// File: frequent.rs
// Foldmine: Frequent structural motif mining
//
// Algorithm overview:
//   1. Collect "frequent seeds": geometric hashes whose document frequency
//      falls in [min_count, max_count].  Rare hashes (below min_count) are
//      infrequent by the anti-monotone property; extremely common hashes
//      (above max_count) encode generic structural features (helix-helix
//      contacts, etc.) and are excluded.
//
//   2. Grow abstract motif graphs via DFS.  A motif is a connected graph
//      whose nodes are abstract residue positions (0, 1, 2, …) and whose
//      edges carry a geometric hash value.  The support of a motif is
//      estimated as the intersection of the structure-ID lists for all edges.
//      Anti-monotone pruning: once the intersection drops below min_count,
//      all extensions are also infrequent.
//
//   3. Canonical deduplication: every motif is normalised to its
//      lexicographically smallest edge-list over all node-relabelings.
//      The DashMap `seen` set prevents the same canonical motif from being
//      reported by multiple DFS branches.
//
// Performance optimisations (A + B):
//   A. Work-stealing DFS: at shallow depths (< PARALLEL_DEPTH) candidates
//      are spawned as rayon tasks so that unbalanced DFS trees distribute
//      naturally across idle threads.
//   B. Bitset intersection: when the database has ≥ BITSET_THRESHOLD
//      structures, structure-ID lists are stored as dense bit-vectors.
//      Intersection becomes a bitwise AND across u64 words (~64× faster),
//      with SIMD auto-vectorisation exploited by the compiler.

use std::collections::HashMap;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Mutex,
};

use dashmap::DashMap;
use rayon::prelude::*;

use crate::geometry::core::{GeometricHash, HashType};
use crate::index::indextable::FolddiscoIndex;
use crate::utils::convert::map_u8_to_aa;

// ---------------------------------------------------------------------------
// Tuning constants
// ---------------------------------------------------------------------------

/// Use bitsets when the database has at least this many structures.
/// Below this threshold the bitset memory overhead (N/8 bytes per seed)
/// is not worth the cost, and sorted-list intersection is faster.
const BITSET_THRESHOLD: usize = 500_000;

/// Spawn rayon tasks for DFS candidates at depths 0..PARALLEL_DEPTH.
/// Deeper levels run sequentially to avoid excessive task-spawn overhead.
const PARALLEL_DEPTH: usize = 2;

const TOP_K_STRUCTURES: usize = 5;

// ---------------------------------------------------------------------------
// Public data types
// ---------------------------------------------------------------------------

/// A single directed edge in an abstract motif graph.
/// `node_a` → `node_b` with geometric hash `hash`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MotifEdge {
    pub node_a: u8,
    pub node_b: u8,
    pub hash: u32,
}

/// An abstract motif graph: `num_nodes` labeled residue positions connected
/// by a sorted list of `MotifEdge` entries.
///
/// Invariant: `edges` is always kept in sorted order.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MotifGraph {
    pub num_nodes: u8,
    pub edges: Vec<MotifEdge>,
}

/// Type alias that signals a `MotifGraph` has been canonicalized.
pub type CanonicalMotif = MotifGraph;

/// A discovered frequent motif with its support statistics.
pub struct FrequentMotif {
    pub motif_id: usize,
    pub motif: CanonicalMotif,
    pub support: f32,
    pub count: usize,
    /// Up to TOP_K_STRUCTURES structure IDs.
    pub top_structure_ids: Vec<usize>,
    /// Inverse document frequency: Σ log2(N / hash_count_i) over all edges.
    pub motif_idf: f32,
    /// Mean residue count of matching structures (estimated from top_structure_ids).
    pub mean_nres: f32,
    /// Length-adjusted IDF: motif_idf × (mean_nres + 1)^(-0.5).
    pub adj_idf: f32,
}

/// Configuration for the mining algorithm.
pub struct MiningConfig {
    /// Minimum support fraction in [0, 1].  Default: 0.01
    pub min_support: f32,
    /// Maximum frequency ratio in [0, 1].  Default: 0.5
    pub max_freq_ratio: f32,
    /// Maximum number of residue positions per motif.  Default: 6
    pub max_nodes: u8,
    /// Maximum number of edges per motif.  0 = auto (complete graph).
    pub max_edges: usize,
    /// Maximum number of seeds to consider.  0 = no limit.
    pub max_seeds: usize,
    /// Stop after finding this many motifs.  0 = no limit.
    pub max_results: usize,
    /// Number of threads.
    pub threads: usize,
    /// Minimum adj_idf score to include in output.  0.0 = no filter.
    pub min_idf: f32,
}

// ---------------------------------------------------------------------------
// Internal seed representation
// ---------------------------------------------------------------------------

/// A frequent single-edge seed with optional bitset representation.
struct Seed {
    hash: u32,
    /// Sorted structure-ID list.  Always present (needed for small N and
    /// for materialising the top-K structure names on output).
    ids: Vec<usize>,
    /// Dense bitset over structure IDs.  Present only when
    /// `total >= BITSET_THRESHOLD`.
    bits: Option<Vec<u64>>,
}

// ---------------------------------------------------------------------------
// Bitset helpers  (Option B)
// ---------------------------------------------------------------------------

/// Convert a sorted list of structure IDs to a dense bitset.
/// `n_words = ceil(total / 64)`.
fn ids_to_bitset(ids: &[usize], n_words: usize) -> Vec<u64> {
    let mut bits = vec![0u64; n_words];
    for &id in ids {
        bits[id / 64] |= 1u64 << (id % 64);
    }
    bits
}

/// Bitwise AND of two bitsets into `out` (which is resized to `n_words`).
/// Returns the popcount of the result — the number of set bits.
///
/// The compiler auto-vectorises this loop to AVX2/AVX-512 when the target
/// supports it.
#[inline]
fn bitset_and_popcount(a: &[u64], b: &[u64], out: &mut Vec<u64>) -> usize {
    debug_assert_eq!(a.len(), b.len());
    out.clear();
    out.extend(a.iter().zip(b.iter()).map(|(&x, &y)| x & y));
    out.iter().map(|w| w.count_ones() as usize).sum()
}

/// Extract the first `k` set-bit positions from a bitset.
fn bitset_first_k(bits: &[u64], k: usize) -> Vec<usize> {
    let mut out = Vec::new();
    'outer: for (word_idx, &w) in bits.iter().enumerate() {
        let mut w = w;
        while w != 0 {
            let bit = w.trailing_zeros() as usize;
            out.push(word_idx * 64 + bit);
            if out.len() == k {
                break 'outer;
            }
            w &= w - 1; // clear lowest set bit
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Canonicalization helpers
// ---------------------------------------------------------------------------

fn next_permutation(arr: &mut Vec<u8>) -> bool {
    let n = arr.len();
    if n <= 1 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && arr[i - 1] >= arr[i] {
        i -= 1;
    }
    if i == 0 {
        return false;
    }
    let mut j = n - 1;
    while arr[j] <= arr[i - 1] {
        j -= 1;
    }
    arr.swap(i - 1, j);
    arr[i..].reverse();
    true
}

pub fn canonicalize(motif: &MotifGraph) -> CanonicalMotif {
    let n = motif.num_nodes as usize;
    let mut perm: Vec<u8> = (0..n as u8).collect();
    let mut best: Option<Vec<MotifEdge>> = None;

    loop {
        let mut inv = [0u8; 16];
        for (new_label, &old_label) in perm.iter().enumerate() {
            inv[old_label as usize] = new_label as u8;
        }
        let mut relabeled: Vec<MotifEdge> = motif
            .edges
            .iter()
            .map(|e| MotifEdge {
                node_a: inv[e.node_a as usize],
                node_b: inv[e.node_b as usize],
                hash: e.hash,
            })
            .collect();
        relabeled.sort();
        if best.is_none() || relabeled < *best.as_ref().unwrap() {
            best = Some(relabeled);
        }
        if !next_permutation(&mut perm) {
            break;
        }
    }
    CanonicalMotif {
        num_nodes: motif.num_nodes,
        edges: best.unwrap(),
    }
}

// ---------------------------------------------------------------------------
// Sorted-list intersection  (fallback for small N)
// ---------------------------------------------------------------------------

#[inline]
fn sorted_intersect_into(a: &[usize], b: &[usize], out: &mut Vec<usize>) {
    out.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
}

// ---------------------------------------------------------------------------
// Seed collection
// ---------------------------------------------------------------------------

fn collect_seeds(
    index: &FolddiscoIndex,
    min_count: usize,
    max_count: usize,
    total: usize,
) -> Vec<Seed> {
    let n = index.total_hashes;
    let n_words = total.div_ceil(64);
    let use_bits = total >= BITSET_THRESHOLD;

    (0..n)
        .into_par_iter()
        .filter_map(|i| {
            let h = index.loaded_hashes[i];
            let byte_span = if i + 1 < index.loaded_offsets.len() {
                index.loaded_offsets[i + 1] - index.loaded_offsets[i]
            } else {
                0
            };
            if byte_span < min_count {
                return None;
            }
            let ids = index.get_entries(h);
            let df = ids.len();
            if df >= min_count && df <= max_count {
                let bits = use_bits.then(|| ids_to_bitset(&ids, n_words));
                Some(Seed { hash: h, ids, bits })
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// DFS motif growth  (Option A + B combined)
// ---------------------------------------------------------------------------

/// Candidate produced inside the DFS inner loop before spawning/recursing.
struct Candidate {
    motif: MotifGraph,
    ids: Vec<usize>,
    bits: Option<Vec<u64>>,
}

#[allow(clippy::too_many_arguments)]
fn dfs_grow(
    motif: &MotifGraph,
    current_ids: &[usize],
    current_bits: Option<&[u64]>,
    seeds: &[Seed],
    min_count: usize,
    max_nodes: u8,
    max_edges: usize,
    max_results: usize,
    total: usize,
    n_words: usize,
    depth: usize,
    seen: &DashMap<CanonicalMotif, ()>,
    results: &Mutex<Vec<FrequentMotif>>,
    counter: &AtomicUsize,
    buf_ids: &mut Vec<usize>,
    buf_bits: &mut Vec<u64>,
) {
    // Early exit if result limit reached.
    if max_results > 0 && counter.load(Ordering::Relaxed) >= max_results {
        return;
    }

    let canonical = canonicalize(motif);
    if seen.insert(canonical.clone(), ()).is_some() {
        return; // already processed by another thread
    }

    // Record this motif.
    let count = current_ids.len();
    let support = count as f32 / total as f32;
    let top_structure_ids: Vec<usize> = if let Some(cb) = current_bits {
        bitset_first_k(cb, TOP_K_STRUCTURES)
    } else {
        current_ids.iter().take(TOP_K_STRUCTURES).copied().collect()
    };
    let motif_id = counter.fetch_add(1, Ordering::Relaxed);
    results.lock().unwrap().push(FrequentMotif {
        motif_id,
        motif: canonical,
        support,
        count,
        top_structure_ids,
        motif_idf: 0.0,
        mean_nres: 0.0,
        adj_idf: 0.0,
    });

    // Check extension limits.
    let at_node_limit = motif.num_nodes >= max_nodes;
    let at_edge_limit = motif.edges.len() >= max_edges;
    if at_node_limit && at_edge_limit {
        return;
    }

    let next_node = motif.num_nodes;
    // At shallow depths we collect candidates and spawn them in parallel.
    // At deeper depths we recurse sequentially to avoid spawn overhead.
    let use_parallel = depth < PARALLEL_DEPTH;
    let mut candidates: Vec<Candidate> = if use_parallel { Vec::new() } else { Vec::new() };

    for seed in seeds {
        // --- Intersection (Option B: bitsets when available) ---
        let new_count = if let (Some(cb), Some(sb)) = (current_bits, seed.bits.as_deref()) {
            bitset_and_popcount(cb, sb, buf_bits)
        } else {
            sorted_intersect_into(current_ids, &seed.ids, buf_ids);
            buf_ids.len()
        };
        if new_count < min_count {
            continue;
        }

        // Materialise the new sorted IDs and optional bitset for this candidate.
        // We need the sorted IDs for recursion (and for output); we need bits
        // for further bitset intersections in child calls.
        let new_ids: Vec<usize> = if current_bits.is_some() && seed.bits.is_some() {
            // Derive sorted IDs from the intersection bitset.
            bitset_first_k(buf_bits, total)
        } else {
            buf_ids.clone()
        };
        let new_bits: Option<Vec<u64>> = if current_bits.is_some() && seed.bits.is_some() {
            Some(buf_bits.clone())
        } else {
            None
        };

        // --- Case A: new edge from existing node → brand-new node ---
        if !at_node_limit && !at_edge_limit {
            for node_a in 0..motif.num_nodes {
                let mut edges = motif.edges.clone();
                edges.push(MotifEdge { node_a, node_b: next_node, hash: seed.hash });
                edges.sort();
                let candidate = MotifGraph { num_nodes: next_node + 1, edges };
                let cand_canonical = canonicalize(&candidate);
                if !seen.contains_key(&cand_canonical) {
                    if use_parallel {
                        candidates.push(Candidate {
                            motif: cand_canonical,
                            ids: new_ids.clone(),
                            bits: new_bits.clone(),
                        });
                    } else {
                        dfs_grow(
                            &cand_canonical,
                            &new_ids,
                            new_bits.as_deref(),
                            seeds,
                            min_count,
                            max_nodes,
                            max_edges,
                            max_results,
                            total,
                            n_words,
                            depth + 1,
                            seen,
                            results,
                            counter,
                            buf_ids,
                            buf_bits,
                        );
                    }
                }
            }
        }

        // --- Case B: chord edge between two existing nodes ---
        if !at_edge_limit {
            for node_a in 0..motif.num_nodes {
                for node_b in 0..motif.num_nodes {
                    if node_a == node_b {
                        continue;
                    }
                    let edge = MotifEdge { node_a, node_b, hash: seed.hash };
                    if motif.edges.contains(&edge) {
                        continue;
                    }
                    let mut edges = motif.edges.clone();
                    edges.push(edge);
                    edges.sort();
                    let candidate = MotifGraph { num_nodes: motif.num_nodes, edges };
                    let cand_canonical = canonicalize(&candidate);
                    if !seen.contains_key(&cand_canonical) {
                        if use_parallel {
                            candidates.push(Candidate {
                                motif: cand_canonical,
                                ids: new_ids.clone(),
                                bits: new_bits.clone(),
                            });
                        } else {
                            dfs_grow(
                                &cand_canonical,
                                &new_ids,
                                new_bits.as_deref(),
                                seeds,
                                min_count,
                                max_nodes,
                                max_edges,
                                max_results,
                                total,
                                n_words,
                                depth + 1,
                                seen,
                                results,
                                counter,
                                buf_ids,
                                buf_bits,
                            );
                        }
                    }
                }
            }
        }
    }

    // --- Option A: spawn collected candidates in parallel ---
    if use_parallel && !candidates.is_empty() {
        rayon::scope(|s| {
            for cand in candidates {
                s.spawn(move |_| {
                    let mut local_buf_ids: Vec<usize> = Vec::new();
                    let mut local_buf_bits: Vec<u64> = if n_words > 0 {
                        vec![0u64; n_words]
                    } else {
                        Vec::new()
                    };
                    dfs_grow(
                        &cand.motif,
                        &cand.ids,
                        cand.bits.as_deref(),
                        seeds,
                        min_count,
                        max_nodes,
                        max_edges,
                        max_results,
                        total,
                        n_words,
                        depth + 1,
                        seen,
                        results,
                        counter,
                        &mut local_buf_ids,
                        &mut local_buf_bits,
                    );
                });
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub fn mine_frequent_motifs(
    index: &FolddiscoIndex,
    lookup: &[(String, usize, usize, f32, usize)],
    config: &MiningConfig,
) -> Vec<FrequentMotif> {
    let total = lookup.len();
    if total == 0 {
        return vec![];
    }

    let min_count = ((config.min_support * total as f32).ceil() as usize).max(1);
    let max_count = (config.max_freq_ratio * total as f32).floor() as usize;

    if min_count > max_count {
        eprintln!(
            "Warning: min_count ({}) > max_count ({}). No seeds will be found.",
            min_count, max_count
        );
        return vec![];
    }

    // Phase 1: collect frequent single-edge seeds.
    let mut seeds = collect_seeds(index, min_count, max_count, total);
    seeds.sort_by_key(|s| s.ids.len()); // rarest first

    if seeds.is_empty() {
        eprintln!("No frequent seeds found. Try lowering --min-support.");
        return vec![];
    }

    // Apply max_seeds limit (rarest seeds first, already sorted).
    if config.max_seeds > 0 && seeds.len() > config.max_seeds {
        seeds.truncate(config.max_seeds);
    }

    let n_words = if total >= BITSET_THRESHOLD { total.div_ceil(64) } else { 0 };

    let max_edges = if config.max_edges > 0 {
        config.max_edges
    } else {
        let k = config.max_nodes as usize;
        k * (k - 1) / 2
    };
    let max_results = config.max_results;

    // Phase 2: parallel per-seed DFS (seed-level parallelism via par_iter;
    // within each DFS tree, work-stealing via rayon::scope at shallow depth).
    let seen: DashMap<CanonicalMotif, ()> = DashMap::new();
    let results: Mutex<Vec<FrequentMotif>> = Mutex::new(Vec::new());
    let counter = AtomicUsize::new(0);

    seeds.par_iter().for_each(|seed| {
        let root = MotifGraph {
            num_nodes: 2,
            edges: vec![MotifEdge { node_a: 0, node_b: 1, hash: seed.hash }],
        };
        let mut buf_ids: Vec<usize> = Vec::with_capacity(seed.ids.len());
        let mut buf_bits: Vec<u64> = if n_words > 0 { vec![0u64; n_words] } else { Vec::new() };
        dfs_grow(
            &root,
            &seed.ids,
            seed.bits.as_deref(),
            &seeds,
            min_count,
            config.max_nodes,
            max_edges,
            max_results,
            total,
            n_words,
            0, // depth = 0 at seed root
            &seen,
            &results,
            &counter,
            &mut buf_ids,
            &mut buf_bits,
        );
    });

    let mut output = results.into_inner().unwrap();

    // Phase 3: compute IDF scores for every motif.
    // Build hash → document-count map from seeds (already have ids.len()).
    let hash_count_map: HashMap<u32, usize> =
        seeds.iter().map(|s| (s.hash, s.ids.len())).collect();

    for m in &mut output {
        // motif_idf: Σ log2(N / hash_count_i) over all edges.
        m.motif_idf = m
            .motif
            .edges
            .iter()
            .map(|e| {
                let c = hash_count_map.get(&e.hash).copied().unwrap_or(1);
                (total as f32 / c as f32).log2().max(0.0)
            })
            .sum();

        // mean_nres: estimated from top_structure_ids (≤ TOP_K_STRUCTURES sample).
        if !m.top_structure_ids.is_empty() {
            let sum_nres: f32 = m
                .top_structure_ids
                .iter()
                .filter_map(|&id| lookup.get(id).map(|e| e.2 as f32))
                .sum();
            m.mean_nres = sum_nres / m.top_structure_ids.len() as f32;
        }

        // adj_idf: length-penalized score (length penalty exponent = 0.5).
        m.adj_idf = m.motif_idf * (m.mean_nres + 1.0).powf(-0.5);
    }

    // Apply min_idf filter.
    if config.min_idf > 0.0 {
        output.retain(|m| m.adj_idf >= config.min_idf);
    }

    // Sort: adj_idf descending → support descending → fewer edges first.
    output.sort_by(|a, b| {
        b.adj_idf
            .partial_cmp(&a.adj_idf)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.support.partial_cmp(&a.support).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.motif.edges.len().cmp(&b.motif.edges.len()))
    });

    // Reassign motif_id in final sorted order.
    for (i, r) in output.iter_mut().enumerate() {
        r.motif_id = i;
    }
    output
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn decode_hash(hash: u32, hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> String {
    let gh = GeometricHash::from_u32(hash, hash_type);
    let mut feat = vec![0.0f32; 16];
    gh.reverse_hash(nbin_dist, nbin_angle, &mut feat);
    let aa1 = map_u8_to_aa(feat[0] as u8);
    let aa2 = map_u8_to_aa(feat[1] as u8);
    let ca_dist = feat[2];
    let cb_dist = feat[3];
    let angle = feat[4];
    if feat[5] != 0.0 || feat[6] != 0.0 {
        format!("{}/{}/{:.2}/{:.2}/{:.1}/{:.1}/{:.1}", aa1, aa2, ca_dist, cb_dist, angle, feat[5], feat[6])
    } else {
        format!("{}/{}/{:.2}/{:.2}/{:.1}", aa1, aa2, ca_dist, cb_dist, angle)
    }
}

pub fn format_motif_tsv_row(
    motif: &FrequentMotif,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
    lookup: &[(String, usize, usize, f32, usize)],
) -> String {
    let edge_strs: Vec<String> = motif
        .motif
        .edges
        .iter()
        .map(|e| {
            let decoded = decode_hash(e.hash, hash_type, nbin_dist, nbin_angle);
            format!("{}-{}:{:08x}:{}", e.node_a, e.node_b, e.hash, decoded)
        })
        .collect();
    let top_names: Vec<String> = motif
        .top_structure_ids
        .iter()
        .filter_map(|&id| lookup.get(id).map(|(name, ..)| name.clone()))
        .collect();
    format!(
        "{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{:.1}\t{:.4}\t{}\t{}",
        motif.motif_id,
        motif.motif.num_nodes,
        motif.motif.edges.len(),
        motif.support,
        motif.count,
        motif.motif_idf,
        motif.mean_nres,
        motif.adj_idf,
        edge_strs.join(";"),
        top_names.join(","),
    )
}

pub fn write_results<W: std::io::Write>(
    results: &[FrequentMotif],
    writer: &mut W,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
    lookup: &[(String, usize, usize, f32, usize)],
) -> std::io::Result<()> {
    writeln!(writer, "motif_id\tnum_residues\tnum_edges\tsupport\tcount\tmotif_idf\tmean_nres\tadj_idf\tedges\ttop_structures")?;
    for motif in results {
        writeln!(writer, "{}", format_motif_tsv_row(motif, hash_type, nbin_dist, nbin_angle, lookup))?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalize_single_edge() {
        let m = MotifGraph {
            num_nodes: 2,
            edges: vec![MotifEdge { node_a: 1, node_b: 0, hash: 42 }],
        };
        let c = canonicalize(&m);
        assert_eq!(c.edges, vec![MotifEdge { node_a: 0, node_b: 1, hash: 42 }]);
    }

    #[test]
    fn test_canonicalize_two_edges() {
        let m = MotifGraph {
            num_nodes: 3,
            edges: vec![
                MotifEdge { node_a: 0, node_b: 1, hash: 10 },
                MotifEdge { node_a: 1, node_b: 2, hash: 20 },
            ],
        };
        let c = canonicalize(&m);
        assert_eq!(c.num_nodes, 3);
        assert_eq!(c.edges.len(), 2);
    }

    #[test]
    fn test_sorted_intersect() {
        let a = vec![1, 3, 5, 7, 9];
        let b = vec![2, 3, 5, 8, 9];
        let mut out = Vec::new();
        sorted_intersect_into(&a, &b, &mut out);
        assert_eq!(out, vec![3, 5, 9]);
    }

    #[test]
    fn test_next_permutation() {
        let mut v = vec![0u8, 1, 2];
        assert!(next_permutation(&mut v));
        assert_eq!(v, vec![0, 2, 1]);
        assert!(next_permutation(&mut v));
        assert_eq!(v, vec![1, 0, 2]);
        let mut count = 3;
        while next_permutation(&mut v) { count += 1; }
        assert_eq!(count, 6);
    }

    #[test]
    fn test_bitset_roundtrip() {
        let ids = vec![0usize, 1, 63, 64, 65, 127, 200];
        let n_words = (200 / 64) + 1 + 1; // enough words to cover id=200
        let bits = ids_to_bitset(&ids, n_words);
        let recovered = bitset_first_k(&bits, 1000);
        assert_eq!(recovered, ids);
    }

    #[test]
    fn test_bitset_and_popcount() {
        let a_ids = vec![1usize, 3, 5, 7, 9];
        let b_ids = vec![2usize, 3, 5, 8, 9];
        let n_words = 1;
        let a = ids_to_bitset(&a_ids, n_words);
        let b = ids_to_bitset(&b_ids, n_words);
        let mut out = vec![0u64; n_words];
        let count = bitset_and_popcount(&a, &b, &mut out);
        assert_eq!(count, 3); // {3, 5, 9}
        let result_ids = bitset_first_k(&out, 1000);
        assert_eq!(result_ids, vec![3usize, 5, 9]);
    }

    #[test]
    fn test_bitset_first_k() {
        let ids = vec![0usize, 1, 2, 100, 200];
        let n_words = 200usize.div_ceil(64) + 1;
        let bits = ids_to_bitset(&ids, n_words);
        let top3 = bitset_first_k(&bits, 3);
        assert_eq!(top3, vec![0usize, 1, 2]);
    }
}
