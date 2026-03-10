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

use std::collections::{HashMap, HashSet};
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Mutex,
};

use dashmap::DashMap;
use rayon::prelude::*;

use crate::controller::feature::get_single_feature;
use crate::controller::io::read_compact_structure;
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::indextable::FolddiscoIndex;
use crate::structure::core::CompactStructure;
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
    /// Residue indices (PDB serial numbers) for each top structure, one per motif node.
    /// `top_structure_residues[i]` corresponds to `top_structure_ids[i]`.
    /// An empty inner Vec means annotation failed for that structure.
    pub top_structure_residues: Vec<Vec<usize>>,
    /// Inverse document frequency: Σ log2(N / hash_count_i) over all edges.
    pub motif_idf: f32,
    /// Mean residue count of matching structures (estimated from top_structure_ids).
    pub mean_nres: f32,
    /// Length-adjusted IDF: motif_idf × (mean_nres + 1)^(-0.5).
    pub adj_idf: f32,
    /// Full sorted structure-ID support list.
    /// Populated during DFS only when merge_iso_threshold is Some; cleared after merge.
    pub all_structure_ids: Vec<usize>,
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
    /// Only report motifs where every unordered node pair has ≥1 directed edge
    /// (fully-connected undirected view).  Default: false.
    pub require_complete: bool,
    /// Distance deviation (Å) for fuzzy seed expansion.  0.0 = disabled.
    /// Paper default: 0.5 Å.
    pub fuzzy_dist: f32,
    /// Angle deviation (degrees) for fuzzy seed expansion.  0.0 = disabled.
    /// Paper default: 5.0°.  Internally converted to radians for sin/cos hash types.
    pub fuzzy_angle: f32,
    /// Geometric hash type (copied from the index configuration).
    pub hash_type: HashType,
    /// Number of distance bins (copied from the index configuration).
    pub num_bin_dist: usize,
    /// Number of angle bins (copied from the index configuration).
    pub num_bin_angle: usize,
    /// If Some(t), merge motifs whose structure-support sets have Jaccard ≥ t after DFS.
    /// None = disabled.  Recommended starting value: 0.7.
    pub merge_iso_threshold: Option<f32>,
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
// Sorted-list union  (fuzzy seed expansion)
// ---------------------------------------------------------------------------

#[inline]
fn sorted_union_into(a: &[usize], b: &[usize], out: &mut Vec<usize>) {
    out.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal   => { out.push(a[i]); i += 1; j += 1; }
            std::cmp::Ordering::Less    => { out.push(a[i]); i += 1; }
            std::cmp::Ordering::Greater => { out.push(b[j]); j += 1; }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
}

// ---------------------------------------------------------------------------
// Non-allocating sorted-list intersection count  (for Jaccard merging)
// ---------------------------------------------------------------------------

/// Returns |A ∩ B| for two sorted slices, without allocating.
/// Used by `merge_isomorphic_motifs` to compute Jaccard similarity.
#[inline]
fn sorted_intersect_count(a: &[usize], b: &[usize]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal   => { n += 1; i += 1; j += 1; }
            std::cmp::Ordering::Less    => { i += 1; }
            std::cmp::Ordering::Greater => { j += 1; }
        }
    }
    n
}

// ---------------------------------------------------------------------------
// Fully-connected motif check  (--require-complete)
// ---------------------------------------------------------------------------

/// Returns `true` iff every unordered node pair {a, b} in `motif` has at
/// least one directed edge (a→b or b→a).
///
/// A 2-node motif is trivially satisfied.
fn is_fully_connected(motif: &MotifGraph) -> bool {
    let n = motif.num_nodes as usize;
    for a in 0..n {
        for b in (a + 1)..n {
            if !motif.edges.iter().any(|e| {
                (e.node_a as usize == a && e.node_b as usize == b)
                    || (e.node_a as usize == b && e.node_b as usize == a)
            }) {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Isomorphic motif merging  (--merge-iso)
// ---------------------------------------------------------------------------

/// Merge motifs whose structure-support sets have Jaccard similarity ≥ `threshold`.
///
/// Algorithm:
/// 1. Group motifs by (num_nodes, num_edges) — only same-shape motifs can be isomorphic.
/// 2. Within each group, run union-find: union(i, j) when Jaccard(i, j) ≥ threshold.
/// 3. For each cluster, keep the motif with the highest `count` as representative.
///    Union all `all_structure_ids`, update count/support/top_structure_ids.
///    Reset `top_structure_residues` so annotate_top_structures re-annotates correctly.
/// 4. Retain only cluster representatives; clear `all_structure_ids`.
fn merge_isomorphic_motifs(results: &mut Vec<FrequentMotif>, total: usize, threshold: f32) {
    if results.is_empty() {
        return;
    }

    // --- Simple path-compressed union-find ---
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    fn union(parent: &mut Vec<usize>, x: usize, y: usize) {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx != ry {
            parent[ry] = rx; // merge ry → rx
        }
    }

    let n = results.len();
    let mut parent: Vec<usize> = (0..n).collect();

    // Group indices by (num_nodes, num_edges).
    let mut groups: HashMap<(u8, usize), Vec<usize>> = HashMap::new();
    for (idx, m) in results.iter().enumerate() {
        groups
            .entry((m.motif.num_nodes, m.motif.edges.len()))
            .or_default()
            .push(idx);
    }

    // Pairwise Jaccard within each group.
    for group in groups.values() {
        for (ii, &i) in group.iter().enumerate() {
            for &j in &group[ii + 1..] {
                let a = &results[i].all_structure_ids;
                let b = &results[j].all_structure_ids;
                if a.is_empty() || b.is_empty() {
                    continue;
                }
                let inter = sorted_intersect_count(a, b);
                let union_sz = a.len() + b.len() - inter;
                if union_sz > 0 && inter as f32 / union_sz as f32 >= threshold {
                    union(&mut parent, i, j);
                }
            }
        }
    }

    // Collect clusters: root → list of member indices.
    let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        clusters.entry(root).or_default().push(i);
    }

    // For singleton clusters, just clear all_structure_ids and move on.
    // For multi-member clusters, merge into the member with the highest count.
    let mut keep: Vec<bool> = vec![false; n];
    let mut merge_buf: Vec<usize> = Vec::new();
    let mut tmp_buf: Vec<usize> = Vec::new();

    for members in clusters.values() {
        if members.len() == 1 {
            let rep = members[0];
            keep[rep] = true;
            results[rep].all_structure_ids.clear();
            continue;
        }

        // Pick representative = member with highest count.
        let rep = *members
            .iter()
            .max_by_key(|&&m| results[m].count)
            .unwrap();
        keep[rep] = true;

        // Union all structure ID lists into merge_buf.
        // sorted_union_into(a, b, out): we alternate merge_buf / tmp_buf as output.
        merge_buf.clear();
        for &m in members {
            // SAFETY of borrow: merge_buf and tmp_buf are two distinct Vecs.
            sorted_union_into(merge_buf.as_slice(), &results[m].all_structure_ids, &mut tmp_buf);
            std::mem::swap(&mut merge_buf, &mut tmp_buf);
        }

        // Update representative.
        let merged_count = merge_buf.len();
        results[rep].count = merged_count;
        results[rep].support = merged_count as f32 / total as f32;
        results[rep].top_structure_ids = merge_buf
            .iter()
            .take(TOP_K_STRUCTURES)
            .copied()
            .collect();
        results[rep].top_structure_residues = Vec::new(); // re-annotated later
        results[rep].all_structure_ids.clear();
    }

    // Retain only representatives.
    let mut i = 0;
    results.retain(|_| { let r = keep[i]; i += 1; r });
}

// ---------------------------------------------------------------------------
// Fuzzy neighbor hash computation  (--fuzzy-dist / --fuzzy-angle)
// ---------------------------------------------------------------------------

/// Given a hash value, compute all neighbor hashes obtained by independently
/// applying ±`dev_dist` (Å) to each distance feature and ±`dev_angle` (°) to
/// each angle feature.
///
/// The approach: decode hash → continuous feature values → perturb → re-encode.
///
/// **Unit note**: `GeometricHash::reverse_hash` returns angle features in
/// degrees.  For hash types that encode angles as sin/cos pairs (PDBTrRosetta,
/// TrRosetta, PDBMotifSinCos), `perfect_hash_as_u32` expects those features
/// back in *radians*, so we convert before re-encoding.  For types that store
/// raw degree values (PDBMotif, FolddiscoAngle, FolddiscoDist) no conversion
/// is necessary.
pub fn neighbor_hashes(
    hash: u32,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
    dev_dist: f32,
    dev_angle_deg: f32,
) -> Vec<u32> {
    // Normalise 0 → type-specific defaults so reverse_hash / perfect_hash_as_u32
    // both operate with the same bin counts (0 means "use default" in the index
    // config, but reverse_hash lacks the fallback that perfect_hash has).
    let (nbin_dist, nbin_angle) = effective_nbins(hash_type, nbin_dist, nbin_angle);
    let gh = GeometricHash::from_u32(hash, hash_type);
    let mut decoded = vec![0.0f32; 16];
    gh.reverse_hash(nbin_dist, nbin_angle, &mut decoded);
    // decoded[0,1] = amino acid indices (categorical — do not perturb).
    // decoded[2,3] = distance features (Å).
    // decoded[4..] = angle features (degrees from reverse_hash).

    // Hash types that use sin/cos encoding need angles in radians for re-encoding.
    let angle_to_rad = matches!(
        hash_type,
        HashType::PDBTrRosetta | HashType::TrRosetta | HashType::PDBMotifSinCos
    );

    // Build the base feature vector in the units expected by perfect_hash_as_u32.
    let mut feat_base = decoded.clone();
    if angle_to_rad {
        for slot in feat_base[4..7].iter_mut() {
            *slot = slot.to_radians();
        }
    }

    // Deviation in the units used by perfect_hash_as_u32.
    let dev_angle_enc = if angle_to_rad {
        dev_angle_deg.to_radians()
    } else {
        dev_angle_deg
    };

    let mut out = vec![hash];

    // Perturb distance features (indices 2, 3).
    if dev_dist > 0.0 {
        for idx in [2usize, 3] {
            for &dev in &[dev_dist, -dev_dist] {
                let mut f = feat_base.clone();
                f[idx] += dev;
                out.push(GeometricHash::perfect_hash_as_u32(
                    &f, hash_type, nbin_dist, nbin_angle,
                ));
            }
        }
    }

    // Perturb angle features (indices 4, 5, 6 — skip if both decoded and
    // encoded values are zero, indicating an unused feature slot).
    if dev_angle_enc > 0.0 {
        for idx in [4usize, 5, 6] {
            if decoded[idx] == 0.0 && feat_base[idx] == 0.0 {
                continue;
            }
            for &dev in &[dev_angle_enc, -dev_angle_enc] {
                let mut f = feat_base.clone();
                f[idx] += dev;
                out.push(GeometricHash::perfect_hash_as_u32(
                    &f, hash_type, nbin_dist, nbin_angle,
                ));
            }
        }
    }

    out.sort_unstable();
    out.dedup();
    out
}

// ---------------------------------------------------------------------------
// Bin-count normalisation
// ---------------------------------------------------------------------------

/// Return the effective (nbin_dist, nbin_angle) for `hash_type`, replacing
/// zero with the type-specific default that `perfect_hash_as_u32` uses
/// internally.  This ensures `reverse_hash` — which lacks the same fallback
/// logic — receives valid bin counts and produces physically sensible feature
/// values.
///
/// Defaults match the guards inside each hash type's `perfect_hash` function:
/// | Type              | dist | angle |
/// |-------------------|------|-------|
/// | PDBMotif          |  18  |   9   |
/// | PDBMotifSinCos    |   8  |   3   |
/// | TrRosetta         |   8  |   3   |
/// | PDBTrRosetta      |  16  |   4   |
/// | PointPairFeature  |   8  |   3   |
/// | TertiaryInteraction|  8  |   3   |
/// | Hybrid            |  16  |   4   |
/// | FolddiscoAngle    |   8  |  32   |
/// | FolddiscoDist     |  32  |   8   |
#[inline]
pub fn effective_nbins(hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> (usize, usize) {
    let d = if nbin_dist == 0 {
        match hash_type {
            HashType::PDBMotif                => 18,
            HashType::PDBMotifSinCos          =>  8,
            HashType::TrRosetta               =>  8,
            HashType::PDBTrRosetta            => 16,
            HashType::PointPairFeature        =>  8,
            HashType::TertiaryInteraction     =>  8,
            HashType::Hybrid                  => 16,
            HashType::FolddiscoAngle          =>  8,
            HashType::FolddiscoDist           => 32,
            #[allow(unreachable_patterns)]
            _                                 => 16, // safe fallback
        }
    } else {
        nbin_dist
    };
    let a = if nbin_angle == 0 {
        match hash_type {
            HashType::PDBMotif                =>  9,
            HashType::PDBMotifSinCos          =>  3,
            HashType::TrRosetta               =>  3,
            HashType::PDBTrRosetta            =>  4,
            HashType::PointPairFeature        =>  3,
            HashType::TertiaryInteraction     =>  3,
            HashType::Hybrid                  =>  4,
            HashType::FolddiscoAngle          => 32,
            HashType::FolddiscoDist           =>  8,
            #[allow(unreachable_patterns)]
            _                                 => 16, // safe fallback
        }
    } else {
        nbin_angle
    };
    (d, a)
}

// ---------------------------------------------------------------------------
// Seed collection
// ---------------------------------------------------------------------------

fn collect_seeds(
    index: &FolddiscoIndex,
    min_count: usize,
    max_count: usize,
    total: usize,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
    fuzzy_dist: f32,
    fuzzy_angle: f32,
) -> Vec<Seed> {
    let n = index.total_hashes;
    let n_words = total.div_ceil(64);
    let use_bits = total >= BITSET_THRESHOLD;
    let do_fuzzy = fuzzy_dist > 0.0 || fuzzy_angle > 0.0;

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
            let mut ids = index.get_entries(h);

            // Fuzzy expansion: union in structure IDs from neighboring hashes.
            if do_fuzzy {
                let nbrs = neighbor_hashes(
                    h, hash_type, nbin_dist, nbin_angle, fuzzy_dist, fuzzy_angle,
                );
                let mut buf = Vec::new();
                for nbr_h in nbrs {
                    if nbr_h == h {
                        continue;
                    }
                    let nbr_ids = index.get_entries(nbr_h);
                    if !nbr_ids.is_empty() {
                        sorted_union_into(&ids, &nbr_ids, &mut buf);
                        std::mem::swap(&mut ids, &mut buf);
                    }
                }
            }

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
    require_complete: bool,
    store_all_ids: bool,
    buf_ids: &mut Vec<usize>,
    buf_bits: &mut Vec<u64>,
) {
    // Early exit if result limit reached (counts only motifs actually recorded).
    if max_results > 0 && counter.load(Ordering::Relaxed) >= max_results {
        return;
    }

    let canonical = canonicalize(motif);
    if seen.insert(canonical.clone(), ()).is_some() {
        return; // already processed by another thread
    }

    // Record this motif — skipped when --require-complete is on and the
    // motif graph is not fully connected (every unordered pair has ≥1 edge).
    if !require_complete || is_fully_connected(&canonical) {
        let count = current_ids.len();
        let support = count as f32 / total as f32;
        let top_structure_ids: Vec<usize> = if let Some(cb) = current_bits {
            bitset_first_k(cb, TOP_K_STRUCTURES)
        } else {
            current_ids.iter().take(TOP_K_STRUCTURES).copied().collect()
        };
        let all_structure_ids = if store_all_ids {
            current_ids.to_vec()
        } else {
            Vec::new()
        };
        let motif_id = counter.fetch_add(1, Ordering::Relaxed);
        results.lock().unwrap().push(FrequentMotif {
            motif_id,
            motif: canonical,
            support,
            count,
            top_structure_ids,
            top_structure_residues: Vec::new(), // filled in by annotate_top_structures()
            motif_idf: 0.0,
            mean_nres: 0.0,
            adj_idf: 0.0,
            all_structure_ids,
        });
    }

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
                            require_complete,
                            store_all_ids,
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
                    // Prevent adding a second edge in the same direction between
                    // a node pair that already has one (regardless of hash).
                    // Without this, the same "0-1" pair can appear multiple times
                    // in a motif's edge list with different hash values — an
                    // over-specified, artifact constraint.
                    if motif.edges.iter().any(|e| e.node_a == edge.node_a && e.node_b == edge.node_b) {
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
                                require_complete,
                                store_all_ids,
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
                        require_complete,
                        store_all_ids,
                        &mut local_buf_ids,
                        &mut local_buf_bits,
                    );
                });
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Residue annotation helpers
// ---------------------------------------------------------------------------

/// CA–CA distance cutoff used when scanning a structure for hash matches.
const HASH_DIST_CUTOFF: f32 = 20.0;

/// Maximum candidates per edge before truncation (bounds backtracking cost).
const MAX_CANDIDATES: usize = 500;

/// Given a loaded `CompactStructure` and a `MotifGraph`, find the first
/// consistent assignment of motif nodes to actual residues in the structure.
///
/// Returns a `Vec<usize>` of PDB residue serial numbers, one per abstract
/// node (index 0..motif.num_nodes-1), or an empty Vec if no match is found.
fn find_motif_residues(
    compact: &CompactStructure,
    motif: &MotifGraph,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
) -> Vec<usize> {
    let n = motif.num_nodes as usize;
    let target_hashes: HashSet<u32> = motif.edges.iter().map(|e| e.hash).collect();

    // Single pass over all ordered residue pairs to collect candidates.
    let mut hash_to_pairs: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    let mut feature = vec![0.0f32; 16];

    for i in 0..compact.num_residues {
        for j in 0..compact.num_residues {
            if i == j {
                continue;
            }
            let Some(ca_dist) = compact.get_ca_distance(i, j) else {
                continue;
            };
            if ca_dist > HASH_DIST_CUTOFF {
                continue;
            }
            if !get_single_feature(i, j, compact, hash_type, HASH_DIST_CUTOFF, &mut feature) {
                continue;
            }
            let h = GeometricHash::perfect_hash_as_u32(&feature, hash_type, nbin_dist, nbin_angle);
            if target_hashes.contains(&h) {
                let pairs = hash_to_pairs.entry(h).or_default();
                if pairs.len() < MAX_CANDIDATES {
                    pairs.push((i, j));
                }
            }
        }
    }

    // Sort edge processing order: rarest hash first (fewest candidates) for
    // maximum pruning efficiency in the backtracking search.
    let empty: Vec<(usize, usize)> = Vec::new();
    let mut order: Vec<usize> = (0..motif.edges.len()).collect();
    order.sort_by_key(|&ei| {
        hash_to_pairs.get(&motif.edges[ei].hash).map(|v| v.len()).unwrap_or(0)
    });

    let ordered_edges: Vec<&MotifEdge> =
        order.iter().map(|&ei| &motif.edges[ei]).collect();
    let ordered_cands: Vec<&[(usize, usize)]> = order
        .iter()
        .map(|&ei| {
            hash_to_pairs
                .get(&motif.edges[ei].hash)
                .map(|v| v.as_slice())
                .unwrap_or(empty.as_slice())
        })
        .collect();

    // Backtracking search — assignment[i] = structure residue index for node i.
    let mut assignment = vec![usize::MAX; n];
    if bt_assign(&ordered_edges, &ordered_cands, &mut assignment, 0, n) {
        // Convert 0-based residue array indices to PDB serial numbers.
        assignment
            .iter()
            .map(|&idx| compact.residue_serial[idx] as usize)
            .collect()
    } else {
        vec![]
    }
}

/// Backtracking constraint solver: assigns a unique structure residue to each
/// motif node such that every edge constraint is satisfied.
///
/// `edges` and `cands` are parallel slices, both ordered rarest-first.
/// `assignment[k] == usize::MAX` means node k is not yet assigned.
fn bt_assign(
    edges: &[&MotifEdge],
    cands: &[&[(usize, usize)]],
    assignment: &mut Vec<usize>,
    ei: usize,
    n: usize,
) -> bool {
    if ei == edges.len() {
        // All edges processed; check every node has been assigned.
        return assignment.iter().take(n).all(|&x| x != usize::MAX);
    }

    let na = edges[ei].node_a as usize;
    let nb = edges[ei].node_b as usize;
    let was_na = assignment[na];
    let was_nb = assignment[nb];

    for &(ri, rj) in cands[ei] {
        // Must be consistent with already-fixed nodes.
        if was_na != usize::MAX && was_na != ri { continue; }
        if was_nb != usize::MAX && was_nb != rj { continue; }
        // Residues must be distinct and not already used by another node.
        if ri == rj { continue; }
        if (0..n).any(|k| k != na && assignment[k] == ri) { continue; }
        if (0..n).any(|k| k != nb && assignment[k] == rj) { continue; }

        assignment[na] = ri;
        assignment[nb] = rj;

        if bt_assign(edges, cands, assignment, ei + 1, n) {
            return true;
        }

        // Restore only what this call changed.
        assignment[na] = was_na;
        assignment[nb] = was_nb;
    }
    false
}

/// Post-process mined motifs: for each top structure, load the structure file
/// and find the actual residue positions that realise the motif.
///
/// Structures are cached by file path so each unique structure is loaded at
/// most once.  If a structure file is missing or unreadable, the motif's
/// `top_structure_residues` entry for that structure is left empty and the
/// bare structure name is shown in the output (no error is printed).
///
/// **Residue annotation requires the original structure files to be
/// accessible at the paths stored in the index lookup table.**  When the
/// index was built on a different machine or the files have been moved,
/// annotation is silently skipped and only bare names are shown.
pub fn annotate_top_structures(
    motifs: &mut [FrequentMotif],
    lookup: &[(String, usize, usize, f32, usize)],
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
) {
    // Cache keyed on file path: Some(structure) if loaded OK, None if unavailable.
    let mut cache: HashMap<String, Option<CompactStructure>> = HashMap::new();

    for motif in motifs.iter_mut() {
        motif.top_structure_residues = motif
            .top_structure_ids
            .iter()
            .map(|&id| {
                if id >= lookup.len() {
                    return vec![];
                }
                let path = lookup[id].0.clone();
                let entry = cache.entry(path.clone()).or_insert_with(|| {
                    // Guard: avoid calling read_compact_structure when the file
                    // does not exist — that function uses .expect() internally
                    // and would panic rather than returning Err.
                    if !std::path::Path::new(&path).exists() {
                        return None;
                    }
                    // std::panic::catch_unwind guards against unexpected panics
                    // (e.g., corrupted files, unsupported format) so that a
                    // single bad file does not abort the entire annotation pass.
                    std::panic::catch_unwind(|| {
                        read_compact_structure(&path).ok().map(|(c, _)| c)
                    })
                    .unwrap_or(None)
                });
                match entry {
                    Some(compact) => find_motif_residues(
                        compact,
                        &motif.motif,
                        hash_type,
                        nbin_dist,
                        nbin_angle,
                    ),
                    None => vec![],
                }
            })
            .collect();
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
    let mut seeds = collect_seeds(
        index,
        min_count,
        max_count,
        total,
        config.hash_type,
        config.num_bin_dist,
        config.num_bin_angle,
        config.fuzzy_dist,
        config.fuzzy_angle,
    );
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
            config.require_complete,
            config.merge_iso_threshold.is_some(),
            &mut buf_ids,
            &mut buf_bits,
        );
    });

    let mut output = results.into_inner().unwrap();

    // Phase 2.5: merge isomorphic motifs (--merge-iso).
    // Done before IDF scoring so we only compute IDF for surviving merged motifs.
    if let Some(threshold) = config.merge_iso_threshold {
        merge_isomorphic_motifs(&mut output, total, threshold);
    }

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
    // Normalise 0 → type-specific defaults before calling reverse_hash.
    // Indices built with "use default" bins store 0 in their .type config file,
    // but reverse_hash lacks the automatic fallback that perfect_hash has,
    // which would otherwise produce physically invalid (e.g. negative) values.
    let (nbin_dist, nbin_angle) = effective_nbins(hash_type, nbin_dist, nbin_angle);
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
        .enumerate()
        .filter_map(|(pos, &id)| {
            lookup.get(id).map(|(name, ..)| {
                match motif.top_structure_residues.get(pos) {
                    Some(residues) if !residues.is_empty() => {
                        let res_str = residues
                            .iter()
                            .map(|r| r.to_string())
                            .collect::<Vec<_>>()
                            .join("-");
                        format!("{}:{}", name, res_str)
                    }
                    _ => name.clone(),
                }
            })
        })
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

    #[test]
    fn test_sorted_union_into() {
        let a = vec![1usize, 3, 5, 7];
        let b = vec![2usize, 3, 6, 7, 8];
        let mut out = Vec::new();
        sorted_union_into(&a, &b, &mut out);
        assert_eq!(out, vec![1, 2, 3, 5, 6, 7, 8]);

        // Union with empty
        let mut out2 = Vec::new();
        sorted_union_into(&a, &[], &mut out2);
        assert_eq!(out2, a);

        sorted_union_into(&[], &b, &mut out2);
        assert_eq!(out2, b);
    }

    #[test]
    fn test_is_fully_connected() {
        // 2-node motif: always connected
        let m2 = MotifGraph {
            num_nodes: 2,
            edges: vec![MotifEdge { node_a: 0, node_b: 1, hash: 1 }],
        };
        assert!(is_fully_connected(&m2));

        // 3-node motif: triangle (all pairs covered)
        let m3_full = MotifGraph {
            num_nodes: 3,
            edges: vec![
                MotifEdge { node_a: 0, node_b: 1, hash: 1 },
                MotifEdge { node_a: 1, node_b: 2, hash: 2 },
                MotifEdge { node_a: 0, node_b: 2, hash: 3 },
            ],
        };
        assert!(is_fully_connected(&m3_full));

        // 3-node motif: path (pair {0,2} missing)
        let m3_path = MotifGraph {
            num_nodes: 3,
            edges: vec![
                MotifEdge { node_a: 0, node_b: 1, hash: 1 },
                MotifEdge { node_a: 1, node_b: 2, hash: 2 },
            ],
        };
        assert!(!is_fully_connected(&m3_path));

        // 3-node motif: reverse direction counts as covering the pair
        let m3_rev = MotifGraph {
            num_nodes: 3,
            edges: vec![
                MotifEdge { node_a: 0, node_b: 1, hash: 1 },
                MotifEdge { node_a: 1, node_b: 2, hash: 2 },
                MotifEdge { node_a: 2, node_b: 0, hash: 3 }, // covers {0,2} via 2→0
            ],
        };
        assert!(is_fully_connected(&m3_rev));
    }

    /// Regression test for the --require-complete bug report.
    ///
    /// The user observed a 4-node motif in the output with edges
    ///   0→1, 0→2, 0→3, 1→0, 1→2, 1→3
    /// and claimed --require-complete did not filter it.  The pair {2,3} has no
    /// edge, so is_fully_connected MUST return false for this graph.
    ///
    /// If this test fails, the is_fully_connected implementation is buggy.
    /// If this test passes but the user still sees the motif in output, the
    /// likely cause is either an old binary or require_complete not reaching
    /// dfs_grow (both of which are runtime / deployment issues, not logic bugs).
    #[test]
    fn test_is_fully_connected_4node_missing_pair() {
        // Exact hashes from the user's report
        let m = MotifGraph {
            num_nodes: 4,
            edges: vec![
                MotifEdge { node_a: 0, node_b: 1, hash: 0x009dfd21 },
                MotifEdge { node_a: 0, node_b: 2, hash: 0x009dfd21 },
                MotifEdge { node_a: 0, node_b: 3, hash: 0x009dfd21 },
                MotifEdge { node_a: 1, node_b: 0, hash: 0x120dfd12 },
                MotifEdge { node_a: 1, node_b: 2, hash: 0x009dfd21 },
                MotifEdge { node_a: 1, node_b: 3, hash: 0x120dfd12 },
            ],
        };
        // Pair {2,3} has no covering edge — must be rejected.
        assert!(
            !is_fully_connected(&m),
            "4-node motif with missing pair {{2,3}} must NOT pass is_fully_connected"
        );

        // Adding a 2→3 edge completes all pairs — must be accepted.
        let mut edges_complete = m.edges.clone();
        edges_complete.push(MotifEdge { node_a: 2, node_b: 3, hash: 0x11111111 });
        let m_complete = MotifGraph { num_nodes: 4, edges: edges_complete };
        assert!(
            is_fully_connected(&m_complete),
            "4-node motif with all pairs covered MUST pass is_fully_connected"
        );
    }

    /// Test that effective_nbins fills in sensible defaults for the two most
    /// common hash types used in practice.
    #[test]
    fn test_effective_nbins_defaults() {
        use crate::geometry::core::HashType;
        // PDBTrRosetta ("default") — distance bin = 16, sin/cos bin = 4.
        let (d, a) = effective_nbins(HashType::PDBTrRosetta, 0, 0);
        assert_eq!(d, 16, "PDBTrRosetta default dist bins");
        assert_eq!(a,  4, "PDBTrRosetta default angle bins");

        // PDBMotif — distance bin = 18, angle bin = 9.
        let (d2, a2) = effective_nbins(HashType::PDBMotif, 0, 0);
        assert_eq!(d2, 18, "PDBMotif default dist bins");
        assert_eq!(a2,  9, "PDBMotif default angle bins");

        // Non-zero values must be passed through unchanged.
        let (d3, a3) = effective_nbins(HashType::PDBTrRosetta, 12, 3);
        assert_eq!(d3, 12);
        assert_eq!(a3,  3);
    }

    /// Smoke-test that decode_hash produces physically reasonable distances
    /// (2–20 Å) when called with nbin_dist = 0 (the "use default" sentinel
    /// stored in some index config files).  Before the effective_nbins fix this
    /// produced large negative values (e.g. -232 Å) because reverse_hash
    /// divided by (0 − 1) = −1 instead of the actual bin count.
    #[test]
    fn test_decode_hash_nbin_zero_no_negative_distances() {
        use crate::geometry::core::HashType;
        // Hash from user's report (PDBTrRosetta index built with default bins).
        let hash: u32 = 0x009dfd21;
        let s = decode_hash(hash, HashType::PDBTrRosetta, 0, 0);
        // Parse ca_dist (3rd '/'-separated field after AA1/AA2).
        let parts: Vec<&str> = s.split('/').collect();
        assert!(parts.len() >= 4, "unexpected edge format: {}", s);
        let ca: f32 = parts[2].parse().expect("ca_dist not a float");
        let cb: f32 = parts[3].parse().expect("cb_dist not a float");
        assert!(ca >= 0.0 && ca <= 25.0,
            "ca_dist out of range [0,25]: {}  (full string: {})", ca, s);
        assert!(cb >= 0.0 && cb <= 25.0,
            "cb_dist out of range [0,25]: {}  (full string: {})", cb, s);
    }

    #[test]
    fn test_neighbor_hashes_returns_original() {
        use crate::geometry::core::HashType;
        // With zero deviations, neighbor_hashes should return only the original hash.
        let hash: u32 = 0x00013d8b; // arbitrary test hash
        let result = neighbor_hashes(hash, HashType::PDBTrRosetta, 16, 4, 0.0, 0.0);
        assert_eq!(result, vec![hash]);
    }

    #[test]
    fn test_neighbor_hashes_produces_neighbors() {
        use crate::geometry::core::HashType;
        // With a large enough deviation (> bin_width/2 ≈ 0.6 Å for 16-bin distance
        // over [2,20] Å), we are guaranteed to cross at least one bin boundary.
        // Use 1.5 Å / 45° to ensure crossing regardless of where the hash center lands.
        let hash: u32 = 0x00013d8b;
        let result = neighbor_hashes(hash, HashType::PDBTrRosetta, 16, 4, 1.5, 45.0);
        // Original hash must be included
        assert!(result.contains(&hash), "original hash missing");
        // At least some neighbors expected with this large deviation
        assert!(result.len() > 1, "expected neighbors with large dev, got {:?}", result);
        // All entries sorted and deduplicated
        for i in 1..result.len() {
            assert!(result[i] > result[i - 1], "not sorted/deduped");
        }
    }
}
