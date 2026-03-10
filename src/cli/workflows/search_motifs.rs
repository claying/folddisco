// File: search_motifs.rs
// Foldmine search: check whether a query motif (PDB + residues) appears in a
// mined results TSV produced by `folddisco mine`.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::cli::{AppArgs, print_logo};
use crate::cli::config::read_index_config_from_file;
use crate::controller::feature::get_single_feature;
use crate::controller::frequent::{effective_nbins, neighbor_hashes};
use crate::controller::io::{get_lookup_and_type, read_compact_structure};
use crate::controller::query::parse_query_string;
use crate::geometry::core::{GeometricHash, HashType};
use crate::prelude::{print_log_msg, INFO};

pub const HELP_SEARCH: &str = "\
usage: folddisco search -f <results_tsv> -p <pdb> -q <residues> [options]

Check whether a query structural motif (defined by a PDB file and residue
positions) is present among the motifs in a `folddisco mine` results file.

required:
    -f, --results <PATH>         Mining results TSV (from `folddisco mine`)
    -p, --pdb <PATH>             Query PDB file
    -q, --query <STRING>         Residue specification, same format as
                                 `folddisco query` (e.g. F207,F212,F225,F229)

hash configuration (must match the index used for mining):
    -i, --index <PREFIX>         Read hash type + bin counts from index config
                                 (overrides --hash-type / --nbin-dist / --nbin-angle)
    --hash-type <TYPE>           Hash type [default: PDBTrRosetta]
                                 Values: PDBMotif, PDBMotifSinCos, TrRosetta,
                                 PDBTrRosetta, PointPairFeature, Hybrid,
                                 FolddiscoAngle, FolddiscoDist
    --nbin-dist <INT>            Number of distance bins (0 = type default) [default: 0]
    --nbin-angle <INT>           Number of angle bins (0 = type default) [default: 0]

matching options:
    --fuzzy-dist <FLOAT>         Distance tolerance (Å) for fuzzy hash expansion.
                                 0 = exact only. Paper default: 0.5 [default: 0]
    --fuzzy-angle <FLOAT>        Angle tolerance (degrees) for fuzzy hash expansion.
                                 0 = exact only. Paper default: 5.0 [default: 0]
    --min-match <FLOAT>          Minimum fraction of motif edges that must be found
                                 in the query hash set [default: 1.0 = all edges]

output:
    -o, --output <PATH>          Output TSV path (default: stdout)

general:
    -v, --verbose                Print verbose messages
    -h, --help                   Print this help menu

output columns:
    All original columns from the results TSV, plus two appended columns:
    match_fraction  fraction of motif edges whose hash appears in query set
    matched_edges   semicolon-separated node-pair labels of matching edges (e.g. 0-1;1-2)
";

// ---------------------------------------------------------------------------
// Helper: build the directed query hash set from a PDB file + residue spec.
// Returns None if the PDB could not be read or no residues could be located.
// ---------------------------------------------------------------------------

fn build_query_hashes(
    pdb_path: &str,
    query_string: &str,
    hash_type: HashType,
    nbin_dist: usize,
    nbin_angle: usize,
    fuzzy_dist: f32,
    fuzzy_angle: f32,
) -> Option<HashSet<u32>> {
    // Resolve effective bin counts (0 → type-specific defaults).
    let (eff_dist, eff_angle) = effective_nbins(hash_type, nbin_dist, nbin_angle);

    // Parse residue specification.
    let (residue_specs, _subs) = parse_query_string(query_string, b'A');
    if residue_specs.is_empty() {
        eprintln!("Warning: no residues parsed from query string '{}'", query_string);
        return None;
    }

    // Load PDB → CompactStructure.
    let (compact, _) = read_compact_structure(pdb_path)
        .map_err(|_| eprintln!("[FAIL] Failed to read PDB file: {}", pdb_path))
        .ok()?;

    // Convert (chain byte, serial u64) → structural indices.
    let indices: Vec<usize> = residue_specs
        .iter()
        .filter_map(|(chain, serial)| compact.get_index(chain, serial))
        .collect();

    if indices.is_empty() {
        eprintln!(
            "Warning: none of the query residues could be located in '{}'",
            pdb_path
        );
        return None;
    }

    // Compute hashes for ALL ordered pairs (both i→j and j→i).
    let mut query_hashes: HashSet<u32> = HashSet::new();
    let mut feat = vec![0.0f32; 9];
    let dist_cutoff = 20.0_f32;

    for i in 0..indices.len() {
        for j in 0..indices.len() {
            if i == j {
                continue;
            }
            if get_single_feature(indices[i], indices[j], &compact, hash_type, dist_cutoff, &mut feat) {
                let h = GeometricHash::perfect_hash_as_u32(&feat, hash_type, eff_dist, eff_angle);
                query_hashes.insert(h);
                // Fuzzy expansion: add neighbor hashes.
                if fuzzy_dist > 0.0 || fuzzy_angle > 0.0 {
                    for nbr in neighbor_hashes(h, hash_type, eff_dist, eff_angle, fuzzy_dist, fuzzy_angle) {
                        query_hashes.insert(nbr);
                    }
                }
            }
        }
    }

    Some(query_hashes)
}

// ---------------------------------------------------------------------------
// Helper: for one TSV data row, extract edge hashes and compute match score.
// Returns None if the row cannot be parsed (too few columns).
// Returns Some((match_fraction, matched_edge_labels)).
// ---------------------------------------------------------------------------

fn match_row(
    row: &str,
    query_hashes: &HashSet<u32>,
    min_match: f32,
) -> Option<(f32, String)> {
    let fields: Vec<&str> = row.split('\t').collect();
    // edges column is the 9th column (index 8, 0-based).
    if fields.len() < 9 {
        return None;
    }
    let edges_field = fields[8];
    if edges_field.is_empty() {
        return None;
    }

    let edge_parts: Vec<&str> = edges_field.split(';').collect();
    let total_edges = edge_parts.len();
    if total_edges == 0 {
        return None;
    }

    let mut matched_labels: Vec<&str> = Vec::new();

    for edge_part in &edge_parts {
        // Format: "nodeA-nodeB:hash_hex:decoded_string"
        let mut cols = edge_part.splitn(3, ':');
        let node_pair = cols.next().unwrap_or("");  // "nodeA-nodeB"
        let hash_hex  = cols.next().unwrap_or("");  // 8-char hex u32
        if hash_hex.is_empty() {
            continue;
        }
        if let Ok(h) = u32::from_str_radix(hash_hex, 16) {
            if query_hashes.contains(&h) {
                matched_labels.push(node_pair);
            }
        }
    }

    let match_fraction = matched_labels.len() as f32 / total_edges as f32;
    if match_fraction >= min_match {
        Some((match_fraction, matched_labels.join(";")))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Main workflow entry point.
// ---------------------------------------------------------------------------

pub fn search_motifs(env: AppArgs) {
    match env {
        AppArgs::Search {
            results_path,
            pdb_path,
            query_string,
            index_path,
            hash_type_str,
            nbin_dist,
            nbin_angle,
            fuzzy_dist,
            fuzzy_angle,
            min_match,
            output,
            verbose,
            help: _,
        } => {
            if verbose {
                print_logo();
            }

            let results_path = match results_path {
                Some(p) => p,
                None => {
                    eprintln!("{}", HELP_SEARCH);
                    std::process::exit(1);
                }
            };
            let pdb_path = match pdb_path {
                Some(p) => p,
                None => {
                    eprintln!("Error: -p / --pdb is required.\n\n{}", HELP_SEARCH);
                    std::process::exit(1);
                }
            };
            let query_string = match query_string {
                Some(q) => q,
                None => {
                    eprintln!("Error: -q / --query is required.\n\n{}", HELP_SEARCH);
                    std::process::exit(1);
                }
            };

            // Determine hash type and bin counts.
            let (hash_type, eff_nbin_dist, eff_nbin_angle) = if let Some(ref ipath) = index_path {
                let (_, config_path) = get_lookup_and_type(ipath);
                let cfg = read_index_config_from_file(&config_path);
                if verbose {
                    print_log_msg(
                        INFO,
                        &format!(
                            "Loaded config from index: hash_type={}, nbin_dist={}, nbin_angle={}",
                            cfg.hash_type.to_string(), cfg.num_bin_dist, cfg.num_bin_angle
                        ),
                    );
                }
                (cfg.hash_type, cfg.num_bin_dist, cfg.num_bin_angle)
            } else {
                let ht = HashType::get_with_str(&hash_type_str);
                if verbose {
                    print_log_msg(
                        INFO,
                        &format!(
                            "Using hash_type={} (nbin_dist={}, nbin_angle={})",
                            ht.to_string(), nbin_dist, nbin_angle
                        ),
                    );
                }
                (ht, nbin_dist, nbin_angle)
            };

            // Build query hash set.
            if verbose {
                print_log_msg(
                    INFO,
                    &format!(
                        "Computing query hashes: pdb={}, residues={}",
                        pdb_path, query_string
                    ),
                );
            }
            let query_hashes = match build_query_hashes(
                &pdb_path,
                &query_string,
                hash_type,
                eff_nbin_dist,
                eff_nbin_angle,
                fuzzy_dist,
                fuzzy_angle,
            ) {
                Some(h) => h,
                None => {
                    eprintln!("Error: could not build query hash set. Check PDB path and residue spec.");
                    std::process::exit(1);
                }
            };

            if verbose {
                print_log_msg(
                    INFO,
                    &format!(
                        "Query hash set: {} unique hashes (fuzzy_dist={}, fuzzy_angle={})",
                        query_hashes.len(), fuzzy_dist, fuzzy_angle
                    ),
                );
            }

            // Open the results TSV for reading.
            let tsv_file = File::open(&results_path).unwrap_or_else(|e| {
                eprintln!("Error opening results file '{}': {}", results_path, e);
                std::process::exit(1);
            });
            let reader = BufReader::new(tsv_file);
            let mut lines = reader.lines();

            // Set up writer.
            let stdout;
            let file_out;
            let mut writer: Box<dyn Write> = match output {
                Some(ref path) => {
                    file_out = File::create(path).unwrap_or_else(|e| {
                        eprintln!("Error creating output file '{}': {}", path, e);
                        std::process::exit(1);
                    });
                    Box::new(BufWriter::new(file_out))
                }
                None => {
                    stdout = std::io::stdout();
                    Box::new(BufWriter::new(stdout.lock()))
                }
            };

            // Write header: original header + two new columns.
            if let Some(Ok(header)) = lines.next() {
                writeln!(writer, "{}\tmatch_fraction\tmatched_edges", header)
                    .unwrap_or_else(|e| eprintln!("Write error: {}", e));
            }

            // Process each data row.
            let mut n_matched = 0usize;
            let mut n_total = 0usize;
            for line in lines {
                let row = match line {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Warning: read error: {}", e);
                        continue;
                    }
                };
                if row.trim().is_empty() {
                    continue;
                }
                n_total += 1;
                if let Some((frac, matched)) = match_row(&row, &query_hashes, min_match) {
                    writeln!(writer, "{}\t{:.4}\t{}", row, frac, matched)
                        .unwrap_or_else(|e| eprintln!("Write error: {}", e));
                    n_matched += 1;
                }
            }

            if verbose {
                print_log_msg(
                    INFO,
                    &format!(
                        "Scanned {} motifs, {} matched (min_match={:.2})",
                        n_total, n_matched, min_match
                    ),
                );
            }
        }
        _ => {
            eprintln!("{}", HELP_SEARCH);
            std::process::exit(1);
        }
    }
}
