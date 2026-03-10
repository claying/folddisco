// File: mine_motifs.rs
// Foldmine workflow: mine frequent structural motifs from an index.

use std::fs::File;
use std::io::BufWriter;

use crate::cli::{AppArgs, print_logo};
use crate::cli::config::read_index_config_from_file;
use crate::controller::frequent::{annotate_top_structures, mine_frequent_motifs, write_results, MiningConfig};
use crate::controller::io::get_lookup_and_type;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::{print_log_msg, INFO};

pub const HELP_MINE: &str = "\
usage: folddisco mine -i <index_prefix> [options]

Foldmine: discover frequent structural motifs in a Folddisco index.

input:
    -i, --index <INDEX_PREFIX>   Path prefix to index files (required)

mining options:
    --min-support <FLOAT>        Minimum support fraction in [0,1] [default: 0.01]
    --max-freq <FLOAT>           Maximum frequency ratio in [0,1] [default: 0.5]
    --max-residues <INT>         Maximum residues per motif [default: 6]
    --max-seeds <INT>            Limit number of seed edges considered (0=unlimited) [default: 0]
    --max-results <INT>          Stop after finding this many motifs (0=unlimited) [default: 0]
    --min-idf <FLOAT>            Minimum length-adjusted IDF score to include (0=no filter) [default: 0]
    --require-complete           Only output motifs where every node pair has at least one
                                 directed edge (fully-connected undirected graph) [default: off]
    --fuzzy-dist <FLOAT>         Distance deviation (Å) for fuzzy seed expansion.
                                 Unions structure IDs from ±deviation neighbor hashes.
                                 0 = exact matching only. Paper default: 0.5 [default: 0]
    --fuzzy-angle <FLOAT>        Angle deviation (degrees) for fuzzy seed expansion.
                                 0 = exact matching only. Paper default: 5.0 [default: 0]
    --merge-iso <FLOAT>          Merge motifs whose structure-support sets have Jaccard
                                 similarity ≥ FLOAT. Reduces duplicate/isomorphic entries
                                 (e.g., same geometric contact encoded from opposite edge
                                 directions). Recommended: 0.7. [default: off]

output:
    -o, --output <PATH>          Output TSV path (default: stdout)

general options:
    -t, --threads <INT>          Number of threads [default: 1]
    -v, --verbose                Print verbose messages
    -h, --help                   Print this help menu

output columns:
    motif_id      sequential integer (sorted by adj_idf descending)
    num_residues  number of abstract residue positions in the motif
    num_edges     number of pairwise geometric constraints
    support       fraction of database structures containing the motif
    count         absolute count of matching structures
    motif_idf     sum of per-edge IDF scores: Σ log2(N / hash_count)
    mean_nres     mean residue count of matching structures (estimated)
    adj_idf       length-adjusted IDF: motif_idf × (mean_nres + 1)^(-0.5)
    edges         semicolon-separated edge descriptions:
                  nodeA-nodeB:hash_hex:AA1/AA2/ca_dist/cb_dist/angle[/phi1/phi2]
    top_structures top-5 example structures with residue positions:
                  name:res0-res1-res2-... (PDB serial numbers per motif node)
                  bare name shown when structure file is unavailable
";

pub fn mine_motifs(env: AppArgs) {
    match env {
        AppArgs::Mine {
            index_path,
            min_support,
            max_freq,
            max_residues,
            max_seeds,
            max_results,
            min_idf,
            require_complete,
            fuzzy_dist,
            fuzzy_angle,
            merge_iso_threshold,
            output,
            threads,
            verbose,
            help: _,
        } => {
            if verbose {
                print_logo();
            }

            let index_path = match index_path {
                Some(p) => p,
                None => {
                    eprintln!("{}", HELP_MINE);
                    std::process::exit(1);
                }
            };

            if verbose {
                print_log_msg(
                    INFO,
                    &format!(
                        "Foldmine: mining frequent motifs in index {} with {threads} threads",
                        index_path
                    ),
                );
                print_log_msg(
                    INFO,
                    &format!(
                        "Parameters: min_support={min_support:.4}, max_freq={max_freq:.4}, max_residues={max_residues}"
                    ),
                );
            }

            // Set up thread pool (matches pattern in analyze.rs).
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap_or_else(|e| {
                    eprintln!("Warning: failed to set thread pool: {e}");
                });

            // Load index.
            if verbose {
                print_log_msg(INFO, &format!("Loading index from {}", index_path));
            }
            let (index, _offset_mmap) = load_folddisco_index(&index_path);

            // Load lookup and configuration.
            let (lookup_path, config_path) = get_lookup_and_type(&index_path);
            let config = read_index_config_from_file(&config_path);
            let lookup = load_lookup_from_file(&lookup_path);

            if verbose {
                print_log_msg(
                    INFO,
                    &format!(
                        "Database: {} structures, hash type: {}",
                        lookup.len(),
                        config.hash_type.to_string()
                    ),
                );
            }

            // Build mining config.
            let max_nodes = max_residues as u8;
            let mining_config = MiningConfig {
                min_support,
                max_freq_ratio: max_freq,
                max_nodes,
                // Default: all edges in an undirected complete graph on max_nodes.
                max_edges: 0,
                max_seeds,
                max_results,
                threads,
                min_idf,
                require_complete,
                fuzzy_dist,
                fuzzy_angle,
                merge_iso_threshold,
                hash_type: config.hash_type,
                num_bin_dist: config.num_bin_dist,
                num_bin_angle: config.num_bin_angle,
            };

            // Run the algorithm.
            if verbose {
                print_log_msg(INFO, "Starting Foldmine...");
            }
            let mut results = mine_frequent_motifs(&index, &lookup, &mining_config);

            if verbose {
                print_log_msg(
                    INFO,
                    &format!("Found {} frequent motifs", results.len()),
                );
                print_log_msg(INFO, "Annotating top structures with residue positions...");
            }

            annotate_top_structures(
                &mut results,
                &lookup,
                config.hash_type,
                config.num_bin_dist,
                config.num_bin_angle,
            );

            // Write output.
            match output {
                Some(path) => {
                    let file = File::create(&path).unwrap_or_else(|e| {
                        eprintln!("Error creating output file '{}': {}", path, e);
                        std::process::exit(1);
                    });
                    let mut writer = BufWriter::new(file);
                    write_results(
                        &results,
                        &mut writer,
                        config.hash_type,
                        config.num_bin_dist,
                        config.num_bin_angle,
                        &lookup,
                    )
                    .unwrap_or_else(|e| {
                        eprintln!("Error writing results: {}", e);
                        std::process::exit(1);
                    });
                    if verbose {
                        print_log_msg(INFO, &format!("Results written to {}", path));
                    }
                }
                None => {
                    let stdout = std::io::stdout();
                    let mut writer = BufWriter::new(stdout.lock());
                    write_results(
                        &results,
                        &mut writer,
                        config.hash_type,
                        config.num_bin_dist,
                        config.num_bin_angle,
                        &lookup,
                    )
                    .unwrap_or_else(|e| {
                        eprintln!("Error writing results: {}", e);
                        std::process::exit(1);
                    });
                }
            }
        }
        _ => {
            eprintln!("{}", HELP_MINE);
            std::process::exit(1);
        }
    }
}
