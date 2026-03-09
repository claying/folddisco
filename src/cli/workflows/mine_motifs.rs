// File: mine_motifs.rs
// Foldmine workflow: mine frequent structural motifs from an index.

use std::fs::File;
use std::io::BufWriter;

use crate::cli::{AppArgs, print_logo};
use crate::cli::config::read_index_config_from_file;
use crate::controller::frequent::{mine_frequent_motifs, write_results, MiningConfig};
use crate::controller::io::get_lookup_and_type;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::{print_log_msg, INFO};

pub const HELP_MINE: &str = "\
usage: folddisco foldmine -i <index_prefix> [options]

Foldmine: discover frequent structural motifs in a Folddisco index.

input:
    -i, --index <INDEX_PREFIX>   Path prefix to index files (required)

mining options:
    --min-support <FLOAT>        Minimum support fraction in [0,1] [default: 0.01]
    --max-freq <FLOAT>           Maximum frequency ratio in [0,1] [default: 0.5]
    --max-residues <INT>         Maximum residues per motif [default: 6]
    --max-seeds <INT>            Limit number of seed edges considered (0=unlimited) [default: 0]
    --max-results <INT>          Stop after finding this many motifs (0=unlimited) [default: 0]

output:
    -o, --output <PATH>          Output TSV path (default: stdout)

general options:
    -t, --threads <INT>          Number of threads [default: 1]
    -v, --verbose                Print verbose messages
    -h, --help                   Print this help menu

output columns:
    motif_id      sequential integer
    num_residues  number of abstract residue positions in the motif
    num_edges     number of pairwise geometric constraints
    support       fraction of database structures containing the motif
    count         absolute count of matching structures
    edges         semicolon-separated edge descriptions:
                  nodeA-nodeB:hash_hex:AA1/AA2/ca_dist/cb_dist/angle[/phi1/phi2]
    top_structures top-5 example structure names
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
            };

            // Run the algorithm.
            if verbose {
                print_log_msg(INFO, "Starting Foldmine...");
            }
            let results = mine_frequent_motifs(&index, &lookup, &mining_config);

            if verbose {
                print_log_msg(
                    INFO,
                    &format!("Found {} frequent motifs", results.len()),
                );
            }

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
