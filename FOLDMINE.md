# Foldmine

**Foldmine** is a frequent structural motif discovery algorithm built on top of the Folddisco index. Instead of searching for a *known* query motif, Foldmine inverts the problem: it mines *all* structural motifs that appear in at least a user-specified fraction of the database structures.

The algorithm is analogous to Apriori / gSpan frequent-subgraph mining, but operates entirely in Folddisco's geometric hash space, making it efficient over large protein databases.

---

## Algorithm overview

Foldmine works in three phases:

1. **Seed collection** — Scan the inverted index and collect all geometric hashes whose document frequency falls in `[min_support × N, max_freq × N]`.  Hashes below the minimum are infrequent by the anti-monotone property; hashes above the maximum encode generic structural features (e.g., helix–helix contacts) and are excluded.

2. **DFS motif growth** — Starting from each frequent single-edge seed, build abstract motif graphs by adding edges one at a time. A motif is a connected graph whose nodes are abstract residue positions (0, 1, 2, …) and whose edges carry a geometric hash encoding the pairwise geometry. The support of a multi-edge motif is bounded above by the sorted intersection of all per-edge structure-ID lists. Any motif whose intersection falls below `min_support × N` is pruned along with all its extensions (anti-monotone pruning).

3. **Canonical deduplication** — Every motif is reduced to its lexicographically smallest edge list over all node-label permutations. A concurrent hash map prevents the same canonical motif from being reported more than once across parallel DFS branches.

---

## Usage

```
folddisco foldmine -i <index_prefix> [options]
```

### Required

| Flag | Description |
|------|-------------|
| `-i, --index <PREFIX>` | Path prefix of the Folddisco index |

### Mining options

| Flag | Default | Description |
|------|---------|-------------|
| `--min-support <FLOAT>` | `0.01` | Minimum support fraction in [0, 1].  A motif must appear in at least this fraction of the database structures. |
| `--max-freq <FLOAT>` | `0.5` | Maximum frequency ratio in [0, 1].  Hashes present in more than this fraction are excluded (too generic). |
| `--max-residues <INT>` | `6` | Maximum number of abstract residue positions per motif. |
| `--max-seeds <INT>` | `0` (unlimited) | Limit the number of seed edges explored (rarest first). Useful for very large or dense indices where full exploration is expensive. |
| `--max-results <INT>` | `0` (unlimited) | Stop after discovering this many motifs. |

### Output options

| Flag | Default | Description |
|------|---------|-------------|
| `-o, --output <PATH>` | stdout | Write TSV output to file instead of stdout. |

### General options

| Flag | Default | Description |
|------|---------|-------------|
| `-t, --threads <INT>` | `1` | Number of parallel threads (seed-level parallelism). |
| `-v, --verbose` | off | Print progress messages to stderr. |
| `-h, --help` | — | Print help. |

---

## Output format

Tab-separated with header:

```
motif_id  num_residues  num_edges  support  count  edges  top_structures
```

| Column | Description |
|--------|-------------|
| `motif_id` | Sequential integer (sorted by support descending, then edges ascending) |
| `num_residues` | Number of abstract residue positions in the motif |
| `num_edges` | Number of pairwise geometric constraints |
| `support` | Fraction of database structures containing the motif |
| `count` | Absolute count of matching structures |
| `edges` | Semicolon-separated edge descriptions: `nodeA-nodeB:hash_hex:AA1/AA2/ca_dist_Å/cb_dist_Å/angle°[/phi1°/phi2°]` |
| `top_structures` | Comma-separated names of up to 5 representative matching structures |

### Example row (2-residue, 1-edge motif)

```
0	2	1	0.4000	2	0-1:00013d8b:SER/HIS/6.82/8.23/-113.2/-101.3/-144.5	data/sp/1pq5.pdb,data/sp/4cha.pdb
```

---

## Compatibility with Folddisco indices

Foldmine reads the same index format used by `folddisco index` and `folddisco query`. **All pre-built indices** listed in the Folddisco README work directly:

| Database | URL |
|----------|-----|
| Human proteome | https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4 |
| *E. coli* proteome | https://opendata.mmseqs.org/folddisco/e_coli_folddisco.tar.lz4 |
| AFDB 16 model organisms | https://opendata.mmseqs.org/folddisco/afdb_proteome_v4_folddisco.tar.lz4 |
| Swiss-Prot | https://opendata.mmseqs.org/folddisco/afdb_swissprot_v4_folddisco.tar.lz4 |
| AFDB50 | https://opendata.mmseqs.org/folddisco/afdb50_v4_folddisco.tar.lz4 |
| ESM30 | https://opendata.mmseqs.org/folddisco/highquality_clust30_folddisco.tar.lz4 |

The index prefix is the filename without any extension (e.g., `index/h_sapiens_folddisco`). The tool expects three sidecar files to be present alongside the main index:

```
<prefix>          # main value file (structure-ID lists)
<prefix>.offset   # hash → byte-offset table
<prefix>.lookup   # structure name / metadata lookup
<prefix>.type     # hash type and bin-count configuration
```

---

## Examples

### Small custom index

Build an index from local PDB files and run Foldmine:

```bash
# Index a set of serine peptidase structures
folddisco index -p data/serine_peptidases -i index/serine_peptidases_folddisco

# Mine motifs present in ≥ 40 % of the structures, up to 4 residues
folddisco foldmine -i index/serine_peptidases_folddisco \
    --min-support 0.4 --max-freq 0.9 \
    --max-residues 4 \
    -o motifs.tsv -v
```

### Pre-built Swiss-Prot index

```bash
# Download and extract
mkdir -p index && cd index
aria2c https://opendata.mmseqs.org/folddisco/afdb_swissprot_v4_folddisco.tar.lz4
lz4 -dc afdb_swissprot_v4_folddisco.tar.lz4 | tar -xvf -
cd ..

# Mine motifs present in ≥ 1 % of Swiss-Prot structures, ≤ 3 residues,
# using 8 threads, limited to top 10 000 motifs
folddisco foldmine -i index/afdb_swissprot_v4_folddisco \
    --min-support 0.01 --max-freq 0.5 \
    --max-residues 3 \
    --max-seeds 5000 \
    --max-results 10000 \
    -t 8 \
    -o swissprot_motifs.tsv -v
```

### AFDB50 (very large index)

For databases with millions of structures, use aggressive limits to keep runtime manageable:

```bash
folddisco foldmine -i index/afdb50_v4_folddisco \
    --min-support 0.001 --max-freq 0.3 \
    --max-residues 3 \
    --max-seeds 2000 \
    --max-results 50000 \
    -t 16 \
    -o afdb50_motifs.tsv -v
```

---

## Parameter tuning guide

| Scenario | Recommendation |
|----------|---------------|
| Small index (< 1 000 structures) | Increase `--min-support` (e.g., 0.3–0.5) to avoid noise; omit `--max-seeds`/`--max-results` |
| Large index (> 100 000 structures) | Use `--max-seeds 2000–5000` and `--max-results 10000–50000` to control runtime |
| Catalytic / binding site discovery | `--max-residues 3–4`, `--min-support 0.001–0.01`, `--max-freq 0.3` |
| Scaffold / fold discovery | `--max-residues 4–6`, `--min-support 0.01–0.05`, `--max-freq 0.5` |
| Fast exploration | `--max-residues 2`, no limits on seeds/results |

**Note on `min_count > max_count` errors**: if `--min-support` and `--max-freq` are set such that the required minimum count exceeds the maximum allowed count (possible on tiny datasets), the tool will warn and exit cleanly. Increase `--max-freq` or decrease `--min-support`.

---

## Algorithm details

### Anti-monotone pruning

If a single-edge motif with hash *h* appears in fewer than `min_support × N` structures, then no multi-edge motif containing *h* can be frequent either. This is the core efficiency guarantee (Apriori principle applied to geometric hash intersections).

### Canonical form

Two abstract motif graphs that differ only by relabeling their residue positions (nodes) represent the same structural pattern. Foldmine normalizes every motif by trying all *n*! node-label permutations (at most 6! = 720 for `--max-residues 6`) and keeping the lexicographically smallest edge list. A concurrent DashMap tracks seen canonical forms so each unique motif is reported exactly once.

### Parallelism

- **Seed-level**: each frequent seed edge launches an independent DFS branch (`rayon::par_iter`).
- **Shared state**: the `seen` DashMap and result vector are shared across threads with lock-free / mutex-based access respectively.
- The `-t` flag controls the rayon thread pool size.
