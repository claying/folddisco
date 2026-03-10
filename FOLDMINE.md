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
folddisco mine -i <index_prefix> [options]
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
| `--min-idf <FLOAT>` | `0.0` (no filter) | Minimum length-adjusted IDF score (`adj_idf`).  Filters out generic/trivial motifs such as alpha-helix contacts. Typical values: 0.5–2.0. |
| `--require-complete` | off | Only output motifs where **every unordered node pair has at least one directed edge**.  Ensures all pairwise geometric relationships are captured.  Intermediate partial graphs are still explored internally for pruning efficiency, but not reported. |
| `--fuzzy-dist <FLOAT>` | `0.0` (exact) | Distance deviation (Å) for fuzzy seed expansion.  Each seed hash is expanded by ±deviation on each distance feature; structure IDs from all neighbor hashes are unioned.  Improves recall for structures with minor conformational variation.  Paper default: **0.5 Å**. |
| `--fuzzy-angle <FLOAT>` | `0.0` (exact) | Angle deviation (degrees) for fuzzy seed expansion.  Applied independently to each angle feature.  Paper default: **5.0°**.  Note: internally converted to radians for sin/cos hash types (PDBTrRosetta, TrRosetta, PDBMotifSinCos). |

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
motif_id  num_residues  num_edges  support  count  motif_idf  mean_nres  adj_idf  edges  top_structures
```

| Column | Description |
|--------|-------------|
| `motif_id` | Sequential integer (sorted by `adj_idf` descending, then support descending, then edges ascending) |
| `num_residues` | Number of abstract residue positions in the motif |
| `num_edges` | Number of pairwise geometric constraints |
| `support` | Fraction of database structures containing the motif |
| `count` | Absolute count of matching structures |
| `motif_idf` | Inverse document frequency: Σ log₂(N / hash\_count) over all edges.  Higher = rarer hash combination. |
| `mean_nres` | Mean residue count of matching structures (estimated from up to 5 representative hits). |
| `adj_idf` | Length-adjusted IDF: `motif_idf × (mean_nres + 1)^(−0.5)`.  Penalizes motifs found only in large multi-domain proteins.  Use this column to rank and filter. |
| `edges` | Semicolon-separated edge descriptions: `nodeA-nodeB:hash_hex:AA1/AA2/ca_dist_Å/cb_dist_Å/angle°[/phi1°/phi2°]` |
| `top_structures` | Comma-separated names of up to 5 representative matching structures, with residue positions when available (see below). |

### `top_structures` column format

When the original structure files are accessible at the paths stored in the index, Foldmine annotates each representative structure with the actual PDB residue serial numbers of the residues that realise the motif:

```
name:res0-res1-res2-...
```

- `name` — structure name (as stored in the index lookup table, typically the file path or PDB/AF accession)
- `res0`, `res1`, … — PDB serial numbers of the residue assigned to abstract motif node 0, 1, … respectively

Multiple structures are comma-separated:

```
1abc.pdb:57-102-195,2xyz.pdb:34-89-176,3def.pdb:61-108-201
```

**When structure files are not accessible** (e.g., the index was built on a different machine or the files have been moved), only the bare name is shown:

```
1abc.pdb,2xyz.pdb,3def.pdb
```

No error is printed in this case — residue annotation silently falls back to bare names.

> **Note:** Residue annotation requires the original structure files (`.pdb`, `.cif`, or `.pdb.gz` / `.cif.gz`) to be readable at the paths recorded in the index.  The index itself stores only structure IDs — coordinates are not stored and cannot be recovered from the index alone.

### Example row (2-residue, 1-edge motif, with residue annotation)

```
0	2	1	0.4000	2	8.3214	312.5	0.4706	0-1:00013d8b:SER/HIS/6.82/8.23/-113.2/-101.3/-144.5	data/sp/1pq5.pdb:57-102,data/sp/4cha.pdb:195-40
```

### Example row (bare names, no structure files available)

```
0	2	1	0.4000	2	8.3214	312.5	0.4706	0-1:00013d8b:SER/HIS/6.82/8.23/-113.2/-101.3/-144.5	1pq5.pdb,4cha.pdb
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

> **Pre-built indices and residue annotation:** Pre-built indices reference the AlphaFold CIF files that were present on the server where the index was built.  Those paths are not available on your local machine, so `top_structures` will show bare accession names without residue positions.  To get residue annotations, build a local index with `folddisco index` pointing to your own copy of the structure files.

---

## Examples

### Small custom index

Build an index from local PDB files and run Foldmine:

```bash
# Index a set of serine peptidase structures
folddisco index -p data/serine_peptidases -i index/serine_peptidases_folddisco

# Mine motifs present in >= 40% of the structures, up to 4 residues
# --require-complete ensures all pairwise relationships are captured
folddisco mine -i index/serine_peptidases_folddisco \
    --min-support 0.4 --max-freq 0.9 \
    --max-residues 4 \
    --require-complete \
    -o motifs.tsv -v

# Same search with fuzzy matching (paper defaults)
folddisco mine -i index/serine_peptidases_folddisco \
    --min-support 0.4 --max-freq 0.9 \
    --max-residues 4 \
    --require-complete \
    --fuzzy-dist 0.5 --fuzzy-angle 5.0 \
    -o motifs_fuzzy.tsv -v
```

### Pre-built Swiss-Prot index

```bash
# Download and extract
mkdir -p index && cd index
aria2c https://opendata.mmseqs.org/folddisco/afdb_swissprot_v4_folddisco.tar.lz4
lz4 -dc afdb_swissprot_v4_folddisco.tar.lz4 | tar -xvf -
cd ..

# Mine motifs present in >= 1% of Swiss-Prot structures, <= 3 residues,
# using 8 threads, limited to top 10 000 motifs, filtering trivial helices/sheets
folddisco mine -i index/afdb_swissprot_v4_folddisco \
    --min-support 0.01 --max-freq 0.5 \
    --max-residues 3 \
    --max-seeds 5000 \
    --max-results 10000 \
    --min-idf 1.0 \
    -t 8 \
    -o swissprot_motifs.tsv -v
```

### AFDB50 (very large index)

For databases with millions of structures, use aggressive limits to keep runtime manageable:

```bash
folddisco mine -i index/afdb50_v4_folddisco \
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

### Fully-connected motifs (`--require-complete`)

By default, Foldmine reports every sub-graph discovered during DFS growth, including sparse intermediate forms (e.g., a 3-node motif with only 2 edges).  When `--require-complete` is on, a motif is only reported if **every unordered node pair {i, j} has at least one directed edge** (either i→j or j→i).

| Nodes | Minimum edges required | Example |
|-------|----------------------|---------|
| 2 | 1 | trivially satisfied by any seed |
| 3 | 3 | a directed triangle |
| 4 | 6 | all pairs covered |
| k | k×(k−1)/2 | |

The DFS growth and anti-monotone pruning continue as usual — partial motifs are still explored, just not reported.  This flag is most useful when you want all pairwise geometric constraints explicitly captured rather than relying on the absence of an edge to mean "unconstrained".

### Fuzzy hash matching (`--fuzzy-dist`, `--fuzzy-angle`)

Folddisco geometric hashes discretize distances and angles into bins.  A pair of residues whose geometry lies exactly on a bin boundary may be assigned to a different bin in a slightly different conformation, making the motif miss those structures under exact matching.

The fuzzy expansion (described in the paper as "Extended search") addresses this:

1. For each frequent seed hash *h*, decode it to continuous feature values.
2. Apply ±`fuzzy_dist` (Å) independently to each distance feature and ±`fuzzy_angle` (°) independently to each angle feature — up to ~10 neighbor hashes.
3. Look up the structure ID lists for all neighbor hashes in the index and **union** them with the exact-match list.
4. Use this expanded list as the seed's support set for all subsequent intersection operations.

This makes the effective support test "does the structure contain geometry approximately equal to this hash" rather than "exactly equal".

> **Support inflation warning**: fuzzy expansion increases support counts (more structures match each seed), so effective support fractions are higher than under exact matching.  When using `--fuzzy-dist`/`--fuzzy-angle`, consider raising `--min-support` by 2–3× to compensate.

### Edge uniqueness

Each ordered node pair (nodeA → nodeB) can appear at most once in a motif's edge list. A second edge between the same ordered pair (with a different hash) would represent an over-specified, artifact constraint and is suppressed during DFS growth.

### IDF scoring

After DFS, Foldmine computes three scores for each discovered motif:

| Score | Formula | Meaning |
|-------|---------|---------|
| `motif_idf` | Σ log₂(N / count_i) over all edges | Sum of per-edge inverse document frequencies.  Rarer hash combinations → higher score. |
| `mean_nres` | mean(residues of matching structures) | Average protein length among matching hits.  Estimated from up to 5 representative structures. |
| `adj_idf` | `motif_idf × (mean_nres + 1)^(−0.5)` | Length-penalized IDF.  Penalizes motifs found primarily in large multi-domain proteins (which incidentally contain many structural patterns). |

The output is sorted by `adj_idf` descending (highest specificity first).  Use `--min-idf` to hard-filter low-scoring, generic motifs (alpha-helix contacts typically score < 0.5; specific active-site geometry typically scores > 2.0).

### Residue annotation

After mining, Foldmine attempts to annotate each representative ("top") structure with the actual residue positions (PDB serial numbers) that realise the motif.  The process:

1. Load the structure file from the path stored in the index lookup table.
2. Scan all ordered residue pairs within 20 Å Cα–Cα distance and compute their geometric hashes.
3. Collect candidate (residue_i, residue_j) pairs for each motif edge hash.
4. Run a backtracking constraint solver (rarest edges first for pruning efficiency) to find a consistent, injective assignment of abstract motif nodes to actual residues.
5. Report PDB serial numbers for the assigned residues.

This step requires the original structure files to be accessible.  If they are not, annotation is silently skipped and only bare structure names are shown.

### Parallelism

- **Seed-level**: each frequent seed edge launches an independent DFS branch (`rayon::par_iter`).
- **Work-stealing DFS**: at shallow depths (< 2), candidates are spawned as rayon tasks so that unbalanced DFS trees distribute naturally across idle threads.
- **Shared state**: the `seen` DashMap and result vector are shared across threads with lock-free / mutex-based access respectively.
- The `-t` flag controls the rayon thread pool size.
