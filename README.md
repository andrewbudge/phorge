<div align="center">

<img src="https://raw.githubusercontent.com/andrewbudge/phorge/main/docs/phorge-assets/phorge-mark-ember.svg" alt="phorge logo" width="140">

# phorge

**A composable CLI toolkit for phylogenetics — from NCBI acquisition to supermatrix, in one binary.**

[![Crates.io](https://img.shields.io/crates/v/phorge.svg)](https://crates.io/crates/phorge)
[![Downloads](https://img.shields.io/crates/d/phorge.svg)](https://crates.io/crates/phorge)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![Rust](https://img.shields.io/badge/rust-2024_edition-orange.svg)

</div>

Phorge is a lightweight, composable CLI phylogenetics toolkit. A single binary, Phorge provides many subcommands that replace common chains of bash commands or a collection of individual programs in phylogenetic pipelines. Examples include NCBI sequence acquisition, homology-based gene extraction, concatenation, and alignment quality control.

Phorge has two layers in one binary: lean file tools (`getheaders`, `concat`, `stats`, `coverage`, `convert`, `filter`, `curate`, `align`) and an acquisition layer (`query`, `fetch`, `extract`, `clean`) that pulls and curates sequences from NCBI. A broken external tool or network never affects the self-contained subcommands.

> **Note:** Phorge is under active development. Subcommands may change or be added as the project matures.

## Contents

- [Install](#install)
- [External Dependencies](#external-dependencies)
- **Subcommands**
  - File tools: [getheaders](#getheaders-ghd) · [concat](#concat-liger) · [stats](#stats) · [coverage](#coverage) · [convert](#convert) · [filter](#filter) · [curate](#curate) · [align](#align-aln)
  - Acquisition: [query](#query) · [fetch](#fetch) · [extract](#extract) · [clean](#clean)
- [Pipeline example](#pipeline-example)
- [Planned Subcommands](#planned-subcommands)

## Install

Requires [Rust](https://www.rust-lang.org/tools/install).

```bash
cargo install phorge
```

This builds the binary and adds `phorge` to your PATH. Running it again later upgrades to the latest published version.

To build from the latest (unreleased) source instead:

```bash
cargo install --git https://github.com/andrewbudge/phorge
```

Updating a git install requires `--force` (add `--force` to either command to reinstall the same version).

## External Dependencies

Most phorge subcommands are fully self-contained. The `extract` and `align` subcommands orchestrate external tools and require them to be in your PATH.

| Subcommand | Requires |
|---|---|
| `extract` | [MMseqs2](https://github.com/soedinglab/MMseqs2) |
| `align` | [MAFFT](https://mafft.cbrc.jp/alignment/software/) or [MUSCLE](https://drive5.com/muscle/) |

If you don't have these installed, the easiest path is conda:

```bash
conda env create -f environment.yml
conda activate phorge-tools
```

Or install individually:

```bash
conda install -c bioconda mmseqs2 mafft muscle
```

The `query` and `fetch` subcommands need an internet connection and a valid email address (required by NCBI's Terms of Service for automated E-utilities access), but no external binary.

All other subcommands (`getheaders`, `concat`, `stats`, `coverage`, `convert`, `filter`, `curate`, `clean`) have no external dependencies.

## Subcommands
---
### getheaders (ghd)

Extract headers from FASTA files.

**Example:**

```bash
$ phorge getheaders testdata/test_good.fasta
Cat
Dog
Cat

$ phorge getheaders -u testdata/test_good.fasta
Cat
Dog
```
---
### concat (liger)

Concatenate multiple gene alignments into a supermatrix. Unlike other tools, input files can live anywhere and globs are accepted.

Concat runs in two modes:

- **Exact match (default):** headers must match exactly across files, like FASconCAT and AMAS.
- **Smart match (`-a alias.txt`):** pass an alias list — a file of clean output names (one per line, e.g. `Mus_musculus`) that get matched to messy input headers via case-insensitive substring search. Longer aliases match first to prevent partial collisions (e.g. `Mus musculus domesticus` claims before `Mus musculus`). Once a header is claimed it cannot be matched again. The alias list doubles as a rename map — input headers stay messy, output gets clean names. Requires `-l` for a provenance TSV that records exactly which original header matched each alias.

Concat auto-detects DNA vs amino acid data per gene and adjusts missing characters and partition labels accordingly. FASTA output goes to stdout, partition boundaries to stderr in RAxML/IQ-TREE format by default. NEXUS bundles everything into one file.

**Exact match — clean headers:**

```bash
$ phorge concat gene1.fasta gene2.fasta > supermatrix.fasta
DNA, gene1.fasta = 1-4
DNA, gene2.fasta = 5-8
```

**Smart match — messy headers with an alias list:**

```bash
$ cat alias.txt
Mus_musculus
Rattus_rattus
Xenopus_laevis

$ phorge concat -a alias.txt -l prov.tsv gene1.fasta gene2.fasta > supermatrix.fasta
DNA, gene1.fasta = 1-4
DNA, gene2.fasta = 5-8

$ cat supermatrix.fasta
>Mus_musculus
ATCGATCG
>Rattus_rattus
ATCGNNNN
>Xenopus_laevis
NNNNATCG

$ cat prov.tsv
alias.txt	gene1.fasta	gene2.fasta
Mus_musculus	AB123.1 Mus musculus gene1 cds	XM456.1 Mus musculus gene2 cds
Rattus_rattus	AB124.1 Rattus rattus gene1 cds	MISSING
Xenopus_laevis	MISSING	XM789.1 Xenopus laevis gene2 cds
```

**NEXUS output:**

```bash
$ phorge concat -a alias.txt -l prov.tsv -f nexus gene1.fasta gene2.fasta
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=N GAP=-;
  MATRIX
  Mus_musculus    ATCGATCG
  Rattus_rattus    ATCGNNNN
  Xenopus_laevis    NNNNATCG
;
END;
BEGIN SETS;
  CHARSET gene1.fasta = 1-4;
  CHARSET gene2.fasta = 5-8;
END;
```

**Flags:**
- `-a, --alias` — alias list for smart matching (clean output names that map to messy input headers)
- `-l, --log` — provenance TSV output file (required with `-a`)
- `-f, --format` — output format: fasta (default), nexus (also accepts `n` or `nex`)
- `-m, --missing` — override missing data character (default: auto per data type — N for DNA, X for amino acid, ? for mixed)
- `-p, --partitions` — partition format: raxml (default, also used by IQ-TREE) or nexus
- `--dry-run` — show a matching summary (per-gene match counts and per-taxon coverage) without building the supermatrix
---
### stats

Get basic alignment statistics from FASTA files. Accepts multiple files via globs. Automatically detects DNA vs amino acid sequences.

**Columns:**
- **file** — filename (path stripped)
- **sequences** — number of sequences
- **length** — alignment length (NA if unaligned)
- **type** — `DNA` or `AA` (auto-detected, supports IUPAC ambiguity codes)
- **gc_pct** — GC content as a percentage of real bases (NA for amino acid data)
- **missing_pct** — percentage of gaps and unknown characters
- **variable** — sites with at least 2 different residues (excluding gaps/unknowns)
- **variable_pct** — variable sites as a percentage of alignment length
- **informative** — parsimony-informative sites (at least 2 residues each appearing 2+ times)
- **informative_pct** — informative sites as a percentage of alignment length

**Example:**

```bash
$ phorge stats supermatrix.fasta proteins.fasta
file	sequences	length	type	gc_pct	missing_pct	variable	variable_pct	informative	informative_pct
supermatrix.fasta	3	8	DNA	50.0	33.3	0	0.0	0	0.0
proteins.fasta	4	20	AA	NA	0.0	3	15.0	2	10.0
```

**Flags:**
- `-d, --detailed` — per-sequence statistics (header, length, GC%, missingness)
- `-p, --pretty` — column-aligned output for readability
---
### coverage

Summarize taxa and loci coverage from a concat provenance TSV. Shows how many loci each taxon appears in, or how many taxa each locus has.

**Example:**

```bash
$ phorge coverage -t prov.tsv
taxa	loci_present	loci_missing	pct_missing
Mus_musculus	5/5	0/5	0.0%
Smilodon_populator	2/5	3/5	60.0%

$ phorge coverage -l -p prov.tsv
loci          appearance_count  missing_pct
12S_aln.fas   6/8               25.0%
COX1_aln.fas  6/8               25.0%
```

**Flags:**
- `-t, --taxa` — show per-taxon coverage (how many loci each taxon has)
- `-l, --loci` — show per-loci coverage (how many taxa each locus has)
- `-p, --pretty` — column-aligned output for readability
---
### convert

Convert between common sequence data file types. Auto-detects the input format from file contents.

**Supported formats:**
- FASTA (`f`)
- NEXUS (`n` / `nex` / `nexus`)
- Relaxed PHYLIP (`rp` / `phylip`)
- Strict PHYLIP (`sp`)

**Example:**

```bash
$ phorge convert -o n alignment.fasta
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=N GAP=-;
  MATRIX
  Taxon_A    ATCGATCG
  Taxon_B    ATCGATCG
  Taxon_C    ATCGNNNN
;
END;

$ phorge convert -o rp alignment.nex
3 8
Taxon_A    ATCGATCG
Taxon_B    ATCGATCG
Taxon_C    ATCGNNNN
```

**Flags:**
- `-o, --output_format` — output format: `f` (fasta), `n` (nexus), `rp` (relaxed phylip), `sp` (strict phylip)
---
### query

Search NCBI's `nuccore` database for one or more taxa and write a `query_results.json` manifest — the metadata spine that `fetch` and `clean` read. No sequences are downloaded at this stage; this only collects accessions and their TaxID/name/length. Each TaxID is expanded to its full subtree (`txidNNN[Organism:exp]`), excluding environmental samples.

Requires an internet connection and an email address (NCBI Terms of Service).

**Example:**

```bash
$ phorge query --ingroup 89829 --outgroup 241031 309676 -o run/ --email you@example.org
querying nuccore for txid89829[Organism:exp]
Leptophlebiidae (89829): 3437 records found; retrieving metadata...
...
query complete
  total accessions: 3586
  written to:       run/query_results.json
```

**Flags:**
- `--ingroup` — one or more ingroup TaxIDs (required)
- `--outgroup` — one or more outgroup TaxIDs
- `-o, --out` — output directory (writes `query_results.json`)
- `--email` — email address required by NCBI ToS (required)
- `--api-key` — NCBI API key (optional; raises the rate limit from 3 to 10 req/s)
---
### fetch

Download the sequences for a `query_results.json` manifest, writing raw NCBI FASTA shards to `<out>/raw/`. The download is resumable — a manifest tracks completed shards, so an interrupted run picks up where it left off. Headers are written verbatim; rewriting them is `clean`'s job.

```bash
$ phorge fetch -q run/query_results.json -o run/ --email you@example.org --yes
preflight ready to download  records=3586  chunks=8  est_mb=2.4
shard written  chunk=0  records=500
...
fetch complete  records=3586  shard_dir=run/raw
```

**Flags:**
- `-q, --query` — path to `query_results.json` (from `query`)
- `-o, --out` — output directory (shards go to `<out>/raw/`)
- `--min-length` / `--max-length` — drop records outside a length range before downloading
- `--email` — email address required by NCBI ToS (required)
- `--api-key` — NCBI API key (optional)
- `--yes` — skip the download-size confirmation prompt (for non-interactive use)
---
### extract

Extract gene regions from target organism sequences using homology search. Takes reference gene sequences and one or more target FASTAs (or a directory), runs MMseqs2 `easy-search`, and writes one output FASTA per gene containing the extracted region from each organism that had a hit. The extracted hit region is cut at the MMseqs2 coordinates; the original target header is preserved so downstream tools (`clean`) can recover the accession.

References come in two forms (one is required):
- `-r, --reference` — a single FASTA where each record header is the gene name (`>COX1`, `>ND2`). Convenient for ad-hoc use.
- `--refs` — one FASTA per gene, where the filename stem is the gene name (`COI.fasta` → COI). Each file may hold several sequences to cover divergence across taxa. This is the pipeline form.

Requires [MMseqs2](https://github.com/soedinglab/MMseqs2) installed and in your PATH.

**Example:**

```bash
# refs/ has one file per gene: COI.fasta, 16S.fasta, 28S.fasta, ...
# run/raw/ contains the FASTA shards written by fetch

$ phorge extract --refs refs/*.fasta -t run/raw/ -o run/genes/
Pooled 19 reference sequence(s).
Pooling 8 target files...
Parsing results...
Done. Extracted 7 gene(s) from 3069 hits.

$ ls run/genes/
12S.fasta  16S.fasta  18S.fasta  28S.fasta  COI.fasta  cytb.fasta  H3.fasta
```

**Flags:**
- `-r, --reference` — single reference FASTA, gene name = each record header (`>COX1`)
- `--refs` — per-gene reference FASTAs, gene name = filename stem (`COI.fasta` → COI)
- `-t, --targets` — target organism FASTA files or a directory containing them
- `-o, --output` — output directory for per-gene FASTAs
- `--min-identity` — minimum MMseqs2 sequence identity to keep a hit, 0.0–1.0 (default: 0.7); the sole quality gate, so choose references that cover your taxa
- `-s, --sensitivity` — MMseqs2 sensitivity, 1.0 (fast) to 7.5 (max); default 5.7
- `--flank` — extra bases to grab on either side of each hit (default: 0)
- `--keep-intermediates` — keep the temp directory with pooled targets and raw MMseqs2 output
---
### clean

Join `extract`'s per-gene output back to `query_results.json`, rewrite headers to `TaxID|Name|Accession|Gene`, and deduplicate to one sequence per taxon per gene. This recovers the TaxID and clean taxon name that homology search alone doesn't carry, and collapses the many accessions NCBI holds per taxon down to a single best representative per gene.

Dedup keeps the longest sequence per TaxID, breaking ties by extract identity. Records whose accession isn't found in `query_results.json` are dropped and reported (broken provenance is useless to `concat`).

```bash
$ phorge clean --genes-dir run/genes/ -q run/query_results.json -o run/clean/
Done. Wrote 591 cleaned sequence(s) across 7 gene file(s); dropped 2478 duplicate(s).
```

Use `--prefer` to favour particular records during dedup — for example your own museum vouchers — even when they aren't the longest. A record is preferred if the substring appears in its extract header or its GenBank title:

```bash
$ phorge clean --genes-dir run/genes/ -q run/query_results.json -o run/clean/ --prefer MyLab
Done. Wrote 591 cleaned sequence(s) across 7 gene file(s); dropped 2478 duplicate(s).
  11 kept record(s) matched --prefer ["MyLab"].
```

**Flags:**
- `--genes-dir` — directory of per-gene FASTAs from `extract`
- `-q, --query` — `query_results.json` (the accession → TaxID/name table)
- `-o, --out` — output directory
- `--prefer` — prefer records whose extract header or GenBank title contains this substring during dedup; repeatable, and overrides the longest-sequence rule
---
### align (aln)

Batch align multiple FASTA files using MAFFT or MUSCLE. Runs the aligner on each input file and writes output to a directory with a consistent naming convention. Aligner stderr is captured to `align.log` in the output directory.

**Example:**

```bash
$ phorge align -p mafft -i genes/*.fasta -e _aln -o aligned/
Aligning COI...done
Aligning ND2...done
Aligning 12S...done
Done. Aligned 3 files.
```

Pass custom flags to the aligner after `--` (replaces the default flag):

```bash
# mafft — replace --auto
$ phorge align -p mafft -i genes/*.fasta -e _aln -o aligned/ -- --thread 4 --maxiterate 1000

# muscle — replace -align with -super5 for large datasets
$ phorge align -p muscle -i genes/*.fasta -e _aln -o aligned/ -- -super5
```

**Flags:**
- `-p, --program` — alignment program: `mafft` or `muscle` (name or full path)
- `-i, --input` — input unaligned FASTA files (glob or list)
- `-e, --extension` — suffix to append to output filenames (default: `_aln`)
- `-o, --output` — output directory for aligned files
- `--` — extra flags passed verbatim to the aligner; replaces the default (`--auto` for mafft, `-align` for muscle)
---
### filter

Remove taxa from an alignment that exceed a missingness threshold, have too few loci in a supermatrix, or both. Filters can be used independently or combined — a taxon must pass all applied filters to be kept. Output goes to stdout, summary to stderr.

**Example:**

```bash
# drop taxa with more than 50% gaps in the supermatrix
$ phorge filter supermatrix.fasta --max-missing 0.5 > filtered.fasta
Total taxa: 8
Kept taxa: 6
Dropped taxa: 2

# drop taxa present in fewer than 3 loci (requires the provenance TSV from concat -l)
$ phorge filter supermatrix.fasta --min-loci 3 -l prov.tsv > filtered.fasta

# both filters at once
$ phorge filter supermatrix.fasta --max-missing 0.5 --min-loci 3 -l prov.tsv > filtered.fasta
```

**Flags:**
- `-m, --max-missing` — maximum allowed missingness fraction per taxon (0.0–1.0)
- `-n, --min-loci` — minimum number of loci a taxon must be present in
- `-l, --log` — provenance TSV from `phorge concat -l` (required with `--min-loci`)
---
### curate

Trim alignment columns by parsimony-informativeness and gappiness. A native Rust port of [ClipKIT](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007) (Steenwyk et al. 2020). Accepts multiple files and a glob. Output goes to a directory or stdout, summary to stderr.

Compose keep conditions with `-k` by combining letters: `p` = parsimony-informative, `c` = constant, `s` = smart-gap (auto-threshold), `g` = gappy (fixed threshold).

**Example:**

```bash
# single file to stdout (smart-gap + parsimony filter by default)
$ phorge curate alignment.fasta > trimmed.fasta
Sites in:      10000
Sites kept:    4321
Sites removed: 5679

# batch with glob, output to directory
$ phorge curate aligned/*.fasta -o curated/
COI: 8000 → 3201 sites (4799 removed)
ND2: 6000 → 2874 sites (3126 removed)
12S: 4000 → 1950 sites (2050 removed)
Done. Curated 3 files.

# fixed gap threshold instead of smart-gap
$ phorge curate aligned/*.fasta -k pg --gap-threshold 0.5 -o curated/
```

**Flags:**
- `-k, --keep` — column properties to keep (combine letters); default `ps`
- `--gap-threshold` — max gappiness per column (0.0–1.0), used when `g` is in `--keep` (default: 0.9)
- `-e, --extension` — suffix to append to output filenames (default: `_curated`)
- `-o, --output` — output directory for trimmed files (if omitted, writes to stdout)

**Mode reference:**

| Flag | Description | ClipKIT equivalent |
|---|---|---|
| `-k ps` | parsimony-informative + smart-gap filter | `kpi-smart-gap` (default) |
| `-k pcs` | parsimony-informative + constant + smart-gap filter | `kpic-smart-gap` |
| `-k pg` | parsimony-informative + fixed gap filter | `kpi-gappy` |
| `-k pcg` | parsimony-informative + constant + fixed gap filter | `kpic-gappy` |
| `-k p` | parsimony-informative only | `kpi` |
| `-k pc` | parsimony-informative + constant | `kpic` |
| `-k s` | smart-gap filter only | `smart-gap` |
| `-k g` | fixed gap filter only | `gappy` |

Smart-gap (`s`) automatically determines the gap threshold from the distribution of per-site gappiness values, rather than requiring a fixed cutoff. This is the primary algorithm from ClipKIT and generally produces better results than a hardcoded threshold.

---
## Pipeline example

From taxon IDs to a supermatrix. The acquisition layer (`query → fetch → extract → clean`) turns a list of TaxIDs into deduplicated, provenance-labeled per-gene FASTAs; the self-contained tools (`align → curate → concat`) turn those into a supermatrix.

```bash
# 1. Acquire — TaxIDs in, one curated FASTA per gene out
phorge query --ingroup 89829 --outgroup 241031 309676 -o run/ --email you@example.org
phorge fetch   -q run/query_results.json -o run/ --email you@example.org --yes
phorge extract --refs refs/*.fasta -t run/raw/ -o run/genes/
phorge clean   --genes-dir run/genes/ -q run/query_results.json -o run/clean/ --prefer MyLab

# 2. Build — align, trim, and concatenate into a supermatrix
phorge align  -p mafft -i run/clean/*.fasta -e _aln -o run/aligned/
phorge curate run/aligned/*.fasta -o run/curated/
phorge concat run/curated/*.fasta > supermatrix.fasta
```

`refs/` holds one reference FASTA per gene (e.g. `COI.fasta`, `16S.fasta`); each may contain several sequences spanning your taxa to catch divergent hits.

---
## Planned Subcommands
- **scrub** — alignment outlier detection via pairwise p-distances
- **drafttree** — quick neighbor-joining tree from an MSA for sanity-checking alignments before committing to ML/Bayesian methods
- **view** — in-terminal alignment viewer
- **slice** — cut out or extract sections of an alignment

## Development Note

Phorge is being built as both a real research tool and a vehicle for learning Rust. Development is assisted by Claude (Anthropic), which serves as a teaching aid and coding partner. The design, domain knowledge, and direction are the author's own.

## Author

Andrew Budge
