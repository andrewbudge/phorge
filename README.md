<div align="center">

<img src="https://raw.githubusercontent.com/andrewbudge/phorge/main/docs/phorge-assets/phorge-mark-ember.svg" alt="phorge logo" width="200">

# Phorge

**A composable CLI toolkit for phylogenetics — from local sequence files to a supermatrix.**

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

Phorge is a single-binary toolkit for common steps in phylogenetic pipelines.
Its subcommands replace ad-hoc chains of shell commands with small, independent
tools that read files, write files, and compose in Unix pipelines. Sequence
acquisition is left to NCBI Datasets and GenBank; Phorge takes the files those
tools produce.

Every subcommand works on its own. A missing external tool affects only the
subcommands that need it.

> **Note:** Phorge is under active development. Subcommands may change or be
> added as the project matures.

## Contents

- [Install](#install)
- [External dependencies](#external-dependencies)
- **Subcommands:** [getheaders](#getheaders-ghd) · [stats](#stats) · [filter](#filter) · [convert](#convert) · [curate](#curate) · [concat](#concat-liger) · [coverage](#coverage) · [align](#align-aln) · [extract](#extract) · [unpack](#unpack) · [harvest](#harvest) · [report](#report) · [orthologs](#orthologs-groups) · [rip](#rip)

## Install

From a release tarball:

```bash
curl -LO https://github.com/andrewbudge/phorge/releases/download/v0.6.0/phorge-0.6.0-linux-x86_64.tar.gz
tar -xzf phorge-0.6.0-linux-x86_64.tar.gz
cp phorge-0.6.0-linux-x86_64/bin/phorge ~/.local/bin/
```

Release assets are built for Linux x86_64, macOS x86_64 and macOS arm64.

From source (requires a C11 compiler, zlib and minizip):

```bash
git clone https://github.com/andrewbudge/phorge.git
cd phorge
make
make install PREFIX=$HOME/.local
```

Run `phorge --help` for the command list, or `phorge <command> --help` for a
specific subcommand.

## External dependencies

Most subcommands are fully self-contained. A few orchestrate external programs:

| Subcommand | Requires |
|---|---|
| `extract` | [MMseqs2](https://github.com/soedinglab/MMseqs2) on PATH |
| `align` | [MAFFT](https://mafft.cbrc.jp/alignment/software/) or [MUSCLE](https://drive5.com/muscle/) on PATH |
| `unpack`, `report`, `orthologs` | An NCBI Datasets zip, or an unpacked Datasets directory |
| `rip` | Mitochondrial GenBank flat files (`.gb` / `.gbff`) |

 Install the alignment/search tools with conda if you do not have them:

```bash
conda install -c bioconda mmseqs2 mafft muscle
```

All file inputs — FASTA, taxa lists, query TSVs, NEXUS/PHYLIP, and Datasets
archives — accept gzip-compressed data transparently.

## Subcommands

### getheaders (ghd)

Extract FASTA headers in file order. `-u` keeps only unique headers. With no
files, reads stdin.

```bash
$ phorge getheaders alignment.fa
alpha first
beta second
alpha first

$ phorge getheaders -u alignment.fa
alpha first
beta second
```

### stats

Summarize one or more FASTA files: sequence count, length, DNA/AA type, GC%,
missingness, and (for equal-length alignments) variable and informative sites.

```bash
$ phorge stats alignment.fa
$ phorge stats -d -p alignment.fa     # per-sequence rows, aligned columns
$ phorge stats *.fa
```

### filter

Remove taxa from an alignment that exceed a missing-data threshold, or that
appear in too few loci according to a provenance TSV from `concat -l`.

```bash
phorge filter -m 0.5 alignment.fa > filtered.fa
phorge filter -n 20 -l provenance.tsv alignment.fa > filtered.fa
```

### convert

Detect FASTA, NEXUS or PHYLIP input from its contents and write FASTA, NEXUS,
strict PHYLIP or relaxed PHYLIP.

```bash
phorge convert -o nexus alignment.fa > alignment.nex
phorge convert -o fasta alignment.phy > alignment.fa
```

### curate

Trim alignment columns using ClipKIT-compatible modes: parsimony
informativeness, gappiness, smart-gap threshold, or combinations.

```bash
phorge curate alignment.fa                         # default kpi-smart-gap
phorge curate -k s alignment.fa                    # gap-aware trimming only
phorge curate -k pcg --gap-threshold 0.5 aln.fa    # fixed threshold
phorge curate -o trimmed/ alignment.fa             # write to directory
```

### concat (liger)

Concatenate per-gene alignments into a supermatrix. With no taxa list, taxa
are matched exactly by header; with `-a`, headers are matched to a taxa list
by case-insensitive substring. Writes optional provenance TSV.

```bash
phorge concat gene1.fa gene2.fa gene3.fa > supermatrix.fa
phorge concat -a taxa.txt -l provenance.tsv genes/*.fa > supermatrix.fa
phorge concat -f nexus -p nexus genes/*.fa > supermatrix.nex
```

### coverage

Summarize per-taxon or per-locus coverage from the provenance TSV produced by
`concat -l`.

```bash
phorge coverage -t provenance.tsv
phorge coverage -l -p provenance.tsv
```

### align (aln)

Batch-align FASTA files with MAFFT (default) or MUSCLE, writing results to an
output directory. Arguments after `--` are passed through to the aligner.

```bash
phorge align -o aligned/ unaligned/*.fa
phorge align -p muscle -o aligned/ unaligned/*.fa
phorge align -o aligned/ gene.fa -- --auto
```

### extract

Find gene regions in target FASTAs with MMseqs2 homology search. Minus-strand
hits are reverse-complemented into reference orientation.

```bash
phorge extract --reference refs.fna -o genes/ targets/*.fa
phorge extract --refs refs/ -o genes/ --flank 10 target.fa
```

### unpack

Extract an NCBI Datasets zip into a directory. Unsafe paths are skipped, and a
sorted manifest of entries is written to stdout.

```bash
phorge unpack --list ncbi_dataset.zip
phorge unpack ncbi_dataset.zip -o data/
```

### harvest

Collect all FASTA files from an unpacked Datasets directory into one flat
directory with collision-safe names. `--decompress` writes plain FASTA from
gzipped inputs.

```bash
phorge harvest data/ -o data/fasta --decompress
```

### report

Build a Phorge query TSV from a Datasets zip or unpacked directory. It joins
FASTA header tags with `data_report.jsonl` / `annotation_report.jsonl` and
handles both eukaryotic and prokaryotic packages.

```bash
phorge report ncbi_dataset.zip -o query.tsv
```

### orthologs (groups)

Split a Datasets ortholog package into one FASTA per orthogroup, named by
taxon and ready for alignment. The longest isoform per taxon is kept by
default. Writes an `orthologs.tsv` manifest.

```bash
phorge orthologs ncbi_dataset.zip -o genes/
phorge orthologs ncbi_dataset.zip -o genes/ --keep-all --id taxid
```

### rip

Extract selected loci from mitochondrial GenBank genomes. This tool is
specific to mitochondrial annotations and does not locate loci in general
GenBank files. Handles common mt gene synonyms, `complement()`/`join()`
locations, and minus-strand correction. Writes one FASTA per gene plus a TSV
log.

```bash
phorge rip -g COI,CYTB,12S,16S -o loci/ mitogenomes/*.gb
phorge rip --all --no-trna -o loci/ mitogenomes/*.gb
```

## License

MIT
