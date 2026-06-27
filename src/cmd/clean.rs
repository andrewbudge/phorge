use crate::models::QueryResult;
use anyhow::{Context, Result};
use clap::Args;
use phorge::parse_fasta;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Args)]
pub struct CleanArgs {
    /// Directory of per-gene FASTAs written by `extract`
    #[arg(long)]
    pub genes_dir: PathBuf,

    /// Path to query_results.json (the accession -> TaxID/Name join table)
    #[arg(long, short = 'q')]
    pub query: PathBuf,

    /// Output directory
    #[arg(long, short = 'o')]
    pub out: PathBuf,

    /// Prefer records whose extract header or GenBank title contains this
    /// substring when deduplicating (e.g. `--prefer LabCode` to favour the lab's
    /// own vouchers). Repeatable; a match on any wins. Overrides the
    /// longest-sequence rule.
    #[arg(long)]
    pub prefer: Vec<String>,
}

/// Per-record provenance recovered from an extract output header, paired with
/// the join result. `taxid`/`name` come from query_results.json; `accession`
/// and `ident` come from the extract header itself.
struct CleanRecord {
    taxid: u64,
    name: String,
    accession: String,
    /// Extract's reported identity for this hit. Used only as the dedup
    /// tiebreaker when two records for the same TaxID are equally long.
    ident: f64,
    /// True when a `--prefer` substring matched this record's header or
    /// GenBank title. Outranks length in dedup so favoured vouchers always win.
    preferred: bool,
    seq: String,
}

/// The join target for one accession: everything `clean` needs from
/// query_results.json. `annotation` is the GenBank title, checked (alongside the
/// extract header) for `--prefer` substrings since vouchers like "BYU:IGCEP153"
/// often live only in the title, not the efetch defline.
struct Taxon {
    taxid: u64,
    name: String,
    annotation: String,
}

pub async fn run(args: CleanArgs) -> Result<()> {
    // Build the accession -> (TaxID, Name) join table from query_results.json.
    // The JSON is a flat array of QueryResult, each carrying its own accessions;
    // we flatten across all of them since the join key is the accession alone.
    let json = fs::read_to_string(&args.query)
        .with_context(|| format!("reading {}", args.query.display()))?;
    let results: Vec<QueryResult> =
        serde_json::from_str(&json).with_context(|| format!("parsing {}", args.query.display()))?;

    let mut join: HashMap<String, Taxon> = HashMap::new();
    for result in &results {
        for acc in &result.accessions {
            join.insert(
                acc.accession.clone(),
                Taxon {
                    taxid: acc.taxid,
                    name: acc.taxon_name.clone(),
                    annotation: acc.gene_annotation.clone(),
                },
            );
        }
    }

    fs::create_dir_all(&args.out)
        .with_context(|| format!("creating output directory {}", args.out.display()))?;

    // Collect the per-gene FASTAs extract emitted, in a stable order so the run
    // is deterministic regardless of how the filesystem hands them back.
    let mut gene_files: Vec<PathBuf> = fs::read_dir(&args.genes_dir)
        .with_context(|| format!("reading {}", args.genes_dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension()
                .is_some_and(|ext| ext == "fasta" || ext == "fa")
        })
        .collect();
    gene_files.sort();

    let mut total_kept = 0usize;
    let mut total_dropped = 0usize;
    let mut total_preferred = 0usize;
    let mut missing: Vec<String> = Vec::new();

    for path in &gene_files {
        // The filename stem is the authoritative gene label (extract names files
        // after the gene). The header's gene= tag agrees, but the stem is cleaner.
        let gene = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let path_str = path.to_string_lossy();
        let (records, _len) = parse_fasta(&path_str, false).map_err(|e| anyhow::anyhow!(e))?;

        // Join every record, dropping (and reporting) any whose accession is not
        // in query_results.json — that means a broken provenance chain, and an
        // unlabeled sequence is useless to concat's TaxID matching downstream.
        let mut joined: Vec<CleanRecord> = Vec::new();
        for (header, seq) in records {
            let accession = first_token(&header);
            let ident = parse_ident(&header);
            match join.get(accession) {
                Some(taxon) => {
                    // A --prefer substring may sit in the extract header (the
                    // preserved defline) or only in the GenBank title; check both.
                    let preferred = args.prefer.iter().any(|p| {
                        header.contains(p.as_str()) || taxon.annotation.contains(p.as_str())
                    });
                    joined.push(CleanRecord {
                        taxid: taxon.taxid,
                        name: taxon.name.clone(),
                        accession: accession.to_string(),
                        ident,
                        preferred,
                        seq,
                    });
                }
                None => missing.push(accession.to_string()),
            }
        }

        // Dedup per TaxID within this gene: keep the longest sequence, breaking
        // ties by highest extract identity.
        let before = joined.len();
        let kept = dedup_by_taxid(joined);
        total_dropped += before - kept.len();
        total_kept += kept.len();
        total_preferred += kept.iter().filter(|r| r.preferred).count();

        write_gene(&args.out, &gene, kept)?;
    }

    if !missing.is_empty() {
        missing.sort();
        missing.dedup();
        eprintln!(
            "warning: {} accession(s) had no match in {} and were dropped: {}",
            missing.len(),
            args.query.display(),
            missing.join(", ")
        );
    }
    eprintln!(
        "Done. Wrote {} cleaned sequence(s) across {} gene file(s); dropped {} duplicate(s).",
        total_kept,
        gene_files.len(),
        total_dropped
    );
    if !args.prefer.is_empty() {
        eprintln!(
            "  {} kept record(s) matched --prefer {:?}.",
            total_preferred, args.prefer
        );
    }

    Ok(())
}

/// First whitespace-delimited token of an extract header — the original NCBI
/// defline's accession (e.g. "MN908947.3" from ">MN908947.3 Severe acute...").
fn first_token(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or(header)
}

/// Pull the `ident=` value out of extract's trailing metadata block
/// (`[gene=COI ident=0.987 src=raw.fasta 10-660]`). Returns 0.0 if absent so a
/// malformed header simply sorts last in the dedup tiebreak rather than erroring.
fn parse_ident(header: &str) -> f64 {
    header
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix("ident="))
        .and_then(|v| v.trim_end_matches(']').parse().ok())
        .unwrap_or(0.0)
}

/// Keep one record per TaxID. Ranking, highest first: `--prefer` match, then
/// sequence length, then extract identity. Output is sorted by TaxID then
/// accession for deterministic files.
fn dedup_by_taxid(records: Vec<CleanRecord>) -> Vec<CleanRecord> {
    let mut best: HashMap<u64, CleanRecord> = HashMap::new();
    for rec in records {
        match best.get(&rec.taxid) {
            // Keep the incumbent only if it ranks at least as high as the newcomer.
            Some(cur) if outranks(cur, &rec) => {}
            _ => {
                best.insert(rec.taxid, rec);
            }
        }
    }
    let mut kept: Vec<CleanRecord> = best.into_values().collect();
    kept.sort_by(|a, b| a.taxid.cmp(&b.taxid).then(a.accession.cmp(&b.accession)));
    kept
}

/// True if `a` should be kept over `b` for the same TaxID: prefer-match wins,
/// then longest, then highest identity (identity only settles exact length ties).
fn outranks(a: &CleanRecord, b: &CleanRecord) -> bool {
    let rank = |r: &CleanRecord| (r.preferred, r.seq.len());
    rank(a) > rank(b) || (rank(a) == rank(b) && a.ident >= b.ident)
}

/// Write the cleaned records for one gene with rewritten `TaxID|Name|Accession|Gene`
/// headers. Spaces in the taxon name become underscores so each header stays a
/// single token (concat keys on the leading TaxID field).
fn write_gene(out_dir: &Path, gene: &str, records: Vec<CleanRecord>) -> Result<()> {
    let out_path = out_dir.join(format!("{gene}.fasta"));
    let file =
        File::create(&out_path).with_context(|| format!("creating {}", out_path.display()))?;
    let mut writer = BufWriter::new(file);
    for rec in records {
        let name = rec.name.replace(char::is_whitespace, "_");
        writeln!(writer, ">{}|{}|{}|{}", rec.taxid, name, rec.accession, gene)?;
        writeln!(writer, "{}", rec.seq)?;
    }
    Ok(())
}
