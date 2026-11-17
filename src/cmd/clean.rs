use crate::models::read_query_tsv;
use anyhow::{Context, Result};
use clap::Args;
use phorge::parse_fasta;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Args)]
pub struct CleanArgs {
    /// Input FASTA files written by `extract` (glob or list), or a directory
    /// containing them, e.g. `phorge clean genes/*.fasta -q ... -o ...`
    #[arg(required = true, num_args = 1..)]
    pub input: Vec<String>,

    /// Path to the query TSV (the accession -> TaxID/Name join table)
    #[arg(long, short = 'q')]
    pub query: PathBuf,

    /// Output directory
    #[arg(long, short = 'o')]
    pub output: PathBuf,

    /// Suffix to append to cleaned output filenames, for tracking pipeline
    /// stage (e.g. `COI.fasta` -> `COI_std.fasta`)
    #[arg(short, long, default_value = "_std")]
    pub extension: String,

    /// Prefer records whose extract header or GenBank title contains this
    /// substring when deduplicating (e.g. `--prefer LabCode` to favour the lab's
    /// own vouchers). Repeatable; a match on any wins. Overrides the
    /// longest-sequence rule.
    #[arg(long)]
    pub prefer: Vec<String>,
}

/// Per-record provenance recovered from an extract output header, paired with
/// the join result. `taxid`/`name` come from the query TSV; `accession`
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
/// the query TSV. `annotation` is the GenBank title, checked (alongside the
/// extract header) for `--prefer` substrings since vouchers like "BYU:IGCEP153"
/// often live only in the title, not the efetch defline.
struct Taxon {
    taxid: u64,
    name: String,
    annotation: String,
}

pub async fn run(args: CleanArgs) -> Result<()> {
    // Build the accession -> (TaxID, Name) join table from the query TSV. The
    // table is a flat list of accessions; the join key is the accession alone.
    let records = read_query_tsv(&args.query)?;
    let mut join: HashMap<String, Taxon> = HashMap::new();
    for acc in &records {
        join.insert(
            acc.accession.clone(),
            Taxon {
                taxid: acc.taxid,
                name: acc.taxon_name.clone(),
                annotation: acc.gene_annotation.clone(),
            },
        );
    }

    fs::create_dir_all(&args.output)
        .with_context(|| format!("creating output directory {}", args.output.display()))?;

    // Collect the per-gene FASTAs to clean, in a stable order so the run is
    // deterministic regardless of how the filesystem hands them back.
    let mut gene_files = collect_inputs(&args.input);
    gene_files.sort();

    let mut total_kept = 0usize;
    let mut total_dropped = 0usize;
    let mut total_preferred = 0usize;
    let mut unknown_taxid = 0usize;
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
        // in the query TSV — that means a broken provenance chain, and an
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
                    if taxon.taxid == 0 {
                        unknown_taxid += 1;
                    }
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

        write_gene(&args.output, &gene, &args.extension, kept)?;
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
    if unknown_taxid > 0 {
        // These records are kept one-per-accession rather than collapsed by
        // taxon, because a TaxID of 0 means the query file did not name one.
        eprintln!(
            "  warning: {} record(s) had no TaxID in {}; they were kept \
             per-accession instead of being deduplicated by taxon.",
            unknown_taxid,
            args.query.display()
        );
    }
    if !args.prefer.is_empty() {
        eprintln!(
            "  {} kept record(s) matched --prefer {:?}.",
            total_preferred, args.prefer
        );
    }

    Ok(())
}

/// Collect the FASTA files to clean. Each input may be a file, or a directory
/// containing them (the shell has already expanded any glob before we see it) —
/// the same "file or directory, mixed freely" convention `align` and `extract`
/// use. A directory is scanned for alignment-extension files; a file named
/// explicitly is used regardless of its extension.
fn collect_inputs(inputs: &[String]) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = Vec::new();

    for input in inputs {
        let path = Path::new(input);
        if path.is_dir() {
            for entry in path.read_dir().unwrap_or_else(|e| {
                eprintln!("Error: could not read directory '{}': {}", input, e);
                std::process::exit(1);
            }) {
                let entry = entry.expect("Failed to read directory entry");
                let p = entry.path();
                if p.extension()
                    .is_some_and(|e| e == "fasta" || e == "fa" || e == "fna" || e == "fas")
                {
                    files.push(p);
                }
            }
        } else if path.is_file() {
            files.push(path.to_path_buf());
        } else {
            eprintln!("Warning: '{}' is not a file or directory, skipping.", input);
        }
    }

    files
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
///
/// Records whose TaxID is unknown (0) are keyed by accession instead — see
/// [`dedup_key`].
fn dedup_by_taxid(records: Vec<CleanRecord>) -> Vec<CleanRecord> {
    let mut best: HashMap<DedupKey, CleanRecord> = HashMap::new();
    for rec in records {
        let key = dedup_key(&rec);
        match best.get(&key) {
            // Keep the incumbent only if it ranks at least as high as the newcomer.
            Some(cur) if outranks(cur, &rec) => {}
            _ => {
                best.insert(key, rec);
            }
        }
    }
    let mut kept: Vec<(DedupKey, CleanRecord)> = best.into_iter().collect();
    kept.sort_by(|(ka, a), (kb, b)| ka.cmp(kb).then(a.accession.cmp(&b.accession)));
    kept.into_iter().map(|(_, rec)| rec).collect()
}

/// The identity a record is deduplicated by: its TaxID, or its accession when
/// the TaxID is unknown.
///
/// A TaxID of 0 means "not known" — `read_query_tsv` fills it in when the
/// query file has no `taxid` column, which happens whenever the user filters
/// the TSV down with `cut`/`awk`. Deduplicating those records by TaxID would
/// collapse every unnamed accession into a single survivor and silently delete
/// the rest, so they are kept one-per-accession instead.
///
/// Variant order matters: `Unknown` sorts before `Taxid` so records with an
/// unknown TaxID keep the position the old numeric sort gave them (0 before
/// every real TaxID).
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash)]
enum DedupKey {
    Unknown(String),
    Taxid(u64),
}

fn dedup_key(rec: &CleanRecord) -> DedupKey {
    if rec.taxid == 0 {
        DedupKey::Unknown(rec.accession.clone())
    } else {
        DedupKey::Taxid(rec.taxid)
    }
}

/// True if `a` should be kept over `b` for the same TaxID: prefer-match wins,
/// then longest, then highest identity (identity only settles exact length ties).
fn outranks(a: &CleanRecord, b: &CleanRecord) -> bool {
    let rank = |r: &CleanRecord| (r.preferred, r.seq.len());
    rank(a) > rank(b) || (rank(a) == rank(b) && a.ident >= b.ident)
}

/// Write the cleaned records for one gene with rewritten `TaxID|Name|Accession|Gene`
/// headers. Spaces in the taxon name become underscores so each header stays a
/// single token (concat keys on the leading TaxID field). `extension` is appended
/// to the gene name so the output filename carries its pipeline stage, matching
/// how `align` names its own output (e.g. `COI.fasta` -> `COI_std.fasta`).
fn write_gene(
    out_dir: &Path,
    gene: &str,
    extension: &str,
    records: Vec<CleanRecord>,
) -> Result<()> {
    let out_path = out_dir.join(format!("{gene}{extension}.fasta"));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn record(taxid: u64, accession: &str, seq_len: usize) -> CleanRecord {
        CleanRecord {
            taxid,
            name: "Felis catus".to_string(),
            accession: accession.to_string(),
            ident: 0.99,
            preferred: false,
            seq: "A".repeat(seq_len),
        }
    }

    #[test]
    fn known_taxids_deduplicate_to_one_record_per_taxon() {
        let kept = dedup_by_taxid(vec![
            record(9901, "ACC1", 100),
            record(9901, "ACC2", 200),
            record(9902, "ACC3", 300),
        ]);
        assert_eq!(kept.len(), 2, "one record per TaxID");
        let accessions: Vec<&str> = kept.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(accessions, ["ACC2", "ACC3"], "longest sequence wins");
    }

    #[test]
    fn unknown_taxids_are_not_collapsed_into_each_other() {
        // A query TSV with no `taxid` column fills the field with 0. Keying
        // those records by TaxID would silently delete all but one of them.
        let kept = dedup_by_taxid(vec![
            record(0, "ACC1", 100),
            record(0, "ACC2", 200),
            record(0, "ACC3", 300),
        ]);
        assert_eq!(
            kept.len(),
            3,
            "records with an unknown TaxID are kept one per accession"
        );
        let accessions: Vec<&str> = kept.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(accessions, ["ACC1", "ACC2", "ACC3"]);
    }

    #[test]
    fn unknown_and_known_taxids_do_not_interfere() {
        let kept = dedup_by_taxid(vec![
            record(0, "ACC1", 100),
            record(0, "ACC2", 200),
            record(9901, "ACC3", 300),
            record(9901, "ACC4", 400),
        ]);
        assert_eq!(kept.len(), 3);
        let accessions: Vec<&str> = kept.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(accessions, ["ACC1", "ACC2", "ACC4"]);
    }
}
