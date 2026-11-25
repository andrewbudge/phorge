use crate::models::TaxonGroup;
use crate::ncbi::EutilsClient;
use anyhow::{Context, Result, bail};
use clap::Args;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{IsTerminal, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use tracing::{info, warn};

#[derive(Args)]
pub struct FetchArgs {
    /// Input: either a query TSV written by `query`, or a bare list of accessions
    /// (one per line, blank lines and `#` comments skipped) — e.g. from your own
    /// curation or piped in from another tool. Detected automatically from the
    /// first line. Pass `-` to read from stdin, e.g.
    /// `cut -f1 query.tsv | phorge fetch - -o run/ ...`. A bare list carries no
    /// length or ingroup/outgroup metadata, so --min-length/--max-length and
    /// cross-group dedup only take effect for a query TSV.
    #[arg(required = true)]
    pub query: PathBuf,

    /// Output directory. Shards download here, then collapse into combined.fasta on success.
    #[arg(long, short = 'o')]
    pub output: PathBuf,

    /// Drop records shorter than this before downloading (preflight trim). Only
    /// applies when the input is a query TSV with a length column — a bare
    /// accession list has no length to filter on, and the filter is skipped with
    /// a warning.
    #[arg(long)]
    pub min_length: Option<usize>,

    /// Drop records longer than this before downloading (preflight trim). Same
    /// query-TSV-only limitation as --min-length.
    #[arg(long)]
    pub max_length: Option<usize>,

    /// Email address required by NCBI ToS for automated access
    #[arg(long, short = 'e')]
    pub email: String,

    /// NCBI API key (optional; raises the NCBI rate limit from 3 to 10 req/s)
    #[arg(long)]
    pub api_key: Option<String>,

    /// Skip confirmation prompt (for non-interactive use)
    #[arg(long)]
    pub yes: bool,
}

/// Nucleotide database to fetch from (GenBank/RefSeq nucleotide set).
const DB: &str = "nuccore";

/// Accessions per efetch POST. This is also the resumable unit: one shard file
/// per chunk, so a crash costs at most one chunk's worth of re-download.
const CHUNK_SIZE: usize = 500;

/// Rough FASTA header + line-wrap overhead per record, added to `slen` for the
/// pre-download size estimate. Deliberately approximate — it only sizes the
/// confirmation prompt, nothing downstream depends on it.
const FASTA_OVERHEAD: u64 = 80;

/// Inline retry attempts per chunk before it is left `Failed` for a later resume.
const MAX_CHUNK_RETRIES: u32 = 3;

/// One record to download, reduced from either a query TSV row or a bare
/// accession line. `length` and `taxon_group` are only ever populated together,
/// from a query TSV's `length`/`taxon_group` columns — a bare accession list has
/// neither, which is why the ingroup-wins overlap rule and length filtering both
/// degrade to no-ops for it (see [`load_and_preflight`]).
struct FetchRecord {
    accession: String,
    length: Option<usize>,
    taxon_group: Option<TaxonGroup>,
}

/// Persistent, resumable record of a download. The authoritative resume signal is
/// which shard files exist on disk (see [`Manifest::reconcile`]); this document
/// carries the explicit accession-per-chunk mapping and provenance. Written
/// atomically (temp + rename) after every chunk so an interrupted run never
/// corrupts it.
#[derive(Serialize, Deserialize, Debug)]
struct Manifest {
    run_id: String,
    total_records: usize,
    /// `None` when any input record's length is unknown (i.e. a bare accession
    /// list) — an estimate built from a mix of real and assumed-zero lengths
    /// would be actively misleading, so the confirmation prompt reports "unknown"
    /// instead.
    est_bytes: Option<u64>,
    chunks: Vec<Chunk>,
}

#[derive(Serialize, Deserialize, Debug)]
struct Chunk {
    /// Drives the shard filename (`shard_0003.fasta`) and never changes.
    index: usize,
    accessions: Vec<String>,
    state: ChunkState,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum ChunkState {
    Pending,
    Done,
    Failed,
}

pub async fn run(args: FetchArgs) -> Result<()> {
    let out_dir = args.output.clone();
    let combined_path = out_dir.join("combined.fasta");

    // A combined file exists only once a prior run fully succeeded and collapsed
    // its shards. Treat its presence as "already done" so a re-run is a no-op;
    // delete it to force a fresh fetch.
    if combined_path.exists() {
        info!(output = %combined_path.display(), "combined output already present; nothing to do");
        return Ok(());
    }

    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;
    let manifest_path = out_dir.join("download_manifest.json");

    // Resume an existing run, or build a fresh manifest from query results. The
    // confirmation gate only fires on a fresh run — a resume was already approved.
    let mut manifest = if manifest_path.exists() {
        let m = Manifest::load(&manifest_path)
            .with_context(|| format!("reading {}", manifest_path.display()))?;
        info!(run_id = %m.run_id, chunks = m.chunks.len(), "resuming from existing manifest");
        m
    } else {
        let records = load_and_preflight(&args)?;
        let manifest = Manifest::build(records);
        confirm(&manifest, args.yes)?;
        manifest
    };

    // Shards on disk are the source of truth: mark present shards Done, and redo
    // any chunk whose shard has vanished since the manifest was last written.
    manifest.reconcile(&out_dir);
    manifest
        .save(&manifest_path)
        .with_context(|| format!("writing {}", manifest_path.display()))?;

    let client = EutilsClient::new(args.api_key, args.email).context("building NCBI client")?;
    download(&client, &mut manifest, &out_dir, &manifest_path).await?;

    // download() only returns Ok when every chunk is Done, so all shards are on
    // disk. Collapse them into one multifasta and drop the now-redundant shards.
    consolidate(&out_dir, &combined_path, &manifest)?;

    info!(
        chunks = manifest.chunks.len(),
        records = manifest.total_records,
        output = %combined_path.display(),
        "fetch complete"
    );
    Ok(())
}

/// Load the query input and reduce it to the concrete set of accessions to
/// download. No network is touched here. Order of operations matters:
/// ingroup-wins overlap resolution runs before dedup so a cross-group duplicate
/// is always resolved in the ingroup's favour.
fn load_and_preflight(args: &FetchArgs) -> Result<Vec<FetchRecord>> {
    let content = read_query_input(&args.query)?;
    let mut records =
        parse_query_input(&content).with_context(|| format!("parsing {}", args.query.display()))?;
    let records_in = records.len();

    // Ingroup wins: a sequence cannot honestly be both ingroup and outgroup, so
    // drop the outgroup copy of any accession that also appears in the ingroup.
    // A bare accession list has no taxon_group at all, so both sides of this
    // comparison are always `None` for it and the retain is a no-op — no
    // special-casing needed. Owned (not borrowed from `records`) so the set
    // outlives the retain below.
    let ingroup_ids: HashSet<String> = records
        .iter()
        .filter(|a| a.taxon_group == Some(TaxonGroup::Ingroup))
        .map(|a| a.accession.clone())
        .collect();
    let mut dropped_overlap = 0usize;
    records.retain(|a| {
        let drop =
            a.taxon_group == Some(TaxonGroup::Outgroup) && ingroup_ids.contains(&a.accession);
        dropped_overlap += usize::from(drop);
        !drop
    });

    // Dedup by accession string, first-seen wins. The query TSV lists ingroup
    // taxa first, so first-seen preserves ingroup provenance for within-set dupes.
    let mut seen: HashSet<String> = HashSet::with_capacity(records.len());
    let mut dropped_dup = 0usize;
    records.retain(|a| {
        let fresh = seen.insert(a.accession.clone());
        dropped_dup += usize::from(!fresh);
        fresh
    });

    // Optional, off-by-default length trim. The byte-gate below is the primary
    // cost control; these bounds are an opt-in trim for obvious genomes/fragments.
    // Only a query TSV carries lengths; a bare accession list has none to filter
    // on, so we warn and fetch everything rather than guessing.
    let mut dropped_len = 0usize;
    if args.min_length.is_some() || args.max_length.is_some() {
        if records.iter().any(|a| a.length.is_some()) {
            let min = args.min_length.unwrap_or(0);
            let max = args.max_length.unwrap_or(usize::MAX);
            records.retain(|a| {
                let keep = a.length.is_none_or(|len| len >= min && len <= max);
                dropped_len += usize::from(!keep);
                keep
            });
        } else {
            warn!(
                "--min-length/--max-length need a query TSV's length column; \
                 ignoring for this bare accession list"
            );
        }
    }

    info!(
        records_in,
        dropped_overlap,
        dropped_dup,
        dropped_len,
        records_out = records.len(),
        "preflight complete"
    );

    if records.is_empty() {
        bail!("no records left to fetch after preflight");
    }
    Ok(records)
}

/// Read the query argument: a real path, or stdin when it is exactly `-`.
/// This is what makes `cut -f1 query.tsv | phorge fetch - ...` work — the
/// sniffing in [`parse_query_input`] then decides TSV vs bare-list on whatever
/// came through, same as it would for a file.
fn read_query_input(path: &Path) -> Result<String> {
    if path == Path::new("-") {
        let mut buf = String::new();
        std::io::stdin()
            .read_to_string(&mut buf)
            .context("reading query input from stdin")?;
        Ok(buf)
    } else {
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))
    }
}

/// Sniff the query input and parse it into fetch records. A query TSV is
/// detected by its header row: the first non-blank line contains a tab. Anything
/// else is treated as a bare accession list — one accession per line, blank
/// lines and `#` comments skipped — which is how a user's own hand-curated or
/// tool-generated accession list gets in.
fn parse_query_input(content: &str) -> Result<Vec<FetchRecord>> {
    let mut lines = content.lines().map(str::trim).filter(|l| !l.is_empty());
    let Some(first) = lines.next() else {
        bail!("query input is empty");
    };

    if first.contains('\t') {
        parse_query_tsv(first, lines)
    } else {
        Ok(std::iter::once(first)
            .chain(lines)
            .filter(|l| !l.starts_with('#'))
            .map(|accession| FetchRecord {
                accession: accession.to_string(),
                length: None,
                taxon_group: None,
            })
            .collect())
    }
}

/// Parse a tab-separated query file: `header` is the already-consumed first
/// line, `rows` the remaining lines. Columns are matched by name
/// (case-insensitive) rather than position, so column order — and which extra
/// columns a user's own filtering with `cut`/`awk` leaves in — doesn't matter.
/// Only `accession` is required; `length` and `taxon_group` are read when
/// present so length filtering and ingroup-wins overlap resolution can run.
fn parse_query_tsv<'a>(
    header: &str,
    rows: impl Iterator<Item = &'a str>,
) -> Result<Vec<FetchRecord>> {
    let columns: Vec<&str> = header.split('\t').map(str::trim).collect();
    let col_index = |name: &str| columns.iter().position(|c| c.eq_ignore_ascii_case(name));
    let accession_col = col_index("accession").context("query TSV has no 'accession' column")?;
    let length_col = col_index("length");
    let group_col = col_index("taxon_group");

    let mut records = Vec::new();
    for row in rows {
        let fields: Vec<&str> = row.split('\t').collect();
        // A ragged row (e.g. a trailing blank field trimmed by hand) is skipped
        // rather than panicking — the file might have been hand-edited.
        let Some(accession) = fields.get(accession_col) else {
            continue;
        };
        let length = length_col
            .and_then(|i| fields.get(i))
            .and_then(|s| s.trim().parse::<usize>().ok());
        let taxon_group = group_col.and_then(|i| fields.get(i)).and_then(|s| {
            match s.trim().to_ascii_lowercase().as_str() {
                "ingroup" => Some(TaxonGroup::Ingroup),
                "outgroup" => Some(TaxonGroup::Outgroup),
                _ => None,
            }
        });
        records.push(FetchRecord {
            accession: accession.trim().to_string(),
            length,
            taxon_group,
        });
    }
    Ok(records)
}

impl Manifest {
    fn build(records: Vec<FetchRecord>) -> Self {
        let total_records = records.len();
        // `Option<u64>: Sum<Option<u64>>` collapses to `None` as soon as any
        // record's length is unknown, which is exactly the "all or nothing"
        // behaviour we want: a bare accession list has no lengths at all, so the
        // whole estimate is reported as unknown rather than partially guessed.
        let est_bytes: Option<u64> = records
            .iter()
            .map(|a| a.length.map(|len| len as u64 + FASTA_OVERHEAD))
            .sum();
        let chunks = records
            .chunks(CHUNK_SIZE)
            .enumerate()
            .map(|(index, recs)| Chunk {
                index,
                accessions: recs.iter().map(|a| a.accession.clone()).collect(),
                state: ChunkState::Pending,
            })
            .collect();
        Self {
            run_id: run_id(),
            total_records,
            est_bytes,
            chunks,
        }
    }

    fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)?;
        Ok(serde_json::from_str(&content)?)
    }

    /// Write atomically: a temp file in the same directory, then rename. An
    /// interrupted write leaves the previous manifest (or just a stray temp)
    /// intact rather than a half-written one.
    fn save(&self, path: &Path) -> Result<()> {
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let mut tmp = tempfile::NamedTempFile::new_in(dir).context("creating temp manifest")?;
        tmp.write_all(serde_json::to_string_pretty(self)?.as_bytes())
            .context("writing temp manifest")?;
        tmp.as_file().sync_all().context("flushing temp manifest")?;
        tmp.persist(path)
            .map_err(|e| anyhow::anyhow!("persisting manifest: {e}"))?;
        Ok(())
    }

    /// Bring chunk states in line with what is actually on disk. Present shard =>
    /// Done; a Done chunk whose shard disappeared => back to Pending.
    fn reconcile(&mut self, out_dir: &Path) {
        for chunk in &mut self.chunks {
            let exists = out_dir.join(shard_name(chunk.index)).exists();
            chunk.state = match (chunk.state, exists) {
                (_, true) => ChunkState::Done,
                (ChunkState::Done, false) => ChunkState::Pending,
                (other, false) => other,
            };
        }
    }
}

/// Fetch every not-yet-Done chunk. Downloads what it can: a chunk that fails
/// after retries is marked `Failed` and the loop continues, so one bad chunk
/// doesn't strand the rest. If anything failed, the command errors at the end
/// with a resume hint.
async fn download(
    client: &EutilsClient,
    manifest: &mut Manifest,
    out_dir: &Path,
    manifest_path: &Path,
) -> Result<()> {
    let pending: Vec<usize> = (0..manifest.chunks.len())
        .filter(|&i| manifest.chunks[i].state != ChunkState::Done)
        .collect();
    if pending.is_empty() {
        info!("all shards already present; nothing to download");
        return Ok(());
    }
    info!(
        pending = pending.len(),
        total = manifest.chunks.len(),
        "starting download"
    );

    let mut failed = 0usize;
    for i in pending {
        let index = manifest.chunks[i].index;
        let ids = manifest.chunks[i].accessions.clone();
        let id_refs: Vec<&str> = ids.iter().map(String::as_str).collect();

        let mut attempt = 0u32;
        let outcome = loop {
            attempt += 1;
            match client.efetch_fasta(DB, &id_refs).await {
                Ok(body) => break Ok(body),
                Err(e) if attempt < MAX_CHUNK_RETRIES => {
                    let backoff = Duration::from_secs(1u64 << attempt);
                    warn!(
                        chunk = index,
                        attempt,
                        error = %e,
                        backoff_s = backoff.as_secs(),
                        "chunk fetch failed; retrying"
                    );
                    tokio::time::sleep(backoff).await;
                }
                Err(e) => break Err(e),
            }
        };

        match outcome {
            Ok(body) => {
                write_shard(out_dir, index, &body)
                    .with_context(|| format!("writing shard {index}"))?;
                manifest.chunks[i].state = ChunkState::Done;
                info!(chunk = index, records = ids.len(), "shard written");
            }
            Err(e) => {
                manifest.chunks[i].state = ChunkState::Failed;
                failed += 1;
                warn!(chunk = index, error = %e, "chunk failed after retries; left for resume");
            }
        }
        manifest
            .save(manifest_path)
            .with_context(|| format!("updating {}", manifest_path.display()))?;
    }

    if failed > 0 {
        bail!("{failed} chunk(s) failed; re-run the same command to resume the remaining work");
    }
    Ok(())
}

/// Write a shard atomically: temp file in `out_dir`, fsync, then rename onto the
/// final name. A crash before the rename leaves only the temp file, which a
/// resume ignores (it keys off the final shard name) — so no duplicates and no
/// partial-record corruption are ever possible.
fn write_shard(out_dir: &Path, index: usize, body: &str) -> Result<()> {
    let mut tmp = tempfile::NamedTempFile::new_in(out_dir).context("creating temp shard")?;
    tmp.write_all(body.as_bytes())
        .context("writing temp shard")?;
    tmp.as_file().sync_all().context("flushing temp shard")?;
    tmp.persist(out_dir.join(shard_name(index)))
        .map_err(|e| anyhow::anyhow!("persisting shard {index}: {e}"))?;
    Ok(())
}

/// Collapse every shard into a single multifasta, then delete the shards. Only
/// called once the download fully succeeds, so every shard named in the manifest
/// is present. The combined file is written atomically (temp + rename) so a crash
/// mid-merge never leaves a half-built file a re-run would mistake for done; the
/// shards are removed only after that rename lands.
fn consolidate(out_dir: &Path, combined_path: &Path, manifest: &Manifest) -> Result<()> {
    let mut tmp =
        tempfile::NamedTempFile::new_in(out_dir).context("creating temp combined file")?;
    {
        let mut writer = std::io::BufWriter::new(tmp.as_file_mut());
        for chunk in &manifest.chunks {
            let shard = out_dir.join(shard_name(chunk.index));
            let body = std::fs::read_to_string(&shard)
                .with_context(|| format!("reading shard {}", shard.display()))?;
            // efetch separates records with a blank line; copied verbatim those
            // would leave gaps between sequences in the combined file. Drop blank
            // lines so the output is gap-free FASTA — headers and sequence lines
            // pass through unchanged.
            for line in body.lines().filter(|l| !l.trim().is_empty()) {
                writeln!(writer, "{line}")
                    .with_context(|| format!("appending shard {}", shard.display()))?;
            }
        }
        writer.flush().context("flushing combined file")?;
    }
    tmp.as_file().sync_all().context("syncing combined file")?;
    tmp.persist(combined_path)
        .map_err(|e| anyhow::anyhow!("persisting combined file: {e}"))?;

    // Combined file is durably on disk; the shards are now redundant.
    for chunk in &manifest.chunks {
        let shard = out_dir.join(shard_name(chunk.index));
        std::fs::remove_file(&shard)
            .with_context(|| format!("removing shard {}", shard.display()))?;
    }
    Ok(())
}

fn shard_name(index: usize) -> String {
    format!("shard_{index:04}.fasta")
}

fn run_id() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    format!("fetch-{secs}")
}

/// Show the download size and require confirmation before any bytes move. This
/// is the cost gate: `--yes` skips it; a non-interactive shell without `--yes`
/// errors rather than silently downloading. The prompt is interactive I/O, so it
/// goes straight to stderr rather than through the log.
fn confirm(manifest: &Manifest, yes: bool) -> Result<()> {
    // A bare accession list carries no lengths, so there's nothing to estimate
    // size from — say so plainly rather than showing a fabricated number.
    let size_desc = match manifest.est_bytes {
        Some(bytes) => format!("~{:.1} MB", bytes as f64 / 1_048_576.0),
        None => "unknown size (bare accession list; no length data)".to_string(),
    };
    info!(
        records = manifest.total_records,
        chunks = manifest.chunks.len(),
        est_size = %size_desc,
        "preflight ready to download"
    );
    if yes {
        return Ok(());
    }
    if !std::io::stdin().is_terminal() {
        bail!(
            "refusing to download non-interactively without confirmation: re-run with --yes \
             (this shell is not a TTY)"
        );
    }
    eprint!(
        "About to download {} sequences in {} chunk(s) ({size_desc}). Continue? [y/N] ",
        manifest.total_records,
        manifest.chunks.len()
    );
    std::io::stderr().flush().ok();
    let mut input = String::new();
    std::io::stdin()
        .read_line(&mut input)
        .context("reading confirmation")?;
    if !matches!(input.trim(), "y" | "Y" | "yes" | "Yes") {
        bail!("aborted by user");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bare_accession_list_skips_blanks_and_comments() {
        let input = "# my curated accessions\nAB123456.1\nCD654321.1\n\nEF999999.1\n";
        let records = parse_query_input(input).unwrap();
        let accessions: Vec<&str> = records.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(accessions, ["AB123456.1", "CD654321.1", "EF999999.1"]);
        // No metadata comes from a bare list — length filtering and ingroup-wins
        // overlap resolution both rely on this being None.
        assert!(records.iter().all(|r| r.length.is_none()));
        assert!(records.iter().all(|r| r.taxon_group.is_none()));
    }

    #[test]
    fn query_tsv_is_detected_and_columns_read_by_name() {
        // Column order deliberately doesn't match the field declaration order,
        // to prove lookup is by name, not position.
        let input = "taxon_group\taccession\tlength\n\
                      ingroup\tAB123456.1\t650\n\
                      outgroup\tCD654321.1\t720\n";
        let records = parse_query_input(input).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].accession, "AB123456.1");
        assert_eq!(records[0].length, Some(650));
        assert_eq!(records[0].taxon_group, Some(TaxonGroup::Ingroup));
        assert_eq!(records[1].taxon_group, Some(TaxonGroup::Outgroup));
    }

    #[test]
    fn tsv_without_accession_column_is_rejected() {
        let input = "taxon_name\tlength\nfelis catus\t650\n";
        assert!(parse_query_input(input).is_err());
    }
}
