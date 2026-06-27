use crate::models::{Accession, QueryResult, TaxonGroup};
use crate::ncbi::EutilsClient;
use anyhow::{Context, Result, bail};
use clap::Args;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use tracing::{info, warn};

#[derive(Args)]
pub struct FetchArgs {
    /// Path to query_results.json (written by `query`)
    #[arg(long, short = 'q')]
    pub query: PathBuf,

    /// Output directory. Shards download here, then collapse into combined.fasta on success.
    #[arg(long, short = 'o')]
    pub out: PathBuf,

    /// Write the JSON log here instead of alongside the output (e.g. fast scratch).
    #[arg(long)]
    pub log_dir: Option<PathBuf>,

    /// Drop records shorter than this before downloading (preflight trim)
    #[arg(long)]
    pub min_length: Option<usize>,

    /// Drop records longer than this before downloading (preflight trim)
    #[arg(long)]
    pub max_length: Option<usize>,

    /// Email address required by NCBI ToS for automated access
    #[arg(long)]
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

/// Persistent, resumable record of a download. The authoritative resume signal is
/// which shard files exist on disk (see [`Manifest::reconcile`]); this document
/// carries the explicit accession-per-chunk mapping and provenance. Written
/// atomically (temp + rename) after every chunk so an interrupted run never
/// corrupts it.
#[derive(Serialize, Deserialize, Debug)]
struct Manifest {
    run_id: String,
    total_records: usize,
    est_bytes: u64,
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
    let out_dir = args.out.clone();
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

/// Load `query_results.json` and reduce it to the concrete set of accessions to
/// download. Metadata only — `slen` was captured at query time, so no network is
/// touched here. Order of operations matters: ingroup-wins overlap resolution
/// runs before dedup so a cross-group duplicate is always resolved in the
/// ingroup's favour.
fn load_and_preflight(args: &FetchArgs) -> Result<Vec<Accession>> {
    let content = std::fs::read_to_string(&args.query)
        .with_context(|| format!("reading {}", args.query.display()))?;
    let results: Vec<QueryResult> = serde_json::from_str(&content)
        .with_context(|| format!("parsing {}", args.query.display()))?;

    let mut records: Vec<Accession> = results.into_iter().flat_map(|r| r.accessions).collect();
    let records_in = records.len();

    // Ingroup wins: a sequence cannot honestly be both ingroup and outgroup, so
    // drop the outgroup copy of any accession that also appears in the ingroup.
    // Owned (not borrowed from `records`) so the set outlives the retain below.
    let ingroup_ids: HashSet<String> = records
        .iter()
        .filter(|a| a.taxon_group == TaxonGroup::Ingroup)
        .map(|a| a.accession.clone())
        .collect();
    let mut dropped_overlap = 0usize;
    records.retain(|a| {
        let drop = a.taxon_group == TaxonGroup::Outgroup && ingroup_ids.contains(&a.accession);
        dropped_overlap += usize::from(drop);
        !drop
    });

    // Dedup by accession string, first-seen wins. The query array lists ingroup
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
    let mut dropped_len = 0usize;
    if args.min_length.is_some() || args.max_length.is_some() {
        let min = args.min_length.unwrap_or(0);
        let max = args.max_length.unwrap_or(usize::MAX);
        records.retain(|a| {
            let keep = a.length >= min && a.length <= max;
            dropped_len += usize::from(!keep);
            keep
        });
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

impl Manifest {
    fn build(records: Vec<Accession>) -> Self {
        let total_records = records.len();
        let est_bytes = records
            .iter()
            .map(|a| a.length as u64 + FASTA_OVERHEAD)
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
            let mut reader = std::fs::File::open(&shard)
                .with_context(|| format!("reading shard {}", shard.display()))?;
            std::io::copy(&mut reader, &mut writer)
                .with_context(|| format!("appending shard {}", shard.display()))?;
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
    let mb = manifest.est_bytes as f64 / 1_048_576.0;
    info!(
        records = manifest.total_records,
        chunks = manifest.chunks.len(),
        est_mb = format!("{mb:.1}"),
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
        "About to download {} sequences in {} chunk(s) (~{mb:.1} MB). Continue? [y/N] ",
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
