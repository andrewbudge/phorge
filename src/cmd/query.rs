use crate::models::{Accession, QueryResult, TaxonGroup};
use crate::ncbi::EutilsClient;
use anyhow::{Context, Result};
use clap::Args;
use serde_json::Value;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

#[derive(Args)]
pub struct QueryArgs {
    /// One or more ingroup TaxIDs (e.g. 7088 for Lepidoptera)
    #[arg(long, num_args = 1.., required = true)]
    pub ingroup: Vec<u64>,

    /// One or more outgroup TaxIDs (e.g. a few representatives of sister genera)
    #[arg(long, num_args = 1..)]
    pub outgroup: Vec<u64>,

    /// Output directory
    #[arg(long, short = 'o')]
    pub out: PathBuf,

    /// Write the JSON log here instead of alongside the output (e.g. fast scratch).
    #[arg(long)]
    pub log_dir: Option<PathBuf>,

    /// Email address required by NCBI ToS for automated access
    #[arg(long)]
    pub email: String,

    /// NCBI API key (optional; raises the NCBI rate limit from 3 to 10 req/s)
    #[arg(long)]
    pub api_key: Option<String>,
}

/// Nucleotide database to search. `nuccore` is the canonical name for the
/// GenBank/RefSeq nucleotide set.
const DB: &str = "nuccore";

/// esummary docsums fetched per request. NCBI tolerates more, but 500 keeps each
/// response modest and matches the fetch-stage batch size.
const PAGE_SIZE: usize = 500;

pub async fn run(args: QueryArgs) -> Result<()> {
    let client = EutilsClient::new(args.api_key, args.email).context("building NCBI client")?;

    // One QueryResult per taxon, ingroup first then outgroup. Both groups are
    // queried identically; the only difference is the TaxonGroup tag, which
    // downstream stages (e.g. fetch's ingroup-wins overlap rule) rely on.
    let mut results: Vec<QueryResult> =
        Vec::with_capacity(args.ingroup.len() + args.outgroup.len());
    for &taxid in &args.ingroup {
        results.push(query_taxon(&client, taxid, TaxonGroup::Ingroup).await?);
    }
    for &taxid in &args.outgroup {
        results.push(query_taxon(&client, taxid, TaxonGroup::Outgroup).await?);
    }

    warn_cross_group_overlap(&results);

    std::fs::create_dir_all(&args.out)
        .with_context(|| format!("creating output directory {}", args.out.display()))?;
    let out_path = args.out.join("query_results.json");
    // Top-level array: one element per queried taxon across both groups.
    let json = serde_json::to_string_pretty(&results)?;
    std::fs::write(&out_path, json).with_context(|| format!("writing {}", out_path.display()))?;

    print_summary(&results, &out_path);
    Ok(())
}

/// Run one esearch + esummary sweep over a single taxon and collect its
/// accessions. `[Organism:exp]` excludes environmental samples and expands the
/// TaxID to its full subtree. No gene-name filtering — homology is MMseqs2's job
/// in `extract`.
async fn query_taxon(client: &EutilsClient, taxid: u64, group: TaxonGroup) -> Result<QueryResult> {
    let term = format!("txid{taxid}[Organism:exp]");
    eprintln!("querying {DB} for {term}");

    let handle = client
        .esearch_history(DB, &term)
        .await
        .with_context(|| format!("esearch failed for {term}"))?;

    // A missing taxonomy name is not worth aborting a successful search over —
    // degrade to the bare TaxID and warn.
    let taxon_name = match client.taxonomy_name(taxid).await {
        Ok(name) => name,
        Err(e) => {
            eprintln!("warning: could not resolve a name for txid {taxid} ({e}); using the TaxID");
            format!("txid{taxid}")
        }
    };

    eprintln!(
        "{taxon_name} ({taxid}): {} records found; retrieving metadata...",
        handle.count
    );

    let mut accessions: Vec<Accession> = Vec::with_capacity(handle.count);
    let mut retstart = 0;
    while retstart < handle.count {
        let page = client
            .esummary_page(DB, &handle, retstart, PAGE_SIZE)
            .await
            .with_context(|| format!("esummary failed at offset {retstart}"))?;
        let parsed = parse_docsums(&page, taxid, group, &mut accessions)
            .with_context(|| format!("parsing esummary page at offset {retstart}"))?;
        // Guard against an empty page stalling the loop if NCBI returns fewer
        // records than advertised.
        if parsed == 0 {
            break;
        }
        retstart += PAGE_SIZE;
    }

    Ok(QueryResult {
        taxid,
        taxon_name,
        taxon_group: group,
        total_accessions: accessions.len(),
        accessions,
    })
}

/// Pull the fields phorge needs out of an esummary docsum page, appending an
/// [`Accession`] per record. Each record is stamped with the `query_taxid` and
/// `group` that surfaced it, so provenance survives later flattening. Returns the
/// number parsed. Records missing an accession (or carrying an error) are skipped
/// rather than failing the page.
fn parse_docsums(
    page: &Value,
    query_taxid: u64,
    group: TaxonGroup,
    out: &mut Vec<Accession>,
) -> Result<usize> {
    let result = page
        .get("result")
        .ok_or_else(|| anyhow::anyhow!("esummary response missing 'result'"))?;
    let uids = result
        .get("uids")
        .and_then(|v| v.as_array())
        .ok_or_else(|| anyhow::anyhow!("esummary response missing 'result.uids'"))?;

    let mut parsed = 0;
    for uid in uids {
        let Some(uid) = uid.as_str() else { continue };
        let Some(doc) = result.get(uid) else { continue };
        if doc.get("error").is_some() {
            continue;
        }

        let accession = doc
            .get("accessionversion")
            .and_then(|v| v.as_str())
            .or_else(|| doc.get("caption").and_then(|v| v.as_str()))
            .unwrap_or_default()
            .to_string();
        if accession.is_empty() {
            continue;
        }

        let sourcedb = doc.get("sourcedb").and_then(|v| v.as_str()).unwrap_or("");
        out.push(Accession {
            accession,
            taxon_name: doc
                .get("organism")
                .and_then(|v| v.as_str())
                .unwrap_or_default()
                .to_string(),
            taxid: doc.get("taxid").and_then(|v| v.as_u64()).unwrap_or(0),
            length: doc.get("slen").and_then(|v| v.as_u64()).unwrap_or(0) as usize,
            gene_annotation: doc
                .get("title")
                .and_then(|v| v.as_str())
                .unwrap_or_default()
                .to_string(),
            refseq: sourcedb.eq_ignore_ascii_case("refseq"),
            source_db: match sourcedb {
                "insd" => "GenBank".to_string(),
                "refseq" => "RefSeq".to_string(),
                other => other.to_string(),
            },
            query_taxid,
            taxon_group: group,
            taxonomic_outlier: false,
        });
        parsed += 1;
    }
    Ok(parsed)
}

/// Human-readable summary to stderr (data goes to the JSON file). Keeps the
/// command scriptable: stdout stays clean for piping.
fn print_summary(results: &[QueryResult], out_path: &Path) {
    eprintln!();
    eprintln!("query complete");
    for result in results {
        let distinct_taxa: HashSet<u64> = result.accessions.iter().map(|a| a.taxid).collect();
        let refseq = result.accessions.iter().filter(|a| a.refseq).count();
        let group = match result.taxon_group {
            TaxonGroup::Ingroup => "ingroup",
            TaxonGroup::Outgroup => "outgroup",
        };
        eprintln!(
            "  [{group}] {} ({}): {} accessions, {} distinct taxa, {refseq} RefSeq",
            result.taxon_name,
            result.taxid,
            result.total_accessions,
            distinct_taxa.len()
        );
    }
    let total: usize = results.iter().map(|r| r.total_accessions).sum();
    eprintln!("  total accessions: {total}");
    eprintln!("  written to:       {}", out_path.display());
}

/// Warn — never drop — when the same accession is returned by both an ingroup and
/// an outgroup query. A sequence cannot honestly be both, so this almost always
/// means overlapping or mis-chosen TaxIDs the user should know about. Resolving
/// it (ingroup wins) happens in fetch's preflight.
fn warn_cross_group_overlap(results: &[QueryResult]) {
    let mut ingroup: HashSet<&str> = HashSet::new();
    let mut outgroup: HashSet<&str> = HashSet::new();
    for result in results {
        let set = match result.taxon_group {
            TaxonGroup::Ingroup => &mut ingroup,
            TaxonGroup::Outgroup => &mut outgroup,
        };
        set.extend(result.accessions.iter().map(|a| a.accession.as_str()));
    }

    let mut overlap: Vec<&str> = ingroup.intersection(&outgroup).copied().collect();
    if overlap.is_empty() {
        return;
    }
    overlap.sort_unstable();

    eprintln!(
        "warning: {} accession(s) returned by both ingroup and outgroup queries; \
         check for overlapping TaxIDs (fetch will resolve these, ingroup wins):",
        overlap.len()
    );
    for accession in &overlap {
        eprintln!("  {accession}");
    }
}
