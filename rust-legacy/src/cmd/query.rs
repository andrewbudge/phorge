use crate::models::{Accession, TaxonGroup, render_query_tsv, write_query_tsv};
use crate::ncbi::EutilsClient;
use anyhow::{Context, Result, bail};
use clap::Args;
use serde_json::Value;
use std::collections::HashSet;
use std::io::{self, Write};
use std::path::PathBuf;

#[derive(Args)]
pub struct QueryArgs {
    /// Ingroup taxa: TaxIDs (e.g. 7088) or names (e.g. Felidae). Names are
    /// resolved against NCBI taxonomy and confirmed. Comma-separated, or repeat
    /// the flag.
    #[arg(long, short = 'i', value_delimiter = ',', required = true)]
    pub ingroup: Vec<String>,

    /// Outgroup taxa: TaxIDs or names, same rules as --ingroup.
    #[arg(long, short = 'o', value_delimiter = ',')]
    pub outgroup: Vec<String>,

    /// Search term(s) restricting results to loci/genetic data of interest
    /// (e.g. COX1,12S). OR'd together, then AND'd onto each organism query.
    #[arg(long, short = 't', value_delimiter = ',')]
    pub term: Vec<String>,

    /// Output query file (TSV). Omit, or pass `-`, to write to stdout instead
    /// (pipe into `awk`/`fetch`; use `-y` when piping non-interactively).
    #[arg(long, short = 'q')]
    pub query: Option<PathBuf>,

    /// Email address required by NCBI ToS for automated access.
    #[arg(long, short = 'e')]
    pub email: String,

    /// NCBI API key (optional; raises the NCBI rate limit from 3 to 10 req/s).
    #[arg(long)]
    pub api_key: Option<String>,

    /// Skip name-confirmation prompts (non-interactive). Aborts on an ambiguous
    /// name rather than guessing.
    #[arg(long, short = 'y')]
    pub yes: bool,
}

/// Nucleotide database to search. `nuccore` is the canonical name for the
/// GenBank/RefSeq nucleotide set.
const DB: &str = "nuccore";

/// esummary docsums fetched per request. NCBI tolerates more, but 500 keeps each
/// response modest and matches the fetch-stage batch size.
const PAGE_SIZE: usize = 500;

/// One esearch sweep over a single root taxon. The output TSV is a flat table,
/// but keeping results grouped in memory drives the per-taxon summary and the
/// cross-group overlap check before they are flattened for writing.
struct TaxonQuery {
    taxid: u64,
    taxon_name: String,
    taxon_group: TaxonGroup,
    total_accessions: usize,
    accessions: Vec<Accession>,
}

pub async fn run(args: QueryArgs) -> Result<()> {
    let client = EutilsClient::new(args.api_key, args.email).context("building NCBI client")?;

    // Resolve every ingroup/outgroup entry to a concrete TaxID first (prompting
    // for names), so a typo aborts before any nuccore querying begins.
    let mut roots: Vec<(u64, TaxonGroup)> = Vec::new();
    for raw in &args.ingroup {
        roots.push((
            resolve_input(&client, raw, args.yes).await?,
            TaxonGroup::Ingroup,
        ));
    }
    for raw in &args.outgroup {
        roots.push((
            resolve_input(&client, raw, args.yes).await?,
            TaxonGroup::Outgroup,
        ));
    }

    // One TaxonQuery per root, ingroup first then outgroup. Both groups are
    // queried identically; the only difference is the TaxonGroup tag, which
    // downstream stages (e.g. fetch's ingroup-wins overlap rule) rely on.
    let mut queries: Vec<TaxonQuery> = Vec::with_capacity(roots.len());
    for (taxid, group) in &roots {
        queries.push(query_taxon(&client, *taxid, *group, &args.term).await?);
    }

    warn_cross_group_overlap(&queries);

    // A `-q -` is an explicit request for stdout; no path (or none at all) also
    // means stdout. Only a real file path routes to disk.
    let file_path = args.query.filter(|p| p.as_os_str() != "-");
    print_summary(&queries, file_path.as_deref());

    let records: Vec<Accession> = queries.into_iter().flat_map(|q| q.accessions).collect();
    match &file_path {
        Some(path) => {
            // Ensure the parent directory exists, then write one row per accession.
            if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating output directory {}", parent.display()))?;
            }
            write_query_tsv(path, &records)?;
        }
        None => {
            // Data to stdout (all diagnostics already went to stderr), so the
            // table can be piped straight into `awk`/`fetch`.
            let tsv = render_query_tsv(&records);
            io::stdout()
                .write_all(tsv.as_bytes())
                .context("writing query TSV to stdout")?;
        }
    }
    Ok(())
}

/// Turn one user-supplied ingroup/outgroup entry into a concrete TaxID. A purely
/// numeric entry is taken as a TaxID verbatim (no lookup, no prompt); anything
/// else is a name, resolved against NCBI taxonomy and confirmed with the user
/// unless `--yes` is set.
async fn resolve_input(client: &EutilsClient, raw: &str, yes: bool) -> Result<u64> {
    if let Ok(taxid) = raw.parse::<u64>() {
        return Ok(taxid);
    }

    let matches = client
        .resolve_taxon_name(raw)
        .await
        .with_context(|| format!("resolving taxon name {raw:?}"))?;

    match matches.as_slice() {
        [] => bail!("no NCBI taxon found for name {raw:?}; check spelling or pass a TaxID"),
        [only] => {
            eprintln!(
                "  {raw:?} -> {} (txid{}), {} in {}",
                only.name, only.taxid, only.rank, only.division
            );
            if yes || confirm("  use this taxon?")? {
                Ok(only.taxid)
            } else {
                bail!("aborted: {raw:?} not confirmed")
            }
        }
        many => {
            // Homonym: the same name across divisions. Non-interactive can't guess.
            if yes {
                bail!(
                    "name {raw:?} is ambiguous ({} matches); pass a TaxID instead",
                    many.len()
                );
            }
            eprintln!("  {raw:?} is ambiguous — {} matches:", many.len());
            for (i, m) in many.iter().enumerate() {
                eprintln!(
                    "    {}) {} (txid{}), {} in {}",
                    i + 1,
                    m.name,
                    m.taxid,
                    m.rank,
                    m.division
                );
            }
            Ok(many[pick(many.len())?].taxid)
        }
    }
}

/// Ask a yes/no question on stderr, read the answer from stdin. Anything but an
/// explicit yes is No (safe default). A closed stdin (piped/non-interactive)
/// aborts rather than silently answering, so scripts must pass `--yes`.
fn confirm(prompt: &str) -> Result<bool> {
    eprint!("{prompt} [y/N] ");
    io::stderr().flush()?;
    let mut input = String::new();
    if io::stdin().read_line(&mut input)? == 0 {
        bail!("no input on stdin; re-run with --yes or supply a TaxID");
    }
    Ok(matches!(input.trim(), "y" | "Y" | "yes" | "Yes"))
}

/// Prompt for a 1-based selection in `[1, n]`, returning the 0-based index.
/// Re-prompts on garbage; `q` or EOF aborts.
fn pick(n: usize) -> Result<usize> {
    loop {
        eprint!("  select [1-{n}] (or q to abort): ");
        io::stderr().flush()?;
        let mut input = String::new();
        if io::stdin().read_line(&mut input)? == 0 {
            bail!("no input on stdin; re-run with --yes or supply a TaxID");
        }
        let s = input.trim();
        if s.eq_ignore_ascii_case("q") {
            bail!("aborted at taxon selection");
        }
        match s.parse::<usize>() {
            Ok(i) if (1..=n).contains(&i) => return Ok(i - 1),
            _ => eprintln!("  enter a number between 1 and {n}"),
        }
    }
}

/// Run one esearch + esummary sweep over a single taxon and collect its
/// accessions. `[Organism:exp]` excludes environmental samples and expands the
/// TaxID to its full subtree. Optional `terms` restrict the pull to loci of
/// interest — OR'd together and AND'd onto the organism clause.
async fn query_taxon(
    client: &EutilsClient,
    taxid: u64,
    group: TaxonGroup,
    terms: &[String],
) -> Result<TaxonQuery> {
    let mut term = format!("txid{taxid}[Organism:exp]");
    if !terms.is_empty() {
        let clause = terms
            .iter()
            .map(|t| format!("{t}[All Fields]"))
            .collect::<Vec<_>>()
            .join(" OR ");
        term.push_str(&format!(" AND ({clause})"));
    }
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

    Ok(TaxonQuery {
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

/// Human-readable summary to stderr (data goes to the TSV sink). Keeps the
/// command scriptable: stdout stays clean for piping.
fn print_summary(queries: &[TaxonQuery], out_path: Option<&std::path::Path>) {
    eprintln!();
    eprintln!("query complete");
    for q in queries {
        let distinct_taxa: HashSet<u64> = q.accessions.iter().map(|a| a.taxid).collect();
        let refseq = q.accessions.iter().filter(|a| a.refseq).count();
        eprintln!(
            "  [{}] {} ({}): {} accessions, {} distinct taxa, {refseq} RefSeq",
            q.taxon_group.as_str(),
            q.taxon_name,
            q.taxid,
            q.total_accessions,
            distinct_taxa.len()
        );
    }
    let total: usize = queries.iter().map(|q| q.total_accessions).sum();
    eprintln!("  total accessions: {total}");
    match out_path {
        Some(path) => eprintln!("  written to:       {}", path.display()),
        None => eprintln!("  written to:       stdout"),
    }
}

/// Warn — never drop — when the same accession is returned by both an ingroup and
/// an outgroup query. A sequence cannot honestly be both, so this almost always
/// means overlapping or mis-chosen TaxIDs the user should know about. Resolving
/// it (ingroup wins) happens in fetch's preflight.
fn warn_cross_group_overlap(queries: &[TaxonQuery]) {
    let mut ingroup: HashSet<&str> = HashSet::new();
    let mut outgroup: HashSet<&str> = HashSet::new();
    for q in queries {
        let set = match q.taxon_group {
            TaxonGroup::Ingroup => &mut ingroup,
            TaxonGroup::Outgroup => &mut outgroup,
        };
        set.extend(q.accessions.iter().map(|a| a.accession.as_str()));
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
