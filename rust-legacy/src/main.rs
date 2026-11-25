use clap::{Parser, Subcommand};
use std::time::{SystemTime, UNIX_EPOCH};
use tracing::Instrument;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;

// Exit silently on broken pipe (e.g., piping to head/tail)
#[cfg(unix)]
fn reset_sigpipe() {
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

mod cmd;
mod models;
mod ncbi;

#[derive(Parser)]
#[command(name = "phorge", version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    // --- Lean file tools (sync, no network) ---
    /// Extract headers from FASTA files
    #[command(visible_alias = "ghd")]
    Getheaders(cmd::getheaders::GetheadersArgs),

    /// Concatenate alignments into a supermatrix
    #[command(visible_alias = "liger")]
    Concat(cmd::concat::ConcatArgs),

    /// Alignment summary statistics
    Stats(cmd::stats::StatsArgs),

    /// Summarize taxa coverage from a concat provenance TSV
    Coverage(cmd::coverage::CoverageArgs),

    /// Convert between common sequence data file types
    Convert(cmd::convert::ConvertArgs),

    /// Batch align FASTA files using an external alignment program
    #[command(visible_alias = "aln")]
    Align(cmd::align::AlignArgs),

    /// Trim alignment columns by parsimony-informativeness and gappiness
    Curate(cmd::curate::CurateArgs),

    /// Extract gene regions from target sequences using homology search
    Extract(cmd::extract::ExtractArgs),

    /// Remove taxa exceeding a missingness threshold from an alignment
    Filter(cmd::filter::FilterArgs),

    // --- Acquisition layer (async, talks to NCBI / orchestrates MMseqs2) ---
    /// Retrieve sequence metadata from NCBI for a taxonomic group
    Query(cmd::query::QueryArgs),

    /// Download sequences for a query result set
    Fetch(cmd::fetch::FetchArgs),

    /// Standardize headers and write pipeline-ready FASTAs
    Clean(cmd::clean::CleanArgs),
}

fn command_name(command: &Commands) -> &'static str {
    match command {
        Commands::Getheaders(_) => "getheaders",
        Commands::Concat(_) => "concat",
        Commands::Stats(_) => "stats",
        Commands::Coverage(_) => "coverage",
        Commands::Convert(_) => "convert",
        Commands::Align(_) => "align",
        Commands::Curate(_) => "curate",
        Commands::Extract(_) => "extract",
        Commands::Filter(_) => "filter",
        Commands::Query(_) => "query",
        Commands::Fetch(_) => "fetch",
        Commands::Clean(_) => "clean",
    }
}

/// Single-layer tracing: human-readable narration on stderr. Every subcommand
/// logs the same way — capture it with `2>` if you want a record.
fn init_tracing() -> anyhow::Result<()> {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::registry()
        .with(filter)
        .with(
            tracing_subscriber::fmt::layer()
                .with_writer(std::io::stderr)
                .with_target(false),
        )
        .init();
    Ok(())
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    #[cfg(unix)]
    reset_sigpipe();

    let cli = Cli::parse();
    init_tracing()?;

    // Tag every event in this invocation with a run-id so a single run can be
    // sliced out of the append-only log.
    let name = command_name(&cli.command);
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let run_id = format!("{name}-{secs}");
    let span = tracing::info_span!("invocation", run_id = %run_id, command = name);

    async move {
        match cli.command {
            // Lean file tools: sync, handle their own errors/exit internally.
            Commands::Getheaders(args) => cmd::getheaders::run(args),
            Commands::Concat(args) => cmd::concat::run(args),
            Commands::Stats(args) => cmd::stats::run(args),
            Commands::Coverage(args) => cmd::coverage::run(args),
            Commands::Convert(args) => cmd::convert::run(args),
            Commands::Align(args) => cmd::align::run(args),
            Commands::Curate(args) => cmd::curate::run(args),
            Commands::Extract(args) => cmd::extract::run(args),
            Commands::Filter(args) => cmd::filter::run(args),
            // Acquisition layer: async, propagate errors via anyhow.
            Commands::Query(args) => return cmd::query::run(args).await,
            Commands::Fetch(args) => return cmd::fetch::run(args).await,
            Commands::Clean(args) => return cmd::clean::run(args).await,
        }
        Ok(())
    }
    .instrument(span)
    .await
}
