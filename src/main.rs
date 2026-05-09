use clap::{Parser, Subcommand};

// Exit silently on broken pipe (e.g., piping to head/tail)
#[cfg(unix)]
fn reset_sigpipe() {
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

mod cmd;

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
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

    /// Extract gene regions from target sequences using homology search
    Extract(cmd::extract::ExtractArgs),

    /// Remove taxa exceeding a missingness threshold from an alignment
    Filter(cmd::filter::FilterArgs),

    // Cluster(cmd::cluster::ClusterArgs), // WIP: metagenomics use case, deprioritized
}

fn main() {
    #[cfg(unix)]
    reset_sigpipe();

    let cli = Cli::parse();
    match cli.command {
        Commands::Getheaders(args) => cmd::getheaders::run(args),
        Commands::Concat(args) => cmd::concat::run(args),
        Commands::Stats(args) => cmd::stats::run(args),
        Commands::Coverage(args) => cmd::coverage::run(args),
        Commands::Convert(args) => cmd::convert::run(args),
        Commands::Align(args) => cmd::align::run(args),
        Commands::Extract(args) => cmd::extract::run(args),
        Commands::Filter(args) => cmd::filter::run(args),
        // Commands::Cluster(args) => cmd::cluster::run(args), // WIP: deprioritized
    }
}
