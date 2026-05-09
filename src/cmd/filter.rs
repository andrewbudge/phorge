use cladekit::parse_fasta;
use clap::Args;

#[derive(Args)]
pub struct FilterArgs {
    /// Aligned FASTA file to filter
    pub input: String,

    /// Maximum allowed missingness per taxon (0.0–1.0)
    #[arg(short, long)]
    pub max_missing: Option<f64>,

    /// Provenance TSV from concat -l (required with --min-loci)
    #[arg(short, long)]
    pub log: Option<String>,

    /// Minimum number of loci a taxon must be present in
    #[arg(long)]
    pub min_loci: Option<usize>,
}

pub fn run(args: FilterArgs) {
    // TODO: parse the input FASTA (use parse_fasta from crate::)
    //   - if args.input is Some(path), parse that file
    //   - if None, we'll handle stdin later — for now just expect a file
    let (sequences, _) = parse_fasta(&args.input, false).expect("Could not open file");
    let mut kept = Vec::new();

    for (header, seq) in &sequences {
        let mut keep = true;
        if let Some(threshold) = args.max_missing {
            let mut missing = 0;
            for ch in seq.chars() {
                match ch {
                    '-' | 'N' | '?' | 'X' => missing += 1,
                    _ => {}
                }
            }
            let fraction = missing as f64 / seq.len() as f64;
            if fraction > threshold {
                keep = false;
            }
        }

        // TODO: if let Some(min_loci) = args.min_loci { ... check provenance TSV ... }

        if keep {
            kept.push((header.clone(), seq.clone()));
        }
    }

    // TODO: print kept sequences to stdout in FASTA format (>header\nsequence\n)

    // TODO: print a summary to stderr:
    //   - how many taxa were in the input
    //   - how many passed
    //   - how many were dropped
}
