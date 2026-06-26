use phorge::parse_fasta;
use clap::Args;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

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
    #[arg(short = 'n', long)]
    pub min_loci: Option<usize>,
}

pub fn run(args: FilterArgs) {
    let (sequences, _) = parse_fasta(&args.input, true)
        .expect("Could not read alignment (file not found or sequences are not the same length)");

    // Parse the provenance TSV once before the loop
    let loci_counts: HashMap<String, usize> = if let Some(log) = &args.log {
        let file = File::open(log).expect("Could not open provenance TSV");
        let reader = BufReader::new(file);
        let mut counts = HashMap::new();
        for line in reader.lines().skip(1) {
            let line = line.expect("Failed to read provenance TSV");
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }
            let taxon = fields[0].to_string();
            let count: usize = fields[1..].iter().filter(|&&f| f != "MISSING").count();
            counts.insert(taxon, count);
        }
        counts
    } else {
        HashMap::new()
    };

    let total = sequences.len();
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

        if let Some(min_loci) = args.min_loci
            && let Some(&count) = loci_counts.get(header)
            && count < min_loci
        {
            keep = false;
        }

        if keep {
            kept.push((header.clone(), seq.clone()));
        }
    }

    for (header, seq) in &kept {
        println!(">{}", header);
        println!("{}", seq);
    }

    eprintln!("Total taxa: {}", total);
    eprintln!("Kept taxa: {}", kept.len());
    eprintln!("Dropped taxa: {}", total - kept.len());
}
