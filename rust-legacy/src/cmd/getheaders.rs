use clap::Args;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Args)]
pub struct GetheadersArgs {
    /// Input FASTA files (accepts multiple files and globs); reads stdin if none given
    pub input: Vec<String>,

    /// Print only unique headers (deduplicate)
    #[arg(short, long)]
    pub unique: bool,
}

pub fn run(args: GetheadersArgs) {
    let mut seen = HashSet::new();

    if args.input.is_empty() {
        process_reader(
            BufReader::new(std::io::stdin().lock()),
            args.unique,
            &mut seen,
        );
    } else {
        for filename in &args.input {
            let file = File::open(filename).unwrap_or_else(|e| {
                eprintln!("Error: could not open '{}': {}", filename, e);
                std::process::exit(1);
            });
            process_reader(BufReader::new(file), args.unique, &mut seen);
        }
    }
}

fn process_reader(reader: impl BufRead, unique: bool, seen: &mut HashSet<String>) {
    for line in reader.lines() {
        let line = line.expect("Could not read line");
        if let Some(header) = line.strip_prefix('>') {
            if unique {
                if seen.insert(header.to_string()) {
                    println!("{}", header);
                }
            } else {
                println!("{}", header);
            }
        }
    }
}
