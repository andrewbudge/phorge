use phorge::{is_dna, parse_fasta};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use clap::Args;

#[derive(Args)]
pub struct ConvertArgs {
    /// Input sequence file (format auto-detected from contents: FASTA, NEXUS, PHYLIP)
    pub input_file: String,

    /// Output format: f (fasta), n (nexus), sp (strict phylip), rp (relaxed phylip); writes to stdout
    #[arg(short = 'o', long = "output_format")]
    pub output_format: String,
}

/// Parse a NEXUS file into a list of (name, sequence) pairs in file order.
/// Finds the MATRIX block and reads taxon/sequence pairs until `;`.
fn parse_nexus(lines: impl Iterator<Item = std::io::Result<String>>) -> Vec<(String, String)> {
    let mut seqs = Vec::new();
    let mut in_matrix = false;

    for line in lines {
        let line = line.expect("Failed to read line");
        let trimmed = line.trim().to_string();

        if trimmed.is_empty() {
            continue;
        }

        if trimmed.to_uppercase().starts_with("MATRIX") {
            in_matrix = true;
            continue;
        }

        if in_matrix {
            if trimmed.starts_with(';') {
                break;
            }
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() >= 2 {
                seqs.push((fields[0].to_string(), fields[1].to_uppercase()));
            }
        }
    }
    seqs
}

/// Write sequences as FASTA to stdout.
fn write_fasta(sequences: &[(String, String)]) {
    for (name, seq) in sequences {
        println!(">{}", name);
        println!("{}", seq);
    }
}

/// Write sequences as NEXUS to stdout.
fn write_nexus(sequences: &[(String, String)]) {
    let dna = is_dna(sequences);
    let datatype = if dna { "DNA" } else { "PROTEIN" };
    let missing = if dna { "N" } else { "X" };
    let length = sequences.first().map_or(0, |(_, s)| s.len());

    println!("#NEXUS");
    println!("BEGIN DATA;");
    println!("  DIMENSIONS NTAX={} NCHAR={};", sequences.len(), length);
    println!("  FORMAT DATATYPE={} MISSING={} GAP=-;", datatype, missing);
    println!("  MATRIX");
    for (name, seq) in sequences {
        println!("  {}    {}", name, seq);
    }
    println!(";");
    println!("END;");
}

/// Write sequences as relaxed PHYLIP to stdout.
fn write_relaxed_phylip(sequences: &[(String, String)]) {
    let length = sequences.first().map_or(0, |(_, s)| s.len());
    println!("{} {}", sequences.len(), length);
    for (name, seq) in sequences {
        println!("{}    {}", name, seq);
    }
}

/// Write sequences as strict PHYLIP to stdout.
/// Taxon names are padded/truncated to exactly 10 characters.
fn write_strict_phylip(sequences: &[(String, String)]) {
    let length = sequences.first().map_or(0, |(_, s)| s.len());
    println!("{} {}", sequences.len(), length);
    for (name, seq) in sequences {
        // Truncate to 10 chars by character (not byte) so non-ASCII names can't panic.
        let label: String = name.chars().take(10).collect();
        print!("{:<10}", label);
        println!("{}", seq);
    }
}

pub fn run(args: ConvertArgs) {
    let file = File::open(&args.input_file).expect("Failed to open file");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let first_line = lines
        .next()
        .expect("empty file")
        .expect("Failed to read line");

    let sequences = if first_line.starts_with('>') {
        let (seqs, _) = parse_fasta(&args.input_file, false).expect("Failed to parse fasta file");
        seqs
    } else if first_line.starts_with('#') {
        parse_nexus(lines)
    } else if first_line
        .chars()
        .next()
        .expect("empty line")
        .is_ascii_digit()
    {
        let mut seqs = Vec::new();
        for line in lines {
            let line = line.expect("Failed to read line");
            let trimmed = line.trim().to_string();
            if trimmed.is_empty() {
                continue;
            }
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() >= 2 {
                seqs.push((fields[0].to_string(), fields[1].to_uppercase()));
            }
        }
        seqs
    } else {
        eprintln!("Error: could not detect input format");
        std::process::exit(1);
    };

    match args.output_format.as_str() {
        "f" | "fasta" => write_fasta(&sequences),
        "n" | "nexus" | "nex" => write_nexus(&sequences),
        "rp" | "phylip" => write_relaxed_phylip(&sequences),
        "sp" => write_strict_phylip(&sequences),
        _ => {
            eprintln!(
                "Error: unknown output format '{}'. Use: f, n, sp, or rp",
                args.output_format
            );
            std::process::exit(1);
        }
    }
}
