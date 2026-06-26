use phorge::{is_dna, parse_fasta, print_table};
use clap::Args;
use std::collections::HashMap;
use std::path::Path;

#[derive(Args)]
pub struct StatsArgs {
    /// FASTA files to summarize (accepts multiple files and globs)
    pub files: Vec<String>,

    /// Per-sequence statistics (header, length, GC%, missingness)
    #[arg(short = 'd', long = "detailed")]
    pub detailed: bool,

    /// Column-aligned output for readability
    #[arg(short = 'p', long = "pretty")]
    pub pretty: bool,
}

struct SequenceStats {
    gc_count: usize,
    missing_count: usize,
    total_chrs: usize,
    length: usize,
}

/// Count variable and parsimony-informative sites in an alignment.
/// Works for both DNA and protein — counts any non-gap, non-unknown character.
/// Variable: column has at least 2 different residues (ignoring gaps/X/N/?).
/// Parsimony-informative: column has at least 2 different residues, each appearing at least twice.
fn count_informative_sites(sequences: &[(String, String)], length: usize) -> (usize, usize) {
    let mut variable = 0;
    let mut informative = 0;
    let seqs: Vec<&String> = sequences.iter().map(|(_, s)| s).collect();

    for i in 0..length {
        let mut counts: HashMap<char, usize> = HashMap::new();
        for seq in &seqs {
            let ch = seq.as_bytes()[i] as char;
            match ch {
                '-' | 'N' | 'X' | '?' => {} // ignore gaps and unknowns
                _ => *counts.entry(ch).or_insert(0) += 1,
            }
        }

        let bases_present = counts.len();
        let bases_twice = counts.values().filter(|&&c| c >= 2).count();

        if bases_present > 1 {
            variable += 1;
        }
        if bases_twice >= 2 {
            informative += 1;
        }
    }

    (variable, informative)
}

// Helper func to calc individual seq data
fn calc_seq_stats(seq: &str) -> SequenceStats {
    let mut gc_count = 0;
    let mut missing_count = 0;
    let mut total_chrs = 0;
    for ch in seq.chars() {
        total_chrs += 1;
        match ch {
            'G' | 'C' => gc_count += 1,
            '-' | 'N' | 'X' | '?' => missing_count += 1,
            _ => {}
        }
    }

    SequenceStats {
        gc_count,
        missing_count,
        total_chrs,
        length: seq.len(),
    }
}

pub fn run(args: StatsArgs) {
    let mut rows: Vec<Vec<String>> = Vec::new();

    if args.detailed {
        rows.push(vec![
            "file".into(),
            "header".into(),
            "length".into(),
            "gc_pct".into(),
            "missing_pct".into(),
        ]);
        for file in &args.files {
            let (sequences, _) = parse_fasta(file, false).unwrap_or_else(|e| {
                eprintln!("Error: failed to parse '{}': {}", file, e);
                std::process::exit(1);
            });
            let dna = is_dna(&sequences);
            let filename = Path::new(file).file_name().unwrap().to_str().unwrap();

            for (header, seq) in &sequences {
                let stats = calc_seq_stats(seq);
                let missing_pct = stats.missing_count as f64 / stats.total_chrs as f64 * 100.0;
                let gc_str = if dna {
                    let gc_pct = stats.gc_count as f64
                        / (stats.total_chrs - stats.missing_count) as f64
                        * 100.0;
                    format!("{:.1}", gc_pct)
                } else {
                    "NA".to_string()
                };
                rows.push(vec![
                    filename.to_string(),
                    header.to_string(),
                    stats.length.to_string(),
                    gc_str,
                    format!("{:.1}", missing_pct),
                ]);
            }
        }
    } else {
        rows.push(vec![
            "file".into(),
            "sequences".into(),
            "length".into(),
            "type".into(),
            "gc_pct".into(),
            "missing_pct".into(),
            "variable".into(),
            "variable_pct".into(),
            "informative".into(),
            "informative_pct".into(),
        ]);
        for file in &args.files {
            let (sequences, length) = parse_fasta(file, false).unwrap_or_else(|e| {
                eprintln!("Error: failed to parse '{}': {}", file, e);
                std::process::exit(1);
            });
            let num_sequences = sequences.len();
            let dna = is_dna(&sequences);

            let mut gc_count = 0;
            let mut missing_count = 0;
            let mut total_chars = 0;
            for (_, seq) in &sequences {
                let stats = calc_seq_stats(seq);
                gc_count += stats.gc_count;
                missing_count += stats.missing_count;
                total_chars += stats.total_chrs;
            }

            let missing_pct = missing_count as f64 / total_chars as f64 * 100.0;
            let gc_str = if dna {
                let gc_pct = gc_count as f64 / (total_chars - missing_count) as f64 * 100.0;
                format!("{:.1}", gc_pct)
            } else {
                "NA".to_string()
            };

            let data_type = if dna { "DNA" } else { "AA" };
            let all_equal = sequences.iter().all(|(_, s)| s.len() == length);
            let filename = Path::new(file).file_name().unwrap().to_str().unwrap();

            if all_equal {
                let (variable, informative) = count_informative_sites(&sequences, length);
                let variable_pct = variable as f64 / length as f64 * 100.0;
                let informative_pct = informative as f64 / length as f64 * 100.0;
                rows.push(vec![
                    filename.to_string(),
                    num_sequences.to_string(),
                    length.to_string(),
                    data_type.to_string(),
                    gc_str,
                    format!("{:.1}", missing_pct),
                    variable.to_string(),
                    format!("{:.1}", variable_pct),
                    informative.to_string(),
                    format!("{:.1}", informative_pct),
                ]);
            } else {
                rows.push(vec![
                    filename.to_string(),
                    num_sequences.to_string(),
                    "NA".to_string(),
                    data_type.to_string(),
                    gc_str,
                    format!("{:.1}", missing_pct),
                    "NA".to_string(),
                    "NA".to_string(),
                    "NA".to_string(),
                    "NA".to_string(),
                ]);
            }
        }
    }

    print_table(&rows, args.pretty);
}
