use phorge::{load_taxa_list, parse_fasta};
use clap::Args;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Args)]
pub struct ConcatArgs {
    /// FASTA alignment files (accepts multiple files and globs)
    pub files: Vec<String>,

    /// Alias list for smart matching — output names that map to messy input headers
    #[arg(short, long)]
    pub alias: Option<String>,

    /// Output format: fasta (default) or nexus
    #[arg(short, long, default_value = "FASTA")]
    pub format: String,

    /// Override missing data character (default: N for DNA, X for amino acid)
    #[arg(short, long)]
    pub missing: Option<String>,

    /// Partition format: raxml (default) or nexus
    #[arg(short = 'p', long = "partitions", default_value = "raxml")]
    pub partitions: String,

    /// Provenance TSV output file (required with -a)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,

    /// Dry run: show matching summary without building the supermatrix
    #[arg(long)]
    pub dry_run: bool,
}

/// Per-gene match result: taxon name -> (original FASTA header, sequence).
/// The header is kept alongside the sequence so we can emit the provenance TSV.
type MatchedGene = HashMap<String, (String, String)>;

/// Match taxa names to FASTA headers using case-insensitive substring search.
/// Longer taxa names match first to prevent partial match collisions
/// (e.g., "Mus musculus domesticus" claims before "Mus musculus").
/// Once a header is claimed, no other taxon can match it.
/// Returns taxon -> (original header, sequence) so we can track provenance.
fn match_taxa(taxa: &[String], sequences: &[(String, String)]) -> MatchedGene {
    let mut sorted_taxa = taxa.to_vec();
    sorted_taxa.sort_by_key(|t| std::cmp::Reverse(t.len()));

    let mut claimed_headers = HashSet::new();
    let mut results = HashMap::new();

    for taxon in &sorted_taxa {
        for (header, sequence) in sequences {
            if claimed_headers.contains(header) {
                continue;
            }
            if header.to_lowercase().contains(&taxon.to_lowercase()) {
                claimed_headers.insert(header.clone());
                results.insert(taxon.clone(), (header.clone(), sequence.clone()));
                break;
            }
        }
    }
    results
}

enum DataType {
    Dna,
    AminoAcid,
}

fn detect_data_type(sequence: &str) -> DataType {
    if sequence
        .chars()
        .any(|c| matches!(c, 'E' | 'F' | 'I' | 'L' | 'P' | 'Q'))
    {
        DataType::AminoAcid
    } else {
        DataType::Dna
    }
}

pub fn run(args: ConcatArgs) {
    // Validate the alias file early, before touching gene files
    if let Some(ref taxa_file) = args.alias {
        let alias_path = Path::new(taxa_file);
        if alias_path.is_dir() {
            eprintln!(
                "Error: '{}' is a directory. -a/--alias expects a plain text file with one taxon name per line.",
                taxa_file
            );
            std::process::exit(1);
        }
        if let Ok(f) = File::open(alias_path) {
            let mut first_line = String::new();
            BufReader::new(f).read_line(&mut first_line).ok();
            if first_line.trim_start().starts_with('>') {
                eprintln!(
                    "Error: '{}' looks like a FASTA file. -a/--alias expects a plain text taxa list (one name per line, no '>' headers).",
                    taxa_file
                );
                std::process::exit(1);
            }
        }
    }

    // Parse all gene files
    let mut gene_data = Vec::new();
    for file in &args.files {
        let (sequences, length) = parse_fasta(file, true).expect("Failed to parse fasta file");
        gene_data.push((file, sequences, length));
    }

    let smart_matching = args.alias.is_some();

    // Determine taxa list and build matched_genes
    let taxa: Vec<String>;
    let mut matched_genes: Vec<(&String, MatchedGene, &usize)> = Vec::new();

    if let Some(ref taxa_file) = args.alias {
        // Smart matching mode: load taxa list, match by substring
        taxa = load_taxa_list(taxa_file).expect("Failed to load taxa list");
        for (file, sequences, length) in &gene_data {
            let matched = match_taxa(&taxa, sequences);
            matched_genes.push((file, matched, length));
        }
    } else {
        // Exact match mode: union of all headers across files becomes the taxa list
        let mut seen = HashSet::new();
        let mut taxa_order = Vec::new();
        for (_file, sequences, _length) in &gene_data {
            for (header, _) in sequences {
                if seen.insert(header.clone()) {
                    taxa_order.push(header.clone());
                }
            }
        }
        taxa = taxa_order;
        // Direct lookup — header IS the taxon name, no matching needed
        for (file, sequences, length) in &gene_data {
            let mut matched = HashMap::new();
            for (header, sequence) in sequences {
                matched.insert(header.clone(), (header.clone(), sequence.clone()));
            }
            matched_genes.push((file, matched, length));
        }
    }

    // Non-dry-run smart matching requires a log file to capture provenance
    if smart_matching && !args.dry_run && args.log.is_none() {
        eprintln!("Error: -a/--alias requires -l/--log to write the provenance TSV.");
        eprintln!("  To preview matching without writing output, use --dry-run.");
        std::process::exit(1);
    }

    // Dry run: show matching summary and output tentative provenance TSV, then exit
    if args.dry_run {
        eprintln!("=== Dry Run: Matching Summary ===\n");

        // Per-gene match counts
        eprintln!(
            "{:<40} {:>8} {:>8} {:>8}",
            "Gene", "Seqs", "Matched", "Missing"
        );
        eprintln!("{}", "-".repeat(68));
        for (file, matched, _length) in &matched_genes {
            let name = Path::new(file).file_name().unwrap().to_str().unwrap();
            let matched_count = matched.len();
            let missing_count = taxa.len() - matched_count;
            let seq_count = gene_data
                .iter()
                .find(|(f, _, _)| f == file)
                .map(|(_, s, _)| s.len())
                .unwrap_or(0);
            eprintln!(
                "{:<40} {:>8} {:>8} {:>8}",
                name, seq_count, matched_count, missing_count
            );
        }

        // Per-taxon coverage
        eprintln!("\n{:<40} {:>8}", "Taxon", "Genes");
        eprintln!("{}", "-".repeat(50));
        let mut taxa_coverage: Vec<(&String, usize)> = taxa
            .iter()
            .map(|taxon| {
                let count = matched_genes
                    .iter()
                    .filter(|(_, matched, _)| matched.contains_key(taxon))
                    .count();
                (taxon, count)
            })
            .collect();
        taxa_coverage.sort_by(|a, b| b.1.cmp(&a.1));

        for (taxon, count) in &taxa_coverage {
            eprintln!("{:<40} {:>8}/{}", taxon, count, matched_genes.len());
        }

        // Summary
        let total_genes = matched_genes.len();
        let full_coverage = taxa_coverage
            .iter()
            .filter(|(_, c)| *c == total_genes)
            .count();
        let any_coverage = taxa_coverage.iter().filter(|(_, c)| *c > 0).count();
        let no_coverage = taxa_coverage.iter().filter(|(_, c)| *c == 0).count();
        eprintln!("\n=== Summary ===");
        eprintln!(
            "Taxa: {} total, {} with all genes, {} with some, {} with none",
            taxa.len(),
            full_coverage,
            any_coverage,
            no_coverage
        );
        eprintln!("Genes: {}", total_genes);

        // Output tentative provenance TSV (only meaningful in smart matching mode)
        if smart_matching {
            let taxa_file = args.alias.as_ref().unwrap();
            let gene_names: Vec<String> = matched_genes
                .iter()
                .map(|(file, _, _)| {
                    Path::new(file)
                        .file_name()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .to_string()
                })
                .collect();

            let mut rows: Vec<String> = Vec::new();
            rows.push(format!("{}\t{}", taxa_file, gene_names.join("\t")));
            for taxon in &taxa {
                let mut row = vec![taxon.clone()];
                for (_file, matched, _length) in &matched_genes {
                    if matched.contains_key(taxon) {
                        row.push(matched[taxon].0.clone());
                    } else {
                        row.push("MISSING".to_string());
                    }
                }
                rows.push(row.join("\t"));
            }

            if let Some(log_path) = &args.log {
                let mut f = File::create(log_path).expect("Failed to create provenance log file");
                for row in &rows {
                    writeln!(f, "{}", row).unwrap();
                }
                eprintln!("\nTentative provenance log written to: {}", log_path);
            } else {
                eprintln!("\n=== Tentative Provenance TSV ===");
                for row in &rows {
                    println!("{}", row);
                }
            }
        }

        return;
    }

    // Detect data type per gene using first available sequence
    let gene_types: Vec<DataType> = matched_genes
        .iter()
        .map(|(_, matched, _)| {
            matched
                .values()
                .next()
                .map(|(_, seq)| detect_data_type(seq))
                .unwrap_or(DataType::Dna)
        })
        .collect();

    // Determine overall data type mix
    let has_dna = gene_types.iter().any(|t| matches!(t, DataType::Dna));
    let has_aa = gene_types.iter().any(|t| matches!(t, DataType::AminoAcid));
    let is_mixed = has_dna && has_aa;

    // Build supermatrix: concatenate matched sequences per taxon, fill gaps with missing char
    let mut supermatrix: HashMap<String, String> = HashMap::new();
    for taxon in &taxa {
        for (i, (_file, matched, length)) in matched_genes.iter().enumerate() {
            let entry = supermatrix.entry(taxon.clone()).or_default();
            if matched.contains_key(taxon) {
                entry.push_str(&matched[taxon].1);
            } else {
                let missing_char = match &args.missing {
                    Some(m) => m.clone(),
                    None => {
                        if is_mixed {
                            "?".to_string()
                        } else {
                            match gene_types[i] {
                                DataType::Dna => "N".to_string(),
                                DataType::AminoAcid => "X".to_string(),
                            }
                        }
                    }
                };
                entry.push_str(&missing_char.repeat(**length));
            }
        }
    }

    // Build partition boundaries with data type (used by both output formats)
    let mut partitions = Vec::new();
    let mut position = 1;
    for (i, (file, _matched, length)) in matched_genes.iter().enumerate() {
        let name = Path::new(file).file_name().unwrap().to_str().unwrap();
        partitions.push((
            name.to_string(),
            position,
            position + *length - 1,
            &gene_types[i],
        ));
        position += *length;
    }
    let total_length = position - 1;

    // Determine overall data type label for NEXUS output
    let nexus_datatype = match (has_dna, has_aa) {
        (true, true) => "MIXED",
        (false, true) => "PROTEIN",
        _ => "DNA",
    };

    // Default missing char for NEXUS format line
    let nexus_missing = match &args.missing {
        Some(m) => m.clone(),
        None => match (has_dna, has_aa) {
            (true, true) => "?".to_string(),
            (false, true) => "X".to_string(),
            _ => "N".to_string(),
        },
    };

    let fmt = args.format.to_lowercase();
    if fmt == "nexus" || fmt == "n" || fmt == "nex" {
        // NEXUS: complete file to stdout (data + partitions in one)
        println!("#NEXUS");
        println!("BEGIN DATA;");
        println!("  DIMENSIONS NTAX={} NCHAR={};", taxa.len(), total_length);
        println!(
            "  FORMAT DATATYPE={} MISSING={} GAP=-;",
            nexus_datatype, nexus_missing
        );
        println!("  MATRIX");
        for taxon in &taxa {
            println!("  {}    {}", taxon, supermatrix[taxon]);
        }
        println!(";");
        println!("END;");
        println!("BEGIN SETS;");
        for (name, start, end, _) in &partitions {
            println!("  CHARSET {} = {}-{};", name, start, end);
        }
        println!("END;");
    } else {
        // FASTA: supermatrix to stdout
        for taxon in &taxa {
            println!(">{}", taxon);
            println!("{}", supermatrix[taxon]);
        }
    }

    // Partitions to stderr
    let part_fmt = args.partitions.to_lowercase();
    if part_fmt == "nexus" || part_fmt == "n" || part_fmt == "nex" {
        for (name, start, end, _) in &partitions {
            eprintln!("CHARSET {} = {}-{};", name, start, end);
        }
    } else {
        // RAxML/IQ-TREE format (default)
        for (name, start, end, dtype) in &partitions {
            let model = match dtype {
                DataType::Dna => "DNA",
                DataType::AminoAcid => "WAG",
            };
            eprintln!("{}, {} = {}-{}", model, name, start, end);
        }
    }

    // Write provenance TSV (only in smart matching mode, -l is required with -t)
    if smart_matching {
        let log_path = args.log.as_ref().unwrap();
        let taxa_file = args.alias.as_ref().unwrap();
        let mut log_file = File::create(log_path).expect("Failed to create provenance log file");
        // Header row: taxa list filename, then each gene filename
        let gene_names: Vec<String> = matched_genes
            .iter()
            .map(|(file, _, _)| {
                Path::new(file)
                    .file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .to_string()
            })
            .collect();
        writeln!(log_file, "{}\t{}", taxa_file, gene_names.join("\t"))
            .expect("Failed to write to log file");
        // One row per taxon: taxon name, then matched header or MISSING
        for taxon in &taxa {
            let mut row = vec![taxon.clone()];
            for (_file, matched, _length) in &matched_genes {
                if matched.contains_key(taxon) {
                    row.push(matched[taxon].0.clone());
                } else {
                    row.push("MISSING".to_string());
                }
            }
            writeln!(log_file, "{}", row.join("\t")).expect("Failed to write to log file");
        }
    }
}
