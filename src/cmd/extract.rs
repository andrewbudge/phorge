use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

use clap::Args;

use cladekit::parse_fasta;

#[derive(Args)]
pub struct ExtractArgs {
    /// Reference FASTA with labeled gene sequences (e.g., >COI, >ND2)
    #[arg(short, long)]
    pub reference: String,

    /// Target organism FASTA files or a directory containing them
    #[arg(short, long, required = true, num_args = 1..)]
    pub targets: Vec<String>,

    /// Output directory for per-gene FASTAs
    #[arg(short, long)]
    pub output: String,

    /// Extra bases to grab on either side of each hit
    #[arg(long, default_value_t = 0)]
    pub flank: usize,

    /// Minimum fraction of the reference gene that must be covered (0.0–1.0)
    #[arg(long, default_value_t = 0.5)]
    pub min_coverage: f64,

    /// MMseqs2 sensitivity (1.0=fast, 7.5=max). Higher finds more divergent hits.
    #[arg(short, long, default_value_t = 5.7)]
    pub sensitivity: f64,

    /// Keep intermediate files instead of deleting them after the search
    #[arg(long, default_value_t = false)]
    pub keep_intermediates: bool,
}

fn collect_targets(targets: &[String]) -> Vec<String> {
    let mut files: Vec<String> = Vec::new();

    for target in targets {
        let path = Path::new(target);
        if path.is_dir() {
            for entry in path.read_dir().unwrap() {
                let entry = entry.unwrap();
                let p = entry.path();
                if p.extension().map_or(false, |e| {
                    e == "fasta" || e == "fa" || e == "fna" || e == "fas"
                }) {
                    files.push(p.to_string_lossy().into_owned());
                }
            }
        } else if path.is_file() {
            files.push(target.clone());
        } else {
            eprintln!(
                "Warning: '{}' is not a file or directory, skipping.",
                target
            );
        }
    }

    files
}

// Writes a pooled FASTA to disk with "organism::seq_id" headers.
// Returns a lookup map from that key to (original_header, full_sequence, original_filename)
fn pool_targets(target_files: &[String], pooled_path: &Path) -> HashMap<String, (String, String, String)> {
    let mut writer = File::create(pooled_path).expect("Could not create pooled targets file");
    let mut lookup: HashMap<String, (String, String, String)> = HashMap::new();

    for file in target_files {
        let path = Path::new(file);
        let filename = path.file_name().unwrap().to_str().unwrap().to_string();
        let organism = path
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .replace(' ', "_");

        let (seqs, _) = parse_fasta(file, false).expect("Failed to read target FASTA");

        for (header, seq) in &seqs {
            // MMseqs2 truncates headers at the first whitespace, so we only
            // use the first token as the sequence ID in the pooled key.
            let seq_id = header.split_whitespace().next().unwrap_or(header);
            let pooled_key = format!("{}::{}", organism, seq_id);

            writeln!(writer, ">{}", pooled_key).unwrap();
            writeln!(writer, "{}", seq).unwrap();
            lookup.insert(pooled_key, (header.clone(), seq.clone(), filename.clone()));
        }
    }

    lookup
}

struct Hit {
    query: String,
    target: String,
    tstart: usize,
    tend: usize,
}

// Parses MMseqs2 tabular output. Filters hits by minimum query coverage.
// Expected --format-output: query,target,tstart,tend,qstart,qend,qlen
fn parse_hits(tsv_path: &Path, min_coverage: f64) -> Vec<Hit> {
    let file = File::open(tsv_path).expect("Could not open MMseqs2 output");
    let reader = BufReader::new(file);
    let mut hits = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Error reading MMseqs2 output");
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            continue;
        }

        let qstart: usize = f[4].parse().unwrap_or(0);
        let qend: usize = f[5].parse().unwrap_or(0);
        let qlen: usize = f[6].parse().unwrap_or(1);
        let coverage = (qend.saturating_sub(qstart) + 1) as f64 / qlen as f64;

        if coverage >= min_coverage {
            hits.push(Hit {
                query: f[0].to_string(),
                target: f[1].to_string(),
                tstart: f[2].parse().unwrap_or(1),
                tend: f[3].parse().unwrap_or(1),
            });
        }
    }

    hits
}

pub fn run(args: ExtractArgs) {
    match Command::new("mmseqs").arg("--version").output() {
        Ok(_) => {}
        Err(_) => {
            eprintln!("Error: mmseqs not found. Make sure it is installed and in your PATH.");
            std::process::exit(1);
        }
    }

    let target_files = collect_targets(&args.targets);
    if target_files.is_empty() {
        eprintln!("Error: no target FASTA files found.");
        std::process::exit(1);
    }

    // Unique temp dir per process so parallel runs don't collide
    let tmp_dir = std::env::temp_dir().join(format!("cladekit_extract_{}", std::process::id()));
    fs::create_dir_all(&tmp_dir).expect("Could not create temp directory");

    let pooled_path = tmp_dir.join("pooled_targets.fasta");
    let results_path = tmp_dir.join("results.tsv");
    let mmseqs_tmp = tmp_dir.join("mmseqs_tmp");

    eprintln!("Pooling {} target files...", target_files.len());
    let lookup = pool_targets(&target_files, &pooled_path);

    let log_path = Path::new(&args.output).join("mmseqs.log");
    let log_file = File::create(&log_path).expect("Could not create mmseqs.log");
    let log_file2 = log_file
        .try_clone()
        .expect("Could not clone log file handle");

    eprintln!("Running MMseqs2 easy-search...");
    let status = Command::new("mmseqs")
        .args([
            "easy-search",
            &args.reference,
            pooled_path.to_str().unwrap(),
            results_path.to_str().unwrap(),
            mmseqs_tmp.to_str().unwrap(),
            "--search-type",
            "3", // nucleotide-vs-nucleotide
            "-s",
            &args.sensitivity.to_string(),
            "--format-output",
            "query,target,tstart,tend,qstart,qend,qlen",
        ])
        .stdout(log_file)
        .stderr(log_file2)
        .status()
        .expect("Failed to run mmseqs");

    if !status.success() {
        eprintln!(
            "Error: mmseqs easy-search failed. See {}",
            log_path.display()
        );
        std::process::exit(1);
    }

    eprintln!("Parsing results...");
    let hits = parse_hits(&results_path, args.min_coverage);

    fs::create_dir_all(&args.output).expect("Could not create output directory");

    // One output file per gene, opened lazily as we encounter each gene name
    let mut gene_writers: HashMap<String, File> = HashMap::new();

    for hit in &hits {
        let (original_header, seq, filename) = match lookup.get(&hit.target) {
            Some(t) => t,
            None => {
                eprintln!("Warning: '{}' not found in lookup, skipping.", hit.target);
                continue;
            }
        };

        // MMseqs2 coordinates are 1-based inclusive. Convert to 0-based for Rust slicing.
        // tstart may be > tend on minus-strand hits — take min/max to always get a valid range.
        let raw_start = hit.tstart.min(hit.tend) - 1;
        let raw_end = hit.tstart.max(hit.tend);
        let start = raw_start.saturating_sub(args.flank);
        let end = (raw_end + args.flank).min(seq.len());
        let extracted = &seq[start..end];

        let writer = gene_writers.entry(hit.query.clone()).or_insert_with(|| {
            let out_path = Path::new(&args.output).join(format!("{}.fasta", hit.query));
            File::create(&out_path).expect("Could not create output file")
        });

        writeln!(writer, ">{} [ref={} hit={} {}-{}]", original_header, hit.query, filename, start, end).unwrap();
        writeln!(writer, "{}", extracted).unwrap();
    }

    eprintln!(
        "Done. Extracted {} gene(s) from {} hits.",
        gene_writers.len(),
        hits.len()
    );

    if !args.keep_intermediates {
        fs::remove_dir_all(&tmp_dir).expect("Could not remove temp directory");
    } else {
        eprintln!("Intermediates kept at: {}", tmp_dir.display());
    }
}
