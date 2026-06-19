use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

use clap::{ArgGroup, Args};

use cladekit::parse_fasta;

#[derive(Args)]
// Exactly one reference form is required: a single multi-gene file, or one
// file per gene. They are mutually exclusive.
#[command(group(
    ArgGroup::new("refsource")
        .required(true)
        .args(["reference", "refs"])
))]
pub struct ExtractArgs {
    /// Single reference FASTA with labeled gene records (gene name = each
    /// record header, e.g. >COI, >ND2). For ad-hoc / standalone use.
    #[arg(short, long)]
    pub reference: Option<String>,

    /// Per-gene reference FASTAs (gene name = filename stem, e.g. COI.fasta ->
    /// COI). Each file may hold many sequences to cover divergence. Pipeline form.
    #[arg(long, num_args = 1..)]
    pub refs: Option<Vec<String>>,

    /// Target FASTA files or a directory containing them (e.g. fetch's raw/ dir)
    #[arg(short, long, required = true, num_args = 1..)]
    pub targets: Vec<String>,

    /// Output directory for per-gene FASTAs
    #[arg(short, long)]
    pub output: String,

    /// Minimum MMseqs2 sequence identity for a hit to be kept (0.0–1.0).
    /// This is the sole quality gate; pick references that cover your taxa.
    #[arg(long, default_value_t = 0.7)]
    pub min_identity: f64,

    /// Extra bases to grab on either side of each hit
    #[arg(long, default_value_t = 0)]
    pub flank: usize,

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
            for entry in path.read_dir().unwrap_or_else(|e| {
                eprintln!("Error: could not read directory '{}': {}", target, e);
                std::process::exit(1);
            }) {
                let entry = entry.expect("Failed to read directory entry");
                let p = entry.path();
                if p.extension()
                    .is_some_and(|e| e == "fasta" || e == "fa" || e == "fna" || e == "fas")
                {
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

// Writes a pooled reference FASTA whose record IDs encode the gene name as
// `gene::N`, so MMseqs2 reports the gene in its `query` column. The gene name
// comes from the record header (single --reference file) or the filename stem
// (per-gene --refs files). Returns the number of reference records written.
fn pool_references(args: &ExtractArgs, pooled_path: &Path) -> usize {
    let mut writer = File::create(pooled_path).expect("Could not create pooled reference file");
    let mut counter = 0usize;

    if let Some(reference) = &args.reference {
        let (seqs, _) = parse_fasta(reference, false).expect("Failed to read reference FASTA");
        for (header, seq) in &seqs {
            // gene = first whitespace-delimited token of the header
            let gene = header.split_whitespace().next().unwrap_or(header);
            writeln!(writer, ">{}::{}", gene, counter).unwrap();
            writeln!(writer, "{}", seq).unwrap();
            counter += 1;
        }
    } else if let Some(refs) = &args.refs {
        for file in refs {
            let path = Path::new(file);
            let gene = path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .replace(' ', "_");
            let (seqs, _) = parse_fasta(file, false).expect("Failed to read reference FASTA");
            for (_header, seq) in &seqs {
                writeln!(writer, ">{}::{}", gene, counter).unwrap();
                writeln!(writer, "{}", seq).unwrap();
                counter += 1;
            }
        }
    }

    counter
}

// Writes a pooled FASTA to disk with "organism::seq_id" headers.
// Returns a lookup map from that key to (original_header, full_sequence, original_filename)
fn pool_targets(
    target_files: &[String],
    pooled_path: &Path,
) -> HashMap<String, (String, String, String)> {
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
    gene: String,
    target: String,
    identity: f64,
    tstart: usize,
    tend: usize,
}

// Parses MMseqs2 tabular output. The identity gate is applied in-engine via
// --min-seq-id, so every row here already passed; we just carry fident through
// for traceability in the output header.
// Expected --format-output: query,target,fident,tstart,tend
fn parse_hits(tsv_path: &Path) -> Vec<Hit> {
    let file = File::open(tsv_path).expect("Could not open MMseqs2 output");
    let reader = BufReader::new(file);
    let mut hits = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Error reading MMseqs2 output");
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }

        // query is `gene::N`; recover the gene name before the separator.
        let gene = f[0].split("::").next().unwrap_or(f[0]).to_string();

        hits.push(Hit {
            gene,
            target: f[1].to_string(),
            identity: f[2].parse().unwrap_or(0.0),
            tstart: f[3].parse().unwrap_or(1),
            tend: f[4].parse().unwrap_or(1),
        });
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

    // Output dir holds the per-gene FASTAs and the mmseqs log; create it up front.
    fs::create_dir_all(&args.output).expect("Could not create output directory");

    // Unique temp dir per process so parallel runs don't collide
    let tmp_dir = std::env::temp_dir().join(format!("cladekit_extract_{}", std::process::id()));
    fs::create_dir_all(&tmp_dir).expect("Could not create temp directory");

    let pooled_ref_path = tmp_dir.join("pooled_reference.fasta");
    let pooled_path = tmp_dir.join("pooled_targets.fasta");
    let results_path = tmp_dir.join("results.tsv");
    let mmseqs_tmp = tmp_dir.join("mmseqs_tmp");

    let n_refs = pool_references(&args, &pooled_ref_path);
    eprintln!("Pooled {} reference sequence(s).", n_refs);
    eprintln!("Pooling {} target files...", target_files.len());
    let lookup = pool_targets(&target_files, &pooled_path);

    let log_path = Path::new(&args.output).join("mmseqs.log");
    let log_file = File::create(&log_path).expect("Could not create mmseqs.log");
    let log_file2 = log_file
        .try_clone()
        .expect("Could not clone log file handle");

    eprintln!(
        "Running MMseqs2 easy-search (min identity {})...",
        args.min_identity
    );
    let status = Command::new("mmseqs")
        .args([
            "easy-search",
            pooled_ref_path.to_str().unwrap(),
            pooled_path.to_str().unwrap(),
            results_path.to_str().unwrap(),
            mmseqs_tmp.to_str().unwrap(),
            "--search-type",
            "3", // nucleotide-vs-nucleotide
            "-s",
            &args.sensitivity.to_string(),
            "--min-seq-id",
            &args.min_identity.to_string(),
            "--format-output",
            "query,target,fident,tstart,tend",
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
    let hits = parse_hits(&results_path);

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

        let writer = gene_writers.entry(hit.gene.clone()).or_insert_with(|| {
            let out_path = Path::new(&args.output).join(format!("{}.fasta", hit.gene));
            File::create(&out_path).expect("Could not create output file")
        });

        writeln!(
            writer,
            ">{} [gene={} ident={:.3} src={} {}-{}]",
            original_header, hit.gene, hit.identity, filename, start, end
        )
        .unwrap();
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
