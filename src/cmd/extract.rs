use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

use clap::{ArgGroup, Args};

use phorge::{parse_fasta, reverse_complement};

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

    /// Directory of per-gene reference FASTAs, or a single such file (gene name
    /// = filename stem, e.g. COI.fasta -> COI). Each file may hold many
    /// sequences to cover divergence. Pipeline form. Repeat the flag to draw
    /// from several directories.
    #[arg(long, num_args = 1)]
    pub refs: Option<Vec<String>>,

    /// Target FASTA files, or a directory containing them (e.g. fetch's output dir)
    #[arg(required = true, num_args = 1..)]
    pub targets: Vec<String>,

    /// Output directory for per-gene FASTAs
    #[arg(short, long)]
    pub output: String,

    /// Minimum MMseqs2 sequence identity for a hit to be kept (0.0–1.0)
    #[arg(short, long, default_value_t = 0.7)]
    pub min_identity: f64,

    /// Minimum fraction of the REFERENCE gene the hit must span (0.0–1.0).
    /// Identity says the bases that aligned matched; coverage says enough of
    /// the gene aligned at all. Set 0.0 to keep partial hits of any length.
    #[arg(short, long, default_value_t = 0.5)]
    pub coverage: f64,

    /// Extra bases to grab on either side of each hit
    #[arg(long, default_value_t = 0)]
    pub flank: usize,

    /// MMseqs2 sensitivity (1.0=fast, 7.5=max). Higher finds more divergent hits.
    #[arg(short, long, default_value_t = 5.7)]
    pub sensitivity: f64,

    /// Max targets each reference is aligned against (MMseqs2 --max-seqs). The
    /// prefilter keeps only the top N targets per gene by k-mer score, so if you
    /// have MORE target sequences than this, the least-similar ones silently get
    /// no hit. Raise above your target count for large runs.
    #[arg(long, default_value_t = 300)]
    pub max_seqs: usize,

    /// Cap MMseqs2 RAM (e.g. 8G); splits the search to stay under it. Default: unlimited.
    #[arg(long)]
    pub max_memory_limit: Option<String>,

    /// MMseqs2 k-mer length (MMseqs2 -k). Left to MMseqs2 to choose when unset.
    /// Its automatic choice for nucleotide search is 15, whose index table needs
    /// more RAM than a 16 GB machine has; 13 costs a few hundred MB and finds
    /// the same hits. Lower k is slower and less specific, so only reach for
    /// this when the search dies on memory.
    #[arg(long)]
    pub kmer_length: Option<usize>,

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

// `--refs` takes one value per occurrence, so a shell glob (`--refs refs/*.fa`)
// silently splits: the first file becomes the reference and the rest land in
// TARGETS. That costs you every gene but the first, with no error — so when a
// target sits in the same directory as a `--refs` file, say so. A warning, not
// an error: keeping references and targets in one directory is unusual but not
// forbidden.
fn warn_refs_leaked_into_targets(args: &ExtractArgs, target_files: &[String]) {
    let Some(refs) = &args.refs else { return };
    let ref_dirs: Vec<&Path> = refs
        .iter()
        .map(Path::new)
        .filter(|p| p.is_file())
        .filter_map(|p| p.parent())
        .collect();

    let leaked: Vec<&String> = target_files
        .iter()
        .filter(|t| Path::new(t).parent().is_some_and(|d| ref_dirs.contains(&d)))
        .collect();

    if !leaked.is_empty() {
        eprintln!(
            "Warning: {} target file(s) sit in a --refs directory: {}",
            leaked.len(),
            leaked
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        eprintln!(
            "  --refs takes ONE value; a glob like `--refs refs/*.fasta` leaves the rest as targets."
        );
        eprintln!("  Pass the directory instead: --refs {}", {
            let d = ref_dirs[0].to_string_lossy();
            if d.is_empty() {
                ".".to_string()
            } else {
                d.into_owned()
            }
        });
    }
}

// Writes a pooled reference FASTA whose record IDs encode the gene name as
// `gene::N`, so MMseqs2 reports the gene in its `query` column. The gene name
// comes from the record header (single --reference file) or the filename stem
// (--refs, whose directories are expanded to their FASTAs first). Returns the
// number of reference records written.
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
        // Each --refs value is a directory of per-gene FASTAs or a single such
        // file; collect_targets flattens both to a plain file list, so a gene's
        // name is always the stem of a real file by the time we get here.
        for file in collect_targets(refs) {
            let path = Path::new(&file);
            let gene = path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .replace(' ', "_");
            let (seqs, _) = parse_fasta(&file, false).expect("Failed to read reference FASTA");
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
    coverage: f64,
    /// Always the lower coordinate, whatever strand the hit was on.
    tstart: usize,
    /// Always the higher coordinate.
    tend: usize,
    /// True when MMseqs2 reported the hit on the minus strand (it emits
    /// qstart > qend for those). The extracted slice must be
    /// reverse-complemented before it is written, or the record comes out in
    /// the opposite orientation to every plus-strand record for the same gene
    /// and no aligner can reconcile them.
    minus_strand: bool,
}

// Parses MMseqs2 tabular output. Both quality gates are applied in-engine
// (--min-seq-id and -c), so every row here already passed; we carry fident and
// qcov through only so they can be recorded in the output header.
// Expected --format-output: query,target,fident,qcov,qstart,qend,tstart,tend
fn parse_hits(tsv_path: &Path) -> Vec<Hit> {
    let file = File::open(tsv_path).expect("Could not open MMseqs2 output");
    let reader = BufReader::new(file);
    let mut hits = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Error reading MMseqs2 output");
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }

        // query is `gene::N`; recover the gene name before the separator.
        let gene = f[0].split("::").next().unwrap_or(f[0]).to_string();

        // MMseqs2 signals a minus-strand hit on the QUERY coordinates, not the
        // target: it reports qstart > qend and leaves the target range forward
        // (e.g. qstart=660 qend=1, tstart=201 tend=860). Reading the strand off
        // tstart/tend instead would never fire.
        let qstart: usize = f[4].parse().unwrap_or(1);
        let qend: usize = f[5].parse().unwrap_or(1);
        let tstart: usize = f[6].parse().unwrap_or(1);
        let tend: usize = f[7].parse().unwrap_or(1);

        hits.push(Hit {
            gene,
            target: f[1].to_string(),
            identity: f[2].parse().unwrap_or(0.0),
            coverage: f[3].parse().unwrap_or(0.0),
            // Normalised defensively; the target range is already forward.
            tstart: tstart.min(tend),
            tend: tstart.max(tend),
            minus_strand: qstart > qend,
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
    warn_refs_leaked_into_targets(&args, &target_files);

    // Output dir holds the per-gene FASTAs and the mmseqs log; create it up front.
    fs::create_dir_all(&args.output).expect("Could not create output directory");

    // Unique temp dir per process so parallel runs don't collide
    let tmp_dir = std::env::temp_dir().join(format!("phorge_extract_{}", std::process::id()));
    fs::create_dir_all(&tmp_dir).expect("Could not create temp directory");

    let pooled_ref_path = tmp_dir.join("pooled_reference.fasta");
    let pooled_path = tmp_dir.join("pooled_targets.fasta");
    let results_path = tmp_dir.join("results.tsv");
    let mmseqs_tmp = tmp_dir.join("mmseqs_tmp");

    let n_refs = pool_references(&args, &pooled_ref_path);
    if n_refs == 0 {
        // Nothing to search with: an empty --refs directory, or a reference file
        // holding no records. Fail here rather than handing MMseqs2 an empty
        // query set and reporting zero hits as though the search had run.
        eprintln!("Error: no reference sequences found.");
        std::process::exit(1);
    }
    eprintln!("Pooled {} reference sequence(s).", n_refs);
    eprintln!("Pooling {} target files...", target_files.len());
    let lookup = pool_targets(&target_files, &pooled_path);

    let log_path = Path::new(&args.output).join("mmseqs.log");
    let log_file = File::create(&log_path).expect("Could not create mmseqs.log");
    let log_file2 = log_file
        .try_clone()
        .expect("Could not clone log file handle");

    eprintln!(
        "Running MMseqs2 easy-search (min identity {}, min coverage {})...",
        args.min_identity, args.coverage
    );
    let mut cmd = Command::new("mmseqs");
    cmd.args([
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
        "-c",
        &args.coverage.to_string(),
        // cov-mode 2 = fraction of the QUERY covered. The reference is the
        // query here, so this asks "how much of the gene did I recover?".
        // Mode 0 (query AND target) would be useless: a 650 bp COI hit inside
        // a 16 kb mitogenome covers ~4% of the target and every hit would die.
        "--cov-mode",
        "2",
        "--max-seqs",
        &args.max_seqs.to_string(),
        "--format-output",
        "query,target,fident,qcov,qstart,qend,tstart,tend",
    ]);

    // Only cap memory when the user asks; otherwise let MMseqs2 use what it wants.
    if let Some(limit) = &args.max_memory_limit {
        cmd.args(["--split-memory-limit", limit]);
    }

    // Same rule for k: no flag means MMseqs2 picks, which is right on a machine
    // with the RAM for it.
    if let Some(k) = args.kmer_length {
        cmd.args(["-k", &k.to_string()]);
    }

    let status = cmd
        .stdout(log_file)
        .stderr(log_file2)
        .status()
        .expect("Failed to run mmseqs");

    if !status.success() {
        eprintln!(
            "Error: mmseqs easy-search failed. See {}",
            log_path.display()
        );
        // MMseqs2 reports an out-of-memory prefilter as a bare "Cannot fit
        // databases into NG" buried in a very long log. The cause is almost
        // always the k-mer index table, whose size is set by k rather than by
        // how big the database is, so --max-memory-limit cannot rescue it.
        if fs::read_to_string(&log_path)
            .is_ok_and(|log| log.contains("Cannot fit databases into"))
        {
            eprintln!(
                "\nThe search ran out of memory building its k-mer index. \
                 Splitting cannot help: the index table is sized by k, not by \
                 the database. Retry with a shorter k, e.g. --kmer-length 13."
            );
        }
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

        // MMseqs2 coordinates are 1-based inclusive; parse_hits already put them
        // in (low, high) order, so this only has to shift to 0-based for slicing.
        let start = hit.tstart.saturating_sub(1).saturating_sub(args.flank);
        let end = (hit.tend + args.flank).min(seq.len());

        // A minus-strand hit sits in the target in the opposite orientation to
        // the reference, so the raw slice reads backwards. Flip it, and the
        // flank comes along correctly: reverse-complementing the widened slice
        // puts the upstream flank where a reader expects it.
        let extracted = if hit.minus_strand {
            reverse_complement(&seq[start..end])
        } else {
            seq[start..end].to_string()
        };
        let strand = if hit.minus_strand { '-' } else { '+' };

        let writer = gene_writers.entry(hit.gene.clone()).or_insert_with(|| {
            let out_path = Path::new(&args.output).join(format!("{}.fasta", hit.gene));
            File::create(&out_path).expect("Could not create output file")
        });

        writeln!(
            writer,
            ">{} [gene={} ident={:.3} cov={:.3} strand={} src={} {}-{}]",
            original_header, hit.gene, hit.identity, hit.coverage, strand, filename, start, end
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
