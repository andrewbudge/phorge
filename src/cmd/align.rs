use std::fs::{File, OpenOptions};
use std::path::Path;
use std::process::Command;

use clap::Args;

#[derive(Args)]
pub struct AlignArgs {
    /// Alignment program name or path. Supported: mafft, muscle.
    /// Defaults: mafft uses --auto, muscle uses -align.
    #[arg(short, long)]
    pub program: String,

    /// Input unaligned FASTA files (glob or list)
    #[arg(short, long, required = true, num_args = 1..)]
    pub input: Vec<String>,

    /// Suffix to append to aligned output filenames (default: _aln)
    #[arg(short, long, default_value = "_aln")]
    pub extension: String,

    /// Output directory for aligned files (also writes align.log here)
    #[arg(short, long)]
    pub output: String,

    /// Extra flags passed verbatim to the aligner (after --).
    /// Replaces the default flag: --auto for mafft, -align for muscle.
    /// Example: phorge align -p mafft ... -- --thread 4 --maxiterate 1000
    /// Example: phorge align -p muscle ... -- -super5
    #[arg(last = true)]
    pub passthrough: Vec<String>,
}

pub fn run(args: AlignArgs) {
    // Check if program exists / is callable
    match Command::new(&args.program).arg("--version").output() {
        Ok(_) => {}
        Err(_) => {
            eprintln!(
                "Error: '{}' not found. Make sure it is installed and in your PATH.",
                args.program
            );
            std::process::exit(1);
        }
    }

    std::fs::create_dir_all(&args.output).expect("Could not create output directory");

    // Create a log file for aligner stderr output
    let log_path = Path::new(&args.output).join("align.log");
    File::create(&log_path).expect("Could not create log file");

    let program_name = Path::new(&args.program)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(&args.program)
        .to_ascii_lowercase();
    let is_muscle = program_name.contains("muscle");

    for input_path in &args.input {
        let path = Path::new(input_path);
        let stem = path
            .file_stem()
            .expect("Could not get filename")
            .to_str()
            .unwrap();
        let ext = path
            .extension()
            .map(|e| e.to_str().unwrap())
            .unwrap_or("fasta");
        let output_filename = format!("{}{}.{}", stem, args.extension, ext);
        let output_path = Path::new(&args.output).join(&output_filename);

        eprint!("Aligning {}...", stem);

        let log_file = OpenOptions::new()
            .append(true)
            .open(&log_path)
            .expect("Could not open log file");

        let mut cmd = Command::new(&args.program);

        let status = if is_muscle {
            // MUSCLE v5: muscle [-super5] -align input.fa -output out.fa
            if args.passthrough.is_empty() {
                cmd.arg("-align");
            } else {
                cmd.args(&args.passthrough);
            }
            cmd.arg(input_path)
                .arg("-output")
                .arg(&output_path)
                .stderr(log_file)
                .status()
                .expect("Failed to run aligner")
        } else {
            // MAFFT (and others): mafft [--auto] input.fa > output.fa
            let output_file = File::create(&output_path).expect("Could not create output file");
            if args.passthrough.is_empty() {
                cmd.arg("--auto");
            } else {
                cmd.args(&args.passthrough);
            }
            cmd.arg(input_path)
                .stdout(output_file)
                .stderr(log_file)
                .status()
                .expect("Failed to run aligner")
        };

        if !status.success() {
            eprintln!("FAILED (see {})", log_path.display());
            std::process::exit(1);
        }
        eprintln!("done");
    }

    eprintln!("Done. Aligned {} files.", args.input.len());
}
