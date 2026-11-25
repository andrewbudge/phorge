//! End-to-end check that `extract` returns every hit in the reference's
//! orientation.
//!
//! MMseqs2 searches both strands, so a gene encoded on the minus strand of a
//! target reports its coordinates backwards. Slicing those coordinates without
//! reverse-complementing yields a record that is a perfect hit by every metric
//! in its own header and still unusable: it aligns against its plus-strand
//! counterparts at roughly random-sequence identity, so the aligner produces
//! garbage columns and the taxon is later dropped as an outlier. Nothing about
//! that failure is loud, which is why it is worth an integration test.
//!
//! Requires MMseqs2 on PATH; skipped with a printed note when it is absent.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Mirrors extract.rs's own detection: `mmseqs --version` is not a valid
/// subcommand and exits non-zero, so only the spawn is checked. Testing
/// `status.success()` here would make this test skip itself on a machine where
/// MMseqs2 is installed and working.
fn mmseqs_available() -> bool {
    Command::new("mmseqs").arg("--version").output().is_ok()
}

fn phorge_bin() -> PathBuf {
    // Cargo tells the test where the binary it just built actually is. Deriving
    // it from current_exe() by hand breaks whenever Cargo changes its artifact
    // layout, and a broken path reads as "phorge missing" rather than "test
    // harness out of date" — which is how this test managed to never run.
    PathBuf::from(env!("CARGO_BIN_EXE_phorge"))
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// Deterministic pseudo-random ACGT, so a failure is reproducible.
fn dna(len: usize, seed: u64) -> String {
    let mut state = seed;
    (0..len)
        .map(|_| {
            // xorshift64: enough randomness for a fixture, no dependency.
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ['A', 'C', 'G', 'T'][(state % 4) as usize]
        })
        .collect()
}

/// Reads a FASTA into (first header token, sequence) pairs.
fn read_fasta(path: &Path) -> Vec<(String, String)> {
    let text = fs::read_to_string(path).expect("reading extract output");
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            let id = h.split_whitespace().next().unwrap_or(h).to_string();
            out.push((id, String::new()));
        } else if !line.is_empty() {
            out.last_mut().expect("sequence before any header").1 += line;
        }
    }
    out
}

#[test]
fn extract_returns_minus_strand_hits_in_reference_orientation() {
    if !mmseqs_available() {
        eprintln!("skipping: mmseqs not found on PATH");
        return;
    }

    let dir = std::env::temp_dir().join(format!("phorge_strand_test_{}", std::process::id()));
    let refs = dir.join("refs");
    fs::create_dir_all(&refs).expect("creating fixture directories");

    // One gene, placed in two organisms: forward in the first, reverse-
    // complemented in the second. Both are otherwise identical, so any
    // difference in the output is orientation and nothing else.
    let gene = dna(660, 0x5EED);
    let flank_a = dna(200, 0xA1);
    let flank_b = dna(200, 0xB2);

    fs::write(refs.join("COI.fasta"), format!(">coi_ref\n{gene}\n")).unwrap();
    fs::write(
        dir.join("plus.fna"),
        format!(">AB001.1 Forward organism\n{flank_a}{gene}{flank_b}\n"),
    )
    .unwrap();
    fs::write(
        dir.join("minus.fna"),
        format!(
            ">AB002.1 Reverse organism\n{flank_a}{}{flank_b}\n",
            reverse_complement(&gene)
        ),
    )
    .unwrap();

    let out = dir.join("genes");
    let status = Command::new(phorge_bin())
        .arg("extract")
        .args(["--refs", refs.to_str().unwrap()])
        // MMseqs2's k-mer index table is sized by k, not by the database, so
        // the default k=15 would reserve several GB to search two 1 KB
        // sequences and can OOM-kill the machine. Pin it small: this test
        // should cost milliseconds and a few hundred MB, not a crash.
        .args(["--kmer-length", "13"])
        .args(["-o", out.to_str().unwrap()])
        .arg(dir.join("plus.fna"))
        .arg(dir.join("minus.fna"))
        .output()
        .expect("running phorge extract");
    assert!(
        status.status.success(),
        "extract failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );

    let records = read_fasta(&out.join("COI.fasta"));
    assert_eq!(records.len(), 2, "expected one hit per organism");

    let plus = &records
        .iter()
        .find(|(id, _)| id == "AB001.1")
        .expect("plus-strand record missing")
        .1;
    let minus = &records
        .iter()
        .find(|(id, _)| id == "AB002.1")
        .expect("minus-strand record missing")
        .1;

    // The real assertion: both records come back as the reference gene. Before
    // the fix the minus record was its reverse complement instead.
    assert_eq!(plus, &gene, "plus-strand hit should match the reference");
    assert_eq!(
        minus, &gene,
        "minus-strand hit should be reverse-complemented back into the \
         reference's orientation"
    );

    // Stated as the biologist would see it: two copies of one gene must not
    // look like unrelated sequence to the aligner.
    let identity = plus
        .chars()
        .zip(minus.chars())
        .filter(|(a, b)| a == b)
        .count() as f64
        / plus.len() as f64;
    assert!(
        identity > 0.99,
        "the two records are only {:.1}% identical; a minus-strand hit is \
         being written in the wrong orientation",
        identity * 100.0
    );

    fs::remove_dir_all(&dir).ok();
}
