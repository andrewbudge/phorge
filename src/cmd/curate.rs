use phorge::parse_fasta;
use clap::Args;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Args)]
pub struct CurateArgs {
    /// Aligned FASTA files to trim
    #[arg(required = true, num_args = 1..)]
    pub input: Vec<String>,

    /// Which column properties to keep (combine letters):
    ///   p  parsimony-informative (at least 2 chars each appearing ≥2 times)
    ///   c  constant (only one character, appearing ≥2 times)
    ///   s  smart-gap: auto-determine gap threshold via slope analysis (default)
    ///   g  gappy: trim sites at or above --gap-threshold (fixed)
    ///
    /// Common combinations:
    ///   ps   keep parsimony-informative, smart-gap filter        [default]
    ///   pcs  keep parsimony-informative + constant, smart-gap filter
    ///   pg   keep parsimony-informative, fixed gap filter
    ///   pcg  keep parsimony-informative + constant, fixed gap filter
    ///   p    keep parsimony-informative only (no gap filter)
    ///   pc   keep parsimony-informative + constant (no gap filter)
    ///   s    smart-gap filter only (no parsimony filter)
    ///   g    fixed gap filter only
    #[arg(short, long, default_value = "ps", verbatim_doc_comment)]
    pub keep: String,

    /// Max gappiness per site (0.0–1.0); used when 'g' is in --keep
    #[arg(long, default_value_t = 0.9)]
    pub gap_threshold: f64,

    /// Suffix to append to output filenames
    #[arg(short, long, default_value = "_curated")]
    pub extension: String,

    /// Output directory (if omitted, writes to stdout)
    #[arg(short, long)]
    pub output: Option<String>,
}

enum SeqType {
    Nt,
    Aa,
}

fn detect_seq_type(seqs: &[(String, Vec<u8>)]) -> SeqType {
    let mut unique = std::collections::HashSet::new();
    let mut count = 0;
    for &b in &seqs[0].1 {
        if b == b'-' || b == b'*' {
            continue;
        }
        unique.insert(b.to_ascii_uppercase());
        count += 1;
        if count >= 100 {
            break;
        }
    }
    if unique.len() > 5 {
        SeqType::Aa
    } else {
        SeqType::Nt
    }
}

fn gap_chars_for(seq_type: &SeqType) -> &'static [u8] {
    match seq_type {
        SeqType::Nt => b"-?*XxNn",
        SeqType::Aa => b"-?*Xx",
    }
}

// Finds the gap fraction threshold where the slope of the gap distribution curve
// changes most sharply — the ClipKIT smart-gap heuristic. Returns 1.0 as a safe
// fallback (no trimming) when the distribution is too uniform to find a split.
fn smart_gap_threshold(site_gappiness: &[f64]) -> f64 {
    if site_gappiness.is_empty() {
        return 1.0;
    }

    // Count frequency of each unique gappiness value (rounded to 4 decimal places).
    // Use integer keys (multiply by 10000) to avoid f64 hashing issues.
    let mut counts: std::collections::BTreeMap<u32, usize> = std::collections::BTreeMap::new();
    for &g in site_gappiness {
        let key = (g * 10_000.0).round() as u32;
        *counts.entry(key).or_insert(0) += 1;
    }

    // Sort descending by gap value.
    let mut gaps_arr: Vec<(f64, f64)> = counts
        .iter()
        .map(|(&k, &v)| (k as f64 / 10_000.0, v as f64))
        .collect();
    gaps_arr.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    if gaps_arr.len() < 2 {
        return gaps_arr.first().map(|&(g, _)| g).unwrap_or(1.0);
    }

    let total = site_gappiness.len() as f64;

    // Build cumulative sum of normalized counts.
    let mut cumsum: Vec<f64> = Vec::with_capacity(gaps_arr.len());
    let mut running = 0.0;
    for &(_, c) in &gaps_arr {
        running += c / total;
        cumsum.push(running);
    }

    // Slope between each consecutive pair of gap values; only use first half.
    let mut slopes: Vec<f64> = Vec::new();
    for i in 1..gaps_arr.len() {
        let delta_gap = gaps_arr[i].0 - gaps_arr[i - 1].0;
        if delta_gap.abs() > f64::EPSILON {
            slopes.push(((cumsum[i] - cumsum[i - 1]) / delta_gap).abs());
        }
    }
    let slopes = &slopes[..slopes.len() / 2];

    if slopes.len() < 2 {
        return if slopes.is_empty() {
            1.0
        } else {
            gaps_arr[0].0
        };
    }

    // Return the gap value at the index where consecutive slope difference is largest.
    let max_diff_idx = slopes
        .windows(2)
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            (a[1] - a[0])
                .abs()
                .partial_cmp(&(b[1] - b[0]).abs())
                .unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    gaps_arr[max_diff_idx].0
}

struct SiteStats {
    gappiness: f64,
    is_parsimony_informative: bool,
    is_constant: bool,
}

pub fn run(args: CurateArgs) {
    if let Some(ref out_dir) = args.output {
        std::fs::create_dir_all(out_dir).expect("Could not create output directory");
    }

    let total_files = args.input.len();

    for input_path in &args.input {
        let (sequences, _) = parse_fasta(input_path, true).expect(
            "Could not read alignment (file not found or sequences are not the same length)",
        );

        let seqs: Vec<(String, Vec<u8>)> = sequences
            .into_iter()
            .map(|(h, s)| (h, s.into_bytes()))
            .collect();

        let n_taxa = seqs.len();
        let aln_len = seqs[0].1.len();
        let gap_chars = gap_chars_for(&detect_seq_type(&seqs));

        let keep_p = args.keep.contains('p');
        let keep_c = args.keep.contains('c');
        let keep_s = args.keep.contains('s');
        let keep_g = args.keep.contains('g');

        // First pass: compute per-site stats.
        let site_stats: Vec<SiteStats> = (0..aln_len)
            .map(|col| {
                let mut gap_count = 0usize;
                let mut freq = HashMap::<u8, usize>::new();

                for (_, seq) in &seqs {
                    let base = seq[col];
                    if gap_chars.contains(&base) {
                        gap_count += 1;
                    } else {
                        *freq.entry(base).or_insert(0) += 1;
                    }
                }

                // Round to 4 decimals to match ClipKIT, which compares the
                // rounded site gappiness (np.around(..., 4)) against the
                // threshold. (np.around is round-half-to-even; f64::round is
                // round-half-away-from-zero — they differ only for exact .5
                // ties at the 5th decimal, which k/n gappiness rarely hits.)
                let gappiness = ((gap_count as f64 / n_taxa as f64) * 10_000.0).round() / 10_000.0;
                let qualifying = freq.values().filter(|&&c| c >= 2).count();
                SiteStats {
                    gappiness,
                    is_parsimony_informative: qualifying >= 2,
                    is_constant: freq.len() == 1 && qualifying == 1,
                }
            })
            .collect();

        // Smart-gap determines its threshold from the full gappiness distribution.
        let smart_threshold = if keep_s {
            let vals: Vec<f64> = site_stats.iter().map(|s| s.gappiness).collect();
            smart_gap_threshold(&vals)
        } else {
            args.gap_threshold
        };

        let kept_cols: Vec<usize> = site_stats
            .iter()
            .enumerate()
            .filter_map(|(col, s)| {
                // ClipKIT uses a mode-dependent boundary: kpi-* modes trim
                // strictly ('> threshold', so keep '<='), while smart-gap /
                // gappy / kpic-* trim on '>=' (keep '<'). Mirror that exactly.
                let kpi_only = keep_p && !keep_c;
                let gap_ok = if keep_s {
                    if kpi_only { s.gappiness <= smart_threshold } else { s.gappiness < smart_threshold }
                } else if keep_g {
                    if kpi_only { s.gappiness <= args.gap_threshold } else { s.gappiness < args.gap_threshold }
                } else {
                    true
                };

                let class_ok = if keep_p && keep_c {
                    s.is_parsimony_informative || s.is_constant
                } else if keep_p {
                    s.is_parsimony_informative
                } else {
                    true
                };

                if gap_ok && class_ok { Some(col) } else { None }
            })
            .collect();

        let kept = kept_cols.len();
        let path = Path::new(input_path);
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("fasta");

        if let Some(ref out_dir) = args.output {
            let out_filename = format!("{}{}.{}", stem, args.extension, ext);
            let out_path = Path::new(out_dir).join(&out_filename);
            let file = File::create(&out_path).expect("Could not create output file");
            let mut writer = BufWriter::new(file);

            for (header, seq) in &seqs {
                let trimmed: String = kept_cols.iter().map(|&col| seq[col] as char).collect();
                writeln!(writer, ">{}", header).unwrap();
                writeln!(writer, "{}", trimmed).unwrap();
            }

            eprintln!(
                "{}: {} → {} sites ({} removed)",
                stem,
                aln_len,
                kept,
                aln_len - kept
            );
        } else {
            for (header, seq) in &seqs {
                let trimmed: String = kept_cols.iter().map(|&col| seq[col] as char).collect();
                println!(">{}", header);
                println!("{}", trimmed);
            }

            eprintln!("Sites in:      {}", aln_len);
            eprintln!("Sites kept:    {}", kept);
            eprintln!("Sites removed: {}", aln_len - kept);
        }
    }

    if args.output.is_some() && total_files > 1 {
        eprintln!("Done. Curated {} files.", total_files);
    }
}
