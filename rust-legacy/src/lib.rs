use std::fs::File;
use std::io::{BufRead, BufReader};

/// Parse a FASTA file into a list of (header, sequence) pairs in file order.
/// A Vec (not a map) is used so output order is deterministic and matches the
/// input — every caller iterates these in order; the ones that need keyed
/// lookups (concat's matching, filter's loci counts) build their own maps.
/// When validate_equal is true, all sequences must be the same length (for alignments).
pub fn parse_fasta(
    filename: &str,
    validate_equal: bool,
) -> Result<(Vec<(String, String)>, usize), String> {
    let file = File::open(filename).map_err(|e| format!("Could not open {}: {}", filename, e))?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();
    let mut expected_length: Option<usize> = None;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading file: {}", e))?;
        let line = line.trim().to_string();

        if let Some(header) = line.strip_prefix('>') {
            // Save the previous sequence before starting a new one
            if !current_header.is_empty() {
                current_seq = current_seq.to_uppercase();
                sequences.push((current_header.clone(), current_seq.clone()));
                if validate_equal {
                    match expected_length {
                        None => expected_length = Some(current_seq.len()),
                        Some(len) => {
                            if current_seq.len() != len {
                                return Err(format!(
                                    "Error: Sequence length mismatch in {} : {}",
                                    filename, current_header
                                ));
                            }
                        }
                    }
                }
            }
            current_header = header.to_string();
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.push_str(&line);
        }
    }

    // Save the last sequence (loop only saves when it hits the next '>')
    if !current_header.is_empty() {
        current_seq = current_seq.to_uppercase();
        // None means this is the only sequence — nothing to compare against
        if validate_equal
            && let Some(len) = expected_length
            && current_seq.len() != len
        {
            return Err(format!(
                "Error: Sequence length mismatch in {} : {}",
                filename, current_header
            ));
        }
        sequences.push((current_header, current_seq));
    }

    let length = sequences.first().map_or(0, |(_, s)| s.len());
    Ok((sequences, length))
}

/// Detect whether sequences are DNA or protein.
/// Allows IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V) in addition to A/T/C/G/N/-.
pub fn is_dna(sequences: &[(String, String)]) -> bool {
    for (_, seq) in sequences {
        for ch in seq.chars() {
            match ch {
                'A' | 'T' | 'C' | 'G' | 'N' | '-' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B'
                | 'D' | 'H' | 'V' | '?' => {}
                _ => return false,
            }
        }
    }
    true
}

/// Reverse-complement a nucleotide sequence.
///
/// Covers the same alphabet `is_dna` accepts, so anything that passes that
/// check round-trips: the IUPAC ambiguity codes complement to their mirror
/// (R<->Y, K<->M, B<->V, D<->H; S/W/N are self-complementary), and gap and
/// missing characters pass through so an aligned sequence keeps its columns.
/// Case is not preserved — `parse_fasta` already uppercases everything, and
/// every caller works on its output.
///
/// Anything outside that alphabet is left as-is rather than silently turned
/// into an N: a caller that hands this protein has a bug, and mangling the
/// residues would hide it.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'U' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            other => other,
        })
        .collect()
}

/// Pretty print a table of rows. When pretty is false, prints tab-separated.
/// When pretty is true, pads columns to align.
pub fn print_table(rows: &[Vec<String>], pretty: bool) {
    if !pretty {
        for row in rows {
            println!("{}", row.join("\t"));
        }
        return;
    }
    let num_cols = rows[0].len();
    let mut widths = vec![0usize; num_cols];
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if cell.len() > widths[i] {
                widths[i] = cell.len();
            }
        }
    }
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if i > 0 {
                print!("  ");
            }
            print!("{:<width$}", cell, width = widths[i]);
        }
        println!();
    }
}

/// Read a taxa list file (one name per line) into a Vec.
pub fn load_taxa_list(filename: &str) -> Result<Vec<String>, String> {
    let file = File::open(filename).map_err(|e| format!("Could not open {}: {}", filename, e))?;
    let reader = BufReader::new(file);

    let mut taxa = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading file: {}", e))?;
        let line = line.trim().to_string();
        if !line.is_empty() {
            taxa.push(line);
        }
    }

    Ok(taxa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_round_trips() {
        let seq = "ATCGGGCTA";
        assert_eq!(reverse_complement(seq), "TAGCCCGAT");
        // Applying it twice must return the original, or extract would silently
        // flip records that were already on the plus strand.
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq);
    }

    #[test]
    fn reverse_complement_handles_iupac_and_gaps() {
        // R<->Y, K<->M, B<->V, D<->H; S/W/N self-complement; gaps pass through
        // so an aligned sequence keeps its column count.
        assert_eq!(reverse_complement("RYKMBVDH"), "DHBVKMRY");
        assert_eq!(reverse_complement("SWN"), "NWS");
        assert_eq!(reverse_complement("AT-CG"), "CG-AT");
        assert_eq!(reverse_complement("AT?CG"), "CG?AT");
        // Length is preserved for every accepted character, which is what makes
        // it safe to apply to an alignment row.
        let iupac = "ATCGNRYSWKMBDHV-?";
        assert_eq!(reverse_complement(iupac).len(), iupac.len());
    }

    #[test]
    fn reverse_complement_agrees_with_is_dna_alphabet() {
        // Anything is_dna accepts must survive a round trip unchanged, so the
        // two functions cannot drift apart as the alphabet grows.
        let seq = "ATCGNRYSWKMBDHV-?".to_string();
        assert!(is_dna(&[("h".to_string(), seq.clone())]));
        assert_eq!(reverse_complement(&reverse_complement(&seq)), seq);
    }
}
