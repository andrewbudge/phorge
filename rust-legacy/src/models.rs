use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// A single nucleotide record returned by `query`, populated from an NCBI
/// esummary docsum. No homology/locus information is present at this stage —
/// that is determined later by MMseqs2 in `extract`.
#[derive(Serialize, Deserialize, Debug)]
pub struct Accession {
    pub accession: String,
    pub taxon_name: String,
    /// The record's own source-organism TaxID (from the docsum), e.g. a single
    /// species. Distinct from [`Accession::query_taxid`].
    pub taxid: u64,
    pub length: usize,
    /// GenBank title string (e.g. "...cytochrome c oxidase subunit I..."). Recorded
    /// for traceability only; deliberately NOT used as a homology/quality filter.
    pub gene_annotation: String,
    pub refseq: bool,
    pub source_db: String,
    /// The higher-level TaxID that was queried to surface this record (the
    /// ingroup/outgroup root the user supplied). Stamped here so provenance
    /// survives flattening, and so the lineage check has the root to compare
    /// [`Accession::taxid`] against.
    pub query_taxid: u64,
    /// Which group the querying root belongs to. Stamped per-record (rather than
    /// only on the enclosing query result) so it travels with the accession once
    /// the table is flattened downstream.
    pub taxon_group: TaxonGroup,
    /// Flagged (never silently dropped) when a record's source TaxID does not fall
    /// within the queried taxonomic root. Always `false` until the lineage check is
    /// implemented.
    pub taxonomic_outlier: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum TaxonGroup {
    Ingroup,
    Outgroup,
}

impl TaxonGroup {
    pub fn as_str(self) -> &'static str {
        match self {
            TaxonGroup::Ingroup => "ingroup",
            TaxonGroup::Outgroup => "outgroup",
        }
    }

    fn parse(s: &str) -> Result<Self> {
        match s {
            "ingroup" => Ok(TaxonGroup::Ingroup),
            "outgroup" => Ok(TaxonGroup::Outgroup),
            other => bail!("unknown taxon_group {other:?} (expected ingroup/outgroup)"),
        }
    }
}

/// Column order for the query TSV. The writer emits these left-to-right; the
/// reader maps by name so a user can filter/reorder columns with `awk`/`cut`
/// before handing the file to `fetch`. Keep in sync with the field reads below.
const TSV_COLUMNS: [&str; 10] = [
    "accession",
    "taxon_name",
    "taxid",
    "length",
    "gene_annotation",
    "refseq",
    "source_db",
    "query_taxid",
    "taxon_group",
    "taxonomic_outlier",
];

/// Render accessions as a tab-separated table with a header row. One row per
/// accession — the nested per-taxon grouping from the old JSON is flattened,
/// since `query_taxid`/`taxon_group` on each row already carry that provenance.
/// Pure (no I/O) so the caller can send it to a file or stdout.
pub fn render_query_tsv(records: &[Accession]) -> String {
    let mut out = String::new();
    out.push_str(&TSV_COLUMNS.join("\t"));
    out.push('\n');
    for r in records {
        // A tab/newline in a free-text field would corrupt the row layout. GenBank
        // titles never contain them in practice, but sanitize defensively.
        let fields = [
            clean(&r.accession),
            clean(&r.taxon_name),
            r.taxid.to_string(),
            r.length.to_string(),
            clean(&r.gene_annotation),
            r.refseq.to_string(),
            clean(&r.source_db),
            r.query_taxid.to_string(),
            r.taxon_group.as_str().to_string(),
            r.taxonomic_outlier.to_string(),
        ];
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}

/// Write the query TSV to a file.
pub fn write_query_tsv(path: &Path, records: &[Accession]) -> Result<()> {
    std::fs::write(path, render_query_tsv(records))
        .with_context(|| format!("writing {}", path.display()))
}

/// Read a query TSV back into accessions. Columns are located by header name, so
/// only `accession` is truly required; any missing optional column falls back to
/// a sensible default (0 / empty / ingroup / false).
pub fn read_query_tsv(path: &Path) -> Result<Vec<Accession>> {
    let content =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = content.lines();

    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("{}: empty query file (no header row)", path.display()))?;
    // Map column name -> position for this file's actual layout.
    let index: HashMap<&str, usize> = header
        .split('\t')
        .enumerate()
        .map(|(i, c)| (c, i))
        .collect();
    if !index.contains_key("accession") {
        bail!(
            "{}: query TSV missing required 'accession' column",
            path.display()
        );
    }

    let mut records = Vec::new();
    for (n, line) in lines.enumerate() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        // Closure over this row: fetch a column's value by name, or "" if the
        // column is absent / the row is short.
        let get = |name: &str| {
            index
                .get(name)
                .and_then(|&i| cols.get(i))
                .copied()
                .unwrap_or("")
        };

        let line_no = n + 2; // +1 for the header row, +1 for 1-based counting
        let parse_u64 = |name: &str| -> Result<u64> {
            let v = get(name);
            if v.is_empty() {
                return Ok(0);
            }
            v.parse()
                .with_context(|| format!("{}:{line_no}: bad {name} value {v:?}", path.display()))
        };

        records.push(Accession {
            accession: get("accession").to_string(),
            taxon_name: get("taxon_name").to_string(),
            taxid: parse_u64("taxid")?,
            length: parse_u64("length")? as usize,
            gene_annotation: get("gene_annotation").to_string(),
            refseq: get("refseq") == "true",
            source_db: get("source_db").to_string(),
            query_taxid: parse_u64("query_taxid")?,
            taxon_group: match get("taxon_group") {
                "" => TaxonGroup::Ingroup,
                g => {
                    TaxonGroup::parse(g).with_context(|| format!("{}:{line_no}", path.display()))?
                }
            },
            taxonomic_outlier: get("taxonomic_outlier") == "true",
        });
    }
    Ok(records)
}

/// Replace any tab/newline in a field with a space so it can't break the TSV row
/// structure.
fn clean(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Accession {
        Accession {
            accession: "AB123.1".to_string(),
            taxon_name: "Mus musculus".to_string(),
            taxid: 10090,
            length: 658,
            gene_annotation: "cytochrome c oxidase subunit I".to_string(),
            refseq: false,
            source_db: "GenBank".to_string(),
            query_taxid: 10088,
            taxon_group: TaxonGroup::Outgroup,
            taxonomic_outlier: false,
        }
    }

    #[test]
    fn tsv_round_trips() {
        let dir = std::env::temp_dir().join(format!("phorge_tsv_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.tsv");

        let a = sample();
        write_query_tsv(&path, std::slice::from_ref(&a)).unwrap();
        let back = read_query_tsv(&path).unwrap();

        assert_eq!(back.len(), 1);
        let b = &back[0];
        assert_eq!(b.accession, a.accession);
        assert_eq!(b.taxid, a.taxid);
        assert_eq!(b.length, a.length);
        assert_eq!(b.gene_annotation, a.gene_annotation);
        assert_eq!(b.refseq, a.refseq);
        assert_eq!(b.source_db, a.source_db);
        assert_eq!(b.query_taxid, a.query_taxid);
        assert_eq!(b.taxon_group, a.taxon_group);
    }

    #[test]
    fn reads_reordered_and_partial_columns() {
        // Simulates a user filtering with `cut`/`awk`: columns reordered, some
        // dropped. Only `accession` is required; missing columns default.
        let dir = std::env::temp_dir().join(format!("phorge_tsv_partial_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.tsv");
        std::fs::write(
            &path,
            "taxid\taccession\ttaxon_group\n10090\tAB999.1\toutgroup\n",
        )
        .unwrap();

        let recs = read_query_tsv(&path).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].accession, "AB999.1");
        assert_eq!(recs[0].taxon_group, TaxonGroup::Outgroup);
        assert_eq!(recs[0].taxid, 10090); // located by name despite reordering
        assert_eq!(recs[0].length, 0); // absent column -> default
    }

    #[test]
    fn missing_accession_column_errors() {
        let dir = std::env::temp_dir().join(format!("phorge_tsv_noacc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.tsv");
        std::fs::write(&path, "taxid\tlength\n10090\t658\n").unwrap();
        assert!(read_query_tsv(&path).is_err());
    }
}
