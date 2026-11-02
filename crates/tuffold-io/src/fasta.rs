//! Protein FASTA discovery and parsing.
//!
//! Conventions (deterministic, no inference):
//! - Discovery scans a directory's top level only, accepting `.fas`, `.fasta`,
//!   and `.fa` case-insensitively; candidates sort by case-insensitive file
//!   name and the first wins. A missing FASTA must surface as "no sequence
//!   information", never as fabricated sequences.
//! - Records are `>`-delimited; wrapped sequence lines are joined, all
//!   whitespace is removed, and sequences are stored upper-cased.
//! - Identifiers are the first whitespace-delimited header token plus each of
//!   its `|`-separated segments, so `sp|P00001|P1_HUMAN` registers all four
//!   spellings.
//! - Entries whose first token starts with a reversed-database label (`rev`,
//!   `decoy`, `xxx`, `dec`, `reversed`) followed by a non-alphanumeric
//!   separator are decoys; a forward entry wins any identifier it shares with a
//!   decoy.
//!
//! Peptide coverage of a sequence is `tuffold_core::coverage`.

use std::collections::HashMap;
use std::io;
use std::path::{Path, PathBuf};

use tuffold_core::decoy::is_decoy_accession;

pub const FASTA_EXTENSIONS: &[&str] = &["fas", "fasta", "fa"];
/// Experiment metadata key recording the selected database's file name.
pub const FASTA_METADATA_KEY: &str = "fasta_file";

/// Every top-level FASTA candidate in `dir`, in selection order.
pub fn find_fasta_candidates(dir: &Path) -> io::Result<Vec<PathBuf>> {
    let mut candidates = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        let is_fasta = path.extension().is_some_and(|ext| {
            let ext = ext.to_string_lossy();
            FASTA_EXTENSIONS
                .iter()
                .any(|known| ext.eq_ignore_ascii_case(known))
        });
        if is_fasta && path.is_file() {
            candidates.push(path);
        }
    }
    candidates.sort_by_cached_key(|path| {
        let name = crate::tsv::file_name(path);
        (name.to_lowercase(), name)
    });
    Ok(candidates)
}

/// The selected FASTA in `dir`, or `None` when it has none.
pub fn find_fasta(dir: &Path) -> io::Result<Option<PathBuf>> {
    Ok(find_fasta_candidates(dir)?.into_iter().next())
}

/// One FASTA entry: its identifier spellings, sequence, and decoy flag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub identifiers: Vec<String>,
    /// The header after the first token (`OS=...` and all), when present.
    pub description: Option<String>,
    pub sequence: String,
    pub is_decoy: bool,
}

impl FastaRecord {
    fn from_parts(header: &str, chunks: &[&str]) -> Self {
        let identifiers = header_identifiers(header);
        let description = header
            .split_once(char::is_whitespace)
            .map(|(_, rest)| rest.trim().to_string())
            .filter(|rest| !rest.is_empty());
        let sequence: String = chunks
            .iter()
            .flat_map(|chunk| chunk.chars())
            .filter(|ch| !ch.is_whitespace())
            .collect::<String>()
            .to_uppercase();
        let is_decoy = identifiers
            .first()
            .is_some_and(|first| is_decoy_accession(first));
        Self {
            identifiers,
            description,
            sequence,
            is_decoy,
        }
    }
}

fn header_identifiers(header: &str) -> Vec<String> {
    let Some(first) = header.split_whitespace().next() else {
        return Vec::new();
    };
    let mut identifiers = vec![first.to_string()];
    for segment in first.split('|').filter(|segment| !segment.is_empty()) {
        if !identifiers.iter().any(|known| known == segment) {
            identifiers.push(segment.to_string());
        }
    }
    identifiers
}

/// Parse FASTA text, tolerating wrapped lines and junk before the first header.
pub fn parse_fasta_str(text: &str) -> Vec<FastaRecord> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut records = Vec::new();
    let mut header: Option<&str> = None;
    let mut chunks: Vec<&str> = Vec::new();
    for line in text.lines().map(str::trim) {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(previous) = header.replace(rest.trim()) {
                records.push(FastaRecord::from_parts(previous, &chunks));
            }
            chunks.clear();
        } else if !line.is_empty() && header.is_some() {
            chunks.push(line);
        }
    }
    if let Some(last) = header {
        records.push(FastaRecord::from_parts(last, &chunks));
    }
    records
}

/// Read and parse a FASTA file; invalid UTF-8 is an `InvalidData` error.
pub fn parse_fasta(path: &Path) -> io::Result<Vec<FastaRecord>> {
    Ok(parse_fasta_str(&std::fs::read_to_string(path)?))
}

/// Map every identifier spelling to its record position, forward preferred.
///
/// First-seen wins, except that a forward record takes over an identifier
/// currently held by a decoy.
pub fn build_identifier_index(records: &[FastaRecord]) -> HashMap<String, usize> {
    let mut index: HashMap<String, usize> = HashMap::new();
    for (position, record) in records.iter().enumerate() {
        for identifier in &record.identifiers {
            let takes_over = match index.get(identifier) {
                None => true,
                Some(&current) => records[current].is_decoy && !record.is_decoy,
            };
            if takes_over {
                index.insert(identifier.clone(), position);
            }
        }
    }
    index
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(record: &FastaRecord) -> Vec<&str> {
        record.identifiers.iter().map(String::as_str).collect()
    }

    #[test]
    fn joins_wrapped_lines_and_folds_case() {
        let records =
            parse_fasta_str(">sp|P1|NAMES_HUMAN Protein one\nMAGTA\nSTPGVQ  \n\nAGAPEP\n");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence, "MAGTASTPGVQAGAPEP");
        assert_eq!(records[0].description.as_deref(), Some("Protein one"));
        assert!(!records[0].is_decoy);
        assert_eq!(
            ids(&records[0]),
            ["sp|P1|NAMES_HUMAN", "sp", "P1", "NAMES_HUMAN"]
        );
    }

    #[test]
    fn handles_plain_and_pipe_headers() {
        let records =
            parse_fasta_str(">P00001 plain accession\nMS\n>P00002\nGT\n>sp|P3|N3_HUMAN desc\nAA\n");
        let descriptions: Vec<Option<&str>> =
            records.iter().map(|r| r.description.as_deref()).collect();
        assert_eq!(descriptions, [Some("plain accession"), None, Some("desc")]);
        let all: Vec<Vec<&str>> = records.iter().map(ids).collect();
        assert_eq!(
            all,
            [
                vec!["P00001"],
                vec!["P00002"],
                vec!["sp|P3|N3_HUMAN", "sp", "P3", "N3_HUMAN"]
            ]
        );
    }

    #[test]
    fn ignores_leading_junk_and_keeps_empty_records() {
        let records = parse_fasta_str("not a header\n>P1\nMS\n>\n>P2\nGT\n");
        let all: Vec<Vec<&str>> = records.iter().map(ids).collect();
        assert_eq!(all, [vec!["P1"], vec![], vec!["P2"]]);
        let sequences: Vec<&str> = records.iter().map(|r| r.sequence.as_str()).collect();
        assert_eq!(sequences, ["MS", "", "GT"]);
    }

    #[test]
    fn byte_order_mark_does_not_hide_the_first_record() {
        let records = parse_fasta_str("\u{feff}>sp|P1|N1 first\nMS\n>P2\nGT\n");
        assert_eq!(records.len(), 2);
        assert_eq!(ids(&records[0]), ["sp|P1|N1", "sp", "P1", "N1"]);
        assert_eq!(records[0].sequence, "MS");
    }

    #[test]
    fn records_follow_the_one_decoy_rule() {
        // tuffold_core::decoy is the rule; this checks the parser uses it,
        // TPP's numbered `reverse-` decoys and Scaffold's suffix included.
        let records = parse_fasta_str(
            ">rev_sp|P1|N1 a\nVV\n>reverse-258\nWW\n>P15975-DECOY\nKK\n>revolver_P1\nMS\n>P00001\nGT\n",
        );
        let decoys: Vec<bool> = records.iter().map(|r| r.is_decoy).collect();
        assert_eq!(decoys, [true, true, true, false, false]);
    }

    #[test]
    fn identifier_index_prefers_forward_over_decoy() {
        // Decoy listed first and sharing the accession, as in FragPipe databases.
        let records = parse_fasta_str(
            ">rev_sp|P1|N1_HUMAN decoy\nVV\n>sp|P1|N1_HUMAN forward\nMS\n>rev_sp|P2|N2_HUMAN decoy\nWW\n",
        );
        let index = build_identifier_index(&records);
        assert_eq!(records[index["P1"]].sequence, "MS");
        assert!(!records[index["P1"]].is_decoy);
        assert!(records[index["P2"]].is_decoy);
        assert_eq!(records[index["N1_HUMAN"]].sequence, "MS");
        // The decoy's own full first token still resolves to the decoy.
        assert!(records[index["rev_sp|P2|N2_HUMAN"]].is_decoy);
    }

    #[test]
    fn identifier_index_first_seen_wins_between_forward_records() {
        let index = build_identifier_index(&parse_fasta_str(">P1 alpha\nMS\n>P1 beta\nGT\n"));
        assert_eq!(index["P1"], 0);
    }

    #[test]
    fn discovery_orders_case_insensitively_and_ignores_subfolders() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        for name in ["z_run.fa", "a_database.FAS", "m.fasta"] {
            std::fs::write(root.join(name), ">A\nMS\n").unwrap();
        }
        std::fs::write(root.join("notes.txt"), "not a fasta").unwrap();
        std::fs::create_dir(root.join("control-1")).unwrap();
        std::fs::write(root.join("control-1/nested.fas"), ">A\nMS\n").unwrap();
        assert_eq!(find_fasta(root).unwrap(), Some(root.join("a_database.FAS")));
    }

    #[test]
    fn discovery_returns_none_without_candidates() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("psm.tsv"), "Spectrum\n").unwrap();
        std::fs::write(dir.path().join("protein.txt"), ">A\nMS\n").unwrap();
        assert_eq!(find_fasta(dir.path()).unwrap(), None);
    }
}
