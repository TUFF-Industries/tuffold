//! Where peptides fall on a protein sequence, and how much they cover.
//!
//! Peptides are located by exact, ASCII case-insensitive substring search
//! (residue letters are ASCII; folding never changes a sequence's length): every
//! occurrence becomes a 1-based inclusive range, overlapping occurrences
//! included; absent peptides are skipped. No digestion, variants, or
//! modification interpretation.

use serde::{Deserialize, Serialize};

/// One 1-based, inclusive occurrence of a peptide on a protein sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceCoverage {
    pub start: u32,
    pub end: u32,
    pub peptide_id: String,
}

/// Locate every occurrence of each `(peptide_id, sequence)` pair.
///
/// Output follows input peptide order, each peptide's occurrences left to
/// right; positions count characters, 1-based and inclusive.
pub fn map_peptide_coverage<'a>(
    sequence: &str,
    peptides: impl IntoIterator<Item = (&'a str, &'a str)>,
) -> Vec<SequenceCoverage> {
    let folded = sequence.to_ascii_uppercase();
    let char_position = |byte: usize| -> u32 {
        let chars = if folded.is_ascii() {
            byte
        } else {
            folded[..byte].chars().count()
        };
        u32::try_from(chars).expect("protein sequences are shorter than u32::MAX residues")
    };
    let mut coverage = Vec::new();
    for (peptide_id, peptide) in peptides {
        let folded_peptide = peptide.to_ascii_uppercase();
        if folded_peptide.is_empty() {
            continue;
        }
        let length = u32::try_from(folded_peptide.chars().count())
            .expect("peptides are shorter than u32::MAX residues");
        let mut from = 0;
        while let Some(offset) = folded[from..].find(&folded_peptide) {
            let byte = from + offset;
            let start = char_position(byte) + 1;
            coverage.push(SequenceCoverage {
                start,
                end: start + length - 1,
                peptide_id: peptide_id.to_string(),
            });
            // Advance one character so overlapping occurrences are found.
            from = byte + folded[byte..].chars().next().map_or(1, char::len_utf8);
        }
    }
    coverage
}

/// Which residues (by character) at least one span covers.
pub fn covered_mask(length: usize, spans: &[SequenceCoverage]) -> Vec<bool> {
    let mut mask = vec![false; length];
    for span in spans {
        let start = usize::try_from(span.start)
            .unwrap_or(usize::MAX)
            .saturating_sub(1);
        let end = usize::try_from(span.end).unwrap_or(usize::MAX).min(length);
        for covered in mask.iter_mut().take(end).skip(start) {
            *covered = true;
        }
    }
    mask
}

/// Fraction of residues covered, 0 for an empty sequence.
pub fn covered_fraction(length: usize, spans: &[SequenceCoverage]) -> f64 {
    if length == 0 {
        return 0.0;
    }
    let covered = covered_mask(length, spans).iter().filter(|&&c| c).count();
    covered as f64 / length as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(coverage: &[SequenceCoverage]) -> Vec<(u32, u32, &str)> {
        coverage
            .iter()
            .map(|entry| (entry.start, entry.end, entry.peptide_id.as_str()))
            .collect()
    }

    #[test]
    fn coverage_finds_all_occurrences_left_to_right() {
        let coverage = map_peptide_coverage("MPEPTIDEKAGPEPTIDEKG", [("pep-PEPTIDEK", "PEPTIDEK")]);
        assert_eq!(
            spans(&coverage),
            [(2, 9, "pep-PEPTIDEK"), (12, 19, "pep-PEPTIDEK")]
        );
    }

    #[test]
    fn coverage_includes_overlapping_occurrences() {
        let coverage = map_peptide_coverage("AAAAA", [("pep-AAA", "AAA")]);
        assert_eq!(
            spans(&coverage),
            [(1, 3, "pep-AAA"), (2, 4, "pep-AAA"), (3, 5, "pep-AAA")]
        );
    }

    #[test]
    fn coverage_is_case_insensitive_and_follows_input_order() {
        let coverage =
            map_peptide_coverage("mpeptidekgag", [("pep-B", "GAG"), ("pep-A", "peptidek")]);
        assert_eq!(spans(&coverage), [(10, 12, "pep-B"), (2, 9, "pep-A")]);
    }

    #[test]
    fn coverage_skips_absent_and_empty_peptides() {
        assert!(map_peptide_coverage("MAGTAST", [("pep-X", "PEPTIDEK")]).is_empty());
        assert!(map_peptide_coverage("MAGTAST", [("pep-X", "")]).is_empty());
    }

    #[test]
    fn coverage_counts_characters_not_bytes() {
        let coverage = map_peptide_coverage("ÄAPEP", [("pep-PEP", "PEP")]);
        assert_eq!(spans(&coverage), [(3, 5, "pep-PEP")]);
    }

    #[test]
    fn coverage_counts_each_residue_once() {
        // MPEPTIDEKAGPEPTIDEKG (20): PEPTIDEK at 2-9 and 12-19, plus AGP at
        // 10-12 overlapping the second: covered 2..=19 = 18 residues.
        let sequence = "MPEPTIDEKAGPEPTIDEKG";
        let spans = map_peptide_coverage(sequence, [("a", "PEPTIDEK"), ("b", "AGP")]);
        let mask = covered_mask(sequence.len(), &spans);
        assert!(!mask[0] && mask[1] && mask[18] && !mask[19]);
        assert!((covered_fraction(sequence.len(), &spans) - 18.0 / 20.0).abs() < 1e-12);
        assert_eq!(covered_fraction(0, &[]), 0.0);
    }
}
