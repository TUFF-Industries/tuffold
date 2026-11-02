//! In-silico tryptic digestion for abundance indices.
//!
//! emPAI (Ishihama et al. 2005) divides the peptides a protein showed by the
//! peptides it could show. "Could show" is approximated, as is common, by
//! the fully tryptic peptides without missed cleavages (cleaving after K or R
//! unless P follows) of [`OBSERVABLE_LENGTH`] residues: shorter ones are
//! rarely unique or retained, longer ones rarely fragment well enough.

use std::ops::RangeInclusive;

/// Lengths of tryptic peptides that count as observable.
pub const OBSERVABLE_LENGTH: RangeInclusive<usize> = 6..=30;

/// The fully tryptic peptides of `sequence`, without missed cleavages.
pub fn tryptic_peptides(sequence: &str) -> Vec<&str> {
    let residues: Vec<(usize, char)> = sequence.char_indices().collect();
    let mut peptides = Vec::new();
    let mut start = 0;
    for (i, &(at, residue)) in residues.iter().enumerate() {
        let next = residues.get(i + 1).map(|&(_, next)| next);
        if matches!(residue, 'K' | 'R') && next.is_some_and(|next| next != 'P') {
            let end = at + residue.len_utf8();
            peptides.push(&sequence[start..end]);
            start = end;
        }
    }
    if start < sequence.len() {
        peptides.push(&sequence[start..]);
    }
    peptides
}

/// How many distinct tryptic peptides of `sequence` are observable.
pub fn observable_peptides(sequence: &str) -> usize {
    let mut observable: Vec<&str> = tryptic_peptides(sequence)
        .into_iter()
        .filter(|peptide| OBSERVABLE_LENGTH.contains(&peptide.chars().count()))
        .collect();
    observable.sort_unstable();
    observable.dedup();
    observable.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trypsin_cleaves_after_k_and_r_except_before_p() {
        // M A K | P E P T I D E R | G G K | (end): K before P stays joined.
        assert_eq!(tryptic_peptides("MAKPEPTIDERGGK"), ["MAKPEPTIDER", "GGK"]);
        // A trailing K or R is the end, not an empty peptide.
        assert_eq!(tryptic_peptides("PEPTIDEK"), ["PEPTIDEK"]);
        assert_eq!(tryptic_peptides("KR"), ["K", "R"]);
        assert_eq!(tryptic_peptides(""), Vec::<&str>::new());
    }

    #[test]
    fn observable_peptides_are_distinct_and_six_to_thirty_long() {
        // MAKPEPTIDER (11) | GGK (3) | SAMPLEK (7) | SAMPLEK (7, repeat) |
        // 31 × A, ending in R (32) => observable: MAKPEPTIDER, SAMPLEK.
        let sequence = format!("MAKPEPTIDERGGKSAMPLEKSAMPLEK{}R", "A".repeat(31));
        assert_eq!(observable_peptides(&sequence), 2);
        // Exactly 6 and exactly 30 residues count.
        let edges = format!("ABCDEK{}R", "A".repeat(29));
        assert_eq!(observable_peptides(&edges), 2);
    }
}
