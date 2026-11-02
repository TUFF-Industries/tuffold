//! Peptide forms: a stripped sequence plus its modifications.
//!
//! Two PSMs identify the same peptide form when the sequence and every
//! modification (site and mass delta) agree. The canonical key is a
//! ProForma-style string (`PEPT[+79.9663]IDEK`, `[+42.0106]-PEPTIDEK`,
//! `PEPTIDEK-[+0.9840]`). In the key only, deltas within
//! [`unimod::SNAP_TOLERANCE`] of a known Unimod modification take its exact
//! value, so engines reporting 15.995 and 15.994915 share a key; every other
//! delta is rounded to [`DELTA_DECIMALS`]. The reported delta itself is kept
//! unchanged on the [`Modification`].

use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::unimod;

/// Mass deltas are compared and keyed at this many decimals (0.1 mDa).
pub const DELTA_DECIMALS: usize = 4;

/// Below this magnitude a reported "modification" is no modification (e.g.
/// pepXML terminal masses equal to the unmodified group).
pub const NULL_DELTA: f64 = 0.005;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind", content = "index")]
pub enum ModSite {
    NTerm,
    /// 1-based residue position.
    Residue(u32),
    CTerm,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Modification {
    pub site: ModSite,
    /// Monoisotopic mass difference from the unmodified residue or terminus.
    pub mass_delta: f64,
    /// Name as reported by the source (e.g. "Oxidation"), when it gave one.
    pub name: Option<String>,
    /// Unimod or PSI-MS accession, when the source gave one.
    pub accession: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PeptideForm {
    /// Upper-case stripped sequence.
    pub sequence: String,
    /// Sorted by site (N-term, residues left to right, C-term).
    pub modifications: Vec<Modification>,
}

/// A modification that cannot belong to the peptide it was reported on.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("modification at residue {position} is outside peptide {sequence} ({length} residues)")]
pub struct InvalidModification {
    pub sequence: String,
    pub position: u32,
    pub length: usize,
}

impl PeptideForm {
    /// Build a form, dropping null deltas and ordering modifications by site.
    ///
    /// A residue modification outside `1..=length` is rejected: it would
    /// otherwise vanish from the key and merge with the unmodified form.
    pub fn new(
        sequence: &str,
        modifications: Vec<Modification>,
    ) -> Result<Self, InvalidModification> {
        let sequence = sequence.to_ascii_uppercase();
        let length = sequence.chars().count();
        for modification in &modifications {
            if let ModSite::Residue(position) = modification.site {
                let inside = usize::try_from(position).is_ok_and(|p| (1..=length).contains(&p));
                if !inside {
                    return Err(InvalidModification {
                        sequence,
                        position,
                        length,
                    });
                }
            }
        }
        let mut modifications: Vec<Modification> = modifications
            .into_iter()
            .filter(|modification| modification.mass_delta.abs() >= NULL_DELTA)
            .collect();
        modifications.sort_by(|a, b| {
            a.site
                .cmp(&b.site)
                .then(a.mass_delta.total_cmp(&b.mass_delta))
        });
        Ok(Self {
            sequence,
            modifications,
        })
    }

    pub fn unmodified(sequence: &str) -> Self {
        Self {
            sequence: sequence.to_ascii_uppercase(),
            modifications: Vec::new(),
        }
    }

    /// The canonical ProForma-style key.
    pub fn proforma(&self) -> String {
        let deltas = |site: ModSite| -> String {
            self.modifications
                .iter()
                .filter(|modification| modification.site == site)
                .map(|modification| format_delta(modification.mass_delta))
                .collect()
        };
        let mut key = String::with_capacity(self.sequence.len() + 12);
        let n_term = deltas(ModSite::NTerm);
        if !n_term.is_empty() {
            let _ = write!(key, "{n_term}-");
        }
        for (index, residue) in self.sequence.chars().enumerate() {
            key.push(residue);
            let position = u32::try_from(index + 1).expect("peptides are short");
            key.push_str(&deltas(ModSite::Residue(position)));
        }
        let c_term = deltas(ModSite::CTerm);
        if !c_term.is_empty() {
            let _ = write!(key, "-{c_term}");
        }
        key
    }
}

fn format_delta(delta: f64) -> String {
    let keyed = unimod::snap(delta).unwrap_or(delta);
    format!("[{keyed:+.prec$}]", prec = DELTA_DECIMALS)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn modification(site: ModSite, mass_delta: f64) -> Modification {
        Modification {
            site,
            mass_delta,
            name: None,
            accession: None,
        }
    }

    #[test]
    fn unmodified_key_is_the_sequence() {
        assert_eq!(PeptideForm::unmodified("peptidek").proforma(), "PEPTIDEK");
    }

    #[test]
    fn key_places_residue_and_terminal_deltas() {
        let form = PeptideForm::new(
            "PEPTIDEK",
            vec![
                modification(ModSite::CTerm, 0.984_016),
                modification(ModSite::Residue(4), 79.966_331),
                modification(ModSite::NTerm, 42.010_565),
            ],
        )
        .unwrap();
        assert_eq!(form.proforma(), "[+42.0106]-PEPT[+79.9663]IDEK-[+0.9840]");
    }

    #[test]
    fn engines_with_different_precision_share_a_key() {
        // 3-decimal (15.995), 5-decimal, and 7-decimal reports of Oxidation.
        let keys: Vec<String> = [15.995, 15.99491, 15.994_915_1]
            .map(|delta| {
                PeptideForm::new("MK", vec![modification(ModSite::Residue(1), delta)])
                    .unwrap()
                    .proforma()
            })
            .to_vec();
        assert_eq!(keys, ["M[+15.9949]K", "M[+15.9949]K", "M[+15.9949]K"]);
        // The reported delta itself is not rewritten.
        let form = PeptideForm::new("MK", vec![modification(ModSite::Residue(1), 15.995)]).unwrap();
        assert_eq!(form.modifications[0].mass_delta, 15.995);
    }

    #[test]
    fn null_deltas_are_dropped() {
        // pepXML reports an unmodified N-terminus as mass 1.0079 (delta ~0).
        let form =
            PeptideForm::new("PEPTIDE", vec![modification(ModSite::NTerm, 0.000_075)]).unwrap();
        assert!(form.modifications.is_empty());
        assert_eq!(form.proforma(), "PEPTIDE");
    }

    #[test]
    fn negative_deltas_keep_their_sign() {
        let form = PeptideForm::new(
            "QPEPTIDE",
            vec![modification(ModSite::Residue(1), -17.026_549)],
        )
        .unwrap();
        assert_eq!(form.proforma(), "Q[-17.0265]PEPTIDE");
    }

    #[test]
    fn residue_positions_outside_the_peptide_are_rejected() {
        for position in [0, 9, 100] {
            let error = PeptideForm::new(
                "PEPTIDEK",
                vec![modification(ModSite::Residue(position), 42.0106)],
            )
            .unwrap_err();
            assert_eq!((error.position, error.length), (position, 8));
        }
        assert!(
            PeptideForm::new("PEPTIDEK", vec![modification(ModSite::Residue(8), 42.0106)]).is_ok()
        );
    }

    #[test]
    fn unknown_deltas_round_to_four_decimals() {
        let form =
            PeptideForm::new("PEPK", vec![modification(ModSite::Residue(2), 123.456_789)]).unwrap();
        assert_eq!(form.proforma(), "PE[+123.4568]PK");
    }
}
