//! Peptides as the models read them: residue indices and, per position,
//! the element composition of its modifications (AlphaPeptDeep's input).

use tuffold_core::unimod;
use tuffold_core::{ModSite, Modification, PeptideForm};

/// Elements of the composition vector, in the models' order; the last slot
/// counts any other element.
const ELEMENTS: [&str; 109] = [
    "C", "H", "N", "O", "P", "S", "B", "F", "I", "K", "U", "V", "W", "X", "Y", "Ac", "Ag", "Al",
    "Am", "Ar", "As", "At", "Au", "Ba", "Be", "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf", "Cl", "Cm",
    "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Es", "Eu", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf",
    "Hg", "Ho", "In", "Ir", "Kr", "La", "Li", "Lr", "Lu", "Md", "Mg", "Mn", "Mo", "Na", "Nb", "Nd",
    "Ne", "Ni", "No", "Np", "Os", "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re",
    "Rh", "Rn", "Ru", "Sb", "Sc", "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti",
    "Tl", "Tm", "Xe", "Yb", "Zn", "Zr", "2H", "13C", "15N", "18O", "?",
];
pub(crate) const MOD_FEATURES: usize = ELEMENTS.len();

/// One peptide, padded with an empty position at each end (N- and
/// C-terminal modifications sit there).
#[derive(Debug, Clone, PartialEq)]
pub struct Encoded {
    /// `A` = 1 … `Z` = 26, 0 at the ends.
    pub residues: Vec<u32>,
    /// `residues.len() × MOD_FEATURES`, row-major.
    pub modifications: Vec<f32>,
}

impl Encoded {
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
}

/// `None` for a residue outside A–Z or a modification whose composition
/// isn't known (the models need the elements, not just the mass).
pub fn encode(form: &PeptideForm) -> Option<Encoded> {
    let length = form.sequence.len();
    let mut residues = Vec::with_capacity(length + 2);
    residues.push(0);
    for residue in form.sequence.bytes() {
        if !residue.is_ascii_uppercase() {
            return None;
        }
        residues.push(u32::from(residue - b'A') + 1);
    }
    residues.push(0);
    let mut modifications = vec![0.0; residues.len() * MOD_FEATURES];
    for modification in &form.modifications {
        let position = match modification.site {
            ModSite::NTerm => 0,
            ModSite::Residue(position) => usize::try_from(position).ok()?,
            ModSite::CTerm => length + 1,
        };
        let row = &mut modifications[position * MOD_FEATURES..(position + 1) * MOD_FEATURES];
        for (element, count) in known(modification)?.elements() {
            let slot = ELEMENTS
                .iter()
                .position(|known| *known == element)
                .unwrap_or(MOD_FEATURES - 1);
            row[slot] += count as f32;
        }
    }
    Some(Encoded {
        residues,
        modifications,
    })
}

/// The Unimod entry a modification is: by name when its mass agrees, else
/// by mass alone (engines often report only the mass).
fn known(modification: &Modification) -> Option<unimod::Known> {
    modification
        .name
        .as_deref()
        .and_then(unimod::by_name)
        .filter(|known| (known.mass - modification.mass_delta).abs() <= 0.01)
        .or_else(|| unimod::by_mass(modification.mass_delta))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn form(sequence: &str, modifications: Vec<Modification>) -> PeptideForm {
        PeptideForm {
            sequence: sequence.into(),
            modifications,
        }
    }

    fn modification(site: ModSite, mass_delta: f64) -> Modification {
        Modification {
            site,
            mass_delta,
            name: None,
            accession: None,
        }
    }

    #[test]
    fn residues_and_modification_elements_are_encoded_by_position() {
        let encoded = encode(&form(
            "ACK",
            vec![
                modification(ModSite::NTerm, 42.010_565),
                modification(ModSite::Residue(2), 57.021_464),
            ],
        ))
        .unwrap();
        assert_eq!(encoded.residues, [0, 1, 3, 11, 0]);
        let row = |position: usize| &encoded.modifications[position * 109..(position + 1) * 109];
        // Acetyl H(2) C(2) O on the N-terminus; carbamidomethyl H(3) C(2) N O on C.
        assert_eq!(row(0)[..4], [2.0, 2.0, 0.0, 1.0]);
        assert_eq!(row(2)[..4], [2.0, 3.0, 1.0, 1.0]);
        assert!(row(1).iter().all(|&v| v == 0.0));
    }

    #[test]
    fn unknown_modifications_and_residues_cannot_be_encoded() {
        assert!(encode(&form("PEPK", vec![modification(ModSite::Residue(1), 12.3)])).is_none());
        assert!(encode(&form("PEPk", vec![])).is_none());
    }
}
