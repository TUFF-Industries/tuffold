//! Mass constants for peptide chemistry.
//!
//! Residue masses are of the amino acid *residue* (free amino acid minus
//! water), from the Unimod / NIST monoisotopic and average element masses.

/// Mass of a proton.
pub const PROTON: f64 = 1.007_276_466_621;
/// Monoisotopic mass of water (added once per peptide).
pub const WATER_MONO: f64 = 18.010_564_683_7;
/// Monoisotopic mass of the unmodified N-terminal group (H).
pub const N_TERM_MONO: f64 = 1.007_825_032_07;
/// Monoisotopic mass of the unmodified C-terminal group (OH).
pub const C_TERM_MONO: f64 = 17.002_739_651_6;
/// Average mass of the unmodified N-terminal group (H).
pub const N_TERM_AVERAGE: f64 = 1.007_94;
/// Average mass of the unmodified C-terminal group (OH).
pub const C_TERM_AVERAGE: f64 = 17.007_34;

/// Whether reported masses are monoisotopic or average.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MassType {
    #[default]
    Monoisotopic,
    Average,
}

/// Residue mass of a standard (or common non-standard) amino acid.
pub fn residue_mass(residue: char, mass_type: MassType) -> Option<f64> {
    let (mono, average) = match residue.to_ascii_uppercase() {
        'G' => (57.021_463_72, 57.051_3),
        'A' => (71.037_113_79, 71.078_8),
        'S' => (87.032_028_41, 87.078_2),
        'P' => (97.052_763_85, 97.115_2),
        'V' => (99.068_413_91, 99.132_6),
        'T' => (101.047_678_5, 101.104_7),
        'C' => (103.009_184_5, 103.143_9),
        'L' | 'I' => (113.084_064_0, 113.159_4),
        'N' => (114.042_927_4, 114.103_8),
        'D' => (115.026_943_1, 115.088_6),
        'Q' => (128.058_577_5, 128.130_7),
        'K' => (128.094_963_0, 128.174_1),
        'E' => (129.042_593_1, 129.115_5),
        'M' => (131.040_484_6, 131.192_6),
        'H' => (137.058_911_9, 137.141_1),
        'F' => (147.068_413_9, 147.176_6),
        'U' => (150.953_633_4, 150.038_6),
        'R' => (156.101_111_0, 156.187_5),
        'Y' => (163.063_328_6, 163.176_0),
        'W' => (186.079_312_9, 186.213_2),
        'O' => (237.147_726_9, 237.300_3),
        _ => return None,
    };
    Some(match mass_type {
        MassType::Monoisotopic => mono,
        MassType::Average => average,
    })
}

/// Monoisotopic neutral mass of an unmodified peptide, or `None` when it
/// contains a residue without a defined mass (e.g. `X`, `B`, `Z`).
pub fn peptide_mass(sequence: &str) -> Option<f64> {
    sequence
        .chars()
        .map(|residue| residue_mass(residue, MassType::Monoisotopic))
        .sum::<Option<f64>>()
        .map(|residues| residues + WATER_MONO)
}

/// m/z of a neutral mass at the given (positive) charge.
pub fn mz(neutral_mass: f64, charge: u8) -> f64 {
    let charge = f64::from(charge.max(1));
    (neutral_mass + charge * PROTON) / charge
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_mass_matches_reference() {
        // PEPTIDE monoisotopic neutral mass: 799.359964 Da (Unimod/Expasy).
        let mass = peptide_mass("PEPTIDE").unwrap();
        assert!((mass - 799.359_964).abs() < 1e-5, "{mass}");
        // SLNGEWR by hand: S 87.03203 + L 113.08406 + N 114.04293 + G 57.02146
        // + E 129.04259 + W 186.07931 + R 156.10111 + H2O 18.01056 = 860.41405;
        // + proton 1.00728 = 861.42133 at 1+.
        assert!((mz(peptide_mass("SLNGEWR").unwrap(), 1) - 861.421_33).abs() < 1e-4);
    }

    #[test]
    fn unknown_residues_have_no_mass() {
        assert_eq!(residue_mass('X', MassType::Monoisotopic), None);
        assert_eq!(peptide_mass("PEPXIDE"), None);
        assert_eq!(
            residue_mass('c', MassType::Monoisotopic),
            residue_mass('C', MassType::Monoisotopic)
        );
    }
}
