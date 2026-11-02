//! Monoisotopic mass deltas for common Unimod modifications.
//!
//! Used when a source names a modification by accession without stating its
//! mass (mzIdentML allows this). Unknown accessions return `None` and the
//! importer reports them instead of guessing. Values are Unimod's
//! monoisotopic deltas (unimod.org).

/// (accession number, name, monoisotopic delta, Unimod delta composition)
const TABLE: &[(u32, &str, f64, &str)] = &[
    (1, "Acetyl", 42.010_565, "H(2) C(2) O"),
    (4, "Carbamidomethyl", 57.021_464, "H(3) C(2) N O"),
    (5, "Carbamyl", 43.005_814, "H C N O"),
    (7, "Deamidated", 0.984_016, "H(-1) N(-1) O"),
    (21, "Phospho", 79.966_331, "H O(3) P"),
    (27, "Glu->pyro-Glu", -18.010_565, "H(-2) O(-1)"),
    (28, "Gln->pyro-Glu", -17.026_549, "H(-3) N(-1)"),
    (34, "Methyl", 14.015_650, "H(2) C"),
    (35, "Oxidation", 15.994_915, "O"),
    (36, "Dimethyl", 28.031_300, "H(4) C(2)"),
    (39, "Methylthio", 45.987_721, "H(2) C S"),
    (121, "GG", 114.042_927, "H(6) C(4) N(2) O(2)"),
    (214, "iTRAQ4plex", 144.102_063, "H(12) C(4) 13C(3) N 15N O"),
    (
        730,
        "iTRAQ8plex",
        304.205_360,
        "H(24) C(7) 13C(7) N(3) 15N O(3)",
    ),
    (737, "TMT6plex", 229.162_932, "H(20) C(8) 13C(4) N 15N O(2)"),
    (
        2016,
        "TMTpro",
        304.207_146,
        "H(25) C(8) 13C(7) N 15N(2) O(3)",
    ),
];

/// Deltas this close to a table entry are the same modification reported at
/// lower precision (3-decimal engines are within 0.0005). Small enough to keep
/// TMTpro (304.2071) and iTRAQ8plex (304.2054) apart.
pub const SNAP_TOLERANCE: f64 = 0.001;

/// The table delta nearest to `delta` within [`SNAP_TOLERANCE`], if any.
pub fn snap(delta: f64) -> Option<f64> {
    TABLE
        .iter()
        .map(|(_, _, known, _)| *known)
        .filter(|known| (known - delta).abs() <= SNAP_TOLERANCE)
        .min_by(|a, b| (a - delta).abs().total_cmp(&(b - delta).abs()))
}

/// Mass delta for an accession like `UNIMOD:35`.
pub fn mass_delta(accession: &str) -> Option<f64> {
    let number: u32 = accession
        .strip_prefix("UNIMOD:")
        .or_else(|| accession.strip_prefix("unimod:"))?
        .parse()
        .ok()?;
    TABLE
        .iter()
        .find(|(id, _, _, _)| *id == number)
        .map(|(_, _, delta, _)| *delta)
}

/// A modification in the table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Known {
    pub id: u32,
    pub name: &'static str,
    pub mass: f64,
    /// Elements gained or lost, Unimod style: `H(3) C(2) N O`.
    pub composition: &'static str,
}

impl Known {
    /// `UNIMOD:35`.
    pub fn accession(&self) -> String {
        format!("UNIMOD:{}", self.id)
    }
}

fn known(&(id, name, mass, composition): &(u32, &'static str, f64, &'static str)) -> Known {
    Known {
        id,
        name,
        mass,
        composition,
    }
}

impl Known {
    /// The composition as (element, count) pairs: `H(-1) N(-1) O` gives
    /// `[("H", -1), ("N", -1), ("O", 1)]`; isotopes keep their mass number
    /// (`13C`).
    pub fn elements(&self) -> Vec<(&'static str, i32)> {
        self.composition
            .split_whitespace()
            .filter_map(|part| match part.split_once('(') {
                Some((element, count)) => Some((element, count.strip_suffix(')')?.parse().ok()?)),
                None => Some((part, 1)),
            })
            .collect()
    }
}

/// The modification named `name` (case-insensitive: `Oxidation`, `tmt6plex`).
pub fn by_name(name: &str) -> Option<Known> {
    TABLE
        .iter()
        .find(|(_, known, _, _)| known.eq_ignore_ascii_case(name))
        .map(known)
}

/// The modification whose delta is within [`SNAP_TOLERANCE`] of `delta`.
pub fn by_mass(delta: f64) -> Option<Known> {
    let mass = snap(delta)?;
    TABLE
        .iter()
        .find(|(_, _, known, _)| *known == mass)
        .map(known)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modifications_are_found_by_name_and_by_mass() {
        let tmt = by_name("tmt6plex").unwrap();
        assert_eq!(
            (tmt.name, tmt.accession()),
            ("TMT6plex", "UNIMOD:737".to_string())
        );
        assert_eq!(by_mass(15.9949).unwrap().name, "Oxidation");
        // TMTpro and iTRAQ8plex differ by 0.0018 and stay apart.
        assert_eq!(by_mass(304.2071).unwrap().name, "TMTpro");
        assert_eq!(by_mass(304.2054).unwrap().name, "iTRAQ8plex");
        assert_eq!(by_name("Nonsense"), None);
        assert_eq!(by_mass(1.5), None);
    }

    #[test]
    fn compositions_list_elements_with_counts() {
        assert_eq!(
            by_name("Deamidated").unwrap().elements(),
            [("H", -1), ("N", -1), ("O", 1)]
        );
        assert_eq!(by_name("TMT6plex").unwrap().elements()[2], ("13C", 4));
    }

    #[test]
    fn oxidation_is_one_oxygen() {
        // 16O monoisotopic mass: 15.9949146.
        assert!((mass_delta("UNIMOD:35").unwrap() - 15.994_914_6).abs() < 1e-6);
    }

    #[test]
    fn phospho_is_hpo3() {
        // H 1.00782503 + P 30.97376199 + 3 x O 15.99491462 = 79.96633088.
        assert!((mass_delta("UNIMOD:21").unwrap() - 79.966_330_9).abs() < 1e-6);
    }

    #[test]
    fn snapping_unifies_precision_but_keeps_close_labels_apart() {
        assert_eq!(snap(15.995), Some(15.994_915));
        assert_eq!(snap(15.994_915_1), Some(15.994_915));
        assert_eq!(snap(304.206), Some(304.205_360), "nearer to iTRAQ8plex");
        assert_eq!(snap(304.2071), Some(304.207_146), "TMTpro");
        assert_eq!(snap(15.99), None, "beyond tolerance");
        assert_eq!(snap(123.456), None);
    }

    #[test]
    fn unknown_or_malformed_accessions_have_no_delta() {
        assert_eq!(mass_delta("UNIMOD:999999"), None);
        assert_eq!(mass_delta("MS:1001460"), None);
        assert_eq!(mass_delta("UNIMOD:abc"), None);
    }
}
