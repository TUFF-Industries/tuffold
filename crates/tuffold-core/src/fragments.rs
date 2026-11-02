//! Theoretical b and y fragment ions, and matching them to measured peaks.
//!
//! A b ion holds the first `n` residues (plus N-terminal modifications); a y
//! ion the last `n` (plus C-terminal modifications and water). At charge `z`:
//! `b = (Σ residues + z · H⁺) / z`, `y = (Σ residues + H₂O + z · H⁺) / z`,
//! monoisotopic, with each residue's modification deltas added to it.
//! Fragments are computed up to charge 2, never above the precursor's charge
//! minus one (a 2+ precursor gives 1+ fragments only).

use crate::chemistry::{MassType, PROTON, WATER_MONO, residue_mass};
use crate::peptide::{ModSite, PeptideForm};

/// Highest fragment charge considered.
pub const MAX_FRAGMENT_CHARGE: u8 = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum IonSeries {
    B,
    Y,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentIon {
    pub series: IonSeries,
    /// Residues in the fragment (b3 holds the first three).
    pub ordinal: u32,
    pub charge: u8,
    pub mz: f64,
}

impl FragmentIon {
    /// `b3`, `y7++`.
    pub fn label(&self) -> String {
        let series = match self.series {
            IonSeries::B => 'b',
            IonSeries::Y => 'y',
        };
        let charges = if self.charge > 1 {
            "+".repeat(usize::from(self.charge))
        } else {
            String::new()
        };
        format!("{series}{}{charges}", self.ordinal)
    }
}

/// Fragment charges for a precursor charge: `1..=min(z - 1, 2)`, at least 1.
pub fn fragment_charges(precursor_charge: Option<u8>) -> std::ops::RangeInclusive<u8> {
    let highest = precursor_charge
        .map_or(1, |z| z.saturating_sub(1))
        .clamp(1, MAX_FRAGMENT_CHARGE);
    1..=highest
}

/// Every b and y ion (ordinals `1..len`, charges from [`fragment_charges`]),
/// b before y, by ordinal then charge. `None` when a residue has no mass.
pub fn fragment_ions(form: &PeptideForm, precursor_charge: Option<u8>) -> Option<Vec<FragmentIon>> {
    let residues = residue_masses(form)?;
    let (n_term, c_term) = terminal_deltas(form);
    let length = residues.len();
    let mut ions = Vec::new();
    for (series, offset) in [(IonSeries::B, n_term), (IonSeries::Y, c_term + WATER_MONO)] {
        for ordinal in 1..length {
            let neutral: f64 = match series {
                IonSeries::B => residues[..ordinal].iter().sum(),
                IonSeries::Y => residues[length - ordinal..].iter().sum(),
            };
            for charge in fragment_charges(precursor_charge) {
                let z = f64::from(charge);
                ions.push(FragmentIon {
                    series,
                    ordinal: u32::try_from(ordinal).ok()?,
                    charge,
                    mz: (neutral + offset + z * PROTON) / z,
                });
            }
        }
    }
    Some(ions)
}

/// Each residue's monoisotopic mass with its modification deltas.
fn residue_masses(form: &PeptideForm) -> Option<Vec<f64>> {
    let mut masses: Vec<f64> = form
        .sequence
        .chars()
        .map(|residue| residue_mass(residue, MassType::Monoisotopic))
        .collect::<Option<_>>()?;
    for modification in &form.modifications {
        if let ModSite::Residue(position) = modification.site {
            let index = usize::try_from(position).ok()?.checked_sub(1)?;
            *masses.get_mut(index)? += modification.mass_delta;
        }
    }
    Some(masses)
}

fn terminal_deltas(form: &PeptideForm) -> (f64, f64) {
    form.modifications
        .iter()
        .fold((0.0, 0.0), |(n, c), modification| match modification.site {
            ModSite::NTerm => (n + modification.mass_delta, c),
            ModSite::CTerm => (n, c + modification.mass_delta),
            ModSite::Residue(_) => (n, c),
        })
}

/// How close a peak must be to a theoretical m/z.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    Ppm(f64),
    Da(f64),
}

impl Tolerance {
    pub fn window(self, mz: f64) -> f64 {
        match self {
            Self::Ppm(ppm) => mz * ppm * 1e-6,
            Self::Da(da) => da,
        }
    }
}

/// A theoretical ion matched to a measured peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Match {
    pub ion: FragmentIon,
    /// Index of the peak in the spectrum.
    pub peak: usize,
    /// `(observed - theoretical) / theoretical`, in ppm.
    pub error_ppm: f64,
}

/// Match each ion to the most intense peak within tolerance (ties: the
/// closer one). `mz` must be sorted ascending; `intensity` is parallel.
pub fn annotate(
    ions: &[FragmentIon],
    mz: &[f64],
    intensity: &[f32],
    tolerance: Tolerance,
) -> Vec<Match> {
    assert_eq!(
        mz.len(),
        intensity.len(),
        "m/z and intensity differ in length"
    );
    ions.iter()
        .filter_map(|ion| {
            let window = tolerance.window(ion.mz);
            let first = mz.partition_point(|&peak| peak < ion.mz - window);
            let last = mz.partition_point(|&peak| peak <= ion.mz + window);
            (first..last)
                .max_by(|&a, &b| {
                    intensity[a]
                        .total_cmp(&intensity[b])
                        .then_with(|| (mz[b] - ion.mz).abs().total_cmp(&(mz[a] - ion.mz).abs()))
                })
                .map(|peak| Match {
                    ion: *ion,
                    peak,
                    error_ppm: (mz[peak] - ion.mz) / ion.mz * 1e6,
                })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::peptide::Modification;

    fn close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-5, "{actual} != {expected}");
    }

    fn ion(ions: &[FragmentIon], series: IonSeries, ordinal: u32, charge: u8) -> f64 {
        ions.iter()
            .find(|i| i.series == series && i.ordinal == ordinal && i.charge == charge)
            .unwrap()
            .mz
    }

    #[test]
    fn b_and_y_ions_of_peptide() {
        // P 97.052764, E 129.042593, D 115.026943; H+ 1.007276; H2O 18.010565.
        let ions = fragment_ions(&PeptideForm::unmodified("PEPTIDE"), Some(2)).unwrap();
        // 2+ precursor: 1+ fragments only; 6 b + 6 y.
        assert_eq!(ions.len(), 12);
        close(
            ion(&ions, IonSeries::B, 2, 1),
            97.052_764 + 129.042_593 + 1.007_276,
        );
        close(
            ion(&ions, IonSeries::Y, 1, 1),
            129.042_593 + 18.010_565 + 1.007_276,
        );
        close(
            ion(&ions, IonSeries::Y, 2, 1),
            115.026_943 + 129.042_593 + 18.010_565 + 1.007_276,
        );
    }

    #[test]
    fn modifications_shift_the_fragments_that_hold_them() {
        // Oxidized M1 and acetylated N-terminus: b1 = M 131.040485 + 15.994915
        // + 42.010565 + H+; y1 (K) is untouched: 128.094963 + H2O + H+.
        let form = PeptideForm::new(
            "MPEPTIDEK",
            vec![
                Modification {
                    site: ModSite::NTerm,
                    mass_delta: 42.010_565,
                    name: None,
                    accession: None,
                },
                Modification {
                    site: ModSite::Residue(1),
                    mass_delta: 15.994_915,
                    name: None,
                    accession: None,
                },
            ],
        )
        .unwrap();
        let ions = fragment_ions(&form, Some(3)).unwrap();
        close(
            ion(&ions, IonSeries::B, 1, 1),
            131.040_485 + 15.994_915 + 42.010_565 + 1.007_276,
        );
        close(
            ion(&ions, IonSeries::Y, 1, 1),
            128.094_963 + 18.010_565 + 1.007_276,
        );
        // 3+ precursor: 2+ fragments too; y1++ = (128.094963 + 18.010565 +
        // 2 x 1.007276) / 2.
        close(
            ion(&ions, IonSeries::Y, 1, 2),
            (128.094_963 + 18.010_565 + 2.0 * 1.007_276) / 2.0,
        );
        assert_eq!(ions.len(), 8 * 2 * 2);
    }

    #[test]
    fn fragment_charges_follow_the_precursor() {
        assert_eq!(fragment_charges(None), 1..=1);
        assert_eq!(fragment_charges(Some(1)), 1..=1);
        assert_eq!(fragment_charges(Some(2)), 1..=1);
        assert_eq!(fragment_charges(Some(3)), 1..=2);
        assert_eq!(fragment_charges(Some(5)), 1..=2);
    }

    #[test]
    fn unknown_residues_have_no_fragments() {
        assert!(fragment_ions(&PeptideForm::unmodified("PEPXIDE"), Some(2)).is_none());
    }

    #[test]
    fn labels_show_series_ordinal_and_charge() {
        let ion = |series, ordinal, charge| FragmentIon {
            series,
            ordinal,
            charge,
            mz: 0.0,
        };
        assert_eq!(ion(IonSeries::B, 3, 1).label(), "b3");
        assert_eq!(ion(IonSeries::Y, 7, 2).label(), "y7++");
    }

    #[test]
    fn annotation_takes_the_most_intense_peak_in_the_window() {
        let ions = fragment_ions(&PeptideForm::unmodified("PEPTIDE"), Some(2)).unwrap();
        let b2 = ion(&ions, IonSeries::B, 2, 1);
        let y1 = ion(&ions, IonSeries::Y, 1, 1);
        // Two peaks near b2 (+2 ppm weak, -5 ppm strong), one 30 ppm off y1,
        // and noise.
        let mz = [
            y1 * (1.0 + 30e-6),
            b2 * (1.0 - 5e-6),
            b2 * (1.0 + 2e-6),
            500.0,
        ];
        let intensity = [100.0, 80.0, 10.0, 50.0];
        let matches = annotate(&ions, &mz, &intensity, Tolerance::Ppm(20.0));
        assert_eq!(matches.len(), 1, "y1 is outside 20 ppm");
        assert_eq!(matches[0].ion.ordinal, 2);
        assert_eq!(matches[0].peak, 1);
        assert!((matches[0].error_ppm + 5.0).abs() < 1e-6);
        // 0.5 Da reaches y1 too.
        assert_eq!(
            annotate(&ions, &mz, &intensity, Tolerance::Da(0.5)).len(),
            2
        );
    }
}
