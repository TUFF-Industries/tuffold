//! The MS/MS peaks a Scaffold project stores (`docs/sf3-format.md`).
//!
//! A spectrum (`Spectrum` 4) names its `PeakList` record (180 bytes): the
//! record's id, 20 (intensity, m/z) pairs, and at 164 the head of a chain of
//! `PeakListMassIntensities` blocks (232 bytes: the chain's id, 28 more
//! pairs, the next block's index in its part) holding the rest, NaN after
//! the last pair; 0 means no more. The two never share a peak: the 20 are
//! Scaffold's picks (no isotope peaks among them), the chain everything else,
//! by falling intensity, so a spectrum's peaks are both together. Scaffold
//! keeps at most 300. Unused slots of the 20 hold (0, 0), or (−1, −1) in
//! 5.2.2. All little-endian `f32`.
//!
//! Checked on 19 sample files (Scaffold 4.0–5.3.3): no record's 20 pairs
//! reappear in its chain, and fragment ions of Scaffold's confident PSMs
//! match the two together better than either alone (e.g. 59% of b/y ions
//! within 0.02 Da against 39% for the 20 alone).

use std::collections::HashSet;
use std::sync::Arc;

use tuffold_core::{PeakSpectrum, SpectrumKey};

use super::Sf3Error;
use super::records::{Container, Records, read_f32, read_i32};

pub(crate) const SUMMARY_BLOCK: usize = 180;
const PEAK_BLOCK: usize = 232;
const SUMMARY_PAIRS: usize = 20;
const CHAIN_HEAD: usize = 164;
/// Unused slots of the 20: (0, 0), or (−1, −1) in 5.2.2.
const EMPTY_SLOTS: [(f32, f32); 2] = [(0.0, 0.0), (-1.0, -1.0)];

/// Where one identified spectrum's peaks are, and what it is.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PeakRef {
    pub run: String,
    pub key: SpectrumKey,
    pub native_id: String,
    pub peak_list: i32,
    pub rt_seconds: Option<f64>,
    pub precursor_mz: Option<f64>,
    pub precursor_charge: Option<i16>,
}

/// The peaks of a project's identified spectra, read one spectrum at a
/// time (a large project stores hundreds of megabytes of them).
#[derive(Clone)]
pub struct Sf3Peaks {
    inner: Arc<Stored>,
}

struct Stored {
    summaries: Records,
    chains: Records,
    spectra: Vec<PeakRef>,
}

impl std::fmt::Debug for Sf3Peaks {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sf3Peaks({} spectra)", self.inner.spectra.len())
    }
}

/// Two handles are equal when they hold the same stored peaks.
impl PartialEq for Sf3Peaks {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
    }
}

/// Why one spectrum's peaks couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("{native_id}: {problem}")]
pub struct UnreadablePeaks {
    pub native_id: String,
    pub problem: &'static str,
}

impl Sf3Peaks {
    /// The runs holding spectra with peaks, in first-seen order.
    pub fn runs(&self) -> Vec<&str> {
        let mut seen = HashSet::new();
        self.inner
            .spectra
            .iter()
            .map(|spectrum| spectrum.run.as_str())
            .filter(|run| seen.insert(*run))
            .collect()
    }

    /// Identified spectra (any run) that name a peak list.
    pub fn len(&self) -> usize {
        self.inner.spectra.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.spectra.is_empty()
    }

    /// The run's spectra with their peaks, sorted by m/z. A spectrum whose
    /// stored list is empty is left out: it has no peaks to show.
    pub fn spectra<'a>(
        &'a self,
        run: &'a str,
    ) -> impl Iterator<Item = Result<PeakSpectrum, UnreadablePeaks>> + 'a {
        self.inner
            .spectra
            .iter()
            .filter(move |spectrum| spectrum.run == run)
            .map(|spectrum| self.inner.read(spectrum))
            .filter(|read| !matches!(read, Ok(spectrum) if spectrum.mz.is_empty()))
    }
}

impl Stored {
    fn read(&self, spectrum: &PeakRef) -> Result<PeakSpectrum, UnreadablePeaks> {
        let unreadable = |problem| UnreadablePeaks {
            native_id: spectrum.native_id.clone(),
            problem,
        };
        let summary = self
            .summaries
            .get(spectrum.peak_list)
            .ok_or_else(|| unreadable("its PeakList record is missing"))?;
        let mut pairs: Vec<(f32, f32)> = (0..SUMMARY_PAIRS)
            .filter_map(|k| pair(summary, 4 + 8 * k))
            .filter(|&pair| !EMPTY_SLOTS.contains(&pair))
            .collect();
        let head = read_i32(summary, CHAIN_HEAD).unwrap_or(0);
        if head != 0 {
            let blocks = self
                .chains
                .chain(head)
                .ok_or_else(|| unreadable("its chain of further peaks is broken"))?;
            'blocks: for block in blocks {
                for k in 0..(PEAK_BLOCK - 8) / 8 {
                    match pair(block, 4 + 8 * k) {
                        Some(found) => pairs.push(found),
                        None => break 'blocks, // NaN ends the list
                    }
                }
            }
        }
        if pairs
            .iter()
            .any(|&(intensity, mz)| !(mz.is_finite() && mz > 0.0 && intensity.is_finite()))
        {
            return Err(unreadable("a peak has an invalid m/z or intensity"));
        }
        pairs.sort_by(|a, b| a.1.total_cmp(&b.1));
        Ok(PeakSpectrum {
            key: spectrum.key.clone(),
            native_id: spectrum.native_id.clone(),
            ms_level: 2,
            rt_seconds: spectrum.rt_seconds,
            precursor_mz: spectrum.precursor_mz,
            precursor_charge: spectrum.precursor_charge,
            mz: pairs.iter().map(|&(_, mz)| f64::from(mz)).collect(),
            intensity: pairs.iter().map(|&(intensity, _)| intensity).collect(),
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        })
    }
}

/// An (intensity, m/z) pair at `offset`; `None` at a NaN, which ends a list.
fn pair(block: &[u8], offset: usize) -> Option<(f32, f32)> {
    let intensity = read_f32(block, offset)?;
    let mz = read_f32(block, offset + 4)?;
    (!intensity.is_nan() && !mz.is_nan()).then_some((intensity, mz))
}

/// The peaks of `spectra`, or `None` when the project stores none.
pub(crate) fn load(
    container: &Container,
    spectra: Vec<PeakRef>,
) -> Result<Option<Sf3Peaks>, Sf3Error> {
    if container.optional_stream("PeakList")?.is_none() {
        return Ok(None);
    }
    let summaries = container.records("PeakList", SUMMARY_BLOCK)?;
    let chains = container.chains("PeakListMassIntensities", PEAK_BLOCK)?;
    Ok(Some(Sf3Peaks {
        inner: Arc::new(Stored {
            summaries,
            chains,
            spectra,
        }),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sf3::records::test_records;

    fn summary(id: i32, pairs: &[(f32, f32)], head: i32) -> Vec<u8> {
        let mut block = id.to_le_bytes().to_vec();
        for k in 0..SUMMARY_PAIRS {
            let (intensity, mz) = pairs.get(k).copied().unwrap_or((0.0, 0.0));
            block.extend_from_slice(&intensity.to_le_bytes());
            block.extend_from_slice(&mz.to_le_bytes());
        }
        block.extend_from_slice(&head.to_le_bytes());
        block.resize(SUMMARY_BLOCK, 0);
        block
    }

    fn chain_block(id: i32, pairs: &[(f32, f32)], next: i32) -> Vec<u8> {
        let mut block = id.to_le_bytes().to_vec();
        for k in 0..28 {
            let (intensity, mz) = pairs.get(k).copied().unwrap_or((f32::NAN, f32::NAN));
            block.extend_from_slice(&intensity.to_le_bytes());
            block.extend_from_slice(&mz.to_le_bytes());
        }
        block.extend_from_slice(&next.to_le_bytes());
        block
    }

    fn reference(peak_list: i32) -> PeakRef {
        PeakRef {
            run: "r1".into(),
            key: SpectrumKey::Scan(u32::try_from(peak_list).unwrap()),
            native_id: format!("spectrum {peak_list}"),
            peak_list,
            rt_seconds: Some(60.0),
            precursor_mz: Some(500.0),
            precursor_charge: Some(2),
        }
    }

    #[test]
    fn a_spectrum_is_its_twenty_peaks_and_its_chain_sorted_by_mz() {
        // Record 1: two picked peaks, then a chain of two blocks (1 → 2)
        // with 28 + 1 more; record 2: picked peaks only (head 0).
        let full: Vec<(f32, f32)> = (0..28).map(|k| (10.0, 300.0 + k as f32)).collect();
        let stored = Stored {
            summaries: test_records(
                SUMMARY_BLOCK,
                &[
                    vec![0; SUMMARY_BLOCK],
                    summary(1, &[(900.0, 750.25), (500.0, 120.5)], 1),
                    summary(2, &[(7.0, 410.0), (-1.0, -1.0), (-1.0, -1.0)], 0),
                ],
            ),
            chains: test_records(
                PEAK_BLOCK,
                &[
                    vec![0; PEAK_BLOCK],
                    chain_block(1, &full, 2),
                    chain_block(1, &[(3.0, 99.0)], -1),
                ],
            ),
            spectra: vec![reference(1), reference(2)],
        };
        let one = stored.read(&stored.spectra[0]).unwrap();
        assert_eq!(one.mz.len(), 2 + 28 + 1);
        assert_eq!(one.mz[..3], [99.0, 120.5, 300.0]);
        assert_eq!(one.intensity[..3], [3.0, 500.0, 10.0]);
        assert_eq!(one.mz.last(), Some(&750.25));
        assert!(one.mz.windows(2).all(|w| w[0] <= w[1]));
        let two = stored.read(&stored.spectra[1]).unwrap();
        assert_eq!((two.mz, two.intensity), (vec![410.0], vec![7.0]));
    }

    #[test]
    fn a_broken_chain_or_bad_peak_is_reported_not_guessed() {
        let stored = Stored {
            summaries: test_records(
                SUMMARY_BLOCK,
                &[
                    vec![0; SUMMARY_BLOCK],
                    summary(1, &[(9.0, 100.0)], 5),
                    summary(2, &[(9.0, -3.0)], 0),
                ],
            ),
            chains: test_records(PEAK_BLOCK, &[vec![0; PEAK_BLOCK]]),
            spectra: vec![reference(1), reference(2), reference(3)],
        };
        let problem = |index: usize| stored.read(&stored.spectra[index]).unwrap_err().problem;
        assert_eq!(problem(0), "its chain of further peaks is broken");
        assert_eq!(problem(1), "a peak has an invalid m/z or intensity");
        assert_eq!(problem(2), "its PeakList record is missing");
    }
}
