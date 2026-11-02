//! Quantitation from isobaric reporter ions (TMT, TMTpro, iTRAQ).
//!
//! One run (or the fractions of one plex) holds every sample at once, each
//! labelled with its own reagent; a peptide's MS/MS spectrum shows one
//! reporter ion per channel, in proportion to the peptide in that sample.
//!
//! - **Reporter intensities** are read from each passing PSM's peaks: the
//!   most intense peak within the tolerance of the channel's m/z
//!   ([`extract`]). A PSM counts only when every assigned channel has a
//!   reporter; incomplete ones are counted and left out, since summing a
//!   missing channel as 0 would read as absence.
//! - **Proteins** sum their unique and razor PSMs' reporters per channel, as
//!   Proteome Discoverer does: summing weights the intense PSMs, whose
//!   ratios are the least noisy.
//! - **Normalization** assumes most proteins don't change: each channel is
//!   shifted by the median, over proteins, of its log2 value minus the
//!   protein's mean over channels. That corrects unequal loading.
//! - The result is an [`IntensityQuantities`] with one "sample" per assigned
//!   channel, so comparisons, ANOVA, and the heat map work as for MS1
//!   intensities. Values are relative within the plex.
//!
//! - **Isotope impurities** (when the plex has its lot's certificate,
//!   [`tuffold_core::reporter::Impurities`]): each PSM's reporters over the
//!   reagent set's channels are unmixed by non-negative least squares
//!   ([`correct`]) before anything else, as OpenMS and Proteome Discoverer do.
//!
//! Not corrected (yet): co-isolated peptides, which compress ratios toward
//! 1:1.

mod bridge;

pub use bridge::{Bridged, PlexQuantities, bridge};

use crate::intensity::{IntensityQuantities, IntensityRow};
use crate::{QuantError, Sample};

/// The most intense peak within `tolerance_ppm` of each channel's m/z, or
/// `None` (peaks sorted by m/z).
pub fn extract(
    mz: &[f64],
    intensity: &[f32],
    channels: &[f64],
    tolerance_ppm: f64,
) -> Vec<Option<f64>> {
    channels
        .iter()
        .map(|&target| {
            let window = target * tolerance_ppm * 1e-6;
            let start = mz.partition_point(|&m| m < target - window);
            mz[start..]
                .iter()
                .zip(&intensity[start..])
                .take_while(|(m, _)| **m <= target + window)
                .map(|(_, i)| f64::from(*i))
                .filter(|&i| i > 0.0)
                .reduce(f64::max)
        })
        .collect()
}

/// Every peak within `tolerance_ppm` of each channel's m/z, summed, or `None`
/// (peaks sorted by m/z). For iTRAQ: an Orbitrap can resolve a channel from
/// its neighbour's isotope 55 ppm away, and the certificate's correction
/// counts that isotope as read at this channel, so it must be.
pub fn extract_summed(
    mz: &[f64],
    intensity: &[f32],
    channels: &[f64],
    tolerance_ppm: f64,
) -> Vec<Option<f64>> {
    channels
        .iter()
        .map(|&target| {
            let window = target * tolerance_ppm * 1e-6;
            let start = mz.partition_point(|&m| m < target - window);
            mz[start..]
                .iter()
                .zip(&intensity[start..])
                .take_while(|(m, _)| **m <= target + window)
                .map(|(_, i)| f64::from(*i))
                .filter(|&i| i > 0.0)
                .reduce(|a, b| a + b)
        })
        .collect()
}

/// Sweeps of [`correct`] before it gives up converging (it needs a handful
/// on a certificate's nearly diagonal matrix).
const MAX_SWEEPS: usize = 1_000;

/// The intensities, none negative, whose reporters spread by `matrix`
/// (`matrix[observed][true]`, from [`Impurities::matrix`]) come closest to
/// `observed` in least squares. Solving the linear system exactly can give a
/// channel a negative intensity (a weak channel next to a strong one, with
/// noise); those are 0 here and the rest refitted, rather than clipped
/// afterwards. Coordinate descent on the normal equations, from `observed`,
/// until no channel moves by more than 1e-12 of the strongest.
///
/// [`Impurities::matrix`]: tuffold_core::reporter::Impurities::matrix
pub fn correct(observed: &[f64], matrix: &[Vec<f64>]) -> Vec<f64> {
    let n = observed.len();
    let gram: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|k| (0..n).map(|r| matrix[r][i] * matrix[r][k]).sum())
                .collect()
        })
        .collect();
    let projected: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|r| matrix[r][i] * observed[r]).sum())
        .collect();
    let scale = observed.iter().copied().fold(0.0, f64::max);
    let mut value: Vec<f64> = observed.iter().map(|&o| o.max(0.0)).collect();
    for _ in 0..MAX_SWEEPS {
        let mut moved: f64 = 0.0;
        for i in (0..n).filter(|&i| gram[i][i] > 0.0) {
            let others: f64 = (0..n)
                .filter(|&k| k != i)
                .map(|k| gram[i][k] * value[k])
                .sum();
            let next = ((projected[i] - others) / gram[i][i]).max(0.0);
            moved = moved.max((next - value[i]).abs());
            value[i] = next;
        }
        if moved <= 1e-12 * scale {
            break;
        }
    }
    value
}

/// A protein group's passing PSMs, each with its reporter intensities
/// (parallel to the channels).
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinReporters {
    pub group_id: i64,
    pub accession: String,
    pub psms: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReporterQuantities {
    /// Samples by channel, comparisons and all (see the module notes).
    pub quantities: IntensityQuantities,
    /// The channel each sample was read from.
    pub channels: Vec<String>,
    /// PSMs with every channel, which the proteins are built from.
    pub psms_used: usize,
    /// PSMs missing a reporter in some channel, left out.
    pub psms_incomplete: usize,
}

/// Sum each protein's complete PSMs per channel and normalize the channels.
/// `channels` names each sample's channel (parallel to `samples`).
pub fn quantify_reporters(
    samples: Vec<Sample>,
    channels: Vec<String>,
    proteins: Vec<ProteinReporters>,
) -> Result<ReporterQuantities, QuantError> {
    let n = samples.len();
    if let Some(bad) = proteins
        .iter()
        .find(|p| p.psms.iter().any(|psm| psm.len() != n))
    {
        return Err(QuantError::Misaligned {
            accession: bad.accession.clone(),
            found: bad.psms.iter().map(Vec::len).find(|&l| l != n).unwrap_or(0),
            expected: n,
        });
    }
    let (mut used, mut incomplete) = (0, 0);
    let summed: Vec<(i64, String, Vec<f64>, usize)> = proteins
        .into_iter()
        .filter_map(|protein| {
            let mut sum = vec![0.0; n];
            let mut count = 0;
            for psm in &protein.psms {
                if psm.iter().any(Option::is_none) {
                    incomplete += 1;
                    continue;
                }
                for (total, value) in sum.iter_mut().zip(psm) {
                    *total += value.unwrap_or_default();
                }
                count += 1;
            }
            used += count;
            (count > 0).then_some((protein.group_id, protein.accession, sum, count))
        })
        .collect();
    let shifts = channel_shifts(
        &summed
            .iter()
            .map(|(_, _, s, _)| s.clone())
            .collect::<Vec<_>>(),
        n,
    );
    let rows = summed
        .into_iter()
        .map(|(group_id, accession, sum, count)| {
            let log2: Vec<Option<f64>> = sum
                .iter()
                .zip(&shifts)
                .map(|(&v, shift)| (v > 0.0).then(|| v.log2() - shift))
                .collect();
            IntensityRow {
                group_id,
                accession,
                intensities: log2.iter().map(|l| l.map(f64::exp2)).collect(),
                log2,
                peptides: vec![count; n],
                top3: vec![None; n],
            }
        })
        .collect();
    Ok(ReporterQuantities {
        quantities: IntensityQuantities {
            samples,
            normalized: vec![!shifts.is_empty(); n],
            shifts,
            shared_peptides: used,
            rows,
        },
        channels,
        psms_used: used,
        psms_incomplete: incomplete,
    })
}

/// Per channel, the median over proteins of its log2 value minus the
/// protein's mean (proteins with every channel above 0), centred on 0.
fn channel_shifts(sums: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut deviations: Vec<Vec<f64>> = vec![Vec::new(); n];
    for sum in sums.iter().filter(|s| s.iter().all(|&v| v > 0.0)) {
        let logs: Vec<f64> = sum.iter().map(|v| v.log2()).collect();
        let mean = logs.iter().sum::<f64>() / n as f64;
        for (d, l) in deviations.iter_mut().zip(&logs) {
            d.push(l - mean);
        }
    }
    if deviations.first().is_none_or(Vec::is_empty) {
        return vec![0.0; n];
    }
    deviations
        .into_iter()
        .map(|mut d| {
            d.sort_by(f64::total_cmp);
            let mid = d.len() / 2;
            if d.len() % 2 == 1 {
                d[mid]
            } else {
                (d[mid - 1] + d[mid]) / 2.0
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str, category: &str) -> Sample {
        Sample {
            name: name.into(),
            category: Some(category.into()),
        }
    }

    #[test]
    fn reporters_are_the_strongest_peak_within_tolerance() {
        let mz = [126.1275, 126.1278, 127.1248, 128.20];
        let intensity = [10.0, 30.0, 5.0, 99.0];
        // 20 ppm of 126.127726 is 0.0025: both 126 peaks qualify; 128.134436
        // has none within 0.0026.
        let summed = extract_summed(&mz, &intensity, &[126.127726, 128.134436], 20.0);
        assert_eq!(summed[0], Some(f64::from(intensity[0] + intensity[1])));
        assert_eq!(summed[1], None);
        let found = extract(&mz, &intensity, &[126.127726, 127.124761, 128.134436], 20.0);
        assert_eq!(found, [Some(30.0), Some(5.0), None]);
    }

    #[test]
    fn correction_unmixes_and_keeps_intensities_non_negative() {
        // Three channels; each spills into its neighbours.
        let matrix = vec![
            vec![0.90, 0.02, 0.00],
            vec![0.05, 0.93, 0.03],
            vec![0.01, 0.04, 0.95],
        ];
        let truth = [1000.0, 50.0, 400.0];
        let observed: Vec<f64> = matrix
            .iter()
            .map(|row| row.iter().zip(truth).map(|(m, t)| m * t).sum())
            .collect();
        let corrected = correct(&observed, &matrix);
        assert!(
            corrected
                .iter()
                .zip(truth)
                .all(|(c, t)| (c - t).abs() < 1e-6),
            "{corrected:?}"
        );

        // Channel 0 spills 5% into channel 1, which holds only 2: solving
        // exactly gives 1 a negative intensity (2 − 0.05 · 100/0.95 = −3.26).
        // With 1 at 0, channel 0 minimizes (0.95 t − 100)² + (0.05 t − 2)²:
        // t = (0.95 · 100 + 0.05 · 2)/(0.95² + 0.05²) = 95.1/0.905.
        let matrix = vec![vec![0.95, 0.0], vec![0.05, 1.0]];
        let corrected = correct(&[100.0, 2.0], &matrix);
        assert!((corrected[0] - 95.1 / 0.905).abs() < 1e-9, "{corrected:?}");
        assert_eq!(corrected[1], 0.0);
    }

    #[test]
    fn proteins_sum_complete_psms_and_channels_are_loading_corrected() {
        let samples = vec![sample("a1", "a"), sample("a2", "a"), sample("b1", "b")];
        let channels = vec!["126".into(), "127".into(), "128".into()];
        // Channel 128 was loaded at twice the amount: every unchanged protein
        // reads 2× there. Protein Y is truly up 4× in b1 (8× raw).
        let protein = |id: i64, acc: &str, psms: Vec<Vec<Option<f64>>>| ProteinReporters {
            group_id: id,
            accession: acc.into(),
            psms,
        };
        let proteins = vec![
            protein(1, "X1", vec![vec![Some(100.0), Some(100.0), Some(200.0)]]),
            protein(
                2,
                "X2",
                vec![
                    vec![Some(50.0), Some(50.0), Some(100.0)],
                    vec![Some(10.0), None, Some(1.0)],
                ],
            ),
            protein(3, "X3", vec![vec![Some(7.0), Some(7.0), Some(14.0)]]),
            protein(4, "Y", vec![vec![Some(10.0), Some(10.0), Some(80.0)]]),
        ];
        let result = quantify_reporters(samples, channels, proteins).unwrap();
        assert_eq!((result.psms_used, result.psms_incomplete), (4, 1));
        let y = result
            .quantities
            .rows
            .iter()
            .find(|r| r.accession == "Y")
            .unwrap();
        let fold = y.log2[2].unwrap() - y.log2[0].unwrap();
        assert!(
            (fold - 2.0).abs() < 1e-9,
            "Y up 4× after correction, got 2^{fold}"
        );
        let x = result
            .quantities
            .rows
            .iter()
            .find(|r| r.accession == "X2")
            .unwrap();
        assert!((x.log2[2].unwrap() - x.log2[0].unwrap()).abs() < 1e-9);
    }
}
