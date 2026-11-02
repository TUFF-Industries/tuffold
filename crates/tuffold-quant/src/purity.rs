//! Precursor purity: how much of what an MS/MS scan isolated was its
//! precursor.
//!
//! An isolation window (often 0.7–2 Th) lets in whatever else elutes at that
//! m/z, and with isobaric labels every co-isolated peptide adds its own
//! reporters, pulling the measured ratios toward 1:1 (ratio compression).
//! Purity is the precursor's isotope peaks over all the signal inside the
//! window in an MS1 scan: the precursor ion fraction (MaxQuant's PIF,
//! OpenMS IsobaricAnalyzer's purity, the complement of Proteome
//! Discoverer's isolation interference).
//!
//! The precursor elutes while it is measured, so the MS1 scans before and
//! after an MS/MS scan see it at different heights; purity is computed in
//! both and interpolated to the MS/MS scan's time ([`PurityTracker`]), as
//! IsobaricAnalyzer does. With only one side (the run's first or last
//! scans), that side's value is used.

use std::collections::VecDeque;

use crate::ms1::ISOTOPE_SPACING;

/// m/z tolerance for each isotope peak, in parts per million. Wider than
/// MS1 tracing's 10: a Velos run's precursors sat at a median +9.7 ppm, so
/// at 10 ppm a precursor plainly in the scan read as absent (purity 0) for
/// 3.4% of passing spectra; at Orbitrap resolution an unrelated peak within
/// 20 ppm of an isotope is rare.
pub const PURITY_PPM: f64 = 20.0;
/// Isotope peaks of the precursor counted (M to M+5), where inside the
/// window.
const ISOTOPES: usize = 6;

/// The precursor ion fraction of the window `low..=high` in one MS1 scan
/// (peaks sorted by m/z, as [`PurityTracker::add_ms1`] keeps them): the
/// precursor's isotope peaks, each the most intense peak in the window within
/// `ppm` of M + k·1.00336/z, over all peak intensity in the window. `None`
/// when the window holds no signal at all, or the arrays don't pair up.
pub fn purity(
    mz: &[f64],
    intensity: &[f32],
    precursor_mz: f64,
    charge: u8,
    (low, high): (f64, f64),
    ppm: f64,
) -> Option<f64> {
    if mz.len() != intensity.len() {
        return None;
    }
    let start = mz.partition_point(|&m| m < low);
    let end = mz.partition_point(|&m| m <= high);
    let total: f64 = intensity[start..end].iter().map(|&i| f64::from(i)).sum();
    if total <= 0.0 || charge == 0 {
        return None;
    }
    let spacing = ISOTOPE_SPACING / f64::from(charge);
    let precursor: f64 = (0..ISOTOPES)
        .map(|k| precursor_mz + k as f64 * spacing)
        .take_while(|&target| target <= high)
        .filter(|&target| target >= low)
        .filter_map(|target| {
            let window = target * ppm * 1e-6;
            // Only peaks inside the window, which the total counts.
            let from = mz.partition_point(|&m| m < target - window).max(start);
            let to = mz.partition_point(|&m| m <= target + window).min(end);
            intensity[from..to]
                .iter()
                .map(|&i| f64::from(i))
                .reduce(f64::max)
        })
        .sum();
    let fraction = precursor / total;
    fraction.is_finite().then_some(fraction.min(1.0))
}

/// An MS/MS scan waiting for the MS1 scan after it.
#[derive(Debug, Clone)]
struct Pending<K> {
    key: K,
    rt: f64,
    precursor_mz: f64,
    charge: u8,
    window: (f64, f64),
    /// Purity in the MS1 scan before, and that scan's time.
    before: Option<(f64, f64)>,
}

/// Streams a run's scans in file order and gives each registered MS/MS scan
/// its purity, interpolated between the MS1 scans around it.
#[derive(Debug)]
pub struct PurityTracker<K> {
    ppm: f64,
    /// The last MS1 scan: time, m/z, intensity.
    last: Option<(f64, Vec<f64>, Vec<f32>)>,
    pending: VecDeque<Pending<K>>,
    done: Vec<(K, Option<f64>)>,
}

impl<K> PurityTracker<K> {
    pub fn new(ppm: f64) -> Self {
        Self {
            ppm,
            last: None,
            pending: VecDeque::new(),
            done: Vec::new(),
        }
    }

    /// An MS1 scan: settles every MS/MS scan waiting since the last one.
    /// Peaks with a non-finite or negative value are left out and the rest
    /// sorted by m/z, so a malformed scan can't make a purity up.
    pub fn add_ms1(&mut self, rt: f64, mz: &[f64], intensity: &[f32]) {
        if !rt.is_finite() {
            return;
        }
        let mut peaks: Vec<(f64, f32)> = mz
            .iter()
            .zip(intensity)
            .map(|(&m, &i)| (m, i))
            .filter(|&(m, i)| m.is_finite() && i.is_finite() && i >= 0.0)
            .collect();
        if !peaks.is_sorted_by(|a, b| a.0 <= b.0) {
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
        let (mz, intensity): (Vec<f64>, Vec<f32>) = peaks.into_iter().unzip();
        let (mz, intensity) = (mz.as_slice(), intensity.as_slice());
        while let Some(p) = self.pending.pop_front() {
            let after = purity(mz, intensity, p.precursor_mz, p.charge, p.window, self.ppm);
            let value = match (p.before, after) {
                (Some((before, at)), Some(after)) if rt > at => {
                    before + (after - before) * ((p.rt - at) / (rt - at)).clamp(0.0, 1.0)
                }
                (Some((before, _)), _) => before,
                (None, after) => after.unwrap_or(f64::NAN),
            };
            self.done.push((p.key, value.is_finite().then_some(value)));
        }
        self.last = Some((rt, mz.to_vec(), intensity.to_vec()));
    }

    /// An identified MS/MS scan at `rt` that isolated `window` (absolute m/z
    /// bounds) around a precursor at `precursor_mz` and `charge`.
    pub fn add_msms(&mut self, key: K, rt: f64, precursor_mz: f64, charge: u8, window: (f64, f64)) {
        let before = self.last.as_ref().and_then(|(at, mz, intensity)| {
            purity(mz, intensity, precursor_mz, charge, window, self.ppm).map(|p| (p, *at))
        });
        self.pending.push_back(Pending {
            key,
            rt,
            precursor_mz,
            charge,
            window,
            before,
        });
    }

    /// Purities settled so far (taken out).
    pub fn take(&mut self) -> Vec<(K, Option<f64>)> {
        std::mem::take(&mut self.done)
    }

    /// The end of the run: scans after the last MS1 scan keep the purity
    /// before them.
    pub fn finish(mut self) -> Vec<(K, Option<f64>)> {
        let rest = self
            .pending
            .drain(..)
            .map(|p| (p.key, p.before.map(|(v, _)| v)));
        self.done.extend(rest);
        self.done
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn purity_is_the_precursor_isotopes_share_of_the_window() {
        // A 2+ precursor at 500.0: isotopes at 500.0, 500.501678, 501.003355.
        // Window 499.5..501.5 holds those three (100, 50, 20), a contaminant
        // at 500.9 (30) and noise at 499.6 (10); 502.0 is outside.
        let mz = [499.6, 500.0, 500.501_678, 500.9, 501.003_355, 502.0];
        let intensity = [10.0, 100.0, 50.0, 30.0, 20.0, 999.0];
        let p = purity(&mz, &intensity, 500.0, 2, (499.5, 501.5), PURITY_PPM).unwrap();
        // (100 + 50 + 20) / (10 + 100 + 50 + 30 + 20) = 170 / 210.
        assert!((p - 170.0 / 210.0).abs() < 1e-12, "{p}");
        // Nothing in the window: unknown, not 0.
        assert_eq!(
            purity(&mz, &intensity, 600.0, 2, (599.5, 600.5), 10.0),
            None
        );
        // No peak at all between 499.7 and 499.9: unknown again.
        assert_eq!(
            purity(&mz, &intensity, 499.8, 2, (499.7, 499.9), 10.0),
            None
        );
        // The precursor absent from a window with signal: 0.
        assert_eq!(
            purity(&mz, &intensity, 499.55, 1, (499.5, 499.65), 10.0),
            Some(0.0)
        );
    }

    #[test]
    fn a_malformed_scan_cannot_make_a_purity_up() {
        // Unsorted, with a NaN intensity and an infinite m/z: the NaN and
        // infinite peaks are dropped and the rest sorted, leaving 500.0 (100)
        // and 499.8 (100) in the window: 0.5.
        let mut tracker = PurityTracker::new(PURITY_PPM);
        tracker.add_msms("a", 11.0, 500.0, 1, (499.5, 500.5));
        tracker.add_ms1(
            12.0,
            &[500.0, f64::INFINITY, 499.8, 500.2],
            &[100.0, 50.0, 100.0, f32::NAN],
        );
        assert_eq!(tracker.take(), [("a", Some(0.5))]);
        // Arrays that don't pair up give no purity rather than a panic.
        assert_eq!(
            purity(&[500.0, 501.0], &[1.0], 500.0, 1, (499.0, 502.0), 10.0),
            None
        );
    }

    #[test]
    fn purity_is_interpolated_between_the_ms1_scans_around() {
        let window = (499.5, 500.5);
        // Before (t = 10): precursor 100 of 200 (0.5); after (t = 14): 300 of
        // 400 (0.75). An MS/MS scan at t = 11 is a quarter of the way: 0.5625.
        let mut tracker = PurityTracker::new(PURITY_PPM);
        tracker.add_ms1(10.0, &[499.8, 500.0], &[100.0, 100.0]);
        tracker.add_msms("a", 11.0, 500.0, 1, window);
        assert!(tracker.take().is_empty(), "waits for the next MS1 scan");
        tracker.add_ms1(14.0, &[499.8, 500.0], &[100.0, 300.0]);
        let done = tracker.take();
        assert_eq!(done.len(), 1);
        assert!((done[0].1.unwrap() - 0.5625).abs() < 1e-12, "{done:?}");
        // One after the last MS1 scan keeps the one before (0.75); one before
        // any MS1 scan takes the one after.
        tracker.add_msms("b", 15.0, 500.0, 1, window);
        let mut first = PurityTracker::new(PURITY_PPM);
        first.add_msms("c", 1.0, 500.0, 1, window);
        first.add_ms1(2.0, &[500.0], &[50.0]);
        assert_eq!(first.take(), [("c", Some(1.0))]);
        assert_eq!(tracker.finish(), [("b", Some(0.75))]);
    }
}
