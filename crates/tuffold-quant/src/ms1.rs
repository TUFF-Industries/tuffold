//! Precursor intensity from MS1 scans. For each identified precursor (a
//! peptide at one charge in one run), the tracer follows its monoisotopic
//! and M+1 peaks across the MS1 scans around its identifications and
//! integrates the elution peak. M+2 is left out on purpose: it falls under
//! the noise floor first at low abundance, so including it would shrink a
//! weak sample's area by more than a strong one's and compress fold changes.
//!
//! It streams. MS1 scans and identified MS/MS scans arrive in file order, and
//! only the last `margin` seconds of MS1 scans are kept, so a precursor
//! identified halfway down its peak still sees where the peak began. Nothing
//! needs to be known about a precursor before its first identification
//! passes by.
//!
//! The peak: start from the scan nearest the identifications (their median
//! retention time) and climb the 3-point smoothed trace to its apex. Then walk
//! out each side until the signal falls below 5% of the apex, more than
//! `MAX_GAP` scans in a row miss the monoisotopic peak, or the signal rises
//! again after falling below a fifth of the apex (a neighbouring peak). The
//! area is the trapezoid integral over retention time (seconds) of the summed
//! isotope intensities.
//! A scan counts only when the monoisotopic peak is present. Up to `MAX_GAP`
//! scans missing it inside a peak (common at the edge of the m/z tolerance,
//! and near the detection limit, where a peak's scans are detected one in
//! two) are bridged linearly between the scans that have it, which is
//! exactly what integrating across them does. They are never counted as
//! points or chosen as the apex. Without bridging, a weak peak was cut at
//! its first run of two misses and a dip to half the apex, so a precursor
//! twice as abundant in one sample had a peak several times longer there,
//! and folds between samples were exaggerated.
//!
//! A precursor whose apex lacks the monoisotopic or the M+1 peak, or whose
//! peak spans fewer than 3 scans, gets no feature (`None`). It was measured
//! and nothing was found, which is different from zero and is never reported
//! as one.
//!
//! Ion mobility (timsTOF): an MS1 frame binned by mobility holds an ion once
//! per bin it spans, each peak with its bin's 1/K0. A precursor whose
//! mobility is known is measured only within `mobility_window` of it, the
//! sum of its isotope's peaks over those bins; ions at the same m/z but
//! another mobility are left out, which collapsing the frame would have
//! added in. A precursor without a mobility gets the sum over every bin,
//! which is the collapsed frame.

use std::collections::{HashMap, VecDeque};
use std::hash::Hash;

/// ¹³C − ¹²C, the spacing of a peptide's isotope peaks (in Da).
pub const ISOTOPE_SPACING: f64 = 1.003_355;
/// Isotope peaks traced and summed: M, M+1.
pub const ISOTOPES: usize = 2;
/// A peak ends where the smoothed signal falls below this share of the apex.
const TAIL: f64 = 0.05;
/// Rising again after falling below this share of the apex means a neighbour.
/// A weak peak's smoothed trace dips below half its apex between scans of
/// the same peak, so only a deeper valley separates two.
const VALLEY: f64 = 0.2;
/// Most consecutive scans without the monoisotopic peak that are bridged when
/// the peak resumes; a longer run ends the peak.
const MAX_GAP: usize = 3;
/// Fewer scans than this with signal is not a peak.
const MIN_POINTS: usize = 3;
/// Default half-width of a precursor's mobility window (1/K0).
pub const MOBILITY_WINDOW: f64 = 0.015;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceSettings {
    /// m/z tolerance for each isotope peak, in parts per million.
    pub ppm: f64,
    /// How far (seconds) before the first and after the last identification
    /// a precursor's trace reaches.
    pub margin_seconds: f64,
    /// Half-width (1/K0, V·s/cm²) of the ion mobility window around a
    /// precursor's own mobility, in frames binned by mobility.
    pub mobility_window: f64,
}

impl Default for TraceSettings {
    /// 10 ppm suits Orbitrap and TOF MS1; 90 s covers typical peaks.
    fn default() -> Self {
        Self {
            ppm: 10.0,
            margin_seconds: 90.0,
            mobility_window: MOBILITY_WINDOW,
        }
    }
}

/// A precursor's integrated elution peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Feature {
    /// Σ isotope intensity × seconds.
    pub area: f64,
    pub apex_rt: f64,
    /// Summed isotope intensity at the apex scan.
    pub apex_intensity: f64,
    pub start_rt: f64,
    pub end_rt: f64,
    /// Scans in the peak with the monoisotopic peak present.
    pub points: usize,
    /// At the apex: the monoisotopic peak's m/z error from the traced m/z
    /// (ppm), and the M+1 / M intensity ratio. Match-between-runs judges a
    /// transferred peak by them.
    pub apex_ppm: f64,
    pub apex_ratio: f64,
}

/// One MS1 scan, peaks sorted by m/z, with each peak's ion mobility when
/// the scan is binned by mobility (else `mobility` is empty).
struct Scan {
    rt: f64,
    mz: Vec<f64>,
    intensity: Vec<f32>,
    mobility: Vec<f32>,
}

impl Scan {
    /// Peaks with a non-finite m/z or intensity are dropped. Mobilities are
    /// used only when there is one per peak.
    fn new(rt: f64, mz: &[f64], intensity: &[f32], mobility: &[f32]) -> Self {
        let count = mz.len().min(intensity.len());
        let mobility = if mobility.len() == mz.len() {
            &mobility[..count]
        } else {
            &[]
        };
        let finite = (0..count).all(|i| mz[i].is_finite() && intensity[i].is_finite());
        if finite && mz[..count].is_sorted() {
            return Self {
                rt,
                mz: mz[..count].to_vec(),
                intensity: intensity[..count].to_vec(),
                mobility: mobility.to_vec(),
            };
        }
        let mut order: Vec<usize> = (0..count)
            .filter(|&i| mz[i].is_finite() && intensity[i].is_finite())
            .collect();
        order.sort_by(|&a, &b| mz[a].total_cmp(&mz[b]));
        Self {
            rt,
            mz: order.iter().map(|&i| mz[i]).collect(),
            intensity: order.iter().map(|&i| intensity[i]).collect(),
            mobility: if mobility.is_empty() {
                Vec::new()
            } else {
                order.iter().map(|&i| mobility[i]).collect()
            },
        }
    }

    /// The signal within `ppm` of `target` and its m/z, or 0 at `target`:
    /// the most intense peak, or in a scan binned by mobility the sum of the
    /// peaks in `window` (1/K0 bounds; every bin without one) at their
    /// intensity-weighted m/z.
    fn peak_at(&self, target: f64, ppm: f64, window: Option<(f64, f64)>) -> (f32, f64) {
        let tolerance = target * ppm * 1e-6;
        let start = self.mz.partition_point(|&mz| mz < target - tolerance);
        let end = self.mz.partition_point(|&mz| mz <= target + tolerance);
        if self.mobility.is_empty() {
            return (start..end).fold((0.0, target), |best, i| {
                if self.intensity[i] > best.0 {
                    (self.intensity[i], self.mz[i])
                } else {
                    best
                }
            });
        }
        let inside = |i: &usize| {
            window.is_none_or(|(low, high)| (low..=high).contains(&f64::from(self.mobility[*i])))
        };
        let (total, weighted) = (start..end).filter(inside).fold((0.0, 0.0), |(t, w), i| {
            let intensity = f64::from(self.intensity[i]);
            (t + intensity, w + intensity * self.mz[i])
        });
        if total > 0.0 {
            #[allow(clippy::cast_possible_truncation)]
            (total as f32, weighted / total)
        } else {
            (0.0, target)
        }
    }
}

/// An open trace of one precursor.
struct Trace {
    mz: f64,
    charge: u8,
    /// The precursor's mobility window (1/K0), when its mobility is known.
    window: Option<(f64, f64)>,
    /// Retention times of its identifications.
    identified_at: Vec<f64>,
    /// Scans after this retention time are past the trace.
    until: f64,
    /// Each scan: time, isotope intensities, and the monoisotopic peak's m/z.
    points: Vec<(f64, [f32; ISOTOPES], f64)>,
}

impl Trace {
    fn sample(&mut self, scan: &Scan, ppm: f64) {
        let mut isotopes = [0.0; ISOTOPES];
        let mut mono_mz = self.mz;
        for (n, slot) in isotopes.iter_mut().enumerate() {
            let target = self.mz + n as f64 * ISOTOPE_SPACING / f64::from(self.charge);
            let (intensity, mz) = scan.peak_at(target, ppm, self.window);
            *slot = intensity;
            if n == 0 {
                mono_mz = mz;
            }
        }
        self.points.push((scan.rt, isotopes, mono_mz));
    }

    fn feature(&self) -> Option<Feature> {
        let signal: Vec<f64> = self
            .points
            .iter()
            .map(|(_, iso, _)| {
                if iso[0] > 0.0 {
                    iso.iter().map(|&i| f64::from(i)).sum()
                } else {
                    0.0
                }
            })
            .collect();
        // Dropouts of up to `MAX_GAP` scans between measured ones bridged;
        // `signal` keeps what was measured.
        let bridged = bridge(&signal);
        let smooth: Vec<f64> = (0..bridged.len())
            .map(|i| {
                let window = &bridged[i.saturating_sub(1)..(i + 2).min(bridged.len())];
                window.iter().sum::<f64>() / window.len() as f64
            })
            .collect();
        let mut apex = self.apex(&smooth)?;
        if signal[apex] == 0.0 {
            // Climbed onto a bridged scan: the strongest measured scan within
            // the bridged run's reach.
            apex = (0..signal.len())
                .filter(|&i| signal[i] > 0.0 && i.abs_diff(apex) <= MAX_GAP)
                .max_by(|&a, &b| signal[a].total_cmp(&signal[b]))?;
        }
        let (_, isotopes, apex_mz) = self.points[apex];
        if isotopes[0] <= 0.0 || isotopes[1] <= 0.0 {
            return None;
        }
        let top = smooth[apex];
        let stop = |from: usize, to: usize| {
            // Stepping from `from` to its neighbour `to`.
            let missed = |i: usize| bridged[i] == 0.0;
            let gap = missed(to)
                && if to < from {
                    to == 0 || missed(to - 1)
                } else {
                    to + 1 >= signal.len() || missed(to + 1)
                };
            gap || smooth[to] < TAIL * top
                || (smooth[to] > smooth[from] && smooth[from] < VALLEY * top)
        };
        let mut start = apex;
        while start > 0 && !stop(start, start - 1) {
            start -= 1;
        }
        let mut end = apex;
        while end + 1 < signal.len() && !stop(end, end + 1) {
            end += 1;
        }
        let points = signal[start..=end].iter().filter(|&&s| s > 0.0).count();
        if points < MIN_POINTS {
            return None;
        }
        let area = (start..end)
            .map(|i| {
                let width = self.points[i + 1].0 - self.points[i].0;
                width * (bridged[i] + bridged[i + 1]) / 2.0
            })
            .sum();
        Some(Feature {
            area,
            apex_rt: self.points[apex].0,
            apex_intensity: signal[apex],
            start_rt: self.points[start].0,
            end_rt: self.points[end].0,
            points,
            apex_ppm: (apex_mz - self.mz) / self.mz * 1e6,
            apex_ratio: f64::from(isotopes[1]) / f64::from(isotopes[0]),
        })
    }

    /// The local maximum of `smooth` reached by climbing from the scan
    /// nearest the identifications (or the nearest scan with signal, within
    /// two scans of it).
    fn apex(&self, smooth: &[f64]) -> Option<usize> {
        let center = median(&self.identified_at)?;
        let nearest = (0..self.points.len()).min_by(|&a, &b| {
            (self.points[a].0 - center)
                .abs()
                .total_cmp(&(self.points[b].0 - center).abs())
        })?;
        let mut at = [0, 1, -1, 2, -2]
            .into_iter()
            .filter_map(|offset: isize| nearest.checked_add_signed(offset))
            .find(|&i| smooth.get(i).is_some_and(|&s| s > 0.0))?;
        loop {
            let left = at.checked_sub(1).filter(|&i| smooth[i] > smooth[at]);
            let right = Some(at + 1).filter(|&i| i < smooth.len() && smooth[i] > smooth[at]);
            at = match (left, right) {
                (Some(l), Some(r)) => {
                    if smooth[l] >= smooth[r] {
                        l
                    } else {
                        r
                    }
                }
                (Some(next), None) | (None, Some(next)) => next,
                (None, None) => return Some(at),
            };
        }
    }
}

/// `signal` with each run of at most `MAX_GAP` zeros between two measured
/// scans replaced by the straight line between them.
fn bridge(signal: &[f64]) -> Vec<f64> {
    let mut bridged = signal.to_vec();
    let mut i = 1;
    while i < signal.len() {
        if signal[i] != 0.0 || signal[i - 1] == 0.0 {
            i += 1;
            continue;
        }
        let end = (i..signal.len())
            .find(|&j| signal[j] != 0.0)
            .unwrap_or(signal.len());
        if end < signal.len() && end - i <= MAX_GAP {
            let (before, after) = (signal[i - 1], signal[end]);
            let steps = (end - i + 1) as f64;
            for (k, slot) in bridged[i..end].iter_mut().enumerate() {
                *slot = before + (after - before) * (k + 1) as f64 / steps;
            }
        }
        i = end;
    }
    bridged
}

fn median(values: &[f64]) -> Option<f64> {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    match n {
        0 => None,
        _ if n % 2 == 1 => Some(sorted[n / 2]),
        _ => Some((sorted[n / 2 - 1] + sorted[n / 2]) / 2.0),
    }
}

/// Traces precursors through one run's scans, fed in file order.
pub struct Tracer<K> {
    settings: TraceSettings,
    recent: VecDeque<Scan>,
    open: HashMap<K, Trace>,
    closed: Vec<(K, Option<Feature>)>,
    ms1_scans: usize,
}

impl<K: Eq + Hash + Clone> Tracer<K> {
    pub fn new(settings: TraceSettings) -> Self {
        Self {
            settings,
            recent: VecDeque::new(),
            open: HashMap::new(),
            closed: Vec::new(),
            ms1_scans: 0,
        }
    }

    /// MS1 scans seen so far.
    pub fn ms1_scans(&self) -> usize {
        self.ms1_scans
    }

    /// The next MS1 scan (peaks centroided; any order), with each peak's ion
    /// mobility when the scan is binned by mobility (else empty).
    pub fn add_ms1(&mut self, rt: f64, mz: &[f64], intensity: &[f32], mobility: &[f32]) {
        if !rt.is_finite() {
            return;
        }
        self.ms1_scans += 1;
        let scan = Scan::new(rt, mz, intensity, mobility);
        // Only traces near this retention time are open (about ±margin of
        // their identifications), so these passes are short.
        let ended: Vec<(K, Trace)> = self.open.extract_if(|_, trace| rt > trace.until).collect();
        self.closed
            .extend(ended.into_iter().map(|(key, trace)| (key, trace.feature())));
        for trace in self.open.values_mut() {
            trace.sample(&scan, self.settings.ppm);
        }
        self.recent.push_back(scan);
        let oldest = rt - self.settings.margin_seconds;
        while self.recent.front().is_some_and(|scan| scan.rt < oldest) {
            self.recent.pop_front();
        }
    }

    /// An MS/MS scan at `rt` identified as precursor `key`, whose
    /// monoisotopic peak is at `mz` and ion mobility (1/K0) `mobility` when
    /// known. A precursor's trace keeps its first identification's mobility.
    pub fn add_identification(
        &mut self,
        key: K,
        mz: f64,
        charge: u8,
        rt: f64,
        mobility: Option<f64>,
    ) {
        if charge == 0 || !mz.is_finite() || !rt.is_finite() {
            return;
        }
        let margin = self.settings.margin_seconds;
        if let Some(trace) = self.open.get_mut(&key) {
            trace.identified_at.push(rt);
            trace.until = trace.until.max(rt + margin);
            return;
        }
        let half = self.settings.mobility_window;
        let mut trace = Trace {
            mz,
            charge,
            window: mobility
                .filter(|m| m.is_finite() && *m > 0.0)
                .map(|m| (m - half, m + half)),
            identified_at: vec![rt],
            until: rt + margin,
            points: Vec::new(),
        };
        for scan in self.recent.iter().filter(|scan| scan.rt >= rt - margin) {
            trace.sample(scan, self.settings.ppm);
        }
        self.open.insert(key, trace);
    }

    /// Every traced precursor's feature (`None`: traced, no peak). A
    /// precursor identified again after its trace closed was traced twice;
    /// the larger peak stands.
    pub fn finish(mut self) -> HashMap<K, Option<Feature>> {
        let open: Vec<(K, Trace)> = self.open.drain().collect();
        self.closed
            .extend(open.into_iter().map(|(key, trace)| (key, trace.feature())));
        let mut features: HashMap<K, Option<Feature>> = HashMap::new();
        for (key, feature) in self.closed {
            let best = features.entry(key).or_insert(None);
            if feature.map(|f| f.area) > best.map(|f| f.area) {
                *best = feature;
            }
        }
        features
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MZ: f64 = 500.0;
    const RATIOS: [f64; ISOTOPES] = [1.0, 0.6];

    /// Gaussian elution peaks `(center, height)` of a 2+ precursor at `MZ`,
    /// sampled once a second over 0..=240 s.
    fn scans(peaks: &[(f64, f64)]) -> Vec<(f64, Vec<f64>, Vec<f32>)> {
        (0..=240)
            .map(|second| {
                let rt = f64::from(second);
                let height: f64 = peaks
                    .iter()
                    .map(|&(center, h)| h * (-(rt - center).powi(2) / (2.0 * 25.0)).exp())
                    .sum();
                let mz = (0..ISOTOPES)
                    .map(|n| MZ + n as f64 * ISOTOPE_SPACING / 2.0)
                    .collect();
                let intensity = RATIOS.iter().map(|r| (r * height) as f32).collect();
                (rt, mz, intensity)
            })
            .collect()
    }

    /// Feed scans in order, with identifications of key 1 at `ids`.
    fn trace(
        scans: &[(f64, Vec<f64>, Vec<f32>)],
        ids: &[f64],
        settings: TraceSettings,
    ) -> Option<Feature> {
        let mut tracer = Tracer::new(settings);
        let mut pending: Vec<f64> = ids.to_vec();
        for (rt, mz, intensity) in scans {
            pending.retain(|&id| {
                if id <= *rt {
                    tracer.add_identification(1, MZ, 2, id, None);
                    false
                } else {
                    true
                }
            });
            tracer.add_ms1(*rt, mz, intensity, &[]);
        }
        tracer.finish()[&1]
    }

    /// Height × Σ ratios × σ √(2π), with σ = 5 s.
    fn gaussian_area(height: f64) -> f64 {
        height * RATIOS.iter().sum::<f64>() * 5.0 * (2.0 * std::f64::consts::PI).sqrt()
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        let error = (actual - expected).abs() / expected;
        assert!(error < tolerance, "{actual} vs {expected} ({error})");
    }

    #[test]
    fn integrates_a_peak_identified_on_its_tail() {
        // Identified 8 s past the apex: the scans before it come from the
        // buffer. The 5% cutoff trims the tails (~1.4% of a Gaussian).
        let feature = trace(&scans(&[(100.0, 1e6)]), &[108.0], TraceSettings::default()).unwrap();
        assert_close(feature.area, gaussian_area(1e6), 0.02);
        assert_eq!(feature.apex_rt, 100.0);
        // The isotopes sit exactly at the traced m/z, in RATIOS (M+1 = 0.6 M;
        // both rounded to f32, so compared to 1e-6).
        assert_eq!(feature.apex_ppm, 0.0);
        assert!(
            (feature.apex_ratio - 0.6).abs() < 1e-6,
            "{}",
            feature.apex_ratio
        );
        assert_close(feature.apex_intensity, 1.6e6, 1e-6);
        assert!(feature.start_rt < 90.0 && feature.end_rt > 110.0);
    }

    #[test]
    fn a_single_missed_scan_inside_a_peak_is_bridged() {
        // Scans 97 and 104 lose the monoisotopic peak (just outside the
        // tolerance), one on each side of the apex.
        let dropped: Vec<_> = scans(&[(100.0, 1e6)])
            .into_iter()
            .map(|(rt, mut mz, intensity)| {
                if rt == 97.0 || rt == 104.0 {
                    mz[0] += 0.1;
                }
                (rt, mz, intensity)
            })
            .collect();
        let feature = trace(&dropped, &[100.0], TraceSettings::default()).unwrap();
        assert!(
            feature.start_rt < 90.0 && feature.end_rt > 110.0,
            "{feature:?}"
        );
        assert_close(feature.area, gaussian_area(1e6), 0.02);
        // Bridged scans aren't measurements.
        let whole = trace(&scans(&[(100.0, 1e6)]), &[100.0], TraceSettings::default()).unwrap();
        assert_eq!(feature.points + 2, whole.points);
    }

    /// `scans` with the monoisotopic peak moved out of tolerance at `lost`.
    fn without_monoisotopic(
        scans: Vec<(f64, Vec<f64>, Vec<f32>)>,
        lost: &[f64],
    ) -> Vec<(f64, Vec<f64>, Vec<f32>)> {
        scans
            .into_iter()
            .map(|(rt, mut mz, intensity)| {
                if lost.contains(&rt) {
                    mz[0] += 0.1;
                }
                (rt, mz, intensity)
            })
            .collect()
    }

    #[test]
    fn a_run_of_up_to_three_missed_scans_is_bridged() {
        // 96..=98 lost, three in a row: the line between 95 and 99 replaces
        // them. Exact over [95, 99] is h sigma sqrt(2 pi) (Phi(-0.2) -
        // Phi(-1)) = 3.285 h, the chord 4 (0.6065 + 0.9802) / 2 h = 3.173 h:
        // 0.9% of the area, on top of the 1.4% the 5% cutoff trims.
        let lost = without_monoisotopic(scans(&[(100.0, 1e6)]), &[96.0, 97.0, 98.0]);
        let feature = trace(&lost, &[100.0], TraceSettings::default()).unwrap();
        assert_close(feature.area, gaussian_area(1e6), 0.035);
        let whole = trace(&scans(&[(100.0, 1e6)]), &[100.0], TraceSettings::default()).unwrap();
        assert_eq!(feature.points + 3, whole.points);
        assert!(feature.start_rt < 90.0, "{feature:?}");
    }

    #[test]
    fn a_run_of_four_missed_scans_ends_the_peak() {
        // 96..=99 lost: more than MAX_GAP, so the peak starts at the apex's
        // side of the gap and about half of it is left out.
        let lost = without_monoisotopic(scans(&[(100.0, 1e6)]), &[96.0, 97.0, 98.0, 99.0]);
        let feature = trace(&lost, &[100.0], TraceSettings::default()).unwrap();
        assert!(feature.start_rt >= 100.0, "{feature:?}");
        assert_close(feature.area, gaussian_area(1e6) / 2.0, 0.05);
    }

    #[test]
    fn a_shallow_dip_does_not_separate_a_peak() {
        // Peaks of equal height 20 s apart (sigma 5 s): halfway each is
        // exp(-100/50) = 13.5% of its height, so the sum dips to 27% of the
        // apex. A valley deeper than half would separate them; this one is
        // above a fifth, so it is one peak of twice the area.
        let feature = trace(
            &scans(&[(100.0, 1e6), (120.0, 1e6)]),
            &[100.0],
            TraceSettings::default(),
        )
        .unwrap();
        assert_close(feature.area, 2.0 * gaussian_area(1e6), 0.03);
        assert!(feature.end_rt > 130.0, "{feature:?}");
    }

    #[test]
    fn a_neighbouring_peak_is_left_out() {
        // Two peaks 40 s apart; identified on the second.
        let feature = trace(
            &scans(&[(100.0, 4e6), (140.0, 1e6)]),
            &[141.0],
            TraceSettings::default(),
        )
        .unwrap();
        assert_eq!(feature.apex_rt, 140.0);
        assert_close(feature.area, gaussian_area(1e6), 0.03);
    }

    #[test]
    fn identifications_extend_the_trace() {
        // A margin of 10 s from an identification at 70 s ends the trace at
        // 80 s, before the apex; one at 100 s extends it past the peak.
        let short = TraceSettings {
            margin_seconds: 10.0,
            ..TraceSettings::default()
        };
        let peak = scans(&[(100.0, 1e6)]);
        let cut = trace(&peak, &[70.0], short);
        assert!(cut.is_none_or(|f| f.area < gaussian_area(1e6) / 2.0));
        let extended = trace(&peak, &[90.0, 100.0], short).unwrap();
        assert_eq!(extended.apex_rt, 100.0);
    }

    #[test]
    fn isotope_peaks_must_be_within_tolerance() {
        let shifted = |ppm: f64| -> Vec<(f64, Vec<f64>, Vec<f32>)> {
            scans(&[(100.0, 1e6)])
                .into_iter()
                .map(|(rt, mz, intensity)| {
                    (
                        rt,
                        mz.iter().map(|m| m * (1.0 + ppm * 1e-6)).collect(),
                        intensity,
                    )
                })
                .collect()
        };
        assert!(trace(&shifted(8.0), &[100.0], TraceSettings::default()).is_some());
        assert_eq!(
            trace(&shifted(30.0), &[100.0], TraceSettings::default()),
            None
        );
    }

    #[test]
    fn without_an_m_plus_1_peak_there_is_no_feature() {
        let lonely: Vec<_> = scans(&[(100.0, 1e6)])
            .into_iter()
            .map(|(rt, mz, intensity)| (rt, mz[..1].to_vec(), intensity[..1].to_vec()))
            .collect();
        assert_eq!(trace(&lonely, &[100.0], TraceSettings::default()), None);
    }

    #[test]
    fn unsorted_scans_and_other_charges_are_handled() {
        let mut tracer = Tracer::new(TraceSettings::default());
        for (rt, mut mz, mut intensity) in scans(&[(100.0, 1e6)]) {
            if rt == 100.0 {
                tracer.add_identification((1, 2), MZ, 2, rt, None);
                tracer.add_identification((1, 3), MZ, 3, rt, None);
            }
            mz.reverse();
            intensity.reverse();
            tracer.add_ms1(rt, &mz, &intensity, &[]);
        }
        assert_eq!(tracer.ms1_scans(), 241);
        let features = tracer.finish();
        assert_close(features[&(1, 2)].unwrap().area, gaussian_area(1e6), 0.02);
        // At 3+ the isotope spacing differs: nothing at those m/z values.
        assert_eq!(features[&(1, 3)], None);
    }

    #[test]
    fn a_known_mobility_leaves_out_ions_at_the_same_mz_elsewhere_in_mobility() {
        // Scans binned by mobility: the precursor at 1/K0 1.00 spread over
        // three bins (25%, 50%, 25%), and another ion at the same m/z and
        // time at 1.10, twice as intense.
        let mut tracer = Tracer::new(TraceSettings::default());
        for (rt, mz, intensity) in scans(&[(100.0, 1e6)]) {
            if rt == 100.0 {
                tracer.add_identification("known", MZ, 2, rt, Some(1.0));
                tracer.add_identification("unknown", MZ, 2, rt, None);
            }
            let mut peaks: Vec<(f64, f32, f32)> = Vec::new();
            for (&m, &i) in mz.iter().zip(&intensity) {
                for (mobility, share) in [(0.995, 0.25), (1.0, 0.5), (1.005, 0.25), (1.1, 2.0)] {
                    peaks.push((m, i * share, mobility));
                }
            }
            let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
            let intensity: Vec<f32> = peaks.iter().map(|p| p.1).collect();
            let mobility: Vec<f32> = peaks.iter().map(|p| p.2).collect();
            tracer.add_ms1(rt, &mz, &intensity, &mobility);
        }
        let features = tracer.finish();
        let known = features["known"].unwrap();
        assert_close(known.area, gaussian_area(1e6), 0.02);
        assert_eq!(known.apex_ppm, 0.0);
        // Without a mobility every bin counts: the collapsed frame.
        assert_close(
            features["unknown"].unwrap().area,
            3.0 * gaussian_area(1e6),
            0.02,
        );
    }
}
