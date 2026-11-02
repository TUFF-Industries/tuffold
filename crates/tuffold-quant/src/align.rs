//! Retention-time alignment between two runs, for match-between-runs. Given
//! the precursors identified in both a donor run and an acceptor run (one
//! retention time each), it learns how the acceptor's gradient drifts from
//! the donor's, so a precursor identified only in the donor can be looked
//! for where it should elute in the acceptor.
//!
//! The drift is a shift (acceptor − donor) that changes slowly along the
//! gradient, so it is a running median of the shifts in donor-RT order. A
//! median, not a mean or a fitted curve, because some pairs are wrong (the
//! same sequence identified from a different analyte); a window that holds
//! mostly right pairs ignores the wrong ones however far off they are. No
//! shape is assumed, so a column that sped up halfway through is followed as
//! well as a steady drift.
//!
//! The window is the larger of 25 pairs and 5% of the pairs, centred on each
//! pair and shrinking to what exists at the ends. Pairs sharing a donor RT
//! are merged into one knot (the median of their running medians) so the
//! mapping is a function.

use std::cmp::Ordering;

/// Fewest finite pairs an alignment is learned from. With fewer, a median
/// window would be most of the data and the shift untrustworthy; no
/// alignment is better than a guessed one.
pub const MIN_PAIRS: usize = 50;

/// Narrowest running-median window, in pairs.
const MIN_WINDOW: usize = 25;

/// A window holds one pair in this many of all pairs, when that exceeds
/// [`MIN_WINDOW`].
const WINDOW_DIVISOR: usize = 20;

/// Scale turning a median absolute deviation into a standard deviation for
/// normally distributed values.
const MAD_TO_SD: f64 = 1.4826;

/// How one run's retention times map onto another's.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// (donor_rt, shift) knots in donor-RT order; shift = acceptor − donor.
    knots: Vec<(f64, f64)>,
    /// Robust spread of the residuals (1.4826 × their median absolute value), seconds.
    pub spread: f64,
    /// Pairs the alignment rests on.
    pub pairs: usize,
}

impl Alignment {
    /// Acceptor RT for a donor RT: linear interpolation between knots;
    /// constant shift beyond the ends (the gradient's edges are not
    /// extrapolated). A NaN donor RT gives NaN.
    pub fn map(&self, donor_rt: f64) -> f64 {
        if donor_rt.is_nan() {
            return f64::NAN;
        }
        let after = self.knots.partition_point(|&(x, _)| x <= donor_rt);
        let shift = match (after.checked_sub(1), self.knots.get(after)) {
            (None, _) => self.knots.first().map_or(0.0, |k| k.1),
            (Some(before), None) => self.knots[before].1,
            (Some(before), Some(&(x1, s1))) => {
                let (x0, s0) = self.knots[before];
                // `x0 <= donor_rt < x1`, so the span is positive.
                s0 + (donor_rt - x0) / (x1 - x0) * (s1 - s0)
            }
        };
        donor_rt + shift
    }
}

/// Learn the alignment from `(donor_rt, acceptor_rt)` pairs in seconds, one
/// per precursor identified in both runs. Pairs with a non-finite time are
/// dropped. `None` when fewer than [`MIN_PAIRS`] remain.
pub fn align(pairs: &[(f64, f64)]) -> Option<Alignment> {
    let mut sorted: Vec<(f64, f64)> = pairs
        .iter()
        .copied()
        .filter(|(d, a)| d.is_finite() && a.is_finite())
        .collect();
    if sorted.len() < MIN_PAIRS {
        return None;
    }
    sorted.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));

    let n = sorted.len();
    let shifts: Vec<f64> = sorted.iter().map(|&(d, a)| a - d).collect();
    let half = MIN_WINDOW.max(n / WINDOW_DIVISOR) / 2;
    let smoothed: Vec<f64> = (0..n)
        .map(|i| median(&shifts[i.saturating_sub(half)..=(i + half).min(n - 1)]))
        .collect();

    let mut knots: Vec<(f64, f64)> = Vec::with_capacity(n);
    let mut start = 0;
    while start < n {
        let donor = sorted[start].0;
        let end = start
            + sorted[start..]
                .iter()
                .take_while(|p| p.0.total_cmp(&donor) == Ordering::Equal)
                .count();
        knots.push((donor, median(&smoothed[start..end])));
        start = end;
    }

    let mut alignment = Alignment {
        knots,
        spread: 0.0,
        pairs: n,
    };
    let residuals: Vec<f64> = sorted
        .iter()
        .map(|&(d, a)| (a - alignment.map(d)).abs())
        .collect();
    alignment.spread = MAD_TO_SD * median(&residuals);
    Some(alignment)
}

/// Median of a non-empty slice of finite values; the mean of the middle two
/// when the count is even.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[mid]
    } else {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 60 pairs at donor 0, 10, …, 590 s, each shifted by `shift(donor)`.
    fn pairs_with(shift: impl Fn(f64) -> f64) -> Vec<(f64, f64)> {
        (0..60)
            .map(|i| {
                let donor = 10.0 * i as f64;
                (donor, donor + shift(donor))
            })
            .collect()
    }

    #[test]
    fn constant_shift_maps_everywhere() {
        // Every shift is 30, so every running median is 30 and every knot is
        // (donor, 30). Interpolating between equal shifts gives 30 and beyond
        // the ends the end shift (30) applies: map(d) = d + 30 for any d.
        // Every residual is acceptor − (donor + 30) = 0, so the spread is
        // 1.4826 × 0 = 0.
        let a = align(&pairs_with(|_| 30.0)).unwrap();
        assert_eq!(a.map(255.0), 285.0); // inside, between knots
        assert_eq!(a.map(-100.0), -70.0); // below the range
        assert_eq!(a.map(1000.0), 1030.0); // above the range
        assert_eq!(a.spread, 0.0);
        assert_eq!(a.pairs, 60);
    }

    #[test]
    fn piecewise_shift_is_followed() {
        // n = 60, so w = max(25, 60 / 20 = 3) = 25 and the window is
        // [i − 12, i + 12]. Pair i has donor 10i; its shift is 0 for i < 30
        // (donor < 300) and 60 for i ≥ 30.
        // Donor 100 is i = 10: window [0, 22], all i < 30, so every shift is
        // 0, the median is 0, and map(100) = 100 + 0.
        // Donor 500 is i = 50: window [38, 62] clipped to [38, 59], all
        // i ≥ 30, so every shift is 60 and map(500) = 500 + 60.
        let a = align(&pairs_with(|d| if d < 300.0 { 0.0 } else { 60.0 })).unwrap();
        assert_eq!(a.map(100.0), 100.0);
        assert_eq!(a.map(500.0), 560.0);
    }

    #[test]
    fn outliers_do_not_move_the_alignment() {
        // Pairs with shift 30, except i = 5, 17, 29, 41, 53 (spaced 12
        // apart), whose acceptor is donor + 600. Windows hold 25 pairs
        // (fewer at the ends, never under 13: [0, 12] at i = 0). Two
        // outliers are 12 apart, so a window holding 3 of them spans at
        // least 24 + 1 = 25 pairs and holds 22 right ones; a window of 13 to
        // 24 pairs holds at most 2 outliers. In every window the right
        // pairs (shift 30) are the clear majority, so the middle value (or
        // both middle values, for an even count) is 30 and every knot is
        // (donor, 30). map(d) = d + 30 everywhere.
        // Residuals: 0 for the 55 right pairs and 600 − 30 = 570 for the 5
        // outliers; the median absolute residual is the 30th/31st of 60
        // sorted values, both 0, so the spread is 0.
        let outliers = [5, 17, 29, 41, 53];
        let pairs: Vec<(f64, f64)> = (0..60)
            .map(|i| {
                let donor = 10.0 * i as f64;
                let shift = if outliers.contains(&i) { 600.0 } else { 30.0 };
                (donor, donor + shift)
            })
            .collect();
        let a = align(&pairs).unwrap();
        for i in 0..60 {
            let donor = 10.0 * i as f64;
            assert_eq!(a.map(donor), donor + 30.0, "knot {i}");
        }
        assert_eq!(a.map(255.0), 285.0);
        assert_eq!(a.spread, 0.0);
        assert_eq!(a.pairs, 60);
    }

    #[test]
    fn too_few_pairs_give_none() {
        // 49 pairs are under MIN_PAIRS = 50.
        let few: Vec<(f64, f64)> = pairs_with(|_| 30.0).into_iter().take(49).collect();
        assert!(align(&few).is_none());
        // 60 pairs with 11 NaN donor RTs leave 60 − 11 = 49 finite pairs.
        let mut nan = pairs_with(|_| 30.0);
        for p in nan.iter_mut().take(11) {
            p.0 = f64::NAN;
        }
        assert!(align(&nan).is_none());
        // One NaN fewer (10) leaves 50, which is enough.
        nan[10].0 = 100.0;
        assert_eq!(align(&nan).unwrap().pairs, 50);
    }

    #[test]
    fn duplicate_donor_times_stay_a_function() {
        // 60 pairs with shift 10 at donors 0, 10, …, 590, plus a 61st at the
        // same donor as i = 30 (300 s) with shift 20: (300, 320). Sorted by
        // donor then acceptor, (300, 310) is index 30 and (300, 320) index
        // 31; n = 61, w = max(25, 61 / 20 = 3) = 25, window ±12.
        // Index 30: [18, 42] holds 25 pairs, one with shift 20 and 24 with
        // shift 10, so the median is 10. Index 31: [19, 43], likewise 10.
        // No other window holds more than that one 20, so every smoothed
        // shift is 10. The two knots at 300 merge to the median of (10, 10)
        // = 10, so map(300) = 310 and map(295) = 305.
        // Residuals: 0 for 60 pairs and 320 − 310 = 10 for one, so the
        // median absolute residual is 0 and the spread is 0.
        let mut pairs = pairs_with(|_| 10.0);
        pairs.push((300.0, 320.0));
        let a = align(&pairs).unwrap();
        assert_eq!(a.map(300.0), 310.0);
        assert_eq!(a.map(295.0), 305.0);
        assert_eq!(a.map(305.0), 315.0);
        assert_eq!(a.spread, 0.0);
        assert_eq!(a.pairs, 61);
    }
}
