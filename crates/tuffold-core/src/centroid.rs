//! Profile spectra to centroids: each local maximum becomes one peak, at the
//! vertex of a parabola through the logs of the apex and its two neighbours
//! (exact for Gaussian peaks), or, for a flat top or an apex at the edge, the
//! intensity-weighted m/z of the points above half its height. Its intensity
//! is the summed intensity of the whole peak, out to the valleys either side. Simple and deterministic; enough for a
//! search engine that wants centroids, not a substitute for a vendor's
//! peak picker.

/// Centroid one profile spectrum. `mz` must be finite, sorted ascending, and
/// as long as `intensity` (finite too); anything else returns no peaks.
pub fn centroid(mz: &[f64], intensity: &[f32]) -> (Vec<f64>, Vec<f32>) {
    let finite = mz.iter().all(|v| v.is_finite()) && intensity.iter().all(|v| v.is_finite());
    if !finite || mz.len() != intensity.len() || mz.windows(2).any(|w| w[0] > w[1]) {
        return (Vec::new(), Vec::new());
    }
    let n = intensity.len();
    let mut peaks = Vec::new();
    let mut i = 0;
    while i < n {
        // A plateau of equal values counts as one apex (its weighted centre).
        let mut end = i;
        while end + 1 < n && intensity[end + 1] == intensity[i] {
            end += 1;
        }
        let rising = i == 0 || intensity[i - 1] < intensity[i];
        let falling = end + 1 == n || intensity[end + 1] < intensity[i];
        if intensity[i] > 0.0 && rising && falling {
            peaks.push(peak(mz, intensity, i, end));
        }
        i = end + 1;
    }
    peaks.into_iter().unzip()
}

/// The peak whose apex spans `first..=last`.
fn peak(mz: &[f64], intensity: &[f32], first: usize, last: usize) -> (f64, f32) {
    // Out to the valley (or the zero) on each side.
    let mut left = first;
    while left > 0 && intensity[left - 1] > 0.0 && intensity[left - 1] <= intensity[left] {
        left -= 1;
    }
    // The valley between two peaks belongs to the right-hand one (whose left
    // side takes it above), so no point is counted twice.
    let n = intensity.len();
    let mut right = last;
    while right + 1 < n
        && intensity[right + 1] > 0.0
        && intensity[right + 1] <= intensity[right]
        && (right + 2 == n || intensity[right + 2] <= intensity[right + 1])
    {
        right += 1;
    }
    if first == last && first > 0 && last + 1 < intensity.len() {
        let (a, b, c) = (intensity[first - 1], intensity[first], intensity[first + 1]);
        if a > 0.0 && c > 0.0 {
            let (la, lb, lc) = (f64::from(a).ln(), f64::from(b).ln(), f64::from(c).ln());
            let curvature = la - 2.0 * lb + lc;
            // Even spacing around the apex is assumed by the formula.
            let spacing = (mz[first + 1] - mz[first - 1]) / 2.0;
            if curvature < 0.0 {
                let offset = 0.5 * (la - lc) / curvature;
                let area: f32 = intensity[left..=right].iter().sum();
                return (mz[first] + offset * spacing, area);
            }
        }
    }
    let half = intensity[first] / 2.0;
    let (weighted, weight) =
        (left..=right)
            .filter(|&j| intensity[j] >= half)
            .fold((0.0, 0.0), |(sum, total), j| {
                let w = f64::from(intensity[j]);
                (sum + mz[j] * w, total + w)
            });
    let area: f32 = intensity[left..=right].iter().sum();
    (weighted / weight, area)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Gaussian peak sampled every 0.002 m/z.
    fn profile(peaks: &[(f64, f32)], from: f64, to: f64) -> (Vec<f64>, Vec<f32>) {
        let step = 0.002;
        let points = ((to - from) / step) as usize;
        let sigma = 0.004;
        (0..points)
            .map(|k| {
                let x = from + k as f64 * step;
                let y: f64 = peaks
                    .iter()
                    .map(|&(center, height)| {
                        f64::from(height) * (-((x - center) / sigma).powi(2) / 2.0).exp()
                    })
                    .sum();
                (x, y as f32)
            })
            .unzip()
    }

    #[test]
    fn gaussian_peaks_centroid_to_their_centers() {
        // Centers between sample points, so the weighting has to work.
        let (mz, intensity) = profile(&[(500.2731, 1000.0), (501.2765, 400.0)], 500.0, 502.0);
        let (centroids, heights) = centroid(&mz, &intensity);
        assert_eq!(centroids.len(), 2, "{centroids:?}");
        assert!((centroids[0] - 500.2731).abs() < 1e-5, "{}", centroids[0]);
        assert!((centroids[1] - 501.2765).abs() < 1e-5, "{}", centroids[1]);
        // Areas keep the ratio of the heights (same width).
        assert!((heights[0] / heights[1] - 2.5).abs() < 0.01);
    }

    #[test]
    fn overlapping_peaks_split_at_the_valley() {
        let (mz, intensity) = profile(&[(600.000, 800.0), (600.020, 500.0)], 599.95, 600.07);
        let (centroids, _) = centroid(&mz, &intensity);
        assert_eq!(centroids.len(), 2, "{centroids:?}");
    }

    #[test]
    fn plateaus_and_edges() {
        // A flat top is one peak, centred on the plateau.
        let (mz, intensity) = centroid(&[1.0, 2.0, 3.0, 4.0, 5.0], &[0.0, 5.0, 5.0, 5.0, 0.0]);
        assert_eq!((mz, intensity), (vec![3.0], vec![15.0]));
        // A peak at the very edge still counts.
        assert_eq!(centroid(&[1.0, 2.0], &[3.0, 1.0]).0, [1.0]);
        // Nothing, all zero, or malformed input: no peaks.
        assert!(centroid(&[], &[]).0.is_empty());
        assert!(centroid(&[1.0, 2.0], &[0.0, 0.0]).0.is_empty());
        assert!(centroid(&[2.0, 1.0], &[1.0, 2.0]).0.is_empty());
        assert!(centroid(&[1.0], &[1.0, 2.0]).0.is_empty());
        assert!(centroid(&[1.0, f64::NAN], &[1.0, 2.0]).0.is_empty());
        assert!(centroid(&[1.0, 2.0], &[1.0, f32::NAN]).0.is_empty());
    }

    #[test]
    fn a_valley_is_counted_once() {
        // Two triangles sharing the valley point (2.0 at m/z 3).
        let (mz, area) = centroid(&[1.0, 2.0, 3.0, 4.0, 5.0], &[1.0, 5.0, 2.0, 6.0, 1.0]);
        assert_eq!(mz.len(), 2);
        // Left peak: 1 + 5; right peak: 2 + 6 + 1. Total = the whole profile.
        assert_eq!(area, [6.0, 9.0]);
    }
}
