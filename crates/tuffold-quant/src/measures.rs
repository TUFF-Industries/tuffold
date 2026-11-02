//! Abundance indices besides counts and MS1 intensity, as Scaffold reports them.
//!
//! - **emPAI** (Ishihama et al. 2005): `10^(observed / observable) − 1`,
//!   where observed is the protein's distinct peptide sequences with a
//!   passing spectrum in the sample and observable its tryptic peptides
//!   ([`tuffold_core::digest::observable_peptides`]). Peptides with missed
//!   cleavages or non-tryptic ends count as observed, so observed can pass
//!   observable; the index is reported as computed, not capped.
//! - **Top 3** (Silva et al. 2006): the mean of the three most intense
//!   precursors' normalized MS1 areas in the sample, or of the one or two
//!   there are. Relative only here: Silva's absolute calibration needs a
//!   spiked standard Tuffold doesn't know about.

/// Precursors averaged by Top 3.
pub const TOP_N: usize = 3;

/// emPAI of one protein in one sample; `None` without observable peptides.
pub fn empai(observed: usize, observable: usize) -> Option<f64> {
    (observable > 0).then(|| 10_f64.powf(observed as f64 / observable as f64) - 1.0)
}

/// The mean of the [`TOP_N`] largest areas (linear scale); `None` without
/// any. Non-finite or non-positive areas are ignored.
pub fn top3(areas: impl IntoIterator<Item = f64>) -> Option<f64> {
    let mut areas: Vec<f64> = areas
        .into_iter()
        .filter(|area| area.is_finite() && *area > 0.0)
        .collect();
    if areas.is_empty() {
        return None;
    }
    areas.sort_by(|a, b| b.total_cmp(a));
    let top = &areas[..areas.len().min(TOP_N)];
    Some(top.iter().sum::<f64>() / top.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empai_follows_ishihama() {
        // 10^(2/4) − 1 = √10 − 1 = 2.16227766...
        assert!((empai(2, 4).unwrap() - 2.162_277_660_168_379_5).abs() < 1e-12);
        assert_eq!(empai(0, 4), Some(0.0));
        // 10^(3/3) − 1 = 9: every observable peptide seen.
        assert!((empai(3, 3).unwrap() - 9.0).abs() < 1e-12);
        assert_eq!(empai(1, 0), None);
    }

    #[test]
    fn top3_averages_the_three_most_intense() {
        // 50, 40, 30 of {10, 50, 30, 40}: 40.
        assert_eq!(top3([10.0, 50.0, 30.0, 40.0]), Some(40.0));
        // Fewer than three: their mean.
        assert_eq!(top3([8.0, 4.0]), Some(6.0));
        assert_eq!(top3([f64::NAN, 0.0, -1.0]), None);
        assert_eq!(top3([]), None);
    }
}
