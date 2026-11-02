//! Protein intensities by Tukey's median polish (Tukey 1977, *Exploratory
//! Data Analysis*), from one protein's peptide intensities across samples.
//! It is how FlashLFQ and MSstats summarize proteins.
//!
//! The protein's peptides × samples table of log2 intensities is split into
//! an overall level, an effect for each peptide (how well it ionizes) and an
//! effect for each sample (the protein's abundance there), plus what is left
//! over. Medians do the splitting, by sweeping row medians and column medians
//! out of the table in turn until nothing moves. A sample's value is the
//! overall level plus its effect, and its level is set so the profile's total
//! matches the protein's summed intensity, as it is for MaxLFQ.
//!
//! The fit is robust: one peptide interfered with in one sample (a
//! co-eluting ion, a mis-traced peak) lands in that cell's left-over and
//! moves no effect, because a median ignores one extreme value among several.
//! MaxLFQ fits a line through pairwise medians by least squares, and when two
//! samples share only one or two peptides, a single bad peptide sets their
//! ratio outright and drags every sample connected to them. Median polish
//! looks at all the peptides of a sample at once.
//!
//! Only samples connected by shared peptides (directly or through a chain)
//! can be compared, so a protein gets values in its largest such group of
//! samples (ties: the larger total intensity) and none elsewhere, chosen as
//! MaxLFQ chooses it. If one run caught peptide A and another caught peptide
//! B, comparing them would compare two peptides' ionization efficiencies
//! rather than the protein's abundance. A protein seen in one sample only
//! keeps its summed intensity there. Missing cells are skipped, never filled.

use crate::pairwise::{largest_component, median, pairwise_ratios, summed_intensities};

/// Sweeps of rows then columns before the fit is accepted as it stands.
const MAX_ITERATIONS: usize = 10;
/// An iteration that moves no effect by more than this (log2) has converged.
const CONVERGED: f64 = 1e-9;

/// Log2 protein intensities from `peptides[p][s]`, the log2 intensity of
/// peptide `p` in sample `s` (`None` where not measured).
///
/// Samples get values only in the largest connected group of samples. There,
/// peptides measured in fewer than two of the group's samples are left out
/// (they say nothing between samples) unless that would leave none, the
/// table is polished (rows first, then columns, at most [`MAX_ITERATIONS`]
/// times), and the profile is shifted so Σ 2^x over the group equals the
/// group's summed intensity over all the protein's peptides.
pub(crate) fn median_polish(peptides: &[Vec<Option<f64>>], samples: usize) -> Vec<Option<f64>> {
    debug_assert!(peptides.iter().all(|p| p.len() == samples));
    // Only whether two samples share a peptide matters here, not the ratio.
    let ratio = pairwise_ratios(peptides, samples, median, 1);
    let summed = summed_intensities(peptides, samples);
    let mut profile = vec![None; samples];
    if let Some((component, target)) = largest_component(&ratio, &summed) {
        let table = component_table(peptides, &component);
        let fitted = polish(&table, component.len());
        // Level: Σ 2^x over the component equals Σ of the summed intensities.
        let current: f64 = fitted.iter().map(|&x| x.exp2()).sum();
        let shift = target.log2() - current.log2();
        for (&s, x) in component.iter().zip(fitted) {
            profile[s] = Some(x + shift);
        }
    }
    profile
}

/// The peptides × component-samples table to polish: peptides measured in at
/// least two of the component's samples, or all it measures at all when none
/// is.
fn component_table(peptides: &[Vec<Option<f64>>], component: &[usize]) -> Vec<Vec<Option<f64>>> {
    let rows: Vec<Vec<Option<f64>>> = peptides
        .iter()
        .map(|p| component.iter().map(|&s| p[s]).collect::<Vec<_>>())
        .filter(|row| row.iter().any(Option::is_some))
        .collect();
    let repeated: Vec<Vec<Option<f64>>> = rows
        .iter()
        .filter(|row| row.iter().flatten().count() >= 2)
        .cloned()
        .collect();
    if repeated.is_empty() { rows } else { repeated }
}

/// Tukey's median polish of `table` (`columns` wide, `None` cells skipped),
/// like R's `medpolish`: row medians out, then column medians out, repeated.
/// After each column sweep the median of the column effects moves into the
/// overall level. Returns overall + column effect for each column.
fn polish(table: &[Vec<Option<f64>>], columns: usize) -> Vec<f64> {
    let mut left_over = table.to_vec();
    let mut column_effect = vec![0.0; columns];
    let mut overall = 0.0;
    for _ in 0..MAX_ITERATIONS {
        // Only the sizes of the moves matter for convergence; a row's own
        // effect (how well its peptide ionizes) is never needed again.
        let mut largest_move = 0.0_f64;
        for row in left_over.iter_mut() {
            let delta = typical(row.iter().flatten().copied());
            row.iter_mut().flatten().for_each(|cell| *cell -= delta);
            largest_move = largest_move.max(delta.abs());
        }
        for (j, effect) in column_effect.iter_mut().enumerate() {
            let delta = typical(left_over.iter().filter_map(|row| row[j]));
            left_over
                .iter_mut()
                .filter_map(|row| row[j].as_mut())
                .for_each(|cell| *cell -= delta);
            *effect += delta;
            largest_move = largest_move.max(delta.abs());
        }
        let shift = typical(column_effect.iter().copied());
        column_effect.iter_mut().for_each(|effect| *effect -= shift);
        overall += shift;
        largest_move = largest_move.max(shift.abs());
        if largest_move <= CONVERGED {
            break;
        }
    }
    column_effect
        .iter()
        .map(|effect| overall + effect)
        .collect()
}

/// The median of some values (even counts: the mean of the two middle ones);
/// 0 for none, so an empty row or column moves nothing.
fn typical(values: impl Iterator<Item = f64>) -> f64 {
    let mut sorted: Vec<f64> = values.collect();
    if sorted.is_empty() {
        return 0.0;
    }
    sorted.sort_by(f64::total_cmp);
    median(&sorted)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn log(values: &[Option<f64>]) -> Vec<Option<f64>> {
        values.iter().map(|v| v.map(f64::log2)).collect()
    }

    fn linear(profile: &[Option<f64>]) -> Vec<Option<f64>> {
        profile.iter().map(|v| v.map(f64::exp2)).collect()
    }

    fn assert_close(actual: &[Option<f64>], expected: &[Option<f64>]) {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected) {
            match (a, e) {
                (Some(a), Some(e)) => {
                    assert!((a - e).abs() < 1e-9 * e.abs().max(1.0), "{a} vs {e}")
                }
                _ => assert_eq!(a, e),
            }
        }
    }

    /// Σ 2^x of the profile's values.
    fn total(profile: &[Option<f64>]) -> f64 {
        profile.iter().flatten().map(|x| x.exp2()).sum()
    }

    #[test]
    fn constant_row_offsets_leave_the_sample_differences_exact() {
        // Rows are already log2: p1 = [10, 11, 12], p2 = [20, 21, 22].
        // Row sweep: medians 11 and 21, leaving [-1, 0, 1] in both rows.
        // Column sweep: the medians are -1, 0, 1, which empties the table;
        // column effects [-1, 0, 1], whose median 0 moves nothing into the
        // overall level. The next iteration moves nothing: converged.
        // So the profile steps up by exactly 1 then 1 (2 over the two).
        let peptides = vec![
            vec![Some(10.0), Some(11.0), Some(12.0)],
            vec![Some(20.0), Some(21.0), Some(22.0)],
        ];
        let profile = median_polish(&peptides, 3);
        let (a, b, c) = (
            profile[0].unwrap(),
            profile[1].unwrap(),
            profile[2].unwrap(),
        );
        assert!((b - a - 1.0).abs() < 1e-12 && (c - a - 2.0).abs() < 1e-12);
        // Level: sample j sums 2^(10+j) + 2^(20+j) = 2^j · 1,049,600, so the
        // total is 7 · 1,049,600 = 7,347,200 and, with parts 1 : 2 : 4, the
        // profile reproduces each sample's own sum.
        assert!((total(&profile) - 7_347_200.0).abs() < 1e-6);
        assert_close(
            &linear(&profile),
            &[Some(1_049_600.0), Some(2_099_200.0), Some(4_198_400.0)],
        );
    }

    #[test]
    fn one_interfered_peptide_in_one_sample_moves_nothing() {
        // p1 = [10, 11, 12], p2 = [20, 21, 22], p3 = [30, 31, 42]: p3's third
        // value is 10 too high (interference).
        // Row sweep: medians 11, 21, 31, leaving [-1, 0, 1], [-1, 0, 1] and
        // [-1, 0, 11].
        // Column sweep: column medians median(-1,-1,-1) = -1,
        // median(0,0,0) = 0, median(1,1,11) = 1, leaving every cell 0 except
        // p3's third, 10. Column effects [-1, 0, 1], median 0.
        // Second iteration: p3's row median of [0, 0, 10] is 0 and every
        // column median is 0, so nothing moves: converged, and the profile
        // steps by 1 and 1 as p1 and p2 say.
        // (Means instead of medians: the row-centred third column is
        // 1, 1, 11, mean 4.33 against -1 and 0 in the others, so the third
        // sample would sit 5.33 above the first instead of 2.)
        let peptides = vec![
            vec![Some(10.0), Some(11.0), Some(12.0)],
            vec![Some(20.0), Some(21.0), Some(22.0)],
            vec![Some(30.0), Some(31.0), Some(42.0)],
        ];
        let profile = median_polish(&peptides, 3);
        let (a, b, c) = (
            profile[0].unwrap(),
            profile[1].unwrap(),
            profile[2].unwrap(),
        );
        assert!((b - a - 1.0).abs() < 1e-12 && (c - b - 1.0).abs() < 1e-12);
        // Level: still the summed intensity, outlier included (it is there).
        let summed: f64 = peptides.iter().flatten().flatten().map(|x| x.exp2()).sum();
        assert!((total(&profile) - summed).abs() < 1e-6 * summed);
    }

    #[test]
    fn missing_cells_are_skipped_and_a_single_sample_peptide_is_dropped() {
        // p1 = [10, 11, 12], p2 = [-, 21, 22], p3 = [30, -, -] (one sample:
        // dropped, though its 2^30 stays in sample 0's summed intensity).
        // Row sweep: p1 median 11 -> [-1, 0, 1]; p2 median of (21, 22) is
        // 21.5 -> [-, -0.5, 0.5].
        // Column sweep: medians -1; median(0, -0.5) = -0.25;
        // median(1, 0.5) = 0.75. Left over: p1 [0, 0.25, 0.25],
        // p2 [-, -0.25, -0.25]; column effects [-1, -0.25, 0.75].
        // Median of those is -0.25: effects [-0.75, 0, 1], overall -0.25.
        // Second iteration: row medians 0.25 and -0.25 leave p1
        // [-0.25, 0, 0] and p2 [-, 0, 0]; column medians -0.25, 0, 0 give
        // effects [-1, 0, 1], whose median is 0. Third iteration moves
        // nothing. Profile = -0.25 + [-1, 0, 1]: steps of 1 and 1 (2 over
        // the two), as p1 says; p2's missing cell changed nothing.
        let peptides = vec![
            vec![Some(10.0), Some(11.0), Some(12.0)],
            vec![None, Some(21.0), Some(22.0)],
            vec![Some(30.0), None, None],
        ];
        let profile = median_polish(&peptides, 3);
        let (a, b, c) = (
            profile[0].unwrap(),
            profile[1].unwrap(),
            profile[2].unwrap(),
        );
        assert!((b - a - 1.0).abs() < 1e-12 && (c - a - 2.0).abs() < 1e-12);
        // Level: sums 2^10 + 2^30, 2^11 + 2^21, 2^12 + 2^22.
        let target = (1u64 << 10) + (1 << 30) + (1 << 11) + (1 << 21) + (1 << 12) + (1 << 22);
        assert!((total(&profile) - target as f64).abs() < 1e-6 * target as f64);
    }

    #[test]
    fn only_the_largest_connected_group_of_samples_gets_values() {
        // A and B share p1; C and D share p2; the groups share nothing, so
        // comparing them would compare p1 with p2. Both have two samples, so
        // the tie goes to the larger total: A + B = 2 + 6 = 8, C + D =
        // 5 + 9 = 14.
        // The group {C, D} has one row [log2 5, log2 9]; its median leaves
        // ±(log2 9 - log2 5)/2, which the column sweep takes as the
        // effects. So the profile is in ratio 5 : 9, and the level (sum 14)
        // makes it exactly 5 and 9.
        let peptides = vec![
            log(&[Some(2.0), Some(6.0), None, None]),
            log(&[None, None, Some(5.0), Some(9.0)]),
        ];
        let profile = linear(&median_polish(&peptides, 4));
        assert_close(&profile, &[None, None, Some(5.0), Some(9.0)]);
        // A larger group beats a larger total: {A, B, C} (3 samples, via p1
        // and p2) wins over D alone.
        let peptides = vec![
            log(&[Some(2.0), Some(4.0), None, None]),
            log(&[None, Some(4.0), Some(8.0), None]),
            log(&[None, None, None, Some(1000.0)]),
        ];
        let profile = median_polish(&peptides, 4);
        assert!(profile[..3].iter().all(Option::is_some) && profile[3].is_none());
    }

    #[test]
    fn a_protein_in_one_sample_keeps_its_summed_intensity() {
        // Peptides of 3 and 5 in the second sample only: every peptide is in
        // one sample, so none is dropped (that would leave none); a single
        // column polishes to 0, and the level makes it 3 + 5 = 8.
        let peptides = vec![log(&[None, Some(3.0), None]), log(&[None, Some(5.0), None])];
        let profile = linear(&median_polish(&peptides, 3));
        assert_close(&profile, &[None, Some(8.0), None]);
        // And a protein with no measurements has no values.
        assert_eq!(median_polish(&[], 2), vec![None, None]);
    }
}
