//! The statistics spectral-count comparisons need, from first principles:
//! log-gamma (Lanczos), the regularized incomplete beta function (continued
//! fraction, Lentz), Welch's t-test, Fisher's exact test, and
//! Benjamini-Hochberg q-values. Checked against scipy (see the tests).

/// ln Γ(x) for x > 0 (Lanczos, g = 7, nine coefficients; ~1e-15 relative).
pub fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection: Γ(x) Γ(1 - x) = π / sin(πx).
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + 7.5;
    let series = COEFFICIENTS[1..]
        .iter()
        .enumerate()
        .fold(COEFFICIENTS[0], |sum, (i, c)| {
            sum + c / (x + i as f64 + 1.0)
        });
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// The regularized incomplete beta function I_x(a, b).
pub fn incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    // The continued fraction converges fast for x < (a + 1) / (a + b + 2).
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_fraction(b, a, 1.0 - x) / b
    }
}

fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    const EPSILON: f64 = 1e-15;
    let mut c = 1.0;
    let mut d = 1.0 - (a + b) * x / (a + 1.0);
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut fraction = d;
    for m in 1..=300 {
        let m = f64::from(m);
        for numerator in [
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)),
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0)),
        ] {
            d = 1.0 + numerator * d;
            if d.abs() < TINY {
                d = TINY;
            }
            c = 1.0 + numerator / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            fraction *= d * c;
        }
        if (d * c - 1.0).abs() < EPSILON {
            break;
        }
    }
    fraction
}

/// Two-sided p-value of Student's t at `df` degrees of freedom.
pub fn t_two_sided(t: f64, df: f64) -> f64 {
    incomplete_beta(df / 2.0, 0.5, df / (df + t * t))
}

/// Welch's two-sample t-test: `(t, df, p)`. `None` with fewer than two
/// values in a group, or when both groups have no variance (t undefined).
pub fn welch(a: &[f64], b: &[f64]) -> Option<(f64, f64, f64)> {
    welch_floored(a, b, 0.0, 0.0)
}

/// Welch's t-test with each group's variance at least its floor (see
/// [`counting_floor`]).
pub fn welch_floored(a: &[f64], b: &[f64], floor_a: f64, floor_b: f64) -> Option<(f64, f64, f64)> {
    if a.len() < 2 || b.len() < 2 {
        return None;
    }
    let summary = |values: &[f64], floor: f64| {
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (n, mean, variance.max(floor) / n)
    };
    let (na, mean_a, se_a) = summary(a, floor_a);
    let (nb, mean_b, se_b) = summary(b, floor_b);
    let se = se_a + se_b;
    // No spread at all, up to rounding: nothing to test (a relative bound, so
    // float noise around equal values can't pass for a tiny variance).
    if se.sqrt() <= 1e-12 * mean_a.abs().max(mean_b.abs()).max(1.0) {
        return None;
    }
    let t = (mean_a - mean_b) / se.sqrt();
    let df = se.powi(2) / (se_a.powi(2) / (na - 1.0) + se_b.powi(2) / (nb - 1.0));
    Some((t, df, t_two_sided(t, df)))
}

/// The least variance replicate log2(count + 1) values can have: Poisson
/// counting noise, by the delta method, at the values' mean count `m`:
/// `m / ((m + 1)² ln² 2)`. Spectral counts of 1, 1 against 0, 0 have no
/// spread at all, but that says nothing about their counting noise, and
/// without a floor such rows test as highly significant.
pub fn counting_floor(log2_values: &[f64]) -> f64 {
    if log2_values.is_empty() {
        return 0.0;
    }
    let mean = log2_values.iter().map(|v| v.exp2() - 1.0).sum::<f64>() / log2_values.len() as f64;
    let ln2 = std::f64::consts::LN_2;
    mean.max(0.0) / ((mean.max(0.0) + 1.0).powi(2) * ln2 * ln2)
}

/// Two-sided Fisher's exact test of the 2 x 2 table `[[a, b], [c, d]]`:
/// the summed probability of every table with the same margins that is no
/// more likely than the observed one (scipy's definition).
///
/// Table probabilities are built by the ratio recurrence
/// `P(x + 1) / P(x) = (col1 - x)(row1 - x) / ((x + 1)(total - col1 - row1 + x + 1))`
/// and normalized by their sum, which stays exact for very large counts
/// (no cancellation between huge log-factorials).
pub fn fisher_exact(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let (row1, col1, total) = (a + b, a + c, a + b + c + d);
    let low = row1.saturating_sub(total - col1);
    let high = row1.min(col1);
    let mut ln_p = Vec::with_capacity(usize::try_from(high - low + 1).unwrap_or(0));
    let mut current = 0.0_f64;
    for x in low..=high {
        ln_p.push(current);
        if x < high {
            let numerator = (col1 - x) as f64 * (row1 - x) as f64;
            // x >= low guarantees total + x + 1 >= col1 + row1 (no underflow).
            let denominator = (x + 1) as f64 * (total + x + 1 - col1 - row1) as f64;
            current += (numerator / denominator).ln();
        }
    }
    let peak = ln_p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let weights: Vec<f64> = ln_p.iter().map(|lp| (lp - peak).exp()).collect();
    let sum: f64 = weights.iter().sum();
    let observed = weights[usize::try_from(a - low).unwrap_or(0)];
    // Relative slack so tables tied with the observed one count (as scipy).
    let threshold = observed * (1.0 + 1e-7);
    (weights.iter().filter(|&&w| w <= threshold).sum::<f64>() / sum).min(1.0)
}

/// Benjamini-Hochberg q-values; `None` p-values stay `None` and do not count
/// towards the number of tests.
pub fn benjamini_hochberg(p: &[Option<f64>]) -> Vec<Option<f64>> {
    let mut tested: Vec<(usize, f64)> = p
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.map(|p| (i, p)))
        .collect();
    tested.sort_by(|x, y| x.1.total_cmp(&y.1));
    let m = tested.len() as f64;
    let mut q = vec![None; p.len()];
    let mut lowest = 1.0_f64;
    for (rank, (index, value)) in tested.iter().enumerate().rev() {
        lowest = lowest.min(value * m / (rank as f64 + 1.0));
        q[*index] = Some(lowest);
    }
    q
}

/// Splits of the pooled samples enumerated for an exact Mann–Whitney
/// p-value; beyond this the normal approximation is used.
const EXACT_SPLITS: u64 = 200_000;

/// Two-sided Mann–Whitney U test (Wilcoxon rank-sum): the p-value that
/// `a` and `b` come from the same distribution. Ties get midranks. Exact
/// (every split of the pooled values into groups of these sizes, ties and
/// all) when there are at most [`EXACT_SPLITS`] splits, else the normal
/// approximation with tie and continuity corrections, as scipy's
/// `mannwhitneyu`. `None` when a group is empty. With few replicates the
/// smallest p it can give is large (1/3 for 2 vs 2, 0.1 for 3 vs 3).
pub fn mann_whitney(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.is_empty() || b.is_empty() || a.iter().chain(b).any(|v| !v.is_finite()) {
        return None;
    }
    let (n1, n2) = (a.len(), b.len());
    let n = n1 + n2;
    let ranks = midranks(&a.iter().chain(b).copied().collect::<Vec<_>>());
    let rank_sum_a: f64 = ranks[..n1].iter().sum();
    let offset = (n1 * (n1 + 1)) as f64 / 2.0;
    let mean = (n1 * n2) as f64 / 2.0;
    let observed = (rank_sum_a - offset - mean).abs();
    if splits(n, n1) <= EXACT_SPLITS {
        // Every choice of n1 positions for group a.
        let mut chosen: Vec<usize> = (0..n1).collect();
        let (mut extreme, mut total) = (0_u64, 0_u64);
        loop {
            let sum: f64 = chosen.iter().map(|&i| ranks[i]).sum();
            total += 1;
            if (sum - offset - mean).abs() >= observed - 1e-9 {
                extreme += 1;
            }
            // Next combination in lexicographic order.
            let Some(i) = (0..n1).rev().find(|&i| chosen[i] != i + n - n1) else {
                break;
            };
            chosen[i] += 1;
            for j in i + 1..n1 {
                chosen[j] = chosen[j - 1] + 1;
            }
        }
        return Some(extreme as f64 / total as f64);
    }
    // Normal approximation, variance corrected for ties.
    let mut sorted = ranks.clone();
    sorted.sort_by(f64::total_cmp);
    let mut ties = 0.0;
    let mut start = 0;
    while start < n {
        let end = (start..n)
            .find(|&j| sorted[j] != sorted[start])
            .unwrap_or(n);
        let t = (end - start) as f64;
        ties += t * t * t - t;
        start = end;
    }
    let (n1f, n2f, nf) = (n1 as f64, n2 as f64, n as f64);
    let variance = n1f * n2f / 12.0 * ((nf + 1.0) - ties / (nf * (nf - 1.0)));
    if variance <= 0.0 {
        return Some(1.0);
    }
    let z = ((observed - 0.5).max(0.0)) / variance.sqrt();
    Some(erfc(z / std::f64::consts::SQRT_2).min(1.0))
}

/// The smallest two-sided p-value Mann–Whitney can give for groups of
/// these sizes: both complete separations out of every split. 1/3 for 2 vs
/// 2, 0.1 for 3 vs 3, 0.029 for 4 vs 4.
pub fn mann_whitney_floor(n1: usize, n2: usize) -> f64 {
    if n1 == 0 || n2 == 0 {
        return 1.0;
    }
    (2.0 / splits(n1 + n2, n1) as f64).min(1.0)
}

/// C(n, k), saturating.
fn splits(n: usize, k: usize) -> u64 {
    let k = k.min(n - k);
    (0..k).fold(1_u64, |c, i| {
        c.saturating_mul((n - i) as u64) / (i as u64 + 1)
    })
}

/// 1-based ranks, ties sharing their mean rank.
fn midranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&x, &y| values[x].total_cmp(&values[y]));
    let mut ranks = vec![0.0; values.len()];
    let mut start = 0;
    while start < order.len() {
        let end = (start..order.len())
            .find(|&j| values[order[j]] != values[order[start]])
            .unwrap_or(order.len());
        let rank = (start + end + 1) as f64 / 2.0;
        for &i in &order[start..end] {
            ranks[i] = rank;
        }
        start = end;
    }
    ranks
}

/// The complementary error function (Numerical Recipes' Chebyshev fit,
/// fractional error below 1.2e-7: plenty for a p-value).
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let polynomial = -z * z - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
    let value = t * polynomial.exp();
    if x >= 0.0 { value } else { 2.0 - value }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mann_whitney_matches_hand_and_scipy() {
        // Complete separation: 2 of C(6,3) = 20 splits are as extreme
        // (scipy exact: 0.1); 2 vs 2: 2 of 6 (scipy: 1/3), the best it
        // can ever do with two replicates.
        close(
            mann_whitney(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap(),
            0.1,
            1e-12,
        );
        close(
            mann_whitney(&[1.0, 2.0], &[3.0, 4.0]).unwrap(),
            1.0 / 3.0,
            1e-12,
        );
        // Ties: pooled 1,2,2,2,3,4 rank 1,3,3,3,5,6; a's ranks sum to 7, so
        // U = 1, |U − 4.5| = 3.5. Splits as extreme: rank sum ≤ 7 ({1,3,3},
        // 3 ways) or ≥ 14 ({3,5,6}, 3 ways): 6 of 20 = 0.3 (scipy
        // permutation_test on U: 0.3).
        close(
            mann_whitney(&[1.0, 2.0, 2.0], &[2.0, 3.0, 4.0]).unwrap(),
            0.3,
            1e-12,
        );
        // 12 vs 12 with ties (C(24,12) splits: normal approximation); scipy
        // mannwhitneyu(method="asymptotic"): 0.002890129724154294.
        let a = [1.0, 2.0, 2.0, 3.0, 5.0, 5.0, 6.0, 7.0, 8.0, 8.0, 9.0, 10.0];
        let b = [
            4.0, 6.0, 7.0, 9.0, 10.0, 11.0, 11.0, 12.0, 13.0, 14.0, 14.0, 15.0,
        ];
        let p = mann_whitney(&a, &b).unwrap();
        assert!((p / 0.002_890_129_724_154_294 - 1.0).abs() < 1e-6, "{p}");
        assert_eq!(mann_whitney(&[2.0, 2.0], &[2.0, 2.0]), Some(1.0));
        // 2 / C(n1 + n2, n1).
        close(mann_whitney_floor(2, 2), 1.0 / 3.0, 1e-12);
        close(mann_whitney_floor(3, 3), 0.1, 1e-12);
        close(mann_whitney_floor(4, 4), 2.0 / 70.0, 1e-12);
        assert_eq!(mann_whitney_floor(1, 1), 1.0);
        assert_eq!(mann_whitney(&[], &[1.0]), None);
    }

    const LN2SQ: f64 = std::f64::consts::LN_2 * std::f64::consts::LN_2;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn ln_gamma_matches_known_values() {
        // Γ(0.5) = √π; Γ(10) = 9! = 362880; scipy gammaln(1e6).
        close(ln_gamma(0.5), std::f64::consts::PI.sqrt().ln(), 1e-13);
        close(ln_gamma(10.0), 362_880_f64.ln(), 1e-13);
        close(ln_gamma(1e6), 12_815_504.569_147_611, 1e-13);
    }

    #[test]
    fn incomplete_beta_matches_the_binomial_sum() {
        // I_0.4(2, 3) = Σ_{j=2..4} C(4,j) 0.4^j 0.6^(4-j)
        //             = 6·0.16·0.36 + 4·0.064·0.6 + 0.0256 = 0.5248.
        close(incomplete_beta(2.0, 3.0, 0.4), 0.5248, 1e-12);
        assert_eq!(incomplete_beta(2.0, 3.0, 0.0), 0.0);
        assert_eq!(incomplete_beta(2.0, 3.0, 1.0), 1.0);
    }

    #[test]
    fn welch_matches_scipy() {
        // [1,2,3] vs [4,5,6]: means 2 and 5, variances 1, se² = 2/3,
        // t = -3/√(2/3) = -3.674235, df = (2/3)² / (2 · (1/3)² / 2) = 4;
        // scipy ttest_ind(equal_var=False): p = 0.0213116411.
        let (t, df, p) = welch(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap();
        close(t, -3.674_234_614_174_767, 1e-12);
        close(df, 4.0, 1e-12);
        close(p, 0.021_311_641_128_756_713, 1e-9);
        // Unequal sizes and variances: scipy p = 0.3943345634897475.
        let (_, _, p) = welch(&[1.0, 2.5, 2.9, 4.1], &[3.0, 5.5]).unwrap();
        close(p, 0.394_334_563_489_747_5, 1e-9);
        assert_eq!(welch(&[1.0], &[2.0, 3.0]), None);
        assert_eq!(welch(&[2.0, 2.0], &[2.0, 2.0]), None);
        // A floor lifts a variance that's smaller, and leaves a larger one.
        let (t, _, _) = welch_floored(&[1.0, 1.0], &[0.0, 0.0], 0.5, 0.0).unwrap();
        assert!((t - 1.0 / (0.5_f64 / 2.0).sqrt()).abs() < 1e-12);
        let plain = welch(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap();
        assert_eq!(
            welch_floored(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], 0.1, 0.1),
            Some(plain)
        );
    }

    #[test]
    fn the_counting_floor_follows_the_delta_method() {
        // One count on average: 1 / (4 ln² 2) ≈ 0.5203.
        assert!((counting_floor(&[1.0, 1.0]) - 1.0 / (4.0 * LN2SQ)).abs() < 1e-12);
        // Zero counts carry no counting noise; large counts little.
        assert_eq!(counting_floor(&[0.0, 0.0]), 0.0);
        // m = 100: 100 / (101² ln² 2) ≈ 0.0204.
        assert!((counting_floor(&[(101.0_f64).log2()]) - 100.0 / (10_201.0 * LN2SQ)).abs() < 1e-12);
        assert_eq!(counting_floor(&[]), 0.0);
    }

    #[test]
    fn fisher_matches_the_hand_count_and_scipy() {
        // Margins 4/4 of 8: P(x) ∝ C(4,x)C(4,4-x) = 1, 16, 36, 16, 1 (/70).
        // Observed x = 3 (16/70); tables no likelier: x = 0, 1, 3, 4 -> 34/70.
        close(fisher_exact(3, 1, 1, 3), 34.0 / 70.0, 1e-12);
        close(
            fisher_exact(10, 2, 3, 15),
            0.000_536_724_119_143_435_7,
            1e-9,
        );
        close(fisher_exact(0, 5, 100, 95), 0.059_383_210_468_586_585, 1e-9);
        // Huge margins stay exact: (0, 1, 1e6, 1e6) is as likely as it gets
        // either way, p = 1 (scipy 1.0).
        assert_eq!(fisher_exact(0, 1, 1_000_000, 1_000_000), 1.0);
        assert_eq!(fisher_exact(0, 0, 0, 0), 1.0);
        // Margins over half the total (low > 0): [[5, 0], [1, 0]], scipy 1.0.
        assert_eq!(fisher_exact(5, 0, 1, 0), 1.0);
        // [[9, 1], [1, 9]]: margins 10/10 of 20; P(x) ∝ C(10,x)C(10,10-x);
        // x = 0, 1, 9, 10 are no likelier: (1 + 100 + 100 + 1) / C(20,10) =
        // 202 / 184756 (scipy 0.001093333910671372).
        close(fisher_exact(9, 1, 1, 9), 202.0 / 184_756.0, 1e-12);
    }

    #[test]
    fn benjamini_hochberg_steps_down_from_the_largest() {
        // Sorted 0.01, 0.03, 0.04, 0.2 -> p·4/rank = 0.04, 0.06, 0.0533, 0.2;
        // monotone from the end: 0.04, 0.0533, 0.0533, 0.2 (scipy agrees).
        let q = benjamini_hochberg(&[Some(0.01), Some(0.04), Some(0.03), None, Some(0.2)]);
        let q: Vec<f64> = q.into_iter().map(|q| q.unwrap_or(-1.0)).collect();
        for (actual, expected) in q.iter().zip([0.04, 0.16 / 3.0, 0.16 / 3.0, -1.0, 0.2]) {
            close(*actual, expected, 1e-12);
        }
    }
}
