//! Pairwise ratios between samples, and what is built on them.
//!
//! The idea is MaxLFQ's (Cox et al. 2014, *Mol Cell Proteomics* 13:2513):
//! compare two samples only on the peptides both measured, take a typical
//! log2 ratio of those peptides, and fit one level per sample to every such
//! pairwise ratio by least squares. Summing peptide intensities misleads
//! whenever samples identify different peptides (a sample that caught three
//! more peptides looks more abundant); a ratio over shared peptides doesn't.
//!
//! This module holds the pieces: [`pairwise_ratios`] (the median or
//! half-sample mode of each pair's differences), the connected components of
//! samples those ratios link, [`largest_component`] (the rule protein
//! intensities use: only samples connected by shared peptides, directly or
//! through a chain of pairs, can be compared, so a protein gets values in its
//! largest such group, ties to the larger total intensity, and none elsewhere;
//! one shared peptide is enough to connect two samples), and the
//! least-squares fit behind [`relative_levels`], which normalization uses
//! over all peptides at once. Protein intensities themselves are built by
//! median polish in `polish.rs`.

/// Log2 of each sample's summed linear intensity over all peptides (`None`
/// where it measured none): what a protein keeps when nothing can be compared.
pub(crate) fn summed_intensities(
    peptides: &[Vec<Option<f64>>],
    samples: usize,
) -> Vec<Option<f64>> {
    (0..samples)
        .map(|s| {
            let total: f64 = peptides.iter().filter_map(|p| p[s]).map(f64::exp2).sum();
            (total > 0.0).then(|| total.log2())
        })
        .collect()
}

/// The connected component of samples (by `ratio`) to give values to, with
/// the summed linear intensity of its samples (the level its profile is set
/// to). The largest component wins; a tie goes to the larger total.
pub(crate) fn largest_component(
    ratio: &[Vec<Option<f64>>],
    summed: &[Option<f64>],
) -> Option<(Vec<usize>, f64)> {
    let total = |component: &[usize]| -> f64 {
        component
            .iter()
            .filter_map(|&s| summed[s])
            .map(f64::exp2)
            .sum()
    };
    components(ratio, summed)
        .into_iter()
        .max_by(|a, b| a.len().cmp(&b.len()).then(total(a).total_cmp(&total(b))))
        .map(|component| {
            let target = total(&component);
            (component, target)
        })
}

/// Each sample's log2 level relative to the others in its connected
/// component (centred on 0 there), from the half-sample mode of pairwise
/// ratios over the peptides two samples share (pairs sharing fewer than
/// `min_shared` don't count), fitted by least squares. `None` for samples
/// connected to no other. Normalization uses this over all peptides at once.
pub(crate) fn relative_levels(
    peptides: &[Vec<Option<f64>>],
    samples: usize,
    min_shared: usize,
) -> Vec<Option<f64>> {
    debug_assert!(peptides.iter().all(|p| p.len() == samples));
    let ratio = pairwise_ratios(peptides, samples, half_sample_mode, min_shared);
    let present: Vec<Option<f64>> = (0..samples)
        .map(|s| peptides.iter().any(|p| p[s].is_some()).then_some(0.0))
        .collect();
    let mut levels = vec![None; samples];
    for component in components(&ratio, &present) {
        if component.len() < 2 {
            continue;
        }
        let fitted = fit(&component, &ratio);
        let mean = fitted.iter().sum::<f64>() / fitted.len() as f64;
        for (&s, x) in component.iter().zip(fitted) {
            levels[s] = Some(x - mean);
        }
    }
    levels
}

/// `ratio[j][k]`: median of log2(peptide in j) − log2(peptide in k) over
/// peptides measured in both; `None` when they share none.
pub(crate) fn pairwise_ratios(
    peptides: &[Vec<Option<f64>>],
    samples: usize,
    summary: fn(&[f64]) -> f64,
    min_shared: usize,
) -> Vec<Vec<Option<f64>>> {
    let mut ratio = vec![vec![None; samples]; samples];
    for j in 0..samples {
        for k in j + 1..samples {
            let mut differences: Vec<f64> = peptides
                .iter()
                .filter_map(|p| Some(p[j]? - p[k]?))
                .collect();
            if differences.is_empty() || differences.len() < min_shared {
                continue;
            }
            differences.sort_by(f64::total_cmp);
            let typical = summary(&differences);
            ratio[j][k] = Some(typical);
            ratio[k][j] = Some(-typical);
        }
    }
    ratio
}

/// The median of sorted values.
pub(crate) fn median(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// The half-sample mode of sorted values (Bickel & Frühwirth 2006): keep
/// the densest half (the run of ⌈n/2⌉ values with the smallest range) until
/// three or fewer remain. It finds the largest cluster even when the rest
/// lean to one side, where a median is pulled along.
pub(crate) fn half_sample_mode(sorted: &[f64]) -> f64 {
    let mut values = sorted;
    while values.len() > 3 {
        let half = values.len().div_ceil(2);
        let start = (0..=values.len() - half)
            .min_by(|&a, &b| {
                (values[a + half - 1] - values[a]).total_cmp(&(values[b + half - 1] - values[b]))
            })
            .unwrap_or(0);
        values = &values[start..start + half];
    }
    match values {
        [a, b, c] if b - a < c - b => (a + b) / 2.0,
        [a, b, c] if b - a > c - b => (b + c) / 2.0,
        [_, b, _] => *b,
        [a, b] => (a + b) / 2.0,
        [a] => *a,
        _ => 0.0,
    }
}

/// Samples with any peptide, grouped into connected components by shared
/// peptides.
fn components(ratio: &[Vec<Option<f64>>], summed: &[Option<f64>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; summed.len()];
    let mut found = Vec::new();
    for start in 0..summed.len() {
        if seen[start] || summed[start].is_none() {
            continue;
        }
        seen[start] = true;
        let mut component = vec![start];
        let mut next = 0;
        while next < component.len() {
            let j = component[next];
            next += 1;
            for k in 0..summed.len() {
                if !seen[k] && ratio[j][k].is_some() {
                    seen[k] = true;
                    component.push(k);
                }
            }
        }
        component.sort_unstable();
        found.push(component);
    }
    found
}

/// Least-squares log2 levels of a component, up to a constant: minimize
/// Σ (x_j − x_k − r_jk)² over its pairs. The normal equations are the graph
/// Laplacian; fixing the first sample at 0 makes them positive definite.
fn fit(component: &[usize], ratio: &[Vec<Option<f64>>]) -> Vec<f64> {
    let n = component.len();
    if n == 1 {
        return vec![0.0];
    }
    // Unknowns: x_1..x_{n-1} (x_0 = 0).
    let m = n - 1;
    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];
    for (row, &j) in component.iter().enumerate().skip(1) {
        for (col, &k) in component.iter().enumerate() {
            if let Some(r) = ratio[j][k] {
                a[row - 1][row - 1] += 1.0;
                b[row - 1] += r;
                if col > 0 {
                    a[row - 1][col - 1] -= 1.0;
                }
            }
        }
    }
    let solved = solve(a, b);
    std::iter::once(0.0).chain(solved).collect()
}

/// Gaussian elimination with partial pivoting (small, well-posed systems).
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))
            .unwrap_or(col);
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (above, below) = a.split_at_mut(col + 1);
        let pivot_row = &above[col];
        for (offset, row) in below.iter_mut().enumerate() {
            let factor = row[col] / pivot_row[col];
            for (value, &p) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *value -= factor * p;
            }
            b[col + 1 + offset] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let rest: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - rest) / a[row][row];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn log(values: &[Option<f64>]) -> Vec<Option<f64>> {
        values.iter().map(|v| v.map(f64::log2)).collect()
    }

    fn assert_close(actual: &[Option<f64>], expected: &[Option<f64>]) {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected) {
            match (a, e) {
                (Some(a), Some(e)) => assert!((a - e).abs() < 1e-9, "{a} vs {e}"),
                _ => assert_eq!(a, e),
            }
        }
    }

    #[test]
    fn levels_follow_the_ratios_over_shared_peptides_only() {
        // Shared peptides say B = 2 × A (log2 −1 on both); p3, caught only
        // by B, takes no part. Levels are centred on 0.
        let peptides = vec![
            log(&[Some(8.0), Some(16.0)]),
            log(&[Some(4.0), Some(8.0)]),
            log(&[None, Some(1000.0)]),
        ];
        assert_close(&relative_levels(&peptides, 2, 1), &[Some(-0.5), Some(0.5)]);
    }

    #[test]
    fn samples_connect_through_a_chain_of_shared_peptides() {
        // A–B share p1 (B = 2A), B–C share p2 (C = 4B); A and C share none.
        // Levels 0, 1, 3 before centring on their mean 4/3.
        let peptides = vec![
            log(&[Some(10.0), Some(20.0), None]),
            log(&[None, Some(5.0), Some(20.0)]),
        ];
        assert_close(
            &relative_levels(&peptides, 3, 1),
            &[Some(-4.0 / 3.0), Some(-1.0 / 3.0), Some(5.0 / 3.0)],
        );
    }

    #[test]
    fn pairs_sharing_too_few_peptides_do_not_connect_samples() {
        // A and B share two peptides, A and C one, D none. With at least two
        // shared peptides required only A and B are linked (B = 2A).
        let peptides = vec![
            log(&[Some(8.0), Some(16.0), None, None]),
            log(&[Some(4.0), Some(8.0), None, None]),
            log(&[Some(2.0), None, Some(2.0), None]),
            log(&[None, None, None, None]),
        ];
        assert_close(
            &relative_levels(&peptides, 4, 2),
            &[Some(-0.5), Some(0.5), None, None],
        );
    }

    #[test]
    fn the_half_sample_mode_finds_the_largest_cluster() {
        // Five values near 0 and four near 2: the median (0.1) leans toward
        // the 2s; the mode stays in the cluster at 0.
        let values = [-0.1, -0.05, 0.0, 0.05, 0.1, 1.9, 2.0, 2.1, 2.2];
        let mode = half_sample_mode(&values);
        assert!(mode.abs() <= 0.05, "{mode}");
        assert_eq!(half_sample_mode(&[1.0, 2.0, 2.5]), 2.25);
        assert_eq!(half_sample_mode(&[3.0]), 3.0);
    }
}
