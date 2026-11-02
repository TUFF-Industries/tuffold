//! Sample QC in the Quant tab: a PCA of the samples, their correlations, and
//! a line per sample.

use gpui_kit::{
    Bounds, Div, Hsla, ParentElement as _, Styled as _, canvas, div, fill, point, px, relative,
    size,
};
use tuffold_quant::qc::{MIN_GROUPS, MIN_PCA_SAMPLES, Pca, SampleQc};

use super::{group_digits, panel_shell, panel_title};
use crate::quant::QcShown;
use crate::theme::{MONO, Palette};

const PCA_HEIGHT: f32 = 300.;
const CELL_MAX: f32 = 34.;
/// More samples than these and the labels would overlap.
const PCA_LABELS: usize = 24;
const CELL_TEXT: usize = 12;
/// Past this the cells would be too thin to see.
const GRID_SAMPLES: usize = 100;
const LABEL_WIDTH: f32 = 96.;
const LABEL_HEIGHT: f32 = 14.;
/// Labels nearer than this share of the plot's width could overlap.
const LABEL_SHARE: f32 = 0.2;
/// Above this share of the colour scale, a cell takes light-on-strong text.
const STRONG_CELL: f64 = 0.6;

pub(super) fn qc_body(shown: &QcShown, palette: &Palette) -> Div {
    let qc = &shown.qc;
    div()
        .flex()
        .flex_col()
        .gap_3()
        .child(
            div()
                .text_size(px(12.))
                .text_color(palette.ink_muted)
                .child(format!(
                    "{} samples · {} protein groups · {}",
                    qc.samples.len(),
                    group_digits(qc.groups),
                    shown.values
                )),
        )
        .child(
            div()
                .flex()
                .gap_4()
                .child(
                    panel_shell(palette)
                        .flex_1()
                        .min_w_0()
                        .child(panel_title("PCA"))
                        .child(pca_view(qc, palette)),
                )
                .child(
                    panel_shell(palette)
                        .flex_1()
                        .min_w_0()
                        .child(panel_title("Correlation"))
                        .child(correlation_view(qc, palette)),
                ),
        )
        .child(
            panel_shell(palette)
                .child(panel_title("Samples"))
                .child(sample_rows(qc, palette)),
        )
}

/// A colour per category, in first-seen order; uncategorized samples are
/// the axis grey.
fn category_colors(qc: &SampleQc, palette: &Palette) -> Vec<(Option<String>, Hsla)> {
    let series = [
        palette.signal,
        palette.ion_y,
        palette.ok,
        palette.caution,
        palette.fault,
        palette.ink_muted,
    ];
    let mut seen: Vec<(Option<String>, Hsla)> = Vec::new();
    for sample in &qc.samples {
        if seen.iter().any(|(c, _)| *c == sample.category) {
            continue;
        }
        let color = match sample.category {
            None => palette.rule_strong,
            Some(_) => {
                let named = seen.iter().filter(|(c, _)| c.is_some()).count();
                series[named % series.len()]
            }
        };
        seen.push((sample.category.clone(), color));
    }
    seen
}

fn pca_view(qc: &SampleQc, palette: &Palette) -> Div {
    let Some(pca) = &qc.pca else {
        return div()
            .h(px(PCA_HEIGHT))
            .text_size(px(13.))
            .text_color(palette.ink_faint)
            .child(format!(
                "Needs {MIN_PCA_SAMPLES} samples and {MIN_GROUPS} groups measured in all of them, \
                 with some variation."
            ));
    };
    let colors = category_colors(qc, palette);
    let color_of = |category: &Option<String>| {
        colors
            .iter()
            .find(|(c, _)| c == category)
            .map_or(palette.rule_strong, |(_, color)| *color)
    };
    let (x_of, y_of) = scales(pca);
    let dots: Vec<(f32, f32, Hsla)> = pca
        .scores
        .iter()
        .zip(&qc.samples)
        .map(|(&(x, y), sample)| (x_of(x), y_of(y), color_of(&sample.category)))
        .collect();
    let tops = label_tops(&dots);
    let labels = (qc.samples.len() <= PCA_LABELS).then(|| {
        pca.scores
            .iter()
            .zip(&qc.samples)
            .zip(tops)
            .map(|((&(x, _), sample), top)| {
                div()
                    .absolute()
                    .left(relative(x_of(x)))
                    .top(px(top))
                    .ml(px(8.))
                    .font_family(MONO)
                    .text_size(px(10.))
                    .text_color(palette.ink_muted)
                    .child(sample.name.clone())
            })
    });
    let axis = palette.rule_strong;
    let caption = |text: String| {
        div()
            .font_family(MONO)
            .text_size(px(11.))
            .text_color(palette.ink_faint)
            .child(text)
    };
    let legend = colors
        .iter()
        .filter(|_| colors.len() > 1)
        .map(|(category, color)| {
            div()
                .flex()
                .items_center()
                .gap_1()
                .child(div().size(px(8.)).rounded_full().bg(*color))
                .child(caption(
                    category.clone().unwrap_or_else(|| "no category".into()),
                ))
        });
    div()
        .flex()
        .flex_col()
        .gap_1()
        .child(
            div()
                .relative()
                .w_full()
                .h(px(PCA_HEIGHT))
                .child(
                    canvas(
                        |_, _, _| {},
                        move |bounds: Bounds<_>, (), window, _| {
                            let (left, top) = (bounds.left(), bounds.top());
                            let (width, height) = (bounds.size.width, bounds.size.height);
                            window.paint_quad(fill(
                                Bounds::new(point(left, top + height), size(width, px(1.))),
                                axis,
                            ));
                            window.paint_quad(fill(
                                Bounds::new(point(left, top), size(px(1.), height)),
                                axis,
                            ));
                            for (dx, dy, color) in &dots {
                                let center = point(left + width * *dx, top + height * *dy);
                                let radius = px(4.);
                                window.paint_quad(
                                    fill(
                                        Bounds::new(
                                            point(center.x - radius, center.y - radius),
                                            size(radius * 2., radius * 2.),
                                        ),
                                        *color,
                                    )
                                    .corner_radii(radius),
                                );
                            }
                        },
                    )
                    .absolute()
                    .top_0()
                    .left_0()
                    .size_full(),
                )
                .children(labels.into_iter().flatten()),
        )
        .child(
            div()
                .flex()
                .flex_wrap()
                .items_center()
                .gap_3()
                .child(caption(format!(
                    "PC1 {:.1}% · PC2 {:.1}% · {} groups in every sample",
                    pca.explained[0] * 100.0,
                    pca.explained[1] * 100.0,
                    group_digits(pca.groups)
                )))
                .child(div().flex_1())
                .children(legend),
        )
}

/// Where each dot's label goes, in px from the plot's top: beside its dot,
/// moved down from any label above it that it would overprint. Replicates
/// sitting together are the usual case.
fn label_tops(dots: &[(f32, f32, Hsla)]) -> Vec<f32> {
    let mut order: Vec<usize> = (0..dots.len()).collect();
    order.sort_by(|&a, &b| dots[a].1.total_cmp(&dots[b].1));
    let mut tops = vec![0.; dots.len()];
    let mut placed: Vec<(f32, f32)> = Vec::new();
    for at in order {
        let (x, y, _) = dots[at];
        let mut top = PCA_HEIGHT * y - LABEL_HEIGHT / 2.;
        for &(other_x, other_top) in &placed {
            if (x - other_x).abs() < LABEL_SHARE && top < other_top + LABEL_HEIGHT {
                top = other_top + LABEL_HEIGHT;
            }
        }
        placed.push((x, top));
        tops[at] = top;
    }
    tops
}

/// Score to a share of the plot, leaving the labels room on the right. PC2
/// gets no more room than its size earns: stretched to fill the height, a
/// component with 1% of the variance would look as far apart as PC1. The plot
/// is about twice as wide as high, so half PC1's span fills the height at
/// equal units.
fn scales(pca: &Pca) -> (impl Fn(f64) -> f32 + use<>, impl Fn(f64) -> f32 + use<>) {
    let range = |values: Vec<f64>| {
        let low = values.iter().copied().fold(f64::INFINITY, f64::min);
        let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (low, (high - low).max(1e-9))
    };
    let (x_low, x_span) = range(pca.scores.iter().map(|p| p.0).collect());
    let (y_low, y_span) = range(pca.scores.iter().map(|p| p.1).collect());
    let y_center = y_low + y_span / 2.0;
    let y_span = y_span.max(x_span / 2.0);
    let y_low = y_center - y_span / 2.0;
    (
        move |x: f64| (0.05 + 0.75 * (x - x_low) / x_span) as f32,
        move |y: f64| (0.92 - 0.84 * (y - y_low) / y_span) as f32,
    )
}

fn correlation_view(qc: &SampleQc, palette: &Palette) -> Div {
    let n = qc.samples.len();
    if n > GRID_SAMPLES {
        return div()
            .text_size(px(13.))
            .text_color(palette.ink_faint)
            .child(format!(
                "{n} samples are too many to draw; Export TSV has every pair."
            ));
    }
    let others = || {
        qc.correlation.iter().enumerate().flat_map(|(a, row)| {
            row.iter()
                .enumerate()
                .filter(move |(b, _)| a != *b)
                .filter_map(|(_, r)| *r)
        })
    };
    let lowest = others().fold(1.0, f64::min);
    let highest = others().fold(-1.0, f64::max);
    let cell = (PCA_HEIGHT / n.max(1) as f32).min(CELL_MAX);
    let with_text = n <= CELL_TEXT;
    let span = (highest - lowest).max(1e-9);
    let label = |name: &str, width: f32| {
        div()
            .w(px(width))
            .flex_none()
            .truncate()
            .font_family(MONO)
            .text_size(px(10.))
            .text_color(palette.ink_muted)
            .child(name.to_string())
    };
    let columns = with_text.then(|| {
        div()
            .flex()
            .gap_px()
            .child(div().w(px(LABEL_WIDTH)).flex_none())
            .children(qc.samples.iter().map(|sample| {
                div()
                    .flex_1()
                    .min_w_0()
                    .truncate()
                    .text_center()
                    .font_family(MONO)
                    .text_size(px(10.))
                    .text_color(palette.ink_muted)
                    .child(sample.name.clone())
            }))
    });
    let rows = qc
        .correlation
        .iter()
        .enumerate()
        .zip(&qc.samples)
        .map(|((a, row), sample)| {
            let cells = row.iter().enumerate().map(move |(b, r)| {
                // A sample against itself says nothing, and the scale leaves it out.
                let r = if a == b { None } else { *r };
                let share = r.map(|r| ((r - lowest) / span).clamp(0.0, 1.0));
                let background = share.map_or(palette.paper_sunken, |t| palette.heat(t, 1.0));
                let text = if share.is_some_and(|t| t > STRONG_CELL) {
                    palette.on_signal
                } else {
                    palette.ink
                };
                let shown = match r {
                    Some(r) => format!("{r:.3}"),
                    None if a == b => String::new(),
                    None => "–".to_string(),
                };
                div()
                    .flex_1()
                    .min_w_0()
                    .h(px(cell - 1.))
                    .flex()
                    .items_center()
                    .justify_center()
                    .bg(background)
                    .font_family(MONO)
                    .text_size(px(10.))
                    .text_color(text)
                    .children(with_text.then_some(shown))
            });
            div()
                .flex()
                .items_center()
                .gap_px()
                .child(label(&sample.name, LABEL_WIDTH))
                .children(cells)
        });
    div()
        .flex()
        .flex_col()
        .gap_px()
        .children(columns)
        .children(rows)
        .children((n >= 2 && lowest <= highest).then(|| {
            div()
                .mt_2()
                .font_family(MONO)
                .text_size(px(11.))
                .text_color(palette.ink_faint)
                .child(format!("Pearson r, {lowest:.3} to {highest:.3}"))
        }))
}

fn sample_rows(qc: &SampleQc, palette: &Palette) -> Div {
    let line = |cells: [String; 5], color: Hsla| {
        let [name, category, measured, median, mean] = cells;
        let right = |text: String| div().w(px(96.)).flex_none().text_right().child(text);
        div()
            .flex()
            .items_center()
            .gap_3()
            .font_family(MONO)
            .text_size(px(12.))
            .text_color(color)
            .child(div().flex_1().min_w_0().truncate().child(name))
            .child(div().w(px(140.)).flex_none().truncate().child(category))
            .child(right(measured))
            .child(right(median))
            .child(right(mean))
    };
    let rows = qc.samples.iter().enumerate().map(|(s, sample)| {
        let others: Vec<f64> = qc.correlation[s]
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != s)
            .filter_map(|(_, r)| *r)
            .collect();
        let mean = (!others.is_empty()).then(|| others.iter().sum::<f64>() / others.len() as f64);
        line(
            [
                sample.name.clone(),
                sample.category.clone().unwrap_or_else(|| "–".into()),
                group_digits(sample.measured),
                sample.median_log2.map_or("–".into(), |v| format!("{v:.2}")),
                mean.map_or("–".into(), |v| format!("{v:.3}")),
            ],
            palette.ink,
        )
    });
    div()
        .flex()
        .flex_col()
        .gap_1()
        .child(line(
            [
                "Sample".into(),
                "Category".into(),
                "Groups".into(),
                "Median log2".into(),
                "Mean r".into(),
            ],
            palette.ink_faint,
        ))
        .children(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_of_samples_together_do_not_overprint() {
        let color = Hsla::default();
        let dots = [(0.8, 0.5, color), (0.81, 0.5, color), (0.1, 0.5, color)];
        let tops = label_tops(&dots);
        assert!((tops[0] - tops[1]).abs() >= LABEL_HEIGHT);
        assert_eq!(
            tops[2],
            PCA_HEIGHT * 0.5 - LABEL_HEIGHT / 2.,
            "far away, stays"
        );
    }
}
