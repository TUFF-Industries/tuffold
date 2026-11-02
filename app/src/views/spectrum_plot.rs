//! The annotated MS/MS spectrum: sticks painted through `canvas` (b ions in
//! one color, y in another, unexplained peaks in rule gray), ion labels that
//! never overlap, an m/z axis, and the fragment ladder above.

use gpui_kit::{
    Bounds, Div, FontWeight, Hsla, IntoElement as _, ParentElement, Styled, canvas, div, fill,
    point, px, relative, size,
};
use tuffold_core::fragments::IonSeries;

use crate::browse::{CoverageLine, LadderStep, SpectrumPlot};
use crate::theme::{MONO, Palette};

const PLOT_HEIGHT: f32 = 230.;
/// Labels closer than this (fraction of the width) to a taller one are hidden.
const LABEL_GAP: f32 = 0.028;
const MAX_LABELS: usize = 40;

fn ion_color(series: IonSeries, palette: &Palette) -> Hsla {
    match series {
        IonSeries::B => palette.ion_b,
        IonSeries::Y => palette.ion_y,
    }
}

/// Indices of annotated peaks to label: tallest first, skipping any too
/// close to one already labelled.
pub fn label_choice(plot: &SpectrumPlot) -> Vec<usize> {
    let span = (plot.max_mz - plot.min_mz).max(1.0);
    let mut annotated: Vec<usize> = (0..plot.peaks.len())
        .filter(|&i| plot.peaks[i].ion.is_some())
        .collect();
    annotated.sort_by(|&a, &b| plot.peaks[b].relative.total_cmp(&plot.peaks[a].relative));
    let mut chosen: Vec<usize> = Vec::new();
    for index in annotated {
        let x = ((plot.peaks[index].mz - plot.min_mz) / span) as f32;
        let clear = chosen.iter().all(|&other| {
            (((plot.peaks[other].mz - plot.min_mz) / span) as f32 - x).abs() > LABEL_GAP
        });
        if clear {
            chosen.push(index);
        }
        if chosen.len() == MAX_LABELS {
            break;
        }
    }
    chosen
}

/// Round m/z steps for about six axis ticks.
fn ticks(min: f64, max: f64) -> Vec<f64> {
    let span = (max - min).max(1.0);
    let raw = span / 6.0;
    let magnitude = 10_f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|factor| factor * magnitude)
        .find(|step| *step >= raw)
        .unwrap_or(10.0 * magnitude);
    let first = (min / step).ceil() * step;
    (0..)
        .map(|i| first + f64::from(i) * step)
        .take_while(|tick| *tick <= max)
        .collect()
}

pub fn spectrum_plot_view(plot: &SpectrumPlot, palette: &Palette) -> Div {
    let span = (plot.max_mz - plot.min_mz).max(1.0);
    let x_of = move |mz: f64| ((mz - plot.min_mz) / span) as f32;
    let sticks: Vec<(f32, f32, Hsla, bool)> = plot
        .peaks
        .iter()
        .map(|peak| {
            let color = peak
                .ion
                .as_ref()
                .map_or(palette.rule_strong, |(series, _)| {
                    ion_color(*series, palette)
                });
            (x_of(peak.mz), peak.relative, color, peak.ion.is_some())
        })
        .collect();
    let axis = palette.rule_strong;
    let labels = label_choice(plot).into_iter().map(|index| {
        let peak = &plot.peaks[index];
        let (series, label) = peak.ion.as_ref().expect("labelled peaks are annotated");
        let top = 1.0 - 0.86 * peak.relative;
        div()
            .absolute()
            .left(relative(x_of(peak.mz)))
            .top(px(PLOT_HEIGHT * top - 16.))
            .ml(px(-10.))
            .font_family(MONO)
            .text_size(px(11.))
            .font_weight(FontWeight::MEDIUM)
            .text_color(ion_color(*series, palette))
            .child(label.clone())
    });
    let tick_labels = ticks(plot.min_mz, plot.max_mz).into_iter().map(|tick| {
        div()
            .absolute()
            .left(relative(x_of(tick)))
            .top(px(PLOT_HEIGHT + 4.))
            .ml(px(-14.))
            .font_family(MONO)
            .text_size(px(11.))
            .text_color(palette.ink_faint)
            .child(format!("{tick:.0}"))
    });
    div()
        .relative()
        .w_full()
        .h(px(PLOT_HEIGHT + 22.))
        .child(
            canvas(
                |_, _, _| {},
                move |bounds: Bounds<_>, (), window, _| {
                    let height = px(PLOT_HEIGHT) * 0.86;
                    let floor = bounds.top() + px(PLOT_HEIGHT);
                    window.paint_quad(fill(
                        Bounds::new(point(bounds.left(), floor), size(bounds.size.width, px(1.))),
                        axis,
                    ));
                    // Unexplained peaks first, so annotated ones paint on top.
                    for annotated in [false, true] {
                        for (x, relative, color, is_annotated) in &sticks {
                            if *is_annotated != annotated {
                                continue;
                            }
                            let width = if annotated { px(2.) } else { px(1.) };
                            let stick = height * *relative;
                            let left = bounds.left() + bounds.size.width * *x;
                            window.paint_quad(fill(
                                Bounds::new(point(left, floor - stick), size(width, stick)),
                                *color,
                            ));
                        }
                    }
                },
            )
            .absolute()
            .top_0()
            .left_0()
            .w_full()
            .h(px(PLOT_HEIGHT)),
        )
        .children(labels)
        .children(tick_labels)
}

/// The peptide with b ions marked above the cleavage and y ions below.
pub fn ladder_view(steps: &[LadderStep], palette: &Palette) -> Div {
    let cells = steps.iter().enumerate().flat_map(|(index, step)| {
        let mut residue = div()
            .w(px(18.))
            .flex()
            .justify_center()
            .rounded(px(3.))
            .font_family(MONO)
            .text_size(px(16.))
            .font_weight(if step.localized {
                FontWeight::BOLD
            } else if step.modified {
                FontWeight::SEMIBOLD
            } else {
                FontWeight::NORMAL
            })
            .text_color(if step.modified {
                palette.signal
            } else {
                palette.ink
            })
            .child(step.residue.to_string());
        if step.localized {
            residue = residue
                .bg(palette.signal_soft)
                .border_1()
                .border_color(palette.signal);
        }
        let last = index + 1 == steps.len();
        let mark = |on: bool, color: Hsla| {
            div()
                .w(px(2.))
                .h(px(10.))
                .bg(if on { color } else { palette.rule })
        };
        let cleavage = (!last).then(|| {
            div()
                .w(px(8.))
                .flex()
                .flex_col()
                .items_center()
                .justify_between()
                .h(px(28.))
                .child(mark(step.b, palette.ion_b))
                .child(mark(step.y, palette.ion_y))
        });
        std::iter::once(residue.into_any_element()).chain(cleavage.map(|c| c.into_any_element()))
    });
    div().flex().items_center().children(cells)
}

/// The leader's sequence with covered residues highlighted.
pub fn coverage_view(lines: &[CoverageLine], palette: &Palette) -> Div {
    let rows = lines.iter().map(|line| {
        // Runs of equal coverage render as one span each.
        let mut runs: Vec<(bool, String)> = Vec::new();
        for (residue, covered) in &line.residues {
            match runs.last_mut() {
                Some((state, text)) if state == covered => text.push(*residue),
                _ => runs.push((*covered, residue.to_string())),
            }
        }
        div()
            .flex()
            .gap_3()
            .child(
                div()
                    .w(px(34.))
                    .flex()
                    .justify_end()
                    .text_color(palette.ink_faint)
                    .child(line.start.to_string()),
            )
            .child(
                div()
                    .flex()
                    .children(runs.into_iter().map(|(covered, text)| {
                        if covered {
                            div()
                                .bg(palette.signal_soft)
                                .text_color(palette.ink)
                                .font_weight(FontWeight::MEDIUM)
                                .child(text)
                        } else {
                            div().text_color(palette.ink_faint).child(text)
                        }
                    })),
            )
    });
    div()
        .flex()
        .flex_col()
        .gap_0p5()
        .font_family(MONO)
        .text_size(px(12.))
        .children(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::browse::PlotPeak;

    fn plot(peaks: Vec<(f64, f32, bool)>) -> SpectrumPlot {
        SpectrumPlot {
            peaks: peaks
                .into_iter()
                .map(|(mz, relative, annotated)| PlotPeak {
                    mz,
                    relative,
                    ion: annotated.then(|| (IonSeries::B, "b2".into())),
                })
                .collect(),
            min_mz: 0.0,
            max_mz: 1000.0,
            ladder: Vec::new(),
            matched: 0,
            theoretical: 0,
            explained: 0.0,
        }
    }

    #[test]
    fn labels_go_to_the_tallest_annotated_peaks_without_overlap() {
        // 500 (0.9) and 510 (1.0) are 1% apart: only the taller keeps a label.
        // 800 is far enough; 300 is not annotated.
        let chosen = label_choice(&plot(vec![
            (300.0, 1.0, false),
            (500.0, 0.9, true),
            (510.0, 1.0, true),
            (800.0, 0.2, true),
        ]));
        assert_eq!(chosen, [2, 3]);
    }

    #[test]
    fn ticks_land_on_round_values() {
        assert_eq!(
            ticks(0.0, 1050.0),
            [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
        );
        assert_eq!(
            ticks(0.0, 300.0),
            [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
        );
    }
}
