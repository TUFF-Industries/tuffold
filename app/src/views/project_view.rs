use gpui_kit::{Div, FontWeight, Hsla, IntoElement, ParentElement, Styled, div, px, relative};

use super::validation_panel::{ValidationControls, validation_panel};
use super::{callout, eyebrow, group_digits, panel, panel_shell, panel_title};
use crate::project::{Finding, FindingLevel, ProjectReport, StatusKind};
use crate::theme::{MONO, Palette};

pub fn project_view(
    report: &ProjectReport,
    controls: &ValidationControls,
    palette: &Palette,
) -> impl IntoElement {
    let has_data = report.summary.psms > 0;
    div()
        .flex()
        .flex_col()
        .w_full()
        .max_w(px(1216.))
        .mx_auto()
        .px(px(48.))
        .pt_8()
        .pb(px(64.))
        .child(header(report, palette))
        .children(has_data.then(|| stat_strip(report, palette)))
        .children(validation_panel(&report.validation, controls, palette))
        .child(
            div()
                .mt_6()
                .flex()
                .items_start()
                .gap_6()
                .child(
                    div()
                        .flex()
                        .flex_col()
                        .gap_6()
                        .flex_basis(relative(0.5))
                        .min_w_0()
                        .child(imports_panel(report, palette))
                        .children(has_data.then(|| runs_panel(report, palette))),
                )
                .child(
                    div()
                        .flex_1()
                        .min_w_0()
                        .child(findings_panel(&report.findings, palette)),
                ),
        )
}

/// The project's state, one line (its name is in the title bar).
fn header(report: &ProjectReport, palette: &Palette) -> impl IntoElement {
    let status = report.status();
    let color = match status.kind {
        StatusKind::Ready => palette.ok,
        StatusKind::Limited => palette.caution,
        StatusKind::Invalid => palette.fault,
    };
    div()
        .flex()
        .items_center()
        .gap_2()
        .pb_4()
        .text_size(px(12.))
        .text_color(color)
        .child(div().size(px(6.)).rounded_full().bg(color))
        .child(status.label)
}

/// What's in the project, as a line: the validation below is the headline.
fn stat_strip(report: &ProjectReport, palette: &Palette) -> impl IntoElement {
    let summary = &report.summary;
    let count = group_digits;
    let line = [
        format!(
            "{} PSMs ({} decoy)",
            count(summary.psms),
            count(summary.decoy_psms)
        ),
        format!("{} peptides", count(summary.peptides)),
        format!(
            "{} proteins ({} decoy, {} with sequences)",
            count(summary.proteins),
            count(summary.decoy_proteins),
            count(summary.proteins_with_sequence)
        ),
        format!("{} spectra", count(summary.spectra)),
        format!("{} MS runs", count(summary.ms_runs)),
    ]
    .join(" · ");
    div()
        .pb_6()
        .font_family(MONO)
        .text_size(px(12.))
        .text_color(palette.ink_muted)
        .child(line)
}

fn imports_panel(report: &ProjectReport, palette: &Palette) -> Div {
    if report.imports.is_empty() {
        return panel("Imports", palette).child(
            div()
                .text_size(px(13.))
                .text_color(palette.ink_muted)
                .child("Nothing imported yet. Add search results with Add data…"),
        );
    }
    let rows = report.imports.iter().map(|record| {
        let file = record
            .source_path
            .rsplit(['/', '\\'])
            .next()
            .unwrap_or(&record.source_path)
            .to_string();
        div()
            .flex()
            .items_center()
            .gap_3()
            .py_2()
            .border_b_1()
            .border_color(palette.rule)
            .text_size(px(13.))
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .font_family(MONO)
                    .truncate()
                    .child(file),
            )
            .child(
                div()
                    .px_2()
                    .rounded(px(3.))
                    .bg(palette.paper_sunken)
                    .font_family(MONO)
                    .text_size(px(11.))
                    .text_color(palette.ink_muted)
                    .child(record.format.clone()),
            )
            .child(
                div()
                    .w(px(140.))
                    .truncate()
                    .text_color(palette.ink_muted)
                    .child(record.software.clone()),
            )
            .child(
                div()
                    .w(px(64.))
                    .flex()
                    .justify_end()
                    .font_family(MONO)
                    .child(group_digits(record.psms)),
            )
    });
    panel("Imports", palette).child(div().flex().flex_col().children(rows))
}

fn runs_panel(report: &ProjectReport, palette: &Palette) -> impl IntoElement {
    let runs = &report.summary.runs;
    let max = runs.iter().map(|run| run.psms).max().unwrap_or(0).max(1);
    let bars = runs.iter().map(|run| {
        let label_color = if run.psms == 0 {
            palette.caution
        } else {
            palette.ink
        };
        div()
            .flex()
            .items_center()
            .gap_3()
            .text_size(px(13.))
            .child(
                div()
                    .w(px(176.))
                    .flex()
                    .flex_col()
                    .min_w_0()
                    .child(
                        div()
                            .font_family(MONO)
                            .text_color(label_color)
                            .truncate()
                            .child(run.name.clone()),
                    )
                    .children((run.biosample != run.name).then(|| {
                        div()
                            .text_size(px(11.))
                            .text_color(palette.ink_faint)
                            .truncate()
                            .child(format!("biosample {}", run.biosample))
                    })),
            )
            .child(
                div()
                    .flex_1()
                    .h(px(8.))
                    .rounded_full()
                    .bg(palette.paper_sunken)
                    .child(
                        div()
                            .h_full()
                            .rounded_full()
                            .bg(palette.signal)
                            .w(relative(run.psms as f32 / max as f32)),
                    ),
            )
            .child(
                div()
                    .w(px(56.))
                    .flex()
                    .justify_end()
                    .font_family(MONO)
                    .text_color(label_color)
                    .child(group_digits(run.psms)),
            )
    });
    panel("PSMs per MS run", palette).child(div().flex().flex_col().gap_3().children(bars))
}

fn findings_panel(findings: &[Finding], palette: &Palette) -> impl IntoElement {
    let count = div()
        .px_2()
        .rounded(px(3.))
        .bg(palette.paper_sunken)
        .font_family(MONO)
        .text_size(px(11.))
        .text_color(palette.ink_muted)
        .child(findings.len().to_string());
    let group = |level: FindingLevel, title: &'static str, color: Hsla, background: Hsla| {
        let items: Vec<&Finding> = findings.iter().filter(|f| f.level == level).collect();
        (!items.is_empty()).then(|| {
            div()
                .flex()
                .flex_col()
                .gap_2()
                .child(eyebrow(title, palette))
                .children(items.into_iter().map(move |finding| {
                    let heading = match &finding.code {
                        Some(code) => format!("{} · {code}", finding.source),
                        None => finding.source.clone(),
                    };
                    callout(color, background)
                        .flex()
                        .flex_col()
                        .gap_1()
                        .child(
                            div()
                                .font_family(MONO)
                                .text_size(px(11.))
                                .font_weight(FontWeight::MEDIUM)
                                .text_color(color)
                                .child(heading),
                        )
                        .child(div().text_size(px(13.)).child(finding.message.clone()))
                }))
        })
    };
    panel_shell(palette)
        .child(
            div()
                .flex()
                .items_center()
                .gap_2()
                .child(panel_title("Findings"))
                .child(count),
        )
        .children(findings.is_empty().then(|| {
            div()
                .text_size(px(13.))
                .text_color(palette.ok)
                .child("Nothing to report for the last import.")
        }))
        .children(group(
            FindingLevel::Error,
            "Not imported",
            palette.fault,
            palette.fault_soft,
        ))
        .children(group(
            FindingLevel::Warning,
            "Imported with a limitation",
            palette.caution,
            palette.caution_soft,
        ))
}
