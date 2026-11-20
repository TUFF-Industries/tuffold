use gpui_kit::{Div, FontWeight, Hsla, IntoElement, ParentElement, Styled, div, px, relative};

use super::validation_panel::{ValidationControls, validation_panel};
use super::{
    callout, fact, facts, group_digits, section, section_header, section_shell, section_title,
    subsection_title,
};
use crate::project::{Finding, FindingLevel, ProjectReport, StatusKind};
use crate::theme::Palette;

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
                .mt(px(48.))
                .flex()
                .items_start()
                .gap(px(48.))
                .child(
                    div()
                        .flex()
                        .flex_col()
                        .gap(px(40.))
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
        .pb_3()
        .text_size(px(13.))
        .font_weight(FontWeight::MEDIUM)
        .text_color(color)
        .child(div().size(px(6.)).rounded_full().bg(color))
        .child(status.label)
}

/// What's in the project, as a row of counts: the validation below is the
/// headline.
fn stat_strip(report: &ProjectReport, palette: &Palette) -> impl IntoElement {
    let summary = &report.summary;
    let count = group_digits;
    let of_which = |text: String| {
        div()
            .text_color(palette.ink_faint)
            .child(format!("({text})"))
    };
    facts()
        .pb(px(40.))
        .text_size(px(13.))
        .child(
            fact(count(summary.psms), "PSMs", palette)
                .child(of_which(format!("{} decoy", count(summary.decoy_psms)))),
        )
        .child(fact(count(summary.peptides), "peptides", palette))
        .child(
            fact(count(summary.proteins), "proteins", palette).child(of_which(format!(
                "{} decoy, {} with sequences",
                count(summary.decoy_proteins),
                count(summary.proteins_with_sequence)
            ))),
        )
        .child(fact(count(summary.spectra), "spectra", palette))
        .child(fact(count(summary.ms_runs), "MS runs", palette))
}

fn imports_panel(report: &ProjectReport, palette: &Palette) -> Div {
    if report.imports.is_empty() {
        return section("Imports", palette).child(
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
                    .flex()
                    .flex_col()
                    .child(div().font_weight(FontWeight::MEDIUM).truncate().child(file))
                    .child(
                        div()
                            .text_size(px(12.))
                            .text_color(palette.ink_faint)
                            .truncate()
                            .child(format!("{}, {}", record.software, record.format)),
                    ),
            )
            .child(
                div()
                    .w(px(96.))
                    .flex()
                    .justify_end()
                    .gap_1()
                    .child(group_digits(record.psms))
                    .child(div().text_color(palette.ink_faint).child("PSMs")),
            )
    });
    section("Imports", palette).child(div().flex().flex_col().children(rows))
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
                    .h(px(6.))
                    .border_l_1()
                    .border_color(palette.rule_strong)
                    .child(
                        div()
                            .h_full()
                            .bg(palette.stain)
                            .w(relative(run.psms as f32 / max as f32)),
                    ),
            )
            .child(
                div()
                    .w(px(56.))
                    .flex()
                    .justify_end()
                    .text_color(label_color)
                    .child(group_digits(run.psms)),
            )
    });
    section("PSMs per MS run", palette).child(div().flex().flex_col().gap_3().children(bars))
}

fn findings_panel(findings: &[Finding], palette: &Palette) -> impl IntoElement {
    let group = |level: FindingLevel, title: &'static str, color: Hsla, background: Hsla| {
        let items: Vec<&Finding> = findings.iter().filter(|f| f.level == level).collect();
        (!items.is_empty()).then(|| {
            div()
                .flex()
                .flex_col()
                .gap_2()
                .child(subsection_title(title, palette))
                .children(items.into_iter().map(move |finding| {
                    callout(color, background)
                        .flex()
                        .flex_col()
                        .gap_1()
                        .child(
                            div()
                                .flex()
                                .items_baseline()
                                .gap_2()
                                .text_size(px(13.))
                                .child(
                                    div()
                                        .font_weight(FontWeight::SEMIBOLD)
                                        .text_color(color)
                                        .child(finding.source.clone()),
                                )
                                .children(finding.code.clone().map(|code| {
                                    div()
                                        .text_size(px(12.))
                                        .text_color(palette.ink_muted)
                                        .child(code)
                                })),
                        )
                        .child(
                            div()
                                .text_size(px(13.))
                                .line_height(px(20.))
                                .child(finding.message.clone()),
                        )
                }))
        })
    };
    let count = (!findings.is_empty()).then(|| findings.len().to_string().into());
    section_shell()
        .child(section_header(palette).child(section_title("Findings", count, palette)))
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
