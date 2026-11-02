//! The validation panel: what passes at the current thresholds, one-click
//! threshold changes, how each engine was ranked, and what to know.

use std::rc::Rc;

use gpui_kit::component::Selectable as _;
use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::component::{Disableable as _, Sizable as _};
use gpui_kit::{App, Div, FontWeight, ParentElement, Styled, Window, div, px};
use tuffold_store::ValidationSummary;
use tuffold_validate::{NoteCode, ThresholdMode, Thresholds, format_rate};

use super::{callout, eyebrow, group_digits, panel_shell, panel_title, track};
use crate::project::{
    Level, LevelView, MIN_PEPTIDE_CHOICES, ValidationState, choices, level_views, switched,
};
use crate::theme::{MONO, Palette};

/// Called with the thresholds the user picked.
pub type ThresholdHandler = Rc<dyn Fn(Thresholds, &mut Window, &mut App)>;
/// Called when the user asks to validate.
pub type ValidateHandler = Rc<dyn Fn(&mut Window, &mut App)>;

/// What the panel may do; `None` handlers render read-only controls.
#[derive(Clone)]
pub struct ValidationControls {
    pub on_thresholds: Option<ThresholdHandler>,
    pub on_validate: Option<ValidateHandler>,
    /// A job is running; controls wait for it.
    pub busy: bool,
}

pub fn validation_panel(
    state: &ValidationState,
    controls: &ValidationControls,
    palette: &Palette,
) -> Option<Div> {
    match state {
        ValidationState::NoData => None,
        ValidationState::Pending { stale } => Some(pending(*stale, controls, palette)),
        ValidationState::Failed(message) => Some(
            panel_shell(palette).child(panel_title("Validation")).child(
                callout(palette.fault, palette.fault_soft)
                    .text_size(px(13.))
                    .text_color(palette.fault)
                    .child(format!(
                        "Validation could not run: {message}. The imported data is saved."
                    )),
            ),
        ),
        ValidationState::Done(summary) => Some(done(summary, controls, palette)),
    }
}

fn pending(stale: bool, controls: &ValidationControls, palette: &Palette) -> Div {
    let message = if stale {
        "Data was imported since this project was last validated, so its results are out of \
         date. Validate to rank PSMs, estimate FDR, and group proteins again."
    } else {
        "This project has not been validated yet. Validate to rank PSMs, estimate FDR, and \
         group proteins."
    };
    let button = Button::new("validate-now")
        .primary()
        .small()
        .label("Validate now")
        .disabled(controls.busy || controls.on_validate.is_none());
    let button = match &controls.on_validate {
        Some(handler) => {
            let handler = Rc::clone(handler);
            button.on_click(move |_, window, cx| handler(window, cx))
        }
        None => button,
    };
    panel_shell(palette).child(panel_title("Validation")).child(
        div()
            .flex()
            .items_center()
            .gap_4()
            .child(
                div()
                    .flex_1()
                    .text_size(px(13.))
                    .text_color(palette.ink_muted)
                    .child(message),
            )
            .child(button),
    )
}

fn done(summary: &ValidationSummary, controls: &ValidationControls, palette: &Palette) -> Div {
    // Only the unusual setting is worth a word.
    let rescoring = (!summary.settings.rescore).then_some("rescoring off");
    let thresholds = summary.settings.thresholds;
    // Threshold by FDR (q-values) or by probability (1 − PEP).
    let modes = [
        (ThresholdMode::Fdr, "FDR"),
        (ThresholdMode::Probability, "Probability"),
    ]
    .map(|(mode, label)| {
        chip(
            format!("mode-{}", mode.as_str()),
            label.to_string(),
            thresholds.mode == mode,
            switched(thresholds, mode),
            controls,
        )
    });
    let tiles = level_views(summary)
        .into_iter()
        .enumerate()
        .map(|(index, view)| {
            let tile = level_tile(&view, summary, controls, palette);
            if index == 0 {
                tile
            } else {
                tile.pl_6().border_l_1().border_color(palette.rule)
            }
        });
    panel_shell(palette)
        .child(
            div()
                .flex()
                .items_center()
                .justify_between()
                .child(panel_title("Validation"))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap_3()
                        .children(rescoring.map(|rescoring| {
                            div()
                                .text_size(px(12.))
                                .text_color(palette.caution)
                                .child(rescoring)
                        }))
                        .child(track(palette).children(modes)),
                ),
        )
        .child(div().flex().items_start().children(tiles))
        .children(
            (summary.localization.scored_psms > 0)
                .then(|| localization_row(&summary.localization, palette)),
        )
        // Which fragment tolerance localization used: a fact, not a caution.
        .children(
            summary
                .notes
                .iter()
                .filter(|note| {
                    note.code == NoteCode::Localization && summary.localization.scored_psms > 0
                })
                .map(|note| {
                    div()
                        .text_size(px(11.))
                        .text_color(palette.ink_faint)
                        .child(note.message.clone())
                }),
        )
        .children((!summary.engines.is_empty()).then(|| engines(summary, palette)))
        // Engine notes repeat what each engine row already says; the
        // localization note is shown under its row.
        .children(
            summary
                .notes
                .iter()
                .filter(|note| !is_engine_note(note.code) && note.code != NoteCode::Localization)
                .map(|note| {
                    callout(palette.caution, palette.caution_soft)
                        .flex()
                        .flex_col()
                        .gap_1()
                        .child(
                            div()
                                .font_family(MONO)
                                .text_size(px(11.))
                                .font_weight(FontWeight::MEDIUM)
                                .text_color(palette.caution)
                                .child(note.code.as_str()),
                        )
                        .child(div().text_size(px(13.)).child(note.message.clone()))
                }),
        )
}

fn level_tile(
    view: &LevelView,
    summary: &ValidationSummary,
    controls: &ValidationControls,
    palette: &Palette,
) -> Div {
    let thresholds = summary.settings.thresholds;
    let current = view.level.threshold(&thresholds);
    let fdr_choices = choices(thresholds.mode).into_iter().map(|choice| {
        chip(
            format!("{}-{:?}-{choice}", thresholds.mode.as_str(), view.level),
            format_rate(choice),
            choice == current,
            view.level.with(&thresholds, choice),
            controls,
        )
    });
    let peptide_choices = (view.level == Level::ProteinGroups)
        .then(|| {
            MIN_PEPTIDE_CHOICES.into_iter().map(move |choice| {
                chip(
                    format!("min-peptides-{choice}"),
                    format!("≥{choice} pep"),
                    choice == thresholds.min_peptides,
                    Thresholds {
                        min_peptides: choice,
                        ..thresholds
                    },
                    controls,
                )
            })
        })
        .into_iter()
        .flatten();
    div()
        .flex_1()
        .flex()
        .flex_col()
        .pr_4()
        .min_w_0()
        .child(
            div()
                .text_size(px(13.))
                .text_color(palette.ink_muted)
                .child(view.level.label()),
        )
        .child(
            div()
                .mt_1()
                .font_family(MONO)
                .text_size(px(30.))
                .line_height(px(36.))
                .font_weight(FontWeight::MEDIUM)
                .text_color(if view.targets > 0 {
                    palette.ink
                } else {
                    palette.ink_faint
                })
                .child(group_digits(view.targets)),
        )
        .child(
            div()
                .mt_1()
                .text_size(px(12.))
                .text_color(palette.ink_faint)
                .child(view.caption.clone()),
        )
        .child(
            div()
                .mt_3()
                .flex()
                .flex_wrap()
                .gap_2()
                .child(track(palette).children(fdr_choices))
                .children(
                    (view.level == Level::ProteinGroups)
                        .then(|| track(palette).children(peptide_choices)),
                ),
        )
        .children(view.hint.clone().map(|hint| {
            div()
                .mt_2()
                .text_size(px(12.))
                .text_color(palette.caution)
                .child(hint)
        }))
}

fn chip(
    id: String,
    label: String,
    selected: bool,
    thresholds: Thresholds,
    controls: &ValidationControls,
) -> Button {
    let button = Button::new(id).xsmall().label(label);
    let button = button.ghost().selected(selected);
    match &controls.on_thresholds {
        Some(handler) if !selected => {
            let handler = Rc::clone(handler);
            button
                .disabled(controls.busy)
                .on_click(move |_, window, cx| handler(thresholds, window, cx))
        }
        _ => button,
    }
}

fn engines(summary: &ValidationSummary, palette: &Palette) -> Div {
    let rows = summary.engines.iter().map(|engine| {
        let ranked = engine.method != "unscored";
        let color = match engine.method.as_str() {
            "rescored" | "best_score" => palette.ink_muted,
            _ => palette.caution,
        };
        div()
            .flex()
            .items_start()
            .gap_3()
            .py_2()
            .border_b_1()
            .border_color(palette.rule)
            .text_size(px(13.))
            .child(
                div()
                    .w(px(176.))
                    .flex_none()
                    .flex()
                    .flex_col()
                    .child(
                        div()
                            .font_family(MONO)
                            .truncate()
                            .child(engine.name.clone()),
                    )
                    .child(
                        div()
                            .text_size(px(11.))
                            .text_color(palette.ink_faint)
                            .child(format!(
                                "{} PSMs · {} {}",
                                group_digits(engine.psms),
                                group_digits(engine.decoys),
                                if engine.decoys == 1 {
                                    "decoy"
                                } else {
                                    "decoys"
                                }
                            )),
                    ),
            )
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .text_color(color)
                    .child(engine.detail.clone()),
            )
            .child(
                div()
                    .w(px(96.))
                    .flex_none()
                    .flex()
                    .justify_end()
                    .font_family(MONO)
                    .text_color(if ranked {
                        palette.ink
                    } else {
                        palette.ink_faint
                    })
                    .child(if ranked {
                        format!("{} @1%", group_digits(engine.accepted))
                    } else {
                        "—".to_string()
                    }),
            )
    });
    div()
        .mt_2()
        .flex()
        .flex_col()
        .child(eyebrow("Engines", palette))
        .children(rows)
}

fn localization_row(loc: &tuffold_store::LocalizationSummary, palette: &Palette) -> Div {
    let flr = loc
        .estimated_flr
        .map_or_else(|| "—".into(), |f| format!("{:.2}%", f * 100.0));
    let mut detail = format!(
        "{} of {} ambiguous-PTM spectra localized at ≥{:.0}% site probability · est. FLR {flr}",
        group_digits(loc.localized_psms),
        group_digits(loc.scored_psms),
        loc.probability_threshold * 100.0,
    );
    if loc.missing_peaks_psms > 0 {
        detail.push_str(&format!(
            " · {} missing peaks",
            group_digits(loc.missing_peaks_psms)
        ));
    }
    if loc.capped_psms > 0 {
        detail.push_str(&format!(
            " · {} capped at 256 isoforms",
            group_digits(loc.capped_psms)
        ));
    }
    if loc.unscorable_psms > 0 {
        detail.push_str(&format!(
            " · {} unscorable",
            group_digits(loc.unscorable_psms)
        ));
    }
    div()
        .mt_1()
        .pt_2()
        .border_t_1()
        .border_color(palette.rule)
        .flex()
        .items_center()
        .justify_between()
        .gap_3()
        .text_size(px(12.))
        .child(eyebrow("PTM localization", palette))
        .child(
            div()
                .flex_1()
                .min_w_0()
                .truncate()
                .text_color(palette.ink_muted)
                .child(detail),
        )
}

fn is_engine_note(code: NoteCode) -> bool {
    matches!(code, NoteCode::EngineUnscored | NoteCode::RescoringFallback)
}
