mod browser;
mod empty_state;
mod models_offer;
mod project_view;
mod qc_view;
mod quant_view;
mod spectrum_plot;
mod tables;
mod validation_panel;

pub use browser::{Browser, BrowserMode, ShowProtein};
pub use empty_state::{EmptyStateProps, empty_state};
pub use models_offer::{ModelsOfferProps, models_offer};
pub use project_view::project_view;
pub use quant_view::{OpenGroup, QuantView};
pub use validation_panel::{ThresholdHandler, ValidateHandler, ValidationControls};

use gpui_kit::{Div, FontWeight, IntoElement, ParentElement, SharedString, Styled, div, px};

use crate::theme::Palette;

/// The sunken track a set of exclusive choices sits in, like the tabs.
pub(crate) fn track(palette: &Palette) -> Div {
    div()
        .flex()
        .items_center()
        .p(px(2.))
        .gap(px(2.))
        .rounded(px(6.))
        .bg(palette.paper_sunken)
}

/// A section of a page: no box. Lists of facts sit on the page itself, and
/// their title row carries the rule (see `section_header`).
pub(crate) fn section_shell() -> Div {
    div().flex().flex_col().gap_3().min_w_0()
}

/// A section's title row, ruled underneath; callers add the title and any
/// controls (spread to the right).
pub(crate) fn section_header(palette: &Palette) -> Div {
    div()
        .flex()
        .items_center()
        .justify_between()
        .gap_3()
        .min_h(px(28.))
        .pb_2()
        .border_b_1()
        .border_color(palette.rule_strong)
}

/// A title, with what it counts beside it in a lighter ink (no pill).
pub(crate) fn section_title(
    title: impl Into<SharedString>,
    count: Option<SharedString>,
    palette: &Palette,
) -> Div {
    div()
        .flex()
        .items_baseline()
        .gap_2()
        .text_size(px(15.))
        .line_height(px(20.))
        .child(
            div()
                .font_weight(FontWeight::SEMIBOLD)
                .text_color(palette.ink)
                .child(title.into()),
        )
        .children(count.map(|count| div().text_color(palette.ink_faint).child(count)))
}

/// A section with a plain title.
pub(crate) fn section(title: impl Into<SharedString>, palette: &Palette) -> Div {
    section_shell().child(section_header(palette).child(section_title(title, None, palette)))
}

/// A heading inside a section (a group of findings, the engines).
pub(crate) fn subsection_title(title: impl Into<SharedString>, palette: &Palette) -> Div {
    div()
        .text_size(px(13.))
        .font_weight(FontWeight::SEMIBOLD)
        .text_color(palette.ink_muted)
        .child(title.into())
}

/// The one boxed surface: a raised plate for plots you look at (a spectrum,
/// the volcano, a heat map). Lists and forms sit on the page instead.
pub(crate) fn plate(palette: &Palette) -> Div {
    div()
        .flex()
        .flex_col()
        .gap_3()
        .p_5()
        .rounded(px(4.))
        .border_1()
        .border_color(palette.rule)
        .bg(palette.paper_raised)
}

/// A number and what it counts, for a row of facts: "255 spectra".
pub(crate) fn fact(
    value: impl Into<SharedString>,
    label: impl Into<SharedString>,
    palette: &Palette,
) -> Div {
    div()
        .flex()
        .flex_none()
        .gap_1()
        .child(
            div()
                .font_weight(FontWeight::MEDIUM)
                .text_color(palette.ink)
                .child(value.into()),
        )
        .child(div().text_color(palette.ink_muted).child(label.into()))
}

/// Facts side by side, set apart by space rather than joined by separators.
pub(crate) fn facts() -> Div {
    div()
        .flex()
        .flex_wrap()
        .items_baseline()
        .gap_x(px(20.))
        .gap_y_1()
}

/// A left-ruled message block for errors and findings.
pub(crate) fn callout(color: gpui_kit::Hsla, background: gpui_kit::Hsla) -> Div {
    div()
        .px_4()
        .py_3()
        .rounded(px(2.))
        .border_l_2()
        .border_color(color)
        .bg(background)
}

/// A full-width failure message above the project (e.g. a failed Add data).
pub(crate) fn failure_banner(message: &str, palette: &Palette) -> impl IntoElement {
    div().px(px(48.)).pt_6().child(
        callout(palette.fault, palette.fault_soft)
            .text_size(px(13.))
            .text_color(palette.fault)
            .child(message.to_string()),
    )
}

/// A full-width note above the project about what just happened (an undo).
pub(crate) fn notice_banner(message: &str, palette: &Palette) -> impl IntoElement {
    div().px(px(48.)).pt_6().child(
        callout(palette.signal, palette.signal_soft)
            .text_size(px(13.))
            .text_color(palette.ink)
            .child(message.to_string()),
    )
}

/// `12345` -> `12,345`.
pub(crate) fn group_digits(value: usize) -> String {
    let digits = value.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

#[cfg(test)]
mod tests {
    use super::group_digits;

    #[test]
    fn groups_thousands() {
        assert_eq!(group_digits(0), "0");
        assert_eq!(group_digits(999), "999");
        assert_eq!(group_digits(1_000), "1,000");
        assert_eq!(group_digits(1_234_567), "1,234,567");
    }
}
