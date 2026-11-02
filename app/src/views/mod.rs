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

use crate::theme::{MONO, Palette};

/// A small label above a heading.
pub(crate) fn eyebrow(text: impl Into<SharedString>, palette: &Palette) -> Div {
    div()
        .text_size(px(12.))
        .text_color(palette.ink_faint)
        .child(text.into())
}

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

/// A raised surface; callers add their own title row.
pub(crate) fn panel_shell(palette: &Palette) -> Div {
    div()
        .flex()
        .flex_col()
        .gap_3()
        .px_6()
        .pt_4()
        .pb_6()
        .rounded(px(6.))
        .border_1()
        .border_color(palette.rule)
        .bg(palette.paper_raised)
        .shadow_sm()
}

pub(crate) fn panel_title(title: impl Into<SharedString>) -> Div {
    div()
        .text_size(px(13.))
        .font_weight(FontWeight::SEMIBOLD)
        .child(title.into())
}

/// A raised surface with a title row.
pub(crate) fn panel(title: impl Into<SharedString>, palette: &Palette) -> Div {
    panel_shell(palette).child(panel_title(title))
}

/// A left-ruled message block for errors and findings.
pub(crate) fn callout(color: gpui_kit::Hsla, background: gpui_kit::Hsla) -> Div {
    div()
        .px_4()
        .py_3()
        .rounded(px(3.))
        .border_l_3()
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

pub(crate) fn mono(text: impl Into<SharedString>) -> impl IntoElement {
    div().font_family(MONO).child(text.into())
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
