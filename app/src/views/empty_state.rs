use gpui_kit::component::Disableable as _;
use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::{
    App, ClickEvent, FontWeight, InteractiveElement as _, IntoElement, ParentElement, Styled,
    Window, div, px,
};

use super::{callout, mono};
use crate::theme::{MONO, Palette};

type Handler = Box<dyn Fn(&ClickEvent, &mut Window, &mut App) + 'static>;

pub struct EmptyStateProps<'a> {
    /// What is running, when something is.
    pub busy: Option<&'a str>,
    pub failure: Option<&'a str>,
    pub on_new: Handler,
    pub on_open: Handler,
    pub on_demo: Handler,
}

pub fn empty_state(palette: &Palette, props: EmptyStateProps<'_>) -> impl IntoElement {
    let is_busy = props.busy.is_some();
    let (new_key, open_key) = shortcuts();
    let action = |button: Button, key: &'static str| {
        div()
            .flex()
            .items_center()
            .justify_between()
            .gap_4()
            .child(button)
            .child(
                div()
                    .font_family(MONO)
                    .text_size(px(11.))
                    .text_color(palette.ink_faint)
                    .child(key),
            )
    };
    let column = div()
        .id("empty-copy")
        .w(px(340.))
        .flex()
        .flex_col()
        .gap_2()
        .child(
            div()
                .mb_4()
                .text_size(px(20.))
                .font_weight(FontWeight::SEMIBOLD)
                .child("Tuffold"),
        )
        .child(action(
            Button::new("new-project")
                .primary()
                .label(props.busy.unwrap_or("New project from data…"))
                .loading(is_busy)
                .disabled(is_busy)
                .on_click(props.on_new),
            new_key,
        ))
        .child(action(
            Button::new("open-project")
                .ghost()
                .label("Open project…")
                .disabled(is_busy)
                .on_click(props.on_open),
            open_key,
        ))
        .child(action(
            Button::new("open-demo")
                .ghost()
                .label("Open the demo project")
                .disabled(is_busy)
                .on_click(props.on_demo),
            "",
        ))
        .child(
            div()
                .mt_6()
                .pt_3()
                .border_t_1()
                .border_color(palette.rule)
                .flex()
                .flex_col()
                .gap_1()
                .text_size(px(12.))
                .text_color(palette.ink_faint)
                .child(mono(
                    "FragPipe · MaxQuant · pepXML · mzIdentML · Mascot · PD · Scaffold",
                ))
                .child(mono("mzML · MGF · timsTOF .d")),
        )
        .children(props.failure.map(|message| {
            callout(palette.fault, palette.fault_soft)
                .mt_4()
                .text_size(px(12.))
                .text_color(palette.fault)
                .child(message.to_string())
        }));

    div()
        .size_full()
        .flex()
        .items_center()
        .justify_center()
        .child(column)
}

fn shortcuts() -> (&'static str, &'static str) {
    if cfg!(target_os = "macos") {
        ("⌘N", "⌘O")
    } else {
        ("Ctrl+N", "Ctrl+O")
    }
}
