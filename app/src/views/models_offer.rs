//! The offer of deep-learning models, above a project when none are installed.

use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::component::{Disableable as _, Sizable as _};
use gpui_kit::{App, ClickEvent, IntoElement, ParentElement, Styled, Window, div, px};

use super::callout;
use crate::theme::Palette;

type Handler = Box<dyn Fn(&ClickEvent, &mut Window, &mut App) + 'static>;

pub struct ModelsOfferProps {
    pub on_download: Handler,
    pub on_install: Handler,
    pub on_dismiss: Handler,
}

pub fn models_offer(palette: &Palette, busy: bool, props: ModelsOfferProps) -> impl IntoElement {
    div().px(px(48.)).pt_6().child(
        callout(palette.signal, palette.signal_soft)
            .flex()
            .items_center()
            .gap_2()
            .text_size(px(13.))
            .text_color(palette.ink)
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .child("Validation is running without deep-learning models (29 MB)."),
            )
            .child(
                Button::new("models-download")
                    .small()
                    .primary()
                    .label("Download")
                    .on_click(props.on_download),
            )
            .child(
                Button::new("models-install")
                    .small()
                    .ghost()
                    .label("Install file…")
                    .disabled(busy)
                    .on_click(props.on_install),
            )
            .child(
                Button::new("models-dismiss")
                    .small()
                    .ghost()
                    .label("Not now")
                    .on_click(props.on_dismiss),
            ),
    )
}
