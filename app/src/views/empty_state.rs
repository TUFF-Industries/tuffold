use gpui_kit::component::Disableable as _;
use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::{
    App, ClickEvent, FontWeight, InteractiveElement as _, IntoElement, ParentElement, Styled,
    Window, div, px,
};

use super::callout;
use super::spectrum_plot::{LadderScale, cleavage_marks};
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

/// The wordmark: TUFFOLD is a valid peptide (U is selenocysteine, O
/// pyrrolysine), so the name is set as one, with a fragment ladder like the
/// spectrum view's. The series have gaps, as real ones do.
const WORDMARK: &str = "TUFFOLD";
const WORDMARK_B: [bool; 6] = [true, true, true, false, true, false];
const WORDMARK_Y: [bool; 6] = [false, true, true, true, true, true];
const WORDMARK_SCALE: LadderScale = LadderScale {
    letter: 40.,
    cell: 30.,
    gap: 16.,
    tick: 16.,
    stroke: 3.,
};

fn wordmark(palette: &Palette) -> impl IntoElement {
    let scale = WORDMARK_SCALE;
    let residues: Vec<char> = WORDMARK.chars().collect();
    let cells = residues.iter().enumerate().flat_map(|(index, residue)| {
        let letter = div()
            .w(px(scale.cell))
            .flex()
            .justify_center()
            .font_family(MONO)
            .font_weight(FontWeight::MEDIUM)
            .text_size(px(scale.letter))
            .line_height(px(scale.letter * 1.25))
            .text_color(palette.ink)
            .child(residue.to_string());
        let bond = (index + 1 < residues.len())
            .then(|| cleavage_marks(WORDMARK_B[index], WORDMARK_Y[index], scale, palette));
        std::iter::once(letter.into_any_element()).chain(bond.map(|b| b.into_any_element()))
    });
    div().flex().items_center().children(cells)
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
                    .text_size(px(12.))
                    .text_color(palette.ink_faint)
                    .child(key),
            )
    };
    let formats = |label: &'static str, list: &'static str| {
        div()
            .flex()
            .gap_4()
            .child(
                div()
                    .w(px(112.))
                    .flex_none()
                    .text_color(palette.ink_muted)
                    .child(label),
            )
            .child(div().flex_1().min_w_0().text_color(palette.ink).child(list))
    };
    let column = div()
        .id("empty-copy")
        .w(px(520.))
        .flex()
        .flex_col()
        .child(wordmark(palette))
        .child(
            div()
                .mt_5()
                .mb_8()
                .text_size(px(15.))
                .line_height(px(22.))
                .text_color(palette.ink_muted)
                .child("Validate MS/MS search results, infer proteins, and compare samples."),
        )
        .child(
            div()
                .flex()
                .flex_col()
                .gap_1()
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
                )),
        )
        .child(
            div()
                .mt_8()
                .pt_4()
                .border_t_1()
                .border_color(palette.rule)
                .flex()
                .flex_col()
                .gap_2()
                .text_size(px(12.))
                .line_height(px(18.))
                .child(formats(
                    "Search results",
                    "FragPipe, MaxQuant, pepXML, mzIdentML, Mascot, Proteome Discoverer, \
                     Scaffold",
                ))
                .child(formats("Spectra", "mzML, MGF, timsTOF .d")),
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
        .pb(px(64.))
        .child(column)
}

fn shortcuts() -> (&'static str, &'static str) {
    if cfg!(target_os = "macos") {
        ("⌘N", "⌘O")
    } else {
        ("Ctrl+N", "Ctrl+O")
    }
}
