//! Design tokens. Direction: colour only ever means data. The chrome is
//! graphite ink on a cool acrylamide-white (gel, not paper), so the colours
//! that remain carry meaning: the fragment-ion pair (b in Coomassie blue, y in
//! carmine) that also stains every quantity (heat maps, coverage, bars), and
//! ok / caution / fault for data state. Real actions are ink-filled; a choice
//! is a selected ghost button. Plex Sans sets everything, its default figures
//! are tabular so numbers align without a monospace; Plex Mono is reserved for
//! amino-acid sequences, where residues must line up.
//!
//! Tokens are authored in OKLCH (perceptually even steps across themes) and
//! converted once to GPUI colors. The palette follows the system appearance
//! and is also pushed into gpui-component's theme so stock components match.

use std::borrow::Cow;

use gpui_kit::component::theme::{Theme, ThemeMode};
use gpui_kit::{App, Global, Hsla, Rgba, SharedString, Window, px};

pub const SANS: &str = "IBM Plex Sans";
pub const MONO: &str = "IBM Plex Mono";

/// Convert OKLCH (lightness in percent, chroma, hue in degrees) to a GPUI color.
pub fn oklch(lightness_pct: f32, chroma: f32, hue_deg: f32) -> Hsla {
    oklcha(lightness_pct, chroma, hue_deg, 1.0)
}

pub fn oklcha(lightness_pct: f32, chroma: f32, hue_deg: f32, alpha: f32) -> Hsla {
    let l = lightness_pct / 100.0;
    let (sin, cos) = hue_deg.to_radians().sin_cos();
    let (a, b) = (chroma * cos, chroma * sin);
    let l_ = (l + 0.396_337_78 * a + 0.215_803_76 * b).powi(3);
    let m_ = (l - 0.105_561_346 * a - 0.063_854_17 * b).powi(3);
    let s_ = (l - 0.089_484_18 * a - 1.291_485_5 * b).powi(3);
    let linear = [
        4.076_741_7 * l_ - 3.307_711_6 * m_ + 0.230_969_94 * s_,
        -1.268_438 * l_ + 2.609_757_4 * m_ - 0.341_319_38 * s_,
        -0.004_196_086 * l_ - 0.703_418_6 * m_ + 1.707_614_7 * s_,
    ];
    let [r, g, b] = linear.map(|channel| {
        let channel = channel.clamp(0.0, 1.0);
        if channel <= 0.003_130_8 {
            12.92 * channel
        } else {
            1.055 * channel.powf(1.0 / 2.4) - 0.055
        }
    });
    Rgba { r, g, b, a: alpha }.into()
}

#[derive(Debug, Clone, Copy)]
pub struct Palette {
    pub paper: Hsla,
    pub paper_raised: Hsla,
    pub paper_sunken: Hsla,
    pub ink: Hsla,
    pub ink_muted: Hsla,
    pub ink_faint: Hsla,
    pub rule: Hsla,
    pub rule_strong: Hsla,
    /// Interface emphasis (primary buttons, focus, selection marks): ink, not
    /// a hue, so it never competes with the data colours.
    pub signal: Hsla,
    pub signal_hover: Hsla,
    /// Selected rows and notices.
    pub signal_soft: Hsla,
    pub on_signal: Hsla,
    /// Data emphasis, Coomassie blue: quantities (heat maps, bars, coverage),
    /// modified residues, "up" in comparisons. Same hue as the b ions.
    pub stain: Hsla,
    pub stain_soft: Hsla,
    /// Annotated b and y fragment ions in spectra (data series, not state).
    pub ion_b: Hsla,
    pub ion_y: Hsla,
    pub ok: Hsla,
    pub ok_soft: Hsla,
    pub caution: Hsla,
    pub caution_soft: Hsla,
    pub fault: Hsla,
    pub fault_soft: Hsla,
    /// Heat map lightness at z = 0 and at the extremes (OKLCH L, percent);
    /// highs take the signal hue, lows the ion-y hue, as on the volcano.
    pub heat_mid: f32,
    pub heat_end: f32,
}

impl Global for Palette {}

/// Coomassie blue (the b ions, the stain on quantities) and carmine (the y
/// ions, "down"), as OKLCH hues.
const STAIN_HUE: f32 = 272.0;
const CARMINE_HUE: f32 = 358.0;

/// z-scores at or beyond this take the heat map's end colours.
pub const HEAT_CLAMP: f64 = 2.5;

/// The colour scale's end for rows of `n` values: a z-score over `n` values
/// can't exceed (n − 1)/√n (1.5 for four), so few samples would otherwise
/// never reach the strong colours.
pub fn heat_clamp(n: usize) -> f64 {
    let n = n.max(2) as f64;
    ((n - 1.0) / n.sqrt()).min(HEAT_CLAMP)
}

impl Palette {
    /// A heat map cell's colour for z-score `z` (diverging, through a
    /// neutral middle), reaching full colour at `clamp`.
    pub fn heat(&self, z: f64, clamp: f64) -> Hsla {
        let t = (z.abs() / clamp).min(1.0) as f32;
        let lightness = self.heat_mid + (self.heat_end - self.heat_mid) * t;
        let hue = if z >= 0.0 { STAIN_HUE } else { CARMINE_HUE };
        oklch(lightness, 0.19 * t, hue)
    }

    /// A quantity's cell colour, like stain density on a gel: from barely
    /// tinted to mid Coomassie at `t` = 1, never so dark that ink text on it
    /// stops being readable.
    pub fn stain_density(&self, t: f32) -> Hsla {
        let t = t.clamp(0.0, 1.0);
        let (from, to) = if self.paper.l > 0.5 {
            (97.0, 70.0)
        } else {
            (24.0, 44.0)
        };
        oklch(from + (to - from) * t, 0.13 * t, STAIN_HUE)
    }

    /// Ink or paper, whichever reads better on `background` (a heat cell).
    pub fn text_on(&self, background: Hsla) -> Hsla {
        let dark_ink = contrast(self.ink, background) >= contrast(self.paper_raised, background);
        if dark_ink {
            self.ink
        } else {
            self.paper_raised
        }
    }

    pub fn light() -> Self {
        let stain = oklch(45.0, 0.2, STAIN_HUE);
        Self {
            paper: oklch(97.8, 0.004, 230.0),
            paper_raised: oklch(99.7, 0.002, 230.0),
            paper_sunken: oklch(94.8, 0.006, 235.0),
            ink: oklch(23.0, 0.02, 255.0),
            ink_muted: oklch(45.0, 0.018, 255.0),
            ink_faint: oklch(52.0, 0.014, 255.0),
            rule: oklch(90.5, 0.006, 240.0),
            rule_strong: oklch(80.0, 0.01, 245.0),
            signal: oklch(25.0, 0.02, 255.0),
            signal_hover: oklch(35.0, 0.022, 255.0),
            signal_soft: oklch(92.0, 0.012, 250.0),
            on_signal: oklch(98.5, 0.003, 230.0),
            stain,
            stain_soft: oklch(92.5, 0.04, STAIN_HUE),
            ion_b: stain,
            ion_y: oklch(54.0, 0.21, CARMINE_HUE),
            ok: oklch(52.0, 0.12, 158.0),
            ok_soft: oklch(94.5, 0.04, 158.0),
            caution: oklch(52.0, 0.12, 70.0),
            caution_soft: oklch(95.0, 0.05, 88.0),
            fault: oklch(52.0, 0.19, 28.0),
            fault_soft: oklch(95.0, 0.03, 28.0),
            heat_mid: 96.0,
            heat_end: 48.0,
        }
    }

    pub fn dark() -> Self {
        let stain = oklch(71.0, 0.15, STAIN_HUE);
        Self {
            paper: oklch(18.0, 0.008, 250.0),
            paper_raised: oklch(21.5, 0.009, 250.0),
            paper_sunken: oklch(15.5, 0.008, 250.0),
            ink: oklch(93.5, 0.005, 230.0),
            ink_muted: oklch(75.0, 0.01, 240.0),
            ink_faint: oklch(61.0, 0.012, 245.0),
            rule: oklch(28.5, 0.01, 250.0),
            rule_strong: oklch(38.0, 0.012, 250.0),
            signal: oklch(91.0, 0.006, 230.0),
            signal_hover: oklch(83.0, 0.008, 230.0),
            signal_soft: oklch(30.0, 0.014, 250.0),
            on_signal: oklch(18.0, 0.008, 250.0),
            stain,
            stain_soft: oklch(30.0, 0.07, STAIN_HUE),
            ion_b: stain,
            ion_y: oklch(72.0, 0.17, CARMINE_HUE),
            ok: oklch(75.0, 0.13, 158.0),
            ok_soft: oklch(28.0, 0.05, 158.0),
            caution: oklch(80.0, 0.13, 80.0),
            caution_soft: oklch(29.0, 0.05, 75.0),
            fault: oklch(72.0, 0.16, 28.0),
            fault_soft: oklch(29.0, 0.06, 28.0),
            heat_mid: 25.0,
            heat_end: 70.0,
        }
    }
}

/// WCAG relative luminance.
fn luminance(color: Hsla) -> f32 {
    let rgba = color.to_rgb();
    let linear = |channel: f32| {
        if channel <= 0.040_45 {
            channel / 12.92
        } else {
            ((channel + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * linear(rgba.r) + 0.7152 * linear(rgba.g) + 0.0722 * linear(rgba.b)
}

/// WCAG contrast ratio between two colours (1 to 21).
fn contrast(a: Hsla, b: Hsla) -> f32 {
    let (a, b) = (luminance(a), luminance(b));
    (a.max(b) + 0.05) / (a.min(b) + 0.05)
}

/// Register the bundled IBM Plex fonts (SIL OFL 1.1, see assets/fonts/OFL.txt).
pub fn load_fonts(cx: &mut App) -> anyhow::Result<()> {
    let fonts: Vec<Cow<'static, [u8]>> = vec![
        Cow::Borrowed(include_bytes!("../assets/fonts/IBMPlexSans-Regular.ttf")),
        Cow::Borrowed(include_bytes!("../assets/fonts/IBMPlexSans-Medium.ttf")),
        Cow::Borrowed(include_bytes!("../assets/fonts/IBMPlexSans-SemiBold.ttf")),
        Cow::Borrowed(include_bytes!("../assets/fonts/IBMPlexMono-Regular.ttf")),
        Cow::Borrowed(include_bytes!("../assets/fonts/IBMPlexMono-Medium.ttf")),
    ];
    cx.text_system().add_fonts(fonts)
}

/// Match the system appearance and apply the Tuffold palette on top of the
/// component theme. Call at startup and whenever the appearance changes.
pub fn sync_with_system(window: &mut Window, cx: &mut App) {
    Theme::sync_system_appearance(Some(window), cx);
    let mode = cx.global::<Theme>().mode;
    apply(mode, cx);
}

/// Apply the Tuffold palette for an explicit mode (snapshots, tests).
pub fn apply(mode: ThemeMode, cx: &mut App) {
    Theme::change(mode, None, cx);
    let palette = match mode {
        ThemeMode::Dark => Palette::dark(),
        ThemeMode::Light => Palette::light(),
    };
    Theme::update(cx, |theme| {
        theme.font_family = SharedString::new_static(SANS);
        theme.mono_font_family = SharedString::new_static(MONO);
        theme.font_size = px(13.);
        theme.radius = px(4.);
        let colors = &mut theme.colors;
        colors.background = palette.paper;
        colors.foreground = palette.ink;
        colors.border = palette.rule;
        colors.input = palette.rule_strong;
        colors.muted = palette.paper_sunken;
        colors.muted_foreground = palette.ink_muted;
        colors.primary = palette.signal;
        colors.primary_hover = palette.signal_hover;
        colors.primary_foreground = palette.on_signal;
        colors.button_primary = palette.signal;
        colors.button_primary_hover = palette.signal_hover;
        colors.button_primary_active = palette.signal_hover;
        colors.button_primary_foreground = palette.on_signal;
        colors.button = palette.paper_raised;
        colors.button_hover = palette.paper_sunken;
        colors.button_active = palette.paper_sunken;
        colors.button_foreground = palette.ink;
        colors.secondary = palette.paper_raised;
        colors.secondary_hover = palette.paper_sunken;
        // Selected ghost buttons: the views' segmented choices.
        colors.secondary_active = palette.rule;
        colors.ring = palette.signal;
        colors.success = palette.ok;
        colors.warning = palette.caution;
        colors.danger = palette.fault;
        colors.table = palette.paper_raised;
        colors.table_even = palette.paper_raised;
        colors.table_head = palette.paper_sunken;
        colors.table_head_foreground = palette.ink_muted;
        colors.table_hover = palette.paper_sunken;
        colors.table_active = palette.signal_soft;
        colors.table_active_border = palette.ink_faint;
        colors.table_row_border = palette.rule;
        colors.scrollbar_thumb = palette.rule_strong;
    });
    cx.set_global(palette);
}

#[cfg(test)]
mod heat_tests {
    use super::*;

    #[test]
    fn the_heat_scale_fits_the_possible_z_scores() {
        assert!((heat_clamp(4) - 1.5).abs() < 1e-12);
        assert_eq!(heat_clamp(100), HEAT_CLAMP);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rgb(color: Hsla) -> [u8; 3] {
        let rgba = color.to_rgb();
        [rgba.r, rgba.g, rgba.b].map(|channel| (channel * 255.0).round() as u8)
    }

    #[test]
    fn heat_cells_take_the_more_readable_text() {
        // Mid-ramp cells (OKLCH L near 60) sit between ink and paper, so no
        // text reaches 4.5:1 there; 3:1 is the floor anywhere on the ramp.
        for palette in [Palette::light(), Palette::dark()] {
            for step in 0..=20 {
                let t = f64::from(step) / 20.0;
                for z in [t, -t] {
                    let background = palette.heat(z, 1.0);
                    let text = palette.text_on(background);
                    assert!(
                        contrast(text, background) >= 3.0,
                        "z {z}: contrast {}",
                        contrast(text, background)
                    );
                }
            }
            for step in 0..=20 {
                let background = palette.stain_density(step as f32 / 20.0);
                assert!(contrast(palette.ink, background) >= 4.5);
            }
        }
    }

    #[test]
    fn oklch_matches_reference_srgb_values() {
        assert_eq!(rgb(oklch(100.0, 0.0, 0.0)), [255, 255, 255]);
        assert_eq!(rgb(oklch(0.0, 0.0, 0.0)), [0, 0, 0]);
        // CSS Color 4 reference: oklch(62.8% 0.2577 29.23) is pure sRGB red.
        assert_eq!(rgb(oklch(62.8, 0.2577, 29.23)), [255, 0, 0]);
        // oklch(45.2% 0.313 264.05) is pure sRGB blue.
        assert_eq!(rgb(oklch(45.2, 0.313, 264.05)), [0, 0, 255]);
    }

    #[test]
    fn stain_is_blue_and_paper_is_light_in_the_light_theme() {
        let palette = Palette::light();
        let [r, g, b] = rgb(palette.stain);
        assert!(b > r && b > g, "stain should read blue: {r},{g},{b}");
        assert!(palette.paper.l > 0.9);
        assert!(Palette::dark().paper.l < 0.2);
    }
}
