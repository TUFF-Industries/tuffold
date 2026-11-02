//! Design tokens. Direction: an instrument panel printed on lab paper — warm
//! paper surfaces, ink typography, one instrument-blue signal color, and
//! semantic ok / caution / fault colors that only ever mean data state. Blue
//! fill is for real actions; a choice is a selected ghost button.
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
    pub signal: Hsla,
    pub signal_hover: Hsla,
    /// Selection and highlight backgrounds (selected table rows, coverage).
    pub signal_soft: Hsla,
    pub on_signal: Hsla,
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
        let hue = if z >= 0.0 { 263.0 } else { 350.0 };
        oklch(lightness, 0.19 * t, hue)
    }

    pub fn light() -> Self {
        Self {
            paper: oklch(97.6, 0.007, 85.0),
            paper_raised: oklch(99.4, 0.003, 85.0),
            paper_sunken: oklch(95.0, 0.009, 85.0),
            ink: oklch(22.0, 0.014, 262.0),
            ink_muted: oklch(46.0, 0.016, 262.0),
            ink_faint: oklch(64.0, 0.012, 262.0),
            rule: oklch(89.0, 0.009, 85.0),
            rule_strong: oklch(79.0, 0.012, 85.0),
            signal: oklch(52.0, 0.2, 263.0),
            signal_hover: oklch(46.0, 0.2, 263.0),
            signal_soft: oklch(93.0, 0.035, 263.0),
            on_signal: oklch(99.0, 0.005, 263.0),
            ion_b: oklch(52.0, 0.2, 263.0),
            ion_y: oklch(55.0, 0.2, 350.0),
            ok: oklch(52.0, 0.12, 158.0),
            ok_soft: oklch(94.5, 0.04, 158.0),
            caution: oklch(52.0, 0.12, 68.0),
            caution_soft: oklch(95.0, 0.055, 85.0),
            fault: oklch(52.0, 0.19, 28.0),
            fault_soft: oklch(95.0, 0.03, 28.0),
            heat_mid: 95.0,
            heat_end: 50.0,
        }
    }

    pub fn dark() -> Self {
        Self {
            paper: oklch(17.5, 0.012, 262.0),
            paper_raised: oklch(21.0, 0.014, 262.0),
            paper_sunken: oklch(15.0, 0.011, 262.0),
            ink: oklch(93.0, 0.008, 85.0),
            ink_muted: oklch(74.0, 0.012, 85.0),
            ink_faint: oklch(56.0, 0.012, 262.0),
            rule: oklch(28.0, 0.014, 262.0),
            rule_strong: oklch(38.0, 0.016, 262.0),
            signal: oklch(70.0, 0.15, 263.0),
            signal_hover: oklch(76.0, 0.13, 263.0),
            signal_soft: oklch(30.0, 0.07, 263.0),
            on_signal: oklch(16.0, 0.03, 263.0),
            ion_b: oklch(72.0, 0.14, 263.0),
            ion_y: oklch(73.0, 0.15, 350.0),
            ok: oklch(75.0, 0.13, 158.0),
            ok_soft: oklch(28.0, 0.05, 158.0),
            caution: oklch(80.0, 0.13, 80.0),
            caution_soft: oklch(29.0, 0.05, 75.0),
            fault: oklch(72.0, 0.16, 28.0),
            fault_soft: oklch(29.0, 0.06, 28.0),
            heat_mid: 26.0,
            heat_end: 72.0,
        }
    }
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
        colors.table_active_border = palette.signal;
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
    fn oklch_matches_reference_srgb_values() {
        assert_eq!(rgb(oklch(100.0, 0.0, 0.0)), [255, 255, 255]);
        assert_eq!(rgb(oklch(0.0, 0.0, 0.0)), [0, 0, 0]);
        // CSS Color 4 reference: oklch(62.8% 0.2577 29.23) is pure sRGB red.
        assert_eq!(rgb(oklch(62.8, 0.2577, 29.23)), [255, 0, 0]);
        // oklch(45.2% 0.313 264.05) is pure sRGB blue.
        assert_eq!(rgb(oklch(45.2, 0.313, 264.05)), [0, 0, 255]);
    }

    #[test]
    fn signal_is_blue_and_paper_is_light_in_the_light_theme() {
        let palette = Palette::light();
        let [r, g, b] = rgb(palette.signal);
        assert!(b > r && b > g, "signal should read blue: {r},{g},{b}");
        assert!(palette.paper.l > 0.9);
        assert!(Palette::dark().paper.l < 0.2);
    }
}
