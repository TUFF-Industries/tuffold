//! Table delegates for the browsing views (gpui-component's virtualized
//! `DataTable`): protein groups (as a list or a biosample matrix), a group's
//! peptides, and a peptide's spectra. Rows are plain store records; cells
//! format them. Every column sorts (click its header); a selection follows
//! its row.

use gpui_kit::component::table::{Column, ColumnSort, TableDelegate, TableState};
use gpui_kit::{
    AnyElement, App, Context, FontWeight, IntoElement, ParentElement, SharedString, Styled, Window,
    div, px,
};
use tuffold_store::{Biosample, GroupPeptide, ProteinListRow, SpectrumListRow};
use tuffold_validate::{ThresholdMode, format_rate};

use crate::browse::{SampleMeasure, sample_value};
use crate::sorted::{Compared, Direction, SortedRows, missing_last, natural};
use crate::theme::{MONO, Palette};
use crate::views::group_digits;

/// Orders two rows ascending by a column.
pub type Compare<R> = Box<dyn Fn(usize, &R, &R) -> Compared>;

/// A table delegate whose rows sort by column.
pub trait Sortable: TableDelegate + Sized + 'static {
    type Row;
    fn rows(&self) -> &SortedRows<Self::Row>;
    fn rows_mut(&mut self) -> &mut SortedRows<Self::Row>;
    fn comparator(&self) -> Compare<Self::Row>;

    /// Replace the rows, keeping the current sort.
    fn set_rows(&mut self, rows: Vec<Self::Row>) {
        let compare = self.comparator();
        self.rows_mut().replace(rows, compare);
    }

    /// The column with its sort indicator (every column is sortable).
    fn mark(&self, col_ix: usize, column: Column) -> Column {
        match self.rows().direction(col_ix) {
            Some(Direction::Ascending) => column.sort(ColumnSort::Ascending),
            Some(Direction::Descending) => column.sort(ColumnSort::Descending),
            None => column.sortable(),
        }
    }
}

/// Sort after the table finishes handling the header click (the delegate
/// cannot reach the selection from inside `perform_sort`), then move the
/// selection to where its row went.
pub(crate) fn sort_later<D: Sortable>(
    col_ix: usize,
    sort: ColumnSort,
    window: &mut Window,
    cx: &mut Context<TableState<D>>,
) {
    let sort = match sort {
        ColumnSort::Ascending => Some((col_ix, Direction::Ascending)),
        ColumnSort::Descending => Some((col_ix, Direction::Descending)),
        ColumnSort::Default => None,
    };
    cx.defer_in(window, move |state, _, cx| {
        let selected = state.selected_row();
        let compare = state.delegate().comparator();
        let moved = state
            .delegate_mut()
            .rows_mut()
            .sort_keeping(sort, compare, selected);
        // Moving the selection scrolls to it; leave the view alone when the
        // selected row stays put.
        if let Some(row) = moved.filter(|&row| Some(row) != selected) {
            state.set_selected_row(row, cx);
        }
        cx.notify();
    });
}

/// A Samples-view cell: counts as whole numbers, NSAF in scientific
/// notation, emPAI to two decimals, "—" when unknown.
fn measure_text(measure: SampleMeasure, value: Option<f64>) -> String {
    match value {
        None => "—".into(),
        Some(v) if measure.is_count() => group_digits(v.round() as usize),
        Some(0.0) => "0".into(),
        Some(v) if measure == SampleMeasure::Nsaf => format!("{v:.2e}"),
        Some(v) => format!("{v:.2}"),
    }
}

fn q_order(a: Option<f64>, b: Option<f64>) -> Compared {
    missing_last(a, b, f64::total_cmp)
}

fn text(value: impl Into<SharedString>) -> AnyElement {
    div().truncate().child(value.into()).into_any_element()
}

fn mono(value: impl Into<SharedString>) -> AnyElement {
    div()
        .font_family(MONO)
        .truncate()
        .child(value.into())
        .into_any_element()
}

fn number(value: impl Into<SharedString>) -> AnyElement {
    div()
        .w_full()
        .flex()
        .justify_end()
        .font_family(MONO)
        .child(value.into())
        .into_any_element()
}

fn q_text(q: Option<f64>) -> String {
    q.map_or_else(|| "—".to_string(), format_rate)
}

/// "99.8%", rounded down so nothing reads as more certain than it is (or
/// seems to pass a threshold by rounding); "—" when unknown.
pub(crate) fn probability_text(probability: Option<f64>) -> String {
    probability.map_or_else(
        || "—".to_string(),
        |p| {
            let tenths = (p * 1000.0 + 1e-9).floor() / 10.0;
            let text = format!("{tenths:.1}");
            format!("{}%", text.strip_suffix(".0").unwrap_or(&text))
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probabilities_round_down() {
        assert_eq!(probability_text(Some(0.999_6)), "99.9%");
        assert_eq!(probability_text(Some(0.989_9)), "98.9%");
        assert_eq!(probability_text(Some(0.95)), "95%");
        assert_eq!(probability_text(Some(1.0)), "100%");
        assert_eq!(probability_text(None), "—");
    }
}

/// The score column: q-values when thresholds are FDRs, probabilities
/// when they are probabilities, so the column matches what decides passing.
fn score_column(mode: ThresholdMode, width: f32) -> Column {
    match mode {
        ThresholdMode::Fdr => Column::new("q", "q"),
        ThresholdMode::Probability => Column::new("probability", "Prob."),
    }
    .width(px(width))
    .text_right()
}

fn score_text(mode: ThresholdMode, q: Option<f64>, probability: Option<f64>) -> String {
    match mode {
        ThresholdMode::Fdr => q_text(q),
        ThresholdMode::Probability => probability_text(probability),
    }
}

fn score_order(
    mode: ThresholdMode,
    (q_a, p_a): (Option<f64>, Option<f64>),
    (q_b, p_b): (Option<f64>, Option<f64>),
) -> Compared {
    match mode {
        ThresholdMode::Fdr => q_order(q_a, q_b),
        ThresholdMode::Probability => q_order(p_a, p_b),
    }
}

/// How protein groups are laid out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinLayout {
    /// Accession, description, peptides, spectra, coverage, q-value.
    List,
    /// Accession, description, spectral counts per biosample, the total.
    Samples,
}

pub struct ProteinTable {
    pub layout: ProteinLayout,
    /// Which score column to show (from the thresholds' mode).
    pub mode: ThresholdMode,
    pub rows: SortedRows<ProteinListRow>,
    pub biosamples: Vec<Biosample>,
    /// What the Samples layout shows per biosample.
    pub measure: SampleMeasure,
    /// Each biosample's NSAF denominator ([`crate::browse::saf_totals`]).
    pub saf_totals: Vec<f64>,
    /// The largest value of the measure, for the heat map's scale.
    pub highest: f64,
}

impl ProteinTable {
    pub fn new(layout: ProteinLayout) -> Self {
        Self {
            layout,
            mode: ThresholdMode::Fdr,
            rows: SortedRows::default(),
            biosamples: Vec::new(),
            measure: SampleMeasure::default(),
            saf_totals: Vec::new(),
            highest: 1.0,
        }
    }
}

impl Sortable for ProteinTable {
    type Row = ProteinListRow;

    fn rows(&self) -> &SortedRows<ProteinListRow> {
        &self.rows
    }

    fn rows_mut(&mut self) -> &mut SortedRows<ProteinListRow> {
        &mut self.rows
    }

    fn comparator(&self) -> Compare<ProteinListRow> {
        let (layout, samples, mode) = (self.layout, self.biosamples.len(), self.mode);
        let (measure, saf_totals) = (self.measure, self.saf_totals.clone());
        Box::new(move |col, a, b| match (layout, col) {
            (_, 0) => natural(&a.accession, &b.accession).into(),
            (_, 1) => missing_last(a.description.as_ref(), b.description.as_ref(), Ord::cmp),
            (ProteinLayout::List, 2) => a.peptides.cmp(&b.peptides).into(),
            (ProteinLayout::List, 3) => a.spectra.cmp(&b.spectra).into(),
            (ProteinLayout::List, 4) => q_order(a.coverage, b.coverage),
            (ProteinLayout::List, _) => {
                score_order(mode, (a.q_value, a.probability), (b.q_value, b.probability))
            }
            (ProteinLayout::Samples, ix) if ix <= samples + 1 => {
                let value = |row: &ProteinListRow| sample_value(row, ix - 2, measure, &saf_totals);
                q_order(value(a), value(b))
            }
            (ProteinLayout::Samples, _) => a.spectra.cmp(&b.spectra).into(),
        })
    }
}

impl TableDelegate for ProteinTable {
    fn columns_count(&self, _: &App) -> usize {
        match self.layout {
            ProteinLayout::List => 6,
            ProteinLayout::Samples => self.biosamples.len() + 3,
        }
    }

    fn rows_count(&self, _: &App) -> usize {
        self.rows.len()
    }

    fn perform_sort(
        &mut self,
        col_ix: usize,
        sort: ColumnSort,
        window: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) {
        sort_later(col_ix, sort, window, cx);
    }

    fn column(&self, col_ix: usize, _: &App) -> Column {
        let column = match (self.layout, col_ix) {
            (_, 0) => Column::new("protein", "Protein").width(px(172.)),
            (ProteinLayout::List, 1) => Column::new("description", "Description").width(px(104.)),
            (ProteinLayout::Samples, 1) => {
                Column::new("description", "Description").width(px(220.))
            }
            (ProteinLayout::List, 2) => Column::new("peptides", "Pep.").width(px(52.)).text_right(),
            (ProteinLayout::List, 3) => Column::new("spectra", "Spectra")
                .width(px(72.))
                .text_right(),
            (ProteinLayout::List, 4) => Column::new("coverage", "Cov.").width(px(58.)).text_right(),
            (ProteinLayout::List, _) => score_column(self.mode, 62.),
            (ProteinLayout::Samples, ix) if ix <= self.biosamples.len() + 1 => {
                let sample = &self.biosamples[ix - 2];
                Column::new(format!("sample-{}", sample.id), sample.name.clone())
                    .width(px(104.))
                    .text_right()
            }
            (ProteinLayout::Samples, _) => Column::new("total", "Total spectra")
                .width(px(104.))
                .text_right(),
        };
        self.mark(col_ix, column)
    }

    fn render_td(
        &mut self,
        row_ix: usize,
        col_ix: usize,
        _: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) -> impl IntoElement {
        let Some(row) = self.rows.get(row_ix) else {
            return div().into_any_element();
        };
        match (self.layout, col_ix) {
            (_, 0) if row.starred => mono(format!("★ {}", row.accession)),
            (_, 0) => mono(row.accession.clone()),
            (_, 1) => text(row.description.clone().unwrap_or_default()),
            (ProteinLayout::List, 2) => number(row.peptides.to_string()),
            (ProteinLayout::List, 3) => number(group_digits(row.spectra)),
            (ProteinLayout::List, 4) => number(
                row.coverage
                    .map_or_else(|| "—".to_string(), |c| format!("{:.0}%", c * 100.0)),
            ),
            (ProteinLayout::List, _) => number(score_text(self.mode, row.q_value, row.probability)),
            (ProteinLayout::Samples, ix) if ix <= self.biosamples.len() + 1 => {
                let palette = *cx.global::<Palette>();
                let value = sample_value(row, ix - 2, self.measure, &self.saf_totals);
                let heat = value.map_or(0.0, |v| (v / self.highest.max(f64::MIN_POSITIVE)) as f32);
                div()
                    .w_full()
                    .h_full()
                    .px_2()
                    .flex()
                    .items_center()
                    .justify_end()
                    .rounded(px(3.))
                    .bg(palette.signal.opacity(0.06 + 0.5 * heat.sqrt()))
                    .font_family(MONO)
                    .text_color(if heat > 0.45 {
                        palette.on_signal
                    } else {
                        palette.ink
                    })
                    .child(measure_text(self.measure, value))
                    .into_any_element()
            }
            (ProteinLayout::Samples, _) => number(group_digits(row.spectra)),
        }
    }

    fn cell_text(&self, row_ix: usize, col_ix: usize, _: &App) -> String {
        let Some(row) = self.rows.get(row_ix) else {
            return String::new();
        };
        match (self.layout, col_ix) {
            (_, 0) => row.accession.clone(),
            (_, 1) => row.description.clone().unwrap_or_default(),
            (ProteinLayout::List, 2) => row.peptides.to_string(),
            (ProteinLayout::List, 3) => row.spectra.to_string(),
            (ProteinLayout::List, 4) => row.coverage.map(|c| format!("{c:.4}")).unwrap_or_default(),
            (ProteinLayout::List, _) => match self.mode {
                ThresholdMode::Fdr => row.q_value,
                ThresholdMode::Probability => row.probability,
            }
            .map(|value| value.to_string())
            .unwrap_or_default(),
            (ProteinLayout::Samples, ix) if ix <= self.biosamples.len() + 1 => {
                sample_value(row, ix - 2, self.measure, &self.saf_totals)
                    .map(|v| v.to_string())
                    .unwrap_or_default()
            }
            (ProteinLayout::Samples, _) => row.spectra.to_string(),
        }
    }
}

#[derive(Default)]
pub struct PeptideTable {
    pub mode: ThresholdMode,
    pub rows: SortedRows<GroupPeptide>,
}

impl Sortable for PeptideTable {
    type Row = GroupPeptide;

    fn rows(&self) -> &SortedRows<GroupPeptide> {
        &self.rows
    }

    fn rows_mut(&mut self) -> &mut SortedRows<GroupPeptide> {
        &mut self.rows
    }

    fn comparator(&self) -> Compare<GroupPeptide> {
        let mode = self.mode;
        Box::new(move |col, a, b| match col {
            0 => a
                .localized_peptide
                .as_ref()
                .unwrap_or(&a.proforma)
                .cmp(b.localized_peptide.as_ref().unwrap_or(&b.proforma))
                .into(),
            1 => q_order(a.min_site_probability, b.min_site_probability),
            2 => a.role.cmp(&b.role).into(),
            3 => a.spectra.cmp(&b.spectra).into(),
            _ => score_order(mode, (a.q_value, a.probability), (b.q_value, b.probability)),
        })
    }
}

impl TableDelegate for PeptideTable {
    fn columns_count(&self, _: &App) -> usize {
        5
    }

    fn rows_count(&self, _: &App) -> usize {
        self.rows.len()
    }

    fn perform_sort(
        &mut self,
        col_ix: usize,
        sort: ColumnSort,
        window: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) {
        sort_later(col_ix, sort, window, cx);
    }

    fn column(&self, col_ix: usize, _: &App) -> Column {
        let column = match col_ix {
            0 => Column::new("peptide", "Peptide").width(px(124.)),
            1 => Column::new("sites", "Sites").width(px(74.)),
            2 => Column::new("role", "Role").width(px(48.)),
            3 => Column::new("spectra", "Spectra")
                .width(px(58.))
                .text_right(),
            _ => score_column(self.mode, 54.),
        };
        self.mark(col_ix, column)
    }

    fn render_td(
        &mut self,
        row_ix: usize,
        col_ix: usize,
        _: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) -> impl IntoElement {
        let Some(row) = self.rows.get(row_ix) else {
            return div().into_any_element();
        };
        match col_ix {
            0 => mono(
                row.localized_peptide
                    .as_ref()
                    .unwrap_or(&row.proforma)
                    .clone(),
            ),
            1 => {
                let palette = *cx.global::<Palette>();
                let text = row.site_probabilities.as_deref().unwrap_or("—").to_string();
                let color = if row.site_probabilities.is_some() {
                    palette.ink_muted
                } else {
                    palette.ink_faint
                };
                div()
                    .w_full()
                    .truncate()
                    .font_family(MONO)
                    .text_size(px(11.))
                    .text_color(color)
                    .child(text)
                    .into_any_element()
            }
            2 => {
                let palette = *cx.global::<Palette>();
                let color = if row.role == "shared" {
                    palette.ink_faint
                } else {
                    palette.ink_muted
                };
                div()
                    .text_size(px(12.))
                    .text_color(color)
                    .child(row.role.clone())
                    .into_any_element()
            }
            3 => number(group_digits(row.spectra)),
            _ => number(score_text(self.mode, row.q_value, row.probability)),
        }
    }

    fn cell_text(&self, row_ix: usize, col_ix: usize, _: &App) -> String {
        self.rows
            .get(row_ix)
            .map_or_else(String::new, |row| match col_ix {
                0 => row
                    .localized_peptide
                    .as_ref()
                    .unwrap_or(&row.proforma)
                    .clone(),
                1 => row.site_probabilities.clone().unwrap_or_default(),
                2 => row.role.clone(),
                3 => row.spectra.to_string(),
                _ => match self.mode {
                    ThresholdMode::Fdr => row.q_value,
                    ThresholdMode::Probability => row.probability,
                }
                .map(|value| value.to_string())
                .unwrap_or_default(),
            })
    }
}

#[derive(Default)]
pub struct SpectrumTable {
    pub mode: ThresholdMode,
    pub rows: SortedRows<SpectrumListRow>,
}

impl Sortable for SpectrumTable {
    type Row = SpectrumListRow;

    fn rows(&self) -> &SortedRows<SpectrumListRow> {
        &self.rows
    }

    fn rows_mut(&mut self) -> &mut SortedRows<SpectrumListRow> {
        &mut self.rows
    }

    fn comparator(&self) -> Compare<SpectrumListRow> {
        let mode = self.mode;
        Box::new(move |col, a, b| match col {
            0 => natural(&a.run, &b.run).into(),
            1 => natural(&a.key, &b.key).into(),
            2 => missing_last(a.charge, b.charge, i64::cmp),
            3 => q_order(a.precursor_mz, b.precursor_mz),
            _ => score_order(
                mode,
                (Some(a.q_value), Some(a.probability)),
                (Some(b.q_value), Some(b.probability)),
            ),
        })
    }
}

impl TableDelegate for SpectrumTable {
    fn columns_count(&self, _: &App) -> usize {
        5
    }

    fn rows_count(&self, _: &App) -> usize {
        self.rows.len()
    }

    fn perform_sort(
        &mut self,
        col_ix: usize,
        sort: ColumnSort,
        window: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) {
        sort_later(col_ix, sort, window, cx);
    }

    fn column(&self, col_ix: usize, _: &App) -> Column {
        let column = match col_ix {
            0 => Column::new("run", "Run").width(px(88.)),
            1 => Column::new("spectrum", "Spectrum").width(px(86.)),
            2 => Column::new("charge", "z").width(px(36.)).text_right(),
            3 => Column::new("mz", "m/z").width(px(78.)).text_right(),
            _ => score_column(self.mode, 58.),
        };
        self.mark(col_ix, column)
    }

    fn render_td(
        &mut self,
        row_ix: usize,
        col_ix: usize,
        _: &mut Window,
        cx: &mut Context<TableState<Self>>,
    ) -> impl IntoElement {
        let Some(row) = self.rows.get(row_ix) else {
            return div().into_any_element();
        };
        let palette = *cx.global::<Palette>();
        let cell = match col_ix {
            0 => mono(row.run.clone()),
            1 => mono(row.key.clone()),
            2 => number(row.charge.map_or_else(String::new, |z| format!("{z}+"))),
            3 => number(
                row.precursor_mz
                    .map_or_else(String::new, |mz| format!("{mz:.4}")),
            ),
            _ => number(score_text(
                self.mode,
                Some(row.q_value),
                Some(row.probability),
            )),
        };
        // Spectra failing the PSM threshold stay listed, dimmed.
        let color = if row.passes {
            palette.ink
        } else {
            palette.ink_faint
        };
        div()
            .w_full()
            .text_color(color)
            .font_weight(FontWeight::NORMAL)
            .child(cell)
            .into_any_element()
    }

    fn cell_text(&self, row_ix: usize, col_ix: usize, _: &App) -> String {
        self.rows
            .get(row_ix)
            .map_or_else(String::new, |row| match col_ix {
                0 => row.run.clone(),
                1 => row.key.clone(),
                2 => row.charge.map(|z| z.to_string()).unwrap_or_default(),
                3 => row
                    .precursor_mz
                    .map(|mz| mz.to_string())
                    .unwrap_or_default(),
                _ => match self.mode {
                    ThresholdMode::Fdr => row.q_value,
                    ThresholdMode::Probability => row.probability,
                }
                .to_string(),
            })
    }
}
