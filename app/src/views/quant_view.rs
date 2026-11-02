//! The Quant tab: compare two categories of biosamples (volcano plot and
//! table) by spectral counts (Welch or Fisher) or by MS1 intensity (median
//! polish, Welch), edit which category each biosample belongs to, and export the
//! table. Loading and saving run off the UI thread.

use std::path::PathBuf;

use gpui_kit::component::Selectable as _;
use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::component::input::{Input, InputState};
use gpui_kit::component::table::{
    Column, ColumnSort, DataTable, TableDelegate, TableEvent, TableState,
};
use gpui_kit::component::{Disableable as _, Sizable as _};
use gpui_kit::{
    AnyElement, App, AppContext as _, Bounds, Context, Div, Entity, EventEmitter, Hsla,
    InteractiveElement as _, IntoElement, ParentElement, PathPromptOptions, Render,
    StatefulInteractiveElement as _, Styled, Subscription, Window, canvas, div, fill, point, px,
    relative, size,
};
use tuffold_quant::intensity::IntensityComparison;
use tuffold_quant::{Comparison, PrimaryTest};

use super::qc_view::qc_body;
use super::tables::{Compare, Sortable, sort_later};
use super::{callout, group_digits, panel_shell, panel_title, track};
use crate::quant::{
    self, AllShown, Certificate, CompareRow, Export, Method, QcShown, QuantData, SIGNIFICANCE,
    Scope, Shown, Volcano, ms1_status, plex_title, reporter_reagent_hint, reporter_status,
    shown_all, shown_counts, shown_intensity, shown_qc, shown_reporters, volcano,
};
use crate::sorted::{SortedRows, missing_last, natural};
use crate::theme::{MONO, Palette};
use tuffold_core::reporter::Plex;
use tuffold_quant::anova::AnovaRow;
use tuffold_quant::heatmap::Heatmap;

const PLOT_HEIGHT: f32 = 320.;

/// The user asked to open a protein group in the Proteins view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenGroup(pub i64);

pub struct ComparisonTable {
    rows: SortedRows<CompareRow>,
    names: (String, String),
}

fn number(value: String) -> AnyElement {
    div()
        .w_full()
        .flex()
        .justify_end()
        .font_family(MONO)
        .child(value)
        .into_any_element()
}

/// A value, or "—" when unknown.
fn shown(value: Option<f64>, format: impl Fn(f64) -> String) -> String {
    value.map_or_else(|| "—".to_string(), format)
}

impl Sortable for ComparisonTable {
    type Row = CompareRow;

    fn rows(&self) -> &SortedRows<CompareRow> {
        &self.rows
    }

    fn rows_mut(&mut self) -> &mut SortedRows<CompareRow> {
        &mut self.rows
    }

    fn comparator(&self) -> Compare<CompareRow> {
        Box::new(move |col, a, b| match col {
            0 => natural(&a.accession, &b.accession).into(),
            1 => missing_last(a.mean_a, b.mean_a, f64::total_cmp),
            2 => missing_last(a.mean_b, b.mean_b, f64::total_cmp),
            3 => missing_last(a.log2_fold, b.log2_fold, f64::total_cmp),
            _ => missing_last(a.q, b.q, f64::total_cmp),
        })
    }
}

impl TableDelegate for ComparisonTable {
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
            0 => Column::new("protein", "Protein").width(px(176.)),
            1 => Column::new("a", self.names.0.clone())
                .width(px(78.))
                .text_right(),
            2 => Column::new("b", self.names.1.clone())
                .width(px(78.))
                .text_right(),
            3 => Column::new("fold", "log2 fold").width(px(76.)).text_right(),
            _ => Column::new("q", "q").width(px(72.)).text_right(),
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
        let q = row.q;
        match col_ix {
            0 => div()
                .font_family(MONO)
                .truncate()
                .child(row.accession.clone())
                .into_any_element(),
            1 => number(shown(row.mean_a, |v| format!("{v:.1}"))),
            2 => number(shown(row.mean_b, |v| format!("{v:.1}"))),
            3 => number(shown(row.log2_fold, |v| format!("{v:+.2}"))),
            _ => div()
                .w_full()
                .flex()
                .justify_end()
                .font_family(MONO)
                .text_color(if q.is_some_and(|q| q <= SIGNIFICANCE) {
                    palette.signal
                } else {
                    palette.ink_muted
                })
                .child(q.map_or_else(|| "—".to_string(), |q| format!("{q:.1e}")))
                .into_any_element(),
        }
    }

    fn cell_text(&self, row_ix: usize, col_ix: usize, _: &App) -> String {
        self.rows
            .get(row_ix)
            .map_or_else(String::new, |row| match col_ix {
                0 => row.accession.clone(),
                1 => row.mean_a.map(|v| v.to_string()).unwrap_or_default(),
                2 => row.mean_b.map(|v| v.to_string()).unwrap_or_default(),
                3 => row.log2_fold.map(|v| v.to_string()).unwrap_or_default(),
                _ => row.q.map(|q| q.to_string()).unwrap_or_default(),
            })
    }
}

/// The ANOVA table: every category compared at once.
pub struct AnovaTable {
    rows: SortedRows<AnovaRow>,
}

impl Sortable for AnovaTable {
    type Row = AnovaRow;

    fn rows(&self) -> &SortedRows<AnovaRow> {
        &self.rows
    }

    fn rows_mut(&mut self) -> &mut SortedRows<AnovaRow> {
        &mut self.rows
    }

    fn comparator(&self) -> Compare<AnovaRow> {
        Box::new(move |col, a, b| match col {
            0 => natural(&a.accession, &b.accession).into(),
            1 => missing_last(a.range, b.range, f64::total_cmp),
            2 => missing_last(a.f, b.f, f64::total_cmp),
            _ => missing_last(a.q, b.q, f64::total_cmp),
        })
    }
}

impl TableDelegate for AnovaTable {
    fn columns_count(&self, _: &App) -> usize {
        4
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
            0 => Column::new("protein", "Protein").width(px(196.)),
            1 => Column::new("range", "log2 range")
                .width(px(86.))
                .text_right(),
            2 => Column::new("f", "F").width(px(76.)).text_right(),
            _ => Column::new("q", "q").width(px(72.)).text_right(),
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
        match col_ix {
            0 => div()
                .font_family(MONO)
                .truncate()
                .child(row.accession.clone())
                .into_any_element(),
            1 => number(shown(row.range, |v| format!("{v:.2}"))),
            2 => number(shown(row.f, |v| format!("{v:.1}"))),
            _ => div()
                .w_full()
                .flex()
                .justify_end()
                .font_family(MONO)
                .text_color(if row.q.is_some_and(|q| q <= SIGNIFICANCE) {
                    palette.signal
                } else {
                    palette.ink_muted
                })
                .child(shown(row.q, |q| format!("{q:.1e}")))
                .into_any_element(),
        }
    }

    fn cell_text(&self, row_ix: usize, col_ix: usize, _: &App) -> String {
        self.rows
            .get(row_ix)
            .map_or_else(String::new, |row| match col_ix {
                0 => row.accession.clone(),
                1 => row.range.map(|v| v.to_string()).unwrap_or_default(),
                2 => row.f.map(|v| v.to_string()).unwrap_or_default(),
                _ => row.q.map(|q| q.to_string()).unwrap_or_default(),
            })
    }
}

pub struct QuantView {
    path: PathBuf,
    data: Option<QuantData>,
    pair: Option<(String, String)>,
    method: Method,
    test: PrimaryTest,
    comparison: Option<Comparison>,
    intensity_comparison: Option<IntensityComparison>,
    reporter_comparison: Option<IntensityComparison>,
    /// Why the intensity comparison can't be shown (no MS1 data, ...).
    intensity_problem: Option<String>,
    /// Why the reporter comparison can't be shown (no plex, ...).
    reporter_problem: Option<String>,
    /// Why the two-category counts comparison can't be shown.
    pair_problem: Option<String>,
    shown: Option<Shown>,
    table: Entity<TableState<ComparisonTable>>,
    scope: Scope,
    all: Option<AllShown>,
    /// Why every category can't be compared (fewer than two, no MS1, ...).
    all_problem: Option<String>,
    qc: Option<QcShown>,
    qc_problem: Option<String>,
    anova_table: Entity<TableState<AnovaTable>>,
    /// One category input per biosample, rebuilt when the samples change.
    inputs: Vec<(String, Entity<InputState>)>,
    /// Name of the stored plex selected in the TMT plex editor (`None` when
    /// defining a new plex or when the project has none yet).
    plex_selected: Option<String>,
    /// Whether the panel is defining a new plex alongside existing ones.
    plex_creating: bool,
    /// Name input for a new plex.
    plex_name_input: Entity<InputState>,
    /// Reagent set selected in the TMT plex editor.
    plex_reagents: Plex,
    plex_synced: bool,
    /// MS runs selected for the plex in the editor.
    plex_runs: Vec<String>,
    /// One biosample input per reporter channel in `plex_reagents`.
    plex_inputs: Vec<(String, Entity<InputState>)>,
    /// Channels marked as holding the pooled reference sample.
    plex_references: Vec<String>,
    /// Whether the user has modified `plex_references` in the editor since
    /// selecting the plex or reagent set.
    plex_references_touched: bool,
    /// Whether the first click on Remove plex has happened and awaits
    /// confirmation.
    confirm_remove_plex: bool,
    busy: bool,
    /// A save is running: saves are never dropped, and a reload asked for
    /// meanwhile runs after it.
    saving: bool,
    reload_after_save: bool,
    exporting: bool,
    /// The spectra-file dialog for matching between runs is open.
    picking: bool,
    failure: Option<String>,
    notice: Option<String>,
    generation: u64,
    _subscriptions: Vec<Subscription>,
}

impl EventEmitter<OpenGroup> for QuantView {}

impl QuantView {
    pub fn new(
        path: PathBuf,
        data: Option<QuantData>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let table = cx.new(|cx| {
            TableState::new(
                ComparisonTable {
                    rows: SortedRows::default(),
                    names: (String::new(), String::new()),
                },
                window,
                cx,
            )
        });
        let subscription =
            cx.subscribe_in(&table, window, |_, table, event: &TableEvent, _, cx| {
                if let TableEvent::DoubleClickedRow(row) = event
                    && let Some(id) = table.read(cx).delegate().rows.get(*row).map(|r| r.group_id)
                {
                    cx.emit(OpenGroup(id));
                }
            });
        let anova_table = cx.new(|cx| {
            TableState::new(
                AnovaTable {
                    rows: SortedRows::default(),
                },
                window,
                cx,
            )
        });
        let anova_subscription = cx.subscribe_in(
            &anova_table,
            window,
            |_, table, event: &TableEvent, _, cx| {
                if let TableEvent::DoubleClickedRow(row) = event
                    && let Some(id) = table.read(cx).delegate().rows.get(*row).map(|r| r.group_id)
                {
                    cx.emit(OpenGroup(id));
                }
            },
        );
        let plex_name_input = cx.new(|cx| InputState::new(window, cx).placeholder("plex name"));
        let mut view = Self {
            path,
            data: None,
            pair: None,
            method: Method::Counts,
            test: PrimaryTest::Fisher,
            comparison: None,
            intensity_comparison: None,
            reporter_comparison: None,
            intensity_problem: None,
            reporter_problem: None,
            pair_problem: None,
            shown: None,
            table,
            scope: Scope::Pair,
            all: None,
            all_problem: None,
            qc: None,
            qc_problem: None,
            anova_table,
            inputs: Vec::new(),
            plex_selected: None,
            plex_creating: false,
            plex_name_input,
            plex_reagents: Plex::Tmt6,
            plex_synced: false,
            plex_runs: Vec::new(),
            plex_inputs: Vec::new(),
            plex_references: Vec::new(),
            plex_references_touched: false,
            confirm_remove_plex: false,
            busy: false,
            saving: false,
            reload_after_save: false,
            exporting: false,
            picking: false,
            failure: None,
            notice: None,
            generation: 0,
            _subscriptions: vec![subscription, anova_subscription],
        };
        match data {
            Some(data) => view.apply(data, cx),
            None => view.reload(cx),
        }
        view
    }

    /// Quantify again (after validation, a threshold change, or new data).
    pub fn reload(&mut self, cx: &mut Context<Self>) {
        if self.saving {
            self.reload_after_save = true;
            return;
        }
        let path = self.path.clone();
        self.run(
            cx,
            move || quant::load(&path),
            |this, data, cx| this.apply(data, cx),
        );
    }

    fn run<T: Send + 'static>(
        &mut self,
        cx: &mut Context<Self>,
        work: impl FnOnce() -> Result<T, String> + Send + 'static,
        done: impl FnOnce(&mut Self, T, &mut Context<Self>) + 'static,
    ) {
        self.generation += 1;
        let generation = self.generation;
        self.busy = true;
        self.notice = None;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let result = cx.background_executor().spawn(async move { work() }).await;
            let _ = this.update(cx, |this, cx| {
                if this.generation != generation {
                    return;
                }
                this.busy = false;
                match result {
                    Ok(value) => {
                        this.failure = None;
                        done(this, value, cx);
                    }
                    Err(message) => this.failure = Some(message),
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn apply(&mut self, data: QuantData, cx: &mut Context<Self>) {
        self.confirm_remove_plex = false;
        // Keep a half-typed "New plex" draft across reloads; saving or
        // removing a plex clears `plex_creating` before calling `apply`.
        if self.plex_creating {
            let free: Vec<String> = data
                .runs_for_plex(None)
                .into_iter()
                .filter(|r| r.in_other_plex.is_none())
                .map(|r| r.name)
                .collect();
            self.plex_runs.retain(|r| free.contains(r));
        } else {
            let chosen = self
                .plex_selected
                .as_deref()
                .and_then(|name| data.plex_named(name))
                .or_else(|| data.plexes.first())
                .cloned();
            if let Some(plex) = &chosen {
                self.plex_selected = Some(plex.name.clone());
                self.plex_reagents = plex.plex;
                self.plex_runs = plex.runs.clone();
                self.plex_references = plex.references.clone();
            } else {
                self.plex_selected = None;
                self.plex_runs.clear();
                self.plex_references.clear();
            }
            self.plex_references_touched = false;
            self.plex_synced = false;
        }
        self.inputs.clear();
        let categories = data.categories_for(self.method);
        let still_valid = self
            .pair
            .as_ref()
            .is_some_and(|(a, b)| categories.contains(a) && categories.contains(b));
        if !still_valid {
            self.pair = data.default_pair_for(self.method);
        }
        self.data = Some(data);
        // Keep the user's test choice unless the pair itself changed.
        self.compare(!still_valid, cx);
    }

    /// Recompute comparisons; `reset_test` picks the counts' primary test.
    fn compare(&mut self, reset_test: bool, cx: &mut Context<Self>) {
        self.comparison = None;
        self.intensity_comparison = None;
        self.reporter_comparison = None;
        self.intensity_problem = None;
        self.reporter_problem = None;
        self.pair_problem = None;
        if let Some(data) = &self.data {
            if data.plexes.is_empty() {
                self.reporter_problem = Some(
                    "This project has no plex yet; choose the reagent set, name the \
                     biosample in each channel below, and click Save plex."
                        .into(),
                );
            } else if let Err(problem) = &data.reporters {
                self.reporter_problem = Some(problem.clone());
            }
            if let Some((a, b)) = &self.pair {
                match data.counts.compare(a, b) {
                    Ok(comparison) => {
                        if reset_test {
                            self.test = comparison.primary;
                        }
                        self.comparison = Some(comparison);
                    }
                    Err(error) => self.pair_problem = Some(error.to_string()),
                }
                match &data.intensity {
                    Ok(quant) => match quant.quantities.compare(a, b) {
                        Ok(comparison) => self.intensity_comparison = Some(comparison),
                        Err(error) => self.intensity_problem = Some(error.to_string()),
                    },
                    Err(problem) => self.intensity_problem = Some(problem.clone()),
                }
                if let Ok(quant) = &data.reporters {
                    match quant.quantities.quantities.compare(a, b) {
                        Ok(comparison) => self.reporter_comparison = Some(comparison),
                        Err(error) => self.reporter_problem = Some(error.to_string()),
                    }
                }
            }
        }
        self.refresh_table(cx);
    }

    fn refresh_table(&mut self, cx: &mut Context<Self>) {
        self.shown = match self.method {
            Method::Counts => self.comparison.as_ref().map(|c| shown_counts(c, self.test)),
            Method::Intensity => self.intensity_comparison.as_ref().map(shown_intensity),
            Method::Reporters => self.reporter_comparison.as_ref().map(shown_reporters),
        };
        let rows = self
            .shown
            .as_ref()
            .map(|s| s.rows.clone())
            .unwrap_or_default();
        let names = self.pair.clone().unwrap_or_default();
        self.table.update(cx, |table, cx| {
            let delegate = table.delegate_mut();
            delegate.names = names;
            delegate.set_rows(rows);
            table.refresh(cx);
            // Rows changed: a highlight by index would land on another protein.
            table.clear_selection(cx);
        });
        self.all = None;
        self.all_problem = None;
        if self.scope == Scope::All
            && let Some(data) = &self.data
        {
            match shown_all(data, self.method) {
                Ok(all) => self.all = Some(all),
                Err(problem) => self.all_problem = Some(problem),
            }
        }
        self.qc = None;
        self.qc_problem = None;
        if self.scope == Scope::Qc
            && let Some(data) = &self.data
        {
            match shown_qc(data, self.method) {
                Ok(qc) => self.qc = Some(qc),
                Err(problem) => self.qc_problem = Some(problem),
            }
        }
        let rows = self
            .all
            .as_ref()
            .map(|all| all.anova.rows.clone())
            .unwrap_or_default();
        self.anova_table.update(cx, |table, cx| {
            table.delegate_mut().set_rows(rows);
            table.refresh(cx);
            table.clear_selection(cx);
        });
        cx.notify();
    }

    fn set_side(&mut self, baseline: bool, category: String, cx: &mut Context<Self>) {
        let Some((mut a, mut b)) = self.pair.clone() else {
            return; // fewer than two categories: nothing to compare
        };
        // Picking the other side's category swaps the sides.
        if baseline {
            if category == b {
                b = a.clone();
            }
            a = category;
        } else {
            if category == a {
                a = b.clone();
            }
            b = category;
        }
        self.pair = Some((a, b));
        self.failure = None;
        self.compare(true, cx);
    }

    pub fn set_test(&mut self, test: PrimaryTest, cx: &mut Context<Self>) {
        self.test = test;
        self.refresh_table(cx);
    }

    pub fn set_method(&mut self, method: Method, cx: &mut Context<Self>) {
        self.method = method;
        if let Some(data) = &self.data {
            let categories = data.categories_for(method);
            let valid = self
                .pair
                .as_ref()
                .is_some_and(|(a, b)| categories.contains(a) && categories.contains(b));
            if !valid {
                self.pair = data.default_pair_for(method);
                self.compare(false, cx);
                return;
            }
        }
        self.refresh_table(cx);
    }

    pub fn set_scope(&mut self, scope: Scope, cx: &mut Context<Self>) {
        self.scope = scope;
        self.refresh_table(cx);
    }

    pub fn select_plex(&mut self, name: &str, cx: &mut Context<Self>) {
        let Some(plex) = self.data.as_ref().and_then(|d| d.plex_named(name)).cloned() else {
            return;
        };
        self.plex_creating = false;
        self.confirm_remove_plex = false;
        self.plex_selected = Some(plex.name);
        self.plex_reagents = plex.plex;
        self.plex_runs = plex.runs;
        self.plex_references = plex.references;
        self.plex_references_touched = false;
        self.plex_synced = false;
        self.failure = None;
        cx.notify();
    }

    pub fn start_new_plex(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        self.plex_creating = true;
        self.confirm_remove_plex = false;
        self.plex_runs.clear();
        self.plex_references.clear();
        self.plex_references_touched = false;
        self.plex_synced = false;
        self.plex_name_input = cx.new(|cx| InputState::new(window, cx).placeholder("plex name"));
        self.failure = None;
        cx.notify();
    }

    pub fn set_plex_reagents(&mut self, plex: Plex, cx: &mut Context<Self>) {
        if self.plex_reagents != plex {
            self.plex_reagents = plex;
            self.confirm_remove_plex = false;
            cx.notify();
        }
    }

    fn toggle_run(&mut self, run: &str, cx: &mut Context<Self>) {
        self.confirm_remove_plex = false;
        if let Some(at) = self.plex_runs.iter().position(|r| r == run) {
            self.plex_runs.remove(at);
        } else {
            self.plex_runs.push(run.to_string());
        }
        cx.notify();
    }

    /// Ask for the lot's certificate file, then save the plex with it.
    fn load_certificate(&mut self, cx: &mut Context<Self>) {
        self.confirm_remove_plex = false;
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Use certificate".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(file) = picked
                .await
                .ok()
                .and_then(Result::ok)
                .flatten()
                .and_then(|paths| paths.into_iter().next())
            else {
                return;
            };
            let _ = this.update(cx, |this, cx| {
                this.save_plex(Certificate::Load(file), cx);
            });
        })
        .detach();
    }

    /// A save or matching is writing the project (reloads wait for it, and
    /// the workspace holds back its own project-changing actions).
    pub fn is_writing(&self) -> bool {
        self.saving
    }

    /// Ask for the runs' spectra files (or timsTOF `.d` folders), then match
    /// between runs in them.
    fn match_between_runs(&mut self, cx: &mut Context<Self>) {
        if self.busy || self.saving || self.picking {
            return;
        }
        self.picking = true;
        cx.notify();
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: true,
            multiple: true,
            prompt: Some("Match".into()),
        });
        cx.spawn(async move |this, cx| {
            let files = picked
                .await
                .ok()
                .and_then(Result::ok)
                .flatten()
                .filter(|paths| !paths.is_empty());
            let _ = this.update(cx, |this, cx| {
                this.picking = false;
                // A save may have started while the dialog was open.
                match files {
                    Some(files) if !this.busy && !this.saving => this.run_match(files, cx),
                    _ => cx.notify(),
                }
            });
        })
        .detach();
    }

    /// Match between runs off the UI thread, then show the quantities that
    /// include the transfers and a summary (or why it stopped). Like a save,
    /// it holds back reloads while it runs.
    fn run_match(&mut self, files: Vec<PathBuf>, cx: &mut Context<Self>) {
        let path = self.path.clone();
        self.saving = true;
        self.busy = true;
        self.notice = None;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let outcome = cx
                .background_executor()
                .spawn(async move { quant::match_runs(&path, &files) })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.saving = false;
                this.busy = false;
                this.generation += 1;
                this.failure = None;
                let mut problems = Vec::new();
                match outcome.result {
                    Ok(summary) => this.notice = Some(summary),
                    Err(message) => problems.push(message),
                }
                match outcome.data {
                    Ok(data) => this.apply(data, cx),
                    Err(message) => problems.push(message),
                }
                problems.dedup();
                if !problems.is_empty() {
                    this.failure = Some(problems.join("\n"));
                }
                // The matching quantified again itself; a reload queued
                // meanwhile would only clear the message shown here.
                this.reload_after_save = false;
                cx.notify();
            });
        })
        .detach();
    }

    fn toggle_reference(&mut self, channel: &str, cx: &mut Context<Self>) {
        self.confirm_remove_plex = false;
        self.plex_references_touched = true;
        if let Some(at) = self.plex_references.iter().position(|c| c == channel) {
            self.plex_references.remove(at);
        } else {
            self.plex_references.push(channel.to_string());
        }
        cx.notify();
    }

    fn edited_plex(&self) -> Option<&tuffold_store::PlexDefinition> {
        if self.plex_creating {
            return None;
        }
        let data = self.data.as_ref()?;
        self.plex_selected
            .as_deref()
            .and_then(|name| data.plex_named(name))
            .or_else(|| data.plexes.first())
    }

    fn save_plex(&mut self, certificate: Certificate, cx: &mut Context<Self>) {
        self.confirm_remove_plex = false;
        let has_plexes = self.data.as_ref().is_some_and(|d| !d.plexes.is_empty());
        let is_new = self.plex_creating || !has_plexes;
        let typed_name = self.plex_name_input.read(cx).value().to_string();
        let name = if self.plex_creating {
            typed_name
        } else if let Some(selected) = &self.plex_selected {
            selected.clone()
        } else if typed_name.trim().is_empty() {
            "plex".to_string()
        } else {
            typed_name
        };
        let channels: Vec<(String, String)> = self
            .plex_inputs
            .iter()
            .map(|(channel, input)| (channel.clone(), input.read(cx).value().to_string()))
            .collect();
        let references = quant::references_for_save(
            self.edited_plex(),
            self.plex_reagents,
            &channels,
            &self.plex_references,
            self.plex_references_touched,
        );
        let runs = self.plex_runs.clone();
        let plex = self.plex_reagents;
        let path = self.path.clone();
        let selected_after = name.trim().to_string();
        self.saving = true;
        self.busy = true;
        self.notice = None;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move {
                    quant::save_plex_named(
                        &path,
                        &name,
                        is_new,
                        plex,
                        &runs,
                        &channels,
                        references.as_deref(),
                        certificate,
                    )
                })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.saving = false;
                this.busy = false;
                this.generation += 1;
                match result {
                    Ok(data) => {
                        this.failure = None;
                        this.plex_creating = false;
                        this.plex_selected = Some(selected_after);
                        this.apply(data, cx);
                        this.notice = Some("Plex saved.".into());
                    }
                    Err(message) => this.failure = Some(message),
                }
                if std::mem::take(&mut this.reload_after_save) {
                    this.reload(cx);
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn remove_plex(&mut self, cx: &mut Context<Self>) {
        let Some(name) = self.edited_plex().map(|p| p.name.clone()) else {
            return;
        };
        if !self.confirm_remove_plex {
            self.confirm_remove_plex = true;
            cx.notify();
            return;
        }
        self.confirm_remove_plex = false;
        let path = self.path.clone();
        let removed_name = name.clone();
        self.saving = true;
        self.busy = true;
        self.notice = None;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { quant::remove_plex(&path, &name) })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.saving = false;
                this.busy = false;
                this.generation += 1;
                match result {
                    Ok(data) => {
                        this.failure = None;
                        this.plex_creating = false;
                        this.plex_selected = None;
                        this.apply(data, cx);
                        this.notice = Some(format!("Removed plex {removed_name}."));
                    }
                    Err(message) => this.failure = Some(message),
                }
                if std::mem::take(&mut this.reload_after_save) {
                    this.reload(cx);
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn save_groups(&mut self, cx: &mut Context<Self>) {
        let assignments: Vec<(String, String)> = self
            .inputs
            .iter()
            .map(|(name, input)| (name.clone(), input.read(cx).value().to_string()))
            .collect();
        let path = self.path.clone();
        self.saving = true;
        self.busy = true;
        self.notice = None;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { quant::save_categories(&path, &assignments) })
                .await;
            let _ = this.update(cx, |this, cx| {
                // A save always applies, and invalidates any load that may
                // have read the categories before it committed.
                this.saving = false;
                this.busy = false;
                this.generation += 1;
                match result {
                    Ok(data) => {
                        this.failure = None;
                        this.apply(data, cx);
                        this.notice = Some("Sample groups saved.".into());
                    }
                    Err(message) => this.failure = Some(message),
                }
                if std::mem::take(&mut this.reload_after_save) {
                    this.reload(cx);
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn export(&mut self, cx: &mut Context<Self>) {
        if self.exporting || self.data.is_none() {
            return;
        }
        self.exporting = true;
        let directory = self
            .path
            .parent()
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("."));
        let save = cx.prompt_for_new_path(&directory, Some("quantitation.tsv"));
        cx.spawn(async move |this, cx| {
            let target = save.await.ok().and_then(Result::ok).flatten();
            // Snapshot what is shown once the file is chosen, not before.
            let Ok(Some((target, what))) = this.update(cx, |this, _| {
                let Some(target) = target else {
                    this.exporting = false;
                    return None;
                };
                let data = this.data.as_ref()?;
                let what = match (this.scope, this.method) {
                    (Scope::All, _) => Export::Anova(this.all.as_ref()?.anova.clone()),
                    (Scope::Qc, _) => Export::Qc(this.qc.as_ref()?.qc.clone()),
                    (Scope::Pair, Method::Counts) => {
                        Export::Counts(data.counts.clone(), this.comparison.clone())
                    }
                    (Scope::Pair, Method::Intensity) => Export::Intensity(
                        data.intensity.clone().ok()?,
                        this.intensity_comparison.clone(),
                    ),
                    (Scope::Pair, Method::Reporters) => Export::Reporters(
                        data.reporters.clone().ok()?,
                        this.reporter_comparison.clone(),
                    ),
                };
                Some((target, what))
            }) else {
                let _ = this.update(cx, |this, _| this.exporting = false);
                return;
            };
            let shown = target.display().to_string();
            let result = cx
                .background_executor()
                .spawn(async move { quant::export(&target, &what) })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.exporting = false;
                match result {
                    Ok(()) => this.notice = Some(format!("Exported {shown}.")),
                    Err(message) => this.failure = Some(message),
                }
                cx.notify();
            });
        })
        .detach();
    }

    /// Keep one input per biosample (prefilled with its category) and one
    /// input per reporter channel in the selected reagent set.
    fn sync_inputs(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let Some(data) = &self.data else {
            return;
        };
        let quantities = &data.counts;
        let names: Vec<&str> = quantities.samples.iter().map(|s| s.name.as_str()).collect();
        let current: Vec<&str> = self.inputs.iter().map(|(name, _)| name.as_str()).collect();
        if names != current {
            let samples = quantities.samples.clone();
            self.inputs = samples
                .into_iter()
                .map(|sample| {
                    let value = sample.category.unwrap_or_default();
                    let input = cx.new(|cx| {
                        InputState::new(window, cx)
                            .placeholder("category")
                            .default_value(value)
                    });
                    (sample.name, input)
                })
                .collect();
        }
        let wanted: Vec<String> = self
            .plex_reagents
            .channels()
            .into_iter()
            .map(|c| c.name)
            .collect();
        let current_channels: Vec<&str> = self
            .plex_inputs
            .iter()
            .map(|(channel, _)| channel.as_str())
            .collect();
        if !self.plex_synced || wanted.iter().map(String::as_str).ne(current_channels) {
            let typed: std::collections::HashMap<String, String> = if self.plex_synced {
                self.plex_inputs
                    .iter()
                    .map(|(c, input)| (c.clone(), input.read(cx).value().to_string()))
                    .collect()
            } else {
                std::collections::HashMap::new()
            };
            let stored: std::collections::HashMap<String, String> = self
                .edited_plex()
                .map(|p| p.channels.iter().cloned().collect())
                .unwrap_or_default();
            self.plex_inputs = wanted
                .into_iter()
                .map(|channel| {
                    let value = typed
                        .get(&channel)
                        .filter(|v| !v.is_empty())
                        .or_else(|| stored.get(&channel))
                        .cloned()
                        .unwrap_or_default();
                    let input = cx.new(|cx| {
                        InputState::new(window, cx)
                            .placeholder("biosample")
                            .default_value(value)
                    });
                    (channel, input)
                })
                .collect();
            self.plex_synced = true;
        }
    }

    fn toolbar(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let categories = self
            .data
            .as_ref()
            .map(|d| d.categories_for(self.method))
            .unwrap_or_default();
        let (a, b) = self.pair.clone().unwrap_or_default();
        let side = |baseline: bool, chosen: &str, cx: &mut Context<Self>| {
            categories
                .iter()
                .enumerate()
                .map(|(index, category)| {
                    let id = (if baseline { "baseline" } else { "versus" }, index);
                    let button = Button::new(id).xsmall().label(category.clone());
                    let button = button.ghost().selected(category == chosen);
                    let category = category.clone();
                    button.on_click(cx.listener(move |this, _, _, cx| {
                        this.set_side(baseline, category.clone(), cx);
                    }))
                })
                .collect::<Vec<_>>()
        };
        // Welch needs replicates in both categories.
        let replicated = self
            .comparison
            .as_ref()
            .is_some_and(|c| c.primary == PrimaryTest::Welch);
        let tests = [
            (PrimaryTest::Welch, "Welch"),
            (PrimaryTest::Fisher, "Fisher"),
            (PrimaryTest::MannWhitney, "Mann–Whitney"),
        ]
        .map(|(test, label)| {
            let button = Button::new(label)
                .xsmall()
                .label(label)
                .disabled(test != PrimaryTest::Fisher && !replicated);
            let button = button.ghost().selected(test == self.test);
            button.on_click(cx.listener(move |this, _, _, cx| this.set_test(test, cx)))
        });
        let methods = [
            (Method::Counts, "measure-counts", "Spectral counts"),
            (Method::Intensity, "measure-intensity", "MS1 intensity"),
            (Method::Reporters, "measure-reporters", "Reporter ions"),
        ]
        .map(|(method, id, label)| {
            let button = Button::new(id).xsmall().label(label);
            let button = button.ghost().selected(method == self.method);
            button.on_click(cx.listener(move |this, _, _, cx| this.set_method(method, cx)))
        });
        let scopes = [
            (Scope::Pair, "scope-pair", "Two categories"),
            (Scope::All, "scope-all", "All categories"),
            (Scope::Qc, "scope-qc", "Sample QC"),
        ]
        .map(|(scope, id, label)| {
            let button = Button::new(id).xsmall().label(label);
            let button = button.ghost().selected(scope == self.scope);
            button.on_click(cx.listener(move |this, _, _, cx| this.set_scope(scope, cx)))
        });
        let group = |label: &str, children: Vec<Button>| {
            div()
                .flex()
                .items_center()
                .gap_1()
                .child(
                    div()
                        .mr_1()
                        .text_size(px(12.))
                        .text_color(palette.ink_faint)
                        .child(label.to_string()),
                )
                .child(track(palette).children(children))
        };
        div()
            .flex()
            .flex_wrap()
            .items_center()
            .gap_6()
            .child(group("Compare", scopes.into_iter().collect()))
            .children(
                (self.scope == Scope::Pair && categories.len() >= 2)
                    .then(|| group("Baseline", side(true, &a, cx))),
            )
            .children(
                (self.scope == Scope::Pair && categories.len() >= 2)
                    .then(|| group("Versus", side(false, &b, cx))),
            )
            .child(group("Measure", methods.into_iter().collect()))
            .children(
                (self.scope == Scope::Pair && self.method == Method::Counts)
                    .then(|| group("Test", tests.into_iter().collect())),
            )
            .child(div().flex_1())
            .children((self.method == Method::Intensity).then(|| {
                // Matching needs a validated project (any loaded data is) and
                // MS1 data to put the transfers beside.
                Button::new("match-between-runs")
                    .ghost()
                    .small()
                    .label("Match between runs…")
                    .disabled(
                        self.busy
                            || self.picking
                            || self.exporting
                            || self
                                .data
                                .as_ref()
                                .is_none_or(|data| data.intensity.is_err()),
                    )
                    .on_click(cx.listener(|this, _, _, cx| this.match_between_runs(cx)))
            }))
            .child(
                Button::new("export-quant")
                    .ghost()
                    .small()
                    .label("Export TSV…")
                    .disabled(
                        self.data.is_none()
                            || self.busy
                            || self.exporting
                            || (self.scope == Scope::All && self.all.is_none())
                            || (self.scope == Scope::Qc && self.qc.is_none())
                            || (self.scope == Scope::Pair
                                && self.method == Method::Intensity
                                && self.intensity_comparison.is_none())
                            || (self.scope == Scope::Pair
                                && self.method == Method::Reporters
                                && self.reporter_comparison.is_none()),
                    )
                    .on_click(cx.listener(|this, _, _, cx| this.export(cx))),
            )
    }

    fn plex_panel(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let plexes = self
            .data
            .as_ref()
            .map(|d| d.plexes.as_slice())
            .unwrap_or_default();
        let edited = self.edited_plex();
        let selector = (!plexes.is_empty()).then(|| {
            let mut buttons: Vec<Button> = plexes
                .iter()
                .enumerate()
                .map(|(index, plex)| {
                    let active = !self.plex_creating && edited.is_some_and(|e| e.name == plex.name);
                    let button = Button::new(("select-plex", index))
                        .xsmall()
                        .label(plex.name.clone())
                        .disabled(self.busy);
                    let button = button.ghost().selected(active);
                    let name = plex.name.clone();
                    button.on_click(cx.listener(move |this, _, _, cx| {
                        this.select_plex(&name, cx);
                    }))
                })
                .collect();
            let new_button = Button::new("new-plex")
                .xsmall()
                .label("New plex")
                .disabled(self.busy);
            let new_button = new_button.ghost().selected(self.plex_creating);
            buttons.push(new_button.on_click(cx.listener(|this, _, window, cx| {
                this.start_new_plex(window, cx);
            })));
            div()
                .flex()
                .flex_wrap()
                .items_center()
                .gap_1()
                .child(
                    div()
                        .mr_1()
                        .text_size(px(12.))
                        .text_color(palette.ink_faint)
                        .child("Plex"),
                )
                .children(buttons)
        });
        let name_row = self.plex_creating.then(|| {
            div()
                .flex()
                .items_center()
                .gap_2()
                .child(
                    div()
                        .text_size(px(12.))
                        .text_color(palette.ink_faint)
                        .child("Name"),
                )
                .child(
                    div()
                        .w(px(180.))
                        .child(Input::new(&self.plex_name_input).small()),
                )
        });
        let reagents = Plex::ALL.map(|plex| {
            let button = Button::new(plex.name()).xsmall().label(plex_title(plex));
            let button = button.ghost().selected(plex == self.plex_reagents);
            button.on_click(cx.listener(move |this, _, _, cx| this.set_plex_reagents(plex, cx)))
        });
        let runs = self
            .data
            .as_ref()
            .map(|d| d.runs_for_plex(edited.map(|p| p.name.as_str())))
            .unwrap_or_default();
        let first_plex_only =
            plexes.is_empty() || (plexes.len() == 1 && !self.plex_creating && edited.is_some());
        let runs_row = (!runs.is_empty()).then(|| {
            let buttons = runs.iter().enumerate().map(|(index, run)| {
                let selected = self.plex_runs.contains(&run.name);
                let label = match &run.in_other_plex {
                    Some(other) => format!("{} ({other})", run.name),
                    None => run.name.clone(),
                };
                let button = Button::new(("plex-run", index))
                    .xsmall()
                    .label(label)
                    .disabled(self.busy || run.in_other_plex.is_some());
                let button = button
                    .ghost()
                    .selected(selected && run.in_other_plex.is_none());
                let name = run.name.clone();
                button.on_click(cx.listener(move |this, _, _, cx| {
                    this.toggle_run(&name, cx);
                }))
            });
            let note = if first_plex_only {
                "For the first plex only, leaving every run unselected uses all runs."
            } else {
                "Pick the runs of this pool; runs in another plex are disabled."
            };
            div()
                .flex()
                .flex_col()
                .gap_1()
                .child(
                    div()
                        .flex()
                        .flex_wrap()
                        .items_center()
                        .gap_1()
                        .child(
                            div()
                                .mr_1()
                                .text_size(px(12.))
                                .text_color(palette.ink_faint)
                                .child("Runs"),
                        )
                        .children(buttons),
                )
                .child(
                    div()
                        .text_size(px(12.))
                        .text_color(palette.ink_faint)
                        .child(note),
                )
        });
        let cells: Vec<(usize, &(String, Entity<InputState>))> =
            self.plex_inputs.iter().enumerate().collect();
        let rows = cells.chunks(3).map(|chunk| {
            div()
                .flex()
                .items_center()
                .gap_3()
                .children(chunk.iter().map(|&(index, (channel, input))| {
                    let reference = self.plex_references.contains(channel);
                    let toggle = Button::new(("plex-reference", index))
                        .xsmall()
                        .label("ref")
                        .disabled(self.busy);
                    let toggle = toggle.ghost().selected(reference);
                    let channel_name = channel.clone();
                    div()
                        .flex_1()
                        .min_w_0()
                        .flex()
                        .items_center()
                        .gap_1p5()
                        .child(
                            div()
                                .w(px(38.))
                                .font_family(MONO)
                                .text_size(px(12.))
                                .text_color(palette.ink_muted)
                                .child(channel.clone()),
                        )
                        .child(div().flex_1().min_w_0().child(Input::new(input).small()))
                        .child(toggle.on_click(cx.listener(move |this, _, _, cx| {
                            this.toggle_reference(&channel_name, cx);
                        })))
                }))
        });
        let remove_button = edited.map(|_| {
            let button = Button::new("remove-plex")
                .small()
                .label(if self.confirm_remove_plex {
                    "Confirm remove"
                } else {
                    "Remove plex"
                })
                .disabled(self.busy);
            let button = if self.confirm_remove_plex {
                button.danger()
            } else {
                button.ghost()
            };
            button.on_click(cx.listener(|this, _, _, cx| this.remove_plex(cx)))
        });
        panel_shell(palette)
            .child(
                div()
                    .flex()
                    .items_center()
                    .justify_between()
                    .gap_2()
                    .child(panel_title("Plex"))
                    .child(
                        div()
                            .flex()
                            .items_center()
                            .gap_2()
                            .children(remove_button)
                            .child(
                                Button::new("save-plex")
                                    .primary()
                                    .small()
                                    .label("Save plex")
                                    .disabled(self.busy || self.plex_inputs.is_empty())
                                    .on_click(cx.listener(|this, _, _, cx| {
                                        this.save_plex(Certificate::Keep, cx);
                                    })),
                            ),
                    ),
            )
            .child(
                div()
                    .text_size(px(12.))
                    .text_color(palette.ink_faint)
                    .child(
                        "Choose the reagent set and name the biosample in each reporter channel; \
                         leave unused channels empty. Mark the pooled reference channel (ref) \
                         to combine plexes.",
                    ),
            )
            .child(
                div()
                    .flex()
                    .flex_col()
                    .gap_1p5()
                    .children(selector)
                    .children(name_row)
                    .child(
                        div()
                            .flex()
                            .flex_wrap()
                            .items_center()
                            .gap_1()
                            .child(
                                div()
                                    .mr_1()
                                    .text_size(px(12.))
                                    .text_color(palette.ink_faint)
                                    .child("Reagents"),
                            )
                            .children(reagents),
                    )
                    .children(runs_row),
            )
            .child(div().flex().flex_col().gap_2().children(rows))
            .child(self.certificate_row(palette, cx))
    }

    /// Whether reporters are corrected for the lot's isotope impurities, with
    /// buttons to load a certificate or drop the plex's.
    fn certificate_row(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let saved = self.edited_plex();
        let corrected =
            saved.is_some_and(|p| p.impurities.is_some() && p.plex == self.plex_reagents);
        let text = if corrected {
            "Isotope impurities: corrected with the lot's certificate."
        } else {
            "Isotope impurities: not corrected. Load the certificate of your reagent lot (a text \
             file, one line per channel: its name, then 4 percentages for −2/−1/+1/+2 Da, or 8 \
             for TMTpro) to correct them."
        };
        let disabled = self.busy || self.plex_inputs.is_empty();
        div()
            .flex()
            .items_center()
            .gap_3()
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .text_size(px(12.))
                    .text_color(palette.ink_faint)
                    .child(text),
            )
            .children(corrected.then(|| {
                Button::new("remove-certificate")
                    .ghost()
                    .small()
                    .label("Remove")
                    .disabled(disabled)
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.save_plex(Certificate::Remove, cx);
                    }))
            }))
            .child(
                Button::new("load-certificate")
                    .small()
                    .label("Certificate…")
                    .disabled(disabled)
                    .on_click(cx.listener(|this, _, _, cx| this.load_certificate(cx))),
            )
    }

    fn groups_panel(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let totals = self
            .data
            .as_ref()
            .map(|d| d.counts.totals.clone())
            .unwrap_or_default();
        let several_plexes = self.data.as_ref().is_some_and(|d| d.plexes.len() > 1);
        let channel_of: std::collections::HashMap<&str, String> = self
            .data
            .as_ref()
            .map(|d| {
                d.plexes
                    .iter()
                    .flat_map(|p| {
                        p.channels.iter().map(move |(c, s)| {
                            let label = if several_plexes {
                                format!("{} {c}", p.name)
                            } else {
                                format!("channel {c}")
                            };
                            (s.as_str(), label)
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();
        let only_plex_samples = self.method == Method::Reporters && !channel_of.is_empty();
        let rows = self
            .inputs
            .iter()
            .enumerate()
            .filter(|(_, (name, _))| !only_plex_samples || channel_of.contains_key(name.as_str()))
            .map(|(index, (name, input))| {
                let spectra = totals.get(index).copied().unwrap_or(0);
                let note = match (channel_of.get(name.as_str()), spectra) {
                    (Some(channel), 0) => channel.clone(),
                    (Some(channel), n) if self.method == Method::Reporters => {
                        format!("{channel} · {} spectra", group_digits(n))
                    }
                    (_, n) => format!("{} spectra", group_digits(n)),
                };
                div()
                    .flex()
                    .items_center()
                    .gap_3()
                    .child(
                        div()
                            .w(px(150.))
                            .font_family(MONO)
                            .text_size(px(13.))
                            .truncate()
                            .child(name.clone()),
                    )
                    .child(div().w(px(180.)).child(Input::new(input).small()))
                    .child(
                        div()
                            .flex_1()
                            .min_w_0()
                            .truncate()
                            .text_size(px(12.))
                            .text_color(palette.ink_faint)
                            .child(note),
                    )
            });
        panel_shell(palette)
            .child(
                div()
                    .flex()
                    .items_center()
                    .justify_between()
                    .child(panel_title("Sample groups"))
                    .child(
                        Button::new("save-groups")
                            .primary()
                            .small()
                            .label("Save groups")
                            .disabled(self.busy || self.inputs.is_empty())
                            .on_click(cx.listener(|this, _, _, cx| this.save_groups(cx))),
                    ),
            )
            .child(div().flex().flex_col().gap_2().children(rows))
    }
}

/// x ticks at whole log2 folds, y ticks at whole -log10 q.
fn whole_ticks(max: f64) -> Vec<f64> {
    let step = if max > 12.0 {
        4.0
    } else if max > 6.0 {
        2.0
    } else {
        1.0
    };
    (0..)
        .map(|i| f64::from(i) * step)
        .take_while(|t| *t <= max)
        .collect()
}

fn volcano_view(volcano: &Volcano, palette: &Palette) -> Div {
    let (x_max, y_max) = (volcano.x_max, volcano.y_max);
    let x_of = move |fold: f64| ((fold + x_max) / (2.0 * x_max)) as f32;
    let y_of = move |height: f64| (1.0 - height / y_max) as f32;
    let dots: Vec<(f32, f32, Hsla)> = volcano
        .points
        .iter()
        .map(|p| {
            let color = match (p.significant, p.log2_fold >= 0.0) {
                (false, _) => palette.rule_strong,
                (true, true) => palette.signal,
                (true, false) => palette.ion_y,
            };
            (x_of(p.log2_fold), y_of(p.height), color)
        })
        .collect();
    let (axis, threshold_y) = (palette.rule_strong, y_of(volcano.threshold));
    let labels = volcano.points.iter().filter(|p| p.labelled).map(|p| {
        div()
            .absolute()
            .left(relative(x_of(p.log2_fold)))
            .top(px(PLOT_HEIGHT * y_of(p.height) - 18.))
            .ml(px(6.))
            .font_family(MONO)
            .text_size(px(10.))
            .text_color(palette.ink_muted)
            .child(
                p.accession
                    .split('|')
                    .nth(1)
                    .unwrap_or(&p.accession)
                    .to_string(),
            )
    });
    let x_ticks = whole_ticks(x_max)
        .into_iter()
        .flat_map(|t| if t == 0.0 { vec![0.0] } else { vec![-t, t] })
        .map(|t| {
            div()
                .absolute()
                .left(relative(x_of(t)))
                .top(px(PLOT_HEIGHT + 4.))
                .ml(px(-6.))
                .font_family(MONO)
                .text_size(px(11.))
                .text_color(palette.ink_faint)
                .child(format!("{t:.0}"))
        });
    let y_ticks = whole_ticks(y_max).into_iter().skip(1).map(|t| {
        div()
            .absolute()
            .left(px(4.))
            .top(px(PLOT_HEIGHT * y_of(t) - 7.))
            .font_family(MONO)
            .text_size(px(11.))
            .text_color(palette.ink_faint)
            .child(format!("{t:.0}"))
    });
    // Names the dashed line where the fewest points sit on it.
    let threshold_label = div()
        .absolute()
        .left(relative(volcano.threshold_label_x as f32))
        .top(px(PLOT_HEIGHT * threshold_y - 16.))
        .font_family(MONO)
        .text_size(px(10.))
        .text_color(palette.ink_faint)
        .child(format!("q = {:.0}%", SIGNIFICANCE * 100.0));
    div()
        .relative()
        .w_full()
        .h(px(PLOT_HEIGHT + 22.))
        .child(
            canvas(
                |_, _, _| {},
                move |bounds: Bounds<_>, (), window, _| {
                    let (left, top, width) = (bounds.left(), bounds.top(), bounds.size.width);
                    let height = px(PLOT_HEIGHT);
                    window.paint_quad(fill(
                        Bounds::new(point(left, top + height), size(width, px(1.))),
                        axis,
                    ));
                    window.paint_quad(fill(
                        Bounds::new(point(left + width * 0.5, top), size(px(1.), height)),
                        axis,
                    ));
                    // The significance threshold, dashed.
                    let y = top + height * threshold_y;
                    let mut x = px(0.);
                    while x < width {
                        window.paint_quad(fill(
                            Bounds::new(point(left + x, y), size(px(6.), px(1.))),
                            axis,
                        ));
                        x += px(10.);
                    }
                    for (dx, dy, color) in &dots {
                        let center = point(left + width * *dx, top + height * *dy);
                        let radius = px(3.5);
                        window.paint_quad(
                            fill(
                                Bounds::new(
                                    point(center.x - radius, center.y - radius),
                                    size(radius * 2., radius * 2.),
                                ),
                                *color,
                            )
                            .corner_radii(radius),
                        );
                    }
                },
            )
            .absolute()
            .top_0()
            .left_0()
            .w_full()
            .h(px(PLOT_HEIGHT)),
        )
        .children(labels)
        .children(x_ticks)
        .children(y_ticks)
        .child(threshold_label)
}

impl QuantView {
    /// Every category at once: summary, heat map, and the ANOVA table.
    fn all_body(&self, all: &AllShown, palette: &Palette) -> AnyElement {
        let tested = all.anova.rows.iter().filter(|r| r.p.is_some()).count();
        let significant = all
            .anova
            .rows
            .iter()
            .filter(|r| r.q.is_some_and(|q| q <= SIGNIFICANCE))
            .count();
        let summary = format!(
            "{significant} of {tested} tested protein groups differ at {}% FDR ({} not testable) · \
             {}",
            SIGNIFICANCE * 100.0,
            all.anova.rows.len() - tested,
            all.test,
        );
        div()
            .flex()
            .flex_col()
            .gap_3()
            .child(
                div()
                    .text_size(px(12.))
                    .text_color(palette.ink_muted)
                    .child(summary),
            )
            .child(
                div()
                    .flex()
                    .gap_4()
                    .child(
                        panel_shell(palette)
                            .flex_1()
                            .min_w_0()
                            .child(
                                div()
                                    .text_size(px(12.))
                                    .text_color(palette.ink_faint)
                                    .child(all.caption.clone()),
                            )
                            .child(heatmap_view(&all.heatmap, palette)),
                    )
                    .child(
                        div()
                            .w(px(460.))
                            .flex_none()
                            .h(px(PLOT_HEIGHT + 90.))
                            .flex()
                            .flex_col()
                            .gap_2()
                            .child(
                                div()
                                    .flex_1()
                                    .min_h_0()
                                    .child(DataTable::new(&self.anova_table).small()),
                            ),
                    ),
            )
            .into_any_element()
    }
}

/// Row labels need this much height per row to be drawn.
const LABEL_ROW: f32 = 9.;

fn heatmap_view(map: &Heatmap, palette: &Palette) -> Div {
    if map.rows.is_empty() {
        return div()
            .text_size(px(13.))
            .text_color(palette.ink_muted)
            .child("Nothing to draw: no protein group could be tested.");
    }
    let row_height = (PLOT_HEIGHT / map.rows.len() as f32).clamp(4., 18.);
    let height = row_height * map.rows.len() as f32;
    let columns = map.columns.len().max(1);
    let label_width = if row_height >= LABEL_ROW { 96. } else { 0. };
    let clamp = crate::theme::heat_clamp(map.columns.len());
    let cells: Vec<(usize, usize, Hsla)> = map
        .rows
        .iter()
        .enumerate()
        .flat_map(|(r, row)| {
            row.z.iter().enumerate().map(move |(c, z)| {
                (
                    r,
                    c,
                    z.map_or(palette.paper_sunken, |z| palette.heat(z, clamp)),
                )
            })
        })
        .collect();
    // Category bands above the columns: a label where each category starts.
    let mut bands: Vec<(usize, usize, String)> = Vec::new();
    for (c, category) in map.column_categories.iter().enumerate() {
        let name = category.clone().unwrap_or_else(|| "(none)".into());
        match bands.last_mut() {
            Some((_, end, last)) if *last == name => *end = c + 1,
            _ => bands.push((c, c + 1, name)),
        }
    }
    let band_labels = bands.into_iter().map(|(start, end, name)| {
        div()
            .absolute()
            .left(relative(start as f32 / columns as f32))
            .w(relative((end - start) as f32 / columns as f32))
            .top_0()
            .px_1()
            .border_b_1()
            .border_color(palette.rule_strong)
            .text_size(px(11.))
            .text_color(palette.ink_muted)
            .truncate()
            .child(name)
    });
    let column_labels = map.columns.iter().enumerate().map(|(c, name)| {
        div()
            .absolute()
            .left(relative(c as f32 / columns as f32))
            .w(relative(1. / columns as f32))
            .top(px(height + 4.))
            .flex()
            .justify_center()
            .font_family(MONO)
            .text_size(px(10.))
            .text_color(palette.ink_faint)
            .truncate()
            .child(name.clone())
    });
    let row_labels = (label_width > 0.).then(|| {
        map.rows
            .iter()
            .enumerate()
            .map(|(r, row)| {
                div()
                    .absolute()
                    .left_0()
                    .top(px(r as f32 * row_height))
                    .h(px(row_height))
                    .flex()
                    .items_center()
                    .font_family(MONO)
                    .text_size(px(10.))
                    .text_color(palette.ink_muted)
                    .child(
                        row.accession
                            .split('|')
                            .nth(1)
                            .unwrap_or(&row.accession)
                            .to_string(),
                    )
            })
            .collect::<Vec<_>>()
    });
    div()
        .flex()
        .gap_2()
        .children(row_labels.map(|labels| {
            div()
                .relative()
                .flex_none()
                .w(px(label_width))
                .mt(px(20.))
                .h(px(height))
                .children(labels)
        }))
        .child(
            div()
                .relative()
                .flex_1()
                .min_w_0()
                .h(px(height + 40.))
                .children(band_labels)
                .child(
                    div()
                        .absolute()
                        .top(px(20.))
                        .left_0()
                        .w_full()
                        .h(px(height))
                        .children(column_labels)
                        .child(
                            canvas(
                                |_, _, _| {},
                                move |bounds: Bounds<_>, (), window, _| {
                                    let width = bounds.size.width / columns as f32;
                                    let gap = px(1.);
                                    for (r, c, color) in &cells {
                                        let origin = point(
                                            bounds.left() + width * *c as f32,
                                            bounds.top() + px(row_height * *r as f32),
                                        );
                                        window.paint_quad(fill(
                                            Bounds::new(
                                                origin,
                                                size(
                                                    width - gap,
                                                    px(row_height) - gap.min(px(row_height / 4.)),
                                                ),
                                            ),
                                            *color,
                                        ));
                                    }
                                },
                            )
                            .absolute()
                            .top_0()
                            .left_0()
                            .size_full(),
                        ),
                ),
        )
}

impl Render for QuantView {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        self.sync_inputs(window, cx);
        let palette = *cx.global::<Palette>();
        let categories = self
            .data
            .as_ref()
            .map(|d| d.categories_for(self.method))
            .unwrap_or_default();
        let has_plex = self.data.as_ref().is_some_and(|d| !d.plexes.is_empty());
        let messages = div()
            .flex()
            .flex_col()
            .gap_2()
            .children(self.failure.clone().map(|message| {
                callout(palette.fault, palette.fault_soft)
                    .text_size(px(13.))
                    .text_color(palette.fault)
                    .child(message)
            }))
            .children(self.notice.clone().map(|message| {
                div()
                    .text_size(px(13.))
                    .text_color(palette.ok)
                    .child(message)
            }));
        let status = match self.method {
            Method::Counts => None,
            Method::Intensity => self
                .data
                .as_ref()
                .and_then(|d| d.intensity.as_ref().ok())
                .map(ms1_status),
            Method::Reporters => self
                .data
                .as_ref()
                .and_then(|d| d.reporters.as_ref().ok())
                .map(reporter_status),
        }
        .map(|status| {
            div()
                .text_size(px(12.))
                .text_color(palette.ink_faint)
                .child(status)
        });
        let reagent_hint = (self.method == Method::Reporters)
            .then(|| {
                self.data
                    .as_ref()
                    .and_then(|d| d.reporters.as_ref().ok())
                    .and_then(reporter_reagent_hint)
            })
            .flatten()
            .map(|hint| {
                callout(palette.caution, palette.caution_soft)
                    .text_size(px(13.))
                    .text_color(palette.caution)
                    .child(hint)
            });
        let problem = match (self.scope, self.method) {
            (Scope::All, _) => self.all_problem.clone(),
            (Scope::Qc, _) => self.qc_problem.clone(),
            (Scope::Pair, Method::Counts) => self.pair_problem.clone(),
            (Scope::Pair, Method::Intensity) => (categories.len() >= 2)
                .then(|| self.intensity_problem.clone())
                .flatten(),
            (Scope::Pair, Method::Reporters) => (!has_plex || categories.len() >= 2)
                .then(|| self.reporter_problem.clone())
                .flatten(),
        }
        .map(|message| {
            callout(palette.caution, palette.caution_soft)
                .text_size(px(13.))
                .text_color(palette.caution)
                .child(message)
        });
        let body = match (&self.scope, &self.all, &self.shown) {
            (Scope::All, Some(all), _) => self.all_body(all, &palette),
            (Scope::All, None, _) => div().into_any_element(),
            (Scope::Qc, _, _) => match &self.qc {
                Some(qc) => qc_body(qc, &palette).into_any_element(),
                None => div().into_any_element(),
            },
            (Scope::Pair, _, shown) => match shown {
                Some(shown) => {
                    let plot = volcano(shown);
                    let caveat = match (self.method, self.test, self.comparison.as_ref()) {
                        (Method::Counts, PrimaryTest::Fisher, Some(c))
                            if c.primary == PrimaryTest::Welch =>
                        {
                            Some(
                                "Fisher overstates significance with replicates; Welch is the primary test \
                                 here."
                                    .to_string(),
                            )
                        }
                        (Method::Counts, PrimaryTest::MannWhitney, Some(c)) => {
                            mann_whitney_caveat(c.samples_a.len(), c.samples_b.len())
                        }
                        _ => None,
                    }
                    .map(|text| {
                        div()
                            .text_size(px(12.))
                            .text_color(palette.caution)
                            .child(text)
                    });
                    let untested = match (plot.untested, self.method) {
                        (0, _) => String::new(),
                        (n, Method::Intensity | Method::Reporters) => {
                            format!(" ({n} not testable: they need two values in each category)")
                        }
                        (n, Method::Counts) => format!(" ({n} not testable)"),
                    };
                    let summary = format!(
                        "{} of {} protein groups differ at {}% FDR{untested} · {}",
                        plot.significant,
                        shown.rows.len(),
                        SIGNIFICANCE * 100.0,
                        shown.test,
                    );
                    div()
                        .flex()
                        .flex_col()
                        .gap_3()
                        .child(
                            div()
                                .text_size(px(12.))
                                .text_color(palette.ink_muted)
                                .child(summary),
                        )
                        .children(caveat)
                        .child(
                            div()
                                .flex()
                                .gap_4()
                                .child(
                                    panel_shell(&palette)
                                        .flex_1()
                                        .min_w_0()
                                        .child(volcano_view(&plot, &palette)),
                                )
                                .child(
                                    div()
                                        .w(px(500.))
                                        .flex_none()
                                        .h(px(PLOT_HEIGHT + 90.))
                                        .flex()
                                        .flex_col()
                                        .gap_2()
                                        .child(
                                            div()
                                                .text_size(px(12.))
                                                .text_color(palette.ink_faint)
                                                .child(shown.means),
                                        )
                                        .child(
                                            div()
                                                .flex_1()
                                                .min_h_0()
                                                .child(DataTable::new(&self.table).small()),
                                        ),
                                ),
                        )
                        .into_any_element()
                }
                None if self.method == Method::Reporters && !has_plex => div().into_any_element(),
                None if categories.len() < 2 => div()
                    .text_size(px(14.))
                    .text_color(palette.ink_muted)
                    .child(
                        "Comparisons need biosamples in at least two categories. Name a category \
                         for each biosample below (for example control and treated) and save.",
                    )
                    .into_any_element(),
                None => div().into_any_element(),
            },
        };
        let panels = if self.method == Method::Reporters {
            div()
                .flex()
                .gap_4()
                .items_start()
                .child(self.plex_panel(&palette, cx).flex_1().min_w_0())
                .child(self.groups_panel(&palette, cx).flex_1().min_w_0())
                .into_any_element()
        } else {
            self.groups_panel(&palette, cx).into_any_element()
        };
        div()
            .id("quant")
            .size_full()
            .overflow_y_scroll()
            .p_4()
            .flex()
            .flex_col()
            .gap_4()
            .child(self.toolbar(&palette, cx))
            .child(messages)
            .children(status)
            .children(reagent_hint)
            .children(problem)
            .child(body)
            .child(panels)
    }
}

/// Mann–Whitney's limit with few replicates, when it rules out any call.
fn mann_whitney_caveat(a: usize, b: usize) -> Option<String> {
    let floor = tuffold_quant::stats::mann_whitney_floor(a, b);
    (floor > SIGNIFICANCE).then(|| {
        format!(
            "With {a} and {b} samples the smallest p Mann–Whitney can give is {floor:.2}, so \
             nothing can pass {:.0}% FDR; it needs more replicates (4 and 4 or more).",
            SIGNIFICANCE * 100.0
        )
    })
}
