//! The linked browsing views. Selecting a protein group loads its peptides,
//! selecting a peptide its spectra, selecting a spectrum its annotated peaks;
//! each step reads the project on a background thread.
//!
//! Loads are ordered by level (list, group, peptide, spectrum). Starting one
//! cancels any pending load at its level and below, never above, so a quick
//! click on a spectrum cannot drop the peptide load it depends on, and a
//! reload is never dropped by a selection. Each level remembers the id it is
//! showing or loading; a click on that id again is ignored.

use std::path::PathBuf;

use gpui_kit::component::Selectable as _;
use gpui_kit::component::Sizable as _;
use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::component::input::{Input, InputEvent, InputState};
use gpui_kit::component::select::{SearchableVec, Select, SelectEvent, SelectItem, SelectState};
use gpui_kit::component::table::{DataTable, TableEvent, TableState};
use gpui_kit::{
    AppContext as _, Context, Div, Entity, EventEmitter, FontWeight, InteractiveElement as _,
    IntoElement, ParentElement, Render, SharedString, StatefulInteractiveElement as _, Styled,
    Subscription, Window, div, px, relative,
};
use tuffold_core::fragments::Tolerance;
use tuffold_store::{GroupDetail, ProteinListRow, SpectrumDetail};
use tuffold_validate::{ThresholdMode, format_rate};

use super::spectrum_plot::{coverage_view, ladder_view, spectrum_plot_view};
use super::tables::{
    PeptideTable, ProteinLayout, ProteinTable, Sortable as _, SpectrumTable, probability_text,
};
use super::{callout, eyebrow, group_digits, panel_shell, panel_title, track};
use crate::browse::{
    self, BrowserData, COVERAGE_LINE, Drill, NoPlot, ProteinFilter, SampleMeasure, TOLERANCES,
    coverage_lines, facets, highest_value, saf_totals, spectrum_plot, tolerance_label,
};
use crate::theme::{MONO, Palette};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrowserMode {
    Proteins,
    Samples,
}

/// The user asked to see the selected group in the Proteins view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShowProtein;

/// What a background load fills, outermost first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Level {
    List = 0,
    Group = 1,
    Peptide = 2,
    Spectrum = 3,
}

/// One choice in a filter menu: its label shows how many groups have it.
#[derive(Debug, Clone, PartialEq)]
struct Choice {
    label: SharedString,
    value: String,
}

impl SelectItem for Choice {
    type Value = String;

    fn title(&self) -> SharedString {
        self.label.clone()
    }

    fn value(&self) -> &String {
        &self.value
    }
}

/// A filter menu's choices; searchable menus filter them as the user types.
type Menu = SearchableVec<Choice>;

fn choices(counted: Vec<(String, usize)>) -> Menu {
    counted
        .into_iter()
        .map(|(value, n)| Choice {
            label: format!("{value} · {}", group_digits(n)).into(),
            value,
        })
        .collect::<Vec<_>>()
        .into()
}

pub struct Browser {
    path: PathBuf,
    pub mode: BrowserMode,
    proteins: Entity<TableState<ProteinTable>>,
    samples: Entity<TableState<ProteinTable>>,
    peptides: Entity<TableState<PeptideTable>>,
    spectra: Entity<TableState<SpectrumTable>>,
    drill: Drill,
    tolerance: Tolerance,
    failure: Option<String>,
    /// Per [`Level`]: bumped by a load at that level or above; a finished
    /// load applies only if its level's generation is unchanged.
    generations: [u64; 4],
    /// Per [`Level`]: whether a load is in flight.
    pending: [bool; 4],
    /// The group, peptide, and spectrum shown or being loaded.
    wanted: [Option<i64>; 3],
    /// Every group passing the thresholds; the tables list those `filter`
    /// keeps.
    all_rows: Vec<ProteinListRow>,
    filter: ProteinFilter,
    search: Entity<InputState>,
    min_spectra: Entity<InputState>,
    organisms: Entity<SelectState<Menu>>,
    modifications: Entity<SelectState<Menu>>,
    go_terms: Entity<SelectState<Menu>>,
    /// The menus' choices changed since they were last given to the menus
    /// (which needs a window, so happens on render).
    facets_changed: bool,
    /// Whether any listed group names an organism / carries a modification
    /// (the menus are left out otherwise).
    has_organisms: bool,
    has_modifications: bool,
    has_go: bool,
    _subscriptions: Vec<Subscription>,
}

impl EventEmitter<ShowProtein> for Browser {}

impl Browser {
    /// A browser over `path`; `data` shows immediately (snapshots), otherwise
    /// it loads in the background.
    pub fn new(
        path: PathBuf,
        mode: BrowserMode,
        data: Option<BrowserData>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let proteins =
            cx.new(|cx| TableState::new(ProteinTable::new(ProteinLayout::List), window, cx));
        let samples =
            cx.new(|cx| TableState::new(ProteinTable::new(ProteinLayout::Samples), window, cx));
        let peptides = cx.new(|cx| TableState::new(PeptideTable::default(), window, cx));
        let spectra = cx.new(|cx| TableState::new(SpectrumTable::default(), window, cx));
        let search =
            cx.new(|cx| InputState::new(window, cx).placeholder("Accession or description"));
        let min_spectra = cx.new(|cx| InputState::new(window, cx).placeholder("0"));
        let organisms = cx
            .new(|cx| SelectState::new(Menu::from(Vec::new()), None, window, cx).searchable(true));
        let modifications = cx.new(|cx| SelectState::new(Menu::from(Vec::new()), None, window, cx));
        let go_terms = cx
            .new(|cx| SelectState::new(Menu::from(Vec::new()), None, window, cx).searchable(true));
        let subscriptions = vec![
            cx.subscribe(&search, |this, input, event: &InputEvent, cx| {
                if matches!(event, InputEvent::Change) {
                    this.filter.text = input.read(cx).value().to_string();
                    this.refilter(cx);
                }
            }),
            cx.subscribe(&min_spectra, |this, input, event: &InputEvent, cx| {
                if matches!(event, InputEvent::Change) {
                    // Anything that isn't a whole number counts as no minimum.
                    this.filter.min_spectra = input.read(cx).value().trim().parse().unwrap_or(0);
                    this.refilter(cx);
                }
            }),
            cx.subscribe(&organisms, |this, _, event: &SelectEvent<Menu>, cx| {
                let SelectEvent::Confirm(value) = event;
                this.filter.organism = value.clone();
                this.refilter(cx);
            }),
            cx.subscribe(&go_terms, |this, _, event: &SelectEvent<Menu>, cx| {
                let SelectEvent::Confirm(value) = event;
                this.filter.go = value.clone();
                this.refilter(cx);
            }),
            cx.subscribe(&modifications, |this, _, event: &SelectEvent<Menu>, cx| {
                let SelectEvent::Confirm(value) = event;
                this.filter.modification = value.clone();
                this.refilter(cx);
            }),
            cx.subscribe_in(
                &proteins,
                window,
                |this, table, event: &TableEvent, _, cx| {
                    if let TableEvent::SelectRow(row) = event {
                        let id = table.read(cx).delegate().rows.get(*row).map(|r| r.group_id);
                        this.select_group(id, cx);
                    }
                },
            ),
            cx.subscribe_in(
                &samples,
                window,
                |this, table, event: &TableEvent, _, cx| {
                    let row = match event {
                        TableEvent::SelectRow(row) | TableEvent::DoubleClickedRow(row) => *row,
                        _ => return,
                    };
                    let id = table.read(cx).delegate().rows.get(row).map(|r| r.group_id);
                    this.select_group(id, cx);
                    if matches!(event, TableEvent::DoubleClickedRow(_)) {
                        cx.emit(ShowProtein);
                    }
                },
            ),
            cx.subscribe_in(
                &peptides,
                window,
                |this, table, event: &TableEvent, _, cx| {
                    if let TableEvent::SelectRow(row) = event {
                        let id = table
                            .read(cx)
                            .delegate()
                            .rows
                            .get(*row)
                            .map(|r| r.peptide_id);
                        this.select_peptide(id, cx);
                    }
                },
            ),
            cx.subscribe_in(
                &spectra,
                window,
                |this, table, event: &TableEvent, _, cx| {
                    if let TableEvent::SelectRow(row) = event {
                        let id = table
                            .read(cx)
                            .delegate()
                            .rows
                            .get(*row)
                            .map(|r| r.spectrum_id);
                        this.select_spectrum(id, cx);
                    }
                },
            ),
        ];
        let mut browser = Self {
            path,
            mode,
            proteins,
            samples,
            peptides,
            spectra,
            drill: Drill::default(),
            tolerance: TOLERANCES[0],
            failure: None,
            generations: [0; 4],
            pending: [false; 4],
            wanted: [None; 3],
            all_rows: Vec::new(),
            filter: ProteinFilter::default(),
            search,
            min_spectra,
            organisms,
            modifications,
            go_terms,
            facets_changed: false,
            has_organisms: false,
            has_modifications: false,
            has_go: false,
            _subscriptions: subscriptions,
        };
        match data {
            Some(data) => browser.apply(data, cx),
            None => browser.reload(cx),
        }
        browser
    }

    /// Select `group_id` in the protein list (when it is listed) and show it.
    pub fn open_group(&mut self, group_id: i64, cx: &mut Context<Self>) {
        let row = self
            .proteins
            .read(cx)
            .delegate()
            .rows
            .iter()
            .position(|r| r.group_id == group_id);
        if let Some(row) = row {
            self.proteins
                .update(cx, |table, cx| table.set_selected_row(row, cx));
        }
        self.select_group(Some(group_id), cx);
    }

    /// Reload everything (after validation or a threshold change), keeping
    /// the selected group when it still passes.
    pub fn reload(&mut self, cx: &mut Context<Self>) {
        let path = self.path.clone();
        let prefer = self.wanted[0];
        self.start(
            Level::List,
            cx,
            move || browse::load(&path, prefer),
            |this, data, cx| this.apply(data, cx),
        );
    }

    fn loading(&self) -> bool {
        self.pending.iter().any(|&pending| pending)
    }

    /// Run `load` in the background, then `done` unless a load at this level
    /// or above started meanwhile.
    fn start<T: Send + 'static>(
        &mut self,
        level: Level,
        cx: &mut Context<Self>,
        load: impl FnOnce() -> Result<T, String> + Send + 'static,
        done: impl FnOnce(&mut Self, T, &mut Context<Self>) + 'static,
    ) {
        let index = level as usize;
        for deeper in index..self.generations.len() {
            self.generations[deeper] += 1;
            self.pending[deeper] = false;
        }
        self.pending[index] = true;
        let generation = self.generations[index];
        cx.notify();
        cx.spawn(async move |this, cx| {
            let result = cx.background_executor().spawn(async move { load() }).await;
            let _ = this.update(cx, |this, cx| {
                if this.generations[index] != generation {
                    return;
                }
                this.pending[index] = false;
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

    fn apply(&mut self, data: BrowserData, cx: &mut Context<Self>) {
        let selected = data.drill.group.as_ref().map(|g| g.group_id);
        // NSAF's denominators come from every group, not the filtered ones.
        let saf_totals = saf_totals(&data.list);
        self.all_rows = data.list.rows.clone();
        self.facets_changed = true;
        let shown = self.filter.apply(&self.all_rows);
        let measure = self.samples.read(cx).delegate().measure;
        let highest = highest_value(&shown, data.list.biosamples.len(), measure, &saf_totals);
        // Every table shows the score the thresholds use (q or probability).
        let mode = data.list.thresholds.mode;
        self.peptides
            .update(cx, |table, _| table.delegate_mut().mode = mode);
        self.spectra
            .update(cx, |table, _| table.delegate_mut().mode = mode);
        for table in [&self.proteins, &self.samples] {
            let (rows, biosamples) = (shown.clone(), data.list.biosamples.clone());
            table.update(cx, |table, cx| {
                // Mode and biosamples first: the sort may use either.
                table.delegate_mut().mode = mode;
                table.delegate_mut().biosamples = biosamples;
                table.delegate_mut().highest = highest;
                table.delegate_mut().saf_totals = saf_totals.clone();
                table.delegate_mut().set_rows(rows);
                table.refresh(cx);
                let row =
                    selected.and_then(|id| table.delegate().rows.position(|r| r.group_id == id));
                match row {
                    Some(row) => table.set_selected_row(row, cx),
                    None => table.clear_selection(cx),
                }
            });
        }
        self.show_drill(data.drill, cx);
    }

    /// Show a loaded drill-down, selecting its first peptide and spectrum.
    fn show_drill(&mut self, drill: Drill, cx: &mut Context<Self>) {
        let peptides = drill
            .group
            .as_ref()
            .map(|g| g.peptides.clone())
            .unwrap_or_default();
        let spectra = drill.spectra.rows.clone();
        self.wanted = [
            drill.group.as_ref().map(|g| g.group_id),
            drill.peptide_id,
            drill.spectrum.as_ref().map(|s| s.spectrum_id),
        ];
        self.drill = drill;
        let peptide = self.wanted[1];
        self.peptides.update(cx, |table, cx| {
            table.delegate_mut().set_rows(peptides);
            table.refresh(cx);
            // Select the loaded peptide wherever the current sort shows it.
            let row = peptide.and_then(|id| table.delegate().rows.position(|r| r.peptide_id == id));
            match row {
                Some(row) => table.set_selected_row(row, cx),
                None => table.clear_selection(cx),
            }
        });
        self.set_spectra(spectra, cx);
    }

    fn set_spectra(&mut self, rows: Vec<tuffold_store::SpectrumListRow>, cx: &mut Context<Self>) {
        let spectrum = self.wanted[2];
        self.spectra.update(cx, |table, cx| {
            table.delegate_mut().set_rows(rows);
            table.refresh(cx);
            let row =
                spectrum.and_then(|id| table.delegate().rows.position(|r| r.spectrum_id == id));
            match row {
                Some(row) => table.set_selected_row(row, cx),
                None => table.clear_selection(cx),
            }
        });
    }

    fn select_group(&mut self, group_id: Option<i64>, cx: &mut Context<Self>) {
        let Some(group_id) = group_id else { return };
        if self.wanted[0] == Some(group_id) {
            return;
        }
        self.wanted = [Some(group_id), None, None];
        let path = self.path.clone();
        self.start(
            Level::Group,
            cx,
            move || browse::load_group(&path, group_id),
            |this, drill, cx| this.show_drill(drill, cx),
        );
    }

    fn select_peptide(&mut self, peptide_id: Option<i64>, cx: &mut Context<Self>) {
        let Some(peptide_id) = peptide_id else { return };
        if self.wanted[1] == Some(peptide_id) {
            return;
        }
        self.wanted[1] = Some(peptide_id);
        self.wanted[2] = None;
        let path = self.path.clone();
        self.start(
            Level::Peptide,
            cx,
            move || browse::load_peptide(&path, peptide_id),
            move |this, (spectra, spectrum), cx| {
                this.wanted[2] = spectrum.as_ref().map(|s| s.spectrum_id);
                this.drill.peptide_id = Some(peptide_id);
                let rows = spectra.rows.clone();
                this.drill.spectra = spectra;
                this.drill.spectrum = spectrum;
                this.set_spectra(rows, cx);
            },
        );
    }

    fn select_spectrum(&mut self, spectrum_id: Option<i64>, cx: &mut Context<Self>) {
        let Some(spectrum_id) = spectrum_id else {
            return;
        };
        if self.wanted[2] == Some(spectrum_id) {
            return;
        }
        self.wanted[2] = Some(spectrum_id);
        let path = self.path.clone();
        self.start(
            Level::Spectrum,
            cx,
            move || browse::load_spectrum(&path, spectrum_id),
            |this, spectrum, _| this.drill.spectrum = spectrum,
        );
    }

    fn proteins_view(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let count = self.proteins.read(cx).delegate().rows.len();
        div()
            .flex()
            .size_full()
            .gap_4()
            .p_4()
            .child(
                div()
                    .flex()
                    .flex_col()
                    .gap_2()
                    .w(relative(0.42))
                    .h_full()
                    .child(list_heading("Protein groups", self.listed(count), palette))
                    .child(self.filter_bar(palette, cx))
                    .child(
                        div()
                            .flex_1()
                            .min_h_0()
                            .child(DataTable::new(&self.proteins).small()),
                    ),
            )
            .child(
                div()
                    .id("drill")
                    .flex_1()
                    .min_w_0()
                    .h_full()
                    .overflow_y_scroll()
                    .flex()
                    .flex_col()
                    .gap_4()
                    .pr_2()
                    .children(self.failure.as_deref().map(|message| {
                        callout(palette.fault, palette.fault_soft)
                            .text_size(px(13.))
                            .text_color(palette.fault)
                            .child(message.to_string())
                    }))
                    .children(self.drill.group.as_ref().map(|group| {
                        let starred = self
                            .all_rows
                            .iter()
                            .any(|row| row.group_id == group.group_id && row.starred);
                        let star = Button::new("star-group").xsmall().label(if starred {
                            "★ Starred"
                        } else {
                            "☆ Star"
                        });
                        let star = star.ghost().selected(starred);
                        group_header(group, palette).child(
                            star.on_click(cx.listener(|this, _, _, cx| this.toggle_star(cx))),
                        )
                    }))
                    .children(self.drill.group.as_ref().and_then(|group| {
                        let row = self
                            .all_rows
                            .iter()
                            .find(|r| r.group_id == group.group_id)?;
                        (!row.go.is_empty()).then(|| {
                            div()
                                .flex_none()
                                .text_size(px(12.))
                                .text_color(palette.ink_muted)
                                .child(format!("GO slim: {}", row.go.join(" · ")))
                        })
                    }))
                    .children(
                        self.drill
                            .group
                            .as_ref()
                            .map(|group| similarity_view(group, palette, cx)),
                    )
                    .child(
                        div()
                            .flex()
                            .gap_4()
                            .h(px(220.))
                            .flex_none()
                            .child(
                                div()
                                    .flex_1()
                                    .min_w_0()
                                    .flex()
                                    .flex_col()
                                    .gap_2()
                                    .child(list_heading(
                                        "Peptides",
                                        group_digits(self.peptides.read(cx).delegate().rows.len()),
                                        palette,
                                    ))
                                    .child(
                                        div()
                                            .flex_1()
                                            .min_h_0()
                                            .child(DataTable::new(&self.peptides).small()),
                                    ),
                            )
                            .child(
                                div()
                                    .flex_1()
                                    .min_w_0()
                                    .flex()
                                    .flex_col()
                                    .gap_2()
                                    .child(list_heading(
                                        "Spectra",
                                        page_label(&self.drill.spectra),
                                        palette,
                                    ))
                                    .child(
                                        div()
                                            .flex_1()
                                            .min_h_0()
                                            .child(DataTable::new(&self.spectra).small()),
                                    ),
                            ),
                    )
                    .children(
                        self.drill
                            .spectrum
                            .as_ref()
                            .map(|spectrum| self.spectrum_view(spectrum, palette, cx)),
                    )
                    .children(
                        self.drill
                            .group
                            .as_ref()
                            .map(|group| coverage_panel(group, palette)),
                    ),
            )
    }

    fn spectrum_view(
        &self,
        detail: &SpectrumDetail,
        palette: &Palette,
        cx: &mut Context<Self>,
    ) -> Div {
        let mode = self.spectra.read(cx).delegate().mode;
        let heading = detail.assigned.as_ref().map_or_else(
            || "No peptide assigned".to_string(),
            |assigned| {
                let decoy = if assigned.is_decoy { " · decoy" } else { "" };
                let score = match mode {
                    ThresholdMode::Fdr => format!("q {}", format_rate(assigned.q_value)),
                    ThresholdMode::Probability => {
                        format!("prob. {}", probability_text(Some(assigned.probability)))
                    }
                };
                let peptide = assigned
                    .localized_peptide
                    .clone()
                    .unwrap_or_else(|| proforma(&assigned.form));
                format!("{peptide} · {score}{decoy}")
            },
        );
        let site_line = detail.assigned.as_ref().and_then(|assigned| {
            let probs = assigned.site_probabilities.as_deref()?;
            let ascore = assigned
                .best_ascore
                .map_or_else(String::new, |a| format!(" · Ascore {a:.1}"));
            Some(format!("Sites: {probs}{ascore}"))
        });
        let facts = format!(
            "{} · {}{}{}",
            detail.run,
            detail.key,
            detail
                .precursor_charge
                .map_or_else(String::new, |z| format!(" · {z}+")),
            detail
                .precursor_mz
                .map_or_else(String::new, |mz| format!(" · m/z {mz:.4}")),
        );
        let chips = TOLERANCES.iter().enumerate().map(|(index, &tolerance)| {
            let selected = tolerance == self.tolerance;
            let button = Button::new(("tolerance", index))
                .xsmall()
                .label(tolerance_label(tolerance));
            let button = button.ghost().selected(selected);
            button.on_click(cx.listener(move |this, _, _, cx| {
                this.tolerance = tolerance;
                cx.notify();
            }))
        });
        let body = match spectrum_plot(detail, self.tolerance) {
            Ok(plot) => div()
                .flex()
                .flex_col()
                .gap_3()
                .child(
                    div()
                        .overflow_hidden()
                        .child(ladder_view(&plot.ladder, palette)),
                )
                .child(
                    div()
                        .font_family(MONO)
                        .text_size(px(12.))
                        .text_color(palette.ink_muted)
                        .child(format!(
                            "{} of {} b/y ions · {:.0}% of intensity explained",
                            plot.matched,
                            plot.theoretical,
                            plot.explained * 100.0
                        )),
                )
                .child(spectrum_plot_view(&plot, palette)),
            Err(why) => div()
                .text_size(px(13.))
                .text_color(palette.ink_muted)
                .child(match why {
                    NoPlot::NoPeaks => {
                        "No peaks for this spectrum: add the run's .mzML or .mgf with Add data… \
                         to see it annotated."
                    }
                    NoPlot::Unreadable => "This spectrum's stored peaks could not be read.",
                }),
        };
        panel_shell(palette)
            .flex_none()
            .child(
                div()
                    .flex()
                    .items_start()
                    .justify_between()
                    .gap_4()
                    .child(
                        div()
                            .flex()
                            .flex_col()
                            .min_w_0()
                            .child(eyebrow("Spectrum", palette))
                            .child(
                                div()
                                    .font_family(MONO)
                                    .text_size(px(14.))
                                    .font_weight(FontWeight::MEDIUM)
                                    .truncate()
                                    .child(heading),
                            )
                            .child(
                                div()
                                    .font_family(MONO)
                                    .text_size(px(12.))
                                    .text_color(palette.ink_faint)
                                    .child(facts),
                            )
                            .children(site_line.map(|line| {
                                div()
                                    .font_family(MONO)
                                    .text_size(px(12.))
                                    .text_color(palette.signal)
                                    .truncate()
                                    .child(line)
                            })),
                    )
                    .child(track(palette).flex_none().children(chips)),
            )
            .child(body)
            .children((!detail.psms.is_empty()).then(|| engines_view(detail, palette)))
    }

    /// List the groups the filter keeps, keeping the selection when it stays.
    fn refilter(&mut self, cx: &mut Context<Self>) {
        let shown = self.filter.apply(&self.all_rows);
        for table in [&self.proteins, &self.samples] {
            let rows = shown.clone();
            table.update(cx, |table, cx| {
                let selected = table
                    .selected_row()
                    .and_then(|row| table.delegate().rows.get(row))
                    .map(|row| row.group_id);
                let delegate = table.delegate_mut();
                delegate.highest = highest_value(
                    &rows,
                    delegate.biosamples.len(),
                    delegate.measure,
                    &delegate.saf_totals,
                );
                delegate.set_rows(rows);
                table.refresh(cx);
                match selected.and_then(|id| table.delegate().rows.position(|r| r.group_id == id)) {
                    Some(row) => table.set_selected_row(row, cx),
                    None => table.clear_selection(cx),
                }
            });
        }
        // The group on show was filtered out: show the first one listed
        // instead. (When nothing is listed, the detail stays as it was.)
        let shown_group = self.wanted[0];
        if !shown.iter().any(|r| Some(r.group_id) == shown_group)
            && let Some(first) = shown.first().map(|r| r.group_id)
        {
            self.open_group(first, cx);
        }
        cx.notify();
    }

    /// Keep only the groups whose leader reaches GO slim term `term`
    /// (`name (aspect)`, as the menu lists it); `None` clears it.
    pub fn set_go_filter(&mut self, term: Option<String>, cx: &mut Context<Self>) {
        self.filter.go = term;
        self.facets_changed = true;
        self.refilter(cx);
    }

    fn toggle_starred_only(&mut self, cx: &mut Context<Self>) {
        self.filter.starred_only = !self.filter.starred_only;
        self.refilter(cx);
    }

    /// Star or unstar the shown group's leader, in the project and here.
    fn toggle_star(&mut self, cx: &mut Context<Self>) {
        let Some(group) = self.drill.group.as_ref() else {
            return;
        };
        let Some(row) = self.all_rows.iter().find(|r| r.group_id == group.group_id) else {
            return;
        };
        let (accession, starred) = (row.accession.clone(), !row.starred);
        let path = self.path.clone();
        cx.spawn(async move |this, cx| {
            let saved = cx
                .background_executor()
                .spawn(async move { crate::project::set_starred(&path, &accession, starred) })
                .await;
            let _ = this.update(cx, |this, cx| {
                match saved {
                    Ok(()) => {
                        let group_id = this.drill.group.as_ref().map(|g| g.group_id);
                        for row in this
                            .all_rows
                            .iter_mut()
                            .filter(|r| Some(r.group_id) == group_id)
                        {
                            row.starred = starred;
                        }
                        this.failure = None;
                        this.refilter(cx);
                    }
                    Err(message) => this.failure = Some(message),
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn clear_filter(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        self.filter = ProteinFilter::default();
        for input in [&self.search, &self.min_spectra] {
            input.update(cx, |input, cx| input.set_value("", window, cx));
        }
        for select in [&self.organisms, &self.modifications, &self.go_terms] {
            select.update(cx, |select, cx| select.set_selected_index(None, window, cx));
        }
        self.refilter(cx);
    }

    /// Search, minimum spectra, organism, modification, and starred.
    fn filter_bar(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let starred = Button::new("filter-starred").xsmall().label("★ Starred");
        let starred = starred.ghost().selected(self.filter.starred_only);
        div()
            .flex()
            .flex_wrap()
            .items_center()
            .gap_2()
            .child(div().w(px(200.)).child(Input::new(&self.search).small()))
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap_1()
                    .text_size(px(12.))
                    .text_color(palette.ink_faint)
                    .child("Spectra ≥")
                    .child(
                        div()
                            .w(px(48.))
                            .child(Input::new(&self.min_spectra).small()),
                    ),
            )
            .children(self.has_organisms.then(|| {
                div().w(px(180.)).child(
                    Select::new(&self.organisms)
                        .small()
                        .placeholder("Any organism")
                        .cleanable(true)
                        .menu_width(px(320.)),
                )
            }))
            .children(self.has_modifications.then(|| {
                div().w(px(160.)).child(
                    Select::new(&self.modifications)
                        .small()
                        .placeholder("Any modification")
                        .cleanable(true)
                        .menu_width(px(240.)),
                )
            }))
            .children(self.has_go.then(|| {
                div().w(px(200.)).child(
                    Select::new(&self.go_terms)
                        .small()
                        .placeholder("Any GO term")
                        .cleanable(true)
                        .menu_width(px(340.)),
                )
            }))
            .child(starred.on_click(cx.listener(|this, _, _, cx| this.toggle_starred_only(cx))))
            .children((!self.filter.is_empty()).then(|| {
                Button::new("filter-clear")
                    .xsmall()
                    .ghost()
                    .label("Clear")
                    .on_click(cx.listener(|this, _, window, cx| this.clear_filter(window, cx)))
            }))
    }

    /// "133", or "12 of 133" when the filter hides some.
    fn listed(&self, shown: usize) -> String {
        if shown == self.all_rows.len() {
            group_digits(shown)
        } else {
            format!(
                "{} of {}",
                group_digits(shown),
                group_digits(self.all_rows.len())
            )
        }
    }

    /// Show another measure in the Samples view, rescaling its heat map.
    pub fn set_measure(&mut self, measure: SampleMeasure, cx: &mut Context<Self>) {
        self.samples.update(cx, |table, cx| {
            let delegate = table.delegate_mut();
            delegate.measure = measure;
            delegate.highest = highest_value(
                delegate.rows.iter(),
                delegate.biosamples.len(),
                measure,
                &delegate.saf_totals,
            );
            // Re-sort by the same column with the new values.
            let rows: Vec<_> = delegate.rows.iter().cloned().collect();
            delegate.set_rows(rows);
            table.refresh(cx);
        });
        cx.notify();
    }

    fn samples_view(&self, palette: &Palette, cx: &mut Context<Self>) -> Div {
        let table = self.samples.read(cx).delegate();
        let (count, measure) = (table.rows.len(), table.measure);
        let buttons = SampleMeasure::ALL
            .into_iter()
            .enumerate()
            .map(|(index, option)| {
                let button = Button::new(("sample-measure", index))
                    .xsmall()
                    .label(option.label());
                let button = button.ghost().selected(option == measure);
                button.on_click(cx.listener(move |this, _, _, cx| this.set_measure(option, cx)))
            });
        div()
            .flex()
            .flex_col()
            .size_full()
            .gap_2()
            .p_4()
            .child(
                div()
                    .flex()
                    .items_end()
                    .justify_between()
                    .child(list_heading(
                        &format!("{} by biosample", measure.label()),
                        self.listed(count),
                        palette,
                    ))
                    .child(track(palette).children(buttons)),
            )
            .child(self.filter_bar(palette, cx))
            // The default measure needs no explaining.
            .children((measure != SampleMeasure::Spectra).then(|| {
                div()
                    .text_size(px(12.))
                    .text_color(palette.ink_faint)
                    .child(measure.explain())
            }))
            .child(
                div()
                    .flex_1()
                    .min_h_0()
                    .child(DataTable::new(&self.samples).small()),
            )
    }
}

/// Similar groups shown as matrix columns, most shared first.
const SIMILAR_COLUMNS: usize = 8;

/// Which other passing groups share this group's peptides (Scaffold's
/// similarity view): the shared peptides down, this group and the others
/// across, a dot where a group holds the peptide. Clicking a column opens
/// that group.
fn similarity_view(group: &GroupDetail, palette: &Palette, cx: &mut Context<Browser>) -> Div {
    let header = div()
        .flex()
        .items_center()
        .gap_2()
        .child(panel_title("Similarity"))
        .child(
            div()
                .text_size(px(12.))
                .text_color(palette.ink_faint)
                .child(if group.similar.is_empty() {
                    "No other passing group shares its peptides.".to_string()
                } else {
                    format!(
                        "Shares peptides with {} group{} · cluster of {} groups",
                        group.similar.len(),
                        if group.similar.len() == 1 { "" } else { "s" },
                        group.cluster_size
                    )
                }),
        );
    if group.similar.is_empty() {
        return div().flex_none().child(header);
    }
    let columns = &group.similar[..group.similar.len().min(SIMILAR_COLUMNS)];
    let shared: Vec<&tuffold_store::GroupPeptide> = group
        .peptides
        .iter()
        .filter(|p| columns.iter().any(|c| c.shared.contains(&p.peptide_id)))
        .collect();
    let cell = |width: f32| div().w(px(width)).flex_none().flex().justify_center();
    let leader = group
        .members
        .first()
        .map_or_else(String::new, |m| m.accession.clone());
    let heading = div()
        .flex()
        .items_end()
        .gap_1()
        .text_size(px(11.))
        .font_family(MONO)
        .child(
            div()
                .w(px(220.))
                .flex_none()
                .text_color(palette.ink_faint)
                .child("peptide"),
        )
        .child(
            cell(84.).child(
                div()
                    .truncate()
                    .text_color(palette.ink)
                    .child(short_accession(&leader)),
            ),
        )
        .children(columns.iter().map(|column| {
            let id = column.group_id;
            cell(84.).child(
                div()
                    .id(("similar", usize::try_from(id).unwrap_or_default()))
                    .truncate()
                    .text_color(palette.signal)
                    .cursor_pointer()
                    .hover(|style| style.underline())
                    .child(short_accession(&column.accession))
                    .on_click(cx.listener(move |this, _, _, cx| this.open_group(id, cx))),
            )
        }));
    let rows = shared.iter().map(|peptide| {
        let dot = |held: bool| {
            cell(84.).child(div().size(px(8.)).rounded_full().bg(if held {
                palette.signal
            } else {
                palette.rule
            }))
        };
        div()
            .flex()
            .items_center()
            .gap_1()
            .h(px(20.))
            .text_size(px(11.))
            .font_family(MONO)
            .child(
                div()
                    .w(px(220.))
                    .flex_none()
                    .truncate()
                    .text_color(palette.ink)
                    .child(peptide.proforma.clone()),
            )
            .child(dot(true))
            .children(
                columns
                    .iter()
                    .map(|column| dot(column.shared.contains(&peptide.peptide_id))),
            )
    });
    let unique = group.peptides.len() - shared.len();
    div()
        .flex()
        .flex_col()
        .flex_none()
        .gap_1()
        .child(header)
        .child(heading)
        .children(rows)
        .child(
            div()
                .text_size(px(11.))
                .text_color(palette.ink_faint)
                .child(format!(
                    "{unique} peptide{} only this group holds not shown{}",
                    if unique == 1 { "" } else { "s" },
                    if group.similar.len() > SIMILAR_COLUMNS {
                        format!(
                            " · {} more groups not shown",
                            group.similar.len() - SIMILAR_COLUMNS
                        )
                    } else {
                        String::new()
                    }
                )),
        )
}

/// `P02768` from `sp|P02768|ALBU_HUMAN`; other accessions as they are.
fn short_accession(accession: &str) -> String {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        [_, id, _] => (*id).to_string(),
        _ => accession.to_string(),
    }
}

/// "12", or "500 of 1,234" when only the best are listed.
fn page_label(page: &tuffold_store::SpectrumPage) -> String {
    if page.is_capped() {
        format!(
            "best {} of {}",
            group_digits(page.rows.len()),
            group_digits(page.total)
        )
    } else {
        group_digits(page.total)
    }
}

fn list_heading(title: &str, count: impl Into<SharedString>, palette: &Palette) -> Div {
    div()
        .flex()
        .items_center()
        .gap_2()
        .child(panel_title(title.to_string()))
        .child(
            div()
                .px_2()
                .rounded(px(3.))
                .bg(palette.paper_sunken)
                .font_family(MONO)
                .text_size(px(11.))
                .text_color(palette.ink_muted)
                .child(count.into()),
        )
}

fn proforma(form: &tuffold_core::PeptideForm) -> String {
    form.proforma()
}

/// The group's leader, description, members, and headline numbers.
fn group_header(group: &GroupDetail, palette: &Palette) -> Div {
    let leader = group.members.first();
    let others: Vec<String> = group
        .members
        .iter()
        .skip(1)
        .map(|m| format!("{} ({})", m.accession, m.role))
        .collect();
    // Unique and razor peptides only, like the protein list: shared spectra
    // are credited to another group.
    let spectra: usize = group
        .peptides
        .iter()
        .filter(|p| p.role != "shared")
        .map(|p| p.spectra)
        .sum();
    let coverage = group.sequence.as_deref().map(|sequence| {
        tuffold_core::coverage::covered_fraction(sequence.chars().count(), &group.coverage)
    });
    // Counted like the protein list's "Pep." column: distinct unique and razor
    // sequences (the table below lists every modified form, shared ones too).
    let sequences = group
        .peptides
        .iter()
        .filter(|p| p.role != "shared")
        .map(|p| p.sequence.as_str())
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    let facts = format!(
        "{} peptide sequences · {} spectra{}",
        group_digits(sequences),
        group_digits(spectra),
        coverage.map_or_else(String::new, |c| format!(" · {:.0}% coverage", c * 100.0)),
    );
    div()
        .flex()
        .flex_none()
        .items_end()
        .justify_between()
        .gap_4()
        .pb_2()
        .border_b_1()
        .border_color(palette.rule)
        .child(
            div()
                .flex()
                .flex_col()
                .min_w_0()
                .child(eyebrow("Protein group", palette))
                .child(
                    div()
                        .font_family(MONO)
                        .text_size(px(18.))
                        .font_weight(FontWeight::SEMIBOLD)
                        .truncate()
                        .child(leader.map_or_else(String::new, |m| m.accession.clone())),
                )
                .children(
                    leader
                        .and_then(|m| m.description.clone())
                        .map(|description| {
                            div()
                                .text_size(px(13.))
                                .text_color(palette.ink_muted)
                                .truncate()
                                .child(description)
                        }),
                )
                .children((!others.is_empty()).then(|| {
                    div()
                        .text_size(px(12.))
                        .text_color(palette.ink_faint)
                        .truncate()
                        .child(format!("Also: {}", others.join(", ")))
                })),
        )
        .child(
            div()
                .flex_none()
                .font_family(MONO)
                .text_size(px(12.))
                .text_color(palette.ink_muted)
                .child(facts),
        )
}

/// The leader's sequence with its peptides highlighted.
fn coverage_panel(group: &GroupDetail, palette: &Palette) -> Div {
    let body = match group.sequence.as_deref() {
        Some(sequence) => coverage_view(
            &coverage_lines(sequence, &group.coverage, COVERAGE_LINE),
            palette,
        ),
        None => div()
            .text_size(px(13.))
            .text_color(palette.ink_muted)
            .child(
                "No sequence for this protein, so coverage cannot be shown. Import results that \
                 carry sequences (a FragPipe folder with its FASTA, or mzIdentML) to see it.",
            ),
    };
    panel_shell(palette)
        .flex_none()
        .child(eyebrow("Sequence coverage", palette))
        .child(body)
}

fn engines_view(detail: &SpectrumDetail, palette: &Palette) -> Div {
    let rows = detail.psms.iter().map(|psm| {
        let scores = psm
            .scores
            .iter()
            .map(|(name, value)| format!("{name} {}", short_number(*value)))
            .collect::<Vec<_>>()
            .join(" · ");
        div()
            .flex()
            .gap_3()
            .py_1()
            .border_t_1()
            .border_color(palette.rule)
            .text_size(px(12.))
            .child(
                div()
                    .w(px(150.))
                    .flex_none()
                    .truncate()
                    .child(psm.engine.clone()),
            )
            .child(
                div()
                    .w(px(180.))
                    .flex_none()
                    .font_family(MONO)
                    .truncate()
                    .child(psm.proforma.clone()),
            )
            .child(
                div()
                    .w(px(64.))
                    .flex_none()
                    .font_family(MONO)
                    .child(psm.q_value.map_or_else(|| "—".into(), format_rate)),
            )
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .truncate()
                    .text_color(palette.ink_muted)
                    .child(scores),
            )
    });
    div()
        .flex()
        .flex_col()
        .child(eyebrow("Engines", palette))
        .children(rows)
}

/// Compact numbers for score lists: `0.00012` -> `1.2e-4`, `23.456` -> `23.46`.
pub fn short_number(value: f64) -> String {
    if value != 0.0 && (value.abs() < 0.01 || value.abs() >= 1e5) {
        format!("{value:.1e}")
    } else {
        let text = format!("{value:.2}");
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

impl Render for Browser {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let palette = *cx.global::<Palette>();
        if self.facets_changed {
            self.facets_changed = false;
            let found = facets(&self.all_rows);
            self.has_organisms = !found.organisms.is_empty();
            self.has_modifications = !found.modifications.is_empty();
            self.has_go = !found.go.is_empty();
            self.organisms.update(cx, |select, cx| {
                select.set_items(choices(found.organisms), window, cx)
            });
            self.modifications.update(cx, |select, cx| {
                select.set_items(choices(found.modifications), window, cx)
            });
            self.go_terms.update(cx, |select, cx| {
                select.set_items(choices(found.go), window, cx)
            });
            // New items, so point each menu at the filter's value again (its
            // count, and its place in the list, may have changed).
            for (select, value) in [
                (&self.organisms, &self.filter.organism),
                (&self.modifications, &self.filter.modification),
                (&self.go_terms, &self.filter.go),
            ] {
                if let Some(value) = value {
                    select.update(cx, |select, cx| {
                        select.set_selected_value(value, window, cx)
                    });
                }
            }
        }
        let empty = self.all_rows.is_empty();
        if empty && !self.loading() {
            return div()
                .p(px(48.))
                .text_size(px(15.))
                .text_color(palette.ink_muted)
                .child(
                    "No protein group passes the current thresholds. Loosen them on the Overview, \
                     or check the validation notes there.",
                )
                .into_any_element();
        }
        match self.mode {
            BrowserMode::Proteins => self.proteins_view(&palette, cx).into_any_element(),
            BrowserMode::Samples => self.samples_view(&palette, cx).into_any_element(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::short_number;

    #[test]
    fn numbers_stay_short() {
        assert_eq!(short_number(23.456), "23.46");
        assert_eq!(short_number(12.0), "12");
        assert_eq!(short_number(0.00012), "1.2e-4");
        assert_eq!(short_number(0.0), "0");
    }
}
