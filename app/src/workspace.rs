//! The root view: creates, opens, and extends `.tuff` projects, running every
//! import off the UI thread, and shows the empty state or the project — its
//! Overview, and once validated, the linked Proteins and Samples views.

use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::Arc;

use gpui_kit::component::button::{Button, ButtonVariants as _};
use gpui_kit::component::menu::DropdownMenu as _;
use gpui_kit::component::{Disableable as _, IconName, Sizable as _};
use gpui_kit::{
    Anchor, App, AppContext as _, AsyncApp, Context, Div, Entity, FocusHandle, FontWeight,
    InteractiveElement as _, IntoElement, ParentElement, PathPromptOptions, Render,
    StatefulInteractiveElement as _, Styled, Subscription, WeakEntity, Window, actions, div, px,
};

use tuffold_validate::Thresholds;

use crate::browse::{BrowserData, SampleMeasure};
use crate::project::{self, ProjectReport, ValidationState};
use crate::theme::{MONO, Palette, SANS};
use crate::views::{
    Browser, BrowserMode, EmptyStateProps, ModelsOfferProps, OpenGroup, QuantView, ShowProtein,
    ThresholdHandler, ValidateHandler, ValidationControls, empty_state, failure_banner,
    models_offer, notice_banner, project_view,
};

actions!(
    tuffold,
    [
        NewProject,
        OpenProject,
        OpenDemo,
        AddData,
        SearchSpectra,
        ExportWorkbook,
        ExportLibrary,
        ExportMzid,
        AddGoAnnotations,
        SaveDiagnostics,
        UndoLastChange,
        ShowOverview,
        ShowProteins,
        ShowSamples,
        ShowQuant,
        Quit
    ]
);

const SUGGESTED_NAME: &str = "study.tuff";
const DEMO_NAME: &str = "tuffold-demo.tuff";

/// The views of an open project.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Overview,
    Proteins,
    Samples,
    Quant,
}

pub struct Workspace {
    project: Option<Arc<ProjectReport>>,
    tab: Tab,
    /// The linked views, once the project is validated.
    browser: Option<Entity<Browser>>,
    _browser_events: Option<Subscription>,
    /// The Quant tab, created on first visit.
    quant: Option<Entity<QuantView>>,
    _quant_events: Option<Subscription>,
    /// Redraws the toolbar when the Quant view starts or ends a write.
    _quant_writes: Option<Subscription>,
    /// What is running; also the guard against starting a second job.
    busy: Option<&'static str>,
    /// The last failure (a failed open, create, or add).
    failure: Option<String>,
    /// The last export written, for the banner.
    exported: Option<String>,
    /// Offer the deep-learning models: none are installed and this wasn't
    /// dismissed.
    offer_models: bool,
    focus_handle: FocusHandle,
}

impl Workspace {
    pub fn new(window: &mut Window, cx: &mut Context<Self>) -> Self {
        let focus_handle = cx.focus_handle();
        window.focus(&focus_handle, cx);
        Self {
            project: None,
            tab: Tab::Overview,
            browser: None,
            _browser_events: None,
            quant: None,
            _quant_events: None,
            _quant_writes: None,
            busy: None,
            failure: None,
            exported: None,
            offer_models: false,
            focus_handle,
        }
    }

    /// Show `measure` in the Samples view (snapshots and tests).
    pub fn set_sample_measure(&mut self, measure: SampleMeasure, cx: &mut Context<Self>) {
        if let Some(browser) = &self.browser {
            browser.update(cx, |browser, cx| browser.set_measure(measure, cx));
        }
    }

    pub fn set_go_filter(&mut self, term: Option<String>, cx: &mut Context<Self>) {
        if let Some(browser) = &self.browser {
            browser.update(cx, |browser, cx| browser.set_go_filter(term, cx));
        }
    }

    /// A workspace showing `report` on `tab` with browser data already loaded
    /// (snapshots and tests).
    pub fn with_browser(
        report: ProjectReport,
        data: BrowserData,
        tab: Tab,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let mut workspace = Self::with_project(report, window, cx);
        workspace.attach_browser(Some(data), window, cx);
        workspace.set_tab(tab, cx);
        workspace
    }

    fn validated(&self) -> bool {
        self.project
            .as_ref()
            .is_some_and(|report| matches!(report.validation, ValidationState::Done(_)))
    }

    /// Create the browser for the open, validated project.
    fn attach_browser(
        &mut self,
        data: Option<BrowserData>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        let mode = browser_mode(self.tab);
        let browser = cx.new(|cx| Browser::new(path, mode, data, window, cx));
        self._browser_events = Some(cx.subscribe(&browser, |this, _, _: &ShowProtein, cx| {
            this.set_tab(Tab::Proteins, cx);
        }));
        self.browser = Some(browser);
    }

    fn set_tab(&mut self, tab: Tab, cx: &mut Context<Self>) {
        self.tab = tab;
        self.exported = None;
        let mode = browser_mode(tab);
        if let Some(browser) = &self.browser {
            browser.update(cx, |browser, cx| {
                browser.mode = mode;
                cx.notify();
            });
        }
        cx.notify();
    }

    fn show_overview(&mut self, _: &ShowOverview, _: &mut Window, cx: &mut Context<Self>) {
        self.set_tab(Tab::Overview, cx);
    }

    fn show_proteins(&mut self, _: &ShowProteins, _: &mut Window, cx: &mut Context<Self>) {
        if self.validated() {
            self.set_tab(Tab::Proteins, cx);
        }
    }

    fn show_samples(&mut self, _: &ShowSamples, _: &mut Window, cx: &mut Context<Self>) {
        if self.validated() {
            self.set_tab(Tab::Samples, cx);
        }
    }

    fn show_quant(&mut self, _: &ShowQuant, _: &mut Window, cx: &mut Context<Self>) {
        if self.validated() {
            self.set_tab(Tab::Quant, cx);
        }
    }

    /// Create the Quant view for the open project (`data` shows immediately).
    fn attach_quant(
        &mut self,
        data: Option<crate::quant::QuantData>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        let quant = cx.new(|cx| QuantView::new(path, data, window, cx));
        self._quant_events = Some(cx.subscribe(&quant, |this, _, event: &OpenGroup, cx| {
            if let Some(browser) = &this.browser {
                browser.update(cx, |browser, cx| browser.open_group(event.0, cx));
            }
            this.set_tab(Tab::Proteins, cx);
        }));
        self._quant_writes = Some(cx.observe(&quant, |_, _, cx| cx.notify()));
        self.quant = Some(quant);
    }

    /// A workspace on the Quant tab with its data loaded, showing `choice`
    /// (snapshots).
    pub fn with_quant(
        report: ProjectReport,
        data: BrowserData,
        quant_data: crate::quant::QuantData,
        choice: crate::quant::QuantChoice,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let mut workspace = Self::with_browser(report, data, Tab::Quant, window, cx);
        workspace.attach_quant(Some(quant_data), window, cx);
        if let Some(quant) = &workspace.quant {
            quant.update(cx, |quant, cx| {
                quant.set_test(choice.test, cx);
                quant.set_method(choice.method, cx);
                quant.set_scope(choice.scope, cx);
            });
        }
        workspace
    }

    fn open_demo(&mut self, _: &OpenDemo, _: &mut Window, cx: &mut Context<Self>) {
        if !self.begin("Building the demo…", cx) {
            return;
        }
        let directory = std::env::var_os("HOME").map_or_else(|| PathBuf::from("."), PathBuf::from);
        let save = cx.prompt_for_new_path(&directory, Some(DEMO_NAME));
        cx.spawn(async move |this, cx| {
            let Some(path) = save.await.ok().and_then(Result::ok).flatten() else {
                return finish(&this, cx, None);
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::create_demo(&path) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Look, off the UI thread, for the deep-learning models, to offer them
    /// when there are none.
    pub fn check_models(&mut self, cx: &mut Context<Self>) {
        cx.spawn(async move |this, cx| {
            let status = cx
                .background_executor()
                .spawn(async { tuffold_engine::models_status() })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.offer_models = !status.installed;
                cx.notify();
            });
        })
        .detach();
    }

    /// Show the offer of models (snapshots).
    pub fn show_models_offer(&mut self) {
        self.offer_models = true;
    }

    /// Pick the downloaded archive and install the models from it.
    fn install_models(&mut self, cx: &mut Context<Self>) {
        if !self.begin("Installing models…", cx) {
            return;
        }
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Install".into()),
        });
        cx.spawn(async move |this, cx| {
            let archive = picked_paths(picked.await).and_then(|paths| paths.into_iter().next());
            let result = match archive {
                None => None,
                Some(archive) => Some(
                    cx.background_executor()
                        .spawn(async move { tuffold_engine::install_models(&archive) })
                        .await,
                ),
            };
            let _ = this.update(cx, |this, cx| {
                this.busy = None;
                match result {
                    Some(Ok(_)) => {
                        this.failure = None;
                        this.offer_models = false;
                        this.exported =
                            Some("Models installed. Validate again to use them.".into());
                    }
                    Some(Err(message)) => this.failure = Some(message),
                    None => {}
                }
                cx.notify();
            });
        })
        .detach();
    }

    /// Save the report to send with a bug. Works with no project open, since
    /// that may be what's broken.
    fn save_diagnostics(&mut self, _: &SaveDiagnostics, _: &mut Window, cx: &mut Context<Self>) {
        if !self.begin("Saving diagnostics…", cx) {
            return;
        }
        let source = self.project.as_ref().map(|report| report.path.clone());
        let directory = source
            .as_deref()
            .and_then(Path::parent)
            .map(Path::to_path_buf)
            .or_else(|| std::env::var_os("HOME").map(PathBuf::from))
            .unwrap_or_else(|| PathBuf::from("."));
        let save = cx.prompt_for_new_path(&directory, Some("tuffold-diagnostics.txt"));
        cx.spawn(async move |this, cx| {
            let target = save.await.ok().and_then(Result::ok).flatten();
            let result = match target {
                None => None,
                Some(target) => Some(
                    cx.background_executor()
                        .spawn(async move {
                            let report = tuffold_engine::diagnostics_report(source.as_deref());
                            std::fs::write(&target, report)
                                .map(|()| {
                                    format!(
                                        "Saved {}. Read it before sending it.",
                                        target.display()
                                    )
                                })
                                .map_err(|e| format!("{}: {e}", target.display()))
                        })
                        .await,
                ),
            };
            let _ = this.update(cx, |this, cx| {
                this.busy = None;
                match result {
                    Some(Ok(message)) => {
                        this.failure = None;
                        this.exported = Some(message);
                    }
                    Some(Err(message)) => this.failure = Some(message),
                    None => {}
                }
                cx.notify();
            });
        })
        .detach();
    }

    /// A workspace already showing `report` (snapshots and tests).
    pub fn with_project(
        report: ProjectReport,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        Self {
            project: Some(Arc::new(report)),
            ..Self::new(window, cx)
        }
    }

    /// A workspace showing a failure on the empty state (snapshots and tests).
    pub fn with_failure(message: String, window: &mut Window, cx: &mut Context<Self>) -> Self {
        Self {
            failure: Some(message),
            ..Self::new(window, cx)
        }
    }

    /// A job of ours is running, or the Quant view is writing the project
    /// (matching between runs, a plex or group save): nothing else may start.
    fn blocked(&self, cx: &App) -> bool {
        self.busy.is_some() || self.quant.as_ref().is_some_and(|q| q.read(cx).is_writing())
    }

    fn begin(&mut self, label: &'static str, cx: &mut Context<Self>) -> bool {
        if self.blocked(cx) {
            return false;
        }
        self.busy = Some(label);
        // The last export's message is old news once something else starts.
        self.exported = None;
        cx.notify();
        true
    }

    fn new_project(&mut self, _: &NewProject, _: &mut Window, cx: &mut Context<Self>) {
        if !self.begin("Importing…", cx) {
            return;
        }
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: true,
            multiple: true,
            prompt: Some("Import".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(inputs) = picked_paths(picked.await) else {
                return finish(&this, cx, None);
            };
            let directory = inputs
                .first()
                .and_then(|input| input.parent())
                .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
            let Ok(save) = this.update(cx, |_, cx| {
                cx.prompt_for_new_path(&directory, Some(SUGGESTED_NAME))
            }) else {
                return;
            };
            let Some(path) = save.await.ok().and_then(Result::ok).flatten() else {
                return finish(&this, cx, None);
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::create(&path, &inputs) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    fn open_project(&mut self, _: &OpenProject, _: &mut Window, cx: &mut Context<Self>) {
        if !self.begin("Opening…", cx) {
            return;
        }
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Open".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(path) = picked_paths(picked.await).and_then(|mut paths| paths.pop()) else {
                return finish(&this, cx, None);
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::open(&path) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Write the validated project's Excel workbook where the user picks.
    fn export_workbook(&mut self, _: &ExportWorkbook, _: &mut Window, cx: &mut Context<Self>) {
        self.export("xlsx", project::export_workbook, cx);
    }

    /// Write the validated project's spectral library (.msp) where the user
    /// picks.
    fn export_library(&mut self, _: &ExportLibrary, _: &mut Window, cx: &mut Context<Self>) {
        self.export("msp", project::export_library, cx);
    }

    /// Write the validated project's results as mzIdentML where the user
    /// picks.
    fn export_mzid(&mut self, _: &ExportMzid, _: &mut Window, cx: &mut Context<Self>) {
        self.export("mzid", project::export_mzid, cx);
    }

    /// Ask where to save an export with `extension`, then `write` it from the
    /// project on a background thread.
    fn export(
        &mut self,
        extension: &'static str,
        write: fn(&Path, &Path) -> Result<String, String>,
        cx: &mut Context<Self>,
    ) {
        let Some(source) = self.project.as_ref().map(|r| r.path.clone()) else {
            return;
        };
        if !self.validated() || !self.begin("Exporting…", cx) {
            return;
        }
        let directory = source
            .parent()
            .map_or_else(|| PathBuf::from("."), PathBuf::from);
        let name = source.file_stem().map_or_else(
            || format!("tuffold.{extension}"),
            |stem| format!("{}.{extension}", stem.to_string_lossy()),
        );
        let save = cx.prompt_for_new_path(&directory, Some(&name));
        cx.spawn(async move |this, cx| {
            let target = save.await.ok().and_then(Result::ok).flatten();
            let result = match target {
                None => None,
                Some(target) => Some(
                    cx.background_executor()
                        .spawn(async move { write(&source, &target) })
                        .await,
                ),
            };
            let _ = this.update(cx, |this, cx| {
                this.busy = None;
                match result {
                    Some(Ok(message)) => {
                        this.failure = None;
                        this.exported = Some(message);
                    }
                    Some(Err(message)) => this.failure = Some(message),
                    None => {}
                }
                cx.notify();
            });
        })
        .detach();
    }

    /// Pick the GO ontology, a slim, and annotation files, and add the
    /// project's GO slim terms.
    fn add_go(&mut self, _: &AddGoAnnotations, _: &mut Window, cx: &mut Context<Self>) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        if !self.begin("Adding GO annotations…", cx) {
            return;
        }
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: true,
            prompt: Some("Add GO".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(inputs) = picked_paths(picked.await) else {
                let _ = this.update(cx, |this, cx| {
                    this.busy = None;
                    cx.notify();
                });
                return;
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::add_go(&path, &inputs) })
                .await;
            let _ = this.update(cx, |this, cx| {
                this.busy = None;
                match result {
                    Ok(message) => {
                        this.failure = None;
                        this.exported = Some(message);
                        if let Some(browser) = &this.browser {
                            browser.update(cx, |browser, cx| browser.reload(cx));
                        }
                    }
                    Err(message) => this.failure = Some(message),
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn add_data(&mut self, _: &AddData, _: &mut Window, cx: &mut Context<Self>) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        if !self.begin("Importing…", cx) {
            return;
        }
        let picked = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: true,
            multiple: true,
            prompt: Some("Add".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(inputs) = picked_paths(picked.await) else {
                return finish(&this, cx, None);
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::add(&path, &inputs) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Pick a FASTA, then spectra files, search them with the built-in
    /// search, and validate.
    fn search_spectra(&mut self, _: &SearchSpectra, _: &mut Window, cx: &mut Context<Self>) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        if !self.begin("Searching…", cx) {
            return;
        }
        let fasta = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Use This FASTA".into()),
        });
        cx.spawn(async move |this, cx| {
            let Some(fasta) = picked_paths(fasta.await).and_then(|paths| paths.into_iter().next())
            else {
                return finish(&this, cx, None);
            };
            let Ok(spectra) = this.update(cx, |_, cx| {
                cx.prompt_for_paths(PathPromptOptions {
                    files: true,
                    // Bruker timsTOF runs are `.d` folders.
                    directories: true,
                    multiple: true,
                    prompt: Some("Search These Spectra".into()),
                })
            }) else {
                return;
            };
            let Some(spectra) = picked_paths(spectra.await) else {
                return finish(&this, cx, None);
            };
            let result = cx
                .background_executor()
                .spawn(async move { project::search(&path, &fasta, &spectra) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Undo the last import, search, or sample-group change (again: redo).
    fn undo_last_change(&mut self, _: &UndoLastChange, _: &mut Window, cx: &mut Context<Self>) {
        let Some(path) = self.project.as_ref().map(|report| report.path.clone()) else {
            return;
        };
        if !self.begin("Undoing…", cx) {
            return;
        }
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { project::undo(&path) })
                .await;
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Validate the open project; the last import's findings stay.
    fn validate(&mut self, cx: &mut Context<Self>) {
        let Some(current) = self.project.clone() else {
            return;
        };
        if !self.begin("Validating…", cx) {
            return;
        }
        let path = current.path.clone();
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { project::validate(&path) })
                .await
                .map(|report| ProjectReport {
                    findings: current.findings.clone(),
                    ..report
                });
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// Apply thresholds the user picked; the last import's findings stay.
    fn change_thresholds(&mut self, thresholds: Thresholds, cx: &mut Context<Self>) {
        let Some(current) = self.project.clone() else {
            return;
        };
        if !self.begin("Applying thresholds…", cx) {
            return;
        }
        let path = current.path.clone();
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { project::set_thresholds(&path, &thresholds) })
                .await
                .map(|report| ProjectReport {
                    findings: current.findings.clone(),
                    ..report
                });
            finish(&this, cx, Some(result));
        })
        .detach();
    }

    /// One bar: the project's name, the view tabs, and the project's
    /// actions as icon buttons (their names are in the tooltips and the
    /// menus).
    fn title_bar(&self, palette: &Palette, cx: &mut Context<Self>) -> impl IntoElement {
        let is_busy = self.blocked(cx);
        let validated = self.validated();
        let name = self.project.as_ref().map(|report| {
            report
                .path
                .file_stem()
                .map(|stem| stem.to_string_lossy().into_owned())
                .unwrap_or_else(|| report.path.display().to_string())
        });
        let focus = self.focus_handle.clone();
        div()
            .flex()
            .flex_none()
            .items_center()
            .gap_3()
            .h(px(40.))
            .px_3()
            .border_b_1()
            .border_color(palette.rule)
            .bg(palette.paper_raised)
            .child(
                div()
                    .text_size(px(13.))
                    .font_weight(FontWeight::SEMIBOLD)
                    .text_color(if name.is_some() {
                        palette.ink_faint
                    } else {
                        palette.ink
                    })
                    .child("Tuffold"),
            )
            .children(name.map(|name| {
                div()
                    .max_w(px(280.))
                    .overflow_hidden()
                    .whitespace_nowrap()
                    .text_ellipsis()
                    .text_size(px(13.))
                    .font_weight(FontWeight::MEDIUM)
                    .child(name)
            }))
            .children(self.browser.is_some().then(|| self.tabs(palette, cx)))
            .child(div().flex_1())
            .children(self.synthetic().then(|| {
                div()
                    .px_2()
                    .py(px(1.))
                    .rounded(px(3.))
                    .bg(palette.caution_soft)
                    .text_color(palette.caution)
                    .font_family(MONO)
                    .text_size(px(11.))
                    .child("SYNTHETIC DATA")
            }))
            .children(self.busy.filter(|_| self.project.is_some()).map(|label| {
                div()
                    .text_size(px(12.))
                    .text_color(palette.ink_muted)
                    .child(label)
            }))
            .children(self.project.is_some().then(|| {
                div()
                    .flex()
                    .items_center()
                    .gap_1()
                    .child(
                        Button::new("add-data")
                            .ghost()
                            .small()
                            .icon(IconName::Plus)
                            .tooltip_with_action("Add data", &AddData, None)
                            .disabled(is_busy)
                            .on_click(cx.listener(|this, _, window, cx| {
                                this.add_data(&AddData, window, cx)
                            })),
                    )
                    .child(
                        Button::new("search")
                            .ghost()
                            .small()
                            .icon(IconName::Search)
                            .tooltip_with_action("Search spectra", &SearchSpectra, None)
                            .disabled(is_busy)
                            .on_click(cx.listener(|this, _, window, cx| {
                                this.search_spectra(&SearchSpectra, window, cx)
                            })),
                    )
                    .child(
                        Button::new("export")
                            .ghost()
                            .small()
                            .label("Export")
                            .dropdown_caret(true)
                            .disabled(is_busy || !validated)
                            .dropdown_menu_with_anchor(Anchor::TopRight, move |menu, _, _| {
                                menu.action_context(focus.clone())
                                    .menu("Excel workbook…", Box::new(ExportWorkbook))
                                    .menu("Spectral library (.msp)…", Box::new(ExportLibrary))
                                    .menu("mzIdentML…", Box::new(ExportMzid))
                            }),
                    )
                    .child(
                        Button::new("open-other")
                            .ghost()
                            .small()
                            .icon(IconName::FolderOpen)
                            .tooltip_with_action("Open project", &OpenProject, None)
                            .disabled(is_busy)
                            .on_click(cx.listener(|this, _, window, cx| {
                                this.open_project(&OpenProject, window, cx)
                            })),
                    )
            }))
    }
}

impl Workspace {
    /// The views as a segmented control.
    fn tabs(&self, palette: &Palette, cx: &mut Context<Self>) -> impl IntoElement {
        let tabs = [
            (Tab::Overview, "Overview"),
            (Tab::Proteins, "Proteins"),
            (Tab::Samples, "Samples"),
            (Tab::Quant, "Quant"),
        ];
        let (raised, muted, ink) = (palette.paper_raised, palette.ink_muted, palette.ink);
        div()
            .flex()
            .items_center()
            .ml_2()
            .p(px(2.))
            .gap(px(2.))
            .rounded(px(6.))
            .bg(palette.paper_sunken)
            .children(tabs.into_iter().map(|(tab, label)| {
                let item = div()
                    .id(label)
                    .px(px(10.))
                    .h(px(24.))
                    .flex()
                    .items_center()
                    .rounded(px(4.))
                    .text_size(px(12.5))
                    .cursor_pointer()
                    .child(label);
                let item = if tab == self.tab {
                    item.bg(raised).font_weight(FontWeight::MEDIUM).shadow_xs()
                } else {
                    item.text_color(muted)
                        .hover(move |style| style.text_color(ink))
                };
                item.on_click(cx.listener(move |this, _, _, cx| this.set_tab(tab, cx)))
            }))
    }

    /// Whether the open project holds generated (demo) data.
    fn synthetic(&self) -> bool {
        self.project.as_ref().is_some_and(|report| {
            report
                .imports
                .iter()
                .any(|record| record.format == tuffold_core::SourceFormat::Synthetic.as_str())
        })
    }
}

fn browser_mode(tab: Tab) -> BrowserMode {
    match tab {
        Tab::Samples => BrowserMode::Samples,
        Tab::Overview | Tab::Proteins | Tab::Quant => BrowserMode::Proteins,
    }
}

/// The chosen paths, or `None` when the dialog was cancelled or failed.
fn picked_paths<E>(
    result: Result<anyhow::Result<Option<Vec<PathBuf>>>, E>,
) -> Option<Vec<PathBuf>> {
    result.ok()?.ok()?.filter(|paths| !paths.is_empty())
}

/// End a job; `None` means the user cancelled, which leaves the screen as it was.
fn finish(
    this: &WeakEntity<Workspace>,
    cx: &mut AsyncApp,
    result: Option<Result<ProjectReport, String>>,
) {
    // The window may have closed meanwhile; then there is nothing to update.
    let _ = this.update(cx, |this, cx| {
        this.busy = None;
        match result {
            None => {}
            Some(Ok(report)) => {
                let same_project = this
                    .project
                    .as_ref()
                    .is_some_and(|current| current.path == report.path);
                let validated = matches!(report.validation, ValidationState::Done(_));
                this.project = Some(Arc::new(report));
                this.failure = None;
                this.exported = None;
                match (&this.browser, same_project && validated) {
                    (Some(browser), true) => {
                        browser.update(cx, |browser, cx| browser.reload(cx));
                        if let Some(quant) = &this.quant {
                            quant.update(cx, |quant, cx| quant.reload(cx));
                        }
                    }
                    // Rebuilt at the next frame for a new or validated project.
                    _ => {
                        this.browser = None;
                        this._browser_events = None;
                        this.quant = None;
                        this._quant_events = None;
                        this._quant_writes = None;
                    }
                }
                // Show the overview when there is something to read there.
                if !validated || this.project.as_ref().is_some_and(|p| p.notice.is_some()) {
                    this.tab = Tab::Overview;
                }
            }
            // An open project stays on screen with the failure above it.
            Some(Err(message)) => this.failure = Some(message),
        }
        cx.notify();
    });
}

impl Workspace {
    /// A tab's content under the failure and export messages, which the
    /// overview shows in its own body; without this, exporting or failing
    /// from another tab would say nothing.
    fn with_banners(&self, content: Div, palette: &Palette) -> Div {
        div()
            .flex_1()
            .min_h_0()
            .flex()
            .flex_col()
            .children(
                self.failure
                    .as_deref()
                    .map(|message| div().pb_4().child(failure_banner(message, palette))),
            )
            .children(
                self.exported
                    .as_deref()
                    .map(|message| div().pb_4().child(notice_banner(message, palette))),
            )
            .child(div().flex_1().min_h_0().child(content))
    }
}

impl Render for Workspace {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let palette = *cx.global::<Palette>();
        if self.validated() && self.browser.is_none() {
            self.attach_browser(None, window, cx);
        }
        if self.validated() && self.tab == Tab::Quant && self.quant.is_none() {
            self.attach_quant(None, window, cx);
        }
        let workspace = cx.entity().downgrade();
        let on_thresholds: ThresholdHandler = Rc::new(move |thresholds, _, cx| {
            // The window may be closing; then there is nothing to change.
            let _ = workspace.update(cx, |this, cx| this.change_thresholds(thresholds, cx));
        });
        let workspace = cx.entity().downgrade();
        let on_validate: ValidateHandler = Rc::new(move |_, cx| {
            let _ = workspace.update(cx, |this, cx| this.validate(cx));
        });
        let controls = ValidationControls {
            on_thresholds: Some(on_thresholds),
            on_validate: Some(on_validate),
            busy: self.blocked(cx),
        };
        let body = match &self.project {
            Some(report) => div()
                .flex()
                .flex_col()
                .children(
                    self.failure
                        .as_deref()
                        .map(|message| failure_banner(message, &palette)),
                )
                .children(
                    report
                        .notice
                        .as_deref()
                        .map(|message| notice_banner(message, &palette)),
                )
                .children(
                    self.exported
                        .as_deref()
                        .map(|message| notice_banner(message, &palette)),
                )
                .children(self.offer_models.then(|| {
                    models_offer(
                        &palette,
                        self.blocked(cx),
                        ModelsOfferProps {
                            on_download: Box::new(|_, _, cx| {
                                cx.open_url(tuffold_engine::MODELS_URL);
                            }),
                            on_install: Box::new(
                                cx.listener(|this, _, _, cx| this.install_models(cx)),
                            ),
                            on_dismiss: Box::new(cx.listener(|this, _, _, cx| {
                                this.offer_models = false;
                                cx.notify();
                            })),
                        },
                    )
                }))
                .child(project_view(report, &controls, &palette))
                .into_any_element(),
            None => empty_state(
                &palette,
                EmptyStateProps {
                    busy: self.busy,
                    failure: self.failure.as_deref(),
                    on_new: Box::new(
                        cx.listener(|this, _, window, cx| {
                            this.new_project(&NewProject, window, cx)
                        }),
                    ),
                    on_open: Box::new(cx.listener(|this, _, window, cx| {
                        this.open_project(&OpenProject, window, cx)
                    })),
                    on_demo: Box::new(
                        cx.listener(|this, _, window, cx| this.open_demo(&OpenDemo, window, cx)),
                    ),
                },
            )
            .into_any_element(),
        };

        div()
            .id("workspace")
            .key_context("Workspace")
            .track_focus(&self.focus_handle)
            .on_action(cx.listener(Self::new_project))
            .on_action(cx.listener(Self::open_project))
            .on_action(cx.listener(Self::add_data))
            .on_action(cx.listener(Self::export_workbook))
            .on_action(cx.listener(Self::export_library))
            .on_action(cx.listener(Self::export_mzid))
            .on_action(cx.listener(Self::add_go))
            .on_action(cx.listener(Self::save_diagnostics))
            .on_action(cx.listener(Self::search_spectra))
            .on_action(cx.listener(Self::undo_last_change))
            .on_action(cx.listener(Self::open_demo))
            .on_action(cx.listener(Self::show_overview))
            .on_action(cx.listener(Self::show_proteins))
            .on_action(cx.listener(Self::show_samples))
            .on_action(cx.listener(Self::show_quant))
            .size_full()
            .flex()
            .flex_col()
            .bg(palette.paper)
            .text_color(palette.ink)
            .font_family(SANS)
            .text_size(px(15.))
            .child(self.title_bar(&palette, cx))
            .child(match (&self.browser, self.tab) {
                (Some(_), Tab::Quant) if self.quant.is_some() => self
                    .with_banners(div().size_full().children(self.quant.clone()), &palette)
                    .into_any_element(),
                (Some(browser), Tab::Proteins | Tab::Samples) => self
                    .with_banners(div().size_full().child(browser.clone()), &palette)
                    .into_any_element(),
                _ => div()
                    .id("main")
                    .flex_1()
                    .overflow_y_scroll()
                    .child(body)
                    .into_any_element(),
            })
    }
}
