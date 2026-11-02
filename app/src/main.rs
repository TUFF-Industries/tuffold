//! Tuffold desktop application (GPUI).

use gpui_kit::{
    App, AppContext as _, Bounds, KeyBinding, Menu, MenuItem, WindowBounds, WindowOptions, px, size,
};

use tuffold_app::theme;
use tuffold_app::workspace::{
    AddData, AddGoAnnotations, ExportLibrary, ExportMzid, ExportWorkbook, NewProject, OpenDemo,
    OpenProject, Quit, SaveDiagnostics, SearchSpectra, ShowOverview, ShowProteins, ShowQuant,
    ShowSamples, UndoLastChange, Workspace,
};

fn main() {
    tuffold_engine::install_panic_log("tuffold-app");
    gpui_kit::application().run(|cx: &mut App| {
        gpui_kit::init(cx);
        if let Err(error) = theme::load_fonts(cx) {
            eprintln!("warning: bundled fonts failed to load ({error}); using system fonts");
        }
        cx.bind_keys([
            KeyBinding::new("secondary-n", NewProject, None),
            KeyBinding::new("secondary-o", OpenProject, None),
            KeyBinding::new("secondary-i", AddData, None),
            KeyBinding::new("secondary-shift-f", SearchSpectra, None),
            KeyBinding::new("secondary-z", UndoLastChange, None),
            KeyBinding::new("secondary-1", ShowOverview, None),
            KeyBinding::new("secondary-2", ShowProteins, None),
            KeyBinding::new("secondary-3", ShowSamples, None),
            KeyBinding::new("secondary-4", ShowQuant, None),
            KeyBinding::new("secondary-q", Quit, None),
        ]);
        cx.on_action(|_: &Quit, cx| cx.quit());
        cx.set_menus([Menu::new("Tuffold").items([
            MenuItem::action("New Project from Data…", NewProject),
            MenuItem::action("Open Project…", OpenProject),
            MenuItem::action("Add Data…", AddData),
            MenuItem::action("Search Spectra…", SearchSpectra),
            MenuItem::action("Add GO Annotations…", AddGoAnnotations),
            MenuItem::action("Export Excel Workbook…", ExportWorkbook),
            MenuItem::action("Export Spectral Library…", ExportLibrary),
            MenuItem::action("Export mzIdentML…", ExportMzid),
            MenuItem::action("Undo Last Change", UndoLastChange),
            MenuItem::action("Open Demo Project…", OpenDemo),
            MenuItem::action("Save Diagnostics…", SaveDiagnostics),
            MenuItem::separator(),
            MenuItem::action("Overview", ShowOverview),
            MenuItem::action("Proteins", ShowProteins),
            MenuItem::action("Samples", ShowSamples),
            MenuItem::action("Quant", ShowQuant),
            MenuItem::separator(),
            MenuItem::action("Quit Tuffold", Quit),
        ])]);

        let options = WindowOptions {
            window_bounds: Some(WindowBounds::Windowed(Bounds::centered(
                None,
                size(px(1280.), px(820.)),
                cx,
            ))),
            window_min_size: Some(size(px(900.), px(600.))),
            ..Default::default()
        };
        let opened = gpui_kit::open_window(options, cx, |window, cx| {
            window.set_window_title("Tuffold");
            theme::sync_with_system(window, cx);
            window
                .observe_window_appearance(theme::sync_with_system)
                .detach();
            cx.new(|cx| {
                let mut workspace = Workspace::new(window, cx);
                workspace.check_models(cx);
                workspace
            })
        });
        if let Err(error) = opened {
            eprintln!("error: could not open the Tuffold window: {error}");
            cx.quit();
        }
        cx.activate(true);
    });
}
