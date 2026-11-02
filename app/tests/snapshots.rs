//! Headless visual snapshots of the main screens in both themes.
//!
//! Renders through GPUI's real renderer (Metal on macOS), so it checks what a
//! user would actually see: the bundled fonts, the palette, and the layout.
//! PNGs are written to `target/tmp/snapshots/` for review. Platforms without a
//! headless renderer skip, like gpui-kit's own rendering tests.

fn main() {
    #[cfg(target_os = "macos")]
    macos::run();
    #[cfg(not(target_os = "macos"))]
    println!("snapshots: skipped; GPUI supplies a headless renderer only on macOS");
}

#[cfg(target_os = "macos")]
mod macos {
    use std::path::{Path, PathBuf};
    use std::sync::Arc;

    use gpui_kit::component::ThemeMode;
    use gpui_kit::test::TestWindowExt as _;
    use gpui_kit::{
        AppContext as _, Bounds, HeadlessAppContext, WindowBounds, WindowOptions, px, size,
    };
    use image::RgbaImage;
    use tuffold_app::browse;
    use tuffold_app::project;
    use tuffold_app::quant;
    use tuffold_app::theme;
    use tuffold_app::workspace::{Tab, Workspace};

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Screen {
        Empty,
        Project,
        /// The fixture project with an illustrative validated summary, so the
        /// panel is seen populated (the fixtures themselves cannot be ranked).
        Validated,
        /// With no deep-learning models installed.
        ModelsOffer,
        /// The fixture project with results out of date (data imported since).
        Stale,
        /// Just after an undo: the notice above the overview.
        Undone,
        /// The demo project's linked views.
        Proteins,
        /// The same, thresholded by probability (Scaffold-style).
        ProteinsByProbability,
        /// A group sharing peptides with others: the similarity matrix.
        ProteinsSimilar,
        /// Filtered by a GO slim term, a group's terms under its header.
        ProteinsGo,
        /// A group with an oxidized multi-M peptide selected: site probabilities
        /// in the Peptides table and Spectrum header, and the localized residue
        /// badged on the fragment ladder.
        ProteinsLocalized,
        Samples,
        /// The Samples view showing emPAI.
        SamplesEmpai,
        Quant,
        /// The counts comparison by Mann–Whitney (its small-n caveat).
        QuantMannWhitney,
        QuantIntensity,
        QuantReporters,
        QuantTwoPlexes,
        QuantAll,
        QuantQc,
        OpenFailed,
    }

    /// The demo project, built once for every render that needs it.
    fn demo_path() -> &'static Path {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf)> = std::sync::OnceLock::new();
        let (_, path) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo.tuff");
            project::create_demo(&path).expect("demo builds");
            (dir, path)
        });
        path
    }

    fn demo_workspace(
        tab: Tab,
        window: &mut gpui_kit::Window,
        cx: &mut gpui_kit::Context<Workspace>,
    ) -> Workspace {
        let path = demo_path();
        let report = project::open(path).expect("demo opens");
        let data = browse::load(path, None).expect("demo browses");
        Workspace::with_browser(report, data, tab, window, cx)
    }

    /// A copy of the demo with one more passing group, led by a protein no
    /// group holds, sharing three peptides of the top group (the demo's
    /// isoforms are subsets, so none share peptides): the similarity matrix.
    fn similar_demo() -> (PathBuf, i64) {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf, i64)> =
            std::sync::OnceLock::new();
        let (_, path, group) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo-similar.tuff");
            std::fs::copy(demo_path(), &path).expect("demo copies");
            let top = tuffold_store::Project::open(&path)
                .expect("demo opens")
                .protein_list()
                .expect("demo lists")
                .rows[0]
                .group_id;
            let conn = rusqlite::Connection::open(&path).expect("demo opens");
            conn.execute_batch(&format!(
                "INSERT INTO protein_groups (id, leader_protein_id, is_decoy, score, q_value,
                                             peptides, pep)
                 SELECT 9000, id, 0, 1.0, 0.001, 3, 0.002 FROM proteins
                  WHERE is_decoy = 0 AND id NOT IN (SELECT protein_id FROM protein_group_members)
                  ORDER BY id LIMIT 1;
                 INSERT INTO protein_group_members (group_id, protein_id, role)
                 SELECT 9000, leader_protein_id, 'leader' FROM protein_groups WHERE id = 9000;
                 INSERT INTO protein_group_peptides (group_id, peptide_id, role)
                 SELECT 9000, peptide_id, 'shared' FROM protein_group_peptides
                  WHERE group_id = {top} ORDER BY peptide_id LIMIT 3;"
            ))
            .expect("the extra group is added");
            (dir, path, top)
        });
        (path.clone(), *group)
    }

    /// A copy of the demo whose proteins have UniProt-style accessions and GO
    /// slim terms (the demo's own accessions aren't UniProt's, so get none):
    /// every third target is cytosolic, the top group among them.
    fn go_demo() -> (PathBuf, i64) {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf, i64)> =
            std::sync::OnceLock::new();
        let (_, path, group) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo-go.tuff");
            std::fs::copy(demo_path(), &path).expect("demo copies");
            let conn = rusqlite::Connection::open(&path).expect("demo opens");
            // sp|DEMO0007|DEMO7_SYNTH -> sp|Q00007|DEMO7_SYNTH
            conn.execute_batch(
                "UPDATE proteins SET accession = replace(accession, 'sp|DEMO', 'sp|Q0')",
            )
            .expect("accessions rename");
            drop(conn);
            let mut project = tuffold_store::Project::open(&path).expect("demo opens");
            let top = project.protein_list().expect("demo lists").rows[0].group_id;
            let terms = [
                ("GO:0005829", "cytosol", "component"),
                ("GO:0005634", "nucleus", "component"),
                ("GO:0003824", "catalytic activity", "function"),
                (
                    "GO:0006091",
                    "generation of precursor metabolites and energy",
                    "process",
                ),
            ]
            .map(|(id, name, aspect)| (id.to_string(), name.to_string(), aspect.to_string()));
            let top_uniprot = tuffold_core::uniprot::uniprot_accession(
                &project.protein_list().expect("demo lists").rows[0].accession,
            )
            .expect("a UniProt-style accession")
            .to_string();
            let mut annotations = Vec::new();
            for (index, uniprot) in project
                .uniprot_accessions()
                .expect("ids")
                .iter()
                .enumerate()
            {
                let picked: &[usize] = match index % 3 {
                    _ if *uniprot == top_uniprot => &[0, 2, 3],
                    0 => &[0, 2, 3],
                    1 => &[1],
                    _ => &[1, 2],
                };
                annotations.extend(
                    picked
                        .iter()
                        .map(|&t| (uniprot.clone(), terms[t].0.clone())),
                );
            }
            project.set_go(&terms, &annotations).expect("GO saves");
            (dir, path, top)
        });
        (path.clone(), *group)
    }

    /// A copy of the demo thresholded at the probability defaults.
    fn probability_demo_path() -> &'static Path {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf)> = std::sync::OnceLock::new();
        let (_, path) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo-probability.tuff");
            std::fs::copy(demo_path(), &path).expect("demo copies");
            project::set_thresholds(&path, &tuffold_validate::Thresholds::probability_default())
                .expect("thresholds apply");
            (dir, path)
        });
        path
    }

    /// A copy of the demo with TMT6 reporter peaks on its identified spectra,
    /// a defined plex (`c1..c3` control, `t1..t3` treated), and every fifth
    /// protein group regulated 4× in treated.
    fn reporters_demo_path() -> &'static Path {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf)> = std::sync::OnceLock::new();
        let (_, path) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo-reporters.tuff");
            std::fs::copy(demo_path(), &path).expect("demo copies");
            let conn = rusqlite::Connection::open(&path).expect("demo opens");
            conn.execute_batch("UPDATE spectra SET peaks_partial = 1")
                .expect("peaks marked replaceable");
            let spectra: Vec<(String, String, i64)> = conn
                .prepare(
                    "SELECT m.name, s.spectrum_key, MIN(gp.group_id)
                     FROM spectra s
                     JOIN ms_runs m ON m.id = s.ms_run_id
                     JOIN spectrum_results sr ON sr.spectrum_id = s.id
                     JOIN protein_group_peptides gp ON gp.peptide_id = sr.peptide_id
                     GROUP BY m.name, s.spectrum_key
                     ORDER BY m.name, s.id",
                )
                .expect("query")
                .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
                .expect("rows")
                .collect::<Result<_, _>>()
                .expect("spectra");
            drop(conn);
            let mut store = tuffold_store::Project::open(&path).expect("demo opens");
            let channels = tuffold_core::reporter::Plex::Tmt6.channels();
            let mut by_run: std::collections::BTreeMap<String, Vec<(String, i64)>> =
                std::collections::BTreeMap::new();
            for (run, key, group_id) in spectra {
                by_run.entry(run).or_default().push((key, group_id));
            }
            for (run, list) in by_run {
                let mut attacher = store.attach_spectra(&run).expect("attacher");
                for (i, (key, group_id)) in list.into_iter().enumerate() {
                    let fold = match group_id % 6 {
                        0 => 4.0,
                        1 => 0.25,
                        _ => 1.0,
                    };
                    let mut mz = Vec::with_capacity(channels.len());
                    let mut intensity = Vec::with_capacity(channels.len());
                    for (c, channel) in channels.iter().enumerate() {
                        let truth = if c >= 3 { fold } else { 1.0 };
                        let loading = if channel.name == "131" { 1.5 } else { 1.0 };
                        let noise = 1.0 + 0.04 * ((i * 7 + c * 13) as f64).sin();
                        mz.push(channel.mz);
                        intensity.push((10_000.0 * truth * loading * noise) as f32);
                    }
                    attacher
                        .add(&tuffold_core::PeakSpectrum {
                            key: tuffold_core::SpectrumKey::from_native_id(&key),
                            native_id: key,
                            ms_level: 2,
                            rt_seconds: None,
                            precursor_mz: None,
                            precursor_charge: None,
                            mz,
                            intensity,
                            profile: false,
                            parent: None,
                            isolation: None,
                            precursor_mobility: None,
                            mobility: Vec::new(),
                        })
                        .expect("spectrum attaches");
                }
                attacher.finish().expect("attach finishes");
            }
            drop(store);
            let samples = ["c1", "c2", "c3", "t1", "t2", "t3"];
            quant::save_plex(
                &path,
                tuffold_core::reporter::Plex::Tmt6,
                &channels
                    .iter()
                    .zip(samples)
                    .map(|(c, s)| (c.name.clone(), (*s).to_string()))
                    .collect::<Vec<_>>(),
                &[],
                quant::Certificate::Keep,
            )
            .expect("plex saves");
            quant::save_categories(
                &path,
                &samples
                    .iter()
                    .map(|s| {
                        (
                            (*s).to_string(),
                            if s.starts_with('c') {
                                "control".into()
                            } else {
                                "treated".into()
                            },
                        )
                    })
                    .collect::<Vec<_>>(),
            )
            .expect("categories save");
            (dir, path)
        });
        path
    }

    /// The reporter demo split into two TMT6 plexes (`pool-1` over the
    /// control runs, `pool-2` over the treated runs), bridged through their
    /// `126` reference channel.
    fn two_plexes_demo_path() -> &'static Path {
        static DEMO: std::sync::OnceLock<(tempfile::TempDir, PathBuf)> = std::sync::OnceLock::new();
        let (_, path) = DEMO.get_or_init(|| {
            let dir = tempfile::tempdir().expect("temp dir");
            let path = dir.path().join("demo-two-plexes.tuff");
            std::fs::copy(reporters_demo_path(), &path).expect("demo copies");
            quant::remove_plex(&path, "plex").expect("default plex removes");
            let channels = tuffold_core::reporter::Plex::Tmt6.channels();
            let mut categories = Vec::new();
            for (plex, runs, n) in [
                ("pool-1", ["control-A", "control-B"], [1, 2]),
                ("pool-2", ["treated-A", "treated-B"], [3, 4]),
            ] {
                let samples = [
                    format!("{plex}-ref"),
                    format!("c{}", n[0]),
                    format!("c{}", n[1]),
                    format!("t{}", n[0]),
                    format!("t{}", n[1]),
                ];
                let assigned: Vec<(String, String)> = channels
                    .iter()
                    .zip(&samples)
                    .map(|(c, s)| (c.name.clone(), s.clone()))
                    .collect();
                quant::save_plex_named(
                    &path,
                    plex,
                    true,
                    tuffold_core::reporter::Plex::Tmt6,
                    &runs.map(str::to_string),
                    &assigned,
                    Some(&["126".into()]),
                    quant::Certificate::Keep,
                )
                .expect("plex saves");
                for sample in &samples[1..] {
                    categories.push((
                        sample.clone(),
                        if sample.starts_with('c') {
                            "control".into()
                        } else {
                            "treated".into()
                        },
                    ));
                }
            }
            quant::save_categories(&path, &categories).expect("categories save");
            (dir, path)
        });
        path
    }

    fn validated_project(dir: &Path) -> project::ProjectReport {
        use tuffold_store::{EngineRecord, LevelCounts, ValidationSummary};
        use tuffold_validate::{Settings, Thresholds};
        let engine = |name: &str, psms, decoys, accepted| EngineRecord {
            name: name.into(),
            psms,
            decoys,
            method: "rescored".into(),
            detail: "Rescored by cross-validated logistic regression over 9 features.".into(),
            features: Vec::new(),
            accepted,
        };
        let summary = ValidationSummary {
            settings: Settings {
                rescore: true,
                thresholds: Thresholds {
                    min_peptides: 2,
                    ..Thresholds::default()
                },
            },
            validated_at: Some(0),
            engines: vec![
                engine("MSFragger", 48_210, 9_102, 21_877),
                engine("Comet", 47_985, 9_844, 19_406),
            ],
            peptide_scoring: None,
            notes: Vec::new(),
            spectra: LevelCounts {
                targets: 23_512,
                decoys: 231,
                lowest_q: Some(0.000_04),
            },
            peptides: LevelCounts {
                targets: 14_087,
                decoys: 139,
                lowest_q: Some(0.000_07),
            },
            protein_groups: LevelCounts {
                targets: 2_964,
                decoys: 28,
                lowest_q: Some(0.000_3),
            },
            localization: tuffold_store::LocalizationSummary {
                scored_psms: 1_842,
                localized_psms: 1_519,
                missing_peaks_psms: 0,
                capped_psms: 0,
                unscorable_psms: 0,
                localized_sites: 1_519,
                total_target_sites: 1_842,
                probability_threshold: 0.75,
                estimated_flr: Some(0.008_4),
            },
        };
        project::ProjectReport {
            validation: project::ValidationState::Done(Box::new(summary)),
            ..fixture_project(dir)
        }
    }

    /// Drill into the first demo protein group that has a localized multi-M
    /// oxidized peptide, with that peptide selected so its localized residue
    /// and site probabilities render in both the table and the spectrum view.
    fn localized_demo_data(path: &Path) -> browse::BrowserData {
        let store = tuffold_store::Project::open(path).expect("demo opens");
        let list = store.protein_list().expect("demo lists");
        for row in &list.rows {
            let Some(mut group) = store.group_detail(row.group_id).expect("group loads") else {
                continue;
            };
            if let Some(pos) = group
                .peptides
                .iter()
                .position(|p| p.site_probabilities.is_some())
            {
                let pep = group.peptides.remove(pos);
                let peptide_id = pep.peptide_id;
                group.peptides.insert(0, pep);
                let (spectra, spectrum) =
                    browse::load_peptide(path, peptide_id).expect("peptide loads");
                return browse::BrowserData {
                    list,
                    drill: browse::Drill {
                        peptide_id: Some(peptide_id),
                        group: Some(group),
                        spectra,
                        spectrum,
                    },
                };
            }
        }
        panic!("the demo has localized multi-M oxidized peptides");
    }

    /// A project built from every fixture: FragPipe, pepXML, mzIdentML, and
    /// an MGF whose run ("test") none of them created, so it fails visibly.
    fn fixture_project(dir: &Path) -> project::ProjectReport {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures");
        let inputs = [
            "fragpipe_synthetic",
            "third_party/pyteomics/test.pep.xml",
            "third_party/pyteomics/test.mzid",
            "third_party/pyteomics/test.mgf",
        ]
        .map(|relative| {
            fixtures
                .join(relative)
                .canonicalize()
                .expect("fixture exists")
        });
        let mut report = project::create(&dir.join("study.tuff"), &inputs).expect("project builds");
        // Show a stable path in the snapshot regardless of the temp dir.
        report.path = PathBuf::from("/projects/study.tuff");
        report
    }

    fn render(screen: Screen, mode: ThemeMode) -> RgbaImage {
        let mut cx = HeadlessAppContext::with_platform(
            gpui_kit::platform::current_platform(true).text_system(),
            Arc::new(gpui_kit::assets::Assets),
            gpui_kit::platform::current_headless_renderer,
        );
        cx.update(|cx| {
            gpui_kit::init(cx);
            theme::load_fonts(cx).expect("bundled fonts load");
            theme::apply(mode, cx);
        });
        let scratch = tempfile::tempdir().expect("temp dir");
        let options = WindowOptions {
            window_bounds: Some(WindowBounds::Windowed(Bounds {
                origin: Default::default(),
                size: size(px(1280.), px(820.)),
            })),
            focus: false,
            show: false,
            ..Default::default()
        };
        let (handle, _) = cx
            .update(|cx| {
                gpui_kit::open_window(options, cx, |window, cx| {
                    cx.new(|cx| match screen {
                        Screen::Empty => Workspace::new(window, cx),
                        Screen::Project => {
                            Workspace::with_project(fixture_project(scratch.path()), window, cx)
                        }
                        Screen::Stale => Workspace::with_project(
                            project::ProjectReport {
                                validation: project::ValidationState::Pending { stale: true },
                                ..fixture_project(scratch.path())
                            },
                            window,
                            cx,
                        ),
                        Screen::Proteins => demo_workspace(Tab::Proteins, window, cx),
                        Screen::ProteinsByProbability => {
                            let path = probability_demo_path();
                            Workspace::with_browser(
                                project::open(path).expect("demo opens"),
                                browse::load(path, None).expect("demo browses"),
                                Tab::Proteins,
                                window,
                                cx,
                            )
                        }
                        Screen::Samples => demo_workspace(Tab::Samples, window, cx),
                        Screen::ProteinsSimilar => {
                            let (path, group) = similar_demo();
                            Workspace::with_browser(
                                project::open(&path).expect("demo opens"),
                                browse::load(&path, Some(group)).expect("demo browses"),
                                Tab::Proteins,
                                window,
                                cx,
                            )
                        }
                        Screen::ProteinsGo => {
                            let (path, group) = go_demo();
                            let mut workspace = Workspace::with_browser(
                                project::open(&path).expect("demo opens"),
                                browse::load(&path, Some(group)).expect("demo browses"),
                                Tab::Proteins,
                                window,
                                cx,
                            );
                            workspace.set_go_filter(Some("cytosol (component)".into()), cx);
                            workspace
                        }
                        Screen::ProteinsLocalized => {
                            let path = demo_path();
                            Workspace::with_browser(
                                project::open(path).expect("demo opens"),
                                localized_demo_data(path),
                                Tab::Proteins,
                                window,
                                cx,
                            )
                        }
                        Screen::SamplesEmpai => {
                            let mut workspace = demo_workspace(Tab::Samples, window, cx);
                            workspace.set_sample_measure(browse::SampleMeasure::Empai, cx);
                            workspace
                        }
                        Screen::Quant
                        | Screen::QuantMannWhitney
                        | Screen::QuantIntensity
                        | Screen::QuantReporters
                        | Screen::QuantTwoPlexes
                        | Screen::QuantAll
                        | Screen::QuantQc => {
                            let path = match screen {
                                Screen::QuantReporters => reporters_demo_path(),
                                Screen::QuantTwoPlexes => two_plexes_demo_path(),
                                _ => demo_path(),
                            };
                            let method = match screen {
                                Screen::Quant | Screen::QuantMannWhitney => quant::Method::Counts,
                                Screen::QuantReporters | Screen::QuantTwoPlexes => {
                                    quant::Method::Reporters
                                }
                                _ => quant::Method::Intensity,
                            };
                            let scope = match screen {
                                Screen::QuantAll => quant::Scope::All,
                                Screen::QuantQc => quant::Scope::Qc,
                                _ => quant::Scope::Pair,
                            };
                            Workspace::with_quant(
                                project::open(path).expect("demo opens"),
                                browse::load(path, None).expect("demo browses"),
                                quant::load(path).expect("demo quantifies"),
                                quant::QuantChoice {
                                    scope,
                                    method,
                                    // Fisher, so significant points and labels render.
                                    test: if screen == Screen::QuantMannWhitney {
                                        tuffold_quant::PrimaryTest::MannWhitney
                                    } else {
                                        tuffold_quant::PrimaryTest::Fisher
                                    },
                                },
                                window,
                                cx,
                            )
                        }
                        Screen::ModelsOffer => {
                            let mut workspace = Workspace::with_project(
                                validated_project(scratch.path()),
                                window,
                                cx,
                            );
                            workspace.show_models_offer();
                            workspace
                        }
                        Screen::Validated => {
                            Workspace::with_project(validated_project(scratch.path()), window, cx)
                        }
                        Screen::Undone => Workspace::with_project(
                            project::ProjectReport {
                                notice: Some(
                                    "Undid “import test.mgf”. Undo again to redo it.".into(),
                                ),
                                ..validated_project(scratch.path())
                            },
                            window,
                            cx,
                        ),
                        Screen::OpenFailed => Workspace::with_failure(
                            "notes.tuff is not a Tuffold project (.tuff) file".into(),
                            window,
                            cx,
                        ),
                    })
                })
            })
            .expect("window opens");
        cx.update_window(handle, |_, window, cx| window.render_frame(cx))
            .expect("frame renders");
        cx.capture_screenshot(handle)
            .expect("Metal headless rendering is available")
    }

    fn out_dir() -> PathBuf {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("snapshots");
        std::fs::create_dir_all(&dir).expect("snapshot dir");
        dir
    }

    /// Fraction of pixels that differ from the top-left pixel: a blank or
    /// failed render is (nearly) uniform.
    fn non_background_fraction(image: &RgbaImage) -> f64 {
        let background = image.get_pixel(0, image.height() - 1);
        let differing = image.pixels().filter(|pixel| *pixel != background).count();
        differing as f64 / (image.width() * image.height()) as f64
    }

    pub fn run() {
        let screens = [
            ("empty", Screen::Empty),
            ("project", Screen::Project),
            ("validated", Screen::Validated),
            ("models-offer", Screen::ModelsOffer),
            ("undone", Screen::Undone),
            ("stale", Screen::Stale),
            ("proteins", Screen::Proteins),
            ("proteins-probability", Screen::ProteinsByProbability),
            ("proteins-similar", Screen::ProteinsSimilar),
            ("proteins-go", Screen::ProteinsGo),
            ("proteins-localized", Screen::ProteinsLocalized),
            ("samples", Screen::Samples),
            ("samples-empai", Screen::SamplesEmpai),
            ("quant", Screen::Quant),
            ("quant-mann-whitney", Screen::QuantMannWhitney),
            ("quant-intensity", Screen::QuantIntensity),
            ("quant-reporters", Screen::QuantReporters),
            ("quant-two-plexes", Screen::QuantTwoPlexes),
            ("quant-all", Screen::QuantAll),
            ("quant-qc", Screen::QuantQc),
            ("open-failed", Screen::OpenFailed),
        ];
        for (name, screen) in screens {
            let light = render(screen, ThemeMode::Light);
            let dark = render(screen, ThemeMode::Dark);
            for (mode, image) in [("light", &light), ("dark", &dark)] {
                let path = out_dir().join(format!("{name}-{mode}.png"));
                image.save(&path).expect("snapshot saves");
                let ink = non_background_fraction(image);
                assert!(
                    ink > 0.005,
                    "{name}-{mode} looks blank ({ink:.4} non-background)"
                );
                println!(
                    "snapshot {name}-{mode}: {}x{} -> {}",
                    image.width(),
                    image.height(),
                    path.display()
                );
            }
            // Plain asserts: a failing assert_eq! would print megabytes of pixels.
            assert!(
                light.as_raw() != dark.as_raw(),
                "{name}: light and dark render identically"
            );
            assert!(
                render(screen, ThemeMode::Light).as_raw() == light.as_raw(),
                "{name}: rendering is not deterministic"
            );
        }
        println!("snapshots: ok");
    }
}
