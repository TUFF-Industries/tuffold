//! Undo and redo of imports and sample-group edits.

use std::path::{Path, PathBuf};

use tuffold_engine::{UndoError, import_inputs, last_change, set_categories, undo, undo_path};
use tuffold_store::Project;

fn fixture(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(relative)
}

fn imports(path: &Path) -> usize {
    Project::open(path).unwrap().imports().unwrap().len()
}

#[test]
fn an_import_undoes_and_redoes() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    drop(Project::create(&path).unwrap());
    assert!(matches!(undo(&path), Err(UndoError::NothingToUndo(_))));

    let pepxml = fixture("third_party/pyteomics/test.pep.xml");
    import_inputs(
        &mut Project::open(&path).unwrap(),
        std::slice::from_ref(&pepxml),
        None,
    );
    assert_eq!(imports(&path), 1);
    let pending = last_change(&path).unwrap().unwrap();
    assert_eq!(
        (pending.label.as_str(), pending.redo),
        ("import test.pep.xml", false)
    );

    let undone = undo(&path).unwrap();
    assert_eq!(undone.label, "import test.pep.xml");
    assert_eq!(imports(&path), 0);
    // The project itself carries no undo note; the copy now offers a redo.
    assert_eq!(Project::open(&path).unwrap().undo_note().unwrap(), None);
    assert!(last_change(&path).unwrap().unwrap().redo);

    let redone = undo(&path).unwrap();
    assert!(redone.redo);
    assert_eq!(imports(&path), 1);
    assert!(!last_change(&path).unwrap().unwrap().redo);
    // No temporary files are left behind.
    let mut names: Vec<String> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(names, ["p.tuff", "p.tuff.undo"]);
}

#[test]
fn a_failed_import_keeps_the_previous_undo_point() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    drop(Project::create(&path).unwrap());
    let pepxml = fixture("third_party/pyteomics/test.pep.xml");
    import_inputs(
        &mut Project::open(&path).unwrap(),
        std::slice::from_ref(&pepxml),
        None,
    );
    // The same file again is refused: nothing changed, nothing to undo.
    let again = import_inputs(&mut Project::open(&path).unwrap(), &[pepxml], None);
    assert!(again.iter().all(|outcome| outcome.is_failure()));
    let pending = last_change(&path).unwrap().unwrap();
    assert_eq!(pending.label, "import test.pep.xml");
    undo(&path).unwrap();
    assert_eq!(imports(&path), 0);
}

#[test]
fn sample_groups_undo() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    tuffold_engine::create_demo_project(&path).unwrap();
    let category = |name: &str| {
        Project::open(&path)
            .unwrap()
            .protein_list()
            .unwrap()
            .biosamples
            .into_iter()
            .find(|b| b.name == name)
            .unwrap()
            .category
    };
    assert_eq!(category("treated-B").as_deref(), Some("treated"));
    set_categories(&mut Project::open(&path).unwrap(), &[("treated-B", None)]).unwrap();
    assert_eq!(category("treated-B"), None);
    assert_eq!(undo(&path).unwrap().label, "change sample groups");
    assert_eq!(category("treated-B").as_deref(), Some("treated"));
    // An unknown biosample changes nothing, and keeps the redo.
    let unknown = set_categories(&mut Project::open(&path).unwrap(), &[("nobody", None)]);
    assert!(unknown.is_err());
    assert!(last_change(&path).unwrap().unwrap().redo);
}

#[test]
fn an_interrupted_project_is_not_swapped() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    drop(Project::create(&path).unwrap());
    import_inputs(
        &mut Project::open(&path).unwrap(),
        &[fixture("third_party/pyteomics/test.pep.xml")],
        None,
    );
    let journal = dir.path().join("p.tuff-journal");
    std::fs::write(&journal, b"").unwrap();
    assert!(matches!(undo(&path), Err(UndoError::Interrupted(_))));
    assert_eq!(imports(&path), 1);
    assert!(undo_path(&path).exists());
}

#[test]
fn an_undo_copy_from_a_replaced_project_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    drop(Project::create(&path).unwrap());
    import_inputs(
        &mut Project::open(&path).unwrap(),
        &[fixture("third_party/pyteomics/test.pep.xml")],
        None,
    );
    // The project is deleted and a new one made under the same name.
    std::fs::remove_file(&path).unwrap();
    drop(Project::create(&path).unwrap());
    let refused = undo(&path).unwrap_err();
    assert!(matches!(refused, UndoError::Unrelated { .. }), "{refused}");
    assert!(refused.to_string().contains("delete"), "{refused}");
}

#[test]
fn an_interrupted_undo_is_reported_not_repeated() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    drop(Project::create(&path).unwrap());
    import_inputs(
        &mut Project::open(&path).unwrap(),
        &[fixture("third_party/pyteomics/test.pep.xml")],
        None,
    );
    std::fs::copy(&path, dir.path().join("p.tuff.undo.swap")).unwrap();
    let refused = undo(&path).unwrap_err();
    assert!(matches!(refused, UndoError::Unfinished { .. }), "{refused}");
    assert_eq!(imports(&path), 1);
}

#[test]
fn changes_that_change_nothing_keep_the_undo_point() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    tuffold_engine::create_demo_project(&path).unwrap();
    set_categories(&mut Project::open(&path).unwrap(), &[("treated-B", None)]).unwrap();
    // Already so: no new checkpoint, the real change stays undoable.
    set_categories(
        &mut Project::open(&path).unwrap(),
        &[("treated-B", Some(" "))],
    )
    .unwrap();
    assert_eq!(
        last_change(&path).unwrap().unwrap().label,
        "change sample groups"
    );
    let saved = std::fs::metadata(undo_path(&path))
        .unwrap()
        .modified()
        .unwrap();
    // Spectra for runs the project doesn't have attach nothing.
    let outcomes = import_inputs(
        &mut Project::open(&path).unwrap(),
        &[fixture("third_party/pyteomics/test.mgf")],
        None,
    );
    assert!(outcomes.iter().all(|o| o.is_failure()), "{outcomes:?}");
    let still = std::fs::metadata(undo_path(&path))
        .unwrap()
        .modified()
        .unwrap();
    assert_eq!(saved, still);
    undo(&path).unwrap();
}

#[test]
fn removing_a_plex_undoes_and_an_unknown_plex_keeps_the_undo_point() {
    use tuffold_core::reporter::Plex;
    use tuffold_engine::{remove_plex, set_plex};
    use tuffold_store::PlexDefinition;

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    tuffold_engine::create_demo_project(&path).unwrap();
    let definition = PlexDefinition {
        name: "pool-1".into(),
        plex: Plex::Tmt6,
        runs: vec!["control-A".into(), "control-B".into()],
        channels: vec![("126".into(), "c1".into()), ("127".into(), "t1".into())],
        impurities: None,
        references: vec!["126".into()],
    };
    set_plex(&mut Project::open(&path).unwrap(), &definition).unwrap();
    assert_eq!(last_change(&path).unwrap().unwrap().label, "define a plex");

    // An unknown name fails without overwriting the undo point.
    assert!(remove_plex(&mut Project::open(&path).unwrap(), "missing").is_err());
    assert_eq!(last_change(&path).unwrap().unwrap().label, "define a plex");

    remove_plex(&mut Project::open(&path).unwrap(), " pool-1 ").unwrap();
    assert!(Project::open(&path).unwrap().plexes().unwrap().is_empty());
    assert_eq!(last_change(&path).unwrap().unwrap().label, "remove a plex");

    // Undo brings the plex back, and undo again removes it again.
    let undone = undo(&path).unwrap();
    assert_eq!(undone.label, "remove a plex");
    assert_eq!(
        Project::open(&path).unwrap().plexes().unwrap(),
        [definition]
    );
    let redone = undo(&path).unwrap();
    assert!(redone.redo);
    assert!(Project::open(&path).unwrap().plexes().unwrap().is_empty());
}
