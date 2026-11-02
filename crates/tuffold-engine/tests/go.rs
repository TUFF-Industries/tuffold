//! GO annotations on a project: a demo copy whose top group's leader gets a
//! UniProt accession, and hand-written ontology, slim, and GAF files.

use std::path::PathBuf;

use tuffold_engine::{GoFiles, annotate_go, create_demo_project};
use tuffold_store::Project;

const ONTOLOGY: &str = "format-version: 1.2

[Term]
id: GO:0007049
name: cell cycle
namespace: biological_process

[Term]
id: GO:0000082
name: G1/S transition of mitotic cell cycle
namespace: biological_process
relationship: part_of GO:0007049 ! cell cycle

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
";

const SLIM: &str = "format-version: 1.2
subsetdef: goslim_generic \"Generic GO slim\"

[Term]
id: GO:0007049
name: cell cycle
namespace: biological_process

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
";

// P12345: G1/S transition (reaches cell cycle); NOT in cytoplasm (skipped);
// an unknown term; Q99999 is not in the project.
const GAF: &str = "!gaf-version: 2.2
UniProtKB\tP12345\tCDK1\tinvolved_in\tGO:0000082\tPMID:1\tIDA\t\tP
UniProtKB\tP12345\tCDK1\tNOT|located_in\tGO:0005737\tPMID:1\tIDA\t\tC
UniProtKB\tP12345\tCDK1\tenables\tGO:9999999\tPMID:1\tIDA\t\tF
UniProtKB\tQ99999\tXYZ\tlocated_in\tGO:0005737\tPMID:1\tIDA\t\tC
";

#[test]
fn proteins_get_the_slim_terms_their_annotations_reach() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let leader = Project::open(&path).unwrap().protein_list().unwrap().rows[0]
        .accession
        .clone();
    rusqlite::Connection::open(&path)
        .unwrap()
        .execute(
            "UPDATE proteins SET accession = 'sp|P12345|CDK1_HUMAN' WHERE accession = ?1",
            [&leader],
        )
        .unwrap();
    let write = |name: &str, text: &str| -> PathBuf {
        let file = dir.path().join(name);
        std::fs::write(&file, text).unwrap();
        file
    };
    let files = GoFiles::classify(&[
        write("annotations.gaf", GAF),
        write("go-basic.obo", ONTOLOGY),
        write("goslim_generic.obo", SLIM),
    ])
    .unwrap();
    let mut project = Project::open(&path).unwrap();
    let report = annotate_go(&mut project, &files).unwrap();
    assert_eq!(report.uniprot_proteins, 1);
    assert_eq!(report.annotated, 1);
    assert_eq!(report.terms, 1);
    assert_eq!(report.unknown_terms, 1);
    let rows = project.protein_list().unwrap().rows;
    let top = rows
        .iter()
        .find(|r| r.accession == "sp|P12345|CDK1_HUMAN")
        .unwrap();
    assert_eq!(top.go, ["cell cycle (process)"]);
    assert!(
        rows.iter()
            .filter(|r| r.group_id != top.group_id)
            .all(|r| r.go.is_empty())
    );
}

#[test]
fn picked_files_sort_into_ontology_slim_and_annotations() {
    let paths = |names: &[&str]| names.iter().map(PathBuf::from).collect::<Vec<_>>();
    let files = GoFiles::classify(&paths(&[
        "goa_human.gaf.gz",
        "go-basic.obo",
        "goslim_generic.obo",
        "goa_yeast.gaf.gz",
    ]))
    .unwrap();
    assert_eq!(files.ontology, PathBuf::from("go-basic.obo"));
    assert_eq!(files.slim, PathBuf::from("goslim_generic.obo"));
    assert_eq!(files.annotations.len(), 2);
    let missing = GoFiles::classify(&paths(&["goa_human.gaf.gz", "go-basic.obo"])).unwrap_err();
    assert!(missing.contains("GO slim"), "{missing}");
    let stray = GoFiles::classify(&paths(&["notes.txt"])).unwrap_err();
    assert!(stray.contains("notes.txt"), "{stray}");
}
