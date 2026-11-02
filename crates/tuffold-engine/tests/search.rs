//! The built-in search end to end: the demo's proteins as a FASTA and its
//! spectra as MGF files, searched with Sage, imported, and validated.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use tuffold_core::PeakSpectrum;
use tuffold_engine::demo::generate;
use tuffold_engine::{SearchSettings, search_files, validate_project};
use tuffold_store::Project;
use tuffold_validate::Settings;

fn write_fasta(dir: &Path, demo: &tuffold_engine::demo::Demo) -> PathBuf {
    let text = demo
        .batch
        .proteins
        .iter()
        .filter_map(|p| Some(format!(">{}\n{}\n", p.accession, p.sequence.as_ref()?)))
        .collect::<String>();
    let path = dir.join("demo.fasta");
    std::fs::write(&path, text).unwrap();
    path
}

fn write_mgf(dir: &Path, run: &str, spectra: &[PeakSpectrum]) -> PathBuf {
    let mut text = String::new();
    for s in spectra {
        let _ = writeln!(text, "BEGIN IONS\nTITLE={}", s.native_id);
        if let Some(mz) = s.precursor_mz {
            let _ = writeln!(text, "PEPMASS={mz}");
        }
        if let Some(z) = s.precursor_charge {
            let _ = writeln!(text, "CHARGE={z}+");
        }
        if let Some(rt) = s.rt_seconds {
            let _ = writeln!(text, "RTINSECONDS={rt}");
        }
        for (mz, intensity) in s.mz.iter().zip(&s.intensity) {
            let _ = writeln!(text, "{mz} {intensity}");
        }
        text.push_str("END IONS\n");
    }
    let path = dir.join(format!("{run}.mgf"));
    std::fs::write(&path, text).unwrap();
    path
}

#[test]
fn searching_the_demo_spectra_finds_the_expressed_proteins() {
    let dir = tempfile::tempdir().unwrap();
    let demo = generate(7);
    let fasta = write_fasta(dir.path(), &demo);
    // Two of the four runs keep the test quick.
    let files: Vec<PathBuf> = demo
        .peaks
        .iter()
        .take(2)
        .map(|(run, scans)| {
            let msms: Vec<PeakSpectrum> =
                scans.iter().filter(|s| s.ms_level >= 2).cloned().collect();
            write_mgf(dir.path(), run, &msms)
        })
        .collect();

    let path = dir.path().join("searched.tuff");
    let mut project = Project::create(&path).unwrap();
    let report = search_files(&mut project, &fasta, &files, SearchSettings::default()).unwrap();
    assert_eq!(report.failures(), 0, "{report:?}");
    // The FASTA's decoys are replaced by Sage's own.
    assert!(report.dropped_decoys > 0);
    for run in &report.runs {
        let searched = run.outcome.as_ref().unwrap();
        assert!(
            searched.searched > 1000,
            "{}: {}",
            run.run,
            searched.searched
        );
        let imported = searched.imported.as_ref().unwrap();
        assert!(imported.import.psms > 0 && imported.attach.attached > 0);
    }

    let summary = validate_project(&mut project, &Settings::default()).unwrap();
    assert_eq!(summary.engines.len(), 1);
    assert_eq!(summary.engines[0].name, "Sage");
    // DEMO0000..0149 are expressed; a passing group led by anything else is
    // a false discovery (an isoform sharing a parent's peptides can be).
    let list = project.protein_list().unwrap();
    let expressed = |accession: &str| {
        accession
            .get(7..11)
            .and_then(|digits| digits.parse::<usize>().ok())
            .is_some_and(|index| index < 150)
    };
    let wrong = list
        .rows
        .iter()
        .filter(|r| !expressed(&r.accession))
        .count();
    assert!(list.rows.len() >= 80, "{} groups", list.rows.len());
    assert!(
        wrong * 100 <= list.rows.len() * 3,
        "{wrong} of {}",
        list.rows.len()
    );

    // Searching the same file again is refused, not duplicated.
    let again = search_files(&mut project, &fasta, &files[..1], SearchSettings::default()).unwrap();
    assert_eq!(again.failures(), 1);
}
