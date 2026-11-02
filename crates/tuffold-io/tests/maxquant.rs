//! MaxQuant `txt` import from a hand-written folder; expected values are
//! derived by hand below. Real output runs only when
//! `TUFFOLD_MAXQUANT_SAMPLES` names a folder of `txt` folders.

use std::path::Path;

use tuffold_core::{ModSite, SourceFormat, SpectrumKey};
use tuffold_io::{MaxQuantError, detect_format, import_maxquant};

const HEADER: &str = "Raw file\tScan number\tSequence\tModified sequence\tProteins\tCharge\tm/z\t\
                      Mass\tRetention time\tScore\tDelta score\tPEP\tReverse";

// Masses by hand (monoisotopic residues + water 18.010565):
// MPEPTIDE = 131.04049 + 781.34938 (PEPTIDE) = 912.38987; + water
//   = 930.400435; + acetyl 42.010565 + oxidation 15.994915 = 988.405915.
// ACDCK = 71.03711 + 2×103.00919 + 115.02694 + 128.09496 = 520.17739;
//   + water = 538.187955; + 2 × carbamidomethyl 57.021464 = 652.230883.
// EDITPEPK = 781.34938 + 128.09496 = 909.44434; + water = 927.454905.
const ROWS: &[&str] = &[
    "runA\t100\tMPEPTIDE\t_(Acetyl (Protein N-term))M(Oxidation (M))PEPTIDE_\tP1;CON__P2\t2\t\
     495.21\t988.405915\t10.5\t120.5\t80.25\t0.0001\t",
    "runA\t100\tACDCK\t_ACDCK_\tP3\t2\t327.12\t652.230883\t10.5\t40\t5\t0.2\t",
    "runA\t200\tEDITPEPK\t_EDITPEPK_\tREV__P4\t2\t464.73\t927.454905\t12\t35\t1\t0.5\t+",
    "runA\t300\tPEPTIDE\t_PEPT(Sulfo (T))IDE_\tP1\t2\t440.67\t879.3\t13\t50\t2\t0.1\t",
];

fn txt_folder(dir: &Path, rows: &[&str]) -> std::path::PathBuf {
    let folder = dir.join("txt");
    std::fs::create_dir(&folder).unwrap();
    let mut msms = format!("{HEADER}\n");
    for row in rows {
        msms.push_str(row);
        msms.push('\n');
    }
    std::fs::write(folder.join("msms.txt"), msms).unwrap();
    std::fs::write(
        folder.join("parameters.txt"),
        "Parameter\tValue\nVersion\t2.4.2.0\nFixed modifications\tCarbamidomethyl (C)\n",
    )
    .unwrap();
    std::fs::write(
        folder.join("summary.txt"),
        "Raw file\tExperiment\nrunA\tcontrol\nTotal\t\n",
    )
    .unwrap();
    folder
}

#[test]
fn psms_ranks_modifications_and_experiments_import() {
    let dir = tempfile::tempdir().unwrap();
    let folder = txt_folder(dir.path(), ROWS);
    assert_eq!(detect_format(&folder), Some(SourceFormat::MaxQuant));
    assert_eq!(
        detect_format(&folder.join("msms.txt")),
        Some(SourceFormat::MaxQuant)
    );

    let imported = import_maxquant(&folder).unwrap();
    let batch = &imported.batch;
    assert_eq!(batch.software.version.as_deref(), Some("2.4.2.0"));
    assert_eq!(batch.biosamples.len(), 1);
    assert_eq!(batch.biosamples[0].name, "control");
    assert_eq!(batch.ms_runs[0].biosample, "control");
    // Scan 300's unknown modification leaves it out.
    assert_eq!(batch.spectra.len(), 2);
    assert_eq!(batch.spectra[0].rt_seconds, Some(630.0));
    assert_eq!(batch.psms.len(), 3);

    let best = &batch.psms[0];
    assert_eq!(best.spectrum, SpectrumKey::Scan(100));
    assert_eq!(best.rank, 1);
    assert_eq!(best.peptide.proforma(), "[+42.0106]-M[+15.9949]PEPTIDE");
    let names: Vec<&str> = best.scores.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, ["Andromeda:score", "MaxQuant:delta score"]); // no PEP
    assert_eq!(best.evidence.len(), 2);

    let second = &batch.psms[1];
    assert_eq!(
        (second.rank, second.peptide.sequence.as_str()),
        (2, "ACDCK")
    );
    assert_eq!(second.peptide.modifications[0].site, ModSite::Residue(2));

    let reverse = &batch.psms[2];
    assert!(batch.decoy_index().is_decoy(reverse));
    assert!(!batch.decoy_index().is_decoy(best)); // CON__ is a target
    assert!(
        imported.notes.iter().any(|note| note.contains("Sulfo (T)")),
        "{:?}",
        imported.notes
    );
}

#[test]
fn masses_that_disagree_refuse_the_file() {
    let dir = tempfile::tempdir().unwrap();
    // SILAC-labelled masses: K +8.014199 that the sequence doesn't show.
    let rows = [
        "runA\t1\tEDITPEPK\t_EDITPEPK_\tP4\t2\t468.74\t935.469104\t12\t35\t1\t0.5\t",
        "runA\t2\tEDITPEPK\t_EDITPEPK_\tP4\t2\t468.74\t935.469104\t12\t36\t1\t0.5\t",
    ];
    let folder = txt_folder(dir.path(), &rows);
    let error = import_maxquant(&folder).unwrap_err();
    assert!(
        matches!(
            error,
            MaxQuantError::Masses {
                failures: 2,
                checked: 2,
                ..
            }
        ),
        "{error}"
    );
}

#[test]
fn a_folder_without_msms_is_not_maxquant() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        detect_format(dir.path()),
        Some(SourceFormat::FragPipeBundle)
    );
    assert!(matches!(
        import_maxquant(dir.path()),
        Err(MaxQuantError::NoMsms(_))
    ));
}

#[test]
#[ignore = "needs MaxQuant output (TUFFOLD_MAXQUANT_SAMPLES)"]
fn sample_folders_import() {
    let folder = std::env::var("TUFFOLD_MAXQUANT_SAMPLES").expect("set TUFFOLD_MAXQUANT_SAMPLES");
    for entry in std::fs::read_dir(folder).unwrap().flatten() {
        if detect_format(&entry.path()) != Some(SourceFormat::MaxQuant) {
            continue;
        }
        let imported = import_maxquant(&entry.path()).unwrap_or_else(|e| panic!("{e}"));
        let batch = &imported.batch;
        let decoys = batch.decoy_index();
        let rank_one = batch.psms.iter().filter(|psm| psm.rank == 1).count();
        let decoy = batch
            .psms
            .iter()
            .filter(|psm| psm.rank == 1 && decoys.is_decoy(psm))
            .count();
        println!(
            "== {}: {} runs, {} PSMs ({rank_one} rank 1, {decoy} decoy), {} proteins",
            entry.path().display(),
            batch.ms_runs.len(),
            batch.psms.len(),
            batch.proteins.len()
        );
        for note in &imported.notes {
            println!("   note: {note}");
        }
    }
}
