//! Proteome Discoverer import. Tiny PD 2/3- and 1.x-shaped databases are
//! written by hand here (only the tables and columns the importer reads);
//! the expected values are derived by hand below. Real files run only when
//! `TUFFOLD_PD_SAMPLES` names a folder (`cargo test -p tuffold-io --test pd
//! -- --ignored --nocapture`).

use std::path::{Path, PathBuf};

use rusqlite::Connection;
use tuffold_core::{ModSite, SpectrumKey};
use tuffold_io::{PdError, import_pd};

const SCHEMA: &str = "
CREATE TABLE SchemaInfo (Version INT, Kind TEXT, Date TEXT, SoftwareVersion TEXT, Comment TEXT);
CREATE TABLE WorkflowInputFiles (WorkflowID INT, FileID INT, FileName TEXT);
CREATE TABLE ProcessingNodeScores (ScoreID INT, ScoreName TEXT, DisplayName TEXT, IsMainScore INT);
CREATE TABLE FoundModifications (ModificationID INT, Name TEXT, DeltaMonoisotopicMass REAL,
    PositionType INT, UnimodAccession INT);
CREATE TABLE TargetPsms (WorkflowID INT, PeptideID INT, IdentifyingNodeName TEXT,
    IdentifyingNodeTypeName TEXT, Sequence TEXT, Charge INT, MassOverCharge REAL,
    TheoreticalMass REAL, RetentionTime REAL, FirstScan INT, LastScan INT, SpectrumFileId INT,
    SearchEngineRank INT, XCorr REAL, DeltaCn REAL, MatchConfidence INT);
CREATE TABLE DecoyPsms (WorkflowID INT, PeptideID INT, IdentifyingNodeName TEXT,
    IdentifyingNodeTypeName TEXT, Sequence TEXT, Charge INT, MassOverCharge REAL,
    TheoreticalMass REAL, RetentionTime REAL, FirstScan INT, LastScan INT, SpectrumFileId INT,
    SearchEngineRank INT, XCorr REAL, DeltaCn REAL);
CREATE TABLE TargetPsmsFoundModifications (TargetPsmsWorkflowID INT, TargetPsmsPeptideID INT,
    FoundModificationsModificationID INT, Position INT);
CREATE TABLE DecoyPsmsFoundModifications (DecoyPsmsWorkflowID INT, DecoyPsmsPeptideID INT,
    FoundModificationsModificationID INT, Position INT);
CREATE TABLE TargetProteins (UniqueSequenceID INT, Accession TEXT, Description TEXT,
    Sequence TEXT, FastaTitleLines TEXT);
CREATE TABLE DecoyProteins (UniqueSequenceID INT, Accession TEXT, Description TEXT,
    Sequence TEXT, FastaTitleLines TEXT);
CREATE TABLE TargetProteinsTargetPsms (TargetProteinsUniqueSequenceID INT,
    TargetPsmsWorkflowID INT, TargetPsmsPeptideID INT);
CREATE TABLE DecoyProteinsDecoyPsms (DecoyProteinsUniqueSequenceID INT,
    DecoyPsmsWorkflowID INT, DecoyPsmsPeptideID INT);
";

// Theoretical MH+ by hand (monoisotopic residues, + water 18.010565,
// + proton 1.007276):
// PEPTIDE: 2×97.05276 + 2×129.04259 + 101.04768 + 113.08406 + 115.02694
//   = 781.34938; + water + proton = 800.367221; + phospho 79.966331
//   = 880.333552.
// AGTSK: 71.03711 + 57.02146 + 101.04768 + 87.03203 + 128.09496
//   = 444.23324; + water + proton = 463.251081; + acetyl 42.010565
//   = 505.261646.
// EDITPEP (the decoy) has PEPTIDE's composition: 800.367221.
const ROWS: &str = r"
INSERT INTO SchemaInfo VALUES (113, 'Result', '', '3.1.0.638', 'Created');
INSERT INTO WorkflowInputFiles VALUES (1, 9054, 'D:\raw\QS-1.raw');
INSERT INTO ProcessingNodeScores VALUES (653, 'XCorr', 'XCorr', 1);
INSERT INTO FoundModifications VALUES (1, 'Phospho', 79.966331, 0, 21),
    (2, 'Acetyl', 42.010565, 3, 1);
INSERT INTO TargetPsms VALUES
    (1, 10, 'Sequest HT (6)', 'Sequest HT', 'PEPTIDE', 2, 440.67, 880.333552, 0.5, 6, 6,
        9054, 1, 2.5, 0.0, 3),
    (1, 11, 'Sequest HT (6)', 'Sequest HT', 'AGTSK', 1, 505.26, 505.261646, 1.0, 9, 11,
        9054, 2, 1.25, 0.1, 1);
INSERT INTO DecoyPsms VALUES
    (1, 20, 'Sequest HT (6)', 'Sequest HT', 'EDITPEP', 2, 400.69, 800.367221, 0.5, 6, 6,
        9054, 1, 0.75, 0.0);
INSERT INTO TargetPsmsFoundModifications VALUES (1, 10, 1, 4), (1, 11, 2, 1);
INSERT INTO TargetProteins VALUES
    (100, 'protein', '1 sequence', 'MPEPTIDEK', '>protein 1 sequence'),
    (200, 'protein', '2 sequence', 'AGTSKR', '>protein 2 sequence');
INSERT INTO DecoyProteins VALUES (300, 'protein', '1 sequence', 'KEDITPEPM', '>protein 1 sequence');
INSERT INTO TargetProteinsTargetPsms VALUES (100, 1, 10), (200, 1, 11);
INSERT INTO DecoyProteinsDecoyPsms VALUES (300, 1, 20);
";

fn database(dir: &Path, name: &str, extra: &str) -> PathBuf {
    let path = dir.join(name);
    let conn = Connection::open(&path).unwrap();
    conn.execute_batch(SCHEMA).unwrap();
    conn.execute_batch(ROWS).unwrap();
    conn.execute_batch(extra).unwrap();
    path
}

#[test]
fn psms_modifications_proteins_and_decoys_import() {
    let dir = tempfile::tempdir().unwrap();
    let imported = import_pd(&database(dir.path(), "study.pdResult", "")).unwrap();

    assert_eq!(imported.batches.len(), 1);
    let batch = &imported.batches[0];
    // One node of its type: named by the type.
    assert_eq!(batch.software.name, "Sequest HT");
    assert_eq!(
        batch.software.version.as_deref(),
        Some("Proteome Discoverer 3.1.0.638")
    );
    assert_eq!(batch.ms_runs.len(), 1);
    assert_eq!(batch.ms_runs[0].name, "QS-1");
    assert_eq!(batch.psms.len(), 3);

    let phospho = &batch.psms[0];
    assert_eq!(phospho.peptide.proforma(), "PEPT[+79.9663]IDE");
    assert_eq!(phospho.peptide.modifications[0].site, ModSite::Residue(4));
    assert_eq!(
        phospho.peptide.modifications[0].accession.as_deref(),
        Some("UNIMOD:21")
    );
    assert_eq!(phospho.spectrum, SpectrumKey::Scan(6));
    assert_eq!((phospho.rank, phospho.charge), (1, 2));
    // (880.333552 − proton) / 2 + proton
    assert!((phospho.calculated_mz.unwrap() - 440.670414).abs() < 1e-6);
    let names: Vec<&str> = phospho.scores.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, ["XCorr", "ΔCn"]); // PD's Confidence is not imported
    // Both targets share the accession "protein": their title lines name them.
    assert_eq!(phospho.evidence[0].accession, "protein 1 sequence");

    let acetyl = &batch.psms[1];
    assert_eq!(acetyl.peptide.proforma(), "[+42.0106]-AGTSK");
    assert_eq!(acetyl.peptide.modifications[0].site, ModSite::NTerm);
    assert_eq!(acetyl.spectrum, SpectrumKey::Native("scans 9-11".into()));
    assert_eq!(acetyl.rank, 2);
    assert_eq!(acetyl.evidence[0].accession, "protein 2 sequence");

    let decoy = &batch.psms[2];
    assert_eq!(decoy.peptide.sequence, "EDITPEP");
    // Scan 6's decoy is rank 1 of the decoy search (its SearchEngineRank)
    // though the target outscores it (XCorr 2.5 vs 0.75): the two searches
    // keep their own ranks, and validation competes them per spectrum.
    assert_eq!(
        (decoy.spectrum.clone(), decoy.rank),
        (SpectrumKey::Scan(6), 1)
    );
    assert_eq!(decoy.evidence[0].accession, "DECOY_protein");
    let decoys = batch.decoy_index();
    assert!(decoys.is_decoy(decoy) && !decoys.is_decoy(phospho));
    let decoy_protein = batch
        .proteins
        .iter()
        .find(|p| p.accession == "DECOY_protein")
        .unwrap();
    assert!(decoy_protein.is_decoy && decoy_protein.decoy_explicit);

    // Spectrum 6 has a target and a decoy match: one spectrum, RT in seconds.
    assert_eq!(batch.spectra.len(), 2);
    assert_eq!(batch.spectra[0].rt_seconds, Some(30.0));
}

#[test]
fn a_misread_modification_refuses_the_file() {
    let dir = tempfile::tempdir().unwrap();
    // Phospho moved onto the decoy without changing its theoretical mass.
    let path = database(
        dir.path(),
        "bad.pdResult",
        "INSERT INTO DecoyPsmsFoundModifications VALUES (1, 20, 1, 4);",
    );
    let error = import_pd(&path).unwrap_err();
    assert!(matches!(error, PdError::Layout(_)), "{error}");
    assert!(error.to_string().contains("theoretical"), "{error}");
}

/// A PD 1.x `.msf` (only what the importer reads). Protein 10 is
/// `MKPEPTIDEKSAMPLERGGG`. Monoisotopic residues: S 87.03203, A 71.03711,
/// M 131.04049, P 97.05276, L 113.08406, E 129.04259, R 156.10111,
/// D 115.02694, I 113.08406, T 101.04768, K 128.09496, W 186.07931;
/// water 18.01056, proton 1.00728.
/// - Spectrum 1, target `SAMPLER`, acetylated N-terminus (terminal table)
///   and oxidized M at 0-based position 2: residues 784.39015 + water +
///   42.01057 + 15.99491 + proton = 861.41347 MH+. Linked to protein 10.
/// - Spectrum 2, decoy `EDITPEPK` (`PEPTIDE` reversed, K kept): residues
///   909.44434 → 928.46218 MH+. Not linked; `PEPTIDE` reversed but its K
///   is in protein 10, so it is protein 10's decoy.
/// - Spectrum 3, decoy `WWWWK`: 4 × 186.07931 + 128.09496 = 872.4122 →
///   891.43004 MH+. Not linked, and in no stored protein: a placeholder.
const PD1_SCHEMA: &str = "
CREATE TABLE SchemaInfo (Version INT, Kind TEXT, Date TEXT, SoftwareVersion TEXT, Comment TEXT);
CREATE TABLE FileInfos (FileID INT, FileName TEXT);
CREATE TABLE SpectrumHeaders (SpectrumID INT, MassPeakID INT, FirstScan INT, LastScan INT,
    RetentionTime REAL, Charge INT, Mass REAL);
CREATE TABLE MassPeaks (MassPeakID INT, Mass REAL, FileID INT);
CREATE TABLE ProcessingNodes (ProcessingNodeNumber INT, FriendlyName TEXT);
CREATE TABLE ProcessingNodeScores (ScoreID INT, ScoreName TEXT, FriendlyName TEXT);
CREATE TABLE AminoAcidModifications (AminoAcidModificationID INT, ModificationName TEXT,
    DeltaMass REAL, PositionType INT, UnimodAccession INT);
CREATE TABLE Proteins (ProteinID INT, Sequence TEXT);
CREATE TABLE ProteinAnnotations (ProteinID INT, Description TEXT);
CREATE TABLE Peptides (ProcessingNodeNumber INT, PeptideID INT, SpectrumID INT,
    SearchEngineRank INT, Sequence TEXT);
CREATE TABLE Peptides_decoy (ProcessingNodeNumber INT, PeptideID INT, SpectrumID INT,
    SearchEngineRank INT, Sequence TEXT);
CREATE TABLE PeptideScores (ProcessingNodeNumber INT, PeptideID INT, ScoreID INT,
    ScoreValue REAL);
CREATE TABLE PeptideScores_decoy (ProcessingNodeNumber INT, PeptideID INT, ScoreID INT,
    ScoreValue REAL);
CREATE TABLE PeptidesProteins (ProcessingNodeNumber INT, PeptideID INT, ProteinID INT);
CREATE TABLE PeptidesProteins_decoy (ProcessingNodeNumber INT, PeptideID INT, ProteinID INT);
CREATE TABLE PeptidesAminoAcidModifications (ProcessingNodeNumber INT, PeptideID INT,
    AminoAcidModificationID INT, Position INT);
CREATE TABLE PeptidesAminoAcidModifications_decoy (ProcessingNodeNumber INT, PeptideID INT,
    AminoAcidModificationID INT, Position INT);
CREATE TABLE PeptidesTerminalModifications (ProcessingNodeNumber INT, PeptideID INT,
    TerminalModificationID INT);
CREATE TABLE PeptidesTerminalModifications_decoy (ProcessingNodeNumber INT, PeptideID INT,
    TerminalModificationID INT);

INSERT INTO SchemaInfo VALUES (79, 'Result', '2013-11-06', '1.4.0.288', 'Created');
INSERT INTO FileInfos VALUES (7, 'Z:\\131026\\secretome_1.raw');
INSERT INTO SpectrumHeaders VALUES
    (1, 1, 100, 100, 10.0, 2, 861.41347),
    (2, 2, 200, 200, 20.0, 2, 928.46218),
    (3, 3, 300, 300, 30.0, 2, 891.43004);
INSERT INTO MassPeaks VALUES (1, 431.21, 7), (2, 464.73, 7), (3, 446.22, 7);
INSERT INTO ProcessingNodes VALUES (2, 'Spectrum Files'), (4, 'Sequest HT');
INSERT INTO ProcessingNodeScores VALUES (4, 'XCorr', 'XCorr');
INSERT INTO AminoAcidModifications VALUES
    (3, 'Acetyl', 42.010565, 1, 1), (45, 'Oxidation', 15.994915, 0, 35);
INSERT INTO Proteins VALUES (10, 'MKPEPTIDEKSAMPLERGGG');
INSERT INTO ProteinAnnotations VALUES (10, '>sp|P00001|TEST_HUMAN Test protein OS=Homo sapiens');
INSERT INTO Peptides VALUES (4, 1, 1, 1, 'SAMPLER');
INSERT INTO PeptideScores VALUES (4, 1, 4, 3.5);
INSERT INTO PeptidesProteins VALUES (4, 1, 10);
INSERT INTO PeptidesAminoAcidModifications VALUES (4, 1, 45, 2);
INSERT INTO PeptidesTerminalModifications VALUES (4, 1, 3);
INSERT INTO Peptides_decoy VALUES (4, 1, 2, 1, 'EDITPEPK'), (4, 2, 3, 1, 'WWWWK');
INSERT INTO PeptideScores_decoy VALUES (4, 1, 4, 1.2), (4, 2, 4, 0.8);
";

fn pd1_database(dir: &Path, extra: &str) -> PathBuf {
    let path = dir.join("old.msf");
    let conn = Connection::open(&path).unwrap();
    conn.execute_batch(PD1_SCHEMA).unwrap();
    conn.execute_batch(extra).unwrap();
    path
}

#[test]
fn pd1_msf_files_import_with_every_decoy_counted() {
    let dir = tempfile::tempdir().unwrap();
    let imported = import_pd(&pd1_database(dir.path(), "")).unwrap();
    assert_eq!(imported.batches.len(), 1);
    let batch = &imported.batches[0];
    assert_eq!(batch.software.name, "Sequest HT");
    assert_eq!(
        batch.software.version.as_deref(),
        Some("Proteome Discoverer 1.4.0.288")
    );
    assert_eq!(batch.ms_runs[0].name, "secretome_1");
    assert_eq!(batch.psms.len(), 3);

    let target = &batch.psms[0];
    assert_eq!(target.spectrum, SpectrumKey::Scan(100));
    assert_eq!(target.peptide.sequence, "SAMPLER");
    let mut sites: Vec<ModSite> = target
        .peptide
        .modifications
        .iter()
        .map(|m| m.site)
        .collect();
    sites.sort();
    assert_eq!(sites, [ModSite::NTerm, ModSite::Residue(3)]);
    assert_eq!(target.scores[0].name, "XCorr");
    assert_eq!(target.scores[0].value, 3.5);
    assert_eq!(target.evidence[0].accession, "sp|P00001|TEST_HUMAN");
    // (861.41347 − 1.00728) / 2 + 1.00728
    assert!((target.calculated_mz.unwrap() - 431.21038).abs() < 1e-4);

    let accessions: Vec<&str> = batch.psms[1..]
        .iter()
        .map(|psm| psm.evidence[0].accession.as_str())
        .collect();
    assert_eq!(
        accessions,
        ["DECOY_sp|P00001|TEST_HUMAN", "DECOY_unstored:WWWWK"]
    );
    assert!(
        batch
            .proteins
            .iter()
            .all(|p| p.is_decoy == p.accession.starts_with("DECOY_"))
    );
    let protein = batch
        .proteins
        .iter()
        .find(|p| p.accession == "sp|P00001|TEST_HUMAN")
        .unwrap();
    assert_eq!(
        protein.description.as_deref(),
        Some("Test protein OS=Homo sapiens")
    );
    assert_eq!(protein.sequence.as_deref(), Some("MKPEPTIDEKSAMPLERGGG"));
    let notes = imported.notes.join(" ");
    assert!(
        notes.contains("decoys 1) that Proteome Discoverer linked"),
        "{notes}"
    );
    assert!(notes.contains("1 decoy PSMs come from proteins"), "{notes}");
}

#[test]
fn a_pd1_psm_missing_its_precursor_mass_refuses_the_file() {
    let dir = tempfile::tempdir().unwrap();
    // The oxidation dropped: SAMPLER weighs 15.99491 less than spectrum 1,
    // and one PSM in three misses (over 5%).
    let path = pd1_database(dir.path(), "DELETE FROM PeptidesAminoAcidModifications;");
    let error = import_pd(&path).unwrap_err();
    assert!(matches!(error, PdError::Layout(_)), "{error}");
    assert!(error.to_string().contains("1 of 3 PSMs miss"), "{error}");
}

#[test]
fn other_databases_are_named_not_misread() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.msf");
    Connection::open(&other)
        .unwrap()
        .execute_batch("CREATE TABLE Foo (x INT);")
        .unwrap();
    assert!(matches!(import_pd(&other), Err(PdError::NotPd(_))));
}

#[test]
#[ignore = "needs Proteome Discoverer sample files (TUFFOLD_PD_SAMPLES)"]
fn sample_files_import() {
    let Ok(folder) = std::env::var("TUFFOLD_PD_SAMPLES") else {
        panic!("set TUFFOLD_PD_SAMPLES");
    };
    let mut files: Vec<PathBuf> = walk(Path::new(&folder));
    files.sort();
    assert!(!files.is_empty(), "no .msf/.pdResult files");
    for file in files {
        let started = std::time::Instant::now();
        let imported =
            import_pd(&file).unwrap_or_else(|error| panic!("{}: {error}", file.display()));
        println!(
            "== {} ({:.1} s)",
            file.display(),
            started.elapsed().as_secs_f64()
        );
        for note in &imported.notes {
            println!("   note: {note}");
        }
        for batch in &imported.batches {
            let decoys = batch.decoy_index();
            let rank_one: Vec<_> = batch.psms.iter().filter(|psm| psm.rank == 1).collect();
            let decoy = rank_one.iter().filter(|psm| decoys.is_decoy(psm)).count();
            println!(
                "   {}: {} runs, {} spectra, {} PSMs ({} rank 1, {decoy} decoy), {} proteins",
                batch.software.name,
                batch.ms_runs.len(),
                batch.spectra.len(),
                batch.psms.len(),
                rank_one.len(),
                batch.proteins.len()
            );
            assert!(batch.psms.iter().all(|psm| !psm.scores.is_empty()));
        }
    }
}

fn walk(folder: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for entry in std::fs::read_dir(folder)
        .expect("readable folder")
        .flatten()
    {
        let path = entry.path();
        if path.is_dir() {
            out.extend(walk(&path));
        } else if path.extension().is_some_and(|ext| {
            ext.eq_ignore_ascii_case("msf") || ext.eq_ignore_ascii_case("pdresult")
        }) {
            out.push(path);
        }
    }
    out
}
