//! Importing a set of inputs into a project.
//!
//! Identification results are imported first (each atomically), then spectra
//! files attach their peaks, so spectra always find the runs they belong to.
//! Attaching a file also traces its identified precursors through its MS1
//! scans (MS1 intensity; see `tuffold_quant::ms1`) in the same pass.
//! One failing input never stops the others: every input gets an [`Outcome`].

use std::path::{Path, PathBuf};

use crate::undo::{Checkpoint, import_count};
use tuffold_core::centroid::centroid;
use tuffold_core::{PeakSpectrum, SourceFormat};
use tuffold_io::{
    InspectionIssue, SpectraError, StoredPeaks, THERMO_RAW, detect_format, find_search_files,
    import_path, inspect_bundle, is_spectra_file, is_thermo_raw, read_scans, spectra_run_name,
};
use tuffold_quant::ms1::{TraceSettings, Tracer};
use tuffold_quant::purity::{PURITY_PPM, PurityTracker};
use tuffold_store::{AttachReport, ImportReport, ImportSource, Project, StoreError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputKind {
    Identifications(SourceFormat),
    Spectra,
}

/// What kind of input a path is, from its name (and whether it is a folder).
pub fn classify(path: &Path) -> Option<InputKind> {
    if is_spectra_file(path) {
        Some(InputKind::Spectra)
    } else {
        detect_format(path).map(InputKind::Identifications)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum OutcomeKind {
    Imported {
        format: SourceFormat,
        software: String,
        report: ImportReport,
        /// Inspection warnings (FragPipe folders); the data still imported.
        warnings: Vec<InspectionIssue>,
        /// What the importer left out or wants known (Scaffold projects).
        notes: Vec<String>,
    },
    Attached(AttachReport),
    Failed(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Outcome {
    pub input: PathBuf,
    pub kind: OutcomeKind,
}

impl Outcome {
    pub fn is_failure(&self) -> bool {
        matches!(self.kind, OutcomeKind::Failed(_))
    }
}

/// A sentence describing repeated scores or evidence an import kept once,
/// for front ends to show; `None` when there were none.
pub fn duplicate_note(report: &ImportReport) -> Option<String> {
    let mut parts = Vec::new();
    if report.duplicate_scores > 0 {
        parts.push(format!(
            "{} repeated score name(s) on the same PSM (the first value was kept)",
            report.duplicate_scores
        ));
    }
    if report.duplicate_evidence > 0 {
        parts.push(format!(
            "{} repeated protein locus/loci on the same PSM (kept once)",
            report.duplicate_evidence
        ));
    }
    (!parts.is_empty()).then(|| format!("The source listed {}.", parts.join(" and ")))
}

/// A sentence on the MS3 scans (SPS-MS3) an attached file held, if any.
pub fn ms3_note(report: &AttachReport) -> Option<String> {
    (report.ms3_scans > 0).then(|| {
        let skipped = report.ms3_scans - report.ms3_attached;
        let mut note = format!(
            "{} of {} MS3 scans stored with the identified MS/MS spectrum they were acquired \
             from; TMT reporters are read there",
            report.ms3_attached, report.ms3_scans
        );
        if skipped > 0 {
            note.push_str(&format!(
                " ({skipped} name no identified spectrum, or one that already has its MS3)"
            ));
        }
        note
    })
}

/// A sentence on what an attached file contributed to MS1 intensity.
pub fn ms1_note(report: &AttachReport) -> String {
    if report.ms1_scans == 0 {
        let damaged = match report.ms1_scans_damaged {
            0 => String::new(),
            n => format!(" ({n} damaged MS1 scans left out)"),
        };
        return match report.ms1_scans_without_rt {
            0 => format!(
                "no MS1 scans read{damaged}, so no precursor intensities (those need an \
                 mzML with its MS1 scans, or a timsTOF .d)"
            ),
            n => format!(
                "{n} MS1 scans but none has a retention time{damaged}, so no precursor \
                 intensities"
            ),
        };
    }
    let mut extra = Vec::new();
    if report.ms1_scans_without_rt > 0 {
        extra.push(format!(
            "{} MS1 scans without a retention time skipped",
            report.ms1_scans_without_rt
        ));
    }
    if report.ms1_scans_damaged > 0 {
        extra.push(format!(
            "{} damaged MS1 scans (data that can't be decoded, or peaks no instrument \
             could measure) left out",
            report.ms1_scans_damaged
        ));
    }
    if report.precursors_unplaceable > 0 {
        extra.push(format!(
            "{} identifications without a known precursor m/z (unknown residues or no \
             charge) not traced",
            report.precursors_unplaceable
        ));
    }
    let extra: String = extra.iter().map(|e| format!("; {e}")).collect();
    let replaced = match report.precursors_replaced {
        0 => String::new(),
        n => format!(" (replacing the {n} traced from this run's earlier spectra file)"),
    };
    format!(
        "{} identified precursors traced through {} MS1 scans, {} with an elution \
         peak{replaced}{extra}",
        report.precursors_traced, report.ms1_scans, report.precursors_measured
    )
}

/// Import `inputs` into `project`: identifications first, then spectra.
///
/// `run` names the run of a spectra file whose file name differs from it;
/// it applies to every spectra input.
///
/// The project is copied aside first, and the copy kept as the undo point
/// when anything was imported (see [`crate::undo`]).
pub fn import_inputs(project: &mut Project, inputs: &[PathBuf], run: Option<&str>) -> Vec<Outcome> {
    let label = match inputs {
        [one] => format!(
            "import {}",
            one.file_name().map_or_else(
                || one.display().to_string(),
                |n| n.to_string_lossy().into_owned()
            )
        ),
        many => format!("import {} inputs", many.len()),
    };
    let checkpoint = match Checkpoint::take(project, &label) {
        Ok(checkpoint) => checkpoint,
        Err(error) => {
            // Without an undo copy, import nothing: the user expects undo.
            return inputs
                .iter()
                .map(|input| Outcome {
                    input: input.clone(),
                    kind: OutcomeKind::Failed(format!(
                        "not imported: an undo copy could not be saved first ({error})"
                    )),
                })
                .collect();
        }
    };
    let before = import_count(project);
    let mut outcomes = import_each(project, inputs, run);
    // Kept only for a real change: a new import, peaks (or MS3 peaks) newly
    // attached, or a run's first tracing. Re-attaching a file that was
    // already attached changes nothing, so it must not replace the real
    // change's undo copy.
    let imported = import_count(project) > before
        || outcomes.iter().any(|outcome| {
            matches!(&outcome.kind, OutcomeKind::Attached(report)
                if report.attached > 0
                    || report.ms3_attached > 0
                    || (report.precursors_traced > 0 && report.precursors_replaced == 0))
        });
    if let (true, Some(checkpoint)) = (imported, checkpoint)
        && let Err(error) = checkpoint.keep()
    {
        outcomes.push(Outcome {
            input: project.path().map(Path::to_path_buf).unwrap_or_default(),
            kind: OutcomeKind::Failed(format!(
                "imported, but the undo copy could not be kept, so this can't be undone \
                     ({error})"
            )),
        });
    }
    outcomes
}

fn import_each(project: &mut Project, inputs: &[PathBuf], run: Option<&str>) -> Vec<Outcome> {
    let mut identifications = Vec::new();
    let mut spectra = Vec::new();
    let mut outcomes = Vec::new();
    for input in inputs {
        // The canonical path identifies a file for duplicate-import checks, so
        // a symlink or `..` alias cannot import the same file twice.
        let absolute = std::fs::canonicalize(input)
            .or_else(|_| std::path::absolute(input))
            .unwrap_or_else(|_| input.clone());
        match classify(&absolute) {
            Some(InputKind::Spectra) => spectra.push(absolute),
            Some(InputKind::Identifications(_)) => identifications.push(absolute),
            None => {
                let why = if is_thermo_raw(&absolute) {
                    THERMO_RAW.into()
                } else {
                    "unrecognized input; expected a FragPipe or MaxQuant folder, a .pep.xml/.pepXML \
                     file, an .mzid, .dat, .pdResult/.msf, .sf3 or .sfdb file, or a .mzML/.mgf \
                     spectra file"
                        .into()
                };
                outcomes.push(Outcome {
                    input: absolute,
                    kind: OutcomeKind::Failed(why),
                });
            }
        }
    }
    for input in identifications {
        let kind = import_identifications(project, &input);
        outcomes.push(Outcome { input, kind });
    }
    for input in spectra {
        let run_name = run.map_or_else(|| spectra_run_name(&input), str::to_string);
        let kind = match attach_file(project, &input, &run_name) {
            Ok(report) => OutcomeKind::Attached(report),
            Err(message) => OutcomeKind::Failed(message),
        };
        outcomes.push(Outcome { input, kind });
    }
    outcomes
}

fn import_identifications(project: &mut Project, input: &Path) -> OutcomeKind {
    // A FragPipe folder holding the engine's pepXML is imported from those,
    // so the psm.tsv inspection does not apply.
    let reads_psm_table = detect_format(input) == Some(SourceFormat::FragPipeBundle)
        && find_search_files(input).is_ok_and(|files| files.is_empty());
    let warnings = if reads_psm_table {
        inspect_bundle(input).warnings
    } else {
        Vec::new()
    };
    let imported = match import_path(input) {
        Ok(imported) => imported,
        Err(error) => return OutcomeKind::Failed(error.to_string()),
    };
    if imported.batches.is_empty() {
        return OutcomeKind::Failed(format!("nothing to import. {}", imported.notes.join(" ")));
    }
    let source = ImportSource {
        path: input.display().to_string(),
        format: imported.format,
    };
    match project.import_all(&source, &imported.batches) {
        Ok(reports) => {
            let mut report = reports[0].clone();
            for other in &reports[1..] {
                report.absorb(other);
            }
            let mut notes = imported.notes;
            if let Some(peaks) = &imported.peaks {
                match attach_stored_peaks(project, peaks) {
                    Ok(note) => notes.push(note),
                    Err(error) => {
                        notes.push(format!("The project's peaks were not attached: {error}"))
                    }
                }
            }
            OutcomeKind::Imported {
                format: imported.format,
                software: imported
                    .batches
                    .iter()
                    .map(|batch| batch.software.name.as_str())
                    .collect::<Vec<_>>()
                    .join(" + "),
                report,
                warnings,
                notes,
            }
        }
        Err(error) => OutcomeKind::Failed(error.to_string()),
    }
}

/// Attach the peaks a Scaffold project stores to its identified spectra, run
/// by run, as a spectra file's would be (never replacing stored peaks). A
/// spectrum whose peaks can't be read is left without and counted.
fn attach_stored_peaks(project: &mut Project, peaks: &StoredPeaks) -> Result<String, String> {
    let (mut attached, mut already, mut unreadable) = (0, 0, 0);
    let mut example = None;
    for run in peaks.runs() {
        let spectra = peaks.spectra(run).filter_map(|spectrum| match spectrum {
            Ok(spectrum) => Some(Ok(spectrum)),
            Err(problem) => {
                unreadable += 1;
                example.get_or_insert(problem.to_string());
                None
            }
        });
        let report = attach_peaks(project, run, spectra, true)?;
        attached += report.attached;
        already += report.already_attached;
    }
    let mut note = format!(
        "Peaks of {attached} of the {} identified spectra were read from the project",
        peaks.len()
    );
    if already > 0 {
        note.push_str(&format!(" ({already} already had peaks, which were kept)"));
    }
    if attached > 0 {
        note.push_str(
            "; Scaffold keeps at most 300 per spectrum, so attaching the runs' spectra files \
             replaces them with complete ones",
        );
    }
    if unreadable > 0 {
        note.push_str(&format!(
            "; {unreadable} spectra's peaks couldn't be read and were left out (e.g. {})",
            example.unwrap_or_default()
        ));
    }
    note.push('.');
    Ok(note)
}

/// Attach a spectra file's MS/MS peaks to `run` and trace the run's
/// identified precursors through its MS1 scans, in one pass and one
/// transaction.
pub(crate) fn attach_file(
    project: &mut Project,
    input: &Path,
    run: &str,
) -> Result<AttachReport, String> {
    let scans = read_scans(input).map_err(|error| error.to_string())?;
    attach_read(project, run, scans)
}

/// [`attach_scans`] over a reader's scans: a damaged MS1 scan (a timsTOF
/// frame, or an mzML scan with impossible peaks) is left out and counted,
/// and its MS/MS spectra are attached as usual; any other error stops.
fn attach_read(
    project: &mut Project,
    run: &str,
    scans: impl Iterator<Item = Result<PeakSpectrum, SpectraError>>,
) -> Result<AttachReport, String> {
    let mut damaged = 0;
    let scans = scans.filter_map(|scan| match scan {
        Err(SpectraError::DamagedMs1 { .. }) => {
            damaged += 1;
            None
        }
        other => Some(other.map_err(|e| e.to_string())),
    });
    let mut report = attach_scans(project, run, scans)?;
    report.ms1_scans_damaged = damaged;
    Ok(report)
}

/// Attach scans (MS1 and MS/MS, in file order) to `run`. Profile MS1 scans
/// are centroided first; scans without a retention time can't be placed in
/// time and don't take part in tracing.
pub(crate) fn attach_scans(
    project: &mut Project,
    run: &str,
    scans: impl Iterator<Item = Result<PeakSpectrum, String>>,
) -> Result<AttachReport, String> {
    attach_peaks(project, run, scans, false)
}

/// [`attach_scans`], for complete peaks (a spectra file's) or `partial` ones
/// (a project file's, which a spectra file later replaces).
fn attach_peaks(
    project: &mut Project,
    run: &str,
    scans: impl Iterator<Item = Result<PeakSpectrum, String>>,
    partial: bool,
) -> Result<AttachReport, String> {
    let mut attacher = if partial {
        project.attach_partial_peaks(run)
    } else {
        project.attach_spectra(run)
    }
    .map_err(|error: StoreError| error.to_string())?;
    let mut tracer = Tracer::new(TraceSettings::default());
    let mut purity = PurityTracker::<String>::new(PURITY_PPM);
    let mut without_rt = 0;
    for scan in scans {
        let scan = scan?;
        if scan.ms_level == 1 {
            let Some(rt) = scan.rt_seconds else {
                without_rt += 1;
                continue;
            };
            let centroided;
            let (mz, intensity) = if scan.profile {
                centroided = centroid(&scan.mz, &scan.intensity);
                (&centroided.0, &centroided.1)
            } else {
                (&scan.mz, &scan.intensity)
            };
            tracer.add_ms1(rt, mz, intensity, &scan.mobility);
            // A frame binned by mobility holds an ion once per bin, and
            // purity reads one peak per isotope, so it isn't measured on
            // them (timsTOF MS/MS spectra carry no isolation window anyway).
            if scan.mobility.is_empty() {
                purity.add_ms1(rt, mz, intensity);
                for (key, value) in purity.take() {
                    attacher
                        .record_purity(&key, value)
                        .map_err(|error| error.to_string())?;
                }
            }
            continue;
        }
        if let Some(rt) = scan.rt_seconds {
            let precursors = attacher.precursors(&scan.key);
            for precursor in precursors {
                let key = (precursor.peptide_id, precursor.charge);
                tracer.add_identification(
                    key,
                    precursor.mz,
                    precursor.charge,
                    rt,
                    scan.precursor_mobility,
                );
            }
            // Purity of the MS/MS scan's own window around its (first)
            // identified precursor; an MS3 scan's precursor is its MS/MS's.
            if let (2, Some(window), Some(precursor)) =
                (scan.ms_level, scan.isolation, precursors.first())
            {
                purity.add_msms(
                    scan.key.canonical(),
                    rt,
                    precursor.mz,
                    precursor.charge,
                    window,
                );
            }
        }
        attacher.add(&scan).map_err(|error| error.to_string())?;
    }
    for (key, value) in purity.finish() {
        attacher
            .record_purity(&key, value)
            .map_err(|error| error.to_string())?;
    }
    let ms1_scans = tracer.ms1_scans();
    attacher
        .record_ms1(ms1_scans, &tracer.finish())
        .map_err(|error| error.to_string())?;
    let mut report = attacher.finish().map_err(|error| error.to_string())?;
    report.ms1_scans_without_rt = without_rt;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixtures() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures")
    }

    #[test]
    fn a_damaged_ms1_frame_is_left_out_and_counted() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        crate::create_demo_project(&path).unwrap();
        let mut project = Project::open(&path).unwrap();
        let damaged = || SpectraError::DamagedMs1 {
            file: "run.d".into(),
            scan: "frame=7".into(),
            problem: "Decompression fails".into(),
        };
        let ms1 = PeakSpectrum {
            key: tuffold_core::SpectrumKey::Native("frame=8".into()),
            native_id: "frame=8".into(),
            ms_level: 1,
            rt_seconds: Some(1.0),
            precursor_mz: None,
            precursor_charge: None,
            mz: vec![500.0],
            intensity: vec![10.0],
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        };

        let report = attach_read(
            &mut project,
            "control-A",
            [Err(damaged()), Ok(ms1)].into_iter(),
        )
        .unwrap();

        assert_eq!((report.ms1_scans, report.ms1_scans_damaged), (1, 1));
        assert!(ms1_note(&report).contains("1 damaged MS1 scans"));
        // Any other unreadable spectrum still stops the attach.
        let invalid = SpectraError::Invalid {
            file: "run.d".into(),
            id: "spectrum 3".into(),
            problem: "bad".into(),
        };
        assert!(attach_read(&mut project, "control-A", [Err(invalid)].into_iter()).is_err());
    }

    /// Three scans in run control-A's mzML, windows 300-1600: scan 1 is MS1
    /// (m/z 400/500/600, intensities 10/20/30, RT 1 min), scan 2 is MS1
    /// with an m/z of 1.5788030582872526e-19 first (RT 2 min), scan 3 is
    /// MS/MS (m/z 200/300, intensities 5/6, RT 3 min) of scan 2. The
    /// arrays are base64 of little-endian doubles (Python's struct).
    const DAMAGED_PRECURSOR_SCAN: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0" id="control-A">
<cvList count="2"><cv id="MS" fullName="MS" version="4.1" URI="x"/><cv id="UO" fullName="UO" version="1" URI="x"/></cvList>
<run id="control-A"><spectrumList count="3">
<spectrum index="0" id="scan=1" defaultArrayLength="3">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/><scan>
<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
<scanWindowList count="1"><scanWindow><cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="300"/><cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1600"/></scanWindow></scanWindowList>
</scan></scanList>
<binaryDataArrayList count="2">
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/><binary>AAAAAAAAeUAAAAAAAEB/QAAAAAAAwIJA</binary></binaryDataArray>
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/><binary>AAAAAAAAJEAAAAAAAAA0QAAAAAAAAD5A</binary></binaryDataArray>
</binaryDataArrayList></spectrum>
<spectrum index="1" id="scan=2" defaultArrayLength="3">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/><scan>
<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="2" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
<scanWindowList count="1"><scanWindow><cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="300"/><cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1600"/></scanWindow></scanWindowList>
</scan></scanList>
<binaryDataArrayList count="2">
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/><binary>AAAAoIxMBzwAAAAAAEB/QAAAAAAAwIJA</binary></binaryDataArray>
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/><binary>AAAAAAAAJEAAAAAAAAA0QAAAAAAAAD5A</binary></binaryDataArray>
</binaryDataArrayList></spectrum>
<spectrum index="2" id="scan=3" defaultArrayLength="2">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/><scan>
<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="3" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
</scan></scanList>
<precursorList count="1"><precursor spectrumRef="scan=2"><selectedIonList count="1"><selectedIon>
<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="500"/>
<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/>
</selectedIon></selectedIonList></precursor></precursorList>
<binaryDataArrayList count="2">
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/><binary>AAAAAAAAaUAAAAAAAMByQA==</binary></binaryDataArray>
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/><cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/><cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/><binary>AAAAAAAAFEAAAAAAAAAYQA==</binary></binaryDataArray>
</binaryDataArrayList></spectrum>
</spectrumList></run></mzML>"#;

    #[test]
    fn an_mzml_ms1_scan_with_garbage_peaks_is_left_out_counted_and_its_msms_still_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        crate::create_demo_project(&path).unwrap();
        let mzml = dir.path().join("control-A.mzML");
        std::fs::write(&mzml, DAMAGED_PRECURSOR_SCAN).unwrap();
        let mut project = Project::open(&path).unwrap();

        let report = attach_file(&mut project, &mzml, "control-A").unwrap();

        // Scan 1 is traced; scan 2 is left out and counted; scan 3, the
        // MS/MS of scan 2, is read (it names no demo spectrum, so isn't stored).
        assert_eq!(
            (report.ms1_scans, report.ms1_scans_damaged, report.read),
            (1, 1, 1)
        );
        let note = ms1_note(&report);
        assert!(note.contains("through 1 MS1 scans"), "{note}");
        assert!(note.contains("1 damaged MS1 scans"), "{note}");
    }

    #[test]
    fn precursor_purity_is_measured_around_each_identified_spectrum() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        crate::create_demo_project(&path).unwrap();
        // An identified spectrum of control-A and its precursor.
        let keys: Vec<String> = rusqlite::Connection::open(&path)
            .unwrap()
            .prepare(
                "SELECT s.spectrum_key FROM spectra s JOIN ms_runs r ON r.id = s.ms_run_id
                 WHERE r.name = 'control-A' ORDER BY s.id",
            )
            .unwrap()
            .query_map([], |row| row.get(0))
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        let mut project = Project::open(&path).unwrap();
        let (key, precursor) = {
            let attacher = project.attach_spectra("control-A").unwrap();
            keys.iter()
                .map(|k| tuffold_core::SpectrumKey::from_native_id(k))
                .find_map(|k| attacher.precursors(&k).first().map(|p| (k, *p)))
                .unwrap()
        };
        let (mono, z) = (precursor.mz, precursor.charge);
        let m1 = mono + tuffold_quant::ms1::ISOTOPE_SPACING / f64::from(z);
        // Before (t = 10): M 100, M+1 50, a contaminant 0.2 above M at 50:
        // 150 / 200 = 0.75. After (t = 12): contaminant 150: 150 / 300 = 0.5.
        // The MS/MS scan at t = 11 is halfway: 0.625.
        let ms1 = |rt: f64, contaminant: f32| PeakSpectrum {
            key: tuffold_core::SpectrumKey::Native(format!("ms1-{rt}")),
            native_id: format!("ms1-{rt}"),
            ms_level: 1,
            rt_seconds: Some(rt),
            precursor_mz: None,
            precursor_charge: None,
            mz: vec![mono, mono + 0.2, m1],
            intensity: vec![100.0, contaminant, 50.0],
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        };
        let msms = PeakSpectrum {
            key: key.clone(),
            native_id: key.canonical(),
            ms_level: 2,
            rt_seconds: Some(11.0),
            precursor_mz: Some(mono),
            precursor_charge: Some(i16::from(z)),
            mz: vec![200.0],
            intensity: vec![1.0],
            profile: false,
            parent: None,
            isolation: Some((mono - 1.0, mono + 1.0)),
            precursor_mobility: None,
            mobility: Vec::new(),
        };
        let report = attach_scans(
            &mut project,
            "control-A",
            [ms1(10.0, 50.0), msms, ms1(12.0, 150.0)]
                .into_iter()
                .map(Ok),
        )
        .unwrap();
        assert_eq!(report.purity_measured, 1);
        drop(project);
        let purity: f64 = rusqlite::Connection::open(&path)
            .unwrap()
            .query_row(
                "SELECT precursor_purity FROM spectra WHERE spectrum_key = ?1",
                [key.canonical()],
                |row| row.get(0),
            )
            .unwrap();
        assert!((purity - 0.625).abs() < 1e-9, "{purity}");
    }

    #[test]
    fn inputs_are_classified_by_name() {
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(
            classify(dir.path()),
            Some(InputKind::Identifications(SourceFormat::FragPipeBundle))
        );
        assert_eq!(classify(Path::new("a.mzML")), Some(InputKind::Spectra));
        assert_eq!(
            classify(Path::new("a.mzid")),
            Some(InputKind::Identifications(SourceFormat::MzIdentMl))
        );
        assert_eq!(classify(Path::new("a.txt")), None);
    }

    #[test]
    fn one_failure_does_not_stop_the_others() {
        let dir = tempfile::tempdir().unwrap();
        let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
        let notes = dir.path().join("notes.txt");
        std::fs::write(&notes, "hi").unwrap();
        let inputs = vec![
            notes,
            fixtures().join("third_party/pyteomics/test.mgf"), // run "test" is unknown
            fixtures().join("third_party/pyteomics/test.mzid"),
        ];

        let outcomes = import_inputs(&mut project, &inputs, None);

        // Order: unrecognized inputs, then identifications, then spectra.
        let failed: Vec<bool> = outcomes.iter().map(Outcome::is_failure).collect();
        assert_eq!(failed, [true, false, true]);
        let OutcomeKind::Imported {
            software, report, ..
        } = &outcomes[1].kind
        else {
            panic!("{:?}", outcomes[1]);
        };
        assert_eq!((software.as_str(), report.psms), ("Sequest + Mascot", 22));
        assert_eq!(project.summary().unwrap().psms, 22);
    }

    #[test]
    fn fragpipe_warnings_travel_with_the_outcome() {
        let dir = tempfile::tempdir().unwrap();
        let bundle = dir.path().join("bundle");
        copy_dir(&fixtures().join("fragpipe_synthetic"), &bundle);
        std::fs::remove_file(bundle.join("protein.fas")).unwrap();
        let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();

        let outcomes = import_inputs(&mut project, &[bundle], None);

        let OutcomeKind::Imported { warnings, .. } = &outcomes[0].kind else {
            panic!("{:?}", outcomes[0]);
        };
        let codes: Vec<&str> = warnings.iter().map(|w| w.code.as_str()).collect();
        assert_eq!(codes, ["filtered-psms", "missing-fasta"]);
    }

    #[test]
    fn fragpipe_folders_with_engine_results_skip_the_psm_table_warnings() {
        let dir = tempfile::tempdir().unwrap();
        let bundle = dir.path().join("bundle");
        copy_dir(&fixtures().join("fragpipe_synthetic"), &bundle);
        // The pyteomics sample stands in for an engine's own pepXML.
        std::fs::copy(
            fixtures().join("third_party/pyteomics/test.pep.xml"),
            bundle.join("control-1/run.pep.xml"),
        )
        .unwrap();
        let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();

        let outcomes = import_inputs(&mut project, &[bundle], None);

        let OutcomeKind::Imported {
            software,
            warnings,
            notes,
            ..
        } = &outcomes[0].kind
        else {
            panic!("{:?}", outcomes[0]);
        };
        assert_eq!(software, "SEQUEST");
        assert!(warnings.is_empty(), "{warnings:?}");
        assert!(
            notes[0].contains("SEQUEST's own results (1 pepXML file,"),
            "{notes:?}"
        );
    }

    fn copy_dir(from: &Path, to: &Path) {
        std::fs::create_dir_all(to).unwrap();
        for entry in std::fs::read_dir(from).unwrap() {
            let path = entry.unwrap().path();
            let target = to.join(path.file_name().unwrap());
            if path.is_dir() {
                copy_dir(&path, &target);
            } else {
                std::fs::copy(&path, &target).unwrap();
            }
        }
    }
}
