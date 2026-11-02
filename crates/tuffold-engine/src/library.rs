//! Spectral-library export (NIST `.msp`, which SpectraST, Skyline, and NIST
//! MS Search read): one entry per peptide ion that passes the saved
//! thresholds, its best-scoring passing spectrum (best replicate, not a
//! consensus), annotated with b and y ions.
//!
//! - `Name: SEQUENCE/charge`; modifications go in the comment as NIST's
//!   `Mods=count/position,residue,UnimodName/...` with 0-based positions.
//!   Terminal modifications sit on the end residue, as NIST writes them. A
//!   peptide carrying a modification without a Unimod name can't be written
//!   that way, so it is left out and counted.
//! - `Fullname=K.PEPTIDE.R/2` gives the residues around the peptide (from
//!   the leading protein's sequence, else the engine's report; left out
//!   when neither says), which SpectraST reads to count enzymatic termini.
//!   SpectraST ignores a `Prob=` in `.msp` (it sets 1), so none is written;
//!   `Qvalue=` carries Tuffold's.
//! - `MW` is the precursor m/z times the charge (as SpectraST writes it);
//!   `PrecursorMZ` and `Parent` are Tuffold's theoretical m/z.
//! - Peaks keep their measured m/z and intensity. Each is annotated with the
//!   b or y ion it matches within 20 ppm (`y4/0.0012`: ion, then observed −
//!   theoretical in Da; `^2` for doubly charged), else `?`.

use std::io::{self, Write};

use tuffold_core::chemistry::{mz, peptide_mass};
use tuffold_core::fragments::{Tolerance, annotate, fragment_ions};
use tuffold_core::{ModSite, PeptideForm};
use tuffold_store::{LibraryEntry, Project, StoreError};

/// Fragment tolerance for annotating library peaks.
const ANNOTATION_TOLERANCE: Tolerance = Tolerance::Ppm(20.0);

#[derive(Debug, thiserror::Error)]
pub enum LibraryError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error("this project has no current validation results; validate it first")]
    NotValidated,
}

/// What went into the library and what was left out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct LibraryReport {
    pub entries: usize,
    /// Passing peptide ions without any spectrum with peaks attached.
    pub without_peaks: usize,
    /// Peptide ions carrying a modification without a Unimod name.
    pub unnamed_modifications: usize,
}

/// Write the project's spectral library as NIST `.msp`.
pub fn write_msp(project: &Project, out: &mut impl Write) -> Result<LibraryReport, LibraryError> {
    let current = project
        .validation_summary()?
        .is_some_and(|summary| summary.validated_at.is_some());
    if !current {
        return Err(LibraryError::NotValidated);
    }
    let library = project.library()?;
    let mut report = LibraryReport {
        without_peaks: library.without_peaks,
        ..LibraryReport::default()
    };
    for entry in &library.entries {
        match entry_text(entry) {
            Some(text) => {
                out.write_all(text.as_bytes())?;
                report.entries += 1;
            }
            None => report.unnamed_modifications += 1,
        }
    }
    Ok(report)
}

/// One `.msp` entry, blank line included; `None` when a modification has no
/// name (or a residue no mass).
fn entry_text(entry: &LibraryEntry) -> Option<String> {
    let form = &entry.form;
    let mods = mods_field(form)?;
    let charge = u8::try_from(entry.charge).ok().filter(|&z| z > 0)?;
    let neutral = peptide_mass(&form.sequence)?
        + form.modifications.iter().map(|m| m.mass_delta).sum::<f64>();
    let precursor = mz(neutral, charge);
    let mut text = format!(
        "Name: {}/{charge}\nMW: {:.4}\nPrecursorMZ: {precursor:.4}\n",
        form.sequence,
        precursor * f64::from(charge)
    );
    let mut comment = vec![
        "Spec=BestReplicate".to_string(),
        format!("Mods={mods}"),
        format!("Parent={precursor:.4}"),
        format!("Charge={charge}"),
        format!("Nreps={}", entry.replicates),
        format!("Qvalue={:.6}", entry.q_value),
        format!("Score={:.4}", entry.score),
        format!("Source={}:{}", no_spaces(&entry.run), no_spaces(&entry.key)),
    ];
    if let Some(protein) = &entry.protein {
        comment.push(format!("Protein={}", no_spaces(protein)));
    }
    if let Some((before, after)) = entry.flanks {
        comment.push(format!(
            "Fullname={before}.{}.{after}/{charge}",
            form.sequence
        ));
    }
    if let Some(rt) = entry.rt_seconds {
        comment.push(format!("RetentionTime={rt:.2}"));
    }
    comment.push("Origin=Tuffold".to_string());
    text.push_str(&format!("Comment: {}\n", comment.join(" ")));
    let peaks = &entry.peaks;
    let annotations = annotations(form, charge, &peaks.mz, &peaks.intensity);
    text.push_str(&format!("Num peaks: {}\n", peaks.mz.len()));
    for ((mz, intensity), annotation) in peaks.mz.iter().zip(&peaks.intensity).zip(annotations) {
        text.push_str(&format!("{mz:.4}\t{intensity}\t\"{annotation}\"\n"));
    }
    text.push('\n');
    Some(text)
}

/// NIST's `Mods=` value: `0`, or `count/position,residue,name/...`.
fn mods_field(form: &PeptideForm) -> Option<String> {
    if form.modifications.is_empty() {
        return Some("0".into());
    }
    let residues: Vec<char> = form.sequence.chars().collect();
    let last = residues.len().checked_sub(1)?;
    let items: Option<Vec<String>> = form
        .modifications
        .iter()
        .map(|modification| {
            let name = modification
                .name
                .as_deref()
                .filter(|n| !n.contains([' ', ',', '/']))?;
            let position = match modification.site {
                ModSite::NTerm => 0,
                ModSite::CTerm => last,
                ModSite::Residue(position) => usize::try_from(position).ok()?.checked_sub(1)?,
            };
            Some(format!("{position},{},{name}", residues.get(position)?))
        })
        .collect();
    let items = items?;
    Some(format!("{}/{}", items.len(), items.join("/")))
}

/// One annotation per peak (`?` when no b or y ion matches it).
fn annotations(form: &PeptideForm, charge: u8, mz: &[f64], intensity: &[f32]) -> Vec<String> {
    let mut labels = vec!["?".to_string(); mz.len()];
    let Some(ions) = fragment_ions(form, Some(charge)) else {
        return labels;
    };
    for found in annotate(&ions, mz, intensity, ANNOTATION_TOLERANCE) {
        let ion = found.ion;
        let series = match ion.series {
            tuffold_core::fragments::IonSeries::B => 'b',
            tuffold_core::fragments::IonSeries::Y => 'y',
        };
        let charge = if ion.charge > 1 {
            format!("^{}", ion.charge)
        } else {
            String::new()
        };
        let label = format!(
            "{series}{}{charge}/{:.4}",
            ion.ordinal,
            mz[found.peak] - ion.mz
        );
        // A peak matching several ions lists them all, comma-separated.
        let slot = &mut labels[found.peak];
        if slot == "?" {
            *slot = label;
        } else {
            slot.push(',');
            slot.push_str(&label);
        }
    }
    labels
}

fn no_spaces(text: &str) -> String {
    text.replace(char::is_whitespace, "_")
}

#[cfg(test)]
mod tests {
    use tuffold_core::Modification;
    use tuffold_store::Peaks;

    use super::*;

    fn modification(site: ModSite, delta: f64, name: Option<&str>) -> Modification {
        Modification {
            site,
            mass_delta: delta,
            name: name.map(String::from),
            accession: None,
        }
    }

    #[test]
    fn mods_follow_nist_positions() {
        let form = PeptideForm::new(
            "MPEPTIDEK",
            vec![
                modification(ModSite::NTerm, 42.010565, Some("Acetyl")),
                modification(ModSite::Residue(1), 15.994915, Some("Oxidation")),
                modification(ModSite::CTerm, -0.984016, Some("Amidated")),
            ],
        )
        .unwrap();
        // N-term on residue 0 (M), oxidation on M at 0, C-term on K at 8.
        assert_eq!(
            mods_field(&form).unwrap(),
            "3/0,M,Acetyl/0,M,Oxidation/8,K,Amidated"
        );
        assert_eq!(
            mods_field(&PeptideForm::unmodified("PEPTIDEK")).unwrap(),
            "0"
        );
        let unnamed = PeptideForm::new(
            "PEPTIDEK",
            vec![modification(ModSite::Residue(2), 1.0, None)],
        )
        .unwrap();
        assert_eq!(mods_field(&unnamed), None);
    }

    #[test]
    fn an_entry_carries_mass_comment_and_annotated_peaks() {
        // PEPTIDEK: residues P 97.05276 E 129.04259 P T 101.04768 I 113.08406
        // D 115.02694 E K 128.09496, + H2O 18.01056 = 927.45493 Da;
        // 2+: (927.45493 + 2 × 1.00728) / 2 = 464.73474.
        // y1 = K + H2O + H+ = 128.09496 + 18.01056 + 1.00728 = 147.11280.
        let entry = LibraryEntry {
            form: PeptideForm::unmodified("PEPTIDEK"),
            charge: 2,
            run: "run 1".into(),
            key: "scan=7".into(),
            rt_seconds: Some(61.5),
            q_value: 0.001,
            score: 3.25,
            replicates: 4,
            protein: Some("sp|P1|ONE".into()),
            flanks: Some(('K', 'A')),
            peaks: Peaks {
                mz: vec![147.1130, 300.0],
                intensity: vec![100.0, 5.0],
            },
        };
        let text = entry_text(&entry).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "Name: PEPTIDEK/2");
        assert_eq!(lines[1], "MW: 929.4695");
        assert_eq!(lines[2], "PrecursorMZ: 464.7347");
        assert_eq!(
            lines[3],
            "Comment: Spec=BestReplicate Mods=0 Parent=464.7347 Charge=2 Nreps=4 \
             Qvalue=0.001000 Score=3.2500 Source=run_1:scan=7 Protein=sp|P1|ONE \
             Fullname=K.PEPTIDEK.A/2 RetentionTime=61.50 Origin=Tuffold"
        );
        assert_eq!(lines[4], "Num peaks: 2");
        // 147.1130 − 147.1128 = +0.0002.
        assert_eq!(lines[5], "147.1130\t100\t\"y1/0.0002\"");
        assert_eq!(lines[6], "300.0000\t5\t\"?\"");
        assert_eq!(lines[7], "");
    }
}
