//! Tuffold workflows over the importers and the project store.
//!
//! The CLI, the desktop app, and (later) the MCP server all call these, so
//! every front end imports and validates data the same way.

mod deep;
pub mod demo;
mod diagnostics;
pub mod export;
pub mod go;
pub mod import;
pub mod library;
pub mod mbr;
pub mod methods;
pub mod mzid;
pub mod quant;
pub mod search;
pub mod undo;
pub mod validation;

pub use deep::{ModelsStatus, install_models, models_status};
pub use demo::{DemoError, create_demo_project};
pub use diagnostics::{diagnostics_report, install_panic_log, log_dir};
pub use export::{ExportError, MAX_ROWS, WorkbookReport, write_localization_tsv, write_workbook};
pub use go::{GoFiles, GoImportError, GoReport, annotate_go};
pub use import::{
    InputKind, Outcome, OutcomeKind, classify, duplicate_note, import_inputs, ms1_note, ms3_note,
};
pub use library::{LibraryError, LibraryReport, write_msp};
pub use mbr::{MbrError, MbrRunReport, match_between_runs};
pub use methods::{Methods, MethodsError, methods};
pub use mzid::{MzidError, MzidReport, write_mzid};
pub use quant::{
    IntensityQuant, QuantifyError, REPORTER_MIN_PURITY, ReporterQuant, quantify_intensity,
    quantify_project, quantify_reporters, quantify_reporters_with, write_anova_tsv,
    write_intensity_tsv, write_qc_tsv, write_quant_tsv, write_reporter_tsv,
};
pub use search::{
    Imported, RunReport, RunSearched, SearchFailure, SearchReport, prepare_search, search_files,
    search_prepared,
};
pub use tuffold_predict::MODELS_URL;
pub use tuffold_search::{
    FragmentTolerance, MIN_PEAKS, SearchOptions, SearchSettings, Searcher, Skipped,
};
pub use undo::{Checkpoint, UndoError, last_change, undo, undo_path};
pub use validation::{
    DeepLearning, ValidateError, apply_thresholds, validate_project, validate_project_with,
};

use tuffold_store::Project;

/// Assign biosamples to categories (`None` clears one), all or nothing, as an
/// undoable change.
pub fn set_categories(
    project: &mut Project,
    assignments: &[(&str, Option<&str>)],
) -> Result<(), UndoError> {
    let current: Vec<(String, Option<String>)> = project
        .protein_list()?
        .biosamples
        .into_iter()
        .map(|b| (b.name, b.category))
        .collect();
    let unchanged = assignments.iter().all(|(name, category)| {
        let wanted = category.map(str::trim).filter(|c| !c.is_empty());
        current
            .iter()
            .any(|(known, now)| known == name && now.as_deref() == wanted)
    });
    if unchanged {
        // Nothing to change, so nothing to undo: keep the last undo point.
        return Ok(());
    }
    let checkpoint = Checkpoint::take(project, "change sample groups")?;
    project.set_categories(assignments)?;
    if let Some(checkpoint) = checkpoint {
        checkpoint.keep()?;
    }
    Ok(())
}

/// Define (or redefine) an isobaric plex, with an undo point like sample
/// groups. Names are trimmed, blank channel assignments dropped, and
/// reference channels put in channel order here, so front ends pass what the
/// user typed. When `definition.runs` is empty,
/// every run in the project that is not in another plex belongs to this one.
pub fn set_plex(
    project: &mut Project,
    definition: &tuffold_store::PlexDefinition,
) -> Result<(), UndoError> {
    let mut definition = definition.clone();
    definition.name = definition.name.trim().to_string();
    if definition.runs.is_empty() {
        let taken: Vec<String> = project
            .plexes()?
            .into_iter()
            .filter(|p| p.name != definition.name)
            .flat_map(|p| p.runs)
            .collect();
        definition.runs = project
            .summary()?
            .runs
            .into_iter()
            .map(|r| r.name)
            .filter(|r| !taken.contains(r))
            .collect();
    }
    let order: Vec<String> = definition
        .plex
        .channels()
        .into_iter()
        .map(|c| c.name)
        .collect();
    definition.channels = definition
        .channels
        .into_iter()
        .map(|(c, s)| (c.trim().to_string(), s.trim().to_string()))
        .filter(|(c, s)| !c.is_empty() && !s.is_empty())
        .collect();
    definition
        .channels
        .sort_by_key(|(c, _)| order.iter().position(|o| o == c).unwrap_or(usize::MAX));
    let mut references: Vec<String> = definition
        .references
        .iter()
        .map(|r| r.trim().to_string())
        .filter(|r| !r.is_empty())
        .collect();
    references.sort_by_key(|c| order.iter().position(|o| o == c).unwrap_or(usize::MAX));
    references.dedup();
    definition.references = references;
    let mut canonical = definition.clone();
    canonical.runs.sort();
    if project.plexes()?.contains(&canonical) {
        return Ok(());
    }
    let checkpoint = Checkpoint::take(project, "define a plex")?;
    project.set_plex(&definition)?;
    if let Some(checkpoint) = checkpoint {
        checkpoint.keep()?;
    }
    Ok(())
}

/// Remove a named plex, with an undo point like [`set_plex`]. Its channel
/// biosamples stay; an unknown name is [`tuffold_store::StoreError::UnknownPlex`].
pub fn remove_plex(project: &mut Project, name: &str) -> Result<(), UndoError> {
    // An unknown name is refused before the undo copy (a whole project) is
    // made for nothing.
    let known: Vec<String> = project.plexes()?.into_iter().map(|p| p.name).collect();
    if !known.iter().any(|n| n == name.trim()) {
        return Err(tuffold_store::StoreError::UnknownPlex {
            name: name.trim().to_string(),
            known,
        }
        .into());
    }
    let checkpoint = Checkpoint::take(project, "remove a plex")?;
    project.remove_plex(name.trim())?;
    if let Some(checkpoint) = checkpoint {
        checkpoint.keep()?;
    }
    Ok(())
}

/// The lot certificate a plex redefined without one keeps: the stored plex's
/// of that name, if it has the same reagent set (a certificate lists one
/// set's channels). Front ends that don't edit the certificate pass this, so
/// fixing a channel never drops it.
pub fn kept_impurities(
    project: &Project,
    name: &str,
    plex: tuffold_core::reporter::Plex,
) -> Result<Option<tuffold_core::reporter::Impurities>, tuffold_store::StoreError> {
    Ok(stored_plex(project, name, plex)?.and_then(|p| p.impurities))
}

/// The reference channels a plex redefined without naming them keeps, as
/// [`kept_impurities`] does its certificate.
pub fn kept_references(
    project: &Project,
    name: &str,
    plex: tuffold_core::reporter::Plex,
) -> Result<Vec<String>, tuffold_store::StoreError> {
    Ok(stored_plex(project, name, plex)?.map_or_else(Vec::new, |p| p.references))
}

/// The stored plex of that name, if it has that reagent set.
fn stored_plex(
    project: &Project,
    name: &str,
    plex: tuffold_core::reporter::Plex,
) -> Result<Option<tuffold_store::PlexDefinition>, tuffold_store::StoreError> {
    Ok(project
        .plexes()?
        .into_iter()
        .find(|p| p.name == name.trim() && p.plex == plex))
}
