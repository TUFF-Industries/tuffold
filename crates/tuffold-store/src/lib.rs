//! Tuffold project files (`.tuff`).
//!
//! A [`Project`] wraps one `.tuff` SQLite file: create it, import
//! [`ImportBatch`]es into it (deduplicating runs, spectra, peptides, and
//! proteins across imports), store validation results, and query summaries.
//! See [`schema`] for the table layout and versioning rules.

mod attach;
mod browse;
mod copy;
mod exchange;
mod import;
mod library;
mod mbr;
mod ms1;
mod peaks;
mod plex;
mod predict;
mod report;
pub use peaks::Peaks;
pub use predict::PsmPeaks;
pub mod schema;
mod summary;
mod validation;

use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use rusqlite::{Connection, OpenFlags, OptionalExtension};
use tuffold_core::{ImportBatch, SourceFormat};

pub use attach::{AttachReport, Attacher};
pub use browse::{
    AssignedPeptide, Biosample, EnginePsm, GroupDetail, GroupMember, GroupPeptide, ProteinList,
    ProteinListRow, SimilarGroup, SpectrumDetail, SpectrumListRow, SpectrumPage, StoredPeaks,
};
pub use copy::UndoNote;
pub use exchange::{
    Exchange, ExchangeEvidence, ExchangeGroup, ExchangePeptide, ExchangeProtein, ExchangeRun,
    ExchangeSpectrum,
};
pub use import::ImportReport;
pub use library::{Library, LibraryEntry};
pub use mbr::{MbrPrecursor, MbrRun, TRANSFER_FDR, Transfer, TransferSummary};
pub use ms1::{IntensityList, IntensityListRow, Ms1Run, PeptideIntensity, Precursor};
pub use plex::{PlexDefinition, ReporterList, ReporterListRow};
pub use report::{PeptideReportRow, Report, SpectrumReportRow};
pub use summary::{ImportRecord, ProjectNames, ProjectSummary, RunSummary};
pub use tuffold_validate::LocalizationSummary;
pub use validation::{
    EngineRecord, InferenceRows, LevelCounts, PeptideScoring, StoredPeptide, ValidationSummary,
};

#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    #[error("{0} already exists; choose a new project file name")]
    AlreadyExists(String),
    #[error("{0} is not a Tuffold project (.tuff) file")]
    NotAProject(String),
    #[error(
        "this project was written by a newer Tuffold (schema v{0}); \
         update Tuffold to open it"
    )]
    UnsupportedVersion(i32),
    #[error("{path} was already imported into this project as {format}")]
    AlreadyImported { path: String, format: &'static str },
    #[error(
        "no MS run named '{run}' in this project (known runs: {}); import the search \
         results first, or name the run explicitly",
        known.join(", ")
    )]
    UnknownRun { run: String, known: Vec<String> },
    #[error(
        "no biosample named '{name}' in this project (biosamples: {})",
        known.join(", ")
    )]
    UnknownBiosample { name: String, known: Vec<String> },
    #[error(
        "no plex named '{name}' in this project (plexes: {})",
        if known.is_empty() { "none; define one with `tuffold plex`".to_string() } else { known.join(", ") }
    )]
    UnknownPlex { name: String, known: Vec<String> },
    #[error("{plex} has no channel '{channel}' (channels: {})", known.join(", "))]
    UnknownChannel {
        channel: String,
        plex: &'static str,
        known: Vec<String>,
    },
    #[error("channel {0} is assigned twice")]
    DuplicateChannel(String),
    #[error("name the plex")]
    UnnamedPlex,
    #[error("reference channel {0} holds no biosample in this plex; assign it first")]
    UnassignedReference(String),
    #[error("the impurities given are for {given}, but the plex is {plex}")]
    ImpuritiesForOtherPlex {
        given: &'static str,
        plex: &'static str,
    },
    #[error("assign at least one channel to a biosample")]
    NoChannels,
    #[error("run {run} already belongs to plex {plex}")]
    RunInPlex { run: String, plex: String },
    #[error("{0} is in use by another program (or another Tuffold window); try again")]
    Busy(String),
    #[error("this project file is damaged: {0}")]
    Damaged(String),
    #[error("project database error: {0}")]
    Sqlite(#[from] rusqlite::Error),
}

/// Identifies the file an [`ImportBatch`] came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImportSource {
    /// Absolute path of the imported file or bundle directory.
    pub path: String,
    pub format: SourceFormat,
}

pub struct Project {
    conn: Connection,
}

impl Project {
    /// Create a new, empty project file. Refuses to overwrite.
    pub fn create(path: &Path) -> Result<Self, StoreError> {
        if path.exists() {
            return Err(StoreError::AlreadyExists(path.display().to_string()));
        }
        let conn = Connection::open(path)?;
        configure(&conn)?;
        schema::initialize(&conn)?;
        conn.execute(
            "INSERT INTO meta (key, value) VALUES ('created_with', ?1)",
            [concat!("tuffold ", env!("CARGO_PKG_VERSION"))],
        )?;
        Ok(Self { conn })
    }

    /// Open an existing project, migrating an older schema forward.
    pub fn open(path: &Path) -> Result<Self, StoreError> {
        let not_a_project = || StoreError::NotAProject(path.display().to_string());
        if !path.is_file() {
            return Err(not_a_project());
        }
        let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_WRITE)?;
        let application_id: i32 = conn
            .pragma_query_value(None, "application_id", |row| row.get(0))
            .map_err(|error| match error.sqlite_error_code() {
                // Another program is writing to it; it is still a project.
                Some(rusqlite::ErrorCode::DatabaseBusy | rusqlite::ErrorCode::DatabaseLocked) => {
                    StoreError::Busy(path.display().to_string())
                }
                _ => not_a_project(),
            })?;
        if application_id != schema::APPLICATION_ID {
            return Err(not_a_project());
        }
        let version: i32 = conn.pragma_query_value(None, "user_version", |row| row.get(0))?;
        if version > schema::SCHEMA_VERSION {
            return Err(StoreError::UnsupportedVersion(version));
        }
        configure(&conn)?;
        schema::migrate(&conn, version)?;
        Ok(Self { conn })
    }

    /// Open a project only to read it: nothing is migrated or written, so an
    /// older schema may fail its queries rather than be upgraded.
    pub fn open_read_only(path: &Path) -> Result<Self, StoreError> {
        let not_a_project = || StoreError::NotAProject(path.display().to_string());
        if !path.is_file() {
            return Err(not_a_project());
        }
        let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
        let application_id: i32 = conn
            .pragma_query_value(None, "application_id", |row| row.get(0))
            .map_err(|_| not_a_project())?;
        if application_id != schema::APPLICATION_ID {
            return Err(not_a_project());
        }
        let version: i32 = conn.pragma_query_value(None, "user_version", |row| row.get(0))?;
        if version > schema::SCHEMA_VERSION {
            return Err(StoreError::UnsupportedVersion(version));
        }
        Ok(Self { conn })
    }

    /// Import one batch atomically: either everything lands or nothing does.
    pub fn import(
        &mut self,
        source: &ImportSource,
        batch: &ImportBatch,
    ) -> Result<ImportReport, StoreError> {
        let mut reports = self.import_all(source, std::slice::from_ref(batch))?;
        Ok(reports.remove(0))
    }

    /// Import several batches from one source (a Scaffold project's search
    /// engines) in one transaction. Validation groups PSMs by import, so
    /// each batch is its own import, recorded as `<path>#<software>` when
    /// there are several.
    pub fn import_all(
        &mut self,
        source: &ImportSource,
        batches: &[ImportBatch],
    ) -> Result<Vec<ImportReport>, StoreError> {
        if batches.is_empty() {
            return Ok(Vec::new()); // nothing to record, nothing to invalidate
        }
        let already: Option<i64> = self
            .conn
            .query_row(
                "SELECT id FROM imports WHERE format = ?2 AND (source_path = ?1
                     OR substr(source_path, 1, length(?1) + 1) = ?1 || '#')",
                (&source.path, source.format.as_str()),
                |row| row.get(0),
            )
            .optional()?;
        if already.is_some() {
            return Err(StoreError::AlreadyImported {
                path: source.path.clone(),
                format: source.format.as_str(),
            });
        }
        let tx = self.conn.transaction()?;
        let now = unix_now();
        let mut reports = Vec::with_capacity(batches.len());
        for batch in batches {
            let source = if batches.len() > 1 {
                ImportSource {
                    path: format!("{}#{}", source.path, batch.software.name),
                    format: source.format,
                }
            } else {
                source.clone()
            };
            reports.push(import::write(&tx, &source, batch, now)?);
        }
        // Results computed before this import no longer describe the project.
        validation::invalidate(&tx)?;
        tx.commit()?;
        Ok(reports)
    }

    /// Start attaching measured spectra to the identified spectra of `run`.
    /// Nothing is written until [`Attacher::finish`].
    pub fn attach_spectra(&mut self, run: &str) -> Result<Attacher<'_>, StoreError> {
        Attacher::begin(self.conn.transaction()?, run, false)
    }

    /// Like [`Project::attach_spectra`], for partial peaks a project file
    /// stores (a Scaffold project keeps at most 300 per spectrum): they fill
    /// spectra without peaks, and a spectra file attached later replaces
    /// them.
    pub fn attach_partial_peaks(&mut self, run: &str) -> Result<Attacher<'_>, StoreError> {
        Attacher::begin(self.conn.transaction()?, run, true)
    }

    /// Protein groups passing the saved thresholds with their MS1 precursor
    /// intensities per biosample, and how much MS1 data each run has.
    pub fn intensity_list(&self) -> Result<IntensityList, StoreError> {
        ms1::intensity_list(&self.conn)
    }

    /// Define (or redefine) an isobaric plex: its runs and channel samples.
    /// Channel biosamples are created when new.
    pub fn set_plex(&mut self, definition: &PlexDefinition) -> Result<(), StoreError> {
        plex::set_plex(&mut self.conn, definition)
    }

    pub fn plexes(&self) -> Result<Vec<PlexDefinition>, StoreError> {
        plex::plexes(&self.conn)
    }

    /// Remove a plex (its runs are free for another plex again). Its
    /// channel biosamples stay; an unknown name is
    /// [`StoreError::UnknownPlex`].
    pub fn remove_plex(&mut self, name: &str) -> Result<(), StoreError> {
        plex::remove_plex(&self.conn, name)
    }

    /// The plex's passing PSMs per protein group with their reporter ions,
    /// leaving out spectra whose measured precursor purity is below
    /// `min_purity` (counted).
    pub fn reporter_list(
        &self,
        plex: &str,
        min_purity: Option<f64>,
    ) -> Result<ReporterList, StoreError> {
        plex::reporter_list(&self.conn, plex, min_purity)
    }

    /// Each run's identified precursors with their retention times and m/z,
    /// for match-between-runs.
    pub fn mbr_inputs(&mut self) -> Result<Vec<MbrRun>, StoreError> {
        mbr::inputs(&mut self.conn)
    }

    /// Replace a run's match-between-runs transfers.
    pub fn save_transfers(
        &mut self,
        run_id: i64,
        transfers: &[Transfer],
        summary: &TransferSummary,
    ) -> Result<(), StoreError> {
        mbr::save(&mut self.conn, run_id, transfers, summary)
    }

    /// Every peptide and spectrum passing the saved thresholds (targets only),
    /// for reports and exports.
    pub fn report(&self) -> Result<Report, StoreError> {
        report::report(&self.conn)
    }

    /// The stored peaks of one spectrum, if attached.
    pub fn peaks(
        &self,
        run: &str,
        key: &tuffold_core::SpectrumKey,
    ) -> Result<Option<Peaks>, StoreError> {
        attach::stored_peaks(&self.conn, run, &key.canonical())
    }

    pub fn summary(&self) -> Result<ProjectSummary, StoreError> {
        summary::project_summary(&self.conn)
    }

    pub fn imports(&self) -> Result<Vec<ImportRecord>, StoreError> {
        summary::imports(&self.conn)
    }

    /// Every sample, category, run, file, and protein name in the project.
    pub fn names(&self) -> Result<ProjectNames, StoreError> {
        summary::names(&self.conn)
    }

    /// Every validated spectrum, peptide, protein, and group, for results
    /// exports (see [`Exchange`]).
    pub fn exchange(&self) -> Result<Exchange, StoreError> {
        exchange::exchange(&self.conn)
    }

    /// The spectra a spectral library is built from (see [`Library`]).
    pub fn library(&self) -> Result<Library, StoreError> {
        library::library(&self.conn)
    }

    /// The rank-1 PSMs, peptides, and proteins validation works on.
    pub fn validation_input(&self) -> Result<tuffold_validate::Input, StoreError> {
        validation::input(&self.conn)
    }

    /// The peptides of the rank-1 PSMs, by id, for prediction.
    pub fn rank1_peptides(
        &self,
    ) -> Result<std::collections::HashMap<i64, tuffold_core::PeptideForm>, StoreError> {
        predict::rank1_peptides(&self.conn)
    }

    /// Whether any rank-1 PSM has peaks.
    pub fn has_rank1_peaks(&self) -> Result<bool, StoreError> {
        predict::has_rank1_peaks(&self.conn)
    }

    /// The rank-1 PSMs whose spectra have peaks: (PSM id, spectrum id).
    pub fn rank1_with_peaks(&self) -> Result<Vec<(i64, i64)>, StoreError> {
        predict::rank1_with_peaks(&self.conn)
    }

    /// Visit the rank-1 PSMs that have peaks, `chunk` at a time.
    pub fn for_each_rank1_spectrum<E: From<StoreError>>(
        &self,
        chunk: usize,
        visit: impl FnMut(&[PsmPeaks]) -> Result<(), E>,
    ) -> Result<(), E> {
        predict::for_each_rank1_spectrum(&self.conn, chunk, visit)
    }

    /// Visit the stored peaks (if any) of the requested spectra one at a time.
    pub fn for_each_spectrum_peaks(
        &self,
        spectrum_ids: &[i64],
        visit: impl FnMut(i64, Option<&Peaks>),
    ) -> Result<(), StoreError> {
        self.snapshot(|conn| predict::for_each_spectrum_peaks(conn, spectrum_ids, visit))
    }

    /// Replace all validation results (and save the settings) atomically.
    pub fn save_validation(
        &mut self,
        settings: &tuffold_validate::Settings,
        results: &tuffold_validate::Validation,
    ) -> Result<(), StoreError> {
        let tx = self.conn.transaction()?;
        validation::save(&tx, settings, results, unix_now())?;
        tx.commit()?;
        Ok(())
    }

    /// Validated peptides and the protein map, for inference at a new
    /// peptide threshold.
    pub fn inference_rows(&self) -> Result<InferenceRows, StoreError> {
        validation::inference_rows(&self.conn)
    }

    /// Replace the protein groups and save new thresholds atomically.
    pub fn save_inference(
        &mut self,
        thresholds: &tuffold_validate::Thresholds,
        inference: &tuffold_validate::ProteinInference,
    ) -> Result<(), StoreError> {
        let tx = self.conn.transaction()?;
        validation::save_inference(&tx, thresholds, inference)?;
        tx.commit()?;
        Ok(())
    }

    /// Save thresholds that need no recomputation (PSM and protein FDR,
    /// minimum peptides): reported counts follow them directly.
    pub fn save_thresholds(
        &mut self,
        thresholds: &tuffold_validate::Thresholds,
    ) -> Result<(), StoreError> {
        let tx = self.conn.transaction()?;
        validation::save_thresholds(&tx, thresholds)?;
        tx.commit()?;
        Ok(())
    }

    /// Protein groups passing the saved thresholds, with spectral counts per
    /// biosample and sequence coverage.
    pub fn protein_list(&self) -> Result<ProteinList, StoreError> {
        self.snapshot(browse::protein_list)
    }

    /// A protein group's members, peptides, and coverage; `None` when the id
    /// is not a current group.
    pub fn group_detail(&self, group_id: i64) -> Result<Option<GroupDetail>, StoreError> {
        self.snapshot(|conn| browse::group_detail(conn, group_id))
    }

    /// The validated spectra assigned a peptide, best first: at most `limit`
    /// of them, with the total so a capped list can say so.
    pub fn peptide_spectra(
        &self,
        peptide_id: i64,
        limit: usize,
    ) -> Result<SpectrumPage, StoreError> {
        self.snapshot(|conn| browse::peptide_spectra(conn, peptide_id, limit))
    }

    /// One spectrum: peaks, the assigned peptide, and every engine's PSM.
    pub fn spectrum_detail(&self, spectrum_id: i64) -> Result<Option<SpectrumDetail>, StoreError> {
        self.snapshot(|conn| browse::spectrum_detail(conn, spectrum_id))
    }

    /// Assign a biosample to a category (e.g. "control") for comparisons, or
    /// clear it with `None`. Results are unaffected; only grouping changes.
    pub fn set_category(
        &mut self,
        biosample: &str,
        category: Option<&str>,
    ) -> Result<(), StoreError> {
        self.set_categories(&[(biosample, category)])
    }

    /// Replace the project's GO slim annotations: `terms` are (id, name,
    /// aspect) and `annotations` (UniProt accession, term id). Like stars,
    /// these describe proteins, not data: validation stays as it is.
    pub fn set_go(
        &mut self,
        terms: &[(String, String, String)],
        annotations: &[(String, String)],
    ) -> Result<(), StoreError> {
        let tx = self.conn.transaction()?;
        tx.execute_batch("DELETE FROM go_annotations; DELETE FROM go_terms;")?;
        {
            let mut term =
                tx.prepare("INSERT INTO go_terms (id, name, aspect) VALUES (?1, ?2, ?3)")?;
            for (id, name, aspect) in terms {
                term.execute((id, name, aspect))?;
            }
            let mut annotation = tx
                .prepare("INSERT OR IGNORE INTO go_annotations (uniprot, go_id) VALUES (?1, ?2)")?;
            for (uniprot, go_id) in annotations {
                annotation.execute((uniprot, go_id))?;
            }
        }
        tx.commit()?;
        Ok(())
    }

    /// The UniProt accessions of the project's target proteins.
    pub fn uniprot_accessions(&self) -> Result<Vec<String>, StoreError> {
        let mut statement = self
            .conn
            .prepare("SELECT accession FROM proteins WHERE is_decoy = 0")?;
        let accessions = statement
            .query_map([], |row| row.get::<_, String>(0))?
            .collect::<Result<Vec<_>, _>>()?;
        let mut found: Vec<String> = accessions
            .iter()
            .filter_map(|a| tuffold_core::uniprot::uniprot_accession(a).map(str::to_string))
            .collect();
        found.sort();
        found.dedup();
        Ok(found)
    }

    /// Star or unstar a protein (by accession). Stars are the user's marks,
    /// not data: they survive validation and imports.
    pub fn set_starred(&mut self, accession: &str, starred: bool) -> Result<(), StoreError> {
        if starred {
            self.conn.execute(
                "INSERT OR IGNORE INTO starred_proteins (accession) VALUES (?1)",
                [accession],
            )?;
        } else {
            self.conn.execute(
                "DELETE FROM starred_proteins WHERE accession = ?1",
                [accession],
            )?;
        }
        Ok(())
    }

    /// The file this project is stored in.
    pub fn path(&self) -> Option<&Path> {
        self.conn
            .path()
            .map(Path::new)
            .filter(|path| !path.as_os_str().is_empty())
    }

    /// Write a consistent copy of the whole project to `to` (a new file).
    pub fn save_copy(&self, to: &Path) -> Result<(), StoreError> {
        copy::save_copy(&self.conn, to)
    }

    /// This project's identity, shared by its copies (created on first use).
    pub fn identity(&mut self) -> Result<String, StoreError> {
        copy::identity(&self.conn)
    }

    /// Run `work` while holding the project's exclusive lock.
    pub fn exclusive<T>(&self, work: impl FnOnce() -> T) -> Result<T, StoreError> {
        copy::exclusive(&self.conn, work)
    }

    /// The undo note of a saved copy, if this is one.
    pub fn undo_note(&self) -> Result<Option<UndoNote>, StoreError> {
        copy::undo_note(&self.conn)
    }

    /// Set (or clear) the undo note.
    pub fn set_undo_note(&mut self, note: Option<&UndoNote>) -> Result<(), StoreError> {
        copy::set_undo_note(&mut self.conn, note)
    }

    /// Several assignments at once: all apply, or (on an unknown biosample)
    /// none do.
    pub fn set_categories(
        &mut self,
        assignments: &[(&str, Option<&str>)],
    ) -> Result<(), StoreError> {
        let tx = self.conn.transaction()?;
        for (biosample, category) in assignments {
            let category = category.map(str::trim).filter(|c| !c.is_empty());
            let changed = tx.execute(
                "UPDATE biosamples SET category = ?2 WHERE name = ?1",
                rusqlite::params![biosample, category],
            )?;
            if changed == 0 {
                let known = tx
                    .prepare("SELECT name FROM biosamples ORDER BY name")?
                    .query_map([], |row| row.get(0))?
                    .collect::<Result<Vec<String>, _>>()?;
                return Err(StoreError::UnknownBiosample {
                    name: biosample.to_string(),
                    known,
                });
            }
        }
        tx.commit()?;
        Ok(())
    }

    /// Run several reads against one consistent snapshot, so a write from
    /// another connection (a validation or threshold change) can never land
    /// between them.
    fn snapshot<T>(
        &self,
        read: impl FnOnce(&Connection) -> Result<T, StoreError>,
    ) -> Result<T, StoreError> {
        let tx = self.conn.unchecked_transaction()?;
        let value = read(&tx)?;
        tx.commit()?;
        Ok(value)
    }

    /// Settings, engines, notes, and counts at the saved thresholds; `None`
    /// when the project was never validated.
    pub fn validation_summary(&self) -> Result<Option<ValidationSummary>, StoreError> {
        validation::summary(&self.conn)
    }

    /// PTM site localization summary at `probability_threshold` (for spectra
    /// passing the saved PSM threshold); `None` when not validated.
    pub fn localization_summary(
        &self,
        probability_threshold: f64,
    ) -> Result<Option<tuffold_validate::LocalizationSummary>, StoreError> {
        self.snapshot(|conn| {
            let Some((settings, validated_at)) = validation::settings(conn)? else {
                return Ok(None);
            };
            if validated_at.is_none() {
                return Ok(None);
            }
            validation::localization_summary(conn, &settings.thresholds, probability_threshold)
                .map(Some)
        })
    }
}

fn configure(conn: &Connection) -> Result<(), StoreError> {
    conn.pragma_update(None, "foreign_keys", true)?;
    conn.pragma_update(None, "journal_mode", "DELETE")?;
    Ok(())
}

fn unix_now() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| {
            i64::try_from(elapsed.as_secs()).unwrap_or(i64::MAX)
        })
}
