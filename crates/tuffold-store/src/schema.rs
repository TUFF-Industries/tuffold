//! The `.tuff` schema and its migrations.
//!
//! A `.tuff` file is a SQLite database marked with `application_id`
//! [`APPLICATION_ID`] ("TUFF") and versioned through `user_version`. Opening a
//! file migrates it forward one version at a time; a file written by a newer
//! Tuffold is refused rather than guessed at.
//!
//! Tables (schema v1):
//! - `meta(key, value)`: project-level facts (`created_with`, `name`).
//! - `imports`: one row per imported file (path, format, software, time).
//! - `biosamples` 1—n `ms_runs` 1—n `spectra`: the sample hierarchy. Spectra are
//!   unique per run by `spectrum_key` (`scan=N`, `index=N`, or a native id),
//!   which is how searches of the same run by different engines meet.
//! - `peptides` (unique ProForma key) 1—n `modifications`.
//! - `proteins` (unique accession, decoy flag, optional sequence).
//! - `psms`: one row per reported hit per import, linked to its spectrum and
//!   peptide; `psm_scores` holds every score verbatim; `psm_proteins` holds the
//!   protein evidence with positions and flanking residues (one row per locus).
//! - `psm_decoys` (view): a PSM is a decoy when it has evidence and every
//!   protein it maps to is a decoy. Derived, never stored, so decoy
//!   information from a later import applies to earlier PSMs too.
//!
//! Schema v2 adds validation, kept apart so imported source data is never
//! overwritten:
//! - `peptide_decoys` (view): a peptide is a decoy when it has evidence and
//!   every protein any of its PSMs maps to is a decoy; `psm_decoys` now
//!   follows its peptide, so one sequence is never both.
//! - `validation`: the single row of settings (rescoring, thresholds) and when
//!   results were computed; `validated_at` is NULL while they are stale.
//! - `validation_engines`, `validation_notes`: how each engine was ranked and
//!   what the user should know (`stage` = `psms` or `proteins`).
//! - `psm_results` (rank-1 PSMs within their engine), `spectrum_results` (the
//!   peptide each spectrum is assigned after engines vote), `peptide_results`.
//! - `protein_groups` 1—n `protein_group_members` (leader, equivalent,
//!   subset) and `protein_group_peptides` (unique, razor, shared).
//!
//! Importing new data clears the results (they no longer describe the
//! project) but keeps the settings.
//!
//! Schema v3 adds probability thresholds and protein probabilities:
//! - `validation.threshold_mode`: `fdr` (the default, as before) or
//!   `probability`. The `psm_fdr`, `peptide_fdr`, and `protein_fdr` columns
//!   hold the threshold for that mode: an FDR, or a probability.
//! - `protein_groups.pep`: the group's posterior error probability (NULL
//!   when it lost its picked competition, like `q_value`).
//!
//! Schema v4 adds MS1 precursor intensity, measured when a spectra file with
//! MS1 scans is attached:
//! - `ms_runs.ms1_scans`: MS1 scans read when the run's precursors were last
//!   traced; NULL when no MS1 data was ever read for the run.
//! - `precursor_features`: one row per traced precursor (peptide and charge)
//!   per run, with its integrated elution peak. `area` is NULL when the
//!   precursor was traced and no peak was found. Rows are measurements
//!   derived from the spectra file, replaced whole for the run whenever a
//!   file with MS1 scans is attached again; imports don't touch them.
//!
//! Schema v9 adds isobaric plexes (TMT): `plexes` (name, reagent set),
//! `plex_runs` (each run in at most one plex), and `plex_channels` (the
//! biosample each reporter channel holds). See `plex.rs`.
//!
//! Schema v16 adds PTM site localization results, computed during validation
//! and cleared by imports like other validation results:
//! - `ptm_spectrum_summary`: one row per validated spectrum with ambiguous
//!   variable modification sites (`scored`, `missing_peaks`, `capped`, or
//!   `unscorable`).
//! - `ptm_site_scores`: per-candidate-site probabilities, pairwise Ascores,
//!   and decoy-amino-acid control sites (`is_decoy = 1`) for FLR estimation.
//!
//! Schema v17 adds `spectra.ion_mobility`: the precursor's ion mobility
//! (1/K0) of an identified spectrum, from a timsTOF spectra file, recorded
//! when its peaks are attached; match-between-runs looks for a transfer at
//! its donors' mobility.

use rusqlite::Connection;

use crate::StoreError;

/// "TUFF" as a big-endian 32-bit integer.
pub const APPLICATION_ID: i32 = 0x5455_4646;
/// The schema version this build writes.
pub const SCHEMA_VERSION: i32 = 17;

const V1: &str = "
CREATE TABLE meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
) WITHOUT ROWID;

CREATE TABLE imports (
    id               INTEGER PRIMARY KEY,
    source_path      TEXT NOT NULL,
    format           TEXT NOT NULL,
    software         TEXT NOT NULL,
    software_version TEXT,
    imported_at      INTEGER NOT NULL, -- unix seconds
    UNIQUE (source_path, format)
);

CREATE TABLE biosamples (
    id       INTEGER PRIMARY KEY,
    name     TEXT NOT NULL UNIQUE,
    category TEXT
);

CREATE TABLE ms_runs (
    id           INTEGER PRIMARY KEY,
    biosample_id INTEGER NOT NULL REFERENCES biosamples(id),
    name         TEXT NOT NULL UNIQUE,
    spectra_file TEXT
);

CREATE TABLE spectra (
    id               INTEGER PRIMARY KEY,
    ms_run_id        INTEGER NOT NULL REFERENCES ms_runs(id),
    spectrum_key     TEXT NOT NULL,
    source_id        TEXT,
    precursor_mz     REAL,
    precursor_charge INTEGER,
    rt_seconds       REAL,
    peaks            BLOB,
    UNIQUE (ms_run_id, spectrum_key)
);

CREATE TABLE peptides (
    id       INTEGER PRIMARY KEY,
    proforma TEXT NOT NULL UNIQUE,
    sequence TEXT NOT NULL
);
CREATE INDEX peptides_sequence ON peptides(sequence);

CREATE TABLE modifications (
    peptide_id INTEGER NOT NULL REFERENCES peptides(id),
    ordinal    INTEGER NOT NULL,
    site_kind  TEXT NOT NULL CHECK (site_kind IN ('n_term', 'residue', 'c_term')),
    position   INTEGER, -- 1-based residue position for site_kind = 'residue'
    mass_delta REAL NOT NULL,
    name       TEXT,
    accession  TEXT,
    PRIMARY KEY (peptide_id, ordinal)
) WITHOUT ROWID;

CREATE TABLE proteins (
    id          INTEGER PRIMARY KEY,
    accession   TEXT NOT NULL UNIQUE,
    description TEXT,
    sequence    TEXT,
    is_decoy    INTEGER NOT NULL CHECK (is_decoy IN (0, 1)),
    -- 1 when a source flagged decoy status explicitly (mzIdentML isDecoy);
    -- explicit flags win over status inferred from the accession.
    decoy_explicit INTEGER NOT NULL DEFAULT 0 CHECK (decoy_explicit IN (0, 1))
);

CREATE TABLE psms (
    id              INTEGER PRIMARY KEY,
    import_id       INTEGER NOT NULL REFERENCES imports(id),
    spectrum_id     INTEGER NOT NULL REFERENCES spectra(id),
    peptide_id      INTEGER NOT NULL REFERENCES peptides(id),
    rank            INTEGER NOT NULL,
    charge          INTEGER NOT NULL,
    calculated_mz   REAL,
    experimental_mz REAL,
    pass_threshold  INTEGER CHECK (pass_threshold IN (0, 1))
);
CREATE INDEX psms_spectrum ON psms(spectrum_id);
CREATE INDEX psms_peptide ON psms(peptide_id);
CREATE INDEX psms_import ON psms(import_id);

CREATE TABLE psm_scores (
    psm_id    INTEGER NOT NULL REFERENCES psms(id),
    name      TEXT NOT NULL,
    accession TEXT,
    value     REAL NOT NULL,
    PRIMARY KEY (psm_id, name)
) WITHOUT ROWID;

CREATE TABLE psm_proteins (
    id         INTEGER PRIMARY KEY,
    psm_id     INTEGER NOT NULL REFERENCES psms(id),
    protein_id INTEGER NOT NULL REFERENCES proteins(id),
    start_pos  INTEGER,
    end_pos    INTEGER,
    pre        TEXT,
    post       TEXT,
    -- A peptide occurring twice in one protein is two loci (two rows).
    UNIQUE (psm_id, protein_id, start_pos)
);
CREATE INDEX psm_proteins_protein ON psm_proteins(protein_id);

CREATE VIEW psm_decoys AS
SELECT p.id AS psm_id,
       EXISTS (SELECT 1 FROM psm_proteins pp WHERE pp.psm_id = p.id)
       AND NOT EXISTS (
           SELECT 1 FROM psm_proteins pp JOIN proteins pr ON pr.id = pp.protein_id
           WHERE pp.psm_id = p.id AND pr.is_decoy = 0
       ) AS is_decoy
FROM psms p;
";

const V2: &str = "
DROP VIEW psm_decoys;

CREATE VIEW peptide_decoys AS
SELECT pe.id AS peptide_id,
       EXISTS (
           SELECT 1 FROM psms p JOIN psm_proteins pp ON pp.psm_id = p.id
           WHERE p.peptide_id = pe.id
       )
       AND NOT EXISTS (
           SELECT 1 FROM psms p
           JOIN psm_proteins pp ON pp.psm_id = p.id
           JOIN proteins pr ON pr.id = pp.protein_id
           WHERE p.peptide_id = pe.id AND pr.is_decoy = 0
       ) AS is_decoy
FROM peptides pe;

CREATE VIEW psm_decoys AS
SELECT p.id AS psm_id, d.is_decoy
FROM psms p JOIN peptide_decoys d ON d.peptide_id = p.peptide_id;

CREATE TABLE validation (
    id             INTEGER PRIMARY KEY CHECK (id = 1),
    rescore        INTEGER NOT NULL CHECK (rescore IN (0, 1)),
    psm_fdr        REAL NOT NULL,
    peptide_fdr    REAL NOT NULL,
    protein_fdr    REAL NOT NULL,
    min_peptides   INTEGER NOT NULL,
    validated_at   INTEGER, -- unix seconds; NULL while results are stale
    validated_with TEXT
);

CREATE TABLE validation_engines (
    ordinal  INTEGER PRIMARY KEY,
    name     TEXT NOT NULL,
    psms     INTEGER NOT NULL,
    decoys   INTEGER NOT NULL,
    method   TEXT NOT NULL CHECK (method IN ('rescored', 'best_score', 'unscored')),
    detail   TEXT NOT NULL,
    features TEXT NOT NULL, -- one name per line
    accepted INTEGER NOT NULL
);

CREATE TABLE validation_notes (
    id      INTEGER PRIMARY KEY,
    stage   TEXT NOT NULL CHECK (stage IN ('psms', 'proteins')),
    code    TEXT NOT NULL,
    message TEXT NOT NULL
);

CREATE TABLE psm_results (
    psm_id  INTEGER PRIMARY KEY REFERENCES psms(id),
    score   REAL NOT NULL,
    q_value REAL, -- NULL when another PSM of the engine won the spectrum
    pep     REAL
);

CREATE TABLE spectrum_results (
    spectrum_id INTEGER PRIMARY KEY REFERENCES spectra(id),
    peptide_id  INTEGER NOT NULL REFERENCES peptides(id),
    is_decoy    INTEGER NOT NULL CHECK (is_decoy IN (0, 1)),
    score       REAL NOT NULL,
    q_value     REAL NOT NULL,
    pep         REAL NOT NULL,
    engines     INTEGER NOT NULL
);
CREATE INDEX spectrum_results_peptide ON spectrum_results(peptide_id);

CREATE TABLE peptide_results (
    peptide_id       INTEGER PRIMARY KEY REFERENCES peptides(id),
    is_decoy         INTEGER NOT NULL CHECK (is_decoy IN (0, 1)),
    score            REAL NOT NULL,
    q_value          REAL NOT NULL,
    pep              REAL NOT NULL,
    best_spectrum_id INTEGER NOT NULL REFERENCES spectra(id)
);

CREATE TABLE protein_groups (
    id                INTEGER PRIMARY KEY,
    leader_protein_id INTEGER NOT NULL REFERENCES proteins(id),
    is_decoy          INTEGER NOT NULL CHECK (is_decoy IN (0, 1)),
    score             REAL NOT NULL,
    q_value           REAL, -- NULL when the group lost to its target-decoy partner
    peptides          INTEGER NOT NULL -- distinct unique + razor sequences
);

CREATE TABLE protein_group_members (
    group_id   INTEGER NOT NULL REFERENCES protein_groups(id),
    protein_id INTEGER NOT NULL REFERENCES proteins(id),
    role       TEXT NOT NULL CHECK (role IN ('leader', 'equivalent', 'subset')),
    PRIMARY KEY (group_id, protein_id)
) WITHOUT ROWID;
CREATE INDEX protein_group_members_protein ON protein_group_members(protein_id);

CREATE TABLE protein_group_peptides (
    group_id   INTEGER NOT NULL REFERENCES protein_groups(id),
    peptide_id INTEGER NOT NULL REFERENCES peptides(id),
    role       TEXT NOT NULL CHECK (role IN ('unique', 'razor', 'shared')),
    PRIMARY KEY (group_id, peptide_id)
) WITHOUT ROWID;
CREATE INDEX protein_group_peptides_peptide ON protein_group_peptides(peptide_id);
";

const V3: &str = "
ALTER TABLE validation ADD COLUMN threshold_mode TEXT NOT NULL DEFAULT 'fdr'
    CHECK (threshold_mode IN ('fdr', 'probability'));
ALTER TABLE protein_groups ADD COLUMN pep REAL;
";

const V4: &str = "
ALTER TABLE ms_runs ADD COLUMN ms1_scans INTEGER;

CREATE TABLE precursor_features (
    ms_run_id      INTEGER NOT NULL REFERENCES ms_runs(id),
    peptide_id     INTEGER NOT NULL REFERENCES peptides(id),
    charge         INTEGER NOT NULL,
    area           REAL, -- NULL: traced, no peak found
    apex_rt        REAL,
    apex_intensity REAL,
    start_rt       REAL,
    end_rt         REAL,
    points         INTEGER,
    PRIMARY KEY (ms_run_id, peptide_id, charge)
) WITHOUT ROWID;
CREATE INDEX precursor_features_peptide ON precursor_features(peptide_id);
";

/// Proteins the user starred, by accession (so a star survives
/// validating again when the protein still leads a group).
const V5: &str = "
CREATE TABLE starred_proteins (
    accession TEXT PRIMARY KEY
) WITHOUT ROWID;
";

/// GO slim terms (`aspect`: process, function, component) and which UniProt
/// accessions reach them, replaced whole by each GO import.
const V6: &str = "
CREATE TABLE go_terms (
    id     TEXT PRIMARY KEY,
    name   TEXT NOT NULL,
    aspect TEXT NOT NULL CHECK (aspect IN ('process', 'function', 'component'))
) WITHOUT ROWID;
CREATE TABLE go_annotations (
    uniprot TEXT NOT NULL,
    go_id   TEXT NOT NULL REFERENCES go_terms(id),
    PRIMARY KEY (uniprot, go_id)
) WITHOUT ROWID;
";

/// How peptides were scored: `pooled` over their spectra (with the
/// features used), or by their best spectrum (`detail` says why); NULL
/// until validated.
const V7: &str = "
ALTER TABLE validation ADD COLUMN peptide_method TEXT
    CHECK (peptide_method IN ('pooled', 'best_spectrum'));
ALTER TABLE validation ADD COLUMN peptide_detail TEXT;
ALTER TABLE validation ADD COLUMN peptide_features TEXT; -- one name per line
";

/// Peaks a project file stored for its spectra (Scaffold keeps at most 300):
/// a spectra file attached later replaces them, unlike complete peaks.
const V8: &str = "
ALTER TABLE spectra ADD COLUMN peaks_partial INTEGER NOT NULL DEFAULT 0
    CHECK (peaks_partial IN (0, 1));
";

/// Isobaric plexes: the runs of one TMT pool and the biosample each
/// reporter channel holds.
const V9: &str = "
CREATE TABLE plexes (
    id    INTEGER PRIMARY KEY,
    name  TEXT NOT NULL UNIQUE,
    label TEXT NOT NULL -- the reagent set, e.g. 'tmt10'
);
CREATE TABLE plex_runs (
    ms_run_id INTEGER PRIMARY KEY REFERENCES ms_runs(id),
    plex_id   INTEGER NOT NULL REFERENCES plexes(id) ON DELETE CASCADE
);
CREATE TABLE plex_channels (
    plex_id      INTEGER NOT NULL REFERENCES plexes(id) ON DELETE CASCADE,
    channel      TEXT NOT NULL,
    biosample_id INTEGER NOT NULL REFERENCES biosamples(id),
    PRIMARY KEY (plex_id, channel)
) WITHOUT ROWID;
";

/// The score a peptide's protein group ranks by, when it isn't the peptide's
/// score (decoys' sibling counts imputed); NULL in older results means the
/// score.
const V10: &str = "
ALTER TABLE peptide_results ADD COLUMN group_score REAL;
";

/// A plex's reagent lot certificate (`Impurities::to_text`), when given:
/// reporters are corrected for isotope impurities with it.
const V11: &str = "
ALTER TABLE plexes ADD COLUMN impurities TEXT;
";

/// Channels holding the pooled reference sample that bridges plexes.
const V12: &str = "
ALTER TABLE plex_channels ADD COLUMN is_reference INTEGER NOT NULL DEFAULT 0
    CHECK (is_reference IN (0, 1));
";

/// Peaks of the MS3 scan acquired from an identified MS/MS spectrum
/// (SPS-MS3: the TMT reporters are read there).
const V13: &str = "
ALTER TABLE spectra ADD COLUMN ms3_peaks BLOB;
";

/// The share of the isolation window that was the identified precursor, in
/// the MS1 scans around the spectrum (NULL: not measured).
const V14: &str = "
ALTER TABLE spectra ADD COLUMN precursor_purity REAL;
";

/// Match-between-runs: precursors identified elsewhere, measured where they
/// should elute in a run that didn't identify them, with the q-value of the
/// target-decoy competition that judged them; and what matching did per run.
const V15: &str = "
CREATE TABLE transfers (
    ms_run_id    INTEGER NOT NULL REFERENCES ms_runs(id),
    peptide_id   INTEGER NOT NULL REFERENCES peptides(id),
    charge       INTEGER NOT NULL,
    predicted_rt REAL NOT NULL,
    area         REAL,
    apex_rt      REAL,
    q_value      REAL,
    PRIMARY KEY (ms_run_id, peptide_id, charge)
) WITHOUT ROWID;
CREATE TABLE transfer_runs (
    ms_run_id   INTEGER PRIMARY KEY REFERENCES ms_runs(id),
    candidates  INTEGER NOT NULL,
    known       INTEGER NOT NULL,
    trained     INTEGER NOT NULL CHECK (trained IN (0, 1)),
    targets_won INTEGER NOT NULL,
    decoys_won  INTEGER NOT NULL,
    accepted    INTEGER NOT NULL
);
";

/// PTM site localization per spectrum and candidate site, replaced by
/// validation and cleared by imports.
const V16: &str = "
CREATE TABLE ptm_spectrum_summary (
    spectrum_id        INTEGER PRIMARY KEY REFERENCES spectra(id),
    status             TEXT NOT NULL CHECK (status IN ('scored', 'missing_peaks', 'capped', 'unscorable')),
    isoform_count      INTEGER NOT NULL,
    localized_peptide  TEXT,
    best_ascore        REAL,
    min_probability    REAL,
    false_target_sites INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE ptm_site_scores (
    spectrum_id         INTEGER NOT NULL REFERENCES spectra(id),
    is_decoy            INTEGER NOT NULL DEFAULT 0 CHECK (is_decoy IN (0, 1)),
    position            INTEGER NOT NULL,
    residue             TEXT NOT NULL,
    mass_delta          REAL NOT NULL,
    modification_name   TEXT,
    probability         REAL NOT NULL,
    ascore              REAL NOT NULL,
    determining_matched INTEGER NOT NULL,
    determining_total   INTEGER NOT NULL,
    assigned            INTEGER NOT NULL CHECK (assigned IN (0, 1)),
    PRIMARY KEY (spectrum_id, is_decoy, position, mass_delta)
) WITHOUT ROWID;
CREATE INDEX ptm_site_scores_spectrum ON ptm_site_scores(spectrum_id);
";

/// The precursor's ion mobility (1/K0) of an identified spectrum, when its
/// spectra file has one (timsTOF).
const V17: &str = "
ALTER TABLE spectra ADD COLUMN ion_mobility REAL;
";

/// Migration `n` upgrades a file from version `n` to `n + 1`.
const MIGRATIONS: &[&str] = &[
    V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17,
];

pub(crate) fn initialize(conn: &Connection) -> Result<(), StoreError> {
    conn.pragma_update(None, "application_id", APPLICATION_ID)?;
    migrate(conn, 0)
}

pub(crate) fn migrate(conn: &Connection, from: i32) -> Result<(), StoreError> {
    let from = usize::try_from(from).map_err(|_| StoreError::UnsupportedVersion(from))?;
    for (version, sql) in MIGRATIONS.iter().enumerate().skip(from) {
        let tx = conn.unchecked_transaction()?;
        tx.execute_batch(sql)?;
        let next = i32::try_from(version + 1).expect("few migrations");
        tx.pragma_update(None, "user_version", next)?;
        tx.commit()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_version_1_file_migrates_and_keeps_its_data() {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute_batch(V1).unwrap();
        conn.pragma_update(None, "user_version", 1).unwrap();
        conn.execute_batch(
            "INSERT INTO imports VALUES (1, 'a.mzid', 'mzidentml', 'X', NULL, 0);
             INSERT INTO biosamples VALUES (1, 'b', NULL);
             INSERT INTO ms_runs VALUES (1, 1, 'r', NULL);
             INSERT INTO spectra (id, ms_run_id, spectrum_key) VALUES (1, 1, 'scan=1');
             INSERT INTO peptides VALUES (1, 'PEPTIDEK', 'PEPTIDEK');
             INSERT INTO proteins (id, accession, is_decoy) VALUES (1, 'rev_P1', 1);
             INSERT INTO psms (id, import_id, spectrum_id, peptide_id, rank, charge)
                 VALUES (1, 1, 1, 1, 1, 2);
             INSERT INTO psm_proteins (psm_id, protein_id) VALUES (1, 1);",
        )
        .unwrap();

        migrate(&conn, 1).unwrap();

        let version: i32 = conn
            .pragma_query_value(None, "user_version", |row| row.get(0))
            .unwrap();
        assert_eq!(version, SCHEMA_VERSION);
        let decoy: bool = conn
            .query_row(
                "SELECT is_decoy FROM psm_decoys WHERE psm_id = 1",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert!(decoy, "the rebuilt view still reads the old data");
        let results: i64 = conn
            .query_row("SELECT COUNT(*) FROM protein_groups", [], |row| row.get(0))
            .unwrap();
        assert_eq!(results, 0);
    }

    #[test]
    fn version_2_settings_keep_reading_as_fdr() {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute_batch(V1).unwrap();
        conn.execute_batch(V2).unwrap();
        conn.pragma_update(None, "user_version", 2).unwrap();
        conn.execute(
            "INSERT INTO validation (id, rescore, psm_fdr, peptide_fdr, protein_fdr, min_peptides)
             VALUES (1, 1, 0.01, 0.02, 0.05, 2)",
            [],
        )
        .unwrap();
        migrate(&conn, 2).unwrap();
        let (mode, protein): (String, f64) = conn
            .query_row(
                "SELECT threshold_mode, protein_fdr FROM validation",
                [],
                |row| Ok((row.get(0)?, row.get(1)?)),
            )
            .unwrap();
        assert_eq!((mode.as_str(), protein), ("fdr", 0.05));
    }
}
