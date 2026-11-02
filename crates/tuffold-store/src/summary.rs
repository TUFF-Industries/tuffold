//! Read-side summaries of a project.

use rusqlite::Connection;

use crate::StoreError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunSummary {
    pub name: String,
    pub biosample: String,
    pub spectra: usize,
    pub psms: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ProjectSummary {
    pub imports: usize,
    pub biosamples: usize,
    pub ms_runs: usize,
    pub spectra: usize,
    pub psms: usize,
    pub decoy_psms: usize,
    pub peptides: usize,
    pub proteins: usize,
    pub decoy_proteins: usize,
    /// Proteins with a stored sequence (coverage is possible for these).
    pub proteins_with_sequence: usize,
    pub runs: Vec<RunSummary>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImportRecord {
    pub id: i64,
    pub source_path: String,
    pub format: String,
    pub software: String,
    pub software_version: Option<String>,
    pub imported_at: i64,
    pub psms: usize,
}

fn count(conn: &Connection, sql: &str) -> Result<usize, StoreError> {
    let value: i64 = conn.query_row(sql, [], |row| row.get(0))?;
    Ok(usize::try_from(value).unwrap_or_default())
}

pub(crate) fn project_summary(conn: &Connection) -> Result<ProjectSummary, StoreError> {
    let mut statement = conn.prepare(
        "SELECT r.name, b.name,
                (SELECT COUNT(*) FROM spectra s WHERE s.ms_run_id = r.id),
                (SELECT COUNT(*) FROM psms p JOIN spectra s ON s.id = p.spectrum_id
                  WHERE s.ms_run_id = r.id)
         FROM ms_runs r JOIN biosamples b ON b.id = r.biosample_id
         ORDER BY r.name",
    )?;
    let runs = statement
        .query_map([], |row| {
            Ok(RunSummary {
                name: row.get(0)?,
                biosample: row.get(1)?,
                spectra: row.get::<_, i64>(2)?.try_into().unwrap_or_default(),
                psms: row.get::<_, i64>(3)?.try_into().unwrap_or_default(),
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ProjectSummary {
        imports: count(conn, "SELECT COUNT(*) FROM imports")?,
        biosamples: count(conn, "SELECT COUNT(*) FROM biosamples")?,
        ms_runs: count(conn, "SELECT COUNT(*) FROM ms_runs")?,
        spectra: count(conn, "SELECT COUNT(*) FROM spectra")?,
        psms: count(conn, "SELECT COUNT(*) FROM psms")?,
        decoy_psms: count(conn, "SELECT COUNT(*) FROM psm_decoys WHERE is_decoy = 1")?,
        peptides: count(conn, "SELECT COUNT(*) FROM peptides")?,
        proteins: count(conn, "SELECT COUNT(*) FROM proteins")?,
        decoy_proteins: count(conn, "SELECT COUNT(*) FROM proteins WHERE is_decoy = 1")?,
        proteins_with_sequence: count(
            conn,
            "SELECT COUNT(*) FROM proteins WHERE sequence IS NOT NULL",
        )?,
        runs,
    })
}

pub(crate) fn imports(conn: &Connection) -> Result<Vec<ImportRecord>, StoreError> {
    let mut statement = conn.prepare(
        "SELECT i.id, i.source_path, i.format, i.software, i.software_version, i.imported_at,
                (SELECT COUNT(*) FROM psms p WHERE p.import_id = i.id)
         FROM imports i ORDER BY i.id",
    )?;
    let records = statement
        .query_map([], |row| {
            Ok(ImportRecord {
                id: row.get(0)?,
                source_path: row.get(1)?,
                format: row.get(2)?,
                software: row.get(3)?,
                software_version: row.get(4)?,
                imported_at: row.get(5)?,
                psms: row.get::<_, i64>(6)?.try_into().unwrap_or_default(),
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    Ok(records)
}

/// Every name a project holds that could say what or whose the data is:
/// for privacy aliases, which must know them all before a reply leaves.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ProjectNames {
    /// Biosample names, by id.
    pub biosamples: Vec<String>,
    /// Distinct categories, sorted.
    pub categories: Vec<String>,
    /// Plex names, by id.
    pub plexes: Vec<String>,
    /// Run names, by id.
    pub runs: Vec<String>,
    /// Spectra files attached to runs, by run id.
    pub spectra_files: Vec<String>,
    /// Import source paths, by id.
    pub sources: Vec<String>,
    /// Protein accessions, decoys included, by id.
    pub accessions: Vec<String>,
}

pub(crate) fn names(conn: &Connection) -> Result<ProjectNames, StoreError> {
    let column = |sql: &str| -> Result<Vec<String>, StoreError> {
        let mut statement = conn.prepare(sql)?;
        let values = statement
            .query_map([], |row| row.get(0))?
            .collect::<Result<Vec<String>, _>>()?;
        Ok(values)
    };
    Ok(ProjectNames {
        biosamples: column("SELECT name FROM biosamples ORDER BY id")?,
        categories: column(
            "SELECT DISTINCT category FROM biosamples WHERE category IS NOT NULL
             ORDER BY category",
        )?,
        plexes: column("SELECT name FROM plexes ORDER BY id")?,
        runs: column("SELECT name FROM ms_runs ORDER BY id")?,
        spectra_files: column(
            "SELECT spectra_file FROM ms_runs WHERE spectra_file IS NOT NULL ORDER BY id",
        )?,
        sources: column("SELECT source_path FROM imports ORDER BY id")?,
        accessions: column("SELECT accession FROM proteins ORDER BY id")?,
    })
}
