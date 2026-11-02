//! The Excel workbook export: everything a validated project reports, one
//! sheet per kind of thing, at the saved thresholds.
//!
//! - **Summary**: the project, the thresholds, what passed at each level,
//!   each engine, and each run.
//! - **Methods**: the publication paragraph and its references
//!   ([`crate::methods`]).
//! - **Proteins**, **Peptides**, **Spectra**: every passing target row.
//! - **ANOVA (counts)** and **ANOVA (MS1)**, when two or more categories have
//!   data (MS1 only with MS1 data).
//!
//! Unknown values are empty cells, never zeros. Headers are frozen and
//! filterable. Excel holds 1,048,575 data rows per sheet; a longer list is
//! cut there and the Summary says so (rows are never dropped silently).

use std::path::Path;

use rust_xlsxwriter::{Format, Workbook, Worksheet, XlsxError};
use tuffold_quant::anova::AnovaComparison;
use tuffold_store::{Project, StoreError};
use tuffold_validate::Level;

use crate::methods::{MethodsError, methods};
use crate::quant::{QuantifyError, quantify_intensity, quantify_project};

/// Data rows a worksheet can hold below its header.
pub const MAX_ROWS: usize = 1_048_575;

#[derive(Debug, thiserror::Error)]
pub enum ExportError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Quantify(#[from] QuantifyError),
    #[error(transparent)]
    Methods(#[from] MethodsError),
    #[error("{path}: {message}")]
    Write { path: String, message: String },
}

/// What was written: each sheet and its data rows, and any sheet cut at
/// Excel's row limit (name, rows the project had).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct WorkbookReport {
    pub sheets: Vec<(String, usize)>,
    pub cut: Vec<(String, usize)>,
}

/// Cell formats used across sheets.
struct Formats {
    header: Format,
    wrap: Format,
    scientific: Format,
    percent: Format,
    fixed2: Format,
    fixed4: Format,
}

impl Formats {
    fn new() -> Self {
        Self {
            header: Format::new().set_bold(),
            wrap: Format::new().set_text_wrap(),
            scientific: Format::new().set_num_format("0.00E+00"),
            percent: Format::new().set_num_format("0.0%"),
            fixed2: Format::new().set_num_format("0.00"),
            fixed4: Format::new().set_num_format("0.0000"),
        }
    }
}

/// One cell of a data row.
enum Cell<'a> {
    Text(&'a str),
    Owned(String),
    Count(usize),
    Number(Option<f64>, Option<&'a Format>),
}

/// Writes a header row, then data rows, into one sheet.
struct Sheet<'w> {
    sheet: &'w mut Worksheet,
    row: u32,
    width: u16,
}

impl<'w> Sheet<'w> {
    fn new(
        workbook: &'w mut Workbook,
        name: &str,
        header: &[String],
        formats: &Formats,
    ) -> Result<Self, XlsxError> {
        let sheet = workbook.add_worksheet();
        sheet.set_name(name)?;
        for (col, title) in header.iter().enumerate() {
            let col = u16::try_from(col).unwrap_or(u16::MAX);
            sheet.write_string_with_format(0, col, title, &formats.header)?;
            let width = (title.chars().count() as f64 + 2.0).clamp(9.0, 40.0);
            sheet.set_column_width(col, width)?;
        }
        sheet.set_freeze_panes(1, 0)?;
        Ok(Self {
            sheet,
            row: 1,
            width: u16::try_from(header.len()).unwrap_or(u16::MAX),
        })
    }

    fn push(&mut self, cells: &[Cell<'_>]) -> Result<(), XlsxError> {
        for (col, cell) in cells.iter().enumerate() {
            let col = u16::try_from(col).unwrap_or(u16::MAX);
            match cell {
                Cell::Text(text) => {
                    self.sheet.write_string(self.row, col, *text)?;
                }
                Cell::Owned(text) => {
                    self.sheet.write_string(self.row, col, text)?;
                }
                Cell::Count(n) => {
                    self.sheet.write_number(self.row, col, *n as f64)?;
                }
                Cell::Number(Some(value), Some(format)) if value.is_finite() => {
                    self.sheet
                        .write_number_with_format(self.row, col, *value, format)?;
                }
                Cell::Number(Some(value), None) if value.is_finite() => {
                    self.sheet.write_number(self.row, col, *value)?;
                }
                Cell::Number(..) => {} // unknown: an empty cell
            }
        }
        self.row += 1;
        Ok(())
    }

    /// Filter buttons over the header; returns the data rows written.
    fn finish(self) -> Result<usize, XlsxError> {
        if self.width > 0 {
            self.sheet
                .autofilter(0, 0, self.row.saturating_sub(1).max(1), self.width - 1)?;
        }
        Ok(self.row as usize - 1)
    }
}

/// Write the workbook of a validated project to `path`.
pub fn write_workbook(project: &Project, path: &Path) -> Result<WorkbookReport, ExportError> {
    let write_error = |error: XlsxError| ExportError::Write {
        path: path.display().to_string(),
        message: error.to_string(),
    };
    let counts = quantify_project(project)?;
    let intensity = quantify_intensity(project).ok();
    let methods = methods(project)?;
    let list = project.protein_list()?;
    let report = project.report()?;
    let formats = Formats::new();
    let mut workbook = Workbook::new();
    let mut written = WorkbookReport::default();
    let mut record = |name: &str, total: usize, rows: usize| {
        written.sheets.push((name.to_string(), rows));
        if rows < total {
            written.cut.push((name.to_string(), total));
        }
    };

    // Summary first, filled in last (it reports the cuts); reserve it.
    workbook
        .add_worksheet()
        .set_name("Summary")
        .map_err(write_error)?;

    // Methods.
    {
        let sheet = workbook.add_worksheet();
        sheet.set_name("Methods").map_err(write_error)?;
        sheet.set_column_width(0, 120).map_err(write_error)?;
        let mut row = 0;
        for paragraph in &methods.paragraphs {
            sheet
                .write_string_with_format(row, 0, paragraph, &formats.wrap)
                .map_err(write_error)?;
            row += 2;
        }
        sheet
            .write_string_with_format(row, 0, "References", &formats.header)
            .map_err(write_error)?;
        for (n, reference) in methods.references.iter().enumerate() {
            row += 1;
            sheet
                .write_string_with_format(row, 0, format!("{}. {reference}", n + 1), &formats.wrap)
                .map_err(write_error)?;
        }
    }

    // Proteins.
    let samples: Vec<String> = list.biosamples.iter().map(|b| b.name.clone()).collect();
    let ms1 = intensity
        .as_ref()
        .filter(|q| q.runs.iter().any(|r| r.ms1_scans.unwrap_or(0) > 0));
    let mut header: Vec<String> = [
        "Accession",
        "Description",
        "Equivalent proteins",
        "Peptides",
        "Spectra",
        "Coverage",
        "Length",
        "q-value",
        "Probability",
    ]
    .map(String::from)
    .to_vec();
    for prefix in [
        "Spectra",
        "Exclusive spectra",
        "Exclusive peptides",
        "Normalized spectra",
        "NSAF",
        "emPAI",
    ] {
        header.extend(samples.iter().map(|s| format!("{prefix}: {s}")));
    }
    if ms1.is_some() {
        for prefix in ["MS1 intensity", "Top 3"] {
            header.extend(samples.iter().map(|s| format!("{prefix}: {s}")));
        }
    }
    {
        let mut sheet =
            Sheet::new(&mut workbook, "Proteins", &header, &formats).map_err(write_error)?;
        sheet.sheet.set_column_width(1, 48).map_err(write_error)?;
        for row in list.rows.iter().take(MAX_ROWS) {
            let quant = counts.rows.iter().find(|q| q.group_id == row.group_id);
            let lfq = ms1.and_then(|q| {
                q.quantities
                    .rows
                    .iter()
                    .find(|r| r.group_id == row.group_id)
            });
            let mut cells = vec![
                Cell::Text(&row.accession),
                Cell::Text(row.description.as_deref().unwrap_or("")),
                Cell::Count(row.equivalents),
                Cell::Count(row.peptides),
                Cell::Count(row.spectra),
                Cell::Number(row.coverage, Some(&formats.percent)),
                Cell::Number(row.length.map(|l| l as f64), None),
                Cell::Number(row.q_value, Some(&formats.scientific)),
                Cell::Number(row.probability, Some(&formats.percent)),
            ];
            cells.extend(row.spectra_by_sample.iter().map(|&n| Cell::Count(n)));
            cells.extend(
                row.exclusive_spectra_by_sample
                    .iter()
                    .map(|&n| Cell::Count(n)),
            );
            cells.extend(
                row.exclusive_peptides_by_sample
                    .iter()
                    .map(|&n| Cell::Count(n)),
            );
            let n = samples.len();
            cells.extend(
                (0..n).map(|s| Cell::Number(quant.map(|q| q.normalized[s]), Some(&formats.fixed2))),
            );
            cells.extend((0..n).map(|s| {
                Cell::Number(
                    quant.and_then(|q| q.nsaf.as_ref().map(|v| v[s])),
                    Some(&formats.scientific),
                )
            }));
            cells.extend((0..n).map(|s| {
                Cell::Number(
                    quant.and_then(|q| q.empai.as_ref().map(|v| v[s])),
                    Some(&formats.fixed2),
                )
            }));
            if ms1.is_some() {
                cells.extend((0..n).map(|s| {
                    Cell::Number(
                        lfq.and_then(|r| r.intensities[s]),
                        Some(&formats.scientific),
                    )
                }));
                cells.extend(
                    (0..n).map(|s| {
                        Cell::Number(lfq.and_then(|r| r.top3[s]), Some(&formats.scientific))
                    }),
                );
            }
            sheet.push(&cells).map_err(write_error)?;
        }
        let rows = sheet.finish().map_err(write_error)?;
        record("Proteins", list.rows.len(), rows);
    }

    // Peptides.
    let mut header: Vec<String> = [
        "Peptide",
        "Sequence",
        "Localized peptide",
        "Site probabilities",
        "Best Ascore",
        "q-value",
        "Probability",
        "Score",
        "Spectra",
    ]
    .map(String::from)
    .to_vec();
    header.extend(samples.iter().map(|s| format!("Spectra: {s}")));
    header.push("Protein groups".into());
    {
        let mut sheet =
            Sheet::new(&mut workbook, "Peptides", &header, &formats).map_err(write_error)?;
        sheet.sheet.set_column_width(0, 32).map_err(write_error)?;
        sheet.sheet.set_column_width(2, 32).map_err(write_error)?;
        for row in report.peptides.iter().take(MAX_ROWS) {
            let groups: Vec<String> = row
                .groups
                .iter()
                .map(|(accession, role)| format!("{accession} ({role})"))
                .collect();
            let mut cells = vec![
                Cell::Text(&row.proforma),
                Cell::Text(&row.sequence),
                Cell::Text(row.localized_peptide.as_deref().unwrap_or("")),
                Cell::Text(row.site_probabilities.as_deref().unwrap_or("")),
                Cell::Number(row.best_ascore, Some(&formats.fixed2)),
                Cell::Number(Some(row.q_value), Some(&formats.scientific)),
                Cell::Number(Some(row.probability), Some(&formats.percent)),
                Cell::Number(Some(row.score), Some(&formats.fixed4)),
                Cell::Count(row.spectra),
            ];
            cells.extend(row.spectra_by_sample.iter().map(|&n| Cell::Count(n)));
            cells.push(Cell::Owned(groups.join("; ")));
            sheet.push(&cells).map_err(write_error)?;
        }
        let rows = sheet.finish().map_err(write_error)?;
        record("Peptides", report.peptides.len(), rows);
    }

    // Spectra.
    let header: Vec<String> = [
        "Run",
        "Biosample",
        "Spectrum",
        "Charge",
        "Precursor m/z",
        "Retention time (min)",
        "Peptide",
        "Localized peptide",
        "Site probabilities",
        "Best Ascore",
        "q-value",
        "Probability",
        "Score",
        "Engines",
    ]
    .map(String::from)
    .to_vec();
    {
        let mut sheet =
            Sheet::new(&mut workbook, "Spectra", &header, &formats).map_err(write_error)?;
        sheet.sheet.set_column_width(6, 32).map_err(write_error)?;
        sheet.sheet.set_column_width(7, 32).map_err(write_error)?;
        for row in report.spectra.iter().take(MAX_ROWS) {
            sheet
                .push(&[
                    Cell::Text(&row.run),
                    Cell::Text(&row.biosample),
                    Cell::Text(&row.key),
                    Cell::Number(row.charge.map(|c| c as f64), None),
                    Cell::Number(row.precursor_mz, Some(&formats.fixed4)),
                    Cell::Number(row.rt_seconds.map(|s| s / 60.0), Some(&formats.fixed2)),
                    Cell::Text(&row.proforma),
                    Cell::Text(row.localized_peptide.as_deref().unwrap_or("")),
                    Cell::Text(row.site_probabilities.as_deref().unwrap_or("")),
                    Cell::Number(row.best_ascore, Some(&formats.fixed2)),
                    Cell::Number(Some(row.q_value), Some(&formats.scientific)),
                    Cell::Number(Some(row.probability), Some(&formats.percent)),
                    Cell::Number(Some(row.score), Some(&formats.fixed4)),
                    Cell::Count(row.engines as usize),
                ])
                .map_err(write_error)?;
        }
        let rows = sheet.finish().map_err(write_error)?;
        record("Spectra", report.spectra.len(), rows);
    }

    // ANOVA sheets, where two or more categories have data.
    let mut anovas: Vec<(&str, AnovaComparison)> = Vec::new();
    if let Ok(anova) = counts.anova() {
        anovas.push(("ANOVA (counts)", anova));
    }
    if let Some(Ok(anova)) = ms1.map(|q| q.quantities.anova()) {
        anovas.push(("ANOVA (MS1)", anova));
    }
    for (name, anova) in &anovas {
        let mut header = vec!["Accession".to_string()];
        header.extend(anova.categories.iter().map(|c| format!("Mean log2: {c}")));
        header.extend(anova.categories.iter().map(|c| format!("Measured: {c}")));
        header.extend(["log2 range", "F", "p-value", "q-value"].map(String::from));
        let mut sheet = Sheet::new(&mut workbook, name, &header, &formats).map_err(write_error)?;
        for row in anova.rows.iter().take(MAX_ROWS) {
            let mut cells = vec![Cell::Text(&row.accession)];
            cells.extend(
                row.means
                    .iter()
                    .map(|&m| Cell::Number(m, Some(&formats.fixed2))),
            );
            cells.extend(row.measured.iter().map(|&n| Cell::Count(n)));
            cells.extend([
                Cell::Number(row.range, Some(&formats.fixed2)),
                Cell::Number(row.f, Some(&formats.fixed2)),
                Cell::Number(row.p, Some(&formats.scientific)),
                Cell::Number(row.q, Some(&formats.scientific)),
            ]);
            sheet.push(&cells).map_err(write_error)?;
        }
        let rows = sheet.finish().map_err(write_error)?;
        record(name, anova.rows.len(), rows);
    }

    write_summary(project, &mut workbook, &formats, &written).map_err(|e| match e {
        SummaryError::Xlsx(error) => write_error(error),
        SummaryError::Store(error) => ExportError::Store(error),
    })?;
    workbook.save(path).map_err(write_error)?;
    Ok(written)
}

enum SummaryError {
    Xlsx(XlsxError),
    Store(StoreError),
}

impl From<XlsxError> for SummaryError {
    fn from(error: XlsxError) -> Self {
        Self::Xlsx(error)
    }
}

impl From<StoreError> for SummaryError {
    fn from(error: StoreError) -> Self {
        Self::Store(error)
    }
}

fn write_summary(
    project: &Project,
    workbook: &mut Workbook,
    formats: &Formats,
    written: &WorkbookReport,
) -> Result<(), SummaryError> {
    let summary = project.summary()?;
    let validation = project.validation_summary()?;
    let sheet = workbook.worksheet_from_name("Summary")?;
    sheet.set_column_width(0, 34)?;
    sheet.set_column_width(1, 90)?;
    let mut row = 0;
    let mut line = |sheet: &mut Worksheet, key: &str, value: String| -> Result<(), XlsxError> {
        sheet.write_string_with_format(row, 0, key, &formats.header)?;
        sheet.write_string(row, 1, value)?;
        row += 1;
        Ok(())
    };
    let synthetic = project
        .imports()?
        .iter()
        .any(|import| import.format == "synthetic");
    if synthetic {
        line(
            sheet,
            "SYNTHETIC DEMO DATA",
            "Nothing in this workbook was measured; it was generated to explore Tuffold.".into(),
        )?;
    }
    line(
        sheet,
        "Written by",
        format!("Tuffold {}", env!("CARGO_PKG_VERSION")),
    )?;
    line(
        sheet,
        "Project",
        format!(
            "{} biosamples, {} MS runs, {} spectra, {} PSMs",
            summary.biosamples, summary.ms_runs, summary.spectra, summary.psms
        ),
    )?;
    if let Some(validation) = &validation {
        let thresholds = &validation.settings.thresholds;
        line(
            sheet,
            "Thresholds",
            format!(
                "spectra {}, peptides {}, protein groups {}, at least {} peptide(s) per group",
                thresholds.describe(Level::Psm),
                thresholds.describe(Level::Peptide),
                thresholds.describe(Level::Protein),
                thresholds.min_peptides
            ),
        )?;
        for (name, counts) in [
            ("Spectra passing", validation.spectra),
            ("Peptides passing", validation.peptides),
            ("Protein groups passing", validation.protein_groups),
        ] {
            let fdr = counts.estimated_fdr().map_or_else(String::new, |f| {
                format!(", decoy-estimated FDR {}", tuffold_validate::format_rate(f))
            });
            line(
                sheet,
                name,
                format!("{} targets ({} decoys{fdr})", counts.targets, counts.decoys),
            )?;
        }
        let loc = &validation.localization;
        if loc.scored_psms > 0
            || loc.missing_peaks_psms > 0
            || loc.capped_psms > 0
            || loc.unscorable_psms > 0
        {
            let flr = loc.estimated_flr.map_or_else(String::new, |f| {
                format!(
                    ", decoy-AA estimated FLR {}",
                    tuffold_validate::format_rate(f)
                )
            });
            let mut value = format!(
                "{} of {} ambiguous passing PSMs ({} of {} sites) localized at >= {:.0}% site \
                 probability{flr}",
                loc.localized_psms,
                loc.scored_psms,
                loc.localized_sites,
                loc.total_target_sites,
                loc.probability_threshold * 100.0,
            );
            if loc.missing_peaks_psms > 0 {
                value.push_str(&format!("; {} without peaks", loc.missing_peaks_psms));
            }
            if loc.capped_psms > 0 {
                value.push_str(&format!("; {} over 256 isoforms", loc.capped_psms));
            }
            if loc.unscorable_psms > 0 {
                value.push_str(&format!("; {} unscorable", loc.unscorable_psms));
            }
            line(sheet, "PTM localization", value)?;
        }
        for engine in &validation.engines {
            line(
                sheet,
                &format!("Engine: {}", engine.name),
                format!(
                    "{} rank-1 PSMs ({} decoy), {}: {}",
                    engine.psms, engine.decoys, engine.method, engine.detail
                ),
            )?;
        }
    }
    for run in &summary.runs {
        line(
            sheet,
            &format!("Run: {}", run.name),
            format!(
                "biosample {}, {} spectra, {} PSMs",
                run.biosample, run.spectra, run.psms
            ),
        )?;
    }
    for (name, rows) in &written.sheets {
        line(sheet, &format!("Sheet: {name}"), format!("{rows} rows"))?;
    }
    for (name, total) in &written.cut {
        line(
            sheet,
            &format!("CUT: {name}"),
            format!(
                "the project has {total} rows; Excel holds {MAX_ROWS}, so the sheet has the \
                 first {MAX_ROWS}. Export the TSV tables for all of them."
            ),
        )?;
    }
    Ok(())
}

/// Write a tab-separated table of every passing target spectrum with
/// ambiguous modification localization (`report().spectra` rows where
/// `localized_peptide` is present). Returns the number of data rows written.
pub fn write_localization_tsv(
    out: &mut impl std::io::Write,
    project: &Project,
) -> Result<usize, ExportError> {
    let report = project.report()?;
    let write_err = |e: std::io::Error| ExportError::Write {
        path: "tsv".into(),
        message: e.to_string(),
    };
    writeln!(
        out,
        "run\tbiosample\tspectrum\tcharge\tprecursor_mz\trt_seconds\tpeptide\tlocalized_peptide\tsite_probabilities\tbest_ascore\tq_value\tprobability"
    )
    .map_err(write_err)?;
    let mut written = 0;
    for row in &report.spectra {
        let Some(localized) = row.localized_peptide.as_deref() else {
            continue;
        };
        let charge = row.charge.map_or_else(String::new, |c| c.to_string());
        let mz = row
            .precursor_mz
            .map_or_else(String::new, |m| format!("{m:.4}"));
        let rt = row
            .rt_seconds
            .map_or_else(String::new, |s| format!("{s:.2}"));
        let probs = row.site_probabilities.as_deref().unwrap_or("");
        let ascore = row
            .best_ascore
            .map_or_else(String::new, |a| format!("{a:.2}"));
        writeln!(
            out,
            "{}\t{}\t{}\t{charge}\t{mz}\t{rt}\t{}\t{localized}\t{probs}\t{ascore}\t{:.6e}\t{:.4}",
            row.run, row.biosample, row.key, row.proforma, row.q_value, row.probability,
        )
        .map_err(write_err)?;
        written += 1;
    }
    Ok(written)
}
