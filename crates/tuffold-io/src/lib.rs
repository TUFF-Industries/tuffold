//! Tuffold importers.
//!
//! - [`fasta`]: protein-database parsing, identifier indexing, and exact-match
//!   peptide coverage.
//! - [`fragpipe`]: FragPipe/Philosopher results bundles — read-only
//!   [`inspect_bundle`] diagnostics and the [`import_bundle`] mapping.
//! - [`go`]: Gene Ontology terms, slims, and annotations (OBO, GAF).
//! - [`mascot`]: Mascot results files (`.dat`).
//! - [`maxquant`]: MaxQuant `txt` folders (`msms.txt`).
//! - [`pd`]: Proteome Discoverer 2.x/3.x results (`.pdResult`, `.msf`).
//! - [`sf3`]: Scaffold `.sf3` projects (read-only), one batch per engine.
//! - [`sfdb`]: Scaffold `.sfdb` files (Scaffold 5, Scaffold DDA; read-only).
//! - [`diagnostics`]: the typed issue vocabulary shared by every importer.

pub mod diagnostics;
pub mod fasta;
pub mod fragpipe;
pub mod go;
pub mod mascot;
pub mod maxquant;
pub mod mzid;
pub mod pd;
pub mod pepxml;
pub mod sf3;
pub mod sfdb;
mod source;
pub mod spectra;
mod tsv;
mod xml;

pub use diagnostics::{InspectionIssue, IssueCode, Severity};
pub use fragpipe::{
    BundleInspection, ParseError, SearchFile, find_search_files, import_bundle,
    import_search_files, inspect_bundle,
};
pub use mascot::{MascotError, MascotImport, import_mascot};
pub use maxquant::{MaxQuantError, MaxQuantImport, import_maxquant};
pub use mzid::import_mzid;
pub use pd::{PdError, PdImport, import_pd};
pub use pepxml::import_pepxml;
pub use sf3::{Sf3Error, Sf3Import, Sf3Peaks, UnreadablePeaks, import_sf3};
pub use sfdb::{SfdbError, SfdbImport, SfdbPeaks, import_sfdb};
pub use source::{ImportError, Imported, StoredPeaks, detect_format, import_path};
pub use spectra::{
    SpectraError, THERMO_RAW, is_spectra_file, is_thermo_raw, read_scans, read_spectra,
    spectra_run_name,
};
pub use xml::XmlImportError;
