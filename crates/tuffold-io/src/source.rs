//! One entry point for every supported input: detect the format, import it.

use std::path::Path;

use tuffold_core::{ImportBatch, PeakSpectrum, SourceFormat};

use crate::fragpipe;

#[derive(Debug, thiserror::Error)]
pub enum ImportError {
    #[error(
        "{0}: unrecognized input; expected a FragPipe results folder, a .pep.xml/.pepXML file, \
         an .mzid file, a Mascot .dat file, a MaxQuant txt folder, a Proteome Discoverer \
         .pdResult/.msf file, a Scaffold .sf3 project, or a Scaffold .sfdb file"
    )]
    UnknownFormat(String),
    #[error("{0}: {raw}", raw = crate::spectra::THERMO_RAW)]
    ThermoRaw(String),
    #[error(transparent)]
    FragPipe(#[from] fragpipe::ParseError),
    #[error(transparent)]
    Xml(#[from] crate::xml::XmlImportError),
    #[error(transparent)]
    Mascot(#[from] crate::mascot::MascotError),
    #[error(transparent)]
    MaxQuant(#[from] crate::maxquant::MaxQuantError),
    #[error("{file}: {source}")]
    Pd {
        file: String,
        #[source]
        source: crate::pd::PdError,
    },
    #[error("{file}: {source}")]
    Sf3 {
        file: String,
        #[source]
        source: crate::sf3::Sf3Error,
    },
    #[error("{file}: {source}")]
    Sfdb {
        file: String,
        #[source]
        source: crate::sfdb::SfdbError,
    },
}

/// Peaks a Scaffold file stores for its identified spectra.
#[derive(Debug, Clone, PartialEq)]
pub enum StoredPeaks {
    Sf3(crate::Sf3Peaks),
    Sfdb(crate::SfdbPeaks),
}

impl StoredPeaks {
    /// The runs holding spectra with peaks, in first-seen order.
    pub fn runs(&self) -> Vec<&str> {
        match self {
            Self::Sf3(peaks) => peaks.runs(),
            Self::Sfdb(peaks) => peaks.runs(),
        }
    }

    /// Identified spectra (any run) with peaks stored.
    pub fn len(&self) -> usize {
        match self {
            Self::Sf3(peaks) => peaks.len(),
            Self::Sfdb(peaks) => peaks.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The run's spectra with their peaks, sorted by m/z.
    pub fn spectra<'a>(
        &'a self,
        run: &'a str,
    ) -> Box<dyn Iterator<Item = Result<PeakSpectrum, crate::UnreadablePeaks>> + 'a> {
        match self {
            Self::Sf3(peaks) => Box::new(peaks.spectra(run)),
            Self::Sfdb(peaks) => Box::new(peaks.spectra(run).into_iter()),
        }
    }
}

/// What one input produced: usually one batch; a Scaffold project gives one
/// per search engine, with notes on what was left out.
#[derive(Debug, Clone, PartialEq)]
pub struct Imported {
    pub format: SourceFormat,
    pub batches: Vec<ImportBatch>,
    pub notes: Vec<String>,
    /// Peaks the source itself stores (Scaffold files), for the
    /// identified spectra; attached after the batches are imported.
    pub peaks: Option<StoredPeaks>,
}

/// Decide the format from the path alone: folders are FragPipe bundles, files
/// go by extension (case-insensitive).
pub fn detect_format(path: &Path) -> Option<SourceFormat> {
    if crate::maxquant::is_maxquant_folder(path) {
        return Some(SourceFormat::MaxQuant);
    }
    if path.is_dir() {
        return Some(SourceFormat::FragPipeBundle);
    }
    let name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    if name == "msms.txt" || name == "evidence.txt" {
        Some(SourceFormat::MaxQuant)
    } else if name.ends_with(".pep.xml") || name.ends_with(".pepxml") {
        Some(SourceFormat::PepXml)
    } else if name.ends_with(".mzid") {
        Some(SourceFormat::MzIdentMl)
    } else if name.ends_with(".dat") {
        Some(SourceFormat::Mascot)
    } else if name.ends_with(".sf3") {
        Some(SourceFormat::Sf3)
    } else if name.ends_with(".sfdb") {
        Some(SourceFormat::Sfdb)
    } else if name.ends_with(".pdresult") || name.ends_with(".msf") {
        Some(SourceFormat::ProteomeDiscoverer)
    } else {
        None
    }
}

/// Import any supported input.
pub fn import_path(path: &Path) -> Result<Imported, ImportError> {
    let format = detect_format(path).ok_or_else(|| {
        let shown = path.display().to_string();
        if crate::spectra::is_thermo_raw(path) {
            ImportError::ThermoRaw(shown)
        } else {
            ImportError::UnknownFormat(shown)
        }
    })?;
    let one = |batch: ImportBatch| Imported {
        format,
        batches: vec![batch],
        notes: Vec::new(),
        peaks: None,
    };
    Ok(match format {
        SourceFormat::FragPipeBundle => {
            let searches = fragpipe::find_search_files(path)?;
            if searches.is_empty() {
                one(fragpipe::import_bundle(path)?)
            } else {
                let batch = fragpipe::import_search_files(path, &searches)?;
                let note = format!(
                    "read {}'s own results ({} pepXML file{}, every match and decoy); \
                     FragPipe's psm.tsv and its PeptideProphet/Percolator verdicts were \
                     not used",
                    batch.software.name,
                    searches.len(),
                    if searches.len() == 1 { "" } else { "s" }
                );
                Imported {
                    format,
                    batches: vec![batch],
                    notes: vec![note],
                    peaks: None,
                }
            }
        }
        SourceFormat::PepXml => one(crate::pepxml::import_pepxml(path)?),
        SourceFormat::MzIdentMl => one(crate::mzid::import_mzid(path)?),
        SourceFormat::Sf3 => {
            let imported = crate::sf3::import_sf3(path).map_err(|source| ImportError::Sf3 {
                file: path.display().to_string(),
                source,
            })?;
            Imported {
                format,
                batches: imported.batches,
                notes: imported.notes,
                peaks: imported.peaks.map(StoredPeaks::Sf3),
            }
        }
        SourceFormat::Sfdb => {
            let imported = crate::sfdb::import_sfdb(path).map_err(|source| ImportError::Sfdb {
                file: path.display().to_string(),
                source,
            })?;
            Imported {
                format,
                batches: imported.batches,
                notes: imported.notes,
                peaks: imported.peaks.map(StoredPeaks::Sfdb),
            }
        }
        SourceFormat::Mascot => {
            let imported = crate::mascot::import_mascot(path)?;
            Imported {
                format,
                batches: vec![imported.batch],
                notes: imported.notes,
                peaks: None,
            }
        }
        SourceFormat::MaxQuant => {
            let imported = crate::maxquant::import_maxquant(path)?;
            Imported {
                format,
                batches: vec![imported.batch],
                notes: imported.notes,
                peaks: None,
            }
        }
        SourceFormat::ProteomeDiscoverer => {
            let imported = crate::pd::import_pd(path).map_err(|source| ImportError::Pd {
                file: path.display().to_string(),
                source,
            })?;
            Imported {
                format,
                batches: imported.batches,
                notes: imported.notes,
                peaks: None,
            }
        }
        SourceFormat::Synthetic | SourceFormat::SageSearch => {
            return Err(ImportError::UnknownFormat(path.display().to_string()));
        }
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn thermo_raw_files_are_refused_with_how_to_convert_them() {
        let dir = tempfile::tempdir().unwrap();
        let raw = dir.path().join("run01.RAW");
        std::fs::write(&raw, b"not read").unwrap();
        let import = import_path(&raw).err().unwrap();
        assert!(matches!(import, ImportError::ThermoRaw(_)), "{import}");
        assert!(
            import.to_string().contains("ThermoRawFileParser"),
            "{import}"
        );
        let spectra = crate::spectra::read_spectra(&raw).err().unwrap();
        assert!(matches!(
            spectra,
            crate::spectra::SpectraError::ThermoRaw(_)
        ));
        // Any other unknown file keeps the general message.
        let other = dir.path().join("notes.txt");
        std::fs::write(&other, b"x").unwrap();
        assert!(matches!(
            import_path(&other),
            Err(ImportError::UnknownFormat(_))
        ));
    }

    use super::*;

    #[test]
    fn formats_are_detected_from_the_path() {
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(
            detect_format(dir.path()),
            Some(SourceFormat::FragPipeBundle)
        );
        for (name, format) in [
            ("run.pep.xml", Some(SourceFormat::PepXml)),
            ("RUN.pepXML", Some(SourceFormat::PepXml)),
            ("interact.pep.xml.gz", Some(SourceFormat::PepXml)),
            ("search.mzid", Some(SourceFormat::MzIdentMl)),
            ("search.MZID.gz", Some(SourceFormat::MzIdentMl)),
            ("project.sf3", Some(SourceFormat::Sf3)),
            ("P706_090622.SFDB", Some(SourceFormat::Sfdb)),
            ("Study.pdResult", Some(SourceFormat::ProteomeDiscoverer)),
            ("search.msf", Some(SourceFormat::ProteomeDiscoverer)),
            ("F012345.dat", Some(SourceFormat::Mascot)),
            ("notes.txt", None),
        ] {
            assert_eq!(detect_format(&dir.path().join(name)), format, "{name}");
        }
    }
}
