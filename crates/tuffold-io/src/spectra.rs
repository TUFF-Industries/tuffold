//! Measured spectra from mzML and MGF files (via `mzdata`) and Bruker
//! timsTOF `.d` folders (via `timsrust`).
//!
//! [`read_spectra`] yields only MS/MS (level ≥ 2) spectra: identifications
//! never point at MS1 scans. [`read_scans`] also yields mzML's MS1 scans
//! (level 1), in file order, for precursor intensity. Each spectrum gets the same [`SpectrumKey`] the identification
//! importers derive, so peaks attach to the right PSMs:
//! - mzML: from the native id (`... scan=N`, `index=N`, else verbatim). An
//!   MS1 scan holding peaks no instrument could have measured is yielded as
//!   [`SpectraError::DamagedMs1`] instead (see `impossible_peaks`).
//! - MGF: `SCANS=`, else a `scan=N` or TPP-style `run.start.end.z` title,
//!   else the 0-based position in the file (`index=N`, as mzIdentML uses).
//! - timsTOF `.d` (DDA-PASEF): one MS/MS spectrum per precursor, its PASEF
//!   frames summed by timsrust (as Sage and MSFragger search them), keyed by
//!   the precursor's id in the `Precursors` table (`Scan(id)`, native id
//!   `precursor=id`), with the precursor's ion mobility (1/K0). [`read_scans`]
//!   also yields each MS1 frame binned by mobility (`tdf::bin`), every peak
//!   with its bin's 1/K0, so the MS1 tracer can measure a precursor at its
//!   own mobility.
//!
//! The run name is the file name without its spectra extension, matching
//! how the identification importers name runs. A retention time of exactly
//! zero is treated as absent (mzdata reports missing times as 0).

use std::path::Path;

use mzdata::io::{MGFReader, MzMLReader};
use mzdata::prelude::*;
use mzdata::spectrum::{MultiLayerSpectrum, RefPeakDataLevel, SignalContinuity};
use mzdata::utils::mass_charge_ratio;
use tuffold_core::{PeakSpectrum, SpectrumKey};

/// What a Thermo `.raw` file needs before Tuffold can read it.
pub const THERMO_RAW: &str = "a Thermo .raw file, which Tuffold doesn't read directly; \
     convert it to mzML first (ThermoRawFileParser `-f 2`, or ProteoWizard msconvert) and add \
     the .mzML (see the user guide, Thermo .raw files)";

/// A file named `*.raw` (any case): Thermo's format, not read here.
pub fn is_thermo_raw(path: &Path) -> bool {
    path.is_file()
        && path
            .extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("raw"))
}

#[derive(Debug, thiserror::Error)]
pub enum SpectraError {
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: std::io::Error,
    },
    #[error("{file}: spectrum '{id}': {problem}")]
    Invalid {
        file: String,
        id: String,
        problem: String,
    },
    #[error("{0}: not a spectra file; expected .mzML, .mgf, or a Bruker .d folder")]
    UnknownFormat(String),
    #[error("{0}: {raw}", raw = THERMO_RAW)]
    ThermoRaw(String),
    #[error("{file}: {problem}")]
    Tdf { file: String, problem: String },
    /// An MS1 scan that can't be trusted: a timsTOF frame whose stored data
    /// can't be decoded, or an mzML scan holding peaks no instrument could
    /// have measured (see [`impossible_peaks`]). Readers of [`read_scans`]
    /// may leave it out (counted) and go on: the scans around it are intact,
    /// and MS/MS spectra are never skipped this way.
    #[error("{file}: MS1 {scan} is damaged ({problem})")]
    DamagedMs1 {
        file: String,
        /// Its native id (`scan=1914`, `frame=7`).
        scan: String,
        problem: String,
    },
    #[error(
        "{file}: only {read} of the {declared} spectra the file declares could be read; \
         it is truncated or corrupt from spectrum {read} on"
    )]
    Truncated {
        file: String,
        read: u64,
        declared: u64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    MzMl,
    Mgf,
    /// A Bruker timsTOF `.d` folder holding `analysis.tdf`.
    Tdf,
}

fn format_of(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if name.ends_with(".mzml") {
        Some(Format::MzMl)
    } else if name.ends_with(".mgf") {
        Some(Format::Mgf)
    } else if name.ends_with(".d") && path.join("analysis.tdf").is_file() {
        Some(Format::Tdf)
    } else {
        None
    }
}

/// Whether the path names a spectra file this module reads.
pub fn is_spectra_file(path: &Path) -> bool {
    format_of(path).is_some()
}

/// `HeLa_01.mzML` -> `HeLa_01`.
pub fn spectra_run_name(path: &Path) -> String {
    let name = crate::tsv::file_name(path);
    let lowered = name.to_ascii_lowercase();
    [".mzml", ".mgf", ".d"]
        .iter()
        .find(|extension| lowered.ends_with(*extension))
        .map_or(name.clone(), |extension| {
            name[..name.len() - extension.len()].to_string()
        })
}

type Spectra = Box<dyn Iterator<Item = Result<PeakSpectrum, SpectraError>>>;

/// Every MS/MS spectrum in the file, in file order.
pub fn read_spectra(path: &Path) -> Result<Spectra, SpectraError> {
    open(path, false)
}

/// Every spectrum in the file, MS1 scans included (`ms_level` 1), in file
/// order. MGF holds no MS1 scans, so it yields the same as [`read_spectra`].
/// A timsTOF `.d` yields each MS1 frame binned by ion mobility and
/// centroided (peaks carry their 1/K0), merged with its MS/MS spectra in
/// retention-time order.
pub fn read_scans(path: &Path) -> Result<Spectra, SpectraError> {
    open(path, true)
}

fn open(path: &Path, with_ms1: bool) -> Result<Spectra, SpectraError> {
    let format = format_of(path).ok_or_else(|| {
        let shown = path.display().to_string();
        if is_thermo_raw(path) {
            SpectraError::ThermoRaw(shown)
        } else {
            SpectraError::UnknownFormat(shown)
        }
    })?;
    let file_label = crate::tsv::file_name(path);
    if format == Format::Tdf {
        let spectra = tdf::TdfSpectra::open(path, file_label.clone())?;
        return Ok(if with_ms1 {
            Box::new(tdf::TdfScans::open(path, file_label, spectra)?)
        } else {
            Box::new(spectra)
        });
    }
    let io_error = |source| SpectraError::Io {
        file: file_label.clone(),
        source,
    };
    // mzdata stops quietly at a parse error, so a truncated mzML would look
    // complete; the spectrumList count it records lets us notice. (MGF has
    // no declared count.)
    let (spectra, declared): (Box<dyn Iterator<Item = MultiLayerSpectrum>>, Option<u64>) =
        match format {
            Format::MzMl => {
                let reader = MzMLReader::open_path(path).map_err(io_error)?;
                let declared = reader.spectrum_count_hint();
                (Box::new(reader), declared)
            }
            Format::Mgf => (
                Box::new(MGFReader::open_path(path).map_err(io_error)?),
                None,
            ),
            Format::Tdf => unreachable!("read above"),
        };
    Ok(Box::new(Checked {
        spectra,
        declared,
        read: 0,
        file: file_label,
        format,
        with_ms1,
        done: false,
    }))
}

/// Yields converted MS/MS spectra, then an error if fewer spectra (of any
/// level) were read than the file declares.
struct Checked {
    spectra: Box<dyn Iterator<Item = MultiLayerSpectrum>>,
    declared: Option<u64>,
    read: u64,
    file: String,
    format: Format,
    with_ms1: bool,
    done: bool,
}

impl Iterator for Checked {
    type Item = Result<PeakSpectrum, SpectraError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        for spectrum in self.spectra.by_ref() {
            self.read += 1;
            if self.with_ms1 || spectrum.ms_level() >= 2 || self.format == Format::Mgf {
                return Some(convert(&spectrum, self.format, &self.file));
            }
        }
        self.done = true;
        match self.declared {
            Some(declared) if declared != self.read => Some(Err(SpectraError::Truncated {
                file: self.file.clone(),
                read: self.read,
                declared,
            })),
            _ => None,
        }
    }
}

fn convert(
    spectrum: &MultiLayerSpectrum,
    format: Format,
    file: &str,
) -> Result<PeakSpectrum, SpectraError> {
    let native_id = spectrum.id().to_string();
    let invalid = |problem: String| SpectraError::Invalid {
        file: file.to_string(),
        id: native_id.clone(),
        problem,
    };
    let (mz, intensity) = match spectrum.peaks() {
        RefPeakDataLevel::Missing => (Vec::new(), Vec::new()),
        RefPeakDataLevel::RawData(arrays) => {
            let mz = arrays
                .mzs()
                .map_err(|error| invalid(format!("m/z array: {error}")))?;
            let intensity = arrays
                .intensities()
                .map_err(|error| invalid(format!("intensity array: {error}")))?;
            if mz.len() != intensity.len() {
                return Err(invalid(format!(
                    "{} m/z values but {} intensities",
                    mz.len(),
                    intensity.len()
                )));
            }
            (mz.into_owned(), intensity.into_owned())
        }
        RefPeakDataLevel::Centroid(peaks) => {
            peaks.iter().map(|peak| (peak.mz, peak.intensity)).unzip()
        }
        // mzdata reads an MGF peak line with a charge column ("345.1 237 3") as
        // a deconvoluted peak, converting m/z to neutral mass; MGF still means
        // m/z, so undo that conversion exactly.
        // mzdata also orders those peaks by neutral mass; restore m/z order.
        RefPeakDataLevel::Deconvoluted(peaks) if format == Format::Mgf => {
            let mut pairs: Vec<(f64, f32)> = peaks
                .iter()
                .map(|peak| {
                    (
                        mass_charge_ratio(peak.neutral_mass, peak.charge),
                        peak.intensity,
                    )
                })
                .collect();
            pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
            pairs.into_iter().unzip()
        }
        RefPeakDataLevel::Deconvoluted(_) => {
            return Err(invalid(
                "holds deconvoluted (neutral mass) peaks, not m/z peaks".into(),
            ));
        }
    };
    // Some mzML converters have written scans whose arrays hold garbage in
    // places; such a scan can't be trusted in part, so an MS1 scan is left
    // out whole. (MS/MS spectra are read as they are.)
    if format == Format::MzMl
        && spectrum.ms_level() == 1
        && let Some(problem) = impossible_peaks(&mz, &intensity, scan_window(spectrum))
    {
        return Err(SpectraError::DamagedMs1 {
            file: file.to_string(),
            scan: native_id,
            problem,
        });
    }
    // The first precursor that names an ion: when an MGF block gives CHARGE
    // before PEPMASS, mzdata records an empty placeholder precursor first.
    let ion = spectrum
        .description()
        .precursor
        .iter()
        .find_map(|precursor| precursor.ion());
    // The scan it was acquired from: ThermoRawFileParser gives every
    // precursor of an MSⁿ scan the n−1 scan's id ("Master Scan Number"), so
    // the first one's is enough.
    let parent = spectrum
        .description()
        .precursor
        .iter()
        .find_map(|precursor| precursor.precursor_id.clone())
        .filter(|id| !id.is_empty());
    // The isolation window of that precursor, when the file gives one.
    let isolation = spectrum
        .description()
        .precursor
        .iter()
        .find(|precursor| precursor.ion().is_some())
        .map(|precursor| &precursor.isolation_window)
        .map(|w| (f64::from(w.lower_bound), f64::from(w.upper_bound)))
        .filter(|&(low, high)| low > 0.0 && high > low);
    let minutes = spectrum.start_time();
    Ok(PeakSpectrum {
        key: match format {
            Format::MzMl | Format::Tdf => SpectrumKey::from_native_id(&native_id),
            Format::Mgf => mgf_key(spectrum, &native_id),
        },
        // MGF has no levels (mzdata may say 1); everything in it is MS/MS.
        ms_level: match format {
            Format::MzMl | Format::Tdf => spectrum.ms_level(),
            Format::Mgf => spectrum.ms_level().max(2),
        },
        rt_seconds: (minutes > 0.0).then_some(minutes * 60.0),
        precursor_mz: ion.map(|ion| ion.mz),
        precursor_charge: ion
            .and_then(|ion| ion.charge)
            .and_then(|charge| i16::try_from(charge).ok()),
        mz,
        intensity,
        profile: spectrum.signal_continuity() == SignalContinuity::Profile,
        native_id,
        parent,
        isolation,
        // timsTOF mobility comes from `.d` folders; mzML converted from them
        // isn't read for it.
        precursor_mobility: None,
        mobility: Vec::new(),
    })
}

/// No ion has an m/z below a proton's (1.007).
const MIN_MZ: f64 = 1.0;
/// Above any analyzer's range; the ceiling for a scan that declares no
/// window.
const MAX_MZ: f64 = 1e6;
/// How far outside its declared scan window a peak may lie. Calibration
/// puts Orbitrap peaks within a few mTh of it (at most 2.1 mTh below the
/// lower limit in six runs, none above); 5 Th also leaves room for profile
/// grids that run past a limit (not checked on an ion trap).
const WINDOW_SLACK: f64 = 5.0;
/// No detector reports a peak this intense (the strongest MS1 peak of six
/// Q Exactive HF runs is 7e9; f32 reaches 3e38, where garbage lands).
const MAX_INTENSITY: f32 = 1e15;

/// The m/z range a spectrum declares it scanned (the union of its scan
/// windows), when it declares a usable one.
fn scan_window(spectrum: &MultiLayerSpectrum) -> Option<(f64, f64)> {
    let (low, high) = spectrum
        .description()
        .acquisition
        .scans
        .iter()
        .flat_map(|scan| scan.scan_windows.iter())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), w| {
            (
                low.min(f64::from(w.lower_bound)),
                high.max(f64::from(w.upper_bound)),
            )
        });
    (low > 0.0 && high > low && high.is_finite()).then_some((low, high))
}

/// Whether an MS1 scan's peaks include ones no instrument could have
/// measured, described if so: an m/z that isn't a number, is below a
/// proton's, or lies more than [`WINDOW_SLACK`] outside the scan's window
/// (without one, above [`MAX_MZ`]); an intensity that isn't a number, is
/// negative, or exceeds [`MAX_INTENSITY`]. Normal data can't trip it. A scan
/// that does has arrays that are garbage in places (ThermoRawFileParser
/// 2.0.0-dev wrote some), and the good-looking peaks beside the garbage are
/// no more trustworthy: in two of the three MS1 scans seen, the peaks that
/// pass sum to 19 and 3.9 times the scan's own total ion current.
fn impossible_peaks(mz: &[f64], intensity: &[f32], window: Option<(f64, f64)>) -> Option<String> {
    let (low, high) = window.map_or((MIN_MZ, MAX_MZ), |(low, high)| {
        ((low - WINDOW_SLACK).max(MIN_MZ), high + WINDOW_SLACK)
    });
    let impossible = |&(_, (&m, &y)): &(usize, (&f64, &f32))| {
        let mz_possible = m.is_finite() && (low..=high).contains(&m);
        let intensity_possible = y.is_finite() && (0.0..=MAX_INTENSITY).contains(&y);
        !(mz_possible && intensity_possible)
    };
    let mut peaks = mz.iter().zip(intensity).enumerate().filter(impossible);
    let (first, (&m, &y)) = peaks.next()?;
    Some(format!(
        "{} of {} peaks are impossible, the first (peak {first}) at m/z {m:e} with intensity \
         {y:e}; m/z must lie within {low}-{high} and intensity within 0-{MAX_INTENSITY:e}",
        peaks.count() + 1,
        mz.len(),
    ))
}

fn mgf_key(spectrum: &MultiLayerSpectrum, title: &str) -> SpectrumKey {
    let scans = spectrum
        .description()
        .get_param_by_name("scans") // mzdata lower-cases unrecognized MGF headers
        .and_then(|param| {
            param
                .value
                .to_string()
                .split(['-', ','])
                .next()?
                .trim()
                .parse()
                .ok()
        });
    if let Some(scan) = scans {
        return SpectrumKey::Scan(scan);
    }
    if let SpectrumKey::Scan(scan) = SpectrumKey::from_native_id(title) {
        return SpectrumKey::Scan(scan);
    }
    let mut parts = title.trim().rsplitn(4, '.');
    let (charge, end, start) = (parts.next(), parts.next(), parts.next());
    let is_tpp = charge.is_some_and(|c| c.parse::<u8>().is_ok())
        && end.is_some_and(|e| e.parse::<u32>().is_ok());
    if let (true, Some(Ok(scan))) = (is_tpp, start.map(str::parse::<u32>)) {
        return SpectrumKey::Scan(scan);
    }
    SpectrumKey::Index(u32::try_from(spectrum.index()).unwrap_or(u32::MAX))
}

mod tdf {
    use std::iter::Peekable;
    use std::path::Path;

    use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};
    use timsrust::readers::{FrameReader, MetadataReader, SpectrumReader};
    use timsrust::{Frame, MSLevel, Spectrum};
    use tuffold_core::centroid::centroid;
    use tuffold_core::{PeakSpectrum, SpectrumKey};

    use super::SpectraError;

    /// Width (1/K0, V·s/cm²) of the mobility bins MS1 frames are summed in:
    /// about a quarter of a peptide's mobility peak on a timsTOF Pro (FWHM
    /// ~0.02), fine enough for the tracer's window around a precursor, and
    /// ~6× the collapsed frame's peaks rather than every scan's.
    const MOBILITY_BIN: f64 = 0.005;

    /// Runs a timsrust call, turning a panic on malformed data into an error
    /// (as the search does for Sage), so a bad file fails or skips only
    /// what it must instead of taking the program down.
    fn guarded<T, E: std::fmt::Display>(call: impl FnOnce() -> Result<T, E>) -> Result<T, String> {
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(call)) {
            Ok(result) => result.map_err(|error| error.to_string()),
            Err(_) => Err("the timsTOF reader failed on this data".into()),
        }
    }

    /// A `.d` folder's MS/MS spectra, one per PASEF precursor, in precursor
    /// order.
    pub(super) struct TdfSpectra {
        reader: SpectrumReader,
        next: usize,
        file: String,
    }

    impl TdfSpectra {
        pub(super) fn open(path: &Path, file: String) -> Result<Self, SpectraError> {
            let reader =
                guarded(|| SpectrumReader::new(path)).map_err(|error| SpectraError::Tdf {
                    file: file.clone(),
                    problem: format!("cannot read the timsTOF data: {error}"),
                })?;
            Ok(Self {
                reader,
                next: 0,
                file,
            })
        }
    }

    impl Iterator for TdfSpectra {
        type Item = Result<PeakSpectrum, SpectraError>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.next >= self.reader.len() {
                return None;
            }
            let index = self.next;
            self.next += 1;
            Some(
                guarded(|| self.reader.get(index))
                    .map(|spectrum| convert(&spectrum, index))
                    .map_err(|problem| SpectraError::Invalid {
                        file: self.file.clone(),
                        id: format!("spectrum {index}"),
                        problem,
                    }),
            )
        }
    }

    /// One summed PASEF spectrum. Without a precursor (not seen in DDA
    /// data) it is keyed by its position, as MGF spectra without scans are.
    pub(super) fn convert(spectrum: &Spectrum, index: usize) -> PeakSpectrum {
        let precursor = spectrum.precursor.as_ref();
        let (key, native_id) = match precursor.and_then(|p| u32::try_from(p.index).ok()) {
            Some(id) => (SpectrumKey::Scan(id), format!("precursor={id}")),
            None => {
                let index = u32::try_from(index).unwrap_or(u32::MAX);
                (SpectrumKey::Index(index), format!("index={index}"))
            }
        };
        PeakSpectrum {
            key,
            native_id,
            ms_level: 2,
            rt_seconds: precursor
                .map(|p| p.rt)
                .filter(|rt| rt.is_finite() && *rt > 0.0),
            precursor_mz: precursor.map(|p| p.mz).filter(|mz| *mz > 0.0),
            precursor_charge: precursor
                .and_then(|p| p.charge)
                .and_then(|charge| i16::try_from(charge).ok())
                .filter(|&charge| charge > 0),
            precursor_mobility: precursor
                .map(|p| p.im)
                .filter(|im| im.is_finite() && *im > 0.0),
            mz: spectrum.mz_values.clone(),
            // timsTOF intensities are detector counts; f32 holds them exactly
            // to 2^24.
            #[allow(clippy::cast_possible_truncation)]
            intensity: spectrum.intensities.iter().map(|&i| i as f32).collect(),
            profile: false,
            parent: None,
            isolation: None,
            mobility: Vec::new(),
        }
    }

    /// A `.d` folder's MS1 frames and MS/MS spectra, in retention-time order,
    /// so the MS1 tracer sees each identification among the frames around it.
    pub(super) struct TdfScans {
        frames: FrameReader,
        mz: Tof2MzConverter,
        im: Scan2ImConverter,
        /// Mobility scans summed per bin; `None` when the run has no
        /// mobility range, and frames are summed whole.
        scans_per_bin: Option<usize>,
        /// Positions (in the frame reader) of the MS1 frames, in time order.
        ms1: Vec<usize>,
        next_ms1: usize,
        spectra: Peekable<TdfSpectra>,
        file: String,
    }

    impl TdfScans {
        pub(super) fn open(
            path: &Path,
            file: String,
            spectra: TdfSpectra,
        ) -> Result<Self, SpectraError> {
            let unreadable = |error: String| SpectraError::Tdf {
                file: file.clone(),
                problem: format!("cannot read the timsTOF MS1 frames: {error}"),
            };
            let frames = guarded(|| FrameReader::new(path)).map_err(unreadable)?;
            let metadata = guarded(|| MetadataReader::new(path)).map_err(unreadable)?;
            let per_scan =
                (metadata.im_converter.convert(1_u32) - metadata.im_converter.convert(0_u32)).abs();
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let scans_per_bin = (per_scan.is_finite() && per_scan > 0.0)
                .then(|| ((MOBILITY_BIN / per_scan).round() as usize).max(1));
            let mut ms1: Vec<(f64, usize)> = (0..frames.len())
                .filter_map(|i| frames.get_frame_without_coordinates(i).ok().map(|f| (i, f)))
                .filter(|(_, f)| f.ms_level == MSLevel::MS1)
                .map(|(i, f)| (f.rt_in_seconds, i))
                .collect();
            ms1.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            Ok(Self {
                frames,
                mz: metadata.mz_converter,
                im: metadata.im_converter,
                scans_per_bin,
                ms1: ms1.into_iter().map(|(_, i)| i).collect(),
                next_ms1: 0,
                spectra: spectra.peekable(),
                file,
            })
        }

        fn ms1_rt(&self) -> Option<f64> {
            let &position = self.ms1.get(self.next_ms1)?;
            self.frames
                .get_frame_without_coordinates(position)
                .ok()
                .map(|f| f.rt_in_seconds)
        }

        fn read_ms1(&mut self) -> Result<PeakSpectrum, SpectraError> {
            let position = self.ms1[self.next_ms1];
            self.next_ms1 += 1;
            let frame = guarded(|| self.frames.get(position)).map_err(|problem| {
                // Bruker numbers frames from 1, as `Frames.Id` does.
                let frame = self
                    .frames
                    .get_frame_without_coordinates(position)
                    .map_or(position + 1, |f| f.index);
                SpectraError::DamagedMs1 {
                    file: self.file.clone(),
                    scan: format!("frame={frame}"),
                    problem,
                }
            })?;
            let (mz, intensity, mobility) = match self.scans_per_bin {
                Some(scans) => bin(&frame, &self.mz, &self.im, scans),
                None => {
                    let (mz, intensity) = collapse(&frame, &self.mz);
                    (mz, intensity, Vec::new())
                }
            };
            Ok(PeakSpectrum {
                key: SpectrumKey::Native(format!("frame={}", frame.index)),
                native_id: format!("frame={}", frame.index),
                ms_level: 1,
                rt_seconds: Some(frame.rt_in_seconds).filter(|rt| rt.is_finite()),
                precursor_mz: None,
                precursor_charge: None,
                mz,
                intensity,
                profile: false,
                parent: None,
                isolation: None,
                precursor_mobility: None,
                mobility,
            })
        }
    }

    impl Iterator for TdfScans {
        type Item = Result<PeakSpectrum, SpectraError>;

        /// The MS1 frame comes first at equal times: an MS/MS spectrum's
        /// precursor time is that of the MS1 frame it was picked from.
        fn next(&mut self) -> Option<Self::Item> {
            let ms1 = self.ms1_rt();
            let ms2 = match self.spectra.peek() {
                Some(Ok(spectrum)) => Some(spectrum.rt_seconds.unwrap_or(f64::NEG_INFINITY)),
                Some(Err(_)) => Some(f64::NEG_INFINITY),
                None => None,
            };
            match (ms1, ms2) {
                (Some(a), Some(b)) if b < a => self.spectra.next(),
                (Some(_), _) => Some(self.read_ms1()),
                (None, _) => self.spectra.next(),
            }
        }
    }

    /// One MS1 frame summed over its ion-mobility scans and centroided.
    pub(super) fn collapse(frame: &Frame, mz: &Tof2MzConverter) -> (Vec<f64>, Vec<f32>) {
        centroid_counts(
            &frame.tof_indices,
            &frame.intensities,
            frame.intensity_correction_factor,
            mz,
        )
    }

    /// One MS1 frame binned by ion mobility and centroided: its scans are
    /// summed `scans_per_bin` at a time, each bin centroided on its own, and
    /// each peak tagged with its bin's mobility (1/K0 at the bin's middle
    /// scan). Peaks are returned sorted by m/z.
    pub(super) fn bin(
        frame: &Frame,
        mz: &Tof2MzConverter,
        im: &Scan2ImConverter,
        scans_per_bin: usize,
    ) -> (Vec<f64>, Vec<f32>, Vec<f32>) {
        let scans = frame.scan_offsets.len().saturating_sub(1);
        let mut peaks: Vec<(f64, f32, f32)> = Vec::new();
        let mut first = 0;
        while first < scans {
            let last = (first + scans_per_bin).min(scans);
            let (start, end) = (frame.scan_offsets[first], frame.scan_offsets[last]);
            let (Some(tof), Some(counts)) = (
                frame.tof_indices.get(start..end),
                frame.intensities.get(start..end),
            ) else {
                break;
            };
            #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
            let mobility = im.convert((first + last - 1) as f64 / 2.0) as f32;
            let (bin_mz, bin_intensity) =
                centroid_counts(tof, counts, frame.intensity_correction_factor, mz);
            peaks.extend(
                bin_mz
                    .into_iter()
                    .zip(bin_intensity)
                    .map(|(m, i)| (m, i, mobility)),
            );
            first = last;
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        (
            peaks.iter().map(|p| p.0).collect(),
            peaks.iter().map(|p| p.1).collect(),
            peaks.iter().map(|p| p.2).collect(),
        )
    }

    /// Detector counts at the same TOF index add up; a gap of one or more
    /// empty indices separates peaks (zeros are put there so the centroider
    /// sees each run on its own), and each local maximum becomes a peak.
    fn centroid_counts(
        tof_indices: &[u32],
        intensities: &[u32],
        correction: f64,
        mz: &Tof2MzConverter,
    ) -> (Vec<f64>, Vec<f32>) {
        let mut counts: Vec<(u32, u64)> = tof_indices
            .iter()
            .zip(intensities)
            .map(|(&tof, &count)| (tof, u64::from(count)))
            .collect();
        counts.sort_unstable_by_key(|&(tof, _)| tof);
        let mut profile_mz = Vec::with_capacity(counts.len());
        let mut profile_intensity = Vec::with_capacity(counts.len());
        let mut last: Option<u32> = None;
        let factor = if correction > 0.0 { correction } else { 1.0 };
        let mut push = |tof: u32, value: f64| {
            profile_mz.push(mz.convert(tof));
            #[allow(clippy::cast_possible_truncation)]
            profile_intensity.push(value as f32);
        };
        let mut i = 0;
        while i < counts.len() {
            let tof = counts[i].0;
            let mut total = 0;
            while i < counts.len() && counts[i].0 == tof {
                total += counts[i].1;
                i += 1;
            }
            match last {
                Some(previous) if tof > previous + 1 => {
                    push(previous + 1, 0.0);
                    if tof > previous + 2 {
                        push(tof - 1, 0.0);
                    }
                }
                _ => {}
            }
            #[allow(clippy::cast_precision_loss)]
            push(tof, total as f64 * factor);
            last = Some(tof);
        }
        centroid(&profile_mz, &profile_intensity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peaks_inside_the_window_and_its_slack_are_possible() {
        let window = Some((300.0, 1600.0));
        // The slack is 5 on each side, inclusive: 295 and 1605 pass.
        let mz = [295.0, 300.0, 1600.0, 1605.0];
        assert_eq!(impossible_peaks(&mz, &[1.0; 4], window), None);
        // Intensity 0 and the 1e15 ceiling pass; so does a scan with no peaks.
        assert_eq!(
            impossible_peaks(&[400.0, 500.0], &[0.0, MAX_INTENSITY], window),
            None
        );
        assert_eq!(impossible_peaks(&[], &[], window), None);
        // With no window: from 1 to 1e6.
        assert_eq!(impossible_peaks(&[1.0, 1e6], &[1.0; 2], None), None);
    }

    #[test]
    fn peaks_no_instrument_could_measure_are_counted_and_the_first_described() {
        let window = Some((300.0, 1600.0));
        // Peak 1 is below the window less its slack, peak 3 isn't a number,
        // peak 4 is past the top; peak 2 has too intense a signal.
        let mz = [400.0, 294.999, 500.0, f64::NAN, 1605.001, 600.0];
        let intensity = [1.0, 1.0, 1.1e15, 1.0, 1.0, 1.0];

        let problem = impossible_peaks(&mz, &intensity, window).unwrap();

        assert!(
            problem.starts_with("4 of 6 peaks are impossible, the first (peak 1) at m/z 2.94999e2"),
            "{problem}"
        );
        // Without a window the proton's m/z and 1e6 bound it; 298.999 passes.
        assert_eq!(
            impossible_peaks(&[0.999_9, 298.999], &[1.0; 2], None)
                .unwrap()
                .split(',')
                .next(),
            Some("1 of 2 peaks are impossible")
        );
        // An intensity that isn't a number, or is negative.
        assert!(impossible_peaks(&[400.0], &[f32::NAN], window).is_some());
        assert!(impossible_peaks(&[400.0], &[f32::INFINITY], window).is_some());
        assert!(impossible_peaks(&[400.0], &[-0.5], window).is_some());
        assert!(impossible_peaks(&[400.0, 1_000_001.0], &[1.0; 2], None).is_some());
    }

    #[test]
    fn timstof_folders_are_spectra_when_they_hold_analysis_tdf() {
        let dir = tempfile::tempdir().unwrap();
        let run = dir.path().join("Sample35_2087.d");
        std::fs::create_dir(&run).unwrap();
        // A `.d` folder without analysis.tdf (e.g. Agilent's) isn't ours.
        assert!(!is_spectra_file(&run));
        std::fs::write(run.join("analysis.tdf"), b"").unwrap();
        assert!(is_spectra_file(&run));
        assert_eq!(spectra_run_name(&run), "Sample35_2087");
    }

    #[test]
    fn an_ms1_frame_is_summed_over_mobility_and_centroided() {
        // m/z = (10 + 0.01 * tof)^2: tof 499, 500, 501 -> 224.7001, 225,
        // 225.3001; tof 900 -> 361.
        let mz_converter =
            timsrust::converters::Tof2MzConverter::from_boundaries(100.0, 400.0, 1000);
        // Two mobility scans; tof 500 appears in both.
        let frame = timsrust::Frame {
            scan_offsets: vec![0, 3, 5],
            tof_indices: vec![499, 500, 900, 500, 501],
            intensities: vec![2, 10, 5, 6, 3],
            index: 7,
            rt_in_seconds: 12.5,
            ms_level: timsrust::MSLevel::MS1,
            intensity_correction_factor: 2.0,
            ..timsrust::Frame::default()
        };

        let (mz, intensity) = tdf::collapse(&frame, &mz_converter);

        // Summed counts 2, 16, 3 at tof 499..501 (x2 correction: 4, 32, 6),
        // then a gap, then 5 (x2: 10) at tof 900. The first peak's apex is
        // the parabola through ln 4, ln 32, ln 6: offset
        // 0.5 (ln 4 - ln 6) / (ln 4 - 2 ln 32 + ln 6) = 0.054013 steps of
        // (225.3001 - 224.7001) / 2 = 0.3, so 225.016204; its area is the
        // whole run, 42. The lone point at tof 900 stays where it is.
        assert_eq!(intensity, [42.0, 10.0]);
        assert!((mz[0] - 225.016_204).abs() < 1e-5, "{}", mz[0]);
        assert!((mz[1] - 361.0).abs() < 1e-9, "{}", mz[1]);

        // Binned two scans at a time it is the same frame, at the bin's
        // mobility: scan 0.5 of 1.3 - 0.5 * scan / 1000, 1.29975.
        let im = timsrust::converters::Scan2ImConverter::from_boundaries(0.8, 1.3, 1000);
        let (whole_mz, whole, mobility) = tdf::bin(&frame, &mz_converter, &im, 2);
        assert_eq!((whole_mz, whole), (mz, intensity));
        assert_eq!(mobility, [1.299_75_f32; 2]);

        // One scan per bin: scan 0 holds tof 499, 500 (2, 10; x2: 24) and
        // 900 (x2: 10), scan 1 tof 500, 501 (6, 3; x2: 18), each its own
        // peak, sorted by m/z, at 1.3 and 1.2995.
        let (by_scan_mz, by_scan, mobility) = tdf::bin(&frame, &mz_converter, &im, 1);
        assert_eq!(by_scan, [24.0, 18.0, 10.0]);
        assert!(by_scan_mz.is_sorted());
        assert_eq!(mobility, [1.3_f32, 1.2995, 1.3]);
    }

    #[test]
    fn a_pasef_spectrum_is_keyed_by_its_precursor() {
        let spectrum = timsrust::Spectrum {
            mz_values: vec![226.2136, 241.5291],
            intensities: vec![9.0, 31.0],
            precursor: Some(timsrust::Precursor {
                mz: 973.0248,
                rt: 1544.83,
                im: 1.18,
                charge: Some(2),
                intensity: Some(4627.0),
                index: 69656,
                frame_index: 14560,
            }),
            index: 69655,
            collision_energy: 42.0,
            isolation_mz: 973.0,
            isolation_width: 2.0,
        };
        let peaks = tdf::convert(&spectrum, 69655);
        assert_eq!(peaks.key, SpectrumKey::Scan(69656));
        assert_eq!(peaks.native_id, "precursor=69656");
        assert_eq!(peaks.ms_level, 2);
        assert_eq!(peaks.rt_seconds, Some(1544.83));
        assert_eq!(peaks.precursor_mz, Some(973.0248));
        assert_eq!(peaks.precursor_charge, Some(2));
        assert_eq!(peaks.intensity, [9.0, 31.0]);
        assert_eq!(peaks.precursor_mobility, Some(1.18));
        let orphan = timsrust::Spectrum {
            precursor: None,
            ..spectrum
        };
        assert_eq!(tdf::convert(&orphan, 7).key, SpectrumKey::Index(7));
    }

    #[test]
    fn run_names_drop_the_spectra_extension() {
        assert_eq!(spectra_run_name(Path::new("/d/HeLa_01.mzML")), "HeLa_01");
        assert_eq!(spectra_run_name(Path::new("x.MGF")), "x");
        assert!(is_spectra_file(Path::new("a.mzml")));
        assert!(!is_spectra_file(Path::new("a.mzid")));
    }
}
