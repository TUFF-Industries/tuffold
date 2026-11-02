//! The MS/MS peaks a Scaffold `.sfdb` stores (`docs/sfdb-format.md`).
//!
//! Two layouts. Scaffold 5 writes them as text on the `Spectrum` row
//! (`MASS_ARRAY`, `INTENSITY_ARRAY`: comma-separated, at most 300). Scaffold
//! DDA 6.5 leaves those empty and writes one `SpectrumData` row per spectrum:
//! a protobuf message `{1: format name, 2: payload}` whose payload, gzip'd
//! (`FORMAT` = `gzip`), is `{1: {2: {1: scan, 2: title}}}` with the scan
//! `{1: precursor m/z (f64), 2: time (f32), 3: precursor intensity (f32),
//! 4: {1: peak…}, 5: …}` and each peak `{1: m/z (f64), 2: time (f32),
//! 3: intensity (f32)}`. Worked out from one file (PXD071927); the decoder
//! takes only that path and refuses anything else.

use std::collections::HashSet;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use rusqlite::{Connection, OpenFlags, OptionalExtension};
use tuffold_core::{PeakSpectrum, SpectrumKey};

use crate::sf3::UnreadablePeaks;

/// A decompressed spectrum larger than this is refused (a peak is 21 bytes
/// of protobuf; 16 MB is ~800,000 peaks, the real ones hold hundreds).
const MAX_INFLATED: u64 = 16 << 20;
/// A message with more fields than this is refused: fields are held as
/// they are parsed, so this bounds the memory a hostile blob can take.
const MAX_FIELDS: usize = 1 << 20;

/// Where one identified spectrum's peaks are, and what it is.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PeakRef {
    pub run: String,
    pub key: SpectrumKey,
    pub sample: i64,
    pub spectrum_id: String,
    pub native_id: String,
    pub rt_seconds: Option<f64>,
    pub precursor_mz: Option<f64>,
    pub precursor_charge: Option<i16>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Layout {
    /// `Spectrum.MASS_ARRAY` / `INTENSITY_ARRAY` text (Scaffold 5).
    Text,
    /// `SpectrumData` blobs (Scaffold DDA).
    Blob,
}

/// The peaks of a file's identified spectra, read run by run from the file.
#[derive(Clone)]
pub struct SfdbPeaks {
    inner: Arc<Stored>,
}

struct Stored {
    path: PathBuf,
    layout: Layout,
    spectra: Vec<PeakRef>,
}

impl std::fmt::Debug for SfdbPeaks {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SfdbPeaks({} spectra)", self.inner.spectra.len())
    }
}

/// Two handles are equal when they hold the same stored peaks.
impl PartialEq for SfdbPeaks {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
    }
}

impl SfdbPeaks {
    pub(crate) fn new(path: &Path, layout: Layout, spectra: Vec<PeakRef>) -> Self {
        Self {
            inner: Arc::new(Stored {
                path: path.to_path_buf(),
                layout,
                spectra,
            }),
        }
    }

    /// The runs holding spectra with peaks, in first-seen order.
    pub fn runs(&self) -> Vec<&str> {
        let mut seen = HashSet::new();
        self.inner
            .spectra
            .iter()
            .map(|spectrum| spectrum.run.as_str())
            .filter(|run| seen.insert(*run))
            .collect()
    }

    /// Identified spectra (any run).
    pub fn len(&self) -> usize {
        self.inner.spectra.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.spectra.is_empty()
    }

    /// The run's spectra with their peaks, sorted by m/z. A spectrum without
    /// stored peaks is left out: it has none to show.
    pub fn spectra(&self, run: &str) -> Vec<Result<PeakSpectrum, UnreadablePeaks>> {
        let refs = self
            .inner
            .spectra
            .iter()
            .filter(|spectrum| spectrum.run == run);
        let conn = match Connection::open_with_flags(
            &self.inner.path,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        ) {
            Ok(conn) => conn,
            Err(_) => {
                return refs
                    .map(|spectrum| Err(unreadable(spectrum, "the file can't be opened again")))
                    .collect();
            }
        };
        refs.filter_map(|spectrum| match read(&conn, self.inner.layout, spectrum) {
            Ok(Some(peaks)) if peaks.mz.is_empty() => None,
            Ok(Some(peaks)) => Some(Ok(peaks)),
            Ok(None) => None,
            Err(problem) => Some(Err(unreadable(spectrum, problem))),
        })
        .collect()
    }
}

fn unreadable(spectrum: &PeakRef, problem: &'static str) -> UnreadablePeaks {
    UnreadablePeaks {
        native_id: spectrum.native_id.clone(),
        problem,
    }
}

fn read(
    conn: &Connection,
    layout: Layout,
    spectrum: &PeakRef,
) -> Result<Option<PeakSpectrum>, &'static str> {
    let pairs = match layout {
        Layout::Text => {
            let arrays: Option<(Option<String>, Option<String>)> = conn
                .prepare_cached(
                    "SELECT MASS_ARRAY, INTENSITY_ARRAY FROM Spectrum \
                     WHERE MSSAMPLE_ID = ?1 AND SPECTRUM_ID = ?2",
                )
                .and_then(|mut statement| {
                    statement.query_row((spectrum.sample, &spectrum.spectrum_id), |row| {
                        Ok((row.get(0)?, row.get(1)?))
                    })
                })
                .optional()
                .map_err(|_| "its Spectrum row can't be read")?;
            let Some((Some(mz), Some(intensity))) = arrays else {
                return Ok(None);
            };
            text_pairs(&mz, &intensity)?
        }
        Layout::Blob => {
            let blob: Option<(Vec<u8>, Option<String>)> = conn
                .prepare_cached(
                    "SELECT DATA, FORMAT FROM SpectrumData \
                     WHERE MSSAMPLE_ID = ?1 AND SPECTRUM_ID = ?2",
                )
                .and_then(|mut statement| {
                    statement.query_row((spectrum.sample, &spectrum.spectrum_id), |row| {
                        Ok((row.get(0)?, row.get(1)?))
                    })
                })
                .optional()
                .map_err(|_| "its SpectrumData row can't be read")?;
            let Some((data, format)) = blob else {
                return Ok(None);
            };
            if format.as_deref() != Some("gzip") {
                return Err("its SpectrumData isn't in the gzip format Tuffold knows");
            }
            blob_pairs(&data)?
        }
    };
    Ok(Some(peak_spectrum(spectrum, pairs)))
}

fn peak_spectrum(spectrum: &PeakRef, mut pairs: Vec<(f64, f32)>) -> PeakSpectrum {
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (mz, intensity) = pairs.into_iter().unzip();
    PeakSpectrum {
        key: spectrum.key.clone(),
        native_id: spectrum.native_id.clone(),
        ms_level: 2,
        rt_seconds: spectrum.rt_seconds,
        precursor_mz: spectrum.precursor_mz,
        precursor_charge: spectrum.precursor_charge,
        mz,
        intensity,
        profile: false,
        parent: None,
        isolation: None,
        precursor_mobility: None,
        mobility: Vec::new(),
    }
}

/// Peaks from Scaffold 5's two comma-separated arrays, which must pair up.
pub(crate) fn text_pairs(mz: &str, intensity: &str) -> Result<Vec<(f64, f32)>, &'static str> {
    let values = |text: &str| -> Result<Vec<f64>, &'static str> {
        text.split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| {
                value
                    .parse::<f64>()
                    .map_err(|_| "its peak arrays hold a non-number")
            })
            .collect()
    };
    let (mz, intensity) = (values(mz)?, values(intensity)?);
    if mz.len() != intensity.len() {
        return Err("its m/z and intensity arrays differ in length");
    }
    checked(mz.into_iter().zip(intensity).map(|(mz, i)| (mz, i as f32)))
}

/// Peaks from a Scaffold DDA `SpectrumData` blob.
pub(crate) fn blob_pairs(data: &[u8]) -> Result<Vec<(f64, f32)>, &'static str> {
    const BAD: &str = "its SpectrumData doesn't have the layout Tuffold knows";
    let outer = fields(data).ok_or(BAD)?;
    let payload = only(&outer, 2).and_then(Field::bytes).ok_or(BAD)?;
    let mut inflated = Vec::new();
    flate2::read::GzDecoder::new(payload)
        .take(MAX_INFLATED + 1)
        .read_to_end(&mut inflated)
        .map_err(|_| "its SpectrumData doesn't decompress")?;
    if inflated.len() as u64 > MAX_INFLATED {
        return Err("its SpectrumData decompresses past Tuffold's limit");
    }
    let top = fields(&inflated).ok_or(BAD)?;
    let container = fields(only(&top, 1).and_then(Field::bytes).ok_or(BAD)?).ok_or(BAD)?;
    let group = fields(only(&container, 2).and_then(Field::bytes).ok_or(BAD)?).ok_or(BAD)?;
    let scan_bytes = group
        .iter()
        .find(|field| field.number == 1)
        .and_then(Field::bytes)
        .ok_or(BAD)?;
    let scan = fields(scan_bytes).ok_or(BAD)?;
    let mut pairs = Vec::new();
    for list in scan.iter().filter(|field| field.number == 4) {
        let list = fields(list.bytes().ok_or(BAD)?).ok_or(BAD)?;
        for peak in list.iter().filter(|field| field.number == 1) {
            let peak = fields(peak.bytes().ok_or(BAD)?).ok_or(BAD)?;
            let mz = only(&peak, 1).and_then(Field::double).ok_or(BAD)?;
            let intensity = only(&peak, 3).and_then(Field::float).ok_or(BAD)?;
            pairs.push((mz, intensity));
        }
    }
    checked(pairs.into_iter())
}

/// Every peak must be one an instrument could measure.
fn checked(pairs: impl Iterator<Item = (f64, f32)>) -> Result<Vec<(f64, f32)>, &'static str> {
    pairs
        .map(|(mz, intensity)| {
            if mz.is_finite() && mz > 0.0 && intensity.is_finite() && intensity >= 0.0 {
                Ok((mz, intensity))
            } else {
                Err("it holds a peak with an impossible m/z or intensity")
            }
        })
        .collect()
}

/// One protobuf field: its number and wire value.
struct Field<'a> {
    number: u64,
    value: Wire<'a>,
}

enum Wire<'a> {
    Varint,
    Fixed64([u8; 8]),
    Bytes(&'a [u8]),
    Fixed32([u8; 4]),
}

impl<'a> Field<'a> {
    fn bytes(&self) -> Option<&'a [u8]> {
        match self.value {
            Wire::Bytes(bytes) => Some(bytes),
            _ => None,
        }
    }

    fn double(&self) -> Option<f64> {
        match self.value {
            Wire::Fixed64(bytes) => Some(f64::from_le_bytes(bytes)),
            _ => None,
        }
    }

    fn float(&self) -> Option<f32> {
        match self.value {
            Wire::Fixed32(bytes) => Some(f32::from_le_bytes(bytes)),
            _ => None,
        }
    }
}

/// The field with this number, when there is exactly one.
fn only<'a, 'b>(fields: &'b [Field<'a>], number: u64) -> Option<&'b Field<'a>> {
    let mut found = fields.iter().filter(|field| field.number == number);
    let first = found.next()?;
    found.next().is_none().then_some(first)
}

/// A message's fields, or `None` when the bytes aren't one (bounds checked;
/// groups and unknown wire types are refused).
fn fields(mut bytes: &[u8]) -> Option<Vec<Field<'_>>> {
    let mut out = Vec::new();
    while !bytes.is_empty() {
        let (key, rest) = varint(bytes)?;
        bytes = rest;
        let number = key >> 3;
        if number == 0 {
            return None;
        }
        let value = match key & 7 {
            0 => {
                let (_, rest) = varint(bytes)?;
                bytes = rest;
                Wire::Varint
            }
            1 => {
                let (value, rest) = bytes.split_first_chunk::<8>()?;
                bytes = rest;
                Wire::Fixed64(*value)
            }
            2 => {
                let (length, rest) = varint(bytes)?;
                let length = usize::try_from(length).ok()?;
                if length > rest.len() {
                    return None;
                }
                let (value, rest) = rest.split_at(length);
                bytes = rest;
                Wire::Bytes(value)
            }
            5 => {
                let (value, rest) = bytes.split_first_chunk::<4>()?;
                bytes = rest;
                Wire::Fixed32(*value)
            }
            _ => return None,
        };
        if out.len() == MAX_FIELDS {
            return None;
        }
        out.push(Field { number, value });
    }
    Some(out)
}

fn varint(bytes: &[u8]) -> Option<(u64, &[u8])> {
    let mut value = 0u64;
    for (i, &byte) in bytes.iter().enumerate().take(10) {
        value |= u64::from(byte & 0x7f) << (7 * i);
        if byte < 0x80 {
            return Some((value, &bytes[i + 1..]));
        }
    }
    None
}

#[cfg(test)]
pub(crate) mod tests {
    use std::io::Write;

    use super::*;

    fn key(number: u64, wire: u64) -> Vec<u8> {
        encode_varint((number << 3) | wire)
    }

    fn encode_varint(mut value: u64) -> Vec<u8> {
        let mut out = Vec::new();
        loop {
            let byte = (value & 0x7f) as u8;
            value >>= 7;
            if value == 0 {
                out.push(byte);
                return out;
            }
            out.push(byte | 0x80);
        }
    }

    fn message(number: u64, body: &[u8]) -> Vec<u8> {
        let mut out = key(number, 2);
        out.extend(encode_varint(body.len() as u64));
        out.extend_from_slice(body);
        out
    }

    /// A `SpectrumData` blob as Scaffold DDA 6.5 writes one, for tests.
    pub(crate) fn blob(precursor: f64, peaks: &[(f64, f32)]) -> Vec<u8> {
        let mut list = Vec::new();
        for &(mz, intensity) in peaks {
            let mut peak = key(1, 1);
            peak.extend(mz.to_le_bytes());
            peak.extend(key(2, 5));
            peak.extend(1.5f32.to_le_bytes());
            peak.extend(key(3, 5));
            peak.extend(intensity.to_le_bytes());
            list.extend(message(1, &peak));
        }
        let mut scan = key(1, 1);
        scan.extend(precursor.to_le_bytes());
        scan.extend(key(3, 5));
        scan.extend(f32::NAN.to_le_bytes());
        scan.extend(message(4, &list));
        scan.extend(key(5, 0));
        scan.push(1);
        let mut group = message(1, &scan);
        group.extend(message(2, b"-7_7"));
        let top = message(1, &message(2, &group));
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&top).unwrap();
        let mut outer = message(1, b"gzip");
        outer.extend(message(2, &gz.finish().unwrap()));
        outer
    }

    #[test]
    fn blobs_decode_to_their_peaks() {
        let data = blob(397.21347, &[(126.127, 15313.37), (113.109, 1375.56)]);
        assert_eq!(
            blob_pairs(&data).unwrap(),
            vec![(126.127, 15313.37), (113.109, 1375.56)]
        );
    }

    #[test]
    fn damaged_blobs_are_refused_not_misread() {
        let data = blob(397.21347, &[(126.127, 15313.37)]);
        assert!(blob_pairs(&data[..data.len() - 3]).is_err());
        assert!(blob_pairs(b"not a blob").is_err());
        assert!(blob_pairs(&blob(400.0, &[(f64::NAN, 1.0)])).is_err());
        assert!(blob_pairs(&blob(400.0, &[(-5.0, 1.0)])).is_err());
    }

    #[test]
    fn text_arrays_pair_up() {
        assert_eq!(
            text_pairs("106.30591,110.07138,", "3984.0,23680.0,").unwrap(),
            vec![(106.30591, 3984.0), (110.07138, 23680.0)]
        );
        assert!(text_pairs("1.0,2.0", "3.0").is_err());
        assert!(text_pairs("1.0,x", "3.0,4.0").is_err());
    }

    #[test]
    fn varints_are_bounded() {
        assert_eq!(varint(&[0x96, 0x01]), Some((150, &[][..])));
        assert_eq!(varint(&[0xff; 11]), None);
        assert_eq!(varint(&[0x80]), None);
    }
}
