//! Peak storage for identified spectra.
//!
//! Blob layout (format 1), zlib-compressed as a whole:
//! `[format: u8 = 1][count: u32 LE][m/z: count × f64 LE][intensity: count × f32 LE]`.
//! m/z keeps full double precision (fragment annotation needs it); intensity
//! is single precision, as every common spectra format stores it. Peaks are
//! stored sorted by m/z (annotation searches them by bisection), whatever
//! order the spectra file recorded.

use std::io::{Read, Write};

use flate2::Compression;
use flate2::read::ZlibDecoder;
use flate2::write::ZlibEncoder;

const FORMAT: u8 = 1;

/// A stored spectrum's peaks, sorted by m/z.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Peaks {
    pub mz: Vec<f64>,
    pub intensity: Vec<f32>,
}

pub(crate) fn encode(mz: &[f64], intensity: &[f32]) -> Vec<u8> {
    let count = mz.len().min(intensity.len());
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| mz[a].total_cmp(&mz[b]));
    let mut raw = Vec::with_capacity(5 + count * 12);
    raw.push(FORMAT);
    raw.extend_from_slice(&u32::try_from(count).unwrap_or(u32::MAX).to_le_bytes());
    for &index in &order {
        raw.extend_from_slice(&mz[index].to_le_bytes());
    }
    for &index in &order {
        raw.extend_from_slice(&intensity[index].to_le_bytes());
    }
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    // Writing into a Vec cannot fail.
    encoder.write_all(&raw).expect("in-memory write");
    encoder.finish().expect("in-memory write")
}

/// Decode a stored blob; `None` when it is corrupt or of an unknown format.
pub(crate) fn decode(blob: &[u8]) -> Option<Peaks> {
    let mut raw = Vec::new();
    ZlibDecoder::new(blob).read_to_end(&mut raw).ok()?;
    let (&format, rest) = raw.split_first()?;
    if format != FORMAT || rest.len() < 4 {
        return None;
    }
    let count = usize::try_from(u32::from_le_bytes(rest[..4].try_into().ok()?)).ok()?;
    let body = &rest[4..];
    if body.len() != count * 12 {
        return None;
    }
    let (mz_bytes, intensity_bytes) = body.split_at(count * 8);
    let mz = mz_bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|chunk| f64::from_le_bytes(*chunk))
        .collect();
    let intensity = intensity_bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|chunk| f32::from_le_bytes(*chunk))
        .collect();
    Some(Peaks { mz, intensity })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_exactly() {
        let mz = vec![100.0, 200.123_456_789_012, 1_999.999_999];
        let intensity = vec![1.5, 0.0, 3.25e6];
        assert_eq!(
            decode(&encode(&mz, &intensity)),
            Some(Peaks { mz, intensity })
        );
    }

    #[test]
    fn peaks_are_stored_sorted_by_mz() {
        let decoded = decode(&encode(&[300.0, 100.0, 200.0], &[3.0, 1.0, 2.0])).unwrap();
        assert_eq!(decoded.mz, [100.0, 200.0, 300.0]);
        assert_eq!(decoded.intensity, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn empty_spectra_round_trip() {
        assert_eq!(decode(&encode(&[], &[])), Some(Peaks::default()));
    }

    #[test]
    fn corrupt_blobs_decode_to_none() {
        assert_eq!(decode(b"not zlib"), None);
        let mut blob = encode(&[1.0], &[2.0]);
        blob.truncate(blob.len() - 3);
        assert_eq!(decode(&blob), None);
    }
}
