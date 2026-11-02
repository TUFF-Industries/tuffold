//! Spectra reader tests.
//!
//! Expectations come from reading the fixture files directly: the mzML arrays
//! were decoded independently (base64 -> little-endian doubles with Python's
//! standard library), never by running this reader.
//! - tiny.pwiz.1.1.mzML: 4 spectra; only "scan=20" is MS2 (RT 5.9905 min,
//!   selected ion m/z 445.34, charge 2, m/z 0,2,...,18, intensity 20,18,...,2).
//!   "scan=19" is MS1 (RT 5.8905 min, m/z 0,1,...,14, intensity 15,14,...,1),
//!   "scan=21" is MS1 with no peaks, and the last (a WIFF-style id) is MS1.
//!   scan=20's precursor has spectrumRef="scan=19" (line 174).
//! - test.mgf: "Spectrum 1" (CHARGE=2+, PEPMASS=983.6, 6 peaks) and
//!   "Spectrum 2" (SCANS=3, PEPMASS=1084.9, RTINSECONDS=25, 6 peaks).

use std::path::{Path, PathBuf};

use tuffold_core::{PeakSpectrum, SpectrumKey};
use tuffold_io::{SpectraError, read_scans, read_spectra, spectra_run_name};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/third_party/pyteomics")
        .join(name)
}

fn read(name: &str) -> Vec<PeakSpectrum> {
    read_spectra(&fixture(name))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
}

#[test]
fn mzml_yields_only_msms_with_peaks_and_precursor() {
    let spectra = read("tiny.pwiz.1.1.mzML");

    assert_eq!(spectra.len(), 1);
    let scan20 = &spectra[0];
    assert_eq!(scan20.key, SpectrumKey::Scan(20));
    assert_eq!(scan20.native_id, "scan=20");
    assert_eq!(scan20.ms_level, 2);
    // 5.9905 minutes = 359.43 seconds.
    assert!((scan20.rt_seconds.unwrap() - 359.43).abs() < 1e-6);
    assert!((scan20.precursor_mz.unwrap() - 445.34).abs() < 1e-9);
    assert_eq!(scan20.precursor_charge, Some(2));
    let mz: Vec<f64> = (0..10).map(|i| f64::from(i) * 2.0).collect();
    let intensity: Vec<f32> = (0..10).map(|i| 20.0 - i as f32 * 2.0).collect();
    assert_eq!((&scan20.mz, &scan20.intensity), (&mz, &intensity));
    // scan=20 declares MS:1000128 "profile spectrum".
    assert!(scan20.profile);
    assert_eq!(scan20.parent.as_deref(), Some("scan=19"));
}

#[test]
fn scans_include_ms1_in_file_order() {
    let scans: Vec<_> = read_scans(&fixture("tiny.pwiz.1.1.mzML"))
        .unwrap()
        .collect();

    // The file is hand-made: its two MS1 scans with peaks (scan=19 and the
    // WIFF-style one) hold m/z 0,1,...,14 against declared windows of
    // 400-1800 and 100-1000, which no instrument produces, so they are
    // reported damaged. scan=20 (MS/MS) and scan=21 (MS1, no peaks) are read.
    let outcome: Vec<String> = scans
        .iter()
        .map(|scan| match scan {
            Ok(spectrum) => format!("ok {} level {}", spectrum.native_id, spectrum.ms_level),
            Err(SpectraError::DamagedMs1 { scan, .. }) => format!("damaged {scan}"),
            Err(other) => panic!("{other}"),
        })
        .collect();
    assert_eq!(
        outcome,
        [
            "damaged scan=19",
            "ok scan=20 level 2",
            "ok scan=21 level 1",
            "damaged sample=1 period=1 cycle=22 experiment=1",
        ]
    );
    assert!(scans[2].as_ref().unwrap().mz.is_empty());
    // MGF has no MS1: the same spectra as read_spectra, all MS/MS.
    let mgf: Vec<u8> = read_scans(&fixture("test.mgf"))
        .unwrap()
        .map(|s| s.unwrap().ms_level)
        .collect();
    assert_eq!(mgf, [2, 2]);
}

#[test]
fn mgf_keys_follow_scans_then_file_position() {
    let spectra = read("test.mgf");

    assert_eq!(spectra.len(), 2);
    let first = &spectra[0];
    assert_eq!(first.key, SpectrumKey::Index(0));
    assert_eq!(first.native_id, "Spectrum 1");
    assert_eq!(
        (first.precursor_mz, first.precursor_charge),
        (Some(983.6), Some(2))
    );
    assert_eq!(first.rt_seconds, None);
    assert_eq!(first.mz, [846.6, 846.8, 847.6, 1640.1, 1640.6, 1895.5]);
    assert_eq!(first.intensity, [73.0, 44.0, 67.0, 291.0, 54.0, 49.0]);

    let second = &spectra[1];
    assert_eq!(second.key, SpectrumKey::Scan(3));
    assert_eq!(second.precursor_mz, Some(1084.9));
    assert!((second.rt_seconds.unwrap() - 25.0).abs() < 1e-9);
    // The charge column (3, 2, 1, ...) does not change the m/z values.
    let expected = [345.1, 370.2, 460.2, 1673.3, 1674.0, 1675.3];
    assert_eq!(second.mz.len(), expected.len());
    for (actual, expected) in second.mz.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-9, "{actual} vs {expected}");
    }
    assert_eq!(second.intensity, [237.0, 128.0, 108.0, 1007.0, 974.0, 79.0]);
}

#[test]
fn run_names_match_identification_conventions() {
    assert_eq!(
        spectra_run_name(&fixture("tiny.pwiz.1.1.mzML")),
        "tiny.pwiz.1.1"
    );
}

#[test]
fn non_spectra_files_are_rejected() {
    let error = read_spectra(&fixture("test.mzid")).err().unwrap();
    assert!(error.to_string().contains("not a spectra file"), "{error}");
}

#[test]
fn truncated_mzml_is_an_error_not_a_short_success() {
    // Cut tiny.pwiz.1.1.mzML inside its third spectrum (scan=21): the file
    // still declares <spectrumList count="4">.
    let text = std::fs::read_to_string(fixture("tiny.pwiz.1.1.mzML")).unwrap();
    let cut = text.find("<spectrum index=\"2\"").unwrap() + 40;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("cut.mzML");
    std::fs::write(&path, &text[..cut]).unwrap();

    let results: Vec<_> = read_spectra(&path).unwrap().collect();

    // scan=20 (the only MS/MS before the cut) still arrives, then the error.
    assert!(matches!(&results[0], Ok(spectrum) if spectrum.key == SpectrumKey::Scan(20)));
    let error = results.last().unwrap().as_ref().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("of the 4 spectra the file declares"),
        "{error}"
    );
}

// --- damaged MS1 scans -----------------------------------------------------
//
// A small mzML written here, so every expectation follows from the values put
// in it. Arrays are uncompressed 64-bit floats; `base64` is checked against
// the encoding of 0, 1, 2 in tiny.pwiz.1.1.mzML (decoded independently).
// Scan windows are 300-1600 unless a scan says otherwise. Per scan:
//   1 MS1 400/500/600, intensities 10/20/30                      -> kept
//   2 MS1 with an m/z of 1.58e-19 (below any proton)             -> damaged
//   3 MS1 with an intensity of 7.4e34 (above 1e15)               -> damaged
//   4 MS1 with a NaN m/z                                         -> damaged
//   5 MS1 with an m/z of 2000, past 1600 + 5                     -> damaged
//   6 MS1 with m/z 299.998 and 1600.5, calibration, in the slack -> kept
//   7 MS1 with no window and an m/z of 1e7 (> 1e6)               -> damaged
//   8 MS1 with no window and m/z 400, 5000                       -> kept
//   9 MS2 whose precursor is scan 2, clean peaks                 -> kept
//  10 MS2 holding the garbage of scan 2                          -> kept
//  11 MS1 with a negative intensity                              -> damaged

fn base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut text = String::new();
    for chunk in bytes.chunks(3) {
        let word = chunk
            .iter()
            .enumerate()
            .fold(0u32, |word, (i, &b)| word | (u32::from(b) << (16 - 8 * i)));
        for i in 0..4 {
            if i <= chunk.len() {
                text.push(ALPHABET[((word >> (18 - 6 * i)) & 63) as usize] as char);
            } else {
                text.push('=');
            }
        }
    }
    text
}

fn encode(values: &[f64]) -> String {
    base64(
        &values
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<_>>(),
    )
}

fn spectrum(
    scan: u32,
    level: u8,
    window: Option<(f64, f64)>,
    mz: &[f64],
    intensity: &[f64],
) -> String {
    let window = window.map_or(String::new(), |(low, high)| {
        format!(
            r#"<scanWindowList count="1"><scanWindow>
<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="{low}"/>
<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="{high}"/>
</scanWindow></scanWindowList>"#
        )
    });
    let precursor = if level == 2 {
        r#"<precursorList count="1"><precursor spectrumRef="scan=2"><selectedIonList count="1">
<selectedIon><cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="500"/>
<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/></selectedIon>
</selectedIonList></precursor></precursorList>"#
    } else {
        ""
    };
    let array = |accession: &str, name: &str, values: &[f64]| {
        format!(
            r#"<binaryDataArray encodedLength="0">
<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
<cvParam cvRef="MS" accession="{accession}" name="{name}" value=""/>
<binary>{}</binary></binaryDataArray>"#,
            encode(values)
        )
    };
    format!(
        r#"<spectrum index="{index}" id="scan={scan}" defaultArrayLength="{n}">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
<scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{scan}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
{window}</scan></scanList>
{precursor}
<binaryDataArrayList count="2">{}{}</binaryDataArrayList></spectrum>"#,
        array("MS:1000514", "m/z array", mz),
        array("MS:1000515", "intensity array", intensity),
        index = scan - 1,
        n = mz.len(),
    )
}

/// Garbage seen in ThermoRawFileParser 2.0.0-dev's output (mixtureArep2 scan
/// 1914 holds both).
const TINY_MZ: f64 = 1.578_803_058_287_252_6e-19;
const HUGE_INTENSITY: f64 = 7.413_142_472_324_203e34;

fn damaged_scans_mzml() -> String {
    let w = Some((300.0, 1600.0));
    let spectra = [
        spectrum(1, 1, w, &[400.0, 500.0, 600.0], &[10.0, 20.0, 30.0]),
        spectrum(2, 1, w, &[TINY_MZ, 500.0, 600.0], &[10.0, 20.0, 30.0]),
        spectrum(
            3,
            1,
            w,
            &[400.0, 500.0, 600.0],
            &[10.0, HUGE_INTENSITY, 30.0],
        ),
        spectrum(4, 1, w, &[400.0, f64::NAN, 600.0], &[10.0, 20.0, 30.0]),
        spectrum(5, 1, w, &[400.0, 2000.0, 600.0], &[10.0, 20.0, 30.0]),
        spectrum(6, 1, w, &[299.998, 1600.5], &[10.0, 20.0]),
        spectrum(7, 1, None, &[400.0, 1e7], &[10.0, 20.0]),
        spectrum(8, 1, None, &[400.0, 5000.0], &[10.0, 20.0]),
        spectrum(9, 2, w, &[200.0, 300.0], &[5.0, 6.0]),
        spectrum(10, 2, w, &[TINY_MZ, 300.0], &[HUGE_INTENSITY, 6.0]),
        spectrum(11, 1, w, &[400.0, 500.0, 600.0], &[10.0, -5.0, 30.0]),
    ];
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0" id="damaged">
<cvList count="2"><cv id="MS" fullName="MS" version="4.1" URI="x"/><cv id="UO" fullName="UO" version="1" URI="x"/></cvList>
<run id="damaged"><spectrumList count="{}">{}</spectrumList></run></mzML>"#,
        spectra.len(),
        spectra.concat()
    )
}

fn write_damaged() -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("damaged.mzML");
    std::fs::write(&path, damaged_scans_mzml()).unwrap();
    (dir, path)
}

#[test]
fn the_test_base64_matches_the_fixtures_encoding() {
    assert_eq!(encode(&[0.0, 1.0, 2.0]), "AAAAAAAAAAAAAAAAAADwPwAAAAAAAABA");
    // RFC 4648 vectors, with padding.
    assert_eq!(base64(b"Ma"), "TWE=");
    assert_eq!(base64(b"M"), "TQ==");
}

#[test]
fn mzml_ms1_scans_with_impossible_peaks_are_damaged_and_msms_is_not_judged() {
    let (_dir, path) = write_damaged();

    let scans: Vec<_> = read_scans(&path).unwrap().collect();

    // Scans arrive in file order, all 11 of them (the damaged ones as
    // errors), and no truncation error follows.
    assert_eq!(scans.len(), 11);
    let outcome: Vec<String> = scans
        .iter()
        .map(|scan| match scan {
            Ok(spectrum) => format!("ok {}", spectrum.native_id),
            Err(SpectraError::DamagedMs1 { scan, .. }) => format!("damaged {scan}"),
            Err(other) => panic!("{other}"),
        })
        .collect();
    assert_eq!(
        outcome,
        [
            "ok scan=1",
            "damaged scan=2",
            "damaged scan=3",
            "damaged scan=4",
            "damaged scan=5",
            "ok scan=6",
            "damaged scan=7",
            "ok scan=8",
            "ok scan=9",
            "ok scan=10",
            "damaged scan=11",
        ]
    );
    // The message names the file, the scan, how many peaks, and the first.
    let message = scans[1].as_ref().unwrap_err().to_string();
    assert!(
        message.starts_with("damaged.mzML: MS1 scan=2 is damaged (1 of 3 peaks are impossible"),
        "{message}"
    );
    assert!(
        message.contains("the first (peak 0) at m/z 1.5788030582872526e-19"),
        "{message}"
    );
    // A scan inside the calibration slack arrives whole, and a kept MS1 scan
    // is converted as before: level, key, no parent, minutes as seconds,
    // intensities as f32.
    assert_eq!(scans[5].as_ref().unwrap().mz, [299.998, 1600.5]);
    let first = scans[0].as_ref().unwrap();
    assert_eq!(
        (first.ms_level, &first.key, &first.parent, first.rt_seconds),
        (1, &SpectrumKey::Scan(1), &None, Some(60.0))
    );
    assert_eq!(first.mz, [400.0, 500.0, 600.0]);
    assert_eq!(first.intensity, [10.0, 20.0, 30.0]);
}

#[test]
fn msms_whose_precursor_scan_is_damaged_is_still_read_as_it_is() {
    let (_dir, path) = write_damaged();

    let spectra: Vec<PeakSpectrum> = read_spectra(&path)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();

    // Scans 9 and 10 are MS/MS of the damaged scan 2: both arrive, and scan
    // 10 keeps its own garbage peaks (MS/MS is not judged).
    let ids: Vec<&str> = spectra.iter().map(|s| s.native_id.as_str()).collect();
    assert_eq!(ids, ["scan=9", "scan=10"]);
    assert_eq!(spectra[0].parent.as_deref(), Some("scan=2"));
    assert_eq!(spectra[1].mz, [TINY_MZ, 300.0]);
}
