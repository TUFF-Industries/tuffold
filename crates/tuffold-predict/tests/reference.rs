//! The models against a PyTorch implementation of AlphaPeptDeep's layers,
//! on its released `generic` weights. Runs when `TUFFOLD_MODELS` names a
//! folder holding them (`tuffold models --install`); weights never enter the
//! repository.

use tuffold_core::{ModSite, Modification, PeptideForm};
use tuffold_predict::{Instrument, Models, encode};

fn models() -> Option<Models> {
    let dir = std::env::var_os("TUFFOLD_MODELS")?;
    Models::load(std::path::Path::new(&dir)).expect("models load")
}

fn form(sequence: &str, modifications: &[(ModSite, f64)]) -> PeptideForm {
    PeptideForm {
        sequence: sequence.into(),
        modifications: modifications
            .iter()
            .map(|&(site, mass_delta)| Modification {
                site,
                mass_delta,
                name: None,
                accession: None,
            })
            .collect(),
    }
}

/// (peptide, charge, pretrained RT, first and last rows of b⁺ b²⁺ y⁺ y²⁺ at
/// NCE 28 on a Q Exactive), from PyTorch.
#[allow(clippy::type_complexity)]
fn cases() -> Vec<(PeptideForm, u8, f32, [f32; 4], [f32; 4])> {
    vec![
        (
            form("LDHESFGGPFK", &[]),
            2,
            0.381_138,
            [0.0, 0.0, 0.049_454, 0.300_859],
            [0.016_661, 0.0, 0.496_695, 0.0],
        ),
        (
            form("AAAAATTMK", &[(ModSite::Residue(8), 15.994_915)]),
            2,
            -0.000_524,
            [0.0, 0.0, 0.028_786, 0.0],
            [0.0, 0.0, 0.128_118, 0.0],
        ),
        (
            form(
                "ACDEFGHIKR",
                &[
                    (ModSite::NTerm, 42.010_565),
                    (ModSite::Residue(2), 57.021_464),
                ],
            ),
            3,
            0.282_450,
            [0.402_452, 0.0, 0.0, 0.235_863],
            [0.0, 0.0, 0.808_72, 0.0],
        ),
    ]
}

#[test]
fn predictions_match_pytorch() {
    let Some(models) = models() else {
        eprintln!("TUFFOLD_MODELS not set; skipping");
        return;
    };
    for (form, charge, rt, first, last) in cases() {
        let encoded = encode(&form).expect("encodes");
        let rows = models
            .fragment_intensities(&[&encoded], charge, 28.0, Instrument::QExactive)
            .unwrap()
            .remove(0);
        assert_eq!(rows.len(), form.sequence.len() - 1);
        for (got, want) in [(rows[0], first), (rows[rows.len() - 1], last)] {
            for (g, w) in got.iter().zip(want) {
                // The reference rounded to 6 decimals and printed tiny
                // negatives as -0.0.
                assert!(
                    (g.max(0.0) - w).abs() < 2e-5,
                    "{}: {got:?} vs {want:?}",
                    form.sequence
                );
            }
        }
        let predicted = models.pretrained_retention(&[&encoded]).unwrap()[0];
        assert!(
            (predicted - rt).abs() < 2e-5,
            "{}: {predicted} vs {rt}",
            form.sequence
        );
    }
}
