//! Throughput of the MS2 model: `TUFFOLD_MODELS=… cargo run --release --example bench`.
use std::time::Instant;
use tuffold_core::PeptideForm;
use tuffold_predict::{Instrument, Models, encode};

fn main() {
    let dir = std::env::var_os("TUFFOLD_MODELS").expect("TUFFOLD_MODELS");
    let models = Models::load(std::path::Path::new(&dir))
        .unwrap()
        .expect("installed");
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    let forms: Vec<PeptideForm> = (0..40960)
        .map(|i| PeptideForm {
            sequence: (0..12)
                .map(|j| letters[(i * 7 + j * 13) % 20] as char)
                .collect(),
            modifications: vec![],
        })
        .collect();
    let encoded: Vec<_> = forms.iter().map(|f| encode(f).unwrap()).collect();
    let refs: Vec<_> = encoded.iter().collect();
    let start = Instant::now();
    for chunk in refs.chunks(256) {
        models
            .fragment_intensities(chunk, 2, 28.0, Instrument::QExactive)
            .unwrap();
    }
    let seconds = start.elapsed().as_secs_f64();
    println!("{:.0} peptides/s", refs.len() as f64 / seconds);
}
