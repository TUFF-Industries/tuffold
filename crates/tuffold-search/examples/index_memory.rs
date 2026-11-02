//! Build only the index for a FASTA and print what was estimated, to check
//! the estimate against the real peak (`/usr/bin/time -l` on macOS,
//! `/usr/bin/time -v` on Linux):
//!
//! `cargo run --release -p tuffold-search --example index_memory -- FASTA [--semi] [--budget-gb N]`

use std::path::Path;

use tuffold_search::{SearchOptions, SearchSettings, Searcher};

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let fasta = args
        .first()
        .ok_or("usage: index_memory FASTA [--semi] [--budget-gb N]")?;
    let budget = args
        .iter()
        .position(|arg| arg == "--budget-gb")
        .and_then(|at| args.get(at + 1))
        .map(|gb| gb.parse::<f64>().map_err(|e| e.to_string()))
        .transpose()?;
    let settings = SearchSettings::from_options(&SearchOptions {
        semi_specific: args.iter().any(|arg| arg == "--semi"),
        max_memory_gb: budget,
        ..SearchOptions::default()
    })?;
    let records = tuffold_io::fasta::parse_fasta(Path::new(fasta)).map_err(|e| e.to_string())?;
    let started = std::time::Instant::now();
    let searcher = Searcher::new(&records, settings).map_err(|e| e.to_string())?;
    let (estimate, budget) = searcher.index_estimate();
    println!(
        "built in {:.1} s: {} peptides indexed; estimated {} digests, {} forms, {} fragments, \
         peak {:.2} GB (budget {:.2} GB)",
        started.elapsed().as_secs_f64(),
        searcher.peptides(),
        estimate.digests,
        estimate.forms,
        estimate.fragments,
        estimate.peak_bytes as f64 / 1e9,
        budget as f64 / 1e9
    );
    Ok(())
}
