//! Scaffold `.sf3` import against real projects. Scaffold's files can't be
//! redistributed, so these run only when `TUFFOLD_SF3_SAMPLES` names a
//! folder of `.sf3` files (`cargo test -p tuffold-io --test sf3 --
//! --ignored --nocapture`). The checks hold for any valid project.

use std::collections::HashSet;
use std::path::PathBuf;

use tuffold_core::chemistry::peptide_mass;
use tuffold_io::import_sf3;

fn samples() -> Vec<PathBuf> {
    let Ok(folder) = std::env::var("TUFFOLD_SF3_SAMPLES") else {
        return Vec::new();
    };
    let mut files: Vec<PathBuf> = std::fs::read_dir(folder)
        .expect("readable sample folder")
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "sf3"))
        .collect();
    files.sort();
    files
}

#[test]
#[ignore = "needs Scaffold sample files (TUFFOLD_SF3_SAMPLES)"]
fn sample_projects_import_consistently() {
    let files = samples();
    assert!(!files.is_empty(), "no .sf3 files in TUFFOLD_SF3_SAMPLES");
    for file in files {
        let started = std::time::Instant::now();
        let imported =
            import_sf3(&file).unwrap_or_else(|error| panic!("{}: {error}", file.display()));
        println!(
            "== {} ({:.1} s)",
            file.display(),
            started.elapsed().as_secs_f64()
        );
        for note in &imported.notes {
            println!("   note: {note}");
        }
        assert!(
            !imported.batches.is_empty(),
            "{}: no engines",
            file.display()
        );
        for batch in &imported.batches {
            let runs: HashSet<&str> = batch.ms_runs.iter().map(|run| run.name.as_str()).collect();
            let spectra: HashSet<(&str, String)> = batch
                .spectra
                .iter()
                .map(|spectrum| (spectrum.ms_run.as_str(), spectrum.key.canonical()))
                .collect();
            assert_eq!(spectra.len(), batch.spectra.len(), "spectrum keys repeat");
            let rank_one = batch.psms.iter().filter(|psm| psm.rank == 1).count();
            let decoys = batch.decoy_index();
            let decoy_rank_one = batch
                .psms
                .iter()
                .filter(|psm| psm.rank == 1 && decoys.is_decoy(psm))
                .count();
            let untitled = batch
                .spectra
                .iter()
                .filter(|s| s.source_id.is_none())
                .count();
            println!("   {untitled} spectra without a title");
            let with_sequence = batch
                .proteins
                .iter()
                .filter(|p| p.sequence.is_some())
                .count();
            println!(
                "   {} {:?}: {} runs, {} spectra, {} PSMs ({rank_one} rank 1, {decoy_rank_one} decoy), \
                 {} proteins ({with_sequence} with sequence)",
                batch.software.name,
                batch.software.version,
                batch.ms_runs.len(),
                batch.spectra.len(),
                batch.psms.len(),
                batch.proteins.len(),
            );
            for psm in &batch.psms {
                assert!(runs.contains(psm.ms_run.as_str()));
                assert!(spectra.contains(&(psm.ms_run.as_str(), psm.spectrum.canonical())));
                assert!(
                    !psm.evidence.is_empty(),
                    "{} has no protein",
                    psm.peptide.proforma()
                );
                assert!(
                    !psm.scores.is_empty(),
                    "{} has no score",
                    psm.peptide.proforma()
                );
                if let (Some(calculated), Some(experimental)) =
                    (psm.calculated_mz, psm.experimental_mz)
                {
                    let error = (calculated - experimental).abs() * f64::from(psm.charge);
                    assert!(error < 4.1, "{}: {error} Da", psm.peptide.proforma());
                }
                assert!(
                    peptide_mass(&psm.peptide.sequence).is_some() || psm.calculated_mz.is_none()
                );
            }
            for protein in &batch.proteins {
                if let Some(sequence) = &protein.sequence {
                    let sequence = sequence.replace('I', "L");
                    for psm in batch.psms.iter().filter(|psm| {
                        psm.evidence
                            .iter()
                            .any(|e| e.accession == protein.accession)
                    }) {
                        assert!(sequence.contains(&psm.peptide.sequence.replace('I', "L")));
                    }
                }
            }
        }
        // Stored peaks: every identified spectrum naming a peak list reads,
        // within Scaffold's 300, sorted by m/z.
        if let Some(peaks) = &imported.peaks {
            let identified: HashSet<(String, String)> = imported
                .batches
                .iter()
                .flat_map(|batch| &batch.spectra)
                .map(|spectrum| (spectrum.ms_run.clone(), spectrum.key.canonical()))
                .collect();
            let mut counts = Vec::new();
            for run in peaks.runs() {
                for spectrum in peaks.spectra(run) {
                    let spectrum =
                        spectrum.unwrap_or_else(|error| panic!("{}: {error}", file.display()));
                    assert!(identified.contains(&(run.to_string(), spectrum.key.canonical())));
                    assert!(spectrum.mz.windows(2).all(|pair| pair[0] <= pair[1]));
                    assert!(spectrum.mz.len() <= 300, "{} peaks", spectrum.mz.len());
                    counts.push(spectrum.mz.len());
                }
            }
            assert!(
                counts.iter().all(|&n| n > 0),
                "empty peak lists are left out"
            );
            counts.sort_unstable();
            println!(
                "   peaks: {} of {} spectra, median {} peaks",
                counts.len(),
                peaks.len(),
                counts.get(counts.len() / 2).copied().unwrap_or(0),
            );
        }
    }
}
