//! The demo project: generated, imported, validated, and browsable.

use tuffold_engine::create_demo_project;
use tuffold_store::{Project, StoredPeaks};

#[test]
fn the_demo_project_validates_and_finds_the_expressed_proteins() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    let summary = create_demo_project(&path).unwrap();
    assert_eq!(summary.engines.len(), 1);
    assert_eq!(
        summary.engines[0].method, "rescored",
        "{}",
        summary.engines[0].detail
    );

    let project = Project::open(&path).unwrap();
    let list = project.protein_list().unwrap();
    let samples: Vec<&str> = list.biosamples.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(
        samples,
        ["control-A", "control-B", "treated-A", "treated-B"]
    );
    // Expressed proteins are DEMO0000..DEMO0149; anything else passing is a
    // false discovery (or an isoform sharing a parent's peptides).
    let expressed = |accession: &str| {
        accession[7..11]
            .parse::<usize>()
            .is_ok_and(|index| index < 150)
    };
    let wrong = list
        .rows
        .iter()
        .filter(|row| !expressed(&row.accession))
        .count();
    assert!(list.rows.len() >= 100, "{}", list.rows.len());
    assert!(wrong * 100 <= list.rows.len() * 3);
    assert!(list.rows.iter().all(|row| row.coverage.is_some()));

    // DEMO0000..0029 are 4x up in treated and DEMO0030..0049 4x down; summed
    // over each set, spectral counts follow (single proteins are too sparse).
    let totals = |range: std::ops::Range<usize>| {
        list.rows
            .iter()
            .filter(|row| range.contains(&row.accession[7..11].parse::<usize>().unwrap_or(999)))
            .fold((0, 0), |(control, treated), row| {
                let counts = &row.spectra_by_sample;
                (
                    control + counts[0] + counts[1],
                    treated + counts[2] + counts[3],
                )
            })
    };
    let (control, treated) = totals(0..30);
    assert!(
        treated > 2 * control,
        "up: control {control}, treated {treated}"
    );
    let (control, treated) = totals(30..50);
    assert!(
        control > 2 * treated,
        "down: control {control}, treated {treated}"
    );

    let detail = project
        .group_detail(list.rows[0].group_id)
        .unwrap()
        .unwrap();
    let spectra = project
        .peptide_spectra(detail.peptides[0].peptide_id, 50)
        .unwrap()
        .rows;
    let spectrum = project
        .spectrum_detail(spectra[0].spectrum_id)
        .unwrap()
        .unwrap();
    assert!(matches!(spectrum.peaks, StoredPeaks::Stored(_)));
    assert!(spectrum.assigned.is_some());
}

#[test]
fn the_demo_quantifies_and_finds_its_regulated_proteins() {
    use tuffold_engine::{quantify_project, write_quant_tsv};
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();
    let quantities = quantify_project(&project).unwrap();
    assert_eq!(quantities.categories(), ["control", "treated"]);
    let comparison = quantities.compare("control", "treated").unwrap();
    // DEMO0000..0029 are 4x up and DEMO0030..0049 4x down in treated. With
    // ~5,000 spectra in all, a handful pass 5% FDR (Fisher), nearly all among
    // those 50 with the right direction. (An unchanged protein can look down:
    // the up-regulated ones take a larger share of treated spectra, the
    // compositional effect of total-count normalization.)
    let index = |accession: &str| accession[7..11].parse::<usize>().unwrap_or(999);
    let significant: Vec<_> = comparison
        .rows
        .iter()
        .filter(|r| r.fisher_q <= 0.05)
        .collect();
    assert!(significant.len() >= 5, "{} significant", significant.len());
    let right = significant
        .iter()
        .filter(|r| match index(&r.accession) {
            0..30 => r.log2_fold > 0.0,
            30..50 => r.log2_fold < 0.0,
            _ => false,
        })
        .count();
    assert!(
        right * 4 >= significant.len() * 3,
        "{right} of {}",
        significant.len()
    );
    // Two samples per category: Welch tests run.
    assert!(comparison.rows.iter().any(|r| r.welch_p.is_some()));
    // Mann–Whitney runs too, but 2 vs 2 can't go below p = 2/C(4,2) = 1/3:
    // complete separation is the best it sees.
    let lowest = comparison
        .rows
        .iter()
        .filter_map(|r| r.mann_whitney_p)
        .fold(1.0, f64::min);
    assert!((lowest - 1.0 / 3.0).abs() < 1e-12, "{lowest}");

    let mut tsv = Vec::new();
    write_quant_tsv(&mut tsv, &quantities, Some(&comparison)).unwrap();
    let text = String::from_utf8(tsv).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(header[0], "accession");
    assert!(header.contains(&"log2_fold:treated/control"));
    assert_eq!(lines.count(), quantities.rows.len());
}

#[test]
fn the_demo_measures_ms1_intensity_and_recovers_its_fold_changes() {
    use tuffold_engine::{quantify_intensity, write_intensity_tsv};
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();
    let quant = quantify_intensity(&project).unwrap();

    // Every run's MS1 scans were read and its precursors traced. About half
    // the identifications are wrong, with no signal at their m/z.
    for run in &quant.runs {
        assert_eq!(run.ms1_scans, Some(1200), "{run:?}");
        assert_eq!(run.untraced, 0);
        assert!(
            run.measured * 3 > run.traced && run.measured < run.traced,
            "{run:?}"
        );
    }
    // Loading factors 1.0, 1.4, 0.7, 1.1: normalization shifts each run by
    // log2(geometric mean / factor) = 0.03, −0.46, 0.54, −0.11.
    for (shift, expected) in quant
        .quantities
        .shifts
        .iter()
        .zip([0.027, -0.458, 0.542, -0.11])
    {
        assert!((shift - expected).abs() < 0.15, "{shift} vs {expected}");
    }

    let comparison = quant.quantities.compare("control", "treated").unwrap();
    let index = |accession: &str| accession[7..11].parse::<usize>().unwrap_or(999);
    let median = |mut values: Vec<f64>| {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    };
    // DEMO0000..0029 are 4× up (log2 +2), DEMO0030..0049 4× down.
    let folds = |range: std::ops::Range<usize>| -> Vec<f64> {
        comparison
            .rows
            .iter()
            .filter(|r| range.contains(&index(&r.accession)))
            .filter_map(|r| r.log2_fold)
            .collect()
    };
    let (up, down, unchanged) = (
        median(folds(0..30)),
        median(folds(30..50)),
        median(folds(50..150)),
    );
    assert!((up - 2.0).abs() < 0.3, "up {up}");
    assert!((down + 2.0).abs() < 0.3, "down {down}");
    assert!(unchanged.abs() < 0.25, "unchanged {unchanged}");
    let significant: Vec<_> = comparison
        .rows
        .iter()
        .filter(|r| r.q.is_some_and(|q| q <= 0.05))
        .collect();
    assert!(significant.len() >= 5, "{} significant", significant.len());
    // Planted changes are called in their direction; at 5% FDR an unchanged
    // group may slip in (1 of 11 here: DEMO0066 at q 0.0496).
    assert!(significant.iter().all(|r| match index(&r.accession) {
        0..30 => r.log2_fold > Some(0.0),
        30..50 => r.log2_fold < Some(0.0),
        _ => true,
    }));
    let unchanged_called = significant
        .iter()
        .filter(|r| index(&r.accession) >= 50)
        .count();
    assert!(
        unchanged_called * 5 <= significant.len(),
        "{unchanged_called} of {} unchanged",
        significant.len()
    );

    // Every category at once: with two categories, ANOVA's calls are mostly
    // the planted ones too. Its pooled variance also tests groups with three
    // values (2 + 1), which Welch can't, so it calls a few more; at 5% FDR
    // an unchanged one may slip in (1 of 13 here).
    let anova = quant.quantities.anova().unwrap();
    assert_eq!(anova.categories, ["control", "treated"]);
    let called: Vec<_> = anova
        .rows
        .iter()
        .filter(|r| r.q.is_some_and(|q| q <= 0.05))
        .collect();
    let wrong = called.iter().filter(|r| index(&r.accession) >= 50).count();
    assert!(called.len() >= 5, "{} called", called.len());
    assert!(
        wrong * 5 <= called.len(),
        "{wrong} of {} wrong",
        called.len()
    );
    let mut table = Vec::new();
    tuffold_engine::write_anova_tsv(&mut table, &anova).unwrap();
    let table = String::from_utf8(table).unwrap();
    assert!(table.starts_with("accession\tmean_log2:control\tmean_log2:treated\t"));
    assert_eq!(table.lines().count(), anova.rows.len() + 1);

    let mut tsv = Vec::new();
    write_intensity_tsv(&mut tsv, &quant.quantities, Some(&comparison)).unwrap();
    let text = String::from_utf8(tsv).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(header[0], "accession");
    assert!(header.contains(&"intensity:control-A") && header.contains(&"welch_q"));
    assert_eq!(lines.count(), quant.quantities.rows.len());
}

#[test]
fn the_demo_exports_a_workbook_and_methods() {
    use tuffold_engine::{methods, write_localization_tsv, write_workbook};
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();

    let text = methods(&project).unwrap().text();
    // Every number is the project's: 4 runs, the default 1% FDR thresholds.
    assert!(text.contains("from 4 MS runs (4 biosamples)"), "{text}");
    assert!(text.starts_with("THESE DATA ARE SYNTHETIC"), "{text}");
    assert!(
        text.contains("Thresholds were 1% FDR for spectra"),
        "{text}"
    );
    assert!(
        text.contains("Tukey's median polish") && text.contains("half-sample mode"),
        "the demo has MS1 data: {text}"
    );
    assert!(
        text.contains(". Tukey JW. Exploratory Data Analysis"),
        "median polish is cited: {text}"
    );
    assert!(
        text.contains("MaxLFQ ["),
        "the normalization idea is cited: {text}"
    );
    // References are numbered in the order they are first cited.
    assert!(
        text.contains("[1]") && text.contains("\n1. Elias JE"),
        "{text}"
    );
    assert!(!text.contains("Sage"), "the demo wasn't searched: {text}");
    // Peptide pooling is cited exactly when it scored the peptides.
    let validation = project.validation_summary().unwrap().unwrap();
    let scoring = validation.peptide_scoring.unwrap();
    assert_eq!(
        text.contains("Peptides were scored from all their spectra"),
        scoring.method == "pooled",
        "{scoring:?}: {text}"
    );
    // Oxidized peptides with two or more methionines are localized and cited.
    assert!(
        validation.localization.scored_psms > 0 && validation.localization.localized_psms > 0,
        "{:?}",
        validation.localization
    );
    assert!(
        text.contains("Variable modifications with more than one candidate residue")
            && text.contains("Taus T")
            && text.contains("Beausoleil SA"),
        "{text}"
    );
    // The decoy-amino-acid reference is cited only with an FLR to report.
    assert_eq!(
        text.contains("Ramsbottom KA"),
        validation.localization.estimated_flr.is_some(),
        "{text}"
    );

    let mut loc_tsv = Vec::new();
    let loc_rows = write_localization_tsv(&mut loc_tsv, &project).unwrap();
    assert_eq!(loc_rows, validation.localization.scored_psms);
    let loc_text = String::from_utf8(loc_tsv).unwrap();
    assert!(loc_text.starts_with(
        "run\tbiosample\tspectrum\tcharge\tprecursor_mz\trt_seconds\tpeptide\tlocalized_peptide\tsite_probabilities\tbest_ascore\tq_value\tprobability\n"
    ));

    let workbook = dir.path().join("demo.xlsx");
    let written = write_workbook(&project, &workbook).unwrap();
    let report = project.report().unwrap();
    assert!(
        report
            .peptides
            .iter()
            .any(|p| p.localized_peptide.is_some() && p.site_probabilities.is_some())
    );
    let sheets: Vec<(&str, usize)> = written
        .sheets
        .iter()
        .map(|(name, rows)| (name.as_str(), *rows))
        .collect();
    assert_eq!(
        sheets[..3],
        [
            ("Proteins", project.protein_list().unwrap().rows.len()),
            ("Peptides", report.peptides.len()),
            ("Spectra", report.spectra.len()),
        ]
    );
    let names: Vec<&str> = sheets.iter().map(|(name, _)| *name).collect();
    assert_eq!(names[3..], ["ANOVA (counts)", "ANOVA (MS1)"]);
    assert!(written.cut.is_empty());
    let bytes = std::fs::read(&workbook).unwrap();
    assert_eq!(&bytes[..2], b"PK", "an xlsx is a zip archive");
}
