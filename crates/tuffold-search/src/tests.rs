//! Searches of spectra built from known peptides' theoretical ions.

use tuffold_core::fragments::fragment_ions;
use tuffold_core::{ModSite, Modification, MsRun, PeakSpectrum, PeptideForm, SpectrumKey};
use tuffold_io::fasta::parse_fasta_str;

use super::*;

const CARBAMIDOMETHYL: f64 = 57.021_464;

/// Three proteins; the peptides below are tryptic pieces of the first two.
const FASTA: &str = ">sp|P1|ONE_TEST First protein\n\
MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK\n\
>sp|P2|TWO_TEST Second protein\n\
MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKR\n\
>sp|P3|THREE_TEST Third protein\n\
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMR\n\
>rev_sp|P3|THREE_TEST\nRMASYEEQGATDLIDLLCTEGDIVVQK\n";

fn form(sequence: &str) -> PeptideForm {
    // Carbamidomethyl on every C, as the default search fixes it.
    let modifications = sequence
        .char_indices()
        .filter(|(_, r)| *r == 'C')
        .map(|(ix, _)| Modification {
            site: ModSite::Residue(u32::try_from(ix + 1).unwrap()),
            mass_delta: CARBAMIDOMETHYL,
            name: None,
            accession: None,
        })
        .collect();
    PeptideForm::new(sequence, modifications).unwrap()
}

/// A clean spectrum of `sequence` at charge 2: every b/y ion plus a little
/// evenly spread noise.
fn spectrum(scan: u32, sequence: &str) -> PeakSpectrum {
    spectrum_of(scan, &form(sequence))
}

fn spectrum_of(scan: u32, form: &PeptideForm) -> PeakSpectrum {
    let deltas: f64 = form.modifications.iter().map(|m| m.mass_delta).sum();
    let precursor = mz(peptide_mass(&form.sequence).unwrap() + deltas, 2);
    let mut peaks: Vec<(f64, f32)> = fragment_ions(form, Some(2))
        .unwrap()
        .into_iter()
        .map(|ion| (ion.mz, 1000.0))
        .chain((0..20).map(|i| (160.0 + 61.3 * f64::from(i), 50.0)))
        .collect();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (mz, intensity) = peaks.into_iter().unzip();
    PeakSpectrum {
        key: SpectrumKey::Scan(scan),
        native_id: format!("scan={scan}"),
        ms_level: 2,
        rt_seconds: Some(60.0 * f64::from(scan)),
        precursor_mz: Some(precursor),
        precursor_charge: Some(2),
        mz,
        intensity,
        profile: false,
        parent: None,
        isolation: None,
        precursor_mobility: None,
        mobility: Vec::new(),
    }
}

fn run() -> MsRun {
    MsRun {
        name: "r1".into(),
        biosample: "s1".into(),
        spectra_file: Some("r1.mzML".into()),
    }
}

#[test]
fn theoretical_spectra_find_their_peptides() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    assert_eq!(searcher.targets(), 3);
    assert_eq!(searcher.dropped_decoys(), 1, "the FASTA's own decoy");
    assert!(searcher.peptides() > 10);

    let peptides = [
        "GLVLIAFSQYLQQCPFDEHVK",
        "AGFAGDDAPR",
        "SALTIQLIQNHFVDEYDPTIEDSYR",
    ];
    let spectra: Vec<PeakSpectrum> = peptides
        .iter()
        .zip(1..)
        .map(|(sequence, scan)| spectrum(scan, sequence))
        .chain([PeakSpectrum {
            precursor_mz: None,
            ..spectrum(9, "AGFAGDDAPR")
        }])
        .collect();
    let found = searcher.search_run(&run(), &spectra);
    assert_eq!(found.searched, 3);
    assert_eq!(
        found.skipped,
        Skipped {
            no_precursor: 1,
            ..Skipped::default()
        }
    );

    let batch = &found.batch;
    assert_eq!(batch.software.name, "Sage");
    assert_eq!(batch.ms_runs, [run()]);
    assert_eq!(batch.biosamples[0].name, "s1");
    let identified: Vec<(&SpectrumKey, &str)> = batch
        .psms
        .iter()
        .map(|psm| (&psm.spectrum, psm.peptide.sequence.as_str()))
        .collect();
    assert_eq!(
        identified,
        [
            (&SpectrumKey::Scan(1), peptides[0]),
            (&SpectrumKey::Scan(2), peptides[1]),
            (&SpectrumKey::Scan(3), peptides[2]),
        ]
    );
    let first = &batch.psms[0];
    // The fixed modification comes back with its exact mass and name.
    let cam = &first.peptide.modifications[0];
    assert_eq!(cam.site, ModSite::Residue(14));
    assert_eq!(
        (cam.mass_delta, cam.name.as_deref()),
        (CARBAMIDOMETHYL, Some("Carbamidomethyl"))
    );
    assert_eq!(first.evidence[0].accession, "sp|P1|ONE_TEST");
    let calculated = first.calculated_mz.unwrap();
    assert!((calculated - spectra[0].precursor_mz.unwrap()).abs() < 1e-9);
    let score = |name: &str| first.scores.iter().find(|s| s.name == name).unwrap().value;
    assert!(score("hyperscore") > 0.0 && score("matched_peaks") >= 10.0);

    let protein = batch
        .proteins
        .iter()
        .find(|p| p.accession == "sp|P1|ONE_TEST")
        .unwrap();
    assert_eq!(protein.description.as_deref(), Some("First protein"));
    assert!(protein.sequence.is_some() && !protein.is_decoy);
    assert_eq!(batch.spectra.len(), 3);
}

#[test]
fn decoy_matches_name_decoy_proteins() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    // A reversed tryptic peptide of P2 (APADDGAFGR, the decoy of
    // AGFAGDDAPR) is in the index; its spectrum matches it as a decoy.
    let found = searcher
        .search_run(&run(), &[spectrum(1, "APADDGAFGR")])
        .batch;
    let psm = &found.psms[0];
    assert_eq!(psm.peptide.sequence, "APADDGAFGR");
    assert_eq!(psm.evidence[0].accession, "rev_sp|P2|TWO_TEST");
    let protein = &found.proteins[0];
    assert!(protein.is_decoy && protein.sequence.is_none());
    assert!(found.psm_is_decoy(psm));
}

#[test]
fn searches_are_deterministic() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    let spectra: Vec<PeakSpectrum> = (1..=6)
        .map(|scan| {
            spectrum(
                scan,
                ["AGFAGDDAPR", "GLVLIAFSQYLQQCPFDEHVK"][scan as usize % 2],
            )
        })
        .collect();
    let a = searcher.search_run(&run(), &spectra).batch;
    let b = searcher.search_run(&run(), &spectra).batch;
    assert_eq!(a, b);
}

#[test]
fn unsearchable_spectra_are_counted_by_reason() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    let good = spectrum(1, "AGFAGDDAPR");
    let spectra = [
        good.clone(),
        // The same scan again: the first is searched, this one counted.
        spectrum(1, "GLVLIAFSQYLQQCPFDEHVK"),
        PeakSpectrum {
            key: SpectrumKey::Scan(3),
            precursor_mz: Some(f64::NAN),
            ..good.clone()
        },
        PeakSpectrum {
            key: SpectrumKey::Scan(4),
            mz: good.mz[..MIN_PEAKS - 1].to_vec(),
            intensity: good.intensity[..MIN_PEAKS - 1].to_vec(),
            ..good.clone()
        },
        PeakSpectrum {
            key: SpectrumKey::Scan(5),
            ms_level: 1,
            ..good.clone()
        },
    ];
    let found = searcher.search_run(&run(), &spectra);
    let expected = Skipped {
        no_precursor: 1,
        too_few_peaks: 1,
        repeated_key: 1,
    };
    assert_eq!((found.searched, found.skipped), (1, expected));
    assert_eq!(found.batch.psms.len(), 1);
    assert_eq!(found.batch.psms[0].peptide.sequence, "AGFAGDDAPR");
    let described = expected.describe().unwrap();
    assert!(described.starts_with("1 without a precursor m/z, 1 with fewer than 15 peaks"));
    assert_eq!(Skipped::default().describe(), None);
}

#[test]
fn targets_sage_would_misread_are_reported() {
    let text = format!("{FASTA}>contam_rev_X\nMKAGFAGDDAPRK\n>sp|P1|ONE_TEST again\nMKK\n");
    let searcher = Searcher::new(&parse_fasta_str(&text), SearchSettings::default()).unwrap();
    assert_eq!(searcher.targets(), 3);
    assert_eq!(
        searcher.excluded_targets(),
        ["contam_rev_X", "sp|P1|ONE_TEST"]
    );
}

#[test]
fn bad_inputs_are_refused() {
    let decoys_only = parse_fasta_str(">rev_P1\nPEPTIDEK\n");
    assert!(matches!(
        Searcher::new(&decoys_only, SearchSettings::default()),
        Err(SearchError::NoTargets)
    ));
    let records = parse_fasta_str(FASTA);
    for settings in [
        SearchSettings {
            precursor_ppm: 0.0,
            ..SearchSettings::default()
        },
        SearchSettings {
            min_length: 9,
            max_length: 8,
            ..SearchSettings::default()
        },
        SearchSettings {
            charges: (3, 2),
            ..SearchSettings::default()
        },
        SearchSettings {
            max_variable_mods: 0,
            ..SearchSettings::default()
        },
    ] {
        assert!(matches!(
            Searcher::new(&records, settings),
            Err(SearchError::Settings(_))
        ));
    }
}

#[test]
fn tmt_labelled_peptides_are_found_with_the_label() {
    // DLGEEHFK from P1, TMT6plex on the N-terminus and on K8.
    let tmt = |site| Modification {
        site,
        mass_delta: 229.162_932,
        name: None,
        accession: None,
    };
    let labelled = PeptideForm::new(
        "DLGEEHFK",
        vec![tmt(ModSite::NTerm), tmt(ModSite::Residue(8))],
    )
    .unwrap();
    let records = parse_fasta_str(FASTA);
    let settings = SearchSettings::default().with_label(Label::Tmt);
    let searcher = Searcher::new(&records, settings).unwrap();
    let found = searcher
        .search_run(&run(), &[spectrum_of(1, &labelled)])
        .batch;
    let psm = &found.psms[0];
    assert_eq!(psm.peptide.sequence, "DLGEEHFK");
    let mods: Vec<(&ModSite, Option<&str>)> = psm
        .peptide
        .modifications
        .iter()
        .map(|m| (&m.site, m.name.as_deref()))
        .collect();
    assert_eq!(
        mods,
        [
            (&ModSite::NTerm, Some("TMT6plex")),
            (&ModSite::Residue(8), Some("TMT6plex")),
        ]
    );
    // Without the label the same spectrum is not that peptide.
    let plain = Searcher::new(&records, SearchSettings::default()).unwrap();
    let unlabelled = plain.search_run(&run(), &[spectrum_of(1, &labelled)]).batch;
    assert!(unlabelled.psms.iter().all(|p| p.peptide != psm.peptide));
}

#[test]
fn the_enzyme_decides_which_peptides_exist() {
    // WVTFISLLLLFSSAYSRGVFRRDTHK: Lys-C cuts after K only. Trypsin would
    // have to miss three cleavages (after SR, FR, and the first R of RR).
    let lys_c = "WVTFISLLLLFSSAYSRGVFRRDTHK";
    let records = parse_fasta_str(FASTA);
    let search = |enzyme| {
        let settings = SearchSettings {
            enzyme,
            ..SearchSettings::default()
        };
        let searcher = Searcher::new(&records, settings).unwrap();
        searcher.search_run(&run(), &[spectrum(1, lys_c)]).batch
    };
    let found = search(Enzyme::LysC);
    assert_eq!(found.psms[0].peptide.sequence, lys_c);
    let tryptic = search(Enzyme::Trypsin);
    assert!(tryptic.psms.iter().all(|p| p.peptide.sequence != lys_c));
}

/// `spectrum` as an instrument would record it in profile mode: every
/// peak a Gaussian (sigma 0.004) sampled every 0.002 m/z.
fn as_profile(spectrum: &PeakSpectrum) -> PeakSpectrum {
    let mut points: Vec<(f64, f32)> = spectrum
        .mz
        .iter()
        .zip(&spectrum.intensity)
        .flat_map(|(&center, &height)| {
            // Sample on a shared grid so neighbouring peaks don't interleave.
            let first = ((center - 0.02) / 0.002).ceil() as i64;
            (first..first + 20).map(move |k| {
                let x = k as f64 * 0.002;
                let y = f64::from(height) * (-((x - center) / 0.004).powi(2) / 2.0).exp();
                (x, y as f32)
            })
        })
        .collect();
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Peaks closer than the window overlap on the grid: sum them.
    points.dedup_by(|later, kept| {
        let same = (later.0 - kept.0).abs() < 1e-9;
        if same {
            kept.1 += later.1;
        }
        same
    });
    let (mz, intensity) = points.into_iter().unzip();
    PeakSpectrum {
        mz,
        intensity,
        profile: true,
        ..spectrum.clone()
    }
}

#[test]
fn profile_spectra_are_centroided_then_searched() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    let profile = as_profile(&spectrum(1, "GLVLIAFSQYLQQCPFDEHVK"));
    assert!(profile.mz.len() > 20 * 30, "{} points", profile.mz.len());
    let found = searcher.search_run(&run(), &[profile]);
    assert_eq!((found.searched, found.centroided), (1, 1));
    assert_eq!(
        found.batch.psms[0].peptide.sequence,
        "GLVLIAFSQYLQQCPFDEHVK"
    );
    let matched = found.batch.psms[0]
        .scores
        .iter()
        .find(|s| s.name == "matched_peaks")
        .unwrap()
        .value;
    assert!(
        matched >= 20.0,
        "{matched} fragments matched after centroiding"
    );
}

#[test]
fn the_index_estimate_bounds_what_is_built() {
    let records = parse_fasta_str(FASTA);
    let searcher = Searcher::new(&records, SearchSettings::default()).unwrap();
    let (estimate, budget) = searcher.index_estimate();
    assert!(estimate.complete);
    // Duplicates and decoys equal to a target are merged after counting, so
    // the count is an upper bound of what the index holds.
    assert!(estimate.forms >= searcher.peptides() as u64);
    assert!(estimate.fragments > estimate.forms);
    assert!(0 < estimate.peak_bytes && estimate.peak_bytes <= budget);

    let semi = SearchSettings {
        semi_specific: true,
        ..SearchSettings::default()
    };
    let (semi_estimate, _) = Searcher::new(&records, semi).unwrap().index_estimate();
    assert!(semi_estimate.forms > estimate.forms && semi_estimate.peak_bytes > estimate.peak_bytes);
}

#[test]
fn searches_that_would_not_fit_are_refused_before_building() {
    let records = parse_fasta_str(FASTA);
    let tight = SearchSettings {
        semi_specific: true,
        memory_budget: Some(1_000),
        ..SearchSettings::default()
    };
    let Err(error) = Searcher::new(&records, tight) else {
        panic!("a 1 kB budget fits no index");
    };
    let SearchError::TooLarge {
        estimate, budget, ..
    } = &error
    else {
        panic!("{error}");
    };
    assert_eq!(*budget, 1_000);
    // Counting stopped a few multiples past the budget.
    assert!(!estimate.complete && estimate.peak_bytes > 4 * 1_000);
    let message = error.to_string();
    assert!(message.contains("fully specific"), "{message}");
    assert!(message.contains("budget is 0.0 GB"), "{message}");

    // The same search with room for it is built.
    let roomy = SearchSettings {
        memory_budget: Some(10_000_000_000),
        ..SearchSettings::default()
    };
    assert!(Searcher::new(&records, roomy).is_ok());
}
