//! The project graph: biosamples, MS runs, which spectra belong to which
//! run, search programs, the modification table, proteins, and the
//! protein → peptide map.

use std::collections::{HashMap, HashSet};

use super::blockbased::{self, BbValue};
use super::java::{ArrayData, Stream, Value};
use super::records::Container;
use super::{Sf3Error, layout};

/// Scaffold's label for a biosample without a category.
const UNCATEGORIZED: &str = "Uncategorized Sample";

pub(crate) struct Run {
    pub name: String,
    pub biosample: String,
    pub category: Option<String>,
    /// Spectrum analyses Scaffold lists for the run.
    pub analyses: Vec<i32>,
}

pub(crate) struct ModificationType {
    pub name: String,
    pub mono_mass: f64,
    /// Residues it may sit on; `n` / `c` are the peptide termini.
    pub residues: String,
}

pub(crate) struct ProteinEntry {
    pub description: Option<String>,
    pub sequence: Option<String>,
}

/// Biosamples and runs from `Experiment`, with each run's spectrum
/// analyses from its `SpectrumCollation`.
pub(crate) fn runs(container: &Container) -> Result<Vec<Run>, Sf3Error> {
    let experiment = container.stream("Experiment")?;
    let collations: HashMap<i64, Vec<i32>> = blockbased::read(container, "SpectrumCollation")?
        .into_iter()
        .filter_map(|collation| {
            let id = collation.int("id")?;
            Some((
                id,
                collation_analyses(collation.field("spectrumGroupIndex")?),
            ))
        })
        .collect();
    let samples = experiment.objects_of("TandemMassSpectrometrySample");
    if samples.is_empty() {
        return Err(layout("Experiment", "no MS runs"));
    }
    let mut runs = Vec::with_capacity(samples.len());
    for sample in samples {
        let fraction = experiment
            .str(sample, "fractionName")
            .filter(|name| !name.trim().is_empty())
            .ok_or_else(|| layout("Experiment", "an MS run has no name"))?
            .to_string();
        let bio = experiment
            .object(sample, "bioSample")
            .ok_or_else(|| layout("Experiment", &format!("run {fraction} has no biosample")))?;
        let biosample = experiment
            .str(bio, "sampleName")
            .filter(|name| !name.trim().is_empty())
            .ok_or_else(|| {
                layout(
                    "Experiment",
                    &format!("run {fraction}'s biosample has no name"),
                )
            })?
            .to_string();
        let category = experiment
            .str(bio, "category")
            .filter(|category| !category.is_empty() && *category != UNCATEGORIZED)
            .map(str::to_string);
        let collation = experiment.int(sample, "spectrumCollationId").unwrap_or(-1);
        let analyses = collations.get(&collation).cloned().ok_or_else(|| {
            layout(
                "SpectrumCollation",
                &format!("no collation {collation} for run {fraction}"),
            )
        })?;
        let expected = experiment
            .int(sample, "numSpectrumAnalyses")
            .ok_or_else(|| {
                layout(
                    "Experiment",
                    &format!("run {fraction} has no spectrum count"),
                )
            })?;
        if usize::try_from(expected) != Ok(analyses.len()) {
            return Err(layout(
                "SpectrumCollation",
                &format!(
                    "run {fraction} lists {} spectra but the project says {expected}",
                    analyses.len()
                ),
            ));
        }
        if let Some(other) = runs
            .iter()
            .find(|run: &&Run| run.biosample == biosample && run.category != category)
        {
            return Err(layout(
                "Experiment",
                &format!(
                    "two biosamples are named {biosample} ({:?} and {category:?})",
                    other.category
                ),
            ));
        }
        runs.push(Run {
            name: fraction,
            biosample,
            category,
            analyses,
        });
    }
    unique_names(&mut runs);
    Ok(runs)
}

/// The spectrum-analysis ids of a collation's `spectrumGroupIndex`: an
/// array of TreeMaps (precursor mass → TIntArrayList).
fn collation_analyses(index: &BbValue) -> Vec<i32> {
    let BbValue::Array(trees) = index else {
        return Vec::new();
    };
    let mut ids = Vec::new();
    for tree in trees {
        let BbValue::Object(tree) = tree else {
            continue;
        };
        for value in tree.custom.iter().skip(1).step_by(2) {
            if let BbValue::Object(list) = value {
                ids.extend(list.custom.iter().filter_map(|id| match id {
                    BbValue::Int(id) => i32::try_from(*id).ok(),
                    _ => None,
                }));
            }
        }
    }
    ids
}

/// Run names are natural keys: a fraction name used by several runs gets
/// its biosample's name in front (Scaffold names fractions "Ctrl" in
/// several biosamples).
fn unique_names(runs: &mut [Run]) {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for run in runs.iter() {
        *counts.entry(run.name.clone()).or_default() += 1;
    }
    for run in runs.iter_mut() {
        if counts[&run.name] > 1 {
            run.name = format!("{}: {}", run.biosample, run.name);
        }
    }
    let mut seen: HashSet<String> = HashSet::new();
    for (index, run) in runs.iter_mut().enumerate() {
        if !seen.insert(run.name.clone()) {
            run.name = format!("{} ({})", run.name, index + 1);
            seen.insert(run.name.clone());
        }
    }
}

/// Search programs by id: (name, version).
pub(crate) fn programs(
    container: &Container,
) -> Result<HashMap<u8, (String, Option<String>)>, Sf3Error> {
    let search = container.stream("DatabaseSearch")?;
    Ok(search
        .objects_of("ProgramDefinition")
        .into_iter()
        .filter_map(|program| {
            let id = u8::try_from(search.int(program, "id")?).ok()?;
            let name = search.str(program, "programName")?.trim().to_string();
            let version = search
                .str(program, "programVersion")
                .map(str::trim)
                .filter(|version| !version.is_empty() && !version.eq_ignore_ascii_case("unknown"))
                .map(str::to_string);
            Some((id, (name, version)))
        })
        .collect())
}

pub(crate) fn modification_types(
    container: &Container,
) -> Result<HashMap<i32, ModificationType>, Sf3Error> {
    let stream = container.stream("ModificationTypes")?;
    let mut types = HashMap::new();
    for object in stream.objects_of("ModificationType") {
        let (Some(id), Some(mono_mass)) =
            (stream.int(object, "id"), stream.float(object, "monoMass"))
        else {
            continue;
        };
        if !mono_mass.is_finite() {
            return Err(layout(
                "ModificationTypes",
                &format!("type {id} has no usable mass"),
            ));
        }
        let name = stream
            .str(object, "referenceName")
            .or_else(|| stream.str(object, "fullName"))
            .unwrap_or_default()
            .to_string();
        let residues = match stream.array(object, "aminoAcidsModified") {
            Some(ArrayData::Chars(residues)) => residues.clone(),
            _ => String::new(),
        };
        if let Ok(id) = i32::try_from(id) {
            types.insert(
                id,
                ModificationType {
                    name,
                    mono_mass,
                    residues,
                },
            );
        }
    }
    if types.is_empty() {
        return Err(layout("ModificationTypes", "no modification types"));
    }
    Ok(types)
}

/// Peptide sequence → accessions of the proteins Scaffold maps it to.
pub(crate) fn protein_map(container: &Container) -> Result<HashMap<String, Vec<String>>, Sf3Error> {
    let Some(stream) = container.optional_stream("IdentificationCacheLookup")? else {
        return Ok(HashMap::new());
    };
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    for object in stream.objects_of("IdentificationCacheLookup") {
        let Some(accession) = stream.str(object, "proteinAccession") else {
            continue;
        };
        for sequence in strings(&stream, stream.array(object, "sequences")) {
            let accessions = map.entry(sequence.to_string()).or_default();
            if !accessions.iter().any(|known| known == accession) {
                accessions.push(accession.to_string());
            }
        }
    }
    Ok(map)
}

fn strings<'a>(stream: &'a Stream, array: Option<&'a ArrayData>) -> Vec<&'a str> {
    match array {
        Some(ArrayData::Values(values)) => values
            .iter()
            .filter_map(|value| match value {
                Value::Str(index) => Some(stream.strings[*index].as_str()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Stored protein sequences and descriptions, by accession.
pub(crate) fn proteins(container: &Container) -> Result<HashMap<String, ProteinEntry>, Sf3Error> {
    if container.file("BBProteinSequence")?.is_none() {
        return Ok(HashMap::new());
    }
    Ok(blockbased::read(container, "BBProteinSequence")?
        .into_iter()
        .filter_map(|protein| {
            let accession = protein.str("accessionNumber")?.to_string();
            let description = protein
                .str("proteinName")
                .or_else(|| protein.str("annotation"))
                .filter(|text| !text.is_empty())
                .map(str::to_string);
            let sequence = protein
                .str("sequence")
                .filter(|sequence| !sequence.is_empty())
                .map(str::to_ascii_uppercase);
            Some((
                accession,
                ProteinEntry {
                    description,
                    sequence,
                },
            ))
        })
        .collect())
}
