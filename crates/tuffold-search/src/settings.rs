//! What to search for: enzyme, modifications, tolerances. Parsed from the
//! short forms people type (`--enzyme lys-c`, `K:TMT6plex`, `^:229.1629`) so
//! the CLI, the MCP server, and the app all read them the same way.

use std::fmt;

use sage_core::database::Builder;
use sage_core::mass::Tolerance;
use tuffold_core::unimod;

use crate::{DECOY_TAG, SearchError};

/// Fragment m/z range indexed and scored (Sage's defaults).
pub(crate) const FRAGMENT_MZ: (f32, f32) = (150.0, 2000.0);
/// Peptide neutral masses digested (Sage's defaults).
const PEPTIDE_MASS: (f64, f64) = (500.0, 5000.0);

/// The protease, with the cleavage rules search engines conventionally use
/// (MSFragger's table).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Enzyme {
    /// After K or R, not before P.
    Trypsin,
    /// After K or R, before P too.
    TrypsinP,
    /// After K, not before P.
    LysC,
    /// Before K.
    LysN,
    /// After R, not before P.
    ArgC,
    /// After D or E, not before P.
    GluC,
    /// Before D.
    AspN,
    /// After F, W, Y, or L, not before P.
    Chymotrypsin,
}

impl Enzyme {
    pub const ALL: [Enzyme; 8] = [
        Self::Trypsin,
        Self::TrypsinP,
        Self::LysC,
        Self::LysN,
        Self::ArgC,
        Self::GluC,
        Self::AspN,
        Self::Chymotrypsin,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Self::Trypsin => "trypsin",
            Self::TrypsinP => "trypsin/p",
            Self::LysC => "lys-c",
            Self::LysN => "lys-n",
            Self::ArgC => "arg-c",
            Self::GluC => "glu-c",
            Self::AspN => "asp-n",
            Self::Chymotrypsin => "chymotrypsin",
        }
    }

    /// `trypsin`, `Lys-C`, `lysc`, `trypsin/p`: case and dashes don't matter.
    pub fn parse(text: &str) -> Result<Self, String> {
        let squash = |s: &str| s.to_lowercase().replace(['-', '_', ' '], "");
        Self::ALL
            .into_iter()
            .find(|enzyme| squash(enzyme.name()) == squash(text))
            .ok_or_else(|| {
                let names: Vec<&str> = Self::ALL.iter().map(|e| e.name()).collect();
                format!("unknown enzyme '{text}' (known: {})", names.join(", "))
            })
    }

    /// Sage's (cleave_at, restrict, c_terminal).
    fn rule(self) -> (&'static str, Option<char>, bool) {
        match self {
            Self::Trypsin => ("KR", Some('P'), true),
            Self::TrypsinP => ("KR", None, true),
            Self::LysC => ("K", Some('P'), true),
            Self::LysN => ("K", None, false),
            Self::ArgC => ("R", Some('P'), true),
            Self::GluC => ("DE", Some('P'), true),
            Self::AspN => ("D", None, false),
            Self::Chymotrypsin => ("FWYL", Some('P'), true),
        }
    }
}

/// Where a modification can sit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Site {
    Residue(char),
    PeptideN,
    PeptideC,
    ProteinN,
    ProteinC,
}

impl Site {
    /// `C`, or `^` / `$` for a peptide's N / C terminus, `[` / `]` for a
    /// protein's (Sage's notation).
    fn parse(text: &str) -> Result<Self, String> {
        match text {
            "^" | "n-term" | "nterm" => Ok(Self::PeptideN),
            "$" | "c-term" | "cterm" => Ok(Self::PeptideC),
            "[" | "protein-n-term" => Ok(Self::ProteinN),
            "]" | "protein-c-term" => Ok(Self::ProteinC),
            _ => {
                let mut chars = text.chars();
                match (chars.next(), chars.next()) {
                    (Some(residue), None) if "ACDEFGHIKLMNPQRSTVWY".contains(residue) => {
                        Ok(Self::Residue(residue))
                    }
                    _ => Err(format!(
                        "'{text}' is not a site: give one residue letter, or ^ $ for a \
                         peptide's N/C terminus, [ ] for a protein's"
                    )),
                }
            }
        }
    }

    fn sage(self) -> String {
        match self {
            Self::Residue(residue) => residue.to_string(),
            Self::PeptideN => "^".into(),
            Self::PeptideC => "$".into(),
            Self::ProteinN => "[".into(),
            Self::ProteinC => "]".into(),
        }
    }

    /// Terminal sites: Sage records these on the peptide's ends.
    pub(crate) fn is_n_terminal(self) -> bool {
        matches!(self, Self::PeptideN | Self::ProteinN)
    }

    pub(crate) fn is_c_terminal(self) -> bool {
        matches!(self, Self::PeptideC | Self::ProteinC)
    }
}

/// One modification to search: where, how heavy, and its Unimod identity
/// when it has one.
#[derive(Debug, Clone, PartialEq)]
pub struct ModRule {
    pub site: Site,
    pub mass: f64,
    pub name: Option<String>,
    pub accession: Option<String>,
}

impl ModRule {
    /// `SITE:MOD`, where MOD is a Unimod name from Tuffold's table
    /// (`K:TMT6plex`, `M:Oxidation`) or a mass delta (`C:57.021464`).
    pub fn parse(text: &str) -> Result<Self, String> {
        let (site, modification) = text
            .split_once(':')
            .ok_or_else(|| format!("'{text}': write SITE:MOD, e.g. K:TMT6plex or C:57.021464"))?;
        let site = Site::parse(site.trim())?;
        let modification = modification.trim();
        let known = match modification.parse::<f64>() {
            Ok(mass) if mass.is_finite() && mass != 0.0 => {
                return Ok(Self::mass(site, mass));
            }
            Ok(_) => return Err(format!("'{text}': the mass must be a non-zero number")),
            Err(_) => unimod::by_name(modification).ok_or_else(|| {
                format!("'{text}': unknown modification '{modification}'; give its mass instead")
            })?,
        };
        Ok(Self::known(site, known))
    }

    /// A mass at a site, named when the mass is in the Unimod table.
    pub fn mass(site: Site, mass: f64) -> Self {
        match unimod::by_mass(mass) {
            Some(known) if (known.mass - mass).abs() < 1e-6 => Self::known(site, known),
            _ => Self {
                site,
                mass,
                name: None,
                accession: None,
            },
        }
    }

    fn known(site: Site, known: unimod::Known) -> Self {
        Self {
            site,
            mass: known.mass,
            name: Some(known.name.to_string()),
            accession: Some(known.accession()),
        }
    }
}

impl fmt::Display for ModRule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.name {
            Some(name) => write!(f, "{}:{name}", self.site.sage()),
            None => write!(f, "{}:{}", self.site.sage(), self.mass),
        }
    }
}

/// Isobaric labels: fixed on K and every peptide N-terminus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Label {
    Tmt,
    TmtPro,
    Itraq4,
    Itraq8,
}

impl Label {
    /// `tmt` (any TMT 6/10/11-plex), `tmtpro` (16/18-plex), `itraq4`, or
    /// `itraq8`.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text.to_lowercase().replace(['-', ' '], "").as_str() {
            "tmt" | "tmt6plex" | "tmt10plex" | "tmt11plex" => Ok(Self::Tmt),
            "tmtpro" | "tmt16plex" | "tmt18plex" => Ok(Self::TmtPro),
            "itraq4" | "itraq4plex" => Ok(Self::Itraq4),
            "itraq8" | "itraq8plex" => Ok(Self::Itraq8),
            _ => Err(format!(
                "unknown label '{text}' (known: tmt, tmtpro, itraq4, itraq8)"
            )),
        }
    }

    pub fn rules(self) -> Vec<ModRule> {
        let name = match self {
            Self::Tmt => "TMT6plex",
            Self::TmtPro => "TMTpro",
            Self::Itraq4 => "iTRAQ4plex",
            Self::Itraq8 => "iTRAQ8plex",
        };
        [Site::Residue('K'), Site::PeptideN]
            .into_iter()
            .filter_map(|site| unimod::by_name(name).map(|known| ModRule::known(site, known)))
            .collect()
    }
}

/// How far a fragment peak may be from a theoretical fragment's m/z, either
/// way.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FragmentTolerance {
    /// Parts per million of the fragment's m/z (high-resolution MS/MS).
    Ppm(f32),
    /// Daltons at every m/z (low-resolution MS/MS, such as an ion trap's,
    /// whose peak widths don't grow with m/z).
    Da(f32),
}

impl FragmentTolerance {
    fn value(self) -> f32 {
        match self {
            Self::Ppm(value) | Self::Da(value) => value,
        }
    }

    /// Sage's symmetric window.
    pub(crate) fn sage(self) -> Tolerance {
        match self {
            Self::Ppm(ppm) => Tolerance::Ppm(-ppm, ppm),
            Self::Da(da) => Tolerance::Da(-da, da),
        }
    }
}

impl fmt::Display for FragmentTolerance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ppm(ppm) => write!(f, "{ppm} ppm"),
            Self::Da(da) => write!(f, "{da} Da"),
        }
    }
}

/// Search choices as typed by a person or sent by an assistant: everything
/// optional, defaults filled in by [`SearchSettings::from_options`]. The CLI
/// and the MCP server both build these, so they read options the same way.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SearchOptions {
    pub precursor_ppm: Option<f32>,
    pub fragment_ppm: Option<f32>,
    /// Fragment tolerance in Da instead of ppm; giving both is an error.
    pub fragment_da: Option<f32>,
    pub enzyme: Option<String>,
    pub semi_specific: bool,
    pub missed_cleavages: Option<u8>,
    /// `SITE:MOD` rules added to the defaults.
    pub fixed: Vec<String>,
    pub variable: Vec<String>,
    pub label: Option<String>,
    pub max_variable_mods: Option<usize>,
    /// Drop the default fixed carbamidomethyl C.
    pub no_carbamidomethyl: bool,
    /// Drop the default variable oxidized M.
    pub no_oxidation: bool,
    /// Most memory the index may take, in GB (10^9 bytes); `None` sizes it
    /// from the memory free when the search starts.
    pub max_memory_gb: Option<f64>,
}

/// What to search for.
#[derive(Debug, Clone, PartialEq)]
pub struct SearchSettings {
    /// Precursor tolerance, ± ppm.
    pub precursor_ppm: f32,
    /// Fragment tolerance, ± ppm or Da.
    pub fragment: FragmentTolerance,
    pub enzyme: Enzyme,
    /// One end may be non-enzymatic.
    pub semi_specific: bool,
    pub missed_cleavages: u8,
    pub min_length: usize,
    pub max_length: usize,
    /// On every matching site.
    pub fixed: Vec<ModRule>,
    /// On any matching site, up to `max_variable_mods` per peptide.
    pub variable: Vec<ModRule>,
    pub max_variable_mods: usize,
    /// Charges tried when a spectrum does not give its precursor's.
    pub charges: (u8, u8),
    /// Most bytes building the index may take; `None` uses
    /// [`crate::default_memory_budget`] when the search starts.
    pub memory_budget: Option<u64>,
}

impl Default for SearchSettings {
    /// Tryptic, 2 missed cleavages, carbamidomethyl C fixed, oxidized M.
    fn default() -> Self {
        let rule = |site, name| ModRule::parse(&format!("{site}:{name}")).ok();
        Self {
            precursor_ppm: 20.0,
            fragment: FragmentTolerance::Ppm(20.0),
            enzyme: Enzyme::Trypsin,
            semi_specific: false,
            missed_cleavages: 2,
            min_length: 7,
            max_length: 50,
            fixed: rule('C', "Carbamidomethyl").into_iter().collect(),
            variable: rule('M', "Oxidation").into_iter().collect(),
            max_variable_mods: 2,
            charges: (2, 4),
            memory_budget: None,
        }
    }
}

/// Bytes in `gb` gigabytes (10^9 bytes), for a memory budget.
fn gigabytes(gb: f64) -> Result<u64, String> {
    if gb.is_finite() && gb > 0.0 {
        Ok((gb * 1e9) as u64)
    } else {
        Err(format!(
            "the memory budget must be a positive number of GB, not {gb}"
        ))
    }
}

impl SearchSettings {
    /// The defaults, changed as `options` asks. Every rule is checked here, so
    /// a typo is an error, never a silently different search.
    pub fn from_options(options: &SearchOptions) -> Result<Self, String> {
        let defaults = Self::default();
        let parse_all = |texts: &[String]| -> Result<Vec<ModRule>, String> {
            texts.iter().map(|text| ModRule::parse(text)).collect()
        };
        let keep = |rules: Vec<ModRule>, drop: bool| if drop { Vec::new() } else { rules };
        let fragment = match (options.fragment_ppm, options.fragment_da) {
            (Some(_), Some(_)) => {
                return Err("give the fragment tolerance in ppm or in Da, not both".to_string());
            }
            (Some(ppm), None) => FragmentTolerance::Ppm(ppm),
            (None, Some(da)) => FragmentTolerance::Da(da),
            (None, None) => defaults.fragment,
        };
        let mut settings = Self {
            precursor_ppm: options.precursor_ppm.unwrap_or(defaults.precursor_ppm),
            fragment,
            enzyme: options
                .enzyme
                .as_deref()
                .map_or(Ok(defaults.enzyme), Enzyme::parse)?,
            semi_specific: options.semi_specific,
            missed_cleavages: options
                .missed_cleavages
                .unwrap_or(defaults.missed_cleavages),
            fixed: keep(defaults.fixed, options.no_carbamidomethyl),
            variable: keep(defaults.variable, options.no_oxidation),
            max_variable_mods: options
                .max_variable_mods
                .unwrap_or(defaults.max_variable_mods),
            memory_budget: options.max_memory_gb.map(gigabytes).transpose()?,
            ..defaults
        };
        settings.fixed.extend(parse_all(&options.fixed)?);
        settings.variable.extend(parse_all(&options.variable)?);
        if let Some(label) = &options.label {
            settings = settings.with_label(Label::parse(label)?);
        }
        settings.check().map_err(|error| error.to_string())?;
        Ok(settings)
    }

    /// One line for reports: `trypsin, up to 2 missed, precursor ±20 ppm,
    /// fragment ±20 ppm, fixed C:Carbamidomethyl, ...` (or `fragment ±0.5 Da`).
    pub fn describe(&self) -> String {
        let list = |rules: &[ModRule]| {
            if rules.is_empty() {
                "none".to_string()
            } else {
                rules
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(" ")
            }
        };
        format!(
            "{}{}, up to {} missed, precursor ±{} ppm, fragment ±{}, fixed {}, variable {}",
            self.enzyme.name(),
            if self.semi_specific { " (semi)" } else { "" },
            self.missed_cleavages,
            self.precursor_ppm,
            self.fragment,
            list(&self.fixed),
            list(&self.variable),
        )
    }

    /// Add one fixed modification.
    pub fn with_fixed(mut self, rule: ModRule) -> Self {
        self.fixed.push(rule);
        self
    }

    /// Add an isobaric label's fixed modifications.
    pub fn with_label(mut self, label: Label) -> Self {
        self.fixed.extend(label.rules());
        self
    }

    pub(crate) fn check(&self) -> Result<(), SearchError> {
        if self.memory_budget == Some(0) {
            return Err(SearchError::Settings(
                "the memory budget must be above 0".into(),
            ));
        }
        let mut fixed_sites: Vec<Site> = self.fixed.iter().map(|rule| rule.site).collect();
        fixed_sites.sort_unstable();
        let doubled = fixed_sites.windows(2).find(|pair| pair[0] == pair[1]);
        // Sage applies one fixed terminal modification per end, in no fixed
        // order: a peptide and a protein rule on the same end would make the
        // search differ from run to run.
        let both_ends = |a: Site, b: Site| fixed_sites.contains(&a) && fixed_sites.contains(&b);
        let usable = |tolerance: f32| tolerance.is_finite() && tolerance > 0.0;
        let problem = if !(usable(self.precursor_ppm) && usable(self.fragment.value())) {
            Some("tolerances must be positive, finite numbers".to_string())
        } else if self.min_length == 0 || self.min_length > self.max_length {
            Some("peptide lengths must run from at least 1 up to the maximum".into())
        } else if self.charges.0 == 0 || self.charges.0 > self.charges.1 {
            Some("charges must run from at least 1 up to the maximum".into())
        } else if !self.variable.is_empty() && self.max_variable_mods == 0 {
            Some("variable modifications need max_variable_mods of 1 or more".into())
        } else if both_ends(Site::PeptideN, Site::ProteinN)
            || both_ends(Site::PeptideC, Site::ProteinC)
        {
            Some(
                "fixed modifications on both a peptide terminus (^ $) and the protein terminus \
                 at the same end ([ ]) clash: keep one, or make the protein one variable"
                    .into(),
            )
        } else {
            doubled.map(|pair| {
                format!(
                    "two fixed modifications on {}: a site takes one fixed modification",
                    pair[0].sage()
                )
            })
        };
        problem.map_or(Ok(()), |p| Err(SearchError::Settings(p)))
    }

    /// Sage's database parameters, through its own (serde) configuration so
    /// unset values keep Sage's defaults.
    pub(crate) fn builder(&self) -> Result<Builder, SearchError> {
        let fixed: serde_json::Map<String, serde_json::Value> = self
            .fixed
            .iter()
            .map(|rule| (rule.site.sage(), serde_json::json!(rule.mass)))
            .collect();
        let mut variable: serde_json::Map<String, serde_json::Value> = serde_json::Map::new();
        for rule in &self.variable {
            let masses = variable
                .entry(rule.site.sage())
                .or_insert_with(|| serde_json::json!([]));
            if let Some(list) = masses.as_array_mut() {
                list.push(serde_json::json!(rule.mass));
            }
        }
        let (cleave_at, restrict, c_terminal) = self.enzyme.rule();
        serde_json::from_value(serde_json::json!({
            "enzyme": {
                "missed_cleavages": self.missed_cleavages,
                "min_len": self.min_length,
                "max_len": self.max_length,
                "cleave_at": cleave_at,
                "restrict": restrict,
                "c_terminal": c_terminal,
                "semi_enzymatic": self.semi_specific,
            },
            "fragment_min_mz": FRAGMENT_MZ.0,
            "fragment_max_mz": FRAGMENT_MZ.1,
            "peptide_min_mass": PEPTIDE_MASS.0,
            "peptide_max_mass": PEPTIDE_MASS.1,
            "static_mods": fixed,
            "variable_mods": variable,
            "max_variable_mods": self.max_variable_mods,
            "decoy_tag": DECOY_TAG,
            "generate_decoys": true,
            // Only a label: the FASTA itself is handed over parsed.
            "fasta": "in-memory",
        }))
        .map_err(|error| SearchError::Settings(error.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enzymes_parse_loosely() {
        assert_eq!(Enzyme::parse("Lys-C"), Ok(Enzyme::LysC));
        assert_eq!(Enzyme::parse("lysc"), Ok(Enzyme::LysC));
        assert_eq!(Enzyme::parse("TRYPSIN/P"), Ok(Enzyme::TrypsinP));
        let error = Enzyme::parse("pepsin").unwrap_err();
        assert!(error.contains("chymotrypsin"), "{error}");
    }

    #[test]
    fn modification_rules_parse_by_name_or_mass() {
        let tmt = ModRule::parse("K:TMT6plex").unwrap();
        assert_eq!(tmt.site, Site::Residue('K'));
        assert_eq!(
            (tmt.mass, tmt.accession.as_deref()),
            (229.162_932, Some("UNIMOD:737"))
        );
        // A mass that is exactly a table entry gets its name.
        let cam = ModRule::parse("C:57.021464").unwrap();
        assert_eq!(cam.name.as_deref(), Some("Carbamidomethyl"));
        // A mass that isn't stays a bare mass, never a guessed name.
        let odd = ModRule::parse("^:12.5").unwrap();
        assert_eq!((odd.site, odd.name), (Site::PeptideN, None));
        assert_eq!(ModRule::parse("[:Acetyl").unwrap().site, Site::ProteinN);
        for bad in ["K", "Z:15.99", "KR:15.99", "M:Mystery", "M:0", "M:nan"] {
            assert!(ModRule::parse(bad).is_err(), "{bad}");
        }
        assert_eq!(tmt.to_string(), "K:TMT6plex");
    }

    #[test]
    fn fixed_terminal_rules_at_one_end_clash() {
        let clash = SearchSettings::default()
            .with_label(Label::Tmt)
            .with_fixed(ModRule::parse("[:Acetyl").unwrap());
        let error = clash.check().unwrap_err().to_string();
        assert!(error.contains("clash"), "{error}");
        // As a variable modification it is fine.
        let mut variable = SearchSettings::default().with_label(Label::Tmt);
        variable.variable.push(ModRule::parse("[:Acetyl").unwrap());
        assert!(variable.check().is_ok());
    }

    #[test]
    fn labels_and_conflicts() {
        let tmt = SearchSettings::default().with_label(Label::Tmt);
        let sites: Vec<Site> = tmt.fixed.iter().map(|r| r.site).collect();
        assert_eq!(
            sites,
            [Site::Residue('C'), Site::Residue('K'), Site::PeptideN]
        );
        assert!(tmt.check().is_ok());
        let doubled = tmt.with_label(Label::TmtPro);
        let error = doubled.check().unwrap_err().to_string();
        assert!(error.contains("two fixed modifications on K"), "{error}");
    }

    #[test]
    fn options_change_the_defaults_and_typos_are_errors() {
        let options = SearchOptions {
            enzyme: Some("Lys-C".into()),
            label: Some("tmt10plex".into()),
            variable: vec!["[:Acetyl".into()],
            no_oxidation: true,
            ..SearchOptions::default()
        };
        let settings = SearchSettings::from_options(&options).unwrap();
        assert_eq!(settings.enzyme, Enzyme::LysC);
        assert_eq!(
            settings.describe(),
            "lys-c, up to 2 missed, precursor ±20 ppm, fragment ±20 ppm, \
             fixed C:Carbamidomethyl K:TMT6plex ^:TMT6plex, variable [:Acetyl"
        );
        assert_eq!(
            SearchSettings::from_options(&SearchOptions::default()).unwrap(),
            SearchSettings::default()
        );
        for bad in [
            SearchOptions {
                enzyme: Some("pepsin".into()),
                ..SearchOptions::default()
            },
            SearchOptions {
                label: Some("itraq".into()),
                ..SearchOptions::default()
            },
            SearchOptions {
                fixed: vec!["C:Mystery".into()],
                ..SearchOptions::default()
            },
            SearchOptions {
                fixed: vec!["C:42".into()],
                ..SearchOptions::default()
            },
            SearchOptions {
                precursor_ppm: Some(-1.0),
                ..SearchOptions::default()
            },
            SearchOptions {
                max_memory_gb: Some(0.0),
                ..SearchOptions::default()
            },
            SearchOptions {
                max_memory_gb: Some(f64::NAN),
                ..SearchOptions::default()
            },
        ] {
            assert!(SearchSettings::from_options(&bad).is_err(), "{bad:?}");
        }
        let budgeted = SearchOptions {
            max_memory_gb: Some(2.5),
            ..SearchOptions::default()
        };
        assert_eq!(
            SearchSettings::from_options(&budgeted)
                .unwrap()
                .memory_budget,
            Some(2_500_000_000)
        );
    }

    #[test]
    fn fragment_tolerance_is_ppm_or_da_never_both() {
        let ppm = SearchOptions {
            fragment_ppm: Some(10.0),
            ..SearchOptions::default()
        };
        let ppm = SearchSettings::from_options(&ppm).unwrap();
        assert_eq!(ppm.fragment, FragmentTolerance::Ppm(10.0));
        assert!(
            ppm.describe().contains(", fragment ±10 ppm,"),
            "{}",
            ppm.describe()
        );
        let da = SearchOptions {
            fragment_da: Some(0.5),
            ..SearchOptions::default()
        };
        let da = SearchSettings::from_options(&da).unwrap();
        assert_eq!(da.fragment, FragmentTolerance::Da(0.5));
        assert_eq!(
            da.describe(),
            "trypsin, up to 2 missed, precursor ±20 ppm, fragment ±0.5 Da, \
             fixed C:Carbamidomethyl, variable M:Oxidation"
        );
        // Sage gets the same window either way of the fragment.
        assert!(matches!(da.fragment.sage(), Tolerance::Da(lo, hi) if lo == -0.5 && hi == 0.5));
        assert!(matches!(ppm.fragment.sage(), Tolerance::Ppm(lo, hi) if lo == -10.0 && hi == 10.0));
        let both = SearchOptions {
            fragment_ppm: Some(20.0),
            fragment_da: Some(0.5),
            ..SearchOptions::default()
        };
        let error = SearchSettings::from_options(&both).unwrap_err();
        assert!(error.contains("not both"), "{error}");
        for bad in [0.0, -0.5, f32::NAN, f32::INFINITY] {
            let options = SearchOptions {
                fragment_da: Some(bad),
                ..SearchOptions::default()
            };
            let error = SearchSettings::from_options(&options).unwrap_err();
            assert!(
                error.contains("tolerances must be positive"),
                "{bad}: {error}"
            );
        }
    }

    #[test]
    fn the_builder_carries_enzyme_and_modifications() {
        let settings = SearchSettings {
            enzyme: Enzyme::LysN,
            semi_specific: true,
            variable: vec![
                ModRule::parse("M:Oxidation").unwrap(),
                ModRule::parse("M:15.0").unwrap(),
                ModRule::parse("[:Acetyl").unwrap(),
            ],
            ..SearchSettings::default()
        };
        let parameters = settings.builder().unwrap().make_parameters();
        assert!(!parameters.enzyme.c_terminal.unwrap_or(true));
        assert_eq!(parameters.enzyme.semi_enzymatic, Some(true));
        assert_eq!(parameters.static_mods.len(), 1);
        assert_eq!(parameters.variable_mods.len(), 2, "M (two masses) and [");
    }
}
