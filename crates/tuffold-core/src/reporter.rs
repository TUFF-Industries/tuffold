//! Isobaric labels (TMT, TMTpro, iTRAQ) and their reporter ions.
//!
//! Every TMT and TMTpro reporter ion is C₈H₁₆N⁺; channels differ only in
//! which atoms are heavy. After 126, channel *i* carries ⌊*i*/2⌋ ¹³C and, when
//! *i* is odd, one ¹⁵N (its name ends in `N`; ¹³C-only channels end in `C`),
//! so 127N, 127C, 128N, … 135N. The m/z values are computed from that
//! composition rather than copied from tables.
//!
//! TMT6 has two reagent sets with the same channel names: today's (127N and
//! 129N, a subset of TMT10) and the original one sold until about 2012,
//! whose 127 and 129 carried ¹³C instead (127C, 129C), 6.3 mDa heavier.
//!
//! Every iTRAQ reporter is C₆H₁₃N₂⁺ (N-methylpiperazine), 113 at its lightest,
//! with ¹³C and ¹⁵N making up the rest: nominal masses 114–117 (4-plex) and
//! 113–119 and 121 (8-plex; 120 would sit on phenylalanine's immonium ion).
//! Published tables list 0.55 mDa more, the electron left in; MaxQuant's
//! values match the composition. Channels are a whole Da apart.

mod impurities;

pub use impurities::{IMPURITY_COLUMNS, Impurities};

/// ¹²C₈ ¹H₁₆ ¹⁴N, less an electron: the 126 reporter ion.
const BASE_MZ: f64 = 8.0 * 12.0 + 16.0 * 1.007_825_032_07 + 14.003_074_004_8 - 0.000_548_579_909;
const C13_SHIFT: f64 = 1.003_354_837_8;
const N15_SHIFT: f64 = 0.997_034_893_4;
/// ¹²C₆ ¹H₁₃ ¹⁴N₂, less an electron: the 113 iTRAQ reporter ion.
const ITRAQ_BASE_MZ: f64 =
    6.0 * 12.0 + 13.0 * 1.007_825_032_07 + 2.0 * 14.003_074_004_8 - 0.000_548_579_909;
/// (nominal mass, ¹³C, ¹⁵N) of the iTRAQ reporters.
const ITRAQ: [(usize, usize, usize); 8] = [
    (113, 0, 0),
    (114, 1, 0),
    (115, 1, 1),
    (116, 2, 1),
    (117, 3, 1),
    (118, 3, 2),
    (119, 4, 2),
    (121, 6, 2),
];

/// A labelling reagent set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Plex {
    Tmt6,
    /// The original TMT6 reagents (127 and 129 are today's 127C and 129C).
    Tmt6Original,
    Tmt10,
    Tmt11,
    TmtPro16,
    TmtPro18,
    Itraq4,
    Itraq8,
}

/// One reporter channel.
#[derive(Debug, Clone, PartialEq)]
pub struct Channel {
    pub name: String,
    pub mz: f64,
    /// TMT: position in the 126, 127N, 127C, … series, the reporter carrying
    /// ⌊index/2⌋ ¹³C and index mod 2 ¹⁵N more than 126. iTRAQ: Da above 113.
    pub index: usize,
}

impl Plex {
    pub const ALL: [Plex; 8] = [
        Plex::Tmt6,
        Plex::Tmt6Original,
        Plex::Tmt10,
        Plex::Tmt11,
        Plex::TmtPro16,
        Plex::TmtPro18,
        Plex::Itraq4,
        Plex::Itraq8,
    ];

    pub fn parse(text: &str) -> Result<Self, String> {
        match text
            .to_ascii_lowercase()
            .replace(['-', '_', ' '], "")
            .as_str()
        {
            "tmt6" | "tmt6plex" => Ok(Self::Tmt6),
            "tmt6original" | "tmt6plexoriginal" | "tmt6old" => Ok(Self::Tmt6Original),
            "tmt10" | "tmt10plex" => Ok(Self::Tmt10),
            "tmt11" | "tmt11plex" => Ok(Self::Tmt11),
            "tmtpro16" | "tmtpro16plex" | "tmt16" | "tmt16plex" => Ok(Self::TmtPro16),
            "tmtpro18" | "tmtpro18plex" | "tmt18" | "tmt18plex" => Ok(Self::TmtPro18),
            "itraq4" | "itraq4plex" => Ok(Self::Itraq4),
            "itraq8" | "itraq8plex" => Ok(Self::Itraq8),
            _ => Err(format!(
                "unknown plex '{text}' (known: tmt6, tmt6-original, tmt10, tmt11, tmtpro16, \
                 tmtpro18, itraq4, itraq8)"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Tmt6 => "tmt6",
            Self::Tmt6Original => "tmt6-original",
            Self::Tmt10 => "tmt10",
            Self::Tmt11 => "tmt11",
            Self::TmtPro16 => "tmtpro16",
            Self::TmtPro18 => "tmtpro18",
            Self::Itraq4 => "itraq4",
            Self::Itraq8 => "itraq8",
        }
    }

    /// The Unimod name of the label the peptides carry.
    pub fn label(self) -> &'static str {
        match self {
            Self::Tmt6 | Self::Tmt6Original | Self::Tmt10 | Self::Tmt11 => "TMT6plex",
            Self::TmtPro16 | Self::TmtPro18 => "TMTpro",
            Self::Itraq4 => "iTRAQ4plex",
            Self::Itraq8 => "iTRAQ8plex",
        }
    }

    /// Whether this is an iTRAQ reagent set, whose certificates give nominal
    /// shifts on channels a whole Da apart.
    pub fn is_itraq(self) -> bool {
        matches!(self, Self::Itraq4 | Self::Itraq8)
    }

    /// How far from a channel's m/z a reporter is read. TMT's must split the
    /// N and C forms of a nominal mass, 6.3 mDa (50 ppm) apart. iTRAQ's have
    /// no such neighbour, and its certificates count a heavier or lighter
    /// reporter against the channel of that nominal mass, up to 55 ppm off
    /// it, so the window takes those in (and sums them): what the correction
    /// assumes.
    pub fn tolerance_ppm(self) -> f64 {
        if self.is_itraq() { 60.0 } else { 20.0 }
    }

    /// Channels in reagent order.
    pub fn channels(self) -> Vec<Channel> {
        match self {
            // The original six: 126, 127N, 128C, 129N, 130C, 131N.
            Self::Tmt6 => sixplex([0, 1, 4, 5, 8, 9]),
            // 126, 127C, 128C, 129C, 130C, 131N.
            Self::Tmt6Original => sixplex([0, 2, 4, 6, 8, 9]),
            Self::Tmt10 => {
                let mut channels = series(10);
                channels[9].name = "131".into();
                channels
            }
            Self::Tmt11 => series(11),
            Self::TmtPro16 => series(16),
            Self::TmtPro18 => series(18),
            Self::Itraq4 => itraq(&ITRAQ[1..5]),
            Self::Itraq8 => itraq(&ITRAQ),
        }
    }

    /// Reagent sets a plex of this one could have been labelled with instead:
    /// the same label and channel names at other m/z (the two TMT6 sets).
    pub fn alternatives(self) -> Vec<Plex> {
        let names =
            |plex: Plex| -> Vec<String> { plex.channels().into_iter().map(|c| c.name).collect() };
        Self::ALL
            .into_iter()
            .filter(|&other| {
                other != self && other.label() == self.label() && names(other) == names(self)
            })
            .collect()
    }
}

impl std::fmt::Display for Plex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

fn channel_mz(i: usize) -> f64 {
    BASE_MZ + (i / 2) as f64 * C13_SHIFT + (i % 2) as f64 * N15_SHIFT
}

fn sixplex(indices: [usize; 6]) -> Vec<Channel> {
    indices
        .into_iter()
        .map(|i| Channel {
            name: (126 + i / 2 + i % 2).to_string(),
            mz: channel_mz(i),
            index: i,
        })
        .collect()
}

fn itraq(reporters: &[(usize, usize, usize)]) -> Vec<Channel> {
    reporters
        .iter()
        .map(|&(nominal, c13, n15)| Channel {
            name: nominal.to_string(),
            mz: ITRAQ_BASE_MZ + c13 as f64 * C13_SHIFT + n15 as f64 * N15_SHIFT,
            index: nominal - 113,
        })
        .collect()
}

fn series(count: usize) -> Vec<Channel> {
    (0..count)
        .map(|i| Channel {
            name: if i == 0 {
                "126".into()
            } else {
                format!(
                    "{}{}",
                    126 + i / 2 + i % 2,
                    if i % 2 == 1 { 'N' } else { 'C' }
                )
            },
            mz: channel_mz(i),
            index: i,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn find(plex: Plex, name: &str) -> f64 {
        plex.channels()
            .into_iter()
            .find(|c| c.name == name)
            .unwrap_or_else(|| panic!("{plex} has no {name}"))
            .mz
    }

    #[test]
    fn reporter_masses_match_the_published_ones() {
        // 126–128N: Tephra's TMTpro table (6 decimals); 134C, 135N: Thermo's
        // TMTpro brochure (5 decimals); TMT6: Thermo's TMT6plex values.
        for (name, mz) in [
            ("126", 126.127726),
            ("127N", 127.124761),
            ("127C", 127.131081),
            ("128N", 128.128116),
        ] {
            assert!((find(Plex::TmtPro18, name) - mz).abs() < 1e-6, "{name}");
        }
        for (name, mz) in [("134C", 134.15457), ("135N", 135.15160)] {
            assert!((find(Plex::TmtPro18, name) - mz).abs() < 1e-5, "{name}");
        }
        for (name, mz) in [
            ("126", 126.127726),
            ("127", 127.124761),
            ("128", 128.134436),
            ("129", 129.131471),
            ("130", 130.141145),
            ("131", 131.138180),
        ] {
            assert!((find(Plex::Tmt6, name) - mz).abs() < 1e-6, "{name}");
        }
        // The original TMT6's 127 and 129 are Thermo's 127C and 129C.
        for (name, mz) in [("127", 127.131081), ("129", 129.137790)] {
            assert!((find(Plex::Tmt6Original, name) - mz).abs() < 1e-6, "{name}");
        }
        assert_eq!(Plex::Tmt6.alternatives(), [Plex::Tmt6Original]);
        assert_eq!(Plex::Tmt6Original.alternatives(), [Plex::Tmt6]);
        assert!(Plex::Tmt10.alternatives().is_empty());
    }

    #[test]
    fn itraq_masses_match_the_published_ones_less_the_electron() {
        // MaxQuant's 4-plex values; the 8-plex's others are the usual table's
        // (113.1078, 118.1120, 119.1153, 121.1220) less 0.55 mDa.
        for (name, mz) in [
            ("114", 114.11068),
            ("115", 115.107715),
            ("116", 116.111069),
            ("117", 117.114424),
        ] {
            assert!((find(Plex::Itraq4, name) - mz).abs() < 1e-5, "{name}");
            assert!((find(Plex::Itraq8, name) - mz).abs() < 1e-5, "{name}");
        }
        for (name, mz) in [
            ("113", 113.1078),
            ("118", 118.1120),
            ("119", 119.1153),
            ("121", 121.1220),
        ] {
            let electron = 0.000_548_579_909;
            assert!(
                (find(Plex::Itraq8, name) - (mz - electron)).abs() < 1e-4,
                "{name}"
            );
        }
        assert_eq!(Plex::Itraq4.channels().len(), 4);
        assert_eq!(Plex::Itraq8.channels().len(), 8);
        assert_eq!(Plex::parse("iTRAQ-8plex"), Ok(Plex::Itraq8));
        assert!(Plex::Itraq4.alternatives().is_empty());
        for label in [Plex::Itraq4.label(), Plex::Itraq8.label()] {
            assert!(crate::unimod::by_name(label).is_some(), "{label}");
        }
    }

    #[test]
    fn plexes_have_their_channels_in_order() {
        let names =
            |plex: Plex| -> Vec<String> { plex.channels().into_iter().map(|c| c.name).collect() };
        assert_eq!(names(Plex::Tmt6).len(), 6);
        assert_eq!(names(Plex::Tmt10).last().unwrap(), "131");
        assert_eq!(names(Plex::Tmt11).last().unwrap(), "131C");
        assert_eq!(names(Plex::TmtPro16).last().unwrap(), "134N");
        assert_eq!(names(Plex::TmtPro18).len(), 18);
        assert_eq!(names(Plex::TmtPro18)[1..3], ["127N", "127C"]);
        for plex in Plex::ALL {
            let mz: Vec<f64> = plex.channels().iter().map(|c| c.mz).collect();
            assert!(mz.windows(2).all(|w| w[0] < w[1]), "{plex}");
            assert_eq!(Plex::parse(plex.name()), Ok(plex));
        }
    }
}
