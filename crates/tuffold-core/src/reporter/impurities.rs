//! Reagent isotope impurities, from the certificate of analysis that comes
//! with each TMT or iTRAQ lot.
//!
//! Every reporter is C₈H₁₆N⁺ with some of its carbons ¹³C and maybe its
//! nitrogen ¹⁵N ([`Channel::index`]). A reagent is never perfectly enriched,
//! so a few percent of a channel's reporter ions carry one or two ¹³C or ¹⁵N
//! more or less than they should, and appear at another m/z: on another
//! channel when one has that composition, else nowhere Tuffold reads. The
//! certificate gives those percentages per channel, of all the channel's
//! reporter ions (what OpenMS assumes too), so the channel keeps 100 minus
//! their sum.
//!
//! Where a TMT impurity lands follows from composition alone: Tuffold reads
//! those reporters within 20 ppm, which separates the N and C forms of a nominal
//! mass (6.3 mDa apart), so an impurity counts against another channel only
//! when it has exactly that channel's composition. Certificates come in two
//! layouts:
//! - eight columns (TMTpro): −2 ¹³C, −1 ¹⁵N −1 ¹³C, −1 ¹³C, −1 ¹⁵N, +1 ¹⁵N,
//!   +1 ¹³C, +1 ¹⁵N +1 ¹³C, +2 ¹³C ([`IMPURITY_COLUMNS`]), which name the
//!   isotope; a nonzero value the channel's composition rules out (127C
//!   losing a ¹⁵N it doesn't have) is refused, since it means the columns
//!   were misread;
//! - four columns (TMT6, TMT10, TMT11): −2, −1, +1, +2 Da. A shift is read
//!   as ¹³C where the reporter can take it, else as ¹⁵N (127N one Da lighter
//!   has lost its ¹⁵N and sits on 126), else refused.
//!
//! iTRAQ certificates have the four columns, and a shift lands on the channel
//! of that nominal mass whichever isotope it is, as MaxQuant and OpenMS
//! treat them; the reporters are read wide enough to include them
//! ([`Plex::tolerance_ppm`]), and every peak in the window is summed. It's kept in the
//! ¹³C columns, and written back as four columns.

use super::Plex;

/// The shifts of an eight-column certificate, in its column order: (¹³C,
/// ¹⁵N) added to the reporter.
pub const IMPURITY_COLUMNS: [(i32, i32); 8] = [
    (-2, 0),
    (-1, -1),
    (-1, 0),
    (0, -1),
    (0, 1),
    (1, 0),
    (1, 1),
    (2, 0),
];

/// A reporter has eight carbons and one nitrogen.
const CARBONS: i32 = 8;

/// One lot's impurities for every channel of its reagent set.
#[derive(Debug, Clone, PartialEq)]
pub struct Impurities {
    plex: Plex,
    /// Per channel of the reagent set, in channel order: percent of its
    /// reporter ions at each [`IMPURITY_COLUMNS`] shift.
    percent: Vec<[f64; 8]>,
}

impl Impurities {
    /// Read a certificate: one line per channel of `plex`, the channel's name
    /// and then four or eight percentages (`NA` or `-` for none), separated
    /// by spaces, tabs, commas, semicolons, or slashes. `#` starts a comment;
    /// a line starting with `channel` is a header. Every channel of the
    /// reagent set must be listed (all zeros for a pure one), so nothing is
    /// assumed.
    pub fn parse(text: &str, plex: Plex) -> Result<Self, String> {
        let channels = plex.channels();
        let mut percent: Vec<Option<[f64; 8]>> = vec![None; channels.len()];
        for (number, raw) in text.lines().enumerate() {
            let line = number + 1;
            let content = raw.split('#').next().unwrap_or_default();
            let tokens: Vec<&str> = content
                .split(|c: char| c.is_whitespace() || "/,;".contains(c))
                .filter(|t| !t.is_empty())
                .collect();
            let Some((&first, values)) = tokens.split_first() else {
                continue;
            };
            let first = first.trim_end_matches(':');
            if first.eq_ignore_ascii_case("channel") {
                continue;
            }
            let Some(at) = channels
                .iter()
                .position(|c| c.name.eq_ignore_ascii_case(first))
            else {
                return Err(format!(
                    "line {line}: {plex} has no channel '{first}' (channels: {})",
                    names(plex)
                ));
            };
            let name = &channels[at].name;
            if percent[at].is_some() {
                return Err(format!("line {line}: channel {name} is listed twice"));
            }
            let values = values
                .iter()
                .map(|v| value(v).map_err(|e| format!("line {line}, channel {name}: {e}")))
                .collect::<Result<Vec<f64>, String>>()?;
            let index = channels[at].index;
            let row = match values.len() {
                4 if plex.is_itraq() => Ok(nominal_columns(&values)),
                8 if plex.is_itraq() => {
                    Err("8 values; an iTRAQ certificate has 4 (−2, −1, +1, +2 Da)".to_string())
                }
                8 => eight_columns(&values, index),
                4 => four_columns(&values, index),
                n => Err(format!(
                    "{n} values; give 4 (−2, −1, +1, +2 Da) or 8 (−2C, −N−C, −C, −N, +N, +C, \
                     +N+C, +2C)"
                )),
            }
            .map_err(|e| format!("line {line}, channel {name}: {e}"))?;
            let total: f64 = row.iter().sum();
            if total >= 100.0 {
                return Err(format!(
                    "line {line}, channel {name}: impurities add up to {total}%, leaving \
                     nothing for the channel itself"
                ));
            }
            percent[at] = Some(row);
        }
        let missing: Vec<&str> = channels
            .iter()
            .zip(&percent)
            .filter(|(_, row)| row.is_none())
            .map(|(c, _)| c.name.as_str())
            .collect();
        if !missing.is_empty() {
            return Err(format!(
                "no impurities for {}: list every channel of the {plex} kit (all zeros for a \
                 pure one)",
                missing.join(", ")
            ));
        }
        Ok(Self {
            plex,
            percent: percent.into_iter().flatten().collect(),
        })
    }

    pub fn plex(&self) -> Plex {
        self.plex
    }

    /// Each channel's name and its eight-column percentages.
    pub fn rows(&self) -> impl Iterator<Item = (String, [f64; 8])> + '_ {
        self.plex
            .channels()
            .into_iter()
            .zip(&self.percent)
            .map(|(c, row)| (c.name, *row))
    }

    /// The certificate (eight columns; iTRAQ's four), which
    /// [`Impurities::parse`] reads back.
    pub fn to_text(&self) -> String {
        let itraq = self.plex.is_itraq();
        self.rows()
            .map(|(name, row)| {
                let values: Vec<String> = row
                    .iter()
                    .zip(IMPURITY_COLUMNS)
                    .filter(|(_, (_, n15))| !itraq || *n15 == 0)
                    .map(|(p, _)| p.to_string())
                    .collect();
                format!("{name}\t{}\n", values.join("/"))
            })
            .collect()
    }

    /// The share of each channel's reporter ions read at each channel, in
    /// channel order: `matrix[observed][true]`. Columns sum to 1 less what
    /// lands on no channel of the set.
    pub fn matrix(&self) -> Vec<Vec<f64>> {
        let channels = self.plex.channels();
        let n = channels.len();
        let mut matrix = vec![vec![0.0; n]; n];
        for (j, (channel, row)) in channels.iter().zip(&self.percent).enumerate() {
            matrix[j][j] = 1.0 - row.iter().sum::<f64>() / 100.0;
            for (&(c13, n15), &p) in IMPURITY_COLUMNS.iter().zip(row) {
                let target = if self.plex.is_itraq() {
                    usize::try_from(channel.index as i32 + c13 + n15).ok()
                } else {
                    shifted(channel.index, c13, n15)
                };
                let Some(target) = target else {
                    continue;
                };
                if let Some(i) = channels.iter().position(|c| c.index == target) {
                    matrix[i][j] += p / 100.0;
                }
            }
        }
        matrix
    }
}

fn names(plex: Plex) -> String {
    plex.channels()
        .into_iter()
        .map(|c| c.name)
        .collect::<Vec<_>>()
        .join(", ")
}

fn value(text: &str) -> Result<f64, String> {
    if text.eq_ignore_ascii_case("na") || text == "-" || text == "—" {
        return Ok(0.0);
    }
    let v: f64 = text
        .trim_end_matches('%')
        .parse()
        .map_err(|_| format!("'{text}' is not a percentage"))?;
    if !(0.0..=100.0).contains(&v) {
        return Err(format!("{v} is not a percentage between 0 and 100"));
    }
    Ok(v)
}

/// The series index of a reporter `index` with `c13` ¹³C and `n15` ¹⁵N
/// more, if a reporter can have that composition.
fn shifted(index: usize, c13: i32, n15: i32) -> Option<usize> {
    let carbons = i32::try_from(index / 2).ok()? + c13;
    let nitrogen = i32::try_from(index % 2).ok()? + n15;
    ((0..=CARBONS).contains(&carbons) && (0..=1).contains(&nitrogen))
        .then(|| usize::try_from(2 * carbons + nitrogen).ok())
        .flatten()
}

fn describe((c13, n15): (i32, i32)) -> String {
    let part = |n: i32, isotope: &str| match n {
        0 => None,
        n => Some(format!("{n:+} {isotope}")),
    };
    [part(n15, "¹⁵N"), part(c13, "¹³C")]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .join(" ")
}

fn eight_columns(values: &[f64], index: usize) -> Result<[f64; 8], String> {
    let mut row = [0.0; 8];
    for ((&shift, &p), slot) in IMPURITY_COLUMNS.iter().zip(values).zip(&mut row) {
        if p > 0.0 && shifted(index, shift.0, shift.1).is_none() {
            return Err(format!(
                "{p}% at {}, which this reporter's composition rules out (that column should \
                 be NA or 0; are the columns in the order −2C, −N−C, −C, −N, +N, +C, +N+C, \
                 +2C?)",
                describe(shift)
            ));
        }
        *slot = p;
    }
    Ok(row)
}

/// −2, −1, +1, +2 Da into the ¹³C columns.
fn nominal_columns(values: &[f64]) -> [f64; 8] {
    let mut row = [0.0; 8];
    for (&p, da) in values.iter().zip([-2, -1, 1, 2]) {
        if let Some(column) = IMPURITY_COLUMNS.iter().position(|&c| c == (da, 0)) {
            row[column] = p;
        }
    }
    row
}

fn four_columns(values: &[f64], index: usize) -> Result<[f64; 8], String> {
    let mut row = [0.0; 8];
    for (&p, da) in values.iter().zip([-2, -1, 1, 2]) {
        if p == 0.0 {
            continue;
        }
        // ¹³C where possible, else one of the Da as ¹⁵N.
        let shift = [(da, 0), (da + 1, -1), (da - 1, 1)]
            .into_iter()
            .find(|&(c13, n15)| shifted(index, c13, n15).is_some())
            .ok_or_else(|| {
                format!("{p}% at {da:+} Da, which no reporter of this channel can be")
            })?;
        let column = IMPURITY_COLUMNS
            .iter()
            .position(|&c| c == shift)
            .ok_or_else(|| format!("{p}% at {da:+} Da has no column for {}", describe(shift)))?;
        row[column] += p;
    }
    Ok(row)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OpenMS's default TMT10 certificate (its `correction_matrix`), four
    /// columns per channel.
    const TMT10: &str = "\
channel  -2  -1  +1  +2
126   0.0/0.0/5.09/0.0
127N  0.0/0.25/5.27/0.0
127C  0.0/0.37/5.36/0.15
128N  0.0/0.65/4.17/0.1
128C  0.08/0.49/3.06/0.0
129N  0.01/0.71/3.07/0.0
129C  0.0/1.32/2.62/0.0
130N  0.02/1.28/2.75/2.53
130C  0.03/2.08/2.23/0.0
131   0.08/1.99/1.65/0.0  # last channel
";

    fn column(matrix: &[Vec<f64>], j: usize) -> Vec<f64> {
        matrix.iter().map(|row| row[j]).collect()
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-12)
    }

    #[test]
    fn four_columns_land_by_composition() {
        let matrix = Impurities::parse(TMT10, Plex::Tmt10).unwrap().matrix();
        // Channels: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131.
        // 127C (one ¹³C): −1 → 126, +1 → 128C, +2 → 129C; keeps
        // 100 − 0.37 − 5.36 − 0.15 = 94.12%.
        assert!(close(
            &column(&matrix, 2),
            &[0.0037, 0.0, 0.9412, 0.0, 0.0536, 0.0, 0.0015, 0.0, 0.0, 0.0]
        ));
        // 127N (no ¹³C): one Da lighter it has lost its ¹⁵N, i.e. 126; one
        // heavier it gained a ¹³C, 128N.
        assert!(close(
            &column(&matrix, 1),
            &[0.0025, 0.9448, 0.0, 0.0527, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ));
        // 131 (four ¹³C, one ¹⁵N): −2 → 129N, −1 → 130N, +1 lands on 132N,
        // which TMT10 lacks, so it only lowers 131's own share.
        assert!(close(
            &column(&matrix, 9),
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0008, 0.0, 0.0199, 0.0, 0.9628]
        ));
    }

    #[test]
    fn itraq_shifts_land_by_nominal_mass_and_round_trip() {
        let text = "113 0/0/6/0.2\n114 0/1/5/0.1\n115 0/2/4/0\n116 0/3/3/0\n\
                    117 0.1/4/2/0\n118 0/5/1/0\n119 0/6/1/0\n121 0.2/7/0.5/0\n";
        let impurities = Impurities::parse(text, Plex::Itraq8).unwrap();
        let m = impurities.matrix();
        // 113's +1 Da reads at 114 and its +2 at 115; it keeps the rest.
        assert!(close(&column(&m, 0)[..3], &[0.938, 0.06, 0.002]));
        // 119's +1 is at 120, no channel; 121's −1 too, and its −2 is at 119.
        assert!((column(&m, 6).iter().sum::<f64>() - 0.99).abs() < 1e-9);
        assert!((m[6][7] - 0.002).abs() < 1e-9 && m[7][6] == 0.0);
        assert!((column(&m, 7).iter().sum::<f64>() - 0.925).abs() < 1e-9);
        let again = Impurities::parse(&impurities.to_text(), Plex::Itraq8).unwrap();
        assert_eq!(again, impurities);
        let four = "114 0/1/5.9/0.2\n115 0/2/5.6/0.1\n116 0/3/4.5/0.1\n117 0.1/4/3.5/0\n";
        let m = Impurities::parse(four, Plex::Itraq4).unwrap().matrix();
        // 114's −1 is at 113, not in the 4-plex.
        assert!((column(&m, 0).iter().sum::<f64>() - 0.99).abs() < 1e-9);
        assert!((m[1][0] - 0.059).abs() < 1e-9);
        let eight = "114 0/0/0/0/0/0/0/0\n115 0/0/0/0/0/0/0/0\n\
                     116 0/0/0/0/0/0/0/0\n117 0/0/0/0/0/0/0/0\n";
        assert!(Impurities::parse(eight, Plex::Itraq4).is_err());
    }

    #[test]
    fn eight_columns_name_the_isotope_and_refuse_what_cannot_be() {
        // OpenMS's TMTpro 16 defaults for the first three channels; the rest
        // are pure.
        let mut text = String::from(
            "126  NA/NA / NA/NA / 0.31/9.09 / 0.02/0.32\n\
             127N NA/NA / NA/0.78 / NA/9.41 / NA/0.33\n\
             127C NA/NA / 0.93/NA / 0.35/8.63 / 0.01/0.27\n",
        );
        for c in &Plex::TmtPro16.channels()[3..] {
            text.push_str(&format!("{} 0 0 0 0 0 0 0 0\n", c.name));
        }
        let impurities = Impurities::parse(&text, Plex::TmtPro16).unwrap();
        let matrix = impurities.matrix();
        // 126: +¹⁵N → 127N, +¹³C → 127C, +¹⁵N+¹³C → 128N, +2¹³C → 128C.
        let mut expected = vec![0.0; 16];
        expected[..5].copy_from_slice(&[0.9026, 0.0031, 0.0909, 0.0002, 0.0032]);
        assert!(close(&column(&matrix, 0), &expected));
        // Read back from its own text.
        assert_eq!(
            Impurities::parse(&impurities.to_text(), Plex::TmtPro16).unwrap(),
            impurities
        );

        // 127C has no ¹⁵N to lose.
        let wrong = text.replace("127C NA/NA / 0.93/NA", "127C NA/NA / 0.93/0.5");
        let error = Impurities::parse(&wrong, Plex::TmtPro16).unwrap_err();
        assert!(
            error.starts_with("line 3, channel 127C: 0.5% at -1 ¹⁵N"),
            "{error}"
        );
    }

    #[test]
    fn certificates_must_be_whole_and_well_formed() {
        let without_131 = TMT10
            .lines()
            .filter(|l| !l.starts_with("131"))
            .collect::<Vec<_>>();
        assert_eq!(
            Impurities::parse(&without_131.join("\n"), Plex::Tmt10).unwrap_err(),
            "no impurities for 131: list every channel of the tmt10 kit (all zeros for a pure \
             one)"
        );
        let twice = format!("{TMT10}126 0 0 0 0\n");
        assert_eq!(
            Impurities::parse(&twice, Plex::Tmt10).unwrap_err(),
            "line 12: channel 126 is listed twice"
        );
        let three = TMT10.replace("126   0.0/0.0/5.09/0.0", "126 0 0 5");
        assert!(
            Impurities::parse(&three, Plex::Tmt10)
                .unwrap_err()
                .starts_with("line 2, channel 126: 3 values")
        );
        // 126 has no heavy atom to lose.
        let lighter = TMT10.replace("126   0.0/0.0/5.09/0.0", "126 0 1 5 0");
        assert_eq!(
            Impurities::parse(&lighter, Plex::Tmt10).unwrap_err(),
            "line 2, channel 126: 1% at -1 Da, which no reporter of this channel can be"
        );
        assert!(
            Impurities::parse(&TMT10.replace("5.09", "101"), Plex::Tmt10)
                .unwrap_err()
                .contains("not a percentage between 0 and 100")
        );
        assert!(
            Impurities::parse(TMT10, Plex::Tmt11)
                .unwrap_err()
                .contains("no channel '131'")
        );
    }
}
