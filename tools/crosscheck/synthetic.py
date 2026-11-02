"""Run mokapot on Tuffold's synthetic validation data.

Ports the generator in crates/tuffold-validate/tests/support/mod.rs draw for
draw (SplitMix64, Box-Muller), hands engine 0 to mokapot as a PSM dataset, and
reports how many target PSMs mokapot accepts at 1% FDR. The Rust test
`one_engine_matches_mokapot_on_the_same_data` records that number.

Usage (needs mokapot; Python's standard library otherwise):
    python tools/crosscheck/synthetic.py [--spectra 6000] [--seed 20261008]
"""

import argparse
import math

MASK = (1 << 64) - 1
PROTEINS = 2000
PRESENT = 200
PEPTIDES_PER_PROTEIN = 5


class SplitMix:
    def __init__(self, seed):
        self.state = seed & MASK

    def next(self):
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def uniform(self):
        return (self.next() >> 11) * (1.0 / (1 << 53))

    def normal(self):
        u1 = self.uniform()
        u2 = self.uniform()
        return math.sqrt(-2.0 * math.log(1.0 - u1)) * math.cos(2.0 * math.pi * u2)

    def index(self, n):
        return self.next() % n


def generate(seed, spectra, engines):
    rng = SplitMix(seed)
    rows = []
    for spectrum in range(1, spectra + 1):
        correct = rng.uniform() < 0.5
        true_peptide = (rng.index(PRESENT), rng.index(PEPTIDES_PER_PROTEIN), False)
        for engine in range(engines):
            if correct:
                peptide = true_peptide
            else:
                decoy = rng.uniform() < 0.5
                peptide = (rng.index(PROTEINS), rng.index(PEPTIDES_PER_PROTEIN), decoy)
            shift = 1.0 if correct else 0.0
            f1 = rng.normal() + 1.6 * shift
            f2 = rng.normal() + 1.2 * shift
            expect = 10.0 ** -(rng.normal() + 1.0 * shift + 2.0)
            rows.append((spectrum, engine, peptide, peptide == true_peptide, f1, f2, expect))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--spectra", type=int, default=6000)
    parser.add_argument("--runs", type=int, default=10)
    args = parser.parse_args()

    rows = [r for r in generate(args.seed, args.spectra, 1) if r[1] == 0]
    print(f"checksum: f1 sum {sum(r[4] for r in rows):.9f}, correct {sum(r[3] for r in rows)}")

    import numpy as np
    import pandas as pd
    import mokapot
    import mokapot.qvalues

    # Built directly rather than through mokapot.read_pin, whose parser
    # (mokapot 0.10) does not work with pandas 3.
    frame = pd.DataFrame(
        {
            "SpecId": [r[0] for r in rows],
            "Label": [not r[2][2] for r in rows],
            "ScanNr": [r[0] for r in rows],
            "f1": [r[4] for r in rows],
            "f2": [r[5] for r in rows],
            "neglogexpect": [-math.log10(max(r[6], 1e-300)) for r in rows],
            "Peptide": [f"PSM{r[0]}" for r in rows],
        }
    )
    # The best possible score knows the generator: each feature's shift for
    # correct matches, in standard deviations (1.6, 1.2, 1.0).
    oracle = 1.6 * frame["f1"] + 1.2 * frame["f2"] + 1.0 * frame["neglogexpect"]
    q = mokapot.qvalues.tdc(oracle.to_numpy(), frame["Label"].to_numpy())
    print(f"oracle score: {int(((q <= 0.01) & frame['Label']).sum())} target PSMs at 1% FDR")

    # mokapot's SVM training is randomized (its own generator, not numpy's
    # global seed) and varies a lot between runs, so report several seeds.
    counts = []
    for seed in range(1, args.runs + 1):
        psms = mokapot.LinearPsmDataset(
            psms=frame,
            target_column="Label",
            spectrum_columns=["ScanNr"],
            peptide_column="Peptide",
            feature_columns=["f1", "f2", "neglogexpect"],
            rng=seed,
        )
        results, _ = mokapot.brew(psms, rng=seed)
        table = results.confidence_estimates["psms"]
        counts.append(int(((table["mokapot q-value"] <= 0.01) & table["Label"]).sum()))
    counts.sort()
    print(f"mokapot {mokapot.__version__}, target PSMs at 1% FDR over seeds 1..{args.runs}:")
    print(f"  {counts}; median {counts[len(counts) // 2]}")


if __name__ == "__main__":
    main()
