"""Are the extra spectra a validation change accepts right? Checked by
retention time, which no Tuffold feature fixes for spectra of a peptide
already accepted in the run (read-only).

BEFORE and AFTER are the same data validated two ways (copies of one
project). For spectra AFTER accepts at 1% and BEFORE doesn't, on peptides
BEFORE accepts in that run, it prints the share within 1 min of those
peptides' median RT in BEFORE, against three references:
- consensus: BEFORE's own spectra of those peptides, leave-one-out;
- marginal: BEFORE's spectra of the same peptides with q > 0.5%, the ones
  BEFORE barely accepts (weighted per peptide like the extras): the matched
  control, since extras are marginal too;
- shuffled: the extras against a random accepted peptide's median.
Extras at least as consistent as the marginal ones show no sign of wrong
spectra on present peptides, the error target-decoy counts can miss when
spectra are pooled. Abundant peptides smeared over
the gradient lower every share, so compare shares, not 100%.

Usage (Python's standard library only):
    python3 -I tools/crosscheck/rt_consistency.py BEFORE.tuff AFTER.tuff [--window 60]
"""

import argparse
import random
import sqlite3
import statistics
from collections import defaultdict


def accepted(path):
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    return {
        spectrum: (peptide, run, rt, q)
        for spectrum, peptide, run, rt, q in db.execute(
            """SELECT r.spectrum_id, r.peptide_id, s.ms_run_id, s.rt_seconds, r.q_value
               FROM spectrum_results r JOIN spectra s ON s.id = r.spectrum_id
               WHERE r.q_value <= 0.01 AND r.is_decoy = 0 AND s.rt_seconds IS NOT NULL"""
        )
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("before")
    parser.add_argument("after")
    parser.add_argument("--window", type=float, default=60, help="seconds")
    args = parser.parse_args()
    before, after = accepted(args.before), accepted(args.after)
    consensus = defaultdict(list)
    for spectrum, (peptide, run, rt, q) in before.items():
        consensus[(peptide, run)].append((spectrum, rt, q))
    extras = defaultdict(list)
    new_peptides = 0
    for spectrum, (peptide, run, rt, _) in after.items():
        if spectrum in before:
            continue
        if (peptide, run) in consensus:
            extras[(peptide, run)].append(rt)
        else:
            new_peptides += 1
    near = lambda rt, ref: abs(rt - ref) <= args.window
    median = {key: statistics.median(rt for _, rt, _ in items) for key, items in consensus.items()}

    def leave_one_out(items, spectrum, rt):
        return near(rt, statistics.median(r for s, r, _ in items if s != spectrum))

    n_extra = sum(len(rts) for rts in extras.values())
    extra_near = sum(near(rt, median[key]) for key, rts in extras.items() for rt in rts)
    loo = [leave_one_out(items, s, rt) for items in consensus.values() if len(items) > 1 for s, rt, _ in items]
    marginal_near = marginal_weight = 0.0
    for key, rts in extras.items():
        items = consensus[key]
        if len(items) < 2:
            continue
        marginal = [leave_one_out(items, s, rt) for s, rt, q in items if q > 0.005]
        if marginal:
            marginal_near += len(rts) * sum(marginal) / len(marginal)
            marginal_weight += len(rts)
    random.seed(1)
    keys = sorted(consensus)
    shuffled = sum(near(rt, median[random.choice(keys)]) for rts in extras.values() for rt in rts)
    print(f"accepted at 1%: before {len(before)}, after {len(after)}")
    print(f"extra spectra on peptides accepted before in the run: {n_extra} (plus {new_peptides} on others)")
    if not n_extra:
        return
    print(f"within {args.window:g} s of the peptide's median RT:")
    print(f"  extras     {extra_near / n_extra:.1%}")
    print(f"  marginal   {marginal_near / max(marginal_weight, 1):.1%} (same peptides, q > 0.5% before)")
    print(f"  consensus  {sum(loo) / max(len(loo), 1):.1%} (leave-one-out)")
    print(f"  shuffled   {shuffled / n_extra:.1%}")


if __name__ == "__main__":
    main()
