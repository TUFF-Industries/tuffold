"""Sibling counts in another rescorer: rebuild Tuffold's peptide-level pooling inputs from a
validated project and give them to mokapot, with the decoys' sibling counts done four ways:
none (no sibling column), unpaired (decoys get 0, as plain pooling does), paired (a decoy counts
its named target protein's identified peptides, `rev_X` -> X), and imputed (each decoy, in the
order of a hash of its id, takes the next quantile of the counts of targets scoring at or below the
decoys' median). Tells whether pairing, not Tuffold's own model, is what keeps the decoys honest.

Peptides are Tuffold's (one per modified form). Features: best and second-best spectrum score
(Tuffold's spectrum scores; the lowest score when there is one spectrum), log(1 + other spectra),
other charges, and log(1 + siblings). Identified = best spectrum passing 1% by target-decoy on the
best spectrum score; a peptide's own sequence never counts. Writes OUT/<variant>.tsv (mokapot
score, is_decoy, sequence) for entrapment_seq.py's tsv mode, and prints the peptides each variant
passes at 1%. Read-only on the project.

usage (mokapot 0.10, pandas): python tools/crosscheck/mokapot_siblings.py PROJECT.tuff OUT
"""

import collections
import hashlib
import math
import os
import sqlite3
import sys

import numpy as np
import pandas as pd

np.float_ = np.float64  # mokapot 0.10 predates NumPy 2
import mokapot  # noqa: E402

DECOY = ("reverse_", "rev_", "REV_", "DECOY_")


def q_values(scores, decoys):
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    q, d, t = [0.0] * len(scores), 0, 0
    for i in order:
        d += decoys[i]
        t += 1 - decoys[i]
        q[i] = min(1.0, (d + 1) / max(t, 1))
    low = 1.0
    for i in reversed(order):
        low = min(low, q[i])
        q[i] = low
    return q


def mix(i):
    return int.from_bytes(hashlib.sha256(str(i).encode()).digest()[:8], "little")


def main():
    path, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    spectra = collections.defaultdict(list)
    for pid, score, charge in db.execute(
        "SELECT r.peptide_id, r.score, p.charge FROM spectrum_results r "
        "JOIN psms p ON p.spectrum_id = r.spectrum_id AND p.peptide_id = r.peptide_id AND p.rank = 1"
    ):
        spectra[pid].append((score, charge))
    info = {pid: (seq, dec) for pid, seq, dec in db.execute(
        "SELECT r.peptide_id, pe.sequence, r.is_decoy FROM peptide_results r JOIN peptides pe ON pe.id = r.peptide_id")}
    ids = sorted(set(spectra) & set(info))
    proteins = {i: (a, d) for i, a, d in db.execute("SELECT id, accession, is_decoy FROM proteins")}
    by_name = {a: i for i, (a, d) in proteins.items() if not d}
    pair = {}
    for i, (a, d) in proteins.items():
        if d:
            for label in DECOY:
                if a.startswith(label) and a[len(label):] in by_name:
                    pair[i] = by_name[a[len(label):]]
    proteins_of = collections.defaultdict(set)
    for pid, prot in db.execute(
        "SELECT DISTINCT p.peptide_id, pp.protein_id FROM psm_proteins pp JOIN psms p ON p.id = pp.psm_id"
    ):
        proteins_of[pid].add(prot)

    lowest = min(s for v in spectra.values() for s, _ in v)
    best = [max(s for s, _ in spectra[i]) for i in ids]
    decoys = [info[i][1] for i in ids]
    identified = collections.defaultdict(set)
    for i, q in zip(ids, q_values(best, decoys)):
        if not info[i][1] and q <= 0.01:
            for prot in proteins_of[i]:
                if not proteins[prot][1]:
                    identified[prot].add(info[i][0])

    def siblings(prot, own):
        seqs = identified.get(prot, ())
        return len(seqs) - (own in seqs)

    target_counts, paired = [], []
    for i in ids:
        seq, dec = info[i]
        if dec:
            target_counts.append(0)
            paired.append(max((siblings(pair[p], seq) for p in proteins_of[i] if p in pair), default=0))
        else:
            c = max((siblings(p, seq) for p in proteins_of[i] if not proteins[p][1]), default=0)
            target_counts.append(c)
            paired.append(c)
    median = sorted(b for b, d in zip(best, decoys) if d)[sum(decoys) // 2]
    low = sorted(c for c, b, d in zip(target_counts, best, decoys) if not d and b <= median)
    order = sorted((mix(i), k) for k, i in enumerate(ids) if decoys[k])
    imputed = list(target_counts)
    for n, (_, k) in enumerate(order):
        imputed[k] = low[int((n + 0.5) * len(low) / len(order))]

    base = pd.DataFrame({
        "specid": ids,
        "label": [not d for d in decoys],
        "peptide": [info[i][0] for i in ids],
        "best": best,
        "second": [sorted((s for s, _ in spectra[i]), reverse=True)[1] if len(spectra[i]) > 1 else lowest for i in ids],
        "others": [math.log1p(len(spectra[i]) - 1) for i in ids],
        "charges": [len({c for _, c in spectra[i]}) - 1 for i in ids],
    })
    variants = {"none": None, "unpaired": target_counts, "paired": paired, "imputed": imputed}
    for name, counts in variants.items():
        df = base.copy()
        feats = ["best", "second", "others", "charges"]
        if counts is not None:
            df["siblings"] = [math.log1p(c) for c in counts]
            feats.append("siblings")
        data = mokapot.LinearPsmDataset(
            psms=df, target_column="label", spectrum_columns="specid",
            peptide_column="peptide", feature_columns=feats, copy_data=True)
        np.random.seed(1)
        results, _ = mokapot.brew(data)
        psms = results.confidence_estimates["psms"]
        with open(f"{out}/{name}.tsv", "w") as f:
            for table, dec in ((psms, 0), (results.decoy_confidence_estimates["psms"], 1)):
                for s, p in zip(table["mokapot score"], table["peptide"]):
                    f.write(f"{s}\t{dec}\t{p}\n")
        print(f"{name}\t{int((psms['mokapot q-value'] <= 0.01).sum())} peptides at 1% (mokapot)", flush=True)


if __name__ == "__main__":
    main()
