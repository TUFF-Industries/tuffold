"""How many correct matches does a validation leave below its 1% cut? Under target-decoy
competition a wrong match is a target or a decoy equally often, so below the cut targets minus
decoys estimates the correct ones left (all of them, and those with q 1-10%). Near 0 means rescoring
is saturated: no feature can add much without new candidates (search space) or new information.
usage: headroom.py PROJECT.tuff ...   (validated projects; spectra and peptides by their scores)"""
import sqlite3
import sys


def remaining(rows):
    rows = sorted(rows, key=lambda r: -r[0])
    d = t = 0
    qs = []
    for _, dec in rows:
        d += dec; t += not dec; qs.append(min(1.0, (d + 1) / max(t, 1)))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    accepted = sum(q <= 0.01 and not dec for (_, dec), q in zip(rows, qs))
    left = sum((1 - 2 * dec) for (_, dec), q in zip(rows, qs) if q > 0.01)
    near = sum((1 - 2 * dec) for (_, dec), q in zip(rows, qs) if 0.01 < q <= 0.10)
    return accepted, left, near


for path in sys.argv[1:]:
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    out = [path.rsplit("/", 1)[-1]]
    for level, table in (("spectra", "spectrum_results"), ("peptides", "peptide_results")):
        rows = db.execute(f"SELECT score, is_decoy FROM {table}").fetchall()
        if not rows:
            out.append(f"{level}: not validated")
            continue
        a, left, near = remaining(rows)
        out.append(f"{level} {a}: ~{left} correct left ({left / max(a, 1):.1%}), ~{near} with q 1-10%")
    print("\t".join(out), flush=True)
