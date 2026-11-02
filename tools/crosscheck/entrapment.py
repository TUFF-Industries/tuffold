"""Entrapment estimates of FDR for validated Tuffold projects and TPP pepXML,
searched against a database from `entrapment_fasta.py` (read-only).

A peptide is entrapment when every non-decoy protein it maps to is
`ENTRAP_`; decoys are counted by Tuffold's own labels (`tuff`) or every
protein `rev_`/`reverse_` (`tpp`). The estimate is Wen et al.'s combined
one, N_E (1 + 1/r) / (N_T + N_E), with r from `entrapment_fasta.py`.

`tuff` uses Tuffold's own q-values: spectra, peptide sequences (I = L, a
sequence is entrapment when all its forms are) and protein groups (leader
`ENTRAP_`), at each --q. `tpp` ranks one best rank-1 hit per spectrum by
iProphet's or PeptideProphet's probability and applies target-decoy
competition, (d + 1)/t; peptide sequences get their best spectrum's
probability. Use `tpp` mode's yardstick on a Tuffold project with `--tdc`
to count both the same way.

Usage (Python's standard library only):
    python3 -I tools/crosscheck/entrapment.py tuff R PROJECT.tuff... [--q 0.005 0.01 0.02] [--tdc]
    python3 -I tools/crosscheck/entrapment.py tpp R FILE.pep.xml... [--analysis interprophet]
"""

import argparse
import sqlite3
import xml.etree.ElementTree as ET

PEPXML = "{http://regis-web.systemsbiology.net/pepXML}"


def estimate(entrapment, total, r):
    return entrapment * (1 + 1 / r) / max(total, 1)


def line(label, level, passing, r):
    """passing: one bool per accepted target, True for entrapment."""
    ne = sum(passing)
    print(f"{label}\t{level}\t{len(passing)}\t{ne}\t{estimate(ne, len(passing), r):.2%}")


def kind(proteins):
    real = [p for p in proteins if not p.startswith(("reverse_", "rev_"))]
    if not real:
        return "decoy"
    return "entrap" if all(p.startswith("ENTRAP_") for p in real) else "real"


def tdc(items, q=0.01):
    """items: (score, kind). Accepted targets at (d + 1)/t ≤ q, as is_entrapment."""
    items = sorted(items, key=lambda x: -x[0])
    decoys = targets = 0
    qs = []
    for _, k in items:
        decoys += k == "decoy"
        targets += k != "decoy"
        qs.append(min(1.0, (decoys + 1) / max(targets, 1)))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    return [k == "entrap" for (_, k), qv in zip(items, qs) if qv <= q and k != "decoy"]


def tuff(path, r, levels, use_tdc):
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    label = path.rsplit("/", 1)[-1]
    proteins = {}
    for peptide, accession in db.execute(
        """SELECT DISTINCT ps.peptide_id, pr.accession FROM psms ps
           JOIN psm_proteins x ON x.psm_id = ps.id JOIN proteins pr ON pr.id = x.protein_id
           WHERE ps.rank = 1"""
    ):
        proteins.setdefault(peptide, []).append(accession)
    entrap = {p for p, accs in proteins.items() if kind(accs) == "entrap"}
    spectra = db.execute(
        """SELECT r.peptide_id, r.q_value, r.score, r.is_decoy, pe.sequence
           FROM spectrum_results r JOIN peptides pe ON pe.id = r.peptide_id"""
    ).fetchall()
    peptides = db.execute(
        """SELECT r.peptide_id, r.q_value, r.score, r.is_decoy, pe.sequence
           FROM peptide_results r JOIN peptides pe ON pe.id = r.peptide_id"""
    ).fetchall()
    groups = db.execute(
        """SELECT g.q_value, pr.accession FROM protein_groups g
           JOIN proteins pr ON pr.id = g.leader_protein_id WHERE g.is_decoy = 0"""
    ).fetchall()
    if use_tdc:
        line(label, "spectra (TDC)", tdc([(s, "decoy" if d else kind(proteins[p])) for p, _, s, d, _ in spectra]), r)
        best = {}
        for p, _, s, d, seq in spectra:
            key = seq.replace("I", "L")
            if key not in best or s > best[key][0]:
                best[key] = (s, "decoy" if d else kind(proteins[p]))
        line(label, "peptide sequences (TDC, best spectrum)", tdc(best.values()), r)
        return
    for q in levels:
        line(label, f"spectra q≤{q:g}", [p in entrap for p, qv, _, d, _ in spectra if not d and qv is not None and qv <= q], r)
        forms = {}
        for p, qv, _, d, seq in peptides:
            if not d and qv is not None and qv <= q:
                forms.setdefault(seq.replace("I", "L"), []).append(p in entrap)
        line(label, f"peptide sequences q≤{q:g}", [all(v) for v in forms.values()], r)
        line(label, f"protein groups q≤{q:g}", [a.startswith("ENTRAP_") for qv, a in groups if qv is not None and qv <= q], r)


def tpp(path, r, analysis):
    best = {}
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag != PEPXML + "spectrum_query":
            continue
        key = (el.get("spectrum").rsplit(".", 3)[0], el.get("start_scan"))
        for hit in el.iter(PEPXML + "search_hit"):
            if hit.get("hit_rank") != "1":
                continue
            proteins = [hit.get("protein")] + [a.get("protein") for a in hit.findall(PEPXML + "alternative_protein")]
            for result in hit.findall(PEPXML + "analysis_result"):
                if result.get("analysis") == analysis:
                    p = float(result[0].get("probability"))
                    if key not in best or p > best[key][0]:
                        best[key] = (p, kind(proteins), hit.get("peptide").replace("I", "L"))
            break
        el.clear()
    label = path.rsplit("/", 1)[-1]
    line(label, "spectra (TDC)", tdc([(p, k) for p, k, _ in best.values()]), r)
    peptides = {}
    for p, k, seq in best.values():
        if seq not in peptides or p > peptides[seq][0]:
            peptides[seq] = (p, k)
    line(label, "peptide sequences (TDC, best spectrum)", tdc(peptides.values()), r)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=["tuff", "tpp"])
    parser.add_argument("r", type=float, help="twin / real peptides, from entrapment_fasta.py")
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--q", type=float, nargs="+", default=[0.01])
    parser.add_argument("--tdc", action="store_true", help="tuff: count as TPP is counted")
    parser.add_argument("--analysis", default="interprophet", help="tpp: or peptideprophet")
    args = parser.parse_args()
    print("file\tlevel\taccepted\tentrapment\testimate")
    for path in args.paths:
        if args.mode == "tuff":
            tuff(path, args.r, args.q, args.tdc)
        else:
            tpp(path, args.r, args.analysis)


if __name__ == "__main__":
    main()
