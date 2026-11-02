"""Entrapment estimates by sequence, for databases whose twins may carry their real protein's name
(embedded entrapment, see embed_pepxml.py). A passing target peptide is entrapment when its sequence
(I = L) occurs in some ENTRAP_ entry of the FASTA and in no real one; decoys are reverse_/rev_.
Spectra: best rank-1 hit per spectrum; peptides: best score per sequence; target-decoy (d + 1)/t
<= 1%. Prints targets, twins, the lower bound N_E/(N_T + N_E) (twins are certainly wrong) and the
combined estimate N_E(1 + 1/r)/(N_T + N_E). tpp-model cuts instead where iProphet's own mean
(1 - p) over target sequences reaches 1%. Read-only, stdlib only.

Then, to show which twins the decoys fail to count: decoys passing, and passing twins per decoy as a
multiple of r/(1 + r), the value if twins and decoys were exchangeable (1.00x; above it decoys
undercount what twins stand for). Twins keep their real protein's tryptic segments
(entrapment_fasta.py), so each twin peptide pairs with the real one at the same offset, and that
ratio is split twice: near-copies (at most 4 positions differ) vs far twins, and shadow twins (their
real counterpart passes, so they can take its spectra, which no wrong real peptide can) vs free
twins, each against its own share of the twin peptides.

usage: python3 -I entrapment_seq.py FASTA tpp FILE.pep.xml interprophet|peptideprophet LABEL R
       python3 -I entrapment_seq.py FASTA tpp-model FILE.pep.xml LABEL R
       python3 -I entrapment_seq.py FASTA tuff PROJECT.tuff LABEL R
       python3 -I entrapment_seq.py FASTA tsv PEPTIDES.tsv LABEL R
(tsv: score, is_decoy, sequence per peptide, as mokapot_siblings.py writes them)
"""
import sqlite3
import sys
import xml.etree.ElementTree as ET

NS = "{http://regis-web.systemsbiology.net/pepXML}"
DECOY = ("reverse_", "rev_")


def read(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        yield name, "".join(seq)


MIN_LEN, MAX_LEN, MISSED = 7, 50, 2  # entrapment_fasta.py's digest
NEAR = 4  # positions a near-copy twin differs in, at most


def spans(seq):
    """(start, end) of the digest's peptides: cut after K/R unless P follows."""
    cuts = [0] + [i + 1 for i, a in enumerate(seq[:-1]) if a in "KR" and seq[i + 1] != "P"] + [len(seq)]
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + MISSED + 2, len(cuts))):
            if MIN_LEN <= cuts[j] - cuts[i] <= MAX_LEN:
                yield cuts[i], cuts[j]


class Twins:
    def __init__(self, fasta):
        twin, real = {}, {}
        for name, seq in read(fasta):
            if not name.startswith(DECOY):
                seq = seq.replace("I", "L")
                if name.startswith("ENTRAP_"):
                    twin[name[len("ENTRAP_"):]] = seq
                else:
                    real[name] = seq
        self.twin_all, self.real_all, self.cache = "#".join(twin.values()), "#".join(real.values()), {}
        self.real_peptides = {s[a:b] for s in real.values() for a, b in spans(s)}
        self.diff, self.counterparts = {}, {}
        for name, t in twin.items():
            s = real.get(name)
            if s is None or len(s) != len(t):
                continue
            for a, b in spans(t):
                p = t[a:b]
                if p in self.real_peptides:
                    continue
                d = sum(x != y for x, y in zip(p, s[a:b]))
                self.diff[p] = min(d, self.diff.get(p, d))
                self.counterparts.setdefault(p, set()).add(s[a:b])

    def entrap(self, seq):
        if seq not in self.cache:
            self.cache[seq] = seq in self.twin_all and seq not in self.real_all
        return self.cache[seq]

    def near(self, seq):
        return self.diff.get(seq, NEAR + 1) <= NEAR

    def shadow(self, seq, passing):
        return bool(self.counterparts.get(seq, set()) & passing)

    def share(self, keep):
        """Twin peptides passing `keep`, per real peptide (r for that kind of twin)."""
        return sum(1 for p in self.counterparts if keep(p)) / len(self.real_peptides)


def tdc(items):
    """(score, is_decoy, sequence) -> (sequences of targets passing 1%, decoys passing)."""
    items = sorted(items, key=lambda x: -x[0])
    d = t = 0
    qs = []
    for _, dec, _ in items:
        d += dec
        t += 1 - dec
        qs.append(min(1.0, (d + 1) / max(t, 1)))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    passing = [(dec, seq) for (_, dec, seq), q in zip(items, qs) if q <= 0.01]
    return [seq for dec, seq in passing if not dec], sum(dec for dec, _ in passing)


def report(twins, label, level, passing, decoys, r):
    entrap = [s for s in passing if twins.entrap(s)]
    ne, n = len(entrap), max(len(passing), 1)
    line = f"{label}\t{level}\t{len(passing) - ne}\t{ne}\t{ne / n:.2%}\t{ne * (1 + 1 / r) / n:.2%}"
    if decoys:
        real = {s for s in passing if not twins.entrap(s)}

        def per_decoy(count, share):
            return f"{count / decoys / (share / (1 + r)):.2f}x"

        far = [s for s in entrap if not twins.near(s)]
        free = [s for s in entrap if not twins.shadow(s, real)]
        r_far = twins.share(lambda p: not twins.near(p))
        r_free = twins.share(lambda p: not twins.shadow(p, real))
        line += (
            f"\tdecoys {decoys}\ttwins/decoy {per_decoy(ne, r)}"
            f"\tnear {ne - len(far)}\tfar {len(far)} {per_decoy(len(far), r_far)}"
            f"\tshadow {ne - len(free)}\tfree {len(free)} {per_decoy(len(free), r_free)}"
        )
    print(line, flush=True)


def tpp_hits(path, analysis):
    """Best rank-1 (probability, is_decoy, sequence) per spectrum."""
    best = {}
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag != NS + "spectrum_query":
            continue
        key = (el.get("spectrum").rsplit(".", 3)[0], el.get("start_scan"))
        for hit in el.iter(NS + "search_hit"):
            if hit.get("hit_rank") != "1":
                continue
            prots = [hit.get("protein")] + [a.get("protein") for a in hit.findall(NS + "alternative_protein")]
            dec = int(all(p.startswith(DECOY) for p in prots))
            for ar in hit.findall(NS + "analysis_result"):
                if ar.get("analysis") == analysis:
                    p = float(ar[0].get("probability"))
                    if key not in best or p > best[key][0]:
                        best[key] = (p, dec, hit.get("peptide").replace("I", "L"))
            break
        el.clear()
    return best


def best_per_sequence(rows):
    out = {}
    for score, dec, seq in rows:
        if seq not in out or score > out[seq][0]:
            out[seq] = (score, dec)
    return [(s, d, q) for q, (s, d) in out.items()]


def main():
    fasta, mode, path = sys.argv[1:4]
    label, r = sys.argv[-2], float(sys.argv[-1])
    twins = Twins(fasta)
    if mode == "tpp":
        spectra = list(tpp_hits(path, sys.argv[4]).values())
        report(twins, label, "spectra", *tdc(spectra), r)
        report(twins, label, "peptides", *tdc(best_per_sequence(spectra)), r)
    elif mode == "tpp-model":
        peptides = sorted(best_per_sequence(tpp_hits(path, "interprophet").values()), key=lambda x: -x[0])
        accepted, false = [], 0.0
        for p, dec, seq in peptides:
            if dec:
                continue
            if (false + 1 - p) / (len(accepted) + 1) > 0.01:
                break
            false += 1 - p
            accepted.append(seq)
        report(twins, label, "peptides (own model)", accepted, None, r)
    elif mode == "tsv":
        with open(path) as f:
            rows = [(float(s), int(d), q.replace("I", "L")) for s, d, q in (line.split() for line in f)]
        report(twins, label, "peptides", *tdc(best_per_sequence(rows)), r)
    else:
        db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        spectra = [
            (s, d, q.replace("I", "L"))
            for s, d, q in db.execute(
                "SELECT r.score, r.is_decoy, pe.sequence FROM spectrum_results r "
                "JOIN peptides pe ON pe.id = r.peptide_id"
            )
        ]
        report(twins, label, "spectra", *tdc(spectra), r)
        peptides = [
            (s, d, q.replace("I", "L"))
            for s, d, q in db.execute(
                "SELECT r.score, r.is_decoy, pe.sequence FROM peptide_results r "
                "JOIN peptides pe ON pe.id = r.peptide_id"
            )
        ]
        report(twins, label, "peptides (Tuffold peptide score)", *tdc(best_per_sequence(peptides)), r)


if __name__ == "__main__":
    main()
