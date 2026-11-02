"""Ground-truth FDP on ProteomeTools synthetic peptide pools (PRIDE PXD004732).

Each pool holds ~1,000 known synthetic peptides. The truth list is the pool's
peptides as the authors' MaxQuant search of the pool alone found them
(`peptides.txt` from the pool's `-tryptic.zip`; ~98% of each pool's sequence),
so a few real analytes (carryover, peptides MaxQuant missed) are missing from
it and read as wrong: compare methods on one search, and see `bins` for where
the floor comes from. Peptides are compared with I = L.

  truth  FASTA OUT.tsv POOL=peptides.txt ...   pool peptides and the FASTA entries holding them
  score  TRUTH.tsv FASTA tpp FILE.pep.xml peptideprophet|interprophet LABEL
  score  TRUTH.tsv FASTA model FILE.pep.xml LABEL    iProphet's own model (mean 1 - p <= 1%)
  score  TRUTH.tsv FASTA tuff PROJECT.tuff LABEL
  bins   TRUTH.tsv FILE.pep.xml                       wrong targets vs decoys by PeptideProphet probability

`score` prints, at 1% by target-decoy competition ((d + 1)/t; spectra by best
rank-1 hit, peptides by best score per sequence), the passing targets, those in
their run's pool (spectra) or any pool (peptides), other pools' peptides,
wrong ones on proteins holding a pool peptide, and related ones (sub- or
superstrings of a pool peptide: synthesis by-products, cleavage variants).
The pool of a run is read from its name (`pool_61_`). FASTA is the target
database searched; decoys in it (`reverse_`) are ignored."""
import bisect
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET

NS = "{http://regis-web.systemsbiology.net/pepXML}"


def il(seq):
    return seq.replace("I", "L")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name and not name.startswith("reverse_"):
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name and not name.startswith("reverse_"):
        yield name, "".join(seq)


def truth(fasta, out, pools):
    names, starts, parts, pos = [], [], [], 0
    for name, seq in read_fasta(fasta):
        names.append(name); starts.append(pos); parts.append(il(seq)); pos += len(seq) + 1
    big = "#".join(parts)
    with open(out, "w") as f:
        for arg in pools:
            pool, path = arg.split("=", 1)
            lines = open(path)
            header = lines.readline().rstrip("\n").split("\t")
            proteins = header.index("Proteins")
            peptides = {il(r.split("\t")[0]) for r in lines if "TUM_" in r.split("\t")[proteins]}
            found = 0
            for pep in sorted(peptides):
                hits, i = set(), big.find(pep)
                while i >= 0:
                    hits.add(names[bisect.bisect_right(starts, i) - 1]); i = big.find(pep, i + 1)
                found += bool(hits)
                f.write(f"{pool}\t{pep}\t{';'.join(sorted(hits))}\n")
            print(f"pool {pool}: {len(peptides)} peptides, {found} in the database")


def load_truth(path):
    pools, present = {}, set()
    for line in open(path):
        pool, pep, prots = line.rstrip("\n").split("\t")
        pools.setdefault(pool, set()).add(pep)
        present.update(p for p in prots.split(";") if p)
    return pools, present


def pool_of(run):
    return re.search(r"pool_(\d+)_", run).group(1)


def tdc(items):
    """(score, is_decoy, sequence, run) -> passing targets as (sequence, run)."""
    items = sorted(items, key=lambda x: -x[0])
    d = t = 0
    qs = []
    for _, dec, _, _ in items:
        d += dec; t += not dec; qs.append(min(1.0, (d + 1) / max(t, 1)))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    return [(seq, run) for (_, dec, seq, run), q in zip(items, qs) if q <= 0.01 and not dec]


def pepxml(path, analysis):
    """Best hit per spectrum: (probability, is_decoy, sequence, run)."""
    best = {}
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag != NS + "spectrum_query":
            continue
        run = el.get("spectrum").rsplit(".", 3)[0]
        for hit in el.iter(NS + "search_hit"):
            if hit.get("hit_rank") != "1":
                continue
            prots = [hit.get("protein")] + [a.get("protein") for a in hit.findall(NS + "alternative_protein")]
            dec = all(p.startswith("reverse_") for p in prots)
            for ar in hit.findall(NS + "analysis_result"):
                if ar.get("analysis") == analysis:
                    p = float(ar[0].get("probability"))
                    key = (run, el.get("start_scan"))
                    if key not in best or p > best[key][0]:
                        best[key] = (p, dec, il(hit.get("peptide")), run)
            break
        el.clear()
    return list(best.values())


def by_peptide(rows):
    best = {}
    for score, dec, seq, _ in rows:
        if seq not in best or score > best[seq][0]:
            best[seq] = (score, dec)
    return [(s, d, seq, None) for seq, (s, d) in best.items()]


def score(truth_path, fasta, mode, path, rest):
    pools, present = load_truth(truth_path)
    union = set().union(*pools.values())
    union_seq = "#".join(sorted(union))
    present_seq = "#".join(il(s) for n, s in read_fasta(fasta) if n in present)
    label = rest[-1]

    def report(level, passing):
        right = sum(s in pools[pool_of(r)] if r else s in union for s, r in passing)
        other = sum(1 for s, r in passing if r and s not in pools[pool_of(r)] and s in union)
        wrong = [s for s, _ in passing if s not in union]
        related = sum(s in union_seq or any(p in s for p in union) for s in wrong)
        on_present = sum(s in present_seq for s in wrong)
        n = len(passing)
        print(f"{label}\t{level}\t{n}\tright {right}\twrong {n - right} ({(n - right) / max(n, 1):.2%})"
              f"\tother-pool {other}\twrong on present proteins {on_present}\trelated {related}"
              f"\twrong less related and other-pool {(n - right - other - related) / max(n, 1):.2%}",
              flush=True)

    if mode == "tpp":
        rows = pepxml(path, rest[0])
        report("spectra", tdc(rows))
        report("peptides", tdc(by_peptide(rows)))
    elif mode == "model":
        accepted, false = [], 0.0
        for p, dec, seq, _ in sorted(by_peptide(pepxml(path, "interprophet")), key=lambda x: -x[0]):
            if dec:
                continue
            if (false + 1 - p) / (len(accepted) + 1) > 0.01:
                break
            false += 1 - p
            accepted.append((seq, None))
        report("peptides, own model 1%", accepted)
    else:
        db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        rows = [(s, d, il(q), run) for s, d, q, run in db.execute(
            "SELECT r.score, r.is_decoy, pe.sequence, m.name FROM spectrum_results r"
            " JOIN peptides pe ON pe.id = r.peptide_id JOIN spectra sp ON sp.id = r.spectrum_id"
            " JOIN ms_runs m ON m.id = sp.ms_run_id")]
        report("spectra", tdc(rows))
        peptides = [(s, d, il(q), None) for s, d, q in db.execute(
            "SELECT r.score, r.is_decoy, pe.sequence FROM peptide_results r"
            " JOIN peptides pe ON pe.id = r.peptide_id")]
        report("peptides (Tuffold peptide score)", tdc(by_peptide(peptides)))


def bins(truth_path, path):
    """Unlisted, unrelated targets vs decoys by PeptideProphet probability over all rank-1 hits:
    where they agree, the decoys model the errors; an excess of unlisted targets at the top is
    real analytes missing from the truth list."""
    pools, _ = load_truth(truth_path)
    union = set().union(*pools.values())
    rows = pepxml(path, "peptideprophet")
    print("probability\ttargets\tright\tother pool\trelated\tunlisted\tdecoys")
    for lo, hi in [(0.999, 2), (0.99, 0.999), (0.95, 0.99), (0.8, 0.95), (0.5, 0.8), (0.2, 0.5), (0, 0.2)]:
        b = [r for r in rows if lo <= r[0] < hi]
        targets = [r for r in b if not r[1]]
        right = sum(s in pools[pool_of(run)] for _, _, s, run in targets)
        other = sum(s in union and s not in pools[pool_of(run)] for _, _, s, run in targets)
        related = sum(s not in union and any(s in u or u in s for u in union) for _, _, s, _ in targets)
        unlisted = len(targets) - right - other - related
        print(f"{lo}-{min(hi, 1)}\t{len(targets)}\t{right}\t{other}\t{related}\t{unlisted}\t{len(b) - len(targets)}")


if __name__ == "__main__":
    command, args = sys.argv[1], sys.argv[2:]
    if command == "truth":
        truth(args[0], args[1], args[2:])
    elif command == "score":
        score(args[0], args[1], args[2], args[3], args[4:])
    elif command == "bins":
        bins(args[0], args[1])
    else:
        sys.exit(__doc__)
