"""Build a one-proteome entrapment database (after Wen et al., Nat Methods 2025).

Every protein whose header contains TEXT (all proteins without --only) gets a
shuffled twin named ENTRAP_<accession>. K/R, and a P right after them, stay
in place and residues are shuffled only within each tryptic segment, so every
tryptic peptide keeps its length and composition (and so its mass). A twin
peptide that also occurs among the real peptides (I = L) has its segment
redrawn; any left are reported, and count as real when scoring. Decoys in the
input (any `rev`/`reverse`/`decoy` label) are dropped: let the search engine
or Tuffold's built-in search make them for the combined database. Seeded.

Prints r, the ratio of twin to real peptides, which `entrapment.py` needs.

Usage (Python's standard library only):
    python3 -I tools/crosscheck/entrapment_fasta.py IN.fasta OUT.fasta \
        [--only "OS=Saccharomyces"] [--seed 20251011]
"""

import argparse
import random
import re

MIN_LEN, MAX_LEN, MISSED = 7, 50, 2  # Tuffold's search defaults
DECOY = re.compile(r"^(?:[a-z]{2}\|)?(?:rev|reverse|decoy)[_\-|]", re.IGNORECASE)


def read(path):
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name:
        yield name, "".join(seq)


def segments(seq):
    """Tryptic pieces: cut after K/R unless P follows."""
    out, start = [], 0
    for i, a in enumerate(seq):
        if a in "KR" and (i + 1 == len(seq) or seq[i + 1] != "P"):
            out.append((start, i + 1))
            start = i + 1
    if start < len(seq):
        out.append((start, len(seq)))
    return out


def peptides(seq):
    segs = segments(seq)
    for i in range(len(segs)):
        for j in range(i, min(i + MISSED + 1, len(segs))):
            p = seq[segs[i][0] : segs[j][1]]
            if MIN_LEN <= len(p) <= MAX_LEN:
                yield p


def fixed(seq, k):
    return k == 0 or seq[k] in "KR" or (seq[k] == "P" and seq[k - 1] in "KR")


def shuffle_segment(seq, a, b, rng):
    free = [k for k in range(a, b) if not fixed(seq, k)]
    residues = [seq[k] for k in free]
    for _ in range(100):
        rng.shuffle(residues)
        out = list(seq[a:b])
        for k, r in zip(free, residues):
            out[k - a] = r
        # No new P right after a K/R: it would block a cleavage.
        if not any(out[k - a] == "P" and seq[k - 1] in "KR" for k in free if k > 0):
            return "".join(out)
    return seq[a:b]


def twin(seq, real, rng):
    out = list(seq)
    for a, b in segments(seq):
        out[a:b] = shuffle_segment(seq, a, b, rng)
    twin_seq, redrawn = "".join(out), 0
    for _ in range(20):
        hits = {p for p in peptides(twin_seq) if p.replace("I", "L") in real}
        if not hits:
            break
        redrawn += len(hits)
        for p in sorted(hits):
            i = twin_seq.find(p)
            for a, b in segments(twin_seq):
                if i <= a < i + len(p):
                    twin_seq = twin_seq[:a] + shuffle_segment(seq, a, b, rng) + twin_seq[b:]
    return twin_seq, redrawn


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("source")
    parser.add_argument("out")
    parser.add_argument("--only", help="twin only proteins whose header contains this")
    parser.add_argument("--seed", type=int, default=20251011)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    entries = [(h, s) for h, s in read(args.source) if not DECOY.match(h)]
    real = {p.replace("I", "L") for _, s in entries for p in peptides(s)}
    chosen = [(h, s) for h, s in entries if args.only is None or args.only in h]
    twins, redrawn = [], 0
    for header, seq in chosen:
        twin_seq, n = twin(seq, real, rng)
        redrawn += n
        accession = header.split()[0]
        twins.append((f"ENTRAP_{accession} shuffled entrapment copy of {accession}", twin_seq))
    twin_peptides = {p.replace("I", "L") for _, s in twins for p in peptides(s)}
    left = len(twin_peptides & real)
    with open(args.out, "w") as out:
        for header, seq in entries + twins:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
    print(f"real entries {len(entries)}, twins {len(twins)}")
    print(f"real peptides {len(real)}, twin peptides {len(twin_peptides - real)}")
    print(f"r = {len(twin_peptides - real) / len(real):.4f}")
    print(f"collisions redrawn {redrawn}, still shared {left} (scored as real)")


if __name__ == "__main__":
    main()
