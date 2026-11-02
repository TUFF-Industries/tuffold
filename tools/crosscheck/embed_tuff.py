"""Embedded entrapment inside a Tuffold project: re-point the standard design's twins (ENTRAP_X)
and their decoys (reverse_ENTRAP_X) to real proteins, in a copy of the project. 'same' renames
them to X / reverse_X, as embed_pepxml.py does to a pepXML (and gives the same results). 'cross'
renames them to the next real protein by length (cyclic) and its reverse: twins still land in
present proteins at the natural rate, but a twin that takes its own real peptide's spectra (same
composition) no longer gets that peptide's protein's siblings. Validate the copy, then score it
with entrapment_seq.py. Stdlib only.
usage: python3 -I embed_tuff.py same|cross FASTA IN.tuff OUT.tuff"""
import shutil
import sqlite3
import sys


def main():
    mode, fasta, src, dst = sys.argv[1:5]
    seqs, name = {}, None
    for line in open(fasta):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    length = {k: sum(map(len, v)) for k, v in seqs.items()}
    real = sorted(
        (k for k in seqs if not k.startswith(("ENTRAP_", "reverse_", "rev_")) and "ENTRAP_" + k in seqs),
        key=lambda k: (length[k], k),
    )
    if mode == "cross":
        to = {real[i]: real[(i + 1) % len(real)] for i in range(len(real))}
    else:
        to = {k: k for k in real}
    shutil.copy(src, dst)
    db = sqlite3.connect(dst)
    ids = {a: i for i, a in db.execute("SELECT id, accession FROM proteins")}

    def protein(accession, decoy):
        if accession not in ids:
            cur = db.execute("INSERT INTO proteins (accession, is_decoy) VALUES (?, ?)", (accession, decoy))
            ids[accession] = cur.lastrowid
        return ids[accession]

    moved = 0
    for accession, old in list(ids.items()):
        if accession.startswith("ENTRAP_"):
            target, decoy = to.get(accession[7:]), 0
        elif accession.startswith("reverse_ENTRAP_"):
            target, decoy = to.get(accession[15:]), 1
            target = target and "reverse_" + target
        else:
            continue
        if target is None:
            print(f"left alone (no real protein in the FASTA): {accession}")
            continue
        new = protein(target, decoy)
        # A PSM already linked to the new protein at the same position keeps that link only.
        moved += db.execute("UPDATE OR IGNORE psm_proteins SET protein_id = ? WHERE protein_id = ?", (new, old)).rowcount
        db.execute("DELETE FROM psm_proteins WHERE protein_id = ?", (old,))
    db.commit()
    print(f"{mode}: {moved} protein links moved; {len(real)} real proteins have twins")


if __name__ == "__main__":
    main()
