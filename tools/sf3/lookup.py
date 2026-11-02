"""IdentificationCacheLookup: protein -> peptide sequences. Are they in the protein's sequence?"""
import sys, collections
sys.path.insert(0, sys.argv[1])
from jser import load, Obj, Arr
from bbload import bb_records
seqs = {o.fields["accessionNumber"]: o.fields["sequence"].replace("I", "L") for _, _, o in bb_records("BBProteinSequence").values() if not isinstance(o, Exception)}
c = collections.Counter(); peps_per = []
stack = list(load("IdentificationCacheLookup")); seen = set()
while stack:
    x = stack.pop()
    if id(x) in seen: continue
    seen.add(id(x))
    if isinstance(x, Obj):
        if "proteinAccession" in x.fields:
            acc = x.fields["proteinAccession"]; ps = x.fields["sequences"].items if x.fields["sequences"] else []
            c["proteins"] += 1; peps_per.append(len(ps))
            if acc in seqs:
                for p in ps: c["in sequence" if p.replace("I", "L") in seqs[acc] else "NOT in sequence"] += 1
            else: c["protein without sequence"] += 1
            c["decoy proteins"] += "DECOY" in acc.upper()
        stack.extend(x.fields.values())
        for a in x.annotations.values(): stack.extend(a)
    elif isinstance(x, Arr) and isinstance(x.items, list): stack.extend(x.items)
print(dict(c), "peptides per protein: max", max(peps_per or [0]), "mean", round(sum(peps_per) / max(1, len(peps_per)), 1))
