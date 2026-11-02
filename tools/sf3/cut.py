"""ProteinIdentification inline accessions vs IdentificationCacheLookup.

The inline field is cut at 22 characters; this counts, for the candidates
Scaffold uses (offset 147 = 0), how a 22-character name completes from the
peptide -> protein map. Run from inside a dump.py OUTDIR:
    python3 -I $TOOLS/cut.py $TOOLS
"""
import sys, collections; sys.path.insert(0, sys.argv[1])
from modlib import *
from jser import load, Obj, Arr
import mods as M
i32 = M.i32
pi = records("ProteinIdentification", 103)
acc = lambda b: b[8:31].split(b'\x0a')[0].decode('latin-1')
cache = collections.defaultdict(set)
stack = list(load("IdentificationCacheLookup")); seen = set()
while stack:
    x = stack.pop()
    if id(x) in seen: continue
    seen.add(id(x))
    if isinstance(x, Obj):
        if "proteinAccession" in x.fields:
            for p in (x.fields["sequences"].items if x.fields["sequences"] else []):
                cache[p].add(x.fields["proteinAccession"])
        stack.extend(x.fields.values())
        for a in x.annotations.values(): stack.extend(a)
    elif isinstance(x, Arr) and isinstance(x.items, list): stack.extend(x.items)
lengths = collections.Counter(len(acc(b)) for b in pi.values())
print("inline lengths (longest five):", sorted(lengths.items())[-5:])
r = collections.Counter(); ex = collections.defaultdict(list)
for rid, b in M.pep.items():
    if b[147] != 0: continue
    s = b[4:34].split(b'\x0a')[0].decode()
    if i32(b, 35): s += M.tails[i32(b, 35)]
    p = i32(b, 122)
    if p not in pi: continue
    a = acc(pi[p])
    if len(a) < 22: r["under 22, whole"] += 1; continue
    listed = cache.get(s)
    if listed is None: k = "peptide not in the map"
    else:
        cands = [x for x in listed if x.startswith(a)]
        k = f"in the map, {min(len(cands), 2)}{'+' if len(cands) > 1 else ''} completion(s)"
        if len(cands) == 1: k += " to a -DECOY name" if cands[0].endswith("-DECOY") else " to a target"
    r[k] += 1
    if len(ex[k]) < 2: ex[k].append((s, a, sorted(listed or [])[:2]))
for k, v in sorted(r.items()): print(f"{v:8} {k}  e.g. {ex[k]}")
