"""Is byte 147 == 0 the top-scoring match per spectrum (per engine's primary score)?"""
import sys, struct, collections
sys.path.insert(0, sys.argv[1])
from modlib import records, load, Obj
PRIMARY = {"Mascot": 70, "Sequest": 74, "MaxQuant": 66, "X! Tandem": 70, "Spectrum Mill": 70}
progs = {}
stack = list(load("DatabaseSearch")); seen = set()
while stack:
    x = stack.pop()
    if id(x) in seen: continue
    seen.add(id(x))
    if isinstance(x, Obj):
        if "programName" in x.fields: progs[x.fields["id"]] = x.fields["programName"]
        stack.extend(x.fields.values())
        for a in x.annotations.values(): stack.extend(a)
sa = records("SpectrumAnalysis", 180); pep = records("PeptideIdentification", 160)
f = lambda b, o: struct.unpack('<f', b[o:o+4])[0]
i = lambda b, o: struct.unpack('<i', b[o:o+4])[0]
groups = collections.defaultdict(list)
for rid, b in pep.items():
    a = sa.get(i(b, 57))
    if a is None: continue
    groups[(i(a, 4), a[8])].append(b)
c = collections.Counter()
for (spec, prog), bs in groups.items():
    name = progs.get(prog, prog)
    chosen = [b for b in bs if b[147] in (0, 1)]
    if not chosen or len(bs) < 2: continue
    for o in (49, 66, 70, 74, 78, 82, 152):
        best, worst = max(f(b, o) for b in bs), min(f(b, o) for b in bs)
        c[(name, o, "max")] += all(f(b, o) == best for b in chosen)
        c[(name, o, "min")] += all(f(b, o) == worst for b in chosen)
    c[(name, "n")] += 1
names = sorted({k[0] for k in c})
for n in names:
    tot = c[(n, "n")]
    print(n, tot, {f"{o}{m}": round(c[(n, o, m)] / tot, 3) for o in (49, 66, 70, 74, 78, 82, 152) for m in ("max", "min") if tot and c[(n, o, m)] / tot > 0.6})
