"""For each field: take each spectrum's top match by that field, then target-vs-decoy AUC and decoy share."""
import sys, struct, collections, bisect
sys.path.insert(0, sys.argv[1])
from modlib import records, load, Obj
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
sa = records("SpectrumAnalysis", 180); pep = records("PeptideIdentification", 160); pi = records("ProteinIdentification", 103)
f = lambda b, o: struct.unpack('<f', b[o:o+4])[0]
i = lambda b, o: struct.unpack('<i', b[o:o+4])[0]
groups = collections.defaultdict(list)
for rid, b in pep.items():
    a = sa.get(i(b, 57)); p = pi.get(i(b, 122))
    if a is None or p is None: continue
    groups[(i(a, 4), a[8])].append((b, b"DECOY" in p[8:28].upper()))
def auc(t, d):
    d = sorted(d)
    return sum(bisect.bisect_left(d, v) + 0.5 * (bisect.bisect_right(d, v) - bisect.bisect_left(d, v)) for v in t) / (len(t) * len(d))
out = collections.defaultdict(dict)
for o in (49, 66, 70, 74, 78, 82, 152):
    for sign in (1, -1):
        per = collections.defaultdict(lambda: ([], []))
        for (spec, prog), bs in groups.items():
            b, dec = max(bs, key=lambda x: sign * f(x[0], o))
            per[progs.get(prog, prog)][dec].append(sign * f(b, o))
        for name, (t, d) in per.items():
            if t and d: out[name][f"{o}{'+' if sign > 0 else '-'}"] = (round(auc(t, d), 3), round(len(d) / (len(t) + len(d)), 3))
for name, r in out.items():
    print(name, sorted(r.items(), key=lambda kv: -kv[1][0])[:6])
