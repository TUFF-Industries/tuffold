"""Target-vs-decoy AUC of each PeptideIdentification float, per engine (0.5 = no signal)."""
import sys, struct, collections
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
sa = records("SpectrumAnalysis", 180); pep = records("PeptideIdentification", 160)
pi = records("ProteinIdentification", 103)
OFFS = [49, 53, 66, 70, 74, 78, 82, 152]
data = collections.defaultdict(lambda: collections.defaultdict(list))
for rid, b in pep.items():
    a = sa.get(struct.unpack('<i', b[57:61])[0]); p = pi.get(struct.unpack('<i', b[122:126])[0])
    if a is None or p is None: continue
    decoy = b"DECOY" in p[8:28].upper()
    top = struct.unpack('<i', a[20:24])[0] == rid
    name = progs.get(a[8], a[8])
    for o in OFFS:
        v = struct.unpack('<f', b[o:o+4])[0]
        if v == v and abs(v) < 1e30: data[(name, top)][o].append((v, decoy))
def auc(pairs):
    t = sorted(v for v, d in pairs if not d); dd = sorted(v for v, d in pairs if d)
    if not t or not dd: return None
    import bisect
    s = sum(bisect.bisect_left(dd, v) + 0.5 * (bisect.bisect_right(dd, v) - bisect.bisect_left(dd, v)) for v in t)
    return s / (len(t) * len(dd))
for (name, top), offs in sorted(data.items(), key=str):
    n = len(offs[49]); dec = sum(d for _, d in offs[49])
    print(f"{name} {'first-listed' if top else 'other'}: {n} PSMs, {dec} decoy", {o: (round(a, 3) if (a := auc(v)) is not None else None) for o, v in offs.items()})
