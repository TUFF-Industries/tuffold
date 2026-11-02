"""Profile PeptideIdentification score fields per engine."""
import sys, struct, collections, statistics
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
OFFS = [45, 49, 53, 66, 70, 74, 78, 82, 126, 130, 134, 138, 142, 146, 152, 156]
by = collections.defaultdict(lambda: collections.defaultdict(list)); bytes_ = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
for rid, b in pep.items():
    a = sa.get(struct.unpack('<i', b[57:61])[0])
    name = progs.get(a[8], a[8]) if a else "?"
    for o in (39, 40, 41, 42, 43, 44, 65): bytes_[name][o][b[o]] += 1
    for o in OFFS:
        raw = b[o:o+4]
        if raw in (b'\xff\xff\xff\x7f',): continue
        by[name][o].append(struct.unpack('<f', raw)[0])
for name in by:
    print(f"-- {name}: {len(by[name][49])} PSMs")
    print("   bytes:", {o: dict(c.most_common(4)) for o, c in bytes_[name].items()})
    for o in OFFS:
        v = [x for x in by[name][o] if x == x and abs(x) < 1e30]
        if not v: continue
        nz = sum(1 for x in v if x != 0)
        if nz == 0: continue
        print(f"   {o:3}: n={len(v)} nonzero={nz} min={min(v):.4g} med={statistics.median(v):.4g} max={max(v):.4g}")
