"""Which SpectrumAnalysis float slots each engine fills, with examples."""
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
        if "programName" in x.fields: progs[x.fields["id"]] = f'{x.fields["programName"]} {x.fields["programVersion"]}'
        stack.extend(x.fields.values())
        for a in x.annotations.values(): stack.extend(a)
sa = records("SpectrumAnalysis", 180)
by = collections.defaultdict(lambda: collections.defaultdict(list)); count = collections.Counter()
for rid, b in sa.items():
    p = b[8]; count[p] += 1
    for o in range(24, 180, 4):
        raw = b[o:o+4]
        if raw in (b'\xff\xff\xff\x7f', b'\x00\x00\x00\x00'): continue
        v = struct.unpack('<f', raw)[0]; iv = struct.unpack('<i', raw)[0]
        by[p][o].append(v if abs(v) > 1e-6 and abs(v) < 1e9 else f"i{iv}")
print("programs", progs)
for p in sorted(by):
    print(f"-- program {p} ({progs.get(p)}): {count[p]} analyses")
    for o in sorted(by[p]):
        vals = by[p][o]
        print(f"   {o:3}: {len(vals):7} e.g. {[round(v,3) if isinstance(v,float) else v for v in vals[:5]]}")
