import sys, struct, itertools
sys.path.insert(0, sys.argv[1])
from modlib import records
pep = records("PeptideIdentification", 160)
f = lambda b, o: struct.unpack('<f', b[o:o+4])[0]
O = (66, 70, 74, 78, 82, 152)
rows = [[f(b, o) for o in O] for b in list(pep.values())[:20000]]
rows = [r for r in rows if all(abs(x) < 1e30 for x in r)]
for (i, a), (j, b), (k, c) in itertools.permutations(list(enumerate(O)), 3):
    if j > k: continue
    for sign in (1, -1):
        ok = sum(abs(r[i] - (r[j] - sign * r[k] if sign == 1 else r[j] + r[k])) < 1e-3 * max(1, abs(r[i])) for r in rows)
        if ok > 0.95 * len(rows): print(f"{a} = {b} {'-' if sign == 1 else '+'} {c}  ({ok}/{len(rows)})")
for (i, a), (j, b) in itertools.permutations(list(enumerate(O)), 2):
    ok = sum(abs(r[i] - r[j]) < 1e-6 for r in rows)
    if ok > 0.95 * len(rows) and i < j: print(f"{a} == {b}  ({ok}/{len(rows)})")
