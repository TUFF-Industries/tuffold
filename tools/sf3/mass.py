"""Peptide + modification mass vs precursor m/z x charge, for every PSM."""
import sys, struct, collections
sys.path.insert(0, sys.argv[1])
import mods as M  # runs the mod check too; reuse its tables
from modlib import records
AA = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=103.00919,
          L=113.08406, I=113.08406, N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259,
          M=131.04049, H=137.05891, F=147.06841, R=156.10111, Y=163.06333, W=186.07931, U=150.95364, O=237.14773)
H2O, P = 18.010565, 1.007276
import os
DEDUP = os.environ.get('DEDUP', '1') == '1'  # later modification at a position wins
sa = records("SpectrumAnalysis", 180); sp = records("Spectrum", 48)
i32 = M.i32
errs = collections.Counter(); worst = []; n = 0; skipped = collections.Counter()
for rid, b in M.pep.items():
    s = b[4:34].split(b'\x0a')[0].decode()
    if i32(b, 35): s += M.tails[i32(b, 35)]
    if any(c not in AA for c in s): skipped["odd residue"] += 1; continue
    ms = [(i32(b, 106 + 4*k), i32(b, 90 + 4*k)) for k in range(3) if i32(b, 90 + 4*k) != M.NONE]
    ms += list(zip(M.locs.get(i32(b, 118), []), M.types.get(i32(b, 102), [])))
    if DEDUP:
        seen = {}
        for p, t in ms: seen[p] = t  # later entry wins
        ms = list(seen.items())
    mass = sum(AA[c] for c in s) + H2O + sum(M.mt[t][1] for _, t in ms)
    a = sa.get(i32(b, 57))
    if a is None: skipped["no analysis"] += 1; continue
    spec = sp.get(i32(a, 4))
    if spec is None: skipped["no spectrum"] += 1; continue
    mz = struct.unpack('<d', spec[12:20])[0]; z = i32(spec, 20)
    if z <= 0: skipped["no charge"] += 1; continue
    obs = mz - P  # stored value is [M+H]+
    best = min(((obs - mass - k * 1.003355) / mass * 1e6 for k in range(-1, 4)), key=abs)
    n += 1
    bare = sum(1 for i, c in enumerate(s) if c == 'C' and i not in {p for p, _ in ms})
    if bare: errs[('bareC', abs(best) < 50)] += 1
    bucket = "<10ppm" if abs(best) < 10 else "<50ppm" if abs(best) < 50 else "<1000ppm" if abs(best) < 1000 else ">1000ppm"
    errs[bucket] += 1
    if abs(best) >= 1000: errs[('prog', a[8], s.count('C'), round((obs-mass)/57.02146, 2))] += 1
    if abs(best) >= 1000 and len(worst) < 5: worst.append((rid, s, ms, round(obs, 4), round(mass, 4), z))
print(f"PSMs checked {n}: {dict(errs)}; skipped {dict(skipped)}")
for w in worst: print("  ", w)
