"""Check the modification layout of PeptideIdentification in one unpacked sf3."""
import sys, struct, collections
sys.path.insert(0, sys.argv[1])
from modlib import *
NONE = 0x7fffffff
def i32(b, o): return struct.unpack('<i', b[o:o+4])[0]

def strings(name):
    recs = records(name, 68); out = {}
    for rid, b in recs.items():
        out.setdefault(i32(b, 0), []).append(b[4:64].split(b'\x02')[0])
    return {k: b"".join(v).decode() for k, v in out.items()}

def int_lists(name):
    out = {}
    for rid, b in records(name, 256).items():
        ints = struct.unpack('<64i', b)
        vals = []
        for v in ints[1:]:
            if v == INT_MIN: break
            vals.append(v)
        out[ints[0]] = vals
    return out

mt = mod_types(); tails = strings("PeptideSequence")
types, locs = int_lists("ModificationTypeIds"), int_lists("ModificationLocations")
pep = records("PeptideIdentification", 160)
bad_id = sum(1 for rid, b in pep.items() if i32(b, 0) != rid)
stats = collections.Counter(); used = collections.Counter(); problems = []
for rid, b in pep.items():
    s = b[4:34].split(b'\x0a')[0].decode()
    tail = i32(b, 35)
    if tail:
        if tail not in tails: problems.append((rid, "missing tail", tail)); continue
        s += tails[tail]
    ms = [(i32(b, 106 + 4*k), i32(b, 90 + 4*k)) for k in range(3) if i32(b, 90 + 4*k) != NONE]
    ot, ol = i32(b, 102), i32(b, 118)
    if ot or ol:
        stats["overflow"] += 1
        if len(types.get(ot, [])) != len(locs.get(ol, [1])): problems.append((rid, "overflow lengths", ot, ol))
        ms += list(zip(locs.get(ol, []), types.get(ot, [])))
    if ms: stats["modified"] += 1
    for pos, t in ms:
        if t not in mt: problems.append((rid, "unknown type", t)); continue
        name, mass, aa, where = mt[t]; used[(t, name, aa)] += 1
        if not (0 <= pos < len(s)): problems.append((rid, s, "pos out of range", pos, name)); continue
        ok = s[pos] in aa or ("n" in aa and pos == 0) or ("c" in aa and pos == len(s) - 1)
        stats["ok" if ok else "mismatch"] += 1
        if not ok: problems.append((rid, s, pos, s[pos], name, aa))
print(f"peptides {len(pep)}, bad ids {bad_id}, {dict(stats)}")
print("types used:", used.most_common(12))
print("problems", len(problems), problems[:8])
