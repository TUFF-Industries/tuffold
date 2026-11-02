import sys, struct, collections; sys.path.insert(0, sys.argv[1])
from modlib import *
from bbload import *
import mods as M
i32 = M.i32
seqs = {}
bad = 0
for _, _, o in bb_records("BBProteinSequence").values():
    if isinstance(o, Exception): bad += 1; continue
    seqs[o.fields["accessionNumber"]] = o.fields["sequence"].replace("I", "L")
pi = records("ProteinIdentification", 103)
acc = lambda b: b[8:28].split(b'\x0a')[0].decode('latin-1')
c = collections.Counter(); miss = []; noseq = collections.Counter()
for rid, b in M.pep.items():
    s = b[4:34].split(b'\x0a')[0].decode()
    if i32(b, 35): s += M.tails[i32(b, 35)]
    p = i32(b, 122)
    if p not in pi: c["no protein record"] += 1; continue
    a = acc(pi[p])
    if a not in seqs: c["no sequence"] += 1; noseq[a[-6:] if "DECOY" in a.upper() or a.startswith(("REV", "rev", "##", "DECOY")) else "target"] += 1; continue
    if s.replace("I", "L") in seqs[a]: c["in sequence"] += 1
    else: c["NOT in sequence"] += 1; miss.append((rid, s, a)) if len(miss) < 4 else None
print(f"proteins {len(pi)} sequences {len(seqs)} (unparsed {bad}); PSM->protein {dict(c)}; no-seq kinds {dict(noseq)}")
decoys = collections.Counter(("decoy" if "DECOY" in acc(b).upper() or acc(b).startswith(("REV_", "rev_", "##")) else "target") for b in pi.values())
print("protein records", dict(decoys), "examples", [acc(b) for b in list(pi.values())[:3]])
print("misses", miss)
