import sys, struct
from jser import load, Obj
INT_MIN = -2**31

def mod_types(path="ModificationTypes"):
    out = {}
    for r in load(path):
        stack = [r]; seen = set()
        while stack:
            x = stack.pop()
            if id(x) in seen: continue
            seen.add(id(x))
            if isinstance(x, Obj):
                if x.cls.name.endswith("ModificationType") and "monoMass" in x.fields:
                    f = x.fields
                    out[f["id"]] = (f["referenceName"], f["monoMass"], "".join(f["aminoAcidsModified"].items) if f["aminoAcidsModified"] else "", f["location"])
                stack.extend(x.fields.values())
                for a in x.annotations.values(): stack.extend(a)
    return out

def int_records(path, size=256):
    d = open(path, 'rb').read()
    recs = {}
    for b in range(len(d) // size):
        ints = struct.unpack(f'<{size//4}i', d[b*size:(b+1)*size])
        if b == 0 and not any(ints): continue
        vals = []
        for v in ints[1:]:
            if v == INT_MIN: break
            vals.append(v)
        recs[b] = (ints[0], vals, ints[-1])
    return recs

def blocks(path, size):
    d = open(path, 'rb').read()
    return [d[i*size:(i+1)*size] for i in range(len(d)//size)]

def seq(blk):
    return blk[4:35].split(b'\x0a')[0].decode('ascii', 'replace')

import os
def table(name):
    """ENTRIES_PER_FILE, FILE_SHIFT, blocks used per part, from a record file's descriptor."""
    stack = list(load(name)); seen = set()
    while stack:
        x = stack.pop()
        if id(x) in seen: continue
        seen.add(id(x))
        if isinstance(x, Obj):
            if x.cls.name.endswith("MathmaticTranslationTable"):
                f = x.fields
                ann = [a for v in f["files"].annotations.values() for a in v]
                raw = b"".join(bytes(a) for a in ann)
                ints = list(struct.unpack(f">{len(raw)//4}i", raw))
                return f["ENTRIES_PER_FILE"], f["FILE_SHIFT"], ints[1:1 + ints[0]]
            stack.extend(x.fields.values())
            for a in x.annotations.values(): stack.extend(a)
    raise ValueError(name)

def records(name, size):
    """id -> block bytes over every part, checking each block's own id."""
    per, shift, used = table(name)
    out = {}
    for part, n in enumerate(used):
        path = f"{name}{part}"
        if not os.path.exists(path): continue
        d = open(path, 'rb').read()
        for blk in range(1 if part == 0 else 0, min(n + 1, len(d) // size)):
            b = d[blk*size:(blk+1)*size]
            rid = (part << shift) | blk
            out[rid] = b
    return out
