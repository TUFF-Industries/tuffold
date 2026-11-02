import sys
sys.path.insert(0, sys.argv[1])
from jser import load, Obj, Arr, Block, ClassDesc, short
maxd = int(sys.argv[3]) if len(sys.argv) > 3 else 6
seen = set()
def show(v, d, name=""):
    pad = "  " * d
    if isinstance(v, Obj):
        if id(v) in seen:
            print(f"{pad}{name}<ref {short(v.cls.name)}>"); return
        seen.add(id(v))
        print(f"{pad}{name}{short(v.cls.name)}")
        if d >= maxd: return
        for k, x in v.fields.items(): show(x, d + 1, k + "=")
        for k, a in v.annotations.items():
            for x in a: show(x, d + 1, f"@{short(k)}:")
    elif isinstance(v, Arr):
        it = v.items
        if isinstance(it, (bytes, str)) or (it and not isinstance(it[0], (Obj, Arr))):
            print(f"{pad}{name}{v.cls.name}[{len(it)}] {it[:12]!r}")
        else:
            print(f"{pad}{name}{v.cls.name}[{len(it)}]")
            if d < maxd:
                for x in it[:40]: show(x, d + 1)
    elif isinstance(v, Block):
        print(f"{pad}{name}block[{len(v)}] {bytes(v[:48]).hex()}")
    elif isinstance(v, ClassDesc):
        print(f"{pad}{name}classdesc {v.name}")
    else:
        print(f"{pad}{name}{v!r}"[:200])
for r in load(sys.argv[2]): show(r, 0)
