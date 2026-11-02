import sys, os, struct
sys.path.insert(0, sys.argv[1])
from jser import load, Obj, Arr, Block, short
d = sys.argv[2]
def val(v, depth=0):
    if isinstance(v, Obj):
        f = {k: val(x, depth+1) for k, x in v.fields.items()} if depth < 2 else "…"
        ann = []
        for c, a in v.annotations.items():
            for x in a:
                if isinstance(x, Block):
                    n=len(x)//4
                    ann.append(list(struct.unpack(f">{n}i", x[:4*n]))[:12])
                else: ann.append(val(x, depth+1))
        return {"cls": short(v.cls.name), **({"f": f} if f else {}), **({"ann": ann} if ann else {})}
    if isinstance(v, Arr): return f"{v.cls.name}[{len(v.items)}]"
    return v
for name in sorted(os.listdir(d)):
    p = os.path.join(d, name)
    data = open(p,'rb').read()
    if not data.startswith(b"\xac\xed") or name == "Experiment": continue
    try:
        r = load(p)
    except Exception as e:
        print(name, "ERR", e); continue
    o = r[0]
    f = {k: val(x) for k, x in o.fields.items() if k not in ("cache",)}
    siblings = sorted(x for x in os.listdir(d) if x.startswith(name) and x != name)
    sizes = [(s, os.path.getsize(os.path.join(d, s))) for s in siblings]
    print(f"== {name} [{short(o.cls.name)}] parts={sizes[:4]}{'…' if len(sizes)>4 else ''}")
    print("   ", str(f)[:700])
