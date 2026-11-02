import sys
sys.path.insert(0, sys.argv[1])
from jser import load, Obj, Arr, Block, short
roots = load(sys.argv[2]); want = sys.argv[3].split(",")
def summ(v, depth=0):
    if isinstance(v, Obj):
        if depth < 1 and v.cls.name.startswith(("java.lang", "java.util.Date","java.io.File")):
            return f"{short(v.cls.name)}{ {k: summ(x, depth+1) for k,x in v.fields.items()} }"[:160]
        return f"<{short(v.cls.name)}>"
    if isinstance(v, Arr):
        it = v.items
        return f"{v.cls.name}[{len(it)}]" + (f" {list(it[:6]) if not isinstance(it,(bytes,str)) else it[:40]!r}" if len(it) else "")
    if isinstance(v, str): return repr(v[:100])
    return repr(v)[:100]
seen=set(); found={}
stack=list(roots)
while stack:
    x=stack.pop()
    if id(x) in seen: continue
    seen.add(id(x))
    if isinstance(x, Obj):
        if short(x.cls.name) in want and short(x.cls.name) not in found: found[short(x.cls.name)]=x
        stack.extend(x.fields.values())
        for a in x.annotations.values(): stack.extend(a)
    elif isinstance(x, Arr) and isinstance(x.items, list): stack.extend(x.items)
    elif isinstance(x, list): stack.extend(x)
for name in want:
    o=found.get(name)
    print("=====", name)
    if not o: continue
    for k,v in o.fields.items(): print(f"   {k}: {summ(v)}")
    for c,a in o.annotations.items(): print(f"   [annotation {short(c)}]: {[summ(x) if not isinstance(x,Block) else 'block'+str(len(x)) for x in a][:8]}")
