import sys, collections
sys.path.insert(0, sys.argv[1])
from jser import load, Obj, Arr, short
roots = load(sys.argv[2])
hist = collections.Counter()
seen = set()
def walk(x):
    stack=[x]
    while stack:
        x=stack.pop()
        if id(x) in seen: continue
        seen.add(id(x))
        if isinstance(x, Obj):
            hist[x.cls.name]+=1
            stack.extend(x.fields.values())
            for a in x.annotations.values(): stack.extend(a)
        elif isinstance(x, Arr) and isinstance(x.items, list):
            stack.extend(x.items)
        elif isinstance(x, list): stack.extend(x)
for r in roots: walk(r)
print("roots:", roots[:3])
for name, n in hist.most_common(70): print(f"{n:7} {name}")
