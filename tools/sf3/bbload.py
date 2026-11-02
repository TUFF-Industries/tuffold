import os, struct
from jser import load
from bbser import BBReader, class_list, BBObj

def bb_records(name):
    """record id -> parsed object, for a block-based object file."""
    desc = load(name)[0]
    classes = class_list(desc)
    idm = desc.fields["serialization"].fields["idToLocation"].fields
    files, locs = idm["file"].items, idm["location"].items
    cache, out = {}, {}
    for rid, (f, loc) in enumerate(zip(files, locs)):
        if f < 0 or loc < 0: continue
        path = f"{name}{f}"
        if path not in cache:
            if not os.path.exists(path): continue
            cache[path] = open(path, "rb").read()
        d = cache[path]
        if d[loc] == 0x73:
            hid, nxt, start = None, None, loc
        else:
            (hid, nxt), start = struct.unpack("<ii", d[loc:loc + 8]), loc + 8
        r = BBReader(d, classes, start)
        try:
            out[rid] = (hid, nxt, r.obj())
        except Exception as e:
            out[rid] = (hid, nxt, e)
    return out

def plain(v, depth=0):
    if isinstance(v, BBObj):
        if depth > 4: return f"<{v.cls.name}>"
        return {"_": v.cls.name.rsplit(".", 1)[-1], **{k: plain(x, depth + 1) for k, x in v.fields.items()}, **({"custom": v.custom} if v.custom else {})}
    return v
