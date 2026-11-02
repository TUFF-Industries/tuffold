"""Profile every byte offset of a fixed-block record file across records."""
import sys, struct, math
path, size = sys.argv[1], int(sys.argv[2])
skip = set(range(int(sys.argv[3]), int(sys.argv[4]))) if len(sys.argv) > 4 else set()
data = open(path, 'rb').read()
blocks = [data[i*size:(i+1)*size] for i in range(1, len(data)//size)]
# keep only head blocks whose id equals position (first-run records)
recs = [b for b in blocks if struct.unpack('<i', b[:4])[0] > 0]
print(f"{len(recs)} records of {size} bytes")
def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals: return None
    return min(vals), max(vals), len(set(vals))
for off in range(4, size):
    if off in skip: continue
    out = []
    for fmt, n, kind in (('<i',4,'i32'), ('<f',4,'f32'), ('<d',8,'f64'), ('<h',2,'i16')):
        if off + n > size: continue
        vals = [struct.unpack(fmt, r[off:off+n])[0] for r in recs]
        if kind.startswith('f'):
            vals = [v if math.isfinite(v) and abs(v) < 1e7 and (v == 0 or abs(v) > 1e-6) else None for v in vals]
            good = sum(v is not None for v in vals)
            if good < len(recs) * 0.9: continue
            s = stats(vals)
            if s and s[2] > 3: out.append(f"{kind}[{s[0]:.4g},{s[1]:.4g}]#{s[2]}")
        else:
            s = stats(vals)
            if s and -1e6 < s[0] and s[1] < 1e7 and s[2] > 1: out.append(f"{kind}[{s[0]},{s[1]}]#{s[2]}")
    b = [r[off] for r in recs]
    out.append(f"u8[{min(b)},{max(b)}]#{len(set(b))}")
    print(off, " ".join(out))
