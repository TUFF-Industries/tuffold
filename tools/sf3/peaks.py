"""Check an unpacked .sf3's stored peaks (read-only, standard library).

Run from inside the dump folder: python3 -I peaks.py TOOLS

A spectrum's peaks are its PeakList record's 20 (intensity, m/z) pairs plus
the PeakListMassIntensities chain named at offset 164 (232-byte blocks: the
chain's id, 28 pairs, the next block's slot in the same part; NaN after the
last pair). Reports, over every PeakList record: chains whose blocks carry
another chain's id, records whose 20 pairs reappear in their chain (never
seen), and peaks per spectrum.
"""
import math
import os
import struct
import sys

sys.path.insert(0, sys.argv[1])
from modlib import records, table  # noqa: E402

i32 = lambda b, o: struct.unpack_from("<i", b, o)[0]  # noqa: E731
EMPTY = {(0.0, 0.0), (-1.0, -1.0)}
per, shift, used = table("PeakListMassIntensities")
parts = {k: open(f"PeakListMassIntensities{k}", "rb").read()
         for k in range(len(used)) if os.path.exists(f"PeakListMassIntensities{k}")}


def chain(head):
    """Pairs of the chain starting at global id `head`, or None if broken."""
    part, slot, pairs = head >> shift, head & ((1 << shift) - 1), []
    for _ in range(100000):
        block = parts.get(part, b"")[slot * 232:(slot + 1) * 232]
        if len(block) < 232 or i32(block, 0) != head:
            return None
        values = struct.unpack_from("<56f", block, 4)
        for k in range(0, 56, 2):
            if math.isnan(values[k]) or math.isnan(values[k + 1]):
                return pairs
            pairs.append((values[k], values[k + 1]))
        slot = i32(block, 228)
        if slot == -1:
            return pairs
    return None


broken = shared = 0
sizes = []
for rid, block in records("PeakList", 180).items():
    picked = [p for p in (struct.unpack_from("<2f", block, 4 + 8 * k) for k in range(20)) if p not in EMPTY]
    head = i32(block, 164)
    rest = chain(head) if head else []
    if rest is None:
        broken += 1
        continue
    if {round(m, 4) for _, m in picked} & {round(m, 4) for _, m in rest}:
        shared += 1
    sizes.append(len(picked) + len(rest))
sizes.sort()
print(f"PeakList records {len(sizes) + broken}: broken chains {broken}, 20 pairs reappearing in their chain {shared}; "
      f"peaks per spectrum min {sizes[0]} median {sizes[len(sizes) // 2]} max {sizes[-1]}")
