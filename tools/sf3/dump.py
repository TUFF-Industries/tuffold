import sqlite3, gzip, sys, os
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
con = sqlite3.connect(f"file:{src}?mode=ro", uri=True)
tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")]
if "SplitFileIdentifiers" in tables:
    print("SplitFileIdentifiers schema:", con.execute("SELECT sql FROM sqlite_master WHERE name='SplitFileIdentifiers'").fetchone()[0])
    for r in con.execute("SELECT * FROM SplitFileIdentifiers LIMIT 5"): print("  ", r)
print("ProteinGroup sample:", con.execute("SELECT * FROM ProteinGroup LIMIT 3").fetchall(), con.execute("SELECT COUNT(*) FROM ProteinGroup").fetchone())
total = 0
for name, blob in con.execute("SELECT filename, file FROM Files"):
    data = gzip.decompress(blob) if blob and blob[:2] == b"\x1f\x8b" else (blob or b"")
    total += len(data)
    open(os.path.join(out, name), "wb").write(data)
print("decompressed total", total)
