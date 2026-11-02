"""Embedded entrapment: rename every ENTRAP_X twin (and reverse_ENTRAP_X decoy) in a pepXML to its
real protein X, so twin peptides belong to present proteins. Standard entrapment twins are absent
proteins, so it cannot see errors that evidence pooling adds to present proteins (siblings, NSP);
embedded twins can. Alternatives that become duplicates within a hit are dropped. Score the result
with entrapment_seq.py (names no longer tell twins apart). Stdlib only.
usage: python3 -I embed_pepxml.py IN.pep.xml OUT.pep.xml"""
import re, sys
name = re.compile(r'protein="([^"]*)"')
seen = set(); dropped = renamed = 0
with open(sys.argv[1]) as src, open(sys.argv[2], "w") as out:
    for line in src:
        if "ENTRAP_" in line:
            line = line.replace("ENTRAP_", ""); renamed += 1
        s = line.lstrip()
        if s.startswith("<search_hit"):
            m = name.search(line); seen = {m.group(1)} if m else set()
        elif s.startswith("<alternative_protein"):
            m = name.search(line)
            if m:
                if m.group(1) in seen: dropped += 1; continue
                seen.add(m.group(1))
        out.write(line)
print(f"{sys.argv[1]}: {renamed} lines renamed, {dropped} duplicate alternatives dropped")
