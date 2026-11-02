"""Validate projects several ways and tabulate what passes (copies only; the
projects given are never written).

Each project is copied into a scratch folder, validated once per variant
with the given `tuffold` binary, and the copy kept as `<name>-<variant>.tuff`
so `entrapment.py` and `rt_consistency.py` can read it afterwards. A
variant is a label and the extra `tuffold validate` arguments, `plain` and
`deep` by default (`--no-deep` and none). Notes Tuffold printed are listed
by code under each row.

Usage (Python's standard library only):
    python3 -I tools/crosscheck/sweep.py --tuffold target/release/tuffold \
        --out SCRATCH PROJECT.tuff... [--variant plain=--no-deep --variant deep=]

Deep variants need the models (`tuffold models`, or TUFFOLD_MODELS).
"""

import argparse
import re
import shutil
import subprocess
import time
from pathlib import Path

LEVELS = ("spectra", "peptides", "protein groups")


def summary(output):
    counts = {}
    for level in LEVELS:
        match = re.search(rf"^  {level}: (\d+) targets pass \((\d+) decoys\)", output, re.M)
        counts[level] = f"{match[1]} ({match[2]})" if match else "—"
    notes = sorted(set(re.findall(r"^note \[([a-z-]+)\]", output, re.M)))
    return counts, notes


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("projects", nargs="+", type=Path)
    parser.add_argument("--tuffold", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--variant", action="append", metavar="LABEL=ARGS")
    args = parser.parse_args()
    variants = [v.split("=", 1) for v in (args.variant or ["plain=--no-deep", "deep="])]
    args.out.mkdir(parents=True, exist_ok=True)
    print("project\tvariant\tspectra\tpeptides\tprotein groups\tseconds\tnotes")
    for project in args.projects:
        for label, extra in variants:
            copy = args.out / f"{project.stem}-{label}.tuff"
            shutil.copyfile(project, copy)
            started = time.monotonic()
            run = subprocess.run(
                [args.tuffold, "validate", *extra.split(), str(copy)],
                capture_output=True,
                text=True,
            )
            seconds = time.monotonic() - started
            if run.returncode != 0:
                print(f"{project.stem}\t{label}\tfailed: {run.stderr.strip()[:200]}")
                continue
            counts, notes = summary(run.stdout + run.stderr)
            row = "\t".join(counts[level] for level in LEVELS)
            print(f"{project.stem}\t{label}\t{row}\t{seconds:.0f}\t{', '.join(notes)}", flush=True)


if __name__ == "__main__":
    main()
