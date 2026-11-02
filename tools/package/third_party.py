#!/usr/bin/env python3
"""Write THIRD-PARTY-LICENSES for a package: every crate the shipped binaries
link, on this platform, with its license text and any NOTICE. Reads what cargo
already downloaded, so run it after a build. A crate shipping no license text
is listed with its SPDX expression; `--check` fails on a crate with neither."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

BINARIES = {"tuffold-app", "tuffold-cli", "tuffold-mcp"}
TEXT_PREFIXES = ("license", "licence", "copying", "notice", "unlicense", "copyright")


def metadata() -> dict:
    output = subprocess.run(
        ["cargo", "metadata", "--locked", "--format-version", "1", "--filter-platform", host()],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return json.loads(output)


def host() -> str:
    version = subprocess.run(["rustc", "-vV"], check=True, capture_output=True, text=True).stdout
    return next(line.split()[1] for line in version.splitlines() if line.startswith("host:"))


def linked(data: dict) -> list[dict]:
    """Packages reachable from the binaries over normal dependencies."""
    packages = {package["id"]: package for package in data["metadata_packages"]}
    nodes = {node["id"]: node for node in data["resolve"]["nodes"]}
    stack = [id for id, package in packages.items() if package["name"] in BINARIES]
    seen = set()
    while stack:
        at = stack.pop()
        if at in seen:
            continue
        seen.add(at)
        for dep in nodes[at]["deps"]:
            if any(kind["kind"] is None for kind in dep["dep_kinds"]):
                stack.append(dep["pkg"])
    return [packages[id] for id in seen]


def texts_in(folder: Path) -> list[Path]:
    return sorted(
        path
        for path in folder.iterdir()
        if path.is_file() and path.name.lower().startswith(TEXT_PREFIXES)
    )


def license_texts(package: dict) -> list[Path]:
    """The crate's own, or for a crate from a git workspace, its repository's."""
    root = Path(package["manifest_path"]).parent
    texts = texts_in(root)
    if texts or not (package.get("source") or "").startswith("git+"):
        return texts
    for folder in root.parents:
        texts = texts_in(folder)
        if texts or (folder / ".git").exists():
            return texts
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = metadata()
    data["metadata_packages"] = data.pop("packages")
    workspace = set(data["workspace_members"])
    packages = sorted(
        (p for p in linked(data) if p["id"] not in workspace),
        key=lambda p: (p["name"], p["version"]),
    )
    unlicensed = []
    sections = []
    for package in packages:
        name = f'{package["name"]} {package["version"]}'
        spdx = package.get("license") or package.get("license_file")
        if not spdx:
            unlicensed.append(name)
        lines = [f"{'=' * 78}\n{name}\nLicense: {spdx or 'not stated'}"]
        if package.get("repository"):
            lines.append(f'Source: {package["repository"]}')
        texts = license_texts(package)
        if not texts:
            lines.append("(the crate ships no license text; see the license named above)")
        for path in texts:
            body = path.read_text(encoding="utf-8", errors="replace").strip()
            lines.append(f"\n--- {path.name}\n{body}")
        sections.append("\n".join(lines))
    header = (
        "Third-party software in Tuffold's binaries, as built for "
        f"{host()}: {len(packages)} crates.\n\n"
    )
    args.output.write_text(header + "\n\n".join(sections) + "\n", encoding="utf-8")
    print(f"{args.output}: {len(packages)} crates")
    if unlicensed:
        print(f"no license stated: {', '.join(unlicensed)}", file=sys.stderr)
        if args.check:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
