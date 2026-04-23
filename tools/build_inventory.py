#!/usr/bin/env python3
"""
build_inventory.py — walk a list of roots and emit INVENTORY.md

Collects every script (.sh, .R, .py, .c, and unsuffixed executables) and
extracts a one-line description from its header comment. Flags files with
no description as ??? so you can see them at a glance.

Usage:
  python3 build_inventory.py <root1> <root2> ... -o INVENTORY.md
  python3 build_inventory.py --config roots.txt -o INVENTORY.md

Example:
  python3 build_inventory.py \
    /scratch/.../inversion-popgen-toolkit \
    /scratch/.../inversion_modules \
    /scratch/.../het_roh \
    -o ~/INVENTORY.md
"""

import argparse
import os
import re
import stat
import sys
from pathlib import Path
from collections import defaultdict

# File types we consider "scripts"
SCRIPT_EXT = {".sh", ".R", ".r", ".py", ".c", ".cpp", ".pl", ".awk"}
# Directories to skip entirely
SKIP_DIRS = {".git", "node_modules", "__pycache__", ".Rproj.user",
             ".ipynb_checkpoints", "renv", ".snakemake",
             # Common output dirs we don't want to walk into
             "results", "figures", "logs", "tracks", "popstats_groups",
             "local_Q_multiscale", "tmp", "cache"}
# Files to skip
SKIP_FILES = {".DS_Store", "README"}


def is_executable_binary(path: Path) -> bool:
    """File with no known extension but executable by user."""
    if path.suffix.lower() in SCRIPT_EXT or path.suffix == "":
        try:
            mode = path.stat().st_mode
            return (mode & stat.S_IXUSR) != 0 and not path.is_symlink()
        except OSError:
            return False
    return False


def extract_description(path: Path, max_lines: int = 80) -> tuple[str, str]:
    """
    Return (description, detection_method).
    Detection methods: HEADER_BLOCK / FIRST_COMMENT / DOCSTRING / NONE.
    """
    try:
        with path.open("r", errors="replace") as f:
            lines = [next(f, "") for _ in range(max_lines)]
    except Exception as e:
        return (f"?? unreadable: {e}", "ERROR")

    # Skip shebang
    start = 1 if lines and lines[0].startswith("#!") else 0

    # ---------------------------------------------------------------------
    # Pattern 1: `# === ... ===` header block with a title line right below
    # Looks for two `====` or `----` lines with a description between them.
    # ---------------------------------------------------------------------
    rule_re = re.compile(r"^\s*#+\s*={5,}\s*$")
    for i in range(start, min(len(lines) - 2, max_lines - 2)):
        if rule_re.match(lines[i]):
            # The title line is usually the filename; the one after is the desc
            # Skip any filename lines; take the first non-empty comment that isn't a rule
            for j in range(i + 1, min(len(lines), i + 10)):
                if rule_re.match(lines[j]):
                    break
                raw = lines[j].strip()
                if not raw.startswith("#"):
                    continue
                content = raw.lstrip("#").strip()
                if not content:
                    continue
                # Skip if it looks like the file's own name
                if content.rstrip(".sh .R .py .c") == path.name.rstrip(".sh .R .py .c"):
                    continue
                if "===" in content or "---" in content:
                    continue
                if content.lower().startswith(("usage:", "inputs:", "outputs:",
                                               "example:", "dependencies:",
                                               "input:", "output:")):
                    continue
                return (content, "HEADER_BLOCK")

    # ---------------------------------------------------------------------
    # Pattern 2: first non-shebang comment line with real content
    # ---------------------------------------------------------------------
    comment_prefix = "#"
    for ln in lines[start:]:
        raw = ln.strip()
        if not raw:
            continue
        if raw.startswith(comment_prefix):
            content = raw.lstrip("#").strip()
            if not content or rule_re.match(ln):
                continue
            if len(content) < 4:
                continue
            return (content, "FIRST_COMMENT")
        # Non-comment line — stop looking
        break

    # ---------------------------------------------------------------------
    # Pattern 3 (Python): module docstring
    # ---------------------------------------------------------------------
    if path.suffix == ".py":
        joined = "".join(lines[start:])
        m = re.search(r'^\s*"""(.*?)"""', joined, re.S | re.M)
        if m:
            first_line = m.group(1).strip().split("\n")[0].strip()
            if first_line:
                return (first_line, "DOCSTRING")

    return ("???", "NONE")


def count_lines(path: Path) -> int:
    try:
        with path.open("rb") as f:
            return sum(1 for _ in f)
    except Exception:
        return 0


def collect(roots: list[Path]) -> list[dict]:
    items = []
    for root in roots:
        root = root.expanduser().resolve()
        if not root.exists():
            print(f"WARN: root not found: {root}", file=sys.stderr)
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            # Prune
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS
                           and not d.startswith(".") or d == ".github"]
            for fn in filenames:
                if fn in SKIP_FILES or fn.endswith((".bak", ".swp", ".log", ".pyc")):
                    continue
                p = Path(dirpath) / fn
                # Skip output files that happen to be under a root
                if p.suffix.lower() in SCRIPT_EXT or is_executable_binary(p):
                    try:
                        size = p.stat().st_size
                    except OSError:
                        continue
                    if size == 0 or size > 5_000_000:  # skip empty or huge
                        continue
                    desc, method = extract_description(p)
                    items.append({
                        "root": str(root),
                        "path": str(p),
                        "relpath": str(p.relative_to(root)),
                        "dir": str(p.parent.relative_to(root)),
                        "name": p.name,
                        "ext": p.suffix,
                        "lines": count_lines(p),
                        "size_kb": round(size / 1024, 1),
                        "desc": desc,
                        "method": method,
                    })
    return items


def render_markdown(items: list[dict], roots: list[Path]) -> str:
    out = []
    out.append("# Pipeline Inventory")
    out.append("")
    out.append("*Auto-generated by `build_inventory.py`. "
               "Shows every script in the searched roots, grouped by directory, "
               "with an auto-extracted one-line description.*")
    out.append("")
    out.append(f"**Searched roots** ({len(roots)}):")
    for r in roots:
        out.append(f"- `{r}`")
    out.append("")

    # Summary stats
    total = len(items)
    with_desc = sum(1 for it in items if it["desc"] != "???")
    no_desc = total - with_desc
    out.append(f"**Totals:** {total} scripts  "
               f"· **{with_desc}** with description  "
               f"· **{no_desc}** need attention (marked ???)")
    out.append("")

    # Things most worth flagging up top
    flagged = [it for it in items if it["desc"] == "???"]
    if flagged:
        out.append("## ??? Scripts without extractable description")
        out.append("")
        out.append("These have no header comment the walker could parse. "
                   "Review and either add a header or mark as deprecated.")
        out.append("")
        out.append("| script | lines | path |")
        out.append("|---|---|---|")
        for it in sorted(flagged, key=lambda x: (x["root"], x["path"])):
            out.append(f"| `{it['name']}` | {it['lines']} | `{it['path']}` |")
        out.append("")

    # Group by root → dir
    by_root = defaultdict(lambda: defaultdict(list))
    for it in items:
        by_root[it["root"]][it["dir"]].append(it)

    for root in sorted(by_root):
        root_short = os.path.basename(root) or root
        out.append(f"## {root_short}")
        out.append("")
        out.append(f"*Full path:* `{root}`  ·  "
                   f"*Scripts:* {sum(len(v) for v in by_root[root].values())}")
        out.append("")
        dirs_sorted = sorted(by_root[root].keys(),
                             key=lambda d: (d == "." and "" or d, d))
        for d in dirs_sorted:
            items_here = sorted(by_root[root][d], key=lambda x: x["name"])
            label = d if d != "." else "(root)"
            out.append(f"### `{label}/`  ({len(items_here)} scripts)")
            out.append("")
            out.append("| file | ln | description |")
            out.append("|---|---:|---|")
            for it in items_here:
                desc = it["desc"].replace("|", "\\|").replace("\n", " ")
                # Truncate for readability
                if len(desc) > 140:
                    desc = desc[:137] + "…"
                out.append(f"| `{it['name']}` | {it['lines']} | {desc} |")
            out.append("")

    # Directory index at the end for quick navigation
    out.append("## Quick directory index")
    out.append("")
    all_dirs = sorted({(it["root"], it["dir"]) for it in items})
    for root, d in all_dirs:
        count = sum(1 for it in items if it["root"] == root and it["dir"] == d)
        label = f"{os.path.basename(root)}/{d}" if d != "." else os.path.basename(root)
        out.append(f"- `{label}` — {count} scripts")
    out.append("")

    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(
        description="Build a pipeline inventory markdown from a list of roots.")
    ap.add_argument("roots", nargs="*", help="Root paths to walk")
    ap.add_argument("--config", help="File with one root per line (blanks / #-comments ignored)")
    ap.add_argument("-o", "--output", default="INVENTORY.md",
                    help="Output markdown path (default: INVENTORY.md)")
    ap.add_argument("--tsv", help="Also emit a raw TSV with all rows")
    args = ap.parse_args()

    roots = [Path(r) for r in args.roots]
    if args.config:
        with open(args.config) as f:
            for ln in f:
                ln = ln.strip()
                if ln and not ln.startswith("#"):
                    roots.append(Path(ln))

    if not roots:
        ap.error("at least one root required (positional or --config)")

    print(f"[build_inventory] walking {len(roots)} roots", file=sys.stderr)
    items = collect(roots)
    print(f"[build_inventory] collected {len(items)} scripts", file=sys.stderr)

    md = render_markdown(items, roots)
    Path(args.output).write_text(md)
    print(f"[build_inventory] wrote {args.output}", file=sys.stderr)

    if args.tsv:
        import csv
        with open(args.tsv, "w", newline="") as f:
            fieldnames = ["root", "dir", "name", "ext", "lines",
                          "size_kb", "method", "desc", "path"]
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for it in items:
                w.writerow({k: it[k] for k in fieldnames})
        print(f"[build_inventory] wrote {args.tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
