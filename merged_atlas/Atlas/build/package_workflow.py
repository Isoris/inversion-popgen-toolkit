#!/usr/bin/env python3
"""
build/package_workflow.py
=========================

Release-time splitter. For daily work, the Atlas is one combined repo
with all workflows mixed at root. At submission / release time, this
script extracts ONE workflow into a self-contained drop-in pack.

Usage:
    python3 build/package_workflow.py inversion
    python3 build/package_workflow.py diversity --output dist/Atlas-diversity-2026-05.tar.gz

Output: a tarball that the recipient can untar into any folder. The
extracted contents drop in next to other workflow packs without
collision (filenames are namespaced by workflow prefix).

What goes in each pack:

    Always (the "core" infrastructure every pack needs):
        Atlas/README.md
        Atlas/shared/                       (with SHARED_VERSION baked in)
        Atlas/build/flatten.py              (so users can re-flatten if they edit)
        Atlas/data/README.md                (the layout doc)
        Atlas/data/precomp/.gitkeep
        Atlas/data/cohort/.gitkeep
        Atlas/data/candidates/.gitkeep
        Atlas/data/comparative/.gitkeep
        Atlas/data/review/<workflow>/.gitkeep    (writable folder for this workflow)

    Workflow-specific:
        Atlas/<workflow>_*.html             (all top-level HTMLs starting with workflow_)
        Atlas/<workflow>_*/                 (matching JS module folders)

    Cohort/precomp/candidates/comparative DATA files: the user supplies
    these from their local pipeline output. Empty pack ships with empty
    folders. (Including 17 MB LG28.json in every pack would defeat the
    point.)
"""

from __future__ import annotations
import argparse, os, re, sys, tarfile
from datetime import date
from pathlib import Path
from typing import List


# ---------------------------------------------------------------------
# What ships in every pack — the shared scaffolding.
# ---------------------------------------------------------------------
CORE_FILES = [
    'README.md',
    'shared/state_io.js',
    'build/flatten.py',
    'build/package_workflow.py',     # so the recipient can re-pack if they want
    'data/README.md',
    'tests/test_modular_smoke.js',
]

CORE_DIRS_KEEP = [
    'data/precomp',
    'data/cohort',
    'data/candidates',
    'data/comparative',
]

# Filename patterns that match workflow-specific files at repo root.
# `inversion_*.html` and `inversion_*/` for the inversion workflow.
def workflow_html_glob(workflow: str) -> re.Pattern:
    return re.compile(rf'^{re.escape(workflow)}_[a-z0-9_]+\.html$')

def workflow_dir_glob(workflow: str) -> re.Pattern:
    return re.compile(rf'^{re.escape(workflow)}_[a-z0-9_]+$')


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------

def collect_files(repo_root: Path, workflow: str) -> List[Path]:
    """Walk the repo and return all files that belong in the pack."""
    out: List[Path] = []

    # Core files — must exist
    for rel in CORE_FILES:
        p = repo_root / rel
        if not p.exists():
            raise SystemExit(f"package: missing core file: {rel}")
        out.append(p)

    # Workflow HTMLs at repo root
    html_re = workflow_html_glob(workflow)
    for p in sorted(repo_root.iterdir()):
        if p.is_file() and html_re.match(p.name):
            out.append(p)

    # Workflow JS module folders
    dir_re = workflow_dir_glob(workflow)
    for p in sorted(repo_root.iterdir()):
        if p.is_dir() and dir_re.match(p.name):
            for sub in p.rglob('*'):
                if sub.is_file():
                    out.append(sub)

    if not out:
        raise SystemExit(f"package: no files found for workflow '{workflow}'. "
                         f"Expected at least one '{workflow}_*.html' at repo root.")
    return out


def make_tarball(repo_root: Path, files: List[Path], workflow: str, out_path: Path) -> None:
    """Build the tarball with an Atlas/ prefix; recipients untar with --strip-components=0
    or rename the top dir as they like."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    arc_root = f"Atlas-{workflow}"
    seen_dirs = set()

    with tarfile.open(out_path, 'w:gz') as tar:
        for f in files:
            rel = f.relative_to(repo_root)
            arcname = f"{arc_root}/{rel.as_posix()}"
            tar.add(f, arcname=arcname)
            # Track parent dirs we've added files to
            for parent in rel.parents:
                seen_dirs.add(parent.as_posix())

        # Add empty .gitkeep entries for the read-only data folders so
        # the pack arrives with the canonical folder skeleton even when
        # the user has no data yet.
        for d in CORE_DIRS_KEEP:
            arcname = f"{arc_root}/{d}/.gitkeep"
            info = tarfile.TarInfo(name=arcname)
            info.size = 0
            tar.addfile(info)

        # Workflow-specific writable folder under data/review/<workflow>/
        rev_keep = f"{arc_root}/data/review/{workflow}/.gitkeep"
        info = tarfile.TarInfo(name=rev_keep)
        info.size = 0
        tar.addfile(info)


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description="Package one workflow into a self-contained tarball")
    ap.add_argument('workflow', help='workflow name, e.g. inversion / diversity / population / assembly')
    ap.add_argument('--output', help='output tarball path (default: dist/Atlas-<workflow>-YYYY-MM.tar.gz)')
    ap.add_argument('--root', default=None, help='repo root (default: parent of this script)')
    args = ap.parse_args(argv)

    if not re.match(r'^[a-z][a-z0-9_]*$', args.workflow):
        ap.error('workflow name must be lowercase letters/digits/underscore')

    if args.root:
        root = Path(args.root).resolve()
    else:
        # Default: parent of build/
        root = Path(__file__).resolve().parent.parent

    if args.output:
        out = Path(args.output).resolve()
    else:
        ym = date.today().strftime('%Y-%m')
        out = root / 'dist' / f'Atlas-{args.workflow}-{ym}.tar.gz'

    files = collect_files(root, args.workflow)
    make_tarball(root, files, args.workflow, out)

    # Report
    total_bytes = sum(f.stat().st_size for f in files)
    print(f"package: workflow '{args.workflow}'")
    print(f"         {len(files)} files, {total_bytes:,} bytes source")
    print(f"         tarball: {out} ({out.stat().st_size:,} bytes compressed)")
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
