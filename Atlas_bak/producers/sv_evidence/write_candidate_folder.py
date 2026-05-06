"""
write_candidate_folder.py
-------------------------
Shared library used by STEP_SV_GT_AGG, STEP_SV_EVID_COMB, and any future
SV-evidence producer. Writes JSON layers into the per-candidate folder
that the SV evidence atlas page expects:

    data/<chrom>/candidates/<candidate_id>/
      ├── sv_genotype_counts.json
      ├── sv_evidence_combinations.json
      └── sv_support_by_sample.json    (future, step 6)

Why this lives here, not in each producer:
  - Folder layout is shared; encoding it once prevents drift between
    producers.
  - File-naming, atomic-write, and manifest-update are non-trivial; one
    well-tested code path is safer than three nearly-identical ones.
  - Manifest aggregation (which layers exist for each candidate) needs
    to coordinate writes anyway.

Atomic writes: each JSON is written to <path>.tmp first, then os.rename
to its final name. This prevents the atlas page from reading a
half-written JSON if the producer is interrupted.

Author:  Quentin Andres
Project: MS_Inversions_North_african_catfish (atlas SV evidence page)
Cohort:  226-sample pure C. gariepinus hatchery (NOT F1 hybrid; NOT C. mac)
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

# Mapping from format_version → canonical filename (in the per-cand folder).
LAYER_FILENAME = {
    "sv_genotype_counts_v1":         "sv_genotype_counts.json",
    "sv_evidence_combinations_v1":   "sv_evidence_combinations.json",
    "sv_support_by_sample_v1":       "sv_support_by_sample.json",  # future
}


def candidate_dir(out_root: str | Path, chrom: str, candidate_id: str) -> Path:
    """
    Resolve the per-candidate output directory:
      <out_root>/<chrom>/candidates/<candidate_id>/

    The chrom argument lets the producer optionally group folders by
    pseudochromosome (matches Quentin's mental model of LG-by-LG work).
    Pass chrom="" or None to flatten:
      <out_root>/candidates/<candidate_id>/

    Always creates the directory if it doesn't exist.
    """
    p = Path(out_root)
    if chrom:
        p = p / chrom / "candidates" / candidate_id
    else:
        p = p / "candidates" / candidate_id
    p.mkdir(parents=True, exist_ok=True)
    return p


def atomic_write_json(payload: dict[str, Any], path: str | Path,
                      indent: int | None = None) -> None:
    """
    Atomically write a JSON payload to path. Strategy:
      1. Validate payload by json.dumps before touching disk.
      2. Write to a sibling .tmp file.
      3. os.replace (atomic on POSIX) to the final name.

    indent=None gives a compact one-line file (smaller, faster to load
    in the browser). Pass indent=2 for human-debuggable output.
    """
    path = Path(path)
    # Validate first — fail before we touch the filesystem.
    text = json.dumps(payload, separators=(",", ":") if indent is None else (",", ": "),
                      indent=indent, ensure_ascii=False, sort_keys=False)
    # Write into the same dir so the rename is atomic (cross-fs renames
    # are not atomic on most platforms).
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
        f.write("\n")
    os.replace(tmp, path)


def write_candidate_layer(payload: dict[str, Any], out_root: str | Path,
                          chrom: str, candidate_id: str,
                          indent: int | None = None) -> Path:
    """
    Write a single layer payload into the per-candidate folder. Looks up
    the canonical filename from `payload["format_version"]`. Returns the
    path that was written.

    Raises ValueError if the format_version is unknown.
    """
    fv = payload.get("format_version")
    if fv not in LAYER_FILENAME:
        raise ValueError(f"unknown format_version: {fv!r}; "
                         f"expected one of {list(LAYER_FILENAME)}")
    # Embed candidate_id if it isn't there already (consistency check).
    if payload.get("candidate_id") and payload["candidate_id"] != candidate_id:
        raise ValueError(
            f"payload candidate_id {payload['candidate_id']!r} != "
            f"argument candidate_id {candidate_id!r}"
        )
    payload.setdefault("candidate_id", candidate_id)

    dest = candidate_dir(out_root, chrom, candidate_id) / LAYER_FILENAME[fv]
    atomic_write_json(payload, dest, indent=indent)
    return dest


def update_chrom_manifest(out_root: str | Path, chrom: str) -> Path:
    """
    Re-scan <out_root>/<chrom>/candidates/ and rewrite a manifest at
    <out_root>/<chrom>/manifest.json. The manifest tells the atlas (or a
    standalone client) which candidate folders exist and which layers
    each candidate has.

    Format:
      {
        "format_version": "sv_evidence_chrom_manifest_v1",
        "chrom": "C_gar_LG28",
        "candidates": [
          {
            "candidate_id": "INV_LG28_001",
            "layers": ["sv_genotype_counts_v1", "sv_evidence_combinations_v1"]
          },
          ...
        ]
      }

    Idempotent — safe to re-run after each producer write.
    """
    p = Path(out_root) / chrom if chrom else Path(out_root)
    cand_root = p / "candidates"
    if not cand_root.is_dir():
        return p / "manifest.json"

    # Build reverse map filename → format_version
    file_to_fv = {v: k for k, v in LAYER_FILENAME.items()}

    cands = []
    for cand_path in sorted(cand_root.iterdir()):
        if not cand_path.is_dir():
            continue
        layers = []
        for fname, fv in file_to_fv.items():
            if (cand_path / fname).is_file():
                layers.append(fv)
        cands.append({"candidate_id": cand_path.name, "layers": layers})

    manifest = {
        "format_version": "sv_evidence_chrom_manifest_v1",
        "chrom": chrom or "",
        "candidates": cands,
    }
    out = p / "manifest.json"
    atomic_write_json(manifest, out, indent=2)
    return out


# Convenience used by the producer SLURM scripts: log the write to stderr
# in a way that's grep-able from the SLURM job log.
def log_write(path: Path) -> None:
    print(f"[write_candidate_folder] wrote {path}", file=sys.stderr)
