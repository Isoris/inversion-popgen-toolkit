#!/usr/bin/env python3
# =============================================================================
# 00_manifest/build_te_manifest.py
# =============================================================================
"""
Walks a root directory full of messy TE / repeat annotation outputs from
multiple catfish species and emits a clean TSV manifest:

    species_id  assembly_name  file_path  file_format  priority  notes

Detection rules (priority lower = preferred):

  format               heuristic                                                 priority
  --------------       --------------------------------------------------------- --------
  edta_intact          filename ends with .EDTA.intact.gff3 (or .intact.gff3)    1
  edta_anno            filename ends with .EDTA.TEanno.gff3 (or .TEanno.gff3)    2
  repeatmasker_out     filename ends with .out and first lines look like RM      3
  gff3                 .gff3 / .gff with feature lines mentioning repeat/TE      4
  bed                  .bed / .bed.gz with looks-like-TE class column            5
  unknown              everything else (logged & skipped at normalize step)      99

Species inference:
  - first matched directory component against the species-map TSV
  - fallback: first looks-like-species token in the path (Gar / Mac / Pangasius /
    Silurus / Ictalurus / Clarias / etc.)
  - manual override via --species_map (TSV: path_substring \\t species_id \\t assembly_name)

Anything not classifiable is written with file_format=unknown so the row stays
visible to the user; nothing is silently dropped.

Usage:
    python build_te_manifest.py \
        --te_root /scratch/lt200308-agbsci/.../TE_outputs \
        --out config/te_file_manifest.tsv \
        [--species_map config/species_map.tsv] \
        [--log_dir output/logs]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Allow `from common import ...` from a sibling package directory
_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    get_logger,
    open_maybe_gz,
    read_tsv,
    write_run_manifest,
    write_tsv,
)


# -----------------------------------------------------------------------------
# format detection
# -----------------------------------------------------------------------------
RM_HEADER_RE = re.compile(r"^\s*SW\s+perc div\.|^\s*score\s+div", re.I)
GFF_REPEAT_RE = re.compile(r"\b(repeat_region|transposable_element|target_site_dup|"
                           r"long_terminal_repeat|TIR|MITE|helitron|LTR_retrotransposon|"
                           r"DNA_transposon|SINE|LINE)\b", re.I)


def sniff_format(path: Path, *, peek_bytes: int = 8192) -> str:
    """Return one of: edta_intact, edta_anno, repeatmasker_out, gff3, bed, unknown."""
    name = path.name.lower()

    # filename-based fast paths
    if name.endswith(".intact.gff3") or name.endswith(".intact.gff3.gz"):
        return "edta_intact"
    if name.endswith(".teanno.gff3") or name.endswith(".teanno.gff3.gz"):
        return "edta_anno"

    # peek
    try:
        with open_maybe_gz(path, "rt") as fh:
            head = fh.read(peek_bytes)
    except Exception:
        return "unknown"

    head_low = head.lower()

    if name.endswith(".out"):
        # RM .out: header carries "SW  perc div." or starts with "score div"
        if RM_HEADER_RE.search(head):
            return "repeatmasker_out"
        # some RM outputs have no header but have the col-19 layout; weak fallback
        if "left" in head_low and "matching" in head_low:
            return "repeatmasker_out"
        return "unknown"

    if name.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        if GFF_REPEAT_RE.search(head):
            return "gff3"
        # still call it gff3; normalize step will skip non-repeat features
        return "gff3"

    if name.endswith((".bed", ".bed.gz")):
        # accept; normalize step will detect a TE-class column or skip
        return "bed"

    return "unknown"


# -----------------------------------------------------------------------------
# species inference
# -----------------------------------------------------------------------------
DEFAULT_SPECIES_TOKENS = {
    "gariepinus": ("C_gariepinus", "C_gariepinus_assembly"),
    "c_gar":      ("C_gariepinus", "C_gariepinus_assembly"),
    "clarias_gar":("C_gariepinus", "C_gariepinus_assembly"),
    "macrocephalus": ("C_macrocephalus", "C_macrocephalus_assembly"),
    "c_mac":      ("C_macrocephalus", "C_macrocephalus_assembly"),
    "pangasius":  ("Pangasius_hypophthalmus", "Pangasius_assembly"),
    "silurus":    ("Silurus_glanis", "Silurus_assembly"),
    "ictalurus":  ("Ictalurus_punctatus", "Ictalurus_assembly"),
    "channel_cat":("Ictalurus_punctatus", "Ictalurus_assembly"),
    "clarias_macrocephalus": ("C_macrocephalus", "C_macrocephalus_assembly"),
}


def infer_species(path: Path, overrides: list[dict]) -> tuple[str, str]:
    """Return (species_id, assembly_name). Manual overrides take priority."""
    p_str = str(path)
    p_low = p_str.lower()
    # manual override first — first match wins
    for ov in overrides:
        sub = (ov.get("path_substring") or "").lower()
        if sub and sub in p_low:
            return (ov["species_id"], ov.get("assembly_name") or f"{ov['species_id']}_assembly")
    # token match
    for token, (sid, asm) in DEFAULT_SPECIES_TOKENS.items():
        if token in p_low:
            return sid, asm
    return ("UNKNOWN_SPECIES", "UNKNOWN_ASSEMBLY")


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--te_root", required=True, type=Path,
                    help="root directory holding messy TE outputs (recursively scanned)")
    ap.add_argument("--out", required=True, type=Path,
                    help="output manifest TSV (e.g. config/te_file_manifest.tsv)")
    ap.add_argument("--species_map", type=Path, default=None,
                    help="optional TSV with path_substring, species_id, assembly_name")
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    ap.add_argument("--include_unknown", action="store_true",
                    help="still write rows for files with file_format=unknown")
    args = ap.parse_args()

    log = get_logger("00_build_te_manifest", args.log_dir)
    log.info("scanning %s", args.te_root)

    overrides: list[dict] = []
    if args.species_map and args.species_map.exists():
        overrides = read_tsv(args.species_map,
                             required_cols=["path_substring", "species_id"])
        log.info("loaded %d species_map override rows", len(overrides))

    cands = [
        ".gff", ".gff3", ".gff.gz", ".gff3.gz",
        ".out", ".bed", ".bed.gz",
    ]

    rows: list[dict] = []
    n_seen = n_kept = n_unknown = 0
    for p in args.te_root.rglob("*"):
        if not p.is_file():
            continue
        if not any(p.name.lower().endswith(s) for s in cands):
            continue
        n_seen += 1

        fmt = sniff_format(p)
        if fmt == "unknown":
            n_unknown += 1
            log.info("unknown format: %s", p)
            if not args.include_unknown:
                continue

        sid, asm = infer_species(p, overrides)
        prio = {"edta_intact": 1, "edta_anno": 2, "repeatmasker_out": 3,
                "gff3": 4, "bed": 5, "unknown": 99}[fmt]

        rows.append({
            "species_id": sid,
            "assembly_name": asm,
            "file_path": str(p.resolve()),
            "file_format": fmt,
            "priority": prio,
            "notes": "" if sid != "UNKNOWN_SPECIES"
                     else "species could not be inferred from path; please edit",
        })
        n_kept += 1

    rows.sort(key=lambda r: (r["species_id"], int(r["priority"]), r["file_path"]))

    cols = ["species_id", "assembly_name", "file_path",
            "file_format", "priority", "notes"]
    write_tsv(args.out, rows, cols)
    write_run_manifest(args.out,
                       te_root=str(args.te_root.resolve()),
                       n_files_seen=n_seen,
                       n_files_kept=n_kept,
                       n_files_unknown=n_unknown,
                       species_seen=sorted({r["species_id"] for r in rows}))

    log.info("seen=%d kept=%d unknown=%d -> %s",
             n_seen, n_kept, n_unknown, args.out)
    log.info("species observed: %s",
             sorted({r["species_id"] for r in rows}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
