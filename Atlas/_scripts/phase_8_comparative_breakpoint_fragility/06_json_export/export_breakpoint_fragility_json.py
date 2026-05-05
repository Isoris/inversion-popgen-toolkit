#!/usr/bin/env python3
# =============================================================================
# 06_json_export/export_breakpoint_fragility_json.py
# =============================================================================
"""
Final assembly step. Reads:

  - config/candidate_breakpoints.tsv
  - output/density/focal/<cid>__<focal>.json
  - output/density/comparative/<cid>__<sp>.json   (any number)
  - output/comparative/map/<cid>.homologs.json
  - output/classification/<cid>.classification.json

Writes:

  - output/per_candidate_json/<cid>.json
        the canonical atlas-loadable JSON described in README §4
  - output/per_chromosome_json/<chrom>.json
        chromosome-level rollup: list of per-candidate objects + chrom summary
  - output/summary_tables/breakpoint_fragility_summary.tsv
        one row per candidate; the columns specified in the original prompt

The exporter is the single place where the canonical schema is materialised.
Everything upstream writes intermediate JSONs; this step is the one to trust
when reading downstream.

Usage:
    python export_breakpoint_fragility_json.py \
        --candidates config/candidate_breakpoints.tsv \
        --species_manifest config/species_manifest.tsv \
        --focal_density_dir output/density/focal \
        --comparative_density_dir output/density/comparative \
        --homologs_dir output/comparative/map \
        --classification_dir output/classification \
        --out_root output
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    get_logger,
    normalize_chrom,
    read_tsv,
    write_run_manifest,
    write_tsv,
)


SCHEMA_VERSION = "comparative_breakpoint_fragility_v0.1"


# columns of the TSV summary, as specified in the original prompt
SUMMARY_COLS = [
    "candidate_id", "chrom", "start", "end", "left_bp", "right_bp",
    "gar_polymorphic",
    "left_TE_density_100kb", "right_TE_density_100kb",
    "left_fold_vs_chr", "right_fold_vs_chr",
    "left_percentile_chr", "right_percentile_chr",
    "both_breakpoints_TE_enriched",
    "inside_TE_density", "chrom_TE_density", "local_TE_density",
    "mac_boundary_present", "mac_orientation",
    "mac_left_TE_density_100kb", "mac_right_TE_density_100kb",
    "pangasius_boundary_present", "pangasius_orientation",
    "pangasius_left_TE_density_100kb", "pangasius_right_TE_density_100kb",
    "num_species_with_boundary", "num_species_TE_enriched",
    "architecture_class", "architecture_label",
    "age_model", "prediction_label", "confidence", "hotspot_tier",
    "safe_conclusion",
]


def _safe(d: dict | None, *path, default=None):
    cur = d
    for k in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
    return cur if cur is not None else default


def _first_target_for_species(per_target: list, species_token: str) -> dict | None:
    """species_token: substring matched case-insensitively against species_id."""
    tok = species_token.lower()
    for t in per_target:
        if tok in t["species_id"].lower():
            return t
    return None


def _get_100kb(side_block: dict | None) -> dict | None:
    if not side_block:
        return None
    return side_block.get("window_100kb")


def export_one(
    *,
    cand_row: dict,
    focal_density: dict,
    homologs: dict,
    classification: dict,
    focal_species_id: str,
    focal_assembly: str,
    log,
) -> tuple[dict, dict]:
    """Return (final_json, summary_row)."""
    cid   = cand_row["candidate_id"]
    chrom = normalize_chrom(cand_row["chrom"])

    L100 = _get_100kb(focal_density.get("left_breakpoint"))
    R100 = _get_100kb(focal_density.get("right_breakpoint"))

    left_te   = _safe(L100, "te_density")
    right_te  = _safe(R100, "te_density")
    left_fold = _safe(L100, "fold_vs_chromosome")
    right_fold= _safe(R100, "fold_vs_chromosome")
    left_pct  = _safe(L100, "percentile_chr")
    right_pct = _safe(R100, "percentile_chr")

    def _enriched(blk: dict | None) -> bool:
        if not blk:
            return False
        f = blk.get("fold_vs_chromosome")
        p = blk.get("percentile_chr")
        return (f is not None and f >= 2.0) or (p is not None and p >= 0.90)

    both_enriched = bool(_enriched(L100) and _enriched(R100))

    inside_te = _safe(focal_density, "inside_interval", "te_density")
    chrom_te  = _safe(focal_density, "chromosome_background", "te_density")
    local_te  = _safe(focal_density, "local_background_2mb", "te_density")

    per_target = classification.get("comparative_targets", [])
    n_with_boundary = sum(1 for t in per_target
                          if t["boundary_status"]["left_boundary_present"]
                          or t["boundary_status"]["right_boundary_present"])
    n_te_enriched   = sum(1 for t in per_target if t.get("te_enriched_at_homolog"))

    # convenience columns for Mac and Pangasius (project-specific)
    mac = _first_target_for_species(per_target, "macrocephalus")
    pan = _first_target_for_species(per_target, "pangasius")

    def _target_left_100kb(t: dict | None) -> dict | None:
        return _get_100kb(_safe(t, "repeat_density", "left_boundary"))

    def _target_right_100kb(t: dict | None) -> dict | None:
        return _get_100kb(_safe(t, "repeat_density", "right_boundary"))

    # --- final JSON --------------------------------------------------------
    out = {
        "schema_version": SCHEMA_VERSION,
        "candidate_id": cid,
        "focal_species": focal_species_id,
        "focal_assembly": focal_assembly,
        "focal_interval": {
            "chrom": chrom,
            "start": int(cand_row["start"]),
            "end":   int(cand_row["end"]),
            "left_breakpoint":  int(cand_row["left_breakpoint"]),
            "right_breakpoint": int(cand_row["right_breakpoint"]),
            "confidence_tier":  cand_row.get("confidence_tier") or None,
            "evidence_layers":  [s.strip() for s in
                                 (cand_row.get("evidence_layers") or "").split(",")
                                 if s.strip()],
            "arrangement_label": cand_row.get("arrangement_label") or None,
        },
        "gar_population_status": {
            "polymorphic_in_gar":
                str(cand_row.get("gar_polymorphic", "true")).lower() != "false",
            "evidence_layers":  [s.strip() for s in
                                 (cand_row.get("evidence_layers") or "").split(",")
                                 if s.strip()],
            "confidence":       cand_row.get("confidence_tier") or "unspecified",
        },
        "polymorphism_confirmed_in_species":
            classification.get("polymorphism_confirmed_in_species",
                               [focal_species_id]),
        "polymorphism_unknown_in_species":
            classification.get("polymorphism_unknown_in_species", []),

        "repeat_density": {
            "left_breakpoint":   focal_density.get("left_breakpoint", {}),
            "right_breakpoint":  focal_density.get("right_breakpoint", {}),
            "inside_interval":   focal_density.get("inside_interval"),
            "left_flank_outside":  focal_density.get("left_flank_outside"),
            "right_flank_outside": focal_density.get("right_flank_outside"),
            "chromosome_background": focal_density.get("chromosome_background", {}),
            "local_background_2mb":  focal_density.get("local_background_2mb", {}),
        },
        "comparative_species": [
            {
                "species_id":     t["species_id"],
                "assembly_name":  t.get("assembly_name") or "",
                "homologous_region": t["homologous_region"],
                "boundary_status":   t["boundary_status"],
                "repeat_density":    t["repeat_density"],
                "te_enriched_at_homolog": t.get("te_enriched_at_homolog", False),
                "interpretation":    t["interpretation"],
            }
            for t in per_target
        ],
        "classification": classification.get("classification", {}),
        "safe_conclusion": classification.get("safe_conclusion", ""),
    }

    # --- summary row -------------------------------------------------------
    cls = classification.get("classification", {})
    summary = {
        "candidate_id": cid,
        "chrom":        chrom,
        "start":        cand_row["start"],
        "end":          cand_row["end"],
        "left_bp":      cand_row["left_breakpoint"],
        "right_bp":     cand_row["right_breakpoint"],
        "gar_polymorphic": (
            "true" if str(cand_row.get("gar_polymorphic", "true")).lower() != "false"
            else "false"
        ),
        "left_TE_density_100kb":  left_te,
        "right_TE_density_100kb": right_te,
        "left_fold_vs_chr":       left_fold,
        "right_fold_vs_chr":      right_fold,
        "left_percentile_chr":    left_pct,
        "right_percentile_chr":   right_pct,
        "both_breakpoints_TE_enriched": "true" if both_enriched else "false",
        "inside_TE_density": inside_te,
        "chrom_TE_density":  chrom_te,
        "local_TE_density":  local_te,
        "mac_boundary_present": ("true" if mac and (
            mac["boundary_status"]["left_boundary_present"] or
            mac["boundary_status"]["right_boundary_present"]) else "false"),
        "mac_orientation": _safe(mac, "homologous_region",
                                 "orientation_relative_to_focal", default=""),
        "mac_left_TE_density_100kb":  _safe(_target_left_100kb(mac),  "te_density"),
        "mac_right_TE_density_100kb": _safe(_target_right_100kb(mac), "te_density"),
        "pangasius_boundary_present": ("true" if pan and (
            pan["boundary_status"]["left_boundary_present"] or
            pan["boundary_status"]["right_boundary_present"]) else "false"),
        "pangasius_orientation": _safe(pan, "homologous_region",
                                       "orientation_relative_to_focal", default=""),
        "pangasius_left_TE_density_100kb":  _safe(_target_left_100kb(pan),  "te_density"),
        "pangasius_right_TE_density_100kb": _safe(_target_right_100kb(pan), "te_density"),
        "num_species_with_boundary": n_with_boundary,
        "num_species_TE_enriched":   n_te_enriched,
        "architecture_class": cls.get("architecture_class"),
        "architecture_label": cls.get("architecture_label"),
        "age_model":          cls.get("age_model"),
        "prediction_label":   cls.get("prediction_label"),
        "confidence":         cls.get("confidence"),
        "hotspot_tier":       classification.get("hotspot_tier"),
        "safe_conclusion":    classification.get("safe_conclusion", ""),
    }

    return out, summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--species_manifest", required=True, type=Path)
    ap.add_argument("--focal_density_dir", required=True, type=Path)
    ap.add_argument("--comparative_density_dir", required=True, type=Path)
    ap.add_argument("--homologs_dir", required=True, type=Path)
    ap.add_argument("--classification_dir", required=True, type=Path)
    ap.add_argument("--out_root", required=True, type=Path)
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("06_export", args.log_dir)
    per_cand_dir = args.out_root / "per_candidate_json"
    per_chrom_dir = args.out_root / "per_chromosome_json"
    summary_dir = args.out_root / "summary_tables"
    per_cand_dir.mkdir(parents=True, exist_ok=True)
    per_chrom_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(args.candidates,
                     required_cols=["candidate_id", "chrom",
                                    "start", "end",
                                    "left_breakpoint", "right_breakpoint"])
    species = read_tsv(args.species_manifest,
                       required_cols=["species_id", "assembly_name", "is_focal"])
    focal_row = next(r for r in species
                     if str(r["is_focal"]).lower() == "true")
    focal_id  = focal_row["species_id"]
    focal_asm = focal_row.get("assembly_name") or ""

    summary_rows: list[dict] = []
    by_chrom: dict[str, list[dict]] = defaultdict(list)

    for c in cands:
        cid = c["candidate_id"]
        focal_density_path = args.focal_density_dir / f"{cid}__{focal_id}.json"
        homolog_path = args.homologs_dir / f"{cid}.homologs.json"
        cls_path = args.classification_dir / f"{cid}.classification.json"

        if not focal_density_path.exists():
            log.warning("missing focal density for %s; skipping", cid)
            continue
        if not cls_path.exists():
            log.warning("missing classification for %s; skipping", cid)
            continue

        focal_density  = json.loads(focal_density_path.read_text())
        homologs       = json.loads(homolog_path.read_text()) if homolog_path.exists() \
                         else {"targets": []}
        classification = json.loads(cls_path.read_text())

        final, summary = export_one(
            cand_row=c,
            focal_density=focal_density,
            homologs=homologs,
            classification=classification,
            focal_species_id=focal_id,
            focal_assembly=focal_asm,
            log=log,
        )

        (per_cand_dir / f"{cid}.json").write_text(
            json.dumps(final, indent=2, sort_keys=False)
        )
        summary_rows.append(summary)
        by_chrom[final["focal_interval"]["chrom"]].append(final)

    # per-chromosome rollups
    for chrom, items in by_chrom.items():
        items.sort(key=lambda x: x["focal_interval"]["start"])
        roll = {
            "schema_version": SCHEMA_VERSION,
            "chrom": chrom,
            "n_candidates": len(items),
            "summary": {
                "by_class":
                    {cls: sum(1 for it in items
                              if (it.get("classification") or {})
                                  .get("architecture_class") == cls)
                     for cls in "ABCDEFG"},
                "n_high_hotspot": sum(
                    1 for it in items
                    if (it.get("classification") or {}).get("confidence") == "high"
                ),
            },
            "candidates": items,
        }
        (per_chrom_dir / f"{chrom}.json").write_text(
            json.dumps(roll, indent=2, sort_keys=False)
        )

    # summary TSV
    summary_path = summary_dir / "breakpoint_fragility_summary.tsv"
    write_tsv(summary_path, summary_rows, SUMMARY_COLS)

    write_run_manifest(args.out_root / "_export.run.json",
                       n_candidates=len(summary_rows),
                       n_chromosomes=len(by_chrom))
    log.info("wrote %d per-candidate JSONs, %d per-chrom JSONs, summary %s",
             len(summary_rows), len(by_chrom), summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
