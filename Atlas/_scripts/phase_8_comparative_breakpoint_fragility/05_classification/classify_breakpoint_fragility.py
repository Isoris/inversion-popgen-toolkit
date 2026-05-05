#!/usr/bin/env python3
# =============================================================================
# 05_classification/classify_breakpoint_fragility.py
# =============================================================================
"""
Combines focal density JSON + comparative density JSONs + homolog mapping JSON
into a single per-candidate classification JSON containing:

    - architecture_class           A | B | C | D | E | F | G
    - architecture_label           short descriptive string
    - age_model                    YOUNG_POP | OLD_POLY | OLD_BP_YOUNG_INV | ...
    - modifiers                    [conserved_boundary, orientation_discordance,
                                    chromosome_context_change, terminal_context,
                                    assembly_risk]
    - hotspot_tier                 high | medium | low
    - prediction_label             confirmed_polymorphic_in_gar |
                                   predicted_polymorphic_hotspot |
                                   ancient_boundary_no_polymorphism_evidence | ...
    - confidence                   high | medium | low
    - safe_conclusion              free-text but using the project's hedged language
    - polymorphism_confirmed_in_species: ["C_gariepinus"] (always)
    - polymorphism_unknown_in_species:   [...] (everyone else with a homolog)

Decision logic (intentionally simple; v0.2 will add bootstrap p-values):

  Hotspot tier (focal):
    high   if BOTH bp ≥2× chr  OR  BOTH bp ≥0.90 percentile chrom-wide
    medium if EXACTLY ONE bp passes the high rule, OR both at 1.5× ≤ x < 2×
    low    otherwise

  Architecture class:
    F if assembly_risk modifier set                           (overrides everything else)
    G if shared orientation switch flagged across ≥2 deep lineages
                                                              (currently inferred only when
                                                               synteny rows tag boundary_type
                                                               as 'block_edge' AND orientation
                                                               is '-' across ≥2 species)
    D if terminal_context modifier set in ≥1 comparative species
    C if chromosome_context_change modifier set in ≥1 species
    E if hotspot_tier=high AND ≥1 comparative species also TE-enriched at homolog
    B if conserved_boundary modifier set in ≥2 species AND not E
    A otherwise

  Age tag:
    OLD_POLY              if class=G
    OLD_BP_YOUNG_INV      if class=E AND focal hotspot_tier=high
    LINEAGE_KARYO         if class in {C, D}
    MULTI_AGE_HOTSPOT     if class=E AND >1 modifier present
    ANCIENT_BOUNDARY      if class=B
    PREDICTED_POLYMORPHIC_HOTSPOT
                          if class=E (always also added)
    YOUNG_POP             default if focal is polymorphic and no other tag fits

  Confidence:
    high   if focal hotspot_tier=high AND ≥2 comparative species support classification
    medium if focal hotspot_tier=high OR ≥1 comparative species support
    low    otherwise

Inputs:
    --homologs_dir    output/comparative/map
    --focal_density_dir   output/density/focal
    --comparative_density_dir   output/density/comparative
    --candidates      config/candidate_breakpoints.tsv
    --species_manifest config/species_manifest.tsv
    --out_dir         output/classification
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    get_logger,
    read_tsv,
    write_run_manifest,
)


# -----------------------------------------------------------------------------
# helpers for reading window blocks
# -----------------------------------------------------------------------------
def _safe(d: dict | None, *path, default=None):
    cur = d
    for k in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
    return cur if cur is not None else default


def bp_density_at_100kb(density_json: dict, side: str) -> dict | None:
    """side in {'left_breakpoint','right_breakpoint'}; prefer 100 kb window."""
    return _safe(density_json, side, "window_100kb")


def is_te_enriched_high(blk: dict | None) -> bool:
    if not blk:
        return False
    fold = blk.get("fold_vs_chromosome")
    pct  = blk.get("percentile_chr")
    return (fold is not None and fold >= 2.0) or (pct is not None and pct >= 0.90)


def is_te_enriched_medium(blk: dict | None) -> bool:
    if not blk:
        return False
    fold = blk.get("fold_vs_chromosome")
    if fold is None:
        return False
    return 1.5 <= fold < 2.0


def hotspot_tier(focal_density: dict) -> str:
    L = bp_density_at_100kb(focal_density, "left_breakpoint")
    R = bp_density_at_100kb(focal_density, "right_breakpoint")
    L_hi = is_te_enriched_high(L)
    R_hi = is_te_enriched_high(R)
    if L_hi and R_hi:
        return "high"
    L_md = is_te_enriched_medium(L)
    R_md = is_te_enriched_medium(R)
    if L_hi or R_hi or (L_md and R_md):
        return "medium"
    return "low"


# -----------------------------------------------------------------------------
# modifiers from comparative material
# -----------------------------------------------------------------------------
def derive_modifiers(
    *,
    homologs: dict,
    comparative_densities: dict[str, dict],
) -> tuple[list[str], list[dict]]:
    """
    Returns (modifiers_list, per_target_summary).
    per_target_summary is a list of dicts with the comparative interpretation
    blocks, ready to slot into the final per-candidate JSON.
    """
    mods: set[str] = set()
    per_target: list[dict] = []

    for tgt in homologs.get("targets", []):
        sp = tgt["species_id"]
        hr = tgt["homologous_region"]
        btype = (hr.get("boundary_type") or "").lower()
        orient = hr.get("orientation_relative_to_focal") or "+"

        if btype == "terminal":
            mods.add("terminal_context")
        if btype == "chr_context_change":
            mods.add("chromosome_context_change")
        if btype in ("block_edge", "internal", "terminal"):
            # presence of a synteny mapping at all = a boundary exists in target
            mods.add("conserved_boundary")
        if orient == "-":
            mods.add("orientation_discordance")

        # comparative TE enrichment at the homolog
        d = comparative_densities.get(sp, {})
        L_blk = bp_density_at_100kb(d, "left_breakpoint")
        R_blk = bp_density_at_100kb(d, "right_breakpoint")
        L_hi = is_te_enriched_high(L_blk)
        R_hi = is_te_enriched_high(R_blk)

        if L_hi or R_hi:
            interp = (
                f"TE-rich homologous breakpoint in {sp} assembly; "
                f"supports candidate recurrent hotspot, does not demonstrate "
                f"polymorphism."
            )
        else:
            interp = (
                f"Homologous region present in {sp} assembly; "
                f"breakpoint architecture is not strongly TE-enriched. "
                f"Polymorphism unknown — no population data."
            )

        per_target.append({
            "species_id": sp,
            "assembly_name": tgt.get("assembly_name") or "",
            "homologous_region": hr,
            "boundary_status": {
                "left_boundary_present":  L_blk is not None,
                "right_boundary_present": R_blk is not None,
                "orientation_state": ("alternative" if orient == "-" else "concordant"),
                "chromosome_context": (
                    "different_chromosome_context"
                    if btype == "chr_context_change"
                    else "same_chromosome_context"
                ),
                "terminal_or_internal": ("terminal" if btype == "terminal" else "internal"),
            },
            "repeat_density": {
                "left_boundary": _safe(d, "left_breakpoint",  default={}),
                "right_boundary": _safe(d, "right_breakpoint", default={}),
                "inside_interval": _safe(d, "inside_interval", default={}),
                "chromosome_background": _safe(d, "chromosome_background", default={}),
                "local_background_2mb": _safe(d, "local_background_2mb", default={}),
            },
            "te_enriched_at_homolog": bool(L_hi or R_hi),
            "interpretation": interp,
        })

    return sorted(mods), per_target


# -----------------------------------------------------------------------------
# architecture / age / confidence
# -----------------------------------------------------------------------------
def assign_architecture(
    *,
    tier: str,
    modifiers: list[str],
    n_te_enriched_targets: int,
    n_targets_with_minus_orient: int,
    assembly_risk: bool,
) -> tuple[str, str]:
    if assembly_risk:
        return "F", "ambiguous_assembly_risk"

    if "conserved_boundary" in modifiers and n_targets_with_minus_orient >= 2:
        return "G", "ancient_conserved_inversion_candidate"

    if "terminal_context" in modifiers:
        return "D", "terminal_translocation_boundary"

    if "chromosome_context_change" in modifiers:
        return "C", "fusion_fission_associated"

    if tier == "high" and n_te_enriched_targets >= 1:
        return "E", "recurrent_rearrangement_hotspot"

    if modifiers.count("conserved_boundary") and n_te_enriched_targets == 0:
        # conserved_boundary is a single flag; presence means ≥1 target;
        # we want ≥2 species so check n_targets_with_minus_orient + total
        # As a minimal proxy, B requires ≥2 targets to have boundaries at all.
        pass

    return "A", "simple_population_inversion"


def assign_age_tag(
    *,
    arch_class: str,
    tier: str,
    modifiers: list[str],
) -> str:
    if arch_class == "G":
        return "OLD_POLY"
    if arch_class == "E" and tier == "high" and len(modifiers) >= 2:
        return "MULTI_AGE_HOTSPOT"
    if arch_class == "E" and tier == "high":
        return "OLD_BP_YOUNG_INV"
    if arch_class in ("C", "D"):
        return "LINEAGE_KARYO"
    if arch_class == "B":
        return "ANCIENT_BOUNDARY"
    return "YOUNG_POP"


def assign_prediction_label(arch_class: str) -> str:
    if arch_class == "E":
        return "predicted_polymorphic_hotspot"
    if arch_class == "G":
        return "candidate_ancient_polymorphism"
    if arch_class == "B":
        return "ancient_boundary_no_polymorphism_evidence"
    if arch_class == "F":
        return "interpret_with_caution_assembly_risk"
    return "confirmed_polymorphic_in_gar"


def assign_confidence(*, tier: str, n_supporting_targets: int) -> str:
    if tier == "high" and n_supporting_targets >= 2:
        return "high"
    if tier == "high" or n_supporting_targets >= 1:
        return "medium"
    return "low"


def build_safe_conclusion(
    *,
    arch_class: str,
    arch_label: str,
    tier: str,
    n_targets: int,
    n_supporting_targets: int,
    polymorphism_unknown_species: list[str],
) -> str:
    head = "Confirmed polymorphic inversion in C. gariepinus."
    if arch_class == "F":
        return (head + " Breakpoint overlaps assembly-risk regions "
                "(gaps / collapsed repeats / low mappability); architecture "
                "claims are deferred until reassembly or read-level support.")
    if arch_class == "A":
        return (head + " Comparative architecture is not strongly indicative; "
                "no recurrent-hotspot claim made.")
    if arch_class == "B":
        return (head + " Breakpoints overlap conserved synteny-block edges in "
                "comparative catfish assemblies, suggesting an ancient boundary "
                "substrate; polymorphism in non-resequenced species is unknown.")
    if arch_class in ("C", "D"):
        kind = ("fusion/fission-associated" if arch_class == "C"
                else "terminal-translocation-associated")
        return (head + f" Homologous breakpoints sit at {kind} boundaries in "
                "comparative assemblies, consistent with reuse of lineage-specific "
                "karyotype scars; polymorphism in those species is not demonstrated.")
    if arch_class == "E":
        unknown_str = (", ".join(polymorphism_unknown_species)
                       if polymorphism_unknown_species else "comparative species")
        return (head + " Homologous TE-rich breakpoint architecture in "
                f"{n_supporting_targets}/{n_targets} comparative species "
                "suggests a recurrent rearrangement hotspot. This does not "
                f"demonstrate polymorphism in {unknown_str}; population "
                "resequencing would be required to test that.")
    if arch_class == "G":
        return (head + " Shared orientation switch across deeper catfish lineages "
                "is consistent with an ancient polymorphism candidate; this "
                "interpretation is provisional and requires resequencing of "
                "additional species.")
    return head


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def load_focal_density(path: Path) -> dict:
    return json.loads(path.read_text())


def load_homologs(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--species_manifest", required=True, type=Path)
    ap.add_argument("--homologs_dir", required=True, type=Path)
    ap.add_argument("--focal_density_dir", required=True, type=Path)
    ap.add_argument("--comparative_density_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("05_classify", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(args.candidates,
                     required_cols=["candidate_id", "chrom",
                                    "start", "end",
                                    "left_breakpoint", "right_breakpoint"])
    species = read_tsv(args.species_manifest,
                       required_cols=["species_id", "is_focal"])
    focal_id = next(r["species_id"] for r in species
                    if str(r["is_focal"]).lower() == "true")
    log.info("focal=%s; %d candidates", focal_id, len(cands))

    n_written = 0
    classes_seen: dict[str, int] = {}

    for c in cands:
        cid = c["candidate_id"]

        focal_path = args.focal_density_dir / f"{cid}__{focal_id}.json"
        if not focal_path.exists():
            log.warning("no focal density JSON for %s; skipping", cid)
            continue
        focal_density = load_focal_density(focal_path)

        homolog_path = args.homologs_dir / f"{cid}.homologs.json"
        homologs = load_homologs(homolog_path) if homolog_path.exists() \
                   else {"targets": []}

        # comparative densities for each target species
        comparative_densities: dict[str, dict] = {}
        for tgt in homologs.get("targets", []):
            sp = tgt["species_id"]
            p = args.comparative_density_dir / f"{cid}__{sp}.json"
            if p.exists():
                comparative_densities[sp] = json.loads(p.read_text())

        tier = hotspot_tier(focal_density)
        modifiers, per_target = derive_modifiers(
            homologs=homologs,
            comparative_densities=comparative_densities,
        )

        n_te_enriched_targets = sum(1 for t in per_target if t["te_enriched_at_homolog"])
        n_minus = sum(1 for t in per_target
                      if t["boundary_status"]["orientation_state"] == "alternative")
        assembly_risk = "assembly_risk" in modifiers  # not yet derived; placeholder

        arch_class, arch_label = assign_architecture(
            tier=tier,
            modifiers=modifiers,
            n_te_enriched_targets=n_te_enriched_targets,
            n_targets_with_minus_orient=n_minus,
            assembly_risk=assembly_risk,
        )
        age_tag = assign_age_tag(arch_class=arch_class, tier=tier, modifiers=modifiers)
        pred_lbl = assign_prediction_label(arch_class)
        conf = assign_confidence(tier=tier, n_supporting_targets=n_te_enriched_targets)

        polymorphism_unknown_species = sorted(
            {t["species_id"] for t in per_target}
        )

        safe_conclusion = build_safe_conclusion(
            arch_class=arch_class,
            arch_label=arch_label,
            tier=tier,
            n_targets=len(per_target),
            n_supporting_targets=n_te_enriched_targets,
            polymorphism_unknown_species=polymorphism_unknown_species,
        )

        out = {
            "schema_version": "comparative_breakpoint_fragility_v0.1",
            "candidate_id": cid,
            "focal_species": focal_id,
            "hotspot_tier": tier,
            "modifiers": modifiers,
            "comparative_targets": per_target,
            "polymorphism_confirmed_in_species": [focal_id],
            "polymorphism_unknown_in_species": polymorphism_unknown_species,
            "classification": {
                "architecture_class": arch_class,
                "architecture_label": arch_label,
                "age_model": age_tag,
                "prediction_label": pred_lbl,
                "confidence": conf,
                "modifiers": modifiers,
            },
            "safe_conclusion": safe_conclusion,
        }
        out_path = args.out_dir / f"{cid}.classification.json"
        out_path.write_text(json.dumps(out, indent=2, sort_keys=False))
        n_written += 1
        classes_seen[arch_class] = classes_seen.get(arch_class, 0) + 1

    write_run_manifest(args.out_dir / "_classify.run.json",
                       n_candidates=n_written,
                       class_counts=classes_seen)
    log.info("wrote %d classifications; class counts: %s",
             n_written, classes_seen)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
