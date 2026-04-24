#!/usr/bin/env python3
"""
STEP_03A_cross_species_bridge.py — cross-species conservation (STEP 3A of 4)

=============================================================================
PIPELINE POSITION
=============================================================================
  STEP 1  assembled-junction forensics       (q4_mechanism)
  STEP 2  single-sided BND support           (q7_existence_audit)
→ STEP 3A cross-species synteny bridge       (cross_species)  ← THIS SCRIPT
  STEP 3B bp-pipeline bridge                 (bp_bridge)       (parallel to 3A)
  STEP 4  structural-class assignment        (phase_9)

Renamed from cross_species_bridge_v6.py in pass 18 (2026-04-24). Runs
in parallel with STEP 3B — both read independent inputs (external
annotation vs phase_6 TSVs) and both write evidence blocks for STEP 4.
Neither depends on the other.

Output block (synteny_v6.json) is consumed by STEP 4 primarily for
the `supported_shared_between_species_inversion` label decision, which
fires when q5_conservation_class ∈ {shared_ancestral, ASEAN_diagnostic}
and existence is confirmed/likely.

=============================================================================
UPGRADE OVER v5
=============================================================================
v5 treated Dollo as primary and tree-polarization equally.
v6 restructures:
  PRIMARY   — direct between-species SV-event overlap (Kuang 2026 style)
              + flank gene-family coherence per side (STEP_09c).
  SECONDARY — tree parsimony polarization (STEP_11), reported with
              confidence.
  CROSS-CHECK — Dollo (if computed), shown only for agreement/disagreement
                flag. Never as primary call.

Rationale: Dollo parsimony forbids re-gain after loss, which is exactly
the brittle assumption for real inversions (recurrence is common, and
sister species rearrangements are often complex rather than simple
binary characters).

=============================================================================
KEYS WRITTEN (v6)
=============================================================================
PRIMARY DIRECT ANNOTATION:
  q5_bs_event_overlap              bool
  q5_bs_event_type                 INVERSION | FUSION | FISSION | NONE | NA
  q5_bs_species_pair_observed      list of "A↔B" pairs where the event is seen
  q5_bs_event_symmetry             same_orientation | inverted_between_species | NA

FLANK COHERENCE (reuses v5):
  q3_left_flank_coherence_score    float
  q3_right_flank_coherence_score   float
  q3_left_flank_n_families         int
  q3_right_flank_n_families        int
  q3_flank_coherence_class         both_coherent | left_only | right_only |
                                    neither | untested

SECONDARY TREE POLARIZATION:
  q5_tree_polarization_attempted    bool
  q5_tree_polarization_direction    inverted | reference | ambiguous | NA
  q5_tree_polarization_confidence   high | medium | low | NA

DOLLO AS CROSS-CHECK (optional):
  q5_dollo_state                    existing key, still populated if Dollo run
  q5_dollo_vs_tree_concordance      agree | disagree | dollo_not_run |
                                     tree_not_run

SYNTHESIS:
  q5_conservation_class             species_specific_polymorphism |
                                    species_specific_rearranged |
                                    shared_ancestral |
                                    ASEAN_diagnostic |
                                    unclassified_complex |
                                    untested

=============================================================================
REGISTRY_CONTRACT
  BLOCKS_WRITTEN:
    - synteny_v6: registries/schemas/structured_block_schemas/synteny_v6.schema.json
      keys: q5_bs_event_overlap, q5_bs_event_type,
            q5_tree_polarization_direction, q5_tree_polarization_confidence,
            q5_dollo_vs_tree_concordance, q5_conservation_class,
            q3_left_flank_coherence_score, q3_right_flank_coherence_score,
            q3_left_flank_n_families, q3_right_flank_n_families,
            q3_flank_coherence_class
      status: WIRED
      note: Flank-coherence keys folded in from superseded
            flank_coherence.schema.json (archived 2026-04-24). Writer
            unchanged; schema's keys_extracted list extended.
  KEYS_IN: none
=============================================================================
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


def load_candidates(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_between_species(path):
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                r["q_start"] = int(r.get("q_start", 0))
                r["q_end"] = int(r.get("q_end", 0))
            except (ValueError, TypeError):
                continue
            rows.append(r)
    return rows


def load_flank_coherence(path):
    if not path or not os.path.exists(path):
        return {}
    out = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                bp = int(r["bp_pos"])
                score = float(r.get("coherence_score", "nan"))
                n_fam = int(r.get("n_families_matched", 0)) \
                    if r.get("n_families_matched") else 0
            except (KeyError, ValueError, TypeError):
                continue
            out[(r["chrom"], r.get("side", "unknown"))].append({
                "bp": bp, "score": score, "n_families": n_fam
            })
    return out


def load_polarization(path):
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                r["bp_pos"] = int(r.get("bp_pos", 0))
            except (ValueError, TypeError):
                continue
            rows.append(r)
    return rows


def load_dollo(path):
    """Optional — a per-candidate Dollo state if previously computed."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r.get("candidate_id")
            if cid:
                out[cid] = r.get("dollo_state", "unresolved")
    return out


# ---------------------------------------------------------------------------
# Primary: between-species event overlap
# ---------------------------------------------------------------------------
def overlap_events(cand_row, bs_rows, overlap_bp):
    """
    Find all between-species events overlapping this candidate's boundaries.
    Return list of {species_pair, event_type, side}.
    """
    try:
        chrom = cand_row["chrom"]
        left = int(cand_row["left_bp"])
        right = int(cand_row["right_bp"])
    except (KeyError, ValueError):
        return []

    hits = []
    for r in bs_rows:
        rchrom = r.get("query_chrom", "")
        if rchrom != chrom:
            continue
        mid = (r["q_start"] + r["q_end"]) // 2
        side = None
        if abs(mid - left) <= overlap_bp:
            side = "left"
        elif abs(mid - right) <= overlap_bp:
            side = "right"
        if side is None:
            continue
        hits.append({
            "species_pair": f"{r.get('query_species')}↔{r.get('target_species')}",
            "event_type": r.get("event_type", "unknown"),
            "side": side,
            "q_start": r["q_start"],
            "q_end": r["q_end"],
        })
    return hits


# ---------------------------------------------------------------------------
# Secondary: flank coherence
# ---------------------------------------------------------------------------
def flank_coherence(flank_idx, chrom, left_bp, right_bp, window=10000):
    def best_near(hits, target_bp):
        if not hits:
            return None, None
        nearest = min(hits, key=lambda h: abs(h["bp"] - target_bp))
        if abs(nearest["bp"] - target_bp) > window:
            return None, None
        return nearest["score"], nearest.get("n_families")

    left_s, left_n = best_near(flank_idx.get((chrom, "left"), []), left_bp)
    right_s, right_n = best_near(flank_idx.get((chrom, "right"), []), right_bp)

    thr = 0.5
    left_ok = left_s is not None and left_s > thr
    right_ok = right_s is not None and right_s > thr

    if left_s is None and right_s is None:
        label = "untested"
    elif left_ok and right_ok:
        label = "both_coherent"
    elif left_ok:
        label = "left_only"
    elif right_ok:
        label = "right_only"
    else:
        label = "neither"

    return left_s, right_s, left_n, right_n, label


def polarize_lookup(pol_rows, chrom, bp, window=10000):
    hits = [r for r in pol_rows
            if r.get("chrom") == chrom and abs(r["bp_pos"] - bp) <= window]
    if not hits:
        return None
    return min(hits, key=lambda r: abs(r["bp_pos"] - bp))


# ---------------------------------------------------------------------------
# Synthesis (v6 rules)
# ---------------------------------------------------------------------------
def derive_conservation_class(
    overlap_hits, flank_class, tree_direction, tree_conf
):
    """
    Six classes. v6 rules.
    """
    # Rule A: no cross-species inputs at all
    if not overlap_hits and flank_class == "untested" and not tree_direction:
        return "untested"

    # Rule B: species_specific_polymorphism — no rearrangement between species
    if not overlap_hits and flank_class in ("both_coherent", "left_only",
                                             "right_only"):
        return "species_specific_polymorphism"

    # Rule C: rearrangement between species but only in query — needs
    # at least one overlap AND coherent flanks in sister species
    if overlap_hits and flank_class == "both_coherent":
        species_pairs = {h["species_pair"] for h in overlap_hits}
        # if the event is observed in >1 species pair → shared_ancestral
        # e.g. Cgar↔Cmac AND Cgar↔Outgroup
        if len(species_pairs) >= 2:
            return "shared_ancestral"
        # Single pair — it's a two-species difference
        # If tree confidence is high and derived in one lineage → ASEAN_diagnostic
        if tree_direction in ("inverted", "reference") and tree_conf == "high":
            return "ASEAN_diagnostic"
        return "species_specific_rearranged"

    # Rule D: overlap but flanks disagree
    if overlap_hits and flank_class in ("left_only", "right_only", "neither"):
        return "unclassified_complex"

    # Default
    return "untested"


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
def try_registry(registries_root):
    if not registries_root:
        return None
    api_path = os.path.join(registries_root, "api", "python")
    if not os.path.isdir(api_path):
        return None
    if api_path not in sys.path:
        sys.path.insert(0, api_path)
    try:
        from registry_loader import load_registry
        return load_registry()
    except Exception as e:
        print(f"[cross_species_v6] registry unavailable: {e}", file=sys.stderr)
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0])
    p.add_argument("--candidates", required=True)
    p.add_argument("--between_species_bed", default="")
    p.add_argument("--flank_coherence_tsv", default="")
    p.add_argument("--polarized_tsv", default="")
    p.add_argument("--dollo_tsv", default="",
                   help="Optional cross-check: per-candidate dollo state")
    p.add_argument("--outdir", required=True)
    p.add_argument("--registries_root", default="")
    p.add_argument("--overlap_bp", type=int, default=100000)
    p.add_argument("--flank_window", type=int, default=10000)
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    cand_rows = load_candidates(args.candidates)
    print(f"[cross_species_v6] {len(cand_rows)} candidates loaded")

    bs_rows = load_between_species(args.between_species_bed)
    flank_idx = load_flank_coherence(args.flank_coherence_tsv)
    pol_rows = load_polarization(args.polarized_tsv)
    dollo = load_dollo(args.dollo_tsv)

    print(f"[cross_species_v6] between-species events: {len(bs_rows)}")
    print(f"[cross_species_v6] flank entries: "
          f"{sum(len(v) for v in flank_idx.values())}")
    print(f"[cross_species_v6] polarization calls: {len(pol_rows)}")
    print(f"[cross_species_v6] dollo calls: {len(dollo)}")

    reg = try_registry(args.registries_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = defaultdict(int)
    for row in cand_rows:
        try:
            chrom = row["chrom"]
            left = int(row["left_bp"])
            right = int(row["right_bp"])
        except (KeyError, ValueError):
            continue
        cid = row["candidate_id"]

        overlap_hits = overlap_events(row, bs_rows, args.overlap_bp)

        left_s, right_s, left_n, right_n, flank_cls = flank_coherence(
            flank_idx, chrom, left, right, args.flank_window
        )

        pol_left = polarize_lookup(pol_rows, chrom, left, args.flank_window)
        pol_right = polarize_lookup(pol_rows, chrom, right, args.flank_window)

        # Tree-polarization synthesis: pick left if available, else right
        pol_ref = pol_left or pol_right
        tree_direction = pol_ref.get("direction") if pol_ref else None
        tree_conf = pol_ref.get("confidence") if pol_ref else None

        # Dollo cross-check
        dollo_state = dollo.get(cid) if dollo else None
        if dollo_state and tree_direction:
            if tree_direction == "inverted" and dollo_state in ("inverted", "derived_inverted"):
                dollo_concord = "agree"
            elif tree_direction == "reference" and dollo_state == "reference":
                dollo_concord = "agree"
            else:
                dollo_concord = "disagree"
        elif dollo_state:
            dollo_concord = "tree_not_run"
        else:
            dollo_concord = "dollo_not_run"

        conservation = derive_conservation_class(
            overlap_hits, flank_cls, tree_direction, tree_conf
        )

        event_types = {h["event_type"] for h in overlap_hits}
        species_pairs = sorted({h["species_pair"] for h in overlap_hits})

        data_direct = {
            "q5_bs_event_overlap": bool(overlap_hits),
            "q5_bs_event_type": (
                "INVERSION" if "INVERSION" in event_types else
                "FUSION" if "FUSION" in event_types else
                "FISSION" if "FISSION" in event_types else
                "NONE" if not overlap_hits else "mixed"
            ),
            "q5_bs_species_pair_observed": species_pairs,
            "q5_bs_event_n_overlaps": len(overlap_hits),
            "q3_left_flank_coherence_score": left_s,
            "q3_right_flank_coherence_score": right_s,
            "q3_left_flank_n_families": left_n,
            "q3_right_flank_n_families": right_n,
            "q3_flank_coherence_class": flank_cls,
            "q5_tree_polarization_attempted": pol_ref is not None,
            "q5_tree_polarization_direction": tree_direction,
            "q5_tree_polarization_confidence": tree_conf,
            "q5_dollo_state": dollo_state,
            "q5_dollo_vs_tree_concordance": dollo_concord,
            "q5_conservation_class": conservation,
        }

        summary[conservation] += 1

        cand_out = outdir / cid / "structured"
        cand_out.mkdir(parents=True, exist_ok=True)
        if not args.dry_run:
            (cand_out / "synteny_v6.json").write_text(json.dumps({
                "block_type": "synteny_v6",
                "candidate_id": cid,
                "source_script": "STEP_03A_cross_species_bridge.py",
                "data": data_direct,
            }, indent=2, default=str))
            if reg is not None:
                try:
                    reg.evidence.write_block(
                        candidate_id=cid, block_type="synteny_v6",
                        data=data_direct,
                        source_script="STEP_03A_cross_species_bridge.py",
                    )
                except Exception as e:
                    print(f"[cross_species_v6] {cid}: registry write "
                          f"failed: {e}", file=sys.stderr)

    print(f"[cross_species_v6] summary: {dict(summary)}")


if __name__ == "__main__":
    main()
