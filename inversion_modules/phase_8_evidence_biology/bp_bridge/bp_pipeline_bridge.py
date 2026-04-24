#!/usr/bin/env python3
"""
bp_pipeline_bridge.py — reads user's existing breakpoint_pipeline outputs
and writes derived keys to the registry.

=============================================================================
ROLE
=============================================================================
The user's breakpoint_pipeline (/followup/breakpoint_pipeline/) already
produces the primary breakpoint signal: per-carrier ancestral fragments,
KDE-mode breakpoints with bootstrap CIs, and a 7-method consensus merge.

This bridge reads the three output TSVs and writes derived Phase 4 keys
to the registry so the v6 structural-class assigner can consume them.

**Supersedes my earlier proposal** of a separate interior-structure
diagnostic, which was re-deriving a coarser version of the same signal
the user had already designed. This script is the wiring the user asked
for.

=============================================================================
INPUTS — three TSVs produced by breakpoint_pipeline
=============================================================================
  --consensus_tsv    candidate_breakpoints_consensus.tsv
                      One row per candidate. Contains final_left_bp,
                      final_right_bp, CI bounds, n_methods_agreeing,
                      sources_available, primary_source,
                      inversion_size_bp, refined_vs_input_shift.

  --fragments_tsv    candidate_ancestral_fragments.tsv
                      One row per (candidate, sample) for each carrier.
                      Contains frag_start_bp, frag_end_bp, frag_length_bp,
                      extended_left_snps, extended_right_snps,
                      left_cor_at_boundary, right_cor_at_boundary.

  --per_method_tsv   candidate_breakpoints_per_method.tsv
                      One row per (candidate, method, side). Contains
                      method name, side, breakpoint_bp, weight, CI,
                      distance_to_consensus_kb, within_agreement_band.

=============================================================================
DERIVED KEYS WRITTEN (per candidate)
=============================================================================
BREAKPOINT CONSENSUS (from consensus_tsv):
  q3_final_left_bp                     int
  q3_final_right_bp                    int
  q3_left_ci_low_bp                    int
  q3_left_ci_high_bp                   int
  q3_right_ci_low_bp                   int
  q3_right_ci_high_bp                  int
  q3_left_ci_width_kb                  float
  q3_right_ci_width_kb                 float
  q3_n_methods_agreeing_left           int
  q3_n_methods_agreeing_right          int
  q3_breakpoint_precision_class        snp_resolution | sub_5kb |
                                        sub_30kb | coarse
  q3_primary_source                    method name
  q3_sources_available                 comma-separated list

FRAGMENT DISTRIBUTION (from fragments_tsv — the interior signal):
  q2_n_fragment_carriers               int
  q2_fragment_interior_recomb_fraction float in [0,1]
  q2_fragment_left_distribution_shape  unimodal | right_skewed | bimodal |
                                        broad | degenerate
  q2_fragment_right_distribution_shape analogous
  q2_interior_class                    clean_simple |
                                        edge_recombinants_only |
                                        bimodal_boundary_signal |
                                        deep_interior_recombinants |
                                        complex_rearrangement_out_of_scope

CONSENSUS QC (from per_method_tsv):
  q3_method_disagreement_kb            largest distance_to_consensus among
                                        methods in agreement band
  q3_sv_method_flagged                 bool — was SV method (weight 0.5)
                                        outside agreement band?

=============================================================================
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path


# -----------------------------------------------------------------------------
# Minimal stats
# -----------------------------------------------------------------------------
def median(xs):
    xs = sorted(x for x in xs if x is not None and x == x)
    if not xs:
        return float("nan")
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def mad(xs):
    xs = [x for x in xs if x is not None and x == x]
    if len(xs) < 2:
        return float("nan")
    m = median(xs)
    return 1.4826 * median([abs(x - m) for x in xs])


def skewness_direction(xs):
    """
    Return 'right_skewed' if the tail extends right of the mode.
    Using: mean > median implies right skew. Approximate check only.
    """
    xs = [x for x in xs if x is not None and x == x]
    if len(xs) < 5:
        return "degenerate"
    m_mean = sum(xs) / len(xs)
    m_med = median(xs)
    m_mad = mad(xs)
    if m_mad == 0 or not (m_mad == m_mad):
        return "degenerate"
    z = (m_mean - m_med) / m_mad
    if abs(z) < 0.2:
        return "unimodal"  # approximately symmetric
    return "right_skewed" if z > 0 else "left_skewed"


def hartigan_dip_quick(xs):
    """
    Lightweight dip/bimodality proxy. Returns a p-value-like score.
    Uses the same ECDF-vs-linear approximation as v6's earlier version.
    Low p = bimodal.
    """
    xs = sorted(x for x in xs if x is not None and x == x)
    n = len(xs)
    if n < 10 or xs[0] == xs[-1]:
        return 1.0
    x_min, x_max = xs[0], xs[-1]
    expected_y = [(x - x_min) / (x_max - x_min) for x in xs]
    ecdf_y = [(i + 1) / n for i in range(n)]
    dip = max(abs(oy - ey) for oy, ey in zip(ecdf_y, expected_y))
    k = dip * math.sqrt(n)
    return min(1.0, max(0.0, 2 * math.exp(-2 * k * k)))


# -----------------------------------------------------------------------------
# Loaders
# -----------------------------------------------------------------------------
def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_int(s, default=None):
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return default


def to_float(s, default=None):
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


# -----------------------------------------------------------------------------
# Derivation: precision class
# -----------------------------------------------------------------------------
def precision_class(left_ci_kb, right_ci_kb):
    """
    Worst of the two sides determines the class.
    """
    if left_ci_kb is None or right_ci_kb is None:
        return "unknown"
    worst = max(left_ci_kb, right_ci_kb)
    if worst < 5:
        return "snp_resolution"
    if worst < 20:
        return "sub_5kb"
    if worst < 100:
        return "sub_30kb"
    return "coarse"


# -----------------------------------------------------------------------------
# Derivation: fragment distribution shape → interior class
# -----------------------------------------------------------------------------
def interior_class(
    left_shape, right_shape,
    left_bimod_p, right_bimod_p,
    recomb_frac, n_carriers,
):
    """
    Apply rules in order to decide interior class.
    """
    if n_carriers is None or n_carriers < 10:
        return "degenerate"

    # Check order matters. Recombinant fraction is the primary signal;
    # fragment-distribution shape is a secondary diagnostic of those
    # recombinants' ancestry.
    #
    # Note: the stdlib dip proxy used here is NOT a real Hartigan dip
    # test — it's an ECDF-vs-linear approximation. When a rigorous
    # bimodality test is required (manuscript), run R's `diptest`
    # package on the per-side fragment distributions. Conservative
    # threshold (< 0.01) avoids false positives from the proxy.

    # Rule 1: out-of-scope gate — most carriers are recombinant
    if recomb_frac is not None and recomb_frac > 0.70:
        return "complex_rearrangement_out_of_scope"

    # Rule 2: deep recombinants — substantial but recoverable
    if recomb_frac is not None and recomb_frac > 0.40:
        return "deep_interior_recombinants"

    # Rule 3: bimodal boundary — conservative dip threshold
    if left_bimod_p < 0.01 or right_bimod_p < 0.01:
        return "bimodal_boundary_signal"

    # Rule 4: edge recombinants
    if recomb_frac is not None and 0.05 <= recomb_frac <= 0.40:
        return "edge_recombinants_only"

    # Rule 5: clean
    return "clean_simple"


# -----------------------------------------------------------------------------
# Core derivation per candidate
# -----------------------------------------------------------------------------
def derive_keys(cid, consensus_row, frag_rows, per_method_rows):
    """
    Compute the derived keys for one candidate.
    """
    data = {}

    # --- Consensus (from consensus_tsv) ---
    if consensus_row:
        fl = to_int(consensus_row.get("final_left_bp"))
        fr = to_int(consensus_row.get("final_right_bp"))
        lcl = to_int(consensus_row.get("left_ci_low"))
        lch = to_int(consensus_row.get("left_ci_high"))
        rcl = to_int(consensus_row.get("right_ci_low"))
        rch = to_int(consensus_row.get("right_ci_high"))
        lcw = to_float(consensus_row.get("left_ci_width_kb"))
        rcw = to_float(consensus_row.get("right_ci_width_kb"))

        data.update({
            "q3_final_left_bp": fl,
            "q3_final_right_bp": fr,
            "q3_left_ci_low_bp": lcl,
            "q3_left_ci_high_bp": lch,
            "q3_right_ci_low_bp": rcl,
            "q3_right_ci_high_bp": rch,
            "q3_left_ci_width_kb": lcw,
            "q3_right_ci_width_kb": rcw,
            "q3_n_methods_agreeing_left": to_int(
                consensus_row.get("n_methods_agreeing_left"), 0),
            "q3_n_methods_agreeing_right": to_int(
                consensus_row.get("n_methods_agreeing_right"), 0),
            "q3_primary_source": consensus_row.get("primary_source"),
            "q3_sources_available": consensus_row.get("sources_available"),
            "q3_inversion_size_bp": to_int(
                consensus_row.get("inversion_size_bp")),
            "q3_breakpoint_precision_class": precision_class(lcw, rcw),
        })
    else:
        # No consensus — bridge can't provide breakpoints
        data["q3_breakpoint_precision_class"] = "unknown"

    # --- Fragment distribution (from fragments_tsv) ---
    if frag_rows:
        starts = [to_int(r.get("frag_start_bp")) for r in frag_rows]
        ends = [to_int(r.get("frag_end_bp")) for r in frag_rows]
        lengths = [to_int(r.get("frag_length_bp")) for r in frag_rows]
        starts = [x for x in starts if x is not None]
        ends = [x for x in ends if x is not None]
        lengths = [x for x in lengths if x is not None]

        # Compute block width from consensus if available; else use
        # max fragment end - min fragment start as a proxy
        block_width = None
        if consensus_row:
            fl = to_int(consensus_row.get("final_left_bp"))
            fr = to_int(consensus_row.get("final_right_bp"))
            if fl is not None and fr is not None:
                block_width = fr - fl
        if block_width is None and lengths:
            block_width = max(lengths)

        recomb_frac = None
        if block_width and lengths:
            short_threshold = 0.80 * block_width
            recomb_frac = sum(1 for l in lengths if l < short_threshold) / len(lengths)

        left_shape = skewness_direction(starts)
        right_shape = skewness_direction(ends)
        left_bimod_p = hartigan_dip_quick(starts)
        right_bimod_p = hartigan_dip_quick(ends)
        n_carriers = len(frag_rows)

        ic = interior_class(
            left_shape, right_shape, left_bimod_p, right_bimod_p,
            recomb_frac, n_carriers
        )

        data.update({
            "q2_n_fragment_carriers": n_carriers,
            "q2_fragment_interior_recomb_fraction": (
                round(recomb_frac, 4) if recomb_frac is not None else None
            ),
            "q2_fragment_left_distribution_shape": left_shape,
            "q2_fragment_right_distribution_shape": right_shape,
            "q2_fragment_left_bimodality_p": (
                round(left_bimod_p, 4) if left_bimod_p is not None else None
            ),
            "q2_fragment_right_bimodality_p": (
                round(right_bimod_p, 4) if right_bimod_p is not None else None
            ),
            "q2_interior_class": ic,
        })
    else:
        data["q2_interior_class"] = "degenerate"

    # --- Consensus QC (from per_method_tsv) ---
    if per_method_rows:
        in_band = [r for r in per_method_rows
                   if str(r.get("within_agreement_band", "")).upper() in
                   ("TRUE", "T", "1", "YES")]
        distances_in_band = [
            abs(to_float(r.get("distance_to_consensus_kb"), 0))
            for r in in_band
        ]
        distances_in_band = [d for d in distances_in_band if d == d]
        if distances_in_band:
            data["q3_method_disagreement_kb"] = round(
                max(distances_in_band), 2)

        # SV method flag — was the SV method outside agreement band?
        sv_methods = [r for r in per_method_rows
                      if "sv" in (r.get("method") or "").lower() or
                      "step37" in (r.get("method") or "").lower()]
        if sv_methods:
            sv_outside = any(
                str(r.get("within_agreement_band", "")).upper() not in
                ("TRUE", "T", "1", "YES")
                for r in sv_methods)
            data["q3_sv_method_flagged"] = sv_outside

    return data


# -----------------------------------------------------------------------------
# Registry
# -----------------------------------------------------------------------------
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
        print(f"[bp_bridge] registry unavailable: {e}", file=sys.stderr)
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0])
    p.add_argument("--consensus_tsv", required=True,
                   help="candidate_breakpoints_consensus.tsv")
    p.add_argument("--fragments_tsv", required=True,
                   help="candidate_ancestral_fragments.tsv")
    p.add_argument("--per_method_tsv", default="",
                   help="candidate_breakpoints_per_method.tsv (optional)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--registries_root", default="")
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    consensus_rows = load_tsv(args.consensus_tsv)
    frag_rows = load_tsv(args.fragments_tsv)
    per_method_rows = load_tsv(args.per_method_tsv) if args.per_method_tsv else []

    print(f"[bp_bridge] consensus rows : {len(consensus_rows)}")
    print(f"[bp_bridge] fragment rows  : {len(frag_rows)}")
    print(f"[bp_bridge] per-method rows: {len(per_method_rows)}")

    # Index by candidate_id
    consensus_by_cid = {r["candidate_id"]: r for r in consensus_rows
                        if r.get("candidate_id")}
    frags_by_cid = defaultdict(list)
    for r in frag_rows:
        if r.get("candidate_id"):
            frags_by_cid[r["candidate_id"]].append(r)
    per_method_by_cid = defaultdict(list)
    for r in per_method_rows:
        if r.get("candidate_id"):
            per_method_by_cid[r["candidate_id"]].append(r)

    all_cids = set(consensus_by_cid.keys()) | set(frags_by_cid.keys())
    print(f"[bp_bridge] total candidates: {len(all_cids)}")

    reg = try_registry(args.registries_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_interior = defaultdict(int)
    summary_precision = defaultdict(int)

    for cid in sorted(all_cids):
        consensus = consensus_by_cid.get(cid)
        frags = frags_by_cid.get(cid, [])
        per_method = per_method_by_cid.get(cid, [])

        data = derive_keys(cid, consensus, frags, per_method)

        summary_interior[data.get("q2_interior_class", "unknown")] += 1
        summary_precision[data.get("q3_breakpoint_precision_class",
                                     "unknown")] += 1

        cand_out = outdir / cid / "structured"
        cand_out.mkdir(parents=True, exist_ok=True)
        if not args.dry_run:
            (cand_out / "fragment_distribution.json").write_text(json.dumps({
                "block_type": "fragment_distribution",
                "candidate_id": cid,
                "source_script": "bp_pipeline_bridge.py",
                "data": data,
            }, indent=2, default=str))
            if reg is not None:
                try:
                    reg.evidence.write_block(
                        candidate_id=cid,
                        block_type="fragment_distribution",
                        data=data,
                        source_script="bp_pipeline_bridge.py",
                    )
                except Exception as e:
                    print(f"[bp_bridge] {cid}: registry write failed: {e}",
                          file=sys.stderr)

    print(f"\n[bp_bridge] interior class distribution:")
    for cls, n in sorted(summary_interior.items(), key=lambda kv: -kv[1]):
        print(f"  {n:5d}  {cls}")
    print(f"\n[bp_bridge] precision class distribution:")
    for cls, n in sorted(summary_precision.items(), key=lambda kv: -kv[1]):
        print(f"  {n:5d}  {cls}")


if __name__ == "__main__":
    main()
