#!/usr/bin/env python3
"""
STEP_C01i_c_nested_composition.py — Phase 4b.3 of v10.1 Phase 4b rewrite

Thin wrapper around nested_composition_core (the vendored MODULE_2B helper).
Runs the composition analysis and writes results as Tier-2 blocks via the
v10 registry library. Produces the internal_ancestry_composition.json block
and derives composite_flag ∈ {clean, maybe_composite, likely_composite,
unknown_no_engine_b}.

This is the script that catches "same interval, different sample groups"
cases — it reads Engine B's per-window ancestry labels and measures whether
the candidate interval has one coherent grouping or is internally composite.

INPUTS:
  --candidates     TSV with candidate_id, chrom, start_bp, end_bp (+ tier)
  --q_cache_dir    Engine B local-Q cache directory (.local_Q_samples.tsv.gz)
  --outdir         fallback output root
  --K              number of ancestry components (default 8)
  --tier_max       default 3

BEHAVIOR:
  - If q_cache_dir is absent or empty: write stub block with
    composite_flag = "unknown_no_engine_b", exit 0. Downstream treats this
    as clean for validation-gating purposes but records the missing info.
  - Otherwise: for each candidate, call nested_composition_core on the
    interval, then derive composite_flag from structure_breakdown.

OUTPUTS:
  per candidate: internal_ancestry_composition.json Tier-2 block
  flat keys extracted: q1_ancestry_composite_flag, q2_pct_samples_two_block,
                       q2_pct_samples_multi_block, q2_mean_fragmentation,
                       q2_mean_internal_entropy

THRESHOLDS for composite_flag (derived from per-sample structure_type
counts):
  clean             ≥80% of samples are homogeneous OR dominant_plus_secondary
  likely_composite  >20% two_block_composite OR >20% multi_block_fragmented
  maybe_composite   everything in between
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

# Vendor path: nested_composition_core.py is next to this script
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Registry library path
REGISTRIES_API = HERE.parent.parent.parent / "phase4_v10" / "registries" / "api" / "python"
if not REGISTRIES_API.exists():
    # Try other sibling paths (when deployed, registry_loader.py will be in
    # ${REGISTRIES}/api/python/; this file is in phase4b_rewrite/python/)
    for candidate in [
        HERE.parent.parent / "phase4_v10" / "registries" / "api" / "python",
        Path(os.environ.get("REGISTRIES", "")) / "api" / "python",
        HERE.parent / "registries" / "api" / "python",
    ]:
        if candidate.exists():
            REGISTRIES_API = candidate
            break

if REGISTRIES_API.exists():
    sys.path.insert(0, str(REGISTRIES_API))

try:
    from registry_loader import load_registry
except ImportError:
    load_registry = None

from nested_composition_core import (
    classify_structure,
    block_analysis,
    load_q_samples,
    load_parents,
)


# =============================================================================
# Composite flag derivation
# =============================================================================
def derive_composite_flag(structure_counts: Counter, n_samples: int) -> Dict:
    """Convert per-sample structure counts into a per-candidate flag.

    Returns a dict with composite_flag, pct_two_block, pct_multi_block,
    pct_homogeneous, pct_gradient, pct_diffuse, rationale.
    """
    if n_samples == 0:
        return {
            "composite_flag": "unknown_no_data",
            "rationale": "no samples with ancestry data in interval",
        }

    pct = {k: structure_counts.get(k, 0) / n_samples for k in [
        "homogeneous", "dominant_plus_secondary", "two_block_composite",
        "continuous_gradient", "multi_block_fragmented", "diffuse_mixed",
        "too_sparse",
    ]}

    pct_homog_like = pct["homogeneous"] + pct["dominant_plus_secondary"]
    pct_two_block = pct["two_block_composite"]
    pct_multi_block = pct["multi_block_fragmented"]
    pct_gradient = pct["continuous_gradient"]
    pct_diffuse = pct["diffuse_mixed"]

    # Decision tree (order matters)
    if pct_two_block > 0.20 or pct_multi_block > 0.20:
        flag = "likely_composite"
        rationale = (f"{pct_two_block:.0%} two_block_composite, "
                     f"{pct_multi_block:.0%} multi_block_fragmented samples")
    elif pct_two_block >= 0.05 or pct_multi_block >= 0.05 or pct_gradient >= 0.20:
        flag = "maybe_composite"
        rationale = (f"{pct_two_block:.0%} two_block, "
                     f"{pct_multi_block:.0%} multi_block, "
                     f"{pct_gradient:.0%} gradient")
    elif pct_homog_like >= 0.80:
        flag = "clean"
        rationale = f"{pct_homog_like:.0%} homogeneous or dominant_plus_secondary"
    else:
        # Fallback: if nothing dominates strongly, call maybe_composite
        flag = "maybe_composite"
        rationale = "no dominant pattern, defaulting to maybe_composite"

    return {
        "composite_flag": flag,
        "rationale": rationale,
        "pct_homogeneous_like": round(pct_homog_like, 3),
        "pct_two_block_composite": round(pct_two_block, 3),
        "pct_multi_block_fragmented": round(pct_multi_block, 3),
        "pct_continuous_gradient": round(pct_gradient, 3),
        "pct_diffuse_mixed": round(pct_diffuse, 3),
    }


# =============================================================================
# Per-parent analysis (rewrite of main() from nested_composition_core, but
# returning structured output for one parent and emitting registry blocks)
# =============================================================================
import math
from collections import defaultdict


def analyze_parent(parent: Dict, q_rows: List[Dict],
                     K: int = 8) -> Optional[Dict]:
    """Run composition analysis on one parent; return block data dict."""
    pid = parent["parent_id"]
    chrom = parent["chrom"]
    pstart = parent["start"]
    pend = parent["end"]

    # Filter windows overlapping the parent region
    overlap = [r for r in q_rows
               if int(r.get("start_bp", 0)) < pend
               and int(r.get("end_bp", 0)) > pstart]
    if not overlap:
        return None

    # Group by sample, sort by position
    by_sample = defaultdict(list)
    for row in overlap:
        sid = row.get("sample_id", row.get("sample_idx", ""))
        by_sample[sid].append(row)

    per_sample = []
    for sid, rows in by_sample.items():
        rows.sort(key=lambda r: int(r.get("start_bp", 0)))
        labels = [str(r.get("assigned_pop", "")) for r in rows]
        n_sub = len(rows)
        if n_sub < 3:
            continue
        blocks, switches = block_analysis(labels)
        label_counts = Counter(labels)
        mc = label_counts.most_common()
        dom_label = mc[0][0]
        dom_frac = mc[0][1] / n_sub
        probs = [c / n_sub for c in label_counts.values()]
        H = -sum(p * math.log(p + 1e-15) for p in probs)
        frag = switches / max(1, n_sub - 1)
        struct = classify_structure(labels)
        sorted_blocks = sorted(blocks, key=lambda x: -x[1])
        largest = sorted_blocks[0][1] if sorted_blocks else 0
        second = sorted_blocks[1][1] if len(sorted_blocks) > 1 else 0

        mean_d12 = sum(float(r.get("delta12", 0)) for r in rows) / n_sub
        mean_H = sum(float(r.get("entropy", 0)) for r in rows) / n_sub
        mean_ena = sum(float(r.get("ena", 0)) for r in rows) / n_sub

        per_sample.append({
            "sample_id": sid,
            "n_windows": n_sub,
            "n_labels_unique": len(label_counts),
            "dominant_label": dom_label,
            "dominant_fraction": round(dom_frac, 4),
            "n_blocks": len(blocks),
            "n_switches": switches,
            "largest_block": largest,
            "second_block": second,
            "fragmentation_score": round(frag, 4),
            "internal_entropy": round(H, 4),
            "structure_type": struct,
            "mean_delta12": round(mean_d12, 4),
            "mean_entropy": round(mean_H, 4),
            "mean_ena": round(mean_ena, 4),
        })

    if not per_sample:
        return None

    # Aggregate
    struct_counts = Counter(r["structure_type"] for r in per_sample)
    struct_breakdown = "; ".join(f"{k}:{v}" for k, v in struct_counts.most_common())
    mean_frag = sum(r["fragmentation_score"] for r in per_sample) / len(per_sample)
    mean_int_H = sum(r["internal_entropy"] for r in per_sample) / len(per_sample)
    mean_switches = sum(r["n_switches"] for r in per_sample) / len(per_sample)

    flag_info = derive_composite_flag(struct_counts, len(per_sample))

    return {
        "parent_id": pid,
        "chrom": chrom,
        "start_bp": pstart,
        "end_bp": pend,
        "K_used": K,
        "n_samples_analyzed": len(per_sample),
        "n_windows_in_parent": len({int(r["start_bp"]) for r in overlap}),
        "dominant_structure_type": struct_counts.most_common(1)[0][0],
        "dominant_structure_count": struct_counts.most_common(1)[0][1],
        "structure_breakdown": struct_breakdown,
        "structure_counts": dict(struct_counts),
        "mean_fragmentation": round(mean_frag, 4),
        "mean_internal_entropy": round(mean_int_H, 4),
        "mean_switches": round(mean_switches, 2),
        "composite_flag": flag_info["composite_flag"],
        "composite_rationale": flag_info["rationale"],
        "pct_homogeneous_like": flag_info["pct_homogeneous_like"],
        "pct_two_block_composite": flag_info["pct_two_block_composite"],
        "pct_multi_block_fragmented": flag_info["pct_multi_block_fragmented"],
        "pct_continuous_gradient": flag_info["pct_continuous_gradient"],
        "pct_diffuse_mixed": flag_info["pct_diffuse_mixed"],
        "per_sample": per_sample,
    }


# =============================================================================
# Stub block for missing Engine B
# =============================================================================
def write_stub_block(reg, candidate_id: str, parent: Dict, outdir: Path,
                      reason: str) -> None:
    data = {
        "parent_id": candidate_id,
        "chrom": parent["chrom"],
        "start_bp": parent["start"],
        "end_bp": parent["end"],
        "K_used": 0,
        "n_samples_analyzed": 0,
        "dominant_structure_type": "unknown",
        "structure_breakdown": "",
        "composite_flag": "unknown_no_engine_b",
        "composite_rationale": reason,
        "pct_homogeneous_like": None,
        "pct_two_block_composite": None,
        "pct_multi_block_fragmented": None,
        "per_sample": [],
    }

    if reg is not None and hasattr(reg.evidence, "write_block"):
        os.environ["CURRENT_SCRIPT"] = "STEP_C01i_c_nested_composition.py"
        reg.evidence.write_block(
            candidate_id=candidate_id,
            block_type="internal_ancestry_composition",
            data=data,
        )
    else:
        # Fallback: write to outdir
        import json
        cand_outdir = outdir / candidate_id
        cand_outdir.mkdir(parents=True, exist_ok=True)
        with (cand_outdir / "internal_ancestry_composition.json").open("w") as f:
            json.dump({
                "block_type": "internal_ancestry_composition",
                "candidate_id": candidate_id,
                "source_script": "STEP_C01i_c_nested_composition.py",
                "data": data,
            }, f, indent=2, default=str)


# =============================================================================
# Main
# =============================================================================
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="Candidate TSV from C01d")
    ap.add_argument("--q_cache_dir", default="",
                    help="Engine B local Q cache directory")
    ap.add_argument("--outdir", default="nested_comp_out",
                    help="Fallback output root")
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--tier_max", type=int, default=3)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load registry
    reg = None
    if load_registry is not None:
        try:
            reg = load_registry()
            print(f"[nested_comp] registry loaded at {reg.root}")
        except Exception as e:
            print(f"[nested_comp] registry load failed: {e}")

    # Load candidates
    parents = load_parents(args.candidates)
    # Filter by tier if present in candidate table
    filtered = []
    with open(args.candidates) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tier = int(row.get("tier", 2))
            if tier <= args.tier_max:
                filtered.append(row)
    # Rebuild parent list aligned with filtered rows
    if filtered:
        parents = [
            {
                "parent_id": row.get("candidate_id") or row.get("interval_id")
                             or f"{row.get('chrom','')}_{row.get('start_bp','')}_{row.get('end_bp','')}",
                "chrom": row.get("chrom", row.get("chr", "")),
                "start": int(row.get("start_bp", row.get("start", 0))),
                "end": int(row.get("end_bp", row.get("end", 0))),
            }
            for row in filtered
        ]

    print(f"[nested_comp] processing {len(parents)} candidates")

    # Engine B cache check — 5-second smoke test before long compute
    q_cache = args.q_cache_dir or os.environ.get("Q_CACHE_DIR", "")
    q_cache_available = bool(q_cache) and os.path.isdir(q_cache)

    if q_cache_available and not os.environ.get("SKIP_SMOKE_TEST"):
        # Try reading the first parent's chromosome to confirm cache format is right.
        # A common failure: cache dir exists but contains leftover files from
        # a previous run, or assigned_pop column is missing.
        import gzip
        sample_chr = parents[0]["chrom"] if parents else "LG01"
        smoke_ok = False
        for ext in [".local_Q_samples.tsv.gz", ".local_Q_samples.tsv"]:
            path = os.path.join(q_cache, sample_chr + ext)
            if os.path.isfile(path):
                opener = gzip.open if path.endswith(".gz") else open
                try:
                    with opener(path, "rt") as f:
                        header = f.readline()
                        if "assigned_pop" not in header:
                            print(f"[smoke_test] ✗ {path} missing assigned_pop column "
                                  "(run ancestry_bridge.R --prepare)")
                            q_cache_available = False
                            break
                        first = f.readline()
                        if not first:
                            print(f"[smoke_test] ✗ {path} is empty")
                            q_cache_available = False
                            break
                        smoke_ok = True
                        print(f"[smoke_test] ✓ Q cache format OK ({sample_chr})")
                        break
                except Exception as e:
                    print(f"[smoke_test] ✗ cannot read {path}: {e}")
                    q_cache_available = False
                    break
        if not smoke_ok and q_cache_available:
            # No file for sample_chr but cache dir exists — warn and fall through
            # (other chromosomes might have data)
            print(f"[smoke_test] warning: no {sample_chr}.local_Q_samples in cache; "
                  "will attempt per-chromosome reads anyway")

    if not q_cache_available:
        print(f"[nested_comp] Q cache directory not found or invalid: {q_cache or '<unset>'}")
        print("[nested_comp] writing stub blocks with composite_flag=unknown_no_engine_b")
        print("[nested_comp] To enable composite detection: run")
        print("    Rscript utils/ancestry_bridge.R --prepare --K 8")
        for parent in parents:
            write_stub_block(reg, parent["parent_id"], parent, outdir,
                              "Engine B local-Q cache not available")
        return 0

    # Cache Q data per chromosome (expensive file reads)
    q_cache_per_chr: Dict[str, List[Dict]] = {}

    for parent in parents:
        pid = parent["parent_id"]
        chrom = parent["chrom"]
        if chrom not in q_cache_per_chr:
            q_cache_per_chr[chrom] = load_q_samples(q_cache, chrom)
            if not q_cache_per_chr[chrom]:
                print(f"[nested_comp] no Q data for {chrom}")

        q_rows = q_cache_per_chr[chrom]
        if not q_rows:
            write_stub_block(reg, pid, parent, outdir,
                              f"no local-Q data for chromosome {chrom}")
            continue

        result = analyze_parent(parent, q_rows, K=args.K)
        if result is None:
            write_stub_block(reg, pid, parent, outdir,
                              "insufficient per-window data in interval")
            continue

        print(f"[nested_comp] {pid}: {result['composite_flag']} "
              f"({result['structure_breakdown']})")

        # Write the full block
        if reg is not None and hasattr(reg.evidence, "write_block"):
            os.environ["CURRENT_SCRIPT"] = "STEP_C01i_c_nested_composition.py"
            reg.evidence.write_block(
                candidate_id=pid,
                block_type="internal_ancestry_composition",
                data=result,
            )
        else:
            import json
            cand_outdir = outdir / pid
            cand_outdir.mkdir(parents=True, exist_ok=True)
            with (cand_outdir / "internal_ancestry_composition.json").open("w") as f:
                json.dump({
                    "block_type": "internal_ancestry_composition",
                    "candidate_id": pid,
                    "source_script": "STEP_C01i_c_nested_composition.py",
                    "data": result,
                }, f, indent=2, default=str)

    print("\n[nested_comp] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
