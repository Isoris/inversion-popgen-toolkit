#!/usr/bin/env python3
"""
test_phase4b_integration.py — end-to-end test for Phase 4b rewrite.

Builds a synthetic Engine B Q cache + candidate table + minimal decomp
outputs, runs STEP_C01i_c_nested_composition.py against them, and checks
that:
  - a clean interval (all samples homogeneous) → composite_flag = clean
  - a composite interval (40% two_block) → composite_flag = likely_composite
  - missing Q cache → composite_flag = unknown_no_engine_b
  - the composite_flag -> promotion_cap mapping works in
    compute_group_validation (mirror of C01f patch)

Does NOT require R. Exercises the Python wrapper + the cap-aware
validation logic.
"""
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY_SCRIPT = HERE.parent / "4b_group_proposal" / "STEP_C01i_c_nested_composition.py"
CORE = HERE.parent / "4b_group_proposal" / "nested_composition_core.py"
SCHEMAS = HERE.parent / "schemas"

# Path to v10 registry_loader.py (needed by the wrapper). Test copies it
# into a temp registries/ root.
V10_LOADER = HERE.parent.parent.parent / "phase4_v10" / "registries" / "api" / "python" / "registry_loader.py"


# =============================================================================
# Compute_group_validation with promotion_cap — Python port of the C01f patch
# =============================================================================
def compute_group_validation(
    current_level="UNCERTAIN",
    t8_concordance=None,
    t9_jackknife_status=None,
    t9_max_delta=None,
    t10_theta_concordance=None,
    layer_d_fisher_p=None,
    layer_d_fisher_or=None,
    promotion_cap=None,
):
    def finite(x):
        return x is not None and isinstance(x, (int, float)) and math.isfinite(x)

    if not current_level:
        current_level = "UNCERTAIN"

    is_fragile = (t9_jackknife_status == "fragile"
                  or (finite(t9_max_delta) and t9_max_delta > 0.3))
    if is_fragile:
        return "SUSPECT"

    candidates = []
    if (finite(layer_d_fisher_p) and finite(layer_d_fisher_or)
            and layer_d_fisher_p < 0.05 and layer_d_fisher_or > 5):
        candidates.append("VALIDATED")
    if (finite(t8_concordance) and t8_concordance >= 0.70
            and t9_jackknife_status == "robust"):
        candidates.append("SUPPORTED")

    rank = {"NONE": 0, "UNCERTAIN": 1, "SUPPORTED": 2, "VALIDATED": 3, "SUSPECT": -1}
    best = max([current_level] + candidates, key=lambda l: rank.get(l, 0))

    if promotion_cap and promotion_cap in rank:
        if rank[best] > rank[promotion_cap]:
            best = promotion_cap
    return best


# =============================================================================
# Synthetic Q cache builders
# =============================================================================
def write_q_cache_clean(path, chrom, samples, K=8):
    """All samples homogeneous: single dominant ancestry label."""
    rows = []
    dom_label = 1
    for sid in samples:
        for i, start in enumerate(range(8000000, 11000000, 500000)):
            end = start + 500000
            q = [0.02] * K
            q[dom_label - 1] = 0.86
            rows.append([chrom, start, end, sid, *q, 0.86, 0.80, 0.4, 1.2, dom_label])
    _write_q_file(path, rows, K)


def write_q_cache_composite(path, chrom, samples, composite_frac=0.4, K=8):
    """composite_frac of samples are two_block_composite, rest homogeneous."""
    n_composite = int(round(len(samples) * composite_frac))
    composite_samples = set(samples[:n_composite])
    rows = []
    for sid in samples:
        for i, start in enumerate(range(8000000, 11000000, 500000)):
            end = start + 500000
            q = [0.02] * K
            if sid in composite_samples:
                # Two-block: first half label 1, second half label 2
                midpoint = 9500000
                label = 1 if start < midpoint else 2
            else:
                label = 1
            q[label - 1] = 0.86
            rows.append([chrom, start, end, sid, *q, 0.86, 0.80, 0.4, 1.2, label])
    _write_q_file(path, rows, K)


def _write_q_file(path, rows, K):
    header = ["chrom", "start_bp", "end_bp", "sample_id"] + \
             [f"Q{i+1}" for i in range(K)] + \
             ["max_q", "delta12", "entropy", "ena", "assigned_pop"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


# =============================================================================
# Test runner
# =============================================================================
def run_wrapper(candidates_tsv, q_cache_dir, outdir, extra_env=None):
    """Invoke STEP_C01i_c_nested_composition.py via subprocess."""
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(
        [sys.executable, str(PY_SCRIPT),
         "--candidates", str(candidates_tsv),
         "--q_cache_dir", str(q_cache_dir),
         "--outdir", str(outdir),
         "--tier_max", "3"],
        env=env, capture_output=True, text=True, timeout=60,
    )
    return result


def read_block(outdir, cid):
    """Read the output JSON block (fallback, since we may not have registry set up)."""
    # First try registry layout
    reg_path = Path(outdir).parent / "registries" / "data" / "evidence_registry" \
               / "per_candidate" / cid / "structured" / "internal_ancestry_composition.json"
    if reg_path.exists():
        return json.loads(reg_path.read_text())
    # Fallback: wrapper writes to outdir/<cid>/internal_ancestry_composition.json
    fb_path = Path(outdir) / cid / "internal_ancestry_composition.json"
    if fb_path.exists():
        return json.loads(fb_path.read_text())
    return None


def main():
    tmp = Path(tempfile.mkdtemp(prefix="phase4b_test_"))
    print(f"[test] tmpdir = {tmp}")
    failures = []

    def check(name, cond, detail=""):
        if cond:
            print(f"  ✓ {name}")
        else:
            failures.append((name, detail))
            print(f"  ✗ {name}: {detail}")

    # ─── Setup: copy v10 registry_loader.py + schemas ──
    registries_root = tmp / "registries"
    (registries_root / "api" / "python").mkdir(parents=True)
    if V10_LOADER.exists():
        shutil.copy(V10_LOADER, registries_root / "api" / "python" / "registry_loader.py")
    # Schemas
    shutil.copytree(SCHEMAS, registries_root / "schemas" / "structured_block_schemas")

    env_base = {
        "REGISTRIES": str(registries_root),
        "PYTHONPATH": f"{registries_root / 'api' / 'python'}:{HERE.parent / 'python'}",
    }

    # ─── Test 1: clean interval ─────────────────────────────────────────
    print("\n[test 1] clean interval — all samples homogeneous")
    q_dir = tmp / "q_cache_clean"
    q_dir.mkdir()
    samples = [f"CGA{i:03d}" for i in range(20)]
    write_q_cache_clean(q_dir / "LG01.local_Q_samples.tsv", "LG01", samples)

    candidates = tmp / "candidates_clean.tsv"
    with candidates.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["candidate_id", "chrom", "start_bp", "end_bp", "tier"])
        w.writerow(["LG01_clean", "LG01", "8500000", "10500000", "1"])

    outdir = tmp / "out_clean"
    r = run_wrapper(candidates, q_dir, outdir, env_base)
    check("wrapper exits 0 on clean", r.returncode == 0,
          detail=f"returncode={r.returncode} stderr={r.stderr[:500]}")

    block = read_block(outdir, "LG01_clean")
    check("clean block was written", block is not None)
    if block:
        flag = block.get("data", {}).get("composite_flag")
        check("clean → composite_flag = clean", flag == "clean",
              detail=f"got {flag}")
        bd = block.get("data", {}).get("structure_breakdown", "")
        check("clean structure_breakdown shows homogeneous dominance",
              "homogeneous" in bd, detail=f"bd={bd}")

    # ─── Test 2: composite interval ─────────────────────────────────────
    print("\n[test 2] composite interval — 40% samples two_block")
    q_dir2 = tmp / "q_cache_comp"
    q_dir2.mkdir()
    write_q_cache_composite(q_dir2 / "LG01.local_Q_samples.tsv",
                              "LG01", samples, composite_frac=0.4)

    candidates2 = tmp / "candidates_comp.tsv"
    with candidates2.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["candidate_id", "chrom", "start_bp", "end_bp", "tier"])
        w.writerow(["LG01_comp", "LG01", "8500000", "10500000", "1"])

    outdir2 = tmp / "out_comp"
    r2 = run_wrapper(candidates2, q_dir2, outdir2, env_base)
    check("wrapper exits 0 on composite", r2.returncode == 0,
          detail=f"stderr={r2.stderr[:500]}")

    block2 = read_block(outdir2, "LG01_comp")
    check("composite block was written", block2 is not None)
    if block2:
        flag2 = block2.get("data", {}).get("composite_flag")
        check("composite 40% → likely_composite",
              flag2 == "likely_composite",
              detail=f"got {flag2}")
        pct = block2.get("data", {}).get("pct_two_block_composite", 0)
        check("pct_two_block_composite ≥ 0.20",
              pct and pct >= 0.20, detail=f"pct={pct}")

    # ─── Test 3: missing Q cache → stub ─────────────────────────────────
    print("\n[test 3] missing Q cache → unknown_no_engine_b stub")
    no_q_dir = tmp / "q_cache_missing"  # does not exist
    candidates3 = tmp / "candidates_stub.tsv"
    with candidates3.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["candidate_id", "chrom", "start_bp", "end_bp", "tier"])
        w.writerow(["LG01_stub", "LG01", "8500000", "10500000", "1"])

    outdir3 = tmp / "out_stub"
    r3 = run_wrapper(candidates3, no_q_dir, outdir3, env_base)
    check("wrapper exits 0 on missing cache", r3.returncode == 0,
          detail=f"stderr={r3.stderr[:500]}")

    block3 = read_block(outdir3, "LG01_stub")
    check("stub block was written", block3 is not None)
    if block3:
        flag3 = block3.get("data", {}).get("composite_flag")
        check("missing cache → composite_flag = unknown_no_engine_b",
              flag3 == "unknown_no_engine_b",
              detail=f"got {flag3}")

    # ─── Test 4: promotion_cap enforcement ──────────────────────────────
    print("\n[test 4] promotion_cap enforcement (C01f patch)")
    # Strong evidence would normally promote to VALIDATED
    strong = dict(
        current_level="UNCERTAIN",
        t8_concordance=0.90, t9_jackknife_status="robust",
        layer_d_fisher_p=0.001, layer_d_fisher_or=20,
    )
    check("no cap + strong evidence → VALIDATED",
          compute_group_validation(**strong) == "VALIDATED",
          detail=compute_group_validation(**strong))

    check("cap=UNCERTAIN + strong evidence → UNCERTAIN (capped)",
          compute_group_validation(**strong, promotion_cap="UNCERTAIN") == "UNCERTAIN")

    check("cap=SUPPORTED + strong evidence → SUPPORTED (capped)",
          compute_group_validation(**strong, promotion_cap="SUPPORTED") == "SUPPORTED")

    check("cap=VALIDATED + strong evidence → VALIDATED (no effective cap)",
          compute_group_validation(**strong, promotion_cap="VALIDATED") == "VALIDATED")

    # Cap only acts downward
    weak = dict(current_level="UNCERTAIN")
    check("cap=SUPPORTED + no promotion evidence → stays UNCERTAIN",
          compute_group_validation(**weak, promotion_cap="SUPPORTED") == "UNCERTAIN")

    # SUSPECT demotion bypasses cap
    check("cap=VALIDATED + fragile jackknife → SUSPECT (demotion wins)",
          compute_group_validation(
              current_level="UNCERTAIN", t9_jackknife_status="fragile",
              promotion_cap="VALIDATED") == "SUSPECT")

    # ─── Report ─────────────────────────────────────────────────────────
    if failures:
        print(f"\n✗ {len(failures)} failure(s)")
        shutil.rmtree(tmp)
        return 1
    print(f"\n✓ all integration tests passed")
    shutil.rmtree(tmp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
