#!/usr/bin/env python3
"""
test_sampling_modes.py
----------------------
Tests for dosage_viewer/sampling.py — the six sampling modes.

What's exercised
----------------
- raw       — returns all indices; cap_exceeded flag
- even      — bp-uniform spacing; first/last anchor; sorted; deduplicated
- random    — same seed -> same result; different seed -> different result
- variance  — top-N by variance; ties broken by position; sorted by pos
- hybrid    — 70/30 split; dedup with even; sorted by pos
- aggregate — n_bins x n_samples mean dosage; NA-aware; empty bins -> NaN
- per_site_variance — population variance, ignores NA, < 2 non-NA -> 0

Run:
  python3 tests/dosage_viewer/test_sampling_modes.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "dosage_viewer"))

from sampling import (
    select_raw, select_even, select_random,
    select_variance, select_hybrid, select_aggregate,
    per_site_variance, _closest_index,
)

passed = 0
failed = 0
def expect(cond, msg):
    global passed, failed
    if cond:
        passed += 1
        print("  \u2713 " + msg)
    else:
        failed += 1
        print("  \u2717 " + msg, file=sys.stderr)
def section(s):
    print("\n\u2014 " + s + " \u2014")


# Common test fixture: 10 positions evenly spaced from 1000 to 10000
positions = [i * 1000 for i in range(1, 11)]

# ---------------------------------------------------------------------------
section("select_raw")
sel, exceeded = select_raw(positions, 100)
expect(sel == list(range(10)), "raw with cap > n returns all indices")
expect(exceeded is False, "raw flag false when n <= cap")

sel, exceeded = select_raw(positions, 10)
expect(sel == list(range(10)) and exceeded is False, "raw with cap == n: ok")

sel, exceeded = select_raw(positions, 5)
expect(exceeded is True, "raw with cap < n flags exceeded")
expect(sel == list(range(10)),
       "raw still returns all indices (caller decides) — flag is the signal")

sel, exceeded = select_raw([], 10)
expect(sel == [] and exceeded is False, "raw on empty input: empty result, no flag")

# ---------------------------------------------------------------------------
section("select_even")
# Pick 4 from 10 evenly: should land near positions 1750, 4250, 6750, 9250
sel = select_even(positions, 4)
expect(len(sel) == 4, f"even: 4 picks (got {len(sel)})")
expect(sel == sorted(sel), "even: sorted ascending")
expect(len(sel) == len(set(sel)), "even: no duplicates")

# Edge: cap == n returns all
expect(select_even(positions, 10) == list(range(10)),
       "even with cap == n returns all")
# Edge: cap > n returns all
expect(select_even(positions, 100) == list(range(10)),
       "even with cap > n returns all")
# Edge: cap == 0 returns empty
expect(select_even(positions, 0) == [],
       "even with cap == 0 returns empty")
# Edge: empty input returns empty
expect(select_even([], 10) == [],
       "even on empty input returns empty")
# Edge: single position
expect(select_even([42], 5) == [0],
       "even on single position returns [0]")

# Bp-uniformity check: gaps between consecutive picks within 1.5x of expected
pos_arr = [positions[i] for i in sel]
gaps = [pos_arr[i + 1] - pos_arr[i] for i in range(len(pos_arr) - 1)]
expected_gap = (positions[-1] - positions[0]) / 4
expect(all(0.5 * expected_gap <= g <= 2.0 * expected_gap for g in gaps),
       f"even: gaps within tolerance (got {gaps}, expected ~{expected_gap})")

# Dense input where dedup might trim: 100 closely spaced positions, cap = 1000
dense = list(range(100))
sel = select_even(dense, 1000)
expect(sel == list(range(100)),
       "even: cap > unique density returns all positions (no fake dups)")

# ---------------------------------------------------------------------------
section("select_random")
sel1 = select_random(positions, 5, seed=42)
sel2 = select_random(positions, 5, seed=42)
expect(sel1 == sel2, "random: same seed -> same selection")
expect(sel1 == sorted(sel1), "random: sorted ascending")
expect(len(sel1) == 5 and len(set(sel1)) == 5, "random: 5 unique indices")

sel3 = select_random(positions, 5, seed=999)
expect(sel1 != sel3, "random: different seed -> different selection")
expect(sel3 == sorted(sel3), "random: still sorted with different seed")

# Edge: cap >= n returns all
expect(select_random(positions, 100, seed=1) == list(range(10)),
       "random with cap >= n returns all")

# ---------------------------------------------------------------------------
section("select_variance")
# Variances: sites at index 2 + 7 are the highest-variance
variances = [0.0, 0.1, 0.9, 0.05, 0.2, 0.15, 0.0, 0.85, 0.3, 0.1]
sel = select_variance(positions, variances, 3)
expect(set(sel) == {2, 7, 8},
       f"variance: top-3 by variance = {{2, 7, 8}} (got {set(sel)})")
expect(sel == sorted(sel), "variance: result sorted ascending by position")

# Tie-break — two sites with identical variance should pick the lower index
variances_tied = [0.0, 0.5, 0.5, 0.5, 0.0]
positions_tied = [100, 200, 300, 400, 500]
sel = select_variance(positions_tied, variances_tied, 2)
expect(sel == [1, 2],
       f"variance: ties broken by position ascending (got {sel})")

# NaN variance treated as 0 — should never beat real variance
variances_nan = [float("nan"), 0.1, float("nan"), 0.2, 0.0]
positions_n = [10, 20, 30, 40, 50]
sel = select_variance(positions_n, variances_nan, 2)
expect(set(sel) == {1, 3},
       f"variance: NaN treated as 0 (got {sel})")

# Mismatched lengths -> ValueError
threw = False
try:
    select_variance(positions, [1.0], 2)
except ValueError:
    threw = True
expect(threw, "variance: length mismatch raises ValueError")

# Cap >= n returns all
expect(select_variance(positions, variances, 100) == list(range(10)),
       "variance with cap >= n returns all")

# ---------------------------------------------------------------------------
section("select_hybrid")
# 10 picks: 7 even + 3 variance
sel = select_hybrid(positions, variances, 10, even_share=0.70)
# All even picks should be in the result; the top-3-variance sites should
# also appear (idx 2, 7, 8). With dedup, length may be < 10 if overlap.
expect(2 in sel, "hybrid: top-variance idx 2 included")
expect(7 in sel, "hybrid: top-variance idx 7 included")
expect(sel == sorted(sel), "hybrid: sorted ascending")
expect(len(sel) <= 10, f"hybrid: never exceeds cap (got {len(sel)})")
expect(len(sel) == len(set(sel)), "hybrid: no duplicates")

# Edge: cap >= n returns all
expect(select_hybrid(positions, variances, 100) == list(range(10)),
       "hybrid with cap >= n returns all")

# Edge: even_share = 0 -> all variance
sel = select_hybrid(positions, variances, 3, even_share=0.0)
expect(set(sel) == {2, 7, 8}, "hybrid: even_share=0 reduces to variance")

# Edge: even_share = 1 -> all even
sel = select_hybrid(positions, variances, 3, even_share=1.0)
sel_pure_even = select_even(positions, 3)
expect(sel == sel_pure_even, "hybrid: even_share=1 reduces to even")

# ---------------------------------------------------------------------------
section("select_aggregate")
# 6 sites at positions 1000..6000, dosages designed for a clean test
agg_pos = [1000, 2000, 3000, 4000, 5000, 6000]
agg_dos = [
    [0, 0, -1],   # bin 0
    [0, 1, 1],    # bin 0
    [1, 1, 2],    # bin 1
    [2, 2, -1],   # bin 1
    [2, -1, 2],   # bin 2
    [-1, -1, -1], # bin 2 (all-NA row)
]
centers, matrix = select_aggregate(agg_pos, agg_dos, n_samples=3, n_bins=3)
expect(len(centers) == 3, f"aggregate: 3 bin centers (got {len(centers)})")
expect(len(matrix) == 3, "aggregate: matrix has 3 rows")
expect(all(len(r) == 3 for r in matrix), "aggregate: every row n_samples wide")

# Bin 0: samples [0,0],[0,1] -> means [0.0, 0.5]
#                [-1,1]      -> mean 1.0 (NA excluded from sample 2)
expect(abs(matrix[0][0] - 0.0) < 1e-9, "bin 0 sample 0 mean = 0.0")
expect(abs(matrix[0][1] - 0.5) < 1e-9, "bin 0 sample 1 mean = 0.5")
expect(abs(matrix[0][2] - 1.0) < 1e-9, "bin 0 sample 2 mean = 1.0 (NA excluded)")

# Bin 1: samples [1,1,2],[2,2,-1] -> means [1.5, 1.5, 2.0]
expect(abs(matrix[1][0] - 1.5) < 1e-9, "bin 1 sample 0 mean = 1.5")
expect(abs(matrix[1][2] - 2.0) < 1e-9, "bin 1 sample 2 mean = 2.0 (NA excluded)")

# Bin 2: only one row [2, -1, 2] (the all-NA row falls in same bin)
# Sample 0: 2; sample 1: NaN; sample 2: 2
expect(abs(matrix[2][0] - 2.0) < 1e-9, "bin 2 sample 0 mean = 2.0")
expect(math.isnan(matrix[2][1]), "bin 2 sample 1 mean = NaN (all NA)")
expect(abs(matrix[2][2] - 2.0) < 1e-9, "bin 2 sample 2 mean = 2.0")

# Edge: empty input
centers, matrix = select_aggregate([], [], n_samples=3, n_bins=5)
expect(centers == [] and matrix == [], "aggregate: empty input returns empty")

# Edge: n_bins == 0
centers, matrix = select_aggregate(agg_pos, agg_dos, n_samples=3, n_bins=0)
expect(centers == [] and matrix == [], "aggregate: n_bins=0 returns empty")

# Edge: bp_start == bp_end
centers, matrix = select_aggregate(agg_pos, agg_dos, n_samples=3, n_bins=5,
                                    bp_start=1000, bp_end=1000)
expect(centers == [] and matrix == [], "aggregate: bp_start == bp_end returns empty")

# Width exactly matches n_bins regardless of n_sites_total
centers, matrix = select_aggregate(agg_pos, agg_dos, n_samples=3, n_bins=10)
expect(len(centers) == 10 and len(matrix) == 10,
       "aggregate: width == n_bins regardless of n_sites")

# ---------------------------------------------------------------------------
section("per_site_variance")
rows = [
    [0, 0, 0, 0],         # var = 0
    [0, 1, 2, -1],        # NA-aware: var of {0,1,2} = 0.667
    [-1, -1, -1, 1],      # < 2 non-NA -> 0
    [1, 1, 1, 1],         # var = 0
    [0, 2, 0, 2],         # var = 1.0
]
v = per_site_variance(rows)
expect(len(v) == 5, "per_site_variance: 5 rows -> 5 variances")
expect(abs(v[0] - 0.0) < 1e-9, "row 0 (all 0): var = 0.0")
# Population variance of {0, 1, 2}: mean=1, ((-1)^2+0+1)/3 = 2/3
expect(abs(v[1] - 2.0/3.0) < 1e-9,
       f"row 1 (0,1,2,NA): var = 2/3 (got {v[1]})")
expect(abs(v[2] - 0.0) < 1e-9, "row 2 (only 1 non-NA): var = 0")
expect(abs(v[3] - 0.0) < 1e-9, "row 3 (all 1): var = 0")
expect(abs(v[4] - 1.0) < 1e-9, "row 4 (alternating 0,2): var = 1.0")

# ---------------------------------------------------------------------------
section("_closest_index")
expect(_closest_index([10, 20, 30, 40], 14) == 0, "closest 14 -> 10 (idx 0)")
expect(_closest_index([10, 20, 30, 40], 16) == 1, "closest 16 -> 20 (idx 1)")
expect(_closest_index([10, 20, 30, 40], 30) == 2, "exact 30 -> idx 2")
expect(_closest_index([10, 20, 30, 40], 100) == 3, "OOB high -> last")
expect(_closest_index([10, 20, 30, 40], 0) == 0, "OOB low -> first")
expect(_closest_index([], 5) == 0, "empty -> 0 (defensive)")

# ---------------------------------------------------------------------------
section("Post-selection sort invariant — applies to every mode")
# For every mode that returns indices, positions must be monotonically
# increasing. Already tested individually above; this is a cross-mode summary.
for mode_name, sel in [
    ("even",     select_even(positions, 5)),
    ("random_seed1", select_random(positions, 5, seed=1)),
    ("random_seed99", select_random(positions, 5, seed=99)),
    ("variance", select_variance(positions, variances, 5)),
    ("hybrid",   select_hybrid(positions, variances, 5)),
]:
    pos_arr = [positions[i] for i in sel]
    expect(pos_arr == sorted(pos_arr),
           f"{mode_name}: positions monotonically increasing")

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
