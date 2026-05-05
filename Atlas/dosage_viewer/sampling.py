"""
sampling.py — site-selection modes for the dosage viewer.

All six modes from S6 spec §Sampling modes implemented as pure functions:

    raw       — return all sites; if n_sites > max_sites, signal cap_exceeded
    even      — pick max_sites positions evenly across the bp range
    random    — reproducible random subset using a seed
    variance  — top max_sites sites by per-site dosage variance
    hybrid    — 70 % even + 30 % variance, dedup, sort by pos
    aggregate — bin into n_bins equal-bp bins, mean dosage per bin per sample

The contract for every mode:

    select_*  takes positions: list[int] (ascending) and returns a
              list[int] of selected indices INTO `positions`. Variance/
              hybrid additionally take `variances`. Aggregate is special —
              it returns binned means rather than indices.

After selection, EVERY mode resorts the result by genomic position. The
spec's "Critical contract (applies to all modes)" rule. The variance
and random modes naturally produce out-of-order indices; we sort.

These functions never raise — callers decide whether cap_exceeded is
fatal. The handler / endpoint module wraps with HTTPException.
"""
from __future__ import annotations

import math
import random as _random
from typing import List, Optional, Tuple


def select_raw(positions: List[int], max_sites: int) -> Tuple[List[int], bool]:
    """Return all indices; second value = True if cap was exceeded."""
    n = len(positions)
    if n <= max_sites:
        return list(range(n)), False
    return list(range(n)), True   # caller decides what to do with cap_exceeded


def select_even(positions: List[int], max_sites: int) -> List[int]:
    """Pick max_sites positions evenly across [positions[0], positions[-1]].

    "Even" here is in genomic-bp space (not array-index space). Bins the
    range into max_sites bins and picks the closest position to each bin's
    midpoint. Ties resolve to the lower index. Duplicates are dropped, so
    the returned list may be shorter than max_sites if the position
    density is uneven (e.g. a long gap with few sites).

    Always returns a sorted, deduplicated list of indices.
    """
    n = len(positions)
    if n <= max_sites:
        return list(range(n))
    if max_sites <= 0:
        return []
    lo = positions[0]
    hi = positions[-1]
    if hi <= lo:
        return list(range(min(max_sites, n)))
    span = hi - lo
    out: List[int] = []
    last_idx = -1
    for k in range(max_sites):
        target = lo + (k + 0.5) * span / max_sites
        idx = _closest_index(positions, target)
        if idx == last_idx:
            idx = idx + 1 if idx + 1 < n else idx
        out.append(idx)
        last_idx = idx
    seen = set()
    dedup: List[int] = []
    for i in out:
        if i not in seen and 0 <= i < n:
            seen.add(i)
            dedup.append(i)
    dedup.sort()
    return dedup


def select_random(positions: List[int], max_sites: int,
                  seed: int = 1) -> List[int]:
    """Reproducible random subset. Same seed -> same selection.

    Always returns indices sorted ascending — the spec's post-selection
    sort invariant applies.
    """
    n = len(positions)
    if n <= max_sites:
        return list(range(n))
    rng = _random.Random(seed)
    indices = rng.sample(range(n), max_sites)
    indices.sort()
    return indices


def select_variance(positions: List[int], variances: List[float],
                    max_sites: int) -> List[int]:
    """Top max_sites by variance descending; tie-break by position ascending.
    Re-sorted by position before return.

    `variances` must be len(positions) — the per-site dosage variance.
    Caller computes this. NaN variance is treated as 0 (uninformative).
    """
    n = len(positions)
    if len(variances) != n:
        raise ValueError(
            f"variance/positions length mismatch: {len(variances)} vs {n}"
        )
    if n <= max_sites:
        return list(range(n))
    # Replace NaN with 0 so they never displace real-variance sites
    safe_var = [(0.0 if (v is None or math.isnan(v)) else v) for v in variances]
    # Sort by (-variance, pos) so highest variance wins, with positional
    # tie-break for determinism.
    order = sorted(range(n),
                   key=lambda i: (-safe_var[i], positions[i]))[:max_sites]
    order.sort()    # final positional sort
    return order


def select_hybrid(positions: List[int], variances: List[float],
                  max_sites: int, seed: int = 1,
                  even_share: float = 0.70) -> List[int]:
    """70 % even + 30 % variance, deduplicated, sorted by position.

    The 70/30 split is documented in the spec table; configurable via
    `even_share` for future tweaking. After dedup the returned list MAY
    be shorter than max_sites if the variance picks already overlap with
    the even picks (which is common — high-variance sites cluster, and
    the even mesh might already cover them).
    """
    n = len(positions)
    if n <= max_sites:
        return list(range(n))
    n_even = int(round(max_sites * even_share))
    n_var = max_sites - n_even
    if n_even <= 0:
        return select_variance(positions, variances, max_sites)
    if n_var <= 0:
        return select_even(positions, max_sites)
    even_picks = set(select_even(positions, n_even))
    var_picks = set(select_variance(positions, variances, n_var))
    merged = sorted(even_picks | var_picks)
    return merged


def select_aggregate(positions: List[int], dosage_rows: List[List[int]],
                      n_samples: int, n_bins: int,
                      bp_start: Optional[int] = None,
                      bp_end: Optional[int] = None
                      ) -> Tuple[List[int], List[List[float]]]:
    """Bin the region into n_bins equal-bp bins; per bin per sample,
    return mean dosage (NA-aware: -1 cells are excluded from the mean).

    Returns (bin_centers, matrix), where matrix is n_bins x n_samples.
    Empty bins yield NaN cells (the standalone viewer's frontend renders
    NaN as a missing-tile colour).

    bp_start/bp_end default to positions[0]/positions[-1]. Empty input
    returns ([], []).
    """
    if not positions:
        return [], []
    if bp_start is None:
        bp_start = positions[0]
    if bp_end is None:
        bp_end = positions[-1]
    if bp_end <= bp_start or n_bins <= 0:
        return [], []

    span = bp_end - bp_start
    bin_w = span / n_bins
    # bin_centers — placed at the midpoint of each bin
    centers = [int(round(bp_start + (k + 0.5) * bin_w))
               for k in range(n_bins)]

    sums = [[0.0] * n_samples for _ in range(n_bins)]
    counts = [[0] * n_samples for _ in range(n_bins)]

    for i, pos in enumerate(positions):
        if pos < bp_start or pos > bp_end:
            continue
        bin_idx = int((pos - bp_start) / bin_w) if bin_w > 0 else 0
        if bin_idx >= n_bins:
            bin_idx = n_bins - 1
        if bin_idx < 0:
            bin_idx = 0
        row = dosage_rows[i]
        for j in range(n_samples):
            v = row[j] if j < len(row) else -1
            if v < 0:
                continue   # NA — skip
            sums[bin_idx][j] += v
            counts[bin_idx][j] += 1

    matrix: List[List[float]] = []
    for k in range(n_bins):
        binrow: List[float] = []
        for j in range(n_samples):
            if counts[k][j] == 0:
                binrow.append(float("nan"))
            else:
                binrow.append(sums[k][j] / counts[k][j])
        matrix.append(binrow)
    return centers, matrix


def _closest_index(sorted_vals: List[int], target: float) -> int:
    """Index of the value in `sorted_vals` closest to `target`.
    Identical to dosage_bridge._closest_index — duplicated to avoid a
    cross-package import (keeps this file standalone)."""
    n = len(sorted_vals)
    if n == 0:
        return 0
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        return 0
    if lo == n:
        return n - 1
    a = abs(sorted_vals[lo - 1] - target)
    b = abs(sorted_vals[lo] - target)
    return (lo - 1) if a <= b else lo


# =============================================================================
# Per-site variance (utility used by variance / hybrid modes)
# =============================================================================

def per_site_variance(dosage_rows: List[List[int]]) -> List[float]:
    """Compute the per-site variance of the dosage values, ignoring NA (-1).

    Used by select_variance / select_hybrid. Pure function — no numpy
    dependency so this stays portable.

    Returns 0.0 for sites with < 2 non-NA samples (insufficient data).
    """
    out: List[float] = []
    for row in dosage_rows:
        n = 0
        s = 0.0
        ss = 0.0
        for v in row:
            if v < 0:
                continue
            n += 1
            s += v
            ss += v * v
        if n < 2:
            out.append(0.0)
            continue
        mean = s / n
        # Population variance (we're describing the cohort, not estimating
        # a population from a sample — every fish in the cohort is observed)
        var = ss / n - mean * mean
        if var < 0:
            var = 0.0
        out.append(var)
    return out
