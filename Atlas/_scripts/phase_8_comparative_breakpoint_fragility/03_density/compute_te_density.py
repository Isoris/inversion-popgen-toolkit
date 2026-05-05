#!/usr/bin/env python3
# =============================================================================
# 03_density/compute_te_density.py
# =============================================================================
"""
Compute TE/repeat density for each candidate window in a species.

Inputs:
  --windows    output/breakpoint_windows/<species>.windows.bed
  --te_bed     output/normalized_te/<species>.te.bed.gz
  --fai        path to species .fa.fai (for chrom-wide percentile sliding windows)
  --species_id e.g. C_gariepinus
  --out_json   one JSON per (species, candidate) pair, written under output/density/

For each window we compute:
  - te_bp_overlap         (bp of the window covered by any TE)
  - callable_bp           (currently == window_bp; v0.2 will subtract a callable mask)
  - te_density            te_bp_overlap / callable_bp
  - te_density_by_class   per-bucket densities (DNA, LTR, ..., unknown)
  - fold_vs_chromosome    te_density / chromosome_te_density
  - fold_vs_local_2mb     te_density / local_bg_te_density
  - percentile_chr        rank of te_density among 100-kb sliding windows
                          on the same chromosome (step 50 kb)
  - breakpoint_TE_hotspot_score
                          = mean(fold_vs_chromosome, fold_vs_local_2mb)
                            * percentile_chr   (NaN-safe, capped at 4.0)

Output JSON (one per candidate, per species role: focal or comparative):
{
  "species_id": "C_gariepinus",
  "candidate_id": "LG28_INV_001",
  "is_focal": true,
  "chromosome_background": {...},
  "local_background_2mb": {...},
  "left_breakpoint":  {"window_50kb": {...}, "window_100kb": {...}, "window_500kb": {...}},
  "right_breakpoint": {"window_50kb": {...}, "window_100kb": {...}, "window_500kb": {...}},
  "inside_interval": {...},
  "left_flank_outside": {...},
  "right_flank_outside": {...}
}

Implementation notes:
  - Pure Python; sweep-line overlap. No pyranges/bedtools dependency.
  - Chrom-wide percentile uses a 100-kb / 50-kb step sliding scan once per chrom.
  - All intervals are 0-based half-open BED.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    TE_CLASS_BUCKETS,
    get_logger,
    normalize_chrom,
    open_maybe_gz,
    write_run_manifest,
)


# -----------------------------------------------------------------------------
# data structures
# -----------------------------------------------------------------------------
def load_te_bed(path: Path) -> dict[str, dict]:
    """
    Returns:
        {chrom: {
            'starts': [..],            # sorted starts, all classes combined
            'ends':   [..],            # parallel ends
            'class':  [..],            # parallel class buckets
            'cum_total': [..],         # prefix-sum of total TE bp by record idx
        }}
    For per-class densities we keep parallel arrays and re-scan; cheap because
    chrom record counts are < 1e6 typically.
    """
    by_chr: dict[str, dict] = defaultdict(lambda: {
        "starts": [], "ends": [], "class": []
    })
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            try:
                chrom = parts[0]
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            cls = parts[6] if len(parts) > 6 else "unknown"
            by_chr[chrom]["starts"].append(s)
            by_chr[chrom]["ends"].append(e)
            by_chr[chrom]["class"].append(cls)

    # ensure sorted (input should be sorted by step 01, but defensive)
    for chrom, d in by_chr.items():
        if not d["starts"]:
            continue
        order = sorted(range(len(d["starts"])),
                       key=lambda i: (d["starts"][i], d["ends"][i]))
        d["starts"] = [d["starts"][i] for i in order]
        d["ends"]   = [d["ends"][i]   for i in order]
        d["class"]  = [d["class"][i]  for i in order]
    return by_chr


def load_fai(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[parts[0]] = int(parts[1])
            except ValueError:
                continue
    return out


# -----------------------------------------------------------------------------
# overlap math
# -----------------------------------------------------------------------------
def overlap_bp_in_window(
    *,
    chrom_idx: dict,
    win_s: int,
    win_e: int,
) -> tuple[int, dict[str, int]]:
    """Return (total_overlap_bp, per_class_bp)."""
    starts = chrom_idx["starts"]
    ends   = chrom_idx["ends"]
    cls    = chrom_idx["class"]
    n = len(starts)
    if n == 0 or win_e <= win_s:
        return 0, {b: 0 for b in TE_CLASS_BUCKETS}

    # records with start < win_e
    hi = bisect_left(starts, win_e)
    # we still need to scan back for records that start before win_s but extend in.
    # use end-based find: linear scan backwards while ends[i] could overlap.
    # In practice, fastest and simplest: scan starts[0..hi] and short-circuit on
    # ends[i] <= win_s.
    total = 0
    per_cls: dict[str, int] = {b: 0 for b in TE_CLASS_BUCKETS}

    # Heuristic: use end-sorted list? We didn't build one. Linear scan over [0, hi)
    # is fine for typical N (a few hundred-thousand records max) and few windows.
    # For tighter perf, switch to interval tree later.
    for i in range(hi):
        e = ends[i]
        if e <= win_s:
            continue
        s = starts[i]
        ov = min(e, win_e) - max(s, win_s)
        if ov <= 0:
            continue
        total += ov
        b = cls[i] if cls[i] in per_cls else "unknown"
        per_cls[b] += ov

    return total, per_cls


def union_overlap_bp_multiwindow(
    *,
    chrom_idx: dict,
    intervals: list[tuple[int, int]],
) -> tuple[int, int, dict[str, int]]:
    """
    Compute (union_window_bp, union_overlap_bp, per_class_overlap_bp) for a
    set of windows. Used for `local_bg_2mb` which is split into pieces.
    Windows are merged before measurement.
    """
    if not intervals:
        return 0, 0, {b: 0 for b in TE_CLASS_BUCKETS}

    # merge windows
    intervals = sorted(intervals)
    merged: list[tuple[int, int]] = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))

    win_bp = sum(e - s for s, e in merged)
    total = 0
    per_cls: dict[str, int] = {b: 0 for b in TE_CLASS_BUCKETS}
    for s, e in merged:
        t, p = overlap_bp_in_window(chrom_idx=chrom_idx, win_s=s, win_e=e)
        total += t
        for k, v in p.items():
            per_cls[k] += v
    return win_bp, total, per_cls


# -----------------------------------------------------------------------------
# chromosome-wide percentile via sliding 100-kb windows (step 50 kb)
# -----------------------------------------------------------------------------
def chrom_density_distribution(
    *,
    chrom_idx: dict,
    chr_len: int,
    win: int = 100_000,
    step: int = 50_000,
) -> list[float]:
    if chr_len < win:
        return []
    out: list[float] = []
    s = 0
    while s + win <= chr_len:
        t, _ = overlap_bp_in_window(chrom_idx=chrom_idx, win_s=s, win_e=s + win)
        out.append(t / win)
        s += step
    return sorted(out)


def percentile_of(value: float, sorted_dist: list[float]) -> float:
    if not sorted_dist:
        return float("nan")
    # fraction of values <= value
    idx = bisect_right(sorted_dist, value)
    return idx / len(sorted_dist)


# -----------------------------------------------------------------------------
# block builders
# -----------------------------------------------------------------------------
def make_window_block(
    *,
    chrom_idx: dict,
    win_s: int,
    win_e: int,
    chrom_density: float,
    local_density: float,
    sorted_chr_dist: list[float] | None,
) -> dict:
    win_bp = max(0, win_e - win_s)
    total, per_cls = overlap_bp_in_window(chrom_idx=chrom_idx,
                                          win_s=win_s, win_e=win_e)
    callable_bp = win_bp  # v0.2 will subtract a callable mask
    density = (total / callable_bp) if callable_bp > 0 else 0.0

    fold_chr = (density / chrom_density) if chrom_density > 0 else float("nan")
    fold_loc = (density / local_density) if local_density > 0 else float("nan")
    pct_chr  = percentile_of(density, sorted_chr_dist) if sorted_chr_dist else float("nan")

    score = float("nan")
    folds = [x for x in (fold_chr, fold_loc) if not _isnan(x)]
    if folds and not _isnan(pct_chr):
        score = min(4.0, (sum(folds) / len(folds)) * pct_chr)

    return {
        "window_bp": win_bp,
        "callable_bp": callable_bp,
        "te_bp_overlap": total,
        "te_density": round(density, 6),
        "te_density_by_class": {b: round(v / callable_bp, 6) if callable_bp else 0.0
                                for b, v in per_cls.items()},
        "fold_vs_chromosome": _round_or_none(fold_chr),
        "fold_vs_local_2mb":  _round_or_none(fold_loc),
        "percentile_chr": _round_or_none(pct_chr),
        "breakpoint_TE_hotspot_score": _round_or_none(score),
    }


def _isnan(x: float) -> bool:
    return x != x


def _round_or_none(x: float, ndigits: int = 4):
    if x is None or _isnan(x):
        return None
    return round(x, ndigits)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--te_bed", required=True, type=Path)
    ap.add_argument("--fai", required=True, type=Path)
    ap.add_argument("--species_id", required=True)
    ap.add_argument("--is_focal", action="store_true")
    ap.add_argument("--out_dir", required=True, type=Path,
                    help="e.g. output/density/focal or output/density/comparative")
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("03_compute_te_density", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    chr_lens = load_fai(args.fai)
    log.info("loaded %d chromosomes from .fai", len(chr_lens))

    te_idx = load_te_bed(args.te_bed)
    log.info("loaded TE BED for %d chromosomes", len(te_idx))

    # group window rows by candidate
    cand_windows: dict[str, list[tuple]] = defaultdict(list)
    chrom_bg_intervals: dict[str, tuple[int, int]] = {}
    with args.windows.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, s, e, cand, wtype, side = parts[:6]
            s, e = int(s), int(e)
            chrom_n = normalize_chrom(chrom)
            if cand == "_CHROM_BG_":
                chrom_bg_intervals[chrom_n] = (s, e)
                continue
            cand_windows[cand].append((chrom_n, s, e, wtype, side))

    log.info("read windows for %d candidates (+ %d chrom_bg rows)",
             len(cand_windows), len(chrom_bg_intervals))

    # per-chromosome cache
    chrom_density_cache: dict[str, float] = {}
    chrom_dist_cache: dict[str, list[float]] = {}

    def get_chrom_density(chrom: str) -> float:
        if chrom in chrom_density_cache:
            return chrom_density_cache[chrom]
        if chrom not in chr_lens or chrom not in te_idx:
            chrom_density_cache[chrom] = 0.0
            return 0.0
        s, e = 0, chr_lens[chrom]
        total, _ = overlap_bp_in_window(chrom_idx=te_idx[chrom],
                                        win_s=s, win_e=e)
        d = total / max(1, e - s)
        chrom_density_cache[chrom] = d
        return d

    def get_chrom_dist(chrom: str) -> list[float]:
        if chrom in chrom_dist_cache:
            return chrom_dist_cache[chrom]
        if chrom not in chr_lens or chrom not in te_idx:
            chrom_dist_cache[chrom] = []
            return []
        d = chrom_density_distribution(chrom_idx=te_idx[chrom],
                                       chr_len=chr_lens[chrom])
        chrom_dist_cache[chrom] = d
        return d

    # process each candidate
    n_written = 0
    for cand_id, rows in cand_windows.items():
        # what chromosome?
        chroms = sorted({r[0] for r in rows})
        if len(chroms) != 1:
            log.warning("candidate %s spans multiple chromosomes: %s; using first",
                        cand_id, chroms)
        chrom = chroms[0]

        if chrom not in te_idx:
            log.warning("no TE BED for chrom %s (cand %s); writing nulls",
                        chrom, cand_id)
            te_idx[chrom] = {"starts": [], "ends": [], "class": []}

        chrom_density = get_chrom_density(chrom)
        sorted_dist = get_chrom_dist(chrom)

        # local 2-Mb: union of all local_bg_2mb pieces for this candidate
        local_pieces = [(s, e) for c, s, e, w, side in rows
                        if w == "local_bg_2mb" and c == chrom]
        local_win_bp, local_te_bp, _local_per_cls = union_overlap_bp_multiwindow(
            chrom_idx=te_idx[chrom], intervals=local_pieces,
        )
        local_density = (local_te_bp / local_win_bp) if local_win_bp > 0 else 0.0

        out: dict = {
            "schema_version": "comparative_breakpoint_fragility_v0.1",
            "species_id": args.species_id,
            "candidate_id": cand_id,
            "chrom": chrom,
            "is_focal": bool(args.is_focal),
            "chromosome_background": {
                "chrom_length_bp": chr_lens.get(chrom, 0),
                "te_density": round(chrom_density, 6),
            },
            "local_background_2mb": {
                "window_bp": local_win_bp,
                "te_bp_overlap": local_te_bp,
                "te_density": round(local_density, 6),
            },
            "left_breakpoint":  {},
            "right_breakpoint": {},
            "inside_interval":  None,
            "left_flank_outside": None,
            "right_flank_outside": None,
        }

        # populate per-window blocks
        for c, s, e, wtype, side in rows:
            blk = make_window_block(
                chrom_idx=te_idx[chrom],
                win_s=s, win_e=e,
                chrom_density=chrom_density,
                local_density=local_density,
                sorted_chr_dist=sorted_dist,
            )
            if wtype.startswith("left_bp_"):
                out["left_breakpoint"][f"window_{wtype.split('_')[-1]}"] = blk
            elif wtype.startswith("right_bp_"):
                out["right_breakpoint"][f"window_{wtype.split('_')[-1]}"] = blk
            elif wtype == "inside_interval":
                out["inside_interval"] = blk
            elif wtype == "left_flank_outside":
                out["left_flank_outside"] = blk
            elif wtype == "right_flank_outside":
                out["right_flank_outside"] = blk
            elif wtype == "local_bg_2mb":
                pass  # already merged into local_background_2mb above

        out_path = args.out_dir / f"{cand_id}__{args.species_id}.json"
        out_path.write_text(json.dumps(out, indent=2, sort_keys=False))
        n_written += 1

    write_run_manifest(args.out_dir / "_density.run.json",
                       species_id=args.species_id,
                       is_focal=bool(args.is_focal),
                       n_candidates=n_written)
    log.info("wrote %d density JSONs to %s", n_written, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
