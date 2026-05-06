#!/usr/bin/env python3
# =============================================================================
# 02_breakpoint_windows/make_breakpoint_windows.py
# =============================================================================
"""
Builds the breakpoint-centred BED windows that the density step will overlap
against the normalised TE BED. Output is one BED per species (focal first,
plus any species with synteny mappings).

Window types per candidate:
    left_bp_50kb        left_breakpoint  +/- 25_000   (window 50 kb wide)
    left_bp_100kb       +/- 50_000
    left_bp_500kb       +/- 250_000
    right_bp_50kb       right_breakpoint +/- 25_000
    right_bp_100kb
    right_bp_500kb
    inside_interval     [start, end]
    left_flank_outside  [max(0, start - flank), start]      flank = inv_len capped at 1 Mb
    right_flank_outside [end, min(chr_len, end + flank)]
    local_bg_2mb        center = midpoint(left_bp, right_bp), +/- 1_000_000,
                          MINUS the union of the breakpoint windows above
    chrom_bg            entire chromosome  (one row per chromosome touched)

Coordinates are clipped to [0, chr_len) using the species `.fai`.

BED columns:
    chrom  start0  end  candidate_id  window_type  side

`side` is "left", "right", "inside", "flank_left", "flank_right",
"local_bg", or "chrom_bg".

Usage (focal-only, default):
    python make_breakpoint_windows.py \
        --candidates config/candidate_breakpoints.tsv \
        --species_manifest config/species_manifest.tsv \
        --out_dir output/breakpoint_windows

With synteny (also writes per-target-species window files):
    python make_breakpoint_windows.py \
        --candidates config/candidate_breakpoints.tsv \
        --species_manifest config/species_manifest.tsv \
        --synteny config/synteny_mapping.tsv \
        --out_dir output/breakpoint_windows
"""
from __future__ import annotations

import argparse
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
)


WINDOW_RADII = {  # half-widths (bp) for the symmetric breakpoint windows
    "50kb":  25_000,
    "100kb": 50_000,
    "500kb": 250_000,
}
LOCAL_BG_RADIUS = 1_000_000  # ±1 Mb around inversion midpoint


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
def load_chr_lens(fai_path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with fai_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[normalize_chrom(parts[0])] = int(parts[1])
            except ValueError:
                continue
    return out


def clip(s: int, e: int, chr_len: int) -> tuple[int, int]:
    return max(0, s), min(chr_len, e)


def subtract_intervals(a: tuple[int, int], blockers: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return a minus the union of blockers. All half-open [s, e)."""
    if a[0] >= a[1]:
        return []
    pieces = [a]
    for bs, be in sorted(blockers):
        if bs >= be:
            continue
        new_pieces: list[tuple[int, int]] = []
        for ps, pe in pieces:
            if be <= ps or bs >= pe:
                new_pieces.append((ps, pe))
                continue
            if bs > ps:
                new_pieces.append((ps, bs))
            if be < pe:
                new_pieces.append((be, pe))
        pieces = new_pieces
        if not pieces:
            break
    return pieces


# -----------------------------------------------------------------------------
# core
# -----------------------------------------------------------------------------
def windows_for_interval(
    *,
    candidate_id: str,
    chrom: str,
    chr_len: int,
    start_1based: int,
    end_1based: int,
    left_bp_1based: int,
    right_bp_1based: int,
) -> list[tuple]:
    """Return list of (chrom, s0, e0, candidate_id, window_type, side) tuples."""
    rows: list[tuple] = []

    # 1-based inclusive -> 0-based half-open
    iv_s, iv_e = start_1based - 1, end_1based
    lb, rb = left_bp_1based - 1, right_bp_1based  # treat bp as a base; window centred on it

    # symmetric BP windows
    bp_window_intervals: list[tuple[int, int]] = []
    for label, half in WINDOW_RADII.items():
        ls, le = clip(lb - half, lb + half, chr_len)
        rs, re = clip(rb - half, rb + half, chr_len)
        rows.append((chrom, ls, le, candidate_id, f"left_bp_{label}",  "left"))
        rows.append((chrom, rs, re, candidate_id, f"right_bp_{label}", "right"))
        bp_window_intervals.append((ls, le))
        bp_window_intervals.append((rs, re))

    # inside interval
    iv_s_c, iv_e_c = clip(iv_s, iv_e, chr_len)
    rows.append((chrom, iv_s_c, iv_e_c, candidate_id, "inside_interval", "inside"))

    # outside flanks: same length as inversion, capped at 1 Mb
    flank = max(1, min(iv_e - iv_s, 1_000_000))
    fls, fle = clip(iv_s - flank, iv_s, chr_len)
    frs, fre = clip(iv_e, iv_e + flank, chr_len)
    rows.append((chrom, fls, fle, candidate_id, "left_flank_outside",  "flank_left"))
    rows.append((chrom, frs, fre, candidate_id, "right_flank_outside", "flank_right"))

    # local 2 Mb background, MINUS the breakpoint windows (use 100 kb as the
    # "exclusion" set; using all radii would zero out big intervals)
    mid = (lb + rb) // 2
    lbgs, lbge = clip(mid - LOCAL_BG_RADIUS, mid + LOCAL_BG_RADIUS, chr_len)
    excl = []
    half100 = WINDOW_RADII["100kb"]
    excl.append(clip(lb - half100, lb + half100, chr_len))
    excl.append(clip(rb - half100, rb + half100, chr_len))
    for ps, pe in subtract_intervals((lbgs, lbge), excl):
        rows.append((chrom, ps, pe, candidate_id, "local_bg_2mb", "local_bg"))

    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--species_manifest", required=True, type=Path)
    ap.add_argument("--synteny", type=Path, default=None,
                    help="optional cross-species synteny mapping")
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("02_make_breakpoint_windows", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(args.species_manifest,
                            required_cols=["species_id", "fai_path", "is_focal"])
    species_idx = {r["species_id"]: r for r in species_rows}

    # focal species
    focal_rows = [r for r in species_rows if str(r["is_focal"]).lower() == "true"]
    if len(focal_rows) != 1:
        log.error("expected exactly one focal species; got %d", len(focal_rows))
        return 1
    focal_id = focal_rows[0]["species_id"]
    log.info("focal species: %s", focal_id)

    cand_cols = ["candidate_id", "chrom", "start", "end",
                 "left_breakpoint", "right_breakpoint"]
    candidates = read_tsv(args.candidates, required_cols=cand_cols)
    log.info("loaded %d candidates", len(candidates))

    # focal windows
    focal_chr_lens = load_chr_lens(Path(focal_rows[0]["fai_path"]))
    focal_recs: list[tuple] = []
    n_chr_seen: set[str] = set()
    for c in candidates:
        chrom = normalize_chrom(c["chrom"])
        if chrom not in focal_chr_lens:
            log.warning("candidate %s on chrom %s not in focal .fai; skipping",
                        c["candidate_id"], chrom)
            continue
        n_chr_seen.add(chrom)
        focal_recs.extend(windows_for_interval(
            candidate_id=c["candidate_id"],
            chrom=chrom,
            chr_len=focal_chr_lens[chrom],
            start_1based=int(c["start"]),
            end_1based=int(c["end"]),
            left_bp_1based=int(c["left_breakpoint"]),
            right_bp_1based=int(c["right_breakpoint"]),
        ))

    # add chrom_bg row(s) — one per chromosome touched
    for chrom in sorted(n_chr_seen):
        focal_recs.append((chrom, 0, focal_chr_lens[chrom],
                           "_CHROM_BG_", "chrom_bg", "chrom_bg"))

    focal_out = args.out_dir / f"{focal_id}.windows.bed"
    focal_recs.sort(key=lambda t: (t[0], t[1], t[2]))
    with focal_out.open("w") as fh:
        for t in focal_recs:
            fh.write("\t".join(str(x) for x in t) + "\n")
    write_run_manifest(focal_out, species_id=focal_id,
                       n_records=len(focal_recs),
                       n_candidates=len(candidates))
    log.info("focal windows -> %s (%d rows)", focal_out, len(focal_recs))

    # -------------------------------------------------------------------------
    # comparative windows from synteny mapping (optional)
    # -------------------------------------------------------------------------
    if args.synteny and args.synteny.exists():
        synt = read_tsv(args.synteny, required_cols=[
            "candidate_id", "target_species", "target_chrom",
            "target_start", "target_end", "orientation",
        ])
        log.info("loaded %d synteny rows for %d candidates",
                 len(synt), len({r['candidate_id'] for r in synt}))

        by_sp: dict[str, list[dict]] = defaultdict(list)
        for r in synt:
            by_sp[r["target_species"]].append(r)

        for sp, rows in by_sp.items():
            if sp not in species_idx:
                log.warning("synteny target species %s not in species_manifest; skipping", sp)
                continue
            chr_lens = load_chr_lens(Path(species_idx[sp]["fai_path"]))
            recs: list[tuple] = []
            chrs_seen: set[str] = set()
            for r in rows:
                chrom = normalize_chrom(r["target_chrom"])
                if chrom not in chr_lens:
                    log.warning("synteny target %s/%s not in fai; skipping",
                                sp, chrom)
                    continue
                ts = int(r["target_start"])
                te = int(r["target_end"])
                # in target coordinates the "breakpoints" are the boundaries of
                # the homologous block; orientation does not change coordinates
                lb_target = min(ts, te)
                rb_target = max(ts, te)
                recs.extend(windows_for_interval(
                    candidate_id=r["candidate_id"],
                    chrom=chrom,
                    chr_len=chr_lens[chrom],
                    start_1based=lb_target,
                    end_1based=rb_target,
                    left_bp_1based=lb_target,
                    right_bp_1based=rb_target,
                ))
                chrs_seen.add(chrom)
            for chrom in sorted(chrs_seen):
                recs.append((chrom, 0, chr_lens[chrom],
                             "_CHROM_BG_", "chrom_bg", "chrom_bg"))

            out = args.out_dir / f"{sp}.windows.bed"
            recs.sort(key=lambda t: (t[0], t[1], t[2]))
            with out.open("w") as fh:
                for t in recs:
                    fh.write("\t".join(str(x) for x in t) + "\n")
            write_run_manifest(out, species_id=sp, n_records=len(recs),
                               n_candidates=len({r['candidate_id'] for r in rows}))
            log.info("comparative windows %s -> %s (%d rows)", sp, out, len(recs))
    else:
        log.info("no synteny mapping provided — focal-only run")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
