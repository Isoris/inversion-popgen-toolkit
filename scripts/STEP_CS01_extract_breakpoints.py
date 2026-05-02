#!/usr/bin/env python3
# =============================================================================
# STEP_CS01_extract_breakpoints.py
# -----------------------------------------------------------------------------
# Read a wfmash PAF (Cgar haplotype = query, Cmac haplotype = target) and
# emit a cs_breakpoints_v1 JSON cataloguing chromosome-scale rearrangements
# between the two species. Each breakpoint is annotated with flanking repeat
# density from the TEfull layer (already computed per-chromosome upstream by
# aggregate_repeat_density.R) so the manuscript can argue "Spalax-style"
# repeat-mediated mechanism for the rearrangements that are flanked by
# elevated all_TE.
#
# Three event types detected:
#   - inversion         : block on Mac chrom X with strand flipped relative
#                         to its neighbours on the Gar walk
#   - translocation     : block jumps from Mac chrom X to Mac chrom Y mid-Gar-chrom
#                         (or the other way round)
#   - fission_or_fusion : a Gar chrom mostly covered by ONE Mac chrom but
#                         with a sub-segment from a different Mac chrom
#                         (gar27 ↔ mac28 difference forces at least one
#                         fission/fusion event in the lineage)
#
# Usage:
#   python3 STEP_CS01_extract_breakpoints.py \
#     --paf /scratch/.../01_wfmash/gar_vs_mac.paf \
#     --te-dir /scratch/.../Cgar_TE_density \
#     --mac-te-dir /scratch/.../Cmac_TE_density \
#     --out /scratch/.../02_breakpoints/cs_breakpoints_v1.json \
#     --min-mapq 1 --min-block-bp 50000 --cluster-radius-bp 50000 \
#     --flank-bp 100000
#
# Note on --min-mapq: default is 1 because wfmash mapq values are typically
# 1-5 (different scoring scheme than minimap2 which produces 0-60). Use
# --min-mapq 40 only if your PAF was produced by minimap2.
#
# Validation hooks (NOT run by this script — Quentin's hand QC):
#   - 5 random breakpoints → IGV spot-check on Cgar BAM
#   - cross-reference Kuang et al. 2024 catfish synteny figures
#   - overlap with the 226-sample polymorphic inversion candidates
#     (see cs_breakpoints_v1.candidate_overlap field, populated when the
#      candidate TSV is passed via --candidates)
#
# =============================================================================

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone

# pandas import is deferred to main() so a syntax error doesn't prevent
# the --help text from rendering on a misconfigured env.

PAF_COLUMNS = [
    "query_name",          # 0   query sequence name (Gar chrom)
    "query_length",        # 1
    "query_start",         # 2   0-based, inclusive
    "query_end",           # 3   0-based, exclusive
    "strand",              # 4   '+' or '-' (mapping orientation)
    "target_name",         # 5   target sequence name (Mac chrom)
    "target_length",       # 6
    "target_start",        # 7
    "target_end",          # 8
    "n_residue_matches",   # 9   #identical bases in alignment
    "alignment_block_length",  # 10  total alignment length incl. indels
    "mapping_quality",     # 11  0..255
]

SCHEMA_VERSION = 2
TOOL_NAME      = "cross_species_breakpoints_v1"   # tool name unchanged (atlas detector)
# v2 additions vs v1:
#   - synteny_blocks: list of {gar_chr, mac_chr, gar_start, gar_end, mac_start,
#                              mac_end, strand, block_size_bp, mapping_quality}
#     i.e. the post-filter PAF rows in compact form. The atlas uses these to
#     derive macro-synteny (one-to-one / fusion / fission / many-to-many)
#     and to overlay inversion candidates against synteny block edges.
#   - chrom_lengths_query / chrom_lengths_target: per-chrom bp lengths (for
#     normalising "% chromosome covered" in the synteny tier).
# The cs_breakpoints_v1 schema (v1) is forward-compatible: an atlas without
# v2 awareness sees synteny_blocks as an unrecognised key and ignores it.


# =============================================================================
# I/O helpers
# =============================================================================

def _open(path: str, mode: str = "rt"):
    """Transparent gzip open for .gz; passthrough otherwise."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_paf(path: str):
    """Load PAF, return DataFrame with named columns + computed length."""
    import pandas as pd
    print(f"[STEP_CS01] loading PAF: {path}", file=sys.stderr)
    # PAF rows can have optional tags after column 12 ("tp:A:P", "cm:i:...",
    # etc.). usecols truncates so we don't pay the parse cost for tags we
    # don't use. dtype hints keep memory predictable for ~100k-row PAFs.
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=range(12),
        names=PAF_COLUMNS,
        dtype={
            "query_name": str, "target_name": str, "strand": str,
            "query_length": "int64", "query_start": "int64", "query_end": "int64",
            "target_length": "int64", "target_start": "int64", "target_end": "int64",
            "n_residue_matches": "int64", "alignment_block_length": "int64",
            "mapping_quality": "int64",
        },
    )
    print(f"[STEP_CS01] PAF rows: {len(df):,}", file=sys.stderr)
    return df


def load_te_density(te_dir: str) -> dict:
    """Load all *_repeat_density_TEfull.json from a directory; key by chrom.

    Each entry is the un-modified parsed JSON, so the atlas-side data shape
    is identical and we can read the all_TE class for flanking densities
    by indexing chromosomes[0].by_class['all_TE'].densities + window_centers_mb.
    Returns {} if the directory doesn't exist (TEfull layer optional —
    breakpoints still emit, just without flanking annotation).
    """
    out: dict = {}
    if not te_dir or not os.path.isdir(te_dir):
        print(f"[STEP_CS01] no TE dir at {te_dir} — flanking annotation skipped", file=sys.stderr)
        return out
    for fn in sorted(os.listdir(te_dir)):
        if not fn.endswith("_repeat_density_TEfull.json"):
            continue
        path = os.path.join(te_dir, fn)
        try:
            with open(path) as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError) as e:
            print(f"[STEP_CS01] WARN: skipping unreadable {fn}: {e}", file=sys.stderr)
            continue
        if not isinstance(data, dict) or "chromosomes" not in data:
            continue
        if not data["chromosomes"]:
            continue
        block = data["chromosomes"][0]
        chrom = block.get("chrom")
        if chrom:
            out[chrom] = data
    print(f"[STEP_CS01] TE density loaded for {len(out)} chrom(s) from {te_dir}", file=sys.stderr)
    return out


# =============================================================================
# Core: walk PAF rows, flag transitions
# =============================================================================

def filter_paf(df, min_mapq: int, min_block_bp: int):
    """Drop low-quality and tiny alignments. Returns new DataFrame."""
    n0 = len(df)
    df = df[df["mapping_quality"] >= min_mapq]
    df = df[df["alignment_block_length"] >= min_block_bp]
    n1 = len(df)
    print(f"[STEP_CS01] filter: mapq>={min_mapq}, block>={min_block_bp:,} bp "
          f"-> kept {n1:,}/{n0:,}", file=sys.stderr)
    return df.reset_index(drop=True)


def detect_breakpoints(df) -> list:
    """Walk each Gar chrom in query_start order; emit breakpoints at transitions.

    Single-pass logic per Gar chrom:
      For each consecutive pair of blocks (prev, cur):
        - if cur.target_name != prev.target_name:
            -> 'translocation' or 'fission_or_fusion' event, depending on
               whether one of the two Mac chroms dominates the Gar chrom
               (handled in the post-pass that aggregates per-Gar-chrom)
        - elif cur.strand != prev.strand AND cur.target_name == prev.target_name:
            -> 'inversion' event (orientation flip on the same Mac chrom)
        - else: no event (continuous syntenic stretch)

    The breakpoint position on Gar is the midpoint between prev.query_end and
    cur.query_start (the inter-block gap centre — typical wfmash gaps are
    small so this is well-defined). On Mac, we report both the prev block
    end and the cur block start because for inversions the breakpoint
    crosses the orientation flip and there isn't a single "mac_pos".
    """
    breakpoints = []
    bp_counter = 0
    for gar_chr, sub in df.groupby("query_name", sort=False):
        sub = sub.sort_values(["query_start", "query_end"]).reset_index(drop=True)
        if len(sub) < 2:
            continue
        for i in range(1, len(sub)):
            prev = sub.iloc[i - 1]
            cur  = sub.iloc[i]

            # Transition flags
            tgt_changed = cur["target_name"] != prev["target_name"]
            ori_changed = cur["strand"]      != prev["strand"]

            if not tgt_changed and not ori_changed:
                continue

            # Position on Gar: midpoint of the inter-block gap
            gap_start = int(prev["query_end"])
            gap_end   = int(cur["query_start"])
            if gap_end < gap_start:
                # Overlapping blocks (rare but possible at boundaries) —
                # use the overlap centre instead
                gap_start, gap_end = gap_end, gap_start
            gar_pos_start = gap_start
            gar_pos_end   = max(gap_end, gap_start + 1)
            gar_pos_mb    = (gar_pos_start + gar_pos_end) / 2 / 1e6

            # Event type discrimination
            if tgt_changed:
                event_type = "translocation_or_fission"
            else:
                event_type = "inversion"

            bp_counter += 1
            bp = {
                "id":             f"cs_bp_{bp_counter:04d}",
                "event_type":     event_type,
                "gar_chr":        str(gar_chr),
                "gar_pos_start":  int(gar_pos_start),
                "gar_pos_end":    int(gar_pos_end),
                "gar_pos_mb":     round(gar_pos_mb, 4),
                # We report both adjacent blocks in full so the atlas can
                # draw the syntenic-block lines without re-deriving them.
                "prev_block": {
                    "mac_chr":          str(prev["target_name"]),
                    "mac_start_bp":     int(prev["target_start"]),
                    "mac_end_bp":       int(prev["target_end"]),
                    "strand":           str(prev["strand"]),
                    "block_size_bp":    int(prev["alignment_block_length"]),
                    "mapping_quality":  int(prev["mapping_quality"]),
                },
                "next_block": {
                    "mac_chr":          str(cur["target_name"]),
                    "mac_start_bp":     int(cur["target_start"]),
                    "mac_end_bp":       int(cur["target_end"]),
                    "strand":           str(cur["strand"]),
                    "block_size_bp":    int(cur["alignment_block_length"]),
                    "mapping_quality":  int(cur["mapping_quality"]),
                },
            }
            breakpoints.append(bp)
    print(f"[STEP_CS01] raw breakpoints: {len(breakpoints):,}", file=sys.stderr)
    return breakpoints


def cluster_breakpoints(bps: list, radius_bp: int) -> list:
    """Merge adjacent breakpoints (same Gar chrom, within `radius_bp`) into one.

    Breakpoint clusters happen when wfmash returns multiple short alignments
    near a true breakpoint (orientation noise from short syntenic islands
    inside an inverted region, or alignment ends that aren't exactly aligned
    bp-for-bp on both sides). Clustering condenses these into one event.

    Strategy:
      - sort by (gar_chr, gar_pos_start)
      - greedy merge: if the next breakpoint is within radius_bp on the same
        gar_chr, merge into the running cluster
      - cluster's representative metadata: take the first breakpoint's
        prev_block + last breakpoint's next_block (this is the wider span
        that the cluster actually encompasses)
      - cluster's event_type: 'inversion' if all members are inversions on
        the same Mac chrom; 'translocation_or_fission' if any target switch
        happens; 'mixed' otherwise
    """
    if not bps:
        return []
    bps = sorted(bps, key=lambda b: (b["gar_chr"], b["gar_pos_start"]))
    out = []
    cluster: list = [bps[0]]
    def _flush():
        if not cluster:
            return
        first = cluster[0]
        last  = cluster[-1]
        types = {b["event_type"] for b in cluster}
        if len(types) == 1:
            event_type = next(iter(types))
        else:
            event_type = "mixed"
        merged = {
            "id":             first["id"] if len(cluster) == 1 else f"{first['id']}__{last['id']}",
            "event_type":     event_type,
            "gar_chr":        first["gar_chr"],
            "gar_pos_start":  first["gar_pos_start"],
            "gar_pos_end":    last["gar_pos_end"],
            "gar_pos_mb":     round((first["gar_pos_start"] + last["gar_pos_end"]) / 2 / 1e6, 4),
            "n_member_breakpoints": len(cluster),
            "prev_block": first["prev_block"],
            "next_block": last["next_block"],
        }
        out.append(merged)
    for bp in bps[1:]:
        head = cluster[-1]
        same_chr = (bp["gar_chr"] == head["gar_chr"])
        within   = (bp["gar_pos_start"] - head["gar_pos_end"]) <= radius_bp
        if same_chr and within:
            cluster.append(bp)
        else:
            _flush()
            cluster = [bp]
    _flush()
    print(f"[STEP_CS01] clustered breakpoints: {len(out):,} "
          f"(radius={radius_bp:,} bp)", file=sys.stderr)
    return out


def reclassify_post_cluster(bps: list, df) -> list:
    """Promote 'translocation_or_fission' to a sharper label using whole-chrom context.

    For each Gar chrom, look at the mix of Mac chroms covered:
      - if exactly 2 Mac chroms cover the Gar chrom AND one of them
        contributes <50% of the aligned bp -> 'fission_or_fusion'
        (one species lineage gained or lost a chrom-fission boundary here)
      - if >=3 Mac chroms or the split is more even -> 'translocation'
    Inversions stay as 'inversion' regardless. 'mixed' clusters keep that
    label but get a `dominant_event` field for the most common member.
    """
    # Per-Gar-chrom Mac coverage map (in aligned bp)
    cov = defaultdict(lambda: defaultdict(int))
    for _, row in df.iterrows():
        cov[row["query_name"]][row["target_name"]] += int(row["alignment_block_length"])
    for bp in bps:
        if bp["event_type"] == "inversion":
            continue
        if bp["event_type"] not in ("translocation_or_fission", "mixed"):
            continue
        gar = bp["gar_chr"]
        macs = cov.get(gar, {})
        if not macs:
            bp["event_type_refined"] = bp["event_type"]
            continue
        total = sum(macs.values())
        sorted_macs = sorted(macs.items(), key=lambda x: -x[1])
        dominant_frac = sorted_macs[0][1] / total if total > 0 else 0.0
        n_macs = len([m for m, v in sorted_macs if v >= 0.05 * total])
        if n_macs <= 2 and dominant_frac >= 0.5:
            refined = "fission_or_fusion"
        else:
            refined = "translocation"
        bp["event_type_refined"] = refined
        bp["gar_chr_mac_cov"] = {m: int(v) for m, v in sorted_macs[:5]}
    return bps


# =============================================================================
# Repeat-density flanking annotation (TEfull layer)
# =============================================================================

def _density_in_window(te_block: dict, cls: str, start_bp: int, end_bp: int):
    """Mean density of `cls` over windows whose centre lies in [start, end].

    Returns (mean, max, n_windows) tuple; (None, None, 0) if no overlap or
    the class is missing. Honors the by_class layout from aggregate_repeat_
    density.R: by_class[cls].densities[i] is parallel to window_centers_mb[i].
    """
    if not te_block or "by_class" not in te_block:
        return None, None, 0
    if cls not in te_block["by_class"]:
        return None, None, 0
    centers_bp = [int(c * 1e6) for c in te_block.get("window_centers_mb", [])]
    densities  = te_block["by_class"][cls].get("densities", [])
    if not centers_bp or not densities:
        return None, None, 0
    vals = []
    for cb, d in zip(centers_bp, densities):
        if d is None:
            continue
        if not isinstance(d, (int, float)):
            continue
        if cb < start_bp or cb > end_bp:
            continue
        vals.append(float(d))
    if not vals:
        return None, None, 0
    return (sum(vals) / len(vals), max(vals), len(vals))


def annotate_flanking_repeats(bps: list, gar_te: dict, mac_te: dict,
                              flank_bp: int, classes_to_annotate: list) -> list:
    """For each breakpoint, compute mean+max density in the ±flank window
    on both Gar and Mac, for the requested classes.

    Output shape per breakpoint:
      flanking_repeat_density_gar: { class: {mean, max, n_windows} }
      flanking_repeat_density_mac: { class: {mean, max, n_windows} }

    The flank is symmetric around gar_pos (or the prev/next block edges on
    Mac). For the Mac side, we use the closer of (prev_block.mac_end_bp,
    next_block.mac_start_bp) as the anchor on each Mac chrom, since the
    breakpoint is split across two Mac positions for translocations.
    """
    for bp in bps:
        gar_chr = bp["gar_chr"]
        gar_anchor = (bp["gar_pos_start"] + bp["gar_pos_end"]) // 2

        # ---- Gar side ----
        gar_block = (gar_te.get(gar_chr) or {}).get("chromosomes", [{}])[0] \
                    if gar_te.get(gar_chr) else {}
        gar_out = {}
        for cls in classes_to_annotate:
            mean, mx, n = _density_in_window(
                gar_block, cls,
                gar_anchor - flank_bp, gar_anchor + flank_bp,
            )
            if mean is not None:
                gar_out[cls] = {
                    "mean": round(mean, 4),
                    "max":  round(mx, 4),
                    "n_windows": n,
                }
        bp["flanking_repeat_density_gar"] = gar_out

        # ---- Mac side ----
        # For inversions the prev/next blocks share a Mac chrom: anchor at
        # the boundary between them (prev.end ≈ next.start). For translocations
        # the two blocks are on different Mac chroms; report each.
        prev = bp["prev_block"]
        nxt  = bp["next_block"]
        mac_out = {}
        for side, blk in (("prev", prev), ("next", nxt)):
            mac_chr = blk["mac_chr"]
            mac_anchor = blk["mac_end_bp"] if side == "prev" else blk["mac_start_bp"]
            mac_block = (mac_te.get(mac_chr) or {}).get("chromosomes", [{}])[0] \
                        if mac_te.get(mac_chr) else {}
            side_out = {}
            for cls in classes_to_annotate:
                mean, mx, n = _density_in_window(
                    mac_block, cls,
                    mac_anchor - flank_bp, mac_anchor + flank_bp,
                )
                if mean is not None:
                    side_out[cls] = {
                        "mean": round(mean, 4),
                        "max":  round(mx, 4),
                        "n_windows": n,
                    }
            mac_out[side] = {
                "mac_chr": mac_chr,
                "anchor_bp": int(mac_anchor),
                "by_class": side_out,
            }
        bp["flanking_repeat_density_mac"] = mac_out

        # ---- Manuscript note: flag Spalax-style enrichment on all_TE ----
        gar_all_te = gar_out.get("all_TE", {})
        if gar_all_te:
            # Compare the local mean to the chrom-wide mean (computed once
            # below — store it per-chrom in a memo)
            chrom_mean = _chrom_all_te_mean(gar_block) if gar_block else None
            if chrom_mean and gar_all_te["mean"] >= 1.5 * chrom_mean:
                bp["manuscript_note"] = (
                    f"Spalax-style enrichment: gar_flank all_TE mean "
                    f"{gar_all_te['mean']:.2f} vs chrom mean {chrom_mean:.2f} "
                    f"({gar_all_te['mean'] / chrom_mean:.1f}x)"
                )
    return bps


_chrom_all_te_mean_cache: dict = {}

def _chrom_all_te_mean(block: dict):
    """Memoized chrom-wide mean of all_TE density."""
    if not block or "by_class" not in block:
        return None
    chrom = block.get("chrom") or id(block)
    if chrom in _chrom_all_te_mean_cache:
        return _chrom_all_te_mean_cache[chrom]
    densities = (block.get("by_class") or {}).get("all_TE", {}).get("densities", [])
    vals = [float(d) for d in densities if isinstance(d, (int, float))]
    mean = sum(vals) / len(vals) if vals else None
    if mean is not None:
        mean = round(mean, 4)
    _chrom_all_te_mean_cache[chrom] = mean
    return mean


# =============================================================================
# Optional: overlap with the 226-sample polymorphic inversion candidates
# =============================================================================

def overlap_polymorphic_candidates(bps: list, candidates_path: str) -> list:
    """For each breakpoint, flag whether it falls within an inversion candidate
    interval from the 226-sample cohort. The candidates TSV has at minimum
    the columns: candidate_id, chrom, start_bp, end_bp.

    This is the manuscript's "Angle 3" claim — locations that produce
    species-level rearrangements over evolutionary time are CURRENTLY
    POLYMORPHIC in the hatchery cohort, supporting the recurrent-hotspot
    mechanism. We annotate per breakpoint: candidate_overlap = list of
    candidate_id strings whose [start,end] contains gar_pos.
    """
    if not candidates_path or not os.path.exists(candidates_path):
        return bps
    import pandas as pd
    df = pd.read_csv(candidates_path, sep="\t")
    needed = {"candidate_id", "chrom", "start_bp", "end_bp"}
    if not needed.issubset(df.columns):
        print(f"[STEP_CS01] WARN: candidate TSV missing cols {needed - set(df.columns)} — skipping overlap",
              file=sys.stderr)
        return bps
    by_chrom = defaultdict(list)
    for _, row in df.iterrows():
        by_chrom[str(row["chrom"])].append((
            int(row["start_bp"]), int(row["end_bp"]), str(row["candidate_id"]),
        ))
    n_overlap = 0
    for bp in bps:
        gar = bp["gar_chr"]
        anchor = (bp["gar_pos_start"] + bp["gar_pos_end"]) // 2
        hits = []
        for s, e, cid in by_chrom.get(gar, []):
            if s <= anchor <= e:
                hits.append(cid)
        if hits:
            bp["candidate_overlap"] = hits
            n_overlap += 1
    print(f"[STEP_CS01] candidate overlap: {n_overlap}/{len(bps)} breakpoints "
          f"coincide with 226-sample polymorphic candidate intervals",
          file=sys.stderr)
    return bps


# =============================================================================
# Main
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Cross-species breakpoint extraction from a wfmash PAF "
                    "(Cgar query, Cmac target). Emits cs_breakpoints_v1 JSON.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--paf", required=True, help="wfmash PAF (Gar query, Mac target)")
    p.add_argument("--out", required=True, help="output cs_breakpoints_v1 JSON path")
    p.add_argument("--te-dir", default=None,
                   help="Cgar TE-density JSON dir (one *_repeat_density_TEfull.json per chrom)")
    p.add_argument("--mac-te-dir", default=None,
                   help="Cmac TE-density JSON dir (one *_repeat_density_TEfull.json per chrom)")
    p.add_argument("--candidates", default=None,
                   help="Optional TSV of 226-sample inversion candidates "
                        "(cols: candidate_id, chrom, start_bp, end_bp). When given, each "
                        "breakpoint is annotated with which candidates it overlaps.")
    p.add_argument("--min-mapq", type=int, default=1,
                   help="drop PAF rows below this mapping quality (default 1, "
                        "appropriate for wfmash whose mapq is typically 1-5; "
                        "use --min-mapq 40 for minimap2-style PAFs)")
    p.add_argument("--min-block-bp", type=int, default=50000,
                   help="drop PAF rows shorter than this in alignment_block_length (default 50000)")
    p.add_argument("--cluster-radius-bp", type=int, default=50000,
                   help="merge adjacent breakpoints within this distance on the same Gar chrom (default 50000)")
    p.add_argument("--flank-bp", type=int, default=100000,
                   help="symmetric flank around each breakpoint for repeat-density mean/max (default 100000)")
    p.add_argument("--classes", default="all_TE,young_TE_all,old_TE_all,Gypsy_LTR_retrotransposon,CACTA_TIR_transposon",
                   help="comma-separated TE classes to annotate flanking density for")
    p.add_argument("--species-query", default="Clarias gariepinus")
    p.add_argument("--species-target", default="Clarias macrocephalus")
    p.add_argument("--haplotype-query", default="fClaHyb_Gar_LG")
    p.add_argument("--haplotype-target", default="fClaHyb_Mac_LG")
    args = p.parse_args()

    classes = [c.strip() for c in args.classes.split(",") if c.strip()]
    t_start = time.time()

    df = load_paf(args.paf)
    df = filter_paf(df, args.min_mapq, args.min_block_bp)

    bps = detect_breakpoints(df)
    bps = cluster_breakpoints(bps, args.cluster_radius_bp)
    bps = reclassify_post_cluster(bps, df)

    gar_te = load_te_density(args.te_dir) if args.te_dir else {}
    mac_te = load_te_density(args.mac_te_dir) if args.mac_te_dir else {}
    bps = annotate_flanking_repeats(bps, gar_te, mac_te, args.flank_bp, classes)

    if args.candidates:
        bps = overlap_polymorphic_candidates(bps, args.candidates)

    # Per-event-type counts for the summary
    counts = defaultdict(int)
    for bp in bps:
        et = bp.get("event_type_refined", bp["event_type"])
        counts[et] += 1

    # ----- v2: synteny_blocks + chrom lengths --------------------------------
    # Compact per-row representation of the post-filter PAF for downstream
    # macro-synteny analysis on the atlas side. We do NOT include match-tag
    # fields (n_residue_matches etc.) — the atlas only needs interval geometry
    # plus strand and mapq for filtering. Sorted by (gar_chr, gar_start) for
    # deterministic output and so the atlas can stream-render without resort.
    synteny_blocks = []
    df_sorted = df.sort_values(["query_name", "query_start", "query_end"]).reset_index(drop=True)
    for _, row in df_sorted.iterrows():
        synteny_blocks.append({
            "gar_chr":         str(row["query_name"]),
            "gar_start":       int(row["query_start"]),
            "gar_end":         int(row["query_end"]),
            "mac_chr":         str(row["target_name"]),
            "mac_start":       int(row["target_start"]),
            "mac_end":         int(row["target_end"]),
            "strand":          str(row["strand"]),
            "block_size_bp":   int(row["alignment_block_length"]),
            "mapping_quality": int(row["mapping_quality"]),
        })

    # Chrom lengths (max query_length / target_length seen per chrom name).
    # Honest about uniqueness: a chrom only appears in the PAF if at least one
    # block aligned, so query_length is the assembly-side full length. We
    # take the max across rows defensively in case wfmash truncated anything.
    chrom_lengths_query = {}
    for _, row in df_sorted.iterrows():
        c = str(row["query_name"]); v = int(row["query_length"])
        if v > chrom_lengths_query.get(c, 0):
            chrom_lengths_query[c] = v
    chrom_lengths_target = {}
    for _, row in df_sorted.iterrows():
        c = str(row["target_name"]); v = int(row["target_length"])
        if v > chrom_lengths_target.get(c, 0):
            chrom_lengths_target[c] = v

    out = {
        "tool":           TOOL_NAME,
        "schema_version": SCHEMA_VERSION,
        "generated_at":   _utc_now_iso(),
        "elapsed_seconds": round(time.time() - t_start, 2),
        "species_query":  {"name": args.species_query, "haplotype": args.haplotype_query},
        "species_target": {"name": args.species_target, "haplotype": args.haplotype_target},
        "input_paf": {
            "path":   os.path.abspath(args.paf),
            "sha256": _sha256(args.paf),
        },
        "params": {
            "min_mapq":          args.min_mapq,
            "min_block_bp":      args.min_block_bp,
            "cluster_radius_bp": args.cluster_radius_bp,
            "flank_bp":          args.flank_bp,
            "classes_annotated": classes,
        },
        "n_breakpoints": len(bps),
        "n_by_event_type": dict(counts),
        "breakpoints": bps,
        # v2 additions:
        "n_synteny_blocks": len(synteny_blocks),
        "synteny_blocks": synteny_blocks,
        "chrom_lengths_query": chrom_lengths_query,
        "chrom_lengths_target": chrom_lengths_target,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=False)
    print(f"[STEP_CS01] wrote {args.out} ({len(bps)} breakpoints, "
          f"{len(synteny_blocks)} synteny blocks, "
          f"event types: {dict(counts)})", file=sys.stderr)


if __name__ == "__main__":
    main()
