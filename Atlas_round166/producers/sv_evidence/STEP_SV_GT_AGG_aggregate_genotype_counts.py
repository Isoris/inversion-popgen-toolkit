#!/usr/bin/env python3
"""
STEP_SV_GT_AGG_aggregate_genotype_counts.py
-------------------------------------------
Aggregate per-SV genotype counts within a candidate inversion's flanking +
boundary + body zones. Emits sv_genotype_counts_v1.json into the per-candidate
folder.

For each SV call within the candidate window:
  - Bin samples into the three karyotype groups H1/H1, H1/H2, H2/H2 (from
    locked_labels), compute the AA / AB / BB / miss counts per group.
  - Run a Fisher's exact test on the H1/H1 vs H2/H2 contingency to get
    OR + raw p-value. Repeat across all SVs in the candidate, then apply
    Benjamini-Hochberg correction to get FDR_BH.
  - Apply a pattern-label decision rule (zone × OR direction × het pattern)
    to classify each SV as canonical_breakpoint_marker /
    dominant_presence_marker / het_specific_marker / sub_haplotype_marker /
    internal_linked_marker / uninformative.

Inputs:
  --vcf         path to merged DELLY+Manta VCF (or any VCF with GT field)
                or a tabular TSV with sv_id, sample_id, GT
  --candidate   path to candidate JSON (must include candidate_id, chrom,
                boundary_left_bp, boundary_right_bp, optional zone_definitions_bp)
  --karyotype   path to karyotype labels TSV (sample_id, label) where label
                is one of HOMO_1 / HET / HOMO_2 (the page-21 lock format)
  --out-root    output root directory; layers written to
                <out>/<chrom>/candidates/<candidate_id>/sv_genotype_counts.json
  --fdr-cutoff  optional FDR threshold used to classify "associated" (default 0.05)
  --indent      JSON indent (default: compact one-line; use 2 for debug)

Cohort note: this producer is for the **226-sample pure C. gariepinus
hatchery cohort** only. K clusters reflect hatchery broodlines, NOT
species admixture. Do NOT use against the F1 hybrid cohort or the
C. macrocephalus wild cohort.

Author: Quentin Andres
Project: MS_Inversions_North_african_catfish
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

# Local sibling import
sys.path.insert(0, str(Path(__file__).parent))
from write_candidate_folder import write_candidate_layer, update_chrom_manifest, log_write


# Karyotype label remapping: page 21 / scrubber's locked_labels use
# HOMO_1 / HET / HOMO_2; the JSON layer uses H1/H1, H1/H2, H2/H2 to match
# Quentin's manuscript-figure conventions.
KARYOTYPE_REMAP = {
    "HOMO_1": "H1/H1",
    "HET":    "H1/H2",
    "HOMO_2": "H2/H2",
    # accept either form
    "H1/H1":  "H1/H2",
    "H1/H2":  "H1/H2",
    "H2/H2":  "H2/H2",
}

# Genotype-coding for VCF GT strings. We accept the standard "/" or "|"
# separators, "." for missing.
def _gt_to_code(gt: str | None) -> str:
    """
    Map a VCF GT string to one of: 'AA' (0/0), 'AB' (0/1 or 1/0),
    'BB' (1/1), 'miss' (./. or ambiguous).
    """
    if not gt or gt in (".", "./.", ".|."):
        return "miss"
    g = gt.replace("|", "/")
    parts = g.split("/")
    if len(parts) != 2 or any(p == "." for p in parts):
        return "miss"
    try:
        a, b = int(parts[0]), int(parts[1])
    except ValueError:
        return "miss"
    if a == 0 and b == 0:
        return "AA"
    if a == 1 and b == 1:
        return "BB"
    return "AB"


def _fisher_2x2_logspace(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """
    Fisher's exact test on a 2×2 table:
        |     | car | non |
        | grp1|  a  |  b  |
        | grp2|  c  |  d  |

    Returns (odds_ratio, two_sided_p_value).

    Uses log-gamma for stability on the larger n we get with 226-sample
    cohorts. Two-sided p is the sum of all tables with probability ≤ the
    observed table's, computed exactly.
    """
    # OR = (a*d) / (b*c); guard against zero with Haldane-Anscombe correction
    if b * c == 0 and a * d == 0:
        odds_ratio = 1.0
    elif b * c == 0:
        odds_ratio = float("inf")
    else:
        odds_ratio = (a * d) / (b * c) if b != 0 and c != 0 else float("inf")

    # Hypergeometric over the row/col margins
    n = a + b + c + d
    r1 = a + b
    c1 = a + c
    if n == 0 or r1 == 0 or c1 == 0:
        return (1.0, 1.0)

    def lnchoose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))

    # logP at observed (a, b, c, d):
    def lp(x: int) -> float:
        # x = number of carriers in group1 (formerly a)
        # b' = r1 - x ; c' = c1 - x ; d' = n - r1 - c1 + x
        if x < max(0, r1 + c1 - n) or x > min(r1, c1):
            return float("-inf")
        return (lnchoose(c1, x) + lnchoose(n - c1, r1 - x) - lnchoose(n, r1))

    log_obs = lp(a)
    if log_obs == float("-inf"):
        return (odds_ratio, 1.0)

    # Sum probabilities of all tables x ∈ [max(0, r1+c1-n), min(r1, c1)]
    # whose log-probability ≤ log_obs (two-sided convention).
    lo = max(0, r1 + c1 - n)
    hi = min(r1, c1)
    log_ps = []
    for x in range(lo, hi + 1):
        lpv = lp(x)
        if lpv != float("-inf") and lpv <= log_obs + 1e-12:
            log_ps.append(lpv)
    if not log_ps:
        return (odds_ratio, 1.0)
    # logsumexp
    m = max(log_ps)
    p = math.exp(m) * sum(math.exp(v - m) for v in log_ps)
    return (odds_ratio, min(1.0, p))


def _benjamini_hochberg(pvals: list[float]) -> list[float]:
    """
    Benjamini-Hochberg FDR-adjusted p-values. Returns a list of the same
    length as `pvals`, with each entry the BH-adjusted q-value.
    """
    n = len(pvals)
    if n == 0:
        return []
    # Sort by p ascending, keeping original indices
    order = sorted(range(n), key=lambda i: pvals[i])
    adj_sorted = [0.0] * n
    cumulative_min = 1.0
    # Walk from largest p (highest rank) backwards, applying the BH formula
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        q = pvals[i] * n / rank
        if q > 1.0:
            q = 1.0
        if q < cumulative_min:
            cumulative_min = q
        adj_sorted[rank - 1] = cumulative_min
    # Map back to original index order
    adj = [0.0] * n
    for rank, i in enumerate(order):
        adj[i] = adj_sorted[rank]
    return adj


def _classify_zone(pos_bp: int, zones: dict | None,
                   bL: int, bR: int) -> str:
    """
    Decide which zone a position falls into:
      left_flank | left_boundary | inversion_body | right_boundary | right_flank

    If zone_definitions_bp is provided in the candidate, use it directly.
    Otherwise infer from boundary_left/right with ± 500 kb defaults.
    """
    if zones:
        for zname in ("left_flank", "left_boundary", "inversion_body",
                      "right_boundary", "right_flank"):
            zr = zones.get(zname)
            if zr and len(zr) == 2 and zr[0] <= pos_bp <= zr[1]:
                return zname
    # Fallback inference
    if pos_bp < bL - 500_000:                  return "left_flank"
    if pos_bp < bL:                             return "left_flank"
    if pos_bp < bL + 500_000:                   return "left_boundary"
    if pos_bp < bR - 500_000:                   return "inversion_body"
    if pos_bp < bR:                             return "right_boundary"
    if pos_bp < bR + 500_000:                   return "right_boundary"
    return "right_flank"


# Pattern-label decision rule. Spec §3.3 of SPEC_sv_evidence_page.md.
# Inputs: zone (str), gt_counts (dict of group → AA/AB/BB/miss), OR, FDR.
# Output: one of canonical_breakpoint_marker / dominant_presence_marker /
# het_specific_marker / sub_haplotype_marker / internal_linked_marker /
# uninformative.
def _classify_pattern(zone: str, gtc: dict, odds_ratio: float, fdr: float,
                      fdr_cutoff: float) -> str:
    h11 = gtc.get("H1/H1", {})
    h12 = gtc.get("H1/H2", {})
    h22 = gtc.get("H2/H2", {})
    # n in each group with non-missing genotype
    def total(g): return (g.get("AA",0) + g.get("AB",0) + g.get("BB",0))
    # carriers = AB or BB (i.e. has at least one ALT)
    def carriers(g): return (g.get("AB",0) + g.get("BB",0))
    n11, n12, n22 = total(h11), total(h12), total(h22)
    c11, c12, c22 = carriers(h11), carriers(h12), carriers(h22)
    f11 = (c11 / n11) if n11 else 0.0
    f12 = (c12 / n12) if n12 else 0.0
    f22 = (c22 / n22) if n22 else 0.0

    # Het-specific: AB rate in heterozygotes >> in either homozygote.
    # This rule is checked BEFORE the FDR gate because by construction a
    # het-specific marker is invisible to the H1/H1-vs-H2/H2 Fisher test
    # (both homozygote groups carry at the same rate, near zero), so its
    # FDR will be ~1 even though the marker is real and important.
    # Recognized by group-fraction structure instead of the OR.
    if n12 >= 5 and (h12.get("AB",0) / n12) >= 0.5 and \
       f11 <= 0.1 and f22 <= 0.1:
        return "het_specific_marker"

    # Everything else gates on FDR (pattern is meaningful only if the
    # H1/H1 vs H2/H2 contrast is statistically supported).
    if fdr >= fdr_cutoff:
        return "uninformative"

    # Boundary + strong OR + clear majority in one homozygous group:
    # canonical breakpoint marker.
    if zone in ("left_boundary", "right_boundary") and odds_ratio >= 10 \
       and ((f22 >= 0.7 and f11 <= 0.1) or (f11 >= 0.7 and f22 <= 0.1)):
        return "canonical_breakpoint_marker"

    # Dominant presence in one group, anywhere — high in one homozygote
    # but absent or near-zero in the other.
    if (f22 >= 0.6 and f11 <= 0.05) or (f11 >= 0.6 and f22 <= 0.05):
        return "dominant_presence_marker"

    # Sub-haplotype: bimodal in one homozygous group.
    for g, n in [(h11, n11), (h22, n22)]:
        if n >= 10:
            f = carriers(g) / n
            if 0.2 <= f <= 0.5:
                return "sub_haplotype_marker"

    # Internal linked: in inversion body with significant OR.
    if zone == "inversion_body" and odds_ratio >= 3:
        return "internal_linked_marker"

    return "uninformative"


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _read_karyotype(path: str | Path) -> dict[str, str]:
    """
    Read karyotype TSV (sample_id, label). Returns dict of
    sample_id → group ('H1/H1', 'H1/H2', or 'H2/H2'). Skips samples
    with unrecognized labels (logs a warning).
    """
    out = {}
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            sid, label = parts[0], parts[1]
            grp = KARYOTYPE_REMAP.get(label)
            if grp is None:
                print(f"[STEP_SV_GT_AGG] line {ln}: skipping unknown label {label!r} for {sid}",
                      file=sys.stderr)
                continue
            out[sid] = grp
    return out


def _read_calls_tsv(path: str | Path):
    """
    Read a tabular TSV of SV calls. Required columns:
      sv_id, chrom, position_bp, sv_type, sample_id, GT
    Optional: end_bp, quality, callers (semicolon-sep)

    Yields one record per SV (grouped by sv_id):
      {sv_id, chrom, position_bp, end_bp, sv_type, quality, callers,
       per_sample: dict(sample_id → 'AA'|'AB'|'BB'|'miss')}
    """
    by_sv: dict[str, dict] = {}
    header = None
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            cells = line.split("\t")
            row = dict(zip(header, cells))
            sv_id = row["sv_id"]
            rec = by_sv.get(sv_id)
            if rec is None:
                rec = {
                    "sv_id":       sv_id,
                    "chrom":       row["chrom"],
                    "position_bp": int(row["position_bp"]),
                    "end_bp":      int(row["end_bp"]) if row.get("end_bp") not in (None, "", ".") else None,
                    "sv_type":     row["sv_type"],
                    "quality":     row.get("quality", "PASS"),
                    "callers":     [c for c in (row.get("callers") or "").split(";") if c],
                    "per_sample":  {},
                }
                by_sv[sv_id] = rec
            rec["per_sample"][row["sample_id"]] = _gt_to_code(row.get("GT"))
    return list(by_sv.values())


def _read_calls_vcf(path: str | Path):
    """
    Minimal VCF reader for the SV-calls case. Looks at FORMAT/GT only;
    builds the same per-SV records as _read_calls_tsv. Designed for
    single-sample-per-line OR multi-sample VCF with a fixed sample-name
    header. This is a deliberately lightweight reader (no external deps);
    for production VCFs, swap in pysam VariantFile.

    Yields per-SV dicts; sv_id is built from CHROM:POS:SVTYPE if no ID.
    """
    samples: list[str] = []
    by_sv: dict[str, dict] = {}
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            if line.startswith("##"):
                continue
            line = line.rstrip("\n")
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                samples = fields[9:]   # everything after FORMAT
                continue
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, _ref, _alt, _qual, fltr, info, fmt = parts[:9]
            gts = parts[9:]
            # Pull SVTYPE / END from INFO
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
            sv_type = info_d.get("SVTYPE", "BND")
            end_bp = int(info_d["END"]) if "END" in info_d else None
            sv_id = vid if (vid and vid != ".") else f"{chrom}:{pos}:{sv_type}"
            # Find GT index in FORMAT
            fmt_keys = fmt.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                continue
            rec = {
                "sv_id":       sv_id,
                "chrom":       chrom,
                "position_bp": int(pos),
                "end_bp":      end_bp,
                "sv_type":     sv_type,
                "quality":     "PASS" if fltr in ("PASS", ".") else fltr,
                "callers":     [],   # populate from INFO/SOURCE if present
                "per_sample":  {},
            }
            for sname, sval in zip(samples, gts):
                if sval == ".":
                    rec["per_sample"][sname] = "miss"
                else:
                    sub = sval.split(":")
                    rec["per_sample"][sname] = _gt_to_code(sub[gt_idx])
            by_sv[sv_id] = rec
    return list(by_sv.values())


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------

def aggregate(calls: list[dict], karyotype: dict[str, str], cand: dict,
              fdr_cutoff: float = 0.05) -> dict:
    """
    Build the sv_genotype_counts_v1 payload from raw SV calls + karyotype
    labels + candidate metadata. Returns a dict ready to be written via
    write_candidate_layer().
    """
    bL = cand["boundary_left_bp"]
    bR = cand["boundary_right_bp"]
    zones = cand.get("zone_definitions_bp")

    # Group counts (n per group)
    grp_n = defaultdict(int)
    for sid, grp in karyotype.items():
        grp_n[grp] += 1

    # First pass: per-SV genotype counts + raw Fisher p
    raw = []
    for c in calls:
        # Filter to SV calls within candidate window (with flanking)
        # — accept anything within [bL - 1Mb, bR + 1Mb]
        if not (bL - 1_000_000 <= c["position_bp"] <= bR + 1_000_000):
            continue
        gtc = {"H1/H1": {"AA":0,"AB":0,"BB":0,"miss":0},
               "H1/H2": {"AA":0,"AB":0,"BB":0,"miss":0},
               "H2/H2": {"AA":0,"AB":0,"BB":0,"miss":0}}
        for sid, code in c["per_sample"].items():
            grp = karyotype.get(sid)
            if not grp:
                continue
            gtc[grp][code] = gtc[grp].get(code, 0) + 1

        # Fisher 2x2: H1/H1 carriers vs H2/H2 carriers
        h11_carrier = gtc["H1/H1"]["AB"] + gtc["H1/H1"]["BB"]
        h11_non     = gtc["H1/H1"]["AA"]
        h22_carrier = gtc["H2/H2"]["AB"] + gtc["H2/H2"]["BB"]
        h22_non     = gtc["H2/H2"]["AA"]
        odds, p = _fisher_2x2_logspace(h11_carrier, h11_non, h22_carrier, h22_non)

        n_with_call = sum(1 for v in c["per_sample"].values()
                          if v in ("AA","AB","BB"))
        zone = _classify_zone(c["position_bp"], zones, bL, bR)
        dist_to_edge = (c["position_bp"] - bL) if abs(c["position_bp"] - bL) < abs(c["position_bp"] - bR) \
                       else (c["position_bp"] - bR)

        raw.append({
            "_call":     c,
            "gtc":       gtc,
            "odds":      odds,
            "p":         p,
            "n_with":    n_with_call,
            "zone":      zone,
            "dist":      dist_to_edge,
        })

    # FDR-BH adjustment over all SVs in the candidate window
    fdrs = _benjamini_hochberg([r["p"] for r in raw])

    sv_calls = []
    for r, fdr in zip(raw, fdrs):
        c = r["_call"]
        pattern = _classify_pattern(r["zone"], r["gtc"], r["odds"], fdr, fdr_cutoff)
        sv_calls.append({
            "sv_id":              c["sv_id"],
            "sv_type":            c["sv_type"],
            "chrom":              c["chrom"],
            "position_bp":        c["position_bp"],
            "end_bp":             c["end_bp"],
            "zone":               r["zone"],
            "distance_to_edge_bp": r["dist"],
            "n_samples_with_call": r["n_with"],
            "quality":            c["quality"],
            "callers":            c["callers"],
            "genotype_counts":    r["gtc"],
            "fisher": {
                "comparison":  "H1/H1_vs_H2/H2",
                "odds_ratio":  None if r["odds"] == float("inf") else r["odds"],
                "p_value":     r["p"],
                "fdr_bh":      fdr,
            },
            "pattern_label":      pattern,
            "notes":              "",
        })

    # Boundary summary (per-side, per-SV-type counts)
    def _summary_for_zone(zone_keys: list[str]) -> dict:
        out = {"BND":  {"n_total":0,"n_associated_fdr_lt_0_05":0},
               "INV":  {"n_total":0,"n_associated_fdr_lt_0_05":0},
               "DEL":  {"n_total":0,"n_associated_fdr_lt_0_05":0},
               "DUP":  {"n_total":0,"n_associated_fdr_lt_0_05":0},
               "Other":{"n_total":0,"n_associated_fdr_lt_0_05":0}}
        for s in sv_calls:
            if s["zone"] not in zone_keys:
                continue
            t = s["sv_type"] if s["sv_type"] in out else "Other"
            out[t]["n_total"] += 1
            if s["fisher"]["fdr_bh"] is not None and s["fisher"]["fdr_bh"] < fdr_cutoff:
                out[t]["n_associated_fdr_lt_0_05"] += 1
        return out

    interval_left  = (zones or {}).get("left_boundary",  [bL - 500_000, bL + 500_000])
    interval_right = (zones or {}).get("right_boundary", [bR - 500_000, bR + 500_000])

    payload = {
        "format_version":      "sv_genotype_counts_v1",
        "candidate_id":        cand["candidate_id"],
        "chrom":               cand["chrom"],
        "boundary_left_bp":    bL,
        "boundary_right_bp":   bR,
        "zone_definitions_bp": zones or {},
        "groups_used": {
            grp: {"n": grp_n[grp], "members": []}
            for grp in ("H1/H1","H1/H2","H2/H2")
        },
        "sv_calls":            sv_calls,
        "boundary_summary": {
            "left":  {"interval_bp": interval_left,
                      "by_sv_type":  _summary_for_zone(["left_boundary"])},
            "right": {"interval_bp": interval_right,
                      "by_sv_type":  _summary_for_zone(["right_boundary"])},
        },
        # upset_top_combinations is now produced by STEP_SV_EVID_COMB; leave
        # empty here so the atlas falls back gracefully if only this layer
        # is loaded.
        "upset_top_combinations": [],
    }
    return payload


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", help="merged VCF with GT field")
    ap.add_argument("--tsv", help="tabular TSV (alternative to --vcf)")
    ap.add_argument("--candidate", required=True, help="candidate JSON")
    ap.add_argument("--karyotype", required=True, help="karyotype TSV")
    ap.add_argument("--out-root", required=True, help="output root")
    ap.add_argument("--fdr-cutoff", type=float, default=0.05)
    ap.add_argument("--indent", type=int, default=None,
                    help="JSON indent (default: compact)")
    ap.add_argument("--no-manifest", action="store_true",
                    help="skip rebuilding the chrom-level manifest.json")
    args = ap.parse_args()

    if not args.vcf and not args.tsv:
        ap.error("must give --vcf or --tsv")

    cand = json.loads(Path(args.candidate).read_text())
    karyotype = _read_karyotype(args.karyotype)
    print(f"[STEP_SV_GT_AGG] {len(karyotype)} karyotyped samples", file=sys.stderr)

    if args.vcf:
        calls = _read_calls_vcf(args.vcf)
    else:
        calls = _read_calls_tsv(args.tsv)
    print(f"[STEP_SV_GT_AGG] {len(calls)} SV records read", file=sys.stderr)

    payload = aggregate(calls, karyotype, cand, fdr_cutoff=args.fdr_cutoff)
    print(f"[STEP_SV_GT_AGG] {len(payload['sv_calls'])} SVs in window after filtering",
          file=sys.stderr)

    out = write_candidate_layer(payload, args.out_root, cand["chrom"],
                                cand["candidate_id"], indent=args.indent)
    log_write(out)

    if not args.no_manifest:
        m = update_chrom_manifest(args.out_root, cand["chrom"])
        log_write(m)


if __name__ == "__main__":
    main()
