#!/usr/bin/env python3
"""
STEP01 – Parse one Clair3 discovery VCF and produce a master annotated table.

Outputs:
  all_variants_annotated.tsv   – one row per VCF record with all parsed fields
                                 plus local neighborhood and overlap annotations.

Usage:
  python STEP01_parse_and_annotate_clair3_vcf.py \
      --vcf  <input.vcf.gz> \
      --outdir <output_directory> \
      --sample SAMPLE_NAME \
      [--ref_fai  ref.fa.fai] \
      [--bed  calling_region.bed]
"""

import os, sys, argparse, math
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam


# ── helpers ──────────────────────────────────────────────────────────────────

def classify_variant(ref, alt):
    if alt is None or alt == ".":
        return "OTHER"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(alt) > len(ref):
        return "INS"
    if len(alt) < len(ref):
        return "DEL"
    return "OTHER"

def classify_gt(gt_str):
    if gt_str is None or gt_str == ".":
        return "OTHER"
    gt = gt_str.replace("|", "/")
    alleles = gt.split("/")
    if "." in alleles:
        return "OTHER"
    uniq = set(alleles)
    if len(uniq) == 1:
        return "HOM_REF" if "0" in uniq else "HOM_VAR"
    return "HET"

def parse_info_source(rec):
    has_p = "P" in rec.info
    has_f = "F" in rec.info
    if has_p and has_f:
        return "PF"
    if has_p:
        return "P"
    if has_f:
        return "F"
    return "NONE"

def safe_field(sample, key):
    try:
        return sample[key]
    except Exception:
        return None

def read_bed_span(bed_path):
    if bed_path is None:
        return {}
    spans = defaultdict(int)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")[:3]
            spans[parts[0]] += max(0, int(parts[2]) - int(parts[1]))
    return dict(spans)

def read_fai_lengths(fai_path):
    if fai_path is None:
        return {}
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            lengths[f[0]] = int(f[1])
    return lengths


# ── main parser ──────────────────────────────────────────────────────────────

def extract_records(vcf_path, sample_name):
    """Parse every record into a dict list."""
    vf = pysam.VariantFile(vcf_path)
    samples = list(vf.header.samples)
    if sample_name not in samples:
        # try first sample
        if len(samples) == 1:
            sample_name = samples[0]
            print(f"[WARN] Requested sample not found, using sole sample: {sample_name}", file=sys.stderr)
        else:
            raise ValueError(f"Sample '{sample_name}' not in VCF. Available: {samples}")

    rows = []
    for rec in vf.fetch():
        # FILTER
        filt_keys = list(rec.filter.keys())
        if len(filt_keys) == 0:
            filt, status = "PASS", "PASS"
        else:
            filt = ";".join(filt_keys)
            status = "PASS" if filt == "PASS" else "FILTERED"

        s = rec.samples[sample_name]

        # GT
        gt_tuple = safe_field(s, "GT")
        gt = "/".join("." if x is None else str(x) for x in gt_tuple) if gt_tuple else "."

        gq = safe_field(s, "GQ")
        dp = safe_field(s, "DP")
        ad = safe_field(s, "AD")
        af = safe_field(s, "AF")
        pl = safe_field(s, "PL")

        alts = rec.alts if rec.alts else []
        alt1 = alts[0] if len(alts) >= 1 else "."
        var_type = classify_variant(rec.ref, alt1)

        ad_ref = ad_alt = ad_alt2 = ad_sum = np.nan
        if ad is not None:
            ad = list(ad)
            if len(ad) >= 1 and ad[0] is not None: ad_ref = float(ad[0])
            if len(ad) >= 2 and ad[1] is not None: ad_alt = float(ad[1])
            if len(ad) >= 3 and ad[2] is not None: ad_alt2 = float(ad[2])
            valid = [x for x in ad if x is not None]
            if valid: ad_sum = float(sum(valid))

        af1 = af2 = np.nan
        if af is not None:
            af = list(af) if isinstance(af, tuple) else ([af] if not isinstance(af, list) else af)
            if len(af) >= 1 and af[0] is not None: af1 = float(af[0])
            if len(af) >= 2 and af[1] is not None: af2 = float(af[1])

        indel_len = np.nan
        if var_type in ("INS", "DEL"):
            indel_len = len(alt1) - len(rec.ref)

        rows.append({
            "CHROM": rec.chrom,
            "POS": int(rec.pos),
            "END": int(rec.stop),
            "REF": rec.ref,
            "ALT1": alt1,
            "ALT_ALL": ",".join(alts) if alts else ".",
            "N_ALT": len(alts),
            "QUAL": float(rec.qual) if rec.qual is not None else np.nan,
            "FILTER": filt,
            "STATUS": status,
            "CALLER_SOURCE": parse_info_source(rec),
            "VAR_TYPE": var_type,
            "IS_SNP": 1 if var_type == "SNP" else 0,
            "IS_INDEL": 1 if var_type in ("INS", "DEL") else 0,
            "IS_SIMPLE_BIALLELIC": 1 if len(alts) == 1 else 0,
            "IS_MULTIALLELIC": 1 if len(alts) > 1 else 0,
            "INDEL_LEN": indel_len,
            "ABS_INDEL_LEN": abs(indel_len) if pd.notnull(indel_len) else np.nan,
            "GT": gt,
            "GT_CLASS": classify_gt(gt),
            "GQ": float(gq) if gq is not None else np.nan,
            "DP": float(dp) if dp is not None else np.nan,
            "AD_REF": ad_ref,
            "AD_ALT": ad_alt,
            "AD_ALT2": ad_alt2,
            "AD_SUM": ad_sum,
            "AF1": af1,
            "AF2": af2,
            "PL": ",".join(map(str, pl)) if pl is not None else "",
        })

    return pd.DataFrame(rows)


# ── local neighborhood annotation ───────────────────────────────────────────

def annotate_local_neighborhood(df):
    """Add distance-to-neighbour and nearby-variant flags per chromosome."""
    if df.empty:
        for c in ("DIST_PREV", "DIST_NEXT",
                   "HAS_NEARBY_VAR_1BP", "HAS_NEARBY_VAR_5BP",
                   "HAS_NEARBY_VAR_10BP", "HAS_NEARBY_VAR_20BP",
                   "N_LOCAL_VAR_10BP", "N_LOCAL_VAR_20BP"):
            df[c] = np.nan
        return df

    out_parts = []
    for chrom, grp in df.groupby("CHROM", sort=False):
        g = grp.sort_values("POS").copy()
        positions = g["POS"].values

        dist_prev = np.full(len(positions), np.nan)
        dist_next = np.full(len(positions), np.nan)
        if len(positions) > 1:
            dist_prev[1:] = positions[1:] - positions[:-1]
            dist_next[:-1] = positions[1:] - positions[:-1]

        g["DIST_PREV"] = dist_prev
        g["DIST_NEXT"] = dist_next

        # nearby flags using vectorised approach
        for win, col_flag, col_count in [
            (1,  "HAS_NEARBY_VAR_1BP",  None),
            (5,  "HAS_NEARBY_VAR_5BP",  None),
            (10, "HAS_NEARBY_VAR_10BP", "N_LOCAL_VAR_10BP"),
            (20, "HAS_NEARBY_VAR_20BP", "N_LOCAL_VAR_20BP"),
        ]:
            flags = np.zeros(len(positions), dtype=int)
            counts = np.zeros(len(positions), dtype=int)
            for i in range(len(positions)):
                p = positions[i]
                lo = np.searchsorted(positions, p - win, side='left')
                hi = np.searchsorted(positions, p + win, side='right')
                n_nearby = (hi - lo) - 1  # exclude self
                flags[i] = 1 if n_nearby > 0 else 0
                counts[i] = n_nearby
            g[col_flag] = flags
            if col_count:
                g[col_count] = counts

        out_parts.append(g)

    return pd.concat(out_parts, ignore_index=True)


# ── overlap / cluster detection ──────────────────────────────────────────────

def annotate_overlaps(df):
    """Detect overlapping records and assign overlap cluster IDs."""
    df["IS_OVERLAPPING"] = 0
    df["OVERLAP_CLUSTER_ID"] = -1

    if df.empty:
        return df

    cluster_id = 0
    for chrom, grp in df.groupby("CHROM", sort=False):
        idx = grp.sort_values("POS").index.tolist()
        if len(idx) == 0:
            continue

        i = 0
        while i < len(idx):
            # start a potential cluster
            cluster_members = [idx[i]]
            cluster_end = df.loc[idx[i], "END"]
            j = i + 1
            while j < len(idx):
                if df.loc[idx[j], "POS"] <= cluster_end:
                    cluster_members.append(idx[j])
                    cluster_end = max(cluster_end, df.loc[idx[j], "END"])
                    j += 1
                else:
                    break

            if len(cluster_members) > 1:
                for m in cluster_members:
                    df.loc[m, "IS_OVERLAPPING"] = 1
                    df.loc[m, "OVERLAP_CLUSTER_ID"] = cluster_id
                cluster_id += 1

            i = j

    return df


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="STEP01: Parse & annotate Clair3 VCF")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--bed", default=None)
    ap.add_argument("--ref_fai", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[STEP01] Parsing VCF: {args.vcf}")
    df = extract_records(args.vcf, args.sample)
    print(f"[STEP01] Records parsed: {len(df)}")

    print("[STEP01] Annotating local neighborhood …")
    df = annotate_local_neighborhood(df)

    print("[STEP01] Detecting overlapping records …")
    df = annotate_overlaps(df)

    out_path = os.path.join(args.outdir, "all_variants_annotated.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[STEP01] Written: {out_path}  ({len(df)} rows)")


if __name__ == "__main__":
    main()
