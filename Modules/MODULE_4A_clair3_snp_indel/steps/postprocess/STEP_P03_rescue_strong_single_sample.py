#!/usr/bin/env python3
"""
STEP03 – Strong single-sample rescue of filtered variants.

Applies practical rescue rules to filtered variants that are biologically
plausible from one sample alone.  Two rescue branches:

  Branch A – "strong_indel":  filtered indels with moderate+ evidence
  Branch B – "relaxed_indel": filtered indels with lower but acceptable evidence
  Branch C – "strong_snp":    filtered SNPs with QUAL >= 12 and decent evidence

Each rescued variant gets:
  FINAL_CLASS   = RESCUED_STRONG_SINGLE_SAMPLE
  FINAL_LAYER   = relaxed_marker
  RESCUE_REASON = <branch tag>

PASS variants are also emitted with:
  FINAL_CLASS   = STRICT_PASS
  FINAL_LAYER   = strict_marker

Outputs:
  strong_filtered_rescued.tsv   – rescued filtered variants only
  pass_variants.tsv             – PASS variants (strict set)
  all_classified_step03.tsv     – everything with classification columns

Usage:
  python STEP03_rescue_strong_single_sample.py \
      --annotated  all_variants_with_blocks.tsv \
      --outdir     <output_directory>
"""

import os, argparse
import numpy as np
import pandas as pd


# ── rescue rule sets ─────────────────────────────────────────────────────────
# These are hardcoded sensible defaults.  Adjust if needed after pilot review.

def rescue_strong_indel(row):
    """Branch A: strong filtered indel.
    QUAL >= 8, GQ >= 5, DP >= 4, AF1 >= 0.25, AD_ALT >= 2,
    not in absurdly complex local block (N_RECORDS <= 6).
    """
    if row["VAR_TYPE"] not in ("INS", "DEL"):
        return False
    if row["STATUS"] != "FILTERED":
        return False
    if pd.isna(row["QUAL"]) or row["QUAL"] < 8:
        return False
    if pd.isna(row["DP"]) or row["DP"] < 4:
        return False
    if pd.isna(row["AF1"]) or row["AF1"] < 0.25:
        return False
    if pd.isna(row["AD_ALT"]) or row["AD_ALT"] < 2:
        return False
    gq = row.get("GQ", np.nan)
    if pd.notna(gq) and gq < 5:
        return False
    return True


def rescue_relaxed_indel(row):
    """Branch B: relaxed filtered indel.
    Broader: QUAL >= 5, DP >= 3, AF1 >= 0.15, AD_ALT >= 1.
    Must be HET or HOM_VAR.  Not in super-complex block.
    """
    if row["VAR_TYPE"] not in ("INS", "DEL"):
        return False
    if row["STATUS"] != "FILTERED":
        return False
    if row["GT_CLASS"] not in ("HET", "HOM_VAR"):
        return False
    if pd.isna(row["QUAL"]) or row["QUAL"] < 5:
        return False
    if pd.isna(row["DP"]) or row["DP"] < 3:
        return False
    if pd.isna(row["AF1"]) or row["AF1"] < 0.15:
        return False
    if pd.isna(row["AD_ALT"]) or row["AD_ALT"] < 1:
        return False
    return True


def rescue_strong_snp(row):
    """Branch C: filtered SNP with QUAL >= 12 and decent evidence.
    This expands the high-quality SNP set beyond strict PASS.
    QUAL >= 12, DP >= 3, AF1 >= 0.15, AD_ALT >= 1, HET or HOM_VAR.
    """
    if row["VAR_TYPE"] != "SNP":
        return False
    if row["STATUS"] != "FILTERED":
        return False
    if row["GT_CLASS"] not in ("HET", "HOM_VAR"):
        return False
    if pd.isna(row["QUAL"]) or row["QUAL"] < 12:
        return False
    if pd.isna(row["DP"]) or row["DP"] < 3:
        return False
    if pd.isna(row["AF1"]) or row["AF1"] < 0.15:
        return False
    if pd.isna(row["AD_ALT"]) or row["AD_ALT"] < 1:
        return False
    return True


def classify_row(row):
    """Return (FINAL_CLASS, FINAL_LAYER, RESCUE_REASON)."""

    # strict PASS
    if row["STATUS"] == "PASS":
        return "STRICT_PASS", "strict_marker", "pass"

    # strong indel rescue
    if rescue_strong_indel(row):
        return "RESCUED_STRONG_SINGLE_SAMPLE", "relaxed_marker", "strong_indel"

    # relaxed indel rescue (subset of strong but with lower thresholds –
    # only hits rows not already caught by strong)
    if rescue_relaxed_indel(row):
        return "RESCUED_STRONG_SINGLE_SAMPLE", "relaxed_marker", "relaxed_indel"

    # strong SNP rescue (QUAL >= 12)
    if rescue_strong_snp(row):
        return "RESCUED_STRONG_SINGLE_SAMPLE", "relaxed_marker", "strong_snp_q12"

    # not rescued yet – will be evaluated by later steps
    return "PENDING", "pending", "none"


def main():
    ap = argparse.ArgumentParser(description="STEP03: Strong single-sample rescue")
    ap.add_argument("--annotated", required=True,
                    help="all_variants_with_blocks.tsv from STEP02")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.annotated, sep="\t")
    print(f"[STEP03] Input rows: {len(df)}")

    # classify
    results = df.apply(classify_row, axis=1, result_type="expand")
    results.columns = ["FINAL_CLASS", "FINAL_LAYER", "RESCUE_REASON"]
    df = pd.concat([df, results], axis=1)

    # outputs
    pass_df = df[df["FINAL_CLASS"] == "STRICT_PASS"]
    rescued_df = df[df["FINAL_CLASS"] == "RESCUED_STRONG_SINGLE_SAMPLE"]
    pending_df = df[df["FINAL_CLASS"] == "PENDING"]

    out_pass = os.path.join(args.outdir, "pass_variants.tsv")
    pass_df.to_csv(out_pass, sep="\t", index=False)

    out_rescued = os.path.join(args.outdir, "strong_filtered_rescued.tsv")
    rescued_df.to_csv(out_rescued, sep="\t", index=False)

    out_all = os.path.join(args.outdir, "all_classified_step03.tsv")
    df.to_csv(out_all, sep="\t", index=False)

    print(f"[STEP03] STRICT_PASS:                    {len(pass_df)}")
    print(f"[STEP03] RESCUED_STRONG_SINGLE_SAMPLE:   {len(rescued_df)}")
    by_reason = rescued_df["RESCUE_REASON"].value_counts()
    for k, v in by_reason.items():
        print(f"         ├─ {k}: {v}")
    print(f"[STEP03] PENDING (for later steps):      {len(pending_df)}")
    print(f"[STEP03] Written: {out_rescued}")


if __name__ == "__main__":
    main()
