#!/usr/bin/env python3
"""
STEP08 – Final variant classification and summary tables.

Merges all classification layers:
  1. STRICT_PASS                  – original Clair3 PASS
  2. RESCUED_STRONG_SINGLE_SAMPLE – from STEP03 rescue
  3. RESCUED_POPULATION_SIGNATURE – weak indels confirmed by population regenotyping
  4. WEAK_SIGNATURE_CANDIDATE     – unconfirmed weak candidates
  5. FILTERED_DISCARD             – filtered variants with no rescue
  6. COMPLEX_LOCAL_BLOCK_ONLY     – variants in complex blocks (informational)

Broader layers:
  strict_marker / relaxed_marker / rescued_population_marker /
  weak_signature_candidate / complex_block_only / discard

Outputs:
  final_variant_classification.tsv
  summary_counts_by_final_class.tsv
  summary_counts_by_final_layer.tsv
  summary_counts_by_type_status_gt.tsv
  summary_local_cluster_stats.tsv
  summary_rescue_stats.tsv
  summary_metric_by_class.tsv

Usage:
  python STEP08_final_classification.py \
      --step03     all_classified_step03.tsv \
      --weak       weak_indel_candidates.tsv \
      --regen      regenotyped_rescued_indels.tsv  (optional) \
      --regen_summary regenotype_cohort_summary.tsv (optional) \
      --blocks     local_overlap_blocks.tsv \
      --outdir     <output_directory> \
      --sample_id  SAMPLE_ID \
      [--min_regen_alt_samples 2]
"""

import os, argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="STEP08: Final classification")
    ap.add_argument("--step03", required=True, help="all_classified_step03.tsv")
    ap.add_argument("--weak", default=None, help="weak_indel_candidates.tsv")
    ap.add_argument("--regen", default=None, help="regenotyped_rescued_indels.tsv (optional)")
    ap.add_argument("--regen_summary", default=None, help="regenotype_cohort_summary.tsv (optional)")
    ap.add_argument("--blocks", default=None, help="local_overlap_blocks.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample_id", default="SAMPLE")
    ap.add_argument("--min_regen_alt_samples", type=int, default=2,
                    help="Min samples with alt call for population rescue confirmation [2]")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load base classification
    df = pd.read_csv(args.step03, sep="\t")
    print(f"[STEP08] Base records: {len(df)}")

    # ── integrate population regenotyping results ────────────────────────────
    # If regenotype results are available, upgrade PENDING weak indels
    # that were confirmed across the population
    regen_confirmed = set()
    if args.regen_summary and os.path.isfile(args.regen_summary):
        regen_sum = pd.read_csv(args.regen_summary, sep="\t")
        confirmed = regen_sum[regen_sum["N_CALLED_ALT"] >= args.min_regen_alt_samples]
        for _, row in confirmed.iterrows():
            key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
            regen_confirmed.add(key)
        print(f"[STEP08] Population-confirmed candidates: {len(regen_confirmed)}")

    # Load weak candidates to cross-reference
    weak_keys = set()
    if args.weak and os.path.isfile(args.weak):
        weak_df = pd.read_csv(args.weak, sep="\t")
        if len(weak_df) > 0:
            for _, row in weak_df.iterrows():
                key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
                weak_keys.add(key)

    # ── final classification ─────────────────────────────────────────────────
    def final_classify(row):
        current = row.get("FINAL_CLASS", "PENDING")
        if current == "STRICT_PASS":
            return "STRICT_PASS", "strict_marker", "pass"
        if current == "RESCUED_STRONG_SINGLE_SAMPLE":
            return current, "relaxed_marker", row.get("RESCUE_REASON", "strong_single_sample")

        # check population rescue for PENDING indels
        if current == "PENDING" and row.get("IS_INDEL", 0) == 1:
            key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT1']}"
            if key in regen_confirmed:
                return "RESCUED_POPULATION_SIGNATURE", "rescued_population_marker", "population_regenotype"
            if key in weak_keys:
                return "WEAK_SIGNATURE_CANDIDATE", "weak_signature_candidate", "weak_candidate"

        # complex block check
        block_id = row.get("BLOCK_ID", -1)
        if block_id >= 0 and row.get("IS_OVERLAPPING", 0) == 1:
            if current == "PENDING":
                return "COMPLEX_LOCAL_BLOCK_ONLY", "complex_block_only", "complex_block"

        if current == "PENDING":
            return "FILTERED_DISCARD", "discard", "none"

        return current, row.get("FINAL_LAYER", "discard"), row.get("RESCUE_REASON", "none")

    results = df.apply(final_classify, axis=1, result_type="expand")
    results.columns = ["FINAL_CLASS_F", "FINAL_LAYER_F", "RESCUE_REASON_F"]
    df["FINAL_CLASS"] = results["FINAL_CLASS_F"]
    df["FINAL_LAYER"] = results["FINAL_LAYER_F"]
    df["RESCUE_REASON"] = results["RESCUE_REASON_F"]

    # output main table
    out_main = os.path.join(args.outdir, "final_variant_classification.tsv")
    df.to_csv(out_main, sep="\t", index=False)
    print(f"[STEP08] Final classification → {out_main}")

    # ── summary tables ───────────────────────────────────────────────────────

    # 1. counts by final class
    tab1 = df.groupby("FINAL_CLASS").size().reset_index(name="N")
    tab1.to_csv(os.path.join(args.outdir, "summary_counts_by_final_class.tsv"),
                sep="\t", index=False)

    # 2. counts by final layer
    tab2 = df.groupby("FINAL_LAYER").size().reset_index(name="N")
    tab2.to_csv(os.path.join(args.outdir, "summary_counts_by_final_layer.tsv"),
                sep="\t", index=False)

    # 3. counts by type × status × GT
    tab3 = (df.groupby(["VAR_TYPE", "STATUS", "GT_CLASS", "FINAL_CLASS"])
              .size().reset_index(name="N"))
    tab3.to_csv(os.path.join(args.outdir, "summary_counts_by_type_status_gt.tsv"),
                sep="\t", index=False)

    # 4. rescue stats
    rescue_tab = (df[df["FINAL_CLASS"].str.startswith("RESCUED")]
                    .groupby(["FINAL_CLASS", "RESCUE_REASON", "VAR_TYPE"])
                    .size().reset_index(name="N"))
    rescue_tab.to_csv(os.path.join(args.outdir, "summary_rescue_stats.tsv"),
                      sep="\t", index=False)

    # 5. local cluster stats
    if args.blocks and os.path.isfile(args.blocks):
        blocks = pd.read_csv(args.blocks, sep="\t")
        blocks.to_csv(os.path.join(args.outdir, "summary_local_cluster_stats.tsv"),
                      sep="\t", index=False)

    # 6. metric summaries by class
    metrics = ["QUAL", "GQ", "DP", "AF1", "AD_ALT", "ABS_INDEL_LEN"]
    stat_rows = []
    for metric in metrics:
        if metric not in df.columns:
            continue
        for fc in df["FINAL_CLASS"].dropna().unique():
            for vt in df["VAR_TYPE"].dropna().unique():
                for gt in df["GT_CLASS"].dropna().unique():
                    vals = df.loc[
                        (df["FINAL_CLASS"] == fc) &
                        (df["VAR_TYPE"] == vt) &
                        (df["GT_CLASS"] == gt),
                        metric
                    ].dropna()
                    if len(vals) == 0:
                        continue
                    stat_rows.append({
                        "METRIC": metric,
                        "FINAL_CLASS": fc,
                        "VAR_TYPE": vt,
                        "GT_CLASS": gt,
                        "N": len(vals),
                        "MEAN": round(float(vals.mean()), 4),
                        "MEDIAN": round(float(vals.median()), 4),
                        "MIN": float(vals.min()),
                        "MAX": float(vals.max()),
                        "Q05": round(float(vals.quantile(0.05)), 4),
                        "Q25": round(float(vals.quantile(0.25)), 4),
                        "Q75": round(float(vals.quantile(0.75)), 4),
                        "Q95": round(float(vals.quantile(0.95)), 4),
                    })

    pd.DataFrame(stat_rows).to_csv(
        os.path.join(args.outdir, "summary_metric_by_class.tsv"),
        sep="\t", index=False)

    # ── print summary ────────────────────────────────────────────────────────
    print("\n[STEP08] ═══════ FINAL CLASSIFICATION SUMMARY ═══════")
    for _, row in tab1.iterrows():
        print(f"  {row['FINAL_CLASS']:40s}  {row['N']:>8d}")
    print(f"  {'TOTAL':40s}  {len(df):>8d}")
    print("[STEP08] Done.")


if __name__ == "__main__":
    main()
