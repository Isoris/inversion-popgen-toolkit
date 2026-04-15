#!/usr/bin/env python3
"""
STEP07B – Merge per-sample regenotyping results + cohort summary (streaming)

This version avoids loading all per-sample TSVs into memory at once.

Usage:
  python STEP07B_merge_regenotype.py \
      --indir   <population_output_dir> \
      --outdir  <population_output_dir>
"""

import os
import glob
import time
import argparse
from collections import defaultdict

import pandas as pd


def _now():
    return time.time()


def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"


def safe_int(x, default=0):
    try:
        if pd.isna(x):
            return default
        return int(x)
    except Exception:
        return default


def safe_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def update_summary(summary, row):
    regen_id = row["REGEN_ID"]

    if regen_id not in summary:
        summary[regen_id] = {
            "REGEN_ID": regen_id,
            "CHROM": row.get("CHROM", ""),
            "POS": row.get("POS", ""),
            "REF": row.get("REF", ""),
            "ALT": row.get("ALT", ""),
            "VAR_TYPE": row.get("VAR_TYPE", ""),
            "N_SAMPLES_TOTAL": 0,
            "N_CALLED_ALT": 0,
            "N_HET": 0,
            "N_HOM_ALT": 0,
            "N_HOM_REF": 0,
            "N_NO_DATA": 0,
            "SUM_SUPPORT_FRACTION": 0.0,
            "TOTAL_READS_SUPPORT": 0,
        }

    rec = summary[regen_id]
    gt = str(row.get("REGEN_GT", "./."))
    sf = safe_float(row.get("SUPPORT_FRACTION", 0.0), 0.0)
    nrs = safe_int(row.get("N_READS_SUPPORT", 0), 0)

    rec["N_SAMPLES_TOTAL"] += 1
    rec["SUM_SUPPORT_FRACTION"] += sf
    rec["TOTAL_READS_SUPPORT"] += nrs

    if gt in {"0/1", "1/1"}:
        rec["N_CALLED_ALT"] += 1
    if gt == "0/1":
        rec["N_HET"] += 1
    elif gt == "1/1":
        rec["N_HOM_ALT"] += 1
    elif gt == "0/0":
        rec["N_HOM_REF"] += 1
    elif gt == "./.":
        rec["N_NO_DATA"] += 1


def finalize_summary(summary_dict):
    rows = []
    for _, rec in summary_dict.items():
        n_total = rec["N_SAMPLES_TOTAL"]
        mean_sf = rec["SUM_SUPPORT_FRACTION"] / n_total if n_total > 0 else 0.0
        alt_freq = rec["N_CALLED_ALT"] / n_total if n_total > 0 else 0.0

        rows.append({
            "REGEN_ID": rec["REGEN_ID"],
            "CHROM": rec["CHROM"],
            "POS": rec["POS"],
            "REF": rec["REF"],
            "ALT": rec["ALT"],
            "VAR_TYPE": rec["VAR_TYPE"],
            "N_SAMPLES_TOTAL": rec["N_SAMPLES_TOTAL"],
            "N_CALLED_ALT": rec["N_CALLED_ALT"],
            "N_HET": rec["N_HET"],
            "N_HOM_ALT": rec["N_HOM_ALT"],
            "N_HOM_REF": rec["N_HOM_REF"],
            "N_NO_DATA": rec["N_NO_DATA"],
            "ALT_FREQ": round(alt_freq, 4),
            "MEAN_SUPPORT_FRACTION": round(mean_sf, 4),
            "TOTAL_READS_SUPPORT": int(rec["TOTAL_READS_SUPPORT"]),
        })

    if not rows:
        return pd.DataFrame()

    summary_df = pd.DataFrame(rows)

    # stable sort if columns exist
    sort_cols = [c for c in ["CHROM", "POS", "REGEN_ID"] if c in summary_df.columns]
    if sort_cols:
        summary_df = summary_df.sort_values(sort_cols, kind="stable").reset_index(drop=True)

    return summary_df


def main():
    t0 = _now()

    ap = argparse.ArgumentParser(description="STEP07B: Merge per-sample regenotype results (streaming)")
    ap.add_argument("--indir", required=True, help="Directory with per-sample TSVs")
    ap.add_argument("--outdir", required=True, help="Output directory (can be same as indir)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pattern = os.path.join(args.indir, "regenotyped_rescued_indels.*.tsv")
    files = sorted(glob.glob(pattern))
    files = [f for f in files if not f.endswith("regenotyped_rescued_indels.tsv")]

    print(f"[STEP07B] Found {len(files)} per-sample regenotype files")

    if len(files) == 0:
        print("[STEP07B] No per-sample files found. Nothing to merge.")
        return

    out_detail = os.path.join(args.outdir, "regenotyped_rescued_indels.tsv")
    out_summary = os.path.join(args.outdir, "regenotype_cohort_summary.tsv")

    # remove previous outputs to avoid accidental appends across reruns
    if os.path.exists(out_detail):
        os.remove(out_detail)
    if os.path.exists(out_summary):
        os.remove(out_summary)

    summary = {}
    wrote_header = False
    n_nonempty_files = 0
    total_rows = 0

    print("[STEP07B] Streaming merge + summary ...")

    usecols = [
        "REGEN_ID", "CHROM", "POS", "REF", "ALT", "VAR_TYPE",
        "REGEN_GT", "SUPPORT_FRACTION", "N_READS_SUPPORT"
    ]

    for i, f in enumerate(files, start=1):
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception as e:
            print(f"[STEP07B] WARN: cannot read {f}: {e}")
            continue

        if df.empty:
            continue

        n_nonempty_files += 1
        total_rows += len(df)

        # append merged detail without keeping all files in memory
        df.to_csv(
            out_detail,
            sep="\t",
            index=False,
            mode="a",
            header=not wrote_header
        )
        wrote_header = True

        # summary update using only needed columns
        missing = [c for c in usecols if c not in df.columns]
        if missing:
            print(f"[STEP07B] WARN: {os.path.basename(f)} missing columns {missing}; filling with defaults")
            for c in missing:
                df[c] = pd.NA

        for row in df[usecols].to_dict(orient="records"):
            update_summary(summary, row)

        if i % 25 == 0 or i == len(files):
            print(
                f"[STEP07B]   Processed {i}/{len(files)} files | "
                f"nonempty={n_nonempty_files} | rows={total_rows} | unique_REGEN_ID={len(summary)}"
            )

        del df

    if not wrote_header:
        print("[STEP07B] No data in any file.")
        return

    print(f"[STEP07B] Merged detail -> {out_detail}")
    print("[STEP07B] Finalizing cohort summary ...")

    summary_df = finalize_summary(summary)
    summary_df.to_csv(out_summary, sep="\t", index=False)

    print(f"[STEP07B] Cohort summary -> {out_summary} ({len(summary_df)} rows)")

    if len(summary_df) > 0:
        n_poly = int((summary_df["N_CALLED_ALT"] > 0).sum())
        n_common = int((summary_df["N_CALLED_ALT"] >= 5).sum())
        n_rare = int(((summary_df["N_CALLED_ALT"] > 0) & (summary_df["N_CALLED_ALT"] < 5)).sum())
        n_mono = int((summary_df["N_CALLED_ALT"] == 0).sum())

        print("[STEP07B] ======= REGENOTYPE SUMMARY =======")
        print(f"  Total candidates:         {len(summary_df)}")
        print(f"  Polymorphic (ALT >= 1):   {n_poly}")
        print(f"    Common (ALT >= 5):      {n_common}")
        print(f"    Rare (1-4 ALT):         {n_rare}")
        print(f"  Monomorphic (ALT = 0):    {n_mono}")

    print(f"\n[STEP07B] Total time: {_elapsed(t0)}")


if __name__ == "__main__":
    main()
