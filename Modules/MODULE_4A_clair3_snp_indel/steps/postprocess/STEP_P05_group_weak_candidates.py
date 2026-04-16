#!/usr/bin/env python3
"""
STEP05 – Group weak indel candidates across samples (population rescue).

Takes many weak_indel_candidates.tsv files (one per sample) and clusters
them by local event similarity.  An event cluster is a group of candidates
from different samples that likely represent the same biological indel.

Grouping strategy (hierarchical):
  Level 1 – Exact: same chrom:pos:ref:alt (left-normalized)
  Level 2 – Local: same chrom, overlapping span, same var_type, same |indel_len|
  Level 3 – Motif: Level 2 + matching inserted/deleted motif

Each cluster gets a POPULATION_CLUSTER_ID and sample-count stats.

Outputs:
  weak_indel_population_clusters.tsv – one row per sample×candidate with cluster IDs
  weak_indel_population_catalog.tsv  – one row per cluster with summary stats

Usage:
  python STEP05_group_weak_candidates_across_samples.py \
      --input_dir   <directory_with_per_sample_weak_indel_candidates> \
      --outdir      <output_directory> \
      [--span_slop  10] \
      [--min_samples 2]

  OR provide a manifest:
  python STEP05_group_weak_candidates_across_samples.py \
      --manifest    <manifest.tsv>   (one weak_indel_candidates.tsv path per line) \
      --outdir      <output_directory>
"""

import os, sys, glob, argparse
from collections import defaultdict

import numpy as np
import pandas as pd


def load_all_candidates(input_dir=None, manifest=None):
    """Load and concatenate weak_indel_candidates.tsv from all samples."""
    files = []
    if manifest:
        with open(manifest) as fh:
            for line in fh:
                p = line.strip()
                if p and os.path.isfile(p):
                    files.append(p)
    elif input_dir:
        # Recursively find all weak_indel_candidates.tsv
        pattern = os.path.join(input_dir, "**", "weak_indel_candidates.tsv")
        files = sorted(glob.glob(pattern, recursive=True))
    else:
        raise ValueError("Need --input_dir or --manifest")

    print(f"[STEP05] Found {len(files)} candidate files")

    dfs = []
    for f in files:
        try:
            d = pd.read_csv(f, sep="\t")
            if len(d) > 0:
                dfs.append(d)
                print(f"  {f}: {len(d)} candidates")
        except Exception as e:
            print(f"[WARN] Cannot read {f}: {e}", file=sys.stderr)

    if not dfs:
        print("[STEP05] No candidates loaded!")
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    print(f"[STEP05] Total candidates loaded: {len(combined)} from {len(dfs)} samples")
    return combined


def cluster_exact(df):
    """Level 1: Exact chrom:pos:ref:alt grouping."""
    df["EXACT_KEY"] = (
        df["CHROM"].astype(str) + ":" +
        df["POS"].astype(str) + ":" +
        df["REF"].astype(str) + ":" +
        df["ALT"].astype(str)
    )
    return df


def cluster_local(df, span_slop=10):
    """Level 2: Group by chrom + overlapping span + same type + same |indel_len|.

    Uses a greedy merge within each (chrom, var_type, abs_indel_len) group.
    """
    df["LOCAL_CLUSTER_ID"] = -1
    cluster_id = 0

    for (chrom, vtype, alen), grp in df.groupby(
        ["CHROM", "VAR_TYPE", "ABS_INDEL_LEN"], sort=False
    ):
        if pd.isna(alen):
            continue
        g = grp.sort_values("POS")
        idx_list = g.index.tolist()
        if not idx_list:
            continue

        i = 0
        while i < len(idx_list):
            members = [idx_list[i]]
            merge_end = int(df.loc[idx_list[i], "END"]) + span_slop

            j = i + 1
            while j < len(idx_list):
                pos_j = int(df.loc[idx_list[j], "POS"])
                if pos_j <= merge_end:
                    members.append(idx_list[j])
                    merge_end = max(merge_end,
                                    int(df.loc[idx_list[j], "END"]) + span_slop)
                    j += 1
                else:
                    break

            for m in members:
                df.loc[m, "LOCAL_CLUSTER_ID"] = cluster_id
            cluster_id += 1
            i = j

    return df


def cluster_motif(df):
    """Level 3: Within each LOCAL_CLUSTER_ID, sub-group by INDEL_MOTIF identity."""
    df["MOTIF_CLUSTER_ID"] = -1
    cluster_id = 0

    for lcid, grp in df.groupby("LOCAL_CLUSTER_ID"):
        if lcid == -1:
            continue
        motif_groups = grp.groupby("INDEL_MOTIF", sort=False)
        for motif, mgrp in motif_groups:
            for idx in mgrp.index:
                df.loc[idx, "MOTIF_CLUSTER_ID"] = cluster_id
            cluster_id += 1

    return df


def build_catalog(df, min_samples=2):
    """Build population catalog from clustered candidates.

    Uses MOTIF_CLUSTER_ID as the finest cluster level.
    Falls back to LOCAL_CLUSTER_ID if motif is not available.
    """
    # Use MOTIF_CLUSTER_ID where available, else LOCAL_CLUSTER_ID
    df["POP_CLUSTER_ID"] = df["MOTIF_CLUSTER_ID"]
    mask_no_motif = df["MOTIF_CLUSTER_ID"] == -1
    if mask_no_motif.any():
        # offset local cluster IDs to avoid collision
        offset = df["MOTIF_CLUSTER_ID"].max() + 1 if df["MOTIF_CLUSTER_ID"].max() >= 0 else 0
        df.loc[mask_no_motif, "POP_CLUSTER_ID"] = (
            df.loc[mask_no_motif, "LOCAL_CLUSTER_ID"] + offset
        )

    catalog_rows = []
    for cid, grp in df.groupby("POP_CLUSTER_ID"):
        if cid == -1:
            continue

        samples = sorted(grp["SAMPLE_ID"].dropna().unique())
        n_samples = len(samples)

        # representative event = row with highest QUAL
        rep = grp.loc[grp["QUAL"].idxmax()] if not grp["QUAL"].isna().all() else grp.iloc[0]

        catalog_rows.append({
            "POP_CLUSTER_ID": cid,
            "CHROM": rep["CHROM"],
            "POS_REPRESENTATIVE": int(rep["POS"]),
            "END_REPRESENTATIVE": int(rep["END"]),
            "REF": rep["REF"],
            "ALT": rep["ALT"],
            "VAR_TYPE": rep["VAR_TYPE"],
            "INDEL_LEN": rep["INDEL_LEN"],
            "ABS_INDEL_LEN": rep.get("ABS_INDEL_LEN", np.nan),
            "INDEL_MOTIF": rep.get("INDEL_MOTIF", ""),
            "EVENT_SIGNATURE": rep.get("EVENT_SIGNATURE", ""),
            "N_SAMPLES": n_samples,
            "SAMPLE_IDS": ",".join(samples),
            "N_TOTAL_OBS": len(grp),
            "MEAN_QUAL": grp["QUAL"].mean(),
            "MAX_QUAL": grp["QUAL"].max(),
            "MEAN_AF1": grp["AF1"].mean(),
            "MEAN_DP": grp["DP"].mean(),
            "MEAN_AD_ALT": grp["AD_ALT"].mean(),
            "TOTAL_READS_SUPPORT": grp["N_READS_SUPPORT_INDEL"].sum(),
            "LEFT_FLANK": rep.get("LEFT_FLANK_20BP", ""),
            "RIGHT_FLANK": rep.get("RIGHT_FLANK_20BP", ""),
            "PASSES_POP_RESCUE": 1 if n_samples >= min_samples else 0,
        })

    catalog = pd.DataFrame(catalog_rows)
    return catalog


def main():
    ap = argparse.ArgumentParser(description="STEP05: Group weak candidates across samples")
    ap.add_argument("--input_dir", default=None,
                    help="Directory tree containing per-sample weak_indel_candidates.tsv files")
    ap.add_argument("--manifest", default=None,
                    help="File listing paths to weak_indel_candidates.tsv (one per line)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--span_slop", type=int, default=10)
    ap.add_argument("--min_samples", type=int, default=2,
                    help="Minimum samples for population rescue [2]")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = load_all_candidates(input_dir=args.input_dir, manifest=args.manifest)
    if df.empty:
        print("[STEP05] Nothing to cluster.")
        return

    print("[STEP05] Clustering Level 1 (exact) …")
    df = cluster_exact(df)

    print("[STEP05] Clustering Level 2 (local span) …")
    df = cluster_local(df, span_slop=args.span_slop)

    print("[STEP05] Clustering Level 3 (motif) …")
    df = cluster_motif(df)

    print("[STEP05] Building population catalog …")
    catalog = build_catalog(df, min_samples=args.min_samples)

    # outputs
    out_clusters = os.path.join(args.outdir, "weak_indel_population_clusters.tsv")
    df.to_csv(out_clusters, sep="\t", index=False)
    print(f"[STEP05] Clustered candidates → {out_clusters}")

    out_catalog = os.path.join(args.outdir, "weak_indel_population_catalog.tsv")
    catalog.to_csv(out_catalog, sep="\t", index=False)
    print(f"[STEP05] Population catalog → {out_catalog}")

    # summary
    n_pass = int(catalog["PASSES_POP_RESCUE"].sum()) if len(catalog) > 0 else 0
    print(f"[STEP05] Total clusters: {len(catalog)}")
    print(f"[STEP05] Clusters with >= {args.min_samples} samples (rescue-eligible): {n_pass}")


if __name__ == "__main__":
    main()
