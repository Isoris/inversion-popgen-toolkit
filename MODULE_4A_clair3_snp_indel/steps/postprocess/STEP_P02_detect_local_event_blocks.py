#!/usr/bin/env python3
"""
STEP02 – Define local event / overlap blocks.

Groups nearby or overlapping primitive VCF records into local event blocks.
A block merges records whose spans overlap OR that fall within --merge_dist bp.

Outputs:
  local_overlap_blocks.tsv  – one row per block with summary stats
  all_variants_with_blocks.tsv – input table with BLOCK_ID column added

Usage:
  python STEP02_detect_local_event_blocks.py \
      --annotated  all_variants_annotated.tsv \
      --outdir     <output_directory> \
      [--merge_dist 20]
"""

import os, argparse
import numpy as np
import pandas as pd


def build_blocks(df, merge_dist=20):
    """Greedy interval-merge per chromosome.
    Two records are in the same block if:
      - their [POS, END] intervals overlap, OR
      - the gap between them is <= merge_dist bp.
    """
    df = df.copy()
    df["BLOCK_ID"] = -1

    block_rows = []
    block_id = 0

    for chrom, grp in df.groupby("CHROM", sort=False):
        g = grp.sort_values("POS")
        idx = g.index.tolist()
        if len(idx) == 0:
            continue

        i = 0
        while i < len(idx):
            members = [idx[i]]
            blk_start = int(df.loc[idx[i], "POS"])
            blk_end = int(df.loc[idx[i], "END"])

            j = i + 1
            while j < len(idx):
                pos_j = int(df.loc[idx[j], "POS"])
                end_j = int(df.loc[idx[j], "END"])
                # merge if overlapping or within merge_dist
                if pos_j <= blk_end + merge_dist:
                    members.append(idx[j])
                    blk_end = max(blk_end, end_j)
                    j += 1
                else:
                    break

            # assign block id
            for m in members:
                df.loc[m, "BLOCK_ID"] = block_id

            # summarize block
            sub = df.loc[members]
            types_present = ",".join(sorted(sub["VAR_TYPE"].dropna().unique()))
            n_pass = int((sub["STATUS"] == "PASS").sum())
            n_filt = int((sub["STATUS"] == "FILTERED").sum())
            max_qual = sub["QUAL"].max() if not sub["QUAL"].isna().all() else np.nan
            gt_classes = ",".join(sorted(sub["GT_CLASS"].dropna().unique()))

            block_rows.append({
                "BLOCK_ID": block_id,
                "CHROM": chrom,
                "BLOCK_START": blk_start,
                "BLOCK_END": blk_end,
                "BLOCK_SPAN": blk_end - blk_start,
                "N_RECORDS": len(members),
                "TYPES_PRESENT": types_present,
                "N_PASS": n_pass,
                "N_FILTERED": n_filt,
                "MAX_QUAL": max_qual,
                "GT_CLASSES": gt_classes,
                "IS_COMPLEX": 1 if len(members) > 1 else 0,
            })

            block_id += 1
            i = j

    blocks_df = pd.DataFrame(block_rows)
    return df, blocks_df


def main():
    ap = argparse.ArgumentParser(description="STEP02: Local event / overlap blocks")
    ap.add_argument("--annotated", required=True, help="all_variants_annotated.tsv from STEP01")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--merge_dist", type=int, default=20,
                    help="Maximum gap (bp) to merge adjacent records into a block [20]")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.annotated, sep="\t")
    print(f"[STEP02] Input rows: {len(df)}")

    df, blocks = build_blocks(df, merge_dist=args.merge_dist)

    out_blocks = os.path.join(args.outdir, "local_overlap_blocks.tsv")
    blocks.to_csv(out_blocks, sep="\t", index=False)
    print(f"[STEP02] Blocks: {len(blocks)} → {out_blocks}")

    out_annot = os.path.join(args.outdir, "all_variants_with_blocks.tsv")
    df.to_csv(out_annot, sep="\t", index=False)
    print(f"[STEP02] Annotated variants with BLOCK_ID → {out_annot}")

    # quick stats
    n_complex = int(blocks["IS_COMPLEX"].sum())
    n_singleton = int((blocks["N_RECORDS"] == 1).sum())
    print(f"[STEP02] Singleton blocks: {n_singleton}, Multi-record blocks: {n_complex}")


if __name__ == "__main__":
    main()
