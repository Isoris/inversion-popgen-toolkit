#!/usr/bin/env python3
"""
STEP06 – Prepare regenotype catalog.

Takes the population catalog from STEP05, selects clusters that pass the
population rescue threshold (>= min_samples), and produces:
  - A TSV catalog of candidate indels to regenotype
  - A BED file of target regions (with padding for pysam fetch)
  - A simple normalized VCF-like representation

Outputs:
  rescued_indel_regenotype_catalog.tsv  – candidate indels for regenotyping
  regenotype_targets.bed                – BED regions for BAM extraction
  regenotype_candidates.vcf             – pseudo-VCF for reference

Usage:
  python STEP06_prepare_regenotype_catalog.py \
      --catalog   weak_indel_population_catalog.tsv \
      --outdir    <output_directory> \
      [--min_samples  2] \
      [--pad_bp       50]
"""

import os, argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="STEP06: Prepare regenotype catalog")
    ap.add_argument("--catalog", required=True,
                    help="weak_indel_population_catalog.tsv from STEP05")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min_samples", type=int, default=2)
    ap.add_argument("--pad_bp", type=int, default=50,
                    help="Padding around target positions for BAM extraction [50]")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    cat = pd.read_csv(args.catalog, sep="\t")
    print(f"[STEP06] Input catalog: {len(cat)} clusters")

    # Filter to rescue-eligible
    rescued = cat[cat["PASSES_POP_RESCUE"] == 1].copy()
    print(f"[STEP06] Clusters with >= {args.min_samples} samples: {len(rescued)}")

    if len(rescued) == 0:
        print("[STEP06] No clusters pass population rescue threshold.")
        # write empty outputs
        pd.DataFrame().to_csv(
            os.path.join(args.outdir, "rescued_indel_regenotype_catalog.tsv"),
            sep="\t", index=False)
        return

    # assign regenotype IDs
    rescued = rescued.reset_index(drop=True)
    rescued["REGEN_ID"] = ["REGEN_" + str(i).zfill(6) for i in range(len(rescued))]

    # TSV catalog
    out_tsv = os.path.join(args.outdir, "rescued_indel_regenotype_catalog.tsv")
    rescued.to_csv(out_tsv, sep="\t", index=False)
    print(f"[STEP06] Catalog → {out_tsv}")

    # BED file
    bed_rows = []
    for _, row in rescued.iterrows():
        chrom = row["CHROM"]
        pos = int(row["POS_REPRESENTATIVE"])
        end = int(row["END_REPRESENTATIVE"])
        start_bed = max(0, pos - 1 - args.pad_bp)
        end_bed = end + args.pad_bp
        name = row["REGEN_ID"]
        bed_rows.append(f"{chrom}\t{start_bed}\t{end_bed}\t{name}")

    out_bed = os.path.join(args.outdir, "regenotype_targets.bed")
    with open(out_bed, "w") as fh:
        fh.write("\n".join(bed_rows) + "\n")
    print(f"[STEP06] BED → {out_bed}")

    # pseudo-VCF
    out_vcf = os.path.join(args.outdir, "regenotype_candidates.vcf")
    with open(out_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=REGEN_ID,Number=1,Type=String,Description="Regenotype catalog ID">\n')
        fh.write('##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description="Number of samples in discovery">\n')
        fh.write('##INFO=<ID=MOTIF,Number=1,Type=String,Description="Indel motif">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in rescued.iterrows():
            info = (f"REGEN_ID={row['REGEN_ID']};"
                    f"N_SAMPLES={row['N_SAMPLES']};"
                    f"MOTIF={row.get('INDEL_MOTIF', '.')}")
            fh.write(f"{row['CHROM']}\t{int(row['POS_REPRESENTATIVE'])}\t"
                      f"{row['REGEN_ID']}\t{row['REF']}\t{row['ALT']}\t"
                      f".\t.\t{info}\n")
    print(f"[STEP06] VCF → {out_vcf}")


if __name__ == "__main__":
    main()
