#!/usr/bin/env python3
"""
01_build_variant_catalog.py — Build population-level variant catalogs + GT matrices

Scans all per-sample postprocess results for a given chromosome and builds:
  - Population catalog: one row per unique variant with frequency stats
  - GT matrix: sample × variant genotype matrix
  - Binary matrix: 0/1 presence for PCA
  - Separate outputs for SNP, INDEL, and COMBINED

Input:
  Per-sample directories containing:
    all_variants_with_phase.tsv  (or all_variants_with_blocks.tsv)
    strong_filtered_rescued.tsv  (STEP03 output, if available)
    final_variant_classification.tsv (STEP08 output, if available)

Output (per variant type):
  catalogs/{TYPE}_catalog.tsv
  matrices/GT_matrix.{TYPE}.tsv
  matrices/binary_genotype_matrix.{TYPE}.tsv
  catalogs/{TYPE}_frequency_spectrum.tsv
"""

import os, sys, argparse, time, glob
from collections import defaultdict, Counter
import numpy as np

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pp_results", required=True, help="postprocess_results/<CHROM>/ directory")
    p.add_argument("--chrom", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--samples_unrelated", default="")
    p.add_argument("--min_qual", type=float, default=0, help="Min QUAL to include [0]")
    return p.parse_args()

def load_sample_variants(sample_dir, sample_id):
    """Load variants from best available per-sample file."""
    # Priority: final_variant_classification > all_variants_with_phase > all_variants_with_blocks
    for fname in ["final_variant_classification.tsv",
                   "all_variants_with_phase.tsv",
                   "all_variants_with_blocks.tsv",
                   "all_variants_annotated.tsv"]:
        path = os.path.join(sample_dir, fname)
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            rows = []
            with open(path) as f:
                header = f.readline().strip().split('\t')
                col = {h: i for i, h in enumerate(header)}
                for line in f:
                    p = line.strip().split('\t')
                    if len(p) < len(header):
                        continue
                    rows.append({h: p[i] for i, h in enumerate(header)})
            return rows, fname
    return [], None

def classify_variant_type(ref, alt, var_type_col=""):
    """Classify as SNP or INDEL."""
    if var_type_col:
        vt = var_type_col.upper()
        if vt == "SNP" or vt == "SNV":
            return "SNP"
        if vt in ("INS", "DEL", "INDEL", "MNP"):
            return "INDEL"
    # Fallback: length-based
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"

def make_variant_key(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}"

def main():
    args = parse_args()
    t0 = _now()
    os.makedirs(os.path.join(args.outdir, "catalogs"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "matrices"), exist_ok=True)

    # Find sample directories
    chrom_dir = args.pp_results
    sample_dirs = []
    for d in sorted(os.listdir(chrom_dir)):
        full = os.path.join(chrom_dir, d)
        if os.path.isdir(full) and not d.startswith("_"):
            sample_dirs.append((d, full))

    print(f"[01] Found {len(sample_dirs)} sample directories in {chrom_dir}")

    # Load unrelated samples
    unrel = set()
    if args.samples_unrelated and os.path.isfile(args.samples_unrelated):
        with open(args.samples_unrelated) as f:
            for line in f:
                s = line.strip()
                if s:
                    unrel.add(s)
        print(f"[01] Unrelated samples: {len(unrel)}")

    # ── Collect all variants across samples ──
    # variant_key → {info dict, set of carrier samples, per-sample GT}
    variant_info = {}       # key → first-seen info dict
    variant_carriers = defaultdict(set)   # key → set of sample_ids
    variant_gts = defaultdict(dict)       # key → {sample: gt_string}
    variant_phase = defaultdict(dict)     # key → {sample: {tier, block_id, phase_gt}}
    sample_list = []

    for si, (sample_id, sdir) in enumerate(sample_dirs):
        rows, source_file = load_sample_variants(sdir, sample_id)
        if not rows:
            print(f"[01]   WARN: no data for {sample_id}")
            continue
        sample_list.append(sample_id)

        for r in rows:
            chrom = r.get("CHROM", args.chrom)
            pos = r.get("POS", "0")
            ref = r.get("REF", ".")
            alt = r.get("ALT1", r.get("ALT", "."))
            gt = r.get("GT", r.get("PHASE_GT", "./."))
            qual = float(r.get("QUAL", 0)) if r.get("QUAL", ".") != "." else 0

            if args.min_qual > 0 and qual < args.min_qual:
                continue

            key = make_variant_key(chrom, pos, ref, alt)

            # Store first-seen info
            if key not in variant_info:
                vtype = classify_variant_type(ref, alt, r.get("VAR_TYPE", ""))
                variant_info[key] = {
                    "CHROM": chrom,
                    "POS": int(pos),
                    "REF": ref,
                    "ALT": alt,
                    "VAR_TYPE": vtype,
                    "INDEL_LEN": r.get("INDEL_LEN", 0),
                    "ABS_INDEL_LEN": r.get("ABS_INDEL_LEN", 0),
                }

            # Is this sample a carrier? (non-ref, non-missing)
            gt_clean = gt.replace("|", "/")
            alleles = gt_clean.split("/")
            is_carrier = any(a not in ("0", ".", "") for a in alleles)
            is_missing = all(a in (".", "") for a in alleles)

            if is_carrier:
                variant_carriers[key].add(sample_id)

            if not is_missing:
                variant_gts[key][sample_id] = gt

            # Phase info
            phase_tier = r.get("PHASE_TIER", "UNPHASED")
            phase_block = r.get("PHASE_BLOCK_ID", "-1")
            phase_gt = r.get("PHASE_GT", "")
            if phase_tier != "UNPHASED":
                variant_phase[key][sample_id] = {
                    "tier": phase_tier,
                    "block_id": phase_block,
                    "phase_gt": phase_gt,
                }

        if (si + 1) % 50 == 0:
            print(f"[01]   Loaded {si+1}/{len(sample_dirs)} samples, "
                  f"{len(variant_info)} unique variants so far …")

    print(f"[01] Total unique variants: {len(variant_info)}")
    print(f"[01] Samples with data: {len(sample_list)}")
    n_total = len(sample_list)
    n_unrel = len([s for s in sample_list if s in unrel])

    # ── Build catalogs and matrices per type ──
    for vtype in ["SNP", "INDEL", "COMBINED"]:
        print(f"\n[01] === Building {vtype} outputs ===")

        # Filter variants by type
        if vtype == "COMBINED":
            keys = sorted(variant_info.keys(),
                          key=lambda k: (variant_info[k]["CHROM"], variant_info[k]["POS"]))
        else:
            keys = sorted(
                [k for k, v in variant_info.items() if v["VAR_TYPE"] == vtype],
                key=lambda k: (variant_info[k]["CHROM"], variant_info[k]["POS"])
            )

        print(f"[01]   Variants: {len(keys)}")

        if not keys:
            print(f"[01]   No {vtype} variants, skipping")
            continue

        # ── Catalog ──
        cat_path = os.path.join(args.outdir, "catalogs", f"{vtype}_catalog.tsv")
        freq_counts = Counter()

        with open(cat_path, 'w') as o:
            o.write("\t".join([
                "VAR_KEY", "CHROM", "POS", "REF", "ALT", "VAR_TYPE",
                "INDEL_LEN", "ABS_INDEL_LEN",
                "N_CARRIERS_ALL", "N_CARRIERS_UNREL",
                "N_GENOTYPED_ALL", "MISSING_FRAC",
                "CARRIER_FREQ_ALL", "CARRIER_FREQ_UNREL",
                "FREQ_CLASS", "N_PHASED_SAMPLES",
            ]) + "\n")

            for key in keys:
                vi = variant_info[key]
                carriers = variant_carriers.get(key, set())
                n_carrier = len(carriers)
                n_carrier_unrel = len(carriers & unrel)
                n_geno = len(variant_gts.get(key, {}))
                n_missing = n_total - n_geno
                miss_frac = round(n_missing / n_total, 4) if n_total > 0 else 0
                cf_all = round(n_carrier / n_total, 4) if n_total > 0 else 0
                cf_unrel = round(n_carrier_unrel / n_unrel, 4) if n_unrel > 0 else 0

                # Frequency class
                if n_carrier == 0:
                    fc = "monomorphic"
                elif n_carrier == 1:
                    fc = "singleton"
                elif n_carrier <= 5:
                    fc = "rare"
                elif n_carrier <= 20:
                    fc = "low_freq"
                elif n_carrier <= int(n_total * 0.8):
                    fc = "common"
                else:
                    fc = "high_freq"
                freq_counts[fc] += 1

                n_phased = len(variant_phase.get(key, {}))

                o.write("\t".join(str(x) for x in [
                    key, vi["CHROM"], vi["POS"], vi["REF"], vi["ALT"], vi["VAR_TYPE"],
                    vi["INDEL_LEN"], vi["ABS_INDEL_LEN"],
                    n_carrier, n_carrier_unrel,
                    n_geno, miss_frac,
                    cf_all, cf_unrel,
                    fc, n_phased,
                ]) + "\n")

        print(f"[01]   Catalog → {cat_path}")
        for fc, cnt in sorted(freq_counts.items()):
            print(f"[01]     {fc}: {cnt}")

        # ── Frequency spectrum ──
        freq_path = os.path.join(args.outdir, "catalogs", f"{vtype}_frequency_spectrum.tsv")
        with open(freq_path, 'w') as o:
            o.write("n_carriers\tcount\n")
            carrier_hist = Counter(len(variant_carriers.get(k, set())) for k in keys)
            for nc in sorted(carrier_hist.keys()):
                o.write(f"{nc}\t{carrier_hist[nc]}\n")
        print(f"[01]   Frequency spectrum → {freq_path}")

        # ── GT matrix ──
        gt_path = os.path.join(args.outdir, "matrices", f"GT_matrix.{vtype}.tsv")
        with open(gt_path, 'w') as o:
            o.write("\t".join(["CHROM", "POS", "REF", "ALT", "VAR_KEY"] + sample_list) + "\n")
            for key in keys:
                vi = variant_info[key]
                gts_row = variant_gts.get(key, {})
                gt_cells = [gts_row.get(s, "./.") for s in sample_list]
                o.write("\t".join([
                    vi["CHROM"], str(vi["POS"]), vi["REF"], vi["ALT"], key
                ] + gt_cells) + "\n")
        print(f"[01]   GT matrix → {gt_path}")

        # ── Binary matrix ──
        bin_path = os.path.join(args.outdir, "matrices", f"binary_genotype_matrix.{vtype}.tsv")
        with open(bin_path, 'w') as o:
            o.write("\t".join(["VAR_KEY"] + sample_list) + "\n")
            for key in keys:
                carriers = variant_carriers.get(key, set())
                cells = ["1" if s in carriers else "0" for s in sample_list]
                o.write(key + "\t" + "\t".join(cells) + "\n")
        print(f"[01]   Binary matrix → {bin_path}")

    print(f"\n[01] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
