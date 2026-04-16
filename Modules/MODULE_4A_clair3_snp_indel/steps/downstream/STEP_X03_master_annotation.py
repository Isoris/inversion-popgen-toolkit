#!/usr/bin/env python3
"""
03_master_annotation.py — Build master annotation tables (DELLY-equivalent)

Merges catalog data with functional annotation, repeat overlap, phase stats.
Produces per-type master tables with frequency class, size class (INDEL),
functional class, and all summary columns needed for downstream plotting.

Input:
  - catalogs/{TYPE}_catalog.tsv (from 01)
  - Annotation BEDs (gene/exon/CDS, repeats) if available
  - Phase prep data (from 02)

Output:
  annotation/master_{TYPE}_annotation.tsv      (SNP, INDEL, COMBINED)
  annotation/frequency_class_summary.{TYPE}.tsv
  annotation/size_class_summary.INDEL.tsv
"""

import os, sys, argparse, time
from collections import defaultdict, Counter

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--catalogs_dir", required=True)
    p.add_argument("--phase_prep_dir", default="")
    p.add_argument("--annot_beds", default="", help="Directory with GENE/EXON/CDS/REPEATS sorted BED")
    p.add_argument("--chrom", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--ref_fai", default="")
    return p.parse_args()

def load_catalog(path):
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= len(header):
                rows.append(dict(zip(header, p)))
    return rows

def size_class_indel(abs_len):
    try:
        l = int(abs_len)
    except (ValueError, TypeError):
        return "unknown"
    if l == 0:
        return "unknown"
    if l <= 5:
        return "micro"
    if l <= 50:
        return "small"
    if l <= 500:
        return "medium"
    return "large"

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "annotation")
    os.makedirs(outdir, exist_ok=True)

    # Load phase summary per variant if available
    phase_per_var = {}
    pmap_path = os.path.join(args.phase_prep_dir, "phase_block_variant_map.tsv") if args.phase_prep_dir else ""
    if pmap_path and os.path.isfile(pmap_path):
        print("[03] Loading phase variant map …")
        with open(pmap_path) as f:
            header = f.readline().strip().split('\t')
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                p = line.strip().split('\t')
                vk = p[col["VAR_KEY"]]
                tier = p[col["PHASE_TIER"]]
                if vk not in phase_per_var:
                    phase_per_var[vk] = {"n_t1": 0, "n_t2": 0, "n_un": 0}
                if tier == "TIER_1_WHATSHAP":
                    phase_per_var[vk]["n_t1"] += 1
                elif tier == "TIER_2_READPAIR":
                    phase_per_var[vk]["n_t2"] += 1
                else:
                    phase_per_var[vk]["n_un"] += 1
        print(f"[03]   Phase data for {len(phase_per_var)} variants")

    # Load chromosome length
    chr_len = 0
    if args.ref_fai and os.path.isfile(args.ref_fai):
        with open(args.ref_fai) as f:
            for line in f:
                p = line.strip().split('\t')
                if p[0] == args.chrom:
                    chr_len = int(p[1])
                    break

    for vtype in ["SNP", "INDEL", "COMBINED"]:
        cat_path = os.path.join(args.catalogs_dir, f"{vtype}_catalog.tsv")
        rows = load_catalog(cat_path)
        if not rows:
            print(f"[03] No {vtype} catalog found, skipping")
            continue

        print(f"[03] === {vtype}: {len(rows)} variants ===")

        # Annotate
        freq_counts = Counter()
        size_counts = Counter()

        out_path = os.path.join(outdir, f"master_{vtype}_annotation.tsv")

        master_cols = [
            "VAR_KEY", "CHROM", "POS", "REF", "ALT", "VAR_TYPE",
            "INDEL_LEN", "ABS_INDEL_LEN", "SIZE_CLASS",
            "N_CARRIERS_ALL", "N_CARRIERS_UNREL",
            "CARRIER_FREQ_ALL", "CARRIER_FREQ_UNREL",
            "FREQ_CLASS", "MISSING_FRAC",
            "N_PHASED_SAMPLES", "N_WHATSHAP_CARRIERS", "N_READPAIR_CARRIERS",
            "WINDOW_1MB",
        ]

        with open(out_path, 'w') as o:
            o.write("\t".join(master_cols) + "\n")

            for r in rows:
                pos = int(r["POS"])
                abs_ilen = r.get("ABS_INDEL_LEN", "0")
                sc = size_class_indel(abs_ilen) if r.get("VAR_TYPE", "") != "SNP" else "SNP"
                fc = r.get("FREQ_CLASS", "unknown")
                freq_counts[fc] += 1
                if sc != "SNP":
                    size_counts[sc] += 1

                vk = r.get("VAR_KEY", "")
                pinfo = phase_per_var.get(vk, {"n_t1": 0, "n_t2": 0, "n_un": 0})

                window = (pos // 1000000) * 1000000

                o.write("\t".join(str(x) for x in [
                    vk, r["CHROM"], pos, r["REF"], r["ALT"], r.get("VAR_TYPE", "."),
                    r.get("INDEL_LEN", 0), abs_ilen, sc,
                    r.get("N_CARRIERS_ALL", 0), r.get("N_CARRIERS_UNREL", 0),
                    r.get("CARRIER_FREQ_ALL", 0), r.get("CARRIER_FREQ_UNREL", 0),
                    fc, r.get("MISSING_FRAC", 0),
                    r.get("N_PHASED_SAMPLES", 0), pinfo["n_t1"], pinfo["n_t2"],
                    window,
                ]) + "\n")

        print(f"[03]   Master annotation → {out_path}")

        # Frequency summary
        fq_path = os.path.join(outdir, f"frequency_class_summary.{vtype}.tsv")
        with open(fq_path, 'w') as o:
            o.write("frequency_class\tcount\n")
            for fc in ["monomorphic", "singleton", "rare", "low_freq", "common", "high_freq"]:
                o.write(f"{fc}\t{freq_counts.get(fc, 0)}\n")
        print(f"[03]   Frequency summary → {fq_path}")

        # Size summary (INDEL only)
        if vtype in ("INDEL", "COMBINED") and size_counts:
            sz_path = os.path.join(outdir, f"size_class_summary.{vtype}.tsv")
            with open(sz_path, 'w') as o:
                o.write("size_class\tcount\n")
                for sc in ["micro", "small", "medium", "large", "unknown"]:
                    o.write(f"{sc}\t{size_counts.get(sc, 0)}\n")
            print(f"[03]   Size class summary → {sz_path}")

    print(f"\n[03] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
