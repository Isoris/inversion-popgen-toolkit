#!/usr/bin/env python3
"""
07_gene_chr_tables.py — Per-chromosome and gene-overlap summary tables

Outputs:
  gene_tables/per_chromosome_summary.{TYPE}.tsv
  gene_tables/variant_density_windows.{TYPE}.tsv  (for genomic landscape plots)
"""

import os, argparse, time
from collections import defaultdict, Counter

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--annotation_dir", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--ref_fai", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--window_size", type=int, default=100000, help="Window size for density [100kb]")
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "gene_tables")
    os.makedirs(outdir, exist_ok=True)

    # Load chromosome lengths
    chr_lengths = {}
    chr_order = []
    with open(args.ref_fai) as f:
        for line in f:
            p = line.strip().split('\t')
            chr_lengths[p[0]] = int(p[1])
            chr_order.append(p[0])

    for vtype in ["SNP", "INDEL", "COMBINED"]:
        annot_path = os.path.join(args.annotation_dir, f"master_{vtype}_annotation.tsv")
        if not os.path.isfile(annot_path):
            continue

        print(f"[07] === {vtype} ===")

        # Load variants
        variants = []
        with open(annot_path) as f:
            header = f.readline().strip().split('\t')
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                p = line.strip().split('\t')
                variants.append(dict(zip(header, p)))

        print(f"[07]   Loaded {len(variants)} variants")

        # ── Per-chromosome summary ──
        chr_data = defaultdict(lambda: {
            "n_var": 0, "n_singleton": 0, "n_rare": 0,
            "n_common": 0, "n_high_freq": 0,
        })
        for v in variants:
            c = v["CHROM"]
            chr_data[c]["n_var"] += 1
            fc = v.get("FREQ_CLASS", "unknown")
            if fc == "singleton":
                chr_data[c]["n_singleton"] += 1
            elif fc == "rare":
                chr_data[c]["n_rare"] += 1
            elif fc == "common":
                chr_data[c]["n_common"] += 1
            elif fc == "high_freq":
                chr_data[c]["n_high_freq"] += 1

        chr_path = os.path.join(outdir, f"per_chromosome_summary.{vtype}.tsv")
        with open(chr_path, 'w') as o:
            o.write("\t".join([
                "CHROM", "LENGTH_MB", "N_VARIANTS", "VAR_PER_MB",
                "N_SINGLETON", "N_RARE", "N_COMMON", "N_HIGH_FREQ",
            ]) + "\n")
            for c in chr_order:
                if c not in chr_data and c != args.chrom:
                    continue
                cd = chr_data.get(c, {"n_var": 0, "n_singleton": 0, "n_rare": 0, "n_common": 0, "n_high_freq": 0})
                length_mb = chr_lengths.get(c, 0) / 1e6
                vpm = round(cd["n_var"] / length_mb, 2) if length_mb > 0 else 0
                o.write("\t".join(str(x) for x in [
                    c, round(length_mb, 2), cd["n_var"], vpm,
                    cd["n_singleton"], cd["n_rare"], cd["n_common"], cd["n_high_freq"],
                ]) + "\n")
        print(f"[07]   Per-chromosome → {chr_path}")

        # ── Variant density in windows ──
        chrom_len = chr_lengths.get(args.chrom, 0)
        ws = args.window_size
        n_windows = (chrom_len // ws) + 1 if chrom_len > 0 else 0

        window_total = [0] * n_windows
        window_snp = [0] * n_windows
        window_indel = [0] * n_windows
        window_singleton = [0] * n_windows
        window_common = [0] * n_windows

        for v in variants:
            if v["CHROM"] != args.chrom:
                continue
            try:
                pos = int(v["POS"])
            except:
                continue
            wi = pos // ws
            if wi >= n_windows:
                continue
            window_total[wi] += 1
            vt = v.get("VAR_TYPE", "")
            if vt == "SNP":
                window_snp[wi] += 1
            else:
                window_indel[wi] += 1
            fc = v.get("FREQ_CLASS", "")
            if fc == "singleton":
                window_singleton[wi] += 1
            elif fc in ("common", "high_freq"):
                window_common[wi] += 1

        wd_path = os.path.join(outdir, f"variant_density_windows.{vtype}.tsv")
        with open(wd_path, 'w') as o:
            o.write("\t".join([
                "CHROM", "WINDOW_START", "WINDOW_END",
                "N_TOTAL", "N_SNP", "N_INDEL",
                "N_SINGLETON", "N_COMMON",
            ]) + "\n")
            for wi in range(n_windows):
                o.write("\t".join(str(x) for x in [
                    args.chrom, wi * ws, (wi + 1) * ws,
                    window_total[wi], window_snp[wi], window_indel[wi],
                    window_singleton[wi], window_common[wi],
                ]) + "\n")
        print(f"[07]   Density windows → {wd_path}")

    print(f"\n[07] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
