#!/usr/bin/env python3
"""
04_per_sample_summary.py — Per-sample variant burden tables

For each sample, counts: total variants, by type, by frequency class,
missingness, phase stats, ancestry. Runs PCA from binary genotype matrix.

Output:
  per_sample/per_sample_{TYPE}_summary.tsv   (SNP, INDEL, COMBINED)
"""

import os, sys, argparse, time
import numpy as np
from collections import defaultdict

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--matrices_dir", required=True)
    p.add_argument("--annotation_dir", required=True)
    p.add_argument("--phase_prep_dir", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--samples_unrelated", default="")
    p.add_argument("--qopt", default="")
    p.add_argument("--qopt_samples", default="")
    return p.parse_args()

def load_ancestry(qopt_path, qopt_samples_path):
    ancestry = {}
    if not (qopt_path and os.path.isfile(qopt_path) and
            qopt_samples_path and os.path.isfile(qopt_samples_path)):
        return ancestry
    with open(qopt_samples_path) as f:
        qsamp = [l.strip() for l in f if l.strip()]
    with open(qopt_path) as f:
        for i, line in enumerate(f):
            if i >= len(qsamp):
                break
            vals = [float(x) for x in line.strip().split()]
            maxq = max(vals)
            cluster = vals.index(maxq) + 1
            ancestry[qsamp[i]] = (f"Q{cluster}", round(maxq, 4))
    return ancestry

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "per_sample")
    os.makedirs(outdir, exist_ok=True)

    unrel = set()
    if args.samples_unrelated and os.path.isfile(args.samples_unrelated):
        with open(args.samples_unrelated) as f:
            for line in f:
                s = line.strip()
                if s:
                    unrel.add(s)

    ancestry = load_ancestry(args.qopt, args.qopt_samples)

    # Load phase summary per sample
    phase_stats = {}
    ps_path = os.path.join(args.phase_prep_dir, "phase_block_summary_per_sample.tsv") if args.phase_prep_dir else ""
    if ps_path and os.path.isfile(ps_path):
        with open(ps_path) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                p = line.strip().split('\t')
                d = dict(zip(header, p))
                phase_stats[d["SAMPLE"]] = d

    for vtype in ["SNP", "INDEL", "COMBINED"]:
        gt_path = os.path.join(args.matrices_dir, f"GT_matrix.{vtype}.tsv")
        bin_path = os.path.join(args.matrices_dir, f"binary_genotype_matrix.{vtype}.tsv")
        annot_path = os.path.join(args.annotation_dir, f"master_{vtype}_annotation.tsv")

        if not os.path.isfile(gt_path):
            print(f"[04] No GT matrix for {vtype}, skipping")
            continue

        print(f"[04] === {vtype} ===")

        # Load annotation for freq class lookup
        var_freqclass = {}
        if os.path.isfile(annot_path):
            with open(annot_path) as f:
                header = f.readline().strip().split('\t')
                col = {h: i for i, h in enumerate(header)}
                for line in f:
                    p = line.strip().split('\t')
                    vk = p[col["VAR_KEY"]]
                    fc = p[col.get("FREQ_CLASS", -1)] if "FREQ_CLASS" in col else "unknown"
                    var_freqclass[vk] = fc

        # Parse GT matrix
        with open(gt_path) as f:
            gtheader = f.readline().strip().split('\t')
            samples = gtheader[5:]

            counters = {s: defaultdict(int) for s in samples}
            n_missing = {s: 0 for s in samples}
            n_total_sites = 0

            for line in f:
                p = line.strip().split('\t')
                var_key = p[4]
                gts = p[5:]
                n_total_sites += 1

                fc = var_freqclass.get(var_key, "unknown")

                for i, gt in enumerate(gts):
                    s = samples[i]
                    gt_clean = gt.replace("|", "/")
                    alleles = gt_clean.split("/")
                    is_missing = all(a in (".", "") for a in alleles)
                    is_carrier = any(a not in ("0", ".", "") for a in alleles)

                    if is_missing:
                        n_missing[s] += 1
                    elif is_carrier:
                        counters[s]["total"] += 1
                        counters[s][f"freq_{fc}"] += 1

        # PCA
        pc_dict = {}
        var_pct = [0, 0]
        if os.path.isfile(bin_path):
            print(f"[04]   Running PCA …")
            with open(bin_path) as f:
                bheader = f.readline().strip().split('\t')
                bsamples = bheader[1:]
                rows = []
                for line in f:
                    p = line.strip().split('\t')
                    rows.append([int(x) for x in p[1:]])

            if rows:
                mat = np.array(rows, dtype=np.float64).T
                col_var = mat.var(axis=0)
                mat = mat[:, col_var > 0]
                if mat.shape[1] > 2:
                    mat_c = mat - mat.mean(axis=0)
                    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
                    pc1 = U[:, 0] * S[0]
                    pc2 = U[:, 1] * S[1]
                    total_var = np.sum(S**2)
                    var_pct = [S[0]**2 / total_var * 100, S[1]**2 / total_var * 100]
                    for i, s in enumerate(bsamples):
                        pc_dict[s] = (round(pc1[i], 4), round(pc2[i], 4))
                    print(f"[04]   PCA: PC1={var_pct[0]:.1f}%, PC2={var_pct[1]:.1f}%")

        # Write output
        out_path = os.path.join(outdir, f"per_sample_{vtype}_summary.tsv")
        with open(out_path, 'w') as o:
            cols = [
                "SAMPLE", "SUBSET", "TOTAL_VARIANTS",
                "N_SINGLETON", "N_RARE", "N_LOW_FREQ", "N_COMMON", "N_HIGH_FREQ",
                "MISSINGNESS",
                "N_TIER1_VARIANTS", "N_TIER2_VARIANTS", "N_UNPHASED",
                "ANCESTRY_CLUSTER", "MAX_Q",
                "PC1", "PC2",
            ]
            o.write("\t".join(cols) + "\n")

            for s in samples:
                c = counters[s]
                miss = round(n_missing[s] / n_total_sites, 4) if n_total_sites > 0 else 0
                subset = "unrelated" if s in unrel else "related"
                anc = ancestry.get(s, (".", "."))
                pc = pc_dict.get(s, (".", "."))
                ps = phase_stats.get(s, {})

                o.write("\t".join(str(x) for x in [
                    s, subset, c["total"],
                    c.get("freq_singleton", 0), c.get("freq_rare", 0),
                    c.get("freq_low_freq", 0), c.get("freq_common", 0),
                    c.get("freq_high_freq", 0),
                    miss,
                    ps.get("N_TIER1_VARIANTS", 0), ps.get("N_TIER2_VARIANTS", 0),
                    ps.get("N_UNPHASED", 0),
                    anc[0], anc[1],
                    pc[0], pc[1],
                ]) + "\n")

        print(f"[04]   {len(samples)} samples → {out_path}")

    print(f"\n[04] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
