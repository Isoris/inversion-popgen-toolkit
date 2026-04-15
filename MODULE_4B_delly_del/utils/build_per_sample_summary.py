#!/usr/bin/env python3
"""
08_build_per_sample_summary.py — Build per-sample DEL summary table
Columns: sample, total_DEL, burden_small, burden_medium, burden_large,
         burden_repeat, burden_nonrepeat, burden_CDS, burden_exon, burden_intronic,
         burden_intergenic, missingness, ancestry_cluster, max_Q, PC1, PC2

Inputs:
  - master_DEL_annotation_226.tsv (from 07)
  - catalog_226.DEL.GT_matrix.tsv
  - samples_unrelated_81.txt
  - Optional: qopt + sample list for ancestry
  - DEL_binary_genotype_matrix.tsv (for PCA)

Output:
  - 11_summary/per_sample_DEL_summary.tsv
"""
import os, sys, argparse
import numpy as np
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--master_annot", required=True)
    p.add_argument("--gt_matrix", required=True)
    p.add_argument("--binary_gt", required=True)
    p.add_argument("--samples_unrelated", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--qopt", default="")
    p.add_argument("--qopt_samples", default="")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load unrelated
    unrel = set()
    with open(args.samples_unrelated) as f:
        for line in f:
            s = line.strip()
            if s: unrel.add(s)

    # Load ancestry
    ancestry = {}
    if args.qopt and os.path.exists(args.qopt) and args.qopt_samples and os.path.exists(args.qopt_samples):
        with open(args.qopt_samples) as f:
            qsamp = [l.strip() for l in f if l.strip()]
        with open(args.qopt) as f:
            for i, line in enumerate(f):
                if i >= len(qsamp): break
                vals = [float(x) for x in line.strip().split()]
                maxq = max(vals)
                cluster = vals.index(maxq) + 1
                ancestry[qsamp[i]] = (f"Q{cluster}", round(maxq, 4))

    # Load master annotation into dict by del_id
    print("Loading master annotation...")
    annot = {}
    with open(args.master_annot) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            annot[d['del_id']] = d

    # Parse GT matrix
    print("Parsing GT matrix...")
    with open(args.gt_matrix) as f:
        gtheader = f.readline().strip().split('\t')
        samples = gtheader[5:]
        n_samples = len(samples)

        # Initialize per-sample counters
        counters = {s: defaultdict(int) for s in samples}
        n_missing = {s: 0 for s in samples}
        n_total_sites = 0

        for line in f:
            p = line.strip().split('\t')
            del_id = p[3]
            gts = p[5:]
            n_total_sites += 1

            a = annot.get(del_id, {})
            sc = a.get('size_class', 'unknown')
            fc = a.get('func_class', 'intergenic')
            rep = int(a.get('repeat_50pct', 0))

            for i, gt in enumerate(gts):
                s = samples[i]
                if gt in ('./.', '.', './0', '0/.'):
                    n_missing[s] += 1
                    continue
                if gt in ('0/0', '0|0'):
                    continue
                # Has alt
                counters[s]['total'] += 1
                counters[s][f'size_{sc}'] += 1
                counters[s][f'func_{fc}'] += 1
                if rep:
                    counters[s]['repeat'] += 1
                else:
                    counters[s]['nonrepeat'] += 1

    # PCA from binary genotype matrix
    print("Running PCA...")
    with open(args.binary_gt) as f:
        bheader = f.readline().strip().split('\t')
        bsamples = bheader[1:]
        rows = []
        for line in f:
            p = line.strip().split('\t')
            rows.append([int(x) for x in p[1:]])

    mat = np.array(rows, dtype=np.float64).T  # samples x DELs
    # Remove zero-variance columns
    col_var = mat.var(axis=0)
    mat = mat[:, col_var > 0]
    # Center
    mat_c = mat - mat.mean(axis=0)
    # SVD for PCA
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]
    total_var = np.sum(S**2)
    var_pct = S**2 / total_var * 100

    pc_dict = {}
    for i, s in enumerate(bsamples):
        pc_dict[s] = (round(pc1[i], 4), round(pc2[i], 4))

    # Write output
    out_file = os.path.join(args.outdir, "per_sample_DEL_summary.tsv")
    with open(out_file, 'w') as o:
        cols = [
            "sample", "subset", "total_DEL",
            "burden_small", "burden_medium", "burden_large",
            "burden_repeat", "burden_nonrepeat",
            "burden_CDS", "burden_exon", "burden_intronic", "burden_intergenic",
            "missingness",
            "ancestry_cluster", "max_Q",
            "PC1", "PC2",
        ]
        o.write('\t'.join(cols) + '\n')
        for s in samples:
            c = counters[s]
            miss = round(n_missing[s] / n_total_sites, 4) if n_total_sites > 0 else 0
            anc = ancestry.get(s, ('.', '.'))
            pc = pc_dict.get(s, ('.', '.'))
            subset = "unrelated_81" if s in unrel else "related"
            row = [
                s, subset, c['total'],
                c.get('size_small', 0), c.get('size_medium', 0), c.get('size_large', 0),
                c.get('repeat', 0), c.get('nonrepeat', 0),
                c.get('func_CDS_overlap', 0), c.get('func_exon_overlap', 0),
                c.get('func_intronic', 0), c.get('func_intergenic', 0),
                miss,
                anc[0], anc[1],
                pc[0], pc[1],
            ]
            o.write('\t'.join(str(x) for x in row) + '\n')

    print(f"Per-sample summary: {len(samples)} samples → {out_file}")
    print(f"  PCA: PC1={var_pct[0]:.1f}%, PC2={var_pct[1]:.1f}%")

if __name__ == "__main__":
    main()
