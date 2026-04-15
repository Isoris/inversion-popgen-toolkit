#!/usr/bin/env python3
import os
import argparse
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

    unrel = set()
    with open(args.samples_unrelated) as f:
        for line in f:
            s = line.strip()
            if s:
                unrel.add(s)

    ancestry = {}
    if args.qopt and os.path.exists(args.qopt) and args.qopt_samples and os.path.exists(args.qopt_samples):
        with open(args.qopt_samples) as f:
            qsamp = [x.strip() for x in f if x.strip()]
        with open(args.qopt) as f:
            for i, line in enumerate(f):
                if i >= len(qsamp):
                    break
                vals = [float(x) for x in line.strip().split()]
                maxq = max(vals)
                cluster = vals.index(maxq) + 1
                ancestry[qsamp[i]] = (f"Q{cluster}", round(maxq, 4))

    annot = {}
    with open(args.master_annot) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d = dict(zip(header, p))
            annot[d["dup_id"]] = d

    with open(args.gt_matrix) as f:
        gthead = f.readline().rstrip("\n").split("\t")
        samples = gthead[5:]

        counters = {s: defaultdict(int) for s in samples}
        missing = {s: 0 for s in samples}
        n_total_sites = 0

        for line in f:
            p = line.rstrip("\n").split("\t")
            dup_id = p[3]
            gts = p[5:]
            n_total_sites += 1

            a = annot.get(dup_id, {})
            sc = a.get("size_class", "unknown")
            fc = a.get("func_class", "intergenic")
            rep = int(a.get("repeat_50pct", 0))

            for i, gt in enumerate(gts):
                s = samples[i]
                if gt in ("./.", ".", "./0", "0/."):
                    missing[s] += 1
                    continue
                if gt in ("0/0", "0|0"):
                    continue

                counters[s]["total"] += 1
                counters[s][f"size_{sc}"] += 1
                counters[s][f"func_{fc}"] += 1
                if rep:
                    counters[s]["repeat"] += 1
                else:
                    counters[s]["nonrepeat"] += 1

    with open(args.binary_gt) as f:
        bhead = f.readline().rstrip("\n").split("\t")
        bsamples = bhead[1:]
        rows = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append([int(x) for x in p[1:]])

    mat = np.array(rows, dtype=np.float64).T
    col_var = mat.var(axis=0)
    mat = mat[:, col_var > 0]
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]
    total_var = np.sum(S ** 2)
    var_pct = (S ** 2 / total_var * 100) if total_var > 0 else [0, 0]

    pc_dict = {}
    for i, s in enumerate(bsamples):
        pc_dict[s] = (round(pc1[i], 4), round(pc2[i], 4))

    out = os.path.join(args.outdir, "per_sample_DUP_summary.tsv")
    with open(out, "w") as o:
        cols = [
            "sample", "subset", "total_DUP",
            "burden_small", "burden_medium", "burden_large",
            "burden_repeat", "burden_nonrepeat",
            "burden_CDS", "burden_exon", "burden_intronic", "burden_intergenic",
            "missingness", "ancestry_cluster", "max_Q", "PC1", "PC2"
        ]
        o.write("\t".join(cols) + "\n")
        for s in samples:
            c = counters[s]
            miss = round(missing[s] / n_total_sites, 4) if n_total_sites else 0
            anc = ancestry.get(s, (".", "."))
            pc = pc_dict.get(s, (".", "."))
            subset = "unrelated_81" if s in unrel else "related"
            row = [
                s, subset, c["total"],
                c.get("size_small", 0), c.get("size_medium", 0), c.get("size_large", 0),
                c.get("repeat", 0), c.get("nonrepeat", 0),
                c.get("func_CDS_overlap", 0), c.get("func_exon_overlap", 0),
                c.get("func_intronic", 0), c.get("func_intergenic", 0),
                miss, anc[0], anc[1], pc[0], pc[1]
            ]
            o.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wrote {out}")
    if len(var_pct) >= 2:
        print(f"PCA: PC1={var_pct[0]:.1f}%, PC2={var_pct[1]:.1f}%")

if __name__ == "__main__":
    main()
