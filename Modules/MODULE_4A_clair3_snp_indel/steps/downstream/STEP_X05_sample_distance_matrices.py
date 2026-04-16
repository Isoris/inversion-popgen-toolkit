#!/usr/bin/env python3
"""
05_sample_distance_matrices.py — Pairwise sharing + signature distances + rare network

Outputs:
  distances/pairwise_shared.{TYPE}.tsv          (sample_i × sample_j shared counts)
  distances/signature_distance_matrix.tsv        (Jaccard distance on phase block signatures)
  distances/rare_sharing_network.tsv             (edges: rare variants shared between sample pairs)
  distances/rare_sharing_nodes.tsv               (node attributes for network viz)
  distances/window_density.{TYPE}.tsv            (1-Mb window variant counts for heatmap)
"""

import os, sys, argparse, time
import numpy as np
from collections import defaultdict, Counter

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--matrices_dir", required=True)
    p.add_argument("--catalogs_dir", required=True)
    p.add_argument("--phase_prep_dir", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--ref_fai", default="")
    p.add_argument("--max_rare_carriers", type=int, default=10,
                   help="Max carriers for rare-sharing network [10]")
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "distances")
    os.makedirs(outdir, exist_ok=True)

    # Chromosome length for windowing
    chr_len = 0
    if args.ref_fai and os.path.isfile(args.ref_fai):
        with open(args.ref_fai) as f:
            for line in f:
                p = line.strip().split('\t')
                if p[0] == args.chrom:
                    chr_len = int(p[1])
                    break

    for vtype in ["SNP", "INDEL", "COMBINED"]:
        bin_path = os.path.join(args.matrices_dir, f"binary_genotype_matrix.{vtype}.tsv")
        cat_path = os.path.join(args.catalogs_dir, f"{vtype}_catalog.tsv")

        if not os.path.isfile(bin_path):
            continue

        print(f"[05] === {vtype} ===")

        # Load binary matrix
        with open(bin_path) as f:
            header = f.readline().strip().split('\t')
            samples = header[1:]
            ns = len(samples)
            var_keys = []
            rows = []
            for line in f:
                p = line.strip().split('\t')
                var_keys.append(p[0])
                rows.append([int(x) for x in p[1:]])

        mat = np.array(rows, dtype=np.uint8)  # (n_vars, n_samples)
        nv = mat.shape[0]
        print(f"[05]   Loaded: {nv} variants × {ns} samples")

        # ── Pairwise shared ──
        print(f"[05]   Computing pairwise sharing …")
        t1 = _now()
        shared = mat.T @ mat   # (ns, ns)

        pw_path = os.path.join(outdir, f"pairwise_shared.{vtype}.tsv")
        with open(pw_path, 'w') as o:
            o.write("sample_i\tsample_j\tn_shared\n")
            for i in range(ns):
                for j in range(i, ns):
                    o.write(f"{samples[i]}\t{samples[j]}\t{shared[i, j]}\n")
        print(f"[05]   Pairwise sharing → {pw_path} ({_elapsed(t1)})")

        # ── 1-Mb window density ──
        # Parse positions from var_keys
        window_counts = defaultdict(lambda: np.zeros(ns, dtype=int))
        for vi, vk in enumerate(var_keys):
            parts = vk.split(":")
            if len(parts) >= 2:
                pos = int(parts[1])
                w = (pos // 1000000) * 1000000
                window_counts[w] += mat[vi, :]

        wd_path = os.path.join(outdir, f"window_density.{vtype}.tsv")
        with open(wd_path, 'w') as o:
            o.write("\t".join(["CHROM", "WINDOW_START", "WINDOW_END"] + samples) + "\n")
            for w in sorted(window_counts.keys()):
                o.write("\t".join(
                    [args.chrom, str(w), str(w + 1000000)] +
                    [str(x) for x in window_counts[w]]
                ) + "\n")
        print(f"[05]   Window density → {wd_path}")

        # ── Rare sharing network ──
        print(f"[05]   Building rare sharing network …")
        t2 = _now()

        # Load carrier counts from catalog
        carrier_counts = {}
        if os.path.isfile(cat_path):
            with open(cat_path) as f:
                cheader = f.readline().strip().split('\t')
                ccol = {h: i for i, h in enumerate(cheader)}
                for line in f:
                    p = line.strip().split('\t')
                    vk = p[ccol["VAR_KEY"]]
                    nc = int(p[ccol["N_CARRIERS_ALL"]])
                    carrier_counts[vk] = nc

        # Build network edges from rare variants
        edge_counts = defaultdict(int)   # (si, sj) → n_shared_rare
        edge_variants = defaultdict(list)

        for vi, vk in enumerate(var_keys):
            nc = carrier_counts.get(vk, int(mat[vi, :].sum()))
            if nc < 2 or nc > args.max_rare_carriers:
                continue
            # Find carriers
            carriers = np.where(mat[vi, :] > 0)[0]
            for ii in range(len(carriers)):
                for jj in range(ii + 1, len(carriers)):
                    si, sj = carriers[ii], carriers[jj]
                    key = (min(si, sj), max(si, sj))
                    edge_counts[key] += 1
                    if len(edge_variants[key]) < 20:
                        edge_variants[key].append(vk)

        net_path = os.path.join(outdir, f"rare_sharing_network.{vtype}.tsv")
        with open(net_path, 'w') as o:
            o.write("sample_i\tsample_j\tn_shared_rare\texample_variants\n")
            for (si, sj), cnt in sorted(edge_counts.items(), key=lambda x: -x[1]):
                examples = ";".join(edge_variants[(si, sj)][:5])
                o.write(f"{samples[si]}\t{samples[sj]}\t{cnt}\t{examples}\n")
        print(f"[05]   Rare network: {len(edge_counts)} edges → {net_path} ({_elapsed(t2)})")

        # Node attributes
        node_path = os.path.join(outdir, f"rare_sharing_nodes.{vtype}.tsv")
        node_degree = Counter()
        node_total_shared = Counter()
        for (si, sj), cnt in edge_counts.items():
            node_degree[si] += 1
            node_degree[sj] += 1
            node_total_shared[si] += cnt
            node_total_shared[sj] += cnt

        with open(node_path, 'w') as o:
            o.write("sample\tn_variants_total\tn_rare_partners\ttotal_rare_shared\n")
            for i, s in enumerate(samples):
                n_total = int(mat[:, i].sum())
                o.write(f"{s}\t{n_total}\t{node_degree.get(i, 0)}\t{node_total_shared.get(i, 0)}\n")
        print(f"[05]   Rare network nodes → {node_path}")

    # ── Phase block signature distance ──
    sig_path = os.path.join(args.phase_prep_dir, "phase_block_variant_map.tsv") if args.phase_prep_dir else ""
    if sig_path and os.path.isfile(sig_path):
        print(f"\n[05] === Phase block signature distances ===")
        # Build: for each sample, the set of (pos, phase_gt) for TIER_1 variants
        sample_signatures = defaultdict(set)
        with open(sig_path) as f:
            header = f.readline().strip().split('\t')
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                p = line.strip().split('\t')
                if p[col["PHASE_TIER"]] == "TIER_1_WHATSHAP":
                    s = p[col["SAMPLE"]]
                    pos = p[col["POS"]]
                    pgt = p[col["PHASE_GT"]]
                    sample_signatures[s].add(f"{pos}:{pgt}")

        sig_samples = sorted(sample_signatures.keys())
        ns_sig = len(sig_samples)

        if ns_sig > 1:
            dist_path = os.path.join(outdir, "signature_distance_matrix.tsv")
            with open(dist_path, 'w') as o:
                o.write("\t".join([""] + sig_samples) + "\n")
                for i, si in enumerate(sig_samples):
                    dists = []
                    for j, sj in enumerate(sig_samples):
                        if i == j:
                            dists.append("0.0000")
                        else:
                            a = sample_signatures[si]
                            b = sample_signatures[sj]
                            union = len(a | b)
                            inter = len(a & b)
                            jd = round(1 - inter / union, 4) if union > 0 else 1.0
                            dists.append(str(jd))
                    o.write(si + "\t" + "\t".join(dists) + "\n")
            print(f"[05]   Signature distance: {ns_sig} samples → {dist_path}")

    print(f"\n[05] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
