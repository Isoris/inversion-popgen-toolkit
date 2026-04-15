#!/usr/bin/env python3
import os
import argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--master_annot", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min_qual", type=float, default=500)
    p.add_argument("--min_pe", type=int, default=5)
    p.add_argument("--max_missing_frac", type=float, default=0.10)
    p.add_argument("--min_carriers", type=int, default=2)
    p.add_argument("--max_carrier_frac", type=float, default=0.95)
    p.add_argument("--max_length", type=int, default=50000)
    p.add_argument("--per_sample_summary", default="")
    p.add_argument("--gt_matrix", default="")
    return p.parse_args()

def write_table(path, rows, cols):
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for d in rows:
            o.write("\t".join(str(d.get(c, ".")) for c in cols) + "\n")

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    dups = []
    with open(args.master_annot) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d = dict(zip(header, p))
            for k in ["pos", "end", "span_bp", "pe", "sr", "precise",
                      "n_carriers_226", "n_carriers_81", "n_genes", "n_exons", "n_cds"]:
                try:
                    d[k] = int(d[k])
                except:
                    d[k] = 0
            for k in ["qual", "missing_frac_226"]:
                try:
                    d[k] = float(d[k])
                except:
                    d[k] = 0.0
            try:
                d["repeat_50pct"] = int(d["repeat_50pct"])
            except:
                d["repeat_50pct"] = 0
            dups.append(d)

    sample_group = {}
    if args.per_sample_summary and os.path.exists(args.per_sample_summary):
        with open(args.per_sample_summary) as f:
            shead = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                sd = dict(zip(shead, p))
                sample_group[sd["sample"]] = sd.get("ancestry_cluster", ".")

    gt_data = {}
    gt_samples = []
    if args.gt_matrix and os.path.exists(args.gt_matrix):
        with open(args.gt_matrix) as f:
            gthead = f.readline().rstrip("\n").split("\t")
            gt_samples = gthead[5:]
            for line in f:
                p = line.rstrip("\n").split("\t")
                dup_id = p[3]
                gts = p[5:]
                gt_data[dup_id] = {}
                for i, gt in enumerate(gts):
                    gt_data[dup_id][gt_samples[i]] = gt not in ("./.", ".", "./0", "0/.", "0/0", "0|0")

    n_total = 226

    tier1 = []
    for d in dups:
        if d["filter"] != "PASS": continue
        if d["precise"] != 1: continue
        if d["qual"] < args.min_qual: continue
        if d["pe"] < args.min_pe: continue
        if d["n_carriers_226"] < args.min_carriers: continue
        if d["n_carriers_226"] > args.max_carrier_frac * n_total: continue
        if d["missing_frac_226"] > args.max_missing_frac: continue
        if d["repeat_50pct"] >= 1: continue
        if d["span_bp"] > args.max_length: continue
        if d["span_bp"] < 50: continue
        mq = d.get("mate_qc_flag", "normal")
        if mq in ("extreme_large",):
            continue
        tier1.append(d)

    tier1.sort(key=lambda x: -x["qual"])

    tier2 = []
    for d in dups:
        if d["func_class"] not in ("CDS_overlap", "exon_overlap"): continue
        if d["qual"] < 200: continue
        if d["pe"] < 3: continue
        if d["n_carriers_226"] < 1: continue
        if d["gene_names"] == ".": continue
        tier2.append(d)

    tier2.sort(key=lambda x: (-x["n_cds"], -x["qual"]))

    tier3 = []
    if gt_data and sample_group:
        groups = defaultdict(list)
        for s, g in sample_group.items():
            if g != ".":
                groups[g].append(s)

        group_names = sorted(groups.keys())
        for d in dups:
            if d["qual"] < 300: continue
            if d["pe"] < 3: continue
            if d["missing_frac_226"] > 0.15: continue

            dup_id = d["dup_id"]
            if dup_id not in gt_data:
                continue

            group_freqs = {}
            for g in group_names:
                mem = groups[g]
                n_alt = sum(1 for s in mem if gt_data[dup_id].get(s, False))
                group_freqs[g] = n_alt / len(mem) if mem else 0

            max_delta = 0
            best_pair = ("", "")
            for i, g1 in enumerate(group_names):
                for g2 in group_names[i+1:]:
                    delta = abs(group_freqs[g1] - group_freqs[g2])
                    if delta > max_delta:
                        max_delta = delta
                        best_pair = (g1, g2)

            if max_delta >= 0.3:
                d["max_delta_freq"] = round(max_delta, 4)
                d["delta_pair"] = f"{best_pair[0]}_vs_{best_pair[1]}"
                d["group_freqs"] = ";".join(f"{g}:{group_freqs[g]:.3f}" for g in group_names)
                tier3.append(d)

    tier3.sort(key=lambda x: -x.get("max_delta_freq", 0))

    base_cols = header
    write_table(os.path.join(args.outdir, "markers_tier1_strict.tsv"), tier1, base_cols)
    write_table(os.path.join(args.outdir, "markers_tier2_gene.tsv"), tier2, base_cols)
    write_table(os.path.join(args.outdir, "markers_tier3_group_informative.tsv"), tier3, base_cols + ["max_delta_freq", "delta_pair", "group_freqs"])

    summary = os.path.join(args.outdir, "marker_selection_summary.tsv")
    with open(summary, "w") as o:
        o.write("tier\tcount\tdescription\n")
        o.write(f"tier1_strict\t{len(tier1)}\tPASS+PRECISE+QUAL>={args.min_qual}+PE>={args.min_pe}+non-repeat+low-missing\n")
        o.write(f"tier2_gene\t{len(tier2)}\tCDS/exon overlap, reasonably supported\n")
        o.write(f"tier3_group\t{len(tier3)}\tGroup frequency delta >= 0.3\n")

    print("Wrote marker tables")

if __name__ == "__main__":
    main()
