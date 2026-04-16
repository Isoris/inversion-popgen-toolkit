#!/usr/bin/env python3
import gzip
import os
import argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf_226", required=True)
    p.add_argument("--gt_matrix_226", required=True)
    p.add_argument("--functional_226", required=True)
    p.add_argument("--repeats_in_bed", required=True)
    p.add_argument("--depth_support", required=True)
    p.add_argument("--mate_distance_qc", required=True)
    p.add_argument("--samples_unrelated", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--qopt", default="")
    p.add_argument("--qopt_samples", default="")
    return p.parse_args()

def parse_info(info_str):
    d = {}
    for x in info_str.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
        else:
            d[x] = True
    return d

def freq_class(n_carriers, n_total):
    if n_carriers == 0: return "absent"
    if n_carriers == 1: return "singleton"
    if n_carriers <= 5: return "rare"
    if n_carriers <= 20: return "low_frequency"
    if n_carriers <= 80: return "common"
    if n_carriers <= 200: return "widespread"
    return "near_fixed"

def size_class(span_bp):
    if span_bp < 50: return "sub_SV"
    if span_bp <= 500: return "small"
    if span_bp <= 5000: return "medium"
    return "large"

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

    vcf_info = {}
    with gzip.open(args.vcf_226, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            chrom, pos, dup_id, qual, flt = p[0], int(p[1]), p[2], p[5], p[6]
            info = parse_info(p[7])
            end = int(info.get("END", pos))
            span = abs(end - pos)
            vcf_info[dup_id] = {
                "chrom": chrom,
                "pos": pos,
                "end": end,
                "span_bp": span,
                "qual": float(qual) if qual != "." else 0.0,
                "filter": flt,
                "pe": int(info.get("PE", 0)),
                "sr": int(info.get("SR", 0)),
                "precise": 1 if "PRECISE" in info else 0,
                "ct": info.get("CT", "."),
            }

    carrier_data = {}
    gt_samples = []
    with open(args.gt_matrix_226) as f:
        header = f.readline().rstrip("\n").split("\t")
        gt_samples = header[5:]
        n_total = len(gt_samples)

        for line in f:
            p = line.rstrip("\n").split("\t")
            dup_id = p[3]
            gts = p[5:]

            n_alt_226 = 0
            n_alt_81 = 0
            n_miss_226 = 0

            for i, gt in enumerate(gts):
                if gt in ("./.", ".", "./0", "0/."):
                    n_miss_226 += 1
                elif gt not in ("0/0", "0|0"):
                    n_alt_226 += 1
                    if gt_samples[i] in unrel:
                        n_alt_81 += 1

            carrier_data[dup_id] = {
                "n_carriers_226": n_alt_226,
                "n_carriers_81": n_alt_81,
                "missing_frac_226": round(n_miss_226 / n_total, 4) if n_total else 0
            }

    func = {}
    with open(args.functional_226) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            dup_key = p[0]
            dup_id = p[4]
            rec = {
                "n_genes": int(p[6]),
                "n_exons": int(p[7]),
                "n_cds": int(p[8]),
                "func_class": p[9],
                "gene_names": p[10],
            }
            func[dup_key] = rec
            func[dup_id] = rec

    in_repeat = set()
    if os.path.exists(args.repeats_in_bed):
        with open(args.repeats_in_bed) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                key = f"{p[0]}:{p[1]}-{p[2]}"
                in_repeat.add(key)
                if len(p) > 3:
                    in_repeat.add(p[3])

    depth = {}
    if os.path.exists(args.depth_support):
        with open(args.depth_support) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                depth[p[0]] = {
                    "depth_n_samples": p[1],
                    "depth_median_ratio": p[2],
                    "depth_label": p[4]
                }

    mate_qc = {}
    if os.path.exists(args.mate_distance_qc):
        with open(args.mate_distance_qc) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                key = f"{p[0]}:{p[1]}-{p[2]}"
                mate_qc[key] = p[6]
                mate_qc[p[3]] = p[6]

    n_total_226 = len(gt_samples)
    n_total_81 = len(unrel)

    header_cols = [
        "dup_id", "chrom", "pos", "end", "span_bp", "size_class",
        "qual", "filter", "pe", "sr", "precise",
        "n_carriers_226", "n_carriers_81", "missing_frac_226",
        "freq_class_226", "freq_class_81",
        "repeat_50pct", "func_class", "n_genes", "n_exons", "n_cds", "gene_names",
        "depth_median_ratio", "depth_label", "mate_qc_flag"
    ]

    rows = []
    freq_counts = defaultdict(int)
    size_counts = defaultdict(int)

    for dup_id, vi in vcf_info.items():
        span = vi["span_bp"]
        sc = size_class(span)
        size_counts[sc] += 1

        cd = carrier_data.get(dup_id, {})
        n226 = cd.get("n_carriers_226", 0)
        n81 = cd.get("n_carriers_81", 0)
        mf = cd.get("missing_frac_226", 0)

        fc226 = freq_class(n226, n_total_226)
        fc81 = freq_class(n81, n_total_81)
        freq_counts[fc226] += 1

        dup_key = f"{vi['chrom']}:{vi['pos']}-{vi['end']}"
        fd = func.get(dup_id, func.get(dup_key, {}))
        rep = 1 if (dup_id in in_repeat or dup_key in in_repeat) else 0
        dd = depth.get(dup_id, depth.get(dup_key, {}))
        mq = mate_qc.get(dup_id, mate_qc.get(dup_key, "."))

        rows.append([
            dup_id, vi["chrom"], vi["pos"], vi["end"], span, sc,
            vi["qual"], vi["filter"], vi["pe"], vi["sr"], vi["precise"],
            n226, n81, mf,
            fc226, fc81,
            rep, fd.get("func_class", "intergenic"),
            fd.get("n_genes", 0), fd.get("n_exons", 0), fd.get("n_cds", 0),
            fd.get("gene_names", "."),
            dd.get("depth_median_ratio", "NA"), dd.get("depth_label", "."),
            mq
        ])

    chr_order = {}
    for dup_id, vi in vcf_info.items():
        if vi["chrom"] not in chr_order:
            chr_order[vi["chrom"]] = len(chr_order)

    rows.sort(key=lambda r: (chr_order.get(r[1], 999), r[2]))

    out_226 = os.path.join(args.outdir, "master_DUP_annotation_226.tsv")
    with open(out_226, "w") as o:
        o.write("\t".join(header_cols) + "\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")

    out_81 = os.path.join(args.outdir, "master_DUP_annotation_81.tsv")
    with open(out_81, "w") as o:
        o.write("\t".join(header_cols) + "\n")
        for r in rows:
            if int(r[12]) > 0:
                o.write("\t".join(str(x) for x in r) + "\n")

    freq_out = os.path.join(args.outdir, "frequency_class_summary.tsv")
    with open(freq_out, "w") as o:
        o.write("frequency_class\tcount\n")
        for fc in ["singleton", "rare", "low_frequency", "common", "widespread", "near_fixed"]:
            o.write(f"{fc}\t{freq_counts.get(fc, 0)}\n")

    size_out = os.path.join(args.outdir, "size_class_summary.tsv")
    with open(size_out, "w") as o:
        o.write("size_class\tcount\n")
        for sc in ["sub_SV", "small", "medium", "large"]:
            o.write(f"{sc}\t{size_counts.get(sc, 0)}\n")

    print("Wrote:", out_226, out_81, freq_out, size_out)

if __name__ == "__main__":
    main()
