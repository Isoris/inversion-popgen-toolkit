#!/usr/bin/env python3
import os
import argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--master_annot", required=True)
    p.add_argument("--ref_fai", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    chr_lengths = {}
    chr_order = []
    with open(args.ref_fai) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            chr_lengths[p[0]] = int(p[1])
            chr_order.append(p[0])

    dups = []
    with open(args.master_annot) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            dups.append(dict(zip(header, p)))

    gene_dup_rows = []
    for d in dups:
        genes = d.get("gene_names", ".")
        if genes == "." or not genes:
            continue
        for g in [x.strip() for x in genes.split(",") if x.strip() and x.strip() != "."]:
            gene_dup_rows.append({
                "dup_id": d["dup_id"],
                "gene": g,
                "func_class": d.get("func_class", "."),
                "chrom": d["chrom"],
                "pos": d["pos"],
                "end": d["end"],
                "span_bp": d.get("span_bp", "."),
                "n_carriers_226": d.get("n_carriers_226", "."),
                "repeat_50pct": d.get("repeat_50pct", "."),
                "qual": d.get("qual", "."),
            })

    out1 = os.path.join(args.outdir, "gene_DUP_overlap_table.tsv")
    with open(out1, "w") as o:
        cols = ["dup_id", "gene", "func_class", "chrom", "pos", "end", "span_bp", "n_carriers_226", "repeat_50pct", "qual"]
        o.write("\t".join(cols) + "\n")
        for r in gene_dup_rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    gene_info = defaultdict(lambda: {"n_dups": 0, "total_carriers": 0, "has_cds": False, "has_exon": False, "has_intronic": False, "chroms": set(), "repeat_count": 0})
    for r in gene_dup_rows:
        g = r["gene"]
        gene_info[g]["n_dups"] += 1
        try:
            gene_info[g]["total_carriers"] += int(r["n_carriers_226"])
        except:
            pass
        fc = r["func_class"]
        if fc == "CDS_overlap":
            gene_info[g]["has_cds"] = True
        elif fc == "exon_overlap":
            gene_info[g]["has_exon"] = True
        elif fc == "intronic":
            gene_info[g]["has_intronic"] = True
        gene_info[g]["chroms"].add(r["chrom"])
        try:
            if int(r["repeat_50pct"]) > 0:
                gene_info[g]["repeat_count"] += 1
        except:
            pass

    out2 = os.path.join(args.outdir, "gene_summary_collapsed.tsv")
    gene_rows = []
    for gene, info in gene_info.items():
        highest = "CDS_overlap" if info["has_cds"] else ("exon_overlap" if info["has_exon"] else "intronic")
        gene_rows.append({
            "gene": gene,
            "n_overlapping_DUPs": info["n_dups"],
            "total_carriers": info["total_carriers"],
            "highest_class": highest,
            "n_repeat_DUPs": info["repeat_count"],
            "chromosomes": ",".join(sorted(info["chroms"])),
        })
    gene_rows.sort(key=lambda x: (-x["n_overlapping_DUPs"], -x["total_carriers"]))

    with open(out2, "w") as o:
        cols = ["gene", "n_overlapping_DUPs", "total_carriers", "highest_class", "n_repeat_DUPs", "chromosomes"]
        o.write("\t".join(cols) + "\n")
        for r in gene_rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    out3 = os.path.join(args.outdir, "top_genes_report.tsv")
    cds_genes = [r for r in gene_rows if r["highest_class"] == "CDS_overlap"]
    top_by_count = gene_rows[:50]
    top_by_cds = cds_genes[:50]
    combined = {r["gene"]: r for r in top_by_count}
    combined.update({r["gene"]: r for r in top_by_cds})
    final_top = sorted(combined.values(), key=lambda x: (-x["n_overlapping_DUPs"]))

    with open(out3, "w") as o:
        cols = ["gene", "n_overlapping_DUPs", "total_carriers", "highest_class", "n_repeat_DUPs", "chromosomes"]
        o.write("\t".join(cols) + "\n")
        for r in final_top:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    chr_data = defaultdict(lambda: {"n_dup": 0, "sizes": [], "repeat": 0, "nonrepeat": 0, "cds": 0, "exon": 0, "intronic": 0, "intergenic": 0})
    for d in dups:
        c = d["chrom"]
        chr_data[c]["n_dup"] += 1
        try:
            chr_data[c]["sizes"].append(int(d["span_bp"]))
        except:
            pass
        try:
            if int(d.get("repeat_50pct", 0)) > 0:
                chr_data[c]["repeat"] += 1
            else:
                chr_data[c]["nonrepeat"] += 1
        except:
            chr_data[c]["nonrepeat"] += 1

        fc = d.get("func_class", "intergenic")
        if fc == "CDS_overlap":
            chr_data[c]["cds"] += 1
        elif fc == "exon_overlap":
            chr_data[c]["exon"] += 1
        elif fc == "intronic":
            chr_data[c]["intronic"] += 1
        else:
            chr_data[c]["intergenic"] += 1

    out4 = os.path.join(args.outdir, "per_chromosome_summary.tsv")
    with open(out4, "w") as o:
        o.write("\t".join(["chrom", "length_Mb", "n_DUP", "DUP_per_Mb", "median_span_bp", "repeat", "nonrepeat", "CDS_overlap", "exon_overlap", "intronic", "intergenic"]) + "\n")
        for c in chr_order:
            if c not in chr_data:
                continue
            cd = chr_data[c]
            length_mb = chr_lengths[c] / 1e6
            dup_per_mb = round(cd["n_dup"] / length_mb, 2) if length_mb > 0 else 0
            med = sorted(cd["sizes"])[len(cd["sizes"]) // 2] if cd["sizes"] else 0
            o.write("\t".join(str(x) for x in [
                c, round(length_mb, 2), cd["n_dup"], dup_per_mb, med,
                cd["repeat"], cd["nonrepeat"], cd["cds"], cd["exon"], cd["intronic"], cd["intergenic"]
            ]) + "\n")

    print("Wrote gene/chromosome summaries")

if __name__ == "__main__":
    main()
