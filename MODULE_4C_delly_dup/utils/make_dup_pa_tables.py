#!/usr/bin/env python3
import gzip
import os

VCF = "07_final_catalogs/catalog_226.DUP.vcf.gz"
OUTDIR = "11_summary"
PREFIX = os.path.join(OUTDIR, "catalog_226.DUP")

def parse_info(info_str):
    d = {}
    for x in info_str.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
        else:
            d[x] = True
    return d

def gt_to_binary(gt):
    if gt in {"./.", ".", "./0", "0/.", ".|.", "./1", "1/.", ".|1", "1|.", ".|0", "0|."}:
        return "NA"
    if gt in {"0/0", "0|0"}:
        return "0"
    if gt in {"0/1", "1/0", "0|1", "1|0", "1/1", "1|1"}:
        return "1"
    return "NA"

def gt_to_dosage(gt):
    if gt in {"./.", ".", "./0", "0/.", ".|.", "./1", "1/.", ".|1", "1|.", ".|0", "0|."}:
        return "NA"
    if gt in {"0/0", "0|0"}:
        return "0"
    if gt in {"0/1", "1/0", "0|1", "1|0"}:
        return "1"
    if gt in {"1/1", "1|1"}:
        return "2"
    return "NA"

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    with gzip.open(VCF, "rt") as f, \
         open(PREFIX + ".PA_binary.tsv", "w") as out_pa, \
         open(PREFIX + ".GT_dosage.tsv", "w") as out_gt, \
         open(PREFIX + ".metadata.tsv", "w") as out_meta, \
         open(PREFIX + ".binary_gt_for_pca.tsv", "w") as out_pca:

        samples = None

        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                out_pa.write("dup_id\t" + "\t".join(samples) + "\n")
                out_gt.write("dup_id\t" + "\t".join(samples) + "\n")
                out_pca.write("dup_id\t" + "\t".join(samples) + "\n")
                out_meta.write("dup_id\tchrom\tpos\tend\tspan_bp\tqual\tfilter\tpe\tsr\tmapq\tn_carriers\tn_het\tn_homalt\tn_ref\tn_missing\n")
                continue

            p = line.rstrip("\n").split("\t")
            chrom, pos, dup_id, qual, flt, info_str, fmt = p[0], int(p[1]), p[2], p[5], p[6], p[7], p[8]
            sample_fields = p[9:]
            info = parse_info(info_str)
            end = int(info.get("END", pos))
            span = abs(end - pos)

            fmt_keys = fmt.split(":")
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else None

            pa_vals = []
            dosage_vals = []
            pca_vals = []

            n_carriers = 0
            n_het = 0
            n_homalt = 0
            n_ref = 0
            n_missing = 0

            for sf in sample_fields:
                vals = sf.split(":")
                gt = vals[gt_idx] if gt_idx is not None and gt_idx < len(vals) else "./."

                pa = gt_to_binary(gt)
                dosage = gt_to_dosage(gt)

                pa_vals.append(pa)
                dosage_vals.append(dosage)
                pca_vals.append("1" if pa == "1" else "0")

                if dosage == "NA":
                    n_missing += 1
                elif dosage == "0":
                    n_ref += 1
                elif dosage == "1":
                    n_carriers += 1
                    n_het += 1
                elif dosage == "2":
                    n_carriers += 1
                    n_homalt += 1

            out_pa.write(dup_id + "\t" + "\t".join(pa_vals) + "\n")
            out_gt.write(dup_id + "\t" + "\t".join(dosage_vals) + "\n")
            out_pca.write(dup_id + "\t" + "\t".join(pca_vals) + "\n")
            out_meta.write(
                f"{dup_id}\t{chrom}\t{pos}\t{end}\t{span}\t{qual}\t{flt}\t"
                f"{info.get('PE','.')}\t{info.get('SR','.')}\t{info.get('MAPQ','.')}\t"
                f"{n_carriers}\t{n_het}\t{n_homalt}\t{n_ref}\t{n_missing}\n"
            )

    print("Wrote DUP PA tables in", OUTDIR)

if __name__ == "__main__":
    main()
