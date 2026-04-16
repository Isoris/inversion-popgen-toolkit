#!/usr/bin/env python3
"""
07_build_master_annotation.py — Build master DEL annotation table
Merges: VCF info, functional class, repeat overlap, depth support, mate-distance QC,
        carrier counts, frequency class, size class

Inputs (from 00_delly_config.sh paths):
  - catalog_226.DEL.vcf.gz          → VCF fields (QUAL, PE, SR, PRECISE, FILTER)
  - catalog_226.DEL.GT_matrix.tsv   → carrier counts + per-sample genotypes
  - catalog_226.functional_class.tsv → gene/exon/CDS overlap
  - catalog_226.DELs_in_repeats.bed → repeat status
  - depth_support_226.tsv           → depth ratio
  - mate_distance_qc_226.tsv        → SVLEN QC flag

Outputs:
  - 11_summary/master_DEL_annotation_226.tsv
  - 11_summary/master_DEL_annotation_81.tsv
  - 11_summary/frequency_class_summary.tsv
  - 11_summary/size_class_summary.tsv
"""
import gzip, os, sys, argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf_226", required=True)
    p.add_argument("--gt_matrix_226", required=True)
    p.add_argument("--functional_226", required=True)
    p.add_argument("--repeats_in_bed", required=True, help="DELs_in_repeats.bed")
    p.add_argument("--depth_support", required=True)
    p.add_argument("--mate_distance_qc", required=True)
    p.add_argument("--samples_unrelated", required=True)
    p.add_argument("--outdir", required=True)
    # Optional: ancestry
    p.add_argument("--qopt", default="", help="NGSadmix .qopt file for ancestry assignment")
    p.add_argument("--qopt_samples", default="", help="Sample list matching qopt row order")
    return p.parse_args()

# ── Frequency class ─────────────────────────────────────────────────────────
def freq_class(n_carriers, n_total):
    if n_carriers == 0:      return "absent"
    if n_carriers == 1:      return "singleton"
    if n_carriers <= 5:      return "rare"
    if n_carriers <= 20:     return "low_frequency"
    if n_carriers <= 80:     return "common"
    if n_carriers <= 200:    return "widespread"
    if n_carriers > 200:     return "near_fixed"
    return "unknown"

def size_class(length_bp):
    if length_bp < 50:       return "sub_SV"
    if length_bp <= 500:     return "small"
    if length_bp <= 5000:    return "medium"
    return "large"

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── 1. Load unrelated samples ───────────────────────────────────────────
    unrel = set()
    with open(args.samples_unrelated) as f:
        for line in f:
            s = line.strip()
            if s: unrel.add(s)
    print(f"Unrelated samples: {len(unrel)}")

    # ── 2. Load ancestry if available ───────────────────────────────────────
    ancestry = {}  # sample -> (cluster_id, max_q)
    if args.qopt and os.path.exists(args.qopt) and args.qopt_samples and os.path.exists(args.qopt_samples):
        with open(args.qopt_samples) as f:
            qsamp = [l.strip() for l in f if l.strip()]
        with open(args.qopt) as f:
            for i, line in enumerate(f):
                if i >= len(qsamp): break
                vals = [float(x) for x in line.strip().split()]
                maxq = max(vals)
                cluster = vals.index(maxq) + 1  # 1-indexed
                ancestry[qsamp[i]] = (f"Q{cluster}", round(maxq, 4))
        print(f"Ancestry loaded for {len(ancestry)} samples")

    # ── 3. Parse VCF for per-DEL INFO fields ────────────────────────────────
    print("Parsing VCF INFO fields...")
    vcf_info = {}  # del_id -> dict
    opener = gzip.open if args.vcf_226.endswith('.gz') else open
    with opener(args.vcf_226, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            chrom, pos, del_id, filt = p[0], p[1], p[2], p[6]
            info = dict(kv.split('=', 1) for kv in p[7].split(';') if '=' in kv)
            flags = set(kv for kv in p[7].split(';') if '=' not in kv)

            vcf_info[del_id] = {
                'chrom': chrom,
                'pos': int(pos),
                'end': int(info.get('END', pos)),
                'qual': float(p[5]) if p[5] != '.' else 0,
                'filter': filt,
                'pe': int(info.get('PE', 0)),
                'sr': int(info.get('SR', 0)),
                'precise': 1 if 'PRECISE' in flags else 0,
                'svlen': abs(int(info.get('SVLEN', 0))),
                'ct': info.get('CT', '.'),
            }
    print(f"  VCF DELs: {len(vcf_info)}")

    # ── 4. Parse GT matrix for carrier counts ───────────────────────────────
    print("Parsing GT matrix for carrier counts...")
    carrier_data = {}  # del_key (chrom:pos-end) -> {n_226, n_81, missing_226, samples_alt}
    gt_samples = []
    with open(args.gt_matrix_226) as f:
        header = f.readline().strip().split('\t')
        gt_samples = header[5:]
        n_total = len(gt_samples)

        for line in f:
            p = line.strip().split('\t')
            chrom, pos, end, del_id = p[0], int(p[1]), int(p[2]), p[3]
            gts = p[5:]

            n_alt_226 = 0
            n_alt_81 = 0
            n_miss_226 = 0
            for i, gt in enumerate(gts):
                if gt in ('./.', '.', './0', '0/.'):
                    n_miss_226 += 1
                elif gt not in ('0/0', '0|0'):
                    n_alt_226 += 1
                    if gt_samples[i] in unrel:
                        n_alt_81 += 1

            carrier_data[del_id] = {
                'n_carriers_226': n_alt_226,
                'n_carriers_81': n_alt_81,
                'n_missing_226': n_miss_226,
                'missing_frac_226': round(n_miss_226 / n_total, 4) if n_total > 0 else 0,
            }
    print(f"  GT entries: {len(carrier_data)}")

    # ── 5. Load functional class ────────────────────────────────────────────
    print("Loading functional class...")
    func = {}  # del_key -> dict
    with open(args.functional_226) as f:
        fheader = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            # del_key, chrom, start, end, id, svlen, n_genes, n_exons, n_cds, class, gene_names
            del_key = p[0]  # "chrom:start-end"
            del_id = p[4] if len(p) > 4 else '.'
            func[del_key] = {
                'n_genes': int(p[6]) if len(p) > 6 else 0,
                'n_exons': int(p[7]) if len(p) > 7 else 0,
                'n_cds': int(p[8]) if len(p) > 8 else 0,
                'func_class': p[9] if len(p) > 9 else 'intergenic',
                'gene_names': p[10] if len(p) > 10 else '.',
            }
            # Also index by del_id for join
            if del_id != '.':
                func[del_id] = func[del_key]
    print(f"  Functional entries: {len(func)}")

    # ── 6. Load repeat overlap ──────────────────────────────────────────────
    print("Loading repeat overlap...")
    in_repeat = set()
    if os.path.exists(args.repeats_in_bed):
        with open(args.repeats_in_bed) as f:
            for line in f:
                p = line.strip().split('\t')
                key = f"{p[0]}:{p[1]}-{p[2]}"
                in_repeat.add(key)
                if len(p) > 3:
                    in_repeat.add(p[3])
    print(f"  DELs with >=50% repeat overlap: {len(in_repeat)}")

    # ── 7. Load depth support ───────────────────────────────────────────────
    print("Loading depth support...")
    depth = {}
    if os.path.exists(args.depth_support):
        with open(args.depth_support) as f:
            f.readline()
            for line in f:
                p = line.strip().split('\t')
                depth[p[0]] = {
                    'depth_n_samples': p[1],
                    'depth_median_ratio': p[2],
                    'depth_label': p[4] if len(p) > 4 else '.',
                }
    print(f"  Depth entries: {len(depth)}")

    # ── 8. Load mate distance QC ────────────────────────────────────────────
    print("Loading mate distance QC...")
    mate_qc = {}
    if os.path.exists(args.mate_distance_qc):
        with open(args.mate_distance_qc) as f:
            f.readline()
            for line in f:
                p = line.strip().split('\t')
                key = f"{p[0]}:{p[1]}-{p[2]}"
                mate_qc[key] = p[6] if len(p) > 6 else '.'
                if len(p) > 3:
                    mate_qc[p[3]] = mate_qc[key]
    print(f"  Mate QC entries: {len(mate_qc)}")

    # ── 9. Build master table ───────────────────────────────────────────────
    print("Building master annotation table...")
    n_total_226 = len(gt_samples)
    n_total_81 = len(unrel)

    out_226 = os.path.join(args.outdir, "master_DEL_annotation_226.tsv")
    freq_counts = defaultdict(int)
    size_counts = defaultdict(int)

    header_cols = [
        "del_id", "chrom", "pos", "end", "svlen", "size_class",
        "qual", "filter", "pe", "sr", "precise",
        "n_carriers_226", "n_carriers_81", "missing_frac_226",
        "freq_class_226", "freq_class_81",
        "repeat_50pct", "func_class", "n_genes", "n_exons", "n_cds", "gene_names",
        "depth_median_ratio", "depth_label",
        "mate_qc_flag",
    ]

    rows = []
    for del_id, vi in vcf_info.items():
        length = vi['end'] - vi['pos']
        if length <= 0:
            length = vi['svlen'] if vi['svlen'] > 0 else abs(vi['end'] - vi['pos'])
        sc = size_class(length)
        size_counts[sc] += 1

        # Carrier data
        cd = carrier_data.get(del_id, {})
        n226 = cd.get('n_carriers_226', 0)
        n81 = cd.get('n_carriers_81', 0)
        mf = cd.get('missing_frac_226', 0)

        fc226 = freq_class(n226, n_total_226)
        fc81 = freq_class(n81, n_total_81)
        freq_counts[fc226] += 1

        # Functional
        del_key = f"{vi['chrom']}:{vi['pos']}-{vi['end']}"
        fd = func.get(del_id, func.get(del_key, {}))

        # Repeat
        rep = 1 if (del_id in in_repeat or del_key in in_repeat) else 0

        # Depth
        dd = depth.get(del_id, depth.get(del_key, {}))

        # Mate QC
        mq = mate_qc.get(del_id, mate_qc.get(del_key, '.'))

        row = [
            del_id, vi['chrom'], vi['pos'], vi['end'], length, sc,
            vi['qual'], vi['filter'], vi['pe'], vi['sr'], vi['precise'],
            n226, n81, mf,
            fc226, fc81,
            rep, fd.get('func_class', 'intergenic'),
            fd.get('n_genes', 0), fd.get('n_exons', 0), fd.get('n_cds', 0),
            fd.get('gene_names', '.'),
            dd.get('depth_median_ratio', 'NA'), dd.get('depth_label', '.'),
            mq,
        ]
        rows.append(row)

    # Sort by chrom, pos
    chr_order = {}
    for del_id, vi in vcf_info.items():
        c = vi['chrom']
        if c not in chr_order:
            chr_order[c] = len(chr_order)
    rows.sort(key=lambda r: (chr_order.get(r[1], 999), r[2]))

    with open(out_226, 'w') as o:
        o.write('\t'.join(header_cols) + '\n')
        for r in rows:
            o.write('\t'.join(str(x) for x in r) + '\n')
    print(f"  Master 226: {len(rows)} DELs → {out_226}")

    # ── 10. 81-only subset ──────────────────────────────────────────────────
    out_81 = os.path.join(args.outdir, "master_DEL_annotation_81.tsv")
    rows_81 = [r for r in rows if r[12] > 0]  # n_carriers_81 > 0
    with open(out_81, 'w') as o:
        o.write('\t'.join(header_cols) + '\n')
        for r in rows_81:
            o.write('\t'.join(str(x) for x in r) + '\n')
    print(f"  Master 81: {len(rows_81)} DELs → {out_81}")

    # ── 11. Summary tables ──────────────────────────────────────────────────
    freq_out = os.path.join(args.outdir, "frequency_class_summary.tsv")
    with open(freq_out, 'w') as o:
        o.write("frequency_class\tcount\n")
        for fc in ["singleton", "rare", "low_frequency", "common", "widespread", "near_fixed"]:
            o.write(f"{fc}\t{freq_counts.get(fc, 0)}\n")
    print(f"  Frequency summary → {freq_out}")

    size_out = os.path.join(args.outdir, "size_class_summary.tsv")
    with open(size_out, 'w') as o:
        o.write("size_class\tcount\n")
        for sc in ["sub_SV", "small", "medium", "large"]:
            o.write(f"{sc}\t{size_counts.get(sc, 0)}\n")
    print(f"  Size summary → {size_out}")

    print("Done.")

if __name__ == "__main__":
    main()
