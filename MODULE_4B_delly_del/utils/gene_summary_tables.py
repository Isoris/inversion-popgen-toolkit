#!/usr/bin/env python3
"""
10_gene_summary_tables.py — Gene-overlap and per-chromosome summary tables

Outputs:
  - 11_summary/gene_DEL_overlap_table.tsv      (one row per DEL-gene pair)
  - 11_summary/gene_summary_collapsed.tsv       (one row per gene)
  - 11_summary/per_chromosome_summary.tsv        (per-chr with functional breakdown)
  - 11_summary/top_genes_report.tsv              (ranked by impact)
"""
import os, argparse
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

    # Load chromosome lengths
    chr_lengths = {}
    chr_order = []
    with open(args.ref_fai) as f:
        for line in f:
            p = line.strip().split('\t')
            chr_lengths[p[0]] = int(p[1])
            chr_order.append(p[0])

    # Load master annotation
    dels = []
    with open(args.master_annot) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            dels.append(d)

    # ── 1. Gene-DEL overlap table ───────────────────────────────────────────
    gene_del_rows = []
    for d in dels:
        genes_str = d.get('gene_names', '.')
        if genes_str == '.' or not genes_str:
            continue
        genes = [g.strip() for g in genes_str.split(',') if g.strip() and g.strip() != '.']
        for gene in genes:
            gene_del_rows.append({
                'del_id': d['del_id'],
                'gene': gene,
                'func_class': d.get('func_class', '.'),
                'chrom': d['chrom'],
                'pos': d['pos'],
                'end': d['end'],
                'svlen': d.get('svlen', '.'),
                'n_carriers_226': d.get('n_carriers_226', '.'),
                'repeat_50pct': d.get('repeat_50pct', '.'),
                'qual': d.get('qual', '.'),
            })

    out_gd = os.path.join(args.outdir, "gene_DEL_overlap_table.tsv")
    with open(out_gd, 'w') as o:
        cols = ['del_id', 'gene', 'func_class', 'chrom', 'pos', 'end', 'svlen',
                'n_carriers_226', 'repeat_50pct', 'qual']
        o.write('\t'.join(cols) + '\n')
        for r in gene_del_rows:
            o.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f"Gene-DEL overlap: {len(gene_del_rows)} pairs → {out_gd}")

    # ── 2. Gene summary (collapsed) ────────────────────────────────────────
    gene_info = defaultdict(lambda: {
        'n_dels': 0, 'total_carriers': 0,
        'has_cds': False, 'has_exon': False, 'has_intronic': False,
        'chroms': set(), 'repeat_count': 0,
    })
    for r in gene_del_rows:
        g = r['gene']
        gene_info[g]['n_dels'] += 1
        try:
            gene_info[g]['total_carriers'] += int(r['n_carriers_226'])
        except:
            pass
        fc = r['func_class']
        if fc == 'CDS_overlap': gene_info[g]['has_cds'] = True
        elif fc == 'exon_overlap': gene_info[g]['has_exon'] = True
        elif fc == 'intronic': gene_info[g]['has_intronic'] = True
        gene_info[g]['chroms'].add(r['chrom'])
        try:
            if int(r['repeat_50pct']) > 0:
                gene_info[g]['repeat_count'] += 1
        except: pass

    out_gs = os.path.join(args.outdir, "gene_summary_collapsed.tsv")
    gene_rows = []
    for gene, info in gene_info.items():
        highest = 'CDS_overlap' if info['has_cds'] else 'exon_overlap' if info['has_exon'] else 'intronic'
        gene_rows.append({
            'gene': gene,
            'n_overlapping_DELs': info['n_dels'],
            'total_carriers': info['total_carriers'],
            'highest_class': highest,
            'n_repeat_DELs': info['repeat_count'],
            'chromosomes': ','.join(sorted(info['chroms'])),
        })
    gene_rows.sort(key=lambda x: (-x['n_overlapping_DELs'], -x['total_carriers']))

    with open(out_gs, 'w') as o:
        cols = ['gene', 'n_overlapping_DELs', 'total_carriers', 'highest_class',
                'n_repeat_DELs', 'chromosomes']
        o.write('\t'.join(cols) + '\n')
        for r in gene_rows:
            o.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f"Gene summary: {len(gene_rows)} genes → {out_gs}")

    # ── 3. Top genes report ─────────────────────────────────────────────────
    out_top = os.path.join(args.outdir, "top_genes_report.tsv")
    # Top 50 by DEL count, then top 50 by CDS overlap
    cds_genes = [r for r in gene_rows if r['highest_class'] == 'CDS_overlap']
    top_by_count = gene_rows[:50]
    top_by_cds = cds_genes[:50]
    combined = {r['gene']: r for r in top_by_count}
    combined.update({r['gene']: r for r in top_by_cds})
    final_top = sorted(combined.values(), key=lambda x: (-x['n_overlapping_DELs']))

    with open(out_top, 'w') as o:
        cols = ['gene', 'n_overlapping_DELs', 'total_carriers', 'highest_class',
                'n_repeat_DELs', 'chromosomes']
        o.write('\t'.join(cols) + '\n')
        for r in final_top:
            o.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f"Top genes: {len(final_top)} → {out_top}")

    # ── 4. Per-chromosome summary ───────────────────────────────────────────
    chr_data = defaultdict(lambda: {
        'n_del': 0, 'sizes': [], 'repeat': 0, 'nonrepeat': 0,
        'cds': 0, 'exon': 0, 'intronic': 0, 'intergenic': 0,
    })
    for d in dels:
        c = d['chrom']
        chr_data[c]['n_del'] += 1
        try: chr_data[c]['sizes'].append(int(d['svlen']))
        except: pass
        try:
            if int(d.get('repeat_50pct', 0)) > 0:
                chr_data[c]['repeat'] += 1
            else:
                chr_data[c]['nonrepeat'] += 1
        except:
            chr_data[c]['nonrepeat'] += 1
        fc = d.get('func_class', 'intergenic')
        if fc == 'CDS_overlap': chr_data[c]['cds'] += 1
        elif fc == 'exon_overlap': chr_data[c]['exon'] += 1
        elif fc == 'intronic': chr_data[c]['intronic'] += 1
        else: chr_data[c]['intergenic'] += 1

    out_chr = os.path.join(args.outdir, "per_chromosome_summary.tsv")
    with open(out_chr, 'w') as o:
        o.write('\t'.join([
            'chrom', 'length_Mb', 'n_DEL', 'DEL_per_Mb',
            'median_svlen', 'repeat', 'nonrepeat',
            'CDS_overlap', 'exon_overlap', 'intronic', 'intergenic',
        ]) + '\n')
        for c in chr_order:
            if c not in chr_data: continue
            cd = chr_data[c]
            length_mb = chr_lengths[c] / 1e6
            del_per_mb = round(cd['n_del'] / length_mb, 2) if length_mb > 0 else 0
            med = sorted(cd['sizes'])[len(cd['sizes'])//2] if cd['sizes'] else 0
            o.write('\t'.join(str(x) for x in [
                c, round(length_mb, 2), cd['n_del'], del_per_mb,
                med, cd['repeat'], cd['nonrepeat'],
                cd['cds'], cd['exon'], cd['intronic'], cd['intergenic'],
            ]) + '\n')
    print(f"Per-chromosome summary → {out_chr}")

if __name__ == "__main__":
    main()
