#!/usr/bin/env python3
"""
09_select_markers.py — Select high-quality DEL markers from master annotation

Tier 1: Strict marker-grade (validation-ready)
Tier 2: Biologically interesting (gene-overlap candidates)
Tier 3: Population-informative (group-discriminating)

Input:
  - master_DEL_annotation_226.tsv

Output:
  - 11_summary/markers_tier1_strict.tsv
  - 11_summary/markers_tier2_gene.tsv
  - 11_summary/markers_tier3_group_informative.tsv
  - 11_summary/marker_selection_summary.tsv
"""
import os, argparse
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--master_annot", required=True)
    p.add_argument("--outdir", required=True)
    # Tier 1 thresholds
    p.add_argument("--min_qual", type=float, default=500)
    p.add_argument("--min_pe", type=int, default=5)
    p.add_argument("--max_missing_frac", type=float, default=0.10)
    p.add_argument("--max_repeat_frac", type=float, default=0.5)  # 0 or 1 in our binary
    p.add_argument("--min_carriers", type=int, default=2)
    p.add_argument("--max_carrier_frac", type=float, default=0.95)
    p.add_argument("--max_length", type=int, default=50000)
    # Ancestry info for tier 3
    p.add_argument("--per_sample_summary", default="", help="per_sample_DEL_summary.tsv")
    p.add_argument("--gt_matrix", default="", help="GT matrix for group freq calculation")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load master annotation
    print("Loading master annotation...")
    dels = []
    with open(args.master_annot) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            d = dict(zip(header, p))
            # Type conversions
            for k in ['pos', 'end', 'svlen', 'pe', 'sr', 'precise',
                       'n_carriers_226', 'n_carriers_81', 'n_genes', 'n_exons', 'n_cds']:
                try: d[k] = int(d[k])
                except: d[k] = 0
            for k in ['qual', 'missing_frac_226']:
                try: d[k] = float(d[k])
                except: d[k] = 0.0
            try: d['repeat_50pct'] = int(d['repeat_50pct'])
            except: d['repeat_50pct'] = 0
            dels.append(d)
    print(f"  Total DELs: {len(dels)}")

    # Load sample-to-group mapping for tier 3
    sample_group = {}
    if args.per_sample_summary and os.path.exists(args.per_sample_summary):
        with open(args.per_sample_summary) as f:
            sheader = f.readline().strip().split('\t')
            for line in f:
                p = line.strip().split('\t')
                sd = dict(zip(sheader, p))
                sample_group[sd['sample']] = sd.get('ancestry_cluster', '.')

    # Load GT matrix for group frequency calculation (tier 3)
    gt_data = {}  # del_id -> {sample: has_alt}
    gt_samples = []
    if args.gt_matrix and os.path.exists(args.gt_matrix):
        with open(args.gt_matrix) as f:
            gtheader = f.readline().strip().split('\t')
            gt_samples = gtheader[5:]
            for line in f:
                p = line.strip().split('\t')
                del_id = p[3]
                gts = p[5:]
                gt_data[del_id] = {}
                for i, gt in enumerate(gts):
                    has_alt = gt not in ('./.', '0/0', './0', '0/.', '.', '0|0')
                    gt_data[del_id][gt_samples[i]] = has_alt

    n_total = 226  # hardcoded for now

    # ── TIER 1: Strict marker-grade ─────────────────────────────────────────
    print("Selecting Tier 1 markers...")
    tier1 = []
    for d in dels:
        if d['filter'] != 'PASS': continue
        if d['precise'] != 1: continue
        if d['qual'] < args.min_qual: continue
        if d['pe'] < args.min_pe: continue
        if d['n_carriers_226'] < args.min_carriers: continue
        if d['n_carriers_226'] > args.max_carrier_frac * n_total: continue
        if d['missing_frac_226'] > args.max_missing_frac: continue
        if d['repeat_50pct'] >= 1: continue  # no repeat overlap >= 50%
        if d['svlen'] > args.max_length: continue
        if d['svlen'] < 50: continue
        # Exclude mate-distance extreme artifacts
        mate_qc = d.get('mate_qc_flag', 'normal')
        if mate_qc in ('extreme_artifact_candidate', 'very_suspicious'): continue
        tier1.append(d)

    tier1.sort(key=lambda x: -x['qual'])
    out1 = os.path.join(args.outdir, "markers_tier1_strict.tsv")
    write_marker_table(out1, tier1, header)
    print(f"  Tier 1: {len(tier1)} strict markers → {out1}")

    # ── TIER 2: Biologically interesting (gene overlap) ─────────────────────
    print("Selecting Tier 2 markers...")
    tier2 = []
    for d in dels:
        if d['func_class'] not in ('CDS_overlap', 'exon_overlap'): continue
        if d['qual'] < 200: continue  # relaxed but still supported
        if d['pe'] < 3: continue
        if d['n_carriers_226'] < 1: continue
        if d['gene_names'] == '.': continue
        tier2.append(d)

    tier2.sort(key=lambda x: (-x['n_cds'], -x['qual']))
    out2 = os.path.join(args.outdir, "markers_tier2_gene.tsv")
    write_marker_table(out2, tier2, header)
    print(f"  Tier 2: {len(tier2)} gene-overlap candidates → {out2}")

    # ── TIER 3: Population-informative ──────────────────────────────────────
    print("Selecting Tier 3 markers...")
    tier3 = []

    if gt_data and sample_group:
        groups = defaultdict(list)
        for s, g in sample_group.items():
            if g != '.':
                groups[g].append(s)

        group_names = sorted(groups.keys())
        print(f"  Groups for tier 3: {', '.join(f'{g}({len(groups[g])})' for g in group_names)}")

        for d in dels:
            if d['qual'] < 300: continue
            if d['pe'] < 3: continue
            if d['missing_frac_226'] > 0.15: continue

            del_id = d['del_id']
            if del_id not in gt_data: continue
            gt = gt_data[del_id]

            # Compute per-group frequency
            group_freqs = {}
            for g in group_names:
                members = groups[g]
                n_alt = sum(1 for s in members if gt.get(s, False))
                group_freqs[g] = n_alt / len(members) if members else 0

            # Find max delta between any two groups
            max_delta = 0
            best_pair = ('', '')
            for i, g1 in enumerate(group_names):
                for g2 in group_names[i+1:]:
                    delta = abs(group_freqs[g1] - group_freqs[g2])
                    if delta > max_delta:
                        max_delta = delta
                        best_pair = (g1, g2)

            if max_delta >= 0.3:  # at least 30% frequency difference
                d['max_delta_freq'] = round(max_delta, 4)
                d['delta_pair'] = f"{best_pair[0]}_vs_{best_pair[1]}"
                d['group_freqs'] = ';'.join(f"{g}:{group_freqs[g]:.3f}" for g in group_names)
                tier3.append(d)

    tier3.sort(key=lambda x: -x.get('max_delta_freq', 0))
    out3 = os.path.join(args.outdir, "markers_tier3_group_informative.tsv")
    extra_cols = ['max_delta_freq', 'delta_pair', 'group_freqs']
    write_marker_table(out3, tier3, header, extra_cols)
    print(f"  Tier 3: {len(tier3)} group-informative markers → {out3}")

    # ── Summary ─────────────────────────────────────────────────────────────
    summary = os.path.join(args.outdir, "marker_selection_summary.tsv")
    with open(summary, 'w') as o:
        o.write("tier\tcount\tdescription\n")
        o.write(f"tier1_strict\t{len(tier1)}\tPASS+PRECISE+QUAL>={args.min_qual}+PE>={args.min_pe}+non-repeat+low-missing\n")
        o.write(f"tier2_gene\t{len(tier2)}\tCDS/exon overlap, reasonably supported\n")
        o.write(f"tier3_group\t{len(tier3)}\tGroup frequency delta >= 0.3\n")
    print(f"  Summary → {summary}")

def write_marker_table(path, rows, all_cols, extra_cols=None):
    cols = all_cols if extra_cols is None else all_cols + extra_cols
    with open(path, 'w') as o:
        o.write('\t'.join(cols) + '\n')
        for d in rows:
            o.write('\t'.join(str(d.get(c, '.')) for c in cols) + '\n')

if __name__ == "__main__":
    main()
