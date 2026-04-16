#!/usr/bin/env python3
"""
06_marker_selection.py — Select tiered variant markers (DELLY 09 equivalent)

Tier 1: Strict high-quality (low missing, common, well-genotyped)
Tier 2: Gene-overlapping candidates (if annotation available)
Tier 3: Group-informative (ancestry-discriminating)

Output:
  markers/markers_tier1_strict.{TYPE}.tsv
  markers/markers_tier3_group_informative.{TYPE}.tsv
  markers/marker_selection_summary.tsv
"""

import os, argparse, time
from collections import defaultdict, Counter

def _now(): return time.time()
def _elapsed(t0):
    dt = time.time() - t0
    return f"{dt:.1f}s" if dt < 60 else f"{dt/60:.1f}min"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--annotation_dir", required=True)
    p.add_argument("--matrices_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--per_sample_dir", default="")
    # Tier 1 thresholds
    p.add_argument("--min_carriers", type=int, default=5)
    p.add_argument("--max_missing_frac", type=float, default=0.15)
    p.add_argument("--max_carrier_frac", type=float, default=0.95)
    # Tier 3
    p.add_argument("--min_delta_freq", type=float, default=0.3)
    return p.parse_args()

def main():
    args = parse_args()
    t0 = _now()
    outdir = os.path.join(args.outdir, "markers")
    os.makedirs(outdir, exist_ok=True)

    # Load sample group assignments for tier 3
    sample_group = {}
    for vtype in ["SNP", "INDEL", "COMBINED"]:
        ps_path = os.path.join(args.per_sample_dir, f"per_sample_{vtype}_summary.tsv")
        if os.path.isfile(ps_path):
            with open(ps_path) as f:
                header = f.readline().strip().split('\t')
                col = {h: i for i, h in enumerate(header)}
                for line in f:
                    p = line.strip().split('\t')
                    s = p[col["SAMPLE"]]
                    anc = p[col.get("ANCESTRY_CLUSTER", -1)] if "ANCESTRY_CLUSTER" in col else "."
                    if anc != ".":
                        sample_group[s] = anc
            break  # only need one

    summary_rows = []

    for vtype in ["SNP", "INDEL", "COMBINED"]:
        annot_path = os.path.join(args.annotation_dir, f"master_{vtype}_annotation.tsv")
        gt_path = os.path.join(args.matrices_dir, f"GT_matrix.{vtype}.tsv")

        if not os.path.isfile(annot_path):
            continue

        print(f"[06] === {vtype} ===")

        # Load annotation
        variants = []
        with open(annot_path) as f:
            header = f.readline().strip().split('\t')
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                p = line.strip().split('\t')
                d = dict(zip(header, p))
                # Type conversions
                for k in ["N_CARRIERS_ALL", "N_CARRIERS_UNREL", "N_PHASED_SAMPLES"]:
                    try:
                        d[k] = int(d[k])
                    except:
                        d[k] = 0
                for k in ["CARRIER_FREQ_ALL", "CARRIER_FREQ_UNREL", "MISSING_FRAC"]:
                    try:
                        d[k] = float(d[k])
                    except:
                        d[k] = 0.0
                variants.append(d)

        n_total_samples = 226  # will be refined from GT matrix

        # ── TIER 1: Strict ──
        tier1 = []
        for d in variants:
            if d["N_CARRIERS_ALL"] < args.min_carriers:
                continue
            if d["MISSING_FRAC"] > args.max_missing_frac:
                continue
            if d["CARRIER_FREQ_ALL"] > args.max_carrier_frac:
                continue
            fc = d.get("FREQ_CLASS", "")
            if fc in ("monomorphic", "singleton"):
                continue
            tier1.append(d)

        tier1.sort(key=lambda x: -x["N_CARRIERS_ALL"])

        out1 = os.path.join(outdir, f"markers_tier1_strict.{vtype}.tsv")
        with open(out1, 'w') as o:
            o.write("\t".join(header) + "\n")
            for d in tier1:
                o.write("\t".join(str(d.get(h, ".")) for h in header) + "\n")
        print(f"[06]   Tier 1: {len(tier1)} → {out1}")

        # ── TIER 3: Group-informative ──
        tier3 = []
        if sample_group and os.path.isfile(gt_path):
            groups = defaultdict(list)
            for s, g in sample_group.items():
                groups[g].append(s)
            group_names = sorted(groups.keys())

            if len(group_names) >= 2:
                # Load GT matrix for group freq calculation
                with open(gt_path) as f:
                    gtheader = f.readline().strip().split('\t')
                    gt_samples = gtheader[5:]
                    n_total_samples = len(gt_samples)

                    gt_data = {}
                    for line in f:
                        p = line.strip().split('\t')
                        vk = p[4]
                        gts = p[5:]
                        has_alt = {}
                        for i, gt in enumerate(gts):
                            gt_clean = gt.replace("|", "/")
                            alleles = gt_clean.split("/")
                            has_alt[gt_samples[i]] = any(a not in ("0", ".", "") for a in alleles)
                        gt_data[vk] = has_alt

                for d in variants:
                    if d["MISSING_FRAC"] > 0.20:
                        continue
                    vk = d["VAR_KEY"]
                    if vk not in gt_data:
                        continue
                    gt = gt_data[vk]

                    group_freqs = {}
                    for g in group_names:
                        members = groups[g]
                        n_alt = sum(1 for s in members if gt.get(s, False))
                        group_freqs[g] = n_alt / len(members) if members else 0

                    max_delta = 0
                    best_pair = ("", "")
                    for i, g1 in enumerate(group_names):
                        for g2 in group_names[i + 1:]:
                            delta = abs(group_freqs[g1] - group_freqs[g2])
                            if delta > max_delta:
                                max_delta = delta
                                best_pair = (g1, g2)

                    if max_delta >= args.min_delta_freq:
                        d["MAX_DELTA_FREQ"] = round(max_delta, 4)
                        d["DELTA_PAIR"] = f"{best_pair[0]}_vs_{best_pair[1]}"
                        d["GROUP_FREQS"] = ";".join(f"{g}:{group_freqs[g]:.3f}" for g in group_names)
                        tier3.append(d)

        tier3.sort(key=lambda x: -x.get("MAX_DELTA_FREQ", 0))
        out3 = os.path.join(outdir, f"markers_tier3_group_informative.{vtype}.tsv")
        extra = ["MAX_DELTA_FREQ", "DELTA_PAIR", "GROUP_FREQS"]
        with open(out3, 'w') as o:
            o.write("\t".join(header + extra) + "\n")
            for d in tier3:
                o.write("\t".join(str(d.get(h, ".")) for h in header + extra) + "\n")
        print(f"[06]   Tier 3: {len(tier3)} → {out3}")

        summary_rows.append({"TYPE": vtype, "TIER1": len(tier1), "TIER3": len(tier3)})

    # Summary
    sum_path = os.path.join(outdir, "marker_selection_summary.tsv")
    with open(sum_path, 'w') as o:
        o.write("TYPE\tTIER1_STRICT\tTIER3_GROUP_INFORMATIVE\n")
        for r in summary_rows:
            o.write(f"{r['TYPE']}\t{r['TIER1']}\t{r['TIER3']}\n")
    print(f"[06] Summary → {sum_path}")
    print(f"[06] Total time: {_elapsed(t0)}")

if __name__ == "__main__":
    main()
