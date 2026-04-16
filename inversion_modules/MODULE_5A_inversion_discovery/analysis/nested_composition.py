#!/usr/bin/env python3
"""
nested_composition.py — Internal ancestry structure classifier.

For each parent region (chromosome, candidate, inversion), finds overlapping
Engine B windows and characterizes the internal Q composition:
  - Dominant ancestry label blocks and switches
  - Fragmentation score
  - Internal entropy
  - Structure classification:
    homogeneous | dominant_plus_secondary | two_block_composite |
    continuous_gradient | multi_block_fragmented | diffuse_mixed

Consumes: Engine B per-sample local Q cache (*.local_Q_samples.tsv.gz)
Writes:   nested_composition.tsv (per-parent × per-sample)
          nested_composition_summary.tsv (per-parent, sample-averaged)

Usage:
  python3 nested_composition.py \
    --q_cache_dir local_Q/ \
    --parents candidates.tsv \
    --outdir nested_composition/ \
    --K 8

Origin: Adapted from MODULE_2B_Step2/helpers/nested_composition.py
"""

import os, sys, csv, gzip, math, argparse
from collections import defaultdict, Counter

def classify_structure(labels):
    """Classify internal Q structure from a sequence of assigned_pop labels."""
    if not labels:
        return "too_sparse"

    unique = set(labels)
    n_unique = len(unique)
    n = len(labels)

    if n < 3:
        return "too_sparse"

    switches = sum(1 for i in range(1, n) if labels[i] != labels[i-1])
    mc = Counter(labels).most_common()
    dom_frac = mc[0][1] / n

    if n_unique == 1:
        return "homogeneous"
    elif dom_frac >= 0.8 and switches <= 2:
        return "dominant_plus_secondary"
    elif n_unique == 2 and 0.3 <= dom_frac <= 0.7:
        return "two_block_composite"
    elif switches >= n * 0.4:
        return "multi_block_fragmented"
    elif switches <= 2 and n_unique <= 3:
        # Check monotonic gradient
        first_half_dom = Counter(labels[:n//2]).most_common(1)[0][0]
        second_half_dom = Counter(labels[n//2:]).most_common(1)[0][0]
        if first_half_dom != second_half_dom:
            return "continuous_gradient"
        return "dominant_plus_secondary"
    else:
        return "diffuse_mixed"


def block_analysis(labels):
    """Decompose label sequence into contiguous blocks."""
    if not labels:
        return [], 0

    blocks = []
    current = labels[0]
    length = 1
    for i in range(1, len(labels)):
        if labels[i] == current:
            length += 1
        else:
            blocks.append((current, length))
            current = labels[i]
            length = 1
    blocks.append((current, length))

    switches = len(blocks) - 1
    return blocks, switches


def load_q_samples(q_cache_dir, chrom):
    """Load per-window × per-sample Q from Engine B cache."""
    for ext in [".local_Q_samples.tsv.gz", ".local_Q_samples.tsv"]:
        path = os.path.join(q_cache_dir, chrom + ext)
        if os.path.isfile(path):
            opener = gzip.open if path.endswith(".gz") else open
            rows = []
            with opener(path, "rt") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    rows.append(row)
            return rows
    return []


def load_parents(path):
    """Load parent intervals (candidates, chromosomes, etc.)."""
    parents = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", row.get("chr", ""))
            start = int(row.get("start_bp", row.get("start", 0)))
            end = int(row.get("end_bp", row.get("end", 0)))
            pid = row.get("candidate_id", row.get("interval_id",
                  row.get("parent_id", f"{chrom}_{start}_{end}")))
            parents.append({"parent_id": pid, "chrom": chrom,
                           "start": start, "end": end})
    return parents


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q_cache_dir", required=True, help="Engine B cache directory")
    ap.add_argument("--parents", required=True, help="Parent intervals TSV")
    ap.add_argument("--outdir", default="nested_composition/")
    ap.add_argument("--K", type=int, default=8)
    ap.add_argument("--sample_list", default=None, help="Sample IDs file")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    parents = load_parents(args.parents)
    K = args.K

    # Load sample names
    sample_names = None
    if args.sample_list and os.path.isfile(args.sample_list):
        with open(args.sample_list) as f:
            sample_names = [l.strip() for l in f if l.strip()]

    per_sample_rows = []
    summary_rows = []

    for parent in parents:
        pid = parent["parent_id"]
        chrom = parent["chrom"]
        pstart = parent["start"]
        pend = parent["end"]

        q_rows = load_q_samples(args.q_cache_dir, chrom)
        if not q_rows:
            print(f"[SKIP] No Q data for {chrom}")
            continue

        # Filter to overlapping windows
        overlap = [r for r in q_rows
                   if int(r.get("start_bp", 0)) < pend
                   and int(r.get("end_bp", 0)) > pstart]

        if not overlap:
            continue

        # Group by sample
        by_sample = defaultdict(list)
        for row in overlap:
            sid = row.get("sample_id", row.get("sample_idx", ""))
            by_sample[sid].append(row)

        # Per-sample composition
        sample_summaries = []
        for sid, rows in by_sample.items():
            rows.sort(key=lambda r: int(r.get("start_bp", 0)))

            labels = [str(r.get("assigned_pop", "")) for r in rows]
            n_sub = len(rows)
            blocks, switches = block_analysis(labels)

            label_counts = Counter(labels)
            mc = label_counts.most_common()
            dom_label = mc[0][0]
            dom_frac = mc[0][1] / n_sub
            sec_label = mc[1][0] if len(mc) > 1 else ""
            sec_frac = mc[1][1] / n_sub if len(mc) > 1 else 0

            # Internal entropy
            probs = [c / n_sub for c in label_counts.values()]
            internal_H = -sum(p * math.log(p + 1e-15) for p in probs)

            frag_score = switches / max(1, n_sub - 1)
            struct_type = classify_structure(labels)

            # Largest blocks
            sorted_blocks = sorted(blocks, key=lambda x: -x[1])
            largest = sorted_blocks[0][1] if sorted_blocks else 0
            second = sorted_blocks[1][1] if len(sorted_blocks) > 1 else 0

            # Mean Q metrics from the rows
            mean_d12 = sum(float(r.get("delta12", 0)) for r in rows) / n_sub
            mean_H = sum(float(r.get("entropy", 0)) for r in rows) / n_sub
            mean_ena = sum(float(r.get("ena", 0)) for r in rows) / n_sub

            rec = {
                "parent_id": pid,
                "chrom": chrom,
                "parent_start": pstart,
                "parent_end": pend,
                "sample_id": sid,
                "n_windows": n_sub,
                "n_labels_unique": len(label_counts),
                "dominant_label": dom_label,
                "dominant_fraction": f"{dom_frac:.4f}",
                "second_label": sec_label,
                "second_fraction": f"{sec_frac:.4f}",
                "n_blocks": len(blocks),
                "n_switches": switches,
                "largest_block": largest,
                "second_block": second,
                "fragmentation_score": f"{frag_score:.4f}",
                "internal_entropy": f"{internal_H:.4f}",
                "structure_type": struct_type,
                "mean_delta12": f"{mean_d12:.4f}",
                "mean_entropy": f"{mean_H:.4f}",
                "mean_ena": f"{mean_ena:.4f}",
            }
            per_sample_rows.append(rec)
            sample_summaries.append(rec)

        # Summary across samples for this parent
        if sample_summaries:
            struct_counts = Counter(r["structure_type"] for r in sample_summaries)
            dom_struct = struct_counts.most_common(1)[0]

            mean_frag = sum(float(r["fragmentation_score"]) for r in sample_summaries) / len(sample_summaries)
            mean_int_H = sum(float(r["internal_entropy"]) for r in sample_summaries) / len(sample_summaries)
            mean_switches = sum(int(r["n_switches"]) for r in sample_summaries) / len(sample_summaries)

            summary_rows.append({
                "parent_id": pid,
                "chrom": chrom,
                "parent_start": pstart,
                "parent_end": pend,
                "n_samples": len(sample_summaries),
                "dominant_structure_type": dom_struct[0],
                "dominant_structure_count": dom_struct[1],
                "mean_fragmentation": f"{mean_frag:.4f}",
                "mean_internal_entropy": f"{mean_int_H:.4f}",
                "mean_switches": f"{mean_switches:.1f}",
                "structure_breakdown": "; ".join(f"{k}:{v}" for k, v in struct_counts.most_common()),
            })

    # Write per-sample
    if per_sample_rows:
        out1 = os.path.join(args.outdir, "nested_composition.tsv")
        fields = list(per_sample_rows[0].keys())
        with open(out1, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(per_sample_rows)
        print(f"[OK] {out1} ({len(per_sample_rows)} rows)")

    # Write summary
    if summary_rows:
        out2 = os.path.join(args.outdir, "nested_composition_summary.tsv")
        fields = list(summary_rows[0].keys())
        with open(out2, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(summary_rows)
        print(f"[OK] {out2} ({len(summary_rows)} parents)")

    print(f"\n[DONE] {len(per_sample_rows)} per-sample rows, {len(summary_rows)} parent summaries")


if __name__ == "__main__":
    main()
