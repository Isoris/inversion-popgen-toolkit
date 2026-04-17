#!/usr/bin/env python3
"""
05_delly_manta_concordance.py — v2: reads from unified candidate table
======================================================================
The concordance between DELLY2 and Manta is now computed in STEP01
(cross_caller_match / concordance_class columns). This script:

  1. Reads the unified matched candidate table
  2. Loads test results from STEP03 (Fisher P, OR, etc.)
  3. Produces a clean concordance report with evidence from both callers
  4. Generates manuscript-ready sentences

Evidence field handling:
  DELLY2: pe = INFO/PE, sr = INFO/SR (global counts)
  Manta:  pe = sum(FORMAT/PR alt), sr = sum(FORMAT/SR alt) (aggregated)
  Both stored as 'pe' and 'sr' in the unified table.
"""
import argparse
import os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True,
                   help="matched_inv_candidates.tsv from STEP01")
    p.add_argument("--tests", default="",
                   help="all_candidates_tests.tsv from STEP03")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load candidates
    cands = []
    with open(args.candidates) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            cands.append(dict(zip(header, line.strip().split('\t'))))
    print(f"Candidates: {len(cands)}")

    # Load test results if available
    tests = {}
    if args.tests and os.path.exists(args.tests):
        with open(args.tests) as f:
            theader = f.readline().strip().split('\t')
            for line in f:
                p = line.strip().split('\t')
                d = dict(zip(theader, p))
                tests[d['inv_id']] = d

    # Classify
    delly_only = [c for c in cands if c['caller'] == 'delly' and c['concordance_class'] != 'both_callers']
    manta_only = [c for c in cands if c['caller'] == 'manta' and c['concordance_class'] != 'both_callers']
    # For both_callers, only count DELLY side to avoid double-counting
    both_callers = [c for c in cands if c['caller'] == 'delly' and c['concordance_class'] == 'both_callers']

    print(f"BOTH_CALLERS: {len(both_callers)}")
    print(f"DELLY_ONLY:   {len(delly_only)}")
    print(f"MANTA_ONLY:   {len(manta_only)}")

    # Build concordance table
    conc_rows = []
    for c in cands:
        t = tests.get(c['inv_id'], {})
        conc_rows.append({
            'inv_id': c['inv_id'],
            'chrom': c['chrom'],
            'bp1': c['bp1_pos'],
            'bp2': c['bp2_pos'],
            'svlen': c['svlen_bp'],
            'caller': c['caller'],
            'qual': c['qual'],
            'pe': c['pe'],
            'sr': c['sr'],
            'precise': c['precise'],
            'concordance_class': c['concordance_class'],
            'cross_caller_match': c['cross_caller_match'],
            'ct_or_junction': c.get('ct', c.get('junction_type', '.')),
            'snake_match': c.get('snake_match', '.'),
            'fisher_OR': t.get('fisher_OR', '.'),
            'fisher_P': t.get('fisher_P', '.'),
            'CA_Z': t.get('CA_Z', '.'),
            'CA_P': t.get('CA_P', '.'),
            'inv_support_frac': t.get('inv_support_frac', '.'),
        })

    # Write concordance table
    conc_file = os.path.join(args.outdir, "delly_manta_concordance.tsv")
    if conc_rows:
        cols = list(conc_rows[0].keys())
        with open(conc_file, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in conc_rows:
                f.write('\t'.join(str(r.get(c, '.')) for c in cols) + '\n')
    print(f"Concordance: {conc_file}")

    # Summary
    summary_file = os.path.join(args.outdir, "concordance_summary.txt")
    n_delly_total = len(both_callers) + len(delly_only)
    n_manta_total = len(both_callers) + len(manta_only)
    n_both = len(both_callers)

    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("DELLY2 × Manta INV Concordance Summary (v2)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"  DELLY2 INV candidates:   {n_delly_total}\n")
        f.write(f"  Manta INV candidates:    {n_manta_total}\n")
        f.write(f"  BOTH_CALLERS:            {n_both}\n")
        f.write(f"  DELLY_ONLY:              {len(delly_only)}\n")
        f.write(f"  MANTA_ONLY:              {len(manta_only)}\n")
        f.write(f"  Unique events:           {n_both + len(delly_only) + len(manta_only)}\n\n")

        if n_delly_total > 0:
            f.write(f"  DELLY concordance rate:  {n_both}/{n_delly_total} ({n_both/n_delly_total*100:.0f}%)\n")
        if n_manta_total > 0:
            f.write(f"  Manta concordance rate:  {n_both}/{n_manta_total} ({n_both/n_manta_total*100:.0f}%)\n")

        # Both-callers with significant test
        both_sig = [c for c in both_callers
                    if c['inv_id'] in tests
                    and tests[c['inv_id']].get('fisher_P', '1') != '.'
                    and float(tests[c['inv_id']].get('fisher_P', '1')) < 0.05]

        f.write(f"\n  Both-callers + significant Fisher P<0.05: {len(both_sig)}\n")

        f.write("\n\nMANUSCRIPT SENTENCE:\n\n")
        # BUGFIX 2026-04-17 (chat 5, FIX 28): guard against zero Manta calls.
        # Previously the sentence unconditionally divided by n_manta_total
        # inside a `if n_delly_total > 0` block, crashing if DELLY had
        # calls but Manta had none.
        if n_delly_total > 0 and n_manta_total > 0:
            f.write(f"  Of {n_delly_total} DELLY2 inversions and {n_manta_total} Manta inversions "
                    f"passing quality filters (DELLY: QUAL≥{os.environ.get('MIN_DELLY_QUAL','300')}, "
                    f"PE≥{os.environ.get('MIN_DELLY_PE','3')}; Manta: QUAL≥{os.environ.get('MIN_MANTA_QUAL','50')}, "
                    f"total alt evidence≥{os.environ.get('MIN_MANTA_TOTAL_ALT','6')}), "
                    f"{n_both} ({n_both/n_delly_total*100:.0f}% of DELLY, "
                    f"{n_both/n_manta_total*100:.0f}% of Manta) were independently detected "
                    f"by both callers with ≥50% reciprocal overlap.")
            if both_sig:
                f.write(f" Among concordant calls, {len(both_sig)} showed statistically "
                        f"significant enrichment of breakpoint-supporting reads "
                        f"(Fisher's exact P < 0.05).")
            f.write("\n")
        elif n_delly_total > 0:
            f.write(f"  {n_delly_total} DELLY2 inversions passed quality filters; "
                    f"no Manta inversions found for comparison.\n")
        elif n_manta_total > 0:
            f.write(f"  {n_manta_total} Manta inversions passed quality filters; "
                    f"no DELLY2 inversions found for comparison.\n")
        else:
            f.write("  No inversions passed quality filters from either caller.\n")

    print(f"Summary: {summary_file}")


if __name__ == "__main__":
    main()
