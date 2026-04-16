#!/usr/bin/env python3
"""
04_validation_plots.py — Publication-quality breakpoint validation figures

For each candidate with a test result, generates:
  - Panel A: Support proportion barplot (REF blue / HET grey / INV red)
  - Panel B: Per-sample evidence score strip
  - Panel C: Evidence breakdown heatmap
  - Panel D: Support frequency trend with CA p-value
  - Panel E: Statistical results text

Also generates a genome-wide summary figure.

Usage:
  python3 04_validation_plots.py \\
    --evidence_dir 01_per_sample_evidence \\
    --stats_dir 03_statistical_tests \\
    --plot_dir 05_plots
"""
import argparse, glob, os, sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--evidence_dir", required=True)
    p.add_argument("--stats_dir", required=True)
    p.add_argument("--plot_dir", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("ERROR: matplotlib/numpy required for plotting", file=sys.stderr)
        sys.exit(1)

    GROUP_COLORS = {'REF': '#4169E1', 'HET': '#A0A0A0', 'INV': '#DC143C'}
    GROUP_ORDER = ['REF', 'HET', 'INV']

    ev_files = sorted(glob.glob(os.path.join(args.evidence_dir, "*_evidence.tsv")))
    test_files = {os.path.basename(f).replace("_tests.tsv", ""): f
                  for f in glob.glob(os.path.join(args.stats_dir, "*_tests.tsv"))}

    print(f"Evidence files: {len(ev_files)}, Test files: {len(test_files)}")

    for ev_file in ev_files:
        inv_id = os.path.basename(ev_file).replace("_evidence.tsv", "")

        # Load evidence
        rows = []
        with open(ev_file) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                p = line.strip().split('\t')
                rows.append(dict(zip(header, p)))
        if not rows:
            continue

        # Load test results
        tests = {}
        if inv_id in test_files:
            with open(test_files[inv_id]) as f:
                th = f.readline().strip().split('\t')
                tp = f.readline().strip().split('\t')
                tests = dict(zip(th, tp))

        # Build counts
        counts = {g: {'yes': 0, 'no': 0} for g in GROUP_ORDER}
        for r in rows:
            g = r.get('group', 'UNKNOWN')
            s = r.get('support', 'no')
            if g in counts:
                counts[g][s] += 1

        totals = {g: counts[g]['yes'] + counts[g]['no'] for g in GROUP_ORDER}
        if all(t == 0 for t in totals.values()):
            continue

        # ── Figure ───────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 4, figsize=(18, 5),
                                 gridspec_kw={'width_ratios': [1, 1, 1.2, 1.5]})
        fig.suptitle(f"Breakpoint Validation: {inv_id}", fontsize=13, fontweight='bold')

        # Panel A: Support proportion
        ax = axes[0]
        x = np.arange(3)
        yes_vals = [counts[g]['yes'] for g in GROUP_ORDER]
        no_vals = [counts[g]['no'] for g in GROUP_ORDER]
        ax.bar(x, yes_vals, color=[GROUP_COLORS[g] for g in GROUP_ORDER],
               edgecolor='black', linewidth=0.5, label='Support')
        ax.bar(x, no_vals, bottom=yes_vals, color='white',
               edgecolor='black', linewidth=0.5, hatch='///')
        for i, g in enumerate(GROUP_ORDER):
            t = totals[g]
            if t > 0:
                frac = counts[g]['yes'] / t * 100
                ax.text(i, t + 0.3, f"{frac:.0f}%", ha='center', fontsize=10,
                        fontweight='bold', color=GROUP_COLORS[g])
        ax.set_xticks(x)
        ax.set_xticklabels([f"{g}\n(n={totals[g]})" for g in GROUP_ORDER])
        ax.set_ylabel("Samples")
        ax.set_title("A  Support by genotype", fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Panel B: Score strip
        ax = axes[1]
        for i, g in enumerate(GROUP_ORDER):
            scores = [int(r.get('total_score', 0)) for r in rows if r.get('group') == g]
            if scores:
                jitter = np.random.uniform(-0.15, 0.15, len(scores))
                ax.scatter(np.full(len(scores), i) + jitter, scores,
                          c=GROUP_COLORS[g], s=30, alpha=0.7, edgecolors='black', linewidths=0.3)
                med = sorted(scores)[len(scores)//2]
                ax.hlines(med, i-0.25, i+0.25, color='black', linewidth=2)
        ax.set_xticks(range(3))
        ax.set_xticklabels(GROUP_ORDER)
        ax.set_ylabel("Evidence score")
        ax.set_title("B  Per-sample score", fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Panel C: Trend
        ax = axes[2]
        fracs = [counts[g]['yes']/totals[g]*100 if totals[g] > 0 else 0 for g in GROUP_ORDER]
        ax.bar(range(3), fracs, color=[GROUP_COLORS[g] for g in GROUP_ORDER],
               edgecolor='black', linewidth=0.5, width=0.6)
        ax.plot(range(3), fracs, 'k--', linewidth=1.5, marker='D', markersize=5,
                markerfacecolor='black')
        for i, f in enumerate(fracs):
            ax.text(i, f+2, f"{f:.0f}%", ha='center', fontsize=10, fontweight='bold')
        ax.set_xticks(range(3))
        ax.set_xticklabels(GROUP_ORDER)
        ax.set_ylabel("% with support")
        ax.set_ylim(0, 115)
        ax.set_title("C  Support trend", fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if tests:
            ca_p = tests.get('CA_P', '?')
            ax.text(1, 108, f"CA P = {ca_p}", ha='center', fontsize=8, fontstyle='italic')

        # Panel D: Stats text
        ax = axes[3]
        ax.axis('off')
        if tests:
            lines = [
                f"Fisher (INV vs non-INV):",
                f"  OR = {tests.get('fisher_OR','?')} ({tests.get('fisher_CI_lo','?')}–{tests.get('fisher_CI_hi','?')})",
                f"  P = {tests.get('fisher_P','?')} {tests.get('fisher_sig','')}",
                "",
                f"Cochran–Armitage trend:",
                f"  Z = {tests.get('CA_Z','?')}, P = {tests.get('CA_P','?')} {tests.get('CA_sig','')}",
                "",
                f"Concordance: {tests.get('concordance','?')}",
                f"INV support: {tests.get('inv_support_frac','?')}",
            ]
            y = 0.92
            for line in lines:
                weight = 'bold' if line and not line.startswith(' ') else 'normal'
                color = '#DC143C' if '***' in line or '**' in line else '#333333'
                ax.text(0.05, y, line, fontsize=9, fontweight=weight, color=color,
                        transform=ax.transAxes, va='top', fontfamily='monospace')
                y -= 0.09
        ax.set_title("D  Statistical tests", fontweight='bold', loc='left')

        plt.tight_layout()
        for ext in ['pdf', 'png']:
            fig.savefig(os.path.join(args.plot_dir, f"{inv_id}_validation.{ext}"),
                       dpi=350, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  {inv_id}: plotted")

    # ── Genome-wide summary ──────────────────────────────────────────────
    master_file = os.path.join(args.stats_dir, "all_candidates_tests.tsv")
    if os.path.exists(master_file):
        results = []
        with open(master_file) as f:
            h = f.readline().strip().split('\t')
            for line in f:
                results.append(dict(zip(h, line.strip().split('\t'))))

        if results:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

            # Left: OR forest plot
            ids = [r['inv_id'] for r in results]
            ors = [float(r['fisher_OR']) for r in results]
            ci_lo = [float(r['fisher_CI_lo']) for r in results]
            ci_hi = [float(r['fisher_CI_hi']) for r in results]

            y_pos = range(len(ids))
            for i, (o, lo, hi) in enumerate(zip(ors, ci_lo, ci_hi)):
                hi_clip = min(hi, 500)
                color = '#DC143C' if float(results[i].get('fisher_P', 1)) < 0.05 else '#999999'
                ax1.plot([lo, hi_clip], [i, i], color=color, linewidth=2)
                ax1.plot(o, i, 'D', color=color, markersize=6)
            ax1.axvline(1, color='black', linestyle='--', linewidth=0.8)
            ax1.set_yticks(list(y_pos))
            ax1.set_yticklabels(ids, fontsize=7)
            ax1.set_xlabel("Odds Ratio (log scale)")
            ax1.set_xscale('log')
            ax1.set_title("Fisher OR (INV vs non-INV)", fontweight='bold')
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

            # Right: CA Z-scores
            zs = [float(r['CA_Z']) for r in results]
            colors = ['#DC143C' if float(r.get('CA_P', 1)) < 0.05 else '#999999' for r in results]
            ax2.barh(list(y_pos), zs, color=colors, edgecolor='black', linewidth=0.3)
            ax2.axvline(0, color='black', linewidth=0.8)
            ax2.set_yticks(list(y_pos))
            ax2.set_yticklabels(ids, fontsize=7)
            ax2.set_xlabel("Cochran–Armitage Z")
            ax2.set_title("Trend test (REF→HET→INV)", fontweight='bold')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

            plt.tight_layout()
            for ext in ['pdf', 'png']:
                fig.savefig(os.path.join(args.plot_dir, f"genome_wide_validation_summary.{ext}"),
                           dpi=350, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            print(f"\nGenome-wide summary: {args.plot_dir}/genome_wide_validation_summary.{{pdf,png}}")

    print("\nAll plots complete.")

if __name__ == "__main__":
    main()
