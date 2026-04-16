#!/usr/bin/env python3
###############################################################################
# S06_plot_species_assign.py
#
# Plots species assignment results: scatter, histogram, and suspicious-sample
# ranking from Mash distance tables.
#
# Usage:
#   python3 S06_plot_species_assign.py \
#     --per_cga results.per_CGA.tsv \
#     [--per_run results.tsv] \
#     [--outprefix mash_species_assign] \
#     [--outdir Figures/]
#
# Outputs (in outdir):
#   <prefix>.scatter_perCGA.{png,pdf}
#   <prefix>.hist_delta_perCGA.{png,pdf}
#   <prefix>.top_suspicious_perCGA.{png,pdf}
#   <prefix>.most_suspicious_perCGA.tsv     (Tables/)
#   [optional] <prefix>.scatter_perRUN.{png,pdf}
#   [optional] <prefix>.most_suspicious_perRUN.tsv
#   S06_plot_species_assign.arg
#   S06_plot_species_assign.results
###############################################################################
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def write_arg(args, outdir):
    path = os.path.join(outdir, "S06_plot_species_assign.arg")
    with open(path, "w") as f:
        f.write("key\tvalue\n")
        f.write(f"step\tS06_plot_species_assign\n")
        f.write(f"per_cga\t{args.per_cga}\n")
        f.write(f"per_run\t{args.per_run or 'NA'}\n")
        f.write(f"outprefix\t{args.outprefix}\n")
        f.write(f"outdir\t{outdir}\n")
        f.write(f"margin\t{args.margin}\n")
        f.write(f"label_top\t{args.label_top}\n")
    return path


def write_results(outdir, files_written, arg_path):
    path = os.path.join(outdir, "S06_plot_species_assign.results")
    with open(path, "w") as f:
        f.write("key\tvalue\tdescription\n")
        f.write(f"arg_file\t{arg_path}\tParameter record\n")
        for fp in files_written:
            f.write(f"output\t{fp}\t{os.path.basename(fp)}\n")


def main():
    ap = argparse.ArgumentParser(description="Plot Mash species assignment results")
    ap.add_argument("--per_cga", required=True, help="results.per_CGA.tsv")
    ap.add_argument("--per_run", default=None, help="results.tsv (optional)")
    ap.add_argument("--outprefix", default="mash_species_assign")
    ap.add_argument("--outdir", default=".", help="Output directory for figures")
    ap.add_argument("--tabledir", default=None, help="Output directory for TSV tables (default: outdir)")
    ap.add_argument("--margin", type=float, default=0.002)
    ap.add_argument("--label_top", type=int, default=15)
    args = ap.parse_args()

    outdir = args.outdir
    tabledir = args.tabledir or outdir
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tabledir, exist_ok=True)

    arg_path = write_arg(args, outdir)
    files_written = []

    df = pd.read_csv(args.per_cga, sep="\t")
    df["median_dist_gar"] = pd.to_numeric(df["median_dist_gar"], errors="coerce")
    df["median_dist_mac"] = pd.to_numeric(df["median_dist_mac"], errors="coerce")
    df["delta_mac_minus_gar"] = df["median_dist_mac"] - df["median_dist_gar"]
    df["absdiff"] = df["delta_mac_minus_gar"].abs()

    # ---- 1) Scatter: dist_gar vs dist_mac ----
    x = df["median_dist_gar"].to_numpy()
    y = df["median_dist_mac"].to_numpy()
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    df_sc = df.loc[finite].copy()

    fig, ax = plt.subplots()
    ax.scatter(x, y, s=18, alpha=0.7)
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.4, linewidth=0.8)
    ax.set_xlabel("Mash dist to Gar")
    ax.set_ylabel("Mash dist to Mac")
    ax.set_title("Species assignment: Gar vs Mac (per CGA)")

    df_sc = df_sc.sort_values("delta_mac_minus_gar", ascending=True)
    for _, r in df_sc.head(args.label_top).iterrows():
        ax.annotate(r["CGA_sample"],
                     (r["median_dist_gar"], r["median_dist_mac"]),
                     fontsize=7, xytext=(3, 3), textcoords="offset points")

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        p = os.path.join(outdir, f"{args.outprefix}.scatter_perCGA.{ext}")
        fig.savefig(p, dpi=200)
        files_written.append(p)
    plt.close(fig)

    # ---- 2) Histogram of delta ----
    d = df["delta_mac_minus_gar"].dropna().to_numpy()
    fig, ax = plt.subplots()
    ax.hist(d, bins=40, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("delta = dist_mac - dist_gar")
    ax.set_ylabel("Count")
    ax.set_title("Reference separation (per CGA)")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        p = os.path.join(outdir, f"{args.outprefix}.hist_delta_perCGA.{ext}")
        fig.savefig(p, dpi=200)
        files_written.append(p)
    plt.close(fig)

    # ---- 3) Ranked suspicious table + plot ----
    sus = df.sort_values("delta_mac_minus_gar", ascending=True).copy()
    sus_cols = ["CGA_sample", "n_runs", "median_dist_gar", "median_dist_mac",
                "delta_mac_minus_gar", "majority_call", "calls_breakdown"]
    sus_path = os.path.join(tabledir, f"{args.outprefix}.most_suspicious_perCGA.tsv")
    sus[sus_cols].to_csv(sus_path, sep="\t", index=False)
    files_written.append(sus_path)

    topn = min(30, len(sus))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.arange(topn), sus["delta_mac_minus_gar"].head(topn).to_numpy(),
            marker="o", linestyle="-", markersize=4)
    ax.set_xticks(np.arange(topn))
    ax.set_xticklabels(sus["CGA_sample"].head(topn).to_list(), rotation=90, fontsize=7)
    ax.set_ylabel("delta (smaller = more suspicious)")
    ax.set_title(f"Top {topn} most suspicious CGAs")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        p = os.path.join(outdir, f"{args.outprefix}.top_suspicious_perCGA.{ext}")
        fig.savefig(p, dpi=200)
        files_written.append(p)
    plt.close(fig)

    # ---- Optional: per-run scatter ----
    if args.per_run:
        dr = pd.read_csv(args.per_run, sep="\t")
        dr["dist_gar"] = pd.to_numeric(dr.get("dist_gar", pd.Series(dtype=float)), errors="coerce")
        dr["dist_mac"] = pd.to_numeric(dr.get("dist_mac", pd.Series(dtype=float)), errors="coerce")

        xr = dr["dist_gar"].to_numpy()
        yr = dr["dist_mac"].to_numpy()
        fr = np.isfinite(xr) & np.isfinite(yr)

        if fr.sum() > 0:
            fig, ax = plt.subplots()
            ax.scatter(xr[fr], yr[fr], s=10, alpha=0.5)
            lo = min(xr[fr].min(), yr[fr].min())
            hi = max(xr[fr].max(), yr[fr].max())
            ax.plot([lo, hi], [lo, hi], "k--", alpha=0.4, linewidth=0.8)
            ax.set_xlabel("dist_gar (per run)")
            ax.set_ylabel("dist_mac (per run)")
            ax.set_title("Mash distance: Gar vs Mac (per run-lane)")
            fig.tight_layout()
            for ext in ["png", "pdf"]:
                p = os.path.join(outdir, f"{args.outprefix}.scatter_perRUN.{ext}")
                fig.savefig(p, dpi=200)
                files_written.append(p)
            plt.close(fig)

        dr["delta_mac_minus_gar"] = dr["dist_mac"] - dr["dist_gar"]
        sus_run_path = os.path.join(tabledir, f"{args.outprefix}.most_suspicious_perRUN.tsv")
        dr.sort_values("delta_mac_minus_gar", ascending=True).to_csv(sus_run_path, sep="\t", index=False)
        files_written.append(sus_run_path)

    write_results(outdir, files_written, arg_path)

    print("Wrote:")
    for f in files_written:
        print(f"  {f}")


if __name__ == "__main__":
    main()
