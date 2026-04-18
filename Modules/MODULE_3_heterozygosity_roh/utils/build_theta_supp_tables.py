#!/usr/bin/env python3
"""
build_theta_supp_tables.py  (v2)

Build supplementary tables for windowed nucleotide diversity (theta_pi) from
ANGSD thetaStat .pestPG output files.

FIX vs v1: correct pestPG column layout. ANGSD thetaStat do_stat writes 14
columns, not 11. The first column is a composite parens string like
   (indexStart,indexStop)(firstPos,lastPos)(WinStart,WinStop)
followed by Chr, WinCenter, tW, tP, tF, tH, tL, Tajima, fuf, fud, fayh, zeng, nSites.

Per-site theta values = tP/nSites (for θ_π) or tW/nSites (for Watterson).

USAGE:
  python build_theta_supp_tables.py \
      --theta-dir /path/to/pestPG/ \
      --window 500000 --step 500000 \
      --out-dir /path/to/out/ \
      [--samples sample_list.txt] \
      [--min-sites 1000]
"""

import argparse
import os
import sys
import re
import glob
import pandas as pd
import numpy as np

# Correct pestPG column layout (14 columns).
# Column 1 is a composite parens string we keep but usually don't need.
PESTPG_COLS = [
    "WinInfo", "Chr", "WinCenter",
    "tW", "tP", "tF", "tH", "tL",
    "Tajima", "fuf", "fud", "fayh", "zeng", "nSites",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--theta-dir", required=True)
    p.add_argument("--window", type=int, required=True)
    p.add_argument("--step", type=int, required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--samples", default=None,
                   help="Optional file with sample IDs to restrict to (one per line)")
    p.add_argument("--min-sites", type=int, default=1000)
    return p.parse_args()


def discover_files(theta_dir, win, step, samples=None):
    pattern = os.path.join(theta_dir, f"*.win{win}.step{step}.pestPG")
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"No files matching {pattern}")
    rx = re.compile(rf"(.+)\.win{win}\.step{step}\.pestPG$")
    sample_files = {}
    for f in files:
        m = rx.search(os.path.basename(f))
        if m:
            sample_files[m.group(1)] = f
    if samples is not None:
        with open(samples) as fh:
            keep = {ln.strip() for ln in fh if ln.strip()}
        sample_files = {s: f for s, f in sample_files.items() if s in keep}
        if not sample_files:
            sys.exit(f"No files match the sample list {samples}")
    print(f"Found {len(sample_files)} samples", file=sys.stderr)
    return sample_files


def read_pestpg(path, sample, min_sites):
    """Read one pestPG file using the correct 14-column layout."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n")
    if not header.startswith("#"):
        sys.exit(f"{path}: expected header line starting with #. Got: {header[:80]!r}")
    n_header_fields = len(header.split("\t"))
    if n_header_fields != len(PESTPG_COLS):
        sys.stderr.write(
            f"  WARNING: {os.path.basename(path)} has {n_header_fields} header fields "
            f"but expected {len(PESTPG_COLS)}. Reading anyway — check output carefully.\n"
        )

    df = pd.read_csv(path, sep="\t", comment=None,
                     skiprows=1, names=PESTPG_COLS)
    for col in ["WinCenter", "tW", "tP", "Tajima", "nSites"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    n_before = len(df)
    df = df.dropna(subset=["nSites", "tP"])
    df = df[df["nSites"] >= min_sites].copy()
    if len(df) == 0:
        sys.stderr.write(
            f"  WARNING: {os.path.basename(path)}: 0 windows passed min-sites filter "
            f"(had {n_before} total)\n"
        )
        return None
    df["theta_pi"] = df["tP"] / df["nSites"]
    df["theta_W"]  = df["tW"] / df["nSites"]
    df["sample"]   = sample
    return df[["sample", "Chr", "WinCenter", "nSites", "theta_pi", "theta_W", "Tajima"]]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sample_files = discover_files(args.theta_dir, args.window, args.step,
                                  samples=args.samples)

    dfs = []
    for sample, path in sample_files.items():
        try:
            d = read_pestpg(path, sample, args.min_sites)
            if d is not None:
                dfs.append(d)
        except Exception as e:
            print(f"  WARNING: failed to read {path}: {e}", file=sys.stderr)
    if not dfs:
        sys.exit("No usable pestPG files read.")
    long = pd.concat(dfs, ignore_index=True)

    # SANITY CHECK
    mean_tp = long["theta_pi"].mean()
    if not np.isfinite(mean_tp) or mean_tp == 0:
        sys.exit(
            f"ERROR: mean θ_π across all windows is {mean_tp}. "
            f"This usually means the pestPG column layout doesn't match. "
            f"Inspect one file with:\n"
            f"   head -2 {list(sample_files.values())[0]}\n"
            f"and compare to the expected 14-column layout "
            f"(WinInfo, Chr, WinCenter, tW, tP, tF, tH, tL, Tajima, "
            f"fuf, fud, fayh, zeng, nSites)."
        )
    print(f"Total (sample, window) rows: {len(long):,}", file=sys.stderr)
    print(f"  sanity:  mean θ_π = {mean_tp:.3e}  median = {long['theta_pi'].median():.3e}",
          file=sys.stderr)

    # Per-window table
    long_out = long.rename(columns={
        "Chr": "Chromosome", "WinCenter": "Window centre (bp)",
        "nSites": "n sites", "theta_pi": "θ_π (per site)",
        "theta_W": "θ_W (per site)", "Tajima": "Tajima D"
    })
    path1 = os.path.join(args.out_dir, "ST_theta_pi_per_window.tsv")
    long_out.to_csv(path1, sep="\t", index=False)
    print(f"  wrote {path1}", file=sys.stderr)

    # Per-chromosome summary: collapse samples→mean per window first,
    # then summarise across windows. Correct approach for a population
    # diversity statistic.
    per_window = (long.groupby(["Chr", "WinCenter"])
                  .agg(theta_pi=("theta_pi", "mean"),
                       theta_W=("theta_W", "mean"),
                       n_samples=("theta_pi", "size"))
                  .reset_index())

    per_chr = (per_window.groupby("Chr")["theta_pi"]
               .agg(["count", "mean", "median",
                     lambda x: x.quantile(0.25), lambda x: x.quantile(0.75),
                     "min", "max"])
               .reset_index())
    per_chr.columns = ["Chromosome", "n windows", "Mean θ_π", "Median θ_π",
                       "Q25 θ_π", "Q75 θ_π", "Min θ_π", "Max θ_π"]
    per_chr["IQR θ_π"] = per_chr["Q75 θ_π"] - per_chr["Q25 θ_π"]
    per_chr = per_chr[["Chromosome", "n windows", "Mean θ_π", "Median θ_π",
                       "IQR θ_π", "Q25 θ_π", "Q75 θ_π", "Min θ_π", "Max θ_π"]]
    per_chr["_lg"] = per_chr["Chromosome"].str.extract(r"LG(\d+)").astype(float)
    per_chr = per_chr.sort_values("_lg").drop(columns="_lg")
    path2 = os.path.join(args.out_dir, "ST_theta_pi_per_chr_summary.tsv")
    per_chr.to_csv(path2, sep="\t", index=False)
    print(f"  wrote {path2}", file=sys.stderr)

    # Per-sample summary
    per_sample = (long.groupby("sample")["theta_pi"]
                  .agg(["count", "mean", "median",
                        lambda x: x.quantile(0.25), lambda x: x.quantile(0.75),
                        "min", "max"])
                  .reset_index())
    per_sample.columns = ["Sample", "n windows", "Mean θ_π", "Median θ_π",
                          "Q25 θ_π", "Q75 θ_π", "Min θ_π", "Max θ_π"]
    per_sample["IQR θ_π"] = per_sample["Q75 θ_π"] - per_sample["Q25 θ_π"]
    per_sample = per_sample[["Sample", "n windows", "Mean θ_π", "Median θ_π",
                             "IQR θ_π", "Q25 θ_π", "Q75 θ_π", "Min θ_π", "Max θ_π"]]
    path3 = os.path.join(args.out_dir, "ST_theta_pi_per_sample_summary.tsv")
    per_sample.to_csv(path3, sep="\t", index=False)
    print(f"  wrote {path3}", file=sys.stderr)

    # Outlier windows
    p99 = per_window["theta_pi"].quantile(0.99)
    p995 = per_window["theta_pi"].quantile(0.995)
    genome_mean = per_window["theta_pi"].mean()
    outliers = per_window[per_window["theta_pi"] > p99].copy()
    outliers["Percentile"] = outliers["theta_pi"].apply(
        lambda x: ">99.5" if x > p995 else "99-99.5"
    )
    outliers["Fold over genome mean"] = outliers["theta_pi"] / genome_mean
    outliers = outliers.rename(columns={
        "Chr": "Chromosome", "WinCenter": "Window centre (bp)",
        "theta_pi": "Mean θ_π (across samples)",
        "n_samples": "n samples",
    })[["Chromosome", "Window centre (bp)", "Mean θ_π (across samples)",
        "Percentile", "Fold over genome mean", "n samples"]]
    outliers = outliers.sort_values("Mean θ_π (across samples)", ascending=False)
    path4 = os.path.join(args.out_dir, "ST_theta_pi_outlier_windows.tsv")
    outliers.to_csv(path4, sep="\t", index=False)
    print(f"  wrote {path4}  ({len(outliers)} windows above 99th percentile)",
          file=sys.stderr)
    print(f"  genome mean θ_π:    {genome_mean:.3e}", file=sys.stderr)
    print(f"  99th percentile:    {p99:.3e}", file=sys.stderr)
    print(f"  99.5th percentile:  {p995:.3e}", file=sys.stderr)


if __name__ == "__main__":
    main()
