#!/usr/bin/env python3
"""
plot_triangle_insulation.py

Hi-C-style rotated triangle heatmaps of the sim_mat (local PCA similarity)
with triangle interval boundaries overlaid and colored by type/nesting.

Produces:
  1. Whole-chromosome triangle-down heatmap with all intervals outlined
  2. Per-interval zoomed views (optional)
  3. Nested sub-regime coloring within intervals

Adapted from MODULE_5C STEP_LD_02_plot_heatmaps.py (same rotated-coord approach).

Usage:
  python plot_triangle_insulation.py \
    --precomp-dir precomp/ \
    --triangle-dir triangles_v2/ \
    --outdir plots_triangles/ \
    [--chrom C_gar_LG01] \
    [--zoom-top 20] \
    [--dpi 400]
"""

import argparse, os, sys, gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
import csv

# ─────────────────────────────────────────────────────────────────────
# PALETTES
# ─────────────────────────────────────────────────────────────────────

# Sim_mat heatmap: dark-to-hot (like your LD plots)
CMAP_SIM = LinearSegmentedColormap.from_list("sim_apex", [
    "#080808", "#1a0a2e", "#3d1261", "#6a1b8a",
    "#9c2e7a", "#cc4466", "#e8734a", "#f5a623"
], N=256)
CMAP_SIM.set_bad(color="white")

# Interval type → outline color
TYPE_COLORS = {
    "strong_triangle":     "#FF1744",   # vivid red
    "moderate_triangle":   "#FF9100",   # amber/orange
    "sharp_but_not_square":"#FFEA00",   # electric yellow
    "patchy_signal":       "#00E676",   # green
    "diffuse_zone":        "#00B0FF",   # cyan
    "weak_zone":           "#B388FF",   # light purple
    "transition":          "#78909C",   # blue-grey
}

# Sub-regime level → fill alpha overlay
SUBREG_COLORS = {
    "hot":  "#FF1744",
    "warm": "#FF9100",
    "cool": "#4FC3F7",
    "cold": "#90A4AE",
}

# ─────────────────────────────────────────────────────────────────────
# I/O HELPERS
# ─────────────────────────────────────────────────────────────────────

def read_tsv_gz(path):
    """Read a gzipped TSV into a list of dicts."""
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows

def load_precomp_simmat(rds_path):
    """
    Load sim_mat from precomp RDS via multiple fallback paths:
    1. .npy file (numpy native)
    2. .bin + .meta file (raw float64 binary, exported by R helper)
    3. rpy2 direct RDS loading
    """
    npy_path = rds_path.replace(".precomp.rds", ".sim_mat.npy")
    npz_path = rds_path.replace(".precomp.rds", ".sim_mat.npz")
    bin_path = rds_path.replace(".precomp.rds", ".sim_mat.bin")
    meta_path = rds_path.replace(".precomp.rds", ".sim_mat.meta")

    if os.path.exists(npy_path):
        return np.load(npy_path)
    if os.path.exists(npz_path):
        d = np.load(npz_path)
        return d[list(d.keys())[0]]

    # Binary format from R helper
    if os.path.exists(bin_path) and os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = dict(line.strip().split("=", 1) for line in f if "=" in line)
        shape = tuple(int(x) for x in meta["shape"].split(","))
        mat = np.fromfile(bin_path, dtype=np.float64).reshape(shape)
        print(f"[plot] Loaded sim_mat from binary: {shape}", file=sys.stderr)
        # Cache as npy for next time
        np.save(npy_path, mat)
        return mat

    # Try rpy2
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
        pc = ro.r(f'readRDS("{rds_path}")')
        sim_mat = np.array(pc.rx2("sim_mat"))
        np.save(npy_path, sim_mat)
        print(f"[plot] Loaded sim_mat via rpy2: {sim_mat.shape}", file=sys.stderr)
        return sim_mat
    except Exception as e:
        print(f"[plot] WARNING: Cannot load {rds_path}: {e}", file=sys.stderr)
        print(f"[plot] Export sim_mat first with:", file=sys.stderr)
        print(f"  Rscript export_simmat_for_plotting.R {os.path.dirname(rds_path)}",
              file=sys.stderr)
        return None

def load_window_grid(rds_path):
    """Load window start/end from precomp RDS → dt$start_bp, dt$end_bp."""
    npy_path = rds_path.replace(".precomp.rds", ".window_grid.npz")
    tsv_path = rds_path.replace(".precomp.rds", ".window_grid.tsv.gz")

    if os.path.exists(npy_path):
        d = np.load(npy_path)
        return d["start_bp"], d["end_bp"]

    # TSV format from R helper
    if os.path.exists(tsv_path):
        rows = read_tsv_gz(tsv_path)
        start_bp = np.array([int(r["start_bp"]) for r in rows])
        end_bp = np.array([int(r["end_bp"]) for r in rows])
        np.savez(npy_path, start_bp=start_bp, end_bp=end_bp)
        print(f"[plot] Loaded window grid from TSV: {len(start_bp)} windows",
              file=sys.stderr)
        return start_bp, end_bp

    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
        pc = ro.r(f'readRDS("{rds_path}")')
        dt = pc.rx2("dt")
        start_bp = np.array(dt.rx2("start_bp"))
        end_bp = np.array(dt.rx2("end_bp"))
        np.savez(npy_path, start_bp=start_bp, end_bp=end_bp)
        return start_bp, end_bp
    except Exception as e:
        print(f"[plot] WARNING: Cannot load window grid: {e}", file=sys.stderr)
        return None, None


# ─────────────────────────────────────────────────────────────────────
# BIN sim_mat (same as R code: average blocks of bin_size windows)
# ─────────────────────────────────────────────────────────────────────

def bin_matrix(mat, bin_size):
    n = mat.shape[0]
    nb = n // bin_size
    bmat = np.zeros((nb, nb))
    for bi in range(nb):
        ri = slice(bi * bin_size, (bi + 1) * bin_size)
        for bj in range(nb):
            rj = slice(bj * bin_size, (bj + 1) * bin_size)
            block = mat[ri, rj]
            bmat[bi, bj] = np.nanmean(block)
    return bmat


# ─────────────────────────────────────────────────────────────────────
# ROTATED TRIANGLE PLOTTER (adapted from MODULE_5C)
# ─────────────────────────────────────────────────────────────────────

def plot_triangle_heatmap(
    sim_mat, chrom, start_bp, end_bp,
    intervals, subregimes, bridges,
    outpath, title=None,
    bin_size=2,
    region_start_mb=None, region_end_mb=None,
    dpi=400, width=14, height=None,
    power=0.5, vmax_quantile=0.995,
    show_subreg=True, show_bridges=True,
    label_intervals=True
):
    """
    Draw a rotated triangle-down heatmap of sim_mat with interval outlines.

    sim_mat:     N×N similarity matrix (raw window-level)
    intervals:   list of dicts from triangle_intervals.tsv.gz
    subregimes:  list of dicts from triangle_subregimes.tsv.gz
    bridges:     list of dicts from triangle_bridges.tsv.gz
    """
    n_raw = sim_mat.shape[0]

    # Bin the matrix
    bmat = bin_matrix(sim_mat, bin_size) if bin_size > 1 else sim_mat.copy()
    n = bmat.shape[0]

    # Window positions in Mb (bin centers)
    if start_bp is not None and len(start_bp) >= n_raw:
        bin_pos_mb = np.array([
            (start_bp[i * bin_size] + end_bp[min((i + 1) * bin_size - 1, n_raw - 1)]) / 2e6
            for i in range(n)
        ])
        dx_mb = np.median(np.diff(bin_pos_mb)) if n > 1 else 0.01
    else:
        dx_mb = 0.01
        bin_pos_mb = np.arange(n) * dx_mb

    x0 = bin_pos_mb[0] - dx_mb / 2
    x_max = bin_pos_mb[-1] + dx_mb / 2

    # Region crop
    if region_start_mb is not None and region_end_mb is not None:
        mask = (bin_pos_mb >= region_start_mb) & (bin_pos_mb <= region_end_mb)
        idx = np.where(mask)[0]
        if len(idx) < 5:
            print(f"[plot] Region too small after crop, skipping", file=sys.stderr)
            return
        bmat = bmat[np.ix_(idx, idx)]
        bin_pos_mb = bin_pos_mb[idx]
        n = len(idx)
        x0 = bin_pos_mb[0] - dx_mb / 2
        x_max = bin_pos_mb[-1] + dx_mb / 2

    # Upper triangle only
    mat = bmat.copy()
    mat[np.tril_indices(n, k=-1)] = np.nan

    # Transform — use power < 0.5 for more contrast in the mid-range
    display = np.ma.masked_invalid(mat)
    display = np.power(np.clip(display, 0, None), power)
    vals_valid = display.compressed()
    if len(vals_valid) > 0:
        vmax = float(np.quantile(vals_valid, vmax_quantile))
    else:
        vmax = 1.0

    if height is None:
        height = width * 0.45

    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Build rotated coordinate mesh
    xs = np.zeros(n + 1)
    xs[0] = x0
    for i in range(n):
        xs[i + 1] = bin_pos_mb[i] + dx_mb / 2

    X = np.zeros((n + 1, n + 1))
    Y = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(n + 1):
            X[i, j] = (xs[i] + xs[j]) / 2
            Y[i, j] = (xs[j] - xs[i]) / 2

    pcm = ax.pcolormesh(X, Y, display, cmap=CMAP_SIM, vmin=0, vmax=vmax,
                        shading='flat', edgecolors='none', linewidth=0,
                        rasterized=True)

    y_max_val = (x_max - x0) / 2
    ax.set_xlim(x0, x_max)
    ax.set_ylim(y_max_val * 1.02, -y_max_val * 0.02)  # triangle_down
    ax.set_aspect("equal")

    # ── Overlay interval boundaries ──
    for iv in intervals:
        iv_start = float(iv["start_mb"])
        iv_end = float(iv["end_mb"])
        iv_type = iv["interval_type"]
        sq = float(iv.get("squareness", 0))

        # Skip if outside view
        if iv_end < bin_pos_mb[0] or iv_start > bin_pos_mb[-1]:
            continue

        color = TYPE_COLORS.get(iv_type, "#AAAAAA")
        # Thicker line for higher squareness
        lw = 1.0 + 2.0 * min(sq, 1.0)  # 1.0 to 3.0
        alpha = 0.4 + 0.5 * min(sq, 1.0)  # 0.4 to 0.9

        # In rotated coordinates, the interval [s, e] forms a diamond/triangle:
        # Bottom-left corner: (s, 0)
        # Bottom-right corner: (e, 0)
        # Top apex: ((s+e)/2, (e-s)/2)
        s, e = iv_start, iv_end
        mid = (s + e) / 2
        h = (e - s) / 2

        triangle = plt.Polygon(
            [(s, 0), (mid, h), (e, 0)],
            closed=True, fill=False,
            edgecolor=color, linewidth=lw, alpha=alpha,
            linestyle="-", zorder=10
        )
        ax.add_patch(triangle)

        # Label with interval ID
        if label_intervals and h > y_max_val * 0.02:
            iid = iv.get("interval_id", "")
            label_y = h * 0.85
            fontsize = max(4, min(6, h / y_max_val * 60))
            ax.text(mid, label_y, f"I{iid}",
                    ha="center", va="center", fontsize=fontsize,
                    color=color, alpha=alpha * 0.9, fontweight="bold",
                    zorder=15)

    # ── Sub-regime fills (subtle colored bands at the base) ──
    if show_subreg and subregimes:
        bar_h = y_max_val * 0.015  # thin bar at y=0
        for sr in subregimes:
            sr_start = float(sr["start_mb"])
            sr_end = float(sr["end_mb"])
            level = sr.get("sub_level", "cool")
            if sr_end < bin_pos_mb[0] or sr_start > bin_pos_mb[-1]:
                continue
            fc = SUBREG_COLORS.get(level, "#90A4AE")
            rect = plt.Rectangle(
                (sr_start, -bar_h * 1.5), sr_end - sr_start, bar_h,
                facecolor=fc, edgecolor="none", alpha=0.7,
                clip_on=False, zorder=20
            )
            ax.add_patch(rect)

    # ── Bridges (connecting lines between intervals) ──
    if show_bridges and bridges:
        for br in bridges:
            # Bridges connect two intervals — draw as thin lines
            # We'd need interval coordinates; skip for now if not in data
            pass

    # ── Colorbar ──
    sm = plt.cm.ScalarMappable(cmap=CMAP_SIM, norm=Normalize(0, vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, shrink=0.6)
    tlab = f"sim^{power}" if power != 1 else "similarity"
    cbar.set_label(tlab, fontsize=9)

    ax.set_xlabel(f"{chrom} position (Mb)", fontsize=10)
    ax.set_ylabel("Pairwise distance (Mb)", fontsize=10)
    ax.tick_params(labelsize=9)

    # ── Legend for interval types ──
    legend_items = []
    types_present = set(iv["interval_type"] for iv in intervals
                        if float(iv["end_mb"]) >= bin_pos_mb[0]
                        and float(iv["start_mb"]) <= bin_pos_mb[-1])
    for t in ["strong_triangle", "moderate_triangle", "sharp_but_not_square",
              "patchy_signal", "diffuse_zone", "weak_zone"]:
        if t in types_present:
            legend_items.append(plt.Line2D([0], [0], color=TYPE_COLORS[t],
                                            lw=2, label=t.replace("_", " ")))
    if legend_items:
        ax.legend(handles=legend_items, loc="lower right", fontsize=6,
                  framealpha=0.85, fancybox=True, edgecolor="grey70")

    # ── Title ──
    if title is None:
        n_iv = len([iv for iv in intervals
                     if float(iv["end_mb"]) >= bin_pos_mb[0]
                     and float(iv["start_mb"]) <= bin_pos_mb[-1]])
        title = f"Triangle Insulation — {chrom}"
    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)

    # Subtitle
    n_strong = len([iv for iv in intervals
                     if iv["interval_type"] == "strong_triangle"
                     and float(iv["end_mb"]) >= bin_pos_mb[0]
                     and float(iv["start_mb"]) <= bin_pos_mb[-1]])
    n_mod = len([iv for iv in intervals
                  if iv["interval_type"] == "moderate_triangle"
                  and float(iv["end_mb"]) >= bin_pos_mb[0]
                  and float(iv["start_mb"]) <= bin_pos_mb[-1]])
    n_total = len([iv for iv in intervals
                    if float(iv["end_mb"]) >= bin_pos_mb[0]
                    and float(iv["start_mb"]) <= bin_pos_mb[-1]])
    subtitle = (f"{n_total} intervals ({n_strong} strong, {n_mod} moderate) | "
                f"bin={bin_size} windows | power={power}")
    ax.text(0.5, 1.005, subtitle, transform=ax.transAxes,
            fontsize=7, ha="center", color="grey")

    fig.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[plot] Saved: {outpath}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser(description="Triangle insulation heatmap plotter")
    ap.add_argument("--precomp-dir", required=True,
                    help="Directory with .precomp.rds files")
    ap.add_argument("--triangle-dir", required=True,
                    help="Directory with triangle_intervals.tsv.gz etc.")
    ap.add_argument("--outdir", default="plots_triangles",
                    help="Output directory for plots")
    ap.add_argument("--chrom", default=None,
                    help="Single chromosome to plot (default: all)")
    ap.add_argument("--bin-size", type=int, default=2,
                    help="Bin size for sim_mat (default: 2)")
    ap.add_argument("--zoom-top", type=int, default=0,
                    help="Generate zoomed views for top N intervals by squareness")
    ap.add_argument("--zoom-pad-mb", type=float, default=2.0,
                    help="Padding around interval for zoom view (Mb)")
    ap.add_argument("--power", type=float, default=0.5,
                    help="Power transform for similarity values")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--width", type=float, default=14)
    ap.add_argument("--no-subreg", action="store_true")
    ap.add_argument("--no-labels", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load triangle results
    iv_path = os.path.join(args.triangle_dir, "triangle_intervals.tsv.gz")
    sr_path = os.path.join(args.triangle_dir, "triangle_subregimes.tsv.gz")
    br_path = os.path.join(args.triangle_dir, "triangle_bridges.tsv.gz")

    if not os.path.exists(iv_path):
        print(f"[ERROR] Not found: {iv_path}", file=sys.stderr)
        sys.exit(1)

    intervals = read_tsv_gz(iv_path)
    subregimes = read_tsv_gz(sr_path) if os.path.exists(sr_path) else []
    bridges = read_tsv_gz(br_path) if os.path.exists(br_path) else []

    print(f"[plot] Loaded {len(intervals)} intervals, "
          f"{len(subregimes)} sub-regimes, {len(bridges)} bridges",
          file=sys.stderr)

    # Get chromosomes
    chroms = sorted(set(iv["chrom"] for iv in intervals))
    if args.chrom:
        chroms = [c for c in chroms if c == args.chrom]
    if not chroms:
        print("[ERROR] No matching chromosomes", file=sys.stderr)
        sys.exit(1)

    # Find precomp files
    precomp_files = {}
    for f in sorted(os.listdir(args.precomp_dir)):
        if f.endswith(".precomp.rds"):
            chrom_name = f.replace(".precomp.rds", "")
            precomp_files[chrom_name] = os.path.join(args.precomp_dir, f)

    for chrom in chroms:
        if chrom not in precomp_files:
            print(f"[plot] SKIP {chrom}: no precomp RDS", file=sys.stderr)
            continue

        print(f"\n[plot] === {chrom} ===", file=sys.stderr)

        # Load sim_mat
        sim_mat = load_precomp_simmat(precomp_files[chrom])
        if sim_mat is None:
            continue

        start_bp, end_bp = load_window_grid(precomp_files[chrom])

        # Filter to this chromosome
        chr_iv = [iv for iv in intervals if iv["chrom"] == chrom]
        chr_sr = [sr for sr in subregimes if sr["chrom"] == chrom]
        chr_br = [br for br in bridges if br.get("chrom", "") == chrom]

        print(f"[plot] sim_mat: {sim_mat.shape[0]}×{sim_mat.shape[0]}, "
              f"{len(chr_iv)} intervals", file=sys.stderr)

        # 1. Whole-chromosome plot
        outpath = os.path.join(args.outdir,
                               f"{chrom}_triangle_insulation_wholechr.png")
        plot_triangle_heatmap(
            sim_mat, chrom, start_bp, end_bp,
            chr_iv, chr_sr, chr_br,
            outpath=outpath,
            bin_size=args.bin_size,
            dpi=args.dpi, width=args.width,
            power=args.power,
            show_subreg=not args.no_subreg,
            label_intervals=not args.no_labels,
        )

        # 2. Zoomed views for top intervals by squareness
        if args.zoom_top > 0:
            sorted_iv = sorted(chr_iv,
                               key=lambda x: float(x.get("squareness", 0)),
                               reverse=True)
            zoom_dir = os.path.join(args.outdir, "zooms")
            os.makedirs(zoom_dir, exist_ok=True)

            for rank, iv in enumerate(sorted_iv[:args.zoom_top], 1):
                iv_start = float(iv["start_mb"])
                iv_end = float(iv["end_mb"])
                iv_id = iv.get("interval_id", rank)
                iv_type = iv["interval_type"]
                sq = float(iv.get("squareness", 0))

                pad = args.zoom_pad_mb
                r_start = max(0, iv_start - pad)
                r_end = iv_end + pad

                zoom_path = os.path.join(
                    zoom_dir,
                    f"{chrom}_zoom_I{iv_id}_{iv_type}_sq{sq:.2f}.png"
                )
                plot_triangle_heatmap(
                    sim_mat, chrom, start_bp, end_bp,
                    chr_iv, chr_sr, chr_br,
                    outpath=zoom_path,
                    bin_size=max(1, args.bin_size // 2),  # finer binning for zoom
                    region_start_mb=r_start, region_end_mb=r_end,
                    dpi=args.dpi, width=12,
                    power=args.power,
                    show_subreg=not args.no_subreg,
                    label_intervals=not args.no_labels,
                    title=f"I{iv_id} — {iv_type} (sq={sq:.2f}) — {chrom}",
                )

    print(f"\n[plot] All done. Outputs in: {args.outdir}/", file=sys.stderr)


if __name__ == "__main__":
    main()
