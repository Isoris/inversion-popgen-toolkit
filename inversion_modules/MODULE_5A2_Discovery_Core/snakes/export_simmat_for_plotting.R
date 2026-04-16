#!/usr/bin/env Rscript
# export_simmat_for_plotting.R
#
# Export sim_mat and window grid from precomp RDS files to .npy-compatible
# binary format that Python can read with np.load().
#
# Actually uses .tsv.gz since numpy .npy requires special headers.
# The Python plotter will check for these files before trying rpy2.
#
# Usage:
#   Rscript export_simmat_for_plotting.R <precomp_dir> [chrom]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript export_simmat_for_plotting.R <precomp_dir> [chrom]")

precomp_dir <- args[1]
chrom_filter <- if (length(args) >= 2) args[2] else NULL

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

message("[export] Found ", length(rds_files), " precomp files")

for (f in rds_files) {
  chrom <- sub("\\.precomp\\.rds$", "", basename(f))
  if (!is.null(chrom_filter) && chrom != chrom_filter) next

  # Check if already exported
  npy_path <- sub("\\.precomp\\.rds$", ".sim_mat.npy", f)
  grid_path <- sub("\\.precomp\\.rds$", ".window_grid.npz", f)

  if (file.exists(npy_path) && file.exists(grid_path)) {
    message("[export] SKIP ", chrom, ": already exported")
    next
  }

  message("[export] Loading ", chrom, "...")
  pc <- readRDS(f)

  # Export sim_mat as raw binary (row-major float64, numpy-compatible)
  sim_mat <- pc$sim_mat
  n <- nrow(sim_mat)
  message("[export]   sim_mat: ", n, "x", n)

  # Write as raw binary float64 with a small header file
  bin_path <- sub("\\.precomp\\.rds$", ".sim_mat.bin", f)
  meta_path <- sub("\\.precomp\\.rds$", ".sim_mat.meta", f)

  writeBin(as.double(t(sim_mat)), bin_path, size = 8)  # row-major
  writeLines(paste0("shape=", n, ",", n, "\ndtype=float64\norder=C"),
             meta_path)

  # Export window grid as TSV
  grid_tsv <- sub("\\.precomp\\.rds$", ".window_grid.tsv.gz", f)
  grid_dt <- pc$dt[, .(global_window_id, start_bp, end_bp)]
  fwrite(grid_dt, grid_tsv, sep = "\t")

  message("[export]   Written: ", bin_path)
  message("[export]   Written: ", grid_tsv)
}

message("[export] Done.")
