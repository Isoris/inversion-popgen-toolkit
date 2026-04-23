#!/usr/bin/env Rscript
# =============================================================================
# q03_coverage_collapse.R
#
# Collapse per-sample mosdepth regions.bed.gz files into a single per-bin
# coverage track. Handles variable region sizes by re-binning to uniform
# BIN_MB windows via coverage-weighted mean.
#
# Input: a directory containing per-sample *.regions.bed.gz files, each with
#        4 columns (chrom, start, end, mean_cov_in_region)
# Output: chrom bin_start_bp bin_end_bp bin_mid_mb
#         mean_cov median_cov cv_across_samples n_samples_low_cov
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
FILES_DIR <- get_arg("--files_glob")
CHROM     <- get_arg("--chrom")
BIN_MB    <- as.numeric(get_arg("--bin_mb", "0.05"))
OUT       <- get_arg("--out")
LOW_FRAC  <- as.numeric(get_arg("--low_frac", "0.5"))  # samples < LOW_FRAC * chrom median
stopifnot(!is.null(FILES_DIR), !is.null(CHROM), !is.null(OUT))

BIN_BP <- as.integer(BIN_MB * 1e6)

files <- list.files(FILES_DIR, pattern = "\\.regions\\.bed\\.gz$",
                    recursive = TRUE, full.names = TRUE)
message("[q03] Found ", length(files), " mosdepth files in ", FILES_DIR)
if (length(files) == 0) stop("No *.regions.bed.gz files found")

# Determine genome extent for this chrom by peeking the first file
peek <- fread(cmd = sprintf("zcat %s | awk '$1==\"%s\"'",
                            shQuote(files[1]), CHROM),
              col.names = c("chr", "start", "end", "cov"))
if (nrow(peek) == 0) {
  stop("Chrom ", CHROM, " not found in first mosdepth file ", files[1])
}
max_bp <- max(peek$end)
n_bins <- as.integer(ceiling(max_bp / BIN_BP))
message("[q03] ", CHROM, ": ", n_bins, " bins to ", round(max_bp/1e6, 2), " Mb")

# Per-bin per-sample mean coverage matrix
cov_mat <- matrix(NA_real_, nrow = n_bins, ncol = length(files))

for (i in seq_along(files)) {
  b <- tryCatch(
    fread(cmd = sprintf("zcat %s | awk '$1==\"%s\"'",
                        shQuote(files[i]), CHROM),
          col.names = c("chr", "start", "end", "cov")),
    error = function(e) NULL)
  if (is.null(b) || nrow(b) == 0) next

  b[, mid := (start + end) / 2]
  b[, bin := as.integer(floor((mid - 1) / BIN_BP)) + 1L]
  b[, w   := (end - start)]
  # Weighted mean coverage per bin
  agg <- b[bin >= 1 & bin <= n_bins,
           .(cov = sum(cov * w, na.rm = TRUE) / sum(w, na.rm = TRUE)),
           by = bin]
  cov_mat[agg$bin, i] <- agg$cov

  if (i %% 20 == 0) message("[q03]   ", i, " / ", length(files), " samples")
}

# Per-bin summaries across samples
mean_cov   <- rowMeans(cov_mat, na.rm = TRUE)
median_cov <- apply(cov_mat, 1, stats::median, na.rm = TRUE)
cv_across  <- apply(cov_mat, 1, function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x); if (m == 0) return(NA_real_)
  stats::sd(x) / m
})
# Chromosome-wide per-sample median coverage → threshold per sample
per_sample_med <- apply(cov_mat, 2, stats::median, na.rm = TRUE)
thresh_per_sample <- LOW_FRAC * per_sample_med
n_low <- vapply(seq_len(n_bins), function(bi) {
  row <- cov_mat[bi, ]
  sum(row < thresh_per_sample, na.rm = TRUE)
}, integer(1))

bin_start_bp <- (seq_len(n_bins) - 1L) * BIN_BP + 1L
bin_end_bp   <- seq_len(n_bins) * BIN_BP
bin_mid_mb   <- round(((bin_start_bp + bin_end_bp) / 2) / 1e6, 4)

keep <- is.finite(mean_cov)
out_dt <- data.table(
  chrom              = CHROM,
  bin_start_bp       = bin_start_bp[keep],
  bin_end_bp         = bin_end_bp[keep],
  bin_mid_mb         = bin_mid_mb[keep],
  mean_cov           = round(mean_cov[keep], 3),
  median_cov         = round(median_cov[keep], 3),
  cv_across_samples  = round(cv_across[keep], 4),
  n_samples_low_cov  = as.integer(n_low[keep])
)

fwrite(out_dt, OUT, sep = "\t", quote = FALSE)
message("[q03] Wrote ", nrow(out_dt), " bins to ", OUT)
