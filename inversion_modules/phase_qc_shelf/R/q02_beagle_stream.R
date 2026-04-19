#!/usr/bin/env Rscript
# =============================================================================
# q02_beagle_stream.R
#
# Stream a BEAGLE genotype-likelihood file in chunks and compute per-bin
# posterior uncertainty metrics. Handles large files (multi-GB) without
# loading the whole thing into RAM.
#
# BEAGLE format (post-ANGSD):
#   marker  allele1  allele2  P00_s1 P01_s1 P11_s1  P00_s2 P01_s2 P11_s2 ...
#   Header row present. The `marker` column in main_qcpass has format
#   "<chrom>_<pos>" for the combined-genome file, but per-chrom files
#   typically keep marker as just a chromosome-unique ID. We join on the
#   companion .pos file for authoritative positions.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- CLI ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
BEAGLE <- get_arg("--beagle")
POS    <- get_arg("--pos")
CHROM  <- get_arg("--chrom")
BIN_MB <- as.numeric(get_arg("--bin_mb", "0.05"))
THRESH <- as.numeric(get_arg("--thresh", "0.9"))
OUT    <- get_arg("--out")
CHUNK  <- as.integer(get_arg("--chunk", "20000"))  # sites per chunk

stopifnot(!is.null(BEAGLE), file.exists(BEAGLE))
stopifnot(!is.null(POS),    file.exists(POS))
stopifnot(!is.null(CHROM), !is.null(OUT))

BIN_BP <- as.integer(BIN_MB * 1e6)
message("[q02] ", CHROM, "   bin=", BIN_MB, "Mb   thresh=", THRESH,
        "   chunk=", CHUNK)

# ---- Load position file ------------------------------------------------------
# .pos files: 2 cols typically (chrom, pos) but sometimes tab-separated with header
pos_dt <- fread(POS, header = FALSE, col.names = c("chrom_pos", "position"),
                fill = TRUE)
# First row may be a header; detect by non-numeric position
if (!is.numeric(pos_dt$position) || is.na(pos_dt$position[1])) {
  pos_dt <- fread(POS, header = TRUE)
  # Try to normalize column names
  if ("pos" %in% names(pos_dt)) setnames(pos_dt, "pos", "position")
  if (!"position" %in% names(pos_dt)) {
    # Try column 2 as position
    setnames(pos_dt, names(pos_dt)[2], "position")
  }
}
pos_dt[, position := as.integer(position)]
# Keep only chrom rows if the file has multiple chroms (combined-genome case)
if ("chrom_pos" %in% names(pos_dt)) {
  pos_dt <- pos_dt[chrom_pos == CHROM | grepl(paste0("^", CHROM), chrom_pos)]
}
message("[q02] Pos rows: ", nrow(pos_dt))

# ---- Peek BEAGLE header to get sample count ---------------------------------
con <- gzfile(BEAGLE, "r")
hdr <- readLines(con, n = 1)
close(con)
hdr_fields <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
n_triplet_cols <- length(hdr_fields) - 3
if (n_triplet_cols %% 3 != 0) {
  stop("Header column count not divisible by 3 after leading metadata cols: ",
       length(hdr_fields))
}
n_samples <- n_triplet_cols / 3
message("[q02] Samples detected: ", n_samples)

# ---- Stream in chunks --------------------------------------------------------
# Accumulators keyed by bin id
bin_accum <- list()  # name = bin_id (as char), value = list()

# Per-sample running counts (uncertainty per bin per sample)
# We store incremental sums and counts, then summarize.
# Memory: n_bins_max × n_samples doubles. For LG01 ~120 Mb / 0.05 = 2400 bins × 226 = 542k doubles. Fine.

# Determine max_bin upfront
max_pos <- max(pos_dt$position, na.rm = TRUE)
n_bins  <- as.integer(ceiling(max_pos / BIN_BP)) + 1L
message("[q02] n_bins: ", n_bins)

# Wide matrices: rows=bins, cols=samples
site_count <- integer(n_bins)
sum_max_post <- numeric(n_bins)                             # sum of max posterior, summed over sites x samples
sum_uncertain <- numeric(n_bins)                            # count of (site x sample) below thresh

# Per-sample per-bin uncertain counts and site counts (for CV + per-sample rate)
uncertain_mat <- matrix(0L, nrow = n_bins, ncol = n_samples)
# We reuse site_count for per-bin site totals (same for all samples per site)

con <- gzfile(BEAGLE, "r")
invisible(readLines(con, n = 1))  # skip header

site_idx <- 0L
repeat {
  lines <- readLines(con, n = CHUNK)
  if (length(lines) == 0) break
  chunk <- fread(text = lines, sep = "\t", header = FALSE)
  # Columns: marker allele1 allele2 [P00 P01 P11]*n_samples
  # Position: rows 1..n of this chunk correspond to pos_dt rows (site_idx+1 .. site_idx+n)
  n_here <- nrow(chunk)

  # Map to positions
  these_positions <- pos_dt$position[(site_idx + 1L):(site_idx + n_here)]
  bin_ids <- as.integer(floor((these_positions - 1L) / BIN_BP)) + 1L  # 1-based

  # Extract posterior matrix: drop 3 metadata cols, convert to numeric matrix
  post <- as.matrix(chunk[, -(1:3)])
  # Reshape: n_sites × (3 * n_samples) -> n_sites × 3 × n_samples
  # Compute per-site per-sample max posterior
  # Fastest: for each sample i, cols = 1+(i-1)*3 : i*3
  for (si in seq_len(n_samples)) {
    c_off <- (si - 1L) * 3L
    mx <- pmax.int(post[, c_off + 1L], post[, c_off + 2L], post[, c_off + 3L])
    unc <- mx < THRESH

    # Accumulate per bin
    # bin_ids is integer vector; use rowsum for speed? For chunk size 20k it's negligible.
    # Use tabulate for counts; for sums, use tapply or data.table
    tb_unc <- tabulate(bin_ids[unc], nbins = n_bins)
    uncertain_mat[, si] <- uncertain_mat[, si] + tb_unc

    # Global sums
    # Use vapply over unique bins? Faster: use .colSums via split — but stick to simple
    sum_max_post   <- sum_max_post   + as.numeric(tapply(mx,  bin_ids, sum, default = 0))[seq_len(n_bins)]
    sum_uncertain  <- sum_uncertain  + as.numeric(tb_unc)
  }

  # Global site count per bin (same across samples since we process all n_samples sites)
  tb_sites <- tabulate(bin_ids, nbins = n_bins)
  site_count <- site_count + tb_sites

  site_idx <- site_idx + n_here
  if (site_idx %% (CHUNK * 5L) == 0L) {
    message("[q02]   streamed ", site_idx, " sites")
  }
}
close(con)
message("[q02] Total sites streamed: ", site_idx)
if (site_idx != nrow(pos_dt)) {
  warning("site count mismatch: beagle=", site_idx, " pos=", nrow(pos_dt))
}

# ---- Summarize per bin -------------------------------------------------------
# Mean max_post per (site x sample) = sum_max_post / (n_sites * n_samples)
# Mean uncertain frac = sum_uncertain / (n_sites * n_samples)
# Per-sample uncertain rate = uncertain_mat[bin, si] / site_count[bin]
# CV across samples = sd(rates) / mean(rates)
# n_samples_high = count of samples with rate > 2 * bin mean

bin_mid_mb <- ((seq_len(n_bins) - 1L) * BIN_BP + BIN_BP / 2) / 1e6
bin_start_bp <- (seq_len(n_bins) - 1L) * BIN_BP + 1L
bin_end_bp <- seq_len(n_bins) * BIN_BP

keep <- site_count > 0
out_dt <- data.table(
  chrom         = CHROM,
  bin_start_bp  = bin_start_bp[keep],
  bin_end_bp    = bin_end_bp[keep],
  bin_mid_mb    = round(bin_mid_mb[keep], 4),
  n_sites       = site_count[keep]
)

denom <- site_count[keep] * n_samples
out_dt[, mean_max_post        := round(sum_max_post[keep] / denom, 4)]
out_dt[, mean_uncertain_frac  := round(sum_uncertain[keep] / denom, 4)]

# Per-sample rates in these bins
rates <- sweep(uncertain_mat[keep, , drop = FALSE],
               1, site_count[keep], "/")
out_dt[, cv_across_samples           := round(apply(rates, 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  sd(x, na.rm = TRUE) / m
}), 4)]
out_dt[, n_samples_with_high_uncertain := apply(rates, 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(0L)
  sum(x > 2 * m, na.rm = TRUE)
})]

fwrite(out_dt, OUT, sep = "\t", quote = FALSE)
message("[q02] Wrote ", nrow(out_dt), " non-empty bins to ", OUT)
