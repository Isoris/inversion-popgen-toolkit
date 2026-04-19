#!/usr/bin/env Rscript
# =============================================================================
# q06_ancestry_tracks.R
#
# Parse the Engine B local_Q cache files for one chromosome and emit two
# compact tracks: per-window (always) and per-sample (if available).
#
# Summary cache columns (from instant_q / ancestry_bridge):
#   chrom, window_id or window_mid_bp or start/end, delta12, entropy, ena,
#   Q1..QK (means across samples in the window), optionally maxQ_label.
#   Column naming varies slightly across versions, so we auto-detect.
#
# Samples cache columns:
#   chrom, window_mid_bp, sample, Q1..QK, maxQ, maxQ_label (optional),
#   secondQ, delta12, delta13, entropy, effective_num_ancestries
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
SUMMARY  <- get_arg("--summary")
SAMPLES  <- get_arg("--samples", NULL)
CHROM    <- get_arg("--chrom")
OUT_WIN  <- get_arg("--out_win")
OUT_SAMP <- get_arg("--out_samp")
stopifnot(!is.null(SUMMARY), !is.null(CHROM), !is.null(OUT_WIN))

# ---- Summary (per-window, always produced) -----------------------------------
message("[q06] Loading window summary: ", SUMMARY)
s <- fread(SUMMARY)
# Normalize names
setnames(s, tolower(names(s)))

# Figure out coordinates
if ("chrom" %in% names(s) && CHROM %in% unique(s$chrom)) {
  s <- s[chrom == CHROM]
}
# Accept several coordinate schemas
start_col <- intersect(c("window_start_bp", "start_bp", "start", "winstart"),
                       names(s))[1]
end_col   <- intersect(c("window_end_bp", "end_bp", "end", "winstop"),
                       names(s))[1]
mid_col   <- intersect(c("window_mid_bp", "mid_bp", "wincenter"),
                       names(s))[1]
if (is.null(start_col) && !is.null(mid_col)) {
  s[, window_start_bp := get(mid_col)]
  start_col <- "window_start_bp"
}
if (is.null(end_col) && !is.null(mid_col)) {
  s[, window_end_bp := get(mid_col)]
  end_col <- "window_end_bp"
}
if (is.null(mid_col)) {
  s[, window_mid_bp := as.integer((get(start_col) + get(end_col)) / 2)]
  mid_col <- "window_mid_bp"
}
# Required metrics
for (col in c("delta12", "entropy")) {
  if (!col %in% names(s)) stop("Summary missing column: ", col)
}
if (!"ena" %in% names(s)) s[, ena := exp(entropy)]
if (!"maxq_label" %in% names(s)) {
  # Derive from Q columns if present
  qcols <- grep("^q[0-9]+$", names(s), value = TRUE)
  if (length(qcols) >= 2) {
    qmat <- as.matrix(s[, ..qcols])
    s[, maxq_label := paste0("K", max.col(qmat, ties.method = "first"))]
  } else {
    s[, maxq_label := NA_character_]
  }
}

# CV of delta12 across samples — requires samples cache
cv_delta <- NA_real_
sample_dt <- NULL
if (!is.null(SAMPLES) && file.exists(SAMPLES)) {
  message("[q06] Loading per-sample cache: ", SAMPLES)
  sample_dt <- fread(SAMPLES)
  setnames(sample_dt, tolower(names(sample_dt)))
  if ("chrom" %in% names(sample_dt) && CHROM %in% unique(sample_dt$chrom)) {
    sample_dt <- sample_dt[chrom == CHROM]
  }
  # Normalize coords on sample table
  if (!"window_mid_bp" %in% names(sample_dt)) {
    mc <- intersect(c("mid_bp", "wincenter", "window_start_bp", "start_bp"),
                    names(sample_dt))[1]
    if (!is.null(mc)) setnames(sample_dt, mc, "window_mid_bp")
  }
  # Per-window CV of delta12 across samples
  cv_dt <- sample_dt[, .(
    cv_delta12 = {
      x <- delta12
      m <- mean(x, na.rm = TRUE)
      if (is.na(m) || m == 0) NA_real_ else stats::sd(x, na.rm = TRUE) / m
    }
  ), by = window_mid_bp]
  # Merge
  s <- merge(s, cv_dt, by.x = mid_col, by.y = "window_mid_bp",
             all.x = TRUE, sort = FALSE)
  setnames(s, "cv_delta12", "cv_delta12_across_samples")
} else {
  s[, cv_delta12_across_samples := NA_real_]
}

# ---- Emit per-window track ---------------------------------------------------
out_win_dt <- data.table(
  chrom           = CHROM,
  window_start_bp = as.integer(s[[start_col]]),
  window_end_bp   = as.integer(s[[end_col]]),
  window_mid_mb   = round(s[[mid_col]] / 1e6, 4),
  delta12         = round(s$delta12, 5),
  entropy         = round(s$entropy, 5),
  ena             = round(s$ena, 4),
  maxQ_label      = s$maxq_label,
  cv_delta12_across_samples = round(s$cv_delta12_across_samples, 4)
)
setorder(out_win_dt, window_start_bp)
fwrite(out_win_dt, OUT_WIN, sep = "\t", quote = FALSE)
message("[q06] Wrote ", nrow(out_win_dt), " window rows -> ", OUT_WIN)

# ---- Emit per-sample track (if samples cache was available) -----------------
if (!is.null(sample_dt) && nrow(sample_dt) > 0 && !is.null(OUT_SAMP)) {
  # Trim to useful columns
  keep_cols <- intersect(c("sample", "chrom", "window_mid_bp",
                           "maxq", "maxq_label", "secondq",
                           "delta12", "delta13", "entropy",
                           "effective_num_ancestries", "ena"),
                         names(sample_dt))
  out_samp <- sample_dt[, ..keep_cols]
  setnames(out_samp, names(out_samp), sub("maxq", "maxQ", names(out_samp)))
  setnames(out_samp, names(out_samp), sub("secondq", "secondQ", names(out_samp)))
  setnames(out_samp, names(out_samp), sub("maxQ_label", "maxQ_label", names(out_samp)))
  setorder(out_samp, sample, window_mid_bp)
  fwrite(out_samp, OUT_SAMP, sep = "\t", quote = FALSE, compress = "gzip")
  message("[q06] Wrote ", nrow(out_samp), " per-sample rows -> ", OUT_SAMP)

  # Also write a compact wide matrix: one row per window, columns = samples,
  # values = maxQ_label (K1..K8). This is what the scrubber consumes directly.
  if ("maxQ_label" %in% names(out_samp)) {
    # Pivot
    wide <- dcast(out_samp, window_mid_bp ~ sample, value.var = "maxQ_label")
    setorder(wide, window_mid_bp)
    out_wide <- sub("\\.tsv\\.gz$", ".maxQ_wide.tsv.gz", OUT_SAMP)
    fwrite(wide, out_wide, sep = "\t", quote = FALSE, compress = "gzip")
    message("[q06] Wrote maxQ wide matrix (", nrow(wide), " windows x ",
            ncol(wide) - 1, " samples) -> ", out_wide)
  }
}
