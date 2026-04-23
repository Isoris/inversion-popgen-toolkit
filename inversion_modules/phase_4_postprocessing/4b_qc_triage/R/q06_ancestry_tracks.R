#!/usr/bin/env Rscript
# =============================================================================
# q06_ancestry_tracks.R  (Engine B native columns, K-agnostic)
# =============================================================================
# Reads Engine B per-chromosome local_Q caches and emits compact per-window +
# per-sample tracks for Q04 / Q07 / Q08 / scrubber consumption.
#
# Engine B summary columns:
#   window_id chrom start_bp end_bp n_sites
#   mean_delta12 mean_entropy mean_ena sd_delta12
# Engine B per-sample columns:
#   window_id chrom start_bp end_bp sample_id sample_idx
#   Q1..QK max_q delta12 entropy ena assigned_pop nQ_above_005 nQ_above_010
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

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

message("[q06] Loading window summary: ", SUMMARY)
s <- fread(SUMMARY)
setnames(s, tolower(names(s)))
if ("chrom" %in% names(s)) s <- s[chrom == CHROM]

required <- c("start_bp", "end_bp", "mean_delta12", "mean_entropy", "mean_ena")
missing  <- setdiff(required, names(s))
if (length(missing) > 0) {
  stop("Summary missing required columns: ", paste(missing, collapse = ", "),
       "\nAvailable: ", paste(names(s), collapse = ", "))
}
s[, window_mid_bp := as.integer((start_bp + end_bp) / 2)]

# ---- Per-sample cache -------------------------------------------------------
sample_dt <- NULL
if (!is.null(SAMPLES) && file.exists(SAMPLES)) {
  message("[q06] Loading per-sample cache: ", SAMPLES)
  sample_dt <- fread(SAMPLES)
  setnames(sample_dt, tolower(names(sample_dt)))
  if ("chrom" %in% names(sample_dt)) sample_dt <- sample_dt[chrom == CHROM]

  # Normalize sample-ID column name to "sample"
  sample_col_aliases <- c("sample", "sample_id", "sampleid", "ind", "individual", "id")
  found_sample <- intersect(sample_col_aliases, names(sample_dt))
  if (length(found_sample) > 0 && found_sample[1] != "sample") {
    setnames(sample_dt, found_sample[1], "sample")
  }

  # Derive maxq_label from assigned_pop (Engine B writes integer K index)
  if (!"maxq_label" %in% names(sample_dt)) {
    if ("assigned_pop" %in% names(sample_dt)) {
      sample_dt[, maxq_label := paste0("K", assigned_pop)]
    } else {
      qcols <- grep("^q[0-9]+$", names(sample_dt), value = TRUE)
      if (length(qcols) >= 2) {
        qm <- as.matrix(sample_dt[, ..qcols])
        sample_dt[, maxq_label := paste0("K", max.col(qm, ties.method = "first"))]
      }
    }
  }
}

# ---- CV of delta12 across samples, per window -------------------------------
if (!is.null(sample_dt) && all(c("window_id", "delta12") %in% names(sample_dt))) {
  cv_dt <- sample_dt[, .(
    cv_delta12_across_samples = {
      x <- delta12; m <- mean(x, na.rm = TRUE)
      if (is.na(m) || m == 0) NA_real_ else stats::sd(x, na.rm = TRUE) / m
    }
  ), by = window_id]
  s <- merge(s, cv_dt, by = "window_id", all.x = TRUE, sort = FALSE)
}
if (!"cv_delta12_across_samples" %in% names(s)) s[, cv_delta12_across_samples := NA_real_]

# ---- Majority-vote maxQ_label per window ------------------------------------
if (!is.null(sample_dt) && all(c("window_id", "maxq_label") %in% names(sample_dt))) {
  dom <- sample_dt[, .(
    maxQ_label = {
      tt <- table(maxq_label, useNA = "no")
      if (length(tt) == 0) NA_character_ else names(tt)[which.max(tt)]
    }
  ), by = window_id]
  s <- merge(s, dom, by = "window_id", all.x = TRUE, sort = FALSE)
}
if (!"maxQ_label" %in% names(s)) s[, maxQ_label := NA_character_]

# ---- Emit per-window track --------------------------------------------------
out_win_dt <- data.table(
  chrom                     = CHROM,
  window_id                 = s$window_id,
  window_start_bp           = as.integer(s$start_bp),
  window_end_bp             = as.integer(s$end_bp),
  window_mid_mb             = round(s$window_mid_bp / 1e6, 4),
  delta12                   = round(s$mean_delta12, 5),
  entropy                   = round(s$mean_entropy, 5),
  ena                       = round(s$mean_ena, 4),
  maxQ_label                = s$maxQ_label,
  cv_delta12_across_samples = round(s$cv_delta12_across_samples, 4)
)
setorder(out_win_dt, window_start_bp)
fwrite(out_win_dt, OUT_WIN, sep = "\t", quote = FALSE)
message("[q06] Wrote ", nrow(out_win_dt), " window rows -> ", OUT_WIN)

# ---- Emit per-sample track + wide maxQ matrix -------------------------------
if (!is.null(sample_dt) && "sample" %in% names(sample_dt) &&
    nrow(sample_dt) > 0 && !is.null(OUT_SAMP)) {
  keep <- intersect(
    c("sample", "chrom", "window_id", "max_q", "maxq", "maxq_label",
      "assigned_pop", "secondq", "delta12", "delta13",
      "entropy", "effective_num_ancestries", "ena"),
    names(sample_dt))
  out_samp <- sample_dt[, ..keep]

  rn <- names(out_samp)
  rn <- sub("^max_q$",      "maxQ",       rn)
  rn <- sub("^maxq$",       "maxQ",       rn)
  rn <- sub("^maxq_label$", "maxQ_label", rn)
  rn <- sub("^secondq$",    "secondQ",    rn)
  setnames(out_samp, names(out_samp), rn)

  mids <- s[, .(window_id, window_mid_bp)]
  out_samp <- merge(out_samp, mids, by = "window_id", all.x = TRUE, sort = FALSE)
  setorder(out_samp, sample, window_mid_bp)
  fwrite(out_samp, OUT_SAMP, sep = "\t", quote = FALSE, compress = "gzip")
  message("[q06] Wrote ", nrow(out_samp), " per-sample rows -> ", OUT_SAMP)

  if ("maxQ_label" %in% names(out_samp)) {
    wide <- dcast(out_samp, window_mid_bp ~ sample, value.var = "maxQ_label")
    setorder(wide, window_mid_bp)
    out_wide <- sub("\\.tsv\\.gz$", ".maxQ_wide.tsv.gz", OUT_SAMP)
    fwrite(wide, out_wide, sep = "\t", quote = FALSE, compress = "gzip")
    message("[q06] Wrote maxQ wide (", nrow(wide), "x", ncol(wide)-1, ") -> ", out_wide)
  }
}

message("[q06] DONE.")
