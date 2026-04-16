#!/usr/bin/env Rscript

# =============================================================================
# STEP_HOBS_confirmation_module.R
#
# SECONDARY CONFIRMATION MINI-MODULE for inversions.
# Adapted from Claire Mérot's Hobs/HWE sliding-window logic
# (https://github.com/clairemerot/angsd_pipeline) and made subset-aware
# and multiscale for structured hatchery catfish data.
#
# IMPORTANT: This is NOT primary inversion discovery.
# Use only as secondary confirmation / interpretation / support.
# See README for rationale.
#
# Computes:
#   1. Site-level Hexp, Hobs, F from ANGSD -doHWE output
#   2. Multiscale window summaries with mean, median, outlier burden
#   3. Candidate interval overlays
#
# Usage:
#   Rscript STEP_HOBS_confirmation_module.R \
#     <hwe_file> <subset_label> <outdir> \
#     [candidate_table=NONE] [chrom_col=Chromo] [pos_col=Position] \
#     [freq_col=Freq] [f_col=F] [hwefreq_col=hweFreq]
#
# Input: ANGSD -doHWE output (tab-separated with header)
#
# Citation note:
#   Hobs/HWE sliding-window logic: Claire Mérot et al.
#   ANGSD/thetaStat: Korneliussen et al.
#   Adapted for hatchery catfish workflow.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP_HOBS_confirmation_module.R <hwe_file> <subset_label> <outdir> [candidate_table=NONE] ...")
}

hwe_file       <- args[1]
subset_label   <- args[2]
outdir         <- args[3]
candidate_file <- if (length(args) >= 4 && args[4] != "NONE") args[4] else NA_character_
chrom_col      <- if (length(args) >= 5) args[5] else "Chromo"
pos_col        <- if (length(args) >= 6) args[6] else "Position"
freq_col       <- if (length(args) >= 7) args[7] else "Freq"
f_col          <- if (length(args) >= 8) args[8] else "F"
hwefreq_col    <- if (length(args) >= 9) args[9] else "hweFreq"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Read HWE data ──────────────────────────────────────────────────────────
message("[INFO] Reading HWE file: ", hwe_file)
hwe <- fread(hwe_file)

# Standardize column names
if (chrom_col != "chrom" && chrom_col %in% names(hwe)) setnames(hwe, chrom_col, "chrom")
if (pos_col != "position" && pos_col %in% names(hwe)) setnames(hwe, pos_col, "position")
if (f_col != "F_val" && f_col %in% names(hwe)) setnames(hwe, f_col, "F_val")
if (hwefreq_col %in% names(hwe)) setnames(hwe, hwefreq_col, "hweFreq")
if (freq_col %in% names(hwe) && freq_col != "Freq") setnames(hwe, freq_col, "Freq")

stopifnot(all(c("chrom", "position") %in% names(hwe)))

message("[INFO] Sites: ", nrow(hwe), " across ", uniqueN(hwe$chrom), " chromosomes")

# ── Compute site-level Hexp, Hobs ────────────────────────────────────────
# Following Claire Mérot logic:
# Hexp = 2*p*(1-p)  where p = hweFreq
# Hobs = Hexp * (1 - F)
if ("hweFreq" %in% names(hwe)) {
  hwe[, Hexp := 2 * hweFreq * (1 - hweFreq)]
} else if ("Freq" %in% names(hwe)) {
  hwe[, Hexp := 2 * Freq * (1 - Freq)]
} else {
  stop("No frequency column found (expected hweFreq or Freq)")
}

if ("F_val" %in% names(hwe)) {
  hwe[, Hobs := Hexp * (1 - F_val)]
} else {
  stop("No F column found")
}

# Optional: minor/major homozygote proportions
if ("Freq" %in% names(hwe) && "F_val" %in% names(hwe)) {
  hwe[, Hminor := Freq^2 + Freq * (1 - Freq) * F_val]
  hwe[, Hmajor := (1 - Freq)^2 + Freq * (1 - Freq) * F_val]
}

hwe[, subset := subset_label]

# ── Write site-level output ──────────────────────────────────────────────
site_out <- file.path(outdir, paste0(subset_label, ".site_level_Hobs_F.tsv.gz"))
fwrite(hwe, site_out, sep = "\t")
message("[INFO] Wrote site-level: ", site_out)

# ── Multiscale window summaries ──────────────────────────────────────────
# Window definitions: size_bp, step_bp, label
windows <- list(
  list(5000,    1000,   "5kb"),
  list(10000,   2000,   "10kb"),
  list(50000,   10000,  "50kb"),
  list(100000,  20000,  "100kb"),
  list(250000,  50000,  "250kb"),
  list(500000,  100000, "500kb"),
  list(1000000, 200000, "1Mb")
)

# Genome-wide robust thresholds (within subset)
genome_median_Hobs <- median(hwe$Hobs, na.rm = TRUE)
genome_mad_Hobs    <- mad(hwe$Hobs, na.rm = TRUE)
genome_median_F    <- median(hwe$F_val, na.rm = TRUE)
genome_mad_F       <- mad(hwe$F_val, na.rm = TRUE)

# Define outlier thresholds: median ± 3*MAD
thresh_low_Hobs  <- genome_median_Hobs - 3 * genome_mad_Hobs
thresh_high_Hobs <- genome_median_Hobs + 3 * genome_mad_Hobs
thresh_low_F     <- genome_median_F - 3 * genome_mad_F
thresh_high_F    <- genome_median_F + 3 * genome_mad_F

# Tag outlier sites
hwe[, outlier_low_Hobs  := Hobs < thresh_low_Hobs]
hwe[, outlier_high_Hobs := Hobs > thresh_high_Hobs]
hwe[, outlier_low_F     := F_val < thresh_low_F]
hwe[, outlier_high_F    := F_val > thresh_high_F]

for (wdef in windows) {
  win_size <- wdef[[1]]
  win_step <- wdef[[2]]
  win_label <- wdef[[3]]

  message("[INFO] Computing windows: ", win_label, " (size=", win_size, ", step=", win_step, ")")

  # Assign each site to its window
  hwe[, win_start := floor(position / win_step) * win_step]

  win_dt <- hwe[, .(
    window_start = win_start[1],
    window_end   = win_start[1] + win_size,
    window_center = win_start[1] + win_size / 2,
    n_sites      = .N,
    mean_Hobs    = round(mean(Hobs, na.rm = TRUE), 6),
    median_Hobs  = round(median(Hobs, na.rm = TRUE), 6),
    sd_Hobs      = round(sd(Hobs, na.rm = TRUE), 6),
    min_Hobs     = round(min(Hobs, na.rm = TRUE), 6),
    max_Hobs     = round(max(Hobs, na.rm = TRUE), 6),
    mean_F       = round(mean(F_val, na.rm = TRUE), 6),
    median_F     = round(median(F_val, na.rm = TRUE), 6),
    sd_F         = round(sd(F_val, na.rm = TRUE), 6),
    min_F        = round(min(F_val, na.rm = TRUE), 6),
    max_F        = round(max(F_val, na.rm = TRUE), 6),
    n_low_Hobs_outlier  = sum(outlier_low_Hobs, na.rm = TRUE),
    n_high_Hobs_outlier = sum(outlier_high_Hobs, na.rm = TRUE),
    n_low_F_outlier     = sum(outlier_low_F, na.rm = TRUE),
    n_high_F_outlier    = sum(outlier_high_F, na.rm = TRUE)
  ), by = .(chrom, win_start)]

  win_dt[, frac_low_Hobs_outlier  := round(n_low_Hobs_outlier / n_sites, 4)]
  win_dt[, frac_high_Hobs_outlier := round(n_high_Hobs_outlier / n_sites, 4)]
  win_dt[, frac_low_F_outlier     := round(n_low_F_outlier / n_sites, 4)]
  win_dt[, frac_high_F_outlier    := round(n_high_F_outlier / n_sites, 4)]
  win_dt[, subset := subset_label]

  win_out <- file.path(outdir, paste0(subset_label, ".win", win_label, ".step",
                                       format(win_step, scientific = FALSE), ".tsv.gz"))
  fwrite(win_dt, win_out, sep = "\t")

  # Clean up temp column
  hwe[, win_start := NULL]
}

# ── Candidate interval overlay ──────────────────────────────────────────
if (!is.na(candidate_file) && file.exists(candidate_file)) {
  message("[INFO] Computing candidate interval overlay")
  cand_dt <- fread(candidate_file)
  stopifnot(all(c("candidate_id", "chrom", "start_bp", "end_bp") %in% names(cand_dt)))

  overlay_list <- list()
  for (i in seq_len(nrow(cand_dt))) {
    cr <- cand_dt[i]
    sub <- hwe[chrom == cr$chrom & position >= cr$start_bp & position <= cr$end_bp]
    if (nrow(sub) == 0) next

    overlay_list[[length(overlay_list) + 1]] <- data.table(
      candidate_id   = cr$candidate_id,
      chrom          = cr$chrom,
      start_bp       = cr$start_bp,
      end_bp         = cr$end_bp,
      subset         = subset_label,
      n_sites        = nrow(sub),
      mean_F         = round(mean(sub$F_val, na.rm = TRUE), 6),
      median_F       = round(median(sub$F_val, na.rm = TRUE), 6),
      mean_Hobs      = round(mean(sub$Hobs, na.rm = TRUE), 6),
      median_Hobs    = round(median(sub$Hobs, na.rm = TRUE), 6),
      frac_low_Hobs  = round(mean(sub$outlier_low_Hobs, na.rm = TRUE), 4),
      frac_high_Hobs = round(mean(sub$outlier_high_Hobs, na.rm = TRUE), 4),
      frac_high_F    = round(mean(sub$outlier_high_F, na.rm = TRUE), 4)
    )
  }

  if (length(overlay_list) > 0) {
    overlay_dt <- rbindlist(overlay_list)
    overlay_out <- file.path(outdir, paste0(subset_label, ".candidate_Hobs_overlay.tsv"))
    fwrite(overlay_dt, overlay_out, sep = "\t")
    message("[INFO] Wrote candidate overlay: ", overlay_out)
  }
}

message("[DONE] Hobs/HWE confirmation module for subset: ", subset_label)
message("  Methods note: HWE/Hobs-based scans were treated as a secondary,")
message("  subset-based diagnostic, not primary inversion detection.")
