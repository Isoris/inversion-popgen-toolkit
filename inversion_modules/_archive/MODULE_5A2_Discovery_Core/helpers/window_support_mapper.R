#!/usr/bin/env Rscript

# =============================================================================
# window_support_mapper.R  (v8.0)
#
# TRACK-SPECIFIC SUPPORT SPAN BUILDER
#
# Uses the focal dense master window grid (from STEP09b windows_master.tsv.gz)
# as the positional backbone. For each focal window and each track
# (snake1/snake2/snake3), determines the effective support span.
#
# Snake 1: uses the exact focal SNP-defined window (no expansion needed)
# Snake 2: expands around focal center until enough grouping markers exist
# Snake 3: expands around focal center until enough phased bp / blocks exist
#
# PASS / WEAK / FAIL QC after support aggregation.
# Over-expanded windows get confidence downgrade even if thresholds are met.
#
# INPUTS:
#   <windows_master.tsv.gz>  — master window registry from STEP09b
#   <dosage_dir>             — per-chr dosage (for Snake 2 marker counting)
#   <phased_dir>             — per-chr phased summaries (for Snake 3 block counting)
#
# OUTPUT:
#   window_support_map.tsv.gz — one row per focal_window × track combination
#
# Usage:
#   Rscript window_support_mapper.R <windows_master> <outdir> \
#     [--dosage_dir <path>] [--phased_dir <path>] \
#     [--s2_min_markers 15] [--s3_min_phased_bp 5000] \
#     [--max_expansion_bp 500000]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript window_support_mapper.R <windows_master> <outdir> ...")
}

windows_file <- args[1]
outdir       <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dosage_dir <- NULL
phased_dir <- NULL
S2_MIN_MARKERS   <- 15L
S3_MIN_PHASED_BP <- 5000L
MAX_EXPANSION_BP <- 500000L
EXPANSION_STEP   <- 25000L
CONFIDENCE_DOWNGRADE_AT_BP <- 200000L

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--dosage_dir" && i < length(args)) {
    dosage_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--phased_dir" && i < length(args)) {
    phased_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--s2_min_markers" && i < length(args)) {
    S2_MIN_MARKERS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--s3_min_phased_bp" && i < length(args)) {
    S3_MIN_PHASED_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--max_expansion_bp" && i < length(args)) {
    MAX_EXPANSION_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

master <- fread(windows_file)
message("[SUPPORT_MAP] Master windows: ", nrow(master))
message("[SUPPORT_MAP] Expansion step: ", EXPANSION_STEP, " bp, max: ", MAX_EXPANSION_BP, " bp")

# =============================================================================
# BUILD SUPPORT MAP
# =============================================================================

support_rows <- list()
support_id <- 0L

for (chr in unique(master$chrom)) {
  chr_windows <- master[chrom == chr]
  n_win <- nrow(chr_windows)
  message("[SUPPORT_MAP] ", chr, ": ", n_win, " windows")

  # Load dosage sites for Snake 2 marker counting
  sites <- NULL
  if (!is.null(dosage_dir)) {
    sf <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))
    if (file.exists(sf)) sites <- fread(sf)
  }

  # Load phased summaries for Snake 3
  phased <- NULL
  if (!is.null(phased_dir)) {
    pf <- file.path(phased_dir, paste0(chr, ".phased_het_summary.tsv.gz"))
    if (file.exists(pf)) phased <- fread(pf)
  }

  for (wi in seq_len(n_win)) {
    wid        <- chr_windows$window_id[wi]
    focal_start <- chr_windows$start_bp[wi]
    focal_end   <- chr_windows$end_bp[wi]
    focal_center <- chr_windows$center_bp[wi]

    # ── Snake 1: exact focal span (no expansion) ─────────────────────
    support_id <- support_id + 1L
    support_rows[[support_id]] <- data.table(
      support_id = support_id, focal_window_id = wid,
      support_track = "snake1", chrom = chr,
      focal_start_bp = focal_start, focal_end_bp = focal_end,
      effective_start_bp = focal_start, effective_end_bp = focal_end,
      effective_span_bp = focal_end - focal_start,
      support_mode = "exact", expansion_steps = 0L,
      n_markers_focal = chr_windows$n_snps[wi],
      n_markers_effective = chr_windows$n_snps[wi],
      support_qc = "PASS",
      confidence_downgrade = 0.0
    )

    # ── Snake 2: expand until enough markers ─────────────────────────
    if (!is.null(sites) && nrow(sites) > 0) {
      eff_start <- focal_start
      eff_end   <- focal_end
      n_markers_focal <- sum(sites$pos >= focal_start & sites$pos <= focal_end)
      n_markers <- n_markers_focal
      exp_steps <- 0L

      while (n_markers < S2_MIN_MARKERS && (eff_end - eff_start) < MAX_EXPANSION_BP) {
        eff_start <- eff_start - EXPANSION_STEP
        eff_end   <- eff_end + EXPANSION_STEP
        n_markers <- sum(sites$pos >= eff_start & sites$pos <= eff_end)
        exp_steps <- exp_steps + 1L
      }

      span <- eff_end - eff_start
      if (n_markers >= S2_MIN_MARKERS) {
        qc <- if (span <= CONFIDENCE_DOWNGRADE_AT_BP) "PASS" else "WEAK"
      } else {
        qc <- "FAIL"
      }
      downgrade <- max(0, min(1, (span - (focal_end - focal_start)) / MAX_EXPANSION_BP))

      support_id <- support_id + 1L
      support_rows[[support_id]] <- data.table(
        support_id = support_id, focal_window_id = wid,
        support_track = "snake2", chrom = chr,
        focal_start_bp = focal_start, focal_end_bp = focal_end,
        effective_start_bp = eff_start, effective_end_bp = eff_end,
        effective_span_bp = span,
        support_mode = if (exp_steps == 0) "exact" else paste0("expanded_", exp_steps, "x"),
        expansion_steps = exp_steps,
        n_markers_focal = n_markers_focal,
        n_markers_effective = n_markers,
        support_qc = qc,
        confidence_downgrade = round(downgrade, 4)
      )
    }

    # ── Snake 3: expand until enough phased bp ───────────────────────
    if (!is.null(phased) && nrow(phased) > 0) {
      eff_start <- focal_start
      eff_end   <- focal_end

      blocks <- unique(phased[, .(block_start, block_end)])
      focal_phased <- sum(pmin(blocks$block_end, focal_end) -
                          pmax(blocks$block_start, focal_start))
      focal_phased <- max(0L, focal_phased)
      phased_bp <- focal_phased
      exp_steps <- 0L

      while (phased_bp < S3_MIN_PHASED_BP && (eff_end - eff_start) < MAX_EXPANSION_BP) {
        eff_start <- eff_start - EXPANSION_STEP
        eff_end   <- eff_end + EXPANSION_STEP
        phased_bp <- sum(pmax(0L, pmin(blocks$block_end, eff_end) -
                                   pmax(blocks$block_start, eff_start)))
        exp_steps <- exp_steps + 1L
      }

      span <- eff_end - eff_start
      if (phased_bp >= S3_MIN_PHASED_BP) {
        qc <- if (span <= CONFIDENCE_DOWNGRADE_AT_BP) "PASS" else "WEAK"
      } else {
        qc <- "FAIL"
      }
      downgrade <- max(0, min(1, (span - (focal_end - focal_start)) / MAX_EXPANSION_BP))

      support_id <- support_id + 1L
      support_rows[[support_id]] <- data.table(
        support_id = support_id, focal_window_id = wid,
        support_track = "snake3", chrom = chr,
        focal_start_bp = focal_start, focal_end_bp = focal_end,
        effective_start_bp = eff_start, effective_end_bp = eff_end,
        effective_span_bp = span,
        support_mode = if (exp_steps == 0) "exact" else paste0("expanded_", exp_steps, "x"),
        expansion_steps = exp_steps,
        n_markers_focal = as.integer(focal_phased),
        n_markers_effective = as.integer(phased_bp),
        support_qc = qc,
        confidence_downgrade = round(downgrade, 4)
      )
    }
  }
}

# =============================================================================
# WRITE
# =============================================================================

support_map <- if (length(support_rows) > 0) rbindlist(support_rows) else {
  data.table(support_id = integer(), focal_window_id = integer(),
             support_track = character(), support_qc = character())
}

f1 <- file.path(outdir, "window_support_map.tsv.gz")
fwrite(support_map, f1, sep = "\t")

message("\n[DONE] Window support mapper complete")
message("  ", f1, " (", nrow(support_map), " entries)")

if (nrow(support_map) > 0) {
  qc_tab <- support_map[, .N, by = .(support_track, support_qc)]
  for (i in seq_len(nrow(qc_tab))) {
    message("  ", qc_tab$support_track[i], " ", qc_tab$support_qc[i], ": ", qc_tab$N[i])
  }
}
