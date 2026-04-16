#!/usr/bin/env Rscript
# =============================================================================
# STEP_D12_bridge_to_C01d.R — Format adapter to C01d triangle-intervals schema
# =============================================================================
#
# WHAT THIS DOES:
#   Converts inv_detect_v9.3 output (staircase blocks + block scores + NN
#   persistence + peel diagnostic) into the exact file format that
#   STEP_C01d_candidate_scoring.R expects from STEP_C01c_triangle_regimes.R.
#
#   This is the BRIDGE that lets you replace C01c with inv_detect while
#   keeping C01d → C01e → ... → C01m unchanged.
#
# INPUT → OUTPUT MAPPING:
#   inv_detect_v9.3 file            → codebase file (C01c format)
#   ────────────────────────────────   ──────────────────────────────
#   blocks_<chr>.tsv                → triangle_intervals.tsv.gz
#   block_scores_<chr>.tsv           → (columns merged into intervals)
#   nn_persistence_<chr>.tsv        → (columns merged into intervals)
#   peel_diagnostic_<chr>.tsv       → (columns merged into intervals)
#   consensus_<chr>.tsv             → triangle_interval_comparison.tsv.gz
#   (staircase profiles)            → triangle_sample_composition.tsv.gz
#   (not applicable)                → triangle_bridges.tsv.gz (empty OK)
#   (not applicable)                → triangle_offdiag_linkage.tsv.gz (empty OK)
#   (bloc parent/child)             → triangle_subtriangles.tsv.gz
#   (not applicable)                → triangle_subregimes.tsv.gz (empty OK)
#
# HOW TO USE:
#   1. Run inv_detect_v9.3 (phases 1-9) → inv_detect_out_v9.3/
#   2. Run this bridge script → triangle_from_detector/
#   3. Point C01d at triangle_from_detector/ instead of the old C01c output
#
#   Rscript 12_bridge_to_codebase.R \
#     --detector_dir inv_detect_out_v9.3/ \
#     --precomp_dir precomp/ \
#     --outdir triangle_from_detector/
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
detector_dir <- NULL; precomp_dir <- NULL; outdir <- "triangle_from_detector"

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--detector_dir" && i < length(args)) { detector_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp_dir" && i < length(args)) { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))  { outdir <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

if (is.null(detector_dir)) stop("--detector_dir required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Find all chromosome block files ----
block_files <- list.files(detector_dir, pattern = "^blocks_.*\\.tsv$", full.names = TRUE)
if (length(block_files) == 0) stop("No blocks_*.tsv files in ", detector_dir)

message("[BRIDGE] Found ", length(block_files), " chromosome block files")

# ---- Process each chromosome ----
all_intervals <- list()
all_subtriangles <- list()
all_composition <- list()
all_comparison <- list()

for (bf in sort(block_files)) {
  chr <- sub("^blocks_", "", sub("\\.tsv$", "", basename(bf)))
  message("\n[BRIDGE] ", chr)

  blocks <- fread(bf)
  if (nrow(blocks) == 0) { message("  No blocks"); next }

  # ---- Load supplementary files if they exist ----
  score_file <- file.path(detector_dir, paste0("block_scores_", chr, ".tsv"))
  nn_file    <- file.path(detector_dir, paste0("nn_persistence_", chr, ".tsv"))
  peel_file  <- file.path(detector_dir, paste0("peel_diagnostic_", chr, ".tsv"))
  cons_file  <- file.path(detector_dir, paste0("consensus_", chr, ".tsv"))

  scores <- if (file.exists(score_file)) fread(score_file) else data.table()
  nn_ev  <- if (file.exists(nn_file))    fread(nn_file)    else data.table()
  peel   <- if (file.exists(peel_file))  fread(peel_file)  else data.table()
  cons   <- if (file.exists(cons_file))  fread(cons_file)  else data.table()

  # Raw-variant scores
  raw_scores <- if (nrow(scores) > 0 && "variant" %in% names(scores)) {
    scores[variant == "raw"]
  } else data.table()

  # ---- Build triangle_intervals format ----
  for (r in seq_len(nrow(blocks))) {
    b <- blocks[r]
    bid <- b$block_id

    # Merge scores
    sc <- if (nrow(raw_scores) > 0 && bid %in% raw_scores$block_id) {
      raw_scores[block_id == bid]
    } else data.table()

    inside_mean <- if (nrow(sc) > 0) sc$inside_mean[1] else b$height
    contrast    <- if (nrow(sc) > 0) sc$contrast[1] else b$step_down_left
    sharpness   <- if (nrow(sc) > 0) sc$sharpness[1] else NA_real_
    occupancy   <- if (nrow(sc) > 0) sc$occupancy[1] else NA_real_
    shape       <- if (nrow(sc) > 0) sc$shape_class[1] else "unknown"
    squareness  <- if (nrow(sc) > 0) sc$squareness[1] else NA_real_

    # Map shape_class → interval_type (C01c naming convention)
    interval_type <- switch(shape %||% "unknown",
      strong_square   = "strong_triangle",
      diffuse_square  = "moderate_triangle",
      diagonal_band   = "sharp_but_not_square",
      noise           = "weak_zone",
      ambiguous       = "diffuse_zone",
      "diffuse_zone"
    )

    # NN persistence
    nn_row <- if (nrow(nn_ev) > 0 && bid %in% nn_ev$candidate_id) {
      nn_ev[candidate_id == bid]
    } else data.table()

    survives_nn40 <- if (nrow(nn_row) > 0 && "survives_nn40" %in% names(nn_row)) {
      nn_row$survives_nn40[1]
    } else NA

    # Peel diagnostic (L1b)
    peel_row <- if (nrow(peel) > 0) {
      peel[block_id == bid & peel_mode == "L1b_chrlocal_kin"]
    } else data.table()

    l1b_effect <- if (nrow(peel_row) > 0) peel_row$effect_class[1] else NA_character_

    # Compute homogeneity (1 - patchiness if available)
    homogeneity <- if (nrow(sc) > 0 && "patchiness" %in% names(sc)) {
      max(0, 1 - sc$patchiness[1])
    } else NA_real_

    # Background median
    bg_median <- if (nrow(sc) > 0) sc$bg_median[1] else NA_real_

    all_intervals[[length(all_intervals) + 1]] <- data.table(
      chrom         = chr,
      interval_id   = bid,
      apex_bin      = round((b$start + b$end) / 2),
      L_bin         = b$start,
      R_bin         = b$end,
      n_bins        = b$width,
      start_bp      = b$start_bp,
      end_bp        = b$end_bp,
      start_mb      = b$start_mb,
      end_mb        = b$end_mb,
      inside_mean   = round(inside_mean %||% NA_real_, 4),
      bg_median     = round(bg_median %||% NA_real_, 4),
      contrast      = round(contrast %||% NA_real_, 4),
      sharpness     = round(sharpness %||% NA_real_, 4),
      homogeneity   = round(homogeneity %||% NA_real_, 3),
      bin_support   = round(occupancy %||% NA_real_, 3),
      apex_score    = round(inside_mean %||% NA_real_, 4),
      interval_type = interval_type,
      n_hot_subregimes = 0L,
      frac_hot      = 0,
      # Extra columns from detector (C01d ignores unknown columns gracefully)
      squareness    = round(squareness %||% NA_real_, 4),
      shape_class   = shape %||% "unknown",
      survives_nn40 = survives_nn40,
      l1b_peel_effect = l1b_effect,
      source        = "inv_detect_v9.3"
    )

    # ---- Build subtriangles (parent/child) ----
    if (!is.na(b$parent_id)) {
      all_subtriangles[[length(all_subtriangles) + 1]] <- data.table(
        chrom       = chr,
        interval_id = b$parent_id,
        sub_id      = bid,
        start_mb    = b$start_mb,
        end_mb      = b$end_mb,
        n_bins      = b$width,
        sub_type    = "child_block"
      )
    }
  }

  message("  Intervals: ", nrow(blocks), " | Subtriangles: ",
          sum(!is.na(blocks$parent_id)))
}

# ---- Write outputs in C01c format ----
iv_dt <- if (length(all_intervals) > 0) rbindlist(all_intervals, fill = TRUE) else data.table()
st_dt <- if (length(all_subtriangles) > 0) rbindlist(all_subtriangles, fill = TRUE) else data.table()

fwrite(iv_dt, file.path(outdir, "triangle_intervals.tsv.gz"), sep = "\t")
fwrite(st_dt, file.path(outdir, "triangle_subtriangles.tsv.gz"), sep = "\t")

# Sample composition: if precomp available, compute basic k=3 PC1 band composition
# This fills the fields C01d needs for scoring dimensions D3-D5
comp_dt <- data.table()
if (!is.null(precomp_dir)) {
  for (rds_file in list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE)) {
    pc <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(pc)) next
    chr <- pc$chrom; dt <- pc$dt
    pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
    if (length(pc1_cols) < 20) next

    chr_blocks <- iv_dt[chrom == chr]
    if (nrow(chr_blocks) == 0) next

    for (r in seq_len(nrow(chr_blocks))) {
      bi <- chr_blocks$L_bin[r]; be <- chr_blocks$R_bin[r]
      if (is.na(bi) || is.na(be) || be > nrow(dt)) next

      # Mean PC1 per sample across block windows
      sub <- as.matrix(dt[bi:be, ..pc1_cols])
      avg_pc1 <- colMeans(sub, na.rm = TRUE)
      avg_pc1 <- avg_pc1[is.finite(avg_pc1)]
      if (length(avg_pc1) < 20) next

      # k=3 clustering
      km <- tryCatch(kmeans(avg_pc1, centers = 3, nstart = 5), error = function(e) NULL)
      if (is.null(km)) next
      co <- order(km$centers[,1])
      band_sizes <- table(km$cluster)

      comp_dt <- rbind(comp_dt, data.table(
        chrom = chr, interval_id = chr_blocks$interval_id[r],
        n_band1 = as.integer(band_sizes[co[1]]),
        n_band2 = as.integer(band_sizes[co[2]]),
        n_band3 = as.integer(band_sizes[co[3]]),
        eff_K = length(unique(km$cluster))
      ), fill = TRUE)
    }
  }
}
fwrite(comp_dt, file.path(outdir, "triangle_sample_composition.tsv.gz"), sep = "\t")
fwrite(data.table(), file.path(outdir, "triangle_bridges.tsv.gz"), sep = "\t")
fwrite(data.table(), file.path(outdir, "triangle_offdiag_linkage.tsv.gz"), sep = "\t")
fwrite(data.table(), file.path(outdir, "triangle_subregimes.tsv.gz"), sep = "\t")

# Comparison table from consensus
if (nrow(cons) > 0) {
  fwrite(cons, file.path(outdir, "triangle_interval_comparison.tsv.gz"), sep = "\t")
} else {
  fwrite(data.table(), file.path(outdir, "triangle_interval_comparison.tsv.gz"), sep = "\t")
}

message("\n=== BRIDGE SUMMARY ===")
message("  Intervals: ", nrow(iv_dt), " across ", length(unique(iv_dt$chrom)), " chromosomes")
message("  Subtriangles (parent/child): ", nrow(st_dt))
message("  Output: ", outdir)
message("")
message("  To use with C01d:")
message("    Rscript STEP_C01d_candidate_scoring.R ", outdir, " <scoring_outdir>")
message("")
message("  This REPLACES the old C01c → C01d flow.")
message("  Everything downstream (C01e-C01m) works unchanged.")
message("\n[DONE]")
