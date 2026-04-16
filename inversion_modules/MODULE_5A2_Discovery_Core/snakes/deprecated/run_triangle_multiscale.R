#!/usr/bin/env Rscript

# =============================================================================
# run_triangle_multiscale.R
#
# WRAPPER: Runs STEP_C01c_triangle_regimes.R at multiple NN smoothing scales.
#
# For each k in {0, 20, 40, 80}:
#   1. Compute NN-smoothed sim_mat from precomp dmat
#   2. Write temporary precomp RDS with smoothed sim_mat
#   3. Run triangle detection
#   4. Prefix all output files with scale tag for gallery sorting
#
# NN smoothing at scale k:
#   For each window i, replace its MDS coordinates with the average of its
#   k nearest neighbors (by dmat distance). Then recompute the distance matrix
#   and similarity matrix from the smoothed coordinates.
#   k=0 is the raw (no smoothing). Higher k = coarser, picks up larger structures.
#
# Usage:
#   Rscript run_triangle_multiscale.R <precomp_dir> <outdir> [--chrom C_gar_LG01] [--scales 0,20,40,80]
#
# Output structure:
#   <outdir>/nn0/    — raw sim_mat results (same as running C01c directly)
#   <outdir>/nn20/   — NN k=20 smoothed
#   <outdir>/nn40/   — NN k=40 smoothed
#   <outdir>/nn80/   — NN k=80 smoothed
#
# Each subdirectory has the full triangle output set.
# Plot files are prefixed: C_gar_LG01_nn0_triangles.png, C_gar_LG01_nn20_triangles.png, etc.
# for easy left-right gallery scrolling.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_triangle_multiscale.R <precomp_dir> <outdir> [opts]")

precomp_dir <- args[1]
base_outdir <- args[2]
chrom_filter <- NULL
scales_str <- Sys.getenv("NN_SIM_SCALES", "20,40,80")
# Always include 0 (raw) as first scale
scales_str_full <- paste0("0,", scales_str)
triangle_script <- NULL

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args))  { chrom_filter <- args[i+1]; i <- i+2L }
  else if (a == "--scales" && i < length(args)) { scales_str_full <- args[i+1]; i <- i+2L }
  else if (a == "--script" && i < length(args)) { triangle_script <- args[i+1]; i <- i+2L }
  else { i <- i+1L }
}

nn_scales <- as.integer(strsplit(scales_str_full, ",")[[1]])
message("[multiscale] Scales: ", paste(nn_scales, collapse = ", "))

# Find triangle script
if (is.null(triangle_script)) {
  # Try relative to this script
  this_script <- commandArgs(trailingOnly = FALSE)
  this_script <- this_script[grep("--file=", this_script)]
  if (length(this_script) > 0) {
    this_dir <- dirname(sub("^--file=", "", this_script[1]))
    triangle_script <- file.path(this_dir, "STEP_C01c_triangle_regimes.R")
  }
}
if (is.null(triangle_script) || !file.exists(triangle_script)) {
  # Try from env
  codebase <- Sys.getenv("DISCOVERY2DIR", "")
  if (nzchar(codebase)) {
    triangle_script <- file.path(codebase, "snakes", "STEP_C01c_triangle_regimes.R")
  }
}
if (is.null(triangle_script) || !file.exists(triangle_script)) {
  stop("[multiscale] Cannot find STEP_C01c_triangle_regimes.R. Use --script <path>")
}
message("[multiscale] Triangle script: ", triangle_script)

RSCRIPT_BIN <- Sys.getenv("RSCRIPT_BIN", "Rscript")

# =============================================================================
# NN SMOOTHING FUNCTION
# =============================================================================

nn_smooth_simmat <- function(pc, k) {
  # Smooth MDS coordinates by averaging k nearest neighbors,
  # then recompute distance matrix and similarity matrix.
  #
  # k=0: return raw sim_mat unchanged
  # k>0: for each window, average MDS coords of its k nearest neighbors
  #
  # Works with old precomp (no dmat) by computing distances from mds_mat.
  
  if (k == 0) return(pc$sim_mat)
  
  # Get MDS coordinates
  mds_cols <- grep("^MDS[0-9]+$", names(pc$dt), value = TRUE)
  if (length(mds_cols) == 0) {
    message("  [nn] No MDS columns found, using raw sim_mat")
    return(pc$sim_mat)
  }
  mds_mat <- as.matrix(pc$dt[, ..mds_cols])
  n <- nrow(mds_mat)
  k_use <- min(k, n - 1L)
  
  # Get or compute distance matrix for NN lookup
  if (!is.null(pc$dmat) && nrow(pc$dmat) == n) {
    dmat <- pc$dmat
    message("  [nn] Using stored dmat for NN lookup")
  } else {
    # Compute from MDS coordinates (Euclidean distance in MDS space)
    message("  [nn] Computing dmat from mds_mat (", n, " x ", ncol(mds_mat), " dims)...")
    dmat <- as.matrix(dist(mds_mat))
  }
  
  # For each window, find k nearest neighbors and average their MDS coords
  smoothed <- matrix(0, nrow = n, ncol = ncol(mds_mat))
  for (i in seq_len(n)) {
    d <- dmat[i, ]
    d[i] <- Inf
    nn_idx <- order(d)[seq_len(k_use)]
    # Include self in average
    all_idx <- c(i, nn_idx)
    smoothed[i, ] <- colMeans(mds_mat[all_idx, , drop = FALSE], na.rm = TRUE)
  }
  
  # Recompute distance matrix from smoothed coordinates
  new_dmat <- as.matrix(dist(smoothed))
  
  # Convert to similarity (same formula as make_sim_mat in precomp)
  finite_vals <- new_dmat[is.finite(new_dmat)]
  dmax <- if (length(finite_vals) > 0) quantile(finite_vals, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(new_dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  
  sim
}

# =============================================================================
# LOAD PRECOMP FILES
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

chroms <- sub("\\.precomp\\.rds$", "", basename(rds_files))
if (!is.null(chrom_filter)) {
  keep <- chroms %in% chrom_filter
  rds_files <- rds_files[keep]
  chroms <- chroms[keep]
}
message("[multiscale] ", length(chroms), " chromosomes")

# =============================================================================
# RUN EACH SCALE
# =============================================================================

for (k in nn_scales) {
  scale_tag <- paste0("nn", k)
  scale_outdir <- file.path(base_outdir, scale_tag)
  temp_precomp <- file.path(base_outdir, paste0(".tmp_precomp_", scale_tag))
  dir.create(scale_outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_precomp, recursive = TRUE, showWarnings = FALSE)
  
  message("\n================================================================")
  message("[multiscale] === Scale: ", scale_tag, " (k=", k, ") ===")
  message("================================================================")
  
  # Create temporary precomp RDS files with the right sim_mat for this scale
  for (fi in seq_along(rds_files)) {
    chr <- chroms[fi]
    message("[multiscale] ", chr, " k=", k, ": loading...")
    pc <- readRDS(rds_files[fi])
    
    # Try to load sim_mat from separate file (v8.5.3 format)
    sim_mat_rds <- file.path(precomp_dir, "sim_mats",
                              paste0(chr, ".sim_mat_nn", k, ".rds"))
    
    if (file.exists(sim_mat_rds)) {
      # v8.5.3: separate sim_mat file — fast, just load it
      pc$sim_mat <- readRDS(sim_mat_rds)
      message("[multiscale]   Loaded ", basename(sim_mat_rds),
              " (", nrow(pc$sim_mat), "x", ncol(pc$sim_mat), ")")
    } else if (k == 0 && !is.null(pc$sim_mat)) {
      # Old precomp with embedded sim_mat — use as-is for nn0
      message("[multiscale]   Using embedded sim_mat (old precomp)")
    } else if (k > 0) {
      # Check if NN sim_mat is embedded in precomp (intermediate format)
      nn_field <- paste0("sim_mat_nn", k)
      if (!is.null(pc[[nn_field]])) {
        pc$sim_mat <- pc[[nn_field]]
        message("[multiscale]   Using embedded ", nn_field)
      } else {
        # Last resort: compute on the fly
        message("[multiscale]   Computing on the fly (k=", k, ")...")
        message("[multiscale]   (Hint: rerun precomp v8.5.3 to pre-store)")
        t0 <- proc.time()
        pc$sim_mat <- nn_smooth_simmat(pc, k)
        elapsed <- round((proc.time() - t0)[3], 1)
        message("[multiscale]   Computed in ", elapsed, "s")
      }
    } else {
      stop("[multiscale] No sim_mat available for ", chr, " k=", k)
    }
    
    # Save temporary precomp with the selected sim_mat
    tmp_rds <- file.path(temp_precomp, paste0(chr, ".precomp.rds"))
    saveRDS(pc, tmp_rds)
  }
  
  # Run triangle detection
  chrom_opts <- if (!is.null(chrom_filter)) paste("--chrom", chrom_filter) else ""
  cmd <- paste(RSCRIPT_BIN, shQuote(triangle_script),
               shQuote(temp_precomp), shQuote(scale_outdir),
               chrom_opts)
  message("[multiscale] Running: ", cmd)
  ret <- system(cmd)
  if (ret != 0) message("[multiscale] WARNING: triangle script returned ", ret)
  
  # Rename plot files with scale prefix for gallery sorting
  plot_dir <- file.path(scale_outdir, "plots")
  if (dir.exists(plot_dir)) {
    plot_files <- list.files(plot_dir, pattern = "\\.png$", full.names = FALSE)
    for (pf in plot_files) {
      # C_gar_LG01_triangles.png → C_gar_LG01_nn0_triangles.png
      # Insert scale tag after chromosome name
      parts <- strsplit(pf, "_")[[1]]
      # Find where the chr name ends (after LG##)
      lg_idx <- grep("^LG[0-9]+$", parts)
      if (length(lg_idx) > 0) {
        new_name <- paste0(
          paste(parts[1:max(lg_idx)], collapse = "_"), "_",
          scale_tag, "_",
          paste(parts[(max(lg_idx)+1):length(parts)], collapse = "_"))
        file.rename(file.path(plot_dir, pf), file.path(plot_dir, new_name))
      }
    }
  }
  
  # Clean up temp precomp
  unlink(temp_precomp, recursive = TRUE)
  
  message("[multiscale] ", scale_tag, " DONE -> ", scale_outdir)
}

message("\n================================================================")
message("[multiscale] All scales complete.")
message("  Output: ", base_outdir, "/")
message("  Scales: ", paste(paste0("nn", nn_scales), collapse = ", "))
message("  Gallery order: sort by filename prefix (nn0, nn20, nn40, nn80)")
message("================================================================")
