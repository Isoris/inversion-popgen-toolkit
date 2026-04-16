#!/usr/bin/env Rscript
# inspect_precomp_schema.R
#
# Shows what's in precomp RDS files: schema version, NN sim_mat scales,
# matrix dimensions, and summarizes triangle results if available.
#
# Usage:
#   Rscript inspect_precomp_schema.R <precomp_dir> [triangle_multiscale_dir]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript inspect_precomp_schema.R <precomp_dir> [triangle_dir]")

precomp_dir <- args[1]
triangle_dir <- if (length(args) >= 2) args[2] else NULL

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

cat("================================================================\n")
cat("  Precomp Schema Inspection\n")
cat("  Dir:", precomp_dir, "\n")
cat("  Files:", length(rds_files), "\n")
cat("================================================================\n\n")

summary_rows <- list()

for (f in rds_files) {
  chrom <- sub("\\.precomp\\.rds$", "", basename(f))
  pc <- readRDS(f)
  
  fields <- names(pc)
  nn_fields <- grep("^sim_mat_nn", fields, value = TRUE)
  nn_scales <- if (length(nn_fields) > 0) {
    as.integer(sub("^sim_mat_nn", "", nn_fields))
  } else integer(0)
  
  stored_scales <- pc$nn_sim_scales
  
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    chrom = chrom,
    schema = pc$schema_version %||% "unknown",
    n_windows = pc$n_windows,
    n_fields = length(fields),
    sim_mat = !is.null(pc$sim_mat),
    sim_mat_dim = if (!is.null(pc$sim_mat)) paste0(nrow(pc$sim_mat), "x", ncol(pc$sim_mat)) else "NA",
    has_dmat = !is.null(pc$dmat),
    has_mds = !is.null(pc$mds_mat),
    nn_stored = paste(sort(nn_scales), collapse = ","),
    nn_config = paste(sort(stored_scales), collapse = ","),
    n_mds_cols = sum(grepl("^MDS[0-9]+$", names(pc$dt))),
    rds_mb = round(file.info(f)$size / 1e6, 1)
  )
}

summ <- rbindlist(summary_rows)
cat("── Precomp Summary ──\n\n")
print(summ, nrows = 50)

# Check if NN sim_mats are available
nn_present <- unique(summ$nn_stored)
cat("\n── NN Sim_mat Status ──\n")
if (all(nn_present == "")) {
  cat("  NO NN sim_mats stored (schema < v8.5.3)\n")
  cat("  Rerun precomp with NN_SIM_SCALES env var to pre-store:\n")
  cat("    export NN_SIM_SCALES=\"20,40,80\"\n")
  cat("    Rscript STEP_C01a_snake1_precompute.R ...\n")
} else {
  cat("  NN scales stored:", unique(nn_present), "\n")
  cat("  Ready for multiscale triangle detection.\n")
}

# Triangle results summary
if (!is.null(triangle_dir) && dir.exists(triangle_dir)) {
  cat("\n── Triangle Multiscale Results ──\n")
  cat("  Dir:", triangle_dir, "\n\n")
  
  scale_dirs <- list.dirs(triangle_dir, recursive = FALSE, full.names = TRUE)
  scale_dirs <- scale_dirs[grepl("nn[0-9]+$", basename(scale_dirs))]
  
  if (length(scale_dirs) == 0) {
    cat("  No nn* subdirectories found.\n")
  } else {
    tri_summary <- list()
    for (sd in sort(scale_dirs)) {
      scale_tag <- basename(sd)
      iv_file <- file.path(sd, "triangle_intervals.tsv.gz")
      if (!file.exists(iv_file)) next
      
      iv <- fread(iv_file)
      n_total <- nrow(iv)
      n_strong <- sum(iv$interval_type == "strong_triangle", na.rm = TRUE)
      n_moderate <- sum(iv$interval_type == "moderate_triangle", na.rm = TRUE)
      n_patchy <- sum(iv$interval_type == "patchy_signal", na.rm = TRUE)
      sq_med <- round(median(as.numeric(iv$squareness), na.rm = TRUE), 3)
      sq_q90 <- round(quantile(as.numeric(iv$squareness), 0.90, na.rm = TRUE), 3)
      
      tri_summary[[length(tri_summary) + 1]] <- data.table(
        scale = scale_tag,
        n_intervals = n_total,
        strong = n_strong,
        moderate = n_moderate,
        patchy = n_patchy,
        sq_median = sq_med,
        sq_q90 = sq_q90
      )
    }
    
    if (length(tri_summary) > 0) {
      tri_dt <- rbindlist(tri_summary)
      print(tri_dt)
      cat("\n  Scale comparison: fewer intervals + higher squareness at higher k = real structure\n")
    }
  }
}

cat("\n================================================================\n")
cat("  Done.\n")
cat("================================================================\n")
