#!/usr/bin/env Rscript

# =============================================================================
# STEP20c_gradient_marker_experiments.R
#
# Per-marker association with anchor-derived continuous gradients.
#
# For each candidate:
#   1. Load anchor projection from STEP20b (u, center_distance, branch_balance)
#   2. Load regional dosage
#   3. For each marker: Spearman correlation with each gradient
#   4. Two modes:
#      MODE A: all samples
#      MODE B: HET-only samples
#   5. Classify markers by dominant gradient association
#
# Marker class hints:
#   - core_inversion_axis: strongest on broad_axis_position (u)
#   - het_center_distance: strongest on center_distance, especially in HET-only
#   - branch_specific: strongest on branch_balance, especially in HET-only
#   - uninformative: weak or inconsistent
#
# Usage:
#   Rscript STEP20c_gradient_marker_experiments.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

has_ggplot <- suppressWarnings(require(ggplot2, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
message("[INFO] Processing ", nrow(cand), " candidates")

# ── Helper: safe Spearman correlation ────────────────────────────────────────
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 5) return(list(rho = NA_real_, p = NA_real_, n = n))
  test <- tryCatch(
    cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) return(list(rho = NA_real_, p = NA_real_, n = n))
  list(rho = test$estimate, p = test$p.value, n = n)
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  message("\n[INFO] Candidate ", cid, " — gradient-marker experiments")

  # ── Load anchor projection ───────────────────────────────────────────────
  anchor_file <- file.path(cand_dir, "candidate_anchor_projection.tsv")
  if (!file.exists(anchor_file)) { message("[SKIP] No anchor file"); next }
  anchor <- fread(anchor_file)
  if (nrow(anchor) < 10) { message("[SKIP] Too few samples in anchor"); next }

  # ── Load dosage ──────────────────────────────────────────────────────────
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[SKIP] Missing dosage"); next
  }

  dos <- fread(dos_file)
  sites <- fread(sites_file)
  dos_sample_cols <- setdiff(names(dos), "marker")

  # Map Ind-style if needed
  if (all(grepl("^Ind", dos_sample_cols)) && length(dos_sample_cols) == nrow(anchor)) {
    setnames(dos, old = dos_sample_cols, new = anchor$sample)
    dos_sample_cols <- anchor$sample
  }

  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) { message("[SKIP] Too few markers"); next }

  sites_reg <- sites[keep]
  X <- as.matrix(dos[keep, ..dos_sample_cols])
  storage.mode(X) <- "double"
  n_markers <- nrow(X)

  # Align anchor samples with dosage columns
  sample_match <- match(dos_sample_cols, anchor$sample)
  if (any(is.na(sample_match))) {
    # Try to find common subset
    common <- intersect(dos_sample_cols, anchor$sample)
    if (length(common) < 10) { message("[SKIP] Sample mismatch"); next }
    X <- X[, match(common, dos_sample_cols), drop = FALSE]
    dos_sample_cols <- common
    sample_match <- match(common, anchor$sample)
  }

  u_vals <- anchor$u[sample_match]
  cd_vals <- anchor$center_distance[sample_match]
  bb_vals <- anchor$branch_balance[sample_match]
  coarse_group <- anchor$coarse_group[sample_match]

  # ── Identify HET samples ────────────────────────────────────────────────
  het_labels <- c("Het", "HET", "G2")
  is_het <- coarse_group %in% het_labels
  n_het <- sum(is_het)

  message("[INFO] Markers: ", n_markers, "  Samples: ", ncol(X),
          "  HET: ", n_het)

  # ═════════════════════════════════════════════════════════════════════════
  # MODE A: ALL SAMPLES
  # ═════════════════════════════════════════════════════════════════════════

  message("[INFO] Mode A: all samples")
  results_all <- vector("list", n_markers)

  for (mi in seq_len(n_markers)) {
    dosage_vec <- X[mi, ]

    sp_u  <- safe_spearman(dosage_vec, u_vals)
    sp_cd <- safe_spearman(dosage_vec, cd_vals)
    sp_bb <- safe_spearman(dosage_vec, bb_vals)

    # Dominant gradient
    rhos <- c(u = abs(sp_u$rho), cd = abs(sp_cd$rho), bb = abs(sp_bb$rho))
    rhos[is.na(rhos)] <- 0
    dominant <- names(which.max(rhos))

    # Marker class hint
    marker_class <- "uninformative"
    if (rhos[dominant] >= 0.3) {
      marker_class <- switch(dominant,
        u = "core_inversion_axis",
        cd = "het_center_distance",
        bb = "branch_specific",
        "uninformative"
      )
    }

    results_all[[mi]] <- data.table(
      candidate_id = cid,
      marker = sites_reg$marker[mi],
      chrom = chr,
      pos = sites_reg$pos[mi],
      n_samples = sp_u$n,
      corr_u = round(sp_u$rho, 6),
      p_u = sp_u$p,
      corr_center_distance = round(sp_cd$rho, 6),
      p_center_distance = sp_cd$p,
      corr_branch_balance = round(sp_bb$rho, 6),
      p_branch_balance = sp_bb$p,
      dominant_gradient = dominant,
      marker_class_hint = marker_class
    )
  }

  all_dt <- rbindlist(results_all)
  fwrite(all_dt,
         file.path(cand_dir, "candidate_marker_gradient_all.tsv"), sep = "\t")

  # ═════════════════════════════════════════════════════════════════════════
  # MODE B: HET-ONLY SAMPLES
  # ═════════════════════════════════════════════════════════════════════════

  if (n_het >= 10) {
    message("[INFO] Mode B: HET-only (", n_het, " samples)")
    het_idx <- which(is_het)
    results_het <- vector("list", n_markers)

    for (mi in seq_len(n_markers)) {
      dosage_het <- X[mi, het_idx]
      u_het <- u_vals[het_idx]
      cd_het <- cd_vals[het_idx]
      bb_het <- bb_vals[het_idx]

      sp_u  <- safe_spearman(dosage_het, u_het)
      sp_cd <- safe_spearman(dosage_het, cd_het)
      sp_bb <- safe_spearman(dosage_het, bb_het)

      rhos <- c(u = abs(sp_u$rho), cd = abs(sp_cd$rho), bb = abs(sp_bb$rho))
      rhos[is.na(rhos)] <- 0
      dominant <- names(which.max(rhos))

      marker_class_het <- "uninformative"
      if (rhos[dominant] >= 0.25) {
        marker_class_het <- switch(dominant,
          u = "het_axis_residual",
          cd = "het_center_distance",
          bb = "het_branch_specific",
          "uninformative"
        )
      }

      results_het[[mi]] <- data.table(
        candidate_id = cid,
        marker = sites_reg$marker[mi],
        chrom = chr,
        pos = sites_reg$pos[mi],
        n_het_samples = sp_u$n,
        corr_u_het = round(sp_u$rho, 6),
        p_u_het = sp_u$p,
        corr_center_distance_het = round(sp_cd$rho, 6),
        p_center_distance_het = sp_cd$p,
        corr_branch_balance_het = round(sp_bb$rho, 6),
        p_branch_balance_het = sp_bb$p,
        dominant_gradient_het = dominant,
        marker_class_hint_het = marker_class_het
      )
    }

    het_dt <- rbindlist(results_het)
    fwrite(het_dt,
           file.path(cand_dir, "candidate_marker_gradient_het_only.tsv"), sep = "\t")
  } else {
    message("[INFO] Too few HET samples (", n_het, "); skipping MODE B")
  }

  # ═════════════════════════════════════════════════════════════════════════
  # TOP GRADIENT MARKERS
  # ═════════════════════════════════════════════════════════════════════════

  # Rank markers by maximum absolute correlation across all modes
  top_markers <- copy(all_dt)
  top_markers[, max_abs_rho := pmax(abs(corr_u), abs(corr_center_distance),
                                     abs(corr_branch_balance), na.rm = TRUE)]
  top_markers <- top_markers[order(-max_abs_rho)]

  # Add HET-only info if available
  if (exists("het_dt") && nrow(het_dt) > 0) {
    top_markers <- merge(top_markers,
                          het_dt[, .(marker, corr_u_het, corr_center_distance_het,
                                     corr_branch_balance_het, marker_class_hint_het)],
                          by = "marker", all.x = TRUE)
  }

  top_markers[, score_rank := seq_len(.N)]

  # Keep top 100 or all if fewer
  n_top <- min(100, nrow(top_markers))
  top_markers_out <- top_markers[seq_len(n_top)]
  top_markers_out[, candidate_id := cid]

  fwrite(top_markers_out,
         file.path(cand_dir, "candidate_top_gradient_markers.tsv"), sep = "\t")

  # ── Diagnostic plots for top markers ─────────────────────────────────────
  if (has_ggplot && n_top >= 5) {
    top5 <- top_markers_out[seq_len(min(5, n_top))]

    gradient_plot_dir <- file.path(plot_dir, "gradient_marker_plots")
    ensure_dir(file.path(gradient_plot_dir, "all_samples"))
    if (n_het >= 10) ensure_dir(file.path(gradient_plot_dir, "het_only"))

    for (ri in seq_len(nrow(top5))) {
      mk <- top5$marker[ri]
      mk_idx <- which(sites_reg$marker == mk)
      if (length(mk_idx) != 1) next

      dos_vec <- X[mk_idx, ]

      plot_dt <- data.table(
        sample = dos_sample_cols,
        dosage = dos_vec,
        u = u_vals,
        center_distance = cd_vals,
        branch_balance = bb_vals,
        group = coarse_group
      )

      # All-sample: dosage vs u
      p <- ggplot(plot_dt, aes(x = u, y = dosage, color = group)) +
        geom_point(size = 1.8, alpha = 0.7) +
        geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.6) +
        scale_color_manual(values = GROUP_COLORS, na.value = "grey60") +
        theme_bw(base_size = 10) +
        labs(title = paste0("Marker ", mk, " vs u (rank ", ri, ")"),
             subtitle = paste0("rho=", top5$corr_u[ri], " class=", top5$marker_class_hint[ri]),
             x = "Broad axis position (u)", y = "Dosage")
      ggsave(file.path(gradient_plot_dir, "all_samples",
                        paste0("rank", ri, "_", gsub(":", "_", mk), "_vs_u.png")),
             p, width = 6, height = 4, dpi = 200)

      # HET-only: dosage vs branch_balance
      if (n_het >= 10) {
        het_plot_dt <- plot_dt[group %in% het_labels]
        if (nrow(het_plot_dt) >= 5) {
          p_het <- ggplot(het_plot_dt, aes(x = branch_balance, y = dosage)) +
            geom_point(size = 2, alpha = 0.7, color = "#D6604D") +
            geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.6) +
            theme_bw(base_size = 10) +
            labs(title = paste0("HET-only: ", mk, " vs branch_balance"),
                 x = "Branch balance", y = "Dosage")
          ggsave(file.path(gradient_plot_dir, "het_only",
                            paste0("rank", ri, "_", gsub(":", "_", mk), "_het_bb.png")),
                 p_het, width = 6, height = 4, dpi = 200)
        }
      }
    }
  }

  message("[INFO] Candidate ", cid, " — gradient markers complete")
}

message("\n[DONE] STEP20c gradient-marker experiments complete")
