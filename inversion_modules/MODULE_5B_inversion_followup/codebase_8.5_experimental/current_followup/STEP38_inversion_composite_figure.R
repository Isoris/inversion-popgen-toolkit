#!/usr/bin/env Rscript

# =============================================================================
# STEP38_inversion_composite_figure.R  (v1.0)
#
# COMPOSITE INVERSION FIGURE — one publication-quality PDF per candidate.
#
# Combines breakpoint SV evidence with population-level inversion signal
# into a single multi-panel figure.
#
# PANELS (top to bottom):
#
#   A — SV BREAKPOINT RIBBON
#       Horizontal ribbon showing the candidate region with left/right
#       boundary zones. SV calls (INV + BND) plotted as jittered points
#       colored by caller (DELLY/Manta) and shaped by SV type. Breakpoint
#       clusters shown as shaded zones. Per-group support bars below.
#
#   B — REGIONAL PCA (3-BAND SCATTER)
#       Classic inversion PCA: PC1 vs PC2 of the candidate region,
#       colored by inferred group (HOMO_1 / HET / HOMO_2 / AMBIGUOUS).
#       Sample counts annotated per group.
#
#   C — THETA (tP) DIVERSITY TRACK
#       Per-group population theta (pairwise θ_π) across the candidate
#       + flanking region. Shows elevated diversity in heterokaryotypes.
#
#   D — SNAKE SUPPORT TRACKS
#       3 horizontal tracks showing per-window PASS/WEAK/FAIL status
#       for Snake 1 (PCA), Snake 2 (community), Snake 3 (GHSL).
#       Consensus support signature below.
#
#   E — BREAKPOINT STATISTICS
#       Summary panel: Fisher test enrichment per cluster, concordance
#       counts, group × boundary support heatmap.
#
# INPUTS:
#   <config.R>     — standard followup config
#   <cid>          — candidate ID
#   Uses existing outputs from:
#     STEP12 (regional PCA), STEP10b (theta), STEP10e/g/h/i (snakes),
#     STEP37 (breakpoint support)
#
# OUTPUTS:
#   candidate_composite_figure.pdf — multi-panel figure
#
# Usage:
#   Rscript STEP38_inversion_composite_figure.R <config.R> <candidate_id>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP38_inversion_composite_figure.R <config.R> <candidate_id>")
}

config_file <- args[1]
cid         <- as.integer(args[2])

source(config_file)

cand <- fread(CANDIDATE_TABLE)
row  <- cand[candidate_id == cid]
if (nrow(row) == 0) stop("Candidate ", cid, " not found")

chr      <- row$chrom
c_start  <- as.numeric(row$start_bp)
c_end    <- as.numeric(row$end_bp)
span_mb  <- (c_end - c_start) / 1e6
flank_bp <- max(100000, (c_end - c_start) * 0.3)

cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))
cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
plots_dir <- file.path(PLOTS_DIR, cand_prefix)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

message("[STEP38] Composite figure for candidate ", cid,
        " (", chr, ":", round(c_start/1e6, 2), "–", round(c_end/1e6, 2), " Mb)")

# =============================================================================
# HELPER: safe file loader
# =============================================================================

safe_load <- function(path) {
  if (file.exists(path)) tryCatch(fread(path), error = function(e) NULL) else NULL
}

fmt_mb <- function(bp) sprintf("%.2f", bp / 1e6)

# =============================================================================
# LOAD ALL DATA SOURCES
# =============================================================================

# STEP12: Regional PCA
pca_dt <- safe_load(file.path(REGIONAL_DIR, "inversion_localpca",
                              paste0("inversion_localpca.candidate_", cid,
                                     ".regional_pca_samples.tsv.gz")))

# STEP21: Group assignments
group_dt <- safe_load(file.path(cand_dir, "candidate_pca_rotated.tsv"))
if (is.null(group_dt)) {
  group_dt <- safe_load(file.path(cand_dir, "candidate_group_assignments.tsv"))
}

# STEP10b: Theta windows
theta_dt <- safe_load(file.path(THETA_OVERLAP_DIR, "het_bridge.sample_tP_windows.tsv.gz"))

# STEP37: Breakpoint support
bp_calls    <- safe_load(file.path(cand_dir, "candidate_bp_sv_calls.tsv.gz"))
bp_clusters <- safe_load(file.path(cand_dir, "candidate_bp_clusters.tsv"))
bp_tests    <- safe_load(file.path(cand_dir, "candidate_bp_group_tests.tsv"))
bp_conc     <- safe_load(file.path(cand_dir, "candidate_bp_concordance.tsv"))
bp_summary  <- safe_load(file.path(cand_dir, "candidate_bp_summary.tsv"))

# Snake tracks
snake1_dir <- file.path(INV_ROOT, "06_mds_candidates", "snake_regions")
snake2_dir <- file.path(INV_ROOT, "06_mds_candidates", "snake2_community")
snake3_dir <- file.path(INV_ROOT, "06_mds_candidates", "snake3_ghsl")
consensus_dir <- file.path(INV_ROOT, "06_mds_candidates", "three_snake_consensus")

s1_track <- safe_load(file.path(snake1_dir, "snake_windows.tsv.gz"))
s2_track <- safe_load(file.path(snake2_dir, "snake2_track.tsv.gz"))
s3_track <- safe_load(file.path(snake3_dir, "snake3_track.tsv.gz"))
cons_track <- safe_load(file.path(consensus_dir, "consensus_windows.tsv.gz"))

# =============================================================================
# PANEL A — SV BREAKPOINT RIBBON
# =============================================================================

build_panel_A <- function() {
  # Define plot range
  x_min <- c_start - flank_bp
  x_max <- c_end + flank_bp

  # Base ribbon showing candidate region
  p <- ggplot() +
    # Candidate region shading
    annotate("rect", xmin = c_start, xmax = c_end, ymin = -Inf, ymax = Inf,
             fill = "#e8e0f0", alpha = 0.5) +
    # Boundary zones
    annotate("rect", xmin = c_start - 50000, xmax = c_start + 50000,
             ymin = -Inf, ymax = Inf, fill = "#fde68a", alpha = 0.3) +
    annotate("rect", xmin = c_end - 50000, xmax = c_end + 50000,
             ymin = -Inf, ymax = Inf, fill = "#fde68a", alpha = 0.3) +
    # Boundary lines
    geom_vline(xintercept = c(c_start, c_end), linetype = "dashed",
               color = "#7c3aed", linewidth = 0.5) +
    annotate("text", x = c_start, y = 1.05, label = paste0("L: ", fmt_mb(c_start), " Mb"),
             hjust = 0.5, vjust = 0, size = 2.8, color = "#7c3aed") +
    annotate("text", x = c_end, y = 1.05, label = paste0("R: ", fmt_mb(c_end), " Mb"),
             hjust = 0.5, vjust = 0, size = 2.8, color = "#7c3aed")

  if (!is.null(bp_calls) && nrow(bp_calls) > 0) {
    # Filter to plot range
    bp_plot <- bp_calls[pos1 >= x_min & pos1 <= x_max]
    if (nrow(bp_plot) > 0) {
      # Jitter y by caller and type
      bp_plot[, y_jitter := runif(.N, 0.1, 0.9)]
      bp_plot[, svtype_upper := toupper(svtype)]

      p <- p +
        geom_point(data = bp_plot, aes(x = pos1, y = y_jitter,
                                        color = caller, shape = svtype_upper),
                   size = 2, alpha = 0.7) +
        scale_color_manual(values = c("DELLY" = "#2563eb", "Manta" = "#dc2626"),
                           name = "Caller") +
        scale_shape_manual(values = c("INV" = 17, "BND" = 4), name = "SV Type")
    }

    # Cluster zones
    if (!is.null(bp_clusters) && nrow(bp_clusters) > 0) {
      cl_plot <- bp_clusters[cluster_start >= x_min & cluster_end <= x_max]
      if (nrow(cl_plot) > 0) {
        p <- p +
          geom_rect(data = cl_plot,
                    aes(xmin = cluster_start, xmax = cluster_end,
                        ymin = -0.15, ymax = -0.02),
                    fill = "#f59e0b", alpha = 0.6, color = "#92400e", linewidth = 0.3) +
          geom_text(data = cl_plot,
                    aes(x = (cluster_start + cluster_end) / 2, y = -0.22,
                        label = paste0("n=", n_unique_samples)),
                    size = 2.2, color = "#92400e")
      }
    }
  }

  p <- p +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6),
                       limits = c(x_min, x_max)) +
    scale_y_continuous(limits = c(-0.3, 1.15), breaks = NULL) +
    labs(x = NULL, y = NULL,
         title = paste0("A — Structural Variant Breakpoint Evidence"),
         subtitle = paste0(chr, ":", fmt_mb(c_start), "–", fmt_mb(c_end),
                           " Mb (", round(span_mb, 2), " Mb span)")) +
    theme_inversion(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "right", legend.key.size = unit(0.35, "cm"),
          plot.margin = margin(5, 10, 2, 10))

  p
}

# =============================================================================
# PANEL B — REGIONAL PCA (3-BAND SCATTER)
# =============================================================================

build_panel_B <- function() {
  if (is.null(pca_dt) || nrow(pca_dt) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No PCA data") +
             theme_void() + labs(title = "B — Regional PCA"))
  }

  # Find PC columns
  pc1_col <- intersect(names(pca_dt), c("PC1", "rPC1", "regional_PC1"))
  pc2_col <- intersect(names(pca_dt), c("PC2", "rPC2", "regional_PC2"))
  grp_col <- intersect(names(pca_dt), c("group", "cluster", "assignment"))

  if (length(pc1_col) == 0 || length(pc2_col) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Missing PC columns") +
             theme_void() + labs(title = "B — Regional PCA"))
  }

  pca_dt[, pc1 := .SD[[1]], .SDcols = pc1_col[1]]
  pca_dt[, pc2 := .SD[[1]], .SDcols = pc2_col[1]]

  if (length(grp_col) > 0) {
    pca_dt[, grp := as.character(.SD[[1]]), .SDcols = grp_col[1]]
  } else if (!is.null(group_dt) && "group" %in% names(group_dt)) {
    samp_col <- intersect(names(pca_dt), c("sample", "sample_id"))
    if (length(samp_col) > 0) {
      pca_dt <- merge(pca_dt, group_dt[, .(sample, group)],
                       by.x = samp_col[1], by.y = "sample", all.x = TRUE)
      pca_dt[, grp := as.character(group)]
    } else {
      pca_dt[, grp := "UNKNOWN"]
    }
  } else {
    pca_dt[, grp := "UNKNOWN"]
  }

  grp_counts <- pca_dt[, .N, by = grp]
  pca_dt <- merge(pca_dt, grp_counts, by = "grp", suffixes = c("", "_count"))

  colors <- safe_group_colors(pca_dt$grp)
  shapes <- safe_shape_values(pca_dt$grp)

  # Group count labels
  count_label <- paste(grp_counts$grp, "=", grp_counts$N, collapse = "  ")

  p <- ggplot(pca_dt, aes(x = pc1, y = pc2, color = grp, shape = grp)) +
    geom_point(size = 1.8, alpha = 0.75) +
    scale_color_manual(values = colors, name = "Group") +
    scale_shape_manual(values = shapes, name = "Group") +
    labs(x = "PC1", y = "PC2",
         title = "B — Regional PCA (3-band inversion pattern)",
         subtitle = count_label) +
    theme_inversion(base_size = 10) +
    theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
          plot.margin = margin(5, 10, 2, 10))

  p
}

# =============================================================================
# PANEL C — THETA (tP) DIVERSITY TRACK
# =============================================================================

build_panel_C <- function() {
  if (is.null(theta_dt) || nrow(theta_dt) == 0 || is.null(group_dt)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No theta data") +
             theme_void() + labs(title = "C — Population Theta"))
  }

  x_min <- c_start - flank_bp
  x_max <- c_end + flank_bp

  # Theta has windows; filter to our region
  theta_cols <- intersect(names(theta_dt), c("chrom", "WinStart", "WinStop",
                                              "WinCenter", "tP", "sample"))
  if (length(theta_cols) < 4) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Theta format mismatch") +
             theme_void() + labs(title = "C — Population Theta"))
  }

  theta_chr <- theta_dt[chrom == chr & WinCenter >= x_min & WinCenter <= x_max]
  if (nrow(theta_chr) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No theta windows in range") +
             theme_void() + labs(title = "C — Population Theta"))
  }

  # Merge with groups
  samp_col <- intersect(names(group_dt), c("sample", "sample_id"))
  if (length(samp_col) > 0 && "sample" %in% names(theta_chr)) {
    theta_grp <- merge(theta_chr, group_dt[, .SD, .SDcols = c(samp_col[1], "group")],
                        by.x = "sample", by.y = samp_col[1], all.x = TRUE)
    theta_grp[is.na(group), group := "UNKNOWN"]

    # Per-group mean theta per window
    theta_agg <- theta_grp[, .(mean_tP = mean(tP, na.rm = TRUE),
                                 se_tP = sd(tP, na.rm = TRUE) / sqrt(.N)),
                             by = .(WinCenter, group)]

    colors <- safe_group_colors(theta_agg$group)

    p <- ggplot(theta_agg, aes(x = WinCenter, y = mean_tP, color = group)) +
      annotate("rect", xmin = c_start, xmax = c_end, ymin = -Inf, ymax = Inf,
               fill = "#e8e0f0", alpha = 0.3) +
      geom_vline(xintercept = c(c_start, c_end), linetype = "dashed",
                 color = "#7c3aed", linewidth = 0.3) +
      geom_line(linewidth = 0.6, alpha = 0.8) +
      geom_ribbon(aes(ymin = mean_tP - se_tP, ymax = mean_tP + se_tP, fill = group),
                  alpha = 0.15, color = NA) +
      scale_color_manual(values = colors, name = "Group") +
      scale_fill_manual(values = colors, guide = "none") +
      scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) +
      labs(x = NULL, y = expression(theta[pi]),
           title = "C — Pairwise Theta Diversity by Group",
           subtitle = "Elevated θπ in heterokaryotypes expected within inversion region") +
      theme_inversion(base_size = 10) +
      theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
            plot.margin = margin(5, 10, 2, 10))
  } else {
    p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Cannot merge theta + groups") +
      theme_void() + labs(title = "C — Population Theta")
  }

  p
}

# =============================================================================
# PANEL D — SNAKE SUPPORT TRACKS
# =============================================================================

build_panel_D <- function() {
  x_min <- c_start - flank_bp
  x_max <- c_end + flank_bp

  # Build a long-form table: window_center, track, status
  tracks <- list()

  if (!is.null(s1_track) && nrow(s1_track) > 0) {
    s1_chr <- s1_track[chrom == chr]
    if (nrow(s1_chr) > 0) {
      s1_chr[, center := (start_bp + end_bp) / 2]
      s1_chr[, status := "PASS"]  # collected windows are all PASS for snake 1
      tracks[["Snake 1\n(PCA)"]] <- s1_chr[center >= x_min & center <= x_max,
                                             .(center, status)]
    }
  }

  if (!is.null(s2_track) && nrow(s2_track) > 0) {
    s2_chr <- s2_track[chrom == chr]
    if (nrow(s2_chr) > 0) {
      s2_chr[, center := (start_bp + end_bp) / 2]
      s2_chr[, status := snake2_status]
      tracks[["Snake 2\n(Community)"]] <- s2_chr[center >= x_min & center <= x_max,
                                                   .(center, status)]
    }
  }

  if (!is.null(s3_track) && nrow(s3_track) > 0) {
    s3_chr <- s3_track[chrom == chr]
    if (nrow(s3_chr) > 0) {
      s3_chr[, center := (start_bp + end_bp) / 2]
      s3_chr[, status := snake3_status]
      tracks[["Snake 3\n(GHSL)"]] <- s3_chr[center >= x_min & center <= x_max,
                                              .(center, status)]
    }
  }

  if (length(tracks) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No snake tracks") +
             theme_void() + labs(title = "D — Snake Support"))
  }

  # Combine
  track_dt <- rbindlist(lapply(names(tracks), function(tn) {
    d <- tracks[[tn]]
    d[, track := tn]
    d
  }))

  status_colors <- c("PASS" = "#16a34a", "WEAK" = "#f59e0b", "FAIL" = "#e5e7eb")

  p <- ggplot(track_dt, aes(x = center, y = track, fill = status)) +
    annotate("rect", xmin = c_start, xmax = c_end,
             ymin = -Inf, ymax = Inf, fill = "#e8e0f0", alpha = 0.3) +
    geom_vline(xintercept = c(c_start, c_end), linetype = "dashed",
               color = "#7c3aed", linewidth = 0.3) +
    geom_tile(height = 0.7, width = (x_max - x_min) / max(50, nrow(track_dt) / 3)) +
    scale_fill_manual(values = status_colors, name = "Status") +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6)) +
    labs(x = paste0("Position on ", chr, " (Mb)"), y = NULL,
         title = "D — Snake Support Tracks",
         subtitle = "Per-window PASS/WEAK/FAIL across 3 continuity snakes") +
    theme_inversion(base_size = 10) +
    theme(legend.position = "right", legend.key.size = unit(0.35, "cm"),
          axis.text.y = element_text(size = 8),
          plot.margin = margin(5, 10, 2, 10))

  # Add consensus if available
  if (!is.null(cons_track) && nrow(cons_track) > 0) {
    cons_chr <- cons_track[chrom == chr]
    if (nrow(cons_chr) > 0) {
      cons_chr[, center := (start_bp + end_bp) / 2]
      cons_chr <- cons_chr[center >= x_min & center <= x_max]
      if (nrow(cons_chr) > 0) {
        sig_colors <- c("PCA+GROUP+GHSL" = "#7c3aed", "PCA+GROUP" = "#2563eb",
                        "PCA+GHSL" = "#16a34a", "GROUP+GHSL" = "#d97706",
                        "PCA_ONLY" = "#93c5fd", "GROUP_ONLY" = "#fcd34d",
                        "GHSL_ONLY" = "#86efac", "NONE" = "#f3f4f6")
        cons_chr[, track := "Consensus"]
        p <- p +
          geom_tile(data = cons_chr,
                    aes(x = center, y = track, fill = support_signature),
                    height = 0.7,
                    width = (x_max - x_min) / max(50, nrow(cons_chr))) +
          scale_fill_manual(values = c(status_colors, sig_colors),
                            name = "Status / Signature", drop = TRUE)
      }
    }
  }

  p
}

# =============================================================================
# PANEL E — BREAKPOINT STATISTICS
# =============================================================================

build_panel_E <- function() {
  stats_parts <- list()

  # Concordance summary
  if (!is.null(bp_conc) && nrow(bp_conc) > 0) {
    n_both   <- sum(bp_conc$concordant_both, na.rm = TRUE)
    n_left   <- sum(bp_conc$left_only, na.rm = TRUE)
    n_right  <- sum(bp_conc$right_only, na.rm = TRUE)
    n_none   <- sum(bp_conc$neither, na.rm = TRUE)
    stats_parts <- c(stats_parts,
                     paste0("Concordance: both=", n_both, " L=", n_left,
                            " R=", n_right, " none=", n_none))
  }

  # Fisher test results
  if (!is.null(bp_tests) && nrow(bp_tests) > 0) {
    sig_tests <- bp_tests[fisher_signif == TRUE]
    if (nrow(sig_tests) > 0) {
      for (i in seq_len(min(3, nrow(sig_tests)))) {
        r <- sig_tests[i]
        stats_parts <- c(stats_parts,
                         paste0(r$cluster_id, " (", r$boundary_side, "): ",
                                r$group, " enriched ",
                                round(r$prop_in_group * 100), "% vs ",
                                round(r$prop_out_group * 100), "%",
                                " (OR=", r$odds_ratio, ", p=", signif(r$fisher_p, 3), ")"))
      }
    } else {
      stats_parts <- c(stats_parts, "No significant group enrichment (p>0.05)")
    }
  }

  # SV call summary
  if (!is.null(bp_summary) && nrow(bp_summary) > 0) {
    stats_parts <- c(stats_parts,
                     paste0("SV calls: ", bp_summary$n_sv_calls,
                            " (INV=", bp_summary$n_inv_calls,
                            " BND=", bp_summary$n_bnd_calls, ")"),
                     paste0("Clusters: ", bp_summary$n_clusters,
                            " (L=", bp_summary$n_left_clusters,
                            " R=", bp_summary$n_right_clusters, ")"),
                     paste0("Mode: ", bp_summary$primary_mode %||% "N/A"))
  }

  if (length(stats_parts) == 0) {
    stats_parts <- "No breakpoint support data available"
  }

  stat_text <- paste(stats_parts, collapse = "\n")

  p <- ggplot() +
    annotate("text", x = 0, y = 1, label = stat_text,
             hjust = 0, vjust = 1, size = 3, family = "mono",
             color = "#334155") +
    xlim(-0.1, 1) + ylim(0, 1.1) +
    labs(title = "E — Breakpoint Support Statistics") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
          plot.margin = margin(5, 10, 5, 10))

  p
}

# =============================================================================
# BUILD COMPOSITE FIGURE
# =============================================================================

message("[STEP38] Building panels...")

panel_A <- build_panel_A()
panel_B <- build_panel_B()
panel_C <- build_panel_C()
panel_D <- build_panel_D()
panel_E <- build_panel_E()

# Assemble with patchwork — ribbon layout
composite <- panel_A / panel_B / panel_C / panel_D / panel_E +
  plot_layout(heights = c(2, 3, 2.5, 2, 1.5)) +
  plot_annotation(
    title = candidate_header(cid, chr, c_start, c_end),
    subtitle = paste0("Composite Inversion Evidence — ",
                      if (!is.null(bp_summary) && bp_summary$has_breakpoint_support)
                        "Breakpoint support detected" else "No breakpoint support"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "#475569")
    )
  )

# =============================================================================
# SAVE
# =============================================================================

out_pdf <- file.path(plots_dir, paste0(cand_prefix, "_composite_figure.pdf"))
ggsave(out_pdf, composite, width = 10, height = 16, device = "pdf")
message("[STEP38] Saved: ", out_pdf)

# Also save PNG for quick viewing
out_png <- file.path(plots_dir, paste0(cand_prefix, "_composite_figure.png"))
ggsave(out_png, composite, width = 10, height = 16, dpi = 300, device = "png")
message("[STEP38] Saved: ", out_png)

message("\n[DONE] STEP38 composite figure for candidate ", cid)
