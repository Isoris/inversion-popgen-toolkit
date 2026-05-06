#!/usr/bin/env Rscript
# =============================================================================
# 07_plots/plot_breakpoint_te_density.R
# =============================================================================
#
# Reads:
#   output/summary_tables/breakpoint_fragility_summary.tsv
#   output/per_candidate_json/<cid>.json   (for one candidate detail figure)
#
# Produces (PDF + PNG @ 300 dpi):
#   output/plots/01_summary_hotspot_distribution.pdf
#       histogram of left & right breakpoint TE_density_100kb across all
#       candidates, with chrom-bg and local-bg reference lines.
#
#   output/plots/02_class_breakdown_by_candidate.pdf
#       boxplot of fold_vs_chr split by architecture_class A..G.
#
#   output/plots/03_focal_vs_homolog_scatter.pdf
#       Gar breakpoint TE density vs homolog TE density per (candidate, target_species).
#
#   output/plots/04_per_candidate/<cid>__breakpoint_TE_density.pdf
#       per-candidate detail: 50kb / 100kb / 500kb stacked panels with
#       per-class colours and a chrom-background ribbon. ONE per candidate.
#
# This script is intentionally self-contained (no project-wide R utilities).
# Uses base R + ggplot2 + jsonlite. No data.table dependency.
#
# Usage:
#   Rscript plot_breakpoint_te_density.R \
#     --summary  output/summary_tables/breakpoint_fragility_summary.tsv \
#     --json_dir output/per_candidate_json \
#     --out_dir  output/plots
#     [--candidate_id LG28_INV_001]      # produce only that one detail figure
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defs <- list(
    summary      = NA_character_,
    json_dir     = NA_character_,
    out_dir      = NA_character_,
    candidate_id = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key %in% names(defs)) {
      defs[[key]] <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  for (k in c("summary", "json_dir", "out_dir")) {
    if (is.na(defs[[k]])) stop(sprintf("--%s is required", k))
  }
  dir.create(defs$out_dir, recursive = TRUE, showWarnings = FALSE)
  defs
}


## ----------------------------------------------------------------------------
## summary-level plots
## ----------------------------------------------------------------------------
plot_summary <- function(summary_path, out_dir) {
  s <- read.delim(summary_path, sep = "\t", stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "None"))
  if (nrow(s) == 0L) {
    message("[plots] empty summary; skipping summary plots")
    return(invisible(NULL))
  }

  num_cols <- c("left_TE_density_100kb", "right_TE_density_100kb",
                "left_fold_vs_chr",      "right_fold_vs_chr",
                "left_percentile_chr",   "right_percentile_chr",
                "inside_TE_density",     "chrom_TE_density",     "local_TE_density",
                "num_species_with_boundary", "num_species_TE_enriched")
  for (c in intersect(num_cols, colnames(s))) {
    s[[c]] <- suppressWarnings(as.numeric(s[[c]]))
  }

  ## 01 — distribution of breakpoint densities
  long <- data.frame(
    candidate_id = rep(s$candidate_id, 2L),
    side         = rep(c("left", "right"), each = nrow(s)),
    density      = c(s$left_TE_density_100kb, s$right_TE_density_100kb)
  )
  long <- long[is.finite(long$density), , drop = FALSE]

  p1 <- ggplot(long, aes(x = density, fill = side)) +
    geom_histogram(bins = 30L, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = mean(s$chrom_TE_density, na.rm = TRUE),
               linetype = "dashed", colour = "grey20") +
    annotate("text", x = mean(s$chrom_TE_density, na.rm = TRUE),
             y = Inf, vjust = 1.5, label = "mean chrom bg",
             colour = "grey20", size = 3) +
    labs(title = "Breakpoint TE density (100 kb windows)",
         subtitle = sprintf("%d candidates", nrow(s)),
         x = "TE density", y = "candidates", fill = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "01_summary_hotspot_distribution.pdf"),
         p1, width = 8, height = 5)

  ## 02 — fold_vs_chr split by class
  if ("architecture_class" %in% colnames(s)) {
    long2 <- data.frame(
      arch = rep(s$architecture_class, 2L),
      side = rep(c("left", "right"), each = nrow(s)),
      fold = c(s$left_fold_vs_chr, s$right_fold_vs_chr)
    )
    long2 <- long2[is.finite(long2$fold) & !is.na(long2$arch), , drop = FALSE]
    if (nrow(long2) > 0L) {
      long2$arch <- factor(long2$arch, levels = c("A", "B", "C", "D", "E", "F", "G"))
      p2 <- ggplot(long2, aes(x = arch, y = fold, fill = side)) +
        geom_boxplot(outlier.size = 0.7, alpha = 0.7) +
        geom_hline(yintercept = c(1, 2), linetype = "dashed", colour = "grey50") +
        labs(title = "Breakpoint enrichment by architecture class",
             subtitle = "fold vs chromosome background (1 = no enrichment, 2 = high-tier rule)",
             x = "architecture class", y = "fold vs chr", fill = NULL) +
        theme_minimal(base_size = 11)
      ggsave(file.path(out_dir, "02_class_breakdown_by_candidate.pdf"),
             p2, width = 8, height = 5)
    }
  }

  ## 03 — focal vs homolog scatter (Mac + Pangasius columns from summary)
  pairs <- list()
  if (all(c("left_TE_density_100kb", "mac_left_TE_density_100kb") %in% colnames(s))) {
    pairs[[length(pairs) + 1L]] <- data.frame(
      species = "C_macrocephalus",
      gar_left  = s$left_TE_density_100kb,
      gar_right = s$right_TE_density_100kb,
      hom_left  = suppressWarnings(as.numeric(s$mac_left_TE_density_100kb)),
      hom_right = suppressWarnings(as.numeric(s$mac_right_TE_density_100kb))
    )
  }
  if (all(c("pangasius_left_TE_density_100kb") %in% colnames(s))) {
    pairs[[length(pairs) + 1L]] <- data.frame(
      species = "Pangasius_hypophthalmus",
      gar_left  = s$left_TE_density_100kb,
      gar_right = s$right_TE_density_100kb,
      hom_left  = suppressWarnings(as.numeric(s$pangasius_left_TE_density_100kb)),
      hom_right = suppressWarnings(as.numeric(s$pangasius_right_TE_density_100kb))
    )
  }
  if (length(pairs) > 0L) {
    df3 <- do.call(rbind, pairs)
    df3 <- df3[is.finite(df3$hom_left) | is.finite(df3$hom_right), , drop = FALSE]
    long3 <- data.frame(
      species = rep(df3$species, 2L),
      side    = rep(c("left", "right"), each = nrow(df3)),
      gar     = c(df3$gar_left, df3$gar_right),
      hom     = c(df3$hom_left, df3$hom_right)
    )
    long3 <- long3[is.finite(long3$gar) & is.finite(long3$hom), , drop = FALSE]
    if (nrow(long3) > 0L) {
      p3 <- ggplot(long3, aes(x = gar, y = hom, colour = side)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        geom_point(alpha = 0.7, size = 1.6) +
        facet_wrap(~ species) +
        labs(title = "Gar breakpoint TE density vs homologous boundary",
             x = "Gar breakpoint TE density (100 kb)",
             y = "Homologous boundary TE density (100 kb)",
             colour = NULL) +
        theme_minimal(base_size = 11)
      ggsave(file.path(out_dir, "03_focal_vs_homolog_scatter.pdf"),
             p3, width = 9, height = 4.5)
    }
  }

  invisible(NULL)
}


## ----------------------------------------------------------------------------
## per-candidate detail
## ----------------------------------------------------------------------------
candidate_detail <- function(json_path, out_path) {
  j <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  cid <- j$candidate_id

  pull <- function(side, win) {
    blk <- j$repeat_density[[side]][[win]]
    if (is.null(blk)) return(NULL)
    cls <- blk$te_density_by_class
    if (is.null(cls) || length(cls) == 0L) return(NULL)
    data.frame(
      side  = side, window = win,
      class = names(cls), density = unlist(cls, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  for (side in c("left_breakpoint", "right_breakpoint")) {
    for (win in c("window_50kb", "window_100kb", "window_500kb")) {
      d <- pull(side, win)
      if (!is.null(d)) rows[[length(rows) + 1L]] <- d
    }
  }
  if (length(rows) == 0L) {
    message("[plots] no per-class density for ", cid, "; skip")
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)
  df$window <- factor(df$window, levels = c("window_50kb", "window_100kb", "window_500kb"))
  df$side   <- factor(df$side,   levels = c("left_breakpoint", "right_breakpoint"))

  chrom_bg <- j$repeat_density$chromosome_background$te_density
  if (is.null(chrom_bg) || !is.finite(chrom_bg)) chrom_bg <- NA_real_

  p <- ggplot(df, aes(x = window, y = density, fill = class)) +
    geom_col(position = "stack") +
    geom_hline(yintercept = chrom_bg, linetype = "dashed", colour = "grey20") +
    facet_wrap(~ side) +
    labs(title = sprintf("%s — breakpoint TE density (per-class stacked)", cid),
         subtitle = sprintf("dashed = chromosome background TE density (%.3f)", chrom_bg),
         x = NULL, y = "TE density (fraction of bp)", fill = "TE class") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 25, hjust = 1))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = 9, height = 5)
  invisible(NULL)
}


## ----------------------------------------------------------------------------
## main
## ----------------------------------------------------------------------------
main <- function() {
  a <- parse_args()

  message("[plots] summary plots from ", a$summary)
  plot_summary(a$summary, a$out_dir)

  per_cand_dir <- file.path(a$out_dir, "04_per_candidate")
  dir.create(per_cand_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.na(a$candidate_id)) {
    p <- file.path(a$json_dir, paste0(a$candidate_id, ".json"))
    if (file.exists(p)) {
      candidate_detail(p, file.path(per_cand_dir,
                                    paste0(a$candidate_id, "__breakpoint_TE_density.pdf")))
    } else {
      message("[plots] candidate JSON not found: ", p)
    }
  } else {
    files <- list.files(a$json_dir, pattern = "\\.json$", full.names = TRUE)
    for (p in files) {
      cid <- sub("\\.json$", "", basename(p))
      candidate_detail(p, file.path(per_cand_dir,
                                    paste0(cid, "__breakpoint_TE_density.pdf")))
    }
  }
  message("[plots] done")
}


if (sys.nframe() == 0L) main()
