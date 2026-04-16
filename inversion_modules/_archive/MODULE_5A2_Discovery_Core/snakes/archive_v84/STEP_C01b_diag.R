#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01b_diag.R  (v8.4)
#
# DIAGNOSTIC PLOTS for C01b core + merge results.
# Standalone -- reads snake_*.tsv.gz + .precomp.rds
#
# Per-chromosome (PDF with one page per chr):
#   01_ideogram_tracks.pdf           -- layered tracks: z/inv/S/M/L/mergeA/B/C
#   02_mds_by_core_family.pdf        -- MDS1x2 colored by S1S/S1M/S1L
#   03_mds_by_merge_region.pdf       -- MDS1x2 colored by merge_A region ID
#   05_nested_boundary_ladder.pdf    -- interval nesting: S/M/L + merge A/B/C
#   07_merge_density_profile.pdf     -- density + spread_class along chromosome
#   08_fuzzy_merge_scores.pdf        -- per-pair fuzzy scores: M, G, combined
#   09_score_decomposition.pdf       -- membership vs geometric scatter per chr
#
# Genome-wide summary (PNG):
#   S01_core_census.png              -- cores per family per chr (stacked bar)
#   S02_merge_census.png             -- merge regions per tier per chr
#   S03_agreement_heatmap.png        -- three-way comparison: agree/composite/weak
#   S04_core_size_violin.png         -- core size distribution by family (genome)
#   S05_density_vs_coherence.png     -- scatter: density vs coherence per region
#   S06_stop_reason_breakdown.png    -- why cores stopped (stacked bar)
#   S07_spread_class_dist.png        -- spread class distribution per merge tier
#   S08_coverage_barplot.png         -- % chromosome covered by cores per chr
#   S09_fuzzy_score_dist.png         -- fuzzy score histogram by tier + merge/reject
#   S10_membership_vs_geometric.png  -- M vs G scatter, genome-wide, by decision
#
# Usage:
#   Rscript STEP_C01b_diag.R <precomp_dir> <snake_outdir> <outdir> [chrom]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript STEP_C01b_diag.R <precomp_dir> <snake_outdir> <outdir> [chrom]")

precomp_dir <- args[1]
snake_dir   <- args[2]
outdir      <- args[3]
chrom_filter <- if (length(args) >= 4 && args[4] != "all") args[4] else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DPI <- 350
THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555"),
        plot.caption = element_text(size = 6, color = "#888", hjust = 0),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

safe_load <- function(f) {
  if (file.exists(f)) tryCatch(fread(f), error = function(e) data.table()) else data.table()
}

# Core/merge color palettes
PAL_FAMILY <- c("1S" = "#1d4ed8", "1M" = "#16a34a", "1L" = "#d97706",
                "merge_A" = "#dc2626", "merge_B" = "#7c3aed", "merge_C" = "#059669")
PAL_STATUS <- c("seed" = "#000000", "accepted" = "#2563eb", "tolerated" = "#f59e0b",
                "bridge_band" = "#06b6d4", "bridge_fuzzy" = "#06b6d4",
                "bridge_skip" = "#67e8f9", "background" = "#e5e7eb")
PAL_SPREAD <- c("dense_continuous" = "#16a34a", "moderate_gaps" = "#f59e0b",
                "sparse_scattered" = "#dc2626", "tiny_fragment" = "#9ca3af")
PAL_AGREE <- c("all_agree" = "#16a34a", "AB_agree_C_split_COMPOSITE" = "#f59e0b",
               "AB_agree_C_none" = "#93c5fd", "A_only_C_exists" = "#c084fc",
               "A_only_weak" = "#dc2626")

# Region ID palette (cycling)
REGION_PAL <- c("#dc2626", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
                "#059669", "#e11d48", "#0284c7", "#65a30d", "#c026d3",
                "#0891b2", "#ca8a04", "#4f46e5", "#15803d", "#be185d")

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01b_diag] Loading data...")

# Snake outputs
win_dt   <- safe_load(file.path(snake_dir, "snake_windows.tsv.gz"))
reg_dt   <- safe_load(file.path(snake_dir, "snake_regions.tsv.gz"))
hier_dt  <- safe_load(file.path(snake_dir, "snake_hierarchy.tsv.gz"))
state_dt <- safe_load(file.path(snake_dir, "snake_window_states.tsv.gz"))
summ_dt  <- safe_load(file.path(snake_dir, "snake_summary.tsv"))
log_dt   <- safe_load(file.path(snake_dir, "snake_decision_log.tsv.gz"))
comp_dt  <- safe_load(file.path(snake_dir, "snake_multiscale_comparison.tsv.gz"))
diag_dt  <- safe_load(file.path(snake_dir, "snake_diagnostics.tsv.gz"))

# Merge outputs (separate merge script)
mreg_dt  <- safe_load(file.path(snake_dir, "snake_merge_regions.tsv.gz"))
mcomp_dt <- safe_load(file.path(snake_dir, "snake_merge_comparison.tsv.gz"))
mscore_dt <- safe_load(file.path(snake_dir, "snake_merge_scores.tsv.gz"))

if (nrow(mscore_dt) > 0) message("[C01b_diag] Fuzzy merge scores: ", nrow(mscore_dt), " pair evaluations")

# Use merge regions if available, otherwise fall back to snake_regions
if (nrow(mreg_dt) > 0) {
  message("[C01b_diag] Using merge regions from standalone merge script")
  merge_reg <- mreg_dt
} else if (nrow(reg_dt) > 0) {
  merge_reg <- reg_dt[snake_phase == "merged"]
} else {
  merge_reg <- data.table()
}

core_reg <- if (nrow(reg_dt) > 0) reg_dt[snake_phase == "core"] else data.table()

# Precomp (for sim_mat, MDS)
rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
precomp_list <- list()
for (f in rds_files) {
  obj <- readRDS(f)
  precomp_list[[obj$chrom]] <- obj
}

chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

message("[C01b_diag] ", length(chroms), " chromosomes | ",
        nrow(core_reg), " core regions | ", nrow(merge_reg), " merge regions")

# =============================================================================
# PER-CHROMOSOME PLOTS
# =============================================================================

message("[C01b_diag] Generating per-chromosome plots...")

# -- 01: Ideogram tracks ---------------------------------------------
build_ideogram <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  chr_len_mb <- max(dt$end_bp) / 1e6

  track_rows <- list()

  # Track 1: robust z (background heatmap)
  if ("robust_z" %in% names(dt)) {
    z_class <- fifelse(abs(dt$robust_z) >= 5, "z>=5",
               fifelse(abs(dt$robust_z) >= 3, "z>=3",
               fifelse(abs(dt$robust_z) >= 2, "z>=2", "z<2")))
    track_rows[[1]] <- data.table(
      pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
      track = "01_robust_z", status = z_class)
  }

  # Track 2: inv_likeness
  if ("inv_likeness" %in% names(dt)) {
    il_class <- fifelse(dt$inv_likeness >= 0.90, "inv>=0.90",
                fifelse(dt$inv_likeness >= 0.75, "inv>=0.75",
                fifelse(dt$inv_likeness >= 0.50, "inv>=0.50", "inv<0.50")))
    track_rows[[length(track_rows) + 1]] <- data.table(
      pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
      track = "02_inv_likeness", status = il_class)
  }

  # Tracks 3-5: cores by family
  if (nrow(core_reg) > 0) {
    for (fam in c("1S", "1M", "1L")) {
      fam_reg <- core_reg[chrom == chr & core_family == fam]
      if (nrow(fam_reg) > 0) {
        fam_wins <- win_dt[chrom == chr & core_family == fam]
        if (nrow(fam_wins) > 0) {
          track_rows[[length(track_rows) + 1]] <- data.table(
            pos_mb = (fam_wins$start_bp + fam_wins$end_bp) / 2e6,
            track = paste0("0", 2 + which(c("1S","1M","1L") == fam), "_core_", fam),
            status = fam_wins$inclusion_status)
        }
      }
    }
  }

  # Tracks 6-8: merge regions A/B/C (shaded blocks)
  if (nrow(merge_reg) > 0) {
    for (mfam in c("merge_A", "merge_B", "merge_C")) {
      mf_reg <- merge_reg[chrom == chr & merge_family == mfam]
      if (nrow(mf_reg) > 0) {
        tnum <- switch(mfam, "merge_A" = "06", "merge_B" = "07", "merge_C" = "08")
        for (ri in seq_len(nrow(mf_reg))) {
          r <- mf_reg[ri]
          # Create points spanning the region
          span_points <- seq(r$start_bp, r$end_bp, length.out = max(3, r$n_windows %||% 5))
          spread <- r$spread_class %||% "unknown"
          track_rows[[length(track_rows) + 1]] <- data.table(
            pos_mb = span_points / 1e6,
            track = paste0(tnum, "_", mfam),
            status = spread)
        }
      }
    }
  }

  if (length(track_rows) == 0) return(NULL)
  track_dt <- rbindlist(track_rows, fill = TRUE)

  all_colors <- c(
    PAL_STATUS, PAL_SPREAD,
    "z>=5" = "#1e3a5f", "z>=3" = "#3b82f6", "z>=2" = "#93c5fd", "z<2" = "#f3f4f6",
    "inv>=0.90" = "#dc2626", "inv>=0.75" = "#f59e0b", "inv>=0.50" = "#fbbf24", "inv<0.50" = "#f3f4f6",
    "unknown" = "#d4d4d4"
  )

  tile_w <- chr_len_mb / n * 0.9

  p <- ggplot(track_dt, aes(x = pos_mb, y = track, fill = status)) +
    geom_tile(height = 0.7, width = tile_w) +
    scale_fill_manual(values = all_colors, na.value = "#f3f4f6", name = "Status") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " -- Snake 1 Ideogram"),
         subtitle = paste0(n, " windows, ", round(chr_len_mb, 1), " Mb"),
         caption = paste0("Tracks: z-score | inv-likeness | S1S/M/L cores | merge A/B/C\n",
                         "Core colors: seed/accepted/tolerated | Merge colors: spread class")) +
    THEME_BASE +
    theme(axis.text.y = element_text(size = 6))
  p
}

# -- 02: MDS by core family ------------------------------------------
build_mds_cores <- function(chr) {
  mds_mat <- precomp_list[[chr]]$mds_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     gid = dt$global_window_id, family = "background")

  chr_win <- win_dt[chrom == chr & snake_phase == "core"]
  if (nrow(chr_win) > 0) {
    for (fam in c("1S", "1M", "1L")) {
      fam_gids <- chr_win[core_family == fam]$global_window_id
      pdt[gid %in% fam_gids, family := fam]
    }
  }
  pdt[, family := factor(family, levels = c("background", "1L", "1M", "1S"))]

  p <- ggplot(pdt[order(family)], aes(x = MDS1, y = MDS2, color = family)) +
    geom_point(size = ifelse(pdt$family[order(pdt$family)] == "background", 0.5, 1.5),
               alpha = ifelse(pdt$family[order(pdt$family)] == "background", 0.2, 0.8)) +
    scale_color_manual(values = c("background" = "#e5e7eb", PAL_FAMILY[c("1S","1M","1L")])) +
    labs(title = paste0(chr, " -- MDS by Core Family"),
         subtitle = paste0("S:", sum(pdt$family == "1S"),
                          " M:", sum(pdt$family == "1M"),
                          " L:", sum(pdt$family == "1L"),
                          " bg:", sum(pdt$family == "background")),
         caption = "Blue: strict (S1S) | Green: moderate (S1M) | Orange: broad (S1L)") +
    THEME_BASE
  p
}

# -- 03: MDS by merge region ID --------------------------------------
build_mds_merge <- function(chr) {
  mds_mat <- precomp_list[[chr]]$mds_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     gid = dt$global_window_id, region = "background")

  # Color by merge_A region
  ma_regs <- merge_reg[chrom == chr & merge_family == "merge_A"]
  ma_wins <- if (nrow(win_dt) > 0) {
    win_dt[chrom == chr & merge_family == "merge_A"]
  } else data.table()

  if (nrow(ma_wins) > 0) {
    for (sid in unique(ma_wins$snake_id)) {
      gids <- ma_wins[snake_id == sid]$global_window_id
      pdt[gid %in% gids, region := paste0("R", sid)]
    }
  }

  # Generate colors for regions
  region_ids <- setdiff(unique(pdt$region), "background")
  if (length(region_ids) == 0) return(NULL)
  rcols <- setNames(
    REGION_PAL[(seq_along(region_ids) - 1L) %% length(REGION_PAL) + 1L],
    region_ids
  )
  rcols <- c("background" = "#f3f4f6", rcols)

  p <- ggplot(pdt[order(region == "background", decreasing = TRUE)],
              aes(x = MDS1, y = MDS2, color = region)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_manual(values = rcols) +
    labs(title = paste0(chr, " -- MDS by Merge-A Region"),
         subtitle = paste0(length(region_ids), " merge_A regions"),
         caption = "Each color = one merge_A region\nGrey = unclaimed background") +
    THEME_BASE +
    guides(color = guide_legend(ncol = 4, override.aes = list(size = 3)))
  p
}

# -- 05: Nested boundary ladder --------------------------------------
build_ladder <- function(chr) {
  chr_core <- core_reg[chrom == chr]
  chr_merge <- merge_reg[chrom == chr]
  if (nrow(chr_core) == 0 && nrow(chr_merge) == 0) return(NULL)

  rows <- list()
  # Cores
  if (nrow(chr_core) > 0) {
    for (ri in seq_len(nrow(chr_core))) {
      r <- chr_core[ri]
      rows[[length(rows) + 1]] <- data.table(
        start_mb = r$start_bp / 1e6, end_mb = r$end_bp / 1e6,
        level = paste0("core_", r$core_family),
        id = r$snake_id, n_win = r$n_windows)
    }
  }
  # Merge regions
  if (nrow(chr_merge) > 0) {
    for (ri in seq_len(nrow(chr_merge))) {
      r <- chr_merge[ri]
      rows[[length(rows) + 1]] <- data.table(
        start_mb = r$start_bp / 1e6, end_mb = r$end_bp / 1e6,
        level = r$merge_family,
        id = r$snake_id, n_win = r$n_windows)
    }
  }

  ldt <- rbindlist(rows)
  level_order <- c("core_1S", "core_1M", "core_1L", "merge_C", "merge_B", "merge_A")
  ldt[, level := factor(level, levels = rev(level_order))]

  pal <- c("core_1S" = "#1d4ed8", "core_1M" = "#16a34a", "core_1L" = "#d97706",
           "merge_A" = "#dc2626", "merge_B" = "#7c3aed", "merge_C" = "#059669")

  p <- ggplot(ldt, aes(xmin = start_mb, xmax = end_mb, ymin = as.numeric(level) - 0.4,
                        ymax = as.numeric(level) + 0.4, fill = level)) +
    geom_rect(alpha = 0.7, color = "#333", linewidth = 0.15) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(breaks = seq_along(level_order), labels = rev(level_order)) +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " -- Nested Boundary Ladder"),
         subtitle = paste0("Cores: S=", sum(ldt$level == "core_1S"),
                          " M=", sum(ldt$level == "core_1M"),
                          " L=", sum(ldt$level == "core_1L"),
                          " | Merge: A=", sum(ldt$level == "merge_A"),
                          " B=", sum(ldt$level == "merge_B"),
                          " C=", sum(ldt$level == "merge_C")),
         caption = "Top: finest (S1S cores) -> Bottom: coarsest (merge_A)\nNesting shows how cores consolidate into merge regions") +
    THEME_BASE
  p
}

# -- 07: Merge density profile ---------------------------------------
build_density_profile <- function(chr) {
  chr_merge <- merge_reg[chrom == chr & merge_family == "merge_A"]
  if (nrow(chr_merge) == 0) return(NULL)
  if (!"density" %in% names(chr_merge)) return(NULL)

  chr_merge[, mid_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(chr_merge, aes(x = mid_mb, y = density, fill = spread_class)) +
    geom_col(width = chr_merge$span_mb * 0.9) +
    scale_fill_manual(values = PAL_SPREAD) +
    geom_hline(yintercept = c(0.30, 0.70), linetype = "dashed", color = "#666") +
    labs(x = paste0(chr, " (Mb)"), y = "Density",
         title = paste0(chr, " -- Merge-A Density Profile"),
         subtitle = paste0(nrow(chr_merge), " regions | ",
                          "dense:", sum(chr_merge$spread_class == "dense_continuous"),
                          " mod:", sum(chr_merge$spread_class == "moderate_gaps"),
                          " sparse:", sum(chr_merge$spread_class == "sparse_scattered"),
                          " tiny:", sum(chr_merge$spread_class == "tiny_fragment")),
         caption = "Density = windows_in_region / possible_windows\nDashed: 0.30 (sparse/moderate) and 0.70 (moderate/dense) thresholds") +
    THEME_BASE
  p
}

# -- 08: Fuzzy merge score profile ------------------------------------
build_fuzzy_profile <- function(chr) {
  if (nrow(mscore_dt) == 0) return(NULL)
  chr_scores <- mscore_dt[chrom == chr]
  if (nrow(chr_scores) == 0) return(NULL)

  chr_scores[, pos_mb := (core_end_bp + next_start_bp) / 2e6]

  p <- ggplot(chr_scores, aes(x = pos_mb)) +
    # Combined fuzzy score
    geom_point(aes(y = fuzzy_score, color = decision), size = 1.2, alpha = 0.7) +
    # Membership and geometric as separate layers
    geom_point(aes(y = membership_score), shape = 1, size = 0.8,
               color = "#7c3aed", alpha = 0.4) +
    geom_point(aes(y = geometric_score), shape = 2, size = 0.8,
               color = "#0891b2", alpha = 0.4) +
    scale_color_manual(values = c("merge" = "#16a34a", "reject" = "#dc2626")) +
    facet_wrap(~ merge_family, ncol = 1) +
    geom_hline(data = data.frame(
      merge_family = c("merge_A", "merge_B", "merge_C"),
      thresh = c(0.30, 0.50, 0.70)
    ), aes(yintercept = thresh), linetype = "dashed", color = "#666", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "Score",
         title = paste0(chr, " -- Fuzzy Merge Scores"),
         subtitle = paste0(nrow(chr_scores), " core-pair evaluations\n",
                          "Circle: membership | Triangle: geometric | Filled: combined"),
         caption = paste0("Green=merged | Red=rejected | Dashed=accept threshold\n",
                         "merge_A: MAX(M,G)>=0.30 | merge_B: mean(M,G)>=0.50 | ",
                         "merge_C: MIN(M,G)>=0.70")) +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold"))
  p
}

# -- 09: Merge score decomposition scatter ----------------------------
build_score_scatter <- function(chr) {
  if (nrow(mscore_dt) == 0) return(NULL)
  chr_scores <- mscore_dt[chrom == chr]
  if (nrow(chr_scores) == 0) return(NULL)

  chr_scores <- chr_scores[is.finite(membership_score) & is.finite(geometric_score)]
  if (nrow(chr_scores) == 0) return(NULL)

  p <- ggplot(chr_scores, aes(x = membership_score, y = geometric_score,
                               color = decision, shape = merge_family)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = c("merge" = "#16a34a", "reject" = "#dc2626")) +
    # Threshold regions
    annotate("rect", xmin = 0, xmax = 0.30, ymin = 0, ymax = 0.30,
             fill = "#fee2e2", alpha = 0.3) +
    annotate("rect", xmin = 0.70, xmax = 1, ymin = 0.70, ymax = 1,
             fill = "#dcfce7", alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#888") +
    labs(x = "Membership score (fuzzy composition)",
         y = "Geometric score (sim_mat continuity)",
         title = paste0(chr, " -- Merge Score Decomposition"),
         subtitle = paste0("Each point = one core pair evaluation\n",
                          sum(chr_scores$decision == "merge"), " merged, ",
                          sum(chr_scores$decision == "reject"), " rejected"),
         caption = paste0("Red zone: both scores low -> always reject\n",
                         "Green zone: both scores high -> always merge\n",
                         "Diagonal: where M and G contribute equally")) +
    THEME_BASE
  p
}

# =============================================================================
# RENDER PER-CHROMOSOME PLOTS
# =============================================================================

plot_builders <- list(
  "01_ideogram_tracks" = build_ideogram,
  "02_mds_by_core_family" = build_mds_cores,
  "03_mds_by_merge_region" = build_mds_merge,
  "05_nested_boundary_ladder" = build_ladder,
  "07_merge_density_profile" = build_density_profile,
  "08_fuzzy_merge_scores" = build_fuzzy_profile,
  "09_score_decomposition" = build_score_scatter
)

plot_dims <- list(
  "01_ideogram_tracks" = c(w = 16, h = 7),
  "02_mds_by_core_family" = c(w = 8, h = 7),
  "03_mds_by_merge_region" = c(w = 8, h = 7),
  "05_nested_boundary_ladder" = c(w = 14, h = 5),
  "07_merge_density_profile" = c(w = 14, h = 5),
  "08_fuzzy_merge_scores" = c(w = 14, h = 10),
  "09_score_decomposition" = c(w = 8, h = 7)
)

for (pname in names(plot_builders)) {
  builder <- plot_builders[[pname]]
  dims <- plot_dims[[pname]]

  plots <- list()
  for (chr in chroms) {
    p <- tryCatch(builder(chr), error = function(e) {
      message("  [WARN] ", pname, " / ", chr, ": ", e$message); NULL
    })
    if (!is.null(p)) plots[[chr]] <- p
  }

  if (length(plots) > 0) {
    f_pdf <- file.path(outdir, paste0(pname, ".pdf"))
    pdf(f_pdf, width = dims["w"], height = dims["h"])
    for (p in plots) print(p)
    dev.off()

    png_dir <- file.path(outdir, "png", pname)
    dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
    for (chr in names(plots)) {
      ggsave(file.path(png_dir, paste0(chr, ".png")),
             plots[[chr]], width = dims["w"], height = dims["h"], dpi = DPI)
    }
    message("[C01b_diag] ", pname, ": ", length(plots), " pages -> ", f_pdf)
  }
}

# =============================================================================
# GENOME-WIDE SUMMARY PLOTS
# =============================================================================

message("[C01b_diag] Generating genome-wide summaries...")
summ_dir <- file.path(outdir, "summaries")
dir.create(summ_dir, recursive = TRUE, showWarnings = FALSE)

# -- S01: Core census (stacked bar) ----------------------------------
if (nrow(core_reg) > 0) {
  cc <- core_reg[, .N, by = .(chrom, core_family)]
  cc[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(cc, aes(x = chrom, y = N, fill = core_family)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = PAL_FAMILY) +
    labs(x = NULL, y = "Number of cores",
         title = "Core Census per Chromosome",
         subtitle = paste0("Total: ", nrow(core_reg), " cores | ",
                          "S:", sum(core_reg$core_family == "1S"),
                          " M:", sum(core_reg$core_family == "1M"),
                          " L:", sum(core_reg$core_family == "1L")),
         caption = "Stacked: S1S (blue) + S1M (green) + S1L (orange)") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S01_core_census.png"), p, width = 14, height = 6, dpi = DPI)
  message("[C01b_diag] S01_core_census.png")
}

# -- S02: Merge census -----------------------------------------------
if (nrow(merge_reg) > 0) {
  mc <- merge_reg[, .N, by = .(chrom, merge_family)]
  mc[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(mc, aes(x = chrom, y = N, fill = merge_family)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = PAL_FAMILY) +
    labs(x = NULL, y = "Merge regions",
         title = "Merge Region Census per Chromosome",
         subtitle = paste0("A:", sum(merge_reg$merge_family == "merge_A"),
                          " B:", sum(merge_reg$merge_family == "merge_B"),
                          " C:", sum(merge_reg$merge_family == "merge_C")),
         caption = "A should have FEWEST (generous), C should have MOST (strict)\nIf inverted -> merge not working properly") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S02_merge_census.png"), p, width = 14, height = 6, dpi = DPI)
  message("[C01b_diag] S02_merge_census.png")
}

# -- S03: Three-way agreement heatmap --------------------------------
use_comp <- if (nrow(mcomp_dt) > 0) mcomp_dt else comp_dt
if (nrow(use_comp) > 0 && "agreement" %in% names(use_comp)) {
  ag <- use_comp[, .N, by = .(chrom, agreement)]
  ag[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(ag, aes(x = chrom, y = N, fill = agreement)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = PAL_AGREE) +
    labs(x = NULL, y = "Regions",
         title = "Three-Way Merge Agreement per Chromosome",
         subtitle = paste0("Total: ", nrow(use_comp), " merge_A regions compared"),
         caption = paste0("Green=all_agree (high confidence) | ",
                         "Orange=COMPOSITE (Snake2 needed) | ",
                         "Red=A_only_weak")) +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S03_agreement_heatmap.png"), p, width = 14, height = 6, dpi = DPI)
  message("[C01b_diag] S03_agreement_heatmap.png")
}

# -- S04: Core size violin -------------------------------------------
if (nrow(core_reg) > 0) {
  p <- ggplot(core_reg, aes(x = core_family, y = pmin(n_windows, 50), fill = core_family)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_fill_manual(values = PAL_FAMILY) +
    labs(x = "Family", y = "Core size (windows, capped 50)",
         title = "Core Size Distribution by Family",
         subtitle = paste0("S median=", median(core_reg[core_family == "1S"]$n_windows),
                          " | M median=", median(core_reg[core_family == "1M"]$n_windows),
                          " | L median=", median(core_reg[core_family == "1L"]$n_windows)),
         caption = "Larger cores = stronger signal\nS1S expected smallest (strict), S1L largest (broad)") +
    THEME_BASE

  ggsave(file.path(summ_dir, "S04_core_size_violin.png"), p, width = 8, height = 6, dpi = DPI)
  message("[C01b_diag] S04_core_size_violin.png")
}

# -- S05: Density vs coherence scatter --------------------------------
if (nrow(merge_reg) > 0 && all(c("density", "coherence") %in% names(merge_reg))) {
  p <- ggplot(merge_reg, aes(x = density, y = coherence, color = merge_family, size = n_windows)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = PAL_FAMILY) +
    scale_size_continuous(range = c(0.5, 4)) +
    geom_vline(xintercept = c(0.30, 0.70), linetype = "dashed", color = "#888") +
    geom_hline(yintercept = c(0.15, 0.25, 0.40), linetype = "dotted", color = "#888") +
    labs(x = "Density (windows / possible)", y = "Coherence",
         title = "Merge Region: Density vs Coherence",
         subtitle = paste0(nrow(merge_reg), " regions across all chromosomes"),
         caption = paste0("Vertical dashed: density thresholds (sparse/moderate/dense)\n",
                         "Horizontal dotted: coherence floors (A=0.15, B=0.25, C=0.40)\n",
                         "Size = number of windows in region")) +
    THEME_BASE

  ggsave(file.path(summ_dir, "S05_density_vs_coherence.png"), p, width = 10, height = 7, dpi = DPI)
  message("[C01b_diag] S05_density_vs_coherence.png")
}

# -- S06: Stop reason breakdown --------------------------------------
if (nrow(log_dt) > 0 && "direction" %in% names(log_dt)) {
  stop_reasons <- log_dt[phase == "core" & grepl("right|left", direction, ignore.case = TRUE)]
  if (nrow(stop_reasons) > 0 && "decision" %in% names(stop_reasons)) {
    sr <- stop_reasons[, .N, by = .(decision)]

    p <- ggplot(sr, aes(x = reorder(decision, -N), y = N)) +
      geom_col(fill = "#3b82f6") +
      coord_flip() +
      labs(x = NULL, y = "Count",
           title = "Core Stop Reasons (All Chromosomes)",
           subtitle = paste0(nrow(stop_reasons), " stop events total"),
           caption = "gap_damage = gap too large or too much accumulated damage\nrejected = score below threshold") +
      THEME_BASE

    ggsave(file.path(summ_dir, "S06_stop_reason_breakdown.png"), p, width = 10, height = 6, dpi = DPI)
    message("[C01b_diag] S06_stop_reason_breakdown.png")
  }
}

# -- S07: Spread class distribution ----------------------------------
if (nrow(merge_reg) > 0 && "spread_class" %in% names(merge_reg)) {
  sp <- merge_reg[, .N, by = .(merge_family, spread_class)]

  p <- ggplot(sp, aes(x = merge_family, y = N, fill = spread_class)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = PAL_SPREAD) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Merge family", y = "Proportion",
         title = "Spread Class Distribution by Merge Tier",
         subtitle = paste0("Total regions: A=", sum(merge_reg$merge_family == "merge_A"),
                          " B=", sum(merge_reg$merge_family == "merge_B"),
                          " C=", sum(merge_reg$merge_family == "merge_C")),
         caption = "dense_continuous (green) = best candidates\nsparse_scattered (red) = likely noise or very fragmented") +
    THEME_BASE

  ggsave(file.path(summ_dir, "S07_spread_class_dist.png"), p, width = 8, height = 6, dpi = DPI)
  message("[C01b_diag] S07_spread_class_dist.png")
}

# -- S08: Coverage barplot --------------------------------------------
cov_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  n_total <- nrow(dt)
  chr_cores <- core_reg[chrom == chr]
  if (nrow(chr_cores) > 0 && nrow(win_dt) > 0) {
    core_gids <- unique(win_dt[chrom == chr & snake_phase == "core"]$global_window_id)
    n_covered <- length(core_gids)
  } else {
    n_covered <- 0L
  }
  cov_rows[[chr]] <- data.table(
    chrom = chr, n_total = n_total, n_covered = n_covered,
    pct = round(100 * n_covered / n_total, 1)
  )
}
if (length(cov_rows) > 0) {
  cov_dt <- rbindlist(cov_rows)
  cov_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(cov_dt, aes(x = chrom, y = pct)) +
    geom_col(fill = "#3b82f6") +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.3, size = 2.5) +
    labs(x = NULL, y = "Coverage (%)",
         title = "Chromosome Coverage by Snake 1 Cores",
         subtitle = paste0("Genome-wide: ",
                          round(100 * sum(cov_dt$n_covered) / sum(cov_dt$n_total), 1), "%"),
         caption = "High coverage may indicate too-permissive thresholds\nExpect 5-30% for chromosomes with inversions") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S08_coverage_barplot.png"), p, width = 14, height = 6, dpi = DPI)
  message("[C01b_diag] S08_coverage_barplot.png")
}

# -- S09: Fuzzy score distribution by tier ----------------------------
if (nrow(mscore_dt) > 0) {
  p <- ggplot(mscore_dt, aes(x = fuzzy_score, fill = decision)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "stack") +
    scale_fill_manual(values = c("merge" = "#16a34a", "reject" = "#dc2626")) +
    facet_wrap(~ merge_family, ncol = 1, scales = "free_y") +
    geom_vline(data = data.frame(
      merge_family = c("merge_A", "merge_B", "merge_C"),
      thresh = c(0.30, 0.50, 0.70)
    ), aes(xintercept = thresh), linetype = "dashed", linewidth = 0.5) +
    labs(x = "Fuzzy merge score", y = "Core pair count",
         title = "Fuzzy Merge Score Distribution by Tier",
         subtitle = paste0(nrow(mscore_dt), " core-pair evaluations across all chromosomes\n",
                          "Merged: ", sum(mscore_dt$decision == "merge"),
                          " | Rejected: ", sum(mscore_dt$decision == "reject")),
         caption = paste0("Dashed = accept threshold per tier\n",
                         "Scores near threshold = boundary cases for tuning")) +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold"))

  ggsave(file.path(summ_dir, "S09_fuzzy_score_dist.png"), p, width = 10, height = 10, dpi = DPI)
  message("[C01b_diag] S09_fuzzy_score_dist.png")
}

# -- S10: Membership vs geometric score scatter (genome-wide) --------
if (nrow(mscore_dt) > 0) {
  ms_valid <- mscore_dt[is.finite(membership_score) & is.finite(geometric_score)]
  if (nrow(ms_valid) > 0) {
    p <- ggplot(ms_valid, aes(x = membership_score, y = geometric_score,
                               color = merge_family)) +
      geom_point(size = 0.8, alpha = 0.4) +
      scale_color_manual(values = PAL_FAMILY) +
      facet_wrap(~ decision, ncol = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#888") +
      labs(x = "Membership score (fuzzy composition)",
           y = "Geometric score (sim_mat)",
           title = "Membership vs Geometric Score (Genome-Wide)",
           subtitle = paste0(nrow(ms_valid), " evaluations | Left: merged | Right: rejected"),
           caption = paste0("Points above diagonal: membership drives the merge\n",
                           "Points below diagonal: geometry drives the merge\n",
                           "Corner clusters reveal which relation is more informative")) +
      THEME_BASE

    ggsave(file.path(summ_dir, "S10_membership_vs_geometric.png"), p, width = 12, height = 6, dpi = DPI)
    message("[C01b_diag] S10_membership_vs_geometric.png")
  }
}

message("\n[DONE] STEP_C01b_diag complete -> ", outdir)
