#!/usr/bin/env Rscript

# =============================================================================
# STEP10f_v3_snake_diagnostics.R  (v8.3.1)
#
# COMPREHENSIVE DIAGNOSTIC PLOTS + TABLES for the multi-scale snake system.
#
# THREE LEVELS:
#   Level 1 — WINDOW:  what each snake sees and decides
#   Level 2 — REGION:  how snakes build, merge, split, fail
#   Level 3 — SUMMARY: one multi-head table per candidate + census + stop-reasons
#
# OUTPUTS (numbered for clean organization):
#   01_chromosome_track_ideogram.pdf     — layered horizontal tracks (MOST USEFUL PLOT)
#   02_mds_scatter_by_snake.pdf          — MDS1×2 colored by snake membership
#   02b_mds_scatter_by_merged_id.pdf     — MDS1×2 colored by merge_A region ID
#   02c_mds_scatter_by_regime.pdf        — MDS1×2 colored by Snake 4 regime
#   03_continuity_score_profiles.pdf     — per-snake-run score vs threshold
#   04_seed_support_diagnostic.pdf       — seed NN distance vs position
#   05_nested_boundary_ladders.pdf       — S/M/L + mergeA/B interval nesting
#   06_region_coherence_heatmaps.pdf     — within-region similarity matrices
#   07_merge_bridge_audits.pdf           — per-bridge continuity strips
#   08_candidate_summary_table.tsv.gz    — multi-head summary (THE TABLE)
#   09_chromosome_census.tsv             — per-chromosome counts
#   10_stop_reason_table.tsv             — why each snake stopped
#   10_stop_reason_barplot.pdf           — stacked bars by snake type
#   11_profile_diversity_<chr>.pdf       — per-stripe Hamming + concordance curves
#
# NEW TRACKS IN IDEOGRAM (when data available):
#   09_Snake2  — Snake 2 PASS/WEAK/FAIL per window
#   10_Snake3  — Snake 3 PASS/WEAK/FAIL per window
#   11_Snake4  — Snake 4 PASS/PARTIAL/WEAK/FAIL per window
#   12_S4_regime — Snake 4 middle-core regime identity (colored by regime)
#   13_control — 500kb control candidates (background comparison)
#
# INPUTS (from STEP10e_v3 + STEP10_v2):
#   <snake_dir>/snake_windows.tsv.gz
#   <snake_dir>/snake_regions.tsv.gz
#   <snake_dir>/snake_hierarchy.tsv.gz
#   <snake_dir>/snake_decision_log.tsv.gz
#   <snake_dir>/snake_window_states.tsv.gz
#   <snake_dir>/snake_diagnostics.tsv.gz
#   <snake_dir>/snake_multiscale_comparison.tsv.gz
#   <snake_dir>/snake_summary.tsv
#   <step10_outprefix>.mds.rds
#
# OPTIONAL (auto-detected):
#   <snake2_dir>/snake2_track.tsv.gz
#   <snake3_dir>/snake3_track.tsv.gz
#   <snake4_dir>/snake4_track.tsv.gz + snake4_regimes.tsv.gz
#   <step10_outprefix>.candidate_regions.tsv.gz  (500kb control)
#
# Usage:
#   Rscript STEP10f_v3_snake_diagnostics.R \
#     <step10_outprefix> <snake_dir> <outdir> [chrom=all] \
#     [--snake2_dir <path>] [--snake3_dir <path>] [--snake4_dir <path>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP10f_v3_snake_diagnostics.R ",
       "<step10_outprefix> <snake_dir> <outdir> [chrom=all] ",
       "[--snake2_dir <path>] [--snake3_dir <path>] [--snake4_dir <path>]")
}

step10_prefix <- args[1]
snake_dir     <- args[2]
outdir        <- args[3]
chrom_filter  <- NULL
snake2_dir    <- NULL
snake3_dir    <- NULL
snake4_dir    <- NULL

# Parse positional + named args
i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--snake2_dir" && i < length(args)) {
    snake2_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--snake3_dir" && i < length(args)) {
    snake3_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--snake4_dir" && i < length(args)) {
    snake4_dir <- args[i + 1]; i <- i + 2L
  } else if (!grepl("^--", a) && a != "all") {
    chrom_filter <- a; i <- i + 1L
  } else {
    i <- i + 1L
  }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[DIAG] Loading data...")

# =============================================================================
# LOAD ALL INPUTS
# =============================================================================

safe_load <- function(path) {
  if (file.exists(path)) tryCatch(fread(path), error = function(e) {
    message("[WARN] Failed to read: ", path); NULL
  }) else { message("[WARN] Missing: ", path); NULL }
}

win_dt   <- safe_load(file.path(snake_dir, "snake_windows.tsv.gz"))
reg_dt   <- safe_load(file.path(snake_dir, "snake_regions.tsv.gz"))
hier_dt  <- safe_load(file.path(snake_dir, "snake_hierarchy.tsv.gz"))
log_dt   <- safe_load(file.path(snake_dir, "snake_decision_log.tsv.gz"))
state_dt <- safe_load(file.path(snake_dir, "snake_window_states.tsv.gz"))
diag_dt  <- safe_load(file.path(snake_dir, "snake_diagnostics.tsv.gz"))
comp_dt  <- safe_load(file.path(snake_dir, "snake_multiscale_comparison.tsv.gz"))
summ_dt  <- safe_load(file.path(snake_dir, "snake_summary.tsv"))

mds_rds <- paste0(step10_prefix, ".mds.rds")
mds_obj <- if (file.exists(mds_rds)) readRDS(mds_rds) else NULL

# ── Optional: Snake 2/3/4 tracks (auto-detect from sibling dirs if not specified)
parent_dir <- dirname(snake_dir)
if (is.null(snake2_dir)) snake2_dir <- file.path(parent_dir, "snake2_community")
if (is.null(snake3_dir)) snake3_dir <- file.path(parent_dir, "snake3_ghsl")
if (is.null(snake4_dir)) snake4_dir <- file.path(parent_dir, "snake4_threeband")

s2_track_dt <- safe_load(file.path(snake2_dir, "snake2_track.tsv.gz"))
s3_track_dt <- safe_load(file.path(snake3_dir, "snake3_track.tsv.gz"))
s4_track_dt <- safe_load(file.path(snake4_dir, "snake4_track.tsv.gz"))
s4_regime_dt <- safe_load(file.path(snake4_dir, "snake4_regimes.tsv.gz"))

# Inversion-likeness scores (from Snake 1 output)
inv_like_dt <- safe_load(file.path(snake_dir, "snake_inv_likeness.tsv.gz"))

# ── Optional: 500kb control candidates
control_file <- paste0(step10_prefix, ".candidate_regions.tsv.gz")
control_dt <- safe_load(control_file)

message("[DIAG] Optional tracks loaded:",
        if (!is.null(s2_track_dt)) paste0(" S2(", nrow(s2_track_dt), ")") else "",
        if (!is.null(s3_track_dt)) paste0(" S3(", nrow(s3_track_dt), ")") else "",
        if (!is.null(s4_track_dt)) paste0(" S4(", nrow(s4_track_dt), ")") else "",
        if (!is.null(s4_regime_dt)) paste0(" regimes(", nrow(s4_regime_dt), ")") else "",
        if (!is.null(control_dt)) paste0(" control(", nrow(control_dt), ")") else "")

if (!is.null(chrom_filter)) {
  for (nm in c("win_dt", "reg_dt", "log_dt", "state_dt", "diag_dt")) {
    dt <- get(nm)
    if (!is.null(dt) && "chrom" %in% names(dt)) assign(nm, dt[chrom == chrom_filter])
  }
}

chroms <- if (!is.null(state_dt)) sort(unique(state_dt$chrom)) else character(0)
message("[DIAG] Chromosomes: ", length(chroms))

# Colors
FAMILY_COLORS <- c(
  "1S" = "#1d4ed8", "1M" = "#16a34a", "1L" = "#d97706",
  "merge_A" = "#dc2626", "merge_B" = "#9333ea",
  "seed" = "#000000", "unused" = "#e5e7eb", "collected" = "#93c5fd"
)
STATUS_COLORS <- c(
  "accepted" = "#16a34a", "tolerated" = "#f59e0b", "seed" = "#1d4ed8",
  "bridge_accepted" = "#06b6d4", "bridge_tolerated" = "#f97316",
  "rejected" = "#ef4444", "stopped" = "#991b1b", "skipped" = "#d1d5db",
  "blocked" = "#7f1d1d"
)

# Snake 2/3/4 PASS/WEAK/FAIL/PARTIAL colors (shared)
SNAKE_STATUS_COLORS <- c(
  "PASS" = "#16a34a", "WEAK" = "#f59e0b", "FAIL" = "#e5e7eb",
  "PARTIAL" = "#f97316"
)

# Regime palette: up to 12 regimes per chromosome + unassigned
REGIME_PALETTE <- c(
  "1" = "#2563eb", "2" = "#16a34a", "3" = "#dc2626",
  "4" = "#9333ea", "5" = "#d97706", "6" = "#0891b2",
  "7" = "#be185d", "8" = "#4f46e5", "9" = "#059669",
  "10" = "#b91c1c", "11" = "#7c3aed", "12" = "#c2410c",
  "0" = "#d1d5db"  # unassigned
)

# 500kb control color
CONTROL_COLOR <- "#a78bfa"

# =============================================================================
# LEVEL 1A — CHROMOSOME LAYERED IDEOGRAM TRACKS (MOST IMPORTANT)
# =============================================================================

build_track_plot <- function(chr) {
  chr_state <- state_dt[chrom == chr]
  chr_win   <- win_dt[chrom == chr]
  chr_reg   <- reg_dt[chrom == chr]
  if (nrow(chr_state) == 0) return(NULL)

  # Build track data: each row = one window × one track
  tracks <- list()

  # Track 1: raw windows (all)
  tracks[[1]] <- chr_state[, .(
    pos = (start_bp + end_bp) / 2, track = "01_all_windows",
    status = collector_state, family = NA_character_
  )]

  # Track 2: seeds
  seeds <- chr_state[collector_state == "seed" | grepl("seed", collector_state, ignore.case = TRUE)]
  if (is.null(seeds) || nrow(seeds) == 0) {
    seeds <- chr_state[max_abs_z >= 2.0]
    seeds[, collector_state := fifelse(max_abs_z >= 3.0, "strong_seed",
                                        fifelse(max_abs_z >= 2.5, "mod_seed", "weak_seed"))]
  }
  if (nrow(seeds) > 0) {
    tracks[[2]] <- seeds[, .(pos = (start_bp + end_bp) / 2, track = "02_seeds",
                              status = collector_state, family = NA_character_)]
  }

  # Tracks 3-5: S1-small, S1-medium, S1-large
  for (fam in c("1S", "1M", "1L")) {
    fam_win <- chr_win[core_family == fam & snake_phase == "core"]
    if (nrow(fam_win) > 0) {
      tnum <- switch(fam, "1S" = "03", "1M" = "04", "1L" = "05")
      tracks[[length(tracks) + 1]] <- fam_win[, .(
        pos = (start_bp + end_bp) / 2,
        track = paste0(tnum, "_S1_", fam),
        status = inclusion_status, family = fam
      )]
    }
  }

  # Tracks 6-7: Merge A, Merge B
  for (mfam in c("merge_A", "merge_B")) {
    mfam_win <- chr_win[merge_family == mfam & snake_phase == "merged"]
    if (nrow(mfam_win) > 0) {
      tnum <- switch(mfam, "merge_A" = "06", "merge_B" = "07")
      tracks[[length(tracks) + 1]] <- mfam_win[, .(
        pos = (start_bp + end_bp) / 2,
        track = paste0(tnum, "_", mfam),
        status = inclusion_status, family = NA_character_
      )]
    }
  }

  # Track 8: QC flags
  if (!is.null(diag_dt) && nrow(diag_dt[chrom == chr]) > 0) {
    chr_diag <- diag_dt[chrom == chr]
    # Map QC flags to regions
    flagged_regs <- chr_diag[qc_flags != "clean"]
    if (nrow(flagged_regs) > 0) {
      # Get the windows for flagged regions
      flagged_sids <- flagged_regs$snake_id
      qc_win <- chr_win[snake_id %in% flagged_sids]
      if (nrow(qc_win) > 0) {
        qc_flags_str <- chr_diag[match(qc_win$snake_id, chr_diag$snake_id)]$qc_flags
        tracks[[length(tracks) + 1]] <- data.table(
          pos = (qc_win$start_bp + qc_win$end_bp) / 2,
          track = "08_QC_flags",
          status = qc_flags_str,
          family = NA_character_
        )
      }
    }
  }

  if (length(tracks) == 0) return(NULL)
  track_dt <- rbindlist(tracks, fill = TRUE)

  # ── ADDITIONAL TRACKS: Snake 2, 3, 4, regime, control ──────────────
  # These are added separately because they come from optional inputs

  # Track 09: Snake 2 status
  if (!is.null(s2_track_dt) && nrow(s2_track_dt[chrom == chr]) > 0) {
    s2_chr <- s2_track_dt[chrom == chr]
    track_dt <- rbind(track_dt, data.table(
      pos = (s2_chr$start_bp + s2_chr$end_bp) / 2,
      track = "09_Snake2", status = s2_chr$snake2_status,
      family = NA_character_
    ), fill = TRUE)
  }

  # Track 09b: Snake 2 profile concordance (diagnostic layer)
  if (!is.null(s2_track_dt) && "profile_concordance" %in% names(s2_track_dt)) {
    s2_chr <- s2_track_dt[chrom == chr]
    s2_prof <- s2_chr[is.finite(profile_concordance)]
    if (nrow(s2_prof) > 0) {
      prof_status <- fifelse(s2_prof$profile_concordance >= 0.5, "PASS",
                     fifelse(s2_prof$profile_concordance >= 0.2, "WEAK", "FAIL"))
      track_dt <- rbind(track_dt, data.table(
        pos = (s2_prof$start_bp + s2_prof$end_bp) / 2,
        track = "09b_Profile", status = prof_status,
        family = NA_character_
      ), fill = TRUE)
    }
  }

  # Track 09c: Inversion-likeness (from local PCA eigenvalues, no MDS)
  if (!is.null(inv_like_dt) && nrow(inv_like_dt[chrom == chr]) > 0) {
    il_chr <- merge(inv_like_dt[chrom == chr],
                    as.data.table(per_chr[[chr]]$out_dt)[, .(global_window_id, start_bp, end_bp)],
                    by = "global_window_id")
    il_chr <- il_chr[is.finite(inv_likeness)]
    if (nrow(il_chr) > 0) {
      il_status <- fifelse(il_chr$inv_likeness >= 0.5, "PASS",
                   fifelse(il_chr$inv_likeness >= 0.25, "WEAK", "FAIL"))
      track_dt <- rbind(track_dt, data.table(
        pos = (il_chr$start_bp + il_chr$end_bp) / 2,
        track = "09c_InvLike", status = il_status,
        family = NA_character_
      ), fill = TRUE)
    }
  }

  # Track 10: Snake 3 status
  if (!is.null(s3_track_dt) && nrow(s3_track_dt[chrom == chr]) > 0) {
    s3_chr <- s3_track_dt[chrom == chr]
    track_dt <- rbind(track_dt, data.table(
      pos = (s3_chr$start_bp + s3_chr$end_bp) / 2,
      track = "10_Snake3", status = s3_chr$snake3_status,
      family = NA_character_
    ), fill = TRUE)
  }

  # Track 11: Snake 4 status
  if (!is.null(s4_track_dt) && nrow(s4_track_dt[chrom == chr]) > 0) {
    s4_chr <- s4_track_dt[chrom == chr]
    track_dt <- rbind(track_dt, data.table(
      pos = (s4_chr$start_bp + s4_chr$end_bp) / 2,
      track = "11_Snake4", status = s4_chr$snake4_status,
      family = NA_character_
    ), fill = TRUE)
  }

  # Track 12: Snake 4 regime identity
  if (!is.null(s4_track_dt) && "middle_regime_id" %in% names(s4_track_dt)) {
    s4_chr <- s4_track_dt[chrom == chr & !is.na(middle_regime_id)]
    if (nrow(s4_chr) > 0) {
      track_dt <- rbind(track_dt, data.table(
        pos = (s4_chr$start_bp + s4_chr$end_bp) / 2,
        track = "12_S4_regime",
        status = paste0("regime_", s4_chr$middle_regime_id),
        family = NA_character_
      ), fill = TRUE)
    }
  }

  # Track 13: 500kb control candidates (shaded regions)
  if (!is.null(control_dt) && nrow(control_dt[chrom == chr]) > 0) {
    ctrl_chr <- control_dt[chrom == chr]
    ctrl_rows <- list()
    for (ci in seq_len(nrow(ctrl_chr))) {
      cr <- ctrl_chr[ci]
      # Create a few points spanning the candidate region
      ctrl_rows[[ci]] <- data.table(
        pos = seq(cr$start_bp, cr$end_bp, length.out = max(3, cr$n_windows %||% 3)) / 1e6 * 1e6,
        track = "13_CONTROL_500kb",
        status = paste0("ctrl_", cr$candidate_id),
        family = NA_character_
      )
    }
    if (length(ctrl_rows) > 0) {
      track_dt <- rbind(track_dt, rbindlist(ctrl_rows), fill = TRUE)
    }
  }

  # Plot
  p <- ggplot(track_dt, aes(x = pos / 1e6, y = track, fill = status)) +
    geom_tile(height = 0.7, width = diff(range(chr_state$start_bp)) / nrow(chr_state) / 1e6 * 0.9) +
    scale_fill_manual(values = c(STATUS_COLORS, FAMILY_COLORS, SNAKE_STATUS_COLORS,
                                  "collected" = "#93c5fd", "unused" = "#f3f4f6",
                                  "seed" = "#1d4ed8",
                                  "strong_seed" = "#1e3a5f", "mod_seed" = "#3b82f6",
                                  "weak_seed" = "#93c5fd",
                                  "clean" = "#d1fae5",
                                  # Regime colors (regime_1 through regime_12 + regime_0)
                                  setNames(REGIME_PALETTE, paste0("regime_", names(REGIME_PALETTE))),
                                  # Control candidate colors (auto-generate)
                                  setNames(rep(CONTROL_COLOR, 100),
                                           paste0("ctrl_", seq_len(100)))),
                      na.value = "#f3f4f6", name = "Status") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " — Multi-Scale Snake Ideogram"),
         subtitle = paste0(nrow(chr_state), " windows, ",
                           sum(chr_state$collector_state != "unused"), " collected")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_text(size = 7, hjust = 1),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          panel.grid.major.y = element_blank())

  # Add region boundary annotations
  if (nrow(chr_reg) > 0) {
    for (ri in seq_len(min(50, nrow(chr_reg)))) {
      p <- p + annotate("segment",
        x = chr_reg$start_bp[ri] / 1e6, xend = chr_reg$end_bp[ri] / 1e6,
        y = -0.2, yend = -0.2,
        color = FAMILY_COLORS[chr_reg$core_family[ri] %||% chr_reg$merge_family[ri]] %||% "#333",
        linewidth = 1.5, alpha = 0.6)
    }
  }

  p
}

# =============================================================================
# LEVEL 1B — MDS SCATTER COLORED BY SNAKE STATUS
# =============================================================================

build_mds_scatter <- function(chr) {
  if (is.null(mds_obj) || is.null(mds_obj$per_chr[[chr]])) return(NULL)

  dt <- as.data.table(mds_obj$per_chr[[chr]]$out_dt)
  if (!("MDS1" %in% names(dt)) || !("MDS2" %in% names(dt))) return(NULL)

  # Merge snake membership
  chr_win <- win_dt[chrom == chr]
  chr_state <- state_dt[chrom == chr]

  dt[, snake_status := "background"]
  if (nrow(chr_state) > 0) {
    dt[global_window_id %in% chr_state[collector_state == "seed"]$global_window_id,
       snake_status := "seed"]
  }
  if (nrow(chr_win) > 0) {
    for (fam in c("1S", "1M", "1L")) {
      fam_wids <- chr_win[core_family == fam]$global_window_id
      dt[global_window_id %in% fam_wids, snake_status := paste0("S1_", fam)]
    }
    bridge_wids <- chr_win[grepl("bridge", inclusion_status)]$global_window_id
    dt[global_window_id %in% bridge_wids, snake_status := "bridge"]
  }

  mds_colors <- c("background" = "#e5e7eb", "seed" = "#000000",
                   "S1_1S" = "#1d4ed8", "S1_1M" = "#16a34a", "S1_1L" = "#d97706",
                   "bridge" = "#06b6d4")

  p <- ggplot(dt, aes(x = MDS1, y = MDS2, color = snake_status)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = mds_colors, name = "Status") +
    labs(title = paste0(chr, " — MDS Scatter by Snake Membership"),
         subtitle = paste0(nrow(dt), " windows")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")

  p
}

# =============================================================================
# LEVEL 1B2 — MDS SCATTER COLORED BY MERGED REGION ID
# =============================================================================

build_mds_scatter_by_merged <- function(chr) {
  if (is.null(mds_obj) || is.null(mds_obj$per_chr[[chr]])) return(NULL)

  dt <- as.data.table(mds_obj$per_chr[[chr]]$out_dt)
  if (!("MDS1" %in% names(dt)) || !("MDS2" %in% names(dt))) return(NULL)

  chr_win <- win_dt[chrom == chr & snake_phase == "merged"]
  if (is.null(chr_win) || nrow(chr_win) == 0) return(NULL)

  dt[, merged_id := "background"]

  # Color by merge_A snake_id (each merged region gets a unique color)
  merge_a_win <- chr_win[merge_family == "merge_A"]
  if (nrow(merge_a_win) > 0) {
    for (sid in unique(merge_a_win$snake_id)) {
      wids <- merge_a_win[snake_id == sid]$global_window_id
      dt[global_window_id %in% wids, merged_id := paste0("M", sid)]
    }
  }

  # Generate colors for merged IDs
  merged_ids <- setdiff(unique(dt$merged_id), "background")
  n_merged <- length(merged_ids)
  if (n_merged == 0) return(NULL)

  merged_pal <- c("background" = "#e5e7eb")
  if (n_merged <= 12) {
    hues <- rainbow(n_merged, s = 0.7, v = 0.8)
  } else {
    hues <- rainbow(n_merged, s = 0.65, v = 0.75)
  }
  for (mi in seq_along(merged_ids)) {
    merged_pal[merged_ids[mi]] <- hues[mi]
  }

  p <- ggplot(dt, aes(x = MDS1, y = MDS2, color = merged_id)) +
    geom_point(data = dt[merged_id == "background"], size = 1, alpha = 0.3) +
    geom_point(data = dt[merged_id != "background"], size = 2, alpha = 0.8) +
    scale_color_manual(values = merged_pal, name = "Merged ID") +
    labs(title = paste0(chr, " — MDS Scatter by Merge A Region ID"),
         subtitle = paste0(n_merged, " merged regions, ", nrow(dt), " windows")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")

  p
}

# =============================================================================
# LEVEL 1B3 — MDS SCATTER COLORED BY SNAKE 4 REGIME
# =============================================================================

build_mds_scatter_by_regime <- function(chr) {
  if (is.null(mds_obj) || is.null(mds_obj$per_chr[[chr]])) return(NULL)
  if (is.null(s4_track_dt) || !("middle_regime_id" %in% names(s4_track_dt))) return(NULL)

  dt <- as.data.table(mds_obj$per_chr[[chr]]$out_dt)
  if (!("MDS1" %in% names(dt)) || !("MDS2" %in% names(dt))) return(NULL)

  s4_chr <- s4_track_dt[chrom == chr & !is.na(middle_regime_id) & middle_regime_id > 0]
  if (nrow(s4_chr) == 0) return(NULL)

  dt[, regime := "unassigned"]
  for (ri in unique(s4_chr$middle_regime_id)) {
    wids <- s4_chr[middle_regime_id == ri]$global_window_id
    dt[global_window_id %in% wids, regime := paste0("regime_", ri)]
  }

  regime_ids <- setdiff(unique(dt$regime), "unassigned")
  n_reg <- length(regime_ids)
  if (n_reg == 0) return(NULL)

  pal <- c("unassigned" = "#e5e7eb")
  for (ri in seq_along(regime_ids)) {
    rid <- sub("regime_", "", regime_ids[ri])
    pal[regime_ids[ri]] <- REGIME_PALETTE[rid] %||% rainbow(1, s = 0.7)
  }

  p <- ggplot(dt, aes(x = MDS1, y = MDS2, color = regime)) +
    geom_point(data = dt[regime == "unassigned"], size = 1, alpha = 0.2) +
    geom_point(data = dt[regime != "unassigned"], size = 2, alpha = 0.8) +
    scale_color_manual(values = pal, name = "S4 Regime") +
    labs(title = paste0(chr, " — MDS Scatter by Snake 4 Middle-Core Regime"),
         subtitle = paste0(n_reg, " regimes, ", sum(dt$regime != "unassigned"),
                           " windows assigned")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")

  p
}

# =============================================================================
# LEVEL 1C — CONTINUITY SCORE PROFILES
# =============================================================================

build_continuity_profile <- function(chr) {
  chr_log <- log_dt[chrom == chr & phase == "core"]
  if (is.null(chr_log) || nrow(chr_log) == 0) return(NULL)

  # Get window positions for x-axis
  chr_state <- state_dt[chrom == chr]
  pos_map <- setNames(chr_state$start_bp, as.character(chr_state$global_window_id))

  chr_log[, to_pos := pos_map[as.character(to_window_id)] / 1e6]
  chr_log <- chr_log[!is.na(to_pos)]

  if (nrow(chr_log) == 0) return(NULL)

  p <- ggplot(chr_log, aes(x = to_pos, y = continuity_score, color = decision)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(aes(yintercept = accept_thresh), linetype = "dashed", color = "#16a34a", linewidth = 0.3) +
    geom_hline(aes(yintercept = toler_thresh), linetype = "dotted", color = "#f59e0b", linewidth = 0.3) +
    facet_wrap(~ core_family, ncol = 1, scales = "free_y") +
    scale_color_manual(values = STATUS_COLORS, name = "Decision") +
    labs(x = paste0(chr, " (Mb)"), y = "Continuity Score",
         title = paste0(chr, " — Continuity Score Profiles by Core Family"),
         subtitle = "Dashed=accept, dotted=tolerate threshold") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

  p
}

# =============================================================================
# LEVEL 2E — NESTED BOUNDARY LADDERS
# =============================================================================

build_nested_ladders <- function(chr) {
  chr_reg <- reg_dt[chrom == chr]
  if (is.null(chr_reg) || nrow(chr_reg) == 0) return(NULL)

  # Build ladder: one row per region, y = family track
  family_order <- c("1S", "1M", "1L", "merge_A", "merge_B")
  chr_reg[, family_label := fifelse(!is.na(core_family), core_family, merge_family)]
  chr_reg <- chr_reg[family_label %in% family_order]
  chr_reg[, y_pos := match(family_label, family_order)]

  if (nrow(chr_reg) == 0) return(NULL)

  p <- ggplot(chr_reg, aes(xmin = start_bp / 1e6, xmax = end_bp / 1e6,
                             ymin = y_pos - 0.35, ymax = y_pos + 0.35,
                             fill = family_label)) +
    geom_rect(alpha = 0.7, color = "#333", linewidth = 0.2) +
    scale_fill_manual(values = FAMILY_COLORS, name = "Family") +
    scale_y_continuous(breaks = seq_along(family_order), labels = family_order) +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " — Nested Boundary Ladders"),
         subtitle = "Core S/M/L → Merge A/B interval nesting") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  p
}

# =============================================================================
# LEVEL 3 — MULTI-HEAD CANDIDATE SUMMARY TABLE (THE TABLE)
# =============================================================================

build_candidate_summary_table <- function() {
  if (is.null(reg_dt) || nrow(reg_dt) == 0) return(data.table())

  # Group by seed/core clusters: find overlapping regions across families
  # Use hierarchy to map cores → merges
  cores <- reg_dt[snake_phase == "core"]
  merges <- reg_dt[snake_phase == "merged"]

  if (nrow(cores) == 0) return(data.table())

  # Build one row per core (or per seed-core group)
  rows <- list()
  for (ci in seq_len(nrow(cores))) {
    cr <- cores[ci]
    sid <- cr$snake_id
    chr <- cr$chrom

    # Block 1: seed stats
    seed_windows <- win_dt[snake_id == sid & inclusion_status == "seed"]
    seed_z <- if (!is.null(state_dt)) {
      state_dt[global_window_id %in% seed_windows$global_window_id]$max_abs_z
    } else NA_real_

    # Block 2-4: per core family
    fam <- cr$core_family

    # Block 5-6: merge info from hierarchy
    if (!is.null(hier_dt) && nrow(hier_dt) > 0) {
      merge_a_links <- hier_dt[core_snake_id == sid & merge_family == "merge_A"]
      merge_b_links <- hier_dt[core_snake_id == sid & merge_family == "merge_B"]
    } else {
      merge_a_links <- data.table()
      merge_b_links <- data.table()
    }

    merge_a_reg <- if (nrow(merge_a_links) > 0) {
      merges[snake_id %in% merge_a_links$merged_snake_id]
    } else data.table()

    merge_b_reg <- if (nrow(merge_b_links) > 0) {
      merges[snake_id %in% merge_b_links$merged_snake_id]
    } else data.table()

    # Block 7: QC
    qc_info <- if (!is.null(diag_dt) && nrow(diag_dt) > 0) {
      # Find QC for merge regions containing this core
      all_merge_sids <- c(merge_a_links$merged_snake_id, merge_b_links$merged_snake_id)
      diag_dt[snake_id %in% all_merge_sids]
    } else data.table()

    # Block 8: stop reason from decision log
    stop_entries <- if (!is.null(log_dt)) {
      log_dt[snake_id == sid & decision %in% c("stopped", "rejected")]
    } else data.table()
    primary_stop <- if (nrow(stop_entries) > 0) stop_entries$reason[1] else "unknown"

    rows[[ci]] <- data.table(
      # Block 1: seed
      chrom = chr,
      core_snake_id = sid,
      core_family = fam,
      core_start_bp = cr$start_bp,
      core_end_bp = cr$end_bp,
      core_n_windows = cr$n_windows,
      seed_max_z = round(max(seed_z, na.rm = TRUE), 2),
      n_seeds = cr$n_seeds %||% 0L,

      # Block 2: core stats
      core_mean_score = cr$mean_score,
      core_min_score = cr$min_score,
      core_coherence = cr$coherence,
      core_n_tolerated = cr$n_tolerated %||% 0L,
      core_stop_reason = primary_stop,

      # Block 5: merge A
      merge_A_merged = nrow(merge_a_reg) > 0,
      merge_A_n_windows = if (nrow(merge_a_reg) > 0) sum(merge_a_reg$n_windows) else 0L,
      merge_A_coherence = if (nrow(merge_a_reg) > 0) round(mean(merge_a_reg$coherence), 4) else NA_real_,
      merge_A_start = if (nrow(merge_a_reg) > 0) min(merge_a_reg$start_bp) else NA_real_,
      merge_A_end = if (nrow(merge_a_reg) > 0) max(merge_a_reg$end_bp) else NA_real_,

      # Block 6: merge B
      merge_B_merged = nrow(merge_b_reg) > 0,
      merge_B_n_windows = if (nrow(merge_b_reg) > 0) sum(merge_b_reg$n_windows) else 0L,
      merge_B_coherence = if (nrow(merge_b_reg) > 0) round(mean(merge_b_reg$coherence), 4) else NA_real_,
      merge_B_start = if (nrow(merge_b_reg) > 0) min(merge_b_reg$start_bp) else NA_real_,
      merge_B_end = if (nrow(merge_b_reg) > 0) max(merge_b_reg$end_bp) else NA_real_,

      # Block 7: QC
      qc_flags = if (nrow(qc_info) > 0) paste(unique(qc_info$qc_flags), collapse = ";") else "no_merge_qc",

      # Block 8: stability
      merge_agreement = if (nrow(merge_a_reg) > 0 && nrow(merge_b_reg) > 0) "both"
                        else if (nrow(merge_b_reg) > 0) "B_only_sensitive"
                        else if (nrow(merge_a_reg) > 0) "A_only"
                        else "no_merge",
      boundary_spread_bp = if (nrow(merge_b_reg) > 0 && nrow(merge_a_reg) > 0)
        max(merge_b_reg$end_bp) - min(merge_b_reg$start_bp) -
        (max(merge_a_reg$end_bp) - min(merge_a_reg$start_bp))
        else NA_real_
    )
  }

  rbindlist(rows, fill = TRUE)
}

# =============================================================================
# LEVEL 3B — SNAKE CENSUS TABLE
# =============================================================================

build_census_table <- function() {
  if (is.null(summ_dt) || nrow(summ_dt) == 0) {
    if (is.null(state_dt) || is.null(reg_dt)) return(data.table())
  }

  # Use existing summary if available
  if (!is.null(summ_dt) && nrow(summ_dt) > 0) {
    census <- copy(summ_dt)
  } else {
    census <- state_dt[, .(n_windows = .N, n_collected = sum(collector_state != "unused")),
                        by = chrom]
  }

  # Add merge counts from regions
  if (!is.null(reg_dt) && nrow(reg_dt) > 0) {
    merge_counts <- reg_dt[snake_phase == "merged", .(
      n_merge_A = sum(merge_family == "merge_A", na.rm = TRUE),
      n_merge_B = sum(merge_family == "merge_B", na.rm = TRUE)
    ), by = chrom]
    census <- merge(census, merge_counts, by = "chrom", all.x = TRUE)
  }

  # Add QC counts
  if (!is.null(diag_dt) && nrow(diag_dt) > 0) {
    qc_counts <- diag_dt[, .(
      n_qc_flagged = sum(qc_flags != "clean"),
      n_split_suggested = sum(grepl("split", qc_flags)),
      n_overmerge = sum(grepl("overmerge", qc_flags))
    ), by = chrom]
    census <- merge(census, qc_counts, by = "chrom", all.x = TRUE)
  }

  census
}

# =============================================================================
# LEVEL 3C — STOP REASON TABLE
# =============================================================================

build_stop_reason_table <- function() {
  if (is.null(log_dt) || nrow(log_dt) == 0) return(data.table())

  # Count stop/reject reasons by family
  stops <- log_dt[decision %in% c("stopped", "rejected", "blocked")]
  if (nrow(stops) == 0) return(data.table())

  reason_tab <- stops[, .N, by = .(core_family, reason)]
  setorder(reason_tab, core_family, -N)
  reason_tab
}

build_stop_reason_plot <- function(reason_tab) {
  if (is.null(reason_tab) || nrow(reason_tab) == 0) return(NULL)

  p <- ggplot(reason_tab, aes(x = core_family, y = N, fill = reason)) +
    geom_col(position = "stack", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", name = "Stop Reason") +
    labs(x = "Snake Family", y = "Count",
         title = "Stop Reasons by Snake Family",
         subtitle = "Why each snake stopped extending") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")

  p
}

# =============================================================================
# GENERATE ALL OUTPUTS
# =============================================================================

message("[DIAG] Generating Level 1 — Window-level plots...")

# 01: Chromosome track ideograms
track_plots <- list()
for (chr in chroms) {
  p <- build_track_plot(chr)
  if (!is.null(p)) track_plots[[chr]] <- p
}
if (length(track_plots) > 0) {
  f01 <- file.path(outdir, "01_chromosome_track_ideogram.pdf")
  pdf(f01, width = 14, height = 6)
  for (p in track_plots) print(p)
  dev.off()
  message("  ", f01, " (", length(track_plots), " pages)")
}

# 02: MDS scatter
mds_plots <- list()
for (chr in chroms) {
  p <- build_mds_scatter(chr)
  if (!is.null(p)) mds_plots[[chr]] <- p
}
if (length(mds_plots) > 0) {
  f02 <- file.path(outdir, "02_mds_scatter_by_snake.pdf")
  pdf(f02, width = 8, height = 7)
  for (p in mds_plots) print(p)
  dev.off()
  message("  ", f02)
}

# 02b: MDS scatter by merged region ID
mds_merged_plots <- list()
for (chr in chroms) {
  p <- build_mds_scatter_by_merged(chr)
  if (!is.null(p)) mds_merged_plots[[chr]] <- p
}
if (length(mds_merged_plots) > 0) {
  f02b <- file.path(outdir, "02b_mds_scatter_by_merged_id.pdf")
  pdf(f02b, width = 8, height = 7)
  for (p in mds_merged_plots) print(p)
  dev.off()
  message("  ", f02b)
}

# 02c: MDS scatter by Snake 4 regime
mds_regime_plots <- list()
for (chr in chroms) {
  p <- build_mds_scatter_by_regime(chr)
  if (!is.null(p)) mds_regime_plots[[chr]] <- p
}
if (length(mds_regime_plots) > 0) {
  f02c <- file.path(outdir, "02c_mds_scatter_by_regime.pdf")
  pdf(f02c, width = 8, height = 7)
  for (p in mds_regime_plots) print(p)
  dev.off()
  message("  ", f02c)
}

# 03: Continuity score profiles
cont_plots <- list()
for (chr in chroms) {
  p <- build_continuity_profile(chr)
  if (!is.null(p)) cont_plots[[chr]] <- p
}
if (length(cont_plots) > 0) {
  f03 <- file.path(outdir, "03_continuity_score_profiles.pdf")
  pdf(f03, width = 12, height = 8)
  for (p in cont_plots) print(p)
  dev.off()
  message("  ", f03)
}

message("[DIAG] Generating Level 2 — Region-level plots...")

# 05: Nested boundary ladders
ladder_plots <- list()
for (chr in chroms) {
  p <- build_nested_ladders(chr)
  if (!is.null(p)) ladder_plots[[chr]] <- p
}
if (length(ladder_plots) > 0) {
  f05 <- file.path(outdir, "05_nested_boundary_ladders.pdf")
  pdf(f05, width = 12, height = 5)
  for (p in ladder_plots) print(p)
  dev.off()
  message("  ", f05)
}

message("[DIAG] Generating Level 3 — Summary tables...")

# 08: Candidate summary table
cand_table <- build_candidate_summary_table()
if (nrow(cand_table) > 0) {
  f08 <- file.path(outdir, "08_candidate_summary_table.tsv.gz")
  fwrite(cand_table, f08, sep = "\t")
  message("  ", f08, " (", nrow(cand_table), " rows)")
}

# 09: Census table
census <- build_census_table()
if (nrow(census) > 0) {
  f09 <- file.path(outdir, "09_chromosome_census.tsv")
  fwrite(census, f09, sep = "\t")
  message("  ", f09)
}

# 10: Stop reason table + plot
reason_tab <- build_stop_reason_table()
if (nrow(reason_tab) > 0) {
  f10t <- file.path(outdir, "10_stop_reason_table.tsv")
  fwrite(reason_tab, f10t, sep = "\t")
  message("  ", f10t)

  reason_plot <- build_stop_reason_plot(reason_tab)
  if (!is.null(reason_plot)) {
    f10p <- file.path(outdir, "10_stop_reason_barplot.pdf")
    ggsave(f10p, reason_plot, width = 8, height = 5)
    message("  ", f10p)
  }
}

# =============================================================================
# 11: Profile diversity plot (per-chromosome, from Snake 2 profile layer)
#
# Four stacked curves per chromosome:
#   Panel A — tailA within-stripe Hamming (higher = more diverse = sub-clusters?)
#   Panel B — tailB within-stripe Hamming
#   Panel C — middle within-stripe Hamming
#   Panel D — profile concordance between adjacent windows (lower = structure changing)
#
# Windows where Snake 1 detected signal are highlighted with background shading.
# This plot answers: "where along the chromosome do same-stripe samples diverge?"
# =============================================================================

if (!is.null(s2_track_dt) && "profile_concordance" %in% names(s2_track_dt) &&
    any(is.finite(s2_track_dt$profile_concordance))) {

  message("\n[DIAG] Building profile diversity plots...")

  for (chr in chroms) {
    s2_chr <- s2_track_dt[chrom == chr]
    if (nrow(s2_chr) == 0 || !any(is.finite(s2_chr$profile_concordance))) next

    s2_chr[, pos_mb := (start_bp + end_bp) / 2e6]

    # Reshape to long for faceting
    prof_long <- rbindlist(list(
      s2_chr[, .(pos_mb, value = profile_tailA_hamming, panel = "A: tailA Hamming")],
      s2_chr[, .(pos_mb, value = profile_tailB_hamming, panel = "B: tailB Hamming")],
      s2_chr[, .(pos_mb, value = profile_middle_hamming, panel = "C: middle Hamming")],
      s2_chr[, .(pos_mb, value = profile_concordance, panel = "D: concordance")]
    ))
    prof_long <- prof_long[is.finite(value)]
    if (nrow(prof_long) == 0) next

    # Snake 1 regions for background shading
    snake1_bg <- data.table()
    if (!is.null(win_dt) && nrow(win_dt[chrom == chr]) > 0) {
      s1_chr <- win_dt[chrom == chr]
      if (nrow(s1_chr) > 0) {
        s1_regs <- s1_chr[, .(xmin = min(start_bp) / 1e6, xmax = max(end_bp) / 1e6),
                          by = snake_id]
        # Replicate across all panels
        snake1_bg <- rbindlist(lapply(unique(prof_long$panel), function(p) {
          cbind(s1_regs, panel = p)
        }))
      }
    }

    panel_order <- c("A: tailA Hamming", "B: tailB Hamming",
                     "C: middle Hamming", "D: concordance")
    prof_long[, panel := factor(panel, levels = panel_order)]
    if (nrow(snake1_bg) > 0) snake1_bg[, panel := factor(panel, levels = panel_order)]

    g <- ggplot(prof_long, aes(x = pos_mb, y = value))

    if (nrow(snake1_bg) > 0) {
      g <- g + geom_rect(data = snake1_bg, inherit.aes = FALSE,
                          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                          fill = "#3b82f6", alpha = 0.10)
    }

    g <- g +
      geom_line(color = "#334155", linewidth = 0.4) +
      geom_point(aes(color = value), size = 0.6) +
      scale_color_gradient2(low = "#2563eb", mid = "#f59e0b", high = "#dc2626",
                            midpoint = 0.3, guide = "none") +
      facet_wrap(~ panel, ncol = 1, scales = "free_y") +
      labs(x = paste0(chr, " position (Mb)"),
           y = "Score",
           title = paste0(chr, " — within-stripe profile diversity"),
           subtitle = "Blue shading = Snake 1 regions. High Hamming = sub-haplotype diversity. Low concordance = structure shift.") +
      theme_minimal(base_size = 9) +
      theme(strip.text = element_text(face = "bold", hjust = 0),
            panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = 7, color = "#64748b"))

    f11 <- file.path(outdir, paste0("11_profile_diversity_", chr, ".pdf"))
    ggsave(f11, g, width = 12, height = 7)
    message("  ", f11)
  }
}

message("\n[DONE] STEP_C02 snake diagnostics complete")
message("  Output: ", outdir)
