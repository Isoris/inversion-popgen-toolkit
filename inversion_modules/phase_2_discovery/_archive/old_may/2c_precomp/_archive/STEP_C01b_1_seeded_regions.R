#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01b_1_seeded_regions.R
#
# Phase 2 / 2c — seeded region-growing from MDS z-outlier seeds.
# Performs seed selection + extension per chromosome, writes per-chr regions.
# Consumed directly by phase_4/4a/STEP_C01d via --cores_dir.
#
# Codebase:    inversion_modules v8.5 / script v9.3.2
# Upstream:    STEP_C01a output — <precomp_dir>/precomp/<chr>.precomp.rds
#              <precomp_dir>/window_dt.tsv.gz
#              (optional) <pruned_samples>.ind for test_26 kin-pruned retention
# Downstream:  phase_4/4a/STEP_C01d_candidate_scoring --cores_dir
#              (seeded regions enter C01d's D2/D5 scoring dimensions at
#              catalog-birth time — no intermediate merge step)
# Formerly:    STEP_C01b_1_cores_wired_registry.R
#              The pipeline's "snake" / "core" terminology has been retired
#              in favour of "seeded region-growing" / "seeded region".
#              See RENAMING.md for the full map. The retired 1D fuzzy merge
#              (STEP_C01b_2_merge.R) is at
#              ../../_archive_superseded/fuzzy_merge_abandoned/.
#
# Algorithm
# ---------
# For each chromosome, three parameter "scale tiers" (1S small / 1M medium /
# 1L large) run independently to capture different-width candidates:
#
#   1. SEED SELECTION. A window is a seed candidate iff:
#        |max_abs_z|     >= seed_z_min     (MDS z-score outlier — primary)
#          OR
#        inv_likeness    >= 0.90           (very high inv-like signal)
#          OR
#        adaptive_seed   == TRUE           (test_07 Beta-distribution gate)
#      AND seed_nn_dist is finite and < seed_nn_max.
#      Z-score is the primary eligibility criterion. inv_likeness and the
#      Beta adaptive flag act as OR-gates for cases where z alone misses
#      a real signal (e.g. compact inversions with modest per-axis z).
#
#   2. SEED ORDERING. Seeds are processed highest-priority first.
#      priority = max_abs_z + 5 * inv_likeness_boost + 2 * spiky_inv_score
#      The spiky boost gives compact candidates priority so they are
#      seeded before larger extensions can swallow them.
#
#   3. EXTENSION (seeded region-growing). From each seed, extend
#      bidirectionally (right then left) accepting adjacent windows whose
#      continuity score exceeds a per-tier accept threshold. Windows below
#      the tolerate threshold accumulate damage (a scalar budget). The
#      extension halts when damage exceeds DMG_MAX or when a sharp drop
#      (> 2 * sd below rolling mean) is detected.
#
#      Morphology soft modifiers (v9.3.2):
#        flat_inv_score > 0.3         reduces damage (coherent plateau)
#        fragmentation_score > 0.3    increases damage (messy territory)
#
# Optional: test_26 kin-pruned retention
# --------------------------------------
# After all regions are found, test each one's sim_mat block contrast on a
# kin-pruned subset of samples. Regions whose signal collapses under kin
# pruning are flagged (NOT dropped — the hatchery-founder caveat means
# collapse under pruning is itself informative, not necessarily
# disqualifying).
#
# Per-chromosome execution
# ------------------------
# The script can process ONE chromosome (--chr or --chr_idx, for SLURM
# array parallelism) or ALL chromosomes serially (default, for rapid
# iteration).
#
# Outputs (per chromosome, in <outdir>/)
# --------------------------------------
#   seeded_regions_<chr>.rds                   list of regions (R object)
#   seeded_regions_windows_<chr>.tsv.gz        per-window region membership
#   seeded_regions_summary_<chr>.tsv.gz        per-region summary
#   seeded_regions_decision_log_<chr>.tsv.gz   per-step accept/tolerate/halt log
#   seeded_regions_window_states_<chr>.tsv.gz  per-window status across tiers
#   seeded_regions_summary_per_chr_<chr>.tsv   single-row chr summary
#
# Usage
# -----
#   Rscript STEP_C01b_1_seeded_regions.R <precomp_dir> <outdir>           # all chr
#   Rscript STEP_C01b_1_seeded_regions.R <precomp_dir> <outdir> --chr C_gar_LG01
#   Rscript STEP_C01b_1_seeded_regions.R <precomp_dir> <outdir> --chr_idx 1
#   Rscript STEP_C01b_1_seeded_regions.R <precomp_dir> <outdir> --tag relaxed
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ── Source load_bridge.R (provides reg for evidence registry) ──
.bridge_file <- Sys.getenv("LOAD_BRIDGE", "")
if (!nzchar(.bridge_file)) {
  for (.bp in c("utils/load_bridge.R", "../utils/load_bridge.R",
                 file.path(Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"),
                           "inversion_codebase_v8.5/utils/load_bridge.R"))) {
    if (file.exists(.bp)) { .bridge_file <- .bp; break }
  }
}
.bridge_available <- FALSE
if (nzchar(.bridge_file) && file.exists(.bridge_file)) {
  tryCatch({ source(.bridge_file); .bridge_available <- TRUE
    message("[seeded_regions] load_bridge.R sourced")
  }, error = function(e) message("[seeded_regions] load_bridge.R failed: ", conditionMessage(e)))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01b_1_cores.R <precomp_dir> <outdir> [--chr <name>] [--chr_idx <N>] [--tag <label>]")
}

precomp_dir <- args[1]
outdir      <- args[2]

run_tag   <- NULL
chr_name  <- NULL
chr_idx   <- NULL
pruned_file <- NULL
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--tag" && i < length(args))     { run_tag <- args[i+1]; i <- i + 2L }
  else if (a == "--chr" && i < length(args)) { chr_name <- args[i+1]; i <- i + 2L }
  else if (a == "--chr_idx" && i < length(args)) { chr_idx <- as.integer(args[i+1]); i <- i + 2L }
  else if (a == "--pruned_samples" && i < length(args)) { pruned_file <- args[i+1]; i <- i + 2L }
  else { i <- i + 1L }
}

if (!is.null(run_tag)) outdir <- file.path(outdir, run_tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# CORE FAMILY PARAMETERS
# =============================================================================

SEED_MDS_AXES    <- 5L
SEED_MIN_Z       <- 1.2
SEED_NEIGHBOR_K  <- 3L
SEED_MAX_NN_DIST <- 0.80

W_PREV <- 0.40; W_ROLL <- 0.40; W_MDS <- 0.20

S1S <- list(
  name = "1S", label = "strict",
  max_gap = 0L, max_bad = 0L, min_windows = 4L,
  accept = 0.60, tolerate = 0.35, roll_size = 3L,
  seed_z_min = 2.5, seed_nn_max = 0.70
)

S1M <- list(
  name = "1M", label = "moderate",
  max_gap = 1L, max_bad = 1L, min_windows = 3L,
  accept = 0.50, tolerate = 0.28, roll_size = 3L,
  seed_z_min = 1.8, seed_nn_max = 0.80
)

S1L <- list(
  name = "1L", label = "broad",
  max_gap = 2L, max_bad = 2L, min_windows = 3L,
  accept = 0.40, tolerate = 0.22, roll_size = 5L,
  seed_z_min = 1.2, seed_nn_max = 0.90
)

CORE_FAMILIES <- list(S1S, S1M, S1L)

# =============================================================================
# LOAD PRECOMPUTED DATA
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

message("[seeded_regions] Loading precomputed data...")
t_load <- proc.time()

# If single chr requested, only load that one
if (!is.null(chr_name)) {
  target_file <- grep(chr_name, rds_files, value = TRUE)
  if (length(target_file) == 0) stop("No precomp file matching: ", chr_name)
  rds_files <- target_file
} else if (!is.null(chr_idx)) {
  if (chr_idx < 1 || chr_idx > length(rds_files)) stop("chr_idx out of range: ", chr_idx)
  rds_files <- rds_files[chr_idx]
}

obj_list <- lapply(rds_files, readRDS)
precomp_list <- list()
chroms <- character()
for (obj in obj_list) {
  precomp_list[[obj$chrom]] <- obj
  chroms <- c(chroms, obj$chrom)
}
message("[seeded_regions] ", length(chroms), " chromosome(s) in ", round((proc.time() - t_load)[3], 1), "s")

# Sample names
sample_names_snake <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { sample_names_snake <- sub("^PC_1_", "", pc1_cols); break }
}
if (!is.null(sample_names_snake)) message("[seeded_regions] Samples: ", length(sample_names_snake))

# Inv-likeness (for threshold modulation)
inv_like_file <- file.path(dirname(precomp_dir), "window_dt.tsv.gz")
window_dt <- if (file.exists(inv_like_file)) fread(inv_like_file) else data.table()

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

mds_set_sim <- function(mds_mat, idxs_a, idx_b) {
  if (length(idxs_a) == 0 || ncol(mds_mat) == 0) return(0.5)
  center_a <- colMeans(mds_mat[idxs_a, , drop = FALSE], na.rm = TRUE)
  pt_b <- mds_mat[idx_b, ]
  d <- sqrt(sum((center_a - pt_b)^2, na.rm = TRUE))
  spreads <- apply(mds_mat, 2, function(x) {
    iq <- IQR(x, na.rm = TRUE)
    if (is.finite(iq) && iq > 0) iq * 2.5 else diff(range(x, na.rm = TRUE))
  })
  max_sp <- max(spreads, na.rm = TRUE)
  if (!is.finite(max_sp) || max_sp == 0) return(0.5)
  max(0, min(1, 1 - d / max_sp))
}

compute_continuity <- function(j, prev, roll_idxs, sim_mat, mds_mat) {
  s_prev <- if (is.finite(sim_mat[prev, j])) sim_mat[prev, j] else 0
  s_roll <- if (length(roll_idxs) > 0) {
    mean(vapply(roll_idxs, function(r) {
      v <- sim_mat[r, j]; if (is.finite(v)) v else 0
    }, numeric(1)))
  } else s_prev
  s_mds <- mds_set_sim(mds_mat, c(roll_idxs, prev), j)
  W_PREV * s_prev + W_ROLL * s_roll + W_MDS * s_mds
}

region_coherence <- function(idxs, sim_mat) {
  if (length(idxs) < 2) return(1.0)
  vals <- sim_mat[idxs, idxs][upper.tri(sim_mat[idxs, idxs])]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(0)
  mean(vals)
}

decision_log <- list()
log_decision <- function(chr, sid, phase, direction, from_idx, to_idx,
                         from_wid, to_wid, score, accept_thresh, toler_thresh,
                         decision, reason, gap_count, bad_count, roll_size,
                         family_name) {
  decision_log[[length(decision_log) + 1L]] <<- data.table(
    chrom = chr, region_id = sid, phase = phase, direction = direction,
    scale_tier = family_name,
    from_window_idx = from_idx, to_window_idx = to_idx,
    from_window_id = from_wid, to_window_id = to_wid,
    continuity_score = round(score, 4),
    accept_thresh = accept_thresh, toler_thresh = toler_thresh,
    decision = decision, reason = reason,
    gap_count = gap_count, bad_count = bad_count, roll_size = roll_size
  )
}

# =============================================================================
# CORE ENGINE (from STEP_C01b_snake1_run.R, unchanged)
# =============================================================================

run_extension_scale <- function(dt, sim_mat, mds_mat, chr, region_id_start,
                            used_global, params) {
  n <- nrow(dt)
  if (n < params$min_windows) return(list(regions = list(), next_id = region_id_start))

  has_inv <- "inv_likeness" %in% names(dt)
  has_adaptive <- "adaptive_seed" %in% names(dt)  # test_07
  has_spiky <- "spiky_inv_score" %in% names(dt)   # v9 morphology

  # Seed eligibility: z-score OR inv_likeness OR test_07 adaptive_seed
  dt[, fam_seed := (max_abs_z >= params$seed_z_min |
                    (has_inv & is.finite(inv_likeness) & inv_likeness >= 0.90) |
                    (has_adaptive & adaptive_seed == TRUE)) &
       is.finite(seed_nn_dist) &
       seed_nn_dist < params$seed_nn_max]

  # Seed priority: z-score + inv_likeness boost + spiky boost
  dt[, seed_priority := max_abs_z]
  if (has_inv) dt[, seed_priority := pmax(seed_priority, inv_likeness * 5, na.rm = TRUE)]
  if (has_spiky) {
    # Compact candidates get priority boost — small real inversions should
    # be seeded early so they're not swallowed by larger greedy regions.
    dt[, seed_priority := seed_priority +
         fifelse(is.finite(spiky_inv_score) & spiky_inv_score > 0.3,
                 spiky_inv_score * 2, 0)]
  }
  seed_order <- order(-dt$seed_priority)

  n_eligible <- sum(dt$fam_seed, na.rm = TRUE)
  n_z_pass <- sum(dt$max_abs_z >= params$seed_z_min, na.rm = TRUE)
  n_nn_pass <- sum(is.finite(dt$seed_nn_dist) & dt$seed_nn_dist < params$seed_nn_max, na.rm = TRUE)

  # Inv-likeness breakdown at multiple thresholds
  inv_counts <- ""
  if (has_inv) {
    il <- dt$inv_likeness[is.finite(dt$inv_likeness)]
    inv_counts <- paste0(
      " inv>=0.90:", sum(il >= 0.90),
      " >=0.85:", sum(il >= 0.85),
      " >=0.80:", sum(il >= 0.80),
      " >=0.75:", sum(il >= 0.75),
      " >=0.70:", sum(il >= 0.70),
      " >=0.60:", sum(il >= 0.60),
      " <0.50:", sum(il < 0.50))
  }

  # Z-score breakdown
  z_counts <- paste0(
    " z>=5:", sum(dt$max_abs_z >= 5, na.rm = TRUE),
    " >=4:", sum(dt$max_abs_z >= 4, na.rm = TRUE),
    " >=3:", sum(dt$max_abs_z >= 3, na.rm = TRUE),
    " >=2:", sum(dt$max_abs_z >= 2, na.rm = TRUE),
    " >=1.5:", sum(dt$max_abs_z >= 1.5, na.rm = TRUE))

  message("    [seeds] ", params$name, ": ", n_eligible, " eligible / ", n, " windows",
          "  (z>=", params$seed_z_min, ": ", n_z_pass,
          ", nn<", round(params$seed_nn_max, 3), ": ", n_nn_pass, ")")
  message("      z: ", z_counts)
  message("      inv:", inv_counts)

  regions <- list()
  sid <- region_id_start
  claimed <- rep(FALSE, n)

  DMG_RECOVER_ACCEPTED  <- params$dmg_recover  %||% -0.05
  DMG_COST_TOLERATED    <- params$dmg_tolerate %||%  0.04
  DMG_COST_GAP          <- params$dmg_gap      %||%  0.08
  DMG_COST_CONSEC_GAP   <- params$dmg_consec   %||%  0.05
  DMG_MAX               <- params$dmg_max      %||%  0.30

  for (si in seed_order) {
    if (claimed[si] || !dt$fam_seed[si]) next

    sid <- sid + 1L
    region <- si; status <- "seed"; scores <- 1.0
    claimed[si] <- TRUE
    seed_z <- dt$max_abs_z[si]
    seed_pos_mb <- round((dt$start_bp[si] + dt$end_bp[si]) / 2e6, 2)
    seed_inv <- if (has_inv && is.finite(dt$inv_likeness[si])) round(dt$inv_likeness[si], 3) else NA
    right_stop_reason <- "chr_end"; left_stop_reason <- "chr_start"
    right_stop_pos <- NA_real_; left_stop_pos <- NA_real_

    # -- EXTEND RIGHT --
    damage <- 0; consec_gap <- 0L; consec_inv_fail <- 0L; pos <- si
    first_right_score <- NA_real_
    recent_scores <- numeric(0)  # track recent scores for sharp-drop detection
    while (pos < n) {
      pos <- pos + 1L
      if (claimed[pos]) {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        if (damage > DMG_MAX) { right_stop_reason <- "gap_damage"; right_stop_pos <- round((dt$start_bp[pos]+dt$end_bp[pos])/2e6,2); break }
        next
      }

      pos_inv <- if (has_inv && is.finite(dt$inv_likeness[pos])) dt$inv_likeness[pos] else chr_inv_median
      inv_thresh_shift <- 0
      if (has_inv && is.finite(chr_inv_median) && is.finite(chr_inv_q75)) {
        spread <- max(chr_inv_q75 - chr_inv_median, 0.02)
        inv_rel <- (pos_inv - chr_inv_median) / spread
        inv_thresh_shift <- -pmax(-1, pmin(1, inv_rel)) * 0.15
      }

      # ── MORPHOLOGY SOFT DAMAGE MODIFIERS (v9.3.2) ──
      # flat_inv_score: reduce damage when traversing a coherent plateau
      # fragmentation_score: increase damage when entering messy territory
      morph_shift <- 0
      has_flat <- "flat_inv_score" %in% names(dt)
      has_frag <- "fragmentation_score" %in% names(dt)
      if (has_flat && is.finite(dt$flat_inv_score[pos]) && dt$flat_inv_score[pos] > 0.3) {
        morph_shift <- morph_shift - dt$flat_inv_score[pos] * 0.08  # reduce damage (lower bar)
      }
      if (has_frag && is.finite(dt$fragmentation_score[pos]) && dt$fragmentation_score[pos] > 0.3) {
        morph_shift <- morph_shift + dt$fragmentation_score[pos] * 0.10  # increase damage (raise bar)
      }

      if (pos_inv < chr_inv_median) { consec_inv_fail <- consec_inv_fail + 1L } else { consec_inv_fail <- 0L }
      consec_penalty <- if (consec_inv_fail > 1L) (consec_inv_fail - 1L) * 0.03 else 0

      eff_accept   <- params$accept + damage + inv_thresh_shift + morph_shift + consec_penalty
      eff_tolerate <- params$tolerate + damage + inv_thresh_shift + morph_shift + consec_penalty

      roll <- tail(region, params$roll_size)
      prev <- region[length(region)]
      sc <- compute_continuity(pos, prev, roll, sim_mat, mds_mat)
      if (is.na(first_right_score)) first_right_score <- sc

      # SHARP-DROP DETECTION: if score drops significantly below rolling mean,
      # this is a real boundary — instant high damage.
      # A gradual decline accumulates damage slowly. A sharp drop kills the region.
      sharp_penalty <- 0
      if (length(recent_scores) >= 3) {
        roll_mean <- mean(tail(recent_scores, 5))
        drop <- roll_mean - sc
        if (drop > 0.20) {
          sharp_penalty <- drop * 0.5  # proportional to how sharp the drop is
        }
      }
      damage <- damage + sharp_penalty

      if (sc >= eff_accept && sharp_penalty < 0.10) {
        region <- c(region, pos); status <- c(status, "accepted"); scores <- c(scores, sc)
        claimed[pos] <- TRUE; damage <- max(0, damage + DMG_RECOVER_ACCEPTED)
        consec_gap <- 0L; consec_inv_fail <- 0L
        recent_scores <- c(recent_scores, sc)
        log_decision(chr, sid, "extend", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "accepted", paste0("dmg=", round(damage, 3),
                                        " inv_shift=", round(inv_thresh_shift, 3),
                                        " sharp=", round(sharp_penalty, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else if (sc >= eff_tolerate && sharp_penalty < 0.10) {
        region <- c(region, pos); status <- c(status, "tolerated"); scores <- c(scores, sc)
        claimed[pos] <- TRUE; damage <- damage + DMG_COST_TOLERATED
        consec_gap <- 0L; consec_inv_fail <- 0L
        recent_scores <- c(recent_scores, sc)
        log_decision(chr, sid, "extend", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "tolerated", paste0("dmg=", round(damage, 3),
                                         " inv_shift=", round(inv_thresh_shift, 3),
                                         " sharp=", round(sharp_penalty, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        reason_str <- paste0("dmg=", round(damage, 3),
                            " inv_shift=", round(inv_thresh_shift, 3),
                            " consec_inv=", consec_inv_fail,
                            " sharp=", round(sharp_penalty, 3))
        log_decision(chr, sid, "extend", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     if (damage > DMG_MAX) "stopped" else "rejected",
                     reason_str,
                     consec_gap, 0L, length(roll), params$name)
        if (damage > DMG_MAX) {
          right_stop_reason <- paste0("rejected_dmg=",round(damage,3),"_sc=",round(sc,3),
                                      "_inv=",round(pos_inv,3),"_sharp=",round(sharp_penalty,3))
          right_stop_pos <- round((dt$start_bp[pos]+dt$end_bp[pos])/2e6,2)
          break
        }
      }
    }

    # -- EXTEND LEFT --
    damage <- 0; consec_gap <- 0L; consec_inv_fail <- 0L; pos <- si
    recent_scores_l <- numeric(0)
    while (pos > 1L) {
      pos <- pos - 1L
      if (claimed[pos]) {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        if (damage > DMG_MAX) { left_stop_reason <- "gap_damage"; left_stop_pos <- round((dt$start_bp[pos]+dt$end_bp[pos])/2e6,2); break }
        next
      }

      pos_inv_l <- if (has_inv && is.finite(dt$inv_likeness[pos])) dt$inv_likeness[pos] else chr_inv_median
      inv_thresh_shift_l <- 0
      if (has_inv && is.finite(chr_inv_median) && is.finite(chr_inv_q75)) {
        spread <- max(chr_inv_q75 - chr_inv_median, 0.02)
        inv_rel_l <- (pos_inv_l - chr_inv_median) / spread
        inv_thresh_shift_l <- -pmax(-1, pmin(1, inv_rel_l)) * 0.15
      }

      # ── MORPHOLOGY SOFT DAMAGE MODIFIERS (v9.3.2) — same as right ──
      morph_shift_l <- 0
      if (has_flat && is.finite(dt$flat_inv_score[pos]) && dt$flat_inv_score[pos] > 0.3) {
        morph_shift_l <- morph_shift_l - dt$flat_inv_score[pos] * 0.08
      }
      if (has_frag && is.finite(dt$fragmentation_score[pos]) && dt$fragmentation_score[pos] > 0.3) {
        morph_shift_l <- morph_shift_l + dt$fragmentation_score[pos] * 0.10
      }

      if (pos_inv_l < chr_inv_median) { consec_inv_fail <- consec_inv_fail + 1L } else { consec_inv_fail <- 0L }
      consec_penalty_l <- if (consec_inv_fail > 1L) (consec_inv_fail - 1L) * 0.03 else 0

      eff_accept   <- params$accept + damage + inv_thresh_shift_l + morph_shift_l + consec_penalty_l
      eff_tolerate <- params$tolerate + damage + inv_thresh_shift_l + morph_shift_l + consec_penalty_l

      roll <- head(region, params$roll_size)
      prev <- region[1]
      sc <- compute_continuity(pos, prev, roll, sim_mat, mds_mat)

      # SHARP-DROP DETECTION (same as right extension)
      sharp_penalty_l <- 0
      if (length(recent_scores_l) >= 3) {
        roll_mean_l <- mean(tail(recent_scores_l, 5))
        drop_l <- roll_mean_l - sc
        if (drop_l > 0.20) {
          sharp_penalty_l <- drop_l * 0.5
        }
      }
      damage <- damage + sharp_penalty_l

      if (sc >= eff_accept && sharp_penalty_l < 0.10) {
        region <- c(pos, region); status <- c("accepted", status); scores <- c(sc, scores)
        claimed[pos] <- TRUE; damage <- max(0, damage + DMG_RECOVER_ACCEPTED)
        consec_gap <- 0L; consec_inv_fail <- 0L
        recent_scores_l <- c(recent_scores_l, sc)
        log_decision(chr, sid, "extend", "left", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "accepted", paste0("dmg=", round(damage, 3),
                                        " inv_shift=", round(inv_thresh_shift_l, 3),
                                        " sharp=", round(sharp_penalty_l, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else if (sc >= eff_tolerate && sharp_penalty_l < 0.10) {
        region <- c(pos, region); status <- c("tolerated", status); scores <- c(sc, scores)
        claimed[pos] <- TRUE; damage <- damage + DMG_COST_TOLERATED
        consec_gap <- 0L; consec_inv_fail <- 0L
        recent_scores_l <- c(recent_scores_l, sc)
        log_decision(chr, sid, "extend", "left", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "tolerated", paste0("dmg=", round(damage, 3),
                                         " inv_shift=", round(inv_thresh_shift_l, 3),
                                         " sharp=", round(sharp_penalty_l, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        reason_str_l <- paste0("dmg=", round(damage, 3),
                              " inv_shift=", round(inv_thresh_shift_l, 3),
                              " consec_inv=", consec_inv_fail,
                              " sharp=", round(sharp_penalty_l, 3))
        log_decision(chr, sid, "extend", "left", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     if (damage > DMG_MAX) "stopped" else "rejected",
                     reason_str_l,
                     consec_gap, 0L, length(roll), params$name)
        if (damage > DMG_MAX) {
          left_stop_reason <- paste0("rejected_dmg=",round(damage,3),"_sc=",round(sc,3),
                                     "_inv=",round(pos_inv_l,3),"_sharp=",round(sharp_penalty_l,3))
          left_stop_pos <- round((dt$start_bp[pos]+dt$end_bp[pos])/2e6,2)
          break
        }
      }
    }

    # -- Emit or kill --
    region_span_mb <- round((dt$end_bp[max(region)] - dt$start_bp[min(region)]) / 1e6, 2)
    n_acc <- sum(status == "accepted"); n_tol <- sum(status == "tolerated")
    mean_inv <- if (has_inv) round(mean(dt$inv_likeness[region], na.rm = TRUE), 3) else NA
    min_inv <- if (has_inv) round(min(dt$inv_likeness[region], na.rm = TRUE), 3) else NA

    if (length(region) >= params$min_windows) {
      regions[[length(regions) + 1]] <- list(
        idxs = region, statuses = status, scores = scores,
        region_id = sid, scale_tier = params$name
      )
      message("    [ALIVE] #", sid, " seed@", seed_pos_mb, "Mb z=", round(seed_z, 2),
              " inv=", seed_inv, " | ", length(region), " windows (",
              round(dt$start_bp[min(region)]/1e6,2), "-", round(dt$end_bp[max(region)]/1e6,2),
              "Mb, ", region_span_mb, "Mb span) | acc=", n_acc, " tol=", n_tol,
              " | mean_inv=", mean_inv, " min_inv=", min_inv,
              " | R_stop: ", right_stop_reason,
              if (!is.na(right_stop_pos)) paste0(" @", right_stop_pos, "Mb") else "",
              " | L_stop: ", left_stop_reason,
              if (!is.na(left_stop_pos)) paste0(" @", left_stop_pos, "Mb") else "")
    } else {
      n_dead_so_far <- sid - region_id_start - length(regions)
      if (n_dead_so_far <= 30) {
        message("    [DEAD]  #", sid, " seed@", seed_pos_mb, "Mb z=", round(seed_z, 2),
                " inv=", seed_inv, " | grew to ", length(region), " windows (need ",
                params$min_windows, ") | first_R=", round(first_right_score, 3),
                " | R_stop: ", right_stop_reason, " | L_stop: ", left_stop_reason)
      }
      claimed[region] <- FALSE
    }
  }

  # Post-family summary
  n_spawned <- sid - region_id_start
  n_survived <- length(regions)
  total_windows <- sum(vapply(regions, function(r) length(r$idxs), integer(1)))
  if (n_survived > 0) {
    sizes <- vapply(regions, function(r) length(r$idxs), integer(1))
    message("    [summary] ", params$name, ": ", n_spawned, " spawned, ",
            n_survived, " survived (", n_spawned - n_survived, " died), ",
            total_windows, " windows total | sizes: min=", min(sizes),
            " median=", median(sizes), " max=", max(sizes))
  }

  list(regions = regions, next_id = sid)
}

# =============================================================================
# MAIN: PROCESS EACH CHROMOSOME
# =============================================================================

region_id <- 0L

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 3) { message("[SKIP] ", chr); next }

  message("\n[seeded_regions] ======= ", chr, " =======")
  dt <- pc$dt; sim_mat <- pc$sim_mat; mds_mat <- pc$mds_mat
  n <- nrow(dt)

  # Per-chr inv_likeness baseline
  chr_inv_median <- NA_real_; chr_inv_q75 <- NA_real_
  if ("inv_likeness" %in% names(dt)) {
    il_vals <- dt$inv_likeness[is.finite(dt$inv_likeness)]
    if (length(il_vals) > 0) {
      chr_inv_median <- median(il_vals)
      chr_inv_q75 <- quantile(il_vals, 0.75)
    }
  }

  # Adaptive accept thresholds
  bg_q <- pc$bg_continuity_quantiles
  bg_accept_S <- bg_accept_M <- bg_accept_L <- NULL
  if (!is.null(bg_q) && length(bg_q) >= 6) {
    bg_accept_S <- as.numeric(bg_q["95%"])
    bg_accept_M <- as.numeric(bg_q["90%"])
    bg_accept_L <- as.numeric(bg_q["85%"])
    message("    [adaptive] bg: q85=", round(bg_q["85%"], 3),
            " q90=", round(bg_q["90%"], 3), " q95=", round(bg_q["95%"], 3))
  }

  if (!is.na(chr_inv_median)) {
    message("    [inv_like] median=", round(chr_inv_median, 3), " q75=", round(chr_inv_q75, 3))
  }

  # Adaptive NN quantiles
  nn_q <- NULL
  if ("seed_nn_dist" %in% names(dt)) {
    nn_vals <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
    if (length(nn_vals) > 20) nn_q <- quantile(nn_vals, c(0.25, 0.50, 0.75), na.rm = TRUE)
  }

  # Run three scale tiers (1S, 1M, 1L)
  all_regions <- list()
  all_window_rows <- list()
  all_region_rows <- list()

  for (fam in CORE_FAMILIES) {
    fam_adapted <- fam

    # Adaptive accept override
    bg_val <- switch(fam$name, "1S" = bg_accept_S, "1M" = bg_accept_M, "1L" = bg_accept_L, NULL)
    if (!is.null(bg_val) && is.finite(bg_val)) {
      fam_adapted$accept   <- bg_val
      fam_adapted$tolerate <- bg_val - 0.15
      message("    [adaptive] ", fam$name, " accept: ", round(fam$accept, 3),
              " -> ", round(fam_adapted$accept, 3))
    }

    # Adaptive NN override
    if (!is.null(nn_q)) {
      nn_adaptive <- switch(fam$name,
        "1S" = as.numeric(nn_q["25%"]),
        "1M" = as.numeric(nn_q["50%"]),
        "1L" = as.numeric(nn_q["75%"]),
        NULL)
      if (!is.null(nn_adaptive) && is.finite(nn_adaptive)) {
        fam_adapted$seed_nn_max <- nn_adaptive
        message("    [adaptive] ", fam$name, " nn_max: ", round(fam$seed_nn_max, 3),
                " -> ", round(nn_adaptive, 3))
      }
    }

    result <- run_extension_scale(dt, sim_mat, mds_mat, chr, region_id, rep(FALSE, n), fam_adapted)
    region_id <- result$next_id
    message("[seeded_regions] ", fam$name, " (", fam$label, "): ", length(result$regions), " seeded regions")

    for (reg in result$regions) {
      all_regions[[length(all_regions) + 1]] <- reg
      idxs <- reg$idxs
      for (k in seq_along(idxs)) {
        all_window_rows[[length(all_window_rows) + 1]] <- data.table(
          chrom = chr, global_window_id = dt$global_window_id[idxs[k]],
          start_bp = dt$start_bp[idxs[k]], end_bp = dt$end_bp[idxs[k]],
          region_id = reg$region_id, extension_phase = "extend",
          scale_tier = fam$name, merge_family = NA_character_,
          inclusion_status = reg$statuses[k],
          continuity_score = round(reg$scores[k], 4)
        )
      }
      all_region_rows[[length(all_region_rows) + 1]] <- data.table(
        region_id = reg$region_id, extension_phase = "extend",
        scale_tier = fam$name, merge_family = NA_character_, chrom = chr,
        start_bp = min(dt$start_bp[idxs]), end_bp = max(dt$end_bp[idxs]),
        n_windows = length(idxs),
        n_seeds = sum(reg$statuses == "seed"),
        n_tolerated = sum(reg$statuses == "tolerated"),
        mean_score = round(mean(reg$scores), 4),
        min_score = round(min(reg$scores), 4),
        coherence = round(region_coherence(idxs, sim_mat), 4),
        test26_status = reg$test26_status %||% NA_character_,
        test26_retention = reg$test26_retention %||% NA_real_
      )
    }
  }

  message("[seeded_regions] Total: ", length(all_regions), " seeded regions")

  # ═══════════════════════════════════════════════════════════════════
  # TEST 26: KIN-PRUNED RETENTION TEST (v9.3.2)
  # ═══════════════════════════════════════════════════════════════════
  # For each seeded region, test if the sim_mat block contrast survives after
  # removing kin-pruned samples. Cores where signal collapses get
  # FLAGGED (not killed) — in a hatchery, collapse can mean either
  # "family artifact" OR "real founder-linked inversion."
  #
  # The flag is stored in the region output and window states for
  # downstream scoring to use as one piece of evidence.
  # ═══════════════════════════════════════════════════════════════════

  pruned_ind_names <- NULL
  if (!is.null(pruned_file) && file.exists(pruned_file)) {
    pruned_raw <- as.character(fread(pruned_file, header = FALSE)[[1]])
    # Map real names to Ind names if needed
    if (length(pruned_raw) > 0 && !grepl("^Ind", pruned_raw[1])) {
      # Try to build mapping from sample_names
      pc1_c <- grep("^PC_1_", names(dt), value = TRUE)
      ind_names <- sub("^PC_1_", "", pc1_c)
      # Assume BAM list order: real_names[1] = Ind0, etc.
      # Load sample list if available
      samples_ind_file <- Sys.getenv("SAMPLES_IND", "")
      if (nzchar(samples_ind_file) && file.exists(samples_ind_file)) {
        real_names <- trimws(readLines(samples_ind_file))
        real_names <- real_names[nzchar(real_names)]
        if (length(real_names) == length(ind_names)) {
          idx <- match(pruned_raw, real_names)
          pruned_ind_names <- paste0("Ind", idx[!is.na(idx)] - 1L)
        }
      }
    } else {
      pruned_ind_names <- pruned_raw
    }
    if (!is.null(pruned_ind_names)) {
      message("[test_26] Pruned samples loaded: ", length(pruned_ind_names))
    }
  }

  if (!is.null(pruned_ind_names) && length(pruned_ind_names) >= 15 && length(all_regions) > 0) {
    message("[test_26] Testing kin-pruned retention for ", length(all_regions), " seeded regions...")
    pruned_pc1_cols <- paste0("PC_1_", pruned_ind_names)
    avail_pruned <- intersect(pruned_pc1_cols, names(dt))

    if (length(avail_pruned) >= 15) {
      for (ci in seq_along(all_regions)) {
        reg <- all_regions[[ci]]
        idxs <- reg$idxs
        if (length(idxs) < 5) { reg$test26_status <- "too_small"; next }

        # Full block contrast (all 226 samples)
        out_idxs <- setdiff(seq_len(n), idxs)
        if (length(out_idxs) < 10) out_idxs <- seq_len(min(50, n))
        full_inside <- mean(sim_mat[idxs, idxs], na.rm = TRUE)
        full_outside <- mean(sim_mat[out_idxs, out_idxs], na.rm = TRUE)
        full_contrast <- full_inside - full_outside

        # Pruned block contrast (recompute sim from pruned PC loadings)
        sub_idxs <- sort(c(sample(out_idxs, min(100, length(out_idxs))),
                            idxs))
        p_mat <- as.matrix(dt[sub_idxs, ..avail_pruned])
        ns <- nrow(p_mat)
        p_sim <- matrix(0, ns, ns)
        for (ii in seq_len(ns)) for (jj in ii:ns) {
          v <- cor(p_mat[ii, ], p_mat[jj, ], use = "pairwise.complete.obs")
          if (is.finite(v)) { p_sim[ii, jj] <- v; p_sim[jj, ii] <- v }
        }
        win_in <- which(sub_idxs %in% idxs)
        out_in <- which(!sub_idxs %in% idxs)

        if (length(win_in) >= 3 && length(out_in) >= 3) {
          p_inside <- mean(p_sim[win_in, win_in], na.rm = TRUE)
          p_outside <- mean(p_sim[out_in, out_in], na.rm = TRUE)
          p_contrast <- p_inside - p_outside
          retention <- if (full_contrast > 0.005) p_contrast / full_contrast else NA

          reg$test26_retention <- round(retention, 3)
          reg$test26_status <- if (!is.finite(retention)) "cannot_compute"
            else if (retention > 0.7) "persists"
            else if (retention > 0.4) "weakened"
            else if (retention > 0.1) "mostly_lost"
            else "collapsed"
        } else {
          reg$test26_status <- "insufficient_pruned"
          reg$test26_retention <- NA_real_
        }

        all_regions[[ci]] <- reg
        if (ci <= 5 || ci %% 20 == 0) {
          message("  region ", ci, " (", reg$scale_tier, " region_id=", reg$region_id,
                  "): test26=", reg$test26_status,
                  " retention=", reg$test26_retention)
        }
      }

      # Summary
      c26_tab <- table(vapply(all_regions, function(r) r$test26_status %||% "untested", character(1)))
      message("[test_26] Results: ", paste(names(c26_tab), c26_tab, sep = "=", collapse = " "))
    } else {
      message("[test_26] Too few pruned PC columns (", length(avail_pruned), ") — skipping")
    }
  } else if (is.null(pruned_file)) {
    message("[test_26] No --pruned_samples provided — skipping retention test")
  }

  # =================================================================
  # ROARY-STYLE WINDOW × CORE PRESENCE/ABSENCE MATRIX
  # =================================================================
  # Rows = ALL windows on this chromosome
  # Fixed columns: chrom, global_window_id, start_bp, end_bp, pos_mb,
  #                max_abs_z, inv_likeness, seed_nn_dist, is_seed_like
  # Per-family columns: in_1S, in_1M, in_1L (TRUE/FALSE)
  # Per-scale region ID: region_1S, region_1M, region_1L (region_id or NA)
  # Per-family status: status_1S, status_1M, status_1L (seed/accepted/tolerated/NA)
  # Per-family score: score_1S, score_1M, score_1L (continuity score or NA)
  # Derived columns:
  #   n_families_claimed: how many families claimed this window (0-3)
  #   all_families_agree: TRUE if all claiming families agree (or only 1 claims)
  #   pa_pattern: e.g. "SML" if all three, "SM" if strict+moderate, "L" if broad only
  #   collector_state: collected / seed_uncollected / unused (backward compatible)

  # Build per-scale lookup: for each scale tier, which windows → which region → what status
  fam_lookup <- list()
  for (fam_name in c("1S", "1M", "1L")) {
    fam_regions <- all_regions[vapply(all_regions, function(r) r$scale_tier == fam_name, logical(1))]
    wid_to_core <- list()
    for (reg in fam_regions) {
      for (k in seq_along(reg$idxs)) {
        wi <- reg$idxs[k]
        wid <- dt$global_window_id[wi]
        wid_to_core[[as.character(wid)]] <- list(
          region_id = reg$region_id,
          status = reg$statuses[k],
          score = reg$scores[k]
        )
      }
    }
    fam_lookup[[fam_name]] <- wid_to_core
  }

  # Build the full matrix
  all_state_rows <- list()
  for (wi in seq_len(n)) {
    wid <- dt$global_window_id[wi]
    wid_str <- as.character(wid)
    has_inv_col <- "inv_likeness" %in% names(dt)
    inv_l <- if (has_inv_col && is.finite(dt$inv_likeness[wi])) dt$inv_likeness[wi] else 0
    nn_d <- if (is.finite(dt$seed_nn_dist[wi])) dt$seed_nn_dist[wi] else NA_real_
    is_seed_like <- (dt$max_abs_z[wi] >= 1.2 | inv_l >= 0.90) &
                    is.finite(nn_d) & nn_d < 0.90

    # Per-family info
    s_info <- fam_lookup[["1S"]][[wid_str]]
    m_info <- fam_lookup[["1M"]][[wid_str]]
    l_info <- fam_lookup[["1L"]][[wid_str]]

    in_1S <- !is.null(s_info); in_1M <- !is.null(m_info); in_1L <- !is.null(l_info)
    n_fam <- sum(c(in_1S, in_1M, in_1L))

    # PA pattern string
    pa <- paste0(if (in_1S) "S" else "", if (in_1M) "M" else "", if (in_1L) "L" else "")
    if (nchar(pa) == 0) pa <- "none"

    # Backward-compatible collector_state
    cstate <- if (n_fam > 0) "collected"
              else if (is_seed_like) "seed_uncollected"
              else "unused"

    all_state_rows[[wi]] <- data.table(
      chrom = chr,
      global_window_id = wid,
      start_bp = dt$start_bp[wi],
      end_bp = dt$end_bp[wi],
      pos_mb = round((dt$start_bp[wi] + dt$end_bp[wi]) / 2e6, 4),
      max_abs_z = round(dt$max_abs_z[wi], 4),
      inv_likeness = round(inv_l, 4),
      seed_nn_dist = round(nn_d, 4),
      is_seed_like = is_seed_like,
      # Per-family PA
      in_1S = in_1S, in_1M = in_1M, in_1L = in_1L,
      region_1S = if (in_1S) s_info$region_id else NA_integer_,
      region_1M = if (in_1M) m_info$region_id else NA_integer_,
      region_1L = if (in_1L) l_info$region_id else NA_integer_,
      status_1S = if (in_1S) s_info$status else NA_character_,
      status_1M = if (in_1M) m_info$status else NA_character_,
      status_1L = if (in_1L) l_info$status else NA_character_,
      score_1S = if (in_1S) round(s_info$score, 4) else NA_real_,
      score_1M = if (in_1M) round(m_info$score, 4) else NA_real_,
      score_1L = if (in_1L) round(l_info$score, 4) else NA_real_,
      # Derived
      n_families_claimed = n_fam,
      pa_pattern = pa,
      collector_state = cstate
    )
  }

  # -- SUMMARY ROW --
  n_s <- sum(vapply(all_regions, function(r) r$scale_tier == "1S", logical(1)))
  n_m <- sum(vapply(all_regions, function(r) r$scale_tier == "1M", logical(1)))
  n_l <- sum(vapply(all_regions, function(r) r$scale_tier == "1L", logical(1)))

  state_dt <- rbindlist(all_state_rows)
  n_collected <- sum(state_dt$collector_state == "collected")
  n_seed_unc  <- sum(state_dt$collector_state == "seed_uncollected")
  n_all3      <- sum(state_dt$n_families_claimed == 3)
  n_2fam      <- sum(state_dt$n_families_claimed == 2)
  n_1fam      <- sum(state_dt$n_families_claimed == 1)

  summ_row <- data.table(
    chrom = chr, n_windows = n,
    n_regions_1S = n_s, n_regions_1M = n_m, n_regions_1L = n_l,
    n_regions_total = length(all_regions),
    windows_collected = n_collected,
    windows_seed_uncollected = n_seed_unc,
    windows_all3 = n_all3, windows_2fam = n_2fam, windows_1fam = n_1fam,
    pct_collected = round(100 * n_collected / n, 1),
    pct_seed_uncollected = round(100 * n_seed_unc / n, 1)
  )

  message("[seeded_regions] PA summary: collected=", n_collected,
          " (all3=", n_all3, " 2fam=", n_2fam, " 1fam=", n_1fam,
          ") seed_uncollected=", n_seed_unc, " unused=", n - n_collected - n_seed_unc)

  # -- SAVE --
  saveRDS(list(chrom = chr, regions = all_regions), file.path(outdir, paste0("seeded_regions_", chr, ".rds")))

  win_dt   <- if (length(all_window_rows) > 0) rbindlist(all_window_rows, fill = TRUE) else data.table()
  reg_dt   <- if (length(all_region_rows) > 0) rbindlist(all_region_rows, fill = TRUE) else data.table()
  log_dt   <- if (length(decision_log) > 0) rbindlist(decision_log, fill = TRUE) else data.table()

  fwrite(win_dt,    file.path(outdir, paste0("seeded_regions_windows_",  chr, ".tsv.gz")), sep = "\t")
  fwrite(reg_dt,    file.path(outdir, paste0("seeded_regions_summary_",  chr, ".tsv.gz")), sep = "\t")
  fwrite(log_dt,    file.path(outdir, paste0("seeded_regions_decision_log_",  chr, ".tsv.gz")), sep = "\t")
  fwrite(state_dt,  file.path(outdir, paste0("seeded_regions_window_states_", chr, ".tsv.gz")), sep = "\t")
  fwrite(summ_row,  file.path(outdir, paste0("seeded_regions_summary_per_chr_",       chr, ".tsv")),    sep = "\t")

  # Reset decision log for next chr
  decision_log <<- list()

  message("[seeded_regions] Saved: ", chr, " (", length(all_regions), " seeded regions, ",
          nrow(win_dt), " windows, ", nrow(state_dt), " states, ",
          nrow(log_dt), " decisions)")
}

message("\n[DONE] STEP_C01b_1_cores complete")
