#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01b_2_region_merge.R
#
# Phase 2 / 2d — seeded-region merge engine.
# Reads per-chromosome seeded regions (from C01b_1) and fuses adjacent /
# overlapping regions into merged candidates using fuzzy max-min composition
# of two parallel relations, adjusted by the landscape scores from
# PHASE_01C (section "landscape integration" below).
#
# Codebase:    inversion_modules v8.5 / script v9.3.x
# Upstream:    phase_2/2c STEP_C01b_1 output — <outdir>/seeded_regions_<chr>.rds
#              phase_2/2c STEP_C01a output — <precomp_dir>/precomp/<chr>.precomp.rds
#              phase_2/2c PHASE_01C output — <landscape_dir>/
#                block_concordance_<chr>.tsv.gz
#                blue_cross_verdicts_<chr>.tsv.gz
#                boundary_catalog_<chr>.tsv.gz
# Downstream:  phase_4/4a STEP_C01d candidate scoring
#              (consumes merged regions via the bridge interface — see
#              phase_2/2e STEP_D12_bridge_to_C01d.R for the matrix track,
#              and legacy triangle_intervals.tsv.gz format)
# Formerly:    STEP_C01b_2_merge.R
#              The pipeline's "snake" / "core" terminology has been retired;
#              see 2c_precomp/RENAMING.md for the full map.
#
# Algorithm — fuzzy max-min composition
# -------------------------------------
# For every adjacent pair of seeded regions (ordered by chromosomal start),
# compute two compatibility scores:
#
#   Relation M (membership compatibility). Soft cluster-membership for
#     samples in each region is computed from PC1 scores (k=3 cluster).
#     The pairwise membership compatibility is max-min composition of
#     the two soft-assignment matrices: T = R ∘ S.
#
#   Relation G (geometric continuity). Mean sim_mat similarity of the
#     tail windows of region A to the head windows of region B (plus
#     the bridge range between them if within max_gap_windows).
#
# The two relations are combined under one of three operators per tier:
#   accept_threshold controls how high the combined score must be;
#   max_gap_windows controls proximity;
#   operator ∈ {MAX, weighted_mean, MIN} controls strictness.
#
# Landscape integration (v8.5)
# ----------------------------
# PHASE_01C produces three landscape tracks per chromosome. These are
# consulted to adjust the fuzzy merge score:
#
#   block_concordance — do the two regions share carriers under PHASE_01C's
#                        block-level sample partitioning? (+ boost if yes)
#   blue_cross_verdicts — is the boundary between them tagged as an
#                        assembly error (bridgeable) or a real transition
#                        (do not bridge)?
#   boundary_catalog — boundary type between the two regions. A
#                        clear_hard boundary penalises the merge score;
#                        soft / inner / diffuse boundaries are neutral or
#                        mildly positive.
#
# The adjustment is applied in landscape_adjust(), not as a hard gate.
#
# Three scale tiers
# -----------------
# As with C01b_1, three parameter tiers run independently:
#   1S (small)  — tight accept_threshold, small max_gap, MAX operator
#   1M (medium) — default
#   1L (large)  — relaxed accept_threshold, larger max_gap, MIN operator
#
# Usage
# -----
#   Rscript STEP_C01b_2_region_merge.R <precomp_dir> <outdir>
#     [--tag <label>] [--landscape_dir <dir>]
#
# Inputs (found under <outdir>/):
#   seeded_regions_<chr>.rds         — per-chr regions from C01b_1
# Inputs (optional, found under <landscape_dir>/):
#   block_concordance_<chr>.tsv.gz   — from PHASE_01C
#   blue_cross_verdicts_<chr>.tsv.gz — from PHASE_01C
#   boundary_catalog_<chr>.tsv.gz    — from PHASE_01C
#
# Outputs (in <outdir>/)
# ----------------------
#   region_merge_<tier>_<chr>.tsv.gz       per-tier merged regions
#   region_candidate_regions.tsv.gz        union across tiers
#   region_merge_scores.tsv.gz             per-pair fuzzy score log
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01b_2_merge.R <precomp_dir> <outdir> [--tag <label>]")
}

precomp_dir <- args[1]
outdir      <- args[2]

run_tag <- NULL
landscape_dir <- NULL
i <- 3L
while (i <= length(args)) {
  if (args[i] == "--tag" && i < length(args)) {
    run_tag <- args[i + 1]; i <- i + 2L
  } else if (args[i] == "--landscape_dir" && i < length(args)) {
    landscape_dir <- args[i + 1]; i <- i + 2L
  } else { i <- i + 1L }
}

if (!is.null(run_tag)) outdir <- file.path(outdir, run_tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# MERGE PARAMETERS
# =============================================================================
#
# Three tiers controlled by accept_threshold + combination operator:
#   merge_A: low threshold + MAX operator (either signal is enough)
#   merge_B: medium threshold + weighted mean (both contribute)
#   merge_C: high threshold + MIN operator (both must be good)
#
# w_membership / w_geometric control the weighted mean blend for merge_B.

MA <- list(
  name = "merge_A", label = "generous",
  max_gap_windows = 15L,
  accept_threshold = 0.30,
  min_windows = 5L,
  w_membership = 0.6, w_geometric = 0.4
)

MB <- list(
  name = "merge_B", label = "moderate",
  max_gap_windows = 8L,
  accept_threshold = 0.50,
  min_windows = 4L,
  w_membership = 0.6, w_geometric = 0.4
)

MC <- list(
  name = "merge_C", label = "strict",
  max_gap_windows = 3L,
  accept_threshold = 0.70,
  min_windows = 6L,
  w_membership = 0.6, w_geometric = 0.4
)

MERGE_FAMILIES <- list(MA, MB, MC)

# QC parameters
QC_INTERNAL_DROP  <- 0.25
QC_DROP_FRAC      <- 0.15

# Landscape integration parameters (v8.5)
# These adjust the fuzzy merge score based on PHASE_01C outputs.
CONC_ALLOW_THRESH  <- 0.70   # block concordance > this → allow merge
CONC_BLOCK_THRESH  <- 0.40   # block concordance < this → block merge
BLUE_CROSS_BOOST   <- 0.15   # assembly_error / same_system → boost score
BLUE_CROSS_PENALTY <- -0.20  # different_systems → penalize score
BOUNDARY_HARD_PEN  <- -0.15  # clear_hard boundary between regions → penalize
BOUNDARY_INNER_BOOST <- 0.10 # inner_hard + same_system → boost

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

region_coherence <- function(idxs, sim_mat) {
  if (length(idxs) < 2) return(1.0)
  vals <- numeric(length(idxs) - 1L)
  for (k in seq_along(vals)) {
    v <- sim_mat[idxs[k], idxs[k + 1L]]
    vals[k] <- if (is.finite(v)) v else 0
  }
  mean(vals, na.rm = TRUE)
}

# ── SOFT BAND MEMBERSHIP ─────────────────────────────────────────────
# k=3 cluster on average PC1 scores, then compute fuzzy membership
# using inverse-distance to each center. Returns N_samples x 3 matrix
# with rows summing to 1. Epsilon is adaptive to center spread.

soft_band_membership <- function(dt, idxs, sample_names) {
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  n_samples <- length(available)
  if (n_samples < 10) return(NULL)

  mat <- as.matrix(dt[idxs, ..available])
  avg_loadings <- if (nrow(mat) > 1) colMeans(mat, na.rm = TRUE) else mat[1, ]

  valid <- is.finite(avg_loadings)
  if (sum(valid) < 10) return(NULL)

  vals <- avg_loadings[valid]
  snames <- sub("^PC_1_", "", names(vals))

  km <- tryCatch(kmeans(vals, centers = 3, nstart = 5), error = function(e) NULL)
  if (is.null(km)) return(NULL)

  # Order centers low -> mid -> high
  center_order <- order(km$centers[, 1])
  centers <- km$centers[center_order, 1]

  # Adaptive epsilon: proportional to spread between centers
  center_range <- diff(range(centers))
  eps <- max(0.01, center_range * 0.02)

  n_valid <- length(vals)
  membership <- matrix(0, nrow = n_valid, ncol = 3)
  for (bi in 1:3) {
    membership[, bi] <- 1 / (abs(vals - centers[bi]) + eps)
  }

  row_sums <- rowSums(membership)
  membership <- membership / row_sums

  colnames(membership) <- paste0("band", 1:3)
  rownames(membership) <- snames
  membership
}

# ── GEOMETRIC CONTINUITY ─────────────────────────────────────────────

geometric_continuity <- function(region_a_idxs, region_b_idxs, sim_mat,
                                  gap_idxs = integer(0)) {
  tail_a <- tail(region_a_idxs, 3L)
  head_b <- head(region_b_idxs, 3L)
  scores <- numeric(0)

  if (length(gap_idxs) > 0) {
    for (gi in head(gap_idxs, 3L)) {
      for (ai in tail_a) {
        s <- sim_mat[ai, gi]
        if (is.finite(s)) scores <- c(scores, s)
      }
    }
    for (gi in tail(gap_idxs, 3L)) {
      for (bi in head_b) {
        s <- sim_mat[gi, bi]
        if (is.finite(s)) scores <- c(scores, s)
      }
    }
  }

  for (ai in tail_a) {
    for (bi in head_b) {
      s <- sim_mat[ai, bi]
      if (is.finite(s)) scores <- c(scores, s)
    }
  }

  if (length(scores) == 0) return(0)
  mean(scores, na.rm = TRUE)
}

# ── FUZZY MAX-MIN COMPOSITION ────────────────────────────────────────
# T[band_a, band_b] = max_sample { min(R[sample,band_a], S[sample,band_b]) }

fuzzy_compose <- function(membership_a, membership_b) {
  shared <- intersect(rownames(membership_a), rownames(membership_b))
  if (length(shared) < 10) return(NULL)

  R <- membership_a[shared, , drop = FALSE]
  S <- membership_b[shared, , drop = FALSE]

  T_mat <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      T_mat[i, j] <- max(pmin(R[, i], S[, j]), na.rm = TRUE)
    }
  }
  colnames(T_mat) <- rownames(T_mat) <- paste0("band", 1:3)
  T_mat
}

# ── COMPOSITION SCORE ────────────────────────────────────────────────
# Extract [0,1] from the 3x3 T matrix.
# High = bands map cleanly (same inversion). Low = diffuse (different).

composition_score <- function(T_mat) {
  if (is.null(T_mat)) return(NA_real_)

  row_max <- apply(T_mat, 1, max)

  best_cols <- apply(T_mat, 1, which.max)
  n_unique <- length(unique(best_cols))
  perm_score <- if (n_unique == 3) 1.0 else if (n_unique == 2) 0.5 else 0.0

  contrast <- mean(apply(T_mat, 1, function(r) {
    sorted <- sort(r, decreasing = TRUE)
    if (length(sorted) >= 2 && sorted[1] > 0) {
      (sorted[1] - sorted[2]) / sorted[1]
    } else 0
  }))

  score <- 0.5 * mean(row_max) + 0.3 * contrast + 0.2 * perm_score
  max(0, min(1, score))
}

# ── COMBINED FUZZY MERGE SCORE ───────────────────────────────────────

compute_merge_score <- function(dt, region_a_idxs, region_b_idxs, sample_names,
                                 sim_mat, combination, w_m = 0.6, w_g = 0.4) {
  gap_start <- max(region_a_idxs) + 1L
  gap_end   <- min(region_b_idxs) - 1L
  gap_idxs <- if (gap_end >= gap_start) seq(gap_start, gap_end) else integer(0)
  gap_idxs <- gap_idxs[gap_idxs >= 1 & gap_idxs <= nrow(dt)]

  # Relation M
  memb_a <- soft_band_membership(dt, region_a_idxs, sample_names)
  memb_b <- soft_band_membership(dt, region_b_idxs, sample_names)

  if (!is.null(memb_a) && !is.null(memb_b)) {
    T_mat <- fuzzy_compose(memb_a, memb_b)
    m_score <- composition_score(T_mat)
  } else {
    m_score <- NA_real_
  }

  # Relation G
  g_score <- geometric_continuity(region_a_idxs, region_b_idxs, sim_mat, gap_idxs)

  # Combine
  if (is.na(m_score) && g_score == 0) {
    final <- 0
  } else if (is.na(m_score)) {
    final <- g_score
  } else if (g_score == 0) {
    final <- m_score
  } else {
    final <- switch(combination,
      "max"  = max(m_score, g_score),
      "min"  = min(m_score, g_score),
      "mean" = w_m * m_score + w_g * g_score
    )
  }

  list(
    score = max(0, min(1, final)),
    membership_score = m_score %||% NA_real_,
    geometric_score = g_score,
    gap_size = length(gap_idxs)
  )
}

# ── LANDSCAPE ADJUSTMENT (v8.5) ──────────────────────────────────────
# Adjusts the raw fuzzy merge score using PHASE_01C landscape data:
#   1. Block concordance: do the two regions share carriers?
#   2. Blue-cross verdicts: is the gap an assembly error or real boundary?
#   3. Boundary catalog: what type of boundary sits between the regions?
# Returns: list(adjusted_score, adjustment, reason)

landscape_adjust <- function(raw_score, dt, region_a_idxs, region_b_idxs,
                              conc_dt_chr, cross_dt_chr, bound_dt_chr,
                              block_pa_chr) {
  adj <- 0
  reasons <- character(0)


  # Skip if no landscape data
  if (is.null(conc_dt_chr) && is.null(cross_dt_chr) && is.null(bound_dt_chr)) {
    return(list(adjusted_score = raw_score, adjustment = 0, reason = "no_landscape"))
  }

  a_end_bp   <- max(dt$end_bp[region_a_idxs])
  b_start_bp <- min(dt$start_bp[region_b_idxs])
  a_mid_bp   <- (min(dt$start_bp[region_a_idxs]) + a_end_bp) / 2
  b_mid_bp   <- (b_start_bp + max(dt$end_bp[region_b_idxs])) / 2

  # --- 1. Block concordance ---
  # Find which blocks the two regions fall in
  if (!is.null(block_pa_chr) && nrow(block_pa_chr) > 0) {
    a_blocks <- unique(block_pa_chr[start_bp >= min(dt$start_bp[region_a_idxs]) &
                                     end_bp <= a_end_bp &
                                     block != "bg", block])
    b_blocks <- unique(block_pa_chr[start_bp >= b_start_bp &
                                     end_bp <= max(dt$end_bp[region_b_idxs]) &
                                     block != "bg", block])

    if (length(a_blocks) > 0 && length(b_blocks) > 0 && !is.null(conc_dt_chr) && nrow(conc_dt_chr) > 0) {
      # Find the best concordance between any block pair
      best_conc <- NA_real_
      for (ba in a_blocks) {
        for (bb in b_blocks) {
          cc_row <- conc_dt_chr[block_a == ba & block_b == bb]
          if (nrow(cc_row) > 0 && is.finite(cc_row$concordance[1])) {
            if (is.na(best_conc) || cc_row$concordance[1] > best_conc)
              best_conc <- cc_row$concordance[1]
          }
        }
      }
      if (!is.na(best_conc)) {
        if (best_conc >= CONC_ALLOW_THRESH) {
          adj <- adj + 0.05
          reasons <- c(reasons, paste0("conc_high=", round(best_conc, 2)))
        } else if (best_conc < CONC_BLOCK_THRESH) {
          adj <- adj - 0.25
          reasons <- c(reasons, paste0("conc_block=", round(best_conc, 2)))
        }
      }
    }
  }

  # --- 2. Blue-cross verdicts in the gap ---
  if (!is.null(cross_dt_chr) && nrow(cross_dt_chr) > 0) {
    gap_crosses <- cross_dt_chr[start_bp >= a_end_bp & end_bp <= b_start_bp]
    if (nrow(gap_crosses) > 0 && "inner_type" %in% names(gap_crosses)) {
      for (ci in seq_len(nrow(gap_crosses))) {
        it <- gap_crosses$inner_type[ci]
        if (is.na(it)) next
        if (grepl("assembly|same_system", it)) {
          adj <- adj + BLUE_CROSS_BOOST
          reasons <- c(reasons, paste0("blue_cross_boost:", it))
        } else if (grepl("different_systems", it)) {
          adj <- adj + BLUE_CROSS_PENALTY
          reasons <- c(reasons, paste0("blue_cross_penalty:", it))
        }
      }
    }
  }

  # --- 3. Boundary catalog between regions ---
  if (!is.null(bound_dt_chr) && nrow(bound_dt_chr) > 0) {
    # Find boundaries whose position falls in the gap
    gap_bounds <- bound_dt_chr[pos_bp >= a_end_bp & pos_bp <= b_start_bp]
    if (nrow(gap_bounds) > 0) {
      for (bi in seq_len(nrow(gap_bounds))) {
        bt <- gap_bounds$boundary_type[bi]
        if (is.na(bt)) next
        if (bt == "clear_hard") {
          adj <- adj + BOUNDARY_HARD_PEN
          reasons <- c(reasons, "boundary_clear_hard")
        } else if (grepl("^inner_hard", bt) && grepl("same_system", bt)) {
          adj <- adj + BOUNDARY_INNER_BOOST
          reasons <- c(reasons, "boundary_inner_same_sys")
        }
      }
    }
  }

  final <- max(0, min(1, raw_score + adj))
  list(
    adjusted_score = final,
    adjustment = round(adj, 4),
    reason = if (length(reasons) > 0) paste(reasons, collapse = ";") else "none"
  )
}

# ── FUZZY MERGE ENGINE ───────────────────────────────────────────────

run_merge_fuzzy <- function(all_regions, dt, sample_names, sim_mat, chr, params,
                            conc_dt_chr = NULL, cross_dt_chr = NULL,
                            bound_dt_chr = NULL, block_pa_chr = NULL) {
  if (length(all_regions) == 0) return(list(regions = list(), score_log = list()))

  region_starts <- vapply(all_regions, function(c) min(c$idxs), integer(1))
  all_regions <- all_regions[order(region_starts)]

  combination <- if (params$accept_threshold < 0.40) "max"
                 else if (params$accept_threshold < 0.65) "mean"
                 else "min"

  w_m <- params$w_membership %||% 0.6
  w_g <- params$w_geometric  %||% 0.4

  merged <- list()
  score_log <- list()
  i <- 1L

  while (i <= length(all_regions)) {
    current <- all_regions[[i]]
    cur_idxs <- current$idxs
    cur_stat <- current$statuses
    cur_scor <- current$scores

    while (i < length(all_regions)) {
      nxt <- all_regions[[i + 1]]
      nxt_start <- min(nxt$idxs)
      cur_end   <- max(cur_idxs)
      gap_size  <- nxt_start - cur_end - 1L

      # Gate 1: proximity
      if (gap_size > params$max_gap_windows || gap_size < 0) break

      # Gate 2: fuzzy score (rolling reference from tail of growing region)
      ref_idxs <- tail(cur_idxs, 30L)
      ms <- compute_merge_score(dt, ref_idxs, nxt$idxs, sample_names,
                                 sim_mat, combination, w_m, w_g)

      # Gate 3 (v8.5): landscape adjustment from PHASE_01C
      la <- landscape_adjust(ms$score, dt, ref_idxs, nxt$idxs,
                              conc_dt_chr, cross_dt_chr, bound_dt_chr,
                              block_pa_chr)
      final_score <- la$adjusted_score

      score_log[[length(score_log) + 1]] <- data.table(
        chrom = chr, merge_family = params$name,
        cur_end_bp = dt$end_bp[cur_end],
        nxt_start_bp = dt$start_bp[nxt_start],
        gap_windows = gap_size,
        fuzzy_score = round(ms$score, 4),
        membership_score = round(ms$membership_score, 4),
        geometric_score = round(ms$geometric_score, 4),
        landscape_adj = la$adjustment,
        landscape_reason = la$reason,
        adjusted_score = round(final_score, 4),
        threshold = params$accept_threshold,
        combination = combination,
        decision = if (final_score >= params$accept_threshold) "merge" else "reject"
      )

      if (final_score < params$accept_threshold) break

      # Merge
      bridge_range <- if (gap_size > 0) seq(cur_end + 1L, nxt_start - 1L) else integer(0)
      bridge_range <- bridge_range[bridge_range >= 1 & bridge_range <= nrow(dt)]

      cur_idxs <- c(cur_idxs, bridge_range, nxt$idxs)
      cur_stat <- c(cur_stat, rep("bridge_fuzzy", length(bridge_range)), nxt$statuses)
      cur_scor <- c(cur_scor, rep(NA_real_, length(bridge_range)), nxt$scores)
      i <- i + 1L
    }

    if (length(cur_idxs) >= params$min_windows) {
      merged[[length(merged) + 1]] <- list(
        idxs = cur_idxs, statuses = cur_stat, scores = cur_scor,
        merge_family = params$name, coherence = NULL
      )
    }

    i <- i + 1L
  }

  list(regions = merged, score_log = score_log)
}

# ── QC ENGINE ────────────────────────────────────────────────────────

run_qc <- function(merged, sim_mat) {
  if (length(merged) == 0) return(data.table())
  diag_rows <- list()
  for (mi in seq_along(merged)) {
    reg <- merged[[mi]]; idxs <- reg$idxs; n <- length(idxs)
    flags <- character(0)

    isims <- numeric(n); isims[1] <- 1.0
    for (k in 2:n) {
      s <- sim_mat[idxs[k-1], idxs[k]]
      isims[k] <- if (is.finite(s)) s else 0
    }
    if (sum(isims < QC_INTERNAL_DROP) / n > QC_DROP_FRAC) flags <- c(flags, "internal_disc")

    n_bridge <- sum(grepl("bridge", reg$statuses), na.rm = TRUE)
    n_original <- sum(reg$statuses %in% c("seed", "accepted"), na.rm = TRUE)
    if (n_original > 0 && n_bridge / n_original > 0.5) flags <- c(flags, "overmerge")

    coh <- reg$coherence %||% region_coherence(idxs, sim_mat)
    if (coh < 0.25) flags <- c(flags, "low_coherence")
    if (length(flags) == 0) flags <- "clean"

    diag_rows[[mi]] <- data.table(
      merged_idx = mi, n_windows = n, coherence = round(coh, 4),
      mean_isim = round(mean(isims), 4), min_isim = round(min(isims), 4),
      n_bridge = n_bridge, bridge_frac = round(n_bridge / n, 3),
      qc_flags = paste(flags, collapse = ";"),
      merge_family = reg$merge_family %||% "unknown"
    )
  }
  rbindlist(diag_rows)
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("[merge] Loading precomputed data...")
t_load <- proc.time()

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

N_LOAD_CORES <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), length(rds_files))
if (N_LOAD_CORES > 1) {
  obj_list <- parallel::mclapply(rds_files, readRDS, mc.cores = N_LOAD_CORES)
} else {
  obj_list <- lapply(rds_files, readRDS)
}

precomp_list <- list()
chroms <- character()
for (obj in obj_list) {
  precomp_list[[obj$chrom]] <- obj
  chroms <- c(chroms, obj$chrom)
}
message("[merge] ", length(chroms), " chromosomes in ",
        round((proc.time() - t_load)[3], 1), "s")

sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) {
    sample_names <- sub("^PC_1_", "", pc1_cols)
    break
  }
}
if (!is.null(sample_names)) message("[merge] Samples: ", length(sample_names))

# Load regions
region_files <- sort(list.files(outdir, pattern = "^seeded_regions_.*\\.rds$", full.names = TRUE))
if (length(region_files) == 0) {
  all_regions_file <- file.path(outdir, "seeded_regions_all.rds")
  if (!file.exists(all_regions_file)) {
    stop("No region files found. Run STEP_C01b_1 first.")
  }
  all_chr_regions <- readRDS(all_regions_file)
} else {
  all_chr_regions <- list()
  for (cf in region_files) {
    obj <- readRDS(cf)
    all_chr_regions[[obj$chrom]] <- obj
  }
}
message("[merge] Core data loaded for ", length(all_chr_regions), " chromosomes")

# =============================================================================
# LOAD LANDSCAPE DATA (v8.5 — from PHASE_01C)
# =============================================================================

landscape_conc   <- list()  # block_concordance per chr
landscape_cross  <- list()  # blue_cross_verdicts per chr
landscape_bound  <- list()  # boundary_catalog per chr
landscape_pa     <- list()  # 01C_window_pa per chr (block membership)
has_landscape    <- FALSE

if (!is.null(landscape_dir) && dir.exists(landscape_dir)) {
  message("[merge] Loading landscape data from: ", landscape_dir)

  # Load per-chromosome files
  for (chr in chroms) {
    f_conc  <- file.path(landscape_dir, paste0("block_concordance_", chr, ".tsv.gz"))
    f_cross <- file.path(landscape_dir, paste0("blue_cross_verdicts_", chr, ".tsv.gz"))
    f_bound <- file.path(landscape_dir, paste0("boundary_catalog_", chr, ".tsv.gz"))
    if (file.exists(f_conc))  landscape_conc[[chr]]  <- fread(f_conc)
    if (file.exists(f_cross)) landscape_cross[[chr]] <- fread(f_cross)
    if (file.exists(f_bound)) landscape_bound[[chr]] <- fread(f_bound)
  }

  # Load window PA (single file, all chromosomes)
  f_pa <- file.path(landscape_dir, "01C_window_pa.tsv.gz")
  if (file.exists(f_pa)) {
    pa_all <- fread(f_pa)
    # Fix column names: landscape may have pos_mb but merge needs start_bp
    if ("pos_mb" %in% names(pa_all) && !"start_bp" %in% names(pa_all)) {
      pa_all[, start_bp := as.integer(pos_mb * 1e6)]
      pa_all[, end_bp := start_bp]
    }
    for (chr in unique(pa_all$chrom)) landscape_pa[[chr]] <- pa_all[chrom == chr]
    rm(pa_all)
  }

  n_conc  <- sum(vapply(landscape_conc,  nrow, integer(1)))
  n_cross <- sum(vapply(landscape_cross, nrow, integer(1)))
  n_bound <- sum(vapply(landscape_bound, nrow, integer(1)))
  n_pa    <- sum(vapply(landscape_pa,    nrow, integer(1)))
  has_landscape <- (n_conc + n_cross + n_bound) > 0

  message("[merge] Landscape: concordance=", n_conc, " crosses=", n_cross,
          " boundaries=", n_bound, " pa_windows=", n_pa)
  if (!has_landscape) message("[merge] WARNING: landscape_dir provided but no data found")
} else {
  if (!is.null(landscape_dir)) message("[merge] WARNING: landscape_dir does not exist: ", landscape_dir)
  message("[merge] Running WITHOUT landscape integration (v8.4 mode)")
}

# =============================================================================
# MAIN MERGE LOOP
# =============================================================================

all_region_rows  <- list()
all_hier_rows    <- list()
all_window_rows  <- list()
all_diag_rows    <- list()
all_compare_rows <- list()
all_summary_rows <- list()
all_score_logs   <- list()
all_merge_pa     <- list()

region_id <- max(0L, unlist(lapply(all_chr_regions, function(x) {
  max(0L, vapply(x$regions, function(c) c$region_id, integer(1)))
})))

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  chr_regions_obj <- all_chr_regions[[chr]]
  if (is.null(pc) || is.null(chr_regions_obj)) { message("[SKIP] ", chr); next }

  dt <- pc$dt; sim_mat <- pc$sim_mat
  n <- nrow(dt)
  all_regions <- chr_regions_obj$regions

  message("\n[merge] ======= ", chr, " =======")
  message("[merge] ", length(all_regions), " regions from C01b_1")
  if (length(all_regions) == 0) next

  # ═══════════════════════════════════════════════════════════════════
  # RUN THREE MERGE FAMILIES (fuzzy composition)
  # ═══════════════════════════════════════════════════════════════════

  for (mfam in MERGE_FAMILIES) {
    result <- run_merge_fuzzy(all_regions, dt, sample_names, sim_mat, chr, mfam,
                               conc_dt_chr  = landscape_conc[[chr]],
                               cross_dt_chr = landscape_cross[[chr]],
                               bound_dt_chr = landscape_bound[[chr]],
                               block_pa_chr = landscape_pa[[chr]])
    merged <- result$regions

    if (length(result$score_log) > 0) {
      all_score_logs <- c(all_score_logs, result$score_log)
    }

    # Spread classes
    n_dense <- 0L; n_mod <- 0L; n_sparse <- 0L; n_tiny <- 0L
    for (reg in merged) {
      sp <- max(dt$end_bp[reg$idxs]) - min(dt$start_bp[reg$idxs])
      np <- sum(dt$start_bp >= min(dt$start_bp[reg$idxs]) &
                dt$end_bp <= max(dt$end_bp[reg$idxs]))
      dens <- if (np > 0) length(reg$idxs) / np else 0
      if (sp / 1e6 < 0.02) n_tiny <- n_tiny + 1L
      else if (dens >= 0.70) n_dense <- n_dense + 1L
      else if (dens >= 0.30) n_mod <- n_mod + 1L
      else n_sparse <- n_sparse + 1L
    }
    message("[merge] ", mfam$name, " (", mfam$label, "): ", length(merged),
            " regions | dense:", n_dense, " mod:", n_mod,
            " sparse:", n_sparse, " tiny:", n_tiny)

    # Store regions
    merged_ids <- integer(length(merged))
    for (mi in seq_along(merged)) {
      region_id <- region_id + 1L
      merged_ids[mi] <- region_id
      reg <- merged[[mi]]; idxs <- reg$idxs

      for (k in seq_along(idxs)) {
        all_window_rows[[length(all_window_rows) + 1]] <- data.table(
          chrom = chr, global_window_id = dt$global_window_id[idxs[k]],
          start_bp = dt$start_bp[idxs[k]], end_bp = dt$end_bp[idxs[k]],
          region_id = region_id, extension_phase = "merged",
          scale_tier = NA_character_, merge_family = mfam$name,
          inclusion_status = reg$statuses[k],
          continuity_score = if (is.finite(reg$scores[k])) round(reg$scores[k], 4) else NA_real_
        )
      }

      span_bp <- max(dt$end_bp[idxs]) - min(dt$start_bp[idxs])
      n_possible <- sum(dt$start_bp >= min(dt$start_bp[idxs]) &
                        dt$end_bp <= max(dt$end_bp[idxs]))
      density <- if (n_possible > 0) length(idxs) / n_possible else 0
      span_mb <- span_bp / 1e6
      mean_inv <- if ("inv_likeness" %in% names(dt))
        round(mean(dt$inv_likeness[idxs], na.rm = TRUE), 4) else NA_real_
      spread_class <- if (span_mb < 0.02) "tiny_fragment"
                      else if (density >= 0.70) "dense_continuous"
                      else if (density >= 0.30) "moderate_gaps"
                      else "sparse_scattered"

      coh <- region_coherence(idxs, sim_mat)

      all_region_rows[[length(all_region_rows) + 1]] <- data.table(
        region_id = region_id, extension_phase = "merged",
        scale_tier = NA_character_, merge_family = mfam$name, chrom = chr,
        start_bp = min(dt$start_bp[idxs]), end_bp = max(dt$end_bp[idxs]),
        n_windows = length(idxs), n_possible_windows = n_possible,
        density = round(density, 3), span_mb = round(span_mb, 3),
        spread_class = spread_class, mean_inv = mean_inv,
        n_seeds = sum(reg$statuses == "seed"),
        n_tolerated = sum(grepl("tolerate", reg$statuses)),
        n_bridge = sum(grepl("bridge", reg$statuses)),
        mean_score = round(mean(reg$scores[is.finite(reg$scores)], na.rm = TRUE), 4),
        min_score = round(min(reg$scores[is.finite(reg$scores)], na.rm = TRUE), 4),
        coherence = round(coh, 4)
      )
    }

    # Hierarchy
    for (mi in seq_along(merged)) {
      ms <- merged[[mi]]$idxs
      for (ci in seq_along(all_regions)) {
        ol <- length(intersect(all_regions[[ci]]$idxs, ms))
        if (ol > 0) {
          all_hier_rows[[length(all_hier_rows) + 1]] <- data.table(
            source_region_id = all_regions[[ci]]$region_id,
            scale_tier = all_regions[[ci]]$scale_tier,
            merged_snake_id = merged_ids[mi],
            merge_family = mfam$name, chrom = chr,
            overlap_windows = ol,
            overlap_frac = round(ol / length(all_regions[[ci]]$idxs), 4)
          )
        }
      }
    }

    # QC
    qc_dt <- run_qc(merged, sim_mat)
    if (nrow(qc_dt) > 0) all_diag_rows[[length(all_diag_rows) + 1]] <- qc_dt
  }

  # ═══════════════════════════════════════════════════════════════════
  # THREE-WAY COMPARISON
  # ═══════════════════════════════════════════════════════════════════

  reg_dt_chr <- rbindlist(all_region_rows[vapply(all_region_rows,
    function(x) is.data.table(x) && nrow(x) > 0 && x$chrom[1] == chr &&
                !is.na(x$extension_phase[1]) && x$extension_phase[1] == "merged",
    logical(1))], fill = TRUE)

  if (nrow(reg_dt_chr) > 0) {
    rA <- reg_dt_chr[merge_family == "merge_A"]
    rB <- reg_dt_chr[merge_family == "merge_B"]
    rC <- reg_dt_chr[merge_family == "merge_C"]

    for (ai in seq_len(nrow(rA))) {
      a_start <- rA$start_bp[ai]; a_end <- rA$end_bp[ai]
      b_ov <- if (nrow(rB) > 0) rB[start_bp <= a_end & end_bp >= a_start] else data.table()
      c_ov <- if (nrow(rC) > 0) rC[start_bp <= a_end & end_bp >= a_start] else data.table()
      n_b <- nrow(b_ov); n_c <- nrow(c_ov)

      agreement <- if (n_b >= 1 && n_c == 1) "all_agree"
                   else if (n_b >= 1 && n_c > 1) "AB_agree_C_split_COMPOSITE"
                   else if (n_b >= 1 && n_c == 0) "AB_agree_C_none"
                   else if (n_b == 0 && n_c >= 1) "A_only_C_exists"
                   else "A_only_weak"

      all_compare_rows[[length(all_compare_rows) + 1]] <- data.table(
        chrom = chr, merge_A_id = rA$region_id[ai],
        merge_A_start = a_start, merge_A_end = a_end,
        merge_A_windows = rA$n_windows[ai], merge_A_density = rA$density[ai],
        merge_A_coherence = rA$coherence[ai],
        n_merge_B = n_b, n_merge_C = n_c,
        merge_B_ids = if (n_b > 0) paste(b_ov$region_id, collapse = ";") else NA,
        merge_C_ids = if (n_c > 0) paste(c_ov$region_id, collapse = ";") else NA,
        agreement = agreement
      )
    }

    if (nrow(rA) > 0) {
      comp_chr <- rbindlist(tail(all_compare_rows, nrow(rA)))
      message("    [merge3] A:", nrow(rA), " B:", nrow(rB), " C:", nrow(rC),
              " | agree:", sum(comp_chr$agreement == "all_agree"),
              " composite:", sum(grepl("COMPOSITE", comp_chr$agreement)),
              " weak:", sum(comp_chr$agreement == "A_only_weak"))
    }
  }

  # Summary
  n_s <- sum(vapply(all_regions, function(r) r$scale_tier == "1S", logical(1)))
  n_m <- sum(vapply(all_regions, function(r) r$scale_tier == "1M", logical(1)))
  n_l <- sum(vapply(all_regions, function(r) r$scale_tier == "1L", logical(1)))
  all_summary_rows[[length(all_summary_rows) + 1]] <- data.table(
    chrom = chr, n_windows = n,
    n_regions_1S = n_s, n_regions_1M = n_m, n_regions_1L = n_l,
    n_regions_total = length(all_regions)
  )

  # =================================================================
  # MERGE ROARY-STYLE WINDOW PA MATRIX
  # =================================================================
  # Like the region PA but for merge regions. For each window:
  # - Which merge tier (A/B/C) claimed it?
  # - What merge region ID?
  # - What inclusion status (region_original / bridge_fuzzy)?
  # - What was the fuzzy membership score for its region?
  # Also loads C01b_1 window states to integrate gap info.

  # Build per-merge-tier lookup
  merge_lookup <- list()
  for (mfam_name in c("merge_A", "merge_B", "merge_C")) {
    fam_regions <- all_region_rows[vapply(all_region_rows, function(x)
      is.data.table(x) && nrow(x) > 0 && x$chrom[1] == chr &&
      !is.na(x$merge_family[1]) && x$merge_family[1] == mfam_name,
      logical(1))]
    wid_map <- list()
    # Also need the window rows for statuses
    fam_wins <- all_window_rows[vapply(all_window_rows, function(x)
      is.data.table(x) && nrow(x) > 0 && x$chrom[1] == chr &&
      !is.na(x$merge_family[1]) && x$merge_family[1] == mfam_name,
      logical(1))]
    if (length(fam_wins) > 0) {
      fw_dt <- rbindlist(fam_wins, fill = TRUE)
      for (ri in seq_len(nrow(fw_dt))) {
        wid_map[[as.character(fw_dt$global_window_id[ri])]] <- list(
          region_id = fw_dt$region_id[ri],
          status = fw_dt$inclusion_status[ri],
          score = fw_dt$continuity_score[ri]
        )
      }
    }
    merge_lookup[[mfam_name]] <- wid_map
  }

  # Load C01b_1 window states for gap info
  region_states_file <- file.path(outdir, paste0("seeded_regions_window_states_", chr, ".tsv.gz"))
  region_states <- if (file.exists(region_states_file)) fread(region_states_file) else data.table()

  # Build the merge PA matrix
  merge_pa_rows <- list()
  for (wi in seq_len(n)) {
    wid <- dt$global_window_id[wi]
    wid_str <- as.character(wid)

    a_info <- merge_lookup[["merge_A"]][[wid_str]]
    b_info <- merge_lookup[["merge_B"]][[wid_str]]
    c_info <- merge_lookup[["merge_C"]][[wid_str]]

    in_A <- !is.null(a_info); in_B <- !is.null(b_info); in_C <- !is.null(c_info)
    n_tiers <- sum(c(in_A, in_B, in_C))

    # Merge PA pattern
    mpa <- paste0(if (in_A) "A" else "", if (in_B) "B" else "", if (in_C) "C" else "")
    if (nchar(mpa) == 0) mpa <- "none"

    # Core-phase info
    region_row <- if (nrow(region_states) > 0) region_states[global_window_id == wid] else data.table()
    region_pa <- if (nrow(region_row) > 0) region_row$pa_pattern[1] else NA_character_
    region_state <- if (nrow(region_row) > 0) region_row$collector_state[1] else NA_character_
    is_seed_like <- if (nrow(region_row) > 0) region_row$is_seed_like[1] else FALSE
    n_region_scale_tiers <- if (nrow(region_row) > 0) region_row$n_families_claimed[1] else 0L

    # Determine role in merge
    merge_role <- if (n_tiers > 0) {
      # Check if it was a bridge or original region
      statuses <- c(
        if (in_A) a_info$status else NULL,
        if (in_B) b_info$status else NULL,
        if (in_C) c_info$status else NULL
      )
      if (any(grepl("bridge", statuses))) "bridge"
      else if (any(statuses %in% c("seed", "accepted"))) "region_window"
      else if (any(statuses == "tolerated")) "tolerated_window"
      else "merged_other"
    } else if (is_seed_like) "gap_with_signal"
    else "background"

    merge_pa_rows[[wi]] <- data.table(
      chrom = chr,
      global_window_id = wid,
      start_bp = dt$start_bp[wi],
      end_bp = dt$end_bp[wi],
      pos_mb = round((dt$start_bp[wi] + dt$end_bp[wi]) / 2e6, 4),
      # Core phase info
      region_pa = region_pa,
      region_state = region_state,
      n_core_families = n_region_scale_tiers,
      is_seed_like = is_seed_like,
      # Merge phase PA
      in_merge_A = in_A, in_merge_B = in_B, in_merge_C = in_C,
      region_A = if (in_A) a_info$region_id else NA_integer_,
      region_B = if (in_B) b_info$region_id else NA_integer_,
      region_C = if (in_C) c_info$region_id else NA_integer_,
      status_A = if (in_A) a_info$status else NA_character_,
      status_B = if (in_B) b_info$status else NA_character_,
      status_C = if (in_C) c_info$status else NA_character_,
      score_A = if (in_A && is.finite(a_info$score)) round(a_info$score, 4) else NA_real_,
      score_B = if (in_B && is.finite(b_info$score)) round(b_info$score, 4) else NA_real_,
      score_C = if (in_C && is.finite(c_info$score)) round(c_info$score, 4) else NA_real_,
      n_merge_tiers = n_tiers,
      merge_pa = mpa,
      merge_role = merge_role
    )
  }

  merge_pa_dt <- rbindlist(merge_pa_rows)
  all_merge_pa <- c(all_merge_pa, list(merge_pa_dt))

  # Log merge PA summary
  n_abc <- sum(merge_pa_dt$merge_pa == "ABC")
  n_ab <- sum(merge_pa_dt$merge_pa == "AB")
  n_a_only <- sum(merge_pa_dt$merge_pa == "A")
  n_bridge <- sum(merge_pa_dt$merge_role == "bridge")
  n_gap_sig <- sum(merge_pa_dt$merge_role == "gap_with_signal")
  message("[merge PA] ", chr, ": ABC=", n_abc, " AB=", n_ab, " A_only=", n_a_only,
          " bridge=", n_bridge, " gap_with_signal=", n_gap_sig)

  # =================================================================
  # SPLIT HEATMAP: upper = sim_mat, lower = merge PA
  # =================================================================
  if (requireNamespace("ggplot2", quietly = TRUE) && n >= 50) {
    library(ggplot2)
    dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

    step_p <- max(1L, n %/% 500)
    idx_p <- seq(1, n, by = step_p)
    ns <- length(idx_p)
    pos_mb <- (dt$start_bp[idx_p] + dt$end_bp[idx_p]) / 2e6

    hm_rows <- list()
    for (ii in seq_len(ns)) {
      for (jj in seq_len(ns)) {
        wi_i <- idx_p[ii]; wi_j <- idx_p[jj]
        if (ii <= jj) {
          v <- sim_mat[wi_i, wi_j]
          hm_rows[[length(hm_rows) + 1]] <- data.table(
            x = pos_mb[ii], y = pos_mb[jj],
            value = if (is.finite(v)) v else 0,
            layer = "sim_mat")
        }
        if (ii >= jj) {
          # Lower: merge PA
          # 0 = background, 0.2 = gap_with_signal, 0.4 = bridge,
          # 0.6 = A only, 0.8 = AB, 1.0 = ABC
          mp_i <- merge_pa_dt[global_window_id == dt$global_window_id[wi_i]]
          mp_j <- merge_pa_dt[global_window_id == dt$global_window_id[wi_j]]
          if (nrow(mp_i) > 0 && nrow(mp_j) > 0) {
            nt <- min(mp_i$n_merge_tiers, mp_j$n_merge_tiers)
            role_i <- mp_i$merge_role[1]; role_j <- mp_j$merge_role[1]
            is_bridge <- role_i == "bridge" | role_j == "bridge"
            is_gap <- role_i == "gap_with_signal" | role_j == "gap_with_signal"
            v <- if (nt == 3) 1.0 else if (nt == 2) 0.8
                 else if (nt == 1 && !is_bridge) 0.6
                 else if (is_bridge) 0.4 else if (is_gap) 0.2 else 0
          } else v <- 0
          hm_rows[[length(hm_rows) + 1]] <- data.table(
            x = pos_mb[ii], y = pos_mb[jj],
            value = v, layer = "merge_pa")
        }
      }
    }
    hm_dt <- rbindlist(hm_rows)

    p_split <- ggplot(hm_dt, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("grey90", "gold3", "cyan3", "steelblue", "darkorange", "red3"),
        name = "Value") +
      coord_fixed() +
      labs(title = paste0(chr, " -- Split Heatmap: simmat (upper) vs merge (lower)"),
           subtitle = paste0(n, " windows | merge_A:", sum(merge_pa_dt$in_merge_A),
                            " B:", sum(merge_pa_dt$in_merge_B),
                            " C:", sum(merge_pa_dt$in_merge_C)),
           x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
           caption = paste0("Upper: raw similarity | Lower: merge PA\n",
                           "red=ABC, orange=AB, blue=A, cyan=bridge, ",
                           "gold=gap_with_signal, grey=background")) +
      theme_minimal(base_size = 8) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.caption = element_text(size = 5, color = "grey60"))

    tryCatch(
      ggsave(file.path(outdir, "plots",
                         paste0(chr, "_A2_simmat_vs_merge.png")),
             p_split, width = 12, height = 12, dpi = 300),
      error = function(e) message("  [PLOT] ", e$message))
    message("[merge] Split heatmap: ", chr)
  }
}

# =============================================================================
# WRITE
# =============================================================================

empty <- function(...) data.table()
reg_dt    <- if (length(all_region_rows) > 0) rbindlist(all_region_rows, fill = TRUE) else empty()
hier_dt   <- if (length(all_hier_rows) > 0) rbindlist(all_hier_rows, fill = TRUE) else empty()
win_dt    <- if (length(all_window_rows) > 0) rbindlist(all_window_rows, fill = TRUE) else empty()
diag_dt   <- if (length(all_diag_rows) > 0) rbindlist(all_diag_rows, fill = TRUE) else empty()
summ_dt   <- if (length(all_summary_rows) > 0) rbindlist(all_summary_rows) else empty()
comp_dt   <- if (length(all_compare_rows) > 0) rbindlist(all_compare_rows, fill = TRUE) else empty()
scores_dt <- if (length(all_score_logs) > 0) rbindlist(all_score_logs, fill = TRUE) else empty()
merge_pa_all <- if (length(all_merge_pa) > 0) rbindlist(all_merge_pa, fill = TRUE) else empty()

files <- c(
  region_merge_regions    = "region_merge_regions.tsv.gz",
  region_merge_hierarchy  = "region_merge_hierarchy.tsv.gz",
  region_merge_windows    = "region_merge_windows.tsv.gz",
  region_merge_qc         = "region_merge_qc.tsv.gz",
  region_merge_summary    = "region_merge_summary.tsv",
  region_merge_comparison = "region_merge_comparison.tsv.gz",
  region_merge_scores     = "region_merge_scores.tsv.gz",
  region_merge_pa         = "region_merge_window_pa.tsv.gz"
)

dts <- list(reg_dt, hier_dt, win_dt, diag_dt, summ_dt, comp_dt, scores_dt, merge_pa_all)
for (i in seq_along(files)) {
  fwrite(dts[[i]], file.path(outdir, files[i]), sep = "\t")
}

message("\n[DONE] STEP_C01b_2_merge (fuzzy composition) complete")
for (f in files) message("  ", file.path(outdir, f))

# Score log summary
if (nrow(scores_dt) > 0) {
  for (mf in unique(scores_dt$merge_family)) {
    mf_dt <- scores_dt[merge_family == mf]
    land_adj_col <- if ("landscape_adj" %in% names(mf_dt)) mf_dt$landscape_adj else rep(0, nrow(mf_dt))
    n_adjusted <- sum(abs(land_adj_col) > 0.001, na.rm = TRUE)
    message("[scores] ", mf, ": ", nrow(mf_dt), " pairs, ",
            sum(mf_dt$decision == "merge"), " merged, ",
            sum(mf_dt$decision == "reject"), " rejected | ",
            "mean fuzzy=", round(mean(mf_dt$fuzzy_score, na.rm = TRUE), 3),
            " memb=", round(mean(mf_dt$membership_score, na.rm = TRUE), 3),
            " geom=", round(mean(mf_dt$geometric_score, na.rm = TRUE), 3),
            " | landscape_adjusted=", n_adjusted)
  }
}

# =============================================================================
# EXPORT: region_candidate_regions.tsv.gz
# =============================================================================

region_cand_rows <- list()
scid <- 0L
if (nrow(reg_dt) > 0) {
  merge_a <- reg_dt[extension_phase == "merged" & merge_family == "merge_A"]
  if (nrow(merge_a) > 0) {
    for (ri in seq_len(nrow(merge_a))) {
      scid <- scid + 1L; r <- merge_a[ri]
      region_cand_rows[[scid]] <- data.table(
        candidate_id = scid, source = "region_merge_A",
        region_id = r$region_id, chrom = r$chrom,
        start_bp = r$start_bp, end_bp = r$end_bp,
        n_windows = r$n_windows, density = r$density,
        span_mb = r$span_mb, spread_class = r$spread_class,
        mean_inv = r$mean_inv, coherence = r$coherence,
        merge_family = "merge_A"
      )
    }
  }
}

region_cand_dt <- if (length(region_cand_rows) > 0) rbindlist(region_cand_rows, fill = TRUE) else data.table()
fwrite(region_cand_dt, file.path(outdir, "region_candidate_regions.tsv.gz"), sep = "\t")
message("[merge] Candidates: ", nrow(region_cand_dt), " -> ", file.path(outdir, "region_candidate_regions.tsv.gz"))

# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })
plot_dir <- file.path(outdir, "merge_diagnostics")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Plot 1: Fuzzy score distributions per merge family ───────────────
# Shows the distribution of fuzzy scores for merged vs rejected pairs.
# Helps calibrate accept_threshold.

if (nrow(scores_dt) > 0) {

  g1 <- ggplot(scores_dt, aes(x = fuzzy_score, fill = decision)) +
    geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.6) +
    facet_wrap(~ merge_family, ncol = 1, scales = "free_y") +
    geom_vline(data = data.table(
      merge_family = c("merge_A", "merge_B", "merge_C"),
      thresh = c(MA$accept_threshold, MB$accept_threshold, MC$accept_threshold)
    ), aes(xintercept = thresh), linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = c(merge = "#2563eb", reject = "#dc2626")) +
    labs(title = "Fuzzy merge score distribution",
         subtitle = "Dashed line = accept threshold",
         x = "Fuzzy score", y = "Region pairs") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")

  ggsave(file.path(plot_dir, "01_fuzzy_score_distribution.pdf"),
         g1, width = 8, height = 10)

  # ── Plot 2: Membership vs Geometric scatter ────────────────────────
  # Shows the two relations against each other. Diagonal = both agree.
  # Off-diagonal = one relation compensates for the other.

  g2 <- ggplot(scores_dt[is.finite(membership_score) & is.finite(geometric_score)],
               aes(x = geometric_score, y = membership_score, color = decision)) +
    geom_point(size = 0.6, alpha = 0.5) +
    facet_wrap(~ merge_family, ncol = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c(merge = "#2563eb", reject = "#dc2626")) +
    labs(title = "Membership vs Geometric scores",
         subtitle = "Diagonal = equal contribution from both relations",
         x = "Geometric score (sim_mat continuity)",
         y = "Membership score (fuzzy band composition)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")

  ggsave(file.path(plot_dir, "02_membership_vs_geometric.pdf"),
         g2, width = 12, height = 5)

  # ── Plot 3: Score vs gap size ──────────────────────────────────────
  # Shows how fuzzy score decays with gap distance.

  g3 <- ggplot(scores_dt, aes(x = gap_windows, y = fuzzy_score, color = decision)) +
    geom_jitter(size = 0.5, alpha = 0.4, width = 0.3) +
    facet_wrap(~ merge_family, ncol = 3) +
    geom_hline(data = data.table(
      merge_family = c("merge_A", "merge_B", "merge_C"),
      thresh = c(MA$accept_threshold, MB$accept_threshold, MC$accept_threshold)
    ), aes(yintercept = thresh), linetype = "dashed", color = "red") +
    scale_color_manual(values = c(merge = "#2563eb", reject = "#dc2626")) +
    labs(title = "Fuzzy score vs gap size",
         x = "Gap (windows between regions)", y = "Fuzzy merge score") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")

  ggsave(file.path(plot_dir, "03_score_vs_gap.pdf"),
         g3, width = 12, height = 5)

  message("[plots] Score diagnostics: 3 PDFs -> ", plot_dir)
}

# ── Plot 4: Per-chromosome merge ideogram ────────────────────────────
# One page per chromosome: shows where merge_A/B/C regions land,
# with density shading and spread class coloring.

if (nrow(reg_dt) > 0) {

  spread_pal <- c(
    dense_continuous = "#16a34a",
    moderate_gaps    = "#d97706",
    sparse_scattered = "#dc2626",
    tiny_fragment    = "#9ca3af"
  )

  for (chr in unique(reg_dt$chrom)) {
    chr_reg <- reg_dt[chrom == chr & extension_phase == "merged"]
    if (nrow(chr_reg) == 0) next

    chr_len <- max(precomp_list[[chr]]$dt$end_bp, na.rm = TRUE)

    # Build rectangle data for each merge family
    chr_reg[, y_min := fifelse(merge_family == "merge_A", 2,
                       fifelse(merge_family == "merge_B", 1, 0))]
    chr_reg[, y_max := y_min + 0.8]
    chr_reg[, start_mb := start_bp / 1e6]
    chr_reg[, end_mb := end_bp / 1e6]

    # Ensure spread_class is a factor for proper coloring
    chr_reg[, spread_class := factor(spread_class,
      levels = c("dense_continuous", "moderate_gaps", "sparse_scattered", "tiny_fragment"))]

    g4 <- ggplot(chr_reg) +
      geom_rect(aes(xmin = start_mb, xmax = end_mb,
                     ymin = y_min, ymax = y_max,
                     fill = spread_class),
                color = "gray30", linewidth = 0.15, alpha = 0.7) +
      scale_fill_manual(values = spread_pal, drop = FALSE) +
      scale_y_continuous(
        breaks = c(0.4, 1.4, 2.4),
        labels = c("C (strict)", "B (moderate)", "A (generous)"),
        limits = c(-0.2, 3.2)
      ) +
      labs(title = paste0(chr, " — Fuzzy merge regions (", nrow(chr_reg), " total)"),
           subtitle = paste0("Chr length: ", round(chr_len / 1e6, 1), " Mb"),
           x = "Position (Mb)", y = "Merge tier", fill = "Spread class") +
      theme_minimal(base_size = 9) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )

    # Add comparison annotations if available
    chr_comp <- comp_dt[chrom == chr]
    if (nrow(chr_comp) > 0) {
      n_agree <- sum(chr_comp$agreement == "all_agree")
      n_comp  <- sum(grepl("COMPOSITE", chr_comp$agreement))
      n_weak  <- sum(chr_comp$agreement == "A_only_weak")
      g4 <- g4 + labs(caption = paste0(
        "agree:", n_agree, " | composite:", n_comp, " | weak:", n_weak))
    }

    ggsave(file.path(plot_dir, paste0("04_ideogram_", chr, ".pdf")),
           g4, width = 14, height = 3.5)
  }
  message("[plots] Ideograms: ", length(unique(reg_dt$chrom)), " PDFs -> ", plot_dir)
}

# ── Plot 5: Fuzzy composition T-matrix heatmaps (selected regions) ───
# For the top merge_A regions by span, show the 3×3 T matrix as a
# heatmap. Diagonal-dominant = same inversion, diffuse = composite.

if (nrow(reg_dt) > 0 && !is.null(sample_names)) {

  merge_a_big <- reg_dt[merge_family == "merge_A" & n_windows >= 20]
  merge_a_big <- merge_a_big[order(-span_mb)]
  n_to_plot <- min(20L, nrow(merge_a_big))

  if (n_to_plot > 0) {
    tmat_plots <- list()

    for (pi in seq_len(n_to_plot)) {
      r <- merge_a_big[pi]
      chr <- r$chrom
      pc <- precomp_list[[chr]]
      if (is.null(pc)) next

      dt_chr <- pc$dt

      # Get window indices for this region
      reg_idxs <- which(dt_chr$start_bp >= r$start_bp & dt_chr$end_bp <= r$end_bp)
      if (length(reg_idxs) < 10) next

      # Split into first half / second half and compute T matrix
      mid <- length(reg_idxs) %/% 2
      half_a <- reg_idxs[1:mid]
      half_b <- reg_idxs[(mid + 1):length(reg_idxs)]

      memb_a <- soft_band_membership(dt_chr, half_a, sample_names)
      memb_b <- soft_band_membership(dt_chr, half_b, sample_names)

      if (is.null(memb_a) || is.null(memb_b)) next

      T_mat <- fuzzy_compose(memb_a, memb_b)
      if (is.null(T_mat)) next

      # Build heatmap data
      t_score <- composition_score(T_mat)
      hm_dt <- data.table(
        band_a = rep(paste0("A:", colnames(T_mat)), each = 3),
        band_b = rep(paste0("B:", colnames(T_mat)), 3),
        value = as.vector(T_mat)
      )

      tmat_plots[[pi]] <- ggplot(hm_dt, aes(x = band_b, y = band_a, fill = value)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = round(value, 2)), size = 3.5) +
        scale_fill_gradient2(low = "#f1f5f9", mid = "#93c5fd", high = "#1e40af",
                             midpoint = 0.5, limits = c(0, 1)) +
        labs(title = paste0(r$chrom, " ", round(r$start_bp/1e6, 1), "-",
                            round(r$end_bp/1e6, 1), " Mb (",
                            r$n_windows, " win, T=", round(t_score, 2), ")"),
             x = "Second half", y = "First half") +
        theme_minimal(base_size = 8) +
        theme(legend.position = "none", axis.text = element_text(size = 7))
    }

    tmat_plots <- tmat_plots[!vapply(tmat_plots, is.null, logical(1))]
    if (length(tmat_plots) > 0) {
      # Arrange in grid
      n_plots <- length(tmat_plots)
      ncol_g <- min(4L, n_plots)
      nrow_g <- ceiling(n_plots / ncol_g)

      g5 <- do.call(gridExtra::arrangeGrob,
                     c(tmat_plots, ncol = ncol_g,
                       top = "Fuzzy composition T matrices (first half vs second half)"))

      ggsave(file.path(plot_dir, "05_T_matrix_heatmaps.pdf"),
             g5, width = ncol_g * 3.5, height = nrow_g * 3)
      message("[plots] T-matrix heatmaps: ", length(tmat_plots), " regions -> ",
              file.path(plot_dir, "05_T_matrix_heatmaps.pdf"))
    }
  }
}

# ── Plot 6: Genome-wide merge summary ───────────────────────────────
# Barplot: per-chromosome count of merge_A regions, colored by spread class.

if (nrow(reg_dt) > 0) {
  merge_a_all <- reg_dt[merge_family == "merge_A"]
  if (nrow(merge_a_all) > 0) {
    merge_a_all[, spread_class := factor(spread_class,
      levels = c("dense_continuous", "moderate_gaps", "sparse_scattered", "tiny_fragment"))]

    # Order chromosomes by name
    chr_order <- unique(merge_a_all$chrom)
    chr_order <- chr_order[order(as.integer(gsub("\\D", "", chr_order)))]
    merge_a_all[, chrom := factor(chrom, levels = chr_order)]

    g6 <- ggplot(merge_a_all, aes(x = chrom, fill = spread_class)) +
      geom_bar() +
      scale_fill_manual(values = c(
        dense_continuous = "#16a34a", moderate_gaps = "#d97706",
        sparse_scattered = "#dc2626", tiny_fragment = "#9ca3af"
      ), drop = FALSE) +
      labs(title = "Merge A regions per chromosome",
           subtitle = paste0("Total: ", nrow(merge_a_all), " regions across ",
                             length(unique(merge_a_all$chrom)), " chromosomes"),
           x = "Chromosome", y = "Merged regions", fill = "Spread class") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")

    ggsave(file.path(plot_dir, "06_genome_merge_summary.pdf"),
           g6, width = 12, height = 5)
    message("[plots] Genome summary -> ", file.path(plot_dir, "06_genome_merge_summary.pdf"))
  }
}

message("\n[DONE] All merge diagnostics complete")
