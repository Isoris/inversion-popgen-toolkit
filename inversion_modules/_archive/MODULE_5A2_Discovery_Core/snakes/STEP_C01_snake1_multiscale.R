#!/usr/bin/env Rscript

# =============================================================================
# STEP10e_v3_multiscale_snake.R  (v8.3)
#
# MULTI-SCALE SNAKE 1 with 3 core families and 2 merge families.
#
# CORE FAMILIES (3 Snake-1 variants):
#   Snake 1S — STRICT core
#     High accept threshold (0.60), 0 gap, 0 bad budget, min 4 windows
#     seed_z_min = 2.5 (strict but reachable with chunked_2x)
#     Produces: safest possible nuclei
#
#   Snake 1M — MODERATE body
#     Standard thresholds (0.50 accept, 0.28 tolerate), gap=1, bad=1
#     seed_z_min = 1.8
#     Produces: plausible fuller bodies including shoulders
#
#   Snake 1L — BROAD exploratory body
#     Lower thresholds (0.40 accept, 0.22 tolerate), gap=2, bad=2
#     seed_z_min = 1.2 (exploratory, top ~11% of distribution)
#     Longer rolling memory (size=5)
#     Produces: diffuse but structured extensions
#
# MERGE FAMILIES (2 Snake-2 variants):
#   Merge A — CONSERVATIVE stitcher
#     Short gap (3 windows), high bridge coherence (0.40), strict coherence floor
#     Produces: high-confidence merged candidates
#
#   Merge B — PERMISSIVE envelope stitcher
#     Larger gap (8 windows), softer bridge threshold (0.25), lower coherence floor
#     Produces: broader candidate envelopes, flagged as less certain
#
# WHY THIS HELPS:
#   - Snake 1 becomes a multi-scale inventory, not just "tiny islands"
#   - Merge A gives high-confidence view, Merge B gives exploratory view
#   - If Merge A says "split" and Merge B says "connected", the region is
#     merge-sensitive — that's informative, not a failure
#   - Reduces pressure on downstream Snake 2/3 to reconstruct from dust
#

# EXECUTION ORDER:
#   1. Snake 1S on all windows (strictest, goes first, claims best seeds)
#   2. Snake 1M on unclaimed windows (moderate, fills shoulders)
#   3. Snake 1L on remaining unclaimed (broadest, captures diffuse structure)
#   4. Merge A on all cores from S+M+L (conservative stitching)
#   5. Merge B on all cores from S+M+L (permissive stitching)
#   6. QC on both merge outputs
#
# OUTPUTS (extends v2):
#   snake_windows.tsv.gz       — per-window with core_family (1S/1M/1L) label
#   snake_regions.tsv.gz       — per-region with core_family + merge_family
#   snake_hierarchy.tsv.gz     — core→merged relationships
#   snake_summary.tsv          — per-chromosome multi-family counts
#   snake_diagnostics.tsv.gz   — QC for both merge families
#   snake_decision_log.tsv.gz  — every extension attempt
#   snake_window_states.tsv.gz — per-window collector states
#   snake_multiscale_comparison.tsv.gz — merge A vs B comparison
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP10e_v3_multiscale_snake.R <step10_outprefix> <outdir>")
}

step10_prefix <- args[1]
outdir        <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PARAMETER SETS — 3 core families + 2 merge families
# =============================================================================

# Shared seed parameters
SEED_MDS_AXES    <- 5L
SEED_MIN_Z       <- 1.2
SEED_NEIGHBOR_K  <- 3L
SEED_MAX_NN_DIST <- 0.80

# Continuity weights (shared across all families)
W_PREV <- 0.40; W_ROLL <- 0.40; W_MDS <- 0.20

# ── Snake 1S: STRICT ─────────────────────────────────────────────────
S1S <- list(
  name = "1S", label = "strict",
  max_gap = 0L, max_bad = 0L, min_windows = 4L,
  accept = 0.60, tolerate = 0.35, roll_size = 3L,
  seed_z_min = 2.5,          # strict but reachable with chunked_2x
  seed_nn_max = 0.70         # tighter NN distance
)

# ── Snake 1M: MODERATE ───────────────────────────────────────────────
S1M <- list(
  name = "1M", label = "moderate",
  max_gap = 1L, max_bad = 1L, min_windows = 3L,
  accept = 0.50, tolerate = 0.28, roll_size = 3L,
  seed_z_min = 1.8,          # moderate: fills gap between S1S(2.5) and S1L(1.2)
  seed_nn_max = 0.80
)

# ── Snake 1L: BROAD ─────────────────────────────────────────────────
S1L <- list(
  name = "1L", label = "broad",
  max_gap = 2L, max_bad = 2L, min_windows = 3L,
  accept = 0.40, tolerate = 0.22, roll_size = 5L,
  seed_z_min = 1.2,          # exploratory: top ~11% of distribution
  seed_nn_max = 0.90
)

CORE_FAMILIES <- list(S1S, S1M, S1L)

# ── Merge A: CONSERVATIVE ────────────────────────────────────────────
MA <- list(
  name = "merge_A", label = "conservative",
  max_gap_between = 3L,
  accept = 0.45, tolerate = 0.30,
  max_bad_bridge = 1L,
  coherence_min = 0.40,
  roll_size = 4L, min_windows = 5L,
  min_sample_jaccard = 0.40   # strict: outlier sets must overlap ≥40%
)

# ── Merge B: PERMISSIVE ─────────────────────────────────────────────
MB <- list(
  name = "merge_B", label = "permissive",
  max_gap_between = 8L,
  accept = 0.30, tolerate = 0.20,
  max_bad_bridge = 4L,
  coherence_min = 0.22,
  roll_size = 6L, min_windows = 4L,
  min_sample_jaccard = 0.20   # permissive: still requires SOME overlap
)

MERGE_FAMILIES <- list(MA, MB)

# QC parameters (shared)
QC_BOUNDARY_WEAK  <- 0.35
QC_INTERNAL_DROP  <- 0.25
QC_DROP_FRAC      <- 0.15
QC_SPLIT_THRESH   <- 0.20

# =============================================================================
# LOAD
# =============================================================================

mds_rds_file <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds_file)) stop("Missing: ", mds_rds_file)
mds_obj <- readRDS(mds_rds_file)
per_chr <- mds_obj$per_chr
if (is.null(per_chr) || length(per_chr) == 0) stop("No per_chr data")
chroms <- names(per_chr)
message("[v3] Multi-scale snake: ", length(chroms), " chromosomes")

# Extract sample names from first chromosome's data (PC columns)
sample_names_snake <- NULL
for (chr_tmp in chroms) {
  dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
  pc1_cols <- grep("^PC_1_", names(dt_tmp), value = TRUE)
  if (length(pc1_cols) > 0) {
    sample_names_snake <- sub("^PC_1_", "", pc1_cols)
    break
  }
}
if (!is.null(sample_names_snake)) {
  message("[v3] Sample names extracted: ", length(sample_names_snake), " samples")
} else {
  message("[v3] WARNING: No PC_1_* columns found — sample-identity gate disabled")
}

# =============================================================================
# INVERSION-LIKENESS SCORE (per-window, no MDS, no background)
#
# Computed directly from local PCA eigenvalues + PC1 score distribution.
# Three components:
#   1. pve1 = λ1 / (λ1 + λ2) — how dominant is PC1?
#   2. eigen_ratio = λ1 / max(λ2, ε) — how much stronger is PC1 vs PC2?
#   3. trimodality = dip test p-value on PC1 scores (low p = multimodal = inv-like)
#
# Combined into a single [0,1] score. Windows scoring high are inversion-like
# regardless of their MDS z-score.
#
# TWO USES:
#   A. As a parallel seed criterion — seed from inv-like windows even if z is low
#   B. As a background filter — exclude inv-like windows from background pool
#      so chunked_2x doesn't accidentally sample other chromosomes' inversions
# =============================================================================

compute_inv_likeness_all <- function(per_chr, chroms, sample_names) {
  all_rows <- list()

  for (chr in chroms) {
    dt <- as.data.table(per_chr[[chr]]$out_dt)
    if (nrow(dt) == 0) next

    # Eigenvalues
    has_lam <- "lam_1" %in% names(dt) && "lam_2" %in% names(dt)
    has_total <- "total" %in% names(dt)

    # PC1 scores
    pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
    has_pc1 <- length(pc1_cols) > 0

    for (wi in seq_len(nrow(dt))) {
      wid <- dt$global_window_id[wi]

      # Component 1: PVE1
      pve1 <- NA_real_
      eigen_ratio <- NA_real_
      if (has_lam) {
        l1 <- dt$lam_1[wi]; l2 <- dt$lam_2[wi]
        if (is.finite(l1) && is.finite(l2) && (l1 + l2) > 0) {
          pve1 <- l1 / (l1 + l2)
          eigen_ratio <- l1 / max(l2, 1e-10)
        }
      }

      # Component 2: Trimodality (Hartigan's dip test on PC1 scores)
      dip_p <- NA_real_
      if (has_pc1) {
        scores <- as.numeric(dt[wi, ..pc1_cols])
        scores <- scores[is.finite(scores)]
        if (length(scores) >= 20) {
          dip_p <- tryCatch({
            # diptest::dip.test — if not available, use gap heuristic
            if (requireNamespace("diptest", quietly = TRUE)) {
              diptest::dip.test(scores)$p.value
            } else {
              # Fallback: gap ratio between k=3 clusters
              km <- tryCatch(kmeans(scores, centers = 3, nstart = 3), error = function(e) NULL)
              if (!is.null(km)) {
                centers <- sort(km$centers[, 1])
                gaps <- diff(centers)
                within_sd <- sqrt(mean(km$withinss / km$size))
                if (within_sd > 0) {
                  # High gap/within_sd ratio = multimodal
                  gap_ratio <- min(gaps) / within_sd
                  # Convert to pseudo-p: high gap_ratio → low p
                  1 / (1 + gap_ratio^2)
                } else 1
              } else 1
            }
          }, error = function(e) NA_real_)
        }
      }

      # Combined score [0, 1]
      # pve1: rescale [0.3, 0.8] → [0, 1] (below 0.3 = no signal, above 0.8 = saturated)
      s_pve <- if (is.finite(pve1)) pmin(1, pmax(0, (pve1 - 0.3) / 0.5)) else 0
      # eigen_ratio: rescale [1, 10] → [0, 1]
      s_eig <- if (is.finite(eigen_ratio)) pmin(1, pmax(0, (eigen_ratio - 1) / 9)) else 0
      # dip_p: low p = multimodal = inv-like. Rescale p [0.05, 0.5] → [1, 0]
      s_dip <- if (is.finite(dip_p)) pmin(1, pmax(0, (0.5 - dip_p) / 0.45)) else 0

      inv_like <- 0.35 * s_pve + 0.25 * s_eig + 0.40 * s_dip

      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        inv_pve1 = round(pve1, 4),
        inv_eigen_ratio = round(eigen_ratio, 2),
        inv_dip_p = round(dip_p, 4),
        inv_likeness = round(inv_like, 4)
      )
    }
  }
  if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
}

message("[v3] Computing per-window inversion-likeness scores...")
inv_like_dt <- compute_inv_likeness_all(per_chr, chroms, sample_names_snake)
if (nrow(inv_like_dt) > 0) {
  n_high <- sum(inv_like_dt$inv_likeness >= 0.90, na.rm = TRUE)
  message("[v3] Inversion-likeness: ", nrow(inv_like_dt), " windows scored, ",
          n_high, " with score >= 0.5")
}

# =============================================================================
# HELPERS (identical to v2 — robust normalization)
# =============================================================================

make_sim_mat <- function(dmat) {
  finite_vals <- dmat[is.finite(dmat)]
  dmax <- if (length(finite_vals) > 0) quantile(finite_vals, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

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

prepare_chr <- function(chr_obj) {
  dt <- as.data.table(chr_obj$out_dt)
  dt <- dt[order(start_bp)]
  dmat <- chr_obj$dmat
  mds_cols <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  mds_mat <- as.matrix(dt[, ..mds_cols])

  # Merge inv_likeness scores (computed globally before snake loop)
  if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
    dt <- merge(dt, inv_like_dt[, .(global_window_id, inv_likeness)],
                by = "global_window_id", all.x = TRUE)
    dt <- dt[order(start_bp)]
  }

  z_cols <- grep("^MDS[0-9]+_z$", names(dt), value = TRUE)

  # ── ROBUST Z RECOMPUTATION ──────────────────────────────────────────
  # Recompute z-scores from raw MDS coordinates using median/MAD instead
  # of the mean/SD values stored in the .mds.rds. This prevents inversion-
  # heavy chromosomes from compressing their own z-scores.
  # The raw MDS coordinates are in mds_cols; we overwrite the _z columns.
  for (mc in mds_cols) {
    zc <- paste0(mc, "_z")
    vv <- dt[[mc]]
    med <- median(vv, na.rm = TRUE)
    mad_val <- mad(vv, na.rm = TRUE)
    if (is.finite(mad_val) && mad_val > 1e-10) {
      dt[[zc]] <- (vv - med) / mad_val
    } else {
      sdev <- sd(vv, na.rm = TRUE)
      if (is.finite(sdev) && sdev > 0) {
        dt[[zc]] <- (vv - mean(vv, na.rm = TRUE)) / sdev
      } else {
        dt[[zc]] <- 0
      }
    }
  }
  z_cols <- grep("^MDS[0-9]+_z$", names(dt), value = TRUE)

  if (length(z_cols) > SEED_MDS_AXES) z_cols <- z_cols[seq_len(SEED_MDS_AXES)]
  if (length(z_cols) > 0) {
    dt[, max_abs_z := apply(.SD, 1, function(x) max(abs(x), na.rm = TRUE)), .SDcols = z_cols]
  } else dt[, max_abs_z := 0]

  o_cols <- grep("^MDS[0-9]+_outlier$", names(dt), value = TRUE)
  if (length(o_cols) > SEED_MDS_AXES) o_cols <- o_cols[seq_len(SEED_MDS_AXES)]
  if (length(o_cols) > 0) {
    dt[, any_outlier := apply(.SD, 1, any, na.rm = TRUE), .SDcols = o_cols]
  } else dt[, any_outlier := FALSE]

  n_dt <- nrow(dt); n_dm <- nrow(dmat)
  if (n_dm != n_dt) {
    n <- min(n_dm, n_dt)
    dt <- dt[seq_len(n)]; dmat <- dmat[seq_len(n), seq_len(n), drop = FALSE]
    mds_mat <- mds_mat[seq_len(n), , drop = FALSE]
  }
  sim_mat <- make_sim_mat(dmat)

  nn_dists <- vapply(seq_len(nrow(dmat)), function(i) {
    d <- dmat[i, ]; d[i] <- Inf; d <- d[is.finite(d)]
    if (length(d) == 0) return(Inf)
    mean(sort(d)[seq_len(min(SEED_NEIGHBOR_K, length(d)))], na.rm = TRUE)
  }, numeric(1))
  dmax <- quantile(dmat[is.finite(dmat)], 0.95, na.rm = TRUE)
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  dt[, seed_nn_dist := nn_dists / dmax]

  list(dt = dt, dmat = dmat, sim_mat = sim_mat, mds_mat = mds_mat)
}

# =============================================================================
# DECISION LOG
# =============================================================================

decision_log <- list()

log_decision <- function(chr, sid, phase, direction, from_idx, to_idx,
                         from_wid, to_wid, score, accept_thresh, toler_thresh,
                         decision, reason, gap_count, bad_count, roll_size,
                         family_name) {
  decision_log[[length(decision_log) + 1L]] <<- data.table(
    chrom = chr, snake_id = sid, phase = phase, direction = direction,
    core_family = family_name,
    from_window_idx = from_idx, to_window_idx = to_idx,
    from_window_id = from_wid, to_window_id = to_wid,
    continuity_score = round(score, 4),
    accept_thresh = accept_thresh, toler_thresh = toler_thresh,
    decision = decision, reason = reason,
    gap_count = gap_count, bad_count = bad_count, roll_size = roll_size
  )
}

# =============================================================================
# PARAMETERIZED CORE ENGINE
# =============================================================================

run_core_family <- function(dt, sim_mat, mds_mat, chr, snake_id_start,
                            used_global, params) {
  # NOTE on used_global:
  # This is a READ-ONLY reference to what previous families claimed.
  # Each family builds its OWN membership. Windows claimed by earlier families
  # are NOT blocked — they can be re-explored. Deduplication happens LATER
  # when building the final core list.
  #
  # Why: shared used[] prevents true S/M/L nesting from the same region.
  # With independent membership, 1S finds a tight core, 1M finds a fuller
  # body overlapping that core, 1L finds a broad envelope. All three are
  # valid views of the same biological structure.

  n <- nrow(dt)
  if (n < params$min_windows) return(list(regions = list(), next_id = snake_id_start))

  # Per-family seed eligibility
  # NOTE (v8.2 fix): removed `any_outlier &` gate.
  # any_outlier is precomputed at Z_THRESH=3.0 in STEP10_v2, which was
  # clamping all families to z≥3.0 even when seed_z_min was set lower.
  # This made 1M (seed_z_min=2.5) and 1L (seed_z_min=2.0) effectively
  # identical to 1S for seeding. Now each family uses its own z threshold.
  #
  # v8.4: inv_likeness as PARALLEL seed gate. A window can be a seed if:
  #   (a) z >= seed_z_min (classic MDS outlier), OR
  #   (b) inv_likeness >= 0.90 (locally inversion-like from PCA eigenvalues)
  # Both still require reasonable NN distance.
  has_inv <- "inv_likeness" %in% names(dt)
  dt[, fam_seed := (max_abs_z >= params$seed_z_min |
                    (has_inv & is.finite(inv_likeness) & inv_likeness >= 0.90)) &
       is.finite(seed_nn_dist) &
       seed_nn_dist < params$seed_nn_max]

  # Order seeds: prefer high z, but inv-like windows with low z still enter
  dt[, seed_priority := max_abs_z]
  if (has_inv) dt[, seed_priority := pmax(seed_priority, inv_likeness * 5, na.rm = TRUE)]
  seed_order <- order(-dt$seed_priority)
  regions <- list()
  sid <- snake_id_start
  claimed <- rep(FALSE, n)  # this family's own claims (not shared)

  # ── DAMAGE MODEL PARAMETERS ────────────────────────────────────────
  # Cumulative damage replaces simple gap counter.
  # Good windows heal. Bad windows wound. Gaps wound harder.
  # Consecutive gaps wound fastest. Snake dies when damage persists.
  DMG_RECOVER_ACCEPTED  <- params$dmg_recover  %||% -0.05  # heal on good window
  DMG_COST_TOLERATED    <- params$dmg_tolerate %||%  0.04  # mild wound
  DMG_COST_GAP          <- params$dmg_gap      %||%  0.08  # wound per gap window
  DMG_COST_CONSEC_GAP   <- params$dmg_consec   %||%  0.05  # extra per consecutive gap
  DMG_MAX               <- params$dmg_max      %||%  0.30  # death threshold

  for (si in seed_order) {
    if (claimed[si] || !dt$fam_seed[si]) next

    sid <- sid + 1L
    snake <- si; status <- "seed"; scores <- 1.0
    claimed[si] <- TRUE

    # ── EXTEND RIGHT with cumulative damage ──────────────────────────
    damage <- 0; consec_gap <- 0L; pos <- si
    while (pos < n) {
      pos <- pos + 1L
      if (claimed[pos]) {
        # Already claimed by THIS family — skip
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        if (damage > DMG_MAX) break
        next
      }

      # Inv-likeness softener: if this window looks inversion-like from PCA,
      # reduce the effective damage penalty. The snake tolerates low-z windows
      # that are structurally inversion-like. Max reduction = 50% of damage.
      inv_soft <- 0
      if (has_inv && is.finite(dt$inv_likeness[pos])) {
        inv_soft <- dt$inv_likeness[pos] * 0.5  # 0 to 0.5 reduction
      }
      eff_damage <- damage * (1 - inv_soft)

      eff_accept  <- params$accept + eff_damage
      eff_tolerate <- params$tolerate + eff_damage

      roll <- tail(snake, params$roll_size)
      prev <- snake[length(snake)]
      sc <- compute_continuity(pos, prev, roll, sim_mat, mds_mat)

      if (sc >= eff_accept) {
        snake <- c(snake, pos); status <- c(status, "accepted"); scores <- c(scores, sc)
        claimed[pos] <- TRUE
        damage <- max(0, damage + DMG_RECOVER_ACCEPTED)
        consec_gap <- 0L
        log_decision(chr, sid, "core", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "accepted", paste0("dmg=", round(damage, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else if (sc >= eff_tolerate) {
        snake <- c(snake, pos); status <- c(status, "tolerated"); scores <- c(scores, sc)
        claimed[pos] <- TRUE
        damage <- damage + DMG_COST_TOLERATED
        consec_gap <- 0L
        log_decision(chr, sid, "core", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     "tolerated", paste0("dmg=", round(damage, 3)),
                     consec_gap, 0L, length(roll), params$name)
      } else {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        log_decision(chr, sid, "core", "right", prev, pos,
                     dt$global_window_id[prev], dt$global_window_id[pos],
                     sc, eff_accept, eff_tolerate,
                     if (damage > DMG_MAX) "stopped" else "rejected",
                     paste0("dmg=", round(damage, 3), "_consec=", consec_gap),
                     consec_gap, 0L, length(roll), params$name)
        if (damage > DMG_MAX) break
      }
    }

    # ── EXTEND LEFT with cumulative damage ───────────────────────────
    damage <- 0; consec_gap <- 0L; pos <- si
    while (pos > 1L) {
      pos <- pos - 1L
      if (claimed[pos]) {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        if (damage > DMG_MAX) break
        next
      }

      # Inv-likeness softener (same as right extension)
      inv_soft_l <- 0
      if (has_inv && is.finite(dt$inv_likeness[pos])) {
        inv_soft_l <- dt$inv_likeness[pos] * 0.5
      }
      eff_damage_l <- damage * (1 - inv_soft_l)

      eff_accept  <- params$accept + eff_damage_l
      eff_tolerate <- params$tolerate + eff_damage_l

      roll <- head(snake, params$roll_size)
      prev <- snake[1]
      sc <- compute_continuity(pos, prev, roll, sim_mat, mds_mat)

      if (sc >= eff_accept) {
        snake <- c(pos, snake); status <- c("accepted", status); scores <- c(sc, scores)
        claimed[pos] <- TRUE
        damage <- max(0, damage + DMG_RECOVER_ACCEPTED)
        consec_gap <- 0L
      } else if (sc >= eff_tolerate) {
        snake <- c(pos, snake); status <- c("tolerated", status); scores <- c(sc, scores)
        claimed[pos] <- TRUE
        damage <- damage + DMG_COST_TOLERATED
        consec_gap <- 0L
      } else {
        consec_gap <- consec_gap + 1L
        damage <- damage + DMG_COST_GAP + (consec_gap - 1L) * DMG_COST_CONSEC_GAP
        if (damage > DMG_MAX) break
      }
    }

    if (length(snake) >= params$min_windows) {
      regions[[length(regions) + 1]] <- list(
        idxs = snake, statuses = status, scores = scores,
        snake_id = sid, core_family = params$name
      )
    } else {
      claimed[snake] <- FALSE  # release windows from failed snake
    }
  }

  list(regions = regions, next_id = sid)
}

# =============================================================================
# SAMPLE-IDENTITY COMPATIBILITY CHECK
# =============================================================================
# The merge must NOT stitch cores that capture DIFFERENT sample groups.
# Two cores can be close in distance and similar in PCA geometry but driven
# by completely different samples (e.g., different families).
#
# Strategy: for each core, identify which samples are "outlier" (high |PC1|)
# using the per-window sample PC loadings from STEP09. Compare the outlier
# sample sets between cores using Jaccard overlap. If overlap is too low,
# these are different biological systems and should NOT be merged.

# Extract sample PC1 loadings matrix for a set of window indices
# Returns: samples × windows matrix of PC1 loadings
get_core_sample_loadings <- function(dt, idxs, sample_names) {
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) == 0) return(NULL)

  mat <- as.matrix(dt[idxs, ..available])
  # Average across windows in the core → one loading per sample
  if (nrow(mat) > 1) {
    avg_loadings <- colMeans(mat, na.rm = TRUE)
  } else {
    avg_loadings <- mat[1, ]
  }
  avg_loadings
}

# Classify samples into outlier groups based on PC1 loadings
classify_core_samples <- function(loadings) {
  if (is.null(loadings) || all(is.na(loadings))) return(NULL)
  loadings[is.na(loadings)] <- 0

  # Use median ± 1 MAD as split point
  med <- median(loadings)
  m <- mad(loadings)
  if (!is.finite(m) || m < 1e-8) m <- sd(loadings)
  if (!is.finite(m) || m < 1e-8) return(NULL)

  high_idx <- which(loadings > med + 0.8 * m)
  low_idx  <- which(loadings < med - 0.8 * m)

  list(high = high_idx, low = low_idx,
       outlier = sort(union(high_idx, low_idx)))
}

# Compute Jaccard overlap between two outlier sample sets
sample_group_jaccard <- function(set_a, set_b) {
  if (length(set_a) == 0 || length(set_b) == 0) return(0)
  inter <- length(intersect(set_a, set_b))
  union_size <- length(union(set_a, set_b))
  if (union_size == 0) return(0)
  inter / union_size
}

# Check if two cores are driven by compatible sample groups
# Returns: list(compatible, jaccard, reason)
check_core_compatibility <- function(dt, core_a_idxs, core_b_idxs,
                                      sample_names, min_jaccard = 0.30) {
  load_a <- get_core_sample_loadings(dt, core_a_idxs, sample_names)
  load_b <- get_core_sample_loadings(dt, core_b_idxs, sample_names)

  if (is.null(load_a) || is.null(load_b)) {
    return(list(compatible = TRUE, jaccard = NA_real_,
                reason = "no_loadings_available"))
  }

  grp_a <- classify_core_samples(load_a)
  grp_b <- classify_core_samples(load_b)

  if (is.null(grp_a) || is.null(grp_b)) {
    return(list(compatible = TRUE, jaccard = NA_real_,
                reason = "insufficient_variance"))
  }

  # Compare outlier sets (samples that are extreme in PC1)
  jacc <- sample_group_jaccard(grp_a$outlier, grp_b$outlier)

  # Also check directional consistency: do the HIGH samples overlap?
  jacc_high <- sample_group_jaccard(grp_a$high, grp_b$high)
  jacc_low  <- sample_group_jaccard(grp_a$low, grp_b$low)

  # Two cores are compatible if their outlier sets substantially overlap
  # OR if one is a subset of the other (shoulder of same system)
  subset_frac <- if (length(grp_a$outlier) > 0 && length(grp_b$outlier) > 0) {
    max(length(intersect(grp_a$outlier, grp_b$outlier)) / length(grp_a$outlier),
        length(intersect(grp_a$outlier, grp_b$outlier)) / length(grp_b$outlier))
  } else 0

  compatible <- jacc >= min_jaccard || subset_frac >= 0.50

  reason <- if (compatible) "sample_groups_compatible"
            else paste0("different_systems_jacc=", round(jacc, 3))

  list(compatible = compatible, jaccard = jacc,
       jaccard_high = jacc_high, jaccard_low = jacc_low,
       subset_frac = subset_frac, reason = reason)
}

# =============================================================================
# PARAMETERIZED MERGE ENGINE (with sample-identity gate)
# =============================================================================

run_merge_family <- function(all_cores, dt, sim_mat, mds_mat, chr, params,
                             sample_names = NULL) {
  if (length(all_cores) == 0) return(list())

  core_starts <- vapply(all_cores, function(c) min(c$idxs), integer(1))
  all_cores <- all_cores[order(core_starts)]

  # ── BRIDGE DAMAGE MODEL ─────────────────────────────────────────────
  # Same cumulative damage concept as core extension.
  # Bridge windows accumulate damage from gaps and tolerated windows.
  # Good bridge windows heal. Damage raises effective thresholds.
  # The gap between cores is treated as initial damage that the bridge
  # must overcome.
  DMG_BRIDGE_RECOVER    <- params$dmg_recover    %||% -0.04
  DMG_BRIDGE_TOLERATED  <- params$dmg_tolerate   %||%  0.03
  DMG_BRIDGE_GAP        <- params$dmg_gap         %||%  0.06
  DMG_BRIDGE_MAX        <- params$dmg_max          %||%  0.35
  max_gap_hard <- params$max_gap_between + 3L

  merged <- list()
  i <- 1L

  while (i <= length(all_cores)) {
    current <- all_cores[[i]]
    cur_idxs <- current$idxs
    cur_stat <- current$statuses
    cur_scor <- current$scores

    while (i < length(all_cores)) {
      nxt <- all_cores[[i + 1]]
      nxt_start <- min(nxt$idxs)
      cur_end   <- max(cur_idxs)
      gap_size  <- nxt_start - cur_end - 1L

      # Hard ceiling: still reject impossibly large gaps or overlaps
      # Hard ceiling still needed — can't merge infinitely far
      if (gap_size > max_gap_hard || gap_size < 0) break

      # Initial damage from the gap itself (before bridge even starts)
      bridge_damage <- gap_size * DMG_BRIDGE_GAP

      # If gap alone already exceeds damage max, don't attempt bridge
      if (bridge_damage > DMG_BRIDGE_MAX) break

      # Staging vector for bridge (C4 fix)
      bridge_range <- if (gap_size > 0) seq(cur_end + 1L, nxt_start - 1L) else integer(0)
      bridge_range <- bridge_range[bridge_range >= 1 & bridge_range <= nrow(dt)]

      bridge_ok <- TRUE
      staged_idxs <- integer(0)
      staged_stat <- character(0)
      staged_scor <- numeric(0)

      if (length(bridge_range) > 0) {
        eval_ctx <- cur_idxs
        for (bi in bridge_range) {
          prev <- eval_ctx[length(eval_ctx)]
          roll <- tail(eval_ctx, params$roll_size)
          sc <- compute_continuity(bi, prev, roll, sim_mat, mds_mat)

          # Effective thresholds raised by accumulated bridge damage
          eff_accept  <- params$accept + bridge_damage
          eff_tolerate <- params$tolerate + bridge_damage

          if (sc >= eff_accept) {
            staged_idxs <- c(staged_idxs, bi)
            staged_stat <- c(staged_stat, "bridge_accepted")
            staged_scor <- c(staged_scor, sc)
            eval_ctx <- c(eval_ctx, bi)
            bridge_damage <- max(0, bridge_damage + DMG_BRIDGE_RECOVER)
          } else if (sc >= eff_tolerate) {
            staged_idxs <- c(staged_idxs, bi)
            staged_stat <- c(staged_stat, "bridge_tolerated")
            staged_scor <- c(staged_scor, sc)
            eval_ctx <- c(eval_ctx, bi)
            bridge_damage <- bridge_damage + DMG_BRIDGE_TOLERATED
          } else {
            bridge_damage <- bridge_damage + DMG_BRIDGE_GAP
            if (bridge_damage > DMG_BRIDGE_MAX) { bridge_ok <- FALSE; break }
          }

          if (bridge_damage > DMG_BRIDGE_MAX) { bridge_ok <- FALSE; break }
        }
      }

      if (!bridge_ok) break

      # ── SAMPLE-IDENTITY COMPATIBILITY GATE ─────────────────────────
      if (!is.null(sample_names) && length(sample_names) > 0) {
        compat <- check_core_compatibility(
          dt, cur_idxs, nxt$idxs, sample_names,
          min_jaccard = params$min_sample_jaccard %||% 0.30
        )
        if (!compat$compatible) {
          log_decision(chr, NA_integer_, "merge", "sample_gate",
                       max(cur_idxs), min(nxt$idxs),
                       dt$global_window_id[max(cur_idxs)],
                       dt$global_window_id[min(nxt$idxs)],
                       compat$jaccard %||% 0, 0.30, 0, "blocked",
                       compat$reason, 0L, 0L, 0L, params$name)
          break
        }
      }

      proposed <- c(cur_idxs, staged_idxs, nxt$idxs)
      coh <- region_coherence(proposed, sim_mat)

      # Coherence floor raised by residual bridge damage
      effective_coherence_min <- params$coherence_min + bridge_damage * 0.5
      if (coh < effective_coherence_min) {
        log_decision(chr, NA_integer_, "merge", "coherence_damage_penalized",
                     max(cur_idxs), min(nxt$idxs),
                     dt$global_window_id[max(cur_idxs)],
                     dt$global_window_id[min(nxt$idxs)],
                     coh, effective_coherence_min, params$coherence_min,
                     "rejected",
                     paste0("coherence_", round(coh, 3), "_below_dmg_floor_",
                            round(effective_coherence_min, 3),
                            "_bridge_dmg=", round(bridge_damage, 3)),
                     gap_size, 0L, 0L, params$name)
        break
      }

      cur_idxs <- proposed
      cur_stat <- c(cur_stat, staged_stat, nxt$statuses)
      cur_scor <- c(cur_scor, staged_scor, nxt$scores)
      i <- i + 1L
    }

    if (length(cur_idxs) >= params$min_windows) {
      # Compute merge confidence: penalized by total gap accumulated
      total_bridge <- sum(grepl("^bridge", cur_stat))
      total_windows <- length(cur_idxs)
      bridge_frac <- total_bridge / max(1, total_windows)
      merge_confidence <- max(0, 1 - bridge_frac * 0.5) *
                          region_coherence(cur_idxs, sim_mat)

      merged[[length(merged) + 1]] <- list(
        idxs = cur_idxs, statuses = cur_stat, scores = cur_scor,
        coherence = region_coherence(cur_idxs, sim_mat),
        merge_confidence = round(merge_confidence, 4),
        merge_family = params$name
      )
    }
    i <- i + 1L
  }
  merged
}

# =============================================================================
# QC ENGINE (shared, parameterized)
# =============================================================================

run_qc <- function(merged, dt, sim_mat) {
  if (length(merged) == 0) return(data.table())
  diag_rows <- list()
  for (mi in seq_along(merged)) {
    reg <- merged[[mi]]; idxs <- reg$idxs; n <- length(idxs); scores <- reg$scores
    flags <- character(0)
    n_edge <- min(2, floor(n / 2))
    if (any(scores[seq_len(n_edge)] < QC_BOUNDARY_WEAK, na.rm = TRUE)) flags <- c(flags, "weak_left")
    if (any(scores[(n - n_edge + 1):n] < QC_BOUNDARY_WEAK, na.rm = TRUE)) flags <- c(flags, "weak_right")

    isims <- numeric(n); isims[1] <- 1.0
    for (k in 2:n) { s <- sim_mat[idxs[k-1], idxs[k]]; isims[k] <- if (is.finite(s)) s else 0 }
    n_drops <- sum(isims < QC_INTERNAL_DROP)
    if (n_drops / n > QC_DROP_FRAC) flags <- c(flags, "internal_disc")
    if (n >= 6) {
      wk <- which.min(isims[3:(n-2)]) + 2L
      if (isims[wk] < QC_SPLIT_THRESH) flags <- c(flags, "split_suggested")
    }
    n_bt <- sum(reg$statuses == "bridge_tolerated", na.rm = TRUE)
    n_core <- sum(reg$statuses %in% c("seed", "accepted"), na.rm = TRUE)
    if (n_core > 0 && n_bt / n_core > 0.5) flags <- c(flags, "overmerge")
    if (!is.null(reg$coherence) && reg$coherence < 0.25) flags <- c(flags, "low_coherence")
    if (length(flags) == 0) flags <- "clean"

    diag_rows[[mi]] <- data.table(
      merged_idx = mi, n_windows = n,
      coherence = round(reg$coherence %||% region_coherence(idxs, sim_mat), 4),
      mean_isim = round(mean(isims), 4), min_isim = round(min(isims), 4),
      qc_flags = paste(flags, collapse = ";"),
      merge_family = reg$merge_family %||% "unknown"
    )
  }
  rbindlist(diag_rows)
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_window_rows  <- list()
all_region_rows  <- list()
all_hier_rows    <- list()
all_diag_rows    <- list()
all_summary_rows <- list()
all_state_rows   <- list()
all_compare_rows <- list()
snake_id <- 0L

for (chr in chroms) {
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) next

  message("\n[v3] ═══════ ", chr, " ═══════")
  prep <- tryCatch(prepare_chr(chr_obj), error = function(e) {
    message("[WARN] ", chr, ": ", e$message); NULL
  })
  if (is.null(prep) || nrow(prep$dt) < 3) { message("[SKIP] ", chr); next }

  dt <- prep$dt; sim_mat <- prep$sim_mat; mds_mat <- prep$mds_mat
  n <- nrow(dt)

  # ═══════════════════════════════════════════════════════════════════
  # PHASE 1: Three core families (S → M → L, INDEPENDENT membership)
  # ═══════════════════════════════════════════════════════════════════
  # Each family builds its own view. A window claimed by 1S can also be
  # claimed by 1M and 1L — this gives true nested S/M/L structure.
  # Deduplication happens downstream (merge families see all cores,
  # diagnostic tables record which families claim which windows).

  all_cores <- list()
  for (fam in CORE_FAMILIES) {
    # Pass a blank used_global — each family explores independently
    result <- run_core_family(dt, sim_mat, mds_mat, chr, snake_id, rep(FALSE, n), fam)
    snake_id <- result$next_id
    message("[v3] ", fam$name, " (", fam$label, "): ", length(result$regions), " cores")

    for (reg in result$regions) {
      all_cores[[length(all_cores) + 1]] <- reg
      idxs <- reg$idxs
      for (k in seq_along(idxs)) {
        all_window_rows[[length(all_window_rows) + 1]] <- data.table(
          chrom = chr, global_window_id = dt$global_window_id[idxs[k]],
          start_bp = dt$start_bp[idxs[k]], end_bp = dt$end_bp[idxs[k]],
          snake_id = reg$snake_id, snake_phase = "core",
          core_family = fam$name, merge_family = NA_character_,
          inclusion_status = reg$statuses[k],
          continuity_score = round(reg$scores[k], 4)
        )
      }
      all_region_rows[[length(all_region_rows) + 1]] <- data.table(
        snake_id = reg$snake_id, snake_phase = "core",
        core_family = fam$name, merge_family = NA_character_, chrom = chr,
        start_bp = min(dt$start_bp[idxs]), end_bp = max(dt$end_bp[idxs]),
        n_windows = length(idxs),
        n_seeds = sum(reg$statuses == "seed"),
        n_tolerated = sum(reg$statuses == "tolerated"),
        mean_score = round(mean(reg$scores), 4),
        min_score = round(min(reg$scores), 4),
        coherence = round(region_coherence(idxs, sim_mat), 4)
      )
    }
  }

  message("[v3] Total cores: ", length(all_cores))

  # ═══════════════════════════════════════════════════════════════════
  # PHASE 2: Two merge families (A and B, independently on same cores)
  # ═══════════════════════════════════════════════════════════════════

  for (mfam in MERGE_FAMILIES) {
    merged <- run_merge_family(all_cores, dt, sim_mat, mds_mat, chr, mfam,
                               sample_names = sample_names_snake)
    message("[v3] ", mfam$name, " (", mfam$label, "): ", length(merged), " merged regions")

    merged_ids <- integer(length(merged))
    for (mi in seq_along(merged)) {
      snake_id <- snake_id + 1L
      merged_ids[mi] <- snake_id
      reg <- merged[[mi]]; idxs <- reg$idxs

      for (k in seq_along(idxs)) {
        all_window_rows[[length(all_window_rows) + 1]] <- data.table(
          chrom = chr, global_window_id = dt$global_window_id[idxs[k]],
          start_bp = dt$start_bp[idxs[k]], end_bp = dt$end_bp[idxs[k]],
          snake_id = snake_id, snake_phase = "merged",
          core_family = NA_character_, merge_family = mfam$name,
          inclusion_status = reg$statuses[k],
          continuity_score = round(reg$scores[k], 4)
        )
      }

      all_region_rows[[length(all_region_rows) + 1]] <- data.table(
        snake_id = snake_id, snake_phase = "merged",
        core_family = NA_character_, merge_family = mfam$name, chrom = chr,
        start_bp = min(dt$start_bp[idxs]), end_bp = max(dt$end_bp[idxs]),
        n_windows = length(idxs),
        n_seeds = sum(reg$statuses == "seed"),
        n_tolerated = sum(grepl("tolerate", reg$statuses)),
        mean_score = round(mean(reg$scores), 4),
        min_score = round(min(reg$scores), 4),
        coherence = round(reg$coherence %||% 0, 4)
      )
    }

    # Hierarchy
    for (mi in seq_along(merged)) {
      ms <- merged[[mi]]$idxs
      for (ci in seq_along(all_cores)) {
        ol <- length(intersect(all_cores[[ci]]$idxs, ms))
        if (ol > 0) {
          all_hier_rows[[length(all_hier_rows) + 1]] <- data.table(
            core_snake_id = all_cores[[ci]]$snake_id,
            core_family = all_cores[[ci]]$core_family,
            merged_snake_id = merged_ids[mi],
            merge_family = mfam$name, chrom = chr,
            overlap_windows = ol,
            overlap_frac = round(ol / length(all_cores[[ci]]$idxs), 4)
          )
        }
      }
    }

    # QC
    qc_dt <- run_qc(merged, dt, sim_mat)
    if (nrow(qc_dt) > 0) {
      qc_dt[, snake_id := merged_ids[merged_idx]]
      qc_dt[, chrom := chr]
      all_diag_rows[[length(all_diag_rows) + 1]] <- qc_dt
    }
  }

  # ═══════════════════════════════════════════════════════════════════
  # PHASE 3: Merge A vs B comparison
  # ═══════════════════════════════════════════════════════════════════

  # Find where merge_A and merge_B agree or disagree
  reg_dt <- rbindlist(all_region_rows[vapply(all_region_rows,
    function(x) x$chrom[1] == chr && x$snake_phase == "merged", logical(1))], fill = TRUE)

  if (nrow(reg_dt) > 0) {
    rA <- reg_dt[merge_family == "merge_A"]
    rB <- reg_dt[merge_family == "merge_B"]

    # For each merge_B region, check if merge_A also found it (overlapping)
    for (bi in seq_len(nrow(rB))) {
      b_start <- rB$start_bp[bi]; b_end <- rB$end_bp[bi]
      a_overlap <- rA[start_bp <= b_end & end_bp >= b_start]

      all_compare_rows[[length(all_compare_rows) + 1]] <- data.table(
        chrom = chr,
        merge_B_id = rB$snake_id[bi],
        merge_B_start = b_start, merge_B_end = b_end,
        merge_B_windows = rB$n_windows[bi],
        merge_B_coherence = rB$coherence[bi],
        n_merge_A_overlapping = nrow(a_overlap),
        merge_A_ids = if (nrow(a_overlap) > 0) paste(a_overlap$snake_id, collapse = ";") else NA,
        merge_A_total_windows = if (nrow(a_overlap) > 0) sum(a_overlap$n_windows) else 0L,
        agreement = if (nrow(a_overlap) > 0) "both_merge" else "B_only_merge_sensitive"
      )
    }
  }

  # Window states
  collected_wids <- unlist(lapply(all_cores, function(r) dt$global_window_id[r$idxs]))
  for (wi in seq_len(n)) {
    wid <- dt$global_window_id[wi]
    
    # Determine seed status using the broadest family threshold (1L) + inv_likeness
    has_inv_col <- "inv_likeness" %in% names(dt)
    inv_l <- if (has_inv_col && is.finite(dt$inv_likeness[wi])) dt$inv_likeness[wi] else 0
    is_seed_like <- (dt$max_abs_z[wi] >= 1.2 | inv_l >= 0.90) &
                    is.finite(dt$seed_nn_dist[wi]) & dt$seed_nn_dist[wi] < 0.90

    all_state_rows[[length(all_state_rows) + 1]] <- data.table(
      chrom = chr, global_window_id = wid,
      start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
      max_abs_z = round(dt$max_abs_z[wi], 4),
      inv_likeness = round(inv_l, 4),
      is_seed_like = is_seed_like,
      collector_state = if (is_seed_like & !(wid %in% collected_wids)) "seed_uncollected"
                        else if (wid %in% collected_wids) "collected"
                        else "unused"
    )
  }

  # Summary
  n_s <- sum(vapply(all_cores, function(r) r$core_family == "1S", logical(1)))
  n_m <- sum(vapply(all_cores, function(r) r$core_family == "1M", logical(1)))
  n_l <- sum(vapply(all_cores, function(r) r$core_family == "1L", logical(1)))
  all_summary_rows[[length(all_summary_rows) + 1]] <- data.table(
    chrom = chr, n_windows = n,
    cores_1S = n_s, cores_1M = n_m, cores_1L = n_l,
    cores_total = length(all_cores)
  )
}

# =============================================================================
# WRITE
# =============================================================================

empty <- function(...) data.table()
win_dt   <- if (length(all_window_rows) > 0) rbindlist(all_window_rows, fill = TRUE) else empty()
reg_dt   <- if (length(all_region_rows) > 0) rbindlist(all_region_rows, fill = TRUE) else empty()
hier_dt  <- if (length(all_hier_rows) > 0) rbindlist(all_hier_rows, fill = TRUE) else empty()
diag_dt  <- if (length(all_diag_rows) > 0) rbindlist(all_diag_rows, fill = TRUE) else empty()
summ_dt  <- if (length(all_summary_rows) > 0) rbindlist(all_summary_rows) else empty()
log_dt   <- if (length(decision_log) > 0) rbindlist(decision_log, fill = TRUE) else empty()
state_dt <- if (length(all_state_rows) > 0) rbindlist(all_state_rows) else empty()
comp_dt  <- if (length(all_compare_rows) > 0) rbindlist(all_compare_rows, fill = TRUE) else empty()

# =============================================================================
# DEBUG: Per-region MDS facet wraps
#
# For each merge_A region with ≥20 windows, produce a facet-wrapped MDS1×2
# scatter plot showing ALL windows in that region. Each facet = one window,
# colored by the 3-band assignment (tailA/middle/tailB). Strip background
# matches the region color from the ideogram plot.
#
# Output: <outdir>/debug_mds_facets/<chr>_mergeA_<id>.pdf
# =============================================================================

debug_dir <- file.path(outdir, "debug_mds_facets")
if (nrow(reg_dt) > 0 && !is.null(mds_obj)) {
  suppressPackageStartupMessages({ library(ggplot2) })
  dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)

  # Use merge_A regions (primary method) — skip very small ones
  merge_a_regs <- reg_dt[grepl("merge_A", family)]
  if (nrow(merge_a_regs) == 0) merge_a_regs <- reg_dt  # fallback to all

  # Palette for region IDs (cycling through distinguishable colors)
  reg_pal <- c("#dc2626", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
               "#059669", "#e11d48", "#0284c7", "#65a30d", "#c026d3")

  n_plotted <- 0L
  for (ri in seq_len(nrow(merge_a_regs))) {
    r <- merge_a_regs[ri]
    chr <- r$chrom
    chr_obj <- mds_obj$per_chr[[chr]]
    if (is.null(chr_obj)) next

    dt_chr <- as.data.table(chr_obj$out_dt)[order(start_bp)]

    # Get windows in this region
    r_wins <- win_dt[chrom == chr &
                     start_bp >= r$start_bp & end_bp <= r$end_bp]
    if (nrow(r_wins) < 5) next  # skip tiny regions

    # Get the MDS coordinates + PC1 scores for these windows
    wids <- unique(r_wins$global_window_id)
    dt_sub <- dt_chr[global_window_id %in% wids]
    if (nrow(dt_sub) < 5) next

    pc1_cols <- grep("^PC_1_", names(dt_sub), value = TRUE)
    if (length(pc1_cols) == 0) next

    # Build facet data: for each window, get sample PC1 scores + 3-band assignment
    facet_rows <- list()
    for (wi in seq_len(nrow(dt_sub))) {
      scores <- as.numeric(dt_sub[wi, ..pc1_cols])
      if (all(is.na(scores))) next

      # Quick k=3 for coloring
      km <- tryCatch(kmeans(scores[is.finite(scores)], centers = 3, nstart = 3),
                     error = function(e) NULL)
      if (is.null(km)) next
      centers <- sort(km$centers[, 1])
      center_order <- order(km$centers[, 1])
      role_map <- setNames(c("tailA", "middle", "tailB"), as.character(center_order))
      roles <- role_map[as.character(km$cluster)]

      # PC2 if available
      pc2_cols <- grep("^PC_2_", names(dt_sub), value = TRUE)
      pc2_scores <- if (length(pc2_cols) > 0) as.numeric(dt_sub[wi, ..pc2_cols]) else rep(0, length(scores))

      wid <- dt_sub$global_window_id[wi]
      pos_mb <- round((dt_sub$start_bp[wi] + dt_sub$end_bp[wi]) / 2e6, 2)

      valid <- is.finite(scores)
      facet_rows[[length(facet_rows) + 1]] <- data.table(
        window_label = paste0("w", wid, " (", pos_mb, " Mb)"),
        PC1 = scores[valid], PC2 = pc2_scores[valid],
        band = roles[valid],
        window_order = wi
      )
    }

    if (length(facet_rows) == 0) next
    fdt <- rbindlist(facet_rows)

    # Sort facets by genomic order
    lev <- unique(fdt$window_label[order(fdt$window_order)])
    fdt[, window_label := factor(window_label, levels = lev)]

    reg_color <- reg_pal[((ri - 1L) %% length(reg_pal)) + 1L]
    reg_label <- paste0(chr, " ", r$family, " #", r$snake_id,
                        " (", round(r$start_bp/1e6, 1), "–", round(r$end_bp/1e6, 1), " Mb, ",
                        nrow(dt_sub), " windows)")

    # Determine layout
    n_facets <- length(lev)
    ncol_facet <- min(6L, n_facets)
    nrow_facet <- ceiling(n_facets / ncol_facet)
    plot_h <- max(4, nrow_facet * 2.5)
    plot_w <- max(6, ncol_facet * 2.5)

    g <- ggplot(fdt, aes(x = PC1, y = PC2, color = band)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(values = c(tailA = "#2563eb", middle = "#f59e0b", tailB = "#dc2626")) +
      facet_wrap(~ window_label, ncol = ncol_facet, scales = "free") +
      labs(title = reg_label, color = "Band") +
      theme_minimal(base_size = 7) +
      theme(
        strip.background = element_rect(fill = reg_color, color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 5),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold")
      )

    fname <- paste0(chr, "_", r$family, "_", r$snake_id, ".pdf")
    ggsave(file.path(debug_dir, fname), g, width = min(plot_w, 20), height = min(plot_h, 30))
    n_plotted <- n_plotted + 1L
  }

  if (n_plotted > 0) message("[v3] Debug MDS facets: ", n_plotted, " PDFs → ", debug_dir)
}

files <- c(
  snake_windows = "snake_windows.tsv.gz",
  snake_regions = "snake_regions.tsv.gz",
  snake_hierarchy = "snake_hierarchy.tsv.gz",
  snake_summary = "snake_summary.tsv",
  snake_diagnostics = "snake_diagnostics.tsv.gz",
  snake_decision_log = "snake_decision_log.tsv.gz",
  snake_window_states = "snake_window_states.tsv.gz",
  snake_multiscale_comparison = "snake_multiscale_comparison.tsv.gz"
)

dts <- list(win_dt, reg_dt, hier_dt, summ_dt, diag_dt, log_dt, state_dt, comp_dt)
for (i in seq_along(files)) {
  fwrite(dts[[i]], file.path(outdir, files[i]), sep = "\t")
}

message("\n[DONE] STEP10e_v3 multi-scale snake complete")
for (f in files) message("  ", file.path(outdir, f))
message("[v3] Cores: 1S=", sum(summ_dt$cores_1S), " 1M=", sum(summ_dt$cores_1M),
        " 1L=", sum(summ_dt$cores_1L))
message("[v3] Decision log: ", nrow(log_dt), " entries")
message("[v3] Merge comparisons: ", nrow(comp_dt), " (",
        sum(comp_dt$agreement == "B_only_merge_sensitive", na.rm = TRUE),
        " merge-sensitive)")

# =============================================================================
# EXPORT: snake_candidate_regions.tsv.gz
# =============================================================================
# This file bridges Snake 1 v3 output to MODULE_5B followup.
# It exports merge_A regions (conservative, high-confidence) as candidate
# regions in the same format as STEP10_v2's candidate_regions.tsv.gz.
#
# If merge_A finds nothing, falls back to core regions (ungrouped).
# merge_B (permissive) regions are NOT exported here — they are for
# exploratory comparison only.

snake_cand_rows <- list()
scid <- 0L

if (nrow(reg_dt) > 0) {
  # Primary: merge_A merged regions
  merge_a <- reg_dt[snake_phase == "merged" & merge_family == "merge_A"]

  if (nrow(merge_a) > 0) {
    for (ri in seq_len(nrow(merge_a))) {
      scid <- scid + 1L
      r <- merge_a[ri]
      # Get the window IDs for SNP counting
      cand_wins <- win_dt[snake_id == r$snake_id & snake_phase == "merged"]
      snake_cand_rows[[scid]] <- data.table(
        candidate_id   = scid,
        source         = "snake_merge_A",
        snake_id       = r$snake_id,
        chrom          = r$chrom,
        start_bp       = r$start_bp,
        end_bp         = r$end_bp,
        n_windows      = r$n_windows,
        mean_score     = r$mean_score,
        coherence      = r$coherence,
        merge_family   = "merge_A"
      )
    }
  }

  # Fallback: ungrouped cores (if no merge_A regions exist for a chromosome)
  merge_a_chroms <- unique(merge_a$chrom)
  cores_no_merge <- reg_dt[snake_phase == "core" & !(chrom %in% merge_a_chroms)]

  if (nrow(cores_no_merge) > 0) {
    for (ri in seq_len(nrow(cores_no_merge))) {
      scid <- scid + 1L
      r <- cores_no_merge[ri]
      snake_cand_rows[[scid]] <- data.table(
        candidate_id   = scid,
        source         = paste0("snake_core_", r$core_family),
        snake_id       = r$snake_id,
        chrom          = r$chrom,
        start_bp       = r$start_bp,
        end_bp         = r$end_bp,
        n_windows      = r$n_windows,
        mean_score     = r$mean_score,
        coherence      = r$coherence,
        merge_family   = NA_character_
      )
    }
  }
}

snake_cand_dt <- if (length(snake_cand_rows) > 0) {
  rbindlist(snake_cand_rows, fill = TRUE)
} else {
  data.table(candidate_id = integer(), source = character(),
             snake_id = integer(), chrom = character(),
             start_bp = numeric(), end_bp = numeric(),
             n_windows = integer(), mean_score = numeric(),
             coherence = numeric(), merge_family = character())
}

scand_file <- file.path(outdir, "snake_candidate_regions.tsv.gz")
fwrite(snake_cand_dt, scand_file, sep = "\t")
message("[v3] Snake candidate regions: ", nrow(snake_cand_dt), " → ", scand_file)

# Inversion-likeness per-window table (for diagnostic plots + background filtering)
if (nrow(inv_like_dt) > 0) {
  f_inv <- file.path(outdir, "snake_inv_likeness.tsv.gz")
  fwrite(inv_like_dt, f_inv, sep = "\t")
  message("[v3] Inversion-likeness scores: ", nrow(inv_like_dt), " windows → ", f_inv)
}
