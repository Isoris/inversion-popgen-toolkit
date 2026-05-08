#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01a_precompute.R
#
# Phase 2 / 2c — precomputation for seeded region-growing (MDS z-outlier
# seed → extension under similarity-based damage cost). Run ONCE, reuse
# many times by downstream runs of C01b with different parameters.
#
# Codebase:    inversion_modules v8.5 / script v9.3.4
# Upstream:    phase_2/2b MDS output — ${MDS_PREFIX}.mds.rds
#              phase_2/2c STEP_C00 output — ${SV_PRIOR_DIR}/sv_prior_<chr>.rds
#              (optional) dosage files in ${DOSAGE_DIR} (second-pass het upgrade)
# Downstream:  phase_2/2d C01b (seeded region-growing runner)
# Formerly:    STEP_C01a_snake1_precompute_wired7_12358_v934.R
#
# Purpose
# -------
# Loads the large .mds.rds (~12 GB) and computes per-chromosome annotations
# that feed the seeded-extension step. The precompute is deliberately slow
# (once, ~1 hour) so that C01b runs in 2-5 minutes when tuning parameters.
#
# Per-window annotations produced
# -------------------------------
#   - Robust z-scores (median / MAD) on each MDS axis, per-chromosome
#   - Similarity matrices derived from lostruct distance matrices
#   - Seed nearest-neighbour distances (drives seed eligibility)
#   - Inversion-likeness score:
#       inv_likeness = 0.45*het + 0.30*trimodality + 0.25*band_discreteness
#       (PVE1 / structure_likeness removed in v9.3 — was circular)
#   - Family-likeness (diagnostic, NOT a gate)
#   - Band Q-stamps (test_05 precompute): per-window per-band dominant
#     Q component, fraction, cross-band Q entropy
#   - Dosage het rate (independent from PCA — second-pass upgrade when
#     --dosage_dir is provided; dosage CV replaces het_contrast as the
#     het component in inv_likeness for windows that have it)
#
# Seeding model
# -------------
# MDS produces z-scored regime coordinates per window. A window becomes
# a seed candidate iff its |z| exceeds a z-threshold OR its inv_likeness
# is very high OR it passes the per-chromosome Beta-distribution adaptive
# threshold (test_07). These are candidates for C01b's seeded extension
# step; the extension itself (cost-accumulating growth along the chromosome)
# happens in C01b, not here.
#
# SV prior annotation (tests 01/02/03/08 from sv_prior)
# -----------------------------------------------------
# At the end of precompute, stamps SV columns onto inv_like_dt by reading
# the per-chromosome sv_prior RDS built by STEP_C00. These are pure
# annotation — they do not alter inv_likeness or any upstream computation.
# Downstream (merge, scoring, decomposition) consumes them.
#
# Additional diagnostics integrated
# ---------------------------------
#   test_05  Q-axis Fst scan via Engine B dispatcher (requires Q matrix
#            and load_bridge.R; falls back to ANOVA family_likeness if not
#            available)
#   test_07  Beta(alpha,beta) adaptive seed threshold — per-chromosome
#            fit on inv_likeness → p-value per window → adaptive_seed flag
#            as alternative to the fixed inv_likeness threshold
#
# Outputs (in <outdir>/)
# ----------------------
#   precomp/<chr>.precomp.rds       per-chr: dt, sim_mat, mds_mat, dmat
#   window_inv_likeness.tsv.gz      per-window inv-likeness scores
#   precomp_summary.tsv             per-chr window counts + timing
#
# Usage
# -----
#   Rscript STEP_C01a_precompute.R <step10_outprefix> <outdir> \
#     [--qmatrix <qopt>] [--relatedness <pairs>] \
#     [--dosage_dir <dir>] [--sv_prior_dir <dir>] \
#     [--mode hatchery|wild]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ── Source load_bridge.R (provides smap, reg, get_Q, get_region_stats) ───────
.bridge_file <- Sys.getenv("LOAD_BRIDGE", "")
if (!nzchar(.bridge_file)) {
  for (.bp in c("utils/load_bridge.R", "../utils/load_bridge.R",
                 file.path(Sys.getenv("BASE", ""), "inversion_codebase_v8.5/utils/load_bridge.R"))) {
    if (file.exists(.bp)) { .bridge_file <- .bp; break }
  }
}
.bridge_available <- FALSE
if (nzchar(.bridge_file) && file.exists(.bridge_file)) {
  tryCatch({
    source(.bridge_file)
    .bridge_available <- TRUE
    message("[PRECOMP] load_bridge.R sourced — Engine B available for test_05 Fst")
  }, error = function(e) {
    message("[PRECOMP] WARNING: load_bridge.R failed: ", conditionMessage(e))
  })
} else {
  message("[PRECOMP] load_bridge.R not found — test_05 Fst scan will be skipped (ANOVA family_likeness only)")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01a_precompute.R <step10_outprefix> <outdir> [--qmatrix <qopt>] [--relatedness <pairs>] [--mode hatchery|wild]")
}

step10_prefix <- args[1]
outdir        <- args[2]
precomp_dir   <- file.path(outdir, "precomp")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)

# Parse optional args
qmatrix_file <- NULL; relatedness_file <- NULL; pop_mode <- "hatchery"
dosage_dir <- NULL
sv_prior_dir <- NULL
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--qmatrix" && i < length(args))      { qmatrix_file <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args)) { relatedness_file <- args[i+1]; i <- i+2 }
  else if (a == "--dosage_dir" && i < length(args)) { dosage_dir <- args[i+1]; i <- i+2 }
  else if (a == "--sv_prior_dir" && i < length(args)) { sv_prior_dir <- args[i+1]; i <- i+2 }
  else if (a == "--mode" && i < length(args))    { pop_mode <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

# =============================================================================
# PARAMETERS (same as C01 — must match)
# =============================================================================

SEED_MDS_AXES    <- 5L
SEED_NEIGHBOR_K  <- 3L

# =============================================================================
# LOAD
# =============================================================================

mds_rds_file <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds_file)) stop("Missing: ", mds_rds_file)

message("[PRECOMP] Loading ", mds_rds_file, " ...")
t_load <- proc.time()
mds_obj <- readRDS(mds_rds_file)
message("[PRECOMP] Loaded in ", round((proc.time() - t_load)[3], 1), "s")

per_chr <- mds_obj$per_chr
if (is.null(per_chr) || length(per_chr) == 0) stop("No per_chr data")
chroms <- names(per_chr)
message("[PRECOMP] ", length(chroms), " chromosomes")

# Extract sample names
sample_names <- NULL
for (chr_tmp in chroms) {
  dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
  pc1_cols <- grep("^PC_1_", names(dt_tmp), value = TRUE)
  if (length(pc1_cols) > 0) {
    sample_names <- sub("^PC_1_", "", pc1_cols)
    break
  }
}
if (!is.null(sample_names)) {
  message("[PRECOMP] Sample names: ", length(sample_names))
}

# =============================================================================
# DETECT npc — how many PCs are present in the upstream MDS output
# =============================================================================
# STEP_A02 / STEP_A03 emit columns PC_1_<sample>, PC_2_<sample>, ...
# up to the requested --npc. We detect npc by scanning column names of any
# non-empty per-chr out_dt. This is the upper bound: a chromosome with too
# few SNPs to fit the full eigendecomposition may have fewer PCs in practice.
# Per-PC presence is checked at use-time. PC1 is mandatory; PC2/3/4 are
# additive — they propagate through to precomp.rds when present.
npc_detected <- 1L   # default fallback if nothing else is found
for (chr_tmp in chroms) {
  dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
  pc_prefix_re <- "^PC_([0-9]+)_"
  pc_cols_all <- grep(pc_prefix_re, names(dt_tmp), value = TRUE)
  if (length(pc_cols_all) > 0) {
    pc_indices <- as.integer(sub(pc_prefix_re, "\\1",
                                 regmatches(pc_cols_all,
                                            regexpr(pc_prefix_re, pc_cols_all))))
    npc_detected <- max(npc_detected, max(pc_indices, na.rm = TRUE))
    break
  }
}
message("[PRECOMP] Detected npc = ", npc_detected,
        " (PC1 + ", max(0L, npc_detected - 1L), " higher-order PCs propagated)")

# =============================================================================
# LOAD Q MATRIX (optional, for family-likeness track)
# =============================================================================
# Q matrix from NGSadmix: n_samples x K ancestry proportions.
# Used to compute per-window family_likeness: do k=3 PC1 bands correlate
# with ancestry composition? If yes → family-driven grouping.
# If no → structural signal (inversion/demographic).

q_mat <- NULL; K_ancestry <- 0L
if (!is.null(qmatrix_file) && file.exists(qmatrix_file)) {
  q_raw <- tryCatch(as.matrix(fread(qmatrix_file, header = FALSE)), error = function(e) NULL)
  if (!is.null(q_raw) && nrow(q_raw) >= length(sample_names)) {
    # NGSadmix .qopt has no header, rows = samples in BAM list order
    # Assume row order matches sample_names (from BAM list)
    if (nrow(q_raw) == length(sample_names)) {
      rownames(q_raw) <- sample_names
      q_mat <- q_raw
      K_ancestry <- ncol(q_mat)
      message("[PRECOMP] Q matrix loaded: ", nrow(q_mat), " samples x K=", K_ancestry)
    } else {
      message("[PRECOMP] Q matrix row count (", nrow(q_raw),
              ") != sample count (", length(sample_names), ") — skipping")
    }
  } else {
    message("[PRECOMP] Q matrix file not readable: ", qmatrix_file)
  }
} else if (!is.null(qmatrix_file)) {
  message("[PRECOMP] Q matrix file not found: ", qmatrix_file, " — family track disabled")
}

# =============================================================================
# LOAD RELATEDNESS (optional, genome-wide pairwise kinship)
# =============================================================================
# From ngsRelate: pairwise kinship for all 226 × 226.
# NOT per-window — used as global supporting evidence in scoring,
# not as a per-window track. Stored for downstream C01d scoring.

kin_mat <- NULL
if (!is.null(relatedness_file) && file.exists(relatedness_file)) {
  kin_raw <- tryCatch(fread(relatedness_file), error = function(e) NULL)
  if (!is.null(kin_raw)) {
    # Expect columns: a, b, rab (or similar)
    kin_cols <- intersect(names(kin_raw), c("a", "b", "rab", "ida", "idb", "KING"))
    if (length(kin_cols) >= 3) {
      message("[PRECOMP] Relatedness loaded: ", nrow(kin_raw), " pairs")
      # Build symmetric matrix
      ids <- sort(unique(c(kin_raw[[1]], kin_raw[[2]])))
      n_ids <- length(ids)
      kin_mat <- matrix(0, n_ids, n_ids, dimnames = list(ids, ids))
      for (ri in seq_len(nrow(kin_raw))) {
        a <- as.character(kin_raw[[1]][ri])
        b <- as.character(kin_raw[[2]][ri])
        v <- kin_raw[[3]][ri]
        if (a %in% ids && b %in% ids) { kin_mat[a, b] <- v; kin_mat[b, a] <- v }
      }
      message("[PRECOMP] Kinship matrix: ", n_ids, " x ", n_ids)
    }
  }
}

message("[PRECOMP] Mode: ", pop_mode,
        " | Q matrix: ", if (!is.null(q_mat)) paste0("K=", K_ancestry) else "none",
        " | Relatedness: ", if (!is.null(kin_mat)) paste0(nrow(kin_mat), "x", ncol(kin_mat)) else "none")

# =============================================================================
# HELPER FUNCTIONS (defined before compute functions that use them)
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

make_sim_mat <- function(dmat) {
  finite_vals <- dmat[is.finite(dmat)]
  dmax <- if (length(finite_vals) > 0) quantile(finite_vals, 0.95, na.rm = TRUE) else 1
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

# =============================================================================
# DOSAGE-BASED HET RATE (genuinely independent from PCA)
# =============================================================================
# For each window, count sites per sample where dosage is het-like (0.3-1.7).
# This is NOT derived from PCA — it uses the raw allele dosage values.
# A window with many samples showing high het rate = real heterozygosity
# = consistent with inversion het genotype.
#
# Returns: data.table with global_window_id, dosage_het_rate_median, dosage_het_rate_sd,
#          dosage_het_rate_cv (coefficient of variation — high CV = bimodal het pattern)

compute_dosage_het_rates <- function(dosage_dir, per_chr, chroms) {
  if (is.null(dosage_dir) || !dir.exists(dosage_dir)) return(data.table())
  all_rows <- list()

  for (chr in chroms) {
    dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
    if (!file.exists(dos_file)) {
      dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv"))
      if (!file.exists(dos_file)) next
    }

    dt <- as.data.table(per_chr[[chr]]$out_dt)
    if (nrow(dt) == 0 || !"start_bp" %in% names(dt)) next

    # Load dosage matrix: rows = sites (with POS column), cols = samples
    dos <- tryCatch(fread(dos_file, nThread = 4L), error = function(e) NULL)
    if (is.null(dos) || nrow(dos) < 100) {
      message("  [dosage] ", chr, ": cannot load or too few sites")
      next
    }

    # Identify position column
    pos_col <- intersect(c("POS", "pos", "position"), names(dos))
    if (length(pos_col) == 0) { message("  [dosage] ", chr, ": no POS column"); next }
    positions <- dos[[pos_col[1]]]
    dos_mat <- as.matrix(dos[, !..pos_col[1]])

    message("  [dosage] ", chr, ": ", nrow(dos_mat), " sites x ", ncol(dos_mat), " samples")

    for (wi in seq_len(nrow(dt))) {
      w_start <- dt$start_bp[wi]; w_end <- dt$end_bp[wi]
      wid <- dt$global_window_id[wi]

      # Sites in this window
      site_idx <- which(positions >= w_start & positions <= w_end)
      if (length(site_idx) < 10) next

      w_dos <- dos_mat[site_idx, , drop = FALSE]

      # Per-sample het rate: fraction of sites with dosage in [0.3, 1.7]
      per_sample_het <- colMeans(w_dos >= 0.3 & w_dos <= 1.7, na.rm = TRUE)
      per_sample_het <- per_sample_het[is.finite(per_sample_het)]
      if (length(per_sample_het) < 20) next

      het_median <- median(per_sample_het)
      het_sd <- sd(per_sample_het)
      het_cv <- if (het_median > 0.01) het_sd / het_median else NA_real_

      # Bimodality of het rate distribution: are some samples very het
      # and others not? This is the inversion signal.
      # For an inversion: het samples have ~50% het sites, hom samples ~10%
      # → bimodal het rate distribution → high CV
      # For background: all samples similar het rate → low CV

      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        dosage_het_rate_median = round(het_median, 4),
        dosage_het_rate_sd = round(het_sd, 4),
        dosage_het_rate_cv = round(het_cv %||% NA_real_, 4)
      )
    }
    message("  [dosage] ", chr, ": ", length(all_rows), " windows processed")
  }
  if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
}

# =============================================================================
# INVERSION-LIKENESS (genome-wide, before per-chr precompute)
# =============================================================================

compute_inv_likeness_all <- function(per_chr, chroms) {
  all_rows <- list()
  for (chr in chroms) {
    dt <- as.data.table(per_chr[[chr]]$out_dt)
    if (nrow(dt) == 0) next
    has_lam <- "lam_1" %in% names(dt) && "lam_2" %in% names(dt)
    pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
    has_pc1 <- length(pc1_cols) > 0

    # NOTE (v9.3): PVE1 baselines removed. PVE1 is circular at this stage
    # (precomp windows are already MDS outliers, so PVE1 just says "outlier
    # is structured"). lam_1/lam_2 are still in out_dt and pass through to
    # the precomp RDS for anyone who wants to read them, but we no longer
    # derive inv_pve1, inv_eigen_ratio, or inv_pve1_excess.

    # ── CHROMOSOME-WIDE FIXED k=3 BANDS (v8.5) ──
    # Computed ONCE per chromosome. Used as fixed group labels for all
    # per-window entropy/delta12/het metrics. Avoids per-window k-means
    # (which is unstable, causes label switching, produces noisy tracks).
    #
    # Method: average PC1 across ALL windows per sample, then k=3 on
    # those averages. This gives the dominant 3-group structure of the
    # chromosome. Individual windows measure conformance to this fixed grouping.
    #
    # NOTE (v9.4): PC1 is the canonical band axis for inv_likeness. PC2/3/4
    # bands are computed separately in `precompute_one_chr` and saved to
    # the precomp RDS for downstream consumers (sub-cluster detection,
    # secondary-inversion candidate detection, scrubber multi-PC views).
    chr_bands <- NULL
    if (has_pc1 && length(pc1_cols) >= 20) {
      pc1_mat <- as.matrix(dt[, ..pc1_cols])
      avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
      valid <- is.finite(avg_pc1)
      if (sum(valid) >= 20) {
        avg_valid <- avg_pc1[valid]
        chr_km <- tryCatch(kmeans(avg_valid, centers = 3, nstart = 10),
                            error = function(e) NULL)
        if (!is.null(chr_km)) {
          co <- order(chr_km$centers[, 1])
          chr_bands <- integer(length(avg_valid))
          chr_bands[chr_km$cluster == co[1]] <- 1L  # low PC1 (band 1)
          chr_bands[chr_km$cluster == co[2]] <- 2L  # mid PC1 (band 2 = het-like)
          chr_bands[chr_km$cluster == co[3]] <- 3L  # high PC1 (band 3)
          names(chr_bands) <- sub("^PC_1_", "", names(avg_valid))
          message("  [", chr, "] Fixed bands: ",
                  sum(chr_bands == 1), "/", sum(chr_bands == 2), "/", sum(chr_bands == 3),
                  " (band1/band2/band3)")
        }
      }
    }

    for (wi in seq_len(nrow(dt))) {
      wid <- dt$global_window_id[wi]
      # PVE1 variables kept as NA for backward compatibility (v9.3: removed)
      pve1 <- NA_real_; eigen_ratio <- NA_real_; lambda_ratio <- NA_real_
      dip_p <- NA_real_; het_contrast <- NA_real_

      # ── Het and structure metrics (v8.5) ──
      # 
      # THREE SOURCES, clearly separated:
      #
      # SOURCE A: Per-window k-means k=3 on PC1 (het_contrast only)
      #   Used for: het_contrast (between/within variance)
      #   NOT independent from PCA
      #
      # SOURCE B: Chromosome-wide fixed k=3 bands applied to this window's PC1
      #   Used for: het_pc1_gap, het_mid_fraction, het_mid_variance,
      #             pc1_bimodality (Ashman D), local_delta12, local_entropy, local_ena
      #   NOT independent from PCA but STABLE (no label switching)
      #
      # SOURCE C: Raw dosage het counts (computed AFTER this loop)
      #   Used for: dosage_het_rate_median, dosage_het_rate_sd, dosage_het_rate_cv
      #   GENUINELY INDEPENDENT from PCA
      #
      # het_contrast (source A) = per-window k-means, unstable but captures local signal
      # All other het metrics (source B) = fixed bands, stable reference frame
      # dosage_het_rate (source C) = filled by post-loop merge, NA here

      # Source A: per-window k-means (only het_contrast)
      # (computed below in the existing km3 block)

      # Source B: fixed-band metrics
      het_pc1_gap <- NA_real_
      het_mid_fraction <- NA_real_
      het_mid_variance <- NA_real_
      pc1_bimodality <- NA_real_       # renamed from pc1_bimodality (honest name)
      local_delta12 <- NA_real_
      local_entropy <- NA_real_
      local_ena <- NA_real_
      band_sizes <- rep(NA_real_, 3)

      # Source C: dosage-based (filled by post-loop merge)
      dosage_het_rate_median <- NA_real_
      dosage_het_rate_sd <- NA_real_
      dosage_het_rate_cv <- NA_real_

      # NOTE (v9.3): lam_1/lam_2 eigenvalue computation removed.
      # PVE1 = lam_1/(lam_1+lam_2) is circular — these windows are
      # already MDS outliers, so PVE1 just says "this outlier has variance."
      # lam_1/lam_2 remain in out_dt as raw passthrough columns.

      if (has_pc1) {
        scores <- as.numeric(dt[wi, ..pc1_cols])
        scores <- scores[is.finite(scores)]
        if (length(scores) >= 20) {
          # Dip test for trimodality
          dip_p <- tryCatch({
            if (requireNamespace("diptest", quietly = TRUE)) {
              diptest::dip.test(scores)$p.value
            } else {
              km <- tryCatch(kmeans(scores, centers = 3, nstart = 3), error = function(e) NULL)
              if (!is.null(km)) {
                centers <- sort(km$centers[, 1])
                gaps <- diff(centers)
                within_sd <- sqrt(mean(km$withinss / km$size))
                if (within_sd > 0) 1 / (1 + (min(gaps) / within_sd)^2) else 1
              } else 1
            }
          }, error = function(e) NA_real_)

          # HET CONTRAST CHECK: if trimodal, does middle band have different
          # het signature than outer bands?
          # For an inversion: middle (HET) samples have dosage ~1 at informative
          # markers -> higher variance on PC1 within the middle group.
          # For family structure: all 3 groups have similar internal variance.
          km3 <- tryCatch(kmeans(scores, centers = 3, nstart = 5), error = function(e) NULL)
          if (!is.null(km3)) {
            co <- order(km3$centers[, 1])
            var_per_group <- vapply(co, function(g) {
              vals <- scores[km3$cluster == g]
              if (length(vals) >= 3) var(vals) else NA_real_
            }, numeric(1))
            # For inversion: middle group (co[2]) has intermediate position
            # but similar variance. Key check: are the 3 groups DISCRETE
            # (low within-group variance vs between-group gap)?
            sizes <- km3$size[co]
            min_frac <- min(sizes) / length(scores)
            # Separation quality: between-group / within-group variance
            between_var <- var(km3$centers[, 1])
            within_var <- mean(var_per_group, na.rm = TRUE)
            if (is.finite(between_var) && is.finite(within_var) && within_var > 0) {
              het_contrast <- between_var / within_var
            }
            # Penalize if smallest group is too small (< 5%)
            if (min_frac < 0.05) het_contrast <- het_contrast * 0.5

            # ── METRICS USING CHROMOSOME-WIDE FIXED BANDS (v8.5 revised) ──
            # chr_bands is computed ONCE before the per-window loop (see below).
            # Here we use it to measure per-window structure against fixed groups.
            # This is NOT per-window k-means. The groups don't change between windows.

            if (exists("chr_bands") && !is.null(chr_bands)) {
              # Map scores to fixed bands (chr_bands is named vector: sample → band)
              score_names <- pc1_cols  # use column names for matching
              sample_ids <- sub("^PC_1_", "", score_names)
              all_scores <- as.numeric(dt[wi, ..pc1_cols])
              names(all_scores) <- sample_ids

              matched <- intersect(names(all_scores), names(chr_bands))
              if (length(matched) >= 20) {
                s <- all_scores[matched]
                b <- chr_bands[matched]

                # Band sizes
                tab <- table(factor(b, levels = 1:3))
                band_sizes <- as.numeric(tab) / sum(tab)

                # Het middle fraction
                het_mid_fraction <- band_sizes[2]

                # Het middle variance on PC1
                mid_s <- s[b == 2]
                het_mid_variance <- if (length(mid_s) >= 3) var(mid_s) else NA_real_

                # PC1 gap: mean PC1 per band, then min gap between adjacent
                band_means <- tapply(s, b, mean, na.rm = TRUE)
                if (length(band_means) >= 2) {
                  sorted_bm <- sort(band_means)
                  het_pc1_gap <- min(diff(sorted_bm))
                }

                # Ashman's D between outer bands
                outer_a <- s[b == 1]; outer_b <- s[b == 3]
                if (length(outer_a) >= 3 && length(outer_b) >= 3) {
                  mu_diff <- abs(mean(outer_a) - mean(outer_b))
                  pooled_sd <- sqrt((var(outer_a) + var(outer_b)) / 2)
                  pc1_bimodality <- if (pooled_sd > 0) mu_diff / pooled_sd else NA_real_
                }

                # Entropy from fixed bands using per-sample distance to band means
                eps_e <- 1e-15
                bm_vec <- as.numeric(band_means[as.character(1:3)])
                if (all(is.finite(bm_vec)) && diff(range(bm_vec)) > 0) {
                  spread <- diff(range(bm_vec)) * 0.01
                  soft <- matrix(0, nrow = length(s), ncol = 3)
                  for (gi in 1:3) soft[, gi] <- 1 / (abs(s - bm_vec[gi]) + spread)
                  soft <- soft / rowSums(soft)

                  per_H <- -rowSums(soft * log(soft + eps_e))
                  local_entropy <- mean(per_H, na.rm = TRUE)
                  local_ena <- exp(local_entropy)

                  per_d12 <- apply(soft, 1, function(p) {
                    ps <- sort(p, decreasing = TRUE); ps[1] - ps[2]
                  })
                  local_delta12 <- mean(per_d12, na.rm = TRUE)
                }
              }
            } else {
              # Fallback: use per-window km3 (less stable but works if chr_bands not set)
              band_sizes <- sizes / sum(sizes)
              het_mid_fraction <- sizes[2] / sum(sizes)
              mid_vals <- scores[km3$cluster == co[2]]
              het_mid_variance <- if (length(mid_vals) >= 3) var(mid_vals) else NA_real_
              sorted_centers <- sort(km3$centers[, 1])
              het_pc1_gap <- min(diff(sorted_centers))
            }
          }
        }
      }

      # ═══════════════════════════════════════════════════════════════
      # NEW BAND-SHAPE DESCRIPTORS (v9.3)
      # ═══════════════════════════════════════════════════════════════
      # These describe the SHAPE of the local PCA cloud, not whether
      # it's "structured" (which is circular at this stage).

      band_discreteness <- NA_real_
      diffuse_score <- NA_real_
      het_intermediacy <- NA_real_
      n_effective_clusters <- NA_real_

      if (!is.null(km3) && has_pc1) {
        scores_all <- as.numeric(dt[wi, ..pc1_cols])
        scores_all <- scores_all[is.finite(scores_all)]
        if (length(scores_all) >= 20) {
          co_local <- order(km3$centers[, 1])
          sizes_local <- km3$size[co_local]
          centers_sorted <- sort(km3$centers[, 1])
          var_per_g <- vapply(co_local, function(g) {
            vals <- scores_all[km3$cluster == g]
            if (length(vals) >= 3) var(vals) else NA_real_
          }, numeric(1))

          # ── band_discreteness ──
          # Ratio of between-group gap to within-group spread.
          # High = crisp 3 clusters (like Candidate 24).
          # Low = diffuse gradient (like Candidate 34).
          mean_within_sd <- sqrt(mean(var_per_g, na.rm = TRUE))
          min_gap <- min(diff(centers_sorted))
          if (is.finite(mean_within_sd) && mean_within_sd > 0) {
            band_discreteness <- min_gap / mean_within_sd
          }

          # ── diffuse_score ──
          # Inverse of trimodality. How gradient-like is the distribution?
          # Uses the dip test p-value: high p = more unimodal/diffuse.
          # Also penalizes if any band is very small (< 8% of samples).
          min_frac_local <- min(sizes_local) / length(scores_all)
          dip_contrib <- if (is.finite(dip_p)) dip_p else 0.5
          size_penalty <- if (min_frac_local < 0.08) 0.3 else 0
          diffuse_score <- pmin(1, dip_contrib + size_penalty)

          # ── het_intermediacy ──
          # Is the middle band truly intermediate between the outers?
          # Measures |mean_het - midpoint(mean_outer1, mean_outer2)| / gap.
          # 0 = perfect intermediacy, 1 = shifted to one side.
          mid_mean <- centers_sorted[2]
          outer_midpoint <- (centers_sorted[1] + centers_sorted[3]) / 2
          full_gap <- centers_sorted[3] - centers_sorted[1]
          if (is.finite(full_gap) && full_gap > 0) {
            het_intermediacy <- 1 - pmin(1, abs(mid_mean - outer_midpoint) / (full_gap / 2))
          }

          # ── n_effective_clusters ──
          # Gap statistic proxy: how many real groups exist?
          # Try k=2,3,4 and pick the k with best silhouette-like metric.
          # Cheap version: compare BSS/TSS for k=2 vs k=3 vs k=4.
          tss <- sum((scores_all - mean(scores_all))^2)
          if (tss > 0) {
            bss_k <- numeric(4)
            bss_k[1] <- 0  # k=1 has 0 BSS
            bss_k[3] <- sum(km3$betweenss) / tss
            for (kk in c(2L, 4L)) {
              kmk <- tryCatch(kmeans(scores_all, centers = kk, nstart = 3),
                              error = function(e) NULL)
              bss_k[kk] <- if (!is.null(kmk)) sum(kmk$betweenss) / tss else 0
            }
            # Effective K: use elbow — where does adding a cluster stop helping?
            # Gain from k-1 to k
            gains <- diff(bss_k)
            # n_effective_clusters = highest k where gain > 0.05
            n_effective_clusters <- 1
            for (kk in 2:4) {
              if (gains[kk - 1] > 0.05) n_effective_clusters <- kk
            }
          }
        }
      }

      # ═══════════════════════════════════════════════════════════════
      # REVISED INV_LIKENESS (v9.3)
      # ═══════════════════════════════════════════════════════════════
      # PVE1 REMOVED from formula — circular (precomp windows are already
      # MDS outliers, so PVE1 just says "this outlier is structured").
      #
      # New formula: 45% HET + 30% TRIMODALITY + 25% BAND DISCRETENESS
      #
      # HET component (45%): adaptive, uses best available source:
      #   Priority 1: dosage_het_rate_cv (genuinely independent)
      #   Priority 2: het_contrast (PCA-derived but informative)
      #   Fallback:   0
      #
      # TRIMODALITY component (30%): dip test
      #
      # BAND DISCRETENESS component (25%): are the 3 groups crisp?
      #   High = discrete clusters = inversion-like
      #   Low = gradient / diffuse = less inversion-like
      # ═══════════════════════════════════════════════════════════════

      # s_dip: trimodality
      s_dip <- if (is.finite(dip_p)) pmin(1, pmax(0, (0.5 - dip_p) / 0.45)) else 0

      # s_het: use het_contrast (PCA-derived) in first pass.
      # Dosage het CV upgrades this in the second pass after dosage merge.
      s_het <- 0
      het_source <- "none"
      if (is.finite(het_contrast)) {
        s_het <- pmin(1, pmax(0, (het_contrast - 2) / 8))
        het_source <- "pc1_het_contrast"
      }

      # s_discrete: band discreteness (new v9.3)
      s_discrete <- 0
      if (is.finite(band_discreteness)) {
        s_discrete <- pmin(1, pmax(0, (band_discreteness - 1) / 4))
      }

      # Combined: 45% het + 30% trimodality + 25% band discreteness
      inv_like <- 0.45 * s_het + 0.30 * s_dip + 0.25 * s_discrete

      # ═══════════════════════════════════════════════════════════════
      # TWO-TRACK DECOMPOSITION (v9.3.1)
      # ═══════════════════════════════════════════════════════════════
      # Track 1: INV LIKENESS — composite (het + trimodality + discreteness)
      #   Already computed as inv_like above.
      #   Per-chromosome calibration via Beta adaptive threshold (test_07)
      #   is applied after the per-window loop.

      # Track 2: FAMILY LIKENESS — do k=3 bands correlate with Q ancestry?
      #
      #   DIAGNOSTIC ANNOTATION ONLY — not a gate or evidence.
      #
      #   High = bands recapitulate genome-wide family/ancestry groups.
      #   Low  = bands crosscut ancestry groups.
      #
      #   IMPORTANT (hatchery biology): high family_likeness does NOT mean
      #   "this is a family artifact." In a hatchery with ~20 founders, a
      #   real inversion that entered through one founder lineage will show
      #   high family_likeness because the carrier samples ARE one ancestry
      #   component. The metric flags the confound but cannot resolve it.
      #
      #   family_likeness is stored for downstream diagnostic review and
      #   PA summary aggregation. It should NOT be used as:
      #     - a seed gate (would kill founder-linked rare inversions)
      #     - a negative weight in inv_likeness or scoring
      #     - evidence for H1 (family structure) without concordance from
      #       T1 (relatedness) + T3 (kin-pruned retention) + T2 (eff_K)
      #
      # Method: one-way ANOVA on Q vectors grouped by per-window k=3 bands.
      # Score = between-band Q variance / total Q variance
      #       = fraction of Q variance explained by the k=3 bands
      family_like <- NA_real_
      family_q_between <- NA_real_
      family_q_within_entropy <- NA_real_

      if (!is.null(q_mat) && !is.null(km3) && has_pc1) {
        # Get per-window k=3 band assignments
        # Use per-window km3 (NOT fixed bands) because we're asking:
        # "does THIS window's grouping track ancestry?"
        all_scores_fw <- as.numeric(dt[wi, ..pc1_cols])
        names(all_scores_fw) <- sub("^PC_1_", "", pc1_cols)
        valid_fw <- is.finite(all_scores_fw)
        samps_fw <- names(all_scores_fw)[valid_fw]

        matched_q <- intersect(samps_fw, rownames(q_mat))
        if (length(matched_q) >= 20 && length(co) == 3) {
          # Assign bands using per-window km3 cluster assignments
          # km3$cluster is indexed by position in scores vector
          # scores was built from dt[wi, ..pc1_cols] filtered by is.finite
          fw_scores <- all_scores_fw[valid_fw]
          fw_names <- names(fw_scores)

          # Predict band for each sample using nearest km3 center
          band_assign <- integer(length(matched_q))
          names(band_assign) <- matched_q
          for (qi in seq_along(matched_q)) {
            s_val <- all_scores_fw[matched_q[qi]]
            if (!is.finite(s_val)) { band_assign[qi] <- NA; next }
            dists_to_centers <- abs(s_val - km3$centers[co, 1])
            band_assign[qi] <- which.min(dists_to_centers)
          }
          band_assign <- band_assign[!is.na(band_assign)]

          if (length(band_assign) >= 15 && length(unique(band_assign)) >= 2) {
            q_sub <- q_mat[names(band_assign), , drop = FALSE]

            # Between-band Q variance:
            # For each Q component, compute variance of band means
            band_ids <- band_assign
            grand_mean_q <- colMeans(q_sub, na.rm = TRUE)
            total_var_q <- sum(apply(q_sub, 2, var, na.rm = TRUE))

            if (total_var_q > 1e-10) {
              between_var_q <- 0
              for (bi_q in unique(band_ids)) {
                bi_samps <- which(band_ids == bi_q)
                if (length(bi_samps) < 2) next
                bi_mean <- colMeans(q_sub[bi_samps, , drop = FALSE], na.rm = TRUE)
                between_var_q <- between_var_q +
                  length(bi_samps) * sum((bi_mean - grand_mean_q)^2)
              }
              between_var_q <- between_var_q / nrow(q_sub)
              family_q_between <- between_var_q

              # Family likeness = fraction of Q variance explained by bands
              # High → bands track ancestry → family-driven
              family_like <- pmin(1, between_var_q / total_var_q)
            }

            # Within-band Q entropy: for each band, compute entropy of mean Q
            band_entropies <- vapply(unique(band_ids), function(bi_q) {
              bi_samps <- which(band_ids == bi_q)
              if (length(bi_samps) < 2) return(NA_real_)
              bi_mean_q <- colMeans(q_sub[bi_samps, , drop = FALSE], na.rm = TRUE)
              bi_mean_q <- bi_mean_q[bi_mean_q > 0]
              if (length(bi_mean_q) < 2) return(0)
              -sum(bi_mean_q * log(bi_mean_q + 1e-15))
            }, numeric(1))
            family_q_within_entropy <- mean(band_entropies, na.rm = TRUE)
          }
        }
      }

      # ═══════════════════════════════════════════════════════════════
      # TEST 05 Q-STAMP: per-band ancestry composition (v9.3.2)
      # ═══════════════════════════════════════════════════════════════
      # For each window's k=3 bands, stamp which Q components dominate
      # each band. This is a DESCRIPTIVE ANNOTATION — not a score or gate.
      #
      # Useful for: visual inspection, manuscript figures, downstream
      # scripts that want to know "what color is each band" without
      # recomputing from Q-matrix + PC1 loadings.
      #
      # NOTE: genome-wide Q (NGSadmix K=8) is used. On chromosomes with
      # big inversions, Q may partially capture inversion signal (circular).
      # This is accepted — stamps are diagnostic, not evidence.
      #
      # Columns added (8):
      #   band1_dom_Q, band1_dom_Q_frac   — dominant Q component in band 1
      #   band2_dom_Q, band2_dom_Q_frac   — dominant Q component in band 2
      #   band3_dom_Q, band3_dom_Q_frac   — dominant Q component in band 3
      #   bands_same_Q    — do band1 and band3 share the same dominant Q?
      #   band_Q_entropy  — Shannon entropy of Q composition across bands
      #                     (low = one Q dominates = family-like)
      # ═══════════════════════════════════════════════════════════════

      band1_dom_Q <- NA_integer_; band1_dom_Q_frac <- NA_real_
      band2_dom_Q <- NA_integer_; band2_dom_Q_frac <- NA_real_
      band3_dom_Q <- NA_integer_; band3_dom_Q_frac <- NA_real_
      bands_same_Q <- NA
      band_Q_entropy <- NA_real_

      if (!is.null(q_mat) && exists("band_assign") && length(band_assign) >= 15) {
        # q_sub is already computed above (samples × K)
        # band_assign maps sample → band (1/2/3)
        for (bi in 1:3) {
          bi_samps <- names(band_assign)[band_assign == bi]
          bi_samps_q <- intersect(bi_samps, rownames(q_mat))
          if (length(bi_samps_q) >= 2) {
            mean_q <- colMeans(q_mat[bi_samps_q, , drop = FALSE], na.rm = TRUE)
            dom_idx <- which.max(mean_q)
            dom_frac <- mean_q[dom_idx]
            if (bi == 1) { band1_dom_Q <- dom_idx; band1_dom_Q_frac <- dom_frac }
            if (bi == 2) { band2_dom_Q <- dom_idx; band2_dom_Q_frac <- dom_frac }
            if (bi == 3) { band3_dom_Q <- dom_idx; band3_dom_Q_frac <- dom_frac }
          }
        }

        # Do band1 and band3 share the same dominant Q?
        if (is.finite(band1_dom_Q) && is.finite(band3_dom_Q)) {
          bands_same_Q <- band1_dom_Q == band3_dom_Q
        }

        # Cross-band Q entropy: how concentrated is ancestry across the 3 bands?
        # Compute the mean Q vector per band, then measure how different they are.
        # If all 3 bands have the same Q profile → low entropy → family.
        # If bands have different Q profiles → high entropy → structural.
        band_q_vecs <- list()
        for (bi in 1:3) {
          bi_samps <- names(band_assign)[band_assign == bi]
          bi_samps_q <- intersect(bi_samps, rownames(q_mat))
          if (length(bi_samps_q) >= 2) {
            band_q_vecs[[bi]] <- colMeans(q_mat[bi_samps_q, , drop = FALSE], na.rm = TRUE)
          }
        }
        if (length(band_q_vecs) >= 2) {
          # Stack band-mean Q vectors and compute row-wise entropy of the
          # "which band does each Q component belong to most?" distribution
          bq_mat <- do.call(rbind, band_q_vecs)  # n_bands × K
          # Normalize columns to sum to 1 (distribution of each Q across bands)
          col_sums <- colSums(bq_mat)
          col_sums[col_sums == 0] <- 1
          bq_norm <- t(t(bq_mat) / col_sums)
          # Per-Q-component entropy across bands
          eps_q <- 1e-15
          per_q_H <- apply(bq_norm, 2, function(p) {
            p <- p[p > eps_q]
            -sum(p * log(p))
          })
          # Mean entropy: high = Q spread evenly across bands = not family
          # Low = Q concentrated in one band = family-like
          band_Q_entropy <- mean(per_q_H, na.rm = TRUE)
        }
      }

      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
        # PVE1 columns — DEPRECATED v9.3 (circular). Kept as NA for compat.
        inv_pve1 = NA_real_, inv_eigen_ratio = NA_real_,
        inv_dip_p = round(dip_p, 4),
        inv_het_contrast = round(het_contrast, 2),
        inv_pve1_excess = NA_real_,
        inv_likeness = round(inv_like, 4),
        # ── BAND-SHAPE DESCRIPTORS (v9.3) ──
        band_discreteness = round(band_discreteness %||% NA_real_, 4),
        diffuse_score = round(diffuse_score %||% NA_real_, 4),
        het_intermediacy = round(het_intermediacy %||% NA_real_, 4),
        n_effective_clusters = as.integer(n_effective_clusters %||% NA_integer_),
        # ── THREE-TRACK DECOMPOSITION (v9.3) ──
        # structure_likeness removed (v9.3.1) — was circular PVE1, always NA
        family_likeness = round(family_like %||% NA_real_, 4),
        family_q_between = round(family_q_between %||% NA_real_, 6),
        family_q_within_entropy = round(family_q_within_entropy %||% NA_real_, 4),
        # ── TEST 05 Q-STAMP (v9.3.2) ──
        band1_dom_Q = as.integer(band1_dom_Q %||% NA_integer_),
        band1_dom_Q_frac = round(band1_dom_Q_frac %||% NA_real_, 4),
        band2_dom_Q = as.integer(band2_dom_Q %||% NA_integer_),
        band2_dom_Q_frac = round(band2_dom_Q_frac %||% NA_real_, 4),
        band3_dom_Q = as.integer(band3_dom_Q %||% NA_integer_),
        band3_dom_Q_frac = round(band3_dom_Q_frac %||% NA_real_, 4),
        bands_same_Q = as.logical(bands_same_Q %||% NA),
        band_Q_entropy = round(band_Q_entropy %||% NA_real_, 4),
        # Multi-source het assessment
        pc1_bimodality = round(pc1_bimodality %||% NA_real_, 4),
        het_pc1_gap = round(het_pc1_gap %||% NA_real_, 4),
        het_mid_fraction = round(het_mid_fraction %||% NA_real_, 4),
        het_mid_variance = round(het_mid_variance %||% NA_real_, 6),
        # Local structure metrics
        local_delta12 = round(local_delta12 %||% NA_real_, 4),
        local_entropy = round(local_entropy %||% NA_real_, 4),
        local_ena = round(local_ena %||% NA_real_, 4),
        # Band sizes for downstream stratification
        band1_frac = round(band_sizes[1] %||% NA_real_, 4),
        band2_frac = round(band_sizes[2] %||% NA_real_, 4),
        band3_frac = round(band_sizes[3] %||% NA_real_, 4),
        # Dosage-based het rate: genuinely independent from PCA
        dosage_het_rate_median = round(dosage_het_rate_median %||% NA_real_, 4),
        dosage_het_rate_sd = round(dosage_het_rate_sd %||% NA_real_, 4),
        dosage_het_rate_cv = round(dosage_het_rate_cv %||% NA_real_, 4)
      )
    }
  }
  if (length(all_rows) > 0) rbindlist(all_rows) else data.table()
}

message("[PRECOMP] Computing inversion-likeness scores...")
t_il <- proc.time()
inv_like_dt <- compute_inv_likeness_all(per_chr, chroms)
message("[PRECOMP] Inv-likeness: ", nrow(inv_like_dt), " windows, ",
        sum(inv_like_dt$inv_likeness >= 0.5, na.rm = TRUE), " with score >= 0.5 (",
        round((proc.time() - t_il)[3], 1), "s)")

# ── Dosage-based het rate (independent from PCA) ──
dosage_het_dt <- data.table()
if (!is.null(dosage_dir)) {
  message("[PRECOMP] Computing dosage het rates from ", dosage_dir, " ...")
  t_dos <- proc.time()
  dosage_het_dt <- compute_dosage_het_rates(dosage_dir, per_chr, chroms)
  message("[PRECOMP] Dosage het: ", nrow(dosage_het_dt), " windows (",
          round((proc.time() - t_dos)[3], 1), "s)")

  # Merge into inv_like_dt using clean data.table join
  if (nrow(dosage_het_dt) > 0) {
    inv_like_dt[dosage_het_dt, on = c("chrom", "global_window_id"),
      `:=`(dosage_het_rate_median = i.dosage_het_rate_median,
           dosage_het_rate_sd = i.dosage_het_rate_sd,
           dosage_het_rate_cv = i.dosage_het_rate_cv)]

    n_filled <- sum(is.finite(inv_like_dt$dosage_het_rate_cv))
    message("[PRECOMP] Dosage merge: ", n_filled, " windows got dosage values")

    # SECOND PASS: recompute inv_likeness for windows that now have dosage CV
    # Dosage CV replaces het_contrast as the het component (45% weight)
    # v9.3: no PVE1, uses band_discreteness instead
    n_with_dosage <- sum(is.finite(inv_like_dt$dosage_het_rate_cv) &
                          inv_like_dt$dosage_het_rate_cv > 0)
    if (n_with_dosage > 0) {
      message("[PRECOMP] Upgrading inv_likeness for ", n_with_dosage,
              " windows with dosage het CV...")
      inv_like_dt[is.finite(dosage_het_rate_cv) & dosage_het_rate_cv > 0,
        inv_likeness := {
          s_h <- pmin(1, pmax(0, (dosage_het_rate_cv - 0.1) / 0.6))
          s_d <- pmin(1, pmax(0, (0.5 - inv_dip_p) / 0.45))
          s_d <- ifelse(is.finite(s_d), s_d, 0)
          s_bd <- pmin(1, pmax(0, (band_discreteness - 1) / 4))
          s_bd <- ifelse(is.finite(s_bd), s_bd, 0)
          round(0.45 * s_h + 0.30 * s_d + 0.25 * s_bd, 4)
        }]
      message("[PRECOMP] Done: ", n_with_dosage, " windows upgraded to dosage-based inv_likeness")
    }
  }
} else {
  message("[PRECOMP] No --dosage_dir provided — dosage het rate columns will be NA")
}

# ── Track summary (v9.3) ──
n_inv <- sum(inv_like_dt$inv_likeness >= 0.5, na.rm = TRUE)
n_fam <- sum(inv_like_dt$family_likeness >= 0.5, na.rm = TRUE)
n_dos <- sum(is.finite(inv_like_dt$dosage_het_rate_median))
n_discrete <- sum(inv_like_dt$band_discreteness > 2, na.rm = TRUE)
n_diffuse <- sum(inv_like_dt$diffuse_score > 0.5, na.rm = TRUE)
message("[PRECOMP] ── TRACK SUMMARY (v9.3) ──")
message("  Inv-like >= 0.5:       ", n_inv, " / ", nrow(inv_like_dt),
        "  (formula: 45% het + 30% trimodality + 25% band_discreteness)")
message("  Band discrete (>2):    ", n_discrete, " / ", nrow(inv_like_dt))
message("  Diffuse (>0.5):        ", n_diffuse, " / ", nrow(inv_like_dt))
message("  Family >= 0.5:         ", n_fam, " / ", nrow(inv_like_dt),
        if (is.null(q_mat)) " (Q not loaded)" else " (diagnostic only)")
message("  Dosage het computed:   ", n_dos, " / ", nrow(inv_like_dt))

# Q-stamp summary (test_05)
if ("band1_dom_Q" %in% names(inv_like_dt)) {
  n_stamped <- sum(is.finite(inv_like_dt$band1_dom_Q))
  n_same_q <- sum(inv_like_dt$bands_same_Q == TRUE, na.rm = TRUE)
  message("  Q-stamps computed:     ", n_stamped, " / ", nrow(inv_like_dt),
          if (is.null(q_mat)) " (Q not loaded)" else "")
  if (n_stamped > 0) {
    message("  Bands same dominant Q: ", n_same_q, " / ", n_stamped,
            " (", round(100 * n_same_q / n_stamped, 1), "% — high = family structure)")
    med_ent <- median(inv_like_dt$band_Q_entropy[is.finite(inv_like_dt$band_Q_entropy)])
    message("  Band Q entropy median: ", round(med_ent, 3),
            " (low = family-like, high = mixed ancestry per band)")
  }
}

# NOTE: family_likeness is a diagnostic annotation, NOT a gate.
# In a hatchery population, real inversions linked to one founder
# lineage will show high family_likeness. The metric flags the
# confound but cannot distinguish artifact from founder-linked
# inversion. Resolution requires hypothesis tests (T1/T2/T3)
# at the candidate scoring stage.

# ═══════════════════════════════════════════════════════════════════
# TEST 07: Beta adaptive threshold (integrated v9.3.1)
# ═══════════════════════════════════════════════════════════════════
# Fit Beta(α,β) to per-chromosome inv_likeness distribution.
# Background windows follow the bulk of the Beta. Inversion windows
# sit in the extreme right tail. The p-value gives a chromosome-
# specific threshold: noisy chromosomes get higher effective
# thresholds, clean chromosomes catch weaker signals.
#
# Adds columns: beta_pval, adaptive_seed, beta_alpha, beta_beta
# Snake 1 can use adaptive_seed == TRUE instead of fixed threshold.
# ═══════════════════════════════════════════════════════════════════

beta_adaptive_pvalues <- function(scores) {
  valid <- scores[is.finite(scores) & scores > 0.001 & scores < 0.999]
  if (length(valid) < 50) {
    return(data.table(beta_pval = rep(NA_real_, length(scores)),
                       adaptive_seed = rep(FALSE, length(scores)),
                       beta_alpha = NA_real_, beta_beta = NA_real_))
  }
  fit <- tryCatch({
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::fitdistr(valid, "beta", start = list(shape1 = 1, shape2 = 5))
    } else {
      m <- mean(valid); v <- var(valid)
      if (v >= m * (1 - m)) v <- m * (1 - m) * 0.99
      a <- m * (m * (1 - m) / v - 1)
      b <- (1 - m) * (m * (1 - m) / v - 1)
      list(estimate = c(shape1 = max(0.1, a), shape2 = max(0.1, b)))
    }
  }, error = function(e) list(estimate = c(shape1 = 1, shape2 = 5)))
  a <- fit$estimate["shape1"]; b <- fit$estimate["shape2"]
  pvals <- 1 - pbeta(scores, a, b)
  pvals[!is.finite(scores)] <- NA_real_
  data.table(beta_pval = pvals, adaptive_seed = pvals < 0.01 & is.finite(pvals),
             beta_alpha = round(a, 3), beta_beta = round(b, 3))
}

# Apply per chromosome
message("[PRECOMP] ── TEST 07: Beta adaptive threshold ──")
inv_like_dt[, c("beta_pval", "adaptive_seed", "beta_alpha", "beta_beta") :=
              list(NA_real_, FALSE, NA_real_, NA_real_)]

for (chr_i in unique(inv_like_dt$chrom)) {
  idx <- inv_like_dt$chrom == chr_i
  scores_chr <- inv_like_dt$inv_likeness[idx]
  bt <- beta_adaptive_pvalues(scores_chr)
  inv_like_dt[idx, `:=`(
    beta_pval = bt$beta_pval,
    adaptive_seed = bt$adaptive_seed,
    beta_alpha = bt$beta_alpha,
    beta_beta = bt$beta_beta
  )]
  n_seeds <- sum(bt$adaptive_seed, na.rm = TRUE)
  message("  ", chr_i, ": Beta(", bt$beta_alpha[1], ",", bt$beta_beta[1],
          ") → ", n_seeds, " adaptive seeds (p<0.01)")
}

n_total_seeds <- sum(inv_like_dt$adaptive_seed, na.rm = TRUE)
message("[PRECOMP] Total adaptive seeds: ", n_total_seeds, " / ", nrow(inv_like_dt))

f_inv <- file.path(outdir, "window_inv_likeness.tsv.gz")

# ── Merge local Q from ancestry_bridge (if cached) ──
# Local Q provides REAL per-window ancestry proportions (delta12, H, ENA)
# computed via fixed-F EM on BEAGLE genotype likelihoods.
# These replace the PC1-proxy metrics with genuine admixture estimates.
# Run ancestry_bridge.R --prepare FIRST to generate the cache.
local_q_dir <- file.path(outdir, "local_Q")
ancestry_bridge_path <- file.path(dirname(dirname(f_inv)), "utils", "ancestry_bridge.R")
if (!file.exists(ancestry_bridge_path)) {
  ancestry_bridge_path <- Sys.glob(file.path(dirname(f_inv), "..", "inversion_codebase_*/utils/ancestry_bridge.R"))[1]
}

if (!is.null(ancestry_bridge_path) && file.exists(ancestry_bridge_path) &&
    dir.exists(local_q_dir)) {
  tryCatch({
    source(ancestry_bridge_path)
    inv_like_dt <- merge_local_Q_into_invlikeness(inv_like_dt, cache_dir = local_q_dir)
    message("[PRECOMP] Local Q columns added: localQ_delta12, localQ_entropy, localQ_ena")
  }, error = function(e) {
    message("[PRECOMP] ancestry_bridge merge failed: ", e$message, " — continuing without")
  })
} else {
  message("[PRECOMP] No local Q cache found at ", local_q_dir,
          " — run ancestry_bridge.R --prepare first to enable real Q metrics")
}

# ═══════════════════════════════════════════════════════════════════
# TEST 05: REAL Fst SCAN via Engine B dispatcher (v9.3.4)
# ═══════════════════════════════════════════════════════════════════
# Per-window family_fst_ratio: Fst(PC1_high, PC1_low) vs Fst(Q_high, Q_low).
#
# This COMPLEMENTS the ANOVA family_likeness (which measures Q variance
# explained by bands). The Fst ratio measures actual allele frequency
# divergence — a window can have high ANOVA family_likeness but low Fst
# if the ancestry groups aren't differentiated at that locus.
#
# family_fst_ratio ~1 → PC1 bands = ancestry grouping (family-driven)
# family_fst_ratio ~0 → PC1 captures something orthogonal to Q (inversion)
#
# Requires: load_bridge.R → get_region_stats(what="Fst", groups=...)
# Cost: ~2 calls per window × ~10,000 windows/chr ≈ 2-3 min with dosage cache
# ═══════════════════════════════════════════════════════════════════

if (.bridge_available && exists("get_region_stats", mode = "function") &&
    !is.null(q_mat) && K_ancestry >= 2) {
  message("[PRECOMP] ── TEST 05: Real Fst scan (Engine B) ──")

  # Helper: extract Fst from dispatcher result
  extract_fst_c5 <- function(s) {
    if (is.null(s$Fst) || length(s$Fst) == 0) return(NA_real_)
    as.numeric(s$Fst[[1]])
  }

  # Initialize columns
  inv_like_dt[, `:=`(
    test05_fst_pc1 = NA_real_,
    test05_fst_q_best = NA_real_,
    test05_fst_q_best_k = NA_integer_,
    test05_family_fst_ratio = NA_real_
  )]

  Q_SPLIT <- 0.30  # top/bottom 30% for Q-axis grouping

  for (chr_i in unique(inv_like_dt$chrom)) {
    chr_rows <- which(inv_like_dt$chrom == chr_i)
    if (length(chr_rows) < 10) next

    # Get per-window band assignments for this chr from precomp
    pc_obj <- NULL
    pf <- file.path(precomp_dir, paste0(chr_i, ".precomp.rds"))
    if (file.exists(pf)) pc_obj <- readRDS(pf)
    if (is.null(pc_obj)) next

    dt_chr <- pc_obj$dt
    pc1_cols_chr <- grep("^PC_1_", names(dt_chr), value = TRUE)
    ind_names_chr <- sub("^PC_1_", "", pc1_cols_chr)

    # Map Ind→CGA for dispatcher groups
    if (!is.null(smap) && grepl("^Ind[0-9]", ind_names_chr[1])) {
      cga_names <- smap$to_real_vec(ind_names_chr)
    } else if (!is.null(real_names) && grepl("^Ind[0-9]", ind_names_chr[1]) &&
               length(real_names) == length(ind_names_chr)) {
      cga_names <- real_names
    } else {
      cga_names <- ind_names_chr  # already CGA or no mapping needed
    }
    names(cga_names) <- ind_names_chr

    # Subsample windows for speed: every 5th window
    scan_idx <- seq(1, length(chr_rows), by = 5L)
    n_done <- 0L; n_fail <- 0L

    for (si in scan_idx) {
      wi <- chr_rows[si]
      s_bp <- inv_like_dt$start_bp[wi]
      e_bp <- inv_like_dt$end_bp[wi]

      # PC1 bands for this window
      pc1_vals <- as.numeric(dt_chr[si, ..pc1_cols_chr])
      names(pc1_vals) <- ind_names_chr
      pc1_vals <- pc1_vals[is.finite(pc1_vals)]
      if (length(pc1_vals) < 30) next

      # k=3 on PC1
      km3 <- tryCatch(kmeans(pc1_vals, centers = 3, nstart = 5), error = function(e) NULL)
      if (is.null(km3)) next
      co <- order(km3$centers[, 1])
      band_assign <- integer(length(pc1_vals))
      names(band_assign) <- names(pc1_vals)
      for (bi in 1:3) band_assign[km3$cluster == co[bi]] <- bi

      b1_cga <- cga_names[names(band_assign[band_assign == 1])]
      b3_cga <- cga_names[names(band_assign[band_assign == 3])]
      b1_cga <- b1_cga[!is.na(b1_cga)]
      b3_cga <- b3_cga[!is.na(b3_cga)]
      if (length(b1_cga) < 5 || length(b3_cga) < 5) next

      # Fst(PC1_band1, PC1_band3)
      fst_pc1 <- tryCatch({
        s <- get_region_stats(chr_i, s_bp, e_bp, what = "Fst",
                               groups = list(b1 = b1_cga, b3 = b3_cga))
        extract_fst_c5(s)
      }, error = function(e) NA_real_)

      # Fst(Q_high, Q_low) for each Q component — find best
      q_fsts <- numeric(K_ancestry)
      for (k in seq_len(K_ancestry)) {
        q_vals <- q_mat[, k]
        names(q_vals) <- rownames(q_mat)
        q_lo <- quantile(q_vals, Q_SPLIT, na.rm = TRUE)
        q_hi <- quantile(q_vals, 1 - Q_SPLIT, na.rm = TRUE)

        # Map to CGA
        high_ind <- names(q_vals)[q_vals >= q_hi]
        low_ind  <- names(q_vals)[q_vals <= q_lo]
        high_cga <- cga_names[high_ind]; high_cga <- high_cga[!is.na(high_cga)]
        low_cga  <- cga_names[low_ind];  low_cga  <- low_cga[!is.na(low_cga)]

        if (length(high_cga) < 5 || length(low_cga) < 5) { q_fsts[k] <- NA; next }

        q_fsts[k] <- tryCatch({
          s <- get_region_stats(chr_i, s_bp, e_bp, what = "Fst",
                                 groups = list(hi = high_cga, lo = low_cga))
          extract_fst_c5(s)
        }, error = function(e) NA_real_)
      }

      # Best Q Fst and ratio
      fst_q_best <- max(q_fsts, na.rm = TRUE)
      fst_q_best_k <- which.max(q_fsts)
      if (!is.finite(fst_q_best)) fst_q_best <- NA_real_

      ratio <- if (is.finite(fst_pc1) && fst_pc1 > 0 && is.finite(fst_q_best)) {
        pmin(2, pmax(0, fst_q_best / fst_pc1))
      } else NA_real_

      inv_like_dt[wi, `:=`(
        test05_fst_pc1 = round(fst_pc1 %||% NA_real_, 4),
        test05_fst_q_best = round(fst_q_best %||% NA_real_, 4),
        test05_fst_q_best_k = as.integer(fst_q_best_k %||% NA_integer_),
        test05_family_fst_ratio = round(ratio %||% NA_real_, 4)
      )]
      n_done <- n_done + 1L
    }

    # Interpolate to non-scanned windows (nearest-neighbor fill)
    scanned <- chr_rows[scan_idx]
    scanned_valid <- scanned[!is.na(inv_like_dt$test05_fst_pc1[scanned])]
    if (length(scanned_valid) > 0) {
      for (col in c("test05_fst_pc1", "test05_fst_q_best",
                     "test05_fst_q_best_k", "test05_family_fst_ratio")) {
        vals <- inv_like_dt[[col]][scanned_valid]
        mids_scanned <- (inv_like_dt$start_bp[scanned_valid] + inv_like_dt$end_bp[scanned_valid]) / 2
        mids_all <- (inv_like_dt$start_bp[chr_rows] + inv_like_dt$end_bp[chr_rows]) / 2
        nn_idx <- findInterval(mids_all, mids_scanned)
        nn_idx[nn_idx == 0] <- 1L
        nn_idx[nn_idx > length(mids_scanned)] <- length(mids_scanned)
        set(inv_like_dt, i = chr_rows, j = col, value = vals[nn_idx])
      }
    }

    message("  ", chr_i, ": ", n_done, " windows scanned (every 5th), interpolated to ", length(chr_rows))
  }

  n_ratio <- sum(is.finite(inv_like_dt$test05_family_fst_ratio))
  n_inversion_like <- sum(inv_like_dt$test05_family_fst_ratio < 0.3, na.rm = TRUE)
  n_family_like <- sum(inv_like_dt$test05_family_fst_ratio > 0.7, na.rm = TRUE)
  message("[PRECOMP] test_05 Fst scan: ", n_ratio, " windows with ratio, ",
          n_inversion_like, " inversion-like (ratio<0.3), ",
          n_family_like, " family-like (ratio>0.7)")
} else {
  message("[PRECOMP] test_05 Fst scan skipped (need Engine B + Q matrix)")
  # Ensure columns exist even when skipped
  if (!"test05_fst_pc1" %in% names(inv_like_dt)) {
    inv_like_dt[, `:=`(
      test05_fst_pc1 = NA_real_,
      test05_fst_q_best = NA_real_,
      test05_fst_q_best_k = NA_integer_,
      test05_family_fst_ratio = NA_real_
    )]
  }
}

# ═══════════════════════════════════════════════════════════════════
# SV PRIOR ANNOTATION — integrated v9.3.1
# ═══════════════════════════════════════════════════════════════════
# Last annotation step — stamps SV columns onto inv_like_dt.
# Reads per-chromosome sv_prior RDS built by STEP_C00_build_sv_prior.R.
# These are PURE ANNOTATION — they do not alter inv_likeness or any
# upstream computation. Can be added/rerun anytime after inv_like_dt exists.
# Downstream scripts (merge blocker, scoring, decomposition seeds) consume them.
#
# Evidence tests stamped here (see phase_4_catalog/ for the test catalogue):
#   Test 01: SV INV call overlap (caller + confidence + AF)
#   Test 02: het-DEL at ±BP_MATCH_WINDOW around INV breakpoints
#   Test 03: het-DEL internal to the INV
#   Test 08: BND-paired triangulation — rescues INVs misclassified as BND
# ═══════════════════════════════════════════════════════════════════

inv_like_dt[, `:=`(
  sv_inv_overlap    = 0L,
  sv_inv_confidence = NA_character_,
  sv_inv_af         = NA_real_,
  sv_het_del_count  = 0L,
  sv_n_anchors      = 0L
)]

if (!is.null(sv_prior_dir) && dir.exists(sv_prior_dir)) {
  message("[PRECOMP] ── SV PRIOR annotation (tests 01/02/03/08) ──")
  n_annotated <- 0L

  for (chr in unique(inv_like_dt$chrom)) {
    # Try current name, then legacy names for backward compatibility
    fl_file <- file.path(sv_prior_dir, paste0("sv_prior_", chr, ".rds"))
    if (!file.exists(fl_file)) {
      fl_file <- file.path(sv_prior_dir, paste0("sv_flashlight_", chr, ".rds"))
    }
    if (!file.exists(fl_file)) {
      fl_file <- file.path(sv_prior_dir, paste0(chr, "_flashlight.rds"))
    }
    if (!file.exists(fl_file)) next

    fl <- tryCatch(readRDS(fl_file), error = function(e) NULL)
    if (is.null(fl)) next

    chr_idx <- which(inv_like_dt$chrom == chr)
    if (length(chr_idx) == 0) next

    inv_calls <- if ("inv_calls" %in% names(fl)) fl$inv_calls else NULL
    bp_dels   <- if ("breakpoint_dels" %in% names(fl)) fl$breakpoint_dels else NULL
    int_dels  <- if ("internal_dels" %in% names(fl)) fl$internal_dels else NULL
    samp_inv  <- if ("sample_inv_states" %in% names(fl)) fl$sample_inv_states else NULL

    # Test 08: BND-triangulated inversions (extend inv_calls with BND pairs)
    bnd_tri <- if ("bnd_triangulated" %in% names(fl)) fl$bnd_triangulated else NULL
    if (!is.null(bnd_tri) && nrow(bnd_tri) > 0) {
      # Convert BND pairs to inv_calls-like format for unified processing
      bnd_as_inv <- data.table(
        inv_id = paste0("BND_", bnd_tri$bnd_left_id),
        chrom = bnd_tri$chrom,
        bp1 = bnd_tri$bp1, bp2 = bnd_tri$bp2,
        svlen = bnd_tri$svlen,
        caller = "bnd_triangulation",
        confidence_level = fifelse(bnd_tri$carrier_jaccard > 0.6, "HIGH",
                            fifelse(bnd_tri$carrier_jaccard > 0.3, "MEDIUM", "LOW")),
        n_carriers = bnd_tri$n_shared_carriers,
        af = round(bnd_tri$n_shared_carriers / n_samples, 3)
      )
      if (is.null(inv_calls) || nrow(inv_calls) == 0) {
        inv_calls <- bnd_as_inv
      } else {
        # Append BND inversions that don't overlap existing INV calls
        for (bi in seq_len(nrow(bnd_as_inv))) {
          b <- bnd_as_inv[bi]
          overlap <- inv_calls[bp1 <= b$bp2 + 50000 & bp2 >= b$bp1 - 50000]
          if (nrow(overlap) == 0) {
            inv_calls <- rbindlist(list(inv_calls, bnd_as_inv[bi]), fill = TRUE)
          }
        }
      }
      message("  ", chr, ": test_08 added ", nrow(bnd_tri), " BND-triangulated inversions")
    }

    if (is.null(inv_calls) || nrow(inv_calls) == 0) next

    for (wi in chr_idx) {
      w_start <- inv_like_dt$start_bp[wi]
      w_end   <- inv_like_dt$end_bp[wi]

      # Test 01: SV INV overlap
      overlapping <- inv_calls[bp1 <= w_end & bp2 >= w_start]
      if (nrow(overlapping) > 0) {
        conf_order <- c("VERY_HIGH", "HIGH", "MEDIUM", "LOW", "MINIMAL", "UNKNOWN")
        best_idx <- 1L
        if ("confidence_level" %in% names(overlapping)) {
          ranks <- match(overlapping$confidence_level, conf_order)
          ranks[is.na(ranks)] <- length(conf_order)
          best_idx <- which.min(ranks)
        }
        set(inv_like_dt, wi, "sv_inv_overlap", 1L)
        set(inv_like_dt, wi, "sv_inv_confidence", overlapping$confidence_level[best_idx])
        if ("af" %in% names(overlapping)) {
          set(inv_like_dt, wi, "sv_inv_af", overlapping$af[best_idx])
        }
        if (!is.null(samp_inv) && nrow(samp_inv) > 0) {
          anchors <- samp_inv[
            inv_id %in% overlapping$inv_id &
            sv_genotype != "MISSING" &
            sv_confidence %in% c("HIGH", "MEDIUM")
          ]
          set(inv_like_dt, wi, "sv_n_anchors", nrow(anchors))
        }
        n_annotated <- n_annotated + 1L
      }

      # Test 02: Het-DELs at breakpoints
      if (!is.null(bp_dels) && nrow(bp_dels) > 0) {
        bp_dels_here <- bp_dels[bp_pos >= w_start & bp_pos <= w_end]
        if (nrow(bp_dels_here) > 0) {
          set(inv_like_dt, wi, "sv_het_del_count", nrow(bp_dels_here))
        }
      }

      # Test 03: Internal het-DELs (hemizygous segments)
      if (!is.null(int_dels) && nrow(int_dels) > 0) {
        int_dels_here <- int_dels[del_start >= w_start & del_end <= w_end]
        if (nrow(int_dels_here) > 0) {
          old_count <- inv_like_dt$sv_het_del_count[wi]
          set(inv_like_dt, wi, "sv_het_del_count", old_count + nrow(int_dels_here))
        }
      }
    }
    message("  ", chr, ": ", sum(inv_like_dt$sv_inv_overlap[chr_idx] == 1),
            " windows with SV INV overlap")
  }

  n_with_sv <- sum(inv_like_dt$sv_inv_overlap == 1, na.rm = TRUE)
  n_with_del <- sum(inv_like_dt$sv_het_del_count > 0, na.rm = TRUE)
  message("[PRECOMP] SV prior annotation complete:")
  message("  Windows with SV INV overlap: ", n_with_sv, " / ", nrow(inv_like_dt))
  message("  Windows with het-DEL: ", n_with_del, " / ", nrow(inv_like_dt))
  message("  Total anchor assignments: ", sum(inv_like_dt$sv_n_anchors, na.rm = TRUE))
} else {
  message("[PRECOMP] No --sv_prior_dir provided — SV annotation columns will be 0/NA")
}

fwrite(inv_like_dt, f_inv, sep = "\t")
message("[PRECOMP] → ", f_inv)

# =============================================================================
# PER-CHROMOSOME PRECOMPUTE (parallelized with mclapply)
# =============================================================================

N_CORES <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), length(chroms))
if (N_CORES > 1) {
  message("[PRECOMP] Using ", N_CORES, " cores for parallel precompute")
}

precompute_one_chr <- function(chr) {
  t_chr <- proc.time()
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) return(NULL)

  dt <- as.data.table(chr_obj$out_dt)
  dt <- dt[order(start_bp)]
  dmat <- chr_obj$dmat
  mds_cols <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  mds_mat <- as.matrix(dt[, ..mds_cols])

  # Merge ALL inv_likeness + three-track + het metrics into precomp dt
  # This way the RDS has everything — downstream scripts don't need the TSV
  if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
    # Avoid duplicate columns that already exist in dt
    il_cols <- setdiff(names(inv_like_dt), c(names(dt), "chrom"))
    il_cols <- c("global_window_id", il_cols)
    dt <- merge(dt, inv_like_dt[, ..il_cols, with = FALSE],
                by = "global_window_id", all.x = TRUE)
    dt <- dt[order(start_bp)]
    message("[PRECOMP] ", chr, ": merged ", length(il_cols) - 1, " inv_like columns into dt")
  }

  # Robust z-scores (median/MAD)
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

  # Align dimensions
  n_dt <- nrow(dt); n_dm <- nrow(dmat)
  if (n_dm != n_dt) {
    n <- min(n_dm, n_dt)
    dt <- dt[seq_len(n)]; dmat <- dmat[seq_len(n), seq_len(n), drop = FALSE]
    mds_mat <- mds_mat[seq_len(n), , drop = FALSE]
  }

  # Similarity matrix
  sim_mat <- make_sim_mat(dmat)

  # ── BACKGROUND CONTINUITY BASELINE ──────────────────────────────────
  # Compute continuity scores between adjacent window pairs in background
  # territory (both windows below median inv_likeness). This gives the
  # "null" distribution — what continuity looks like when there's no
  # inversion. The extension accept threshold should be calibrated relative
  # to this baseline, not set to an arbitrary fixed number.
  #
  # Saved as quantiles so C01b can set:
  #   S1S accept = bg_q95 (very selective — top 5% of background)
  #   S1M accept = bg_q90 (selective — top 10%)
  #   S1L accept = bg_q85 (moderate — top 15%)

  bg_continuity <- numeric(0)
  chr_inv_med <- median(dt$inv_likeness[is.finite(dt$inv_likeness)])
  n_w <- nrow(dt)

  if (n_w >= 10 && is.finite(chr_inv_med)) {
    # Adjacent pairs where both are below-median inv_likeness = background
    bg_pairs <- which(
      dt$inv_likeness[seq_len(n_w - 1)] < chr_inv_med &
      dt$inv_likeness[seq(2, n_w)] < chr_inv_med
    )

    if (length(bg_pairs) >= 20) {
      # Sample up to 2000 pairs for speed
      if (length(bg_pairs) > 2000) bg_pairs <- sort(sample(bg_pairs, 2000))

      bg_continuity <- vapply(bg_pairs, function(i) {
        j <- i + 1L
        s_prev <- if (is.finite(sim_mat[i, j])) sim_mat[i, j] else 0
        # Simple: just use sim_mat[i,j] as proxy for continuity
        # (full compute_continuity needs roll buffer which doesn't exist yet)
        s_prev
      }, numeric(1))
    }
  }

  # Compute quantiles
  if (length(bg_continuity) >= 10) {
    bg_q <- quantile(bg_continuity, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  } else {
    # Fallback: use global sim_mat adjacent pairs
    adj_sims <- vapply(seq_len(n_w - 1), function(i) sim_mat[i, i + 1L], numeric(1))
    bg_q <- quantile(adj_sims, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  }

  message("[PRECOMP] ", chr, ": bg_continuity baseline: ",
          "q50=", round(bg_q["50%"], 3),
          " q85=", round(bg_q["85%"], 3),
          " q90=", round(bg_q["90%"], 3),
          " q95=", round(bg_q["95%"], 3),
          " (from ", length(bg_continuity), " bg pairs)")

  # Seed NN distances
  nn_dists <- vapply(seq_len(nrow(dmat)), function(i) {
    d <- dmat[i, ]; d[i] <- Inf; d <- d[is.finite(d)]
    if (length(d) == 0) return(Inf)
    mean(sort(d)[seq_len(min(SEED_NEIGHBOR_K, length(d)))], na.rm = TRUE)
  }, numeric(1))
  dmax_nn <- quantile(dmat[is.finite(dmat)], 0.95, na.rm = TRUE)
  if (!is.finite(dmax_nn) || dmax_nn == 0) dmax_nn <- 1
  dt[, seed_nn_dist := nn_dists / dmax_nn]

  # ═══════════════════════════════════════════════════════════════════════
  # MORPHOLOGY FEATURES (v9)
  # ═══════════════════════════════════════════════════════════════════════
  # Compute per-window features that describe the LOCAL MORPHOLOGY of
  # the z-score landscape and similarity matrix. These features let
  # downstream code distinguish:
  #   A. Broad flat coherent inversion-like plateau
  #   B. Compact spiky localized inversion-like peak
  #   C. Fragmented messy family/haplotype structure
  #   D. Isolated noise / background
  #
  # All features are written directly into dt so the .precomp.rds
  # is the one-stop cache for everything.
  # ═══════════════════════════════════════════════════════════════════════

  n_w <- nrow(dt)
  z_vec <- dt$max_abs_z  # primary z-score track

  # Fallback-safe: if z_vec is all zero or missing, skip morphology
  if (n_w >= 5 && any(is.finite(z_vec) & z_vec > 0)) {

    z_vec[!is.finite(z_vec)] <- 0

    # ── 1. Z-PROFILE MORPHOLOGY ──────────────────────────────────────

    # 1a. Local jaggedness: mean |z[i] - z[i-1]| in a ±5 neighborhood
    #     High = unstable / fragmented; Low = smooth / coherent
    JAGGED_HALF <- 5L
    jagged_vec <- rep(NA_real_, n_w)
    for (i in seq_len(n_w)) {
      lo <- max(2L, i - JAGGED_HALF)
      hi <- min(n_w, i + JAGGED_HALF)
      local_z <- z_vec[lo:hi]
      if (length(local_z) >= 3) {
        jagged_vec[i] <- mean(abs(diff(local_z)), na.rm = TRUE)
      }
    }
    dt[, local_jaggedness := round(jagged_vec, 4)]

    # 1b. Elevated run detection
    #     For each window, find the contiguous run of above-threshold
    #     windows it belongs to. Use a chromosome-adaptive threshold:
    #     75th percentile of z_vec (i.e. the top quartile is "elevated")
    z_thresh <- max(quantile(z_vec, 0.75, na.rm = TRUE), 1.5)
    is_elevated <- z_vec >= z_thresh

    # Label connected components of elevated windows
    run_id <- integer(n_w)
    current_run <- 0L
    for (i in seq_len(n_w)) {
      if (is_elevated[i]) {
        if (i == 1L || !is_elevated[i - 1L]) current_run <- current_run + 1L
        run_id[i] <- current_run
      }
    }

    # Per-run stats
    run_len_vec   <- rep(0L, n_w)
    run_mean_z    <- rep(NA_real_, n_w)
    run_sd_z      <- rep(NA_real_, n_w)
    run_start_vec <- rep(NA_integer_, n_w)
    run_end_vec   <- rep(NA_integer_, n_w)

    if (current_run > 0) {
      for (rid in seq_len(current_run)) {
        idx <- which(run_id == rid)
        rlen <- length(idx)
        rmean <- mean(z_vec[idx], na.rm = TRUE)
        rsd   <- if (rlen >= 2) sd(z_vec[idx], na.rm = TRUE) else 0
        rstart <- idx[1]; rend <- idx[rlen]
        run_len_vec[idx]   <- rlen
        run_mean_z[idx]    <- rmean
        run_sd_z[idx]      <- rsd
        run_start_vec[idx] <- rstart
        run_end_vec[idx]   <- rend
      }
    }

    dt[, local_run_len     := run_len_vec]
    dt[, local_run_mean_z  := round(run_mean_z, 4)]
    dt[, local_run_sd_z    := round(run_sd_z, 4)]
    dt[, run_start_idx     := run_start_vec]
    dt[, run_end_idx       := run_end_vec]

    # Candidate length in bp (from run boundaries)
    dt[, candidate_length_windows := local_run_len]
    dt[, candidate_length_bp := fifelse(
      is.finite(run_start_idx) & is.finite(run_end_idx),
      as.numeric(dt$end_bp[pmin(run_end_vec, n_w)] - dt$start_bp[pmax(run_start_vec, 1L)]),
      NA_real_
    )]

    # 1c. Peak prominence: how far this window's z exceeds its local
    #     background (±20 windows, using 25th percentile as background)
    PROM_HALF <- 20L
    prom_vec <- rep(NA_real_, n_w)
    for (i in seq_len(n_w)) {
      lo <- max(1L, i - PROM_HALF)
      hi <- min(n_w, i + PROM_HALF)
      bg_z <- z_vec[lo:hi]
      bg_level <- quantile(bg_z, 0.25, na.rm = TRUE)
      prom_vec[i] <- z_vec[i] - bg_level
    }
    dt[, local_peak_prominence := round(prom_vec, 4)]

    # 1d. Plateau flatness: combines run length, low SD, and low jaggedness
    #     High = broad coherent plateau; Low = spiky or no run
    #     Formula: run_len * (1 / (1 + run_sd)) * (1 / (1 + jaggedness))
    pf_vec <- rep(0, n_w)
    for (i in seq_len(n_w)) {
      rl <- run_len_vec[i]
      if (rl >= 3) {
        sd_term  <- 1 / (1 + (run_sd_z[i] %||% 0))
        jag_term <- 1 / (1 + (jagged_vec[i] %||% 0))
        pf_vec[i] <- rl * sd_term * jag_term
      }
    }
    # Normalize to [0, 1] within chromosome
    pf_max <- max(pf_vec, na.rm = TRUE)
    if (is.finite(pf_max) && pf_max > 0) pf_vec <- pf_vec / pf_max
    dt[, plateau_flatness := round(pf_vec, 4)]

    # ── 2. NEIGHBORHOOD SUPPORT ──────────────────────────────────────
    # Fraction of neighbors above z_thresh at different scales
    NBHOOD_RADII <- c(5L, 10L, 20L)
    for (rad in NBHOOD_RADII) {
      sup_vec <- rep(NA_real_, n_w)
      for (i in seq_len(n_w)) {
        lo <- max(1L, i - rad)
        hi <- min(n_w, i + rad)
        neighbors <- is_elevated[lo:hi]
        sup_vec[i] <- mean(neighbors, na.rm = TRUE)
      }
      col_name <- paste0("nbhood_support_", rad)
      dt[[col_name]] <- round(sup_vec, 4)
    }

    # ── 3. SIMILARITY MATRIX MORPHOLOGY ──────────────────────────────
    # For each window, extract a local neighborhood from sim_mat and
    # compute structural descriptors.

    SIM_HALF <- 10L  # ±10 windows = 21x21 local block

    block_compact  <- rep(NA_real_, n_w)
    block_coherence <- rep(NA_real_, n_w)
    block_frag     <- rep(NA_real_, n_w)
    square_support <- rep(NA_real_, n_w)

    # Chromosome-wide median sim (computed once for square_support)
    chr_sim_med <- median(sim_mat[upper.tri(sim_mat)], na.rm = TRUE)

    for (i in seq_len(n_w)) {
      lo <- max(1L, i - SIM_HALF)
      hi <- min(n_w, i + SIM_HALF)
      bsz <- hi - lo + 1L
      if (bsz < 5) next

      block <- sim_mat[lo:hi, lo:hi, drop = FALSE]

      # 3a. Block compactness: mean of the block (high = dense internal similarity)
      block_vals <- block[upper.tri(block)]
      if (length(block_vals) >= 3) {
        block_compact[i] <- mean(block_vals, na.rm = TRUE)
      }

      # 3b. Block coherence: ratio of within-diagonal-band to off-diagonal
      #     A true block has high sim near the diagonal AND high sim far
      #     from diagonal (square-like). A stripe has high near-diag only.
      diag_dist <- abs(row(block) - col(block))
      near_diag <- block[diag_dist <= 2 & diag_dist > 0]
      far_diag  <- block[diag_dist > max(2, bsz %/% 3)]
      if (length(near_diag) >= 2 && length(far_diag) >= 2) {
        mn_near <- mean(near_diag, na.rm = TRUE)
        mn_far  <- mean(far_diag, na.rm = TRUE)
        # Coherence: if far is close to near, it's a solid block
        block_coherence[i] <- if (mn_near > 0.01) mn_far / mn_near else 0
      }

      # 3c. Block fragmentation: SD of block values / mean
      #     High CV = heterogeneous internal structure = fragmented
      #     Low CV = uniform block = coherent
      if (length(block_vals) >= 3) {
        bm <- mean(block_vals, na.rm = TRUE)
        bsd <- sd(block_vals, na.rm = TRUE)
        block_frag[i] <- if (bm > 0.01) bsd / bm else 0
      }

      # 3d. Square support: fraction of the block above chromosome-wide
      #     median similarity. A real block has more high-sim values.
      if (is.finite(chr_sim_med) && length(block_vals) >= 3) {
        square_support[i] <- mean(block_vals > chr_sim_med, na.rm = TRUE)
      }
    }

    dt[, local_block_compactness  := round(block_compact, 4)]
    dt[, local_block_coherence    := round(block_coherence, 4)]
    dt[, local_block_fragmentation := round(block_frag, 4)]
    dt[, local_square_support     := round(square_support, 4)]

    # ── 4. COMPOSITE MORPHOLOGY SCORES ───────────────────────────────
    #
    # These are interpretable heuristic descriptors, not confidence values.
    # Each rewards the features described in the design spec.

    # Helper: safe clamp to [0,1], NA → 0
    clamp01 <- function(x) { x[!is.finite(x)] <- 0; pmin(1, pmax(0, x)) }

    # Helper: replace NA with 0 in a vector
    na0 <- function(x) { x[!is.finite(x)] <- 0; x }

    # 4a. flat_inv_score
    #     Rewards: long run, moderate/high z, low SD, low jaggedness,
    #     high block compactness, high coherence, low fragmentation.
    #     NOTE: family_likeness is NOT used here. In a hatchery, real
    #     founder-linked inversions show high family_likeness. Using it
    #     as a penalty would kill rare real inversions.
    fl_run   <- clamp01(log2(pmax(run_len_vec, 1)) / 5)  # log2(32)=5 → full credit at 32 wins
    fl_z     <- clamp01(na0(run_mean_z) / 4)
    fl_stab  <- clamp01(1 - na0(run_sd_z) / 2)
    fl_jag   <- clamp01(1 - na0(jagged_vec) / 2)
    fl_block <- clamp01(na0(block_compact))
    fl_coher <- clamp01(na0(block_coherence))
    fl_nofrag <- clamp01(1 - na0(block_frag))

    flat_raw <- (0.22 * fl_run + 0.16 * fl_z + 0.12 * fl_stab +
                 0.10 * fl_jag + 0.16 * fl_block + 0.12 * fl_coher +
                 0.12 * fl_nofrag)
    # Zero out windows not in any elevated run
    flat_raw[run_len_vec < 3] <- 0
    dt[, flat_inv_score := round(flat_raw, 4)]

    # 4b. spiky_inv_score
    #     Rewards: local peak prominence, some short-range support,
    #     compact local block, moderate/low fragmentation.
    #     Does NOT require long run. family_likeness NOT used (see 4a note).
    sp_prom  <- clamp01(na0(prom_vec) / 3)
    sp_sup5  <- clamp01(na0(dt$nbhood_support_5))
    sp_block <- clamp01(na0(block_compact))
    sp_nofrag <- clamp01(1 - na0(block_frag) * 0.5)
    # Bonus: short run (3-10 windows) gets full credit; long run demoted
    sp_short <- clamp01(fifelse(run_len_vec >= 3 & run_len_vec <= 12, 1,
                        fifelse(run_len_vec > 12, 0.5, 0.3)))
    # Anti-feature: if the run is very long AND flat, this is a flat candidate, not spiky
    sp_antiflat <- clamp01(1 - na0(dt$plateau_flatness) * 0.5)

    spiky_raw <- (0.28 * sp_prom + 0.17 * sp_sup5 + 0.17 * sp_block +
                  0.12 * sp_nofrag + 0.12 * sp_short + 0.14 * sp_antiflat)
    # Require minimum z to be nonzero
    spiky_raw[z_vec < z_thresh * 0.5] <- 0
    dt[, spiky_inv_score := round(spiky_raw, 4)]

    # 4c. fragmentation_score
    #     Rewards: high jaggedness, high block fragmentation, low coherence,
    #     unstable neighborhood, many short alternating disruptions
    fg_jag   <- clamp01(na0(jagged_vec) / 2)
    fg_frag  <- clamp01(na0(block_frag))
    fg_nocoher <- clamp01(1 - na0(block_coherence))
    # Neighborhood instability: support at ±5 but NOT at ±20
    # → windows where local structure doesn't persist
    sup5  <- na0(dt$nbhood_support_5)
    sup20 <- na0(dt$nbhood_support_20)
    fg_instab <- clamp01(pmax(0, sup5 - sup20))
    # Short alternating runs: if many very short runs exist nearby
    fg_shortrun <- clamp01(fifelse(run_len_vec >= 1 & run_len_vec <= 2, 1,
                           fifelse(run_len_vec == 3, 0.5, 0)))

    frag_raw <- (0.25 * fg_jag + 0.20 * fg_frag + 0.15 * fg_nocoher +
                 0.20 * fg_instab + 0.20 * fg_shortrun)
    # Only compute for windows that have SOME structure
    frag_raw[z_vec < 1.0] <- 0
    dt[, fragmentation_score := round(frag_raw, 4)]

    message("[PRECOMP] ", chr, ": morphology features computed (",
            "flat>0.3: ", sum(dt$flat_inv_score > 0.3, na.rm = TRUE),
            ", spiky>0.3: ", sum(dt$spiky_inv_score > 0.3, na.rm = TRUE),
            ", frag>0.3: ", sum(dt$fragmentation_score > 0.3, na.rm = TRUE), ")")

  } else {
    # Chromosome too small or no z-score variation — fill defaults
    morph_cols <- c("local_jaggedness", "local_run_len", "local_run_mean_z",
                    "local_run_sd_z", "run_start_idx", "run_end_idx",
                    "candidate_length_windows", "candidate_length_bp",
                    "local_peak_prominence", "plateau_flatness",
                    "nbhood_support_5", "nbhood_support_10", "nbhood_support_20",
                    "local_block_compactness", "local_block_coherence",
                    "local_block_fragmentation", "local_square_support",
                    "flat_inv_score", "spiky_inv_score", "fragmentation_score")
    for (mc in morph_cols) dt[[mc]] <- NA_real_
    message("[PRECOMP] ", chr, ": too few windows for morphology (", n_w, ")")
  }

  # ── PER-WINDOW PVE (v9.4) ──
  # Stamps PVE_1..PVE_<npc> per window so downstream consumers (scree plots,
  # scrubber overview page) can show "how much of the local PCA variance
  # is captured by each PC at each window" without recomputing.
  #
  # Definition: PVE_k = lam_k / sum(lam_1..lam_<npc>).
  # We use the sum of available lam_* columns as the denominator (not
  # `total`, which is sum of squared eigenvalues — different scale).
  # If only PC1 is present, PVE_1 = 1.0 by construction.
  lam_cols <- grep("^lam_[0-9]+$", names(dt), value = TRUE)
  if (length(lam_cols) >= 1) {
    lam_mat <- as.matrix(dt[, ..lam_cols])
    # Per-window denominator = sum of lambdas; avoid /0
    lam_sum <- rowSums(lam_mat, na.rm = TRUE)
    lam_sum[lam_sum <= 0 | !is.finite(lam_sum)] <- NA_real_
    for (lc in lam_cols) {
      pc_idx <- as.integer(sub("^lam_", "", lc))
      pve_col <- paste0("PVE_", pc_idx)
      dt[[pve_col]] <- round(dt[[lc]] / lam_sum, 6)
    }
    message("[PRECOMP] ", chr, ": stamped PVE for ", length(lam_cols),
            " PCs (", paste(lam_cols, collapse = ", "), ")")
  } else {
    message("[PRECOMP] ", chr, ": no lam_ columns found, skipping PVE stamping")
  }

  # ── MULTI-PC CHROMOSOME-WIDE BANDS (v9.4) ──
  # Compute k=3 bands for each available PC (PC2/3/4 in addition to PC1
  # which was already done inside compute_inv_likeness_all). PC1 bands are
  # the canonical signal that drives inv_likeness; higher PCs are exposed
  # for downstream sub-cluster detection, secondary-inversion candidate
  # detection, and scrubber multi-PC visualization. They do NOT alter
  # PC1 metrics or any per-window scoring.
  chr_bands_per_pc <- list()
  chr_band_centers_per_pc <- list()
  for (pc_idx in seq_len(npc_detected)) {
    pc_cols_n <- grep(paste0("^PC_", pc_idx, "_"), names(dt), value = TRUE)
    if (length(pc_cols_n) < 20) next
    pc_mat_n <- as.matrix(dt[, ..pc_cols_n])
    avg_pc_n <- colMeans(pc_mat_n, na.rm = TRUE)
    valid_n  <- is.finite(avg_pc_n)
    if (sum(valid_n) < 20) next
    avg_valid_n <- avg_pc_n[valid_n]
    chr_km_n <- tryCatch(kmeans(avg_valid_n, centers = 3, nstart = 10),
                         error = function(e) NULL)
    if (is.null(chr_km_n)) next
    co_n <- order(chr_km_n$centers[, 1])
    bands_n <- integer(length(avg_valid_n))
    bands_n[chr_km_n$cluster == co_n[1]] <- 1L
    bands_n[chr_km_n$cluster == co_n[2]] <- 2L
    bands_n[chr_km_n$cluster == co_n[3]] <- 3L
    names(bands_n) <- sub(paste0("^PC_", pc_idx, "_"), "", names(avg_valid_n))
    key_n <- paste0("pc", pc_idx)
    chr_bands_per_pc[[key_n]] <- bands_n
    chr_band_centers_per_pc[[key_n]] <- chr_km_n$centers[co_n, 1]
    if (pc_idx > 1L) {
      message("  [", chr, "] Fixed bands (PC", pc_idx, "): ",
              sum(bands_n == 1), "/", sum(bands_n == 2), "/", sum(bands_n == 3),
              " (propagated; not used in inv_likeness)")
    }
  }

  # Save precomputed data
  precomp <- list(
    dt = dt, sim_mat = sim_mat, mds_mat = mds_mat,
    chrom = chr, n_windows = nrow(dt),
    bg_continuity_quantiles = bg_q,
    # v9.4: multi-PC propagation. The dt itself already carries PC_1_*,
    # PC_2_*, PC_3_*, PC_4_* columns (whatever was in upstream out_dt).
    # The chr_bands_per_pc / chr_band_centers_per_pc fields expose the
    # k=3 chromosome-wide grouping per PC, which downstream consumers
    # (sub-cluster detection, scrubber) can use directly without
    # recomputing kmeans. PC1 remains the canonical band axis for
    # inv_likeness; higher PCs are advisory.
    npc = npc_detected,
    chr_bands_per_pc = chr_bands_per_pc,
    chr_band_centers_per_pc = chr_band_centers_per_pc,
    # v9.4: chromosome-level scree summary. pve_per_pc_mean[k] = mean of
    # PVE_k across all windows on this chromosome (variance-weighted PVE
    # estimate at chromosome scale). Used for scree plots in the scrubber.
    # pve_per_pc_median is the same with median instead of mean (more
    # robust to outlier windows like inversion regions).
    pve_per_pc_mean = if (length(grep("^PVE_", names(dt))) > 0) {
      pve_cols <- grep("^PVE_[0-9]+$", names(dt), value = TRUE)
      pve_idx  <- as.integer(sub("^PVE_", "", pve_cols))
      ord      <- order(pve_idx)
      setNames(sapply(pve_cols[ord], function(c) mean(dt[[c]], na.rm = TRUE)),
               paste0("PC", pve_idx[ord]))
    } else NULL,
    pve_per_pc_median = if (length(grep("^PVE_", names(dt))) > 0) {
      pve_cols <- grep("^PVE_[0-9]+$", names(dt), value = TRUE)
      pve_idx  <- as.integer(sub("^PVE_", "", pve_cols))
      ord      <- order(pve_idx)
      setNames(sapply(pve_cols[ord], function(c) median(dt[[c]], na.rm = TRUE)),
               paste0("PC", pve_idx[ord]))
    } else NULL
  )
  rds_out <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  saveRDS(precomp, rds_out)

  # BUGFIX 2026-04-17 (FIX 21, CRASH/DESIGN): port NN-smoothed sim_mat
  # production from archived snake1 precomp. Without this, the live
  # precomp only wrote the raw sim_mat (nn=0) inside <chr>.precomp.rds
  # and did NOT write per-scale sim_mat_nn<k>.rds files. Phase 2d
  # (`load_sim_mat` in run_all.R) looks for sim_mat_nn20/40/80/... and
  # was silently finding only nn=0, so:
  #   - Phase 4 (D02 NN persistence, needs >= 2 scales) was skipped
  #   - Phase 5 (D09 NN sweep tree, needs >= 3 scales) was skipped
  #   - C01d's iv$survives_nn40/nn80 were always NA
  #   - C01d's iv$nn_birth was always NA (tree never built)
  #   - The entire NN-persistence track in the scoring was dead
  # FIX 16 (tree merge into scoring table) was correct but couldn't
  # fire because tree was always NULL upstream.
  #
  # Semantics: NN smoothing is NOT genomic-adjacent window averaging.
  # It's MDS-space k-nearest-neighbor averaging: for each window, find
  # its k most-similar (smallest dmat row) windows anywhere on the
  # chromosome, average their MDS coordinates, recompute distances,
  # build a new sim_mat. Windows inside the same structural block
  # (inversion/family/background) have MDS-space neighbors sharing
  # that block; smoothing preserves block-level structure while
  # dissolving noise. Larger k = more smoothing = smaller blocks
  # dissolve into the ambient structure.
  #
  # WHY THE TREE (D09) NEEDS MULTIPLE SCALES:
  # The tree is a persistence-barcode over NN scale. At the coarsest
  # NN, only the biggest-and-most-persistent blocks survive — those
  # become roots. As NN decreases, blocks either stay stable (small
  # boundary refinements), split into children (nested inversions
  # resolving), disappear (oversmoothing artifact — block only looked
  # real at one scale), or novel blocks pop in (scale-specific real
  # small inversion, or noise).
  #
  # nn_birth = coarsest NN scale at which the block first appears as
  # its own distinct unit. High nn_birth = block survives heavy
  # smoothing = strong persistent structure = real large inversion.
  # Low nn_birth = block is scale-specific = noise or small local
  # signal. This is topological-data-analysis persistence for
  # structural blocks along a chromosome.
  #
  # Scale ladder (expanded from 20,40,80 per Quentin 2026-04-17 for
  # 226 samples / 50kb windows / median chr ~900 windows):
  #   20, 40, 80   = fine (catches ~100kb-1Mb inversions, nn_birth=20-80)
  #   120, 160     = middle (catches 1-5Mb inversions, nn_birth=120-160)
  #   200, 240     = coarse (catches 5-10Mb inversions, nn_birth>=200
  #                  which is C01d's D3 saturation threshold)
  #   320          = ceiling (confirms largest blocks persist; anything
  #                  nn_birth=320 is the strongest class)
  # Top scales auto-clamp to (n_windows - 1) for small chromosomes
  # via min(k, n_mds - 1L) at runtime, so nothing breaks.
  #
  # Going beyond 320 is not scientifically useful for this cohort
  # (largest suspected inversion ~10Mb = NN~200; beyond that only
  # chromosome-scale family-LD lingers, which isn't what D09 is
  # classifying). Diagnostic scales (480/640/800) could be added if
  # needed to confirm the noise floor, at the cost of ~2x precomp time.
  sim_dir <- file.path(precomp_dir, "sim_mats")
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

  # Always save nn0 as a separate file for phase 2d to load
  saveRDS(sim_mat, file.path(sim_dir, paste0(chr, ".sim_mat_nn0.rds")))
  message("[PRECOMP] ", chr, ": saved sim_mat_nn0")

  # Scale list from env. Default: scientifically-calibrated ladder that
  # reaches D09's classifier threshold (nn_birth>=200) without wasting
  # precomp time on oversmoothed scales where nothing survives.
  nn_sim_scales <- as.integer(strsplit(
    Sys.getenv("NN_SIM_SCALES",
               "20,40,80,120,160,200,240,320"), ",")[[1]])
  nn_sim_scales <- nn_sim_scales[nn_sim_scales > 0 & is.finite(nn_sim_scales)]

  mds_cols_nn <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  if (length(mds_cols_nn) > 0 && length(nn_sim_scales) > 0) {
    mds_mat_nn <- as.matrix(dt[, ..mds_cols_nn])
    n_mds <- nrow(mds_mat_nn)
    for (k in nn_sim_scales) {
      k_use <- min(k, n_mds - 1L)
      if (k_use < 2L) {
        message("[PRECOMP] ", chr, ": skipping nn", k,
                " (n_windows=", n_mds, " too small)")
        next
      }
      t_nn <- proc.time()
      smoothed <- matrix(0, nrow = n_mds, ncol = ncol(mds_mat_nn))
      for (wi in seq_len(n_mds)) {
        d <- dmat[wi, ]; d[wi] <- Inf
        nn_idx <- order(d)[seq_len(k_use)]
        smoothed[wi, ] <- colMeans(
          mds_mat_nn[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
      }
      nn_dmat <- as.matrix(dist(smoothed))
      nn_sim <- make_sim_mat(nn_dmat)
      saveRDS(nn_sim,
              file.path(sim_dir, paste0(chr, ".sim_mat_nn", k, ".rds")))
      elapsed_nn <- round((proc.time() - t_nn)[3], 1)
      message("[PRECOMP] ", chr, ": saved sim_mat_nn", k,
              " (k_use=", k_use, ", ", elapsed_nn, "s)")
    }
  } else {
    message("[PRECOMP] ", chr,
            ": skipping NN-smoothed sim_mats — no MDS columns or empty scale list")
  }

  elapsed <- round((proc.time() - t_chr)[3], 1)
  message("[PRECOMP] ", chr, ": ", nrow(dt), " windows (", elapsed, "s)")

  data.table(
    chrom = chr, n_windows = nrow(dt),
    n_inv_like_050 = sum(dt$inv_likeness >= 0.5, na.rm = TRUE),
    n_z_above_2 = sum(dt$max_abs_z >= 2.0, na.rm = TRUE),
    n_z_above_3 = sum(dt$max_abs_z >= 3.0, na.rm = TRUE),
    median_max_z = round(median(dt$max_abs_z, na.rm = TRUE), 3),
    q95_max_z = round(quantile(dt$max_abs_z, 0.95, na.rm = TRUE), 3),
    bg_q50 = round(bg_q["50%"], 4),
    bg_q85 = round(bg_q["85%"], 4),
    bg_q90 = round(bg_q["90%"], 4),
    bg_q95 = round(bg_q["95%"], 4),
    n_bg_pairs = length(bg_continuity),
    n_flat_030 = sum(dt$flat_inv_score > 0.3, na.rm = TRUE),
    n_spiky_030 = sum(dt$spiky_inv_score > 0.3, na.rm = TRUE),
    n_frag_030 = sum(dt$fragmentation_score > 0.3, na.rm = TRUE),
    elapsed_sec = elapsed
  )
}

message("[PRECOMP] Processing ", length(chroms), " chromosomes...")
t_all <- proc.time()

if (N_CORES > 1) {
  summary_list <- parallel::mclapply(chroms, precompute_one_chr, mc.cores = N_CORES)
} else {
  summary_list <- lapply(chroms, precompute_one_chr)
}

summary_list <- summary_list[!vapply(summary_list, is.null, logical(1))]
summary_dt <- rbindlist(summary_list)
elapsed_total <- round((proc.time() - t_all)[3], 1)
fwrite(summary_dt, file.path(outdir, "precomp_summary.tsv"), sep = "\t")

message("\n[DONE] Precompute complete (v9 — morphology features)")
message("  Chromosomes: ", length(chroms), " (", N_CORES, " cores, ", elapsed_total, "s total)")
message("  Total windows: ", sum(summary_dt$n_windows))
message("  Windows with inv_likeness >= 0.5: ", sum(summary_dt$n_inv_like_050))
message("  Windows with robust |z| >= 2.0: ", sum(summary_dt$n_z_above_2))
message("  Windows with robust |z| >= 3.0: ", sum(summary_dt$n_z_above_3))
message("  ── Morphology (v9) ──")
message("    flat_inv_score > 0.3:      ", sum(summary_dt$n_flat_030))
message("    spiky_inv_score > 0.3:     ", sum(summary_dt$n_spiky_030))
message("    fragmentation_score > 0.3: ", sum(summary_dt$n_frag_030))
message("  Precomp dir: ", precomp_dir)
message("  Inv-likeness: ", f_inv)
message("  Summary: ", file.path(outdir, "precomp_summary.tsv"))
message("\n  Next: run STEP_C01b seeded extension with desired parameters")
