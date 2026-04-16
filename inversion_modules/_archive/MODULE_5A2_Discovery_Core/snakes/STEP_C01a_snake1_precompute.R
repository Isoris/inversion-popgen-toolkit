#!/usr/bin/env Rscript
# =============================================================================
# STEP_C01a_snake1_precompute.R  (v9.1 -- merged)
#
# PRECOMPUTE STEP -- run once, reuse many times.
#
# WIRING (from v9.0):
#   - load_bridge.R: smap (Ind<->CGA), reg (sample registry), instant Q
#   - Ind->CGA rename at load time via smap
#   - optparse CLI: --chr, --bridge, --no-ancestry, --dosage-dir, --mode
#
# COMPUTATION (from v8.5.2, all kept):
#   - Robust z-scores (median/MAD)
#   - inv_likeness scores (trimodality gated on het, 3-track decomposition)
#   - Dosage het rates (genuinely independent from PCA)
#   - Distance matrix (dmat) and similarity matrix (sim_mat)
#   - Background continuity baselines
#   - Seed NN distances
#
# STORAGE (v8.5.3):
#   - sim_mat saved as SEPARATE RDS files (not inside precomp)
#   - NN-smoothed sim_mats saved at configurable scales (NN_SIM_SCALES env)
#   - Main precomp RDS has $sim_mat = NULL (saves ~640 MB per chr)
#
# =====================================================================
# PRECOMP RDS SCHEMA (v9.1)
# =====================================================================
#   $chrom, $n_windows, $schema_version = "v9.1"
#   $dt        : data.table -- per-window metrics (all columns, CGA names)
#   $sim_mat   : NULL -- REMOVED, use sim_mats/ directory
#   $mds_mat   : matrix [n_win x n_dims]
#   $dmat      : matrix [n_win x n_win] raw distances
#   $bg_continuity_quantiles : named numeric
#   $chr_bands : named integer (CGA names)
#   $nn_normalization : list
#   $sample_names     : character vector (CGA names)
#   $nn_sim_scales    : integer vector
#   $sim_mat_dir      : character
#   $name_format      : "CGA"
#   $has_instant_q    : logical
#
# SEPARATE SIM_MAT FILES (in precomp/sim_mats/):
#   <chr>.sim_mat_nn0.rds, <chr>.sim_mat_nn20.rds, etc.
# =====================================================================
#
# Usage:
#   Rscript STEP_C01a_snake1_precompute.R <step10_outprefix> <outdir> \
#     [--bridge <load_bridge.R>] [--no-ancestry] \
#     [--qmatrix <qopt>] [--relatedness <pairs>] [--dosage_dir <dir>] \
#     [--mode hatchery|wild]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01a_snake1_precompute.R <step10_outprefix> <outdir> [opts]")
}

step10_prefix <- args[1]
outdir        <- args[2]
precomp_dir   <- file.path(outdir, "precomp")
dir.create(precomp_dir, recursive = TRUE, showWarnings = FALSE)

# Parse optional args
qmatrix_file <- NULL; relatedness_file <- NULL; pop_mode <- "hatchery"
dosage_dir <- NULL; bridge_path <- NULL; no_ancestry <- FALSE
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--qmatrix" && i < length(args))      { qmatrix_file <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args)) { relatedness_file <- args[i+1]; i <- i+2 }
  else if (a == "--dosage_dir" && i < length(args)) { dosage_dir <- args[i+1]; i <- i+2 }
  else if (a == "--mode" && i < length(args))    { pop_mode <- args[i+1]; i <- i+2 }
  else if (a == "--bridge" && i < length(args))  { bridge_path <- args[i+1]; i <- i+2 }
  else if (a == "--no-ancestry")                 { no_ancestry <- TRUE; i <- i+1 }
  else { i <- i+1 }
}

# =============================================================================
# WIRE: load_bridge.R (provides smap, reg, instant Q)
# =============================================================================

# Find load_bridge.R
if (is.null(bridge_path) || !file.exists(bridge_path)) {
  bridge_path <- Sys.getenv("LOAD_BRIDGE", "")
  if (!nzchar(bridge_path) || !file.exists(bridge_path)) {
    for (p in c("utils/load_bridge.R",
                "../utils/load_bridge.R",
                file.path(Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"),
                          "inversion_codebase_v8.5/utils/load_bridge.R"))) {
      if (file.exists(p)) { bridge_path <- p; break }
    }
  }
}

has_bridge <- FALSE
if (!is.null(bridge_path) && file.exists(bridge_path)) {
  tryCatch({
    source(bridge_path)
    has_bridge <- TRUE
    message("[PRECOMP] load_bridge.R sourced: ", bridge_path)
    if (exists("smap") && !is.null(smap)) {
      message("[PRECOMP] smap: ", smap$n, " samples (", smap$detect(smap$ind_names[1]), " -> ", smap$real_names[1], ")")
    }
  }, error = function(e) {
    message("[PRECOMP] WARNING: load_bridge.R failed: ", e$message)
    message("[PRECOMP] Continuing without bridge wiring")
  })
} else {
  message("[PRECOMP] load_bridge.R not found (tried: ", bridge_path %||% "none", ")")
  message("[PRECOMP] Continuing without bridge wiring (old-style sample names)")
}

# =============================================================================
# PARAMETERS (same as C01 -- must match)
# =============================================================================

SEED_MDS_AXES    <- 5L
SEED_NEIGHBOR_K  <- 3L

DOSAGE_HET_LO <- 0.3
DOSAGE_HET_HI <- 1.7

# =============================================================================
# LOAD MDS DATA
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

# Extract sample names and rename Ind->CGA if bridge available
sample_names <- NULL
for (chr_tmp in chroms) {
  dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
  pc1_cols <- grep("^PC_1_", names(dt_tmp), value = TRUE)
  if (length(pc1_cols) > 0) {
    sample_names <- sub("^PC_1_", "", pc1_cols)
    break
  }
}

# === THE ROOT FIX: Rename Ind->CGA in ALL per_chr data ===
name_format <- "unknown"
if (!is.null(sample_names) && has_bridge && exists("smap") && !is.null(smap)) {
  if (smap$detect(sample_names[1]) == "ind") {
    message("[PRECOMP] Renaming Ind->CGA across all chromosomes...")
    real_names <- smap$to_real_vec(sample_names)

    for (chr_tmp in chroms) {
      dt_tmp <- as.data.table(per_chr[[chr_tmp]]$out_dt)
      for (pfx in c("PC_1_", "PC_2_")) {
        old_cols <- paste0(pfx, sample_names)
        new_cols <- paste0(pfx, real_names)
        present <- old_cols %in% names(dt_tmp)
        if (any(present)) {
          setnames(dt_tmp, old_cols[present], new_cols[present])
        }
      }
      per_chr[[chr_tmp]]$out_dt <- dt_tmp
    }
    sample_names <- real_names
    name_format <- "CGA"
    message("[PRECOMP] Renamed: ", sample_names[1], " ... ", sample_names[length(sample_names)])
  } else {
    name_format <- "CGA"
    message("[PRECOMP] Columns already use CGA names")
  }
} else {
  name_format <- if (!is.null(sample_names) && grepl("^CGA", sample_names[1])) "CGA" else "Ind"
}

if (!is.null(sample_names)) {
  message("[PRECOMP] Sample names: ", length(sample_names), " (format: ", name_format, ")")
}

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
    if (nrow(q_raw) == length(sample_names)) {
      rownames(q_raw) <- sample_names
      q_mat <- q_raw
      K_ancestry <- ncol(q_mat)
      message("[PRECOMP] Q matrix loaded: ", nrow(q_mat), " samples x K=", K_ancestry)
    } else {
      message("[PRECOMP] Q matrix row count (", nrow(q_raw),
              ") != sample count (", length(sample_names), ") -- skipping")
    }
  } else if (!is.null(qmatrix_file)) {
    message("[PRECOMP] Q matrix file not readable: ", qmatrix_file)
  }
} else if (!is.null(qmatrix_file)) {
  message("[PRECOMP] Q matrix file not found: ", qmatrix_file, " -- family track disabled")
}

# Try instant Q from bridge (v9.0 wiring) if no qmatrix file provided
has_instant_q <- FALSE
if (!no_ancestry && has_bridge && exists("get_Q_summary", mode = "function")) {
  message("[PRECOMP] Instant Q available via load_bridge")
  has_instant_q <- TRUE
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

      # ISSUE_010: using parameterized thresholds
      per_sample_het <- colMeans(w_dos >= DOSAGE_HET_LO & w_dos <= DOSAGE_HET_HI, na.rm = TRUE)
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

    # Compute chromosome-wide baselines for relative scoring
    chr_pve1_med <- NA_real_; chr_eig_ratio_med <- NA_real_
    if (has_lam) {
      l1_all <- dt$lam_1[is.finite(dt$lam_1)]
      l2_all <- dt$lam_2[is.finite(dt$lam_2)]
      if (length(l1_all) > 20 && length(l2_all) > 20) {
        pve1_all <- l1_all / (l1_all + l2_all)
        chr_pve1_med <- median(pve1_all, na.rm = TRUE)
        chr_eig_ratio_med <- median(l1_all / pmax(l2_all, 1e-10), na.rm = TRUE)
      }
    }

    # ── CHROMOSOME-WIDE FIXED k=3 BANDS (v8.5) ──
    # Computed ONCE per chromosome. Used as fixed group labels for all
    # per-window entropy/delta12/het metrics. Avoids per-window k-means
    # (which is unstable, causes label switching, produces noisy tracks).
    #
    # Method: average PC1 across ALL windows per sample, then k=3 on
    # those averages. This gives the dominant 3-group structure of the
    # chromosome. Individual windows measure conformance to this fixed grouping.
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
      pve1 <- eigen_ratio <- dip_p <- NA_real_
      lambda_ratio <- NA_real_; het_contrast <- NA_real_

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

      if (has_lam) {
        l1 <- dt$lam_1[wi]; l2 <- dt$lam_2[wi]
        if (is.finite(l1) && is.finite(l2) && (l1 + l2) > 0) {
          pve1 <- l1 / (l1 + l2)
          eigen_ratio <- l1 / max(l2, 1e-10)
          lambda_ratio <- l1 / max(l2, 0.01)
        }
      }

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
      # REVISED INV_LIKENESS (v8.5.2)
      # ═══════════════════════════════════════════════════════════════
      # ISSUE_014: Trimodality is a GATED term. Without het evidence,
      #   trimodality alone does not make a window inversion-like.
      #   Generic trimodal outliers (family structure) are thus penalized.
      #
      # ISSUE_016: This is the FIRST PASS (PCA-derived het_contrast).
      #   When --dosage_dir is provided, a SECOND PASS replaces het_contrast
      #   with dosage_het_rate_cv (Source C, genuinely independent).
      #   Both passes use the same formula structure.
      #
      # Formula: 50% HET + 25% STRUCTURE + 25% TRIMODALITY (gated)
      # ═══════════════════════════════════════════════════════════════

      # s_pve: PVE1 RELATIVE to chromosome median
      s_pve <- 0
      if (is.finite(pve1) && is.finite(chr_pve1_med)) {
        pve1_excess <- pve1 - chr_pve1_med
        s_pve <- pmin(1, pmax(0, (pve1_excess + 0.1) / 0.3))
      } else if (is.finite(pve1)) {
        s_pve <- pmin(1, pmax(0, (pve1 - 0.3) / 0.5))
      }

      # s_dip: trimodality (raw score, before gating)
      s_dip_raw <- if (is.finite(dip_p)) pmin(1, pmax(0, (0.5 - dip_p) / 0.45)) else 0

      # s_het: use het_contrast (PCA-derived) in first pass.
      # Dosage het CV upgrades this in the second pass after dosage merge.
      s_het <- 0
      het_source <- "none"
      if (is.finite(het_contrast)) {
        s_het <- pmin(1, pmax(0, (het_contrast - 2) / 8))
        het_source <- "pc1_het_contrast"
      }

      # ISSUE_014: GATE trimodality on het evidence.
      # Trimodality only fully contributes if het signal >= 0.1.
      # Without het, trimodality is heavily damped (x0.25) to prevent
      # family-structure outliers from scoring as inversion-like.
      s_dip <- if (s_het >= 0.1) s_dip_raw else s_dip_raw * 0.25

      # Combined: 50% het + 25% structure + 25% trimodality (gated)
      inv_like <- 0.50 * s_het + 0.25 * s_pve + 0.25 * s_dip

      # ═══════════════════════════════════════════════════════════════
      # THREE-TRACK DECOMPOSITION (v8.5.2)
      # ═══════════════════════════════════════════════════════════════
      # ISSUE_015: Track 1 captures GENERAL structural strength (PVE1).
      #   High PVE1 = PC1 dominates — could be inversion, family, demographic.
      #   Does NOT assess whether structure resolves into discrete banded groups.
      #   Trimodality (in Track 2 via s_dip) adds that distinction:
      #   continuous gradient = family LD; discrete 3-group = inversion.
      # Track 1: STRUCTURE LIKENESS — general structural strength
      structure_like <- if (is.finite(pve1)) pmin(1, pmax(0, (pve1 - 0.5) / 0.4)) else 0

      # Track 2: INV LIKENESS — existing composite (trimodality + het_contrast)
      #   Already computed as inv_like above.

      # Track 3: FAMILY LIKENESS — do k=3 bands correlate with Q ancestry?
      #   High = bands recapitulate family/ancestry groups → family-driven
      #   Low  = bands mix ancestry groups → structural signal (inversion)
      #
      # Method: for each k=3 band (from per-window km3), compute the mean
      # Q vector. Then measure how DIFFERENT the band-mean Q vectors are:
      #   - between-band Q variance (high → bands track ancestry → family)
      #   - within-band Q entropy (low → each band is one ancestry → family)
      # Score: between-band Q variance / total Q variance
      #   = fraction of Q variance explained by the k=3 bands
      #   = how well ancestry predicts band membership
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

      all_rows[[length(all_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = wid,
        start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
        inv_pve1 = round(pve1, 4), inv_eigen_ratio = round(eigen_ratio, 2),
        inv_dip_p = round(dip_p, 4),
        inv_het_contrast = round(het_contrast, 2),
        inv_pve1_excess = round(if (is.finite(pve1) && is.finite(chr_pve1_med))
          pve1 - chr_pve1_med else NA_real_, 4),
        inv_likeness = round(inv_like, 4),
        # ── THREE-TRACK DECOMPOSITION (v8.5) ──
        structure_likeness = round(structure_like, 4),
        family_likeness = round(family_like %||% NA_real_, 4),
        family_q_between = round(family_q_between %||% NA_real_, 6),
        family_q_within_entropy = round(family_q_within_entropy %||% NA_real_, 4),
        # Multi-source het assessment (v8.5)
        pc1_bimodality = round(pc1_bimodality %||% NA_real_, 4),
        het_pc1_gap = round(het_pc1_gap %||% NA_real_, 4),
        het_mid_fraction = round(het_mid_fraction %||% NA_real_, 4),
        het_mid_variance = round(het_mid_variance %||% NA_real_, 6),
        # Local structure metrics (v8.5)
        local_delta12 = round(local_delta12 %||% NA_real_, 4),
        local_entropy = round(local_entropy %||% NA_real_, 4),
        local_ena = round(local_ena %||% NA_real_, 4),
        # Band sizes for downstream stratification
        band1_frac = round(band_sizes[1] %||% NA_real_, 4),
        band2_frac = round(band_sizes[2] %||% NA_real_, 4),
        band3_frac = round(band_sizes[3] %||% NA_real_, 4),
        # Dosage-based het rate: genuinely independent from PCA
        # Computed from raw dosage values, not PC1 projection
        # To be filled by a second pass that reads dosage files,
        # or left NA if dosage_dir not available.
        dosage_het_rate_median = round(dosage_het_rate_median %||% NA_real_, 4),
        dosage_het_rate_sd = round(dosage_het_rate_sd %||% NA_real_, 4)
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

  # Merge into inv_like_dt
  if (nrow(dosage_het_dt) > 0) {
    # dosage_het_dt has: chrom, global_window_id, dosage_het_rate_median, _sd, _cv
    # inv_like_dt has placeholder NA columns with same names
    # Update in place
    setkey(dosage_het_dt, chrom, global_window_id)
    setkey(inv_like_dt, chrom, global_window_id)

    match_idx <- inv_like_dt[dosage_het_dt, which = TRUE, on = c("chrom", "global_window_id")]
    match_idx <- match_idx[!is.na(match_idx)]

    if (length(match_idx) > 0) {
      matched_dos <- dosage_het_dt[inv_like_dt[match_idx, .(chrom, global_window_id)],
                                    on = c("chrom", "global_window_id")]
      inv_like_dt[match_idx, dosage_het_rate_median := matched_dos$dosage_het_rate_median]
      inv_like_dt[match_idx, dosage_het_rate_sd := matched_dos$dosage_het_rate_sd]
      inv_like_dt[match_idx, dosage_het_rate_cv := matched_dos$dosage_het_rate_cv]
    }

    # ISSUE_016: SECOND PASS — dosage CV replaces het_contrast as primary.
    # This upgrades the het component from PCA-derived (Source A) to
    # dosage-derived (Source C). Formula structure identical to first pass.
    n_with_dosage <- sum(is.finite(inv_like_dt$dosage_het_rate_cv) &
                          inv_like_dt$dosage_het_rate_cv > 0)
    if (n_with_dosage > 0) {
      message("[PRECOMP] Upgrading inv_likeness for ", n_with_dosage,
              " windows with dosage het CV (second pass)...")
      inv_like_dt[is.finite(dosage_het_rate_cv) & dosage_het_rate_cv > 0,
        inv_likeness := {
          s_h <- pmin(1, pmax(0, (dosage_het_rate_cv - 0.1) / 0.6))
          s_p <- pmin(1, pmax(0, (inv_pve1_excess + 0.1) / 0.3))
          s_p <- fifelse(is.finite(s_p), s_p, 0)
          s_d_raw <- pmin(1, pmax(0, (0.5 - inv_dip_p) / 0.45))
          s_d_raw <- fifelse(is.finite(s_d_raw), s_d_raw, 0)
          # ISSUE_014: gate trimodality on het evidence in second pass too
          s_d <- fifelse(s_h >= 0.1, s_d_raw, s_d_raw * 0.25)
          round(0.50 * s_h + 0.25 * s_p + 0.25 * s_d, 4)
        }]
      message("[PRECOMP] Done: ", n_with_dosage, " windows upgraded to dosage-based inv_likeness")
    }
  }
} else {
  message("[PRECOMP] No --dosage_dir provided — dosage het rate columns will be NA")
}

# ── ISSUE_017: 3-track summary with adaptive + fixed thresholds ──
n_struct <- sum(inv_like_dt$structure_likeness >= 0.5, na.rm = TRUE)
n_inv <- sum(inv_like_dt$inv_likeness >= 0.5, na.rm = TRUE)
n_fam <- sum(inv_like_dt$family_likeness >= 0.5, na.rm = TRUE)
n_dos <- sum(is.finite(inv_like_dt$dosage_het_rate_median))

# Per-chromosome adaptive: count windows above chr-specific q75
n_inv_adaptive <- inv_like_dt[, sum(inv_likeness >= quantile(inv_likeness, 0.75, na.rm = TRUE),
                                     na.rm = TRUE), by = chrom][, sum(V1)]

message("[PRECOMP] \u2550\u2550 THREE-TRACK SUMMARY \u2550\u2550")
message("  Structure >= 0.5 (fixed): ", n_struct, " / ", nrow(inv_like_dt))
message("  Inv-like  >= 0.5 (fixed): ", n_inv, " / ", nrow(inv_like_dt))
message("  Inv-like  >= q75 (adaptive per-chr): ", n_inv_adaptive, " / ", nrow(inv_like_dt))
message("  Family    >= 0.5 (fixed): ", n_fam, " / ", nrow(inv_like_dt), if (is.null(q_mat)) " (Q not loaded)" else "")
message("  Dosage het computed: ", n_dos, " / ", nrow(inv_like_dt))

# Windows that are structured + NOT family = best inversion candidates
if (n_struct > 0 && n_fam > 0) {
  n_candidate <- sum(inv_like_dt$structure_likeness >= 0.5 &
                      inv_like_dt$family_likeness < 0.3, na.rm = TRUE)
  message("  Structured AND low-family (fixed 0.5/0.3): ", n_candidate,
          " (best pre-snake inversion candidates)")
  message("  NOTE (ISSUE_017): fixed thresholds are reported for comparability;")
  message("    per-chr adaptive thresholds (q75 inv, q25 family) are recommended")
  message("    for production candidate selection.")
}

f_inv <- file.path(outdir, "snake_inv_likeness.tsv.gz")

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

fwrite(inv_like_dt, f_inv, sep = "\t")
message("[PRECOMP] → ", f_inv)

# =============================================================================
# HELPER FUNCTIONS
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

  # Merge inv_likeness
  if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
    dt <- merge(dt, inv_like_dt[, .(global_window_id, inv_likeness)],
                by = "global_window_id", all.x = TRUE)
    dt <- dt[order(start_bp)]
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

  # Similarity matrix (raw = nn0)
  sim_mat <- make_sim_mat(dmat)

  # ── BACKGROUND CONTINUITY BASELINE ──────────────────────────────────
  bg_continuity <- numeric(0)
  chr_inv_med <- median(dt$inv_likeness[is.finite(dt$inv_likeness)])
  n_w <- nrow(dt)

  if (n_w >= 10 && is.finite(chr_inv_med)) {
    bg_pairs <- which(
      dt$inv_likeness[seq_len(n_w - 1)] < chr_inv_med &
      dt$inv_likeness[seq(2, n_w)] < chr_inv_med
    )
    if (length(bg_pairs) >= 20) {
      if (length(bg_pairs) > 2000) bg_pairs <- sort(sample(bg_pairs, 2000))
      bg_continuity <- vapply(bg_pairs, function(i) {
        j <- i + 1L
        s_prev <- if (is.finite(sim_mat[i, j])) sim_mat[i, j] else 0
        s_prev
      }, numeric(1))
    }
  }

  if (length(bg_continuity) >= 10) {
    bg_q <- quantile(bg_continuity, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  } else {
    adj_sims <- vapply(seq_len(n_w - 1), function(i) sim_mat[i, i + 1L], numeric(1))
    bg_q <- quantile(adj_sims, c(0.50, 0.75, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
  }

  message("[PRECOMP] ", chr, ": bg_continuity: q50=", round(bg_q["50%"], 3),
          " q85=", round(bg_q["85%"], 3), " q95=", round(bg_q["95%"], 3),
          " (", length(bg_continuity), " bg pairs)")

  # ── Seed NN distances ──────────────────────────────────────────────
  nn_dists <- vapply(seq_len(nrow(dmat)), function(i) {
    d <- dmat[i, ]; d[i] <- Inf; d <- d[is.finite(d)]
    if (length(d) == 0) return(Inf)
    mean(sort(d)[seq_len(min(SEED_NEIGHBOR_K, length(d)))], na.rm = TRUE)
  }, numeric(1))

  nn_finite <- nn_dists[is.finite(nn_dists)]
  nn_anchor <- if (length(nn_finite) > 0) median(nn_finite) else 1
  if (!is.finite(nn_anchor) || nn_anchor == 0) nn_anchor <- 1
  dt[, seed_nn_dist := nn_dists / nn_anchor]

  nn_norm_info <- list(
    method = "median_nn",
    anchor_value = round(nn_anchor, 6),
    description = "seed_nn_dist = raw_nn / median(all_nn). Low values = tightly clustered."
  )

  # ── Recover chr_bands ──────────────────────────────────────────────
  chr_bands_save <- NULL
  pc1_cols_chr <- grep("^PC_1_", names(dt), value = TRUE)
  if (length(pc1_cols_chr) >= 20) {
    pc1_mat_chr <- as.matrix(dt[, ..pc1_cols_chr])
    avg_pc1_s <- colMeans(pc1_mat_chr, na.rm = TRUE)
    valid_bs <- is.finite(avg_pc1_s)
    if (sum(valid_bs) >= 20) {
      avg_valid_s <- avg_pc1_s[valid_bs]
      chr_km_s <- tryCatch(kmeans(avg_valid_s, centers = 3, nstart = 10), error = function(e) NULL)
      if (!is.null(chr_km_s)) {
        co_s <- order(chr_km_s$centers[, 1])
        chr_bands_save <- integer(length(avg_valid_s))
        chr_bands_save[chr_km_s$cluster == co_s[1]] <- 1L
        chr_bands_save[chr_km_s$cluster == co_s[2]] <- 2L
        chr_bands_save[chr_km_s$cluster == co_s[3]] <- 3L
        names(chr_bands_save) <- sub("^PC_1_", "", names(avg_valid_s))
      }
    }
  }

  # ── Save sim_mat as SEPARATE RDS files ─────────────────────────────
  sim_dir <- file.path(precomp_dir, "sim_mats")
  dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(sim_mat, file.path(sim_dir, paste0(chr, ".sim_mat_nn0.rds")))
  message("[PRECOMP] ", chr, ": saved sim_mat_nn0")

  # ── NN-SMOOTHED SIMILARITY MATRICES ────────────────────────────────
  nn_sim_scales <- as.integer(strsplit(
    Sys.getenv("NN_SIM_SCALES", "20,40,80"), ",")[[1]])
  nn_sim_scales <- nn_sim_scales[nn_sim_scales > 0]

  mds_cols_nn <- grep("^MDS[0-9]+$", names(dt), value = TRUE)
  if (length(mds_cols_nn) > 0 && length(nn_sim_scales) > 0) {
    mds_mat_nn <- as.matrix(dt[, ..mds_cols_nn])
    n_mds <- nrow(mds_mat_nn)
    for (k in nn_sim_scales) {
      k_use <- min(k, n_mds - 1L)
      t_nn <- proc.time()
      smoothed <- matrix(0, nrow = n_mds, ncol = ncol(mds_mat_nn))
      for (wi in seq_len(n_mds)) {
        d <- dmat[wi, ]; d[wi] <- Inf
        nn_idx <- order(d)[seq_len(k_use)]
        smoothed[wi, ] <- colMeans(mds_mat_nn[c(wi, nn_idx), , drop = FALSE], na.rm = TRUE)
      }
      nn_dmat <- as.matrix(dist(smoothed))
      nn_sim <- make_sim_mat(nn_dmat)
      saveRDS(nn_sim, file.path(sim_dir, paste0(chr, ".sim_mat_nn", k, ".rds")))
      elapsed_nn <- round((proc.time() - t_nn)[3], 1)
      message("[PRECOMP] ", chr, ": saved sim_mat_nn", k, " (", elapsed_nn, "s)")
    }
  }
  nn_scales_stored <- c(0L, nn_sim_scales)

  # ── Integrate instant Q (v9.0 wiring) ─────────────────────────────
  if (has_instant_q) {
    tryCatch({
      q_summ <- get_Q_summary(chr)
      if (is.data.table(q_summ) && nrow(q_summ) > 0 && "mean_delta12" %in% names(q_summ)) {
        q_summ[, mid_bp := as.integer((start_bp + end_bp) / 2)]
        dt[, mid_bp := as.integer((start_bp + end_bp) / 2)]
        setkey(dt, mid_bp); setkey(q_summ, mid_bp)
        matched <- q_summ[dt, roll = "nearest", mult = "first", on = "mid_bp"]
        if (nrow(matched) > 0) {
          dt[, localQ_delta12 := matched$mean_delta12]
          dt[, localQ_entropy := matched$mean_entropy]
          dt[, localQ_ena     := matched$mean_ena]
          message("[PRECOMP] ", chr, ": instant Q merged: ",
                  sum(is.finite(dt$localQ_delta12)), "/", nrow(dt))
        }
        dt[, mid_bp := NULL]
      }
    }, error = function(e) message("[PRECOMP] ", chr, ": instant Q: ", e$message))
  }

  # ── Save precomp RDS (v9.1 schema) ────────────────────────────────
  precomp <- list(
    chrom                   = chr,
    n_windows               = nrow(dt),
    schema_version          = "v9.1",
    dt                      = dt,
    sim_mat                 = NULL,
    mds_mat                 = mds_mat,
    dmat                    = dmat,
    bg_continuity_quantiles = bg_q,
    chr_bands               = chr_bands_save,
    nn_normalization        = nn_norm_info,
    sample_names            = sample_names,
    nn_sim_scales           = nn_scales_stored,
    sim_mat_dir             = sim_dir,
    name_format             = name_format,
    has_instant_q           = has_instant_q
  )
  rds_out <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  saveRDS(precomp, rds_out)


  elapsed <- round((proc.time() - t_chr)[3], 1)
  message("[PRECOMP] ", chr, ": ", nrow(dt), " windows, schema=v9.1 (",
          length(names(precomp)), " fields, ", elapsed, "s) -> ", basename(rds_out))

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
    nn_norm_method = "median_nn",
    nn_norm_anchor = round(nn_anchor, 6),
    nn_scales = paste(nn_scales_stored, collapse = ","),
    name_format = name_format,
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

message("\n[DONE] Precompute complete (schema v9.1, bridge=", has_bridge, ", name_format=", name_format, ")")
message("  Chromosomes: ", length(chroms), " (", N_CORES, " cores, ", elapsed_total, "s total)")
message("  Total windows: ", sum(summary_dt$n_windows))
message("  Windows with inv_likeness >= 0.5: ", sum(summary_dt$n_inv_like_050))
message("  Windows with robust |z| >= 2.0: ", sum(summary_dt$n_z_above_2))
message("  Windows with robust |z| >= 3.0: ", sum(summary_dt$n_z_above_3))
message("  NN normalization: median_nn (ISSUE_018)")
message("  Precomp dir: ", precomp_dir)
message("  Inv-likeness: ", f_inv)
message("  Summary: ", file.path(outdir, "precomp_summary.tsv"))
message("\n  Next: run STEP_C01b_snake1_run.R with desired parameters")

# =============================================================================
# REGISTER CANONICAL GROUPS (v9.0 wiring)
# =============================================================================

if (has_bridge && exists("reg") && !is.null(reg)) {
  tryCatch({
    register_group("all_226", sample_names,
                    description = "All 226 QC-pass samples (from precomp v9.1)")
    pruned <- get_pruned_ids()
    if (length(pruned) > 0) {
      register_group("unrelated_81", pruned,
                      description = "NAToRA pruned unrelated set")
    }
    message("[PRECOMP] Registered canonical groups in sample registry")
  }, error = function(e) {
    message("[PRECOMP] Group registration failed: ", e$message)
  })
}
