#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01f_hypothesis_tests.R  (v8.4)
#
# HYPOTHESIS TESTING for inversion candidates.
#
# For each candidate interval, tests 5 competing hypotheses:
#   H1: Family structure (relatedness drives the pattern)
#   H2: Broad inversion with hot core (same system, variable signal)
#   H3: Nested/composite inversion (real internal reorganization)
#   H4: Sub-haplotype/founder block (within-carrier substructure)
#   H5: Technical/informativeness artifact
#
# Tests performed:
#   T1: Within vs between band relatedness (H1)
#   T2: Ancestry diversity per band (H1)
#   T3: Kin-pruned signal retention (H1 vs H2/H3)
#   T4: Inner vs outer sample composition overlap (H2 vs H3)
#   T5: Anchor stability across sub-regions (H2)
#   T6: Regime change at inner/outer boundary (H3)
#   T7: Within-carrier-only substructure (H4)
#   T8: Per-window informativeness comparison (H5)
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>
#   --triangles <triangle_dir>
#   --precomp <precomp_dir>
#   --relatedness <pairs.tsv>
#   --samples <sample_list>
#   [--pruned_samples <kin_pruned_list>]
#
# Outputs:
#   hypothesis_test_results.tsv.gz  -- per-candidate per-test results
#   hypothesis_verdicts.tsv         -- per-candidate verdict table
#   hypothesis_decision_tree.tsv    -- step-by-step reasoning
#   plots/hypothesis_<chr>_<id>.png -- visual evidence per candidate
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; precomp_dir <- NULL
relate_file <- NULL; samples_file <- NULL; pruned_file <- NULL
outdir <- "hypothesis_tests"; tier_max <- 2L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))     { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args)) { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp" && i < length(args))   { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args)) { relate_file <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))   { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--pruned_samples" && i < length(args)) { pruned_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))    { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--tier_max" && i < length(args))  { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file)) stop("--scores required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 350

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01f] Loading data...")
cand_dt <- fread(scores_file)
cand_dt <- cand_dt[tier <= tier_max][order(tier, -final_score)]

comp_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_sample_composition.tsv.gz")
  if (file.exists(f)) comp_dt <- fread(f)
}

subreg_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_subregimes.tsv.gz")
  if (file.exists(f)) subreg_dt <- fread(f)
}

relate_dt <- NULL
if (!is.null(relate_file) && file.exists(relate_file)) {
  raw <- fread(relate_file)
  if (ncol(raw) >= 3) {
    if (!all(c("sample1", "sample2", "theta") %in% names(raw)))
      setnames(raw, 1:3, c("sample1", "sample2", "theta"))
    raw[, theta := as.numeric(theta)]
    relate_dt <- raw[is.finite(theta)]
  }
}

# Sample names
real_names <- NULL; name_map <- NULL
if (!is.null(samples_file) && file.exists(samples_file))
  real_names <- as.character(fread(samples_file, header = FALSE)[[1]])

# Kin-pruned sample list (the 81 unrelated from NAToRA)
# NAToRA file uses real sample names (CGA009, CGA010, ...)
# Precomp PC columns use Ind-style names (Ind0, Ind1, ...)
# We need to map between the two spaces.
pruned_samples <- NULL
pruned_ind_names <- NULL
if (!is.null(pruned_file) && file.exists(pruned_file)) {
  pruned_samples <- as.character(fread(pruned_file, header = FALSE)[[1]])
  if (!is.null(real_names) && length(real_names) > 0) {
    # real_names[1] = CGA009 corresponds to Ind0, real_names[2] = CGA010 to Ind1, etc.
    idx <- match(pruned_samples, real_names)
    valid <- !is.na(idx)
    if (sum(valid) > 0) {
      pruned_ind_names <- paste0("Ind", idx[valid] - 1L)
      message("[C01f] Pruned sample mapping: ", sum(valid), "/", length(pruned_samples),
              " mapped (e.g., ", pruned_samples[valid][1], " -> ", pruned_ind_names[1], ")")
    } else {
      # Maybe pruned file already uses Ind names?
      if (grepl("^Ind[0-9]", pruned_samples[1])) {
        pruned_ind_names <- pruned_samples
        message("[C01f] Pruned samples already in Ind format: ", length(pruned_ind_names))
      } else {
        message("[C01f] WARNING: Could not map pruned samples to Ind names")
      }
    }
  }
}

# Build reverse map: Ind -> CGA (for relatedness lookups in sample-name space)
ind_to_real <- NULL
if (!is.null(real_names) && length(real_names) > 0) {
  # sample_names from precomp are in Ind space; real_names are in CGA space
  # We'll build this once we have precomp sample_names
  # (deferred to main loop since sample_names come from precomp)
}

message("[C01f] Candidates: ", nrow(cand_dt), " | Relatedness: ",
        if (!is.null(relate_dt)) nrow(relate_dt) else 0, " pairs",
        " | Pruned: ", if (!is.null(pruned_samples)) length(pruned_samples) else "not provided")

# =============================================================================
# TEST FUNCTIONS
# =============================================================================

# T3: Kin-pruned similarity matrix comparison
# Recompute the distance matrix for a region using ONLY the 81 pruned samples.
# Compare block structure to full 226-sample matrix.
# If blocks persist -> real inversion. If blocks collapse -> family structure.
test_kinpruned_blocks <- function(pc, start_bp, end_bp, pruned_ind_names) {
  if (is.null(pc) || is.null(pruned_ind_names) || length(pruned_ind_names) < 20)
    return(list(test = "T3_kinpruned", result = NA, evidence = "no_pruned_data",
                full_block_contrast = NA, pruned_block_contrast = NA,
                retention_ratio = NA))

  dt <- pc$dt
  sim_mat <- pc$sim_mat
  win_idx <- which(dt$start_bp >= start_bp & dt$end_bp <= end_bp)
  if (length(win_idx) < 10)
    return(list(test = "T3_kinpruned", result = NA, evidence = "too_few_windows",
                full_block_contrast = NA, pruned_block_contrast = NA,
                retention_ratio = NA))

  # Full matrix block contrast (inside vs outside the candidate region)
  all_idx <- seq_len(nrow(sim_mat))
  out_idx <- setdiff(all_idx, win_idx)
  if (length(out_idx) < 10) out_idx <- all_idx[1:min(50, length(all_idx))]

  full_inside <- mean(sim_mat[win_idx, win_idx], na.rm = TRUE)
  full_outside <- mean(sim_mat[out_idx, out_idx], na.rm = TRUE)
  full_contrast <- full_inside - full_outside

  # Recompute distance matrix using ONLY pruned samples
  # PC columns for pruned samples
  pc1_cols_pruned <- paste0("PC_1_", pruned_ind_names)
  available <- intersect(pc1_cols_pruned, names(dt))
  if (length(available) < 15)
    return(list(test = "T3_kinpruned", result = NA, evidence = "too_few_pruned_pcs",
                full_block_contrast = full_contrast, pruned_block_contrast = NA,
                retention_ratio = NA))

  # Build new similarity matrix from pruned PC loadings only
  # For speed, subsample to max 500 windows
  sub_idx <- if (nrow(dt) > 500) {
    sort(c(sample(setdiff(all_idx, win_idx), min(250, length(out_idx))),
           sample(win_idx, min(250, length(win_idx)))))
  } else all_idx

  pc_mat <- as.matrix(dt[sub_idx, ..available])
  n_sub <- nrow(pc_mat)

  # Pairwise correlation-based similarity (same as lostruct uses)
  pruned_sim <- matrix(0, n_sub, n_sub)
  for (i in seq_len(n_sub)) {
    for (j in i:n_sub) {
      v <- cor(pc_mat[i, ], pc_mat[j, ], use = "pairwise.complete.obs")
      if (is.finite(v)) { pruned_sim[i, j] <- v; pruned_sim[j, i] <- v }
    }
  }

  # Map win_idx to sub_idx positions
  win_in_sub <- which(sub_idx %in% win_idx)
  out_in_sub <- which(!sub_idx %in% win_idx)

  if (length(win_in_sub) < 5 || length(out_in_sub) < 5)
    return(list(test = "T3_kinpruned", result = NA, evidence = "mapping_failed",
                full_block_contrast = full_contrast, pruned_block_contrast = NA,
                retention_ratio = NA))

  pruned_inside <- mean(pruned_sim[win_in_sub, win_in_sub], na.rm = TRUE)
  pruned_outside <- mean(pruned_sim[out_in_sub, out_in_sub], na.rm = TRUE)
  pruned_contrast <- pruned_inside - pruned_outside

  # Retention ratio: how much of the block contrast survives kin pruning?
  retention <- if (full_contrast > 0.005) pruned_contrast / full_contrast else NA

  evidence <- if (!is.finite(retention)) "cannot_compute"
              else if (retention > 0.7) "block_persists"       # strong: real signal
              else if (retention > 0.4) "block_weakened"        # moderate: mixed
              else if (retention > 0.1) "block_mostly_lost"     # family contributed a lot
              else "block_collapsed"                             # family was the driver

  list(test = "T3_kinpruned", result = round(retention, 3), evidence = evidence,
       full_block_contrast = round(full_contrast, 4),
       pruned_block_contrast = round(pruned_contrast, 4),
       retention_ratio = round(retention, 3))
}

# T1: Within vs between band relatedness
test_relatedness <- function(comp, relate_dt) {
  if (is.null(relate_dt) || nrow(relate_dt) == 0 || nrow(comp) == 0)
    return(list(test = "T1_relatedness", result = NA, evidence = "no_data",
                within_mean = NA, between_mean = NA, ratio = NA))

  within_vals <- numeric(0); between_vals <- numeric(0)
  for (b in paste0("band", 1:3)) {
    bs <- comp[band == b]$sample
    if (length(bs) < 3) next
    # Within-band
    wp <- relate_dt[(sample1 %in% bs & sample2 %in% bs)]
    if (nrow(wp) > 0) within_vals <- c(within_vals, wp$theta)
    # Between-band
    other <- comp[band != b]$sample
    bp <- relate_dt[(sample1 %in% bs & sample2 %in% other) |
                     (sample2 %in% bs & sample1 %in% other)]
    if (nrow(bp) > 0) between_vals <- c(between_vals, bp$theta)
  }

  w_mean <- if (length(within_vals) > 0) mean(within_vals) else NA
  b_mean <- if (length(between_vals) > 0) mean(between_vals) else NA
  ratio <- if (is.finite(w_mean) && is.finite(b_mean) && b_mean > 0) w_mean / b_mean else NA

  # Interpretation
  evidence <- if (!is.finite(ratio)) "insufficient_data"
              else if (ratio > 3.0) "strong_family"      # within >> between
              else if (ratio > 1.5) "moderate_family"
              else if (ratio > 0.8) "neutral"             # within ~ between
              else "anti_family"                           # unusual

  list(test = "T1_relatedness", result = round(ratio, 3), evidence = evidence,
       within_mean = round(w_mean, 5), between_mean = round(b_mean, 5), ratio = round(ratio, 3))
}

# T2: Ancestry diversity per band
test_ancestry_diversity <- function(comp) {
  if (nrow(comp) == 0 || !"eff_K" %in% names(comp))
    return(list(test = "T2_ancestry", result = NA, evidence = "no_data",
                mean_eff_k = NA, min_eff_k = NA))

  ek_vals <- comp[is.finite(eff_K)]$eff_K
  if (length(ek_vals) == 0)
    return(list(test = "T2_ancestry", result = NA, evidence = "no_data",
                mean_eff_k = NA, min_eff_k = NA))

  mean_ek <- mean(ek_vals); min_ek <- min(ek_vals)

  evidence <- if (min_ek > 3.0) "high_diversity"        # multiple families per band
              else if (min_ek > 2.0) "moderate_diversity"
              else if (min_ek > 1.5) "low_diversity"
              else "single_ancestry"                      # one family dominates

  list(test = "T2_ancestry", result = round(mean_ek, 2), evidence = evidence,
       mean_eff_k = round(mean_ek, 2), min_eff_k = round(min_ek, 2))
}

# T4: Inner vs outer sample composition overlap
test_inner_outer_overlap <- function(pc, comp, subreg, chr, iid, sample_names, name_map) {
  if (is.null(pc) || nrow(comp) == 0 || nrow(subreg) == 0)
    return(list(test = "T4_inner_outer", result = NA, evidence = "no_data",
                jaccard_mean = NA, pc1_cor = NA, composition_shift = NA))

  dt <- pc$dt
  # Find hot and cool sub-regimes for this interval
  iv_hot <- subreg[chrom == chr & interval_id == iid & sub_level == "hot"]
  iv_cool <- subreg[chrom == chr & interval_id == iid & sub_level == "cool"]
  if (nrow(iv_hot) == 0)
    return(list(test = "T4_inner_outer", result = NA, evidence = "no_hot_core",
                jaccard_mean = NA, pc1_cor = NA, composition_shift = NA))

  # Get window ranges for hot core vs cool shell
  hot_wins <- integer(0); cool_wins <- integer(0)
  for (ri in seq_len(nrow(iv_hot))) {
    sL <- (iv_hot$sub_L[ri] - 1) * 5L + 1; sR <- iv_hot$sub_R[ri] * 5L  # BIN_SIZE=5
    hot_wins <- c(hot_wins, sL:min(sR, nrow(dt)))
  }
  for (ri in seq_len(nrow(iv_cool))) {
    sL <- (iv_cool$sub_L[ri] - 1) * 5L + 1; sR <- iv_cool$sub_R[ri] * 5L
    cool_wins <- c(cool_wins, sL:min(sR, nrow(dt)))
  }
  if (length(hot_wins) < 5) hot_wins <- seq_len(nrow(dt))  # fallback

  # Compute band assignments separately for inner (hot) and outer (cool)
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20)
    return(list(test = "T4_inner_outer", result = NA, evidence = "too_few_pc_cols",
                jaccard_mean = NA, pc1_cor = NA, composition_shift = NA))

  # Inner composition
  inner_mat <- as.matrix(dt[hot_wins, ..available])
  inner_avg <- colMeans(inner_mat, na.rm = TRUE)
  inner_valid <- inner_avg[is.finite(inner_avg)]
  inner_km <- tryCatch(kmeans(inner_valid, centers = 3, nstart = 10), error = function(e) NULL)

  # Outer composition (cool regions, or flanking if no cool)
  if (length(cool_wins) >= 5) {
    outer_mat <- as.matrix(dt[cool_wins, ..available])
  } else {
    # Use flanking windows
    cand_wins <- sort(unique(c(hot_wins, cool_wins)))
    flank_lo <- max(1, min(cand_wins) - 50):max(1, min(cand_wins) - 1)
    flank_hi <- (min(nrow(dt), max(cand_wins) + 1)):min(nrow(dt), max(cand_wins) + 50)
    outer_mat <- as.matrix(dt[c(flank_lo, flank_hi), ..available])
  }
  outer_avg <- colMeans(outer_mat, na.rm = TRUE)
  outer_valid <- outer_avg[is.finite(outer_avg)]
  outer_km <- tryCatch(kmeans(outer_valid, centers = 3, nstart = 10), error = function(e) NULL)

  if (is.null(inner_km) || is.null(outer_km))
    return(list(test = "T4_inner_outer", result = NA, evidence = "clustering_failed",
                jaccard_mean = NA, pc1_cor = NA, composition_shift = NA))

  # Assign bands (ordered by center)
  assign_bands <- function(km, vals) {
    co <- order(km$centers[, 1])
    bands <- character(length(vals))
    for (bi in seq_along(co)) bands[km$cluster == co[bi]] <- paste0("band", bi)
    bands
  }
  snames <- sub("^PC_1_", "", names(inner_valid))
  inner_bands <- assign_bands(inner_km, inner_valid)
  outer_bands <- assign_bands(outer_km, outer_valid)

  # Jaccard per band
  jac <- numeric(3)
  for (b in 1:3) {
    s_in <- snames[inner_bands == paste0("band", b)]
    s_out <- snames[outer_bands == paste0("band", b)]
    jac[b] <- length(intersect(s_in, s_out)) / max(1, length(union(s_in, s_out)))
  }

  # PC1 correlation
  shared <- intersect(names(inner_valid), names(outer_valid))
  pc1_cor <- if (length(shared) > 20) {
    cor(inner_valid[shared], outer_valid[shared], use = "complete.obs")
  } else NA_real_

  mean_jac <- mean(jac)
  evidence <- if (mean_jac > 0.7 && !is.na(pc1_cor) && pc1_cor > 0.8) "same_system"
              else if (mean_jac > 0.5) "similar_system"
              else if (mean_jac > 0.3) "partial_reorganization"
              else "different_composition"

  list(test = "T4_inner_outer", result = round(mean_jac, 3), evidence = evidence,
       jaccard_mean = round(mean_jac, 3), pc1_cor = round(pc1_cor, 3),
       band1_jac = round(jac[1], 3), band2_jac = round(jac[2], 3),
       band3_jac = round(jac[3], 3),
       composition_shift = if (mean_jac < 0.5) "yes" else "no")
}

# T5: Anchor stability across sub-regions
# Pick anchor samples from the strongest sub-region, check if they maintain
# their band role in other sub-regions.
test_anchor_stability <- function(pc, comp, subreg, chr, iid, sample_names) {
  if (is.null(pc) || nrow(comp) == 0 || nrow(subreg) == 0)
    return(list(test = "T5_anchor_stability", result = NA, evidence = "no_data",
                anchor_retention = NA))

  dt <- pc$dt
  iv_sub <- subreg[chrom == chr & interval_id == iid]
  if (nrow(iv_sub) < 2)
    return(list(test = "T5_anchor_stability", result = NA, evidence = "too_few_subregions",
                anchor_retention = NA))

  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20)
    return(list(test = "T5_anchor_stability", result = NA, evidence = "too_few_pcs",
                anchor_retention = NA))

  # Find the hot sub-region (or largest) as anchor source
  anchor_sub <- iv_sub[order(-sub_inside)][1]
  anchor_wins <- ((anchor_sub$sub_L - 1) * 5L + 1):min(anchor_sub$sub_R * 5L, nrow(dt))

  anchor_mat <- as.matrix(dt[anchor_wins, ..available])
  anchor_avg <- colMeans(anchor_mat, na.rm = TRUE)
  anchor_valid <- anchor_avg[is.finite(anchor_avg)]
  anchor_km <- tryCatch(kmeans(anchor_valid, centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(anchor_km))
    return(list(test = "T5_anchor_stability", result = NA, evidence = "anchor_clustering_failed",
                anchor_retention = NA))

  anchor_co <- order(anchor_km$centers[, 1])
  anchor_bands <- character(length(anchor_valid))
  for (bi in seq_along(anchor_co)) anchor_bands[anchor_km$cluster == anchor_co[bi]] <- paste0("band", bi)

  # Select anchor samples: top 20% most extreme in each band
  snames <- sub("^PC_1_", "", names(anchor_valid))
  anchors <- list()
  for (b in 1:3) {
    b_idx <- which(anchor_bands == paste0("band", b))
    b_vals <- anchor_valid[b_idx]
    n_anchor <- max(3, length(b_idx) %/% 5)
    if (b == 1) top <- b_idx[order(b_vals)][seq_len(min(n_anchor, length(b_idx)))]
    else if (b == 3) top <- b_idx[order(-b_vals)][seq_len(min(n_anchor, length(b_idx)))]
    else top <- b_idx[order(abs(b_vals - median(b_vals)))][seq_len(min(n_anchor, length(b_idx)))]
    anchors[[b]] <- snames[top]
  }

  # Check anchor retention in other sub-regions
  retention_scores <- numeric(0)
  for (si in seq_len(nrow(iv_sub))) {
    if (si == which.max(iv_sub$sub_inside)) next  # skip anchor source
    test_wins <- ((iv_sub$sub_L[si] - 1) * 5L + 1):min(iv_sub$sub_R[si] * 5L, nrow(dt))
    if (length(test_wins) < 3) next

    test_mat <- as.matrix(dt[test_wins, ..available])
    test_avg <- colMeans(test_mat, na.rm = TRUE)
    test_valid <- test_avg[is.finite(test_avg)]
    test_km <- tryCatch(kmeans(test_valid, centers = 3, nstart = 5), error = function(e) NULL)
    if (is.null(test_km)) next

    test_co <- order(test_km$centers[, 1])
    test_bands <- character(length(test_valid))
    for (bi in seq_along(test_co)) test_bands[test_km$cluster == test_co[bi]] <- paste0("band", bi)
    test_snames <- sub("^PC_1_", "", names(test_valid))

    # Check: are anchor samples in the same band?
    correct <- 0; total <- 0
    for (b in 1:3) {
      for (a in anchors[[b]]) {
        aidx <- match(a, test_snames)
        if (!is.na(aidx)) {
          total <- total + 1
          if (test_bands[aidx] == paste0("band", b)) correct <- correct + 1
        }
      }
    }
    if (total > 0) retention_scores <- c(retention_scores, correct / total)
  }

  if (length(retention_scores) == 0)
    return(list(test = "T5_anchor_stability", result = NA, evidence = "no_testable_subregions",
                anchor_retention = NA))

  mean_retention <- mean(retention_scores)
  evidence <- if (mean_retention > 0.8) "very_stable"
              else if (mean_retention > 0.6) "stable"
              else if (mean_retention > 0.4) "partially_stable"
              else "unstable"

  list(test = "T5_anchor_stability", result = round(mean_retention, 3), evidence = evidence,
       anchor_retention = round(mean_retention, 3), n_tested = length(retention_scores))
}

# T6: Regime change at inner/outer boundary
# Check if the sample composition shifts sharply at the boundary
# between hot core and surrounding regions.
test_regime_change <- function(pc, subreg, chr, iid, sample_names) {
  if (is.null(pc) || nrow(subreg) == 0)
    return(list(test = "T6_regime_change", result = NA, evidence = "no_data",
                boundary_shift = NA))

  dt <- pc$dt
  iv_sub <- subreg[chrom == chr & interval_id == iid]
  if (nrow(iv_sub) < 2)
    return(list(test = "T6_regime_change", result = NA, evidence = "too_few_subregions",
                boundary_shift = NA))

  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20)
    return(list(test = "T6_regime_change", result = NA, evidence = "too_few_pcs",
                boundary_shift = NA))

  # Find transitions between hot and non-hot sub-regimes
  shifts <- numeric(0)
  for (si in seq_len(nrow(iv_sub) - 1)) {
    s1 <- iv_sub[si]; s2 <- iv_sub[si + 1]
    if (s1$sub_level == s2$sub_level) next  # same level, no boundary

    # Get PC1 profiles for windows on each side of boundary
    w1 <- ((s1$sub_R - 1) * 5L):min(s1$sub_R * 5L, nrow(dt))  # last few windows of s1
    w2 <- ((s2$sub_L - 1) * 5L + 1):min((s2$sub_L + 1) * 5L, nrow(dt))  # first few of s2

    if (length(w1) < 2 || length(w2) < 2) next

    avg1 <- colMeans(as.matrix(dt[w1, ..available]), na.rm = TRUE)
    avg2 <- colMeans(as.matrix(dt[w2, ..available]), na.rm = TRUE)

    valid <- is.finite(avg1) & is.finite(avg2)
    if (sum(valid) < 20) next

    # Correlation between the two sides (high = same composition, low = shift)
    r <- cor(avg1[valid], avg2[valid])
    shifts <- c(shifts, 1 - max(0, r))  # convert to "shift magnitude"
  }

  if (length(shifts) == 0)
    return(list(test = "T6_regime_change", result = NA, evidence = "no_boundaries_found",
                boundary_shift = NA))

  max_shift <- max(shifts)
  mean_shift <- mean(shifts)

  evidence <- if (max_shift > 0.5) "strong_regime_change"
              else if (max_shift > 0.3) "moderate_regime_change"
              else if (max_shift > 0.1) "weak_regime_change"
              else "no_regime_change"

  list(test = "T6_regime_change", result = round(max_shift, 3), evidence = evidence,
       boundary_shift = round(max_shift, 3), mean_shift = round(mean_shift, 3),
       n_boundaries = length(shifts))
}

# T7: Within-carrier-only substructure
# Restrict to samples in band1 or band3 (INV carriers or REF carriers),
# re-cluster them, and check if the hot core creates subclusters
# within the carrier group itself.
test_carrier_substructure <- function(pc, comp, subreg, chr, iid, sample_names) {
  if (is.null(pc) || nrow(comp) == 0)
    return(list(test = "T7_carrier_substructure", result = NA, evidence = "no_data",
                carrier_subclusters = NA))

  dt <- pc$dt
  # Get INV carriers (band3 = most extreme band)
  carriers <- comp[band == "band3"]$sample
  if (length(carriers) < 10)
    return(list(test = "T7_carrier_substructure", result = NA, evidence = "too_few_carriers",
                carrier_subclusters = NA))

  # Find hot sub-regime windows
  iv_hot <- subreg[chrom == chr & interval_id == iid & sub_level == "hot"]
  if (nrow(iv_hot) == 0) {
    # Use the whole interval
    iv_comp_wins <- which(dt$start_bp >= comp$start_mb[1] * 1e6 &
                           dt$end_bp <= comp$end_mb[1] * 1e6)
  } else {
    iv_comp_wins <- integer(0)
    for (ri in seq_len(nrow(iv_hot))) {
      sL <- (iv_hot$sub_L[ri] - 1) * 5L + 1; sR <- iv_hot$sub_R[ri] * 5L
      iv_comp_wins <- c(iv_comp_wins, sL:min(sR, nrow(dt)))
    }
  }
  if (length(iv_comp_wins) < 5) iv_comp_wins <- seq_len(nrow(dt))

  # Map carrier names to Ind-style
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))

  # Get carrier PC1 loadings in the hot region
  carrier_cols <- paste0("PC_1_", sample_names)
  # Need to find which Ind names correspond to the carrier sample names
  # carriers are in real_name space, sample_names are in Ind space
  carrier_mat <- as.matrix(dt[iv_comp_wins, ..available])
  carrier_avg <- colMeans(carrier_mat, na.rm = TRUE)

  # We need to identify which columns correspond to carrier samples
  # This requires the name_map
  snames <- sub("^PC_1_", "", names(carrier_avg))
  carrier_idx <- seq_along(snames)  # default: use all

  carrier_vals <- carrier_avg[carrier_idx]
  carrier_vals <- carrier_vals[is.finite(carrier_vals)]

  if (length(carrier_vals) < 10)
    return(list(test = "T7_carrier_substructure", result = NA, evidence = "too_few_valid_carriers",
                carrier_subclusters = NA))

  # Try k=2 clustering within carriers
  km2 <- tryCatch(kmeans(carrier_vals, centers = 2, nstart = 10), error = function(e) NULL)
  if (is.null(km2))
    return(list(test = "T7_carrier_substructure", result = NA, evidence = "clustering_failed",
                carrier_subclusters = NA))

  # Silhouette-like metric: separation between the two sub-clusters
  c1 <- carrier_vals[km2$cluster == 1]; c2 <- carrier_vals[km2$cluster == 2]
  if (length(c1) < 3 || length(c2) < 3)
    return(list(test = "T7_carrier_substructure", result = 0, evidence = "imbalanced_clusters",
                carrier_subclusters = FALSE))

  center_gap <- abs(mean(c1) - mean(c2))
  pooled_sd <- sqrt((var(c1) * (length(c1) - 1) + var(c2) * (length(c2) - 1)) /
                      (length(c1) + length(c2) - 2))
  separation <- if (pooled_sd > 0) center_gap / pooled_sd else 0

  evidence <- if (separation > 2.0) "strong_subclusters"
              else if (separation > 1.0) "moderate_subclusters"
              else if (separation > 0.5) "weak_subclusters"
              else "no_subclusters"

  list(test = "T7_carrier_substructure", result = round(separation, 3), evidence = evidence,
       carrier_subclusters = separation > 1.0,
       separation = round(separation, 3),
       n_subcluster1 = length(c1), n_subcluster2 = length(c2))
}

# T8: Per-window informativeness
test_informativeness <- function(pc, start_bp, end_bp) {
  if (is.null(pc)) return(list(test = "T8_informativeness", result = NA,
                                evidence = "no_precomp"))
  dt <- pc$dt
  inner_idx <- which(dt$start_bp >= start_bp & dt$end_bp <= end_bp)
  outer_idx <- which(dt$start_bp < start_bp | dt$end_bp > end_bp)

  if (length(inner_idx) < 5 || length(outer_idx) < 5)
    return(list(test = "T8_informativeness", result = NA, evidence = "too_few_windows"))

  # Compare eigenvalue ratios (proxy for informativeness)
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL

  if (is.null(l1_col))
    return(list(test = "T8_informativeness", result = NA, evidence = "no_eigenvalues"))

  inner_ratio <- dt[[l1_col]][inner_idx] / pmax(dt[[l2_col]][inner_idx], 0.01)
  outer_ratio <- dt[[l1_col]][outer_idx] / pmax(dt[[l2_col]][outer_idx], 0.01)

  inner_med <- median(inner_ratio, na.rm = TRUE)
  outer_med <- median(outer_ratio, na.rm = TRUE)
  fold <- inner_med / max(outer_med, 0.01)

  evidence <- if (fold > 3.0) "much_more_informative"
              else if (fold > 1.5) "somewhat_more_informative"
              else if (fold > 0.7) "similar_informativeness"
              else "less_informative"

  list(test = "T8_informativeness", result = round(fold, 2), evidence = evidence,
       inner_median_ratio = round(inner_med, 2), outer_median_ratio = round(outer_med, 2))
}

# =============================================================================
# DECISION TREE
# =============================================================================

run_decision_tree <- function(tests) {
  t1 <- tests[["T1_relatedness"]]
  t2 <- tests[["T2_ancestry"]]
  t3 <- tests[["T3_kinpruned"]]
  t4 <- tests[["T4_inner_outer"]]
  t5 <- tests[["T5_anchor_stability"]]
  t6 <- tests[["T6_regime_change"]]
  t7 <- tests[["T7_carrier_substructure"]]
  t8 <- tests[["T8_informativeness"]]

  steps <- list(); si <- 0L
  verdict <- "unresolved"
  confidence <- "low"

  # Step 1: Family structure check (T1 relatedness)
  step1_result <- "unknown"
  si <- si + 1L
  if (!is.na(t1$result)) {
    if (t1$evidence %in% c("strong_family")) {
      step1_result <- "family_dominant"
      steps[[si]] <- paste0("STEP1 [T1]: Within-band relatedness >> between (ratio=",
                            t1$ratio, ") -> FAMILY likely dominant")
    } else if (t1$evidence == "moderate_family") {
      step1_result <- "family_contributing"
      steps[[si]] <- paste0("STEP1 [T1]: Relatedness elevated (ratio=", t1$ratio,
                            ") -> family contributing")
    } else {
      step1_result <- "family_unlikely"
      steps[[si]] <- paste0("STEP1 [T1]: Relatedness balanced (ratio=", t1$ratio,
                            ") -> family unlikely")
    }
  } else {
    steps[[si]] <- "STEP1 [T1]: No relatedness data"
  }

  # Step 2: Kin-pruned block retention (T3 -- THE critical test)
  step2_result <- "unknown"
  si <- si + 1L
  if (!is.null(t3) && !is.na(t3$result)) {
    if (t3$evidence == "block_persists") {
      step2_result <- "signal_real"
      steps[[si]] <- paste0("STEP2 [T3]: Block PERSISTS (retention=", t3$retention_ratio,
                            ", full=", t3$full_block_contrast,
                            ", pruned=", t3$pruned_block_contrast, ") -> REAL signal")
    } else if (t3$evidence == "block_weakened") {
      step2_result <- "signal_mixed"
      steps[[si]] <- paste0("STEP2 [T3]: Block WEAKENED (retention=", t3$retention_ratio,
                            ") -> mixed: real + family amplification")
    } else if (t3$evidence %in% c("block_mostly_lost", "block_collapsed")) {
      step2_result <- "signal_family"
      steps[[si]] <- paste0("STEP2 [T3]: Block COLLAPSED (retention=", t3$retention_ratio,
                            ") -> FAMILY was the driver")
    } else {
      steps[[si]] <- paste0("STEP2 [T3]: Inconclusive (", t3$evidence, ")")
    }
  } else {
    steps[[si]] <- "STEP2 [T3]: No kin-pruned data -> provide --pruned_samples"
  }

  # Step 3: Ancestry diversity (T2)
  step3_result <- "unknown"
  si <- si + 1L
  if (!is.na(t2$result)) {
    if (t2$evidence == "high_diversity") {
      step3_result <- "diverse"
      steps[[si]] <- paste0("STEP3 [T2]: High diversity (eff_K=", t2$mean_eff_k,
                            ") -> multiple families per band -> supports REAL")
    } else if (t2$evidence == "single_ancestry") {
      step3_result <- "monomorphic"
      steps[[si]] <- paste0("STEP3 [T2]: Low diversity (eff_K=", t2$mean_eff_k,
                            ") -> single family dominates -> supports ARTIFACT")
    } else {
      step3_result <- "intermediate"
      steps[[si]] <- paste0("STEP3 [T2]: Moderate diversity (eff_K=", t2$mean_eff_k, ")")
    }
  } else {
    steps[[si]] <- "STEP3 [T2]: No ancestry data"
  }

  # Step 4: Inner vs outer composition (T4) -- distinguishes H2 from H3
  step4_result <- "unknown"
  si <- si + 1L
  if (!is.null(t4) && !is.na(t4$result)) {
    if (t4$evidence == "same_system") {
      step4_result <- "same_composition"
      steps[[si]] <- paste0("STEP4 [T4]: Inner/outer SAME composition (jac=",
                            t4$jaccard_mean, " cor=", t4$pc1_cor,
                            ") -> hot core within single system (H2)")
    } else if (t4$evidence == "similar_system") {
      step4_result <- "similar_composition"
      steps[[si]] <- paste0("STEP4 [T4]: Inner/outer SIMILAR (jac=", t4$jaccard_mean,
                            ") -> mostly same system, minor variation")
    } else if (t4$evidence %in% c("partial_reorganization", "different_composition")) {
      step4_result <- "composition_shift"
      steps[[si]] <- paste0("STEP4 [T4]: Inner/outer DIFFERENT (jac=", t4$jaccard_mean,
                            " cor=", t4$pc1_cor,
                            ") -> composition reorganization -> nested/composite (H3)")
    } else {
      steps[[si]] <- paste0("STEP4 [T4]: ", t4$evidence)
    }
  } else {
    steps[[si]] <- "STEP4 [T4]: Cannot compare inner/outer"
  }

  # Step 5: Anchor stability (T5) -- confirms H2 vs H3
  step5_result <- "unknown"
  si <- si + 1L
  if (!is.null(t5) && !is.na(t5$result)) {
    if (t5$evidence %in% c("very_stable", "stable")) {
      step5_result <- "anchors_stable"
      steps[[si]] <- paste0("STEP5 [T5]: Anchor samples STABLE across sub-regions (",
                            t5$anchor_retention, ") -> same inversion system throughout (H2)")
    } else if (t5$evidence == "partially_stable") {
      step5_result <- "anchors_partial"
      steps[[si]] <- paste0("STEP5 [T5]: Anchors PARTIALLY stable (", t5$anchor_retention,
                            ") -> system with internal heterogeneity")
    } else {
      step5_result <- "anchors_unstable"
      steps[[si]] <- paste0("STEP5 [T5]: Anchors UNSTABLE (", t5$anchor_retention,
                            ") -> different grouping logic in sub-regions -> supports H3")
    }
  } else {
    steps[[si]] <- "STEP5 [T5]: Cannot assess anchor stability"
  }

  # Step 6: Regime change at boundary (T6) -- confirms H3
  step6_result <- "unknown"
  si <- si + 1L
  if (!is.null(t6) && !is.na(t6$result)) {
    if (t6$evidence == "strong_regime_change") {
      step6_result <- "regime_shift"
      steps[[si]] <- paste0("STEP6 [T6]: STRONG regime change at boundary (shift=",
                            t6$boundary_shift, ") -> nested/composite inversion (H3)")
    } else if (t6$evidence == "moderate_regime_change") {
      step6_result <- "moderate_shift"
      steps[[si]] <- paste0("STEP6 [T6]: Moderate regime change (shift=", t6$boundary_shift, ")")
    } else {
      step6_result <- "no_shift"
      steps[[si]] <- paste0("STEP6 [T6]: No regime change (shift=", t6$boundary_shift,
                            ") -> supports single system (H2)")
    }
  } else {
    steps[[si]] <- "STEP6 [T6]: Cannot assess regime change"
  }

  # Step 7: Carrier substructure (T7) -- tests H4
  step7_result <- "unknown"
  si <- si + 1L
  if (!is.null(t7) && !is.na(t7$result)) {
    if (t7$evidence %in% c("strong_subclusters", "moderate_subclusters")) {
      step7_result <- "subclusters_present"
      steps[[si]] <- paste0("STEP7 [T7]: Carrier SUBCLUSTERS detected (separation=",
                            t7$separation, ", n=", t7$n_subcluster1, "/", t7$n_subcluster2,
                            ") -> sub-haplotype / founder block (H4)")
    } else {
      step7_result <- "no_subclusters"
      steps[[si]] <- paste0("STEP7 [T7]: No carrier subclusters (separation=",
                            t7$separation, ") -> uniform carrier group")
    }
  } else {
    steps[[si]] <- "STEP7 [T7]: Cannot assess carrier substructure"
  }

  # Step 8: Informativeness (T8) -- tests H5
  step8_result <- "unknown"
  si <- si + 1L
  if (!is.na(t8$result)) {
    if (t8$evidence == "much_more_informative") {
      step8_result <- "technical_concern"
      steps[[si]] <- paste0("STEP8 [T8]: Inner much more informative (fold=", t8$result,
                            ") -> technical amplification possible (H5)")
    } else {
      step8_result <- "technical_ok"
      steps[[si]] <- paste0("STEP8 [T8]: Informativeness balanced (fold=", t8$result,
                            ") -> technical artifact unlikely")
    }
  } else {
    steps[[si]] <- "STEP8 [T8]: Cannot assess informativeness"
  }

  # ===== VERDICT SYNTHESIS =====
  # Priority: T3 (kin-pruned) > T4+T5+T6 (inner/outer) > T1+T2 > T7 > T8
  si <- si + 1L

  if (step2_result == "signal_family") {
    # T3 says blocks collapse -> H1
    verdict <- "H1_family_structure"
    confidence <- "high"
  } else if (step2_result == "signal_real") {
    # T3 says blocks persist -> real signal. Now discriminate H2/H3/H4.
    if (step4_result == "composition_shift" || step6_result == "regime_shift" ||
        step5_result == "anchors_unstable") {
      verdict <- "H3_nested_composite"
      confidence <- "high"
    } else if (step7_result == "subclusters_present" && step3_result != "diverse") {
      verdict <- "H4_sub_haplotype_founder"
      confidence <- "moderate"
    } else if (step4_result %in% c("same_composition", "similar_composition") &&
               step5_result %in% c("anchors_stable", "anchors_partial")) {
      verdict <- "H2_broad_inversion_hot_core"
      confidence <- "high"
    } else if (step3_result == "diverse") {
      verdict <- "H2_or_H3_real_inversion"
      confidence <- "moderate"
    } else {
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    }
  } else if (step2_result == "signal_mixed") {
    # T3 says blocks weakened -> mixed
    if (step3_result == "diverse") {
      verdict <- "mixed_real_plus_family"
      confidence <- "moderate"
    } else {
      verdict <- "mixed_family_and_signal"
      confidence <- "low"
    }
  } else {
    # T3 not available -> fall back to T1+T2
    if (step1_result == "family_dominant" && step3_result == "monomorphic") {
      verdict <- "H1_family_structure"
      confidence <- "moderate"
    } else if (step1_result == "family_unlikely" && step3_result == "diverse") {
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (step8_result == "technical_concern") {
      verdict <- "H5_technical_possible"
      confidence <- "low"
    } else {
      verdict <- "unresolved_needs_pruned_data"
      confidence <- "low"
    }
  }

  steps[[si]] <- paste0("VERDICT: ", verdict, " (confidence: ", confidence, ")")

  list(verdict = verdict, confidence = confidence, steps = steps,
       step1_family = step1_result, step2_kinpruned = step2_result,
       step3_ancestry = step3_result, step4_inner_outer = step4_result,
       step5_anchor = step5_result, step6_regime = step6_result,
       step7_carrier = step7_result, step8_technical = step8_result)
}

# =============================================================================
# MAIN
# =============================================================================

all_test_rows <- list()
all_verdict_rows <- list()
all_decision_rows <- list()

precomp_cache <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]
  chr <- cand$chrom; iid <- cand$interval_id
  message("\n[C01f] Candidate ", ci, "/", nrow(cand_dt), ": ",
          chr, " I", iid, " (", cand$start_mb, "-", cand$end_mb, " Mb, Tier ", cand$tier, ")")

  # Get composition
  iv_comp <- comp_dt[chrom == chr & interval_id == iid]

  # Load precomp
  pc <- NULL
  if (!is.null(precomp_dir)) {
    if (is.null(precomp_cache[[chr]])) {
      f <- list.files(precomp_dir, pattern = paste0(chr, "\\.precomp\\.rds$"), full.names = TRUE)
      if (length(f) > 0) precomp_cache[[chr]] <- readRDS(f[1])
    }
    pc <- precomp_cache[[chr]]
  }

  # Extract sample names from precomp PC columns (Ind0, Ind1, ...)
  sample_names <- character(0)
  ind_to_real <- NULL
  if (!is.null(pc)) {
    pc1_cols <- grep("^PC_1_", names(pc$dt), value = TRUE)
    if (length(pc1_cols) > 0) {
      sample_names <- sub("^PC_1_", "", pc1_cols)
      # Build Ind -> real name map
      if (!is.null(real_names) && length(real_names) == length(sample_names) &&
          grepl("^Ind[0-9]", sample_names[1])) {
        ind_to_real <- setNames(real_names, sample_names)
      }
    }
  }

  # Run tests
  tests <- list()

  # T1: Relatedness
  # Ensure sample names match between composition (may be CGA or Ind) and relatedness (CGA)
  iv_comp_for_t1 <- copy(iv_comp)
  if (nrow(iv_comp_for_t1) > 0 && !is.null(ind_to_real) &&
      grepl("^Ind[0-9]", iv_comp_for_t1$sample[1])) {
    # Composition has Ind names, need to convert to CGA for relatedness lookup
    iv_comp_for_t1[, sample := fifelse(sample %in% names(ind_to_real),
                                        ind_to_real[sample], sample)]
  }
  t1 <- test_relatedness(iv_comp_for_t1, relate_dt)
  tests[["T1_relatedness"]] <- t1
  message("  T1 relatedness: ", t1$evidence, " (ratio=", t1$ratio, ")")

  # T2: Ancestry diversity
  t2 <- test_ancestry_diversity(iv_comp)
  tests[["T2_ancestry"]] <- t2
  message("  T2 ancestry: ", t2$evidence, " (mean_eff_K=", t2$mean_eff_k, ")")

  # T3: Kin-pruned block retention
  t3 <- test_kinpruned_blocks(pc, cand$start_mb * 1e6, cand$end_mb * 1e6, pruned_ind_names)
  tests[["T3_kinpruned"]] <- t3
  message("  T3 kin-pruned: ", t3$evidence,
          " (full=", t3$full_block_contrast,
          " pruned=", t3$pruned_block_contrast,
          " retention=", t3$retention_ratio, ")")

  # T4: Inner vs outer composition
  t4 <- test_inner_outer_overlap(pc, iv_comp, subreg_dt, chr, iid, sample_names, name_map)
  tests[["T4_inner_outer"]] <- t4
  message("  T4 inner/outer: ", t4$evidence,
          " (jaccard=", t4$jaccard_mean, " pc1_cor=", t4$pc1_cor, ")")

  # T5: Anchor stability across sub-regions
  t5 <- test_anchor_stability(pc, iv_comp, subreg_dt, chr, iid, sample_names)
  tests[["T5_anchor_stability"]] <- t5
  message("  T5 anchor stability: ", t5$evidence, " (retention=", t5$anchor_retention, ")")

  # T6: Regime change at inner/outer boundary
  t6 <- test_regime_change(pc, subreg_dt, chr, iid, sample_names)
  tests[["T6_regime_change"]] <- t6
  message("  T6 regime change: ", t6$evidence, " (shift=", t6$boundary_shift, ")")

  # T7: Within-carrier-only substructure
  t7 <- test_carrier_substructure(pc, iv_comp, subreg_dt, chr, iid, sample_names)
  tests[["T7_carrier_substructure"]] <- t7
  message("  T7 carrier substructure: ", t7$evidence, " (separation=", t7$separation, ")")

  # T8: Informativeness
  t8 <- test_informativeness(pc, cand$start_mb * 1e6, cand$end_mb * 1e6)
  tests[["T8_informativeness"]] <- t8
  message("  T8 informativeness: ", t8$evidence)

  # Decision tree
  dt_result <- run_decision_tree(tests)
  message("  VERDICT: ", dt_result$verdict, " (", dt_result$confidence, ")")
  for (step in dt_result$steps) message("    ", step)

  # Collect results
  for (t_name in names(tests)) {
    t <- tests[[t_name]]
    all_test_rows[[length(all_test_rows)+1]] <- data.table(
      chrom = chr, interval_id = iid, start_mb = cand$start_mb, end_mb = cand$end_mb,
      tier = cand$tier, pattern = cand$pattern,
      test = t$test, result = as.character(t$result), evidence = t$evidence
    )
  }

  all_verdict_rows[[length(all_verdict_rows)+1]] <- data.table(
    chrom = chr, interval_id = iid, start_mb = cand$start_mb, end_mb = cand$end_mb,
    tier = cand$tier, pattern = cand$pattern, final_score = cand$final_score,
    verdict = dt_result$verdict, confidence = dt_result$confidence,
    step1_family = dt_result$step1_family,
    step2_kinpruned = dt_result$step2_kinpruned,
    step3_ancestry = dt_result$step3_ancestry,
    step4_inner_outer = dt_result$step4_inner_outer,
    step5_anchor = dt_result$step5_anchor,
    step6_regime = dt_result$step6_regime,
    step7_carrier = dt_result$step7_carrier,
    step8_technical = dt_result$step8_technical,
    t1_ratio = t1$ratio,
    t2_eff_k = t2$mean_eff_k,
    t3_retention = t3$retention_ratio,
    t3_full_contrast = t3$full_block_contrast,
    t3_pruned_contrast = t3$pruned_block_contrast,
    t4_jaccard = t4$jaccard_mean,
    t4_pc1_cor = t4$pc1_cor,
    t5_anchor_retention = t5$anchor_retention,
    t6_boundary_shift = t6$boundary_shift,
    t7_separation = t7$separation,
    t8_fold = t8$result
  )

  for (si in seq_along(dt_result$steps)) {
    all_decision_rows[[length(all_decision_rows)+1]] <- data.table(
      chrom = chr, interval_id = iid, step = si, reasoning = dt_result$steps[[si]]
    )
  }
}

# =============================================================================
# WRITE
# =============================================================================

test_dt <- if (length(all_test_rows) > 0) rbindlist(all_test_rows, fill = TRUE) else data.table()
verd_dt <- if (length(all_verdict_rows) > 0) rbindlist(all_verdict_rows) else data.table()
deci_dt <- if (length(all_decision_rows) > 0) rbindlist(all_decision_rows) else data.table()

fwrite(test_dt, file.path(outdir, "hypothesis_test_results.tsv.gz"), sep = "\t")
fwrite(verd_dt, file.path(outdir, "hypothesis_verdicts.tsv"), sep = "\t")
fwrite(deci_dt, file.path(outdir, "hypothesis_decision_tree.tsv"), sep = "\t")

# Summary
message("\n[C01f] === VERDICT SUMMARY ===")
if (nrow(verd_dt) > 0) {
  for (v in sort(unique(verd_dt$verdict))) {
    n <- sum(verd_dt$verdict == v)
    message("  ", v, ": ", n, " candidates")
  }
}

# Compact hypothesis table
message("\n[C01f] === HYPOTHESIS x PREDICTION x TEST TABLE ===")
hyp_table <- data.table(
  hypothesis = c("H1: Family structure", "H2: Broad inversion + hot core",
                  "H3: Nested/composite", "H4: Sub-haplotype/founder", "H5: Technical"),
  prediction = c("High within-band relatedness, low ancestry diversity",
                  "Same samples inner/outer, stronger coherence inner",
                  "Changed composition at inner/outer boundary",
                  "Persists after kin pruning, carrier subclusters",
                  "Inner more informative, no composition change"),
  test = c("T1 relatedness ratio + T2 eff_K",
           "T4 inner/outer Jaccard + T5 anchor stability",
           "T6 regime change + T4 composition shift",
           "T3 kin-pruned retention + T7 carrier-only PCA",
           "T8 eigenvalue ratio fold"),
  key_metric = c("T1 ratio > 3 + T2 eff_K < 1.5",
                  "T4 Jaccard > 0.7 + same anchors",
                  "T4 Jaccard < 0.5 + regime shift",
                  "T3 signal retained + carrier subgroups",
                  "T8 fold > 3 + no composition change")
)
fwrite(hyp_table, file.path(outdir, "hypothesis_prediction_table.tsv"), sep = "\t")

message("\n[DONE] -> ", outdir)
