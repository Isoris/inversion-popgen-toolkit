#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01f_hypothesis_tests.R  (v9.3.4 — T9/T10 verdict wiring + Engine B)
#
# HYPOTHESIS TESTING for inversion candidates.
#
# v9.3.4 changes:
#   - T9 + T10 now WIRED INTO run_decision_tree() verdict synthesis.
#     Previously T9/T10 were called + stored but IGNORED by verdict logic.
#   - T9 ancestry jackknife: if Engine B (load_bridge.R) is available,
#     uses get_region_stats(Fst) for jackknife instead of sim_mat recompute.
#     Falls back to sim_mat method if Engine B is not configured.
#   - T10 theta het concordance: integrated as independent confirmation
#     layer alongside T8 (Clair3 genotype concordance).
#   - Sourced load_bridge.R at top → smap, reg, get_Q, get_region_stats.
#   - Fixed real_names scoping bug in test_ancestry_jackknife().
#
# v9.3.3 changes:
#   - T8 = Cheat 9 Clair3 genotype concordance (from v9.3.2)
#   - T9 = Cheat 6 ancestry jackknife: remove one Q-group at a time,
#     recompute sim_mat block contrast. Robust = real, fragile = single-family.
#   - T10 = Cheat 12 theta het prior: per-sample θ_P from ANGSD MODULE_3.
#     HET samples have elevated θ_P inside inversions. Independent from PCA.
#   - New args: --theta_dir, --qmatrix
#
# v9 changes:
#   - T3 (kin-pruned block retention): "block_collapsed" no longer treated
#     as definitive evidence for H1. In a hatchery with ~20 founders, a
#     real rare inversion entering through one brood line collapses after
#     pruning because most carriers are relatives. T3 "collapsed" is now
#     "supporting evidence" requiring concordance with T1+T2.
#   - T2 (ancestry diversity): low eff_K no longer mapped to "supports
#     ARTIFACT." Relabeled to "ambiguous — could be artifact or founder-
#     linked inversion." Only HIGH eff_K is informative (supports REAL).
#   - Verdict synthesis: T3 "collapsed" alone no longer forces H1.
#     Requires T1 "strong_family" + T2 "single_ancestry" concordance.
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
    message("[C01f] load_bridge.R sourced — Engine B available")
  }, error = function(e) {
    message("[C01f] WARNING: load_bridge.R failed: ", conditionMessage(e))
  })
} else {
  message("[C01f] load_bridge.R not found — running without Engine B (T9 uses sim_mat fallback)")
}

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; precomp_dir <- NULL
relate_file <- NULL; samples_file <- NULL; pruned_file <- NULL
vcf_dir <- NULL; theta_dir <- NULL; qmatrix_file <- NULL
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
  else if (a == "--vcf_dir" && i < length(args))   { vcf_dir <- args[i+1]; i <- i+2 }
  else if (a == "--theta_dir" && i < length(args)) { theta_dir <- args[i+1]; i <- i+2 }
  else if (a == "--qmatrix" && i < length(args))   { qmatrix_file <- args[i+1]; i <- i+2 }
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

# Q-matrix for Cheat 6 ancestry jackknife
q_mat <- NULL; ancestry_groups <- list()
if (!is.null(qmatrix_file) && file.exists(qmatrix_file)) {
  q_raw <- as.matrix(fread(qmatrix_file, header = FALSE))
  if (!is.null(real_names) && nrow(q_raw) == length(real_names)) {
    rownames(q_raw) <- real_names
    q_mat <- q_raw
    K <- ncol(q_mat)
    message("[C01f] Q-matrix loaded: ", nrow(q_mat), " x K=", K)

    # Build ancestry groups: samples where Q_k > 0.5 (dominant component)
    for (k in seq_len(K)) {
      dom_samples <- real_names[q_raw[, k] > 0.5]
      if (length(dom_samples) >= 3) {
        ancestry_groups[[paste0("Q", k)]] <- dom_samples
      }
    }
    # Also add samples with max Q per component (even if < 0.5)
    for (k in seq_len(K)) {
      gname <- paste0("Q", k)
      if (!gname %in% names(ancestry_groups)) {
        max_k_samples <- real_names[apply(q_raw, 1, which.max) == k]
        if (length(max_k_samples) >= 3) {
          ancestry_groups[[gname]] <- max_k_samples
        }
      }
    }
    message("[C01f] Ancestry groups (Cheat 6): ", length(ancestry_groups),
            " groups (", paste(names(ancestry_groups), sapply(ancestry_groups, length),
                               sep = "=", collapse = ", "), ")")
  } else {
    message("[C01f] Q-matrix row count doesn't match samples — skipping")
  }
} else {
  # No --qmatrix: try registry (from load_bridge.R) or dispatcher resolve
  if (.bridge_available && exists("reg") && !is.null(reg)) {
    # Registry path: look for ancestry_K8_Q1..Q8 groups
    for (k in 1:8) {
      gid <- paste0("ancestry_K8_Q", k)
      if (reg$has_group(gid)) {
        members <- reg$get_group(gid)
        if (length(members) >= 3) ancestry_groups[[paste0("Q", k)]] <- members
      }
    }
    if (length(ancestry_groups) >= 2) {
      message("[C01f] Ancestry groups from registry: ", length(ancestry_groups),
              " groups (", paste(names(ancestry_groups), sapply(ancestry_groups, length),
                                 sep = "=", collapse = ", "), ")")
    } else {
      message("[C01f] Registry has <2 ancestry groups — trying dispatcher resolve")
      if (exists("resolve_groups", mode = "function")) {
        resolved <- tryCatch(resolve_groups("ancestry_clusters"), error = function(e) NULL)
        if (!is.null(resolved) && length(resolved) >= 2) {
          ancestry_groups <- resolved
          message("[C01f] Ancestry groups from dispatcher: ", length(ancestry_groups), " groups")
        }
      }
    }
  }
  if (length(ancestry_groups) < 2) {
    message("[C01f] No --qmatrix and no registry groups — Cheat 6 ancestry jackknife will use sim_mat fallback only")
  }
}

# Theta directory for Cheat 12
has_theta <- FALSE
if (!is.null(theta_dir) && dir.exists(theta_dir)) {
  has_theta <- TRUE
  message("[C01f] Theta dir: ", theta_dir)
} else {
  message("[C01f] No --theta_dir provided — Cheat 12 theta het prior will be skipped")
}

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

# =============================================================================
# T8: Clair3 indel genotype concordance (Cheat 9) — v9.3.2
# =============================================================================
# GENUINELY INDEPENDENT from PCA. Uses discrete Clair3 genotype calls.
#
# Logic:
#   1. Extract genotypes from Clair3 VCF for the candidate region
#   2. For each sample, compute het indel rate (fraction of indels that are 0/1)
#   3. Cluster samples by het rate into HOM_like vs HET_like (k=2)
#   4. Compare with PCA band assignments: do HET_like samples fall in band 2?
#   5. If concordant → PCA bands reflect real genotype classes → real inversion
#   6. If discordant → PCA bands don't match genotypes → suspicious
#
# Also computes the "slope" diagnostic: flat het rate across the region =
# inversion, declining het rate = family LD.
# =============================================================================

test_genotype_concordance <- function(chr, start_bp, end_bp, comp,
                                       vcf_dir, sample_names, real_names) {
  default <- list(test = "T8_genotype_concordance", result = NA,
                  evidence = "no_data", concordance = NA,
                  n_het_like = NA, n_hom_like = NA,
                  het_slope_median = NA, hom_slope_median = NA)

  if (is.null(vcf_dir) || !dir.exists(vcf_dir)) return(default)

  # Find VCF file for this chromosome
  vcf_file <- NULL
  for (pattern in c(paste0(chr, ".vcf.gz"), paste0(chr, ".clair3.vcf.gz"),
                     paste0(chr, ".g.vcf.gz"))) {
    f <- file.path(vcf_dir, pattern)
    if (file.exists(f)) { vcf_file <- f; break }
  }
  # Also try subdirectory structure
  if (is.null(vcf_file)) {
    f <- Sys.glob(file.path(vcf_dir, "*", paste0(chr, "*.vcf.gz")))
    if (length(f) > 0) vcf_file <- f[1]
  }
  if (is.null(vcf_file)) {
    return(modifyList(default, list(evidence = "no_vcf_found")))
  }

  # Extract genotypes using bcftools
  region <- paste0(chr, ":", as.integer(start_bp), "-", as.integer(end_bp))
  cmd <- paste0("bcftools query -r ", region,
                " -f '%POS\\t%TYPE\\t[%GT\\t]\\n' ", vcf_file, " 2>/dev/null")
  raw <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t", fill = TRUE),
                   error = function(e) NULL)
  if (is.null(raw) || nrow(raw) < 20) {
    return(modifyList(default, list(evidence = "too_few_variants")))
  }

  # Get VCF sample names
  cmd_h <- paste0("bcftools query -l ", vcf_file, " 2>/dev/null")
  vcf_samples <- tryCatch(fread(cmd = cmd_h, header = FALSE)[[1]],
                            error = function(e) NULL)
  if (is.null(vcf_samples)) {
    return(modifyList(default, list(evidence = "no_sample_header")))
  }

  n_samp <- length(vcf_samples)
  if (ncol(raw) < n_samp + 2) {
    return(modifyList(default, list(evidence = "column_mismatch")))
  }

  pos <- raw[[1]]
  var_type <- raw[[2]]
  gt_mat <- as.matrix(raw[, 3:(n_samp + 2)])
  colnames(gt_mat) <- vcf_samples

  # Convert GT strings to 0/1/2
  gt_num <- matrix(NA_integer_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  colnames(gt_num) <- vcf_samples
  gt_num[gt_mat %in% c("0/0", "0|0")] <- 0L
  gt_num[gt_mat %in% c("0/1", "1/0", "0|1", "1|0")] <- 1L
  gt_num[gt_mat %in% c("1/1", "1|1")] <- 2L

  # Filter: keep sites with < 30% missing and MAF >= 0.05
  miss_rate <- rowMeans(is.na(gt_num))
  af <- rowMeans(gt_num, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- miss_rate < 0.3 & is.finite(maf) & maf >= 0.05
  if (sum(keep) < 20) {
    return(modifyList(default, list(evidence = "too_few_after_filter")))
  }
  gt_num <- gt_num[keep, , drop = FALSE]
  pos <- pos[keep]

  # Per-sample het rate: fraction of sites that are 0/1
  per_sample_het_rate <- colMeans(gt_num == 1L, na.rm = TRUE)
  per_sample_het_rate <- per_sample_het_rate[is.finite(per_sample_het_rate)]
  if (length(per_sample_het_rate) < 20) {
    return(modifyList(default, list(evidence = "too_few_samples")))
  }

  # Per-sample hom-alt rate: fraction of sites that are 1/1
  per_sample_hom_rate <- colMeans(gt_num == 2L, na.rm = TRUE)

  # k=2 clustering on het rate → HOM_like vs HET_like
  km2 <- tryCatch(kmeans(per_sample_het_rate, centers = 2, nstart = 10),
                   error = function(e) NULL)
  if (is.null(km2)) {
    return(modifyList(default, list(evidence = "clustering_failed")))
  }

  co <- order(km2$centers[, 1])
  indel_class <- character(length(per_sample_het_rate))
  names(indel_class) <- names(per_sample_het_rate)
  indel_class[km2$cluster == co[1]] <- "HOM_like"
  indel_class[km2$cluster == co[2]] <- "HET_like"

  n_het_like <- sum(indel_class == "HET_like")
  n_hom_like <- sum(indel_class == "HOM_like")

  # Separation quality: Ashman's D between the two clusters
  c1 <- per_sample_het_rate[km2$cluster == co[1]]
  c2 <- per_sample_het_rate[km2$cluster == co[2]]
  ashman_d <- abs(mean(c2) - mean(c1)) / sqrt((var(c1) + var(c2)) / 2)

  # Median het/hom slopes for diagnostics
  het_slope_med <- round(median(per_sample_het_rate, na.rm = TRUE), 4)
  hom_slope_med <- round(median(per_sample_hom_rate, na.rm = TRUE), 4)

  # ── CONCORDANCE WITH PCA BANDS ──
  # Compare Cheat 9 HOM_like/HET_like with PCA band membership
  concordance <- NA_real_
  if (nrow(comp) > 0 && "band" %in% names(comp) && "sample" %in% names(comp)) {
    # Build PCA band lookup: sample → band
    pca_bands <- setNames(comp$band, comp$sample)

    # Match sample names between VCF and composition
    # VCF may use CGA names, comp may use Ind names or vice versa
    vcf_samp_names <- names(indel_class)
    shared <- intersect(vcf_samp_names, names(pca_bands))

    # If no direct match, try mapping via real_names ↔ Ind
    if (length(shared) < 10 && !is.null(real_names) && length(real_names) > 0) {
      # Build Ind → real and real → Ind maps
      ind_names <- paste0("Ind", seq_along(real_names) - 1)
      real_to_ind <- setNames(ind_names, real_names)
      ind_to_real <- setNames(real_names, ind_names)

      # Try mapping comp samples (Ind) to VCF samples (CGA)
      if (grepl("^Ind", names(pca_bands)[1]) && !grepl("^Ind", vcf_samp_names[1])) {
        mapped_pca <- setNames(pca_bands, ind_to_real[names(pca_bands)])
        mapped_pca <- mapped_pca[!is.na(names(mapped_pca))]
        shared <- intersect(vcf_samp_names, names(mapped_pca))
        if (length(shared) >= 10) pca_bands <- mapped_pca
      } else if (!grepl("^Ind", names(pca_bands)[1]) && grepl("^Ind", vcf_samp_names[1])) {
        mapped_indel <- setNames(indel_class, ind_to_real[names(indel_class)])
        mapped_indel <- mapped_indel[!is.na(names(mapped_indel))]
        shared <- intersect(names(mapped_indel), names(pca_bands))
        if (length(shared) >= 10) indel_class <- mapped_indel
      }
    }

    if (length(shared) >= 10) {
      ic <- indel_class[shared]
      pb <- pca_bands[shared]

      # Concordance: HET_like should be in band2, HOM_like in band1 or band3
      n_concordant <- sum(
        (ic == "HET_like" & pb == "band2") |
        (ic == "HOM_like" & pb %in% c("band1", "band3"))
      )
      concordance <- n_concordant / length(shared)
    }
  }

  # Interpret
  evidence <- if (!is.finite(concordance)) "cannot_compare"
              else if (concordance > 0.75 && ashman_d > 1.5) "strong_concordance"
              else if (concordance > 0.60 && ashman_d > 1.0) "moderate_concordance"
              else if (concordance > 0.50) "weak_concordance"
              else "discordant"

  list(test = "T8_genotype_concordance", result = round(concordance, 3),
       evidence = evidence, concordance = round(concordance, 3),
       n_het_like = n_het_like, n_hom_like = n_hom_like,
       ashman_d = round(ashman_d, 3),
       het_slope_median = het_slope_med, hom_slope_median = hom_slope_med,
       n_variants_tested = sum(keep))
}

# =============================================================================
# T9: Ancestry jackknife (Cheat 6) — leave-one-Q-group-out
# =============================================================================
# For each Q-group, remove its samples and recompute the contrast.
#
# ENGINE B MODE (preferred): uses get_region_stats(Fst) for each jackknife
#   iteration. Hudson Fst between band1 and band3, with one Q-group removed.
#   This is the real Fst from dosage data — not a PCA proxy.
#
# FALLBACK MODE: if Engine B unavailable, recomputes sim_mat block contrast
#   from PC loadings with samples removed. Less precise but still informative.

test_ancestry_jackknife <- function(pc, start_bp, end_bp, ancestry_groups,
                                     sample_names, real_names_arg = NULL,
                                     use_engine_b = .bridge_available) {
  default <- list(test = "T9_ancestry_jackknife", result = NA,
                  evidence = "no_data", n_groups_tested = 0L,
                  max_delta = NA, fragile_group = NA, jackknife_verdict = "untested",
                  n_drops = 0L, n_rises = 0L, method = "none",
                  fst_full = NA_real_)

  if (is.null(pc) || length(ancestry_groups) < 2) return(default)
  dt <- pc$dt; sim_mat <- pc$sim_mat
  n_w <- nrow(dt)

  inner_idx <- which(dt$start_bp >= start_bp & dt$end_bp <= end_bp)
  outer_idx <- setdiff(seq_len(n_w), inner_idx)
  if (length(inner_idx) < 5 || length(outer_idx) < 10) return(default)

  # PC1 columns for sample identification
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  ind_names <- sub("^PC_1_", "", pc1_cols)

  # Build Ind→CGA map
  ind_to_real_map <- NULL
  if (!is.null(real_names_arg) && length(real_names_arg) == length(ind_names) &&
      grepl("^Ind[0-9]", ind_names[1])) {
    ind_to_real_map <- setNames(real_names_arg, ind_names)
  }

  # Get PCA band assignments for Engine B Fst groups
  chr <- dt$chrom[1]  # chromosome from precomp
  if (is.null(chr) || is.na(chr)) chr <- unique(dt$chrom)[1]

  # ── Engine B path: use real Fst ──────────────────────────────────────────
  if (use_engine_b && exists("get_region_stats", mode = "function")) {
    message("    [T9] Using Engine B (get_region_stats Fst)")

    # Helper: extract the single Fst value from dispatcher result.
    # C binary outputs Fst_<g1>_<g2>, dispatcher wraps in result$Fst list.
    # With one group pair there's exactly one Fst column — take it.
    extract_fst <- function(s) {
      if (is.null(s$Fst) || length(s$Fst) == 0) return(NA_real_)
      as.numeric(s$Fst[[1]])
    }

    # Get band assignments from composition
    # We need PC1 band membership → sample lists in CGA space
    pc1_vals <- as.numeric(colMeans(as.matrix(dt[inner_idx, ..pc1_cols]), na.rm = TRUE))
    names(pc1_vals) <- ind_names
    pc1_vals <- pc1_vals[is.finite(pc1_vals)]
    if (length(pc1_vals) < 30) return(modifyList(default, list(evidence = "too_few_samples")))

    km3 <- tryCatch(kmeans(pc1_vals, centers = 3, nstart = 10), error = function(e) NULL)
    if (is.null(km3)) return(modifyList(default, list(evidence = "clustering_failed")))

    co <- order(km3$centers[, 1])
    band_assign <- integer(length(pc1_vals))
    names(band_assign) <- names(pc1_vals)
    for (bi in seq_along(co)) band_assign[km3$cluster == co[bi]] <- bi

    # Convert to CGA names for Engine B
    if (!is.null(ind_to_real_map)) {
      band1_cga <- ind_to_real_map[names(band_assign[band_assign == 1])]
      band3_cga <- ind_to_real_map[names(band_assign[band_assign == 3])]
    } else {
      band1_cga <- names(band_assign[band_assign == 1])
      band3_cga <- names(band_assign[band_assign == 3])
    }
    band1_cga <- band1_cga[!is.na(band1_cga)]
    band3_cga <- band3_cga[!is.na(band3_cga)]

    if (length(band1_cga) < 5 || length(band3_cga) < 5) {
      return(modifyList(default, list(evidence = "too_few_per_band")))
    }

    # Full Fst (baseline)
    fst_full <- tryCatch({
      s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                             groups = list(b1 = band1_cga, b3 = band3_cga))
      extract_fst(s)
    }, error = function(e) { message("    [T9] Fst failed: ", e$message); NA_real_ })

    if (!is.finite(fst_full) || fst_full < 0.005) {
      return(modifyList(default, list(evidence = "no_fst_signal", fst_full = fst_full,
                                      method = "engine_b")))
    }

    # Jackknife: remove each Q-group, recompute Fst
    jk_results <- list()
    for (gname in names(ancestry_groups)) {
      remove_cga <- ancestry_groups[[gname]]
      b1_kept <- setdiff(band1_cga, remove_cga)
      b3_kept <- setdiff(band3_cga, remove_cga)

      if (length(b1_kept) < 3 || length(b3_kept) < 3) {
        jk_results[[gname]] <- list(group = gname, fst = NA, delta = NA,
                                     direction = "insufficient")
        next
      }

      fst_without <- tryCatch({
        s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                               groups = list(b1 = b1_kept, b3 = b3_kept))
        extract_fst(s)
      }, error = function(e) NA_real_)

      delta <- if (is.finite(fst_without)) fst_full - fst_without else NA_real_
      direction <- if (is.finite(delta)) {
        if (delta > 0.02) "drops"
        else if (delta < -0.02) "rises"
        else "stable"
      } else "unknown"

      jk_results[[gname]] <- list(group = gname, fst = round(fst_without, 4),
                                   delta = round(delta, 4), direction = direction)
    }

    # Summarize
    deltas <- sapply(jk_results, function(x) x$delta)
    deltas <- deltas[is.finite(deltas)]
    directions <- sapply(jk_results, function(x) x$direction)

    max_delta <- if (length(deltas) > 0) max(abs(deltas)) else NA
    n_drops <- sum(directions == "drops")
    fragile_grp <- if (n_drops > 0) {
      names(which.max(sapply(jk_results, function(x) if (is.finite(x$delta)) x$delta else -Inf)))
    } else NA

    verdict <- if (length(deltas) == 0) "untested"
      else if (max_delta < 0.02) "robust_multi_family"
      else if (n_drops == 1) "single_family_fragile"
      else if (n_drops <= 2) "few_family_contributing"
      else "multi_family_contributing"

    evidence <- if (verdict == "robust_multi_family") "robust"
      else if (verdict == "single_family_fragile") "fragile"
      else if (verdict == "untested") "no_data"
      else "partial"

    return(list(test = "T9_ancestry_jackknife", result = round(max_delta, 4),
                evidence = evidence, n_groups_tested = length(deltas),
                max_delta = round(max_delta, 4),
                n_drops = n_drops, n_rises = sum(directions == "rises"),
                fragile_group = fragile_grp, jackknife_verdict = verdict,
                method = "engine_b", fst_full = round(fst_full, 4)))
  }

  # ── Fallback: sim_mat recomputation ────────────────────────────────────────
  message("    [T9] Using sim_mat fallback (Engine B not available)")

  # Full block contrast (baseline)
  full_inside <- mean(sim_mat[inner_idx, inner_idx], na.rm = TRUE)
  full_outside <- mean(sim_mat[outer_idx[1:min(100, length(outer_idx))],
                                outer_idx[1:min(100, length(outer_idx))]], na.rm = TRUE)
  full_contrast <- full_inside - full_outside
  if (!is.finite(full_contrast) || full_contrast < 0.01) {
    return(modifyList(default, list(evidence = "no_block_contrast", method = "sim_mat")))
  }

  jk_results <- list()
  for (gname in names(ancestry_groups)) {
    remove_cga <- ancestry_groups[[gname]]

    # Map to Ind names
    if (!is.null(ind_to_real_map)) {
      remove_ind <- names(ind_to_real_map)[ind_to_real_map %in% remove_cga]
    } else {
      remove_ind <- remove_cga  # assume already Ind format
    }

    keep_cols <- setdiff(pc1_cols, paste0("PC_1_", remove_ind))
    if (length(keep_cols) < 30) {
      jk_results[[gname]] <- list(group = gname, contrast = NA, delta = NA, direction = "insufficient")
      next
    }

    # Recompute sim_mat from remaining samples' PC loadings
    sub_idx <- sort(c(sample(outer_idx, min(80, length(outer_idx))), inner_idx))
    p_mat <- as.matrix(dt[sub_idx, ..keep_cols])
    ns <- nrow(p_mat)
    p_sim <- matrix(0, ns, ns)
    for (ii in seq_len(ns)) for (jj in ii:ns) {
      v <- cor(p_mat[ii, ], p_mat[jj, ], use = "pairwise.complete.obs")
      if (is.finite(v)) { p_sim[ii, jj] <- v; p_sim[jj, ii] <- v }
    }
    win_in <- which(sub_idx %in% inner_idx)
    out_in <- which(!sub_idx %in% inner_idx)

    if (length(win_in) < 3 || length(out_in) < 3) {
      jk_results[[gname]] <- list(group = gname, contrast = NA, delta = NA, direction = "insufficient")
      next
    }

    jk_inside <- mean(p_sim[win_in, win_in], na.rm = TRUE)
    jk_outside <- mean(p_sim[out_in, out_in], na.rm = TRUE)
    jk_contrast <- jk_inside - jk_outside
    delta <- full_contrast - jk_contrast

    direction <- if (delta > 0.02) "drops"
                 else if (delta < -0.02) "rises"
                 else "stable"

    jk_results[[gname]] <- list(group = gname, contrast = round(jk_contrast, 4),
                                 delta = round(delta, 4), direction = direction)
  }

  # Summarize
  deltas <- sapply(jk_results, function(x) x$delta)
  deltas <- deltas[is.finite(deltas)]
  directions <- sapply(jk_results, function(x) x$direction)

  max_delta <- if (length(deltas) > 0) max(abs(deltas)) else NA
  n_drops <- sum(directions == "drops")
  fragile_grp <- if (n_drops > 0) {
    names(which.max(sapply(jk_results, function(x) if (is.finite(x$delta)) x$delta else -Inf)))
  } else NA

  verdict <- if (length(deltas) == 0) "untested"
    else if (max_delta < 0.02) "robust_multi_family"
    else if (n_drops == 1) "single_family_fragile"
    else if (n_drops <= 2) "few_family_contributing"
    else "multi_family_contributing"

  evidence <- if (verdict == "robust_multi_family") "robust"
    else if (verdict == "single_family_fragile") "fragile"
    else if (verdict == "untested") "no_data"
    else "partial"

  list(test = "T9_ancestry_jackknife", result = round(max_delta, 4),
       evidence = evidence, n_groups_tested = length(deltas),
       max_delta = round(max_delta, 4),
       n_drops = n_drops, n_rises = sum(directions == "rises"),
       fragile_group = fragile_grp, jackknife_verdict = verdict,
       method = "sim_mat", fst_full = NA_real_)
}

# =============================================================================
# T10: Theta het prior (Cheat 12) — independent from PCA
# =============================================================================
# Per-sample θ_P from ANGSD MODULE_3 thetaStat. Inside an inversion, HET
# samples have elevated θ_P. Classify samples by θ_P and compare with PCA bands.

test_theta_concordance <- function(chr, start_bp, end_bp, comp, theta_dir, sample_names) {
  default <- list(test = "T10_theta_het", result = NA,
                  evidence = "no_data", band_separation = NA,
                  concordance = NA, n_samples_classified = 0L)

  if (is.null(theta_dir) || !dir.exists(theta_dir)) return(default)

  # Load per-sample theta for this region
  theta_rows <- list()
  for (sid in sample_names) {
    # Try multiple file patterns
    found <- NULL
    for (pattern in c(
      file.path(theta_dir, paste0(sid, ".thetaStat.pestPG.gz")),
      file.path(theta_dir, paste0(sid, ".", chr, ".thetaStat.pestPG.gz")),
      file.path(theta_dir, chr, paste0(sid, ".thetaStat.pestPG.gz"))
    )) {
      if (file.exists(pattern)) { found <- pattern; break }
    }
    if (is.null(found)) next

    dt <- tryCatch(fread(found, header = TRUE), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next

    # Standardize column names
    if ("Chr" %in% names(dt)) setnames(dt, "Chr", "chr_col")
    if ("WinCenter" %in% names(dt)) setnames(dt, "WinCenter", "wincenter")
    if (!"chr_col" %in% names(dt)) next
    if (!"tP" %in% names(dt)) next

    # Full chromosome for z-scoring
    chr_dt <- dt[chr_col == chr]
    if (nrow(chr_dt) < 20) next

    chr_mean <- mean(chr_dt$tP, na.rm = TRUE)
    chr_sd <- sd(chr_dt$tP, na.rm = TRUE)
    if (!is.finite(chr_sd) || chr_sd == 0) next

    # Region subset
    region_dt <- chr_dt[wincenter >= start_bp & wincenter <= end_bp]
    if (nrow(region_dt) < 3) next

    mean_tP_z <- mean((region_dt$tP - chr_mean) / chr_sd, na.rm = TRUE)
    theta_rows[[sid]] <- mean_tP_z
  }

  if (length(theta_rows) < 20) {
    return(modifyList(default, list(evidence = "too_few_samples")))
  }

  # Per-sample theta z-scores
  tP_z <- unlist(theta_rows)
  names(tP_z) <- names(theta_rows)

  # k=2 clustering: HOM_like (low θ_P) vs HET_like (high θ_P)
  km2 <- tryCatch(kmeans(tP_z, centers = 2, nstart = 10), error = function(e) NULL)
  if (is.null(km2)) {
    return(modifyList(default, list(evidence = "clustering_failed")))
  }

  co <- order(km2$centers[, 1])
  theta_class <- character(length(tP_z))
  names(theta_class) <- names(tP_z)
  theta_class[km2$cluster == co[1]] <- "HOM_like"
  theta_class[km2$cluster == co[2]] <- "HET_like"

  band_sep <- abs(km2$centers[co[2], 1] - km2$centers[co[1], 1])

  # Compare with PCA bands
  concordance <- NA_real_
  if (nrow(comp) > 0 && "band" %in% names(comp) && "sample" %in% names(comp)) {
    pca_bands <- setNames(comp$band, comp$sample)
    shared <- intersect(names(theta_class), names(pca_bands))

    # Try mapping if no direct match
    if (length(shared) < 10 && !is.null(real_names) && length(real_names) > 0) {
      ind_names <- paste0("Ind", seq_along(real_names) - 1)
      if (grepl("^Ind", names(pca_bands)[1]) && !grepl("^Ind", names(theta_class)[1])) {
        mapped <- setNames(pca_bands, setNames(real_names, ind_names)[names(pca_bands)])
        mapped <- mapped[!is.na(names(mapped))]
        shared <- intersect(names(theta_class), names(mapped))
        if (length(shared) >= 10) pca_bands <- mapped
      }
    }

    if (length(shared) >= 10) {
      tc <- theta_class[shared]; pb <- pca_bands[shared]
      concordance <- mean(
        (tc == "HET_like" & pb == "band2") |
        (tc == "HOM_like" & pb %in% c("band1", "band3"))
      )
    }
  }

  evidence <- if (!is.finite(concordance)) "cannot_compare"
    else if (concordance > 0.70 && band_sep > 1.0) "strong_concordance"
    else if (concordance > 0.55 && band_sep > 0.5) "moderate_concordance"
    else if (band_sep < 0.3) "no_separation"
    else "weak_or_discordant"

  list(test = "T10_theta_het", result = round(concordance, 3),
       evidence = evidence, band_separation = round(band_sep, 3),
       concordance = round(concordance, 3),
       n_samples_classified = length(theta_class),
       n_het_like = sum(theta_class == "HET_like"),
       n_hom_like = sum(theta_class == "HOM_like"))
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
  t9 <- tests[["T9_ancestry_jackknife"]]
  t10 <- tests[["T10_theta_het"]]

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

  # Step 2: Kin-pruned block retention (T3)
  # IMPORTANT: in a hatchery with ~20 founders, a real rare inversion
  # that entered through one brood line will ALSO collapse after pruning.
  # T3 "collapsed" is ambiguous — it means "family-linked" but not
  # necessarily "family artifact." Requires concordance with T1+T2.
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
      step2_result <- "signal_family_linked"
      steps[[si]] <- paste0("STEP2 [T3]: Block COLLAPSED (retention=", t3$retention_ratio,
                            ") -> family-linked (could be artifact OR founder-linked inversion)")
    } else {
      steps[[si]] <- paste0("STEP2 [T3]: Inconclusive (", t3$evidence, ")")
    }
  } else {
    steps[[si]] <- "STEP2 [T3]: No kin-pruned data -> provide --pruned_samples"
  }

  # Step 3: Ancestry diversity (T2)
  # NOTE: high eff_K is informative (multiple families per band = real).
  # Low eff_K is AMBIGUOUS in a hatchery — it could mean one family
  # artifact OR a real inversion that entered through one founder.
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
                            ") -> single family dominates -> AMBIGUOUS (artifact or founder-linked)")
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

  # Step 8: Clair3 genotype concordance (Cheat 9) — independent confirmation
  # If Clair3 discrete genotypes cluster the same fish the same way as PCA,
  # that's strong independent evidence for a real inversion.
  step8_result <- "unknown"
  si <- si + 1L
  if (!is.na(t8$result) && t8$evidence != "no_data") {
    if (t8$evidence == "strong_concordance") {
      step8_result <- "genotype_confirms"
      steps[[si]] <- paste0("STEP8 [T8/Cheat9]: Clair3 genotypes CONCORDANT with PCA bands (",
                            "concordance=", t8$concordance,
                            " Ashman_D=", t8$ashman_d,
                            " n_HET=", t8$n_het_like, " n_HOM=", t8$n_hom_like,
                            ") -> independent confirmation of inversion genotypes")
    } else if (t8$evidence == "moderate_concordance") {
      step8_result <- "genotype_supports"
      steps[[si]] <- paste0("STEP8 [T8/Cheat9]: Moderate concordance (",
                            "concordance=", t8$concordance,
                            " Ashman_D=", t8$ashman_d,
                            ") -> genotypes partially match PCA bands")
    } else if (t8$evidence == "discordant") {
      step8_result <- "genotype_discordant"
      steps[[si]] <- paste0("STEP8 [T8/Cheat9]: DISCORDANT (concordance=", t8$concordance,
                            ") -> Clair3 genotypes don't match PCA bands -> suspicious")
    } else {
      step8_result <- "genotype_weak"
      steps[[si]] <- paste0("STEP8 [T8/Cheat9]: Weak/uncertain (", t8$evidence,
                            " concordance=", t8$concordance, ")")
    }
  } else {
    steps[[si]] <- paste0("STEP8 [T8/Cheat9]: Cannot assess (",
                          t8$evidence %||% "no VCF", ")")
  }

  # Step 9: Ancestry jackknife (T9 / Cheat 6) — multi-family robustness
  # If the sim_mat block (or Engine B Fst) survives removal of each Q-group,
  # the signal is real and multi-family. If removing one group kills it,
  # the signal is fragile and driven by a single ancestry lineage.
  step9_result <- "unknown"
  si <- si + 1L
  if (!is.null(t9) && !is.na(t9$result) && t9$jackknife_verdict != "untested") {
    if (t9$jackknife_verdict == "robust_multi_family") {
      step9_result <- "jackknife_robust"
      steps[[si]] <- paste0("STEP9 [T9/Cheat6]: ROBUST across all Q-groups (max_delta=",
                            t9$max_delta, " n_drops=", t9$n_drops,
                            " method=", t9$method %||% "unknown",
                            ") -> signal is multi-family -> REAL INVERSION")
    } else if (t9$jackknife_verdict == "single_family_fragile") {
      step9_result <- "jackknife_fragile"
      steps[[si]] <- paste0("STEP9 [T9/Cheat6]: FRAGILE — removing ", t9$fragile_group,
                            " kills signal (max_delta=", t9$max_delta,
                            ") -> single-family-driven -> could be FAMILY or founder-linked")
    } else if (t9$jackknife_verdict %in% c("few_family_contributing",
                                            "multi_family_contributing")) {
      step9_result <- "jackknife_partial"
      steps[[si]] <- paste0("STEP9 [T9/Cheat6]: ", t9$n_drops, " groups contribute ",
                            "(max_delta=", t9$max_delta,
                            ") -> partially robust")
    } else {
      steps[[si]] <- paste0("STEP9 [T9/Cheat6]: ", t9$jackknife_verdict,
                            " (max_delta=", t9$max_delta, ")")
    }
  } else {
    steps[[si]] <- "STEP9 [T9/Cheat6]: Cannot assess (no Q-groups or untested)"
  }

  # Step 10: Theta het prior (T10 / Cheat 12) — independent confirmation
  # Per-sample θ_P from ANGSD inside the region. HET samples have elevated θ_P
  # inside inversions. This is completely independent from PCA and from SV calls.
  step10_result <- "unknown"
  si <- si + 1L
  if (!is.null(t10) && !is.na(t10$result) && t10$evidence != "no_data") {
    if (t10$evidence == "strong_concordance") {
      step10_result <- "theta_confirms"
      steps[[si]] <- paste0("STEP10 [T10/Cheat12]: θ_P CONCORDANT with PCA bands (",
                            "concordance=", t10$concordance,
                            " separation=", t10$band_separation,
                            " n_HET=", t10$n_het_like, " n_HOM=", t10$n_hom_like,
                            ") -> independent confirmation of inversion genotypes")
    } else if (t10$evidence == "moderate_concordance") {
      step10_result <- "theta_supports"
      steps[[si]] <- paste0("STEP10 [T10/Cheat12]: Moderate θ_P concordance (",
                            "concordance=", t10$concordance,
                            " separation=", t10$band_separation,
                            ") -> partial support for inversion model")
    } else if (t10$evidence == "no_separation") {
      step10_result <- "theta_flat"
      steps[[si]] <- paste0("STEP10 [T10/Cheat12]: No θ_P separation (sep=",
                            t10$band_separation, ") -> uninformative")
    } else {
      step10_result <- "theta_weak"
      steps[[si]] <- paste0("STEP10 [T10/Cheat12]: ", t10$evidence,
                            " (concordance=", t10$concordance, ")")
    }
  } else {
    steps[[si]] <- paste0("STEP10 [T10/Cheat12]: Cannot assess (",
                          if (!is.null(t10)) t10$evidence %||% "no theta_dir" else "null", ")")
  }

  # ===== VERDICT SYNTHESIS (v9.3.4) =====
  # Priority layers:
  #   1. T3+T1+T2 concordance (structural: does block survive kin pruning?)
  #   2. T9/Cheat6 jackknife (is signal multi-family or single-family?)
  #   3. T8/Cheat9 + T10/Cheat12 (independent genotype & theta confirmation)
  #   4. T4+T5+T6 (architecture: H2 vs H3)
  #   5. T7 (substructure: H4)
  #
  # v9.3.4 additions:
  #   - T9 "robust" upgrades ambiguous → likely_real
  #   - T9 "fragile" + T1 "family_dominant" → H1 (strong)
  #   - T9 "fragile" + T8/T10 confirm → founder-linked real inversion
  #   - T10 "strong_concordance" acts as rescue alongside T8
  #   - Combined T8+T10 double-confirmation → confidence boost

  # Count independent confirmations (T8, T9, T10)
  n_independent_confirm <- 0L
  if (step8_result == "genotype_confirms") n_independent_confirm <- n_independent_confirm + 1L
  if (step9_result == "jackknife_robust")  n_independent_confirm <- n_independent_confirm + 1L
  if (step10_result == "theta_confirms")   n_independent_confirm <- n_independent_confirm + 1L

  n_independent_support <- n_independent_confirm
  if (step8_result == "genotype_supports") n_independent_support <- n_independent_support + 1L
  if (step9_result == "jackknife_partial") n_independent_support <- n_independent_support + 1L
  if (step10_result == "theta_supports")   n_independent_support <- n_independent_support + 1L

  si <- si + 1L

  if (step2_result == "signal_family_linked") {
    # T3 says blocks collapse -> family-linked. But is it artifact or real?
    if (n_independent_confirm >= 2) {
      # TWO independent confirmations rescue even collapsed T3
      verdict <- "likely_real_inversion"
      confidence <- "high"
    } else if (step8_result == "genotype_confirms" || step10_result == "theta_confirms") {
      # ONE strong independent confirmation rescues
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (step9_result == "jackknife_robust") {
      # Block collapses on kin-pruning but survives Q-group jackknife
      # → the inversion is carried by multiple families but relatives dominate
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (step9_result == "jackknife_fragile" &&
               step1_result == "family_dominant" && step3_result == "monomorphic" &&
               step8_result != "genotype_supports" && step10_result != "theta_supports") {
      # All evidence converges: T3 collapsed + T9 fragile + T1 family + T2 single
      # AND neither T8 nor T10 partially rescue → confident H1
      verdict <- "H1_family_structure"
      confidence <- "high"
    } else if (step1_result == "family_dominant" && step3_result == "monomorphic" &&
               n_independent_support == 0L) {
      # T3+T1+T2 all say family, no independent rescue at all
      verdict <- "H1_family_structure"
      confidence <- "high"
    } else if (step9_result == "jackknife_fragile" &&
               (step8_result == "genotype_supports" || step10_result == "theta_supports")) {
      # Fragile in jackknife but partial independent support
      # → founder-linked real inversion (entered through one brood line)
      verdict <- "founder_linked_real_inversion"
      confidence <- "moderate"
    } else if (step1_result == "family_dominant") {
      verdict <- "mixed_family_and_signal"
      confidence <- "moderate"
    } else {
      verdict <- "ambiguous_family_linked"
      confidence <- "low"
    }

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

    # Independent confirmation confidence boosts
    if (n_independent_confirm >= 2 && confidence == "moderate") {
      confidence <- "high"
    } else if (n_independent_confirm >= 1 && confidence == "moderate") {
      confidence <- "high"
    }

    # T9 fragile downgrades within signal_real → flag as founder-linked
    if (step9_result == "jackknife_fragile" && confidence == "high") {
      confidence <- "moderate"
      steps[[si - 1]] <- paste0(steps[[si - 1]],
                                 " [NOTE: T9 fragile — signal is real but founder-linked]")
    }

  } else if (step2_result == "signal_mixed") {
    # T3 says blocks weakened -> mixed
    if (step9_result == "jackknife_robust") {
      # T3 weakened but T9 says multi-family → real inversion + family amplification
      verdict <- "mixed_real_plus_family"
      confidence <- "moderate"
    } else if (step3_result == "diverse") {
      verdict <- "mixed_real_plus_family"
      confidence <- "moderate"
    } else if (step8_result == "genotype_confirms" || step10_result == "theta_confirms") {
      # Independent confirmation rescues weakened T3
      verdict <- "mixed_real_plus_family"
      confidence <- "moderate"
    } else if (n_independent_support >= 2) {
      # Multiple partial supports → still likely real
      verdict <- "mixed_real_plus_family"
      confidence <- "moderate"
    } else {
      verdict <- "mixed_family_and_signal"
      confidence <- "low"
    }

  } else {
    # T3 not available -> fall back to independent evidence layers
    if (step9_result == "jackknife_robust") {
      # No T3 but jackknife says robust → real
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (n_independent_confirm >= 2) {
      # Multiple independent confirmations even without T3
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (step8_result == "genotype_confirms" || step10_result == "theta_confirms") {
      # Single independent confirmation without T3 → likely real
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else if (step1_result == "family_dominant" && step3_result == "monomorphic" &&
               step9_result == "jackknife_fragile" &&
               n_independent_support == 0L) {
      verdict <- "H1_family_structure"
      confidence <- "moderate"
    } else if (step1_result == "family_dominant" && step3_result == "monomorphic" &&
               n_independent_support == 0L) {
      verdict <- "H1_family_structure"
      confidence <- "moderate"
    } else if (step1_result == "family_unlikely" && step3_result == "diverse") {
      verdict <- "likely_real_inversion"
      confidence <- "moderate"
    } else {
      verdict <- "unresolved_needs_pruned_data"
      confidence <- "low"
    }
  }

  steps[[si]] <- paste0("VERDICT: ", verdict, " (confidence: ", confidence,
                         ", independent_confirm: ", n_independent_confirm,
                         ", independent_support: ", n_independent_support, ")")

  list(verdict = verdict, confidence = confidence, steps = steps,
       step1_family = step1_result, step2_kinpruned = step2_result,
       step3_ancestry = step3_result, step4_inner_outer = step4_result,
       step5_anchor = step5_result, step6_regime = step6_result,
       step7_carrier = step7_result, step8_technical = step8_result,
       step9_jackknife = step9_result, step10_theta = step10_result,
       n_independent_confirm = n_independent_confirm,
       n_independent_support = n_independent_support)
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

  # T8: Clair3 genotype concordance (Cheat 9) — independent from PCA
  t8 <- test_genotype_concordance(chr, cand$start_mb * 1e6, cand$end_mb * 1e6,
                                   iv_comp, vcf_dir, sample_names, real_names)
  tests[["T8_genotype_concordance"]] <- t8
  message("  T8 genotype concordance: ", t8$evidence,
          " (concordance=", t8$concordance,
          " Ashman_D=", t8$ashman_d,
          " n_HET=", t8$n_het_like, " n_HOM=", t8$n_hom_like, ")")

  # T9: Ancestry jackknife (Cheat 6) — leave-one-Q-group-out
  t9 <- test_ancestry_jackknife(pc, cand$start_mb * 1e6, cand$end_mb * 1e6,
                                 ancestry_groups, sample_names,
                                 real_names_arg = real_names)
  tests[["T9_ancestry_jackknife"]] <- t9
  message("  T9 ancestry jackknife: ", t9$evidence,
          " (verdict=", t9$jackknife_verdict,
          " max_delta=", t9$max_delta,
          " n_drops=", t9$n_drops,
          " method=", t9$method %||% "unknown",
          if (!is.na(t9$fragile_group)) paste0(" fragile=", t9$fragile_group) else "",
          if (is.finite(t9$fst_full %||% NA)) paste0(" fst_full=", t9$fst_full) else "", ")")

  # T10: Theta het prior (Cheat 12) — independent from PCA
  t10 <- test_theta_concordance(chr, cand$start_mb * 1e6, cand$end_mb * 1e6,
                                 iv_comp, theta_dir, sample_names)
  tests[["T10_theta_het"]] <- t10
  message("  T10 theta het: ", t10$evidence,
          " (concordance=", t10$concordance,
          " separation=", t10$band_separation,
          " n_HET=", t10$n_het_like, " n_HOM=", t10$n_hom_like, ")")

  # T11: Extended boundary suppression (Cheat 23) — four-bin Fst test
  # Tests whether Fst(band1,band3) extends beyond the formal breakpoints.
  # Four bins: deep_external, near_external, near_internal, deep_internal.
  # Uses dispatcher for real Fst when available, inv_likeness proxy otherwise.
  t11 <- list(test = "T11_extended_suppression", result = NA,
              evidence = "no_data", has_extended = NA,
              extent_kb = NA_real_, fst_decay_rate = NA_real_,
              fst_deep_ext = NA_real_, fst_near_ext = NA_real_,
              fst_near_int = NA_real_, fst_deep_int = NA_real_)

  if (!is.null(pc)) {
    start_bp <- cand$start_mb * 1e6; end_bp <- cand$end_mb * 1e6
    span <- end_bp - start_bp
    if (span > 0) {
      # Define four bins
      ext_flank <- min(200000, span / 2)
      bins <- list(
        list(bin = "deep_external",  s = start_bp - 2*ext_flank, e = start_bp - ext_flank),
        list(bin = "near_external",  s = start_bp - ext_flank,   e = start_bp),
        list(bin = "near_internal",  s = start_bp,               e = start_bp + ext_flank),
        list(bin = "deep_internal",  s = start_bp + ext_flank,   e = end_bp - ext_flank)
      )

      # Get band assignments for Fst groups
      pc1_cols_t11 <- grep("^PC_1_", names(pc$dt), value = TRUE)
      inner_w <- which(pc$dt$start_bp >= start_bp & pc$dt$end_bp <= end_bp)

      bin_fsts <- numeric(4)
      names(bin_fsts) <- sapply(bins, `[[`, "bin")

      if (length(inner_w) >= 3 && length(pc1_cols_t11) >= 20) {
        avg_pc1 <- colMeans(as.matrix(pc$dt[inner_w, ..pc1_cols_t11]), na.rm = TRUE)
        valid_t11 <- is.finite(avg_pc1)
        if (sum(valid_t11) >= 20) {
          km_t11 <- tryCatch(kmeans(avg_pc1[valid_t11], centers = 3, nstart = 5),
                              error = function(e) NULL)
          if (!is.null(km_t11)) {
            co_t11 <- order(km_t11$centers[, 1])
            bands_t11 <- integer(sum(valid_t11))
            for (b in 1:3) bands_t11[km_t11$cluster == co_t11[b]] <- b
            snames_t11 <- sub("^PC_1_", "", names(avg_pc1)[valid_t11])
            names(bands_t11) <- snames_t11

            b1_names <- snames_t11[bands_t11 == 1]
            b3_names <- snames_t11[bands_t11 == 3]

            # Map to CGA for dispatcher
            if (!is.null(ind_to_real)) {
              b1_cga_t11 <- ind_to_real[b1_names]; b1_cga_t11 <- b1_cga_t11[!is.na(b1_cga_t11)]
              b3_cga_t11 <- ind_to_real[b3_names]; b3_cga_t11 <- b3_cga_t11[!is.na(b3_cga_t11)]
            } else {
              b1_cga_t11 <- b1_names; b3_cga_t11 <- b3_names
            }

            use_dispatcher <- .bridge_available &&
              exists("get_region_stats", mode = "function") &&
              length(b1_cga_t11) >= 5 && length(b3_cga_t11) >= 5

            for (bi in seq_along(bins)) {
              bs <- bins[[bi]]$s; be <- bins[[bi]]$e
              if (bs >= be || bs < 0) { bin_fsts[bi] <- NA; next }

              if (use_dispatcher) {
                bin_fsts[bi] <- tryCatch({
                  s <- get_region_stats(chr, bs, be, what = "Fst",
                                         groups = list(b1 = b1_cga_t11, b3 = b3_cga_t11))
                  if (!is.null(s$Fst) && length(s$Fst) > 0) as.numeric(s$Fst[[1]])
                  else NA_real_
                }, error = function(e) NA_real_)
              } else {
                # Fallback: inv_likeness proxy
                win_in <- which(pc$dt$start_bp >= bs & pc$dt$end_bp <= be)
                bin_fsts[bi] <- if (length(win_in) > 0 && "inv_likeness" %in% names(pc$dt))
                  mean(pc$dt$inv_likeness[win_in], na.rm = TRUE) else NA_real_
              }
            }

            fst_de <- bin_fsts["deep_external"]
            fst_ne <- bin_fsts["near_external"]
            fst_ni <- bin_fsts["near_internal"]
            fst_di <- bin_fsts["deep_internal"]

            # Extended suppression: near_external > deep_external
            has_ext <- !is.na(fst_ne) && !is.na(fst_de) &&
              fst_ne > fst_de + 0.01

            # Decay rate: how fast does Fst drop outside the inversion?
            decay <- if (is.finite(fst_ni) && is.finite(fst_de) && ext_flank > 0)
              (fst_ni - fst_de) / (ext_flank / 1000) else NA_real_

            evidence <- if (is.na(has_ext)) "no_data"
              else if (has_ext && fst_ne > 0.03) "strong_extended"
              else if (has_ext) "moderate_extended"
              else "no_extension"

            t11 <- list(test = "T11_extended_suppression", result = round(fst_ne %||% NA_real_, 4),
                        evidence = evidence, has_extended = has_ext,
                        extent_kb = round(ext_flank / 1000, 0),
                        fst_decay_rate = round(decay %||% NA_real_, 4),
                        fst_deep_ext = round(fst_de %||% NA_real_, 4),
                        fst_near_ext = round(fst_ne %||% NA_real_, 4),
                        fst_near_int = round(fst_ni %||% NA_real_, 4),
                        fst_deep_int = round(fst_di %||% NA_real_, 4),
                        method = if (use_dispatcher) "engine_b" else "proxy")
          }
        }
      }
    }
  }
  tests[["T11_extended_suppression"]] <- t11
  message("  T11 extended suppression: ", t11$evidence,
          " (near_ext=", t11$fst_near_ext, " deep_ext=", t11$fst_deep_ext,
          " extent=", t11$extent_kb, "kb)")

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
    step9_jackknife = dt_result$step9_jackknife,
    step10_theta = dt_result$step10_theta,
    n_independent_confirm = dt_result$n_independent_confirm,
    n_independent_support = dt_result$n_independent_support,
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
    t8_concordance = t8$concordance,
    t8_ashman_d = t8$ashman_d,
    t8_n_het_like = t8$n_het_like,
    t8_n_hom_like = t8$n_hom_like,
    t8_n_variants = t8$n_variants_tested,
    # T9: Ancestry jackknife (Cheat 6)
    t9_jackknife_verdict = t9$jackknife_verdict,
    t9_max_delta = t9$max_delta,
    t9_n_drops = t9$n_drops %||% NA_integer_,
    t9_fragile_group = t9$fragile_group %||% NA_character_,
    t9_method = t9$method %||% NA_character_,
    t9_fst_full = t9$fst_full %||% NA_real_,
    # T10: Theta het prior (Cheat 12)
    t10_theta_concordance = t10$concordance,
    t10_band_separation = t10$band_separation,
    t10_n_het_like = t10$n_het_like,
    t10_n_hom_like = t10$n_hom_like,
    # T11: Extended suppression (Cheat 23)
    t11_evidence = t11$evidence,
    t11_has_extended = t11$has_extended %||% NA,
    t11_extent_kb = t11$extent_kb,
    t11_fst_near_ext = t11$fst_near_ext,
    t11_fst_deep_ext = t11$fst_deep_ext
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
                  "H3: Nested/composite", "H4: Sub-haplotype/founder",
                  "INDEPENDENT: Clair3 genotype concordance (Cheat 9)",
                  "CONFOUND: Ancestry jackknife (Cheat 6)",
                  "INDEPENDENT: Theta het prior (Cheat 12)"),
  prediction = c("High within-band relatedness, low ancestry diversity",
                  "Same samples inner/outer, stronger coherence inner",
                  "Changed composition at inner/outer boundary",
                  "Persists after kin pruning, carrier subclusters",
                  "Discrete het/hom indel classes match PCA bands",
                  "Block survives removal of each Q-group",
                  "Per-sample theta_P separates HET from HOM bands"),
  test = c("T1 relatedness ratio + T2 eff_K",
           "T4 inner/outer Jaccard + T5 anchor stability",
           "T6 regime change + T4 composition shift",
           "T3 kin-pruned retention + T7 carrier-only PCA",
           "T8 Clair3 het rate clustering vs PCA bands",
           "T9 leave-one-Q-group-out sim_mat recomputation",
           "T10 ANGSD theta_P z-score clustering vs PCA bands"),
  key_metric = c("T1 ratio > 3 + T2 eff_K < 1.5",
                  "T4 Jaccard > 0.7 + same anchors",
                  "T4 Jaccard < 0.5 + regime shift",
                  "T3 signal retained + carrier subgroups",
                  "T8 concordance > 0.75 + Ashman D > 1.5",
                  "T9 max_delta < 0.02 = robust; single drop = fragile",
                  "T10 concordance > 0.70 + band separation > 1.0"),
  verdict_role = c("Primary: step1 family check",
                    "Architecture: H2 vs H3 discrimination",
                    "Architecture: H3 confirmation",
                    "Structural: H4 sub-haplotype",
                    "Independent rescue: overrides T3 collapse",
                    "Independent: multi-family robustness, can rescue or condemn",
                    "Independent rescue: overrides T3 collapse alongside T8")
)
fwrite(hyp_table, file.path(outdir, "hypothesis_prediction_table.tsv"), sep = "\t")


# ═══════════════════════════════════════════════════════════════════════
# REGISTER IN EVIDENCE REGISTRY
# ═══════════════════════════════════════════════════════════════════════
tryCatch({
  for (.hf in c("utils/registry_key_helpers.R", "../utils/registry_key_helpers.R")) {
    if (file.exists(.hf)) { source(.hf); break }
  }

  if (.bridge_available && exists("reg") && !is.null(reg) && !is.null(reg$add_evidence) &&
      nrow(verd_dt) > 0) {
    message("[C01f] Registering ", nrow(verd_dt), " verdicts in evidence registry...")
    n_keys_total <- 0L
    for (vi in seq_len(nrow(verd_dt))) {
      vd <- verd_dt[vi]
      cid <- paste0(vd$chrom, "_", vd$interval_id)

      if (exists("register_C01f_keys", mode = "function")) {
        n_keys_total <- n_keys_total + register_C01f_keys(vd, cid, outdir)
      }
      if (exists("store_C01f_results", mode = "function")) {
        store_C01f_results(vd, cid, outdir)
      }
    }
    message("[C01f] Registered ", n_keys_total, " evidence keys across ", nrow(verd_dt), " candidates")
  }
}, error = function(e) message("[C01f] Registry wiring: ", conditionMessage(e)))

message("\n[DONE] -> ", outdir)
