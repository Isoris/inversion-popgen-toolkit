#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01f_b_negative_controls.R  (v8.4)
#
# NEGATIVE CONTROLS for hypothesis tests.
#
# For each test in C01f (T1-T8), compute the SAME metric on regions where
# there is NO inversion signal. This establishes the null distribution:
#   - What does T1 (relatedness ratio) look like in random background?
#   - What does T3 (kin-pruned retention) look like when there's no block?
#   - What does T8 (informativeness fold) look like genome-wide?
#
# Without this, all thresholds in C01f are arbitrary. With this, we can say:
#   "T3 retention = 0.75 for candidate X; background mean = 0.12 +/- 0.08"
#   That's a real result.
#
# Control regions are selected as:
#   A. Background windows: random non-candidate windows across the genome
#   B. Ancestry-matched controls: regions where ancestry structure is strong
#      but NO triangle/block structure exists (rules out ancestry instability)
#   C. Family-enriched controls: regions where relatedness is high but
#      NO PCA signal exists (rules out family artifact)
#
# Also tests ANCESTRY INSTABILITY:
#   For each candidate, measure how much band membership correlates with
#   genome-wide ancestry (Q-matrix). If high correlation = ancestry drives
#   the bands, not a local inversion.
#
# Inputs: same as C01f + candidate scores from C01d
#
# Outputs:
#   null_distributions.tsv.gz     -- per-metric null distributions
#   control_comparisons.tsv       -- candidate vs null (z-scores / percentiles)
#   ancestry_instability.tsv      -- per-candidate ancestry correlation
#   plots/null_vs_candidate.png   -- visual comparison
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; precomp_dir <- NULL
samples_file <- NULL; relate_file <- NULL; pruned_file <- NULL
ancestry_file <- NULL; outdir <- "negative_controls"
n_controls <- 50L  # number of random control regions per chromosome

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))      { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args))  { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp" && i < length(args))    { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args)) { relate_file <- args[i+1]; i <- i+2 }
  else if (a == "--pruned_samples" && i < length(args)) { pruned_file <- args[i+1]; i <- i+2 }
  else if (a == "--ancestry" && i < length(args))   { ancestry_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--n_controls" && i < length(args)) { n_controls <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file) || is.null(precomp_dir)) stop("--scores and --precomp required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 300

# =============================================================================
# LOAD
# =============================================================================

message("[C01f_b] Loading data...")
cand_dt <- fread(scores_file)

# Intervals (to know where candidates ARE, so controls avoid them)
iv_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_intervals.tsv.gz")
  if (file.exists(f)) iv_dt <- fread(f)
}

# Sample names + mapping
real_names <- NULL
if (!is.null(samples_file) && file.exists(samples_file))
  real_names <- as.character(fread(samples_file, header = FALSE)[[1]])

# Relatedness
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

# Pruned samples
pruned_samples <- NULL; pruned_ind_names <- NULL
if (!is.null(pruned_file) && file.exists(pruned_file)) {
  pruned_samples <- as.character(fread(pruned_file, header = FALSE)[[1]])
  if (!is.null(real_names) && length(real_names) > 0) {
    idx <- match(pruned_samples, real_names)
    pruned_ind_names <- paste0("Ind", idx[!is.na(idx)] - 1L)
  }
}

# Q-matrix for ancestry instability test
q_mat <- NULL
if (!is.null(ancestry_file) && file.exists(ancestry_file)) {
  q_raw <- as.matrix(fread(ancestry_file, header = FALSE))
  if (!is.null(real_names) && nrow(q_raw) == length(real_names)) {
    rownames(q_raw) <- real_names
    q_mat <- q_raw
    message("[C01f_b] Q-matrix: ", nrow(q_mat), " x K=", ncol(q_mat))
  }
}

message("[C01f_b] Candidates: ", nrow(cand_dt), " | Controls per chr: ", n_controls)

# =============================================================================
# GENERATE CONTROL REGIONS
# =============================================================================

generate_controls <- function(precomp_dir, iv_dt, n_per_chr, cand_dt) {
  rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
  controls <- list()

  for (f in rds_files) {
    pc <- readRDS(f)
    chr <- pc$chrom; dt <- pc$dt; n_w <- nrow(dt)
    if (n_w < 100) next

    # Candidate windows to AVOID
    cand_wins <- integer(0)
    for (ci in seq_len(nrow(cand_dt[chrom == chr]))) {
      c <- cand_dt[chrom == chr][ci]
      idx <- which(dt$start_bp >= c$start_mb * 1e6 & dt$end_bp <= c$end_mb * 1e6)
      cand_wins <- c(cand_wins, idx)
    }
    # Also avoid triangle intervals
    for (ti in seq_len(nrow(iv_dt[chrom == chr]))) {
      iv <- iv_dt[chrom == chr][ti]
      idx <- which(dt$start_bp >= iv$start_bp & dt$end_bp <= iv$end_bp)
      cand_wins <- c(cand_wins, idx)
    }
    cand_wins <- unique(cand_wins)
    bg_wins <- setdiff(seq_len(n_w), cand_wins)
    if (length(bg_wins) < 50) next

    # Sample random control regions (same size as median candidate)
    median_span <- if (nrow(cand_dt) > 0) median(cand_dt$span_mb * 1e6 /
                    ((dt$end_bp[2] - dt$start_bp[1]) %||% 1e4)) else 40
    ctrl_size <- max(20L, as.integer(median_span))

    for (ci in seq_len(n_per_chr)) {
      start_idx <- sample(bg_wins[bg_wins <= n_w - ctrl_size], 1)
      end_idx <- start_idx + ctrl_size - 1
      if (end_idx > n_w) next
      # Check no overlap with candidates
      if (any(start_idx:end_idx %in% cand_wins)) next

      controls[[length(controls) + 1]] <- list(
        chrom = chr, start_idx = start_idx, end_idx = end_idx,
        start_bp = dt$start_bp[start_idx], end_bp = dt$end_bp[end_idx],
        start_mb = dt$start_bp[start_idx] / 1e6, end_mb = dt$end_bp[end_idx] / 1e6,
        pc = pc
      )
    }
  }
  controls
}

message("[C01f_b] Generating control regions...")
controls <- generate_controls(precomp_dir, iv_dt, n_controls, cand_dt)
message("[C01f_b] Generated ", length(controls), " control regions")

# =============================================================================
# COMPUTE NULL METRICS ON CONTROL REGIONS
# =============================================================================

compute_control_metrics <- function(ctrl, relate_dt, pruned_ind_names, q_mat, real_names) {
  pc <- ctrl$pc; dt <- pc$dt; sim_mat <- pc$sim_mat
  chr <- ctrl$chrom; wi_s <- ctrl$start_idx; wi_e <- ctrl$end_idx
  n_w <- nrow(sim_mat)

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  sample_names <- sub("^PC_1_", "", pc1_cols)
  available <- intersect(pc1_cols, names(dt))

  # T1-like: k=3 cluster in control region, compute within/between relatedness
  t1_ratio <- NA_real_
  if (length(available) >= 20 && !is.null(relate_dt)) {
    mat <- as.matrix(dt[wi_s:wi_e, ..available])
    avg <- colMeans(mat, na.rm = TRUE); valid <- is.finite(avg)
    if (sum(valid) >= 20) {
      vals <- avg[valid]; snames <- sub("^PC_1_", "", names(vals))
      km <- tryCatch(kmeans(vals, centers = 3, nstart = 5), error = function(e) NULL)
      if (!is.null(km)) {
        co <- order(km$centers[, 1])
        bands <- character(length(vals))
        for (bi in seq_along(co)) bands[km$cluster == co[bi]] <- paste0("band", bi)
        # Map to real names for relatedness
        if (!is.null(real_names) && grepl("^Ind", snames[1])) {
          ind_to_real <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
          snames_r <- ind_to_real[snames]; snames_r[is.na(snames_r)] <- snames[is.na(snames_r)]
        } else snames_r <- snames

        within_v <- numeric(0); between_v <- numeric(0)
        for (b in paste0("band", 1:3)) {
          bs <- snames_r[bands == b]
          wp <- relate_dt[sample1 %in% bs & sample2 %in% bs]
          if (nrow(wp) > 0) within_v <- c(within_v, wp$theta)
          other <- snames_r[bands != b]
          bp <- relate_dt[(sample1 %in% bs & sample2 %in% other) |
                           (sample2 %in% bs & sample1 %in% other)]
          if (nrow(bp) > 0) between_v <- c(between_v, bp$theta)
        }
        w <- mean(within_v); b <- mean(between_v)
        if (is.finite(w) && is.finite(b) && b > 0) t1_ratio <- w / b
      }
    }
  }

  # T3-like: block contrast (should be near 0 for controls)
  t3_retention <- NA_real_
  if (!is.null(pruned_ind_names) && length(pruned_ind_names) >= 15) {
    all_idx <- seq_len(n_w); out_idx <- setdiff(all_idx, wi_s:wi_e)
    if (length(out_idx) >= 10) {
      full_inside <- mean(sim_mat[wi_s:wi_e, wi_s:wi_e], na.rm = TRUE)
      full_outside <- mean(sim_mat[out_idx[1:min(100, length(out_idx))],
                                    out_idx[1:min(100, length(out_idx))]], na.rm = TRUE)
      full_contrast <- full_inside - full_outside
      # Pruned: recompute from PC loadings of pruned samples only
      pruned_cols <- paste0("PC_1_", pruned_ind_names)
      avail_p <- intersect(pruned_cols, names(dt))
      if (length(avail_p) >= 15) {
        sub_idx <- sort(c(sample(out_idx, min(100, length(out_idx))),
                           wi_s:min(wi_e, n_w)))
        p_mat <- as.matrix(dt[sub_idx, ..avail_p])
        ns <- nrow(p_mat); p_sim <- matrix(0, ns, ns)
        for (ii in seq_len(ns)) for (jj in ii:ns) {
          v <- cor(p_mat[ii, ], p_mat[jj, ], use = "pairwise.complete.obs")
          if (is.finite(v)) { p_sim[ii, jj] <- v; p_sim[jj, ii] <- v }
        }
        win_in <- which(sub_idx %in% (wi_s:wi_e))
        out_in <- which(!sub_idx %in% (wi_s:wi_e))
        if (length(win_in) >= 5 && length(out_in) >= 5) {
          p_inside <- mean(p_sim[win_in, win_in], na.rm = TRUE)
          p_outside <- mean(p_sim[out_in, out_in], na.rm = TRUE)
          p_contrast <- p_inside - p_outside
          if (full_contrast > 0.005) t3_retention <- p_contrast / full_contrast
        }
      }
    }
  }

  # T8-like: eigenvalue ratio fold (should be ~1 for controls)
  t8_fold <- NA_real_
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else NULL
  if (!is.null(l1_col)) {
    inner_r <- dt[[l1_col]][wi_s:wi_e] / pmax(dt[[l2_col]][wi_s:wi_e], 0.01)
    outer_r <- dt[[l1_col]][setdiff(seq_len(n_w), wi_s:wi_e)] /
               pmax(dt[[l2_col]][setdiff(seq_len(n_w), wi_s:wi_e)], 0.01)
    im <- median(inner_r, na.rm = TRUE); om <- median(outer_r, na.rm = TRUE)
    if (om > 0) t8_fold <- im / om
  }

  # Ancestry instability: correlation between PC1 band and dominant Q-component
  ancestry_cor <- NA_real_
  if (!is.null(q_mat) && length(available) >= 20) {
    mat <- as.matrix(dt[wi_s:wi_e, ..available])
    avg <- colMeans(mat, na.rm = TRUE); valid <- is.finite(avg)
    if (sum(valid) >= 20) {
      vals <- avg[valid]; snames <- sub("^PC_1_", "", names(vals))
      if (!is.null(real_names) && grepl("^Ind", snames[1])) {
        rn <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
        snames_r <- rn[snames]
      } else snames_r <- snames
      shared <- intersect(snames_r[!is.na(snames_r)], rownames(q_mat))
      if (length(shared) > 20) {
        pc1_shared <- vals[match(shared, snames_r)]
        # Dominant Q-component for each sample
        dom_q <- apply(q_mat[shared, , drop = FALSE], 1, max)
        ancestry_cor <- abs(cor(pc1_shared, dom_q, use = "complete.obs"))
      }
    }
  }

  data.table(
    chrom = chr, start_mb = ctrl$start_mb, end_mb = ctrl$end_mb,
    type = "control",
    t1_ratio = round(t1_ratio, 4),
    t3_retention = round(t3_retention, 4),
    t8_fold = round(t8_fold, 4),
    ancestry_cor = round(ancestry_cor, 4)
  )
}

message("[C01f_b] Computing null metrics on ", length(controls), " control regions...")
null_rows <- list()
for (ci in seq_along(controls)) {
  if (ci %% 50 == 0) message("  Control ", ci, "/", length(controls))
  null_rows[[ci]] <- compute_control_metrics(controls[[ci]], relate_dt, pruned_ind_names,
                                              q_mat, real_names)
}
null_dt <- rbindlist(null_rows, fill = TRUE)
message("[C01f_b] Null distributions computed: ", nrow(null_dt), " control regions")

# =============================================================================
# ANCESTRY INSTABILITY TEST ON CANDIDATES
# =============================================================================

message("[C01f_b] Computing ancestry instability for candidates...")
anc_rows <- list()
precomp_cache <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]; chr <- cand$chrom
  if (is.null(precomp_cache[[chr]])) {
    f <- list.files(precomp_dir, pattern = paste0(chr, "\\.precomp\\.rds$"), full.names = TRUE)
    if (length(f) > 0) precomp_cache[[chr]] <- readRDS(f[1])
  }
  pc <- precomp_cache[[chr]]
  if (is.null(pc) || is.null(q_mat)) next
  dt <- pc$dt
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  available <- intersect(pc1_cols, names(dt))
  sample_names <- sub("^PC_1_", "", pc1_cols)

  win_idx <- which(dt$start_bp >= cand$start_mb * 1e6 & dt$end_bp <= cand$end_mb * 1e6)
  if (length(win_idx) < 5) next

  mat <- as.matrix(dt[win_idx, ..available])
  avg <- colMeans(mat, na.rm = TRUE); valid <- is.finite(avg)
  if (sum(valid) < 20) next
  vals <- avg[valid]; snames <- sub("^PC_1_", "", names(vals))

  if (!is.null(real_names) && grepl("^Ind", snames[1])) {
    rn <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
    snames_r <- rn[snames]
  } else snames_r <- snames

  shared <- intersect(snames_r[!is.na(snames_r)], rownames(q_mat))
  if (length(shared) < 20) next

  pc1_shared <- vals[match(shared, snames_r)]
  # Test correlation with EACH Q-component, not just dominant
  q_cors <- numeric(ncol(q_mat))
  for (ki in seq_len(ncol(q_mat))) {
    q_cors[ki] <- abs(cor(pc1_shared, q_mat[shared, ki], use = "complete.obs"))
  }
  max_q_cor <- max(q_cors, na.rm = TRUE)
  dominant_k <- which.max(q_cors)

  # Null comparison
  null_ancestry <- null_dt[is.finite(ancestry_cor)]$ancestry_cor
  percentile <- if (length(null_ancestry) > 5) {
    mean(null_ancestry <= max_q_cor)
  } else NA_real_

  # Verdict
  anc_verdict <- if (max_q_cor > 0.6 && !is.na(percentile) && percentile > 0.95) {
    "ancestry_driven"  # PC1 bands strongly track ancestry
  } else if (max_q_cor > 0.4) {
    "ancestry_correlated"
  } else {
    "ancestry_independent"
  }

  anc_rows[[length(anc_rows) + 1]] <- data.table(
    chrom = chr, interval_id = cand$interval_id,
    start_mb = cand$start_mb, end_mb = cand$end_mb,
    max_ancestry_cor = round(max_q_cor, 4),
    dominant_K = dominant_k,
    per_K_cors = paste(round(q_cors, 3), collapse = ","),
    null_percentile = round(percentile, 3),
    ancestry_verdict = anc_verdict
  )
}
anc_dt <- if (length(anc_rows) > 0) rbindlist(anc_rows) else data.table()

# =============================================================================
# COMPARE CANDIDATES VS NULL
# =============================================================================

# Load candidate hypothesis test results if available
hyp_file <- NULL
hi <- match("--hyp_results", args)
if (!is.na(hi) && hi < length(args)) hyp_file <- args[hi + 1]
hyp_dt <- data.table()
if (!is.null(hyp_file) && file.exists(hyp_file)) hyp_dt <- fread(hyp_file)

compare_rows <- list()
for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]; chr <- cand$chrom; iid <- cand$interval_id

  # Get candidate metrics from hypothesis tests
  cand_t1 <- NA; cand_t3 <- NA; cand_t8 <- NA
  if (nrow(hyp_dt) > 0) {
    hv <- hyp_dt[chrom == chr & interval_id == iid]
    if (nrow(hv) > 0) {
      if ("t1_ratio" %in% names(hv)) cand_t1 <- hv$t1_ratio[1]
      if ("t3_retention" %in% names(hv)) cand_t3 <- hv$t3_retention[1]
      if ("t8_fold" %in% names(hv)) cand_t8 <- hv$t8_fold[1]
    }
  }

  # Z-scores against null
  null_t1 <- null_dt[is.finite(t1_ratio)]$t1_ratio
  null_t3 <- null_dt[is.finite(t3_retention)]$t3_retention
  null_t8 <- null_dt[is.finite(t8_fold)]$t8_fold

  z_t1 <- if (is.finite(cand_t1) && length(null_t1) > 5)
    (cand_t1 - mean(null_t1)) / max(sd(null_t1), 0.01) else NA
  z_t3 <- if (is.finite(cand_t3) && length(null_t3) > 5)
    (cand_t3 - mean(null_t3)) / max(sd(null_t3), 0.01) else NA
  z_t8 <- if (is.finite(cand_t8) && length(null_t8) > 5)
    (cand_t8 - mean(null_t8)) / max(sd(null_t8), 0.01) else NA

  p_t1 <- if (is.finite(cand_t1) && length(null_t1) > 5) mean(null_t1 >= cand_t1) else NA
  p_t3 <- if (is.finite(cand_t3) && length(null_t3) > 5) mean(null_t3 >= cand_t3) else NA
  p_t8 <- if (is.finite(cand_t8) && length(null_t8) > 5) mean(null_t8 >= cand_t8) else NA

  compare_rows[[ci]] <- data.table(
    chrom = chr, interval_id = iid, tier = cand$tier,
    cand_t1 = cand_t1, null_t1_mean = round(mean(null_t1), 4),
    null_t1_sd = round(sd(null_t1), 4), z_t1 = round(z_t1, 2), p_t1 = round(p_t1, 4),
    cand_t3 = cand_t3, null_t3_mean = round(mean(null_t3), 4),
    null_t3_sd = round(sd(null_t3), 4), z_t3 = round(z_t3, 2), p_t3 = round(p_t3, 4),
    cand_t8 = cand_t8, null_t8_mean = round(mean(null_t8), 4),
    null_t8_sd = round(sd(null_t8), 4), z_t8 = round(z_t8, 2), p_t8 = round(p_t8, 4)
  )
}
compare_dt <- if (length(compare_rows) > 0) rbindlist(compare_rows, fill = TRUE) else data.table()

# =============================================================================
# WRITE
# =============================================================================

message("[C01f_b] Writing outputs...")
fwrite(null_dt, file.path(outdir, "null_distributions.tsv.gz"), sep = "\t")
fwrite(compare_dt, file.path(outdir, "control_comparisons.tsv"), sep = "\t")
fwrite(anc_dt, file.path(outdir, "ancestry_instability.tsv"), sep = "\t")

# Summary
message("\n[C01f_b] === NULL DISTRIBUTION SUMMARY ===")
for (metric in c("t1_ratio", "t3_retention", "t8_fold", "ancestry_cor")) {
  vals <- null_dt[[metric]]
  vals <- vals[is.finite(vals)]
  if (length(vals) > 0)
    message("  ", metric, ": mean=", round(mean(vals), 4),
            " sd=", round(sd(vals), 4),
            " q05=", round(quantile(vals, 0.05), 4),
            " q95=", round(quantile(vals, 0.95), 4))
}

if (nrow(anc_dt) > 0) {
  message("\n[C01f_b] === ANCESTRY INSTABILITY ===")
  for (v in sort(unique(anc_dt$ancestry_verdict)))
    message("  ", v, ": ", sum(anc_dt$ancestry_verdict == v))
}

message("\n[DONE] -> ", outdir)
