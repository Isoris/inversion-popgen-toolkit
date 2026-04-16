#!/usr/bin/env Rscript

# =============================================================================
# STEP41_candidate_membership_trajectories.R  (v1.0)
#
# PER-SAMPLE MEMBERSHIP TRACKING across windows within a candidate.
#
# The "graph of belonging" — for each sample, which structural sub-cluster
# does it belong to in each window? Does it stay put or switch roles?
#
# DISTINGUISHES TWO HYPOTHESES:
#
#   H1 — STABLE SUB-VARIANTS:
#     The B arrangement has internal variants (B1, B2, B3).
#     Each sample carries a fixed combination.
#     Sub-cluster membership is constant across all windows.
#     Evidence: high trajectory stability, low switching rate.
#
#   H2 — POSITION-DEPENDENT SYSTEM CHANGES:
#     Multiple arrangement systems overlap in the candidate interval.
#     Some samples start as Homo_2(system1) but transition to
#     Homo_2(system2) at a different position.
#     Evidence: switching events clustered at specific positions,
#     trajectory instability in certain interval stretches.
#
# METHOD:
#   For each window within the candidate:
#     1. Run local PCA on the window's dosage data
#     2. Identify sub-clusters within each broad group (Homo_1/Het/Homo_2)
#        using k-means on PC2 (the secondary axis WITHIN each group)
#     3. Assign each sample a sub-cluster label
#   Then across windows:
#     4. Track each sample's sub-cluster trajectory
#     5. Compute stability metrics and detect switching events
#     6. Classify as H1 (stable) or H2 (position-dependent)
#
# INPUTS:
#   <config.R>     — standard followup config
#   <cid>          — candidate ID
#   <dosage_dir>   — per-chr dosage files
#
# OUTPUTS:
#   candidate_membership_matrix.tsv.gz — sample × window sub-cluster assignments
#   candidate_switching_events.tsv.gz  — detected sub-cluster switches
#   candidate_trajectory_summary.tsv   — per-sample stability metrics
#   candidate_membership_verdict.tsv   — H1 vs H2 classification
#   candidate_membership_heatmap.pdf   — sample × window heatmap of membership
#
# Usage:
#   Rscript STEP41_candidate_membership_trajectories.R <config.R> <cid>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: STEP41 <config.R> <cid>")

config_file <- args[1]
cid <- as.integer(args[2])
source(config_file)

# Parameters
N_SUB_MAX    <- 4L    # max sub-clusters to look for within each broad group
MIN_SUB_SIZE <- 3L    # min samples in a sub-cluster to count
SWITCH_GAP   <- 2L    # gap tolerance for switch detection (windows)

cand <- fread(CANDIDATE_TABLE)
row <- cand[candidate_id == cid]
if (nrow(row) == 0) stop("Candidate ", cid, " not found")
chr <- row$chrom; c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))
cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)

message("[STEP41] Membership trajectories for candidate ", cid)
message("[STEP41] ", chr, ":", round(c_start/1e6, 2), "–", round(c_end/1e6, 2), " Mb")

# =============================================================================
# LOAD DOSAGE + GROUP ASSIGNMENTS
# =============================================================================

dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
if (!file.exists(dos_file) || !file.exists(sites_file)) stop("Missing dosage for ", chr)

dos <- fread(dos_file); sites <- fread(sites_file)
sample_names <- fread(SAMPLES_IND, header = FALSE)[[1]]
sample_names <- sample_names[nchar(sample_names) > 0]
dos_cols <- setdiff(names(dos), "marker")
if (all(grepl("^Ind", dos_cols)) && length(sample_names) == length(dos_cols)) {
  setnames(dos, old = dos_cols, new = sample_names); dos_cols <- sample_names
}

# Load group assignments (from STEP21 or similar)
group_dt <- NULL
gf <- file.path(cand_dir, "candidate_group_assignments.tsv")
if (!file.exists(gf)) gf <- file.path(cand_dir, "candidate_pca_rotated.tsv")
if (file.exists(gf)) {
  group_dt <- fread(gf)
  # Standardize column names
  gcol <- intersect(names(group_dt), c("group", "coarse_group", "coarse_group_refined"))
  scol <- intersect(names(group_dt), c("sample", "sample_id"))
  if (length(gcol) > 0 && length(scol) > 0) {
    group_dt <- group_dt[, .SD, .SDcols = c(scol[1], gcol[1])]
    setnames(group_dt, c("sample", "broad_group"))
  } else group_dt <- NULL
}

if (is.null(group_dt)) {
  message("[STEP41] No group assignments found — will use PCA-based grouping per window")
}

# =============================================================================
# DEFINE WINDOWS WITHIN CANDIDATE
# =============================================================================

# Load MDS data for window positions
mds_rds <- file.path(INV_ROOT, "06_mds_candidates", "inversion_localpca.mds.rds")
if (!file.exists(mds_rds)) stop("Missing: ", mds_rds)
mds_obj <- readRDS(mds_rds)

chr_dt <- as.data.table(mds_obj$per_chr[[chr]]$out_dt)
chr_dt <- chr_dt[order(start_bp)]
cand_windows <- chr_dt[start_bp >= c_start & end_bp <= c_end]
n_win <- nrow(cand_windows)
message("[STEP41] Windows in candidate: ", n_win)

if (n_win < 5) {
  message("[STEP41] Too few windows — skipping")
  fwrite(data.table(candidate_id = cid, verdict = "TOO_FEW_WINDOWS"),
         file.path(cand_dir, "candidate_membership_verdict.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

# =============================================================================
# PER-WINDOW SUB-CLUSTERING
# =============================================================================

# For each window: run local PCA, split into broad groups, then find
# sub-clusters within each group using PC2/PC3

membership_rows <- list()  # sample × window

for (wi in seq_len(n_win)) {
  w_start <- cand_windows$start_bp[wi]
  w_end   <- cand_windows$end_bp[wi]
  wid     <- cand_windows$global_window_id[wi]

  # Extract markers
  w_keep <- which(sites$pos >= w_start & sites$pos <= w_end)
  if (length(w_keep) < 15) {
    # Not enough markers — assign NA for this window
    for (si in seq_along(sample_names)) {
      membership_rows[[length(membership_rows) + 1]] <- data.table(
        sample = sample_names[si], window_idx = wi, window_id = wid,
        start_bp = w_start, end_bp = w_end,
        broad_group = NA_character_, sub_cluster = NA_character_,
        pc1 = NA_real_, pc2 = NA_real_
      )
    }
    next
  }

  Xw <- as.matrix(dos[w_keep, ..dos_cols])
  storage.mode(Xw) <- "double"

  # Quick local PCA (covariance on centered dosage)
  Xc <- Xw - rowMeans(Xw, na.rm = TRUE)
  covmat <- suppressWarnings(cov(Xc, use = "pairwise.complete.obs"))
  covmat[is.na(covmat)] <- 0

  ee <- tryCatch(eigen(covmat, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ee) || length(ee$values) < 2) next

  pc1 <- ee$vectors[, 1]; pc2 <- ee$vectors[, 2]

  # Broad 3-way grouping on PC1
  km3 <- tryCatch(kmeans(pc1, centers = 3, nstart = 5), error = function(e) NULL)
  if (is.null(km3)) next

  co <- order(km3$centers[, 1])
  broad_map <- c("Homo_1", "Het", "Homo_2")
  broad_labels <- broad_map[match(km3$cluster, co)]

  # Sub-clustering within each broad group on PC2
  sub_labels <- rep(NA_character_, length(pc1))
  for (bg in broad_map) {
    bg_idx <- which(broad_labels == bg)
    if (length(bg_idx) < MIN_SUB_SIZE) {
      sub_labels[bg_idx] <- paste0(bg, "_sub1")
      next
    }

    bg_pc2 <- pc2[bg_idx]

    # Determine optimal k for sub-clusters (up to N_SUB_MAX)
    # Use gap in sorted PC2 values as a simple heuristic
    best_k <- 1L
    if (length(bg_idx) >= 2 * MIN_SUB_SIZE) {
      sorted_pc2 <- sort(bg_pc2)
      gaps <- diff(sorted_pc2)
      # A sub-cluster split if the largest gap is > 2× median gap
      if (length(gaps) > 0) {
        med_gap <- median(gaps)
        big_gaps <- sum(gaps > 2 * med_gap & med_gap > 0)
        best_k <- min(N_SUB_MAX, big_gaps + 1L)
      }
    }

    if (best_k > 1 && length(bg_idx) >= best_k * MIN_SUB_SIZE) {
      km_sub <- tryCatch(kmeans(bg_pc2, centers = best_k, nstart = 5),
                          error = function(e) NULL)
      if (!is.null(km_sub)) {
        sub_order <- order(km_sub$centers[, 1])
        sub_labels[bg_idx] <- paste0(bg, "_sub", match(km_sub$cluster, sub_order))
      } else {
        sub_labels[bg_idx] <- paste0(bg, "_sub1")
      }
    } else {
      sub_labels[bg_idx] <- paste0(bg, "_sub1")
    }
  }

  # Record
  for (si in seq_along(sample_names)) {
    membership_rows[[length(membership_rows) + 1]] <- data.table(
      sample = sample_names[si], window_idx = wi, window_id = wid,
      start_bp = w_start, end_bp = w_end,
      broad_group = broad_labels[si],
      sub_cluster = sub_labels[si],
      pc1 = round(pc1[si], 4), pc2 = round(pc2[si], 4)
    )
  }
}

if (length(membership_rows) == 0) {
  message("[STEP41] No valid windows — skipping")
  quit(save = "no", status = 0)
}

mem_dt <- rbindlist(membership_rows, fill = TRUE)
message("[STEP41] Membership matrix: ", nrow(mem_dt), " entries (",
        uniqueN(mem_dt$sample), " samples × ", n_win, " windows)")

# =============================================================================
# DETECT SWITCHING EVENTS
# =============================================================================

# For each sample, track sub_cluster trajectory across windows
# A "switch" = sub_cluster changes between adjacent (or near-adjacent) windows

switches <- list()
stability <- list()

for (samp in unique(mem_dt$sample)) {
  traj <- mem_dt[sample == samp][order(window_idx)]
  traj <- traj[!is.na(sub_cluster)]
  if (nrow(traj) < 3) next

  # Count switches
  n_switch <- 0L
  switch_positions <- c()
  for (ti in 2:nrow(traj)) {
    if (traj$sub_cluster[ti] != traj$sub_cluster[ti - 1]) {
      n_switch <- n_switch + 1L
      switch_positions <- c(switch_positions, traj$start_bp[ti])
      switches[[length(switches) + 1]] <- data.table(
        sample = samp,
        from_window = traj$window_idx[ti - 1],
        to_window = traj$window_idx[ti],
        switch_bp = traj$start_bp[ti],
        from_cluster = traj$sub_cluster[ti - 1],
        to_cluster = traj$sub_cluster[ti],
        from_broad = traj$broad_group[ti - 1],
        to_broad = traj$broad_group[ti]
      )
    }
  }

  # Stability = fraction of windows with no switch
  n_valid <- nrow(traj)
  frac_stable <- (n_valid - n_switch) / n_valid

  # Dominant sub-cluster
  dom_tab <- table(traj$sub_cluster)
  dominant <- names(dom_tab)[which.max(dom_tab)]
  frac_dominant <- as.numeric(max(dom_tab)) / n_valid

  stability[[length(stability) + 1]] <- data.table(
    sample = samp,
    n_valid_windows = n_valid,
    n_switches = n_switch,
    switch_rate = round(n_switch / max(1, n_valid - 1), 4),
    frac_stable = round(frac_stable, 4),
    dominant_sub_cluster = dominant,
    frac_dominant = round(frac_dominant, 4),
    broad_group = if (!is.null(group_dt) && samp %in% group_dt$sample)
      group_dt[sample == samp]$broad_group[1] else traj$broad_group[1]
  )
}

switch_dt <- if (length(switches) > 0) rbindlist(switches) else data.table()
stab_dt   <- if (length(stability) > 0) rbindlist(stability) else data.table()

message("[STEP41] Switching events: ", nrow(switch_dt))
message("[STEP41] Samples with switches: ", sum(stab_dt$n_switches > 0))

# =============================================================================
# VERDICT: H1 (stable sub-variants) vs H2 (position-dependent changes)
# =============================================================================

if (nrow(stab_dt) > 0) {
  # Overall stability
  mean_stability <- mean(stab_dt$frac_stable)
  frac_switchers <- mean(stab_dt$n_switches > 0)
  mean_switch_rate <- mean(stab_dt$switch_rate)

  # Are switches clustered at specific positions?
  switch_clustering <- "dispersed"
  if (nrow(switch_dt) >= 5) {
    switch_bps <- switch_dt$switch_bp
    # Check if switches cluster at a few hotspot positions
    bp_density <- density(switch_bps, bw = (c_end - c_start) / 20)
    peak_height <- max(bp_density$y)
    baseline_height <- mean(bp_density$y)
    if (peak_height > 3 * baseline_height) {
      switch_clustering <- "hotspot"  # switches concentrated at specific positions
    }
  }

  # Verdict
  if (mean_stability >= 0.85 && frac_switchers < 0.15) {
    h_verdict <- "H1_STABLE_SUBVARIANTS"
    h_detail <- paste0("Sub-cluster membership is stable across windows (mean stability=",
                        round(mean_stability, 2), ", ", round(frac_switchers*100, 1),
                        "% switchers). Sub-clusters likely represent fixed sub-variants ",
                        "within the arrangement, not position-dependent system changes.")
    h_recommendation <- "Analyze with STEP33 within-stripe analysis to characterize sub-variants."
  } else if (switch_clustering == "hotspot") {
    h_verdict <- "H2_POSITIONAL_HOTSPOT"
    h_detail <- paste0("Switches are concentrated at specific genomic positions (mean stability=",
                        round(mean_stability, 2), ", ", round(frac_switchers*100, 1),
                        "% switchers). This suggests the candidate contains multiple arrangement ",
                        "systems that transition at specific breakpoints.")
    h_recommendation <- paste0("Split candidate at switch hotspot positions. ",
                               "Each sub-interval may represent a different system.")
  } else if (frac_switchers >= 0.30) {
    h_verdict <- "H2_WIDESPREAD_INSTABILITY"
    h_detail <- paste0("Many samples switch sub-clusters across windows (mean stability=",
                        round(mean_stability, 2), ", ", round(frac_switchers*100, 1),
                        "% switchers). This suggests a genuinely composite candidate ",
                        "with multiple overlapping systems.")
    h_recommendation <- "Run STEP40 for window-defined sub-candidates, then re-analyze each."
  } else {
    h_verdict <- "MIXED_SIGNAL"
    h_detail <- paste0("Moderate switching (stability=", round(mean_stability, 2),
                        ", ", round(frac_switchers*100, 1), "% switchers). ",
                        "Neither clearly H1 nor H2.")
    h_recommendation <- "Inspect membership heatmap visually to determine pattern."
  }
} else {
  h_verdict <- "INSUFFICIENT_DATA"
  h_detail <- "Not enough valid trajectories to assess"
  h_recommendation <- "Check input data quality"
  mean_stability <- NA_real_; frac_switchers <- NA_real_
  mean_switch_rate <- NA_real_; switch_clustering <- NA_character_
}

verdict_dt <- data.table(
  candidate_id = cid, chrom = chr,
  n_windows = n_win, n_samples = length(sample_names),
  mean_stability = round(mean_stability, 4),
  frac_switchers = round(frac_switchers, 4),
  mean_switch_rate = round(mean_switch_rate, 4),
  switch_clustering = switch_clustering,
  verdict = h_verdict,
  detail = h_detail,
  recommendation = h_recommendation
)

message("[STEP41] Verdict: ", h_verdict)

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

fwrite(mem_dt, file.path(cand_dir, "candidate_membership_matrix.tsv.gz"), sep = "\t")
fwrite(switch_dt, file.path(cand_dir, "candidate_switching_events.tsv.gz"), sep = "\t")
fwrite(stab_dt, file.path(cand_dir, "candidate_trajectory_summary.tsv"), sep = "\t")
fwrite(verdict_dt, file.path(cand_dir, "candidate_membership_verdict.tsv"), sep = "\t")

# =============================================================================
# MEMBERSHIP HEATMAP
# =============================================================================

message("[STEP41] Building membership heatmap...")

# Convert sub_cluster to numeric for heatmap coloring
all_subclusters <- sort(unique(mem_dt$sub_cluster[!is.na(mem_dt$sub_cluster)]))
mem_dt[, sub_cluster_num := match(sub_cluster, all_subclusters)]

# Create wide matrix: samples × windows
mem_wide <- dcast(mem_dt, sample ~ window_idx, value.var = "sub_cluster_num")

# Order samples by broad group and then by dominant sub-cluster
if (nrow(stab_dt) > 0) {
  sample_order <- stab_dt[order(broad_group, dominant_sub_cluster, -frac_dominant)]$sample
  sample_order <- c(sample_order, setdiff(mem_wide$sample, sample_order))
} else {
  sample_order <- mem_wide$sample
}
mem_wide <- mem_wide[match(sample_order, sample)]

mat <- as.matrix(mem_wide[, -1])
rownames(mat) <- mem_wide$sample

# Plot
n_cols <- length(all_subclusters)
pal <- if (n_cols <= 8) {
  c("#1d4ed8", "#3b82f6", "#93c5fd",  # Homo_1 shades
    "#f59e0b", "#fbbf24", "#fde68a",  # Het shades
    "#dc2626", "#f87171", "#fca5a5",  # Homo_2 shades
    "#6b7280")[seq_len(n_cols)]
} else rainbow(n_cols, s = 0.7)

fig_path <- file.path(cand_dir, "candidate_membership_heatmap.pdf")
n_samp <- nrow(mat)
fig_h <- max(6, n_samp * 0.04 + 2)

pdf(fig_path, width = max(8, n_win * 0.15 + 3), height = fig_h)
par(mar = c(4, 8, 3, 6))
image(t(mat[nrow(mat):1, ]), col = pal, axes = FALSE,
      main = paste0("Candidate ", cid, " — Sub-cluster membership (",
                     h_verdict, ")"),
      xlab = "Window position →", ylab = "")
# Sample labels (if not too many)
if (n_samp <= 80) {
  axis(2, at = seq(0, 1, length.out = n_samp),
       labels = rev(rownames(mat)), las = 1, cex.axis = 0.4)
}
# Window position labels
win_labels <- round(cand_windows$start_bp / 1e6, 1)
axis(1, at = seq(0, 1, length.out = n_win),
     labels = win_labels, las = 2, cex.axis = 0.5)
# Legend
legend("right", legend = all_subclusters, fill = pal,
       cex = 0.5, title = "Sub-cluster", inset = -0.15, xpd = TRUE)
dev.off()

message("[STEP41] Heatmap: ", fig_path)
message("\n[DONE] STEP41 membership trajectories for candidate ", cid)
message("  Verdict: ", h_verdict)
