#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01g_layer4_tube_graph.R  (v8.4)
#
# SAMPLE-AWARE TUBE GRAPH  --  Layer 4
#
# For each merged candidate (from C01d scoring + C01c triangles):
#   1. Choose anchor region (hottest subtriangle / strongest sub-regime)
#   2. Define anchor sample roles (LEFT / MIDDLE / RIGHT with confidence)
#   3. Track sample affinities across all windows in the candidate
#   4. Build tube graph: nodes = local regime pieces, edges = lane transitions
#   5. Stage D: annotate confounding / failure diagnosis
#   6. Stage E: annotate complex inversion / recombinant interpretation
#   7. Output labels + split suggestions
#
# The tube graph has TWO annotation layers on the same structure:
#   Layer D: why does the simple inversion model fail? (family, ROH, technical)
#   Layer E: what complex biology explains the signal? (recombinant, double-XO)
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>
#   --triangles <triangle_dir>     (intervals, subregimes, subtriangles, bridges)
#   --precomp <precomp_dir>
#   --samples <sample_list.tsv>
#   --relatedness <natora_pairs.tsv>  (real sample names, NOT Ind names)
#   [--snake2_dir <dir>]           (snake2 band assignments if available)
#   [--pruned_samples <list>]      (81 kin-pruned unrelated)
#   --outdir <dir>
#
# Outputs:
#   tube_nodes.tsv.gz              -- per-node: local regime + composition + confounding
#   tube_edges.tsv.gz              -- per-edge: lane continuity + transition type
#   tube_sample_tracks.tsv.gz      -- per-sample per-window: affinity + lane
#   tube_verdicts.tsv              -- per-candidate: Stage D + E labels
#   plots/<chr>_<id>_tube.png      -- tube graph visualization
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; precomp_dir <- NULL
samples_file <- NULL; relate_file <- NULL; snake2_dir <- NULL
pruned_file <- NULL; ghsl_dir <- NULL; outdir <- "tube_graph"; tier_max <- 2L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))      { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args))  { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp" && i < length(args))    { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args)) { relate_file <- args[i+1]; i <- i+2 }
  else if (a == "--snake2_dir" && i < length(args)) { snake2_dir <- args[i+1]; i <- i+2 }
  else if (a == "--ghsl_dir" && i < length(args))   { ghsl_dir <- args[i+1]; i <- i+2 }
  else if (a == "--pruned_samples" && i < length(args)) { pruned_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--tier_max" && i < length(args))   { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file)) stop("--scores required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 350
WINDOW_SMOOTH <- 3L  # rolling window for affinity smoothing

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01g] Loading data...")
cand_dt <- fread(scores_file)[tier <= tier_max][order(tier, -final_score)]

# Triangle outputs
iv_dt <- fread(file.path(triangle_dir, "triangle_intervals.tsv.gz"))
subreg_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_subregimes.tsv.gz")),
                       error = function(e) data.table())
subtri_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_subtriangles.tsv.gz")),
                       error = function(e) data.table())
comp_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_sample_composition.tsv.gz")),
                     error = function(e) data.table())

# Sample names (real = CGA, ind = Ind0)
real_names <- NULL
if (!is.null(samples_file) && file.exists(samples_file))
  real_names <- as.character(fread(samples_file, header = FALSE)[[1]])

# Relatedness (real sample names)
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
pruned_samples <- NULL
if (!is.null(pruned_file) && file.exists(pruned_file))
  pruned_samples <- as.character(fread(pruned_file, header = FALSE)[[1]])

precomp_cache <- list()

# Load GHSL data once (all chromosomes in one file)
ghsl_all <- NULL
if (!is.null(ghsl_dir)) {
  ghsl_f <- file.path(ghsl_dir, "snake3_sample_window.tsv.gz")
  if (file.exists(ghsl_f)) {
    ghsl_all <- tryCatch(fread(ghsl_f), error = function(e) NULL)
    if (!is.null(ghsl_all)) message("[C01g] GHSL loaded: ", nrow(ghsl_all), " sample-window records")
  }
}

message("[C01g] Candidates: ", nrow(cand_dt))

# =============================================================================
# STEP 1+2: ANCHOR REGION + SAMPLE ROLES
# =============================================================================

choose_anchor_and_roles <- function(dt, subreg, subtri, chr, iid, sample_names) {
  # Find the hottest sub-region or subtriangle as anchor
  hot_sub <- subreg[chrom == chr & interval_id == iid & sub_level == "hot"]
  if (nrow(hot_sub) == 0) hot_sub <- subreg[chrom == chr & interval_id == iid]
  if (nrow(hot_sub) == 0) return(NULL)

  # Pick the largest hot sub-region
  anchor <- hot_sub[order(-sub_n)][1]
  anchor_wins <- ((anchor$sub_L - 1) * 5L + 1):min(anchor$sub_R * 5L, nrow(dt))
  if (length(anchor_wins) < 5) return(NULL)

  # Build average PC1 loadings for anchor region
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20) return(NULL)

  anchor_mat <- as.matrix(dt[anchor_wins, ..available])
  anchor_avg <- colMeans(anchor_mat, na.rm = TRUE)
  valid <- is.finite(anchor_avg)
  if (sum(valid) < 20) return(NULL)

  vals <- anchor_avg[valid]; snames <- sub("^PC_1_", "", names(vals))

  # k=3 clustering for L/M/R roles

  km <- tryCatch(kmeans(vals, centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(km)) return(NULL)

  co <- order(km$centers[, 1])
  roles <- character(length(vals))
  for (bi in seq_along(co)) {
    role <- c("LEFT", "MIDDLE", "RIGHT")[bi]
    roles[km$cluster == co[bi]] <- role
  }

  # Confidence: distance from own center vs distance to nearest other center
  centers <- km$centers[co, 1]
  confidence <- numeric(length(vals))
  for (si in seq_along(vals)) {
    own_center <- centers[match(roles[si], c("LEFT", "MIDDLE", "RIGHT"))]
    other_centers <- centers[c("LEFT", "MIDDLE", "RIGHT") != roles[si]]
    dist_own <- abs(vals[si] - own_center)
    dist_other <- min(abs(vals[si] - other_centers))
    confidence[si] <- dist_other / max(dist_own + dist_other, 0.001)
  }

  # Core anchors: top 80% confidence
  core_threshold <- quantile(confidence, 0.2)

  data.table(
    sample = snames, role = roles, pc1 = round(vals, 4),
    confidence = round(confidence, 3),
    is_core_anchor = confidence >= core_threshold,
    anchor_sub_L = anchor$sub_L, anchor_sub_R = anchor$sub_R
  )
}

# =============================================================================
# STEP 3: SAMPLE AFFINITY TRACKING ACROSS WINDOWS
# =============================================================================

track_sample_affinities <- function(dt, anchors, sample_names, win_start, win_end) {
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20) return(data.table())

  anchor_left <- anchors[role == "LEFT" & is_core_anchor == TRUE]$sample
  anchor_mid <- anchors[role == "MIDDLE" & is_core_anchor == TRUE]$sample
  anchor_right <- anchors[role == "RIGHT" & is_core_anchor == TRUE]$sample

  if (length(anchor_left) < 3 || length(anchor_mid) < 3 || length(anchor_right) < 3)
    return(data.table())

  # For each window, compute each sample's affinity to each anchor group
  results <- list()
  for (wi in seq(win_start, min(win_end, nrow(dt)), by = WINDOW_SMOOTH)) {
    w_range <- wi:min(wi + WINDOW_SMOOTH - 1, nrow(dt))
    w_mat <- as.matrix(dt[w_range, ..available])
    w_avg <- colMeans(w_mat, na.rm = TRUE)
    valid <- is.finite(w_avg)
    if (sum(valid) < 20) next

    w_vals <- w_avg[valid]; w_names <- sub("^PC_1_", "", names(w_vals))

    # Compute mean distance to each anchor group
    left_center <- mean(w_vals[w_names %in% anchor_left], na.rm = TRUE)
    mid_center <- mean(w_vals[w_names %in% anchor_mid], na.rm = TRUE)
    right_center <- mean(w_vals[w_names %in% anchor_right], na.rm = TRUE)

    if (!is.finite(left_center) || !is.finite(mid_center) || !is.finite(right_center)) next

    for (si in seq_along(w_vals)) {
      d_left <- abs(w_vals[si] - left_center)
      d_mid <- abs(w_vals[si] - mid_center)
      d_right <- abs(w_vals[si] - right_center)
      total <- d_left + d_mid + d_right
      if (total < 0.001) total <- 0.001

      # Soft probabilities (inverse distance)
      p_left <- (1/d_left) / (1/d_left + 1/d_mid + 1/d_right)
      p_mid <- (1/d_mid) / (1/d_left + 1/d_mid + 1/d_right)
      p_right <- (1/d_right) / (1/d_left + 1/d_mid + 1/d_right)
      # Handle zero distance
      if (!is.finite(p_left)) { p_left <- 1; p_mid <- 0; p_right <- 0 }
      if (!is.finite(p_mid)) { p_left <- 0; p_mid <- 1; p_right <- 0 }
      if (!is.finite(p_right)) { p_left <- 0; p_mid <- 0; p_right <- 1 }

      entropy <- -sum(c(p_left, p_mid, p_right) *
                        log2(pmax(c(p_left, p_mid, p_right), 1e-10)))
      hard_lane <- c("LEFT", "MIDDLE", "RIGHT")[which.max(c(p_left, p_mid, p_right))]

      results[[length(results)+1]] <- data.table(
        window = wi, sample = w_names[si],
        p_left = round(p_left, 3), p_mid = round(p_mid, 3), p_right = round(p_right, 3),
        entropy = round(entropy, 3), hard_lane = hard_lane,
        pos_mb = round((dt$start_bp[wi] + dt$end_bp[min(wi+2, nrow(dt))]) / 2e6, 3)
      )
    }
  }
  if (length(results) > 0) rbindlist(results) else data.table()
}

# =============================================================================
# STEP 4: WINDOW-LEVEL SUMMARIES
# =============================================================================

compute_window_summaries <- function(tracks, anchors) {
  if (nrow(tracks) == 0) return(data.table())

  anchor_left <- anchors[role == "LEFT" & is_core_anchor == TRUE]$sample
  anchor_mid <- anchors[role == "MIDDLE" & is_core_anchor == TRUE]$sample
  anchor_right <- anchors[role == "RIGHT" & is_core_anchor == TRUE]$sample

  tracks[, .(
    anchor_left_retention = mean(hard_lane[sample %in% anchor_left] == "LEFT", na.rm = TRUE),
    anchor_mid_retention = mean(hard_lane[sample %in% anchor_mid] == "MIDDLE", na.rm = TRUE),
    anchor_right_retention = mean(hard_lane[sample %in% anchor_right] == "RIGHT", na.rm = TRUE),
    mid_entropy = mean(entropy[sample %in% anchor_mid], na.rm = TRUE),
    overall_entropy = mean(entropy, na.rm = TRUE),
    mid_centrality = mean(p_mid[sample %in% anchor_mid], na.rm = TRUE),
    lr_asymmetry = abs(mean(p_left[sample %in% anchor_left], na.rm = TRUE) -
                        mean(p_right[sample %in% anchor_right], na.rm = TRUE)),
    pos_mb = pos_mb[1]
  ), by = window]
}

# =============================================================================
# GHSL HAPLOTYPE CONTRAST (Snake 3 support)
# =============================================================================
# For each sample in a candidate region, compute per-window haplotype contrast
# using Clair3 phased genotypes. The contrast = fraction of het sites where
# the two haplotypes carry different alleles relative to the local PCA grouping.
#
# For recombinant detection: if a HET sample's contrast DROPS in one sub-region
# but stays high in the rest, that sub-region has a recombinant haplotype where
# one copy switched from inversion to reference arrangement.
#
# Input: per-chr GT matrix from Clair3 (samples x positions), phased as 0|1 or 1|0
# The ghsl_dir should contain per-chr files: snake3_sample_window.tsv.gz
# with columns: window, sample, contrast, n_het_sites, phase_block_bp

compute_ghsl_per_node <- function(ghsl_dt, node_wins, anchors, chr) {
  if (is.null(ghsl_dt) || nrow(ghsl_dt) == 0) return(list(
    ghsl_mid_contrast = NA, ghsl_left_contrast = NA, ghsl_right_contrast = NA,
    ghsl_recomb_signal = NA, ghsl_available = FALSE
  ))

  # Filter to windows in this node
  node_ghsl <- ghsl_dt[global_window_id %in% node_wins]
  if (nrow(node_ghsl) == 0) return(list(
    ghsl_mid_contrast = NA, ghsl_left_contrast = NA, ghsl_right_contrast = NA,
    ghsl_recomb_signal = NA, ghsl_available = FALSE
  ))

  mid_samples <- anchors[role == "MIDDLE" & is_core_anchor == TRUE]$sample
  left_samples <- anchors[role == "LEFT" & is_core_anchor == TRUE]$sample
  right_samples <- anchors[role == "RIGHT" & is_core_anchor == TRUE]$sample

  ghsl_mid <- mean(node_ghsl[sample_id %in% mid_samples]$ghsl_total_weighted, na.rm = TRUE)
  ghsl_left <- mean(node_ghsl[sample_id %in% left_samples]$ghsl_total_weighted, na.rm = TRUE)
  ghsl_right <- mean(node_ghsl[sample_id %in% right_samples]$ghsl_total_weighted, na.rm = TRUE)

  # Recombinant signal: MIDDLE samples should have HIGH contrast (heterozygous)
  # LEFT/RIGHT should have LOW contrast (homozygous for ref or inv)
  # If MIDDLE contrast drops, those samples may be recombinant
  expected_mid_high <- ghsl_mid > ghsl_left && ghsl_mid > ghsl_right
  recomb_signal <- if (!is.finite(ghsl_mid)) NA
                   else if (expected_mid_high) 0  # normal: HET has high contrast
                   else ghsl_left + ghsl_right - 2 * ghsl_mid  # positive = recombinant-like

  list(ghsl_mid_contrast = round(ghsl_mid, 4),
       ghsl_left_contrast = round(ghsl_left, 4),
       ghsl_right_contrast = round(ghsl_right, 4),
       ghsl_recomb_signal = round(recomb_signal, 4),
       ghsl_available = TRUE)
}

# GHSL edge comparison: detect contrast changes between consecutive nodes
# A contrast DROP for middle samples = recombinant boundary
# A contrast SWAP (left becomes high, right stays low) = double-crossover boundary
compute_ghsl_edge <- function(ghsl_n1, ghsl_n2) {
  if (!ghsl_n1$ghsl_available || !ghsl_n2$ghsl_available)
    return(list(ghsl_mid_delta = NA, ghsl_left_delta = NA, ghsl_right_delta = NA,
                ghsl_edge_type = "no_ghsl"))

  mid_delta <- ghsl_n2$ghsl_mid_contrast - ghsl_n1$ghsl_mid_contrast
  left_delta <- ghsl_n2$ghsl_left_contrast - ghsl_n1$ghsl_left_contrast
  right_delta <- ghsl_n2$ghsl_right_contrast - ghsl_n1$ghsl_right_contrast

  # Classify the GHSL transition
  ghsl_type <- if (abs(mid_delta) < 0.05 && abs(left_delta) < 0.05 && abs(right_delta) < 0.05) {
    "ghsl_stable"  # no change in haplotype contrast
  } else if (mid_delta < -0.1 && abs(left_delta) < 0.05) {
    "ghsl_mid_drop"  # middle loses contrast = recombinant
  } else if (left_delta > 0.1 && right_delta < -0.05) {
    "ghsl_side_swap"  # left gains contrast, right loses = double-XO
  } else if (right_delta > 0.1 && left_delta < -0.05) {
    "ghsl_side_swap"
  } else if (abs(mid_delta) > 0.1) {
    "ghsl_mid_change"  # middle shifts but unclear direction
  } else {
    "ghsl_gradual"
  }

  list(ghsl_mid_delta = round(mid_delta, 4),
       ghsl_left_delta = round(left_delta, 4),
       ghsl_right_delta = round(right_delta, 4),
       ghsl_edge_type = ghsl_type)
}

# =============================================================================
# STEP 5+6: BUILD TUBE GRAPH (nodes + edges + two annotation layers)
# =============================================================================

build_tube_graph <- function(win_summaries, subreg, anchors, relate_dt, chr, iid,
                              real_names, ind_to_real, ghsl_dt = NULL) {
  if (nrow(win_summaries) == 0) return(list(nodes = data.table(), edges = data.table()))

  # Nodes = sub-regime pieces
  iv_sub <- subreg[chrom == chr & interval_id == iid]
  if (nrow(iv_sub) == 0) {
    # Single node for the whole interval
    iv_sub <- data.table(sub_L = min(win_summaries$window) %/% 5 + 1,
                          sub_R = max(win_summaries$window) %/% 5 + 1,
                          sub_level = "whole", sub_inside = NA, sub_n = nrow(win_summaries),
                          chrom = chr, interval_id = iid)
  }

  nodes <- list()
  for (ni in seq_len(nrow(iv_sub))) {
    sr <- iv_sub[ni]
    node_wins <- win_summaries[window >= (sr$sub_L - 1) * 5 + 1 &
                                window <= sr$sub_R * 5]
    if (nrow(node_wins) == 0) next

    # Node metrics
    mid_ret <- mean(node_wins$anchor_mid_retention, na.rm = TRUE)
    left_ret <- mean(node_wins$anchor_left_retention, na.rm = TRUE)
    right_ret <- mean(node_wins$anchor_right_retention, na.rm = TRUE)
    mean_entropy <- mean(node_wins$overall_entropy, na.rm = TRUE)

    # Relatedness within anchor-middle samples (confounding check)
    mid_relatedness <- NA_real_
    if (!is.null(relate_dt) && nrow(anchors) > 0) {
      mid_samples <- anchors[role == "MIDDLE" & is_core_anchor == TRUE]$sample
      # Map Ind -> real if needed
      if (!is.null(ind_to_real)) {
        mid_real <- ind_to_real[mid_samples]
        mid_real <- mid_real[!is.na(mid_real)]
      } else {
        mid_real <- mid_samples
      }
      mid_pairs <- relate_dt[sample1 %in% mid_real & sample2 %in% mid_real]
      if (nrow(mid_pairs) > 0) mid_relatedness <- mean(mid_pairs$theta, na.rm = TRUE)
    }

    # GHSL haplotype contrast for this node
    node_win_ids <- node_wins$window
    ghsl_node <- compute_ghsl_per_node(ghsl_dt, node_win_ids, anchors, chr)

    nodes[[ni]] <- data.table(
      chrom = chr, interval_id = iid, node_id = ni,
      sub_L = sr$sub_L, sub_R = sr$sub_R, sub_level = sr$sub_level,
      n_windows = nrow(node_wins),
      anchor_mid_retention = round(mid_ret, 3),
      anchor_left_retention = round(left_ret, 3),
      anchor_right_retention = round(right_ret, 3),
      mean_entropy = round(mean_entropy, 3),
      mid_relatedness = round(mid_relatedness, 5),
      ghsl_mid_contrast = ghsl_node$ghsl_mid_contrast,
      ghsl_left_contrast = ghsl_node$ghsl_left_contrast,
      ghsl_right_contrast = ghsl_node$ghsl_right_contrast,
      ghsl_recomb_signal = ghsl_node$ghsl_recomb_signal,
      ghsl_available = ghsl_node$ghsl_available,
      # Stage D: confounding label
      confounding_label = fifelse(
        !is.na(mid_relatedness) & mid_relatedness > 0.1, "family_enriched",
        fifelse(mean_entropy > 1.2, "diffuse_noisy",
        fifelse(mid_ret < 0.5, "broken_middle", "clean"))),
      # Stage E: regime label (preliminary, updated by edges)
      regime_label = sr$sub_level
    )
  }
  node_dt <- if (length(nodes) > 0) rbindlist(nodes) else data.table()

  # Edges = transitions between consecutive nodes
  edges <- list()
  if (nrow(node_dt) >= 2) {
    for (ei in seq_len(nrow(node_dt) - 1)) {
      n1 <- node_dt[ei]; n2 <- node_dt[ei + 1]

      # Lane overlap: do the same samples keep the same roles?
      w1 <- win_summaries[window >= (n1$sub_L - 1) * 5 + 1 & window <= n1$sub_R * 5]
      w2 <- win_summaries[window >= (n2$sub_L - 1) * 5 + 1 & window <= n2$sub_R * 5]

      # Retention comparison
      mid_shift <- abs(n1$anchor_mid_retention - n2$anchor_mid_retention)
      left_shift <- abs(n1$anchor_left_retention - n2$anchor_left_retention)
      right_shift <- abs(n1$anchor_right_retention - n2$anchor_right_retention)
      total_shift <- mid_shift + left_shift + right_shift

      # Asymmetry change
      asym_change <- if (nrow(w1) > 0 && nrow(w2) > 0) {
        abs(mean(w1$lr_asymmetry, na.rm = TRUE) - mean(w2$lr_asymmetry, na.rm = TRUE))
      } else 0

      # Classify edge transition
      edge_type <- if (total_shift < 0.15) {
        "same_tube"
      } else if (total_shift < 0.15 && n2$mean_entropy < n1$mean_entropy) {
        "same_tube_cleaner"
      } else if (mid_shift > 0.3 && left_shift < 0.15 && right_shift < 0.15) {
        "middle_switch"  # middle group changes but sides stay
      } else if (left_shift > 0.3 || right_shift > 0.3) {
        "side_switch"  # one side reorganizes
      } else if (mid_shift > 0.2 && asym_change > 0.2) {
        "lane_switch"
      } else if (total_shift > 0.5) {
        "hard_break"
      } else {
        "gradual_transition"
      }

      # Stage D edge confounding
      d_label <- "clean"
      if (n1$confounding_label == "family_enriched" || n2$confounding_label == "family_enriched")
        d_label <- "confounded_by_family"
      if (n1$confounding_label == "diffuse_noisy" || n2$confounding_label == "diffuse_noisy")
        d_label <- "noisy_transition"

      # GHSL edge comparison (orthogonal haplotype-contrast evidence)
      ghsl_n1 <- list(ghsl_mid_contrast = n1$ghsl_mid_contrast,
                       ghsl_left_contrast = n1$ghsl_left_contrast,
                       ghsl_right_contrast = n1$ghsl_right_contrast,
                       ghsl_available = n1$ghsl_available)
      ghsl_n2 <- list(ghsl_mid_contrast = n2$ghsl_mid_contrast,
                       ghsl_left_contrast = n2$ghsl_left_contrast,
                       ghsl_right_contrast = n2$ghsl_right_contrast,
                       ghsl_available = n2$ghsl_available)
      ghsl_edge <- compute_ghsl_edge(ghsl_n1, ghsl_n2)

      # Stage E: combine PCA-based edge_type with GHSL evidence
      # GHSL confirmation strengthens the call; GHSL contradiction weakens it
      e_label <- if (edge_type == "same_tube" || edge_type == "same_tube_cleaner") {
        if (ghsl_edge$ghsl_edge_type == "ghsl_stable") "continuous_inversion"
        else if (ghsl_edge$ghsl_edge_type == "ghsl_mid_drop") "recombinant_ghsl_confirmed"
        else "continuous_inversion"
      } else if (edge_type == "middle_switch") {
        if (ghsl_edge$ghsl_edge_type == "ghsl_mid_drop") "recombinant_confirmed"
        else if (ghsl_edge$ghsl_edge_type == "ghsl_mid_change") "possible_recombinant"
        else "possible_recombinant"
      } else if (edge_type == "side_switch") {
        if (ghsl_edge$ghsl_edge_type == "ghsl_side_swap") "double_crossover_confirmed"
        else "possible_double_crossover"
      } else if (edge_type == "lane_switch") {
        "regime_change"
      } else if (edge_type == "hard_break") {
        "system_boundary"
      } else {
        "gradual"
      }

      edges[[ei]] <- data.table(
        chrom = chr, interval_id = iid,
        from_node = ei, to_node = ei + 1,
        mid_shift = round(mid_shift, 3),
        left_shift = round(left_shift, 3),
        right_shift = round(right_shift, 3),
        total_shift = round(total_shift, 3),
        asym_change = round(asym_change, 3),
        edge_type = edge_type,
        ghsl_edge_type = ghsl_edge$ghsl_edge_type,
        ghsl_mid_delta = ghsl_edge$ghsl_mid_delta,
        ghsl_left_delta = ghsl_edge$ghsl_left_delta,
        ghsl_right_delta = ghsl_edge$ghsl_right_delta,
        stage_d_label = d_label,
        stage_e_label = e_label
      )
    }
  }
  edge_dt <- if (length(edges) > 0) rbindlist(edges) else data.table()

  list(nodes = node_dt, edges = edge_dt)
}

# =============================================================================
# STEP 7+8: VERDICTS + SPLIT SUGGESTIONS
# =============================================================================

compute_verdict <- function(node_dt, edge_dt) {
  if (nrow(node_dt) == 0) return(list(verdict = "no_data", split = "no_data"))

  # Stage D: overall confounding assessment
  n_family <- sum(node_dt$confounding_label == "family_enriched")
  n_broken <- sum(node_dt$confounding_label == "broken_middle")
  n_clean <- sum(node_dt$confounding_label == "clean")
  frac_clean <- n_clean / nrow(node_dt)

  stage_d <- if (n_family > nrow(node_dt) / 2) "family_dominated"
             else if (n_broken > nrow(node_dt) / 2) "unstable_middle"
             else if (frac_clean > 0.7) "clean"
             else "mixed_confounding"

  # Stage E: overall biological interpretation
  if (nrow(edge_dt) == 0) {
    stage_e <- "single_regime"
  } else {
    n_recomb <- sum(edge_dt$stage_e_label %in% c("possible_recombinant", "recombinant_confirmed",
                                                    "recombinant_ghsl_confirmed"))
    n_dxo <- sum(edge_dt$stage_e_label %in% c("possible_double_crossover", "double_crossover_confirmed"))
    n_regime <- sum(edge_dt$stage_e_label == "regime_change")
    n_boundary <- sum(edge_dt$stage_e_label == "system_boundary")
    n_continuous <- sum(edge_dt$stage_e_label == "continuous_inversion")
    n_ghsl_confirmed <- sum(edge_dt$stage_e_label %in% c("recombinant_confirmed",
                                                           "double_crossover_confirmed",
                                                           "recombinant_ghsl_confirmed"))

    stage_e <- if (n_continuous == nrow(edge_dt)) "broad_clean_inversion"
               else if (n_continuous >= nrow(edge_dt) * 0.7) "hot_core_inversion"
               else if (n_ghsl_confirmed > 0 && n_dxo > 0) "double_crossover_ghsl_confirmed"
               else if (n_ghsl_confirmed > 0 && n_recomb > 0) "recombinant_ghsl_confirmed"
               else if (n_recomb > 0 && n_dxo == 0) "simple_recombinant"
               else if (n_dxo > 0) "double_crossover_like"
               else if (n_regime > 0) "nested_composite"
               else if (n_boundary > 0) "adjacent_systems"
               else "complex_mosaic"
  }

  # Split suggestion (conservative: need convergent evidence)
  split <- "no_split"
  if (nrow(edge_dt) > 0) {
    hard_breaks <- edge_dt[edge_type == "hard_break"]
    regime_changes <- edge_dt[stage_e_label == "regime_change"]
    if (nrow(hard_breaks) > 0 && nrow(regime_changes) > 0) {
      split <- "split_at_hard_break"
    } else if (nrow(hard_breaks) > 0) {
      split <- "inspect_manually"
    } else if (nrow(regime_changes) > 0 && stage_d == "clean") {
      split <- "split_at_regime_change"
    }
  }

  list(stage_d = stage_d, stage_e = stage_e, split = split)
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_nodes <- list(); all_edges <- list(); all_tracks <- list(); all_verdicts <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]
  chr <- cand$chrom; iid <- cand$interval_id
  message("\n[C01g] ", ci, "/", nrow(cand_dt), ": ", chr, " I", iid,
          " (", cand$start_mb, "-", cand$end_mb, " Mb)")

  # Load precomp
  if (is.null(precomp_cache[[chr]])) {
    f <- list.files(precomp_dir, pattern = paste0(chr, "\\.precomp\\.rds$"), full.names = TRUE)
    if (length(f) > 0) precomp_cache[[chr]] <- readRDS(f[1])
  }
  pc <- precomp_cache[[chr]]
  if (is.null(pc)) { message("  No precomp"); next }
  dt <- pc$dt

  # Filter GHSL data for this chromosome
  ghsl_chr <- NULL
  if (!is.null(ghsl_all) && "chrom" %in% names(ghsl_all)) {
    ghsl_chr <- ghsl_all[chrom == chr]
    if (nrow(ghsl_chr) > 0) message("  GHSL: ", nrow(ghsl_chr), " records for ", chr)
    else ghsl_chr <- NULL
  }

  # Sample names (Ind-space for PC columns)
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  sample_names <- sub("^PC_1_", "", pc1_cols)
  ind_to_real <- NULL
  if (!is.null(real_names) && length(real_names) == length(sample_names) &&
      grepl("^Ind[0-9]", sample_names[1])) {
    ind_to_real <- setNames(real_names, sample_names)
  }

  # Step 1+2: Anchor
  anchors <- choose_anchor_and_roles(dt, subreg_dt, subtri_dt, chr, iid, sample_names)
  if (is.null(anchors)) { message("  No anchor"); next }
  n_left <- sum(anchors$role == "LEFT")
  n_mid <- sum(anchors$role == "MIDDLE")
  n_right <- sum(anchors$role == "RIGHT")
  message("  Anchor: L=", n_left, " M=", n_mid, " R=", n_right,
          " (core: ", sum(anchors$is_core_anchor), ")")

  # Step 3: Track affinities
  start_bp <- cand$start_mb * 1e6; end_bp <- cand$end_mb * 1e6
  win_start <- max(1, which.min(abs(dt$start_bp - start_bp)))
  win_end <- min(nrow(dt), which.min(abs(dt$end_bp - end_bp)))
  tracks <- track_sample_affinities(dt, anchors, sample_names, win_start, win_end)
  message("  Tracks: ", nrow(tracks), " sample-window records")

  # Step 4: Window summaries
  win_sum <- compute_window_summaries(tracks, anchors)

  # Step 5+6: Tube graph
  tg <- build_tube_graph(win_sum, subreg_dt, anchors, relate_dt, chr, iid,
                          real_names, ind_to_real, ghsl_dt = ghsl_chr)
  message("  Graph: ", nrow(tg$nodes), " nodes, ", nrow(tg$edges), " edges")

  # Step 7+8: Verdict
  vrd <- compute_verdict(tg$nodes, tg$edges)
  message("  Stage D: ", vrd$stage_d, " | Stage E: ", vrd$stage_e,
          " | Split: ", vrd$split)

  # Collect
  if (nrow(tg$nodes) > 0) all_nodes[[length(all_nodes)+1]] <- tg$nodes
  if (nrow(tg$edges) > 0) all_edges[[length(all_edges)+1]] <- tg$edges
  if (nrow(tracks) > 0 && nrow(tracks) < 500000)
    all_tracks[[length(all_tracks)+1]] <- tracks[, .(chrom = chr, interval_id = iid,
                                                      window, sample, p_left, p_mid, p_right,
                                                      entropy, hard_lane, pos_mb)]

  all_verdicts[[length(all_verdicts)+1]] <- data.table(
    chrom = chr, interval_id = iid,
    start_mb = cand$start_mb, end_mb = cand$end_mb,
    tier = cand$tier, pattern = cand$pattern,
    stage_d = vrd$stage_d, stage_e = vrd$stage_e,
    split_suggestion = vrd$split,
    n_nodes = nrow(tg$nodes), n_edges = nrow(tg$edges),
    n_clean_nodes = sum(tg$nodes$confounding_label == "clean"),
    n_family_nodes = sum(tg$nodes$confounding_label == "family_enriched"),
    n_continuous_edges = if (nrow(tg$edges) > 0) sum(tg$edges$stage_e_label == "continuous_inversion") else 0,
    n_recombinant_edges = if (nrow(tg$edges) > 0) sum(tg$edges$stage_e_label == "possible_recombinant") else 0,
    n_dxo_edges = if (nrow(tg$edges) > 0) sum(tg$edges$stage_e_label == "possible_double_crossover") else 0
  )
}

# =============================================================================
# WRITE
# =============================================================================

message("\n[C01g] Writing...")
nd <- if (length(all_nodes) > 0) rbindlist(all_nodes, fill = TRUE) else data.table()
ed <- if (length(all_edges) > 0) rbindlist(all_edges, fill = TRUE) else data.table()
td <- if (length(all_tracks) > 0) rbindlist(all_tracks, fill = TRUE) else data.table()
vd <- if (length(all_verdicts) > 0) rbindlist(all_verdicts) else data.table()

fwrite(nd, file.path(outdir, "tube_nodes.tsv.gz"), sep = "\t")
fwrite(ed, file.path(outdir, "tube_edges.tsv.gz"), sep = "\t")
fwrite(td, file.path(outdir, "tube_sample_tracks.tsv.gz"), sep = "\t")
fwrite(vd, file.path(outdir, "tube_verdicts.tsv"), sep = "\t")

message("[C01g] Nodes: ", nrow(nd), " | Edges: ", nrow(ed),
        " | Tracks: ", nrow(td), " | Verdicts: ", nrow(vd))

if (nrow(vd) > 0) {
  message("\n[C01g] === STAGE D SUMMARY ===")
  for (v in sort(unique(vd$stage_d))) message("  ", v, ": ", sum(vd$stage_d == v))
  message("\n[C01g] === STAGE E SUMMARY ===")
  for (v in sort(unique(vd$stage_e))) message("  ", v, ": ", sum(vd$stage_e == v))
}

message("\n[DONE] -> ", outdir)
