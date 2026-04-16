#!/usr/bin/env Rscript

# =============================================================================
# STEP10d_layer_comparison_qc.R
#
# AUTOMATIC QC GATE — runs after STEP10c, compares Layer A / B / B2 / C.
#
# For each candidate, reads the distance matrices produced by STEP10c and asks:
#
#   1. GROUPING DISCORDANCE: Which windows are grouped together by Layer A
#      (PCA shape) but separated by Layer B2 (raw-vector induced distances)?
#      → Detects the "same shape, different fish" problem.
#
#   2. RESCUE: Which windows are weak/noisy in Layer A but form coherent
#      groups in Layer B2?
#      → Detects subtle inversion signal that PCA-shape misses.
#
#   3. ENCODING AGREEMENT: Do B2_minor, B2_major, and B2_012 agree on
#      window grouping? Where they disagree, is it discretization artifact
#      or real subgroup signal?
#
#   4. B2 vs C DELTA: Does B2 (full ordered vectors) find structure that
#      Layer C (scalar fingerprints) misses? Large delta = within-band
#      substructure or marker-order patterns that mean dosage erases.
#
#   5. FISH-STRUCTURE CONSISTENCY: For windows grouped by B2, do the
#      sample-belonging vectors (Layer B) also agree? If B2 groups windows
#      but B disagrees, it means the windows share fine-grained profile
#      structure but differ in coarse band membership — suspicious.
#
# OUTPUTS PER CANDIDATE:
#   candidate_layer_qc_report.tsv        — per-window-pair diagnostic table
#   candidate_layer_qc_summary.tsv       — candidate-level summary flags
#   candidate_layer_qc_figure.png        — composite diagnostic panel
#
# CROSS-CANDIDATE OUTPUT:
#   all_candidates_layer_qc_summary.tsv  — one-row-per-candidate for triage
#
# Usage:
#   Rscript STEP10d_layer_comparison_qc.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

has_ggplot   <- suppressWarnings(require(ggplot2, quietly = TRUE))
has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Load candidate table ─────────────────────────────────────────────────────
# Prefer snake-derived candidates (primary method) over 500kb control candidates.
SNAKE_CANDIDATE_TABLE <- file.path(INV_ROOT, "06_mds_candidates",
                                   "snake_regions_multiscale",
                                   "snake_candidate_regions.tsv.gz")

if (file.exists(SNAKE_CANDIDATE_TABLE)) {
  cand <- fread(SNAKE_CANDIDATE_TABLE)
  message("[STEP10d] Using SNAKE candidate table: ", SNAKE_CANDIDATE_TABLE)
} else if (file.exists(CANDIDATE_TABLE)) {
  cand <- fread(CANDIDATE_TABLE)
  message("[STEP10d] Using CONTROL candidate table (500kb): ", CANDIDATE_TABLE)
} else {
  stop("No candidate table found")
}
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
cand <- cand[order(chrom, start_bp)]
message("[STEP10d] Layer comparison QC for ", nrow(cand), " candidates")

# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

# Read a distance matrix TSV (window × window) from STEP10c
read_dist_mat <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  dt <- tryCatch(fread(filepath), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  # First col is "window", last col is "candidate_id"
  wcol <- "window"
  cidcol <- "candidate_id"
  mat_cols <- setdiff(names(dt), c(wcol, cidcol))
  if (length(mat_cols) < 2) return(NULL)
  labels <- dt[[wcol]]
  mat <- as.matrix(dt[, ..mat_cols])
  rownames(mat) <- labels; colnames(mat) <- labels
  storage.mode(mat) <- "double"
  mat
}

# Cluster a distance matrix, return group assignments at height h
cluster_at_h <- function(dmat, h = 0.5) {
  if (is.null(dmat) || nrow(dmat) < 2) return(NULL)
  d <- dmat
  med <- median(d, na.rm = TRUE)
  # v7.3 fix: guard against all-NA distance matrix
  if (is.na(med)) return(NULL)
  d[is.na(d)] <- med; diag(d) <- 0
  d <- (d + t(d)) / 2
  hc <- hclust(as.dist(d), method = "average")
  cutree(hc, h = h)
}

# Adjusted Rand Index between two clustering vectors
adj_rand <- function(a, b) {
  if (length(a) != length(b) || length(a) < 2) return(NA_real_)
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 2) return(NA_real_)
  a <- a[ok]; b <- b[ok]
  tab <- table(a, b); n <- sum(tab)
  if (n < 2) return(NA_real_)
  sc <- sum(choose(tab, 2))
  sa <- sum(choose(rowSums(tab), 2))
  sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2)
  mx <- 0.5 * (sa + sb)
  if (mx == e) return(1)
  (sc - e) / (mx - e)
}

# For two symmetric matrices, compute elementwise correlation (Mantel-like)
mat_cor <- function(m1, m2) {
  if (is.null(m1) || is.null(m2)) return(NA_real_)
  suppressWarnings(cor(as.vector(m1), as.vector(m2), use = "complete.obs"))
}

# For a pair of windows (i, j): are they "same group" in clustering g?
same_group <- function(g, i, j) {
  if (is.null(g)) return(NA)
  !is.na(g[i]) && !is.na(g[j]) && g[i] == g[j]
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════

all_summaries <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp)

  # Must match STEP10c's naming convention
  cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))

  cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
  plot_dir <- file.path(PLOTS_DIR, cand_prefix, "window_belonging")
  ensure_dir(cand_dir); ensure_dir(plot_dir)

  message("\n[STEP10d] ════════════════════════════════════════════════")
  message("[STEP10d] ", cand_prefix, " — layer comparison QC")

  # ─── Load distance matrices (names match STEP10c output) ─────────────
  dist_A          <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_A_snp.tsv")))
  dist_B          <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_B_belonging_hamming.tsv")))
  dist_B2_minor   <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_B2_rawvec_minor.tsv")))
  dist_B2_major   <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_B2_rawvec_major.tsv")))
  dist_B2_012     <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_B2_rawvec_012.tsv")))
  dist_C_minor    <- read_dist_mat(file.path(cand_dir, paste0(cand_prefix, "_dist_C_scalar_minor.tsv")))

  # Also read pre-computed agreement table
  agreement_file <- file.path(cand_dir, paste0(cand_prefix, "_layer_agreement.tsv"))
  agreement_dt <- if (file.exists(agreement_file)) fread(agreement_file) else NULL

  # Check minimum data
  if (is.null(dist_A) || is.null(dist_B2_minor)) {
    message("[SKIP] Missing distance matrices for candidate ", cid)
    next
  }
  nw <- nrow(dist_A)
  wlabels <- rownames(dist_A)
  if (nw < 3) { message("[SKIP] <3 windows"); next }

  message("[STEP10d] ", nw, " windows loaded")

  # ─── Cluster each layer at h=0.5 ────────────────────────────────────
  # h=0.5 means: windows with distance < 0.5 are "same system"
  grp_A        <- cluster_at_h(dist_A, 0.5)
  grp_B        <- cluster_at_h(dist_B, 0.5)
  grp_B2_minor <- cluster_at_h(dist_B2_minor, 0.5)
  grp_B2_major <- cluster_at_h(dist_B2_major, 0.5)
  grp_B2_012   <- cluster_at_h(dist_B2_012, 0.5)
  grp_C_minor  <- cluster_at_h(dist_C_minor, 0.5)

  # ═════════════════════════════════════════════════════════════════════
  # 1. PER-WINDOW-PAIR DIAGNOSTIC TABLE
  # ═════════════════════════════════════════════════════════════════════

  pairs <- list()
  for (i in seq_len(nw - 1)) {
    for (j in (i + 1):nw) {
      pairs[[length(pairs) + 1]] <- data.table(
        window_i = wlabels[i],
        window_j = wlabels[j],

        # Raw distances
        dist_A        = dist_A[i, j],
        dist_B        = if (!is.null(dist_B)) dist_B[i, j] else NA_real_,
        dist_B2_minor = dist_B2_minor[i, j],
        dist_B2_major = if (!is.null(dist_B2_major)) dist_B2_major[i, j] else NA_real_,
        dist_B2_012   = if (!is.null(dist_B2_012)) dist_B2_012[i, j] else NA_real_,
        dist_C_minor  = if (!is.null(dist_C_minor)) dist_C_minor[i, j] else NA_real_,

        # Same-group flags at h=0.5
        same_A        = same_group(grp_A, i, j),
        same_B        = same_group(grp_B, i, j),
        same_B2_minor = same_group(grp_B2_minor, i, j),
        same_B2_major = same_group(grp_B2_major, i, j),
        same_B2_012   = same_group(grp_B2_012, i, j),
        same_C_minor  = same_group(grp_C_minor, i, j)
      )
    }
  }
  pair_dt <- rbindlist(pairs)

  # ─── Derived flags ──────────────────────────────────────────────────
  # v7.3 fix: use explicit NA-safe logic — (x == TRUE) returns NA when x is NA,
  # so wrap with NA→FALSE coercion to ensure flags are always logical.
  na_false <- function(x) fifelse(is.na(x), FALSE, x)

  # DISCORDANCE: grouped by A but separated by B2
  pair_dt[, A_grouped_B2_split := na_false(same_A == TRUE & same_B2_minor == FALSE)]
  # RESCUE: separated by A but grouped by B2
  pair_dt[, A_split_B2_grouped := na_false(same_A == FALSE & same_B2_minor == TRUE)]
  # B2 encoding agreement (all three B2 sub-layers concur)
  pair_dt[, B2_encodings_agree := na_false(same_B2_minor == same_B2_major & same_B2_minor == same_B2_012)]
  # B2 vs C delta: grouped by B2 but split by C (within-band structure)
  pair_dt[, B2_grouped_C_split := na_false(same_B2_minor == TRUE & same_C_minor == FALSE)]
  # Fish-structure check: B2 groups but B disagrees
  pair_dt[, B2_grouped_B_split := na_false(same_B2_minor == TRUE & same_B == FALSE)]
  # Suspicious: shape agrees, coarse bands agree, but raw vectors disagree
  pair_dt[, subtle_divergence := na_false(same_A == TRUE & same_B == TRUE & same_B2_minor == FALSE)]

  pair_dt[, candidate_id := cid]
  fwrite(pair_dt, file.path(cand_dir, paste0(cand_prefix, "_layer_qc_report.tsv")), sep = "\t")

  # ═════════════════════════════════════════════════════════════════════
  # 2. CANDIDATE-LEVEL SUMMARY
  # ═════════════════════════════════════════════════════════════════════

  np <- nrow(pair_dt)

  # Cross-layer Mantel correlations
  mantel_A_B2    <- mat_cor(dist_A, dist_B2_minor)
  mantel_B_B2    <- mat_cor(dist_B, dist_B2_minor)
  mantel_B2m_B2M <- mat_cor(dist_B2_minor, dist_B2_major)
  mantel_B2m_B2d <- mat_cor(dist_B2_minor, dist_B2_012)
  mantel_B2_C    <- mat_cor(dist_B2_minor, dist_C_minor)

  # ARI between clusterings
  ari_A_B2    <- adj_rand(grp_A, grp_B2_minor)
  ari_B_B2    <- adj_rand(grp_B, grp_B2_minor)
  ari_B2m_B2M <- adj_rand(grp_B2_minor, grp_B2_major)
  ari_B2m_B2d <- adj_rand(grp_B2_minor, grp_B2_012)

  # Count diagnostics
  n_A_grouped_B2_split <- sum(pair_dt$A_grouped_B2_split, na.rm = TRUE)
  n_A_split_B2_grouped <- sum(pair_dt$A_split_B2_grouped, na.rm = TRUE)
  n_B2_enc_disagree    <- sum(!pair_dt$B2_encodings_agree, na.rm = TRUE)
  n_B2_grouped_C_split <- sum(pair_dt$B2_grouped_C_split, na.rm = TRUE)
  n_B2_grouped_B_split <- sum(pair_dt$B2_grouped_B_split, na.rm = TRUE)
  n_subtle_divergence  <- sum(pair_dt$subtle_divergence, na.rm = TRUE)

  # Number of distinct groups per layer
  n_groups_A   <- length(unique(grp_A))
  n_groups_B   <- if (!is.null(grp_B)) length(unique(grp_B)) else NA_integer_
  n_groups_B2m <- length(unique(grp_B2_minor))
  n_groups_B2M <- if (!is.null(grp_B2_major)) length(unique(grp_B2_major)) else NA_integer_
  n_groups_B2d <- if (!is.null(grp_B2_012)) length(unique(grp_B2_012)) else NA_integer_

  # ─── QC flags ────────────────────────────────────────────────────────
  #   FLAG_SHAPE_FISH_MISMATCH: A groups windows that B2 separates
  #   FLAG_B2_RESCUES:          B2 finds coherent groups that A misses
  #   FLAG_ENCODING_DISAGREE:   B2 encodings disagree on >20% of pairs
  #   FLAG_WITHIN_BAND_STRUCT:  B2 finds structure that scalar C misses
  #   FLAG_SUBTLE_DIVERGENCE:   A+B agree but B2 disagrees (hidden split)
  #   FLAG_B2_FINER:            B2 finds more groups than A

  flag_shape_fish     <- n_A_grouped_B2_split > 0
  flag_B2_rescues     <- n_A_split_B2_grouped > 0
  flag_enc_disagree   <- n_B2_enc_disagree > (np * 0.2)
  flag_within_band    <- n_B2_grouped_C_split > 0
  flag_subtle_div     <- n_subtle_divergence > 0
  flag_B2_finer       <- n_groups_B2m > n_groups_A

  summary_dt <- data.table(
    candidate_id = cid, chrom = chr, n_windows = nw, n_pairs = np,

    # Group counts
    n_groups_A = n_groups_A, n_groups_B = n_groups_B,
    n_groups_B2_minor = n_groups_B2m, n_groups_B2_major = n_groups_B2M,
    n_groups_B2_012 = n_groups_B2d,

    # Mantel correlations
    mantel_A_vs_B2 = round(mantel_A_B2, 4),
    mantel_B_vs_B2 = round(mantel_B_B2, 4),
    mantel_B2minor_vs_B2major = round(mantel_B2m_B2M, 4),
    mantel_B2minor_vs_B2_012 = round(mantel_B2m_B2d, 4),
    mantel_B2_vs_C = round(mantel_B2_C, 4),

    # ARI
    ari_A_vs_B2 = round(ari_A_B2, 4),
    ari_B_vs_B2 = round(ari_B_B2, 4),
    ari_B2minor_vs_B2major = round(ari_B2m_B2M, 4),
    ari_B2minor_vs_B2_012 = round(ari_B2m_B2d, 4),

    # Discordance counts
    n_A_grouped_B2_split = n_A_grouped_B2_split,
    n_A_split_B2_grouped = n_A_split_B2_grouped,
    n_B2_encoding_disagree = n_B2_enc_disagree,
    n_B2_grouped_C_split = n_B2_grouped_C_split,
    n_B2_grouped_B_split = n_B2_grouped_B_split,
    n_subtle_divergence = n_subtle_divergence,

    # QC flags
    FLAG_SHAPE_FISH_MISMATCH = flag_shape_fish,
    FLAG_B2_RESCUES          = flag_B2_rescues,
    FLAG_ENCODING_DISAGREE   = flag_enc_disagree,
    FLAG_WITHIN_BAND_STRUCT  = flag_within_band,
    FLAG_SUBTLE_DIVERGENCE   = flag_subtle_div,
    FLAG_B2_FINER_THAN_A     = flag_B2_finer
  )
  fwrite(summary_dt, file.path(cand_dir, paste0(cand_prefix, "_layer_qc_summary.tsv")), sep = "\t")

  # Diagnostic log
  message("[STEP10d] Groups: A=", n_groups_A,
          " B=", n_groups_B,
          " B2m=", n_groups_B2m,
          " B2M=", n_groups_B2M,
          " B2d=", n_groups_B2d)
  message("[STEP10d] Mantel A↔B2=", round(mantel_A_B2, 3),
          "  B↔B2=", round(mantel_B_B2, 3),
          "  B2m↔B2M=", round(mantel_B2m_B2M, 3),
          "  B2m↔B2d=", round(mantel_B2m_B2d, 3),
          "  B2↔C=", round(mantel_B2_C, 3))
  message("[STEP10d] ARI A↔B2=", round(ari_A_B2, 3),
          "  B↔B2=", round(ari_B_B2, 3),
          "  B2m↔B2M=", round(ari_B2m_B2M, 3),
          "  B2m↔B2d=", round(ari_B2m_B2d, 3))
  if (flag_shape_fish)
    message("[STEP10d] *** FLAG: ", n_A_grouped_B2_split,
            " window-pairs grouped by A but split by B2 (same shape, different fish)")
  if (flag_B2_rescues)
    message("[STEP10d] *** FLAG: ", n_A_split_B2_grouped,
            " window-pairs split by A but grouped by B2 (B2 rescues weak signal)")
  if (flag_enc_disagree)
    message("[STEP10d] *** FLAG: B2 encodings disagree on ", n_B2_enc_disagree,
            "/", np, " pairs (>20% — check discretization/polarity)")
  if (flag_subtle_div)
    message("[STEP10d] *** FLAG: ", n_subtle_divergence,
            " pairs where A+B agree but B2 disagrees (hidden substructure)")

  all_summaries[[length(all_summaries) + 1]] <- summary_dt

  # ═════════════════════════════════════════════════════════════════════
  # 3. COMPOSITE DIAGNOSTIC FIGURE
  # ═════════════════════════════════════════════════════════════════════

  if (has_ggplot && has_patchwork && nw >= 3) {
    message("[STEP10d] Generating diagnostic figure")

    # ── Panel A: Scatterplot of dist_A vs dist_B2_minor per pair ────
    pair_dt[, discordance_type := fifelse(
      A_grouped_B2_split, "A grouped / B2 split",
      fifelse(A_split_B2_grouped, "A split / B2 grouped",
        fifelse(subtle_divergence, "Subtle divergence", "Concordant"))
    )]

    p1 <- ggplot(pair_dt, aes(x = dist_A, y = dist_B2_minor,
                               color = discordance_type)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", alpha = 0.5) +
      scale_color_manual(values = c(
        "A grouped / B2 split"  = "#D73027",
        "A split / B2 grouped"  = "#4575B4",
        "Subtle divergence"     = "#F46D43",
        "Concordant"            = "grey60"
      )) +
      labs(title = "A: Layer A vs B2 distance",
           subtitle = paste0("Mantel r=", round(mantel_A_B2, 3),
                             "  ARI=", round(ari_A_B2, 3)),
           x = "Layer A (SNP shape)", y = "Layer B2 (raw-vector induced)",
           color = "Discordance") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7))

    # ── Panel B: B2 minor vs B2 major vs B2 012 ──────────────────────
    enc_dt <- rbind(
      pair_dt[, .(x = dist_B2_minor, y = dist_B2_major,
                  pair = paste(window_i, window_j), comp = "minor vs major")],
      pair_dt[, .(x = dist_B2_minor, y = dist_B2_012,
                  pair = paste(window_i, window_j), comp = "minor vs 012")]
    )

    p2 <- ggplot(enc_dt, aes(x = x, y = y)) +
      geom_point(size = 1.8, alpha = 0.6, color = "#2166AC") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ comp, scales = "free") +
      labs(title = "B: B2 encoding agreement",
           subtitle = paste0("minor↔major r=", round(mantel_B2m_B2M, 3),
                             "  minor↔012 r=", round(mantel_B2m_B2d, 3)),
           x = "B2 minor distance", y = "Comparison distance") +
      theme_bw(base_size = 10)

    # ── Panel C: B2 vs C delta (within-band structure detection) ─────
    if (!is.null(dist_C_minor)) {
      p3 <- ggplot(pair_dt, aes(x = dist_C_minor, y = dist_B2_minor)) +
        geom_point(size = 2, alpha = 0.6,
                   color = fifelse(pair_dt$B2_grouped_C_split %in% TRUE, "#D73027", "grey60")) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "C: B2 vs C (scalar fingerprint)",
             subtitle = paste0("r=", round(mantel_B2_C, 3),
                               "  B2-grouped-C-split: ", n_B2_grouped_C_split),
             x = "Layer C (scalar)", y = "Layer B2 (raw-vector)",
             caption = "Red = B2 groups but C splits") +
        theme_bw(base_size = 10)
    } else {
      p3 <- ggplot() + theme_void() +
        labs(title = "C: Layer C not available")
    }

    # ── Panel D: Group count comparison barplot ──────────────────────
    bar_dt <- data.table(
      layer = c("A (shape)", "B (bands)", "B2 minor", "B2 major", "B2 012"),
      n_groups = c(n_groups_A,
                   ifelse(is.na(n_groups_B), 0, n_groups_B),
                   n_groups_B2m,
                   ifelse(is.na(n_groups_B2M), 0, n_groups_B2M),
                   ifelse(is.na(n_groups_B2d), 0, n_groups_B2d))
    )
    bar_dt[, layer := factor(layer, levels = layer)]

    p4 <- ggplot(bar_dt, aes(x = layer, y = n_groups, fill = layer)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = n_groups), vjust = -0.3, size = 3.5) +
      scale_fill_manual(values = c("#4575B4", "#74ADD1", "#D73027",
                                    "#F46D43", "#FDAE61")) +
      labs(title = "D: Window groups per layer",
           subtitle = paste0(nw, " windows, h=0.5 cut"),
           x = NULL, y = "Groups") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # ── Panel E: Flag summary tile ───────────────────────────────────
    flag_dt <- data.table(
      flag = c("Shape/fish\nmismatch", "B2 rescues", "Encoding\ndisagree",
               "Within-band\nstructure", "Subtle\ndivergence", "B2 finer\nthan A"),
      value = c(flag_shape_fish, flag_B2_rescues, flag_enc_disagree,
                flag_within_band, flag_subtle_div, flag_B2_finer),
      count = c(n_A_grouped_B2_split, n_A_split_B2_grouped, n_B2_enc_disagree,
                n_B2_grouped_C_split, n_subtle_divergence,
                max(0, n_groups_B2m - n_groups_A))
    )
    flag_dt[, flag := factor(flag, levels = flag)]

    p5 <- ggplot(flag_dt, aes(x = flag, y = 1, fill = value)) +
      geom_tile(color = "white", linewidth = 1.5) +
      geom_text(aes(label = count), size = 4, fontface = "bold") +
      scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#A6DBA0"),
                        labels = c("TRUE" = "FLAGGED", "FALSE" = "OK")) +
      labs(title = "E: QC flags",
           subtitle = paste0("Candidate ", cid, " — ", nw, " windows"),
           x = NULL, y = NULL, fill = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "right")

    # ── Compose ──────────────────────────────────────────────────────
    fig <- (p1 | p2) / (p3 | p4) / p5 +
      patchwork::plot_annotation(
        title = paste0("STEP10d Layer Comparison QC — ", cand_prefix),
        subtitle = paste0(chr, " | ", nw, " windows | ",
                          sum(c(flag_shape_fish, flag_B2_rescues, flag_enc_disagree,
                                flag_within_band, flag_subtle_div, flag_B2_finer)),
                          " flags raised"),
        theme = theme(plot.title = element_text(face = "bold", size = 14))
      )

    figpath <- file.path(plot_dir, paste0(cand_prefix, "_layer_qc_figure.png"))
    ggsave(figpath, fig, width = 14, height = 14, dpi = 300)
    message("[STEP10d] Figure saved: ", figpath)
  }

  message("[STEP10d] ", cand_prefix, " — QC complete")
}

# ═══════════════════════════════════════════════════════════════════════════
# CROSS-CANDIDATE SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

if (length(all_summaries) > 0) {
  all_dt <- rbindlist(all_summaries, fill = TRUE)

  # Overall triage: sort by number of flags, then by discordance count
  flag_cols <- grep("^FLAG_", names(all_dt), value = TRUE)
  all_dt[, n_flags := rowSums(.SD == TRUE, na.rm = TRUE), .SDcols = flag_cols]
  all_dt[, total_discordance := n_A_grouped_B2_split + n_subtle_divergence]
  setorder(all_dt, -n_flags, -total_discordance)

  outpath <- file.path(FOLLOWUP_DIR, "all_candidates_layer_qc_summary.tsv")
  fwrite(all_dt, outpath, sep = "\t")
  message("\n[STEP10d] Cross-candidate summary: ", outpath)

  # Log triage
  message("[STEP10d] ════════════════════════════════════════════════")
  message("[STEP10d] TRIAGE SUMMARY (", nrow(all_dt), " candidates)")
  message("[STEP10d] ────────────────────────────────────────────────")
  for (ri in seq_len(min(nrow(all_dt), 20))) {
    r <- all_dt[ri]
    flags <- paste(flag_cols[as.logical(r[, ..flag_cols])], collapse = ", ")
    if (nchar(flags) == 0) flags <- "none"
    message(sprintf("[STEP10d]   cand %3d | %s | %d win | %d flags: %s",
                    r$candidate_id, r$chrom, r$n_windows, r$n_flags, flags))
  }
}

message("\n[DONE] STEP10d layer comparison QC complete")