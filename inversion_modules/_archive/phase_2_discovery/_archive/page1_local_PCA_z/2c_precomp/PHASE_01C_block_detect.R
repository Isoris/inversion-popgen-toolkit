#!/usr/bin/env Rscript

# =============================================================================
# PHASE_01C_block_detect.R  (v8.5)
#
# UNIFIED: block detection + boundary classification + blue-cross diagnosis +
# block concordance matrix + sample grouping per block.
#
# Merges v8.4 scripts: C01a3 (blue cross) + C01a4 (block comparison) +
# new boundary classification.
#
# Outputs consumed by PHASE_02 (snake cores + merge + triangles):
#   landscape/block_registry_<chr>.tsv.gz      — blocks with type + metrics
#   landscape/boundary_catalog_<chr>.tsv.gz    — classified boundaries
#   landscape/blue_cross_verdicts_<chr>.tsv.gz — assembly error vs real
#   landscape/block_concordance_<chr>.tsv.gz   — block-pair concordance
#   landscape/01C_window_pa.tsv.gz             — per-window block membership
#   landscape/plots/<chr>_L1_local_sim_blocks.png
#   landscape/plots/<chr>_L2_block_facets.png
#   landscape/plots/<chr>_L3_concordance.png
#   landscape/plots/<chr>_L4_boundary_catalog.png
#
# Usage:
#   Rscript PHASE_01C_block_detect.R <precomp_dir> <outdir> \
#     [--mode hatchery|wild] [--gaps gaps.bed] [--agp scaffold.agp] \
#     [--chrom C_gar_LG28] [--local_range 80]
#
# =============================================================================
# REGISTRY_CONTRACT
#   BLOCKS_WRITTEN:
#     - block_detect: registries/schemas/structured_block_schemas/block_detect.schema.json
#       keys: (none — not currently written from this script)
#       status: BLOCKED_ON_NO_CANDIDATE_JOIN
#       note: PHASE_01C produces per-chromosome blocks (block_registry_<chr>.tsv.gz)
#             BEFORE candidate_ids exist. The schema asks for per-candidate
#             fields; candidates are created by C01d downstream. Proper wiring
#             belongs in C01d (or a new C00b attribution step), where the join
#             block_id → candidate_id is available (iv_dt carries block_id as
#             interval_id, see C01d L174).
#
#             Secondary blocker: even if the join existed, 5 of the 9 schema
#             fields aren't emitted by this script. block_registry_<chr>.tsv.gz
#             columns are: chrom, block_id, start_bp, end_bp, start_mb, end_mb,
#             start_global_window_id, end_global_window_id, n_windows, mean_sim,
#             has_3bands (see script L566-574). Schema fields NOT computed here:
#               - confidence_score       (no composite score exists per block)
#               - inner_boundary_left_bp  (blue-cross inner edges live in the
#                                          boundary_catalog output, not the
#                                          block_registry)
#               - inner_boundary_right_bp (same)
#               - mean_outside_similarity (not computed — would need to scan
#                                          local_sim outside block windows)
#               - contrast_ratio          (not computed)
#
#             Two ways forward, in order of preference:
#               (a) extend this script to compute the 5 missing fields and
#                   add them to block_reg_rows; then add a register_block_detect
#                   helper that C01d calls per-candidate after looking up the
#                   block row by interval_id.
#               (b) drop the 5 missing fields from the schema (weaker — the
#                   inner_boundary fields are the real scientific value here).
#   KEYS_IN: none
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript PHASE_01C_block_detect.R <precomp_dir> <outdir> [opts]")

precomp_dir <- args[1]
outdir      <- args[2]
pop_mode    <- "hatchery"
gaps_file   <- NULL; agp_file <- NULL; chrom_filter <- NULL
LOCAL_RANGE <- 80L; MIN_BLOCK <- 8L; TOP_BLOCKS <- 25L
FLANK_W     <- 15L

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--mode" && i < length(args))         { pop_mode <- args[i+1]; i <- i+2L }
  else if (a == "--gaps" && i < length(args))    { gaps_file <- args[i+1]; i <- i+2L }
  else if (a == "--agp" && i < length(args))     { agp_file <- args[i+1]; i <- i+2L }
  else if (a == "--chrom" && i < length(args))   { chrom_filter <- args[i+1]; i <- i+2L }
  else if (a == "--local_range" && i < length(args)) { LOCAL_RANGE <- as.integer(args[i+1]); i <- i+2L }
  else { i <- i+1L }
}

dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

# Mode-dependent thresholds
if (pop_mode == "hatchery") {
  BLOCK_THRESH_ABOVE <- 0.03   # above bg_median to be "in block"
  DROP_THRESH <- 0.18           # blue-cross drop threshold
  MIN_BLOCK_SIM <- 0.52         # min surrounding sim for blue-cross context
  message("[01C] Mode: HATCHERY")
} else {
  BLOCK_THRESH_ABOVE <- 0.02
  DROP_THRESH <- 0.20
  MIN_BLOCK_SIM <- 0.55
  message("[01C] Mode: WILD")
}

DPI <- 300

BLOCK_PAL <- c(
  "#dc2626", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
  "#059669", "#e11d48", "#0284c7", "#65a30d", "#c026d3",
  "#0891b2", "#ca8a04", "#4f46e5", "#15803d", "#be185d",
  "#6366f1", "#ea580c", "#0d9488", "#a21caf", "#78716c",
  "#b91c1c", "#1d4ed8", "#047857", "#b45309", "#6d28d9",
  "#0e7490", "#a16207", "#4338ca", "#166534", "#9d174d"
)

# =============================================================================
# LOAD
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
precomp_list <- list()
n_rds_total   <- length(rds_files)
n_rds_loaded  <- 0L
n_rds_failed  <- 0L
load_failures <- character(0)
for (f in rds_files) {
  obj <- tryCatch(readRDS(f), error = function(e) {
    load_failures <<- c(load_failures, paste0(basename(f), ": ", conditionMessage(e)))
    NULL
  })
  if (is.null(obj)) { n_rds_failed <- n_rds_failed + 1L; next }
  precomp_list[[obj$chrom]] <- obj
  n_rds_loaded <- n_rds_loaded + 1L
}
chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { sample_names <- sub("^PC_1_", "", pc1_cols); break }
}
n_samples <- length(sample_names)

# ── Loading QC report ──
message("[01C] ── Load summary ──")
message(sprintf("  RDS files found:       %d in %s", n_rds_total, precomp_dir))
message(sprintf("  RDS loaded:            %d", n_rds_loaded))
if (n_rds_failed > 0) {
  message(sprintf("  RDS failed:            %d", n_rds_failed))
  for (lf in head(load_failures, 5)) message("    FAIL: ", lf)
  if (length(load_failures) > 5) message("    ... and ", length(load_failures) - 5, " more")
}
if (!is.null(chrom_filter)) {
  message(sprintf("  Chrom filter applied:  keeping %d of %d chroms matching '%s'",
    length(chroms), n_rds_loaded, paste(chrom_filter, collapse = ",")))
}
message(sprintf("  Samples detected:      %d (from PC_1_* columns)", n_samples))

# Per-chrom window count summary (before filtering chroms with <50 windows)
if (length(chroms) > 0) {
  n_win_per_chr <- sapply(chroms, function(c) precomp_list[[c]]$n_windows %||% nrow(precomp_list[[c]]$dt))
  total_windows <- sum(n_win_per_chr)
  small_chroms  <- chroms[n_win_per_chr < 50]
  message(sprintf("  Total windows:         %d across %d chroms (mean=%d/chr median=%d/chr min=%d max=%d)",
    total_windows, length(chroms),
    as.integer(round(mean(n_win_per_chr))), as.integer(stats::median(n_win_per_chr)),
    min(n_win_per_chr), max(n_win_per_chr)))
  if (length(small_chroms) > 0) {
    message(sprintf("  Chroms <50 windows (will be skipped): %d  [%s]",
      length(small_chroms), paste(small_chroms, collapse = ",")))
  }
  # Mean kb per window (window size QC — helps catch precomp config drift)
  pc_first <- precomp_list[[chroms[1]]]
  if (!is.null(pc_first$dt) && "start_bp" %in% names(pc_first$dt) && "end_bp" %in% names(pc_first$dt)) {
    widths <- pc_first$dt$end_bp - pc_first$dt$start_bp
    if (length(widths) > 0) {
      message(sprintf("  Window size (on %s): mean=%.1f kb median=%.1f kb range=[%d, %d] bp",
        chroms[1], mean(widths) / 1000, stats::median(widths) / 1000,
        min(widths), max(widths)))
    }
  }
}

# Optional: assembly gaps + AGP
gap_dt <- data.table()
if (!is.null(gaps_file) && file.exists(gaps_file)) {
  gap_dt <- fread(gaps_file, header = FALSE)
  if (ncol(gap_dt) >= 3) setnames(gap_dt, 1:3, c("chrom", "start", "end"))
}
agp_dt <- data.table()
if (!is.null(agp_file) && file.exists(agp_file)) {
  agp_raw <- fread(agp_file, header = FALSE, fill = TRUE, sep = "\t")
  if (ncol(agp_raw) >= 6) agp_dt <- agp_raw[V5 == "W", .(chrom = V1, start = V2, end = V3, contig = V6)]
}

# =============================================================================
# HELPERS
# =============================================================================

# Cluster samples into 3 bands on PC1 within a set of windows
get_bands <- function(dt, idxs, available_pc1) {
  if (length(idxs) < 3 || length(available_pc1) < 20) return(NULL)
  mat <- as.matrix(dt[idxs, ..available_pc1])
  avg <- colMeans(mat, na.rm = TRUE)
  valid <- is.finite(avg)
  if (sum(valid) < 20) return(NULL)
  vals <- avg[valid]; snames <- sub("^PC_1_", "", names(vals))
  km <- tryCatch(kmeans(vals, centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(km)) return(NULL)
  co <- order(km$centers[, 1])
  bands <- integer(length(vals))
  for (b in 1:3) bands[km$cluster == co[b]] <- b
  list(bands = setNames(bands, snames), pc1 = setNames(as.numeric(vals), snames))
}

# Concordance between two band assignments
band_concordance <- function(bands_a, bands_b) {
  shared <- intersect(names(bands_a), names(bands_b))
  if (length(shared) < 20) return(list(concordance = NA, jaccard = NA))
  ba <- bands_a[shared]; bb <- bands_b[shared]
  conc <- sum(ba == bb) / length(shared)
  jac <- mean(vapply(1:3, function(b) {
    sa <- which(ba == b); sb <- which(bb == b)
    inter <- length(intersect(sa, sb)); union <- length(union(sa, sb))
    if (union > 0) inter / union else 0
  }, numeric(1)))
  list(concordance = conc, jaccard = jac)
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_blocks <- list(); all_boundaries <- list(); all_crosses <- list()
all_concordances <- list(); all_pa <- list()

# Per-chromosome summary rows — written to disk at the end.
chr_summary_rows <- list()
# Track skip reasons for the QC block at the end.
chr_skip_log <- list()
# Sink for chroms that pass every gate but yield zero blocks (vs chroms skipped upfront).
n_chroms_input <- length(chroms)
n_chroms_too_small <- 0L
n_chroms_no_pc1 <- 0L
n_chroms_no_blocks <- 0L
n_chroms_processed <- 0L

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 50) {
    n_chroms_too_small <- n_chroms_too_small + 1L
    chr_skip_log[[chr]] <- sprintf("skipped: n_windows=%s < 50",
      if (is.null(pc)) "NULL" else as.character(pc$n_windows))
    next
  }
  dt <- pc$dt; sim_mat <- pc$sim_mat; n <- nrow(dt)
  pc1_cols <- paste0("PC_1_", sample_names)
  available_pc1 <- intersect(pc1_cols, names(dt))
  if (length(available_pc1) < 20) {
    n_chroms_no_pc1 <- n_chroms_no_pc1 + 1L
    chr_skip_log[[chr]] <- sprintf("skipped: only %d PC_1_* columns (need >=20)", length(available_pc1))
    next
  }
  n_chroms_processed <- n_chroms_processed + 1L

  message("\n[01C] === ", chr, " === (", n, " windows)")

  # =========================================================================
  # STEP 1: Local similarity profile + block detection
  # =========================================================================
  local_sim <- numeric(n)
  for (wi in seq_len(n)) {
    lo <- max(1, wi - LOCAL_RANGE); hi <- min(n, wi + LOCAL_RANGE)
    nb <- lo:hi; nb <- nb[nb != wi]
    local_sim[wi] <- mean(sim_mat[wi, nb], na.rm = TRUE)
  }

  bg_median <- median(local_sim, na.rm = TRUE)
  in_block <- local_sim > bg_median + BLOCK_THRESH_ABOVE

  # Connected components
  block_id <- integer(n); bid <- 0L
  for (wi in seq_len(n)) {
    if (!in_block[wi]) { block_id[wi] <- 0L; next }
    if (wi == 1 || !in_block[wi - 1]) bid <- bid + 1L
    block_id[wi] <- bid
  }

  # Filter small + take top N
  block_tab <- table(block_id[block_id > 0])
  valid_blocks <- as.integer(names(block_tab[block_tab >= MIN_BLOCK]))
  block_id[!block_id %in% valid_blocks] <- 0L
  if (length(valid_blocks) == 0) {
    n_chroms_no_blocks <- n_chroms_no_blocks + 1L
    chr_skip_log[[chr]] <- sprintf("no blocks found (local_sim bg_median=%.4f, threshold>+%.2f)",
      bg_median, BLOCK_THRESH_ABOVE)
    # Still emit a summary row so the user sees "I processed this chrom but found nothing"
    # Columns mirror the full summary row below (same order, same names) so the
    # per-chrom TSV is a proper rectangular table.
    chr_summary_rows[[length(chr_summary_rows) + 1L]] <- data.table(
      chrom = chr,
      n_windows = n,
      n_samples_with_pc1 = length(available_pc1),
      bg_median_sim = round(bg_median, 4),
      frac_in_block = 0,
      n_blocks = 0L,
      mean_block_windows    = NA_real_,
      median_block_windows  = NA_real_,
      max_block_windows     = NA_integer_,
      total_block_windows   = 0L,
      mean_block_mean_sim   = NA_real_,
      # Boundary-type counts (outer + inner total)
      n_boundaries_clear_hard      = 0L,
      n_boundaries_clear_soft      = 0L,
      n_boundaries_diffuse         = 0L,
      n_boundaries_chromosome_edge = 0L,
      n_boundaries_inner           = 0L,
      # Blue-cross inner_type breakdown
      n_blue_crosses                  = 0L,
      n_inner_hard_same_system        = 0L,
      n_inner_hard_different_systems  = 0L,
      n_inner_hard_assembly           = 0L,
      n_inner_hard_suspect_assembly   = 0L,
      n_inner_soft_boundary_candidate = 0L,
      n_inner_ambiguous               = 0L,
      # Concordance stats (off-diagonal)
      mean_concordance    = NA_real_,
      median_concordance  = NA_real_,
      n_high_concordance_pairs = 0L,
      n_low_concordance_pairs  = 0L,
      n_blocks_with_3bands     = 0L
    )
    message(sprintf("[01C]   %s: no blocks (bg_median=%.4f)", chr, bg_median))
    next
  }

  block_sizes <- sort(block_tab[as.character(valid_blocks)], decreasing = TRUE)
  keep <- as.integer(names(block_sizes))[seq_len(min(TOP_BLOCKS, length(block_sizes)))]
  new_id <- setNames(seq_along(keep), as.character(keep))
  block_label <- character(n)
  for (wi in seq_len(n)) {
    bid_str <- as.character(block_id[wi])
    block_label[wi] <- if (bid_str %in% names(new_id)) paste0("B", new_id[bid_str]) else "bg"
  }
  n_blocks <- length(keep)
  message("  ", n_blocks, " blocks detected")

  # Block metadata
  block_meta <- list()
  for (bi in seq_len(n_blocks)) {
    bname <- paste0("B", bi)
    idxs <- which(block_label == bname)
    block_meta[[bi]] <- list(
      name = bname, idxs = idxs,
      start_bp = min(dt$start_bp[idxs]), end_bp = max(dt$end_bp[idxs]),
      start_mb = round(min(dt$start_bp[idxs]) / 1e6, 3),
      end_mb = round(max(dt$end_bp[idxs]) / 1e6, 3),
      n_windows = length(idxs),
      mean_sim = round(mean(local_sim[idxs]), 4)
    )
  }

  # =========================================================================
  # STEP 2: Sample grouping per block
  # =========================================================================
  block_bands <- list(); block_pc1 <- list()
  for (bi in seq_len(n_blocks)) {
    bm <- block_meta[[bi]]
    bg <- get_bands(dt, bm$idxs, available_pc1)
    if (!is.null(bg)) {
      block_bands[[bm$name]] <- bg$bands
      block_pc1[[bm$name]] <- bg$pc1
    }
  }

  # =========================================================================
  # STEP 3: Boundary classification
  # =========================================================================
  # For each pair of adjacent blocks AND for each block edge, classify
  # the boundary transition shape from the local_sim profile.

  boundary_rows <- list()
  for (bi in seq_len(n_blocks)) {
    bm <- block_meta[[bi]]
    idxs <- bm$idxs

    # LEFT boundary of this block
    left_edge <- min(idxs)
    if (left_edge > 5) {
      # Measure transition: how many windows does it take for sim to drop?
      profile <- local_sim[max(1, left_edge - 20):left_edge]
      inside_mean <- mean(local_sim[idxs[1:min(10, length(idxs))]])
      outside_mean <- mean(local_sim[max(1, left_edge - 20):(left_edge - 1)])
      drop_mag <- inside_mean - outside_mean

      # Count windows for transition (sim > midpoint)
      midpoint <- (inside_mean + outside_mean) / 2
      trans_width <- sum(profile > midpoint & profile < inside_mean - 0.02)

      btype <- if (drop_mag > 0.15 && trans_width <= 3) "clear_hard"
               else if (drop_mag > 0.10 && trans_width <= 8) "clear_soft"
               else if (drop_mag > 0.05) "diffuse"
               else "chromosome_edge"

      boundary_rows[[length(boundary_rows) + 1]] <- data.table(
        chrom = chr, block = bm$name, side = "left",
        global_window_id = if ("global_window_id" %in% names(dt)) dt$global_window_id[left_edge] else NA_integer_,
        window_idx = left_edge,
        pos_bp = dt$start_bp[left_edge],
        pos_mb = round(dt$start_bp[left_edge] / 1e6, 3),
        drop_magnitude = round(drop_mag, 4),
        transition_width = trans_width,
        boundary_type = btype
      )
    }

    # RIGHT boundary
    right_edge <- max(idxs)
    if (right_edge < n - 5) {
      profile <- local_sim[right_edge:min(n, right_edge + 20)]
      inside_mean <- mean(local_sim[idxs[max(1, length(idxs) - 9):length(idxs)]])
      outside_mean <- mean(local_sim[(right_edge + 1):min(n, right_edge + 20)])
      drop_mag <- inside_mean - outside_mean

      midpoint <- (inside_mean + outside_mean) / 2
      trans_width <- sum(profile > midpoint & profile < inside_mean - 0.02)

      btype <- if (drop_mag > 0.15 && trans_width <= 3) "clear_hard"
               else if (drop_mag > 0.10 && trans_width <= 8) "clear_soft"
               else if (drop_mag > 0.05) "diffuse"
               else "chromosome_edge"

      boundary_rows[[length(boundary_rows) + 1]] <- data.table(
        chrom = chr, block = bm$name, side = "right",
        global_window_id = if ("global_window_id" %in% names(dt)) dt$global_window_id[right_edge] else NA_integer_,
        window_idx = right_edge,
        pos_bp = dt$end_bp[right_edge],
        pos_mb = round(dt$end_bp[right_edge] / 1e6, 3),
        drop_magnitude = round(drop_mag, 4),
        transition_width = trans_width,
        boundary_type = btype
      )
    }
  }

  # =========================================================================
  # STEP 4: Blue-cross detection inside blocks
  # =========================================================================
  cross_rows <- list()
  roll_k <- 20L

  for (wi in (roll_k + 1):(n - roll_k)) {
    if (block_label[wi] == "bg") next  # only inside blocks
    surr_left  <- mean(local_sim[max(1, wi - roll_k):(wi - 3)], na.rm = TRUE)
    surr_right <- mean(local_sim[(wi + 3):min(n, wi + roll_k)], na.rm = TRUE)
    surr_mean  <- (surr_left + surr_right) / 2
    if (surr_mean < MIN_BLOCK_SIM) next

    drop <- surr_mean - local_sim[wi]
    if (drop < DROP_THRESH) next

    # Check zone
    zone <- max(1, wi - 2):min(n, wi + 2)
    zone_drop <- surr_mean - mean(local_sim[zone], na.rm = TRUE)
    if (zone_drop < DROP_THRESH * 0.7) next

    cross_rows[[length(cross_rows) + 1]] <- data.table(
      chrom = chr,
      global_window_id = if ("global_window_id" %in% names(dt)) dt$global_window_id[wi] else NA_integer_,
      window_idx = wi, block = block_label[wi],
      pos_mb = round((dt$start_bp[wi] + dt$end_bp[wi]) / 2e6, 4),
      start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
      drop = round(drop, 4), zone_drop = round(zone_drop, 4)
    )
  }

  if (length(cross_rows) > 0) {
    cross_dt <- rbindlist(cross_rows)
    # Merge nearby crosses. The summary picks the strongest (max drop) and
    # keeps its global_window_id as the identifier; chrom is constant across
    # all rows of this chromosome so pass it through directly.
    cross_dt[, group := cumsum(c(TRUE, diff(window_idx) > 10))]
    cross_dt <- cross_dt[, .(
      chrom = chrom[1],
      global_window_id = global_window_id[which.max(drop)],
      window_idx = window_idx[which.max(drop)],
      block = block[1],
      pos_mb = pos_mb[which.max(drop)],
      start_bp = min(start_bp), end_bp = max(end_bp),
      drop = max(drop), zone_drop = max(zone_drop),
      n_cross_windows = .N
    ), by = group][, group := NULL]

    # Test each cross: concordance across it
    for (ci in seq_len(nrow(cross_dt))) {
      cx <- cross_dt[ci]; wi <- cx$window_idx
      left_idx <- max(1, wi - FLANK_W):(wi - 2)
      right_idx <- (wi + 2):min(n, wi + FLANK_W)
      if (length(left_idx) < 5 || length(right_idx) < 5) next

      bg_left <- get_bands(dt, left_idx, available_pc1)
      bg_right <- get_bands(dt, right_idx, available_pc1)
      if (is.null(bg_left) || is.null(bg_right)) next

      cc <- band_concordance(bg_left$bands, bg_right$bands)
      # BUGFIX 2026-04-17: cor() throws a warning and returns NA when the
      # input length is 0 or 1, and NaN when all values are identical.
      # Compute the intersect once, guard on length >= 3.
      shared_s <- intersect(names(bg_left$pc1), names(bg_right$pc1))
      pc1_cor <- if (length(shared_s) >= 3) {
        tryCatch(cor(bg_left$pc1[shared_s], bg_right$pc1[shared_s]),
                 warning = function(w) NA_real_,
                 error = function(e) NA_real_)
      } else NA_real_

      # Gap / contig boundary check
      in_gap <- FALSE; at_ctg <- FALSE
      mid_bp <- (cx$start_bp + cx$end_bp) / 2
      if (nrow(gap_dt) > 0) in_gap <- any(gap_dt$chrom == chr & gap_dt$start <= cx$end_bp & gap_dt$end >= cx$start_bp)
      if (nrow(agp_dt) > 0) {
        chr_agp <- agp_dt[chrom == chr]
        if (nrow(chr_agp) > 0) at_ctg <- min(abs(c(chr_agp$start, chr_agp$end) - mid_bp)) < 50000
      }

      # Classify the cross as inner boundary type
      conc <- cc$concordance; jac <- cc$jaccard
      inner_type <- if (in_gap || at_ctg) {
        if (conc > 0.7) "inner_hard_assembly" else "inner_hard_suspect_assembly"
      } else if (conc > 0.85 && jac > 0.7) "inner_hard_same_system"
      else if (conc < 0.50) "inner_hard_different_systems"
      else if (conc < 0.65) "inner_soft_boundary_candidate"
      else "inner_ambiguous"

      cross_dt[ci, `:=`(
        concordance = round(conc, 3), jaccard = round(jac, 3),
        pc1_cor = round(pc1_cor, 3),
        in_gap = in_gap, at_contig_boundary = at_ctg,
        inner_type = inner_type
      )]

      # Also add to boundary catalog
      boundary_rows[[length(boundary_rows) + 1]] <- data.table(
        chrom = chr, block = cx$block, side = "inner",
        global_window_id = cx$global_window_id,
        window_idx = cx$window_idx,
        pos_bp = as.integer(mid_bp),
        pos_mb = round(mid_bp / 1e6, 3),
        drop_magnitude = cx$drop,
        transition_width = cx$n_cross_windows,
        boundary_type = inner_type
      )
    }
    cross_dt[, chrom := chr]
    all_crosses[[length(all_crosses) + 1]] <- cross_dt
    message("  Blue crosses: ", nrow(cross_dt))
  }

  # =========================================================================
  # STEP 5: Block concordance matrix
  # =========================================================================
  bnames <- names(block_bands)
  n_b <- length(bnames)
  if (n_b >= 2) {
    conc_rows <- list()
    for (bi in seq_len(n_b)) {
      for (bj in bi:n_b) {
        cc <- band_concordance(block_bands[[bi]], block_bands[[bj]])
        if (is.na(cc$concordance)) next
        conc_rows[[length(conc_rows) + 1]] <- data.table(
          chrom = chr, block_a = bnames[bi], block_b = bnames[bj],
          concordance = round(cc$concordance, 3),
          jaccard = round(cc$jaccard, 3)
        )
        if (bi != bj) {
          conc_rows[[length(conc_rows) + 1]] <- data.table(
            chrom = chr, block_a = bnames[bj], block_b = bnames[bi],
            concordance = round(cc$concordance, 3),
            jaccard = round(cc$jaccard, 3)
          )
        }
      }
    }
    all_concordances[[length(all_concordances) + 1]] <- rbindlist(conc_rows)
    message("  Concordance pairs: ", length(conc_rows))
  }

  # =========================================================================
  # STEP 6: Block registry + window PA
  # =========================================================================
  block_reg_rows <- list()
  for (bi in seq_len(n_blocks)) {
    bm <- block_meta[[bi]]
    has_bands <- bm$name %in% names(block_bands)
    # First/last window of the block (in local order). These paired global_window_ids
    # let downstream scripts join back to the precompute dt without needing to match
    # on bp coordinates.
    first_idx <- min(bm$idxs); last_idx <- max(bm$idxs)
    block_reg_rows[[bi]] <- data.table(
      chrom = chr, block_id = bm$name,
      start_bp = bm$start_bp, end_bp = bm$end_bp,
      start_mb = bm$start_mb, end_mb = bm$end_mb,
      start_global_window_id = if ("global_window_id" %in% names(dt)) dt$global_window_id[first_idx] else NA_integer_,
      end_global_window_id   = if ("global_window_id" %in% names(dt)) dt$global_window_id[last_idx]  else NA_integer_,
      n_windows = bm$n_windows, mean_sim = bm$mean_sim,
      has_3bands = has_bands
    )
  }
  all_blocks[[length(all_blocks) + 1]] <- rbindlist(block_reg_rows)

  # Window PA
  pa_rows <- list()
  for (wi in seq_len(n)) {
    pa_rows[[wi]] <- data.table(
      chrom = chr, global_window_id = dt$global_window_id[wi],
      pos_mb = round((dt$start_bp[wi] + dt$end_bp[wi]) / 2e6, 4),
      block = block_label[wi],
      local_sim = round(local_sim[wi], 4),
      in_block = block_label[wi] != "bg"
    )
  }
  all_pa[[length(all_pa) + 1]] <- rbindlist(pa_rows)
  all_boundaries[[length(all_boundaries) + 1]] <- rbindlist(boundary_rows, fill = TRUE)

  # =========================================================================
  # PER-CHROMOSOME SUMMARY ROW
  # =========================================================================
  # Assemble everything useful about what happened on this chromosome for
  # the summary TSV and the end-of-run genome-wide report.

  # Block-size stats
  block_win_counts <- sapply(block_meta, function(bm) bm$n_windows)
  block_mean_sims  <- sapply(block_meta, function(bm) bm$mean_sim)
  total_block_windows <- sum(block_win_counts)
  frac_in_block <- round(total_block_windows / n, 4)
  mean_block_win <- round(mean(block_win_counts), 1)
  median_block_win <- stats::median(block_win_counts)
  max_block_win <- as.integer(max(block_win_counts))
  mean_block_mean_sim <- round(mean(block_mean_sims), 4)
  n_blocks_with_3bands <- sum(sapply(block_meta, function(bm) bm$name %in% names(block_bands)))

  # Boundary-type counts (from this chromosome's boundary_rows)
  bt_counts <- c(clear_hard = 0L, clear_soft = 0L, diffuse = 0L,
                 chromosome_edge = 0L, inner = 0L)
  if (length(boundary_rows) > 0) {
    chr_bounds <- rbindlist(boundary_rows, fill = TRUE)
    if ("boundary_type" %in% names(chr_bounds)) {
      # Outer boundary types
      for (bt in c("clear_hard", "clear_soft", "diffuse", "chromosome_edge")) {
        bt_counts[bt] <- sum(chr_bounds$boundary_type == bt, na.rm = TRUE)
      }
      # Inner boundaries: everything with side == "inner" (inner_hard_*, inner_soft_*, inner_ambiguous)
      if ("side" %in% names(chr_bounds)) {
        bt_counts["inner"] <- sum(chr_bounds$side == "inner", na.rm = TRUE)
      }
    }
  }

  # Blue-cross inner_type counts (only exists when crosses were found)
  it_counts <- c(inner_hard_same_system = 0L, inner_hard_different_systems = 0L,
                 inner_hard_assembly = 0L, inner_hard_suspect_assembly = 0L,
                 inner_soft_boundary_candidate = 0L, inner_ambiguous = 0L)
  n_blue <- 0L
  if (exists("cross_dt") && !is.null(cross_dt) && nrow(cross_dt) > 0) {
    n_blue <- nrow(cross_dt)
    if ("inner_type" %in% names(cross_dt)) {
      for (it in names(it_counts)) {
        it_counts[it] <- sum(cross_dt$inner_type == it, na.rm = TRUE)
      }
    }
  }

  # Concordance stats (exclude self-pairs where block_a == block_b)
  mean_conc <- NA_real_; median_conc <- NA_real_
  n_high_conc <- 0L; n_low_conc <- 0L
  if (exists("conc_rows") && length(conc_rows) > 0) {
    cdt <- rbindlist(conc_rows)
    off_diag <- cdt[block_a != block_b]
    if (nrow(off_diag) > 0) {
      mean_conc   <- round(mean(off_diag$concordance, na.rm = TRUE), 4)
      median_conc <- round(stats::median(off_diag$concordance, na.rm = TRUE), 4)
      n_high_conc <- sum(off_diag$concordance >= 0.85, na.rm = TRUE)
      n_low_conc  <- sum(off_diag$concordance <  0.50, na.rm = TRUE)
    }
  }

  chr_summary_rows[[length(chr_summary_rows) + 1L]] <- data.table(
    chrom = chr,
    n_windows = n,
    n_samples_with_pc1 = length(available_pc1),
    bg_median_sim = round(bg_median, 4),
    frac_in_block = frac_in_block,
    n_blocks = n_blocks,
    mean_block_windows    = mean_block_win,
    median_block_windows  = median_block_win,
    max_block_windows     = max_block_win,
    total_block_windows   = total_block_windows,
    mean_block_mean_sim   = mean_block_mean_sim,
    # Boundary-type counts (outer + inner total)
    n_boundaries_clear_hard      = bt_counts["clear_hard"],
    n_boundaries_clear_soft      = bt_counts["clear_soft"],
    n_boundaries_diffuse         = bt_counts["diffuse"],
    n_boundaries_chromosome_edge = bt_counts["chromosome_edge"],
    n_boundaries_inner           = bt_counts["inner"],
    # Blue-cross inner_type breakdown
    n_blue_crosses                  = n_blue,
    n_inner_hard_same_system        = it_counts["inner_hard_same_system"],
    n_inner_hard_different_systems  = it_counts["inner_hard_different_systems"],
    n_inner_hard_assembly           = it_counts["inner_hard_assembly"],
    n_inner_hard_suspect_assembly   = it_counts["inner_hard_suspect_assembly"],
    n_inner_soft_boundary_candidate = it_counts["inner_soft_boundary_candidate"],
    n_inner_ambiguous               = it_counts["inner_ambiguous"],
    # Concordance stats (off-diagonal, block_a != block_b)
    mean_concordance    = mean_conc,
    median_concordance  = median_conc,
    n_high_concordance_pairs = n_high_conc,
    n_low_concordance_pairs  = n_low_conc,
    n_blocks_with_3bands     = n_blocks_with_3bands
  )

  # Per-chrom log line
  message(sprintf("[01C]   %s: n_win=%d  bg_med=%.3f  blocks=%d (in_block=%.1f%%  mean_sz=%.0f)  boundaries=%d (hard=%d soft=%d diffuse=%d inner=%d)  blue_x=%d  mean_conc=%s",
    chr, n, bg_median, n_blocks,
    100 * frac_in_block, mean_block_win,
    length(boundary_rows),
    bt_counts["clear_hard"], bt_counts["clear_soft"], bt_counts["diffuse"], bt_counts["inner"],
    n_blue,
    if (is.finite(mean_conc)) sprintf("%.3f", mean_conc) else "NA"
  ))

  # Clean up per-chrom objects so the 'exists()' checks work correctly next iter
  if (exists("cross_dt")) rm(cross_dt)
  if (exists("conc_rows")) rm(conc_rows)
}

# =============================================================================
# WRITE
# =============================================================================

message("\n[01C] Writing...")

block_dt <- if (length(all_blocks) > 0) rbindlist(all_blocks, fill = TRUE) else data.table()
bound_dt <- if (length(all_boundaries) > 0) rbindlist(all_boundaries, fill = TRUE) else data.table()
cross_dt_all <- if (length(all_crosses) > 0) rbindlist(all_crosses, fill = TRUE) else data.table()
conc_dt <- if (length(all_concordances) > 0) rbindlist(all_concordances, fill = TRUE) else data.table()
pa_dt <- if (length(all_pa) > 0) rbindlist(all_pa, fill = TRUE) else data.table()

for (chr in unique(block_dt$chrom)) {
  fwrite(block_dt[chrom == chr], file.path(outdir, paste0("block_registry_", chr, ".tsv.gz")), sep = "\t")
  fwrite(bound_dt[chrom == chr], file.path(outdir, paste0("boundary_catalog_", chr, ".tsv.gz")), sep = "\t")
}
if (nrow(cross_dt_all) > 0) {
  for (chr in unique(cross_dt_all$chrom))
    fwrite(cross_dt_all[chrom == chr], file.path(outdir, paste0("blue_cross_verdicts_", chr, ".tsv.gz")), sep = "\t")
}
if (nrow(conc_dt) > 0) {
  for (chr in unique(conc_dt$chrom))
    fwrite(conc_dt[chrom == chr], file.path(outdir, paste0("block_concordance_", chr, ".tsv.gz")), sep = "\t")
}
fwrite(pa_dt, file.path(outdir, "01C_window_pa.tsv.gz"), sep = "\t")

# =============================================================================
# PER-CHROM SUMMARY TSV + GENOME-WIDE REPORT
# =============================================================================

chr_summary_dt <- if (length(chr_summary_rows) > 0) rbindlist(chr_summary_rows) else data.table()
if (nrow(chr_summary_dt) > 0) {
  f_chr_summ <- file.path(outdir, "01C_chrom_summary.tsv")
  fwrite(chr_summary_dt, f_chr_summ, sep = "\t")
  message("[01C] Per-chrom summary: ", f_chr_summ)
}

# ── Processing-stage accounting ──
message("\n[01C] ═══ PROCESSING SUMMARY ═══")
message(sprintf("  Chromosomes input:              %d", n_chroms_input))
message(sprintf("  Processed (blocks found):       %d", nrow(chr_summary_dt) - n_chroms_no_blocks))
message(sprintf("  Processed (no blocks found):    %d", n_chroms_no_blocks))
message(sprintf("  Skipped (<50 windows):          %d", n_chroms_too_small))
message(sprintf("  Skipped (<20 PC_1_* cols):      %d", n_chroms_no_pc1))
if (length(chr_skip_log) > 0) {
  message("")
  message("  Skip reasons:")
  for (chr in names(chr_skip_log)) {
    message(sprintf("    %-20s %s", chr, chr_skip_log[[chr]]))
  }
}

# ── Genome-wide block stats ──
if (nrow(chr_summary_dt) > 0 && any(chr_summary_dt$n_blocks > 0)) {
  chr_with_blocks <- chr_summary_dt[n_blocks > 0]
  message("\n[01C] ═══ BLOCK SUMMARY (genome-wide) ═══")
  message(sprintf("  Total blocks detected:          %d", sum(chr_with_blocks$n_blocks)))
  message(sprintf("  Chroms with blocks:             %d / %d",
    nrow(chr_with_blocks), nrow(chr_summary_dt)))
  message(sprintf("  Blocks per chrom:               mean=%.1f  median=%.0f  min=%d  max=%d",
    mean(chr_with_blocks$n_blocks), stats::median(chr_with_blocks$n_blocks),
    min(chr_with_blocks$n_blocks), max(chr_with_blocks$n_blocks)))
  message(sprintf("  Windows in blocks (frac):       mean=%.2f  median=%.2f  max=%.2f",
    mean(chr_with_blocks$frac_in_block),
    stats::median(chr_with_blocks$frac_in_block),
    max(chr_with_blocks$frac_in_block)))
  message(sprintf("  Mean block size (windows):      mean=%.1f  median=%.0f  max=%d",
    mean(chr_with_blocks$mean_block_windows, na.rm = TRUE),
    stats::median(chr_with_blocks$mean_block_windows, na.rm = TRUE),
    max(chr_with_blocks$max_block_windows, na.rm = TRUE)))
  message(sprintf("  Blocks with valid 3-band split: %d / %d (%.0f%%)",
    sum(chr_with_blocks$n_blocks_with_3bands),
    sum(chr_with_blocks$n_blocks),
    100 * sum(chr_with_blocks$n_blocks_with_3bands) / max(1L, sum(chr_with_blocks$n_blocks))))
}

# ── Boundary-type totals ──
if (nrow(bound_dt) > 0) {
  message("\n[01C] ═══ BOUNDARY SUMMARY (genome-wide) ═══")
  for (bt in sort(unique(bound_dt$boundary_type))) {
    n_bt <- sum(bound_dt$boundary_type == bt, na.rm = TRUE)
    pct <- 100 * n_bt / nrow(bound_dt)
    message(sprintf("  %-32s %6d  (%5.1f%%)", bt, n_bt, pct))
  }
  # Outer vs inner breakdown
  if ("side" %in% names(bound_dt)) {
    n_outer <- sum(bound_dt$side %in% c("left", "right"), na.rm = TRUE)
    n_inner <- sum(bound_dt$side == "inner", na.rm = TRUE)
    message(sprintf("  [outer edges: %d | inner: %d]", n_outer, n_inner))
  }
}

# ── Blue-cross totals ──
if (nrow(cross_dt_all) > 0) {
  message("\n[01C] ═══ BLUE CROSS SUMMARY (genome-wide) ═══")
  message(sprintf("  Total blue crosses detected: %d", nrow(cross_dt_all)))
  for (it in sort(unique(cross_dt_all$inner_type[!is.na(cross_dt_all$inner_type)]))) {
    n_it <- sum(cross_dt_all$inner_type == it, na.rm = TRUE)
    pct <- 100 * n_it / nrow(cross_dt_all)
    message(sprintf("  %-32s %6d  (%5.1f%%)", it, n_it, pct))
  }
  n_assembly <- sum(cross_dt_all$inner_type %in%
    c("inner_hard_assembly", "inner_hard_suspect_assembly"), na.rm = TRUE)
  n_different_systems <- sum(cross_dt_all$inner_type == "inner_hard_different_systems", na.rm = TRUE)
  n_soft_cand <- sum(cross_dt_all$inner_type == "inner_soft_boundary_candidate", na.rm = TRUE)
  message("")
  message(sprintf("  Assembly-error candidates:      %d", n_assembly))
  message(sprintf("  Different-systems (real bdy):   %d", n_different_systems))
  message(sprintf("  Soft-boundary candidates:       %d", n_soft_cand))
}

# ── Concordance totals ──
if (nrow(conc_dt) > 0) {
  off_diag <- conc_dt[block_a != block_b]
  if (nrow(off_diag) > 0) {
    message("\n[01C] ═══ CONCORDANCE SUMMARY (genome-wide, off-diagonal pairs) ═══")
    message(sprintf("  Pairs examined:         %d", nrow(off_diag)))
    message(sprintf("  Concordance:            mean=%.3f  median=%.3f  min=%.3f  max=%.3f",
      mean(off_diag$concordance, na.rm = TRUE),
      stats::median(off_diag$concordance, na.rm = TRUE),
      min(off_diag$concordance, na.rm = TRUE),
      max(off_diag$concordance, na.rm = TRUE)))
    message(sprintf("  Jaccard:                mean=%.3f  median=%.3f",
      mean(off_diag$jaccard, na.rm = TRUE),
      stats::median(off_diag$jaccard, na.rm = TRUE)))
    message(sprintf("  High concordance >=0.85 (likely same system): %d (%.1f%%)",
      sum(off_diag$concordance >= 0.85, na.rm = TRUE),
      100 * sum(off_diag$concordance >= 0.85, na.rm = TRUE) / nrow(off_diag)))
    message(sprintf("  Low  concordance <0.50  (different systems): %d (%.1f%%)",
      sum(off_diag$concordance <  0.50, na.rm = TRUE),
      100 * sum(off_diag$concordance <  0.50, na.rm = TRUE) / nrow(off_diag)))
  }
}

message("\n[DONE] -> ", outdir)
