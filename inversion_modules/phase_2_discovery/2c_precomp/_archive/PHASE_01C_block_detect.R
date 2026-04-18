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
for (f in rds_files) { obj <- readRDS(f); precomp_list[[obj$chrom]] <- obj }
chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { sample_names <- sub("^PC_1_", "", pc1_cols); break }
}
n_samples <- length(sample_names)
message("[01C] ", length(chroms), " chromosomes, ", n_samples, " samples")

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

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 50) next
  dt <- pc$dt; sim_mat <- pc$sim_mat; n <- nrow(dt)
  pc1_cols <- paste0("PC_1_", sample_names)
  available_pc1 <- intersect(pc1_cols, names(dt))
  if (length(available_pc1) < 20) next

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
  if (length(valid_blocks) == 0) { message("  No blocks"); next }

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
      window_idx = wi, block = block_label[wi],
      pos_mb = round((dt$start_bp[wi] + dt$end_bp[wi]) / 2e6, 4),
      start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
      drop = round(drop, 4), zone_drop = round(zone_drop, 4)
    )
  }

  if (length(cross_rows) > 0) {
    cross_dt <- rbindlist(cross_rows)
    # Merge nearby
    cross_dt[, group := cumsum(c(TRUE, diff(window_idx) > 10))]
    cross_dt <- cross_dt[, .(
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
    block_reg_rows[[bi]] <- data.table(
      chrom = chr, block_id = bm$name,
      start_bp = bm$start_bp, end_bp = bm$end_bp,
      start_mb = bm$start_mb, end_mb = bm$end_mb,
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

  message("  Blocks: ", n_blocks, " | Boundaries: ", length(boundary_rows))
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

# Summary
if (nrow(bound_dt) > 0) {
  message("\n[01C] === BOUNDARY SUMMARY ===")
  for (bt in sort(unique(bound_dt$boundary_type)))
    message("  ", bt, ": ", sum(bound_dt$boundary_type == bt))
}
if (nrow(cross_dt_all) > 0) {
  message("\n[01C] === BLUE CROSS SUMMARY ===")
  for (it in sort(unique(cross_dt_all$inner_type[!is.na(cross_dt_all$inner_type)])))
    message("  ", it, ": ", sum(cross_dt_all$inner_type == it, na.rm = TRUE))
}

message("\n[DONE] -> ", outdir)
