#!/usr/bin/env Rscript

# =============================================================================
# 01_dosage_signal.R  (v1.1 — registry-wired)
#
# DOSAGE-BASED INVERSION SIGNAL — first stage of the breakpoint refinement
# workflow.
#
# For each candidate:
#   1. Load dosage matrix + rotated PCA group assignments (from registry)
#   2. Identify informative markers:
#        delta_i = mean(X_i | HOMO_2) - mean(X_i | HOMO_1)
#        informative if |delta_i| >= DELTA_MIN
#   3. Build the marker × marker correlation graph on informative markers
#   4. Extract the largest connected component with correlation >= RHO_BLOCK
#      as the core block
#   5. Compute per-sample core consensus:
#        c_j = mean(X_ij | i in top-tier core)
#   6. Extend block boundaries outward one SNP at a time:
#        for each candidate extension SNP i,
#          r_i = |cor(X_i, c)|
#        stop when r_i < RHO_EXT
#
# REGISTRY WIRING (chat-18):
#   - Candidate coords come from reg$intervals$get_candidate(cid)
#   - Karyotype assignments come from reg$samples$get_groups_for_candidate(cid)
#   - Dosage stays on disk (path via $DOSAGE_DIR env or config)
#   - Output block "dosage_blocks" is written via reg$evidence$write_block
#   - Per-marker table stays as raw TSV under evidence_registry/per_candidate/<cid>/raw/
#
# OUTPUTS:
#   evidence_registry/per_candidate/<cid>/raw/dosage_informative_markers.tsv.gz
#   evidence_registry/per_candidate/<cid>/structured/dosage_blocks.json
#
# References:
#   find_informative_markers     -- ported from STEP_C01h
#   find_coseg_blocks            -- ported from STEP_C01i
#   extend_block_boundaries      -- ported from STEP_C01i
#
# Usage:
#   Rscript 01_dosage_signal.R [cid=all]
#   Rscript 01_dosage_signal.R LG28_1          # single candidate
#   Rscript 01_dosage_signal.R --config my_overrides.R LG28_1
# =============================================================================

# --- Entry: one line sources registry + ancestry stack -----------------------
Sys.setenv(CURRENT_SCRIPT = "01_dosage_signal.R")
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
if (!file.exists(.bridge)) {
  # Fallback search: repo-root /utils/, $BASE/utils/
  for (p in c("utils/registry_bridge.R",
              "../utils/registry_bridge.R",
              file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
    if (file.exists(p)) { .bridge <- p; break }
  }
}
if (!file.exists(.bridge)) {
  stop("[01_dosage_signal] cannot locate utils/registry_bridge.R. ",
       "Set $REGISTRY_BRIDGE or $BASE and retry.")
}
source(.bridge)
# After sourcing: reg, smap, get_Q, get_region_stats, BRIDGE_PATHS available

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- CLI: optional --config + candidate id filter ---------------------------
# New call pattern:
#   Rscript 01_dosage_signal.R [--config <overrides.R>] [<cid>|all]
args <- commandArgs(trailingOnly = TRUE)
config_file <- ""
cid_filter  <- NA_character_
.ai <- 1L
while (.ai <= length(args)) {
  a <- args[.ai]
  if (a == "--config" && .ai < length(args)) {
    config_file <- args[.ai + 1L]
    .ai <- .ai + 2L
  } else if (a != "all") {
    cid_filter <- a
    .ai <- .ai + 1L
  } else {
    .ai <- .ai + 1L
  }
}
if (nzchar(config_file) && file.exists(config_file)) {
  source(config_file)   # config file is now for PARAMETER overrides only, not paths
}

# --- Paths: dosage still on disk; everything else from registry -------------
# DOSAGE_DIR is the only required external path (dosage files are too large
# to live in the registry). The registry handles candidate coords, sample
# groups, output writes.
if (!exists("DOSAGE_DIR")) {
  DOSAGE_DIR <- Sys.getenv("DOSAGE_DIR", "")
  if (!nzchar(DOSAGE_DIR)) {
    DOSAGE_DIR <- file.path(BRIDGE_PATHS$BASE, "popstruct_thin",
                             "04_beagle_byRF_majmin")
  }
}

if (!exists("ensure_dir"))
  ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }

# =============================================================================
# Parameters
# =============================================================================
BP01_PARAMS <- list(
  delta_min            = 1.0,    # |h2_mean - h1_mean| >= this to be informative
  rho_block            = 0.7,    # min correlation for core-block adjacency
  rho_ext              = 0.6,    # min correlation to extend past block edge
  min_core_size        = 20L,    # min SNPs in the core block
  max_extend_snps      = 5000L,  # max SNPs to scan past block edge
  top_core_quantile    = 0.75,   # core markers for consensus: top this quantile
  subsample_cor_n      = 5000L,  # for cor matrix if too many markers
  min_group_n          = 3L      # min samples per karyotype to compute means
)
# config override
if (exists("BP01_PARAMS_OVERRIDE") && is.list(BP01_PARAMS_OVERRIDE)) {
  for (k in names(BP01_PARAMS_OVERRIDE)) BP01_PARAMS[[k]] <- BP01_PARAMS_OVERRIDE[[k]]
}

# =============================================================================
# Core functions
# =============================================================================

# Compute per-group mean dosage and signed informativeness delta.
# X: markers × samples matrix
# groups: character vector of length ncol(X) with HOMO_1/HET/HOMO_2
compute_informative_markers <- function(X, groups, delta_min = 1.0,
                                          min_group_n = 3L) {
  h1_idx <- which(groups == "HOMO_1")
  h2_idx <- which(groups == "HOMO_2")
  he_idx <- which(groups == "HET")

  if (length(h1_idx) < min_group_n || length(h2_idx) < min_group_n) {
    return(list(
      h1_mean = rep(NA_real_, nrow(X)),
      h2_mean = rep(NA_real_, nrow(X)),
      het_mean = rep(NA_real_, nrow(X)),
      delta = rep(NA_real_, nrow(X)),
      informative = logical(nrow(X)),
      n_informative = 0L,
      reason = sprintf("insufficient group sizes (H1=%d, H2=%d)",
                       length(h1_idx), length(h2_idx))
    ))
  }

  h1_mean <- rowMeans(X[, h1_idx, drop = FALSE], na.rm = TRUE)
  h2_mean <- rowMeans(X[, h2_idx, drop = FALSE], na.rm = TRUE)
  het_mean <- if (length(he_idx) >= min_group_n)
    rowMeans(X[, he_idx, drop = FALSE], na.rm = TRUE)
  else rep(NA_real_, nrow(X))

  delta <- h2_mean - h1_mean
  delta[!is.finite(delta)] <- 0
  informative <- abs(delta) >= delta_min

  list(
    h1_mean = h1_mean, h2_mean = h2_mean, het_mean = het_mean,
    delta = delta,
    informative = informative,
    n_informative = sum(informative),
    reason = if (sum(informative) >= BP01_PARAMS$min_core_size)
      "ok" else sprintf("only %d informative markers", sum(informative))
  )
}

# Find the core block of high-correlation informative markers.
# Returns the SNP index vector defining the core block (indices into the
# informative-marker subset), plus the selected top-tier core for consensus.
find_core_block <- function(X, info_idx, rho_block = 0.7,
                             min_core = 20L, subsample_n = 5000L,
                             top_quantile = 0.75) {
  if (length(info_idx) < min_core) {
    return(list(block_idx = integer(0), core_idx = integer(0),
                reason = sprintf("only %d informative markers", length(info_idx))))
  }

  # Restrict to informative markers
  Xi <- X[info_idx, , drop = FALSE]
  n_m <- nrow(Xi)

  # Subsample for correlation if too many markers
  if (n_m > subsample_n) {
    sub_rows <- sort(sample(n_m, subsample_n))
  } else {
    sub_rows <- seq_len(n_m)
  }
  Xi_sub <- Xi[sub_rows, , drop = FALSE]

  # Mean-impute missing per marker for correlation
  for (mi in seq_len(nrow(Xi_sub))) {
    na_idx <- !is.finite(Xi_sub[mi, ])
    if (any(na_idx)) {
      m_val <- mean(Xi_sub[mi, !na_idx])
      if (is.finite(m_val)) Xi_sub[mi, na_idx] <- m_val
      else Xi_sub[mi, na_idx] <- 1
    }
  }
  Xi_sub[!is.finite(Xi_sub)] <- 1

  cor_mat <- suppressWarnings(cor(t(Xi_sub), use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- 0
  abs_cor <- abs(cor_mat)

  # Build adjacency and find connected components via BFS
  adj <- abs_cor >= rho_block
  diag(adj) <- FALSE
  n_sub <- nrow(adj)
  visited <- rep(FALSE, n_sub)
  comp_id <- rep(0L, n_sub)
  next_id <- 0L

  for (seed in seq_len(n_sub)) {
    if (visited[seed] || !any(adj[seed, ])) next
    next_id <- next_id + 1L
    queue <- seed
    while (length(queue) > 0) {
      curr <- queue[1]; queue <- queue[-1]
      if (visited[curr]) next
      visited[curr] <- TRUE
      comp_id[curr] <- next_id
      nbrs <- which(adj[curr, ] & !visited)
      queue <- c(queue, nbrs)
    }
  }

  # Largest component
  if (next_id == 0L) {
    return(list(block_idx = integer(0), core_idx = integer(0),
                reason = sprintf("no component reaches rho_block=%.2f", rho_block)))
  }
  comp_sizes <- tabulate(comp_id[comp_id > 0])
  largest_id <- which.max(comp_sizes)
  if (max(comp_sizes) < min_core) {
    return(list(block_idx = integer(0), core_idx = integer(0),
                reason = sprintf("largest component = %d < %d", max(comp_sizes), min_core)))
  }

  # Map back to full informative-marker indexing
  sub_in_block <- which(comp_id == largest_id)
  block_in_info <- sub_rows[sub_in_block]

  # Top-tier core: markers with highest mean in-block absolute correlation
  # (recompute just for the block)
  block_cor_sub <- abs_cor[sub_in_block, sub_in_block, drop = FALSE]
  mean_cor <- rowMeans(block_cor_sub, na.rm = TRUE)
  cor_thr <- quantile(mean_cor, top_quantile, na.rm = TRUE)
  core_in_sub <- sub_in_block[mean_cor >= cor_thr]
  core_in_info <- sub_rows[core_in_sub]

  list(
    block_idx = info_idx[block_in_info],
    core_idx = info_idx[core_in_info],
    n_block = length(block_in_info),
    n_core = length(core_in_info),
    mean_block_cor = round(mean(block_cor_sub, na.rm = TRUE), 4),
    reason = "ok"
  )
}

# Extend block boundaries outward one SNP at a time.
# Returns extended left/right marker indices (into full X) and the correlations
# at the stopping points.
extend_boundaries <- function(X, block_idx, core_idx,
                               rho_ext = 0.6, max_extend = 5000L) {
  if (length(block_idx) == 0L) return(NULL)

  n_total <- nrow(X)
  orig_left  <- min(block_idx)
  orig_right <- max(block_idx)

  # Per-sample core consensus: mean dosage over core markers
  core_mat <- X[core_idx, , drop = FALSE]
  core_consensus <- colMeans(core_mat, na.rm = TRUE)

  test_cor <- function(i) {
    m_gt <- X[i, ]
    valid <- is.finite(m_gt) & is.finite(core_consensus)
    if (sum(valid) < 20) return(NA_real_)
    r <- suppressWarnings(cor(m_gt[valid], core_consensus[valid]))
    if (!is.finite(r)) NA_real_ else abs(r)
  }

  # Extend LEFT
  new_left <- orig_left
  left_cor_stop <- NA_real_
  if (orig_left > 1) {
    for (mi in seq(orig_left - 1L, max(1L, orig_left - max_extend), by = -1L)) {
      r <- test_cor(mi)
      if (is.na(r) || r < rho_ext) { left_cor_stop <- r; break }
      new_left <- mi
    }
  }

  # Extend RIGHT
  new_right <- orig_right
  right_cor_stop <- NA_real_
  if (orig_right < n_total) {
    for (mi in seq(orig_right + 1L, min(n_total, orig_right + max_extend), by = 1L)) {
      r <- test_cor(mi)
      if (is.na(r) || r < rho_ext) { right_cor_stop <- r; break }
      new_right <- mi
    }
  }

  list(
    new_left = new_left, new_right = new_right,
    orig_left = orig_left, orig_right = orig_right,
    extended_left_snps  = orig_left  - new_left,
    extended_right_snps = new_right  - orig_right,
    left_cor_stop = left_cor_stop,
    right_cor_stop = right_cor_stop,
    core_consensus = core_consensus
  )
}

# =============================================================================
# Main per-candidate routine
# =============================================================================
process_candidate <- function(cid, chr, c_start, c_end) {
  # Ensure the candidate directory exists in evidence_registry (auto via
  # register_candidate if needed). reg$evidence$write_block will create its
  # own subdirs.
  if (!is.null(reg$intervals$get_candidate(cid))) {
    # candidate already registered — good
  } else {
    message("[01_dosage_signal] cid=", cid,
            " not in interval_registry; registering now.")
    reg$intervals$add_candidate(candidate_id = cid, chrom = chr,
                                 start_bp = c_start, end_bp = c_end)
    reg$evidence$register_candidate(cid, chrom = chr,
                                      start_bp = c_start, end_bp = c_end,
                                      tier = NA_integer_, score = NA_real_)
  }

  # Resolve karyotype groups via registry (replaces reading
  # candidate_pca_rotated.tsv from a followup folder).
  grs <- tryCatch(reg$samples$get_groups_for_candidate(cid),
                  error = function(e) NULL)
  if (is.null(grs) ||
      length(grs$HOM_REF) == 0L ||
      length(grs$HET)     == 0L ||
      length(grs$HOM_INV) == 0L) {
    message("[01_dosage_signal] cid=", cid,
            " SKIP — karyotype groups not registered (need inv_",
            cid, "_HOM_REF/HET/HOM_INV)")
    return(invisible(NULL))
  }
  # Build the (sample, group) map in the internal HOMO_1/HET/HOMO_2 convention
  # used by the downstream algorithm code. HOM_REF → HOMO_1, HOM_INV → HOMO_2.
  rot <- rbindlist(list(
    data.table(sample = grs$HOM_REF, coarse_group_refined = "HOMO_1"),
    data.table(sample = grs$HET,     coarse_group_refined = "HET"),
    data.table(sample = grs$HOM_INV, coarse_group_refined = "HOMO_2")
  ))

  dos_file   <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[01_dosage_signal] cid=", cid,
            " SKIP — dosage files missing under ", DOSAGE_DIR)
    return(invisible(NULL))
  }

  dos   <- fread(dos_file)
  sites <- fread(sites_file)

  # Use full region with some flank for extension. Flank = 30% of span on each side.
  span <- c_end - c_start
  flank <- max(500000L, as.integer(span * 0.3))
  scan_start <- max(1L, c_start - flank)
  scan_end   <- c_end + flank

  keep <- which(sites$pos >= scan_start & sites$pos <= scan_end)
  if (length(keep) < 200) {
    message("[01_dosage_signal] cid=", cid, " SKIP — only ", length(keep), " markers in scan region")
    return(invisible(NULL))
  }
  dos_reg   <- dos[keep]
  sites_reg <- sites[keep]

  # Rename Ind0/Ind1/... to real sample names if needed
  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot)) {
    setnames(dos_reg, old = sc, new = rot$sample)
  }
  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 20) {
    message("[01_dosage_signal] cid=", cid, " SKIP — <20 matched samples")
    return(invisible(NULL))
  }

  X <- as.matrix(dos_reg[, ..sample_cols])
  storage.mode(X) <- "double"
  groups <- rot[match(sample_cols, sample), coarse_group_refined]

  message("[01_dosage_signal] cid=", cid, " ", chr, " scan=", scan_start, "-", scan_end,
          " markers=", nrow(X), " samples=", ncol(X))

  # --- Step 1: informative markers ------------------------------------------
  im <- compute_informative_markers(X, groups,
                                      delta_min   = BP01_PARAMS$delta_min,
                                      min_group_n = BP01_PARAMS$min_group_n)
  info_idx <- which(im$informative)
  message("[01_dosage_signal]   informative markers: ", im$n_informative,
          " / ", nrow(X), "  (delta_min=", BP01_PARAMS$delta_min, ")")

  # Write informative markers TSV (one row per marker in the scan region,
  # including non-informative — useful for the diagnostic figure). This stays
  # as a raw artifact under the evidence_registry per-candidate dir, because
  # the per-marker table is too large to live as a JSON block or manifest row.
  cand_raw_dir <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                             "evidence_registry", "per_candidate", cid, "raw")
  ensure_dir(cand_raw_dir)

  markers_dt <- data.table(
    candidate_id = cid, chrom = chr,
    pos          = sites_reg$pos,
    marker_idx   = seq_len(nrow(X)),
    h1_mean      = round(im$h1_mean, 4),
    het_mean     = round(im$het_mean, 4),
    h2_mean      = round(im$h2_mean, 4),
    delta        = round(im$delta, 4),
    abs_delta    = round(abs(im$delta), 4),
    informative  = im$informative
  )
  markers_tsv <- file.path(cand_raw_dir, "dosage_informative_markers.tsv.gz")
  fwrite(markers_dt, markers_tsv, sep = "\t", compress = "gzip")

  # Sidecar RDS with the same informative-SNP info in machine-friendly form
  # (added 2026-04-24 chat D). Lets downstream consumers (GC detector, the
  # audit path that cross-checks phase-6 vs GC-detector diagnostic masks per
  # CODE_REVIEW_TO_FIGURE_OUT.md §5, etc.) read the informative-SNP set
  # directly without re-parsing the TSV. Also carries the per-informative-
  # SNP group means, which the polarity-aware version of the GC detector
  # scan (review §7) needs in order to compare each sample's dosage at SNP i
  # against the per-SNP group means rather than against fixed 0.5/1.5
  # thresholds.
  markers_rds <- file.path(cand_raw_dir, "dosage_informative_markers.rds")
  saveRDS(list(
    candidate_id          = cid,
    chrom                 = chr,
    scan_start_bp         = as.integer(scan_start),
    scan_end_bp           = as.integer(scan_end),
    delta_min             = as.numeric(BP01_PARAMS$delta_min),
    n_markers_in_scan     = as.integer(nrow(X)),
    n_informative         = as.integer(im$n_informative),
    # Indices into the scan-region marker order (1-based, same as marker_idx
    # column in the TSV). Use sites_reg$pos[informative_marker_idx] to
    # recover positions.
    informative_marker_idx = as.integer(info_idx),
    # Parallel bp-position vector for the informative set, ascending
    informative_pos_bp    = as.integer(sites_reg$pos[info_idx]),
    # Per-informative-SNP group means — the same values the polarity-aware
    # comparison needs. Indexed 1:1 with informative_marker_idx.
    h1_mean_informative   = as.numeric(im$h1_mean[info_idx]),
    h2_mean_informative   = as.numeric(im$h2_mean[info_idx]),
    het_mean_informative  = as.numeric(im$het_mean[info_idx]),
    delta_informative     = as.numeric(im$delta[info_idx]),
    # Convenience: sign of delta per informative SNP. +1 means "HOMO_2 > HOMO_1"
    # (so HOMO_1-class samples expected low dosage here); -1 means the other
    # polarity (HOMO_1-class samples expected high dosage here). Review §7's
    # polarity bug is specifically about detectors that assume this is always
    # +1. Saving it so the fix can read it straight.
    polarity_informative  = as.integer(sign(im$delta[info_idx]))
  ), markers_rds)

  if (im$n_informative < BP01_PARAMS$min_core_size) {
    message("[01_dosage_signal]   insufficient informative markers; writing status block")
    reg$evidence$write_block(cid, "dosage_blocks", list(
      status = "insufficient_informative_markers",
      reason = im$reason,
      n_informative = im$n_informative,
      delta_min = BP01_PARAMS$delta_min,
      markers_tsv_path = markers_tsv,
      markers_rds_path = markers_rds
    ))
    return(invisible(NULL))
  }

  # --- Step 2: find core block ---------------------------------------------
  cb <- find_core_block(X, info_idx,
                         rho_block     = BP01_PARAMS$rho_block,
                         min_core      = BP01_PARAMS$min_core_size,
                         subsample_n   = BP01_PARAMS$subsample_cor_n,
                         top_quantile  = BP01_PARAMS$top_core_quantile)
  if (length(cb$block_idx) == 0L) {
    message("[01_dosage_signal]   no core block: ", cb$reason)
    reg$evidence$write_block(cid, "dosage_blocks", list(
      status = "no_core_block",
      reason = cb$reason,
      markers_tsv_path = markers_tsv,
      markers_rds_path = markers_rds
    ))
    return(invisible(NULL))
  }
  message("[01_dosage_signal]   core block: ", cb$n_block, " markers; top-tier core: ",
          cb$n_core, "; mean_block_cor=", cb$mean_block_cor)

  # --- Step 3: extend boundaries -------------------------------------------
  ext <- extend_boundaries(X, cb$block_idx, cb$core_idx,
                            rho_ext     = BP01_PARAMS$rho_ext,
                            max_extend  = BP01_PARAMS$max_extend_snps)

  if (is.null(ext)) {
    message("[01_dosage_signal]   extension failed")
    reg$evidence$write_block(cid, "dosage_blocks", list(
      status = "extension_failed",
      markers_tsv_path = markers_tsv,
      markers_rds_path = markers_rds
    ))
    return(invisible(NULL))
  }

  # --- Step 4: write structured block via registry --------------------------
  block_data <- list(
    status = "ok",
    # before extension (core block)
    core_left_marker_idx   = ext$orig_left,
    core_right_marker_idx  = ext$orig_right,
    core_left_bp           = as.integer(sites_reg$pos[ext$orig_left]),
    core_right_bp          = as.integer(sites_reg$pos[ext$orig_right]),
    core_n_markers         = cb$n_block,
    core_top_tier_n        = cb$n_core,
    core_mean_block_cor    = round(cb$mean_block_cor, 6),
    # after extension — these are the BLOCK-extension breakpoint estimates
    # consumed by 03_consensus_merge (as one of the per-method sources)
    ext_left_marker_idx    = ext$new_left,
    ext_right_marker_idx   = ext$new_right,
    ext_left_bp            = as.integer(sites_reg$pos[ext$new_left]),
    ext_right_bp           = as.integer(sites_reg$pos[ext$new_right]),
    extended_left_snps     = ext$extended_left_snps,
    extended_right_snps    = ext$extended_right_snps,
    left_cor_stop          = round(ext$left_cor_stop, 4),
    right_cor_stop         = round(ext$right_cor_stop, 4),
    # parameters used (for reproducibility)
    delta_min              = BP01_PARAMS$delta_min,
    rho_block              = BP01_PARAMS$rho_block,
    rho_ext                = BP01_PARAMS$rho_ext,
    # input candidate for comparison
    input_start_bp         = as.integer(c_start),
    input_end_bp           = as.integer(c_end),
    shift_left_kb          = round((sites_reg$pos[ext$new_left]  - c_start) / 1000, 2),
    shift_right_kb         = round((sites_reg$pos[ext$new_right] - c_end)   / 1000, 2),
    extended_span_bp       = as.integer(sites_reg$pos[ext$new_right] -
                                          sites_reg$pos[ext$new_left]),
    markers_tsv_path       = markers_tsv,
    markers_rds_path       = markers_rds
  )
  reg$evidence$write_block(cid, "dosage_blocks", block_data)

  message("[01_dosage_signal]   done. Extended block: ",
          format(sites_reg$pos[ext$new_left],  big.mark = ","), " - ",
          format(sites_reg$pos[ext$new_right], big.mark = ","),
          " (", round((sites_reg$pos[ext$new_right] - sites_reg$pos[ext$new_left]) / 1e6, 3),
          " Mb)")

  invisible(list(markers_dt = markers_dt, block_data = block_data, ext = ext))
}

# =============================================================================
# Driver
# =============================================================================
main <- function() {
  # Candidates come from interval_registry, not a flat CANDIDATE_TABLE file.
  all_cands <- reg$intervals$get_windows()  # returns NULL if not windows API
  # Actually, the candidate list lives in interval_registry$candidate_intervals
  # via list_candidates() if available. Fallback to bulk read:
  cand_path <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "interval_registry", "candidate_intervals.tsv")
  if (!file.exists(cand_path)) {
    stop("[01_dosage_signal] candidate_intervals.tsv not found at ", cand_path,
         ". Has interval_registry been populated by phase 4a (STEP_C01d)?")
  }
  cand <- fread(cand_path)
  if (!is.na(cid_filter)) {
    cand <- cand[candidate_id == cid_filter]
  }
  if (nrow(cand) == 0) {
    message("[01_dosage_signal] No candidates to process",
            if (!is.na(cid_filter)) paste0(" (filter: ", cid_filter, ")") else "")
    return(invisible(NULL))
  }
  message("[01_dosage_signal] Processing ", nrow(cand), " candidate(s) from interval_registry")

  for (ci in seq_len(nrow(cand))) {
    row <- cand[ci]
    cid <- as.character(row$candidate_id)
    chr <- as.character(row$chrom)
    c_start <- as.numeric(row$start_bp)
    c_end   <- as.numeric(row$end_bp)
    tryCatch(
      process_candidate(cid, chr, c_start, c_end),
      error = function(e) {
        message("[01_dosage_signal] cid=", cid, " ERROR: ", conditionMessage(e))
      }
    )
  }
  message("[01_dosage_signal] DONE")
}

main()
