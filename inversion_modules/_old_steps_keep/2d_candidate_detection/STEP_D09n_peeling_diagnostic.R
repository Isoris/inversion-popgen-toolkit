#!/usr/bin/env Rscript

# =============================================================================
# STEP_D09n_peeling_diagnostic.R — Blockwise peeling diagnostic (sample zeroing)
#
# SAMPLE-PEELING DIAGNOSTIC COMPANION TO DECOMPOSITION
#
# =============================================================================
#
# IMPORTANT DISTINCTIONS:
#   - Sample peeling != matrix treatment (transforms A1-A6 change VALUES)
#   - Sample peeling zeroes out specific fish, then re-examines local blocks
#   - This step is a DIAGNOSTIC complement to C01i decomposition
#   - Block disappearance after peeling is NOT automatically proof of artifact
#     because: (a) same fish may carry both parent and child inversions
#              (b) power drops when you remove samples
#              (c) signal may be real but sample-dependent
#
# MATHEMATICAL APPROACH — ZERO-MASKING:
#   The precomp dt has PC_<k>_<sample> columns = per-sample PCA loading
#   scores at each window. Each window's loading vector is:
#
#     v_i = [PC_1_fish1, PC_1_fish2, ..., PC_1_fish226]
#
#   To "peel" fish47: multiply its entry by 0 in the loading vector.
#   The angular distance between windows recomputes instantly from the
#   masked vectors. No matrix rebuild ceremony — just multiply and go.
#
#     v_i_peeled = v_i * mask    where mask[fish47] = 0, all others = 1
#     d_angular(i,j) = arccos( v_i_peeled · v_j_peeled / (|v_i_peeled| |v_j_peeled|) )
#     sim(i,j) = 1 - d_angular / dmax
#
#   Optionally, mask can be SOFT (0 < w < 1) for partial downweighting,
#   e.g. w=0.5 for second-degree relatives. But default is hard 0/1.
#
# THREE PEEL LEVELS:
#
#   Level 1 — Genome-wide relatedness peel (ngsRelate):
#     Remove first-degree relatives from whole-genome kinship estimates.
#     Useful baseline but COARSE: genome-wide kinship averages over 28
#     chromosomes. Two fish can be 3 generations apart genome-wide
#     but share a WHOLE CHROMOSOME identical-by-descent. This peel
#     misses chromosome-local family sharing.
#
#   Level 2 — Chromosome-local relatedness peel (from sim_mat):
#     MORE INTERESTING. Identify fish that co-segregate locally on THIS
#     chromosome by clustering sample PC1 profiles within the block.
#     Fish in the same local cluster share haplotype structure here
#     regardless of genome-wide kinship. This catches 3-generation-distant
#     fish that inherited the same chromosome segment by chance.
#     With Ne~20, this happens constantly.
#
#   Level 3 — Block-aware driver peel (from C01i/PC1 leverage):
#     For a specific block, identify the fish driving that local structure
#     (extreme PC1 loading / HOM_INV from C01i decomposition) and test
#     removal specifically for that block.
#
# INPUTS:
#   --precomp_dir <dir>     — precomp/<chr>.precomp.rds files
#   --blocks <file>         — decomposed block table (staircase or C01i)
#   --pruned_list <file>    — pruned_samples.txt OR per-chr pruning table (Level 1)
#   --coseg_dir <dir>       — C01i output (Level 3, optional)
#   --chr <chr>             — chromosome to process
#   --outdir <dir>
#   [--kin_threshold 0.2]   — kinship threshold for Level 1 (unused if pruned list)
#   [--local_kin_thresh 0.7] — PC1 correlation threshold for L1b
#   [--max_peel_frac 0.3]   — never zero out more than 30% of samples
#   [--local_k 5]           — k for local clustering (Level 2)
#   [--n_pcs 2]             — PCs to use from loading vectors
#
# OUTPUTS:
#   peel_diagnostic_<chr>.tsv        — per-block comparison table
#   peel_sample_groups_<chr>.tsv     — which samples were zeroed per block
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# CLI PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
precomp_dir <- NULL; blocks_file <- NULL; pruned_list <- NULL
coseg_dir <- NULL; chr <- NULL; outdir <- "peel_diagnostic"
kin_threshold <- 0.2; local_kin_thresh <- 0.7
max_peel_frac <- 0.3; local_k <- 5L; n_pcs <- 2L
window_size_bp <- 50000L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp_dir" && i < length(args))      { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--blocks" && i < length(args))       { blocks_file <- args[i+1]; i <- i+2 }
  else if (a == "--pruned_list" && i < length(args))  { pruned_list <- args[i+1]; i <- i+2 }
  else if (a == "--relatedness" && i < length(args))  { pruned_list <- args[i+1]; i <- i+2 }  # backward compat
  else if (a == "--coseg_dir" && i < length(args))    { coseg_dir <- args[i+1]; i <- i+2 }
  else if (a == "--chr" && i < length(args))          { chr <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))       { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--kin_threshold" && i < length(args)) { kin_threshold <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--local_kin_thresh" && i < length(args)) { local_kin_thresh <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--max_peel_frac" && i < length(args)) { max_peel_frac <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--local_k" && i < length(args))      { local_k <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--n_pcs" && i < length(args))        { n_pcs <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--window_size_bp" && i < length(args)) { window_size_bp <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# CORE: Build sample mask vector
# =============================================================================
# mask = numeric vector of length n_samples, values in [0, 1].
# 0 = fully peeled, 1 = fully present, 0.5 = half-weighted.

build_mask <- function(sample_names, peel_samples, soft_weight = 0) {
  mask <- rep(1, length(sample_names))
  names(mask) <- sample_names
  peel_idx <- match(peel_samples, sample_names)
  peel_idx <- peel_idx[!is.na(peel_idx)]
  mask[peel_idx] <- soft_weight
  mask
}


# =============================================================================
# CORE: Extract loading matrix for a window range, apply mask, compute sim_mat
# =============================================================================
# This is the whole operation:
#   1. Extract PC_k_<sample> columns for windows [start, end]
#   2. Multiply each sample column by its mask value (0 = zeroed out)
#   3. Compute angular distance between window loading vectors
#   4. Convert to similarity

masked_local_simmat <- function(dt, start_win, end_win, sample_names,
                                 mask, n_pcs = 2L) {
  win_idx <- start_win:end_win
  n_win <- length(win_idx)
  if (n_win < 3) return(NULL)

  # Build loading matrix: rows = windows, cols = samples (across all PCs)
  load_mat <- NULL
  for (k in seq_len(n_pcs)) {
    pc_cols <- paste0("PC_", k, "_", sample_names)
    pc_cols <- intersect(pc_cols, names(dt))
    if (length(pc_cols) == 0) next

    sub <- as.matrix(dt[win_idx, ..pc_cols])

    # Apply mask: multiply each column by its sample's mask weight
    col_samples <- sub(paste0("^PC_", k, "_"), "", pc_cols)
    for (j in seq_along(pc_cols)) {
      sname <- col_samples[j]
      if (sname %in% names(mask)) {
        sub[, j] <- sub[, j] * mask[sname]
      }
    }

    load_mat <- if (is.null(load_mat)) sub else cbind(load_mat, sub)
  }

  if (is.null(load_mat) || ncol(load_mat) < 5) return(NULL)

  # Impute remaining NAs with 0 (masked samples are already 0)
  load_mat[is.na(load_mat)] <- 0

  # Remove zero-variance columns (fully masked or constant)
  col_var <- apply(load_mat, 2, var, na.rm = TRUE)
  keep_cols <- which(col_var > 1e-15)
  if (length(keep_cols) < 3) return(NULL)
  load_mat <- load_mat[, keep_cols, drop = FALSE]

  # Normalize rows to unit vectors
  row_norms <- sqrt(rowSums(load_mat^2))
  row_norms[row_norms == 0] <- 1
  load_norm <- load_mat / row_norms

  # Cosine similarity → angular distance → sim_mat
  cos_sim <- tcrossprod(load_norm)
  cos_sim <- pmin(pmax(cos_sim, -1), 1)
  d_angular <- acos(cos_sim) / pi
  d_angular[!is.finite(d_angular)] <- 0
  diag(d_angular) <- 0

  dmax <- quantile(d_angular[d_angular > 0], 0.95, na.rm = TRUE)
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(d_angular / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1

  sim
}


# =============================================================================
# BLOCK METRICS from local sim_mat
# =============================================================================

compute_block_metrics <- function(sim_local) {
  n <- nrow(sim_local)
  if (is.null(sim_local) || n < 4) {
    return(list(inside_mean = NA, contrast = NA,
                squareness = NA, occupancy = NA))
  }

  ut <- sim_local[upper.tri(sim_local)]
  ut <- ut[is.finite(ut)]
  inside_mean <- if (length(ut) > 0) mean(ut) else NA_real_

  # Squareness (far/near ratio)
  bw <- n
  near_limit <- max(1L, round(bw * 0.2))
  far_start  <- max(1L, round(bw * 0.8))
  near_vals <- c(); far_vals <- c()
  for (ii in 1:(bw-1)) for (jj in (ii+1):bw) {
    d <- jj - ii; v <- sim_local[ii, jj]
    if (!is.finite(v)) next
    if (d <= near_limit) near_vals <- c(near_vals, v)
    if (d >= far_start) far_vals <- c(far_vals, v)
  }
  near_med <- if (length(near_vals) > 0) median(near_vals) else NA
  far_med  <- if (length(far_vals) > 0)  median(far_vals)  else NA
  squareness <- if (!is.na(near_med) && near_med > 0) far_med / near_med else NA

  # Contrast: core vs edge bands within the sub-matrix
  edge_band <- min(3L, n %/% 4)
  if (edge_band >= 1 && n > edge_band * 2) {
    core_vals <- sim_local[(edge_band+1):(n-edge_band), (edge_band+1):(n-edge_band)]
    core_ut <- core_vals[upper.tri(core_vals)]
    core_ut <- core_ut[is.finite(core_ut)]
    edge_vals <- c(sim_local[1:edge_band, (edge_band+1):n],
                   sim_local[(n-edge_band+1):n, 1:(n-edge_band)])
    edge_vals <- edge_vals[is.finite(edge_vals)]
    core_mean <- if (length(core_ut) > 0) mean(core_ut) else inside_mean
    edge_mean <- if (length(edge_vals) > 0) mean(edge_vals) else inside_mean
    contrast <- core_mean - edge_mean
  } else {
    contrast <- NA_real_
  }

  # Occupancy
  if (length(ut) > 10) {
    bg_level <- quantile(ut, 0.25, na.rm = TRUE)
    occupancy <- sum(ut > bg_level + 0.02) / length(ut)
  } else {
    occupancy <- NA_real_
  }

  list(inside_mean = round(inside_mean, 4),
       contrast    = round(contrast, 4),
       squareness  = round(squareness, 4),
       occupancy   = round(occupancy, 4))
}


# =============================================================================
# LEVEL 1: Pruned-list peel
# =============================================================================
# Accepts TWO input formats:
#
# Format A — genome-wide pruned_samples.txt (81 IDs, one per line):
#   The existing file from 06_relatedness/pruned_samples.txt.
#   These 81 are the KEPT samples. Peel = everyone NOT in this list.
#   COARSE: genome-wide ngsRelate averages over 28 chromosomes.
#
# Format B — per-chromosome pruning table (TSV with chr + sample columns):
#   chr          sample         status
#   C_gar_LG01   CGA_001        keep
#   C_gar_LG01   CGA_002        remove
#   C_gar_LG01   CGA_003        keep
#   C_gar_LG02   CGA_001        keep
#   C_gar_LG02   CGA_002        keep    <- kept on LG02 but removed on LG01!
#
#   With Ne~20, two fish can be genome-wide unrelated but share a WHOLE
#   chromosome by IBD. Per-chromosome pruning catches this. Different
#   chromosomes get different "keep" lists (LG01 might keep 90, LG02 only 70).
#
# Auto-detection: if input is a 1-column file (just sample IDs) → Format A.
#   If it has chr/sample columns → Format B, filter to current chromosome.
#   If neither exists → skip L1 entirely.

build_pruned_peel_set <- function(pruned_input, sample_names, chr_name,
                                   max_peel_frac = 0.3) {
  if (is.null(pruned_input) || !file.exists(pruned_input)) {
    message("  [L1] No pruned list — skipping")
    return(NULL)
  }

  # Try reading as TSV first
  raw <- tryCatch(fread(pruned_input, header = TRUE), error = function(e) NULL)

  if (!is.null(raw) && ncol(raw) >= 2) {
    # ---- Format B: per-chromosome table ----
    names(raw) <- tolower(names(raw))

    # Find chr column
    chr_col <- intersect(c("chr", "chrom", "chromosome"), names(raw))[1]
    samp_col <- intersect(c("sample", "sample_id", "id", "ida"), names(raw))[1]

    if (!is.na(chr_col) && !is.na(samp_col)) {
      # Filter to this chromosome
      chr_rows <- raw[get(chr_col) == chr_name]
      if (nrow(chr_rows) == 0) {
        # Try without prefix
        chr_clean <- sub("^C_gar_", "", chr_name)
        chr_rows <- raw[get(chr_col) == chr_clean | get(chr_col) == chr_name]
      }

      # If table has prune_method column (from 11b), prefer ngsRelate rows
      if ("prune_method" %in% names(chr_rows) && nrow(chr_rows) > 0) {
        ngsrel_rows <- chr_rows[prune_method == "ngsRelate"]
        if (nrow(ngsrel_rows) > 0) {
          chr_rows <- ngsrel_rows
          message("  [L1] Using ngsRelate-based pruning for ", chr_name)
        }
        # else fall through to whatever rows exist
      }

      if (nrow(chr_rows) > 0) {
        # If there's a status column, use it
        status_col <- intersect(c("status", "keep", "action"), names(chr_rows))[1]
        if (!is.na(status_col)) {
          kept <- chr_rows[get(status_col) %in% c("keep", "kept", "TRUE", "1")][[samp_col]]
        } else {
          # All listed samples are "kept" (pruned = complement)
          kept <- chr_rows[[samp_col]]
        }

        kept <- intersect(kept, sample_names)
        to_remove <- setdiff(sample_names, kept)

        if (length(to_remove) > length(sample_names) * max_peel_frac) {
          message("  [L1] Per-chr pruning would remove ", length(to_remove),
                  " (>", max_peel_frac*100, "%) — clamping")
          to_remove <- to_remove[seq_len(floor(length(sample_names) * max_peel_frac))]
        }

        message("  [L1] Per-chromosome pruned list: keep ", length(kept),
                ", peel ", length(to_remove), " of ", length(sample_names),
                " for ", chr_name)
        return(to_remove)
      }
    }
  }

  # ---- Format A: genome-wide pruned_samples.txt (one ID per line) ----
  kept_lines <- tryCatch(readLines(pruned_input), error = function(e) NULL)
  if (is.null(kept_lines)) return(NULL)

  # Clean: strip paths, extensions
  kept <- sapply(kept_lines, function(x) {
    x <- trimws(x)
    x <- basename(x)
    x <- sub("\\.sorted\\.markdup\\.bam$", "", x)
    x <- sub("\\.bam$", "", x)
    x
  }, USE.NAMES = FALSE)
  kept <- kept[nchar(kept) > 0]
  kept <- intersect(kept, sample_names)

  if (length(kept) == 0) return(NULL)

  to_remove <- setdiff(sample_names, kept)
  message("  [L1] Genome-wide pruned list: keep ", length(kept),
          ", peel ", length(to_remove), " of ", length(sample_names))
  to_remove
}


# =============================================================================
# LEVEL 1b: Chromosome-local relatedness from PC1 (FAST, no ngsRelate)
# =============================================================================
# Cheap alternative to per-chromosome ngsRelate. Uses the WHOLE chromosome
# PC1 trajectories (not just one block) to identify fish that co-segregate
# across THIS chromosome. No external file needed — computed from precomp dt.
#
# Two fish whose PC1 profiles correlate across ALL windows of this chromosome
# = they share haplotype structure chromosome-wide = local relatives.
# This catches the "3 generations apart genome-wide but share a whole LG" case.
#
# Cost: one cor() call on 226×N_windows matrix. Seconds.

build_chrwide_local_peel_set <- function(dt, sample_names, local_kin_thresh = 0.7,
                                          max_peel_frac = 0.3) {
  pc1_cols <- paste0("PC_1_", sample_names)
  pc1_cols <- intersect(pc1_cols, names(dt))
  if (length(pc1_cols) < 20) return(NULL)

  # Full chromosome PC1 matrix: rows = windows, cols = samples
  # Subsample windows for speed if chromosome is huge
  n_win <- nrow(dt)
  if (n_win > 2000) {
    win_idx <- sort(sample(n_win, 2000))
  } else {
    win_idx <- seq_len(n_win)
  }

  sub <- as.matrix(dt[win_idx, ..pc1_cols])
  sub[is.na(sub)] <- 0
  col_samples <- sub("^PC_1_", "", colnames(sub))

  # Sample × sample correlation across ALL chromosome windows
  samp_cor <- cor(sub, use = "pairwise.complete.obs")
  samp_cor[!is.finite(samp_cor)] <- 0
  diag(samp_cor) <- 0  # exclude self

  # For each sample, check if it has a "local relative" (|cor| > threshold)
  # Greedy removal: remove the most-connected sample per round
  adj <- abs(samp_cor) >= local_kin_thresh
  to_remove <- character(0)
  max_remove <- floor(length(sample_names) * max_peel_frac)

  while (any(adj) && length(to_remove) < max_remove) {
    conn <- colSums(adj)
    if (max(conn) == 0) break
    worst_idx <- which.max(conn)
    worst_name <- col_samples[worst_idx]
    to_remove <- c(to_remove, worst_name)
    adj[worst_idx, ] <- FALSE
    adj[, worst_idx] <- FALSE
  }

  to_remove <- intersect(to_remove, sample_names)
  if (length(to_remove) == 0) return(NULL)

  message("  [L1b] Chr-local PC1 relatedness peel: ", length(to_remove),
          " of ", length(sample_names), " (cor threshold=", local_kin_thresh, ")")
  to_remove
}


# =============================================================================
# LEVEL 2: Chromosome-local co-segregation peel
# =============================================================================
# THE INTERESTING ONE. Genome-wide kinship misses local IBD sharing.
# Two fish 3 generations apart genome-wide can share an entire chromosome
# segment IBD. In Ne~20, this happens constantly — bad luck inheritance.
#
# Method:
#   1. Extract PC1 scores for all samples across the block's windows
#   2. Each sample is a trajectory: [PC1 at window_1, PC1 at window_2, ...]
#   3. Correlate trajectories between all sample pairs
#   4. Fish whose trajectories are highly correlated = locally co-inherited
#   5. Cluster into local_k groups by trajectory similarity
#   6. Peel the dominant cluster (largest local family group)

build_local_peel_set <- function(dt, start_win, end_win, sample_names,
                                  local_k = 5L, max_peel_frac = 0.3) {
  pc1_cols <- paste0("PC_1_", sample_names)
  pc1_cols <- intersect(pc1_cols, names(dt))
  if (length(pc1_cols) < 10) return(NULL)

  # Extract: rows = windows in block, cols = samples
  sub <- as.matrix(dt[start_win:end_win, ..pc1_cols])
  col_samples <- sub("^PC_1_", "", colnames(sub))

  # Transpose: rows = samples, cols = windows
  sub_t <- t(sub)
  sub_t[is.na(sub_t)] <- 0

  # Sample × sample correlation of PC1 trajectories across the block
  samp_cor <- cor(t(sub_t), use = "pairwise.complete.obs")
  samp_cor[!is.finite(samp_cor)] <- 0

  # Cluster samples by local co-segregation
  k_use <- min(local_k, length(col_samples) %/% 3)
  if (k_use < 2) return(NULL)

  d <- as.dist(1 - samp_cor)
  hc <- hclust(d, method = "ward.D2")
  groups <- cutree(hc, k = k_use)
  names(groups) <- col_samples

  # Find the LARGEST group = dominant local family cluster
  group_sizes <- table(groups)
  dominant_group <- as.integer(names(which.max(group_sizes)))
  dominant_samples <- names(groups[groups == dominant_group])

  # Clamp
  max_remove <- floor(length(sample_names) * max_peel_frac)
  if (length(dominant_samples) > max_remove) {
    dominant_samples <- sample(dominant_samples, max_remove)
  }

  message("  [L2] Local co-segregation peel: ", length(dominant_samples),
          " samples from dominant local cluster (k=", k_use,
          ", group ", dominant_group, " of ", length(group_sizes), ")")

  list(peel_samples = dominant_samples,
       group_labels = groups,
       dominant_group = dominant_group)
}


# =============================================================================
# LEVEL 3: Block-aware driver peel (PC1 leverage + C01i karyotype)
# =============================================================================

build_driver_peel_set <- function(dt, start_win, end_win, sample_names,
                                   coseg_samples = NULL, max_peel_frac = 0.3) {
  pc1_cols <- paste0("PC_1_", sample_names)
  pc1_cols <- intersect(pc1_cols, names(dt))
  if (length(pc1_cols) < 10) return(NULL)

  sub <- as.matrix(dt[start_win:end_win, ..pc1_cols])
  col_samples <- sub("^PC_1_", "", colnames(sub))

  # PC1 leverage: mean absolute loading across block windows
  leverage <- colMeans(abs(sub), na.rm = TRUE)
  names(leverage) <- col_samples

  top_n <- min(30L, floor(length(col_samples) * max_peel_frac))
  sorted <- sort(leverage, decreasing = TRUE)
  top_leverage <- names(sorted)[seq_len(top_n)]

  # Combine with C01i HOM_INV if available
  if (!is.null(coseg_samples) && nrow(coseg_samples) > 0) {
    inv_samps <- intersect(
      coseg_samples[group %in% c("HOMO_INV", "RARE_INV")]$sample,
      sample_names)
    drivers <- unique(c(inv_samps, top_leverage))
  } else {
    drivers <- top_leverage
  }

  max_remove <- floor(length(sample_names) * max_peel_frac)
  if (length(drivers) > max_remove) drivers <- drivers[seq_len(max_remove)]

  message("  [L3] Driver peel: ", length(drivers), " samples")
  drivers
}


# =============================================================================
# CLASSIFY PEEL EFFECT
# =============================================================================

classify_peel_effect <- function(before, after) {
  if (is.na(after$contrast) || is.na(before$contrast)) return("ambiguous")
  if (!is.finite(before$contrast) || before$contrast <= 0) return("ambiguous")

  ratio <- after$contrast / before$contrast

  if (after$contrast < 0.01 || ratio < 0.15) return("disappeared")
  if (ratio < 0.50)                           return("weakened")
  if (ratio >= 0.80)                          return("stable")

  # Squareness improvement = hidden sub-block revealed
  if (!is.na(before$squareness) && !is.na(after$squareness) &&
      after$squareness > before$squareness + 0.10 &&
      after$contrast > 0.02) {
    return("revealed_child")
  }

  "weakened"
}


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

run_peel_diagnostic <- function(precomp, blocks, pruned_list = NULL,
                                 coseg_dir = NULL, chr_name = NULL,
                                 kin_threshold = 0.2, local_kin_thresh = 0.7,
                                 max_peel_frac = 0.3,
                                 local_k = 5L, n_pcs_use = 2L) {

  dt <- precomp$dt
  N <- nrow(dt)
  chr_name <- chr_name %||% precomp$chrom

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  sample_names <- sub("^PC_1_", "", pc1_cols)
  n_samples <- length(sample_names)

  message("[PEEL] ", chr_name, " | ", N, " windows | ", n_samples,
          " samples | ", nrow(blocks), " blocks")

  # ---- L1 peel set: pruned list (genome-wide or per-chr) ----
  l1_peel <- build_pruned_peel_set(pruned_list, sample_names, chr_name,
                                    max_peel_frac)

  # ---- L1b peel set: chromosome-local PC1 relatedness (free, no file needed) ----
  l1b_peel <- build_chrwide_local_peel_set(dt, sample_names, local_kin_thresh,
                                            max_peel_frac)

  # ---- Load C01i sample classifications ----
  coseg_dt <- NULL
  if (!is.null(coseg_dir)) {
    cf <- file.path(coseg_dir, "marker_coseg_samples.tsv.gz")
    if (file.exists(cf)) {
      coseg_dt <- fread(cf)
      if ("chrom" %in% names(coseg_dt)) coseg_dt <- coseg_dt[chrom == chr_name]
      message("  Loaded C01i sample classes: ", nrow(coseg_dt))
    }
  }

  # ---- Process each block ----
  results <- list()
  peel_groups <- list()

  for (bi in seq_len(nrow(blocks))) {
    block <- blocks[bi]
    bid <- block$block_id
    bs <- block$start; be <- block$end; bw <- be - bs + 1
    if (bw < 5 || bs < 1 || be > N) next

    parent <- block$parent_id %||% NA_integer_

    message("  Block ", bid, " [", bs, "-", be, "] width=", bw)

    # ---- BEFORE: all-samples sim_mat ----
    mask_full <- build_mask(sample_names, character(0))
    sim_before <- masked_local_simmat(dt, bs, be, sample_names,
                                       mask_full, n_pcs_use)
    if (is.null(sim_before)) next
    m_before <- compute_block_metrics(sim_before)

    # Helper to run one peel comparison
    do_peel <- function(peel_samples, peel_mode, group_label) {
      if (length(peel_samples) < 1) return(NULL)
      if (length(peel_samples) > n_samples * max_peel_frac) return(NULL)

      mask <- build_mask(sample_names, peel_samples)
      sim_after <- masked_local_simmat(dt, bs, be, sample_names,
                                        mask, n_pcs_use)
      if (is.null(sim_after)) return(NULL)
      m_after <- compute_block_metrics(sim_after)
      effect <- classify_peel_effect(m_before, m_after)

      message("    ", peel_mode, ": -", length(peel_samples),
              " fish → ", effect, " (contrast ",
              m_before$contrast, " → ", m_after$contrast, ")")

      peel_groups[[length(peel_groups) + 1]] <<- data.table(
        block_id = bid, peel_mode = peel_mode,
        sample = peel_samples, action = "zeroed"
      )

      data.table(
        block_id          = bid,
        parent_id         = parent,
        chromosome        = chr_name,
        start_win         = bs, end_win = be,
        start_bp          = (bs - 1L) * window_size_bp,
        end_bp            = be * window_size_bp,
        peel_mode         = peel_mode,
        removed_sample_group = group_label,
        n_removed         = length(peel_samples),
        before_inside     = m_before$inside_mean,
        after_inside      = m_after$inside_mean,
        before_contrast   = m_before$contrast,
        after_contrast    = m_after$contrast,
        before_squareness = m_before$squareness,
        after_squareness  = m_after$squareness,
        before_occupancy  = m_before$occupancy,
        after_occupancy   = m_after$occupancy,
        effect_class      = effect
      )
    }

    # ---- L1: pruned list (genome-wide or per-chr) ----
    if (!is.null(l1_peel) && length(l1_peel) > 0) {
      r <- do_peel(l1_peel, "L1_pruned_list", "pruned_relatives")
      if (!is.null(r)) results[[length(results) + 1]] <- r
    }

    # ---- L1b: chromosome-local PC1 relatedness (fast, no external file) ----
    if (!is.null(l1b_peel) && length(l1b_peel) > 0) {
      r <- do_peel(l1b_peel, "L1b_chrlocal_kin", "chrwide_PC1_relatives")
      if (!is.null(r)) results[[length(results) + 1]] <- r
    }

    # ---- L2: chromosome-local co-segregation ----
    local_info <- build_local_peel_set(dt, bs, be, sample_names,
                                        local_k, max_peel_frac)
    if (!is.null(local_info)) {
      r <- do_peel(local_info$peel_samples, "L2_local_coseg",
                    paste0("local_cluster_", local_info$dominant_group))
      if (!is.null(r)) results[[length(results) + 1]] <- r
    }

    # ---- L3: block drivers (PC1 leverage + C01i) ----
    drivers <- build_driver_peel_set(dt, bs, be, sample_names,
                                      coseg_dt, max_peel_frac)
    if (!is.null(drivers) && length(drivers) > 0) {
      r <- do_peel(drivers, "L3_block_drivers", "high_leverage_plus_HOM_INV")
      if (!is.null(r)) results[[length(results) + 1]] <- r
    }

    # ---- L3b: HOM_INV only (if C01i available) ----
    if (!is.null(coseg_dt) && nrow(coseg_dt) > 0) {
      inv_only <- intersect(
        coseg_dt[group %in% c("HOMO_INV", "RARE_INV")]$sample,
        sample_names)
      if (length(inv_only) >= 2) {
        r <- do_peel(inv_only, "L3_hom_inv_only", "C01i_HOMO_INV")
        if (!is.null(r)) results[[length(results) + 1]] <- r
      }
    }
  }

  result_dt <- if (length(results) > 0) rbindlist(results, fill = TRUE) else data.table()
  groups_dt <- if (length(peel_groups) > 0) rbindlist(peel_groups, fill = TRUE) else data.table()

  list(diagnostics = result_dt, sample_groups = groups_dt)
}


# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================

print_interpretation <- function(result_dt) {
  if (nrow(result_dt) == 0) {
    message("\n  No peel diagnostics produced.")
    return(invisible(NULL))
  }

  message("\n=== PEEL DIAGNOSTIC RESULTS ===")
  message("  Total comparisons: ", nrow(result_dt))

  for (mode in unique(result_dt$peel_mode)) {
    sub <- result_dt[peel_mode == mode]
    tab <- table(sub$effect_class)
    message("\n  ", mode, " (", nrow(sub), " blocks):")
    for (cls in names(tab)) message("    ", cls, ": ", tab[cls])
  }

  message("\n  === INTERPRETATION ===")
  message("  L1 pruned_list: uses your existing pruned_samples.txt or per-chr table")
  message("    Accepts genome-wide (81 IDs) OR per-chromosome table (different per LG)")
  message("    stable → not driven by known relatives → structural")
  message("    disappeared → family structure artifact")
  message("")
  message("  L1b chrlocal_kin: FAST, no external file, computed from precomp PC1")
  message("    Correlates ALL sample PC1 trajectories across the whole chromosome.")
  message("    Catches 'bad luck inheritance' — fish unrelated genome-wide but")
  message("    sharing a whole LG segment by chance (constant with Ne~20).")
  message("    LG01 might remove 40 fish, LG02 might remove 25 — different per chr.")
  message("    stable → NOT driven by chr-local family sharing → real structure")
  message("    disappeared → chromosome-local family block, not inversion")
  message("")
  message("  L2 local_coseg: block-level clustering of PC1 trajectories")
  message("    Finer than L1b — clusters within the BLOCK interval only.")
  message("    stable → signal not from dominant local cluster → real")
  message("    revealed_child → hidden sub-block under family signal")
  message("")
  message("  L3 block_drivers / hom_inv_only:")
  message("    weakened → EXPECTED (you removed carriers)")
  message("    stable → signal from many fish → robust")
  message("")
  message("  REMEMBER: removing carriers of a real inversion weakens it too.")
  message("  L1b + L2 together give the strongest family-vs-structure diagnostic.")
}


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

if (!interactive()) {
  if (is.null(precomp_dir) || is.null(blocks_file) || is.null(chr)) {
    stop("Required: --precomp_dir, --blocks, --chr\n",
         "Usage: Rscript STEP_C01n_blockwise_peeling_diagnostic.R \\\n",
         "  --precomp_dir precomp/ --blocks blocks.tsv --chr C_gar_LG01 \\\n",
         "  [--pruned_list pruned_samples.txt OR per_chr_pruning.tsv] \\\n",
         "  [--coseg_dir multi_inversion/] [--outdir peel_diagnostic/] \\\n",
         "  [--local_k 5] [--local_kin_thresh 0.7] [--n_pcs 2]")
  }

  precomp_file <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  if (!file.exists(precomp_file)) stop("Missing: ", precomp_file)
  message("[PEEL] Loading: ", precomp_file)
  precomp <- readRDS(precomp_file)

  blocks <- fread(blocks_file)
  for (col in c("chr", "chromosome", "chrom")) {
    if (col %in% names(blocks)) { blocks <- blocks[get(col) == chr]; break }
  }
  if (nrow(blocks) == 0) { message("[PEEL] No blocks for ", chr); quit(save = "no") }

  out <- run_peel_diagnostic(
    precomp, blocks,
    pruned_list = pruned_list,
    coseg_dir = coseg_dir,
    chr_name = chr,
    kin_threshold = kin_threshold,
    local_kin_thresh = local_kin_thresh,
    max_peel_frac = max_peel_frac,
    local_k = local_k,
    n_pcs_use = n_pcs
  )

  fwrite(out$diagnostics, file.path(outdir, paste0("peel_diagnostic_", chr, ".tsv")), sep = "\t")
  fwrite(out$sample_groups, file.path(outdir, paste0("peel_sample_groups_", chr, ".tsv")), sep = "\t")

  print_interpretation(out$diagnostics)
  message("\n[DONE] -> ", outdir)
}
