#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01i_multi_inversion_decomposition.R  (v8.4)
#
# MULTI-INVERSION DECOMPOSITION via marker co-segregation.
#
# Problem: a candidate region may contain 2+ overlapping inversions.
# PCA only finds the dominant one. Marker peeling is not sample-aware.
#
# Solution: markers that belong to the SAME inversion co-segregate across
# samples. If marker X and marker Y are both inside inversion A, then
# every sample that is HET at X is also HET at Y. But markers from
# different inversions segregate INDEPENDENTLY.
#
# Algorithm:
#   1. Extract genotypes from Clair3 VCF for the candidate region
#   2. Discretize to 012 vectors (per sample, per marker)
#   3. Compute marker x marker correlation matrix (across samples)
#   4. Cluster markers into co-segregating BLOCKS (each block = one system)
#   5. For each block: classify samples using the block's markers
#   6. Cross-reference blocks: overlapping, nested, independent?
#   7. Detect rare inversions: markers co-segregating in <5% of samples
#
# Key insight: this does NOT use PCA at all. It works directly on discrete
# genotype vectors. Hamming distance for sample comparison. Correlation
# for marker co-segregation. No averaging. No collapsing.
#
# Uses the vector engine primitives from STEP20 (encode_string_012,
# hamming_012, consensus_012).
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>
#   --vcf_dir <dir>           -- Clair3 per-chr VCFs
#   --samples <sample_list>
#   --outdir <dir>
#
# Outputs:
#   marker_coseg_blocks.tsv.gz       -- per-block: markers, span, system_id
#   marker_coseg_samples.tsv.gz      -- per-sample per-block: 012 vector, class
#   multi_inv_systems.tsv.gz         -- per-system: span, frequency, n_markers
#   multi_inv_overlaps.tsv           -- pairwise system relationships
#   multi_inv_summary.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; vcf_dir <- NULL; samples_file <- NULL
outdir <- "multi_inversion"; tier_max <- 3L
min_block_markers <- 20L; cor_threshold <- 0.5
pop_mode <- "hatchery"

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))          { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--vcf_dir" && i < length(args))    { vcf_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--tier_max" && i < length(args))   { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--min_block" && i < length(args))  { min_block_markers <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--cor_thresh" && i < length(args)) { cor_threshold <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--mode" && i < length(args))       { pop_mode <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

# Mode-dependent extension parameters
# Hatchery: high background LD from small Ne, so require HIGHER correlation
# to extend (otherwise extends into family LD).
# Wild: lower background LD, can extend with lower correlation.
if (pop_mode == "hatchery") {
  EXT_COR_DROP <- 0.5       # stop extending when cor drops below this
  EXT_MAX_SCAN <- 300L      # don't scan too far (family LD extends far)
  FRAG_COR_MIN <- 0.4       # per-sample fragment: stricter
  FRAG_MAX_SCAN <- 500L
  message("[C01i] Mode: HATCHERY (cor_drop=0.5, conservative extension)")
} else {
  EXT_COR_DROP <- 0.3
  EXT_MAX_SCAN <- 500L
  FRAG_COR_MIN <- 0.2
  FRAG_MAX_SCAN <- 1000L
  message("[C01i] Mode: WILD (cor_drop=0.3, standard extension)")
}

if (is.null(scores_file) || is.null(vcf_dir)) stop("--scores and --vcf_dir required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 300
cand_dt <- fread(scores_file)[tier <= tier_max][order(tier, -final_score)]
message("[C01i] Candidates: ", nrow(cand_dt))

# =============================================================================
# VECTOR ENGINE PRIMITIVES (from STEP20)
# =============================================================================

discretize_gt <- function(gt_str) {
  # VCF GT string -> 0/1/2
  fifelse(gt_str %in% c("0/0", "0|0"), 0L,
  fifelse(gt_str %in% c("0/1", "1/0", "0|1", "1|0"), 1L,
  fifelse(gt_str %in% c("1/1", "1|1"), 2L, NA_integer_)))
}

encode_012 <- function(x) paste(ifelse(is.na(x), "9", as.character(x)), collapse = "")

hamming_012 <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  nv <- sum(ok)
  if (nv == 0) return(NA_real_)
  sum(a[ok] != b[ok]) / nv
}

consensus_012 <- function(mat) {
  apply(mat, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA_integer_)
    as.integer(which.max(tabulate(row + 1L, nbins = 3L)) - 1L)
  })
}

# =============================================================================
# EXTRACT GENOTYPE MATRIX
# =============================================================================

extract_gt_matrix <- function(vcf_file, chr, start_bp, end_bp) {
  if (!file.exists(vcf_file)) return(NULL)
  region <- paste0(chr, ":", as.integer(start_bp), "-", as.integer(end_bp))
  cmd <- paste0("bcftools query -r ", region,
                " -f '%POS\\t[%GT\\t]\\n' ", vcf_file, " 2>/dev/null")
  raw <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"),
                   error = function(e) NULL)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  cmd_h <- paste0("bcftools query -l ", vcf_file, " 2>/dev/null")
  samp <- tryCatch(fread(cmd = cmd_h, header = FALSE)[[1]], error = function(e) NULL)
  if (is.null(samp)) return(NULL)

  ns <- length(samp)
  if (ncol(raw) < ns + 1) return(NULL)
  pos <- raw[[1]]

  gt <- matrix(NA_integer_, nrow(raw), ns)
  colnames(gt) <- samp
  for (j in seq_len(ns)) {
    gt[, j] <- discretize_gt(as.character(raw[[j + 1]]))
  }

  # Filter: MAF >= 0.02, missing < 40%
  miss <- rowMeans(is.na(gt))
  af <- rowMeans(gt, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- miss < 0.4 & maf >= 0.02
  if (sum(keep) < 50) return(NULL)

  list(pos = pos[keep], gt = gt[keep, ], samples = samp)
}

# =============================================================================
# MARKER CO-SEGREGATION CLUSTERING
# =============================================================================

find_coseg_blocks <- function(gt, pos, min_block, cor_thresh) {
  nm <- nrow(gt); ns <- ncol(gt)
  if (nm < min_block * 2) return(list())

  # Subsample markers for correlation if too many (>5000)
  if (nm > 5000) {
    sub_idx <- sort(sample(nm, 5000))
    gt_sub <- gt[sub_idx, ]; pos_sub <- pos[sub_idx]
  } else {
    gt_sub <- gt; pos_sub <- pos; sub_idx <- seq_len(nm)
  }
  nm_sub <- nrow(gt_sub)

  # Compute marker x marker correlation (across samples)
  # Use absolute correlation: markers in same inversion correlate +1 or -1
  # (depending on which allele is REF vs ALT)
  message("    Computing marker correlation (", nm_sub, " x ", nm_sub, ")...")

  # Impute NA with row median for correlation
  gt_imp <- gt_sub
  for (mi in seq_len(nm_sub)) {
    na_idx <- is.na(gt_imp[mi, ])
    if (any(na_idx)) gt_imp[mi, na_idx] <- median(gt_imp[mi, !na_idx], na.rm = TRUE)
  }

  cor_mat <- cor(t(gt_imp), use = "pairwise.complete.obs")
  cor_mat[!is.finite(cor_mat)] <- 0
  abs_cor <- abs(cor_mat)

  # Threshold: build adjacency from high absolute correlation
  adj <- abs_cor >= cor_thresh
  diag(adj) <- FALSE

  # Connected components (BFS)
  visited <- rep(FALSE, nm_sub)
  block_id <- rep(0L, nm_sub)
  bid <- 0L

  for (mi in seq_len(nm_sub)) {
    if (visited[mi] || !any(adj[mi, ])) next
    bid <- bid + 1L
    queue <- mi
    while (length(queue) > 0) {
      curr <- queue[1]; queue <- queue[-1]
      if (visited[curr]) next
      visited[curr] <- TRUE
      block_id[curr] <- bid
      neighbors <- which(adj[curr, ] & !visited)
      queue <- c(queue, neighbors)
    }
  }

  # Filter small blocks
  block_tab <- table(block_id[block_id > 0])
  valid <- as.integer(names(block_tab[block_tab >= min_block]))
  if (length(valid) == 0) return(list())

  # Renumber
  new_ids <- setNames(seq_along(valid), as.character(valid))
  block_id[block_id > 0 & block_id %in% valid] <-
    new_ids[as.character(block_id[block_id > 0 & block_id %in% valid])]
  block_id[!block_id %in% seq_along(valid)] <- 0L

  # Map back to full marker set if subsampled
  full_block_id <- rep(0L, nm)
  full_block_id[sub_idx] <- block_id

  # Extend to non-subsampled markers: assign each unsampled marker
  # to the block of its nearest subsampled neighbor (by correlation)
  if (nm > nm_sub) {
    non_sub <- setdiff(seq_len(nm), sub_idx)
    for (mi in non_sub) {
      # Correlate this marker with block consensus profiles
      m_vals <- gt[mi, ]
      m_imp <- m_vals; m_imp[is.na(m_imp)] <- median(m_imp, na.rm = TRUE)
      best_cor <- 0; best_block <- 0L
      for (bi in seq_along(valid)) {
        block_markers <- sub_idx[block_id == bi]
        if (length(block_markers) == 0) next
        block_consensus <- consensus_012(gt[block_markers, , drop = FALSE])
        r <- abs(cor(m_imp, block_consensus, use = "complete.obs"))
        if (is.finite(r) && r > best_cor) { best_cor <- r; best_block <- bi }
      }
      if (best_cor >= cor_thresh * 0.8) full_block_id[mi] <- best_block
    }
  }

  # Build block info
  blocks <- list()
  for (bi in seq_along(valid)) {
    b_idx <- which(full_block_id == bi)
    if (length(b_idx) < min_block) next
    b_pos <- pos[b_idx]
    blocks[[length(blocks) + 1]] <- list(
      block_id = bi, marker_idx = b_idx, marker_pos = b_pos,
      n_markers = length(b_idx),
      span_start = min(b_pos), span_end = max(b_pos),
      span_mb = round((max(b_pos) - min(b_pos)) / 1e6, 3),
      mean_internal_cor = round(mean(abs_cor[block_id == bi, block_id == bi], na.rm = TRUE), 3)
    )
  }
  blocks
}

# =============================================================================
# CLASSIFY SAMPLES PER BLOCK
# =============================================================================

classify_samples_per_block <- function(gt, block, samples) {
  b_gt <- gt[block$marker_idx, , drop = FALSE]
  ns <- ncol(b_gt)

  # Per-sample: count 0, 1, 2 at block markers
  n0 <- colSums(b_gt == 0, na.rm = TRUE)
  n1 <- colSums(b_gt == 1, na.rm = TRUE)
  n2 <- colSums(b_gt == 2, na.rm = TRUE)
  n_valid <- n0 + n1 + n2

  # Het rate per sample across block markers
  het_rate <- n1 / pmax(n_valid, 1)

  # Dosage mean (0-2 scale)
  dosage_mean <- (n1 + 2 * n2) / pmax(n_valid, 1)

  # Classify: k=3 on dosage_mean
  km <- tryCatch(kmeans(dosage_mean, centers = 3, nstart = 10), error = function(e) NULL)
  if (!is.null(km)) {
    co <- order(km$centers[, 1])
    group <- character(ns)
    group[km$cluster == co[1]] <- "HOMO_REF"
    group[km$cluster == co[2]] <- "HET"
    group[km$cluster == co[3]] <- "HOMO_INV"
  } else {
    # Fallback: threshold-based
    group <- fifelse(dosage_mean < 0.5, "HOMO_REF",
             fifelse(dosage_mean > 1.5, "HOMO_INV", "HET"))
  }

  # Check for rare: if smallest group has < 5 samples, check if it's real
  sizes <- table(group)
  if (min(sizes) <= 2 && max(sizes) > ns * 0.8) {
    # Possibly rare: relabel using k=2
    km2 <- tryCatch(kmeans(dosage_mean, centers = 2, nstart = 10), error = function(e) NULL)
    if (!is.null(km2)) {
      co2 <- order(km2$centers[, 1])
      minor_idx <- which(km2$cluster == co2[2])
      major_idx <- which(km2$cluster == co2[1])
      if (length(minor_idx) <= 5) {
        # Rare: check for het bridge
        group <- rep("HOMO_REF", ns)
        group[minor_idx] <- "RARE_INV"
        # Find intermediates
        mid_dosage <- (km2$centers[co2[1], 1] + km2$centers[co2[2], 1]) / 2
        het_candidates <- which(abs(dosage_mean - mid_dosage) < 0.3 &
                                 !seq_len(ns) %in% minor_idx)
        if (length(het_candidates) >= 2) group[het_candidates] <- "RARE_HET"
      }
    }
  }

  # Profile string per sample
  profiles <- character(ns)
  for (si in seq_len(ns)) profiles[si] <- encode_012(b_gt[, si])

  freq <- (sum(group %in% c("HOMO_INV", "RARE_INV")) * 2 +
            sum(group %in% c("HET", "RARE_HET"))) / (2 * ns)

  data.table(
    sample = samples, group = group,
    dosage_mean = round(dosage_mean, 3), het_rate = round(het_rate, 3),
    n_valid = n_valid, profile = profiles,
    frequency = round(freq, 4)
  )
}

# =============================================================================
# CORE PACKAGE EXTRACTION
# =============================================================================
# From a co-segregation block + sample classifications, extract the "core"
# markers: sites where HOMO_REF and HOMO_INV are perfectly (or near-perfectly)
# separated. These are the inversion's diagnostic SNP package.

extract_core_package <- function(gt, block, sample_class) {
  b_gt <- gt[block$marker_idx, , drop = FALSE]
  ref_idx <- match(sample_class[group == "HOMO_REF"]$sample, colnames(gt))
  inv_idx <- match(sample_class[group %in% c("HOMO_INV", "RARE_INV")]$sample, colnames(gt))
  ref_idx <- ref_idx[!is.na(ref_idx)]; inv_idx <- inv_idx[!is.na(inv_idx)]

  if (length(ref_idx) < 3 || length(inv_idx) < 2) return(NULL)

  core_rows <- list()
  for (mi in seq_len(nrow(b_gt))) {
    g_ref <- b_gt[mi, ref_idx]; g_inv <- b_gt[mi, inv_idx]
    g_ref <- g_ref[!is.na(g_ref)]; g_inv <- g_inv[!is.na(g_inv)]
    if (length(g_ref) < 3 || length(g_inv) < 2) next

    freq0_ref <- mean(g_ref == 0); freq0_inv <- mean(g_inv == 0)
    separation <- abs(freq0_ref - freq0_inv)

    # Core: near-perfect separation (>0.85 difference)
    # Shell: moderate separation (0.5-0.85)
    # Noise: weak (<0.5)
    tier <- if (separation > 0.85) "core"
            else if (separation > 0.5) "shell"
            else "noise"

    if (tier != "noise") {
      ref_allele <- if (freq0_ref > 0.5) 0L else 2L
      core_rows[[length(core_rows) + 1]] <- data.table(
        marker_idx = block$marker_idx[mi],
        pos = block$marker_pos[mi],
        tier = tier, separation = round(separation, 3),
        ref_allele = ref_allele, inv_allele = 2L - ref_allele
      )
    }
  }
  if (length(core_rows) > 0) rbindlist(core_rows) else NULL
}

# =============================================================================
# SEED EXTENSION USING CORE PACKAGE
# =============================================================================
# Given the core package (diagnostic markers inside the triangle), extend
# the inversion boundaries by checking if markers OUTSIDE the triangle
# also correlate with the core.
#
# This is like Snake 1 seed extension but on discrete genotypes:
# grow left/right from the block boundary, checking each marker for
# correlation with the core consensus.

extend_block_boundaries <- function(gt_full, pos_full, core_package,
                                     sample_class, block,
                                     max_extend = EXT_MAX_SCAN,
                                     cor_drop = EXT_COR_DROP) {
  if (is.null(core_package) || nrow(core_package) < 10) return(block)

  # Build core consensus: for each sample, dosage at core markers
  core_idx <- core_package[tier == "core"]$marker_idx
  if (length(core_idx) < 5) core_idx <- core_package$marker_idx[seq_len(min(20, nrow(core_package)))]

  core_gt <- gt_full[core_idx, , drop = FALSE]
  core_consensus <- colMeans(core_gt, na.rm = TRUE)  # per-sample mean dosage at core

  n_total <- nrow(gt_full)
  orig_left <- min(block$marker_idx); orig_right <- max(block$marker_idx)

  # Extend LEFT
  new_left <- orig_left
  for (mi in seq(orig_left - 1, max(1, orig_left - max_extend), by = -1)) {
    m_gt <- gt_full[mi, ]
    m_valid <- !is.na(m_gt) & is.finite(core_consensus)
    if (sum(m_valid) < 20) next
    r <- abs(cor(m_gt[m_valid], core_consensus[m_valid]))
    if (!is.finite(r) || r < cor_drop) break
    new_left <- mi
  }

  # Extend RIGHT
  new_right <- orig_right
  for (mi in seq(orig_right + 1, min(n_total, orig_right + max_extend))) {
    m_gt <- gt_full[mi, ]
    m_valid <- !is.na(m_gt) & is.finite(core_consensus)
    if (sum(m_valid) < 20) next
    r <- abs(cor(m_gt[m_valid], core_consensus[m_valid]))
    if (!is.finite(r) || r < cor_drop) break
    new_right <- mi
  }

  extended <- new_left < orig_left || new_right > orig_right
  if (extended) {
    new_idx <- new_left:new_right
    new_pos <- pos_full[new_idx]
    block$marker_idx <- new_idx
    block$marker_pos <- new_pos
    block$n_markers <- length(new_idx)
    block$span_start <- min(new_pos)
    block$span_end <- max(new_pos)
    block$span_mb <- round((max(new_pos) - min(new_pos)) / 1e6, 3)
    block$extended_left <- orig_left - new_left
    block$extended_right <- new_right - orig_right
  } else {
    block$extended_left <- 0L
    block$extended_right <- 0L
  }
  block
}

# =============================================================================
# EXPORT FOR C01h RECOMBINANT SCANNER
# =============================================================================
# Write per-system informative marker files that C01h can consume directly
# instead of computing its own informative markers from the crude 3-band.

export_for_recombinant_scanner <- function(core_package, block, sample_class,
                                            chr, iid, system_id, outdir) {
  if (is.null(core_package) || nrow(core_package) == 0) return(invisible(NULL))

  recomb_dir <- file.path(outdir, "for_C01h")
  dir.create(recomb_dir, recursive = TRUE, showWarnings = FALSE)

  # Write informative markers
  markers_out <- core_package[, .(pos, tier, separation, ref_allele, inv_allele)]
  markers_out[, `:=`(chrom = chr, interval_id = iid, system_id = system_id)]
  fwrite(markers_out,
         file.path(recomb_dir, paste0(chr, "_I", iid, "_S", system_id, "_markers.tsv.gz")),
         sep = "\t")

  # Write sample classifications
  class_out <- sample_class[, .(sample, group, dosage_mean, het_rate)]
  class_out[, `:=`(chrom = chr, interval_id = iid, system_id = system_id)]
  fwrite(class_out,
         file.path(recomb_dir, paste0(chr, "_I", iid, "_S", system_id, "_samples.tsv.gz")),
         sep = "\t")

  # Write het sample list for C01h to scan
  het_samples <- sample_class[group %in% c("HET", "RARE_HET")]$sample
  if (length(het_samples) > 0) {
    writeLines(het_samples,
               file.path(recomb_dir, paste0(chr, "_I", iid, "_S", system_id, "_het_list.txt")))
  }

  invisible(NULL)
}

# =============================================================================
# DEEP EXTENSION: PER-SAMPLE ANCESTRAL FRAGMENT BOUNDARIES
# =============================================================================
# After basic seed extension (block-level), scan further out per-sample.
# For each carrier (HET or HOMO_INV), check how far their genotype
# correlation with the core consensus persists. Samples whose correlation
# drops earlier had an ancestor that recombined closer to the breakpoint.
# The per-sample extension distance = ancestral fragment length.

compute_ancestral_fragments <- function(gt_full, pos_full, core_package,
                                         sample_class, block,
                                         max_scan = FRAG_MAX_SCAN,
                                         cor_min = FRAG_COR_MIN) {
  if (is.null(core_package) || nrow(core_package) < 10) return(data.table())

  core_idx <- core_package[tier == "core"]$marker_idx
  if (length(core_idx) < 5) core_idx <- core_package$marker_idx[seq_len(min(20, nrow(core_package)))]

  # Core consensus: per-sample mean dosage at core markers
  core_gt <- gt_full[core_idx, , drop = FALSE]
  core_consensus <- colMeans(core_gt, na.rm = TRUE)

  # Only scan carriers (HET + INV)
  carriers <- sample_class[group %in% c("HET", "HOMO_INV", "RARE_HET", "RARE_INV")]$sample
  carrier_idx <- match(carriers, colnames(gt_full))
  carrier_idx <- carrier_idx[!is.na(carrier_idx)]
  if (length(carrier_idx) < 2) return(data.table())

  n_total <- nrow(gt_full)
  orig_left <- min(block$marker_idx); orig_right <- max(block$marker_idx)

  frag_rows <- list()
  for (ci in carrier_idx) {
    samp <- colnames(gt_full)[ci]
    samp_consensus <- core_gt[, ci]

    # Scan LEFT: how far does this sample's correlation persist?
    left_boundary <- orig_left
    for (mi in seq(orig_left - 1, max(1, orig_left - max_scan), by = -1)) {
      m_gt_all <- gt_full[mi, ]
      m_valid <- !is.na(m_gt_all) & is.finite(core_consensus)
      if (sum(m_valid) < 10) next
      # Per-sample: does THIS sample's genotype at this marker correlate with core?
      this_gt <- gt_full[mi, ci]
      if (is.na(this_gt)) next
      # Compare to core consensus for this sample
      r <- abs(cor(gt_full[mi, m_valid], core_consensus[m_valid]))
      if (!is.finite(r) || r < cor_min) break
      left_boundary <- mi
    }

    # Scan RIGHT
    right_boundary <- orig_right
    for (mi in seq(orig_right + 1, min(n_total, orig_right + max_scan))) {
      m_gt_all <- gt_full[mi, ]
      m_valid <- !is.na(m_gt_all) & is.finite(core_consensus)
      if (sum(m_valid) < 10) next
      this_gt <- gt_full[mi, ci]
      if (is.na(this_gt)) next
      r <- abs(cor(gt_full[mi, m_valid], core_consensus[m_valid]))
      if (!is.finite(r) || r < cor_min) break
      right_boundary <- mi
    }

    frag_start <- pos_full[left_boundary]
    frag_end <- pos_full[right_boundary]
    frag_bp <- frag_end - frag_start
    ext_left <- orig_left - left_boundary
    ext_right <- right_boundary - orig_right

    frag_rows[[length(frag_rows) + 1]] <- data.table(
      sample = samp,
      group = sample_class[sample == samp]$group[1],
      frag_start_bp = frag_start,
      frag_end_bp = frag_end,
      frag_length_bp = frag_bp,
      frag_length_mb = round(frag_bp / 1e6, 3),
      extended_left = ext_left,
      extended_right = ext_right,
      block_start_bp = pos_full[orig_left],
      block_end_bp = pos_full[orig_right]
    )
  }

  if (length(frag_rows) > 0) rbindlist(frag_rows) else data.table()
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_blocks <- list(); all_samples <- list()
all_overlaps <- list(); all_summary <- list()
all_cores <- list(); all_fragments <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]; chr <- cand$chrom; iid <- cand$interval_id
  message("\n[C01i] ", ci, "/", nrow(cand_dt), ": ", chr, " I", iid,
          " (", cand$start_mb, "-", cand$end_mb, " Mb)")

  vcf_file <- file.path(vcf_dir, paste0(chr, ".vcf.gz"))
  if (!file.exists(vcf_file)) {
    vcf_file <- list.files(vcf_dir, pattern = paste0(chr, ".*\\.vcf\\.gz$"),
                            full.names = TRUE)[1]
    if (is.na(vcf_file)) { message("  No VCF"); next }
  }

  gt_data <- extract_gt_matrix(vcf_file, chr, cand$start_mb * 1e6, cand$end_mb * 1e6)
  if (is.null(gt_data)) { message("  Extraction failed"); next }
  message("  Markers: ", nrow(gt_data$gt), " x ", ncol(gt_data$gt), " samples")

  # Find co-segregating blocks
  blocks <- find_coseg_blocks(gt_data$gt, gt_data$pos, min_block_markers, cor_threshold)
  if (length(blocks) == 0) { message("  No co-segregating blocks found"); next }
  message("  Blocks: ", length(blocks))

  # Classify samples per block + core extraction + extension + C01h export
  for (bi in seq_along(blocks)) {
    blk <- blocks[[bi]]
    samp_class <- classify_samples_per_block(gt_data$gt, blk, gt_data$samples)
    samp_class[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]

    # Extract core package (diagnostic markers)
    core_pkg <- extract_core_package(gt_data$gt, blk, samp_class)
    n_core <- if (!is.null(core_pkg)) sum(core_pkg$tier == "core") else 0
    n_shell <- if (!is.null(core_pkg)) sum(core_pkg$tier == "shell") else 0

    # Seed extension using core package
    blk_ext <- extend_block_boundaries(gt_data$gt, gt_data$pos, core_pkg,
                                        samp_class, blk)
    # If extended, reclassify samples with expanded marker set
    if (blk_ext$extended_left > 0 || blk_ext$extended_right > 0) {
      samp_class <- classify_samples_per_block(gt_data$gt, blk_ext, gt_data$samples)
      samp_class[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]
      # Re-extract core with extended markers
      core_pkg <- extract_core_package(gt_data$gt, blk_ext, samp_class)
      n_core <- if (!is.null(core_pkg)) sum(core_pkg$tier == "core") else 0
      n_shell <- if (!is.null(core_pkg)) sum(core_pkg$tier == "shell") else 0
    }

    # Export for C01h recombinant scanner
    export_for_recombinant_scanner(core_pkg, blk_ext, samp_class,
                                    chr, iid, bi, outdir)

    # Deep extension: per-sample ancestral fragment boundaries
    frags <- compute_ancestral_fragments(gt_data$gt, gt_data$pos, core_pkg,
                                          samp_class, blk_ext)
    if (nrow(frags) > 0) {
      frags[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]
      all_fragments[[length(all_fragments) + 1]] <- frags
      frag_range <- range(frags$frag_length_mb)
      message("    Ancestral fragments: ", nrow(frags), " carriers, ",
              "range ", frag_range[1], "-", frag_range[2], " Mb")
    }

    n_ref <- sum(samp_class$group == "HOMO_REF")
    n_het <- sum(samp_class$group %in% c("HET", "RARE_HET"))
    n_inv <- sum(samp_class$group %in% c("HOMO_INV", "RARE_INV"))
    is_rare <- n_inv <= 5

    message("    Block ", bi, ": ", blk_ext$n_markers, " markers (",
            n_core, " core + ", n_shell, " shell), ",
            blk_ext$span_mb, " Mb, freq=", samp_class$frequency[1],
            " (REF=", n_ref, " HET=", n_het, " INV=", n_inv,
            if (is_rare) " RARE" else "",
            if (blk_ext$extended_left > 0 || blk_ext$extended_right > 0)
              paste0(" ext:L+", blk_ext$extended_left, "/R+", blk_ext$extended_right)
            else "", ")")

    # Store core package
    if (!is.null(core_pkg)) {
      core_pkg[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]
      all_cores[[length(all_cores) + 1]] <- core_pkg
    }

    all_blocks[[length(all_blocks) + 1]] <- data.table(
      chrom = chr, interval_id = iid, system_id = bi,
      n_markers = blk_ext$n_markers, n_core = n_core, n_shell = n_shell,
      span_start = blk_ext$span_start, span_end = blk_ext$span_end,
      span_mb = blk_ext$span_mb, mean_cor = blk$mean_internal_cor,
      frequency = samp_class$frequency[1],
      n_homo_ref = n_ref, n_het = n_het, n_homo_inv = n_inv,
      is_rare = is_rare,
      extended_left = blk_ext$extended_left %||% 0L,
      extended_right = blk_ext$extended_right %||% 0L
    )
    all_samples[[length(all_samples) + 1]] <- samp_class
    blocks[[bi]] <- blk_ext  # update for overlap calculation
  }

  # =====================================================================
  # MARKER COMPOSITION HEATMAP (pangenome-style)
  # =====================================================================
  # Rows = markers sorted by tier (core → shell → non-core)
  # Columns = samples grouped by karyotype/regime
  # Color = allele dosage (0=REF blue, 1=HET light, 2=INV dark)
  # Core markers show clean blocks (fixed within each group).
  # Shell markers show some within-group variation.
  # Non-core markers are random — no group structure.
  # Shared rows (same color across groups) = not inversion-diagnostic.

  if (length(all_cores) > 0 && nrow(gt_data$gt) > 0) {
    # Collect all core/shell markers for this candidate
    cand_cores <- rbindlist(all_cores[vapply(all_cores, function(x)
      x$chrom[1] == chr & x$interval_id[1] == iid, logical(1))], fill = TRUE)

    if (nrow(cand_cores) > 0) {
      # Get sample classifications for ordering columns
      cand_samps <- rbindlist(all_samples[vapply(all_samples, function(x)
        x$chrom[1] == chr & x$interval_id[1] == iid, logical(1))], fill = TRUE)

      # For multi-system: combine system + group for column ordering
      if ("system_id" %in% names(cand_samps)) {
        cand_samps[, regime_label := paste0("S", system_id, "_", group)]
      } else {
        cand_samps[, regime_label := group]
      }

      # Column order: group by regime, within regime sort by dosage
      regime_order <- cand_samps[, .(med_dos = median(dosage_mean, na.rm = TRUE)),
                                  by = regime_label][order(med_dos)]$regime_label
      cand_samps[, regime_label := factor(regime_label, levels = regime_order)]
      samp_order <- cand_samps[order(regime_label, dosage_mean)]$sample

      # Row order: core first, then shell, then sample of non-core
      core_idx <- cand_cores[tier == "core"]$marker_idx
      shell_idx <- cand_cores[tier == "shell"]$marker_idx
      all_block_idx <- unique(unlist(lapply(blocks, function(b) b$marker_idx)))
      non_core_idx <- setdiff(seq_len(nrow(gt_data$gt)), c(core_idx, shell_idx))
      # Subsample non-core to keep plot manageable
      if (length(non_core_idx) > length(core_idx) + length(shell_idx)) {
        non_core_idx <- sort(sample(non_core_idx,
                                     min(length(non_core_idx),
                                         length(core_idx) + length(shell_idx))))
      }

      row_idx <- c(core_idx, shell_idx, non_core_idx)
      row_tier <- c(rep("core", length(core_idx)),
                     rep("shell", length(shell_idx)),
                     rep("non-core", length(non_core_idx)))

      # Extract genotype submatrix
      col_idx <- match(samp_order, gt_data$samples)
      col_idx <- col_idx[!is.na(col_idx)]
      if (length(col_idx) >= 10 && length(row_idx) >= 10) {
        sub_gt <- gt_data$gt[row_idx, col_idx, drop = FALSE]

        # Build data.table for plotting
        hm_rows <- list()
        for (ri in seq_len(nrow(sub_gt))) {
          for (ci in seq_len(ncol(sub_gt))) {
            hm_rows[[length(hm_rows) + 1]] <- data.table(
              marker = ri, sample = ci, dosage = sub_gt[ri, ci])
          }
        }
        hm_dt <- rbindlist(hm_rows)

        # Tier annotation breaks
        tier_breaks <- c(length(core_idx),
                          length(core_idx) + length(shell_idx))

        # Regime annotation breaks
        regime_sizes <- cand_samps[match(samp_order[col_idx], sample)]
        if (nrow(regime_sizes) > 0 && "regime_label" %in% names(regime_sizes)) {
          regime_breaks <- which(regime_sizes$regime_label[-nrow(regime_sizes)] !=
                                  regime_sizes$regime_label[-1])
        } else regime_breaks <- integer(0)

        pComp <- ggplot(hm_dt, aes(x = sample, y = marker, fill = dosage)) +
          geom_tile() +
          scale_fill_gradient2(low = "#2563eb", mid = "#e5e7eb", high = "#991b1b",
                                midpoint = 1, na.value = "white",
                                name = "Dosage",
                                breaks = c(0, 1, 2),
                                labels = c("0 (REF)", "1 (HET)", "2 (INV)")) +
          # Tier separator lines (horizontal)
          geom_hline(yintercept = tier_breaks + 0.5, color = "black", linewidth = 0.6) +
          # Regime separator lines (vertical)
          geom_vline(xintercept = regime_breaks + 0.5, color = "black", linewidth = 0.6) +
          # Tier labels
          annotate("text", x = -2, y = length(core_idx) / 2,
                    label = "CORE", size = 2, angle = 90, fontface = "bold") +
          annotate("text", x = -2, y = length(core_idx) + length(shell_idx) / 2,
                    label = "SHELL", size = 2, angle = 90) +
          annotate("text", x = -2, y = tier_breaks[2] + length(non_core_idx) / 2,
                    label = "other", size = 2, angle = 90, color = "grey50") +
          labs(title = paste0(chr, " I", iid, " -- Marker Composition by Regime"),
               subtitle = paste0(length(core_idx), " core + ", length(shell_idx),
                                " shell + ", length(non_core_idx), " other markers | ",
                                length(unique(cand_samps$regime_label)), " regimes"),
               x = "Samples (grouped by regime)", y = "Markers (sorted by tier)") +
          theme_minimal(base_size = 8) +
          theme(axis.text = element_blank(), axis.ticks = element_blank(),
                legend.key.size = unit(0.3, "cm"))

        tryCatch(ggsave(file.path(outdir, "plots",
                                   paste0(chr, "_I", iid, "_marker_composition.png")),
                         pComp, width = 12, height = 10, dpi = DPI),
                 error = function(e) message("  [PLOT] ", e$message))
        message("  Marker composition heatmap: ", nrow(sub_gt), " markers x ",
                ncol(sub_gt), " samples")
      }
    }
  }

  # Pairwise block overlaps
  if (length(blocks) >= 2) {
    for (b1 in seq_len(length(blocks) - 1)) {
      for (b2 in (b1 + 1):length(blocks)) {
        blk1 <- blocks[[b1]]; blk2 <- blocks[[b2]]
        ovl_start <- max(blk1$span_start, blk2$span_start)
        ovl_end <- min(blk1$span_end, blk2$span_end)
        ovl_bp <- max(0, ovl_end - ovl_start)

        # Sample carrier overlap
        s1 <- all_samples[[length(all_samples) - length(blocks) + b1]]
        s2 <- all_samples[[length(all_samples) - length(blocks) + b2]]
        carriers1 <- s1[group %in% c("HET", "HOMO_INV", "RARE_HET", "RARE_INV")]$sample
        carriers2 <- s2[group %in% c("HET", "HOMO_INV", "RARE_HET", "RARE_INV")]$sample
        shared <- length(intersect(carriers1, carriers2))
        jac <- shared / max(1, length(union(carriers1, carriers2)))

        # Marker correlation between blocks
        shared_markers <- intersect(blk1$marker_idx, blk2$marker_idx)

        all_overlaps[[length(all_overlaps) + 1]] <- data.table(
          chrom = chr, interval_id = iid,
          system_1 = b1, system_2 = b2,
          overlap_bp = ovl_bp, overlap_mb = round(ovl_bp / 1e6, 3),
          carrier_overlap = shared, carrier_jaccard = round(jac, 3),
          shared_markers = length(shared_markers),
          relationship = fifelse(jac > 0.8, "same_system",
                        fifelse(jac > 0.5, "linked",
                        fifelse(ovl_bp > 0, "overlapping", "independent")))
        )
      }
    }
  }

  all_summary[[length(all_summary) + 1]] <- data.table(
    chrom = chr, interval_id = iid,
    n_systems = length(blocks),
    n_rare = sum(vapply(blocks, function(b) {
      sc <- classify_samples_per_block(gt_data$gt, b, gt_data$samples)
      sum(sc$group %in% c("HOMO_INV", "RARE_INV")) <= 5
    }, logical(1))),
    systems_desc = paste(vapply(seq_along(blocks), function(bi) {
      b <- blocks[[bi]]
      paste0("S", bi, ":", b$span_mb, "Mb,", b$n_markers, "m")
    }, character(1)), collapse = " | ")
  )
}

# =============================================================================
# WRITE
# =============================================================================

message("\n[C01i] Writing...")
blk_dt <- if (length(all_blocks) > 0) rbindlist(all_blocks, fill = TRUE) else data.table()
smp_dt <- if (length(all_samples) > 0) rbindlist(all_samples, fill = TRUE) else data.table()
ovl_dt <- if (length(all_overlaps) > 0) rbindlist(all_overlaps) else data.table()
sum_dt <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()
cor_dt <- if (length(all_cores) > 0) rbindlist(all_cores, fill = TRUE) else data.table()

fwrite(blk_dt, file.path(outdir, "marker_coseg_blocks.tsv.gz"), sep = "\t")
fwrite(smp_dt, file.path(outdir, "marker_coseg_samples.tsv.gz"), sep = "\t")
fwrite(ovl_dt, file.path(outdir, "multi_inv_overlaps.tsv"), sep = "\t")
fwrite(sum_dt, file.path(outdir, "multi_inv_summary.tsv"), sep = "\t")
fwrite(cor_dt, file.path(outdir, "core_packages.tsv.gz"), sep = "\t")

frag_dt <- if (length(all_fragments) > 0) rbindlist(all_fragments, fill = TRUE) else data.table()
fwrite(frag_dt, file.path(outdir, "ancestral_fragments.tsv.gz"), sep = "\t")

if (nrow(frag_dt) > 0) {
  message("\n[C01i] === ANCESTRAL FRAGMENTS ===")
  message("  Total carriers with fragments: ", nrow(frag_dt))
  message("  Fragment length range: ", round(min(frag_dt$frag_length_mb), 2), " - ",
          round(max(frag_dt$frag_length_mb), 2), " Mb")
  message("  Mean fragment: ", round(mean(frag_dt$frag_length_mb), 2), " Mb")
  # Samples with shorter fragments = more ancestral recombination near breakpoints
  short_frags <- frag_dt[frag_length_mb < median(frag_dt$frag_length_mb) * 0.7]
  if (nrow(short_frags) > 0)
    message("  Short fragments (<70% median): ", nrow(short_frags), " carriers",
            " (possible breakpoint-proximal recombinants)")
}

if (nrow(sum_dt) > 0) {
  message("\n[C01i] === DECOMPOSITION ===")
  for (ri in seq_len(nrow(sum_dt)))
    message("  ", sum_dt$chrom[ri], " I", sum_dt$interval_id[ri], ": ",
            sum_dt$n_systems[ri], " systems (", sum_dt$n_rare[ri], " rare) -- ",
            sum_dt$systems_desc[ri])
}
message("\n[DONE] -> ", outdir)
