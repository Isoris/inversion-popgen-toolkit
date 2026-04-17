#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01l_local_structure_segments.R
#
# LOCAL STRUCTURE QUANTIFICATION FOR INVERSION SEGMENTS
#
# For each confirmed inversion candidate, quantifies the local ancestry/
# structure intensity across 5 spatial segments:
#   1. Left flank (outside inversion, upstream)
#   2. Inversion left half
#   3. Inversion core (central third)
#   4. Inversion right half
#   5. Right flank (outside inversion, downstream)
#
# Uses dosage data (NOT global NGSadmix Q) to compute local membership
# proportions via per-segment k-means/PCA, then derives:
#   - Delta_12: dominance of leading local background (maxP - secondP)
#   - Local entropy: H = -sum(P_i * ln(P_i + eps))
#   - Effective number of backgrounds: ENA = exp(H)
#
# These are computed at the LOCAL level per segment, NOT from whole-genome
# admixture, because this hatchery population is a patchwork where ancestry
# changes across chromosomes and even within chromosomes.
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>    — scored candidates (Tier 1-3)
#   --coseg_dir <dir>                     — C01i decomposition output
#   --dosage_dir <dir>                    — per-chr dosage matrices
#   --precomp_dir <dir>                   — precomp RDS (for PC1 bands)
#   --outdir <dir>
#   [--flank_bp 500000]                   — flank size in bp (default 500kb)
#   [--k_local 4]                         — k for local structure clustering
#   [--tier_max 3]                        — max tier to analyze
#
# Outputs:
#   segment_sample_structure.tsv.gz       — per-sample × per-segment metrics
#   segment_summary.tsv.gz               — per-segment aggregated metrics
#   segment_by_inversion_class.tsv.gz    — metrics stratified by HOM/HET/INV
#   cell_structure_map.tsv.gz            — per-window structure metrics
#   plots/                                — segment boxplots, entropy maps
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── chat-13 wiring: registry + write_block_safe helpers ─────────────────────
.try_source <- function(path) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    tryCatch({ source(path); TRUE }, error = function(e) {
      message("[C01l] could not source ", path, ": ", conditionMessage(e))
      FALSE
    })
  } else FALSE
}
.registries_root <- Sys.getenv("REGISTRIES", "")
.pipe_base       <- Sys.getenv("BASE", ".")
.reg_loader <- file.path(ifelse(nzchar(.registries_root),
                                  .registries_root,
                                  file.path(.pipe_base,
                                             "inversion-popgen-toolkit",
                                             "registries")),
                          "api", "R", "registry_loader.R")
.helpers    <- file.path(.pipe_base, "inversion-popgen-toolkit",
                           "inversion_modules", "phase_4_postprocessing",
                           "4b_group_proposal", "lib_decompose_helpers.R")
.try_source(.reg_loader)
.try_source(.helpers)

reg <- if (exists("load_registry", mode = "function")) {
  tryCatch(load_registry(), error = function(e) {
    message("[C01l] load_registry() failed: ", conditionMessage(e),
            " — falling back to JSON sidecar writes")
    NULL
  })
} else NULL

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; coseg_dir <- NULL; dosage_dir <- NULL
precomp_dir <- NULL; outdir <- "local_structure"
FLANK_BP <- 500000L; K_LOCAL <- 4L; tier_max <- 3L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))       { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--coseg_dir" && i < length(args)) { coseg_dir <- args[i+1]; i <- i+2 }
  else if (a == "--dosage_dir" && i < length(args)) { dosage_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp_dir" && i < length(args)) { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))  { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--flank_bp" && i < length(args)) { FLANK_BP <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--k_local" && i < length(args)) { K_LOCAL <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--tier_max" && i < length(args)) { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file)) stop("--scores required")
if (is.null(dosage_dir) && is.null(precomp_dir)) stop("--dosage_dir or --precomp_dir required")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

EPS <- 1e-15

message("[C01l] Local Structure Quantification for Inversion Segments")
message("[C01l] Scores: ", scores_file)
message("[C01l] Coseg: ", coseg_dir %||% "none")
message("[C01l] Dosage: ", dosage_dir %||% "from precomp")
message("[C01l] K_local: ", K_LOCAL, " Flank: ", FLANK_BP / 1e3, " kb")

# =============================================================================
# HELPERS
# =============================================================================

# Compute local membership proportions via k-means on dosage PCA
# Returns: matrix (n_samples × k) of soft membership proportions
compute_local_memberships <- function(dosage_mat, k = K_LOCAL) {
  # dosage_mat: markers × samples
  if (nrow(dosage_mat) < 10 || ncol(dosage_mat) < 20) return(NULL)

  # Remove markers with no variance
  vars <- apply(dosage_mat, 1, var, na.rm = TRUE)
  dosage_mat <- dosage_mat[vars > 0.01 & is.finite(vars), , drop = FALSE]
  if (nrow(dosage_mat) < 5) return(NULL)

  # PCA on transposed dosage (samples as rows)
  t_dos <- t(dosage_mat)
  t_dos[is.na(t_dos)] <- 0  # impute NA as 0 for PCA
  n_pc <- min(k + 1, ncol(t_dos), nrow(t_dos) - 1)
  pca <- tryCatch(prcomp(t_dos, center = TRUE, scale. = FALSE, rank. = n_pc),
                   error = function(e) NULL)
  if (is.null(pca)) return(NULL)

  pc_scores <- pca$x[, seq_len(min(k, ncol(pca$x))), drop = FALSE]

  # k-means on PC scores
  actual_k <- min(k, nrow(pc_scores) %/% 5)  # need at least 5 per cluster
  if (actual_k < 2) actual_k <- 2
  km <- tryCatch(kmeans(pc_scores, centers = actual_k, nstart = 10),
                  error = function(e) NULL)
  if (is.null(km)) return(NULL)

  # Convert hard assignments to soft memberships via inverse distance
  n_s <- nrow(pc_scores)
  memb <- matrix(0, nrow = n_s, ncol = actual_k)
  for (ci in seq_len(actual_k)) {
    dists <- sqrt(rowSums((pc_scores - matrix(km$centers[ci, ], nrow = n_s,
                                                ncol = ncol(pc_scores), byrow = TRUE))^2))
    memb[, ci] <- 1 / (dists + EPS)
  }
  # Normalize rows to sum to 1
  row_sums <- rowSums(memb)
  memb <- memb / row_sums
  rownames(memb) <- rownames(t_dos)
  memb
}

# Compute structure metrics from a membership proportion vector
compute_metrics <- function(p_vec) {
  # p_vec: vector of proportions summing to ~1
  p_vec <- p_vec[is.finite(p_vec) & p_vec > 0]
  if (length(p_vec) < 2) return(list(delta12 = NA_real_, entropy = NA_real_, ena = NA_real_))

  p_sorted <- sort(p_vec, decreasing = TRUE)
  delta12 <- p_sorted[1] - p_sorted[2]
  entropy <- -sum(p_vec * log(p_vec + EPS))
  ena <- exp(entropy)

  list(delta12 = delta12, entropy = entropy, ena = ena)
}

# Define the 5 segments for a candidate inversion
define_segments <- function(inv_start, inv_end, flank_bp = FLANK_BP) {
  inv_span <- inv_end - inv_start
  inv_third <- inv_span / 3

  list(
    left_flank     = c(max(0, inv_start - flank_bp), inv_start),
    inv_left_half  = c(inv_start, inv_start + inv_third),
    inv_core       = c(inv_start + inv_third, inv_end - inv_third),
    inv_right_half = c(inv_end - inv_third, inv_end),
    right_flank    = c(inv_end, inv_end + flank_bp)
  )
}

# =============================================================================
# LOAD DATA
# =============================================================================

# Candidates
cand_dt <- fread(scores_file)
if ("tier" %in% names(cand_dt)) cand_dt <- cand_dt[tier <= tier_max]
message("[C01l] ", nrow(cand_dt), " candidates (tier <= ", tier_max, ")")

# C01i sample classifications (if available)
coseg_samples <- NULL
if (!is.null(coseg_dir)) {
  f_smp <- file.path(coseg_dir, "marker_coseg_samples.tsv.gz")
  if (file.exists(f_smp)) {
    coseg_samples <- fread(f_smp)
    message("[C01l] Loaded coseg sample classifications: ", nrow(coseg_samples))
  }
}

# Precomp (for PC1 dosage proxy when dosage_dir not available)
precomp_list <- list()
if (!is.null(precomp_dir)) {
  rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
  for (f in rds_files) { obj <- readRDS(f); precomp_list[[obj$chrom]] <- obj }
  message("[C01l] Loaded precomp for ", length(precomp_list), " chromosomes")
}

# =============================================================================
# MAIN LOOP: per candidate
# =============================================================================

all_seg_sample  <- list()
all_seg_summary <- list()
all_seg_class   <- list()
all_cell_map    <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]
  chr <- cand$chrom
  inv_start <- cand$start_bp
  inv_end <- cand$end_bp
  cand_id <- if ("candidate_id" %in% names(cand)) cand$candidate_id
             else if ("region_id" %in% names(cand)) cand$region_id
             else ci

  message("\n[C01l] === Candidate ", cand_id, ": ", chr, " ",
          round(inv_start / 1e6, 2), "-", round(inv_end / 1e6, 2), " Mb ===")

  # chat-13 Finding AT: scale flank_bp with candidate span. For a 300 kb
  # candidate the default 500 kb flanks reach 1.67× the inversion length
  # and risk including unrelated structure; for a ~1 Mb inversion, 500 kb
  # (50% of span) is appropriate. Clamp between 200 kb floor (need enough
  # flank to characterise the outside regime) and 500 kb ceiling (default
  # ceiling preserved for large candidates).
  span_bp <- as.integer(inv_end - inv_start)
  cand_flank_bp <- as.integer(max(200000L, min(500000L, span_bp)))
  message("  flank_bp (chat-13 AT scaling): ", cand_flank_bp / 1e3, " kb ",
          "(span = ", round(span_bp / 1e3), " kb)")

  # Define 5 segments (use per-candidate scaled flank)
  segs <- define_segments(inv_start, inv_end, flank_bp = cand_flank_bp)
  seg_names <- c("left_flank", "inv_left_half", "inv_core", "inv_right_half", "right_flank")

  # chat-13 WIRING: remember where in all_seg_summary / all_seg_class this
  # candidate's rows begin, so we can slice them out after the segment
  # loop and write the per-candidate block.
  seg_summary_start <- length(all_seg_summary) + 1L
  seg_class_start   <- length(all_seg_class)   + 1L

  # Get PC1 data from precomp as dosage proxy
  pc <- precomp_list[[chr]]
  if (is.null(pc)) { message("  [SKIP] no precomp for ", chr); next }
  dt <- pc$dt
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  if (length(pc1_cols) < 20) { message("  [SKIP] too few PC1 columns"); next }
  sample_names <- sub("^PC_1_", "", pc1_cols)

  # Get inversion genotype assignments from C01i
  inv_class <- NULL
  if (!is.null(coseg_samples)) {
    cs_chr <- coseg_samples[grepl(chr, sample, fixed = FALSE)]
    # Try matching by candidate region overlap
    # C01i labels: HOMO_REF, HET, HOMO_INV
    if (nrow(cs_chr) > 0 && "group" %in% names(cs_chr)) {
      inv_class <- cs_chr[, .(sample, group)]
      setnames(inv_class, "group", "inv_genotype")
    }
  }

  # Fallback: use PC1 k=3 bands from precomp
  if (is.null(inv_class) || nrow(inv_class) == 0) {
    # Get windows in the inversion region
    inv_idx <- which(dt$start_bp >= inv_start & dt$end_bp <= inv_end)
    if (length(inv_idx) < 5) { message("  [SKIP] too few windows in region"); next }

    pc1_mat <- as.matrix(dt[inv_idx, ..pc1_cols])
    avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
    valid <- is.finite(avg_pc1)
    if (sum(valid) < 20) next

    vals <- avg_pc1[valid]
    km <- tryCatch(kmeans(vals, centers = 3, nstart = 10), error = function(e) NULL)
    if (is.null(km)) next

    co <- order(km$centers[, 1])
    geno <- rep("HET", length(vals))
    geno[km$cluster == co[1]] <- "HOM_STD"
    geno[km$cluster == co[3]] <- "HOM_INV"

    inv_class <- data.table(sample = sub("^PC_1_", "", names(vals)),
                             inv_genotype = geno)
  }

  message("  Inversion classes: ",
          paste(inv_class[, .N, by = inv_genotype][, paste0(inv_genotype, "=", N)], collapse = " "))

  # ── Compute local structure per segment ──
  for (si in seq_along(seg_names)) {
    seg_name <- seg_names[si]
    seg_range <- segs[[si]]

    # Get windows in this segment
    seg_idx <- which(dt$start_bp >= seg_range[1] & dt$end_bp <= seg_range[2])
    if (length(seg_idx) < 3) {
      message("  ", seg_name, ": too few windows (", length(seg_idx), ")")
      next
    }

    # Build dosage-like matrix from PC1 loadings
    # (PC1 captures the dominant axis which is the inversion signal)
    pc1_mat <- as.matrix(dt[seg_idx, ..pc1_cols])

    # Compute local memberships
    memb <- compute_local_memberships(pc1_mat, k = K_LOCAL)
    if (is.null(memb)) {
      message("  ", seg_name, ": membership computation failed")
      next
    }

    # Per-sample metrics
    for (ri in seq_len(nrow(memb))) {
      samp <- rownames(memb)[ri]
      samp_clean <- sub("^PC_1_", "", samp)
      p_vec <- memb[ri, ]
      m <- compute_metrics(p_vec)

      ig <- inv_class[sample == samp_clean]$inv_genotype
      if (length(ig) == 0) ig <- "UNKNOWN"

      all_seg_sample[[length(all_seg_sample) + 1L]] <- data.table(
        candidate_id = cand_id, chrom = chr,
        inv_start = inv_start, inv_end = inv_end,
        segment = seg_name, sample_id = samp_clean,
        inv_genotype = ig,
        delta12 = round(m$delta12, 4),
        entropy = round(m$entropy, 4),
        ena = round(m$ena, 4),
        max_p = round(max(p_vec), 4),
        n_windows = length(seg_idx)
      )
    }

    # Segment-level summary
    seg_metrics <- rbindlist(lapply(seq_len(nrow(memb)), function(ri) {
      m <- compute_metrics(memb[ri, ])
      data.table(delta12 = m$delta12, entropy = m$entropy, ena = m$ena)
    }))

    all_seg_summary[[length(all_seg_summary) + 1L]] <- data.table(
      candidate_id = cand_id, chrom = chr,
      inv_start = inv_start, inv_end = inv_end,
      segment = seg_name,
      n_windows = length(seg_idx),
      n_samples = nrow(memb),
      mean_delta12 = round(mean(seg_metrics$delta12, na.rm = TRUE), 4),
      mean_entropy = round(mean(seg_metrics$entropy, na.rm = TRUE), 4),
      mean_ena = round(mean(seg_metrics$ena, na.rm = TRUE), 4),
      sd_delta12 = round(sd(seg_metrics$delta12, na.rm = TRUE), 4)
    )

    # By inversion class
    for (ig in unique(inv_class$inv_genotype)) {
      ig_samples <- inv_class[inv_genotype == ig]$sample
      ig_idx <- which(sub("^PC_1_", "", rownames(memb)) %in% ig_samples)
      if (length(ig_idx) < 3) next

      ig_metrics <- rbindlist(lapply(ig_idx, function(ri) {
        m <- compute_metrics(memb[ri, ])
        data.table(delta12 = m$delta12, entropy = m$entropy, ena = m$ena)
      }))

      all_seg_class[[length(all_seg_class) + 1L]] <- data.table(
        candidate_id = cand_id, chrom = chr, segment = seg_name,
        inv_genotype = ig, n_samples = length(ig_idx),
        mean_delta12 = round(mean(ig_metrics$delta12, na.rm = TRUE), 4),
        mean_entropy = round(mean(ig_metrics$entropy, na.rm = TRUE), 4),
        mean_ena = round(mean(ig_metrics$ena, na.rm = TRUE), 4)
      )
    }
  }

  # ── Cell-wise map: per-window structure metrics ──
  # Broader range: from left_flank start to right_flank end
  map_start <- segs[["left_flank"]][1]
  map_end <- segs[["right_flank"]][2]
  map_idx <- which(dt$start_bp >= map_start & dt$end_bp <= map_end)

  if (length(map_idx) >= 10) {
    # Sliding window: compute local structure in blocks of 20 windows
    block_size <- 20L
    step_size <- 5L

    for (bi in seq(1, length(map_idx) - block_size + 1, by = step_size)) {
      block_idx <- map_idx[bi:min(bi + block_size - 1, length(map_idx))]
      block_center <- mean(c(dt$start_bp[block_idx[1]], dt$end_bp[tail(block_idx, 1)]))

      pc1_block <- as.matrix(dt[block_idx, ..pc1_cols])
      memb <- compute_local_memberships(pc1_block, k = min(K_LOCAL, 3))
      if (is.null(memb)) next

      # Aggregate metrics across samples
      block_metrics <- rbindlist(lapply(seq_len(nrow(memb)), function(ri) {
        m <- compute_metrics(memb[ri, ])
        data.table(delta12 = m$delta12, entropy = m$entropy, ena = m$ena)
      }))

      # Determine which segment this block falls in
      seg_label <- if (block_center < inv_start) "left_flank"
                   else if (block_center < inv_start + (inv_end - inv_start) / 3) "inv_left"
                   else if (block_center < inv_end - (inv_end - inv_start) / 3) "inv_core"
                   else if (block_center < inv_end) "inv_right"
                   else "right_flank"

      # Dominant group: which cluster has most mass on average
      mean_memb <- colMeans(memb, na.rm = TRUE)
      dom_group <- which.max(mean_memb)
      dom_frac <- max(mean_memb)

      all_cell_map[[length(all_cell_map) + 1L]] <- data.table(
        candidate_id = cand_id, chrom = chr,
        pos_bp = as.integer(block_center),
        pos_mb = round(block_center / 1e6, 4),
        segment = seg_label,
        dominant_group = dom_group,
        dominance_frac = round(dom_frac, 4),
        dominance_gap = round(sort(mean_memb, decreasing = TRUE)[1] -
                               sort(mean_memb, decreasing = TRUE)[2], 4),
        mean_delta12 = round(mean(block_metrics$delta12, na.rm = TRUE), 4),
        mean_entropy = round(mean(block_metrics$entropy, na.rm = TRUE), 4),
        mean_ena = round(mean(block_metrics$ena, na.rm = TRUE), 4)
      )
    }
  }

  # chat-13 WIRING: emit per-candidate local_structure_segments block
  # per registries/schemas/structured_block_schemas/
  # local_structure_segments.schema.json. boundary_sharpness is a small
  # derived object: sharp/moderate/diffuse classifications for left, right,
  # and overall, plus the raw core/flank delta12 values that fed them.
  {
    # Slice out just this candidate's rows from the per-candidate appends
    # above. Any rows beyond seg_summary_start at this point belong to
    # cand_id because each candidate's loop iteration is serial.
    cand_seg_summary <- if (length(all_seg_summary) >= seg_summary_start) {
      rbindlist(all_seg_summary[seq(seg_summary_start, length(all_seg_summary))],
                fill = TRUE)
    } else data.table()
    cand_seg_class <- if (length(all_seg_class) >= seg_class_start) {
      rbindlist(all_seg_class[seq(seg_class_start, length(all_seg_class))],
                fill = TRUE)
    } else data.table()

    # ── Boundary sharpness: compare core delta12 to flanks. Sharp when
    # core >> flank (clear edge); diffuse when core ≈ flank. Thresholds
    # tentative — chat-14 calibration may tune.
    bs_core  <- cand_seg_summary[segment == "inv_core",      mean_delta12][1]
    bs_lflk  <- cand_seg_summary[segment == "left_flank",    mean_delta12][1]
    bs_rflk  <- cand_seg_summary[segment == "right_flank",   mean_delta12][1]
    bs_lhalf <- cand_seg_summary[segment == "inv_left_half", mean_delta12][1]
    bs_rhalf <- cand_seg_summary[segment == "inv_right_half",mean_delta12][1]

    classify_side <- function(core, flank) {
      if (is.na(core) || is.na(flank)) return(NA_character_)
      diff <- core - flank
      if      (diff > 0.30) "sharp"
      else if (diff > 0.10) "moderate"
      else                  "diffuse"
    }
    left_class  <- classify_side(bs_core, bs_lflk)
    right_class <- classify_side(bs_core, bs_rflk)
    overall_class <- if (is.na(left_class) || is.na(right_class)) NA_character_
                     else if (left_class == right_class)           left_class
                     else                                           "asymmetric"

    .as_block_list <- function(dt, pick_cols) {
      if (is.null(dt) || nrow(dt) == 0) return(list())
      present <- intersect(pick_cols, names(dt))
      lapply(seq_len(nrow(dt)), function(i) as.list(dt[i, ..present]))
    }

    lss_block <- list(
      candidate_id = as.character(cand_id),
      segment_summary = .as_block_list(cand_seg_summary,
        c("segment", "n_windows", "n_samples",
          "mean_delta12", "mean_entropy", "mean_ena", "sd_delta12")),
      by_inversion_class = .as_block_list(cand_seg_class,
        c("segment", "inv_genotype", "n_samples",
          "mean_delta12", "mean_entropy", "mean_ena")),
      boundary_sharpness = list(
        left     = left_class,
        right    = right_class,
        overall  = overall_class,
        core_delta12        = if (is.na(bs_core))  NULL else round(bs_core, 4),
        left_flank_delta12  = if (is.na(bs_lflk))  NULL else round(bs_lflk, 4),
        right_flank_delta12 = if (is.na(bs_rflk))  NULL else round(bs_rflk, 4),
        left_half_delta12   = if (is.na(bs_lhalf)) NULL else round(bs_lhalf, 4),
        right_half_delta12  = if (is.na(bs_rhalf)) NULL else round(bs_rhalf, 4)
      )
    )

    cand_outdir <- file.path(outdir, "blocks", as.character(cand_id))
    dir.create(cand_outdir, recursive = TRUE, showWarnings = FALSE)
    if (exists("write_block_safe", mode = "function")) {
      write_block_safe(
        reg             = reg,
        candidate_id    = as.character(cand_id),
        block_type      = "local_structure_segments",
        data            = lss_block,
        source_script   = "STEP_C01l_local_structure_segments.R",
        outdir_fallback = cand_outdir
      )
    }
  }

  message("  Done: ", sum(vapply(all_seg_sample, function(x)
    x$candidate_id[1] == cand_id, logical(1))), " sample-segment records")
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

seg_sample_dt <- if (length(all_seg_sample) > 0) rbindlist(all_seg_sample, fill = TRUE) else data.table()
seg_summ_dt   <- if (length(all_seg_summary) > 0) rbindlist(all_seg_summary, fill = TRUE) else data.table()
seg_class_dt  <- if (length(all_seg_class) > 0) rbindlist(all_seg_class, fill = TRUE) else data.table()
cell_map_dt   <- if (length(all_cell_map) > 0) rbindlist(all_cell_map, fill = TRUE) else data.table()

fwrite(seg_sample_dt, file.path(outdir, "segment_sample_structure.tsv.gz"), sep = "\t")
fwrite(seg_summ_dt, file.path(outdir, "segment_summary.tsv.gz"), sep = "\t")
fwrite(seg_class_dt, file.path(outdir, "segment_by_inversion_class.tsv.gz"), sep = "\t")
fwrite(cell_map_dt, file.path(outdir, "cell_structure_map.tsv.gz"), sep = "\t")

message("\n[C01l] Outputs:")
message("  ", file.path(outdir, "segment_sample_structure.tsv.gz"), " (", nrow(seg_sample_dt), " rows)")
message("  ", file.path(outdir, "segment_summary.tsv.gz"), " (", nrow(seg_summ_dt), " rows)")
message("  ", file.path(outdir, "segment_by_inversion_class.tsv.gz"), " (", nrow(seg_class_dt), " rows)")
message("  ", file.path(outdir, "cell_structure_map.tsv.gz"), " (", nrow(cell_map_dt), " rows)")

# =============================================================================
# DIAGNOSTIC FIGURES
# =============================================================================

if (nrow(seg_sample_dt) > 0) {
  message("\n[C01l] Generating plots...")

  # Factor ordering for segments
  seg_order <- c("left_flank", "inv_left_half", "inv_core", "inv_right_half", "right_flank")
  seg_labels <- c("Left\nflank", "Inv\nleft", "Inv\ncore", "Inv\nright", "Right\nflank")
  seg_sample_dt[, segment := factor(segment, levels = seg_order, labels = seg_labels)]

  seg_pal <- c("HOM_STD" = "#3b82f6", "HET" = "#22c55e", "HOM_INV" = "#ef4444",
               "HOMO_REF" = "#3b82f6", "HOMO_INV" = "#ef4444",
               "RARE_INV" = "#a855f7", "RARE_HET" = "#84cc16", "UNKNOWN" = "grey60")

  for (cid in unique(seg_sample_dt$candidate_id)) {
    cand_data <- seg_sample_dt[candidate_id == cid]
    if (nrow(cand_data) < 50) next

    chr_lab <- cand_data$chrom[1]
    start_mb <- round(cand_data$inv_start[1] / 1e6, 1)
    end_mb <- round(cand_data$inv_end[1] / 1e6, 1)

    # ── Plot L1: Delta12 boxplot by segment × inversion class ──
    p1 <- ggplot(cand_data[is.finite(delta12)],
                  aes(x = segment, y = delta12, fill = inv_genotype)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7) +
      scale_fill_manual(values = seg_pal, name = "Inv. genotype") +
      labs(title = paste0(chr_lab, " ", start_mb, "-", end_mb, " Mb: Local Dominance (\u039412)"),
           subtitle = "Higher = one local background dominates | Lower = mixed structure",
           x = NULL, y = expression(Delta[12] ~ "(max P - second P)"),
           caption = paste0("Source: k=", K_LOCAL, " clustering on PC1 loadings per segment\n",
                           "Computed locally per segment (NOT from global admixture Q)\n",
                           "Eq: \u039412 = max(P) - second(P), where P = local membership proportions")) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.caption = element_text(size = 6, color = "grey50", hjust = 0),
            legend.position = "bottom")

    # ── Plot L2: Entropy boxplot ──
    p2 <- ggplot(cand_data[is.finite(entropy)],
                  aes(x = segment, y = entropy, fill = inv_genotype)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7) +
      scale_fill_manual(values = seg_pal, name = "Inv. genotype") +
      labs(title = paste0(chr_lab, " ", start_mb, "-", end_mb, " Mb: Local Entropy (H)"),
           subtitle = "Higher = more mixed local structure | Lower = cleaner single background",
           x = NULL, y = expression("H = " ~ -sum(P[i] ~ ln(P[i]))),
           caption = paste0("Source: local membership proportions from segment-specific PCA + k-means\n",
                           "Eq: H = -\u03A3 P_i \u00B7 ln(P_i + \u03B5), \u03B5 = 1e-15")) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.caption = element_text(size = 6, color = "grey50", hjust = 0),
            legend.position = "bottom")

    # ── Plot L3: ENA boxplot ──
    p3 <- ggplot(cand_data[is.finite(ena)],
                  aes(x = segment, y = ena, fill = inv_genotype)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7) +
      scale_fill_manual(values = seg_pal, name = "Inv. genotype") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(title = paste0(chr_lab, " ", start_mb, "-", end_mb, " Mb: Effective Number of Backgrounds (ENA)"),
           subtitle = "ENA \u2248 1 = single clean background | ENA > 2 = mixed patchwork",
           x = NULL, y = "ENA = exp(H)",
           caption = paste0("Source: exp(local entropy)\n",
                           "Eq: ENA = exp(H) where H = local Shannon entropy\n",
                           "Biological expectation: inversion core should have lower ENA than flanks")) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.caption = element_text(size = 6, color = "grey50", hjust = 0),
            legend.position = "bottom")

    for (ext in c("png", "pdf")) {
      tryCatch({
        ggsave(file.path(outdir, "plots", paste0(chr_lab, "_I", cid, "_L1_delta12.", ext)),
               p1, width = 10, height = 6, dpi = 300)
        ggsave(file.path(outdir, "plots", paste0(chr_lab, "_I", cid, "_L2_entropy.", ext)),
               p2, width = 10, height = 6, dpi = 300)
        ggsave(file.path(outdir, "plots", paste0(chr_lab, "_I", cid, "_L3_ena.", ext)),
               p3, width = 10, height = 6, dpi = 300)
      }, error = function(e) message("  [PLOT] ", e$message))
    }
  }

  # ── Plot L4: Cell-wise structure map ──
  if (nrow(cell_map_dt) > 0) {
    cell_map_dt[, segment := factor(segment,
      levels = c("left_flank", "inv_left", "inv_core", "inv_right", "right_flank"))]

    for (cid in unique(cell_map_dt$candidate_id)) {
      cm <- cell_map_dt[candidate_id == cid]
      if (nrow(cm) < 5) next
      chr_lab <- cm$chrom[1]

      inv_s <- cand_dt[candidate_id == cid | region_id == cid]$start_bp[1] / 1e6
      inv_e <- cand_dt[candidate_id == cid | region_id == cid]$end_bp[1] / 1e6

      # 3-panel: dominant group + dominance gap + entropy
      p4a <- ggplot(cm, aes(x = pos_mb, y = mean_delta12)) +
        geom_rect(aes(xmin = inv_s, xmax = inv_e, ymin = -Inf, ymax = Inf),
                   fill = "lightblue", alpha = 0.1, inherit.aes = FALSE) +
        geom_line(color = "steelblue", linewidth = 0.5) +
        geom_point(aes(color = segment), size = 1) +
        labs(y = expression(Delta[12]), title = paste0(chr_lab, " I", cid, ": Local Structure Map")) +
        theme_minimal(base_size = 8)

      p4b <- ggplot(cm, aes(x = pos_mb, y = mean_entropy)) +
        geom_rect(aes(xmin = inv_s, xmax = inv_e, ymin = -Inf, ymax = Inf),
                   fill = "lightblue", alpha = 0.1, inherit.aes = FALSE) +
        geom_line(color = "darkorange", linewidth = 0.5) +
        geom_point(aes(color = segment), size = 1) +
        labs(y = "Entropy (H)") +
        theme_minimal(base_size = 8)

      p4c <- ggplot(cm, aes(x = pos_mb, y = mean_ena)) +
        geom_rect(aes(xmin = inv_s, xmax = inv_e, ymin = -Inf, ymax = Inf),
                   fill = "lightblue", alpha = 0.1, inherit.aes = FALSE) +
        geom_line(color = "red3", linewidth = 0.5) +
        geom_point(aes(color = segment), size = 1) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
        labs(x = paste0(chr_lab, " (Mb)"), y = "ENA",
             caption = "Blue shading = inversion region | Dashed = ENA=1 (single background)") +
        theme_minimal(base_size = 8) +
        theme(plot.caption = element_text(size = 6, color = "grey50", hjust = 0))

      if (requireNamespace("patchwork", quietly = TRUE)) {
        library(patchwork)
        p4 <- p4a / p4b / p4c + plot_layout(guides = "collect") &
          theme(legend.position = "bottom")

        for (ext in c("png", "pdf")) {
          tryCatch(ggsave(file.path(outdir, "plots",
                                     paste0(chr_lab, "_I", cid, "_L4_structure_map.", ext)),
                           p4, width = 14, height = 10, dpi = 300),
                   error = function(e) message("  [PLOT] ", e$message))
        }
      }
    }
  }
}

message("\n[DONE] Local structure quantification -> ", outdir)
