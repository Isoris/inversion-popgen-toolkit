#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01h_recombinant_scanner.R  (v8.4)
#
# RECOMBINANT DETECTION within inversion candidates.
#
# For each candidate inversion, uses Clair3 VCF genotypes directly
# (NOT PCA distances) to detect recombinant haplotypes.
#
# Logic:
#   1. From Clair3 VCF, extract genotypes at all variant sites in the region
#   2. Identify inversion-informative markers: sites where Homo_1 group
#      is mostly 0/0 and Homo_2 group is mostly 1/1 (or vice versa)
#   3. For each HET sample, compute rolling het-rate across informative markers
#   4. True HET: flat het-rate (~0.5) across the whole inversion
#      Recombinant: step function, het-rate drops from ~0.5 to ~0 at the
#      recombination breakpoint
#   5. Detect changepoints in the per-sample het-rate trajectory
#   6. Classify: true_het / simple_recombinant / double_crossover / mosaic
#
# Also computes:
#   - Per-sample haplotype consistency score
#   - Pairwise genotype concordance between HET samples
#   - Breakpoint position estimates for recombinants
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>  -- Tier 1/2 candidates
#   --triangles <triangle_dir>          -- composition (band assignments)
#   --vcf_dir <dir>                     -- Clair3 per-chr VCFs (<chr>.vcf.gz)
#   --samples <sample_list.tsv>
#   --outdir <dir>
#   [--min_informative 50]              -- minimum informative markers
#   [--het_window 20]                   -- rolling window size for het-rate
#
# Outputs:
#   recombinant_calls.tsv.gz            -- per-sample classification + breakpoints
#   informative_markers.tsv.gz          -- per-candidate informative sites
#   het_trajectories.tsv.gz             -- per-sample rolling het-rate
#   recombinant_summary.tsv             -- per-candidate counts
#   plots/<cand>_het_trajectories.png   -- visual per candidate
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; vcf_dir <- NULL
samples_file <- NULL; outdir <- "recombinant_scan"
min_informative <- 50L; het_window <- 20L; tier_max <- 2L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))      { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args))  { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--vcf_dir" && i < length(args))    { vcf_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--min_informative" && i < length(args)) { min_informative <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--het_window" && i < length(args)) { het_window <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--tier_max" && i < length(args))   { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file) || is.null(vcf_dir)) stop("--scores and --vcf_dir required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 300

# =============================================================================
# LOAD
# =============================================================================

message("[C01h] Loading...")
cand_dt <- fread(scores_file)[tier <= tier_max][order(tier, -final_score)]

comp_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_sample_composition.tsv.gz")
  if (file.exists(f)) comp_dt <- fread(f)
}

real_names <- NULL
if (!is.null(samples_file) && file.exists(samples_file))
  real_names <- as.character(fread(samples_file, header = FALSE)[[1]])

message("[C01h] Candidates: ", nrow(cand_dt), " | VCF dir: ", vcf_dir)

# =============================================================================
# EXTRACT GENOTYPES FROM VCF (uses bcftools)
# =============================================================================

extract_genotypes <- function(vcf_file, chr, start_bp, end_bp) {
  if (!file.exists(vcf_file)) return(NULL)

  region <- paste0(chr, ":", as.integer(start_bp), "-", as.integer(end_bp))

  # bcftools query: extract POS and all sample genotypes
  cmd <- paste0("bcftools query -r ", region,
                " -f '%POS\\t[%GT\\t]\\n' ", vcf_file, " 2>/dev/null")
  raw <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"),
                   error = function(e) NULL)
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # Get sample names from VCF header
  cmd_h <- paste0("bcftools query -l ", vcf_file, " 2>/dev/null")
  sample_names <- tryCatch(fread(cmd = cmd_h, header = FALSE)[[1]],
                            error = function(e) NULL)
  if (is.null(sample_names)) return(NULL)

  # Parse: first column = POS, rest = genotypes
  n_samp <- length(sample_names)
  if (ncol(raw) < n_samp + 1) return(NULL)

  pos <- raw[[1]]
  gt_mat <- as.matrix(raw[, 2:(n_samp + 1)])
  colnames(gt_mat) <- sample_names

  # Convert GT strings to numeric: 0/0=0, 0/1=1, 1/1=2, ./. or other=NA
  gt_num <- matrix(NA_integer_, nrow = nrow(gt_mat), ncol = ncol(gt_mat))
  colnames(gt_num) <- sample_names
  for (j in seq_len(ncol(gt_mat))) {
    g <- gt_mat[, j]
    gt_num[, j] <- fifelse(g %in% c("0/0", "0|0"), 0L,
                   fifelse(g %in% c("0/1", "1/0", "0|1", "1|0"), 1L,
                   fifelse(g %in% c("1/1", "1|1"), 2L, NA_integer_)))
  }

  list(pos = pos, gt = gt_num, samples = sample_names)
}

# =============================================================================
# FIND INFORMATIVE MARKERS
# =============================================================================

find_informative_markers <- function(gt_data, homo1_samples, homo2_samples,
                                      min_freq_diff = 0.7, max_missing = 0.3) {
  gt <- gt_data$gt; pos <- gt_data$pos; samples <- gt_data$samples

  h1_idx <- match(homo1_samples, samples); h1_idx <- h1_idx[!is.na(h1_idx)]
  h2_idx <- match(homo2_samples, samples); h2_idx <- h2_idx[!is.na(h2_idx)]

  if (length(h1_idx) < 5 || length(h2_idx) < 5) return(NULL)

  informative <- logical(length(pos))
  h1_allele <- integer(length(pos))  # expected allele for Homo_1 at informative sites

  for (si in seq_along(pos)) {
    g1 <- gt[si, h1_idx]; g2 <- gt[si, h2_idx]
    g1 <- g1[!is.na(g1)]; g2 <- g2[!is.na(g2)]
    if (length(g1) < 3 || length(g2) < 3) next

    # Missingness check
    miss1 <- 1 - length(g1) / length(h1_idx)
    miss2 <- 1 - length(g2) / length(h2_idx)
    if (miss1 > max_missing || miss2 > max_missing) next

    # Frequency of allele 0 in each group
    freq0_h1 <- mean(g1 == 0)
    freq0_h2 <- mean(g2 == 0)

    # Informative: one group mostly 0/0, other mostly 1/1
    if (freq0_h1 > (1 - min_freq_diff) && freq0_h2 < min_freq_diff) {
      informative[si] <- TRUE; h1_allele[si] <- 0L
    } else if (freq0_h1 < min_freq_diff && freq0_h2 > (1 - min_freq_diff)) {
      informative[si] <- TRUE; h1_allele[si] <- 2L
    }
  }

  if (sum(informative) < 5) return(NULL)

  data.table(
    pos = pos[informative],
    h1_expected = h1_allele[informative],
    h2_expected = 2L - h1_allele[informative],
    idx = which(informative)
  )
}

# =============================================================================
# PER-SAMPLE HET-RATE TRAJECTORY
# =============================================================================

compute_het_trajectory <- function(gt_data, markers, sample_name, window_size) {
  samples <- gt_data$samples
  sidx <- match(sample_name, samples)
  if (is.na(sidx)) return(NULL)

  # Extract genotypes at informative sites
  gts <- gt_data$gt[markers$idx, sidx]
  pos <- markers$pos
  n_m <- length(gts)
  if (n_m < window_size) return(NULL)

  # Rolling het-rate
  half_w <- window_size %/% 2
  results <- list()
  for (mi in seq(half_w + 1, n_m - half_w)) {
    w_range <- (mi - half_w):(mi + half_w)
    w_gts <- gts[w_range]
    w_valid <- w_gts[!is.na(w_gts)]
    if (length(w_valid) < 5) next

    het_rate <- mean(w_valid == 1)
    hom_ref_rate <- mean(w_valid == 0)
    hom_alt_rate <- mean(w_valid == 2)

    # Concordance with Homo_1 expected alleles
    h1_exp <- markers$h1_expected[w_range]
    h1_concordance <- mean(w_gts == h1_exp, na.rm = TRUE)
    h2_concordance <- mean(w_gts == markers$h2_expected[w_range], na.rm = TRUE)

    results[[length(results) + 1]] <- data.table(
      marker_idx = mi, pos = pos[mi],
      het_rate = round(het_rate, 3),
      hom_ref_rate = round(hom_ref_rate, 3),
      hom_alt_rate = round(hom_alt_rate, 3),
      h1_concordance = round(h1_concordance, 3),
      h2_concordance = round(h2_concordance, 3)
    )
  }
  if (length(results) > 0) rbindlist(results) else NULL
}

# =============================================================================
# DETECT CHANGEPOINTS (simple: max absolute shift in het-rate)
# =============================================================================

detect_changepoints <- function(traj, min_shift = 0.25, min_segment = 10) {
  if (is.null(traj) || nrow(traj) < min_segment * 2) return(list(
    n_changepoints = 0L, changepoint_pos = integer(0),
    class = "insufficient_data", max_shift = 0
  ))

  hr <- traj$het_rate
  n <- length(hr)

  # Binary segmentation: find position where splitting into two segments
  # maximizes the difference in mean het-rate
  best_pos <- 0; best_shift <- 0
  for (cp in min_segment:(n - min_segment)) {
    left_mean <- mean(hr[1:cp])
    right_mean <- mean(hr[(cp+1):n])
    shift <- abs(left_mean - right_mean)
    if (shift > best_shift) { best_shift <- shift; best_pos <- cp }
  }

  changepoints <- integer(0)
  if (best_shift >= min_shift) {
    changepoints <- traj$pos[best_pos]

    # Check for double-crossover: second changepoint in the larger segment
    if (best_pos > min_segment * 2) {
      hr_left <- hr[1:best_pos]
      best2 <- 0; shift2 <- 0
      for (cp2 in min_segment:(length(hr_left) - min_segment)) {
        s <- abs(mean(hr_left[1:cp2]) - mean(hr_left[(cp2+1):length(hr_left)]))
        if (s > shift2) { shift2 <- s; best2 <- cp2 }
      }
      if (shift2 >= min_shift) changepoints <- c(traj$pos[best2], changepoints)
    }
    if (n - best_pos > min_segment * 2) {
      hr_right <- hr[(best_pos+1):n]
      best3 <- 0; shift3 <- 0
      for (cp3 in min_segment:(length(hr_right) - min_segment)) {
        s <- abs(mean(hr_right[1:cp3]) - mean(hr_right[(cp3+1):length(hr_right)]))
        if (s > shift3) { shift3 <- s; best3 <- cp3 }
      }
      if (shift3 >= min_shift) changepoints <- c(changepoints, traj$pos[best_pos + best3])
    }
  }

  # Classify
  n_cp <- length(changepoints)
  cls <- if (n_cp == 0 && mean(hr) > 0.3) "true_het"
         else if (n_cp == 0 && mean(hr) <= 0.3) "true_homo"
         else if (n_cp == 1) "simple_recombinant"
         else if (n_cp == 2) "double_crossover"
         else "complex_mosaic"

  list(n_changepoints = n_cp, changepoint_pos = changepoints,
       class = cls, max_shift = round(best_shift, 3),
       mean_het = round(mean(hr), 3))
}

# =============================================================================
# MAIN
# =============================================================================

all_calls <- list()
all_markers <- list()
all_traj <- list()
all_summary <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]; chr <- cand$chrom; iid <- cand$interval_id
  message("\n[C01h] ", ci, "/", nrow(cand_dt), ": ", chr, " I", iid,
          " (", cand$start_mb, "-", cand$end_mb, " Mb)")

  # Get band assignments
  iv_comp <- comp_dt[chrom == chr & interval_id == iid]
  if (nrow(iv_comp) == 0) { message("  No composition"); next }

  # Map band names to karyotype
  homo1 <- iv_comp[band == "band1"]$sample
  het_samples <- iv_comp[band == "band2"]$sample
  homo2 <- iv_comp[band == "band3"]$sample
  message("  Groups: Homo1=", length(homo1), " Het=", length(het_samples),
          " Homo2=", length(homo2))

  if (length(homo1) < 5 || length(homo2) < 5 || length(het_samples) < 3) {
    message("  Too few samples in groups"); next
  }

  # Find VCF file
  vcf_file <- file.path(vcf_dir, paste0(chr, ".vcf.gz"))
  if (!file.exists(vcf_file)) {
    # Try other naming
    vcf_file <- list.files(vcf_dir, pattern = paste0(chr, ".*\\.vcf\\.gz$"),
                            full.names = TRUE)[1]
    if (is.na(vcf_file) || !file.exists(vcf_file)) {
      message("  No VCF for ", chr); next
    }
  }

  # Extract genotypes
  start_bp <- cand$start_mb * 1e6; end_bp <- cand$end_mb * 1e6
  gt_data <- extract_genotypes(vcf_file, chr, start_bp, end_bp)
  if (is.null(gt_data)) { message("  VCF extraction failed"); next }
  message("  Variants: ", length(gt_data$pos), " | Samples: ", length(gt_data$samples))

  # Find informative markers
  markers <- find_informative_markers(gt_data, homo1, homo2)
  if (is.null(markers) || nrow(markers) < min_informative) {
    message("  Informative markers: ", if (is.null(markers)) 0 else nrow(markers),
            " (need ", min_informative, ")")
    next
  }
  message("  Informative markers: ", nrow(markers))

  markers[, `:=`(chrom = chr, interval_id = iid)]
  all_markers[[length(all_markers) + 1]] <- markers

  # Per-HET-sample trajectory + classification
  for (samp in het_samples) {
    traj <- compute_het_trajectory(gt_data, markers, samp, het_window)
    if (is.null(traj)) next

    cp <- detect_changepoints(traj)

    traj[, `:=`(chrom = chr, interval_id = iid, sample = samp)]
    all_traj[[length(all_traj) + 1]] <- traj

    all_calls[[length(all_calls) + 1]] <- data.table(
      chrom = chr, interval_id = iid, sample = samp,
      class = cp$class,
      n_changepoints = cp$n_changepoints,
      changepoint_positions = paste(cp$changepoint_pos, collapse = ";"),
      max_shift = cp$max_shift,
      mean_het_rate = cp$mean_het,
      n_informative = nrow(markers)
    )
  }

  # Also check Homo samples for misclassification
  for (samp in c(homo1, homo2)) {
    traj <- compute_het_trajectory(gt_data, markers, samp, het_window)
    if (is.null(traj)) next
    cp <- detect_changepoints(traj)
    if (cp$class != "true_homo") {
      all_calls[[length(all_calls) + 1]] <- data.table(
        chrom = chr, interval_id = iid, sample = samp,
        class = paste0("misclassified_homo_", cp$class),
        n_changepoints = cp$n_changepoints,
        changepoint_positions = paste(cp$changepoint_pos, collapse = ";"),
        max_shift = cp$max_shift,
        mean_het_rate = cp$mean_het,
        n_informative = nrow(markers)
      )
    }
  }

  # Summary
  calls_this <- rbindlist(all_calls[vapply(all_calls, function(x)
    x$chrom[1] == chr & x$interval_id[1] == iid, logical(1))])
  if (nrow(calls_this) > 0) {
    all_summary[[length(all_summary) + 1]] <- data.table(
      chrom = chr, interval_id = iid,
      start_mb = cand$start_mb, end_mb = cand$end_mb,
      n_het_tested = sum(calls_this$class != "true_homo"),
      n_true_het = sum(calls_this$class == "true_het"),
      n_simple_recomb = sum(calls_this$class == "simple_recombinant"),
      n_double_xo = sum(calls_this$class == "double_crossover"),
      n_mosaic = sum(calls_this$class == "complex_mosaic"),
      n_misclass = sum(grepl("misclassified", calls_this$class)),
      recomb_fraction = round(sum(calls_this$class %in% c("simple_recombinant",
        "double_crossover", "complex_mosaic")) / max(1, nrow(calls_this)), 3)
    )
    message("  Results: ", sum(calls_this$class == "true_het"), " true_het, ",
            sum(calls_this$class == "simple_recombinant"), " recomb, ",
            sum(calls_this$class == "double_crossover"), " dbl_XO, ",
            sum(grepl("misclassified", calls_this$class)), " misclass")
  }

  # Plot: het-rate trajectories colored by class
  traj_this <- rbindlist(all_traj[vapply(all_traj, function(x)
    x$chrom[1] == chr & x$interval_id[1] == iid, logical(1))], fill = TRUE)
  if (nrow(traj_this) > 0 && nrow(calls_this) > 0) {
    traj_this <- merge(traj_this, calls_this[, .(sample, class)], by = "sample")
    pT <- ggplot(traj_this, aes(x = pos / 1e6, y = het_rate,
                                 group = sample, color = class)) +
      geom_line(alpha = 0.4, linewidth = 0.3) +
      scale_color_manual(values = c("true_het" = "grey50", "true_homo" = "steelblue",
                                     "simple_recombinant" = "red3",
                                     "double_crossover" = "darkred",
                                     "complex_mosaic" = "purple3"),
                          name = "Class") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
      labs(title = paste0(chr, " I", iid, " -- Het-rate Trajectories"),
           subtitle = paste0(nrow(markers), " informative markers, ",
                            nrow(calls_this), " samples classified"),
           x = paste0(chr, " (Mb)"), y = "Rolling het-rate",
           caption = "Flat ~0.5 = true HET | Step down = recombinant | Two steps = double-XO") +
      theme_minimal(base_size = 9)

    tryCatch(ggsave(file.path(outdir, "plots",
                               paste0(chr, "_I", iid, "_het_trajectories.png")),
                     pT, width = 12, height = 5, dpi = DPI),
             error = function(e) message("  [PLOT] ", e$message))
  }
}

# =============================================================================
# WRITE
# =============================================================================

message("\n[C01h] Writing...")
calls_dt <- if (length(all_calls) > 0) rbindlist(all_calls, fill = TRUE) else data.table()
markers_dt <- if (length(all_markers) > 0) rbindlist(all_markers, fill = TRUE) else data.table()
traj_dt <- if (length(all_traj) > 0) rbindlist(all_traj, fill = TRUE) else data.table()
summ_dt <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()

fwrite(calls_dt, file.path(outdir, "recombinant_calls.tsv.gz"), sep = "\t")
fwrite(markers_dt, file.path(outdir, "informative_markers.tsv.gz"), sep = "\t")
fwrite(traj_dt, file.path(outdir, "het_trajectories.tsv.gz"), sep = "\t")
fwrite(summ_dt, file.path(outdir, "recombinant_summary.tsv"), sep = "\t")

if (nrow(summ_dt) > 0) {
  message("\n[C01h] === RECOMBINANT SUMMARY ===")
  message("  Total candidates scanned: ", nrow(summ_dt))
  message("  True HET: ", sum(summ_dt$n_true_het))
  message("  Simple recombinant: ", sum(summ_dt$n_simple_recomb))
  message("  Double crossover: ", sum(summ_dt$n_double_xo))
  message("  Complex mosaic: ", sum(summ_dt$n_mosaic))
  message("  Misclassified homo: ", sum(summ_dt$n_misclass))
}

message("\n[DONE] -> ", outdir)
