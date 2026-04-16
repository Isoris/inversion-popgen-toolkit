#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01j_regime_compatibility_engine.R  (v8.4)
#
# COMPATIBILITY-BASED REGIME ENGINE
#
# Inspired by hifiasm's site-vector compatibility grouping, adapted for
# population-level inversion detection in a hatchery family soup.
#
# At each sliding window position along a candidate region:
#   1. Build per-sample site vectors from discrete genotypes (012)
#   2. Group samples into COMPATIBLE CLASSES (samples whose vectors
#      can come from the same haplotype arrangement)
#   3. Measure regime properties:
#      - n_groups: how many compatible classes exist
#      - thickness: how many markers support each class
#      - membership: which samples are in which class
#      - persistence: do memberships stay stable across windows?
#   4. Detect regime transitions: where properties change = boundaries
#   5. Label each regime segment: background_soup / structured_inversion /
#      nested_regime / recombinant_zone / transition
#
# Key insight: outside inversions, the family soup creates many thin
# fragmented groups. Inside inversions, samples collapse into fewer
# thick groups. The NUMBER and THICKNESS of compatible groups is the
# regime signal, not PCA eigenvalues.
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>
#   --vcf_dir <dir>              -- Clair3 per-chr VCFs
#   --coseg_dir <dir>            -- C01i output (core packages, sample classes)
#   --samples <sample_list>
#   --outdir <dir>
#   [--window_markers 50]        -- markers per sliding window
#   [--step_markers 10]          -- step size
#
# Outputs:
#   regime_windows.tsv.gz        -- per-window: n_groups, thickness, entropy
#   regime_segments.tsv.gz       -- per-segment: type, span, stability
#   regime_memberships.tsv.gz    -- per-sample per-window: group assignment
#   regime_transitions.tsv       -- boundary positions + transition type
#   regime_state_labels.tsv      -- final per-position state label
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; vcf_dir <- NULL; coseg_dir <- NULL
samples_file <- NULL; outdir <- "regime_engine"
window_markers <- 50L; step_markers <- 10L; tier_max <- 3L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))          { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--vcf_dir" && i < length(args))    { vcf_dir <- args[i+1]; i <- i+2 }
  else if (a == "--coseg_dir" && i < length(args))  { coseg_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--window_markers" && i < length(args)) { window_markers <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--step_markers" && i < length(args))   { step_markers <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--tier_max" && i < length(args))   { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file) || is.null(vcf_dir)) stop("--scores and --vcf_dir required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE)

DPI <- 300
cand_dt <- fread(scores_file)[tier <= tier_max][order(tier, -final_score)]
message("[C01j] Candidates: ", nrow(cand_dt))

# =============================================================================
# GENOTYPE EXTRACTION
# =============================================================================

extract_gt <- function(vcf_file, chr, start_bp, end_bp) {
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
    g <- as.character(raw[[j + 1]])
    gt[, j] <- fifelse(g %in% c("0/0", "0|0"), 0L,
               fifelse(g %in% c("0/1", "1/0", "0|1", "1|0"), 1L,
               fifelse(g %in% c("1/1", "1|1"), 2L, NA_integer_)))
  }
  # Filter: MAF >= 0.03, missing < 40%
  miss <- rowMeans(is.na(gt))
  af <- rowMeans(gt, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  keep <- miss < 0.4 & maf >= 0.03
  if (sum(keep) < 30) return(NULL)
  list(pos = pos[keep], gt = gt[keep, ], samples = samp)
}

# =============================================================================
# COMPATIBILITY GROUPING (core algorithm)
# =============================================================================
# For a set of site vectors (samples x markers, coded 0/1/2):
# Group samples whose vectors are "compatible" = could come from the
# same pair of haplotype arrangements.
#
# Two samples are compatible if their genotype vectors are consistent
# with being the same diplotype class. For inversion detection:
# - Two 0/0 samples: compatible (both REF/REF)
# - Two 1/1 samples: compatible (both INV/INV)
# - Two 0/1 samples: compatible (both HET)
# - 0/0 and 1/1: compatible (different homozygotes)
# - 0/0 and 0/1: compatible (one is HET of the other)
#
# Incompatibility arises when samples SHOULD be in the same group
# but their vectors disagree at too many positions = different arrangements.
#
# We use hierarchical clustering on Hamming distance with adaptive k.

compute_compatibility_groups <- function(gt_window, min_group_size = 3L) {
  ns <- ncol(gt_window); nm <- nrow(gt_window)
  if (ns < 10 || nm < 5) return(NULL)

  # Per-sample dosage mean across window markers
  dosage <- colMeans(gt_window, na.rm = TRUE)

  # Pairwise Hamming distance (on discrete 012)
  dist_mat <- matrix(0, ns, ns)
  for (i in seq_len(ns - 1)) {
    for (j in (i + 1):ns) {
      ok <- !is.na(gt_window[, i]) & !is.na(gt_window[, j])
      nv <- sum(ok)
      if (nv < 5) { dist_mat[i, j] <- 1; dist_mat[j, i] <- 1; next }
      d <- sum(gt_window[ok, i] != gt_window[ok, j]) / nv
      dist_mat[i, j] <- d; dist_mat[j, i] <- d
    }
  }

  # Hierarchical clustering
  hc <- tryCatch(hclust(as.dist(dist_mat), method = "ward.D2"),
                  error = function(e) NULL)
  if (is.null(hc)) return(NULL)

  # Find optimal k using silhouette-like criterion
  # Try k = 2..8, pick k that maximizes between/within ratio
  best_k <- 2L; best_score <- -Inf
  for (k in 2:min(8, ns %/% min_group_size)) {
    cl <- cutree(hc, k = k)
    sizes <- table(cl)
    if (min(sizes) < min_group_size) next

    # Between/within ratio
    within_d <- 0; between_d <- 0; n_within <- 0; n_between <- 0
    for (i in seq_len(ns - 1)) {
      for (j in (i + 1):ns) {
        if (cl[i] == cl[j]) { within_d <- within_d + dist_mat[i, j]; n_within <- n_within + 1 }
        else { between_d <- between_d + dist_mat[i, j]; n_between <- n_between + 1 }
      }
    }
    if (n_within == 0 || n_between == 0) next
    ratio <- (between_d / n_between) / max(within_d / n_within, 0.001)
    if (ratio > best_score) { best_score <- ratio; best_k <- k }
  }

  cl <- cutree(hc, k = best_k)

  # Group properties
  groups <- list()
  for (gi in seq_len(best_k)) {
    g_idx <- which(cl == gi)
    g_gt <- gt_window[, g_idx, drop = FALSE]

    # Thickness: how many markers show consensus (majority agrees)
    consensus <- apply(g_gt, 1, function(row) {
      row <- row[!is.na(row)]
      if (length(row) == 0) return(NA_real_)
      max(tabulate(row + 1L, nbins = 3L)) / length(row)
    })
    thickness <- sum(consensus > 0.7, na.rm = TRUE)

    # Mean dosage (identifies if group is REF-like, HET-like, or INV-like)
    mean_dos <- mean(dosage[g_idx], na.rm = TRUE)

    # Internal consistency (mean pairwise Hamming within group)
    internal_d <- if (length(g_idx) > 1) {
      mean(dist_mat[g_idx, g_idx][upper.tri(dist_mat[g_idx, g_idx])])
    } else 0

    groups[[gi]] <- list(
      group_id = gi, n_samples = length(g_idx), sample_idx = g_idx,
      thickness = thickness, mean_dosage = round(mean_dos, 3),
      internal_distance = round(internal_d, 4),
      dosage_class = if (mean_dos < 0.5) "REF_like"
                     else if (mean_dos > 1.5) "INV_like"
                     else "HET_like"
    )
  }

  list(n_groups = best_k, groups = groups, cluster = cl,
       separation_score = round(best_score, 3),
       dist_mat = dist_mat)
}

# =============================================================================
# SLIDING WINDOW REGIME TRACKING
# =============================================================================

track_regimes <- function(gt_data, window_size, step_size) {
  pos <- gt_data$pos; gt <- gt_data$gt; samples <- gt_data$samples
  nm <- nrow(gt); ns <- ncol(gt)
  if (nm < window_size) return(list(windows = data.table(), memberships = data.table()))

  win_rows <- list(); memb_rows <- list()
  prev_cl <- NULL

  for (wi_start in seq(1, nm - window_size + 1, by = step_size)) {
    wi_end <- wi_start + window_size - 1
    w_gt <- gt[wi_start:wi_end, , drop = FALSE]
    w_pos_start <- pos[wi_start]; w_pos_end <- pos[wi_end]
    w_pos_mid <- (w_pos_start + w_pos_end) / 2

    cg <- compute_compatibility_groups(w_gt)
    if (is.null(cg)) next

    # Regime metrics
    max_thickness <- max(vapply(cg$groups, function(g) g$thickness, integer(1)))
    mean_thickness <- mean(vapply(cg$groups, function(g) g$thickness, integer(1)))
    max_group_size <- max(vapply(cg$groups, function(g) g$n_samples, integer(1)))
    min_group_size <- min(vapply(cg$groups, function(g) g$n_samples, integer(1)))

    # Group entropy: how evenly distributed are samples across groups?
    sizes <- vapply(cg$groups, function(g) g$n_samples, integer(1))
    props <- sizes / sum(sizes)
    entropy <- -sum(props * log2(pmax(props, 1e-10)))

    # Persistence: compare to previous window
    persistence <- NA_real_
    if (!is.null(prev_cl) && length(prev_cl) == length(cg$cluster)) {
      # ARI between consecutive windows
      tab <- table(prev_cl, cg$cluster)
      n <- sum(tab)
      if (n >= 2) {
        sc <- sum(choose(tab, 2))
        sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
        e <- sa * sb / choose(n, 2); mx <- 0.5 * (sa + sb)
        persistence <- if (mx == e) 1 else (sc - e) / (mx - e)
      }
    }
    prev_cl <- cg$cluster

    # Dosage class distribution
    n_ref <- sum(vapply(cg$groups, function(g) g$dosage_class == "REF_like", logical(1)) *
                  vapply(cg$groups, function(g) g$n_samples, integer(1)))
    n_het <- sum(vapply(cg$groups, function(g) g$dosage_class == "HET_like", logical(1)) *
                  vapply(cg$groups, function(g) g$n_samples, integer(1)))
    n_inv <- sum(vapply(cg$groups, function(g) g$dosage_class == "INV_like", logical(1)) *
                  vapply(cg$groups, function(g) g$n_samples, integer(1)))

    win_rows[[length(win_rows) + 1]] <- data.table(
      window = length(win_rows) + 1L,
      pos_start = w_pos_start, pos_end = w_pos_end,
      pos_mid_mb = round(w_pos_mid / 1e6, 4),
      n_markers = window_size,
      n_groups = cg$n_groups,
      max_thickness = max_thickness,
      mean_thickness = round(mean_thickness, 1),
      max_group_size = max_group_size,
      min_group_size = min_group_size,
      entropy = round(entropy, 3),
      separation = cg$separation_score,
      persistence = round(persistence, 3),
      n_ref_like = n_ref, n_het_like = n_het, n_inv_like = n_inv
    )

    # Per-sample group membership
    for (si in seq_len(ns)) {
      gi <- cg$cluster[si]
      grp <- cg$groups[[gi]]
      memb_rows[[length(memb_rows) + 1]] <- data.table(
        window = length(win_rows),
        sample = samples[si],
        group_id = gi,
        dosage_class = grp$dosage_class,
        pos_mid_mb = round(w_pos_mid / 1e6, 4)
      )
    }
  }

  list(
    windows = if (length(win_rows) > 0) rbindlist(win_rows) else data.table(),
    memberships = if (length(memb_rows) > 0) rbindlist(memb_rows) else data.table()
  )
}

# =============================================================================
# REGIME SEGMENTATION + STATE LABELING
# =============================================================================

segment_regimes <- function(win_dt) {
  if (nrow(win_dt) == 0) return(list(segments = data.table(), transitions = data.table()))

  # Label each window based on STRUCTURAL QUALITY, not group count.
  # Thick + high separation + persistent = structured signal (inversion)
  # Thin + low separation + unstable = unstructured noise (soup)
  # Group count tells you COMPLEXITY, not quality.
  win_dt[, structure_score := (
    pmin(1, mean_thickness / 25) * 0.4 +     # thick groups = structured
    pmin(1, separation / 5) * 0.35 +           # well-separated = real
    pmin(1, fifelse(is.finite(persistence), persistence, 0)) * 0.25  # stable = real
  )]

  win_dt[, complexity := fifelse(
    n_groups <= 2, "simple",
    fifelse(n_groups <= 4, "moderate",
    fifelse(n_groups <= 6, "complex", "highly_complex"))
  )]

  win_dt[, state := fifelse(
    structure_score > 0.6 & complexity == "simple",
    "clean_inversion",
    fifelse(
      structure_score > 0.6 & complexity == "moderate",
      "structured_moderate",
      fifelse(
        structure_score > 0.6 & complexity %in% c("complex", "highly_complex"),
        "structured_complex",
        fifelse(
          structure_score > 0.35 & structure_score <= 0.6,
          "weak_signal",
          fifelse(
            structure_score <= 0.35 & mean_thickness < 8,
            "background_soup",
            "transition"
          )
        )
      )
    )
  )]

  # Run-length encode states into segments
  rle_state <- rle(win_dt$state)
  seg_rows <- list(); pos <- 0L
  for (ri in seq_along(rle_state$lengths)) {
    len <- rle_state$lengths[ri]; state <- rle_state$values[ri]
    seg_start <- pos + 1; seg_end <- pos + len; pos <- pos + len
    if (seg_start > nrow(win_dt) || seg_end > nrow(win_dt)) next

    seg_rows[[ri]] <- data.table(
      segment_id = ri,
      state = state,
      n_windows = len,
      start_mb = win_dt$pos_mid_mb[seg_start],
      end_mb = win_dt$pos_mid_mb[seg_end],
      mean_n_groups = round(mean(win_dt$n_groups[seg_start:seg_end]), 1),
      mean_thickness = round(mean(win_dt$mean_thickness[seg_start:seg_end]), 1),
      mean_separation = round(mean(win_dt$separation[seg_start:seg_end]), 2),
      mean_persistence = round(mean(win_dt$persistence[seg_start:seg_end], na.rm = TRUE), 3),
      mean_entropy = round(mean(win_dt$entropy[seg_start:seg_end]), 3)
    )
  }
  seg_dt <- if (length(seg_rows) > 0) rbindlist(seg_rows) else data.table()

  # Transitions: where state changes
  trans_rows <- list()
  if (nrow(seg_dt) >= 2) {
    for (ti in seq_len(nrow(seg_dt) - 1)) {
      s1 <- seg_dt[ti]; s2 <- seg_dt[ti + 1]
      group_change <- s2$mean_n_groups - s1$mean_n_groups
      thickness_change <- s2$mean_thickness - s1$mean_thickness

      trans_type <- if (s1$state == "background_soup" & grepl("structured|clean", s2$state)) {
        "enter_structured"
      } else if (grepl("structured|clean", s1$state) & s2$state == "background_soup") {
        "exit_structured"
      } else if (s1$state == "clean_inversion" & s2$state == "structured_complex") {
        "complexity_increase"
      } else if (s1$state == "structured_complex" & s2$state == "clean_inversion") {
        "complexity_decrease"
      } else if (grepl("structured|clean", s1$state) & s2$state == "weak_signal") {
        "signal_weakening"
      } else if (s1$state == "weak_signal" & grepl("structured|clean", s2$state)) {
        "signal_strengthening"
      } else if (s1$state != s2$state) {
        "regime_change"
      } else {
        "continuation"
      }

      trans_rows[[ti]] <- data.table(
        from_segment = ti, to_segment = ti + 1,
        from_state = s1$state, to_state = s2$state,
        boundary_mb = round((s1$end_mb + s2$start_mb) / 2, 4),
        group_change = round(group_change, 1),
        thickness_change = round(thickness_change, 1),
        transition_type = trans_type
      )
    }
  }
  trans_dt <- if (length(trans_rows) > 0) rbindlist(trans_rows) else data.table()

  list(segments = seg_dt, transitions = trans_dt)
}

# =============================================================================
# MEMBERSHIP STABILITY ANALYSIS (who's inside, and do they stay?)
# =============================================================================

analyze_membership_stability <- function(memb_dt, win_dt) {
  if (nrow(memb_dt) == 0) return(data.table())

  # For each sample, compute: how often does it stay in the same group
  # across consecutive windows?
  samples <- unique(memb_dt$sample)
  windows <- sort(unique(memb_dt$window))

  stab_rows <- list()
  for (samp in samples) {
    s_memb <- memb_dt[sample == samp][order(window)]
    if (nrow(s_memb) < 3) next

    # Count group switches
    switches <- sum(s_memb$group_id[-1] != s_memb$group_id[-nrow(s_memb)])
    switch_rate <- switches / (nrow(s_memb) - 1)

    # Dominant group
    dom_group <- names(sort(table(s_memb$group_id), decreasing = TRUE))[1]
    dom_frac <- mean(s_memb$group_id == dom_group)

    # Dominant dosage class
    dom_class <- names(sort(table(s_memb$dosage_class), decreasing = TRUE))[1]
    class_frac <- mean(s_memb$dosage_class == dom_class)

    # Is this sample a recombinant? (switches dosage class mid-region)
    class_switches <- sum(s_memb$dosage_class[-1] != s_memb$dosage_class[-nrow(s_memb)])

    sample_type <- if (class_switches == 0 && dom_class == "HET_like") "stable_het"
                   else if (class_switches == 0 && dom_class == "REF_like") "stable_ref"
                   else if (class_switches == 0 && dom_class == "INV_like") "stable_inv"
                   else if (class_switches == 1) "simple_recombinant"
                   else if (class_switches == 2) "double_crossover"
                   else if (class_switches >= 3) "complex_mosaic"
                   else "unstable"

    stab_rows[[length(stab_rows) + 1]] <- data.table(
      sample = samp,
      n_windows = nrow(s_memb),
      n_switches = switches,
      switch_rate = round(switch_rate, 3),
      dominant_group = dom_group,
      dominant_frac = round(dom_frac, 3),
      dominant_class = dom_class,
      class_frac = round(class_frac, 3),
      class_switches = class_switches,
      sample_type = sample_type
    )
  }
  if (length(stab_rows) > 0) rbindlist(stab_rows) else data.table()
}

# =============================================================================
# MAIN
# =============================================================================

all_windows <- list(); all_segments <- list()
all_transitions <- list(); all_stability <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]; chr <- cand$chrom; iid <- cand$interval_id
  message("\n[C01j] ", ci, "/", nrow(cand_dt), ": ", chr, " I", iid,
          " (", cand$start_mb, "-", cand$end_mb, " Mb)")

  vcf_file <- file.path(vcf_dir, paste0(chr, ".vcf.gz"))
  if (!file.exists(vcf_file)) {
    vcf_file <- list.files(vcf_dir, pattern = paste0(chr, ".*\\.vcf\\.gz$"),
                            full.names = TRUE)[1]
    if (is.na(vcf_file)) { message("  No VCF"); next }
  }

  # Extract with flanks (20% on each side)
  span <- cand$end_mb - cand$start_mb
  flank <- span * 0.2
  gt_data <- extract_gt(vcf_file, chr,
                         (cand$start_mb - flank) * 1e6,
                         (cand$end_mb + flank) * 1e6)
  if (is.null(gt_data)) { message("  Extraction failed"); next }
  message("  Markers: ", nrow(gt_data$gt), " x ", ncol(gt_data$gt), " samples")

  # Track regimes
  regime <- track_regimes(gt_data, window_markers, step_markers)
  if (nrow(regime$windows) == 0) { message("  No windows"); next }
  regime$windows[, `:=`(chrom = chr, interval_id = iid)]
  message("  Windows: ", nrow(regime$windows),
          " | Groups range: ", min(regime$windows$n_groups), "-",
          max(regime$windows$n_groups))

  # Segment + label
  seg <- segment_regimes(regime$windows)
  if (nrow(seg$segments) > 0) {
    seg$segments[, `:=`(chrom = chr, interval_id = iid)]
    message("  Segments: ", nrow(seg$segments), " (",
            paste(names(table(seg$segments$state)),
                  table(seg$segments$state), sep = "=", collapse = " "), ")")
  }
  if (nrow(seg$transitions) > 0) {
    seg$transitions[, `:=`(chrom = chr, interval_id = iid)]
    message("  Transitions: ", nrow(seg$transitions), " (",
            paste(names(table(seg$transitions$transition_type)),
                  table(seg$transitions$transition_type), sep = "=", collapse = " "), ")")
  }

  # Membership stability
  regime$memberships[, `:=`(chrom = chr, interval_id = iid)]
  stab <- analyze_membership_stability(regime$memberships, regime$windows)
  if (nrow(stab) > 0) {
    stab[, `:=`(chrom = chr, interval_id = iid)]
    n_recomb <- sum(stab$sample_type == "simple_recombinant")
    n_dxo <- sum(stab$sample_type == "double_crossover")
    n_mosaic <- sum(stab$sample_type == "complex_mosaic")
    message("  Stability: ", sum(stab$sample_type == "stable_het"), " stable_het, ",
            sum(stab$sample_type == "stable_ref"), " stable_ref, ",
            sum(stab$sample_type == "stable_inv"), " stable_inv, ",
            n_recomb, " recomb, ", n_dxo, " dbl_XO, ", n_mosaic, " mosaic")
  }

  all_windows[[length(all_windows) + 1]] <- regime$windows
  if (nrow(seg$segments) > 0) all_segments[[length(all_segments) + 1]] <- seg$segments
  if (nrow(seg$transitions) > 0) all_transitions[[length(all_transitions) + 1]] <- seg$transitions
  if (nrow(stab) > 0) all_stability[[length(all_stability) + 1]] <- stab

  # =================================================================
  # DIAGNOSTIC PLOTS (5 panels)
  # =================================================================
  if (nrow(regime$windows) > 5 && nrow(regime$memberships) > 0) {
    wd <- regime$windows
    md <- regime$memberships

    # --- Panel 1: Regime state strip ---
    pR <- ggplot(wd) +
      geom_tile(aes(x = pos_mid_mb, y = 1, fill = state), height = 0.8) +
      scale_fill_manual(values = c(
        "clean_inversion" = "#dc2626", "structured_moderate" = "#ea580c",
        "structured_complex" = "#7c3aed", "weak_signal" = "#d97706",
        "background_soup" = "#d1d5db", "transition" = "#f3f4f6"),
        name = "State") +
      labs(title = paste0(chr, " I", iid, " -- Regime Profile"),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 6),
            plot.title = element_text(size = 10, face = "bold"))

    # --- Panel 2: Vector-type mosaic ---
    # Each position has N compatible groups (not just 3).
    # Each sample belongs to one group at each position.
    # The mosaic shows: for each regime interval, what are the
    # distinct vector types and how many samples belong to each?
    #
    # Samples sorted by their DOMINANT group across windows, so
    # samples that stay in the same group are contiguous rows.
    # Color = group_id (not just REF/HET/INV — actual compatible class).
    # When complexity changes (entering overlap zone), new colors appear.
    # When a recombinant exists, its row switches color mid-chromosome.

    # Assign persistent group labels across windows
    # (raw group_id can shuffle between windows, need to align them)
    # Use dosage_class + group_id combination as stable label
    md[, vector_label := paste0(dosage_class, "_g", group_id)]

    # Count unique vector types
    n_vtypes <- length(unique(md$vector_label))

    # Build a color palette for vector types
    # REF-like groups get blue shades, HET-like get grey/green, INV-like get red shades
    vtype_levels <- sort(unique(md$vector_label))
    vtype_colors <- character(length(vtype_levels))
    ref_blues <- c("#1e3a5f", "#2563eb", "#60a5fa", "#93c5fd")
    het_greys <- c("#374151", "#6b7280", "#9ca3af", "#d1d5db")
    inv_reds  <- c("#7f1d1d", "#dc2626", "#f87171", "#fca5a5")
    ri <- 1; hi <- 1; ii <- 1
    for (vi in seq_along(vtype_levels)) {
      vt <- vtype_levels[vi]
      if (grepl("REF_like", vt)) {
        vtype_colors[vi] <- ref_blues[pmin(ri, length(ref_blues))]; ri <- ri + 1
      } else if (grepl("HET_like", vt)) {
        vtype_colors[vi] <- het_greys[pmin(hi, length(het_greys))]; hi <- hi + 1
      } else {
        vtype_colors[vi] <- inv_reds[pmin(ii, length(inv_reds))]; ii <- ii + 1
      }
    }
    names(vtype_colors) <- vtype_levels

    # Order samples by: sample_type then dominant vector label then stability
    if (nrow(stab) > 0) {
      # Find each sample's dominant vector label
      dom_vtype <- md[, .(dom_vtype = names(sort(table(vector_label), decreasing = TRUE))[1]),
                       by = sample]
      stab2 <- merge(stab, dom_vtype, by = "sample", all.x = TRUE)
      type_order <- c("stable_ref", "stable_inv", "stable_het",
                       "simple_recombinant", "double_crossover",
                       "complex_mosaic", "unstable")
      stab2[, type_rank := match(sample_type, type_order, nomatch = 99)]
      samp_order <- stab2[order(type_rank, dom_vtype, -class_frac)]$sample
      md[, sample := factor(sample, levels = samp_order)]
    }

    # Subsample windows for plotting
    n_win_plot <- length(unique(md$window))
    if (n_win_plot > 300) {
      keep_wins <- unique(md$window)[seq(1, n_win_plot, length.out = 300)]
      md_plot <- md[window %in% keep_wins]
    } else md_plot <- copy(md)

    pMosaic <- ggplot(md_plot, aes(x = pos_mid_mb, y = sample, fill = vector_label)) +
      geom_tile() +
      scale_fill_manual(values = vtype_colors, name = "Vector\ntype") +
      labs(x = NULL, y = NULL,
           subtitle = paste0("Sample vector-type mosaic (", n_vtypes,
                            " types; sorted: stable top, recombinants bottom)")) +
      theme_minimal(base_size = 7) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 5))

    # Add separators between sample type groups
    if (nrow(stab) > 0 && exists("stab2")) {
      ordered_types <- stab2[match(samp_order, sample)]$sample_type
      type_breaks <- which(ordered_types[-length(ordered_types)] != ordered_types[-1])
      for (tb in type_breaks) {
        pMosaic <- pMosaic +
          geom_hline(yintercept = tb + 0.5, color = "white", linewidth = 0.4)
      }
    }

    # --- Panel 3: Stacked area — vector-type proportions along chromosome ---
    # Each vector type gets its own area in the stack.
    # When complexity changes, new types appear/disappear.
    # Recombinant zones show type proportions shifting.
    prop_vtype <- md[, .N, by = .(window, pos_mid_mb, vector_label)]
    prop_vtype[, total := sum(N), by = window]
    prop_vtype[, fraction := N / total]
    prop_vtype[, vector_label := factor(vector_label, levels = vtype_levels)]

    pArea <- ggplot(prop_vtype, aes(x = pos_mid_mb, y = fraction, fill = vector_label)) +
      geom_area(alpha = 0.85, color = NA, position = "stack") +
      scale_fill_manual(values = vtype_colors, name = "Vector\ntype") +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
      labs(x = NULL, y = "Proportion",
           subtitle = "Population composition by vector type (new colors = new regimes)") +
      theme_minimal(base_size = 8) +
      theme(legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 5))

    # --- Panel 3: Multi-signal tracks (like enrichment plot) ---
    # Structure score, thickness, separation, entropy all on same axes
    # Smoothed with loess for cleaner visualization
    signal_dt <- wd[, .(pos_mid_mb, structure_score, n_groups,
                         mean_thickness, separation, entropy)]

    pS <- ggplot(signal_dt, aes(x = pos_mid_mb)) +
      geom_ribbon(aes(ymin = 0, ymax = structure_score), fill = "#dc2626", alpha = 0.15) +
      geom_smooth(aes(y = structure_score, color = "Structure score"),
                   method = "loess", span = 0.2, se = FALSE, linewidth = 0.8) +
      geom_smooth(aes(y = pmin(1, mean_thickness / 30), color = "Thickness (norm)"),
                   method = "loess", span = 0.2, se = FALSE, linewidth = 0.7) +
      geom_smooth(aes(y = pmin(1, separation / 6), color = "Separation (norm)"),
                   method = "loess", span = 0.2, se = FALSE, linewidth = 0.7) +
      geom_smooth(aes(y = entropy / max(entropy + 0.01, na.rm = TRUE), color = "Entropy (norm)"),
                   method = "loess", span = 0.2, se = FALSE, linewidth = 0.6,
                   linetype = "dashed") +
      scale_color_manual(values = c("Structure score" = "#dc2626",
                                     "Thickness (norm)" = "#ea580c",
                                     "Separation (norm)" = "#2563eb",
                                     "Entropy (norm)" = "#6b7280"),
                          name = "Signal") +
      scale_y_continuous(limits = c(0, 1.05)) +
      labs(x = NULL, y = "Normalized\nsignal",
           subtitle = "Multi-signal regime tracks") +
      theme_minimal(base_size = 8) +
      theme(legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 6))

    # --- Panel 5: N groups + persistence ---
    pG <- ggplot(wd, aes(x = pos_mid_mb)) +
      geom_step(aes(y = n_groups), color = "#2563eb", linewidth = 0.5) +
      geom_line(data = wd[is.finite(persistence)],
                 aes(y = persistence * max(n_groups, na.rm = TRUE)),
                 color = "#16a34a", linewidth = 0.5, alpha = 0.7) +
      scale_y_continuous(sec.axis = sec_axis(
        ~ . / max(wd$n_groups, na.rm = TRUE), name = "Persistence")) +
      labs(x = paste0(chr, " (Mb)"), y = "N groups",
           subtitle = "Blue = compatible groups | Green = persistence (ARI)") +
      theme_minimal(base_size = 8)

    # Assemble 5-panel figure
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      p_all <- pR / pMosaic / pArea / pS / pG +
        plot_layout(heights = c(0.5, 3, 1.2, 1.5, 1))

      tryCatch(ggsave(file.path(outdir, "plots",
                                 paste0(chr, "_I", iid, "_regimes.png")),
                       p_all, width = 14, height = 18, dpi = DPI),
               error = function(e) message("  [PLOT] ", e$message))
    }
  }
}

# =============================================================================
# WRITE
# =============================================================================

message("\n[C01j] Writing...")
w_dt <- if (length(all_windows) > 0) rbindlist(all_windows, fill = TRUE) else data.table()
s_dt <- if (length(all_segments) > 0) rbindlist(all_segments, fill = TRUE) else data.table()
t_dt <- if (length(all_transitions) > 0) rbindlist(all_transitions, fill = TRUE) else data.table()
st_dt <- if (length(all_stability) > 0) rbindlist(all_stability, fill = TRUE) else data.table()

fwrite(w_dt, file.path(outdir, "regime_windows.tsv.gz"), sep = "\t")
fwrite(s_dt, file.path(outdir, "regime_segments.tsv.gz"), sep = "\t")
fwrite(t_dt, file.path(outdir, "regime_transitions.tsv"), sep = "\t")
fwrite(st_dt, file.path(outdir, "regime_stability.tsv.gz"), sep = "\t")

# State labels: combine segment states into per-position labels
if (nrow(w_dt) > 0) {
  labels <- w_dt[, .(chrom, interval_id, pos_mid_mb, state, complexity,
                      structure_score, n_groups, mean_thickness, separation,
                      persistence, entropy,
                      n_ref_like, n_het_like, n_inv_like)]
  fwrite(labels, file.path(outdir, "regime_state_labels.tsv.gz"), sep = "\t")
}

# =================================================================
# REGIME ROARY-STYLE WINDOW PA MATRIX
# =================================================================
# For every window in each candidate: regime state + membership summary
# This can be joined with the core PA and merge PA by global_window_id

regime_pa_rows <- list()
if (nrow(w_dt) > 0 && nrow(st_dt) > 0) {
  for (chr in unique(w_dt$chrom)) {
    chr_w <- w_dt[chrom == chr]
    chr_st <- st_dt[chrom == chr]
    for (iid in unique(chr_w$interval_id)) {
      iid_w <- chr_w[interval_id == iid]
      iid_st <- chr_st[interval_id == iid]

      # Per-sample type counts
      n_stable_ref <- sum(iid_st$sample_type == "stable_ref")
      n_stable_het <- sum(iid_st$sample_type == "stable_het")
      n_stable_inv <- sum(iid_st$sample_type == "stable_inv")
      n_recomb <- sum(iid_st$sample_type == "simple_recombinant")
      n_dxo <- sum(iid_st$sample_type == "double_crossover")
      n_mosaic <- sum(iid_st$sample_type == "complex_mosaic")

      for (wi in seq_len(nrow(iid_w))) {
        w <- iid_w[wi]
        regime_pa_rows[[length(regime_pa_rows) + 1]] <- data.table(
          chrom = chr,
          interval_id = iid,
          regime_window = w$window,
          pos_mid_mb = w$pos_mid_mb,
          regime_state = w$state,
          regime_complexity = w$complexity,
          regime_structure_score = w$structure_score,
          regime_n_groups = w$n_groups,
          regime_thickness = w$mean_thickness,
          regime_separation = w$separation,
          regime_persistence = w$persistence,
          regime_entropy = w$entropy,
          frac_ref = round(w$n_ref_like / max(1, w$n_ref_like + w$n_het_like + w$n_inv_like), 3),
          frac_het = round(w$n_het_like / max(1, w$n_ref_like + w$n_het_like + w$n_inv_like), 3),
          frac_inv = round(w$n_inv_like / max(1, w$n_ref_like + w$n_het_like + w$n_inv_like), 3),
          # Candidate-level sample type summary
          cand_n_stable_ref = n_stable_ref,
          cand_n_stable_het = n_stable_het,
          cand_n_stable_inv = n_stable_inv,
          cand_n_recomb = n_recomb,
          cand_n_dxo = n_dxo,
          cand_n_mosaic = n_mosaic
        )
      }
    }
  }
}

regime_pa_dt <- if (length(regime_pa_rows) > 0) rbindlist(regime_pa_rows, fill = TRUE) else data.table()
fwrite(regime_pa_dt, file.path(outdir, "regime_window_pa.tsv.gz"), sep = "\t")

if (nrow(st_dt) > 0) {
  message("\n[C01j] === SAMPLE TYPE SUMMARY ===")
  for (v in sort(unique(st_dt$sample_type)))
    message("  ", v, ": ", sum(st_dt$sample_type == v))
}

# =============================================================================
# COMBINED ALL-LAYERS STRIP PLOT
# =============================================================================
# One figure per candidate showing ALL pipeline layers as horizontal strips:
#   Strip 1: raw sim_mat diagonal profile (eigenvalue or mean kNN sim)
#   Strip 2: core PA (SML / SM / L / seed_uncollected / unused)
#   Strip 3: merge PA (ABC / AB / A / bridge / gap_with_signal / background)
#   Strip 4: triangle PA (strong / moderate / patchy / diffuse / outside)
#   Strip 5: regime state (clean_inversion / structured_* / weak / soup)
#   Strip 6: tier assignment (if available)
#
# Loads PA files from other scripts via --core_pa, --merge_pa, --tri_pa

core_pa_file <- NULL; merge_pa_file <- NULL; tri_pa_file <- NULL
ci <- match("--core_pa", args); if (!is.na(ci) && ci < length(args)) core_pa_file <- args[ci + 1]
mi <- match("--merge_pa", args); if (!is.na(mi) && mi < length(args)) merge_pa_file <- args[mi + 1]
tp <- match("--tri_pa", args); if (!is.na(tp) && tp < length(args)) tri_pa_file <- args[tp + 1]
sc <- match("--scores", args)  # already parsed above

core_pa <- if (!is.null(core_pa_file) && file.exists(core_pa_file)) fread(core_pa_file) else data.table()
merge_pa <- if (!is.null(merge_pa_file) && file.exists(merge_pa_file)) fread(merge_pa_file) else data.table()
tri_pa <- if (!is.null(tri_pa_file) && file.exists(tri_pa_file)) fread(tri_pa_file) else data.table()

has_all_pa <- nrow(core_pa) > 0 && nrow(merge_pa) > 0 && nrow(tri_pa) > 0

if (has_all_pa && nrow(regime_pa_dt) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  message("[C01j] Building combined all-layers strip plots...")

  for (chr_iid in unique(paste0(regime_pa_dt$chrom, "_", regime_pa_dt$interval_id))) {
    parts <- strsplit(chr_iid, "_(?=[^_]+$)", perl = TRUE)[[1]]
    chr <- parts[1]; iid <- as.integer(parts[2])

    cand_regime <- regime_pa_dt[chrom == chr & interval_id == iid]
    if (nrow(cand_regime) == 0) next

    pos_range <- range(cand_regime$pos_mid_mb)
    flank <- diff(pos_range) * 0.15

    # Filter all PAs to this region
    c_pa <- core_pa[chrom == chr & pos_mb >= pos_range[1] - flank & pos_mb <= pos_range[2] + flank]
    m_pa <- merge_pa[chrom == chr & pos_mb >= pos_range[1] - flank & pos_mb <= pos_range[2] + flank]
    t_pa <- tri_pa[chrom == chr & pos_mb >= pos_range[1] - flank & pos_mb <= pos_range[2] + flank]

    # Build strip data
    strips <- list()

    # Strip 1: Core PA
    if (nrow(c_pa) > 0 && "pa_pattern" %in% names(c_pa)) {
      core_map <- c("SML" = 5, "SM" = 4, "ML" = 3.5, "SL" = 3.5,
                     "S" = 3, "M" = 3, "L" = 2, "none" = 0)
      c_pa[, strip_val := core_map[pa_pattern] %||% 0]
      c_pa[collector_state == "seed_uncollected", strip_val := 1]
      strips[[1]] <- c_pa[, .(pos_mb, value = strip_val / 5,
                               strip = "1_cores")]
    }

    # Strip 2: Merge PA
    if (nrow(m_pa) > 0 && "merge_pa" %in% names(m_pa)) {
      merge_map <- c("ABC" = 5, "AB" = 4, "AC" = 3.5, "BC" = 3.5,
                      "A" = 3, "B" = 2.5, "C" = 2, "none" = 0)
      m_pa[, strip_val := merge_map[merge_pa] %||% 0]
      m_pa[merge_role == "bridge", strip_val := 1.5]
      m_pa[merge_role == "gap_with_signal", strip_val := 1]
      strips[[2]] <- m_pa[, .(pos_mb, value = strip_val / 5,
                               strip = "2_merge")]
    }

    # Strip 3: Triangle PA
    if (nrow(t_pa) > 0 && "tri_state" %in% names(t_pa)) {
      tri_map <- c("strong_triangle" = 5, "moderate_triangle" = 4,
                    "sharp_but_not_square" = 3, "patchy_signal" = 2.5,
                    "diffuse_zone" = 2, "weak_zone" = 1, "outside" = 0,
                    "transition" = 0.5)
      t_pa[, strip_val := tri_map[tri_state] %||% 0]
      strips[[3]] <- t_pa[, .(pos_mb, value = strip_val / 5,
                               strip = "3_triangles")]
    }

    # Strip 4: Regime state
    regime_map <- c("clean_inversion" = 5, "structured_moderate" = 4,
                     "structured_complex" = 3.5, "weak_signal" = 2,
                     "background_soup" = 0.5, "transition" = 1)
    cand_regime[, strip_val := regime_map[regime_state] %||% 0]
    strips[[4]] <- cand_regime[, .(pos_mb = pos_mid_mb, value = strip_val / 5,
                                    strip = "4_regimes")]

    if (length(strips) >= 2) {
      all_strips <- rbindlist(strips, fill = TRUE)

      p_combined <- ggplot(all_strips, aes(x = pos_mb, y = strip, fill = value)) +
        geom_tile(height = 0.8) +
        scale_fill_gradientn(
          colours = c("grey90", "grey70", "gold3", "darkorange", "red3", "darkred"),
          limits = c(0, 1), name = "Signal") +
        scale_y_discrete(labels = c("1_cores" = "Snake cores",
                                     "2_merge" = "Fuzzy merge",
                                     "3_triangles" = "Triangles",
                                     "4_regimes" = "Regimes")) +
        labs(title = paste0(chr, " I", iid, " -- All Pipeline Layers"),
             subtitle = "Each strip = one detection/analysis layer. Red = strong signal. Grey = nothing.",
             x = paste0(chr, " (Mb)"), y = NULL,
             caption = paste0("Cores: red=SML,orange=2fam,blue=1fam,gold=seed_unc | ",
                             "Merge: red=ABC,orange=AB | ",
                             "Triangles: red=strong,orange=moderate | ",
                             "Regimes: red=clean,orange=structured")) +
        theme_minimal(base_size = 9) +
        theme(plot.title = element_text(size = 10, face = "bold"),
              plot.caption = element_text(size = 5, color = "grey60", hjust = 0))

      tryCatch(
        ggsave(file.path(outdir, "plots",
                           paste0(chr, "_I", iid, "_A5_all_layers.png")),
               p_combined, width = 16, height = 4, dpi = 300),
        error = function(e) message("  [PLOT] ", e$message))
    }
  }
}

message("\n[DONE] -> ", outdir)
