#!/usr/bin/env Rscript

# =============================================================================
# STEP40_candidate_internal_coherence.R  (v1.0)
#
# PER-CANDIDATE INTERNAL COHERENCE DIAGNOSTIC
#
# Asks: is this candidate one coherent system, or a composite of multiple
# regimes that were incorrectly merged?
#
# WORKFLOW (in this order, not reversed):
#   Step 1: Run Snake 2 + Snake 3 metrics across ALL windows of the full
#           candidate. Do NOT pre-split by samples.
#   Step 2: Detect internal coherence breaks — stretches where Snake 2/3
#           continuity collapses
#   Step 3: Propose window-defined sub-candidates if breaks are found
#   Step 4: Only then compare sample roles across sub-candidates
#
# WHY THIS ORDER:
#   - Pre-splitting by samples forces the answer too early
#   - It can cut apart real heterokaryotype structure
#   - It makes Snake 2/3 look cleaner just because you pre-separated
#   - You lose the evidence that the candidate was composite
#
# INPUTS:
#   <config.R>                     — standard followup config
#   <candidate_id>                 — candidate to analyze
#   <snake2_dir>                   — Snake 2 track output
#   <snake3_dir>                   — Snake 3 track output
#   <snake_dir>                    — Snake 1 multiscale output
#   <consensus_dir>                — 3-snake consensus
#
# OUTPUTS:
#   candidate_coherence_profile.tsv.gz   — per-window coherence scores
#   candidate_internal_breaks.tsv        — detected break positions
#   candidate_sub_candidates.tsv         — proposed sub-candidates
#   candidate_role_comparison.tsv        — sample role comparison across subs
#   candidate_coherence_verdict.tsv      — final verdict
#   candidate_coherence_diagnostic.pdf   — multi-panel diagnostic figure
#
# Usage:
#   Rscript STEP40_candidate_internal_coherence.R <config.R> <cid> \
#     [--snake2_dir <path>] [--snake3_dir <path>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: STEP40_candidate_internal_coherence.R <config.R> <cid> ...")

config_file <- args[1]
cid <- as.integer(args[2])

source(config_file)

# Parse optional args
snake2_dir    <- file.path(INV_ROOT, "06_mds_candidates", "snake2_community")
snake3_dir    <- file.path(INV_ROOT, "06_mds_candidates", "snake3_ghsl")
snake_dir     <- file.path(INV_ROOT, "06_mds_candidates", "snake_regions_multiscale")
consensus_dir <- file.path(INV_ROOT, "06_mds_candidates", "three_snake_consensus")

# Configurable thresholds (not hardcoded biological truth)
BREAK_THRESH     <- 0.30   # combined score below this = potential break
MIN_COHERENT_RUN <- 3L     # min coherent windows to define a sub-candidate
MIN_BREAK_RUN    <- 2L     # min incoherent windows to count as a real break
EDGE_FRAC        <- 0.15   # breaks within this fraction of candidate edges = edge-only

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--snake2_dir" && i < length(args)) { snake2_dir <- args[i+1]; i <- i+2L }
  else if (a == "--snake3_dir" && i < length(args)) { snake3_dir <- args[i+1]; i <- i+2L }
  else if (a == "--snake_dir" && i < length(args)) { snake_dir <- args[i+1]; i <- i+2L }
  else if (a == "--break_thresh" && i < length(args)) { BREAK_THRESH <- as.numeric(args[i+1]); i <- i+2L }
  else if (a == "--min_coherent" && i < length(args)) { MIN_COHERENT_RUN <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--min_break" && i < length(args)) { MIN_BREAK_RUN <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--edge_frac" && i < length(args)) { EDGE_FRAC <- as.numeric(args[i+1]); i <- i+2L }
  else i <- i + 1L
}

message("[STEP40] Thresholds: break=", BREAK_THRESH, " min_coherent=", MIN_COHERENT_RUN,
        " min_break=", MIN_BREAK_RUN, " edge_frac=", EDGE_FRAC)

# Load candidate
cand <- fread(CANDIDATE_TABLE)
row <- cand[candidate_id == cid]
if (nrow(row) == 0) stop("Candidate ", cid, " not found")
chr <- row$chrom; c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))
cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)

message("[STEP40] Internal coherence diagnostic for candidate ", cid,
        " (", chr, ":", round(c_start/1e6, 2), "–", round(c_end/1e6, 2), " Mb)")

# =============================================================================
# STEP 1: Load Snake 2/3 tracks for the candidate region
# =============================================================================

safe_load <- function(path) {
  if (file.exists(path)) tryCatch(fread(path), error = function(e) NULL) else NULL
}

s2 <- safe_load(file.path(snake2_dir, "snake2_track.tsv.gz"))
s3 <- safe_load(file.path(snake3_dir, "snake3_track.tsv.gz"))
s1 <- safe_load(file.path(snake_dir, "snake_windows.tsv.gz"))
cons <- safe_load(file.path(consensus_dir, "consensus_windows.tsv.gz"))

# Filter to candidate region
filter_region <- function(dt, chr, start, end) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  dt[chrom == chr & start_bp >= start & end_bp <= end]
}

s2_cand <- filter_region(s2, chr, c_start, c_end)
s3_cand <- filter_region(s3, chr, c_start, c_end)
s1_cand <- filter_region(s1, chr, c_start, c_end)
cons_cand <- filter_region(cons, chr, c_start, c_end)

n_s2 <- if (!is.null(s2_cand)) nrow(s2_cand) else 0
n_s3 <- if (!is.null(s3_cand)) nrow(s3_cand) else 0
message("[STEP40] Snake 2 windows: ", n_s2, ", Snake 3 windows: ", n_s3)

# =============================================================================
# STEP 2: Build coherence profile across the candidate
# =============================================================================

# Merge all snake tracks by window position
profile_windows <- NULL

if (n_s2 > 0) {
  profile_windows <- s2_cand[, .(
    chrom, global_window_id, start_bp, end_bp,
    s2_hard_nn = mean_hard_nn,
    s2_fuzzy_nn = mean_fuzzy_nn,
    s2_middle = mean_middle_stability,
    s2_tail = mean_tail_match,
    s2_combined = combined_score,
    s2_status = snake2_status
  )]
}

if (n_s3 > 0 && !is.null(profile_windows)) {
  s3_slim <- s3_cand[, .(global_window_id,
    s3_status = snake3_status)]
  # Add GHSL columns if available
  ghsl_cols <- intersect(c("ghsl_mean", "ghsl_median", "ghsl_total_weighted"),
                          names(s3_cand))
  if (length(ghsl_cols) > 0) {
    s3_slim <- cbind(s3_slim, s3_cand[, ..ghsl_cols])
  }
  profile_windows <- merge(profile_windows, s3_slim,
                            by = "global_window_id", all.x = TRUE)
} else if (n_s3 > 0) {
  profile_windows <- s3_cand[, .(chrom, global_window_id, start_bp, end_bp,
    s3_status = snake3_status)]
}

if (is.null(profile_windows) || nrow(profile_windows) == 0) {
  message("[STEP40] No snake data for this candidate — skipping")
  # Write empty verdict
  fwrite(data.table(candidate_id = cid, verdict = "NO_DATA"),
         file.path(cand_dir, "candidate_coherence_verdict.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

setorder(profile_windows, start_bp)
n_windows <- nrow(profile_windows)

# Add consensus info
if (!is.null(cons_cand) && nrow(cons_cand) > 0) {
  cons_slim <- cons_cand[, .(global_window_id,
    n_snakes = n_snakes_supporting, support_sig = support_signature)]
  profile_windows <- merge(profile_windows, cons_slim,
                            by = "global_window_id", all.x = TRUE)
}

message("[STEP40] Profile windows: ", n_windows)

# =============================================================================
# STEP 3: Detect internal coherence breaks
# =============================================================================

# A "break" is a stretch of windows where Snake 2 coherence drops below
# a threshold. We look for transitions from coherent → incoherent → coherent.

# A "break" is a stretch where multiple signals agree that coherence collapsed.
# NOT just s2_combined alone — that would be brittle.
# A window is "incoherent" if s2_combined is low AND at least one of:
#   - middle-band stability dropped below 0.25
#   - hard and fuzzy NN diverge strongly (hard low, fuzzy also low)
#   - consensus support dropped to 0 or 1 snake
# This multi-signal requirement prevents false breaks from single-score noise.

detect_breaks <- function(profile) {
  if (!"s2_combined" %in% names(profile)) {
    return(list(breaks = data.table(), sub_candidates = data.table()))
  }

  scores <- profile$s2_combined
  scores[is.na(scores)] <- 0

  # Multi-signal incoherence: combined is low AND at least one confirming signal
  mid_drop <- if ("s2_middle" %in% names(profile)) {
    !is.na(profile$s2_middle) & profile$s2_middle < 0.25
  } else rep(FALSE, nrow(profile))

  nn_both_low <- if ("s2_hard_nn" %in% names(profile) && "s2_fuzzy_nn" %in% names(profile)) {
    (!is.na(profile$s2_hard_nn) & profile$s2_hard_nn < 0.30) &
    (!is.na(profile$s2_fuzzy_nn) & profile$s2_fuzzy_nn < 0.30)
  } else rep(FALSE, nrow(profile))

  weak_consensus <- if ("n_snakes" %in% names(profile)) {
    !is.na(profile$n_snakes) & profile$n_snakes <= 1
  } else rep(FALSE, nrow(profile))

  combined_low <- scores < BREAK_THRESH
  has_confirming_signal <- mid_drop | nn_both_low | weak_consensus

  # Coherent = NOT (combined_low AND confirmed)
  profile[, coherent := !(combined_low & has_confirming_signal)]

  # Find runs of coherent/incoherent
  runs <- rle(profile$coherent)
  n_runs <- length(runs$lengths)

  breaks <- list()
  subs <- list()
  sub_id <- 0L

  pos <- 0L
  for (ri in seq_len(n_runs)) {
    run_start <- pos + 1L
    run_end <- pos + runs$lengths[ri]
    pos <- run_end

    if (!runs$values[ri] && runs$lengths[ri] >= MIN_BREAK_RUN) {
      # This is a break zone
      breaks[[length(breaks) + 1]] <- data.table(
        break_id = length(breaks) + 1L,
        break_start_bp = profile$start_bp[run_start],
        break_end_bp = profile$end_bp[run_end],
        break_n_windows = runs$lengths[ri],
        mean_s2_score = round(mean(scores[run_start:run_end]), 4),
        break_start_idx = run_start,
        break_end_idx = run_end
      )
    }

    if (runs$values[ri] && runs$lengths[ri] >= MIN_COHERENT_RUN) {
      # This is a coherent sub-candidate
      sub_id <- sub_id + 1L

      # Get middle stability stats for this sub
      mid_vals <- profile$s2_middle[run_start:run_end]
      mid_vals <- mid_vals[is.finite(mid_vals)]

      subs[[length(subs) + 1]] <- data.table(
        sub_candidate_id = sub_id,
        parent_candidate_id = cid,
        chrom = chr,
        start_bp = profile$start_bp[run_start],
        end_bp = profile$end_bp[run_end],
        n_windows = runs$lengths[ri],
        mean_s2_combined = round(mean(scores[run_start:run_end]), 4),
        mean_middle_stability = if (length(mid_vals) > 0) round(mean(mid_vals), 4) else NA_real_,
        min_s2_combined = round(min(scores[run_start:run_end]), 4),
        start_idx = run_start,
        end_idx = run_end
      )
    }
  }

  list(
    breaks = if (length(breaks) > 0) rbindlist(breaks) else data.table(),
    sub_candidates = if (length(subs) > 0) rbindlist(subs) else data.table()
  )
}

result <- detect_breaks(profile_windows)
n_breaks <- nrow(result$breaks)
n_subs <- nrow(result$sub_candidates)

message("[STEP40] Breaks detected: ", n_breaks)
message("[STEP40] Sub-candidates proposed: ", n_subs)

# =============================================================================
# STEP 4: Compare sample roles across sub-candidates (if applicable)
# =============================================================================

role_comparison <- data.table()

if (n_subs >= 2 && !is.null(s1_cand) && nrow(s1_cand) > 0) {
  # For each sub-candidate, find which Snake 1 cores overlap
  # and compare their scale_tier membership
  for (si in seq_len(n_subs)) {
    sub <- result$sub_candidates[si]
    s1_sub <- s1_cand[start_bp >= sub$start_bp & end_bp <= sub$end_bp]

    # Count by core family
    n_1S <- sum(s1_sub$scale_tier == "1S", na.rm = TRUE)
    n_1M <- sum(s1_sub$scale_tier == "1M", na.rm = TRUE)
    n_1L <- sum(s1_sub$scale_tier == "1L", na.rm = TRUE)

    # Count by Snake 2 status
    prof_sub <- profile_windows[start_bp >= sub$start_bp & end_bp <= sub$end_bp]
    n_pass <- sum(prof_sub$s2_status == "PASS", na.rm = TRUE)
    n_weak <- sum(prof_sub$s2_status == "WEAK", na.rm = TRUE)

    role_comparison <- rbind(role_comparison, data.table(
      sub_candidate_id = sub$sub_candidate_id,
      start_bp = sub$start_bp, end_bp = sub$end_bp,
      n_windows = sub$n_windows,
      s1_cores_1S = n_1S, s1_cores_1M = n_1M, s1_cores_1L = n_1L,
      s2_pass = n_pass, s2_weak = n_weak,
      mean_middle_stability = sub$mean_middle_stability
    ), fill = TRUE)
  }

  # ── Cross-sub-candidate middle-band comparison ───────────────────────
  # The key question: do the same samples fill the middle band in
  # different sub-candidates, or do different systems take over?
  #
  # If the same middle samples appear across sub-candidates → one system
  # with a break zone in between (could still be one inversion).
  # If different middle samples → truly different systems were merged.
  #
  # Also check: is the middle internally structured (sub-branches)?
  # A middle containing sys1/sys2 AND sys1/sys3 pairings would show
  # internal substructure that a single Jaccard misses.

  mds_rds <- paste0(step10_prefix, ".mds.rds")
  mds_obj_local <- if (file.exists(mds_rds)) readRDS(mds_rds) else NULL
  middle_comparison <- data.table()

  if (n_subs >= 2 && !is.null(mds_obj_local) && chr %in% names(mds_obj_local$per_chr)) {
    chr_obj <- mds_obj_local$per_chr[[chr]]
    dt_full <- as.data.table(chr_obj$out_dt)[order(start_bp)]
    dmat <- chr_obj$dmat

    # For each sub-candidate, identify middle-band samples via quick PC1 split
    sub_middles <- list()
    for (si in seq_len(n_subs)) {
      sub <- result$sub_candidates[si]
      sub_mask <- dt_full$start_bp >= sub$start_bp & dt_full$end_bp <= sub$end_bp
      sub_idx <- which(sub_mask)
      if (length(sub_idx) < 3 || nrow(dmat) < max(sub_idx)) next

      # Average distance across sub-candidate windows → rough sample structure
      sub_dmat <- dmat[sub_idx, sub_idx, drop = FALSE]
      avg_dist <- colMeans(sub_dmat, na.rm = TRUE)

      # Quick MDS dim 1 for middle identification
      mds1 <- tryCatch(cmdscale(as.dist(sub_dmat), k = 1), error = function(e) NULL)
      if (is.null(mds1)) next

      pc1 <- as.numeric(mds1)
      km <- tryCatch(kmeans(pc1, centers = 3, nstart = 5), error = function(e) NULL)
      if (is.null(km)) next

      co <- order(km$centers[,1])
      middle_cluster <- co[2]
      sub_middles[[si]] <- which(km$cluster == middle_cluster)
    }

    # Compare middle sets between all pairs of sub-candidates
    for (a in seq_len(n_subs - 1)) {
      for (b in (a + 1):n_subs) {
        if (is.null(sub_middles[[a]]) || is.null(sub_middles[[b]])) next
        ma <- sub_middles[[a]]; mb <- sub_middles[[b]]
        jacc <- length(intersect(ma, mb)) / max(1, length(union(ma, mb)))
        middle_comparison <- rbind(middle_comparison, data.table(
          sub_A = a, sub_B = b,
          n_middle_A = length(ma), n_middle_B = length(mb),
          middle_overlap = length(intersect(ma, mb)),
          middle_jaccard = round(jacc, 4),
          interpretation = if (jacc >= 0.50) "same_middle_system"
                          else if (jacc >= 0.20) "partially_shared_middle"
                          else "different_middle_systems"
        ))
      }
    }
  }

  if (nrow(middle_comparison) > 0) {
    fwrite(middle_comparison, file.path(cand_dir, "candidate_middle_comparison.tsv"), sep = "\t")
    message("[STEP40] Middle-band comparison: ", nrow(middle_comparison), " pairs")
    for (mi in seq_len(nrow(middle_comparison))) {
      message("  Sub ", middle_comparison$sub_A[mi], " vs ",
              middle_comparison$sub_B[mi], ": Jaccard=",
              middle_comparison$middle_jaccard[mi], " → ",
              middle_comparison$interpretation[mi])
    }
  }
}

# =============================================================================
# STEP 5: Verdict
# =============================================================================

if (n_breaks == 0) {
  verdict <- "COHERENT"
  verdict_detail <- "Snake 2/3 show one continuous system across the full candidate"
  recommendation <- "No splitting needed — candidate appears to be one structural regime"
} else if (n_subs >= 2) {
  # Check whether breaks are internal or edge-only
  # Edge = within EDGE_FRAC of candidate start or end
  edge_zone_bp <- (c_end - c_start) * EDGE_FRAC
  all_breaks_at_edges <- all(
    result$breaks$break_start_bp <= c_start + edge_zone_bp |
    result$breaks$break_end_bp >= c_end - edge_zone_bp
  )
  any_internal_break <- any(
    result$breaks$break_start_bp > c_start + edge_zone_bp &
    result$breaks$break_end_bp < c_end - edge_zone_bp
  )

  if (any_internal_break) {
    # Enrich with middle-band comparison if available
    mid_verdict_extra <- ""
    if (nrow(middle_comparison) > 0) {
      n_diff <- sum(middle_comparison$interpretation == "different_middle_systems")
      n_same <- sum(middle_comparison$interpretation == "same_middle_system")
      if (n_diff > 0) {
        mid_verdict_extra <- paste0(" Middle-band analysis: ", n_diff,
          " sub-candidate pair(s) have different middle systems — ",
          "likely multiple inversion-like arrangements merged.")
      } else if (n_same > 0) {
        mid_verdict_extra <- paste0(" Middle-band analysis: same samples in middle ",
          "across sub-candidates — may be one system with a gap, not truly composite.")
      }
    }
    verdict <- "COMPOSITE"
    verdict_detail <- paste0(n_breaks, " internal break(s) detected, ", n_subs,
                             " sub-candidates proposed.", mid_verdict_extra)
    recommendation <- paste0("Consider splitting into ", n_subs,
                             " window-defined sub-candidates. ",
                             "Then compare sample roles across sub-candidates.")
  } else {
    verdict <- "WEAK_EDGES"
    verdict_detail <- paste0(n_breaks, " break(s) at candidate edges only — ",
                             "core remains coherent, boundaries are weak")
    recommendation <- paste0("Core is coherent but boundaries are uncertain. ",
                             "Consider trimming edges rather than splitting.")
  }
} else {
  # Breaks found but not enough coherent runs for sub-candidates
  verdict <- "WEAK_EDGES"
  verdict_detail <- paste0(n_breaks, " break(s) but no clear sub-candidates above ",
                           MIN_COHERENT_RUN, " windows")
  recommendation <- "Candidate may have noisy edges but no convincing internal split"
}

verdict_dt <- data.table(
  candidate_id = cid, chrom = chr,
  start_bp = c_start, end_bp = c_end,
  n_windows = n_windows,
  n_breaks = n_breaks,
  n_sub_candidates = n_subs,
  verdict = verdict,
  verdict_detail = verdict_detail,
  recommendation = recommendation,
  mean_s2_combined = if ("s2_combined" %in% names(profile_windows))
    round(mean(profile_windows$s2_combined, na.rm = TRUE), 4) else NA_real_,
  mean_middle_stability = if ("s2_middle" %in% names(profile_windows))
    round(mean(profile_windows$s2_middle, na.rm = TRUE), 4) else NA_real_
)

message("[STEP40] Verdict: ", verdict, " — ", verdict_detail)

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

fwrite(profile_windows, file.path(cand_dir, "candidate_coherence_profile.tsv.gz"), sep = "\t")
fwrite(result$breaks, file.path(cand_dir, "candidate_internal_breaks.tsv"), sep = "\t")
fwrite(result$sub_candidates, file.path(cand_dir, "candidate_sub_candidates.tsv"), sep = "\t")
fwrite(role_comparison, file.path(cand_dir, "candidate_role_comparison.tsv"), sep = "\t")
fwrite(verdict_dt, file.path(cand_dir, "candidate_coherence_verdict.tsv"), sep = "\t")

# =============================================================================
# DIAGNOSTIC FIGURE
# =============================================================================

if ("s2_combined" %in% names(profile_windows)) {
  message("[STEP40] Building diagnostic figure...")

  fmt_mb <- function(bp) sprintf("%.2f", bp / 1e6)

  # Panel A: Snake 2 combined score profile
  pA <- ggplot(profile_windows, aes(x = start_bp / 1e6, y = s2_combined)) +
    geom_line(color = "#2563eb", linewidth = 0.5) +
    geom_point(aes(color = s2_status), size = 1.2) +
    geom_hline(yintercept = BREAK_THRESH, linetype = "dashed", color = "red", linewidth = 0.3) +
    scale_color_manual(values = c("PASS" = "#16a34a", "WEAK" = "#f59e0b", "FAIL" = "#e5e7eb"),
                       name = "Status") +
    labs(x = NULL, y = "S2 Combined",
         title = paste0("Candidate ", cid, " — Internal Coherence (", verdict, ")"),
         subtitle = paste0(chr, ":", fmt_mb(c_start), "–", fmt_mb(c_end), " Mb | ",
                           verdict_detail)) +
    theme_minimal(base_size = 9)

  # Add break zones
  if (n_breaks > 0) {
    for (bi in seq_len(n_breaks)) {
      pA <- pA + annotate("rect",
        xmin = result$breaks$break_start_bp[bi] / 1e6,
        xmax = result$breaks$break_end_bp[bi] / 1e6,
        ymin = -Inf, ymax = Inf, fill = "#fca5a5", alpha = 0.3)
    }
  }

  # Panel B: Middle-band stability
  pB <- ggplot(profile_windows, aes(x = start_bp / 1e6, y = s2_middle)) +
    geom_line(color = "#7c3aed", linewidth = 0.5) +
    geom_point(size = 0.8, color = "#7c3aed") +
    labs(x = NULL, y = "Middle Stability",
         subtitle = "Middle-band sample preservation between adjacent windows") +
    theme_minimal(base_size = 9)

  # Panel C: Hard vs fuzzy NN
  pC <- ggplot(profile_windows) +
    geom_line(aes(x = start_bp / 1e6, y = s2_hard_nn), color = "#1d4ed8", linewidth = 0.4) +
    geom_line(aes(x = start_bp / 1e6, y = s2_fuzzy_nn), color = "#f59e0b", linewidth = 0.4) +
    labs(x = NULL, y = "NN Score",
         subtitle = "Blue = hard kNN preservation, Amber = fuzzy kNN preservation") +
    theme_minimal(base_size = 9)

  # Panel D: Consensus support
  pD <- NULL
  if ("n_snakes" %in% names(profile_windows)) {
    pD <- ggplot(profile_windows, aes(x = start_bp / 1e6, y = n_snakes,
                                       fill = factor(n_snakes))) +
      geom_col(width = diff(range(profile_windows$start_bp)) / n_windows / 1e6 * 0.9) +
      scale_fill_manual(values = c("0" = "#f3f4f6", "1" = "#fde68a",
                                    "2" = "#93c5fd", "3" = "#86efac"),
                        name = "Snakes") +
      labs(x = paste0(chr, " (Mb)"), y = "# Snakes",
           subtitle = "Number of snakes supporting each window") +
      theme_minimal(base_size = 9)
  }

  # Panel E: Sub-candidate boundaries
  pE <- NULL
  if (n_subs >= 2) {
    sub_dt <- result$sub_candidates
    pE <- ggplot(sub_dt, aes(xmin = start_bp / 1e6, xmax = end_bp / 1e6,
                              ymin = 0, ymax = 1,
                              fill = factor(sub_candidate_id))) +
      geom_rect(alpha = 0.6, color = "#333", linewidth = 0.3) +
      scale_fill_brewer(palette = "Set2", name = "Sub-candidate") +
      labs(x = paste0(chr, " (Mb)"), y = NULL,
           subtitle = paste0("Proposed sub-candidates (", n_subs, " blocks)")) +
      theme_minimal(base_size = 9) +
      theme(axis.text.y = element_blank())
  }

  # Assemble
  panels <- list(pA, pB, pC)
  if (!is.null(pD)) panels <- c(panels, list(pD))
  if (!is.null(pE)) panels <- c(panels, list(pE))

  n_panels <- length(panels)
  fig_height <- 3 * n_panels

  fig_path <- file.path(cand_dir, "candidate_coherence_diagnostic.pdf")
  pdf(fig_path, width = 12, height = fig_height)
  for (p in panels) print(p)
  dev.off()

  message("[STEP40] Figure: ", fig_path)
}

message("\n[DONE] STEP40 internal coherence diagnostic for candidate ", cid)
message("  Verdict: ", verdict)
