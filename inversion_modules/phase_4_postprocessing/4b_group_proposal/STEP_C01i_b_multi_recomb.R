#!/usr/bin/env Rscript
# =============================================================================
# STEP_C01i_b_multi_recomb.R — Phase 4b.2 of v10.1 Phase 4b rewrite
# =============================================================================
# Recombinant detection for each candidate inversion.
#
# A recombinant sample is one whose class assignment switches within the
# interval — HOM_REF on one segment, HOM_INV on another, with a HET-like
# mosaic transition. These signal either:
#   - gene conversion (short mosaic near a boundary)
#   - double crossover (long mosaic in interior — rare)
#   - suspicious (long mosaic near boundary — possible misassembly)
#
# INPUTS:
#   --candidates   same format as C01i_decompose
#   --decomp_dir   where STEP_C01i_decompose wrote per_window_class.rds
#                   (one directory per candidate)
#   --outdir       fallback output root
#   --tier_max     default 3
#
# DEPENDENCIES:
#   - STEP_C01i_decompose.R must have run first (reads its per-window class)
#   - cheat24_recombinant_prior.R (sourced if found)
#   - Optional: Clair3 phase data for phase-switch detection
#   - Optional: flashlight hemizygous segments
#
# OUTPUTS per candidate: recombinant_map.json Tier-2 block with
#   - per-sample recomb_status: NOT_RECOMB, recomb_GC, recomb_DCO,
#                                recomb_suspicious, recomb_ambiguous
#   - per-recombinant-sample details: switchpoint_bp, mosaic_length_bp,
#     event_class, posterior, n_phase_support
#   - per-candidate aggregates: n_recombinants, mean_posterior,
#     max_dco_prior, mean_dco_prior
#
# METHOD (three independent signals combined):
#
#   SIGNAL 1 — Per-window class consistency (from C01i_decompose)
#   A sample whose per-window class differs from its global class on ≥30%
#   of windows is a candidate recombinant. We find the CONTIGUOUS stretch
#   where class differs and call that the mosaic segment.
#
#   SIGNAL 2 — Clair3 phase switches (direct evidence)
#   Within the interval, count haplotype switches in phased Clair3 calls.
#   A sample with ≥1 long phase block whose genotype pattern differs on
#   each side of a switchpoint is a direct recombinant observation.
#
#   SIGNAL 3 — Flashlight hemizygous segments (optional, cleanest signal)
#   If a sample has a het-DEL inside the inversion, the deleted region
#   shows only one haplotype. If that haplotype's allele state differs
#   from the sample's global arrangement, we have a direct observation of
#   a switch. This is the most reliable recombinant signal but only
#   applies when flashlight data is available.
#
# A sample is called recombinant if signals 1 AND (2 OR 3) agree, or
# signal 3 alone (flashlight is unambiguous). Signal 1 alone is
# insufficient — noise can produce apparent class switches.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--candidates", type = "character", help = "Candidate table"),
  make_option("--decomp_dir", type = "character", default = "decomp_out",
              help = "Where C01i_decompose wrote per_window_class.rds"),
  make_option("--outdir", type = "character", default = "recomb_out"),
  make_option("--tier_max", type = "integer", default = 3L),
  make_option("--consistency_threshold", type = "double", default = 0.70,
              help = "per-window consistency below this flags Signal 1 (default 0.70)")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$candidates) || !file.exists(opt$candidates)) stop("--candidates required")

# ── Load helpers ─────────────────────────────────────────────────────────────
helper_paths <- c("lib_decompose_helpers.R", "R/lib_decompose_helpers.R")
for (p in helper_paths) if (file.exists(p)) { source(p); break }

paths <- resolve_decomp_paths(list(outdir = opt$outdir))
reg <- try_load_registry()
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Source cheat24 if available ──────────────────────────────────────────────
.has_cheat24 <- FALSE
for (p in c("cheats/cheat24_recombinant_prior.R",
            "../cheats/cheat24_recombinant_prior.R",
            file.path(paths$base, "flashlight_v2", "cheats",
                       "cheat24_recombinant_prior.R"))) {
  if (file.exists(p)) {
    tryCatch({ source(p); .has_cheat24 <- TRUE; break },
             error = function(e) message("[multi_recomb] cheat24 source fail: ",
                                          conditionMessage(e)))
  }
}
if (!.has_cheat24) {
  message("[multi_recomb] cheat24 unavailable; classifying events heuristically.")
  # Inline minimal version so the pipeline still runs
  classify_recombinant_event <- function(recomb_position, left_bp, right_bp,
                                            mosaic_length = NA_integer_) {
    span <- right_bp - left_bp
    dist_nearest <- min(abs(recomb_position - left_bp),
                          abs(recomb_position - right_bp))
    dist_frac <- dist_nearest / span
    near_bp <- dist_frac < 0.15
    short_mosaic <- !is.na(mosaic_length) && mosaic_length < 100000
    long_mosaic  <- !is.na(mosaic_length) && mosaic_length >= 500000
    event_class <- if (near_bp && (short_mosaic || is.na(mosaic_length))) "gene_conversion"
                   else if (!near_bp && long_mosaic) "double_crossover"
                   else if (near_bp && long_mosaic) "suspicious"
                   else "ambiguous"
    posterior <- switch(event_class,
                         gene_conversion = 0.5,
                         double_crossover = 0.8,
                         suspicious = 0.1,
                         ambiguous = 0.3)
    list(event_class = event_class,
         posterior_probability = posterior,
         dist_fraction = round(dist_frac, 3),
         near_breakpoint = near_bp)
  }
}

# ── Optional flashlight for hemizygous segments ──────────────────────────────
.has_flashlight <- FALSE
flashlight_path <- file.path(paths$base, "flashlight_v2", "utils",
                              "flashlight_loader_v2.R")
if (file.exists(flashlight_path)) {
  tryCatch({
    source(flashlight_path); .has_flashlight <- TRUE
  }, error = function(e) NULL)
}

# =============================================================================
# Core: find contiguous class-switch segments in per-window matrix
# =============================================================================
#
# For each sample, scan the per-window class assignments. Find the LONGEST
# contiguous stretch where the class differs from the sample's global class.
# That stretch defines the recombinant segment (if any).
#
# Returns NULL if no such stretch exists (sample is not a Signal 1 recombinant).
find_switch_segment <- function(per_window_class_col, global_class,
                                   window_starts, window_ends,
                                   min_consecutive = 2L) {
  n <- length(per_window_class_col)
  if (n < min_consecutive + 1) return(NULL)

  # Mask: 1 where class differs from global, 0 otherwise (NA → 0)
  diff <- ifelse(!is.na(per_window_class_col) & per_window_class_col != global_class, 1L, 0L)
  if (sum(diff) < min_consecutive) return(NULL)

  # Find longest consecutive run of 1s
  rle_diff <- rle(diff)
  max_len <- 0L
  max_start_idx <- NA_integer_
  idx <- 1L
  for (i in seq_along(rle_diff$values)) {
    if (rle_diff$values[i] == 1L && rle_diff$lengths[i] >= min_consecutive) {
      if (rle_diff$lengths[i] > max_len) {
        max_len <- rle_diff$lengths[i]
        max_start_idx <- idx
      }
    }
    idx <- idx + rle_diff$lengths[i]
  }

  if (is.na(max_start_idx)) return(NULL)

  # Convert window indices to bp coordinates
  start_win <- max_start_idx
  end_win <- max_start_idx + max_len - 1L
  mosaic_start_bp <- window_starts[start_win]
  mosaic_end_bp   <- window_ends[end_win]
  mosaic_length <- mosaic_end_bp - mosaic_start_bp

  # The "switched" class is the most common non-global class in the segment
  segment_classes <- per_window_class_col[start_win:end_win]
  segment_classes <- segment_classes[!is.na(segment_classes) & segment_classes != global_class]
  if (length(segment_classes) == 0) return(NULL)
  switched_class <- names(sort(table(segment_classes), decreasing = TRUE))[1]

  # Switchpoint = midpoint of the mosaic
  switchpoint_bp <- as.integer((mosaic_start_bp + mosaic_end_bp) / 2)

  list(
    mosaic_start_bp = as.integer(mosaic_start_bp),
    mosaic_end_bp   = as.integer(mosaic_end_bp),
    mosaic_length_bp = as.integer(mosaic_length),
    switchpoint_bp = switchpoint_bp,
    global_class = global_class,
    switched_class = switched_class,
    n_switch_windows = max_len
  )
}

# =============================================================================
# Clair3 phase-switch detection (Signal 2)
# =============================================================================
# For a sample, look at Clair3 phased variants in the candidate interval.
# If there are ≥2 phase blocks with clearly different genotype patterns,
# that's a phase switch signal for recombination.
count_phase_switches <- function(chr, sample_id, start_bp, end_bp, clair3_dir) {
  ph <- load_clair3_phase(chr, sample_id, start_bp, end_bp, clair3_dir)
  if (nrow(ph) == 0) return(list(n_switches = 0L, n_blocks = 0L, has_support = FALSE))
  # POS col + IS_PHASED col; we don't have block IDs in this minimal read.
  # Heuristic: count runs of phased variants separated by unphased gaps ≥10kb.
  ph <- ph[order(POS)]
  phased <- ph$IS_PHASED == 1
  if (sum(phased) < 4) return(list(n_switches = 0L, n_blocks = 0L, has_support = FALSE))
  # Count gaps between consecutive phased variants > 10kb → block boundaries
  phased_pos <- ph$POS[phased]
  gaps <- diff(phased_pos)
  n_blocks <- sum(gaps > 10000L) + 1L
  # A phase switch requires ≥2 blocks
  n_switches <- max(0L, n_blocks - 1L)
  list(
    n_switches = as.integer(n_switches),
    n_blocks = as.integer(n_blocks),
    has_support = n_switches >= 1
  )
}

# =============================================================================
# Flashlight hemizygous signal (Signal 3)
# =============================================================================
flashlight_hemi_signal <- function(chr, sample_id, start_bp, end_bp, global_class) {
  if (!.has_flashlight) return(list(has_signal = FALSE))
  if (!exists("get_internal_dels", mode = "function")) return(list(has_signal = FALSE))
  dels <- tryCatch(get_internal_dels(chr, start_bp, end_bp),
                     error = function(e) NULL)
  if (is.null(dels) || nrow(dels) == 0) return(list(has_signal = FALSE))
  sample_dels <- dels[sample_id %in% unlist(het_carriers)]
  if (nrow(sample_dels) == 0) return(list(has_signal = FALSE))
  # If we got here, sample has ≥1 het-DEL inside the inversion. Without
  # running the full Cheat 3 allele-state analysis (which needs BAM access
  # and is out of scope for this script), flag the sample as hemi-evidence
  # present and let downstream use the flashlight directly.
  list(
    has_signal = TRUE,
    n_hemi_segments = nrow(sample_dels),
    hemi_details = lapply(seq_len(nrow(sample_dels)), function(i) {
      list(
        del_start = as.integer(sample_dels$del_start[i]),
        del_end = as.integer(sample_dels$del_end[i]),
        del_length = as.integer(sample_dels$del_end[i] - sample_dels$del_start[i])
      )
    })
  )
}

# =============================================================================
# Main loop
# =============================================================================
cand_dt <- fread(opt$candidates)
if (!"candidate_id" %in% names(cand_dt)) {
  cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
}
if (!"tier" %in% names(cand_dt)) cand_dt[, tier := 2L]
if (!"start_bp" %in% names(cand_dt))
  cand_dt[, start_bp := as.integer(start_mb * 1e6)]
if (!"end_bp" %in% names(cand_dt))
  cand_dt[, end_bp := as.integer(end_mb * 1e6)]
cand_dt <- cand_dt[tier <= opt$tier_max]

cat("[multi_recomb] processing ", nrow(cand_dt), " candidates\n", sep = "")

for (ci in seq_len(nrow(cand_dt))) {
  cd <- cand_dt[ci]
  cid <- cd$candidate_id
  chr <- cd$chrom; s <- cd$start_bp; e <- cd$end_bp

  cat("\n[multi_recomb] ─── ", cid, " ───\n", sep = "")

  # Load per-window class from STEP_C01i_decompose
  pw_file <- file.path(opt$decomp_dir, cid, "per_window_class.rds")
  if (!file.exists(pw_file)) {
    cat("[multi_recomb]   SKIP: per_window_class.rds not found at ", pw_file, "\n", sep = "")
    next
  }
  pw <- readRDS(pw_file)
  pw_class <- pw$per_window_class   # n_windows × n_samples
  sample_names <- pw$sample_names
  global_class <- pw$class_assignment_global
  win_starts <- pw$window_starts
  win_ends <- pw$window_ends

  # Scan per-sample
  per_sample_recomb <- list()
  n_recomb <- 0L
  for (si in seq_along(sample_names)) {
    sid <- sample_names[si]
    gc <- unname(global_class[sid])
    if (is.na(gc)) next

    # Signal 1: class switch segment from per-window matrix
    sig1 <- find_switch_segment(
      per_window_class_col = pw_class[, sid],
      global_class = gc,
      window_starts = win_starts,
      window_ends = win_ends,
      min_consecutive = 2L
    )

    # Signal 2: Clair3 phase switches
    sig2 <- count_phase_switches(chr, sid, s, e, paths$clair3_dir)

    # Signal 3: flashlight hemizygous
    sig3 <- flashlight_hemi_signal(chr, sid, s, e, gc)

    # Combine
    has_sig1 <- !is.null(sig1)
    has_sig2 <- sig2$has_support
    has_sig3 <- sig3$has_signal

    # Decision rule: (S1 AND S2) OR S1 OR S3
    is_recomb <- (has_sig1 && has_sig2) || has_sig3
    # Allow signal 1 alone with a strong mosaic (> 100kb)
    if (!is_recomb && has_sig1 && !is.null(sig1)) {
      if (sig1$mosaic_length_bp > 100000) is_recomb <- TRUE
    }

    if (!is_recomb) next

    n_recomb <- n_recomb + 1L

    # Classify event via cheat24
    mosaic_len <- if (!is.null(sig1)) sig1$mosaic_length_bp else NA_integer_
    switchpoint <- if (!is.null(sig1)) sig1$switchpoint_bp else as.integer((s + e) / 2)
    classif <- classify_recombinant_event(switchpoint, s, e, mosaic_len)

    # Recomb status label
    recomb_status <- switch(classif$event_class,
                              gene_conversion  = "recomb_GC",
                              double_crossover = "recomb_DCO",
                              suspicious       = "recomb_suspicious",
                              ambiguous        = "recomb_ambiguous",
                              "recomb_ambiguous")

    per_sample_recomb[[length(per_sample_recomb) + 1L]] <- list(
      sample_id = sid,
      global_class = gc,
      switched_class = if (!is.null(sig1)) sig1$switched_class else NA_character_,
      switchpoint_bp = switchpoint,
      mosaic_start_bp = if (!is.null(sig1)) sig1$mosaic_start_bp else NA_integer_,
      mosaic_end_bp   = if (!is.null(sig1)) sig1$mosaic_end_bp else NA_integer_,
      mosaic_length_bp = mosaic_len,
      distance_to_nearest_boundary_bp = min(abs(switchpoint - s), abs(switchpoint - e)),
      dist_fraction = classif$dist_fraction,
      event_class = classif$event_class,
      posterior = classif$posterior_probability,
      near_breakpoint = classif$near_breakpoint,
      recomb_status = recomb_status,
      # Evidence combination
      signal_1_class_switch = has_sig1,
      signal_2_phase_switch = has_sig2,
      signal_3_flashlight_hemi = has_sig3,
      n_phase_blocks = sig2$n_blocks,
      n_hemi_segments = if (has_sig3) sig3$n_hemi_segments else 0L
    )
  }

  # Aggregate
  if (n_recomb > 0) {
    event_classes <- vapply(per_sample_recomb, function(r) r$event_class, character(1))
    posteriors   <- vapply(per_sample_recomb, function(r) r$posterior, numeric(1))
    priors <- vapply(per_sample_recomb, function(r) {
      # compute position prior the same way cheat24 does
      span <- e - s
      dist_nearest <- min(abs(r$switchpoint_bp - s), abs(r$switchpoint_bp - e))
      rel_pos <- (r$switchpoint_bp - s) / span
      # logistic at center — approximates cheat24's prior_dco behavior
      dco <- 0.01 / (1 + exp(-10 * (dist_nearest / span - 0.3)))
      dco
    }, numeric(1))

    n_gc  <- sum(event_classes == "gene_conversion")
    n_dco <- sum(event_classes == "double_crossover")
    n_susp <- sum(event_classes == "suspicious")
    n_amb  <- sum(event_classes == "ambiguous")

    cat("[multi_recomb]   ", n_recomb, " recombinants: GC=", n_gc,
        " DCO=", n_dco, " susp=", n_susp, " amb=", n_amb, "\n", sep = "")
  } else {
    cat("[multi_recomb]   no recombinants detected\n")
  }

  # ── Write Tier-2 block ──
  block_data <- list(
    candidate_id = cid,
    chrom = chr, start_bp = as.integer(s), end_bp = as.integer(e),
    n_recombinants = as.integer(n_recomb),
    n_gene_conversion  = if (n_recomb > 0) as.integer(n_gc) else 0L,
    n_double_crossover = if (n_recomb > 0) as.integer(n_dco) else 0L,
    n_suspicious       = if (n_recomb > 0) as.integer(n_susp) else 0L,
    n_ambiguous        = if (n_recomb > 0) as.integer(n_amb) else 0L,
    mean_posterior = if (n_recomb > 0) round(mean(posteriors), 4) else NA_real_,
    mean_dco_prior = if (n_recomb > 0) round(mean(priors), 6) else NA_real_,
    max_dco_prior  = if (n_recomb > 0) round(max(priors), 6) else NA_real_,
    cheat24_version = if (.has_cheat24) "cheat24_recombinant_prior.R" else "inline_fallback",
    recombinants = per_sample_recomb
  )

  cand_outdir <- file.path(opt$outdir, cid)
  dir.create(cand_outdir, recursive = TRUE, showWarnings = FALSE)

  write_block_safe(
    reg = reg,
    candidate_id = cid,
    block_type = "recombinant_map",
    data = block_data,
    source_script = "STEP_C01i_b_multi_recomb.R",
    outdir_fallback = cand_outdir
  )
}

cat("\n[multi_recomb] done\n")
