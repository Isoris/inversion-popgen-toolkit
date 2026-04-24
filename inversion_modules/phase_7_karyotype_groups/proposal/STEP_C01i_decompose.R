#!/usr/bin/env Rscript
# =============================================================================
# STEP_C01i_decompose.R — Phase 4b.1 of v10.1 Phase 4b rewrite
# =============================================================================
# Per-candidate sample decomposition into HOM_REF / HET / HOM_INV classes.
#
# INPUTS:
#   --candidates   TSV with at least: candidate_id, chrom, start_bp, end_bp, tier
#   --outdir       per-candidate output root (used as fallback if registry
#                   library is absent)
#   --tier_max     maximum tier to process (default: 3)
#
# INDEPENDENT (no --hyp_dir, no --boundary_dir)
# OUTPUTS per candidate: internal_dynamics.json Tier-2 block containing
#   - per-sample class assignment (HOM_REF/HET/HOM_INV)
#   - per-window dosage track
#   - silhouette_score, bic_gap_k3_vs_k2, phase_support
#   - seeded_by_sv_prior flag (if applicable)
#   - discordant_samples (where sv_prior seeds disagree with clustering)
#   - samples_constrained_cheat2 (het-DEL at breakpoint → not HOM_INV)
#
# DOES NOT register groups yet. STEP_C01i_d_seal does that after all three
# phase 4b scripts complete, using combined evidence.
#
# Sources:
#   - lib_decompose_helpers.R for shared utilities
#   - flashlight_v2/patches/patch_C01i_decomposition_sv_prior.R for cheats 1/2/3
#   - inspired by superseded/STEP_C01i_multi_inversion_decomposition.R structure
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ── Parse args ───────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--candidates", type = "character",
              help = "Candidate table TSV from C01d"),
  make_option("--outdir", type = "character", default = "decomp_out",
              help = "Output root (fallback if no registry)"),
  make_option("--tier_max", type = "integer", default = 3L,
              help = "Only process candidates with tier <= this"),
  make_option("--k_decomp", type = "integer", default = 3L,
              help = "Number of classes (default 3)"),
  make_option("--step03_seeds_dir", type = "character", default = NA_character_,
              help = paste0("Directory containing STEP03 per-candidate seed ",
                            "TSVs (04_deconvolution_seeds/<cid>_seeds.tsv). ",
                            "If set, combined with sv_prior seeds via ",
                            "lib_step03_seed_loader.R::combine_tier1_seeds(). ",
                            "If unset AND sv_prior seeds are unavailable, ",
                            "the candidate is skipped with decomp_status=",
                            "'no_seeding'. The unsupervised k-means fallback ",
                            "has been removed (chat-12, per chat-9 design).")),
  make_option("--min_seeds_per_class", type = "integer", default = 2L,
              help = paste0("Minimum samples per (HOM_REF/HET/HOM_INV) class ",
                            "in the combined seed set (default 2)."))
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$candidates) || !file.exists(opt$candidates)) {
  stop("--candidates must be an existing file")
}

# ── Load helpers + registry ──────────────────────────────────────────────────
helper_paths <- c(
  "lib_decompose_helpers.R",
  "R/lib_decompose_helpers.R",
  file.path(dirname(sys.frame(1)$ofile %||% "."), "lib_decompose_helpers.R")
)
helper_found <- FALSE
for (p in helper_paths) {
  if (file.exists(p)) { source(p); helper_found <- TRUE; break }
}
if (!helper_found) {
  stop("[C01i_decompose] cannot find lib_decompose_helpers.R")
}

# Chat-12: load the Tier-1 seed combiner (sv_prior ∪ STEP03, conflict drop,
# no unsupervised fallback). If missing, seed combination falls back to the
# sv_prior-only path inside the main loop.
seed_loader_paths <- c(
  "lib_step03_seed_loader.R",
  "R/lib_step03_seed_loader.R",
  file.path(dirname(sys.frame(1)$ofile %||% "."), "lib_step03_seed_loader.R")
)
.has_seed_loader <- FALSE
for (p in seed_loader_paths) {
  if (file.exists(p)) {
    tryCatch({ source(p); .has_seed_loader <- TRUE; break },
             error = function(e) message("[C01i_decompose] seed loader source failed: ",
                                          conditionMessage(e)))
  }
}
if (!.has_seed_loader) {
  message("[C01i_decompose] lib_step03_seed_loader.R not found — ",
          "STEP03 seed path disabled, sv_prior-only seeding.")
}

paths <- resolve_decomp_paths(list(outdir = opt$outdir))
reg <- try_load_registry()
if (is.null(reg)) {
  message("[C01i_decompose] registry unavailable — writing to ", opt$outdir)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Optional sv_prior loader (for seeded k-means cheats) ──────────────────
.has_sv_prior <- FALSE
sv_prior_path <- file.path(paths$base, "flashlight_v2", "utils",
                              "flashlight_loader_v2.R")
if (file.exists(sv_prior_path)) {
  tryCatch({
    source(sv_prior_path)
    .has_sv_prior <- TRUE
    message("[C01i_decompose] sv_prior loader sourced")
  }, error = function(e) {
    message("[C01i_decompose] sv_prior source failed: ", conditionMessage(e))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# Seeded k-means (Cheat 1) — uses SV genotypes as cluster initialization
# =============================================================================
#
# If the candidate has ≥10 SV-genotyped anchor samples (from DELLY/Manta
# sv_prior), use their genotypes as k-means seeds instead of random init.
# This is much more reliable than unsupervised clustering at 9x coverage.
seeded_kmeans <- function(x, seeds, k = 3L, nstart_fallback = 25L) {
  seed_ref <- intersect(seeds$HOM_REF, names(x))
  seed_het <- intersect(seeds$HET, names(x))
  seed_inv <- intersect(seeds$HOM_INV, names(x))
  n_seeds <- length(seed_ref) + length(seed_het) + length(seed_inv)

  if (n_seeds < 10 || length(seed_ref) < 3 || length(seed_inv) < 2) {
    # Fallback: standard k-means
    km <- kmeans_with_quality(x, k = k, nstart = nstart_fallback)
    if (!is.null(km)) km$seeded <- FALSE
    return(km)
  }

  # Seed centroids from anchor means
  center_ref <- mean(x[seed_ref], na.rm = TRUE)
  center_inv <- mean(x[seed_inv], na.rm = TRUE)
  center_het <- if (length(seed_het) >= 3) {
    mean(x[seed_het], na.rm = TRUE)
  } else {
    (center_ref + center_inv) / 2
  }
  centers <- sort(c(center_ref, center_het, center_inv))

  # k-means with seeded starts
  km <- tryCatch(
    kmeans(x, centers = matrix(centers, ncol = 1), iter.max = 30),
    error = function(e) NULL
  )

  if (is.null(km)) {
    # Manual assignment fallback
    dists <- matrix(0, length(x), k)
    for (ci in seq_len(k)) dists[, ci] <- abs(x - centers[ci])
    cluster <- apply(dists, 1, which.min)
    names(cluster) <- names(x)
    return(list(
      cluster = cluster, centers = matrix(centers, ncol = 1),
      size = as.integer(table(factor(cluster, levels = 1:k))),
      tot_withinss = NA_real_, betweenss = NA_real_,
      silhouette_mean = compute_silhouette_mean(matrix(x, ncol = 1), cluster),
      seeded = TRUE,
      n_seed_ref = length(seed_ref), n_seed_het = length(seed_het),
      n_seed_inv = length(seed_inv)
    ))
  }

  list(
    cluster = km$cluster, centers = km$centers, size = km$size,
    tot_withinss = km$tot.withinss, betweenss = km$betweenss,
    silhouette_mean = compute_silhouette_mean(matrix(x, ncol = 1), km$cluster),
    seeded = TRUE,
    n_seed_ref = length(seed_ref),
    n_seed_het = length(seed_het),
    n_seed_inv = length(seed_inv)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Per-window class track — does the sample's class stay consistent across
# all windows of the candidate? If not, they might be recombinant, but that's
# settled downstream by STEP_C01i_b_multi_recomb. Here we just compute the
# per-window class and consistency score.
compute_per_window_class <- function(pc1_mat, class_assignment_global, centers_global) {
  # pc1_mat: n_windows × n_samples
  # class_assignment_global: named character vector of global class per sample
  # centers_global: the 3 global kmeans centers (sorted low → high)
  stopifnot(!is.null(colnames(pc1_mat)))
  sample_names <- colnames(pc1_mat)
  n_win <- nrow(pc1_mat)

  # For each window, assign each sample to the nearest of the three global
  # centers. Produces an n_windows × n_samples matrix of class labels.
  class_mat <- matrix(NA_character_, nrow = n_win, ncol = ncol(pc1_mat),
                       dimnames = list(NULL, sample_names))
  labels <- c("HOM_REF", "HET", "HOM_INV")
  for (w in seq_len(n_win)) {
    row <- pc1_mat[w, ]
    for (si in seq_along(row)) {
      if (is.na(row[si])) next
      dists <- abs(row[si] - as.numeric(centers_global))
      class_mat[w, si] <- labels[which.min(dists)]
    }
  }

  # Per-sample consistency: fraction of windows matching the global class
  consistency <- vapply(sample_names, function(s) {
    gc <- class_assignment_global[s]
    if (is.na(gc)) return(NA_real_)
    col <- class_mat[, s]
    col <- col[!is.na(col)]
    if (length(col) == 0) return(NA_real_)
    mean(col == gc)
  }, numeric(1))

  list(per_window_class = class_mat, consistency = consistency)
}

# ─────────────────────────────────────────────────────────────────────────────
# Main loop
# =============================================================================

cat("[C01i_decompose] loading candidates...\n")
cand_dt <- fread(opt$candidates)

# Normalize column names
if (!"candidate_id" %in% names(cand_dt)) {
  if ("interval_id" %in% names(cand_dt) && "chrom" %in% names(cand_dt)) {
    cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
  } else stop("[C01i_decompose] candidates table needs candidate_id or chrom+interval_id")
}
if (!"tier" %in% names(cand_dt)) cand_dt[, tier := 2L]
if (!"start_bp" %in% names(cand_dt)) {
  if ("start_mb" %in% names(cand_dt)) cand_dt[, start_bp := as.integer(start_mb * 1e6)]
  else stop("[C01i_decompose] need start_bp or start_mb")
}
if (!"end_bp" %in% names(cand_dt)) {
  if ("end_mb" %in% names(cand_dt)) cand_dt[, end_bp := as.integer(end_mb * 1e6)]
  else stop("[C01i_decompose] need end_bp or end_mb")
}

cand_dt <- cand_dt[tier <= opt$tier_max]
cat("[C01i_decompose] processing ", nrow(cand_dt), " candidates at tier <= ",
    opt$tier_max, "\n", sep = "")

# Per-candidate outputs (for fallback)
cand_outdirs <- character(nrow(cand_dt))
for (ci in seq_len(nrow(cand_dt))) {
  cand_outdirs[ci] <- file.path(opt$outdir, cand_dt$candidate_id[ci])
  dir.create(cand_outdirs[ci], recursive = TRUE, showWarnings = FALSE)
}

# Process each candidate
for (ci in seq_len(nrow(cand_dt))) {
  cd <- cand_dt[ci]
  cid <- cd$candidate_id
  chr <- cd$chrom
  s <- cd$start_bp
  e <- cd$end_bp

  cat("\n[C01i_decompose] ─── ", cid, " (", chr, ":", round(s/1e6, 2), "-",
      round(e/1e6, 2), " Mb) ───\n", sep = "")

  # ── Load PC loadings from precomp ──
  pc <- load_precomp(chr, paths$precomp_dir)
  if (is.null(pc)) {
    cat("[C01i_decompose]   SKIP: no precomp for ", chr, "\n")
    next
  }

  loadings <- extract_pc_loadings(pc, s, e, min_windows = 5L)
  if (is.null(loadings)) {
    cat("[C01i_decompose]   SKIP: insufficient windows or no PC columns\n")
    next
  }

  sample_names <- loadings$sample_names
  n_samples <- length(sample_names)
  n_windows <- loadings$n_windows

  # Build per-sample mean PC1 (1D clustering — matches v9.3.4 behavior)
  mean_pc1 <- colMeans(loadings$pc1_mat, na.rm = TRUE)
  mean_pc2 <- colMeans(loadings$pc2_mat, na.rm = TRUE)
  # Optional: 2D clustering on (PC1, PC2) for more stable results
  # We use 1D (PC1 only) to preserve compatibility with v9.3.4 and because
  # PC1 is where the arrangement signal lives. PC2 mostly captures ancestry.

  # ── Cheat 1: Flashlight seeded init ──
  # The sv_prior_seeds list is gathered first; then combine_tier1_seeds()
  # (chat-12) unions sv_prior ∪ STEP03 seeds, drops samples where the two
  # sources disagree on class, and emits status="no_seeding" if the result
  # can't meet min_seeds_per_class on all three classes. The unsupervised
  # k-means fallback has been REMOVED (chat-12, per chat-9 design): if
  # Tier-1 has no seeding, the candidate is skipped with
  # decomp_status="no_seeding" and no groups are proposed.
  sv_prior_seeds <- list(HOM_REF = character(0), HET = character(0),
                            HOM_INV = character(0))
  cheat2_constrained <- character(0)
  sv_prior_available <- FALSE

  if (.has_sv_prior && exists("load_sv_prior", mode = "function")) {
    fl <- tryCatch(load_sv_prior(chr), error = function(e) NULL)
    if (!is.null(fl)) {
      if (exists("get_sv_anchors", mode = "function")) {
        anchors <- tryCatch(
          get_sv_anchors(chr, s, e, min_confidence = "MEDIUM"),
          error = function(e) NULL
        )
        if (!is.null(anchors) && nrow(anchors) >= 10) {
          sv_prior_seeds$HOM_REF <- unique(anchors[sv_genotype == "HOM_REF"]$sample_id)
          sv_prior_seeds$HET     <- unique(anchors[sv_genotype == "HET"]$sample_id)
          sv_prior_seeds$HOM_INV <- unique(anchors[sv_genotype == "HOM_INV"]$sample_id)
          sv_prior_available <- TRUE
          cat("[C01i_decompose]   sv_prior seeds: REF=",
              length(sv_prior_seeds$HOM_REF),
              " HET=", length(sv_prior_seeds$HET),
              " INV=", length(sv_prior_seeds$HOM_INV), "\n", sep = "")
        }
      }
      # Cheat 2: Het-DEL at breakpoint → not HOM_INV
      if (exists("get_breakpoint_dels", mode = "function")) {
        bp_dels_left  <- tryCatch(get_breakpoint_dels(chr, s, window = 50000L),
                                    error = function(e) list(het_carriers = list()))
        bp_dels_right <- tryCatch(get_breakpoint_dels(chr, e, window = 50000L),
                                    error = function(e) list(het_carriers = list()))
        cheat2_constrained <- unique(c(
          unlist(bp_dels_left$het_carriers),
          unlist(bp_dels_right$het_carriers)
        ))
        cheat2_constrained <- intersect(cheat2_constrained, sample_names)
        if (length(cheat2_constrained) > 0) {
          cat("[C01i_decompose]   cheat2: ", length(cheat2_constrained),
              " samples constrained (het-DEL at breakpoint)\n", sep = "")
        }
      }
    }
  }

  # ── Combined Tier-1 seeding (chat-12) ──
  step03_dir <- if (!is.na(opt$step03_seeds_dir) && nzchar(opt$step03_seeds_dir))
                    opt$step03_seeds_dir else NULL

  seed_result <- if (.has_seed_loader) {
    combine_tier1_seeds(
      chr = chr, cid = cid, start_bp = s, end_bp = e,
      sv_prior_seeds    = if (sv_prior_available) sv_prior_seeds else NULL,
      step03_seeds_dir    = step03_dir,
      min_seeds_per_class = opt$min_seeds_per_class
    )
  } else if (sv_prior_available) {
    # Seed loader missing — degrade gracefully to sv_prior-only
    list(status = "sv_prior_only", reason = NA_character_,
         seeds = sv_prior_seeds, n_agree = 0L, n_conflict = 0L,
         source = "sv_prior_legacy")
  } else {
    list(status = "no_seeding",
         reason = "neither sv_prior nor STEP03 seeds available",
         seeds = NULL, n_agree = 0L, n_conflict = 0L, source = "none")
  }

  if (identical(seed_result$status, "no_seeding")) {
    cat("[C01i_decompose]   SKIP (decomp_status=no_seeding): ",
        seed_result$reason, "\n", sep = "")
    # Emit a minimal Tier-2 block so downstream scripts don't silently
    # miss the candidate — they can read decomp_status and skip accordingly.
    block_data <- list(
      candidate_id   = cid,
      chrom          = chr,
      start_bp       = as.integer(s),
      end_bp         = as.integer(e),
      n_samples      = as.integer(n_samples),
      n_windows      = as.integer(n_windows),
      decomp_status  = "no_seeding",
      decomp_reason  = seed_result$reason,
      seed_source    = seed_result$source,
      n_seed_conflict= as.integer(seed_result$n_conflict %||% 0L)
    )
    write_block_safe(
      reg = reg,
      candidate_id = cid,
      block_type   = "internal_dynamics",
      data         = block_data,
      source_script= "STEP_C01i_decompose.R"
    )
    next
  }

  seeds <- seed_result$seeds
  sv_prior_mode <- seed_result$source  # e.g. "sv_prior+step03", "step03",
                                          # "sv_prior", "sv_prior_legacy"

  # ── K-means clustering (seeded only; unsupervised fallback removed) ──
  km <- seeded_kmeans(mean_pc1, seeds, k = opt$k_decomp)
  if (is.null(km)) {
    cat("[C01i_decompose]   SKIP: seeded_kmeans failed\n")
    next
  }
  if (!isTRUE(km$seeded)) {
    # seeded_kmeans returned the unsupervised fallback (n_seeds<10 or
    # under-represented classes). Per chat-9 design we do NOT ship this —
    # emit no_seeding and continue.
    cat("[C01i_decompose]   SKIP (decomp_status=no_seeding): ",
        "seeds too sparse for seeded_kmeans (n_ref=", length(seeds$HOM_REF),
        " n_het=", length(seeds$HET),
        " n_inv=", length(seeds$HOM_INV), ")\n", sep = "")
    block_data <- list(
      candidate_id   = cid,
      chrom          = chr,
      start_bp       = as.integer(s),
      end_bp         = as.integer(e),
      n_samples      = as.integer(n_samples),
      n_windows      = as.integer(n_windows),
      decomp_status  = "no_seeding",
      decomp_reason  = "seeded_kmeans fell back to unsupervised — not shipping",
      seed_source    = seed_result$source,
      n_seed_conflict= as.integer(seed_result$n_conflict %||% 0L)
    )
    write_block_safe(
      reg = reg,
      candidate_id = cid,
      block_type   = "internal_dynamics",
      data         = block_data,
      source_script= "STEP_C01i_decompose.R"
    )
    next
  }

  # ── Assign canonical class labels (HOM_REF/HET/HOM_INV by PC1 ordering) ──
  cluster_order <- order(as.numeric(km$centers))
  labels_by_cluster <- setNames(c("HOM_REF", "HET", "HOM_INV"), cluster_order)
  class_assignment <- labels_by_cluster[as.character(km$cluster)]
  names(class_assignment) <- sample_names

  # ── Apply Cheat 2 (het-DEL at breakpoint forces HET, not HOM_INV) ──
  n_cheat2_reclassified <- 0L
  if (length(cheat2_constrained) > 0) {
    idx <- sample_names %in% cheat2_constrained & class_assignment == "HOM_INV"
    n_cheat2_reclassified <- sum(idx)
    if (n_cheat2_reclassified > 0) {
      class_assignment[idx] <- "HET"
      cat("[C01i_decompose]   cheat2 reclassified: ",
          n_cheat2_reclassified, " HOM_INV → HET\n", sep = "")
    }
  }

  # ── Flag sv_prior-discordant samples ──
  discordant <- character(0)
  if (sv_prior_mode == "seeded") {
    for (sid in sample_names) {
      sv_gt <- if (sid %in% seeds$HOM_REF) "HOM_REF"
                else if (sid %in% seeds$HET) "HET"
                else if (sid %in% seeds$HOM_INV) "HOM_INV"
                else NA
      if (!is.na(sv_gt) && sv_gt != class_assignment[sid]) {
        discordant <- c(discordant, sid)
      }
    }
    if (length(discordant) > 0) {
      cat("[C01i_decompose]   discordant (SV != clustering): ",
          length(discordant), "\n", sep = "")
    }
  }

  # ── Per-window class track + consistency score ──
  centers_sorted <- sort(as.numeric(km$centers))
  pw <- compute_per_window_class(loadings$pc1_mat, class_assignment, centers_sorted)

  # ── Quality metrics ──
  silhouette <- km$silhouette_mean %||% NA_real_
  bic_gap <- tryCatch(compute_bic_gap(mean_pc1), error = function(e) NA_real_)

  # ── Phase support (how many HET samples have phased Clair3 in region?) ──
  n_phase_support <- 0L
  het_samples <- names(class_assignment)[class_assignment == "HET"]
  for (sid in het_samples) {
    ph <- load_clair3_phase(chr, sid, s, e, paths$clair3_dir)
    if (nrow(ph) > 0 && mean(ph$IS_PHASED) > 0.5) {
      n_phase_support <- n_phase_support + 1L
    }
  }
  phase_concordance <- if (length(het_samples) > 0) {
    n_phase_support / length(het_samples)
  } else NA_real_

  # ── Summary counts ──
  class_counts <- as.list(table(factor(class_assignment,
                                          levels = c("HOM_REF", "HET", "HOM_INV"))))

  cat("[C01i_decompose]   classes: HOM_REF=", class_counts$HOM_REF,
      " HET=", class_counts$HET,
      " HOM_INV=", class_counts$HOM_INV,
      " | silhouette=", round(silhouette, 3),
      " bic_gap=", round(bic_gap, 3), "\n", sep = "")

  # ── Assemble Tier-2 block data ──
  per_sample <- lapply(seq_along(sample_names), function(si) {
    sid <- sample_names[si]
    list(
      sample_id = sid,
      pca_class = unname(class_assignment[sid]),
      mean_pc1  = round(unname(mean_pc1[sid]), 4),
      mean_pc2  = round(unname(mean_pc2[sid]), 4),
      window_consistency = round(unname(pw$consistency[sid]), 3),
      sv_prior_seeded  = sid %in% c(seeds$HOM_REF, seeds$HET, seeds$HOM_INV),
      discordant = sid %in% discordant,
      cheat2_constrained = sid %in% cheat2_constrained
    )
  })

  block_data <- list(
    candidate_id = cid,
    chrom = chr, start_bp = as.integer(s), end_bp = as.integer(e),
    n_samples = as.integer(n_samples),
    n_windows = as.integer(n_windows),
    k_used = as.integer(opt$k_decomp),
    # Chat-12: explicit status field so downstream readers can distinguish
    # "succeeded with seeds" from "skipped with no_seeding".
    decomp_status = "seeded_ok",
    seed_source   = seed_result$source,
    n_seed_agree    = as.integer(seed_result$n_agree    %||% 0L),
    n_seed_conflict = as.integer(seed_result$n_conflict %||% 0L),
    # Class counts (preliminary — final counts come from C01i_d_seal)
    n_HOM_REF_prelim = as.integer(class_counts$HOM_REF %||% 0),
    n_HET_prelim     = as.integer(class_counts$HET %||% 0),
    n_HOM_INV_prelim = as.integer(class_counts$HOM_INV %||% 0),
    # Quality metrics
    silhouette_score = round(silhouette, 4),
    # chat-13 Finding AV: non-gating decomp_quality annotation derived
    # from silhouette. At 9x with limited samples, silhouette >= 0.40 is
    # "clean" (well-separated clusters), below is "noisy". This is a
    # flag only — the chat-14 HPC run on Quentin's cohort will
    # recalibrate the 0.40 cutoff against the observed distribution.
    # C01d / C01f must NOT gate on this for the first wiring pass.
    decomp_quality = if (is.na(silhouette)) NA_character_
                      else if (silhouette >= 0.40) "clean"
                      else "noisy",
    bic_gap_k3_vs_k2 = round(bic_gap, 4),
    phase_concordance = round(phase_concordance, 3),
    n_phase_supported = as.integer(n_phase_support),
    # Flashlight integration
    sv_prior_mode = sv_prior_mode,
    n_seed_HOM_REF = as.integer(km$n_seed_ref %||% length(seeds$HOM_REF)),
    n_seed_HET     = as.integer(km$n_seed_het %||% length(seeds$HET)),
    n_seed_HOM_INV = as.integer(km$n_seed_inv %||% length(seeds$HOM_INV)),
    n_discordant = as.integer(length(discordant)),
    n_cheat2_reclassified = as.integer(n_cheat2_reclassified),
    discordant_samples = discordant,
    cheat2_constrained_samples = cheat2_constrained,
    # Centers + per-sample + per-window
    cluster_centers_pc1 = round(centers_sorted, 4),
    per_sample = per_sample,
    # Don't include the full per-window × per-sample matrix in JSON
    # — too big. Save it as an RDS file alongside.
    per_window_class_file = "per_window_class.rds"
  )

  # Save the big per-window matrix separately
  rds_path <- file.path(cand_outdirs[ci], "per_window_class.rds")
  saveRDS(list(
    per_window_class = pw$per_window_class,
    window_starts = loadings$win_starts,
    window_ends = loadings$win_ends,
    window_mids = loadings$win_mids,
    sample_names = sample_names,
    class_assignment_global = class_assignment,
    centers = centers_sorted
  ), rds_path)

  # ── Write Tier-2 block ──
  write_block_safe(
    reg = reg,
    candidate_id = cid,
    block_type = "internal_dynamics",
    data = block_data,
    source_script = "STEP_C01i_decompose.R"
  )
}

cat("\n[C01i_decompose] done — ", nrow(cand_dt), " candidates processed\n")
