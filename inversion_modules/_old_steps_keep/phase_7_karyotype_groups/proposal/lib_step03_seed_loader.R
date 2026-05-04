#!/usr/bin/env Rscript
# =============================================================================
# lib_step03_seed_loader.R — Tier 1 k-means seeding for C01i decompose
# =============================================================================
#
# WHAT "SEED" MEANS IN THIS FILE
# ------------------------------
# "Seed" here is a k-MEANS CENTROID SEED. C01i_decompose runs k-means (k=3)
# on per-sample PC1 values across the candidate inversion to cluster samples
# into HOM_REF / HET / HOM_INV groups. k-means needs initial centroids.
# Random initialization produces unstable / arbitrary class labels. This
# loader supplies principled initial centroids by pre-assigning confident
# samples to each of the three classes from external priors. k-means then
# starts from the means of those pre-assigned subsets.
#
# IT IS NOT: a random seed, a PRNG seed, the "seed" from the seed-extender
#            in GHSL, or a genomic seed region. It is the anchor sample set
#            that tells k-means "start HOM_REF cluster here, HET there,
#            HOM_INV there."
#
# The companion function that actually consumes the seed set lives in
# STEP_C01i_decompose.R::seeded_kmeans() — see that function for the exact
# centroid-init math.
#
# THE TWO SEED SOURCES
# --------------------
# This loader combines TWO INDEPENDENT PRIORS, each of which can pre-assign
# samples to HOM_REF / HET / HOM_INV classes based on different evidence:
#
#   1) sv_prior  (Cheat 1) — derived from SV-caller dosage signal.
#      Pre-computed upstream; passed in as the `sv_prior_seeds` arg. This
#      is the already-existing pathway.
#
#   2) STEP03    — short for the producer script
#          inversion_modules/phase_3_refine/STEP_D03_statistical_tests_and_seeds.py
#      Runs Fisher's exact test + Cochran-Armitage trend test on SV support
#      counts per arrangement class, emits per-candidate TSVs at
#          <step03_seeds_dir>/<inv_id>_seeds.tsv
#      with columns (sample_id, assigned_class, p_value, n_supporting_snps).
#      Samples surviving p ≤ 0.05 AND n_supporting_snps ≥ 5 are proposed as
#      seeds for that class. This is the INDEPENDENT statistical-evidence
#      pathway added by chat 9.
#
# The two sources ask different questions:
#   sv_prior: "Does this sample's dosage pattern LOOK LIKE a
#              REF/HET/INV carrier across the candidate?"
#   STEP03:   "Does this sample's per-SV-call support COUNT AS significant
#              for a REF/HET/INV assignment under Fisher/Armitage?"
# Both converging on the same class for a sample is strong evidence.
# Disagreement means neither is trustworthy enough to anchor k-means →
# drop that sample from the seed set (conservative).
#
# =============================================================================
# DESIGN HISTORY (chat-9 / chat-13)
# =============================================================================
# Implements chat-9 design §Tier 1 (Finding S closure):
#
#   "Add STEP03-seeded path. When STEP03's seed file at
#    04_deconvolution_seeds/{inv_id}_seeds.tsv is available, use its REF/HET/
#    INV sample assignments to seed k-means centroids the same way sv_prior
#    does. STEP03 seeds come from Fisher/Armitage significance — independent
#    evidence from SV calls."
#
#   "Combine the two seed sources when both available. Priority:
#    sv_prior > STEP03 if they conflict for a given sample. Union the
#    agreeing subsets; drop the disagreeing ones from the seed set rather
#    than forcing a pick."
#
# NOTE (chat-13 Finding AY): the implemented behaviour is stricter than
# the quoted chat-9 language. `combine_tier1_seeds()` DROPS conflicting
# samples from BOTH sources rather than PICKING sv_prior's assignment
# over STEP03's. Drop-on-conflict is safer than priority-pick because it
# doesn't manufacture a seed from disagreement: if the two sources
# disagree on a sample, neither source is confident enough to anchor
# the k-means, so that sample stays out of the seed set entirely. The
# "priority: sv_prior > STEP03" phrasing above is preserved here for
# traceability back to the chat-9 design doc but should be read as
# "sv_prior wins only in the DEGENERATE agreement case; genuine
# conflict is always resolved by dropping both".
#
#   "Remove the unsupervised k-means fallback. If neither seed source is
#    available for a candidate, emit decomp_status = 'no_seeding' and skip."
#
# =============================================================================
# USAGE in STEP_C01i_decompose.R (inline patch):
#
#   source("lib_step03_seed_loader.R")
#   seeds <- combine_tier1_seeds(
#     chr = chr, cid = cid, start_bp = s, end_bp = e,
#     sv_prior_seeds = sv_prior_seeds,   # existing Cheat 1 output
#     step03_seeds_dir = paths$step03_seeds_dir
#   )
#   if (seeds$status == "no_seeding") {
#     cat("[C01i_decompose]   SKIP: ", seeds$reason, "\n"); next
#   }
#   # seeds$seeds is a list(HOM_REF, HET, HOM_INV) of sample IDs — used as
#   # the anchor set for the three k-means centroids.
#   km <- seeded_kmeans(mean_pc1, seeds$seeds, k = opt$k_decomp)
#
# =============================================================================
# GLOSSARY (quick reference)
# =============================================================================
#   seed             = sample pre-assigned to a class, used as centroid anchor
#                      for the k=3 k-means that decomposes the cohort into
#                      HOM_REF / HET / HOM_INV
#   sv_prior         = seed source #1: Cheat 1 dosage-pattern prior
#   STEP03           = seed source #2: Fisher+Armitage statistical test run
#                      by STEP_D03_statistical_tests_and_seeds.py
#   agreement        = both sources assign the same class to a sample → keep
#   conflict         = sources assign different classes → DROP (both)
#   no_seeding       = fewer than min_seeds_per_class (default 2) samples
#                      survive in ANY of the three classes → skip candidate
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Load STEP03 seed file for a candidate (source #2 of the two seed sources)
# =============================================================================
# STEP03 = inversion_modules/phase_3_refine/STEP_D03_statistical_tests_and_seeds.py
# That script runs Fisher's exact test + Cochran-Armitage trend test per
# candidate and emits one TSV per inversion at:
#   <step03_seeds_dir>/<inv_id>_seeds.tsv
#
# Expected columns:
#   sample_id, assigned_class ∈ {REF, HET, INV} (or HOM_REF/HOM_INV/HET),
#   p_value, odds_ratio, n_supporting_snps
#
# Filters samples where p_value <= 0.05 AND n_supporting_snps >= min_snps
# (default 5). Anything below this is considered unreliable and dropped
# from the seed set — such samples will not anchor any k-means centroid.
#
# Returns a named list(HOM_REF, HET, HOM_INV) of sample_id character vectors.
# Returns NULL if the file doesn't exist (→ caller falls back to sv_prior
# alone, or if that's also absent, no_seeding and the candidate is skipped).
# =============================================================================
load_step03_seeds <- function(inv_id, step03_seeds_dir,
                                 p_value_thresh = 0.05,
                                 min_supporting_snps = 5L) {
  if (is.null(step03_seeds_dir) || !nzchar(step03_seeds_dir)) return(NULL)
  f <- file.path(step03_seeds_dir, paste0(inv_id, "_seeds.tsv"))
  if (!file.exists(f)) {
    # Also try candidate_id with dot-to-underscore normalization
    f2 <- file.path(step03_seeds_dir, paste0(gsub("\\.", "_", inv_id),
                                               "_seeds.tsv"))
    if (file.exists(f2)) f <- f2 else return(NULL)
  }
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (!"sample_id" %in% names(dt) || !"assigned_class" %in% names(dt)) {
    warning("[step03_seeds] ", f, " missing sample_id or assigned_class")
    return(NULL)
  }

  # Apply filters if columns exist
  if ("p_value" %in% names(dt) && !is.null(p_value_thresh)) {
    dt <- dt[!is.na(p_value) & p_value <= p_value_thresh]
  }
  if ("n_supporting_snps" %in% names(dt) && !is.null(min_supporting_snps)) {
    dt <- dt[!is.na(n_supporting_snps) &
             n_supporting_snps >= min_supporting_snps]
  }
  if (!nrow(dt)) return(NULL)

  # Normalize class labels to canonical HOM_REF/HET/HOM_INV
  dt[, class_norm := fifelse(assigned_class %in% c("REF", "HOM_REF"), "HOM_REF",
                      fifelse(assigned_class %in% c("INV", "HOM_INV"), "HOM_INV",
                       fifelse(assigned_class == "HET", "HET",
                               NA_character_)))]
  dt <- dt[!is.na(class_norm)]
  if (!nrow(dt)) return(NULL)

  list(
    HOM_REF = unique(dt[class_norm == "HOM_REF"]$sample_id),
    HET     = unique(dt[class_norm == "HET"]$sample_id),
    HOM_INV = unique(dt[class_norm == "HOM_INV"]$sample_id)
  )
}

# =============================================================================
# Combine sv_prior + STEP03 seeds with priority and disagreement rules
# =============================================================================
# "Tier 1 seeds" = the final sample→class assignments that will anchor the
# three k-means centroids in C01i_decompose::seeded_kmeans(). "Tier 1"
# distinguishes this from any downstream tier that further refines the
# k-means output (e.g. Tier 2 validation in C01f).
#
# chat-9 §Tier 1 rule: "Union the agreeing subsets; drop the disagreeing
# ones from the seed set rather than forcing a pick. Priority: sv_prior >
# STEP03 if they conflict."  See top-of-file NOTE re: chat-13 finding AY
# — genuine conflict drops both sources, priority-pick applies only in
# agreement cases.
#
# Concretely, for each sample:
#   - in BOTH sources with SAME class  → keep in seed set for that class
#   - in BOTH sources with DIFFERENT class → DROP from seed set entirely
#     (conservative: if the two sources disagree, neither can be trusted
#      to anchor k-means for this sample)
#   - in ONLY ONE source → keep with that source's class
# Final seed list per class = union(agreeing samples, single-source samples).
#
# Returns list:
#   status ∈ {"ok", "sv_prior_only", "step03_only", "no_seeding"}
#     ok            — both sources available, seeds populated from the
#                     combined set, all three classes have ≥ min_seeds_per_class
#     sv_prior_only — only sv_prior available, passed through
#     step03_only   — only STEP03 available, passed through
#     no_seeding    — neither available OR combined set too sparse after
#                     conflict drops → caller MUST skip the candidate
#                     (unsupervised k-means fallback was removed per chat-9)
#   reason     — character string if no_seeding, else NA
#   seeds      — list(HOM_REF, HET, HOM_INV) of sample_ids, NULL if no_seeding
#   n_agree    — count of samples where BOTH sources agreed on class
#   n_conflict — count of samples where sources DISAGREED (dropped)
#   source     — character label, useful for logging
# =============================================================================
combine_tier1_seeds <- function(chr, cid, start_bp, end_bp,
                                   sv_prior_seeds = NULL,
                                   step03_seeds_dir = NULL,
                                   min_seeds_per_class = 2L) {
  # Coerce sv_prior to the canonical list shape if needed
  fl <- if (!is.null(sv_prior_seeds) &&
            all(c("HOM_REF", "HET", "HOM_INV") %in% names(sv_prior_seeds))) {
    lapply(sv_prior_seeds[c("HOM_REF","HET","HOM_INV")], function(x) {
      if (is.null(x)) character(0) else as.character(x)
    })
  } else NULL

  st <- load_step03_seeds(cid, step03_seeds_dir)

  has_fl <- !is.null(fl) && sum(lengths(fl)) > 0
  has_st <- !is.null(st) && sum(lengths(st)) > 0

  if (!has_fl && !has_st) {
    return(list(
      status = "no_seeding",
      reason = "neither sv_prior nor STEP03 seeds available",
      seeds  = NULL,
      n_agree = 0L, n_conflict = 0L,
      source = "none"
    ))
  }

  if (has_fl && !has_st) {
    ok <- vapply(fl, function(x) length(x) >= min_seeds_per_class,
                 logical(1))
    if (!all(ok)) {
      return(list(
        status = "no_seeding",
        reason = paste0("sv_prior seeds too sparse (classes with ≥",
                         min_seeds_per_class, ": ",
                         sum(ok), "/3)"),
        seeds = NULL, n_agree = 0L, n_conflict = 0L,
        source = "sv_prior_sparse"
      ))
    }
    return(list(status = "sv_prior_only", reason = NA_character_,
                seeds = fl, n_agree = 0L, n_conflict = 0L,
                source = "sv_prior"))
  }

  if (!has_fl && has_st) {
    ok <- vapply(st, function(x) length(x) >= min_seeds_per_class,
                 logical(1))
    if (!all(ok)) {
      return(list(
        status = "no_seeding",
        reason = paste0("STEP03 seeds too sparse (classes with ≥",
                         min_seeds_per_class, ": ",
                         sum(ok), "/3)"),
        seeds = NULL, n_agree = 0L, n_conflict = 0L,
        source = "step03_sparse"
      ))
    }
    return(list(status = "step03_only", reason = NA_character_,
                seeds = st, n_agree = 0L, n_conflict = 0L,
                source = "step03"))
  }

  # Both available — union agreements, drop conflicts
  # Build a map: sample_id → sv_prior class
  fl_map <- c(
    setNames(rep("HOM_REF", length(fl$HOM_REF)), fl$HOM_REF),
    setNames(rep("HET",     length(fl$HET)),     fl$HET),
    setNames(rep("HOM_INV", length(fl$HOM_INV)), fl$HOM_INV)
  )
  st_map <- c(
    setNames(rep("HOM_REF", length(st$HOM_REF)), st$HOM_REF),
    setNames(rep("HET",     length(st$HET)),     st$HET),
    setNames(rep("HOM_INV", length(st$HOM_INV)), st$HOM_INV)
  )

  all_samples <- union(names(fl_map), names(st_map))
  agreed <- list(HOM_REF = character(0), HET = character(0), HOM_INV = character(0))
  n_agree    <- 0L
  n_conflict <- 0L

  for (sid in all_samples) {
    cls_fl <- fl_map[sid]
    cls_st <- st_map[sid]
    if (is.na(cls_fl) && !is.na(cls_st)) {
      agreed[[cls_st]] <- c(agreed[[cls_st]], sid)
    } else if (!is.na(cls_fl) && is.na(cls_st)) {
      agreed[[cls_fl]] <- c(agreed[[cls_fl]], sid)
    } else if (!is.na(cls_fl) && !is.na(cls_st)) {
      if (cls_fl == cls_st) {
        agreed[[cls_fl]] <- c(agreed[[cls_fl]], sid)
        n_agree <- n_agree + 1L
      } else {
        n_conflict <- n_conflict + 1L
      }
    }
  }

  ok <- vapply(agreed, function(x) length(x) >= min_seeds_per_class,
               logical(1))
  if (!all(ok)) {
    return(list(
      status = "no_seeding",
      reason = paste0("combined seeds too sparse after conflict drops (",
                       n_conflict, " conflicts dropped; classes with ≥",
                       min_seeds_per_class, ": ", sum(ok), "/3)"),
      seeds = NULL, n_agree = n_agree, n_conflict = n_conflict,
      source = "combined_sparse"
    ))
  }

  list(
    status     = "ok",
    reason     = NA_character_,
    seeds      = agreed,
    n_agree    = n_agree,
    n_conflict = n_conflict,
    source     = "sv_prior+step03"
  )
}
