# =============================================================================
# PATCH: STEP_C01f_hypothesis_tests — comp_from_registry + validation promotion
# =============================================================================
# Applies to: STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R
# Purpose:    (1) Replace internal PC1 k-means-based carrier assignment with
#                 registry-lookup of groups written by C01i.
#             (2) After all tests complete, promote q6_group_validation to
#                 SUPPORTED / VALIDATED or demote to SUSPECT.
#             (3) Keep the old k-means logic as a fallback for backward
#                 compatibility (if C01i has not run, or if groups are empty).
#
# Apply as:   Three separate str_replace operations — marked (A), (B), (C).
#             (A) adds two new functions near the top of the script
#             (B) patches the four call sites that build `comp` from k-means
#             (C) adds the validation-promotion block at end of per-candidate loop
# =============================================================================


# =============================================================================
# (A) ADD NEW HELPER FUNCTIONS — insert after the args-parsing block,
#     before the first test_* function definition. Easy marker: look for the
#     line `# T3: Kin-pruned similarity matrix comparison` and insert
#     the block below immediately BEFORE it.
# =============================================================================

# ─── comp_from_registry ───────────────────────────────────────────────────────
# Build the `comp` data.table (columns: sample, band, source) from C01i's
# registered groups for this candidate. Returns NULL if groups are absent
# or insufficient, so the caller can fall back to k-means.
#
# Band mapping convention (matches old k-means ordering by mean PC1):
#   band1 = HOM_REF    (lowest PC1)
#   band2 = HET        (middle PC1)
#   band3 = HOM_INV    (highest PC1)
#   RECOMBINANT samples are EXCLUDED — they violate the
#   "one band per sample" assumption of T1/T2/T7/T8/T10. Their class
#   switches inside the interval, so averaging them into any band
#   pollutes the within-band statistic. They are handled separately in
#   Q2 internal dynamics, not here.
comp_from_registry <- function(cid, sample_names,
                                 exclude_recombinants = TRUE,
                                 min_per_band = 3L) {
  if (!exists("reg") || is.null(reg) || is.null(reg$has_group)) {
    return(NULL)
  }
  # Accept either HOM_REF (canonical, v10) or HOM_STD (alias, v9) group names.
  ref_gid <- if (reg$has_group(paste0("inv_", cid, "_HOM_REF"))) {
    paste0("inv_", cid, "_HOM_REF")
  } else if (reg$has_group(paste0("inv_", cid, "_HOM_STD"))) {
    paste0("inv_", cid, "_HOM_STD")
  } else {
    NULL
  }
  het_gid  <- paste0("inv_", cid, "_HET")
  inv_gid  <- paste0("inv_", cid, "_HOM_INV")
  rec_gid  <- paste0("inv_", cid, "_RECOMBINANT")

  if (is.null(ref_gid) || !reg$has_group(het_gid) || !reg$has_group(inv_gid)) {
    return(NULL)  # need all three core classes to proceed
  }

  ref_samps <- reg$get_group(ref_gid)
  het_samps <- reg$get_group(het_gid)
  inv_samps <- reg$get_group(inv_gid)
  rec_samps <- if (reg$has_group(rec_gid)) reg$get_group(rec_gid) else character(0)

  # Sanity: require min_per_band members in each band or fall back
  if (length(ref_samps) < min_per_band ||
      length(het_samps) < min_per_band ||
      length(inv_samps) < min_per_band) {
    return(NULL)
  }

  # Intersect with sample_names so we return only samples present in this
  # run (e.g., after kin-pruning).
  ref_samps <- intersect(ref_samps, sample_names)
  het_samps <- intersect(het_samps, sample_names)
  inv_samps <- intersect(inv_samps, sample_names)

  # Build comp. Column names match what test_relatedness(), test_ancestry_diversity(),
  # test_carrier_substructure(), test_theta_concordance(), and
  # test_genotype_concordance() currently expect: columns `sample` and `band`.
  comp <- data.table::data.table(
    sample = c(ref_samps, het_samps, inv_samps),
    band   = c(rep("band1", length(ref_samps)),
               rep("band2", length(het_samps)),
               rep("band3", length(inv_samps))),
    source = "C01i_registry"
  )

  # If any extra columns were expected by specific tests (e.g. `eff_K` for T2),
  # they are attached by the caller after this. The base shape (sample, band)
  # is what all tests agree on.
  attr(comp, "recombinants_excluded") <- rec_samps
  comp
}


# ─── build_comp_with_fallback ────────────────────────────────────────────────
# Convenience wrapper: try registry first, fall back to k-means on PC1.
# `km_fn` is a closure that runs the old inline k-means logic for this
# candidate and returns a `comp` data.table in the same format. Each call
# site passes its own closure (the 4 call sites had slightly different
# inline code — see patch B).
build_comp_with_fallback <- function(cid, sample_names, km_fn) {
  comp <- comp_from_registry(cid, sample_names)
  if (!is.null(comp) && nrow(comp) >= 9) {  # at least 3 per band
    message("[C01f] ", cid, ": using C01i groups from registry (",
            nrow(comp), " samples, ",
            length(attr(comp, "recombinants_excluded")), " recombinants excluded)")
    return(comp)
  }
  message("[C01f] ", cid, ": registry groups absent/insufficient — falling back to PC1 k-means")
  km_comp <- km_fn()
  if (!is.null(km_comp)) attr(km_comp, "source") <- "kmeans_fallback"
  km_comp
}


# =============================================================================
# (B) PATCH THE FOUR CALL SITES that currently build `comp` via inline k-means.
# =============================================================================
#
# In v9.3.4 C01f, `comp` is built in four places (grep for `kmeans(...)`):
#
#   1. In the per-candidate main loop, around line ~1800-ish, where `comp`
#      is first constructed for passage to T1/T2/T3/T7/T8/T10.
#   2. In test_inner_outer_overlap() T4 — rebuilds a local comp for
#      the inner/outer comparison. LEAVE AS-IS: T4 needs the REGIONAL
#      composition, not the global genotype groups. Comment added.
#   3. In test_anchor_stability() T5 — builds per-subregion anchor bands.
#      LEAVE AS-IS: T5 tests whether bands are stable across sub-regions,
#      which requires per-subregion re-clustering. Comment added.
#   4. In test_regime_change() T6 — per-boundary regime k-means.
#      LEAVE AS-IS: same reason as T5.
#
# So only call site (1) needs the patch. T4/T5/T6 get explanatory comments.
#
# ─── CALL SITE 1 (main loop, the ONE patch) ───────────────────────────────
# Find the block in the candidate-loop that looks like (paraphrased):
#
#     # Build band assignment from PC1
#     pc_avg <- colMeans(dt[rows, ..available], na.rm = TRUE)
#     km <- kmeans(pc_avg, centers = 3, nstart = 10)
#     co <- order(km$centers[, 1])
#     bands <- character(length(pc_avg))
#     for (bi in seq_along(co)) bands[km$cluster == co[bi]] <- paste0("band", bi)
#     snames <- sub("^PC_1_", "", names(pc_avg))
#     comp <- data.table(sample = snames, band = bands)
#     ... (possibly add eff_K column here)
#
# Replace with:
#
#     # Build band assignment: prefer C01i's registered groups; fall back to
#     # k-means on PC1 if groups are absent (e.g., C01i did not run, or the
#     # candidate is Tier 4 so had no decomposition).
#     km_fallback <- function() {
#       pc_avg <- colMeans(dt[rows, ..available], na.rm = TRUE)
#       km <- kmeans(pc_avg, centers = 3, nstart = 10)
#       co <- order(km$centers[, 1])
#       bands <- character(length(pc_avg))
#       for (bi in seq_along(co)) bands[km$cluster == co[bi]] <- paste0("band", bi)
#       snames <- sub("^PC_1_", "", names(pc_avg))
#       data.table::data.table(sample = snames, band = bands, source = "kmeans_fallback")
#     }
#     comp <- build_comp_with_fallback(cid, sample_names, km_fallback)
#     if (is.null(comp)) {
#       message("[C01f] ", cid, ": no groups and k-means failed — skipping group-based tests")
#       next
#     }
#
# ─── eff_K column ─────────────────────────────────────────────────────────
# T2 (test_ancestry_diversity) wants an `eff_K` column on `comp`. After the
# comp is built (registry OR fallback), attach eff_K the same way as the
# old code does — with the Q-matrix lookup. That block doesn't change.
# Insert eff_K annotation AFTER build_comp_with_fallback returns.
#
# ─── CALL SITES 2, 3, 4 (T4, T5, T6) ──────────────────────────────────────
# Add a comment above each inline kmeans() call:
#
#     # DO NOT replace with registry lookup. T4/T5/T6 require PER-SUBREGION
#     # re-clustering by design — the test IS "do bands change across
#     # sub-regions?". Registry groups are GLOBAL across the candidate
#     # interval; they cannot answer that question.


# =============================================================================
# (C) VALIDATION PROMOTION — computed INSIDE the block, not via add_evidence
# =============================================================================
#
# DESIGN NOTE (learned from v10 sanity test):
# The validation level is a derived property of the hypothesis_verdict block.
# Writing it via reg$add_evidence() separately creates a race with the
# schema's keys_extracted directive (the block's group_validation_after
# field gets flattened to q6_group_validation by the library, and whichever
# write lands last wins). The correct pattern is:
#
#   1. Compute the new level as a local value
#   2. Put it into the block's data dict as `group_validation_after`
#   3. Let reg$evidence$write_block() extract it into q6_group_validation
#      via the schema's keys_extracted rule.
#
# The block is the single source of truth. The flat key is a projection.
#
# This function therefore RETURNS the new level and does not touch the
# registry. The caller puts the returned value into the block data.
#
# Insert AFTER the args-parsing block, near comp_from_registry().

compute_group_validation <- function(current_level,
                                       t8_concordance = NA,
                                       t9_jackknife_status = NA_character_,
                                       t9_max_delta = NA,
                                       t10_theta_concordance = NA,
                                       layer_d_fisher_p = NA,
                                       layer_d_fisher_or = NA) {
  # Pure function: compute the new level from test outcomes.
  # Order of logic:
  #   1. SUSPECT demotion wins (fragile jackknife beats any promotion)
  #   2. Otherwise, take the strongest promotion that triggers
  #   3. Otherwise, keep the current level
  #
  # Monotone upgrade ranking (excluding SUSPECT which is orthogonal):
  #   NONE(0) < UNCERTAIN(1) < SUPPORTED(2) < VALIDATED(3)

  if (is.na(current_level) || !nzchar(current_level)) {
    current_level <- "UNCERTAIN"
  }

  # Demotion check first
  is_fragile <- identical(t9_jackknife_status, "fragile") ||
                (is.finite(t9_max_delta) && t9_max_delta > 0.3)
  if (is_fragile) return("SUSPECT")

  # Promotion: compute candidate levels from independent evidence streams
  candidate_levels <- character(0)

  # Layer D Fisher OR test — strongest evidence
  if (is.finite(layer_d_fisher_p) && is.finite(layer_d_fisher_or) &&
      layer_d_fisher_p < 0.05 && layer_d_fisher_or > 5) {
    candidate_levels <- c(candidate_levels, "VALIDATED")
  }

  # Clair3 concordance + robust jackknife
  if (is.finite(t8_concordance) && t8_concordance >= 0.70 &&
      identical(t9_jackknife_status, "robust")) {
    candidate_levels <- c(candidate_levels, "SUPPORTED")
  }

  # Take the max of current + candidate levels
  rank_map <- c(NONE = 0L, UNCERTAIN = 1L, SUPPORTED = 2L, VALIDATED = 3L,
                 SUSPECT = -1L)
  all_levels <- c(current_level, candidate_levels)
  ranks <- rank_map[all_levels]
  names(all_levels)[which.max(ranks)]
}

# ─── CALL IT at end of per-candidate loop, BEFORE writing the block ──────
#
# Where v9.3.4 currently does something like:
#
#   vd <- data.table(candidate_id = cid, verdict = verdict,
#                    verdict_confidence = conf,
#                    t8_concordance = t8$concordance, ...)
#   fwrite(rbind(verd_dt, vd), file.path(outdir, "hypothesis_verdicts.tsv"))
#   register_C01f_keys(vd, cid, outdir)
#
# Change to:
#
#   # Compute the new validation level from test outcomes
#   current_level <- tryCatch({
#     ev <- reg$get_evidence(cid, "q6_group_validation")
#     if (nrow(ev) > 0) ev$value[1] else "UNCERTAIN"
#   }, error = function(e) "UNCERTAIN")
#
#   # Chat 7 FIX 31: read Layer D OR + p from registry (written by phase_3 STEP03).
#   # Never hard-code NA — that made the VALIDATED promotion gate unreachable.
#   layer_d_fisher_p <- tryCatch({
#     if (exists("reg") && !is.null(reg$evidence) && !is.null(reg$evidence$get_evidence)) {
#       ev <- reg$evidence$get_evidence(cid, "q7_layer_d_fisher_p")
#       if (!is.null(ev) && is.data.frame(ev) && nrow(ev) > 0) {
#         suppressWarnings(as.numeric(ev$value[1]))
#       } else NA_real_
#     } else NA_real_
#   }, error = function(e) NA_real_)
#   layer_d_fisher_or <- tryCatch({
#     if (exists("reg") && !is.null(reg$evidence) && !is.null(reg$evidence$get_evidence)) {
#       ev <- reg$evidence$get_evidence(cid, "q7_layer_d_fisher_or")
#       if (!is.null(ev) && is.data.frame(ev) && nrow(ev) > 0) {
#         suppressWarnings(as.numeric(ev$value[1]))
#       } else NA_real_
#     } else NA_real_
#   }, error = function(e) NA_real_)
#
#   new_level <- compute_group_validation(
#     current_level         = current_level,
#     t8_concordance        = t8$concordance,
#     t9_jackknife_status   = t9$status,
#     t9_max_delta          = t9$max_delta,
#     t10_theta_concordance = t10$concordance,
#     layer_d_fisher_p      = layer_d_fisher_p,   # NA-safe registry read
#     layer_d_fisher_or     = layer_d_fisher_or   # NA-safe registry read
#   )
#
#   message("[C01f] ", cid, ": group validation ", current_level, " → ", new_level)
#
#   # Include in verdict table
#   vd <- data.table(candidate_id = cid,
#                    verdict                  = verdict,
#                    verdict_confidence       = conf,
#                    t8_concordance           = t8$concordance,
#                    step9_jackknife          = t9$status,
#                    t9_max_delta             = t9$max_delta,
#                    t10_theta_concordance    = t10$concordance,
#                    t11_has_extended         = t11$has_extended,
#                    t11_extent_kb            = t11$extent_kb,
#                    t11_fst_decay_rate       = t11$fst_decay_rate,
#                    group_validation_before  = current_level,
#                    group_validation_after   = new_level,
#                    comp_source              = attr(comp, "source") %||% "unknown")
#
#   # v10 path: write Tier-2 block — schema extracts q6_group_validation from
#   # vd$group_validation_after automatically.
#   if (exists("reg", inherits = TRUE) && !is.null(reg$evidence) &&
#       !is.null(reg$evidence$write_block)) {
#     reg$evidence$write_block(
#       candidate_id = cid,
#       block_type   = "hypothesis_verdict",
#       data         = as.list(vd[1])
#     )
#   } else {
#     # v9 fallback: flat add_evidence + register_C01f_keys
#     register_C01f_keys(vd, cid, outdir)
#     tryCatch(reg$add_evidence(cid, "q6_group_validation",
#                                 value = new_level, script = "C01f"),
#              error = function(e) NULL)
#   }
#
# Note: `compute_group_validation()` is a pure function with no side effects.
# This makes it testable in isolation — the test in tests/test_compute_group_validation.R
# exercises all transitions without needing a real registry.


# =============================================================================
# NOTES ON RECOMBINANTS (cross-cutting concern)
# =============================================================================
#
# The new comp_from_registry() EXCLUDES recombinant samples from `comp`. This
# is correct for T1/T2/T3/T7/T8/T10 (these tests assume samples have a
# single class across the interval).
#
# BUT test_genotype_concordance() (T8) uses Clair3 GTs to check PCA bands.
# Because recombinants are excluded from `comp`, T8 will score only the
# ~219 clean samples. That's what we want — Clair3 concordance is a test
# of "can discrete genotypes predict the three clean classes?", not
# "can they also predict the messy recombinants?" The recombinants are
# handled by Q2 internal_dynamics (per-window class tracking).
#
# If you want a sanity check on recombinants, add an optional per-candidate
# attr that records which samples were excluded:
#
#   attr(comp, "recombinants_excluded") <- rec_samps
#
# This is already set by comp_from_registry(), so downstream code that
# cares can access it via attr(comp, "recombinants_excluded").
