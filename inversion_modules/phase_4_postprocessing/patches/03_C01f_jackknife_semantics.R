# =============================================================================
# PATCH: compute_group_validation — richer jackknife semantics
# =============================================================================
# Applies to: the compute_group_validation() function in
#             v10 patches/C01f/C01f_comp_from_registry_patch.R
#             (which v10.1 already extended with promotion_cap).
#
# REASON:
# The old rule `if t9_status == "fragile" → SUSPECT` conflated two things:
#   (a) the PCA groups are driven by family structure, not arrangement
#       → the inversion is probably not real at all OR the groups are wrong
#   (b) the inversion is carried by one family → it IS real, just restricted
#
# These should not be treated the same way. The cheat6 jackknife verdict
# already distinguishes them. This patch propagates that distinction into
# the validation level + quality flags.
#
# jackknife_verdict values (from cheat6_ancestry_jackknife_v934.R):
#   robust_multi_family         inversion segregates in multiple families
#                               → safe to promote, no concerns
#   multi_family_contributing   ≥3 groups contribute but sensitivity varies
#                               → safe to promote but tag "partial_robustness"
#   few_family                  2-3 groups drive signal
#                               → partial; tag "restricted_family_spread"
#                                  do not promote above SUPPORTED
#   single_family_fragile       one family dominates
#                               → tag "family_specific_polymorphism"
#                                  do not promote to VALIDATED but STAY at
#                                  UNCERTAIN or SUPPORTED (real but restricted)
#   (internal fallback) fragile  genuinely suspect, same as before → SUSPECT
#
# Critical insight: "single_family_fragile" is NOT equivalent to "fragile"
# in the confidence sense. It's a biological category — family-restricted
# polymorphism — that deserves its own handling.
# =============================================================================

# ─── REPLACE compute_group_validation() with this version (v2.1) ────────────

compute_group_validation <- function(current_level,
                                       t8_concordance = NA,
                                       t9_jackknife_verdict = NA_character_,
                                       t9_max_delta = NA,
                                       t10_theta_concordance = NA,
                                       layer_d_fisher_p = NA,
                                       layer_d_fisher_or = NA,
                                       promotion_cap = NA_character_) {
  # Returns: list(level = <LEVEL>, quality_flags = <chr vec>, family_linkage = <str>)
  # This version returns a richer object instead of just a string. Callers
  # that only want the level can use result$level.

  if (is.na(current_level) || !nzchar(current_level)) current_level <- "UNCERTAIN"

  quality_flags <- character(0)
  family_linkage <- "unknown"
  internal_cap <- NA_character_  # local cap derived from jackknife

  # ── Classify the jackknife verdict ──
  jk <- tolower(t9_jackknife_verdict %||% "")
  if (jk == "robust_multi_family") {
    family_linkage <- "multi_family"
    # No restriction on promotion
  } else if (jk == "multi_family_contributing") {
    family_linkage <- "multi_family"
    quality_flags <- c(quality_flags, "partial_robustness")
  } else if (jk == "few_family") {
    family_linkage <- "few_family"
    quality_flags <- c(quality_flags, "restricted_family_spread")
    internal_cap <- "SUPPORTED"  # don't validate, but it's real
  } else if (jk == "single_family_fragile") {
    family_linkage <- "single_family"
    quality_flags <- c(quality_flags, "family_specific_polymorphism")
    # Keep at UNCERTAIN or SUPPORTED — this IS a real inversion, just
    # restricted to one family. Do not promote to VALIDATED.
    internal_cap <- "SUPPORTED"
  } else if (jk == "fragile" || jk == "pca_family_driven") {
    # Genuinely suspect: PCA groups track ancestry, not arrangement.
    # This is the old "SUSPECT" path.
    family_linkage <- "pca_family_confounded"
    quality_flags <- c(quality_flags, "pca_groups_family_confounded")
    return(list(level = "SUSPECT",
                quality_flags = quality_flags,
                family_linkage = family_linkage))
  }

  # Fallback: raw t9_max_delta-only signal (when verdict missing)
  if (!nzchar(jk) && is.finite(t9_max_delta) && t9_max_delta > 0.3) {
    # Without a verdict string, we can't tell "family-specific but real"
    # from "pca_family_confounded". Be conservative and treat it as family
    # linkage uncertain + cap.
    quality_flags <- c(quality_flags, "high_jackknife_delta_no_verdict")
    internal_cap <- "SUPPORTED"
    family_linkage <- "uncertain"
  }

  # ── Candidate promotions (same rules as v10) ──
  candidate_levels <- character(0)
  if (is.finite(layer_d_fisher_p) && is.finite(layer_d_fisher_or) &&
      layer_d_fisher_p < 0.05 && layer_d_fisher_or > 5) {
    candidate_levels <- c(candidate_levels, "VALIDATED")
  }
  if (is.finite(t8_concordance) && t8_concordance >= 0.70 &&
      identical(jk, "robust_multi_family")) {
    candidate_levels <- c(candidate_levels, "SUPPORTED")
  }
  # NEW: also allow SUPPORTED promotion on single_family inversion if
  # T8 Clair3 concordance is high — the family-restricted pattern IS
  # the inversion, Clair3 agreeing confirms the groups.
  if (is.finite(t8_concordance) && t8_concordance >= 0.70 &&
      jk %in% c("single_family_fragile", "few_family")) {
    candidate_levels <- c(candidate_levels, "SUPPORTED")
  }

  rank_map <- c(NONE = 0L, UNCERTAIN = 1L, SUPPORTED = 2L, VALIDATED = 3L,
                 SUSPECT = -1L)
  all_levels <- c(current_level, candidate_levels)
  ranks <- rank_map[all_levels]
  best <- names(all_levels)[which.max(ranks)]

  # ── Apply caps (external promotion_cap AND internal jackknife cap) ──
  for (cap in c(promotion_cap, internal_cap)) {
    if (!is.na(cap) && nzchar(cap) && cap %in% names(rank_map)) {
      if (rank_map[best] > rank_map[cap]) {
        best <- cap
      }
    }
  }

  list(
    level = best,
    quality_flags = if (length(quality_flags) == 0) "normal" else quality_flags,
    family_linkage = family_linkage
  )
}

# ─── CALL SITE ───────────────────────────────────────────────────────────────
#
# The new return type is a list, not a string. Update call sites:
#
#   OLD:
#     new_level <- compute_group_validation(...)
#     vd[, group_validation_after := new_level]
#
#   NEW:
#     v <- compute_group_validation(
#       current_level         = current_level,
#       t8_concordance        = t8$concordance,
#       t9_jackknife_verdict  = t9$verdict,             # <-- the string, not status
#       t9_max_delta          = t9$max_delta,
#       t10_theta_concordance = t10$concordance,
#       layer_d_fisher_p      = layer_d_fisher_p,
#       layer_d_fisher_or     = layer_d_fisher_or,
#       promotion_cap         = cap
#     )
#     vd[, group_validation_after := v$level]
#     vd[, quality_flags          := paste(v$quality_flags, collapse = ",")]
#     vd[, family_linkage         := v$family_linkage]
#
# The family_linkage value feeds the frequency.json block's family_linkage
# field (replacing the old binary single_family/multi_family axis with the
# four-class one above). classify_inversions picks this up directly.
#
# IMPORTANT: the hypothesis_verdict schema must expand its family_linkage
# enum to include:
#   {multi_family, few_family, single_family, pca_family_confounded, unknown}
#
# The v10.1 schema update is in:
#   schemas/hypothesis_verdict.schema.json (bumped to v2)
