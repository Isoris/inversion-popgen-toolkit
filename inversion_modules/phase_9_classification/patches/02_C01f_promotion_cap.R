# =============================================================================
# PATCH: C01f — respect q6_validation_promotion_cap from Phase 4b
# =============================================================================
# Applies to: STEP_C01f_hypothesis_tests (already patched in v10 with
#              compute_group_validation).
# Purpose:    C01i_d_seal sets q6_validation_promotion_cap = "UNCERTAIN"
#             for composite candidates (composite_flag = likely_composite).
#             C01f must READ this cap before attempting any promotion.
#             Without this patch, a composite interval could still get
#             promoted to SUPPORTED / VALIDATED on strong T8/T9/Layer D,
#             producing confidently-wrong groups.
#
# This is ~30 lines of modification to the compute_group_validation()
# function from v10. Apply by replacing the function body.
# =============================================================================


# ─── REPLACE compute_group_validation() with this version ───────────────────

compute_group_validation <- function(current_level,
                                       t8_concordance = NA,
                                       t9_jackknife_status = NA_character_,
                                       t9_max_delta = NA,
                                       t10_theta_concordance = NA,
                                       layer_d_fisher_p = NA,
                                       layer_d_fisher_or = NA,
                                       promotion_cap = NA_character_) {
  # Pure function. Same logic as v10 + a new promotion_cap argument.
  # If promotion_cap is set (e.g., "UNCERTAIN"), the function cannot
  # return a level ranked higher than the cap — regardless of how strong
  # the evidence is.
  #
  # promotion_cap is typically set by STEP_C01i_d_seal based on the
  # internal_ancestry_composition block: likely_composite → cap=UNCERTAIN,
  # because groups from a composite interval are unreliable and shouldn't
  # be promoted on the strength of tests that don't know about the
  # compositeness.

  if (is.na(current_level) || !nzchar(current_level)) {
    current_level <- "UNCERTAIN"
  }

  # Demotion check first (same as v10) — SUSPECT bypasses the cap because
  # it's a demotion below any cap anyway.
  is_fragile <- identical(t9_jackknife_status, "fragile") ||
                (is.finite(t9_max_delta) && t9_max_delta > 0.3)
  if (is_fragile) return("SUSPECT")

  # Candidate promotions (same rules as v10)
  candidate_levels <- character(0)
  if (is.finite(layer_d_fisher_p) && is.finite(layer_d_fisher_or) &&
      layer_d_fisher_p < 0.05 && layer_d_fisher_or > 5) {
    candidate_levels <- c(candidate_levels, "VALIDATED")
  }
  if (is.finite(t8_concordance) && t8_concordance >= 0.70 &&
      identical(t9_jackknife_status, "robust")) {
    candidate_levels <- c(candidate_levels, "SUPPORTED")
  }

  rank_map <- c(NONE = 0L, UNCERTAIN = 1L, SUPPORTED = 2L, VALIDATED = 3L,
                 SUSPECT = -1L)
  all_levels <- c(current_level, candidate_levels)
  ranks <- rank_map[all_levels]
  best <- names(all_levels)[which.max(ranks)]

  # NEW: apply promotion_cap if set. best is replaced by min(best, cap)
  # in the ranking order.
  if (!is.na(promotion_cap) && nzchar(promotion_cap) &&
      promotion_cap %in% names(rank_map)) {
    cap_rank <- rank_map[promotion_cap]
    best_rank <- rank_map[best]
    if (!is.na(best_rank) && !is.na(cap_rank) && best_rank > cap_rank) {
      # Log before capping — useful for audit
      message("[C01f] promotion cap ", promotion_cap, " overrides ", best)
      best <- promotion_cap
    }
  }

  best
}

# ─── CALL SITE UPDATE ───────────────────────────────────────────────────────
#
# Where C01f currently calls compute_group_validation(), pass the cap
# from the registry:
#
#   cap_row <- reg$evidence$get_evidence(cid, "q6_validation_promotion_cap")
#   cap <- if (!is.null(cap_row) && length(cap_row) > 0 && nrow(cap_row) > 0) {
#     cap_row$value[1]
#   } else NA_character_
#
#   new_level <- compute_group_validation(
#     current_level         = current_level,
#     t8_concordance        = t8$concordance,
#     t9_jackknife_status   = t9$status,
#     t9_max_delta          = t9$max_delta,
#     t10_theta_concordance = t10$concordance,
#     layer_d_fisher_p      = layer_d_fisher_p,
#     layer_d_fisher_or     = layer_d_fisher_or,
#     promotion_cap         = cap
#   )
#
# That's the only change at the call site. The promotion logic is now
# aware of the composite_flag without C01f needing to read
# internal_ancestry_composition directly — the cap is a single flat key.
