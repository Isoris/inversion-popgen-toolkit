# =============================================================================
# GROUP VALIDATION GATE — Pre-check for characterization convergence rules
# =============================================================================
#
# Before evaluating whether a biological question is ANSWERED, check whether
# the GROUPS used to compute the keys are trustworthy enough for that question.
#
# Group validation levels:
#   VALIDATED:   OR test passed (p < 0.05) + at least one of (GHSL, jackknife)
#   SUPPORTED:   No OR test, but GHSL concordance > 0.7 + jackknife robust
#   UNCERTAIN:   Only PCA, no independent confirmation
#   SUSPECT:     family_likeness > 0.5 OR jackknife fragile
#   NONE:        No groups assigned yet (decomposition hasn't run)
#
# Each question has a MINIMUM group validation level:
#   Q1 (what is it):        NONE      — uses sim_mat structure, not groups
#   Q2 (what's inside):     UNCERTAIN — recombinant detection uses groups but
#                                        is self-validating (switching pattern)
#   Q3 (boundaries):        NONE      — boundary evidence is group-independent
#   Q4 (mechanism):         NONE      — junction sequence doesn't use groups
#   Q5 (age):               SUPPORTED — Fst needs reliable group separation;
#                                        GDS doesn't need groups (pairwise)
#   Q6 (frequency):         SUPPORTED — freq_inv is directly from group counts
#   Q7 (is it real):        NONE      — Layer D uses groups but IS the validation
#
# Special cases:
#   Q5 can partially answer with GDS alone (no group requirement)
#   but Fst/dXY/Tajima D per-arrangement require SUPPORTED groups.
#   So Q5 can be PARTIAL without groups but needs SUPPORTED for ANSWERED.
#
# =============================================================================

group_validation_minimum <- list(
  Q1 = "NONE",
  Q2 = "UNCERTAIN",
  Q3 = "NONE",
  Q4 = "NONE",
  Q5 = "SUPPORTED",   # for full answer; GDS alone can reach PARTIAL
  Q6 = "SUPPORTED",
  Q7 = "NONE"         # Q7 Layer D IS the group validation test
)

# Validation level ordering (for comparison)
validation_order <- c("NONE" = 0, "SUSPECT" = 1, "UNCERTAIN" = 2,
                       "SUPPORTED" = 3, "VALIDATED" = 4)

assess_group_validation <- function(keys) {
  # Determine group validation level from available keys

  # BUGFIX 2026-04-17 (chat 7, FIX 32): add fisher_or > 5 to the OR-pass
  # gate. Schema description requires both "OR>5, p<0.05"; C01f's
  # compute_group_validation() tests both (L525–526). Previously this
  # gate tested p only, which would award VALIDATED on a significant-
  # but-low-OR Fisher test that C01f would withhold. Default OR = 0
  # means a missing value blocks the gate (conservative).
  or_p <- suppressWarnings(as.numeric(keys[["q7_layer_d_fisher_p"]] %||% NA))
  or_val <- suppressWarnings(as.numeric(keys[["q7_layer_d_fisher_or"]] %||% NA))
  or_tested <- !is.na(or_p) && !is.na(or_val)
  or_passed <- or_tested && or_p < 0.05 && or_val > 5

  ghsl_conc <- suppressWarnings(as.numeric(keys[["q7_layer_c_ghsl_pct_pass"]] %||% NA))
  ghsl_high <- !is.na(ghsl_conc) && ghsl_conc > 0.7

  jk_status <- keys[["q7_t9_jackknife_status"]] %||% keys[["q7_t9_jackknife_verdict"]] %||% "unknown"
  jk_robust <- jk_status %in% c("robust_multi_family", "robust")

  # FIX 43 (chat 9): q7_t9_jackknife_* keys are aspirational (see
  # 4e/compute_candidate_status.R build_key_spec). The semantic info is
  # reliably available via q6_family_linkage (written by 4b seal initially,
  # overwritten by 4c C01f after jackknife). Prefer that path.
  family_linkage <- keys[["q6_family_linkage"]] %||% "unknown"
  if (family_linkage %in% c("multi_family")) {
    jk_robust <- TRUE   # family_linkage=multi_family IS the robust case
  }

  family_lik <- suppressWarnings(as.numeric(keys[["q1_family_likeness_mean"]] %||% 0))
  # FIX 43 (chat 9): canonical fragile signal is family_linkage=pca_family_confounded
  # (written by C01f after jackknife, per patch 03 semantics). v9 fallback
  # to jk_status==single_family_fragile kept but secondary.
  jk_fragile <- family_linkage == "pca_family_confounded" ||
                jk_status == "single_family_fragile"

  # FIX 43 (chat 9): has_groups check accept both v10 HOM_REF and v9 HOM_STD
  # aliases, and fall through to the frequency block summary (n_total).
  has_groups <- !is.null(keys[["q6_n_HOM_REF"]]) ||
                !is.null(keys[["q6_n_HOM_STD"]]) ||
                !is.null(keys[["q6_n_total"]])

  if (!has_groups) {
    return(list(level = "NONE", reason = "No groups assigned (decomposition pending)"))
  }

  if (family_lik > 0.5 || jk_fragile) {
    if (or_passed) {
      # OR test overrides suspicion — the breakpoint evidence validates the groups
      return(list(level = "VALIDATED",
                  reason = "OR test validates despite family_likeness/jackknife concern"))
    }
    return(list(level = "SUSPECT",
                reason = paste0("family_likeness=", round(family_lik, 2),
                                if (jk_fragile) ", jackknife_fragile" else "")))
  }

  if (or_passed && (ghsl_high || jk_robust)) {
    return(list(level = "VALIDATED",
                reason = paste0("OR=", signif(or_val, 3),
                                ", p=", signif(or_p, 3),
                                if (ghsl_high) " + GHSL_concordant" else "",
                                if (jk_robust) " + jackknife_robust" else "")))
  }

  if (or_passed) {
    return(list(level = "VALIDATED",
                reason = paste0("OR=", signif(or_val, 3),
                                ", p=", signif(or_p, 3), " (OR test alone)")))
  }

  if (ghsl_high && jk_robust) {
    return(list(level = "SUPPORTED",
                reason = "GHSL concordant + jackknife robust (no OR test)"))
  }

  if (ghsl_high || jk_robust) {
    return(list(level = "SUPPORTED",
                reason = paste0(if (ghsl_high) "GHSL concordant" else "jackknife robust",
                                " (partial confirmation, no OR test)")))
  }

  return(list(level = "UNCERTAIN",
              reason = "PCA groups only, no independent confirmation"))
}


# Gate check: can this question advance beyond MEASURED?
check_group_gate <- function(question, group_val) {
  minimum <- group_validation_minimum[[question]]
  current_level <- validation_order[group_val$level]
  required_level <- validation_order[minimum]
  passes <- current_level >= required_level
  list(passes = passes, required = minimum, actual = group_val$level,
       reason = if (passes) "OK" else paste0("Need ", minimum, ", have ", group_val$level))
}
