#!/usr/bin/env Rscript
# =============================================================================
# characterize_candidate.R — Convergence-based biological characterization
# =============================================================================
#
# THREE SYSTEMS (independent):
#   System 1: DISCOVERY    → confidence tier (from compute_candidate_status.R)
#   System 2: COLLECTION   → key fill percentage (352 keys per v2 spec)
#   System 3: CHARACTER.   → THIS SCRIPT: question-level convergence
#
# For each of 7 questions, evaluates whether the collected keys CONVERGE
# on a biological conclusion. Returns ANSWERED / PARTIAL / MEASURED /
# CONTRADICTED / EMPTY per question.
#
# FIX 44 (chat 9): this script is a FUNCTION LIBRARY — it defines the
# characterize_q1..q7 functions and the master characterize_candidate()
# but has no driver. The companion driver `run_characterize.R` sources
# this file + group_validation_gate.R, reads per-candidate keys, and
# writes outputs. This file is also safe to source() standalone from an
# interactive R session for inspecting individual candidates.
#
# =============================================================================

# Helpers — MUST be defined before the gate-locator block below uses %||%.
# (FIX 46, chat 10: FIX 44's gate-locator used %||% before it was defined.
# Every caller hit "could not find function %||%" before any characterize_qN
# function could run. Also simplified sys.frame(1)$ofile path — that form
# only works when the script is sourced via source(file, chdir=TRUE) and
# crashes when sourced the normal way. The driver run_characterize.R now
# sources the gate file explicitly, so the locator is defensive-only.)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a[1]))) b else a
# FIX 47 (chat 10): safe_num() previously did `if (is.finite(x)) x else d`
# without guarding length. For a missing key, keys[[k]] -> NULL,
# as.numeric(NULL) -> numeric(0), is.finite(numeric(0)) -> logical(0),
# and `if (logical(0))` throws "argument is of length zero". Every
# characterize_qN function hit this on its first safe_num call for an
# absent key. Guard length + NULL first, then scalarize before the
# is.finite branch.
safe_num <- function(x, d = NA_real_) {
  if (is.null(x) || length(x) == 0) return(d)
  x <- suppressWarnings(as.numeric(x[1]))
  if (length(x) == 1 && is.finite(x)) x else d
}
safe_bool <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  if (is.logical(x)) return(isTRUE(x[1]))
  identical(toupper(as.character(x[1])), "TRUE")
}
has_key <- function(keys, k) {
  v <- keys[[k]]
  !is.null(v) && length(v) > 0 && !is.na(v[1]) && v[1] != ""
}

# Locate and source group_validation_gate.R (defines assess_group_validation
# and check_group_gate used below). Support both "characterize_candidate.R
# sourced directly" and "sourced via run_characterize.R" cases.
# ofile is only available when this file was sourced via source(); guard
# with tryCatch so a direct Rscript invocation doesn't blow up.
#
# NOTE (2026-04-24, pass 16): the gate moved to phase_9_classification/
# (co-located with this file). Old phase_7 paths kept as deprecated
# fallback — remove when confirmed no external caller depends on them.
.char_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.char_ofile <- if (is.null(.char_ofile) || !is.character(.char_ofile) ||
                   !nzchar(.char_ofile)) "" else .char_ofile
.gate_file_candidates <- c(
  Sys.getenv("GROUP_VALIDATION_GATE", ""),
  if (nzchar(.char_ofile))
    file.path(dirname(.char_ofile), "group_validation_gate.R")  # co-located (primary)
  else character(),
  "group_validation_gate.R",
  "phase_9_classification/group_validation_gate.R"
)
.gate_sourced <- FALSE
.gate_source_path <- NA_character_
for (.gf in .gate_file_candidates) {
  if (nzchar(.gf) && file.exists(.gf)) {
    tryCatch({ source(.gf); .gate_sourced <- TRUE; .gate_source_path <- .gf; break },
             error = function(e) NULL)
  }
}
if (!.gate_sourced &&
    (!exists("assess_group_validation", mode = "function") ||
     !exists("check_group_gate", mode = "function"))) {
  warning("[characterize_candidate] could not locate group_validation_gate.R - ",
          "set GROUP_VALIDATION_GATE env var or source the gate file before ",
          "calling characterize_candidate(). assess_group_validation + ",
          "check_group_gate must be in scope at call time.")
}

# =============================================================================
# Q1: WHAT IS IT?
# =============================================================================
characterize_q1 <- function(keys) {
  shape <- keys[["q1_shape_class"]] %||% NULL
  if (is.null(shape)) return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                                   reason = "inv_detect hasn't run (q1_shape_class missing)"))

  n_children <- safe_num(keys[["q1_n_children"]], 0)
  landscape  <- keys[["q1_landscape_category"]] %||% "unknown"
  family_lik <- safe_num(keys[["q1_family_likeness_mean"]], 0)
  d8_peel    <- safe_num(keys[["q1_d08_peel_or_hyp"]], 0.5)
  squareness <- safe_num(keys[["q1_squareness"]], 0)
  nn40       <- safe_bool(keys[["q1_nn_survives_40"]])
  inv_lik    <- safe_num(keys[["q1_inv_likeness_mean"]], 0)

  ev_for <- character(); ev_against <- character()

  # Check for contradictions first
  if (shape == "strong_square" && family_lik > 0.5) {
    ev_against <- c(ev_against, "shape=strong_square BUT family_likeness>0.5")
  }
  if (nn40 && d8_peel < 0.3) {
    ev_against <- c(ev_against, "survives_nn40 BUT peel_destroyed (suspicious)")
  }
  if (n_children > 0 && landscape == "standard") {
    ev_against <- c(ev_against, "children_detected BUT landscape=standard (inconsistent)")
  }

  if (length(ev_against) >= 2) {
    return(list(status = "CONTRADICTED", label = NA,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste("Multiple contradictions:", paste(ev_against, collapse = "; "))))
  }

  # Artifact detection
  if (family_lik > 0.5 && d8_peel < 0.3 && inv_lik < 0.4) {
    ev_for <- c("family_likeness>0.5", "peel_destroyed", "inv_likeness<0.4")
    return(list(status = "ANSWERED", label = "likely_artifact",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "Family structure signal, not inversion"))
  }

  # Simple polymorphic
  if (shape == "strong_square" && squareness > 0.7 && n_children == 0) {
    ev_for <- c(ev_for, paste0("shape=strong_square(sq=", round(squareness, 2), ")"))
    if (nn40) ev_for <- c(ev_for, "survives_nn40")
    if (d8_peel >= 0.6) ev_for <- c(ev_for, paste0("peel_stable(", round(d8_peel, 2), ")"))

    if (length(ev_for) >= 2 && length(ev_against) == 0) {
      return(list(status = "ANSWERED", label = "simple_polymorphic",
                  evidence_for = ev_for, evidence_against = ev_against,
                  reason = paste(length(ev_for), "converging lines")))
    }
    return(list(status = "PARTIAL", label = "simple_polymorphic",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "Shape suggests simple but insufficient supporting evidence"))
  }

  # Nested/complex
  if (n_children >= 1 || grepl("nested|complex", landscape)) {
    label <- if (n_children >= 2 || landscape == "complex_system") "complex_system" else "nested_system"
    ev_for <- c(paste0("n_children=", n_children), paste0("landscape=", landscape))
    return(list(status = if (length(ev_for) >= 2) "ANSWERED" else "PARTIAL",
                label = label, evidence_for = ev_for, evidence_against = ev_against,
                reason = paste("Nested/complex:", paste(ev_for, collapse = ", "))))
  }

  # Fallback
  ev_for <- c(paste0("shape=", shape))
  list(status = "MEASURED", label = NA,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = "Shape measured but doesn't clearly resolve type")
}


# =============================================================================
# Q2: WHAT'S HAPPENING INSIDE?
# =============================================================================
characterize_q2 <- function(keys) {
  n_rec <- safe_num(keys[["q2_n_recombinant"]], NA)
  if (is.na(n_rec)) return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                                 reason = "Decomposition hasn't run (q2_n_recombinant missing)"))

  stab_pct <- safe_num(keys[["q2_class_stability_pct"]], NA)
  kw_p     <- safe_num(keys[["q2_switching_kw_pval"]], 1)
  left_n   <- safe_num(keys[["q2_recomb_left_n"]], 0)
  center_n <- safe_num(keys[["q2_recomb_center_n"]], 0)
  right_n  <- safe_num(keys[["q2_recomb_right_n"]], 0)
  entropy  <- safe_num(keys[["q2_class_entropy_mean"]], NA)
  phase_ratio <- safe_num(keys[["q2_phase_block_ratio"]], NA)
  indel_conc <- safe_num(keys[["q2_indel_class_concordance"]], NA)

  # ── chat-13 wiring: pull keys from the four new structured blocks ──
  # These are ANNOTATION-ONLY for the first wiring pass; the audit at
  # §6 and the handoff §"Not in scope" note that gate promotion (the
  # new keys actually changing status/label) comes after chat-14 HPC
  # calibration. Present them as evidence_for / evidence_against lines
  # so they show up in the characterization report without altering
  # the existing decision flow.
  #
  # From regime_segments (C01j):
  q2_regime_n_segments      <- safe_num(keys[["q2_regime_n_segments"]], NA)
  q2_regime_dominant        <- keys[["q2_regime_dominant_state"]] %||% NA_character_
  q2_regime_has_recomb      <- keys[["q2_regime_has_recomb"]] %||% NA
  q2_regime_n_samples_chg   <- safe_num(keys[["q2_regime_n_samples_changed"]], NA)
  # From local_structure_segments (C01l):
  q2_bound_left             <- keys[["q2_boundary_sharpness_left"]]    %||% NA_character_
  q2_bound_right            <- keys[["q2_boundary_sharpness_right"]]   %||% NA_character_
  q2_bound_overall          <- keys[["q2_boundary_sharpness_overall"]] %||% NA_character_
  q2_core_delta12           <- safe_num(keys[["q2_local_struct_core_delta12"]],   NA)
  q2_lflk_delta12           <- safe_num(keys[["q2_local_struct_lflank_delta12"]], NA)
  q2_rflk_delta12           <- safe_num(keys[["q2_local_struct_rflank_delta12"]], NA)
  # From distance_concordance (C01m):
  q2_inv_vs_fam             <- safe_num(keys[["q2_distance_conc_inv_vs_fam_score"]], NA)
  q2_n_persistent           <- safe_num(keys[["q2_distance_conc_n_persistent"]],    NA)
  q2_n_decaying             <- safe_num(keys[["q2_distance_conc_n_decaying"]],      NA)
  # From regime_sample_dag (C01i_b):
  q2_dag_frac_r_fired       <- safe_num(keys[["q2_dag_fraction_samples_R_fired"]],   NA)
  q2_dag_frac_dev           <- safe_num(keys[["q2_dag_fraction_samples_deviating"]], NA)
  q2_dag_dominant_cohort    <- keys[["q2_dag_dominant_regime_cohort"]] %||% NA_character_
  q2_dag_dominant_support   <- safe_num(keys[["q2_dag_dominant_regime_support"]],     NA)
  q2_dag_n_r_fired          <- safe_num(keys[["q2_dag_n_samples_R_fired"]],           NA)

  ev_for <- character(); ev_against <- character()

  # Annotation-only evidence collection from the four new blocks.
  if (!is.na(q2_regime_dominant))
    ev_for <- c(ev_for, paste0("regime_dominant=", q2_regime_dominant))
  if (!is.na(q2_regime_n_segments))
    ev_for <- c(ev_for, paste0("regime_n_segments=", as.integer(q2_regime_n_segments)))
  if (isTRUE(q2_regime_has_recomb) ||
      (is.character(q2_regime_has_recomb) && q2_regime_has_recomb == "TRUE")) {
    ev_for <- c(ev_for, "regime_track_has_recomb_signal")
  }
  if (!is.na(q2_bound_overall))
    ev_for <- c(ev_for, paste0("boundary_sharpness=", q2_bound_overall))
  if (!is.na(q2_inv_vs_fam))
    ev_for <- c(ev_for, paste0("inv_vs_family_score=", signif(q2_inv_vs_fam, 3)))
  if (!is.na(q2_dag_frac_r_fired))
    ev_for <- c(ev_for, paste0("dag_frac_R_fired=", signif(q2_dag_frac_r_fired, 3)))
  if (!is.na(q2_dag_dominant_cohort))
    ev_for <- c(ev_for, paste0("dag_dominant_regime=", q2_dag_dominant_cohort,
                               if (!is.na(q2_dag_dominant_support))
                                 paste0("(support=", signif(q2_dag_dominant_support, 3), ")")
                               else ""))

  # Soft contradictions (annotation-only; don't drive CONTRADICTED yet)
  if (!is.na(q2_inv_vs_fam) && q2_inv_vs_fam < 0.2) {
    ev_against <- c(ev_against,
                    paste0("inv_vs_family_score=", signif(q2_inv_vs_fam, 3),
                           " (low — family LD signal dominates)"))
  }
  if (!is.na(q2_bound_overall) && q2_bound_overall == "diffuse") {
    ev_against <- c(ev_against, "boundary_sharpness=diffuse (no crisp edge)")
  }

  # Contradictions
  if (n_rec == 0 && !is.na(stab_pct) && stab_pct < 90) {
    ev_against <- c(ev_against, "n_recombinant=0 BUT stability<90% (decomp issue?)")
  }
  if (!is.na(phase_ratio) && phase_ratio < 1) {
    ev_against <- c(ev_against, "phase_blocks shorter inside than flanking (opposite of suppression)")
  }
  if (!is.na(indel_conc) && indel_conc < 0.3) {
    ev_against <- c(ev_against, "INDEL classes don't match PCA classes")
  }
  if (length(ev_against) >= 2) {
    return(list(status = "CONTRADICTED", label = NA,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste(ev_against, collapse = "; ")))
  }

  # Clean
  if (n_rec == 0 && !is.na(stab_pct) && stab_pct > 98) {
    ev_for <- c("n_recombinant=0", paste0("stability=", round(stab_pct, 1), "%"))
    if (!is.na(entropy) && entropy < 0.2) ev_for <- c(ev_for, "low_entropy")
    return(list(status = "ANSWERED", label = "clean",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "No recombinants, high class stability"))
  }

  # Edge conversion vs internal crossover
  if (n_rec > 0) {
    edge_n <- left_n + right_n
    ev_for <- c(paste0("n_recombinant=", n_rec))

    if (edge_n > center_n && kw_p < 0.05) {
      ev_for <- c(ev_for, paste0("edge=", edge_n, ">center=", center_n), paste0("kw_p=", signif(kw_p, 3)))
      label <- if (center_n > 0) "edge_conversion" else "edge_conversion"
      # Check if also internal crossover
      n_dco <- safe_num(keys[["q2_n_double_crossover"]], 0)
      if (n_dco > 0) {
        ev_for <- c(ev_for, paste0("double_crossover=", n_dco))
        label <- "internal_crossover"
      }
      return(list(status = "ANSWERED", label = label,
                  evidence_for = ev_for, evidence_against = ev_against,
                  reason = paste(length(ev_for), "converging lines")))
    }

    # Mosaic
    if (!is.na(stab_pct) && stab_pct < 80) {
      ev_for <- c(ev_for, paste0("stability=", round(stab_pct, 1), "%"))
      return(list(status = "ANSWERED", label = "mosaic",
                  evidence_for = ev_for, evidence_against = ev_against,
                  reason = "Low class stability + recombinants = complex mosaic"))
    }

    return(list(status = "PARTIAL", label = "edge_conversion",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "Recombinants detected but position/switching pattern unclear"))
  }

  list(status = "MEASURED", label = NA,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = "Recombination data exists but insufficient for conclusion")
}


# =============================================================================
# Q3: WHAT ARE THE BOUNDARIES DOING?
# =============================================================================
characterize_q3 <- function(keys) {
  lh <- safe_num(keys[["q3_left_hardness"]], NA)
  rh <- safe_num(keys[["q3_right_hardness"]], NA)
  if (is.na(lh) && is.na(rh))
    return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                reason = "Boundary catalog hasn't run"))

  lv <- keys[["q3_left_verdict"]] %||% "unknown"
  rv <- keys[["q3_right_verdict"]] %||% "unknown"
  l_cheats <- safe_num(keys[["q3_left_n_cheats"]], 0)
  r_cheats <- safe_num(keys[["q3_right_n_cheats"]], 0)
  l_fossil <- safe_bool(keys[["q3_left_is_fossil"]])
  r_fossil <- safe_bool(keys[["q3_right_is_fossil"]])
  dropout <- safe_num(keys[["q3_dropout_rate"]], NA)

  ev_for <- character(); ev_against <- character()

  # Contradictions
  if (!is.na(lh) && lh > 0.7 && l_cheats == 0) {
    ev_against <- c(ev_against, "left_hardness>0.7 BUT 0 cheats support it")
  }
  if (!is.na(rh) && rh > 0.7 && r_cheats == 0) {
    ev_against <- c(ev_against, "right_hardness>0.7 BUT 0 cheats support it")
  }
  if (!is.na(dropout) && dropout > 0.5) {
    ev_against <- c(ev_against, paste0("carrier_dropout=", round(dropout*100), "% (SV evidence unreliable)"))
  }

  classify_boundary <- function(hardness, verdict, n_cheats, is_fossil) {
    if (is_fossil) return("fossil")
    if (is.na(hardness)) return("unknown")
    if (hardness > 0.7 && verdict == "confirmed_structural" && n_cheats >= 3) return("hard")
    if (hardness > 0.3 && n_cheats >= 2) return("soft")
    return("eroded")
  }

  left_class <- classify_boundary(lh, lv, l_cheats, l_fossil)
  right_class <- classify_boundary(rh, rv, r_cheats, r_fossil)
  label <- paste0(left_class, "_", right_class)

  ev_for <- c(paste0("left=", left_class, "(h=", round(lh %||% 0, 2), ",cheats=", l_cheats, ")"),
              paste0("right=", right_class, "(h=", round(rh %||% 0, 2), ",cheats=", r_cheats, ")"))

  if (length(ev_against) >= 2) {
    return(list(status = "CONTRADICTED", label = label,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste(ev_against, collapse = "; ")))
  }

  status <- if (left_class != "unknown" && right_class != "unknown" &&
                l_cheats + r_cheats >= 3) "ANSWERED"
            else if (left_class != "unknown" || right_class != "unknown") "PARTIAL"
            else "MEASURED"

  list(status = status, label = label,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = paste(label, "—", l_cheats + r_cheats, "total cheats supporting"))
}


# =============================================================================
# Q4: HOW DID IT FORM?
# =============================================================================
characterize_q4 <- function(keys) {
  has_sd <- keys[["q4_has_inverted_sd"]]
  jtype_l <- keys[["q4_junction_type_left"]]
  te_enrich <- keys[["q4_te_enrichment"]]

  n_available <- sum(!is.null(has_sd), !is.null(jtype_l), !is.null(te_enrich))
  if (n_available == 0)
    return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                reason = "No mechanism analysis run (cheats 14/27/28/29 pending)"))

  has_sd <- safe_bool(has_sd)
  sd_id <- safe_num(keys[["q4_sd_identity_pct"]], 0)
  biser <- keys[["q4_biser2_concordance"]] %||% "untested"
  mh_l <- safe_num(keys[["q4_mh_length_left"]], NA)
  mh_r <- safe_num(keys[["q4_mh_length_right"]], NA)
  te_fam_same <- safe_bool(keys[["q4_te_same_family"]])
  te_orient <- keys[["q4_te_same_orientation"]] %||% "NA"

  ev_for <- character(); ev_against <- character()

  # NAHR
  if (has_sd && sd_id > 90) {
    ev_for <- c(ev_for, paste0("inverted_SD(id=", round(sd_id, 1), "%)"))
    if (biser == "agree_nahr") ev_for <- c(ev_for, "BISER2_agrees")
    if (!is.na(mh_l) && mh_l >= 8) ev_for <- c(ev_for, paste0("MH=", mh_l, "bp"))

    if (has_sd && biser == "disagree") {
      ev_against <- c(ev_against, "SD found BUT BISER2 disagrees")
    }

    conf <- if (length(ev_for) >= 3) "confirmed" else if (length(ev_for) >= 2) "likely" else "probable"
    status <- if (length(ev_against) > 0) "CONTRADICTED"
              else if (conf %in% c("confirmed", "likely")) "ANSWERED"
              else "PARTIAL"
    return(list(status = status, label = "NAHR",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste("NAHR", conf, "—", paste(ev_for, collapse = ", "))))
  }

  # TE-mediated
  if (!is.null(te_enrich) && te_enrich == "ENRICHED" && te_fam_same) {
    ev_for <- c("TE_enriched", "same_TE_family")
    if (te_orient == "inverted") ev_for <- c(ev_for, "inverted_orientation")
    conf <- if (length(ev_for) >= 3) "likely" else "probable"
    return(list(status = if (conf == "likely") "ANSWERED" else "PARTIAL",
                label = "TE_mediated",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste("TE_mediated", conf)))
  }

  # NHEJ / MMEJ (require junction data)
  if (!is.na(mh_l) && !is.na(mh_r)) {
    if (mh_l <= 2 && mh_r <= 2 && !has_sd) {
      ev_for <- c("no_SD", paste0("MH_left=", mh_l), paste0("MH_right=", mh_r))
      return(list(status = "ANSWERED", label = "NHEJ",
                  evidence_for = ev_for, evidence_against = ev_against,
                  reason = "Blunt junctions, no SD substrate"))
    }
    if (mh_l >= 3 && mh_l <= 20 && mh_r >= 3 && !has_sd) {
      ev_for <- c("no_SD", paste0("MH_left=", mh_l), paste0("MH_right=", mh_r))
      return(list(status = "ANSWERED", label = "MMEJ",
                  evidence_for = ev_for, evidence_against = ev_against,
                  reason = "Short microhomology, no SD"))
    }
  }

  if (n_available == 1) {
    return(list(status = "MEASURED", label = NA,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "Only one mechanism evidence line available"))
  }

  list(status = "PARTIAL", label = "unclassified",
       evidence_for = ev_for, evidence_against = ev_against,
       reason = "Mechanism evidence incomplete or ambiguous")
}


# =============================================================================
# Q5: HOW OLD IS IT?
# =============================================================================
characterize_q5 <- function(keys) {
  gds <- safe_num(keys[["q5_gds_gap"]], NA)
  fst <- safe_num(keys[["q5_fst_b1b3"]], NA)

  if (is.na(gds) && is.na(fst))
    return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                reason = "No age proxies computed (cheat30 + popstats pending)"))

  gds_pctl <- safe_num(keys[["q5_gds_gap_percentile"]], NA)
  rho_p    <- safe_num(keys[["q5_gds_fst_spearman_p"]], NA)
  div_shape <- keys[["q5_diversity_shape"]] %||% "unknown"
  dxy_ratio <- safe_num(keys[["q5_dxy_ratio"]], NA)
  fixed_diff <- safe_num(keys[["q5_fixed_differences"]], NA)

  ev_for <- character(); ev_against <- character()

  if (!is.na(gds) && !is.na(gds_pctl)) ev_for <- c(ev_for, paste0("GDS_gap=", round(gds, 4), "(", round(gds_pctl), "pctl)"))
  if (!is.na(fst)) ev_for <- c(ev_for, paste0("Fst=", round(fst, 4)))

  # Need ranking to classify
  if (is.na(gds_pctl)) {
    return(list(status = "MEASURED", label = NA,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "GDS/Fst measured but cross-candidate ranking not done"))
  }

  # Determine age from percentile
  if (gds_pctl >= 67) {
    age_label <- "old"
    if (div_shape %in% c("deep_U", "shallow_U")) ev_for <- c(ev_for, paste0("diversity=", div_shape))
    else if (div_shape == "flat") ev_against <- c(ev_against, "diversity=flat (expected U for old)")
    if (!is.na(fixed_diff) && fixed_diff > 0) ev_for <- c(ev_for, paste0("fixed_diffs=", fixed_diff))
    if (!is.na(fixed_diff) && fixed_diff == 0) ev_against <- c(ev_against, "zero_fixed_diffs (unexpected for old)")
  } else if (gds_pctl <= 33) {
    age_label <- "young"
    if (div_shape == "flat") ev_for <- c(ev_for, "diversity=flat")
    else if (div_shape %in% c("deep_U", "shallow_U")) ev_against <- c(ev_against, paste0("diversity=", div_shape, " (unexpected for young)"))
    if (!is.na(fixed_diff) && fixed_diff > 5) ev_against <- c(ev_against, paste0("fixed_diffs=", fixed_diff, " (unexpected for young)"))
  } else {
    age_label <- "intermediate"
  }

  # Check Fst agreement
  if (!is.na(rho_p) && rho_p > 0.05) {
    ev_against <- c(ev_against, "GDS-Fst correlation not significant (age proxies may disagree)")
  }

  # Gap 2: cheat30 GDS between/within consistency gate.
  # If between-arrangement GDS is meaningfully higher than within-arrangement
  # GDS, the gds_gap is biologically supported (not noise), which corroborates
  # the age_label regardless of direction. Near-zero or inverted separation
  # means the gds_gap percentile shouldn't be trusted as age evidence.
  gds_within_std <- safe_num(keys[["q5_gds_within_std"]], NA)
  gds_within_inv <- safe_num(keys[["q5_gds_within_inv"]], NA)
  gds_between    <- safe_num(keys[["q5_gds_between"]],    NA)
  het_pattern    <- keys[["q5_gds_het_pattern"]] %||% NA_character_
  if (!is.na(gds_between) && !is.na(gds_within_std) && !is.na(gds_within_inv)) {
    bw_gap <- gds_between - max(gds_within_std, gds_within_inv)
    if (bw_gap >= 0.05) {
      ev_for <- c(ev_for, paste0("gds_between_within_consistent(gap=", round(bw_gap, 3), ")"))
    } else {
      ev_against <- c(ev_against, paste0("gds_between≈within (gap=", round(bw_gap, 3), ", age proxy may be noise)"))
    }
  }
  if (!is.na(het_pattern) && nzchar(het_pattern) && het_pattern != "unknown") {
    ev_for <- c(ev_for, paste0("het_pattern=", het_pattern))
  }

  if (length(ev_against) >= 2) {
    return(list(status = "CONTRADICTED", label = age_label,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste("Age evidence contradictory:", paste(ev_against, collapse = "; "))))
  }

  status <- if (length(ev_for) >= 3 && length(ev_against) == 0) "ANSWERED"
            else if (length(ev_for) >= 2) "PARTIAL"
            else "MEASURED"

  list(status = status, label = age_label,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = paste(age_label, "—", length(ev_for), "supporting,", length(ev_against), "against"))
}


# =============================================================================
# Q6: HOW COMMON IS IT?
# =============================================================================
characterize_q6 <- function(keys) {
  freq <- safe_num(keys[["q6_freq_inv"]], NA)
  if (is.na(freq))
    return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                reason = "Decomposition hasn't produced frequency estimate"))

  n_total <- safe_num(keys[["q6_n_total"]], 0)
  cv <- safe_num(keys[["q6_freq_cv_across_qgroups"]], NA)
  jk_delta <- safe_num(keys[["q6_jackknife_max_delta"]], NA)
  family <- keys[["q6_family_linkage"]] %||% "unknown"

  ev_for <- c(paste0("freq=", round(freq, 3)))
  ev_against <- character()

  label <- if (freq < 0.05) "rare"
           else if (freq < 0.15) "low"
           else if (freq < 0.50) "intermediate"
           else if (freq < 0.85) "high"
           else "nearly_fixed"

  if (n_total >= 200) ev_for <- c(ev_for, paste0("n_classified=", n_total, "/226"))
  if (n_total < 200) ev_against <- c(ev_against, paste0("only ", n_total, "/226 classified"))

  if (family == "multi_family") ev_for <- c(ev_for, "multi_family")
  if (!is.na(jk_delta) && jk_delta > 0.10) {
    ev_against <- c(ev_against, paste0("jackknife_delta=", round(jk_delta, 3), " (freq fragile)"))
  }

  # Carrier concordance
  conc <- safe_num(keys[["q3_carrier_concordance"]], NA)
  if (!is.na(conc) && conc < 0.6) {
    ev_against <- c(ev_against, paste0("PCA-SV carrier concordance=", round(conc, 2), " (freq may be biased)"))
  }

  status <- if (n_total >= 200 && length(ev_against) == 0) "ANSWERED"
            else if (length(ev_against) >= 2) "CONTRADICTED"
            else "PARTIAL"

  list(status = status, label = label,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = paste(label, "(", round(freq, 3), ")"))
}


# =============================================================================
# Q7: IS IT REAL?
# =============================================================================
characterize_q7 <- function(keys) {
  layer_a <- safe_bool(keys[["q7_layer_a_detected"]])
  layer_b <- safe_bool(keys[["q7_layer_b_detected"]])
  layer_c <- safe_bool(keys[["q7_layer_c_ghsl_detected"]])
  layer_d <- safe_bool(keys[["q7_layer_d_tested"]])
  fisher_p <- safe_num(keys[["q7_layer_d_fisher_p"]], 1)
  fisher_or <- safe_num(keys[["q7_layer_d_fisher_or"]], 1)
  ghsl_qual <- keys[["q7_layer_c_ghsl_quality"]] %||% "ABSENT"
  # FIX 43 (chat 9): q7_t9_jackknife_status is aspirational (no writer).
  # q6_family_linkage carries the same semantic info reliably (written by
  # seal + overwritten by C01f jackknife), so prefer it.
  family <- keys[["q6_family_linkage"]] %||% "unknown"
  jackknife <- keys[["q7_t9_jackknife_verdict"]] %||%
               keys[["q7_t9_jackknife_status"]] %||% "unknown"
  dropout <- safe_num(keys[["q7b_observed_dropout_pct"]], NA)

  if (!layer_a && !layer_b && !layer_c)
    return(list(status = "EMPTY", label = NA, evidence_for = NULL, evidence_against = NULL,
                reason = "No detection layers have run"))

  ev_for <- character(); ev_against <- character()
  layers_pass <- 0L

  if (layer_a) { ev_for <- c(ev_for, "Layer_A(PCA)"); layers_pass <- layers_pass + 1L }
  if (layer_b) { ev_for <- c(ev_for, "Layer_B(SV)"); layers_pass <- layers_pass + 1L }
  if (layer_c && ghsl_qual %in% c("HIGH", "MODERATE")) {
    ev_for <- c(ev_for, paste0("Layer_C(GHSL_", ghsl_qual, ")")); layers_pass <- layers_pass + 1L
  }
  if (layer_d && fisher_p < 0.05 && fisher_or > 3) {
    ev_for <- c(ev_for, paste0("Layer_D(OR=", round(fisher_or, 1), ",p=", signif(fisher_p, 3), ")"))
    layers_pass <- layers_pass + 1L
  }

  if (family == "multi_family" || jackknife == "robust_multi_family") {
    ev_for <- c(ev_for, "multi_family_robust")
  }

  # Contradictions
  if (layer_a && layer_d && fisher_p > 0.2 && fisher_or < 2) {
    if (is.na(dropout) || dropout < 0.3) {
      ev_against <- c(ev_against, "PCA detects but Fisher association FAILS (not dropout — real disagreement)")
    } else {
      ev_against <- c(ev_against, paste0("Fisher weak BUT dropout=", round(dropout*100), "% (test may be compromised)"))
    }
  }
  if (layer_c && !layer_a) {
    ev_against <- c(ev_against, "GHSL signal without PCA block (unusual — investigate)")
  }

  if (length(ev_against) >= 2) {
    return(list(status = "CONTRADICTED", label = NA,
                evidence_for = ev_for, evidence_against = ev_against,
                reason = paste(ev_against, collapse = "; ")))
  }

  label <- if (layers_pass >= 3) "confirmed"
           else if (layers_pass == 2) "likely"
           else if (layers_pass == 1) "candidate"
           else "insufficient"

  if (grepl("artifact|family_structure|H1_family", keys[["q7_verdict"]] %||% "")) {
    label <- "artifact"
    ev_against <- c(ev_against, paste0("verdict=", keys[["q7_verdict"]]))
  }

  status <- if (label == "confirmed" && length(ev_against) == 0) "ANSWERED"
            else if (label %in% c("likely", "confirmed")) "PARTIAL"
            else if (label == "artifact") "ANSWERED"
            else "MEASURED"

  list(status = status, label = label,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = paste(layers_pass, "/4 layers pass"))
}


# =============================================================================
# Q_QC_SHELF: DATA-QUALITY TRIAGE (pass 13)
# =============================================================================
# Reads q_qc_shelf_* keys appended by _qc_shelf_reader.R (env-gated by
# QC_SHELF_EVIDENCE_DIR in compute_candidate_status.R). This is a SOFT flag
# — it never caps tier or group validation; it just labels the shelf's
# data quality so downstream consumers can de-prioritize messy candidates
# without dropping them from the catalog.
#
# Status semantics:
#   ANSWERED   — ran and got a definitive flag (clean / coverage_artifact /
#                high_uncertain / low_snp / messy)
#   MEASURED   — ran but came back "unknown" (thresholds couldn't decide)
#   EMPTY      — phase_5_qc_triage did not run for this candidate
# =============================================================================
characterize_q_qc_shelf <- function(keys) {
  ran <- safe_bool(keys[["q_qc_shelf_ran"]])

  if (!ran) {
    return(list(status = "EMPTY", label = NA,
                evidence_for = NULL, evidence_against = NULL,
                reason = "phase_5_qc_triage did not run for this candidate"))
  }

  flag <- keys[["q_qc_shelf_flag"]] %||% "unknown"
  fst_fold <- safe_num(keys[["q_qc_shelf_fst_enrichment_fold"]], NA)
  hovere_het <- safe_num(keys[["q_qc_shelf_hovere_het_inside"]], NA)
  snp_r <- safe_num(keys[["q_qc_shelf_snp_ratio"]], NA)
  unc_r <- safe_num(keys[["q_qc_shelf_uncertain_ratio"]], NA)
  cov_r <- safe_num(keys[["q_qc_shelf_coverage_ratio"]], NA)

  ev_for <- character()
  ev_against <- character()

  # Positive biological signals
  if (is.finite(fst_fold) && fst_fold >= 2) {
    ev_for <- c(ev_for, paste0("Fst_enrichment=", round(fst_fold, 2), "x"))
  }
  if (is.finite(hovere_het) && hovere_het >= 1.3) {
    ev_for <- c(ev_for, paste0("HoverE_HET=", round(hovere_het, 2),
                                " (heterokaryotype signature)"))
  }

  # Artifact signals
  if (is.finite(snp_r) && snp_r <= 0.5) {
    ev_against <- c(ev_against, paste0("SNP_ratio=", round(snp_r, 2),
                                        " (shelf ≤ 50% of ref)"))
  }
  if (is.finite(unc_r) && unc_r >= 2) {
    ev_against <- c(ev_against, paste0("BEAGLE_uncertain_ratio=", round(unc_r, 2),
                                        " (shelf ≥ 2× ref)"))
  }
  if (is.finite(cov_r) && cov_r <= 0.6) {
    ev_against <- c(ev_against, paste0("coverage_ratio=", round(cov_r, 2),
                                        " (shelf ≤ 60% of ref)"))
  }

  # Status + label come directly from the flag
  if (flag == "unknown") {
    return(list(status = "MEASURED", label = "unknown",
                evidence_for = ev_for, evidence_against = ev_against,
                reason = "phase_5_qc_triage ran but thresholds could not decide"))
  }

  status <- "ANSWERED"
  list(status = status, label = flag,
       evidence_for = ev_for, evidence_against = ev_against,
       reason = paste0("flag=", flag,
                       if (flag == "clean") " (shelf passes data-quality triage)"
                       else if (flag == "messy") " (multiple artifact signals)"
                       else ""))
}


# =============================================================================
# MASTER FUNCTION: characterize one candidate
# =============================================================================
characterize_candidate <- function(keys) {
  # ── Step 0: Assess group validation level ──
  group_val <- assess_group_validation(keys)

  results <- list()

  # ── Step 1: Evaluate each question WITH group gate ──
  for (q in paste0("Q", 1:7)) {
    q_func <- switch(q,
      Q1 = characterize_q1,
      Q2 = characterize_q2,
      Q3 = characterize_q3,
      Q4 = characterize_q4,
      Q5 = characterize_q5,
      Q6 = characterize_q6,
      Q7 = characterize_q7
    )

    # Run the characterization
    result <- q_func(keys)

    # Apply group gate: if groups aren't good enough, cap at PARTIAL
    gate <- check_group_gate(q, group_val)
    if (!gate$passes && result$status == "ANSWERED") {
      result$status <- "PARTIAL"
      result$reason <- paste0(result$reason,
        " [CAPPED: groups=", gate$actual, ", need=", gate$required, "]")
      if (is.null(result$evidence_against)) result$evidence_against <- character()
      result$evidence_against <- c(result$evidence_against,
        paste0("group_validation=", gate$actual, " (need ", gate$required, ")"))
    }

    results[[q]] <- result
  }

  # ── Step 1b: Q_QC_SHELF — soft triage, NO group gate ──
  # Data-quality flag is evaluated regardless of group validation. A messy
  # shelf is messy whether or not the groups are well-proposed.
  results[["Q_QC_SHELF"]] <- characterize_q_qc_shelf(keys)

  # Count statuses — ONLY over Q1–Q7 (Q_QC_SHELF is a soft flag, not one
  # of the seven scientific questions; including it would dilute the
  # completion-% interpretation).
  q17_statuses <- sapply(results[paste0("Q", 1:7)], function(r) r$status)
  n_answered <- sum(q17_statuses == "ANSWERED")
  n_partial <- sum(q17_statuses == "PARTIAL")
  n_measured <- sum(q17_statuses == "MEASURED")
  n_contradicted <- sum(q17_statuses == "CONTRADICTED")
  n_empty <- sum(q17_statuses == "EMPTY")

  # Q_QC_SHELF reported separately
  qcs <- results[["Q_QC_SHELF"]]
  qcs_note <- if (qcs$status == "EMPTY") "qc_shelf=not_run"
              else sprintf("qc_shelf=%s", qcs$label %||% "unknown")

  list(
    per_question = results,
    group_validation = group_val,
    summary = list(
      answered = n_answered,
      partial = n_partial,
      measured = n_measured,
      contradicted = n_contradicted,
      empty = n_empty,
      group_level = group_val$level,
      qc_shelf_flag = qcs$label %||% NA_character_,
      qc_shelf_status = qcs$status,
      characterization_string = paste0(
        n_answered, "/7 ANSWERED, ",
        n_partial, " PARTIAL, ",
        n_measured, " MEASURED, ",
        n_contradicted, " CONTRADICTED, ",
        n_empty, " EMPTY",
        " [groups=", group_val$level, "; ", qcs_note, "]"
      )
    )
  )
}


# =============================================================================
# FORMAT: human-readable report for one candidate
# =============================================================================
format_characterization <- function(cid, char_result) {
  lines <- c(
    paste(rep("=", 70), collapse = ""),
    paste0("CHARACTERIZATION: ", cid),
    paste(rep("=", 70), collapse = ""),
    paste0("  ", char_result$summary$characterization_string),
    ""
  )

  q_names <- c(
    Q1 = "What is it?",
    Q2 = "What's inside?",
    Q3 = "Boundaries?",
    Q4 = "How formed?",
    Q5 = "How old?",
    Q6 = "How common?",
    Q7 = "Is it real?"
  )

  for (q in paste0("Q", 1:7)) {
    r <- char_result$per_question[[q]]
    status_icon <- switch(r$status,
      ANSWERED = "\u2705",      # green check
      PARTIAL = "\u26A0\uFE0F",  # warning
      MEASURED = "\u2B1C",       # white square
      CONTRADICTED = "\u274C",   # red X
      EMPTY = "\u2B1B",         # black square
      "?"
    )
    label_str <- if (!is.null(r$label) && !is.na(r$label)) paste0(" -> ", r$label) else ""
    lines <- c(lines, sprintf("  %s  %-12s %-20s %s%s",
                               status_icon, r$status, q_names[q], r$reason, label_str))

    if (length(r$evidence_for) > 0) {
      lines <- c(lines, paste0("       FOR: ", paste(r$evidence_for, collapse = ", ")))
    }
    if (length(r$evidence_against) > 0) {
      lines <- c(lines, paste0("       AGAINST: ", paste(r$evidence_against, collapse = ", ")))
    }
  }

  lines <- c(lines, paste(rep("=", 70), collapse = ""))
  paste(lines, collapse = "\n")
}
