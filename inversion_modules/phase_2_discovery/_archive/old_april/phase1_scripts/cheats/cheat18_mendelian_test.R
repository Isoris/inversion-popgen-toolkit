#!/usr/bin/env Rscript
# =============================================================================
# cheat18_mendelian_test.R — Mendelian transmission test
#
# BIOLOGY:
#   If inversion genotypes are real, they must follow Mendelian inheritance.
#   In a hatchery with known family structure (from ngsRelate), we test:
#   do related pairs have compatible genotypes?
#   HOM_INV × HOM_REF → all offspring HET (never HOM_INV or HOM_REF)
#   Violations = wrong genotype call or false inversion.
#
# INPUT:  ngsRelate pairwise kinship, decomposition genotypes
# OUTPUT: per-candidate violation rate, enrichment in relatives
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
THETA_1ST_DEGREE    <- 0.20
THETA_DUPLICATE     <- 0.40
THETA_2ND_DEGREE    <- 0.10
MIN_PAIRS_TESTED    <- 5L
VIOLATION_PASS      <- 0.05
VIOLATION_MARGINAL  <- 0.15

# ── Load kinship pairs from ngsRelate ─────────────────────────────────

load_kinship_pairs <- function(ngsrelate_file, theta_threshold = THETA_1ST_DEGREE) {
  if (!file.exists(ngsrelate_file)) {
    message("[cheat18] ngsRelate file not found: ", ngsrelate_file)
    return(data.table())
  }
  dt <- tryCatch(fread(ngsrelate_file), error = function(e) data.table())
  if (nrow(dt) == 0) return(data.table())

  # Detect column names (ngsRelate output varies)
  if ("ida" %in% names(dt) && "idb" %in% names(dt)) {
    setnames(dt, c("ida","idb"), c("sample_a","sample_b"), skip_absent = TRUE)
  }
  theta_col <- intersect(c("theta","rab","KING"), names(dt))[1]
  if (is.na(theta_col)) {
    message("[cheat18] Cannot find kinship column"); return(data.table())
  }
  setnames(dt, theta_col, "theta", skip_absent = TRUE)

  # Classify pair type
  dt <- dt[theta > theta_threshold & theta < THETA_DUPLICATE]
  dt[, pair_type := fifelse(theta > THETA_1ST_DEGREE, "first_degree",
                             fifelse(theta > THETA_2ND_DEGREE, "second_degree", "other"))]
  dt[, .(sample_a, sample_b, theta, pair_type)]
}

# ── Check Mendelian compatibility ─────────────────────────────────────

check_mendelian_compatibility <- function(gt_a, gt_b) {
  # Impossible parent-offspring combinations
  # HOM_REF × HOM_REF → only HOM_REF possible
  # HOM_REF × HOM_INV → only HET possible
  # HOM_INV × HOM_INV → only HOM_INV possible
  # All other combos have at least one compatible outcome

  if (is.na(gt_a) || is.na(gt_b)) return(NA)

  pair <- paste(sort(c(gt_a, gt_b)), collapse = "_")
  # These pairs cannot produce certain offspring:
  impossible <- list(
    "HOM_INV_HOM_REF" = c("HOM_REF", "HOM_INV"),  # must be HET
    "HOM_REF_HOM_REF" = c("HET", "HOM_INV"),       # must be HOM_REF
    "HOM_INV_HOM_INV" = c("HET", "HOM_REF")        # must be HOM_INV
  )
  # For first-degree pairs (not knowing who is parent/offspring),
  # check if the PAIR itself is compatible
  # i.e., could one member be parent and the other offspring?
  # HOM_REF × HOM_INV → offspring must be HET → this pair is only valid if observed
  # as a pair, since both are "possible" relative states for parent-offspring

  # Simplified: check if pair genotypes are compatible as ANY relative pair
  if (pair == "HOM_INV_HOM_REF") return(FALSE)  # violation: PO must produce HET
  if (pair == "HOM_REF_HOM_REF") return(TRUE)   # compatible
  if (pair == "HOM_INV_HOM_INV") return(TRUE)   # compatible
  if (pair == "HET_HOM_REF") return(TRUE)
  if (pair == "HET_HOM_INV") return(TRUE)
  if (pair == "HET_HET") return(TRUE)
  TRUE  # default compatible
}

# Actually: for parent-offspring, the real impossibles are:
# Parent HOM_REF × Offspring HOM_INV → VIOLATION
# Parent HOM_INV × Offspring HOM_REF → VIOLATION
# Since we don't know direction, HOM_REF ↔ HOM_INV in either direction = violation

check_mendelian_po <- function(gt_a, gt_b) {
  if (is.na(gt_a) || is.na(gt_b)) return(NA)
  pair <- paste(sort(c(gt_a, gt_b)), collapse = "_")
  # The only impossible parent-offspring pair is HOM_REF × HOM_INV
  if (pair == "HOM_INV_HOM_REF") return(FALSE)
  TRUE
}

# ── Run Mendelian test on a candidate ─────────────────────────────────

run_mendelian_test <- function(kinship_pairs, decomp_genotypes,
                                candidate_id = "unknown") {
  if (nrow(kinship_pairs) == 0 || length(decomp_genotypes) == 0) {
    return(list(n_pairs_tested = 0L, n_violations = 0L,
                violation_rate = NA_real_, confidence = NA_real_,
                verdict = "NO_DATA", violation_details = data.table()))
  }

  first_deg <- kinship_pairs[pair_type == "first_degree"]
  results <- list()
  n_tested <- 0L; n_viol <- 0L

  for (i in seq_len(nrow(first_deg))) {
    sa <- first_deg$sample_a[i]
    sb <- first_deg$sample_b[i]
    ga <- decomp_genotypes[sa]
    gb <- decomp_genotypes[sb]
    if (is.na(ga) || is.na(gb)) next

    n_tested <- n_tested + 1L
    compat <- check_mendelian_po(ga, gb)
    if (!is.na(compat) && !compat) {
      n_viol <- n_viol + 1L
      results[[length(results)+1]] <- data.table(
        sample_a = sa, sample_b = sb,
        gt_a = ga, gt_b = gb,
        theta = first_deg$theta[i])
    }
  }

  viol_rate <- if (n_tested > 0) n_viol / n_tested else NA_real_
  confidence <- if (n_tested >= MIN_PAIRS_TESTED) 1 - viol_rate else NA_real_
  verdict <- if (is.na(viol_rate)) "NO_DATA"
    else if (viol_rate < VIOLATION_PASS) "PASSES"
    else if (viol_rate < VIOLATION_MARGINAL) "MARGINAL"
    else "FAILS"

  list(n_pairs_tested = n_tested, n_violations = n_viol,
       violation_rate = round(viol_rate, 4),
       confidence = round(confidence, 4),
       verdict = verdict,
       violation_details = if (length(results) > 0) rbindlist(results) else data.table())
}

# ── Enrichment in relatives ───────────────────────────────────────────

enrichment_in_relatives <- function(kinship_pairs, decomp_genotypes) {
  if (nrow(kinship_pairs) == 0 || length(decomp_genotypes) == 0)
    return(list(enrichment_ratio = NA_real_, fisher_p = NA_real_))

  gt_named <- decomp_genotypes
  first_deg <- kinship_pairs[pair_type == "first_degree"]

  # Related pairs: same genotype rate
  related_same <- 0L; related_total <- 0L
  for (i in seq_len(nrow(first_deg))) {
    ga <- gt_named[first_deg$sample_a[i]]
    gb <- gt_named[first_deg$sample_b[i]]
    if (is.na(ga) || is.na(gb)) next
    related_total <- related_total + 1L
    if (ga == gb) related_same <- related_same + 1L
  }

  # Random pairs: expected same genotype rate from frequencies
  gt_tab <- table(gt_named[!is.na(gt_named)])
  freqs  <- gt_tab / sum(gt_tab)
  random_same_rate <- sum(freqs^2)

  related_same_rate <- if (related_total > 0) related_same / related_total else NA_real_

  enrichment <- if (is.finite(related_same_rate) && random_same_rate > 0)
    round(related_same_rate / random_same_rate, 3) else NA_real_

  # Fisher exact test
  fisher_p <- tryCatch({
    n_random <- sum(gt_tab) * (sum(gt_tab) - 1) / 2  # total possible pairs
    mat <- matrix(c(related_same, related_total - related_same,
                     round(n_random * random_same_rate),
                     round(n_random * (1 - random_same_rate))),
                   nrow = 2)
    fisher.test(mat)$p.value
  }, error = function(e) NA_real_)

  list(enrichment_ratio = enrichment, fisher_p = fisher_p,
       related_same_rate = round(related_same_rate, 4),
       random_same_rate = round(random_same_rate, 4))
}

# ── Search mode ────────────────────────────────────────────────────────

search_mendelian <- function(chr, zone_start, zone_end, ...) {
  data.table(method = "mendelian_test", best_bp = NA_integer_,
             score = 0, is_precise = FALSE,
             detail = "search_not_applicable")
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat18 <- function(chr, candidate_id, decomp_genotypes,
                         ngsrelate_file = Sys.getenv("NGSRELATE_FILE")) {
  message("[cheat18] Candidate: ", candidate_id,
          " | ", length(decomp_genotypes), " genotyped samples")

  kinship <- load_kinship_pairs(ngsrelate_file)
  message("[cheat18] Loaded ", nrow(kinship), " first-degree pairs")

  mendelian <- run_mendelian_test(kinship, decomp_genotypes, candidate_id)
  message("[cheat18] Tested: ", mendelian$n_pairs_tested,
          " | Violations: ", mendelian$n_violations,
          " (", mendelian$violation_rate, ")",
          " → ", mendelian$verdict)

  enrichment <- enrichment_in_relatives(kinship, decomp_genotypes)
  message("[cheat18] Enrichment: ", enrichment$enrichment_ratio,
          " (p=", enrichment$fisher_p, ")")

  list(mendelian = mendelian, enrichment = enrichment,
       search_result = search_mendelian(chr, 0, 0))
}
