#!/usr/bin/env Rscript
# =============================================================================
# STEP_C01i_d_seal.R — Phase 4b.4 (synthesis) of v10.1 Phase 4b rewrite
# =============================================================================
# Final synthesis of Phase 4b outputs. For each candidate, reads:
#   - internal_dynamics.json            (from STEP_C01i_decompose)
#   - recombinant_map.json              (from STEP_C01i_b_multi_recomb)
#   - internal_ancestry_composition.json (from STEP_C01i_c_nested_composition)
#
# Resolves per-sample final class according to the rules in
# docs/PHASE4B_REWRITE_ARCHITECTURE.md §3:
#
#   Rule 1: recomb_status != NOT_RECOMB → FINAL = RECOMBINANT
#   Rule 2: structure_type = multi_block_fragmented AND recomb_status =
#           NOT_RECOMB → FINAL = RECOMBINANT (weak, event_class = ambiguous)
#   Rule 3: structure_type = two_block_composite → keep pca_class + flag
#           sample as in_composite_region (it's correct for whichever
#           system we're calling on, but warn downstream)
#   Rule 4: Otherwise → FINAL = pca_class
#
# Then registers the four groups (HOM_REF, HET, HOM_INV, RECOMBINANT) in
# sample_registry, plus optional GC/DCO subgroups, plus HOM_STD alias.
# Sets q6_group_validation = UNCERTAIN (capped) based on quality metrics:
#
#   composite_flag = likely_composite       → q6_group_validation = UNCERTAIN
#                                              with quality_flag =
#                                              "composite_cap"
#                                              C01f MUST NOT promote above
#                                              UNCERTAIN.
#   silhouette_score < 0.25                  → quality_flag = "low_silhouette"
#                                              C01f may promote but with caution.
#   bic_gap < 0.05                           → quality_flag = "weak_k3_fit"
#   phase_concordance < 0.30                 → quality_flag = "low_phase"
#   everything else                          → q6_group_validation = UNCERTAIN
#                                              with quality_flag = "normal"
#
# INPUTS:
#   --candidates   candidate table
#   --decomp_dir   STEP_C01i_decompose output root (for fallback read)
#   --recomb_dir   STEP_C01i_b_multi_recomb output root
#   --nested_dir   STEP_C01i_c_nested_composition output root
#   --outdir       fallback for this script's outputs
#   --tier_max     default 3
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--candidates", type = "character", help = "Candidate table"),
  make_option("--decomp_dir", type = "character", default = "decomp_out"),
  make_option("--recomb_dir", type = "character", default = "recomb_out"),
  make_option("--nested_dir", type = "character", default = "nested_comp_out"),
  make_option("--outdir",     type = "character", default = "seal_out"),
  make_option("--tier_max",   type = "integer", default = 3L)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$candidates) || !file.exists(opt$candidates)) stop("--candidates required")

# ── Helpers + registry ───────────────────────────────────────────────────────
for (p in c("lib_decompose_helpers.R", "R/lib_decompose_helpers.R")) {
  if (file.exists(p)) { source(p); break }
}
paths <- resolve_decomp_paths(list(outdir = opt$outdir))
reg <- try_load_registry()
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Read all three blocks for one candidate
# =============================================================================
read_all_blocks <- function(cid, decomp_dir, recomb_dir, nested_dir, reg) {
  # Each read_block_safe tries registry first, then JSON fallback
  list(
    dyn = read_block_safe(reg, cid, "internal_dynamics",
                            file.path(decomp_dir, cid)),
    rec = read_block_safe(reg, cid, "recombinant_map",
                            file.path(recomb_dir, cid)),
    nst = read_block_safe(reg, cid, "internal_ancestry_composition",
                            file.path(nested_dir, cid))
  )
}

# =============================================================================
# Resolve final per-sample class
# =============================================================================
resolve_final_classes <- function(blocks) {
  dyn <- blocks$dyn; rec <- blocks$rec; nst <- blocks$nst
  if (is.null(dyn)) {
    warning("[seal] cannot resolve without internal_dynamics block")
    return(NULL)
  }

  # dyn$data$per_sample is a list of per-sample records
  dyn_ps <- dyn$data$per_sample
  if (is.null(dyn_ps) || length(dyn_ps) == 0) return(NULL)

  # Index recombinants by sample_id
  recomb_by_sample <- list()
  if (!is.null(rec) && !is.null(rec$data$recombinants)) {
    for (r in rec$data$recombinants) {
      recomb_by_sample[[as.character(r$sample_id)]] <- r
    }
  }

  # Index nested_composition per-sample structure type
  nst_struct_by_sample <- list()
  if (!is.null(nst) && !is.null(nst$data$per_sample)) {
    for (p in nst$data$per_sample) {
      nst_struct_by_sample[[as.character(p$sample_id)]] <-
        as.character(p$structure_type)
    }
  }

  # Apply rules
  final_records <- list()
  for (s in dyn_ps) {
    sid <- as.character(s$sample_id)
    pca_class <- as.character(s$pca_class)

    rec_info <- recomb_by_sample[[sid]]
    struct <- nst_struct_by_sample[[sid]]

    # Rule 1: recombinant signal from multi_recomb
    is_recomb <- !is.null(rec_info)
    recomb_source <- if (is_recomb) "signal_detection" else NA_character_
    event_class <- if (is_recomb) as.character(rec_info$event_class) else NA_character_
    # chat-13: pass multi_recomb's canonical recomb_subgroup through to
    # register_all_groups. The full split (RECOMBINANT / RECOMBINANT_GC /
    # RECOMBINANT_DCO) lives on the multi_recomb per-sample record.
    recomb_subgroup <- if (is_recomb) {
      as.character(rec_info$recomb_subgroup %||% NA_character_)
    } else NA_character_

    # Rule 2: structure type says fragmented but no recomb detected → flag as weak
    if (!is_recomb && !is.null(struct) && struct == "multi_block_fragmented") {
      is_recomb <- TRUE
      recomb_source <- "nested_comp_fragmented"
      event_class <- "ambiguous"
      recomb_subgroup <- NA_character_
    }

    # Rule 3: two_block_composite flag (keep pca_class, mark sample)
    in_composite <- !is.null(struct) && struct == "two_block_composite"

    final_class <- if (is_recomb) "RECOMBINANT" else pca_class

    final_records[[length(final_records) + 1L]] <- list(
      sample_id = sid,
      pca_class = pca_class,
      structure_type = struct %||% NA_character_,
      # chat-13 Finding BI: C01i_b emits event_class values
      # "gene_conversion_embedded", "double_crossover",
      # "single_crossover_or_ambiguous". Accept both the current-canonical
      # names and the short legacy names ("gene_conversion",
      # "double_crossover", "suspicious"/"ambiguous"). Map to a stable
      # recomb_status token.
      recomb_status = if (is_recomb) {
        ec_norm <- event_class %||% "ambiguous"
        if      (ec_norm %in% c("gene_conversion_embedded", "gene_conversion"))
          "recomb_GC"
        else if (ec_norm == "double_crossover")
          "recomb_DCO"
        else if (ec_norm %in% c("single_crossover_or_ambiguous",
                                  "suspicious", "ambiguous", NA_character_))
          "recomb_ambiguous"
        else
          "recomb_ambiguous"
      } else "NOT_RECOMB",
      recomb_event_class = event_class %||% NA_character_,
      recomb_subgroup = recomb_subgroup %||% NA_character_,
      recomb_source = recomb_source %||% NA_character_,
      in_composite_region = in_composite,
      final_class = final_class,
      # Diagnostics passed through
      window_consistency = s$window_consistency %||% NA_real_,
      sv_prior_seeded = s$sv_prior_seeded %||% FALSE,
      discordant = s$discordant %||% FALSE,
      cheat2_constrained = s$cheat2_constrained %||% FALSE
    )
  }

  final_records
}

# =============================================================================
# Determine quality flag and validation level
# =============================================================================
determine_validation <- function(blocks, n_samples_total) {
  dyn <- blocks$dyn; nst <- blocks$nst

  silhouette <- dyn$data$silhouette_score %||% NA_real_
  bic_gap    <- dyn$data$bic_gap_k3_vs_k2 %||% NA_real_
  phase_conc <- dyn$data$phase_concordance %||% NA_real_

  composite_flag <- if (!is.null(nst)) {
    nst$data$composite_flag %||% "unknown_no_engine_b"
  } else "unknown_no_engine_b"

  # Initial level: always UNCERTAIN after Phase 4b (C01f can promote later)
  level <- "UNCERTAIN"
  quality_flags <- character(0)
  promotion_cap <- NA_character_  # if set, downstream cannot exceed this

  if (composite_flag == "likely_composite") {
    quality_flags <- c(quality_flags, "composite_cap")
    promotion_cap <- "UNCERTAIN"
  }
  if (composite_flag == "maybe_composite") {
    quality_flags <- c(quality_flags, "maybe_composite_flagged")
  }
  if (!is.na(silhouette) && silhouette < 0.25) {
    quality_flags <- c(quality_flags, "low_silhouette")
  }
  if (!is.na(bic_gap) && bic_gap < 0.05) {
    quality_flags <- c(quality_flags, "weak_k3_fit")
  }
  if (!is.na(phase_conc) && phase_conc < 0.30) {
    quality_flags <- c(quality_flags, "low_phase")
  }
  if (length(quality_flags) == 0) quality_flags <- "normal"

  list(
    level = level,
    quality_flags = quality_flags,
    promotion_cap = promotion_cap,
    silhouette = silhouette,
    bic_gap = bic_gap,
    phase_concordance = phase_conc,
    composite_flag = composite_flag
  )
}

# =============================================================================
# Register the four groups (+ subgroups + alias)
# =============================================================================
register_all_groups <- function(cid, chrom, final_records) {
  if (is.null(reg) || is.null(reg$samples$add_group)) {
    message("[seal] registry unavailable — groups not persisted")
    return(list(registered = 0L))
  }

  class_groups <- list(
    HOM_REF     = character(0),
    HET         = character(0),
    HOM_INV     = character(0),
    RECOMBINANT = character(0),
    RECOMBINANT_GC  = character(0),
    RECOMBINANT_DCO = character(0)
  )
  for (r in final_records) {
    fc <- r$final_class
    class_groups[[fc]] <- c(class_groups[[fc]], r$sample_id)
    if (fc == "RECOMBINANT") {
      # chat-13 wiring: prefer the explicit recomb_subgroup from
      # multi_recomb's per-sample record if present. Fall back to
      # deriving from event_class (Finding BI: accept both canonical
      # "gene_conversion_embedded" and legacy "gene_conversion"; same
      # for "double_crossover").
      sg <- r$recomb_subgroup %||% NA_character_
      ec <- r$recomb_event_class %||% "ambiguous"
      target <- if (!is.na(sg) && sg %in% c("RECOMBINANT_GC",
                                              "RECOMBINANT_DCO")) {
        sg
      } else if (ec %in% c("gene_conversion_embedded", "gene_conversion")) {
        "RECOMBINANT_GC"
      } else if (ec == "double_crossover") {
        "RECOMBINANT_DCO"
      } else NA_character_
      if (!is.na(target)) {
        class_groups[[target]] <- c(class_groups[[target]], r$sample_id)
      }
    }
  }

  n_reg <- 0L

  # Core four
  for (cls in c("HOM_REF", "HET", "HOM_INV", "RECOMBINANT")) {
    samps <- class_groups[[cls]]
    if (length(samps) < 2L) next
    grp_id <- paste0("inv_", cid, "_", cls)
    ok <- tryCatch({
      reg$samples$add_group(
        group_id = grp_id, sample_ids = samps,
        chrom = chrom, inv_id = cid, subgroup = cls,
        description = paste0("C01i_d_seal: ", cls, " (n=", length(samps), ")"),
        overwrite = TRUE
      )
      TRUE
    }, error = function(e) { message("[seal] add_group fail ", grp_id, ": ", e$message); FALSE })
    if (ok) n_reg <- n_reg + 1L
  }

  # HOM_STD alias for backward compat
  if (length(class_groups$HOM_REF) >= 2L) {
    tryCatch(reg$samples$add_group(
      group_id = paste0("inv_", cid, "_HOM_STD"),
      sample_ids = class_groups$HOM_REF,
      chrom = chrom, inv_id = cid, subgroup = "HOM_STD",
      description = paste0("ALIAS of inv_", cid, "_HOM_REF (v9 migration)"),
      overwrite = TRUE
    ), error = function(e) NULL)
  }

  # Subgroups
  for (sub in c("RECOMBINANT_GC", "RECOMBINANT_DCO")) {
    samps <- class_groups[[sub]]
    if (length(samps) < 3L) next
    grp_id <- paste0("inv_", cid, "_", sub)
    tryCatch({
      reg$samples$add_group(
        group_id = grp_id, sample_ids = samps,
        chrom = chrom, inv_id = cid, subgroup = sub,
        description = paste0("cheat24: ", sub, " (n=", length(samps), ")"),
        overwrite = TRUE
      )
      n_reg <- n_reg + 1L
    }, error = function(e) NULL)
  }

  list(registered = n_reg, class_counts = lapply(class_groups, length))
}

# =============================================================================
# Main loop
# =============================================================================
cand_dt <- fread(opt$candidates)
if (!"candidate_id" %in% names(cand_dt)) {
  cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
}
if (!"tier" %in% names(cand_dt)) cand_dt[, tier := 2L]
cand_dt <- cand_dt[tier <= opt$tier_max]

cat("[seal] processing ", nrow(cand_dt), " candidates\n", sep = "")

seal_summary <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cd <- cand_dt[ci]
  cid <- cd$candidate_id
  chr <- cd$chrom

  cat("\n[seal] ─── ", cid, " ───\n", sep = "")

  blocks <- read_all_blocks(cid, opt$decomp_dir, opt$recomb_dir,
                              opt$nested_dir, reg)
  if (is.null(blocks$dyn)) {
    cat("[seal]   SKIP: internal_dynamics block not found\n")
    next
  }

  # Resolve per-sample final classes
  final_records <- resolve_final_classes(blocks)
  if (is.null(final_records) || length(final_records) == 0) {
    cat("[seal]   SKIP: no samples to resolve\n")
    next
  }

  class_counts <- table(vapply(final_records, function(r) r$final_class,
                                  character(1)))
  n_ref <- as.integer(class_counts["HOM_REF"] %||% 0)
  n_het <- as.integer(class_counts["HET"] %||% 0)
  n_inv <- as.integer(class_counts["HOM_INV"] %||% 0)
  n_rec <- as.integer(class_counts["RECOMBINANT"] %||% 0)
  n_tot <- length(final_records)
  freq_inv <- if (n_tot > 0) round((2 * n_inv + n_het) / (2 * n_tot), 4) else NA_real_

  # Determine validation level + quality flags
  valid <- determine_validation(blocks, n_tot)

  cat("[seal]   classes: REF=", n_ref, " HET=", n_het, " INV=", n_inv,
      " REC=", n_rec, " | composite=", valid$composite_flag,
      " | quality_flags=", paste(valid$quality_flags, collapse = ","), "\n",
      sep = "")

  # Register groups
  reg_result <- register_all_groups(cid, chr, final_records)
  cat("[seal]   registered ", reg_result$registered, " groups\n", sep = "")

  # Write validation/composite keys to registry (as flat evidence)
  if (!is.null(reg) && !is.null(reg$evidence$add_evidence)) {
    tryCatch(reg$evidence$add_evidence(cid, "q6_group_validation",
                                         value = valid$level,
                                         script = "STEP_C01i_d_seal.R"),
             error = function(e) NULL)
    # FIX 40 (chat 9): q1_composite_flag is written exclusively by the
    # internal_ancestry_composition.schema.json keys_extracted rule (from
    # the Python nested_composition block's `composite_flag` field). The
    # previous explicit add_evidence call here duplicated that write under
    # a different key name (schema wrote q1_ancestry_composite_flag; seal
    # wrote q1_composite_flag — two keys for the same value). Schema key
    # renamed to q1_composite_flag in the same fix; seal's explicit write
    # removed. One canonical key, one canonical writer.
    # See AUDIT_LOG_chat9 Finding X.
    tryCatch(reg$evidence$add_evidence(cid, "q2_decomp_quality_flags",
                                         value = paste(valid$quality_flags, collapse = ","),
                                         script = "STEP_C01i_d_seal.R"),
             error = function(e) NULL)
    if (!is.na(valid$promotion_cap)) {
      tryCatch(reg$evidence$add_evidence(cid, "q6_validation_promotion_cap",
                                           value = valid$promotion_cap,
                                           script = "STEP_C01i_d_seal.R"),
               error = function(e) NULL)
    }
    # ── family_linkage / polymorphism_class placeholders:
    #    REMOVED 2026-04-24. The register_C01i_frequency_block call below
    #    now writes both fields inside the frequency block, and the v3
    #    schema's keys_extracted emits q6_family_linkage and
    #    q6_polymorphism_class as flat keys. The old add_evidence placeholder
    #    writes were redundant duplicates (last-write-wins, so harmless but
    #    unnecessary). C01f later overwrites the two fields via
    #    update_C01f_frequency_linkage after the jackknife runs.

    # ── NEW: frequency block write (Option A architecture, v3 schema) ──
    #    register_C01i_frequency_block lives in utils/registry_key_helpers.R.
    #    Source the helpers if not already loaded (same guard pattern used
    #    by STEP_C01d/STEP_C01f). The helper computes HWE locally from
    #    class counts, derives freq_class from freq_inv, and writes the
    #    frequency block with placeholder family_linkage/polymorphism_class.
    #    C01f later overwrites those two fields via update_C01f_frequency_linkage.
    for (.hf in c("utils/registry_key_helpers.R", "../utils/registry_key_helpers.R",
                  file.path(Sys.getenv("BASE", ""),
                            "inversion-popgen-toolkit/utils/registry_key_helpers.R"))) {
      if (file.exists(.hf)) { source(.hf); break }
    }
    if (exists("register_C01i_frequency_block", mode = "function")) {
      tryCatch(
        register_C01i_frequency_block(
          list(n_HOM_REF = n_ref, n_HET = n_het, n_HOM_INV = n_inv,
               n_RECOMBINANT = n_rec, n_total = n_tot, freq_inv = freq_inv),
          cid, opt$outdir
        ),
        error = function(e) {
          message("[seal]   frequency block write failed for ", cid, ": ",
                  conditionMessage(e))
        }
      )
    }
  }

  # Record for summary
  seal_summary[[length(seal_summary) + 1L]] <- list(
    candidate_id = cid, chrom = chr,
    n_HOM_REF = n_ref, n_HET = n_het, n_HOM_INV = n_inv,
    n_RECOMBINANT = n_rec, n_total = n_tot, freq_inv = freq_inv,
    composite_flag = valid$composite_flag,
    quality_flags = paste(valid$quality_flags, collapse = ","),
    promotion_cap = valid$promotion_cap %||% NA_character_,
    family_linkage = "unknown",             # placeholder; C01f overwrites
    polymorphism_class = "unclassified",    # placeholder; C01f overwrites
    silhouette = valid$silhouette,
    bic_gap = valid$bic_gap,
    phase_concordance = valid$phase_concordance,
    n_groups_registered = reg_result$registered
  )
}

# Write summary table
if (length(seal_summary) > 0) {
  summ_dt <- rbindlist(seal_summary, fill = TRUE)
  summ_path <- file.path(opt$outdir, "seal_summary.tsv")
  fwrite(summ_dt, summ_path, sep = "\t")
  cat("\n[seal] summary written to ", summ_path, "\n", sep = "")
  cat("[seal] composite flag breakdown:\n")
  print(table(summ_dt$composite_flag, useNA = "ifany"))
}

cat("\n[seal] done\n")
