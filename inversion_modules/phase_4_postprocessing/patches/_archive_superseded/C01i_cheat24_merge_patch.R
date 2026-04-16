# =============================================================================
# PATCH: STEP_C01i_decomposition — merge cheat24 event_class into `result`
# =============================================================================
# Applies to: STEP_C01i_decomposition_rewired_24_v934_registry.R
# Purpose:    C01i currently runs cheat24_recombinant_prior.R and produces
#             `prior_dt` as a STANDALONE summary, but does not merge the
#             per-sample event_class back onto the main `result` table.
#             Without this merge, the RECOMBINANT_GC / RECOMBINANT_DCO
#             subgroups that the recombinant-groups patch registers will
#             never fire (the column `recomb_event_class` doesn't exist
#             on `result`, so the subset returns zero samples silently).
#
# This patch is TINY (one merge) but essential. Without it, the
# recombinant-classification work that cheat24 does is discarded.
# =============================================================================

# ─── LOCATE ─────────────────────────────────────────────────────────────────
# In v9.3.4 C01i, cheat24 is sourced and invoked in a block that starts near
# line ~240. The block writes `prior_dt` (a data.table with one row per
# recombinant-sample per candidate) and, at the end, typically does:
#
#     if (nrow(prior_dt) > 0) {
#       fwrite(prior_dt, file.path(DECOMP_DIR, "cheat24_recombinant_prior.tsv.gz"))
#       message("[decomp] Cheat 24: ", nrow(prior_dt), " recombinants classified")
#     }
#
# Insert the merge IMMEDIATELY AFTER that fwrite, inside the same block.


# ─── INSERT THIS ────────────────────────────────────────────────────────────
    # Merge per-sample event_class back into `result` so downstream group
    # registration (see C01i_recombinant_groups_patch.R) can subset on it.
    # prior_dt must contain columns: sample_id, candidate_id, event_class
    # (the cheat24 output uses `event_class` directly per classify_recombinant_event)
    if (exists("prior_dt") && is.data.table(prior_dt) && nrow(prior_dt) > 0 &&
        all(c("sample_id", "candidate_id", "event_class") %in% names(prior_dt))) {
      # If result already has a recomb_event_class column from a previous
      # run, drop it first to avoid the .x / .y suffix gotcha
      if ("recomb_event_class" %in% names(result)) {
        result[, recomb_event_class := NULL]
      }
      result <- merge(
        result,
        prior_dt[, .(sample_id, candidate_id,
                       recomb_event_class = event_class)],
        by = c("sample_id", "candidate_id"),
        all.x = TRUE  # keep every row in result; recomb_event_class is NA
                      # for samples that are not recombinants
      )
      n_rec_classified <- sum(!is.na(result$recomb_event_class))
      message("[decomp] Cheat 24: merged event_class onto result (",
              n_rec_classified, " recombinant rows annotated)")
    } else {
      message("[decomp] Cheat 24: no prior_dt or wrong columns — skipping merge. ",
              "RECOMBINANT_GC/DCO subgroups will not be registered.")
    }


# ─── NOTES ──────────────────────────────────────────────────────────────────
#
# 1) `event_class` values from cheat24 are:
#    "gene_conversion", "double_crossover", "suspicious"
#
#    The recombinant-groups patch filters on these exact strings. Keep them
#    aligned. If cheat24 ever renames a class, update both sides.
#
# 2) This merge is idempotent — the drop-if-exists guard at the top lets
#    C01i be re-run on the same candidate without accumulating duplicate
#    columns.
#
# 3) If cheat24 ever produces multiple rows per (sample_id, candidate_id)
#    pair (e.g., two recombination events in one sample), this merge will
#    duplicate result rows. Current cheat24 writes one row per recombinant-
#    sample-candidate, but add a safety check if that changes:
#
#       stopifnot(!anyDuplicated(prior_dt, by = c("sample_id", "candidate_id")))
