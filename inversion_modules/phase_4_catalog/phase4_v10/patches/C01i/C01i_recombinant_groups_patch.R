# =============================================================================
# PATCH: STEP_C01i_decomposition — add RECOMBINANT group registration
# =============================================================================
# Applies to: STEP_C01i_decomposition_rewired_24_v934_registry.R (v9.3.4)
# Target:     the group-registration block at lines ~356-366
# Purpose:    register the fourth genotype class (RECOMBINANT) and, where
#             cheat24 resolves them, the two subclasses (GC, DCO). Also sets
#             q6_group_validation = UNCERTAIN so downstream can gate on it.
#
# Design:     HOM_STD and HOM_REF are aliased for a migration window. The
#             handoff uses HOM_REF; v9.3.4 code uses HOM_STD. This patch
#             writes BOTH group aliases so either name works during migration.
#             Once all consumers read the canonical HOM_REF name, remove
#             the HOM_STD alias. (See DEPRECATION note at the end.)
#
# Apply with: str_replace in C01i, replacing the inner FOR loop that currently
#             registers only HOM_STD / HET / HOM_INV.
# =============================================================================


# ─── OLD BLOCK (what's there now, C01i v9.3.4 lines ~357-366) ────────────────
#
#         # Register genotype sample groups (for downstream dispatcher group resolution)
#         for (cls in c("HOM_STD", "HET", "HOM_INV")) {
#           samps <- result[get(id_col[1]) == r[[id_col[1]]] &
#                                get(class_col[1]) == cls]$sample_id
#           if (length(samps) > 0) {
#             reg$add_group(paste0("inv_", cid, "_", cls), samps,
#                            chrom = r$chrom, inv_id = cid, subgroup = cls,
#                            description = paste0("C01i decomposition: ", cls))
#           }
#         }

# ─── NEW BLOCK — replace with this ───────────────────────────────────────────

        # Register four genotype groups: HOM_REF, HET, HOM_INV, RECOMBINANT
        # NOTE: `HOM_STD` is kept as an ALIAS for `HOM_REF` during the v9→v10
        # migration. Both group ids are registered. Downstream scripts should
        # prefer `HOM_REF` per the handoff naming convention. Fallback to
        # `HOM_STD` if an older consumer is still around.
        register_four_classes <- function(cid, chrom, decomp_sub) {
          # decomp_sub: data.table for this candidate, columns include
          # sample_id and the class column (HOM_STD|HOM_REF|HET|HOM_INV|Recombinant)
          cls_col <- class_col[1]

          # Harmonize class labels: accept either HOM_STD or HOM_REF from decomp
          # (different code paths in C01i have used both)
          labels_present <- unique(decomp_sub[[cls_col]])
          ref_label <- intersect(c("HOM_REF", "HOM_STD"), labels_present)[1]
          if (is.na(ref_label)) ref_label <- "HOM_REF"  # default if none present

          # Canonical four groups
          class_map <- list(
            HOM_REF     = decomp_sub[get(cls_col) == ref_label]$sample_id,
            HET         = decomp_sub[get(cls_col) == "HET"]$sample_id,
            HOM_INV     = decomp_sub[get(cls_col) == "HOM_INV"]$sample_id,
            RECOMBINANT = decomp_sub[get(cls_col) == "Recombinant"]$sample_id
          )

          n_registered <- 0L
          for (cls in names(class_map)) {
            samps <- class_map[[cls]]
            if (length(samps) == 0) next
            grp_id <- paste0("inv_", cid, "_", cls)
            reg$add_group(grp_id, samps,
                           chrom = chrom, inv_id = cid, subgroup = cls,
                           description = paste0("C01i decomposition: ", cls,
                                                 " (n=", length(samps), ")"))
            n_registered <- n_registered + 1L

            # Write HOM_STD alias alongside HOM_REF for backward compatibility
            # during the migration window. Remove after v10 consumers stabilize.
            if (cls == "HOM_REF") {
              alias_id <- paste0("inv_", cid, "_HOM_STD")
              reg$add_group(alias_id, samps,
                             chrom = chrom, inv_id = cid, subgroup = "HOM_STD",
                             description = paste0("ALIAS of inv_", cid,
                                                   "_HOM_REF (v9→v10 migration)"))
            }
          }

          # Register cheat24 subgroups if available — these split RECOMBINANT
          # into gene-conversion (GC) vs double-crossover (DCO) events. Useful
          # downstream for mechanism Q4 (GC has different junction signature
          # than DCO) and for excluding DCO samples from per-window class
          # purity checks in Q2 internal dynamics.
          if ("recomb_event_class" %in% names(decomp_sub)) {
            gc_samps  <- decomp_sub[get(cls_col) == "Recombinant" &
                                      recomb_event_class == "gene_conversion"]$sample_id
            dco_samps <- decomp_sub[get(cls_col) == "Recombinant" &
                                      recomb_event_class == "double_crossover"]$sample_id
            if (length(gc_samps) > 0) {
              reg$add_group(paste0("inv_", cid, "_RECOMBINANT_GC"), gc_samps,
                             chrom = chrom, inv_id = cid,
                             subgroup = "RECOMBINANT_GC",
                             description = "cheat24: gene-conversion recombinants")
              n_registered <- n_registered + 1L
            }
            if (length(dco_samps) > 0) {
              reg$add_group(paste0("inv_", cid, "_RECOMBINANT_DCO"), dco_samps,
                             chrom = chrom, inv_id = cid,
                             subgroup = "RECOMBINANT_DCO",
                             description = "cheat24: double-crossover recombinants")
              n_registered <- n_registered + 1L
            }
          }

          # Set the per-candidate group validation level to UNCERTAIN.
          # C01f will promote this to SUPPORTED / VALIDATED or demote to
          # SUSPECT based on T8 Clair3, T9 jackknife, and Layer D OR test.
          tryCatch({
            reg$add_evidence(cid, "q6_group_validation", value = "UNCERTAIN",
                              file_path = "", script = "C01i")
          }, error = function(e) {
            message("[C01i] could not set q6_group_validation: ", conditionMessage(e))
          })

          n_registered
        }

        # Slice this candidate's rows and register its groups
        decomp_sub <- result[get(id_col[1]) == r[[id_col[1]]]]
        n_grps <- register_four_classes(cid, r$chrom, decomp_sub)
        if (n_grps > 0) {
          message("[C01i] ", cid, ": registered ", n_grps, " groups + validation=UNCERTAIN")
        }


# ─── NOTES ────────────────────────────────────────────────────────────────────
#
# 1) The `recomb_event_class` column must be added to the C01i result table
#    by the cheat24 block that already runs at lines ~230-280. Currently
#    cheat24 produces `prior_dt` as a separate data.table. Modify the
#    cheat24 block to merge `event_class` back onto `result` per sample,
#    then this patch will pick it up automatically.
#
#    The minimal merge at end of cheat24 block:
#        result <- merge(result, prior_dt[, .(sample_id, candidate_id,
#                                             recomb_event_class = event_class)],
#                         by = c("sample_id", "candidate_id"), all.x = TRUE)
#
# 2) reg$add_group already has an `overwrite = FALSE` default that skips
#    existing groups. For re-runs of C01i on the same candidate, pass
#    `overwrite = TRUE` (add as argument to register_four_classes and
#    thread through). Or call reg$rm_group first. For now the default
#    skip-and-warn behavior is fine for the first run.
#
# 3) DEPRECATION: once all downstream consumers (C01f, cheats, characterize,
#    classify) read `HOM_REF` instead of `HOM_STD`, remove the alias block
#    above. Track usage with:
#        grep -rn 'HOM_STD' inversion-popgen-toolkit/
#    and migrate call sites one at a time. The alias is cheap (two rows
#    in sample_groups.tsv per candidate) but the ambiguity is confusing,
#    so drop it once safe.
