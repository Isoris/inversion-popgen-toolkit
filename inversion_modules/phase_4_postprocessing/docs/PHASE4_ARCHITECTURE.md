# =============================================================================
# PHASE 4 ARCHITECTURE — v10 (Registry-as-Catalog + Group Validation Gates)
# =============================================================================
# Date: 2026-04-16
# Supersedes: ad-hoc C01d→C01g→C01f→C01i ordering in v9.3.4
# Scope: Phase 4 of the inversion pipeline (catalog birth + group proposal +
#        group validation + group-dependent evidence + final re-score)
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0. Rule summary
# ─────────────────────────────────────────────────────────────────────────────
#
# A candidate's evidence is produced in a phased way because different
# evidence streams require different "group validation levels" to be
# meaningful. C01i PROPOSES groups from PCA structure. C01f VALIDATES
# them against independent evidence (Clair3 GHSL, Fisher OR, jackknife,
# theta). Downstream cheats CONSUME validated groups.
#
# Single ordering rule (the ONLY ordering rule):
#
#   If a script writes a Tier-2 block that has `group_validation_required`
#   in its schema, it runs AFTER the script(s) that produce that validation
#   level. Everything else runs as early as possible, in parallel.
#
# This replaces the C01d→C01g→C01f→C01i chain with five phases that the
# orchestrator dispatches as SLURM array jobs with minimal dependencies.


# ─────────────────────────────────────────────────────────────────────────────
# 1. Group semantics — the corrected model
# ─────────────────────────────────────────────────────────────────────────────
#
# A genotype group for a candidate inversion is a set of samples that share
# one arrangement state at that locus. For a polymorphic inversion there
# are four possible states, not three:
#
#   HOM_REF     common-arrangement homozygote
#   HET         heterozygote
#   HOM_INV     rare-arrangement homozygote
#   RECOMBINANT class-switching sample (position-dependent)
#
# Why RECOMBINANT is a fourth class, not a subclass of HET:
#
#   A recombinant sample carries a chimeric haplotype: part of the
#   inversion is inherited from the HOM_REF parent, part from HOM_INV.
#   Averaged across the interval it looks HET, but per-window it is
#   HOM_REF on one side of the switchpoint and HOM_INV on the other
#   (gene conversion) or HOM_REF on both ends and HOM_INV in the middle
#   (double crossover). Treating these samples as HET in Fst/dXY/theta
#   computations pollutes the group estimates.
#
# cheat24 already classifies each recombinant into:
#   gene_conversion   — short mosaic tract near a boundary
#   double_crossover  — long mosaic tract in the interior
#   suspicious        — ambiguous (possible genotyping error)


# ─────────────────────────────────────────────────────────────────────────────
# 2. Group registration — the four registered groups per candidate
# ─────────────────────────────────────────────────────────────────────────────
#
# For each candidate `cid = <chrom>_<interval_id>`, C01i registers FOUR
# groups into sample_registry (was previously three):
#
#   inv_<cid>_HOM_REF       n samples in stripe 1
#   inv_<cid>_HET           n samples in stripe 2 (true heterozygotes only)
#   inv_<cid>_HOM_INV       n samples in stripe 3
#   inv_<cid>_RECOMBINANT   n samples flagged by cheat24 as class-switching
#
# Plus two optional sub-groups when cheat24 resolves the event class
# (useful for mechanism Q4 and frequency Q6 downstream):
#
#   inv_<cid>_RECOMBINANT_GC   gene-conversion recombinants
#   inv_<cid>_RECOMBINANT_DCO  double-crossover recombinants
#
# The `subgroup` field of sample_groups.tsv carries the class label so
# resolve_groups() can retrieve any of them by name.
#
# What the flat table stores vs. what the evidence blocks store:
#
#   sample_groups.tsv                        one row per group per candidate,
#                                            DOMINANT class across the interval
#
#   evidence_registry/<cid>/structured/      per-window class, switchpoints,
#     internal_dynamics.json                 event_class, posterior, mosaic
#                                            length per recombinant sample
#
# So: flat groups give you "who are the 7 recombinants in candidate
# LG12_17?", and the structured block gives you "for sample CGA045 in
# LG12_17, which windows are HOM_REF and which are HOM_INV, where is
# the switchpoint, and is this a GC or DCO event?".


# ─────────────────────────────────────────────────────────────────────────────
# 3. Group validation levels — promoted and demoted by independent evidence
# ─────────────────────────────────────────────────────────────────────────────
#
# Each registered group carries a validation level stored as evidence key
# `q6_group_validation` on the candidate (one level per candidate, applies
# to all four of its groups since they are jointly validated or not).
#
#   NONE        no groups assigned yet (before C01i runs)
#   UNCERTAIN   PCA-clustered groups, no independent confirmation
#               (set by C01i after decomposition)
#   SUPPORTED   GHSL concordance high + jackknife robust OR Clair3 concordant
#               (promoted by C01f Layer C or T9+T8)
#   VALIDATED   Fisher OR test passed p<0.05 on breakpoint read evidence
#               (promoted by C01f Layer D)
#   SUSPECT     jackknife fragile (family drives the signal) or
#               Clair3 strongly discordant (promoted-then-demoted)
#
# Transitions (monotone except SUSPECT which is a demotion):
#
#   NONE ──C01i──→ UNCERTAIN ──┬─C01f T8/T9──→ SUPPORTED ──C01f T4(LayerD)──→ VALIDATED
#                              │
#                              └─C01f T9 fragile──→ SUSPECT
#
# A candidate with SUSPECT groups is NOT removed from the catalog — it
# gets q7_verdict = "family_artifact" and its Q5/Q6 answers stay EMPTY.
# Downstream cheats check group validation before computing group-based
# statistics.


# ─────────────────────────────────────────────────────────────────────────────
# 4. Per-question group requirements (from handoff section 4)
# ─────────────────────────────────────────────────────────────────────────────
#
# Q1 architecture         NONE        uses sim_mat, not groups
# Q2 internal dynamics    UNCERTAIN   recombinant detection needs groups
# Q3 boundaries           NONE        boundary evidence is group-independent
# Q4 mechanism            NONE        junction sequence doesn't use groups
# Q5 age (GDS sub-block)  NONE        GDS works on dosage, no groups needed
# Q5 age (Fst sub-block)  SUPPORTED   Fst needs reliable group separation
# Q6 frequency            SUPPORTED   direct from group counts
# Q6 burden regression    VALIDATED   requires confirmed group assignments
# Q7 existence            NONE        Layer D IS the validation test
#
# This lookup is the authoritative scheduling constraint. A characterization
# script asking "can I answer Q5?" checks `q6_group_validation ≥ SUPPORTED`.


# ─────────────────────────────────────────────────────────────────────────────
# 5. The five phases
# ─────────────────────────────────────────────────────────────────────────────
#
# Phase 4a — GROUP-INDEPENDENT EVIDENCE (parallel, no ordering dependency)
#   ▪ C01d pass-1 scoring (D1..D7, D9, D10, D12; leaves D8=0.5, D11=0)
#     writes: Tier-2 existence_layer_a.json, morphology.json, block_detect.json
#     writes: q1_d01..d07, q1_d09, q1_d10, q1_d12, q7_tier(provisional)
#     registers: candidate entries (catalog birth)
#
#   ▪ C01g boundary catalog
#     writes: Tier-2 boundary_left.json, boundary_right.json
#     writes: q3_* keys (~9 per boundary × 2 = ~18 per candidate)
#     does NOT depend on C01d because it reads its own five boundary sources
#     BUT waits until catalog birth so it can register per-candidate keys
#
#   ▪ Existence Layer B (SV concordance from C00 flashlight)
#     writes: Tier-2 existence_layer_b.json
#     writes: q7_layer_b_* keys
#
#   Phase 4a exit condition: catalog born + Tier-2 blocks for Q1/Q3/Q7-LayerB
#   written for every candidate. `q6_group_validation = NONE` for all.
#
#
# Phase 4b — GROUP PROPOSAL (C01i runs once per candidate, typically Tier ≤ 3)
#   ▪ C01i decomposition
#     reads: candidate list (tier ≤ 3) from registry
#     writes: Tier-2 internal_dynamics.json, band_composition.json
#     writes: q2_* keys, q6_n_HOM_REF, q6_n_HET, q6_n_HOM_INV, q6_n_Recombinant
#     registers: four groups per candidate (HOM_REF, HET, HOM_INV, RECOMBINANT)
#                plus GC/DCO subgroups where cheat24 resolves them
#     sets:     q6_group_validation = UNCERTAIN
#
#   Phase 4b exit condition: groups registered + q6_group_validation=UNCERTAIN
#   for every Tier≤3 candidate. Tier 4 candidates skipped (no decomposition).
#
#
# Phase 4c — GROUP VALIDATION (C01f runs once per Tier≤2 candidate)
#   ▪ C01f hypothesis tests — uses registered groups via comp_from_registry()
#     T1, T2: read `comp` built from registry (was k-means on PC1)
#     T3:     kin-pruned block comparison (group-independent)
#     T4:     inner/outer composition — reads registered groups
#     T5, T6: group-independent (PCA sub-region regime)
#     T7:     within-carrier substructure — reads HOM_INV group
#     T8:     Clair3 genotype concordance — reads all three classes
#     T9:     ancestry jackknife — reads Q-groups (unchanged)
#     T10:    theta het — reads all three classes
#     T11:    extended suppression — group-independent (Fst decay)
#
#     writes: Tier-2 hypothesis_verdict.json, existence_layer_c.json,
#             existence_layer_d.json
#     writes: q7_verdict, q7_verdict_confidence, q7_t1..t11 keys
#     writes: q3_extended_suppression, q6_family_linkage
#
#   ▪ Validation promotion rules (run at end of C01f per candidate):
#
#     IF t8_concordance ≥ 0.70 AND t9_jackknife_status = "robust":
#         q6_group_validation = SUPPORTED
#
#     IF Layer D Fisher OR p < 0.05 AND OR > 5:
#         q6_group_validation = VALIDATED   # strongest level
#
#     IF t9_jackknife_status = "fragile" (max_delta > 0.3):
#         q6_group_validation = SUSPECT     # demotion
#
#     Default if nothing triggers: stay at UNCERTAIN.
#
#   Phase 4c exit condition: q6_group_validation set for every Tier≤2 candidate.
#
#
# Phase 4d — GROUP-DEPENDENT EVIDENCE (parallel, gated by validation level)
#   ▪ cheat30 gds_by_genotype (Q5 age)
#       requires q6_group_validation ≥ NONE    (GDS works on all candidates)
#       but Fst sub-block requires ≥ SUPPORTED
#
#   ▪ cheat6 ancestry jackknife (if not already run in C01f)
#       requires q6_group_validation ≥ UNCERTAIN
#
#   ▪ burden_regression (Q6)
#       requires q6_group_validation = VALIDATED
#
#   ▪ classify_inversions first pass — can now assign labels for questions
#     whose requirements are met. Remaining questions stay unclassified.
#
#   Phase 4d exit condition: all cheats that CAN run given validation levels
#   have run. Tier-2 blocks + Tier-3 keys populated up to their gates.
#
#
# Phase 4e — FINAL CATALOG RE-SCORE (C01d pass-2)
#   ▪ C01d pass-2 re-reads boundary + hypothesis keys from registry
#     (not from --boundary_dir / --hyp_dir TSV files, though the TSV path
#     still works as a fallback)
#     fills: D8 = max(peel_score, hypothesis_verdict_score)
#            D11 = boundary_concordance from registered q3_*_verdict
#            D12 (already filled in pass-1, double-check)
#     recomputes: final_score = weighted sum of D1..D12
#                 dim_positive, tier
#     applies Cheat 25 viability test — DEAD candidates → Tier 4
#     rewrites: candidate_scores.tsv.gz (snapshot) + re-registers q1/q7 keys
#
#   ▪ characterize_candidate runs per candidate
#     walks the 18-block schema, assigns per-question status
#     writes characterization.json
#
#   ▪ classify_inversions final pass
#     assigns 14 axes + short_tag
#     writes classification.json
#
#   Phase 4e exit condition: every candidate has a final tier, verdict,
#   characterization.json, and classification.json. Done.


# ─────────────────────────────────────────────────────────────────────────────
# 6. The registry as source of truth
# ─────────────────────────────────────────────────────────────────────────────
#
# In v10, the "catalog" stops being a single TSV file and becomes a view over
# the registry. Three concrete changes:
#
#   (a) Scripts write Tier-2 JSON blocks via reg$evidence$write_block(). The
#       library extracts flat keys automatically per each block's schema and
#       appends them to per_candidate/<cid>/keys.tsv. Scripts do not write
#       keys directly.
#
#   (b) candidate_scores.tsv.gz is a SNAPSHOT projected from the registry
#       by C01d pass-2 for convenience + legacy downstream. It is no longer
#       the source of truth. The source of truth is the evidence_registry
#       per_candidate/<cid>/structured/ directory + keys.tsv.
#
#   (c) Queries like "all Tier-1 candidates with VALIDATED groups and Fst
#       sub-block populated" become reg$query$* methods, not grep on the
#       catalog file.
#
# The library provides fallback behavior: if any tier-2 block is missing
# its schema or fails validation, the script writes to the old per-module
# output directory and emits a warning. The pipeline runs without the
# library (via the same fallback), just without structured block writes.


# ─────────────────────────────────────────────────────────────────────────────
# 7. Why this lands on Option 3 in practice (but cheaper than it looked)
# ─────────────────────────────────────────────────────────────────────────────
#
# Of the three options originally considered:
#
#   Option 1 (re-run C01d three times): rejected — C01d is already written
#       to be idempotent on enrichment, but three passes hides what's really
#       happening (orthogonal evidence streams, not sequential scoring).
#
#   Option 2 (rewrite C01f internals): subsumed by Option 3. The only
#       rewrite needed in C01f is one helper that builds `comp` from
#       registered groups instead of from PC1 k-means.
#
#   Option 3 (split C01f conceptually into pre/post): adopted. But the
#       "split" is NOT two scripts — it's recognizing that different tests
#       in C01f have different group requirements. T3/T5/T6/T11 need none,
#       T4/T9 need the Q-ancestry groups (already handled), T1/T2/T7/T8/T10
#       need C01i's genotype groups. The orchestrator puts C01f AFTER C01i
#       and C01f reads groups from the registry — that's the whole "split".
#
# The conceptual change that matters: C01f is no longer "the hypothesis
# test script". It's "the group validation script". Its output promotes
# groups from UNCERTAIN to SUPPORTED/VALIDATED (or demotes to SUSPECT).
# The hypothesis verdict is a byproduct.


# ─────────────────────────────────────────────────────────────────────────────
# 8. What changes where (concrete patch map)
# ─────────────────────────────────────────────────────────────────────────────
#
# C01i STEP_C01i_decomposition_rewired_24_v934_registry.R
#   — add RECOMBINANT group registration (was missing)
#   — add RECOMBINANT_GC / RECOMBINANT_DCO subgroups from cheat24 output
#   — name convention: HOM_STD → HOM_REF (align with handoff)
#       Both names kept as aliases for a migration window.
#   — write q6_group_validation = "UNCERTAIN" at end
#
# C01f STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R
#   — add comp_from_registry() helper (new function, ~35 lines)
#   — patch 4 call sites that build `comp` via k-means (each becomes
#     `comp <- comp_from_registry(cid, sample_names, fallback_km = <old inline km>)`)
#   — add validation promotion block at end of per-candidate loop that
#     sets q6_group_validation based on t8/t9/Layer D
#   — optionally consume cheat24 RECOMBINANT group to exclude recombinant
#     samples from Fst group comparisons (cleaner signal)
#
# C01d STEP_C01d_candidate_scoring_wired_25_v934_registry.R
#   — no code change needed. Already supports --boundary_dir and --hyp_dir.
#   — Orchestrator calls it twice: pass-1 without those flags, pass-2 with them.
#   — In v10 with the library, --boundary_dir and --hyp_dir are replaced by
#     registry queries (single flag --from_registry).
#
# C01g STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R
#   — no ordering change. Runs in phase 4a in parallel with C01d pass-1.
#   — optional: migrate boundary_catalog_unified.tsv.gz writes to Tier-2
#     boundary_left.json + boundary_right.json via the library.
#
# Orchestrator run_phase4.sh (new)
#   — five SLURM array jobs with afterok dependencies:
#       4a_scoring_pass1   (array over chromosomes)
#       4a_boundary        (array over chromosomes, parallel with 4a_scoring_pass1)
#       4b_decomposition   (array over candidates, depends on 4a_scoring_pass1)
#       4c_hypothesis      (array over candidates, depends on 4b_decomposition)
#       4d_group_cheats    (array, depends on 4c_hypothesis)
#       4e_final_score     (single job, depends on 4c_hypothesis + 4d_group_cheats)
#
# registry_loader.R (new, Stage 2)
#   — three registry classes (sample / interval / evidence)
#   — write_block() with schema validation
#   — query methods (Layer 2 + Layer 3 per handoff section 7)
#
# schemas/structured_block_schemas/*.schema.json (new, Stage 2)
#   — 18 block schemas, each specifying required fields + keys_extracted
