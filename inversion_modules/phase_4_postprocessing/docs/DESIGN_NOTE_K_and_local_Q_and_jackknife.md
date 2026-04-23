# =============================================================================
# DESIGN NOTE — K choice and local_Q cache
# =============================================================================
# Date: 2026-04-16
# Scope: answers to two real operational questions about the ancestry layer
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 1. local_Q: where it lives, who produces it, who consumes it
# ─────────────────────────────────────────────────────────────────────────────
#
# THERE ARE TWO ARTIFACTS FROM THE ANCESTRY LAYER. They are not the same.
#
# ARTIFACT A — per-sample per-window ancestry track (raw)
#   Path:       ${PRECOMP_OUT}/local_Q/<chr>.local_Q_samples.tsv.gz
#   Producer:   utils/ancestry_bridge.R --prepare  (calls instant_q.cpp)
#   Format:     one row per (window × sample) with columns:
#                 chrom, start_bp, end_bp, sample_id,
#                 Q1..QK, max_q, delta12, entropy, ena, assigned_pop
#   Consumers:  nested_composition.py / STEP_C01i_c_nested_composition.py
#               (walks assigned_pop sequences per sample to detect
#                two_block / multi_block / continuous_gradient patterns)
#
# ARTIFACT B — per-window summary merged into precomp RDS
#   Path:       ${PRECOMP_DIR}/<chr>.precomp.rds → list(dt = ...)
#   Producer:   C01a via ancestry_bridge::merge_local_Q_into_invlikeness()
#   Format:     additional columns on the inv_likeness table:
#                 localQ_delta12, localQ_entropy, localQ_ena
#               These are AVERAGES over samples for each window.
#   Consumers:  downstream scoring that wants per-window ancestry diagnostics
#               (NOT nested_composition — the averages discard the per-sample
#                signal that nested_composition needs).
#
# KEY POINT: nested_composition REQUIRES artifact A. It does not work from
# artifact B. If your precomp was run WITHOUT `ancestry_bridge.R --prepare`
# first, only artifact B exists (and maybe not even that — merging is
# conditional on dir.exists(local_q_dir)).
#
# TO CHECK what you have right now on LANTA:
#   ls ${PRECOMP_OUT}/local_Q/*.local_Q_samples.tsv.gz    # artifact A
#   # vs
#   Rscript -e 'pc <- readRDS("<chr>.precomp.rds"); print(grep("^localQ", names(pc$dt), value=TRUE))'
#   # artifact B: should print "localQ_delta12" "localQ_entropy" "localQ_ena"
#   # if these are all you have, nested_composition will write stubs.
#
# TO PRODUCE ARTIFACT A:
#   1. Make sure Engine B is compiled (instant_q binary present)
#   2. Run: Rscript utils/ancestry_bridge.R --prepare --chr <CHR> --K 8
#   3. This writes ${PRECOMP_OUT}/local_Q/<CHR>.local_Q_samples.tsv.gz
#   4. Optionally then re-run C01a to get artifact B too
#
# (The precomp C01a also includes per-window localQ_* columns from
# ancestry_bridge::merge_local_Q_into_invlikeness(). That merge reads the
# already-cached artifact A and averages it. So artifact B ⊂ artifact A:
# if you have A you can derive B. Not the other way around.)
#
# v10.1 BEHAVIOR if artifact A is absent:
#   STEP_C01i_c_nested_composition.py writes a stub block with
#   composite_flag = "unknown_no_engine_b". Downstream treats this as
#   "clean for gating purposes but record the missing signal".
#   NO composite detection happens without artifact A.
#
# RECOMMENDATION: run ancestry_bridge.R --prepare as a blocking
# precompute BEFORE C01a, not after. Then C01a merges artifact A into its
# RDS, AND nested_composition can read A directly. The code already
# supports this ordering — it's just that the Sep 2026 runs were done
# before ancestry_bridge existed.


# ─────────────────────────────────────────────────────────────────────────────
# 2. Why K=8 for this cohort (and when to deviate)
# ─────────────────────────────────────────────────────────────────────────────
#
# K appears in FOUR places in your current pipeline:
#
#   (a) NGSadmix Q matrix for the whole cohort (once, stored in .qopt file)
#       — this is genome-wide family/ancestry structure from thinned SNPs.
#   (b) Cheat 5 family-Fst scan in C01a (loops over each K, line 1142)
#       — finds which Q component best separates high/low allele frequency.
#   (c) Cheat 6 jackknife in the ancestry groups (hardcoded K=8 at line 120)
#       — leave-one-ancestry-group-out Fst delta.
#   (d) instant_q local Q (per window, same K as genome-wide)
#       — assigned_pop per window per sample, used by nested_composition.
#
# The question is: what's the RIGHT K for each use?
#
# UPFRONT — this cohort is NOT structure-matched by K. The hatchery has
# 20+ small families, so the K that would actually match the ancestry
# structure is ~20. But K is used here as a *scan parameter for inversion
# regime detection*, not as a tool for describing cohort structure, and
# the right K for inversion detection is well below the structure K. The
# ancestry layer routinely scans K=2..12 (and occasionally up to 20) to
# verify stability; K=8 is the default presented here for per-candidate
# downstream steps.
#
# Why K=8 is the default for inversion-regime scans:
#
#   1. An inversion usually has a dosage signal that concentrates in
#      1-2 Q components at moderate K. At K=8 you have enough resolution
#      that an inversion's carriers cluster into 1-2 components cleanly
#      (the rest of the structure lives in the other 6-7). Going higher
#      fragments the inversion signal across too many components.
#
#   2. Statistical tractability. At K=8 the average group size is ~28
#      samples, which is the minimum useful regime for Fst / jackknife
#      tests on a 226-sample cohort. At K=20 (structure-matched), you
#      get ~11 samples per group — Fst estimator variance dominates and
#      "drop-one-group" jackknife barely moves the signal because each
#      group is too small to matter.
#
#   3. Composite-detection alphabet. nested_composition reads the local
#      assigned_pop label per window per sample; at K=8 you have 8
#      possible labels and a two_block composite (2 labels, each ≥30%)
#      is easy to distinguish from noise. At K=20 with noisy local Q,
#      spurious label changes make everything look fragmented →
#      false-positive composites.
#
# So K=8 is a tool choice tuned for the inversion-regime job, NOT a
# claim that the cohort has ~8 deep biological groupings. Per-use
# notes:
#
#   (a) genome-wide Q: K=8 is the default downstream consumers assume,
#       but the ancestry layer produces Q at K=2..12 (K sweep). Any
#       downstream step that wants a different K can pull it from the
#       K-sharded cache.
#   (b) Cheat 5 family-Fst scan: K=8 is the sweet spot between Q
#       component count and per-component sample size. Do not go above
#       K=12 — estimator variance takes over.
#   (c) Cheat 6 jackknife: uses K=8 because removing one of 8 groups
#       removes ~28 samples, big enough to see real signal changes.
#   (d) instant_q local Q: K=8 gives the 8-label alphabet that
#       composite-detection is calibrated for.
#
# CONCLUSION: K=8 IS THE RIGHT DEFAULT FOR INVERSION-REGIME DETECTION in
# this 226-sample cohort. It's a scan-parameter choice that balances
# (signal concentration in 1-2 components) × (per-component sample
# size ≥25) × (composite-detection alphabet sensible for local Q).
# It is NOT chosen to match the hatchery's ancestry structure —
# structure-matching would argue for K=20 or higher and would be the
# wrong choice for this pipeline's statistics.
#
# WHEN TO DEVIATE:
#
#   - Composite-flagged candidates at K=8 → re-run nested_composition at
#     K=12 ONLY for those candidates. If the composite pattern survives,
#     it's real. If it disappears at K=12, K=8 was too coarse.
#     Cheap: <20 candidates × 2 min each.
#
#   - Ancestry-structure questions (how many founder lines, which
#     samples share broodline background) → use the full K=2..12 sweep
#     directly, don't rely on the K=8 default.
#
#   - A new cohort with different sample size → re-pick K to keep
#     samples/component ≥ 25; that's the binding statistical constraint.
#
# IMPLEMENTATION: make K a config variable, not a hardcoded 8.
#   Add to 00_inversion_config.sh:
#       export ANCESTRY_K=${ANCESTRY_K:-8}
#   In cheat6_ancestry_jackknife_v934.R line 120:
#       for (k in 1:8) { ... }
#     →
#       K_ancestry <- as.integer(Sys.getenv("ANCESTRY_K", "8"))
#       for (k in seq_len(K_ancestry)) { ... }
#
# Everything else reads K from the Q matrix's ncol, so it's already
# adaptive. Only cheat6 has the hardcoded 8.


# ─────────────────────────────────────────────────────────────────────────────
# 3. Jackknife semantics — four distinct outcomes, not one axis
# ─────────────────────────────────────────────────────────────────────────────
#
# The old v10 rule "if t9_jackknife = fragile → SUSPECT" was wrong because
# it conflated two biological categories:
#
#   (A) SUSPECT — the PCA groups are confounded with family structure.
#       Removing one ancestry group eliminates the signal because the
#       signal WAS just that group vs. others. The "inversion" is not
#       a real arrangement — it's a family marker masquerading as one.
#       → demote to SUSPECT.
#
#   (B) SINGLE-FAMILY POLYMORPHISM — the inversion IS real, but it's
#       segregating in only one family line. Removing that family makes
#       the signal vanish because without the carriers, there's nothing
#       to see. This is not a problem; family-restricted polymorphisms
#       are common and scientifically important.
#       → do NOT demote. Tag family_linkage = single_family,
#         polymorphism_class = family_restricted.
#
# Cheat6 already distinguishes these:
#   jackknife_verdict = "single_family_fragile"     → case B
#   jackknife_verdict = "pca_family_confounded"     → case A (not currently
#                                                      emitted by cheat6 but
#                                                      can be derived)
#
# The v10.1.1 patch (C01f_jackknife_semantics_patch.R) propagates the
# distinction through compute_group_validation() with three consequences:
#
#   1. "single_family_fragile" → tag family_specific_polymorphism,
#      cap at SUPPORTED (don't promote to VALIDATED — we don't have
#      Fisher OR confirmation), but do NOT demote to SUSPECT.
#
#   2. "pca_family_confounded" → SUSPECT (same as old rule).
#
#   3. High T8 Clair3 concordance on a single_family_fragile candidate
#      is now sufficient for SUPPORTED promotion. Rationale: if Clair3
#      discrete genotypes agree with the PCA groups within the single
#      family, the groups are right. The inversion is real and
#      genotyped correctly. We just can't promote to VALIDATED because
#      Fisher OR test requires cohort-wide carrier/non-carrier comparison,
#      which doesn't apply to a family-restricted variant.
#
# This actually matters for the manuscript. Without this patch, any
# family-restricted inversion would have ended up as SUSPECT and its
# Q5/Q6 would be unclassified. With this patch, family-restricted
# inversions are properly cataloged with polymorphism_class =
# family_restricted and get their age/origin axes filled.
#
# This is also the answer to "can we find inversions that changed the
# karyotype number" — those are often species-restricted or line-restricted,
# which previously looked SUSPECT under the old rule but are actually the
# most interesting biological cases.


# ─────────────────────────────────────────────────────────────────────────────
# 4. What's still hardcoded that should be configurable
# ─────────────────────────────────────────────────────────────────────────────
#
# Current hardcoded values that should move to 00_inversion_config.sh:
#
#   ANCESTRY_K                 default 8
#   DELTA_THRESH               default 0.10 (cheat6, defines "robust")
#   COMPOSITE_TWO_BLOCK_FRAC   default 0.20 (nested_composition flag threshold)
#   COMPOSITE_MULTI_BLOCK_FRAC default 0.20
#   COMPOSITE_HOMOG_FRAC       default 0.80
#   MOSAIC_SHORT_BP            default 100000 (cheat24)
#   MOSAIC_LONG_BP             default 500000
#   BREAKPOINT_WINDOW_BP       default 50000
#
# None of these need changing for the first manuscript run — the defaults
# are sensible. But making them env-var overridable means you can re-run
# one candidate with adjusted thresholds without editing code.
