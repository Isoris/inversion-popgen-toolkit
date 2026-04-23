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
# For (a) genome-wide Q, K should match the number of deep ancestry groups
# in the cohort. For 226 pure *C. gariepinus* catfish from a hatchery with
# ~8 broodlines, K ≈ 8 is biologically justified. [Quentin: verify the
# actual broodline count for the MS_Inversions cohort and adjust this
# number if it's different — the earlier draft of this doc said "~6 founder
# lines + two parent species", which was a cohort-conflation with the F1
# hybrid assembly paper's cohort.] At this K:
#   - K=5: collapses lines together, loses fine family structure
#   - K=8: ~28 samples per group on average — good for Fst stats
#   - K=12: ~19 per group — starts to get noisy for Fst
#   - K=20: ~11 per group — too small; Fst becomes noise-dominated;
#           jackknife "drop-one-group" barely moves the signal because
#           each group is tiny
#
# For (b) Cheat 5 — K determines how many Q components are scanned. Higher
# K gives more chances to find a high-Fst component. But at K=20 with
# 11 samples per group, the Fst estimator variance dominates. K=8 is the
# inflection point where signal starts beating noise. Do NOT go above K=12
# for cheat 5.
#
# For (c) Cheat 6 jackknife — same story. Removing one Q group at K=8
# removes ~28 samples, big enough to see real signal changes; at K=20
# removing ~11 samples is within noise. K=8 is appropriate.
#
# For (d) instant_q local Q — K determines the alphabet size of the
# assigned_pop label. nested_composition relies on label changes within
# a sample. At K=8 you have 8 possible labels; a two_block composite is
# easy to detect (2 out of 8 labels appearing, each ≥30%). At K=20 with
# noisy local Q estimates, spurious label changes make everything look
# fragmented → false-positive composite flags. K=8 is safer than K=20
# for composite detection.
#
# CONCLUSION: K=8 IS THE RIGHT DEFAULT for this 226-sample cohort. It's
# not arbitrary — it matches the biology (~8 deep groupings) AND the
# statistical floor (~28 samples per group for reliable Fst/jackknife).
#
# WHEN TO DEVIATE:
#
#   - Composite-flagged candidates at K=8 → re-run nested_composition at
#     K=12 ONLY for those candidates. If the composite pattern survives,
#     it's real. If it disappears at K=12, K=8 was too coarse.
#     Cheap: <20 candidates × 2 min each.
#
#   - A new cohort with different biology → re-pick K based on
#     expected_deep_groups AND samples/group ≥ 25.
#
#   - The full multi-K approach (evaluate at K=5, 8, 12, report stable
#     results) is principled but costs 3× more NGSadmix time + 3× the
#     downstream analysis. Skip for the manuscript; revisit post-submission.
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
