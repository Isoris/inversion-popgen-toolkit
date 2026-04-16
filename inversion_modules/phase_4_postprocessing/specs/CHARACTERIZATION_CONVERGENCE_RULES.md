# =============================================================================
# CHARACTERIZATION CONVERGENCE RULES
# =============================================================================
#
# For each of the 7 biological questions, this defines:
#   - Which keys are REQUIRED for an answer
#   - What "convergence" means (which keys must agree)
#   - What triggers each status level
#   - What CONTRADICTED looks like (specific failure patterns)
#
# Status levels:
#   ANSWERED      — multiple independent evidence lines converge on one label
#   PARTIAL       — some evidence but not enough for a confident label
#   MEASURED      — raw numbers exist but interpretation/cross-validation pending
#   CONTRADICTED  — evidence lines point in conflicting directions
#   EMPTY         — required keys not yet computed
#
# The label (e.g., "simple_polymorphic") is only assigned at ANSWERED.
# At PARTIAL, a provisional label + caveat is given.
# At MEASURED/EMPTY, no label — just "pending".
# At CONTRADICTED, no label — just "investigate" + the conflict description.
#
# =============================================================================


# ═══════════════════════════════════════════════════════════════════════════════
# Q1: WHAT IS IT?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: simple_polymorphic / rare_polymorphic / fixed_fossil /
#         nested_system / complex_system / ambiguous_signal / likely_artifact
#
# REQUIRED keys (minimum for any assessment):
#   q1_shape_class, q1_n_children, q1_landscape_category
#
# CONVERGENCE for ANSWERED:
#   Test 1: Shape is definitive
#     q1_shape_class = "strong_square" AND q1_squareness > 0.7
#     AND q1_n_children = 0
#     → simple_polymorphic (if genotype counts show 3 classes)
#     → rare_polymorphic (if one homozygote class absent)
#
#   Test 2: Nested detection agrees
#     q1_n_children >= 1 AND q1_landscape_category contains "nested"
#     AND at least 2 of: block_compactness > 0.5, nn_survives_40 = TRUE,
#     consensus_confidence in (HIGH, MEDIUM)
#     → nested_system or complex_system
#
#   Test 3: Artifact detection
#     q1_family_likeness_mean > 0.5 AND d8_peel < 0.3 (peel destroys signal)
#     AND q1_inv_likeness_mean < 0.4
#     → likely_artifact
#
# CONTRADICTED when:
#   - shape_class = "strong_square" BUT family_likeness > 0.5
#     (looks like inversion structurally but driven by family)
#   - n_children > 0 BUT landscape_category = "standard"
#     (detector found children but classifier says simple)
#   - nn_survives_80 = TRUE BUT peel_effect = "disappeared"
#     (robust to smoothing but not to kin removal — suspicious)
#
# PARTIAL when:
#   Shape and NN agree but fewer than 2 supporting lines
#
# MEASURED when:
#   Keys filled but cross-checks not computed (e.g., have shape but no peel)
#
# EMPTY when:
#   q1_shape_class is NULL (inv_detect hasn't run)


# ═══════════════════════════════════════════════════════════════════════════════
# Q2: WHAT'S HAPPENING INSIDE?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: clean / edge_conversion / internal_crossover / mosaic / not_assessed
#
# REQUIRED keys:
#   q2_n_recombinant, q2_class_stability_pct, q2_switching_kw_pval
#
# CONVERGENCE for ANSWERED:
#   Test 1: Clean interior
#     q2_n_recombinant = 0 AND q2_class_stability_pct > 98%
#     AND q2_class_entropy_mean < 0.2
#     → clean
#
#   Test 2: Edge conversion
#     q2_n_recombinant > 0 AND q2_recomb_left_n + q2_recomb_right_n > q2_recomb_center_n
#     AND q2_switching_kw_pval < 0.05 (classes have different switch rates)
#     → edge_conversion
#
#   Test 3: Internal crossover
#     q2_n_double_crossover > 0 AND q2_recomb_center_n > 0
#     AND these center events are > 500kb from both breakpoints
#     → internal_crossover (requires edge_conversion also usually present)
#
#   Test 4: Mosaic
#     q2_mean_switch_rate_std > 0.05 OR q2_n_windows_unimodal > 20% of total
#     AND q2_class_stability_pct < 80%
#     → mosaic (complex, possibly overlapping inversions)
#
# CONVERGENCE BOOST (strengthens answer):
#   q2_diversity_gradient = "deep_U" consistent with edge_conversion or older
#   q2_phase_block_ratio > 2 confirms recombination suppression inside
#   q2_indel_class_concordance > 0.7 confirms PCA classes are real
#
# CONTRADICTED when:
#   - q2_n_recombinant = 0 BUT q2_class_stability_pct < 90%
#     (no recombinants detected but samples are unstable — decomposition problem?)
#   - q2_phase_block_ratio < 1 (longer blocks OUTSIDE than inside —
#     opposite of recombination suppression expectation)
#   - q2_indel_class_concordance < 0.3 (INDEL classes don't match PCA classes —
#     the PCA grouping might not reflect real haplotype structure)


# ═══════════════════════════════════════════════════════════════════════════════
# Q3: WHAT ARE THE BOUNDARIES DOING?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: hard_hard / hard_soft / hard_eroded / soft_soft / eroded_eroded /
#         fossil_fossil / not_assessed
#
# REQUIRED keys:
#   q3_left_hardness, q3_right_hardness, q3_left_verdict, q3_right_verdict
#
# CONVERGENCE for ANSWERED:
#   For EACH boundary independently:
#     "hard": hardness > 0.7 AND verdict = "confirmed_structural"
#             AND n_cheats >= 3 AND at least one of:
#             clip_count > 5, or sv_bp within 10kb, or fst_step > 0.05
#
#     "soft": hardness 0.3-0.7 AND verdict in ("likely_structural", "probable")
#             AND n_cheats >= 2
#
#     "eroded": hardness < 0.3 OR verdict = "unresolved"
#               BUT some evidence exists (clip_count > 0 or depth_anomaly present)
#
#     "fossil": is_fossil = TRUE (Cheat 17 detected orphaned boundary scar)
#
#   Label = left_right (e.g., "hard_soft")
#
# CONVERGENCE BOOST:
#   q3_left_concordance_kb < 20 (all boundary estimates agree within 20kb)
#   q3_left_sv_bp within 5kb of q3_left_bp (SV and PCA agree)
#   Carrier reconciliation: q3_carrier_concordance > 0.8
#
# CONTRADICTED when:
#   - q3_left_hardness > 0.7 BUT q3_left_n_cheats = 0
#     (sim_mat says sharp but no independent evidence)
#   - q3_left_sv_bp exists BUT > 100kb from q3_left_bp
#     (SV and PCA boundaries wildly disagree — wrong match?)
#   - q3_left_fst_step < 0 (Fst DECREASES at boundary — opposite of expected)
#   - q3_dropout_rate > 0.5 (more than half of PCA carriers missing from SV —
#     either massive dropout or PCA is wrong)


# ═══════════════════════════════════════════════════════════════════════════════
# Q4: HOW DID IT FORM?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: NAHR / NHEJ / MMEJ / TE_mediated / MMBIR / unclassified
#
# REQUIRED keys:
#   At least 2 of: q4_has_inverted_sd, q4_junction_type_left, q4_te_enrichment
#   (mechanism requires multiple evidence lines)
#
# CONVERGENCE for ANSWERED (follows Porubsky et al. 2022 decision tree):
#
#   NAHR (3 tests must agree):
#     q4_has_inverted_sd = TRUE (inverted SD pair flanking breakpoints)
#     AND q4_sd_identity_pct > 90% (high identity = recent enough for HR)
#     AND q4_biser2_concordance = "agree_nahr" (independent SD finder agrees)
#     OPTIONAL BOOST: q4_mh_length >= 8 (long microhomology = HR signature)
#     → mechanism_confidence = "confirmed" if all 3, "likely" if 2/3
#
#   NHEJ (2 tests):
#     q4_has_inverted_sd = FALSE (no SD substrate)
#     AND q4_mh_length_left <= 2 AND q4_mh_length_right <= 2 (blunt or tiny MH)
#     → mechanism_confidence = "likely" (NHEJ is a diagnosis of exclusion)
#
#   MMEJ (2 tests):
#     q4_mh_length_left 3-20 bp AND q4_mh_length_right 3-20 bp
#     AND q4_has_inverted_sd = FALSE
#     → mechanism_confidence = "probable"
#
#   TE_mediated (3 tests):
#     q4_te_enrichment = "ENRICHED"
#     AND q4_te_same_family = TRUE (same TE family at both breakpoints)
#     AND q4_te_same_orientation = "inverted"
#     → mechanism_confidence = "likely"
#
#   MMBIR (diagnosis of exclusion):
#     q4_has_inverted_sd = FALSE AND q4_mh_length < 3
#     AND junction shows templated insertions or complex rearrangement
#     → mechanism_confidence = "tentative"
#
# CONTRADICTED when:
#   - q4_has_inverted_sd = TRUE BUT q4_biser2_concordance = "disagree"
#     (minimap2 finds SDs but BISER2 doesn't — check alignment quality)
#   - q4_junction_type_left = "te_mediated" BUT q4_te_enrichment = "DEPLETED"
#     (junction says TE but breakpoint region is TE-depleted — paradox)
#   - q4_mh_length_left > 20 (very long MH — could be SD that minimap2 missed)
#
# PARTIAL when:
#   Only junction OR SD data available, not both


# ═══════════════════════════════════════════════════════════════════════════════
# Q5: HOW OLD IS IT?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: old / intermediate / young / unranked
#
# REQUIRED keys:
#   q5_gds_gap AND q5_fst_b1b3 (minimum: two independent age proxies)
#
# CONVERGENCE for ANSWERED:
#   Test 1: Ranking agreement
#     q5_gds_gap_percentile in same third as Fst rank
#     (both say top third → old, both say bottom third → young)
#     AND q5_gds_fst_spearman_p < 0.05 (correlation is significant across all candidates)
#
#   Test 2: Diversity profile consistent
#     old: q5_diversity_shape = "deep_U" or "shallow_U" (nucleotide diversity
#          elevated near breakpoints from accumulated gene flux)
#     young: q5_diversity_shape = "flat" (no time for gene flux)
#
#   Test 3: Per-arrangement divergence consistent
#     old: q5_dxy_ratio >> 1 (arrangements diverged inside inversion)
#     AND q5_fixed_differences > 0 (some sites fixed for different alleles)
#     young: q5_dxy_ratio ≈ 1, q5_fixed_differences = 0
#
# CONVERGENCE BOOST:
#   q5_dollo_mya provides absolute calibration (if synteny tested)
#   q5_boundary_status contains "eroded" → consistent with old
#
# CONTRADICTED when:
#   - GDS says old (top third) BUT Fst says young (bottom third)
#     → one proxy is biased; check carrier set and sample size per arrangement
#   - q5_diversity_shape = "deep_U" BUT q5_gds_gap_percentile < 33
#     → U-shape without arrangement divergence — possible complex demography
#   - q5_fixed_differences > 10 BUT q5_gds_gap_percentile < 50
#     → fixed differences suggest old, but GDS says not — possible gene conversion
#       resetting GDS while fixed diffs accumulate at non-converted sites
#
# MEASURED when:
#   GDS and Fst both computed but ranking/comparison not done
#   (need classify_inversions.R to run the cross-candidate comparison)


# ═══════════════════════════════════════════════════════════════════════════════
# Q6: HOW COMMON IS IT?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: rare / low / intermediate / high / nearly_fixed
#
# REQUIRED keys:
#   q6_freq_inv, q6_n_HOM_STD, q6_n_HET, q6_n_HOM_INV
#
# CONVERGENCE for ANSWERED:
#   Test 1: Frequency classification (straightforward)
#     freq_inv < 0.05 → rare
#     freq_inv 0.05-0.15 → low
#     freq_inv 0.15-0.50 → intermediate
#     freq_inv 0.50-0.85 → high
#     freq_inv > 0.85 → nearly_fixed (but see FM5: might be species-diagnostic)
#
#   Test 2: HWE consistency
#     q6_hwe_deviation tells you whether genotype ratios are expected
#     CAVEAT: F1 hatchery violates HWE by construction, so hwe_consistent
#     is the expected result for most inversions, not evidence of selection
#
#   Test 3: Family distribution
#     q6_freq_cv_across_qgroups low (< 0.3) → evenly distributed = true polymorphism
#     q6_freq_cv_across_qgroups high (> 0.5) → concentrated in few families
#     q6_family_linkage = "multi_family" → not a single-family artifact
#
# CONVERGENCE is relatively easy for Q6 — frequency is a direct measurement.
# The complications are in INTERPRETATION (is the frequency driven by selection
# or drift? — that's Q6 selection_pattern, which requires Q5 age context).
#
# CONTRADICTED when:
#   - q6_freq_inv differs significantly between PCA-derived and SV-derived counts
#     (q3_carrier_concordance < 0.6 → the two methods disagree on WHO carries it)
#   - q6_n_total < 200 (too many unclassified samples — decomposition incomplete)
#   - q6_jackknife_max_delta > 0.10 (removing one Q-group changes frequency by >10%
#     → frequency estimate is fragile, driven by one family)


# ═══════════════════════════════════════════════════════════════════════════════
# Q7: IS IT REAL?
# ═══════════════════════════════════════════════════════════════════════════════
#
# Labels: confirmed / likely / candidate / artifact
#
# REQUIRED keys:
#   At least q7_layer_a_detected (PCA) — the minimum detection
#
# CONVERGENCE for ANSWERED:
#   "confirmed" requires ≥3 of 4 layers passing:
#     Layer A: q7_layer_a_detected = TRUE AND q7_layer_a_beta_qval < 0.10
#     Layer B: q7_layer_b_detected = TRUE AND (delly OR manta)
#     Layer C: q7_layer_c_ghsl_detected = TRUE AND q7_layer_c_ghsl_quality in (HIGH, MODERATE)
#     Layer D: q7_layer_d_fisher_p < 0.05 AND q7_layer_d_fisher_or > 3
#
#     PLUS: q6_family_linkage = "multi_family" OR q7_t9_jackknife_status = "robust"
#     (at least one of these — rules out single-family artifact)
#
#   "likely" requires 2/4 layers passing + no contradictions
#
#   "candidate" requires 1/4 layers passing
#
#   "artifact" requires:
#     q7_verdict contains "artifact" or "family_structure" or "H1_family"
#     OR (q7_t9_jackknife_status = "single_family_fragile" AND
#         q1_family_likeness_mean > 0.5 AND d8_peel < 0.3)
#
# CONTRADICTED when:
#   - Layer A passes but Layer D fails significantly
#     (PCA says inversion, but breakpoint reads DON'T correlate with PCA carriers)
#     → either PCA is detecting family structure, or SV breakpoints are wrong
#   - Layer B passes but Layer A fails
#     (SV callers see breakpoints but PCA sees no block — fixed inversion?
#      check if all samples are carriers → SV-only tier, not contradiction)
#   - Layer C (GHSL) passes but Layer A fails
#     (haplotype divergence without PCA block — small inversion below PCA resolution?
#      or GHSL detecting something PCA can't see at 50kb windows)
#   - q7b_observed_dropout_pct > 50% AND q7_layer_d_fisher_p > 0.05
#     (massive carrier dropout may have destroyed the Fisher test —
#      don't conclude "no association", conclude "test compromised")


# ═══════════════════════════════════════════════════════════════════════════════
# IMPLEMENTATION NOTES
# ═══════════════════════════════════════════════════════════════════════════════
#
# 1. The characterization engine reads the key registry and applies these rules.
#    It does NOT compute anything new — it only interprets existing values.
#
# 2. CONTRADICTED is the most valuable status. It flags problems that need
#    human investigation. The specific contradiction pattern tells you where
#    to look.
#
# 3. Counter-proof: when evidence AGAINST a conclusion exists, it's recorded
#    alongside the evidence FOR. Example:
#      Q5 age = "old"
#      Evidence FOR: GDS gap 90th percentile, Fst top third, diversity = U-shape
#      Evidence AGAINST: zero fixed differences between arrangements
#      Status: PARTIAL (not ANSWERED, because the counter-evidence matters)
#
# 4. The completion % (System 2) and characterization status (System 3) are
#    displayed together but computed independently:
#      Collection: 78% (274/352 keys)
#      Characterization: 5/7 ANSWERED, 1 PARTIAL, 1 MEASURED
#
#    A candidate can have high collection but low characterization (lots of
#    measurements but they don't converge) or low collection but high
#    characterization (few measurements but they all agree perfectly).
#
# 5. The convergence rules reference specific key names from the registry
#    specification v2 (352 keys). When a key doesn't exist yet, the rule
#    that depends on it simply can't fire, which reduces the maximum
#    achievable status for that question.
#
# =============================================================================
