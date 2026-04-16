# =============================================================================
# INVERSION EVIDENCE REGISTRY v2 — Updated Specification
# =============================================================================
#
# v2 changes (2026-04-16 audit session):
#
#   1. Q7 Layer C → GHSL haplotype contrast (independent from PCA)
#      Q7 Layer D → Genotype-breakpoint association (was Layer C)
#      Total: 4 independence layers (A/B/C/D), not 3
#
#   2. NEW: Q7B — SV caller audit keys (failure mode tracking)
#      Tracks stochastic dropout, BND fragmentation, filter chain losses
#      27 new keys documenting WHERE evidence was lost
#
#   3. NEW: Q3 carrier reconciliation keys
#      PCA vs SV carrier agreement, dropout flags, population prior status
#      8 new keys
#
#   4. Updated Q7 confidence tiers: pathway-based, not count-based
#      Tier 1 = 4/4 or 3/3 layers converge
#      Tier 2 = PCA+GHSL or PCA+SV without association
#      SV-only tiers for orphan inversions
#
#   5. Completion calculation: keys_resolved / (total - not_applicable) × 100
#
# Total keys: 317 (original) + 35 (new) = 352
#
# =============================================================================


# ═══════════════════════════════════════════════════════════════════════════════
# Q1–Q6: UNCHANGED from v1 (see INVERSION_REGISTRY_SPECIFICATION.md)
# ═══════════════════════════════════════════════════════════════════════════════
#
#   Q1: WHAT IS IT?                    49 keys (unchanged)
#   Q2: WHAT'S HAPPENING INSIDE IT?    40 keys (unchanged)
#   Q3: WHAT ARE THE BOUNDARIES DOING? 73 keys + 8 new = 81 keys
#   Q4: HOW DID IT FORM?              47 keys (unchanged)
#   Q5: HOW OLD IS IT?                39 keys (unchanged)
#   Q6: HOW COMMON IS IT?             28 keys (unchanged)
#   Q7: IS IT REAL?                   41 keys + 27 new = 68 keys
#
#   TOTAL: 317 + 35 = 352 keys


# ═══════════════════════════════════════════════════════════════════════════════
# Q3 NEW KEYS: CARRIER RECONCILIATION (PCA vs SV)
# ═══════════════════════════════════════════════════════════════════════════════
#
# These keys track the agreement between PCA-derived carrier assignments
# and SV-derived genotypes, and document where stochastic dropout caused
# discrepancies. The PCA carrier set is treated as truth for matched
# candidates (hundreds of SNPs > a few breakpoint reads).
#
# Added to Q3 because they relate to boundary/breakpoint evidence quality.

#   q3_n_carriers_pca         int       C01i    FALSE   PCA-derived carrier count (HET + HOM_INV)
#   q3_n_carriers_sv          int       C00     TRUE    SV VCF GT carrier count
#   q3_carrier_concordance    numeric   audit   FALSE   fraction of PCA carriers also SV carriers
#   q3_n_pca_carrier_sv_ref   int       audit   FALSE   samples PCA=carrier but SV=0/0 (dropout)
#   q3_n_sv_carrier_pca_ref   int       audit   FALSE   samples SV=carrier but PCA=REF (rare)
#   q3_dropout_rate           numeric   audit   FALSE   q3_n_pca_carrier_sv_ref / q3_n_carriers_pca
#   q3_population_prior_applied  bool   regeno  TRUE    was population regenotyping run?
#   q3_n_rescued_by_prior     int       regeno  TRUE    samples rescued from 0/0 to carrier by prior


# ═══════════════════════════════════════════════════════════════════════════════
# Q7 UPDATED: IS IT REAL? (68 keys, was 41)
# ═══════════════════════════════════════════════════════════════════════════════
#
# FOUR INDEPENDENCE LAYERS (was three):
#
# ┌──────────┐     ┌──────────┐     ┌──────────┐
# │ Layer A  │     │ Layer B  │     │ Layer C  │
# │ Local PCA│     │SV callers│     │  GHSL    │
# │(dosage   │     │(BAM read │     │(phased   │
# │covariance│     │ geometry)│     │ genotype │
# │lostruct) │     │DELLY+    │     │ haplotype│
# │          │     │Manta     │     │ contrast)│
# └────┬─────┘     └────┬─────┘     └────┬─────┘
#      │                 │                │
#      └────────┬────────┴────────┬───────┘
#               ▼                 ▼
#        ┌──────────────────────────────┐
#        │         Layer D              │
#        │  Genotype-breakpoint         │
#        │  association (Fisher OR)     │
#        │  Links A+C carriers to       │
#        │  B breakpoint evidence       │
#        └──────────────────────────────┘
#
# WHY 4 LAYERS:
#   A uses ANGSD dosage → lostruct PCA → sim_mat blocks
#   B uses BAM reads → discordant pairs + split reads → SV calls
#   C uses Clair3 phased VCFs → within-sample haplotype divergence
#   D links them: do the same fish show up as carriers in A, B, and C?
#
#   A and C are BOTH "population structure" but from DIFFERENT data:
#     A = ANGSD dosage (genotype likelihoods, unphased, callable regions only)
#     C = Clair3/WhatsHap phased genotypes (phased, includes some repeats)
#   A false positive from family structure would show HIGH between-band
#   concordance in GHSL (same family = similar haplotypes). A real inversion
#   shows LOW between-band concordance (different arrangements = different haplotypes).
#
# Summary label (stored as "confidence"):
#   confirmed       — ≥3/4 layers + multi-family + Fisher OR significant
#   likely          — 2/4 layers positive (A+B, A+C, or B+C)
#   candidate       — 1/4 layers (PCA only, SV only, or GHSL only)
#   artifact        — hypothesis tests say family structure
#

# Detail keys — LAYER A (local PCA): UNCHANGED from v1
#   q7_layer_a_detected       bool      C01b    FALSE
#   q7_layer_a_inv_likeness   numeric   C01a    FALSE
#   q7_layer_a_beta_pval      numeric   C01a    FALSE
#   q7_layer_a_beta_qval      numeric   fdr     FALSE
#   q7_layer_a_core_family    cat       C01b    FALSE
#   q7_layer_a_pa_pattern     cat       C01b    FALSE

# Detail keys — LAYER B (SV callers): UNCHANGED from v1
#   q7_layer_b_detected       bool      C00     TRUE
#   q7_layer_b_delly          bool      C00     TRUE
#   q7_layer_b_manta          bool      C00     TRUE
#   q7_layer_b_bnd_triang     bool      C00     TRUE
#   q7_layer_b_n_carriers     int       C00     TRUE
#   q7_layer_b_pe_support     int       C00     TRUE
#   q7_layer_b_sr_support     int       C00     TRUE
#   q7_layer_b_cipos          text      C00     TRUE
#   q7_layer_b_ciend          text      C00     TRUE

# Detail keys — LAYER C (GHSL haplotype contrast): NEW
#   q7_layer_c_ghsl_detected  bool      C04     FALSE   GHSL partition found at this region?
#   q7_layer_c_ghsl_contrast  numeric   C04     FALSE   mean GHSL contrast score across core windows
#   q7_layer_c_ghsl_n_pass    int       C04     FALSE   number of windows with GHSL PASS
#   q7_layer_c_ghsl_pct_pass  numeric   C04     FALSE   fraction of core windows with GHSL PASS
#   q7_layer_c_ghsl_quality   cat       C04     FALSE   HIGH (>70% PASS) / MODERATE / LOW / ABSENT
#   q7_layer_c_partition_stable bool    C04     FALSE   partition stability > 0.5?
#   q7_layer_c_ghsl_version   cat       C04     TRUE    v4 / v5 / v6

# Detail keys — LAYER D (genotype-breakpoint association): was Layer C in v1
#   q7_layer_d_tested         bool      C01f    FALSE   Fisher/Armitage test run?
#   q7_layer_d_fisher_or      numeric   C01f    FALSE   Fisher's exact test odds ratio
#   q7_layer_d_fisher_p       numeric   C01f    FALSE   Fisher's exact test p-value
#   q7_layer_d_armitage_z     numeric   C01f    FALSE   Cochran-Armitage trend test Z
#   q7_layer_d_armitage_p     numeric   C01f    FALSE   Cochran-Armitage p-value
#   q7_layer_d_concordance    cat       C01f    FALSE   descriptive concordance statement

# Independence summary: UPDATED for 4 layers
#   q7_n_layers_tested        int       classify mixed  how many layers were testable (0-4)
#   q7_n_layers_passed        int       classify mixed  how many layers passed (0-4)
#   q7_independence_class     cat       classify mixed  4/4=confirmed, 3/4=confirmed,
#                                                        2/4=likely, 1/4=candidate, 0/4=artifact


# ═══════════════════════════════════════════════════════════════════════════════
# Q7B NEW: SV CALLER AUDIT KEYS (failure mode tracking)
# ═══════════════════════════════════════════════════════════════════════════════
#
# These keys document WHERE and HOW evidence was lost in the SV calling
# pipeline. They are diagnostic — they don't contribute to the confidence
# tier directly, but they explain WHY a candidate might have low D7 or
# missing Layer B, and they guide the population regenotyping rescue.
#
# FAILURE MODE 1: Per-sample stochastic dropout
#   q7b_delly_raw_carriers    int       raw_vcf TRUE    carriers in raw VCF (before strict filter)
#   q7b_delly_strict_carriers int       strict  TRUE    carriers in strict catalog (after QUAL>=300)
#   q7b_delly_carrier_loss    int       audit   TRUE    raw - strict (lost to filtering)
#   q7b_manta_raw_carriers    int       raw_vcf TRUE    carriers in Manta raw VCF
#   q7b_manta_pass_carriers   int       pass    TRUE    carriers in Manta PASS catalog
#   q7b_manta_carrier_loss    int       audit   TRUE    raw - PASS
#   q7b_expected_dropout_pct  numeric   binom   TRUE    expected dropout % at 9× (binomial model)
#   q7b_observed_dropout_pct  numeric   audit   FALSE   actual PCA-carrier vs SV-carrier discordance
#
# FAILURE MODE 2: INV vs BND fragmentation
#   q7b_delly_inv_present     bool      4D      TRUE    DELLY INV call exists at this locus?
#   q7b_delly_bnd_3to3        bool      4E      TRUE    DELLY BND with CT=3to3 at left boundary?
#   q7b_delly_bnd_5to5        bool      4E      TRUE    DELLY BND with CT=5to5 at right boundary?
#   q7b_manta_inv_present     bool      4G      TRUE    Manta INV call exists?
#   q7b_manta_bnd_inv3        bool      4G_raw  TRUE    Manta raw BND with INV3 signal?
#   q7b_manta_bnd_inv5        bool      4G_raw  TRUE    Manta raw BND with INV5 signal?
#   q7b_bnd_rescued           bool      STEP06  TRUE    BND pair successfully rescued to INV candidate?
#   q7b_bnd_rescue_concordance numeric  STEP06  FALSE   carrier overlap between paired BND junctions
#
# FAILURE MODE 3: Filter chain losses
#   q7b_delly_site_in_raw     bool      4D      TRUE    site exists in raw merged BCF?
#   q7b_delly_site_passes_strict bool   4D      TRUE    site survives PRECISE+QUAL300+PE3?
#   q7b_delly_site_qual       numeric   4D      TRUE    actual QUAL value
#   q7b_delly_site_pe         int       4D      TRUE    actual PE count
#   q7b_manta_site_in_raw     bool      4G      TRUE    site exists in Manta split catalog?
#   q7b_manta_site_passes_pass bool     4G      TRUE    site survives PASS+QUAL20?
#   q7b_manta_site_qual       numeric   4G      TRUE    actual QUAL value
#
# FAILURE MODE 4: Population prior status
#   q7b_pop_prior_applied     bool      regeno  TRUE    population regenotyping run?
#   q7b_pop_prior_freq_est    numeric   regeno  FALSE   estimated allele frequency from confident calls
#   q7b_pop_prior_n_rescued   int       regeno  FALSE   samples rescued by posterior > 0.5


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIDENCE TIER ASSIGNMENT (pathway-based, replaces count-based)
# ═══════════════════════════════════════════════════════════════════════════════
#
# OLD: count how many of 12 dimensions pass → ≥8 = Tier 1
# NEW: which EVIDENCE PATHWAYS are satisfied?
#
# ┌─────────┬──────────────────────────────┬────────────────────────────────────┐
# │ Tier    │ Pathway                      │ Layers required                    │
# ├─────────┼──────────────────────────────┼────────────────────────────────────┤
# │ 1       │ Full convergence             │ A + B + C + D (4/4 tested, 4/4    │
# │         │                              │ passed, multi-family)              │
# │ 1       │ Strong convergence           │ A + B + D (3/3 tested, 3/3 passed │
# │         │                              │ — C not tested/available)          │
# ├─────────┼──────────────────────────────┼────────────────────────────────────┤
# │ 2       │ PCA + GHSL convergence       │ A + C strong, B absent,            │
# │         │                              │ D not testable                     │
# │ 2       │ PCA + SV partial             │ A + B present, D marginal          │
# │         │                              │ or C not tested                    │
# ├─────────┼──────────────────────────────┼────────────────────────────────────┤
# │ 3       │ Mixed evidence               │ Any 1 of {A,B,C} strong +          │
# │         │                              │ ≥4 D-dimensions positive           │
# ├─────────┼──────────────────────────────┼────────────────────────────────────┤
# │ 4       │ Weak / artifact              │ ≤1 layer, verdict = artifact,      │
# │         │                              │ or Cheat 25 DEAD                   │
# ├─────────┼──────────────────────────────┼────────────────────────────────────┤
# │ SV-1    │ SV-only (both callers)       │ B only (DELLY + Manta), no A or C  │
# │ SV-2    │ SV-only (single caller)      │ B partial, no A or C               │
# └─────────┴──────────────────────────────┴────────────────────────────────────┘
#
# KEY DESIGN PRINCIPLE: D7 = 0 no longer kills a candidate.
# A Tier 2 Pathway "PCA + GHSL" candidate has zero SV support BY DEFINITION.
# The PCA + GHSL agreement (two independent methods using different data
# and algorithms) is sufficient evidence without breakpoints.
#
# The 12 D-dimensions from C01d remain as ANNOTATION (the heatmap,
# the composite score, the per-dimension breakdown) but they no longer
# determine the tier. The tier comes from the 4-layer pathway assessment.


# ═══════════════════════════════════════════════════════════════════════════════
# THREE-AXIS STATUS SYSTEM
# ═══════════════════════════════════════════════════════════════════════════════
#
# Every candidate has three independent status values:
#
# AXIS 1: CONFIDENCE TIER (pathway-based, from above)
#   → "How sure are we this is real?"
#   → Tier 1 / 2 / 3 / 4 / SV-1 / SV-2
#
# AXIS 2: COMPLETION PERCENTAGE (from registry key counts)
#   → "How much do we know about it?"
#   → 0-100%, with per-question breakdown (Q1–Q7)
#   → Formula: resolved / (total - not_applicable) × 100
#   → 352 total keys possible
#
# AXIS 3: EVOLUTIONARY CLASS (from Q3+Q4+Q5+Q6 when sufficient keys)
#   → "What kind of inversion is it?"
#   → young_polymorphic / intermediate / old_polymorphic /
#     fixed_species_diagnostic / complex_nested / unresolved
#   → Only assigned when Q4 ≥ 30% AND Q5 ≥ 30% completion
#
# EXAMPLE OUTPUT:
#   LG12_inv17   Tier 1   78% complete   old_polymorphic     [A_full_convergence]
#   LG22_inv3    Tier 2   31% complete   young_polymorphic   [B_pca_ghsl_convergence]
#   LG01_inv42   Tier 3   44% complete   unresolved          [C_mixed]
#   LG06_inv8    SV-1     12% complete   unresolved          [E_sv_only_both_callers]
#
# The completion percentage is what drives the "37% resolved" display —
# it tells you exactly what's missing and what to run next.


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENTED FAILURE MODES (from 2026-04-16 pipeline audit)
# ═══════════════════════════════════════════════════════════════════════════════
#
# FM1: CALLABLE MASK BOUNDARY SHIFT
#   PCA windows exclude softmasked (atgc) regions. Inversion breakpoints
#   often sit in repeats. PCA boundaries are shifted inward by repeat width.
#   Consequences: recombinant undercount, truncated erosion gradient, age bias.
#   Keys affected: q3_left_bp, q3_right_bp (PCA estimates), q5_* (age).
#   Mitigation: use SV breakpoints (q3_left_sv_bp) for mechanism/age analyses.
#
# FM2: PER-SAMPLE STOCHASTIC DROPOUT AT 9×
#   At Poisson(λ=4.5) for HET carriers, P(≤3 reads) ≈ 34%.
#   DELLY/Manta genotype these as 0/0 without population prior.
#   Consequences: carrier count underestimate, allele freq bias, reduced
#   Fisher OR power, BND junction carrier discordance.
#   Keys affected: q7_layer_b_n_carriers, q6_freq_inv, q7_layer_d_fisher_or.
#   Mitigation: population regenotyping script; PCA carrier set as truth.
#
# FM3: INV vs BND FRAGMENTATION
#   Same inversion appears as INV in some samples, BND in others.
#   convertInversion.py needs both junctions per-sample → stochastic.
#   Post-conversion BND catalog has zero inversion BNDs.
#   Keys affected: q7_layer_b_detected, q7b_bnd_rescued.
#   Mitigation: STEP06 BND rescue from raw pre-conversion VCF.
#
# FM4: STRICT FILTER CHAIN LOSSES
#   DELLY: PRECISE=1 + QUAL≥300 + PE≥3. Manta: PASS + QUAL≥20.
#   Aggregate QUAL deflated by carrier dropout → borderline sites removed.
#   Keys affected: q7b_delly_site_passes_strict.
#   Mitigation: consult raw VCFs for D7 rescue in Tier 3-4 candidates.
#
# FM5: FIXED INVERSIONS INVISIBLE TO PCA
#   If all 226 samples carry the inversion, PCA has zero contrast.
#   SV callers detect breakpoints but no PCA core exists.
#   Keys affected: q7_layer_a_detected = FALSE for fixed inversions.
#   Tier: SV-1 or SV-2 (dedicated SV-only pathway).
#
# FM6: NESTED INVERSIONS
#   PCA sees one block; SV callers may report two calls with different
#   carrier sets. Greedy 1:1 matching misses the second.
#   Keys affected: q1_n_children, q1_landscape_category.
#   Mitigation: many-to-many matching, carrier concordance check.


# ═══════════════════════════════════════════════════════════════════════════════
# KEY COUNT SUMMARY (v2)
# ═══════════════════════════════════════════════════════════════════════════════
#
#   Q1: WHAT IS IT?                    49 keys
#   Q2: WHAT'S HAPPENING INSIDE IT?    40 keys
#   Q3: WHAT ARE THE BOUNDARIES DOING? 81 keys (+8 carrier reconciliation)
#   Q4: HOW DID IT FORM?              47 keys
#   Q5: HOW OLD IS IT?                39 keys
#   Q6: HOW COMMON IS IT?             28 keys
#   Q7: IS IT REAL?                   68 keys (+7 GHSL layer, +20 SV audit)
#   ─────────────────────────────────────────
#   TOTAL:                            352 keys
#
#   + 7 summary labels (one per question)
#   + 3 status axes (confidence tier, completion %, evolutionary class)
#   = 362 total evidence items per candidate
