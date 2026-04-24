# =============================================================================
# PHASE 4B REWRITE — Decomposition + Multi-Recombinant + Nested Composition
# =============================================================================
# Date: 2026-04-16
# Version: v10.1 (addendum to v10)
# Scope: Phase 4b of the inversion pipeline. Replaces the old single-script
#        C01i decomposition with a three-script coherent subsystem that
#        handles recombinants properly, detects composite intervals, and
#        flags sample-grouping-varies-across-subintervals cases.
#
# Depends on: v10 registry_loader, sample_registry, Engine B (instant_q
#             local-Q cache, optional but recommended)
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0. What this replaces and why
# ─────────────────────────────────────────────────────────────────────────────
#
# v9.3.4 had ONE script STEP_C01i_decomposition_rewired_24_v934_registry.R that:
#   - ran one PCA+kmeans on the whole candidate interval
#   - assigned each sample to HOM_REF/HET/HOM_INV/Recombinant (four labels,
#     but the Recombinant class was basically unused)
#   - registered 3 groups (missed RECOMBINANT entirely — see v10 patch)
#   - called cheat24 recombinant prior but didn't merge event_class back
#   - had no awareness of sub-intervals, composite structure, or the
#     "same interval, different sample groups" case
#
# The rewrite splits the work into three coherent scripts that run in
# parallel in Phase 4b and together produce a full decomposition:
#
#   4b.1  STEP_C01i_decompose.R           groups + per-window class refinement
#   4b.2  STEP_C01i_b_multi_recomb.R      recombinant detection with cheat24
#   4b.3  STEP_C01i_c_nested_composition.py  composite-interval flag
#
# All three write their outputs as Tier-2 structured blocks:
#   - internal_dynamics.json         (from C01i_decompose)
#   - recombinant_map.json           (from C01i_b_multi_recomb)
#   - internal_ancestry_composition.json (from C01i_c_nested_composition)
#
# A fourth script, STEP_C01i_d_seal.R, runs AFTER all three complete and
# produces the final group registrations + validation level for the candidate.
# This is where the four groups (HOM_REF, HET, HOM_INV, RECOMBINANT) plus
# the composite_flag get written to sample_registry.


# ─────────────────────────────────────────────────────────────────────────────
# 1. The three problems the rewrite solves
# ─────────────────────────────────────────────────────────────────────────────
#
# PROBLEM 1: Recombinants are a separate class, not bad HETs
#   v9.3.4 silently put recombinants in HET, which poisoned HET's Fst and
#   theta. The rewrite identifies them cleanly using three independent
#   signals (PCA residuals, phase switches from Clair3/WhatsHap, and
#   flashlight hemizygous segments if available) and splits them into
#   their own group. Cheat24 then classifies each recombinant as
#   gene_conversion, double_crossover, suspicious, or ambiguous.
#
# PROBLEM 2: Composite intervals look like one system but aren't
#   The current v10 catalog assumes one interval → one arrangement system.
#   If two inversions sit at overlapping coordinates with different family
#   memberships, the PCA k-means produces a hybrid grouping that is valid
#   for neither. The nested_composition script reads Engine B's local
#   ancestry labels and measures whether the internal ancestry structure
#   is homogeneous (one system), two_block (two systems), fragmented
#   (recombinant-rich or composite), or continuous_gradient (admixture).
#
# PROBLEM 3: Quality of decomposition is unaudited
#   v9.3.4 registered groups with no indication of how reliable they are.
#   The rewrite produces three auditable quality metrics per candidate:
#     - silhouette_score        (how well-separated are the three classes?)
#     - phase_concordance       (what fraction of HET samples have phased Clair3?)
#     - bic_gap                 (how much better is k=3 than k=2?)
#   These feed C01f's validation promotion — a low-quality decomposition
#   gets capped at UNCERTAIN.


# ─────────────────────────────────────────────────────────────────────────────
# 2. Architecture: the four scripts
# ─────────────────────────────────────────────────────────────────────────────
#
#   ┌─────────────────────────────────────────────────────────────────────┐
#   │  Phase 4a complete                                                  │
#   │  Tier ≤ 3 candidates in catalog                                     │
#   └───────┬─────────────────────────────────────────────────────────────┘
#           │
#           ├──────────────────────────┬────────────────────────┐
#           │                          │                        │
#           ▼                          ▼                        ▼
#   ┌───────────────┐        ┌──────────────────┐    ┌────────────────────┐
#   │ 4b.1          │        │ 4b.2             │    │ 4b.3               │
#   │ C01i_decompose│        │ C01i_b_multi_    │    │ C01i_c_nested_     │
#   │               │        │   recomb         │    │   composition      │
#   │ Per-candidate │        │ Per-candidate    │    │ Per-candidate      │
#   │ - PCA+kmeans  │        │ - Phase switches │    │ - Engine B local Q │
#   │   with        │        │   from Clair3    │    │ - Per-window       │
#   │   flashlight  │        │ - Hemi segments  │    │   ancestry labels  │
#   │   seeding     │        │   from flashlight│    │ - Block analysis   │
#   │ - Per-window  │        │ - cheat24        │    │ - Composite flag   │
#   │   class track │        │   classify       │    │                    │
#   │ - Silhouette  │        │ - Per-sample     │    │                    │
#   │ - BIC k=3 vs  │        │   event_class    │    │                    │
#   │   k=2         │        │                  │    │                    │
#   └───────┬───────┘        └────────┬─────────┘    └─────────┬──────────┘
#           │                         │                        │
#           │ internal_dynamics.json  │ recombinant_map.json   │ internal_ancestry
#           │                         │                        │ _composition.json
#           │                         │                        │
#           └───────┬─────────────────┴────────────────────────┘
#                   │
#                   ▼
#           ┌──────────────────┐
#           │ 4b.4             │
#           │ C01i_d_seal      │   SYNTHESIS step
#           │                  │
#           │ Reads all three  │   - merges per-sample labels across the 3 scripts
#           │ JSON blocks      │   - handles conflicts (e.g., PCA says HET,
#           │                  │     phase says recombinant → recombinant wins)
#           │ Registers 4      │   - produces FINAL per-sample class
#           │ groups +         │   - registers 4 groups in sample_registry
#           │ composite_flag + │   - sets q6_group_validation = UNCERTAIN
#           │ validation       │     (downgraded to SUSPECT if composite+fragmented,
#           │ level            │      capped at UNCERTAIN if composite_flag=likely)
#           └──────────────────┘
#
# Why split into three parallel scripts + synthesis rather than one monolithic:
#   1. Each script has one focused responsibility → easier to test + debug
#   2. Scripts 4b.1/4b.2/4b.3 can run in PARALLEL (SLURM array over candidates,
#      three jobs per candidate, ~3× wall-time speedup)
#   3. When something goes wrong for one candidate, you know from which script
#   4. Each produces its own Tier-2 block — observable, auditable, citable in
#      the manuscript methods


# ─────────────────────────────────────────────────────────────────────────────
# 3. Per-sample final class resolution (C01i_d_seal logic)
# ─────────────────────────────────────────────────────────────────────────────
#
# After the three parallel scripts complete, each sample has up to three
# labels:
#
#   pca_class       ∈ {HOM_REF, HET, HOM_INV}  from C01i_decompose
#   recomb_status   ∈ {NOT_RECOMB, recomb_GC, recomb_DCO,
#                        recomb_suspicious, recomb_ambiguous}
#                                                  from C01i_b_multi_recomb
#   structure_type  ∈ {homogeneous, dominant_plus_secondary,
#                        two_block_composite, continuous_gradient,
#                        multi_block_fragmented, diffuse_mixed}
#                                                  from C01i_c_nested_composition
#
# Resolution rules (applied per sample, in order):
#
#   Rule 1: If recomb_status != NOT_RECOMB → FINAL = RECOMBINANT
#           This overrides PCA. A sample with detected phase switch is
#           recombinant even if PCA clusters it as HET — that's the
#           entire reason we detect it.
#
#   Rule 2: If structure_type = multi_block_fragmented AND recomb_status =
#           NOT_RECOMB → FINAL = RECOMBINANT (weak), with
#           recomb_event_class = "ambiguous"
#           Phase data might be missing/noisy; ancestry fragmentation is
#           also a recombinant signal. Flag for manual review.
#
#   Rule 3: If structure_type = two_block_composite → keep pca_class but
#           flag sample as in_composite_region. These samples are in a
#           region where two systems overlap. Their class is correct for
#           whichever system we're calling on; the composite flag warns
#           downstream.
#
#   Rule 4: Otherwise → FINAL = pca_class
#
# Per-candidate composite flag (written to registry as q1_composite_flag):
#
#   clean        ≤ 5% of samples are two_block_composite AND
#                ≤ 5% are multi_block_fragmented
#   maybe_composite   5-20% two_block OR 5-20% multi_block
#   likely_composite  > 20% two_block OR > 20% multi_block
#
# Group registration (in sample_registry):
#
#   inv_<cid>_HOM_REF          samples with FINAL = HOM_REF
#   inv_<cid>_HET              samples with FINAL = HET
#   inv_<cid>_HOM_INV          samples with FINAL = HOM_INV
#   inv_<cid>_RECOMBINANT      samples with FINAL = RECOMBINANT
#
# Subgroups (registered only when ≥3 samples):
#
#   inv_<cid>_RECOMBINANT_GC   cheat24 event_class = gene_conversion
#   inv_<cid>_RECOMBINANT_DCO  cheat24 event_class = double_crossover
#
# Plus one alias for backward compat during migration:
#   inv_<cid>_HOM_STD          alias of inv_<cid>_HOM_REF
#
# Validation level after C01i_d_seal:
#
#   composite_flag = likely_composite  → UNCERTAIN (capped; C01f cannot promote)
#   composite_flag = maybe_composite   → UNCERTAIN (C01f can promote)
#   composite_flag = clean             → UNCERTAIN (C01f can promote)
#
#   silhouette_score < 0.25             → UNCERTAIN + quality_flag
#   bic_gap < 0.05 (k=3 barely beats k=2) → UNCERTAIN + 2-class_suspect


# ─────────────────────────────────────────────────────────────────────────────
# 4. Handling the "same interval, different sample groups" case
# ─────────────────────────────────────────────────────────────────────────────
#
# As discussed in the audit: the v10.1 catalog does NOT split composite
# intervals into separate system_ids. Instead it:
#
#   (a) Detects compositeness via nested_composition's structure_breakdown
#   (b) Flags the candidate with composite_flag = likely_composite
#   (c) Registers the (probably wrong) hybrid groups from C01i
#   (d) Caps q6_group_validation at UNCERTAIN so downstream Fst/burden/age
#       cheats refuse to run on these samples
#   (e) In classify_inversions, the final tag for a likely_composite
#       candidate ends with _composite_undecomposed and Q5/Q6 stay
#       "unclassified"
#
# This is the honest state: the system detects composite intervals, flags
# them explicitly, and prevents downstream analyses from using their
# unreliable groups. Full decomposition into multiple system_ids is a
# post-manuscript task.


# ─────────────────────────────────────────────────────────────────────────────
# 5. Engine B dependency
# ─────────────────────────────────────────────────────────────────────────────
#
# C01i_c_nested_composition REQUIRES Engine B's local-Q cache (instant_q
# output, one .local_Q_samples.tsv.gz per chromosome). This is the only
# part of the rewrite that genuinely depends on Engine B.
#
# Fallback behavior: if the cache is absent, C01i_c_nested_composition
# writes a stub block with composite_flag = "unknown_no_engine_b" and
# exits 0. Downstream treats unknown_no_engine_b the same as clean for
# validation-gating purposes (does NOT cap promotion), but classify_inversions
# records the missing information.
#
# C01i_decompose gets a soft bonus from flashlight (SV seeding + Cheat 2
# het-DEL constraint), but works without it. Falls back to standard
# k-means if flashlight is unavailable.
#
# C01i_b_multi_recomb works with Clair3 phase blocks alone (no Engine B).
# Flashlight hemizygous segments are an additional signal when available.


# ─────────────────────────────────────────────────────────────────────────────
# 6. File manifest
# ─────────────────────────────────────────────────────────────────────────────
#
# R/
#   STEP_C01i_decompose.R              ~600 lines — replaces v9.3.4 C01i
#   STEP_C01i_b_multi_recomb.R         ~450 lines — recombinant detection
#   STEP_C01i_d_seal.R                 ~300 lines — synthesis + register
#   lib_decompose_helpers.R            ~200 lines — shared utilities
#
# python/
#   STEP_C01i_c_nested_composition.py  ~250 lines — thin wrapper around
#                                         MODULE_2B's nested_composition
#                                         that writes to registry blocks
#   nested_composition_core.py         271 lines — the uploaded script,
#                                         vendored unchanged as a library
#
# schemas/
#   internal_dynamics.schema.json         updated with recombinant list
#   recombinant_map.schema.json           new
#   internal_ancestry_composition.schema.json  new
#
# patches/
#   C01f_composite_flag_gate.R         ~30 lines — one rule in
#                                         compute_group_validation
#
# orchestrator/
#   run_phase4b.sh                     SLURM chain for the 4 scripts
#
# tests/
#   test_c01i_d_seal_resolution.py     resolution rules unit tests
#   test_phase4b_integration.py        end-to-end on synthetic data
