# =============================================================================
# INTEGRATION AUDIT — Detector v9.3 + Flashlight v3 + Codebase v8.5.9
# =============================================================================
#
# =============================================================================
# WHAT REPLACES WHAT
# =============================================================================
#
# OLD FLOW (codebase v8.5.9):
#   C01a precompute → C01b snake → C01c triangles → C01d scoring → C01e-m
#                                   ^^^^^^^^^^^
#                                   OUTDATED
#                                   chromosome-wide insulation = noisy
#                                   (Session 6 finding)
#
# NEW FLOW (v9.3 replaces C01c only):
#   C01a precompute → C01b snake → inv_detect_v9.3 → 12_bridge → C01d → C01e-m
#                                  ^^^^^^^^^^^^^^^   ^^^^^^^^^
#                                  REPLACES C01c     format adapter
#                                  9 phases:         produces
#                                  staircase +       triangle_intervals.tsv.gz
#                                  cheats + NN +     in C01c format so C01d
#                                  consensus +       reads it unchanged
#                                  peeling
#
# EVERYTHING AFTER C01d IS UNCHANGED:
#   C01d  candidate_scoring    ← reads triangle_intervals.tsv.gz (same format)
#   C01e  candidate_figures    ← reads scored candidates (same format)
#   C01f  hypothesis_tests     ← unchanged
#   C01g  boundary_refinement  ← unchanged (flashlight)
#   C01h  recombinant_scanner  ← unchanged
#   C01i  decomposition        ← unchanged (reads scored candidates)
#   C01j  regime_engine        ← unchanged
#   C01k  annotated_simmat     ← unchanged
#   C01l  local_structure      ← unchanged
#   C01m  distance_concordance ← unchanged
#   C01n  peeling_diagnostic   ← NEW (runs inside inv_detect or standalone)
#
# =============================================================================
# CONCRETE COMMANDS TO SWITCH
# =============================================================================
#
# BEFORE (old C01c):
#   Rscript STEP_C01c_triangle_regimes.R precomp/ triangle_out/
#   Rscript STEP_C01d_candidate_scoring.R triangle_out/ scoring_out/
#
# AFTER (new inv_detect → bridge → C01d):
#   # Step 1: Run detector (replaces C01c)
#   Rscript inv_detect_v9.3/00_run_all.R \
#     --chr C_gar_LG01 --sim-mat-dir precomp/ --precomp-dir precomp/ \
#     --outdir inv_detect_out_v9.3/
#
#   # Step 2: Bridge to codebase format
#   Rscript inv_detect_v9.3/12_bridge_to_codebase.R \
#     --detector_dir inv_detect_out_v9.3/ --outdir triangle_from_detector/
#
#   # Step 3: Score (unchanged command, different input dir)
#   Rscript STEP_C01d_candidate_scoring.R triangle_from_detector/ scoring_out/
#
# THAT'S IT. C01d reads the same file format. Everything downstream unchanged.
#
# =============================================================================
# FULL PRODUCTION PIPELINE (recommended order)
# =============================================================================
#
# # 0. One-time: build per-chr pruning table (~5 min)
# sbatch inv_detect_v9.3/11a_run_ngsrelate_perchr.sh
# Rscript inv_detect_v9.3/11b_build_perchr_pruning_table.R \
#   --ngsrelate_dir ngsrelate_perchr/ --outfile per_chr_pruned.tsv
#
# # 1. Precompute (if not already done)
# Rscript STEP_C01a_snake1_precompute.R ...
#
# # 2. Snake regions (if not already done)
# Rscript STEP_C01b_1_cores.R ...
#
# # 3. Detector (REPLACES C01c) — all 28 chr
# sbatch inv_detect_v9.3/slurm_run.sh
#
# # 4. Bridge to codebase format
# Rscript inv_detect_v9.3/12_bridge_to_codebase.R \
#   --detector_dir inv_detect_out_v9.3/ --outdir triangle_from_detector/
#
# # 5. Score + downstream (UNCHANGED from here on)
# Rscript STEP_C01d_candidate_scoring.R triangle_from_detector/ scoring_out/
# Rscript STEP_C01i ... --scores scoring_out/...
# Rscript STEP_C01j ... --scores scoring_out/...
#
# =============================================================================
# RELATEDNESS / PEELING CONFOUND (honest discussion)
# =============================================================================
#
# ngsRelate on a chromosome with a big inversion inflates kinship between
# carriers. Two unrelated fish both HOM_INV share ~80K concordant sites →
# ngsRelate sees them as relatives. Per-chr ngsRelate prunes them.
#
# For pop-gen this is a problem. For the peeling diagnostic it's a feature:
#   - L1 (ngsRelate pruning): removes "relatives" which may be carriers.
#     If block disappears → those fish were essential. Family or carrier?
#   - L1b (PC1 chr-local): independent method, same question.
#   - L2 (block-level): tiebreaker. If L2=stable but L1=disappeared →
#     ngsRelate pruned inversion carriers, not family.
#     If L2 also disappeared → genuinely family-driven.
#
# First-degree pruning doesn't have to be perfect. It may remove inversion
# carriers. That's OK — it's a diagnostic view, not a filter.
#
# =============================================================================
# CROSS-VALIDATION TABLE
# =============================================================================
#
# | Source           | Signal                       | Weight |
# |------------------|------------------------------|--------|
# | Detector v9.3    | nn_birth >= 80               | +++    |
# | Detector v9.3    | shape_class = strong_square   | ++     |
# | Detector v9.3    | consensus >= 3 variants       | ++     |
# | Peel diagnostic  | L1b chr-local peel = stable   | +++    |
# | Peel diagnostic  | L2 block coseg peel = stable  | ++     |
# | Flashlight       | cheat25 viability = PASS      | ++     |
# | Flashlight       | cheat26 kin-pruned retention   | +++    |
# | Codebase         | C01i decomposition found      | ++     |
# | Codebase         | C01j regime = structured_inv  | ++     |
# | SV evidence      | DELLY/Manta breakpoint match  | +++    |
#
# A candidate with 5+ positive signals = very high confidence.
#
# =============================================================================
