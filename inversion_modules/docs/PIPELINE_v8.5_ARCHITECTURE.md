# =============================================================================
# INVERSION PIPELINE v8.5 — REORGANIZED ARCHITECTURE
# =============================================================================
#
# RATIONALE: v8.4 accumulated 15+ scripts with cross-dependencies. v8.5
# consolidates into logical phases with clear inputs/outputs.
#
# KEY CHANGES FROM v8.4:
#   - Merged C01a3 (blue cross) + C01a4 (block comparison) into PHASE_01
#   - Added boundary classification (clear/soft/inner) to PHASE_01
#   - Block concordance and blue-cross verdicts consumed by PHASE_02 and PHASE_03
#   - PA matrices standardized: every phase writes <phase>_window_pa.tsv.gz
#   - Split heatmaps standardized: every phase writes <chr>_<phase>_split.png
#   - All-layers strip plot runs as standalone at the end (not embedded in C01j)
#
# v8.5 ADDITIONS (this release):
#   - C01b_2 merge now consumes landscape data (block_concordance,
#     blue_cross_verdicts, boundary_catalog) via --landscape_dir
#   - Landscape adjustment modifies fuzzy merge scores: concordance gates,
#     blue-cross boost/penalty, boundary-type penalty
#   - Sim_mat contrast transforms added to C01a_diag (plots 16 & 17):
#     row-wise Z-score and row quantile-normalization
#   - SLURM launchers for all v8.5 scripts (01C, C01c, C01d, C01j)
#   - All launchers updated to v8.5 codebase path
#   - Config updated: CODEBASE=v8.5, new output dirs added
#
# =============================================================================
# PHASE OVERVIEW
# =============================================================================
#
# PHASE_01: LANDSCAPE (run once, reuse many times)
#   01A_precompute.R         — eigenvalues, z-scores, sim_mat, MDS, inv-likeness
#   01A_diag.R               — 17 diagnostic plots per chr + genome-wide summaries
#                               (v8.5: +plot 16 sim_mat Z-score, +plot 17 quantile)
#   01C_block_detect.R       — detect blocks, classify boundaries, blue-cross
#                               test, block concordance matrix, sample grouping
#   01D_network_diag.R       — sample network per block (ancestry + relatedness)
#
#   Outputs consumed by all downstream phases:
#     precomp/<chr>.precomp.rds
#     landscape/block_registry_<chr>.tsv.gz      — blocks with boundaries + types
#     landscape/block_concordance_<chr>.tsv.gz    — block-vs-block concordance
#     landscape/blue_cross_verdicts_<chr>.tsv.gz  — assembly error vs real boundary
#     landscape/boundary_catalog_<chr>.tsv.gz     — classified boundaries
#     landscape/01C_window_pa.tsv.gz              — block PA per window
#
# PHASE_02: DETECTION (snake cores + merge, uses landscape)
#   02A_snake_cores.R        — S1S/S1M/S1L cores with sharp-drop, PA matrix
#   02B_fuzzy_merge.R        — 3-tier merge, uses block concordance + blue-cross
#                               verdicts to decide merge/split (v8.5: IMPLEMENTED)
#   02C_triangles.R          — insulation boundaries, squareness, hatchery mode
#
#   Outputs:
#     detection/02A_window_pa.tsv.gz  — core PA
#     detection/02B_window_pa.tsv.gz  — merge PA
#     detection/02C_window_pa.tsv.gz  — triangle PA
#     detection/<chr>_D1_split_cores.png
#     detection/<chr>_D2_split_merge.png
#     detection/<chr>_D3_split_triangles.png
#
# PHASE_03: SCORING + GATING (first pass)
#   03A_candidate_scoring.R  — 10 dimensions + PA aggregation, tier assignment
#
#   Outputs:
#     scoring/candidate_scores.tsv.gz
#     scoring/<chr>_S1_evidence_heatmap.png
#
# PHASE_04: CHARACTERIZATION (deep analysis on Tier 1-3)
#   04A_decomposition.R      — marker co-segregation, core packages, extension
#   04B_recombinant_scanner.R — het-rate trajectory from VCF
#   04C_regime_engine.R      — compatibility grouping, sample types
#   04D_hypothesis_tests.R   — 8 statistical tests
#   04E_negative_controls.R  — null distributions + ancestry instability
#   04F_tube_graph.R         — node/edge topology + GHSL integration
#
#   Outputs:
#     characterization/04C_window_pa.tsv.gz  — regime PA
#     characterization/ancestral_fragments.tsv.gz
#     characterization/regime_state_labels.tsv.gz
#
# PHASE_05: CONFIRMATION + FIGURES
#   05A_candidate_scoring.R  — second pass with all 10 dimensions + PA
#   05B_candidate_figures.R  — per-candidate publication panels
#   05C_annotated_simmat.R   — 6-layer overlay on raw sim_mat
#   05D_all_layers_strip.R   — combined strip plot (standalone, loads all PAs)
#
#   Outputs:
#     figures/<chr>_F1_all_layers.png
#     figures/<chr>_<cand>_composite.png
#
# =============================================================================
# RUN ORDER (SLURM)
# =============================================================================
#
#   1. sbatch LAUNCH_C01a_snake1_precompute.slurm       # PHASE_01A (once)
#   2. sbatch LAUNCH_C01b_snake1.slurm --step diag_precomp  # 01A diagnostics
#   3. sbatch LAUNCH_01C_block_detect.slurm             # PHASE_01C landscape
#   4. sbatch LAUNCH_C01b_snake1.slurm --step cores     # PHASE_02A cores
#   5. sbatch LAUNCH_C01b_snake1.slurm --step merge     # PHASE_02B merge (+landscape)
#   6. sbatch LAUNCH_C01c_triangles.slurm               # PHASE_02C triangles
#   7. sbatch LAUNCH_C01d_scoring.slurm                 # PHASE_03A scoring
#   8. sbatch LAUNCH_C01j_regime_engine.slurm           # PHASE_04C regimes
#
# =============================================================================
# BOUNDARY CLASSIFICATION TYPES (computed in 01C)
# =============================================================================
#
# Each boundary between adjacent blocks or at block edges is classified:
#
#   clear_hard    — sim drops from >0.7 to <0.3 in ≤3 bins. Sharp structural
#                   boundary. High confidence inversion breakpoint.
#
#   clear_soft    — sim drops from >0.6 to <0.4 over 4-10 bins. Gradual
#                   transition. Could be: partial recombination zone near
#                   breakpoint, or family LD fadeout.
#
#   inner_hard    — narrow blue stripe (≤5 bins) inside a block. Both sides
#                   have high similarity. Either assembly error (if high
#                   concordance across) or sub-system transition.
#
#   inner_soft    — broader low-similarity zone (6-20 bins) inside a block.
#                   Regime change or overlapping system boundary.
#
#   diffuse       — no clear boundary, similarity gradually changes over
#                   >20 bins. Family LD transition, not structural.
#
#   chromosome_edge — block at start/end of chromosome, no comparison possible.
#
# =============================================================================
# MERGE DECISION RULES (02B uses landscape outputs) — v8.5 IMPLEMENTED
# =============================================================================
#
# When considering merging two adjacent cores:
#
#   1. Check block_concordance: if the two cores fall in blocks with
#      concordance > 0.70, they share carriers → boost score +0.05.
#      If concordance < 0.40, different systems → penalize -0.25.
#
#   2. Check blue_cross_verdicts: if the gap between cores contains a
#      blue cross with verdict "assembly_error" or "same_system",
#      boost merge score by +0.15.
#      If verdict is "different_systems", penalize by -0.20.
#
#   3. Check boundary_catalog: if the boundary between cores is
#      "clear_hard", penalize merge by -0.15 (real breakpoint).
#      If "inner_hard" + "same_system", boost by +0.10
#      (assembly artifact inside one system).
#
#   All adjustments are additive on the raw fuzzy score. Final score
#   is clamped to [0, 1]. Landscape data is optional: without
#   --landscape_dir, merge runs in v8.4 compatibility mode.
#
# =============================================================================
# SIM_MAT CONTRAST TRANSFORMS (v8.5, diagnostic addon in C01a_diag)
# =============================================================================
#
# Two additional views alongside the raw sim_mat:
#
#   Plot 16: Row-wise Z-score
#     For each row i: z[i,j] = (sim[i,j] - mean(row_i)) / sd(row_i)
#     Then symmetrised: (Z + Z') / 2
#     Removes family baseline inflation. Blocks where nearby windows
#     have higher-than-expected similarity stand out in red.
#     Clipped at ±3 SD for display.
#
#   Plot 17: Row quantile-normalization
#     For each row i: replace sim values with within-row rank / n
#     Then symmetrised. Produces [0,1] matrix.
#     Shows relative similarity structure free of global baseline.
#
#   Both are DIAGNOSTIC views, not pipeline modifications. The raw
#   sim_mat remains the input to all downstream computation.
#
# =============================================================================
