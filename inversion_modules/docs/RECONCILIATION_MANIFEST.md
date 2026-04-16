# =============================================================================
# v8.4 → v8.5 RECONCILIATION MANIFEST
# =============================================================================
# Updated: 2026-04-07 (second pass — all v8.5 items complete)
#
# HOW TO USE:
#   1. Scripts marked UPDATED: replace the HPC version with the tarball version
#   2. Scripts marked NEW: add to the snakes/ directory
#   3. Scripts marked UNCHANGED: keep the HPC version as-is
#   4. Scripts marked DELETE: remove .old/.deprecated/.bak files
#   5. Scripts marked SUPERSEDED: C01a3 + C01a4 are merged into PHASE_01C
#
# =============================================================================

# === UPDATED (replace HPC version) ==========================================

# STEP_C01a_snake1_precompute.R  (341 → 418 lines)
#   CHANGES:
#   - Inv-likeness: relative PVE1 scoring (excess over chr median)
#   - Inv-likeness: relative eigenvalue ratio scoring
#   - NEW dimension: het_contrast (between/within variance ratio from k=3)
#   - NEW formula: 0.25*pve + 0.20*eig + 0.25*dip + 0.30*het (was 0.35/0.25/0.40)
#   - Hatchery note: PVE1 weight should be further reduced (see REMAINING)
#   - NEW output columns: inv_het_contrast, inv_pve1_excess
#   STATUS: ✅ READY

# STEP_C01a_diag.R  (1190 → ~1310 lines)
#   CHANGES (v8.5):
#   - NEW plot 16: sim_mat row-wise Z-score contrast transform
#   - NEW plot 17: sim_mat row quantile-normalization contrast transform
#   - Both produce per-chr PDF + PNG, same render loop as other plots
#   - Header updated, version tag v8.5
#   STATUS: ✅ READY

# STEP_C01b_1_cores.R  (605 → 822 lines)
#   CHANGES:
#   - Roary-style PA matrix: per-window per-family
#   - Sharp-drop boundary detector
#   - Split heatmap A1
#   STATUS: ✅ READY

# STEP_C01b_2_merge.R  (1145 → ~1280 lines)
#   CHANGES (v8.5):
#   - NEW: --landscape_dir CLI argument
#   - NEW: loads block_concordance, blue_cross_verdicts, boundary_catalog from 01C
#   - NEW: landscape_adjust() function applies concordance gates, blue-cross
#     boost/penalty, boundary-type penalty to fuzzy merge scores
#   - Score log now includes: landscape_adj, landscape_reason, adjusted_score
#   - Summary log reports landscape-adjusted pair counts
#   - Backward-compatible: runs in v8.4 mode without --landscape_dir
#   - Version tag v8.5
#   STATUS: ✅ READY (was ⚠️ PARTIAL in v8.4)

# STEP_C01c_triangle_regimes.R  (1448 lines)
#   No changes from v8.4→v8.5 session
#   STATUS: ✅ READY

# STEP_C01d_candidate_scoring.R  (603 lines)
#   No changes from v8.4→v8.5 session
#   STATUS: ✅ READY

# STEP_C01i_multi_inversion_decomposition.R  (885 lines)
#   No changes from v8.4→v8.5 session
#   STATUS: ✅ READY

# STEP_C01j_regime_compatibility_engine.R  (927 lines)
#   No changes from v8.4→v8.5 session
#   STATUS: ✅ READY

# === NEW SCRIPTS =============================================================

# PHASE_01C_block_detect.R  (494 lines)
#   Merges C01a3 + C01a4 + boundary classification into one script.
#   STATUS: ✅ READY

# === NEW SLURM LAUNCHERS (v8.5) =============================================

# LAUNCH_01C_block_detect.slurm       — PHASE_01C landscape
# LAUNCH_C01c_triangles.slurm         — PHASE_02C triangles
# LAUNCH_C01d_scoring.slurm           — PHASE_03A scoring
# LAUNCH_C01j_regime_engine.slurm     — PHASE_04C regimes
#
# All launchers:
#   - Use ${RSCRIPT_BIN} from config (not hardcoded)
#   - Source config with set -a / set +a
#   - set -euo pipefail
#   - mamba activate assembly
#   - All existing launchers updated: config path → v8.5
#   STATUS: ✅ READY

# === UPDATED CONFIG ==========================================================

# 00_inversion_config.sh
#   CHANGES (v8.5):
#   - CODEBASE path: inversion_codebase_v8.4 → inversion_codebase_v8.5
#   - NEW dirs: LANDSCAPE_DIR, TRIANGLE_DIR, SCORING_DIR,
#     DECOMPOSITION_DIR, REGIME_DIR
#   - inv_init_dirs() creates new directories
#   STATUS: ✅ READY

# === UNCHANGED (keep HPC version) ===========================================

# STEP_C01a2_network_diag.R, STEP_C01b_diag.R, STEP_C02_snake1_diag.R,
# STEP_C02b_snake_composite_debug.R, STEP_C03_snake2_community.R,
# STEP_C04_snake3_ghsl.R, STEP_C05_snake4_threeband.R,
# STEP_C06_consensus.R, STEP_C07_hmm_regime.R, STEP_C08_rescue_pass2.R
#
# Also keep: STEP_C01e, C01f, C01f_b, C01g, C01h, C01k,
# rare_inversion_score.R, fuzzy_merge_score.R

# === SUPERSEDED (keep as reference) =========================================

# STEP_C01a3_blue_cross_diagnostic.R  → merged into PHASE_01C
# STEP_C01a4_block_comparison_diagnostic.R → merged into PHASE_01C

# === DELETE (old/deprecated files on HPC) ====================================

# All *.old, *.old.old, *.bak*, *.deprecated, *.neededit files

# =============================================================================
# REMAINING WORK (not done yet)
# =============================================================================

# ❌ PVE1 in hatchery mode: reduce weight further (0.10 or 0.0).
#    Het_contrast breaks for rare inversions (<10% carriers) because k=3
#    doesn't produce meaningful groups. Need fallback: when smallest k=3
#    group < 5% of samples, use eigenvalue ratio only.

# ❌ Integrate STEP20-41 (followup scripts) into PHASE_04/05 numbering
#    See LEGACY_FOLLOWUP_BRIDGE.md for priorities.

# ❌ Method name for repo: proposed "SnakeInv" or "Mamba"

# =============================================================================
# COMPLETED IN THIS SESSION (v8.5)
# =============================================================================

# ✅ C01b_2 merge: landscape integration (concordance + blue_cross + boundary)
# ✅ Sim_mat contrast transform: row Z-score + quantile (C01a_diag plots 16/17)
# ✅ SLURM launchers: 01C, C01c, C01d, C01j (4 new launchers)
# ✅ All 18 existing launchers: config path updated to v8.5
# ✅ Config: CODEBASE=v8.5, new output dirs, inv_init_dirs updated
# ✅ Docs: PIPELINE_v8.5_ARCHITECTURE.md + RECONCILIATION_MANIFEST.md updated
# =============================================================================
