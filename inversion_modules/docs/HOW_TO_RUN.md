# =============================================================================
# HOW TO RUN — Inversion Pipeline v8.5
# =============================================================================
#
# Prerequisites:
#   - All 226 QC-passing BAMs mapped to fClaHyb_Gar_LG.fa (28 chr)
#   - Dosage matrices in 04_dosage_by_chr/
#   - MDS precomputed in 06_mds_candidates/ (from STEP_B01)
#   - conda env "assembly" with R 4.x, data.table, ggplot2, etc.
#
# Config: 00_inversion_config.sh
#   - Source at top of every script/launcher
#   - All paths resolve from BASE
#   - RSCRIPT_BIN set to conda env Rscript
#
# =============================================================================
# FULL PIPELINE RUN ORDER
# =============================================================================
#
# ── PHASE_01: LANDSCAPE (run once) ──────────────────────────────────
#
# Step 1: Precompute (eigenvalues, z-scores, sim_mat, MDS, inv-likeness)
sbatch launchers/LAUNCH_C01a_snake1_precompute.slurm

# Step 2: Precompute diagnostics (17 plot types per chr + genome-wide)
sbatch launchers/LAUNCH_C01b_snake1.slurm --step diag_precomp

# Step 3: Block detection + boundary classification + blue-cross
sbatch launchers/LAUNCH_01C_block_detect.slurm
#   Optional: single chromosome:
#   sbatch launchers/LAUNCH_01C_block_detect.slurm C_gar_LG28

# Step 4 (optional): Network diagnostics
sbatch launchers/LAUNCH_C01b_snake1.slurm --step network_diag

#
# ── PHASE_02: DETECTION ─────────────────────────────────────────────
#
# Step 5: Snake cores (S1S/S1M/S1L with sharp-drop + PA)
sbatch launchers/LAUNCH_C01b_snake1.slurm --step cores

# Step 6: Fuzzy merge (uses landscape data from step 3 automatically)
sbatch launchers/LAUNCH_C01b_snake1.slurm --step merge

#
# ── PHASE_03: SCORING ───────────────────────────────────────────────
#
# Step 8: Candidate scoring (10D + PA aggregation)
sbatch launchers/LAUNCH_C01d_scoring.slurm

#
# ── PHASE_04: CHARACTERIZATION ──────────────────────────────────────
#
# Step 9: Regime compatibility engine
sbatch launchers/LAUNCH_C01j_regime_engine.slurm

#
# ── SHORTCUTS ───────────────────────────────────────────────────────
#
# Run cores + merge + diag in one job:
sbatch launchers/LAUNCH_C01b_snake1.slurm --step cores,merge,diag

# Run merge only with a tag (subfolder):
sbatch launchers/LAUNCH_C01b_snake1.slurm --step merge --tag relaxed_v2

#
# ── ITERATING ───────────────────────────────────────────────────────
#
# After changing parameters in C01b_1 or C01b_2:
#   - Steps 1-3 don't need rerunning (landscape is stable)
#   - Rerun from step 5 (cores) onward
#   - Use --tag to save alternative parameter sets
#
# After changing 01C parameters (block thresholds):
#   - Rerun from step 3 onward
#
# =============================================================================
