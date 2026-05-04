# HOW TO RUN — Inversion Detector v9.3

## What changed from v9.1

| Component | v9.1 (Session 7) | v9.3 (Session 8) |
|-----------|-------------------|-------------------|
| Staircase | Global decay profile | Per-window left/right traversal |
| Boundaries | Any drop = boundary | Sustained drops only (persist_n=3) |
| Artifacts | Discarded | Separate artifact registry |
| Matrix cheats | Not implemented | A1-A6 all implemented |
| Bloc scoring | Not implemented | 6 metrics + shape classifier |
| NN module | Check nn0 blocks in nn40 | Independent staircase per NN |
| NN tree | Not implemented | Full interval tree with nn_birth |
| Consensus | Not implemented | Cross-variant matching |

## Architecture overview

```
Phase 1: Staircase on raw sim_mat → block registry + artifact registry
Phase 2: Matrix cheats → 6 treated variants (distcorr, localnorm, denoised, resid_bg, edge, support)
Phase 3: Bloc scoring → inside_mean, squareness, sharpness, occupancy, patchiness per variant
Phase 4: NN persistence → independent staircase at nn0/20/40/80, match registries
Phase 5: NN sweep tree → full interval tree with nn_birth classification (optional)
Phase 6: Consensus → match blocks across variants, confidence = HIGH/MODERATE/LOW
Phase 7: Evidence stack → flatness, CV, GHSL, SV overlap, local PCA/ICA
Phase 8: Final scoring table → everything joined, calibration check
Phase 9: Peeling diagnostic → L1/L1b/L2/L3 sample-zeroing, family vs structure test
```

## Setup

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/snakes

# Unpack (or copy directory)
cp -r inv_detect_v9.3 .
# OR: tar xzf inv_detect_v9.3.tar.gz

conda activate assembly
```

## Step 1: Quick test on LG01 (staircase only)

This is the FIRST thing to run.  No matrix cheats, fast.

```bash
sbatch inv_detect_v9.3/slurm_test_LG01.sh
# ~5-10 min.  Check:
tail -f logs/inv_test_LG01.out
```

Or interactively:

```r
setwd("inv_detect_v9.3")
source("00_config.R")
source("01_staircase_detector.R")

# Load LG01 nn0 sim_mat
pc <- readRDS("../precomp/C_gar_LG01.precomp.rds")
smat <- pc$sim_mat

# Run staircase
result <- detect_blocks_staircase(smat)

# Inspect
print(result$blocks[, .(block_id, start, end, width, start_mb, end_mb,
                         height, step_down_left, step_down_right, n_artifacts)])

# Does it find I57? (24-28 Mb → bins ~480-560)
result$blocks[start_mb >= 23 & end_mb <= 29]

# Does it find I70? (29-33.5 Mb → bins ~580-670)
result$blocks[start_mb >= 28 & end_mb <= 34]

# I86 should be absent or have low height/step_down
result$blocks[start_mb >= 33 & end_mb <= 39]

# Quick diagnostic plot
plot(result$profiles$left_extent + result$profiles$right_extent,
     type="l", main="Total extent (block signal)", xlab="Window", ylab="Extent")
```

### What to look for

- **I57 (24-28 Mb):** Should show as 1 large block with nested sub-blocks.
  Height should be high (>0.15), step_down should be sharp (>0.05).
- **I70 (29-33.5 Mb):** Should show as 1 clean block, homogeneous.
- **I86 (34-38 Mb):** Should be ABSENT, or present as a weak block with
  low height and small step_down.  shape_class should be diagonal_band.
- **Total blocks:** If >200, increase min_block_width or min_drop.
  If <5, decrease them.

### Tuning

```r
# Too many blocks → increase thresholds
CFG$STAIR_MIN_BLOCK_WIDTH <- 10L
CFG$STAIR_MIN_DROP <- 0.05
result <- detect_blocks_staircase(smat)

# Too few blocks → decrease thresholds
CFG$STAIR_MIN_BLOCK_WIDTH <- 3L
CFG$STAIR_MIN_DROP <- 0.02
result <- detect_blocks_staircase(smat)

# Sustained-drop tolerance (how many bins must stay below)
CFG$STAIR_PERSIST_N <- 5L   # stricter: fewer boundaries
CFG$STAIR_PERSIST_N <- 2L   # looser: more boundaries
```


## Step 2: Add matrix cheats + bloc scoring

```r
CHR <- "C_gar_LG01"
SIM_MAT_DIR <- "../precomp"
OUTDIR <- "../inv_detect_out_v9.3"
PHASES <- "1:3"
source("00_run_all.R")

# Check bloc scores
scores <- fread(file.path(OUTDIR, "bloc_scores_C_gar_LG01.tsv"))
print(scores[variant == "raw", .(block_id, contrast, squareness,
                                  occupancy, patchiness, shape_class)])

# Compare across variants
dcast(scores, block_id ~ variant, value.var = "contrast")
```


## Step 3: Add NN persistence

```r
PHASES <- "1:4"
source("00_run_all.R")

nn <- fread(file.path(OUTDIR, "nn_persistence_C_gar_LG01.tsv"))
print(nn[, .(candidate_id, ref_start, ref_end, survives_nn40,
             survives_nn80, nn_topology)])
```


## Step 4: Full pipeline on LG01

```r
PHASES <- "all"
SKIP_CHEATS <- FALSE
SKIP_TREE <- FALSE
source("00_run_all.R")

# Full scoring table
tab <- fread(file.path(OUTDIR, "scoring_table_C_gar_LG01.tsv"))
# Sort by strongest candidates
tab[order(-squareness, -contrast)]
```


## Step 5: Run all 28 chromosomes

Only after LG01 looks good.

```bash
sbatch inv_detect_v9.3/slurm_run.sh
squeue -u $USER
```


## Step 6: Genome-wide merge

```r
files <- list.files("../inv_detect_out_v9.3",
                    pattern="^scoring_table_.*\\.tsv$", full.names=TRUE)
all_tabs <- lapply(files, fread)
genome <- rbindlist(all_tabs, fill=TRUE)

# Sort by strength
genome[order(-squareness, -contrast)][1:30,
  .(chr, start_mb, end_mb, width, shape_class, squareness,
    contrast, survives_nn40, survives_nn80)]
```


## File inventory after full run

```
inv_detect_out_v9.3/
  blocks_C_gar_LG01.tsv          ← Phase 1: block registry
  artifacts_C_gar_LG01.tsv       ← Phase 1: artifact registry
  staircase_C_gar_LG01.rds       ← Phase 1: full result (profiles)
  variants_C_gar_LG01.rds        ← Phase 2: all matrix variants
  bloc_scores_C_gar_LG01.tsv     ← Phase 3: scores per variant
  nn_persistence_C_gar_LG01.tsv  ← Phase 4: NN survival
  nn_tree_C_gar_LG01.tsv         ← Phase 5: interval tree
  consensus_C_gar_LG01.tsv       ← Phase 6: cross-variant consensus
  ev_flatness_C_gar_LG01.tsv     ← Phase 7a: decay shape
  ev_cv_C_gar_LG01.tsv           ← Phase 7b: interior uniformity
  ev_ghsl_C_gar_LG01.tsv         ← Phase 7c: partition stability
  ev_sv_C_gar_LG01.tsv           ← Phase 7d: SV breakpoint matching
  ev_pca_C_gar_LG01.tsv          ← Phase 7e: local PCA/ICA
  scoring_table_C_gar_LG01.tsv   ← Phase 8: everything joined (incl. peel cols)
  peel_diagnostic_C_gar_LG01.tsv ← Phase 9: before/after peeling per block
  peel_groups_C_gar_LG01.tsv     ← Phase 9: which samples zeroed per block
```


## Step 0 (optional): Build per-chromosome pruning table

Only needs to run once. Uses precomp PC1 to identify chromosome-local
relatives without rerunning ngsRelate. Different chromosomes get different
prune sets (LG01 might keep 186, LG02 might keep 200).

```bash
Rscript inv_detect_v9.3/11_build_perchr_pruning_table.R \
  --precomp_dir precomp/ \
  --outfile per_chr_pruned.tsv \
  --cor_threshold 0.7

# Output: per_chr_pruned.tsv with columns: chr | sample | status
# Pass to --pruned-list in the main pipeline
```


## All configurable parameters

All live in `00_config.R`.  Nothing hidden.

| Parameter | Default | Module | Effect |
|-----------|---------|--------|--------|
| STAIR_MIN_BLOCK_WIDTH | 5 | 01 | Min block size in bins |
| STAIR_MIN_DROP | 0.03 | 01 | Min step height for boundary |
| STAIR_PERSIST_N | 3 | 01 | Sustained drop must persist N bins |
| STAIR_SMOOTH_SPAN | 3 | 01 | Running median span |
| CHEAT_LOCALNORM_WINDOW | 50 | 07 | ±bins for local z-score |
| CHEAT_DENOISED_LAMBDA | 1.0 | 07 | Fused lasso / proxy strength |
| CHEAT_RESIDBG_RANK | 3 | 07 | SVD rank for background |
| CHEAT_SUPPORT_QUANTILE | 0.75 | 07 | Adaptive threshold quantile |
| BLOC_NEAR_FRAC | 0.2 | 08 | Near zone for squareness |
| BLOC_FAR_FRAC | 0.2 | 08 | Far zone for squareness |
| BLOC_EDGE_DEPTH | 3 | 08 | Bins for edge contrast |
| NN_OVERLAP_THRESH | 0.70 | 09 | Reciprocal overlap for matching |
| CONSENSUS_MIN_VARIANTS | 3 | 10 | Min variants for HIGH confidence |
| CONSENSUS_OVERLAP | 0.50 | 10 | Overlap for cross-variant matching |
| NN_SURVIVAL_RATIO | 0.5 | 02 | Signal retention threshold |
| SV_MAX_DIST_KB | 500 | 06 | Max SV breakpoint distance |
| WINDOW_SIZE_BP | 50000 | all | Bin-to-bp conversion |


## Known limitations

- **Matrix cheats are slow** on 9000×9000 matrices.  A2 (local contrast) and
  A3 (denoising) are O(N²×window²).  Use `--skip-cheats` for initial tests.
- **GHSL module** uses sim_mat fallback, not real GHSL v5.3 band assignments.
- **Denoising** uses iterated median filter as proxy.  Install `flsa` for real
  fused lasso: `install.packages("flsa")`.
- **ICA** requires `fastICA`: `install.packages("fastICA")`.
- **NN sweep tree** (Phase 5) only uses pre-loaded scales.  The full
  nn20-to-nn4000 sweep requires `run_triangle_multiscale.R` to compute
  sim_mats at each scale first.
- **No plotting yet.**  Add diagnostic plots after seeing real output.


## What to tell Claude next session

"Here's the scoring table for LG01 and the consensus table.
The known inversions at [positions] have these scores.
The known negatives at [positions] have these scores.
Help me: (1) set thresholds, (2) add diagnostic plots,
(3) wire in real GHSL band assignments, (4) add per-stair Q recomputation."
