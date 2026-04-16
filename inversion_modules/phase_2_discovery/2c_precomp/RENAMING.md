# RENAMING.md — inversion pipeline terminology migration

**Status:** incremental. The renames below were applied to the three phase 2/2c scripts first. Every other script in the codebase that references the old names needs the same substitutions applied one by one, using the grep recipes at the bottom of this doc to audit each file.

**Rule of thumb:** every rename in this file corresponds to a search-replace. There is no logic change anywhere — this is pure terminology cleanup.

## 1. Why the rename

The pipeline grew names organically during development: "flashlight" for an SV evidence lookup, "snake" for a sequential region-growing algorithm (borrowed from the Kass–Witkin–Terzopoulos 1988 active contours paper), "cores" for the seeded regions the algorithm produces, "cheats" for add-on validation tests. These are all informal nicknames, not published or scientific terms. The manuscript and the codebase use scientific vocabulary. The rename makes the codebase match.

## 2. Master rename table

### 2.1 Concept terms (primary renames)

| Old | New | What it actually is | Why this name |
|---|---|---|---|
| `flashlight` | `sv_prior` | Per-chromosome RDS packaging SV INV calls, het-DELs at breakpoints, BND-triangulated INVs, and the sample×inversion genotype table | It is an SV-derived informed prior consumed by downstream scripts |
| `snake` (algorithm) | `seeded region-growing`, or short: `seeded extension` | Greedy 1D seed-and-extend along a chromosome with damage-budget early termination | "Region growing" is standard image-segmentation / computational biology terminology. Important clarification: seeds come from **MDS z-score outliers**, not from the algorithm's output — this was not clear under the old name |
| `snake` (output vector) | `region` (or `seeded_region` in docs) | The contiguous stretch of windows grown from a seed | What it is |
| `core` | `seeded_region` (noun), or the column `scale_tier` when referring to 1S/1M/1L | The seeded extension output before merge | "Core" was a historical term meaning "this is the before-merge output"; that role is clearer as "seeded region". The 1S/1M/1L labels refer to parameter scale (small/medium/large), not biological families |
| `core family` / `family` (algorithm param set) | `scale_tier` | The three parameter sets run independently per chromosome | These are parameter scales, not biological families — "family" was misleading |
| `cheat` prefix | `test_NN` where NN is the two-digit cheat number | Add-on validation tests that confirm or contradict inversion identity | Each one performs a statistical or evidence test; "cheat" implied shortcut/hack |
| `snake1_*` file prefix | `seeded_regions_*` | Output files from the region-growing step | Matches new terminology |
| `snake_regions_multiscale/` dir | `regions/` (suggested, **NOT YET APPLIED** — see section 5) | Container dir for all region-level outputs | Shorter, scientific |

### 2.2 Script filename renames

| Old filename | New filename |
|---|---|
| `STEP_C00_build_flashlight_wired8_registry.R` | `STEP_C00_build_sv_prior.R` |
| `STEP_C01a_snake1_precompute_wired7_12358_v934.R` | `STEP_C01a_precompute.R` |
| `STEP_C01b_1_cores_wired_registry.R` | `STEP_C01b_1_seeded_regions.R` |

### 2.3 R identifiers (variables, functions, fields)

| Old | New | Scope |
|---|---|---|
| `flashlight` (list object in C00) | `sv_prior` | R object |
| `FLASH_DIR` | `SV_PRIOR_DIR` | env / config var |
| `FLASHLIGHT_LOADER` | `SV_PRIOR_LOADER` | env var for the loader path |
| `flashlight_dir` (CLI arg) | `sv_prior_dir` | C01a CLI |
| `--flashlight_dir` | `--sv_prior_dir` | C01a CLI flag |
| `snake_id` | `region_id` | column + field in C01b |
| `snake_id_start` | `region_id_start` | function parameter |
| `snake_phase` | `extension_phase` | column in decision log + regions table |
| `snake_span_mb` | `region_span_mb` | local variable |
| `core_family` | `scale_tier` | column + field |
| `core_1S` / `core_1M` / `core_1L` | `region_1S` / `region_1M` / `region_1L` | window-states columns |
| `cores_1S` / `cores_1M` / `cores_1L` | `n_regions_1S` / `n_regions_1M` / `n_regions_1L` | summary columns |
| `cores_total` | `n_regions_total` | summary column |
| `all_cores` | `all_regions` | R list |
| `fam_cores` | `fam_regions` | R list |
| `run_core_family` | `run_extension_scale` | function name |
| Loop variable `snake` (growing index vector) | `region` | local in `run_extension_scale` |
| Phase string `"core"` in decision_log | `"extend"` | literal value |
| `cheat2_*` | `test02_*` | C00 internals |
| `cheat5_*` | `test05_*` | C01a Fst-scan columns |
| `cheat26_*` | `test26_*` | C01b kin-pruning columns |
| `n_cheat2_confirmed` | `n_test02_confirmed` | summary column |

### 2.4 Log prefixes

| Old | New |
|---|---|
| `[flashlight]` | `[sv_prior]` |
| `[cheat8]` | `[sv_prior]` (was the BND sub-module of the flashlight build) |
| `[cheat26]` | `[test_26]` |
| `[cores]` | `[seeded_regions]` |

### 2.5 Output file basenames

| Old basename | New basename | Produced by |
|---|---|---|
| `sv_flashlight_<chr>.rds` | `sv_prior_<chr>.rds` | C00 |
| `flashlight_summary.tsv` | `sv_prior_summary.tsv` | C00 |
| `snake_inv_likeness.tsv.gz` | `window_inv_likeness.tsv.gz` | C01a |
| `snake1_cores_<chr>.rds` | `seeded_regions_<chr>.rds` | C01b |
| `snake1_core_windows_<chr>.tsv.gz` | `seeded_regions_windows_<chr>.tsv.gz` | C01b |
| `snake1_core_regions_<chr>.tsv.gz` | `seeded_regions_summary_<chr>.tsv.gz` | C01b |
| `snake1_decision_log_<chr>.tsv.gz` | `seeded_regions_decision_log_<chr>.tsv.gz` | C01b |
| `snake1_window_states_<chr>.tsv.gz` | `seeded_regions_window_states_<chr>.tsv.gz` | C01b |
| `snake1_summary_<chr>.tsv` | `seeded_regions_summary_per_chr_<chr>.tsv` | C01b |

## 3. Cheat → test_NN file rename map

Each `cheatN_*.R` file in `cheats/` becomes `testNN_*.R`. Two-digit zero-padded for consistent sort order.

| Old filename | New filename | What it measures |
|---|---|---|
| `cheat4_boundary_sharpness.R` | `test04_boundary_sharpness.R` | Fst / dXY / Hobs step function at inversion boundaries |
| `cheat5_qaxis_fst_scan.R` | `test05_qaxis_fst_scan.R` | Q-axis Fst scan (Engine B) |
| `cheat6_ancestry_jackknife.R` | `test06_ancestry_jackknife.R` | Leave-one-family-out Fst robustness |
| `cheat7_beta_adaptive.R` | `test07_beta_adaptive_threshold.R` | Beta(α, β) adaptive seed threshold |
| `cheat8_bnd_triangulation.R` | `test08_bnd_triangulation.R` | BND paired-breakpoint triangulation |
| `cheat9_indel_slope.R` | `test09_cumulative_indel_slope.R` | Cumulative-indel-slope heterokaryotype detector |
| `cheat10_depth_boundary.R` | `test10_depth_boundary.R` | Read-depth anomaly at sim_mat boundaries |
| `cheat11_clip_boundary.R` | `test11_clip_boundary.R` | Soft/hard-clipped read pileup at boundaries |
| `cheat12_theta_het_prior.R` | `test12_theta_het_prior.R` | Per-sample θ_P heterozygosity prior |
| `cheat13_breakpoint_proximity_weight.R` | `test13_breakpoint_proximity_weight.R` | Breakpoint-proximal weighting |
| `cheat14_repeat_architecture.R` | `test14_repeat_architecture.R` | Repeat / SD architecture at breakpoints |
| `cheat14_self_align.sh` | `test14_self_align.sh` | minimap2 self-alignment (test_14 helper) |
| `cheat15_recurrence_test.R` | `test15_recurrence_vs_single_origin.R` | Recurrence vs single-origin test |
| `cheat16_multi_arrangement.R` | `test16_multi_arrangement.R` | Multi-arrangement (k > 3) detector |
| `cheat17_fossil_breakpoints.R` | `test17_fossil_breakpoints.R` | Fossil-breakpoint detection |
| `cheat18_mendelian_test.R` | `test18_mendelian_transmission.R` | Mendelian transmission test |
| `cheat19_diversity_gradient.R` | `test19_diversity_gradient.R` | Within-inversion diversity gradient (U-shape) |
| `cheat20_junction_classifier.R` | `test20_junction_microhomology.R` | Breakpoint junction microhomology classifier |
| `cheat20_extract_junctions.py` | `test20_extract_junctions.py` | Junction extraction (test_20 helper) |
| `cheat21_te_context.R` | `test21_te_enrichment.R` | TE enrichment / depletion at breakpoints |
| `cheat22_secondary_structure.R` | `test22_secondary_structure.R` | Secondary-structure scan within inversion classes |
| `cheat23_extended_suppression.R` | `test23_extended_boundary_suppression.R` | Extended boundary recombination suppression |
| `cheat24_recombinant_prior.R` | `test24_recombinant_prior.R` | Position-aware recombinant prior |
| `cheat25_block_viability.R` | `test25_block_viability.R` | Block-viability 4-test battery |
| `cheat26_kin_pruned_retention.R` | `test26_kin_pruned_retention.R` | Kin-pruned signal retention |

## 4. Per-file migration checklist

Apply this procedure to every script that references any old term.

### Step 1. Filename

```bash
git mv <old>.R <new>.R      # for each file being renamed
```

### Step 2. Internal search-replace

Run this sed sequence against any single file. It is safe to run more than once.

```bash
FILE="$1"

# Concept term -> identifier renames
sed -i \
  -e 's/\bflashlight\b/sv_prior/g' \
  -e 's/\bFLASH_DIR\b/SV_PRIOR_DIR/g' \
  -e 's/\bFLASHLIGHT_LOADER\b/SV_PRIOR_LOADER/g' \
  -e 's/\bflashlight_dir\b/sv_prior_dir/g' \
  -e 's/--flashlight_dir/--sv_prior_dir/g' \
  -e 's/\bsnake_id\b/region_id/g' \
  -e 's/\bsnake_phase\b/extension_phase/g' \
  -e 's/\bsnake_id_start\b/region_id_start/g' \
  -e 's/\bsnake_span_mb\b/region_span_mb/g' \
  -e 's/\bcore_family\b/scale_tier/g' \
  -e 's/\bcore_1S\b/region_1S/g' \
  -e 's/\bcore_1M\b/region_1M/g' \
  -e 's/\bcore_1L\b/region_1L/g' \
  -e 's/\bcores_1S\b/n_regions_1S/g' \
  -e 's/\bcores_1M\b/n_regions_1M/g' \
  -e 's/\bcores_1L\b/n_regions_1L/g' \
  -e 's/\bcores_total\b/n_regions_total/g' \
  -e 's/\ball_cores\b/all_regions/g' \
  -e 's/\bfam_cores\b/fam_regions/g' \
  -e 's/\brun_core_family\b/run_extension_scale/g' \
  "$FILE"

# cheatN_ -> testNN_  (per-file renames — run only for the cheats you have touched)
sed -i \
  -e 's/\bcheat2_verification\b/test02_verification/g' \
  -e 's/\bcheat2_dt\b/test02_dt/g' \
  -e 's/\bcheat2_concordant\b/test02_concordant/g' \
  -e 's/\bcheat2_works\b/test02_works/g' \
  -e 's/\bcheat5_fst_pc1\b/test05_fst_pc1/g' \
  -e 's/\bcheat5_fst_q_best\b/test05_fst_q_best/g' \
  -e 's/\bcheat5_fst_q_best_k\b/test05_fst_q_best_k/g' \
  -e 's/\bcheat5_family_fst_ratio\b/test05_family_fst_ratio/g' \
  -e 's/\bcheat26_status\b/test26_status/g' \
  -e 's/\bcheat26_retention\b/test26_retention/g' \
  "$FILE"

# Output file basenames
sed -i \
  -e 's|\bsv_flashlight_|sv_prior_|g' \
  -e 's|flashlight_summary\.tsv|sv_prior_summary.tsv|g' \
  -e 's|snake_inv_likeness\.tsv\.gz|window_inv_likeness.tsv.gz|g' \
  -e 's|snake1_cores_|seeded_regions_|g' \
  -e 's|snake1_core_windows_|seeded_regions_windows_|g' \
  -e 's|snake1_core_regions_|seeded_regions_summary_|g' \
  -e 's|snake1_decision_log_|seeded_regions_decision_log_|g' \
  -e 's|snake1_window_states_|seeded_regions_window_states_|g' \
  -e 's|snake1_summary_|seeded_regions_summary_per_chr_|g' \
  "$FILE"

# Log prefixes
sed -i \
  -e 's/\[flashlight\]/[sv_prior]/g' \
  -e 's/\[cheat8\]/[sv_prior]/g' \
  -e 's/\[cheat26\]/[test_26]/g' \
  -e 's/\[cores\]/[seeded_regions]/g' \
  "$FILE"
```

### Step 3. Loop variable `snake` → `region`

This one is risky to sed blindly because `snake` appears in comments too and you want the replacement in comments as well — but it can collide with user-facing documentation words like "sneak" or "snakemake". Do this one carefully per-file:

```bash
# Check first what matches
grep -nE "\bsnake\b|\bsnakes\b" <FILE>

# Apply only if the matches are all the loop variable and comment references:
sed -i -e 's/\bsnake\b/region/g' -e 's/\bsnakes\b/regions/g' <FILE>
```

Never run this on files that mention `snakemake` — if such files exist, do the replacement by hand.

### Step 4. Log decision phase literal

```bash
sed -i 's|log_decision(\([^,]*\), \([^,]*\), "core",|log_decision(\1, \2, "extend",|g' <FILE>
```

### Step 5. Header metadata block

Replace by hand using the template below.

### Step 6. Audit

After applying, confirm nothing stale remains:

```bash
grep -nE "\bflashlight|\bsnake|\bcore[^_a-z]|Cheat|CHEAT|\bcheat[0-9]" <FILE>
```

Allowed residues: `Formerly:` provenance lines, explicit "was cheatN" breadcrumb comments, and backward-compat read paths like `sv_flashlight_` in file-existence fallbacks. Everything else is stale.

## 5. Things NOT yet renamed (scoped out)

These require coordinated changes across many files and are tracked here for a future pass:

### 5.1 The `snake_regions_multiscale/` directory

`00_inversion_config.sh` defines:

```
SNAKE1_DIR="${MDS_DIR}/snake_regions_multiscale"
LANDSCAPE_DIR="${MDS_DIR}/snake_regions_multiscale/landscape"
TRIANGLE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles"
SCORING_DIR="${MDS_DIR}/snake_regions_multiscale/scoring"
DECOMPOSITION_DIR="${MDS_DIR}/snake_regions_multiscale/decomposition"
REGIME_DIR="${MDS_DIR}/snake_regions_multiscale/regime_engine"
PRECOMP_DIR="${SNAKE1_DIR}/precomp"
SIM_MATS_DIR="${PRECOMP_DIR}/sim_mats"
TRIANGLE_MULTISCALE_DIR="${MDS_DIR}/snake_regions_multiscale/triangles_multiscale"
SNAKE2_DIR / SNAKE3_DIR / SNAKE4_DIR
SNAKE1_DIAG_DIR / SNAKE_PRIMARY_DIR / SNAKE_SECONDARY_DIR / SNAKE_DIAG_DIR
SNAKE_CAND_FILE
```

Proposed replacements when you do the config pass:

```
REGIONS_DIR="${MDS_DIR}/regions"                                 # was SNAKE1_DIR
REGIONS_DIAG_DIR="${MDS_DIR}/regions_diagnostics"                # was SNAKE1_DIAG_DIR
RG_COMMUNITY_DIR / RG_GHSL_DIR / RG_THREEBAND_DIR                # was SNAKE2/3/4_DIR
LANDSCAPE_DIR="${REGIONS_DIR}/landscape"
TRIANGLE_DIR="${REGIONS_DIR}/triangles"
SCORING_DIR="${REGIONS_DIR}/scoring"
DECOMPOSITION_DIR="${REGIONS_DIR}/decomposition"
REGIME_DIR="${REGIONS_DIR}/regime_engine"
PRECOMP_DIR="${REGIONS_DIR}/precomp"
SIM_MATS_DIR="${PRECOMP_DIR}/sim_mats"
TRIANGLE_MULTISCALE_DIR="${REGIONS_DIR}/triangles_multiscale"
REGIONS_PRIMARY_CAND_FILE                                        # was SNAKE_CAND_FILE
SV_PRIOR_DIR="${REGIONS_DIR}/sv_prior"                           # pattern match
```

This touches every downstream module (5B / 5C / 5D / 5E / MODULE_5A2_breakpoint_validation) plus any already-produced output directories on LANTA. Coordinate with a directory rename and a set of symlinks for backward compatibility.

### 5.2 The snake1 / snake2 / snake3 / snake4 algorithm variants

These are four feature sets of the same region-growing algorithm. Current internal distinctions:

| Old | What it adds | Proposed |
|---|---|---|
| snake1 | Multiscale similarity (base algorithm) | `region_grow_multiscale` or short `rg_multiscale` |
| snake2 | Community detection on sample graph | `region_grow_community` / `rg_community` |
| snake3 | GHSL haplotype contrast | `region_grow_ghsl` / `rg_ghsl` |
| snake4 | Three-band decomposition | `region_grow_threeband` / `rg_threeband` |

When doing that rename, keep the numbering in the short form (snake1 → rg1 or rg_multiscale) so cross-references stay stable.

### 5.3 `N_CORES` (CPU count)

This one is **not** renamed — it's the OS-level CPU count for `parallel::mclapply`, which is genuinely "cores" in the correct CS sense.

## 6. Backward-compatibility read paths

C01a reads the SV prior with this fallback ladder:

```r
fl_file <- file.path(sv_prior_dir, paste0("sv_prior_", chr, ".rds"))             # new
if (!file.exists(fl_file)) {
  fl_file <- file.path(sv_prior_dir, paste0("sv_flashlight_", chr, ".rds"))       # old
}
if (!file.exists(fl_file)) {
  fl_file <- file.path(sv_prior_dir, paste0(chr, "_flashlight.rds"))              # oldest
}
```

This is intentional — it lets C01a read sv_prior RDS files produced by pre-rename C00 without re-running C00. Once every producer script has been re-run, the fallbacks can be dropped.

## 7. Provenance convention

Every renamed script carries a `Formerly:` line in its header so grep against old commits / old docs still lands on the right file:

```
# Formerly: STEP_C00_build_flashlight_wired8_registry.R
```

Do this for every script you rename.

## 8. Audit commands

**Find every stale reference across the codebase.** Run from the repository root:

```bash
grep -rn --include='*.R' --include='*.py' --include='*.sh' --include='*.md' \
  -E "\bflashlight|\bFLASH_DIR|\bsnake(_id|_phase|_span|[0-9])|\bcore(_family|s_total|s_1[SML])|\bcheat[0-9]" \
  .
```

**Count stale references per file** so you can prioritise:

```bash
grep -rl --include='*.R' --include='*.py' --include='*.sh' \
  -E "\bflashlight|\bsnake(_id|_phase|_span|[0-9])|\bcore_family|\bcheat[0-9]" . \
  | xargs -I {} sh -c 'echo "$(grep -cE "\bflashlight|\bsnake(_id|_phase|_span|[0-9])|\bcore_family|\bcheat[0-9]" "{}") {}"' \
  | sort -rn
```

**Confirm a single file is fully migrated:**

```bash
grep -nE "\bflashlight|\bsnake(_id|_phase|_span|[0-9])|\bcore(_family|s_total|s_1[SML])|\bcheat[0-9]" <FILE>
# Allowed residues: `Formerly:`, `was cheatN` breadcrumbs, backward-compat read paths
```

## 9. Files already migrated (phase 2 / 2c and 2d)

**Phase 2 / 2c:**
- `phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R`
- `phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`
- `phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R`

**Phase 2 / 2d (was `inv_detect_v9.3/`):**
- All 22 scripts + 3 launchers, fully renamed. See full map in section 12 below.

## 12. Phase 2 / 2d candidate detection rename map

The `inv_detect_v9.3/` drop has been migrated into
`phase_2_discovery/2d_candidate_detection/`. It replaces only the old
`STEP_C01c_triangle_regimes.R`; everything downstream from C01d is
unchanged, because the bridge script writes the legacy
`triangle_intervals.tsv.gz` format.

### Filenames

| Old | New |
|---|---|
| `inv_detect_v9.3/` | `phase_2_discovery/2d_candidate_detection/` |
| `00_run_all.R` | `run_all.R` |
| `01_staircase_detector.R` | `STEP_D01_staircase_boundaries.R` |
| `02_evidence_nn_persistence.R` | `STEP_D02_nn_persistence.R` |
| `03_evidence_flatness.R` | `STEP_D03_flatness.R` |
| `04_evidence_interior_cv.R` | `STEP_D04_interior_cv.R` |
| `05_evidence_ghsl.R` | `STEP_D05_ghsl_stability.R` |
| `06_evidence_sv_overlap.R` | `STEP_D06_sv_overlap.R` |
| `07_matrix_cheats.R` | `STEP_D07_matrix_transforms.R` |
| `08_bloc_scoring.R` | `STEP_D08_block_scoring.R` |
| `08_evidence_local_pca.R` | `STEP_D08b_local_pca.R` |
| `09_nn_sweep_tree.R` | `STEP_D09_nn_sweep_tree.R` |
| `10_consensus.R` | `STEP_D10_variant_consensus.R` |
| `11a_run_ngsrelate_perchr.sh` | `STEP_D11a_ngsrelate_perchr.sh` |
| `11b_build_perchr_pruning_table.R` | `STEP_D11b_perchr_pruning_table.R` |
| `12_bridge_to_codebase.R` | `STEP_D12_bridge_to_C01d.R` |
| `13_plot_annotated_simmat.R` | `STEP_D13_plot_annotated_simmat.R` |
| `14_landscape_classifier.R` | `STEP_D14_landscape_classifier.R` |
| `15_plot_zoomed_systems.R` | `STEP_D15_plot_zoomed_regions.R` |
| `16_test_staircase.R` | `tests/test_staircase.R` |
| `17_plot_marginal_tracks.R` | `STEP_D17_plot_marginal_tracks.R` |
| `STEP_C01n_blockwise_peeling_diagnostic.R` | `STEP_D09n_peeling_diagnostic.R` |
| `slurm_run.sh` | `LAUNCH_run_all.slurm` |
| `slurm_test_LG01.sh` | `LAUNCH_test_LG01.slurm` |
| `slurm_peel.sh` | `LAUNCH_peel.slurm` |

### Function renames

Run against every R file to migrate references:

```bash
sed -i \
  -e 's/\bcheat_distance_correction\b/mtx_distance_correction/g' \
  -e 's/\bcheat_local_contrast\b/mtx_local_contrast/g' \
  -e 's/\bcheat_denoise_proxy\b/mtx_denoise_proxy/g' \
  -e 's/\bcheat_denoise_fused_lasso\b/mtx_denoise_fused_lasso/g' \
  -e 's/\bcheat_background_residual\b/mtx_background_residual/g' \
  -e 's/\bcheat_background_smooth\b/mtx_background_smooth/g' \
  -e 's/\bcheat_edge_map\b/mtx_edge_map/g' \
  -e 's/\bcheat_support_map\b/mtx_support_map/g' \
  -e 's/\bcompute_bloc_metrics\b/compute_block_metrics/g' \
  "$FILE"
```

### Config variable renames

```bash
sed -i \
  -e 's/\bCHEAT_DISTCORR_ENABLED\b/MTX_DISTCORR_ENABLED/g' \
  -e 's/\bCHEAT_LOCALNORM_ENABLED\b/MTX_LOCALNORM_ENABLED/g' \
  -e 's/\bCHEAT_LOCALNORM_WINDOW\b/MTX_LOCALNORM_WINDOW/g' \
  -e 's/\bCHEAT_DENOISED_ENABLED\b/MTX_DENOISED_ENABLED/g' \
  -e 's/\bCHEAT_DENOISED_LAMBDA\b/MTX_DENOISED_LAMBDA/g' \
  -e 's/\bCHEAT_RESIDBG_ENABLED\b/MTX_RESIDBG_ENABLED/g' \
  -e 's/\bCHEAT_RESIDBG_RANK\b/MTX_RESIDBG_RANK/g' \
  -e 's/\bCHEAT_RESIDBG_KERNEL\b/MTX_RESIDBG_KERNEL/g' \
  -e 's/\bCHEAT_EDGE_ENABLED\b/MTX_EDGE_ENABLED/g' \
  -e 's/\bCHEAT_SUPPORT_ENABLED\b/MTX_SUPPORT_ENABLED/g' \
  -e 's/\bCHEAT_SUPPORT_QUANTILE\b/MTX_SUPPORT_QUANTILE/g' \
  -e 's/\bBLOC_NEAR_FRAC\b/BLOCK_NEAR_FRAC/g' \
  -e 's/\bBLOC_FAR_FRAC\b/BLOCK_FAR_FRAC/g' \
  -e 's/\bBLOC_EDGE_DEPTH\b/BLOCK_EDGE_DEPTH/g' \
  "$FILE"
```

### CLI flag + output basename renames

```bash
sed -i \
  -e 's/--skip-cheats/--skip-transforms/g' \
  -e 's/\bskip_cheats\b/skip_transforms/g' \
  -e 's/\bSKIP_CHEATS\b/SKIP_TRANSFORMS/g' \
  -e 's/\bbloc_scores_/block_scores_/g' \
  "$FILE"
```

### Preserved legacy interfaces (do NOT rename)

- **Output filename `triangle_intervals.tsv.gz`** from
  `STEP_D12_bridge_to_C01d.R`. C01d reads this exact name. Changing it
  would require changing C01d, which is out of scope. Flagged as a
  future-cleanup item.
- **Geometric "triangle" / "upper triangle" / "lower triangle"** in
  `STEP_D17_plot_marginal_tracks.R`. These refer to the geometric
  upper/lower half of a symmetric matrix. Correct math, keep as-is.
- **`STAIR_*`, `NN_*`, `CONSENSUS_*`, `LANDSCAPE_*`, `PEEL_*`** config
  prefixes. These describe the parameter's function, not a historical
  nickname.

## 13. Files now known to still contain old terms

Original section 10 list, updated:

- `00_inversion_config.sh` (path variables — see section 5.1)
- `phase_2_discovery/2c_precomp/patches/*.R` (all 6 flashlight patches)
- `phase_2_discovery/2d_cores/STEP_C01b_2_merge.R` and the other snake*
  files (these may now be superseded by `2d_candidate_detection/` but
  the seed-based merge step still exists — confirm whether C01b_2 is
  kept as the merge step for the seed-based track)
- `phase_2_discovery/2d_cores/STEP_C03_snake2_community.R`
- `phase_2_discovery/2d_cores/STEP_C05_snake4_threeband.R`
- `phase_2_discovery/2d_cores/STEP_C06_consensus.R`
- `phase_2_discovery/2d_cores/STEP_C07_hmm_regime.R`
- `phase_2_discovery/2d_cores/STEP_C08_rescue_pass2.R`
- `phase_2_discovery/2e_ghsl/STEP_C04_snake3_ghsl_v5.R`
- `phase_4_catalog/` — every script that writes `cheatN_*` columns
- `phase_4_catalog/cheats/*` — all 23 cheat files (see section 3 for
  the rename map)
- `MODULE_5A2_breakpoint_validation/` — minor references
- `MODULE_5B / 5C / 5D / 5E` — consume `SNAKE*_DIR` paths, need path
  rename

Prioritise by the count from the audit command in section 8.

**Important clarification:** the folder called `phase_2_discovery/2d_cores/`
in your earlier tree note is separate from the new
`phase_2_discovery/2d_candidate_detection/`. The former is the
seed-based track's merge/community/threeband/consensus/HMM/rescue
scripts (the "snakes" family). The latter is the matrix-based track
from `inv_detect_v9.3/`. Both live under phase 2 discovery, both
produce candidates, both feed C01d.

Consider renaming the seed-based folder to something like
`phase_2_discovery/2d_seeded_extension/` or
`phase_2_discovery/2d_region_growing/` to avoid confusion with the
new `2d_candidate_detection/`. Or keep the two as explicit parallel
sub-phases:

```
phase_2_discovery/
├── 2c_precomp/                 (both tracks read from here)
├── 2d_seeded_regions/          (was 2d_cores — seed-based track)
└── 2d_candidate_detection/     (bumped up from 2d to avoid collision)
```

This is a folder-layout decision for you to make; I've used
`2d_candidate_detection/` in this drop on the assumption that the
seed-based track's folder will be renamed / bumped separately.

## 11. Contract — what does not change

- Logic of any algorithm (seed selection, extension, damage accumulation, decision thresholds, Fst computation, etc.) is identical
- Numerical outputs are identical
- CLI argument structure and downstream RDS schemas are compatible (old reads old, new reads new, and C01a reads both via the fallback ladder)
- Output directory structure is the same until section 5.1 is done
