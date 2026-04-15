# Pipeline Wiring v9.0 — Ind→CGA Root Fix + Cross-Module Registry

## What this delivers

This package rewires the entire pipeline to fix the Ind0/Ind1 vs CGA009/CGA010 sample name mismatch at the root, and connects the sample registry + unified ancestry module across all scripts.

## Architecture

```
pipeline_bridge.sh  ←── bash entry point (SLURM launchers source this)
    ↓ exports env vars
load_bridge.R       ←── R entry point (all R scripts source this)
    ↓ provides: smap, reg, get_Q(), get_region_stats()
    ├── sample_map.R      (Ind↔CGA translation, unchanged)
    ├── sample_registry.R (group management, unchanged)
    ├── instant_q.R       (fixed-F EM engine wrapper)
    └── region_stats_dispatcher.R (unified stats)
```

## The Root Fix

**Problem:** BEAGLE files use `Ind0, Ind1, ...`. Precomp RDS files inherit these names. Every script crossing the BEAGLE↔Clair3 boundary needs mapping.

**Fix:** `STEP_C01a_snake1_precompute.R` now renames columns at BEAGLE load time:

```r
# In STEP_C01a (the root fix):
real_names <- smap$to_real_vec(raw_sample_names)  # Ind0→CGA009
setnames(bgl, old_cols, new_cols)                  # rename in-place
```

After this, precomp RDS has `PC_1_CGA009` instead of `PC_1_Ind0`. Every downstream script just works.

**Cost:** One precomp rerun (~2 hours for 28 chromosomes as SLURM array).

## Files delivered

### Core infrastructure (copy to `inversion_codebase_v8.5/utils/`)

| File | Purpose |
|------|---------|
| `utils/load_bridge.R` | Universal R loader — sources smap, registry, instant_q, dispatcher |
| `utils/pipeline_bridge.sh` | Universal bash loader — exports all env vars for R scripts |

### Rewired inversion scripts (copy to `inversion_codebase_v8.5/`)

| File | Changes |
|------|---------|
| `STEP_C01a_snake1_precompute.R` | **ROOT FIX**: Ind→CGA at load. Instant Q integration. Registry `all_226`. |
| `STEP_C04_snake3_ghsl_v3.R` | Removed sample_map hack. Auto-registers GHSL genotype groups. |
| `STEP_C01d_candidate_scoring.R` | Direct CGA names from precomp. Uses `get_region_stats()` for Fst/theta_pi. |
| `STEP_C01i_multi_inversion_decomposition.R` | Native CGA. Registers `inv_{cid}_HOM_REF/HET/HOM_INV` groups. |
| `STEP_C01f_hypothesis_tests.R` | CGA in all output labels. Q divergence test via instant_q. |

### Updated MODULE_2B (replace `MODULE_2B_Structure/`)

| File | Changes |
|------|---------|
| `00_module2b_config.sh` | Added REGISTRY_DIR, LOAD_BRIDGE, ANCESTRY_CONFIG, SAMPLES_IND |
| `scripts/STEP_A04_best_seed_by_K.R` | Auto-registers `ancestry_K8_Q1..Q8` groups. Adds `Q_major_K8` to master. |
| `launchers/LAUNCH_module2b_figures.sh` | Sources pipeline_bridge.sh. Passes --bridge to A04. |

### Unchanged files (included for completeness)

All B-series viz scripts, slurm workers, compute launcher, A01/A05 — unchanged.
`sample_map.R` and `sample_registry.R` — unchanged (they were already correct).

## Deployment

```bash
BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04

# 1. Copy infrastructure
cp utils/load_bridge.R      ${BASE}/inversion_codebase_v8.5/utils/
cp utils/pipeline_bridge.sh ${BASE}/inversion_codebase_v8.5/utils/

# 2. Copy rewired inversion scripts (backup originals first!)
cd ${BASE}/inversion_codebase_v8.5
for f in STEP_C01a STEP_C04 STEP_C01d STEP_C01i STEP_C01f; do
  cp MODULE_5A2_Discovery_Core/${f}*.R MODULE_5A2_Discovery_Core/${f}.bak.R 2>/dev/null
done
cp inversion_rewired/STEP_C01a_snake1_precompute.R  MODULE_5A2_Discovery_Core/
cp inversion_rewired/STEP_C04_snake3_ghsl_v3.R       MODULE_5A2_Discovery_Core/
cp inversion_rewired/STEP_C01d_candidate_scoring.R   MODULE_5A2_Discovery_Core/
cp inversion_rewired/STEP_C01i_multi_inversion_decomposition.R MODULE_5A2_Discovery_Core/
cp inversion_rewired/STEP_C01f_hypothesis_tests.R    MODULE_5A2_Discovery_Core/

# 3. Replace MODULE_2B
cp -r MODULE_2B_Structure/ ${BASE}/MODULE_2B_Structure/

# 4. Initialize registry (first time only)
source ${BASE}/inversion_codebase_v8.5/utils/pipeline_bridge.sh

# 5. Re-run precomp (THE KEY STEP — fixes everything downstream)
#    Submit as SLURM array:
sbatch --array=1-28 --account=lt200308 --partition=compute \
  --time=4:00:00 --mem=64G --cpus-per-task=4 \
  --wrap='
    source /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/utils/pipeline_bridge.sh
    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BASE}/chrom_list.txt)
    Rscript ${BASE}/inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/STEP_C01a_snake1_precompute.R --chr $CHR
  '
```

## Registry groups auto-registered

| Group ID | When | Source |
|----------|------|--------|
| `all_226` | Precomp (LG01) | STEP_C01a |
| `unrelated_81` | Precomp (LG01) or pipeline_bridge.sh | NAToRA pruned list |
| `ancestry_K8_Q1` .. `ancestry_K8_Q8` | STEP_A04 | NGSadmix best-K clusters |
| `inv_{cid}_HOM_REF` | STEP_C01i | Per-inversion genotype decomposition |
| `inv_{cid}_HET` | STEP_C01i | Per-inversion genotype decomposition |
| `inv_{cid}_HOM_INV` | STEP_C01i | Per-inversion genotype decomposition |
| `ghsl_{cid}_HOM_REF/HET/HOM_INV` | STEP_C04 | Snake3 GHSL classification |

Groups are **NOT** created for:
- Per-window k=3 bands (transient intermediate results)
- Intermediate clustering results
- Module 6 founder packs (deferred to Module 6 rewiring)

## What doesn't need fixing

These scripts read from a single source chain (precomp only) and will automatically get CGA names after the precomp rerun:

- `STEP_C01b_1_cores.R` — reads precomp only
- `STEP_C01b_2_merge.R` — reads precomp + landscape (both fixed)
- `STEP_C01m_distance_concordance.R` — precomp only
- `C01a_diag_*.R` — precomp only (plots now show CGA labels)
- `STEP_C01b_plot_heatmaps.R` — core labels now CGA
- `STEP_C01l_local_structure_segments.R` — output columns now CGA

## Config coexistence

All three configs share `BASE` but use module-prefixed variables:

| Config | Prefix | Example vars |
|--------|--------|-------------|
| `00_inversion_config.sh` | `INVDIR`, `SNAKE*_DIR` | `inv_log`, `inv_die` |
| `00_ancestry_config.sh` | `INSTANT_Q_*`, `POPSTATS_*` | `anc_log`, `anc_die` |
| `00_module2b_config.sh` | `MODULE2B_*` | `m2b_log`, `m2b_die` |

`pipeline_bridge.sh` sources all three safely. `load_bridge.R` parses all three via env vars.
