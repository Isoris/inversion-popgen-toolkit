# MODULE_2B — Population Structure & Visualization

Extended ancestry inference with K=2–20, per-chromosome structure, evalAdmix model diagnostics, K-hierarchy merge tree, and publication-quality figure suite. Consumes MODULE_2A SNP panels and produces the definitive ancestry assignments, kinship figures, and sample registry entries used by all downstream modules.

## Pipeline

```
MODULE_2A outputs (BEAGLE GLs, thin panels, pruned sample list)
  │
  ├─ A. Compute ──────────────────────────────────────────────────────
  │   A01  generate task list (K×seed×thin×set combos)
  │   A02  NGSadmix array (K=2..20 × 5 seeds × panels × sample sets)
  │   A03  evalAdmix array (same task list, dependency on A02)
  │   A04  best-seed-by-K selection + palette + registry export
  │   A05  K-hierarchy merge tree (founder lineage structure)
  │
  └─ B. Visualization ────────────────────────────────────────────────
      B01  admixture faceted composite (K=2..20 barplots + evalAdmix)
      B02  kinship network (ggraph, ancestry-colored)
      B03  kinship heatmap (226×226, PCAngsd-ordered)
      B04  PCA + kinship overlay (PC1 vs PC2 with edges)
      B05  K-hierarchy Sankey (alluvial founder merge tree)
      B06  relatedness summary (4-panel pruning composite)
      B07  evalAdmix diagnostic (per-K residual heatmaps)
  │
  └─→ Ancestry Q tables + palette + figures + sample registry
        consumed by MODULE_3, MODULE_5A, MODULE_6, manuscript
```

## Directory layout

```
MODULE_2B_Structure/
  00_module2b_config.sh                  ← central config (source this everywhere)
  README.md
  docs/
    MODULE_2B_methods.md                 ← manuscript-ready methods prose
  scripts/
    STEP_A01_generate_tasks.sh           ← compute: task list generator
    STEP_A04_best_seed_by_K.R            ← compute: best-seed selection + exports
    STEP_A05_k_hierarchy_tree.R          ← compute: founder merge tree
    STEP_B01_admixture_faceted.R         ← viz: K=2..20 barplots + evalAdmix strips
    STEP_B02_kinship_network.R           ← viz: ggraph kinship network
    STEP_B03_kinship_heatmap.R           ← viz: 226×226 annotated heatmap
    STEP_B04_pca_kinship_overlay.R       ← viz: PCA + kinship edges
    STEP_B05_k_hierarchy_sankey.R        ← viz: alluvial founder merge tree
    STEP_B06_relatedness_summary.R       ← viz: 4-panel pruning composite
    STEP_B07_evaladmix_diagnostic.R      ← viz: per-K residual diagnostics
  slurm/
    SLURM_A02_ngsadmix_worker.sh         ← SLURM array: one NGSadmix per task
    SLURM_A03_evaladmix_worker.sh        ← SLURM array: one evalAdmix per task
  launchers/
    LAUNCH_module2b_compute.sh           ← submits A01→A02→A03 with dependencies
    LAUNCH_module2b_figures.sh           ← runs A04→A05→B01..B07 for each panel
```

## Naming convention

Matches inversion codebase v8.5 and MODULE_2A:

- Config: `00_module2b_config.sh` (like `00_inversion_config.sh`)
- Scripts: `STEP_{PHASE}{NN}_{description}.{R|sh}`
  - Phase A = compute (A01 tasks, A02 NGSadmix, A03 evalAdmix, A04 best-seed, A05 merge tree)
  - Phase B = visualization (B01–B07)
- SLURM: `SLURM_{PHASE}{NN}_{description}.sh`
- Launchers: `LAUNCH_module2b_{purpose}.sh`
- Helper functions: `m2b_log`, `m2b_die`, `m2b_check_file`, `m2b_init_dirs`

## Usage

```bash
# ── Compute (submit to SLURM) ──────────────────────────────────────
bash launchers/LAUNCH_module2b_compute.sh               # whole-genome only
bash launchers/LAUNCH_module2b_compute.sh --per-chr      # + per-chromosome
bash launchers/LAUNCH_module2b_compute.sh --dry-run      # preview without submitting

# ── Figures (after all jobs finish) ─────────────────────────────────
bash launchers/LAUNCH_module2b_figures.sh --best-k 8
```

## Task matrix

| Scope | K values | Seeds | Thin panels | Sample sets | Tasks |
|-------|----------|-------|-------------|-------------|-------|
| Whole-genome | 2–20 | 5 | 500, 1000 | all, pruned | 285 |
| + per-chromosome | 2–20 | 5 | 500 | all | +2,660 |
| **Total with --per-chr** | | | | | **2,945** |

## Key design decisions

**Extended K range.** K=2–20 with 5 seeds (vs MODULE_2A's K=2–12 with 3 seeds) provides comprehensive exploration of population substructure. The hatchery founder population may contain subtle sub-lineages only visible at higher K.

**Two-phase launcher architecture.** Compute (SLURM arrays with dependency chains) and figures (sequential R scripts) are separate launchers. This avoids blocking SLURM resources during figure generation and allows re-running figures without re-computing.

**K-hierarchy merge tree.** For consecutive K values (K, K−1), ancestry components are matched by maximum Q-overlap to determine which components merge. This reveals the nested founder lineage structure in the hatchery population.

**Cross-module wiring (v9.0).** Config exports `REGISTRY_DIR`, `LOAD_BRIDGE`, `ANCESTRY_CONFIG`, and `INVERSION_CONFIG` paths. Best-seed selection (A04) auto-registers ancestry cluster groups to the sample registry for consumption by MODULE_5A inversion analysis.

**Theme system.** All visualization scripts source `theme_systems_plate.R` for consistent publication styling across the entire manuscript figure set.

## Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| K range | 2–20 | A01, A02, A04 |
| Seeds per K | 5 | A01, A02 |
| NGSadmix minMaf | 0.05 | A02 |
| evalAdmix nits | 5 | A03 |
| evalAdmix mistol | 0.05 | A03 |
| PCAngsd eigenvectors | 6 | A04 |
| PCAngsd EM iterations | 250 | A04 |
| Theta first-degree | 0.177 | B02, B03, B06 |
| Theta dup/MZ | 0.354 | B02, B03, B06 |
| Best-seed rule | loglik → mean\|resid\| → max\|resid\| → seed | A04 |
| Palette | 20-color stable palette | A04, B01–B07 |

## Canonical outputs

| File | Description | Used by |
|------|-------------|---------|
| `best_seed_by_K.tsv` | One row per K: selected best seed + metrics | MODULE_5A, MODULE_6 |
| `sample_main_ancestry_by_K.tsv` | Per-sample: dominant cluster, full Q, color per K | All downstream |
| `cluster_palette_by_K.tsv` | Stable palette per K × cluster | All figures |
| `k_hierarchy_tree.tsv` | Merge tree: which components merge K→K−1 | B05, manuscript |
| `combined_best_seed_by_K.tsv` | all + pruned with sample_set column | Manuscript |
| `figures/*.pdf` | Publication figures (B01–B07) | Manuscript |

## Dependencies

R (data.table, ggplot2, ggraph, ComplexHeatmap, ggalluvial, patchwork), NGSadmix, evalAdmix, PCAngsd, theme_systems_plate.R

## Sidecar convention

Every step writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands (written at start)
- **`{step}.results`** — all output file paths + descriptions (written at end)

## Methods

See [`docs/MODULE_2B_methods.md`](docs/MODULE_2B_methods.md) for manuscript-ready methods prose.
