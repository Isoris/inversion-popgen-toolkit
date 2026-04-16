# MODULE_4F — DELLY Translocation (TRA) Calling

DELLY2 translocation discovery, site merging, cohort regenotyping, germline filtering, functional annotation, and publication plots. Follows the same 6-job SLURM chain as MODULE_4B–4E. Reuses shared markdup BAMs and exclusion BED.

> **Note on numbering:** MODULE_4F was originally 4G. Letters were shifted because DELLY2 INS calling was dropped (unreliable at ~5× short-read coverage). The current MODULE_4 series is: 4A (Clair3), 4B (DEL), 4C (DUP), 4D (INV), 4E (BND), 4F (TRA), 4G (Manta).

## Pipeline

```
Shared markdup BAMs (from 00_markdup, built by MODULE_4B A01)
  │
  ├─ A. SV calling (6-job SLURM chain) ──────────────────────────────
  │   A01  prep: validate reused markdup BAMs + exclude BED
  │   A02  DELLY call -t TRA per sample (30 parallel × 4 threads)
  │   A03  merge → regenotype → merge cohort → subset 81 → germline filter
  │   A04  annotation layers (gene overlap, inter-chromosomal pairs)
  │   A05  summary report
  │
  └─ B. Visualization ───────────────────────────────────────────────
      B01  TRA-specific publication plots
  │
  └─→ Final TRA catalog + GT matrices
        consumed by manuscript
```

## Directory layout

```
MODULE_4F_TRA_Delly/
  00_module4f_config.sh              ← central config (SV_TYPE=TRA)
  README.md
  docs/
    MODULE_4F_methods.md             ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh         ← validate shared inputs
    SLURM_A02_delly_discovery.sh     ← DELLY call -t TRA
    SLURM_A03_merge_genotype.sh      ← merge → regenotype → filter
    SLURM_A04_annotation_layers.sh   ← functional annotation
    SLURM_A05_summary_report.sh      ← counts + stats
    SLURM_B01_plot_results.sh        ← publication figures
  launchers/
    LAUNCH_module4f.sh               ← submit 6-job dependency chain
  utils/
    plot_TRA_results.R               ← TRA publication plots
```

## Usage

```bash
bash launchers/LAUNCH_module4f.sh
```

## Shared resources

Reuses from MODULE_4B: markdup BAMs (`00_markdup/`), exclusion BED, sample lists.

## TRA vs BND

DELLY distinguishes TRA (inter-chromosomal translocations) from BND (general breakends including intra-chromosomal). TRA events are a subset of BND where both breakpoint mates map to different chromosomes. Both are called separately because DELLY applies different evidence models and filtering for each type.

## Parameters

Same framework as MODULE_4B–4E with `SV_TYPE=TRA`. See `00_module4f_config.sh`.

## Methods

See [`docs/MODULE_4F_methods.md`](docs/MODULE_4F_methods.md) for manuscript-ready methods prose.
