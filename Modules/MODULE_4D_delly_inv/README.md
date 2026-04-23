# MODULE_4D — DELLY Inversion Calling

DELLY2 inversion discovery, site merging, cohort regenotyping, germline filtering, functional annotation, and publication plots. Follows the same 6-job SLURM chain as MODULE_4B/4C with INV-specific parameters. Reuses shared markdup BAMs and exclusion BED. The DELLY INV catalog provides an independent SV-caller complement to the local-PCA-based inversion detection in MODULE_5A.

## Why this module exists (for the inversion paper)

Direct SV-caller evidence for inversions. Every candidate that MODULE_5A flags from population-level signal (local PCA, GHSL v5 within-sample haplotype divergence, 4-layer consensus combining Layer A dosage and Layer C haplotype contrast) should ideally have a concordant DELLY INV call with split-read support across both breakpoints. Conversely, DELLY INV calls with no MODULE_5A population signal are candidates for low-frequency inversions below the detection limit of local PCA.

The DELLY INV catalog anchors MODULE_5A2 (breakpoint validation), which cross-references every INV call against BAM-level evidence (pysam split-read orientation, INV3 vs INV5 tags) and tests concordance between DELLY and Manta (MODULE_4G) calls using Fisher / χ² / Cochran-Armitage tests.

## Pipeline

```
Shared markdup BAMs (from 00_markdup, built by MODULE_4B A01)
  │
  ├─ A. SV calling (6-job SLURM chain) ──────────────────────────────
  │   A01  prep: validate reused markdup BAMs + exclude BED
  │   A02  DELLY call -t INV per sample (30 parallel × 4 threads)
  │   A03  merge → regenotype → merge cohort → subset 81 → germline filter
  │   A04  annotation layers (gene/exon/CDS, repeats, depth)
  │   A05  summary report
  │
  └─ B. Visualization ───────────────────────────────────────────────
      B01  INV-specific publication plots
  │
  └─→ Final INV catalog + GT matrices
        consumed by MODULE_5A (DELLY anchor for inversion validation),
        MODULE_5B (breakpoint validation), manuscript
```

## Directory layout

```
MODULE_4D_INV_Delly/
  00_module4d_config.sh              ← central config (SV_TYPE=INV)
  README.md
  docs/
    MODULE_4D_methods.md             ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh         ← validate shared inputs
    SLURM_A02_delly_discovery.sh     ← DELLY call -t INV
    SLURM_A03_merge_genotype.sh      ← merge → regenotype → filter
    SLURM_A04_annotation_layers.sh   ← functional annotation
    SLURM_A05_summary_report.sh      ← counts + stats
    SLURM_B01_plot_results.sh        ← publication figures
  launchers/
    LAUNCH_module4d.sh               ← submit 6-job dependency chain
  utils/
    plot_INV_results.R               ← INV publication plots (v2)
    plot_INV_results_v1.R            ← legacy v1 plots
```

## Usage

```bash
bash launchers/LAUNCH_module4d.sh
```

## Shared resources

Reuses from MODULE_4B: markdup BAMs (`00_markdup/`), exclusion BED, sample lists.

## Cross-module significance

The DELLY INV catalog is consumed by MODULE_5A as an independent anchor for inversion validation. Candidates detected by local PCA (MODULE_5A) that also have DELLY INV support receive higher confidence scores. The DELLY INV GT matrix feeds into the hypothesis testing framework (STEP_C01f) as "Cheat 9" genotype concordance evidence.

## Parameters

Same framework as MODULE_4B/4C with `SV_TYPE=INV`. See `00_module4d_config.sh`.

## Methods

See [`docs/MODULE_4D_methods.md`](docs/MODULE_4D_methods.md) for manuscript-ready methods prose.
