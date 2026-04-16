# MODULE_4E — DELLY Breakend (BND) Calling

DELLY2 translocation/breakend discovery, site merging, cohort regenotyping, germline filtering, functional annotation, and publication plots. Follows the same 6-job SLURM chain as MODULE_4B–4D. Reuses shared markdup BAMs and exclusion BED.

## Why this module exists (for the inversion paper)

BND (breakend) records are DELLY's output when the caller sees split-read evidence but cannot assemble the full SV structure. Every inversion creates two BND calls — a `CT=3to3` (right-facing) and a `CT=5to5` (left-facing) junction. When DELLY's INV-typing misses one — because a breakpoint falls in a repeat, or supporting reads are marginal — the two BND records can still recover the inversion.

MODULE_5A2 STEP06 extracts all BND records with inversion-orientation tags, pairs 3to3+5to5 junctions by distance, cross-references against the INV catalog, and flags orphan BNDs: isolated breakpoint evidence where no INV was called. These are the hidden inversions the typer missed and are rescue candidates for MODULE_5A's discovery set.

> **Note on numbering:** MODULE_4E was originally 4F. The letter E was freed because DELLY2 INS calling was dropped — DELLY's insertion detection is unreliable for short-read data at ~5× coverage, and insertions are instead captured by Manta (MODULE_4F) and Clair3 small indels (MODULE_4A).

## Pipeline

```
Shared markdup BAMs (from 00_markdup, built by MODULE_4B A01)
  │
  ├─ A. SV calling (6-job SLURM chain) ──────────────────────────────
  │   A01  prep: validate reused markdup BAMs + exclude BED
  │   A02  DELLY call -t BND per sample (30 parallel × 4 threads)
  │   A03  merge → regenotype → merge cohort → subset 81 → germline filter
  │   A04  annotation layers (gene overlap, inter-/intra-chromosomal)
  │   A05  summary report
  │
  └─ B. Visualization ───────────────────────────────────────────────
      B01  BND-specific publication plots (circos-style, per-chr bars)
  │
  └─→ Final BND catalog + GT matrices
        consumed by MODULE_5A (BND inversion signal), manuscript
```

## Directory layout

```
MODULE_4E_BND_Delly/
  00_module4e_config.sh              ← central config (SV_TYPE=BND)
  README.md
  docs/
    MODULE_4E_methods.md             ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh         ← validate shared inputs
    SLURM_A02_delly_discovery.sh     ← DELLY call -t BND
    SLURM_A03_merge_genotype.sh      ← merge → regenotype → filter
    SLURM_A04_annotation_layers.sh   ← functional annotation
    SLURM_A05_summary_report.sh      ← counts + stats
    SLURM_B01_plot_results.sh        ← publication figures
  launchers/
    LAUNCH_module4e.sh               ← submit 6-job dependency chain
  utils/
    plot_BND_results.R               ← BND publication plots
```

## Usage

```bash
bash launchers/LAUNCH_module4e.sh
```

## Shared resources

Reuses from MODULE_4B: markdup BAMs (`00_markdup/`), exclusion BED, sample lists.

## Cross-module significance

BND calls capture inter-chromosomal translocations and intra-chromosomal rearrangements that may indicate inversion breakpoints. The MODULE_5A inversion detection pipeline consumes BND calls as supporting evidence — paired BND calls flanking an inversion candidate strengthen the hypothesis that the region represents a true polymorphic inversion rather than a population structure artifact.

## Parameters

Same framework as MODULE_4B–4D with `SV_TYPE=BND`. See `00_module4e_config.sh`.

## Methods

See [`docs/MODULE_4E_methods.md`](docs/MODULE_4E_methods.md) for manuscript-ready methods prose.
