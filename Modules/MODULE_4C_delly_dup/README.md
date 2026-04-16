# MODULE_4C — DELLY Duplication Calling

DELLY2 duplication discovery, site merging, cohort regenotyping, germline filtering, functional annotation, and downstream analysis. Follows the same architecture as MODULE_4B (DEL) with DUP-specific parameters. Reuses the shared markdup BAMs and exclusion BED from MODULE_4B.

## Pipeline

```
Shared markdup BAMs (from 00_markdup, built by MODULE_4B A01)
  │
  ├─ A. SV calling (4-job SLURM chain) ──────────────────────────────
  │   A01  prep: validate reused markdup BAMs + exclude BED
  │   A02  DELLY call -t DUP per sample (30 parallel × 4 threads)
  │   A03  merge → regenotype → merge cohort → subset 81 → germline filter
  │   A04  annotation layers (gene/exon/CDS, repeats, depth)
  │   A05  summary report
  │
  └─ B. Downstream ──────────────────────────────────────────────────
      B01  master annotation, per-sample burden, marker selection
      B02  extended DUP plots
  │
  └─→ Final DUP catalog + GT matrices + marker selections
        consumed by MODULE_5A, manuscript
```

## Directory layout

```
MODULE_4C_DUP_Delly/
  00_module4c_config.sh              ← central config
  README.md
  docs/
    MODULE_4C_methods.md             ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh         ← validate shared inputs
    SLURM_A02_delly_discovery.sh     ← DELLY call -t DUP
    SLURM_A03_merge_genotype.sh      ← merge → regenotype → filter
    SLURM_A04_annotation_layers.sh   ← functional annotation
    SLURM_A05_summary_report.sh      ← counts + stats
    SLURM_B01_downstream_analysis.sh ← master annotation + burden
    SLURM_B02_extended_plots.sh      ← extended figures
  launchers/
    LAUNCH_module4c.sh               ← submit SLURM dependency chain
    LAUNCH_downstream.sh             ← downstream analysis + plots
  utils/
    build_master_annotation.py       ← master DUP annotation table
    build_per_sample_summary.py      ← per-sample burden
    select_markers.py                ← tiered marker selection
    gene_summary_tables.py           ← gene/chr overlap tables
    plot_DUP_extended.R              ← extended DUP plots
    make_dup_pa_tables.py            ← DUP presence/absence tables
    check_dup_span.sh                ← DUP span QC helper
```

## Usage

```bash
# Full pipeline
bash launchers/LAUNCH_module4c.sh

# Downstream (after pipeline completes)
bash launchers/LAUNCH_downstream.sh
```

## Shared resources (from MODULE_4B)

This module does NOT re-create markdup BAMs or the exclusion BED. It reuses:
- `${BASE}/delly_sv/00_markdup/` — duplicate-marked BAMs
- `${BASE}/delly_sv/exclude.minimal.bed` — callable-based exclusion BED
- `${BASE}/delly_sv/samples_all_226.txt` / `samples_unrelated_81.txt`

MODULE_4B A01 must have been run at least once before running this module.

## Output structure

```
delly_sv_DUP/
  01_discovery/          Per-sample DUP BCFs
  02_merged_sites/       Shared site list
  03_genotyped/          Per-sample regenotyped BCFs
  04_merged_cohort/      226-sample merged BCF
  05_subset_81/          81-sample unrelated subset
  06_germline_filtered/  Germline-filtered BCF
  07_final_catalogs/     Final VCFs + GT matrices + BEDs
```

## Parameters

Same as MODULE_4B except SV type = DUP. See `00_module4c_config.sh` for all values.

## Dependencies

DELLY2 (v1.7.3), samtools, bcftools, bedtools, python3, R (data.table, ggplot2)

## Methods

See [`docs/MODULE_4C_methods.md`](docs/MODULE_4C_methods.md) for manuscript-ready methods prose.
