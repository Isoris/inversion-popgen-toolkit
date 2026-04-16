# MODULE_4B — DELLY Deletion Calling

DELLY2 deletion discovery, site merging, cohort regenotyping, germline filtering, functional annotation, depth support QC, and publication figure suite. Produces the final DEL catalog (226-sample + 81-unrelated subsets) with GT matrices, marker selections, and per-sample burden tables consumed by downstream SV modules and the manuscript.

## Pipeline

```
Markdup BAMs from MODULE_1 (226 samples, raw minimap2 → markdup)
  │
  ├─ A. SV calling pipeline (6-job SLURM dependency chain) ──────────
  │   A01  prep inputs (markdup BAMs, exclude BED, sample lists)
  │   A02  DELLY call -t DEL per sample (30 parallel × 4 threads)
  │   A03  merge sites → regenotype → merge cohort → subset 81 → germline filter
  │   A04  annotation layers (gene/exon/CDS, repeats, depth, mate QC)
  │   A05  summary report (counts, per-sample stats, spectrum)
  │
  ├─ B. Visualization & downstream ──────────────────────────────────
  │   B01  publication plots (genome heatmap, PCA, sharing, size dist)
  │   B02  downstream analysis (master annotation, burden, markers)
  │   B03  extended plots (ancestry overlays, detailed breakdowns)
  │
  └─→ Final DEL catalog + GT matrices + marker selections + figures
        consumed by MODULE_5A (inversion support), MODULE_6, manuscript
```

## Directory layout

```
MODULE_4B_DEL_Delly/
  00_module4b_config.sh                    ← central config
  README.md
  exclude.minimal.bed                      ← callable-based exclusion BED
  samples_all_226.txt                      ← full sample list
  samples_unrelated_81.txt                 ← pruned sample list
  samples_81.txt                           ← alias
  docs/
    MODULE_4B_methods.md                   ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh               ← markdup BAMs, exclude BED
    SLURM_A02_delly_discovery.sh           ← DELLY call -t DEL per sample
    SLURM_A03_merge_genotype.sh            ← merge → regenotype → filter
    SLURM_A04_annotation_layers.sh         ← functional + repeat + depth
    SLURM_A05_summary_report.sh            ← counts, stats, spectrum
    SLURM_B01_plot_results.sh              ← publication figures
    SLURM_B02_downstream_analysis.sh       ← master annotation + burden
    SLURM_B03_extended_plots.sh            ← extended figure suite
  launchers/
    LAUNCH_module4b.sh                     ← submit 6-job dependency chain
    LAUNCH_downstream.sh                   ← downstream analysis + plots
  utils/
    plot_delly_results.R                   ← core publication plots
    build_master_annotation.py             ← master DEL annotation table
    build_per_sample_summary.py            ← per-sample burden table
    select_markers.py                      ← tiered marker selection
    gene_summary_tables.py                 ← gene/chromosome overlap tables
    plot_DEL_extended.R                    ← extended analysis plots
    summarize_results_for_claude.sh        ← compact summary for debugging
    visualize_results.sh                   ← IGV visualization helper
```

## Naming convention

- Config: `00_module4b_config.sh`
- SLURM: `SLURM_{A|B}{NN}_{description}.sh` — Phase A = SV calling, Phase B = viz/downstream
- Launchers: `LAUNCH_{purpose}.sh`
- Utils: standalone tools (R, Python) called by SLURM scripts
- Helper functions: `dv_log`, `dv_die`, `dv_check_file`, `dv_init_dirs`

## Usage

```bash
# Full pipeline (6-job SLURM chain with dependencies)
bash launchers/LAUNCH_module4b.sh

# Downstream analysis (after main pipeline completes)
bash launchers/LAUNCH_downstream.sh
```

The launcher submits: A01 → A02 → A03 → A04 → A05 → B01, each with `--dependency=afterok` on the previous job.

## Key design decisions

**DELLY on raw markdup BAMs.** DELLY requires unfiltered (not MAPQ/TLEN-filtered) BAMs with duplicate marking. MODULE_1's population-genomics-filtered BAMs are NOT used. Instead, A01 runs `samtools markdup` on the raw minimap2 BAMs (manifest column 2) to produce markdup copies.

**Callable-based exclusion BED.** Rather than scanning for N-blocks, the exclude BED is built empirically from 50-kb bins with callable_bp < 500 (from PA-Roary mosdepth), plus unconditional 50-kb chromosome-end masking on every chromosome.

**Two-tier catalog.** The full 226-sample merged BCF and an 81-sample unrelated subset are both produced. Germline filtering (`delly filter -f germline`) is applied to the 81-sample subset for population-level analyses.

**Strict marker selection.** Three-tier marker grading: Tier 1 (PASS + PRECISE + QUAL>500 + PE>4), Tier 2 (relaxed quality), Tier 3 (all passing germline filter).

## Output directory structure

```
delly_sv/
  00_markdup/            ← duplicate-marked BAMs (from raw minimap2)
  01_discovery/          ← per-sample DEL BCFs
  02_merged_sites/       ← shared site list
  03_genotyped/          ← per-sample regenotyped BCFs
  04_merged_cohort/      ← 226-sample merged BCF
  05_subset_81/          ← 81-sample unrelated subset
  06_germline_filtered/  ← germline-filtered BCF
  07_final_catalogs/     ← final VCFs + GT matrices + BEDs
  08_annotation/         ← functional + repeat annotation
  09_depth_support/      ← mosdepth depth QC
  10_mate_distance_qc/   ← SVLEN-based flags
  11_summary/            ← report + counts
  12_plots/              ← publication figures
```

## Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| DELLY version | 1.7.3 (compiled) | A02 |
| DELLY -t | DEL | A02 |
| DELLY -h | 4 (threads per call) | A02 |
| Parallel calls | 30 | A02 |
| Chr-end mask | 50 kb | A01 |
| Excl min callable_bp | 500 (per 50-kb bin) | A01 |
| Excl min block | 50 kb | A01 |
| Depth window | 500 bp | A04 |
| Depth MAPQ | 30 | A04 |
| Strict DEL QUAL | > 500 | config |
| Strict DEL PE | > 4 | config |
| Mate warn | 20 kb | config |
| Mate suspicious | 50 kb | config |

## Dependencies

DELLY2 (v1.7.3), samtools, bcftools, mosdepth, bedtools, HTSlib, Boost, python3, R (data.table, ggplot2, cowplot, viridis, ComplexHeatmap)

## Sidecar convention

Every step writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands
- **`{step}.results`** — all output file paths + descriptions

## Methods

See [`docs/MODULE_4B_methods.md`](docs/MODULE_4B_methods.md) for manuscript-ready methods prose.
