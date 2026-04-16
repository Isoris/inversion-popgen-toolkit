# MODULE_4G — Manta SV Calling (All Types)

Manta structural variant discovery across all SV types (DEL, DUP, INS, INV/BND) in a single per-sample run, with inversion conversion, cohort merging, per-type splitting, evidence-based tiered filtering, and summary reporting. Complements the type-specific DELLY pipelines (MODULE_4B–4F) with Manta's assembly-based approach and insertion detection.

> **Note on numbering:** MODULE_4G was originally 4H. Letters were shifted because DELLY2 INS calling was dropped. The final MODULE_4 series is: 4A (Clair3), 4B (DEL), 4C (DUP), 4D (INV), 4E (BND), 4F (TRA), 4G (Manta).

## Pipeline

```
Shared markdup BAMs (from 00_markdup, built by MODULE_4B A01)
  │
  ├─ A. SV calling (5-job SLURM chain) ──────────────────────────────
  │   A01  prep: call regions BED, validate markdup BAMs
  │   A02  Manta configManta.py + runWorkflow.py per sample
  │        (30 parallel × 4 threads, minCandidateVariantSize=50)
  │   A03  per-sample INV conversion → bcftools merge → subset 81
  │        → split by type (DEL/DUP/INS_small/INS_large/INV/BND)
  │   A04  annotation + evidence-based tiered filtering
  │        (3 tiers: lenient / publication / strict)
  │   A05  summary report (per-sample burden, size dist, sharing)
  │
  └─ B. Visualization ───────────────────────────────────────────────
      B01  Manta publication plots
  │
  └─→ Per-type catalogs + GT matrices + tiered filtered sets
        consumed by MODULE_5A, MODULE_6, manuscript
```

## Directory layout

```
MODULE_4G_ALL_Manta/
  00_module4g_config.sh                  ← central config
  README.md
  configManta.custom.ini                 ← minCandidateVariantSize = 50
  call_regions.bed(.gz, .gz.tbi)         ← callable regions for Manta
  exclude.minimal.bed                    ← exclusion BED
  exclude.sorted.bed                     ← sorted exclusion BED
  genome.bed                             ← full genome intervals
  samples_all_226.txt / samples_*.txt    ← sample lists
  docs/
    MODULE_4G_methods.md                 ← manuscript-ready methods
  slurm/
    SLURM_A01_prep_inputs.sh             ← validate inputs, call regions
    SLURM_A02_manta_discovery.sh         ← Manta per-sample discovery
    SLURM_A03_merge_and_split.sh         ← convert INV → merge → split
    SLURM_A04_annotation_filter.sh       ← annotation + 3-tier filtering
    SLURM_A04b_filter_all_sv_types.sh    ← standalone filter wrapper
    SLURM_A05_summary_report.sh          ← summary + burden tables
    SLURM_B01_plot_results.sh            ← publication figures
  launchers/
    LAUNCH_module4g.sh                   ← submit 5-job dependency chain
  utils/
    convertInversion_py3.py              ← Manta BND→INV conversion
    plot_manta_results.R                 ← publication plots
    summarize_results_for_claude.sh      ← compact summary for debugging
```

## Usage

```bash
bash launchers/LAUNCH_module4g.sh
```

## Key differences from DELLY modules

**Single-run all types.** Manta calls DEL, DUP, INS, INV, and BND in a single per-sample run, then post-processing splits into per-type catalogs. DELLY runs separate pipelines per SV type.

**No merge/regenotype.** Manta operates per-sample (configManta.py + runWorkflow.py). Cohort-level analysis comes from bcftools merge of per-sample VCFs, not DELLY-style site-list regenotyping.

**Inversion conversion.** Manta encodes inversions as BND pairs with INV3/INV5 tags. `convertInversion_py3.py` converts these to standard INV records. This must run per-sample BEFORE merging — compound IDs in merged VCFs crash the converter.

**Insertion classes.** Manta uniquely detects insertions, split into two classes:
- Small INS: fully assembled, SVLEN present (50–~200 bp)
- Large INS: incompletely assembled, LEFT_SVINSSEQ/RIGHT_SVINSSEQ present, unknown total length

**Three-tier filtering.** Evidence-based filtering on PR_alt + SR_alt (summed across carriers), per-carrier minimum evidence, QUAL, size, and carrier count. IMPRECISE variants are flagged but not removed (real SVs at ~9× are often IMPRECISE due to low coverage).

## Output structure

```
MODULE_4G_ALL_Manta/
  01_per_sample/          ← raw Manta output per sample
  01b_per_sample_converted/ ← INV-converted per-sample VCFs
  02_merged_cohort/       ← 226-sample merged VCF
  03_subset_81/           ← 81-sample unrelated subset
  04_split_by_type/       ← per-type VCFs (DEL/DUP/INS_small/INS_large/INV/BND)
  05_final_catalogs/      ← PASS VCFs + GT matrices
  06_annotation/          ← functional annotation
  07_depth_support/       ← mosdepth QC
  08_summary/             ← summary tables + report
  09_plots/               ← publication figures
  10_filtered_catalogs/   ← tiered filter output
```

## Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| Manta minCandidateVariantSize | 50 bp | A02 |
| Manta threads per sample | 4 | A02 |
| Parallel samples | 30 | A02 |
| Strict min QUAL | 20 | config |
| Tier 1 (lenient) | min evidence, any precision | A04 |
| Tier 2 (publication) | good evidence + carriers + QUAL | A04 |
| Tier 3 (strict) | strong evidence + many carriers | A04 |

## Dependencies

Manta (manta_py2 conda env), samtools, bcftools, bedtools, python3 (for convertInversion), R (data.table, ggplot2)

## Methods

See [`docs/MODULE_4G_methods.md`](docs/MODULE_4G_methods.md) for manuscript-ready methods prose.
