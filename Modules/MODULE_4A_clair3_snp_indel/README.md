# MODULE_4A — Clair3 SNP & INDEL Discovery

Clair3 variant calling with WhatsHap phasing, five-phase post-processing pipeline (parse → phase blocks → rescue → population regenotyping → classification), marker handoff packaging, and downstream cohort-wide analysis suite. Produces phased per-sample VCFs, classified variant tables, marker packages, and population-level catalogs consumed by MODULE_6 (founder packs) and MODULE_CONSERVATION.

## Why this module exists (for the inversion paper)

Per-sample SNP and indel discovery with hard genotypes. This is a **different variant catalog from MODULE_2A**, not a duplicate:

- MODULE_2A uses ANGSD to call **biallelic SNPs** from genotype likelihoods across the cohort. That catalog is designed for low-coverage population genetics — Q estimation, admixture, Fst — where per-sample hard genotypes are never materialized.
- MODULE_4A uses Clair3 to call **per-sample SNPs and indels** with hard phased genotypes. ANGSD does not call indels, and its biSNP filter drops sites (multi-allelic, low-depth-in-individuals, outside callable regions) that matter for per-sample analyses.

The Clair3 catalog feeds two downstream analyses the inversion paper depends on:

1. **MODULE_6 founder-pack fragment classes** — the 100-INDEL sliding window framework is the backbone of founder-pack analysis. Without a dense, phase-aware indel catalog at cohort scale, founder-pack fragment classes cannot be resolved. `[CONFIRM: founder packs use indels only, not SNPs — based on memory context "100-INDEL windows".]`
2. **MODULE_CON deleterious-variant scoring** — SnpEff / SIFT4G / VESM annotations on Clair3 variants feed the breeding-concern classification (BC0–BC4) overlaid on inversion candidate regions to identify inversions harboring elevated deleterious burden.

The three-branch rescue system (STEP_P03 through P05) targets indels near variant hotspots where Clair3's confidence is marginal — disproportionately the breakpoint-proximal indels that matter for inversion analysis.

> Numbering gaps in this module (STEP_D02, STEP_X08–X09, STEP_X13–X19) are reserved slots for future steps; current production uses the steps listed in the scripts directory.

## Pipeline overview

```
Filtered BAMs from MODULE_1 (226 samples, ~9× mean)
  │
  ├─ D. Discovery ────────────────────────────────────────────────────
  │   D01  prepare inputs (per-chr BEDs, manifest, dispatch table)
  │   D02  Clair3 discovery (SLURM array: per sample × per chromosome)
  │   D03  validate outputs (VCF completeness check)
  │   D04  merge phase audit
  │
  ├─ P. Post-processing (5 SLURM phases) ────────────────────────────
  │   Phase 1: array ×226  → P01-P04+P02B  (per-sample)
  │     P01  parse & annotate VCF
  │     P02  local event block detection
  │     P02B phase-aware blocks (WhatsHap PS + read-pair rescue)
  │     P03  strong single-sample rescue
  │     P04  export weak indel candidates
  │
  │   Phase 2: single job → P05-P06  (population)
  │     P05  group weak candidates across samples
  │     P06  prepare regenotype catalog
  │
  │   Phase 3: array ×226  → P07A  (per-sample regenotyping)
  │     P07A regenotype one BAM against shared catalog
  │
  │   Phase 4: single job → P07B  (merge)
  │     P07B merge per-sample regenotype + cohort summary
  │
  │   Phase 5: array ×226  → P08-P10  (per-sample final)
  │     P08  six-class final variant classification
  │     P09  marker handoff package (flanking seqs, BED, master table)
  │     P10  publication figures
  │
  └─ X. Downstream (cohort-wide analysis) ────────────────────────────
      X01  build variant catalog (all chr, all samples)
      X02  phase block catalog
      X03  master annotation (gene/exon/CDS overlap)
      X04  per-sample summary
      X05  sample distance matrices
      X06  marker selection
      X07  gene/chromosome tables
      X10–X12  plots (main, extended, phase signatures)
      X20  annotate variant consequences (VEF)
      X21  score breeding concern
  │
  └─→ Classified variants + phase blocks + markers + population catalog
        consumed by MODULE_6, MODULE_CONSERVATION, manuscript
```

## Directory layout

```
MODULE_4A_SNP_INDEL_Clair3/
  README.md
  docs/
    MODULE_4A_methods.md                         ← manuscript-ready methods
    MODULE_CONSERVATION_RECIPE.sh                ← recipe for conservation pipeline
  steps/
    discovery/
      STEP_D01_prepare_inputs.sh                 ← BEDs, manifest, dispatch table
      STEP_D03_validate_outputs.sh               ← VCF completeness check
      STEP_D04_merge_phase_audit.sh              ← merge phase stats
    postprocess/
      STEP_P01_parse_and_annotate_vcf.py         ← parse Clair3 VCF + annotate
      STEP_P02_detect_local_event_blocks.py      ← local overlap blocks
      STEP_P02B_phase_aware_blocks.py            ← WhatsHap + read-pair phasing
      STEP_P03_rescue_strong_single_sample.py    ← strong filtered rescue
      STEP_P04_export_weak_indel_candidates.py   ← weak indel candidate export
      STEP_P05_group_weak_candidates.py          ← population clustering
      STEP_P06_prepare_regenotype_catalog.py     ← shared regenotype catalog
      STEP_P07A_regenotype_one_sample.py         ← per-sample BAM scan
      STEP_P07B_merge_regenotype.py              ← merge + cohort summary
      STEP_P08_final_classification.py           ← six-class classification
      STEP_P09_marker_handoff_package.py         ← flanking seqs, BED, master table
      STEP_P10_publication_figure.R              ← per-sample figures
    downstream/
      00_downstream_config.sh                    ← config for downstream suite
      STEP_X01_build_variant_catalog.py          ← cohort-wide catalog
      STEP_X02_phase_block_catalog.py            ← phase block catalog
      STEP_X03_master_annotation.py              ← gene/exon/CDS overlap
      STEP_X04_per_sample_summary.py             ← per-sample burden
      STEP_X05_sample_distance_matrices.py       ← distance matrices
      STEP_X06_marker_selection.py               ← marker selection
      STEP_X07_gene_chr_tables.py                ← gene/chr tables
      STEP_X10_plot_main.R                       ← main publication plots
      STEP_X11_plot_extended.R                   ← extended plots
      STEP_X12_plot_phase_signatures.R           ← phase signature plots
      STEP_X20_annotate_variant_consequences.py  ← VEF annotation
      STEP_X21_score_breeding_concern.py         ← breeding concern scoring
  slurm/
    SLURM_D02_clair3_discovery.sh                ← SLURM array: Clair3 per sample
    SLURM_D02b_clair3_discovery_allchr.sh        ← SLURM array: all sample×chr
    SLURM_P01_per_sample_phase1.sh               ← Phase 1: STEP P01-P04+P02B
    SLURM_P02_population_shared.sh               ← Phase 2: STEP P05-P06
    SLURM_P03_regenotype_array.sh                ← Phase 3: STEP P07A array
    SLURM_P04_regenotype_merge.sh                ← Phase 4: STEP P07B merge
    SLURM_P05_classify_array.sh                  ← Phase 5: STEP P08-P10
    SLURM_P06_population_steps_v2.sh             ← alternative population runner
  launchers/
    LAUNCH_module4a.sh                           ← master controller (full pipeline)
    LAUNCH_postprocessing.sh                     ← post-processing with resume
    LAUNCH_downstream.sh                         ← downstream analysis suite
    LAUNCH_downstream_plots.sh                   ← downstream plots only
    LAUNCH_downstream_vef.sh                     ← VEF annotation + scoring
    LAUNCH_batches_4ch.sh                        ← batch helper (4-chr groups)
  utils/
    run_pilot_CGA009.sh                          ← pilot test script
```

## Naming convention

Three phase prefixes reflect the pipeline's three major stages:

- Phase D = Discovery (Clair3 calling): `STEP_D01–D04`
- Phase P = Post-processing (per-sample + population): `STEP_P01–P10`
- Phase X = Downstream (cohort-wide analysis): `STEP_X01–X21`
- SLURM: `SLURM_{D|P}{NN}_{description}.sh`
- Launchers: `LAUNCH_{purpose}.sh`

## Usage

### Full pipeline

```bash
bash launchers/LAUNCH_module4a.sh --all_chroms
```

### Step by step

```bash
# 1. Prepare
bash steps/discovery/STEP_D01_prepare_inputs.sh

# 2. Discover (per-chromosome SLURM arrays)
for CHR in $(cat meta/chromosome_list.txt); do
    sbatch --array=1-226 slurm/SLURM_D02_clair3_discovery.sh --chrom "$CHR"
done

# 3. Validate
bash steps/discovery/STEP_D03_validate_outputs.sh --all_chroms

# 4. Post-processing
bash launchers/LAUNCH_postprocessing.sh --all_chroms
```

### Resume post-processing from any step

```bash
bash launchers/LAUNCH_postprocessing.sh --chrom C_gar_LG01 --step 7    # Phase 3+4+5
bash launchers/LAUNCH_postprocessing.sh --chrom C_gar_LG01 --step 8    # Phase 5 only
```

| `--step` | Runs |
|----------|------|
| 1 | All 5 phases (full run) |
| 5 | Phase 2+3+4+5 (from population clustering) |
| 7 | Phase 3+4+5 (regenotype array + merge + classify) |
| 7b | Phase 4+5 (merge only + classify) |
| 8 | Phase 5 only (classify + markers + figures) |

### Downstream

```bash
bash launchers/LAUNCH_downstream.sh           # full downstream
bash launchers/LAUNCH_downstream_plots.sh     # plots only
bash launchers/LAUNCH_downstream_vef.sh       # VEF annotation + breeding concern
```

## Clair3 parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| --platform | ilmn | Illumina short reads |
| --qual | 20 | PASS/LowQual threshold; post-processing rescues below |
| --snp_min_af | 0.08 | Low floor to capture candidates for rescue |
| --indel_min_af | 0.08 | Same rationale |
| --min_mq | 20 | Minimum mapping quality |
| --min_coverage | 2 | Allow low-coverage sites (~5× data) |
| --gvcf | enabled | For future joint genotyping |
| --use_whatshap_for_final_output_phasing | enabled | Produces PS tags + phased GTs |

## Post-processing architecture

The five SLURM phases implement a staged per-sample → population → per-sample workflow:

**Phase 1 (per-sample):** Parse Clair3 VCF, detect local event blocks, build phase-aware haplotype blocks (WhatsHap PS tags for TIER_1, read-pair union-find for TIER_2), rescue strong filtered variants, export weak indel candidates.

**Phase 2 (population):** Cluster weak indel candidates across all 226 samples, build a shared regenotype catalog of candidate sites.

**Phase 3 (per-sample):** Each sample independently scans its BAM against the shared catalog using a linear-pass CandidateIndex with binary search.

**Phase 4 (merge):** Concatenate per-sample regenotype results and compute cohort-level summary statistics.

**Phase 5 (per-sample):** Final six-class variant classification, marker handoff package (flanking sequences, BED coordinates, master table), and publication figures.

## Key design decisions

**WhatsHap phasing is required.** The `--use_whatshap_for_final_output_phasing` flag produces PS tags and pipe-separated phased genotypes that feed directly into STEP_P02B phase-aware blocks. Phase blocks are short at ~5× (typically 2–5 variants) but locally reliable.

**Three-branch rescue system.** STEP_P03 rescues strong filtered variants using evidence-based criteria. STEP_P04 exports weak candidates. STEP_P05–P07 use population-level regenotyping to validate weak candidates across the cohort.

**Six-class final classification.** STEP_P08 assigns each variant to one of six classes based on individual quality, population support, and regenotype evidence.

## Dependencies

Clair3 (Apptainer/Singularity container), samtools, pysam, python3, R (data.table, ggplot2)

## Sidecar convention

Every step writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands
- **`{step}.results`** — all output file paths + descriptions

## Methods

See [`docs/MODULE_4A_methods.md`](docs/MODULE_4A_methods.md) for manuscript-ready methods prose.
