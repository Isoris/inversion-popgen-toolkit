# MODULE_3 — Heterozygosity & Runs of Homozygosity

Per-sample genome-wide heterozygosity, multiscale local theta diversity tracks, ngsF-HMM runs of homozygosity (ROH), FROH computation, and publication figure suite with ancestry/relatedness overlays. Produces the heterozygosity and inbreeding metrics consumed by the breeding/complementarity module and manuscript figures.

## Why this module exists (for the inversion paper)

Inversions suppress recombination in heterokaryotypes. The direct prediction: carriers of an inversion should show locally elevated observed heterozygosity across the inverted region relative to a Hardy-Weinberg expectation. This module produces the genome-wide heterozygosity baseline and the ngsF-HMM ROH catalog. `phase_qc_shelf/STEP_Q07b` + `STEP_Q07c` (the per-group Hobs pair that superseded MODULE_5E in April 2026) then asks, for each inversion candidate, whether the region shows the expected Hobs pattern against this baseline.

The per-sample F_ROH estimates also enter MODULE_6 (founder-pack analysis) as the inbreeding-coefficient covariate for directional gradient scoring — distinguishing genuinely rare founder-lineage signal from coincidental ROH-driven allele sharing.

## Pipeline

```
Filtered BAMs from MODULE_1 (226 samples, ~9× mean)
+ BEAGLE GLs from MODULE_2A
+ Ancestry labels from MODULE_2B
  │
  ├─ A. Compute ──────────────────────────────────────────────────────
  │   A01  validate inputs, build QC-pass BAM list + sample list
  │   A02  per-sample SAF → SFS → theta (main 500kb + multiscale)
  │        [sequential or SLURM array, one sample at a time]
  │   A03  ngsF-HMM (10 replicates, keep best by likelihood)
  │   A04  parse .ibd → BED ROH, compute FROH, het in/out ROH
  │
  └─ B. Visualization & reporting ────────────────────────────────────
      B01  core het plots, ROH plots, per-chr plots, ideograms,
           scatter stats, metadata overlays, statistics, report
  │
  └─→ Master summary tables + ROH BEDs + theta tracks + figures
        consumed by breeding module, MODULE_5A, manuscript
```

## Directory layout

```
MODULE_3_Heterozygosity_ROH/
  00_module3_config.sh                     ← central config
  README.md
  docs/
    MODULE_3_methods.md                    ← manuscript-ready methods prose
  steps/
    STEP_A01_prep_inputs.sh                ← validate inputs, build sample list
    STEP_A02_run_heterozygosity.sh         ← per-sample SAF → SFS → theta
    STEP_A03_run_ngsF_HMM.sh              ← ngsF-HMM multi-replicate
    STEP_A04_parse_roh_and_het.sh          ← parse ROH, compute FROH
    STEP_B01_run_all_plots.sh              ← all plots + stats + report
  slurm/
    SLURM_A02_heterozygosity_worker.sh     ← SLURM array: one sample per task
  launchers/
    LAUNCH_module3.sh                      ← main wrapper
  utils/
    parse_roh_and_het.py                   ← .ibd → BED ROH + FROH + het in/out
    write_report.py                        ← auto-generate Methods/Results markdown
    plot_heterozygosity_core.R             ← genome-wide het distributions
    plot_roh_core.R                        ← ROH/FROH genome-wide plots
    plot_roh_by_chromosome.R               ← per-chromosome ROH heatmaps
    plot_roh_metadata_overlays.R           ← ancestry/family grouped plots
    plot_scatter_stats.R                   ← het vs ROH/FROH/depth scatters
    plot_theta_ideogram.R                  ← local theta ideogram tracks
    run_stats.R                            ← sample + chromosome level stats
```

## Naming convention

Matches inversion codebase v8.5 and MODULE_2A/2B:

- Config: `00_module3_config.sh`
- Steps: `STEP_{PHASE}{NN}_{description}.sh` — Phase A = compute, Phase B = visualization
- SLURM: `SLURM_{PHASE}{NN}_{description}.sh`
- Launchers: `LAUNCH_module3.sh`
- Helper functions: `hr_log`, `hr_die`, `hr_check_file`, `hr_init_dirs`

## Usage

```bash
# Run everything sequentially
bash launchers/LAUNCH_module3.sh

# Run only step 2
bash launchers/LAUNCH_module3.sh --step 2

# Run from step 3 onward
bash launchers/LAUNCH_module3.sh --from 3

# Submit step 2 as SLURM array (one sample per task)
bash launchers/LAUNCH_module3.sh --step 2 --slurm
# Then after SLURM completes:
bash launchers/LAUNCH_module3.sh --from 3
```

## Practical run order

1. `--step 1` — validate inputs, build QC-pass BAM and sample lists
2. `--step 2 --slurm` — submit per-sample heterozygosity as SLURM array (226 jobs)
3. *(wait for SLURM completion)*
4. `--from 3` — ngsF-HMM → ROH parsing → plots + stats + report

## Key design decisions

**Per-sample SAF/SFS.** Each sample gets its own SAF and folded 1D-SFS via ANGSD/realSFS, from which per-site theta (Watterson's θ_W and pairwise θ_π) is computed. This is the standard approach for individual-level heterozygosity from low-coverage data.

**Multiscale theta windows.** In addition to the main 500 kb non-overlapping windows, three finer scales (5 kb/1 kb step, 10 kb/2 kb step, 50 kb/10 kb step) are computed for local diversity landscape visualization. These feed the ideogram plots and downstream inversion support analyses.

**ngsF-HMM for ROH.** Runs of homozygosity are called using ngsF-HMM with 10 random restarts (seeds 42–51), keeping the replicate with the highest log-likelihood. ngsF-HMM operates on genotype likelihoods directly (from BEAGLE format), avoiding hard-call genotyping bias at low coverage.

**Dual ROH bin system.** ROH are classified into length bins (short/medium/long) reflecting different demographic processes — recent inbreeding produces long ROH, ancestral bottlenecks produce short ROH. Both absolute length and FROH (proportion of genome in ROH) are reported.

**Breeding module separation.** This module produces descriptive het/ROH/FROH metrics only. Breeding value scoring, complementarity analysis, and cross-prediction are handled by a separate downstream module that consumes MODULE_3 outputs.

## Output directory structure

```
het_roh/
  01_inputs_check/     ← QC-pass BAM list, sample list, BEAGLE, .pos, .ind
  02_heterozygosity/   ← per-sample SAF, SFS, theta (main + multiscale)
  03_ngsF_HMM/         ← ngsF-HMM outputs (10 reps, best selected)
  04_roh_summary/      ← parsed ROH BEDs, FROH tables, het in/out ROH
  05_inversion_support/ ← theta tracks for inversion region support
  06_plots_core/       ← het + ROH core figures
  07_plots_metadata/   ← ancestry/relatedness overlay figures
  08_stats/            ← statistical summaries
  09_final_tables/     ← master summary tables
  10_report/           ← auto-generated manuscript sections
```

## Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| ANGSD -minQ | 20 | A02 |
| ANGSD -minMapQ | 30 | A02 |
| ANGSD -C | 50 | A02 |
| realSFS maxIter | 2000 | A02 |
| realSFS tolerance | 1e-16 | A02 |
| thetaStat main window | 500 kb / 500 kb step | A02 |
| thetaStat multiscale | 5kb/1kb, 10kb/2kb, 50kb/10kb | A02 |
| ngsF-HMM replicates | 10 | A03 |
| ngsF-HMM seed base | 42 | A03 |

## Dependencies

ANGSD (doSaf, realSFS, thetaStat), ngsF-HMM, samtools, python3, R (data.table, ggplot2)

## Sidecar convention

Every step writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands (written at start)
- **`{step}.results`** — all output file paths + descriptions (written at end)

## Methods

See [`docs/MODULE_3_methods.md`](docs/MODULE_3_methods.md) for manuscript-ready methods prose.
