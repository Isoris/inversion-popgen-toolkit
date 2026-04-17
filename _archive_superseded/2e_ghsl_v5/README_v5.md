# phase_2_discovery/2e_ghsl — Snake 3 / GHSL v5

Status as of audit 2026-04-17. **FINDING 19 upgraded to FIX 22 in
chat 4 continuation — 2e IS now wired into C01d's D10 scoring.**

## What this module does

`STEP_C04_snake3_ghsl_v5.R` computes within-sample haplotype
divergence per window from phased Clair3 VCFs. For each sample, it
measures how different the two haplotypes are locally; the population
distribution of those divergences is scored for rank stability and
bimodality. Inversion windows are expected to show stable, bimodal
divergence patterns (INV/INV samples consistently low, INV/nonINV
samples consistently high). Output is per-window scores plus
per-sample karyotype calls.

`STEP_C04b_snake3_ghsl_figures.R` is the plotting companion.

## Wiring (FIX 22)

`run_all.R` in 2d_candidate_detection now reads
`<ghsl-dir>/annot/<chr>.ghsl_v5.annot.rds` (if present) during
Phase 8 and aggregates per-block over `global_window_id` in
`[block.start, block.end]`. Five columns land on the scoring table:

- `ghsl_v5_score_max`        — peak composite score within block
- `ghsl_rank_stability_max`  — peak inter-window rank stability
- `ghsl_div_contrast_z_max`  — peak bimodality-contrast z-score
- `ghsl_div_bimodal_frac`    — fraction of windows with >=2 density modes
- `ghsl_pass_frac`           — fraction of windows with score > 0.65
- `ghsl_n_scored_windows`    — non-NA-scored windows in the block

C01d's D10 dimension (partition stability) now blends:
- When Clair3-phased data exists for the block: `0.60 * d10_ghsl +
  0.40 * d10_simmat` where `d10_ghsl = 0.50 * ghsl_v5_score_max +
  0.25 * ghsl_div_bimodal_frac + 0.25 * ghsl_pass_frac`.
- When only sim_mat data: `d10_simmat` alone (pre-FIX-22 behaviour).

C01d emits `d10_source` in its scoring table (values:
`simmat_only`, `ghsl_and_simmat`) so you can see which path fired
per candidate.

## What IS still NOT wired

- C04 doesn't consume the C01i group registry. It runs on phased
  Clair3 VCFs directly. If C04 were rewritten to accept registry
  groups it could produce per-karyotype-class validation. Deferred
  (noted by Quentin 2026-04-17 as "fix later").
- 2c_precomp README's "Layer C = GHSL" framework is now mostly
  accurate — C04's Clair3 signal does reach C01d's scoring
  (FIX 22), though D05's sim_mat-based partition also still feeds
  D10 as a secondary signal. Both paths are present.

## 4-layer independence framework status

- **Layer A (dosage/MDS):** live via C01a precomp → D01-D08
- **Layer B (SV calls):** live via C00 sv_prior → C01g boundary catalog
- **Layer C (GHSL):** live via C04 → run_all Phase 8 → C01d D10 (FIX 22)
- **Layer D (genotype-breakpoint association):** via hypothesis tests
  in C01f (phase 4c)

## Running 2e

    Rscript STEP_C04_snake3_ghsl_v5.R <precomp_dir> <ghsl_prep_dir> <outdir> \
      [--chrom C_gar_LG01] [--test_windows 50] [--ncores 4]

Requires Clair3 phased VCFs per chromosome. As of 2026-04-17,
Clair3 is done through LG12; LG13-16 running; LG17-27 + LG02-06
stragglers over ~3 days. C04 can be run per-chromosome as each
finishes phasing; Phase 8 of run_all.R picks up whatever annot.rds
shards exist at the time of the detection run.
