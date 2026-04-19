# MODULE_QC_ShelfDiagnosis

Systematic QC tracks to diagnose whether a Z plateau on a local-PCA profile
(e.g. LG28 at 15–18 Mb) is a real biological signal or a data-quality artifact.

## Layout

```
inversion_modules/phase_qc_shelf/
├── 00_config.sh                      # paths + tool bins, sourced by every script
├── install.sh                        # one-shot setup + Engine B/F compile
├── run_all.sh                        # driver: Q01..Q07 + Q04
├── export_precomp_to_json.R          # build scrubber JSON (theta + ancestry)
├── STEP_Q01_snp_density.sh           # SNP density per bin from .pos.fixed
├── STEP_Q02_beagle_uncertainty.sh    # BEAGLE posterior uncertainty per bin
├── STEP_Q03_coverage_tracks.sh       # per-bin coverage from mosdepth
├── STEP_Q04_plot_diagnostic.sh       # stacked up-to-10-track diagnostic PDF
├── STEP_Q05_aggregate_theta.sh       # pestPG -> per-sample theta_pi long TSV
├── STEP_Q06_ancestry_tracks.sh       # Engine B local_Q -> per-window ancestry
├── STEP_Q06_multiscale.sh            # Engine B at 1x/5x/10x window scales
├── STEP_Q07_popstats.sh              # Engine F (Fst/dXY/theta/Tajima) per group
├── R/
│   ├── q02_beagle_stream.R
│   ├── q03_coverage_collapse.R
│   ├── q04_compose_plot.R
│   ├── q05_aggregate_theta.R
│   ├── q06_ancestry_tracks.R
│   └── q07_build_groups.R
├── slurm/
│   └── array_28chrom.sh
└── docs/
    └── interpretation.md
```

## What the Q-steps produce

| Step | Compute | Key output | Diagnostic |
|---|---|---|---|
| Q01 | SNPs/window from .pos.fixed | SNP density track | low = mapping gap / low polymorphism |
| Q02 | BEAGLE posterior uncertainty | uncertain_frac track | high = caller lacked confidence |
| Q03 | mosdepth coverage | mean_cov + CV track | low+low CV = uniform mappability drop (repeat) |
| Q05 | ANGSD pestPG aggregation | per-sample θπ (heterozygosity) | low+uniform = low polymorphism; low+stratified = real inversion |
| Q06 | Engine B local_Q | per-window Δ12, ENA, maxQ | stable Δ12 = real structure; spiky = noise |
| Q06-multi | Engine B at 1x/5x/10x | Δ12 at 3 scales | scale-stable = robust; scale-dependent = label switching |
| Q07 | Engine F region_popstats | Fst/dXY/θ/Tajima per group | Hom1-vs-Hom2 Fst spike = inversion signature |
| Q04 | composite plot | 11-panel stacked PDF | shelf vs ref summary + interpretation guide |
| Q08 | sample × SNP genotype heatmap for shelf, 9 multi-scale panels | per-panel PDF + composite | 3-band pattern across all scales = robust real inversion; scale-dependent = noise |

## Quick start — LG28 with shelf

```bash
cd ${BASE}/inversion-popgen-toolkit/inversion_modules/phase_qc_shelf
bash install.sh                              # one-time: paths + compile

SHELF_START_MB=15 SHELF_END_MB=18 \
  bash run_all.sh C_gar_LG28
```

## Engine B (Q06) and Engine F (Q07) explained

**Engine B = `instant_q`**: C++ binary, per-window ancestry Q estimation via
fixed-F EM on BEAGLE GL. Outputs Δ12/entropy/ENA per window + per-sample maxQ.

**Engine F = `region_popstats`**: C binary, per-window Hudson Fst, dXY, dA, θ_π,
θ_W, Tajima's D across 1 or more sample groups. Reads BEAGLE, runs at fixed
windowed scale (default 50kb/10kb to match pestPG).

Q07 uses Engine F with **two grouping strategies** in separate runs:

1. **Ancestry groups** — samples bucketed by their dominant maxQ from Engine B.
   Per-group Fst answers "is the signal driven by ancestry structure?"
2. **Inversion genotypes (Hom1/Het/Hom2)** — k-means(3) on PC1 loadings
   averaged over shelf windows, then Fst across the whole chromosome. The
   classical inversion signature: Fst spikes at the shelf and drops outside.

## Multi-scale ancestry

See INSTALL.md.

## Scrubber (pca_scrubber.html)

Build the JSON with ancestry integration:
```bash
Rscript export_precomp_to_json.R \
  --precomp  ${PRECOMP_DIR}/C_gar_LG28.precomp.rds \
  --samples  samples_with_ancestry.tsv \
  --theta    results/tracks/theta.C_gar_LG28.win50000.step10000.tsv.gz \
  --ancestry results/tracks/ancestry_window.C_gar_LG28.tsv \
  --ancestry_samples results/tracks/ancestry_sample.C_gar_LG28.maxQ_wide.tsv.gz \
  --out      LG28_scrubber.json
```

## Skip patterns

```bash
# Reuse everything, only regenerate the plot:
SKIP_Q01=1 SKIP_Q02=1 SKIP_Q03=1 SKIP_Q05=1 SKIP_Q06=1 SKIP_Q06_MULTI=1 SKIP_Q07=1 \
  SHELF_START_MB=15 SHELF_END_MB=18 \
  bash run_all.sh C_gar_LG28
```

## Dependencies

- R with data.table, ggplot2, patchwork
- mosdepth in PATH (only if RUN_MOSDEPTH=1 in Q03)
- Existing ANGSD BEAGLE + .pos.fixed (from `inversion_localpca_v7/02_snps_beagle`)
- Existing mosdepth output at `het_roh/mosdepth_output` (or set MOSDEPTH_EXISTING)

## What the tracks tell you

See `docs/interpretation.md`.
