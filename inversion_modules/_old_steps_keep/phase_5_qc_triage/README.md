# `phase_5_qc_triage/` — data-quality QC on candidate shelves

(Formerly the top-level `phase_qc_shelf/` module. Folded into phase_4 as
sub-block 4b in pass 12, 2026-04-24. Promoted to top-level phase_5 in
pass 15, 2026-04-24.)

Systematic QC tracks to diagnose whether a Z-plateau on a local-PCA
profile (e.g. LG28 at 15–18 Mb) is a real biological signal or a
data-quality artifact. Produces the soft `q_qc_shelf_flag` (one of
`clean` / `low_snp` / `high_uncertain` / `coverage_artifact` / `messy`)
that rides along on the final candidate catalog via pass-13's reader
in `phase_9_classification/`.

**Soft gate.** A messy candidate still gets characterized and tiered —
the flag just tells downstream readers to de-prioritize that candidate
without dropping it. See `docs/MODULE_MAP.md` for the broader semantics.

---

## Two modes

`phase_5_qc_triage` ships two drivers covering two distinct operating
modes. Pick the one that matches what you're doing:

### Mode 1 — single-candidate (`run_chrom.sh`)

For **exploration and debugging**: you have one chromosome with one
shelf region in mind, you want the full Q01→Q09 pipeline run on it
with explicit breakpoints, and you want all diagnostic plots emitted.

```bash
cd ${BASE}/inversion-popgen-toolkit/inversion_modules/phase_5_qc_triage
bash install.sh                                   # one-time: paths + compile

# Single candidate, explicit shelf + breakpoint coords:
SHELF_START_MB=15 SHELF_END_MB=18 \
BP1_MB=15.115  BP2_MB=18.005 \
  bash run_chrom.sh C_gar_LG28
```

Outputs land in `${QC_OUT}/C_gar_LG28/` with all nine Q-step
sub-directories. Use skip flags (`SKIP_Q01=1 SKIP_Q02=1 …`) to
re-run only specific steps after fixing something.

This is the mode to use when you're investigating a known shelf,
tuning thresholds, or reproducing a figure for one candidate.

### Mode 2 — genome-wide (`run_all_28chrom.sh --resume`)

For **production runs**: you want to sweep the whole genome,
registering every candidate in the shared results registry, with
resume semantics for incremental updates.

```bash
cd ${BASE}/inversion-popgen-toolkit/inversion_modules/phase_5_qc_triage

# 1) populate shelf_coords.tsv from phase_4_catalog's candidate list:
bash scripts/bridge_from_phase4.sh \
  ${PHASE4_OUT}/candidate_summary.tsv \
  --min-tier 3 --min-span-kb 50

# 2) submit the genome-wide SLURM array (Q01→Q10 per chrom, then stitch):
bash run_all_28chrom.sh --resume
```

`--resume` skips chromosomes that already have Q10 registry entries —
useful for re-runs after a pipeline update or a partial crash.
`--query-only` prints the current registry summary without computing.
`--no-array` runs sequentially on the login node (debug only — do
not use for production; ~10h wall for 28 chroms).

Outputs stitch into a global table. The bridge script in step 1 reads
from `phase_4_catalog/candidate_summary.tsv` and writes
`shelf_coords.tsv` which the SLURM array consumes — this is the
wiring between phase_4 and phase_5.

---

## Layout

```
inversion_modules/phase_5_qc_triage/
├── 00_config.sh                      # paths + tool bins, sourced by every script
├── install.sh                        # one-shot setup + Engine B/F compile
├── run_chrom.sh                      # ★ MODE 1 driver: single candidate (Q01→Q09)
├── run_all_28chrom.sh                # ★ MODE 2 driver: genome-wide (Q01→Q10 via SLURM array)
├── run_all.sh                        # legacy driver (Q01–Q08 only; use run_chrom.sh instead)
├── export_precomp_to_json.R          # build scrubber JSON (theta + ancestry)
│
├── STEP_Q01_snp_density.sh           # SNP density per bin from .pos.fixed
├── STEP_Q02_beagle_uncertainty.sh    # BEAGLE posterior uncertainty per bin
├── STEP_Q03_coverage_tracks.sh       # per-bin coverage from mosdepth
├── STEP_Q04_plot_diagnostic.sh       # stacked diagnostic PDF
├── STEP_Q05_aggregate_theta.sh       # pestPG → per-sample theta_pi long TSV
├── STEP_Q06_ancestry_tracks.sh       # Engine B local_Q → per-window ancestry
├── STEP_Q06_multiscale.sh            # Engine B at 1x/5x/10x window scales
├── STEP_Q06_precompute.sh            # Engine B instant_q precompute
├── STEP_Q07_popstats.sh              # Engine F (Fst/dXY/theta/Tajima) per group
├── STEP_Q07b_hobs_per_group.sh       # Hobs per karyotype group (MODULE_5E successor)
├── STEP_Q07c_hobs_windower.sh        # Hobs windowed summariser
├── STEP_Q08_shelf_heatmap.sh         # sample × SNP heatmap, 9-panel composite
├── STEP_Q09_gap_characterization.sh  # interior-structure gap characterization
├── STEP_Q09b_shelf_ld_check.sh       # LD sanity check across the shelf
├── STEP_Q10_register.sh              # write per-candidate summary.json → registry
│
├── R/                                # R workers called by STEP_Q* shell wrappers
│   ├── q02_beagle_stream.R
│   ├── q03_coverage_collapse.R
│   ├── q04_compose_plot.R
│   ├── q05_aggregate_theta.R
│   ├── q06_ancestry_tracks.R
│   └── q07_build_groups.R
├── scripts/
│   ├── bridge_from_phase4.sh         # populate shelf_coords.tsv from phase_4_catalog
│   ├── build_ref_n_bed.sh
│   └── registry_query.sh             # query per-candidate / genome-wide Q10 results
├── slurm/
│   └── array_28chrom.sh              # SLURM array wrapper dispatched by run_all_28chrom
├── registry_bridge/                  # writes to shared results registry
├── shelf_coords.tsv                  # per-chromosome shelf + breakpoint coords (mode 2 input)
└── docs/
    └── interpretation.md             # how to read the Q-step tracks
```

---

## What the Q-steps produce

| Step | Compute | Key output | What it tells you |
|---|---|---|---|
| Q01 | SNPs/window from .pos.fixed | SNP density track | low = mapping gap / low polymorphism |
| Q02 | BEAGLE posterior uncertainty | uncertain_frac track | high = caller lacked confidence |
| Q03 | mosdepth coverage | mean_cov + CV track | low + low CV = uniform mappability drop (repeat) |
| Q04 | composite plot | 11-panel stacked PDF | shelf vs ref summary + interpretation guide |
| Q05 | ANGSD pestPG aggregation | per-sample θπ (heterozygosity) | low + uniform = low polymorphism; low + stratified = real inversion |
| Q06 | Engine B local_Q (single + multi-scale) | per-window Δ12, ENA, maxQ | stable Δ12 = real structure; scale-dependent = noise |
| Q07 | Engine F region_popstats | Fst/dXY/θ/Tajima per group | Hom1-vs-Hom2 Fst spike = inversion signature |
| Q07b/c | Hobs per karyotype group | per-group Hobs track | MODULE_5E successor (Hobs confirmation) |
| Q08 | sample × SNP shelf heatmap (9 scales) | per-panel PDF + composite | 3-band pattern across all scales = robust; scale-dependent = noise |
| Q09 | interior-structure gap characterization | gap annotations | signals nested inversions / double crossovers / GC tracts |
| Q09b | LD sanity check across shelf | shelf LD tracks | confirms cohesive haplotype structure |
| Q10 | per-candidate summary.json | registry entry | **the artifact consumed by phase_9** |

`Q10` is what feeds `phase_9_classification/_qc_shelf_reader.R`. That's
where the `q_qc_shelf_*` keys get attached to the final candidate
table.

---

## The `q_qc_shelf_*` key family (wiring to phase_9)

Q10 writes `summary.json` per candidate with the fields the pass-13
reader expects. In `phase_9_classification/compute_candidate_status.R`,
setting `QC_SHELF_EVIDENCE_DIR` to the Q10 output directory appends
13 keys to every candidate row:

| Key | Semantic |
|---|---|
| `q_qc_shelf_ran` | logical — did phase_5 run for this cid |
| `q_qc_shelf_flag` | one of clean / low_snp / high_uncertain / coverage_artifact / messy / unknown |
| `q_qc_shelf_fst_enrichment_fold` | Fst(shelf)/Fst(flanking) — supporting biological signal |
| `q_qc_shelf_hovere_het_inside` | HoverE in heterokaryotypes — inversion signature |
| `q_qc_shelf_snp_ratio` | SNP density(shelf)/SNP density(flanking) |
| `q_qc_shelf_uncertain_ratio` | BEAGLE uncertainty(shelf)/BEAGLE uncertainty(flanking) |
| `q_qc_shelf_coverage_ratio` | coverage(shelf)/coverage(flanking) |
| … | (six more threshold-decision keys — see `compute_candidate_status.R` build_key_spec) |

All 13 keys are **soft flags**: they never cap tier, never block
characterization, never gate validation. They just add columns
readers can de-prioritize on.

Threshold overrides via env vars (`QC_SHELF_FST_THRESH`,
`QC_SHELF_HOVERE_THRESH`, etc.) so tuning doesn't require code edits.

---

## Engines (Q06, Q07)

**Engine B = `instant_q`**: C++ binary, per-window ancestry Q estimation
via fixed-F EM on BEAGLE GL. Outputs Δ12 / entropy / ENA per window
plus per-sample maxQ. Called by Q06 (single-scale), Q06-multi
(1x/5x/10x), and the precompute step (cache for downstream reads).

**Engine F = `region_popstats`**: C binary, per-window Hudson Fst, dXY,
dA, θ_π, θ_W, Tajima's D across 1+ sample groups. Reads BEAGLE,
runs at fixed windowed scale (default 50kb/10kb to match pestPG).

Q07 uses Engine F with two grouping strategies in separate runs:

1. **Ancestry groups** — samples bucketed by their dominant maxQ from
   Engine B. Per-group Fst answers "is the signal driven by ancestry
   structure?"
2. **Inversion genotypes (Hom1/Het/Hom2)** — k-means(3) on PC1
   loadings averaged over shelf windows, then Fst across the whole
   chromosome. Classical inversion signature: Fst spikes at the
   shelf and drops outside.

---

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

---

## Skip patterns (mode 1)

```bash
# Reuse everything, only regenerate the plot:
SKIP_Q01=1 SKIP_Q02=1 SKIP_Q03=1 SKIP_Q05=1 \
SKIP_Q06=1 SKIP_Q06_MULTI=1 SKIP_Q07=1 SKIP_Q08=1 SKIP_Q09=1 \
  SHELF_START_MB=15 SHELF_END_MB=18 \
  bash run_chrom.sh C_gar_LG28
```

---

## Dependencies

- R with data.table, ggplot2, patchwork
- mosdepth in PATH (only if `RUN_MOSDEPTH=1` in Q03)
- Existing ANGSD BEAGLE + .pos.fixed (from `phase_2_discovery/2a_local_pca`)
- Existing mosdepth output at `het_roh/mosdepth_output` (or set `MOSDEPTH_EXISTING`)
- Engine B (`instant_q`) and Engine F (`region_popstats`) compiled — handled by `install.sh`

---

## Interpretation

See `docs/interpretation.md` for how to read the Q-step tracks and
translate them into a `q_qc_shelf_flag` decision.

---

## Related docs

- `docs/MODULE_MAP.md` — overall phase map (see the phase_5 sub-section)
- `phase_9_classification/_qc_shelf_reader.R` — reads this phase's Q10 outputs
- `phase_9_classification/compute_candidate_status.R` — where `q_qc_shelf_*` keys land on the candidate table
- `CHANGELOG.md` — full history of this module (phase_qc_shelf → 4b_qc_triage → phase_5_qc_triage)
