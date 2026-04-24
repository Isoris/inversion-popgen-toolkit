# CHANGELOG — phase_5_qc_triage

(Module history: `phase_qc_shelf/` → `phase_4_postprocessing/4b_qc_triage/` (pass 12, 2026-04-24) → `phase_5_qc_triage/` (pass 15, 2026-04-24). This CHANGELOG retains entries from all three name eras.)

## Pass 11 task 2 finalization — 2026-04-24 (post pass 15)

- README rewritten to reflect the 2-mode structure: `run_chrom.sh` for single-candidate exploration, `run_all_28chrom.sh --resume` for genome-wide production. Legacy `run_all.sh` flagged as superseded.
- Added "`q_qc_shelf_*` key family" section documenting the 13 soft-flag keys appended to the final candidate catalog via pass 13's reader wiring (`phase_9_classification/_qc_shelf_reader.R`).
- Layout section expanded to show all 16 STEP_Q* scripts including Q06_precompute, Q07b, Q07c, Q09, Q09b, Q10 — previously missing.
- Wiring diagram (phase_4_catalog → bridge → shelf_coords.tsv → Q10 → phase_9) moved into `docs/MODULE_MAP.md` under a new phase_5_qc_triage sub-section so MODULE_MAP is discoverable on its own without jumping here first.

No code change. README + docs only.

## v3.1 — 2026-04-20 (Engine H / Merot Hobs integration)

Wires the pre-built Engine H (hobs_windower) and patched ANGSD
(Isoris/angsd_fixed_HWE) into the pipeline as Q07b + Q07c steps. Adds
group-stratified HWE deviation as a second axis of inversion evidence
alongside Fst.

### Rationale

At an inversion SNP inside the Het karyotype group, every Het sample
is heterozygous (one inverted + one non-inverted haplotype), so
Hobs_Het → 1, Hexp_Het → 0.5, **Hobs/Hexp → 2.0**. In the Hom1 and
Hom2 groups, samples are fixed at their own arrangement, so
Hobs/Hexp → 0. Outside the inversion, all three groups sit at ~1.0
(HWE-consistent). This three-way signature, scientifically validated
by Mérot et al., is orthogonal to Fst (which uses allele frequency
differences) and provides independent confirmation of inversion
polymorphism.

### New steps

- **`STEP_Q07b_hobs_per_group.sh`** — runs patched ANGSD with `-doHWE 1`
  separately for each karyotype group (Hom1 / Het / Hom2 BAM subsets),
  emitting per-group `.hwe.gz`. Three ANGSD calls launched in parallel
  per chromosome. Binary path: `$BASE/angsd_fixed_HWE/angsd`
  (fallback search in common locations).
- **`STEP_Q07c_hobs_windower.sh`** — feeds each per-group `.hwe.gz`
  into the existing `hobs_windower` C binary at 7 multi-scale sliding
  windows (5 kb / 10 kb / 50 kb / 100 kb / 250 kb / 500 kb / 1 Mb).
  Merges the three per-group windowed outputs into a single wide TSV:
  `hobs_merged.<CHR>.tsv.gz`.
- **`R/q07c_merge_hobs.R`** — R worker that merges hobs_windower output
  into the wide format expected by Q04 and Q10.

### Q04 upgrade

- New panel p10 (bottom of diagnostic stack): three lines — HoverE_HOM1
  (blue), HoverE_HET (red), HoverE_HOM2 (green) — with reference lines
  at y=1 (HWE) and y=2 (expected Het excess at inversion).
  At a real inversion, HET line peaks toward 2 inside the shelf, both
  HOM lines drop toward 0 inside the shelf, all three sit near 1
  elsewhere. The cleanest inversion visual available.

### Q10 upgrade

- Engine H outputs registered into `results_registry`:
  - `kind=group_hwe stat=hobs_hexp who=hom1_het_hom2 where=<cid>` —
    full HoverE profile as a pairwise-style file
  - `kind=interval_summary stat=hobs_summary where=<cid>` — scalar
    summary stats (HoverE_Hom1/Het/Hom2 means inside vs outside shelf)

### Runtime cost

Per chromosome on 4 cpus:
- Q07b ANGSD (3 groups parallel, 2 cpus each): ~15 min
- Q07c hobs_windower: ~30 sec
- Q04 panel rendering: +2 sec
Total added to pipeline: ~15 min per chromosome, or ~7 hrs across
28 chromosomes on a SLURM array with 8 concurrent tasks.

### Not yet tested on LANTA

Smoke test before launching the 28-chrom sweep:
```bash
bash STEP_Q07b_hobs_per_group.sh C_gar_LG28
bash STEP_Q07c_hobs_windower.sh C_gar_LG28
ls -la ${QC_TRACKS}/hobs_merged.C_gar_LG28.tsv.gz
# re-render Q04 to see the new panel
bash STEP_Q04_plot_diagnostic.sh C_gar_LG28
```

## v3.0 — 2026-04-20 (registry-wired)

Major milestone: phase_qc_shelf is no longer a standalone QC viewer. It now
**populates all four registries** of the inversion-popgen-toolkit, so every
inversion discovered becomes a queryable database row rather than a PDF on
disk.

### New step: Q10 registry bridge

- **`STEP_Q10_register.sh` + `R/q10_register.R`**: for each chromosome, reads
  phase_qc_shelf outputs and registers them via the existing registry APIs:
  - `interval_registry.add(cid, chrom, start_bp, end_bp, method="phase_qc_shelf")`
  - `sample_registry.add_group("inv_<cid>_HOM1"/HET/HOM2, ids, dimension="karyotype", parent="all_226")`
  - `results_registry.put_pairwise(kind, stat="fst", who="hom1_vs_hom2", where=cid, ...)`
  - `results_registry.put_interval_summary(cid, stat="fst_summary", values=...)`
  - `evidence_registry.add_candidate(cid, files=[pdfs + tsvs + summary.json])`
- Idempotent: re-running with the same inputs updates in place. Use
  `FORCE_OVERWRITE=1` to replace existing cids.
- Provenance: every registry row records script name, timestamp, and sha256
  hashes of the precomp + invgt + popstats inputs.
- ALL mode: `bash STEP_Q10_register.sh ALL` walks every chromosome that has
  phase_qc_shelf outputs.

### Query CLI: registry_query.sh

- **`scripts/registry_query.sh` + `R/registry_query.R`**: command-line
  interrogation of the registries. Commands:
  - `summary` — counts by chrom / kind / dimension
  - `list_candidates [--chrom X] [--method Y]` — all cids with coords
  - `describe <cid>` — interval + karyotype groups + results + evidence
  - `karyotypes <cid>` — print Hom1/Het/Hom2 sample IDs
  - `fst <cid>` — show Fst profile and summary stats

### 28-chromosome orchestration

- **`slurm/array_28chrom.sh`**: SLURM array driver, one task per chrom,
  throttled at 8 concurrent. 4 CPUs, 16 GB, 4 hr walltime per chrom.
  Reads optional `shelf_coords.tsv` for per-chromosome shelf/BP coords.
- **`run_all_28chrom.sh`**: end-to-end orchestrator. Three modes:
  - default: submit SLURM array and exit
  - `--no-array`: run sequentially on login node (debugging)
  - `--query-only`: skip compute, just print registry summary
  - `--resume`: skip chromosomes already in registry
- **`shelf_coords.tsv`**: seed file with LG28 coordinates from this session.
  New inversions found get added here as discovered.

### Integration with existing registries

Pre-existing (from chat-16):
- `registries/data/sample_registry/` — append-only, auto-backup
- `registries/data/interval_registry/` — UUID row_ids, FK enforcement
- `registries/data/evidence_registry/` — per-candidate file bundles
- `registries/data/results_registry/` — pairwise/summary results with
  `ask(kind, stat, who, where)` query plane

phase_qc_shelf writes into all four through their existing APIs. No
schema changes required; the integration uses documented public methods.

### Workflow from "fresh" to "all inversions databased"

```bash
# one time:
cd /scratch/.../phase_qc_shelf
bash scripts/build_ref_n_bed.sh ${REFERENCE_FASTA}  # dropout stipple source
bash install.sh

# once per genome:
sbatch slurm/array_28chrom.sh          # launches 28 tasks, ~2-8 hrs wall
# ... wait for completion ...
bash STEP_Q09_gap_characterization.sh ALL
bash STEP_Q10_register.sh ALL          # idempotent catch-any-missed
bash scripts/registry_query.sh summary
```

### What the registry enables

Any downstream question about inversions is now a one-liner:

```bash
# Every inversion in the cohort
bash scripts/registry_query.sh list_candidates

# All inversions on a specific chromosome
bash scripts/registry_query.sh list_candidates --chrom C_gar_LG01

# Everything about a specific candidate
bash scripts/registry_query.sh describe inv_C_gar_LG28_15_18

# Fst profile for an inversion
bash scripts/registry_query.sh fst inv_C_gar_LG28_15_18

# Sample IDs in each karyotype
bash scripts/registry_query.sh karyotypes inv_C_gar_LG28_15_18
```

And from any R or Python script, the same data is retrievable via
`reg$results$ask(...)`, `reg$samples$get_subgroups(cid)`,
`reg$evidence$get_candidate(cid)`. One source of truth.

## v2.1.1 — 2026-04-20 (audit round)

Audit of v2.1 before wider deployment caught and fixed several integration
gaps. All fixes are mechanical — no scientific reinterpretation required.

### Audit fixes

- **Q08 heatmap breakpoint markers** (`R/q08_shelf_heatmap.R` +
  `STEP_Q08_shelf_heatmap.sh`): added `--breakpoint1_mb` / `--breakpoint2_mb`
  args and `ComplexHeatmap::anno_mark` layer on the top column annotation.
  The Q08 heatmap now labels the exact BP positions above the columns.
- **run_chrom.sh env var naming**: was setting `BP1_MB`/`BP2_MB` but the Q04
  driver reads `BREAKPOINT1_MB`/`BREAKPOINT2_MB`. Fixed the translation.
  Also added `SMOOTH_WIN` / `SNP_DENSITY_SCALE_KB` pass-through.
- **run_chrom.sh Q08 call**: now passes `BP1_MB`/`BP2_MB` through.
- **Q06 multiscale cache discovery** (`STEP_Q06_multiscale.sh`): `resolve_summary`
  and `resolve_samples` now also look at `scale_dense/K<NN>/...` and
  `scale_thin/K<NN>/...` layouts emitted by the v3 multi-K precompute.
- **Q07 popstats cache discovery** (`STEP_Q07_popstats.sh`): was only
  looking in flat and `scale_1x/`. Now prefers `scale_dense` (higher SNP
  resolution → better invgt grouping), then `scale_thin`, then legacy paths.

### Not changed

- The Engine B precompute's symlink-to-canonical-K behavior means K=8
  output at `scale_<res>/K08/` is mirrored as `scale_<res>/...` flat files,
  so `STEP_Q06_ancestry_tracks.sh` finds them via existing lookup. No patch
  needed there, confirmed by audit.

## v2.1 — 2026-04-20 (LG28 follow-up: dropout semantics)

Second wave of fixes after the v2.0 LG28 run uncovered a 360 kb reference
dropout inside the inversion. Adds explicit "no-data" semantics so readers
can't mistake missing data for biology, plus a gap characterization step
and a one-shot per-chromosome wrapper.

### New features

- **No-data background strips** in Q04 panels (`detect_nodata_regions()`
  in `R/q04_compose_plot.R`). Light-gray rectangles render behind every
  panel at contiguous runs of zero-SNP windows, so the reader sees "data
  is absent here" distinct from "value is low here".
- **Two-variant Q04 output**: `--no_nodata_strips` flag renders a clean
  version without the gray bands. Driven by env var `Q04_NO_STRIPS=1`
  and filename suffix `OUT_SUFFIX=.clean` in the Q04 driver.
- **Two-layer ideogram stipple**:
  - Light-gray dot pattern over zero-SNP regions (from precomp span_kb)
  - Dark-navy dot pattern over reference-N regions (from optional
    `--ref_n_bed` / `REF_N_BED` BED file)
  - Staggered printer-style dot grid via a new `make_stipple()` helper
  - Ideogram subtitle auto-reports the number of dropout regions
- **`STEP_Q09_gap_characterization.sh`** (new step). For each chromosome,
  emits per-gap TSV and BED with annotations:
  - `frac_N`, `frac_softmasked`, `frac_uppercase` (reference composition)
  - `frac_repeat_bed` (overlap with a RepeatMasker-style BED if provided)
  - `cov_mean`, `cov_n_low` (from mosdepth if provided)
  - `flag`: ASSEMBLY_GAP / REPEAT_RICH / LOW_MAPPABILITY / DROPOUT
  - `ALL` mode stitches per-chrom tables into genome-wide files
- **`scripts/build_ref_n_bed.sh`** utility: one-shot awk scan of a FASTA
  to extract all N-runs >= MIN_GAP_BP into per-chrom BED files. Feeds
  the ideogram dark-stipple layer and Q09 reference annotation.
- **`run_chrom.sh`** wrapper: runs Q01→Q09 for a single chromosome
  including both Q04 variants (with and without strips). Supports
  per-step `SKIP_*` env vars. Prints a summary at the end listing all
  PDFs and gap tables produced.

### Fixes

- **Z candidate list** (`R/q04_compose_plot.R`): added `max_abs_z` as the
  primary candidate and `MDS1_z` as tertiary fallback, covering the two
  column names Engine B writes in different versions of the precomp.
- **span_kb** now computed immediately after precomp load so `detect_nodata`
  can see it on the first `decorate()` call (previously it was computed
  inside the density panel and only appeared after the first few panels
  had already rendered without the strip).

### Scientific finding (LG28, this round)

- Confirmed: **one inversion**, 15.115 – 18.005 Mb, span 2.88 Mb
- Internal **360 kb dropout** at 15.335 – 15.695 Mb: 38 consecutive windows
  with n_sites=0 — hard reference/mappability gap inside the inversion
- Fst on both sides of the gap matches (0.40 at 15.32 → 0.44 at 15.70),
  confirming single inversion with an unmappable hole, NOT sub-inversions
- Shelf dropout is now rendered as a clean gap rather than a false
  drop-to-zero in the Fst panel

## v2.0 — 2026-04-20 (LG28 field-tested)

All fixes below were diagnosed and applied during a single debugging session
that took the module from broken on LANTA to a publication-quality figure
showing a balanced 2.89 Mb inversion on LG28.

### Scientific validation

- **LG28 shelf (15.115–18.005 Mb) confirmed as real segregating inversion**
  - Fst_Hom1_Hom2 = 0.264 at shelf vs 0.011 genome-wide (24× elevation)
  - 60/106/60 group sizes, Hardy-Weinberg consistent (p=0.5)
  - Flat-Z signature (shelf SD = 0.14 vs 0.97 ref) — novel detection regime
    for inversions missed by peak-Z thresholding
  - Asymmetric π between arrangements: π_Hom1=0.27, π_Hom2=0.35

### Breaking changes

None. v2.0 is fully backward-compatible with existing caches and track files.

### Fixes — bash drivers

- **`install.sh`**: fixed set-e bug where `[[ ... ]] && missing=$((missing+1))`
  would exit the script on false-branch when set -e was active.
- **`STEP_Q03_coverage_tracks.sh`**: soft-skip when `MOSDEPTH_EXISTING` dir
  doesn't exist (was `qc_die` hard-kill, now `qc_log` + `return 0`).
- **`STEP_Q06_precompute.sh`**: preserve `INSTANT_Q_BIN` and `LOCAL_Q_DIR`
  before sourcing `00_ancestry_config.sh`. That config clobbers both to
  paths that match its own assumed layout (`${BASE}/unified_ancestry/…`)
  rather than the actual repo-local layout.
- **`STEP_Q06_precompute.sh` v2 (multi-K × multi-scale)**: emits caches at
  `scale_<res>/K<NN>/<CHR>.local_Q_*.tsv.gz`. Supports `Q06_K_SWEEP` env var
  (e.g. "4,8,10,14"). K-equal-to-canonical gets symlinked to the flat
  `scale_<res>/` path for backward-compatible Q06 readers.

### Fixes — R workers

- **`R/q06_ancestry_tracks.R`** (full rewrite): uses Engine B's native column
  names (`mean_delta12`, `mean_entropy`, `mean_ena`, `start_bp`, `end_bp`,
  `window_id`). Handles the per-sample TSV's `sample_id`, `max_q`,
  `assigned_pop`, with fallback to argmax over Q1..QK for `maxQ_label`
  derivation. Previous version erroneously expected pre-renamed columns
  and died with "Summary missing column: delta12".
- **`R/q04_compose_plot.R`** (four patches):
  1. Z-column candidate list extended: `max_abs_z` is now the default read
     from the precomp RDS; `MDS1_z` is a secondary fallback. Previous list
     missed every column name actually present in the precomp.
  2. Placeholder panels are now `NULL` instead of empty `theme_void`
     plots, and the final assembly filters NULL entries from the layout.
     Previously missing-data panels drew as tall blank gaps in the PDF.
  3. SNP density panel is configurable via `--snp_density_scale_kb`
     (default 10). Previous panel plotted `n_snps` which is constant at
     the window-size config (e.g. always 100) and uninformative.
  4. Popstats panels (θπ per invgt, Fst Hom1_Hom2) filter rows with
     `n_sites_used < Q04_MIN_SITES` (default 20) before plotting. Engine F
     writes exact zeros for under-sampled windows; those were previously
     rendered as fake drops-to-baseline that looked like false boundaries.
- **`R/q08_shelf_heatmap.R`**: fixed `subscript out of bounds` in BEAGLE
  streaming. The named-integer `target_lookup[[key]]` pattern throws rather
  than returning NULL when the key is absent. Replaced with a hashed
  environment for proper O(1) "missing key returns NULL" semantics.

### Fixes — config wiring

- New `config.local.sh` anchors all paths cleanly in one file (BASE,
  UNIFIED_ANCESTRY_DIR, LOCAL_Q_DIR, INSTANT_Q_BIN, BEAGLE_DIR,
  REGION_POPSTATS_BIN, SAMPLE_LIST, SAMPLE_LIST_POPSTATS). Translates
  names from the shared `00_inversion_config.sh` (HETDIR→HET_DIR,
  SAMPLES_IND→SAMPLE_LIST) with fallback defaults.

### New capability

- **Multi-K × multi-scale Engine B precompute**. Default: canonical K=8 at
  thin+dense BEAGLE scales (4 total runs per chromosome). With
  `Q06_K_SWEEP="2,4,8,11"` you get 8 runs per chromosome in ~30 seconds.
  Feeds Q06-multiscale panel for ancestry robustness check.

### Tested

- LG28 end-to-end: 10 Q01 bins → 4302 Engine B windows → 2011 popstats rows
  → 9-panel Q08 heatmap + 10-panel Q04 diagnostic PDF. Total wall time
  ~15 minutes on LANTA login node (excluding mosdepth Q03).

## v1.0 — initial release

Module structure, Q01-Q08 step drivers, R workers, install.sh, run_all.sh.
