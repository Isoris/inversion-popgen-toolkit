# utils/lib_ghsl_panel.R — GHSL v6 per-sample panel query API

**Created:** chat 14 (2026-04-18).
**Backend for:** `reg$compute$ghsl_*`, `lib_ghsl_confirmation.R`'s
panel-based helper, figure scripts, ad-hoc queries.

## What problem this solves

The v6 classifier (`phase_2_discovery/2e_ghsl/STEP_C04b_snake3_ghsl_classify.R`)
emits one dense per-sample panel RDS per chromosome containing:

- per-sample rolling GHSL divergence at every computed scale
  (default `s10, s20, s30, s40, s50, s100`)
- per-sample cohort-relative rank at every scale
- per-sample rank-band {`LOW`, `MID`, `HIGH`} at every scale
- per-sample stable-run membership from Part B karyotype calling
- window-level `ghsl_v6_score`, `ghsl_v6_status`, `rank_stability`,
  `div_contrast_z`, `div_bimodal`

This library is the query layer. Consumers don't read the RDS directly;
they call these functions and get back slices / aggregates / wide
matrices / plots with memoized in-session caching.

## Quick reference

```r
source("utils/lib_ghsl_panel.R")
panel_set_default_dir("/scratch/.../ghsl_v6_out")

# dense slice: one row per (sample, window) in [start, end]
ghsl_panel_range("C_gar_LG12", 10e6, 14e6)

# per-sample summary across range (one row per sample)
ghsl_panel_aggregate("C_gar_LG12", 10e6, 14e6,
                      summaries = c("mean", "rank_mean",
                                    "frac_low", "frac_high",
                                    "longest_low_run_bp",
                                    "longest_high_run_bp"))

# multi-subblock scan (e.g. 6 sub-blocks from 4 soft boundaries)
subs <- data.table(
  subblock_id = 1:6,
  start_bp = c(10e6, 10.5e6, 11e6, 11.7e6, 12.4e6, 13.2e6),
  end_bp   = c(10.5e6, 11e6, 11.7e6, 12.4e6, 13.2e6, 14e6))
ghsl_panel_subblock_scan("C_gar_LG12", subs)

# wide matrix for heatmap plotting: [samples × windows]
ghsl_panel_wide_matrix("C_gar_LG12", 10e6, 14e6,
                         metric = "div_roll_s50")

# plot helpers
plot_ghsl_sample_track("C_gar_LG12", "CGA_042",
                         scales = c("s10", "s50"))
plot_ghsl_heatmap_track("C_gar_LG12", 10e6, 14e6,
                          metric = "rank_band_s50")
```

## Via the registry (preferred in pipeline code)

```r
reg <- load_registry()

reg$compute$ghsl_at_candidate("LG12_17")
reg$compute$ghsl_at_interval("C_gar_LG12", 10e6, 14e6,
                               sample_subset = reg$samples$get_carriers(cid, "HET"))
reg$compute$ghsl_at_subblocks("LG12_17", subs)
reg$compute$ghsl_at_block_subregions("LG12_17", "regime_segments")
reg$compute$ghsl_wide_matrix("LG12_17", metric = "div_roll_s20")
```

## Scale ladder

Heavy-engine default (chat 14):

| key   | rolling window (base windows × 5 kb) | use                                             |
|-------|---------------------------------------|--------------------------------------------------|
| s10   | ~50 kb                                | fine detail; sub-inversion boundaries, short GC tracts |
| s20   | ~100 kb                               | short-inversion bodies, sub-block distinction    |
| s30   | ~150 kb                               | intermediate                                     |
| s40   | ~200 kb                               | intermediate                                     |
| s50   | ~250 kb                               | **default** — karyotype calling, scoring, typical inversion-body display |
| s100  | ~500 kb                               | chromosome-overview plotting                     |

Scale-agreement across adjacent scales is itself signal: a sample whose
rank-band is stable across s10→s50 is a confident call; a sample whose
band only appears at s50 and vanishes at s10 is borderline or near a
boundary.

Choosing a non-default scale:
- Give `scale = "s20"` or `scale = 20` (same thing; resolution helper
  accepts int / "20" / "s20").
- Score/status columns (`ghsl_v6_score`, `ghsl_v6_status`) are computed
  at the classifier's `--scale` (default s50) at emit time. They do not
  change if you pick a different scale for aggregation — the scale
  argument only affects `div_roll_<s>` / `rank_in_cohort_<s>` /
  `rank_band_<s>` selection.

## Function reference

### Data access

- `panel_set_default_dir(ghsl_dir)` — globally set GHSL output dir.
- `panel_get_default_dir()` — read current default (falls back to
  `GHSL_DIR` env var).
- `panel_clear_cache(chrom = NULL)` — invalidate the in-session panel
  cache (all chroms or one).
- `load_ghsl_panel(chrom, ghsl_dir = NULL, refresh = FALSE)` — returns
  the full per-chrom data.table; memoized.
- `ghsl_panel_meta(chrom, ghsl_dir = NULL)` — the panel's metadata
  attribute (sample order, window info, scales available, thresholds,
  timestamp, classifier version string).

### Queries

- `ghsl_panel_range(chrom, start_bp, end_bp, sample_ids, ghsl_dir)` —
  dense slice, one row per (sample, window).
- `ghsl_panel_aggregate(chrom, start_bp, end_bp, sample_ids, scale,
  summaries, ghsl_dir)` — one row per sample. Summaries:
  `mean`, `median`, `sd`, `frac_low`, `frac_mid`, `frac_high`,
  `longest_low_run_bp`, `longest_high_run_bp`, `rank_mean`,
  `rank_at_peak`, `stable_run_call`, `n_windows`.
- `ghsl_panel_subblock_scan(chrom, subblocks, sample_ids, scale,
  summaries, include_transitions, ghsl_dir)` — one row per
  (sample, subblock). With `include_transitions = TRUE` (default)
  attaches a `per_sample_patterns` attribute: one row per sample
  with concatenated `band_pattern` (e.g. `"LOW/LOW/LOW/HIGH/HIGH/HIGH"`)
  and `n_transitions` count. A sample with `n_transitions ≥ 1` is the
  recombinant-detection signal at the sub-block level.
- `ghsl_panel_wide_matrix(chrom, start_bp, end_bp, sample_ids, metric,
  sample_order, ghsl_dir)` — `[samples × windows]` matrix ready for
  `image()`, `pheatmap`, `ComplexHeatmap`, etc. Categorical metrics
  (`rank_band_<s>`, `stable_run_call`) come back as character;
  numeric metrics come back numeric.

### Plot helpers (minimal; grab the underlying data and compose with
`figrid` / `ggplot2` for publication figures)

- `plot_ghsl_sample_track(chrom, sample_id, scales, start_bp, end_bp,
  ghsl_dir, out_png)` — multi-scale divergence line track for one sample.
- `plot_ghsl_heatmap_track(chrom, start_bp, end_bp, metric,
  sample_order, sample_ids, ghsl_dir, out_png)` — cohort heatmap; good
  for the "yellow-track" style figure panel.
- `plot_ghsl_subblock_panel(chrom, subblocks, sample_ids, scale,
  ghsl_dir, out_png)` — sample × subblock rank-band grid.

## Results-aware registry integration

`reg$compute$ghsl_at_block_subregions(cid, block_type)` reads an
existing evidence block for a candidate and queries GHSL across each
declared sub-region. Tries these common sub-region field names in
order: `segments`, `sub_blocks`, `subblocks`, `regions`, `tracts`,
`windows`, `soft_boundaries_blocks`. Accepts either `id` or
`subblock_id`, and start-field aliases `start` / `bp_start` /
`begin_bp` / `start_bp` (same for end).

Supported block types where this is known to work without extra
plumbing:
- `boundary` — if the block exposes left+right + soft-boundaries list
- `regime_segments` — uses the `segments` field from
  `STEP_C01j_regime_compatibility_engine.R`
- `recombinant_map` — uses the tract list from `STEP_C01i_b_multi_recomb.R`

If a block stores its sub-regions under a different key, either extend
the field-name search list in `ghsl_at_block_subregions` or wrap an
ad-hoc caller that builds a `subblocks` data.table explicitly and
calls `ghsl_at_subblocks` directly.

## Caching and memory

`load_ghsl_panel` caches the full per-chrom data.table in an env-level
cache keyed by chrom. A typical chrom panel is 30–80 MB compressed on
disk and 300–500 MB in memory once deserialized; keep this in mind for
scripts that touch many chroms in sequence. Use `panel_clear_cache(chrom)`
to free a specific chrom after you're done with it.

## File layout expected

```
<ghsl_dir>/
  snake3v6_window_track.tsv.gz
  snake3v6_karyotype_calls.tsv.gz
  snake3v6_summary.tsv
  annot/
    C_gar_LG01.ghsl_v6.annot.rds
    C_gar_LG01.ghsl_v6.karyotypes.rds
    ...
  per_sample/
    C_gar_LG01.ghsl_v6.per_sample.rds
    ...
```

Files can also live directly under `<ghsl_dir>` without the `annot/`
or `per_sample/` subdirectories; the loader falls back automatically.

## When to use what

| Task | Function |
|---|---|
| Candidate-level per-sample interval summary | `reg$compute$ghsl_at_candidate()` |
| Free-form bp range query (e.g. a recombinant tract) | `reg$compute$ghsl_at_interval()` |
| Complex candidate with soft boundaries → sub-block comparison | `reg$compute$ghsl_at_subblocks()` or `ghsl_at_block_subregions()` |
| Figure yellow-track / cohort heatmap | `reg$compute$ghsl_wide_matrix()` then plot |
| Exploratory — all windows in a region | `ghsl_panel_range()` |
| Tier-3 SPLIT detection (within-sample karyotype change) | **Not this library.** Use `lib_ghsl_confirmation.R::ghsl_per_sample_in_interval()`. Panel collapses SPLIT. |
| Tier-3 alternative (interval class summary) | `lib_ghsl_confirmation.R::ghsl_per_sample_panel_in_interval()` (delegates here) |

## Known caveats

- **SPLIT is invisible at aggregate level.** A sample whose rank-band
  toggles LOW→HIGH→LOW inside a single interval averages to MID and
  looks uninformative in `ghsl_panel_aggregate`. Use the
  per-stable-run logic in `lib_ghsl_confirmation.R` for that.
- **Rank_in_cohort is recomputed at query time for each scale** in the
  emit block; the panel stores it precomputed, but the cohort is the
  cohort that was live when the classifier ran. If you restrict
  `sample_ids` at query time, rank values are still relative to the
  full cohort, not the subset. If you need within-subset ranks,
  compute them explicitly on `ghsl_panel_range` output.
- **Short intervals (< ~10 base windows ≈ 50 kb)** give very few
  rolling cells at s50+ — aggregation is noisy. Drop to s10 or s20
  for short intervals, or fall back to the window-level `metrics`
  from `annot/*.ghsl_v6.annot.rds`.

## Chat-14 acknowledgements

Scale ladder choice (10/20/30/40/50/100) and the "many tracks as plots"
framing came from the chat-14 design discussion. The library shape
(backend + registry wrapper, same pattern as `pairwise_stat` for FST)
matches the existing toolkit idiom.
