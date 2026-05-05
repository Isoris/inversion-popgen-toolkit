# phase_X_ushape_evolution

**Inversion evolutionary-shape classification layer for the Inversion Atlas.**

Quentin Andres · MS_Inversions_North_african_catfish · LANTA `lt200308`
Codebase v8.5 · Phase X · v0.1 (initial scaffold)

---

## What this module does

For each inversion candidate, the module computes per-window between-arrangement
divergence (dXY, FST), within-arrangement diversity (π_HOMO1, π_HOMO2), and
allele-frequency contrast across the inversion plus flanking regions. It then
summarises the profile **shape** (U vs internal peak vs flat-deep vs young vs
asymmetric) and assigns one of seven shape classes — gated by an oldness check
so that young, weakly-diverged inversions don't get mis-classified.

The output is a JSON layer (`ushape_evolution_v1.json`) loadable by the atlas
in two ways:

- **Server-primary (default):** the popstats server has a new endpoint
  `POST /api/ushape/candidate` that computes one candidate on demand. The
  atlas calls it when the user opens a candidate page; the result is cached
  client-side so opening the same candidate again is instant.
- **Static-fallback:** running `LAUNCH_U00_run_all.sh` produces a single JSON
  for the entire candidate cohort. The atlas can load this at startup and
  treat the file as a frozen snapshot (offline / no-server demos).

Both paths use the same R library (`R/ushape_lib.R` + `ushape_classify.R` +
`ushape_io.R`) and produce identical candidate blocks.

## What this module does NOT do

- It is **not** an adaptation test. It classifies shape, not cause.
  See `docs/ushape_evolution_class_rules.md` for the safe-language section.
- It does **not** include any TE / repeat / SD / mappability content.
  Breakpoint architecture is the job of `phase_8_comparative_breakpoint_fragility/`.
  Mixing the two would conflate "is the breakpoint fragile?" with "what does
  the divergence profile look like?" — they answer different questions.
- It does **not** validate karyotype groups. The module trusts the input
  HOMO_1 / HET / HOMO_2 assignments coming from the upstream PCA pipeline.
  Wrong groups → wrong shape; the upstream PCA is the source of truth.

## Inputs

The CLI/batch path needs three files plus a dosage directory; the server path
takes the same information inline in a POST body.

| input | shape | who consumes it |
|---|---|---|
| `candidates.tsv`        | `candidate_id, chrom, start_bp, end_bp, left_breakpoint, right_breakpoint`              | U01 + U06 |
| `candidate_groups.tsv`  | `candidate_id, sample_id, karyotype_group` (HOMO_1 / HET / HOMO_2)                      | U01       |
| `<chrom>.dosage.tsv.gz` | `chrom, pos, sample1, ...` with 0/1/2/NA dosages                                        | U01       |
| `config_ushape.yaml`    | thresholds + windowing                                                                  | all steps |

Matched-background z-scores (U03 / atlas inline) use **local-flank
resampling** — bootstrap `inside_n` windows from the union of left+right
flank windows, recompute the inside summary statistic per resample, z-score
the observed value against the resampled null. No genomewide windows are
required. The same helper (`flank_resample_zscores` in `R/ushape_lib.R`) is
called by both the offline batch and the popstats-server endpoint, so the
numbers in `output/ushape_evolution_v1.json` match what the atlas shows
when fetching dynamically.

The honest interpretation: this captures **local** variance correctly
(better than a naive "ratio against flank mean" which treats the flank
mean as a zero-variance point estimate). It does **not** establish a
chromosome-wide null. The JSON exposes this as
`matched_background.bg_mode = "flank_resample"` so downstream code can
distinguish it from a true genomewide null if you ever add one.

## Outputs

```
output/
├── window_stats/<candidate_id>.window_stats.tsv         # U01
├── candidate_raw_summary.tsv                             # U01
├── candidate_shape_scores.tsv                            # U02 (+z from U03)
├── candidate_shape_classes.tsv                           # U02
├── classification_summary.tsv                            # U02
├── candidate_matched_background.tsv                      # U03
├── shape_feature_matrix.tsv                              # U04
├── shape_pca_coordinates.tsv                             # U04
├── shape_cluster_assignments.tsv                         # U04
├── shape_cluster_summary.tsv                             # U04
├── shape_cluster_feature_means.tsv                       # U04
├── permanova_results.tsv / anosim_results.tsv / ...      # U05
├── validation_summary.txt                                # U05
├── ushape_evolution_v1.json                              # U06  ← atlas loads this
└── plots/<candidate_id>.{pdf,png}, all_candidates.pdf   # U07
```

## Pipeline DAG

```
                   candidates.tsv          candidate_groups.tsv
                          │                         │
                          ▼                         ▼
                ┌────────────────────────────────────────┐
                │ U01: compute_window_stats              │
                │   (offline path) reads dosage TSVs     │
                │   (server path)  reuses region_popstats│
                │                  output via R lib      │
                └────────────────────────────────────────┘
                          │
                          ▼
            window_stats/<cid>.window_stats.tsv
                  +  candidate_raw_summary.tsv
                          │
                          ▼
              ┌────────────────────────────┐
              │ U02: score + classify      │  age-gated 2-step rules
              └────────────────────────────┘
                          │
                          ▼
              candidate_shape_scores.tsv
              candidate_shape_classes.tsv
                          │
              ┌───────────┴───────────┐
              ▼                       ▼
     U03 matched-bg z-scores   U04 PCA + Ward + PAM
              │                       │
              └───────────┬───────────┘
                          ▼
                    U05 PERMANOVA / ANOSIM
                          │
                          ▼
                    U06 JSON exporter
                          │
                          ▼
                ushape_evolution_v1.json
                          │
                          ▼
                    U07 per-candidate plots
```

## Server endpoint (primary entry)

The popstats server gets one new endpoint. To wire it in, paste this into
`popstats_server.py` next to the existing `/api/popstats/groupwise`:

```python
from ushape_endpoint import UshapeReq, ushape_candidate as _ushape_handler

@app.post("/api/ushape/candidate")
async def ushape_candidate(req: UshapeReq) -> Response:
    return await _ushape_handler(
        req,
        ensure_ready=_ensure_ready, samples=SAMPLES, cache=CACHE,
        engines=ENGINES, run_region_popstats=run_region_popstats,
        resolve_beagle=resolve_beagle, beagle_dir=Path(CFG["beagle_dir"]),
        ushape_r_dir=Path(CFG.get("ushape_r_dir", "R")),
    )
```

Add to `popstats_server.config.yaml`:

```yaml
ushape_r_dir: /scratch/lt200308-agbsci/.../phase_X_ushape_evolution/R
```

The endpoint reuses `region_popstats` for FST/dXY/π and runs the R lib via
`Rscript --vanilla` once per request. Results are content-addressable cached
under `cache/ushape/<sha256>.json` so repeat opens of the same candidate are
free. Cache invalidates automatically when the engine binary or the R lib
files change.

## Atlas integration (client side)

`atlas_js/ushape_renderer.js` is a single drop-in script. Add **one** `<script>`
tag near the bottom of `Inversion_atlas.html`, after the popstats client lib
loads but before the candidate-router fires:

```html
<script src="js/ushape_renderer.js"></script>
<script>
  UShapeRenderer.init({
    server_url:        window.atlasServer && atlasServer.url,
    static_layer_url:  "json/ushape_evolution_v1.json",  // optional fallback
  });
</script>
```

Then add these container elements wherever you want the U-shape surfaces to
appear (the renderer no-ops if a container is absent — safe to deploy
incrementally):

| where | element id |
|---|---|
| candidate page card               | `<div id="ushape_card"></div>` |
| page 6 popstats track group       | `<div id="ushape_popstats_panel"></div>` |
| page 4 axis chip                  | `<div id="ushape_axis_chip"></div>` |
| page 17 stats table               | `<div id="ushape_stats_table"></div>` |
| candidate gallery filter dropdown | `<select id="ushape_filter_select"></select>` |

When the user opens a candidate, call:

```js
// `groups` is already in the atlas as `locked_labels` for the candidate.
// `interval` carries chrom + breakpoints from the candidate manifest.
UShapeRenderer.fetchForCandidate(candidateId, groups, interval)
  .then(() => UShapeRenderer.renderAll(candidateId));
```

This single call hits the server, caches the result, then renders the card,
the SVG popstats panel, the axis chip, the stats table, and refreshes the
gallery filter counts. The classification pill carries a `data-tip` attribute
that pops a tooltip on hover with the reason, key scores, sample sizes, and
flags.

## Class definitions (one-line)

| class | one-liner |
|---|---|
| `neutral_like_U_shape`                 | dXY higher at edges than center — "suspension bridge" |
| `locally_adapted_internal_peak_like`   | dXY/FST internal peak dominates edges                  |
| `locally_adapted_breakpoint_like`      | strong edge signal without center peak                 |
| `flat_high_deep_structural_haplotype`  | broadly elevated dXY, low flatness                      |
| `young_weak_divergence`                | below age gate — shape uninformative                    |
| `asymmetric_edge`                      | one edge much hotter than the other                     |
| `complex_mixed`                        | multiple shape signals firing at once                   |
| `insufficient_data`                    | too few inside windows or non-finite scores             |

Full plain-language rules: `docs/ushape_evolution_class_rules.md`.

## Important caution

Reviewer 2 will catch overclaiming. Use these phrasings:

- "candidates show a U-shaped between-arrangement divergence profile
  consistent with long-suppressed recombination" ✅
- "candidates display an internal peak in differentiation consistent with
  locally-adapted alleles, although our hatchery design cannot test
  adaptation directly" ✅
- "this inversion is adaptive" ❌
- "this inversion is neutral" ❌
- "recombination is suppressed inside this inversion" — only after a
  separate per-karyotype recombination test

The atlas tooltip already shows "Class: shape verdict, not a causal claim" on
every pill. Don't strip that line.

## Roadmap

- **v0.2:** karyotype-purity gating in U01 (drop ambiguous samples below 0.80
  PC1 confidence). Currently relies on the upstream PCA being clean.
- **v0.3:** weight-AFD score by karyotype confidence rather than treating
  every assignment as equal weight.
- **v0.4:** integrate matched-background z-scores into the rule classifier
  (currently descriptive only).
- **v0.5:** TE-age (Kimura) co-layer cross-link in the candidate card —
  references the phase_8 module without duplicating its data.

## Cross-references

- breakpoint architecture, comparative-fragility, TE density →
  `phase_8_comparative_breakpoint_fragility/`
- karyotype calling and PCA-band assignment → upstream `pca_scrubber_v3.html`
  + `STEP_D17` boundary scan
- popstats-server endpoint contract → `Atlas/server_turn1/popstats_server.py`
- the planning bundle that motivated this module →
  `Atlas/handoffs_to_implement/turn129_bundle/`
