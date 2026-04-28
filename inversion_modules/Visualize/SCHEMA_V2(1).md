# SCHEMA_V2 — JSON layers contract for the local PCA scrubber

**Version:** 2.10 (2026-04-28)
**Status:** Specification — scrubber implements detection in v3.20, two-level
candidate promotion in v3.80, K=6 substructure in v3.81–v3.87. Marker /
PCR assay module (phase 13, §10) is specification-only — not yet
implemented. Band diagnostics module (§11) has scrubber-side compute and
UI in v3.92–v3.93; the three source layers (`theta_pi_panel`,
`roh_intervals`, `sample_froh`) are spec-only — emitted-from-R-side not
yet written. Het-shape sub-module (§11.1) shipped in v3.93 with three
in-scrubber placements (inline ridgeline column, candidate-page figure,
hover popover) using existing GHSL panel data. Live dosage heatmap
rendering contract added in v2.5; **v2.6 substantially extends it** with
the FIG_C08-style visual contract from `STEP28_candidate_marker_heatmap.R`
(5 left tracks + 1 top track), two new planned source layers
(`candidate_sample_coherence`, `candidate_marker_polarity`), and a
region-binding contract that supports both static (candidate-bound) and
live (cursor-bound, ±5 windows) views. **v2.7** adds the hover-tooltip
contract — overlay-only updates (no redraws), one-line format
"sample · marker (bp) · dosage=N · group · quality · flipped", with cell-
stable suppression and RAF coalescing. **v2.8** adds the boundary-zone
refinement contract (§12): an optional `boundary_evidence` layer,
per-candidate `boundary_left`/`boundary_right` records with weighted
score + support tracks + support_class, a candidate-level
`breakpoint_status` enum, and a strict terminology contract — "boundary
zone" is the default; "exact breakpoint" is reserved for junction-level
evidence only. **v2.9** flips three planned-status flags to shipped:
the v3.94 dosage heatmap has its R-side emit (STEP_M04 + STEP_M05 in
v3.95); v3.97's `boundary_evidence` layer has its R-side emit
(STEP_M06 in v3.99); the v3.98 NULL_CONCORD K=6 fix is documented in
the band-diagnostics module. No new contract surface in 2.9 — it's a
status-reconciliation bump. **v2.10** adds the within-regime
subclustering contract to §10: per-regime split → SVD-rotate `(u, v)`
→ standardize → DBSCAN/HDBSCAN — never raw `(u, v)` across all
samples. Four `cluster_mode` values defined: `none` (**default** —
subclustering not run, current scrubber behaviour preserved),
`rotated_within_regime` (opt-in production mode), and two named
diagnostic alternatives `raw_uv` and `one_dimensional_rotated_axis`.
Subcluster labels (`HOMO_1_sub2`, `HET_sub_noise`, etc.) ride on the
`candidate_sample_coherence` layer's `coarse_group_refined` column;
new optional `u_rot`/`v_rot`/`cluster_mode` columns. **Subclusters are
diagnostic, not paper-level inversion states** — manuscripts report
HOMO_1 / HET / HOMO_2 only.

This document defines the canonical JSON shape that flows between
`export_precomp_to_json_v3.R` (LANTA cluster), the scrubber (browser),
and the downstream cluster phases 6–12 of the inversion-popgen-toolkit.

---

## 1. Why a versioned schema

Until now, each scrubber feature was added by extending the JSON shape
without a version field. That worked while the scrubber was the only
consumer. Now that the scrubber emits a **candidate registry** that
cluster phases 6+ depend on, the shape becomes a contract — and
contracts need versions.

**Versioning rules:**

- `schema_version: 1` is the legacy shape (no `_layers_present` field).
  The scrubber **infers** which layers are present by checking which
  top-level keys exist.
- `schema_version: 2` is the new shape. It declares `_layers_present`
  explicitly so the scrubber doesn't have to infer.
- **Adding a new layer or column** → bump to 2.1, 2.2, etc. Backwards
  compatible; older readers ignore unknown fields.
- **Renaming or removing** a layer/column → bump to 3.x. Breaking change;
  cluster scripts must be updated in lockstep.

The scrubber accepts both v1 and v2 input. It writes v2 only.

---

## 1.5. File staging convention — the "chapter end" pattern

The scrubber doesn't talk to the cluster directly. It reads JSON files
that the user drag-drops into the browser. To make this manageable
across 28 chromosomes and many phases, all per-chromosome JSONs land in
**one staging area**, organized by chromosome.

### The pattern

```
inversion_modules/
└── scrubber/
    └── data/
        ├── README.md
        ├── LG01/
        │   ├── LG01_phase2_precomp.json        ← end of phase 2
        │   ├── LG01_phase3_proposals.json      ← end of phase 3
        │   ├── LG01_phase4a_ghsl.json          ← end of phase 4a
        │   ├── LG01_phase4b_theta.json         ← end of phase 4b
        │   ├── LG01_phase4c_dip.json           ← end of phase 4c
        │   ├── LG01_phase4d_concordance.json   ← end of phase 4d
        │   ├── LG01_phase4e_subcandidates.json ← end of phase 4e
        │   ├── LG01_phase5_consolidated.json   ← end of phase 5
        │   ├── LG01_phase9_classification.json ← end of phase 9
        │   └── LG01_phase12_cargo.json         ← end of phase 12
        ├── LG02/
        │   └── ... (same structure)
        ...
        └── LG28/
            └── ...
```

### Why one staging area

Reviewing 28 chromosomes is real work. The user shouldn't have to hunt
across `phase_2_discovery/2c_precomp/output/`, `phase_3_catalog_assembly/output/`,
`phase_4_resolution_and_classification/4d_dual_clustering/output/`, etc.
to find the JSONs for one chromosome. Everything for `LG14` lives in
`scrubber/data/LG14/`. Drag the folder's JSONs in, get to work.

### The chapter-end pattern

Each phase is a "chapter." Cluster scripts run, do their work, and
write their *own* outputs (RDS files, intermediate TSVs, audit logs)
to their *own* phase folders during the chapter. At the **end** of the
chapter — when the phase has finished — one designated emit script
distills what the scrubber needs and writes a single JSON to
`scrubber/data/<chrom>/`.

This means each phase has exactly one "emit-for-scrubber" step as its
last action. Not every script writes to `scrubber/data/`. Not even most
scripts. Just the chapter-closing emit script.

| Phase | Emit script (chapter-end)              | Output filename                              |
|-------|----------------------------------------|----------------------------------------------|
| 2     | `export_precomp_to_json_v3.R`          | `<chrom>_phase2_precomp.json`                |
| 3     | `STEP_X02_emit_proposals_for_scrubber.R` | `<chrom>_phase3_proposals.json`            |
| 4a    | `STEP_R03_emit_ghsl_for_scrubber.R`    | `<chrom>_phase4a_ghsl.json`                  |
| 4b    | `STEP_T06_emit_theta_for_scrubber.R`   | `<chrom>_phase4b_theta.json`                 |
| 4c    | `STEP_DT02_emit_dip_for_scrubber.R`    | `<chrom>_phase4c_dip.json`                   |
| 4d    | `STEP_DC05_diagnostic_panel.R`         | `<chrom>_phase4d_concordance.json`           |
| 4e    | `STEP_SC02_emit_subcandidates.R`       | `<chrom>_phase4e_subcandidates.json`         |
| 5     | `STEP_Q99_emit_consolidated.R`         | `<chrom>_phase5_consolidated.json`           |
| 9     | `STEP_CL99_emit_classification.R`      | `<chrom>_phase9_classification.json`         |
| 12    | `STEP_GC99_emit_cargo.R`               | `<chrom>_phase12_cargo.json`                 |
| 13    | `STEP_M99_emit_markers.R`              | `<chrom>_phase13_markers.json`               |

The naming convention is `<chrom>_<phase>_<role>.json`:

- `<chrom>` matches the JSON's `chrom` field exactly. Used for cross-chrom safety in the scrubber.
- `<phase>` is `phase2`, `phase3`, `phase4a`–`phase4e`, `phase5`, `phase9`, `phase12`.
- `<role>` is short and descriptive (`precomp`, `proposals`, `concordance`, `consolidated`).

The scrubber doesn't care what the file is named. It reads
`schema_version` and `_layers_present` from the JSON contents to know
what's inside. The naming convention is for human navigation, not for
the scrubber's parsing.

### Why "chapter end" not "every script"

Two architectural considerations:

**1. Audit trail vs. working area separation.** The phase folder is the
audit trail — it accumulates RDS files, intermediate TSVs, log files,
all the cluster's working outputs. The scrubber doesn't read any of
this. The `scrubber/data/<chrom>/` folder is the *working area* for
human review — only the distilled JSONs the scrubber needs.

You can `rm -rf scrubber/data/LG28/` and rebuild it by re-running each
phase's chapter-end script (no expensive recomputation needed if the
phase folder is preserved). You can't easily do the reverse — phase
folders are the source of truth.

**2. One emit step per phase = one place to change the schema.** When
the scrubber's expected JSON shape changes (a new layer, a renamed
field), only the chapter-end emit script needs updating, not every
script in the phase. This is what makes versioned schema viable.

### What if a phase wants to emit multiple files?

It can. Phase 4 has five sub-phases (4a–4e), each with its own emit
script and its own JSON. That's still "chapter end" — each sub-phase
is its own chapter. The user sees five files for phase 4 in
`scrubber/data/<chrom>/` and can drag them in independently as each
sub-phase finishes on the cluster.

What's not allowed: a single sub-phase emitting multiple JSONs
(e.g., `phase4d_concordance_part1.json`, `phase4d_concordance_part2.json`).
One sub-phase, one emit script, one JSON. If the data is too big for
one file, that's a sign the layer needs splitting into two
sub-phases or two distinct layer names.

---

## 2. The layer model

A v2 JSON file is a collection of optional layers, each owned by one
phase of the inversion-popgen-toolkit pipeline. Layers accumulate
monotonically: phase N's output preserves phase N-1's output and adds
its own keys.

**Top-level required fields (always present):**

```jsonc
{
  "schema_version": 2,
  "chrom": "C_gar_LG28",
  "n_windows": 4302,
  "n_samples": 226,
  "_layers_present": ["windows", "envelopes", "tracks", "sv_evidence",
                      "candidates_registry"],
  "_generated_at": "2026-04-27T14:30:00Z",
  "_generator": "export_precomp_to_json_v3.R"
}
```

`_layers_present` is the truth. If a layer name is in the array, the
corresponding top-level keys are guaranteed to exist and be valid.
If it is absent, the keys may be missing entirely.

**Layer ownership (corrected 2026-04-27 per parallel-chat architecture):**

The pipeline runs in this order. Each phase emits its own enrichment JSON
file per chromosome. Files don't merge on disk — the scrubber unions them
in-session as the user drag-drops them.

| Layer name              | Phase  | Owned by             | Written by                              | Status        |
|-------------------------|--------|----------------------|-----------------------------------------|---------------|
| `windows`               | 1+2    | precompute           | `export_precomp_to_json_v3.R`           | exists (v1)   |
| `envelopes`             | 1+2    | landscape            | `export_precomp_to_json_v3.R`           | exists (v1)   |
| `tracks`                | 1+2    | C++ pop-stats engine | `export_precomp_to_json_v3.R`           | exists (v1)   |
| `samples`               | 1+2    | sample registry      | `export_precomp_to_json_v3.R`           | exists (v1)   |
| `candidate_proposals`   | 3      | cluster GHSL+staging | `STEP_X01_emit_proposals.R` (planned)   | planned       |
| `cluster_labels_ghsl`   | 4a     | GHSL Part C          | `STEP_R01_ghsl_interval_kmeans.R`       | planned       |
| `cusum_ghsl`            | 4a     | GHSL Part D          | `STEP_R02_ghsl_cusum.R`                 | planned       |
| `ghsl_karyotype_runs`   | 4a     | GHSL Part B          | `export_ghsl_to_json_v2.R`              | exists (v1)   |
| `ghsl_heatmap`          | 4a     | GHSL panel (slim)    | `export_ghsl_to_json_v1.R` (deprecated) | superseded    |
| `ghsl_panel`            | 4a     | GHSL dense panel     | `export_ghsl_to_json_v2.R`              | exists        |
| `ghsl_kstripes`         | 4a     | per-K stripe assigns | `export_ghsl_to_json_v2.R`              | exists        |
| `cluster_labels_theta`  | 4b     | θπ Part C1+C2        | `STEP_T03+T04_theta_clustering.R`       | planned       |
| `cusum_theta`           | 4b     | θπ Part D            | `STEP_T05_theta_cusum.R`                | planned       |
| `dosage_dip`            | 4c     | Hartigan dip         | `STEP_DT01_dip_per_band.R`              | planned       |
| `concordance_tables`    | 4d     | dual-clustering      | `STEP_DC05_diagnostic_panel.R`          | planned       |
| `cusum_concordance`     | 4d     | DC06 cross-signal    | `STEP_DC06_cusum_concordance.R`         | planned       |
| `subcandidates_emitted` | 4e     | cell-mean cliff-walk | `STEP_SC01_emit_subcandidates.R`        | planned       |
| `candidates_registry`   | (user) | **scrubber**         | scrubber's "📋 export registry" button  | exists (v3.20)|
| `sv_evidence`           | 5      | Layer D OR test      | `STEP_D03_statistical_tests_*.py`       | exists, moves |
| `bnd_rescue`            | 5      | BND-paired junctions | `STEP_B06_bnd_rescue.py`                | exists, moves |
| `boundaries_refined`    | 5      | Q07d/Q07e cohort     | `STEP_Q07d_fst_boundary_refine.R`       | partial       |
| `qc_flags`              | 5      | Q01–Q11 battery      | `phase_5_qc_triage/`                    | partial       |
| `groups_validated`      | 5      | tier upgrade         | (rolled-up phase 5 output)              | planned       |
| `classification`        | 9      | terminal class       | `assign_structural_class_v7.py`         | exists        |
| `gene_cargo`            | 12     | gene content         | (planned)                               | planned       |
| `marker_catalogue`      | 13     | regime-diagnostic markers | `STEP_M01_emit_marker_catalogue.R` | planned       |
| `marker_primers`        | 13     | primer designs       | `STEP_M02_emit_marker_primers.R`        | planned       |
| `marker_panel_summary`  | 13     | panel rollup         | `STEP_M03_emit_marker_panel_summary.R`  | planned       |
| `theta_pi_panel`        | 4b     | dense θπ matrix      | `STEP_T06_emit_theta_pi_panel.R`        | planned (v3.92) |
| `roh_intervals`         | 3 (MODULE_3 ROH) | per-sample ROH spans | `STEP_R12_emit_roh_intervals.R`         | planned (v3.92) |
| `sample_froh`           | 3 (MODULE_3 ROH) | genome-wide FROH per sample | `STEP_R13_emit_sample_froh.R`     | planned (v3.92) |
| `candidate_sample_coherence` | 13 follow-up | per-(candidate × sample) coherence + `stripe_quality` | `STEP_M05_emit_step29_layers.R` (wraps `STEP29_candidate_coherence_and_polarity.R` TSVs) | shipped v3.95 |
| `candidate_marker_polarity`  | 13 follow-up | per-(candidate × marker) polarity (L1, L2, final flip, support class, block id) | `STEP_M05_emit_step29_layers.R` (wraps `STEP29_candidate_coherence_and_polarity.R` TSVs) | shipped v3.95 |
| `dosage_chunks`         | 13 follow-up | per-region dosage matrix tiles (LRU loaded by scrubber on demand) | `STEP_M04_emit_dosage_chunks.R` | shipped v3.95 |
| `boundary_evidence`     | 12 (boundary refinement) | per-candidate FST + per-regime θπ + discordant_pair_pileup + sv_anchors at 5 kb scan-window resolution | `STEP_M06_emit_boundary_evidence.R` | shipped post-v3.98 |

**Architecture notes:**

- **Layer D (SV/breakpoint) lives in phase 5, not phase 3.** This was clarified in the parallel-chat decision: phase 3 is pure catalogue assembly (no validation), and Layer D needs intervals as input so it violates the phase-3 contract. SV evidence becomes a phase-5 QC validation step.
- **Two distinct candidate layers:**
  - `candidate_proposals` — cluster-authored, automated dedup of phase-2 streams. Read-only from the scrubber's perspective; loaded on phase-3 enrichment import to populate `state.candidateList` as a starting point.
  - `candidates_registry` — scrubber-authored, the user's curated list (proposals + edits + accept/reject). Emitted by the "📋 export registry" button. This is what cluster phases 4+ read as their candidate input.
- **The dual-clustering page (page 6 concordance)** lights up when `concordance_tables` is in `_layers_present` — i.e., after phase 4d output is loaded.
- **The validation banner on page 6** lights up when `sv_evidence` is in `_layers_present` — i.e., after phase 5 output is loaded.
- **CUSUM layers (`cusum_ghsl`, `cusum_theta`, `cusum_concordance`)** carry per-sample changepoint tables, per-candidate boundary distributions (5'/3' median + IQR + spread_class tight/ragged), and sub-system clusters. These are manuscript-grade observations.
- **`ghsl_panel` is the dense back-end for L2-as-triangle aggregation.** It carries the full (sample × window × scale) divergence matrix plus per-sample rank and rank-band at each scale. Browser-side helpers (`ghslDivAt`, `ghslAggregateRange`, `ghslAggregateFocal`, `ghslAggregateInterval`) compute per-L2 aggregates on click without a cluster round-trip. JSON size 50–300 MB per chromosome; loads fine locally.
- **`ghsl_kstripes` carries per-K stripe assignments** (K=2..6) at the primary rolling scale. Stripe 1 is always the lowest-rank stripe; per-stripe mean and median ranks are pre-computed for instant catalogue summaries.
- **Focal-with-flanks aggregation:** for wide L2s, mean GHSL across the full L2 range can dilute the signal (gene-conversion erosion, sub-systems). The focal±N viewer (`ghslAggregateFocal`, `ghslAggregateInterval`) takes a focal window plus a half-flank radius (default 5 windows = 11 windows total) and aggregates only that. Used by tracked-samples view and per-L2 inspection. Range-aggregation (`ghslAggregateRange`) stays available for catalogue stats where the full-L2 view is the right question.

The scrubber renders UI conditional on `_layers_present`. Layers it does
not recognize are silently ignored (so future expansions don't break it).

---

## 2.5. Interval Collector (planned, session 3)

The Interval Collector is a **passive accumulator** of candidate intervals
derived from the interaction between tracked samples and the GHSL panel.
It runs automatically when tracked samples change; the user does not have
to open a page for collection to happen. A page exists for inspection but
collection is not gated on it.

### Data flow

```
tracked_samples (user-driven, changes as user clicks ⭐ on samples)
   +
ghsl_panel (loaded from phase 4a JSON)
   +
ghsl_kstripes (loaded from phase 4a JSON)
   ↓
collectorEngine() — runs in scrubber, no I/O
   ↓
state.collectedIntervals = [
  { interval_id, start_bp, end_bp, start_w, end_w,
    n_supporting_pairs, mean_strength,
    captured_at, source: "tracked_pair_disagreement",
    derived_from_samples: ["CGA042", "CGA107"],
    K_used: 3 }
  ...
]
```

### Algorithm sketch (session 3 implementation)

For each pair of currently-tracked samples (s_i, s_j):

1. Look up `ghsl_kstripes.by_k["3"].stripe_per_sample[i]` and `[j]`.
2. If samples are in different stripes globally, scan windows: at each
   window w, compute `rank_band[s_i][w]` and `rank_band[s_j][w]` (LOW/MID/HIGH).
3. Find consecutive runs where the pair agrees (same band) and disagrees
   (different bands). Boundaries are at the agree↔disagree transitions.
4. Each transition window contributes one boundary signal for this pair.

Then aggregate across all tracked pairs:

5. For each window w, count how many pairs have a boundary signal at w
   (or within ±2 windows for smoothing).
6. Find local maxima in this count vector. Each maximum above a threshold
   becomes a captured boundary.
7. Consecutive captured boundaries define intervals — each interval gets
   an entry in `state.collectedIntervals`.

### Storage and persistence

- `state.collectedIntervals` is per-chromosome (keyed by chrom).
- Persisted to `localStorage` like the candidate list.
- Cleared when chromosome changes (separate per chrom; can be browsed back).
- Exported with the candidate registry as a separate column or as
  derived-from rows.

### What goes into the registry export

Two options for review at session 3 time:

- **Option A.** Collector intervals are auto-promoted into the candidate
  registry as `discovery_method = "scrubber_collector"` candidates, alongside
  user-promoted L2 candidates. The cluster phase 4 picks them up as proposals.
- **Option B.** Collector intervals stay in a separate "candidates_collected"
  layer. The user reviews them on a page and can promote individual ones
  into the registry. More conservative.

Default leans toward (B) — collector signals can be noisy without further
QC and we don't want auto-promotion to pollute the registry. User decides
per interval.

### Why "no need to open the page"

The collector runs whenever:
- A sample is added/removed from tracked_samples
- The user changes K (k-stripes choice)
- A new ghsl_panel/ghsl_kstripes is loaded

The output (`state.collectedIntervals`) is queryable from any page —
it's just state. The collector page exists to visualize the captured
intervals (positions, supporting pairs, strength) but UI rendering is
not on the critical path. Other pages can show "N intervals collected"
in the header.

---

## 3. v1 → v2 inference

When the scrubber loads a v1 file (no `schema_version` field, or
`schema_version: 1`), it **infers** the layers from top-level keys:

```js
function inferLayersFromV1(json) {
  const layers = [];
  if (json.windows && json.n_windows > 0)        layers.push('windows');
  if (json.l1_envelopes || json.l2_envelopes)    layers.push('envelopes');
  if (json.tracks && Object.keys(json.tracks).length > 0) layers.push('tracks');
  if (json.samples && json.samples.length > 0)   layers.push('samples');
  // sv_evidence, breakpoints_refined, etc. were not produced by v1 exporters
  return layers;
}
```

This means a v1 JSON behaves identically to a v2 JSON with
`_layers_present: ['windows', 'envelopes', 'tracks', 'samples']` —
which is exactly the case today.

---

## 4. The candidates_registry layer (the scrubber's contract output)

This is the **canonical candidate list** that cluster phases 6+ consume.
It is written by the scrubber when the user clicks
*"📋 export candidate registry"* on page 4.

Two parallel files are produced:

- `<chrom>_candidate_registry.tsv` — canonical, one row per candidate.
  Easy to read in R / Python / Excel. **This is what cluster phases read.**
- `<chrom>_candidate_registry.json` — same data as JSON. For tools that
  prefer JSON (e.g. embedded in a v2 layer JSON).

### TSV columns (canonical, fixed order)

| Column                  | Type    | Source           | Description                                                                                       |
|-------------------------|---------|------------------|---------------------------------------------------------------------------------------------------|
| `candidate_id`          | string  | scrubber         | Unique id (e.g. `cand_lg28_142_a`); stable across reruns of the scrubber                          |
| `chrom`                 | string  | scrubber         | Chromosome (matches the JSON's `chrom`)                                                           |
| `start_bp`              | integer | scrubber         | Candidate start, 1-based inclusive                                                                |
| `end_bp`                | integer | scrubber         | Candidate end, 1-based inclusive                                                                  |
| `start_w`               | integer | scrubber         | First window index (0-based)                                                                      |
| `end_w`                 | integer | scrubber         | Last window index (0-based, inclusive)                                                            |
| `n_windows`             | integer | scrubber         | end_w - start_w + 1                                                                               |
| `span_mb`               | float   | scrubber         | (end_bp - start_bp) / 1e6, 3 decimals                                                             |
| `K`                     | integer | scrubber         | Number of bands (karyotype groups)                                                                |
| `source`                | string  | scrubber         | `l2_single` / `l2_merge` / `lock_promote` (how the candidate was created)                         |
| `l2_indices`            | string  | scrubber         | Comma-separated L2 indices (e.g. `"3,4,5"`)                                                       |
| `ref_l2`                | integer | scrubber         | The L2 used as the band reference                                                                 |
| `ref_window`            | integer | scrubber         | Window index where bands were anchored                                                            |
| `sigma_verdict`         | string  | scrubber         | `TWO_INVERSIONS` / `CROSSOVER_ARTIFACTS` / `NOISY_REGION` / `NA` / `UNDETERMINED`                  |
| `sigma_q50`             | float   | scrubber         | 50th percentile of per-sample σ across the candidate span                                         |
| `n_drifters`            | integer | scrubber         | Samples with σ > 2·q50                                                                            |
| `top_family_largest_band` | string | scrubber         | Family ID (e.g. `F4`) of the dominant family in the largest band                                  |
| `top_family_pct`        | float   | scrubber         | Fraction of the largest band held by `top_family_largest_band`                                    |
| `n_supporting_intervals` | integer | scrubber (v3.80) | Number of L2 intervals that support this candidate (= length of `l2_indices`)                    |
| `aggregate_concordance` | float   | scrubber (v3.80) | Mean Hungarian-aligned K=3 concord between each supporting L2 and the anchor (`ref_l2`)          |
| `band_continuity_pct`   | float   | scrubber (v3.80) | `same_rate + adj_rate` from `bandContinuity` — fraction of (sample × adjacent-L2-pair) where the K=3 PC1-rank label stays the same or moves by one |
| `band_continuity_verdict` | string | scrubber (v3.80) | `NESTED` / `INTERMEDIATE` / `UNSTABLE` / `NO_DATA` (thresholds: `cont ≥ 0.85 + far ≤ 0.05` → NESTED; `cont ≤ 0.50 OR far ≥ 0.20` → UNSTABLE) |
| `regime_counts`         | string  | scrubber (v3.80) | Per-K=3-regime fish counts: `g0:N0;g1:N1;g2:N2`. Excludes ambiguous fish.                        |
| `n_ambiguous_fish`      | integer | scrubber (v3.80) | Number of fish with `confidence < ambiguous_threshold` (default 0.5)                              |
| `subband_nesting_status` | string | scrubber (v3.81) | K=6 → K=3 nesting verdict: `NESTED` (all K=6 groups pure ≥ 0.80) / `MIXED` (some pure, some not) / `CROSS_CUTTING` (none pure) / `NA` (no K=6 pass) |
| `n_k6_groups`           | integer | scrubber (v3.81) | Count of K=6 groups with at least one supporting observation                                      |
| `n_k6_pure`             | integer | scrubber (v3.81) | Count of K=6 groups with purity ≥ `purity_threshold` (default 0.80)                              |
| `qc_status`             | string  | scrubber         | `accepted` / `rejected` / `review` (set by user; default `accepted` for saved candidates)         |
| `notes`                 | string  | scrubber         | Free-text notes (newlines and tabs replaced with spaces)                                          |
| `created_at`            | string  | scrubber         | ISO 8601 timestamp                                                                                |
| `scrubber_version`      | string  | scrubber         | The scrubber that wrote this row (e.g. `v3.20`)                                                   |

**Stability promise:** column names will not change without a schema bump.
New columns may be appended (with a 2.x bump). Existing columns will not
be reordered.

### JSON mirror

```jsonc
{
  "schema_version": 2,
  "chrom": "C_gar_LG28",
  "scrubber_version": "v3.20",
  "exported_at": "2026-04-27T14:30:00Z",
  "n_candidates": 12,
  "candidates": [
    {
      "candidate_id": "cand_lg28_142_a",
      "chrom": "C_gar_LG28",
      "start_bp": 13240000,
      "end_bp": 17110000,
      // ... all TSV columns as object keys
      "locked_labels": [0, 0, 1, 2, 0, ...]   // also retained for round-trip into scrubber
    }
  ]
}
```

The `locked_labels` field is JSON-only (not in TSV). It lets the user
reload the registry into the scrubber and have the bands restored
exactly as decided.

---

## 5. The downstream contract for cluster phases 6+

A cluster script that wants to enrich the candidate registry must:

1. **Read** `<chrom>_candidate_registry.tsv` (or `.json`).
2. **Compute** its layer's payload, keyed by `candidate_id`.
3. **Emit** an enrichment JSON with at minimum:

```jsonc
{
  "schema_version": 2,
  "chrom": "C_gar_LG28",
  "_layers_present": ["breakpoints_refined"],          // just the new layers added
  "_generator": "STEP_PHASE6_breakpoint_refine.R",
  "_generated_at": "2026-04-27T17:30:00Z",
  "breakpoints_refined": {
    "cand_lg28_142_a": {
      "start_bp_refined": 13241827,
      "end_bp_refined":   17109934,
      "start_ci_low": 13241500, "start_ci_high": 13242100,
      "end_ci_low":   17109800, "end_ci_high":   17110300,
      "method": "DELLY_PR_split_consensus",
      "n_supporting_reads": 47
    },
    // ...
  }
}
```

The scrubber then loads this enrichment file (drag-and-drop, alongside
the original precomp JSON) and detects `breakpoints_refined` in
`_layers_present`. The corresponding UI (column on page 3 / panel on
page 2 / its own page) becomes available.

**Layer files do not need to be merged on disk.** The scrubber is happy
to load multiple JSONs and union their `_layers_present` arrays. This
keeps each cluster phase's output as its own audit-able artifact.

---

## 6. Versioning checklist for the maintainer (you, future-you, me)

When making a change:

- [ ] **Adding a column to candidates_registry** → bump to 2.1, append column at end of TSV.
- [ ] **Renaming or removing a column** → bump to 3.0, document which cluster scripts need updating.
- [ ] **Adding a brand-new layer** (e.g. `gene_cargo`) → bump to 2.1, document in §2 table.
- [ ] **Changing the meaning of an existing field** (e.g. start_bp from 1-based to 0-based) → bump to 3.0, document.
- [ ] **Internal scrubber refactor that doesn't change JSON shape** → no bump needed.

Keep this file updated. The scrubber's `console.log` on JSON load
should also report the schema version it detected.

---

## 7. Two-level candidate promotion (v3.80, schema 2.1)

The candidate registry was originally a flat per-candidate row. v3.80
restructured the export around a **two-level promotion model**:

- **Level 1 (interval → candidate)** — supporting L2 intervals are merged
  into one `candidate_registry` row with aggregate evidence (concordance,
  band continuity, regime counts).
- **Level 2 (fish → regime within candidate)** — each fish gets one row
  per candidate it appears in, with a consensus K=3 regime call,
  confidence, and (when K=6 substructure is computed) a sub-band path.

**Why two levels.** The fish does not "belong to interval L2_03." The fish
belongs to a *regime* (g0/g1/g2) inside a *candidate*. Conflating those two
operations was the original sin of the v3.20 schema — `top_family_pct` and
`locked_labels` were doing duty as both per-candidate aggregates and
per-fish summaries. v3.80 separates them:

```
Chromosome
└── Inversion candidate                              (Level 1 row in candidate_registry.tsv)
    ├── supporting L2/L3 intervals                   (rows in interval_support.tsv)
    ├── evidence summary                             (columns on candidate_registry.tsv)
    └── fish regime calls                            (rows in fish_regime_calls.tsv)
        ├── g0 / g1 / g2     ← primary K=3 regime
        └── g0a/g0b/g1a/...  ← K=6 sub-band annotation (v3.81)
```

### Three-table export

When the user clicks *"📋 export candidate registry"*, the scrubber writes
**four files** (was two in v3.20):

| File                                       | Rows                                  | Level | Added in |
|--------------------------------------------|---------------------------------------|-------|----------|
| `<chrom>_candidate_registry.tsv`           | one per candidate                     | 1     | v3.20    |
| `<chrom>_candidate_registry.json`          | round-trippable mirror + `locked_labels` | 1  | v3.20    |
| `<chrom>_fish_regime_calls.tsv`            | one per (fish × candidate)            | 2     | v3.80    |
| `<chrom>_interval_support.tsv`             | one per (interval × candidate)        | 1     | v3.80    |

### `<chrom>_fish_regime_calls.tsv` (v3.80, extended in v3.81)

One row per `(sample × candidate)` pair. Records the consensus regime
call — what cluster phases will read as the fish's karyotype call for
that candidate.

| Column              | Type    | Description                                                                                       |
|---------------------|---------|---------------------------------------------------------------------------------------------------|
| `sample_id`         | string  | Sample identifier (CGA / individual)                                                              |
| `candidate_id`      | string  | Candidate this row belongs to                                                                      |
| `regime`            | string  | `g0` / `g1` / `g(K-1)` for the consensus regime, or `ambiguous` when `confidence < ambiguous_threshold` |
| `major_regime`      | string  | (v3.81) Alias of `regime` for K=3 majority. Empty string for ambiguous. Made explicit so downstream consumers don't need to remember which is the major. |
| `confidence`        | float   | Fraction `n_supporting / n_intervals` for the modal K=3 regime                                    |
| `n_intervals`       | integer | Number of supporting L2 intervals where the fish had a valid call                                 |
| `n_supporting`      | integer | Count of intervals voting for the modal regime                                                    |
| `ambiguous`         | string  | `TRUE` / `FALSE`                                                                                  |
| `subband_path`      | string  | (v3.81) Pipe-separated K=6 sub-bands across supporting intervals: `g0a|g0b|g0a`. Empty when K=6 unavailable. |
| `subband_stability` | float   | (v3.81) Fraction of intervals where the fish's K=6 sub-band's parent matches the K=3 regime. `1.0` = fish stays inside one major regime; `<1` = K=6 path crosses parents. `NA` when no K=6 pass. |
| `votes`             | string  | Pipe-separated per-interval aligned K=3 labels: `g0|g0|g1` for audit                              |

**Consensus rule.** For each fish, collect its Hungarian-aligned K=3
label at every supporting L2. The consensus regime is the modal label;
confidence is `count(mode) / n_valid`. If the max-fraction drops below
`ambiguous_threshold` (default 0.5), the fish is flagged `ambiguous` and
its regime is set to the empty string in `major_regime` and the literal
`"ambiguous"` in `regime`.

**PC1-rank ordering.** Before Hungarian alignment, every L2's labels are
ranked by PC1 centroid so `g0` always means "lowest-PC1 cluster" and
`g(K-1)` always means "highest-PC1 cluster". This makes regime labels
biologically stable across L2s and across candidates.

### `<chrom>_interval_support.tsv` (v3.80)

One row per `(interval × candidate)` pair. Records which L2 intervals
back each candidate and what role each plays.

| Column            | Type    | Description                                                                          |
|-------------------|---------|--------------------------------------------------------------------------------------|
| `candidate_id`    | string  | Candidate this row belongs to                                                         |
| `l2_idx`          | integer | L2 envelope index in the source `precomp.json`                                       |
| `l2_label`        | string  | L2 envelope's `candidate_id` from precomp (e.g. `L2_0001/03`)                        |
| `role`            | string  | `core` (the anchor L2, used as Hungarian-alignment basis) / `support` (others)        |
| `start_bp`        | integer | L2 envelope `start_bp`                                                               |
| `end_bp`          | integer | L2 envelope `end_bp`                                                                  |
| `n_windows`       | integer | `_e0 - _s0 + 1`                                                                      |
| `concord_to_ref`  | float   | Hungarian-aligned K=3 concord vs the anchor `ref_l2`. `1.0` for the core L2 itself.  |

The `concord_to_ref` column is a per-supporting-L2 read of the same
quantity that's *averaged* into `aggregate_concordance` on the candidate
registry — useful for deciding which L2s are tightly aligned vs loosely
when reviewing why a candidate's aggregate is below 1.

---

## 8. K=6 nested substructure semantics (v3.81, schema 2.1)

K=3 is the **primary regime** call (the karyotype). K=6 is **substructure**
inside K=3. K=6 only gets promoted to a separate inversion call when it
fails to nest cleanly inside K=3.

### The nesting decision rule

For each L2 interval, cluster at K=6 in addition to K=3, then map K=6
groups to their K=3 parents by max-overlap:

```
K=6 group | parent K=3 regime
k0        | g0  (most samples in k0 are g0 in K=3)
k1        | g0
k2        | g1
k3        | g1
k4        | g2
k5        | g2
```

Each K=6 group's **purity** is `max_overlap / total_in_group`. A K=6
group with purity ≥ `purity_threshold` (default 0.80) is *cleanly nested*
in its K=3 parent. The candidate's `subband_nesting_status` is then:

| Verdict          | Definition                                          | Interpretation                              |
|------------------|-----------------------------------------------------|---------------------------------------------|
| `NESTED`         | All K=6 groups with data have purity ≥ threshold    | One inversion with substructure. Keep K=6 as `g0a/g0b/...` annotation. Do NOT promote K=6 groups to separate inversions. |
| `MIXED`          | Some K=6 groups pure, others not                    | Partial nesting. Review per-K6-group purity before deciding.                                |
| `CROSS_CUTTING`  | No K=6 groups pure                                  | Possible second system / independent inversion overlapping this candidate.                  |
| `NA` / `NO_DATA` | K=6 cluster pass unavailable                        | Older candidate written before v3.81 OR cluster failed at K=6                              |

### Sub-band labels: `g{parent}{letter}`

Within each K=3 parent, K=6 groups are ranked by PC1 centroid and
assigned ordinal letters (`a`, `b`, ...). So a fish's per-L2 sub-band
path looks like `g0a → g0b → g0a` (stayed within g0) or
`g0a → g1a → g0b` (crossed parents — `subband_stability < 1`).

### Why K=3 is the primary call

The biology is the K=3 regime: a fish is a homozygous or heterozygous
karyotype call. K=6 captures finer haplotype variation within those
regimes — useful for diagnostic atlas work and the manuscript's "nested
arrangements" framing, but not the karyotype call itself.

The published wording for this:

> *Major K=3 regimes were used as candidate-level arrangement classes.
> Higher-K clusters were retained as within-regime sub-bands and were not
> interpreted as independent inversions unless they cut across the K=3
> structure or showed independent PC-axis support.*

That sentence's preconditions (the cut-across detection) are exactly what
`subband_nesting_status` exposes.

### The `k6_substructure` field on the JSON-mirror

The round-trippable `<chrom>_candidate_registry.json` retains a per-
candidate `k6_substructure` object (in addition to the TSV's flat
columns) so the scrubber can rehydrate a saved candidate's substructure
without re-running clustering:

```jsonc
"k6_substructure": {
  "K6": 6,
  "parent_of_k6": [0, 0, 1, 1, 2, 2],   // K=3 parent for each K=6 group
  "purity": [1.0, 0.95, 1.0, 0.88, 1.0, 0.97],
  "n_pure": 6, "n_with_data": 6,
  "verdict": "NESTED",
  "purity_threshold": 0.80,
  "subband_letter": ["a", "b", "a", "b", "a", "b"]   // ordinal letter per K=6 group
}
```

### Computation cost

Computed once at candidate-commit time (`commitL3Draft` in scrubber).
Two passes per candidate (K=3 + K=6) over the Hungarian-aligned per-L2
labels. Total cost is O(n_intervals × n_samples × K!) — small for the
typical few-L2 candidate. Subsequent loads of saved candidates rehydrate
the stored `k6_substructure` without recomputing.

### When to flag CROSS_CUTTING for biological review

The verdict alone isn't biology. Use `subband_nesting_status =
CROSS_CUTTING` as a **filter** for which candidates need extra scrutiny,
not as evidence by itself. A cross-cutting K=6 in a region where the
breakpoint geometry, σ profile, and concordance tables all support a
single inversion may simply mean K=6 over-clusters at that K — review
the K-mode 'both' L3 panel and the per-fish `subband_path` patterns
before splitting into two candidates.

---

## 9. Versioning history

| Schema | Date       | Scrubber | Changes                                                        |
|--------|------------|----------|----------------------------------------------------------------|
| 1      | (legacy)   | <v3.20   | No version field; layers inferred from key presence             |
| 2      | 2026-04-27 | v3.20    | `_layers_present` declared; candidate registry contract         |
| 2.1    | 2026-04-28 | v3.80–v3.87 | Two-level promotion: candidate (L1) + fish regime (L2). Added 9 columns to candidate_registry; added 2 sibling TSV files (fish_regime_calls.tsv + interval_support.tsv). Added K=6 nested substructure (v3.81) with NESTED/MIXED/CROSS_CUTTING verdict and `g0a/g0b/...` sub-band annotation. |
| 2.2    | 2026-04-28 | (planned)   | Marker / PCR assay module (phase 13). Added three layers (`marker_catalogue`, `marker_primers`, `marker_panel_summary`) and three sibling TSVs for hatchery-deployable diagnostic panels per inversion candidate. |
| 2.3    | 2026-04-28 | v3.92       | Band diagnostics module (§11). Per-band confounder/support metrics computed on focal L2 and stored on candidate at commit time. Added three planned source layers (`theta_pi_panel`, `roh_intervals`, `sample_froh`) and a sibling TSV (`<chrom>_band_diagnostics.tsv`). Scrubber-side compute + render + export are live; R emit scripts for the three new layers are planned. Existing GHSL panel and dosage het reads are wired today. |
| 2.4    | 2026-04-28 | v3.93       | Het-shape sub-module (§11.1). `computeBandDiagnostics` returns a new `het_shape` field with 16-bin histograms per band + support-ratio (K=3 → middle / mean(flanks); else max / min spread). Three in-scrubber placements share one renderer: inline SVG ridgeline column in the L3 contingency tile diagnostics table, FIG_C07-style stacked ridgeline figure on the candidate page, click-popover from the het mini-chip pill. TSV got 4 new columns (`het_support_ratio`, `het_support_kind`, `het_support_passes`, `het_hist_bins`). No new source-layer requirement — uses existing GHSL panel. |
| 2.5    | 2026-04-28 | (planned)   | Live dosage heatmap rendering contract added to §10 (Marker module). Forward spec — no scrubber implementation yet. Defines: visible-region = current 5–10 local PCA windows; default 200 markers, hard cap 500; selection cascade (missingness filter → diagnostic_score → dosage variance fallback); per-region dosage chunk JSON shape (`<chrom>_dosage_chunk_<start>_<end>.json`); LRU chunk cache; explicit redraw triggers (K, row order, visible window, cap, chunk-load); 150 ms debounce on scrub. No implementation, no test impact. |
| 2.6    | 2026-04-28 | (planned)   | Substantial expansion of the dosage heatmap contract (still planned) using `STEP28_candidate_marker_heatmap.R` as the visual reference. Five left-side annotation tracks (K=3 PCA group, K=6 sub-band, Quality, GHSL, het) and one top-side track (Polarity). Sample order = `coarse_group_refined, u`. Dosage ramp `#2166AC / #F7F7F7 / #B2182B`. Two view modes share one renderer: candidate-bound (static) and cursor-bound (live, ±5 windows). Marker cap independent of region source. Stripe-quality button when coherence layer absent (faithful port of STEP29 algorithm). Two new planned source layers: `candidate_sample_coherence`, `candidate_marker_polarity` (both have R-side reference in `STEP29_candidate_coherence_and_polarity.R`). Placement (popup vs docked vs floating) deliberately undecided — to be chosen during v3.94 implementation. |
| 2.7    | 2026-04-28 | v3.96       | Hover-tooltip contract added to §10 ("Hover affordances"). One-line format `sample · marker (bp) · dosage=N · group · quality · flipped` with absent-field omission (no `null` placeholders). NA dosage shown as `dosage=NA`; the `-1` sentinel is never user-facing. Schema-mandated as overlay updates only — never a redraw. Cell-stable suppression and RAF coalescing for performance. Two hover contexts (`hover_static`, `hover_cursor`) on `state.__dosageHm`, both reset on every redraw. Annotation tracks (left + top) do not currently tooltip — only the dosage matrix. R-side unchanged. |
| 2.8    | 2026-04-28 | (planned)   | Boundary zone refinement module (§12). New optional `boundary_evidence` layer (cluster-side R emit, planned, with FST/theta-pi/discordant-pair/SV-anchor tracks). Per-candidate `boundary_left`/`boundary_right` records carrying zone bp range, weighted score 0..1, contributing tracks list, support_class (strong/moderate/weak/ambiguous), source (auto/manual/auto_plus_manual), and SV anchors in zone. Candidate-level `breakpoint_status` enum: boundary_zone_only / SV_supported / junction_supported. Auto-propose algorithm: per-track step magnitude in MAD units → weighted combined score → strongest edges in left/right halves → ±5 window zones. Manual override via E/F keys at cursor; B saves; R resets. New "boundaries" page (page 11 internally, displayed as tab "3 boundaries") wedged between candidate-focus and catalogue. 14 new TSV columns (`left_boundary_zone_start` ... `boundary_notes`); registry width 31 → 45. **Terminology contract enforced**: "boundary zone" is default; "exact breakpoint" reserved for junction-level evidence only. |
| 2.9    | 2026-04-28 | v3.98 (+ R-side STEP_M06 post-v3.98) | Status-reconciliation bump — no new contract surface. Three flags flipped from "planned" to "shipped": (1) v3.94 dosage heatmap layer (`dosage_chunks`, `candidate_sample_coherence`, `candidate_marker_polarity`) has its R-side emit (`STEP_M04_emit_dosage_chunks.R` + `STEP_M05_emit_step29_layers.R`, shipped v3.95); (2) v3.97 boundary refinement layer (`boundary_evidence`) has its R-side emit (`STEP_M06_emit_boundary_evidence.R`, shipped post-v3.98) with FST + per-regime θπ (homo1/het/homo2) + discordant_pair_pileup + sv_anchors; track-by-track availability — each track independently optional based on inputs supplied; (3) v3.98 NULL_CONCORD K=6 → 0.24 (was incorrectly 0.30 fallback) documented in band-diagnostics module. v3.80 UI text rewrite (intervals → promoted → candidates; fish → assigned → regimes) shipped in v3.98 with new help-page Vocabulary section. **Still planned in v3.x**: `STEP_T06_emit_theta_pi_panel.R`, `STEP_R12_emit_roh_intervals.R`, `STEP_R13_emit_sample_froh.R` (band-diagnostics §11 source layers); marker module phase-13 R pipeline (§10). |
| 2.10   | 2026-04-28 | (planned)   | Within-regime subclustering contract added to §10 (dosage heatmap module). DBSCAN/HDBSCAN must NOT run on raw `(u, v)` across all samples — the major-regime axis dominates and within-regime structure gets washed out. Per-regime contract: split by `coarse_group` (HOMO_1 / HET / HOMO_2) → SVD-rotate within-regime `(u, v)` to `(u_rot, v_rot)` → z-score standardize → DBSCAN/HDBSCAN. Subcluster labels are within-regime, dense, 1-based: `HOMO_1_sub1`, `HOMO_1_sub_noise`, `HET_sub2`, etc. **Subclusters are diagnostic within-regime haplotype backgrounds, NOT paper-level inversion states** — manuscript reports `HOMO_1 / HET / HOMO_2` only. New `cluster_mode` parameter with four values: `none` (**default** — subclustering not run, preserves current scrubber behaviour), `rotated_within_regime` (opt-in production mode — the 5-step contract), `raw_uv` (diagnostic alternative — DBSCAN on raw `(u, v)`; almost always recovers the major-regime split), `one_dimensional_rotated_axis` (diagnostic alternative — per-regime split + 1D `u_rot` only, fallback for low-coverage where `v_rot` is noise-dominated). Default may flip to `rotated_within_regime` in a future schema bump after calibration on the 226-sample LG28 cohort. Labels ride on `candidate_sample_coherence` layer's existing `coarse_group_refined` column; new optional `u_rot`/`v_rot`/`cluster_mode` columns. Heatmap K=3 PCA group track gains an inner sub-stripe under the two rotated modes. Recompute trigger: candidate commit time; changing `cluster_mode` invalidates stored labels. |

---

## 10. Marker / PCR assay module (planned, schema 2.2)

The marker module closes the loop from inversion discovery to farm-usable
diagnostic assays. For each high-confidence inversion candidate, it
designs a small panel of PCR markers that distinguish K=3 regimes
(g0/g1/g2), letting hatchery staff assign arrangement classes from
low-cost genotyping rather than from whole-genome sequencing.

### Why "candidate-regime diagnostic markers", not "breakpoint PCR"

Breakpoint PCR requires exact junction sequences on both inversion
arrangements. We do not have those — short-read assemblies in the
226-sample cohort were not built for breakpoint-resolution work, and the
F1 hybrid reference's breakpoint coordinates carry assembly-level
uncertainty. What we **do** have, robustly, is the regime call per fish
per candidate (`fish_regime_calls.tsv`, schema §7) and dense per-sample
allele information from the resequencing.

So the marker logic is variant-frequency-driven, not junction-driven:

```
candidate inversion regime calls
       ↓
find SNPs/indels enriched in each regime (g0 / g1 / g2)
       ↓
select robust markers (specificity, MAF, family-spread, locus quality)
       ↓
design primer pairs around each marker
       ↓
rank and assemble panels (3–10 markers per candidate)
       ↓
export hatchery genotyping sheet
```

The output framing for the manuscript:

> *For each high-confidence candidate structural haplotype regime, we
> designed candidate diagnostic PCR marker panels to support low-cost
> screening of arrangement classes in hatchery populations.*

This phrasing is deliberate. We are not claiming breakpoint detection or
direct structural assays — the panels are *diagnostic of the regime call*,
which is itself a population-genomics product. The panels are validated
by whether they reproduce the K=3 regime call from the genome data, not
by whether they detect a junction.

### Four marker classes

For each candidate the module attempts four marker classes:

| Class                         | What it tags                                  | Required prerequisite                          |
|-------------------------------|-----------------------------------------------|-------------------------------------------------|
| **Regime-specific**           | Variant common in one regime, rare in others  | `freq_target ≥ 0.80` AND `freq_others ≤ 0.10`   |
| **Heterozygote-pattern**      | Variant whose g1 frequency is intermediate    | `0.35 ≤ freq_g1 ≤ 0.65` AND extremes for g0/g2 |
| **Internal diagnostic**       | SNP/indel inside the candidate interval, no breakpoint claim required | inside `[start_bp, end_bp]` |
| **Breakpoint-proximal**       | Marker close to the inferred breakpoint zone (when available) | requires `boundaries_refined` from phase 5 |

The first three are always attempted; the fourth is opportunistic and
only fires when phase 5 has refined boundaries. None of these is
"breakpoint PCR" in the strict sense — even the breakpoint-proximal class
is just a regime-correlated marker that happens to sit near the
breakpoint.

### Panel construction rules

A marker is **not** an assay. A panel is. The minimum and recommended
panel sizes are:

- **Minimum panel**: 3–5 markers per candidate
- **Recommended panel**: 5–10 markers per candidate
- **Multiplex constraint**: per-panel Tm spread ≤ 4 °C, no primer-dimer
  hits at default tolerance, similar amplicon size buckets (small / mid)

A single marker can fail in the field (local mutation, family
background, PCR failure). Panels of 3+ markers tolerate one failure and
return a probabilistic regime call rather than a hard call. The panel's
output is:

```
predicted_regime + confidence
```

### Module placement (phase 13)

The marker module is **phase 13** — runs after phase 12 (gene cargo)
because gene cargo informs which markers to avoid (don't design primers
inside critical gene exons or repetitive regions). Phase 13 reads
`fish_regime_calls.tsv` (schema §7), the source VCFs, and the gene_cargo
JSON. It writes three TSVs and one `<chrom>_phase13_markers.json`.

### Three exported TSVs

#### `<chrom>_candidate_marker_catalogue.tsv`

One row per (candidate × marker). The "marker discovery" output —
all variants that passed filtering as plausible markers, with their
allele frequencies in each regime.

| Column                | Type    | Description                                                                                       |
|-----------------------|---------|---------------------------------------------------------------------------------------------------|
| `candidate_id`        | string  | Joins to `candidate_registry.tsv`                                                                  |
| `marker_id`           | string  | Stable ID (e.g. `cand_lg28_142_a_M001`)                                                          |
| `chr`                 | string  | Chromosome                                                                                        |
| `pos`                 | integer | 1-based position                                                                                  |
| `variant_type`        | string  | `SNP` / `INDEL` / `MNP`                                                                          |
| `ref`                 | string  | Reference allele                                                                                 |
| `alt`                 | string  | Alt allele (or pipe-joined alts)                                                                |
| `marker_class`        | string  | `regime_specific` / `het_pattern` / `internal_diagnostic` / `breakpoint_proximal`               |
| `target_regime`       | string  | `g0` / `g1` / `g2` / `g0_or_g2` (the regime(s) this marker preferentially tags)                  |
| `freq_g0`             | float   | Alt-allele frequency in fish called g0 by `fish_regime_calls.tsv`                                |
| `freq_g1`             | float   | Alt-allele frequency in fish called g1                                                           |
| `freq_g2`             | float   | Alt-allele frequency in fish called g2                                                           |
| `n_g0`                | integer | Number of fish (called g0) genotyped at this site (excludes ambiguous + missing)                 |
| `n_g1`                | integer | Same for g1                                                                                      |
| `n_g2`                | integer | Same for g2                                                                                      |
| `specificity_score`   | float   | `freq_target − max(freq_others)`. Range [-1, 1]. Higher = more discriminative.                   |
| `family_spread_score` | float   | Fraction of breeding families that contain ≥1 fish carrying the alt allele. Penalizes family-private markers. |
| `locus_quality_score` | float   | Composite of read-depth quartile, mappability, and proximity to repeats. Range [0, 1].          |
| `selected_for_panel`  | string  | `TRUE` / `FALSE` — passed all filters AND chosen by panel-assembly step                         |
| `panel_id`            | string  | The panel this marker belongs to (e.g. `cand_lg28_142_a_panel1`); empty if not selected         |

**Specificity rule.** A marker is *discriminative* when
`specificity_score ≥ 0.50` for `regime_specific` class (the alt allele's
frequency in the target regime is at least 0.50 higher than in the next-
highest regime). Lower thresholds may be admitted for `het_pattern` —
those are evaluated by absolute distance from 0.50 in the heterozygous
regime, not by max-of-others margin.

**Family-spread rule.** `family_spread_score ≥ 0.50` (allele present in
at least half of contributing families) avoids family-private markers
that look discriminative on the cohort but won't generalize. This rule
is critical for hatchery cohorts where breeding structure can fake
regime-specificity.

#### `<chrom>_candidate_marker_primers.tsv`

One row per primer pair. A marker can have multiple primer-pair
attempts (different windows around the variant); the best is selected
into a panel.

| Column                  | Type    | Description                                                                                       |
|-------------------------|---------|---------------------------------------------------------------------------------------------------|
| `candidate_id`          | string  | Joins to `candidate_registry.tsv`                                                                  |
| `marker_id`             | string  | Joins to `candidate_marker_catalogue.tsv`                                                          |
| `primer_pair_id`        | string  | Stable ID (e.g. `cand_lg28_142_a_M001_P1`)                                                       |
| `forward_primer`        | string  | 5'→3' DNA sequence                                                                                |
| `reverse_primer`        | string  | 5'→3' DNA sequence                                                                                |
| `forward_pos`           | integer | Forward primer's 5' bp position                                                                   |
| `reverse_pos`           | integer | Reverse primer's 5' bp position                                                                   |
| `tm_forward`            | float   | Predicted melting temperature, °C                                                                 |
| `tm_reverse`            | float   | Predicted melting temperature, °C                                                                 |
| `gc_forward`            | float   | GC fraction                                                                                       |
| `gc_reverse`            | float   | GC fraction                                                                                       |
| `amplicon_size`         | integer | Expected product size in bp                                                                       |
| `amplicon_seq`          | string  | Reference sequence of the amplicon (uppercase)                                                    |
| `target_regime`         | string  | Same as marker's `target_regime`                                                                  |
| `panel_id`              | string  | Panel this primer pair was assigned to                                                            |
| `selected`              | string  | `TRUE` / `FALSE` — chosen as the active pair for this marker in the panel                        |
| `primer_dimer_score`    | float   | Higher = more dimer risk; pairs with score ≥ 8 are rejected                                      |
| `off_target_hits`       | integer | Count of in-silico off-target hits across the genome (≤ 3 mismatches)                            |
| `warning`               | string  | Pipe-separated warnings: `repeat_region`, `near_gene_exon`, `low_mappability`, `tm_mismatch`, `large_amplicon`, ... Empty when clean. |

Empty `warning` is the goal. A primer pair with any warning is still
listed but is unlikely to be `selected = TRUE` unless no clean
alternative exists.

#### `<chrom>_candidate_marker_panel_summary.tsv`

One row per panel — the rollup the hatchery actually uses.

| Column                     | Type    | Description                                                                                       |
|----------------------------|---------|---------------------------------------------------------------------------------------------------|
| `candidate_id`             | string  | Joins to `candidate_registry.tsv`                                                                  |
| `panel_id`                 | string  | Stable ID (e.g. `cand_lg28_142_a_panel1`)                                                        |
| `n_markers`                | integer | Number of markers in the panel (3–10)                                                            |
| `n_g0_diagnostic`          | integer | Markers tagging g0                                                                                |
| `n_g1_diagnostic`          | integer | Markers tagging g1 (heterozygote-pattern)                                                         |
| `n_g2_diagnostic`          | integer | Markers tagging g2                                                                                |
| `panel_class`              | string  | `minimum` (3–4 markers) / `recommended` (5–10) / `extended` (>10, rare)                          |
| `tm_min`                   | float   | Minimum primer Tm in panel                                                                        |
| `tm_max`                   | float   | Maximum primer Tm                                                                                 |
| `tm_spread`                | float   | `tm_max − tm_min`. Multiplex-friendly when ≤ 4 °C.                                              |
| `amplicon_size_min`        | integer | Smallest amplicon, bp                                                                             |
| `amplicon_size_max`        | integer | Largest amplicon, bp                                                                              |
| `panel_specificity_mean`   | float   | Mean of selected markers' `specificity_score`                                                     |
| `panel_specificity_min`    | float   | Worst marker's specificity                                                                        |
| `expected_call_accuracy`   | float   | Probability the panel reproduces the genome-based regime call. Computed by leave-one-fish-out validation against `fish_regime_calls.tsv`. |
| `confidence_tier`          | string  | `HIGH` (≥ 0.95) / `MEDIUM` (0.85 – 0.95) / `LOW` (< 0.85)                                       |
| `warnings`                 | string  | Panel-level warnings: `low_family_coverage`, `single_diagnostic_for_regime`, `tm_spread_high`, `multiplex_unsafe`, `no_g1_marker` |
| `recommended_for_export`   | string  | `TRUE` / `FALSE` — panel is hatchery-deployable as written                                       |

**Confidence tiers.** `expected_call_accuracy` is the clearest
hatchery-facing number — *given the cohort we calibrated on, how often
does the panel agree with the genome-based regime call?*. Computed by
leave-one-fish-out: drop fish *i*, predict its regime from its alt-
allele dosage on the panel's markers (using the remaining cohort's
allele frequencies), compare to the genome call. Repeat for all fish.
The fraction that match is the tier number.

**Why not validate against true karyotype.** True karyotype is
unavailable for the 226-sample cohort (no cytology). The genome-based
regime call (`fish_regime_calls.tsv`) is the ground truth for panel
validation. This is documented honestly in the manuscript: panels
reproduce the genome-based regime call at confidence X, not the
underlying biological karyotype.

### Page 7 — marker page (planned scrubber UI)

The scrubber gains a new tab/page that activates when
`marker_panel_summary` is in `_layers_present`. Per-candidate view shows:

- **Header**: candidate ID, chr / start / end / span_mb, regime counts
  from `candidate_registry`
- **Recommended panel card**: panel ID, marker count, confidence tier,
  expected call accuracy, panel-level warnings
- **Marker table**: per-marker rows with marker_id, position, target
  regime, specificity score, primer pair (forward / reverse), amplicon
  size, Tm, warning chips
- **Interpretation block**: human-readable rules for the panel — e.g.
  *"Marker A supports g0; Marker B supports g2; Marker C separates
  g0/g2 from g1. Combined panel assigns a fish to a candidate regime
  with 96% accuracy on the 226-sample calibration cohort."*
- **Export buttons**: download the three TSVs filtered to this candidate
  (handy for sending one panel to a collaborator without exporting the
  whole atlas)

### Live dosage heatmap rendering (planned, schema 2.5 → revised 2.6)

A per-marker dosage heatmap is part of the §10 module's view.
**v2.6 substantially extends the v2.5 spec** with the FIG_C08-style
visual contract from `STEP28_candidate_marker_heatmap.R` (the existing
R-side reference plot), plus a region-binding contract that supports
both static (per-candidate) and live (cursor-following) views.

Naming convention: this is the "dosage heatmap"; in R-side and
manuscript contexts it is **FIG_C08**.

#### Visual layout (FIG_C08-style, expanded)

The heatmap is a Canvas matrix of **samples (rows) × markers (columns)**.
Cells are colored by per-sample dosage at that marker on a 3-stop
diverging ramp:

```
0 (HOMO_REF)  → #2166AC  (blue)
1 (HET)       → #F7F7F7  (near-white)
2 (HOMO_ALT)  → #B2182B  (red)
NA            → #FFFFFF  (white) or no fill (transparent)
```

Sample order within the matrix is `coarse_group_refined, u` — group
block first (so HOMO_1 / HET / HOMO_2 form clean horizontal bands),
then by anchor-PC1 score within group. v3.94 may approximate `u` with
raw PC1 within band when the anchor-rotated geometry from STEP20b is
unavailable (acceptable as v0; gives qualitatively the same ordering).

#### Within-regime subclustering (schema 2.10, post-v3.99)

When sub-structure exists *inside* a major regime — e.g. two distinct
HOMO_1 haplotype backgrounds carrying different ancestral arrangements
or recombination histories — running density-based clustering naively
on raw `(u, v)` across all samples conflates this signal with the
top-level HOMO_1 / HET / HOMO_2 split. The major-regime axis dominates
distances and the minor within-regime axis gets washed out.

**Contract: subclustering is per-regime, in a within-regime rotated
frame.** Per regime separately:

1. Extract local PCA coordinates `(u, v)` for samples in that regime.
   Source: `cand.fish_calls[si]` for K=3 candidates (regime-aware), or
   the per-window `cl.labels[]` for cursor-mode regions outside any
   committed candidate.
2. **Rotate.** Run PCA / SVD on the within-regime `(u, v)` matrix.
   Project each sample onto the rotated axes:
   - `u_rot` = projection on the within-regime PC1 (the dominant axis
     of variation *inside* this regime)
   - `v_rot` = projection on the within-regime PC2 (orthogonal)
3. **Standardize.** Z-score `u_rot` and `v_rot` independently to unit
   variance so DBSCAN's `eps` parameter behaves the same way along
   both axes. Mean-centering is implicit in the SVD step.
4. **Cluster.** Run DBSCAN (or HDBSCAN when `min_cluster_size` is more
   useful than `eps`) on the standardized `(u_rot, v_rot)`. Recommended
   defaults — DBSCAN: `eps = 0.5`, `min_samples = 3`. HDBSCAN:
   `min_cluster_size = 3`, `min_samples = 1`.
5. **Label.** Subcluster IDs are within-regime, dense, 1-based:
   `HOMO_1_sub1`, `HOMO_1_sub2`, ..., `HET_sub1`, `HET_sub2`, ...,
   `HOMO_2_sub1`, `HOMO_2_sub2`, ... DBSCAN noise points (label `-1`)
   become `<REGIME>_sub_noise` (e.g. `HOMO_1_sub_noise`).

**Subclusters are DIAGNOSTIC, not paper-level.** They are within-regime
haplotype backgrounds — useful for spotting recombinants, hybrid
ancestry, or sub-broodline structure inside what looks like a clean
homozygous block. The paper-level inversion regime call remains
`HOMO_1 / HET / HOMO_2`. Manuscript text should never report
"`HOMO_1_sub2` carries 12 samples" as a primary finding — that is a
description of within-regime broodline structure, not of the inversion
itself. Subclusters surface as a refinement annotation track in the
heatmap (between K=3 PCA group and K=6 sub-band), not as an
inversion-state label anywhere.

**`cluster_mode` parameter.** Four values are defined. The scrubber UI
exposes a primary toggle between the two production modes (`none` ↔
`rotated_within_regime`); the other two are named diagnostic
alternatives that exist in the contract but are not primary user
toggles.

| Mode | Default? | Description | When to use |
|---|---|---|---|
| `none` | **default** | Subclustering is not run. `coarse_group_refined == coarse_group` (no `_sub<N>` suffix). Heatmap K=3 track renders its existing solid-color cells. `u_rot`/`v_rot` are null. | Default until within-regime subclustering has been calibrated against real cohort data. Matches the current scrubber behaviour (no DBSCAN runs anywhere). |
| `rotated_within_regime` | opt-in | The 5-step contract above. Per-regime split → SVD-rotate `(u, v)` → standardize → DBSCAN/HDBSCAN. | The standard subclustering mode once calibrated. Surfaces within-regime haplotype backgrounds while preserving the paper-level regime call. |
| `raw_uv` | diagnostic alternative | DBSCAN on raw `(u, v)` across all samples (no per-regime split, no rotation). Equivalent to "ignore the major-regime split entirely". | Diagnostic only — to see what naive clustering would have produced. Almost always recovers the major-regime split as the dominant clusters and misses within-regime structure. Not a recommended production mode. |
| `one_dimensional_rotated_axis` | diagnostic alternative | Per-regime split → SVD-rotate, but cluster only on `u_rot` (1D). Discard `v_rot`. | When `v_rot` is dominated by sequencing noise rather than haplotype-background signal — produces cleaner clusters at the cost of losing any 2D structure. Useful as a fallback when low coverage makes the second within-regime axis non-informative. |

**Why `none` is the default.** Subclustering is a *new* analytical
contract introduced in schema 2.10. Until it has been calibrated on
the 226-sample LG28 cohort (and any per-regime sample-size lower
bounds tuned), the safe default is to leave the existing heatmap
behaviour unchanged: rows are sorted by `coarse_group_refined, u`
where `coarse_group_refined == coarse_group`, no sub-stripe renders,
no labels ship in the TSV. Once validated, the default may flip to
`rotated_within_regime` in a future schema bump — that change is
explicitly out of scope for 2.10.

**JSON shape.** Subcluster labels ship in the
`candidate_sample_coherence` layer as an additional column
`coarse_group_refined` (already promised by the v3.94 contract;
schema 2.10 makes the value space explicit):

```jsonc
{
  "candidate_id": 48,
  "sample": "CGA037",
  "coarse_group": "HOMO_1",                    // unchanged: paper-level regime
  "coarse_group_refined": "HOMO_1_sub2",       // NEW: within-regime subcluster (== coarse_group when cluster_mode=='none')
  "cluster_mode": "rotated_within_regime",     // NEW: which mode produced the label ('none' | 'rotated_within_regime' | 'raw_uv' | 'one_dimensional_rotated_axis')
  "u_rot": -0.83,                              // NEW: within-regime rotated coords (z-scored); null if cluster_mode in {'none', 'raw_uv'}
  "v_rot":  0.41,                              // NEW; null if cluster_mode in {'none', 'raw_uv', 'one_dimensional_rotated_axis'}
  "agreement_fraction": 0.83,
  "coherence_class": "coherent",
  "dist_to_centroid": 0.42,
  "centroid_z_score": 1.1,
  "stripe_quality": "core",
  "tier_rule": "coherence=coherent;z=1.1;group=HOMO_1",
  "n_informative_markers": 47
}
```

When `cluster_mode == "none"` (the default — subclustering not run),
`coarse_group_refined` falls back to `coarse_group` (no `_sub<N>`
suffix), `u_rot`/`v_rot` are null, and the `cluster_mode` column
records `"none"` so downstream consumers can tell "subclustering was
skipped" from "subclustering ran and produced this label". The dosage
heatmap row sort (`coarse_group_refined, u`) still works — it's just
that the within-regime ordering uses anchor-PC1 `u` instead of `u_rot`.

When `cluster_mode == "raw_uv"`, `coarse_group_refined` carries the
raw-clustering label (e.g. `cluster_3`) and `u_rot`/`v_rot` are null —
the rotation step was skipped.

When `cluster_mode == "one_dimensional_rotated_axis"`, `u_rot` is
populated; `v_rot` is null.

When `cluster_mode == "rotated_within_regime"`, both `u_rot` and
`v_rot` are populated.

**Heatmap rendering.** When subcluster labels are present *and*
`cluster_mode != "none"` *and* `cluster_mode != "raw_uv"` (i.e. only
under `rotated_within_regime` or `one_dimensional_rotated_axis`), the
K=3 PCA group track gains an inner sub-stripe showing the subcluster
assignment (one sub-color per `_sub<N>`, monotone within each major
regime's hue). Otherwise — labels absent, `cluster_mode == "none"`,
or `cluster_mode == "raw_uv"` (which produces global-clustering labels
not within-regime sub-labels) — the K=3 track renders its existing
solid-color cells unchanged.

**Recompute trigger.** Subcluster labels are computed at candidate
commit time (same as `band_diagnostics` and `het_shape`) and stored on
the candidate. Changing `cluster_mode` invalidates the stored labels;
the candidate page surfaces a "recompute subclusters" button.

#### Annotation tracks

**Five left-side annotation tracks**, drawn outermost (left edge of
canvas) to innermost (adjacent to dosage matrix):

| Track | Source | Cell rendering |
|---|---|---|
| **K=3 PCA group** | `cl.labels` per L2 OR `cand.fish_calls[si].regime` per candidate | Solid color from `groupColor(k)` (existing v3 palette: g0 / g1 / g2). |
| **K=6 sub-band** | `cand.k6_substructure.parent_of_k6` + per-fish K=6 path | Solid color from K=6 palette. Empty (background) when region does not overlap an active candidate or no K=6 pass for the candidate. |
| **Quality** (`stripe_quality`) | `candidate_sample_coherence` layer if loaded; otherwise scrubber-computed via the v3.94 button (see "Stripe-quality button" below) | Categorical: `core` `#1B7837` (green), `peripheral` `#F4A582` (peach), `junk` `#D73027` (red), `unknown` `#BDBDBD` (grey). Empty when no source. |
| **GHSL** | `ghsl_panel`, mean across all available scales over the visible region | Continuous viridis ramp. v3.94 ships this and `Het` reading from the same primary-scale matrix; the two diverge automatically when a dedicated multi-scale GHSL aggregate or a separate het matrix ships. |
| **Het** | `ghsl_panel` primary scale only | Continuous viridis ramp, distinct visual treatment (slightly different ramp endpoints or single-color saturation). |

**One top-side annotation track**:

| Track | Source | Cell rendering |
|---|---|---|
| **Polarity** (`final_flip_decision`) | `candidate_marker_polarity` layer if loaded | Binary: `FALSE` `#2166AC` (blue, "matches dominant orientation"), `TRUE` `#D73027` (red, "marker flipped"). Empty when polarity layer not loaded — top track simply doesn't render. |

**Why two GHSL/het tracks today when they share a source.** Two reasons:
(1) future-proofing — when a dedicated het matrix or a multi-scale GHSL
aggregate ships, the two tracks already exist as separate visual
channels and just diverge automatically; (2) the GHSL primary scale by
construction *is* phased-SNP heterozygosity, so calling them the same
thing today is honest, but operationally they answer different
questions (GHSL = "how diverged is this sample's haplotype block from
others?", het = "how heterozygous is this sample at this locus?").
Surfaced as two tracks now to avoid a confusing rename later.

#### Region binding — two view modes share one renderer

The heatmap's region (`{ start_bp, end_bp }`) is one of two sources:

**(a) Candidate-bound (static).** Region = `cand.start_bp..cand.end_bp`.
Recomputed only when the active candidate changes. Equivalent to the R
reference (one heatmap per candidate, full candidate span).

**(b) Cursor-bound (live).** Region = the bp range covered by the 11
local-PCA windows centered on `state.cur` (i.e. `state.cur ± 5`,
inclusive). Recomputed every time `state.cur` changes — including when
arrow keys advance the L3 contingency carousel, since the carousel and
the scrubber share `state.cur` already. The window radius is **fixed
at ±5** in v3.94 (not user-adjustable yet); this matches the existing
`stepMode === 'win5'` arrow-key step.

**Marker cap is independent of region source.** Whether region comes
from (a) or (b), the marker-selection cascade applies the same way
(see "Marker selection within the visible region" below). Region
defines *where*, the cap defines *how many*.

#### Marker selection within the visible region

When the region contains more markers than the user-selected display
cap, the scrubber selects which to show, in this priority order:

1. Collect every marker whose position lies inside the region.
2. Drop markers exceeding a missingness threshold (default `0.20`).
3. Rank the remaining markers by `diagnostic_score` from
   `marker_catalogue.tsv` if loaded.
4. Otherwise, rank by per-marker dosage variance across the cohort
   (computed once per chunk, cached).
5. Take the top `N` markers, where `N` is the user-selected display
   cap. Default `N = 200`; toolbar offers `100 / 200 / 500` (hard cap
   500 — no option for more).

If the region has fewer than `N` markers passing the missingness
filter, show all available — never synthesize, never duplicate.

#### Stripe-quality button

The Quality track requires per-(candidate × sample) `stripe_quality`
classification (`core` / `peripheral` / `junk` / `unknown`). Two
sources:

1. **`candidate_sample_coherence` layer loaded** (preferred). Track
   reads directly from this layer's `stripe_quality` column.
2. **Layer not loaded.** A button on the heatmap toolbar — "Compute
   stripe quality (~1 s)" — runs the scrubber-side equivalent of
   STEP29's coherence algorithm on the currently-loaded dosage chunk.
   Result is cached on `cand.__stripeQuality` for the active
   candidate. Re-runs when the candidate or chunk changes. Computed
   silently for the active candidate; not pre-computed for all
   candidates.

The scrubber-side coherence algorithm is a faithful port of STEP29
lines 200–278: agreement_fraction with group-specific thresholds (HET
uses `0.55 / 0.40` cuts because expected intermediate; HOMO_1/HOMO_2
use `0.70 / 0.45`); `centroid_z` from PC1 score relative to the band's
centroid in MAD units; tier from the conjunction of coherence_class +
centroid_z. Implementation note: the algorithm is O(n_samples ×
n_informative_markers), where n_informative_markers is the top quartile
by |Δ_HOMO|. On 226 samples × ~50 informative markers this runs in
single-digit ms — the "~1 s" message accounts for chunk fetch + DOM
work, not the algorithm itself.

#### Placement (deliberately undecided)

The heatmap's home in the scrubber UI is **not specified by this
schema**. v3.94 will decide between:

- **Popup** (modal-like overlay, dismissable, similar to the v3.93 het
  popover but larger)
- **Docked strip** below the L3 contingency row on page 1 (always-on,
  follows the cursor when in cursor-bound mode)
- **Floating panel** (resizable, draggable, persistent across page
  changes)
- **Section on candidate page 2** for the static view + a separate
  placement for the live view

This is a UI decision that depends on what feels right when seen on
real data. The schema fixes the **visual contract** (tracks, colors,
sample order, region binding) but explicitly leaves the **container**
to implementation.

#### Redraw triggers

The Canvas heatmap MUST redraw on, and only on:

1. User changes `K` (re-clusters → row order changes).
2. User changes the row or column sort.
3. Region changes (candidate switch in static mode, `state.cur`
   change in live mode).
4. User changes display cap `N`.
5. Dosage chunk for a newly-visible region finishes loading.
6. `candidate_sample_coherence` layer finishes loading (Quality track
   appears).
7. `candidate_marker_polarity` layer finishes loading (Polarity top
   track appears).
8. Stripe-quality button completes (Quality track populates from
   scrubber-side compute).

Avoid redraw on hover, tooltip events, or any selection that doesn't
change the marker set or sort order. Hover affordances (e.g. "sample
CGA037, marker M0124, dosage 1") MUST be implemented as overlay
updates on the cached Canvas image data, never as full redraws.

In **cursor-bound mode**, redraws are debounced **150 ms** after the
last `state.cur` change to avoid thrashing during rapid arrow-key
pressing or scrubber drag.

#### Dosage chunk file shape (R emit shipped v3.95)

The R emit script `STEP_M04_emit_dosage_chunks.R` (shipped v3.95) writes
per-region chunks under a predictable URL pattern:

```
<chrom>_dosage_chunk_<region_start_bp>_<region_end_bp>.json
```

Chunks are aligned to L1 boundaries (so chunks line up with the
scrubber's natural windowing). Loaded on demand by the scrubber as the
region changes; LRU-cached in memory (cap ~50 chunks ≈ a few hundred
MB max for the 226-sample cohort).

Each chunk JSON:

```jsonc
{
  "chrom": "C_gar_LG28",
  "region_start_bp": 12000000,
  "region_end_bp": 14000000,
  "n_markers": 412,
  "n_samples": 226,
  "samples": ["CGA001", "CGA002", ...],
  "markers": [
    { "marker_id": "M0001", "pos_bp": 12018473,
      "diagnostic_score": 0.92, "missingness": 0.04 },
    ...
  ],
  "dosage": [
    [0, 1, 2, 0, 1, ...],   // marker 0, length n_samples (0/1/2 or -1 = NA)
    [1, 1, 0, 2, 1, ...],
    ...
  ]
}
```

#### `candidate_sample_coherence` layer (R emit shipped v3.95)

Per-(candidate × sample) row, faithful mirror of STEP29's
`candidate_sample_coherence.tsv` output. Emitted by
`STEP_M05_emit_step29_layers.R`:

```jsonc
{
  "candidate_sample_coherence": [
    {
      "candidate_id": 48,
      "sample": "CGA037",
      "coarse_group": "HOMO_1",          // 'HOMO_1' | 'HET' | 'HOMO_2'
      "agreement_fraction": 0.83,        // 0..1
      "coherence_class": "coherent",     // 'coherent' | 'intermediate' | 'discordant' | 'insufficient'
      "dist_to_centroid": 0.42,
      "centroid_z_score": 1.1,
      "stripe_quality": "core",          // 'core' | 'peripheral' | 'junk' | 'unknown'
      "tier_rule": "coherence=coherent;z=1.1;group=HOMO_1",
      "n_informative_markers": 47
    },
    ...
  ]
}
```

Detection: presence of `data.candidate_sample_coherence` as an array.
Layer name: `candidate_sample_coherence`. Activates the Quality track
in heatmap mode (a) static — the live mode (b) requires per-window
recompute against the live region's dosage, not stored coherence, so
it always uses the button-computed path.

#### `candidate_marker_polarity` layer (R emit shipped v3.95)

Per-(candidate × marker) row, faithful mirror of STEP29's
`candidate_marker_polarity.tsv` output. Emitted by
`STEP_M05_emit_step29_layers.R`:

```jsonc
{
  "candidate_marker_polarity": [
    {
      "candidate_id": 48,
      "chrom": "C_gar_LG28",
      "pos": 13782104,
      "marker_id": "M0124",
      "raw_mean_HOMO_1": 0.21,
      "raw_mean_HOMO_2": 1.78,
      "raw_mean_HET": 0.95,
      "delta_hom_mean": 1.57,
      "abs_delta_hom": 1.57,
      "polarity_L1_sign": 1,
      "polarity_L1_confidence": 1.0,
      "polarity_L1_ambiguous": false,
      "polarity_L1_reversed": false,
      "block_id": 7,
      "polarity_L2_sign": 1,
      "polarity_L2_confidence": 0.95,
      "polarity_L2_reversed": false,
      "final_flip_decision": false,
      "flip_source": "L2_block",          // 'L2_block' | 'L1_marker'
      "support_class": "strong"           // 'strong' | 'moderate' | 'weak' | 'ambiguous'
    },
    ...
  ]
}
```

Detection: presence of `data.candidate_marker_polarity` as an array.
Layer name: `candidate_marker_polarity`. Activates the top Polarity
track when present.

#### Failure modes (must handle gracefully)

- **Region with zero markers** → empty heatmap with a "no markers in
  this region" placeholder. No error.
- **Dosage chunk fetch fails** → toolbar shows "(loading failed)";
  rest of UI unaffected.
- **Scrub faster than chunks load** → 150 ms debounce; cancel
  in-flight fetches no longer for the current region.
- **Display cap larger than markers in region** → show all available.
- **Stripe-quality button pressed but no dosage chunk loaded yet** →
  button disabled with hint "load dosage first".
- **Active candidate switches while heatmap rendered in cursor-bound
  mode** → K=6 sub-band track may change content; redraw normally.
- **Cursor scrubs outside any active candidate** in cursor-bound mode
  → K=6 sub-band track shows empty / background; other tracks render
  normally; no error.

#### Why the specific numbers

- **±5 windows** ≈ 50 kb to 1 Mb of genome at typical local-PCA window
  sizes; matches the existing `stepMode === 'win5'` arrow step. Wider
  blurs fine structure; narrower lacks enough markers to fill the
  matrix.
- **200 default markers** × 226 samples at typical pane width gives
  ~3-px cells — readable per cell, no zoom needed.
- **500 hard cap** is the empirical break-point where K-change redraws
  start dropping frames on the 13" laptop hardware the scrubber runs
  on routinely.
- **150 ms debounce** matches the v3.91 click-debounce idiom and is
  short enough not to feel laggy during arrow-key scrubbing.

This contract is forward spec — the heatmap component does not yet
exist in the scrubber. v3.94 will implement it. The contract is
sufficient for the R emit script (`STEP_M04_emit_dosage_chunks.R`)
and the scrubber-side renderer to be built in parallel without
further design loops.

#### Hover affordances (v3.96, schema 2.7)

A one-line tooltip appears next to the cursor when hovering any cell
of the dosage matrix. Tooltip MUST be implemented as overlay updates
on the cached Canvas image data, never as full redraws. The renderer
already returns the layout fields needed for hit-testing
(`{ matX, matY, matW, matH, cellW, cellH }`); the tooltip pipeline
consumes those.

**Format.** Components space-`·`-space separated; absent fields
omitted (no `null` placeholders, no double separators):

```
sample · marker (bp) · dosage=N · group · quality · flipped
```

Example with all fields populated:

```
CGA037 · M0124 (12,018,473 bp) · dosage=1 · HOMO_1 · core
```

The `flipped` segment appears only when `final_flip_decision === true`.
The `quality` segment is omitted when value is `unknown` or absent.
The `group` segment is omitted when sample-group lookup returns null
(e.g. cursor-bound bare-genome mode without an active candidate).
NA dosage cells render `dosage=NA`; the `-1` sentinel encoding is
never exposed in user-facing strings.

**Annotation tracks** (left + top) do not currently tooltip — only the
dosage matrix proper. Hit-tests outside the matrix's `[matX, matY,
matX+matW, matY+matH)` half-open rectangle return null and the tooltip
hides. This is the v3.96 behaviour; annotation-track tooltips may be
added in a later increment if requested.

**Performance contract.**

- Cell-stable suppression: when the cursor stays in the same cell
  across multiple `mousemove` events, no DOM writes occur. The check
  is on `(marker_idx, sample_idx)` equality.
- RAF coalescing: at most one DOM write per frame. Multiple
  `mousemove` events between frames coalesce into a single tooltip
  update via `requestAnimationFrame`.
- Tooltip element is reused (one per canvas, lazily created on first
  mount). `pointer-events: none` on the tooltip prevents the cursor
  from flickering between hit/no-hit at cell boundaries that overlap
  the tooltip.

**Hover state lifecycle.**

The state object `state.__dosageHm` carries two hover contexts:
`hover_static` (page-2 candidate heatmap) and `hover_cursor` (L3
docked strip), each shaped `{ last_cell, raf_handle }`. Both reset
to `{ null, null }` whenever the heatmap redraws (the layout the
prior `last_cell` referred to is now stale).

**Position clamping.**

Default offset is +12 px right and +4 px below the cursor. When the
tooltip would overflow the canvas wrapper, it flips to left or above
respectively, then clamps to a 4 px margin from the wrapper bounds.

### Cross-references

- Phase 13 reads `fish_regime_calls.tsv` (schema §7) as ground truth for
  per-regime allele-frequency computation — the marker logic is meaningless
  without those calls.
- Phase 13 may consult `gene_cargo` (phase 12) to avoid designing primers
  inside annotated exons.
- Phase 13 may consult `boundaries_refined` (phase 5) to attempt
  breakpoint-proximal markers; falls back to internal-diagnostic when
  unavailable.
- The page activation in the scrubber is conditional on
  `marker_panel_summary` being in `_layers_present`. The other two layers
  are not strictly required for page rendering (the page can show
  "panel summary present, marker detail not loaded" if only the summary
  is dragged in).

---

## 11. Band diagnostics module (v3.92, schema 2.3)

For each PCA band/regime in the focal L2 envelope (or any contiguous
interval), the band diagnostics module computes per-band summary
statistics across up to five data layers. The result is a
**confounder/support check**: it tells the user whether a band is likely
a real structural haplotype regime, or whether it might be influenced
by ROH/inbreeding, unusual heterozygosity, local diversity, or GHSL
haplotype-divergence signal.

### Why "diagnostics", not "validation"

A band that fails one of these checks is not automatically wrong — but
it is worth pausing on. Real inversion regimes typically:

- Show some GHSL spread between bands (one arrangement carries longer
  haplotype blocks than the other).
- Show a θ/π shift if the regimes have different effective population
  sizes, recombination histories, or sample-size imbalances.
- Show elevated dosage heterozygosity in the **middle** band of a K=3
  fit when the inversion is heterozygous-on-the-middle (the classic
  inversion-Hardy-Weinberg pattern).
- Show no excess ROH overlap in any one band — ROH is a per-sample
  property, not a structural one, so a band with much higher ROH
  overlap suggests inbreeding-driven LD, not inversion-driven LD.

The diagnostics panel makes those expectations visible per candidate.

### Two render paths sharing one compute

`computeBandDiagnostics(cl, env, l2idx)` is the single source of truth.
It returns a structured object that drives both render paths:

- **Mini-chips** (always-visible per L2 tile in the L3 contingency
  panel). Per band, a row of four small colored rectangles in the
  band's color, white text inside: `GHSL <value>`, `θπ <value>`,
  `het <value>`, `ROH <pct>`. Quick visual scan across bands without
  expanding anything.
- **Full panel** (`<details>`-collapsible per L2 tile, default
  collapsed). 14-column table with all 12 metrics + flags column.

### The 12 metrics per band

| Metric                          | Source layer       | Notes |
|---------------------------------|--------------------|-------|
| `n_samples`                     | cluster labels     | Per-band sample count from `cl.labels`. |
| `mean_GHSL`, `median_GHSL`      | `ghsl_panel`       | Per-sample mean over panel windows whose midpoint falls inside `env.start_bp..env.end_bp`, then mean/median across samples in the band. Uses `panel.primary_scale`. |
| `mean_theta_pi`, `median_theta_pi` | `theta_pi_panel` | Same shape and aggregation as GHSL. **Layer planned, not yet emitted by R.** |
| `mean_dosage_heterozygosity`, `median_dosage_heterozygosity` | `ghsl_panel` (same scale) | By construction in `STEP_C04_snake3_ghsl_v6.R`, the GHSL primary-scale matrix is per-sample phased-SNP heterozygosity rate. We treat het as a separate metric for clarity even though it currently shares the same source. |
| `roh_overlap_pct`               | `roh_intervals`    | % of band's samples with at least one ROH span overlapping the interval. **Layer planned.** |
| `mean_ROH_fraction_of_interval`, `median_ROH_fraction_of_interval` | `roh_intervals` | Per-sample fraction = clipped ROH bp inside interval / interval length, capped at 1. Samples with no overlapping ROH contribute 0 (not NaN). **Layer planned.** |
| `mean_FROH_of_samples_in_band`, `median_FROH_of_samples_in_band` | `sample_froh` | Genome-wide FROH per sample, mean/median across band's samples. Independent of the interval. **Layer planned.** |

When a source layer is missing, the corresponding metrics are `null`,
which the renderer displays as `?`. **Missing data does not break the
app and does not cause spurious flag firing** — flags are only
evaluated against metrics that are present.

### Six flags

| Flag                  | Fires when                                                                                | Source needed |
|-----------------------|-------------------------------------------------------------------------------------------|---------------|
| `het_high`            | One band's `het_mean` > 1.5 × max of other bands' `het_mean`                              | het           |
| `middle_het_support`  | K=3, middle PC1-rank band has highest het AND > 1.2 × both flanking bands                 | het           |
| `GHSL_support`        | Spread `max(ghsl_mean) / min(ghsl_mean)` > 1.3 across bands. Marked on the max band.      | GHSL          |
| `theta_pi_shift`      | Same logic, on `theta_pi_mean`. Marked on the max band.                                   | θ/π           |
| `ROH_confounded`      | One band's `roh_overlap_pct` > 2 × max of other bands AND > 5 % absolute floor            | ROH           |
| `FROH_high`           | One band's `froh_mean` > 1.5 × max of other bands                                         | FROH          |

Flags attach to the band they describe. Reading a band's row tells you
every flag that fires on that band.

### `<chrom>_band_diagnostics.tsv` — sibling TSV (v3.92)

One row per (candidate × band). Emitted alongside the candidate
registry / fish_regime / interval_support TSVs when at least one
candidate carries band diagnostics.

| Column                   | Type    | Notes |
|--------------------------|---------|-------|
| `candidate_id`           | string  | Matches `candidate_registry.id`. |
| `band`                   | string  | `g0`, `g1`, `g2`, ... (PC1-ranked).
| `n_samples`              | int     | Per-band sample count. |
| `ghsl_mean`              | float   | `NA` if GHSL panel not loaded. |
| `ghsl_median`            | float   | `NA` if missing. |
| `theta_pi_mean`          | float   | `NA` if θπ panel not loaded. |
| `theta_pi_median`        | float   | `NA` if missing. |
| `het_mean`               | float   | `NA` if GHSL panel not loaded (shares source). |
| `het_median`             | float   | `NA` if missing. |
| `roh_overlap_pct`        | float   | 0–100 with two decimal places, `NA` if ROH layer not loaded. |
| `roh_frac_mean`          | float   | 0–1, `NA` if missing. |
| `roh_frac_median`        | float   | 0–1, `NA` if missing. |
| `froh_mean`              | float   | 0–1 (genome-wide), `NA` if FROH layer not loaded. |
| `froh_median`            | float   | 0–1, `NA` if missing. |
| `flags`                  | string  | Pipe-separated flag names, e.g. `het_high\|GHSL_support`. Empty when no flags. |
| `has_ghsl`               | int     | 0 / 1. Reflects `data_status.ghsl` at compute time. Lets consumers tell missing-data from genuine zero. |
| `has_theta_pi`           | int     | 0 / 1. |
| `has_het`                | int     | 0 / 1. (Same as `has_ghsl` until θπ panel ships a separate het channel.) |
| `has_roh`                | int     | 0 / 1. |
| `has_froh`               | int     | 0 / 1. |

**Numeric format note.** Because the metrics span six orders of magnitude
(`theta_pi` ~ 1e-3, `roh_overlap_pct` ~ 1e1), the TSV uses
`Number.toExponential(4)` for all metric values except `roh_overlap_pct`
(plain `.toFixed(2)`). This keeps the TSV machine-parseable without per-metric
locale-specific formatting.

### `band_diagnostics` field on candidate (JSON-mirror)

The full `computeBandDiagnostics` return object (with `bands` array and
`data_status`) is attached to each candidate at commit time, computed
at the candidate's `ref_l2`. Round-trips intact through the
candidate-list JSON export/import path.

### What the diagnostic does NOT replace

- **It is not a structural call.** A band with `GHSL_support` flag is
  more likely to be a structural regime; a band without is not
  necessarily wrong. Use in conjunction with the K=6 nesting verdict
  (§8) and the band continuity (NESTED / MIXED / CROSS_CUTTING /
  UNSTABLE).
- **It is not a fish-level karyotype call.** That is `fish_regime_calls.tsv`
  (§7).
- **It is not a population-level diversity scan.** That is the rolling
  `tracks` block (phase 1+2). Band diagnostics restricts to one
  interval and asks, "do the bands behave plausibly here?"

### Why the source layers are still "planned"

Two of the four source layers wired today are real
(`ghsl_panel` is exists, dosage het rides on the same matrix); the other
three (`theta_pi_panel`, `roh_intervals`, `sample_froh`) are speced and
detected by the scrubber but no R emit script has been written. When
those scripts ship, the panel automatically lights up — no scrubber-side
change needed. The contract here is sufficient for a parallel-chat R
implementation when the user is ready.

---

## 11.1. Het-shape sub-module (v3.93, schema 2.4)

A single-number heterozygosity ratio hides the **shape** of the
distribution. Inversion regimes have qualitatively different
heterozygosity shapes — the heterozygote band shows a tight peak near
0.5 (or wherever the ploidy makes the expected het cluster), while
homozygote bands show broader low-het distributions sometimes with
outlier tails or bimodality from sub-systems. Capturing shape requires
a histogram, not just a mean.

### What the sub-module adds

`computeBandDiagnostics` was extended to return a new `het_shape` field
alongside `bands` and `data_status`. Renderer plumbing has three
placements that all consume this one structure.

### `het_shape` field on the diag

```jsonc
{
  "n_bins": 16,
  "bin_lo": 0.0,
  "bin_hi": 1.0,
  "bin_edges": [0.0, 0.0625, 0.125, ..., 1.0],   // length n_bins+1
  "counts_per_band": [
    [0, 2, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],   // band 0
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 1],   // band 1 (HET-like)
    [0, 0, 4, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]    // band 2
  ],
  "support_ratio": 5.667,
  "support_kind": "middle_vs_flanks",   // or 'max_min_spread' or 'none'
  "support_passes": true,                // ratio >= 1.5
  "support_marginal": false              // ratio in [1.2, 1.5)
}
```

### The support ratio

A single quick-glance summary number that rolls up the K-band
heterozygosity spread.

- **K = 3** uses `support_kind = "middle_vs_flanks"`:
  `median(g1) / mean(median(g0), median(g2))`. This is the
  inversion-Hardy-Weinberg test in disguise — does the middle PC1-rank
  band carry a heterozygosity peak relative to the homozygote flanks?
- **K ≠ 3 (or K=3 fallback when middle/flanks medians missing)** uses
  `support_kind = "max_min_spread"`: `max(median) / min(median)` across
  bands.
- **support_passes** at `ratio ≥ 1.5` (✓ glyph, strong signal)
- **support_marginal** at `1.2 ≤ ratio < 1.5` (?? glyph, weak signal)
- **No glyph** otherwise (or `none` support_kind when fewer than 2 bands
  have het data)

### Three placements

The renderer ships three views; all read from the same `het_shape`:

1. **Inline ridgeline column** (placement A) — small SVG sparkline
   (~80×14 px) in the band-diagnostics `<details>` table inside each
   L3 contingency tile. One row per band. Visible only when the user
   expands the panel. Tiny but enough to see narrow vs broad vs
   bimodal.
2. **FIG_C07-style stacked ridgeline figure** (placement B) — full
   readable figure on the candidate page (page 2), styled like the
   existing `FIG_C07 — Heterozygosity by group` plot. Per-band median
   tick on the baseline; X axis 0..1 with 0.0 / 0.25 / 0.5 / 0.75 /
   1.0 labels. Computed at the candidate's `ref_l2`. Falls back to
   on-the-fly recompute when a candidate predates v3.92 (no
   `band_diagnostics` saved on the object).
3. **Hover-popover from the het mini-chip** (placement C) — clicking
   any `het` pill in the L3 tile mini-chips row pops a floating panel
   with the same full ridgeline figure. Click-outside or `×` button
   dismisses. No persistent screen real-estate cost.

### One-line summary text

`Het: 5.4× ✓ · g0=0.15 | g1=0.84 | g2=0.16` — rendered above the
`<details>` table when the het layer is loaded, and inside the
candidate-page section. Same data the popover uses (no `<details>`
expand needed for the line).

### Shared SVG ridgeline kernel

`_hetRidgelineSvg(binCounts, color, width, height)` is the single
shape-drawing primitive. Builds an area path under the histogram +
top-edge stroke, with band-color fill at 0.55 opacity. Used at 80×14
for placement A, embedded at larger size inside `_hetRidgelineFigureHtml`
for placements B and C. Ratio of width to height is left to the caller —
inline cells stay tiny; figure rows are ~36 px tall with stacked
baselines.

### Histogram round-trip in TSV

The sibling TSV (`<chrom>_band_diagnostics.tsv`) gained 4 columns at the
end so downstream R code can re-render the ridgeline without recomputing
from raw GHSL panel data:

| Column                | Type   | Notes |
|-----------------------|--------|-------|
| `het_support_ratio`   | float  | Same value as `het_shape.support_ratio`. `NA` when not computed. |
| `het_support_kind`    | string | `middle_vs_flanks` \| `max_min_spread` \| `none`. |
| `het_support_passes`  | int    | 0 / 1. Reflects `support_passes`; `support_marginal` is not separately encoded — interpret a `support_passes=0` row alongside the ratio to recover marginal status (`1.2 ≤ ratio < 1.5`). |
| `het_hist_bins`       | string | Pipe-separated 16 integers, `0\|2\|0\|...` covering `[0, 1]` in 16 equal bins. Empty when het not loaded for this candidate. |

**Total TSV columns are now 24** (was 20 in v2.3).

### Numerical conventions

- 16 bins on `[0, 1]` is fixed for v3.93. The heterozygosity domain is
  bounded by construction (it's a rate), so a fixed extent works
  cleanly. Future versions may extend (e.g. 32 bins or auto-binning)
  if the panel resolution warrants — but downstream consumers should
  read `het_shape.n_bins` and `het_shape.bin_edges`, not assume 16.
- All bands share the same bin grid, so per-band counts are directly
  comparable and stackable.
- The Y-axis of the rendered ridgeline normalizes per-band by the
  band's own max bin count (placement A inline) or by the global max
  across bands (placement B / C figure) — the former optimizes for
  shape-comparison-within-row; the latter optimizes for sample-count
  comparison across bands. This is a deliberate dual choice; both have
  legitimate interpretive uses.

### Threshold rationale

The 1.5×/1.2× thresholds for `support_passes` and `support_marginal`
were chosen empirically from FIG_C07-equivalent plots in the
manuscript: candidates that read visually as "real inversion regimes"
sit comfortably above 1.5×, while genuinely ambiguous candidates fall
in the 1.2–1.5× band. This mirrors the existing `middle_het_support`
flag in §11 (which fires at K=3 middle band > 1.2× both flanks and is
preserved unchanged for backwards compatibility).

---

## 12. Boundary zone refinement (v3.97, schema 2.8)

After candidates are promoted on the Candidate Builder (page 2), the
**Boundaries page** (page 11 internally; displayed as tab "3
boundaries", wedged between candidate-focus and catalogue) refines each
promoted candidate's `[start_bp, end_bp]` interval into approximate
left and right *boundary zones* using multiple evidence tracks.

### Terminology contract — CRITICAL, locked

The schema enforces three terms with strict semantic boundaries.
Implementations MUST NOT conflate them.

| Term | Meaning | When used |
|---|---|---|
| **boundary zone** | An *interval* `[zone_start_bp, zone_end_bp]` of plausible breakpoint locations, derived from population-genetic signal transitions (PCA, dosage, FST, GHSL, etc.). | Default. Always emitted when boundary refinement runs. |
| **exact breakpoint** | A single bp position with sub-100bp resolution. | Only when junction-level evidence exists (split reads at a specific position, or an assembly junction). NEVER inferred from population statistics alone. |
| **SV anchor** | A DELLY/Manta/etc. SV call near the boundary zone. Adds confidence; does NOT alone justify "exact". | Bonus signal — when its bp position falls inside the auto-detected zone, confidence increases. When SV anchor is the *only* evidence beyond population signal, status is `SV_supported` not `junction_supported`. |

The status enum on each candidate encodes the strongest available
evidence:

```
breakpoint_status :=
  'boundary_zone_only'   -- population signal only; default for all auto-detected zones
  'SV_supported'         -- SV anchor (DELLY/Manta INV/BND/DEL) falls inside the zone
  'junction_supported'   -- split-read or assembly junction at sub-100bp resolution
```

The candidate-level field is the *minimum* of the two edges' statuses
(the candidate is `junction_supported` only if both edges are). Per-edge
status lives in `boundary_left.support_class` (orthogonal axis: how
much *track* support, not what *kind* of evidence).

### Why "boundary zone" not "exact breakpoint"

Local PCA, dosage, FST, GHSL all report population-genetic signals
that change *across* a transition zone of variable width — typically
5–50 kb depending on recombination suppression strength. They CANNOT
pinpoint a single bp. Reporting an "exact breakpoint" from these
signals would be a category error.

Exact breakpoints come from:

- Split-read alignment at the breakpoint (BWA-MEM soft-clipped reads
  with mate at the inversion's other end)
- Assembly junctions in long-read assemblies (PacBio HiFi, ONT)
- Sanger sequencing of PCR products spanning the junction

When that evidence ships, the schema upgrades `breakpoint_status` to
`junction_supported` and an additional field
`boundary_left.exact_breakpoint_bp` (integer) appears alongside the
zone. Until then, all `_bp` fields are zone-bounded.

This terminology discipline is enforced at two layers:

- TSV column names use `*_boundary_zone_*` (not `*_breakpoint_*`)
- UI strings use "boundary zone" not "breakpoint" — except in the
  `breakpoint_status` field name itself, which describes what *kind*
  of breakpoint evidence we have (including `boundary_zone_only` for
  "none beyond zone")

### `boundary_evidence` layer (optional R-side input — shipped via STEP_M06)

When R-side has precomputed boundary-relevant signals (FST, theta-pi,
discordant-pair pileups, SV anchors), they ship as a new layer.
**Emitted by `STEP_M06_emit_boundary_evidence.R`** (post-v3.98); the
v3.97 scrubber detects and consumes the layer when present, including
optional `theta_pi_het` (third regime track):

```jsonc
{
  "boundary_evidence": [
    {
      "candidate_id": 48,
      "chrom": "C_gar_LG28",
      "scan_start_bp": 11000000,        // candidate ± scan_radius
      "scan_end_bp":   15000000,
      "scan_window_bp": 5000,           // bp resolution of arrays below
      "tracks": {
        "fst": [0.02, 0.03, 0.05, 0.18, 0.21, ...],         // length = (scan_end-scan_start)/scan_window
        "theta_pi_homo1": [0.0014, 0.0013, ...],
        "theta_pi_het":   [0.0021, 0.0019, ...],            // optional; populated when g1 has >= 2 fish
        "theta_pi_homo2": [0.0011, 0.0012, ...],
        "discordant_pair_pileup": [3, 5, 12, 89, 102, ...],  // count per scan window
        "sv_anchors": [
          { "kind": "DELLY_INV", "ct": "3to3", "pos_bp": 12345670, "qual": 80 },
          { "kind": "Manta_INV5", "pos_bp": 14250130, "qual": 95 }
        ]
      }
    }
  ]
}
```

Detection: `Array.isArray(data.boundary_evidence)` AND first row has
`candidate_id` and `tracks`. Layer name: `boundary_evidence`. The
boundaries page works **without** this layer — it falls back to
in-scrubber scoring from existing layers (precomp, GHSL panel, dosage
chunks, candidate_marker_polarity). The layer is purely additive;
presence increases auto-propose confidence.

The `scan_window_bp` resolution is intentionally coarser than the
local-PCA window step (typical 5–10 kb) — boundary scoring benefits
from already-smoothed input that the scrubber doesn't re-blur further.

### Per-candidate boundary fields (extends `candidates_registry`)

Each promoted candidate gains nested boundary records on its JSON
representation (and flattened columns in the TSV mirror):

```jsonc
{
  "candidate_id": 48,
  // ... existing fields ...

  "boundary_left": {
    "zone_start_bp": 11820000,
    "zone_end_bp":   11860000,
    "score":         0.78,            // 0..1; weighted sum of contributing tracks
    "support":       ["pca_drop", "dosage_transition", "ghsl_step", "polarity_change"],
                                       // sorted, deduped; subset of TRACK_NAMES
    "support_class": "strong",        // 'strong'|'moderate'|'weak'|'ambiguous'
    "source":        "auto",          // 'auto'|'manual'|'auto_plus_manual'
    "sv_anchors_in_zone": [           // copies from boundary_evidence.tracks.sv_anchors
      { "kind": "DELLY_INV", "pos_bp": 11842300 }
    ],
    "notes": "",
    "set_at":  "2026-04-29T14:32:11+0700",
    "set_by":  "scrubber_auto"        // or 'scrubber_manual' for E/F overrides
  },

  "boundary_right": { /* same shape */ },

  "breakpoint_status": "boundary_zone_only",   // candidate-level min of both edges
  "boundary_notes":    ""                       // optional free text per candidate
}
```

**TSV mirror columns** (added to `REGISTRY_TSV_COLUMNS` after the
`subband_nesting_status` block, before `qc_status`):

```
'left_boundary_zone_start',          // integer bp; '' if unset
'left_boundary_zone_end',
'left_boundary_score',               // 0..1, 4dp; '' if unset
'left_boundary_support',             // pipe-separated track names
'left_boundary_support_class',       // strong|moderate|weak|ambiguous|''
'left_boundary_source',              // auto|manual|auto_plus_manual|''
'right_boundary_zone_start',
'right_boundary_zone_end',
'right_boundary_score',
'right_boundary_support',
'right_boundary_support_class',
'right_boundary_source',
'breakpoint_status',                 // boundary_zone_only|SV_supported|junction_supported|''
'boundary_notes',
```

That's **14 new columns**. Total registry width was 31 (post-v3.93);
becomes 45. Backward compatibility: missing values emit as `''` (empty
string, NOT `'NA'`) so tools that parse "header row defines shape" don't
break.

### `support_class` derivation

Confidence on each edge comes from how many evidence tracks back the
proposed boundary zone. Tracks whose contribution exceeded a per-track
threshold get listed in `support[]`. `support_class` is derived from
`|support|` (the count) AND the weighted `score`, because an edge with
one very strong track is qualitatively different from an edge with
three moderate tracks:

```
support_class :=
  'strong'      |support| >= 3   AND  score >= 0.60
  'moderate'    |support| == 2   OR  (|support| >= 3 AND 0.40 <= score < 0.60)
  'weak'        |support| == 1   OR  (|support| >= 2 AND score < 0.40)
  'ambiguous'   |support| == 0
```

### `TRACK_NAMES` enum (canonical set)

Values that may appear in `support[]`, with their **polarity** — `+1`
means "rises INSIDE the inversion" (so a rising step at the boundary
is a left-edge-like signal); `-1` means "falls inside" (the step gets
sign-flipped before the combined sum, so a falling step is correctly
treated as edge-like for ALL tracks):

```
'pca_drop'                +1   |Δ PVE1| or |Δ PC1| Z-score across the edge
'dosage_transition'       +1   |Δ mean dosage| between flanking windows
'band_continuity_drop'    -1   precomp band-continuity DROPS inside (regimes recluster)
'polarity_change'         +1   candidate_marker_polarity sign flip near edge
'similarity_edge'         -1   precomp similarity DROPS inside (haplotypes diverge)
'ghsl_step'               +1   |Δ ghsl_div_median| from ghsl_panel
'het_transition'          +1   |Δ het rate| from ghsl_panel primary scale
'fst_edge'                +1   FST between regimes is HIGH inside
'theta_pi_step'           -1   theta-pi within regime DROPS inside (suppressed recombination)
'discordant_pile'         +1   discordant-pair pileup spikes AT boundary
'sv_anchor'               +1   SV anchor density spikes AT boundary
'snp_density_quality'     —    QUALITY FILTER ONLY — never appears in support[]
```

The polarity convention solves a real-data issue: tracks like band
continuity and theta-pi *drop* inside the inversion because of
recombination suppression, while tracks like PCA and GHSL *rise*. The
algorithm sign-flips falling-inside tracks' steps so all four `support[]`
contributions add coherently into the combined-left and combined-right
summations.

`snp_density_quality` is a filter, not a contributor: when SNP density
falls below threshold within the candidate scan region, auto-propose
flags those windows as low-quality and excludes them from edge
detection. It is recorded as a quality flag on the boundary record but
NEVER appears in the `support[]` list.

### Scan region

Default scan region: `[start_bp - 1.5 Mb, end_bp + 1.5 Mb]`, clamped to
chromosome bounds. Configurable via toolbar in v3.97 (1 / 2 / 5 Mb).
The candidate's interior `[start_bp, end_bp]` IS scanned (so an inner
edge that's actually the true boundary can be detected when the
committed candidate was over-large), but flagged with an
"interior_edge" note when the proposed zone falls inside the
candidate's pre-existing interval.

When the candidate spans more than 3 Mb, the scan region expands to
`start_bp - 0.5 * span, end_bp + 0.5 * span` to keep the interior /
flanking ratio sensible.

### Auto-propose algorithm (scrubber-side)

For each candidate `c` with promoted `[start_bp, end_bp]`:

**Step 1 — Build per-track score arrays.** For each available track,
compute a per-window score over the scan region. Local-PCA window
grid is the natural unit (matches `state.data.windows`).

**Step 2 — Smooth scores with rolling median (width 3 windows).**
Suppresses single-window noise. Width 3 chosen because typical
local-PCA window step is 20 markers ≈ 1 kb; smoothing at 3 windows
≈ 3 kb ≪ scan_window_bp from the optional R-side layer.

Smooth the *score* (not the step) — smoothing the step array would
erase clean single-window transitions entirely (median of `[0, spike, 0]`
= 0), which is the opposite of what we want.

**Step 3 — Per-track signed step on the smoothed series:**

```
step[t][w] = smoothed[t][w+1] - smoothed[t][w]
```

**Step 4 — Per-track normalized magnitude (MAD scaling).** Divide each
step by the track's median absolute deviation over the scan region.
Tracks become comparable on a unitless 0..∞ scale where steps > 2 are
"interesting".

**Step 5 — Combined score per window per side:**

```
left_combined[w]  = Σ_t  weight[t] * max(0,  step[t][w])
right_combined[w] = Σ_t  weight[t] * max(0, -step[t][w])
```

Default `weight` (sums to 1.0):

```
pca_drop:               0.20
dosage_transition:      0.18
band_continuity_drop:   0.14
ghsl_step:              0.12
polarity_change:        0.10
het_transition:         0.06
similarity_edge:        0.05
fst_edge:               0.05    (only when boundary_evidence loaded)
theta_pi_step:          0.04    (only when boundary_evidence loaded)
discordant_pile:        0.04    (only when boundary_evidence loaded)
sv_anchor:              0.02    (only when boundary_evidence loaded)
```

When optional tracks are absent, their weights are removed and the
remaining weights are renormalized to sum to 1.0.

**Step 6 — Pick edge windows.** Strongest left-edge window in the left
half of the scan region; strongest right-edge window in the right half.
The two halves are searched independently. Each half's search excludes
the outermost 10% of the scan region (avoids edge artifacts at scan
boundaries).

**Step 7 — Build boundary zones.** Each picked window `w*` becomes a
zone `[windows[w*-radius].start_bp, windows[w*+radius].end_bp]` where
`radius = 5` windows by default. This matches the `state.cur ± 5`
cursor radius used elsewhere — consistent mental anchor.

**Step 8 — Populate `support[]`.** For each track `t` whose contribution
to `combined[w*]` exceeds 0.5 × the median nonzero contribution within
the scan, append `t` to `support[]`. Sort and dedupe.

**Step 9 — Compute `score` and `support_class`** per the formulas above.
Emit `boundary_left` and `boundary_right`. Set `breakpoint_status` to
`boundary_zone_only`. Auto-propose NEVER sets `SV_supported` or
`junction_supported` — the user must explicitly promote to those
statuses via the UI to acknowledge they have reviewed the underlying
evidence.

**Step 10 — SV anchor bonus.** When `boundary_evidence.tracks.sv_anchors`
contains an SV with `pos_bp` inside the auto-detected zone, copy it
into `boundary_left.sv_anchors_in_zone` (or right).

### Manual override (page-scoped hotkeys)

Hotkeys are page-scoped — they fire only when the boundaries page is
active and no text input has focus.

| Key | Action |
|---|---|
| **E** | Set `boundary_left` zone centered on `state.cur`, radius=5 windows. `source` becomes `manual` (or `auto_plus_manual` if right was auto-set). `support[]` is computed fresh against the user-chosen window using the same step-magnitude rules. |
| **F** | Symmetric — `boundary_right`. |
| **B** | Persist current boundary annotation (commit staging state into `state.candidates[i]` so it survives navigation). Auto-detected zones are *staged* but not *committed* until the user explicitly accepts. |
| **R** | Reset both `boundary_left` and `boundary_right`; `breakpoint_status` reverts to `''`. |

The page also provides explicit buttons for users who prefer mouse:

- "Auto-propose"           — runs the algorithm (re-runs if already done)
- "Accept auto boundaries" — stages → commits both edges (= B after auto)
- "Override left at cursor"  — equivalent to E
- "Override right at cursor" — equivalent to F
- "Reset boundaries"       — equivalent to R
- "Save"                   — equivalent to B
- `breakpoint_status` dropdown — `boundary_zone_only` /
  `SV_supported` / `junction_supported`. Defaults to
  `boundary_zone_only`. User manually promotes to stronger statuses
  when they have reviewed the evidence.

### Display

The boundaries page shows the candidate's scan region with stacked
evidence tracks on a shared bp axis:

| Track | Source | Always shown? |
|---|---|---|
| Local PCA (PVE1 + PC1 Z-score) | precomp | yes |
| Dosage heatmap (compact strip) | dosage_chunks | when loaded |
| Band continuity | precomp | yes |
| Marker polarity | candidate_marker_polarity | when loaded |
| GHSL change | ghsl_panel | when loaded |
| Heterozygosity | ghsl_panel | when loaded |
| FST | boundary_evidence | when loaded |
| Theta-pi (per regime) | boundary_evidence | when loaded |
| Discordant-pair pileup | boundary_evidence | when loaded |
| SV anchors | boundary_evidence | when loaded |
| Combined boundary score | computed scrubber-side | yes |

Proposed left/right boundary zones drawn as **shaded vertical bands**
across all tracks, color-coded by `support_class`:

| support_class | Hex | Opacity |
|---|---|---|
| strong    | `#1B7837` (green) | 0.20 |
| moderate  | `#F4A582` (amber) | 0.18 |
| weak      | `#FDDBC7` (peach) | 0.16 |
| ambiguous | `#BDBDBD` (grey)  | 0.10 |

SV anchors inside zones drawn as small triangles above the dosage
heatmap track. The candidate's pre-existing `[start_bp, end_bp]`
drawn as faint dashed vertical lines on each side — the *anchor* the
boundaries refine *from*.

The dosage heatmap on this page is the **same renderer as page 2** —
`drawDosageHeatmap` with the candidate scan region bound, capped at
the user-selected marker count. This means hover tooltips work
identically (v3.96 spec) and the visual idiom is the same across
pages.

### Failure modes (must handle gracefully)

1. **No layers loaded beyond precomp.** Auto-propose runs with the
   minimum tier (PCA drop only, possibly with GHSL/dosage if those
   load). `support_class` will be `weak` or `ambiguous`. UI shows a
   "limited evidence" notice with the list of layers that would
   improve the result.
2. **Scan region clipped at chromosome bounds.** Zones near chrom ends
   may be one-sided — fine, just displayed as such.
3. **Candidate larger than scan region (>3 Mb).** Scan expands to
   `start - 0.5*span, end + 0.5*span`. Toolbar can override.
4. **All track scores flat.** No edge detected; zones not proposed;
   UI shows "no clear boundary signal — manual annotation required".
5. **User presses E or F outside the candidate.** Allowed — the zone
   is created at that window. The candidate's pre-existing
   `[start_bp, end_bp]` is NOT updated (boundaries refine *toward*
   the true edges; the original interval stays as a reference).
6. **Saved boundaries when candidate already had boundaries.** New
   values overwrite. Previous values survive in candidate-history (a
   sibling array maintained for QC; spec out of scope here).
7. **`boundary_evidence` layer present but `candidate_id` not found.**
   Auto-propose falls back to the layer-absent path; no error.

### Cross-references

- §4 (`candidates_registry`) — the boundary fields extend the
  registry shape and add 14 TSV columns
- §10 (Marker module) — the dosage heatmap renderer is reused on
  this page; hover tooltips (§10 Hover affordances, schema 2.7)
  apply unchanged
- §11 (Band diagnostics) — `band_continuity_pct` from band
  diagnostics is the source of `band_continuity_drop` track when
  available; precomp's window-level `band_continuity_score` is the
  fallback
- DELLY/Manta integration: `STEP_M06_emit_boundary_evidence.R` (shipped
  post-v3.98) produces the `boundary_evidence` layer from existing SV
  pipeline outputs (MODULE_4x series; `inversion_codebase_v8.5`).
  Track-by-track availability: `--dosage` + `--fish_regimes` enables
  `fst` + per-regime θπ; `--sv_vcf` enables `sv_anchors`;
  `--discordant_bed` enables `discordant_pair_pileup`. Each track is
  independently optional.

### Open tuning parameters (v3.97 implementation)

These are settings, not contract:

- `track_weight` map — defaults given above; tunable after real-data
  validation
- `radius = 5` windows for zone half-width — matches cursor scrubber
  step; adjustable
- 0.5×median threshold for `support[]` inclusion — empirical default
- `support_class` thresholds — locked in schema, but documenting that
  v3.97 may surface tuning knobs for power users
