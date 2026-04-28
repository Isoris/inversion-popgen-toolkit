# SCHEMA_V2 — JSON layers contract for the local PCA scrubber

**Version:** 2.0 (2026-04-27)
**Status:** Specification — scrubber implements detection in v3.20

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
