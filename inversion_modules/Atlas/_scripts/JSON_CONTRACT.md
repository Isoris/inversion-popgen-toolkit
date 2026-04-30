# Atlas — JSON Contract Reference

What's inside each precomp JSON, what it means, who produces it.

This is the single source of truth for the JSON layer schemas. Both
the R exporters (cluster-side) and the atlas (browser-side) write
against this contract. If they diverge, this document is the
arbiter.

---

## File overview

Three JSON files per chromosome:

| File | Atlas pages | Producer |
|---|---|---|
| `LG##_dosage.json` | page 1 (dosage local PCA) | existing dosage exporter (not built this session) |
| `LG##_ghsl.json` | page 3 (GHSL local PCA) + page 4 catalogue | `export_ghsl_to_json_v3.R` |
| `LG##_theta.json` | page 12 (θπ local PCA) | `STEP_TR_B_classify_theta.R` + `STEP_TR_D_augment_theta_json.R` |

All three are loaded into the atlas by drag-drop or Browse button.
The atlas merges their layers onto a single `state.data` object.
Cross-page features (e.g. coloring by GHSL kstripe in the page-1
PCA scatter) read from the merged layers regardless of which file
contributed them.

---

## Top-level structure (all files)

Every precomp JSON has these top-level keys:

```jsonc
{
  "schema_version": 2,                      // integer, currently 2
  "chrom": "C_gar_LG28",                    // chromosome label
  "n_samples": 226,                         // sample count
  "n_windows": 4302,                        // window count for this stream
  "samples": [...],                         // optional, per-sample metadata
  "_layers_present": ["tracks", ...],       // which layers are populated
  "_generated_at": "2026-04-29T14:00:00",   // ISO8601 generation timestamp
  "_generator": "export_ghsl_to_json_v3.R"  // script that wrote this file

  // Layer payloads at top level — see per-layer specs below
}
```

The `_layers_present` array is the source of truth for what the
atlas should expect to find. The atlas's `detectSchemaAndLayers()`
validates that each declared layer has its top-level key populated
with the expected shape.

---

## GHSL JSON (`LG##_ghsl.json`) — 8 layers

Produced by `export_ghsl_to_json_v3.R` consolidating four cluster-side
sources: STEP_C04 v6 RDS files (annot, per_sample, karyotypes),
STEP_C04c local PCA RDS, STEP_C04d D17 wrapper TSVs.

### Layer: `tracks`

Per-window aggregate values for the cross-page chromosome track strip.

```jsonc
{
  "ghsl_score":          { "values": [...], "pos_bp": [...], "min": ..., "max": ..., "mean": ... },
  "ghsl_div_median":     { ... },
  "ghsl_rank_stability": { ... }
}
```

Each track is `{ values, pos_bp, min, max, mean }`. The `values` and
`pos_bp` arrays have length `n_windows`. The atlas's track strip
displays these as colored stripes across the chromosome.

### Layer: `ghsl_panel`

Per-sample × per-window divergence matrices, multiple smoothing
scales.

```jsonc
{
  "schema_version": 2,
  "primary_scale": "s50",
  "n_samples": 226,
  "n_windows": 4302,
  "sample_names": ["CGA001", ...],          // length n_samples
  "start_bp":     [12345, 17345, ...],      // length n_windows — window grid coordinates
  "end_bp":       [17345, 22345, ...],
  "div_roll": {                             // per-scale matrices
    "s10":  [[/*n_windows*/], ...],         // n_samples arrays, each n_windows long
    "s50":  [[/*n_windows*/], ...],
    "s100": [[/*n_windows*/], ...]
  },
  "div_median": [/*n_windows*/]             // cohort-mean pre-aggregated at primary_scale
}
```

The `div_roll[scale][sample_idx][window_idx]` is the rolling-smoothed
divergence (phased-het / total). Scales `s10` / `s50` / `s100` are
windows of 10, 50, 100 windows for smoothing. Atlas page 3-bis renders
this as a samples × windows heatmap.

### Layer: `ghsl_kstripes`

K=2..6 stripe assignments per sample (cross-page coloring).

```jsonc
{
  "schema_version": 2,
  "primary_scale": "s50",
  "n_samples": 226,
  "sample_mean_rank": [/*n_samples, 0..1*/], // rank-mean-normalized
  "by_k": {
    "2": { "stripe_per_sample": [1,2,1,...], "n_per_stripe": [113, 113] },
    "3": { "stripe_per_sample": [1,2,3,...], "n_per_stripe": [76, 75, 75] },
    "4": { ... },
    "5": { ... },
    "6": { ... }
  }
}
```

`stripe_per_sample[si]` is an integer in `[1..K]` — which stripe
sample `si` belongs to at that K. Used by the atlas's
`_resolveSampleColorByMode("ghsl", ...)` for per-sample coloring on
page-1 lines, page-12 lines, all PCA scatters.

### Layer: `ghsl_karyotype_runs`

Stable LOW/HIGH per-sample contiguous regions from STEP_C04 v6.

```jsonc
{
  "schema_version": 2,
  "n_runs": 487,
  "runs": [
    {
      "sample_name": "CGA001",
      "start_bp": 6420000,
      "end_bp":   6655000,
      "class":    "LOW",                    // "LOW" | "HIGH" | "UNKNOWN"
      "n_windows": 47,
      "mean_score": 0.812
    },
    ...
  ]
}
```

Atlas page 3-bis renders these as horizontal sample-bar tracks.

### Layer: `ghsl_envelopes` (PRIMARY candidate set)

Production candidate envelopes from STEP_C04b PASS-runs (the
calibrated upstream classifier). This is the **PRIMARY** GHSL
candidate set per ADR-2.

```jsonc
{
  "schema_version": 2,
  "layer": "ghsl_envelopes",
  "chrom": "C_gar_LG28",
  "source": "STEP_C04b annot RDS (ghsl_v6_status == 'PASS')",
  "z_threshold_equivalent": null,           // n/a — not |Z|-derived
  "l2_envelopes": [
    {
      "l2_id":      "C_gar_LG28_ghsl_L2_001",
      "win_start":  1284,                   // 0-based window indices
      "win_end":    1331,
      "start_bp":   6420000,                // bp coordinates
      "end_bp":     6655000,
      "span_kb":    235.0,
      "n_windows":  48,
      "n_pass":     45,                     // PASS-status windows in run
      "n_weak":     2,                      // WEAK-status (allowed inside merge gap)
      "n_fail":     1,
      "mean_score": 0.873,
      "peak_score": 0.957
    }
  ],
  "l1_envelopes": [
    {
      "l1_id":      "C_gar_LG28_ghsl_L1_001",
      "win_start":  1280,
      "win_end":    1340,
      "start_bp":   6400000,
      "end_bp":     6700000,
      "span_kb":    300.0,
      "n_l2":       2                        // # of L2 envelopes inside this L1
    }
  ]
}
```

**L2 = tight individual candidates. L1 = merged groups (L2's within
9 windows of each other).** Both ship; the atlas overlays L2 as
gold rectangles on the page-3 sim_mat heatmap.

### Layer: `ghsl_local_pca`

Per-window local PCA outputs (the heart of the page-3 default view).

```jsonc
{
  "schema_version": 2,
  "layer": "ghsl_local_pca",
  "chrom": "C_gar_LG28",
  "scale": "5kb",
  "pad": 1,                                 // local-PCA neighborhood half-width
  "n_samples": 226,
  "n_windows": 4302,
  "sample_order": ["CGA001", ...],

  // Per-window PCA loadings — sign-AMBIGUOUS
  "pc1_loadings": [[/*n_samples*/], ...],   // n_windows arrays
  "pc2_loadings": [[/*n_samples*/], ...],

  // Per-window PCA loadings — sign-ALIGNED to anchor (what the lines panel renders)
  "pc1_loadings_aligned": [[/*n_samples*/], ...],
  "pc2_loadings_aligned": [[/*n_samples*/], ...],

  // Per-window eigenvalues
  "lambda_1":     [/*n_windows*/],
  "lambda_2":     [/*n_windows*/],
  "lambda_ratio": [/*n_windows*/],          // lambda_1 / lambda_2

  // Per-window robust |Z| profile (1D, sign-invariant)
  "z_profile":    [/*n_windows*/],
  "z_top10_mean": [/*n_windows*/],

  // Per-window MDS coordinates (cmdscale k=2 of (1-sim_mat))
  "mds_coords": {
    "mds1": [/*n_windows*/],
    "mds2": [/*n_windows*/]
  },

  // Anchor window for sign alignment (highest |Z| window)
  "anchor_window_idx": 1284,

  // Window×window similarity (sign-invariant via |cor|)
  "sim_mat_format": "upper_triangle_float32",
  "sim_mat_n":      4302,
  "sim_mat":        [/*n*(n+1)/2 floats, row-major upper triangle*/],

  "_compute_meta": {
    "weighting":       "n_phased_het",
    "smoothing_input": "div_roll_s50"
  }
}
```

**Sign alignment is critical.** Per-window local PCA produces sign-
ambiguous PC1/PC2 vectors. Without alignment, scrolling the cursor
across windows shows random sign flips that look like noise. The
algorithm: pick anchor (argmax z_profile), flip every other window's
PC1 to maximize correlation with anchor's PC1. PC2 aligned to
anchor's PC2 after PC1 is fixed.

**Sim_mat packing:** flat array of length `n * (n+1) / 2` storing
the upper triangle row-major. Atlas reader unpacks via
`_unpackGhSimMat()`, caches full n×n on `state._ghSimFull`.

### Layer: `ghsl_secondary_envelopes`

|Z|-threshold envelopes from STEP_C04c local PCA (secondary cross-
check; ADR-2).

```jsonc
{
  "schema_version": 2,
  "layer": "ghsl_secondary_envelopes",
  "chrom": "C_gar_LG28",
  "source": "STEP_C04c local-PCA z_profile threshold",
  "z_threshold":    2.5,
  "min_l2_windows": 5,
  "merge_gap":      3,
  "secondary_l2_envelopes": [
    {
      "sec_id":    "C_gar_LG28_ghsl_sec_001",
      "win_start": 1284,
      "win_end":   1331,
      "start_bp":  6420000,
      "end_bp":    6655000,
      "span_kb":   235.0,
      "n_windows": 48,
      "peak_z":    5.42,
      "mean_z":    3.71
    }
  ],
  "secondary_l1_envelopes": [...]            // similar shape, n_l2 instead of peak_z
}
```

### Layer: `ghsl_d17_envelopes`

D17 cross-block boundary detector outputs (secondary cross-check;
ADR-2).

```jsonc
{
  "schema_version": 2,
  "layer": "ghsl_d17_envelopes",
  "chrom": "C_gar_LG28",
  "source": "STEP_C04d D17 cross-block boundary detector",
  "detector": "STEP_D17_multipass_L1_only_v7.R + STEP_D17_multipass_L2_v8.R",

  "l1_envelopes": [
    {
      "candidate_id": "C_gar_LG28_d17L1_0001",
      "start_w":      1280,                // window indices
      "end_w":        1340,
      "start_bp":     6400000,
      "end_bp":       6700000,
      "n_windows":    61,
      "mean_sim":     0.732,
      "density_p70":  0.812,
      "status":       "STABLE_BLUE"        // STABLE_BLUE | EDGE | FAKE | DEMOTED
    }
  ],
  "l2_envelopes": [
    {
      "candidate_id": "C_gar_LG28_d17L2_0001_03",
      "start_w":      1284,
      "end_w":        1331,
      "start_bp":     6420000,
      "end_bp":       6655000,
      "n_windows":    48,
      "mean_sim":     0.867,
      "density_p70":  0.913,
      "status":       "STABLE_BLUE",
      "parent_l1_id": "C_gar_LG28_d17L1_0001"
    }
  ],

  "l1_boundaries": [
    {
      "boundary_idx":      "C_gar_LG28_d17L1_b001",
      "boundary_w":        1280,
      "boundary_bp":       6400000,
      "validation_status": "STABLE_BLUE",
      "boundary_score":    3.214,
      "grow_max_z":        -0.087           // grow validator output
    }
  ],
  "l2_boundaries": [...]                    // similar shape with parent_l1_id
}
```

**status semantics** (D17's L1 classification):
- `STABLE_BLUE`: cross-block stays bluish (low similarity) at all W
  in the grow series — REAL boundary
- `EDGE`: only the largest W validates — boundary at a chromosome
  edge or near a region transition
- `FAKE`: cross-block drifts toward positive z at large W —
  interior cross of a single inversion, not a true boundary
- `DEMOTED`: was a candidate but failed multiple validators

---

## θπ JSON (`LG##_theta.json`) — 5 layers

Produced by `STEP_TR_B_classify_theta.R` (4 layers) + augmented by
`STEP_TR_D_augment_theta_json.R` (adds 5th layer).

### Layer: `tracks`

Same shape as GHSL's tracks layer. Tracks contributed:
`theta_pi_median`, `theta_pi_z`, `lambda_ratio`.

### Layer: `theta_pi_per_window`

Per-sample × per-window θπ matrix.

```jsonc
{
  "schema_version": 2,
  "scale_label":    "win10000.step2000",
  "n_samples":      226,
  "n_windows":      16500,
  "samples":        ["CGA001", ...],
  "windows": [
    { "window_idx": 0, "start_bp": 0,    "end_bp": 9999,  "pos_bp": 5000 },
    { "window_idx": 1, "start_bp": 2000, "end_bp": 11999, "pos_bp": 7000 },
    ...
  ],
  "values":  [/*n_samples × n_windows, row-major*/],   // Float32Array-shaped
  "n_sites": [/*n_windows*/]                           // callable site count per window
}
```

### Layer: `theta_pi_local_pca`

Per-window local PCA on per-sample θπ. Same shape as GHSL's
`ghsl_local_pca` BUT v4 currently OMITS sim_mat, mds_coords, and
sign-aligned loadings. The v5 retrofit (deferred — ADR-7) will add
those for cross-stream parity.

```jsonc
{
  "schema_version": 2,
  "scale_label":    "win10000.step2000",
  "pad":            1,
  "n_samples":      226,
  "n_windows":      16500,
  "sample_order":   ["CGA001", ...],

  "pc1_loadings":   [[/*226*/], ...],       // sign-AMBIGUOUS, length n_windows
  "pc2_loadings":   [[/*226*/], ...],

  "lambda_1":       [/*n_windows*/],
  "lambda_2":       [/*n_windows*/],
  "lambda_ratio":   [/*n_windows*/],

  "z_profile":      [/*n_windows*/],
  "z_top10_mean":   [/*n_windows*/]

  // Pending v5 retrofit (currently absent):
  // "pc1_loadings_aligned": ...,
  // "pc2_loadings_aligned": ...,
  // "anchor_window_idx":    ...,
  // "mds_coords":           {...},
  // "sim_mat_format":       ...,
  // "sim_mat_n":            ...,
  // "sim_mat":              ...
}
```

### Layer: `theta_pi_envelopes` (PRIMARY |Z|-threshold detector)

|Z|-threshold envelopes from STEP_TR_B's 1D scan. PRIMARY for θπ
because no upstream calibrated detector exists (ADR-2).

```jsonc
{
  "schema_version": 2,
  "layer":          "theta_pi_envelopes",
  "chrom":          "C_gar_LG28",
  "source":         "STEP_TR_B |Z| threshold (PRIMARY)",
  "z_threshold":    2.5,
  "min_l2_windows": 5,
  "merge_gap":      3,
  "l2_envelopes": [
    {
      "l2_id":     "C_gar_LG28_theta_L2_001",
      "win_start": 4220,
      "win_end":   4287,
      "start_bp":  8420000,
      "end_bp":    8540000,
      "span_kb":   120.0,
      "n_windows": 68,
      "peak_z":    5.91,
      "mean_z":    3.42
    }
  ],
  "l1_envelopes": [...]
}
```

### Layer: `theta_d17_envelopes` (PRIMARY D17 cross-block detector)

Added by STEP_TR_D post-processor from STEP_TR_C wrapper TSVs.
**Also PRIMARY for θπ** — both detectors run as parallel candidate
sets without primary/secondary distinction (ADR-2).

Same shape as GHSL's `ghsl_d17_envelopes` layer.

---

## Dosage JSON (`LG##_dosage.json`)

Produced by the existing dosage exporter (not built this session).
Layers consumed by atlas page 1, page 2, page 4 catalogue.

This document does NOT specify the dosage JSON contract — that's
in the atlas's existing `SCHEMA_V2.md` (cluster-side at
`inversion_modules/scrubber/docs/SCHEMA_V2.md`). Key layers known
from atlas code:

- `windows` — per-window dosage local PCA outputs (top-level array)
- `l1_envelopes` / `l2_envelopes` — D17 candidates
- `tracks` — per-window aggregates
- `samples` — per-sample metadata + cluster labels
- `candidate_proposals` — Phase 3 cluster proposals
- ... ~30 more layers (see atlas `detectSchemaAndLayers()` switch
  statement at HTML line ~24400 for the full list)

---

## Cross-cutting fields

### `_layers_present`

Array of strings naming which top-level layer keys are populated.
Atlas's `detectSchemaAndLayers()` reads this to decide what UI to
enable. Layers declared here that aren't actually populated are
flagged as schema mismatches.

### Window indices vs bp coordinates

All layers use **0-based window indices** for `win_start` / `win_end`
/ `start_w` / `end_w` fields. Bp coordinates use **1-based, inclusive**
for `start_bp` / `end_bp` (matching SAM/BED conventions for human
readability, even though they're inconsistent with the half-open
0-based BED standard).

### Sample order

Within a single JSON, the sample order is consistent across layers:
the order in `ghsl_panel.sample_names` (or
`theta_pi_per_window.samples`) matches the order in `pc1_loadings`
inner arrays, `kstripes.by_k.X.stripe_per_sample`, etc. This is
enforced by the exporter (writes from one canonical
`sample_order` source).

Across JSONs (dosage / GHSL / θπ): sample order may differ if the
streams' R pipelines pulled samples from different config files.
**The atlas tolerates this**: cross-page coloring lookups go through
sample names (string identity), not indices. As long as
`sample_names[i]` is the same string in both files, cross-page
features work.

---

## Adding a new field

When adding a field to an existing layer:

1. Update the producer (R script) to emit it
2. Update the consumer (atlas JS) to read it
3. Update this document's per-layer spec
4. Bump `schema_version` if the change is breaking; otherwise leave it

When adding a new layer:

1. Update the producer's exporter to write the layer
2. Add a `case '<layer_name>':` to atlas
   `detectSchemaAndLayers()` switch
3. Add the layer to the merge-fix patch's case statements (so
   drag-drop enrichment merge picks it up)
4. Add an entry to this document's "Layer:" list above
5. If the layer needs cross-page rendering, add a renderer (page-3
   style) or extend `_resolveSampleColorByMode` for coloring

---

## Honest caveats

- **The dosage JSON contract is documented elsewhere**, not here.
  This document covers only what was built this session (GHSL +
  θπ extensions).
- **The θπ v5 retrofit is pending.** When STEP_TR_B v5 ships,
  `theta_pi_local_pca` will gain sim_mat, mds_coords, and
  sign-aligned loadings. The deferred-fields section above will
  become populated.
- **Field names follow R/snake_case convention** because they
  originate cluster-side. The atlas reads them as-is.
- **No test enforces this document.** It's prose. Drift between
  this doc and actual JSON output is possible. If you hit a field
  mismatch, the actual JSON wins; update this doc to match.

---

*Generated end of April 2026 session. Schemas pin the GHSL 8-layer
contract and θπ 5-layer contract.*
