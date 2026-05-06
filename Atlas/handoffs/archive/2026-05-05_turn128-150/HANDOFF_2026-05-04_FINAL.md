# SV Evidence Page + Producers — Final Handoff

**Date:** 2026-05-04
**Project:** MS_Inversions_North_african_catfish manuscript
**Atlas component:** SV evidence page (`5b SV evidence` tab) in `Inversion_atlas.html`
**HPC component:** Producer scripts under `producers/sv_evidence/`
**Working tree at session end:** `/home/claude/Atlas/`

**🟡 Read alongside `STEP6_DESIGN_NOTE.md`** — that note answers three architectural questions Quentin asked at session end (Procrustes/PCA-on-side: nonsense; do dosage and karyotype have the same shape: no, they share only the sample axis; what schema for step 6: reuse the existing `fixture_sv_support_by_sample_v1.json` shape with row-as-string compact form). Step 6 should not start without reading that note first.

---

## What's complete (this session and prior)

### Atlas SV evidence page — 7-step build plan

| step | feature | tests | status |
|---|---|---|---|
| 1 | Skeleton + tab + empty-state mount | 288/288 | ✅ DONE |
| 2 | Loader + main SV table (sort/filter/paginate/TSV export) | 198/198 | ✅ DONE |
| 3 | Locus track strip (zones / axis / SV calls / placeholders) | 57/57 | ✅ DONE |
| 3.5 | Cursor + hover + zoom + E/F/Enter/Esc + ↑↓ + view presets | 48/48 | ✅ DONE |
| 4 | Right-rail boundary summary + interactive legend + row annotations | 60/60 | ✅ DONE |
| 4.5 | Region-select mode (drag → region filter) | 29/29 | ✅ DONE |
| 5 | UpSet panel (the redirect) + folder-walk drag-drop | 50/50 | ✅ DONE |
| **6** | **Sample × SV heatmap** | — | **TODO** |
| 7 | Python emitters: STEP_SV_GT_AGG + STEP_SV_EVID_COMB + manifest + SLURM | 36/36 | ✅ DONE |

### Cross-species page — small fix

Dotplot hover now triggers only when the cursor is over the actual SVG square, not the unfilled left/right margins of the flex-centered panel. Bound to inner `<svg>`, re-bound on each `renderMini()`. Doc-level close-timer also tracks the SVG.

### Test totals

**766 / 766 assertions passing.**

```
tests/sv_evidence/test_step1_skeleton.js     288 PASS
tests/sv_evidence/test_step2_table.js        198 PASS
tests/sv_evidence/test_step3_locus.js         57 PASS
tests/sv_evidence/test_step3_5_cursor.js      48 PASS
tests/sv_evidence/test_step4_rightrail.js     60 PASS
tests/sv_evidence/test_step4_5_select.js      29 PASS
tests/sv_evidence/test_step5_upset.js         50 PASS
tests/producers/test_step7_producers.py       36 PASS  (Python end-to-end)
```

---

## File tree (latest)

```
Atlas/
├── Inversion_atlas.html             ← main HTML file, only 4 surgical edits since project start
├── HANDOFF_2026-05-04_FINAL.md      ← THIS FILE
├── js/
│   ├── atlas_sv_evidence.js          ← SV evidence module (3738 lines, ~129 KB)
│   └── atlas_dotplot.js              ← cross-species dotplot, with hover fix
├── json/sv_genotype_counts/
│   └── cand_LG28_15Mb.json           ← dev/test fixture
├── producers/sv_evidence/             ← Python producers (NEW this turn)
│   ├── README.md                     ← producer documentation
│   ├── write_candidate_folder.py     ← shared library (atomic write, manifest, layout)
│   ├── STEP_SV_GT_AGG_aggregate_genotype_counts.py
│   ├── STEP_SV_EVID_COMB_emit_combinations.py
│   └── run_sv_evidence_pipeline.slurm
├── tests/
│   ├── sv_evidence/
│   │   ├── fixture_sv_genotype_counts_v1.json
│   │   ├── fixture_sv_evidence_combinations_v1.json
│   │   ├── test_step1_skeleton.js   ... test_step5_upset.js   (7 test files)
│   │   └── smoke_test_step1.html    ... smoke_test_step5.html  (7 browser smoke tests)
│   └── producers/
│       └── test_step7_producers.py  ← end-to-end producer test
└── specs_todo/
    ├── SPEC_sv_evidence_page__upset_redirect.md       ← step-5 schema redirect
    └── SPEC_sv_evidence_page__per_candidate_folder.md ← per-candidate folder layout
```

---

## Three deliberate-non-features (blocked by schema; unlock by step 6)

The SV evidence page deliberately does *not* yet do these three things, because they all need a third schema layer (`sv_support_by_sample_v1`) that doesn't exist yet:

1. **Recount table genotype columns within `selectedSamples`.** When the user clicks an UpSet bar, the table still shows full-cohort genotype counts. Recomputing within the selection needs per-sample × per-SV genotype data.
2. **Dim non-selected SV glyphs in locus.** Same blocker.
3. **Per-cell hover tooltip on heatmap.** Step 6 needs to render the heatmap.

All three become trivial once `sv_support_by_sample_v1` is loaded. Step 6 is the natural place to spec + use it.

---

## Schemas (reference)

### `sv_genotype_counts_v1` — the main per-candidate layer

Produced by `STEP_SV_GT_AGG`. One file per candidate. Atlas displays steps 1-4.

```jsonc
{
  "format_version":     "sv_genotype_counts_v1",
  "candidate_id":       "INV_LG28_003",
  "chrom":              "C_gar_LG28",
  "boundary_left_bp":   15142000,
  "boundary_right_bp":  18124000,
  "zone_definitions_bp": {
    "left_flank":     [14140000, 14640000],
    "left_boundary":  [14640000, 15640000],
    "inversion_body": [15640000, 17620000],
    "right_boundary": [17620000, 18620000],
    "right_flank":    [18620000, 19120000]
  },
  "groups_used": {
    "H1/H1": {"n": 61,  "members": []},
    "H1/H2": {"n": 103, "members": []},
    "H2/H2": {"n": 62,  "members": []}
  },
  "sv_calls": [
    {
      "sv_id":               "SV001",
      "sv_type":             "BND" | "INV" | "DEL" | "DUP" | "Other",
      "chrom":               "C_gar_LG28",
      "position_bp":         15142380,
      "end_bp":              null,
      "zone":                "left_boundary",
      "distance_to_edge_bp": 380,
      "n_samples_with_call": 226,
      "quality":             "PASS",
      "callers":             ["delly2", "manta"],
      "genotype_counts": {
        "H1/H1": {"AA":59,"AB":2, "BB":0, "miss":0},
        "H1/H2": {"AA":3, "AB":98,"BB":2, "miss":0},
        "H2/H2": {"AA":0, "AB":3, "BB":58,"miss":1}
      },
      "fisher": {
        "comparison":  "H1/H1_vs_H2/H2",
        "odds_ratio":  42.1,
        "p_value":     1.2e-12,
        "fdr_bh":      2.1e-9
      },
      "pattern_label":      "canonical_breakpoint_marker",
      "notes":              ""
    }
  ],
  "boundary_summary": {
    "left":  {"interval_bp": [14640000, 15640000], "by_sv_type": { /* ... */ }},
    "right": {"interval_bp": [17620000, 18620000], "by_sv_type": { /* ... */ }}
  },
  "upset_top_combinations": []
}
```

### `sv_evidence_combinations_v1` — UpSet panel layer

Produced by `STEP_SV_EVID_COMB`. One file per candidate, sibling of the gt counts file.

```jsonc
{
  "format_version":   "sv_evidence_combinations_v1",
  "candidate_id":     "INV_LG28_003",
  "n_samples_total":  226,
  "evidence_types": [
    {"id":"left_SA",      "label":"Left split-read",  "side":"left",  "kind":"SA",     "tier":1},
    {"id":"right_SA",     "label":"Right split-read", "side":"right", "kind":"SA",     "tier":1},
    {"id":"left_PE",      "label":"Left PE",          "side":"left",  "kind":"PE",     "tier":2},
    {"id":"right_PE",     "label":"Right PE",         "side":"right", "kind":"PE",     "tier":2},
    {"id":"Manta_INV_GT", "label":"Manta INV",        "side":null,    "kind":"caller", "tier":3},
    {"id":"DELLY_INV_GT", "label":"DELLY INV",        "side":null,    "kind":"caller", "tier":3},
    {"id":"MAPQ0_left",   "label":"Left MAPQ0",       "side":"left",  "kind":"mapq0",  "tier":4},
    {"id":"MAPQ0_right",  "label":"Right MAPQ0",      "side":"right", "kind":"mapq0",  "tier":4}
  ],
  "combinations": [
    {"members":["left_SA","right_SA","left_PE","right_PE"], "intersection_size":42, "samples":["s_h2h2_0", ...]},
    {"members":["MAPQ0_left","MAPQ0_right"],                "intersection_size":21, "samples":[...]},
    {"members":["Manta_INV_GT","DELLY_INV_GT"],             "intersection_size":17, "samples":[...]},
    {"members":["left_SA"],                                 "intersection_size": 8, "samples":[...]}
  ],
  "per_evidence_totals": {
    "left_SA":      {"n_samples": 67},
    "right_SA":     {"n_samples": 65},
    /* ... */
  }
}
```

### `sv_support_by_sample_v1` — TODO for step 6

Not yet specified or built. Suggested shape:

```jsonc
{
  "format_version": "sv_support_by_sample_v1",
  "candidate_id":   "INV_LG28_003",
  "samples":        ["FL01_001", "FL01_002", ...],   // canonical row order
  "sv_ids":         ["SV001", "SV002", ...],          // canonical column order
  "encoding":       "AA=0, AB=1, BB=2, miss=-1",
  "matrix":         [[0,1,0,...],[1,0,2,...],...]   // n_samples × n_svs, dense
}
```

Sparse encoding could be considered if the matrix is mostly AA. Producer = `STEP_SV_SUPPORT_emit_support_by_sample.py` (greenfield).

---

## Step 6 plan (next chat's main task)

**Architectural decisions are in `STEP6_DESIGN_NOTE.md`. Read that first.** The short version:

- **Schema**: `sv_support_by_sample_v1` per-candidate JSON, sibling of the other two layers. Use the existing `fixture_sv_support_by_sample_v1.json` shape (row-as-string compact form + `row_groups` karyotype banding index). Already registered in `producers/sv_evidence/write_candidate_folder.py` `LAYER_FILENAME`.
- **No Procrustes panel.** No PCA-on-the-side. The heatmap with karyotype-banded sample ordering is the figure.
- **Karyotype shape ≠ SV-support shape.** They share only the sample axis. Karyotype is the row-annotation strip on the left edge of the heatmap, not a coordinate to align against.
- **Three deferred features unlock automatically** once `_state.supportLayer` is loaded: recount table on UpSet selection, dim non-selected SV glyphs, hover tooltip.

### Implementation order

1. `_validateSupportLayer` (mirrors `_validateCombinationsLayer`).
2. Add third branch in `_ingestJsonText` for `sv_support_by_sample_v1` → `_state.supportLayer`.
3. **Compact heatmap renderer** (right-rail, below UpSet): aggregated as group × SV mean dosage → 3 rows × n_svs cells.
4. **Large heatmap renderer** (main-area, replaces a card under the locus): n_samples × n_svs, hover tooltip, click-to-highlight, default sort = karyotype-banded by `row_groups`.
5. Wire UpSet `selectedSamples` → heatmap dim of non-selected rows.
6. Wire heatmap cell click → `_state.highlightedSvId` (already wired into locus glyph + table).
7. Unlock the three deferred features (now feasible).
8. Producer `STEP_SV_SUPPORT_emit_support_by_sample.py`. Reuse the GT parser from `STEP_SV_GT_AGG`.
9. Extend `run_sv_evidence_pipeline.slurm` to run all three producers.
10. Extend `tests/producers/test_step7_producers.py` with support-layer round-trip.

### Acceptance criteria

- Compact heatmap renders without scroll bars (fits 280px width); rows = three karyotype groups; cols = n_svs in canonical order.
- Main heatmap is interactive: hover shows `(sample_id, sv_id, GT)`; click cell sets `highlightedSvId` (locus glyph + table row both light up).
- UpSet bar click → heatmap dims non-selected sample rows.
- Esc clears all highlights including heatmap.
- Light + dark mode both work (CSS variables only).

---

## State shape (full reference)

```js
_state = {
  // Lifecycle
  rootEl, layer, layerLoading, layerError, activeCandidateId,

  // Filters / table
  filters, sortColumnId, sortDirection, pageSize, pageIndex,

  // Step 3 — locus crosstalk
  highlightedSvId,

  // Step 3.5 — cursor + view
  cursorBp, markerLeftBp, markerRightBp, viewPreset, hotkeysAttached,
  _lastWindow, _customWindow,

  // Step 4 — annotations + legend highlight
  rowAnnotations, highlightPattern,

  // Step 4.5 — region-select
  selectMode,        // 'zoom' | 'select'
  selection,         // {startBp, endBp} | null
  _selectDrag,       // {startBp, currentBp} | null

  // Step 5 — UpSet + folder-walk
  combinationsLayer,         // sv_evidence_combinations_v1 | null
  selectedSamples,           // Set<sample_id> | null
  activeCombinationIndex,    // number | null

  // Step 6 — heatmap (TBD)
  // supportLayer,           // sv_support_by_sample_v1 | null
  // heatmapSelection,       // {sampleId, svId} | null
}
```

---

## Public API surface

```js
window.AtlasSVEvidence = {
  // Lifecycle
  init({ root }),
  loadCandidate(cid?),
  refresh(),

  // Filters
  setFilters(filters),
  exportFilteredTSV(),

  // Step 3.5 — cursor + hotkeys
  attachHotkeys(),
  detachHotkeys(),
  setViewPreset('default'|'left_close'|'right_close'),
  setCursorBp(bp),

  // Step 4.5 — region select
  setSelectMode('zoom'|'select'),
  setSelection(startBp, endBp),

  // Step 5 — UpSet
  onUpSetBarClick(i),
  clearSampleSelection(),
  setCombinationsLayer(obj),

  _state, _const, _internals          // for tests
}
```

---

## Critical constraints (Quentin's non-negotiables)

1. **Three separate catfish cohorts; never conflate:**
   - F₁ hybrid (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
   - **226-sample pure *C. gariepinus* hatchery cohort on LANTA — current inversion work.** K clusters here are hatchery broodline structure, NOT species admixture.
   - Pure *C. macrocephalus* wild cohort — future paper.

2. **The SV evidence page is read-only with respect to candidate boundaries.** Boundary editing belongs to page 11. The page does NOT modify `bs.staging`, `state.candidate.boundary_*`, the karyotype groups, or SV calls. The E/F keys here drop *visual-only pins* as read-only annotation.

3. **Empty/partial states never throw.** Always render whatever subset of data is available. Missing combinations layer → empty hint. Missing per-sample data → fall back to cohort-wide counts.

4. **All colours via CSS variables.** Works in light AND dark mode. No hard-coded `#fff` or `#000`.

5. **Surgical HTML edits only.** `Inversion_atlas.html` has had only 4 edits since step 1, all additive. Maintain this discipline.

6. **Producers write atomically + update the manifest.** Don't shortcut to direct `open()` writes.

7. **The pattern classifier's het-specific rule fires before the FDR gate.** This is intentional — see `STEP_SV_GT_AGG_aggregate_genotype_counts.py` `_classify_pattern` comment for why.

---

## Pickup commands for next chat

```bash
# Verify working tree
cd /home/claude/Atlas
node --check js/atlas_sv_evidence.js
node --check js/atlas_dotplot.js

# Run all atlas tests (730 should pass)
for f in tests/sv_evidence/test_step*.js; do
  echo "=== $(basename $f) ==="
  node "$f" 2>&1 | tail -3
done

# Run producer tests (36 should pass)
python3 tests/producers/test_step7_producers.py

# Browser smoke test for latest step
open tests/sv_evidence/smoke_test_step5.html
```

For step 6 start by:
1. Reading the 3 deferred-features list in this handoff (the heatmap unlocks all three).
2. Spec the `sv_support_by_sample_v1` schema (decide sparse vs dense first — talk to Quentin).
3. Add a loader path in `_ingestJsonText` for the new format_version.
4. Build the heatmap renderer (compact + large).
5. Wire the deferred-feature unlocks: when `supportLayer` is present, recount table on UpSet selection; dim non-selected glyphs; per-cell tooltip.
6. Add `STEP_SV_SUPPORT_emit_support_by_sample.py` producer.
7. Update the SLURM wrapper to optionally run all three producers.
8. End-to-end test in `tests/producers/test_step6_step7_combined.py` or similar.

Stop point in current code: search for `_state.combinationsLayer` and `_state.selectedSamples` to find every place that the third layer would integrate. The TODOs are clearly marked.
