# SV Evidence Page — Handoff for Next Chat

**Date:** 2026-05-04
**Project:** MS_Inversions_North_african_catfish manuscript — SV evidence page in `Inversion_atlas.html`
**Author:** Quentin Andres
**Working directory at session end:** `/home/claude/Atlas/` (mirror of his GitHub Desktop tree)

---

## Status overview

The SV evidence page is **functionally complete through step 5** (out of 7 steps in the build plan). All passing 730 tests. Steps 6 (heatmap) and step 7 (Python emitters + per-candidate folder writer) are the remaining pieces.

Plus one small isolated fix to the cross-species page's dotplot panel (hover-only-on-the-square), shipped this session.

### Steps shipped this session and earlier

| step  | feature                                                     | tests   | status |
|-------|-------------------------------------------------------------|---------|--------|
| 1     | Skeleton + tab + empty-state mount                          | 288/288 | DONE   |
| 2     | Loader + main SV table (sort/filter/paginate/TSV export)    | 198/198 | DONE   |
| 3     | Locus track strip (zones / axis / SV calls / placeholders)  | 57/57   | DONE   |
| 3.5   | Cursor + hover + zoom + E/F/Enter/Esc + ↑↓ + view presets   | 48/48   | DONE   |
| 4     | Right-rail boundary summary + interactive legend + row anno | 60/60   | DONE   |
| 4.5   | Region-select mode (drag → region filter)                   | 29/29   | DONE   |
| 5     | UpSet panel (the redirect) + folder-walk drag-drop          | 50/50   | DONE   |
| 6     | Sample × SV heatmap (right-rail compact + main-area large)  | —       | TODO   |
| 7     | Python emitters: STEP_SV_GT_AGG + STEP_SV_EVID_COMB         | —       | TODO   |
|       |                                                             |         |        |
| dotplot fix | Cross-species page dotplot hover only on SVG square   | manual  | DONE   |

**Total: 730/730 assertions passing across 7 test files.**

---

## Architecture summary

### File structure (in Quentin's GitHub repo)

```
Atlas/
├── Inversion_atlas.html         ← main HTML file (62k lines), SURGICAL EDITS ONLY
├── js/
│   ├── atlas_sv_evidence.js     ← step 1-5 module (3738 lines, ~129 KB)
│   └── atlas_dotplot.js         ← cross-species dotplot (modified for hover fix)
├── json/
│   └── sv_genotype_counts/
│       └── cand_LG28_15Mb.json  ← test/dev fixture
├── tests/sv_evidence/
│   ├── fixture_sv_genotype_counts_v1.json
│   ├── fixture_sv_evidence_combinations_v1.json
│   ├── test_step1_skeleton.js          (288 PASS)
│   ├── test_step2_table.js             (198 PASS)
│   ├── test_step3_locus.js              (57 PASS)
│   ├── test_step3_5_cursor.js           (48 PASS)
│   ├── test_step4_rightrail.js          (60 PASS)
│   ├── test_step4_5_select.js           (29 PASS)
│   ├── test_step5_upset.js              (50 PASS)
│   ├── smoke_test_step1.html            ← browser smoke tests
│   ├── smoke_test_step2.html
│   ├── smoke_test_step3.html
│   ├── smoke_test_step3_5.html
│   ├── smoke_test_step4.html
│   ├── smoke_test_step4_5.html
│   └── smoke_test_step5.html
└── specs_todo/
    ├── SPEC_sv_evidence_page__upset_redirect.md       ← step-5 redirect (new schema)
    └── SPEC_sv_evidence_page__per_candidate_folder.md ← per-cand folder layout
```

### State shape (all in `_state`)

```js
{
  // Step 1-2 core:
  rootEl, layer, layerLoading, layerError, activeCandidateId,
  filters, sortColumnId, sortDirection, pageSize, pageIndex,
  // Step 3 locus crosstalk:
  highlightedSvId,
  // Step 3.5 cursor + view:
  cursorBp, markerLeftBp, markerRightBp, viewPreset, hotkeysAttached,
  _lastWindow, _customWindow,
  // Step 4 annotations + legend highlight:
  rowAnnotations, highlightPattern,
  // Step 4.5 region-select:
  selectMode ('zoom' | 'select'), selection ({startBp, endBp}|null),
  _selectDrag,
  // Step 5 UpSet + folder-walk:
  combinationsLayer, selectedSamples (Set|null), activeCombinationIndex
}
```

### Surgical edits in Inversion_atlas.html (4 total)

Only 4 edits since project start, all surgical and additive:

1. Tab button at line ~4876: `<button class="tab-button" data-tab="page_sv_evidence">5b SV evidence</button>`
2. Page DOM at line ~8498: `<div id="page_sv_evidence" class="page">…<div id="sv_evidence_root"></div>…</div>`
3. Script include at line ~61861: `<script src="js/atlas_sv_evidence.js"></script>`
4. Tab dispatcher at line ~57077: SV evidence branch calls `AtlasSVEvidence.loadCandidate()` and `attachHotkeys`/`detachHotkeys` on activate/deactivate.

No further HTML edits since step 1. All step 2-5 work is purely additive in the JS module.

### Module API surface (public)

```js
window.AtlasSVEvidence = {
  // Lifecycle
  init({ root }),
  loadCandidate(cid?),              // resolves candidate from atlas state
  refresh(),                         // re-render everything

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
  setCombinationsLayer(obj),         // accepts validated layer

  _state, _const, _internals          // for tests
}
```

---

## Step 5 detail (latest, the redirect)

### Schema

`sv_evidence_combinations_v1.json` — sibling of `sv_genotype_counts_v1.json` per candidate:

```json
{
  "format_version": "sv_evidence_combinations_v1",
  "candidate_id": "cand_LG28_15Mb",
  "n_samples_total": 226,
  "evidence_types": [
    { "id": "left_SA",      "label": "Left split-read",  "side": "left",  "kind": "SA",     "tier": 1 },
    { "id": "right_SA",     "label": "Right split-read", "side": "right", "kind": "SA",     "tier": 1 },
    { "id": "left_PE",      "label": "Left PE",          "side": "left",  "kind": "PE",     "tier": 2 },
    { "id": "right_PE",     "label": "Right PE",         "side": "right", "kind": "PE",     "tier": 2 },
    { "id": "Manta_INV_GT", "label": "Manta INV",        "side": null,    "kind": "caller", "tier": 3 },
    { "id": "DELLY_INV_GT", "label": "DELLY INV",        "side": null,    "kind": "caller", "tier": 3 },
    { "id": "MAPQ0_left",   "label": "Left MAPQ0",       "side": "left",  "kind": "mapq0",  "tier": 4 },
    { "id": "MAPQ0_right",  "label": "Right MAPQ0",      "side": "right", "kind": "mapq0",  "tier": 4 }
  ],
  "combinations": [
    { "members": ["left_SA","right_SA","left_PE","right_PE"], "intersection_size": 42, "samples": ["s_h2h2_0", ...] },
    { "members": ["MAPQ0_left","MAPQ0_right"],                "intersection_size": 21, "samples": ["s_h1h1_0", ...] },
    { "members": ["Manta_INV_GT","DELLY_INV_GT"],             "intersection_size": 17, "samples": [...] },
    { "members": ["left_SA"],                                 "intersection_size":  8, "samples": [...] }
  ],
  "per_evidence_totals": {
    "left_SA":      { "n_samples": 67 },
    "right_SA":     { "n_samples": 65 },
    ...
  }
}
```

### Folder-walk drag-drop

`_handlePageDrop` now uses `webkitGetAsEntry()` to detect folder drops. Walks the folder tree, collects all `*.json` files, dispatches each by `format_version`:

- `sv_genotype_counts_v1` → `_state.layer`
- `sv_evidence_combinations_v1` → `_state.combinationsLayer`

Backwards-compatible: single-file drop still works via the same dispatcher. Future formats just need a new branch in `_ingestJsonText`.

### UpSet click handshake (current behaviour)

1. User clicks a top bar in the UpSet panel.
2. `_state.selectedSamples = new Set(combination.samples)`.
3. `_state.activeCombinationIndex = i`.
4. UpSet panel re-renders: clicked bar gets accent stroke + amber fill.
5. Locus readout shows a "samples: 42 fish [left_SA + right_SA + left_PE + right_PE] [clear]" pill.
6. Esc, the UpSet pill, or the readout pill all clear it.
7. Second click on the same bar toggles off.

### What step 5 does NOT do yet (intentional, blocked by schema)

The current `sv_genotype_counts_v1` does **not** carry per-sample carrier lists per SV. So when `selectedSamples` is set, we **can't** yet:

- Recompute the table's H1/H1, H1/H2, H2/H2 genotype-count columns within the selection.
- Dim individual SV glyphs in the locus based on whether their carriers intersect the selection.

Both are spec-extensible. The cleanest route: add a `samples_with_call` field per `sv_calls[]` entry in `sv_genotype_counts_v1` (a list of sample_ids that carry the SV). Step 6 or step 7 (whichever ships first) is the natural place to spec + use that.

---

## What's next — ordered priorities

### Step 6 — Sample × SV heatmap (next obvious step)

**Spec:** `SPEC_sv_evidence_page.md §3.6` originally; `SPEC_sv_evidence_page__upset_redirect.md` may overlap. Two views:

- **Right-rail compact** (below UpSet): N×M heatmap, samples × SVs, cell colour = AA/AB/BB/miss. Intended as a thumbnail.
- **Main-area large**: full-size, click-cell-to-highlight, hover-tooltip with `(sample_id, sv_id, AA/AB/BB)`.

Both need per-sample × per-SV genotype data. **This is the schema extension noted above.** Suggested approach:

1. Spec a new sibling layer `sv_support_by_sample_v1.json` (the third file in Quentin's per-candidate folder). Format: a sparse matrix of (sample_id, sv_id) → AA/AB/BB/miss code. Producer = MODULE_5A2.
2. Atlas-side: load it via the same folder-walk path. Render the heatmap from it.
3. Once `sv_support_by_sample` is loaded, **also** unlock the "dim non-selected SVs" and "recount columns within selection" features for step 5's UpSet click. Two birds, one schema.

### Step 7 — Python emitters

Three writers, all in MODULE_5A2 / MODULE_4 producer code:

1. **`STEP_SV_GT_AGG_aggregate_genotype_counts.py`** — emits `sv_genotype_counts_v1.json`. Inputs: DELLY2 + Manta VCFs, the candidate boundaries, the locked karyotype labels. Computes Fisher OR + FDR_BH + pattern label per SV.
2. **`STEP_SV_EVID_COMB_emit_combinations.py`** — emits `sv_evidence_combinations_v1.json`. Inputs: BAM-evidence track JSON (S7 / `phase_8_comparative_breakpoint_fragility`), DELLY/Manta GT calls, MAPQ0 regions. Computes per-sample 0/1 evidence vector → group by combination → top N.
3. **`STEP_SV_SUPPORT_emit_support_by_sample.py`** — emits `sv_support_by_sample_v1.json` (the new sibling for step 6). Per-sample × per-SV genotype matrix.

All three writers should write into the same per-candidate folder structure:

```
data/LG28/candidates/INV_LG28_003/
  ├── sv_genotype_counts.json
  ├── sv_evidence_combinations.json
  └── sv_support_by_sample.json
```

### Three deferred features to revisit after step 6/7

These are noted in code comments but skipped because they need schema upgrades:

1. **Recount table genotype columns within `selectedSamples`** (currently shows full-cohort counts even when an UpSet bar is active).
2. **Dim non-selected SV glyphs in locus** when `selectedSamples` is set.
3. **Per-cell hover tooltip on heatmap** showing genotype + sample id.

All three become trivial once `sv_support_by_sample_v1` is loaded.

---

## Critical constraints to maintain

These are non-negotiables Quentin reiterated. Any drift = bug:

1. **Three separate catfish cohorts; never conflate:**
   - F₁ hybrid (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
   - **226-sample pure *C. gariepinus* hatchery cohort on LANTA — current inversion work.** K clusters here are hatchery broodline structure, NOT species admixture.
   - Pure *C. macrocephalus* wild cohort — future paper.

2. **Page is read-only with respect to candidate boundaries.** Boundary editing belongs to page 11 (`_bndAttachHotkeys`). The SV evidence page does NOT modify `bs.staging`, `state.candidate.boundary_*`, the karyotype groups, or SV calls. The E/F keys here drop **visual-only pins** as read-only annotation; they look like page-11 muscle-memory mnemonics but cannot edit.

3. **Empty/partial states never throw.** Always render whatever subset of data is available. Missing combinations layer → empty hint. Missing per-sample data → fall back to cohort-wide counts. Etc.

4. **All colours via CSS variables.** Works in light AND dark mode. No hard-coded `#fff` or `#000`.

5. **Surgical HTML edits only.** `Inversion_atlas.html` has had only 4 edits since step 1 — purely additive. All step 2-5 work is in the JS module. Keep this discipline.

---

## Pickup commands for next chat

```bash
# Verify environment
cd /home/claude/Atlas
node --check js/atlas_sv_evidence.js
node --check js/atlas_dotplot.js

# Run all tests (730 should pass)
for f in tests/sv_evidence/test_step*.js; do
  echo "=== $(basename $f) ==="
  node "$f" 2>&1 | tail -3
done

# Browser smoke test for latest step
# open tests/sv_evidence/smoke_test_step5.html
```

For step 6, start by:
1. Reading `SPEC_sv_evidence_page__upset_redirect.md` and `SPEC_sv_evidence_page__per_candidate_folder.md` in `specs_todo/`.
2. Designing `sv_support_by_sample_v1.json` schema (sparse matrix; per-row = sample, per-col = SV; values = AA/AB/BB/miss codes; or per-row = SV with carrier list — either, but pick one).
3. Adding loader path in `_ingestJsonText` for the new format_version.
4. Building the heatmap renderer (compact + large).

For step 7, the Python emitter is greenfield. Match the JSON schemas exactly (they're locked by the atlas-side validators).

---

## Files in `/mnt/user-data/outputs/` from this session and prior

Step bundles are in zip form, plus standalone copies of the latest JS/spec/README files. The most current is `sv_evidence_step5_drop.zip` — full Atlas tree with all step 1-5 work + dotplot fix + handoff.

Earlier per-step bundles preserved for inspection:
- `sv_evidence_step1_drop.zip` through `sv_evidence_step4_5_drop.zip`
