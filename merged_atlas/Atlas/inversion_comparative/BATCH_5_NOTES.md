# Batch 5 — `inversion_comparative` extraction notes

**Pages:** page16 (cross-species overview), page16b (multi-species cockpit), page5 (help)
**Source:** `legacy/Inversion_atlas.html` (75,617 lines)
**Tests added:** `tests/test_comparative_page16.js`, `tests/test_comparative_page16b.js`, `tests/test_comparative_page5.js` (45 assertions across 3 files, all green)

---

## ⚠️ Important — handoff doc / actual content discrepancy

`HANDOFF_BATCH_5.md` describes page5 as the **"Multi-species comparison page (Pairwise dotplots, synteny graph, breakpoint annotations)"** at lines 8159–9316. **This is wrong.**

In legacy/Inversion_atlas.html, lines 8159–9316 are the **static Quick-reference / Help page** (`<button data-page="page5" data-stage="help">` at legacy line 5142).
- It is purely declarative HTML (CSS-styled tables of tabs, hotkeys, schema definitions, pipeline diagrams).
- It has no JS render function in the legacy file.
- The atlas's tab bar labels it as "16 help" — i.e., position 16 in the navbar, named "help".

The **actual** multi-species classification cockpit is **page16b** (`<div id="page16b">` at legacy line 8033, renderer `_renderMultiSpeciesPage` at legacy line 27497). I extracted page16b correctly per the line ranges in the handoff (8033–8100); the multi-species JS lives in `page16b.js`.

Decision: I followed the handoff's line ranges literally (8159–9316 → `page5.html`) and made the JS module a stub with explanatory comments. The merge chat should treat page5 as a static help page.

---

## Extraction summary

| Page | HTML lines | JS lines | Legacy JS source |
|---|---:|---:|---|
| page16 | 84 | 2,556 | `[20971..21114] ∪ [23717..26025] ∪ [28367..28419]` |
| page16b | 68 | 2,403 | `[26023..28366]` (incl. the two `DOTPLOT_MASHMAP_*` constants at 26023–26024 anchored as a prelude block) |
| page5 | 1,158 | 34 (stub) | none — static page |

Function bodies are extracted **verbatim** from the legacy file. No code was modified, refactored, or restructured. Per the handoff: "Don't try to make pages work end-to-end."

---

## TODO_MISSING summary (frequency-sorted, both pages combined)

| Symbol | Count | Legacy line | Notes / suggested resolution |
|---|---:|---:|---|
| `_esc` | **108** | 13925 | HTML-escape helper. Used by every page in the atlas. **TODO_PROMOTE_TO_SHARED** — belongs in `shared/dom.js` or `shared/format.js`. |
| `_getRepeatDensity` | 3 | 14477 | Accessor for `state.repeatDensity[chrom]`. Almost certainly already extracted by batch 2 (boundaries page lives there); merge chat should import from there. |
| `setCur` | 1 | 51747 | Scrubber: jump-to-window. Used only by `_csBpJumpToWindow`. |
| `drawZ` | 1 | 31839 | Scrubber: \|Z\| panel draw. Guarded by `typeof === 'function'`. |
| `drawSim` | 1 | 31331 | Scrubber: similarity-matrix draw. Guarded. |
| `drawLinesPanel` | 1 | 34894 | Scrubber: lines panel draw. Guarded. |
| `drawWinSumStrip` | 1 | (not defined) | Optional hook. Both call sites use `typeof === 'function'` guard, so absence is non-fatal. Resolution: leave as-is, or merge chat removes the dead reference. |
| `window.popgenDotplot` | (lib) | external | Vendor library — popgen dot plot panel renderer. Used by page16's `_renderCrossSpeciesDotplot` and page16b's multi-species dotplot rendering. Loaded as a `<script>` tag in the assembled HTML. |
| `window.popgenFocalVsBg` | (lib) | external | Vendor library — focal-vs-background panel renderer. Used by page16's `_renderCrossSpeciesFocalVsBg`. Same loading mechanism. |

The five scrubber draw helpers (`setCur`, `drawZ`, `drawSim`, `drawLinesPanel`, `drawWinSumStrip`) are only invoked from `_csBpJumpToWindow` to repaint the |Z| panel after a cross-species jump, and are wrapped in `typeof === 'function'` guards. **Page16 will degrade gracefully if these are absent.** They will be provided by batch 1's scrubber module (`page1.js` / `inversion_discovery/`). The merge chat wires them via either a shared scrubber bridge or unchanged `window.*` globals.

---

## Window-mounted exports (no TODO needed — defined and exported by these modules)

These are functions defined in the extracted modules and explicitly mounted on `window.*` for cross-page consumers (notably the page1 scrubber's canvas hover wiring):

**From page16.js:**
`window._csBpHitTestXFromList`, `window._csBpHitTest2D`, `window._csBpJumpToWindow`, `window._csBpHoverEnter`, `window._csBpHoverLeave`, `window._wireCsBpHoverOnCanvas`, `window._CS_BP_HOVER_TOL_PX`

**From page16b.js:**
The complete `_ms*` API surface (~25 functions) plus `window._isCompTEFragilityJSON`, `window._isDxyPerInversionJSON`, `window._isKaryotypeLineageJSON`, `window._isPhyloTreeJSON`, `window._isSyntenyMultispeciesJSON`, `window._storeCompTEFragility`, `window._storeDxyPerInversion`, `window._storeKaryotypeLineage`, `window._MS_DEFAULT_SPECIES`

The merge chat will decide whether to keep the window-global mount pattern or convert to ES module exports. The legacy file uses both — the window mount is the cross-page interop boundary.

---

## State slots touched (for `shared/state.js` SLOT_REGISTRY audit)

**Owned by page16.js (cross-species page):**
- `state.crossSpecies` — primary data, written by `_storeCrossSpecies`
- `state._crossSpeciesUI` — UI scratch (active_id, filters, etc.)
- `state._csSyntenyCache`, `state._csSyntenyEdgesCache`, `state._csInversionContextCache` — derived caches
- `state._csOverlayIndex` — per-window breakpoint overlay index
- `state._csDotplotPanel`, `state._csDotplotPanelFp` — popgenDotplot panel handles
- `state._csHoverActive`, `state._csHoverRaf`, `state._csHeaderScrollWired`, `state._crossSpeciesKeysBound` — UI flags
- `state._focalVsBg` — focal-vs-bg panel state (shared with page11)

**Owned by page16b.js (multi-species cockpit):**
- `state.dotplotMashmap` (loader; consumed also by page16)
- `state.syntenyMultispecies`, `state.phyloTree`, `state.dxyPerInversion`, `state.compTEFragility`, `state.karyotypeLineage` — six new JSON layers
- `state._multiSpeciesUI` — UI scratch (active_species, show_busco)
- `state._msClassifications`, `state.classifications` — classification persistence

**Read but not owned (cross-page):**
- `state.repeatDensity` (boundaries page / page11)
- `state.candidateList` (catalogue / page10)
- `state.cur`, `state.data` (scrubber / page1)

---

## Decisions worth flagging

1. **Where dotplot_mashmap_v1 IO lives.** The IO loader functions for `state.dotplotMashmap` (legacy 26026–26124) are co-located in `page16b.js` because the multi-species page is the primary consumer. Page16 reads `state.dotplotMashmap` only as an optional overlay for `_renderCrossSpeciesDotplot`. **The merge chat may decide to promote the IO loader to `shared/state_io.js` if both pages end up depending on it equally.**

2. **`_csBuildPermResultHtml` placement.** This function lives at legacy line 28367, *after* the entire multi-species block. It belongs logically to the cross-species permutation test (called from `_renderCrossSpeciesSynteny`'s click handler at legacy 25794). I extracted it into `page16.js` as Block C. The position in the legacy file appears to be an accident of patch ordering.

3. **`page5.js` is a stub.** Because legacy page5 is purely declarative HTML with no render function, the JS module is a stub with `renderPage5(state)` (no-op) and a `PAGE5_META = { id, stage, label, num, static: true }` object. This keeps page-loader routing uniform.

---

## Test results

```
$ node tests/test_comparative_page5.js
pass: 11   fail: 0

$ node tests/test_comparative_page16.js
pass: 9    fail: 0

$ node tests/test_comparative_page16b.js
pass: 25   fail: 0

# All shared tests still green (no shared/ files modified):
$ for t in tests/test_shared_*.js tests/test_modular_smoke.js; do node "$t"; done
test_shared_contingency.js:    pass: 45   fail: 0
test_shared_het_rate.js:       pass: 20   fail: 0
test_shared_hungarian.js:      pass: 40   fail: 0
test_shared_kmeans.js:         pass: 35   fail: 0
test_shared_per_l2_cluster.js: pass: 61   fail: 0
test_shared_state.js:          pass: 78   fail: 0
test_modular_smoke.js:         pass: 58   fail: 0

$ LEGACY_ATLAS=$PWD/../legacy/Inversion_atlas.html node tests/test_legacy_parity.js
pass: 118  fail: 0
```

**Total: 500 assertions across 11 test files, 0 failures.**

---

## Files shipped

```
Atlas/inversion_comparative/
├── page16.html       (   84 lines)
├── page16.js         (2,556 lines)
├── page16b.html      (   68 lines)
├── page16b.js        (2,403 lines)
├── page5.html        (1,158 lines)
├── page5.js          (   34 lines)
└── BATCH_5_NOTES.md  (this file)

Atlas/tests/
├── test_comparative_page16.js   (62 lines)
├── test_comparative_page16b.js  (74 lines)
└── test_comparative_page5.js    (60 lines)
```

`shared/` — untouched ✓
Legacy regression tests — all passing ✓
