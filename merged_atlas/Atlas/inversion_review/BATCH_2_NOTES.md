# Batch 2 — extraction notes

`inversion_review/` pages: page11, page_sv_evidence, page4, page6, page7.

All shared regression tests pass (455/455). All per-page module tests pass
(24/24). `shared/` was not modified.

---

## Files extracted

| File | Lines | Source (legacy) | Notes |
|---|---|---|---|
| `inversion_review/page11.html` | 16 | 7906–7921 | DOM mount slot for "4 boundaries" |
| `inversion_review/page11.js` | ~245 | 30175–30314 (`renderBoundariesPage`) + 30318–30345 (hotkey wiring) | Verbatim entry-point + 5-key hotkey handler. ~25 `_bnd*` helpers marked TODO_MISSING. |
| `inversion_review/page_sv_evidence.html` | 3 | 9329–9331 | Mount slot only (`#sv_evidence_root`). |
| `inversion_review/page_sv_evidence.js` | ~135 | (none — external `js/atlas_sv_evidence.js`) | Lifecycle wrapper around `window.AtlasSVEvidence`. |
| `inversion_review/page4.html` | 73 | 7573–7645 | Two-pane DOM: cand-list pane (left) + karyotype/tier pane (right). |
| `inversion_review/page4.js` | ~250 | 62508–62526 (`karyoState`) + 62800–62836 (`renderCandidateKaryotype`) + 62840–62861 (`_refreshSubviewButtonStyles`) + 63411–63492 (`renderCandidateTier`) | Verbatim dispatcher, tier renderer, button styler, page-local state. Body + grid helpers marked TODO_MISSING. |
| `inversion_review/page6.html` | 12 | 7647–7658 | Popstats stack DOM (toolbar + ps-stack + no-chrom + gallery tray). |
| `inversion_review/page6.js` | ~115 | (none — external `js/atlas_page6_wiring.js` + siblings) | Lifecycle wrapper around `window.renderPopstatsPage`. |
| `inversion_review/page7.html` | 11 | 7660–7670 | Ancestry stack DOM. |
| `inversion_review/page7.js` | ~110 | (none — external sibling renderer) | Lifecycle wrapper around `window.renderAncestryPage`. |

| Test | Pass | Fail |
|---|---|---|
| `tests/test_review_page11.js` | 5 | 0 |
| `tests/test_review_page_sv_evidence.js` | 5 | 0 |
| `tests/test_review_page4.js` | 6 | 0 |
| `tests/test_review_page6.js` | 4 | 0 |
| `tests/test_review_page7.js` | 4 | 0 |
| **Total** | **24** | **0** |

Shared regression suite (must keep passing through the migration):

| Test | Pass | Fail |
|---|---|---|
| `tests/test_shared_contingency.js` | 45 | 0 |
| `tests/test_shared_het_rate.js` | 20 | 0 |
| `tests/test_shared_hungarian.js` | 40 | 0 |
| `tests/test_shared_kmeans.js` | 35 | 0 |
| `tests/test_shared_per_l2_cluster.js` | 61 | 0 |
| `tests/test_shared_state.js` | 78 | 0 |
| `tests/test_modular_smoke.js` | 58 | 0 |
| `tests/test_legacy_parity.js` | 118 | 0 |
| **Total** | **455** | **0** |

---

## TODO_MISSING markers — frequency-sorted

Every TODO_MISSING is referenced once at its declaration site, by design
(declaration only — no hand-rolled stubs to chase). The merge chat should
treat this as a flat catalogue of unresolved deps.

### Boundary refinement (page11) — 25 markers, all unique

Boundary algorithm core (likely → `shared/boundaries.js`):
- `_ensureBoundariesState` — legacy 17763
- `_boundaryScanRange` — legacy 17813
- `_findSVAnchorsInZone` — legacy 17895
- `_buildBoundaryTrackScores` — legacy 17942
- `_computeBoundaryEdges` — legacy 18127
- `_buildBoundaryRecord` — legacy 18256

Page-local UI helpers (probably stay in `page11.js`):
- `_bndFindCandidate` — legacy 18297
- `_bndFmtBp` — legacy 18308
- `_bndCloneRecord` — legacy 18317
- `_bndStageFromCandidate` — legacy 18337
- `_bndPopulateCandidateSelect` — legacy 18351
- `_bndUpdateRadiusButtons` — legacy 18383
- `_bndUpdateSaveButton` — legacy 18394
- `_bndUpdateStatusSelect` — legacy 18404
- `_bndSelectCandidate` — legacy 18415
- `_bndAutoPropose` — legacy 18424
- `_bndManualOverride` — legacy 18479
- `_bndOverrideLeft` — legacy 18622
- `_bndOverrideRight` — legacy 18623
- `_bndReset` — legacy 18626
- `_bndSave` — legacy 18640
- `_bndRefreshSummary` — legacy 18658
- `_bndDrawTracks` — legacy 18726
- `_bndRefreshUI` — legacy 18911

Side panels on the boundary page:
- `_renderBndFocalVsBg` — legacy 20507 (focal-vs-bg widget; also used by page16)
- `_wireRepeatDensityEscapeReset` — legacy 20583
- `_wireRepeatDensityClassCycle` — legacy 20856
- `_wireRepeatDensityAllTeToggle` — legacy 20895
- `_wireRepeatDensityArrowNav` — legacy 20941

Cross-page (will be needed by pages 4, 6, 7, 11 → strong shared candidate):
- `_renderCandidateNavInline` — uncatalogued line; produces the prev/next
  bar with idPrefix per page

### Karyotype / tier (page4) — 8 markers

Karyotype body helpers (large copy-port):
- `_renderCandidateKaryotypeBody` — legacy 63091 (~200 LOC)
- `buildKaryotypeRows` — legacy 62703
- `_isKaryoTwoTrack` — uncatalogued; sigma-area helper
- `sigmaProfileCandidate` — uncatalogued; sigma-profile helper
- `getKaryotypeLabel` — legacy 37022 (also used by other pages)
- `getKaryotypeLabelCaveat` — legacy 37040 (also used by other pages)
- `groupColor` — palette helper (almost certainly cross-page → shared)
- `refreshCandidateListUI` — legacy 62528 (every page with a cand-list pane)

Tier-specific:
- `_renderTierAxesGrid` — legacy ~63495 (immediately after
  `renderCandidateTier`, not yet line-confirmed)

### SV evidence (page_sv_evidence) — 3 markers

External-file deps (live in `js/atlas_sv_evidence.js`, not in legacy HTML):
- `AtlasSVEvidence.init`
- `AtlasSVEvidence.loadCandidate`
- `AtlasSVEvidence.destroy`

### Popstats (page6) — 4 markers

External-file deps:
- `renderPopstatsPage` — `js/atlas_page6_wiring.js` (only ever referenced
  via `typeof === 'function'` guards at legacy 59626 and 59729 — never
  defined inside the legacy HTML)
- `popgenLive` — `js/atlas_request_layer.js`
- `popgenPage6` — `js/atlas_page6_wiring.js`
- `popgenGallery` — `js/atlas_track_gallery.js`

### Ancestry (page7) — 2 markers

External-file deps:
- `renderAncestryPage` — sibling of `atlas_page6_wiring.js` (same pattern
  as `renderPopstatsPage`; only `typeof === 'function'` guards at legacy
  59629 and 59732)
- ancestry layer loader (drag-drop for `<chrom>_phase4_ancestry.json`) —
  likely lives in main bootstrap

---

## Cross-page promotion candidates (recommendations for merge chat)

Functions referenced by ≥2 pages in this batch that should probably go to `shared/`:

| Function | Used by | Where it lives now |
|---|---|---|
| `_renderCandidateNavInline` | page4, page6, page7, page11 | inlined in legacy |
| `getKaryotypeLabel` / `getKaryotypeLabelCaveat` | page4 (and pages 1/2 from batch 1) | legacy 37022 / 37040 |
| `groupColor` | every page that draws a band | inlined in legacy |
| `refreshCandidateListUI` | page4 explicitly; cand-list pane appears on multiple pages | legacy 62528 |

The boundary algorithm core (`_ensureBoundariesState` through
`_buildBoundaryRecord`, legacy 17763–18295) is page11-only at present, but
the underlying scoring logic is generic enough that batch 5's comparative
page might want to reuse `_findSVAnchorsInZone`. Worth keeping the door
open during merge.

---

## state.* references not in SLOT_REGISTRY

Catalogued from each page module's footer. Splitting by intent:

### Already in SLOT_REGISTRY (confirmed)

- `state.candidate` (cross_atlas) — sub-fields `.id`, `.chrom`, `.start_bp`,
  `.end_bp`, `.K`, `.locked_labels`, `.tracks`, `.ref_l2` are intra-candidate
  schema, not slot-level
- `state.candidateList` (cross_atlas) — read by `_bndPopulateCandidateSelect`
- `state.data` (transient) — read by every page; sub-layers
  `state.data.final_classification`, `state.data.classification`,
  `state.data.ancestry` are layer-attached, not slot-level

### Page-private (do NOT promote)

- `state.repeatDensity` — TE density layer, page11-only consumer
- `state.ncRNADensity` — ncRNA density layer, page11-only consumer
- `state.popstatsLive` — IndexedDB-cached responses, owned by
  `atlas_request_layer.js`, page6-only
- `state.popstatsTracksOn` — Set of active chip IDs, page6 UI flag
- `state.popstatsGalleryOpen` — gallery tray collapsed flag, page6 UI flag
- `state.ancestryViewChips` — Set of active view chips, page7 UI flag

### Promotion candidate

- `state.activeMode` — referenced inside `_renderCandidateKaryotypeBody`
  (legacy two-mode toggle, "simple" vs "detailed"). Currently NOT in
  SLOT_REGISTRY but it's a persisted user preference. Recommend adding to
  SLOT_REGISTRY as `persisted` with key `inversion_atlas.activeMode`.

---

## Decisions made (pre-empting "should this go in shared/?")

1. **`AtlasSVEvidence` stays external.** The SV evidence renderer lives in
   `js/atlas_sv_evidence.js` (script-tag inventory at legacy line 54981).
   Promoting it would mean porting IndexedDB caching + an UpSet plot lib;
   that's a separate extraction pass, not in scope for batch 2.
2. **`renderPopstatsPage` and `renderAncestryPage` stay external.** They
   are referenced via `typeof === 'function'` guards but have no inline
   body in `Inversion_atlas.html`. They are siblings of `atlas_page6_wiring.js`
   in the original deployment; the page module here is a lifecycle wrapper.
3. **`karyoState` stays page-local.** It's a UI-only object (sortKey,
   sortAsc, filter, bandFilter, subview) used by exactly one page. The
   `subview` field has its own localStorage key (`pca_scrubber_v3.candSubview`)
   and that persistence is handled inside `page4.js` via the new
   `setKaryoSubview()` export.
4. **Hotkey wiring (`_bndKeyHandler`/`_bndAttachHotkeys`/`_bndDetachHotkeys`)
   stays in `page11.js`.** The handler is page11-specific (it gates on
   `#page11.active` and only fires the bnd hotkeys). The merge chat does
   not need to promote it.
5. **`renderBoundariesPage` was extracted verbatim.** All 30+ helper
   references kept as bare global lookups; the merge chat patches them.
   Did NOT inline any helpers I might "guess" at — pure copy-and-mark.
6. **`renderCandidateTier` was extracted verbatim including the entire
   empty-state HTML and value-state HTML blocks.** Both blocks reference
   `_renderTierAxesGrid` which is the immediate next function in the
   legacy (line ~63495); marked TODO_MISSING so the merge chat picks it up.

---

## Hard rules — confirmed compliance

1. ✅ `shared/` not modified (no files added, none changed).
2. ✅ `tests/test_shared_*.js`, `tests/test_modular_smoke.js`,
       `tests/test_legacy_parity.js` not modified.
3. ✅ Every unresolved reference marked `TODO_MISSING(name) — legacy line N`.
4. ✅ No attempt to make pages run end-to-end.
