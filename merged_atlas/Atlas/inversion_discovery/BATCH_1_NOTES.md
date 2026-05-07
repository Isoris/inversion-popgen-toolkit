# Batch 1 notes

**Author**: parallel-chat-1 (inversion_discovery extraction).
**Scope**: pages 1, 2, 8, 12, 15, 19 — extracted LITERALLY from `legacy/Inversion_atlas.html`.
**Status**: 6 modules + 6 HTML fragments + 6 minimal tests shipped. All 455 baseline
shared/ tests still pass; the 6 new module-load tests add 45 checks (all green).

## What was extracted

| Page | Tab name | Stage | HTML lines | Functions extracted | Notes |
|---|---|---|---|---|---|
| page1  | "1 local PCA \|z\|"  | discovery  | 5474–6817 | 26 | The big one. All canonical renderers (drawSim/drawZ/drawLinesPanel/drawPCA/drawAnchorStrip/renderL3Panel + state mutators) extracted. |
| page12 | "2 local PCA θπ"     | discovery  | 6831–7173 | 8  | All `_drawTh*` panel renderers + θπ status helpers. |
| page15 | "2b local PCA GHSL"  | discovery  | 7180–7246 | 1  | Only `_refreshGhslLayerStatus` exists in legacy; the panel renderers are not yet authored. |
| page2  | "3 candidate focus"  | discovery  | 7248–7259 | 2  | `renderCandidateMetadata` (composes ~15 sub-panels via helpers) + `wireCandidateNav`. |
| page8  | "10 windows"         | refinement | 7672–7774 | 0  | **No JS handlers in legacy yet** — pure HTML scaffold. winSumTable / winSumStripCanvas / winSumZFilter handlers will be authored fresh. |
| page19 | "6 negative regions" | discovery  | 7378–7571 | 0  | **No JS handlers in legacy yet** — pure HTML scaffold. nrLoadBtn / nrExportCsvBtn / nrTableSlot handlers will be authored fresh. |

## TODO_MISSING summary

### Functions referenced but not extracted

Sorted by frequency (most-referenced first). Frequency = number of pages each
function is called from. Suggested location is a heuristic — merge chat decides.

| Function | Frequency | Pages referring | Suggested location |
|---|---|---|---|
| `toX` | 2 | page1, page12 | ?? — let merge chat decide |
| `toY` | 2 | page1, page12 | ?? — let merge chat decide |
| `_assignCandidateLanes` | 1 | page1 | shared/candidate_helpers.js |
| `_buildJumpMask` | 1 | page1 | ?? — let merge chat decide |
| `_drawBandTraceStrip` | 1 | page1 | shared/band_trace.js |
| `_drawDiamondOverlay` | 1 | page1 | ?? — let merge chat decide |
| `_drawInheritanceLabelsStrip` | 1 | page1 | shared/band_trace.js |
| `_drawLineageStrip` | 1 | page1 | shared/lines_panel.js |
| `_drawRegimeBreadthStrip` | 1 | page1 | shared/band_trace.js |
| `_drawSnpDensityShade` | 1 | page1 | shared/density_overlays.js |
| `_drawSnpDensityStrip` | 1 | page1 | shared/density_overlays.js |
| `_drawTrackedLinkageStrip` | 1 | page1 | shared/tracked_panel.js |
| `_drawTransitionRateStrip` | 1 | page1 | shared/density_overlays.js |
| `_drawWRow` | 1 | page1 | ?? — let merge chat decide |
| `_drawWinNavLane` | 1 | page1 | ?? — let merge chat decide |
| `_ensureCsOverlayIndex` | 1 | page1 | ?? — let merge chat decide |
| `_navigateToCandidate` | 1 | page2 | shared/candidate_helpers.js |
| `_paintCandidateBands` | 1 | page1 | shared/candidate_helpers.js |
| `_refreshScreeInset` | 1 | page1 | ?? — let merge chat decide |
| `_resolveSampleScopeColor` | 1 | page1 | ?? — let merge chat decide |
| `_vColor` | 1 | page1 | ?? — let merge chat decide |
| `_wRowBand` | 1 | page1 | ?? — let merge chat decide |
| `_winNavBand` | 1 | page1 | ?? — let merge chat decide |
| `_wireCandidateBandClicks` | 1 | page2 | shared/candidate_helpers.js |
| `_wireCandidateBlockChips` | 1 | page2 | shared/candidate_helpers.js |
| `_wireCandidateDosageHeatmap` | 1 | page2 | shared/candidate_helpers.js |
| `_wireCandidateHaplotypeAnnotations` | 1 | page2 | shared/candidate_helpers.js |
| `_wireCandidateRegimeRow` | 1 | page2 | shared/candidate_helpers.js |
| `addCandidateToList` | 1 | page2 | shared/candidate_helpers.js |
| `allSampleIdx` | 1 | page1 | ?? — let merge chat decide |
| `candidateAgeOriginHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateAncestryConfoundHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateBandComposition` | 1 | page2 | shared/candidate_helpers.js |
| `candidateBandsHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateBlockChipsHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateDosageHeatmapHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateFromJSON` | 1 | page2 | shared/candidate_helpers.js |
| `candidateHaplotypeAnnotationsHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateHeaderHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateHetShapeHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateListClosestIndex` | 1 | page2 | shared/candidate_helpers.js |
| `candidateListIndexOf` | 1 | page2 | shared/candidate_helpers.js |
| `candidateListSortedByPos` | 1 | page2 | shared/candidate_helpers.js |
| `candidateNavHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateNotesHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateProfileHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateRegimeRowHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateRichCardHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateSigmaChartHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateSubbandHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateSummaryHtml` | 1 | page2 | shared/candidate_helpers.js |
| `candidateToJSON` | 1 | page2 | shared/candidate_helpers.js |
| `colorFor` | 1 | page12 | ?? — let merge chat decide |
| `currentMbRange` | 1 | page1 | ?? — let merge chat decide |
| `dataset` | 1 | page1 | ?? — let merge chat decide |
| `drawCandGHSLPerBand` | 1 | page2 | shared/ghsl.js |
| `drawCandLinesPanel` | 1 | page2 | shared/lines_panel.js |
| `drawCandLocalPCA` | 1 | page2 | shared/pca_panel.js |
| `drawCandidateBar` | 1 | page1 | shared/candidate_helpers.js |
| `drawCandidateLocationStrip` | 1 | page2 | shared/candidate_helpers.js |
| `drawCandidateSigmaChart` | 1 | page2 | shared/candidate_helpers.js |
| `drawRect` | 1 | page1 | ?? — let merge chat decide |
| `escapeHtml` | 1 | page1 | ?? — let merge chat decide |
| `families` | 1 | page1 | ?? — let merge chat decide |
| `fillFor` | 1 | page12 | ?? — let merge chat decide |
| `fitCanvas` | 1 | page1 | shared/draw_utils.js |
| `flushRun` | 1 | page1 | ?? — let merge chat decide |
| `formatTrackVal` | 1 | page1 | ?? — let merge chat decide |
| `getActiveSimScale` | 1 | page1 | ?? — let merge chat decide |
| `getL2Cluster` | 1 | page1 | ?? — let merge chat decide |
| `getLinesGrid` | 1 | page1 | shared/lines_panel.js |
| `getLinesSignAt` | 1 | page1 | shared/lines_panel.js |
| `getLinesValuesAt` | 1 | page1 | shared/lines_panel.js |
| `getPCRender` | 1 | page1 | ?? — let merge chat decide |
| `getSampleColor` | 1 | page1 | ?? — let merge chat decide |
| `has` | 1 | page12 | ?? — let merge chat decide |
| `hubs` | 1 | page1 | ?? — let merge chat decide |
| `jittered` | 1 | page1 | ?? — let merge chat decide |
| `kColor` | 1 | page12 | ?? — let merge chat decide |
| `layer` | 1 | page1 | ?? — let merge chat decide |
| `mbAt` | 1 | page1 | ?? — let merge chat decide |
| `niceTicks` | 1 | page1 | ?? — let merge chat decide |
| `palette` | 1 | page12 | shared/draw_utils.js |
| `persistCandidateList` | 1 | page2 | shared/candidate_helpers.js |
| `q` | 1 | page12 | ?? — let merge chat decide |
| `recomputeAnchorConcord` | 1 | page1 | ?? — let merge chat decide |
| `refreshCandidateListUI` | 1 | page2 | shared/candidate_helpers.js |
| `refreshCandidateUI` | 1 | page2 | shared/candidate_helpers.js |
| `renderCatalogue` | 1 | page2 | ?? — let merge chat decide |
| `samples` | 1 | page1 | ?? — let merge chat decide |
| `showHide` | 1 | page12 | ?? — let merge chat decide |
| `sigmaProfileCandidate` | 1 | page2 | shared/candidate_helpers.js |
| `simColor` | 1 | page1 | ?? — let merge chat decide |
| `simColorPDF` | 1 | page1 | shared/draw_utils.js |
| `strokeSamplePath` | 1 | page1 | ?? — let merge chat decide |
| `strokeSamplePathStyled` | 1 | page1 | ?? — let merge chat decide |
| `themeColor` | 1 | page1 | ?? — let merge chat decide |
| `toPx` | 1 | page1 | ?? — let merge chat decide |
| `toPy` | 1 | page1 | ?? — let merge chat decide |
| `trackedColor` | 1 | page1 | shared/tracked_panel.js |
| `wireCandidateAncestryConfound` | 1 | page2 | shared/candidate_helpers.js |
| `wireCandidateButtons` | 1 | page2 | shared/candidate_helpers.js |
| `withAlpha` | 1 | page1 | ?? — let merge chat decide |
| `xAt` | 1 | page12 | ?? — let merge chat decide |
| `xOfWin` | 1 | page1 | ?? — let merge chat decide |
| `xToPx` | 1 | page12 | ?? — let merge chat decide |
| `yAt` | 1 | page12 | ?? — let merge chat decide |
| `yToPx` | 1 | page12 | ?? — let merge chat decide |
| `zColorPDF` | 1 | page1 | shared/draw_utils.js |

**Total unique missing functions: 109.**

Regenerate this list any time with:
```bash
grep -rh 'TODO_MISSING(' inversion_discovery/page*.js | sort | uniq -c | sort -rn
```

### State slots referenced but not in shared/state.js SLOT_REGISTRY

| Slot | Pages | Likely class |
|---|---|---|
| `state._simGeom` | page1, page12 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state._l3CacheFp` | page1 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state._l3CacheRendered` | page1 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state._lineageComputeScheduled` | page1 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state._simMinimapGeom` | page1 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state._thSimGeom` | page12 | ad-hoc geometry / cache (legacy added imperatively at first use) |
| `state.ancestryPalette` | page1 | ? |
| `state.cacheKey` | page1 | transient cache |
| `state.candidateMode` | page1 | persisted UI mode (probably needs to be registered) |
| `state.candidatePageMode` | page2 | ? |
| `state.colorMode` | page1 | persisted UI mode (probably needs to be registered) |
| `state.compareUnit` | page1 | persisted UI mode (probably needs to be registered) |
| `state.crossSpecies` | page1 | ? |
| `state.hubFamilies` | page1 | ? |
| `state.l2GroupCache` | page1 | transient cache |
| `state.l3Draft` | page1 | ? |
| `state.l3Mode` | page1 | ? |
| `state.l3ReclusterMode` | page1 | ? |
| `state.l3SecondaryMetric` | page1 | ? |
| `state.scaleStabilityPanes` | page1 | ? |
| `state.schemaVersion` | page1 | ? |
| `state.secondaryL2` | page1 | ? |
| `state.simInMinimap` | page1 | ? |
| `state.singletonFamilyIds` | page1 | ? |
| `state.smallFamilyIds` | page1 | ? |
| `state.trackedN` | page1 | persisted limit (probably needs to be registered) |
| `state.windowToL1` | page1 | ? |
| `state.windowToL2` | page1 | ? |
| `state.zColorMode` | page1 | ? |
| `state.zHighlightThr` | page1 | ? |
| `state.zValueMode` | page1 | ? |

**Total unique missing slots: 31.**

Slots prefixed with `_` (e.g. `_simGeom`, `_l3CacheFp`) are ad-hoc geometry/cache
writes the legacy added at first use; merge chat should either register them in
SLOT_REGISTRY (with `tag: 'transient'`) or refactor the writers to keep them off
the state object.

## Decisions made

- **Brace-matched function slicing.** A naïve "next function start" slicer pulled
  in inter-function IIFEs (e.g. the PCA-lasso wiring block between `onPCAClick`
  and `togglePlay`), which then ran at module-import time and threw
  `ReferenceError: document is not defined` under Node. Switched to a brace-
  counter that respects strings, template literals, and comments. Net effect:
  module bodies are tighter and TODO_MISSING dropped by ~14 in page1.

- **`state` as first parameter.** Every extracted function takes `state` as its
  first arg, matching the `shared/per_l2_cluster.js` `contextFromState(state)`
  pattern. Function bodies are otherwise LITERAL — no body refactoring.

- **Shared imports declared in every page header**, even when unused. ESM is
  fine with unused named imports, and it lets the merge chat add references
  without re-touching the import line.

- **Page8 + page19 have NO extracted JS.** Confirmed by full-file grep that
  `winSum*` and `nr*` IDs are referenced ONLY from CSS / HTML / comments in
  legacy — no `getElementById` calls, no event wiring. They are pure scaffold.
  Their renderers will be authored fresh by the merge chat or a follow-up batch.

- **Page15 has only `_refreshGhslLayerStatus`.** The page is mostly empty-state
  documentation about layers that don't yet exist; only the layer-pill badge
  has live JS.

- **Did not extract page2's ~25 helper functions** (`candidateNavHtml`, `candidateBandsHtml`,
  `sigmaProfileCandidate`, `wireCandidateButtons`, etc). Per Rule 4, these are
  flagged TODO_MISSING and left for the merge chat (most are likely reused by
  other discovery / catalogue pages and belong in `shared/candidate_helpers.js`).

## What's NOT done

- **Tab routing** (which page activates when, how state.activePage flows) —
  left for merge chat.
- **CSS** — left for merge chat. The HTML fragments still reference legacy
  CSS variables (`var(--ink)`, `var(--rule)`, etc) which assume a global
  stylesheet.
- **The 57 + slot-level TODO_MISSING references.** Some will resolve once
  batches 2–5 land; the rest need to be hoisted into `shared/` or authored
  fresh.
- **No behavioural tests.** The new tests only verify that modules load and
  exports are present. Behavioural tests need DOM + state and depend on the
  merge chat's wiring decisions.
- **Page8 and page19 renderers** — pure scaffold in legacy; need fresh
  authoring.
