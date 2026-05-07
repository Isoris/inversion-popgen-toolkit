# Batch 4 — extraction notes

**Pages**: page10 (Marker panels), page_overview (synthesis stub).
**Sub-atlas**: `inversion_catalogue` (part 2).
**Shared/ touched**: none. All 7 foundation regression suites still pass.

## Files produced

```
Atlas/inversion_catalogue/
├── page10.html            (23 lines, legacy 7873–7895)
├── page10.js              (244 lines — wraps legacy 57837–58043)
├── page_overview.html     (1 line, empty stub from legacy 9322)
├── page_overview.js       (37 lines — empty-render stub)
└── BATCH_4_NOTES.md

Atlas/tests/
├── test_catalogue_page10.js          (10 assertions, all pass)
└── test_catalogue_page_overview.js   ( 5 assertions, all pass)
```

## Page 10 — "Marker panels"

**Discrepancy with handoff table**: HANDOFF_BATCH_4.md lists page10 as
"Help — keyboard shortcuts, dependency graph, layer descriptions". The
actual content of legacy lines 7873–7895 is the **Marker panels** page
(diagnostic PCR marker panels per candidate, gated on
`marker_panel_summary` layer presence, see SCHEMA §10). Confirmed by:

- Legacy `<h2>` reads "Marker panels".
- Subtitle reads "(no marker layer loaded)" when gated off.
- The render function is `renderMarkerPage` (legacy line 57972).
- No "Help" panel exists at this `<div id="page10">` location.

The handoff table description was incorrect; the line range was right.

**Render function**: `renderMarkerPage` (legacy 57972–58043, 72 lines).
**Helper**: `_markerPanelCardHtml` (legacy 57837–57970, 134 lines) is
used **only** by `renderMarkerPage` (verified: 2 grep hits — definition
+ single call site). Because it has no other consumer, it was extracted
together with the render function as a module-private helper rather
than being marked `TODO_PROMOTE_TO_SHARED`. The merge chat can hoist it
to shared/ later if any other page turns out to need it.

**External references**: only `state.data._layers_present`,
`state.data.marker_panel_summary`, `state.data.marker_catalogue`,
`state.data.marker_primers`, and `state.candidateList` — all of which
are part of the canonical `state` contract managed by `shared/state.js`.
DOM access is limited to `getElementById('page10Content')` and
`getElementById('page10Subtitle')`.

**No TODO_MISSING markers were emitted for page10** — the function is
self-contained once `state` is in hand.

## page_overview — empty stub

**Legacy state**: line 9322 is literally `<div id="page_overview" class="page"></div>`
with no body. The corresponding tab button at legacy line 5138 declares
`data-page="page_overview"` `data-stage="synthesis"`.

**JS render function**: none exists in legacy. Verified by:

```
grep -niE "(renderOverview|page_overview|renderPageOverview)" legacy/Inversion_atlas.html
```

returning only the tab button (line 5138) and the empty div (line 9322).
No render function, no helpers, no event wiring. The synthesis-stage
overview was scoped in the tab nav but never implemented.

**Shipped**: an empty `wirePageOverview(state)` stub that returns a
no-op `renderPageOverview()`. Kept so the page registry has a
non-throwing entry. The merge chat (or a future design pass) decides
whether to drop the tab or actually populate it.

## TODO inventory

```
TODO_MISSING(*)                       count
─────────────────────────────────────  ─────
synthesis_overview_design                 1   (page_overview.js, design-level — not a missing function)

TODO_PROMOTE_TO_SHARED(*)             count
─────────────────────────────────────  ─────
(none)
```

The single marker is a design hold, not a missing-function trace; the
merge chat can ignore it unless the synthesis overview is being put
into scope.

## Test results

```
test_catalogue_page10.js          pass: 10   fail: 0
test_catalogue_page_overview.js   pass:  5   fail: 0
```

## Foundation regressions

All shared + modular smoke suites still green:

```
test_shared_contingency.js     OK
test_shared_het_rate.js        OK
test_shared_hungarian.js       OK
test_shared_kmeans.js          OK
test_shared_per_l2_cluster.js  OK
test_shared_state.js           OK
test_modular_smoke.js          OK
test_legacy_parity.js          (skipped — legacy not at hardcoded path; not affected by this batch)
```

`shared/` has not been modified.

## Notes for the merge chat

1. **`_markerPanelCardHtml` is local** — if any of batches 1/2/3/5 also
   need this helper (unlikely, the name is highly specific), it can be
   hoisted to `shared/marker_panels.js`. If it stays single-use, leave
   it where it is.
2. **page_overview is empty by design** — don't flag it as a missing
   render. Either delete the tab from the new build or wait for a
   synthesis-stage UX spec.
3. **The handoff table mis-described page10** as "Help"; the page is
   actually "Marker panels". If "Help" still needs to be extracted,
   it lives somewhere else in the legacy file and was not in batch 4's
   line range.
