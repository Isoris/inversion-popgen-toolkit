# Inversion Atlas — chat-A integrated bundle

Patched `Inversion_atlas.html` plus 11 chat-A JavaScript modules,
spliced in via the per-turn `SPLICE_POINTS.md` instructions.

## Layout

Drop everything in this tarball into the same directory. The HTML
references the JS modules by filename (relative path), so they must
sit alongside.

```
Inversion_atlas.html        ← patched host atlas
atlas_group_engine.js       ← turn 2  (popgen group composition engine)
atlas_request_layer.js      ← turn 3  (popgenLive request layer + IDB cache)
atlas_renderers_turn4.js    ← turn 4  (popstats track renderers)
atlas_floating_dock.js      ← turn 5a (floating dock UI)
atlas_dock_tabs.js          ← turn 5b (dock Groups / Lasso / Snapshots tabs)
atlas_page6_wiring.js       ← turn 6  (page 6 popstats wiring)
atlas_track_gallery.js      ← turn 6.5 (track gallery sidebar)
atlas_turn7.js              ← turn 7  (Q09b shelf-LD on page 3 + Δ12 stratifier)
atlas_overview.js           ← turn 9  (candidate overview spreadsheet — new tab "15b")
atlas_compare.js            ← turn 10 (set-arithmetic comparison library)
atlas_ld.js                 ← turn 11c (fast_ld LD split-heatmap on page 3)
```

## Load order (applied to the HTML, do not reorder)

```html
<script src="js/atlas_group_engine.js"></script>      <!-- turn 2 -->
<script src="js/atlas_request_layer.js"></script>     <!-- turn 3 -->
<script src="js/atlas_renderers_turn4.js"></script>   <!-- turn 4 -->
<script src="js/atlas_floating_dock.js"></script>     <!-- turn 5a -->
<script src="js/atlas_dock_tabs.js"></script>         <!-- turn 5b -->
<script src="js/atlas_page6_wiring.js"></script>      <!-- turn 6 -->
<script src="js/atlas_track_gallery.js"></script>     <!-- turn 6.5 -->
<script src="js/atlas_turn7.js"></script>             <!-- turn 7 -->
<script src="js/atlas_overview.js"></script>          <!-- turn 9 -->
<script src="js/atlas_compare.js"></script>           <!-- turn 10 -->
<script src="js/atlas_ld.js"></script>                <!-- turn 11c -->
```

## What the splices changed in `Inversion_atlas.html`

Surgical edits only. Net +565 lines. Highlights:

- **Page 6 (popstats)**: `collectPopstatsTracks()` now routes through
  `popgenPage6.wrapAllTrackDefs` (turn 6) after a turn-6.5 gallery
  union/filter pass. The `drawPopstatsTracks()` loop dispatches through
  `popgenRenderers.dispatch` first (turn 4). Hover tooltip wiring
  added per-canvas. Toolbar gained `mode: static / live / split` chip
  (retro-pass). Multiline tracks gained per-track `curves / fill /
  heatmap` chip in their headers (retro-pass). New right-side
  `#psGalleryTray` (turn 6.5) holding the track gallery.
- **Page 3 (candidate focus)**: two new mount slots —
  `#page3_q09b_slot` (turn 7) and `#page3_ld_slot` (turn 11c). Q09b
  shelf-LD test panel mounts inside, plus the fast_ld split heatmap
  with HOM1/HOM2/HET picker.
- **New page**: `data-page="page_overview"` tab (num "15b") added
  between marker panel and help. Page div mounts the candidate
  spreadsheet on init.
- **Globals exposed**: `window.popgenDosage`, `window.renderPopstatsPage`,
  `window.renderAncestryPage` so the chat-A modules can find them.
- **Retro-pass monkey-patches** at the very bottom of body:
  monkey-patches `renderPopstatsPage` (track-header chip injection),
  `renderTrackedList` (`popgen.recordSelection` on every tracked
  change → dock Lasso tab populates), and `renderCandidateMetadata`
  (`popgen.invalidateDimensions` + `popgenDockTabs.refresh('groups')`
  on focal-candidate switch).

All edits are guarded with feature-detect checks
(`if (window.popgenX && typeof window.popgenX.method === 'function')`)
so removing any chat-A `<script>` tag degrades gracefully back to the
host's pre-splice behavior.

## Server side

This tarball is **atlas-only**. The fast_ld endpoints, popstats
endpoints, ngsLD-replacement engine, and configuration files live in
the companion `server_bundle.tar.gz` plus `DEPLOY_TO_LANTA.md`.

Without a server, the LD panel and any "Compute" actions in the dock
will fail with `error: fetch ...` in the status line. Everything else
(static-mode tracks, page navigation, dock UI, snapshots, track
gallery toggling, candidate overview spreadsheet, Q09b shelf-LD test
which is pure-JS) works offline.

## Verifying

Open `Inversion_atlas.html` directly in a browser.

Expected: page 1 renders normally; small "G" launcher button at
top-left; toolbar on page 6 shows the new mode chip; new "overview"
tab in the tab bar between marker-panel and help; page 3 shows Q09b
+ LD panels below the catalogue table when a focal candidate is set.

Console should be free of errors. `window.popgen`, `window.popgenLive`,
`window.popgenRenderers`, `window.popgenDock`, `window.popgenDockTabs`,
`window.popgenPage6`, `window.popgenGallery`, `window.popgenTurn7`,
`window.popgenOverview`, `window.popgenCompare`, `window.popgenLD`
should all be defined.

Per-module unit tests live in the source bundle (chat A's
`atlas_turn*` directories), each runnable via `node test_*.js`.

## Open items not addressed by splicing

These were either flagged as future work in chat A or require
out-of-band wiring the user owns:

- A multi-chrom `state.atlas_catalogue` is not auto-fetched. The
  overview falls back to `state.data.candidates` (current chrom).
  To populate cross-chrom, `state.atlas_catalogue = [...]` from
  whatever loader you prefer.
- `popgenCompare` is exposed as a function library; no inline UI
  affordance was wired (would have required editing the tested
  `atlas_dock_tabs.js` / `atlas_overview.js` modules, which the
  splice protocol forbids). Use from the browser console:
  `popgenCompare.compareSlotsByName('HOM1', 'HOM2_LG14', { universe_size: 'cohort' })`.
- Stratified Δ12 resolver (`window.popgenPage7DeltaResolver`) is
  installed but the host's `renderAncestryPage` still reads from
  `state.data.tracks.ancestry_delta12`. Wiring the resolver into
  the live tracks path is non-trivial and was deferred — if a
  per-group ancestry compute response is loaded, the resolver is
  ready to be consulted.
