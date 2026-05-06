# HANDOFF — turn 157 — dosage shim load fix + page6/page7 try/catch (PARTIAL)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (73,493 lines, +58 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

**Honest framing up front**: Quentin reported 8 distinct issues in one
message. This turn ships fixes for #1 (dosage shim load) and #2/#3
(diagnostic guard + likely root-cause for "popstats tracks not
computed when server connected"). Issues #4–#8 are real features that
each warrant their own turn — captured in §8 with concrete scoping.

---

## 0. What Quentin reported

> "The dosage.js shim script is not found for some reasons? Can we
> fix it?
>
> When we click ancestry tab or popstats it doesn't allow us to click
> on other tabs, it doesn't refresh the page.
>
> When the popstats server is connected the tracks are not computed.
>
> It should be able to directly use the groups in the G bar and I-g /
> candidate karyotypes.
>
> We separate the tracks on the right hand side selection for QC and
> for popstats.
>
> And turn off popstats track by default but QC tracks are off by
> default.
>
> The LD track doesn't allow to plot double heatmaps for common versus
> rare karyotype groups.
>
> We need to be able to export the plot as well.
>
> If you cannot manage it just add to a handoff."

Triaged into 8 numbered items in §8.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse, pragmatic.
Tarball is the standard handoff format.

---

## 2. What shipped this turn — issue #1 (dosage shim load)

### 2.1 Root cause

`js/atlas_dosage_bridge.js` exists on disk. The file's responsibility
is to install a synthetic `state.data.dosage_chunks` index whose URLs
template into the popstats server's `/api/dosage/chunk` endpoint. The
atlas-side `_resolveChunkUrl` (line 16774, pre-existing) substitutes
the placeholders at fetch time.

But the `<script src="js/atlas_dosage_bridge.js">` tag was missing
from the trailing script-load block (the registry comment at line
52184 even noted "Bundled but not auto-loaded"). So
`window.AtlasDosageBridge` was never defined. So the heatmap, the
stripe-quality compute, AND the turn-156 V-shape diagnostic could not
pull live dosage data — they stayed stuck on whatever static index
was in `state.data.dosage_chunks` (typically nothing in live
sessions until precomp loads).

This is also the most likely root cause of issue #3 ("popstats
server connected but tracks not computed"). Several popstats tracks
ultimately require dosage chunks — if the bridge isn't installed,
the chunks-by-region lookup returns null, downstream renderers
gracefully degrade to the empty-state, and the user sees "tracks
not computed."

### 2.2 The fix

```html
<script src="js/atlas_dosage_bridge.js"></script>     <!-- turn 157 -->
<script>
  // Auto-install on every positive health probe.
  (function _wireDosageBridgeAutoInstall() {
    if (!window.atlasServer || !window.AtlasDosageBridge) return;
    const orig = window.atlasServer.isAvailable;
    if (window.atlasServer._dosageBridgeWired) return;
    window.atlasServer._dosageBridgeWired = true;
    window.atlasServer.isAvailable = async function (forceRefresh) {
      const ok = await orig.call(this, forceRefresh);
      if (ok && window.AtlasDosageBridge) {
        try {
          window.AtlasDosageBridge.install({
            atlasServer: window.atlasServer,
            stateRef: window.state,
          }).catch((e) => console.warn('[turn 157] install failed:', e));
        } catch (e) {
          console.warn('[turn 157] install threw:', e);
        }
      }
      return ok;
    };
  })();
</script>
```

Properties:
- **Idempotent** — `_dosageBridgeWired` flag prevents double-patching.
  `AtlasDosageBridge.install()` itself is also idempotent (its own
  module-level docstring says so).
- **Non-blocking** — `install()` is async; the `await` resolves before
  the install promise; `isAvailable` returns the original `ok` value.
- **Silent on partial-state** — server up but dosage subsystem disabled
  is a legal state (the turn 145 partial-state pattern). The install
  rejects with `reason: 'manifest_unavailable'` and the heatmap falls
  back to whatever static index is in place.
- **Auto-runs on EVERY isAvailable call** — chrom switches, page
  navigation, and the initial health probe all trigger a re-install
  attempt. The bridge tracks active chrom internally and updates the
  index when the active chrom changes.

---

## 3. What shipped this turn — issue #2 (UI freeze on tabs)

### 3.1 Diagnostic guard (not a root-cause fix)

The popstats and ancestry tab click handlers schedule their renders
via `requestAnimationFrame`. If the render throws, the user sees a
half-rendered page with no console output and (depending on what was
in flight at the time) potentially broken DOM that interferes with
subsequent clicks.

Pre-turn-157:
```js
if (target === 'page6') {
  requestAnimationFrame(() => {
    renderPopstatsPage();   // throws → silent half-render
  });
}
```

Post-turn-157:
```js
if (target === 'page6') {
  requestAnimationFrame(() => {
    try { renderPopstatsPage(); }
    catch (e) {
      console.warn('[turn 157] renderPopstatsPage threw:',
        e && e.message ? e.message : e, e && e.stack);
    }
  });
}
```

This does NOT fix the root-cause throw if there is one — but it gives
us a precise console message Quentin can paste back next session, and
it prevents the failed render from poisoning the click handler's
state on subsequent tab clicks (the click handler itself is unchanged
and still completes normally).

### 3.2 What I couldn't fix without browser-side reproduction

The actual freeze cause depends on what `renderPopstatsPage` /
`renderAncestryPage` are throwing on. Most likely candidates:
- Reading from `state.data.tracks.X` when X is missing
- A canvas resize race (`fitCanvas` against zero-dimension hidden
  canvas)
- An installed listener that grabs the click event stream

Without a reproducible case + console error, fixing the root cause is
guesswork. Ask Quentin to paste the `[turn 157] renderPopstatsPage
threw:` warning content next session — that pins it down precisely.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 156)   | 73,435  | 2840            | 60    |
| Post-turn-157 (current)  | 73,493  | 2865            | 61    |
| **Δ session**            | +58     | +25             | +1    |

Full sweep: **2865 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn157_dosage_bridge_load.js` (25 tests) covers:

1. Dosage bridge `<script src>` tag present + ordered after request layer
2. Auto-install hook IIFE: patches `atlasServer.isAvailable`, idempotent
   via `_dosageBridgeWired`, calls `install` with stateRef + atlasServer,
   only on `ok=true`, catches both promise rejections and synchronous
   throws
3. page6/page7 RAF guards: try/catch present with console.warn
4. Existing flow preserved (atlasServer, render functions,
   `_resolveChunkUrl`, turn 155 + 156 contracts)
5. Bridge file exists on disk, exports `AtlasDosageBridge`,
   `install()` is async

No existing test needed inverting.

---

## 5. What Quentin can verify after this turn

1. **Open the popstats server** (`./run_server.sh` or equivalent on the
   atlas server module).
2. **Reload Inversion_atlas.html.** Open browser console.
3. **Look for** the absence of any `[turn 157] dosage bridge install
   failed` warning. If you see one, paste it back — that means the
   bridge module is present but the install rejected (most likely
   reason: server not up at health probe time, or dosage subsystem
   disabled in `popstats_server.config.yaml`).
4. **Click page 2 (candidate focus)** with a candidate selected.
   **Open the dosage heatmap.** It should now populate from the live
   server. If it still says "no chunks," paste any console messages.
5. **Click popstats tab.** If it still freezes, paste the
   `[turn 157] renderPopstatsPage threw:` warning from console — that
   will pin the root cause for the next turn.

---

## 6. Things I did NOT change

- **`renderPopstatsPage` / `renderAncestryPage` bodies** — only
  wrapped their RAF entry-points in try/catch.
- **`atlasServer.isAvailable` body** — only wrapped the original via
  the auto-install hook.
- **`AtlasDosageBridge` source** — not modified.
- **`_resolveChunkUrl`** — not modified.
- **Anything in turns 152–156** — fully preserved.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 157 patched, 73,493 lines.
- `js/atlas_dosage_bridge.js` — pre-existing, now actually loaded.
- `tests/test_turn157_dosage_bridge_load.js` — NEW, 25 tests.
- `HANDOFF_2026-05-05_turn157_dosage_bridge_load.md` — this file.

Plus prior handoffs (carried for reference): turn 156, 155, 154, 153,
152.

---

## 8. Queued: issues #4–#8 — what each one needs

These are not getting done this turn. Each is a real change with real
scope. Listing them with concrete proposals so the next session can
pick one up cold.

### #4 — Use G-bar groups + I·g + candidate karyotypes for popstats groupings

**Status**: Real feature. Medium-large.

**What today does**: Popstats tracks consume static per-group dosage
data derived at precomp time. Groupings come from whatever groups the
upstream R pipeline (`STEP_Q07_popstats.sh`) wrote.

**What Quentin wants**: Live re-grouping. When you change the G-bar
selection (manual / karyotype / inheritance) on page 1, the popstats
tracks on page 8 should re-compute against the new group definitions.

**Where to wire**:
- `state.gPanelInheritanceThreshold` or any G-panel mutation
  → ideally fires the same `_notifyInheritanceConsumers` (turn 153)
  pattern, with one new consumer: a popstats live-recompute path.
- The popstats live-recompute must POST to `/api/popstats/compute` with
  the new group membership (CGA → group_id mapping). The server-side
  endpoint exists per `js/atlas_request_layer.js`'s `popgenLive` global.
- Track defs in `POPSTATS_TRACKS` need a per-track flag like
  `livePushSupported: true` so renderers know to subscribe to the live
  stash vs the static precomp.

**Estimated**: 200-400 LOC + significant test surface.

**Risk**: server endpoint contract may not match what the page expects.
Worth verifying the request shape end-to-end on a small toy candidate
set before scaling.

### #5 — Separate QC vs popstats track sections in the right-hand selector

**Status**: Real feature. Medium.

**What today does**: `POPSTATS_TRACKS` is a single flat array.
`renderPopstatsPage`'s chip toolbar renders all of them as a single row.
There's no QC/popstats distinction.

**What Quentin wants**: Two grouped sections in the chip toolbar — one
for "QC" tracks (BEAGLE uncert, coverage, low-cov count, SNP density,
sim_mat, ideogram) and one for "popstats" tracks (Z, θπ, Fst, Hobs/Hexp,
Δ12).

**Where to wire**:
- Add `category: 'qc' | 'popstats'` to each entry in `POPSTATS_TRACKS`.
- `renderPopstatsPage` builds two chip rows, one per category, with
  section labels.

**Estimated**: ~80 LOC + tests.

### #6 — Default ON/OFF for QC vs popstats tracks

**Status**: Tied to #5. Small once #5 is in.

**What Quentin wants**:
- popstats tracks: **OFF by default** (currently most are ON when data
  is present)
- QC tracks: **OFF by default** (also currently ON when data present)

Implementation: per-category default in the visibility check:
```js
const visible = (t) => {
  if (hiddenStored.has(t.id)) return false;
  if (t.alwaysOn) return true;
  // turn 158: per-category default — both QC and popstats off by default
  if (t.category === 'qc' || t.category === 'popstats') {
    return hiddenStored.has(t.id) ? false : false;   // off until user enables
  }
  if (!t.hasData) return false;
  return true;
};
```

**Caveat**: This inverts the "data present → visible" default that
exists today. Existing user storage state will be honored, so users
who'd already set their preferred chips will keep their state. New
users will start with everything off, which might feel empty. Worth
adding a "show all" / "hide all" pair of buttons next to the toolbar
to make first-time setup quick.

**Estimated**: ~40 LOC.

### #7 — LD double-heatmap: common vs rare karyotype groups

**Status**: Real feature. Medium-large.

**What today does**: `atlas_ld.js` renders one combined heatmap for
the whole cohort.

**What Quentin wants**: Two side-by-side heatmaps:
- Left: LD computed from the **common** karyotype group (typically the
  HOMO band with ≥ N samples)
- Right: LD computed from the **rare** karyotype group (typically the
  HET / minority HOMO band)

Difference reveals LD that's specific to the inversion polymorphism
vs background LD.

**Where to wire**:
- `atlas_ld.js`'s `makeLDSplitPanel` already exists per the registry
  (line 7282 mentions it). Likely needs a new mode flag
  `splitMode: 'common-vs-rare-karyo'`.
- Group definition: cohort split by `cand.locked_labels` (band-size
  threshold to decide common vs rare — propose minor band, default the
  rest as common).
- The existing fast_ld endpoint accepts a sample subset. Two parallel
  POSTs with different subsets, render both, share the same LD scale.

**Estimated**: 300-500 LOC, mostly UI + request orchestration.

### #8 — Plot export

**Status**: Real feature. Medium.

**What today does**: Most plots are canvas-based. Some let the user
right-click → save image; many do not.

**What Quentin wants**: Explicit "Export" button on each plot —
download as PNG (canvas → toDataURL) or SVG (where the renderer is
SVG-based).

**Where to wire**:
- For canvas plots: a single shared helper `_exportCanvasAsPNG(canvas,
  filename)` using `canvas.toDataURL()` + an anchor click. Add an
  export button to each plot's toolbar.
- Plots that need SVG: the inheritance matrix, the V-shape, the LD
  heatmap, the per-sample lines panel. Each has its own canvas — wrap
  them.

**Estimated**: 100-200 LOC of helpers + per-plot toolbar buttons.

---

## 9. Honest framing

**What's solid:**
- Issue #1 (dosage shim) is straightforwardly fixed. The bridge is now
  loaded and installs automatically on server availability.
- Issue #2 (UI freeze) gets a diagnostic guard that won't make it
  worse, will prevent it from poisoning subsequent tab clicks, and
  surfaces the actual error to console for the next pass.
- Issue #3 (tracks not computed) is most likely resolved as a
  side-effect of #1, but I can't verify without a live server run.
  Quentin should test and report back.

**What's open:**
- The actual root cause of issue #2 — needs browser-side reproduction.
- Issues #4–#8 — each is its own turn, scoping in §8.

**What's queued (in priority-order I'd suggest):**
- **Verify #1 + #3 actually fix what Quentin reported** (next session,
  needs Quentin's hands-on)
- **#6 (defaults)** — tiny, decoupled, immediate clarity win.
- **#5 (QC/popstats split)** — also small, makes #6 cleaner.
- **#8 (plot export)** — bounded; useful for the manuscript figures.
- **#4 (live re-grouping)** — biggest, riskiest, most valuable.
  Worth it but should come after #5/#6 land so the toolbar is
  already rationalised.
- **#7 (LD double-heatmap)** — independent of all the above; pick
  whenever the LD path is hot.

End of handoff.
