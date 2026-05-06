# HANDOFF — turn 162 — band-trace tooltip + TSV export

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (75,194 lines, +401 LOC from 74,793)
**Working dir**: `/home/claude/work/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes** items 2 (hover tooltip) and 4 (TSV export) from the turn 161
handoff §6 NEXT list. Items 1 (run brackets) and 3 (combinatorial
enumeration) remain deliberately deferred.

**Picked up from**: post-turn-161 working tree (3197 / 0 baseline,
77 / 0 on `test_turn161_band_trace_ui.js`).

---

## 0. Manuscript framing carried forward

Same axiom as turn 161: 30 family hubs across 226 fish ⇒ effectively
random mixing of founder backgrounds, so co-segregation at distant L2s
is a real physical observation, not population structure. This turn
adds two surfaces that **make the observation more legible without
crossing the interpretation line**:

- **Tooltip** = "what does this stripe mean numerically?"
- **TSV** = "supplementary table for the figure / pipe to R."

Neither surface labels anything as an inversion. The tooltip just
reports `entropy / dominant_band / fraction / regime / n_valid` for the
hovered L2; the TSV is one row per L2 with the same numerics plus
chromosome coordinates. Quentin reads the pattern.

The **brackets-for-detected-runs** feature (turn 161 §6 item 1) is
still deferred because it crosses into "this is an inversion"
interpretation; per turn 161, only add after Quentin has seen real
LG28 data and asks for it.

---

## 1. What this turn ships

### Hover tooltip on the band-trace strip

Six new functions, all in a fresh `// turn 162 …` section just after
the `// end turn 161 UI layer` marker:

- `_btraceTooltipEnsureEl()` — lazy-creates a fixed-position tooltip
  div with id `btracePointTooltip`. Same visual style as `_inhTooltip`
  (turn 122) and `_vsTooltip` (turn 157A) — dark panel, mono font,
  edge-clamped via `_btraceTooltipShow`.
- `_btraceTooltipBuildHtml(hit, opts)` — pure HTML builder. Renders:
  - Header line: `L2 #<idx>` + regime chip in regime-color (green /
    amber / grey / dark-grey / very-dark)
  - Optional chain badge `chain <ci>` only when `opts.n_chains > 1`
  - Body table: `n_valid: K / N`, optional `dominant: <swatch> b<k>
    (xx.x%)` row (skipped when `dominant_band === -1`),
    `entropy: 0.213` (3 dp)
  - Per-band fractions row: one swatch + `b<k> NN%` per non-zero
    fraction (skipped when `regime === 'no_valid'`)
- `_btraceTooltipShow(hit, clientX, clientY, opts)` — populate +
  position with edge-clamping (`vw - 8`, `vh - 8` margins); falls back
  to right+below cursor by default.
- `_btraceTooltipHide()` — display:none on the el if present.
- `_btraceHitTest(hits, x, y)` — returns the first rect containing
  `(x, y)`, null on miss. Hits are mutually disjoint (one per L2
  column), so first-hit is also only-hit.
- `_wireBandTraceTooltip(canvas)` — installs `mousemove` +
  `mouseleave` on the canvas, idempotent via
  `canvas.dataset.btraceTooltipWired === '1'`. Reads `state._btraceHits`
  and `state.bandTraceCache` on every move (no closure capture of
  stale data; re-renders that mutate state are picked up immediately).

### Hits stash on `_drawBandTraceStrip`

One new state slot, populated by the existing strip drawer:

- `state._btraceHits : array | null` — one
  `{ x, y, w, h, l2_idx, entry }` per visible L2 column. Coordinates
  are CSS pixels (the canvas is `setTransform(dpr,…)`-ed by `fitCanvas`
  before drawing, so the strip rectangles drawn at `pad.l + …`,
  `pad.t - 21` are in CSS-pixel space). Each rect spans the full 7 px
  strip height for a generous hover target. `entry` is the trace's
  `per_l2[i]` record, so the tooltip renders without re-deriving
  anything.

`_bandTraceClearCache` was extended to null `state._btraceHits` too,
so chrom switches don't leave stale rectangles behind.

### TSV export

- `_bandTraceToTSV(trace, opts)` — pure serializer. Header row +
  one row per `per_l2` entry. Columns: `chrom, l2_idx, chain_idx,
  chain_position, start_bp, end_bp, n_valid, n_fish_selected, regime,
  dominant_band, dominant_fraction, entropy, band_fraction_0..K-1`.
  K is read from `trace.K` (not from `band_fractions.length` per row).
  Floats formatted to 6 dp; NaN/null → empty. `start_bp`/`end_bp`
  filled from `opts.envelopes[l2_idx]` when provided, blank otherwise.
- `_bandTraceDownloadTSV()` — pulls the cached trace via
  `_bandTraceGetOrCompute`, calls `_bandTraceToTSV` with current
  chrom + envelopes, then performs a Blob/URL/anchor download.
  Filename: `band_trace_<chrom>_n<n_fish>_K<K>.tsv`. Headless-safe:
  when `document` or `Blob` are unavailable, returns the filename
  without performing the download (so tests can verify naming logic
  without simulating a real DOM). Wraps the actual download in
  try/catch so a blocked browser path doesn't crash the lines panel.

### Lines header UI

One new button next to the existing `🔍 trace`:

```
[ ] band trace    🔍 trace    📊 TSV
```

- `linesBandTraceExportBtn` (📊 TSV): same pill style as `🔍 trace`.
- Click handler: calls `_bandTraceDownloadTSV()`. On null return
  (no fish-set / no cached trace), surfaces a friendly alert
  (`"No active band-trace to export. Click 🔍 trace first."`).
  Does not throw.

### `drawLinesPanel` integration

One additional line in the per-source canvas creation loop, right next
to the existing `_wireInheritancePillTooltip` install:

```js
if (src === 'pc1' && typeof _wireBandTraceTooltip === 'function') {
  try { _wireBandTraceTooltip(cv); } catch (_) {}
}
```

PC1-only because the strip itself is PC1-only. Idempotent via the
canvas dataset marker, so re-running `drawLinesPanel` (which destroys +
recreates these canvases on every paint) just attaches handlers to the
fresh nodes without leaks.

### Coordinate-space decision (worth recording)

The codebase has **two patterns** for canvas tooltips:

1. **`_inhTooltip` (turn 122)**: scales CSS cursor coords to
   canvas-internal coords via `canvas.width / rect.width`. Stores hits
   in… also CSS pixels per the source comment. This is **latently
   broken on dpr > 1** — the multiplication doesn't cancel because the
   regions are CSS pixels, not canvas-internal. Works in test envs
   (dpr = 1) and on standard monitors; would miss-fire on a Retina
   display.
2. **`_vsTooltip` (turn 157A)**: CSS pixels throughout, no scaling.
   Comment explicitly notes "hits are produced AFTER setTransform, so
   they're already CSS-pixel coordinates."

I deliberately followed pattern (2). The lines panel canvases get
`setTransform(dpr, 0, 0, dpr, 0, 0)` from `fitCanvas`, so all draw
coordinates (`pad.l`, `pad.t - 21`, `xLo`, `w`, `stripH`) are CSS
pixels. The mousemove handler uses raw `clientX - rect.left` /
`clientY - rect.top` directly without DPR scaling.

This is a known divergence from the older `_inhTooltip` pattern; the
fix to `_inhTooltip` is out of scope for this turn (separate bug, would
need its own change + test).

---

## 2. What this turn does NOT do

Same scope discipline as turn 161:

- **Run brackets.** Still deferred. Turn 161 already gated this on
  Quentin seeing real LG28 data first and asking for it.
- **Combinatorial enumeration button.** ~150 LOC and needs design
  thought; would not have fit alongside two other features.
- **Chain-break tick visualization.** Turn 161 listed this as a
  follow-up; still deferred.
- **Clipboard fallback for the TSV.** Browser `<a download>` is
  sufficient for now. If Quentin reports the download being blocked
  by his browser, we can add `navigator.clipboard.writeText(tsv)` as
  a try/catch fallback in 5 LOC.
- **CSV variant or per-run summary export.** TSV is one shape; if
  Quentin wants a runs-only sheet (one row per `_bandTraceRegimeRuns`
  entry), that's a separate ~30 LOC export function.
- **Band-picker dropdown next to the trace button.** Turn 161 §7
  flagged this as a small follow-up. Still queued.
- **`_inhTooltip` DPR fix.** Out of scope; latent bug noted.

---

## 3. Files touched

```
Inversion_atlas.html                                +401 LOC
  - state._btraceHits slot (written by drawer, cleared by clearCache)
  - _BTRACE_TOOLTIP_EL_ID constant
  - _btraceTooltipEnsureEl
  - _btraceTooltipBuildHtml(hit, opts)
  - _btraceTooltipShow(hit, clientX, clientY, opts)
  - _btraceTooltipHide
  - _btraceHitTest(hits, x, y)
  - _wireBandTraceTooltip(canvas)
  - _bandTraceToTSV(trace, opts)
  - _bandTraceDownloadTSV()
  - 9 window exports
  - linesBandTraceExportBtn ("📊 TSV") in lines header
  - export button click handler in the same wiring block as
    linesBandTraceFromCandBtn
  - _drawBandTraceStrip: hits stash init + push per L2 column
  - _bandTraceClearCache: also nulls _btraceHits
  - drawLinesPanel: invokes _wireBandTraceTooltip on PC1 sub-canvas

tests/test_turn162_band_trace_tooltip_export.js     new (115 assertions)
```

No existing functions modified semantically. `_drawBandTraceStrip`
gained a hit-recording side effect; `_bandTraceClearCache` clears one
extra slot. Both are additive.

---

## 4. Test results

**Single test**: 115 / 0 across 14 sections:

1. Source-pattern checks — function defs (8), exports (5), DOM
   ids (2), tooltip element id constant (1), wiring patterns (3),
   hits init/push (2), clearCache hits-clear (1) — 22 assertions
2-3. `_btraceTooltipBuildHtml` — co_seg / partial / no_valid /
   null-defensive paths, regime chip color, dominant swatch, entropy
   formatting, chain badge gating — 17 assertions
4. `_btraceHitTest` — hit / miss / edge-inclusive / null-safe — 8
5. `_wireBandTraceTooltip` — idempotency via dataset marker — 4
6-8. `_drawBandTraceStrip` hits — shape, ordering, height, off-screen
   filtering, toggle-off no-op, no-fish-set no-op, clearCache — 11
9-11. `_bandTraceToTSV` — header, K=3 columns, per-row fields, 6 dp
   floats, no_valid handling, empty per_l2 → header-only, null trace
   → null, missing envelopes → blank bp cols — 27 assertions
12. `_bandTraceDownloadTSV` — null when no fish-set, headless filename
   format (`band_trace_LG28_n5_K3.tsv`), null on empty envelopes — 7
13. DOM constants + tooltip element id — 3
14. Regression — turn 161/160/130/157A/122 surfaces still wired — 12

**Full sweep**: **3312 / 0** across all `test_turn*.js` files
(was 3197 / 0 at turn 161 close). Zero regressions.

JS-brace balance: clean (12,510 / 12,510). Largest script block parses
under `node --check` with no syntax errors. HTML parser shows 0 errors.

---

## 5. What Quentin should exercise

The two new surfaces:

**Tooltip flow:**
1. Load LG28, focus a candidate, click `🔍 trace` (turn 161 flow).
2. Hover the strip with the cursor.
3. A tooltip appears with: `L2 #<idx>` + a colored regime chip
   (green = co_seg, amber = partial, grey = fanned, dark-grey = sparse,
   very-dark = no_valid). Below: `n_valid: K / N`, optionally
   `dominant: <swatch> b<k> (xx.x%)`, `entropy: 0.213`. Below that, a
   row of small swatches showing each non-zero band's fraction.
4. When chains break (`n_chains > 1` in the cached trace), each
   tooltip header gets an additional `chain <i>` badge so Quentin can
   tell which Hungarian-chain segment the L2 belongs to.

**TSV flow:**
1. After running a trace, click `📊 TSV`.
2. Browser downloads `band_trace_<chrom>_n<n_fish>_K<K>.tsv`.
3. Open in Excel / R / awk:
   ```
   chrom  l2_idx  chain_idx  chain_position  start_bp  end_bp  n_valid  n_fish_selected  regime    dominant_band  dominant_fraction  entropy  band_fraction_0  band_fraction_1  band_fraction_2
   LG28   0       0          0               1234567   2345678 12       12               co_seg    0              0.833333           0.213000 0.833333         0.000000         0.166667
   ...
   ```
4. Filter by `regime == 'co_seg' AND dominant_band != <focal-band>`
   to get the long-range co-segregating L2s in one row each — that's
   the manuscript supplementary table.

**Specific patterns to verify:**

- Hover the focal candidate's own footprint: `regime == co_seg`,
  `dominant_fraction` should be very close to 1, entropy close to 0.
  If not, the trace is not seeing the focal candidate's fish at its
  own position — that would indicate a Hungarian-projection or label
  mismatch.
- Click `📊 TSV` immediately after `🔍 trace`. The export should
  include every L2 on the chromosome (count rows == count L2 envelopes).
  If short, the cache-key invalidation logic dropped some.
- Hover an L2 outside the chromosome's view range (zoom in first).
  The strip clips correctly; tooltip should not fire on rectangles
  that aren't painted.

---

## 6. What's NEXT (relative to the broader queue)

In priority order, mostly inherited from turn 161 §6:

1. **Run brackets for detected runs** (~30 LOC) — opt-in switch in the
   strip header. Draws thick top brackets over the L2 spans where
   `_bandTraceRegimeRuns` returned a run. Crosses into interpretation —
   only add after Quentin sees real data and asks for it.

2. **Combinatorial enumeration button** (~150 LOC) — auto-finder.
   Loops through all subsets of bands at a focal L2 (capped at 2^K - 2),
   produces ranked list of (subset, run) pairs by signal strength.
   Output as a small dialog. This is the "stop me from clicking 🔍 trace
   for every band" surface.

3. **Run the strip on real LG28** — calibrate the regime thresholds
   (`co_seg ≤ 0.40 / fanned ≥ 0.85`). Quentin's job, not ours; once
   he runs it, the calibration constants land in `_BTRACE_COSEG_ENTROPY_MAX`
   / `_BTRACE_FANNED_ENTROPY_MIN`.

4. **Band-picker dropdown** next to `🔍 trace` (~30 LOC) — small select
   that lets the user trace the focal candidate's band 0/1/2 without
   dropping into the console with `_bandTraceFromFocalCandidate({bandIdx: N})`.

5. **Per-run TSV export** (~30 LOC) — second download path that emits
   one row per `_bandTraceRegimeRuns` entry. Useful when Quentin wants
   the runs as a manuscript table without doing the per-L2 → run
   collapse in R.

6. **Chain-break tick visualization** (~20 LOC) — thin tick between
   chain neighbors in the strip when `chain_idx` differs.

7. **Clipboard fallback for TSV** (~10 LOC) — when Blob/anchor path
   fails, copy TSV text to clipboard with a toast.

After 1+2 land, `SPEC_distant_band_concordance_fish_trajectory.md` Slice 4
moves from `specs_todo/` to `specs_done/`. Items 4-7 are nice-to-haves.

There's also one **separate bug** worth flagging:

- **`_inhTooltip` DPR scaling** is latently wrong (CSS hits scaled by
  DPR ratio). On dpr=1 it works; on Retina/HiDPI it would miss-fire.
  Would need its own turn + test to fix without breaking the
  existing inheritance-pill flow.

---

## 7. Honest framing

**What's solid:**
- The tooltip is read-only over data the strip already computed —
  `state._btraceHits` is just a coordinate index into the existing
  `bandTraceCache`. No new compute, no cache-invalidation surface
  area.
- The TSV serializer is a pure function of the trace; it has no
  dependency on DOM or canvas. Easy to call from the console for
  debugging (`copy(_bandTraceToTSV(state.bandTraceCache, {chrom: 'LG28', envelopes: state.data.l2_envelopes}))`).
- The download path is headless-safe (returns filename without
  crashing in environments where Blob/document don't exist), and
  wrapped in try/catch so a blocked browser path doesn't take down
  the lines panel.
- Coordinate-space discipline: CSS pixels throughout, mirroring the
  newer `_vsTooltip` (turn 157A) pattern. Will work on Retina; the
  older `_inhTooltip` pattern would not.
- Click handler on the export button is independent of the toggle
  state — if the user has a fish-set but the strip is off, the TSV
  still exports (the export reads `bandTraceCache`, not
  `bandTraceOn`). This is intentional: the user might want the data
  without the visual.

**What's risky:**
- The tooltip could feel busy on a dense chromosome (28 L2s
  side-by-side at ~20 px each). If the L2 columns get smaller than
  ~6 px on a small viewport, the hit rectangles will overlap visual
  perception even though they don't overlap geometrically. No
  regression — the existing inheritance pills have the same property
  — but worth eyeballing on Quentin's actual screen size.
- TSV K column count is fixed per file at the trace's K. If Quentin
  runs traces at K=3 then K=6 and downloads both, the files will
  have different column counts — that's correct but worth noting in
  any awk/R pipeline.
- Filename uses `n<n_fish>` not the focal candidate id — if Quentin
  runs `🔍 trace` on the same band of the same candidate twice
  (with no other state change), both downloads will have the same
  filename and the second will overwrite the first in Downloads.
  Adding a timestamp suffix would fix this; deferred until requested
  because it adds noise to the typical filename.
- The `applyData` chrom-switch path was already clearing `bandTraceCache`
  via `_bandTraceClearCache`. With turn 162's edit, that same call
  now also clears `_btraceHits`. Verified by test, but if a future
  refactor breaks the `_bandTraceClearCache` call from `applyData`,
  hits could survive a chrom switch and fire false tooltips. The
  test catches this regression at the source-pattern level.

**What's queued:**
- Quentin runs LG28 → tooltip + TSV exercised on real data
- Threshold calibration based on hatchery noise
- Run-brackets only after Quentin asks
- Combinatorial enumeration as its own turn

End of handoff.
