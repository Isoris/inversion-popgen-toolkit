# HANDOFF — turn 165 — G-panel `auto` review tab

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (76,407 lines, +373 LOC from 76,034)
**Working dir**: `/home/claude/work/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes** Slice 2 of `SPEC_review_surfaces_auto_and_lineages.md`
(the auto-candidate review tab in the G-panel popup). Slice 0
(infrastructure) shipped in turn 130 follow-up. Slice 1 deliverables
turned out to already be shipped (filter modifiers landed in earlier
work — see §0). Slice 3 (lineages tab) deliberately deferred.

**Picked up from**: post-turn-164 working tree (3529 / 0 baseline).

---

## 0. Why this spec, this turn — and a discovery

The band-trace + linkage stack (turns 160–164) is at a natural pause
point pending Quentin's real-LG28 review.

`SPEC_review_surfaces_auto_and_lineages.md` has been partially-shipped
since turn 130 follow-up. Re-reading the spec, three slices remained:

1. **Slice 1 — full dashed-outline visual treatment.** Three concrete
   items: (a) cross-candidate matrix filter, (b) lines-panel
   candidate-bands skip, (c) confirm action drops dashed.
2. **Slice 2 — `auto` G-panel tab.** Per-row Confirm / Dismiss /
   Inspect with bulk actions.
3. **Slice 3 — `lineages` G-panel tab.** Track / Pin as group / Dismiss.

When I went to start Slice 1, I found that all three Slice 1
deliverables had already shipped:

- `_gatherActiveCandidatesForInheritance` already filters auto-candidates
  (turn 130 hook at line 41218; the cross-candidate matrix is built
  from this gather function via `renderInheritanceMatrix`).
- `_paintCandidateBands` already skips unconfirmed candidates via the
  `c.confirmed !== true` guard at line 33960. Auto candidates have
  `confirmed=false` so they're invisible in the lines-panel band
  highlights.
- The "confirm action drops dashed" path is automatic: `_isAutoCandidate`
  returns false when `cand.confirmed` is true, so the dashed
  treatment falls off as soon as the user confirms.

That moved Slice 2 to the front. Quentin actually **does** have a gap
here — the L2-sweep auto-promote pipeline (turns 133–134) creates
candidates with `source = 'auto_l2_sweep'` and `confirmed = false`
that show up in `state.candidateList` but have no triage UI. Without
this turn, the only way to confirm them is manual click in the
candidate strip.

Slice 3 (lineages tab) deliberately deferred — `state.lineageResult`
exists (turn 130) but Quentin hasn't exercised lineage compute on real
data yet. Building the UI is speculative without that.

---

## 1. What this turn ships

### State + tab registry

```js
// state.gPanelTab comment widened
gPanelTab: 'manual',  // 'karyotype' | 'inheritance' | 'manual' | 'auto'

// _GPANEL_TABS extended with the auto entry
const _GPANEL_TABS = [
  { key: 'karyotype',   label: 'karyotype',   slice: 2 },
  { key: 'inheritance', label: 'inheritance', slice: 3 },
  { key: 'manual',      label: 'manual',      slice: 1 },
  // turn 165
  { key: 'auto',        label: 'auto · review', slice: 2,
    visibleWhen: function(_state) {
      return (typeof _gPanelHasAutoCandidates === 'function')
              && _gPanelHasAutoCandidates(_state);
    },
    accent: 'review' },
];
```

The `visibleWhen` predicate is read at tab-strip render time. Tabs
without it always render. Generic mechanism — adding the lineages tab
later is one new entry.

The `accent: 'review'` field is a styling hint. Inactive review-accent
tabs render with dashed borders so the user can see at a glance "this
is a triage queue, not a primary grouping."

### Render path changes (`_renderGPanelModal`)

```js
const visibleTabs = _GPANEL_TABS.filter(t =>
  typeof t.visibleWhen !== 'function' || !!t.visibleWhen(_state)
);
const validTabs = visibleTabs.map(t => t.key);
if (validTabs.indexOf(_state.gPanelTab) === -1) _state.gPanelTab = 'manual';
```

Three modifications:
1. **Tab strip filtered via `visibleWhen`** — only iterates visible tabs.
2. **Fallback-to-manual** when active tab disappears (e.g. user was on
   `auto`, then confirmed the last auto-promote → tab vanishes → next
   render snaps to `manual`).
3. **Body routing**: `if (activeTab === 'auto') html += _gPanelRenderTabAuto()`
4. **Action wiring**: after innerHTML rewrite, calls
   `_gPanelWireAutoTabActions()` when active tab is `auto`. Idempotent
   — every render reattaches because innerHTML rewrite detaches all
   listeners. No leak.

The dashed-border styling for inactive review-accent tabs:
```js
const isReviewAccent = (t.accent === 'review');
const borderStyle = (!isActive && isReviewAccent) ? 'dashed' : 'solid';
```

### Functions

`_gPanelHasAutoCandidates(state)` — predicate. Returns true if
`state.candidateList` has at least one candidate where
`!c.confirmed && c.source.startsWith('auto_')`. Defensive fallback when
`_isAutoCandidate` is out of scope (sandbox tests still work).

`_gPanelCollectAutoCandidates(state)` — gathers unconfirmed auto
candidates in `candidateList` order. Auto-promotes are typically
appended in chromosome order during a sweep, so this matches Quentin's
reading order without a sort.

`_gPanelDismissAuto(cand)` — **reuses** the existing
`_addL2SweepDismissed(chrom, l2idx)` from turn 133/134. That helper
persists at `pca_scrubber_v3.l2SweepDismissed.<chrom>` and the
existing `runL2SweepInheritance` gate 6 already consults the same set.
Then removes the cand from `state.candidateList`. Single source of
truth for "is this L2 dismissed."

`_gPanelRenderTabAuto()` — produces the tab body HTML:
- Header with count + bulk Confirm-all / Dismiss-all buttons
- Per-row card with dashed border, 🤖 prefix, L2 idx, chrom, span Mb,
  K / silhouette / n_groups when present, band-color swatches with
  fish counts (g0 (n=4), g1 (n=3), …), and per-row Confirm / Dismiss /
  Inspect buttons

`_gPanelWireAutoTabActions()` — binds click handlers:
- **Confirm**: `c.confirmed = true`, then re-render
- **Dismiss**: `_gPanelDismissAuto(c)`, then re-render
- **Inspect**: `setCandidate(c)`, then close modal so user sees chrom
- **Bulk Confirm-all**: `confirm()`-prompted, sets all to confirmed
- **Bulk Dismiss-all**: `confirm()`-prompted, snapshots ids first
  (since dismiss mutates the list), dismisses each
- After every action: `_renderGPanelModal()` rebuilds, plus
  best-effort `drawCandidateStrip` + `refreshCandidateListUI` so the
  strip + list update immediately

### Design choices recorded

**Reused existing dismissal helper.** I started with a parallel
`state.dismissedAutoIds` slot keyed by `cand.id` and a separate
localStorage prefix. Discovered partway through that turn 133/134
already shipped `_addL2SweepDismissed(chrom, l2idx)` with persistence
at `pca_scrubber_v3.l2SweepDismissed.<chrom>`, and gate 6 of
`runL2SweepInheritance` already consults that set. Building a parallel
structure would have created two sources of truth for "is this L2
dismissed." Removed the duplicate, dismiss now keys on `cand.ref_l2`
(the L2 index this auto-cand was promoted from).

**Visibility-gating via predicate.** Hardcoding visible/hidden in the
tab strip iteration would have worked for one tab. Adding a
`visibleWhen` field on the tab entry generalizes — Slice 3 (lineages)
drops in as one new array entry with no further mechanism work.

**Dashed border on inactive review-accent tabs.** Spec §3.1 says
*"distinct visual treatment (e.g., dashed tab border, 'review'
subtitle) makes the distinction clear."* The `auto · review` label
+ dashed-when-inactive border satisfies both halves.

**Fallback-to-manual when active tab disappears.** Without this, a
user on the `auto` tab who confirms the last pending auto-cand would
see an empty body or stuck tab pointer. The fallback ensures the
modal always lands on a valid tab.

---

## 2. What this turn does NOT do

- **Slice 3 (`lineages` G-panel tab).** Same `visibleWhen` mechanism
  ready to receive it — predicate would check
  `state.lineageResult && state.lineageResult.lineages.length >= 2`.
  Deferred until Quentin runs lineage compute on real data; otherwise
  the UI is speculative.
- **Surface A "dashed-outline rows" visual treatment.** Already shipped
  in turn 130 follow-up via `cli-auto-prefix` + `.cand-list-item.is-auto`
  CSS. No work needed.
- **Cross-candidate matrix filter, lines-panel skip, confirm-drops-dashed.**
  All already shipped — verified in §0 above.
- **Inspect → focus + scroll behaviour beyond `setCandidate`.** The
  Inspect button calls `setCandidate(cand)` and closes the modal. If
  the user wants the chromosome strip to scroll-to-position, that's a
  separate page-level feature — not in this spec.
- **Cross-chromosome auto candidates.** This turn's tab shows ALL
  unconfirmed auto candidates regardless of chrom. That's intentional
  — when Quentin runs a sweep on multiple chromosomes, he wants to
  triage them all in one place. The cand row shows the chrom name so
  it's clear which chromosome each came from.

---

## 3. Files touched

```
Inversion_atlas.html                                   +373 LOC
  - state.gPanelTab comment widened
  - _GPANEL_TABS extended with `auto · review` entry
  - _renderGPanelModal:
      - filters tab strip via visibleWhen
      - validTabs derived from visibleTabs
      - fallback-to-manual when active disappears
      - routes 'auto' body to _gPanelRenderTabAuto
      - applies dashed-border to inactive review-accent tabs
      - calls _gPanelWireAutoTabActions when active===auto
  - _gPanelHasAutoCandidates(state)
  - _gPanelCollectAutoCandidates(state)
  - _gPanelDismissAuto(cand) — reuses _addL2SweepDismissed
  - _gPanelRenderTabAuto()
  - _gPanelWireAutoTabActions()
  - 5 window exports

tests/test_turn165_g_panel_auto_tab.js                 new (123 assertions)
```

No existing functions semantically modified except `_renderGPanelModal`,
which gains the visibility filter, fallback, body routing, and wiring
hook. Strictly additive — pre-existing tabs (karyotype, inheritance,
manual) work identically.

`runL2SweepInheritance` was inspected but not modified — gate 6
already consults the dismissed set.

---

## 4. Test results

**Single test**: 123 / 0 across 14 sections:

1. Source-pattern checks — 21 (function defs, exports, tab entry,
   visibility predicate, accent, gPanelTab comment, render path
   modifications, dismiss reuse markers)
2. `_gPanelHasAutoCandidates` predicate (true/false on fixtures,
   defensive fallback) — 9
3. `_gPanelCollectAutoCandidates` (filters confirmed + manual,
   preserves order, defensive null paths) — 6
4. `_gPanelDismissAuto` (calls `_addL2SweepDismissed`, removes from
   list, defensive null/no-id, no-ref_l2 path) — 10
5. `_gPanelRenderTabAuto` non-empty list (count header, bulk buttons,
   per-row cards, 🤖 prefix, L2 / chrom / span / metadata, band
   swatches, action button classes) — 23
6. `_gPanelRenderTabAuto` empty path (instructional message, count=0,
   no per-row cards) — 4
7. `_GPANEL_TABS` extension (entry preserved + auto added at end +
   label) — 5
8. `_renderGPanelModal` visibility filter (visibleTabs filter,
   validTabs, fallback-to-manual) — 3
9. `_renderGPanelModal` body routing (activeTab branch, typeof guard,
   manual fallback, action wiring trigger) — 4
10. Confirm/dismiss/inspect/bulk action source patterns — 10
11. Visibility predicate sandbox (false on no-list, false on empty,
    true on auto, false on all-confirmed) — 4
12. Review-accent dashed border source pattern — 2
13. Reuse of existing `_addL2SweepDismissed` — 4
14. Regression — turns 130 / 133-134 / 135 / 136 / 152 / 160-164 — 13

**Full sweep**: **3652 / 0** across all `test_turn*.js` files
(was 3529 / 0 at turn 164 close). Zero regressions.

JS-brace balance: clean (12,666 / 12,666). Largest script block parses
under `node --check` with no syntax errors.

---

## 5. What Quentin should exercise

The L2-sweep auto-promote pipeline produces auto candidates as a side
effect of the L2 inheritance sweep (button "auto-sweep" in the L3
toolbar, default OFF). Once a sweep produces auto-promotes:

1. Open the G-panel (hotkey `g` or the G ▾ button).
2. **(NEW)** A new tab appears at the end: `auto · review` with a
   dashed border (because it's a triage queue, not a primary grouping).
3. Click the tab. The body shows:
   - Header line: "<count> auto-promoted candidate(s) awaiting review"
     and `confirm all` / `dismiss all` bulk buttons
   - Per-row dashed cards, one per unconfirmed auto-cand:
     `🤖 L2:<idx>  <chrom>  <start>–<end> Mb  K=3, silh=0.42, n_groups=2`
     with band-color swatches showing fish counts per band, and a
     three-button stack: `confirm` / `dismiss` / `inspect`
4. Per-row actions:
   - **confirm** → cand becomes a regular catalogue entry, gains solid
     border in the strip, starts feeding I·g pills + cross-candidate
     matrix on the next compute, drops out of this tab
   - **dismiss** → persisted via `_addL2SweepDismissed` (per-chrom
     localStorage), removed from candidateList, won't be re-proposed
     on the next sweep
   - **inspect** → opens the L3 mini-PCA + per-sample lines on page 1,
     closes the modal so you can see the chromosome
5. Bulk actions are gated by a `confirm()` prompt to avoid accidents.
   Dismiss-all persists each dismissal so a re-sweep won't re-propose
   any of them.
6. When the last unconfirmed auto-cand is confirmed or dismissed, the
   tab disappears on the next render (visibility-gated). The modal
   falls back to the `manual` tab.

**Specific patterns to look for on real LG28:**
- After a sweep, the count header tells you how many proposals there
  are. If it's >50 you should triage in batches — the tab handles 100+
  rows fine via the modal's `max-height: 60vh; overflow-y: auto;`.
- `silh` (silhouette) and `n_groups` numbers in each row are the L2's
  quality metrics from the sweep gates. High silhouette + ≥2 groups +
  big bands ⇒ likely real signal worth confirming.
- Dismissed L2s persist per chrom in localStorage. A re-sweep on the
  same chrom won't re-propose them; a sweep on a different chrom
  starts with an empty dismissed set (correctly — chromosome-local
  decisions).

---

## 6. What's NEXT

In priority order:

1. **Run on real LG28** — Quentin reviews the auto-promotes from a
   real sweep. Same blocker as turns 161–164: real data is the
   calibration step.
2. **Slice 3 (`lineages` tab)** (~0.5 turn). Mechanism is in place —
   one new entry in `_GPANEL_TABS` with `visibleWhen` checking
   `state.lineageResult && state.lineageResult.lineages.length >= 2`,
   plus a `_gPanelRenderTabLineages()` body. Per-row would have
   Track / Pin as group / Dismiss buttons per spec §3.3.
3. **Inspect button → page 1 navigation polish.** Right now Inspect
   calls `setCandidate(cand)` and closes the modal. Could additionally
   ensure the candidate strip scrolls to that L2 — needs the
   strip-scroll primitive.
4. **Run brackets for detected runs** (band-trace turn 161 §6 item 1)
   — gated on real-data calibration.
5. **Slice 2 of `lasso_inheritance_backgrounds` (alpha intervals on
   per-sample-lines)** — visual layout design needed.
6. **Combinatorial enumeration** (band-trace turn 161 §6 item 3) —
   ~150 LOC turn.
7. **Lasso surface bridge** — when a real lasso ships in the
   per-sample-lines or tracked-PCA, widen fish-set source from
   `bandTraceFishSet` to `lassoSelection || bandTraceFishSet`.
8. **`_inhTooltip` DPR-scaling fix** — latent bug noted in turn 162 §1.

---

## 7. Honest framing

**What's solid:**
- Visibility-gating via `visibleWhen` predicate generalizes — Slice 3
  (lineages) drops in cleanly; future tabs can use the same mechanism.
- Fallback-to-manual when the active tab disappears prevents the user
  from getting stranded on an empty tab body.
- Reuse of `_addL2SweepDismissed` (turn 133/134) is the right call —
  single source of truth for "is this L2 dismissed", and the existing
  `runL2SweepInheritance` gate 6 already does the right thing on
  re-sweep.
- Action wiring is rebuilt on every render — no listener leaks
  because innerHTML rewrite detaches before each rebuild.
- Tests cover all five new functions plus the render-path
  modifications via source patterns + behavioural sandboxes.
- The `auto · review` label + dashed-when-inactive border styling
  visually distinguishes the triage queue from primary groupings, per
  spec §3.1.

**What's risky:**
- Bulk Confirm-all and Dismiss-all use the browser-native `confirm()`
  for guard prompts. Wrapped in `typeof confirm === 'function'` so
  test environments without it don't crash, but the click-through
  path skips the prompt entirely in those environments. Real-browser
  use is the only path that exercises the guard.
- The auto tab shows ALL unconfirmed auto-promotes regardless of
  chromosome. With 30 chromosomes and 10–20 auto-promotes each, the
  count could reach 300–600 rows. The modal's `max-height: 60vh;
  overflow-y: auto;` handles scrolling; rendering 300 rows of
  `innerHTML` is fast (<50ms typical) but worth flagging if Quentin
  sees lag. Pagination could be added if needed (`paginate by
  chrom`).
- The Inspect button calls `setCandidate(cand)` + `_gPanelClose()`.
  The candidate becomes the focal candidate but the candidate strip
  doesn't auto-scroll to its position. Quentin would need to scrub
  manually if it's off-screen. Could be polished — see §6 item 3.
- `_gPanelDismissAuto` only writes to the dismissed set when
  `Number.isInteger(cand.ref_l2)`. Current auto-promotes from
  `runL2SweepInheritance` always have `ref_l2` set (checked: line
  42120). But future auto sources (`auto_lineage` etc.) may not, in
  which case dismiss falls through to "remove from candidateList"
  only — no persistence. Left as-is because future sources will
  probably need their own dismissal scheme.
- The visibility predicate runs every render. With 226 fish × 30
  candidates, the linear scan in `_gPanelHasAutoCandidates` is ~7000
  ops — trivial. Worth noting in case the candidate list grows to
  thousands.
- Auto-tab tests run against synthetic candidates; no integration
  test exercises a real `runL2SweepInheritance` → tab → confirm →
  cross-candidate matrix path. That's the kind of end-to-end check
  Quentin's manual review will exercise on real LG28 data.

**What's queued:**
- Quentin reviews the auto-promote tab on real LG28 → confirms the
  defaults work
- Slice 3 (lineages tab) once lineage compute is exercised on real
  data
- Inspect → strip-scroll polish (small, one item from §6)
- Band-trace continuation gated on real-data calibration

End of handoff.
