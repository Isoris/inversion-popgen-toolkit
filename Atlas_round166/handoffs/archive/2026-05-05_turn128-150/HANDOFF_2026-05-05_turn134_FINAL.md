# HANDOFF — turn 134 — L2-sweep inspector modal + chrom-load hook (Slice 2)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (66,708 lines, 1064/0 tests passing)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: `HANDOFF_2026-05-05_turn133_FINAL.md` (L2-sweep Slice 1).

This turn ships **Slice 2** of `SPEC_l2_sweep_inheritance.md` — the
sweep-inspector modal — plus the **chrom-load hook** that turn 133's
§11 flagged as a Slice 1 gap.

---

## 0. Cohort discipline (NEVER conflate)

Same as turn 133. Three separate cohorts:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent
surname.

---

## 1. What this turn shipped

### 1.1 Chrom-load hook (Slice 1 follow-up)

Turn 133's handoff §11 flagged: *"no `applyData` / `onDataLoad` hook
yet that re-runs the sweep when the user switches chromosome with the
toggle already on."* That gap is now filled.

`applyData()` (line ~44823) now calls `invalidateL2SweepCache()` on
chrom switch, then if `state.l2SweepEnabled === true`:
1. Runs `runL2SweepInheritance({ force: true })` against the new chrom's
   L2 inventory
2. Calls `_autoPromoteFromSweep(result)` to add fresh auto candidates
3. Refreshes the inspector trigger-button availability

Order is correct: candidates from localStorage are loaded BEFORE the
sweep runs (so dedupe gate 4 sees them), and the hook fires AFTER
schema detection + layer-status refresh (so we don't sweep before
data is actually live).

### 1.2 Inspector modal (Slice 2 proper)

Per SPEC §4.2: "a modal that lists every L2 in the sweep with its
inheritance group_id, silhouette, band sizes, and a single-click
'promote' / 'dismiss' pair of buttons."

Plus per turn-133's recommended additions: silhouette + group-count
filters, status filter (all/pending/promoted/dismissed), bulk
promote-visible, bulk dismiss-visible, re-sweep button.

**Trigger**: `🔍 inspect` button in the L3 toolbar, sibling to the
L2-sweep checkbox. Disabled (greyed at 45% opacity) until
`state.l2SweepResult` is non-null.

**Modal**: `#l2SweepInspectorOverlay`, sibling to `#schemaRegistryOverlay`
and `#jsScriptsRegistryOverlay`. Same close pattern: ✕ button, click
outside, or Esc.

**Renderer**: `_renderL2SweepInspectorModal()` builds the full HTML
each call. Re-renders on every interaction (filter change, row
action, bulk action, re-sweep). Cheap given typically <100 rows.

**Per-row action: Promote**. The user explicitly clicked → bypasses
the 5 auto-promote gates. Builds the same auto-candidate shape as
`_autoPromoteFromSweep` (`source: 'auto_l2_sweep'`, `confirmed: false`,
inspector-bypass note in `notes`), pushes via `addCandidateToList`.

**Per-row action: Dismiss**. Adds l2idx to the per-chrom dismissed
set (same set Slice 1 honors). If the L2 currently has an
**auto-promoted** candidate (`source: 'auto_l2_sweep'` AND
`!confirmed`), also drops that candidate via `removeCandidateFully`.
**Manual or confirmed candidates are LEFT ALONE** — dismissing the
sweep view shouldn't nuke user-curated work.

**Bulk actions**: snapshot the current visible (filter-passing) row
set, then iterate. Promote skips already-promoted; dismiss skips
already-dismissed. Re-render after.

**Re-sweep**: calls `runL2SweepInheritance({ force: true })`, then
re-renders. Lets the user invalidate stale results without leaving
the modal (e.g., after K change, kMode switch, or chrom hot-reload).

---

## 2. Files changed

```
Inversion_atlas.html
  ├─ +block at ~line 36213 (after Slice 1 window-exports):
  │    L2-SWEEP INSPECTOR MODAL — Slice 2 (turn 134)
  │    ~580 LOC of:
  │      - 3 _L2SI_DEFAULT_* filter-default constants
  │      - _l2siFilterState() — lazy filter state on
  │        state.l2SweepInspectorFilter
  │      - _l2siRowMatchesFilter(meta, status, filter) — pure predicate
  │      - _l2siGroupsTouching(result, item_idx) — count groups per row
  │      - _l2siBuildRows() — derive renderable rows from sweep result
  │      - _l2siStatusFor / _l2siStatusChip / _l2siFmtBpRange — formatting
  │      - _l2siPromoteRow(l2idx) — gate-bypass per-row promote
  │      - _l2siDismissRow(l2idx) — dismiss + drop-auto-candidate
  │      - _renderL2SweepInspectorModal() — top-level renderer + interaction wiring
  │      - _refreshL2SweepInspectBtnAvailability() — gates trigger button on l2SweepResult
  │      - _wireL2SweepInspectBtn IIFE — DOM-ready binder
  │      - 9 window.X = X exports
  │
  ├─ +DOM trigger button in L3 toolbar (~line 6320):
  │    #l2SweepInspectBtn — "🔍 inspect", starts disabled, sibling to
  │    #l2SweepToggleLabel. Tooltip explains what the inspector does.
  │
  ├─ +DOM overlay div in header (~line 4696):
  │    #l2SweepInspectorOverlay — full-screen modal overlay, z-index
  │    1000, sibling to #schemaRegistryOverlay and
  │    #jsScriptsRegistryOverlay. Starts display:none.
  │
  ├─ +chrom-load hook in applyData (~line 44824):
  │    invalidateL2SweepCache() →
  │    runL2SweepInheritance({force:true}) + _autoPromoteFromSweep
  │      [if state.l2SweepEnabled]
  │    _refreshL2SweepInspectBtnAvailability()
  │
  └─ +button refresh in _l2SweepEnabledChange (~line 36185):
       call _refreshL2SweepInspectBtnAvailability after sweep+promote
       so toggle-on lights up the inspector button immediately.

  ALSO: defensive guard — _wireL2SweepInspectBtn now checks
  `if (btn.dataset && btn.dataset.l2sweepBound)` instead of bare
  property access. Caught during sandbox testing; a real DOM has
  `dataset` but stub elements may not.

tests/test_turn134_l2_sweep_inspector.js  [NEW]
  104 tests across 9 sections + 7 sandbox sub-scenarios:
    1. Source-level: function definitions + DOM elements + filter constants
    2. DOM elements (overlay div, trigger button)
    3. applyData chrom-load hook presence + ordering
    4. _l2siRowMatchesFilter behavior (NaN handling, all 3 filter axes)
    5. _l2siGroupsTouching (alternate "cells" vs "members" shape)
    6. _l2siStatusFor (dismissed-wins-over-promoted edge case)
    7. _l2siFmtBpRange formatting + null guards
    8. Full sandbox: build rows + promote + dismiss + filter state + button gating
    9. Modal renderer DOM output (empty + populated)
```

LOC delta: 66,092 → 66,708 = **+616** (Slice 2 ~580 + chrom-load
hook ~25 + button-refresh hooks ~10).

---

## 3. Tests

```
tests/test_turn134_l2_sweep_inspector.js — 104 / 0
```

Coverage map highlights:

| Sub-scenario | What it locks down |
|---|---|
| 4 — filter predicate | NaN silhouette is filtered out by minSil > 0; minGroups counts only populated bands; status='all' accepts everything |
| 5 — groupsTouching | Both `members` and `cells` member-array shapes; both `item_idx` and `itemIdx` keys; null result returns 0 |
| 6 — statusFor | Dismissed wins over promoted (most-recent-user-action semantics); null inputs don't crash |
| 7 — fmtBpRange | "1.00–2.00 Mb (1.00 Mb)" format; null inputs return "—" |
| 8a — buildRows | Status correctly derived from candidateList AND dismissed-set localStorage |
| 8b — promoteRow | Gate-bypass; idempotent (re-promote = no-op); candidate has correct shape |
| 8c — dismissRow | Adds to dismissed set; drops auto candidates; LEAVES manual candidates intact |
| 8d — filterState | Lazy-init on first call; mutation persists (object ref); defaults match constants |
| 8e — buttonGating | Disabled when no sweep; enabled when sweep populated; cursor reflects state |
| 9a — empty render | Empty-state message; close button; filter inputs; bulk buttons all present |
| 9b — populated render | One promote + one dismiss button per row; bp range formatted; status chips |

| | Tests | Files |
|---|---|---|
| Turn 132 baseline | 879 | 33 |
| Turn 133 Slice 1 | 960 | 34 |
| Turn 134 Slice 2 | 1064 | 35 |
| Δ this turn | +104 | +1 |

Zero regressions.

---

## 4. Sandbox vs. real-data note (still applies)

Same caveat as turn 133: tests run in `vm.createContext` with stubbed
DOM and stubbed action helpers. What the sandbox does NOT exercise:

- The real `addCandidateToList` → `persistCandidateList` →
  `_rebuildCandidateRegistries` chain. Slice 2 calls the same paths
  Slice 1 does, so if Slice 1's first real-data smoke test passes
  (auto candidates appear in BOTH `state.candidateList` AND
  `state.candidates`), Slice 2 will too.
- Real DOM event flow (focus, hover, real keyboard Esc handling).
  The renderer wires `addEventListener` for clicks + filter changes
  but the actual interaction loop is browser-side.
- Real CSS rendering — the modal uses `var(--good)`, `var(--rule)`,
  `var(--ink-dim)`, `var(--accent, #f5a524)` etc. matching the
  schema modal pattern; should look identical visually.

The browser smoke test for Slice 2: load any chrom, toggle L2-sweep
on, click 🔍 inspect, verify the modal opens with the expected row
count + headers. Click promote on a pending row → verify the chip
flips to "🤖 promoted" and the count in the header updates. Click
dismiss on a promoted row → verify both the chip flips to "🚫
dismissed" AND the candidate disappears from the catalogue.

---

## 5. State / window slots added this turn

```
state.l2SweepInspectorFilter   object — { minSil, minGroups, status }
                                lazy-initialized by _l2siFilterState
                                NOT persisted to localStorage (per-session)
```

```
window-exposed functions (added this turn):
_renderL2SweepInspectorModal             — top-level modal renderer
_refreshL2SweepInspectBtnAvailability    — gates trigger button visibility
_l2siBuildRows                           — pure row builder
_l2siRowMatchesFilter                    — pure filter predicate
_l2siGroupsTouching                      — count groups touching an item
_l2siStatusFor                           — promoted/dismissed/pending derivation
_l2siPromoteRow                          — explicit gate-bypass promote
_l2siDismissRow                          — dismiss + drop auto candidate
_l2siFmtBpRange                          — bp range to compact "X.XX Mb" string
_l2siFilterState                         — filter slot lazy initializer
```

```
DOM ids claimed (this turn):
#l2SweepInspectBtn          — trigger button in L3 toolbar
#l2SweepInspectorOverlay    — modal overlay (top-level header sibling)
#l2SweepInspectClose        — close button in modal header
#l2siFilterMinSil           — filter input
#l2siFilterMinGrps          — filter input
#l2siFilterStatus           — filter dropdown
#l2siFilterReset            — reset filter button
#l2siBulkPromote            — bulk promote-visible button
#l2siBulkDismiss            — bulk dismiss-visible button
#l2siResweep                — re-sweep button
._l2siPromoteBtn (class)    — per-row promote button
._l2siDismissBtn (class)    — per-row dismiss button
```

```
localStorage keys claimed:
(none new this turn — Slice 2 reuses Slice 1's
 pca_scrubber_v3.l2SweepDismissed.<chrom> via _addL2SweepDismissed)
```

---

## 6. Backups present

```
Inversion_atlas.html.bak_pre_l2_sweep_slice1     (turn 133 baseline)
Inversion_atlas.html.bak_post_l2_sweep_slice1    (after turn 133)
Inversion_atlas.html.bak_post_l2_sweep_inspector (this turn — current state)
```

`.bak_*` files NOT in bundle.

---

## 7. What this is NOT

- **Not a candidate-promotion sandbox**. Promote in the inspector
  goes straight to `state.candidateList` via the production
  `addCandidateToList` path. There's no "preview" intermediate
  state. If the user clicks promote and immediately clicks dismiss,
  the candidate is created and then removed; both writes go through
  the canonical paths. (Acceptable: it's how the rest of the atlas
  works, and the dismiss path is fast.)
- **Not lasso/multi-select**. Per-row + bulk-visible only. If the
  user wants to promote a specific subset, they use the filters to
  narrow the visible set, then bulk-promote. Lasso would need a row
  click-to-toggle UI; not in this spec.
- **Not column-sortable**. Rows render in `start_bp` order (the same
  order `runL2SweepInheritance` produces). For LG28 with ~30 L2s
  that's fine; if Quentin wants sort-by-silhouette-descending, that's
  a Slice 3 ask.
- **Not row-click navigation**. Clicking a row doesn't navigate the
  scrubber to that L2. Could be added — `state.cur = env._s0` plus
  a `requestRepaint` would do it. Deliberately deferred: the
  inspector is for triage, not for jumping into the L2; the user
  closes the modal first if they want to dive in.
- **Not "auto candidate visualization"**. The catalogue already shows
  auto candidates (turn 130 plumbing); this modal is a separate
  surface for sweep-result triage. They show the same data through
  different lenses.

---

## 8. Where to start the next chat

Same recommendation as turn 133 §10b — **calibrate Slice 1 + 2 on
real data**.

The L2-sweep + inspector is now feature-complete enough that
Quentin can:
1. Toggle L2-sweep on, load LG28
2. Auto-promotes appear (if any pass the gates)
3. Click 🔍 inspect to see the full sweep result
4. Filter + promote / dismiss to triage
5. Verify the LG28 ~2.89 Mb inversion is auto-promoted as ONE
   candidate (not multiple overlapping L2s; gate 5 should collapse
   them)
6. Repeat on LG12

If gates need tuning:
- `_AUTO_PROMOTE_MIN_SILHOUETTE` (currently 0.30)
- `_AUTO_PROMOTE_DEDUPE_BP` (currently 100kb)
- `_AUTO_PROMOTE_MIN_BAND_SIZE` (currently 5)

These are top-of-file constants in the L2-SWEEP block; trivial to
adjust.

### Other options:

**Option 8a — Trajectory matrix viewer Slice 3** (turn 132 §10b
queue item #3). Requires lineage compute to be working; that landed
in turn 130.

**Option 8b — Cross-chromosome lineages Slice 1** (queue item #4).

**Option 8c — Page-12 Slice 8** (state machine + remaining renderers
from turn 132). Tracked-samples + L3 contingency for page 12; needs
state.thCur / thTracked / thAnchorWin / thViewMode.

**Option 8d — G-panel scaffold Slice 1** (Tier 2 review surface).
Was queued behind L2-sweep in turn 132.

---

## 9. Honest framing

**What turn 134 actually delivered:**

- Closed the chrom-load hook gap that Slice 1 left open (now
  switching chrom with toggle on re-runs the sweep against the new
  chrom's L2 inventory, no re-toggling needed).
- A working sweep inspector with per-row + bulk promote/dismiss,
  silhouette + group-count + status filters, re-sweep, defensive
  status derivation that survives the `dismissed wins over promoted`
  edge case.
- 104 new tests covering every action path, filter axis, status
  derivation, formatting helper, and DOM render shape.

**What it deliberately didn't deliver:**

- Column sort. Single sort key (start_bp) for now.
- Row-click scrubber navigation. Triage doesn't need it; if Quentin
  asks, Slice 3.
- Persistent filter state. Resets on reload — reasonable since the
  user's filter intent for one sweep doesn't generally apply to the
  next chrom.
- A "🤖 swept" badge on the lines header. Mentioned in SPEC §4.1 but
  is pure styling on existing chrome; would be a one-liner if Quentin
  asks.

**What this means for the manuscript:**

Same as turn 133: the L2-sweep gives a near-comprehensive auto-recovery
claim for the methods section. The inspector is operationally valuable
for Quentin during the curation pass but doesn't add anything to the
manuscript narrative directly — except possibly a supplementary figure
showing the inspector's filter UI as the "review surface" by which the
candidate atlas was curated.

---

## 10. Bundle contents

When bundled for next turn:

```
Inversion_atlas.html              (current, 66,708 lines)
tests/                            (all *.js)
  ├─ test_turn134_l2_sweep_inspector.js     [NEW, 104/0]
  ├─ test_turn133_l2_sweep_auto_promote.js  (unchanged, 81/0)
  └─ test_turn132_*.js                      (unchanged, 321/0)
specs_todo/                       (active build queue)
specs_new_turn131/                (pending review queue, unchanged)
HANDOFF_2026-05-05_turn134_FINAL.md  (this file)
all previous handoffs             (kept for history)
OBSERVATIONS_TO_FIX.txt
```

Walk the map carefully, respect cohort discipline, don't break the
test suite. L2-sweep Slice 1+2 are in a clean stopping state — pick up
with calibration on real LG28/LG12 data, or pivot to one of the other
priority items.
