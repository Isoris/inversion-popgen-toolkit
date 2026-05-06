# HANDOFF — turn 147 + 147b — Atlas UI fixes #1 + #2

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (71,620 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 146 Turn-D handoff.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok).

---

## 1. What this turn ships

Two of the four screenshot-fixes Quentin flagged in turn 144 are now done.

### 1.1 Fix #1 — L3 toolbar consolidation (turn 147)

Quentin: *"L3 contingency-table toolbar sprawls across 3 rows on page 1,
~120px wasted vertical."*

The page-1 L3 toolbar carried 21 controls across 2 rows by design (a
manual `flex-basis: 100%` line break split rows 1 and 2), but on
narrower viewports row 1 wrapped and produced 3 visual rows.

**Fix**: a `⋯ more` disclosure pattern that hides 14 secondary controls
behind one click while keeping the workhorse controls always-visible.

- 14 secondary controls each carry `class="l3-more-item"`:
  color mode, K-mode, het toggle, L2-sweep toggle, inspect, G-panel,
  band-continuity scope, scale-mode, scale panes, detailed,
  show tracked, clear spot, pin 2nd, live dosage
- 7 workhorse controls stay visible:
  layout group (Focal/+L/+R/+L/R/±2), compare unit (L2/1w/5w/10w/Nw + N
  input), recluster + draft action row + edit row, ★ promote, ⋯ more
  toggle, ▼ collapse
- New CSS: `#l3Panel.l3-more-collapsed .l3-more-item { display: none !important; }`
  (the `!important` beats inline-style display set by per-control logic
  like `#l3ScalePanes` flipping with scale-mode toggle)
- New button: `#l3MoreToggleBtn` — shows `⋯ more` (collapsed default)
  or `⋯ less` (expanded). Accent-colored when expanded (matches `.active`
  convention)
- State persists: `localStorage 'pca_scrubber_v3.l3MoreExpanded'`
- The redundant `<div class="l3-head-break">` was removed — with the
  secondary controls now collapsed by default, the forced row split was
  creating wasteful empty space; `flex-wrap: wrap` on `.l3-head` handles
  narrow viewports without it
- Every hidden control keeps its existing event listeners — only
  visibility changes. No functionality lost.

### 1.2 Fix #2 — Windows-(N) button engagement state (turn 147b)

Quentin (screenshot 2): the "📊 Windows (N)" button at the top-right
showed no visual indication of engagement when N≠1 (i.e., when the user
had selected a custom step size from the sidebar `winN` mode).

Pre-fix:
- `state.stepMode='l2'` → label "📊 Windows (1)", dim. **Misleading** —
  looks like step=1 is active, but arrows actually jump L2 envelopes.
- `state.stepMode='winN'` → label "📊 Windows (1)", dim. **Wrong** —
  the cycler isn't the live setting, but the user IS stepping in
  windows; the button hides this completely.

Post-fix (turn 147b):
- `'l2'` (or `null`) → label "📊 Windows (—)", dim. Em-dash signals
  "I'm not the active step setting." Tooltip explains arrows jump L2.
- `'win1' / 'win5' / 'win10'` → label shows the actual N, **lit**
  (accent border + accent text).
- `'winN'` → label shows the custom N from `state.stepModeN`, **lit**.
  Tooltip explains "set to N windows per step (custom value from sidebar)."
- Defensive: invalid `stepModeN` (e.g. negative) falls back to N=1.
- The sidebar `#stepModeNInput` now refreshes the cycler label
  immediately when the user types a new N value (was silent before).

Plus: turn 128's test was updated to reflect the new (correct) behavior.
Two of its assertions were *codifying the bug Quentin flagged* — those
have been inverted to assert the fix.

---

## 2. Files changed / added

| File | Change |
|---|---|
| `Inversion_atlas.html` | +75 LOC for turn 147 (CSS rule + 14 class additions + toggle button + JS wiring) + ~30 LOC for turn 147b (rewritten `_refreshStepSizeBtn` + extra refresh hook in stepModeNInput handler). Net +108 LOC. |
| `tests/test_turn147_l3_toolbar_disclosure.js` | NEW — 66 / 0 |
| `tests/test_turn147b_windows_n_engagement.js` | NEW — 26 / 0 |
| `tests/test_turn128_window_cycler_lit.js` | UPDATED — 27 / 0 (was 22 / 4 after turn 147b changed cycler behavior; 2 assertions inverted to test the fix, 2 brittle source-level patterns relaxed) |
| `tests/test_turn134_l2_sweep_inspector.js` | UPDATED — 104 / 0 (was 103 / 1 after turn 147 added a class attribute that pushed text past a 600-char lookahead; lookahead bumped to 800) |

---

## 3. Tests

### 3.1 New + updated suites

- `tests/test_turn147_l3_toolbar_disclosure.js` — **66 / 0** across 8 sections:
  CSS rules (4), HTML markers (4), 14 secondary controls each marked (14),
  5 workhorses NOT marked (5), line-break removed (2), JS wiring patterns (7),
  sandboxed click round-trip with mocked localStorage (10),
  functional preservation (20).
- `tests/test_turn147b_windows_n_engagement.js` — **26 / 0**:
  source-level pattern checks (7), sandboxed `_refreshStepSizeBtn` against
  every stepMode value with edge cases (16), CSS rule still wired (2),
  + 1 totals.
- `tests/test_turn128_window_cycler_lit.js` — **27 / 0** (was 22 / 4):
  the 2 inverted assertions now test the corrected behavior (cycler LIT
  for winN, label "(—)" for l2); the 2 brittle source-level patterns now
  use larger lookahead windows.
- `tests/test_turn134_l2_sweep_inspector.js` — **104 / 0** (was 103 / 1):
  a `class="l3-more-item"` attribute pushed the inspector button's "🔍 inspect"
  text past the test's 600-char lookahead; bumped to 800.

### 3.2 Adjacent suites unchanged

- turn 146 breeding-card bulk: **136 / 0**
- turn 145 server-unify: **78 / 0**
- turn 144 breeding-card render: **162 / 0**
- turn 143 breeding-card compute: **196 / 0**
- turn 142 cohort_diversity loader: **133 / 0**
- turn 129 L3 het coloring: **58 / 0** (touched `#l3HetToggleLabel`, no regression)
- turn 133 L2 sweep auto-promote: **81 / 0** (touched `#l2SweepToggleLabel`, no regression)
- turn 135 G-panel slice 1: **73 / 0** (touched `#gPanelOpenBtn`, no regression)

### 3.3 Full sweep

**Across parseable turn-numbered tests: 2262 / 0** (was 2169 at end of
turn 146, +93 net: +66 turn 147 + 26 turn 147b + 1 from the updated
turn 128 test). Same 10 broken-on-baseline (missing fixtures, unrelated).

---

## 4. Atlas state

|                          | LOC     | Tests | Files |
|---                       |---      |---    |---    |
| Pre-turn (turn 146)      | 71,512  | 2169  | 49    |
| Post-turn (this)         | 71,620  | 2262  | 51    |
| Δ                        | +108    | +93   | +2    |

---

## 5. Architecture (turn 147)

### Disclosure pattern (turn 147)

The L3 toolbar uses a single state class on the panel root:

```html
<div class="panel l3-more-collapsed" id="l3Panel">    <!-- ships collapsed -->
  <div class="l3-head">
    <!-- always-visible: title, layout group, compare unit, ... -->
    <div class="l3-layout l3-more-item" id="l3ColorMode">     <!-- hidden -->
    <div class="l3-layout l3-more-item" id="l3KMode">         <!-- hidden -->
    ...14 total secondary controls...
    <button id="l3MoreToggleBtn">⋯ more</button>              <!-- toggle -->
    <button id="l3CollapseBtn">▼</button>                     <!-- panel collapse -->
  </div>
  <div class="l3-body">...</div>
</div>
```

```css
#l3Panel.l3-more-collapsed .l3-more-item {
  display: none !important;       /* beats per-control inline style */
}
```

The `!important` is essential because some hidden controls have their own
visibility logic (e.g. `#l3ScalePanes` flips display on scale-mode toggle);
when our disclosure says "hide", we need to win.

### Step-cycler (turn 147b)

The button reflects ALL five `state.stepMode` values:

| stepMode | label             | lit?  |
|---       |---                |---    |
| `l2`     | `📊 Windows (—)`  | dim   |
| `null`   | `📊 Windows (—)`  | dim   |
| `win1`   | `📊 Windows (1)`  | LIT   |
| `win5`   | `📊 Windows (5)`  | LIT   |
| `win10`  | `📊 Windows (10)` | LIT   |
| `winN`   | `📊 Windows (N)` (where N=`state.stepModeN`) | LIT |

Tooltip varies per state and explains what arrows currently do.

---

## 6. What's NOT done — still queued

Two screenshot-fixes from turn 144 remain:

3. **1w/Nw L3 panel parity with L2 mode** — `renderL3PanelSlab`
   (line ~43230) needs the same rich pane structure as `renderL3Panel`
   (line ~42706). Currently the slab-mode (1w / 5w / 10w / Nw)
   renderer is a downgraded variant of the L2-mode renderer.
   Estimated ~150–250 LOC.

4. **G-popup → popgen page merge** — fold `_gPanelToggle` modal overlay
   into the popgen page as a second tab strip with different background
   shade. Per Quentin: "merge into one page with 2 rows of tabs
   (different bg shade)". Estimated ~400 LOC.

Plus from the docs track:

5. **Tutorial authoring** — eight `data-status="pending"` cards in the
   help-page Tutorials section need actual step-by-step content. The
   four-turn breeding-card build is now complete, so the
   breeding-card tutorial is the highest-value first one to author.

---

## 7. Honest framing

**What turn 147 + 147b actually delivered:**

- The L3 toolbar's row count goes from 3 (or 2 on widest viewports) to
  **1 by default**. Power users who need the secondary controls click
  `⋯ more` once and the cluster reveals; their choice persists across
  reloads.
- The Windows-(N) button now correctly reflects all five step modes,
  including the bug case (`winN`) Quentin flagged.
- Both fixes preserve every event listener and feature; the changes are
  purely visual + state-reflection. No workflow had to change.
- Two pre-existing tests had to be updated: turn 128 (to invert
  assertions that codified the very bug being fixed) and turn 134 (a
  brittle 600-char lookahead window).
- The full test suite grew by 93 tests and stays at 0 failures.

**What it deliberately didn't deliver:**

- A bigger restructure of the L3 toolbar (e.g. moving secondary
  controls to a dedicated settings popover). The disclosure pattern is
  the smallest possible visual fix that addresses Quentin's specific
  complaint without touching muscle memory.
- Anything beyond the cycler button for the `winN` case. The sidebar
  step-mode bar (`#stepModeBar`) and the L3 compare-unit
  (`#l3CompareUnit`) already had their own visual feedback.

**Manuscript impact:** small but real. Page 1 (the discovery page) is
where Quentin spends most of his time. Reclaiming ~120px of vertical
space and giving the cycler honest visual feedback both reduce the
friction of every minute of every session he runs against the atlas.

Walk the map carefully, respect cohort discipline, don't break the
test suite.
