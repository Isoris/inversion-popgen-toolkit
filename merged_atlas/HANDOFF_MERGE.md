# MERGE chat — repair + integrate the 5 batch outputs

**For the merge chat that runs after batches 1–5.**

The 5 parallel batch chats are done. Their outputs have been pre-merged into the bundled `Atlas/` folder, and their TODO inventories are pre-computed in `Atlas/docs/merge_inputs/`. Your job: resolve the TODOs, ship the final repo.

---

## State on arrival (already verified)

```
Foundation tests:  455 / 455 ✓   (Step 0–2.A intact, bit-identical to legacy)
Batch tests:       194 / 194 ✓   (every page module imports cleanly)
TOTAL:             649 / 649 ✓
```

The 5 batch chats did clean work. **No one broke `shared/`** (read-only rule held). Every page module imports without throwing. You are starting from a green tree.

What's left: **resolve TODO markers** and **assemble 4 working sub-atlas HTMLs**.

---

## Bundle structure on arrival

```
.
├── HANDOFF_MERGE.md                     ← this file
├── HANDOFF_SUPPLEMENT.md                ← important corrections + lessons learned
├── Atlas/
│   ├── shared/                          (locked, READ-ONLY for you too unless promoting)
│   │   ├── contingency.js, hungarian.js, het_rate.js, kmeans.js,
│   │   ├── per_l2_cluster.js, state.js, state_io.js
│   ├── inversion_discovery/
│   │   ├── page1.{js,html}              (4361 LOC JS, 89 TODOs — the monster)
│   │   ├── page2.{js,html}              (251 LOC, 46 TODOs)
│   │   ├── page8.{js,html}, page12.{js,html}, page15.{js,html}, page19.{js,html}
│   │   └── BATCH_1_NOTES.md             ← read this for batch-1 context
│   ├── inversion_review/
│   │   ├── page11.{js,html}, page4.{js,html}, page6.{js,html}, page7.{js,html},
│   │   ├── page_sv_evidence.{js,html}
│   │   └── BATCH_2_NOTES.md
│   ├── inversion_catalogue/
│   │   ├── page3.{js,html}, page9.{js,html}, page10.{js,html}, page17.{js,html},
│   │   ├── page18.{js,html}, page21.{js,html}, page_overview.{js,html}
│   │   └── BATCH_3_NOTES.md, BATCH_4_NOTES.md
│   ├── inversion_comparative/
│   │   ├── page5.{js,html} ← STUB! help page, NOT multi-species. See HANDOFF_SUPPLEMENT.md
│   │   ├── page16.{js,html}, page16b.{js,html} ← multi-species lives in page16b
│   │   └── BATCH_5_NOTES.md
│   ├── tests/
│   │   ├── test_shared_*.js (8 files)   ← READ-ONLY foundation regression tests
│   │   ├── test_modular_smoke.js, test_legacy_parity.js  ← READ-ONLY
│   │   └── test_<phase>_page<N>.js (21 files)  ← from the batches
│   └── docs/
│       ├── MIGRATION_INVENTORY.md
│       ├── MIGRATION_LOG.md             ← update at end with what shipped
│       └── merge_inputs/                ← PRE-COMPUTED inventories (use these!)
│           ├── TODO_MISSING_inventory.txt   197 markers, frequency-sorted
│           ├── TODO_PROMOTE_inventory.txt   2 markers
│           ├── TODO_SLOTS_inventory.txt     31 missing slots
│           └── all_todo_locations.txt       Every TODO with file:line
├── legacy/
│   └── Inversion_atlas.html             3.5 MB legacy (turn 165 close binary)
└── handoff_docs/                        Original Quentin handoffs + SPEC BLOCK 2
```

---

## Verify state on first turn

```bash
cd Atlas
# Foundation regression
PASS=0; FAIL=0
for t in tests/test_shared_*.js tests/test_modular_smoke.js tests/test_legacy_parity.js; do
  out=$(LEGACY_ATLAS=$PWD/../legacy/Inversion_atlas.html node "$t" 2>&1)
  if echo "$out" | grep -q "fail: 0"; then
    n=$(echo "$out" | grep -oE "pass: [0-9]+" | head -1 | grep -oE "[0-9]+")
    PASS=$((PASS + n))
  fi
done
echo "FOUNDATION: $PASS pass, expecting 455"
# → expect: FOUNDATION: 455 pass

# Batch tests
PASS=0
for t in tests/test_discovery_*.js tests/test_review_*.js tests/test_catalogue_*.js tests/test_comparative_*.js; do
  out=$(node "$t" 2>&1)
  if echo "$out" | grep -q "fail: 0"; then
    n=$(echo "$out" | grep -oE "pass: [0-9]+" | head -1 | grep -oE "[0-9]+")
    PASS=$((PASS + n))
  fi
done
echo "BATCHES:    $PASS pass, expecting 194"
# → expect: BATCHES: 194 pass
```

If anything mismatches, STOP and investigate before extracting anything new.

---

## TODO bucket distribution (already analyzed)

```
High-frequency (≥3 refs)    → 2 markers   → promote to shared/
Medium (2 refs each)        → 15 markers  → triage case by case
Low (1 ref only)            → 180 markers → keep in the page that needs them
Missing state slots         → 31          → add to SLOT_REGISTRY
```

The work is **bounded and small** — most TODOs are one-off references that just need their function pulled from legacy into the page that calls them. Only ~17 functions need to go to `shared/`.

---

## Step-by-step protocol

### Step 1 — Promote `_esc` to `shared/` (10 min)

This is the highest-impact single action. `_esc` is HTML-escape — used by 4+ batches. Extract from `Inversion_atlas.html` line ~13925 (per BATCH_5_NOTES.md):

```js
// Atlas/shared/dom.js
export function _esc(s) { /* extracted body, paste verbatim from legacy */ }
// Maybe also: escapeHtml (alias used in some legacy code)
export const escapeHtml = _esc;

if (typeof window !== 'undefined') {
  window._esc = _esc;
  window.escapeHtml = escapeHtml;
}
```

Then in every page module that has `// TODO_MISSING(_esc)`, add:
```js
import { _esc } from '../shared/dom.js';
```
…and remove the TODO comment.

Add a parity test in `tests/test_legacy_parity.js` for `_esc` (extend the existing `wanted` list and add a comparison block).

### Step 2 — Resolve the medium-frequency 15 markers

Look at `docs/merge_inputs/TODO_MISSING_inventory.txt`, lines 2–17. Each appears 2–3 times. Decide per item:

- **Used by 2 different sub-atlases** → promote to `shared/<thematic>.js`
- **Used by 2 pages within ONE sub-atlas** → put in a per-sub-atlas helpers file, e.g. `inversion_discovery/_helpers.js`

For each promoted function: extract from legacy, add to module, write a short parity test for it, update imports in the consuming pages.

### Step 3 — Resolve the low-frequency 180 markers

For each, find the function in legacy (`docs/merge_inputs/all_todo_locations.txt` has the page:line for every TODO; grep legacy by function name to find its definition). Extract the body, paste it into the page module that needs it. Don't try to refactor — just paste verbatim and fix syntax if needed.

This is mechanical work — many TODOs are tiny utility functions (`xToPx`, `yToPx`, `withAlpha`, `formatTrackVal`, etc.). Some are larger render functions. Either way: copy from legacy, paste in page, remove TODO comment.

### Step 4 — Add the 31 missing state slots

For each entry in `docs/merge_inputs/TODO_SLOTS_inventory.txt`, add to `SLOT_REGISTRY` in `shared/state.js`. Pick a tag based on usage:

- `state._foo` (single underscore) → typically `transient` (cache, geometry, UI scratch)
- `state.fooMode`, `state.fooEnabled` → typically `persisted` (with persistKey)
- Slots referenced by `candidate*` family → typically `cross_atlas` (already mostly there)

The slot list looks like:
```
state._l3CacheFp, state._l3CacheRendered, state._lineageComputeScheduled,
state._simGeom, state._simMinimapGeom, state._thSimGeom, state.ancestryPalette,
state.cacheKey, state.candidateMode, state.candidatePageMode, ...
```

Most are clearly transient (UI cache, geometry). A few might warrant cross_atlas (e.g., `candidateMode` if it carries between atlases).

Update `tests/test_shared_state.js` if you add cross_atlas slots — round-trip test should still pass.

### Step 5 — Build the 4 sub-atlas HTMLs

For each sub-atlas, create `Atlas/inversion_<phase>.html` containing:

1. **Boilerplate**: doctype, charset, viewport, title.
2. **CSS**: copy from legacy `Inversion_atlas.html` lines 1–~5000 (the entire CSS block — start at the first `<style>` tag, end at `</style>`). All 4 sub-atlases use the SAME CSS for consistency.
3. **Tab bar**: the navigation bar HTML, but **filtered to this sub-atlas's pages only**. Look at legacy lines ~5020–5080 — each `<button data-page="pageN" data-stage="discovery|refinement|...">` maps to a stage. Keep only buttons whose stage is owned by this sub-atlas:

   | Sub-atlas | Stages to keep |
   |---|---|
   | inversion_discovery | discovery, refinement |
   | inversion_review | refinement, classification |
   | inversion_catalogue | synthesis, help |
   | inversion_comparative | compare |

4. **Page divs**: concatenate all `Atlas/inversion_<phase>/page*.html` for that sub-atlas.
5. **Script module**: a single `<script type="module">` block:
   ```js
   // Import all page modules
   import * as page1 from './inversion_discovery/page1.js';
   import * as page2 from './inversion_discovery/page2.js';
   // ... etc
   import { makeState, readPersistedSlots } from './shared/state.js';
   // ... other shared

   // Build state
   const state = makeState();
   Object.assign(state, readPersistedSlots());

   // Tab routing
   document.querySelectorAll('#tabBar button[data-page]').forEach(btn => {
     btn.addEventListener('click', () => {
       const pageId = btn.dataset.page;
       document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
       document.getElementById(pageId)?.classList.add('active');
       document.querySelectorAll('#tabBar button').forEach(b => b.classList.remove('active'));
       btn.classList.add('active');
       // Dispatch render. Each page exports a renderPage<N>(state) by convention.
       const renderFn = ({ page1: page1.renderPage1, page2: page2.renderPage2, /* ... */ })[pageId];
       if (typeof renderFn === 'function') renderFn(state);
     });
   });

   // Initial render — first page
   /* trigger first tab click */
   ```

   Each page's batch-extracted module should already export a `renderPage<N>` (or similar) function. Check `BATCH_*_NOTES.md` for export conventions.

6. **Save**.

### Step 6 — Browser sanity check

Open each `Atlas/inversion_<phase>.html` in a browser:
- The tab bar should appear with only that sub-atlas's tabs.
- Clicking a tab should show that page's content (even if the content is a blank canvas).
- No JavaScript console errors (or only data-load errors, which are expected without precomp JSONs).

You don't need pixel-perfect parity. You need: **structure works, modules load, no syntax errors**.

### Step 7 — Final test run

```bash
cd Atlas/
PASS=0; FAIL=0
for t in tests/*.js; do
  out=$(LEGACY_ATLAS=$PWD/../legacy/Inversion_atlas.html node "$t" 2>&1)
  if echo "$out" | grep -q "fail: 0"; then
    n=$(echo "$out" | grep -oE "pass: [0-9]+" | head -1 | grep -oE "[0-9]+")
    PASS=$((PASS + n))
  else
    n=$(echo "$out" | grep -oE "fail: [0-9]+" | head -1 | grep -oE "[0-9]+")
    FAIL=$((FAIL + n))
    echo "FAILED: $t"
  fi
done
echo "TOTAL: $PASS pass, $FAIL fail"
# Expected: ≥649 pass (foundation 455 + batches 194), 0 fail
# Plus any new shared parity tests you added (e.g., _esc parity)
```

### Step 8 — Update MIGRATION_LOG.md

Append a section documenting:
- Step 2.B/C/D done (extracted via parallel batches + merge)
- Step 3 done (review pages migrated)
- Step 4 done (catalogue pages migrated)
- Step 5 done (comparative pages migrated)
- New shared modules added (`shared/dom.js`, etc.)
- Final test counts
- What's deferred for follow-up (see "deferred" section below)

### Step 9 — Bundle and ship

```bash
cd /tmp
tar czf /mnt/user-data/outputs/Atlas_FINAL_2026-05-05.tar.gz Atlas/
```

---

## What's likely DEFERRED (don't try to do this turn)

These are real but post-merge work:

- **Pixel-perfect parity** with the legacy single-file Atlas. Some interactions (e.g. Shift+click lasso, complex keyboard shortcuts) may behave slightly differently. Document; don't fix in this turn.
- **The dispatcher** — legacy `Inversion_atlas.html` has scattered tab routing (5 different sites that each toggle `.page.active`). The unified dispatcher you write should handle ~80% of cases; the long tail of "promotions auto-jump to page2" etc. is worth deferring.
- **Tab numbers** — legacy tabs are numbered (1, 2, 2b, 3, 4, 5, 5b, 6, 7, 8, 9, 10, 11, 16). The numbers are preserved per page in `data-page="pageN"`. The tab labels (`<span class="num">N</span>`) should be preserved as well. If the numbering looks weird in a particular sub-atlas after filtering, just leave it; Quentin will adjust.
- **CSS scoping** — all 4 sub-atlases sharing ALL legacy CSS is overkill. Some optimization possible but DEFER.
- **Window-mounted globals** — page modules expose some `window.*` globals (e.g., `window._csBpJumpToWindow` per BATCH_5_NOTES.md). These are CROSS-PAGE bridges that matter when discovery and comparative both load. For single-sub-atlas HTMLs they're harmless. Document but don't refactor.

---

## What "done" looks like

- ✅ All 4 sub-atlas HTMLs produced and **load in a browser** (tabs visible, modules import without console errors)
- ✅ All TODOs resolved (or explicitly documented as deferred in `MIGRATION_LOG.md`)
- ✅ All 31 missing slots added to `SLOT_REGISTRY`
- ✅ `_esc` and any other promoted shared modules have parity tests
- ✅ Foundation tests: still 455 pass, zero broken
- ✅ Batch tests: still 194 pass, zero broken
- ✅ MIGRATION_LOG.md updated
- ✅ Final tarball: `Atlas_FINAL_2026-05-05.tar.gz` shipped to `/mnt/user-data/outputs/`

---

## Honest scope estimate

This is **6–10 turns of focused work**. The breakdown:

| Step | Turns |
|---|---|
| Verify state, read inventory | 1 |
| Promote `_esc` to shared/ | 1 |
| Resolve medium-frequency 15 markers | 1–2 |
| Resolve low-frequency 180 markers (mechanical) | 2–3 |
| Add 31 missing slots | 1 |
| Build 4 sub-atlas HTMLs | 1–2 |
| Browser sanity, final tests, bundle | 1 |

The mechanical TODOs (Step 4 above) are the bulk of the time. Many are tiny — paste a 5-line function from legacy into the page that needs it. Some require thought (e.g., what does `_drawBandTraceStrip` need from state?), but you have all the context: legacy lines are searchable, batch notes document conventions, the foundation modules cover the heavy infrastructure.

If you find yourself spending more than 1 turn on a single TODO, consider deferring it — mark `// DEFERRED: not resolved by merge chat — needs attention` and move on. Quentin will iterate.

---

## Final reminder

The 5 batch chats kept their hands off `shared/`. You should too unless you're explicitly promoting a function that has 3+ references. The foundation's 455 parity tests are the contract — they MUST keep passing throughout your work.

If a foundation test fails after one of your edits, undo immediately. The right impulse is: "did I accidentally edit a file in `shared/`?" — yes 90% of the time.

Good luck.
