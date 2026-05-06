# HANDOFF — turn 154 — G-panel inheritance compute UX hardening

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,861 lines, +159 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 153 handoff. Closes the UX hole Quentin reported:
clicking ↻ recompute appeared to do nothing because compute could
return null silently and the tab body couldn't tell the user why.

---

## 0. What Quentin asked for

> "when we click compute for inheritance so far it didn't work because
> the UI wasn't ready as well so let's manage it now."

The compute path itself worked, but:
1. **Tab didn't auto-compute on first open** — every fresh open
   required a manual click even when items were stable.
2. **No visible click feedback** — synchronous compute is fast, the
   button never changed during compute, click felt unregistered.
3. **`runInheritanceCompute` returned `null` silently** for several
   distinct reasons (insufficient items, all-K=0, exception in
   `inheritanceGroupClustering`) — all rendered as the same generic
   "No compute result yet" placeholder, indistinguishable from
   "you haven't clicked yet".
4. **No success indicator** — even when compute succeeded, no status
   line confirmed it. Cards just appeared, or didn't.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse,
pragmatic. Tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 New state slot

```
state.inheritanceLastStatus: {
  ok:        boolean,
  reason:    'computed' | 'cached' | 'insufficient_items' |
             'null_result' | 'exception',
  n_items:   number,
  n_groups:  number,        // when ok=true
  threshold: number,        // when applicable
  error:     string,        // when reason='exception'
  message:   string,        // human-readable, displayed in pill
  at:        number,        // Date.now() of last update
} | null
```

Set on every code path inside `runInheritanceCompute`. Cleared by
`invalidateInheritanceCache`. Read by the inheritance tab to render
the status pill.

### 2.2 `runInheritanceCompute` rewritten with status tracking

Every early-return path now writes a status. The
`inheritanceGroupClustering` call is wrapped in try/catch so
exceptions are captured into the status instead of bubbling. Each
status carries a precise `reason` enum and a human-readable
`message`.

The five outcomes:

- **`computed`** — fresh compute returned a non-null result. Status
  carries `n_groups` and `threshold`.
- **`cached`** — cache key matched, returned the cached result. Status
  carries `n_groups` from cached result.
- **`insufficient_items`** — fewer than 2 candidates have locked
  labels. Status `ok=false`, message explains the gate rule.
- **`null_result`** — `inheritanceGroupClustering` returned null
  (typically all candidates have K=0, so band_index is empty). Status
  message hints at the K=0 cause.
- **`exception`** — `inheritanceGroupClustering` threw. Status
  captures the error message string.

### 2.3 Tab body auto-computes on mount

```
if (items.length >= 2 && typeof runInheritanceCompute === 'function') {
  const expectedKey = _inheritanceCacheKey(items, mode);
  const cacheStale = (expectedKey == null) ||
                     (state.inheritanceCacheKey !== expectedKey) ||
                     !state.inheritanceResult;
  if (cacheStale) {
    runInheritanceCompute({ threshold: threshold });   // no force — cache hit OK
  }
}
```

No `force: true` because the cache hit is a valid no-op outcome (and
still writes a `cached` status). The tab now mirrors what the
lines-strip already did in its draw-time stale-cache guard, so all
inheritance consumers behave consistently.

### 2.4 Status pill in the tab UI

Above the threshold slider row, a colored pill displays the most
recent compute outcome:

- Green left-border + ● icon when `ok=true`
- Grey + ○ when `reason='insufficient_items'` (informational, not error)
- Amber + ⚠ when compute failed (exception or null_result)

The pill text is `status.message`, e.g. `"computed: 4 groups across
12 candidates @ threshold 0.15"` or `"compute returned null (no
usable bands; check K values on items)"`.

### 2.5 Distinct messages for 0-groups vs null-result

Pre-turn-154 the body had two states: result-with-cards or
no-result-placeholder. Now three:

- **`result && n_groups >= 1`** → render cards (unchanged)
- **`result && n_groups === 0`** → "Compute succeeded but produced 0
  inheritance groups. Try lowering the threshold (currently 0.15) ..."
- **`!result && status.ok=false`** → "No groups to display — see
  status above."
- **`!result && !status`** → original placeholder (defensive
  fallback; auto-compute should make this branch unreachable post-turn-154)

### 2.6 Recompute button click feedback

Clicking ↻ recompute now:
1. Captures original button text
2. Sets text to `"computing…"` and `disabled = true`
3. Runs compute synchronously
4. Calls `_renderGPanelModal()` — re-render replaces the entire modal
   so the button is recreated fresh (text + disabled state reset
   automatically)
5. If re-render fails, restores the original button manually so the
   user can retry

Compute is synchronous so the flash is brief on small inputs; on
larger candidate sets (≥ 50 items) it'll be more visible. Either way
the click is visibly acknowledged.

---

## 3. What did NOT change

- **`inheritanceGroupClustering`** — untouched. Still returns null on
  insufficient bands; the difference is `runInheritanceCompute` now
  catches and reports.
- **`_gpInhRenderResult`** — untouched.
- **`_gpInhRenderEmpty`** — untouched.
- **Threshold slider debounce path** — untouched. Slider was already
  visible-feedback through the live value label; the recompute button
  was the silent one.
- **Lines-strip pills** — untouched. They were already using the
  draw-time stale-cache guard.
- **Turn 153 auto-register hook** — untouched. Still fires on
  candidate mutation; turn 154 just makes the UI more honest about
  what compute did or didn't do.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 153)   | 72,702  | 2685            | 57    |
| Post-turn-154 (current)  | 72,861  | 2746            | 58    |
| **Δ session**            | +159    | +61             | +1    |

Full sweep: **2746 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn154_compute_ux_hardening.js` (61 tests) covers:

1. `runInheritanceCompute` writes status on all five paths
   (insufficient_items, cached, exception, null_result, computed)
2. `invalidateInheritanceCache` clears the status
3. Tab body auto-computes on mount (uses cache-key staleness check,
   does NOT pass force: true)
4. Status pill renders with correct color, icon, and message per
   outcome
5. Distinct messages for 0-groups vs null-result vs never-computed
6. Recompute button flash + disable + restore on failure
7. Sandboxed end-to-end status-write behaviour for all five paths,
   including timestamp recency
8. Existing flow preserved — earlier turn contracts (130, 152, 153)
   still in place

No existing test needed inverting this turn. The runtime behaviour
that was codified by previous tests (the placeholder text on first
open, etc.) is now defensive-fallback and still emitted in the
no-status branch — the assertions still pass.

---

## 5. What Quentin can do next session

To exercise:

1. Open the G-panel (`g`) → switch to inheritance tab.
2. **First time** — with ≥ 2 candidates: tab auto-computes and shows
   a green status pill `"computed: N groups across M candidates @
   threshold 0.15"` plus the cards.
3. **Less than 2 candidates** — grey status pill `"○ Need ≥2
   candidates with locked labels (you have 1)"` + the existing
   empty-state body.
4. **Click ↻ recompute** — button briefly shows `"computing…"` then
   re-renders to `"↻ recompute"` with a fresh status pill (success
   or failure depending on data state).
5. **Drag the threshold slider** — value label tracks live; debounce
   fires after 220 ms; status pill updates with new threshold.
6. **Threshold below sensible bound** — if compute produces 0 groups,
   you get a distinct amber-accent body explaining "succeeded but 0
   groups, try lowering the threshold".
7. **Edit candidates while popup open** — turn 153 auto-register
   fires, cache invalidates (which also clears status), turn 154
   tab-body auto-compute on next render rebuilds status. Should feel
   instant.

If anything still feels silent or unclear, that's a real bug — let
me know what you see.

---

## 6. Things I almost broke and fixed

- **`runInheritanceCompute` swallowed exceptions silently.** Pre-turn-154
  if `inheritanceGroupClustering` threw, the exception bubbled up to
  the caller — which in the recompute click handler was caught by
  `catch (e) { console.warn(...) }`, never reaching the user. Now
  caught inside `runInheritanceCompute` itself, status set with the
  error message, UI shows it.
- **Auto-compute used `force: true` initially.** That would re-run on
  every tab open, including ones where the cache was perfectly valid.
  Wasteful. Removed `force` — cache hit is fine, status writes
  `reason: 'cached'` so the user still sees confirmation.
- **Status pill on insufficient_items conflicted with empty-state
  body.** Both wanted to say "you don't have enough candidates".
  Resolved: pill is grey ("informational"), empty-state body
  elaborates with the gate rule. Two complementary surfaces, not a
  duplication.
- **Button flash flickered when re-render was fast.** Synchronous
  compute can finish in < 1 ms for small inputs; the button text
  flipped to "computing…" and back so fast it looked like a no-op.
  Decided to live with this: the flash IS the flash. On larger
  inputs (≥ 50 candidates) it'll be visible. Adding a forced minimum
  delay would be theatre.
- **`inheritanceLastStatus` not cleared on invalidate.** Without the
  clear, after a candidate-mutation the user would see a stale "ok:
  computed" status pill describing data that no longer existed.
  Fixed in `invalidateInheritanceCache`.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 154 patched, 72,861 lines.
- `tests/test_turn154_compute_ux_hardening.js` — NEW, 61 tests.
- `HANDOFF_2026-05-05_turn154_compute_ux_hardening.md` — this file.

Plus prior handoffs (carried for reference):
- `HANDOFF_2026-05-05_turn153_inheritance_auto_register.md`
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md`
- `HANDOFF_2026-05-05_turn151_slab_het_coloring.md`

---

## 8. Honest framing

**What's solid:**
- Every compute outcome is now observable. The user can't click a
  button and see nothing — there's always a status pill change.
- Auto-compute on mount removes the "click first" friction. The tab
  is now a passive view, not a manual harness.
- The `state.inheritanceLastStatus` slot is also available to other
  consumers (lines-strip tooltip, matrix popup) if we want to surface
  the status there too. Not done this turn but trivial later.
- Test coverage for the runtime status-write behaviour is end-to-end
  — every code path has a sandboxed assertion proving it sets the
  expected status.

**What's NOT done (and why that's right):**
- **Status pill on the matrix popup.** Currently the matrix popup
  uses its own `imxStatus` element which says "N candidates" or
  "(need ≥2 candidates)". Not strictly the same as the new compute
  status. Could unify in a follow-up but they serve slightly
  different purposes (matrix shows "I have N items to render", G-panel
  shows "compute succeeded with M groups").
- **No "compute is taking a while" spinner for very large inputs.**
  Compute is synchronous. For 100+ candidates with K=6 each, it could
  take long enough to feel laggy. Real fix would move the compute to
  a worker, which is a larger turn.
- **No history of past computes.** Status only stores the most
  recent. If a user wanted to compare "before / after threshold
  change", they'd need to remember. Probably fine.
- **Threshold still not in `_inheritanceCacheKey`.** Same workaround
  as turn 152 — `_gpInhSetThreshold` invalidates explicitly. ~10 LOC
  honesty fix queued as option E.

**What's queued (from turn 153 §8):**
- B — FIG_C30 V-shape `(u, agreement_fraction)` plot
- E — threshold in `_inheritanceCacheKey` (~10 LOC honesty fix)
- D — pivot based on what het/inheritance reveals
- A (UV refactor), C (lines het), F (per-group fish list)

End of handoff.
