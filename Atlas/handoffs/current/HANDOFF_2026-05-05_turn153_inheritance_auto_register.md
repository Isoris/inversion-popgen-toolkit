# HANDOFF — turn 153 — Auto-register inheritance on candidate-list mutations

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,702 lines, +109 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 152 handoff. This turn closes the consistency gap
flagged in turn 152 §6 — the G-panel inheritance tab and the matrix
popup were not auto-redrawing when the candidate set changed.

---

## 0. What Quentin asked for

> "we must solve the G pill panel first because it must totally
> register inversion candidates automatically"

Pre-turn-153, the G-panel inheritance tab and matrix popup read
`state.inheritanceResult` directly. Lines-strip pills had a draw-time
stale-cache guard that auto-recomputed on cache-key mismatch — the
G-panel and matrix did not. So promoting a candidate elsewhere left
the G-panel showing stale group cards until the user did something to
re-render it.

Fix this turn: move invalidation upstream to the candidate-mutation
funnel. Every add / remove / import / L3-commit / refinement that
funnels through `persistCandidateList` now invalidates the inheritance
cache when the candidate set has actually changed and notifies all
open consumers to re-render.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse,
pragmatic. Tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 New helper — `_autoRegisterInheritanceOnCandidateChange()`

Single-purpose: compare current cache key vs expected, invalidate +
notify on mismatch, no-op on match.

```
_autoRegisterInheritanceOnCandidateChange()
  ├── _gatherActiveCandidatesForInheritance()
  ├── _inheritanceCacheKey(items, mode) → expectedKey
  ├── if (cachedKey === expectedKey) return        ← idempotent
  ├── invalidateInheritanceCache()
  └── debounced via _inheritanceNotifyScheduled →
        queueMicrotask(_notifyInheritanceConsumers)
```

Cheap: one cache-key compute (≤ a few hundred candidates × O(label
fingerprint hash)). Idempotent: no-op when the key matches what's
cached. Sandbox-safe: missing helpers fail fast before any state
mutation.

### 2.2 New helper — `_notifyInheritanceConsumers()`

Re-renders all open inheritance consumers. Each is independently
guarded; failure in one does not block the others.

```
_notifyInheritanceConsumers()
  ├── window.requestRepaint?.()                    ← lines-strip pills
  ├── if (gPanelOpen && gPanelTab === 'inheritance')
  │     _renderGPanelModal()                       ← G-panel tab
  └── if (matrixModal && visible)
        _redrawInheritanceMatrix()                 ← matrix popup
```

Each consumer call is wrapped in try/catch so partial notification is
better than none.

### 2.3 Wired into `persistCandidateList()`

```
function persistCandidateList() {
  if (!state.data) return;
  if (typeof _rebuildCandidateRegistries === 'function') { ... }   // turn 129
  // turn 153: invalidate + notify on candidate set change
  if (typeof _autoRegisterInheritanceOnCandidateChange === 'function') { ... }
  try { localStorage.setItem(...); } catch (e) { ... }
}
```

Hook is positioned **after** `_rebuildCandidateRegistries` (so the
gatherer sees the up-to-date `state.candidates` dict) and **before**
the localStorage write (deterministic order; localStorage is a side
effect, not a dependency). Both calls are guarded with `typeof` and
wrapped in try/catch — non-blocking even in degraded environments.

### 2.4 Why this single change covers all paths

Every candidate mutation in the codebase funnels through
`persistCandidateList()`:

- `addCandidateToList` → `persistCandidateList`
- `removeCandidateFromList` → `persistCandidateList`
- `removeCandidateFully` → `removeCandidateFromList` → `persistCandidateList`
- `commitL3Draft` → `persistCandidateList` (turn 129 contract)
- JSON import handler → `persistCandidateList`
- L2-sweep auto-promote → `persistCandidateList`
- Refinement / regime updates → `persistCandidateList`

Adding the hook here covers all of them at once.

---

## 3. What did NOT change

- **`runInheritanceCompute` / `inheritanceGroupClustering`** — untouched.
- **`_inheritanceCacheKey`** — same shape (still doesn't include
  threshold; that workaround stays in `_gpInhSetThreshold` for now).
- **`_drawInheritanceLabelsStrip` draw-time stale guard** — still
  there. Belt + suspenders: the upstream invalidate is the new primary
  path, the draw-time guard remains as backup.
- **Matrix popup recompute button** — still works for forced refreshes.
- **G-panel Slice 3 implementation (turn 152)** — untouched.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 152)   | 72,592  | 2630            | 56    |
| Post-turn-153 (current)  | 72,702  | 2685            | 57    |
| **Δ session**            | +109    | +55             | +1    |

Full sweep: **2685 / 0** across `tests/test_turn*.js`. JS syntax: clean.
HTML parser: 0 errors.

`tests/test_turn153_inheritance_auto_register.js` (55 tests) covers:

1. Helper declarations + window exports
2. `persistCandidateList` wiring (placement, guards, try/catch)
3. `_autoRegister` body — gather call, cache-key compare, idempotent
   no-op, debounce flag, queueMicrotask + setTimeout fallback,
   sandbox-safety via typeof checks
4. `_notifyInheritanceConsumers` body — three independent consumer
   branches, each guarded, each wrapped in try/catch
5. Sandboxed `_autoRegister` behaviour — no-op on key match, invalidate
   on mismatch, debounce coalesces 3 rapid calls to 1 microtask, flag
   resets after fire, missing-helpers fail-soft
6. Sandboxed `_notifyInheritanceConsumers` behaviour — `requestRepaint`
   calls, G-panel modal re-render conditional on open + on inheritance
   tab, matrix popup re-draw conditional on visibility, all-three-fire
   end-to-end test, partial-failure isolation (thrown error in one
   consumer does not block others)
7. Existing flow preserved

No existing test needed inverting this turn — `persistCandidateList`'s
contract was extended (new hook added before localStorage write), not
changed.

---

## 5. What Quentin can do next session

To exercise:

1. Open the G-panel (`g`), switch to the **inheritance** tab, click
   ↻ recompute → cards appear.
2. Leave the popup open, switch to a different page or stay where you
   are, and **promote a new candidate** (page 1 → lock + promote).
3. The G-panel inheritance tab should **automatically re-render** with
   the new candidate folded into the group cards. No "↻ recompute"
   click needed.
4. Same flow for **removing** a candidate (page 2 ✕, page 4 red ✕,
   etc.) — the inheritance tab should reflect the removal immediately.
5. Same flow for the **matrix popup** — open it, then promote a
   candidate; the matrix should resize and re-draw with the new row
   and column included.
6. Lines-strip pills already worked pre-turn-153 (the existing
   draw-time stale guard). They still work after turn 153; the new
   upstream path runs first and `requestRepaint` is the trigger.

If any consumer doesn't refresh, that's a real bug — let me know what
broke.

---

## 6. Things I almost broke and fixed

- **`_autoRegister` calling `_gatherActiveCandidatesForInheritance` in
  a sandbox without it.** Each helper is now individually guarded with
  `typeof === 'function'` checks at the top of `_autoRegister` so the
  function fails fast and doesn't mutate state when its dependencies
  are missing.
- **Microtask debounce flag could leak.** If `_notifyInheritanceConsumers`
  threw, the flag would never reset and subsequent mutations wouldn't
  schedule new notifies. Fixed by setting the flag to false at the
  *start* of the fire callback (before the try/catch), not the end.
- **G-panel modal re-rendering when on the wrong tab.** Calling
  `_renderGPanelModal` from notify when the user is on the manual or
  karyotype tab would re-render those tabs unnecessarily and visually
  flicker. Now guarded by `gPanelTab === 'inheritance'`.
- **Matrix popup hidden but in the DOM.** First DOM-only check would
  redraw a hidden modal. Now also checks
  `style.display !== 'none'` so the redraw only fires when the modal
  is actually visible.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 153 patched, 72,702 lines.
- `tests/test_turn153_inheritance_auto_register.js` — NEW, 55 tests.
- `HANDOFF_2026-05-05_turn153_inheritance_auto_register.md` — this file.

Plus the prior handoffs (carried for reference):
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md`
- `HANDOFF_2026-05-05_turn151_slab_het_coloring.md`

---

## 8. Honest framing

**What's solid:**
- The fix is a single funnel hook, not a sprawling change. One
  `persistCandidateList` mutation = one auto-register call = one
  notify pass (debounced). Tests cover all the moving parts.
- Idempotent: no-op when nothing changed, so calls from non-list
  paths (e.g. favorites toggles in Save-Session payload) don't churn
  the inheritance cache.
- Sandbox-safe: every helper reference is typeof-guarded. Tests
  prove the function does not throw or mutate state when its
  dependencies are absent.

**What's NOT done (and why that's right):**
- **Threshold still not in cache key.** Same workaround as turn 152.
  `_gpInhSetThreshold` calls `invalidateInheritanceCache()` directly.
  A proper fix would extend `_inheritanceCacheKey` to fold threshold
  in. ~10 LOC, separate turn.
- **No active-page awareness.** The notify pass calls
  `requestRepaint` unconditionally. If the user is on a page that
  doesn't show the lines strip, the call is wasted. Cheap (it's a
  single function call), so probably fine.
- **No notify for closed-but-cached G-panel state.** If the user
  closes the popup, candidates change, then re-opens — the popup
  re-renders fresh on open, which is correct. No work needed.

**What's queued (from turn 152 §8):**
- A — UV refactor on locked groups
- B — FIG_C30 V-shape `(u, agreement_fraction)` plot
- C — Per-sample-lines het coloring
- D — Pivot based on what het / inheritance reveals
- E — Threshold in inheritance cache key (~10 LOC honesty fix)
- F — Per-group "show fish" expand toggle (~50 LOC)

The G-panel auto-register is now done. Quentin should exercise it on
real data and tell me whether B (next bounded scientific deliverable)
or D (a pivot based on what the inheritance + het reveal) makes more
sense next session.

End of handoff.
