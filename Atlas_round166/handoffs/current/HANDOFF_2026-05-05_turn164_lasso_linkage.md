# HANDOFF — turn 164 — fish-set linkage table

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (76,034 lines, +528 LOC from 75,506)
**Working dir**: `/home/claude/work/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes** Slices 1 + 3 of `SPEC_lasso_inheritance_backgrounds.md`
(the purity computation + the linkage-table sidebar). Slices 2 (alpha
intervals on per-sample-lines), 4 (enlarged page), and 5 (candidate-page
integration) deliberately deferred — each needs design conversations
about visual layout and the lasso surface itself isn't yet built in
the codebase.

**Picked up from**: post-turn-163 working tree (3411 / 0 baseline).

---

## 0. Why this spec, this turn

The band-trace work (turns 160–163) hit a natural pause point — the
remaining items (run brackets, threshold calibration) need Quentin's
real-LG28 eyes before being safe to ship. Continuing to push on
band-trace without his calibration data would optimize for synthetic
fixtures.

So I picked the next item using these criteria:
- Self-contained (no JSON pre-compute, no HPC dependency)
- Manuscript-relevant (feeds the inversion paper directly)
- Not gated on Quentin's real-data review
- Right scope for one turn

`SPEC_lasso_inheritance_backgrounds.md` slices 1 + 3 fit perfectly. The
linkage table is the **inverse direction** of the band-trace pipeline:

| Direction | Operates on | Asks |
|---|---|---|
| Band-trace (turns 160-163) | Every L2 envelope along the chromosome | "Where do these fish co-segregate?" |
| **Linkage (this turn)** | Saved candidates (sparse, user-promoted) | "Which band of each promoted candidate do these fish predominantly land in?" |

The two are complementary. They share `state.bandTraceFishSet` (turn
161) as the input, so when the user clicks 🔍 trace, the linkage
table is immediately populated for that fish-set across every
confirmed candidate.

The lasso integration itself (turn 130 trigger) is deliberately
**not** built in this turn — `state.lassoSelection` doesn't exist in
the codebase yet, only `state.linesLassoRect`. Bridging those is its
own scope. For now the fish-set source is just `bandTraceFishSet`,
which Quentin already has a workflow for populating.

---

## 1. What this turn ships

### Pure compute: `_computeLassoLinkage(fishSet, candidateList, opts)`

Given a fish-set (Set or Array of integer sample indices) and a
candidate list, returns a per-candidate purity map:

```js
{
  n_fish_selected: 5,
  n_candidates_seen: 3,
  purity_threshold: 0.7,
  min_band_size: 5,
  per_candidate: {
    'I1': {
      id, chrom, start_bp, end_bp, K,
      best_band, best_purity, n_in_best_band,
      n_lasso_seen,                            // fish-set members assignable in this candidate
      per_band: [{band, n_in_lasso, fraction}],
      is_strong_link,                          // best_purity ≥ threshold AND n_in_best_band ≥ min_band_size
    },
    ...
  },
  strong_links: ['I1', 'I2', 'I3'],            // sorted: purity desc, then n_in_best_band desc, then id asc
}
```

Headless-pure: takes both inputs explicitly so tests can drive it
without building a state object. The state-aware accessor is the
separate `_lassoLinkageGetOrCompute()`.

Skips:
- Unconfirmed candidates (`!c.confirmed`)
- Candidates without `locked_labels`
- Out-of-range fish indices (< 0 or ≥ `locked_labels.length`)
- Fish indices labelled `-1` or ≥ K

Defensive: returns `null` for null fish-set, empty fish-set, null
candidate list, non-array candidate list.

### Cache layer

- `_lassoLinkageCacheKey(chromFilter, fishSet, candidateList)` —
  deterministic key. Order-insensitive in fish-set; sensitive to
  chrom filter, fish-set fingerprint, candidate-list fingerprint
  (computed from the IDs of confirmed candidates with locked_labels,
  so add/remove/promote/demote all bust the cache without needing a
  separate version counter).
- `_lassoLinkageGetOrCompute(opts)` — reads `state.bandTraceFishSet`
  + `state.candidateList`, returns the cached result if the key
  matches; otherwise recomputes. Caches on
  `state.lassoLinkageCache` + `state.lassoLinkageCacheKey`.
- `_invalidateLassoLinkageCache()` — clears both slots.
- **Hooked into `setBandTraceFishSet` (turn 161)** so any fish-set
  change automatically invalidates the linkage cache. Verified by
  test (assertion at the source-pattern level).

### TSV export

- `_lassoLinkageToTSV(result, opts)` — pure serializer. 4 comment
  lines (`# n_fish_selected`, `# n_candidates_seen`, `# purity_threshold`,
  `# min_band_size`) + header row + data rows. Sort: strong links
  first by purity desc, then weak by purity desc, ties broken by
  `n_in_best_band` desc, then `id` asc. Floats to 6 dp; NaN/null →
  empty.
- `_lassoLinkageDownloadTSV()` — Blob/anchor download path with
  same headless-safe contract as turn 162's `_bandTraceDownloadTSV`.
  Filename: `lasso_linkage_<chrom>_n<n_fish>.tsv`. Returns the
  filename in headless environments without trying to call DOM
  APIs.

### Modal popover

`_openLassoLinkagePopover()` opens a fixed-position modal with:
- Title: "Fish-set linkage" + summary stats line
  (`n_fish=… · n_candidates=… · N strong link(s) (purity≥0.7)`)
- Status line explaining sort order
- Table with columns: candidate, chrom, span (Mb), best band (with
  band-color swatch via `_gpKaryoColor`), purity %, n in band /
  n lasso seen, strong?
- Embedded `📊 TSV` button → routes to `_lassoLinkageDownloadTSV()`
- Close button (`close ×`) + close-on-Esc + close-on-backdrop-click
- Footer paragraph reiterating the observation-only framing

When no fish-set is active, the modal still opens but renders an
instructional empty state ("Click 🔍 trace on the lines panel
first") rather than failing silently. This is intentional —
explaining the missing prerequisite is more useful than a no-op.

The modal HTML structure mirrors `openVShapePlot` (turn 156).
Idempotent — first open builds the DOM, subsequent opens reuse it
and just call `_renderLassoLinkageTable()` to refresh.

### Lines header UI

One new button next to `📊 runs`:

```
[ ] band trace    🔍 trace    [largest ▾]    📊 TSV    📊 runs    🔗 linkage
```

`linesBandTraceLinkageBtn` has an informative title and routes its
click to `_openLassoLinkagePopover()` via try/catch.

---

## 2. What this turn does NOT do

Same scope discipline as turns 161–163:

- **Slice 2 (alpha intervals on per-sample-lines).** Visual layout
  on the per-sample-lines panel is its own design pass — already two
  strips (lineage + band-trace) + candidate band highlights occupy
  the upper area. Adding alpha-shaded intervals at every strong-link
  candidate's bp range needs careful thought about Z-order with the
  existing turn 141 candidate-bands feature.
- **Slice 4 (enlarged per-sample-lines page).** A whole new page in
  the atlas — its own multi-turn scope.
- **Slice 5 (candidate-page integration).** Auto-lasso on band click
  + open the modal — depends on slice 2 + 4 being further along.
- **Lasso surface itself.** `state.lassoSelection` doesn't exist;
  only `state.linesLassoRect` and committed-rectangle slots. The
  lasso → fish-set bridge is its own turn. For now,
  `state.bandTraceFishSet` is the canonical input.
- **Genome-wide chromosome navigator** in the modal (spec §2 of the
  enlarged-page slice). Out of scope.
- **Filter-by-purity-threshold UI control** in the modal. The
  threshold is settable via `_computeLassoLinkage(.., {purity_threshold:
  0.5})` from the console; a UI slider is deferred until Quentin
  asks.
- **Click-to-jump on table rows.** Spec §1.4 mentions "Click any row
  → highlight that candidate on the chromosome strip." Deferred
  because the highlight target depends on which page the user is on
  and which chromosome is loaded — needs design.

---

## 3. Files touched

```
Inversion_atlas.html                                +528 LOC
  - linesBandTraceLinkageBtn ("🔗 linkage") in lines header
  - _LASSO_LINKAGE_DEFAULT_PURITY_THRESHOLD = 0.7
  - _LASSO_LINKAGE_DEFAULT_MIN_BAND_SIZE    = 5
  - _LASSO_LINKAGE_MODAL_ID = 'lassoLinkageModal'
  - _computeLassoLinkage(fishSet, candidateList, opts)
  - _lassoLinkageCacheKey(chromFilter, fishSet, candidateList)
  - _lassoLinkageGetOrCompute(opts)
  - _invalidateLassoLinkageCache()
  - _lassoLinkageToTSV(result, opts)
  - _lassoLinkageDownloadTSV()
  - _openLassoLinkagePopover()
  - _closeLassoLinkagePopover()
  - _renderLassoLinkageTable()
  - 12 window exports
  - linkage button click handler in the wiring block
  - setBandTraceFishSet (turn 161) gains an _invalidateLassoLinkageCache
    call so fish-set changes bust the linkage cache too

tests/test_turn164_lasso_linkage.js                 new (118 assertions)
```

No existing functions semantically modified. `setBandTraceFishSet`
gains one new tail call wrapped in try/catch — strictly additive.

---

## 4. Test results

**Single test**: 118 / 0 across 14 sections:

1. Source-pattern checks — 9 function defs, 4 window exports, 4 DOM
   bits, 3 constants, 1 setBandTraceFishSet hook check — 21
2. Clean lasso (purity 1.0 at source band, sort-on-tie) — 13
3. Split lasso (purity < 1.0, threshold gating) — 9
4. Strong-link threshold + min_band_size override — 5
5. chrom_filter — 4
6. Defensive (null/empty/unconfirmed/no-locked_labels/out-of-range) — 9
7. Cache key (order-insensitivity, sensitivity to inputs) — 6
8. `_lassoLinkageGetOrCompute` (cache hit, recompute on change, null paths) — 7
9. `_invalidateLassoLinkageCache` — 2
10. `_lassoLinkageToTSV` (header, sort, NaN-safe, null) — 17
11. `_lassoLinkageDownloadTSV` (filename, null path) — 3
12. setBandTraceFishSet hook check — 1
13. Lines header (button id, label, title, click handler) — 4
14. Regression — turns 163/162/161/160/130/156/122 still wired — 16

**Full sweep**: **3529 / 0** across all `test_turn*.js` files
(was 3411 / 0 at turn 163 close). Zero regressions.

JS-brace balance: clean (12,608 / 12,608). Largest script block
parses under `node --check` with no syntax errors.

---

## 5. What Quentin should exercise

Once a band-trace is active (workflow from turns 161–163):

1. Run a trace as usual: focus a candidate, optionally pick a band,
   click `🔍 trace`. The lines-panel strip lights up with the
   band-trace.
2. **(NEW)** Click `🔗 linkage`. A modal opens.
3. The table shows every confirmed candidate (across all chromosomes
   in `state.candidateList`) and which of its bands the trace's
   fish-set falls into. Strong links (purity ≥ 0.7, n ≥ 5) are at
   the top with a green-tinted background; weak links follow,
   dimmed.
4. Read the patterns:
   - **The focal candidate itself** should be at the top with
     `purity = 100%` and `band = <whatever band you picked>`.
     Sanity check — if it isn't, the linkage compute didn't see
     the focal candidate.
   - **Other strong links** = candidates where the same fish
     predominantly land in one band. Each one is an inheritance
     signal between the focal candidate and that other candidate.
     **This is the manuscript table.**
   - **Weak links** = candidates where the fish split across bands.
     Either the focal candidate's fish-set isn't structured at
     this other candidate, or the other candidate's K-means
     partition is noisy.
5. Click `📊 TSV` inside the modal to download the table.
   Filename: `lasso_linkage_<chrom>_n<n_fish>.tsv`. The first 4
   lines are `#`-prefixed comments documenting the parameters.

**Specific patterns to look for on real LG28:**

- The focal candidate should always be a strong link. Verify by
  pinning a candidate, picking its largest band, tracing → opening
  linkage. Top row should show that candidate at 100% purity.
- A trace with the small HOM_INV band (b2 of K=3, often n=3) might
  produce zero strong links because n_in_best_band < 5 by
  construction. That's the right behaviour — the threshold prevents
  spurious matches on tiny bands.
- If the same fish-set produces strong links at multiple distant
  candidates (e.g. I1 on LG28 ↔ I7 on LG12), that's a **cross-
  chromosome linkage** signal. The current implementation doesn't
  highlight cross-chromosome differently — every confirmed
  candidate, regardless of chromosome, shows up in the table.

---

## 6. What's NEXT

In priority order:

1. **Run on real LG28** — Quentin's review of the linkage table on
   actual data. This is the same blocker as the band-trace turns:
   we need his eyes before doing anything that crosses into
   interpretation.
2. **Slice 2 (alpha intervals on per-sample-lines)** — visual
   companion to the table. ~80 LOC + design conversation about
   Z-order with existing strips.
3. **Click-to-jump on table rows** (~40 LOC). When the user clicks
   a row, if the candidate is on the current chrom, scroll the
   lines panel to its bp range and flash a highlight; otherwise
   show a "switch to LG12" affordance. Needs the cross-chrom
   navigation primitive.
4. **Filter UI in the modal** — slider for purity threshold, input
   for min_band_size. Trivial UX (~30 LOC) but only worth it after
   Quentin says the defaults aren't right for hatchery data.
5. **Lasso surface bridge** — when a lasso surface (PC1 lines or
   tracked PCA) ships, widen the fish-set source from
   `bandTraceFishSet` to `lassoSelection || bandTraceFishSet`.
6. **Combinatorial enumeration** (band-trace turn 161 §6 item 3)
   — its own ~150 LOC turn.
7. **Run brackets** (band-trace turn 161 §6 item 1) — gated on
   Quentin's real-data calibration.

After items 2 + 3 land, the spec moves from `specs_todo/` to a
partially-shipped state alongside the band-trace work; it can fully
graduate once items 4 + 5 (and the enlarged page) ship — those are
multi-turn scopes.

---

## 7. Honest framing

**What's solid:**
- Pure compute is testable in isolation — takes both inputs
  explicitly. 30+ behavioural assertions cover clean, split,
  threshold, chrom_filter, defensive, and out-of-range paths.
- Cache layer integrates with the existing band-trace fish-set
  setter so fish-set changes automatically bust the linkage cache.
  Test verifies the source-level hook + a behavioural recompute.
- The TSV format is self-describing (4 `#`-comment lines explain
  the parameters) and sorts in a deterministic, predictable order
  (strong → weak, then purity desc, then by id asc). Quentin can
  diff the TSVs across runs without spurious row reorderings.
- Headless-safe download path follows the same try/catch pattern as
  turns 162 + 163 — returns the filename in test/headless
  environments, never crashes the panel on Blob/anchor failures.
- Modal HTML structure mirrors the V-shape diagnostic (turn 156),
  so the UX is consistent: close-on-Esc, close-on-backdrop, close
  button with `×`.

**What's risky:**
- The linkage compute walks `state.candidateList` linearly. With 30
  candidates × 226 fish, that's ~6800 operations per call — trivial.
  At 500 candidates (after L2-sweep auto-promote ships fully) it's
  113k ops, still trivial. But the cache key fingerprints **every
  confirmed candidate's id + K + locked_labels.length** which is
  O(n_candidates × avg_id_length). For 500 candidates with 12-char
  ids, that's 6000 char-iterations per cache-key check — still fast,
  but it's the only piece that scales linearly with candidate count.
  If this becomes hot, the right fix is a candidate-list version
  counter.
- The cache key is **order-sensitive on the candidate list** because
  the fingerprint walks in array order. If two passes produce the
  same set of confirmed candidates in different orders (e.g. after
  a sort), the key changes and we recompute. The result is still
  correct; it's just a wasted recompute. A sort-then-fingerprint
  variant would fix this; left as-is because candidate-list reorder
  is rare in practice (only happens during explicit re-sort).
- The modal opens **even when no fish-set is active**, showing an
  instructional empty state. This is friendlier than failing
  silently but it does mean the user can repeatedly click 🔗 linkage
  and get an empty modal each time. Could gate the button visually
  (grey out when `state.bandTraceFishSet` is null) — deferred until
  Quentin reports it as confusing.
- The TSV's `is_strong_link` column is a single-bit flag at the
  default thresholds (0.7 / 5). If Quentin recomputes the table at
  a non-default threshold via the console, the TSV reports the
  thresholds in the comment header but the `is_strong_link` column
  reflects the override. Reading the file requires reading the
  header — clean enough for a supplementary table, but worth
  flagging.
- I didn't add a click-to-jump from the modal table rows to the
  candidate on the chromosome strip (spec §1.4). Right call for a
  one-turn scope, but worth flagging that the modal is read-only
  for now — Quentin can copy a candidate id and search for it but
  the table rows aren't actionable.
- The linkage button click always opens the modal, even on a
  rapid-double-click. The modal display logic is idempotent
  (subsequent opens just refresh the table) so this is fine, but
  if a redraw is racing with a state mutation, the table could
  show stale numbers briefly. The cache invalidation on
  `setBandTraceFishSet` handles the common case.

**What's queued:**
- Quentin reviews the linkage table on real LG28 → confirms the
  defaults work
- Slice 2 (alpha intervals) — visual design conversation
- Click-to-jump from table rows — cross-page navigation primitive
- Lasso surface bridge — when a lasso ships
- Band-trace turn 161 §6 items 1, 3 — gated on real-data
  calibration

End of handoff.
