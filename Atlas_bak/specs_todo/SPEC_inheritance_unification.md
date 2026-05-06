# SPEC — Inheritance system unification (replace candidate-only with universal-input)

**Status**: drafted turn 130. **Decision-pending**, NOT to be implemented
until L2-sweep + sliding-window + fish-trajectory have all shipped and
been used on real data. Estimated ~1.5 turns IF the decision goes this
way.

**Trigger** (Quentin, turn 130):
> *"Replace the candidate-based version with the L2-sweep. Stop requiring
> user-promoted candidates entirely."*

This is option 3 from the turn-130 menu. It's a refactor, not a feature.
The point of writing the spec now is to record the **decision criterion**
so future-Claude (or current-Quentin) knows when to pull the trigger and
when not to.

---

## 1. The current architecture

`runInheritanceCompute()` calls `_gatherActiveCandidatesForInheritance()`
which reads from `state.candidates` (now correctly bridged from
`state.candidateList` thanks to turn 129's `_rebuildCandidateRegistries`).
The function returns a list of items that have:

- `id` — the candidate id
- `labels` — Hungarian-aligned-to-itself K-means labels (locked at promote)
- `K` — number of bands
- `start_bp` / `end_bp`
- `seq_num` — sequential index by start_bp on this chromosome

`inheritanceGroupClustering(items)` then runs Jaccard / agglomerative on
the (item, band) pairs.

## 2. What unification would mean

Replace `_gatherActiveCandidatesForInheritance()` with a dispatcher:

```
function _gatherInheritanceItems() {
  const mode = state.inheritanceItemSource || 'candidates';
  switch (mode) {
    case 'candidates':       return _gatherActiveCandidatesForInheritance();
    case 'l2_sweep':         return _gatherL2SweepItems();
    case 'sliding_window':   return _gatherSlidingTilesItems();
    case 'union':            return _gatherUnionItems();         // candidates + sweep + tiles, dedupe by bp overlap
  }
}
```

Plus `state.inheritanceItemSource` switching all downstream consumers to
read from the same gatherer. The I·g pills, tooltip, suggestions, and
matrix viewer all get the unified result.

## 3. Why we should NOT do this yet

Three risks:

1. **Schema drift in items**: candidate items have a lot of metadata
   (locked_labels_raw, vocab, classifier outputs, parent_l2, fish_calls
   per-sample). L2-sweep and sliding-window items don't. Code that
   assumes `cand.fish_calls` exists will break. Audit needed.

2. **Cache key explosion**: `_inheritanceCacheKey` currently includes
   per-candidate label fingerprints. If items can come from 3 sources,
   the cache key needs to incorporate the source. Otherwise the cache
   gets stale or wrongly-shared.

3. **UI semantics shift**: today the I·g pills mean "user promoted these
   N candidates and the inheritance compute clustered their bands". If
   pills can come from auto-sweep, the user's mental model changes —
   they might trust the pills less if they're auto-generated. The
   `🤖` badge in the L2-sweep spec signals this, but unification would
   blur the line.

## 4. When we SHOULD do this

If, after L2-sweep + sliding-window + fish-trajectory have shipped:

- Quentin reports the candidate-based version is **redundant** with the
  sweep results (the sweep auto-promotes basically everything that
  matters).
- AND the sweep + tile + trajectory results agree more than they
  disagree on a typical chromosome.
- AND the manuscript's Figure / Methods would be cleaner with one
  unified inheritance system.

Then unification is a refactor turn that simplifies the code.

If, instead, Quentin reports the systems give **different** information
that he wants both of:
- candidate-based: "what I confirmed manually"
- sweep / tile / trajectory: "what the algorithm proposes"

Then we keep them parallel forever. No unification.

## 5. Decision criterion

Run for one full chromosome (LG28 — already the validation chromosome)
through every system after Slices 1 land:

1. `runInheritanceCompute()` on confirmed candidates only — record N
   inheritance groups, group composition.
2. `runL2SweepInheritance()` — record N groups + which L2s auto-promoted.
3. `runSlidingWindowInheritance(10)` — record N groups + tile coverage.
4. `runLineageCompute()` — record N lineages per fish + dendrogram.

Compare:
- How many fish change lineage / inheritance group between systems?
- Do the auto-promoted candidates match the user's manual promotions?
- Are the I·g pills the same labels?

If agreement ≥ 80% across systems on LG28 → unify (this spec).
If agreement < 80% → keep parallel (don't ship this spec).

Let the data decide. Don't unify speculatively.

## 6. Implementation sketch (only if decision is GO)

### Slice 1 — dispatcher (~0.5 turn)
- [ ] `_gatherInheritanceItems(mode)` dispatcher
- [ ] `state.inheritanceItemSource` slot + localStorage
- [ ] Audit every downstream consumer for `cand.*` field access; replace
      with safe accessors that handle missing fields.

### Slice 2 — UI picker (~0.5 turn)
- [ ] Dropdown in the lines-panel header: "items: [candidates | sweep |
      tiles | union]"
- [ ] Force-recompute on switch.

### Slice 3 — cache key + invalidation (~0.5 turn)
- [ ] Add source to `_inheritanceCacheKey`
- [ ] Invalidate when source changes
- [ ] Per-source dismiss sets

## 7. What this is NOT

- **Not deletion of `state.candidateList`.** That's the canonical user
  curation, kept regardless of whether it feeds inheritance.
- **Not deletion of the parallel registry bridge.** The bridge stays;
  it just becomes one of multiple input sources.
- **Not a turn we should burn before LG28 measurement.** Premature
  unification risks re-doing it the wrong way.

## 8. Tests (only if shipped)

- Each source independently produces valid items.
- Switching source invalidates cache + re-renders pills.
- Union source dedupes overlapping items by bp range.
- Backward compatibility: existing tests for `runInheritanceCompute()`
  with default source ('candidates') still pass.
