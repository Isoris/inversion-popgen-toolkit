# SPEC — Interval Collector (passive automatic accumulator) [TIER 2]

**Status**: drafted turn 130 final session. Spec'd in chat `819d8454`
("Live PCA visualization with contingency tables") as a planned "session
3" feature for the scrubber. Not yet implemented in atlas. Estimated
~1 turn.

**Trigger** (Quentin, chat `819d8454`):
> *"We could do something like use the tracked samples to capture based
> on contingency table very dynamically the intervals of boundaries and
> save to export. Its the interval collector page (all automatic) we
> don't need to open it."*

**Tier 2** — atlas/UI engineering, no direct lit anchor. Useful but not
critical-path for the master automatic pipeline.

---

## 1. The idea

A **passive accumulator** that runs whenever tracked samples change and
records candidate intervals derived from pairwise stripe-disagreement
patterns. The user does not have to open a page for collection to
happen — the page is for inspection only.

## 2. Algorithm

For each pair of currently-tracked samples (s_i, s_j):

1. Look up `ghsl_kstripes.by_k["3"].stripe_per_sample[i]` and `[j]`.
2. If samples are in different stripes globally, scan windows: at each
   window w, compute their per-window K-band assignment.
3. Find consecutive runs where the pair agrees (same band) and disagrees
   (different bands). Boundaries are at agree↔disagree transitions.
4. Each transition contributes one boundary signal for this pair.

Aggregate across all tracked pairs:

5. For each window w, count how many pairs have a boundary signal at w
   (or within ±2 windows for smoothing).
6. Find local maxima in this count vector. Each maximum above a
   threshold becomes a captured boundary.
7. Consecutive captured boundaries define intervals — each gets an entry
   in `state.collectedIntervals`.

## 3. State

```
state.collectedIntervals: { [chrom]: [
  { interval_id, start_bp, end_bp, start_w, end_w,
    n_supporting_pairs, mean_strength,
    captured_at, source: "tracked_pair_disagreement",
    derived_from_samples: ["CGA042", "CGA107"],
    K_used: 3 }
] }
```

Persisted to localStorage. Cleared per chrom switch (separate per
chrom; can browse back).

## 4. Surface

A page that shows captured intervals (no compute on page open — collector
runs in background regardless). Inspect mode:
- Sortable table: id, span, n_pairs, mean_strength, captured_at
- Click → navigate to that window range
- Promote → manual confirm into `state.candidateList`

## 5. Auto-promote vs review-and-promote

Two options for review at session-3 time:

- **Option A (auto-promote)**: collector intervals join the candidate
  registry as `discovery_method = "scrubber_collector"` candidates,
  alongside user-promoted L2 candidates.
- **Option B (review-and-promote)**: collector intervals stay in their
  own slot; user reviews and promotes the strong ones.

Default: **Option B** (less intrusive). User can flip a setting to
"auto-promote everything above strength X."

## 6. Implementation slices

### Slice 1 — collector engine + state slot (~0.5 turn)
- `collectorEngine(trackedSamples, ghslPanel, kStripes)` returns array
- Recompute on tracked-samples change
- Store in `state.collectedIntervals`
- Persist to localStorage

### Slice 2 — inspection page (~0.3 turn)
- Sortable table view
- Click handlers
- Promote-to-candidate action

### Slice 3 — registry export integration (~0.2 turn)
- Include collected intervals in candidate JSON export
- Tag as `discovery_method = "scrubber_collector"`

## 7. Open questions

1. **Threshold for "boundary signal"**: how many pairs must transition
   at the same window to count as a captured boundary? Default n=3.
2. **Smoothing window**: ±2 windows is loose enough to absorb noise,
   tight enough to preserve resolution. Configurable.
3. **What about small lassos?** Need ≥3 tracked samples to compute
   pair transitions. Fewer → empty collector.
4. **Cross-chrom collector**: stays per-chrom (different K-band
   assignments per chrom).

## 8. Tests

- Synthetic: 5 tracked samples, planted boundary at window 100 →
  collector captures boundary at window 100 ± 2.
- No tracked samples → empty collector.
- Single tracked sample → empty collector (need pairs).

## 9. Cross-references

- `SPEC_lasso_inheritance_backgrounds.md` — uses lassoed samples for
  per-candidate purity; this collector uses them for boundary signal.
- `SPEC_l2_sweep_inheritance.md` — different auto-detection pathway
  (Jaccard between candidates).
- Chat `819d8454` for the original design.
