# SPEC — L2-sweep inheritance mode

**Status**: Slices 1 + 2 SHIPPED — auto-promote in turn 133
(`tests/test_turn133_l2_sweep_auto_promote.js`), inspector in turn 134
(`tests/test_turn134_l2_sweep_inspector.js`). The dispatcher refactor that
makes inheritance compute consume L2 envelopes directly (rather than just
auto-promoting them into the candidate list) is still pending — see also
`SPEC_inheritance_unification.md` for the broader decision context. Drafted
turn 130.

**Trigger** (Quentin, turn 130):
> *"It could auto promote candidates (i feel its right)."*
>
> *"L2-sweep mode: same Jaccard/agglomerative pipeline, but inputs are ALL
> L2 envelopes on the chrom (not just user-promoted candidates)."*

The current `inheritanceGroupClustering` (turn 117) only sees candidates
the user has explicitly promoted. Most users won't promote candidates
densely enough to populate the inheritance matrix on a chromosome with
many small inversions. The L2-sweep mode runs the same algorithm with a
different input set: **every L2 envelope on the chromosome** treated as
an implicit candidate.

---

## 1. Why this is needed

Three failure modes the current candidate-based system has:

1. **Empty matrix**: chromosome has 18 L2s but the user only promoted 2
   candidates → only 2 items reach the inheritance compute → nothing
   useful comes out.
2. **Selection bias**: users tend to promote the obvious candidates and
   miss the weaker but real ones. The matrix then reflects what was
   visually salient, not what's structurally informative.
3. **Cold start**: every new chromosome needs the user to do at least
   3–4 manual promotions before the I·g pills become useful. High
   activation energy.

L2-sweep removes all three by running the inheritance algorithm
unconditionally on the full L2 inventory.

## 2. Algorithm

Identical to `inheritanceGroupClustering` (the per-(item, band) Jaccard
clustering), but the input items are L2 envelopes instead of candidates:

```
items = state.data.l2_envelopes
    .filter(env => isUsableL2(env))
    .map(env => ({
      id: `L2:${env._idx}`,
      labels: getL2Cluster(env._idx).fixedKLabels,   // n_samples Int8Array
      K: state.k || 3,
      start_bp: env.start_bp,
      end_bp: env.end_bp,
      meta: { source: 'l2_sweep', l2_idx: env._idx },
    }));
runInheritanceCompute({ force: true, items });
```

`isUsableL2` filters out L2s where K-means failed (silhouette below
threshold, fewer than `_IGC_MIN_BANDS_FOR_CLUSTERING` populated bands,
etc.) — they would just inject noise into the Jaccard distance matrix.

## 3. Auto-promote candidates from sweep results

This is Quentin's turn-130 addition. After the L2-sweep computes
inheritance groups, the atlas can auto-promote L2s that look like real
candidates by surfacing them in `state.candidateList` with
`confirmed: false` and a tag `auto_promoted_from: 'l2_sweep'`.

### Selection rule (configurable)

An L2 auto-promotes when ALL of:

1. The L2's K-means silhouette ≥ `_AUTO_PROMOTE_MIN_SILHOUETTE`
   (default 0.30 — a calibration knob).
2. At least 2 inheritance groups are populated by this L2's bands
   (otherwise this L2 is "all one ancestry" and not informative).
3. The L2's smallest band has ≥ `_IGC_MIN_BAND_SIZE` fish (default 5).
4. The L2 has not already been confirmed as a candidate by the user.
5. The L2 is not within `_AUTO_PROMOTE_DEDUPE_BP` of an already-saved
   candidate (default 100kb — avoids cluttering the catalogue with
   adjacent L2s that the user already promoted manually).

### Auto-promote mechanics

- Tag: `cand.source = 'auto_l2_sweep'`, `cand.auto_promoted_at = ISO date`,
  `cand.confirmed = false`.
- Visible in catalogue with a distinct chip (e.g., 🤖 robot icon).
- User can confirm (promotes to standard candidate, clears the auto tag)
  or dismiss (drops from list, adds to a per-chrom `dismissedAutoIds`
  set in localStorage so the same L2 isn't re-auto-promoted on next
  sweep).
- `persistCandidateList()` after the sweep so localStorage round-trip
  works.

### Why this is OK to do automatically

The risk of auto-promote is creating spurious candidates the user has
to manually clean up. The five gates above (silhouette + ≥2 groups +
band size + dedupe + not-dismissed) keep the false-positive rate low.
And the user always retains the dismiss override — auto-promotes are
suggestions, not commitments.

The reward is large: the inheritance matrix becomes useful from the
first chromosome load, not after 4–5 manual promotions.

## 4. UI surfacing

### 4.1 Toggle in the L3 toolbar

Checkbox: `[ ] L2-sweep`. Default off (so power users can keep the
manual workflow). When on:

- On chromosome load (`applyData` callback), trigger
  `runL2SweepInheritance()` after a short debounce.
- Auto-promotes propagate to the candidate strip + catalogue.
- I·g pills draw on the per-sample-lines panel reflecting the broader
  group structure.
- A small `🤖 swept` badge appears on the lines header.

State slot: `state.l2SweepEnabled` (boolean, default false). localStorage:
`pca_scrubber_v3.l2SweepEnabled`.

### 4.2 Sweep-result inspector

Optional Slice 2: a modal that lists every L2 in the sweep with its
inheritance group_id, silhouette, band sizes, and a single-click
"promote" / "dismiss" pair of buttons. Lets the user quickly triage all
L2s on the chromosome.

## 5. Implementation slices

### Slice 1 — sweep + auto-promote (~1 turn)
- [ ] `runL2SweepInheritance()` orchestrator
- [ ] `isUsableL2(env)` filter
- [ ] `_autoPromoteFromSweep(result)` honoring the 5 gates + dismiss set
- [ ] localStorage `pca_scrubber_v3.l2SweepDismissed.<chrom>` Set
- [ ] L3 toolbar checkbox + handler
- [ ] Tests against synthetic L2 envelopes with planted inheritance

### Slice 2 — sweep inspector modal (~0.5 turn)
- [ ] Modal listing all swept L2s with promote/dismiss buttons
- [ ] Filter by silhouette / group count
- [ ] Bulk promote / bulk dismiss

## 6. Open design questions

1. **Default ON or default OFF?** Current spec says OFF for caution.
   Could flip to ON once Slice 1 has run on Quentin's real chromosomes
   without nuisance auto-promotes.
2. **Silhouette threshold**: 0.30 is a guess. Calibrate on LG28 + LG12.
3. **Dedupe radius**: 100kb is arbitrary. Chromosome 28's known
   ~2.89 Mb inversion is the test case — the sweep should auto-promote
   ONE candidate covering it, not 5 overlapping L2s.
4. **Detailed mode**: should sweep run on `candidates_detailed` too
   (K=6)? Default: yes, separate cache per mode.

## 7. What this is NOT

- **Not a replacement for manual promotion.** Users still curate.
  Sweep just provides a starting set.
- **Not a recompute on every UI tick.** The sweep result is cached
  per-chromosome; only re-runs on chromosome load or explicit user
  click of a "re-sweep" button.
- **Not the trajectory inheritance system.** That's a separate axis
  (fish-trajectory clustering, see
  `SPEC_distant_band_concordance_fish_trajectory.md`). L2-sweep can
  benefit FROM trajectory lineages (use them to weight Jaccard) but
  doesn't depend on them.

## 8. Tests (Slice 1)

- 6 synthetic L2s with planted inheritance pattern → sweep recovers it.
- Auto-promote gate verification: low-silhouette L2s don't promote;
  L2s with only 1 group don't promote.
- Dedupe verification: 3 adjacent L2s overlapping a single biological
  inversion → only one auto-promote.
- Dismiss persistence: dismiss a sweep result, reload chromosome,
  same L2 not re-promoted.
- localStorage round-trip on `l2SweepEnabled` and `l2SweepDismissed`.
