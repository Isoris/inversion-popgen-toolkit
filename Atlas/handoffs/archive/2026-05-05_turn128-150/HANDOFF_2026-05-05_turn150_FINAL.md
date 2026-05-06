# HANDOFF — turn 150 — Slab-aware U/V cluster modes (closes turn 149's flagged gap)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,090 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 149 handoff.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok).

---

## 1. What this turn ships

Turn 149 closed the slab recluster-dropdown gap by routing kmeans-K3 /
kmeans-K6 to actual slab compare and surfacing a "↩ slab fallback" notice
for the five U/V modes (which had no slab-aware implementation). Turn
150 closes that flagged gap: the U/V modes now run their real algorithms
on slab ranges and the fallback notice is now rare (only fires on
genuine cluster failures).

### 1.1 Architecture — extracted shared cores

The L2 U/V clusterers all followed the same shape:

1. Get rotation cache via `_getOrComputeUVRotation(l2idx)`
2. Use `rot.us` / `rot.vs` / centroids / `baseLabels`
3. Run mode-specific clustering on those
4. Wrap via `_wrapKmeansResultAsCluster(result, l2idx, K, reason)`

Step 1 was the only L2-bound piece (it read `env._s0..env._e0` from the
L2 envelope object). Step 4's `l2idx` parameter was never actually used
inside `_wrapKmeansResultAsCluster` — vestigial.

The refactor splits the rotation pipeline cleanly:

- **`_aggregateWindowRangeForUV(s, e)`** — phase 1, range-keyed. Per-sample
  mean PC1×sign / PC2 across windows. Always uses MEAN regardless of
  `state.aggMethod` (matches the L2 convention; the rotation math doesn't
  need the median variant).
- **`_computeUVRotationCore(xs, ys, nS)`** — phase 2, range-agnostic.
  K-means3 on aggregated (xs, ys) → axis from centroid[0]→centroid[2] →
  rotate to (u, v). Returns the same shape the L2 cache used to return
  directly.
- **`_clusterFromRotation_UV{Rotated,Denoise,DBSCAN,DistRank,DistFuzzy}(rot)`**
  — phase 3, range-agnostic. Each takes a rotation result; runs mode-
  specific clustering; wraps via `_wrapKmeansResultAsCluster(..., null, K, reason)`.

Both L2 and slab paths plug into the same phase 2 and phase 3 cores.
Phase 1 has two callers — L2 uses `env._s0, env._e0`; slab uses the
slab range `s, e` directly.

### 1.2 Slab-side wiring

- **`_getOrComputeUVRotationSlab(s, e)`** — slab variant of the rotation
  cache. Uses `state.slabUVRotationCache` (separate from
  `state.l2UVRotationCache`). Cache key `${s}_${e}`. Same data-identity
  invalidation as the L2 cache. Defends against bad ranges (`BAD_RANGE`
  reason for `s < 0`, `e >= n_windows`, `s > e`).
- **`clusterSlab_UV{Rotated,Denoise,DBSCAN,DistRank,DistFuzzy}(s, e)`** —
  five thin wrappers, each one-liner over the appropriate `_clusterFromRotation_UV*`
  core fed by `_getOrComputeUVRotationSlab(s, e)`.
- **`getSlabClusterByMode(s, e, mode)`** — dispatcher mirroring
  `getL2ClusterByMode`. Routes:
  - `kmeans-K3` (or null/undefined) → `getSlabClusterAt(s, e, state.k)`
  - `kmeans-K6` → `getSlabClusterAt(s, e, 6)`
  - `uv-rotated` / `distance-uv` (alias) → `clusterSlab_UVRotated`
  - `uv-denoise` → `clusterSlab_UVDenoise`
  - `uv-dbscan` → `clusterSlab_UVDBSCAN`
  - `uv-dist-rank` → `clusterSlab_UVDistRank`
  - `uv-dist-fuzzy` → `clusterSlab_UVDistFuzzy`
  - default (unknown future modes like `gmm-3`) → fall back to
    `getSlabClusterAt(s, e, state.k)` (defensive against stale
    localStorage values for retired modes)
  - cache: `state.slabGroupCacheByMode`, key `${s}_${e}_${mode}`,
    same data-identity invalidation as `state.slabGroupCache`.
- **`compareSlabPair_byMode`** — rewritten. The structure now mirrors
  `compareL2Pair_byMode` more closely:
  - K derived once at the top (`kmeans-K6` → 6, else `state.k`)
  - kmeans-K3 / kmeans-K6 / null → short-circuit via `compareSlabPair`
    with `fellBack: false` always
  - U/V modes → fetch via `getSlabClusterByMode`, build the
    contingency directly with `alignLabels` + `chiSquare`/`fisher2x2`,
    return `isSlabPair: true` with `reclusterMode: m`, `fellBack: false`
  - genuine cluster failure (e.g. KMEANS_FAILED in rotation) →
    fall back to kmeans-K3 with `fellBack: true` so the renderer's
    notice still surfaces

### 1.3 L2 path: refactored, behavior preserved

- **`_getOrComputeUVRotation(l2idx)`** — refactored: now calls
  `_aggregateWindowRangeForUV(env._s0, env._e0)` then
  `_computeUVRotationCore(agg.xs, agg.ys, d.n_samples)`. Same cache,
  same return shape.
- **`clusterL2_UV{Rotated,Denoise,DBSCAN,DistRank,DistFuzzy}(l2idx)`** —
  now thin wrappers (2 lines each). Get rotation, call shared core.

This is the one risk surface: L2 mode functionality changed at the
implementation level. Tests confirm same return shapes for clean and
edge cases (degenerate axis, KMEANS_FAILED rotation, etc.). No
adjacent test regressed.

---

## 2. What's now actually working (vs. before)

**Before turn 150 (turn 149 state):**

User picks compare unit `1w` / `5w` / `10w` / `Nw` AND picks recluster
mode `uv-rotated` (or any U/V mode). Result: contingency tables built
with kmeans-K3, accent-bordered notice "↩ slab fallback: uv-rotated
not available for slabs · using kmeans-K3" prepended above each
contingency.

**After turn 150:**

Same user action. Result: contingency tables built with the actual U/V
algorithm — slab rotation cache → mode-specific clustering → real
contingency. No fallback notice. The notice machinery from turn 149
remains in place but only fires on genuine cluster failures (e.g.
KMEANS_FAILED on degenerate rotation), which is what it was designed
to surface anyway.

---

## 3. Files changed / added

| File | Change |
|---|---|
| `Inversion_atlas.html` | +239 LOC. Functions touched/added: `_aggregateWindowRangeForUV` (NEW), `_computeUVRotationCore` (NEW), `_getOrComputeUVRotation` (refactored to use helpers), `_getOrComputeUVRotationSlab` (NEW), `_clusterFromRotation_UV{Rotated,Denoise,DBSCAN,DistRank,DistFuzzy}` (NEW × 5), `clusterL2_UV*` (refactored × 5), `clusterSlab_UV*` (NEW × 5), `getSlabClusterByMode` (NEW), `compareSlabPair_byMode` (rewritten). |
| `tests/test_turn150_slab_uv_cluster_modes.js` | NEW — 109 / 0 |
| `tests/test_turn149_slab_recluster_dispatch.js` | UPDATED — 71 / 0 (was 69 / 0). Two assertions inverted to test the new behavior: U/V `fellBack:true` → `fellBack:false`. Source-pattern checks for the rewritten function structure. |

The turn 149 update follows the same pattern as turn 147b's update of
turn 128 — when the underlying behavior is corrected, the test that
codified the bug must be inverted to assert the fix.

---

## 4. Tests

### 4.1 New suite

**`tests/test_turn150_slab_uv_cluster_modes.js` — 109 / 0** across 14 sections:

1. Aggregation + rotation core extracted (8) — `_aggregateWindowRangeForUV`
   uses MEAN unconditionally, `_computeUVRotationCore` returns ok:true
   with us/vs/angle/centroids, fails cleanly on KMEANS_FAILED, uses
   cx[2]−cx[0] / cy[2]−cy[0] for the Hom1↔Hom2 axis.
2. L2 rotation refactor (4) — uses extracted helpers, no inline
   aggregation loop remaining, still uses `state.l2UVRotationCache`.
3. Slab rotation cache (7) — own cache, key `${s}_${e}`, defends
   against bad ranges, data-identity invalidation.
4. Slab U/V clusterers (15) — five `clusterSlab_UV*` wrappers each
   call `_getOrComputeUVRotationSlab(s, e)` then `_clusterFromRotation_UV*(rot)`.
5. Shared post-rotation cores (6) — five `_clusterFromRotation_UV*`
   declarations + multiline-tolerant count of `_wrapKmeansResultAsCluster(... null, 3, ...)`
   calls (≥5 expected).
6. L2 clusterers preserved (10) — five `clusterL2_UV*` still take
   l2idx, now use the same shared cores via `_getOrComputeUVRotation`.
7. `getSlabClusterByMode` dispatcher (12) — routes all 8 mode names
   correctly (incl. `distance-uv` alias), uses
   `state.slabGroupCacheByMode`, default-case fallback to kmeans-K3.
8. `compareSlabPair_byMode` rewritten (6) — U/V branch fetches via
   `getSlabClusterByMode`, builds contingency via `alignLabels` + `chiSquare`,
   verdict from `mergeThr` / LOW_POWER, sets `isSlabPair: true`,
   handles cluster failure with fellBack:true.
9. **Sandboxed end-to-end rotation core (10)** — synthetic 9-sample 3-band
   data; rotation succeeds with angle ≈ 0° on x-axis-aligned data,
   rotation succeeds with angle ≈ 45° on diagonal data; post-rotation
   axis is u-axis with vs collapsed near zero.
10. **Sandboxed `_clusterFromRotation_UVRotated` (5)** — produces ok:true
    K=3 cluster; degenerate path falls back to `baseLabels` with
    correct reason; rot.ok:false propagates correctly.
11. **Sandboxed `getSlabClusterByMode` dispatch (10)** — routes each of
    8 modes to the correct underlying function; cache hit verified
    (repeat call doesn't re-invoke clusterer); fresh range → fresh call.
12. compareSlabPair_byMode no longer hardcodes fellBack:true (1).
13. Slab renderer fallback notice still wired (2).
14. L2 path preserved (7) — `getL2ClusterByMode`, `compareL2Pair_byMode`,
    five `clusterL2_UV*` wrappers all still take l2idx.

### 4.2 Turn 149 test updated

Source-pattern checks rewritten for the new `compareSlabPair_byMode`
shape (single K derivation, kmeans short-circuit, U/V real-dispatch).
Two sandbox assertions inverted (`fellBack:true` → `fellBack:false`
for U/V modes). New assertion: U/V cluster failure → `fellBack:true`
(genuine fallback path).

### 4.3 Adjacent suites unchanged

- turn 147 L3 toolbar disclosure: 66 / 0
- turn 147b windows-N engagement: 26 / 0
- turn 148 slab L3 parity: 62 / 0

### 4.4 Full sweep

**Across parseable turn-numbered tests: 2504 / 0** (was 2393 baseline,
+111 net: +109 turn 150 + 2 from turn 149 updates).

JS syntax check on the 2.78MB extracted main script block: clean.

---

## 5. Atlas state

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 147)   | 71,620  | 2262            | 51    |
| Post-turn-148            | 71,767  | 2324            | 52    |
| Post-turn-149            | 71,851  | 2393            | 53    |
| Post-turn-150 (this)     | 72,090  | 2504            | 54    |
| Δ session (3 turns)      | +470    | +242            | +3    |

---

## 6. What's still queued

From the turn 147 → 149 sequence, the queue is now:

4. **G-popup → popgen page merge** (~400 LOC) — fold the 3-tab modal
   (karyotype / inheritance / manual) into popgen page (`page12`) as
   a second tab strip with different background shade. Per Quentin's
   brief: *"merge into one page with 2 rows of tabs (different bg
   shade)"*. **Highest-LOC remaining UI fix; would benefit from a
   sketch from Quentin to ground the desired layout.**
5. **Tutorial authoring** — eight `data-status="pending"` cards.
   Launcher infrastructure isn't yet built either; needs ~110 LOC of
   overlay launcher + ~80 LOC content per tutorial. Visual assets
   (screenshots / GIFs) would significantly improve effectiveness;
   text-only tutorials are still useful but not as good as the
   *"30 seconds to figure 3"* pattern they're meant to mirror.
6. ~~**Slab-aware U/V cluster modes**~~ — **CLOSED THIS TURN.**

New from this turn:

7. **Cross-pane spotlight in U/V mode** — pre-existing inconsistency
   that affects both L2 and slab modes equally. When the user picks
   a U/V recluster mode, the contingency table cells are computed
   against U/V labels (via `getL2ClusterByMode` / `getSlabClusterByMode`),
   but `_applySpotlightHighlights` looks up cluster labels via
   `getL2Cluster` / `getSlabClusterAt` (kmeans-K3 / -K6 only). So the
   cell where a spotlighted sample appears highlighted may not match
   the cell where its sample-data row actually lives in the U/V
   contingency. Estimated ~30 LOC: pass `meta.reclusterMode` through
   `__l3_meta` and route the cluster lookup through the by-mode
   dispatcher in `_applySpotlightHighlights`. **Slab and L2 should
   be fixed together for consistency.**

---

## 7. Honest framing

**What turn 150 actually delivered:**

- The U/V degradation notice from turn 149 is now a rare-edge
  diagnostic instead of the every-day reality for slab-mode U/V users.
  All five U/V modes (uv-rotated, uv-denoise, uv-dbscan, uv-dist-rank,
  uv-dist-fuzzy, plus the distance-uv alias) compute real
  contingency tables on slab ranges.
- The L2 implementation was refactored to share the same cores. No
  code duplication: 5 cores serve both L2 and slab. Bug fixes to U/V
  algorithms now propagate to both modes automatically.
- Sandboxed end-to-end test of the rotation pipeline against
  synthetic 3-band data (9 samples, 3 clusters at known centroids),
  including the 45° off-axis case that exercises the rotation step.
  Rotation is verified to align Hom1↔Hom2 with the u-axis and collapse
  perpendicular variation into vs.
- Fail-soft maintained: bad ranges, missing helpers, and KMEANS_FAILED
  rotations all return early without crash. The fallback path from
  turn 149 still works for these genuine failures.
- L2 path explicitly preserved: 5 `clusterL2_UV*` wrappers still take
  l2idx; `getL2ClusterByMode` and `compareL2Pair_byMode` unchanged;
  10 negative tests confirm no slab-mode contamination of L2 logic.

**What it didn't deliver:**

- Spotlight highlighting in U/V mode (item #7 above) — pre-existing
  cross-mode inconsistency not introduced this turn but flagged here
  for completeness.
- An actual end-to-end slab-render integration test against a real
  cohort fixture — the test harness is structurally source-pattern +
  sandboxed-unit, which doesn't trivially extend to the full DOM
  render path. Would require a fixture cohort + DOM emulation; out
  of scope for this turn.

**Manuscript impact:** moderate. U/V modes are research/diagnostic
tools — they're how Quentin stress-tests merge calls when the standard
K-means K=3 partition looks suspicious (e.g. when family clusters
are dragging centroids, when DBSCAN-within-stripe finds hidden
substructure, when fuzzy-distance reveals ambiguous samples). Turn 149
made these tools visible in slab mode but non-functional; turn 150
makes them real. Slab compare units (1w / 5w / 10w / Nw) are
particularly useful at boundary refinement, where Quentin needs the
finest resolution his data can support, and that's exactly where
U/V's noise-handling matters most.

**Bundle**:
`Atlas_full_bundle_2026-05-05_turn150_slab_uv_cluster_modes.tar.gz`.

Walk the map carefully, respect cohort discipline, don't break the
test suite.
