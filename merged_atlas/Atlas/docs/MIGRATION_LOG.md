# Migration log

Running record of what's been extracted from the legacy
`Inversion_atlas.html` into the new modular layout. Each entry
records: what moved, where it landed, whether parity was verified
against the legacy binary, and what's still pending in that area.

This log is the single source of truth for "did we ship X yet?"
The full plan is in `docs/MIGRATION_INVENTORY.md`.

---

## Step 0 — foundation (turn — already shipped)

**Status**: ✅ complete

Atlas folder skeleton, canonical layout, shared infrastructure.

- `shared/state_io.js` (v `SHARED_VERSION = 1`) — JSON I/O, layer
  registry, path helpers, State class, headless TSV/JSON serializers.
- `build/flatten.py` — ES-module → single-file inliner with
  multiline-import support.
- `build/package_workflow.py` — release-time per-workflow splitter.
- `data/` — five-folder canonical layout
  (`precomp/cohort/candidates/comparative/review/`).
- `data/README.md` — locked layout doc.
- `inversion_review.html` + `inversion_review/main.js` — smoke-test
  shell validating the modular pattern.
- `tests/test_modular_smoke.js` — **58/58 pass**.

---

## Step 1 — extract `shared/` pure-compute primitives (turn — complete)

**Status**: ✅ complete, **bit-identical parity verified** against
`Inversion_atlas.html` (turn 165 close binary).

### What moved

| Source (legacy line) | New module | Export |
|---|---|---|
| `_buildContingency` (12250) | `shared/contingency.js` | `buildContingency` |
| `_detectFuseEvents` (12274) | `shared/contingency.js` | `detectFuseEvents` |
| `_computeARI` (12320) | `shared/contingency.js` | `computeARI` |
| `_computeNMI` (12353) | `shared/contingency.js` | `computeNMI` |
| `_scaleStabilityVerdict` (12408) | `shared/contingency.js` | `scaleStabilityVerdict` |
| `_cramersV` (36046) | `shared/contingency.js` | `cramersV` |
| `_chiSqSurvival` (38443) | `shared/contingency.js` | `chiSqSurvival` |
| `_lnGamma` (38482) | `shared/contingency.js` | `lnGamma` |
| `alignLabels` (30843) | `shared/hungarian.js` | `alignLabels` |
| `permutations` (30879) | `shared/hungarian.js` | `permutations` |
| `_hungarianChainProjection` (39120) | `shared/hungarian.js` | `hungarianChainProjection` |
| `_finalizeChain` (39190) | `shared/hungarian.js` | (internal) |
| `_concordanceMatrix` (39212) | `shared/hungarian.js` | `concordanceMatrix` |
| `_LINEAGE_CHAIN_BREAK_AGREEMENT` (39090) | `shared/hungarian.js` | `LINEAGE_CHAIN_BREAK_AGREEMENT` |
| `_HET_RAMP` (16569) | `shared/het_rate.js` | `HET_RAMP` |
| `_hetRateColor` (16575) | `shared/het_rate.js` | `hetRateColor` |

Each module preserves the legacy `window._foo = foo` console-debug
exposures inside an `if (typeof window !== 'undefined')` block.

### What did NOT move (deferred to Step 2 — `inversion_discovery` extraction)

These functions appear pure-ish but read from `state.data`,
`state.__hetRateCache`, or other discovery-side state. They'll move
together with their owning sub-atlas:

| Function | Reason | Will land in |
|---|---|---|
| `_computeHetRateForL2` (16735) | reads `state.data.l2_envelopes`, dosage chunk LRU | `inversion_discovery/het_rate_compute.js` |
| `_computeHetRateForRange` | dosage chunk LRU | same |
| `_computeHetRateForSlab` (16752) | reads `state.data.windows.start_bp/end_bp` | same |
| `_getHetRateCache` (16612) | touches `state.__hetRateCache` | same |
| `_invalidateHetRateCache` (16619) | touches `state.__hetRateCache` | same |
| `getL2Cluster` (10673) / `getL2ClusterAt` (10770) / `getL2ClusterByMode` (11656) / `getSlabClusterByMode` (11958) | the K-means impl reads `state.data.tracks` and the per-L2 cache | `inversion_discovery/k_means.js` |
| `recomputeAnchorConcord` (36081) | thin wrapper around `cramersV` that reads `state.trackingAnchor` | `inversion_discovery/anchor_concord.js` |

### Parity verification

`tests/test_legacy_parity.js`:

- Loads the legacy `Inversion_atlas.html` (turn 165 close binary)
- Extracts each target function source by name + brace-walk
  (works around the truncated `<script>` in the upload)
- Evaluates ONLY those functions in a tiny VM context (no DOM stubs needed)
- Runs a 74-case battery comparing legacy `_foo(...)` to new
  `import { foo }`
- **74/74 pass** including:
  - `cramersV` matches to 14 decimal places
  - `hetRateColor(0.42)` matches `rgb(213,224,235)` exactly
  - `concordanceMatrix` Float32Array bit-identical
  - `chiSqSurvival(1000, 1)` ≈ 1e-219 matches exactly
  - `alignLabels` perm/aligned/concord all bit-identical

### Test totals

```
test_modular_smoke.js          58 pass / 0 fail
test_shared_contingency.js     45 pass / 0 fail
test_shared_hungarian.js       40 pass / 0 fail
test_shared_het_rate.js        20 pass / 0 fail
test_legacy_parity.js          74 pass / 0 fail
                              ----
                              237 total / 0 fail
```

---

## Step 2 — `inversion_discovery.html` + `inversion_discovery/` (in progress)

### Step 2.0 — `shared/state.js` (turn — complete)

**Status**: ✅ complete

The canonical state contract is now `shared/state.js`:

- `SLOT_REGISTRY` — typed inventory of state slots with tags
  (`cross_atlas` / `persisted` / `transient` / `derived`).
  Currently 70+ slots covering the most-touched areas of the legacy
  state object (which has 110+ declared + 120+ runtime-created slots).
  Inventory grows as sub-atlases are extracted.
- `makeState({slots})` — factory returning a fresh state object with
  defaults (deep-cloned so two calls don't share refs)
- `serializeState(state, opts)` / `deserializeState(snapshot)` —
  JSON-friendly snapshot of cross-atlas slots only.
  Round-trips Set, Map, Int8Array, Float32Array, etc.
- `mergeStateFromSession(state, snapshot)` — apply a snapshot from a
  sibling sub-atlas onto a fresh state in place
- `readPersistedSlots(storage)` / `writePersistedSlot(name, value)` —
  localStorage hydration with type coercion

**Cross-atlas slot list** (the ones that round-trip):
`candidate`, `candidateList`, `candidate_detailed`, `candidateList_detailed`,
`lockedLabels`, `lockedRefL2`, `manualGroups`, `tracked`, `l2SweepResult`,
`activeSampleSet`, `activeSampleReasons`, `activeSampleRules`,
`bandTraceFishSet`, `candidate_review_decisions`, `locked_karyotype_groups`.

These are the only slots that need to travel between sub-atlases via
`data/review/inversion/sessions/`. Every other slot is local to one
atlas and lives in localStorage or in-memory.

`tests/test_shared_state.js`: **78/78 pass**, covering registry shape,
makeState defaults + deep-clone, curated subset, JSON round-trip for
all typed forms, transient-slot exclusion, mergeStateFromSession,
localStorage round-trip with coercion, future-version handling, and
empty/malformed inputs.

### Step 2 strategy reframing — discovered after auditing

**Original plan**: extract page1 in one or two turns.

**What the audit showed**: page1 reaches **680 of 1064 top-level
functions** in the legacy atlas. Page1 isn't a clean island — opening
a candidate in the scrubber triggers boundary code, classification,
band-trace, lineage compute, etc. A page-by-page extraction would
duplicate ~80% of the function graph across sub-atlases.

**Revised strategy** — extract by **layer**, not by page:

1. **Layer A — pure compute primitives**  ✅ done in Step 1
   (contingency, hungarian, het_rate; 13 primitives, 74 parity tests)
2. **Layer B — state contract**           ✅ done in Step 2.0
   (`shared/state.js`; 78 tests)
3. **Layer C — data loaders + caches**    *next*
   - `shared/per_l2_cluster.js` — `getL2Cluster`, `getL2ClusterAt`,
     `getL2ClusterByMode`, `clusterL2`, `clusterL2AtK`,
     `kmeans1D`, `kmeans2D`, `silhouette1D`, `adaptiveK1D`
   - `shared/het_rate_compute.js` — `_computeHetRateForL2`,
     `_computeHetRateForRange`, `_computeHetRateForSlab`,
     `_getHetRateCache`, `_invalidateHetRateCache`
   - `shared/dosage_chunks.js` — chunk LRU + index lookup
   - `shared/pc_helpers.js` — `getPC`, `getPCByAxis`, `availablePCs`,
     `computePC1Signs`, `_getOrComputeUVRotation`
4. **Layer D — panel renderers** — split per-sub-atlas:
   - `inversion_discovery/sim_mat_panel.js`
   - `inversion_discovery/z_panel.js`
   - `inversion_discovery/lines_panel.js`
   - `inversion_discovery/pca_panel.js`
   - `shared/l3_panel.js` (used by both discovery and review)
5. **Layer E — page glue + DOM wiring** — the thinnest layer per
   sub-atlas, lives in `inversion_discovery/main.js`,
   `inversion_review/main.js`, etc.

This means **Step 2 will be 4-5 sub-steps** (2.0 / 2.A clusters /
2.B het + chunks / 2.C panels / 2.D page glue), but each one ships a
testable layer. The extraction risk is much lower because tests
verify each layer in isolation against the legacy.

### Step 2.A — `shared/kmeans.js` + `shared/per_l2_cluster.js` (turn — complete)

**Status**: ✅ complete, **bit-identical parity verified** for the pure
K-means primitives.

#### What moved

| Source (legacy line) | New module | Export |
|---|---|---|
| `kmeans1D` (10092) | `shared/kmeans.js` | `kmeans1D` |
| `kmeans2D` (10133) | `shared/kmeans.js` | `kmeans2D` |
| `silhouette1D` (10178) | `shared/kmeans.js` | `silhouette1D` |
| `adaptiveK1D` (10224) | `shared/kmeans.js` | `adaptiveK1D` |
| `aggregateL2` (10249) | `shared/per_l2_cluster.js` | `aggregateL2` |
| `sampleSpreadL2` (10294) | `shared/per_l2_cluster.js` | `sampleSpreadL2` |
| `sampleSpreadRange` (10304) | `shared/per_l2_cluster.js` | `sampleSpreadRange` |
| `sigmaProfileL2` (10342) | `shared/per_l2_cluster.js` | `sigmaProfileL2` |
| `clusterL2` (10503) | `shared/per_l2_cluster.js` | `clusterL2` |
| `clusterL2AtK` (10694) | `shared/per_l2_cluster.js` | `clusterL2AtK` |
| `getL2Cluster` (10673) | replaced by `ClusterCache` class | `ClusterCache` |
| `ensureL2Cache` (10665) | replaced by `clusterCacheKey` + `ClusterCache` | `clusterCacheKey` |

#### Refactor pattern: explicit context

The legacy `clusterL2(l2idx)` reads from a global `state` object. The
new API takes an explicit `ctx` with the relevant fields:

```js
import { contextFromState, clusterL2, ClusterCache } from './shared/per_l2_cluster.js';

const ctx = contextFromState(state);   // bundles state.k, state.aggMethod,
                                        // state.kMode, state.kRange, state.silThreshold,
                                        // state.minNGroup, state.minNWin, state.flipPC1,
                                        // state.pc1Sign + a getPC callback

const cl = clusterL2(ctx, l2idx);       // pure call, no global reads

// cache:
const cache = new ClusterCache();
const cl2 = cache.getOrCompute(ctx, l2idx);  // memoized; auto-invalidates on ctx change
```

Same shape as `hungarianChainProjection` from Step 1: pure functions
with caller-supplied accessors. No hidden state, no globals.

#### Tests

```
tests/test_shared_kmeans.js          35 pass / 0 fail   (synthetic fixtures)
tests/test_shared_per_l2_cluster.js  61 pass / 0 fail   (synthetic fixtures)
tests/test_legacy_parity.js          118 pass / 0 fail  (was 74, +44 for kmeans)
```

The kmeans parity battery covers `kmeans1D`, `kmeans2D`, `silhouette1D`,
`adaptiveK1D` against the legacy implementation:
- `silhouette1D` 3-clean-clusters case: `0.9731919434440444` matches to all 16 decimals
- `adaptiveK1D` always picks the same K, returns the same labels and centers
- All edge cases (empty input, K=1, NaN paths, kMin fallback) match exactly

The per_l2_cluster module isn't covered by the parity test (would
require a complete state.data fixture — large to embed). Synthetic
fixtures in test_shared_per_l2_cluster.js verify the math is correct
end-to-end and that the API contracts hold.

#### Test totals

```
test_modular_smoke.js                58 ✓
test_shared_contingency.js           45 ✓
test_shared_hungarian.js             40 ✓
test_shared_het_rate.js              20 ✓
test_shared_state.js                 78 ✓
test_shared_kmeans.js                35 ✓  NEW this turn
test_shared_per_l2_cluster.js        61 ✓  NEW this turn
test_legacy_parity.js               118 ✓  (extended +44 cases)
                                    ----
                                    455 total / 0 fail
```

### Step 2.B — `shared/het_rate_compute.js` + `shared/dosage_chunks.js` (next)

Plan: extract `_computeHetRateForL2`, `_computeHetRateForRange`,
`_computeHetRateForSlab`, `_getHetRateCache`, `_invalidateHetRateCache`
from legacy line 16735+. State-coupled — the cache lives in `state`,
keyed off the dosage_chunks layer load. Refactor: caller passes
`ctx = { data, hetRateCache, dosageChunks: {get, has, ... } }` and
the new module never touches a global.

Also `shared/dosage_chunks.js` for the chunk LRU + index lookup.

These two unblock the L3 panel renderers (which use het rate to
color the mini-PCA dots in v3.94+ "het coloring" mode).

---

## Step 3 — `inversion_review.html` + `inversion_review/`

**Status**: not started

**Scope**: page11 + page_sv_evidence + page4 + page7 + page6 + the
existing band-trace UI + G-panel auto tab. Behavioural parity only;
new SPEC BLOCK 2 features (auto-promote / bulk actions / sample-
concordance proposals) layer on AFTER migration is complete.

---

## Step 4–6 — catalogue / comparative / parity verification

**Status**: not started

---

## Step 7 — NEW work begins

**Status**: not started

SPEC BLOCK 2 features (review page auto-promote etc.).
SPEC BLOCK 1 R-module (LANTA-side band track extraction).

---

## Notes for next chat

The legacy atlas tarball uploaded to this chat was truncated at
the upload pipeline at ~11 MB (large LG28.json eats the budget).
The truncated file ends mid-`_bundleTierBlock`. None of the
primitives extracted in Step 1 were affected — they all sit
between line 12250 and line 39220, well before the truncation point.

To proceed with Step 2 cleanly, ask the user to bundle the next
tarball **without `data/precomp/<chrom>/<chrom>.json` files**. Drops
the tarball to ~250 KB, uploads cleanly, and the R-pipeline outputs
can be fetched separately on the dev machine.
