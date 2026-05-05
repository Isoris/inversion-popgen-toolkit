# HANDOFF — turn 148 + turn 149 — Slab-mode L3 panel parity (full close)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (71,851 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 147 + 147b handoff and the interim turn 148 handoff.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok).

---

## 1. What this session ships

This session closes the queued screenshot fix #3 from turn 147 — slab-mode
L3 panel parity with L2 mode — across two coherent turns. Turn 148 ports
the seven structural feature gaps; turn 149 closes the recluster-dropdown
gap that turn 148 introduced.

### Turn 148 — Slab-mode → L2-mode structural parity

When Quentin sets the compare unit to anything other than `L2` (i.e.
`1w` / `5w` / `10w` / `Nw`), `renderL3Panel` short-circuits to
`renderL3PanelSlab`. The slab variant rendered a working contingency
table but was missing seven features developed in the L2 renderer
across turns 32–99:

1. **Per-pane toolbar mirrors** (K group, color group, recluster select)
   via `_l3PaneHeaderToolsHtml()` — clicks already wired by
   `_wireL3PaneToolsDelegation`. Interpolated into both focal and
   neighbor `head.innerHTML`. The vestigial `<span>SIM —</span>`
   placeholder in the focal head was removed.
2. **Karyotype chips ABOVE the mini-PCA** on focal panes (turn 99
   placement). `_kSpecificMetaInlineHtml(cl, null)` works off
   `cl.n_per_group` + `cl.centers` — both available on slab clusters.
   `l2idx=null` skips the L2-only sub-band annotation chip cleanly.
3. **Spotlight click setup** on every mini-PCA via `_setupL3MiniClick`.
4. **`canvas.__l3_render` cache** populated by `drawSlabMiniPCA` (with
   `isSlab:true` and `slabRange`) so the click handler can hit-test.
5. **Slab-aware hit-test branch** in `_setupL3MiniClick`. The handler
   now branches on `ctx.isSlab` and uses `aggregateSlab(...)` to get
   slab-mean PC1/PC2 — matching what `drawSlabMiniPCA` actually plots.
   At W=1 the L2 single-window path and the slab-mean path give
   identical results; at W=5 / 10 / N this branch was the difference
   between landing on the right sample and being off by one.
6. **Band-selector strip** on focal pane in candidate mode via
   `_l3BandSelectorHtml(K, bandCounts)` — already K-and-counts-only
   so it ports cleanly.
7. **`col.__l3_meta` stash** on neighbor columns with `isSlab:true` +
   `leftRange`/`rightRange` (instead of L2 mode's `leftIdx`/`rightIdx`)
   + `invPerm`. **`_applySpotlightHighlights`** branches on
   `meta.isSlab` and routes cluster lookups via
   `getSlabClusterAt(range[0], range[1], K)`. Renderer ends with
   `_applySpotlightHighlights(null)` in a try/catch (parity with L2
   line 44823).

### Turn 149 — Slab recluster dropdown dispatch (loose-end close)

Turn 148 made the per-pane recluster dropdown visible in slab mode but
the slab compare silently used kmeans-K3 regardless of
`state.l3ReclusterMode`. That's worse than a missing button — silently
incorrect output is worse than visibly missing functionality. Turn 149
closes this:

- **New `compareSlabPair_byMode(leftRange, rightRange, mode)`** in the
  slab cluster API section. For `kmeans-K3` and `kmeans-K6` it routes
  to existing `compareSlabPair` with the right K. For U/V modes
  (`uv-rotated`, `uv-denoise`, `uv-dbscan`, `uv-dist-rank`,
  `uv-dist-fuzzy`, `distance-uv`) it falls back to kmeans-K3 with
  `fellBack: true` and `requestedMode: '<the U/V mode>'`. Honest
  degradation rather than silent ignore.
- **`renderL3PanelSlab` reMode dispatch** (parity with L2 line 44832):
  reads `state.l3ReclusterMode`; for non-default modes routes through
  `compareSlabPair_byMode`; for default keeps the simple
  `compareSlabPair` path (no behavioral change for unaffected users).
- **Visible fallback notice**: when `cmp.fellBack` is true, the
  contingency content gets prepended with a small accent-bordered
  notice div: "↩ slab fallback: **`<requestedMode>`** not available
  for slabs · using kmeans-K3", with an educational title attribute
  explaining why and what to do (use kmeans-K3/K6 in slab mode, or
  switch back to L2 mode for U/V).
- **`col.__l3_meta.K` uses `cmp.K`** (via `ctK`) instead of section K,
  so kmeans-K6 reMode against `l3KMode='k3'` sections gets the correct
  K=6 contingency dimensions for spotlight resolution. Same
  K-mismatch trade-off L2 mode accepts (mini-PCA palette stays at
  section K).

### Why U/V modes don't yet have slab implementations

The L2 U/V clusterers (`clusterL2_UVRotated`, `_UVDenoise`, `_UVDBSCAN`,
`_UVDistRank`, `_UVDistFuzzy`) all build per-sample (u, v) coordinates
via `_getOrComputeUVRotation(l2idx)`, which reads
`env._s0..env._e0` from the L2 envelope object. The math itself
(K-means3-for-centroids → principal-axis rotation → U/V projection,
then mode-specific clustering on top) is window-range-agnostic and
would port cleanly to slab ranges. The blocker is purely the
window-range plumbing.

A clean refactor would:

1. Extract `_computeUVRotationFromRange(s, e)` core.
2. Build `_getOrComputeUVRotationSlab(s, e)` with a slab cache.
3. Build `clusterSlab_UVRotated(s, e, K)` etc. mirroring the L2
   versions but reading the slab rotation cache.
4. Add five-mode dispatch in `compareSlabPair_byMode` similar to
   `getL2ClusterByMode`.

Estimated ~250 LOC across five clusterers + the rotation helper +
dispatch. Out of scope for this session — flagged for a future turn.

---

## 2. Files changed / added

| File | Change |
|---|---|
| `Inversion_atlas.html` | +231 LOC total (+147 turn 148 + +84 turn 149). Five functions touched: `renderL3PanelSlab` (head HTML restructure, focal chips above, spotlight click, band selector, `__l3_meta`, post-pass call, reMode dispatch with fallback notice), `drawSlabMiniPCA` (`canvas.__l3_render` cache), `_setupL3MiniClick` (slab-aware hit-test branch via `aggregateSlab`), `_applySpotlightHighlights` (slab-aware cluster lookup), `compareSlabPair_byMode` (NEW). |
| `tests/test_turn148_slab_l3_parity.js` | NEW — 62 / 0 |
| `tests/test_turn149_slab_recluster_dispatch.js` | NEW — 69 / 0 |

No pre-existing tests required updating mid-session for turn 148. One
turn 148 test had to be relaxed for turn 149 (the `invPerm = new
Array(K)` pattern became `new Array(ctK)` when K-mismatch was honored).
The relaxation is correct — the new behavior is what we want.

---

## 3. Tests

### 3.1 New suites

**`tests/test_turn148_slab_l3_parity.js` — 62 / 0** across 14 sections:
1. Per-pane toolbar in slab mode (4)
2. Karyotype chips above mini-PCA (3)
3. Spotlight click on slab mini-PCAs (1)
4. `drawSlabMiniPCA __l3_render` context (5)
5. Slab-aware hit-test branch in `_setupL3MiniClick` (5)
6. Band-selector strip in candidate mode (3)
7. `col.__l3_meta` on slab neighbor columns (4)
8. `_applySpotlightHighlights` call from slab renderer (2)
9. `_applySpotlightHighlights` slab branch (5)
10. Sandboxed slab-spotlight execution (7) — runs the actual ported
    function in a vm sandbox against a mock column, verifies
    `getSlabClusterAt` is called (not `getL2Cluster`), confirms the
    correct `(r, c)` cell + axis headers tagged with the right CSS class
11. Sandboxed tracked-only spotlight (1)
12. Defensive paths (2) — null `__l3_meta`, null cluster lookup
13. L2 path preserved (5) — explicit negative tests
14. Pre-existing slab features still wired (8)

**`tests/test_turn149_slab_recluster_dispatch.js` — 69 / 0** across 10 sections:
1. `compareSlabPair_byMode` signature (2)
2. Mode dispatch within byMode (4) — kmeans-K3, kmeans-K6, U/V fallback
3. Output shape annotations (5) — reclusterMode, requestedMode,
   fellBack, null guards
4. Slab renderer dispatch (5)
5. Fallback notice rendering (6) — accent border, title, prefix-not-replace
6. `__l3_meta` uses contingency K via `ctK` (3)
7. **Sandboxed `compareSlabPair_byMode` execution (32)** — every mode
   tested end-to-end: kmeans-K3, kmeans-K6, all six U/V modes,
   null mode, missing ranges, inner-null pass-through, unknown
   future mode (gmm-3 fallback)
8. Renderer integration sanity (2)
9. L2 path preserved (2)
10. Code organization (2)

### 3.2 Adjacent suites unchanged

All L3-touching turn tests run clean against the modified source:

- turn 128 spacebar-cut-in-candidate-mode: 20 / 0
- turn 128c cycler-l3-sync: 28 / 0
- turn 128c resize-dpr: 21 / 0
- turn 129 L3 het coloring: 58 / 0
- turn 132 page12 panel mirror: 78 / 0
- turn 134 L2 sweep inspector: 104 / 0
- turn 136 G-panel karyotype: 100 / 0
- turn 147 L3 toolbar disclosure: 66 / 0
- turn 147b windows-N engagement: 26 / 0

### 3.3 Full sweep

**Across parseable turn-numbered tests: 2393 / 0** (was 2262 baseline,
+131 from this session: +62 turn 148 + +69 turn 149).

Note on the sweep number: the previous "2262" baseline in turn 147
handoff is the same scrape methodology used here (counts both
`PASS:`/`FAIL:` and `N / N tests passed` summary formats). The "1568"
number reported in the interim turn 148 handoff was an undercount from
a less-robust parser, not a regression.

JS syntax check on the 2.78MB extracted main script block: clean.
HTML parser: 0 errors.

---

## 4. Atlas state

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 147)   | 71,620  | 2262 (sweep #)  | 51    |
| Post-turn-148            | 71,767  | 2324            | 52    |
| Post-turn-149 (this)     | 71,851  | 2393            | 53    |
| Δ session                | +231    | +131            | +2    |

---

## 5. What's still queued

From the turn 147 handoff, **two items remain**:

4. **G-popup → popgen page merge** — fold `_gPanelToggle` modal overlay
   (a 3-tab modal: karyotype / inheritance / manual) into the popgen
   page (`page12`) as a second tab strip with a different background
   shade. Per Quentin's brief: *"merge into one page with 2 rows of
   tabs (different bg shade)"*. Estimated ~400 LOC. **Highest-LOC
   remaining UI fix; would benefit from a screenshot or sketch from
   Quentin to ground the desired layout.**

5. **Tutorial authoring** — eight `data-status="pending"` cards in the
   help-page Tutorials section. Per turn 147 handoff: breeding-card
   tutorial highest-value first. **Caveat I discovered this session**:
   the tutorial launcher infrastructure isn't yet built either — the
   cards are static placeholder UI with no click-to-overlay machinery.
   So "authoring" a tutorial today actually means (a) building the
   overlay launcher infrastructure (~110 LOC) AND (b) writing the
   step content (~80 LOC per tutorial) AND (c) the tutorial really
   wants screenshots / animated GIFs for full effect. Quentin should
   decide whether the text-only version is worth shipping or whether
   to wait for visual-asset capability.

New from this session:

6. **Slab-aware U/V cluster modes** — `compareSlabPair_byMode` falls
   back to kmeans-K3 for all five U/V modes today. A future refactor
   could extract `_computeUVRotationFromRange(s, e)` from
   `_getOrComputeUVRotation` and ship slab variants of the five
   `clusterL2_UV*` functions. ~250 LOC.

---

## 6. Honest framing

**What this session actually delivered:**

- Slab mode's L3 panel now matches L2 mode's pane structure across the
  seven feature gaps that mattered for daily use: per-pane controls,
  karyotype chips above the plot, click-to-spotlight, band selector
  in candidate mode, and contingency-cell highlighting from spotlight.
- The hit-test for slab clicks correctly uses slab-mean PCs (matching
  what's plotted), not single-window PCs. At W=1 these are identical;
  at W=5/10/N this matters.
- The recluster dropdown — newly visible in slab mode after turn 148
  — now actually dispatches. Two K-means modes work; U/V modes
  visibly degrade with an in-place notice. No more silent
  "dropdown does nothing" trap.
- The L2 path is preserved unchanged — explicit negative tests
  confirm no `isSlab:true` leaks into L2 `__l3_meta`, no
  `getSlabClusterAt` routing for L2 columns, L2 mode still uses
  `compareL2Pair_byMode` for its dispatch.
- 131 new tests with sandboxed execution coverage of both the spotlight
  highlight path and the slab byMode dispatcher across all six U/V
  modes plus edge cases (null mode, missing ranges, unknown future
  mode, inner-null pass-through).
- Defensive everywhere: missing `__l3_meta`, missing cluster, missing
  helper, edge-case slab ranges — all return early without crash.

**What it didn't deliver:**

- Full focal-content statistics body parity (heterozygosity ridgelines,
  family-purity diagnostics, ROH-overlap, etc. — separate ~400+ LOC
  undertaking).
- Slab-aware U/V cluster implementations (~250 LOC follow-up; honestly
  flagged with both code-level reason and degradation path in place).
- An end-to-end slab-render integration test (would require fixture
  cohort + DOM emulation — out of scope for source-pattern + sandboxed
  unit tests, which is the existing test harness style).

**Manuscript impact:** moderate-to-substantial. Slab mode is the unit
Quentin uses to push from L2-envelope view down to fine
boundary-resolution work — it's where the L2 → L3 transition
narrative lives. Having parity with L2 mode means the same cognitive
primitives (click a sample, see it highlighted across all panes;
toggle K from a pane header without travelling to the global toolbar;
see karyotype info above the plot where the eye goes first) work
without mode-specific muscle-memory gaps. The recluster dropdown now
behaves honestly — kmeans modes work, U/V modes give a clear "not
yet" message.

**Bundle**: tarball
`Atlas_full_bundle_2026-05-05_turn149_slab_recluster_dispatch.tar.gz`.

Walk the map carefully, respect cohort discipline, don't break the
test suite.
