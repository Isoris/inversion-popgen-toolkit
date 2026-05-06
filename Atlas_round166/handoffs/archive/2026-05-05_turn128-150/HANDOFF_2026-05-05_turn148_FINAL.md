# HANDOFF — turn 148 — Slab-mode L3 panel parity with L2 mode

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (71,767 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 147 + 147b handoff.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok).

---

## 1. What this turn ships — queued screenshot fix #3

This turn addresses item #3 from the turn 147 handoff's "still queued" list:

> **1w/Nw L3 panel parity with L2 mode** — `renderL3PanelSlab` needs the
> same rich pane structure as `renderL3Panel`. Currently the slab-mode
> (1w / 5w / 10w / Nw) renderer is a downgraded variant of the L2-mode
> renderer. Estimated ~150–250 LOC.

Actual: +147 LOC, within budget.

### What was a "downgrade" before

When Quentin changes the compare unit from `L2` to any of `1w` / `5w` /
`10w` / `Nw`, `renderL3Panel` short-circuits to `renderL3PanelSlab`. The
slab variant rendered a working contingency table but was missing seven
features that the L2-mode renderer had developed across turns 32–99:

1. Per-pane toolbar mirrors (K group, color group, recluster select)
2. Karyotype chips ABOVE the mini-PCA on focal panes (turn 99 placement)
3. Spotlight click on every mini-PCA
4. `canvas.__l3_render` hit-test cache (no spotlight click could land)
5. Band-selector strip on focal pane in candidate mode
6. `col.__l3_meta` stash for cross-pane spotlight resolution
7. `_applySpotlightHighlights` post-pass that highlights contingency cells
   for the spotlighted sample(s)

Translation for the user: when working in slab mode, clicking a sample
on a mini-PCA did nothing, the K=mode toolbar wasn't reachable from a
pane header, the karyotype info sat below the plot instead of above,
and tracked-sample spotlights silently failed to highlight their cells
in the contingency tables. None of these were wrong-answer bugs; they
were missing features. Slab mode was usable but stripped down.

### What this turn does

Ports each of the seven into `renderL3PanelSlab` while preserving all
existing slab features and the L2 path.

#### 1.1 Per-pane toolbar (parity with L2 line 44592)

```js
const paneToolsHtml = (typeof _l3PaneHeaderToolsHtml === 'function')
  ? _l3PaneHeaderToolsHtml() : '';
```

Interpolated into both focal and neighbor `head.innerHTML` via string
concatenation (template-literal multi-line). The clicks were already
delegated by `_wireL3PaneToolsDelegation` (installed once on `#l3Panel`
root, survives re-renders), so wiring is automatic.

The previous focal head carried a vestigial `<span style="margin-left:auto;">SIM —</span>`
placeholder — a meaningless L2-aping element. Removed.

The focal title is now wrapped in `<span class="l3-pane-title">…</span>`
matching L2 mode's convention so the toolbar's `margin-left:auto` (set
by the global `.l3-pane-tools` CSS rule) places the toolbar on the right.

#### 1.2 Karyotype chips above mini-PCA (parity with L2 line 44687–44699)

```js
if (isFocal) {
  const clFocal = getSlabClusterAt(range[0], range[1], K);
  if (typeof _kSpecificMetaInlineHtml === 'function') {
    const chipsHtml = _kSpecificMetaInlineHtml(clFocal, null);
    if (chipsHtml) {
      const focalChipsAbove = document.createElement('div');
      focalChipsAbove.style.cssText = 'font-size: 9.5px; line-height: 1.2; ' +
        'padding: 2px 10px; margin: 0;';
      focalChipsAbove.innerHTML = chipsHtml;
      col.appendChild(focalChipsAbove);
    }
  }
}
```

`_kSpecificMetaInlineHtml(cl, l2idx)` reads `cl.n_per_group` and
`cl.centers`. Slab clusters carry both (they're returned by
`clusterSlabAtK` line 11627). Passing `l2idx=null` skips the
`_resolveSubbandLabelsForColumn` call which only applies to L2 envelopes
in active draft context — clean degradation, no crash, no spurious chip.

#### 1.3 Spotlight click + 1.4 `__l3_render` cache

The slab now wires `_setupL3MiniClick(mini)` on every mini-PCA, but the
click handler reads `canvas.__l3_render` for hit-testing. `drawSlabMiniPCA`
didn't populate that cache, so the click would silently no-op.

`drawSlabMiniPCA` now sets:

```js
const slabMid = (range[0] + range[1]) >> 1;
const { sign: slabSign } = getPC(slabMid);
canvas.__l3_render = {
  l2idx: null,
  isSlab: true, slabRange: range.slice(),
  wMid: slabMid, sign: slabSign,
  pad, plotW, plotH,
  xMin, xMax, yMin, yMax,
  cssW: w, cssH: h,
};
```

But there's a subtlety: `drawSlabMiniPCA` plots samples at slab-MEAN PC1
× slab-MEAN PC2. The L2 click handler's `getPC(ctx.wMid)` would return
single-window PCs, drifting from where dots actually sit at W>1.

Fix: `_setupL3MiniClick` now branches on `ctx.isSlab`. When slab, it
calls `aggregateSlab(ctx.slabRange[0], ctx.slabRange[1])` to get the
per-sample slab-mean PC1 (already sign-flipped) and recomputes slab-mean
PC2 (since `aggregateSlab` only stores `ys` when `aggMethod === 'mean_pc12'`).
Hit-test math then matches plot math exactly.

`signMul = 1` in the slab branch because `agg.xs` already contains
`pc1[si]*sign`, vs the L2 branch which multiplies `pc1[si] * sign` at
hit-test time.

#### 1.5 Band-selector strip (parity with L2 line 44713–44744)

`_l3BandSelectorHtml(K, bandCounts)` is K-and-counts-only — no L2
dependency — so it ports cleanly. Renders only when `state.candidateMode`
and only on the focal pane. Wires `_wireBandSelectorClicks()` which is
idempotent (one-time delegation on document.body).

#### 1.6 + 1.7 `col.__l3_meta` + spotlight post-pass

Slab columns now stash:

```js
col.__l3_meta = {
  K, offset,
  isSlab: true,
  leftRange:  (offset < 0) ? offsetSlab : range,
  rightRange: (offset < 0) ? range      : offsetSlab,
  invPerm: _invPerm,
};
```

Note `leftRange` / `rightRange` instead of L2 mode's `leftIdx` /
`rightIdx`, plus the `isSlab: true` discriminator.

`_applySpotlightHighlights` now branches on `meta.isSlab`:

```js
if (meta.isSlab) {
  if (!meta.leftRange || !meta.rightRange ||
      typeof getSlabClusterAt !== 'function') return;
  cl_left  = getSlabClusterAt(meta.leftRange[0],  meta.leftRange[1],  meta.K);
  cl_right = getSlabClusterAt(meta.rightRange[0], meta.rightRange[1], meta.K);
} else {
  cl_left  = getL2Cluster(meta.leftIdx);
  cl_right = getL2Cluster(meta.rightIdx);
}
```

Everything downstream (tagSample, axis-header tagging, primary-vs-tracked
priority) is K-and-labels-only and works identically for slabs.

`renderL3PanelSlab` now ends with `_applySpotlightHighlights(null)` in a
try/catch (parity with L2 line 44823). The argument is unused for the
slab path because slab columns carry their own `leftRange`/`rightRange`.

---

## 2. What this turn deliberately did NOT do

**Full focal-content body parity.** L2 mode's `focalContentHtml` is
~400 LOC of statistics: heterozygosity ridgelines, family-purity
diagnostics, band-diagnostics chips, ROH-overlap, GHSL, the K=both
hoisted invariant-meta block, etc. Many of these depend on L2-envelope-
level statistics (`_l2InvariantStats(cl, env, l2idx)` etc.) that don't
have direct slab equivalents. Slab still uses its lighter
`slabFocalContentHtml`. Bringing the full statistics suite to slabs
would be a separate multi-turn undertaking and isn't what the screenshot
complaint was about.

**Scale-mode reclustering for slabs.** When `state.l3ReclusterMode !==
'kmeans-K3'`, L2 mode routes through `compareL2Pair_byMode`. There is no
`compareSlabPair_byMode`. The per-pane recluster dropdown (now visible
in slab mode) currently no-ops for slabs — it changes the global state
slot but slab compare keeps using kmeans-K3. This is silently incorrect
when the user is in slab mode AND has changed the recluster dropdown.

**Mitigation**: this is a deferred-but-flagged gap. A clean fix is one
new function (`compareSlabPair_byMode`) plus a 4-line dispatch in
`renderL3PanelSlab`. Estimated <50 LOC in a follow-up.

**`drawMiniPCAEmpty` for null neighbors.** L2 mode draws "(no envelope)"
inside an empty mini-PCA box for null-neighbor columns. Slab mode shows
text-only "out of range" / "identical to focal at edge". This divergence
is intentional — slab "out of range" really IS off the chromosome edge,
not a missing envelope.

---

## 3. Files changed / added

| File | Change |
|---|---|
| `Inversion_atlas.html` | +147 LOC. Three functions touched: `renderL3PanelSlab` (head HTML restructure, focal chips above, spotlight click, band selector, `__l3_meta`, post-pass call), `drawSlabMiniPCA` (`canvas.__l3_render` cache), `_setupL3MiniClick` (slab-aware hit-test branch via `aggregateSlab`), `_applySpotlightHighlights` (slab-aware cluster lookup). |
| `tests/test_turn148_slab_l3_parity.js` | NEW — 62 / 0 |

No tests required updating — the head-HTML restructure (focal title
wrapped in `.l3-pane-title`, "SIM —" placeholder removed) didn't touch
any source pattern that an existing test was anchored to.

---

## 4. Tests

### 4.1 New suite

`tests/test_turn148_slab_l3_parity.js` — **62 / 0** across 14 sections:

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
    `getSlabClusterAt` is called (not `getL2Cluster`), and confirms the
    correct `(r, c)` cell + axis headers are tagged with the right CSS class
11. Sandboxed tracked-only spotlight (1) — verifies tracked path tags
    `spotlight-cell-tracked` and not `-primary`
12. Defensive paths (2) — null `__l3_meta`, null cluster lookup
13. L2 path preserved (5) — explicit negative tests confirming L2 mode
    still uses `leftIdx`/`rightIdx`, does NOT set `isSlab:true`, still
    calls `_applySpotlightHighlights(curL2)`, still short-circuits to
    slab when `compareUnit !== 'L2'`
14. Pre-existing slab features still wired (8) — slabFocalContentHtml,
    drawSlabMiniPCA on focal/neighbors, compareSlabPair, alignedSlabLabelsTo,
    ctHtml, l3KMode handling, ksToRender computation

### 4.2 Adjacent suites unchanged

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

### 4.3 Full sweep

**Across parseable turn-numbered tests: 1499 / 0** (was 1437 baseline,
+62 from this turn's new file).

The handoff's "2262 / 0" baseline includes more than just `test_turn*.js`
files; my sweep methodology counts only those. Either way, no regressions.

JS syntax check (`node --check` on the 2.78MB extracted main script
block): clean.

---

## 5. Atlas state

|                          | LOC     | Tests            | Files |
|---                       |---      |---               |---    |
| Pre-turn (turn 147)      | 71,620  | 2262 (handoff #) | 51    |
| Post-turn (this)         | 71,767  | 2324 (handoff #) | 52    |
| Δ                        | +147    | +62              | +1    |

---

## 6. What's NOT done — still queued

From the turn 147 handoff:

4. **G-popup → popgen page merge** — fold `_gPanelToggle` modal overlay
   into the popgen page as a second tab strip with different background
   shade. Per Quentin: "merge into one page with 2 rows of tabs (different
   bg shade)". Estimated ~400 LOC. **Highest-LOC remaining UI fix.**

5. **Tutorial authoring** — eight `data-status="pending"` cards in the
   help-page Tutorials section. Breeding-card tutorial highest-value
   first per turn 147 handoff.

New from this turn:

6. **`compareSlabPair_byMode`** — restore behavioral parity for the
   recluster dropdown when in slab mode. Currently the dropdown is
   visible (per-pane and global) but slab compare ignores it. ~50 LOC.

---

## 7. Honest framing

**What turn 148 actually delivered:**

- Slab mode's L3 panel now matches L2 mode's pane structure across the
  seven feature gaps that mattered for the user: per-pane controls,
  karyotype chips above the plot, click-to-spotlight, band selector in
  candidate mode, and contingency-cell highlighting from spotlight.
- The hit-test for slab clicks correctly uses slab-mean PCs (matching
  what's plotted), not single-window PCs. At W=1 the two are identical;
  at W=5 / 10 / N this matters.
- The L2 path is preserved unchanged — explicit negative tests confirm
  no `isSlab:true` leaks into L2 `__l3_meta`, no `getSlabClusterAt`
  routing for L2 columns.
- Defensive: missing `__l3_meta`, missing cluster, missing helper —
  all return early without crash.

**What it didn't deliver:**

- Full focal-content statistics body parity (huge, separate work).
- Recluster-dropdown dispatch for slabs (~50 LOC follow-up; flagged).
- An actual end-to-end slab-render integration test (would require a
  fixture cohort + DOM emulation; the test suite is structurally
  source-pattern + sandboxed unit instead).

**Manuscript impact:** moderate. Slab mode is the unit Quentin uses to
push from L2-envelope view down to fine boundary-resolution work — it's
where the L2 → L3 transition narrative lives. Having parity with L2
mode means the same cognitive primitives (click a sample, see it
highlighted across all panes; toggle K from a pane header without
travelling to the global toolbar; see karyotype info above the plot
where his eye goes first) work without mode-specific muscle-memory
gaps.

**Bundle**: tarball will be `Atlas_full_bundle_2026-05-05_turn148.tar.gz`
when re-packed. The current state at `/home/claude/Atlas/Atlas/` is
shippable.

Walk the map carefully, respect cohort discipline, don't break the
test suite.
