# SPEC — Sliding-window inheritance mode

**Status**: drafted turn 130. Not yet implemented. Estimated ~1 turn after
the L2-sweep mode (`SPEC_l2_sweep_inheritance.md`) ships.

**Trigger** (Quentin, turn 130):
> *"In the per sample lines and the line breadth taking 10 windows at a
> time or using L2 blocs and do it pairwise on the per sample lines lines."*

L2 envelopes have variable widths (some are tiny, some span tens of
windows). For finer resolution and to catch inheritance signals that
fall **between** L2 boundaries, we sweep a fixed-width window across
the chromosome, K-means each window, and feed those windowed clusterings
into the same inheritance machinery.

---

## 1. How this differs from L2-sweep

| | L2-sweep | Sliding-window |
|---|---|---|
| Items | Existing L2 envelopes | Fixed-width tiles (configurable: 5/10/20 windows) |
| Resolution | Variable (data-defined) | Uniform (user-defined) |
| Signal | "What does this L2 carry?" | "What does this stretch carry, regardless of L2 boundaries?" |
| K-means | Already cached per L2 | Recompute per tile (cached) |
| Use case | Confirms / extends user-promoted candidates | Catches inheritance structure between L2s, especially in dead zones |

Sliding-window is **strictly finer + noisier** than L2-sweep. It will
produce more groups but more spurious ones. The two complement each
other; neither replaces the other.

## 2. Algorithm

```
window_size = state.lineageSlideWindowSize    // default 10 windows
step        = state.lineageSlideStep          // default = window_size (non-overlapping)
n_total_w   = state.data.n_windows
tiles       = [[i, i + window_size) for i in range(0, n_total_w, step)]

for each tile (start_w, end_w):
    extract per-sample PC1 + PC2 averaged over tile windows
    K-means with K = state.k (default 3)
    record fixedKLabels[n_samples]
    record silhouette

items = tiles
    .filter(t => isUsableTile(t))
    .map(t => ({
      id: `tile:${t.start_w}-${t.end_w}`,
      labels: t.fixedKLabels,
      K: state.k,
      start_bp: state.data.windows[t.start_w].start_bp,
      end_bp: state.data.windows[t.end_w - 1].end_bp,
      meta: { source: 'sliding_window', start_w: t.start_w, end_w: t.end_w },
    }))

runInheritanceCompute({ force: true, items });
```

`isUsableTile` rejects tiles where:
- Fewer than 2 PCA bands populate (silhouette undefined or sub-threshold).
- Coverage of usable windows in the tile < 50% (most windows are NA /
  failed precomp).
- The tile's overall PVE1 / PVE2 indicates no real structure.

## 3. Compute cost

Per chromosome with 10,000 windows and tile size 10 → 1,000 tiles. Each
tile needs:
- Per-tile PC averaging: O(n_samples × tile_size) = 226 × 10 = 2,260 ops
- K-means with K=3: ~20 iterations × 226 samples = ~5k ops per tile
- Total: ~7k ops per tile × 1,000 tiles = 7M ops for K-means

Plus the inheritance Jaccard matrix on ~1,000 tiles × 3 bands = 3,000
band-items: 9M pairwise distances. This is the bottleneck.

Mitigation: chunked compute via `requestIdleCallback` so the UI never
blocks. Cache result per `(chrom, window_size, step)` in
`state.__lineageSlideCache`.

## 4. UI surfacing

### 4.1 L3 toolbar dropdown

Reuse the existing window-mode cycler vocabulary (`win5`, `win10`,
`winN`). Add a new picker `lineage:` with options:

- `off` (default)
- `L2` (delegate to L2-sweep when that ships)
- `5w` (sliding 5-window tiles)
- `10w` (sliding 10-window tiles)
- `Nw` (custom; pulls from `state.l3CompareUnitN` already wired in turn 128c)

State: `state.lineageScale` ∈ `{'off', 'L2', 'win5', 'win10', 'winN'}`.
localStorage: `pca_scrubber_v3.lineageScale`.

### 4.2 Lineage strip on per-sample-lines

Same strip as `SPEC_distant_band_concordance_fish_trajectory.md` §4.3.
Sliding-window mode produces denser strips (more pips per Mb).

## 5. Implementation slices

### Slice 1 — compute (~0.5 turn)
- [ ] `_buildSlidingTiles(window_size, step)` returns tile metadata
- [ ] `_kmeansForTile(tile, K)` recomputes if cache miss
- [ ] `runSlidingWindowInheritance(window_size)` orchestrator
- [ ] Cache: `state.__lineageSlideCache: Map<key, result>`

### Slice 2 — UI picker (~0.5 turn)
- [ ] Add `lineage:` dropdown to L3 toolbar
- [ ] Wire to `runSlidingWindowInheritance` / `runL2SweepInheritance` /
      `runLineageCompute` (the three engines, dispatched by mode)
- [ ] Re-render the lineage strip on change

## 6. Open design questions

1. **Tile overlap**: non-overlapping (step = window_size) is simplest
   and matches the L3 toolbar's existing `win5` / `win10` semantics.
   50% overlap (step = window_size / 2) would smooth boundaries but
   double the cost. Default: non-overlapping.
2. **Tile K**: matches `state.k` (3 or 6). Detailed mode runs K=6
   automatically.
3. **Auto-promote from sliding tiles?** Spec says NO — too noisy. The
   tile resolution is finer than the natural inversion length scale, so
   auto-promoting tile-as-candidate would produce candidate spam. User
   must confirm manually if they want one of these as a candidate.

## 7. What this is NOT

- **Not the L3 cycler's `win1/win5/win10` modes** (those are about the
  L3 contingency-table comparison unit, not about generating items
  for the inheritance compute).
- **Not real per-window K-means** at every individual window — that's
  the existing per-window precomp PC data. We're aggregating across
  tile windows here, not single-window analysis.

## 8. Tests (Slice 1)

- Synthetic 100-window chromosome with two planted inversions in
  windows 30-50 and 70-90 → sliding 10w produces tiles whose
  inheritance grouping recovers both inversions.
- Tile recompute: changing `window_size` invalidates cache; re-running
  with same size returns cached result.
- Performance: chromosome with 10k windows + tile 10 completes in
  under 2 seconds (chunked).
