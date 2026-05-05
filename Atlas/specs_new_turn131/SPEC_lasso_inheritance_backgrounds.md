# SPEC — Lasso-driven inheritance backgrounds + enlarged per-sample lines page

**Status**: drafted turn 130 final session. Major design ask from chat
`f74cf5d4`. Not yet implemented. Estimated ~1.5–2 turns total.

**Trigger** (Quentin, chat `f74cf5d4`):
> *"In per sample lines but also separately in tracked sample PCA when we
> lasso some bands, example we lasso samples like 40 samples then it
> looks what bands are they belonging to the most like 95% belonging to
> band 5 (because we lasso band 5 but like maybe have some sample that
> we lasso by mistake anyway). So it will track the samples. As usual.
> But it should also draw background alpha intervals on per sample lines
> for the whole inheritance systems that they belong to more far in the
> chromosome, their linkage."*
>
> *"Bc in the end we could reuse the per sample lines in a new page and
> have it large and move cursor left right and choose groups and it
> would annotate everything. But also for candidate page it will show
> like our candidate if we click on band 3 it will show its genome wide
> linked inheritance as I explained in an enlarged per sample lines and
> have the table of linkage of inheritance genome wide."*

A unified "linkage exploration" surface where the user lassos fish (in
either per-sample-lines or tracked-PCA), and the atlas:
1. Shows the lasso's purity score per candidate ("96% of these fish are
   in band 5 of candidate I7")
2. Draws alpha-shaded intervals on per-sample-lines showing which
   *other* candidates these fish predominantly belong to
3. Surfaces a table of "genome-wide linked inheritance" for the lasso

---

## 1. The interaction

### 1.1 Sources of lasso

- **Per-sample lines panel**: lasso on the line trajectories
- **Tracked-PCA scatter**: lasso on dots in the focal mini-PCA
- **Candidate band click**: clicking a band selects all samples in
  that band (equivalent to a clean lasso of just one band)

All three feed into the same `lassoSelection` state slot.

### 1.2 Lasso purity computation

For each candidate `C` and each band `b` of `C`:

```
purity[C][b] = |lasso_fish ∩ band_fish[C][b]| / |lasso_fish|
```

For each candidate `C`:

```
best_band[C] = argmax_b purity[C][b]
best_purity[C] = max_b purity[C][b]
```

A candidate is **strongly linked** to the lasso when `best_purity[C] > 0.7`
and the band `best_band[C]` has at least 5 fish.

### 1.3 Background alpha intervals

For each candidate `C` strongly linked to the lasso:
- Draw an alpha-shaded interval on the per-sample-lines panel spanning
  `C.start_bp` to `C.end_bp`
- Color = the K-means color of `best_band[C]`
- Alpha proportional to `best_purity[C]`

Effect: when the user lassos band 5 of candidate I3 (LG28 at 15-18 Mb),
the per-sample-lines panel shows alpha-shaded boxes at every other
candidate where these same 40 fish predominantly fall in one band.
Visual confirmation of inheritance linkage.

### 1.4 Linkage table

Sortable table:

```
candidate    span_mb       best_band  purity   n_lasso_in_band
I1 (LG28)    15.1–18.0     band 2     1.00     40 (the lasso itself)
I7 (LG28)    21.5–22.8     band 5     0.96     38
I12 (LG12)   4.2–7.1       band 3     0.92     37
I17 (LG14)   18.2–19.6     band 1     0.84     34
...
```

Click any row → highlight that candidate on the chromosome strip.

## 2. Enlarged per-sample-lines page

A dedicated page (probably new top-level page) where:
- Per-sample-lines panel takes the full screen (or most of it)
- Cursor left/right scrubs across the chromosome
- Lasso behavior is the same as the small panel
- The linkage table sits as a sidebar
- Genome-wide chromosome navigator at the top to switch chromosomes

This is essentially a **dedicated linkage-exploration mode**. Useful when
the user is investigating a specific inheritance hypothesis ("does this
broodline carry the same arrangements at all my candidates?").

## 3. State + data flow

```
state.lassoSelection: Set<sample_id>  (existing)
state.lassoLinkage: {                 (new, computed from lasso)
  per_candidate: {
    [cand_id]: {
      best_band: number,
      best_purity: number,
      n_in_band: number
    }
  },
  strong_links: cand_id[]   // best_purity > 0.7
}
```

Recompute `lassoLinkage` whenever `lassoSelection` changes. Cheap if
candidate count ≤ 200; full genome (~500 candidates after L2-sweep
auto-promote) needs O(n_cand × n_lasso) compute, still fast.

## 4. Implementation slices

### Slice 1 — lasso purity computation (~0.5 turn)
- `_computeLassoLinkage(lassoFish)` walks all candidates, returns the
  per-candidate purity map
- Cache invalidation on lasso change
- Tests against synthetic candidates

### Slice 2 — alpha intervals on per-sample-lines (~0.5 turn)
- Render alpha-shaded boxes per strong-link candidate
- Color = best_band's K-means color
- Alpha = purity scaled to [0.1, 0.5]
- Hover → tooltip with purity %

### Slice 3 — linkage table sidebar (~0.3 turn)
- Sortable table of strong-link candidates
- Click → highlight
- Filter by purity threshold

### Slice 4 — enlarged per-sample-lines page (~0.5 turn)
- Full-page layout
- Cursor scrub + lasso reuse
- Sidebar with linkage table + chromosome navigator

### Slice 5 — candidate-page integration (~0.3 turn)
- On the page-2 candidate-focus page, clicking a band auto-lassos
  that band's fish
- The enlarged per-sample-lines page opens as a modal showing the
  band's genome-wide linkage
- Linkage table is per-candidate

## 5. Open questions

1. **Purity threshold**: 0.7 is a reasonable default. User-tunable.
2. **n_in_band threshold**: at least 5 fish in `best_band[C]` to count
   as a strong link. Avoid spurious matches on tiny bands.
3. **Lasso-mistake tolerance**: if 5% of the lasso are "by mistake" (in
   the wrong band), the purity drops to 0.95. Threshold 0.7 is forgiving.
4. **Multi-chrom lasso**: the lasso operates within the active chromosome.
   The linkage table can extend genome-wide once `chromSummary` is
   populated (per `SPEC_multichrom_load_orchestrator.md`).
5. **Visual clutter**: alpha-shading lots of candidates makes the panel
   noisy. Cap at top-10 strong links by purity. Toggle to show all.

## 6. Tests

- Synthetic lasso of all band-2 fish at candidate I3 → purity at I3
  band 2 = 1.0, others lower.
- Lasso with one mistake fish → purity at I3 band 2 = (n-1)/n.
- Alpha intervals render at correct candidate spans.
- Linkage table sorted by purity descending.

## 7. Cross-references

- `SPEC_distant_band_concordance_fish_trajectory.md` — fish-trajectory
  lineage compute is a *non-interactive* version of this same idea
  (which fish co-cluster across all candidates). Lasso is the
  interactive version.
- `SPEC_g_panel_unified_groups.md` — lasso selections can be saved as
  manual groups via the G-panel.
- Chat `f74cf5d4` for the original design conversation.
- Existing lasso infrastructure (state.lassoSelection) shipped earlier.

## 8. What this is NOT

- Not a clustering algorithm. The lasso is user-driven; this spec just
  computes purity given the lasso.
- Not a replacement for the inheritance Jaccard matrix. The matrix is
  the global structure; lasso is local exploration.
- Not a classifier. It doesn't say "these fish are lineage A."
  It says "these fish co-cluster at candidates X, Y, Z."
