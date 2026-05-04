# Turn 120 — PCA scree-plot inset on the tracked-sample plane

**Spec implemented**: `SPEC_pca_scree_inset.md` (drafted this turn).
**Sibling deliverable**: `PATCH_export_precomp_lam_top_k.md` — R-side
patch for `export_precomp_to_json_v3.R` to ship the new `lam_top_k`
field per window. Quentin runs this on LANTA at his convenience; the
atlas already supports both old and new precomp shapes.

## What shipped

A toggleable mini scree plot in the upper-right corner of the PCA
scatter (page 1, `#pcaCanvasWrap`). Renders top-k local-PCA eigenvalues
for the current window as a bar plot. Answers in 2 seconds: *"is this
window's variance dominated by PC1 (clean structure) or smeared across
many PCs (noisy / unstable)?"*

User flow:

1. Click the **scree** checkbox in the tracked-samples panel (next to
   trails / sign-align PC1 / lasso). Toggle persists across reloads.
2. Mini bar plot appears in the upper-right corner of the PCA scatter.
3. Scrub through windows — the bars update each cursor move (cheap, pure
   SVG, no canvas).
4. Each bar's height is `λᵢ / λ₁` (normalized to the dominant
   eigenvalue). Steep first bar = clean inversion structure. Flat bars
   = noisy window.
5. λ₁/λ₂ ratio shown below the bars to one decimal.
6. With current precomp (schema ≤2.15): 2 bars + hint *"k≥3 needs
   precomp ≥2.16"*. With patched precomp (schema 2.16): 5 bars (or up
   to 7 if Quentin chooses to ship more in the patch).

## What was deliberately NOT shipped

- **DA eigenvalues inset** (the second inset in Quentin's reference
  screenshot). Quentin agreed: the atlas doesn't run DAPC, so DA
  eigenvalues don't apply to local PCA. *"The DA no need bc its not
  applicable to PCA."*
- **Per-band silhouette inset** (the natural sibling diagnostic).
  Out of scope for this turn; potential follow-up.
- **Card flip on the candidate page** (mentioned by Quentin in turn
  119). Still deferred until the age-origin panel sits in production.

## Files modified / added

### `Inversion_atlas.html` (+~190 lines: 58,843 → 59,033)

- **HTML**: scree toggle checkbox in `#pcaTrackedAside`'s settings row;
  `#screeInset` div as a child of `#pcaCanvasWrap` (absolutely
  positioned, hidden by default).
- **CSS** (~50 lines): `.scree-inset`, `.scree-inset-label`,
  `.scree-inset-svg`, `.scree-inset-bar`, `.scree-inset-ratio`,
  `.scree-inset-fallback-hint`. Translucent panel background so the
  underlying scatter is still visible at the inset edges.
- **JS** (~140 lines):
  - `_SCREE_BAR_COLORS` palette (PC1 teal, PC2 grey, descending).
  - `state.screePlotEnabled` boolean with localStorage persist/restore.
  - `_renderScreeInsetHTML()` pure-SVG renderer. Reads `w.lam_top_k`
    if present; falls back to `[w.lam1, w.lam2]` + hint when only
    legacy fields are loaded. Defensive: filters NaN/zero/negative
    eigenvalues, sorts descending, caps at 7 bars.
  - `_refreshScreeInset()` — toggles visibility of the inset div and
    writes innerHTML.
  - Toggle handler wired to the new checkbox; persists state on
    change.
  - `_refreshScreeInset()` hooked into `drawPCA()`'s tail (try/catch
    fail-soft) so the inset auto-updates on every cursor move.

### Test added: `test_turn120_scree_inset.js` (33 assertions, all passing)

13 behavioural scenarios in a Node `vm` sandbox + 12 source-level
checks for helper presence, CSS classes, drawPCA hook, etc.

Behavioural scenarios:
- **Test 1**: scree disabled → empty render
- **Test 2**: no data loaded → empty
- **Test 3**: cur < 0 → empty
- **Test 4**: full `lam_top_k` (5 elements) → 5 `<rect>`s, no fallback
  hint, λ₁/λ₂ ratio shown
- **Test 5**: legacy `[lam1, lam2]` → 2 `<rect>`s, fallback hint with
  `precomp ≥2.16` mention
- **Test 6**: window object with no eigenvalues → empty
- **Test 7**: defensive descending sort — unsorted input gets sorted
  before rendering (verified via λ₁/λ₂ ratio computed from sorted
  values, not raw input)
- **Test 8**: cap at 7 bars when `lam_top_k` has 10 elements
- **Test 9**: NaN / zero / negative values filtered out before
  rendering
- **Test 10**: PC1 bar uses the documented `#7ad3db` palette color
- **Test 11**: label includes the current window index (`w0`, `w42`,
  etc.)
- **Test 12**: `cur` out of range (beyond windows.length-1) → empty
- **Test 13**: fewer than 2 valid eigenvalues → empty (need at least
  2 for a scree to make sense)

### Spec added: `specs_todo/SPEC_pca_scree_inset.md`

The full design document — diagnostic question, why precomp must
change, wire format, atlas-side rendering, R-side patch summary,
test strategy, what was deliberately not shipped.

### R-side patch added: `PATCH_export_precomp_lam_top_k.md`

A standalone document Quentin can read on LANTA and apply to
`export_precomp_to_json_v3.R`. Covers two cases: (a) full
eigendecomp already computed (one-line patch), (b) pipeline
short-circuits to top-2 (small refactor of the local-PCA worker call
to bump `k = 2` → `k = 5`). Includes a validation procedure
(regenerate one chrom, verify the atlas inset shows 5 bars).

## Verification

- `node --check` on full inline atlas JS: clean
- 33/33 turn-120 scree inset tests green
- 74/74 turn-119 age & origin tests still green (no regression)
- 82/82 turn-118 enriched bundle tests still green (no regression)
- 15/15 turn-117 focal-vs-bg tests still green (no regression)
- 14/14 turn-117b dot plot hover tests still green (no regression)
- **218/218 total no regression**

## Lines of code

| File | Before | After | Δ |
|---|---|---|---|
| `Inversion_atlas.html` | 58,843 | 59,033 | +190 |
| `test_turn120_scree_inset.js` | (new) | 257 | +257 |
| `specs_todo/SPEC_pca_scree_inset.md` | (new) | 217 | +217 |
| `PATCH_export_precomp_lam_top_k.md` | (new) | 145 | +145 |

## Roll-out

1. **Atlas-side is shipped now.** With current precomp (schema 2.15),
   Quentin can already toggle the scree on and see the 2-bar fallback
   inset with the *"k≥3 needs precomp ≥2.16"* hint. This is useful as a
   "is the visual idiom right?" sanity check before doing the precomp
   re-run.
2. **R-side patch is up to Quentin.** Apply
   `PATCH_export_precomp_lam_top_k.md`, regenerate one chrom (LG28 is
   the canonical proving ground), drop the new JSON in the atlas, see
   5 bars. If the visual idiom feels right, regenerate the rest of the
   chroms.
3. **No urgency on the precomp re-run.** The fallback inset is already
   informative — the λ₁/λ₂ ratio alone tells you "PC1 dominates" vs
   "no structure". The 5-bar version adds the elbow shape, which is
   nicer but not blocking.

## Known caveats / not-yet-addressed

- **Per-window scope, not per-L2.** Decided per-window because the
  cursor diagnostic is what Quentin asked for ("did our window is all
  unstable"). If per-L2 turns out to be more visually stable in
  practice (less flickering as the cursor slides inside one envelope),
  that's a one-line change in `_renderScreeInsetHTML`.
- **No interaction with the inset.** It's a passive diagnostic display
  — no click handlers, no tooltip on hover. Decided this fits the "1
  glance to read" use case better than an interactive panel. If
  Quentin wants click-to-drilldown later (e.g. open the windows page
  filtered to similar-shape windows), it's a small extension.
- **Bottom-right corner is free.** The original screenshot had a DA
  eigenvalues inset there. We dropped it; the corner is empty. If
  Quentin wants a per-band silhouette diagnostic later, it lives there.
