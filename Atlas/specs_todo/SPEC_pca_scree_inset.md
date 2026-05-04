# SPEC — PCA scree-plot inset on the tracked-sample plane

**Status**: drafted 2026-05-03, scope = page 1 PCA scatter (`#pcaPanel` →
`#pcaCanvasWrap`).
**Schema impact**: `SCHEMA_V2.md` 2.15 → 2.16 (add `lam_top_k` per
window).
**Decision basis**: user request — *"At least a mini scree plot it's ok.
The DA no need bc its not applicable to PCA. ... Idk you can try its just
to see like yeah is my window stable or did our window is all unstable
and we have more variations in the next PC its a mini diagnosis"* — and
*"for the first 4-5 PCs if possible"*.

## 1. The diagnostic question

When the user is on a window in the per-sample-lines/PCA scrubber, the
question to answer in 2 seconds: **is this window's local PCA dominated
by PC1 (real biaxis structure → likely inversion signal) or is variance
smeared across many PCs (noisy / unstable / no clean structure)?**

A scree plot of the top-k eigenvalues (k = 5 or 6) gives the elbow shape
that answers exactly this:

- **Steep elbow at PC1** (λ₁ ≫ λ₂ ≈ λ₃ ≈ ...) → real inversion-style structure.
- **Two-step elbow** (λ₁ > λ₂, both > λ₃) → two-axis structure, possibly
  two overlapping inversions or PC2-relevant covariate.
- **Flat / shallow** (λ₁ ≈ λ₂ ≈ ... ≈ λ₅) → noisy window, no dominant
  axis, variance smeared. Don't trust the K-means assignments in this
  window.
- **Sharp drop after PC2 with PC2 substantial** → classic single-inversion-
  with-frequency-imbalance signature.

This is one of the most-information-per-pixel diagnostics possible for a
local-PCA scrubber. Honest call: this is more useful than I initially
suggested.

## 2. Why precomp must change

`Inversion_atlas.html` line 27746–27747 explicitly notes:

> *"the only honest percentage we can give from per-window data alone
> (full-dimensional variance is not in JSON)"*

Today's per-window object carries `w.lam1`, `w.lam2` only. A real scree
needs at least 5 eigenvalues. So `export_precomp_to_json_v3.R` (the
canonical emitter — see line 4974 of the atlas) must learn to write
`lam_top_k`. Per-window cost: 5 doubles ≈ 40 bytes. Per-chrom cost:
40 × ~3000 windows = ~120 KB. For LG28 specifically that's a 0.3% size
increase on a 40 MB precomp. Negligible.

The R code change itself depends on whether the upstream local-PCA
pipeline already computes the full eigendecomposition. Two cases:

**Case A: full eigendecomp already computed.** The local-PCA worker
returns the full eigenvalue spectrum but the emitter only writes the
top-2. Patch is one line: write `eigenvalues[1:5]` instead of
`eigenvalues[1:2]`.

**Case B: pipeline short-circuits to top-2.** Many local-PCA implementations
use power iteration / RSpectra to extract just the top-2 eigenvectors for
speed. Extending to top-5 means switching to a full `eigen()` call or
`RSpectra::eigs_sym(..., k = 5)`. Per-window slowdown ~2-5× on a
nsamples×nsamples covariance matrix (still tens of milliseconds).

**Status**: I don't have the R source uploaded so I can't tell which
case applies. Quentin has it on LANTA at
`/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/` (per the
memory notes). The R-side patch is included as a separate file in this
turn's deliverables — **a unified diff against `export_precomp_to_json_v3.R`**
that handles both cases gracefully (Case B is the safe path: always
recompute with `RSpectra::eigs_sym(k=5)` rather than relying on whatever
the upstream worker returns).

## 3. Wire format — schema 2.16 §<scree>

Per-window addition to the existing `windows[]` array:

```json
{
  "windows": [
    {
      "idx": 0, "start": 1234567, "end": 1334567, "center_mb": 1.28,
      "z": 1.42, "lam1": 0.82, "lam2": 0.31,
      "lam_top_k": [0.82, 0.31, 0.18, 0.12, 0.08],   // NEW (schema 2.16)
      "n_snps": 412, "L1": "L1_3", "L2": "L2_5",
      "...": "..."
    },
    ...
  ]
}
```

Constraints:

- `lam_top_k` is an **array of 5 to 8 doubles** (k flexible to allow 6
  or 7 for cohorts where additional structure matters; default 5).
- Sorted descending: `lam_top_k[0] >= lam_top_k[1] >= ...`.
- `lam_top_k[0]` MUST equal `w.lam1` and `lam_top_k[1]` MUST equal
  `w.lam2` (when both are present); the atlas asserts this on load and
  warns if they disagree.
- Backward-compat: omitting `lam_top_k` is fine; the atlas falls back to
  `[lam1, lam2]` and shows a "spectrum not in precomp yet" hint in the
  inset.

## 4. Atlas-side rendering

### 4.1 Inset placement

Top-right corner of `#pcaCanvasWrap`. Absolutely-positioned div with
`top: 8px; right: 8px;`. Width: 110px. Height: 70px. Translucent panel
background. Renders an inline SVG bar plot.

The screenshot from Quentin (DAPC plot) showed two insets — top-right
(PCA eigenvalues) and bottom-right (DA eigenvalues). Quentin agreed
DAPC is not applicable to local PCA; we ship only the PCA eigenvalues
inset for now. The bottom-right corner stays free for a possible future
"per-band silhouette" inset (separate spec, not this turn).

### 4.2 Toggle UI

A new checkbox in `#pcaTrackedAside`'s settings row:

```
[ ] trails  [ ] sign-align PC1  [ ] lasso  [ ] scree
```

State: `state.screePlotEnabled` (boolean), persisted to
`localStorage.inversion_atlas.screePlotEnabled`. Default OFF (matches
the user's request — *"its not active by default"*).

When toggled on, the inset div appears; toggled off, the inset div is
removed. No layout reflow because the inset is absolutely-positioned.

### 4.3 Bar plot rendering

5 (or 6 / 7) bars in descending height. Each bar is drawn as a fraction
of `lam_top_k[0]` (the dominant eigenvalue), so the first bar is always
full-height. Colors match the existing PC1/PC2 palette where possible:

- bar 0: PC1 colour (the existing PCA scatter's PC1 highlight)
- bar 1: PC2 colour
- bars 2–5: descending grey shades, decreasingly opaque

A small label above the inset reads "PC eigenvalues" (small caps,
`var(--mono)`, 9px). Below the bars, a one-line ratio diagnostic:
`λ₁/λ₂ = X.X` formatted to 1 decimal.

When `lam_top_k` is missing (legacy precomp): show 2 bars only and a
small italic hint "spectrum k≥3 needs precomp ≥2.16". Doesn't break,
just degrades.

When the cursor is in the empty state (no window selected): the inset
hides itself (the toggle is on but there's nothing to plot).

### 4.4 Per-window vs per-L2 scope

Per-window (matches the cursor's scrubbing). Each cursor move updates
the inset. This is cheap (5 bars, ~30 SVG elements). I considered
per-L2 (more visually stable) but per-window is the right diagnostic
target — when the user is parking on a noisy window, the inset
immediately tells them "yeah this one is unstable".

## 5. Atlas implementation checklist

- [ ] `state.screePlotEnabled` (boolean, default false), persist to
  `localStorage.inversion_atlas.screePlotEnabled`.
- [ ] Toggle checkbox in `#pcaTrackedAside` (one new label after `lasso`).
- [ ] `_renderScreeInset(wObj)` — pure SVG renderer reading
  `wObj.lam_top_k` (with the `[lam1, lam2]` fallback). ~60 lines.
- [ ] `#screeInset` div, absolutely-positioned in `#pcaCanvasWrap`'s
  upper-right.
- [ ] Wire to scrubber updates: each time `state.cur` changes,
  re-render the inset if `state.screePlotEnabled` is true. Already a
  hook somewhere — `requestRender` or similar.
- [ ] Hide the inset when no window is selected (`state.cur < 0`) or
  when the toggle is off.
- [ ] CSS: `.scree-inset { ... }`, `.scree-inset-bar`, `.scree-inset-label`,
  `.scree-inset-ratio`, `.scree-inset-fallback-hint`. ~25 lines.

Estimated atlas line count: **~120 lines**.

## 6. R-side patch (separate deliverable)

`PATCH_export_precomp_lam_top_k.R` — a small additive change to
`export_precomp_to_json_v3.R`. Two-line shape:

```r
# Add to the per-window emit loop:
lam_top_k = head(window_eigenvalues_full, 5),
```

Where `window_eigenvalues_full` is the full (or top-k) eigenvalue vector
from the local PCA. If the upstream worker doesn't already return the
full spectrum, the patch includes a small wrapper that recomputes with
`RSpectra::eigs_sym(window_cov_matrix, k = 5, which = "LM")`.

The patch handles the schema bump: `meta$schema_version <- "2.16"` and
adds a release note. Backward compatibility: old precomps still load
(the atlas falls back to the 2-bar inset).

Quentin runs this patch on LANTA, regenerates the precomp for one
chromosome to verify, and the inset auto-extends to 5 bars.

## 7. Tests

- Source-level: helper presence, CSS classes, toggle wired, schema
  version handling.
- Behavioural: `_renderScreeInset` returns valid SVG with 5 bars when
  `lam_top_k` is provided; returns 2 bars + hint when only `lam1/lam2`
  are present; returns null when both missing.
- Toggle: localStorage persists across reloads (test the state-restore
  hook).
- Edge cases: `lam_top_k` with NaN values, `lam_top_k.length < 2`,
  `lam_top_k[0] === 0` (all bars zero-height).

## 8. What this is NOT

- **Not DAPC.** The screenshot inset showed both PCA and DA
  eigenvalues; we deliberately drop DA because the atlas doesn't run
  discriminant analysis. Quentin agreed: *"The DA no need bc its not
  applicable to PCA."*
- **Not a global PCA panel.** Per-window only. There is no
  cohort-level scree (that would be a separate analysis — global PCA on
  the whole genome, which lives in MODULE_3 outputs, not in the
  scrubber).
- **Not a band-silhouette inset.** That's a different diagnostic ("how
  good is the K-means clustering?") and would go in the bottom-right
  corner. Out of scope for this turn; potential follow-up if Quentin
  wants it.

## 9. Estimated effort

- Atlas-side: 1 focused turn (~120 atlas lines + ~30 test lines).
- R-side patch: 5 lines of R + a comment block. Quentin runs on LANTA.
- Schema doc: half a page in `SCHEMA_V2.md`.
- Total: small, additive, fully back-compat.
