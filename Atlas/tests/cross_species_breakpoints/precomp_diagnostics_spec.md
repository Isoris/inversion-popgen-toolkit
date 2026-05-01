# Precomp diagnostic-mode specification (atlas v4 turn 83)

This document specifies what the `STEP_R..._precomp.R` script needs to add
so the Inversion Atlas diagnostic mode (residual coloring, suspicion scoring,
PC3/PC4 axes, eigenvector loading track) shows precomp-derived numbers rather
than atlas-side approximations.

Status:
- The atlas already reads `pc1`, `pc2`, `lam1`, `lam2`, `band` per window.
- The atlas already supports `pc3`, `pc4` in the per-window object — the X/Y
  scatter dropdowns auto-populate when these are present. **Just running
  precomp with `--npc 4` is enough for PC3/PC4 to appear.**
- The atlas already computes per-fish residuals atlas-side from `pc1`+`pc2`
  when precomp doesn't ship them. Day-3 work is to add fields that improve
  the accuracy of the atlas's diagnostic numbers.

---

## What to add to precomp

### Already-supported (just enable them)

```r
# In the precomp R script, when calling the local PCA decomposition:
local_pca <- prcomp(dosage_matrix, center=TRUE, scale.=FALSE, rank.=4)
# rank=4 (or ncomp=4) instead of rank=2 — this gives you PC1..PC4
```

Then in the per-window record:

```r
window_record$pc1 <- local_pca$x[, 1]    # length = n_samples
window_record$pc2 <- local_pca$x[, 2]    # length = n_samples
window_record$pc3 <- local_pca$x[, 3]    # NEW — atlas already supports this
window_record$pc4 <- local_pca$x[, 4]    # NEW — atlas already supports this
```

That's all that's needed for PC3/PC4 to surface in the X/Y scatter dropdowns.
The `availablePCs()` function in the atlas auto-detects them.

### NEW — required for full diagnostic mode

#### 1. Per-window eigenvalues (lam3, lam4)

The atlas already shows `λ₁ XX% · λ₂ YY%` in the PCA scatter axis labels for
PC1×PC2. Extend this rule to PC3+ by emitting `lam3` and `lam4`:

```r
window_record$lam1 <- local_pca$sdev[1]^2
window_record$lam2 <- local_pca$sdev[2]^2
window_record$lam3 <- local_pca$sdev[3]^2    # NEW
window_record$lam4 <- local_pca$sdev[4]^2    # NEW
```

Atlas already reads `lam1`/`lam2` for the existing variance-fraction display;
adding `lam3`/`lam4` will let me extend that display to PC3/PC4 axes when they
are selected.

#### 2. Per-window per-sample residual_z (band_residual_z)

The atlas computes this client-side from `pc1`+`pc2` as a fallback. The
precomp-derived version is more robust because it can use:
- a Gaussian-mixture fit (mclust) instead of K-means, giving the actual
  posterior probability per band rather than a hard assignment
- per-band σ derived from the mixture covariance, not the rough sd we use

```r
# After K-means (or mclust) assignment to k bands
library(mclust)
fit <- Mclust(local_pca$x[, 1:2], G=k)        # k = number of expected bands (3 for diploid)
band_label <- fit$classification              # length n_samples, values in 1..k
mixture_prob <- apply(fit$z, 1, max)          # per-fish max posterior
# Residual: Mahalanobis distance from sample's PC position to its assigned
# band centroid, using the band's covariance
residuals <- numeric(n_samples)
for (i in 1:n_samples) {
  b <- band_label[i]
  pos <- local_pca$x[i, 1:2]
  centroid <- fit$parameters$mean[, b]
  Sigma <- fit$parameters$variance$sigma[, , b]
  residuals[i] <- sqrt(t(pos - centroid) %*% solve(Sigma) %*% (pos - centroid))
}
window_record$band_residual_z <- residuals       # NEW — required field name
window_record$band            <- band_label - 1  # NEW — 0-indexed for atlas
window_record$mixture_prob_max <- mixture_prob   # NEW — optional but useful
```

Atlas behaviour:
- If `band_residual_z` is present, the atlas uses it directly and reports
  `source: 'precomp'` in the diagnostic cache (vs `'atlas'` for the
  client-computed fallback).
- If absent, atlas computes residuals from `pc1`+`pc2` itself.

#### 3. Per-window eigenvector loadings (top-N or full)

Used for the **eigenvector loading track** (a thin track below the |Z| profile
showing where each PC does its work in genomic coordinates).

For L2 windows (smaller, many of them): emit only the top-50 loci per PC
to keep file size down. For L1 windows: full loadings are tractable.

```r
# Top-50 loci per PC by |loading|
top_loci <- function(loadings_vec, n=50) {
  ord <- order(abs(loadings_vec), decreasing=TRUE)[1:min(n, length(loadings_vec))]
  data.frame(
    locus_idx = ord - 1,                           # 0-indexed
    loading   = loadings_vec[ord]
  )
}

window_record$pc_loadings <- list(
  pc1 = top_loci(local_pca$rotation[, 1]),
  pc2 = top_loci(local_pca$rotation[, 2]),
  pc3 = top_loci(local_pca$rotation[, 3]),         # NEW
  pc4 = top_loci(local_pca$rotation[, 4])          # NEW
)
```

Or, even more compact, just per-PC summary statistics:

```r
window_record$pc_loadings_summary <- list(
  pc1 = list(max_abs=max(abs(local_pca$rotation[, 1])),
             rms=sqrt(mean(local_pca$rotation[, 1]^2))),
  pc2 = list(...),
  pc3 = list(...),
  pc4 = list(...)
)
```

The atlas can build the loading track from either shape.

#### 4. Per-candidate per-sample diagnostics block (optional)

Aggregating the per-window residuals across all L2 sub-windows of a candidate
gives the suspicion table on the karyotype/tier page. The atlas can compute
this client-side too (it already does), but precomp can pre-aggregate it
to save browser CPU on big cohorts:

```r
candidate_record$per_sample_diagnostics <- data.frame(
  sample_id              = sample_names,
  band_residual_z_mean   = ...,
  band_residual_z_max    = ...,
  band_consistency       = ...,    # fraction of windows in the modal band
  mixture_prob_max       = ...,    # min over windows of max-posterior
  suspicious_carrier     = ...     # boolean: ≥2 of 3 metrics fail
)
```

When present, the atlas reads it; when absent, atlas computes from the
per-window residuals.

---

## File-size impact

For 226 samples, ~5000 loci per L2 window, ~600 L2 windows:

- `pc3`+`pc4` per window: 226 × 2 × 4 bytes × 600 = ~1.1 MB total. Negligible.
- `lam3`+`lam4` per window: 8 bytes × 600 = ~5 KB total. Free.
- `band_residual_z` per window: 226 × 4 bytes × 600 = ~540 KB. Negligible.
- `mixture_prob_max` per window: same as residual_z. Negligible.
- `pc_loadings` (top-50 per PC × 4 PCs): 200 floats × 600 windows × 4 bytes
  = ~500 KB total. Manageable.
- `per_sample_diagnostics` per candidate: 226 × ~5 fields × 4 bytes × ~10 cands
  = ~50 KB. Free.

**Total addition: ~2 MB** on top of the current precomp file. Well within
budget for a precomp JSON that already runs at ~50 MB for this study.

---

## Migration order (suggested)

1. **Day 3 morning** — re-run precomp with `--npc 4` (or whatever flag you
   use to ask for 4 PCs). PC3/PC4 dropdowns will populate; variance display
   will show λ₁/λ₂ for any selected PC pair as soon as `lam3`/`lam4` ship.

2. **Day 3 afternoon** — add `band_residual_z` + `band` + `mixture_prob_max`
   per window. The atlas's residual coloring switches from `source: 'atlas'`
   to `source: 'precomp'` automatically; numbers become more accurate.

3. **Day 4 (optional)** — add `pc_loadings` for the eigenvector track. Atlas
   builds the track from whichever shape you ship (full / top-N / summary).

4. **Day 5 (optional)** — add `per_sample_diagnostics` per candidate to
   pre-aggregate the suspicion summary; saves browser CPU on big cohorts.

5. **Cross-page cluster registration** (whenever pages 2 & 3 ship clusters)
   — emit the top-level `cross_page_clusters` block (section 6 below).
   Atlas modes `cluster_theta_pi` and `cluster_ghsl` light up automatically.

6. **Q-ancestry block** (when NGSadmix outputs are wrangled) — emit the
   top-level `ngsadmix_q` block (section 7 below). The Q-ancestry color
   mode lights up automatically; legend swaps to figure-style Q-proportion
   bars per K-means cluster. Recommendation: project the pruned samples
   onto the inference-set Q-clusters so all 226 fish get colored.

7. **Detailed-mode classifier inputs** (whenever per-sample heterozygosity
   per candidate is computed) — emit `per_sample_het` per candidate (see
   section 8 below). Atlas-side classifier (`classifyDetailedCandidate`,
   v4 turn 90) immediately switches from `source: 'atlas_proxy'` to
   `source: 'precomp_het'`, producing more accurate per-sample divergence
   classes and band interpretation (hom-like / het-like / mixed).

8. **SNP density per window** (whenever variant counts are loaded) — add
   `n_snps` per window (see section 9 below). The SNP-density strip on
   the per-sample lines panel (v4 turn 95) switches from λ₁ proxy to real
   counts. Useful for interpreting diamond patterns where bands compress
   in low-density windows.

9. **Deferred** — SNP-distance residual (`band_dosage_distance_z`,
   section 5 below). Schema reserved; not implemented unless PC-residual
   proves insufficient for catching rare-arrangement carriers. Decision
   rule: revisit if a known rare-carrier slips through PC-residual.

Atlas is feature-complete on day 3; items 4–9 are optional / deferred.

---

## Atlas API for verification

```js
// In the browser console after loading precomp:
state.data.windows[0].pc3              // should be array of n_samples floats
state.data.windows[0].lam3             // should be a number > 0
state.data.windows[0].band_residual_z  // should be array of n_samples floats
state.data.windows[0].band             // should be array of integers in 0..k-1

// Verify diagnostic cache is using precomp source:
_diagClearCache();
const r = _diagComputeWindowResiduals(0);
console.log(r.source);   // should be 'precomp' (was 'atlas' before this work)

// Trigger residual recoloring:
state.colorMode = 'residual';
drawPCA();   // PCA scatter recolors blue→amber→red by per-fish residual
```

---

## v4 turn 84 additions — SNP-distance mode + cross-page cluster registry

### 5. SNP-distance residual (orthogonal to PC-residual) — **DEFERRED, schema-only**

**Status (v4 turn 85):** spec'd but **not implemented** in precomp or atlas.
The candidate inputs (`dosage.tsv` + `sites.tsv`, ~200 MB) exist on the HPC
but the cost/benefit is unclear:
  - Output JSON size: only ~540 KB extra (one z-score per fish per window).
  - Build-time cost: precomp would have to load and process the 200 MB
    dosage matrix per chromosome (probably already loaded for the PCA step,
    so marginal cost is small if integrated in the same pass).
  - Diagnostic value: the SNP-distance residual catches a **strict subset**
    of cases where PCA misses something — specifically, fish that project
    onto a band cleanly in PC1×PC2 but have an orthogonal SNP signature.
    Most rare-arrangement carriers are already caught by the PC-residual
    diagnostic (turn 83), which is implemented and works today.

**Decision rule for revisiting:** if the PC-residual diagnostic flags fewer
suspicious carriers than expected (e.g., a known rare-carrier fish slips
through PC-residual but you have orthogonal evidence it's atypical),
implementing the SNP-distance mode becomes worthwhile. Until then, leave
the spec here as a reserved field name + algorithm so the implementation is
a drop-in if/when it's needed.

**Question this would answer:** if a fish projects onto the same PCA band as
others but its actual sequence-level dissimilarity is large, it carries
something the PCA missed. The PC-residual catches positional outliers in
PC space; the SNP-distance residual catches sequence-level outliers.

**Not circular** — uses the same input matrix but a different metric. PCA
finds the dominant variance axes, projecting fish into 2D; SNP-distance
measures actual pairwise dissimilarity. Two fish at the same PC1 value can
have very different SNPs along orthogonal directions.

**Reserved field name:** `band_dosage_distance_z` (per-window, length n_samples)

**Reserved atlas mode:** `colorMode === 'residual_dosage'` — alongside the
existing `'residual'` (PC-space). Atlas should keep the mode dispatch as a
no-op fallback to grey until the field is shipped. (Not currently wired in
the code; add a single line to `getSampleColor` when implementing.)

**Reference algorithm** (R-side, when implementing):

```r
# Per window, after K-means assignment to bands:
band_label <- kmeans_band                    # length n_samples, 0..k-1
# Per-band SNP-by-SNP centroid (mean dosage profile)
band_centroids <- list()
for (b in 0:(k-1)) {
  members <- which(band_label == b)
  band_centroids[[b + 1]] <- colMeans(dosage_matrix[members, , drop=FALSE])
}
# Per-fish: Manhattan distance from sample's dosage profile to its assigned
# band centroid, normalised by the within-band Manhattan-distance σ
band_dosage_distance <- numeric(n_samples)
for (i in 1:n_samples) {
  b <- band_label[i] + 1
  band_dosage_distance[i] <- sum(abs(dosage_matrix[i, ] - band_centroids[[b]]))
}
# Z-score relative to within-band mean and sd of distances
within_means <- numeric(k)
within_sds   <- numeric(k)
for (b in 0:(k-1)) {
  ds <- band_dosage_distance[band_label == b]
  within_means[b + 1] <- mean(ds)
  within_sds[b + 1]   <- sd(ds)
}
band_dosage_distance_z <- numeric(n_samples)
for (i in 1:n_samples) {
  b <- band_label[i] + 1
  band_dosage_distance_z[i] <- (band_dosage_distance[i] - within_means[b]) / within_sds[b]
}
window_record$band_dosage_distance_z <- band_dosage_distance_z   # NEW
```

**Output file-size impact (when implemented):** 226 × 4 bytes × 600 windows
= ~540 KB. Same as the existing `band_residual_z` field. **Not** 200 MB —
the 200 MB is the input dosage matrix on the HPC, which precomp consumes
and reduces to a single z-score per fish per window.

**Atlas usage (when implemented):** an additional color mode `'residual_dosage'`
becomes available alongside the existing `'residual'`. The diagnostic value
is in the comparison: a fish where PC-residual is low but SNP-residual is
high is the canonical rare-arrangement signal — PCA missed it but the actual
sequence dissimilarity caught it. The atlas would show both side by side so
the user can spot disagreement.

**Migration path when revisiting:** add `band_dosage_distance_z` to the
per-window record in precomp; add a single line to atlas-side `getSampleColor`:
```js
if (mode === 'residual_dosage') {
  // Same shape as _diagSampleColor but reads w.band_dosage_distance_z
  // instead of computing residuals from PC1+PC2.
  return _diagSampleColorDosage(si, state.cur);
}
```
plus the corresponding `_diagSampleColorDosage` helper (mirror of
`_diagSampleColor` from turn 83). One UI button on `colorModeBar`.
**Total atlas-side work to enable: ~30 lines of JS.** The infrastructure
from turn 83 (residual color ramp, cache, fallback to grey) all reuses.

### 6. Cross-page cluster registration

The atlas's cross-page coloring needs cluster labels from page 1 (dosage),
page 2 (θπ MDS), and page 3 (GHSL). Each page registers its cluster
assignments per candidate; the atlas can then color by any registered
source while looking at any other page.

**Schema** (one entry per candidate per source):

```json
{
  "cross_page_clusters": {
    "dosage": {
      "cand_63": { "labels": [0, 1, 2, 0, 1, ...], "k": 3 },
      "cand_47": { "labels": [...], "k": 3 }
    },
    "theta_pi": {
      "cand_63": { "labels": [...], "k": 4 },
      ...
    },
    "ghsl": {
      "cand_63": { "labels": [...], "k": 3 },
      ...
    }
  }
}
```

This block goes at the top level of the precomp JSON. Atlas reads it on
load and populates `state._crossPageClusters` automatically; cluster modes
become available with no further wiring.

**File-size impact:** 226 ints × 3 sources × ~10 candidates = ~7 KB. Free.

**Note:** the `theta_pi` and `ghsl` sources are populated by the day-3+ work
on pages 2 and 3 respectively. The atlas-side mode dropdowns are already
in place (turn 84) and will gracefully fall through to grey when those
sources aren't yet populated.

---

## v4 turn 86 additions — Q-ancestry block (NGSadmix K=2..20)

For the new "Q ancestry" color mode on page 1's PCA scatter. The atlas
already has the UI (button, K dropdown, hard/blend toggle, per_cluster /
cohort legend toggle) and the rendering logic (color resolution + legend
swap with proportion bars). It greys out gracefully until precomp ships
the Q-vector block.

### Schema — top-level `ngsadmix_q` block

```json
{
  "ngsadmix_q": {
    "K_values": [2, 4, 6, 8, 10, 12],
    "by_K": {
      "2":  { "q_vectors": [0.95, 0.05, 0.55, 0.45, ...],
              "group_names": ["Q1", "Q2"] },
      "4":  { "q_vectors": [...], "group_names": [...] },
      "8":  { "q_vectors": [...], "group_names": [...] }
    },
    "best_K": 4,
    "evaluator": "evalAdmix"
  }
}
```

- `q_vectors` is row-major: `n_samples * K` floats, sample-major order matching
  the atlas's `state.data.samples` order.
- `group_names` is optional; defaults to `["Q1", "Q2", ...]` if absent.
- `best_K` is optional; if present, the atlas auto-selects it as the active K
  on load. Otherwise the smallest K is selected.
- `evaluator` is for documentation; the atlas doesn't use it.

### Atlas-side ingestion (when implementing the precomp loader)

```js
// After loading the precomp JSON:
if (state.data && state.data.ngsadmix_q) {
  const block = state.data.ngsadmix_q;
  for (const K of block.K_values) {
    const entry = block.by_K[String(K)];
    if (!entry || !entry.q_vectors) continue;
    _qaRegisterK(K, entry.q_vectors, entry.group_names);
  }
  if (Number.isInteger(block.best_K)) {
    _qaSetK(block.best_K);
  }
}
```

That's all the wiring needed. The Q-ancestry color mode lights up
automatically as soon as `_qaRegisterK` populates the registry.

### File-size impact

For 226 samples, K values [2, 4, 6, 8, 10, 12] = 42 K-columns total:
226 × 42 × 4 bytes ≈ 38 KB total. Free.

For the full K=2..20 sweep (19 K values, 209 K-columns total): ~190 KB.
Still negligible.

### NAToRA pruning consideration

Quentin's pipeline runs NAToRA pruning to ~81 unrelated samples for the
NGSadmix step. The Q-vectors are computed on the unrelated set, so the
`q_vectors` array is naturally length `n_unrelated * K` = `81 * K`, not
`226 * K`.

For the atlas to color all 226 fish (including pruned-out family members),
two options:

1. **Project the pruned samples onto the unrelated-set Q-clusters.** Each
   pruned sample inherits the Q-vector of its closest kin in the unrelated
   set. This is what NGSadmix does internally if you ask it to — flag the
   pruned samples and re-run inference holding the cluster centroids fixed.

2. **Ship two Q matrices** — one for the 81 unrelated samples (used for the
   NGSadmix biological inference) and one for the full 226 (used for the
   atlas display, projected). Tag the unrelated set so the atlas can mark
   them visually if useful.

Recommendation: **option 1** — produce a 226-row Q-vector by projecting,
ship that, and add a per-sample boolean `is_unrelated` to indicate which
fish were in the inference set. The atlas can use the boolean to mark
inference-set fish with a tiny dot or border (same idea as the FOCAL marker
on the hatchery health page).

```json
{
  "ngsadmix_q": {
    "K_values": [...],
    "by_K": { ... 226-row matrices ... },
    "is_inference_set": [true, false, true, ...]   // length 226
  }
}
```

The atlas-side wiring picks this up naturally; the visual marker for
inference-set vs projected can be added later if needed.

---

## v4 turn 96 additions — detailed-mode classifier inputs + SNP density

### 8. Per-sample heterozygosity per candidate

The detailed-mode classifier (`classifyDetailedCandidate`, v4 turn 90)
operationalises the framework doc: per fish, compute within-sample
divergence in the candidate's window range; per band, summarize the
distribution of within-sample divergence and label the band as
`hom-like` / `het-like` / `mixed` / `ambiguous`; then assign H-system
labels by median PC1 ordering of hom-like bands.

Atlas-side fallback: when this field is absent, the classifier uses a
PC-residual proxy (saturating function of `_diagComputeCandidateSuspicion`
output from v4 turn 83). This is correlated with within-band divergence
but not equivalent to true heterozygosity.

#### Schema — per-candidate `per_sample_het` field

```json
{
  "candidates": {
    "C_gar_LG28_d17L2_0001_03": {
      "K": 3,
      "locked_labels": [...],   // existing
      "start_w": 142,           // existing
      "end_w": 287,             // existing
      "per_sample_het": [       // NEW — Float32Array, length n_samples
        0.12, 0.45, 0.10, 0.42, 0.13, ...
      ]
    }
  }
}
```

#### Definition

For candidate C with windows `[start_w, end_w]`:

```
For each sample s:
  per_sample_het[C][s] = mean over windows w ∈ [start_w, end_w] of:
    fraction of het sites in fish s in window w
  (using the SAME variant set the local PCA ran on)
```

If a fish has no called sites in the candidate's window range, emit
`null` or `NaN`. Atlas-side will treat as "no data" and skip the fish.

Recommended scale: heterozygosity per fish per window in the range
`[0, 1]` where 0 = no het sites among called sites, 0.5 = ~50% het
sites (full heterozygote signature for a biallelic system), 1.0 =
unlikely (would mean every called site is heterozygous).

#### Atlas-side detection

The atlas calls `classifyDetailedCandidate(cand, cand.per_sample_het)`.
When `cand.per_sample_het` is a typed array of length n_samples,
classification source = `'precomp_het'`. When absent, falls back to
`'atlas_proxy'`. The classifier output structure is identical in both
cases; only the underlying numbers differ.

#### File-size impact

For 226 samples × ~50 candidates × Float32 = ~45 KB per chromosome.
Negligible.

---

### 9. Per-window SNP count

The SNP-density strip on the per-sample lines panel (v4 turn 95)
shows where local PCA loses resolution. When `n_snps` is per-window,
the strip is precise; when absent, atlas computes a λ₁ proxy from
PC1 variance which is correlated but not identical (a window can have
many SNPs but no structural signal, or few SNPs with strong signal).

#### Schema — per-window `n_snps` field

```json
{
  "windows": [
    {
      "center_mb": 12.34,
      "pc1": [...],           // existing
      "pc2": [...],           // existing
      "lam1": 18.2,            // existing (we already use this as proxy)
      "lam2": 4.1,             // existing
      "n_snps": 142             // NEW — number of variant sites in this window
    }
  ]
}
```

Alternative names atlas accepts (in priority order):
1. `n_snps` (preferred)
2. `snp_count` (alternate)
3. `lambda1` / `variance_pc1` (existing, used as proxy when no SNP count)

#### Definition

Number of variant sites the local PCA ran on for that window — i.e.,
the row count of the dosage matrix slice for that window. NOT the
total number of records in the input VCF (which would include indels,
multiallelics, low-callability sites that may have been filtered).
Should match exactly the input count to `prcomp()` for that window.

#### Atlas-side detection

`_snpDensityForWindow(wi)` falls through priorities:
- `precomp:n_snps` (priority 1)
- `precomp:snp_count` (priority 1, alt name)
- `precomp:lambda1` (priority 2)
- `precomp:variance_pc1` (priority 2, alt name)
- `proxy:pc1_variance` computed atlas-side from `pc1[]` (priority 3)
- `none` (no source)

The strip's color ramp normalizes per-visible-range, so even rough
proxies still convey relative information.

#### File-size impact

For ~5000 windows per chromosome × int32 = ~20 KB per chromosome.
Negligible.

#### Why this matters

The diamond patterns Quentin observed on real candidates (synchronized
slants in homozygote-like bands while heterozygote-like rails stay
flat) are best explained by per-window SNP-density variation — the
B-vs-not-B contrast survives reduced SNP support, but the A-vs-C and
B-internal contrasts lose resolution. Without `n_snps` the atlas can
only show the symptom (compressed bands) via the λ₁ proxy; with
`n_snps` it can show the cause directly.

---
