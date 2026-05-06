# Atlas-side reference: haplotype grouping, regime separation, inversion-system logic

**Purpose.** This document is a precise technical reference for what the atlas
computes and represents in the domain of haplotype grouping, regime
separation, and inversion-system inference. It is written so that a next
chat — given this document plus phase 4b's README plus the actual phase 4b
scripts — can do a script-by-script overlap analysis and decide which
logic to bring from where without rederiving context.

**This document is descriptive, not prescriptive.** It does not say what
should happen next. It says what the atlas does, with line numbers and
parameters, so a synthesis decision can be made empirically.

**Generated end of turn 113, 2026-05-01.** Reflects atlas state at line
count 49,558. All line numbers below refer to `Inversion_atlas.html` at
that snapshot.

---

## Table of contents

1. [The conceptual frame](#1-the-conceptual-frame)
2. [Atlas data model](#2-atlas-data-model)
3. [K-modes and the parallel candidate registry](#3-k-modes-and-the-parallel-candidate-registry)
4. [Per-window labelling and Hungarian alignment](#4-per-window-labelling-and-hungarian-alignment)
5. [The diamond detector](#5-the-diamond-detector)
6. [Sigma profile + detailed classifier](#6-sigma-profile--detailed-classifier)
7. [Structural-haplotype transition graph](#7-structural-haplotype-transition-graph)
8. [Band reach and regime breadth](#8-band-reach-and-regime-breadth)
9. [Per-sample lines panel and regime separation](#9-per-sample-lines-panel-and-regime-separation)
10. [L1/L2/L3 boundary scan](#10-l1l2l3-boundary-scan)
11. [Heterozygosity histogram and het-shape sub-module](#11-heterozygosity-histogram-and-het-shape-sub-module)
12. [SNP density confounder visualization](#12-snp-density-confounder-visualization)
13. [Q-ancestry display](#13-q-ancestry-display)
14. [Karyotype label vocabulary](#14-karyotype-label-vocabulary)
15. [What the atlas does NOT compute (and why)](#15-what-the-atlas-does-not-compute-and-why)
16. [Atlas state model: what's persisted, what's transient](#16-atlas-state-model-whats-persisted-whats-transient)
17. [Possible export shapes for the atlas → pipeline contract](#17-possible-export-shapes-for-the-atlas--pipeline-contract)

---

## 1. The conceptual frame

### 1.1 The simplification: 3 PCA bands → REF/HET/INV

The simplest model for a polymorphic inversion in a diploid cohort:
- One reference arrangement (REF) and one inverted arrangement (INV)
- Three diploid genotypes: REF/REF (HOMO_1), REF/INV (HET), INV/INV (HOMO_2)
- Local PCA inside the inversion gives 3 PCA bands corresponding to those 3 genotypes
- PC1 ordering: HOMO_1 < HET < HOMO_2 (or reversed; sign is arbitrary)
- The middle band has 2× the spread of the outer bands (heterozygotes mix
  the two haplotypes' signals)

This is what STEP_BP_01..04 currently assumes (`HOMO_1 / HET / HOMO_2`).

### 1.2 The reality: more haplotypes can coexist

Quentin's framing, recorded verbatim from a parallel chat:

> *Observed PCA bands: B1...B6. Inferred haplotype classes: H1, H2, H3.
> Diploid combinations: H1/H1, H1/H2, H1/H3, H2/H2, H2/H3, H3/H3. Why 3+
> bands ≠ simple REF/INV — multiple structural haplotypes can coexist
> with the same orientation, distinguished by sequence background,
> internal SVs, repeat history, etc. And the A/B/C alleles producing up
> to 6 diploid classes, with band order on PCA not necessarily following
> biological order.*

The biological mechanism: an inversion fixed in one lineage doesn't
destroy variation; it suppresses recombination between arrangements but
allows variation to accumulate independently within each arrangement.
After enough time, a single orientation can carry multiple haplotype
classes (H1, H2, H3, ...) that are detectable by local PCA but are NOT
the same as REF/INV.

### 1.3 Naming conventions in active use

Three conventions appear across Quentin's notes:

- **B1..B6** — observed PCA bands. Band index is by PC1 position, NOT
  biology. Band B3 might be H2/H2 in one candidate and H1/H3 in another.
- **H1, H2, H3** — inferred haplotype classes. Independent of orientation.
  An H1/H1 individual has two H1 haplotypes; an H1/H2 individual has one
  of each.
- **AA, AB, BB** — allele-pair shorthand. `AA` = homozygous for allele A,
  `AB` = heterozygous, `BB` = homozygous for allele B. Used when there
  are 2 alleles. With 3 haplotypes, the shorthand extends: AA, AB, AC,
  BB, BC, CC (6 diploid combinations).

### 1.4 What "regime" means

A **regime** is a stretch of the chromosome where the inversion-system
identity is stable. Inside one regime, the same set of haplotypes are
segregating at the same frequencies, and individual fish stay in the
same band across windows. A regime change is a position where the
underlying inversion system shifts: different haplotypes, different
breakpoints, or both.

This is distinct from a **boundary** (an inversion endpoint) and a
**recombinant edge** (a position where one fish's ancestry switches but
the population-level system is unchanged).

The atlas operationalises this distinction with three independent
signals:
- **L1/L2 boundaries** = changes in the population-level system structure
- **Diamonds** = changes in band stability within one system
- **Transition rate** = positions where many fish change band at once
- **Band reach / regime breadth** = how confined fish stay within bands

---

## 2. Atlas data model

The atlas reads precomp JSONs containing per-window data and extracts
per-window per-sample band labels via local-PCA + K-means. Every
algorithm below operates on this layered structure.

### 2.1 Top-level data

```
state.data = {
  chrom, n_windows,
  windows[]: {                             // length = n_windows
    start_bp, end_bp, n_snps,
    pca: {                                 // local PCA at this window
      pc1: Float32Array[n_samples],
      pc2: Float32Array[n_samples],
      lambda1, lambda2,                    // eigenvalues
    },
    z, ...                                 // |z|, lambda ratios, etc.
  },
  samples[]: { id, family, ... },          // length = n_samples = 226
  l1_envelopes[]: { start_w, end_w, peak_w, score },
  l2_envelopes[]: { start_w, end_w, peak_w, score, parent_l1 },
  ghsl_panel: { ... },                     // optional, blocked on Clair3
  theta_pi_panel: { ... },                 // optional, HPC export wiring pending
  ...
}
```

### 2.2 Per-L2 K-means clustering

For each L2 envelope, the atlas can run K-means on the per-window PC1/PC2
loadings to partition samples into K bands. This is computed lazily and
cached. The cache is keyed by `(L2_idx, K)`.

```
getL2Cluster(l2_idx) → {
  l2_idx,
  K: 3 | 6 | 8,                            // the K used
  usedK,                                   // actual K achieved (may downgrade)
  labels: Int8Array[n_samples],            // each sample's band, -1 if invalid
  fixedKLabels: Int8Array[n_samples],      // labels permuted to canonical order
  centers: Array[K] of {pc1, pc2},
  silhouette,                              // cluster quality metric
  ...
}
```

**Note**: `labels` is the raw K-means output (cluster ids 0..K-1 in
arbitrary order). `fixedKLabels` is the same labels permuted so band 0
is always the leftmost-on-PC1 band. The Hungarian alignment in §4 uses
`fixedKLabels` for cross-L2 comparison.

### 2.3 Per-window per-sample band assignment

For a sample at a window, the band is:
1. Find which L2 envelope contains this window (`l2_at_window`).
2. Look up `getL2Cluster(l2_idx).fixedKLabels[sample_idx]`.

So per-window labels exist BUT the K-means is per-L2, not per-window.
Within one L2, every window has the same per-sample labels by construction.

This is why "per-window class vector" output from phase 4b.1 maps cleanly
onto the atlas's L2 partitioning: one L2 = one class vector, repeated
across the windows of that L2.

---

## 3. K-modes and the parallel candidate registry

### 3.1 Three K-modes

The atlas supports K=3, K=6, and K=3+6 modes:

- **K=3**: classic REF/HET/INV partition. Sample falls into one of 3 bands
  per L2.
- **K=6**: finer partition. Used when 3 bands are clearly insufficient
  (e.g. multiple haplotypes per orientation).
- **K=3+6**: hybrid mode. K=3 verdict and K=6 verdict computed in
  parallel; cross-cutting cases (where K=3 says one thing and K=6 disagrees
  on band identity) trigger a `CROSS_CUTTING` verdict in pca_scrubber_v3.

Setting: `state.k = 3 | 6`. Persists in localStorage.

### 3.2 Parallel candidate registry (turns 87-89)

Two parallel candidate registries coexist:
- **default** (`state.candidates`): K=3-derived candidates
- **detailed** (`state.candidates_detailed`): K=6-derived candidates

Both are populated from the same precomp JSON but differ in the K-mode
used for clustering. The active mode is `state.activeMode = 'default' |
'detailed'`. Cross-realm-safe deep-clone: `_pcrDeepCloneCandidate(c)`
(JSDOM-safe; handles typed arrays via `constructor.name`).

Switching modes does NOT recompute clusters — it switches which registry
is active. Both are kept in memory simultaneously.

### 3.3 Merge isolation

`assertSameMode(a, b)` enforces that contingency tables and merge
operations stay within one K-mode. Mixing K=3 and K=6 candidates in one
contingency would corrupt the mode tracking. Helpers:
`buildContingencyForCandidates`, `mergeIsolationAudit`.

---

## 4. Per-window labelling and Hungarian alignment

### 4.1 The cross-L2 alignment problem

K-means assigns cluster ids 0..K-1 to bands in arbitrary order. Two
adjacent L2s can both have K=3 clusters but with cluster 0 = "leftmost
band" in one and cluster 0 = "rightmost band" in the other. To compare
band membership across L2s, the labels must be aligned.

### 4.2 The fix

`alignLabels(labelsA, labelsB, K)` computes a Hungarian-optimal permutation
`perm[K]` such that `labelsB` permuted by `perm` matches `labelsA` as
closely as possible (max overlap). Returns:

```
{
  perm: Int8Array[K],     // perm[k_A] = k_B_aligned
  matchCount: int,         // total matches under best permutation
  agreement: float,        // matchCount / n_samples
}
```

The convention: `perm[k_A] = k_B_aligned` means "band k in A maps to band
perm[k] in B." To project B INTO A's space, apply the inverse permutation.

### 4.3 Where alignment is used

- **Diamond detector** (§5): bands need consistent identity across
  windows of one L2. Within one L2 they already are; the diamond detector
  doesn't use Hungarian.
- **Transition graph** (§7): per-boundary stats use `compareL2Pair(left, right)`
  which calls Hungarian to align right's labels to left's space before
  counting transitions.
- **Band reach** (§8): chains alignments across multiple L2s so each
  fish's band-reach counts the same biological band consistently.

### 4.4 Limitations

Hungarian assumes the same K on both sides. K=3 ↔ K=6 transitions can't
be aligned this way; the K=3+6 mode handles those via the
`CROSS_CUTTING` verdict instead.

---

## 5. The diamond detector

**Defined**: Inversion_atlas.html lines 27225–27429, turn 92.
**Concept**: A "diamond" is a region where one band visibly splits into
sub-bands while at least one other band stays stable. It indicates
multi-haplotype structure within a single orientation.

### 5.1 Parameters

```
_DD_SPLIT_RATIO_THRESHOLD  = 1.5   // band sd > 1.5× baseline = "splitting"
_DD_MIN_DIAMOND_WINDOWS    = 3     // splitting must persist ≥3 windows
_DD_STABLE_VARIANCE_FRAC   = 0.05  // mean drift < 5% of total spread = stable
```

### 5.2 Algorithm

For each band b at each L2:

1. **Compute per-window sd**: for each window in the L2, compute the PC1
   standard deviation of samples assigned to band b.
2. **Compute baseline**: 25th percentile of the per-window sd values
   across the L2 (assumes the splitting is local; baseline reflects the
   stable region).
3. **Detect splitting ranges**: contiguous windows where sd > 1.5 ×
   baseline, persisting for ≥3 windows.
4. **For each splitting range**, check stability of OTHER bands:
   - A band is "stable" in a range if its mean PC1 drifts by < 5% of the
     total candidate spread.
   - Count stable vs slanting other-bands.
5. **Classify**:
   - `strict` if ≥1 other band is stable
   - `strict2` if ≥2 other bands are stable

### 5.3 Output shape

```
detectDiamonds(candidate) → [
  {
    splitting_band: int,       // band index that splits
    diamond_start_w, diamond_end_w,
    stable_bands: int[],       // bands that stay flat
    slanting_bands: int[],     // bands that drift in same range
    strict: bool,              // ≥1 stable
    strict2: bool,             // ≥2 stable
    baseline_spread, peak_spread, peak_spread_ratio,
  },
  ...
]
```

`summarizeDiamonds(candidate)` returns boolean+counts for catalogue
column rendering: `n_loose, n_strict, n_strict2, has_loose, has_strict,
has_strict2`.

### 5.4 Biological interpretation

A strict diamond strongly suggests **two parallel inversion systems** at
the same locus rather than recombination artifacts:
- If one band stays parallel-flat while another splits, the splitting
  band's samples are differentiated by something the stable-band samples
  are NOT differentiated by — i.e., a second segregating system.
- Recombination would tend to drift everyone, not split one band while
  others stay flat.

The detailed classifier in §6 reaches the same conclusion via a different
statistic (sigma profile bimodality).

---

## 6. Sigma profile + detailed classifier

**Defined**: Inversion_atlas.html lines 8055–8200, turn 90.
**Concept**: Per-sample PC1 standard deviation across the candidate's
windows, summarized by quantiles + bimodality coefficient. Different
profiles distinguish two-inversion-systems from crossover-artifacts from
noisy-region.

### 6.1 Per-sample sigma

For each sample s and each window w within the candidate, the sample's
band centroid is `centers[labels[s]].pc1`. The per-sample sigma is:

```
σ_s = sd_w (pc1[s, w] - centers[labels[s, w]].pc1)
```

i.e., how far from its band's centroid the sample drifts across windows.
Implemented as `sampleSpreadL2(l2idx)` and `sampleSpreadRange(start_w,
end_w)`.

### 6.2 Bimodality coefficient

Sarle's BC: `(skew^2 + 1) / kurt`. Threshold > 5/9 ≈ 0.555 = bimodal.

### 6.3 Classifier verdicts

```
sigmaProfileL2(l2idx, usedK) → {
  q50, q90, q95,            // sigma quantiles across samples
  n_high, ratio_high,       // count + fraction of σ > 2× q50
  bimodality_coef, is_bimodal,
  verdict: 'NA' | 'TWO_INVERSIONS' | 'CROSSOVER_ARTIFACTS' | 'NOISY_REGION' | 'UNDETERMINED',
  reason: string,
  top_high: [{si, sigma}, ...],   // sample IDs with highest σ
}
```

Verdict logic:
- **NA**: K ≤ 3 (no extra bands to explain)
- **TWO_INVERSIONS**: `q50 < 0.05 AND ratio_high < 0.05` AND K ≥ 4 →
  everyone is stable, K bands all parallel = stacked inversion systems
- **CROSSOVER_ARTIFACTS**: `is_bimodal AND 2% < ratio_high < 25%` →
  small high-σ tail = a few drifting samples = double-crossover or gene
  conversion in heterozygotes
- **NOISY_REGION**: `q50 > 0.10` → everyone has high σ = low power
- **UNDETERMINED**: anything else

### 6.4 What this tells you that diamond detection doesn't

Diamond detection looks at band shape (does one band split). Sigma
profile looks at sample shape (do samples drift from their centroids).

- TWO_INVERSIONS: bimodality low everywhere, K large → diamonds may NOT
  be present if both systems are individually stable
- CROSSOVER_ARTIFACTS: bimodality high with small tail → diamonds may
  not show because only a few samples are drifting

So §5 and §6 are complementary. Use both for triangulation.

---

## 7. Structural-haplotype transition graph

**Defined**: Inversion_atlas.html lines 27431–27572, turn 101.
**Concept**: Walk along the chromosome, at each L2 boundary count how
many fish change band. High transition rate = many fish switching at the
same coordinate = structural-haplotype regime boundary.

### 7.1 Parameters

```
_SHTG_HOTSPOT_THRESHOLD = 0.30   // transition_rate ≥ 0.30 = "hotspot"
_SHTG_MIN_SAMPLES = 5            // need ≥5 fish in both L2s to count
```

### 7.2 Algorithm

For each adjacent pair of L2 envelopes (i, i+1):

1. Get per-sample labels in each L2: `cl.fixedKLabels` and `cr.fixedKLabels`.
2. Hungarian-align right to left via `compareL2Pair(left, right)`.
3. For each sample s with valid labels in both L2s:
   - `flAligned = perm[fl]` (project left band into right's space)
   - if `flAligned ≠ fr`, the sample changed band → increment `n_changed`
4. `transition_rate = n_changed / n_valid`
5. Build edge map: `{from_band: {to_band: count}}` showing where samples went.

### 7.3 Output

```
computeStructuralHaplotypeTransitionGraph() → {
  boundaries: [
    {
      l2_left, l2_right,
      position_mb,                  // midpoint between L2 envelopes
      n_samples, n_changed,
      transition_rate,              // 0..1
      edges: [{from_band, to_band, count}, ...]   // sorted by count desc
    },
    ...
  ],
  hotspots: [...same shape, filtered to rate ≥ 0.30, sorted desc]
}
```

### 7.4 Biological interpretation

- **Low transition rate** (< 0.10): most fish keep same band across the
  boundary → no regime change at this position
- **Medium transition rate** (0.10–0.30): some fish change → either real
  but partial transition (e.g., recombinant fish at edges of a real
  boundary) or noise
- **High transition rate** (≥ 0.30, "hotspot"): regime boundary likely
  here. The edges reveal HOW fish move (e.g., "B1→B0: 80, B2→B1: 30"
  shows a coherent rotation of band identities)

This is distinct from a per-fish recombinant event — that's a single
fish changing band, which would barely move the population transition
rate. Hotspots are population-level signals.

---

## 8. Band reach and regime breadth

**Defined**: Inversion_atlas.html lines 27575–27780, turns 103-104.
**Concept**: Per-fish, count how many distinct bands the fish visits
across a chain of L2s. Aggregate to per-region "regime breadth": how
disciplined-vs-spread are the fish in this region?

### 8.1 Parameters

```
_BREACH_NARROW_THRESHOLD     = 2     // band_reach ≤ 2 = "narrow"
_BREACH_NARROW_FRAC_HIGH     = 0.70  // ≥70% narrow → regime is narrow
_BREACH_NARROW_FRAC_LOW      = 0.30  // <30% narrow → regime is wide
_BREACH_MIN_FISH_PER_BAND    = 5     // need ≥5 fish to call band populated
```

### 8.2 Algorithm

For an L2 chain `[i, i+1, ..., i+k]`:

1. **Get labels per L2**: K-means cluster labels for each L2.
2. **Chain Hungarian alignment**: leftmost L2 is identity. For each L2 j > 0,
   align L2[j] to L2[j-1]'s projected space. Apply the inverse permutation
   to project L2[j] INTO L2[0]'s space.
3. **Per-sample reach**: for each fish, count distinct projected bands
   visited across the chain (using a 32-bit bitmask per fish; reach =
   popcount of the bitmask).
4. **Per-region summary**:
   - `narrow_fraction = count(reach ≤ 2) / n_valid`
   - `bands_populated = count(bands with ≥5 visitors)`
   - `regime_breadth = 'narrow'` if `narrow_fraction ≥ 0.70` AND ≥2 bands
     populated; `'wide'` if `< 0.30`; `'medium'` otherwise; `'no_signal'`
     if no valid fish or <2 bands populated.

### 8.3 Output

```
computeBandReachAcrossL2s(l2_indices) → {
  l2_indices, n_samples, n_valid, K,
  per_sample_band_reach: Int8Array[n_samples],   // each fish's reach
  per_band_visit_count: Int32Array[K],
  narrow_fraction,
  bands_populated,
  regime_breadth: 'narrow' | 'medium' | 'wide' | 'no_signal',
}
```

`computeWindowedBandReachPerL2()` is the convenience wrapper: for each
L2 i, computes band reach over a 3-L2 sliding window `[i-1, i, i+1]`,
adds `center_l2`, `window_l2_indices`, `center_mb` to each entry.

### 8.4 Biological interpretation

- **Narrow regime**: most fish stay in 1-2 bands across multiple L2s →
  the inversion system is stable, fish carry consistent ancestry.
  Recombinants are rare.
- **Wide regime**: most fish visit 3+ bands across the chain → either
  the system is unstable (multiple regimes squeezed into the chain) or
  many recombinant samples.
- **Medium**: in between. Manuscript text might say "moderately
  recombinant background" or "transition zone".

The regime-breadth strip on the atlas (turn 104, `_drawRegimeBreadthStrip`)
displays this as colored chunks above the per-sample lines panel: green=narrow,
amber=medium, red=wide, grey=no_signal.

---

## 9. Per-sample lines panel and regime separation

The per-sample lines panel renders 226 fish as 226 lines along the
chromosome, one PC1 value per window per fish. This is the atlas's
primary visualization for regime separation.

### 9.1 Coloring modes (lcMode)

```
state.lcMode = 'default' | 'kmeans' | 'family' | 'dosage' | 'ghsl' | 'het' | 'theta_pi' | 'froh' | 'confounder_alert'
```

- **default**: all lines same color (low-alpha, just shape)
- **kmeans**: color by K-means cluster (band) at the focal L2. Shows
  which fish are in which band as the focal L2 changes.
- **family**: color by family (NAToRA-derived clusters). Shows whether
  related fish track together (suggesting common haplotype inheritance).
- **dosage / ghsl / het / theta_pi / froh**: color by an external
  per-sample-per-window scalar layer. Most are stubs awaiting their
  layer JSONs.
- **confounder_alert**: highlight samples flagged as potential
  confounders (low coverage, high missingness, etc.).

### 9.2 Rendering

`_resolveSampleColorByMode(sampleIdx, windowIdx, mode)` returns RGBA. The
lines code calls this per-sample per-frame; offscreen-canvas caches by
mode (cache key includes `|lc=`).

### 9.3 What this tells you that band-reach doesn't

Band-reach is an aggregate statistic. The lines panel shows the actual
trajectories. Two regimes with the same `narrow_fraction = 0.70` can
look very different:
- **Regime A**: 70% of fish stay in band 0; 30% bounce randomly. Lines
  look like a flat carpet plus 30% spaghetti.
- **Regime B**: 70% of fish stay in band 0 in the left half AND band 2
  in the right half (a crossing pattern). Same narrow_fraction but the
  trajectories tell a completely different story.

Visual inspection of the lines panel separates A from B. Band-reach
flattens the difference.

### 9.4 Trail and sign-align controls

```
state.trailOn (bool), state.trailN (int)    // how many windows of trail to show
state.flipPC1 (bool)                          // sign-align PC1 across windows
state.pc1Sign: Float32Array[n_windows]        // +1 or -1 per window
```

`flipPC1=true` flips the sign of PC1 in windows where the cohort
average drifts so all windows show "low" on the bottom. Without this,
PCA's arbitrary sign convention makes adjacent windows look like
opposite regimes when they're actually the same.

---

## 10. L1/L2/L3 boundary scan

The atlas's boundary detection lives in `STEP_D17_multipass_L1_only_v7.R`
and `STEP_D17_multipass_L2_v8.R` (R-side, runs upstream of the atlas).
The atlas READS the resulting envelopes from precomp JSON.

### 10.1 L1 envelopes

L1 = chromosome-scale envelopes. Boundaries detected by perpendicular-ray
boundary validation on the nn80-smoothed sim_mat. Each L1 envelope is a
contiguous window range with a peak (highest score window) and a score.

### 10.2 L2 envelopes

L2 = within-L1 sub-envelopes. Detected by multipass quadrant validation
on the nn40-smoothed sim_mat. Anchor-shift correction baked in (`+ceiling((G+1)/2)`
shift at write time). L2 envelopes are nested inside L1.

### 10.3 L3 sub-blocks (in the atlas)

L3 = candidate-level. Originally an R-side step but now driven by the
atlas's L3 contingency UI. The user can:
- Promote an L2 to a candidate ("promote" button)
- Refine candidate boundaries via window-mode drafts (turn 106)
- Merge or separate adjacent slabs based on contingency table verdict

The atlas's "candidate" object is the L3-level artifact. It has:
```
candidate = {
  candidate_id, chrom,
  start_w, end_w,            // window indices
  start_bp, end_bp,           // genomic coords
  K, K_used,                  // 3 | 6 | 8
  source: 'L1' | 'L2' | 'l3_draft_w' | ...,
  resolution: 'L1' | 'L2' | 'W',
  parent_l2,                  // if promoted from an L2
  ...
}
```

### 10.4 Layered detection philosophy

The hierarchy is intentional: L1 finds gross structure, L2 refines, L3
commits. Each level over-splits relative to the truth (false positives
cheap, false negatives not). The L3 contingency table validation is
where over-splits get merged back if the data says they should.

---

## 11. Heterozygosity histogram and het-shape sub-module

**Defined**: turn 98 in atlas, `_renderDetailedHetHistogram`.
**Concept**: For a candidate, plot the distribution of per-sample
heterozygosity. Three peaks (low/mid/high) indicate a clean REF/HET/INV
structure. Single peak indicates monomorphism. Bimodal indicates
something else.

### 11.1 What's histogrammed

For each sample, compute per-sample heterozygosity inside the candidate.
The per-sample value is currently atlas-derived (proxied from PC1
spread; not a true het count) when the `per_sample_het` layer isn't
loaded. With the layer, it's the real per-site het rate from genotypes.

Source flag: `'precomp_het'` (real layer) or `'atlas_proxy'` (PC1 spread).

### 11.2 Bin shape

16 bins, range 0..1 by default. Visualization shows histogram bars +
mode markers if K=2 detected.

### 11.3 K=2 detection

`detectBimodalHet()` (related to `classifyDetailedCandidate` from turn 90):
- Find local maxima
- Test gap between modes (default `_DBD_GAP_MIN = 0.15`)
- Test z-score separation (default `_DBD_Z_MIN = 2.0`)
- Min 3 samples per mode
- Thresholds: `_DBD_LOW_HET_THRESHOLD = 0.40`, `_DBD_HIGH_HET_THRESHOLD = 0.50`

If bimodal: classifies the candidate as having a heterozygous mid-band
(REF/INV system) vs a homozygous-only structure.

### 11.4 Why this matters for grouping

Three bands on PCA could be three haplotype classes (H1/H2/H3 all
homozygous-ish) OR a REF/HET/INV system (where the middle band is
truly heterozygous). Heterozygosity tells them apart:
- REF/HET/INV: middle band has HIGH het (sites segregating)
- H1/H1, H2/H2, H3/H3: all bands have LOW het if each H is a single
  haplotype background

So het + PCA together do what neither alone can do: distinguish
"3 homozygous classes" from "REF/HET/INV system."

---

## 12. SNP density confounder visualization

**Defined**: turn 95, 3-mode toggle.
**Concept**: SNP density per window varies by 10× along chromosomes.
Bands look more or less spread depending on how many SNPs the local PCA
had. The atlas surfaces this as a confounder.

### 12.1 Modes

- **off**: don't show SNP density
- **strip**: thin colored strip above the lines panel showing SNP-density
  contrast
- **shade**: shade the lines-panel background by SNP density (light =
  low density)

### 12.2 Data source priority

```
n_snps from precomp:n_snps    →    actual SNP count per window
fallback: precomp:snp_count
fallback: precomp:lambda1     →    PCA eigenvalue 1 (correlates w/ informativeness)
fallback: precomp:variance_pc1
fallback: atlas-side proxy from PC1 spread distribution
```

### 12.3 Why this matters for regime calls

Two windows can have the same band structure but very different SNP
density:
- **Window A**: 200 SNPs, clear 3-band PCA
- **Window B**: 20 SNPs, 3 bands but very fuzzy

If you call regime change between A and B based on band fuzziness, you're
calling SNP-density confounder, not real biology.

The 3-mode toggle lets the user check: does the regime breadth track SNP
density? If yes, the regime call is suspect.

---

## 13. Q-ancestry display

The atlas has an ancestry page (page 9) that displays Q matrices when
loaded. This is a DISPLAY feature, not an analysis feature.

### 13.1 Layers consumed

```
ancestry_window: { start_bp, samples × K Q matrix }
ancestry_sample: { samples list, K, q_matrix }
ancestry_q_means: { q_means by region }
ancestry_q_global: per-sample chromosome-wide Q
ancestry_q_chrom: per-chromosome per-sample Q
snp_q_support: per-SNP Q-support
```

### 13.2 What the atlas does NOT do with Q

The atlas does NOT classify per-sample ancestry structure into
single_block / two_block_composite / gradient / fragmented (which is
what phase 4b.3 does). It just paints Q values onto the screen.

### 13.3 Bridge to phase 4b.3

If the atlas were to add the four-category classification, it would do
so by:
1. Reading `ancestry_window.q_matrix[w, sample, k]`
2. For each sample, walking across w and looking at how Q changes
3. Classifying: single_block (Q nearly constant), two_block (one
   discrete jump), gradient (smooth slide), fragmented (many small jumps)

This is straightforward to add. ~100 lines. Not currently in the atlas.

---

## 14. Karyotype label vocabulary

**Defined**: turn 87, `labelVocab` system.
**Concept**: Map K-means cluster ids to human-readable labels. Two
vocabularies: legacy and detailed.

### 14.1 Legacy vocab (K=3)

`g0 → "REF"`, `g1 → "HET"`, `g2 → "INV"`. The classic shorthand.

### 14.2 Detailed vocab (K=6)

Partial nesting: `g0 → "g0a"`, `g1 → "g1a"`, `g2 → "g1b"`, `g3 → "g0a"`,
`g4 → "g1c"`, `g5 → "g0b"`, ... — the actual mapping reflects the
sub-band structure within a 3-band frame. Documented in `KARYO_LABEL_DETAILED`.

### 14.3 What this is NOT

It's NOT a haplotype-class assignment (H1/H2/H3). It's a way to label
clusters. The cluster IS the unit; the label is just text.

To go from clusters → haplotype class, an external decision is needed:
- "g0 in candidate LG28_cand_01 = H1/H1"
- "g1a = H1/H2"
- ...

This decision currently lives nowhere in the atlas. It's the missing
annotation surface (called out in the architectural-friction handoff
section). When the atlas builds the per-candidate annotation surface
described in the handoff (a small future turn), the haplotype-class
labels go there.

### 14.4 Persistence

`state.labelVocab = 'legacy' | 'detailed'`. Persists in localStorage as
`'inversion_atlas.labelVocab'`. Default 'legacy'.

---

## 15. What the atlas does NOT compute (and why)

Listed for the next chat's overlap analysis.

| Phase 4b feature | Why atlas can't do it |
|---|---|
| **SV-anchor seeded K-means** (cheat 1) | Atlas has no access to flashlight/sv_prior data. Would need a new layer JSON. |
| **Het-DEL constrained clustering** (cheat 2) | Same — needs internal_dels from sv_prior. |
| **Clair3 phase switches** (signal 2 of 4b.2) | Needs per-sample read-level haplotype info from Clair3 postprocess. Not in atlas. |
| **Flashlight hemizygous segments** (signal 3 of 4b.2) | Same as above — sv_prior cache is HPC-side. |
| **Engine B local-Q ancestry** (4b.3 input) | Atlas displays Q if loaded but doesn't run NGSadmix or instant_q.cpp. |
| **GC vs DCO classification** (cheat24) | Needs the multi-signal fusion of 4b.2. Atlas has diamond detection (visual analog) but doesn't classify event types. |
| **Registry commit** (4b.4 seal) | Atlas writes to its own state, not to the canonical sample_registry. |
| **q6_ flat keys** (q6_group_validation, q6_validation_promotion_cap, etc.) | Same — registry-side machinery. |
| **Bootstrap CI on breakpoints** (STEP_BP_02) | Atlas displays breakpoints but doesn't compute CIs from BEAGLE dosage. |
| **Per-sample ancestral fragment scan** | Same — needs per-SNP per-sample dosage, atlas only has per-window per-sample PC1. |

The pattern: **atlas can do everything that operates on per-window per-sample
PC1/PC2 loadings**. It cannot do anything that needs per-SNP per-sample
data or external signal layers (Clair3, flashlight, Engine B local-Q,
BCFtools).

---

## 16. Atlas state model: what's persisted, what's transient

Critical for the export contract design.

### 16.1 Persisted across sessions (localStorage / IndexedDB)

- `state.labelVocab` (legacy/detailed)
- `state.k` (3 or 6, "K mode")
- `state.activeMode` (default/detailed)
- View preferences (lcMode, trailN, flipPC1, etc.)
- Repeat density preferences (yMode, viewMode, loess bands, customSpan,
  yModeUserPicked) — turn 111
- Recent files history — turn 110
- Catalogue diamond mode — turn 93
- IndexedDB-backed JSON cache — turn 73f, auto-restores loaded JSONs

### 16.2 Computed at runtime (not persisted)

- `getL2Cluster(idx)` — K-means cluster results, lazily cached in memory
- `detectDiamonds(candidate)` — recomputed on demand
- `sigmaProfileL2(idx, K)` — recomputed on demand
- `computeStructuralHaplotypeTransitionGraph()` — recomputed on demand
- `computeBandReachAcrossL2s(indices)` — recomputed on demand

### 16.3 In `state.data` (loaded from precomp JSON)

The whole §2 data model. Persists as long as the JSON is loaded. The
IndexedDB cache restores it on reload.

### 16.4 Candidates registry

`state.candidates` (default mode) and `state.candidates_detailed`
(detailed mode). Both arrays of candidate objects. Per-session, kept in
memory. Could be persisted via export.

---

## 17. Possible export shapes for the atlas → pipeline contract

Recorded for the next chat, NOT a final design. This is what the atlas
COULD export today if a button were added.

### 17.1 Per-candidate export

```jsonc
{
  "format_version": "atlas_candidate_export_v1",
  "atlas_provenance": {
    "atlas_version": "v4.t113",
    "exported_at": "2026-05-01T...Z",
    "session_state_hash": "..."
  },
  "cohort": {
    "n_samples": 226,
    "species": "Clarias gariepinus",
    "chrom": "C_gar_LG28",
    "n_windows": 4302,
    "smoothing": "nn40",
    ...
  },
  "candidate": {
    "candidate_id": "LG28_cand_01",
    "chrom": "C_gar_LG28",
    "start_w": 3260, "end_w": 3870,
    "start_bp": 15115000, "end_bp": 18005000,
    "K_used": 6,
    "source": "l3_draft_w",
    "parent_l2": 14,
    // Atlas's clustering output:
    "atlas_band_assignments": {
      "per_sample": {
        "sample_001": { "primary_band": "g0", "band_reach": 1, "stable": true },
        "sample_002": { "primary_band": "g1a", "band_reach": 2, "stable": false },
        ...
      },
      "per_window_per_sample": {                  // optional, large
        "sample_001": [0,0,0,0,1,0,0,...],         // band index per window
        ...
      }
    },
    // User's annotations (currently MISSING — need new UI):
    "haplotype_labels": {
      "g0":  "H1/H1",
      "g1a": "H1/H2",
      "g1b": "H1/H3",
      "g2":  "H2/H2",
      "g3":  "H2/H3",
      "g4":  "H3/H3"
    },
    "annotator_notes": "Diamond detected at 16.2 Mb, suggesting two stacked systems",
    // Atlas-computed analytics:
    "diamond_detection": {
      "n_diamonds": 1,
      "n_strict": 1,
      "diamonds": [
        { "splitting_band": 1, "start_w": 3401, "end_w": 3415, "stable_bands": [0, 2] }
      ]
    },
    "sigma_profile": {
      "verdict": "TWO_INVERSIONS",
      "q50": 0.034,
      "ratio_high": 0.02,
      "bimodality_coef": 0.43
    },
    "regime_breadth": {
      "value": "narrow",
      "narrow_fraction": 0.81,
      "bands_populated": 3
    },
    "transition_boundaries": [
      {
        "position_mb": 15.115,
        "transition_rate": 0.42,
        "is_hotspot": true,
        "edges": [{"from_band": 0, "to_band": 1, "count": 80}, ...]
      },
      {
        "position_mb": 18.005,
        "transition_rate": 0.38,
        "is_hotspot": true,
        "edges": [...]
      }
    ],
    // Inferred sample groups (for phase 4b.4 seal compatibility):
    "sample_groups": {
      "HOM_REF": ["sample_007", "sample_023", ...],   // band g0 samples
      "HET":     ["sample_011", "sample_045", ...],   // band g1 samples
      "HOM_INV": ["sample_002", "sample_094", ...],   // band g2 samples
      "RECOMBINANT": ["sample_113"],                  // band-reach > 2 OR diamond-edge
      // Multi-haplotype case (K=6):
      "H1_H1": [...], "H1_H2": [...], "H1_H3": [...],
      "H2_H2": [...], "H2_H3": [...], "H3_H3": [...]
    }
  }
}
```

### 17.2 What the pipeline reader would do

`STEP_BP_00_read_atlas_export.R` would:
1. Parse the JSON
2. Register the candidate via `reg$intervals$add_candidate()`
3. Register the sample groups via `reg$samples$register_groups_for_candidate()`
4. Optionally write atlas's analytics (diamond, sigma, transition rates)
   as evidence blocks for downstream phases
5. Hand off to STEP_BP_01

The methodology question (binary vs multi-haplotype) determines
whether `sample_groups` uses `HOM_REF/HET/HOM_INV` keys (K=3 path) or
the multi-haplotype keys (K=6 path). The export JSON supports both.

---

## End

This document captures atlas-side computation in enough detail that the
next chat, given phase 4b's README + scripts, can:

- Identify per-script overlap with §3-§14 above
- Assess what's redundant (atlas already does it)
- Assess what's missing (atlas can't do it)
- Plan which logic to bring across the seam
- Design the export JSON shape that carries everything the pipeline needs
- Estimate effort: ~50 lines per signal to wire into existing pipeline,
  vs full rebuild

Pair this with `HANDOFF_2026-05-01_session_summary.md` (which has the
project context and the multi-haplotype methodology question) and the
`phase_7_karyotype_groups/proposal/README.md` (which has the pipeline
architecture).
