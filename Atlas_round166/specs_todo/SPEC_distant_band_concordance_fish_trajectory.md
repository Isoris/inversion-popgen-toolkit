# SPEC — Distant-band concordance / fish-trajectory inheritance system

**Status**: drafted turn 130, recovered from old chat
`f74cf5d4-1415-4fa7-8736-5781470af51e`. **Slice 1 + Slice 2 SHIPPED
turn 130.** Slices 3 + 4 still pending. Designed but never built before
this session.

**Trigger** (Quentin, original chat, paraphrased then verified verbatim):
> *"That group like 80% further in the genome it visits only 2 lines and
> not 5-6. So very distantly it tells us about the very distant inversion
> system."*
>
> *"No matter you put your river fish in the mixed hatchery they have
> their chromosome that come from 1 inversion system or the other so in
> the end we can still kind of 'know.'"*
>
> Turn 130 follow-up: *"It was an old chat when we said that when we
> select bands. Then much further in the genome the same group of fish
> visit many bands but sometimes it visit a subset of bands and we
> talked that its basically the inheritance the long range inheritance
> idk how to explain."*

The chat ended with me proposing "Turn 2 plan: build approach 2". The
turn never happened. This spec re-anchors the design and queues it for
implementation.

---

## 1. Why the existing inheritance pipeline isn't this

The atlas already ships `inheritanceGroupClustering()` (turn 117). That
function clusters **(candidate, band) pairs** by their fish-membership
masks via Jaccard distance. Two bands across two candidates that contain
mostly the same fish get merged into one inheritance group. Output:
group_ids that link bands across candidates (the I·g pills).

That's **the orthogonal axis** to what this spec describes:

| Existing (`inheritanceGroupClustering`) | This spec (distant-band concordance) |
|---|---|
| Inputs: user-promoted candidates with `locked_labels` | Inputs: every L2 envelope on the chromosome (or a chosen range) |
| Feature: per-band fish-membership mask (which fish are in band g of cand C?) | Feature: per-fish L2-band trajectory (which band is fish F in at L2 #1, #2, …) |
| Objective: which BANDS group together by shared fish | Objective: which FISH group together by shared band trajectory |
| Output: group_ids per band → I·g pills | Output: lineage_id per fish → coloring mode + co-membership matrix |
| Requires: ≥2 confirmed candidates | Requires: ≥3 L2s on the chromosome (always true) |

Both are valid; they tell you different things. You need both to do
what you've been describing. The current pipeline answers *"do bands at
distant candidates carry the same fish?"*; this spec answers *"do fish
at distant L2s land in the same band?"*.

## 2. The core computation (Approach 2 from the original chat)

Fixed n_L2 windows along the chromosome. Hungarian-aligned per-L2 K-means
labels (the atlas already computes these for band-reach §8). For every
fish pair (i, j), define:

```
concordance(i, j) = (# L2s where Hungarian-projected band of i == band of j) / n_valid_L2
```

This is a 226 × 226 symmetric matrix in [0, 1]. Diagonal = 1 by
definition. High off-diagonal = same lineage. Hierarchical clustering on
`1 - concordance` produces a per-fish lineage label.

Cost: O(n_samples² × n_L2). For 226 samples × 30 L2s = ~2M comparisons
per chromosome. Sub-second in JavaScript.

### 2.1 Hungarian projection

Atlas `alignLabels(labelsA, labelsB, K)` already exists. To project all
L2s into a common space, walk L2s in order and apply the running
Hungarian permutation:

```
ref = getL2Cluster(L2[0]).fixedKLabels
for i in 1..n_L2-1:
    cur = getL2Cluster(L2[i]).fixedKLabels
    perm = alignLabels(ref, cur, K).perm        // best perm of cur into ref space
    projected[i] = perm.applied(cur)            // perm[cur[s]] for each sample s
    ref = projected[i]                          // chain: next L2 aligns to this one
```

The chain breaks if K differs across L2s (some L2s might be K=2 effective
even when K=3 was requested). Handling: when the Hungarian agreement
score drops below 0.5, mark this L2 as a chain-break and start a fresh
chain for downstream L2s. Concordance becomes a chain-local statistic
per fish pair, then averaged across chains weighted by chain length.

### 2.2 Lineage clustering

Average-linkage agglomerative on `1 - concordance`. Cut threshold:
default 0.50 (i.e., fish pairs with ≥50% trajectory agreement land in
the same lineage). User-tunable in the UI later.

Output:

```
state.lineageResult = {
  n_samples,
  n_L2_chains,                                  // > 1 if Hungarian breaks
  concordance_matrix: Float32Array[n*n],        // upper-triangular packed?
  lineage_id_per_sample: Int32Array[n_samples], // 0..L-1
  n_lineages: int,
  dendrogram: { merges, heights },
  chain_meta: [ { L2_indices, n_L2 } ],
}
```

## 3. Combinatorial selection (Quentin's turn-130 addendum)

Original ask:
> *"When I select 2-3 bands sometimes if there are 6 lines in per sample
> lines then they visit different bands so it should not only be 1 by 1
> it could be combinatorial too."*

Resolution: lineage labels from §2 already encode trajectories. To
support "find fish that visit this combination of bands across these
L2s":

1. User lassos a region in the per-sample-lines panel covering N
   consecutive L2s, plus selects bands (say {0, 2}) at each.
2. Filter fish whose Hungarian-projected band at every selected L2 is in
   the chosen band set.
3. Show their lineage_id distribution (which lineages produced these
   trajectories?) plus their band membership at every other candidate
   on the chromosome.

This is a **derived query** on top of the matrix from §2, not a new
algorithm. So §2 ships first; combinatorial selection is a UI surface
on top.

## 4. Surfacing in the atlas

### 4.1 Per-sample-lines coloring mode

Add `state.linesColorMode === 'lineage'`. Each fish line painted in its
lineage's color (palette: golden-angle rotation, like the candidate-band
spec). Default off; user flips it on from the existing color-mode picker
in the lines header.

### 4.2 Co-membership matrix viewer

Modal triggered from a new toolbar button: shows the 226×226 concordance
matrix as a heatmap, ordered by hierarchical clustering. Pixel = fish
pair. Hover = "fish i / fish j: concordance %". Click a row → highlight
that fish's line in the per-sample-lines panel.

### 4.3 Lineage strip on per-sample-lines

A thin strip above the lines panel (sibling to the regime-breadth strip
and the I·g strip) showing **per-L2** which lineages dominate which
bands. Rendered as small colored pips per L2 per band.

## 5. Implementation slices

### Slice 1 — compute (~1 turn)  ← SHIPPED turn 130
- [x] `_hungarianChainProjection(L2_indices, K)` returns
      `{ projected_labels: Int8Array[n_L2 × n_samples], chain_breaks: int[] }`
- [x] `_concordanceMatrix(projected)` returns `Float32Array[n*n]`
- [x] `_lineageClustering(matrix, threshold)` returns lineage_id_per_sample
      via existing `_agglomerativeAverageLinkage` + `_cutDendrogram`
- [x] `runLineageCompute()` orchestrator + `state.lineageResult` + cache
      key (chrom + n_L2 + threshold)
- [x] Tests: synthetic 4-fish 4-L2 case where 2 fish always co-occur and
      2 others always co-occur — should produce 2 lineages.

### Slice 2 — coloring + lineage strip (~0.5 turn)  ← SHIPPED turn 130
- [x] `state.linesColorMode === 'lineage'` resolver in the existing
      `_resolveSampleColorByMode` + `_resolveSampleScopeColor` stubs
- [x] `_drawLineageStrip()` — per-L2 dominant-lineage colored bar above
      the per-sample-lines panel. Hatched grey for chain-breaks.
- [x] `state.linesLineageStripOn` slot (default true) + setter +
      localStorage persistence (`inversion_atlas.linesLineageStripOn`)
- [x] Toggle checkbox `linesLineageStripToggle` in lines header
- [x] Auto-trigger compute on first paint when `state.lineageResult`
      missing (mirrors the `_drawInheritanceLabelsStrip` lazy pattern)
- [x] Cache invalidation hooked into `applyData` chrom-switch path
- [x] Tests: source-level for picker entry + resolver + strip + state
      slot + setter + localStorage; behavioural for golden-angle hue
      math + edge cases (no result, invalid si, lineage_id -1)

### Slice 3 — co-membership matrix viewer (~0.5 turn)
- [ ] Modal canvas with reorderable matrix + dendrogram on the side
- [ ] Click-to-highlight in the lines panel

### Slice 4 — combinatorial selection UI (~1 turn, after Slices 1–3)
- [ ] Multi-L2 lasso on the lines panel
- [ ] Filter fish by combined band selection
- [ ] Show lineage distribution + cross-chromosome implications

## 6. Open design questions

1. **Threshold default**: 0.50 (>50% trajectory agreement). Genome-wide
   expectation in a hatchery with Ne~20 is unclear — Quentin's data on
   LG28 will tell us if 0.50 is too lax or too strict.
2. **Hungarian chain break threshold**: 0.5 agreement seems reasonable
   but is a guess. May need calibration against a chromosome where the
   ground truth is known.
3. **Lineage palette size**: golden-angle rotation works for ≤30
   lineages; beyond that the palette gets confusable. Cap at 30 with a
   "+N more" badge?
4. **Cross-chromosome lineages**: out of scope for Slices 1–4. A fish's
   lineage on LG12 and on LG28 are computed independently; matching
   them across chromosomes needs a different cohort-wide aggregator
   that's a separate spec.

## 7. What this is NOT

- **Not haplotype phasing.** This works at the per-window K-means band
  level, not at the SNP/marker level. The HPC vector engine
  (`ABABABB`/`010101`) does the marker-level version; this is the
  cheap atlas-side cousin.
- **Not a replacement for `inheritanceGroupClustering`.** That clusters
  bands by fish; this clusters fish by trajectories. Both are useful;
  both should run.
- **Not L2-sweep of the existing algorithm.** Those would rerun Jaccard
  on (every-L2-band, every-L2-band) — same axis as today, just more
  inputs. This is the orthogonal axis.
- **Not the staircase / sim_mat brute force.** That's boundary
  detection at the window-pair level. This operates on K-means labels.

## 8. Tests (Slice 1)

- Synthetic 226 × 30 trajectories with two planted lineages → cluster
  recovery within ε.
- Hungarian chain-break injection: two L2s with K=2-effective in the
  middle of a K=3 chain → chain split into 2 segments, concordance
  computed per-segment then weighted-averaged.
- Edge case: all fish identical trajectory → 1 lineage, concordance
  matrix all 1.
- Edge case: every fish unique trajectory → n_lineages == n_samples.
- Cache: re-running with same inputs returns same Float32Array reference.
