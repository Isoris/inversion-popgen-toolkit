# SPEC — Observable-allele H-label system, empirical-first design

**Status**: drafted turn 137, rewritten turn 138 to flip the ordering
(observe first → derive system size → label, instead of assume K →
label → check). Not yet implemented. Replaces the implicit H-label
semantics in `_KARYO_DETAILED_LABELS` (`Inversion_atlas.html` line
~32940).

**Working file**: `Inversion_atlas.html` line ~32940 area
(`_KARYO_DETAILED_LABELS`, `getKaryotypeLabel`, `getKaryotypeLabelCaveat`).

**Trigger** (Quentin, turn 137-138):
> *"H would be the shorthand for Haplotype of an observable allele of
> a system 1/2/3 [...] maybe we never observe AA smth like that."*
>
> *"We should just observe empirically and change the model later it
> make more sense to me since its a finite set of states. [...] First
> separate Homozygous from heterozygous bands then based on how many
> total heterozygous bands we can kind of know this loci seems to have
> max 2 inversion systems because its 2 observable het bands."*

The existing atlas treats `H1/H1`, `H1/H2`, `H2/H2`, etc. as
diploid-class labels with a fixed K-indexed lookup table. This spec
inverts that ordering: **classify each band as HOM / HET / AMBIGUOUS
first**, then use the (`n_hom_bands`, `n_het_bands`) pair to derive
the inferred allele system size, **then** label. Empirical observation
drives model fitting, not the other way around.

---

## 1. The combinatorial model — why hom/het counts encode system size

With `H` observable allele clusters at a candidate locus, the diploid
state space is fully determined:

| H (observable alleles) | Hom classes | Het classes | Total bands at full obs |
|---|---|---|---|
| 1 | 1 (H1/H1) | 0 | 1 |
| 2 | 2 (H1/H1, H2/H2) | 1 (H1/H2) | 3 |
| 3 | 3 | 3 (H1/H2, H1/H3, H2/H3) | 6 |
| 4 | 4 | 6 | 10 |

Number of observable alleles → number of inversion systems segregating:

- 1 allele → no inversion polymorphism observed
- 2 alleles → 1 inversion system
- 3 alleles → 2 inversion systems
- 4+ alleles → 3+ systems (well beyond the K=6 atlas ceiling and
  unrealistic for hatchery timescales per Porubsky 2022 / Corbett-Detig
  2012; mark as suspect)

**The het-band count is a more reliable system-size indicator than
the hom-band count** in finite cohorts because:

- Missing hom bands are common: the founder lineage carrying H1/H1
  may have been lost or never founded. K-means will simply not
  produce a band where there are no samples to populate it.
- Missing het bands also happen but follow a different logic: if H1
  and H3 alleles both segregate but no H1×H3 cross was ever sampled,
  H1/H3 is missing.
- Het bands have larger expected populations under random mating
  (Hardy-Weinberg) and so are more likely to be sampled when present.

So `n_het_bands ≥ 1` is strong evidence of inversion polymorphism,
and the het count narrows the system size:

| n_het_bands observed | Implied min H | Implied min systems | Notes |
|---|---|---|---|
| 0 | 1 | 0 | No polymorphism, OR all samples in one class (unlikely in a cohort of 226) |
| 1 | 2 | 1 | Standard biallelic inversion locus |
| 2 | 3 | 2 | 3 alleles segregating but ONE het cross unobserved (sample size limitation OR rare cross) |
| 3 | 3 | 2 | 3 alleles, all crosses observed |
| 4-6 | 4 | 3 | Beyond atlas ceiling; flag as suspect |

The spec does NOT claim the inferred H is correct. It claims it's the
**minimum H consistent with the observed het bands**. The actual
underlying system could be larger (with more unobserved classes).

### 1.1 Why empirical-first matters

The existing `_KARYO_DETAILED_LABELS[K]` table hard-codes assumptions
like "K=4 means H1/H1, H1/H2, H2/H2, H1/H3" — i.e. one het class plus
one extra-far hom band. That's *one* possible 4-band regime, but the
actual data could equally show:

- `[HOM, HET, HOM, HET]` → 2 hom + 2 het → suggests 3 alleles with
  H3/H3 missing
- `[HET, HOM, HOM, HOM]` → 3 hom + 1 het → suggests something
  unusual, possibly K-means oversplit a homozygous class (sub-haplotype
  divergence within H1 or H2)
- `[HOM, HET, HET, HOM]` → 2 hom + 2 het → similar to first case,
  same inferred regime
- `[HOM, HET, AMBIGUOUS, HOM]` → 2 hom + 1 het + 1 ambiguous → can't
  call the system size cleanly; flag for manual review

Each of these has different biological interpretation. A fixed table
hides the difference; the empirical-first spec surfaces it.

---

## 2. Per-band classification

For each band `b ∈ [0, K)` of `candidate.locked_labels`, compute the
dosage-HET fraction:

```
members_b      = { si : locked_labels[si] == b }
het_called_b   = { si ∈ members_b : dosage[si] == HET }
non_na_b       = { si ∈ members_b : dosage[si] != NA }
het_fraction_b = |het_called_b| / |non_na_b|       (NaN if non_na_b == 0)
```

Classification with two thresholds:

| Condition | Classification |
|---|---|
| `het_fraction_b >= _HLABEL_HET_FRACTION_HIGH` (default 0.70) | **HET** |
| `het_fraction_b <  _HLABEL_HET_FRACTION_LOW`  (default 0.30) | **HOM** |
| `_HLABEL_HET_FRACTION_LOW <= het_fraction_b < _HLABEL_HET_FRACTION_HIGH` | **AMBIGUOUS** |
| `het_fraction_b == NaN` (no dosage data at all in this band's samples) | **NO_DOSAGE** |

Also defensible: too few non-NA samples (< `_HLABEL_MIN_DOSAGE_N`,
default 5) → also **NO_DOSAGE** even if het_fraction is technically
defined. A 0/1 split on a single sample is statistically meaningless.

Asymmetric thresholds because dosage HET is noisier than dosage HOM
in low-coverage regions; HET clusters are expected to be cleaner than
HOM clusters in well-resolved inversions. Defaults are seeds — the
threshold tuning is a Slice 1 calibration step (§9).

---

## 3. Inferred regime from band counts

After classification, count:

```
n_hom_bands       = bands classified HOM
n_het_bands       = bands classified HET
n_ambiguous       = bands classified AMBIGUOUS
n_no_dosage       = bands classified NO_DOSAGE
```

Then derive:

### 3.1 Inferred minimum allele count

```
implied_H_from_het_count =
   { 0:1, 1:2, 2:3, 3:3, 4:4, 5:4, 6:4 }[n_het_bands]   // see §1 table
```

When `n_het_bands == 2`, we can't tell from het count alone whether
H=3 with one missing cross OR H=2 with band-merge (very unusual). The
table chooses H=3 because that's the more biologically plausible
interpretation in hatchery data.

### 3.2 Consistency check against hom band count

If `n_hom_bands` and `implied_H_from_het_count` are consistent:

- `n_hom_bands == implied_H` → fully observed: all H homozygous
  classes present
- `n_hom_bands == implied_H - 1` → one homozygous class missing
  (e.g., founder lineage H1 lost → no H1/H1 carriers)
- `n_hom_bands == implied_H - 2` → two missing
- `n_hom_bands < implied_H - 2` → suspect: more missing than
  expected, OR the het classification is wrong, OR K-means didn't
  resolve all classes

Consistency outcomes:

| n_hom_bands vs implied_H | regime_consistency |
|---|---|
| equal | `clean` |
| -1 from implied_H | `partial_obs_minor` (one hom class missing) |
| -2 from implied_H | `partial_obs_major` (two hom classes missing) |
| more negative | `suspect` |
| > implied_H | `suspect_extra_hom` (probable K-means oversplit) |

### 3.3 Inferred system count

```
n_inferred_systems = max(0, implied_H - 1)
```

This is the **minimum** number of segregating inversion events at this
locus consistent with the het pattern. The actual count could be
higher (if multiple het classes merged into one band by K-means
under-resolution).

---

## 4. H-label assignment, only AFTER counts are known

Once `n_hom_bands`, `n_het_bands`, `implied_H` are settled,
assign labels:

### 4.1 Sort hom bands by median PC1

The HOM-classified bands get the canonical `H1/H1`, `H2/H2`, ...
labels in PC1 order:

```
hom_bands_sorted_pc1[0]   → label "H1/H1"
hom_bands_sorted_pc1[1]   → label "H2/H2"
hom_bands_sorted_pc1[2]   → label "H3/H3"
```

When `n_hom_bands == implied_H` (consistency = `clean`), this is
unambiguous.

When `n_hom_bands < implied_H` (`partial_obs_*`), we need to figure
out which hom class is missing. The het bands' median PC1s constrain
this: a HET band labelled `Hi/Hj` should sit between the median PC1s
of its two parent hom classes. So:

1. Sort observed hom bands by PC1: `pc1_homs = [p_a, p_b]` (for 2
   observed homs)
2. Sort observed het bands by PC1: `pc1_hets = [h_1, h_2, ...]`
3. If any het has PC1 below `p_a`: there's a missing hom below it
   (label that as H1/H1 of the implied system; observed homs become
   H2/H2 and H3/H3)
4. If any het has PC1 above `p_b`: missing hom is above (observed
   homs are H1/H1 and H2/H2; missing is H3/H3)
5. If all hets are between `p_a` and `p_b` AND `implied_H == 3`: the
   missing hom could be either end, can't determine — fall back to
   ordinal labels for this candidate, mark `regime_consistency` as
   `partial_obs_unconstrained`.

When the inference is constrained (cases 3 or 4), output explicit
labels including the missing class in `missing_hom_classes`.

### 4.2 Het band labels by PC1 interpolation

Once hom band labels are pinned (or pinned with a known missing
position), each HET-classified band sits at some median PC1. The
label `Hi/Hj` is determined by which two hom classes (real or
known-missing) the het band falls between in PC1 order:

- Het band PC1 between H1/H1 and H2/H2 medians → label `H1/H2`
- Between H2/H2 and H3/H3 → `H2/H3`
- Below H1/H1 → suspect (het band shouldn't sit below the lowest
  hom); label `band N*` and add to `ambiguous_bands`
- Above the highest hom → similar suspect handling

For the case where one het sits between H1/H1 and H3/H3 but skips
over (above) H2/H2's PC1, label as `H1/H3` (the cross-pair that
spans the missing middle) — this is a valid configuration in the
3-allele system.

### 4.3 Ambiguous + no-dosage bands

These don't get H-labels. They get fallback labels:

- AMBIGUOUS → `band {N+1}*` with classification reason in tooltip
  ("dosage-HET fraction 0.52 — between thresholds")
- NO_DOSAGE → `band {N+1}` with caveat ("no dosage data — falling
  back to ordinal label")

---

## 5. Output JSON shape

Per-candidate annotation:

```json
{
  "candidate_id": "LG28_cand_03",
  "K": 4,
  "vocab_mode": "observable_allele_v1",

  "thresholds": {
    "het_fraction_high": 0.70,
    "het_fraction_low":  0.30,
    "min_dosage_n":      5
  },

  "band_counts": {
    "n_hom":       2,
    "n_het":       2,
    "n_ambiguous": 0,
    "n_no_dosage": 0,
    "n_total":     4
  },

  "implied_regime": {
    "implied_H":            3,
    "n_inferred_systems":   2,
    "consistency":          "partial_obs_minor",
    "missing_hom_classes":  ["H3/H3"],
    "missing_het_classes":  []
  },

  "bands": [
    {
      "band_idx":           0,
      "label":              "H1/H1",
      "classification":     "HOM",
      "median_pc1":         -2.31,
      "n_members":          84,
      "n_dosage_het":       3,
      "n_dosage_non_na":    82,
      "het_fraction":       0.037
    },
    {
      "band_idx":           1,
      "label":              "H1/H2",
      "classification":     "HET",
      "median_pc1":         0.05,
      "n_members":          99,
      "n_dosage_het":       88,
      "n_dosage_non_na":    95,
      "het_fraction":       0.926
    },
    {
      "band_idx":           2,
      "label":              "H2/H2",
      "classification":     "HOM",
      "median_pc1":         1.87,
      "n_members":          35,
      "n_dosage_het":       1,
      "n_dosage_non_na":    33,
      "het_fraction":       0.030
    },
    {
      "band_idx":           3,
      "label":              "H1/H3",
      "classification":     "HET",
      "median_pc1":         2.95,
      "n_members":          8,
      "n_dosage_het":       7,
      "n_dosage_non_na":    8,
      "het_fraction":       0.875
    }
  ],

  "per_sample_verdicts": [
    { "sample_idx": 17, "sample_id": "Cga_017", "band_idx": 1,
      "label": "H1/H2", "dosage": "HET", "verdict": "AGREE" },
    { "sample_idx": 42, "sample_id": "Cga_042", "band_idx": 0,
      "label": "H1/H1", "dosage": "HET", "verdict": "MISMATCH_FALSE_HOM" }
  ],

  "summary": {
    "n_agree":               198,
    "n_mismatch":            12,
    "n_no_call":             16,
    "agree_fraction":        0.943
  }
}
```

Note `implied_regime` is now first-class output. Atlas surfaces use it
to display the regime summary even when full labels can't be assigned.

---

## 6. Slice plan — empirical-first ordering

This is the biggest change from the turn-137 draft. Slice 1 ships
**counting + classification only**, no labels. Quentin runs it on real
data, looks at the histogram of `(n_hom, n_het)` pairs across all
candidates, and tells us what regimes actually appear. Only after that
empirical pass does Slice 2 add labelling using the empirical regime
distribution as a constraint.

### Slice 1 — classification + regime inference, NO labels (~1 turn)

Atlas Slice 1 ships:

- [ ] `_classifyHLabelBands(candidate, state)` — pure JS, performs §2
      classification and §3 regime inference. Sets
      `candidate.h_band_classifications`, `candidate.h_band_counts`,
      `candidate.h_implied_regime`.
- [ ] G-panel karyotype tab additions:
  - [ ] Per-band rows now show classification chip
        (`HOM` / `HET` / `AMBIGUOUS` / `NO_DOSAGE`) with dosage HET
        fraction tooltip
  - [ ] Header gains a regime summary line:
        *"K=4 · 2 hom + 2 het · implies H=3 (2 systems) ·
        partial_obs_minor (H3/H3 missing) · 94.3% agreement"*
  - [ ] No H-labels yet — bands keep their ordinal labels (`band 1
        (lo)`, etc.) for now
- [ ] Tests cover classification edge cases + regime-derivation table.
- [ ] **Calibration output**: a `[ Export regime histogram TSV ]`
      button on the karyotype tab footer that emits one row per saved
      candidate with `(candidate_id, K, n_hom, n_het, n_ambiguous,
      implied_H, n_systems, consistency)`. Quentin runs this on
      LG28+LG12+a few more, eyeballs the histogram, then we can
      decide what regime list to support.

### Slice 2 — H-label assignment from empirical regimes (~1 turn)

After Quentin reports the regime histogram, Slice 2 ships:

- [ ] Empirical regime list defined as a constant (e.g.,
      `[(2 hom + 1 het), (2 hom + 2 het), (3 hom + 3 het), ...]`)
      based on what actually appears
- [ ] `_assignHLabels(candidate)` — §4.1 + §4.2 logic, with
      empirical-regime guard: refuse to label combinations not in the
      observed regime list, fall through to ordinal labels with a
      caveat
- [ ] `getKaryotypeLabel(bandIdx, K, mode, candidate)` extended
      signature — per-candidate `h_labels` override
- [ ] G-panel karyotype tab now shows H-labels (replacing ordinal
      labels) when the candidate's regime is in the empirical list
- [ ] Caveat string update — three states (legacy / detailed-with-
      h_labels / detailed-fallback) per turn-137 §5.2

### Slice 3 — page-1 PCA mismatch highlighting (~0.5 turn)

When user clicks a "⚠ flag samples that disagree" action in the
G-panel, paint MISMATCH samples on page-1 PCA with a distinct ring
(e.g., red 2px outline) without changing the band-color fill. State
slot: `state.hlMismatchHighlight: bool`. Repaint hook in `drawPCA`.

### Slice 4 — inheritance tab H-label propagation (~0.5 turn)

Folds into G-panel Slice 3. When inheritance group breakdown lists
"cand_42 → g1 (n=80)", change `g1` to `H1/H2` (or whatever the
candidate's per-band label is). Pure rendering update once Slice 2
emits per-candidate `h_labels`.

### Slice 5 — R-side STEP_HL01 (separate spec, ~1 turn)

A separate spec (`SPEC_STEP_HL01_h_label_classifier.md`) for the
LANTA-side classifier that ingests dosage_chunks JSON + candidate
atlas TSV and emits `<chrom>_h_labels.json` matching §5 output shape.
Ships when the dosage_chunks R-side pipeline is producing the JSON
layer reliably.

---

## 7. Why the ordering change matters

Concrete example from the existing atlas's K=4 fixed table assumption:

**Old (turn-137) flow**: candidate has K=4, locked_labels split into
4 bands. Apply `_KARYO_DETAILED_LABELS[4]`, label them
`[H1/H1, H1/H2, H2/H2, H1/H3]`. Use dosage to validate. If band 0
turns out to be HET-classified (high dosage HET fraction) instead of
HOM, we have a labelling-vs-classification conflict — what do we do?
Re-label? Trust the table?

**New (turn-138) flow**: candidate has K=4. Classify bands first:
maybe `[HET, HOM, HOM, HET]`. Count: 2 hom + 2 het → implies H=3,
consistency partial_obs_minor (H3/H3 missing — assuming hom band PC1s
are below the highest het band PC1, otherwise H1/H1 missing).
THEN label by PC1: hom bands are H1/H1 and H2/H2 (or similar based
on which one is missing), het bands at PC1 between them is H1/H2,
het band at higher PC1 is H1/H3.

The new flow has an answer. The old flow has a conflict.

---

## 8. The "we figured this out before" question

Quentin's intuition that "max 2-3 overlapping inversions" is
biologically right is well-supported but is *the model* not *the
data*. In any specific candidate from the actual 226-sample data,
the EFFECTIVE number of segregating systems could be 0, 1, or 2
(rarely 3) depending on which alleles got drawn from the founder
population.

The empirical histogram from Slice 1 will tell us:

- If the modal regime is `1 hom + 1 het` (i.e., K=2 effectively): the
  cohort is dominated by 1-allele or near-fixed loci. Inversion
  polymorphism is rare.
- If the modal regime is `2 hom + 1 het`: standard biallelic
  inversion polymorphism (1 system per locus is the norm).
- If the modal regime is `2 hom + 2 het` or `3 hom + 3 het`: some
  loci have multi-system stacking — 2 inversions overlapping.
- If `4+ het` shows up at any locus: extraordinary, deserves a
  dedicated investigation (probably K-means oversplit; flag for
  manual review).

This is the empirical baseline that the manuscript Methods can cite:
*"Across N candidate loci, X% had 1 het band (single inversion
system), Y% had 2 het bands (3-allele system, possibly 2 inversions),
and Z% had no resolvable polymorphism."*

---

## 9. Open questions

1. **Threshold defaults**: 0.70 / 0.30 / min_n=5 are seeds.
   Calibration via Slice 1's regime histogram on LG28+LG12.
2. **Caching**: cache on candidate object, invalidate when dosage
   layer reloads (mirrors `_inheritanceCacheKey` pattern).
3. **Multi-track candidates** (page 4 dual-track, rare): per-track
   classification — simpler; tracks are independent inversions on the
   same locus by construction.
4. **Asterisk vs explicit "uncertain" UI**: asterisk + tooltip
   explaining the classification reason. Less visual noise than an
   explicit "?" affordance.
5. **`_KARYO_DETAILED_LABELS` legacy path**: keep as fallback when no
   dosage data is loaded — used in demo scenarios and partial-data
   sessions. Don't deprecate.
6. **Strict vs advisory missing-class identification**: when the data
   doesn't unambiguously identify which hom class is missing (all het
   bands sit between the observed homs, leaving the missing one's
   PC1 position unconstrained), fall back to ordinal labels with a
   `partial_obs_unconstrained` consistency flag. The implied regime
   summary still reports the count-based inference; only the H-label
   assignment falls back.

---

## 10. Honest framing

**What this spec accomplishes**:

- **Empirical observation drives the model**, not the other way
  around. Slice 1 just classifies and counts; Slice 2 labels using
  whatever regimes the data actually presents.
- The labelling system handles missing classes as a first-class
  output — `missing_hom_classes`, `missing_het_classes` — instead of
  silently mis-labelling a 4-band candidate as if it were a fully
  observed K=4 system.
- The `implied_regime` block surfaces the inference: H, n_systems,
  consistency. The atlas can show "K=4 · 2 hom + 2 het · implies
  H=3 (2 systems) · H3/H3 missing" without committing to specific
  H-labels until Slice 2 ships.
- Counts of (n_hom, n_het) pairs across the whole catalogue produce
  a manuscript-grade empirical regime distribution.

**What it deliberately doesn't do**:

- Doesn't try to identify WHICH inversion event corresponds to which
  H-letter. H1, H2, H3 stay observable-cluster labels.
- Doesn't model within-haplotype divergence (sub-haplotypes — H1 →
  H1', H1''). If K resolves a sub-haplotype, it shows up as an extra
  band that gets classified HOM (since the sub-haplotype carriers are
  homozygous for the inversion); the regime check then likely flags
  the candidate as `suspect_extra_hom`. That's the right behaviour
  — the user investigates.
- Doesn't autoassert biological regimes. H-labels are observables;
  combining them with GHSL haplotype divergence, het, F_ROH, family
  confounders is the regime framework's job.
- Doesn't ship labels in Slice 1. We wait for the empirical histogram
  before committing to a label scheme.

**What this means for the manuscript**:

Methods sentence (after Slice 1): *"Per-candidate band classifications
were computed by cross-referencing K-means cluster membership against
per-sample dosage HET/HOM calls; bands with dosage-HET fraction ≥
0.70 were classified as heterozygous classes, bands < 0.30 as
homozygous classes, and intermediate bands as low-confidence. Across
N candidate loci, the modal observed regime was X (% of candidates),
with multi-system loci (n_het ≥ 2) representing Y%."*

Methods sentence (after Slice 2): *"Heterozygous-class bands were
labelled in the H-system based on PC1 interpolation between adjacent
homozygous-class bands; missing classes (no observed homozygote for
H1, H2, or H3) were retained as first-class output rather than
back-filled by assumption."*

This grounds the labels in observable structure + dosage agreement,
explicitly accommodating missing classes — which is more honest than
either silently dropping them or asserting them via fixed-table
back-fill.

---

## 11. Companion specs

| Spec | Status | Relationship |
|---|---|---|
| `SPEC_g_panel_unified_groups.md` | Slices 1+2 shipped (turns 135-136) | Slice 1 of THIS spec adds classification chips + regime summary to the karyotype tab; Slice 2 adds H-labels |
| `SPEC_l3_het_dosage_coloring.md` | Implemented turn 129+ | Same dosage_chunks layer source, orthogonal use (color L3 dots vs label bands) |
| `SPEC_karyotype_per_interval_intersection.md` | Drafted turn 130 | Operates on GHSL karyotype calls (different signal); intersection produces the per-sample dosage state THIS spec consumes |
| `SPEC_multi_evidence_regime_framework.md` | Existing | This spec is one input layer; the regime framework synthesizes H-labels + GHSL + ROH + family-confounders into final regime calls |
| `SPEC_STEP_HL01_h_label_classifier.md` | NOT YET DRAFTED | R-side producer for the JSON layer this spec consumes (when dosage R-side is ready) |

---

## 12. Build estimate

| Slice | Effort | Gates on |
|---|---|---|
| 1 (classification + regime inference, NO labels) | ~1 turn | dosage_chunks layer loaded |
| 2 (H-label assignment using empirical regime list) | ~1 turn | Slice 1 regime histogram from Quentin |
| 3 (PCA mismatch highlighting) | ~0.5 turn | Slice 1 |
| 4 (inheritance tab H-label propagation) | ~0.5 turn | Slice 2 + G-panel Slice 3 |
| 5 (R-side STEP_HL01) | separate spec, ~1 turn | dosage R-side pipeline ready |

Total atlas-side: ~3 turns across 4 slices, gated on Quentin's
empirical regime histogram between Slices 1 and 2.

---

**End of spec.** Pending Quentin's review before promotion to
`specs_todo/`.
