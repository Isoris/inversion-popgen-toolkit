# SPEC — Multi-evidence regime framework (K-means is geometric, not biological)

**Status**: drafted turn 130 final session. Doctrine introduced by
Quentin in chat `47cd29b9`. Some elements already in atlas (color modes,
per-band views), but the doctrine itself isn't formalized as a spec.

**Trigger** (Quentin, chat `47cd29b9`):
> *"In the scrubber, K-means bands should be treated as provisional
> geometric clusters, not final biological regimes. Add multiple
> interpretation modes: K-means mode (PCA-derived bands), Dosage mode
> (marker dosage pattern by band/regime), GHSL mode (haplotype-divergence
> support by band/regime), Heterozygosity mode (dosage heterozygosity
> per band/regime), ROH/FROH mode (inbreeding/IBD confounder check),
> Family/relatedness mode (family confounding check)."*
>
> *"The workflow should be: 1) K-means proposes local PCA bands. 2)
> Dosage, GHSL, heterozygosity, ROH/FROH, and family/relatedness
> diagnostics validate or reject the biological interpretation. 3)
> Candidate-level regimes should be assigned only when multiple
> evidence layers agree."*
>
> *"Do not interpret K=6 as six biological regimes by default. K=6
> may represent nested substructure inside three major regimes. The
> paper-level call should usually use major regimes; higher K values
> are diagnostic/substructure layers."*

This is a **doctrine spec**, not a feature spec. It governs how every
other regime-related feature in the atlas is built and labeled.

---

## 1. The doctrine in one sentence

> *"K-means bands are geometric proposals. Biological regimes require
> agreement across multiple independent evidence layers."*

K-means on PC1 is fast and cheap and usually right enough for a first-pass
karyotype call. But:

- Family LD can produce K=3 bands that look like REF/HET/INV but are
  actually broodline-structure
- ROH-driven local diversity changes can look like K-means structure
- Recombinants can sit in the wrong band
- K=6 sub-bands may be inversion sub-types, family substructure within
  a broader regime, or noise

Until ≥2 independent evidence layers agree, the K-means call is
**provisional**.

## 2. The independent evidence layers

| Layer | What it sees | Independent of |
|---|---|---|
| **K-means on PC1** | geometric clusters in PCA space | (the prior; everything else validates it) |
| **Dosage pattern** | per-marker 0/1/2 genotype per sample | PCA basis |
| **GHSL within-sample** | haplotype divergence per sample | PCA basis, dosage |
| **Heterozygosity** | per-sample het rate in candidate region | PCA, dosage, GHSL |
| **F_ROH local** | local inbreeding intensity per sample | independent |
| **Family/relatedness** | broodline membership per sample | independent (genome-wide) |
| **Q ancestry (NGSadmix)** | local ancestry composition | independent of PCA in candidate |
| **θπ per sample** | nucleotide diversity per sample in region | independent |

## 3. Atlas-side surface — already partial

The atlas has color modes (turn ~50 from chat `47cd29b9`):
- `kmeans` (current default)
- `dosage` (live or static)
- `ghsl`
- `het`
- `family`
- `froh` (planned)
- `confounder_alert` (planned)

Each color mode applies the corresponding evidence layer's per-sample
value as the visual color. Switch modes → visualize different evidence
layer instantly. Per-sample lines panel, L3 contingency panes, main
PCA scatter all support the mode.

## 4. The doctrine codification

### 4.1 Visual labeling

Wherever K-means bands are displayed, label them as **provisional**.
Concretely:
- Toolbar buttons "K=3" become "K=3 (geometric)"
- Tooltips on band labels: *"K-means proposal — confirm with dosage / GHSL
  / het before treating as biological regime."*
- Per-candidate verdict panel (`SPEC_hypothesis_test_framework_atlas.md`)
  shows whether the K-means call is multi-layer-confirmed.

### 4.2 Multi-layer confirmation flag

Per candidate, atlas computes:

```
n_layers_agree = count of layers where the per-sample band assignment
                 from that layer agrees with K-means (above threshold)
```

Layers compared to K-means:
- **Dosage**: cluster samples by mean dosage; agreement = ARI vs K-means
- **GHSL**: cluster by within-sample phase divergence; ARI vs K-means
- **Het**: cluster by het rate; ARI vs K-means
- **F_ROH**: doesn't cluster the same way; flag as "F_ROH-confounder
  detected if F_ROH varies systematically across K-means bands"
- **Family**: flag as "family-confounder detected if K=8 broodlines
  align with K-means bands at ARI > 0.5"

Confidence label:
- **HIGH**: ≥3 of dosage/GHSL/het agree (ARI > 0.7)
- **MODERATE**: 2 of 3 agree
- **LOW**: 0–1 agree
- **CONFOUNDED**: family or F_ROH confounder detected

### 4.3 K=6 special handling

K=6 is **never** a biological regime by default. It's:
- Either nested sub-structure inside K=3 regimes (each K=3 regime splits
  into 2 K=6 sub-bands)
- Or cross-cutting (K=6 bands cut across K=3 regimes)
- Or mixed (partial nesting)

The K=6 substructure verdict (NESTED / MIXED / CROSS_CUTTING / NO_DATA)
already exists per turn 60-ish. This spec formalizes that **only NESTED
K=6 verdict allows interpreting K=6 as biological**. Other verdicts:
flag for review.

## 5. Scope refinement per layer

Each evidence layer has a **primary scope** — the band where it's most
informative:

| Layer | Primary scope | Secondary scope |
|---|---|---|
| Dosage | all bands | — |
| GHSL phase divergence | HET band (highest divergence expected) | HOM bands (low divergence) |
| Heterozygosity | HET band (high het rate expected) | HOM bands (low het rate) |
| F_ROH | all bands | — |
| Family | all bands | — |
| θπ | HET band (elevated diversity) | HOM bands |

Tooltips in the atlas tell the user which scope the current mode is
informative for, so a "no contrast on HOM bands" reading isn't
misinterpreted as "the mode is broken."

## 6. Implementation status

### Already implemented (mark as COMPLETE)
- [x] Color mode infrastructure (turn ~50)
- [x] K-means / dosage / GHSL / het modes (turn ~50–95)
- [x] K=6 substructure verdict + chip (turn ~82)
- [x] Lineage color mode (turn 130 Slice 2)

### Still needed (Slice 1 of this spec, ~0.5 turn)
- [ ] Provisional labeling on K-means buttons + tooltips
- [ ] Multi-layer confirmation flag per candidate
- [ ] HIGH/MODERATE/LOW/CONFOUNDED labels
- [ ] Scope-tooltip system per mode

### Still needed (Slice 2 of this spec, ~0.3 turn)
- [ ] Family-confounder detection (auto-flag when K=8 broodlines
      align with K-means at ARI > 0.5)
- [ ] F_ROH-confounder detection (flag when F_ROH varies systematically
      across K-means bands)

### Still needed (Slice 3 of this spec, ~0.3 turn)
- [ ] `confounder_alert` color mode (red on samples flagged as
      confounded by family or F_ROH)

## 7. Tests

- Synthetic candidate with planted K-means agreement on dosage + GHSL
  → multi-layer confirmation = HIGH.
- Family-confounder synthetic: K-means bands match K=8 broodline at
  ARI > 0.5 → CONFOUNDED label.
- K=6 with nested verdict → atlas allows K=6 as biological. With
  CROSS_CUTTING → flagged as not-biological.

## 8. Cross-references

- Chat `47cd29b9` — original doctrine.
- `SPEC_hypothesis_test_framework_atlas.md` — T1/T2/T9 detect
  family-confounder; this spec surfaces it as a UI flag.
- `SPEC_l3_het_dosage_coloring.md` — Slice 1 shipped turn 129;
  belongs in this multi-evidence framework.
- `SPEC_recombinant_dosage_changepoint_detector.md` — recombinant
  detection complements multi-evidence band confirmation.

## 9. What this is NOT

- Not a refactor of existing K-means infrastructure. K-means stays as
  the prior; this spec just relabels it as provisional.
- Not a new clustering algorithm. Other layers cluster independently
  using their own logic; this spec defines how to compare those
  clusterings to K-means.
- Not a replacement for the hypothesis test framework
  (`SPEC_hypothesis_test_framework_atlas.md`). The hypothesis tests
  are the rigorous version; this spec is the visual evidence layer
  that helps the user spot problems before running tests.
