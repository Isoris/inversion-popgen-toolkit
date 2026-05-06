# SPEC — Per-band evidence layer (H, θπ, F_ROH, family per K-band)

**Status**: drafted turn 130 final session. Roadmap items from chat
`819d8454` and the scrubber wishlist. Some shipped (per-band het via
turn 93's het-shape module); others pending JSON exporter changes.

**Trigger** (Quentin via wishlist):
> *"Per-sample H per band — inversion vs not inversion."*
> *"π / θ curve per band — nucleotide diversity profile."*
> *"F_ROH per band — inbreeding-driven confounder check."*
> *"Family by band — relatedness confounder check."*

Each K-band of a candidate should be queryable by additional
diagnostics beyond just dosage. Some are shipped, some need new JSON
producers, some are atlas-derivable.

---

## 1. Layer inventory

| Diagnostic | Per-sample? | Status | Producer |
|---|---|---|---|
| **Het-rate per band** | yes (per K-band) | shipped turn 93 (het-shape) | atlas-derived from dosage_chunks |
| **K-means group color** | one value per sample | shipped | atlas-derived |
| **Per-sample θπ in candidate** | yes | spec'd, needs producer | LANTA `STEP_R36_emit_per_sample_theta_pi.R` (planned) |
| **Per-band θπ summary** | summary per band | not yet | LANTA Engine B per-sample θπ → atlas group |
| **Per-sample F_ROH local** | yes | not yet | LANTA Module 3 ngsF-HMM intersected with candidate region |
| **Per-band F_ROH summary** | summary per band | not yet | atlas-derived from per-sample F_ROH |
| **Per-band family composition** | summary per band | shipped (sidebar text) | atlas-derived from family_id metadata |
| **Per-band damaging-load** | summary per band | not yet | LANTA `MODULE_CONSERVATION` → atlas group |

## 2. Why per-band matters

Each K-band represents a putative karyotype class. Properly characterized
inversions show **predictable per-band profiles**:

| Band | Expected H | Expected θπ | Expected F_ROH | Expected family |
|---|---|---|---|---|
| **HOM_REF** | low (homozygous) | low (one arrangement) | variable | mixed |
| **HET** | high (heterozygous) | elevated (two arrangements) | variable | mixed |
| **HOM_INV** | low (homozygous) | low (one arrangement) | variable | mixed |

If, instead, per-band F_ROH is asymmetric (one band has uniformly high
F_ROH) → **possible F_ROH-confounder**. The K-means is detecting
inbreeding intensity, not inversion karyotype.

If per-band family composition is asymmetric (one band is 90%
broodline-3) → **family-confounder**. The K-means is detecting
broodline structure, not inversion karyotype.

If per-band θπ doesn't show the HET-elevation signal → **suspicious**
candidate. Either it's not a real inversion or the karyotype assignment
is wrong.

## 3. Atlas surface

Per-candidate page-2 view extends the existing band-diagnostics table
(turn 93's het-shape) with new columns:

```
Band | n  | het_med | het_IQR | θπ_med | F_ROH_med | family_top | family_pct
g0   | 60 | 0.05    | 0.02    | 0.0034 | 0.18      | broodline2 | 35%
g1   | 106| 0.42    | 0.08    | 0.0089 | 0.15      | mixed      | n/a
g2   | 60 | 0.04    | 0.02    | 0.0028 | 0.21      | broodline2 | 38%
```

Color cells red where:
- F_ROH varies > 2× across bands (F_ROH-confounder flag)
- One family dominates one band > 50% (family-confounder flag)
- HET band's θπ is NOT > HOM bands' θπ (θπ inversion-confirmation
  fail)

Click any cell → opens detail popover.

## 4. Confounder flags

The atlas computes per-candidate confounder flags:

```
F_ROH_confounder: max(band_F_ROH_med) / min(band_F_ROH_med) > 2.0
family_confounder: any band with > 50% same family_id
theta_pi_confirmation_fail: θπ(HET) <= θπ(HOM_REF) AND θπ(HET) <= θπ(HOM_INV)
```

These flags appear in:
- The per-candidate verdict panel (`SPEC_hypothesis_test_framework_atlas.md`)
- The catalogue page-3 view as a column
- The breeding-readiness card
- The manuscript bundle's per-candidate paragraph

## 5. JSON layer schema

```
sample_theta_pi_in_candidate_v1: {
  per_candidate: {
    [cand_id]: {
      [sample_id]: theta_pi_value
    }
  }
}

sample_froh_in_candidate_v1: {
  per_candidate: {
    [cand_id]: {
      [sample_id]: froh_value
    }
  }
}
```

Each file is per-chromosome bundle, layered into the candidate JSON.
Atlas detects via `detectSchemaAndLayers`.

## 6. Implementation slices

### Slice 1 — atlas-side per-band computation (~0.5 turn)
- `_computePerBandH(candidate)` — already exists via het-shape; expose
  as per-band column
- `_computePerBandFamily(candidate)` — straightforward; family_id
  is per-sample metadata
- Per-band table extension on page-2

### Slice 2 — confounder flag computation (~0.3 turn)
- F_ROH-confounder, family-confounder, θπ-fail flags
- Cell-coloring on the per-band table
- Flag aggregation into candidate-level summary

### Slice 3 — LANTA-side producer for per-sample θπ in candidate (~outside atlas)
- `STEP_R36_emit_per_sample_theta_pi.R` — runs Engine B
  region_popstats.c per candidate per sample
- Emits `sample_theta_pi_in_candidate_v1` JSON layer

### Slice 4 — LANTA-side producer for per-sample F_ROH in candidate (~outside atlas)
- Intersect candidate region with each sample's ngsF-HMM ROH calls
- Emit `sample_froh_in_candidate_v1` JSON layer

### Slice 5 — atlas-side rendering once layers ship (~0.3 turn)
- Detect new JSON layers
- Add θπ and F_ROH columns to per-band table
- Confounder flags activate when layers present

## 7. Open questions

1. **Are per-band summaries enough, or do we need per-sample plots?**
   Per-band summaries fit the verdict-panel surface. Per-sample dot
   plots fit the deep-dive surface. Probably both, with the dot plot
   in a click-to-expand drawer.
2. **What's the confounder-flag threshold?** F_ROH ratio >2 is a
   guess. Calibrate on real LG28 once data is available.
3. **Does damaging-load belong here?** Yes for the breeding-readiness
   card. No for the everyday band-diagnostics view. Surface in the
   card only.

## 8. Cross-references

- `SPEC_multi_evidence_regime_framework.md` — per-band evidence
  layers feed the regime-confidence assessment.
- `SPEC_hypothesis_test_framework_atlas.md` — confounder flags
  inform T1/T2/T9 verdicts.
- `SPEC_per_candidate_breeding_readiness_card.md` — consumes per-band
  burden columns.
- Chat `819d8454` for the original wishlist.
- Chat `5b793a68` for the metric inventory.
