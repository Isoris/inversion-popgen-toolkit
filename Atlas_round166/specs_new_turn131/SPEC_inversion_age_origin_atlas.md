# SPEC — Inversion age & origin classification (atlas surface)

**Status**: drafted turn 130 final session. The classifier exists
LANTA-side as `cheat30_ibs_by_genotype.R` and `cheat15_recurrence_test.R`
(per chats `c1a88579`, `2beb520f`). This spec brings the **atlas-side
visualization + integration** into the spec folder.

**Trigger** (Quentin, chat `c1a88579`):
> *"GDS unimodality within INV carriers (Hartigan's dip test on INV/INV
> pairwise GDS) — fail to reject unimodality → single origin.
> Bimodality → recurrent origin (Porubsky 2022 Fig 3C/D)."*
>
> *"GDS gap (arrangement divergence) >0.08 = old, 0.03-0.08 =
> intermediate, <0.03 = young."*

For each confirmed inversion, the atlas should display a per-candidate
**age + origin classification** based on within-arrangement GDS
(Genotype Dosage Similarity).

---

## 1. The science (quick recap)

For each pair of samples (i, j):
```
GDS(i, j) = 1 - mean(|d_i - d_j|) / 2
```
where `d_i, d_j` are dosage vectors at sites within the inversion.

Stratify GDS values by genotype-pair:
- **REF-REF pairs**: high GDS (same arrangement, similar haplotypes)
- **REF-INV pairs**: lower GDS (different arrangements)
- **INV-INV pairs**: how high?

The **INV-INV distribution shape** reveals origin:
- Unimodal, tight, high mean → single origin (one founding mutation
  spread by drift)
- Bimodal or wide → recurrent origin (multiple independent flips)

The **REF-INV gap** (mean REF-REF GDS − mean REF-INV GDS) reveals age:
- Large gap → arrangements diverged for many generations → **OLD**
- Small gap → arrangements diverged recently → **YOUNG**

## 2. Atlas-side surface

### 2.1 Per-candidate age + origin badge

On the page-2 candidate-focus page:

```
Inversion age & origin
─────────────────────────────────────────────
GDS gap:          0.124  → OLD
INV-INV dip test: p = 0.34 → unimodal
Origin class:     CONFIRMED SINGLE ORIGIN
              (3/3 evidence: unimodal + NHEJ + Fst-age agrees)
```

### 2.2 GDS distribution plot

Density plot showing three distributions overlaid:
- REF-REF (blue)
- REF-INV (red)
- INV-INV (green)

If INV-INV is bimodal, two peaks visible — visual confirmation of
recurrence.

### 2.3 Cross-candidate consistency

A scatter plot: x = GDS gap, y = arrangement Fst (Spearman > 0.5
expected for a coherent age signal). Outliers in this plot are
candidates whose age estimate doesn't agree across methods — flag
for review.

## 3. Classification logic (from Cheat 30 / Cheat 15)

### Origin class (Claim 5/6 in `c1a88579`)
| n_clusters in INV | mechanism | classification |
|---|---|---|
| 1 (unimodal) | NHEJ | confirmed_single |
| 1 | UNKNOWN | likely_single |
| 1 | NAHR | substrate_exists_but_single |
| ≥2 (bimodal) | NAHR | confirmed_recurrent |
| ≥2 | UNKNOWN | likely_recurrent |
| ≥2 | NHEJ | unexpected_multi_cluster (review) |

### Age class
| GDS gap | classification |
|---|---|
| >0.08 | OLD |
| 0.03–0.08 | INTERMEDIATE |
| <0.03 | YOUNG |

## 4. JSON layer schema

`origin_age_classification_v1`:

```json
{
  "per_candidate": {
    "LG28_inv5": {
      "gds_ref_ref_mean": 0.821,
      "gds_ref_inv_mean": 0.697,
      "gds_inv_inv_mean": 0.834,
      "gds_gap": 0.124,
      "dip_test_p": 0.34,
      "n_inv_clusters": 1,
      "mechanism_class": "NHEJ",
      "age_class": "OLD",
      "origin_class": "confirmed_single",
      "fst_arrangement": 0.42,
      "evidence_count": 3,
      "evidence_supporting": ["unimodal_dip", "nhej_mechanism", "fst_age_consistent"]
    }
  },
  "cross_candidate_consistency": {
    "spearman_rho_gap_vs_fst": 0.61,
    "p_value": 1.4e-5
  }
}
```

## 5. Implementation slices

### Slice 1 — JSON layer detection + badge (~0.3 turn)
- Detect `origin_age_classification_v1` layer
- Render age + origin badge on page-2

### Slice 2 — GDS distribution plot (~0.5 turn)
- Three-distribution density plot per candidate
- Tooltips with means and p-values

### Slice 3 — cross-candidate consistency view (~0.3 turn)
- Scatter plot + flag outliers for review

## 6. Open questions

1. **Where does the Fst arrangement come from?** Engine B
   region_popstats.c with HOM_REF and HOM_INV as the two groups. Already
   computable per `region_stats_dispatcher.R`.
2. **What about candidates with too few HOM_INV?** GDS test needs
   ≥5 INV-INV pairs. If fewer, mark age + origin as UNKNOWN. Don't
   guess.
3. **NHEJ vs NAHR mechanism**: requires breakpoint sequence analysis
   (microhomology, inverted SDs). LANTA-side via Cheat 29 / Cheat 27.
   Consumed but not computed by atlas.

## 7. Tests

- Synthetic candidate with planted unimodal INV-INV distribution →
  classified as single_origin.
- Synthetic with two distinct INV haplotype clusters → classified as
  recurrent.
- Synthetic with very small REF-INV gap → YOUNG.

## 8. Dependencies

- LANTA-side `cheat30_ibs_by_genotype.R` (existing).
- LANTA-side `cheat15_recurrence_test.R` (existing).
- Engine B Fst between arrangements (existing).
- Per-candidate karyotype calls (existing).

## 9. Cross-references

- Chat `c1a88579` — full classifier design.
- Chat `2beb520f` — Cheat 15 implementation.
- `SPEC_per_candidate_breeding_readiness_card.md` — origin/age fields
  in the breeding card.
- `SPEC_manuscript_bundle_export.md` — boilerplate prose includes
  age/origin per candidate.
- Existing `SPEC_age_origin_panel.md` (already in folder) covers some
  of the same ground; this spec consolidates.

## 10. What this is NOT

- Not a phylogenetic dating method. The age estimate is relative
  ("OLD vs YOUNG"), not absolute Mya.
- Not a phasing-resolved mechanism caller. Mechanism (NHEJ/NAHR/MMEJ)
  comes from upstream sequence analysis.
- Not a substitute for cross-species comparison. For absolute age,
  catfish-synteny-toolkit is the right tool.
