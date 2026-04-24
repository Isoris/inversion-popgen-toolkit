# `phase_9_classification/` — final catalog with 14-axis classification

This sub-phase is **not yet populated**. Upload `classify_inversions.R`
and `characterize_candidate.R` here when ready.

## What runs here

| Script | Purpose |
|---|---|
| `classify_inversions.R` | reads ALL evidence (layers A–D + decomposition + hypothesis verdict + all cheats), assigns 14-axis classification |
| `characterize_candidate.R` | per-candidate summary narrative |

## The 14 classification axes

From the phase 4 architecture doc:

1. `existence_layer_a` — pass / fail (from 4a)
2. `existence_layer_b` — pass / fail (from 4a)
3. `existence_layer_c` — pass / fail (from 4a)
4. `existence_layer_d` — pass / fail (from 4c)
5. `boundary_quality` — sharp / fuzzy (from 4a + phase 3)
6. `group_validation` — NONE / UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT (from 4c)
7. `internal_structure` — clean / gradient / composite_undecomposed (from 4b.3)
8. `recombinant_class` — none / gene_conversion / double_crossover / mixed (from 4b.2)
9. `family_linkage` — multi_family / few_family / single_family / pca_family_confounded / unknown (from 4c)
10. `polymorphism_class` — cohort_wide / lineage_restricted / family_restricted / unclassified (from 4c)
11. `mechanism_class` — NAHR / NHEJ / MMBIR / unknown (from 4d cheat20)
12. `age_class` — young / intermediate / ancient / unknown (from 4d Q5)
13. `burden_class` — enriched / neutral / depleted / unknown (from 4d Q6)
14. `confidence_tier` — T1 (all layers pass + VALIDATED) / T2 / T3 / T4 (SV-only or Layer-C-only)

## Current state

Per the v10 phase 4 README:

> characterize_candidate and classify_inversions not migrated. The
> architecture doc describes where they fit in phase 4e, but they still
> consume the old ad-hoc TSVs. Migration is mechanical once the
> library is in place.

## Output

- `final_catalog.tsv` — one row per candidate with all 14 axes + coords
- Per-candidate registry entry at
  `evidence_registry/per_candidate/<cid>/final_classification.json`

## Upload checklist

- Current source of `classify_inversions.R`
- Current source of `characterize_candidate.R`
- An example `final_catalog.tsv` from a past run, if available
- Whatever lookup tables / rule files they consume
