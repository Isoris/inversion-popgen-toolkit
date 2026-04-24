# Frequency schema — reconciliation archive

**Date:** 2026-04-24 (chat B continuation)

## What happened

`frequency.schema.json` (v1) and `frequency.v2.schema.json` both claimed
`block_type = "frequency"` — a genuine collision. Loading the schema by
block_type was ambiguous. v1 had better documentation; v2 had better
semantics. Neither was complete.

Reconciled into `frequency.v3.schema.json` at
`registries/schemas/structured_block_schemas/frequency.v3.schema.json`.

## What moved

- `frequency.v1.schema.json` (here) = the original
  `frequency.schema.json`. 4-class `freq_class` enum; 2-value
  `family_linkage`; no `jackknife_max_delta` / `polymorphism_class`;
  but richer HWE documentation and `expected_counts_hwe`.
- `frequency.v2.schema.json` (here, unchanged name) = the
  `frequency.v2.schema.json` from the tree. 5-value `family_linkage`
  (missing `uncertain`); added `jackknife_max_delta`,
  `jackknife_n_contributing` (not extracted), `polymorphism_class`
  (not derived). Lost v1's HWE docs.

## What v3 did differently

| Change | Source |
|---|---|
| Merged `expected_counts_hwe`, `note_on_hatchery_hwe`, richer HWE descriptions | from v1 |
| Kept 5-class `family_linkage` but added missing `uncertain` value | v2 + bug fix |
| Kept `jackknife_max_delta`, `polymorphism_class` | from v2 |
| Added `jackknife_n_contributing` to `keys_extracted` (v2 had it as property but never extracted) | bug fix |
| Widened `freq_class` enum from 4-class (rare/low/common/major at 0.05/0.20/0.50) to 5-class (rare/low/intermediate/high/nearly_fixed at 0.05/0.15/0.50/0.85) to match `characterize_candidate.R`'s actual classification | consumer/schema reconciliation |
| Flattened `class_counts` into `q6_n_HOM_REF`, `q6_n_HET`, etc. as flat keys | gap flagged in SPEC_VS_REALITY.md |
| Added `computed_from` provenance field | new |
| Documented that `polymorphism_class` is derived helper-side in `utils/registry_key_helpers.R::update_C01f_frequency_linkage`, not in C01f source | new |

## Why keep them

Audit trail. If someone later asks "what did `family_linkage` look like
before v3?" you can answer without git archaeology.
