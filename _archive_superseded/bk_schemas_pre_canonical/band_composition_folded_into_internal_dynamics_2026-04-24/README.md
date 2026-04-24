# band_composition — superseded

**Date folded:** 2026-04-24 (chat C continuation)
**Folded into:** `internal_dynamics.schema.json`

## Why

`band_composition` was declared as a standalone 2-key schema back when the
Phase 4b rewrite was first scoped. Its claimed writer was
`STEP_C01i_decompose.R + unified_ancestry/instant_q`, and the two keys it
extracted were:

  - `q1_k_used` — K used for decompose clustering (2 or 3).
  - `q1_ancestry_div_hom_ref_vs_hom_inv` — L1 distance between mean Q
    vectors of the HOM_REF vs HOM_INV classes, flagging family-driven
    signals where the two homozygous classes come from different
    ancestral subpopulations.

During the 2026-04-24 audit (chat B / chat C), two things became clear:

1. **`q1_k_used` is already written** by `STEP_C01i_decompose.R` under the
   existing `internal_dynamics` block as `q2_pca_k_used`. Having the
   same data expose two different flat-key names via two different
   schemas is the exact kind of duplication the registry tries to avoid.
   (The `q2_` prefix denotes the "how did decompose cluster samples?"
   interpretation; `q1_` the "how many bands does the block have?"
   interpretation. Same K either way.)

2. **`q1_ancestry_div_hom_ref_vs_hom_inv` is NOT computed anywhere.** It
   requires joining decompose's per-sample class assignment
   (HOM_REF/HET/HOM_INV) with Engine B's per-window local-Q vectors
   inside the candidate interval, averaging Q per class, then taking the
   L1 distance between the HOM_REF and HOM_INV mean vectors. Decompose
   does not currently pull instant_q data and has no ancestry surface.

A standalone schema for two keys that share the decompose writer was
architectural over-engineering. Quentin's preference (2026-04-24) was
option (a) of the audit: promote the keys into `internal_dynamics` and
retire the standalone schema.

## What changed

- Added two `keys_extracted` entries to `internal_dynamics.schema.json`:
  - `q1_decompose_k_used` ← `k_used` (duplicates existing
    `q2_pca_k_used`; dual-key shim for back-compat with any reader of
    the old band_composition vocabulary). The explicit `decompose_`
    prefix in the key name makes the provenance visible in `keys.tsv`
    so a reader scanning the file can tell which script produced the
    fact.
  - `q1_decompose_ancestry_div_hom_ref_vs_hom_inv` ←
    `ancestry_div_hom_ref_vs_hom_inv` (new field in the
    internal_dynamics block; decompose writes NA until the instant_q
    join is added). Same `decompose_` prefix convention — future
    ancestry-divergence values computed by a different script (e.g.
    `q1_unified_ancestry_div_*`) would use that script's prefix.
- Added the `ancestry_div_hom_ref_vs_hom_inv` property to
  `internal_dynamics.schema.json#properties` with description + nullable
  number type.
- Extended `STEP_C01i_decompose.R`'s `block_data` list to emit
  `ancestry_div_hom_ref_vs_hom_inv = NA_real_` in the `seeded_ok` block
  (the two `no_seeding` stub blocks do not include it — consistent with
  their minimal payload).
- Contract header on `STEP_C01i_decompose.R` updated to declare
  `internal_dynamics` WIRED, with `q1_decompose_ancestry_div_hom_ref_vs_hom_inv`
  listed under `keys_na`.
- This file (`band_composition.schema.json`) archived here as a
  historical record; no live code reads it.

## Remaining work: compute the Q divergence

The only outstanding piece is the instant_q join inside decompose, so
`ancestry_div_hom_ref_vs_hom_inv` stops being NA. Rough sketch:

1. Decompose already has `class_assignment` (per-sample HOM_REF / HET /
   HOM_INV labels) by L437.
2. Load the per-window local-Q TSV for the candidate's chromosome from
   the Engine B cache (path convention: `${LOCAL_Q_DIR}/K*/${chr}.local_Q_samples.tsv.gz`).
3. Filter to windows overlapping the candidate interval.
4. Group by sample, average Q across windows, then group samples by
   `class_assignment` and average per class → 3 mean Q vectors.
5. Compute L1 distance between the HOM_REF and HOM_INV mean vectors.
   Populate `block_data$ancestry_div_hom_ref_vs_hom_inv` with the
   scalar.

Estimate: ~30 lines in decompose, reusing `load_q_samples` from
`unified_ancestry/engines/nested_composition/internal_ancestry_composition.py`
as a reference for the file format. An R port of that loader would sit
alongside `lib_decompose_helpers.R`. Deferred until someone needs the key
to land as non-NA; the schema contract is ready when they do.

## No action required from consumers today

The old band_composition schema exposed `q1_k_used` and
`q1_ancestry_div_hom_ref_vs_hom_inv`. The fold renames these to
`q1_decompose_k_used` and `q1_decompose_ancestry_div_hom_ref_vs_hom_inv`
respectively, making the `STEP_C01i_decompose.R` provenance explicit in
the key name (consistent with existing `q2_pca_*` and `q6_*_prelim`
patterns in this block). Both keys appear in `keys.tsv` via the
`internal_dynamics` block. Same semantics for `q1_decompose_k_used`;
NA for `q1_decompose_ancestry_div_hom_ref_vs_hom_inv` until the Q join
lands.

If any consumer was coded against the literal strings `q1_k_used` or
`q1_ancestry_div_hom_ref_vs_hom_inv`, it needs a one-line update to the
new names. Grep across the codebase at fold time found no such consumers
(band_composition had never been wired, so nothing read from it).
