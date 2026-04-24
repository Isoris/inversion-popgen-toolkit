# `patches/` — apply in numeric order to your current C01f

Three patches that modify `STEP_C01f_hypothesis_tests.R` in your current
inversion codebase. Apply in this exact order. Each one depends on the
previous.

## Application order

```
01_C01f_comp_from_registry.R        ← v10 base patch
02_C01f_promotion_cap.R             ← v10.1 phase 4b rewrite
03_C01f_jackknife_semantics.R       ← v10.1.1 addition
```

## What each patch does

### 01 — `C01f_comp_from_registry.R`

Three insertions marked (A), (B), (C):

- **(A)** Adds two new helper functions (`comp_from_registry()` and
  `build_comp_with_fallback()`) near the top of the script, before the
  first `test_*` function. Insert marker: line containing
  `# T3: Kin-pruned similarity matrix comparison`, insert BEFORE.

- **(B)** Patches the main-loop call site that builds `comp` via PC1
  k-means. **Only the main-loop call site.** Leave T4/T5/T6 internal
  k-means calls alone — those are per-subregion by design.

- **(C)** Adds `compute_group_validation()` (pure function) at the end
  of the per-candidate loop. Its return value goes into the verdict
  row's `group_validation_after` field. The schema extracts that into
  `q6_group_validation` automatically when the block is written. Do NOT
  also call `reg$add_evidence` separately — that creates a race with
  the block's `keys_extracted`.

### 02 — `C01f_promotion_cap.R`

Replaces `compute_group_validation()` with a cap-aware version that
reads `q6_validation_promotion_cap` from the evidence registry. If the
cap is set to UNCERTAIN (which 4b.4 `STEP_C01i_d_seal.R` sets when
`composite_flag = likely_composite`), the function refuses to promote
above UNCERTAIN regardless of test evidence strength.

This is the piece that makes composite intervals safe: C01f cannot
promote them to SUPPORTED/VALIDATED, so downstream Fst/burden/age
cheats in phase 4d refuse to run on their groups.

### 03 — `C01f_jackknife_semantics.R`

Another replacement of `compute_group_validation()`, this time
implementing the **four-way jackknife classification** (v10.1.1):

| Cheat6 verdict            | family_linkage        | validation    |
|---------------------------|-----------------------|---------------|
| robust_multi_family       | multi_family          | can VALIDATE  |
| multi_family_contributing | multi_family          | SUPPORTED     |
| few_family                | few_family            | cap SUPPORTED |
| single_family_fragile     | single_family         | cap SUPPORTED |
| pca_family_confounded     | pca_family_confounded | SUSPECT       |

**Key change**: single-family-fragile is NOT SUSPECT. It's a
family-restricted polymorphism — real, legitimate, lineage-specific.
The ONLY verdict that still becomes SUSPECT is
`pca_family_confounded` (where PCA groups reflect ancestry, not
arrangement).

**Return type changed**: the function now returns a list with fields
`$level`, `$quality_flags`, `$family_linkage` instead of a single
string. **Your call sites need a 3-line update** to unpack the list.
See patch 03 comments for the exact call-site change.

## Related doc

`ENGINE_B_SMOKE_TEST_INSERTS.md` — how to wire the 5-second Engine B
self-test (`../../phase_7_karyotype_groups/proposal/engine_b_smoke_test.R`) into C01a,
cheat6, and the Python pipeline. Not a patch in the mechanical sense;
instructions for where to add `source()` + `run_engine_b_smoke_test()`
calls.

## Prerequisite

Your current `STEP_C01f_hypothesis_tests.R`. Upload it to
`phase_7_karyotype_groups/validation/` when you're ready to apply these patches. The
patches will need your current file's exact layout for the (B) section
in patch 01 especially.

## Superseded patches

See `_archive_superseded/README.md` — two C01i patches from the v10
base that are made obsolete by the 4b rewrite.

## Testing

After applying all 3 patches, run:

```bash
cd inversion_modules/phase_9_classification/tests/
python3 test_compute_group_validation.py     # 22 unit tests
python3 test_jackknife_semantics.py          # 15 tests
python3 test_phase4b_integration.py          # 18 tests
```

All should pass. The logic in each test is the pure Python mirror of
the R function; if the 3 patches are applied correctly the R function
returns identical values to the Python mirror.
