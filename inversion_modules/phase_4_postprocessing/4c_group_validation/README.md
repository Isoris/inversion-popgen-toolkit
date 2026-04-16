# `4c_group_validation/` — validate proposed groups

Phase 4c runs hypothesis tests (T1–T11) on the groups proposed by phase
4b and promotes or demotes the group validation level accordingly. Uses
the v9.3.4 hypothesis-test suite with three integrated patches.

## Contents

```
4c_group_validation/
├── STEP_C01f_hypothesis_tests.R   ← integrated (v9.3.4 + patches 01/02/03)
├── _PRISTINE_v9.3.4.R             ← pristine v9.3.4 for diffing
└── README.md                       ← this file
```

## What's integrated

`STEP_C01f_hypothesis_tests.R` = the v9.3.4 source from
`STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R` with
three patches applied in order. See `../patches/` for the originals.

### Patch [01] — `comp_from_registry` + `build_comp_with_fallback`

**Location:** inserted just before `# T3: Kin-pruned similarity matrix comparison`
(the first test function definition, around line 280).

**What it does:** adds two helper functions that let C01f read group
composition from the C01i sample_registry rather than from
`triangle_sample_composition.tsv.gz`. Falls back to the file if the
registry isn't populated.

**Band mapping (registry → file legacy):**
```
HOM_REF     → band1  (lowest PC1)
HET         → band2
HOM_INV     → band3  (highest PC1)
RECOMBINANT → EXCLUDED from comp (violates "one band per sample" assumption
               of T1/T2/T7/T8/T10)
```

### Patch [02] — promotion_cap argument

**Merged into** the `compute_group_validation()` function inserted by
patch [01]. Reads `q6_validation_promotion_cap` from the evidence
registry (set by `4b_group_proposal/STEP_C01i_d_seal.R` when
`composite_flag = likely_composite`). If set to "UNCERTAIN", the
function refuses to promote above UNCERTAIN regardless of test
evidence strength.

### Patch [03] — four-way jackknife semantics

**Merged into** `compute_group_validation()`. Return type is now a list
with fields `$level` / `$quality_flags` / `$family_linkage`. The
jackknife verdict from T9 (cheat6) is classified four ways:

| verdict (T9)              | family_linkage        | validation          |
|---------------------------|-----------------------|---------------------|
| robust_multi_family       | multi_family          | can VALIDATE        |
| multi_family_contributing | multi_family          | SUPPORTED (w/ tag)  |
| few_family                | few_family            | cap at SUPPORTED    |
| single_family_fragile     | single_family         | cap at SUPPORTED    |
| fragile / pca_family_driven | pca_family_confounded | SUSPECT           |

**Key insight:** `single_family_fragile` is NOT SUSPECT. It's a real,
family-restricted polymorphism. Only `pca_family_confounded` becomes
SUSPECT.

## Where each patch was applied in the file

| Patch | Section | Purpose |
|---|---|---|
| (A) Helper + compute_group_validation | around line 282 (before T3) | define new functions |
| (B) Per-candidate comp build | ~line 1951 (start of per-candidate loop) | prefer registry, fall back to file |
| (C) Validation-level compute + write | ~line 2200 (verdict row construction) | read current + cap, call compute_group_validation, write 5 new columns |

## New columns written to `hypothesis_verdicts.tsv`

The verdict row (and therefore the output TSV + evidence registry) now
has 5 additional columns:

| Column | Values | Set by |
|---|---|---|
| `group_validation_before` | UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT | C01i seal initial state |
| `group_validation_after`  | UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT | compute_group_validation result |
| `quality_flags`           | comma-separated (`normal`, `partial_robustness`, `restricted_family_spread`, `family_specific_polymorphism`, `pca_groups_family_confounded`, `composite_cap` from seal, `high_jackknife_delta_no_verdict`) | |
| `family_linkage`          | multi_family / few_family / single_family / pca_family_confounded / uncertain / unknown | T9 verdict classification |
| `promotion_cap_applied`   | name of cap if triggered, empty otherwise | audit trail |
| `comp_source`             | `C01i_registry` or `file_fallback` | how `comp` was built |

## The `reg$evidence$write_block` path

When the registry is available, the verdict row should be written as a
Tier-2 `hypothesis_verdict` block. The schema's `keys_extracted` rule
automatically promotes `group_validation_after` → `q6_group_validation`,
`family_linkage` → `q6_family_linkage`, `quality_flags` → `q6_quality_flags`.

When the registry is NOT available, flat writes via `register_C01f_keys()`
(the fallback path at the bottom of the script) still work because the
new columns are just additional fields on the verdict row.

## Running

Same as v9.3.4 — no CLI args changed:

```bash
Rscript STEP_C01f_hypothesis_tests.R \
    --scores <candidate_scores.tsv.gz> \
    --triangles <triangle_dir> \
    --precomp <precomp_dir> \
    --relatedness <pairs.tsv> \
    --samples <sample_list> \
    [--pruned_samples <kin_pruned_list>] \
    [--theta_dir <theta_dir>] \
    [--qmatrix <q_matrix>]
```

## Testing

The pure-logic portion of `compute_group_validation()` is tested by:

- `../tests/test_compute_group_validation.py` — 22 unit tests (v10 base)
- `../tests/test_jackknife_semantics.py` — 15 tests (v10.1.1 jackknife)
- `../tests/test_phase4b_integration.py` — 18 end-to-end tests

These are Python mirrors of the R function, testing the same decision
table. If the R function drifts from the Python, the tests will flag
it. Run them after any future edit to `compute_group_validation()`:

```bash
cd ../tests/
for t in test_compute_group_validation.py test_jackknife_semantics.py \
         test_phase4b_integration.py; do
    python3 "$t" || break
done
```

All 5 test suites passed as of the last deployment check:
`python3 ../../check_deployment_v10.py` → `✓ ALL GREEN` (77 tests).

## Diffing the pristine source

To see exactly what the 3 patches changed:

```bash
diff -u _PRISTINE_v9.3.4.R STEP_C01f_hypothesis_tests.R
```

The only changes should be:
- Header rewrite (purely documentation)
- Insertion of `comp_from_registry` / `build_comp_with_fallback` /
  `compute_group_validation` before the T3 function
- Replacement of the per-candidate loop opening (lines ~1944-1976 in
  pristine) with the registry-aware version
- Replacement of the verdict row construction (lines ~2200-2249 in
  pristine) with the registry-aware version including new columns

## Pending

- **Schema update:** the `hypothesis_verdict.schema.json` needs to bump
  to v2 to include the new `family_linkage` enum
  (multi_family / few_family / single_family / pca_family_confounded /
  uncertain / unknown) and the new `quality_flags` field. See
  `../patches/03_C01f_jackknife_semantics.R` line 162-167.
- **Real R parse on LANTA:** bracket-balance does not equal parse-valid.
  Run `Rscript --vanilla -e "invisible(parse(file='STEP_C01f_hypothesis_tests.R'))"`
  on LANTA before first submit.
- **One-candidate smoke test** to verify the comp_from_registry path
  actually reaches T1-T11 with the same results as the file-fallback
  path on a candidate that has both (registry groups + file comp).

## See also

- `../README.md` — phase 4 overview
- `../patches/README.md` — the 3 patches, in application order
- `../patches/_archive_superseded/README.md` — 3 patches NOT to apply
  (two C01i + the obsolete flashlight C01f patch)
- `../4b_group_proposal/` — where the groups come from (4b.4 seal
  writes q6_validation_promotion_cap that this phase respects)
- `../docs/PHASE4_ARCHITECTURE.md` — the overall phase 4 design
