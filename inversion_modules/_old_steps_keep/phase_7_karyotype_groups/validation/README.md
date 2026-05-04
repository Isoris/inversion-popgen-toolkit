# `phase_7_karyotype_groups/validation/` — validate proposed groups

C01f runs hypothesis tests T1–T11 on the groups proposed by C01i
(in `../proposal/`) and promotes or demotes the group validation level
accordingly. Uses the v9.3.4 hypothesis-test suite with three integrated
patches.

## Contents

```
phase_7_karyotype_groups/validation/
├── STEP_C01f_hypothesis_tests.R   ← integrated (v9.3.4 + patches 01/02/03)
├── LAUNCH_C01f_hypothesis.sh      ← SLURM launcher
└── README.md                      ← this file
```

### File-location history (pass 16, 2026-04-24)

- `_PRISTINE_v9.3.4.R` (pristine pre-patch source of
  `STEP_C01f_hypothesis_tests.R`) was archived to
  `_archive_superseded/phase_7_validation_pristine_v934/`. See that
  directory's README for diff instructions.
- `group_validation_gate.R` was moved to
  `../../phase_9_classification/group_validation_gate.R` because both
  its callers live there and the `group_validation_minimum` per-Q table
  is a characterization rule (phase_9), not a group-proposal concept
  (phase_7). The deprecated stub at the old path was removed in pass 20
  (2026-04-24) once both callers were confirmed migrated. The
  writer-side counterpart `compute_group_validation()` stays here
  inside `STEP_C01f_hypothesis_tests.R` — see §"Writer / reader split"
  below.

## What's integrated

`STEP_C01f_hypothesis_tests.R` = the v9.3.4 source from
`STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R` with three
patches applied in order. See `../../phase_9_classification/patches/`
for the originals.

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
registry (set by `../proposal/STEP_C01i_d_seal.R` when
`composite_flag = likely_composite`). If set to "UNCERTAIN", the
function refuses to promote above UNCERTAIN regardless of test evidence
strength.

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

## Writer / reader split — where is the validation logic?

The `compute_group_validation()` function inside
`STEP_C01f_hypothesis_tests.R` is the **writer** — given the T1–T11
test outputs, it computes the validation level and writes
`q6_group_validation` / `q6_family_linkage` / `q6_quality_flags` to the
registry.

The companion **reader**, `assess_group_validation()` in
`../../phase_9_classification/group_validation_gate.R`, re-derives a
validation level from those flat keys at characterization time so that
phase_9 can gate per-Q convergence without re-running C01f. The two
functions implement the same decision table but from different sides:

| Function | Lives in | Input | Output |
|---|---|---|---|
| `compute_group_validation()` | phase_7/validation (this dir) | T1–T11 raw outputs | writes `q6_*` flat keys |
| `assess_group_validation()` | phase_9 | `q6_*` flat keys (read) | returns `level` + `reason` |
| `check_group_gate()` | phase_9 | question + level | passes / fails |

The Python test mirrors live in `../../phase_9_classification/tests/`:
- `test_compute_group_validation.py` — 22 unit tests (v10 base)
- `test_jackknife_semantics.py` — 15 tests (v10.1.1 jackknife)
- `test_phase4b_integration.py` — 18 end-to-end tests

These test the R function's decision table. If the R drifts from the
Python, tests flag it. Run after any future edit to
`compute_group_validation()`:

```bash
cd ../../phase_9_classification/tests/
for t in test_compute_group_validation.py test_jackknife_semantics.py \
         test_phase4b_integration.py; do
    python3 "$t" || break
done
```

*(The tests target writer-side logic but currently live in phase_9/tests.
Candidate for a future move to `phase_7/validation/tests/` for
consistency — tracked as a low-priority item, not done in pass 16.)*

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

Or submit via the SLURM launcher:

```bash
sbatch LAUNCH_C01f_hypothesis.sh
```

## Diffing the pristine source

The pre-patch source is archived at
`../../../_archive_superseded/phase_7_validation_pristine_v934/_PRISTINE_v9.3.4.R`.
To see exactly what the 3 patches changed:

```bash
diff -u \
  ../../../_archive_superseded/phase_7_validation_pristine_v934/_PRISTINE_v9.3.4.R \
  STEP_C01f_hypothesis_tests.R
```

The only changes should be:
- Header rewrite (purely documentation)
- Insertion of `comp_from_registry` / `build_comp_with_fallback` /
  `compute_group_validation` before the T3 function
- Replacement of the per-candidate loop opening (lines ~1944–1976 in
  pristine) with the registry-aware version
- Replacement of the verdict row construction (lines ~2200–2249 in
  pristine) with the registry-aware version including new columns

## Pending

- **Schema update:** `hypothesis_verdict.schema.json` needs to bump to v2
  to include the new `family_linkage` enum
  (multi_family / few_family / single_family / pca_family_confounded /
  uncertain / unknown) and the new `quality_flags` field. See
  `../../phase_9_classification/patches/03_C01f_jackknife_semantics.R`
  line 162–167.
- **Real R parse on LANTA:** bracket-balance does not equal parse-valid.
  Run `Rscript --vanilla -e "invisible(parse(file='STEP_C01f_hypothesis_tests.R'))"`
  on LANTA before first submit.
- **One-candidate smoke test** to verify the comp_from_registry path
  actually reaches T1–T11 with the same results as the file-fallback
  path on a candidate that has both (registry groups + file comp).

## See also

- `../README.md` — phase 7 module root overview
- `../proposal/` — where the groups come from (seal writes
  `q6_validation_promotion_cap` that this sub-phase respects)
- `../../phase_9_classification/group_validation_gate.R` — reader-side
  gate (post pass 16)
- `../../phase_9_classification/patches/README.md` — the 3 patches, in
  application order
- `../../phase_9_classification/patches/_archive_superseded/README.md` —
  3 patches NOT to apply (two C01i + the obsolete flashlight C01f patch)
- `../../../docs/PHASE4_ARCHITECTURE.md` — the overall phase 4 → 9 design
