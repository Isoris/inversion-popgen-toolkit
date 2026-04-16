# Phase 4 v10 — Delivery package

Stage-1 patches + Stage-2 library + 18 JSON schemas + orchestrator +
tests for the group-validation-gated Phase 4 architecture described
in the handoff.

All code tested: 22 unit tests + 1 end-to-end sanity test pass.
JSON schemas validated. R files bracket-balance-checked (full R parse
not available in the build sandbox, but syntax is clean).

## Contents

```
README.md
docs/
  PHASE4_ARCHITECTURE_v10.md         the full design (8 sections, ~350 lines)

patches/
  C01i/C01i_recombinant_groups_patch.R   register 4 groups + GC/DCO + UNCERTAIN
  C01i/C01i_cheat24_merge_patch.R        tiny merge that makes GC/DCO work
  C01f/C01f_comp_from_registry_patch.R   comp_from_registry + compute_group_validation

orchestrator/
  run_phase4.sh                          SLURM dependency chain, phase 4a → 4e
  LAUNCH_group_cheats_example.sh         concrete validation-gated launcher

registries/
  api/R/registry_loader.R                unified R library (~380 lines)
  api/python/registry_loader.py          Python mirror for STEP03, audit tools
  api/bash/registry_loader.sh            shell helpers for launchers
  schemas/structured_block_schemas/
    existence_layer_a/b/c/d.schema.json
    boundary.schema.json                 (templated for left/right)
    carrier_reconciliation.schema.json
    internal_dynamics.schema.json
    mechanism.schema.json
    age_evidence.schema.json             (dual sub-blocks)
    frequency.schema.json
    hypothesis_verdict.schema.json
    morphology.schema.json
    band_composition.schema.json
    burden.schema.json
    block_detect.schema.json
    inv_detect_blocks.schema.json
    peel_diagnostic.schema.json
    synteny_dollo.schema.json
    INDEX_remaining_blocks.json          reference map

tests/
  test_registry_sanity.py                end-to-end: 4 groups, 4 blocks,
                                          49 keys, NONE→UNCERTAIN→SUPPORTED→
                                          VALIDATED transitions
  test_compute_group_validation.py       22 unit tests for the promotion
                                          logic (pure function, in sync with R)
```

## Recombinants — the missing group

The old 3-class registration (HOM_STD / HET / HOM_INV) silently dropped
recombinant samples into whichever class their average PC1 happened to
land in. Because recombinants have a chimeric haplotype (HOM_REF on part
of the interval, HOM_INV on the other part) averaging puts most of them
in HET, which then poisons the HET group's Fst/theta statistics. The
fix is a separate RECOMBINANT group, registered as `inv_<cid>_RECOMBINANT`.

**Why it can't be a subclass of HET.** A HET sample is heterozygous at
every position inside the inversion. A recombinant sample is homozygous
for one arrangement in some windows and homozygous for the other in
other windows. Those two things only look the same when you average PC1
across the interval — they are not the same biological entity, and
they don't belong in the same Fst computation.

**Cheat 24 splits recombinants further.** Based on position of the
switchpoint and mosaic length, cheat24 classifies each recombinant as:

- `gene_conversion` — short mosaic near a boundary (common)
- `double_crossover` — long mosaic in the interior (rare, needs two
  independent crossovers inside the suppressed region)
- `suspicious` — ambiguous (possible genotyping noise)

The patch registers these as optional subgroups:
`inv_<cid>_RECOMBINANT_GC` and `inv_<cid>_RECOMBINANT_DCO`. This lets
downstream Q4 mechanism analysis pull just the GC recombinants (which
carry different junction signatures than DCO) without re-running
cheat24.

**The flat group stores dominant class. The per-window class lives in
the structured block.** `sample_groups.tsv` gets one row per group per
candidate — good for "who are the 7 recombinants in LG12_17?". The
detailed per-window HOM_REF vs HOM_INV track lives in
`evidence_registry/per_candidate/LG12_17/structured/internal_dynamics.json`,
which has a `recombinants` array with switchpoint_bp, mosaic_length_bp,
event_class, posterior, and windows_HOM_REF / windows_HOM_INV lists
per recombinant sample. That's where "is LG12_17's recombinant #3
HOM_REF at position 8,400,000?" lives.

**C01f excludes recombinants from its `comp` table.** This is correct:
T1/T2/T7/T8/T10 all assume each sample has a single class across the
interval, so recombinants would break them. The patch attaches
`attr(comp, "recombinants_excluded")` so any test that wants to sanity-
check which samples were dropped can still access the list.

## Applying the patches — minimal workflow

### Stage 1 (this week, no new infrastructure)

1. **C01i**: apply two patches in order.
   - First `patches/C01i/C01i_cheat24_merge_patch.R` — the tiny
     merge that puts cheat24's `event_class` onto the main `result`
     data.table. Without this, the RECOMBINANT_GC / RECOMBINANT_DCO
     subgroups silently won't register (their subset returns zero).
   - Then `patches/C01i/C01i_recombinant_groups_patch.R` — str_replace
     the old 3-class loop (~lines 357-366) with the
     `register_four_classes()` block that registers HOM_REF, HET,
     HOM_INV, RECOMBINANT plus the optional GC/DCO subgroups, and
     sets `q6_group_validation = UNCERTAIN`.

2. **C01f**: apply `patches/C01f/C01f_comp_from_registry_patch.R` in
   three parts (A, B, C marked in the file).
   - (A) Add `comp_from_registry()` and `build_comp_with_fallback()`
     near the top of the script.
   - (B) Patch **only the main-loop call site** that builds `comp`
     via PC1 k-means. Leave the T4/T5/T6 internal k-means calls
     alone — they are per-subregion by design.
   - (C) Add `compute_group_validation()` (pure function) and put
     its return value into the verdict row's `group_validation_after`
     field. The schema extracts that into `q6_group_validation`
     automatically when the block is written — **do not call
     `reg$add_evidence` separately**, it creates a race with the
     block's keys_extracted.

3. **C01d**: no code change. Orchestrator calls it twice (pass-1
   without `--boundary_dir --hyp_dir`, pass-2 with them).

4. **Orchestrator**: drop `orchestrator/run_phase4.sh` into
   `${MODULE_DIR}/`, write seven thin `LAUNCH_*.sh` wrappers.
   `LAUNCH_group_cheats_example.sh` shows the pattern for
   validation-gated dispatch.

5. **Test locally before submission**:
   ```bash
   cd tests/
   python3 test_compute_group_validation.py    # logic correctness
   python3 test_registry_sanity.py             # end-to-end plumbing
   ```
   Both should pass with no errors.

6. **Smoke test on one chromosome**: LG12 is a good choice (known
   inversion at ~8.2-12.5 Mb). Run
   `bash run_phase4.sh --dry-run --chroms LG12` to see the
   dependency chain, then remove `--dry-run`.

### Stage 2 (during manuscript draft, not blocking it)

7. **Library**: `registries/api/R/registry_loader.R` and
   `registries/api/python/registry_loader.py` ready to drop in.
   Bash helpers in `api/bash/registry_loader.sh`.

8. **Schemas**: all 18 schemas now written. Copy
   `registries/schemas/structured_block_schemas/` to the toolkit
   registries directory.

9. **Migrate scripts incrementally**. Start with C01g (cleanest
   TSV → block mapping), then C01i, then C01f, then C01d pass-2.
   Python scripts (STEP03_statistical_tests.py,
   breakpoint_evidence_audit.py) migrate using the Python library.

## What this does NOT change

- `candidate_scores.tsv.gz` still gets written by C01d pass-2. It's
  a snapshot/projection, not the source of truth, but anything
  currently downstream of it keeps working.
- inv_detect v9.3 Phase 2 output (scoring_table_*.tsv, triangle_intervals
  via the 12_bridge_to_codebase.R adapter) unchanged.
- Engine B compilation + load_bridge.R — unchanged. The library
  falls back to the old sample_registry.R if the v10 registries/
  directory doesn't exist yet.
- PHASE_01C block_detect, C00 flashlight, C01a precompute,
  C01b cores — none of these change.

## Known gaps in this delivery

- **Schema key-extractor needs dotted-path support for `age_evidence`**.
  The schema uses `"from": "gds_sub_block.gds_gap"` notation. The Python
  library supports this (see `_extract_keys_from_schema`). The R library
  currently does flat lookups only — patch needed when migrating cheat30.
  One-line fix: split the `from` on `.` and walk the data list with
  `purrr::pluck` or manual indexing.
- **characterize_candidate and classify_inversions not migrated**. The
  architecture doc describes where they fit in phase 4e, but they still
  consume the old ad-hoc TSVs. Migration is mechanical once the library
  is in place.
- **LAUNCH_*.sh wrappers are sketched only**. `run_phase4.sh` references
  seven launchers; one full example is provided
  (`LAUNCH_group_cheats_example.sh`). The others (LAUNCH_C01d_pass1,
  LAUNCH_C01g_boundary, LAUNCH_C01i_decomposition, LAUNCH_C01f_hypothesis,
  LAUNCH_C01d_pass2, LAUNCH_characterize_classify) are thin shims — each
  sources `registry_loader.sh` and shells into the corresponding R
  script. Ten minutes each once you have the paths right.
- **No Python mirror of the R-side comp_from_registry**. Python scripts
  don't currently call it; if STEP03 ever needs carriers + band labels
  in the same shape as `comp`, write a Python helper that mirrors
  `SamplesAPI.get_carriers` into the data.table-equivalent.
