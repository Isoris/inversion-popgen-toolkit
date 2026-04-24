# Phase 4b Rewrite v10.1.1 — Additions

This file documents what was added on top of the v10.1 delivery
(phase4b_rewrite_delivery_2026-04-16.tar.gz). The additions address three
operational questions:

1. How do we fail fast when Engine B is broken?
2. Is ancestry_bridge the right precompute step, and where does it fit?
3. How do we make the jackknife verdict work for family-restricted
   polymorphisms instead of wrongly rejecting them as SUSPECT?

**All 54 tests pass** (21 resolution + 18 integration + 15 jackknife
semantics).

## New files

```
R/engine_b_smoke_test.R                         5-second self-test
                                                 for Engine B reachability
patches/ENGINE_B_SMOKE_TEST_INSERTS.md          how to wire the smoke test
                                                 into C01a, cheat6, Python
patches/C01f_jackknife_semantics_patch.R        richer compute_group_validation
                                                 with four-way jackknife
                                                 classification
schemas/frequency.v2.schema.json                 5-class family_linkage +
                                                 polymorphism_class axis
docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md  design rationale for K=8,
                                                 what local_Q cache is,
                                                 corrected jackknife semantics
tests/test_jackknife_semantics.py                15 new unit tests (all pass)
```

## Modified files

```
python/STEP_C01i_c_nested_composition.py        added Q cache smoke test
                                                 (reads the first parent's
                                                 chromosome, checks
                                                 assigned_pop column,
                                                 warns and falls through
                                                 to stub path on failure)
```

## What to apply in what order

### Stage 1 — wire the smoke test

1. Drop `R/engine_b_smoke_test.R` into your toolkit `utils/` directory.
2. At the top of `STEP_C01a_snake1_precompute_*.R`, right after the
   `source(load_bridge.R)` line, add:

   ```r
   if (!nzchar(Sys.getenv("SKIP_SMOKE_TEST", ""))) {
     source(file.path(Sys.getenv("BASE", ""),
                       "inversion-popgen-toolkit/utils/engine_b_smoke_test.R"))
     if (!run_engine_b_smoke_test()) {
       stop("[C01a] Engine B smoke test failed — aborting before long compute. ",
            "Override with SKIP_SMOKE_TEST=1 to bypass.")
     }
   }
   ```

3. Same for `ancestry_bridge.R --prepare` mode.
4. Python wrapper is already patched — the smoke test runs automatically
   unless `SKIP_SMOKE_TEST=1` is set.
5. Bash orchestrator: add at the top of `run_phase4b.sh`:

   ```bash
   if [[ -z "${SKIP_SMOKE_TEST:-}" ]]; then
     RUN_SMOKE_TEST=1 Rscript "${MODULE_DIR}/utils/engine_b_smoke_test.R" \
       || { echo "Engine B not ready; see DESIGN_NOTE for fixes"; exit 1; }
   fi
   ```

### Stage 2 — apply the jackknife semantics patch

1. Replace `compute_group_validation()` in your patched C01f with the
   version in `patches/C01f_jackknife_semantics_patch.R`.
2. Update the call site — the function now returns a list with
   `$level`, `$quality_flags`, `$family_linkage` instead of just a
   string. See the patch comments for the 3-line call-site change.
3. Replace `schemas/frequency.schema.json` with
   `schemas/frequency.v2.schema.json` (the v2 schema adds the new
   `polymorphism_class` axis and the 5-valued `family_linkage` enum).
4. In `STEP_C01i_d_seal.R`, the family_linkage value from
   `compute_group_validation` gets written into the frequency block.
   This is a ~10-line addition to seal; not in this delivery yet.
   Suggested location: in the same spot seal writes
   `q6_group_validation`, also write `q6_family_linkage` and
   `q6_polymorphism_class`.

### Stage 3 — run ancestry_bridge first

1. In your pipeline dispatcher (the overall run-all, not run_phase4b.sh),
   add a Phase 1a step before C01a:

   ```bash
   Rscript utils/ancestry_bridge.R --prepare --K 8 --chr <CHR> \
       || { echo "ancestry_bridge --prepare failed"; exit 1; }
   ```

   This writes `${PRECOMP_OUT}/local_Q/<CHR>.local_Q_samples.tsv.gz` and
   its companion summary files. C01a picks them up automatically via
   `merge_local_Q_into_invlikeness()`. nested_composition reads the
   same cache directly.

2. The smoke test runs twice in this flow — once for ancestry_bridge
   before --prepare does work, once for C01a. Both are cheap. The
   second instance is ~instant because Engine B is already warm.

## The four-way jackknife in a table

| verdict (cheat6)          | family_linkage         | validation  | rationale |
|---------------------------|------------------------|-------------|-----------|
| robust_multi_family       | multi_family           | can VALIDATE | cohort-wide polymorphism |
| multi_family_contributing | multi_family           | SUPPORTED   | a bit family-driven but real |
| few_family                | few_family             | cap SUPPORTED | lineage-restricted |
| single_family_fragile     | single_family          | cap SUPPORTED | family-specific polymorphism — **real** |
| pca_family_confounded     | pca_family_confounded  | SUSPECT     | PCA groups are ancestry, not arrangement |

The last row is the ONLY one that still becomes SUSPECT. The previous
v10 rule demoted all four "non-robust" verdicts, which wrongly threw
away family-restricted inversions (the same inversions that are most
likely to be karyotype-changing line-specific events).

## What's NOT in this delivery

- **Cross-K stability analysis** (the "run at K=5, K=7, K=8, K=12 and
  check hierarchical coherence" proposal). This is a separate ~200-line
  module (`STEP_C01i_e_cross_K_stability.R`) that reads pre-computed Q
  matrices at multiple K and checks whether K=8 groups map cleanly onto
  K=12 groups. Useful for refining the family_linkage classification of
  family-restricted candidates. Deferred to v10.2.

- **The 10-line addition to STEP_C01i_d_seal.R** that writes the new
  family_linkage + polymorphism_class into the frequency block. Easy
  but not yet written. Do this when applying stage 2.

- **A bash mirror of engine_b_smoke_test.R** (pure shell, for pipelines
  that don't call R). Not needed yet since all your Engine B entry
  points are R.
