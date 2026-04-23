# `4d_group_proposal/` — propose genotype groups (UNCERTAIN)

Four scripts that together replace the old single
`STEP_C01i_decomposition_rewired_24_v934_registry.R`. Each candidate
interval gets decomposed into genotype groups, recombinants are
detected and classified, and the final output is a registered group
set plus flags at `q6_group_validation = UNCERTAIN`. Validation
happens in phase 4c.

## The four scripts

```
4b.1    STEP_C01i_decompose.R              PCA + k-means clustering
                                             ↓ writes per_window_class.rds
                                               + internal_dynamics.json

4b.2    STEP_C01i_b_multi_recomb.R         three-signal recombinant detection
                                             reads per_window_class.rds
                                             ↓ writes recombinant_map.json

4b.3    STEP_C01i_c_nested_composition.py  nested ancestry composition
                                             reads local_Q cache from ancestry_bridge
                                             ↓ writes internal_ancestry_composition.json

4b.4    STEP_C01i_d_seal.R                 synthesis + group registration
                                             reads all 3 blocks above
                                             ↓ registers 4 groups in sample_registry
                                               writes q6_group_validation=UNCERTAIN
                                               writes q6_validation_promotion_cap
```

## Dependencies

```
                ┌─→ 4b.1 decompose ──→ 4b.2 multi_recomb ─┐
                │                                          ├─→ 4b.4 seal
                └─→ 4b.3 nested_composition ──────────────┘
```

4b.1 and 4b.3 run in parallel (Engine B is independent of the PCA
clustering). 4b.2 waits for 4b.1 (reads its `per_window_class.rds`).
4b.4 gates on all three.

See `../orchestrator/run_phase4b.sh` for the SLURM DAG.

## What the three upstream signals are

### Signal 1 — PCA k-means class track (4b.1)

Optional Cheat 1 (flashlight) seeding of k-means from SV anchors, if
the sv_prior cache has `sample_inv_states` for this candidate. Cheat 2
(het-DEL constraint) if internal_dels are present. Silhouette + BIC
reported as cluster quality metrics; the per-window class vector is
saved for 4b.2 to consume.

### Signal 2 — Clair3 phase switches (4b.2)

Reads Clair3 postprocess output at `$BASE/MODULE_4A_SNP_INDEL50_Clair3/`.
A recombinant appears as a phase switch within the candidate interval.
Silent no-op if Clair3 data isn't available — just reduces sensitivity.

### Signal 3 — Flashlight hemizygous segments (4b.2)

Hemizygous (het-DEL) segments from `sv_prior` that span part of the
candidate interval are a third independent recombinant signal. Catches
recombinants that Signals 1 and 2 miss due to PCA ambiguity.

### Combined decision (4b.2)

A sample is flagged as RECOMBINANT if:
  `(Signal 1 AND Signal 2) OR Signal 3 OR (Signal 1 + >100 kb mosaic)`

Each flagged sample gets an `event_class` from cheat24:
- `gene_conversion` — short mosaic (< 200 kb) near a boundary
- `double_crossover` — long mosaic (> 200 kb) in the interior
- `suspicious` — ambiguous / possible genotyping noise
- `ambiguous` — insufficient signal to classify

## Composite detection (4b.3)

Reads Engine B's local-Q cache (per-window × per-sample ancestry
proportions). Classifies each sample's ancestry structure within the
candidate:

| structure_type | meaning |
|---|---|
| `single_block` | one dominant Q component across the interval |
| `two_block_composite` | two distinct Q regions inside the interval |
| `gradient` | smooth Q shift across the interval |
| `fragmented` | many small shifts (noisy or recombinant) |

The per-candidate `composite_flag` comes from the distribution:

| composite_flag | when |
|---|---|
| `clean` | < 10% samples are two_block_composite |
| `maybe_composite` | 10-40% samples |
| `likely_composite` | > 40% samples |
| `unknown_no_engine_b` | local_Q cache absent |

## What 4b.4 seal does

- Reads all 3 blocks above
- Applies the 4 resolution rules (see `docs/PHASE4B_REWRITE_ARCHITECTURE.md`)
- Registers 4 groups (HOM_REF, HET, HOM_INV, RECOMBINANT) in
  `sample_registry/sample_groups.tsv`
- Optional subgroups: RECOMBINANT_GC / RECOMBINANT_DCO based on cheat24
  event_class
- Writes these flat keys to the evidence registry:
  - `q6_group_validation = "UNCERTAIN"`
  - `q6_validation_promotion_cap = "UNCERTAIN"` if
    `composite_flag == "likely_composite"` (else unset)
  - `q1_composite_flag` (clean / maybe / likely / unknown)
  - `q2_decomp_quality_flags` (comma-separated list)
  - `q6_family_linkage = "unknown"` (placeholder; C01f overwrites)
  - `q6_polymorphism_class = "unclassified"` (placeholder; C01f overwrites)
- Writes `seal_summary.tsv` with per-candidate counts and metrics

## Engine B smoke test

`engine_b_smoke_test.R` — 5-second self-test that confirms Engine B is
compiled and reachable before long computes. Usage in any R script:

```r
if (!nzchar(Sys.getenv("SKIP_SMOKE_TEST", ""))) {
  source(file.path(Sys.getenv("BASE"),
                    "inversion-popgen-toolkit/utils/engine_b_smoke_test.R"))
  if (!run_engine_b_smoke_test()) {
    stop("Engine B smoke test failed. Set SKIP_SMOKE_TEST=1 to bypass.")
  }
}
```

See `../patches/ENGINE_B_SMOKE_TEST_INSERTS.md` for where to wire it
into C01a, cheat6, and Python.

## Optional upstream dependencies (graceful fallback)

All of these are optional — if absent, the pipeline still runs but
with reduced sensitivity:

- **Engine B local-Q cache** at `$Q_CACHE_DIR` — 4b.3 writes
  `composite_flag = unknown_no_engine_b` if absent
- **Flashlight (sv_prior) cache** at
  `$BASE/flashlight_v2/cache/` — Cheats 1 and 2 skipped if absent
- **Clair3 postprocess** at
  `$BASE/MODULE_4A_SNP_INDEL50_Clair3/postprocess_results/` — Signal 2
  skipped if absent
- **cheat24** at `$BASE/flashlight_v2/cheats/cheat24_recombinant_prior.R` —
  inline fallback available

## Shared utilities

`lib_decompose_helpers.R` — shared utilities used by 4b.1, 4b.2, 4b.4:
- Registry-aware path resolution
- K-means + silhouette helpers
- Registry write wrappers that handle both standalone and unified mode

## Tests

All tests live in `../tests/`, three of them cover this sub-phase:

- `test_c01i_d_seal_resolution.py` — 21 resolution rule tests
- `test_jackknife_semantics.py` — 15 tests (cross-phase, relevant here
  because 4b.4 writes the cap that 4c respects)
- `test_phase4b_integration.py` — 18 end-to-end tests

All pass as of last deployment check.

## See also

- `../docs/PHASE4B_REWRITE_ARCHITECTURE.md` — full design (400 lines, 6 sections)
- `../docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md` — rationale for
  K=8, what local_Q is, corrected jackknife semantics
- `../docs/ADDITIONS_v10_1_1.md` — what changed from v10.1 to v10.1.1
- `../orchestrator/run_phase4b.sh` — the SLURM DAG
- `../patches/README.md` — the three C01f patches applied in phase 4c
