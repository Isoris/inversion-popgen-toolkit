# `phase_4_postprocessing/` — per-candidate analysis and classification

Phase 4 takes candidate inversions from phase 2 (discovery) + phase 3
(breakpoint refinement) and turns them into a classified catalog. It's
the main per-candidate analysis phase: existence layers are collected,
genotype groups are proposed, groups are validated against independent
evidence, group-dependent cheats are run on validated groups only, and
a final 14-axis classification is assigned.

## Layout

Flat folder, with letter subfolders for each sub-phase:

```
phase_4_postprocessing/
├── README.md                              ← you are here
│
├── 4a_existence_layers/                   existence evidence (no groups)
│                                            ├ C01d pass-1 scoring
│                                            ├ C01e candidate figures
│                                            └ C01g boundary catalog
│                                          [placeholder — upload your scripts]
│
├── 4b_group_proposal/                     propose groups → UNCERTAIN
│   ├─ STEP_C01i_decompose.R               4b.1 PCA k-means
│   ├─ STEP_C01i_b_multi_recomb.R          4b.2 recombinant detection
│   ├─ STEP_C01i_c_nested_composition.py   4b.3 composite detection
│   ├─ nested_composition_core.py          (vendored helper)
│   ├─ STEP_C01i_d_seal.R                  4b.4 synthesis + registration
│   ├─ lib_decompose_helpers.R             (shared utilities)
│   ├─ engine_b_smoke_test.R               (5-s self-test)
│   └─ README.md
│
├── 4c_group_validation/                   validate → SUPPORTED/VALIDATED/SUSPECT
│                                          [placeholder — upload your C01f]
│                                          apply patches/01 → 02 → 03 in order
│
├── 4d_group_dependent/                    gated cheats (≥ SUPPORTED required)
│                                            cheat6, cheat20, cheat30, Q5, Q6
│                                          [placeholder — upload your cheats]
│
├── 4e_final_classification/               14-axis final tag
│                                            classify_inversions, characterize_candidate
│                                          [placeholder — upload your scripts]
│
├── orchestrator/                          SLURM dependency chains
│   ├─ run_phase4.sh                       full phase 4 orchestrator (EDIT: swap 4b section)
│   ├─ run_phase4b.sh                      4b sub-DAG
│   └─ LAUNCH_group_cheats_example.sh      example of validation-gated launcher
│
├── patches/                               apply in numeric order to your C01f
│   ├─ 01_C01f_comp_from_registry.R        v10 base patch
│   ├─ 02_C01f_promotion_cap.R             v10.1 promotion cap
│   ├─ 03_C01f_jackknife_semantics.R       v10.1.1 four-way jackknife
│   ├─ ENGINE_B_SMOKE_TEST_INSERTS.md      how to wire smoke test
│   ├─ README.md                           application order + semantics
│   └─ _archive_superseded/                DO NOT APPLY — see its README
│       ├─ C01i_cheat24_merge_patch.R      superseded by 4b.2
│       └─ C01i_recombinant_groups_patch.R superseded by 4b.4 seal
│
├── schemas/                               Tier-2 JSON schemas for this phase
│   ├─ internal_dynamics.schema.json       (v2 — replaces v10's)
│   ├─ recombinant_map.schema.json         (from 4b.2)
│   ├─ internal_ancestry_composition.schema.json (from 4b.3)
│   └─ frequency.v2.schema.json            (v10.1.1 with family_linkage)
│
├── tests/                                 all 5 test suites, all passing (77 tests)
│   ├─ test_compute_group_validation.py    22 tests (v10)
│   ├─ test_registry_sanity.py             end-to-end plumbing (v10)
│   ├─ test_c01i_d_seal_resolution.py      21 tests (v10.1)
│   ├─ test_jackknife_semantics.py         15 tests (v10.1.1)
│   └─ test_phase4b_integration.py         18 tests (v10.1)
│
└── docs/                                  architecture + design notes
    ├─ PHASE4_ARCHITECTURE.md              ~350 lines, 8 sections (was v10)
    ├─ PHASE4B_REWRITE_ARCHITECTURE.md     ~400 lines, 6 sections (v10.1)
    ├─ DESIGN_NOTE_K_and_local_Q_and_jackknife.md   (v10.1.1)
    └─ ADDITIONS_v10_1_1.md                v10.1.1 changelog
```

## Data flow

```
FROM PHASE 2D:  triangle_intervals.tsv.gz  (candidate coords + basic shape)
FROM PHASE 3:   refined breakpoints + SV anchors + GHSL haplotype contrast

        ▼
┌──────────────────────────────────────────────────────────────────────┐
│ 4a_existence_layers/                                                 │
│   C01d pass-1: aggregate Layers A+B+C → existence_score              │
│   C01e:        per-candidate figures                                 │
│   C01g:        boundary catalog from phase 3                         │
│   Writes: existence_layer_{a,b,c}.json + boundary.json               │
│   State: group_validation = NONE                                     │
└──────────────────────────────────────────────────────────────────────┘
        ▼
┌──────────────────────────────────────────────────────────────────────┐
│ 4b_group_proposal/                                                   │
│   ┌─→ 4b.1 decompose ──→ 4b.2 multi_recomb ─┐                        │
│   │                                          ├─→ 4b.4 seal           │
│   └─→ 4b.3 nested_composition ──────────────┘                        │
│   Writes: internal_dynamics.json, recombinant_map.json,              │
│           internal_ancestry_composition.json, frequency.v2.json      │
│   Registers: HOM_REF, HET, HOM_INV, RECOMBINANT groups               │
│   State: group_validation = UNCERTAIN                                │
│          (composite → promotion_cap = UNCERTAIN)                     │
└──────────────────────────────────────────────────────────────────────┘
        ▼
┌──────────────────────────────────────────────────────────────────────┐
│ 4c_group_validation/                                                 │
│   C01f: run T1-T10 hypothesis tests (Layer D)                        │
│   Four-way jackknife classification → family_linkage                 │
│   Respects promotion_cap from 4b                                     │
│   Writes: existence_layer_d.json, hypothesis_verdict.json            │
│   State: UNCERTAIN → SUPPORTED | VALIDATED | SUSPECT                 │
└──────────────────────────────────────────────────────────────────────┘
        ▼
┌──────────────────────────────────────────────────────────────────────┐
│ 4d_group_dependent/  (ONLY if q6_group_validation ≥ SUPPORTED)       │
│   cheat6, cheat20, cheat24, cheat30, Q5 age, Q6 burden               │
│   Writes: age_evidence.json, burden.json, mechanism.json             │
└──────────────────────────────────────────────────────────────────────┘
        ▼
┌──────────────────────────────────────────────────────────────────────┐
│ 4e_final_classification/                                             │
│   classify_inversions.R: assign 14-axis tag                          │
│   characterize_candidate.R: per-candidate summary                    │
│   Writes: final_catalog.tsv + final_classification.json per candidate│
└──────────────────────────────────────────────────────────────────────┘
        ▼
TO PHASE 5: per-candidate figures, manuscript panels
```

## Deployment order

```bash
# 1. Apply the 3 C01f patches in order (requires your current C01f)
# See patches/README.md

# 2. Run the 5 test suites to confirm nothing regressed
cd tests/
for t in *.py; do python3 "$t" || break; done

# 3. Deploy orchestrator — edit run_phase4.sh to call run_phase4b.sh
#    for the 4b block (the orchestrator currently references an older
#    4b decomposition call that no longer exists)

# 4. Smoke test on one chromosome
bash orchestrator/run_phase4.sh --dry-run --chroms LG12
```

## What's currently ready vs pending

| Sub-phase | Status |
|---|---|
| 4a | placeholder — C01d / C01e / C01g need uploading |
| 4b | **ready** — 4 scripts + 4 schemas + 3 tests |
| 4c | placeholder — C01f needs uploading + 3 patches applied |
| 4d | placeholder — cheats + Q5 + Q6 need uploading |
| 4e | placeholder — classify_inversions + characterize_candidate need uploading |

Patches are ready. Tests pass. Schemas validate. Registry infrastructure
at `../registries/` is ready. The skeleton is in place — just needs
the scripts for the 4 other sub-phases.

## See also

- `../README.md` — top-level synthesis of all phases
- `../registries/` — registry infrastructure (shared)
- `../phase_2_discovery/2d_detect/STEP_D12_bridge_to_C01d.R` — the
  contract between phase 2 and phase 4a
- `../phase_3_refine/` — Breakpoint validation + BND rescue. This is the
  **writer for Layer D** (`q7_layer_d_*` flat keys via the
  `existence_layer_d` registry block) and for **supplementary Layer B**
  (`q7b_bnd_*` keys via `existence_layer_b_bnd_rescue`). The Layer D
  VALIDATED promotion gate (fisher_p<0.05 AND fisher_or>5) in 4c
  `compute_group_validation()` and the Layer B BND-rescue consumer in
  4a Layer B scoring both depend on phase 3 having run. Earlier
  sessions believed the registry contract was inactive; the chat-5
  audit wired STEP03 and STEP06 through the Python registry API to
  make it real.
