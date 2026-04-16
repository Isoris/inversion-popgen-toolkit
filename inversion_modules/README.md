# `inversion_modules/` — pipeline deployment root

Consolidated deployment tree for the catfish F₁ hybrid inversion pipeline.
Combines phase 2 discovery work (this session stream) with the v10 phase 4
catalog infrastructure (the `2026-04-16` bundle). All paths below are
relative to this directory.

## Status

**`✓ ALL GREEN`** from `python3 check_deployment_v10.py`:

- 26 expected files present
- All JSON schemas parse (22 in `registries/schemas/` + 4 in
  `phase_4_catalog/phase4b_rewrite/schemas/`)
- All R files bracket-balanced (12 native + 6 allowlisted, where the
  naive counter false-positives on valid R containing complex string
  literals)
- All Python files compile cleanly
- All schema block_types cross-reference correctly
- Registry tree creation works in tempdir
- All 5 automated test suites pass:
  - `test_compute_group_validation` (v10 — 22 unit tests)
  - `test_registry_sanity` (v10 — end-to-end plumbing)
  - `test_c01i_d_seal_resolution` (v10.1 — 21 resolution rules)
  - `test_jackknife_semantics` (v10.1.1 — 15 four-way jackknife tests)
  - `test_phase4b_integration` (v10.1 — 18 end-to-end phase 4b tests)

## Top-level layout

```
inversion_modules/
├── README.md                           ← you are here
├── CONFIG_ARCHITECTURE.md              how module configs source the master
├── PATH_WIRING.md                      canonical DELLY/Manta paths
├── check_deployment_v10.py             self-test (~10 s)
│
├── registries/                         v10 registry infrastructure (shared)
│   ├── api/
│   │   ├── R/registry_loader.R         unified R library
│   │   ├── python/registry_loader.py   Python mirror
│   │   └── bash/registry_loader.sh     shell helpers
│   └── schemas/structured_block_schemas/
│                                       18 Tier-2 block schemas
│
├── phase_2_discovery/                  genome scan for inversion-like regions
│   ├── 2c_precomp/                     C00 sv_prior + C01a precompute +
│   │                                   C01b_1 seeded regions + patches
│   │   ├── STEP_C00_build_sv_prior.R
│   │   ├── STEP_C01a_precompute.R
│   │   ├── STEP_C01b_1_seeded_regions.R
│   │   ├── patches/                    legacy integration patches (reference)
│   │   ├── README.md
│   │   └── RENAMING.md                 terminology migration tracker
│   │
│   └── 2d_detect/                      inv_detect v9.3 — matrix-based track
│       ├── run_all.R                   9-phase orchestrator
│       ├── STEP_D01…D17                (22 scripts, new naming)
│       ├── LAUNCH_*.slurm              (3 launchers)
│       ├── _archive_old_names/         original v9.3 docs (authoritative
│       │                               for label semantics)
│       └── README.md
│
├── phase_3_refine/                     breakpoint validation
│   └── MODULE_5A2_breakpoint_validation/
│                                       DELLY/Manta concordance, BND signal,
│                                       refined breakpoint coordinates
│
└── phase_4_catalog/                    per-candidate postprocessing
    ├── phase4_v10/                     base patches + 18 schemas + tests
    │   ├── docs/PHASE4_ARCHITECTURE_v10.md
    │   ├── patches/C01d /C01f /C01i/   registry-aware patches
    │   ├── orchestrator/run_phase4.sh  SLURM dependency chain
    │   └── tests/                      2 test suites, both passing
    │
    └── phase4b_rewrite/                v10.1.1 — four-script subsystem
                                        replacing old C01i
        ├── R/
        │   ├── STEP_C01i_decompose.R          4b.1 PCA + kmeans
        │   ├── STEP_C01i_b_multi_recomb.R     4b.2 recombinant detection
        │   ├── STEP_C01i_d_seal.R             4b.4 synthesis
        │   ├── engine_b_smoke_test.R          5-s self-test
        │   └── lib_decompose_helpers.R        shared utilities
        ├── python/
        │   ├── STEP_C01i_c_nested_composition.py  4b.3 composite detection
        │   └── nested_composition_core.py     (vendored)
        ├── schemas/                    4 schemas (v2 internal_dynamics +
        │                               recombinant_map + composition +
        │                               frequency.v2)
        ├── patches/                    C01f promotion_cap + jackknife
        ├── orchestrator/run_phase4b.sh SLURM DAG for the 4-script flow
        ├── tests/                      3 test suites, all passing
        ├── README.md
        ├── ADDITIONS_v10_1_1.md        what changed in v10.1.1
        └── docs/                       architecture + design notes
```

## The pipeline in one page

This follows the canonical phase layout from `SESSION_AUDIT_2026-04-16.md`.

```
phase_1_inputs/       (not in this deployment — lives elsewhere in the project)
  mask, ANGSD, BEAGLE, NGSadmix K=8, ancestry_bridge --prepare

phase_2_discovery/    scan genome for inversion-like regions
  2a_local_pca/       (legacy lives in MODULE_2A — not renamed yet)
  2b_mds/             (legacy lives in MODULE_2B — not renamed yet)
  2c_precomp/         ✓ HERE — C00 sv_prior, C01a precompute, C01b seeded regions
  2d_detect/          ✓ HERE — inv_detect v9.3 (matrix-based detection)
  2e_ghsl/            (not yet renamed — lives as MODULE-* in the legacy tree)

phase_3_refine/       breakpoint validation (bp resolution)
  ✓ HERE — MODULE_5A2_breakpoint_validation

phase_4_catalog/      per-candidate postprocessing — the main pipeline
  ✓ HERE — phase4_v10/ + phase4b_rewrite/
  Sub-phases are runtime states, not folders:
    4a  C01d scoring, C01e figures, C01g boundaries       (group_validation: NONE)
    4b  C01i decompose/multi_recomb/nested_comp/seal      (writes: UNCERTAIN)
    4c  C01f hypothesis tests                             (→ SUPPORTED/VALIDATED or SUSPECT)
    4d  cheat6/20/24, Q5 age, Q6 burden                   (reads: ≥ SUPPORTED)
    4e  classify_inversions, characterize_candidate       (reads: everything)

phase_5_followup/     (not in this deployment)
phase_6_secondary/    (not in this deployment)
```

## The two detection tracks in phase 2

Both tracks detect candidate inversion blocks. Both feed C01d (phase 4a)
via the legacy `triangle_intervals.tsv.gz` format.

| Track | Folder | Signal | Seed / start |
|---|---|---|---|
| Seed-based | `2c_precomp/STEP_C01b_1_seeded_regions.R` | per-window `max_abs_z` + `inv_likeness` + adaptive β-threshold | **MDS z-score outliers** |
| Matrix-based | `2d_detect/run_all.R` (9 phases) | similarity matrix structure + 6 image-processing transforms + NN persistence | step-down votes across row profiles |

They are designed to fail independently on different kinds of noise.
Candidates visible to both tracks are higher confidence than those
visible to one.

## Group validation — the phase 4 spine

Phase 4 is organized around **validation levels that gate downstream
cheats**. Groups (HOM_REF, HET, HOM_INV, RECOMBINANT) are proposed by
C01i (phase 4b), validated by C01f (phase 4c), and consumed by
group-dependent cheats (phase 4d):

```
NONE        phase 4a — scoring / figures / boundaries; no groups required
UNCERTAIN   phase 4b — groups proposed from PCA; composite intervals capped here
SUPPORTED   phase 4c — family-restricted + multi-family-contributing + few-family
VALIDATED   phase 4c — robust multi-family verdicts only
SUSPECT     phase 4c — only for pca_family_confounded
```

The four-way jackknife semantics in v10.1.1 replaces the old "demote
everything that isn't robust_multi_family" rule with a principled
five-way classification where **single-family-fragile = REAL
(family-restricted polymorphism)**, not SUSPECT. See
`phase_4_catalog/phase4b_rewrite/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`
for the rationale.

## The four key axes written by phase 4b

Phase 4b writes four things into each candidate's evidence registry:

```
internal_dynamics.json             per-window class track, silhouette,
                                   BIC, cluster quality
recombinant_map.json               three-signal recombinant detection
                                   + GC/DCO/suspicious classification
internal_ancestry_composition.json composite_flag ∈ {clean, maybe_composite,
                                                     likely_composite,
                                                     unknown_no_engine_b}
frequency.v2.json                  q6_group_validation
                                   q6_validation_promotion_cap
                                   q6_family_linkage        ← see "Still owed" §1
                                   q6_polymorphism_class     ← see "Still owed" §1
```

## Still owed (from the audit's blocking-items list)

1. **10-line edit to `STEP_C01i_d_seal.R`** to write the new
   `family_linkage` + `polymorphism_class` fields into the frequency
   block. The `compute_group_validation()` function (patched in v10.1.1)
   returns them but seal doesn't currently persist them. See
   `ADDITIONS_v10_1_1.md` section "Stage 2" for the exact location.
2. **C01f call-site updates**. `compute_group_validation()` now returns
   a list (`$level`, `$quality_flags`, `$family_linkage`) instead of a
   string. Call sites need the 3-line update.
3. **`ancestry_bridge.R --prepare`** for all 28 chromosomes before the
   first C01a run. This produces the `local_Q/<chr>.local_Q_samples.tsv.gz`
   cache that `nested_composition` reads.
4. **Full `Rscript --vanilla -e "parse(file='...')"`** on LANTA for
   every R file. The deployment check does bracket balance only; 6 of
   our files are allowlisted as "naive counter false-positives". Those
   need a real parse before SLURM submit.
5. **C01f comp-call-site patch** from `phase4_v10/patches/C01f/` —
   only at the main-loop k-means (not T4/T5/T6 internal calls).

## Still owed (from our phase 2 work)

See `phase_2_discovery/2c_precomp/RENAMING.md` section 10 for the full
list of scripts that still contain legacy `flashlight` / `snake` /
`cheatN` / `core_family` terminology. The v10 phase 4 code intentionally
still uses `cheatNN` and `event_class`-style naming — those stay until a
cross-module rename pass (post-manuscript).

Specific pending items:

1. **Five remaining patch files** at
   `phase_2_discovery/2c_precomp/patches/` have the autogen R parse bug
   (statements inside `c(...)`). Need rewriting using the clean pattern
   from integrated C01a.
2. **`PATH_WIRING.md` diff** to centralise DELLY/Manta paths in the
   master `00_inversion_config.sh`. Needs `readlink -f` verification on
   LANTA first to confirm which directory names are real vs symlinks.
3. **`snake_regions_multiscale/` directory rename** on LANTA — touches
   9 config variables + MODULE_5B/5C/5D/5E + existing output dirs.
   Coordinated via symlinks.
4. **Folder conflict**: `phase_2_discovery/2d_cores/` (the seed-based
   track's merge/community/consensus scripts from the legacy tree)
   still needs renaming. `RENAMING.md` section 13 proposes
   `2d_seeded_extension/` + this folder staying `2d_detect/`, OR
   `2d_cores/` → keep and this folder becomes `2e_matrix_detect/`.
   Your call.

## Scope note — what this deployment does NOT do

- **Does not replace `00_inversion_config.sh`** (the master config
  lives in the inversion codebase, not here). This tree is the module
  deployment; the master config is referenced via `INV_CONFIG` env
  var. See `CONFIG_ARCHITECTURE.md`.
- **Does not rebuild phase 1**. That's `ancestry_bridge.R --prepare`
  and the ANGSD/BEAGLE/NGSadmix pipeline, which live outside this tree.
- **Does not replace C01d or C01g**. Phase 4 uses patches (diff-apply
  style) against the existing scripts, not replacements. See each
  patch file's comment block.
- **Does not split composite intervals**. `likely_composite` stays as
  one catalog entry capped at UNCERTAIN. Multi-system catalog is
  v10.2 (post-manuscript).

## How to deploy

```bash
# 1. Drop into your toolkit
cd inversion-popgen-toolkit/inversion_modules/
tar xzf phase_2_plus_v10_bundle.tar.gz --strip-components=1

# 2. Verify
python3 check_deployment_v10.py
# Should print "✓ ALL GREEN"

# 3. Run R parse on LANTA (the 6 allowlisted files especially)
for f in phase_2_discovery/*/STEP_*.R phase_2_discovery/*/run_all.R \
         phase_4_catalog/**/*.R; do
  Rscript --vanilla -e "invisible(parse(file='$f'))" \
    || { echo "PARSE FAIL: $f"; break; }
done

# 4. Run ancestry_bridge --prepare for each chromosome (phase 1 prep)
for chr in C_gar_LG{01..28}; do
  Rscript utils/ancestry_bridge.R --prepare --K 8 --chr $chr
done

# 5. Phase 4 smoke test
bash phase_4_catalog/phase4b_rewrite/orchestrator/run_phase4b.sh \
  --dry-run --chroms LG12
```

## Referenced external docs

- `SESSION_AUDIT_2026-04-16.md` — uploaded earlier, contains the
  phase-structure decisions this deployment implements
- `phase_2_discovery/2c_precomp/RENAMING.md` — every rename already
  applied + full audit commands for remaining work
- `phase_2_discovery/2c_precomp/README.md` — precompute workflow
- `phase_2_discovery/2d_detect/README.md` — matrix-based detection
- `phase_4_catalog/phase4_v10/docs/PHASE4_ARCHITECTURE_v10.md` — the
  authoritative phase 4 design (~350 lines, 8 sections)
- `phase_4_catalog/phase4b_rewrite/docs/PHASE4B_REWRITE_ARCHITECTURE.md` —
  the four-script subsystem design
- `phase_4_catalog/phase4b_rewrite/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md` —
  why K=8, what local_Q is, why single_family_fragile is not SUSPECT
- `CONFIG_ARCHITECTURE.md` — how module configs source the master
- `PATH_WIRING.md` — the DELLY/Manta canonical paths proposal
