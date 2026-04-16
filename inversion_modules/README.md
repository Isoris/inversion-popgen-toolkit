# `inversion_modules/` — pipeline deployment root

Consolidated deployment tree for the catfish F₁ hybrid inversion
pipeline (226 samples, *C. gariepinus × C. macrocephalus*). All paths
below are relative to this directory.

## Top-level layout

```
inversion_modules/
├── README.md                           ← you are here
├── CONFIG_ARCHITECTURE.md              how module configs source the master
├── PATH_WIRING.md                      canonical DELLY/Manta paths
├── check_deployment_v10.py             self-test (~5 s)
│
├── utils/                              shared helpers
│   ├── ancestry_bridge.R               phase 1 prep + local_Q cache
│   ├── local_q_projector.R             Engine B local-Q computer
│   ├── sample_map.R, sample_registry.R sample bookkeeping
│   └── theme_systems_plate.R           ggplot theming
│
├── phase_1_inputs/                     mask + ANGSD SAF/SFS + BEAGLE
│   ├── steps/   STEP01_mask_regions, STEP02_callable_mask, STEP03_mask_stats
│   ├── launchers/ LAUNCH_STEP04..07   ANGSD SAF, merge, SFS, SNP + BEAGLE
│   └── runners/  run_5A1_discovery_inputs.sh
│
├── phase_2_discovery/                  genome scan for inversion-like regions
│   ├── 2a_local_pca/                   dosage + per-chr local PCA windows
│   ├── 2b_mds/                         MDS (stage1 + stage2)
│   ├── 2c_precomp/                     SV prior + precompute + seeded regions + PHASE_01C
│   │   ├── STEP_C00_build_sv_prior.R
│   │   ├── STEP_C01a_precompute.R
│   │   ├── STEP_C01b_1_seeded_regions.R
│   │   ├── PHASE_01C_block_detect.R
│   │   ├── diags/, patches/
│   │   └── README.md + RENAMING.md
│   ├── 2d_candidate_detection/         staircase detector (primary boundary track)
│   │   ├── STEP_D01..D17 + run_all.R
│   │   ├── LAUNCH_*.slurm
│   │   ├── tests/test_staircase.R
│   │   └── README.md
│   └── 2e_ghsl/                        GHSL haplotype contrast (Layer C)
│
├── phase_3_refine/                     bp-resolution breakpoint validation
│   └── MODULE_5A2_breakpoint_validation/
│                                       DELLY/Manta concordance, BND signal
│
├── phase_4_postprocessing/             per-candidate postprocessing — the spine
│   ├── 4a_existence_layers/            ← catalog birth (C01d/C01e/C01g)
│   ├── 4b_group_proposal/              C01i decompose / multi_recomb / nested_comp / seal
│   ├── 4c_group_validation/            C01f hypothesis tests + gate
│   ├── 4d_group_dependent/             Q5 age + Q6 burden + cheat28/29/30
│   ├── 4e_final_classification/        characterize_candidate + compute_candidate_status
│   ├── docs/                           PHASE4 + PHASE4B architecture + design notes
│   ├── orchestrator/                   SLURM DAG (run_phase4b.sh)
│   ├── patches/                        C01f registry-aware patches
│   ├── schemas/                        4 Tier-2 block schemas
│   ├── specs/                          evidence registry specs
│   └── tests/                          5 test suites (all pass)
│
├── phase_5_followup/                   per-candidate deep analysis
├── phase_6_secondary/                  LD / Fst / HOBS secondary analyses
│   └── MODULE_5C_Inversion_LD/, MODULE_5D_Inversion_FST/, MODULE_5E_Inversion_HOBS/
│
└── _archive/                           legacy v8.5 HPC tree (kept for reference)
```

There is also `_archive_superseded/` at the deployment root for things
that were dropped from the pipeline (not legacy — dead):

```
_archive_superseded/
└── fuzzy_merge_abandoned/              STEP_C01b_2 + README explaining why
                                        it was dropped (overmerge problem)
```

## The pipeline in one page

```
phase_1_inputs/       mask + ANGSD SAF/SFS + BEAGLE + NGSadmix K=8
                      ancestry_bridge --prepare (produces local_Q cache)
                           │
                           ▼
phase_2_discovery/    scan genome for inversion-like regions
  2a_local_pca/       lostruct per-window PCA on dosage
  2b_mds/             MDS on lostruct distance matrices
  2c_precomp/         SV prior (C00) + heavy precompute (C01a)
                      + seeded regions (C01b_1) + PHASE_01C landscape
  2d_candidate_detection/  ★ staircase detector — primary boundary track
  2e_ghsl/            Clair3 phased-genotype haplotype contrast (Layer C)
                           │
                           ▼
phase_3_refine/       breakpoint validation (bp resolution)
                      DELLY/Manta concordance + BND signal
                           │
                           ▼
phase_4_postprocessing/  per-candidate postprocessing — the main pipeline
  4a  C01d scoring (★ CATALOG BIRTH), C01e figures, C01g boundaries
      group_validation: NONE
  4b  C01i decompose + multi_recomb + nested_comp + seal
      writes: UNCERTAIN
  4c  C01f hypothesis tests + gate
      → SUPPORTED / VALIDATED / SUSPECT
  4d  Q5 age + Q6 burden + cheat28/29/30
      reads: ≥ SUPPORTED
  4e  classify_inversions + characterize_candidate
      reads: everything
                           │
                           ▼
phase_5_followup/     per-candidate deep analysis (MODULE_5B successor)
phase_6_secondary/    LD / Fst / HOBS secondary analyses
```

## Catalog birth — where candidates become first-class objects

`phase_4/4a/STEP_C01d_candidate_scoring` is the pipeline's pivot. At
the moment it runs:

1. Reads the primary boundary track (staircase blocks from
   `phase_2/2d_candidate_detection/`, bridge-converted to
   `triangle_intervals.tsv.gz`)
2. Reads evidence from whichever of these are available, as scoring
   dimensions (not gates):
   - `--cores_dir` → seeded regions from `phase_2/2c_precomp/STEP_C01b_1`
   - `--boundary_dir` → `boundary_catalog_unified.tsv.gz` from C01g
     (built from 5 boundary sources)
   - `--hyp_dir` → hypothesis verdicts from C01f (when available — a
     second C01d pass after 4c)
3. Computes 12 scoring dimensions D1–D12, assigns Tier 1/2/3/4
4. Writes `candidate_scores.tsv.gz` — **the catalog**

Still owed: a `create_candidate_folders.sh` that materialises
per-candidate folder structure (`01_detection/`, `02_genotypes/`,
`04_breakpoints/`, `05_mechanism/`, `06_evolution/`) from the catalog.
Previous chats describe this running right after C01d; it is not in
the current deployment. Downstream (`phase_4/4b`..`4e`) still works
without it — the scripts read the flat catalog directly — but the
per-candidate folder layout that older launchers expect is absent.

## The two tracks in phase 2

Both tracks produce candidate blocks. They are designed to fail on
different kinds of noise — a candidate visible to both is higher
confidence than one visible to only one.

| Track | Primary signal | Seed / start | Folder |
|---|---|---|---|
| Matrix-based (staircase) | vote-based boundaries across row profiles of sim_mat | step-down votes across row profiles | `2d_candidate_detection/` |
| Seed-based (seeded regions) | MDS z-outlier extension under damage budget | MDS z-score outliers | `2c_precomp/STEP_C01b_1_seeded_regions.R` |

The matrix track is **primary** — it runs `run_all.R` (9 phases) per
chromosome and converts its output to the legacy
`triangle_intervals.tsv.gz` format that C01d reads. The seed track is
internal evidence — C01d consumes its output via `--cores_dir` as one
of its scoring dimensions.

### Historical note — no merge step in phase 2

An earlier design had a `2d_seeded_merge/` folder with
`STEP_C01b_2_merge.R` that used 1D fuzzy max-min composition to
consolidate seeded regions before scoring. It was dropped because it
overmerged across real boundaries visible in the 2D sim_mat (LG19, LG28
during testing). The replacement is: two parallel boundary detectors
(staircase in `2d_candidate_detection/`, PHASE_01C in `2c_precomp/`)
feeding phase 4 directly. The retired script is at
`_archive_superseded/fuzzy_merge_abandoned/` with its own README.

## Phase 4 — validation levels

Phase 4 is organised around validation levels that gate downstream
tests. Groups (HOM_REF, HET, HOM_INV, RECOMBINANT) are proposed by
C01i (phase 4b), validated by C01f (phase 4c), and consumed by
group-dependent tests (phase 4d):

```
NONE        phase 4a — scoring / figures / boundaries; no groups required
UNCERTAIN   phase 4b — groups proposed from PCA; composite intervals cap here
SUPPORTED   phase 4c — family-restricted or few-family verdicts
VALIDATED   phase 4c — robust multi-family verdicts
SUSPECT     phase 4c — pca_family_confounded only
```

The four-way jackknife semantics (v10.1.1) treats
`single_family_fragile` as REAL (family-restricted polymorphism), not
SUSPECT. See `phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`.

## The 4-layer independence framework

| Layer | Source | Algorithm | Produced by |
|---|---|---|---|
| A | dosage (genotype likelihoods) | lostruct local PCA → MDS | `phase_2/2a` + `2b` |
| B | SV caller output | DELLY2 + Manta catalogs | `phase_2/2c/STEP_C00` |
| C | Clair3 phased genotypes | GHSL haplotype contrast | `phase_2/2e_ghsl` |
| D | genotype–breakpoint association | Fisher odds ratio linking A, B, C | `phase_4` |

Layers are designed to fail independently. Tier 1 = multi-layer
convergence; Tier 3–4 = single-layer evidence or layer-specific
artefact.

## Deployment check

```bash
python3 check_deployment_v10.py
```

Reports on files present, JSON schemas parsing, R bracket balance,
Python syntax, schema cross-references, registry loader smoke, and
runs the 5 test suites in `phase_4_postprocessing/tests/`.

The `registries/` API + schema directory lives at the repo root
(`inversion-popgen-toolkit/registries/`), not inside this tree — the
deployment-check errors about missing registry files are expected in
this configuration. The 5 tests all pass.

## Still owed

1. Real `Rscript --vanilla -e 'parse(file=...)'` on LANTA for every R
   file, especially the allowlisted ones (naive bracket counter false-
   positives on ggplot dict-arg + complex string literals).
2. Bump `hypothesis_verdict.schema.json` to v2 (add `family_linkage`
   enum + `quality_flags` field) to match v10.1.1 semantics.
3. One-chromosome smoke test through 2c → 2d → phase_3 → 4a → 4b → 4c
   (suggested target: LG12 with the Tier-1 candidate at 52.4–54.6 Mb
   from earlier manuscript drafts).
4. `ancestry_bridge.R --prepare` for all 28 chromosomes before the
   first C01a run (produces `local_Q/<chr>.local_Q_samples.tsv.gz`
   cache read by `nested_composition`).
5. Wire C01f call sites to the new `compute_group_validation()` list
   return (`$level`, `$quality_flags`, `$family_linkage`) instead of
   the pre-v10.1.1 string return.
6. `create_candidate_folders.sh` — the script that materialises
   per-candidate folder structure (`01_detection/`, `02_genotypes/`,
   `04_breakpoints/`, `05_mechanism/`, `06_evolution/`) from
   `candidate_scores.tsv.gz`. Described in past chats, not present in
   this deployment. Low priority — phase 4b–4e consume the flat catalog
   directly.
7. Cross-module rename pass (post-manuscript) to retire remaining
   `cheatNN`/`snake`/`core` terminology in phase_4 and MODULE_5B–E
   consumers. See `2c_precomp/RENAMING.md` section 13 for the
   up-to-date remaining list. Note: C01d's output column names
   `d12_snake_concordance` and `snake_overlap` are consumed by
   `4e/compute_candidate_status.R` and `test_registry_sanity.py` —
   renaming them is a coordinated change across those three files.

## Referenced docs

- `CONFIG_ARCHITECTURE.md` — module configs source the master
- `PATH_WIRING.md` — DELLY/Manta canonical paths
- `phase_2_discovery/2c_precomp/README.md` — precompute + PHASE_01C
  workflow + "why no merge step"
- `phase_2_discovery/2c_precomp/RENAMING.md` — terminology migration
  tracker
- `phase_2_discovery/2d_candidate_detection/README.md` — staircase
  detector
- `phase_4_postprocessing/docs/PHASE4_ARCHITECTURE.md` — the phase 4
  design (sections + flow)
- `phase_4_postprocessing/docs/PHASE4B_REWRITE_ARCHITECTURE.md` — the
  four-script 4b subsystem
- `phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`
  — why K=8, what local_Q is, four-way jackknife semantics
- `phase_4_postprocessing/4a_existence_layers/README.md` — catalog
  birth (C01d/C01e/C01g contract)
- `_archive_superseded/fuzzy_merge_abandoned/README.md` — what the
  retired merge was, why it was dropped
