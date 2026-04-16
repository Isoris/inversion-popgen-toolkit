# SESSION_AUDIT_2026-04-17_v3.md

Continuation of the 2026-04-17 session (v2 audit in the uploaded tarball).
Scope: phase 2 folder restructure, catalog-birth contract repair, docs refresh.

---

## Starting state (from tarball)

Upload: `inversion-popgen-toolkit.tar` — the real git repo, not a
deliverable snapshot. At session start:

- Phase 2 layout: `2a_local_pca`, `2b_mds`, `2c_precomp`, `2d_cores`
  (one file — duplicate of C01b_1), empty `2e_blocs/`,
  `2f_candidate_detection`, `2g_ghsl`
- Phase 4 flat restructure already applied (`phase_4_postprocessing/`
  with 4a–4e subfolders)
- Phase 4 tests: 5 pass
- Fuzzy merge `STEP_C01b_2_merge.R` in `_archive/.../snakes/` but no
  clear "retired" vs "kept" status

## What was done

### 1. Layout cleanup

- Removed empty `2e_blocs/` (content moved into 2c_precomp previously
  per user's instruction)
- First renumbered `2f_candidate_detection/` → `2e_candidate_detection/`,
  `2g_ghsl/` → `2f_ghsl/`
- After the fuzzy-merge retirement decision (below), renumbered again:
  `2e_candidate_detection/` → `2d_candidate_detection/` and
  `2f_ghsl/` → `2e_ghsl/`
- Final phase 2: `2a_local_pca`, `2b_mds`, `2c_precomp`,
  `2d_candidate_detection`, `2e_ghsl`, `_archive`

### 2. `2d_cores/` content decision

Investigation found `2d_cores/STEP_C01b_1_cores_wired_registry.R` was
the pre-rename duplicate of `2c_precomp/STEP_C01b_1_seeded_regions.R`
— same 900-line script, same algorithm, same outputs, only
`snake`/`core` vs `region`/`seeded_region` terminology differed.
Deleted. Folder recycled for merge step (see below), then retired
entirely after the user's clarification about merge.

### 3. The fuzzy merge retirement

User clarification on this turn — direct quotes:

> "our merge it was bugged and overmerged so thats why we didnt really
> fix it. then we developed the staircase traversal method and the bloc
> method to see like can we merge where its reliable if it needs to be
> merged. maybe its no need"

> "I thought that we could update our catalogue based on multiple
> boundaries and breakpoints and the OR [...] eventually at the end use
> the cores mds windows seeded and extended to like assign them to their
> candidate"

Then asked why a separate assignment step is needed at all — the answer
turned out to be: **it isn't**. `STEP_C01d_candidate_scoring.R` already
accepts `--cores_dir` and reads seeded regions as a scoring dimension
at catalog-birth. Confirmed by `conversation_search`:

> From chat `fe8ed302`: "C01d is where the catalog is born.
> create_candidate_folders.sh runs here."

> From chat `5b793a68`, the legacy launcher: C01d already accepts
> `--cores_dir`.

Decisions made:

- `STEP_C01b_2_merge.R` (copied earlier this session to
  `2d_seeded_merge/` as `STEP_C01b_2_region_merge.R`) moved to
  `_archive_superseded/fuzzy_merge_abandoned/` with a README
  documenting the overmerge failure modes (LG19 8–15 Mb, LG28 7–15 Mb,
  carriers not matching merged PA)
- `2d_seeded_merge/` folder deleted
- `2d_` slot given to staircase detector (the real primary boundary
  track) — renumbered `2e_candidate_detection` → `2d_candidate_detection`
- GHSL renumbered `2f_ghsl` → `2e_ghsl`

### 4. Critical plug-in bug found in C01d (FIXED)

C01d's `--cores_dir` reader was looking for the pre-rename input:

```r
pattern = "^snake1_core_regions_.*\\.tsv\\.gz$"
# fields: cr$core_family, cr$snake_id, cr$cheat26_status
```

Current `STEP_C01b_1_seeded_regions.R` writes
`seeded_regions_summary_<chr>.tsv.gz` with `scale_tier`, `region_id`,
`test26_status`. Without the fix, C01d would silently find zero cores,
log `"No core_regions files found"`, and incorrectly cap staircase
candidates at Tier 4 via the "no snake core" branch.

Surgical fix applied: new filename pattern first, legacy fallback
second, column normalisation at read time. Output column names
(`d12_snake_concordance`, `snake_overlap`) preserved because
`4e/compute_candidate_status.R` and `test_registry_sanity.py` read them.

### 5. Same class of bug found in C01g (FIXED)

`STEP_C01g_boundary_catalog.R` similarly reads two legacy-named
inputs:

- `snake1_core_regions_<chr>.tsv.gz` → fixed with new+legacy fallback
- `sv_flashlight_<chr>.rds` → fixed with new+legacy fallback

Output `source` field value `"snake_core"` → `"seeded_region"`. No
downstream consumers filter on this value (verified via grep), so safe
to change.

### 6. Top-level `inversion_modules/README.md` rewritten

Old README described `phase_4_catalog/phase4_v10/` + `phase4b_rewrite/`
— a layout that stopped existing when phase 4 was flattened into
`phase_4_postprocessing/4a–4e/`. New README describes current layout,
catalog-birth flow at C01d with `--cores_dir` / `--boundary_dir` /
`--hyp_dir` contract, both phase-2 tracks, the retirement of the
fuzzy merge, and the 5 validation levels.

### 7. Other README touches

- `phase_2_discovery/2c_precomp/README.md` — added PHASE_01C to intro,
  updated layout, redrew workflow DAG to show C01b_1 and PHASE_01C as
  parallel readers of C01a feeding phase 4 directly, added
  "Architectural note — why no 2d_seeded_merge" explaining the
  retirement
- `phase_2_discovery/2d_candidate_detection/README.md` — fixed stale
  "Merges: yes, via C01b_2" row to note the retirement
- `phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R` header
  — "Downstream: STEP_C01b_2" → "Downstream:
  phase_4/4a/STEP_C01d_candidate_scoring --cores_dir"
- `phase_2_discovery/2c_precomp/patches/README.md` —
  `patch_C01b2_merge_flashlight.R` marked NO LONGER APPLICABLE
- `phase_4_postprocessing/4a_existence_layers/README.md` — fixed
  "scripts not yet in this delivery" claim (they are); added the
  "C01d input contract" table with all 7 CLI flags/positional args;
  corrected the catalog-birth flow to say
  `create_candidate_folders.sh` is still owed
- `docs/HOW_TO_RUN.md` → `HOW_TO_RUN_outdated.md` with deprecation
  banner (matches existing `PIPELINE_v8.5_ARCHITECTURE_outdated.md`
  convention)
- `check_deployment_v10.py` — updated allowlist to point at
  `2e_candidate_detection/` then `2d_candidate_detection/` through the
  two rename waves

### 8. `_archive_superseded/fuzzy_merge_abandoned/`

New directory created at the deployment root (parallel to `_archive/`)
with:

- `STEP_C01b_2_region_merge.R` — the retired merge script, with
  terminology migration applied (so future readers see the new names)
- `README.md` — documents the overmerge failure modes, the replacement
  architecture (two parallel boundary detectors + seeded regions as
  C01d scoring dimension), and why the file is kept as design history

## Final state

```
inversion_modules/
├── README.md                           ← rewritten this session
├── CONFIG_ARCHITECTURE.md
├── PATH_WIRING.md
├── check_deployment_v10.py             allowlist updated
├── utils/                              ancestry_bridge + local_q_projector present
├── phase_1_inputs/                     populated
├── phase_2_discovery/
│   ├── 2a_local_pca/
│   ├── 2b_mds/
│   ├── 2c_precomp/                     ← SV prior + precompute + seeded + PHASE_01C
│   │   ├── STEP_C00_build_sv_prior.R
│   │   ├── STEP_C01a_precompute.R
│   │   ├── STEP_C01b_1_seeded_regions.R  (header updated)
│   │   ├── PHASE_01C_block_detect.R
│   │   ├── diags/, patches/
│   │   └── README.md + RENAMING.md
│   ├── 2d_candidate_detection/         ← staircase (primary boundary track)
│   │   └── README.md updated
│   ├── 2e_ghsl/
│   └── _archive/                       legacy v8.5 HPC tree (untouched)
├── phase_3_refine/                     populated
├── phase_4_postprocessing/
│   ├── 4a_existence_layers/            ← catalog birth
│   │   ├── STEP_C01d_*.R               fixed (cores_dir reader)
│   │   ├── STEP_C01e_*.R
│   │   ├── STEP_C01g_*.R               fixed (cores + flashlight readers)
│   │   └── README.md updated
│   ├── 4b_group_proposal/, 4c, 4d, 4e
│   ├── docs/, orchestrator/, patches/, schemas/, specs/, tests/
├── phase_5_followup/, phase_6_secondary/  populated
└── _archive/                           legacy v8.5 HPC tree
_archive_superseded/
└── fuzzy_merge_abandoned/              ← new this session
    ├── STEP_C01b_2_region_merge.R
    └── README.md
```

Deployment check: 5/5 phase 4 test suites pass. 74 pre-existing
warnings, same as session start — no regressions. The 74 break down
roughly as: registry loader path mismatch (registries live at repo
root, not `inversion_modules/` — 6 warnings), archived legacy files
with naive bracket-counter false positives (~60 warnings), a few
pre-existing broken utility files (`theme_systems_plate.R`,
`debug_core_pcas.R`, `STEP_C01e_candidate_figures.R`).

## Key architectural points locked in this session

1. **No merge step in phase 2.** The 1D fuzzy merge overmerged on test
   chromosomes; replaced by two parallel boundary detectors (staircase
   in 2d, PHASE_01C in 2c) + seeded regions as internal evidence.
2. **Seeded regions are internal evidence, not detection output.**
   They feed C01d's D12 scoring dimension via `--cores_dir`, they
   don't define candidates themselves.
3. **Candidate boundaries come from the OR of three sources** at C01g:
   PHASE_01C landscape, staircase, seeded-region edges, plus SV prior
   and blue-cross verdicts (5 sources total in C01g). Union is by
   proximity merging into `boundary_catalog_unified.tsv.gz`.
4. **Catalog is born at C01d.** It reads the staircase-derived
   candidate list, consumes cores / boundary / hypothesis as evidence
   dimensions, writes `candidate_scores.tsv.gz`, assigns tiers.
5. **`create_candidate_folders.sh` is aspirational, not current.**
   Past chats describe it; not in this deployment. Phase 4b–4e read
   the flat catalog directly.

## Still owed (in priority order)

### Blocking a LANTA run

1. Real `Rscript --vanilla -e 'parse(file=...)'` on LANTA for every R
   file — the naive bracket counter is necessary but not sufficient.
   Especially check the allowlisted files (C01b_2 merge, C01f, C01b_1,
   C00, ggplot-heavy detection scripts).
2. One-chromosome smoke test through 2c → 2d → phase_3 → 4a → 4b → 4c.
   Suggested target: LG12 with the Tier-1 candidate at 52.4–54.6 Mb
   from earlier manuscript drafts.
3. `ancestry_bridge.R --prepare` for all 28 chromosomes before the
   first C01a run (produces `local_Q/<chr>.local_Q_samples.tsv.gz`
   cache that `nested_composition` reads).

### Schema / contract

4. Bump `hypothesis_verdict.schema.json` to v2 (add `family_linkage`
   enum + `quality_flags` field) to match v10.1.1 semantics.
5. Wire C01f call sites to the new `compute_group_validation()` list
   return (`$level`, `$quality_flags`, `$family_linkage`) instead of
   the pre-v10.1.1 string return.

### Optional / low priority

6. `create_candidate_folders.sh` if per-candidate folder layout is
   wanted. Not blocking downstream 4b–4e.
7. Cross-module rename pass (post-manuscript) — retire remaining
   `cheatNN`/`snake`/`core` terminology in phase 4 and MODULE_5B–E
   consumers. C01d's output columns `d12_snake_concordance` and
   `snake_overlap` are read by `4e/compute_candidate_status.R` and
   `test_registry_sanity.py`; renaming is coordinated across those
   three files.

## For the next session

If the next chat focuses on **running the pipeline**, start with items
1 and 2 above — parse check and smoke test on LG12.

If the next chat focuses on **phase 4b/4c integration**, items 4 and 5
plus reading through C01f to map the call sites of
`compute_group_validation()` after the v10.1.1 return-type change.

If the next chat focuses on **manuscript**, point at
`phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`
and the top-level README.

Deliverable tarball from this session:
`inversion-popgen-toolkit_phase2_restructured_2026-04-17.tar`
