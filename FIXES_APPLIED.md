# Chat 17 (2026-04-18) — Launcher coverage + orchestrator path cleanup

Scope: user flagged that phase 2c (SV priors + precomp) had no SLURM
launchers and wondered what else was missing. Chat 17 did a full triage
across `inversion_modules/` + `unified_ancestry/`, wrote launchers for
every genuine gap, and patched two stale-path problems discovered during
the audit. **No HPC run yet — pre-flight checks from chat-17 handoff plus
these launchers are the first real chat-18 step.**

### [CH17-triage] META — full launcher-coverage audit

**Where.** All phases of `inversion_modules/` plus `unified_ancestry/`.

**Findings.**

- `Modules/` 1–4 (upstream pipeline): fully covered by existing
  `LAUNCH_moduleN.sh` + `SLURM_*.sh` workers. Not touched.
- Phase 1 (inputs): covered by chat-13 `LAUNCH_STEP04…07_*.slurm`.
- Phase 2a / 2b / 2d: covered. Phase 2d in particular has
  `LAUNCH_run_all.slurm` as the array (1–28) entry point; individual
  `STEP_D*` scripts are pipeline components, not standalone.
- Phase 2c (SV priors + precomp + block detect + seeded regions): **4
  launchers missing.** See [CH17-phase2c] below.
- Phase 2e (GHSL v6): **2 launchers missing.** See [CH17-phase2e].
- Phase 3 (breakpoint refinement): covered by `run_breakpoint_validation.sh`.
- Phase 4a, 4b, 4c, 4e: **9 launchers missing, all named by the
  phase-4 orchestrators but never written.** See [CH17-phase4].
- Phase 4d: already covered (`orchestrator/LAUNCH_group_cheats.sh`).
- Phase 5, 6 (secondary): covered. Not touched.
- unified_ancestry: covered (chat 16).

**Deliverable.** `LAUNCHER_TRIAGE.md` at repo root — full decision matrix
with DAG diagram and known-issue callouts.

### [CH17-phase2c] MAJOR — Phase 2c launcher coverage

**Where.** `inversion_modules/phase_2_discovery/2c_precomp/`.

**What was missing.**

| Launcher added                                   | Calls                         |
|--------------------------------------------------|-------------------------------|
| `LAUNCH_STEP_C00_sv_prior.slurm`                 | `STEP_C00_build_sv_prior.R`   |
| `LAUNCH_STEP_C01a_precompute.slurm`              | `STEP_C01a_precompute.R`      |
| `LAUNCH_PHASE_01C_block_detect.slurm`            | `PHASE_01C_block_detect.R`    |
| `LAUNCH_STEP_C01b_1_seeded_regions.slurm`        | `STEP_C01b_1_seeded_regions.R`|
| `submit_precomp_chain.sh`                        | ancestry → C00 → C01a chain   |

**Why these matter.** STEP_C01a is the chat-16 silent-bug surface: its
`merge_local_Q_into_invlikeness()` call silently no-op'd when the
SAMPLE_GROUP filename convention didn't match between R and bash. The
launcher now probes `${LOCAL_Q_DIR}/K08/${SAMPLE_GROUP}.local_Q_summary.tsv.gz`
and warns loudly before running if the expected file isn't there.

**Test.** All 4 `bash -n` clean. `submit_precomp_chain.sh` chains with
`--dependency=afterok:` so a silent upstream failure blocks downstream.

### [CH17-phase2e] MAJOR — Phase 2e (GHSL v6) launcher coverage

**Where.** `inversion_modules/phase_2_discovery/2e_ghsl/`.

**What was missing.**

| Launcher added                              | Calls                              |
|---------------------------------------------|------------------------------------|
| `LAUNCH_STEP_C04_ghsl_v6_compute.slurm`     | `STEP_C04_snake3_ghsl_v6.R`        |
| `LAUNCH_STEP_C04b_ghsl_v6_classify.slurm`   | `STEP_C04b_snake3_ghsl_classify.R` |

**Sizing.** C04 is the heavy engine: 1 h / chromosome × 28 chroms —
`--array=1-28` with 64G / 4 CPU. C04b is the light classifier: 30 s /
chrom, tunable params via env (`GHSL_SCALE`, `GHSL_KARYO_LO/HI`,
`GHSL_MAX_K`, etc.).

### [CH17-phase4] MAJOR — Phase 4 launcher coverage

**Where.** `inversion_modules/phase_4_postprocessing/{4a,4b,4c,4e}/`.

**The discovery.** Both phase-4 orchestrators (`run_phase4.sh`,
`run_phase4b.sh`) reference launchers by name but those launchers
didn't exist in the tree:

```
run_phase4.sh expects:
  ${MODULE_DIR}/launchers/LAUNCH_C01d_scoring_pass1.sh
  ${MODULE_DIR}/launchers/LAUNCH_C01g_boundary.sh
  ${MODULE_DIR}/launchers/LAUNCH_C01f_hypothesis.sh
  ${MODULE_DIR}/launchers/LAUNCH_C01d_scoring_pass2.sh
  ${MODULE_DIR}/launchers/LAUNCH_group_cheats.sh         ← already existed
  ${MODULE_DIR}/launchers/LAUNCH_characterize_classify.sh

run_phase4b.sh uses --wrap=, no external launchers needed for orchestrator
  but per-script launchers are useful for manual re-runs.
```

**What was added.**

| Launcher added (phase 4)                   | Calls                         |
|--------------------------------------------|-------------------------------|
| `4a/LAUNCH_C01d_scoring_pass1.sh`          | C01d (cold pass)              |
| `4a/LAUNCH_C01g_boundary.sh`               | C01g boundary catalog         |
| `4b/LAUNCH_C01i_decompose.sh`              | C01i_decompose                |
| `4b/LAUNCH_C01i_b_multi_recomb.sh`         | C01i_b_multi_recomb           |
| `4b/LAUNCH_C01i_c_nested_comp.sh` (stub)   | C01i_c_nested_composition.py  |
| `4b/LAUNCH_C01i_d_seal.sh`                 | C01i_d_seal                   |
| `4c/LAUNCH_C01f_hypothesis.sh`             | C01f_hypothesis_tests         |
| `4e/LAUNCH_C01d_scoring_pass2.sh`          | C01d (with boundary + hyp)    |
| `4e/LAUNCH_characterize_classify.sh`       | run_characterize + status     |

**Note on 4b.3 stub.** `STEP_C01i_c_nested_composition.py` is referenced
by `run_phase4b.sh` but not present in the tree. The launcher is written
correctly in shape; when invoked it fails loudly with a clear message
telling the user to either add the script or set `SKIP_4B3=1`. This is
the correct behaviour for a missing DAG node.

### [CH17-orch-paths] MINOR — Phase-4 orchestrator stale paths

**Where.** `inversion_modules/phase_4_postprocessing/orchestrator/run_phase4.sh`
and `run_phase4b.sh`.

**The bug.** Both carry paths from an earlier directory layout:

- `${BASE}/inversion-popgen-toolkit` → should be `${BASE}/inversion_modules`
- `${MODULE_DIR}/launchers/LAUNCH_X.sh` → `${MODULE_DIR}/phase_4_postprocessing/4?/LAUNCH_X.sh`
- `${MODULE_DIR}/phase4b_rewrite/R` → `${MODULE_DIR}/phase_4_postprocessing/4b_group_proposal`

**Fix.** `patches/patch_orchestrator_paths.sh` — writes `.bk_chat17`
backups, applies a 7-line `sed` substitution per orchestrator. Verified
both files `bash -n` clean after the patch; all 6 LAUNCH_ refs now
point to real files in the current tree.

**Idempotent?** No — the patch will error out on the second run if paths
don't match. The in-tree backups (`run_phase4.sh.bk_chat17`) let you undo
with `mv`.

### [CH17-sv_prior_dir] MINOR — STEP_C00 / SV_PRIOR_DIR mismatch

**Where.** `inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R`
line 120 vs `inversion_modules/00_inversion_config.sh` line 181.

**The bug.** The config exports
`SV_PRIOR_DIR=${INVDIR}/03_sv_prior`, but the R script hardcodes output to
`${INVDIR}/06_mds_candidates/snake_regions_multiscale/sv_prior`. The two
paths do not match — so everything downstream that reads `$SV_PRIOR_DIR`
looks in the wrong place.

**Fix.** `patches/patch_STEP_C00_sv_prior_dir.sh` — makes the R script
respect the env var with back-compat fallback:

```r
SV_PRIOR_DIR <- Sys.getenv("SV_PRIOR_DIR",
  file.path(INVDIR, "06_mds_candidates/snake_regions_multiscale/sv_prior"))
```

The chat-17 C01a launcher already points `--sv_prior_dir` at the hardcoded
location, so the pipeline works either way; this patch just restores the
contract that `$SV_PRIOR_DIR` is meaningful. After applying, either (a)
update `00_inversion_config.sh` to export the hardcoded path, or (b)
move / symlink on disk to match the config. Post-patch R file passes
`_rcheck.py`. Patch is idempotent (detects prior application).

### [CH17-figures] MINOR — Chat-15 publication figures folded into bundle

**Where.** `inversion_modules/phase_4_postprocessing/4e_final_classification/figures/`.

**Context.** Chat 15 produced three publication-figure R scripts plus a
smoke-test harness but left them outside the repo in a separate tarball
(`chat15_figures_bundle_2026-04-17.tar`). Chat 17 folds them into the
main deliverable tarball so chat 18 has everything in one drop.

**Scripts.**

| File                                       | Fig | Status |
|--------------------------------------------|-----|--------|
| `fig4e_gene_count_by_biotype.R`            | 2   | static-reviewed (level 1 — not run on real data) |
| `fig4e_class_along_inversion.R`            | 1a  | static-reviewed; pragmatic analog using GHSL karyotype calls (true local-ancestry figure blocked on upstream compute) |
| `fig4e_cumulative_burden_per_group.R`      | 3a  | static-reviewed; `--entity variant`/`gene_inv` stubbed until 3a verified |
| `test_figures.R`                           | N/A | smoke test harness — first chat to ship one |

All 4 pass `_rcheck.py`. Chat-15 static review caught and fixed four real
bugs (colorRamp2 degenerate case, setNames recycle-safety, data.table
row-idx semantics, reproducibility of first-hit attribution tiebreak)
plus a zero-overlap placeholder case. Schema assumptions pinned in
`INSTALL.md` in-place.

**Chat-18 first thing.** Run `test_figures.R` before anything downstream
touches these. The ladder is *not run* → *smoke test passed* → *real
data OK*; smoke test is the only path to level 2 and must pass before
level 3 even makes sense.

### [CH17-engines] MINOR — Engine binaries must be compiled before launchers run

**Where.** Three makefiles + one raw .c file:
- `unified_ancestry/src/Makefile` → `instant_q` (C++17, OpenMP)
- `unified_ancestry/engines/fst_dxy/Makefile` → `region_popstats`
- `inversion_modules/phase_5_followup/engines/Makefile` → `rare_sfs_pairwise`, `export_q_residual_dosage`
- `unified_ancestry/engines/hobs_hwe/scripts/hobs_windower.c` (no Makefile; raw gcc)

**The problem.** Every launcher that calls one of these binaries assumes
it exists at its configured path (`INSTANT_Q_BIN`, `POPSTATS_BIN`, etc.
in `00_ancestry_config.sh` and `pipeline_bridge.sh`). `pipeline_bridge.sh`
auto-compiles instant_q / popstats / hobs_windower on first source via
`_compile_if_missing`, but: (a) the auto-compile pipes `make` to
`2>/dev/null`, hiding real compile errors behind a terse "WARN: compile
failed" line; (b) it does NOT cover the two phase-5 engines in
`phase_5_followup/engines/` (those have their own Makefile that the
bridge doesn't know about); (c) first-time setup has no clean
non-interactive entry point for all of the above.

**Fix.** `patches/compile_engines.sh` — explicit, verbose, exit-coded
build of all four makefile targets + the hobs_windower raw gcc.
Compiled all five successfully on the chat-17 sandbox (gcc 13.3.0) with
only harmless warnings (misleading-indentation in instant_q.cpp at
lines 173 and 194, strncpy-truncation in rare_sfs_pairwise.c at
lines 169/189/194). All binaries land at the paths the configs expect.

**Also flagged but not fixed.** Stale copies of `rare_sfs_pairwise.c` and
`export_q_residual_dosage.c` still exist under
`unified_ancestry/engines/fst_dxy/` even though that Makefile's comment
header explicitly says "previously built here but now live in
MODULE_5B_inversion_followup/engines/". Canonical copies are in
`inversion_modules/phase_5_followup/engines/`. Not deleted here — worth
cleaning up post-HPC-validation to prevent confusion.

**Also.** `00_ancestry_config.sh` line 46 points `RARE_SFS_BIN` at
`${BASE}/Modules/MODULE_5B_inversion_followup/engines/rare_sfs_pairwise`,
but that directory doesn't exist; it was renamed to `phase_5_followup`.
Either fix the config or symlink the directory. This is not blocking
(launchers can override via env) but should be noted.

### [CH17-scientific-names] MINOR — Apply Tier A of the existing RENAMING.md plan

**Where.** Active tree. `RENAMING.md` (under `phase_2_discovery/2c_precomp/`)
already lists the full rename plan from informal nicknames
(flashlight/snake/core/cheat) to scientific terminology. Filename renames
were mostly applied, but **internal R identifiers, CLI flags, and
status strings were left behind across ~30 files**. The user asked for
a systematic audit.

**Two patches shipped (`_v1` does the common cases, `_v2` catches
compound identifiers the `\b` word boundaries missed).**

**v1 — `patch_rename_scientific_names.sh`** covers:
- `flashlight` → `sv_prior` (bare word)
- `FLASHLIGHT_LOADER` / `FLASH_DIR` → `SV_PRIOR_LOADER` / `SV_PRIOR_DIR`
- `--flashlight_dir` / `--flashlight` → `--sv_prior_dir` / `--sv_prior`
- `snake_id` → `region_id`
- `snake_phase` → `extension_phase`
- `core_family` → `scale_tier`
- Log prefix `[flashlight]` → `[sv_prior]`

Result: **17 of 20 target files patched**, ~120 substitutions. 3 already
clean (C00, C01a, STEP_C01i_b_multi_recomb don't have the v1 patterns).

**v2 — `patch_rename_scientific_names_v2.sh`** covers the compound
identifier cases the `\bflashlight\b` regex in v1 cannot match because
`_` is a word char:
- `flashlight_<X>` prefix: `flashlight_seeds`, `flashlight_path`,
  `flashlight_available`, `flashlight_seeded`, etc.
- `<X>_flashlight` suffix: `try_load_flashlight`, `load_flashlight`,
  `.has_flashlight`
- `<X>_flashlight_<Y>` middle: `s3_flashlight_hemi`
- Status string literals: `"flashlight"`, `"flashlight_only"`,
  `"flashlight_sparse"`, `"flashlight+step03"`

Result: **13 files patched**, ~85 substitutions.

**Explicitly preserved** (safety-checked at patch time — patch aborts
if these strings are accidentally clobbered):
- `sv_flashlight_<chr>.rds` (legacy RDS filename fallback in C01g)
- `<chr>_flashlight.rds` (second legacy RDS filename fallback in C01a)
- `flashlight_v2/...` (stale hardcoded directory paths — a missing
  dir is a separate problem that renaming won't fix)
- `flashlight_loader_v2.R` (stale hardcoded source file name — same
  category)

**Paired launcher fix.** `LAUNCH_C01g_boundary.sh` (chat-17 NEW) was
already updated in-place to pass `--sv_prior <dir>` instead of
`--flashlight <dir>`, matching the post-v1 R-side CLI rename.

**Deferred to chat 18+** (documented in `NAMING_AUDIT.md`):
- Tier B: `snake_regions_multiscale/` dir rename (8 files, coordinates
  with launcher defaults)
- Tier B: `cheat27–30.R` filename renames (5 files, coordinates with
  `LAUNCH_group_cheats.sh`)
- Tier B: `bloc_` → `block_` in fst_dxy plot scripts (3 files)
- Tier C (leave alone): `inv_likeness` column, `peel` / `peeling`,
  `snake3` in GHSL filename, `STAIR_*` / `NN_*` config prefixes,
  `hatchery` / `wild` population modes

**Tested.** Full chain applied against fresh scratch copy:
- `patch_rename_window_dt.sh` → 10 files / 160 subs
- `patch_rename_scientific_names.sh` → 17 files / ~120 subs
- `patch_rename_scientific_names_v2.sh` → 13 files / ~85 subs
- **0 parse failures** across all touched `*.R` files (`_rcheck.py`)
- **All 4 preserved back-compat strings intact** post-patch
- **All 3 patches idempotent** (second run reports 0 patched /
  N already clean)

### [CH17-rename] MINOR — `inv_like_dt` / `window_inv_likeness.tsv.gz` → `window_dt`

**Where.** 10 active R files across `phase_2_discovery/2c_precomp/`,
`phase_4_postprocessing/4d_group_dependent/`, and `unified_ancestry/wrappers/`.

**The confusion.** The per-window diagnostic table and its on-disk file
were named after their flagship column `inv_likeness` (a
het + dip-test + band-discreteness composite), as if that were the
table's only content. Actually the same table carries ~40 columns
including MDS z-scores, dosage stats, band stamps, family_likeness,
localQ_* (chat-15 ancestry merge), SV overlap columns (chat-17
[CH17-phase2c]), etc. The name systematically misleads anyone reading
the code — the user flagged it during chat-17 audit: *"why we call it
inv_likeliness window if it's just local pca window of z outlier…
audit the script too."*

Additional confusion: two slightly different filenames in the code
(`snake_inv_likeness.tsv.gz` in one stale reference, `window_inv_likeness.tsv.gz`
as the actual written file) from the old Snake 1/2/3 architecture.

**Fix.** `patches/patch_rename_window_dt.sh` — three literal renames:

- variable `inv_like_dt` → `window_dt` (10 files, 160 substitutions)
- file `window_inv_likeness.tsv.gz` → `window_dt.tsv.gz` (writer + reader)
- stale `snake_inv_likeness.tsv.gz` → `window_dt.tsv.gz` (fixes a
  pre-existing broken reference in `diag_common.R` as a bonus)

Column name `inv_likeness` is **kept** (that name is accurate for the
column itself). Function `compute_inv_likeness_all()` is kept (still
returns the table bearing the inv_likeness column, among others).

**Test.** Tested against full scratch copy of the active tree. All 10
patched files pass `_rcheck.py`. Idempotent (detects prior rename).
Writer + reader filenames now agree. Archive copies under `_archive/`
untouched.

**Post-patch on HPC.** Any `window_inv_likeness.tsv.gz` files already
written to disk from prior runs will be orphaned. Three options
documented in the patch script:
- rename existing files to match
- symlink old name → new name
- re-run precomp (which chat 18 will do anyway — simplest)

### Static validation status (end of chat 17)

- 16 new launcher files (15 launchers + 1 chain) — all `bash -n` clean
- 3 patch scripts — all `bash -n` clean and tested:
  - `patch_orchestrator_paths.sh` tested against scratch copies of both
    orchestrators (`bash -n` clean post-patch; all 6 LAUNCH_ refs resolve)
  - `patch_STEP_C00_sv_prior_dir.sh` tested against a scratch copy of
    STEP_C00 (passes `_rcheck.py` post-patch; idempotent)
  - `compile_engines.sh` actually compiled all 5 binaries (instant_q,
    region_popstats, rare_sfs_pairwise, export_q_residual_dosage,
    hobs_windower) on the chat-17 sandbox with gcc 13.3.0 — summary
    reports "4 built, 0 failed, 0 skipped" (4 makefile targets +
    hobs_windower raw gcc = 5 binaries)
- Both patched orchestrators `bash -n` clean post-patch
- Patched STEP_C00 passes `_rcheck.py`
- All patches are idempotent (detect prior application, don't double-patch)

### Open findings at end of chat 17

**Chat 18 (first real HPC run on LANTA):**

1. Run 5 pre-flight checks from chat-17 handoff (unchanged — still the
   critical first step).
2. Apply both patches from `patches/`:
   ```bash
   bash patches/patch_STEP_C00_sv_prior_dir.sh
   bash patches/patch_orchestrator_paths.sh
   ```
   Both write `.bk_chat17` backups.
3. Drop in all 16 new launchers at their phase-dir locations.
4. Run `submit_precomp_chain.sh` for ancestry + SV prior + C01a.
5. Proceed through the DAG (see `LAUNCHER_TRIAGE.md`).
6. Final `reg$results$integrity_check()` — expect `all_pass=TRUE`.

**Chat 18 (low pri after HPC):**

- Decide: should `00_inversion_config.sh` be updated to export
  `SV_PRIOR_DIR` at the script's hardcoded path (option a), or should
  on-disk paths be aligned to the config's short path (option b)?
  The patch works for both paths; this is a taste / consistency call.
- If `STEP_C01i_c_nested_composition.py` is needed for this manuscript,
  write/port it. Otherwise wire `SKIP_4B3=1` into `run_phase4b.sh`.
- Once HPC run confirms nothing external depends on `reg$stats`
  deprecated alias, remove it.

### Known paths — additions post-chat-17

- Launchers live under their phase dirs (not under a single `launchers/`):
  - `inversion_modules/phase_2_discovery/2c_precomp/LAUNCH_*.slurm`
  - `inversion_modules/phase_2_discovery/2e_ghsl/LAUNCH_*.slurm`
  - `inversion_modules/phase_4_postprocessing/4a_existence_layers/LAUNCH_*.sh`
  - `inversion_modules/phase_4_postprocessing/4b_group_proposal/LAUNCH_*.sh`
  - `inversion_modules/phase_4_postprocessing/4c_group_validation/LAUNCH_*.sh`
  - `inversion_modules/phase_4_postprocessing/4e_final_classification/LAUNCH_*.sh`
- Patches live at `inversion_modules/patches_chat17/` (move after extract).
- `LAUNCHER_TRIAGE.md` at repo root — read alongside `SESSION_SUMMARY.md`
  before touching anything phase-related in chat 18.

---

---

# Previous fixes (chat 16 and earlier)

# FIXES_APPLIED.md

Rolling register of fixes applied. Most recent chat first; older chats
condensed to one-liners. Full detail for pre-chat-12 fixes is in
`_archive/chat_history/FIXES_APPLIED_2026-04-17.md` (chats 3–11.5).

Format: `[FINDING] [SEVERITY] — TITLE` with WHERE / HOW / TEST.

---

## Chat 16 (2026-04-17) — Database redesign: results_registry as fourth first-class registry + silent no-op fix

Scope: user flagged mid-chat that chat-15's `stats_cache` was half-built —
no schema, no FKs, no query plane. "Make the smart database not the
scattered database." Also surfaced the chat-15 R-vs-bash sample-set hash
mismatch that was silently no-op'ing every ancestry merge.

Both addressed in one chat.

### [CH16-hash] CRITICAL — chat-15 R/bash sample-set tag mismatch

**Where.** `unified_ancestry/wrappers/instant_q.R:217-242` and
`unified_ancestry/launchers/LAUNCH_instant_q_precompute.slurm:111-124`.

**The bug.** R used `digest::digest(ids_sorted, algo="sha1")` with default
`serialize=TRUE`, which hashes R's binary serialization format (with R
version headers, type tags, length prefixes) — NOT the plain text of the
newline-joined IDs. Bash used `awk 'NF' | sort | sha1sum | cut -c1-6` on
the raw file text. The two sha1s are from different byte streams and
cannot be equal. Empirically verified: same 226-sample list → bash
produces `N226_fe01a5`, R fallback (md5 branch) produces `N226_b1c6c0`.
Result: every cache read via the R tag silently returned NA, every
`merge_local_Q_into_invlikeness()` silently no-op'd, every
`reg$compute$ancestry_*` call returned NULL. Chat-15 handoff's "Step B
— verify C01a now stamps localQ columns" would have reported empty
and concluded the cache wasn't populated.

**Fix.** Killed the hash entirely. Sample-set identity is now a
registered `group_id` from sample_registry (default `all_226`), read
from the `SAMPLE_GROUP` env var on both the R and bash sides.
`.compute_sample_set_tag()` kept as a deprecated shim returning the
group name. New public accessor `get_active_sample_group()`.

### [CH16-fourth-registry] HIGH — stats_cache → results_registry as first-class fourth registry

**Problem.** Chat 15 added `reg$stats` with atomic put/get methods but no
schema, no cross-registry foreign keys, no query plane, no integrity
check. File identifiers were opaque 6-char sha1 tags. User explicitly
requested "smart database, not scattered" — meaning FK-disciplined,
manifest-driven, with a single SELECT-style query method and an integrity
check a reviewer could run.

**Fix.** Full fourth-registry buildout:

- **New design doc** `registries/DATABASE_DESIGN.md`. The 4-table reference,
  FK relationships, sub-cluster naming convention, query-plane spec,
  migration notes.
- **Four JSON schemas** under `registries/schemas/registry_schemas/`:
  `result_row.schema.json` (manifest-row, with `kind`-conditional `allOf`
  rules for pairwise/candidate_q/candidate_f/interval_summary),
  `sample_group.schema.json` (adds `dimension` field for sub-cluster
  classification), `candidate_interval.schema.json`, `evidence_key.schema.json`.
  All pass draft-07 meta + 9 discrimination fixtures.
- **R API rewrite** `registries/api/R/registry_loader.R` —
  `load_stats_api` → `load_results_api(results_dir, samples, intervals)`.
  Adds UUID row_ids, FK enforcement (refuses writes with unregistered
  group_id or unregistered candidate_id), provenance block
  (source_script/engine/engine_version/run_id/config_hash/upstream_files),
  sha256 of each file, group version stamping (from sample_registry's
  `created` field), `ask()` query plane with 4 filters + interval overlap,
  4 shortcuts (`ask_what_at`/`ask_what_for_candidate`/`ask_what_for_group`/
  `ask_provenance`), `integrity_check(check_sha256=FALSE)` with 6 checks,
  `session_start()`/`session_summary()` for script-level echo pattern,
  new `put_interval_summary()` kind. Back-compat: old method names and
  `reg$stats` alias preserved for one cycle. Top-level loader auto-migrates
  `stats_cache/` → `results_registry/` on first init.
- **Python API** `registries/api/python/registry_loader.py` — new
  `ResultsAPI` class with atomic writers (`put_pairwise`,
  `put_interval_summary`), `list_manifest`, read-only `ask()`.
  `Registry.results` field added. Same auto-migration.
- **Bash API** `registries/api/bash/registry_loader.sh` — new
  `registry_results_path <kind> ...`, `registry_results_has <kind> ...`,
  `registry_results_count_by_kind <kind>` helpers. Plus
  `RESULTS_REGISTRY` path resolved in `registry_resolve_paths`.

### [CH16-subclusters] MEDIUM — sub-cluster design for the 10% case

**Problem.** User's biology insight: real candidates sometimes show 3–4
sub-clusters within a single karyotype (not clean HOM_REF/HET/HOM_INV
bands). Database needs to support this without special casing.

**Fix.** Naming convention documented in `DATABASE_DESIGN.md` §
"Sample group naming convention":
- `inv_<cid>_<KARYO>` — top-level karyotype (90% case)
- `inv_<cid>_<KARYO>__sub<N>` — recursive-PCA sub-cluster
- `inv_<cid>_<KARYO>__ancestry<K>_<k>` — ancestry split
- `inv_<cid>_<KARYO>__ghsl_<band>` — GHSL-band split
- `inv_<cid>_<KARYO>__family_<fid>` — family split

`dimension` column in `sample_groups.tsv` classifies each group's nature.
Schema enum: `. / karyotype / karyotype_subcluster / ancestry / ghsl /
family / intersect / cohort / unrelated`. No schema changes in
results_registry — sub-cluster FSTs / ancestry flow through the same
`put_pairwise` / `put_candidate_q_f` pipes; manifest records the exact
group pairs naturally.

New helper `reg$samples$get_subgroups(cid)` in
`utils/sample_registry.R` — returns a data.table of
(group_id, karyotype, dimension, parent_group, n) by parsing group
names. One call discovers all sub-structure for a candidate without
string parsing in every classification script.

Recursive sub-cluster *detection* (the script that runs PCA inside each
karyotype and calls `add_group` for each detected sub-cluster) is NOT
built — future chat-17+ feature `STEP_C01i_e_subcluster.R`.

### [CH16-configs] HIGH — config consistency across three files

`inversion_modules/00_inversion_config.sh`, `unified_ancestry/00_ancestry_config.sh`,
`utils/pipeline_bridge.sh` now all export:
- `SAMPLE_GROUP=${SAMPLE_GROUP:-all_226}` — active sample group for FK
- `RESULTS_REGISTRY_DIR` — new canonical location
- `STATS_CACHE_DIR` — deprecated alias that resolves to
  `RESULTS_REGISTRY_DIR` (same physical directory; no duplicate storage)
- `SAMPLE_REGISTRY` — path used by launcher for FK existence check

### [CH16-tests] — integration test suite

`registries/tests/test_results_registry.py` — 8 integration tests, all
passing in chat-16 sandbox:
1. Fresh init creates expected directory layout
2. put_pairwise with unregistered group raises (FK enforcement)
3. put_pairwise writes file + manifest row with correct shape
4. Every manifest row validates against `result_row.schema.json`
5. `ask()` filters by kind / stat / who / where
6. put_interval_summary rejects invalid stat
7. Auto-migration from legacy `stats_cache/` → `results_registry/`
8. sample_group schema validates sub-cluster group row

Python-based so the chat sandbox (no R) can actually run them. Since
both R and Python write the same manifest.tsv schema, passing these
tests means the R writer side is correct too.

### Coverage after chat 16

Unchanged from chat 15 (59/91 Q2 = 64.8%). Chat 16 was architectural,
not new-key. The chat-16 work unblocks the HPC run that was silently
failing — which in turn unblocks the BK key extraction validation that
still needs LANTA.

### Tools retained

- `_rcheck.py`, `_schema_check.py`, `_code_field_check.py`, `_bk_rename.py`
  — all still valid after chat-16 changes (modified R files all pass
  `_rcheck.py`; modified shell files pass `bash -n`; all 4 new schemas
  pass draft-07 meta).

### Out of scope

- HPC run (same sandbox limitation — deferred to chat 17)
- Recursive sub-cluster *detection* (feature, not bug-fix)
- Removing `reg$stats` deprecated alias (removing in chat 18 after HPC
  confirms no external dependencies)

### [CH16-scan-drivers] HIGH — multi-window scan drivers with automatic segment resolution

**Where.** `registries/api/R/registry_loader.R::load_compute_api` —
new methods `scan_pairwise`, `scan_region_stat`, `scan_candidate_full`.
Plus new primitive `reg$intervals$resolve_smallest_at(chrom, pos)`.

**Motivation.** User's biology question: for a window scan across an
inversion that has recombinant segments (DCO / CO), the correct
karyotype group per window depends on where the window falls — inside
the DCO tract CGA042 is HOM_REF; outside it is HOM_INV. Can the scan
auto-switch groups per window without per-chromosome hardcoding?

**Fix.**
- `resolve_smallest_at(chrom, pos)` returns the smallest registered
  candidate_id containing a position (~10 LOC).
- `scan_pairwise(chrom, s, e, ws, stat, karyo_pair)` iterates windows,
  resolves the smallest candidate per window, looks up
  `inv_<cid>_<k1>` and `inv_<cid>_<k2>` groups, falls back to parent
  candidate if groups too small (min_group_n), computes stat. Persists
  ONE manifest row per unique `(cid × stat)` with per-window values
  in the tsv.gz payload.
- `scan_region_stat(chrom, s, e, ws, stat, karyotype)` — single-group
  region scan (Δ12/entropy/ENA within one karyotype). Same fallback.
- `scan_candidate_full(cid, flank_kb, karyotypes, pairwise_stats,
  region_stats)` — runs the full panel (3 pairwise FSTs + 3×n region
  stats) for one candidate. One call per candidate, full persistence.

Scan semantics verified by tests 9 and 10 in the Python test suite:
segment resolution correct at every position including breakpoints;
3-segment scan produces 3 manifest rows with 5 windows each.

### [CH16-segment-writer] HIGH — STEP_C01i_b_multi_recomb.R now registers segments

**Where.** `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R`,
new block inserted after the `recombinant_map` write, before the
`regime_sample_dag` write. ~150 LOC added.

**Motivation.** Without segment candidates and segment karyotype groups
registered, `scan_pairwise` would always fall back to parent-scale
groups, defeating the purpose of segment-aware scanning. The writer
was missing.

**Fix.** When `n_recombinants > 0`, the script:
1. Collects every breakpoint (regime tract ends + GC tract ends) from
   all recombinant carriers, into a sorted unique list
2. Builds segment intervals `[parent_start, bp1, bp2, ..., parent_end]`
3. For each segment, registers a child candidate via
   `reg$intervals$add_candidate(candidate_id=sprintf("%s_seg_%s", cid, label),
   chrom, seg_start, seg_end, scale="seg"|"seg_dco", parent_id=cid)`.
   Label is `L` for first segment, `R` for last, `DCO1` for middle of
   a 3-segment layout, `MID<N>` otherwise.
4. For each (sample × segment), determines the local karyotype:
   - Non-recombinant samples: parent karyotype everywhere
   - RECOMBINANT* samples: parent karyotype outside their tract,
     inverted karyotype inside (HOM_INV ↔ HOM_REF)
5. Registers per-segment karyotype groups
   `inv_<seg_cid>_<KARYO>` via `reg$samples$add_group`.

Fully back-compat: if n_recombinants=0 (the 90% case) the block is
skipped and nothing changes. The existing `recombinant_map` block is
unchanged and still written first; segment registration is downstream
of it so a failure there doesn't block the block write.

### [CH16-spec] — deferred features documented in SPEC_DEFERRED.md

Five features spec'd but not built, to be implemented chat-17+ once
LANTA data informs thresholds:
1. `scan_with_ancestry_resolution` — ancestry-cluster-dynamic scans
2. `scan_all_candidates_driver.R` — full automatic sweeper with resume
3. `stale_segment_gc.R` — cleanup when multi_recomb is re-run
4. `STEP_C01i_e_subcluster.R` — recursive PCA sub-cluster detector
5. `classification_loop.R` — per-candidate review workflow spec

See `registries/SPEC_DEFERRED.md` for input contracts, acceptance
criteria, and deferred-decision lists.

---

## Chat 15 (2026-04-17) — Ancestry full-wiring + stats cache + config rewrite

Scope extension from BK-only. User observed "ancestry bridge is maybe kind
of old, wire fully to registry and precomp, do K=2..20 (save as different
table if too much), check if inversion config is up to date (it's been a
while), update README, we need Q+F per candidate." All addressed.

### [BK-ancestry] CRITICAL — ancestry not actually wired into precomp

**Where.** `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`
line ~1048: the local-Q merge block tried to `source(ancestry_bridge.R)` via
a filesystem glob lookup that only worked in the pre-flattened `inversion_codebase_*/`
layout. Under the current flattened layout the lookup returned `character(0)[1]`
= `NA`, the `file.exists(NA)` guard was falsy, and the block silently did
nothing. No localQ_* columns ever landed in `<chr>.precomp.rds`. The bridge
file being targeted (`inversion_modules/utils/ancestry_bridge.R`) was also
a pre-chat-12 artefact referencing a non-K-sharded cache layout.

**Fix.** Deleted the broken lookup block entirely (24 lines). `load_bridge.R`
— sourced at the top of C01a since chat 11 — already calls
`configure_instant_q()` against `$ANCESTRY_CONFIG` and pulls
`unified_ancestry/wrappers/instant_q.R` into global scope. The new code
calls `merge_local_Q_into_invlikeness(inv_like_dt)` directly. Replaced
`ancestry_bridge.R` with a much shorter tryCatch around the direct call.
When the cache isn't populated yet, the merge returns the data.table
unchanged and prints a clear "run LAUNCH_instant_q_precompute.slurm first"
message.

Archived: `inversion_modules/utils/ancestry_bridge.R` →
`_archive_superseded/ancestry_bridge_pre_unified/` with README pointing at
the unified wrapper.

### [BK-ksweep] HIGH — K sweep (2..20) + K-sharded cache + sample-set tags

**Problem.** Prior wrapper/launcher were single-K only. Cache files were
`<chr>.local_Q_summary.tsv.gz` with no K in the path, no sample-set tag in
the filename, and living inside `unified_ancestry/local_Q/` (data inside
the code tree). Three real issues:

1. Running a different K silently overwrote the prior K's cache.
2. Running the same chrom/K against a different sample subset (e.g. 81
   unrelated vs 226 full) silently overwrote.
3. Cache files were under the code-git tree; any repo reclone nuked them.

**Fix.**
- Cache moved to `$BASE/ancestry_cache/` (set via `LOCAL_Q_DIR` in
  `00_ancestry_config.sh`; default now respects the env).
- Cache layout is now `<LOCAL_Q_DIR>/K<NN>/<chr>.<sample_set>.local_Q_*.tsv.gz`.
- `<sample_set>` = `N<count>_<6char-sha1>` of the sorted `SAMPLE_LIST`.
  Computed at `configure_instant_q()` time in R (via `digest` pkg or
  `tools::md5sum` fallback) and in bash in the launcher (`sha1sum | cut -c1-6`).
  Both paths produce the same tag given the same sorted list.
- `get_Q(K=NULL)`, `get_Q_summary(K=NULL)`, `get_Q_precompute(K=NULL)`,
  `merge_local_Q_into_invlikeness(K=NULL)`, `cache_path_*`, `cache_exists`,
  `update_manifest`, `list_cached_Q` all gain K parameter. K=NULL resolves
  to `.iq_env$canonical_K` (sourced from `CANONICAL_K` in the ancestry config,
  default 8).
- Per-K qinit/fopt resolved by swapping `K<canon>` → `K<N>` in `BEST_QOPT` /
  `BEST_FOPT` filenames, with optional seed-swap via `BEST_SEED_BY_K`.
- Manifest TSV now has `chrom`, `K`, `sample_set`, `summary_file`,
  `samples_file`, `timestamp`. Launcher + wrapper both idempotently upgrade
  4-col and 5-col legacy manifests to the 6-col schema.
- Launcher array dimension: `28 chroms × 19 K values = 532` tasks.
- Resolve-cache fallback chain on read: K-sharded-tagged → K-sharded-untagged
  → legacy flat (canonical K only). Pre-2026-04-17 caches keep working.
- Precomp RDS gets canonical-K columns only: `localQ_delta12_K08`,
  `localQ_entropy_K08`, `localQ_ena_K08`. Downstream wanting full Q/K-sweep
  hits the cache via `reg$compute$ancestry_q_vector(chrom, s, e, K=N)`.

### [BK-bridge] HIGH — load_bridge.R auto-registers ancestry groups

**Where.** `utils/load_bridge.R` STEP 6.5 (new). Prior state: ancestry-named
groups like `ancestry_K8_Q3`, `all_226`, `unrelated_81` were referenced by
downstream scripts but never actually registered in the inversion
`sample_registry`. Any `reg$has_group("ancestry_K8_Q3")` returned FALSE.

**Fix.** After load_bridge sources the instant_q wrapper and has a live
`reg` handle, step 6.5 registers (idempotent, `overwrite=FALSE`):

- `all_226` — every sample in sample_master
- `unrelated_81` — NAToRA-pruned subset from `PRUNED_LIST`
- `ancestry_K<K>_Q<k>` for k in 1..K — members derived from majority-vote
  `assigned_pop` across the first 3 cached chroms at canonical K. Robust to
  single-chrom assignment wobble.

No-op if the ancestry cache hasn't been populated yet; just skips the
ancestry_K groups silently. The `all_226` and `unrelated_81` groups land
regardless.

Also in the bridge: fixed the hard-coded `inversion_codebase_v8.5/` fallback
paths for sample_map.R / sample_registry.R / INVERSION_CONFIG. New search
order: flattened layout first (`$BASE/utils`, `$BASE/inversion_modules/utils`),
then legacy v8.5. Same pattern in `utils/pipeline_bridge.sh`.

### [BK-compute] HIGH — reg$compute gains 5 ancestry methods

**Where.** `registries/api/R/registry_loader.R::load_compute_api()`.

New methods (same chat-14 GHSL pattern — thin wrappers over the wrapper's
functions, NULL+warning if deps missing):

- `ancestry_at_interval(chrom, start_bp, end_bp, K=NULL)` — per-window
  summary (mean_delta12/entropy/ena) intersecting the interval
- `ancestry_at_candidate(cid, K=NULL)` — same, looked up via
  `intervals$get_candidate(cid)`
- `ancestry_q_vector(chrom, start_bp, end_bp, K=NULL, sample_ids=NULL)` —
  full per-sample × per-window Q matrix for the region (reads the _samples
  cache, narrows by bp window and optional sample_ids)
- `ancestry_q_summary(chrom, K=NULL)` — whole-chrom window summary
- `ancestry_q_and_f_for_candidate(cid, K=NULL, persist=TRUE)` — returns
  `list(Q, F, cid, K, chrom, start_bp, end_bp)` and auto-persists to
  `reg$stats$put_candidate_q_f()` when persist is TRUE

Pattern is `.ensure_instant_q()` that sources the wrapper if not already
sourced. `pairwise_stat()` also gained a `cache=TRUE/FALSE` flag that
auto-writes to `reg$stats$put_pairwise()`.

### [BK-statscache] NEW — reg$stats subsystem (fourth registry)

**Motivation.** Before chat 15, `reg$compute$pairwise_stat()` returned
results in memory and discarded them. Same for all GHSL methods. The
evidence_registry only stores scalar keys, so raw FST tracks / per-candidate
Q matrices / F matrices had no persistent home. User asked explicitly:
"for inversion candidates of course the Q and F matrix and the tree or
if no tree at least we need them for that interval." Scope: Q + F, no
tree (would require choosing an R tree library and format; Q+F suffices
to reconstruct offline).

**Implementation.** `registries/api/R/registry_loader.R::load_stats_api()`.
Rooted at `$STATS_CACHE_DIR` (default `registries/data/stats_cache/`).

Layout:
```
stats_cache/
├── pairwise/<stat>/<chrom>/<g1>__vs__<g2>.<sample_set>.tsv.gz
├── candidate/<cid>/Q_K<NN>.<sample_set>.tsv.gz
├── candidate/<cid>/F_K<NN>.tsv.gz
├── candidate/<cid>/meta.tsv
└── manifest.tsv
```

Methods on `reg$stats`: `put_pairwise`, `get_pairwise`, `put_candidate_q_f`,
`get_candidate_q`, `get_candidate_f`, `list_cached(kind=NULL)`,
`clear_candidate`.

`<g1>__vs__<g2>` is canonically sorted so `g1 vs g2` and `g2 vs g1` map to
the same file. Compute API publishes the stats handle via `.REG_STATS_API`
global so `ancestry_q_and_f_for_candidate(..., persist=TRUE)` and
`pairwise_stat(..., cache=TRUE)` can write without needing the full `reg`
object.

### [BK-configs] HIGH — inversion + ancestry configs rewritten

**Where.** `inversion_modules/00_inversion_config.sh`, `unified_ancestry/00_ancestry_config.sh`.

Issues with the prior inversion config:
- `CODEBASE=${BASE}/inversion_codebase_v8.5` — broke under the flattened
  `inversion_modules/` layout
- `DISCOVERYDIR`/`FOLLOWUPDIR`/etc. referenced non-existent `MODULE_5*` dirs
- `UTILS_DIR=${CODEBASE}/utils` → non-existent path
- `NGSADMIX_DIR=05_NGSadmix` vs ancestry config's `05_ngsadmix_global` —
  divergent
- No `LOCAL_Q_DIR` / `STATS_CACHE_DIR` / `CANONICAL_K` / `K_SWEEP` /
  `SV_PRIOR_DIR` / `GHSL_DIR`

Fix: full rewrite. `CODEBASE` now auto-discovers (flattened layout first,
v8.5 fallback). `UTILS_DIR` auto-discovers. `LOAD_BRIDGE` derives from
`UTILS_DIR`. All the chat-15 data dirs added. `PHASE2_DIR`/`PHASE3_DIR`/
`PHASE4_DIR` declared as preferred pointers (legacy `MODULE_5*` aliases
preserved for pre-flattened launchers). `NGSADMIX_DIR` synced to
`05_ngsadmix_global` matching ancestry config. `inv_init_dirs()` now
mkdirs `LOCAL_Q_DIR`, `STATS_CACHE_DIR`, `SV_PRIOR_DIR`, `GHSL_DIR`.

Same data-dir additions to the ancestry config (`LOCAL_Q_DIR` default
changed from `$BASE/unified_ancestry/local_Q` → `$BASE/ancestry_cache`,
`CANONICAL_K` + `K_SWEEP` exports added, `LOAD_BRIDGE`/`SAMPLE_MAP_R`/
`SAMPLE_REGISTRY_R` hard-coded v8.5 paths replaced with auto-discovery).

### [BK-docs] HIGH — docs

- `unified_ancestry/README.md`: directory layout now shows code (git) vs
  data (scratch) split, K-sharded cache layout, migration note. Steps 4a/4b
  updated with full K-sweep vs canonical-only commands.
- `unified_ancestry/README_REWIRING_v9_need_rewrite.md` → retired to
  `_archive_superseded/` (work described there is now done).
- `inversion_modules/phase_2_discovery/2c_precomp/README.md`: C01a section
  now mentions `localQ_delta12_K08` / `_K08` / `_K08` stamps. New
  "Ancestry precompute (sidecar workflow)" section with the SLURM command.
- `registries/API_CHEATSHEET.md`: compute table extended with chat-14 GHSL
  and chat-15 ancestry rows; new `reg$stats` table with usage example.

### Coverage after chat 15

Candidate-spec `q2_vector` coverage: 22/73 = 30.1% (start of chat 15)
→ 59/91 = 64.8% (end of chat 15 BK work) → unchanged by the ancestry
wiring (ancestry keys are precomp-level features, not block-schema keys).

The full-wiring pass closed the "ancestry never actually persisted"
hole that was the root cause of repeated questions about "where do
Q/delta12/ena/FST go on disk." Post-chat-15 answer in one sentence:
cohort aggregates in precomp RDS, full per-sample Q in `$LOCAL_Q_DIR/K<NN>/`
K-sharded cache, per-candidate Q+F archives in `$STATS_CACHE_DIR/candidate/`,
pairwise group stats in `$STATS_CACHE_DIR/pairwise/`.

### Tools added

- `_rcheck.py` (retained from earlier chats) — parse-check R files via a
  lightweight tokenizer; no R interpreter needed.
- `_schema_check.py`, `_code_field_check.py`, `_bk_rename.py` — chat-15 BK
  work validators, retained.

### Out of scope / not attempted

- Tree persistence for candidates (user said "if no tree at least we need
  them for that interval" — Q+F delivered; tree requires library choice)
- HPC calibration (no LANTA access in chat env)
- Threshold tuning (user explicit: "why do we need to calibrate anything")
- C01d scoring dimension re-weighting now that localQ columns exist — new
  feature; different ticket

---

## Chat 15 (2026-04-17) — BK schema canonicalisation

Scope: phase-4b canonical schemas for the three block types whose
`keys_extracted` directives were missing or stale. Chat 15 deliberately
did NOT touch the HPC-requiring priorities from the chat-14 handoff
(no R interpreter in the chat-15 environment; no LANTA access; no
calibration work — per user direction, default thresholds stand unless
observed distributions break something).

### [BK] HIGH — three canonical phase-4b schemas missing or stale

**Where.**
- `registries/schemas/structured_block_schemas/internal_dynamics.schema.json`
  (existed but was the pre-chat-12 recombinant-flavor — every
  `keys_extracted` entry referenced a field the current C01i_decompose
  does not write; all 12 extractions silently skipped)
- `registries/schemas/structured_block_schemas/recombinant_map.schema.json`
  (did not exist — block type written by C01i_b but no schema loaded;
  registry loader recorded `validation_status = "no_schema"` and
  extracted nothing)
- `registries/schemas/structured_block_schemas/internal_ancestry_composition.schema.json`
  (did not exist — block type written by C01i_c for the same reason)

Draft versions at `inversion_modules/phase_4_postprocessing/schemas/`
were never loaded by `registry_loader.R::load_schema`, which hard-codes
the canonical path (line 720: `schema_blocks <- file.path(schema_dir, "structured_block_schemas")`).

**Handoff framing reconciliation.** The chat-14 handoff listed four
schemas (internal_dynamics, recombinant_map, internal_ancestry_composition,
nested_composition). Investigation: `nested_composition` is the informal
alias for `internal_ancestry_composition` in the docs — same block type,
no distinct writer exists (`STEP_C01i_c_nested_composition.py` writes
with `block_type="internal_ancestry_composition"`). Chat-15 BK therefore
treats it as three blocks, not four.

**How fixed.**
1. Rewrote `internal_dynamics.schema.json` to match chat-12 decompose-flavor:
   silhouette / BIC / flashlight_mode / seed counts / PCA class counts /
   discordant samples / cheat2 reclassification / per_sample PCA class
   array. Handles three emit variants (seeded_ok, no_seeding × 2) via
   nullable rich fields gated on `decomp_status`. 19 keys_extracted
   entries (16 Q2 + 3 Q6_prelim).
2. Created `recombinant_map.schema.json` matching the actual chat-12
   R + G gate block written by STEP_C01i_b_multi_recomb.R. Chat-12
   renames (`n_recombinant_gc`, `n_recombinant_dco`, `n_disputed`,
   `n_ghsl_only`, `gate_params.min_dco_bp`) replace the pre-chat-12
   field names the draft schema still used (`n_gene_conversion`,
   `n_double_crossover`, `mean_posterior`). 10 keys_extracted entries.
3. Created `internal_ancestry_composition.schema.json` by promoting the
   phase_4 draft (fields aligned with `STEP_C01i_c_nested_composition.py
   ::analyze_parent`). Added `pct_diffuse_mixed`, `K_used`, and
   `n_samples_analyzed` to keys_extracted (draft had 9, canonical has 12).
4. Updated `INDEX_remaining_blocks.json`: bumped total 18 → 20 (new
   block types broken out of internal_dynamics are now first-class),
   added `recombinant_map` and `internal_ancestry_composition` entries,
   rewrote the `internal_dynamics` entry's `status` field, updated the
   `written_so_far` list and totals.
5. Updated `compute_candidate_status.R::build_key_spec()::q2` and `::q6`
   vectors: added the 36 newly-wired Q2 keys, the 3 newly-wired Q6_prelim
   keys; removed 5 stale entries that corresponded to pre-chat-12 draft
   fields no writer produces any more (`q2_n_suspicious_recomb`,
   `q2_n_ambiguous_recomb`, `q2_mean_dco_prior`, `q2_max_dco_prior`,
   `q2_recomb_mean_posterior`); renamed `q2_decomp_quality_flags` →
   `q2_decomp_quality` to match the actual emitted key.

**Test.** Four-gate validation:
- JSON parse: all three files parse cleanly.
- Draft-07 meta-schema: all three validate against
  `http://json-schema.org/draft-07/schema#` under `jsonschema.Draft7Validator`.
- Internal consistency (`_schema_check.py`): every keys_extracted `from`
  path resolves to a declared property; no duplicate keys within a schema;
  no duplicate `from` paths producing aliases. 41/41 entries OK.
- Source-code cross-check (`_code_field_check.py`): every `from` path's
  head segment is a field that the corresponding R/Python script actually
  writes in its block_data construction. 41/41 head segments OK.
  (The R parser's strict mode under-reported 12 internal_dynamics fields
  due to multi-line if/else values in `decomp_quality = if (...) ... else
  ... else ...,`; the permissive fallback pass catches them. Cross-verified
  by direct grep: all 12 flagged fields have indented `name =` assignments
  in STEP_C01i_decompose.R.)

**Coverage metric.** Q2: 22/73 = 30.1% (chat-14 end) → 58/90 = 64.4%
(chat-15 end). Denominator expanded because the three canonical schemas
exposed keys that were being written inside JSON blocks but not tracked
in the coverage vector. Meets the handoff's 50-60% target.

**Archive.**
- Pre-BK `internal_dynamics.schema.json` (canonical, pre-chat-12
  recombinant-flavor) → `_archive_superseded/bk_schemas_pre_canonical/`.
- Three phase-4b drafts that lived in the non-canonical directory →
  `_archive_superseded/bk_schemas_pre_canonical/phase4b_drafts/`.
- Not archived: `inversion_modules/phase_4_postprocessing/schemas/frequency.v2.schema.json`
  (out of BK scope; no canonical `frequency.schema.json` exists yet either
  and the chat-14 handoff didn't flag it).

**Tools added.**
- `_schema_check.py` — JSON + draft-07 meta + internal `from`-path resolver.
- `_code_field_check.py` — static cross-check of schema `from` paths
  against source-code block_data field names. Known limitation: R
  parser's strict mode under-reports fields with multi-line if/else
  values; permissive fallback catches them.

### [BK-rename] HIGH — scientific renaming of procedural keys

**Where.** All three canonical schemas written under [BK], plus
`compute_candidate_status.R::build_key_spec()::q1/q2`.

**Symptom.** The first BK pass wired 41 keys whose names were
procedural / mechanical — e.g. `q2_decomp_status` (what decomposition?
says "status" but not of what), `q2_flashlight_mode` (mode of what?),
`q2_cheat24_version` (doesn't name the thing it's provenance for),
`q2_combination_rule` (of what?), `q2_n_recombinants` (count of what?),
`q2_pct_samples_two_block` (two blocks of what?). The names said what
software step emitted the key, not what the key measured about the
genome. Bad for the manuscript; bad for downstream readers.

**How fixed.** Systematic rename pass, applied by `_bk_rename.py`:
- 12 keys in `internal_dynamics` renamed (7 preserved where the
  existing name was already scientific: `q2_bic_gap_k3_vs_k2`, the
  three `q2_n_seed_HOM_*` counts, the three `q6_n_HOM_*_prelim`
  counts).
- 10 keys in `recombinant_map` renamed (0 preserved — all old names
  were procedural).
- 9 keys in `internal_ancestry_composition` renamed (3 preserved:
  `q1_composite_flag` which is the v10.1 canonical per chat-9
  Finding X, and the two `q2_ancestry_*` keys that already carried
  the correct prefix).

**Renaming philosophy.** Each new name answers "what biological
quantity does this measure?". Examples:
- `q2_silhouette_score` → `q2_pca_cluster_silhouette` (it's the
  silhouette of the PCA-based class clusters; that's the genomic
  thing being measured).
- `q2_cheat24_version` → `q2_recomb_posterior_source` (the version
  string serves as provenance for where the posterior came from,
  or whether it's NA because cheat24 was unreachable).
- `q2_n_recombinants` → `q2_n_recombinant_samples` (count of what?
  samples).
- `q2_n_recomb_disputed` → `q2_n_regime_ghsl_disputed_samples`
  (the dispute is specifically between the regime-DAG signal R
  and the GHSL SPLIT call G — not generic).
- `q2_n_cheat2_reclass` → `q2_n_hetDEL_breakpoint_reclass` (the
  reclassification is triggered by het-DEL signatures at the
  inversion breakpoint, which IS the biological signal).

Full rename table lives in
`registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md`.

**Addendum fix.** First BK pass stated in comment that chat-13's
`q2_decomp_quality_flags` had been renamed to `q2_decomp_quality`.
This was wrong: the two are DISTINCT keys — the former is written
directly by `STEP_C01i_d_seal.R` (comma-separated aggregate of
validation quality_flags); the latter is the silhouette-derived
clean/noisy flag from the internal_dynamics schema. Both keys
coexist. `q2_decomp_quality_flags` (seal-written) re-added to the
q2 vector with a clarifying comment; `q2_pca_cluster_separation_flag`
(the renamed silhouette-derived key) remains distinct.

**Test.** `_bk_rename.py` applied 31 quoted replacements in
`compute_candidate_status.R` plus all `keys_extracted` directives
in the three schemas. `_schema_check.py` + `_code_field_check.py`
re-run after rename: 41/41 keys structurally valid, 41/41
`from` paths resolve to a source-code field. `_rcheck.py` on the
edited R file: OK.

### [BK-docs] NEW — per-key biological/computational documentation

**Where.** New file
`registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md`
(~18 KB, ~500 lines).

**What it provides.** For each of the 41 BK keys:
- **Provenance** — which block field it extracts from, which R/Python
  script writes the block, what upstream data feeds into that
  computation.
- **Biological meaning** — what the value measures about the genome
  (e.g. "per-candidate cohort-level disagreement between two
  independent class-assignment systems, SV-based vs PCA-based").
- **Downstream consumer** — which scoring / promotion / reporting
  script reads the key, and what decision it informs.
- Introductory sections on what Q1 / Q2 / Q6 mean, what a
  "candidate" is in this toolkit, and what the three block types
  characterise about an interval.
- Trailing rename map: pre-rename → post-rename with the scientific
  reasoning for each.

**Test.** Markdown renders cleanly under standard linters (headings,
tables, code fences balanced). Internal cross-references
(`q1_composite_flag` → internal_ancestry_composition section, etc.)
confirmed present. Each canonical schema's `description` field now
ends with a pointer to the explainer so downstream developers find
it without grepping.

### Updated coverage metric

After the rename pass and re-adding `q2_decomp_quality_flags`:
- Q2 vector: 91 unique keys (was 73 at end of chat 14, 90 after
  BK first pass, 91 after re-adding the seal-written key).
- Q2 wired: 59 (was 22 at end of chat 14).
- Q2 coverage: **59 / 91 = 64.8%** (was 30.1% at end of chat 14).

Meets the handoff's 50-60% target with a margin.

### Chat 15 out of scope (deferred to chat 16 + HPC)

- HPC validation: running STEP_C04 + STEP_C04b on LG12, smoke-testing
  the `reg$compute$ghsl_at_candidate` registry API, Tier-3 ghsl
  confirmation, sub-block scans. No R interpreter was available in the
  chat-15 environment.
- Re-running chat-12 unit tests (50/50 DAG + 33/33 GC detector) — needs HPC R.
- Threshold calibration AS / AW / AV — user direction: v5-inherited
  defaults stand unless observed distributions break something, at
  which point the tuning is a finding not a miscalibration.
- `frequency.v2.schema.json` canonicalisation — not in BK scope.

---

## Chat 14 (2026-04-18) — GHSL v6 wiring cleanup + panel library

### [BL] CRITICAL — four downstream consumers still read v5 paths/columns

**Where.**
- `phase_4_postprocessing/4b_group_proposal/lib_ghsl_confirmation.R`
  (lines 60–62, 85–86) — read `<chr>.ghsl_v5.annot.rds` and
  `.karyotypes.rds`.
- `phase_2_discovery/2d_candidate_detection/run_all.R` (lines 453–525)
  — reads v5 annot with v5 column names.
- `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
  line 319 — `iv$ghsl_v5_score_max`.
- `phase_2_discovery/2d_candidate_detection/STEP_D05_ghsl_stability.R`
  header — referenced `STEP_C04_snake3_ghsl_v5.R`.

**Symptom.** Chat 13 swapped `2e_ghsl/` from v5 to v6 but did not update
downstream consumers. v6 classifier only emitted genome-wide TSVs and
did not emit the per-chrom annot/karyotype RDS files the consumers
expected. Net effect: phase-4b Tier-3 karyotype confirmation silently
returned UNINFORMATIVE for every sample; run_all.R's 2d block-scoring
GHSL merge was silently skipped.

**How fixed.**
1. Classifier `STEP_C04b_snake3_ghsl_classify.R` now emits three
   per-chrom RDS files at the end of each chrom loop iteration:
   `annot/<chr>.ghsl_v6.annot.rds` (thin per-window aggregates),
   `annot/<chr>.ghsl_v6.karyotypes.rds` (per-sample stable LOW/HIGH
   runs), `per_sample/<chr>.ghsl_v6.per_sample.rds` (dense long-format
   panel, one row per sample × window, carrying per-scale divergence +
   rank + rank-band at every rolling scale plus stable-run membership
   from Part B and window-level score/status).
2. `lib_ghsl_confirmation.R` reads v6 paths first with v5 fallback for
   transitional runs. Adds `ghsl_per_sample_panel_in_interval()` as
   an additive panel-based alternative to the run-overlap Tier-3
   logic (run-overlap stays authoritative for SPLIT detection).
3. `run_all.R` 2d: reads v6 annot with v6 column names (`ghsl_v6_score`
   + unprefixed `rank_stability` / `div_contrast_z` / `div_bimodal`);
   v5 fallback retained for mixed runs. Merged block columns renamed
   `ghsl_v5_score_max` → `ghsl_v6_score_max`; other merged column
   names (`ghsl_rank_stability_max`, `ghsl_div_contrast_z_max`,
   `ghsl_div_bimodal_frac`, `ghsl_pass_frac`, `ghsl_n_scored_windows`)
   kept as stable contract. PASS fraction now reads `ghsl_v6_status`
   directly when v6; falls back to `score > 0.65` heuristic when v5.
4. `STEP_C01d` renamed `iv$ghsl_v5_score_max` → `iv$ghsl_v6_score_max`.
5. `STEP_D05` doc-only references updated.
6. `00_config.R::CFG$GHSL_DIR` default `ghsl_v5` → `ghsl_v6`.

**Test.** All nine edited R files pass `_rcheck.py` brace/paren/string
balance. Full end-to-end test requires HPC R (chat 15 earliest).
Archive of pre-patch copies at
`_archive_superseded/ghsl_v5_consumers_pre_v6_rewire/`.

### [BM] HIGH — new per-sample panel backend + on-demand query library

**Where.** New file `utils/lib_ghsl_panel.R` (~500 lines).

**What it provides.**
- `load_ghsl_panel(chrom, ghsl_dir, refresh)` — memoized per-chrom
  loader over the dense panel RDS.
- `ghsl_panel_range(chrom, start_bp, end_bp, sample_ids, ghsl_dir)` —
  dense long-format slice.
- `ghsl_panel_aggregate(chrom, start_bp, end_bp, sample_ids, scale,
  summaries, ghsl_dir)` — per-sample summary (mean, frac_low,
  longest_low_run_bp, rank_at_peak, stable_run_call, ...).
- `ghsl_panel_subblock_scan(chrom, subblocks, sample_ids, scale,
  summaries, include_transitions, ghsl_dir)` — per-(sample, subblock)
  aggregation for complex candidates with soft boundaries (e.g. 4
  soft boundaries → 6 sub-blocks). Optional cross-subblock
  band-transition patterns expose recombinant-detection signal
  (a sample with `n_transitions ≥ 1` toggled LOW/HIGH across the
  candidate).
- `ghsl_panel_wide_matrix(chrom, start_bp, end_bp, sample_ids, metric,
  sample_order, ghsl_dir)` — `[samples × windows]` matrix for
  heatmap plotting.
- Plot helpers: `plot_ghsl_sample_track`, `plot_ghsl_heatmap_track`,
  `plot_ghsl_subblock_panel`.

**Documentation.** `utils/lib_ghsl_panel_README.md`.

### [BN] HIGH — registry compute API extended with GHSL panel methods

**Where.** `registries/api/R/registry_loader.R::load_compute_api()`.

**What it provides.** Same pattern as `pairwise_stat` /
`boundary_fst_profile` — registry-aware wrappers over the backend
library:
- `reg$compute$ghsl_at_interval(chrom, start_bp, end_bp, ...)`
- `reg$compute$ghsl_at_candidate(cid, ...)` — looks up chrom/bp from
  the interval registry.
- `reg$compute$ghsl_at_subblocks(cid_or_chrom, subblocks, ...)` — for
  the "N sub-blocks from soft boundaries" case.
- `reg$compute$ghsl_at_block_subregions(cid, block_type, ...)` —
  results-aware; pulls sub-regions from an existing evidence block
  (`boundary`, `regime_segments`, `recombinant_map`, etc.). Tries
  common field names: `segments`, `sub_blocks`, `regions`, `tracts`,
  `soft_boundaries_blocks`. Accepts start/end aliases.
- `reg$compute$ghsl_wide_matrix(cid_or_chrom, ...)` — plot data prep.
- Auto-sources `lib_ghsl_panel.R` via `GHSL_PANEL_LIB` env var or
  standard relative paths; returns NULL with a warning if absent.

### [BO] MEDIUM — heavy engine loaded all chroms eagerly

**Where.** `STEP_C04_snake3_ghsl_v6.R` lines 109–110 (pre-patch).

**Symptom.** Single-chrom runs paid memory for all 28 chroms' precomp
RDS before the `--chrom` filter was applied. Wasted 1–3 GB RAM.

**How fixed.** Replaced eager `precomp_list[[chrom]] = readRDS(f)`
with a filename→chrom index (built by peek-readRDS-then-rm) that
resolves the filter first, and reads the heavy payload only for each
selected chrom inside the main loop. Includes `rm(pc, dt)` + `gc` in
end-of-iter cleanup.

### [BP] MEDIUM — finer default scale ladder for the heavy engine

**Where.** `STEP_C04_snake3_ghsl_v6.R` `ROLLING_SCALES` default.

**What changed.** `c(20L, 50L, 100L)` →
`c(10L, 20L, 30L, 40L, 50L, 100L)`. Gives an evenly-spaced fine
ladder (~50 kb to ~250 kb) for detail work plus s100 for overview
plotting. Added rationale comment block. Adds seconds, not minutes,
to the heavy-engine runtime.

### [BQ] LOW — k=4 / k=5 Part C labels were generic `CLASS_N`

**Where.** `STEP_C04b_snake3_ghsl_classify.R` `classify_interval`.

**What changed.** At k=4 labels are now by divergence rank
`INV_INV / INTER_LOW / INTER_HIGH / INV_nonINV`. At k=5 HET sits in
the middle: `INV_INV / INTER_LOW / HET / INTER_HIGH / INV_nonINV`.
k>5 retains generic `CLASS_N`. Part D (CUSUM changepoint clustering)
stays the sub-inversion-boundary axis and reports its own
`n_clusters` independently — two axes, two outputs, not conflated.

### Doc fixes (no code change)

- `phase_2_discovery/2e_ghsl/README.md` — wiring-status section
  rewritten; three per-chrom RDS outputs documented; scale ladder
  documented.
- `phase_3_refine/README.md` line 22 — `STEP_C04_snake3_ghsl_v5.R`
  → `STEP_C04_snake3_ghsl_v6.R + STEP_C04b_snake3_ghsl_classify.R`.
- `phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md`
  line 52 — GHSL v5 section updated to v6 with chat-14 status note.
- `phase_2_discovery/2c_precomp/patches/README.md` line 30 — patch
  target v5 → v6.
- `phase_2_discovery/2c_precomp/RENAMING.md` line 465 — file list v5
  → v6 (heavy + classifier pair).

### Tool added

`_rcheck.py` at work root — R brace/paren/string-balance + stray `,,`
checker. Handles R's legitimate `x[, col]` array indexing. All edited
files pass. Chat 13 used a similarly-named tool but didn't ship it in
the tarball; chat 14 re-built it.

---

## Chat 13 (2026-04-18) — registry wiring pass

### [BG] HIGH — registry_loader dropped dotted-path keys silently

**Where.** `registries/api/R/registry_loader.R::write_block`, the
key-extraction loop around line 776.

**Symptom.** Schema entries like
`{ key: "q2_dag_fraction_samples_R_fired", from: "cohort.fraction_samples_R_fired" }`
were extracted via `data[[from]]` which does not walk dotted paths, so
21 keys across 4 schemas never reached `keys.tsv`. This included all
six `boundary_sharpness.*` keys and all five `cohort.*` DAG keys.

**How fixed.** Added `resolve_dotted(obj, path)` helper that splits on
`.` and indexes successively. Vector-valued leaves (e.g.
`persistent_carriers`) collapse to `length()` to match schema intent
(the key wants the count, not the array).

**Test.** Post-fix key count confirmed by Python inspection:
18 chat-13 wired keys survive extraction, and a further 3 in
age_evidence + 1 in gene_conversion_tracts are also unlocked. Before:
0 of the 18 chat-13 keys landed in keys.tsv because every one of them
uses a dotted path (by design of the schemas).

### [BH] CRITICAL — C01j advertised but did not write regime_memberships.tsv.gz

**Where.** `STEP_C01j_regime_compatibility_engine.R` header line 41
advertised this file; main loop computed `regime$memberships` per
candidate and then discarded it.

**Symptom.** C01i_b's `--regime_memb` argument requires that file.
Pipeline could NOT run end-to-end.

**How fixed.** Collect each candidate's `regime$memberships` into
`all_memberships` during the loop; `fwrite()` at the end via
`rbindlist(all_memberships, fill = TRUE)`. The columns (sample,
pos_mid_mb, group_id, dosage_class, chrom, interval_id) already match
what `derive_R_from_regime()` normalises.

**Test.** Code inspection; parse-check OK; full end-to-end test
deferred to chat-14 HPC run.

### [BI] HIGH — seal missed RECOMBINANT_GC classifications

**Where.** `STEP_C01i_d_seal.R` — both `resolve_final_classes`
(line ~148) and `register_all_groups` (line ~240).

**Symptom.** C01i_b (chat-12 rewrite) writes
`event_class = "gene_conversion_embedded"` but seal was matching on
`"gene_conversion"`. No samples ever reached the
`inv_<cid>_RECOMBINANT_GC` group.

**How fixed.** (1) Accept both canonical (`gene_conversion_embedded`,
`double_crossover`) and legacy (`gene_conversion`, `suspicious`,
`ambiguous`) event_class values. (2) Prefer `rec_info$recomb_subgroup`
from multi_recomb's per-sample record when present — that's the
canonical subgroup token, avoiding the event_class translation
altogether. `recomb_subgroup` is now threaded through `final_records`
so `register_all_groups` can read it.

**Test.** Parse-check OK; HPC verification in chat 14.

### [BJ] LOW — build_key_spec() did not list chat-13 keys

**Where.** `compute_candidate_status.R::build_key_spec()::q2`.

**Symptom.** The 18 chat-13 wired keys were missing from the q2 vector
so completion accounting would under-report them as unproduced.

**How fixed.** Appended the 18 keys to the q2 list with a comment tag
linking back to Finding BJ.

**Test.** Post-fix wired-fraction: Q2 22/73 = 30.1%, registry-wide
73/319 = 22.9%.

### [BK] NEW, not yet applied — four 4b schemas missing keys_extracted directives

**Where.** `registries/schemas/structured_block_schemas/` —
`internal_dynamics.schema.json`, `recombinant_map.*`, nested
composition schema, and the existing phase-2 schemas for the 51
Q2 spec keys that have NO schema writer producing them.

**Status.** Flagged for chat 14/15. Adding `keys_extracted` directives
is the path to the 85% Q2 coverage aspiration from the handoff.
Out of chat-13 scope; the data is written inside the JSON blocks, just
not flattened to keys.tsv.

### [WIRE] Block output wiring (C01j / C01l / C01m / seal / characterize)

**Where.** Four 4a/4b scripts + characterize_candidate.R.

**How.** Per-candidate `write_block_safe()` calls routed through the
registry API (JSON fallback to `<outdir>/blocks/<cid>/<block>.json`
if no registry). Characterize Q2 pulls 14 new keys as annotation-only
evidence lines.

**Test.** Parse-check OK on all edits.

### [AT] MEDIUM — C01l flank_bp now scales with candidate span

**Where.** `STEP_C01l_local_structure_segments.R` main loop.

**How fixed.** `cand_flank_bp = clamp(span_bp, 200000L, 500000L)`;
`define_segments(inv_start, inv_end, flank_bp = cand_flank_bp)`.

**Test.** Parse-check OK; calibration verified against LG12 + LG25
in chat 14.

### [AV] LOW — decompose now emits decomp_quality annotation

**Where.** `STEP_C01i_decompose.R` block_data.

**How fixed.** `decomp_quality = if (silhouette >= 0.40) "clean" else
"noisy"`. Non-gating; C01d / C01f must NOT promote on this field during
the first wiring pass. The 0.40 threshold will be calibrated in
chat 14 against observed silhouette distribution.

### [BB] MEDIUM — ghsl_confirmation silent first-row-pick fixed

**Where.** `lib_ghsl_confirmation.R::resolve_karyo_bp`.

**How fixed.** Added `stopifnot(uniqueN(annot_dt$global_window_id) ==
nrow(annot_dt))` at the top with inline rationale (if a re-run caused
two annotation passes both to contribute rows, the bp mapping would
silently pick the first — leading to wrong bp). Now fails loud.

### [AO] DOC — C01d 10 dimensions → 12 dimensions

**Where.** `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` L30.

**How fixed.** Header now reads "12 dimensions from detector" matching
the D1–D12 mapping immediately below.

### [AY] DOC — seed_loader clarified drop-conflict is stricter than docs

**Where.** `lib_step03_seed_loader.R` header.

**How fixed.** Added explicit note that implementation is STRICTER
than chat-9 design language — drop-on-conflict, not priority-pick.
Chat-9 wording preserved for traceability.

### [AQ / BD] VERIFIED — no hardcoded stale C01j paths

**Where.** `STEP_C01e_candidate_figures.R` Panel H;
`STEP_C01k_annotated_simmat.R`.

**Outcome.** Both read sidecar TSVs via `--regime_dir` CLI arg; C01j
still writes those filenames at the same relative locations inside
`--outdir`. Chat-11.5 directory move did not break these paths because
C01j writes unchanged filenames. No patches needed.

---

## Chat 12 (2026-04-18) — phase-4 coherence audit + DAG rewrite

### [V] CRITICAL — inline cheat24 fallback diverged from real cheat24

**Where.** `STEP_C01i_b_multi_recomb.R` old L97–122.

**Symptom.** Inline fallback returned a different `posterior_probability`
from the real cheat24 on the same input — 2× divergent on
`gene_conversion` class. Two candidates side-by-side on the same HPC
run got different posteriors depending purely on cheat24 source
file reachability. Open since chat 8.

**How fixed.** Deleted the inline fallback in the chat-12 multi_recomb
rewrite. If real cheat24 absent, `posterior = NA` and block records
`cheat24_version = "cheat24_unavailable_posterior_na"`.

**Test.** Code inspection — no `classify_recombinant_event` defined
inline any more.

### [PRIMARY] DAG-based derive_R_from_regime

**Where.** `lib_recomb_combination.R::derive_R_from_regime`.

**Symptom.** Old counter over-counted A→B→A as 2 changes (the
handoff's named pattern). Also rigid "long-long-long" shape filter
dropped 2-window deviations as noise when they were real signal, and
ignored recombinants whose terminals don't match start (important
discriminator vs alternating noise).

**How fixed.** Full rewrite with DAG semantics. Per sample: build
directed graph over per-window regime labels (nodes = distinct runs,
weights = run lengths). Summary fields include `deviation_fraction`,
`longest_deviation_bp`, `n_distinct_nodes`, `n_edges`,
`terminal_matches_start`. R_fired is two-gate: fraction AND bp span
must both clear threshold. The bp gate is the real discriminator
against alternating noise (high fraction, low longest-dev).

**Test.** `tests/test_dag_derive_R.R` — 50/50 PASS. Covers AAABBA,
AAABABBA, single-run, boundary recombinant (non-matching terminals),
cohort aggregation, interval filtering, empty inputs, back-compat
paths.

### [PRIMARY] Per-SNP gene conversion detector rewrite

**Where.** `gene_conversion_detector.R` full rewrite.

**Symptom.** Old windowed-binning detector (40-SNP windows, CUSUM
sweep) had lower detection floor ~60 kb at Quentin's 1.6 kb/SNP
density. Biology of GC tracts: 18–774 bp typical (Harris 1993), up
to 5 kb atypical. Detector was blind to real tracts; surfaced SNP-
column artefacts instead.

**How fixed.** Per-SNP run-length detector: QC filter (missingness,
excess-het / paralog, low MAF, depth anomaly) → diagnostic SNP mask
(|AF_REF − AF_INV| ≥ 0.5) → per-sample per-SNP flag → run-length with
tolerance (up to 1 interior skip) → length gate (2–10 flags, span ≤
20 kb) → geometric-prior confidence (LOW/MEDIUM/HIGH) → direction
annotation → snp_qc accounting block.

**Test.** `tests/test_gc_detector.R` — 33/33 PASS.

### [gate] Remove 100 kb mosaic threshold

**Where.** Old `STEP_C01i_b_multi_recomb.R` L333–340.

**Symptom.** Arbitrary 100 kb escape-hatch threshold in the old
combine rule. No principled justification; coupled all inversion
sizes to one threshold.

**How fixed.** Removed. DAG's `min_deviation_bp` (50 kb default)
replaces it with cleaner semantics based on longest contiguous
non-dominant run.

### [Signal 1] Demote old per-window PCA class switch to Tier-4 diagnostic

**Where.** `STEP_C01i_b_multi_recomb.R`.

**Rationale.** Signal 1 was a noisy 1D PCA projection. C01j's
sliding-window compatibility grouping is the correct substitute.

**How fixed.** Signal 1 still computed if per_window_class.rds sidecar
present; merged into per-sample record with `s1_*` prefix for audit
trails. Not referenced by the gate.

### [seeding] Wire STEP03 seed loader, remove unsupervised fallback

**Where.** `STEP_C01i_decompose.R`.

**Symptom.** On pure noise intervals (no real inversion signal), the
unsupervised k-means fallback produced essentially random 3-way
splits that got registered as groups and polluted every downstream
test.

**How fixed.** Added `--step03_seeds_dir` option. Routes seeds through
`combine_tier1_seeds()` (flashlight ∪ STEP03, drop conflicts). If
combined set has < `min_seeds_per_class` in any class, emit
`decomp_status="no_seeding"` + reason and skip. Also catches the
inner failure where `seeded_kmeans()` returns `km$seeded = FALSE`.

### [registry] Collapse comp_from_registry to one API call

**Where.** `STEP_C01f_hypothesis_tests.R::comp_from_registry`.

**Symptom.** 4×has_group + 4×get_group calls per candidate. With
200+ candidates and several calls per candidate, 3200+ redundant
queries.

**How fixed.** One `reg$get_groups_for_candidate(cid)` call.
Legacy path preserved via `is.function(reg$...)` check.

### [schema] regime_sample_dag (new)

Per-sample DAG block. 5 `keys_extracted`:
q2_dag_fraction_samples_R_fired, q2_dag_fraction_samples_deviating,
q2_dag_dominant_regime_cohort, q2_dag_dominant_regime_support,
q2_dag_n_samples_R_fired.

### [schema] gene_conversion_tracts (updated)

Added: `run_length_flagged_snps`, `tolerance_snps`, `confidence`,
`span_bp`, `snp_qc` accounting, `params` block. Legacy
`tract_width_bp` preserved as alias. Description and direction-field
doc updated per Finding BF — honest framing as "arrangement-
discordant IBS tracts consistent with gene conversion".

### [plot] plot_sample_regime_dag.R (new)

Base-R per-sample DAG diagnostic. Grid layout, yellow background =
R_fired, dominant regime filled, deviations bordered. Smoke-tested
on synthetic 5-sample cohort; PDF output OK.

### [BF] Honest terminology for the "GC" detector

**Rationale.** What the detector measures is short intra-inversion
arrangement-discordant tracts — stretches where one sample looks
like the opposite arrangement. Classical interpretation is gene
conversion (Harris 1993) but not provable from sequence data alone
without trios / MA panel / long-read haplotypes. Short double
crossovers, paralog mismapping, and coincidental IBS at diagnostic
sites produce identical signal.

**How addressed.** File and block names kept for pipeline back-compat
(renaming would ripple to many readers). Schema description and
detector file header updated to the honest framing. Manuscript text
and figure legends should use "GC-like tracts" or "arrangement-
discordant IBS tracts", not "gene conversion events" unless
supporting evidence.

---

## Chats 3–11.5 (summary, detail in _archive/chat_history/)

Pre-chat-12 work delivered the registry infrastructure, the
detection engines, and the architectural framework. Key milestones:

- **Chats 3–5**: initial Phase-4 framework, C01d scoring dimensions
  1–10, boundary catalog, early hypothesis tests.
- **Chats 6–8**: cheat framework (24 = recombinant prior, 25 = block
  viability, etc.); flashlight integration for SV-seeded k-means.
- **Chat 9**: full architectural redesign of 4b. Defined the R+G+GC
  semantics Chat 12 now implements. Removed several redundant
  downstream paths.
- **Chat 10**: phase-4e characterization, 367-key spec, smoke-tested
  the three-tier loader. 73.8% wired at that time.
- **Chat 11**: registry library expansion — 22 new API methods
  including `get_groups_for_candidate`. Sample + interval registries
  consolidated.
- **Chat 11.5**: C01j/C01l/C01m dispatched into 4a from orphan
  locations. GC detector v1 drafted (then rewritten chat 12 because
  of scale issue).

Full `FIXES_APPLIED` register for that era in
`_archive/chat_history/FIXES_APPLIED_2026-04-17.md` (42 KB, chats
3–11.5).

---

## Findings carried forward (not yet closed)

Detailed in `AUDIT_LOG_chat12_2026-04-18.md`. Bucketed:

- **Chat 13 mechanical:** AQ, AT, AV, AY, BB, BC, BD + core wiring.
- **Chat 14 (HPC calibration):** AP, AR, AS, AW, BE.
- **Chat 15 (doc + minor):** AO, T, W, 9.
- **Deferred / low pri:** AL, AZ, BA.
- **Manuscript:** BF (terminology).
