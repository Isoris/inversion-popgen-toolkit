# SESSION_SUMMARY.md

Rolling summary of the inversion popgen toolkit. Current-state only;
historical per-chat narratives live in `_archive/chat_history/`.

**Last updated:** 2026-04-18 (end of chat 17).
**Next chat:** chat 18 — see `HANDOFF_PROMPT_chat18_2026-04-18.md`.

---

## What the toolkit is

Phase-2→5 post-processing pipeline for inversion candidate detection in
a 226-sample North African catfish (*Clarias gariepinus*) cohort (~9×
WGS, 28 chromosomes, LANTA HPC, conda env `assembly`).

Phase-4 is the existence-and-characterization stage, organized in
sub-phases 4a (existence scoring) → 4b (group proposal) → 4c (group
validation) → 4d (group-dependent cheats) → 4e (final classification).
Chat 12 closed the architectural questions; chat 13 did registry wiring
end-to-end; chat 14 cleaned up GHSL v6 and added panel queries; chat 15
added BK schema canonicalisation + the ancestry full-wiring pass +
the half-built `stats_cache`; chat 16 rewrote `stats_cache` as a
first-class fourth registry (`results_registry`) with explicit foreign
keys, JSON schemas, a one-method query plane (`ask()`), and a one-
command integrity check — and fixed the chat-15 R-vs-bash hash mismatch
that was silently no-op'ing every ancestry merge. **Chat 17 did a full
SLURM-launcher-coverage audit, wrote the 15 missing launchers + 1 chain
runner for phases 2c, 2e, 4a, 4b, 4c, 4e, and added two patch scripts
for stale orchestrator paths and the `SV_PRIOR_DIR` config mismatch.
Still no HPC run — that's chat 18.**

## Current architectural state (post-chat-17)

### The four registries (at parity since chat 16)

| Table | Dir | PK | Role |
|---|---|---|---|
| sample_registry   | `data/sample_registry/`   | `group_id` | WHO — samples + named groups |
| interval_registry | `data/interval_registry/` | `candidate_id` | WHERE — windows, candidate intervals, cov paths |
| evidence_registry | `data/evidence_registry/` | `(cid, key)` | WHAT-SCALAR — per-candidate Tier-2 JSON blocks + keys.tsv |
| **results_registry** | `data/results_registry/manifest.tsv` | `row_id (UUID)` | WHAT-NUMERICAL — FST tracks, Q/F matrices, ancestry summaries |

All four have JSON schemas under `registries/schemas/`. All four have
R full-API + Python atomic-writer + bash path-resolver bindings.
`reg$results$integrity_check()` verifies cross-registry FK consistency
in one call.

### Launcher coverage (new in chat 17)

Every phase of the pipeline has a SLURM entry point. Launchers live
under their phase directories (not under a single `launchers/` tree):

- `phase_2_discovery/2c_precomp/`:  4 launchers + `submit_precomp_chain.sh`
- `phase_2_discovery/2e_ghsl/`:     2 launchers (1 array, 1 classify)
- `phase_4_postprocessing/4a_existence_layers/`:     2 launchers (pass-1, boundary)
- `phase_4_postprocessing/4b_group_proposal/`:       4 launchers (decomp, recomb, nested stub, seal)
- `phase_4_postprocessing/4c_group_validation/`:     1 launcher  (C01f)
- `phase_4_postprocessing/4e_final_classification/`: 2 launchers (pass-2, characterize+classify)

All follow the archive template convention: `source ~/.bashrc` →
`mamba activate assembly` → `source 00_inversion_config.sh` with
`set -a` → `source utils/pipeline_bridge.sh` → `inv_init_dirs` → banner
→ `${RSCRIPT_BIN}` call → footer. `-A lt200308 -p compute`. All pass
`bash -n`. See `LAUNCHER_TRIAGE.md` for the full DAG.

### Known path inconsistencies — both patched

1. **Phase-4 orchestrators had stale paths.** `run_phase4.sh` and
   `run_phase4b.sh` pointed at `${BASE}/inversion-popgen-toolkit/…`
   (old tree layout) and `${MODULE_DIR}/phase4b_rewrite/R/…` (old
   subdir name). Patched by `patches/patch_orchestrator_paths.sh`.
2. **`STEP_C00_build_sv_prior.R` ignored `$SV_PRIOR_DIR`.** The R
   script hardcoded output to `${INVDIR}/06_mds_candidates/…/sv_prior`
   while the config exports `SV_PRIOR_DIR=${INVDIR}/03_sv_prior`.
   Patched by `patches/patch_STEP_C00_sv_prior_dir.sh` to honour the
   env var with back-compat fallback.

Both patches are idempotent (detect prior application), write
`.bk_chat17` backups, and have been tested against scratch copies of
the real files.

### The database design (chat 16, unchanged)

**`registries/DATABASE_DESIGN.md`** is the single source of truth.
Summary:

- Every `results_registry` manifest row declares where/who/what/
  how-derived with explicit FK references into sample_registry
  (group_id) and interval_registry (candidate_id).
- No hashes in filenames. Sample-set identity is a registered `group_id`
  (typically `all_226`). Group version = the `created` timestamp in
  sample_registry. Stale caches detected by comparing stored vs current
  version.
- One `reg$results$ask(where, who, what, kind, K, overlap)` method with
  4 filter dimensions + interval-overlap semantics. Shortcuts:
  `ask_what_at`, `ask_what_for_candidate`, `ask_what_for_group`,
  `ask_provenance`.
- `reg$results$integrity_check(check_sha256=FALSE)` — 6 checks
  (FK_group/FK_candidate/group_version_current/file_exists/orphans,
  plus optional sha256_drift). Returns data.table with
  `attr(x, "all_pass")` for script-level gating.
- Sub-cluster naming convention: `inv_<cid>_<KARYO>__sub<N>` /
  `__ancestry<K>_<k>` / `__ghsl_<band>` / `__family_<fid>`. `dimension`
  field in `sample_groups.tsv` classifies each group.

### Registry system (chat 13+ unchanged)

- Sample registry: 22+ methods including `get_groups_for_candidate(cid)`,
  and chat-16's `get_subgroups(cid)`.
- Structured block schemas: every phase-4 block declares
  `keys_extracted`. 41 BK keys across 3 phase-4b schemas from chat 15
  (`internal_dynamics`, `recombinant_map`, `internal_ancestry_composition`)
  — still need HPC validation (chat 18).
- `reg$compute` for GHSL queries (chat 14), ancestry queries (chat 15),
  and multi-window scan drivers (chat 16).

### Ancestry wiring

Chat 16 killed the chat-15 sha1 hash and switched both R and bash to
`SAMPLE_GROUP=all_226` from env. Filenames are
`C_gar_LG12.all_226.local_Q_summary.tsv.gz`. The merge into precomp RDS
now stamps three `localQ_*_K08` columns per window. Chat-17 C01a
launcher probes the expected filename pattern at startup and warns if
absent — first defense against a recurrence of the silent bug.

## Code changed in chat 17

| File                                                      | Change |
|-----------------------------------------------------------|--------|
| `phase_2_discovery/2c_precomp/LAUNCH_STEP_C00_sv_prior.slurm` | **NEW** |
| `phase_2_discovery/2c_precomp/LAUNCH_STEP_C01a_precompute.slurm` | **NEW** — probes localQ files, warns on silent-bug preconditions |
| `phase_2_discovery/2c_precomp/LAUNCH_PHASE_01C_block_detect.slurm` | **NEW** |
| `phase_2_discovery/2c_precomp/LAUNCH_STEP_C01b_1_seeded_regions.slurm` | **NEW** |
| `phase_2_discovery/2c_precomp/submit_precomp_chain.sh` | **NEW** — chains ancestry → C00 → C01a with afterok |
| `phase_2_discovery/2e_ghsl/LAUNCH_STEP_C04_ghsl_v6_compute.slurm` | **NEW** — array 1–28 |
| `phase_2_discovery/2e_ghsl/LAUNCH_STEP_C04b_ghsl_v6_classify.slurm` | **NEW** |
| `phase_4_postprocessing/4a_existence_layers/LAUNCH_C01d_scoring_pass1.sh` | **NEW** |
| `phase_4_postprocessing/4a_existence_layers/LAUNCH_C01g_boundary.sh` | **NEW** |
| `phase_4_postprocessing/4b_group_proposal/LAUNCH_C01i_decompose.sh` | **NEW** |
| `phase_4_postprocessing/4b_group_proposal/LAUNCH_C01i_b_multi_recomb.sh` | **NEW** |
| `phase_4_postprocessing/4b_group_proposal/LAUNCH_C01i_c_nested_comp.sh` | **NEW (stub)** — Python script missing; fails loudly |
| `phase_4_postprocessing/4b_group_proposal/LAUNCH_C01i_d_seal.sh` | **NEW** |
| `phase_4_postprocessing/4c_group_validation/LAUNCH_C01f_hypothesis.sh` | **NEW** |
| `phase_4_postprocessing/4e_final_classification/LAUNCH_C01d_scoring_pass2.sh` | **NEW** |
| `phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh` | **NEW** |
| `phase_4_postprocessing/4e_final_classification/figures/fig4e_*.R` (×3) + `test_figures.R` + `INSTALL.md` | **CARRIED FORWARD** from chat-15 separate bundle; now folded into main tarball |
| `patches_chat17/patch_orchestrator_paths.sh` | **NEW** — fixes run_phase4{,b}.sh stale refs |
| `patches_chat17/patch_STEP_C00_sv_prior_dir.sh` | **NEW** — makes STEP_C00 respect `$SV_PRIOR_DIR` |
| `patches_chat17/patch_rename_window_dt.sh` | **NEW** — `inv_like_dt` → `window_dt` across 10 files (cleans misleading name that mixed "inv_likeness score" with "per-window diagnostic table") |
| `patches_chat17/compile_engines.sh` | **NEW** — explicit build of all 5 C/C++ engine binaries |
| `LAUNCHER_TRIAGE.md` | **NEW** — decision matrix + DAG diagram |
| `FIXES_APPLIED.md` | Chat-17 section prepended |
| `SESSION_SUMMARY.md` | Rolled to end-of-chat-17 state (this file) |

## Static validation status (end of chat 17)

- 16 new launcher files (15 launchers + 1 chain runner) — all `bash -n`
  clean
- 2 patch scripts — both `bash -n` clean and tested against scratch
  copies of the real files
- Both patched orchestrators `bash -n` clean post-patch; all 6 LAUNCH_
  refs in `run_phase4.sh` point to real files in the current tree
- Patched `STEP_C00_build_sv_prior.R` passes `_rcheck.py`
- Both patches are idempotent (detect prior application, don't double-patch)
- No HPC access — all chat-16 + chat-17 work is still static-only

## Open findings at end of chat 17

**Chat 18 (first real HPC run on LANTA) — in order:**

1. 5 pre-flight checks from chat-17 handoff (unchanged, still critical)
2. Extract chat-17 tarball under `${BASE}/inversion_modules/`
3. Apply both patches (write `.bk_chat17` backups):
   `bash patches_chat17/patch_STEP_C00_sv_prior_dir.sh`
   `bash patches_chat17/patch_orchestrator_paths.sh`
4. Verify: `bash -n` on both orchestrators, `_rcheck.py` on C00
5. Run `submit_precomp_chain.sh` — ancestry + SV prior + C01a
6. Verify C01a stamps `localQ_*_K08` with non-NA values (Step 3 from
   chat-17 handoff)
7. BK schema extraction (still pending from chat 15 / chat 16)
8. GHSL v6 compute array + classify (phase 2e)
9. Phase 2d run_all array (existing launcher)
10. Phase 4a / 4b / 4c / 4d / 4e per DAG
11. Final `reg$results$integrity_check()` — expect `all_pass=TRUE`

**Chat 18 (low pri after HPC):**

- Decide `SV_PRIOR_DIR` canonical path (update config or move dir on
  disk — post-patch, either works)
- If needed: write/port `STEP_C01i_c_nested_composition.py` for phase
  4b.3, or wire `SKIP_4B3=1` into orchestrator
- Remove `reg$stats` deprecated alias (once HPC confirms no external
  deps)
- Legacy hash-tagged cache files cleanup (if any survived)
- Recursive sub-cluster detection (`STEP_C01i_e_subcluster.R`) if
  LG12/LG25 PCAs show sub-structure

**Deferred / low pri:** AL (sample_registry shadowing), `frequency.v2`
canonicalisation, terminology finalization (BF).

## Known paths and data conventions (post-chat-17)

- HPC base: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04`
- LANTA SLURM account: `lt200308`
- Results registry: `${BASE}/registries/data/results_registry/`
- Sample group: `SAMPLE_GROUP=all_226` (env default; launcher + R both
  read this)
- Filename convention: `C_gar_LG12.all_226.local_Q_*.tsv.gz`,
  `Q_K08.all_226.tsv.gz`, `all_226__vs__unrelated_81.tsv.gz`
- Manifest: `registries/data/results_registry/manifest.tsv` — append-
  only, one row per artifact
- **Launcher locations (new in chat 17):** per-phase subdirs
  (`phase_2_discovery/2c_precomp/LAUNCH_*.slurm` etc.), not under a
  single `launchers/` tree. See `LAUNCHER_TRIAGE.md` for the full list.

## Next milestones

1. **Chat 18** — first HPC run on LANTA. Pre-flights, apply patches,
   drop in launchers, run precomp chain, verify localQ merge, then
   proceed through the full DAG to `integrity_check()`.
2. **Chat 19+** — if LG12 / LG25 PCAs show sub-structure, build
   `STEP_C01i_e_subcluster.R`. Remove `reg$stats` alias. Continue to
   manuscript.
3. **Manuscript** — figrid figure composition. IF10+ target. The
   fourth registry + integrity check + full launcher coverage is
   a Data Availability paragraph reviewers care about.

## Posture going forward

Keep session summary + fixes-applied rolling. Keep chat-specific handoff
+ audit log (current chat only) at root. Archive per-chat handoff +
audit log from prior chats to `_archive/chat_history/`. Each new chat
reads root-level docs first, then `LAUNCHER_TRIAGE.md` if touching
phase-level code, then `registries/DATABASE_DESIGN.md` if touching
results_registry code.
