# Handoff prompt for chat 18

**Chat 18 priority:** first real HPC run on LANTA, now that chat 17 closed
the last non-HPC blocker (full launcher coverage + two path patches).

## TL;DR — do this in order

1. **Pre-flight checks** (5 of them, unchanged from chat 17 — still the
   critical first gate). See § "Pre-flight checks" below.
2. **Extract chat-17 deliverables** into the tree. See § "Step 0".
3. **Apply the two chat-17 patches** (idempotent, write `.bk_chat17`
   backups). See § "Step 1".
4. **Run the precomp chain** (ancestry → SV prior → C01a). See § "Step 2".
5. **Verify localQ merge** (the chat-16 silent-bug target). See § "Step 3".
6. **Phase 2e GHSL + Phase 2d candidate detection.** See § "Step 4".
7. **Phase 4 DAG** via orchestrator or manually with `--dependency=afterok:`.
   See § "Step 5".
8. **Final integrity check.** See § "Step 6".

**If any pre-flight fails, do not sbatch. Fix and re-check first.**

## Read at the start

1. This file.
2. `SESSION_SUMMARY.md` — current-state snapshot (post-chat-17).
3. `LAUNCHER_TRIAGE.md` — **new in chat 17.** Decision matrix for which
   scripts have launchers, DAG diagram of the full pipeline, known-issue
   callouts. Read before touching any phase-level code.
4. `registries/DATABASE_DESIGN.md` — chat-16 master design doc (unchanged).
5. `FIXES_APPLIED.md` — chat-17 section (CH17-*) explains the launcher
   audit and the two patches. Chat-16 section (CH16-*) covers the
   silent-bug fix that's now in production.

**Do NOT read** unless needed for specific debugging:
- `_archive/chat_history/` pre-chat-17 handoffs
- `_archive_superseded/`
- `HANDOFF_PROMPT_chat17_2026-04-17.md` — now archived, already executed

---

## What chat 17 delivered (non-HPC only; still needs live validation)

- **16 new launcher files.** 15 launchers covering every pipeline gap
  in phases 2c, 2e, 4a, 4b, 4c, 4e, plus 1 chain runner that sbatches
  ancestry → SV prior → C01a with `--dependency=afterok:`. All `bash -n`
  clean. See `LAUNCHER_TRIAGE.md` for the complete table.
- **2 patches for stale paths.** `patch_orchestrator_paths.sh` fixes
  `run_phase4{,b}.sh` references to an old directory layout. `patch_STEP_
  C00_sv_prior_dir.sh` makes the R script honour `$SV_PRIOR_DIR` instead
  of hardcoding its own output path. Both idempotent, both tested
  against scratch copies.
- **The 4b.3 Python script is still missing.** `STEP_C01i_c_nested_
  composition.py` referenced by `run_phase4b.sh` was never ported into
  the current tree. The chat-17 launcher is written correctly in shape
  and fails loudly when invoked without the script. Either add the
  script (chat 18 scope if needed), or set `SKIP_4B3=1` and accept
  that seal (4b.4) will write `unknown_no_engine_b` stubs.

### What chat 17 did NOT do (and why)

- **No HPC run.** Same constraint as chats 14–16 — no LANTA access in
  the sandbox. All launchers, patches, and orchestrators validated
  only by `bash -n` + `_rcheck.py` + scratch-copy tests.
- **No runtime fixes.** Everything fixed in chat 17 was a static
  surface (paths, coverage, config). Any R/Python runtime bug on LANTA
  is chat-18 territory.

---

## Pre-flight checks (same 5 as chat 17)

Still the critical gate. An HPC run that fails silently after 2 hours
costs a day. Pre-flights take 15 minutes.

### Pre-flight 0 — `registry_loader.R` sources cleanly on LANTA R

```bash
cd $BASE
source utils/pipeline_bridge.sh
Rscript -e '
  source(file.path(Sys.getenv("BASE"),
                    "registries/api/R/registry_loader.R"))
  cat("[preflight] registry_loader.R sources cleanly\n")
'
```

### Pre-flight 1 — `all_226` is registered

```bash
Rscript -e '
  source(file.path(Sys.getenv("BASE"), "utils/load_bridge.R"))
  stopifnot(reg$samples$has_group("all_226"))
  grp <- reg$samples$get_group("all_226")
  cat("[preflight] all_226 registered with", length(grp), "samples\n")
  stopifnot(length(grp) == 226)
'
```

### Pre-flight 2 — auto-migration of chat-15 `stats_cache/` (if present)

```bash
ls $BASE/registries/data/stats_cache/ 2>/dev/null && echo "legacy present"
Rscript -e '
  source(file.path(Sys.getenv("BASE"), "utils/load_bridge.R"))
  cat("[preflight] reg$status():\n"); reg$status()
'
ls $BASE/registries/data/stats_cache/ 2>/dev/null && \
  echo "[preflight] WARN: legacy dir still exists post-init"
```

### Pre-flight 3 — segment sanity guard

Same as chat 17 — optional smoke test; first real recombinant candidate
during LG12 run will exercise the guard naturally.

### Pre-flight 4 — old hash-tagged cache files

```bash
find $LOCAL_Q_DIR -name "*.N226_*.tsv.gz" 2>/dev/null | head
# If present, leave in place for now; only delete after new precompute works.
```

### Pre-flight 5 — `reg$status()` baseline

```bash
Rscript -e '
  source(file.path(Sys.getenv("BASE"), "utils/load_bridge.R"))
  reg$status()
'
```

Expected: `samples: 2+ groups`, `intervals: N candidates`,
`evidence: M candidates`, `results: 0 manifest rows` (or whatever
chat-15/16 left behind).

**If all 5 pass**, proceed to Step 0. **If any fails**, fix first.

---

## Step 0 — Extract chat-17 deliverables

```bash
cd $BASE/inversion_modules
tar -xf ~/inversion-popgen-toolkit_chat17_launchers_2026-04-18.tar

# Relocate per-phase launchers to their phase dirs:
mv phase_2c/*.slurm phase_2c/*.sh                  phase_2_discovery/2c_precomp/
mv phase_2e/*.slurm                                 phase_2_discovery/2e_ghsl/
mv phase_4a/*                                       phase_4_postprocessing/4a_existence_layers/
mv phase_4b/*                                       phase_4_postprocessing/4b_group_proposal/
mv phase_4c/*                                       phase_4_postprocessing/4c_group_validation/
mv phase_4e/*                                       phase_4_postprocessing/4e_final_classification/

# Chat-15 figures (carried in the same tarball, still unrun on real data)
mkdir -p phase_4_postprocessing/4e_final_classification/figures
mv figures_chat15/*.R figures_chat15/INSTALL.md \
   phase_4_postprocessing/4e_final_classification/figures/

# Patches go to their own subdir
mkdir -p patches_chat17
mv patches/*.sh patches_chat17/

# Optional: remove now-empty extract dirs
rmdir phase_2c phase_2e phase_4a phase_4b phase_4c phase_4e patches figures_chat15 2>/dev/null

# LAUNCHER_TRIAGE.md goes to repo root alongside SESSION_SUMMARY.md
mv LAUNCHER_TRIAGE.md ../
```

Sanity check that every launcher parses:

```bash
find phase_2_discovery phase_4_postprocessing -name "LAUNCH_*" -newer patches_chat17 -exec bash -n {} \; -print
```

## Step 0b — Run the chat-15 figure smoke test

The three figure scripts in `4e_final_classification/figures/` are at
level 1 (not run on real data). Before running them on the real
candidate table, exercise them against the synthetic data in
`test_figures.R`:

```bash
cd $BASE/inversion_modules/phase_4_postprocessing/4e_final_classification/figures/
Rscript test_figures.R
```

If this fails, the figure scripts will fail on real data too — fix first.
All three parse clean (`_rcheck.py`) so any failure is likely a runtime /
schema issue that the smoke test is there to surface. See `INSTALL.md`
in that dir for schema assumptions and status ladder.

## Step 1 — Apply the two patches

```bash
# Dry-run all first to see what they'll change
bash patches_chat17/patch_orchestrator_paths.sh --dry-run
bash patches_chat17/patch_STEP_C00_sv_prior_dir.sh --dry-run
bash patches_chat17/patch_rename_window_dt.sh --dry-run
bash patches_chat17/patch_rename_scientific_names.sh --dry-run
bash patches_chat17/patch_rename_scientific_names_v2.sh --dry-run

# Apply (order matters: rename_scientific_names v1 before v2)
bash patches_chat17/patch_orchestrator_paths.sh
bash patches_chat17/patch_STEP_C00_sv_prior_dir.sh
bash patches_chat17/patch_rename_window_dt.sh
bash patches_chat17/patch_rename_scientific_names.sh
bash patches_chat17/patch_rename_scientific_names_v2.sh

# Verify — all touched files should still parse clean
bash -n phase_4_postprocessing/orchestrator/run_phase4.sh
bash -n phase_4_postprocessing/orchestrator/run_phase4b.sh
python3 ../_rcheck.py phase_2_discovery/2c_precomp/STEP_C01a_precompute.R
for f in $(find . -name '*.R.bk_chat17*'); do
  orig=${f%.bk_chat17*}
  python3 ../_rcheck.py "$orig" || echo "RCHECK FAIL: $orig"
done
```

## Step 1b — Compile the C/C++ engines

The pipeline depends on 5 compiled binaries that live next to their
Makefiles (or, for hobs_windower, next to a raw .c file). **Every
launcher that needs them assumes they already exist at their configured
paths** (see `00_ancestry_config.sh` and `pipeline_bridge.sh`).

`pipeline_bridge.sh` has an auto-compile hook for instant_q + popstats +
hobs_windower that fires on first source, but it silences make errors
with `2>/dev/null` and does NOT cover the phase-5 engines. Run the
explicit helper once at first setup:

```bash
bash patches_chat17/compile_engines.sh
```

Expected output: `summary: 4 built, 0 failed, 0 skipped` (the phase-5
Makefile builds two binaries in one step, so the "built" count is 4
targets for 5 binaries).

On this sandbox (gcc 13.3.0) all five compile cleanly with only harmless
warnings (misleading-indentation in instant_q.cpp, strncpy-truncation in
rare_sfs_pairwise.c). If LANTA's toolchain is older and any of these
warnings become errors, the most likely culprit is `-std=c++17` on older
g++; pin the g++ version in the assembly conda env if needed.

Binaries produced:
- `unified_ancestry/src/instant_q` (used by ancestry precompute)
- `unified_ancestry/engines/fst_dxy/region_popstats` (FST/dXY tracks)
- `inversion_modules/phase_5_followup/engines/rare_sfs_pairwise`
- `inversion_modules/phase_5_followup/engines/export_q_residual_dosage`
- `unified_ancestry/engines/hobs_hwe/scripts/hobs_windower`

Note: stale `.c` copies of `rare_sfs_pairwise` and
`export_q_residual_dosage` also exist under `unified_ancestry/engines/
fst_dxy/` but that Makefile's comment header explicitly says they are
NOT built there anymore. The compile_engines helper flags this with a
warning. Clean up after HPC validation if desired (chat 19 low-pri).

### Decide: SV_PRIOR_DIR canonical path

After `patch_STEP_C00_sv_prior_dir.sh`, the R script honours `$SV_PRIOR_DIR`
with a back-compat default. Two options:

- **(a) Update config to export the script's hardcoded path** (simpler,
  no data move required):
  ```bash
  # In 00_inversion_config.sh, change line 181:
  export SV_PRIOR_DIR="${SV_PRIOR_DIR:-${MDS_DIR}/snake_regions_multiscale/sv_prior}"
  ```
- **(b) Leave config as-is and point the R script at the short path
  via env** (`$SV_PRIOR_DIR=${INVDIR}/03_sv_prior`). Requires ensuring
  any pre-existing SV prior files live at the same path.

Pick one per taste; the pipeline works for either.

## Step 2 — Precomp chain (the critical path)

```bash
cd phase_2_discovery/2c_precomp
./submit_precomp_chain.sh
# Default: K=8, runs ancestry → SV prior → C01a, each depends on prior.
# Monitor: squeue -u $USER
# Outputs land in standard locations (see submit_precomp_chain.sh header).
```

If you want to skip a step (e.g., ancestry already cached):

```bash
SKIP_ANCESTRY=1 ./submit_precomp_chain.sh    # skip step 1
SKIP_SV=1       ./submit_precomp_chain.sh    # skip step 2
```

## Step 3 — Verify localQ merge (the silent-bug target)

```bash
Rscript -e '
  x <- readRDS("<precomp_outdir>/precomp/C_gar_LG12.precomp.rds")
  cat("columns with localQ:\n")
  print(grep("^localQ_", names(x), value = TRUE))
  cat("n non-NA delta12_K08:", sum(!is.na(x$localQ_delta12_K08)),
      "/", nrow(x), "\n")
'
```

**Expected:** 3 columns (`localQ_delta12_K08`, `localQ_entropy_K08`,
`localQ_ena_K08`), non-NA count == full window count.

**If all NA:** the chat-17 C01a launcher printed a warning at startup
telling you which file it was looking for. Check that file exists.
If missing, the ancestry step didn't land files under the expected
`$SAMPLE_GROUP` filename convention.

## Step 4 — Phase 2e (GHSL) + Phase 2d (candidate detection)

```bash
cd ../2e_ghsl
sbatch LAUNCH_STEP_C04_ghsl_v6_compute.slurm            # array 1–28
# Wait for completion, then:
sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm          # all chroms

cd ../2d_candidate_detection
sbatch LAUNCH_run_all.slurm                              # existing, array 1–28
```

## Step 5 — Phase 4 DAG

Two options:

### Option A: orchestrator (recommended — dependencies handled)

```bash
cd ../../phase_4_postprocessing/orchestrator
bash run_phase4.sh                                # uses afterok chain
# OR with dry-run to see commands first:
bash run_phase4.sh --dry-run
```

Post-patch, this submits:

```
4a_score1  (C01d pass-1) ─┐
4a_bound   (C01g)         ├─→ 4b (run_phase4b.sh) ─→ 4c (C01f) ─→
                           │                                       4d (cheats) ─→ 4e ─→ characterize
```

### Option B: manual submission (for debugging)

```bash
# 4a (parallel)
JID_SCORE1=$(sbatch --parsable 4a_existence_layers/LAUNCH_C01d_scoring_pass1.sh)
JID_BOUND=$(sbatch  --parsable 4a_existence_layers/LAUNCH_C01g_boundary.sh)

# 4b (depends on 4a scoring)
JID_DECOMP=$(sbatch --parsable --dependency=afterok:$JID_SCORE1 \
                    4b_group_proposal/LAUNCH_C01i_decompose.sh)
JID_RECOMB=$(sbatch --parsable --dependency=afterok:$JID_DECOMP \
                    4b_group_proposal/LAUNCH_C01i_b_multi_recomb.sh)
JID_NESTED=$(sbatch --parsable --dependency=afterok:$JID_SCORE1 \
                    4b_group_proposal/LAUNCH_C01i_c_nested_comp.sh)
# ^^ 4b.3 will fail loudly — Python script missing. Either add it, or
#    skip by removing this line and letting seal write unknown_no_engine_b.
JID_SEAL=$(sbatch   --parsable \
                    --dependency=afterok:$JID_DECOMP:$JID_RECOMB:$JID_NESTED \
                    4b_group_proposal/LAUNCH_C01i_d_seal.sh)

# 4c, 4d, 4e in sequence with afterok — see run_phase4.sh for exact DAG
```

## Step 6 — Final integrity check

```bash
Rscript -e '
  source(file.path(Sys.getenv("BASE"), "utils/load_bridge.R"))
  ic <- reg$results$integrity_check()
  print(ic)
  cat("all_pass=", attr(ic, "all_pass"), "\n")
'
```

**Expected:** `all_pass=TRUE`. If not, the 6-check output tells you
exactly which FK/version/file-exists/orphan/sha256 check failed.

---

## Known carry-overs / low priority (unchanged from chat 17)

1. **Recursive sub-cluster detection** (`STEP_C01i_e_subcluster.R`) —
   naming convention in place, pipeline step not built. Tackle if
   LG12 / LG25 PCAs show clear sub-structure.
2. **`reg$stats` deprecated alias** — remove in chat 19 after HPC run
   confirms no external dependencies.
3. **Chat-15 BK schema extractions (BL/BM/BN/BO/BP/BQ)** — still pending
   HPC validation. Chat 16 + 17 didn't touch that code; no regression
   risk. Exercise during chat-18 phase-4b run.
4. **Legacy hash-tagged cache files** — delete once new K=8 precompute
   is done and everything is re-populated under the all_226 convention.
5. **`STEP_C01i_c_nested_composition.py`** — not in the tree. Either
   port from the phase4b_rewrite codebase (if it exists elsewhere) or
   skip by setting `SKIP_4B3=1` for this run.
6. **`SV_PRIOR_DIR` canonicalisation** — post-patch, decide which path
   wins (config or script hardcode) and update the other to match.

---

## Housekeeping at end of chat 18

```bash
mv HANDOFF_PROMPT_chat18_2026-04-18.md _archive/chat_history/
```

Create `HANDOFF_PROMPT_chat19_<date>.md` if needed. Update
`SESSION_SUMMARY.md` + `FIXES_APPLIED.md` in place. Tarball as
`inversion-popgen-toolkit_chat18_<topic>_<date>.tar`.

Expected chat-18 tarball contents: any HPC-discovered fixes, chat-18
handoff, audit log, rolled session summary.
