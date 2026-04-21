# Wiring notes — Stages 1–4 done, Stage 5 (script wiring) pending

**Date:** 2026-04-21
**Chat:** continuation from 2026-04-20 session
**Purpose:** record what was changed in the registry infrastructure so the
next step (wiring the 7 breakpoint scripts) can proceed without
rediscovering it.

---

## Headline

**To use the registry in any R script, source ONE file:**

```r
source("utils/registry_bridge.R")
```

That's the whole entry point. After sourcing, `reg`, `smap`, `get_Q()`,
`get_region_stats()`, and `BRIDGE_PATHS` are all available. Idempotent.

See `registries/HOW_TO_USE.md` for the canonical user-facing docs.

---

## What was touched this turn

### Stage 1 — sample_registry.R canonical location

| File | Action |
|---|---|
| `registries/api/R/sample_registry.R` | **NEW** (copied from `utils/`, 309 lines, identical contents) |
| `utils/sample_registry.R` | **REPLACED** with a compat shim (83 lines) that sources the canonical copy |

### Stage 2 — registry backend for the bridge

| File | Change |
|---|---|
| `registries/api/R/registry_loader.R` | Added 4 top-level method aliases (`list_groups`, `get_master`, `count_groups`, `get_groups_for_candidate`) so flat-API callers work unchanged against the new nested `reg` |
| `utils/load_bridge.R` (pre-rename) | Rewrote STEP 4 to source the full `registries/api/R/registry_loader.R` instead of the flat `utils/sample_registry.R` |
| Same | Added `BRIDGE_PATHS$REGISTRIES_ROOT` with auto-detection |
| Same | Fixed final `rm()` to use `intersect()` (safer with conditional temporaries) |

### Stage 3 — dispatcher symmetric upgrade

| File | Change |
|---|---|
| `unified_ancestry/dispatchers/region_stats_dispatcher.R` | Dispatcher's own registry-acquisition path now tries the full loader first, flat-loader fallback second. Matches `load_bridge.R`'s policy. |

### Stage 4 — rename for discoverability (chat-18)

| File | Action |
|---|---|
| `utils/registry_bridge.R` | **NEW** (canonical bridge, 455 lines — same code as former `load_bridge.R` with an updated header docblock explaining what it is and what it provides) |
| `utils/load_bridge.R` | **REPLACED** with a 3-line compat shim that forwards to `utils/registry_bridge.R`. Existing callers keep working; one-line deprecation message suppressible via `QUIET_LOAD_BRIDGE_SHIM=1` |
| `utils/pipeline_bridge.sh` | Exports `REGISTRY_BRIDGE` (preferred, post-chat-18) alongside legacy `LOAD_BRIDGE` alias. Search order prefers new filename, accepts old. Comments updated. |
| `registries/HOW_TO_USE.md` | **NEW** user-facing docs. The one-line entry point is visible the moment someone `ls` in `registries/` |

### Why rename instead of move

User feedback: "Like now its in the folder utils/ at the root. I think its good. but its not explicit that like its registry." Fix: keep the file in `utils/` (same folder as other R libraries — `sample_map.R`, `lib_ghsl_panel.R`, `theme_systems_plate.R`) but rename it to make the purpose visible. No folder reorganization, no move of dependencies, zero risk of breaking references.

---

## Backward-compatibility invariants (still hold)

**Two alias layers:**

1. **Flat API on `reg`** — the top-level `reg$has_group()`, `reg$add_group()`, `reg$get_group()`, `reg$list_groups()`, `reg$get_master()`, `reg$count_groups()`, `reg$get_groups_for_candidate()`, `reg$register_candidate()`, `reg$add_evidence()`, `reg$get_evidence()`, `reg$has_evidence()` all continue to work on the new nested `reg`. They forward to `reg$samples$*` or `reg$evidence$*`. Used by 7 existing consumer scripts.

2. **Old filename paths** — `source("utils/load_bridge.R")` still works (shim) and `source("utils/sample_registry.R")` still works (shim). The env var `$LOAD_BRIDGE` is still set (to the same value as `$REGISTRY_BRIDGE`). Used by ~22 existing references across inversion_modules/, unified_ancestry/, and Modules/.

**No existing script or launcher needs source-path changes.** Everything keeps running. New scripts should prefer the new names.

---

## What was NOT touched

- **`unified_ancestry/registries/build_registries.py`** — remains named misleadingly. It's not a registry in the DATABASE_DESIGN sense; it's a precomputed TSV index consumed by the `instant_q` / `region_stats` bash CLIs. Rename to `unified_ancestry/job_index/` is a separate task that touches bash paths.
- **No phase-4 script edits.** The 5 phase-4 scripts that use flat API continue to work via aliases.
- **No SLURM launcher edits.** They all use `$LOAD_BRIDGE`, which is now an alias for `$REGISTRY_BRIDGE`. Both env vars are set after `source utils/pipeline_bridge.sh`.

---

## Files changed this tarball

```
staging/
├── utils/
│   ├── registry_bridge.R          NEW (canonical bridge, 455 lines)
│   ├── load_bridge.R              SHIM (forwards to registry_bridge.R)
│   ├── sample_registry.R          SHIM (forwards to registries/api/R/sample_registry.R)
│   └── pipeline_bridge.sh         UPDATED (exports REGISTRY_BRIDGE + LOAD_BRIDGE alias)
├── registries/
│   ├── HOW_TO_USE.md              NEW (user-facing entry-point docs)
│   └── api/R/
│       ├── registry_loader.R      UPDATED (+4 top-level aliases)
│       └── sample_registry.R      NEW (canonical 309-line copy)
├── unified_ancestry/dispatchers/
│   └── region_stats_dispatcher.R  UPDATED (full-loader-first, flat fallback)
└── WIRING_NOTES_STAGE_1_2_3_4.md  this file
```

Deploy on LANTA with:

```bash
tar -xzf wiring_stages_1_2_3_4.tar.gz -C $BASE/
```

Then sanity-check:

```bash
source utils/pipeline_bridge.sh
echo "REGISTRY_BRIDGE: $REGISTRY_BRIDGE"    # should point at utils/registry_bridge.R
echo "LOAD_BRIDGE:     $LOAD_BRIDGE"        # should be the same value

Rscript -e '
  source(Sys.getenv("REGISTRY_BRIDGE"))
  stopifnot(reg$samples$has_group("all_226"))
  cat("[preflight OK] all_226:", length(reg$samples$get_group("all_226")), "samples\n")
  cat("[preflight OK] reg$intervals:", reg$intervals$count_candidates(), "candidates\n")
  cat("[preflight OK] reg$results:", nrow(reg$results$list_cached()), "manifest rows\n")
'
```

If that prints three `[preflight OK]` lines, the infrastructure is
correctly wired and ready for Stage 5 (the 7 breakpoint scripts).

---

## Catalog of still-to-update callers (cosmetic only — they all work via shims)

Files that still reference the old name `load_bridge.R` — everything runs
correctly today via the shim. Batch-update in a later chat when touching
these files for other reasons:

**R scripts (6):**
- `inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R`
- `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`
- `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R`
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/engine_b_smoke_test.R`
- `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`
- `inversion_modules/phase_4_postprocessing/4c_group_validation/_PRISTINE_v9.3.4.R`
- `inversion_modules/phase_4_postprocessing/4d_group_dependent/cheat6_ancestry_jackknife_v934.R`
- `Modules/MODULE_2B_structure/scripts/STEP_A04_best_seed_by_K.R`
- `unified_ancestry/dispatchers/region_stats_dispatcher.R` (self — only in search fallback list)
- `unified_ancestry/wrappers/instant_q.R`

**Bash / SLURM (7):**
- `inversion_modules/00_inversion_config.sh`
- `inversion_modules/phase_2_discovery/2c_precomp/submit_precomp_chain.sh`
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/LAUNCH_C01g_boundary.sh`
- `inversion_modules/phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh`
- `Modules/MODULE_2B_structure/00_module2b_config.sh`
- `unified_ancestry/00_ancestry_config.sh`
- `unified_ancestry/run_ancestry.sh`

**Docs (4):**
- `FIXES_APPLIED.md`
- `HANDOFF_PROMPT_chat18_2026-04-18.md`
- `inversion_modules/CONFIG_ARCHITECTURE.md`
- `inversion_modules/phase_4_postprocessing/patches/ENGINE_B_SMOKE_TEST_INSERTS.md`
- `inversion_modules/phase_4_postprocessing/docs/ADDITIONS_v10_1_1.md`

Grep-and-sed command to batch-update all references when ready:

```bash
# Preview
grep -rln 'utils/load_bridge\.R' \
  --include='*.R' --include='*.sh' --include='*.slurm' --include='*.md' \
  --exclude-dir=_archive --exclude-dir=_archive_superseded .

# Apply (AFTER preview)
grep -rl 'utils/load_bridge\.R' \
  --include='*.R' --include='*.sh' --include='*.slurm' --include='*.md' \
  --exclude-dir=_archive --exclude-dir=_archive_superseded . | \
  xargs sed -i 's|utils/load_bridge\.R|utils/registry_bridge.R|g'
```

---

## What the next step (Stage 5 — script wiring) will use

The 7 breakpoint-pipeline scripts will each start with:

```r
source("utils/registry_bridge.R")
```

After that line, every script is self-contained: it can find its
candidate coordinates via `reg$intervals$get_candidate(cid)`, its
karyotype samples via `reg$samples$get_groups_for_candidate(cid)`,
write results via `reg$results$put_*()` with provenance + sha256,
and check compute-if-missing via `reg$results$ask_what_for_candidate(cid)`.

No chained-source requirements, no hardcoded paths, no order-dependency
between scripts. Each script works standalone.

---

## Validation done in the sandbox

All 6 modified R files + 1 bash file pass:

- Brace/paren/bracket balance check (R-aware tokenizer: respects `"..."`,
  `'...'`, `` `...` `` strings with `\\` escapes, skips comments)
- `bash -n` syntax check (for pipeline_bridge.sh)

No R runtime tests — no R interpreter in the sandbox. Runtime validation
deferred to LANTA. The preflight snippet above is the suggested check.
