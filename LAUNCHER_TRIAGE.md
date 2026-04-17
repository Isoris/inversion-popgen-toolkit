# SLURM launcher coverage — triage

**Scope:** `inversion_modules/` + `unified_ancestry/`. Upstream `Modules/` (1–4) is
already fully launcher-covered and excluded here.

**Convention (style match to archive):** `set -euo pipefail` → `source ~/.bashrc` →
`mamba activate assembly` → source `00_inversion_config.sh` with `set -a` →
source `utils/pipeline_bridge.sh` → `inv_init_dirs` → banner → `${RSCRIPT_BIN}
script.R args` → footer. `-A lt200308 -p compute`. SAMPLE_GROUP / RESULTS_REGISTRY_DIR
propagated via bridge (chat-16 wiring).

---

## Phase 1 — Inputs
**Status: covered.** `launchers/LAUNCH_STEP04…07_*.slurm` exist for the ANGSD
SAF/SFS/beagle chain. `runners/run_5A1_discovery_inputs.sh` orchestrates.
**No new launchers needed.**

## Phase 2a — local PCA
**Status: covered.** `LAUNCH_A01_beagle_to_dosage`, `LAUNCH_A02_local_pca_legacy`,
`LAUNCH_A03_dense_registry_stage1/2`. The STEP_A02/A03 scripts are called by
these. **No new launchers needed.**

## Phase 2b — MDS
**Status: covered.** `LAUNCH_B01_mds_legacy`, `LAUNCH_B01_mds_stage1/2`.
**No new launchers needed.**

## Phase 2c — precomp
**Status: was missing.** Just added (this chat):

| Launcher                                  | Calls                         | Status    |
|-------------------------------------------|-------------------------------|-----------|
| `LAUNCH_STEP_C00_sv_prior.slurm`          | `STEP_C00_build_sv_prior.R`   | ✅ new    |
| `LAUNCH_STEP_C01a_precompute.slurm`       | `STEP_C01a_precompute.R`      | ✅ new    |
| `LAUNCH_PHASE_01C_block_detect.slurm`     | `PHASE_01C_block_detect.R`    | ✅ new    |
| `LAUNCH_STEP_C01b_1_seeded_regions.slurm` | `STEP_C01b_1_seeded_regions.R`| ✅ new    |
| `submit_precomp_chain.sh`                 | chains ancestry→C00→C01a      | ✅ new    |

Ordering: `instant_q_precompute → C00 → C01a → 01C_block_detect → C01b_1_seeded_regions`.

Diag scripts (`diags/STEP_C01a_*`) are post-hoc diagnostics — run interactively,
**no launchers needed**.

## Phase 2d — candidate detection
**Status: covered.** `LAUNCH_run_all.slurm` is the omnibus per-chrom array launcher
(1–28) that calls `run_all.R`, which in turn calls all STEP_D01–D17 scripts.
`LAUNCH_peel.slurm` / `LAUNCH_test_LG01.slurm` for specific sub-modes. Individual
STEP_D scripts are pipeline components, not standalone entry points — no per-step
launchers needed. **No new launchers needed.**

## Phase 2e — GHSL v6
**Status: was missing.** Two scripts, both standalone per-chrom:

| Launcher                               | Calls                          | Status  |
|----------------------------------------|--------------------------------|---------|
| `LAUNCH_STEP_C04_ghsl_v6_compute.slurm`| `STEP_C04_snake3_ghsl_v6.R`    | ✅ new  |
| `LAUNCH_STEP_C04b_ghsl_v6_classify.slurm` | `STEP_C04b_snake3_ghsl_classify.R` | ✅ new |

STEP_C04 is the heavy engine (run ONCE per chromosome, ~77M variants). STEP_C04b
is the light classifier (30 s per chrom, iterate while tuning). C04b depends on C04.

## Phase 3 — breakpoint refinement
**Status: covered.** `run_breakpoint_validation.sh` is the wrapper. The individual
`01–06_*` scripts are steps inside that wrapper. `annotate_population_confidence.sh`
is a utility. **No new launchers needed** unless the individual steps get promoted
to separate SLURM jobs, which isn't how the orchestrator reads.

## Phase 4a — existence layers
**Status: missing (referenced by `run_phase4.sh` orchestrator).** Two launchers
explicitly named by the orchestrator:

| Launcher                           | Calls                                                       | Status  |
|------------------------------------|-------------------------------------------------------------|---------|
| `LAUNCH_C01d_scoring_pass1.sh`     | `STEP_C01d_candidate_scoring_wired_25_v934_registry.R`      | ✅ new  |
| `LAUNCH_C01g_boundary.sh`          | `STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` | ✅ new  |

Other 4a scripts (`C01e`, `C01j`, `C01l`, `C01m`) are internal modules called by
C01d — **no individual launchers needed** based on the orchestrator DAG.
`C01e_candidate_figures` is a post-hoc figure script — runs interactively.

## Phase 4b — group proposal
**Status: missing (referenced by `run_phase4b.sh` orchestrator).** Four jobs,
one per sub-phase:

| Launcher                      | Calls                          | Status  |
|-------------------------------|--------------------------------|---------|
| `LAUNCH_C01i_decompose.sh`    | `STEP_C01i_decompose.R`        | ✅ new  |
| `LAUNCH_C01i_b_multi_recomb.sh` | `STEP_C01i_b_multi_recomb.R` | ✅ new  |
| `LAUNCH_C01i_c_nested_comp.sh` | `STEP_C01i_c_nested_composition.py` (**not in tree**) | ⚠ stub only |
| `LAUNCH_C01i_d_seal.sh`       | `STEP_C01i_d_seal.R`           | ✅ new  |

⚠ `STEP_C01i_c_nested_composition.py` is referenced by `run_phase4b.sh` but not
present in the current tree. Launcher written as a stub that will fail loudly if
the script is absent — this is the correct behaviour (pipeline should not silently
skip a sub-phase that was supposed to run).

`lib_*.R`, `gene_conversion_detector.R`, `engine_b_smoke_test.R`,
`plot_sample_regime_dag.R` are helpers / tests / post-hoc plots — no launchers.

## Phase 4c — group validation
**Status: missing (referenced by `run_phase4.sh`).**

| Launcher                         | Calls                           | Status  |
|----------------------------------|---------------------------------|---------|
| `LAUNCH_C01f_hypothesis.sh`      | `STEP_C01f_hypothesis_tests.R`  | ✅ new  |

`group_validation_gate.R` is a helper (no `commandArgs`) — not standalone.

## Phase 4d — group-dependent cheats
**Status: mostly covered.** `orchestrator/LAUNCH_group_cheats.sh` exists. The
individual `cheat*.R` scripts are called by it. `SLURM_A03b_population_regenotype.sh`
is a genuine SLURM worker (sbatch-able). **No new launchers needed** — the
orchestrator's `LAUNCH_group_cheats.sh` is the entry point.

## Phase 4e — final classification
**Status: missing (referenced by `run_phase4.sh`).**

| Launcher                               | Calls                          | Status  |
|----------------------------------------|--------------------------------|---------|
| `LAUNCH_C01d_scoring_pass2.sh`         | `STEP_C01d_*.R` (same script, pass-2 mode) | ✅ new |
| `LAUNCH_characterize_classify.sh`      | `run_characterize.R` + `compute_candidate_status.R` | ✅ new |

`STEP_C01k_annotated_simmat.R` is a sim-matrix annotator — called internally.
`characterize_candidate.R` is the library; `run_characterize.R` is the driver.

**Figures subdir** (`4e_final_classification/figures/`): contains three
chat-15 publication-figure R scripts (`fig4e_gene_count_by_biotype.R`,
`fig4e_class_along_inversion.R`, `fig4e_cumulative_burden_per_group.R`)
plus `test_figures.R` smoke harness and schema-assumptions doc
(`INSTALL.md`). **These are interactive figure scripts, not SLURM jobs**
— run directly with `Rscript` on a compute node or login node. No
launchers needed. Carried forward from the chat-15 separate bundle;
folded into the main chat-17 tarball so chat 18 has everything in one
drop. Smoke test has never run on real data — first chat-18 action.

## Phase 5 — followup
**Status: already covered for experimental codebase; two new standalone
launchers exist for followup compute:**

- `legacy_followup/STEP14–18_*.slurm` — the v8.5 experimental followup suite.
- `launchers/LAUNCH_q_residual_dosage.slurm`, `LAUNCH_rare_sfs_pairwise.slurm` —
  pop-gen region stats.

**No new launchers needed** unless the candidate-followup gets a new driver.

## Phase 6 — secondary modules (LD / FST / HObs)
**Status: covered.** Each secondary module has its own SLURM files under
`MODULE_5C/5D/5E/`. **No new launchers needed.**

## unified_ancestry
**Status: covered.** `LAUNCH_instant_q_precompute.slurm` (chat-16 rewrite),
`LAUNCH_region_popstats.slurm`, `LAUNCH_rare_sfs_pairwise.slurm`,
`LAUNCH_q_residual_dosage.slurm`, `engines/hobs_hwe/slurm/LAUNCH_hobs_hwe.slurm`.
**No new launchers needed.**

---

## Summary

**New launchers written this session: 16 files (15 launchers + 1 chain runner).**

| Phase | Count | Files |
|-------|-------|-------|
| 2c    | 4 + chain | `LAUNCH_STEP_C00_sv_prior.slurm`, `LAUNCH_STEP_C01a_precompute.slurm`, `LAUNCH_PHASE_01C_block_detect.slurm`, `LAUNCH_STEP_C01b_1_seeded_regions.slurm`, `submit_precomp_chain.sh` |
| 2e    | 2 | `LAUNCH_STEP_C04_ghsl_v6_compute.slurm` (array 1–28), `LAUNCH_STEP_C04b_ghsl_v6_classify.slurm` |
| 4a    | 2 | `LAUNCH_C01d_scoring_pass1.sh`, `LAUNCH_C01g_boundary.sh` |
| 4b    | 4 | `LAUNCH_C01i_decompose.sh`, `LAUNCH_C01i_b_multi_recomb.sh`, `LAUNCH_C01i_c_nested_comp.sh` (stub), `LAUNCH_C01i_d_seal.sh` |
| 4c    | 1 | `LAUNCH_C01f_hypothesis.sh` |
| 4e    | 2 | `LAUNCH_C01d_scoring_pass2.sh`, `LAUNCH_characterize_classify.sh` |

**All 16 files pass `bash -n`.** None run on LANTA yet — static-only validation.

### Pipeline DAG the launchers cover

```
    ┌──── unified_ancestry/LAUNCH_instant_q_precompute.slurm  (existing)
    │             │
    │             ↓
    │     phase_2c/LAUNCH_STEP_C00_sv_prior.slurm
    │             │
    └─────┬───────┘
          ↓
    phase_2c/LAUNCH_STEP_C01a_precompute.slurm
          │
          ↓
    phase_2c/LAUNCH_PHASE_01C_block_detect.slurm
          │
          ↓
    phase_2c/LAUNCH_STEP_C01b_1_seeded_regions.slurm
          │          (plus GHSL track, parallel)
          │
          │     phase_2e/LAUNCH_STEP_C04_ghsl_v6_compute.slurm  (array)
          │             │
          │             ↓
          │     phase_2e/LAUNCH_STEP_C04b_ghsl_v6_classify.slurm
          │
          ↓
    phase_2d/LAUNCH_run_all.slurm  (existing, array 1–28)
          │
          ↓
    phase_4a/LAUNCH_C01d_scoring_pass1.sh   ─┐  (both run in parallel — 4a)
    phase_4a/LAUNCH_C01g_boundary.sh         ─┘
          │
          ↓  (4b DAG — see run_phase4b.sh)
    phase_4b/LAUNCH_C01i_decompose.sh  ──→  LAUNCH_C01i_b_multi_recomb.sh ─┐
                         │                                                 ├─→ LAUNCH_C01i_d_seal.sh
    phase_4b/LAUNCH_C01i_c_nested_comp.sh (parallel with decompose) ──────┘
          │
          ↓
    phase_4c/LAUNCH_C01f_hypothesis.sh
          │
          ↓
    orchestrator/LAUNCH_group_cheats.sh  (existing)
          │
          ↓
    phase_4e/LAUNCH_C01d_scoring_pass2.sh
          │
          ↓
    phase_4e/LAUNCH_characterize_classify.sh
          │
          ↓
    reg$results$integrity_check()  →  figrid manuscript figures
```

### Orchestrator compatibility

The existing `orchestrator/run_phase4.sh` names 6 launchers — this session
provides all 6 (4 new + the 2 that already existed: `LAUNCH_group_cheats.sh`
and the 4b chain via `run_phase4b.sh`). The new launchers accept the
orchestrator's `--chroms` and `--array=...` arguments as no-ops where
appropriate (orchestrator passes them through blindly).

### Known path inconsistencies requiring user attention

1. **STEP_C00 output path vs `SV_PRIOR_DIR` config.** `00_inversion_config.sh`
   exports `SV_PRIOR_DIR=${INVDIR}/03_sv_prior`, but `STEP_C00_build_sv_prior.R`
   hardcodes output to `${INVDIR}/06_mds_candidates/snake_regions_multiscale/sv_prior`.
   Launchers use the script's actual path. Fix one or the other.

2. **`run_phase4.sh` path mismatches.** Script references
   `${MODULE_DIR}/inversion-popgen-toolkit/launchers/*` but the current tree
   uses `inversion_modules/phase_4_postprocessing/{4a,4b,4c,4e}/`. The new
   launchers live in the current tree locations. When running, either:
   - Patch `run_phase4.sh` to point at `inversion_modules/phase_4_postprocessing/4?/`, or
   - `sbatch` the launchers manually in DAG order with `--dependency=afterok:`.

3. **`run_phase4b.sh` path mismatch.** References `phase4b_rewrite/R/`; current
   tree is `4b_group_proposal/`. The new 4b launchers use the current path.
   `run_phase4b.sh` uses `--wrap=` style (no external launcher), so the two
   patterns can coexist — use whichever makes sense.

4. **`STEP_C01i_c_nested_composition.py` missing.** Referenced by
   `run_phase4b.sh`, not present in the tree. The launcher is written
   correctly and will fail loudly with a clear message when the script is
   absent. This is the correct behaviour for a missing DAG node — pipelines
   should not silently skip a sub-phase.

### Pre-flight before first sbatch on LANTA

Run the 5 checks from `HANDOFF_PROMPT_chat17_2026-04-17.md` first:
0. `registry_loader.R` sources cleanly in LANTA R
1. `all_226` is registered in sample_registry
2. Auto-migration of chat-15 `stats_cache/` → `results_registry/` happened
3. No hash-tagged leftover cache files
4. `reg$status()` shows expected baseline

**Then** start with `phase_2c/submit_precomp_chain.sh` — that chains
ancestry → SV prior → C01a precomp with `--dependency=afterok:` so a silent
upstream failure blocks the rest instead of producing junk downstream.
