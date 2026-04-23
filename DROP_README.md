# DROP_README — pass 12: phase_qc_shelf + breakpoint_pipeline folded into phase_4

**Pass:** 12
**Date:** 2026-04-24
**Format:** full-tree tarball (not diff-only — this pass renames 6 directories, safer to ship the whole new tree)
**Scope:** Restructure `inversion_modules/phase_4_postprocessing/` from 5 sub-blocks (4a..4e) to 7 sub-blocks (4a..4g) by folding in `phase_qc_shelf/` and `breakpoint_pipeline/` as 4b and 4c respectively, and renumbering the existing 4b..4e → 4d..4g.

This replaces the architectural story from pass 11 ("phase_qc_shelf is a sibling of phase_4"). Pass 11's shipped code (bridge script) is preserved — it just lives at the new path `4b_qc_triage/scripts/bridge_from_phase4.sh`.

---

## ⚠️ Manual step: delete old directories

Because this is shipped as a full-tree replacement, after drag-drop you need to manually delete the old top-level module directories. GitHub Desktop should display this as a mass rename, but if it doesn't, nuke these two:

```bash
cd /mnt/c/Users/quent/Desktop/inversion-popgen-toolkit

rm -rf inversion_modules/phase_qc_shelf
rm -rf inversion_modules/breakpoint_pipeline
```

The old 4b..4e sub-block directories (`inversion_modules/phase_4_postprocessing/4b_group_proposal/`, etc.) are empty at their old names after the mass rename — GitHub Desktop should show those as moves, but confirm with:

```bash
ls inversion_modules/phase_4_postprocessing/
# Expected: 4a_existence_layers 4b_qc_triage 4c_breakpoint_refinement
#           4d_group_proposal 4e_group_validation 4f_group_dependent
#           4g_final_classification docs orchestrator patches schemas
#           specs tests README.md
```

If you see any of the OLD names (`4b_group_proposal`, `4c_group_validation`, `4d_group_dependent`, `4e_final_classification`) alongside the new names, remove the old ones:

```bash
rm -rf inversion_modules/phase_4_postprocessing/4b_group_proposal
rm -rf inversion_modules/phase_4_postprocessing/4c_group_validation
rm -rf inversion_modules/phase_4_postprocessing/4d_group_dependent
rm -rf inversion_modules/phase_4_postprocessing/4e_final_classification
```

(These should be empty / nonexistent after drag-drop; the commands are safety netting.)

---

## What changed

### Directory moves (6)

| Old path | New path |
|---|---|
| `inversion_modules/phase_qc_shelf/` | `inversion_modules/phase_4_postprocessing/4b_qc_triage/` |
| `inversion_modules/breakpoint_pipeline/` | `inversion_modules/phase_4_postprocessing/4c_breakpoint_refinement/` |
| `inversion_modules/phase_4_postprocessing/4b_group_proposal/` | `.../4d_group_proposal/` |
| `inversion_modules/phase_4_postprocessing/4c_group_validation/` | `.../4e_group_validation/` |
| `inversion_modules/phase_4_postprocessing/4d_group_dependent/` | `.../4f_group_dependent/` |
| `inversion_modules/phase_4_postprocessing/4e_final_classification/` | `.../4g_final_classification/` |

### Docs imported into 4c_breakpoint_refinement

Brought in from a separate upload (Quentin's laptop had additional material that wasn't in the in-repo copy):

- `docs/HANDOFF.md` — the 2026-04-20 session handoff
- `docs/SESSION_AUDIT.md` — same session's audit
- `docs/METHODOLOGY.md` — formal algorithmic description (LaTeX-friendly)
- `docs/MANUSCRIPT_CHUNKS.md` — draft methods/results/discussion paragraphs
- `docs/example_panel_D_output.png` — panel D replica output from STEP 07
- `_archive/architecture_classifier/` — STEP40 classifier + self-test (scaffolding, not core)
- `_archive/heatmap_ordering_upgrade/` — within_group_ordering.R + STEP36/41 (Module 1 scaffolding)
- `_archive/original_workflow_diagrams/` — pre-refactor ASCII diagrams
- `README.md` and `PIPELINE_DIAGRAM.md` at module root (were missing in the in-repo version)

The 8 core `.R` / `.py` / `.sh` scripts in the repo are **newer than the uploaded copies** and were NOT overwritten.

### Sed pass: 48+ files updated

Path references updated in shift order (4e→4g first, then 4d→4f, 4c→4e, 4b→4d; then new 4b=qc_triage, 4c=breakpoint_refinement):

- Orchestrators: `run_phase4.sh`, `run_phase4b.sh`, `LAUNCH_group_cheats.sh`
- All 7 LAUNCH_*.sh inside 4d_group_proposal (4× C01i launchers), 4e_group_validation (C01f), 4g_final_classification (characterize_classify)
- R internal `source()` cross-refs: `plot_sample_regime_dag.R`, `run_characterize.R`, `characterize_candidate.R`, `STEP_C01j_regime_compatibility_engine.R`, `STEP_C01l_local_structure_segments.R`, `STEP_C01m_distance_concordance.py`
- Tests: `test_phase4b_integration.py`, `test_dag_derive_R.R`, `test_gc_detector.R`, `test_figures.R`
- Cross-module: `phase_2_discovery/2e_ghsl/README.md`, `phase_5_followup/README.md`, `phase_6_secondary/README.md`, `phase_7_cargo/extra_plots/compute/_lib_final_labels.R`, `Modules/MODULE_2A_snp_discovery/README.md`, `Modules/MODULE_3_heterozygosity_roh/README.md`
- Tools: `tools/_bk_rename.py`, `tools/_code_field_check.py`, `inversion_modules/check_deployment_v10.py`
- Registry schemas: `BK_KEYS_EXPLAINED.md`, `fragment_distribution.schema.json`
- Root docs: `README.md`, `docs/MODULE_MAP.md`, `docs/ENGINES.md`, `docs/OPEN_TODOS.md`, `docs/STANDARDIZE_FOLDER.md`
- Installers: `tools/install_scripts/*` (old output-path scripts)
- Module root docs: `inversion_modules/README.md`, `CONFIG_ARCHITECTURE.md`, `PATH_WIRING.md`, `phase_4_postprocessing/README.md`

### New docs

- `docs/PHASE4_RENUMBER_PROPOSAL.md` — the design doc laying out the rationale, blast radius, risks, and execution order.
- `docs/MODULE_MAP.md` — updated phase table (removed standalone qc_shelf/bp_pipeline rows), updated sub-block table (4a..4g with all 7 roles), updated oddity #4 (bp pipeline folded in), updated quick-lookup rows.
- `inversion_modules/phase_4_postprocessing/4b_qc_triage/README.md` — title updated to reflect new location.
- `inversion_modules/phase_4_postprocessing/4c_breakpoint_refinement/README.md` — title updated.

### Preserved (no functional change)

- Pass 8 axis-5 wiring (env-gated by `V7_FINAL_DIR`, works via absolute path lookup — insensitive to folder rename). Verified `_axis5_final_label.R` source path in `compute_candidate_status.R` L958 intact.
- Pass 8 Gap 2 GDS block in `characterize_candidate.R::characterize_q5()`. Intact.
- Pass 9 phase_3 STEP_{A,B,D}NN_ naming. Independent of this pass.
- Pass 10 MODULE_5E archive. Independent.
- Pass 11 bridge script — moved from `phase_qc_shelf/scripts/` to `4b_qc_triage/scripts/`. `$here` self-relative resolution means it still works from its new location.
- Registry `METHOD_TAG="phase_qc_shelf"` in `STEP_Q10_register.sh` **kept as-is** — existing registry rows on HPC have this tag and renaming it would orphan them. This is an intentional exception to the path renames.

---

## Pre-ship verification

- [x] All 6 directory moves executed; old names don't exist in working tree
- [x] Round-1 + round-2 sed passes: 48+18 = 66 files updated
- [x] Post-sed grep for old refs: zero remaining (except intentional historical-context refs in MODULE_2A/3 and phase_5 READMEs, which were then separately updated to use new paths)
- [x] `bash -n` passes on 12 critical launchers (orchestrator, LAUNCHers, run_pipeline.sh, run_chrom.sh, run_all_28chrom.sh, bridge, STEP_Q10_register)
- [x] `python3 -m ast.parse` passes on 5 touched `.py` files (check_deployment, test_phase4b_integration, bp_pipeline_bridge, _bk_rename, _code_field_check)
- [x] JSON schema parses cleanly (fragment_distribution.schema.json)
- [x] MD fence balance: 13 major docs checked, all even
- [x] Pass-8 axis-5 wiring intact (V7_FINAL_DIR env var + _axis5_final_label.R source path)
- [x] Pass-11 bridge script still resolves its relative paths correctly (`$here/..` lookup unchanged)
- [x] No cohort-identity regressions

---

## Commit message

```
pass 12: fold phase_qc_shelf + breakpoint_pipeline into phase_4 (4b + 4c)

Restructure phase_4_postprocessing from 5 sub-blocks (4a..4e) to 7
(4a..4g) so the data-dependency sequence is enforced by folder
ordering rather than implied by convention.

Directory moves:
  phase_qc_shelf/              -> phase_4_postprocessing/4b_qc_triage/
  breakpoint_pipeline/         -> phase_4_postprocessing/4c_breakpoint_refinement/
  phase_4/4b_group_proposal/   -> phase_4/4d_group_proposal/
  phase_4/4c_group_validation/ -> phase_4/4e_group_validation/
  phase_4/4d_group_dependent/  -> phase_4/4f_group_dependent/
  phase_4/4e_final_classification/ -> phase_4/4g_final_classification/

Sed pass: 66 files updated across orchestrators, LAUNCH scripts, R
source-refs, Python tests, registry schemas, docs, and tools. Path
shift order was 4e→4g first to prevent collision.

Design rationale + blast radius: docs/PHASE4_RENUMBER_PROPOSAL.md
Updated: docs/MODULE_MAP.md (phase table, sub-block table, oddities)

Also imports docs and _archive into 4c_breakpoint_refinement from
Quentin's laptop copy that were missing from the in-repo version.

Preserves: pass 8 axis-5 wiring (env-gated), pass 8 Gap 2 GDS, pass 9
STEP_{A,B,D}NN_ phase 3, pass 10 MODULE_5E archive, pass 11 bridge
script (moved with the module). METHOD_TAG="phase_qc_shelf" inside
Q10 register kept as-is for registry row compatibility.
```

---

## Deployment note for HPC

On LANTA the layout must be updated to match:

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit

git pull  # after pushing pass 12
# Or if not using git, rsync from laptop post-drag-drop

# Verify the new layout
ls inversion_modules/phase_4_postprocessing/
# Should see: 4a..4g + docs/ orchestrator/ patches/ schemas/ specs/ tests/ README.md
```

Registry rows under `method="phase_qc_shelf"` remain valid (tag
unchanged). Old job output directories under `FOLLOWUP_DIR` are not
affected — they live outside `MODULE_DIR`.

---

## Phase 3-6 reorg progress (updated)

- ✅ **Pass 9**: phase_3_refine Layer-tagged rename
- ✅ **Pass 10**: phase_6 MODULE_5E archive
- 🟡 **Pass 11 (partial)**: phase_qc_shelf scale-out bridge script only
- ✅ **Pass 12 (this pass)**: phase_4 restructure — 4b_qc_triage + 4c_breakpoint_refinement folded in, 4b..4e shifted to 4d..4g
- ⏳ **Pass 11 tasks 2/3/4 (still pending)**: README rewrite of 4b_qc_triage for 2-mode use, MODULE_MAP polishing, 4e reader → 4g reader (q_qc_shelf_* keys in compute_candidate_status.R)

Next options:
- Finish pass 11 tasks 2/3/4 (now that 4g is the target for task 4, not 4e)
- phase_4d (now 4f) regroup — 14 flat files still need sub-grouping
- phase_5 reorg
- OPEN_TODOS items (T2/T3/T7/T12)

Pick whichever you want next.
