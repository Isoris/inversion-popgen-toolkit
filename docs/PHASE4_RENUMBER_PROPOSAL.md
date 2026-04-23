# Phase 4 renumber proposal — insert qc_triage + breakpoint_refinement

**Date:** 2026-04-24
**Status:** proposal — to be executed in pass 12 (or later)
**Rationale chat:** `docs/HANDOFF_2026-04-24_pass11_qc_shelf_wiring.md` lays out why the current sibling layout for `phase_qc_shelf` and `breakpoint_pipeline` mis-represents the real pipeline dependencies.

---

## Why renumber

Current layout has both `phase_qc_shelf` and `breakpoint_pipeline`
sitting at the same level as `phase_4_postprocessing` — as though all
three run in parallel. The real data-dependency chain is **sequential**:

1. `phase_4a` (catalog-birth) produces a candidate list with coarse
   start_mb / end_mb shelf coordinates from local-PCA detection.
2. `phase_qc_shelf` triages each candidate: is the shelf real, or is
   it a data-quality artifact? Writes `q_qc_shelf_flag` (clean /
   low_snp / high_uncertain / coverage_artifact / messy).
3. `breakpoint_pipeline` refines the coarse interval to bp-resolution:
   uses dosage signal + per-carrier ancestral fragment distribution
   to produce `final_left_bp`, `final_right_bp`, CI bounds.
4. Only then do group proposal (C01i), group validation (C01f), group-
   dependent analyses (cheats), and final classification run — all on
   the refined interval + QC flag.

If step 2 or 3 is skipped, step 4's Fisher tests and group dosage
assignments operate on the wrong region. A 200 kb breakpoint error
can mean samples near the shelf edge are assigned to the wrong group,
which cascades through all downstream statistics.

Currently this sequence is **convention**, not **enforcement**. Moving
`phase_qc_shelf` and `breakpoint_pipeline` inside `phase_4_postprocessing`
as ordered subfolders makes the dependency visible in the file tree
and enforceable by the orchestrator.

---

## Proposed layout

```
phase_4_postprocessing/
├── 4a_existence_layers/           (unchanged: catalog-birth from local-PCA)
├── 4b_qc_triage/                  (moved from inversion_modules/phase_qc_shelf/)
├── 4c_breakpoint_refinement/      (moved from inversion_modules/breakpoint_pipeline/)
├── 4d_group_proposal/             (was 4b_group_proposal)
├── 4e_group_validation/           (was 4c_group_validation)
├── 4f_group_dependent/            (was 4d_group_dependent)
├── 4g_final_classification/       (was 4e_final_classification)
├── docs/
├── orchestrator/
├── patches/
├── schemas/
├── specs/
└── tests/
```

No subfolder is removed; the existing 4b→4g shift renumbers but
preserves structure.

---

## Integration principles (preserved)

1. **Soft gates, no hard blocks.** A candidate flagged "messy" by
   4b_qc_triage still passes through 4c..4g. The flag rides as a
   column on the candidate table; downstream scorers may deprioritize
   but never drop.

2. **Fallback to coarse bp.** If 4c_breakpoint_refinement fails (too
   few informative markers, insufficient homozygote carriers, low
   MAF), candidates fall back to the coarse phase_4a interval.
   Tracked via a `bp_refinement_status` column.

3. **Idempotent writes.** Both moved modules already write their own
   registries. Moving them into phase_4 doesn't change the write
   semantics — just the folder location and the orchestrator
   references.

4. **Self-relative paths preserved.** Both modules use
   `$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)` for their own
   internal refs, so internal launcher paths survive the move. Only
   external references (cross-module cites) need updating.

---

## Blast radius (pass 12 scope)

### Directory moves (3 operations)

```bash
cd inversion_modules

mv phase_qc_shelf        phase_4_postprocessing/4b_qc_triage
mv breakpoint_pipeline   phase_4_postprocessing/4c_breakpoint_refinement

cd phase_4_postprocessing
mv 4b_group_proposal        4d_group_proposal
mv 4c_group_validation      4e_group_validation
mv 4d_group_dependent       4f_group_dependent
mv 4e_final_classification  4g_final_classification
```

### Files needing sed updates

**Critical path (orchestrator + launchers):**
- `phase_4_postprocessing/orchestrator/run_phase4.sh` — 10+ refs to 4b/4c/4d/4e
- `phase_4_postprocessing/orchestrator/run_phase4b.sh` — refs to 4b_group_proposal (new name: 4d_group_proposal)
- `phase_4_postprocessing/orchestrator/LAUNCH_group_cheats.sh` — refs to 4d_group_dependent
- `phase_4_postprocessing/4b_group_proposal/LAUNCH_*.sh` × 4 — self-refs to 4b
- `phase_4_postprocessing/4c_group_validation/LAUNCH_C01f_hypothesis.sh`
- `phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh`
- `phase_4_postprocessing/4e_final_classification/_axis5_final_label.R`

**Tests:**
- `phase_4_postprocessing/tests/test_dag_derive_R.R`
- `phase_4_postprocessing/tests/test_gc_detector.R`
- `phase_4_postprocessing/4e_final_classification/figures/test_figures.R`
- `phase_4_postprocessing/4e_final_classification/figures/INSTALL.md`

**Cross-module consumers:**
- `phase_7_cargo/extra_plots/compute/_lib_final_labels.R`
- `phase_2_discovery/2e_ghsl/README.md`
- `phase_6_secondary/README.md`
- `registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md`

**Docs:**
- `docs/MODULE_MAP.md` — restructures phase_4 sub-block table + several
  quick-lookup rows
- `docs/HANDOFF_2026-04-24_pass11_qc_shelf_wiring.md` — historical
  reference to old paths (update with a "superseded by pass 12" note)

**Tools:**
- `tools/_bk_rename.py`
- `tools/_code_field_check.py`
- `inversion_modules/check_deployment_v10.py`

**Root-level docs:**
- `inversion_modules/README.md` — tree diagram
- `inversion_modules/CONFIG_ARCHITECTURE.md` — module list
- `inversion_modules/phase_4_postprocessing/README.md`

### Sed patterns

Apply these, in order, across all files in the blast list:

```bash
# Directory-name renames — apply to path references in configs/scripts/docs.
# Order matters: do the INSERTIONS first (shift 4b..4e → 4d..4g), then the
# new block names (qc_shelf, breakpoint_pipeline). Otherwise the shifted
# names collide with the source names.

# 1. Shift 4e → 4g (do this first — nothing else lands on 4g)
sed -i 's|4e_final_classification|4g_final_classification|g' $FILES
# 2. Shift 4d → 4f
sed -i 's|4d_group_dependent|4f_group_dependent|g' $FILES
# 3. Shift 4c → 4e
sed -i 's|4c_group_validation|4e_group_validation|g' $FILES
# 4. Shift 4b → 4d
sed -i 's|4b_group_proposal|4d_group_proposal|g' $FILES
# 5. Insert new 4b and 4c
sed -i 's|inversion_modules/phase_qc_shelf|inversion_modules/phase_4_postprocessing/4b_qc_triage|g' $FILES
sed -i 's|inversion_modules/breakpoint_pipeline|inversion_modules/phase_4_postprocessing/4c_breakpoint_refinement|g' $FILES
```

For casual `"phase_qc_shelf"` / `"breakpoint_pipeline"` mentions in
docs (without the `inversion_modules/` prefix), handle with separate,
narrower patterns to avoid false positives.

---

## Risks + mitigations

1. **MODULE_MAP phase_qc_shelf sibling language from pass 11.** Pass
   11's HANDOFF explicitly said "sibling, not subfolder." Pass 12 is
   the revision. Update the HANDOFF with a pointer to this renumber
   proposal.

2. **Pass 8's axis 5 wiring** in `compute_candidate_status.R`
   (4g_final_classification after renumber). The wiring uses
   `V7_FINAL_DIR` env var, not hardcoded paths, so it should survive
   the folder rename untouched. Verify nonetheless.

3. **Pass 9's STEP_{A,B,D}NN_ phase 3 rename.** Independent of this.
   No interaction.

4. **Pass 10's 5E archive.** Independent. No interaction.

5. **HPC paths.** User's LANTA tree mirrors the local tree under
   `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/`.
   Renumber requires a coordinated HPC + laptop + registry sweep.
   Existing output directories under the old names stay readable
   (they're under `FOLLOWUP_DIR`, not `MODULE_DIR`), so old results
   aren't orphaned.

6. **Pass 11 bridge script.** `phase_qc_shelf/scripts/bridge_from_phase4.sh`
   moves with the module to `4b_qc_triage/scripts/bridge_from_phase4.sh`.
   Its internal paths use `$here` so the move is safe.

---

## Execution order (for pass 12)

1. Take a checkpoint tarball of the baseline (already have working tree).
2. Directory moves (6 `mv` operations listed above).
3. Sed pass across all blast-list files.
4. Manual check: `compute_candidate_status.R` axis 5 + Gap 2 GDS blocks
   (from pass 8). Both should be path-insensitive via env vars.
5. Manual check: all LAUNCH scripts still chmod +x.
6. Run bash -n / python ast / MD fence checks.
7. Grep for any remaining old path strings (negative check).
8. Update `docs/MODULE_MAP.md` with the new structure.
9. Write new README at `phase_4_postprocessing/README.md` documenting
   the sequence.
10. Package + DROP_README with explicit rm commands for the old empty
    directories.

---

## What stays the same

- Phase names (`phase_4_postprocessing/`) — unchanged.
- All internal launcher names, script names, function names.
- Registry semantics — `reg$evidence`, `reg$samples`, `reg$intervals`,
  `reg$results` APIs unchanged.
- The pass 8 axis-5 wiring, pass 9 STEP_{A,B,D}NN_ phase 3 names,
  pass 10 MODULE_5E archive.
- Pass 11's bridge script logic (moves with the module).
- All per-candidate data products on disk (output directories under
  `FOLLOWUP_DIR` are not under `MODULE_DIR`).

---

## Alternative considered: keep siblings, add orchestrator

The lightweight alternative is to keep `phase_qc_shelf` and
`breakpoint_pipeline` as siblings and just add a master orchestrator
(`run_all_phase_4_and_prereqs.sh`) that chains them in the right
order. Lower blast radius — maybe 3 files.

**Rejected because:**
- Doesn't encode the dependency in the tree; future contributors /
  your future self still see three parallel modules.
- Orchestrator scripts tend to go stale. Tree structure doesn't.
- The renumber is a one-time cost; the clarity compounds.

The full renumber is worth the transient pain.
