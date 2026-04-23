# DROP_README — pass 9: phase_3_refine Layer-tagged rename

**Pass:** 9
**Date:** 2026-04-24
**Format:** diff-only tarball
**Scope:** Rename the 6 numbered phase_3 scripts to `STEP_{A,B,D}NN_…`
tagging each step with the evidence Layer it feeds; rename the config
to match the `MODULE_*/00_moduleX_config.sh` convention used elsewhere;
update all 11 internal + external references.

This is the first pass of the phase 3-6 reorg plan. Phase 4 / 5 / 6
reorgs will ship as separate passes.

---

## ⚠️ Manual `rm` commands (required because diff-only)

```bash
cd /mnt/c/Users/quent/Desktop/inversion-popgen-toolkit/inversion_modules/phase_3_refine

rm 00_breakpoint_validation_config.sh
rm 01_extract_inv_candidates.sh
rm 02_extract_breakpoint_evidence.py
rm 03_statistical_tests_and_seeds.py
rm 04_validation_plots.py
rm 05_delly_manta_concordance.py
rm 06_bnd_inversion_signal.py
```

After these 7 `rm`s + drag-drop, GitHub Desktop will show 7 renames
(old → new filename) plus the modifications to the 4 reference files.

---

## Rename map

| Old | New | Layer |
|---|---|---|
| `00_breakpoint_validation_config.sh` | `00_phase3_config.sh` | shared config |
| `01_extract_inv_candidates.sh` | `STEP_A01_extract_inv_candidates.sh` | A upstream |
| `02_extract_breakpoint_evidence.py` | `STEP_A02_extract_breakpoint_evidence.py` | A upstream (→ D) |
| `03_statistical_tests_and_seeds.py` | `STEP_D03_statistical_tests_and_seeds.py` | **D** writer |
| `04_validation_plots.py` | `STEP_D04_validation_plots.py` | D diagnostic |
| `05_delly_manta_concordance.py` | `STEP_B05_delly_manta_concordance.py` | **B** report |
| `06_bnd_inversion_signal.py` | `STEP_B06_bnd_rescue.py` | **B** writer (rescue) |

No `STEP_C*` in phase 3 — Layer C (GHSL) is produced in
`phase_2_discovery/2e_ghsl/`, not here.

Unchanged (already descriptive, not part of the 01–06 numbered chain):
- `annotate_population_confidence.sh` (optional side-tool)
- `breakpoint_validator_standalone.py` (single-candidate standalone)
- `run_breakpoint_validation.sh` (orchestrator)
- `README.md`

---

## Files modified

### Inside `phase_3_refine/`

- **`00_phase3_config.sh`** — 1 internal header line + 8 `STEP0N` → `STEP_XNN` comment updates
- **`STEP_A01_extract_inv_candidates.sh`** — header comment, config source path, "Next: bash STEP_A02…" hint
- **`STEP_A02_extract_breakpoint_evidence.py`** — docstring header, usage example, argparse help ("STEP01" → "STEP_A01")
- **`STEP_B05_delly_manta_concordance.py`** — docstring header, 3 inline refs ("STEP01"/"STEP03" → "STEP_A01"/"STEP_D03")
- **`STEP_B06_bnd_rescue.py`** — docstring header, 1 inline ref, **registry provenance string** (`source_script="phase_3_refine/STEP_B06_bnd_rescue.py"`)
- **`STEP_D03_statistical_tests_and_seeds.py`** — docstring header, 3 inline refs, **registry provenance string** (`source_script="phase_3_refine/STEP_D03_statistical_tests_and_seeds.py"`)
- **`STEP_D04_validation_plots.py`** — docstring header + usage example
- **`run_breakpoint_validation.sh`** — all 6 `${SCRIPT_DIR}/0N_*` paths updated; 3 `STEP0N` comment mentions updated; 2 `source` lines pointing at config updated
- **`README.md`** — pipeline-steps block rewritten with STEP_XNN names + naming-convention explanation; 2026-04-24 layout note added; all inline step mentions updated (13 replacements via mechanical sweep)

### External references updated

- **`inversion_modules/phase_4_postprocessing/tests/test_registry_sanity.py`** — 2 `source_script="phase_3_refine/..."` strings + 1 comment
- **`inversion_modules/phase_4_postprocessing/specs/STRUCTURED_BLOCK_SCHEMAS.md`** — 1 code-comment example
- **`inversion_modules/phase_4_postprocessing/4a_existence_layers/README.md`** — Layer table (2 rows) + audit note (STEP_B06 / STEP_D03 refs)
- **`inversion_modules/CONFIG_ARCHITECTURE.md`** — 4 refs (TL;DR + tree + code example + See-also)
- **`inversion_modules/PATH_WIRING.md`** — 1 config-name ref
- **`registries/schemas/structured_block_schemas/INDEX_remaining_blocks.json`** — `source_script` field
- **`registries/schemas/structured_block_schemas/existence_layer_b_bnd_rescue.schema.json`** — description text + `source_script` field

### New docs

- **`docs/MODULE_MAP.md`** — added `### phase_3_refine/ scripts` sub-section documenting the step table + Layer-tag convention; updated phase_4 4e description from "axis 5 wiring (pending)" to "axis 5 (wired via V7_FINAL_DIR env)" reflecting pass 8.

---

## Why this layout

From the README's existing Layer table: phase 3 owns Layer D writing
(OR/Fisher test) and contributes supplementary Layer B (BND rescue).
But the old `01_…06_` names gave no hint which step fed which Layer —
you had to read the README table to figure it out. The `_A/_B/_D`
prefix makes each script's role self-identifying:

- `STEP_A01`, `STEP_A02`: upstream candidate discovery and per-sample
  evidence extraction — feeds the downstream Layer D tests.
- `STEP_D03`: the Layer D registry writer (the gate that decides
  VALIDATED tier in phase 4c).
- `STEP_D04`: Layer D diagnostic plots (no registry writes).
- `STEP_B05`: Layer B concordance report (reads unified table from
  STEP_A01).
- `STEP_B06`: Layer B BND rescue registry writer.

Grepping the toolkit for Layer-D-producing steps is now trivial:
`grep -rln "STEP_D" inversion_modules/`.

---

## Pre-ship verification

- [x] All 6 renames applied; old files removed from working tree
- [x] `bash -n` on all 4 `.sh` scripts in phase_3 passes
- [x] `python3 -m ast.parse` on all 6 `.py` scripts in phase_3 passes
- [x] Both edited JSON schemas (INDEX + bnd_rescue) parse as valid JSON
- [x] All markdown fence balances (MODULE_MAP, 4a/README, phase_3/README, PATH_WIRING, CONFIG_ARCHITECTURE): even
- [x] No remaining bare `STEP0N` token anywhere in the live tree
- [x] No remaining references to any of the 7 old filenames in the live tree (except the single intentional mention in the phase_3 README's 2026-04-24 layout note explaining the rename)
- [x] `source_script` fields in registry schemas + test fixtures use new filenames (important for registry audit trails)
- [x] No cohort-identity regressions in any edited file (grep F1/hybrid clean)

Could not runtime-test because R/pytest are not installed in this
sandbox. Static checks are strong. If anything surprises on LANTA,
paste the traceback and I'll fix next round.

---

## Commit message

```
pass 9: phase_3_refine Layer-tagged rename (STEP_{A,B,D}NN_)

Rename the 6 numbered phase_3 scripts from 01_...06_ to
STEP_{A,B,D}NN_ tagging each step with the evidence Layer it feeds
(A = upstream candidate discovery, B = SV caller layer, D = OR
association). No C tag because Layer C/GHSL lives in phase_2e.

Renamed:
  00_breakpoint_validation_config.sh -> 00_phase3_config.sh
  01_extract_inv_candidates.sh       -> STEP_A01_extract_inv_candidates.sh
  02_extract_breakpoint_evidence.py  -> STEP_A02_extract_breakpoint_evidence.py
  03_statistical_tests_and_seeds.py  -> STEP_D03_statistical_tests_and_seeds.py
  04_validation_plots.py             -> STEP_D04_validation_plots.py
  05_delly_manta_concordance.py      -> STEP_B05_delly_manta_concordance.py
  06_bnd_inversion_signal.py         -> STEP_B06_bnd_rescue.py

Updated references:
- run_breakpoint_validation.sh (all 6 step paths + config source)
- All 6 script internal docstrings + usage examples
- 2 registry source_script provenance strings (STEP_D03 and STEP_B06)
- phase_4/tests/test_registry_sanity.py (fixtures)
- phase_4/specs/STRUCTURED_BLOCK_SCHEMAS.md (example)
- phase_4/4a_existence_layers/README.md (Layer table)
- CONFIG_ARCHITECTURE.md (4 refs)
- PATH_WIRING.md (1 ref)
- 2 JSON schema files (INDEX + bnd_rescue source_script)
- docs/MODULE_MAP.md (new phase_3 sub-section + axis5 status update)
- phase_3/README.md (pipeline-steps block rewritten, 2026-04-24 note)

First pass of phase 3-6 reorg. Phase 4 / 5 / 6 to follow.
```

---

## Next up (phase 3-6 reorg continued)

Per the plan: phase_3 done. Options from the inventory I pitched earlier:

- **phase_4d** (14 flat files) — regroup into cheats/, bridges/, regenotype/, wiring/
- **phase_4b** (10 flat files) — regroup by concept: decompose/, multi_recomb/, nested_comp/, seal/
- **phase_6** — decide MODULE_5E fate (keep/consolidate/delete) — fast
- **phase_5** — 2-3 chats of its own

From `docs/OPEN_TODOS.md` (unchanged by this pass):
- T2 / T3 — BREEDING + Q7B wiring (needs spec reconstruction)
- T7 — RENAMING Cat 2 sed pass
- T12 — manuscript MAPQ inconsistency

Pick whichever you want next.
