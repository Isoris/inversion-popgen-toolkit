# OPEN_TODOS.md — live work queue

Clean tracker of work items that are genuinely pending. Extracted from
`docs/HANDOFF.md` (now archived) during the 2026-04-24 docs triage.
Session-artifact narrative stripped; kept only actionable items.

---

## 🟢 Laptop-doable (no HPC needed)

### T2. Step 7 — BREEDING filter snippets

**Status.** `BREEDING_A_broodstock_compatibility.R`,
`BREEDING_C_founder_haplotype_tracker.R`,
`BREEDING_D_recombination_atlas.R` all exist in
`inversion_modules/phase_7_cargo/extra_plots/breeding/`. The prior
session planned an ~8-line snippet per script that reads
`candidate_status.tsv` and filters to supported balanced inversions
only. Not applied.

**Blocker.** None — R edits. Natural to bundle with T3 in one tarball.

### T3. Step 8 gaps

- **Gap 1 (Q7B audit keys).** Extend `build_key_spec()` q7 list in
  `compute_candidate_status.R` to cover the Q7B audit keys.
- **Gap 4 (gene conversion).** Extend q2 list likewise.

**Blocker.** None. Natural to bundle with T2 since both touch BREEDING /
`compute_candidate_status.R` in the same pass.

### T4. Step 6b — GDS consistency (Gap 2)

**Status.** ~6-line addition in `characterize_q5()` inside
`characterize_candidate.R`.

**Blocker.** None. Lower priority per prior handoff.

---

## 🟡 HPC-blocked (wait for LANTA)

### T5. Full-cohort v7 phase_4 classification run

The actual end-to-end validation of the v7 wiring against the 226-sample
cohort. Requires `breakpoint_pipeline/` outputs to exist on HPC.
Expected check afterwards:

```bash
awk -F'\t' 'NR>1 {print $2}' ${BASE}/phase4_v7_blocks/final/final_catalog.tsv \
  | sort | uniq -c | sort -rn
```

Tuning knob if `complex_rearrangement_out_of_scope > 50%`: drop rule 1
threshold in `assign_structural_class_v7.py` from 0.70 → 0.80.

### T6. Engine H validation on LANTA

v3.1 shipped the Hobs-per-group layer untested on HPC.

---

## 🔵 Deferred — broader backlog

These come from `HANDOFF_2026-04-24.md` (top-level handoff):

### T7. RENAMING.md Category 2 — code-identifier renames

Sed recipes already written in
`inversion_modules/phase_2_discovery/2c_precomp/RENAMING.md` §2.3:

- `cheat[N]_*` → `test[NN]_*`  (~100+ hits)
- `snake_id` → `region_id`
- `snake_phase` → `extension_phase`
- `core_family` → `scale_tier`
- ~300 lines of live snake identifier hits

Big sed pass. One tarball per identifier class to keep diffs reviewable.

### T8. RENAMING.md Category 3 — cross-module path renames

`SNAKE*_DIR` env vars; `snake_regions_multiscale/` dir; `flashlight_v2/`
path refs. Deferred per RENAMING.md §5 until HPC available for
coordinated rename.

### T9. `reg$stats` deprecated alias

3 hits in `registries/api/R/registry_loader.R`. Low priority.

### T10. Root cruft — move `_*.py` dev utilities into `tools/`

- `_bk_rename.py`
- `_code_field_check.py`
- `_rcheck.py`
- `_schema_check.py`

### T11. Expand root `README.md` (currently 2 lines)

### T12. Manuscript v19 MAPQ=60 inconsistency

MAPQ=60 justification says "homeologous regions between subgenomes"
but mapping reference is stated as Gar-only 28 pseudochromosomes. Pick
one for v20:
- (a) reword MAPQ justification to "within-subgenome paralogue / repeat
  mismapping"
- (b) confirm mapping was against full 55-pseudochromosome reference
  and update earlier paragraph

`MODULE_1_Methods_Revised.md` §1.6 inherits the same; follow whichever
path user picks.

### T13. `MODULE_1_Methods_Revised.md` TODOs

Has `[TODO: ...]` placeholders throughout because experiments aren't
fully reported yet. Don't touch until explicitly asked.

---

## Conventions

- When a todo gets picked up in a session, it stays here until the
  corresponding code change is committed.
- When done, delete the item — no "✅ completed" markers (those belong
  in git history and audit logs).
- Keep this file short. If it grows past ~150 lines, a cleanup pass
  is overdue.
