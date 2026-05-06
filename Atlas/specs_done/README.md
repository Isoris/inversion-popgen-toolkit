# specs_done/ — shipped atlas + pipeline specs

Specs that have been fully implemented in `Inversion_atlas.html` (or the
companion HPC pipeline scripts) live here. Kept as historical reference:
the atlas code is the canonical implementation, the spec is the design
intent at the time of build.

## Why a top-level folder, not a subdir of `specs_todo/`?

The two original READMEs (`specs_todo/README.md` and
`specs_new_turn131/README.md`) disagreed on this — one said
`specs_shipped/`, the other said `specs_todo/_archive/specs_done/`.
Both READMEs have been updated to point here.

Top-level keeps symmetry with `specs_todo/` (forward-looking) and
`specs_new_turn131/` (pending-review). A spec's lifecycle is then:

```
specs_new_turn131/   ← drafted, awaiting review for build
       ↓
specs_todo/          ← approved, in build queue (or partially built)
       ↓
specs_done/          ← fully shipped, archive only
```

## Inventory (2026-05-05)

| Spec | Shipped where | Test file(s) |
|---|---|---|
| `SPEC_g_panel_unified_groups.md` | turns 135–136 | `test_turn135_g_panel_slice1.js`, `test_turn136_g_panel_karyotype_slice2.js` |
| `SPEC_observable_allele_h_label_system.md` | turns 139–140 | `test_turn139_h_label_classifier_slice1.js`, `test_turn140_h_label_chip.js` |
| `SPEC_sv_evidence_page.md` + `__per_candidate_folder` + `__upset_redirect` | turns leading up to 145 | `test_turn145_server_unify.js`, `tests/sv_evidence/` |
| `SPEC_comparative_te_breakpoint_fragility.md` + `__implementation_status` | v0.1 (turn 123) | `test_turn123_te_fragility.js` |

## Partially-shipped specs that DID NOT move here

These have at least one slice shipped but other slices flagged as
deferred. They live in `specs_todo/` with updated `**Status**:` lines:

- `SPEC_l3_het_dosage_coloring.md` (slice 1 shipped turn 129; slice 2 deferred)
- `SPEC_lines_panel_candidate_bands.md` (slice 1 shipped turn 141; slice 2 deferred)
- `SPEC_age_origin_panel.md` (main panel shipped turn 119; Q5 chip + card-flip deferred)
- `SPEC_l2_sweep_inheritance.md` (slices 1+2 shipped turns 133–134; mode-dispatch refactor pending)
- `SPEC_review_surfaces_auto_and_lineages.md` (slice 0 shipped turn 130; remaining UI slices unblocked but not built)

When the deferred slices land, those specs move here too.
