# HANDOFF — bundle reorganization, 2026-05-05

**Action**: organized handoffs and specs from the 2026-05-05 full bundle.
No atlas code touched. **2991 / 0 tests** still pass.

## What changed

### Handoffs (47 files → 3 buckets)

| Before | After |
|---|---|
| 36 `.md` at root + 11 in `handoffs_to_implement/` | 3 `.md` at root, rest under `handoffs/` |

```
Atlas/
├── HANDOFF_FULL_BUNDLE_2026-05-05.md          ← canonical session pointer (root)
├── DEPLOY_TO_LANTA.md, STEP6_DESIGN_NOTE.md   ← unchanged (root)
├── handoffs/
│   ├── README.md                              ← layout guide (new)
│   ├── current/                               ← 10 files (turn 151-158, this session)
│   │   ├── HANDOFF_2026-05-05_turn151_session_handoff.md
│   │   ├── HANDOFF_2026-05-05_turn151_slab_het_coloring.md
│   │   ├── HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md
│   │   ├── HANDOFF_2026-05-05_turn153_inheritance_auto_register.md
│   │   ├── HANDOFF_2026-05-05_turn154_compute_ux_hardening.md
│   │   ├── HANDOFF_2026-05-05_turn155_threshold_in_cache_key.md
│   │   ├── HANDOFF_2026-05-05_turn156_vshape_diagnostic.md
│   │   ├── HANDOFF_2026-05-05_turn158_popstats_qc_split.md
│   │   ├── turn157A_vshape_tooltip.md         ← renamed for clarity
│   │   └── turn157B_dosage_bridge_load.md     ← renamed for clarity
│   └── archive/
│       ├── 2026-05-05_turn128-150/            ← 25 files (Category B prior session)
│       └── pre_turn128/                       ← 11 files + turn129_review_bundle/
```

`handoffs_to_implement/` is gone (was empty after archiving).

### Specs

| Folder | Before | After |
|---|---|---|
| `specs_todo/` | 27 specs (mix of shipped, partial, forward-looking) | 20 specs (forward-looking + partially-shipped + meta) |
| `specs_done/` | did not exist | 7 spec files (6 specs + 1 paired status doc) |
| `specs_new_turn131/` | 21 specs awaiting review | unchanged |

#### Moved to `specs_done/` (fully shipped):

1. `SPEC_g_panel_unified_groups.md` — turns 135–136
2. `SPEC_observable_allele_h_label_system.md` — turns 139–140
3. `SPEC_sv_evidence_page.md` + `__per_candidate_folder` + `__upset_redirect`
4. `SPEC_comparative_te_breakpoint_fragility.md` + `__implementation_status` — turn 123 v0.1

All five specs had their **Status** lines updated to `SHIPPED in turn N`
with pointers to the corresponding test files.

#### Status-line updates on specs that stay in `specs_todo/`:

5 specs updated to make their partially-shipped state explicit:

- `SPEC_l3_het_dosage_coloring.md` — slice 1 shipped (turn 129); slice 2 deferred (already accurate, no change)
- `SPEC_lines_panel_candidate_bands.md` — slice 1 shipped (turn 141); slice 2 deferred
- `SPEC_age_origin_panel.md` — main panel shipped (turn 119); Q5 chip + card-flip deferred
- `SPEC_l2_sweep_inheritance.md` — slices 1+2 shipped (turns 133–134); dispatcher refactor pending
- `SPEC_review_surfaces_auto_and_lineages.md` — slice 0 shipped (turn 130); remaining UI slices unblocked

### READMEs

- `handoffs/README.md` (new) — layout guide for handoffs
- `specs_done/README.md` (new) — shipped-archive guide with inventory
- `specs_todo/README.md` (rewritten) — was 6-line stale inventory; now 20-spec table organized by lifecycle state
- `specs_new_turn131/README.md` (one-line edit) — promotion rule now points to `specs_done/` (top-level, not the proposed `specs_todo/_archive/specs_done/`)

The two READMEs that disagreed on the archive folder location (`specs_shipped/` vs `specs_todo/_archive/specs_done/`) are now both pointing to a top-level `specs_done/` for symmetry with `specs_todo/` and `specs_new_turn131/`.

## Verification

- `Inversion_atlas.html`: 73,773 LOC unchanged
- Tests: 2991 / 0 (re-run after reorganization, identical to bundle claim)
- 10 environment-broken pre-turn-128 tests (missing fixture files) unchanged — these are not counted in 2991

## What's still open

### Standing queue from FULL_BUNDLE §6

- A — UV refactor on locked groups (~200 LOC)
- C — Per-sample-lines het coloring (~300–500 LOC, perf concern at 226 × 30k)
- F — Per-group "show fish" expand toggle in G-panel inheritance card (~50 LOC)
- H — Cross-candidate V-shape gallery / sparkline matrix
- D — Pivot based on what het / inheritance / V-shape reveal on real data

### Parallel-track items #4 / #7 / #8

- #4 — G-bar / I·g / candidate karyotypes consumed by ancestry & popstats (200–400 LOC)
- #7 — LD track double-heatmap (common vs rare karyotype groups)
- #8 — Plot export PNG / SVG from canvases

### `specs_todo/` remaining work

12 forward-looking specs + 5 partially-shipped + 3 meta = 20 specs.
See `specs_todo/README.md` for the categorized inventory.

`specs_new_turn131/` has 21 specs awaiting Quentin's review-for-build.

End of reorganization handoff.
