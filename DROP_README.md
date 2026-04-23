# DROP_README — pass 8: docs triage + MODULE_MAP + axis5 wire + 4 new READMEs

**Pass:** 8 (after chat N+1 combined pass 5+6+archive)
**Date:** 2026-04-24
**Scope:** Four-part:
1. `docs/` triage (session artifacts archived, HANDOFF TODOs extracted)
2. New `docs/MODULE_MAP.md` (tree documentation, no renames)
3. Axis 5 wired into `compute_candidate_status.R` (T1 closed)
4. Four new READMEs: phase_2/2a, phase_2/2b, phase_5_followup, phase_6_secondary

---

## 1. `docs/` triage — 10 files → 6 files

Archived 4 session-artifact docs (work verified to be applied in tree):

| File | Where it went | Why |
|---|---|---|
| `HANDOFF_2026-04-23_v7_wiring_plan.md` | `_archive/chat_history/` | Deliverables (axis5 helper, assign_structural_class_v7, run_phase4_v7_wiring.sh, breakpoint_pipeline/01–07) all present in tree |
| `SESSION_AUDIT_append.md` | `_archive/chat_history/SESSION_AUDIT_append_2026-04-23.md` | Historical audit |
| `phase4_framework_v7_correction.md` | `_archive/chat_history/phase4_framework_v7_correction_2026-04-23.md` | Correction doc; its mandate ("use breakpoint_pipeline/, don't write new interior-structure code") is baked in |
| `HANDOFF.md` | `_archive/chat_history/HANDOFF_2026-04-23.md` | Priority 1-4 extracted into `docs/OPEN_TODOS.md` first; T1 then closed — see §3 below |

Deleted:

| File | Why |
|---|---|
| `docs/phase_qc_shelf_CHANGELOG.md` | md5-identical duplicate of canonical `inversion_modules/phase_qc_shelf/CHANGELOG.md` |

New:

- **`docs/OPEN_TODOS.md`** — clean todo tracker with three buckets:
  🟢 laptop-doable / 🟡 HPC-blocked / 🔵 deferred backlog. Convention:
  when done, delete the item (no "✅ completed" markers — git history
  records completion).

---

## 2. New `docs/MODULE_MAP.md`

Documents the two parallel module trees:

- `Modules/` (MODULE_1..4, preprocessing + SV callers)
- `inversion_modules/` (phase_1..7, the inversion pipeline)

Includes:
- What each module/phase does + which manuscript section it maps to
- Quick lookup table ("where does X happen?")
- Known numbering oddities (MODULE_5* migrated out of `Modules/`,
  legacy names in `phase_6_secondary/`, `breakpoint_pipeline/` has no
  phase prefix, etc.)

**No renames performed.** Documentation only. Renaming rejected as too
disruptive — too many `source(...)` / config paths would need
coordinated updates.

---

## 3. Axis 5 wiring — T1 closed

### The bug

`_axis5_final_label.R` existed as a standalone helper in
`4e_final_classification/` but was never sourced from
`compute_candidate_status.R`. `add_axis5()` /
`append_axis5_to_rows()` / `read_final_label()` were defined but
called from nowhere in live code.

### The wire

**`compute_candidate_status.R`**: inserted ~48 lines between
`status_dt <- rbindlist(...)` and `fwrite(status_dt, ...)` (around L932).

Env-gated via `V7_FINAL_DIR`:
- If unset → nothing changes, identical output to pre-wire behavior.
- If set + dir exists → sources helper, calls
  `append_axis5_to_rows(status_dt, axis5_dir)`, adds 4 columns to
  `candidate_status.tsv`:
  - `q_overall_structural_class`
  - `axis5_weakest_component`
  - `axis5_justification`
  - `axis5_source`
- If set but dir missing → warning, skips axis 5 (doesn't abort run).

Helper location is resolved robustly — works whether the script is run
via `Rscript --file=...` (uses `commandArgs(trailingOnly=FALSE)`) or
via `source()` (uses `sys.frames()$ofile`), with a CWD fallback.

### LAUNCH script passthrough

**`LAUNCH_characterize_classify.sh`**:
- Added `V7_FINAL_DIR` to the startup echo banner (shows
  `(unset — axis 5 disabled)` when not set so the state is visible).
- Added `export V7_FINAL_DIR="${V7_FINAL_DIR:-}"` before the
  `Rscript` call + short comment explaining what it gates.

### To actually use axis 5

```bash
export V7_FINAL_DIR="${BASE}/path/to/phase4_v7_blocks/final"
sbatch LAUNCH_characterize_classify.sh
```

Or inline:
```bash
V7_FINAL_DIR="${BASE}/.../final" sbatch LAUNCH_characterize_classify.sh
```

If `V7_FINAL_DIR` is unset or empty, behavior is backwards-compatible
with the previous `compute_candidate_status.R` — no axis 5 columns
appear, and the script's downstream consumers see the same schema as
before.

### Pre-ship verification

- Brace/paren balance on whole file: `(){}[]` all match (411/55/99
  matched pairs).
- Block-local: all balanced, no stray `<<-` assignments.
- `_axis5_final_label.R` function names verified: `read_final_label`,
  `append_axis5_to_rows`, `add_axis5` all defined at top level.
- `%||%` operator: helper redefines at its line 82; main script
  defines its own at line 34. Both have semantically equivalent
  definitions (a fallback operator for `is.null || length==0`).
  Sourcing inside the gated block overwrites the main's — benign.
- LAUNCH script: `bash -n` passes.

Note: I could not execute R in this sandbox (no R installed), so
runtime verification is deferred to HPC when you run this. The static
checks (balance, function-name resolution, no-op backwards
compatibility) are strong, but if anything surprises please send the
traceback and I'll fix on the next round.

### What's still undone for v7 end-to-end

The live edits from my pre-upload HANDOFF are T1 (axis 5 wiring —
done here). T2/T3/T4 are separate line-level edits in
`BREEDING_*.R` and `compute_candidate_status.R::build_key_spec()` that
still need reconstruction from whatever spec was in the prior chat.
See `docs/OPEN_TODOS.md`.

---

## 4. New READMEs

| Path | Lines | What's in it |
|---|---|---|
| `inversion_modules/phase_2_discovery/2a_local_pca/README.md` | 53 | BEAGLE→dosage→per-chr local-PCA windows. Workflow DAG, legacy (A02 monolithic) vs parallel (A03 stage1+stage2). Parameters (winsize=100, step=20, npc=2). |
| `inversion_modules/phase_2_discovery/2b_mds/README.md` | 57 | lostruct distances + per-chr MDS → candidate region seeds. MDS modes (per_chr, chunked_Nx). Candidate region assembly params. |
| `inversion_modules/phase_5_followup/README.md` | 105 | Umbrella for per-candidate deep analysis. Menu-style (not strict DAG) — C engines, R plotting, STEP20–41 legacy suite, LEGACY_FOLLOWUP_BRIDGE. |
| `inversion_modules/phase_6_secondary/README.md` | 71 | Umbrella for LD/Fst/Hobs secondary confirmation. Table mapping each sub-module to the biological prediction it tests. Naming-note about MODULE_5* legacy prefixes + pointer to RENAMING Cat 3. |

### Updated

- **`docs/MODULE_MAP.md`**: "Missing READMEs (backlog)" section replaced
  with "None as of 2026-04-24" + convention for future additions.
- **`docs/OPEN_TODOS.md`**: T1 removed; T3 cross-reference updated
  ("bundle with T2" instead of "with T1").

---

## File changes (live tree)

**Added:**
- `docs/OPEN_TODOS.md`
- `docs/MODULE_MAP.md`
- `inversion_modules/phase_2_discovery/2a_local_pca/README.md`
- `inversion_modules/phase_2_discovery/2b_mds/README.md`
- `inversion_modules/phase_5_followup/README.md`
- `inversion_modules/phase_6_secondary/README.md`

**Modified:**
- `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R`
  — added env-gated axis 5 block (~48 lines) between L930 and L982
- `inversion_modules/phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh`
  — added `V7_FINAL_DIR` to banner + export + comment (5 lines)

**Moved:**
- `docs/HANDOFF.md` → `_archive/chat_history/HANDOFF_2026-04-23.md`
- `docs/HANDOFF_2026-04-23_v7_wiring_plan.md` → `_archive/chat_history/`
- `docs/SESSION_AUDIT_append.md` → `_archive/chat_history/SESSION_AUDIT_append_2026-04-23.md`
- `docs/phase4_framework_v7_correction.md` → `_archive/chat_history/phase4_framework_v7_correction_2026-04-23.md`

**Deleted:**
- `docs/phase_qc_shelf_CHANGELOG.md` (duplicate)

---

## Manual `rm` commands

Nothing to remove manually. All changes are either additive, modifying
tracked files, or moves that GitHub Desktop will render as
delete-old + add-new.

---

## Commit message

```
pass 8: docs triage + MODULE_MAP + axis5 wire + 4 new READMEs

Part 1 — docs/ triage (10 → 6 files):
- Archive 4 done session-artifact docs to _archive/chat_history/
  (HANDOFF, SESSION_AUDIT_append, v7_wiring_plan, v7_correction)
- Extract HANDOFF Priority 1-4 into docs/OPEN_TODOS.md
- Delete duplicate docs/phase_qc_shelf_CHANGELOG.md (canonical copy
  at inversion_modules/phase_qc_shelf/CHANGELOG.md)

Part 2 — docs/MODULE_MAP.md:
- Document the Modules/ (MODULE_1..4) + inversion_modules/ (phase_1..7)
  dual-tree layout with manuscript-section mapping; no renames

Part 3 — Axis 5 wired (T1 closed):
- compute_candidate_status.R sources _axis5_final_label.R when
  V7_FINAL_DIR env var is set + valid; appends 4 columns to
  candidate_status.tsv. Backwards-compatible (no-op when V7_FINAL_DIR
  unset)
- LAUNCH_characterize_classify.sh echoes V7_FINAL_DIR state + exports
  it through to Rscript

Part 4 — 4 new READMEs:
- phase_2_discovery/2a_local_pca/     (dosage + local PCA windows)
- phase_2_discovery/2b_mds/           (lostruct + per-chr MDS)
- phase_5_followup/                   (per-cand deep analysis umbrella)
- phase_6_secondary/                  (LD/Fst/Hobs confirmation umbrella)

Missing-READMEs backlog is now empty per MODULE_MAP.md.
```

---

## Pre-ship verification checklist

- [x] `(){}[]` balance on all edited files
- [x] `bash -n` on `LAUNCH_characterize_classify.sh`
- [x] All paths referenced in `MODULE_MAP.md` exist (spot-checked 10+)
- [x] All paths in new phase_5 / phase_6 READMEs exist (verified 26)
- [x] No cohort-identity bugs in new content (grep F1/hybrid — clean)
- [x] No accidental legacy snake vocab (grep `snake_id|snake_phase|core_family` — only live in T7 rename-targets in OPEN_TODOS)
- [x] axis5 helper functions (`append_axis5_to_rows` etc.) exist and are called with correct argument shape
- [x] V7_FINAL_DIR env-gate is backwards-compatible (unset → identical prior behavior)
- [x] OPEN_TODOS updated (T1 removed, T3 cross-ref updated)
- [x] MODULE_MAP backlog section updated

---

## Next up (from `docs/OPEN_TODOS.md`)

🟢 laptop-doable:
- **T2**: BREEDING filter snippets (3 scripts × ~8 lines) — needs spec from a prior chat, probably wants a quick back-and-forth
- **T3**: Q7B / gene-conversion key-spec extensions — same
- **T4**: Step 6b GDS consistency check in `characterize_q5()` (~6 lines)

🔵 deferred:
- **T7**: RENAMING Cat 2 sed pass (big, one identifier class per tarball)
- **T10**: Move `_*.py` dev utils to `tools/`
- **T11**: Expand 2-line `README.md`

Pick whichever you want next.
