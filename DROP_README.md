# DROP_README — pass 8 (diff-only tarball)

**Pass:** 8
**Date:** 2026-04-24
**Format:** diff-only (~267 KB, 19 files) — much smaller than last
pass's full-tree approach per your feedback.
**Scope:** 6 tasks closed this session:
- docs/ triage + MODULE_MAP + OPEN_TODOS tracker
- Axis 5 wired (T1)
- Gap 2 GDS consistency check in `characterize_q5()` (T4)
- 4 new READMEs (phase_2/2a, phase_2/2b, phase_5_followup, phase_6_secondary)
- Root cruft moved into `tools/` (T10)
- Root README expanded from 2 lines to full orientation (T11)

---

## ⚠️ IMPORTANT — manual `rm` commands required

Because this is a **diff-only** tarball, drag-drop won't remove the
files that were deleted or moved-from-source. Run these after dropping:

```bash
cd /mnt/c/Users/quent/Desktop/inversion-popgen-toolkit

# Delete the duplicate CHANGELOG (canonical copy lives with its code)
rm docs/phase_qc_shelf_CHANGELOG.md

# Delete session-artifact docs (archived copies already in the tarball
# under _archive/chat_history/)
rm docs/HANDOFF.md
rm docs/HANDOFF_2026-04-23_v7_wiring_plan.md
rm docs/SESSION_AUDIT_append.md
rm docs/phase4_framework_v7_correction.md

# Delete the 4 _*.py dev utils from root (moved-to copies already in
# the tarball under tools/)
rm _bk_rename.py
rm _code_field_check.py
rm _rcheck.py
rm _schema_check.py
```

After these 9 `rm`s + drag-drop, GitHub Desktop will show:
- 4 renames (root `_*.py` → `tools/_*.py`)
- 4 "renames" of `docs/HANDOFF*.md` / `SESSION_AUDIT_append.md` /
  `phase4_framework_v7_correction.md` → `_archive/chat_history/…`
  (with slight filename tweaks to add date suffix — may appear as
  delete + add instead of rename depending on similarity detection)
- 1 delete of `docs/phase_qc_shelf_CHANGELOG.md`
- 3 modified files in `4e_final_classification/`
- 1 modified root `README.md`
- 6 new files: 4 READMEs + `docs/MODULE_MAP.md` + `docs/OPEN_TODOS.md`

---

## 1. `docs/` triage — 10 files → 6 files

Archived 4 session-artifact docs (work verified to be applied in tree):

| File | Where it went | Why |
|---|---|---|
| `HANDOFF_2026-04-23_v7_wiring_plan.md` | `_archive/chat_history/` | Deliverables (axis5 helper, assign_structural_class_v7, run_phase4_v7_wiring.sh, breakpoint_pipeline/01–07) all present in tree |
| `SESSION_AUDIT_append.md` | `_archive/chat_history/SESSION_AUDIT_append_2026-04-23.md` | Historical audit |
| `phase4_framework_v7_correction.md` | `_archive/chat_history/phase4_framework_v7_correction_2026-04-23.md` | Correction doc; its mandate ("use breakpoint_pipeline/, don't write new interior-structure code") is baked in |
| `HANDOFF.md` | `_archive/chat_history/HANDOFF_2026-04-23.md` | Priority 1-4 extracted into `docs/OPEN_TODOS.md` first — T1 and T4 then closed (see below) |

Deleted:

| File | Why |
|---|---|
| `docs/phase_qc_shelf_CHANGELOG.md` | md5-identical duplicate of canonical `inversion_modules/phase_qc_shelf/CHANGELOG.md` |

New:

- **`docs/OPEN_TODOS.md`** — clean todo tracker, three buckets
  (🟢 laptop-doable / 🟡 HPC-blocked / 🔵 deferred backlog). Down to
  9 open items from 13 at session start.

---

## 2. New `docs/MODULE_MAP.md`

Documents the two parallel module trees:

- `Modules/` (MODULE_1..4, preprocessing + SV callers)
- `inversion_modules/` (phase_1..7, the inversion pipeline)

Includes what each module/phase does, which manuscript section it maps
to, a quick-lookup table, and known numbering oddities. **No renames
performed** — documentation only.

---

## 3. Axis 5 wiring — T1 closed

### Problem

`_axis5_final_label.R` existed as a standalone helper but was never
sourced from `compute_candidate_status.R`.

### Fix

**`compute_candidate_status.R`**: inserted ~48 lines between
`status_dt <- rbindlist(...)` and `fwrite(...)`. Env-gated via
`V7_FINAL_DIR`:
- Unset → no-op, identical to pre-wire behavior.
- Set + dir exists → sources helper, adds 4 columns to
  `candidate_status.tsv` (`q_overall_structural_class`,
  `axis5_weakest_component`, `axis5_justification`, `axis5_source`).
- Set but dir missing → warning, skips axis 5 (doesn't abort).

Helper location resolved robustly (works under `Rscript --file=...`
and `source()`, with CWD fallback).

**`LAUNCH_characterize_classify.sh`**: added `V7_FINAL_DIR` to banner
echo + `export V7_FINAL_DIR="${V7_FINAL_DIR:-}"` before the Rscript call.

### To use

```bash
export V7_FINAL_DIR="${BASE}/path/to/phase4_v7_blocks/final"
sbatch LAUNCH_characterize_classify.sh
```

---

## 4. Gap 2 GDS consistency — T4 closed

### Problem

`cheat30` writes 4 cheat30 GDS keys (`q5_gds_within_std`,
`q5_gds_within_inv`, `q5_gds_between`, `q5_gds_het_pattern`) but
`characterize_q5()` never consumed them. The age inference from
`q5_gds_gap_percentile` wasn't being corroborated (or contradicted)
against the underlying between-vs-within separation that justifies
treating the gap as a real signal.

### Fix

**`characterize_candidate.R::characterize_q5()`**: inserted ~22 lines
between the Fst-agreement block and the convergence block (around
L500). Reads the cheat30 keys; if present:

- `bw_gap = gds_between - max(gds_within_std, gds_within_inv)`
- `bw_gap >= 0.05` → add `gds_between_within_consistent(gap=X)` to
  `ev_for` (supports age inference)
- `bw_gap < 0.05` → add `gds_between≈within (gap=X, age proxy may be
  noise)` to `ev_against` (contradicts: between ≈ within means the
  karyotype split isn't real at the genotype level)

Also: `q5_gds_het_pattern` → `ev_for` when present + meaningful.

Downstream effect: candidates with meaningful cheat30 agreement can
now reach `ANSWERED` status on Q5 (≥3 supporting evidence, 0 against);
candidates where the gds_gap is noise get demoted via the
`ev_against >= 2` convergence rule.

### Biology note

Even for "young" arrangements, between > within is expected — that's
what makes the karyotype assignment biologically real. Near-zero
`bw_gap` isn't "consistent with young"; it's "the inversion itself
isn't well-supported at the genotype-distance level." Hence `ev_against`
in both age regimes.

---

## 5. Four new READMEs

| Path | Lines | What's in it |
|---|---|---|
| `inversion_modules/phase_2_discovery/2a_local_pca/README.md` | 53 | BEAGLE→dosage→per-chr local-PCA windows. DAG, legacy-vs-parallel, params |
| `inversion_modules/phase_2_discovery/2b_mds/README.md` | 57 | lostruct + per-chr MDS → candidate region seeds. MDS modes, region params |
| `inversion_modules/phase_5_followup/README.md` | 105 | Umbrella — menu-style per-candidate deep analysis (C engines, R plotters, STEP20-41) |
| `inversion_modules/phase_6_secondary/README.md` | 71 | Umbrella — LD/Fst/Hobs secondary confirmation, with per-module prediction table |

Missing-READMEs backlog in `docs/MODULE_MAP.md` is now empty.

---

## 6. Root cruft — T10 closed

Moved 4 dev utilities from repo root to `tools/`:

- `_bk_rename.py` → `tools/_bk_rename.py`
- `_code_field_check.py` → `tools/_code_field_check.py`
- `_rcheck.py` → `tools/_rcheck.py`
- `_schema_check.py` → `tools/_schema_check.py`

Docstring paths updated to new locations. Scripts still expect to be
run from repo root (unchanged): `python3 tools/_schema_check.py`.

Verified: no live code sources these scripts. Only self-references
existed (in their own docstrings).

---

## 7. Root README — T11 closed

Expanded `README.md` from 2 lines to a full orientation document:

- What the repo does
- Repo layout in one paragraph + pointer to `docs/MODULE_MAP.md`
- Key docs table
- Dependencies
- Running on LANTA
- Dev utilities (under `tools/`)
- R packages (vendored under `R_packages/`)
- License
- Manuscript

Roughly 85 lines. Links relative so they render both on GitHub and
locally.

---

## File changes (live tree)

**Added (6):**
- `docs/OPEN_TODOS.md`
- `docs/MODULE_MAP.md`
- `inversion_modules/phase_2_discovery/2a_local_pca/README.md`
- `inversion_modules/phase_2_discovery/2b_mds/README.md`
- `inversion_modules/phase_5_followup/README.md`
- `inversion_modules/phase_6_secondary/README.md`

**Modified (4):**
- `README.md` (2 lines → 85 lines)
- `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R` (+48 lines axis 5)
- `inversion_modules/phase_4_postprocessing/4e_final_classification/LAUNCH_characterize_classify.sh` (+5 lines V7_FINAL_DIR plumbing)
- `inversion_modules/phase_4_postprocessing/4e_final_classification/characterize_candidate.R` (+22 lines Gap 2 GDS consistency)

**Moved (4 root → tools/):**
- `_bk_rename.py`
- `_code_field_check.py`
- `_rcheck.py`
- `_schema_check.py`

**Moved to archive (4):**
- `docs/HANDOFF.md` → `_archive/chat_history/HANDOFF_2026-04-23.md`
- `docs/HANDOFF_2026-04-23_v7_wiring_plan.md` → `_archive/chat_history/`
- `docs/SESSION_AUDIT_append.md` → `_archive/chat_history/SESSION_AUDIT_append_2026-04-23.md`
- `docs/phase4_framework_v7_correction.md` → `_archive/chat_history/phase4_framework_v7_correction_2026-04-23.md`

**Deleted (1):**
- `docs/phase_qc_shelf_CHANGELOG.md` (duplicate)

---

## Pre-ship verification

- [x] Brace/paren/bracket balance on `compute_candidate_status.R` (411/55/99 all matched)
- [x] Brace/paren/bracket balance on `characterize_candidate.R` (617/76/181 all matched)
- [x] `bash -n` on `LAUNCH_characterize_classify.sh` passes
- [x] `python3 -m ast.parse` on all 4 moved tools scripts: OK
- [x] All paths referenced in `MODULE_MAP.md` exist (spot-checked 10+)
- [x] All paths in new phase_5/phase_6 READMEs exist (26 verified)
- [x] Gap 2 insertion uses correct key names (`q5_gds_within_std`, `q5_gds_within_inv`, `q5_gds_between`, `q5_gds_het_pattern`) — cross-checked against `build_key_spec()` q5 list in `compute_candidate_status.R` L273
- [x] Gap 2 block biology reviewed (ev_against for bw_gap≈0 even for "young" candidates — see biology note above)
- [x] `safe_num()` and `` `%||%` `` helpers defined at top of `characterize_candidate.R` (lines 31 and 39 — verified)
- [x] No cohort-identity bugs in new content (grep F1/hybrid clean)
- [x] No accidental legacy snake vocab (only in T7 rename-targets)

Could not runtime-test because R is not installed in this sandbox.
Static checks are strong but if anything surprises on LANTA, paste
the traceback and I'll fix next round.

---

## Commit message

```
pass 8: docs triage + MODULE_MAP + axis5 + Gap2 GDS + tools/ move + README

Closed this session: T1 (axis5 wire), T4 (Gap2 GDS consistency in
characterize_q5), T10 (root _*.py → tools/), T11 (root README).

Part 1 — docs/ (10 → 6 files):
- Archive 4 done session-artifact docs to _archive/chat_history/
- Extract HANDOFF Priority 1-4 into docs/OPEN_TODOS.md
- Delete duplicate docs/phase_qc_shelf_CHANGELOG.md

Part 2 — docs/MODULE_MAP.md:
- Document Modules/ + inversion_modules/ trees → manuscript sections
- No renames (documentation only)

Part 3 — Axis 5 wired (T1 closed):
- compute_candidate_status.R sources _axis5_final_label.R when
  V7_FINAL_DIR env var is set; adds 4 columns to candidate_status.tsv
- Backwards-compatible (no-op when V7_FINAL_DIR unset)
- LAUNCH_characterize_classify.sh echoes + exports V7_FINAL_DIR

Part 4 — Gap 2 GDS consistency (T4 closed):
- characterize_q5() now consumes q5_gds_within_std, q5_gds_within_inv,
  q5_gds_between, q5_gds_het_pattern
- bw_gap gate adds corroborating or contradicting evidence to age
  inference; downstream convergence rules unchanged

Part 5 — 4 new READMEs (missing-READMEs backlog empty):
- phase_2_discovery/2a_local_pca/, 2b_mds/
- phase_5_followup/, phase_6_secondary/

Part 6 — root cruft (T10):
- Move _bk_rename.py, _code_field_check.py, _rcheck.py, _schema_check.py
  to tools/. Docstring paths updated. No live callers changed.

Part 7 — root README (T11):
- Expand from 2 lines to full orientation doc with key-docs table,
  dependencies, LANTA usage, pointer to MODULE_MAP.md.
```

---

## Next up (from `docs/OPEN_TODOS.md`)

🟢 laptop-doable (2 items left):
- **T2**: BREEDING filter snippets (3 scripts × ~8 lines) — needs spec from a prior chat, probably wants a back-and-forth
- **T3**: Q7B audit keys + gene-conversion Q2 extensions — same

🔵 deferred:
- **T7**: RENAMING Cat 2 sed pass (big, one identifier class per tarball)
- **T12**: Manuscript MAPQ=60 wording (wait until you pick option a or b)
- **T13**: `MODULE_1_Methods_Revised.md` TODOs — don't touch until asked

🟡 HPC-blocked:
- **T5**, **T6**: full v7 run + Engine H validation on LANTA Monday

Pick whichever you want next.
