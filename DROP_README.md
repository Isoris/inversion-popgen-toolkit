# DROP_README — pass 14: design doc only (12-phase rename proposal)

**Pass:** 14
**Date:** 2026-04-24
**Format:** diff-only, docs only — NO code changes
**Scope:** Full blast-radius analysis + staged execution plan for the 12-phase rename. Execute in a fresh chat.

---

## ⚠️ No manual commands. Zero code changes this pass.

Two new docs, nothing deleted, nothing renamed.

---

## What's in this tarball

1. **`docs/TWELVE_PHASE_RENAME_PROPOSAL.md`**
   - The full design: phase_4 → phase_4/5/6/7/8/9, shift downstream phases 5/6/7 → 10/11/12
   - Blast radius enumeration: 159 files total (77 shell, 22 R, 7 Python, 50 MD, 3 JSON)
   - Staged execution plan: checkpoint → moves → regroup → sed tiers → docs → verify
   - Risks + don't-do list + the preserved-as-is list (pass 8/11/13 wiring)
   - Suggested opening prompt for the fresh chat

2. **`docs/SESSION_AUDIT_2026-04-24_chat-3-end.md`**
   - Full record of what passes 9, 10, 11-partial, 12, 13 shipped
   - Architectural decisions made (soft flag semantics, fold vs sibling, etc.)
   - Things I got wrong and recovered from (pass 12 DROP_README ordering bug, 4f subfolder flailing)
   - Repository state at end of chat 3 — confirmed from filesystem
   - Cohort identity + workflow conventions
   - Open non-pass-14 tasks

3. **This DROP_README.md**

---

## Why design-only this pass

Quentin asked for "Option 5 first" (see full cost) because:

- Chat 3 already shipped 5 passes (9, 10, 11-partial, 12, 13) and is substantially tired
- Pass 14 itself is 2.4× the blast radius of pass 12 (159 vs 66 files)
- Starting a major rename on tired context is how things break silently
- Execution in a fresh chat with the full proposal document in front of it will be cleaner

---

## Why the 12-phase rename

Quentin flagged three related concerns over the course of chat 3:

1. "phase 4 post processing is a strange name its too general"
2. "the biology is not only post process, its also real analysis"
3. "maybe it must be flattened and have many more phases?"

The diagnosis is correct. `phase_4_postprocessing/` currently packs 7 orthogonal jobs (catalog construction, QC triage, breakpoint refinement, group proposal, group validation, evidence biology, final classification) under one misleading label. The proposed 12-phase layout gives each its own top-level phase with a clear single-purpose name.

See `TWELVE_PHASE_RENAME_PROPOSAL.md` for the full rationale.

---

## What NOT to do this pass

- Do not attempt any directory moves
- Do not run any sed passes
- Do not update any launcher paths
- Do not rename the `phase_4_postprocessing/` folder yet

The proposal is the deliverable. Execution is a fresh-chat job.

---

## Commit message

```
pass 14: design doc for 12-phase rename

Adds docs/TWELVE_PHASE_RENAME_PROPOSAL.md with full blast-radius
analysis (159 files) and staged execution plan for flattening
phase_4_postprocessing into 6 dedicated phases (4..9) and shifting
downstream phases (5/6/7 -> 10/11/12).

Supersedes docs/PHASE4_RENUMBER_PROPOSAL.md (pass 12's internal
phase_4 renumber).

Also adds docs/SESSION_AUDIT_2026-04-24_chat-3-end.md documenting
what passes 9, 10, 11-partial, 12, 13 shipped and architectural
decisions made during chat 3.

No code changes. Execute the rename in a fresh chat using the
proposal document.
```

---

## Phase progress

- ✅ Pass 9: phase_3_refine Layer-tagged rename
- ✅ Pass 10: phase_6 MODULE_5E archive
- 🟡 Pass 11 (partial): qc_triage bridge script
- ✅ Pass 12: phase_4 restructure (folded qc_shelf + breakpoint_pipeline in)
- ✅ Pass 13: q_qc_shelf_* reader wired into 4g_final_classification
- 📋 Pass 14 (this pass): design doc for 12-phase rename
- ⏳ Pass 15 (next chat): execute the 12-phase rename per the design doc
