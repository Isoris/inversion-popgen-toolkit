# Session audit — 2026-04-24 chat-3-end (design pass 14)

## What this chat delivered

**5 shipped passes + 1 design doc** for the fresh chat to execute.

### Shipped

1. **Pass 9** — phase_3_refine Layer-tagged rename (`STEP_{A,B,D}NN_` convention). 18 files, 60 KB.
2. **Pass 10** — phase_6 MODULE_5E archive (superseded by phase_qc_shelf Q07b/Q07c). 11 files, 29 KB.
3. **Pass 11 partial** — phase_qc_shelf bridge script from phase_4a catalog. 4 files, 11 KB. Tasks 2/3/4 deferred.
4. **Pass 12** — phase_4 restructure: folded phase_qc_shelf → 4b_qc_triage, breakpoint_pipeline → 4c_breakpoint_refinement, shifted 4b..4e → 4d..4g. 66 files, 998 KB. Major restructure.
5. **Pass 13** — q_qc_shelf_* reader wired into 4g/compute_candidate_status.R; 13 new keys, new `_qc_shelf_reader.R` helper (parallel to `_axis5_final_label.R` from pass 8), new `characterize_q_qc_shelf()` in characterize_candidate.R. Soft flag semantics preserved. 5 files, 33 KB. Actually tested with R installed in sandbox — full end-to-end smoke test passed.

### Design doc only (no moves yet)

6. **Pass 14** (this pass) — `docs/TWELVE_PHASE_RENAME_PROPOSAL.md` — full blast-radius analysis (159 files) and staged execution plan for flattening phase_4_postprocessing into 6 dedicated phases (phase_4 through phase_9) and shifting downstream phases (phase_10, 11, 12). User explicitly chose "Option 1 but after option 5 (see cost first) and in a fresh chat."

## Key architectural decisions made this chat

### Pass 12 decision: fold modules into phase_4, not keep as siblings

Originally in pass 11 I (wrongly) proposed that phase_qc_shelf stay a sibling of phase_4_postprocessing. Quentin pushed back — he wanted the data-dependency chain visible in the tree. Pass 12 executed: 6 directory moves + 66-file sed pass. Final layout was 7 sub-blocks 4a..4g.

### Pass 13 decision: soft flag semantics

Quentin's explicit guidance: "a messy inversion is still worth characterizing, just flagged as messy." Pass 13 honored this — Q_QC_SHELF never caps tier, never blocks characterization, just adds flag columns. Dispatched as "Step 1b" in `characterize_candidate()` with NO group gate. Pass 13's characterization_string kept `/7 ANSWERED` denominator (the 7 scientific questions Q1-Q7) and added `[qc_shelf=<flag>]` as a separate suffix.

### Pass 14 (design) decision: recognize the "phase_4 does too much" problem

Quentin's language: *"phase 4 post processing is a strange name its too general. maybe it must be flattened and have many more phases?"*. Also: *"the biology is not only post process, its also real analysis"*. Correct diagnosis. Phase_4 was doing catalog construction + QC triage + breakpoint refinement + group proposal + group validation + evidence biology + final classification — 7 things under one misleading label.

## Things I got wrong and recovered from

### Pass 11 initial architecture proposal

I proposed phase_qc_shelf stay a sibling of phase_4. User correctly pushed back — the dependency chain should be in the tree. Pass 12 fixed it.

### Pass 12 DROP_README

Initially listed the `rm -rf` commands for old folders at the BOTTOM of the DROP_README. GitHub Desktop's rename detection fails when old and new names coexist, so Quentin ended up with both the old 4b-4e names and the new 4b-4g names simultaneously. Required a fix-up message in the next turn explaining which folders to manually delete. Lesson learned: for any future rename pass, `rm` commands must come BEFORE drag-drop, at the top of the DROP_README, with a prominent warning.

### 4f subfolder design flailing

Took three attempts to get the 4f_group_dependent reorganization right:
1. First proposed: cheats/ sv_audit/ bridges/ regenotype/ (by file type — wrong)
2. Second: mechanism/ evolution/ robustness/ dropout_rescue/ (by biology — better)
3. Third (after Quentin pushed back hard): actually read the schemas to see what Q each script feeds, and realized:
   - Most 4f scripts have `group_validation_required: NONE` — so "group_dependent" name is misleading
   - The real biological classification is by Q1-Q7 (mechanism=Q4, age=Q5, audit=Q7, bridges=Q2+Q3)
   - `population_regenotype.py` doesn't belong in phase_4 at all — author's declared home is MODULE_4D/4E, STEP_A03b
   - `cross_species_bridge_v6.py` writes Q5+Q3 phase_4 keys, so it stays

Key lesson: when the user says "read the classification docs + the wiring in registries" — actually do it. The schemas had the authoritative Q-mapping all along. I was guessing from filenames when I should have been reading the `source_script` and `question` fields in the JSON schemas.

### Pass 14 scope flailing

At the end of the chat I realized my 4f subfolder work was a small piece of a much larger structural problem Quentin had raised three times in different ways:
- "idk maybe they fit after 4a or maybe they should be brought in 4b/c/d/e or maybe we need to make them 4b and change 4b/c/d/e to 4c/d/e/f" (pass 11/12 design)
- "maybe we need to pull it out of phase 4 then?" (pass 14 conversation)
- "phase 4 post processing is a strange name its too general" (pass 14 conversation)

The right response was to stop iterating on 4f and step back to the tree architecture. Quentin got there before I did; I should have paused sooner.

## Repository state at end of chat 3

- Working copy: `/home/claude/work/inversion-popgen-toolkit/`
- User's local: `/mnt/c/Users/quent/Desktop/inversion-popgen-toolkit/`
- HPC: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/`

Phase structure (confirmed from filesystem at end of chat):
```
inversion_modules/
├── phase_1_inputs/
├── phase_2_discovery/
├── phase_3_refine/              (STEP_{A,B,D}NN_ naming from pass 9)
├── phase_4_postprocessing/
│   ├── 4a_existence_layers/
│   ├── 4b_qc_triage/             (folded in from phase_qc_shelf, pass 12)
│   ├── 4c_breakpoint_refinement/ (folded in from breakpoint_pipeline, pass 12)
│   ├── 4d_group_proposal/        (was 4b, pass 12)
│   ├── 4e_group_validation/      (was 4c, pass 12)
│   ├── 4f_group_dependent/       (was 4d, pass 12 — still flat, 14 files)
│   ├── 4g_final_classification/  (was 4e, pass 12; has q_qc_shelf reader from pass 13)
│   └── orchestrator/ patches/ schemas/ specs/ tests/ docs/
├── phase_5_followup/
├── phase_6_secondary/            (MODULE_5E archived in pass 10)
├── phase_7_cargo/
└── utils/
```

All pass 9 / 10 / 12 / 13 shipped tarballs have been applied by the user (pass 11 partial was superseded by pass 12's inclusion of the bridge at the new location).

## Cohort identity (critical, do not conflate)

1. F1 hybrid (*C. gariepinus* × *C. macrocephalus*) — assembly paper ONLY
2. **226 pure *C. gariepinus* hatchery broodstock on LANTA — THIS paper (MS_Inversions_North_african_catfish). K=8 = broodline structure, not hybrid/population structure.**
3. Pure *C. macrocephalus* wild — future paper

**User's full name: Quentin Andres** (NOT Pasquet).

## Workflow conventions

- User drag-drops tarballs via GitHub Desktop (no git CLI on Windows laptop)
- Diff-only tarballs preferred over full-tree (60-270 KB vs 10 MB)
- **`rm` commands go at TOP of DROP_README, BEFORE drag-drop** (learned the hard way in pass 12)
- Every ship needs: `bash -n`, `python ast.parse`, R parse (if R installed), JSON parse, MD fence balance
- Every ship needs verification the previous passes' env-gated wiring survives: axis 5 (V7_FINAL_DIR) from pass 8, q_qc_shelf reader (QC_SHELF_EVIDENCE_DIR) from pass 13
- No emojis. Bilingual EN. "lol" register. Casual pushback expected and welcomed.
- User pushes back productively when Claude is wrong — listen, don't defend.
- When user says "I don't know," that's not a delegation — it's a flag that means "read the code properly before recommending, because I don't want to lock in the wrong thing." Observed twice this chat: on MODULE_5E verdict (where it turned out to be zero-callers + superseded by Q07b/Q07c, clear archive decision), and on 4f structure (where schemas had the authoritative Q-mapping I should have consulted first).

## For the fresh chat

**Suggested opening prompt:**

> Continuing from chat-3-end. Read `docs/TWELVE_PHASE_RENAME_PROPOSAL.md`
> and `docs/SESSION_AUDIT_2026-04-24_chat-3-end.md`. The user chose
> Option 1 (12-phase rename). Execute in stages — don't try to do it
> all at once. Start with stage 1 (checkpoint tarball), stage 2
> (directory moves), stage 3 (evidence_biology internal regroup +
> regenotype move), then stop and verify before stage 4 (sed pass).

**Open non-pass-14 tasks:**

- Pass 11 task 2: rewrite `4b_qc_triage/README.md` with 2-mode structure (single-candidate `run_chrom.sh` vs genome-wide `run_all_28chrom.sh --resume`)
- Pass 11 task 3: polish `docs/MODULE_MAP.md` qc_shelf entries
- OPEN_TODOS T2/T3 (BREEDING + Q7B wiring)
- OPEN_TODOS T7 (RENAMING Cat 2 sed pass)
- OPEN_TODOS T12 (manuscript MAPQ inconsistency)
- HPC-blocked: T5 full v7 run, T6 Engine H validation

These can interleave with pass 14 execution or come later.
