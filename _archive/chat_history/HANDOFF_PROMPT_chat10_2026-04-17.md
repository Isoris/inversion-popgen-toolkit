# Handoff prompt for chat 10

Continuing from chat 9 (2026-04-17). The tarball is the chat-9 output:
`inversion-popgen-toolkit_chat9_fixes_plus_phase4e_partial_2026-04-17.tar`.

## What just happened (chat 9, in one paragraph)

Quentin asked me to understand the full workflow before auditing, then
said "code anything missing" for phase 4e. Scope expanded from the
planned Option-A phase-4c audit into a phase-4 end-to-end connection
audit plus substantial 4e rework. Six fixes landed (FIX 39–44). The
4e driver (`run_characterize.R`) and `SPEC_VS_REALITY.md` writer-audit
doc are the main pending deliverables — chat 9 ran out of budget
before writing them. This chat's job is to finish those plus anything
else that keeps 4e from actually running end-to-end.

## READING ORDER — do not skip

You must understand the full phase-4 workflow BEFORE touching code.
Chat 9's posture (understand-before-fix) is the baseline. Read in this
order; don't jump ahead:

### Step 1: The big picture (≈15 min of reading budget)

1. **`AUDIT_LOG_chat9_2026-04-17.md`** — the six fixes and ten findings
   from the chat you're continuing. This is the most recent context.
2. **`inversion_modules/phase_4_postprocessing/README.md`** — the
   canonical data-flow diagram for phase 4 (4a→4b→4c→4d→4e).
3. **`inversion_modules/phase_4_postprocessing/docs/PHASE4_ARCHITECTURE.md`**
   — the registry-as-catalog design, group validation gate rules,
   per-question group requirements table (section 4 is the
   authoritative Q→group-validation map).
4. **`inversion_modules/phase_4_postprocessing/docs/PHASE4B_REWRITE_ARCHITECTURE.md`**
   — the 4b four-script architecture and seal's resolution rules.
5. **`inversion_modules/phase_4_postprocessing/specs/INVERSION_REGISTRY_SPECIFICATION_v2.md`**
   — the 352-key authoritative spec. This is what 4e's
   `build_key_spec()` targets.
6. **`inversion_modules/phase_4_postprocessing/specs/CHARACTERIZATION_CONVERGENCE_RULES.md`**
   — the per-Q convergence logic characterize_candidate implements.

### Step 2: Phase 2/3 quick pass (high-level only, ~10 min)

Most phase 2/3 code doesn't need reading for this chat — Quentin said
"high level for phase 2–3, focus on phase 4e." The only thing you need
to know about phase 2/3 from a 4e perspective is: they are the
UPSTREAM WRITERS of keys that 4e will consume. Specifically:

- **Phase 2c** (flashlight C00): produces SV priors. Writes Layer B
  keys (`q7_layer_b_*`) and some Q1 mechanism context.
- **Phase 2d** (inv_detect v9.3, triangles): produces candidate coords
  (`triangle_intervals.tsv.gz`). Writes Q1 shape keys (squareness,
  NN persistence, shape_class, d01–d12).
- **Phase 2e** (GHSL C04_snake3_ghsl_v5): produces Layer C haplotype
  contrast. Writes `q7_layer_c_ghsl_*` keys.
- **Phase 3 STEP03** (statistical_tests_and_seeds.py): runs Fisher
  exact + Cochran-Armitage on 3×2 contingency table from per-sample
  PE/SR counts at breakpoints. Writes Layer D keys
  (`q7_layer_d_fisher_or`, `q7_layer_d_fisher_p`,
  `q7_layer_d_armitage_z`, `q7_layer_d_armitage_p`,
  `q7_layer_d_concordance`). This is the "OR test" in chat 9's
  mental model.

For the 4e task: do NOT re-read the phase-2/3 source. Trust the spec.
If you find a missing writer, add it to the `SPEC_VS_REALITY.md`
aspirational list for chat 11+ to wire, not to this chat.

### Step 3: Phase 4 complete audit trail (≈20 min)

The connection-by-connection state of phase 4 as of chat 9:

| Connection | Status |
|---|---|
| 2d → 4a (candidate birth) | Not re-verified chat 9; assumed working |
| Phase 3 → 4c (Layer D read) | INTACT (chat 7 FIX 31 holds) |
| 4b.1/4b.2/4b.3 → 4b.4 (seal synthesis) | INTACT (chat 8 verified) |
| 4b.4 → 4c (groups + promotion_cap) | INTACT (chat 9 verified) |
| 4c → hypothesis_verdicts.tsv | INTACT (chat 9 verified L2355–2425) |
| 4c list-return unpacking → verdict row | INTACT (chat 9 verified 5 cols) |
| T9 verdict string → `compute_group_validation` | **FIXED** FIX 39 |
| seal → 4e via q1_composite_flag | **FIXED** FIX 40 (one writer now) |
| cheat30 internal `compute_pairwise_ibs` call | **FIXED** FIX 41 |
| jackknife_status aspirational → fall back to q6_family_linkage | **FIXED** FIX 43 |
| characterize_candidate → source group_validation_gate | **FIXED** FIX 44 |
| 4e build_key_spec vs v2 spec | **FIXED** FIX 42 (367 keys, aspirational flagged) |
| `characterize_candidate()` has no driver | **STILL BROKEN** → your Priority 1 |
| `run_phase4.sh` calls old-style 4b launcher | **STILL BROKEN** → your Priority 3 |
| 4d launchers stale paths | **STILL BROKEN** → your Priority 3 |

## TASKS FOR THIS CHAT

### PRIORITY 1 — Build `run_characterize.R` (the 4e driver)

**Where:** `inversion_modules/phase_4_postprocessing/4e_final_classification/run_characterize.R`

**What it does:** reads per-candidate registry keys, calls
`characterize_candidate(keys)` from characterize_candidate.R, writes
outputs per candidate.

**CLI contract** (mirror `compute_candidate_status.R`):

```bash
Rscript run_characterize.R <registry_dir> <outdir>
```

**Reuse compute_candidate_status's loader pattern** — the three-tier
fallback at `compute_candidate_status.R` L518–556 (per-candidate
`*.keys.tsv`, then merged `evidence_registry.tsv`, then reconstruct
from `candidate_scores.tsv.gz`). Don't rewrite the loader from scratch;
source or copy the pattern.

**Must source before running:**
- `characterize_candidate.R` (sibling file; gives the characterize_qN
  functions + `characterize_candidate()` + `format_characterization()`)
- `group_validation_gate.R` (in `../4c_group_validation/`; gives
  `assess_group_validation` + `check_group_gate`). FIX 44 made
  characterize_candidate auto-source this, but explicit sourcing in
  the driver is safer for SLURM runs with weird cwd.

**Outputs:**

1. `characterization.tsv` — one row per candidate with columns:
   ```
   candidate_id  group_level  characterization_string
   Q1_status  Q1_label  Q1_reason
   Q2_status  Q2_label  Q2_reason
   ... Q3..Q7 ...
   n_answered  n_partial  n_measured  n_contradicted  n_empty
   ```

2. Per-candidate `per_candidate/<cid>/characterization.txt` — use
   `format_characterization(cid, char_result)` which is already
   defined in characterize_candidate.R L584.

3. `characterization_summary.txt` — genome-wide summary: how many
   candidates at each status cut, distribution of evolutionary
   labels, top CONTRADICTED candidates flagged for review.

**Error handling:** if a candidate has NO keys (registry miss),
produce a row with all Qs as EMPTY rather than crashing. Log to
stderr. Pipeline must be robust to partial registry state because
aspirational keys will be missing for most candidates.

**Integration with `compute_candidate_status.R`:** after producing
`characterization.tsv`, the driver should ALSO produce a combined
file `candidate_full_report.tsv` that merges:
- `compute_candidate_status.R`'s `candidate_status.tsv` (tier +
  completion + evolutionary class)
- `run_characterize.R`'s `characterization.tsv` (per-Q convergence)
on `candidate_id`. This is the "one row, all axes" deliverable.
If `candidate_status.tsv` isn't in `registry_dir`, skip the merge
silently.

### PRIORITY 2 — Write `SPEC_VS_REALITY.md` (the writer-audit doc)

**Where:** `inversion_modules/phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md`

**Structure:**

```markdown
# Spec vs reality — what 4e expects vs what gets written

Generated: <date>
Spec version: v2 (352 keys target; 4e code contains 367 due to dual aliases)

## Summary

| Question | Spec | Wired | Aspirational | % wired |
|---|---|---|---|---|
| Q1 | 49  | XX | YY | ZZ% |
| Q2 | 55  | ...
| Q3 | 74  | ...
| Q4 | 47  | ...
| Q5 | 39  | ...
| Q6 | 28  | ...
| Q7 | 75  | ...
| **Total** | 367 | ... | ... | ... |

## Aspirational keys by module (TODO for wiring work)

### Phase 2 — needs wiring
- Module 2B: ... (keys it produces but doesn't register)
- GHSL v5: q7_layer_c_ghsl_* (all 7 — pipeline exists, no registry writes)
- ...

### Phase 3 — needs wiring
- STEP03 Armitage: q7_layer_d_armitage_z, q7_layer_d_armitage_p, q7_layer_d_concordance

### Phase 4 — needs wiring
- 4d cheat28 tandem repeats: q4_n_tr_at_bp_*, q4_dominant_tr_class_*, ...
- 4d breakpoint audit: q7b_delly_* (partial — some written, some not)
- ...

## Keys with writers (reference; NOT aspirational)

### Phase 4b (all wired per chat 8 verification)
- internal_dynamics.schema.json extracts: q2_silhouette_score, q2_bic_gap_k3_vs_k2, ...
- recombinant_map.schema.json extracts: q2_n_recombinant, q2_n_gene_conversion, ...
- internal_ancestry_composition.schema.json extracts: q1_composite_flag (FIX 40), ...
- frequency.v2.schema.json extracts: q6_freq_inv, q6_family_linkage, q6_polymorphism_class, ...

### Phase 4c (all wired post FIX 31/32/33/39)
- hypothesis_verdicts.tsv columns (via register_C01f_keys helper):
  q6_group_validation, q6_family_linkage, q6_quality_flags, ...
```

**Source of truth:** the `*_aspir` lists in `compute_candidate_status.R::build_key_spec()`
(FIX 42). Walk those lists programmatically — don't reinvent. Use a
small R script or grep to generate the first draft of the doc, then
add prose.

### PRIORITY 3 — Orchestrator + launcher breaks

**Only do this if Priority 1 + 2 are done with budget remaining.**

1. `orchestrator/run_phase4.sh` calls `LAUNCH_C01i_decomposition.sh`
   as one job. 4b is now 4 scripts. Either (a) change run_phase4.sh
   to call `run_phase4b.sh` for the 4b block, or (b) create
   `LAUNCH_C01i_decomposition.sh` as a thin wrapper that calls
   run_phase4b.sh. Option (a) is cleaner.

2. `orchestrator/LAUNCH_group_cheats_example.sh` has stale paths
   (references `cheats/`, `scripts/`, `burden/` subdirs that don't
   exist). Also references scripts `compute_age_fst_subblock.R`
   and `STEP_C01f_c_burden_regression.R` that aren't in the
   tarball. For this chat: update the paths to point at
   `4d_group_dependent/` and rename to `LAUNCH_group_cheats.sh`
   (production). Flag the two missing scripts as aspirational in
   `SPEC_VS_REALITY.md`.

### PRIORITY 4 — Stretch goals (only if budget permits)

- Phase 4a audit (C01d, C01e, C01g) — chat 9 never got to these.
- Layer C (GHSL) read contract: what does C01d do when Clair3 is
  incomplete on some chromosomes?
- The `register_C01f_keys` external helper (Finding W) — lives in
  `utils/registry_key_helpers.R` OUTSIDE the tarball. Verify it
  handles the 5 new verdict columns (group_validation_after,
  quality_flags, family_linkage, promotion_cap_applied, comp_source).

## POSTURE FOR CHAT 10

**Understand before coding.** Read all the docs in Step 1 even if
you think you know the workflow. Chat 9 had to recover from two
near-misses (FIX 38 v1 with wrong impact story, Option A-vs-B on
composite_flag that I got partially wrong initially). The cost of
re-reading is lower than the cost of writing something that doesn't
fit the actual contracts.

**No three-fix ceiling this chat.** Quentin lifted it in chat 9 when
he authorized "code anything missing." But that's not a license for
sloppy work — each code artifact should be (a) bug-backed or
feature-backed by something in the findings list, (b) surgical in
scope, (c) parse-checked before commit. The `run_characterize.R`
driver is probably 150–200 lines and should be written once
cleanly, not incrementally patched.

**Budget planning (chat 10):**
- 15% reading (the docs above + your sources of truth for the driver)
- 10% figuring out the loader contract (how the three fallback paths
  in compute_candidate_status combine; how to reuse without duplicating)
- 40% writing run_characterize.R + testing brace balance / parse
- 15% writing SPEC_VS_REALITY.md (walk the aspirational lists)
- 10% orchestrator/launcher fixes if time permits
- 10% writeup + repack

If budget runs low during writeup: the audit log can be abbreviated,
the session summary can be bullets. But the driver MUST be complete
+ parse-balanced before chat end. Shipping a half-written driver is
worse than shipping no driver.

## OPEN ITEMS CARRIED FORWARD (cumulative; chat 10 doesn't need to fix these)

From chat 5: Finding S — STEP03 seed wiring. Design call needed.
From chat 6: lib design work, blocks registry-block path.
From chat 7 Findings 8 & 9: 4a writes no registry blocks without lib;
  C01g silent skip on missing helper.
From chat 8 Finding T: README structure_type drift.
From chat 8 Finding U: `q2_decomp_quality_flags` forward-declared.
From chat 8 Finding V: cheat24 threshold divergence — design call needed.
From chat 9 Finding W: `register_C01f_keys` external helper needs
  verification for new verdict columns.
From chat 9 Finding Y: 4d production launcher missing.
From chat 9 Finding Z: 4d README stale.
From chat 9 Finding AA: population_regenotype.py misfiled.
From chat 9 Finding AB: 4e key-count comments off by a few (cosmetic,
  partially fixed by FIX 42).

## PARSE-CHECK BACKLOG (cumulative; requires LANTA Rscript)

Thirteen items need LANTA `Rscript -e "parse(file=...)"` + smoke test
before production:

From chat 7:
- phase_3_refine/03_statistical_tests_and_seeds.py (Python AST-parse OK; pipeline smoke OK)
- 4c_group_validation/STEP_C01f_hypothesis_tests.R
- 4c_group_validation/group_validation_gate.R
- 4e_final_classification/compute_candidate_status.R
- 4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R
- 4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R

From chat 8:
- 4b_group_proposal/STEP_C01i_b_multi_recomb.R
- 4b_group_proposal/lib_decompose_helpers.R

From chat 9:
- 4c_group_validation/STEP_C01f_hypothesis_tests.R (RE-PARSE after FIX 39)
- 4c_group_validation/group_validation_gate.R (RE-PARSE after FIX 43)
- 4d_group_dependent/cheat30_gds_by_genotype.R (FIX 41)
- 4e_final_classification/compute_candidate_status.R (RE-PARSE after FIX 42)
- 4e_final_classification/characterize_candidate.R (FIX 43 + FIX 44)
- 4b_group_proposal/STEP_C01i_d_seal.R (FIX 40)
- schemas/internal_ancestry_composition.schema.json (FIX 40; JSON
  verified in-container but schema validator should re-run)

Add to the backlog whatever chat 10 produces.

## END OF CHAT 10

When done, write `AUDIT_LOG_chat10_2026-04-17.md`, append to
`SESSION_SUMMARY_2026-04-17.md` and `FIXES_APPLIED_2026-04-17.md`.
Repack as
`inversion-popgen-toolkit_chat10_phase4e_complete_2026-04-17.tar`
(or `phase4e_plus_orchestrator` if Priority 3 was tackled too).

Remember: the goal is **4e runs end-to-end on real data**. When
Quentin can run both `compute_candidate_status.R` AND
`run_characterize.R` on the LANTA registry and get `characterization.tsv`
+ `candidate_full_report.tsv` out, the finish line is crossed for
4e. Everything else (orchestrator, phase-2/3 wiring, phase-4a audit)
is polish around that finish line.
