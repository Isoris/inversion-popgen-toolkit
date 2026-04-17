# AUDIT LOG — chat 9 of the phase-4 audit series (2026-04-17) — PART F

Scope pivoted mid-chat. Initial plan was Option A (phase 4c audit) with
three-fix ceiling. After Quentin asked for full workflow understanding +
authorization to "code anything missing" for phase 4e, the ceiling was
lifted and scope expanded to an end-to-end phase-4 audit + substantial
4e rework.

Six fixes (39–44) applied. The 4e driver (`run_characterize.R`) and
`SPEC_VS_REALITY.md` doc remain as next-chat work because the chat
budget was nearly exhausted before those could be written.

## Reading done this chat

- Top-level `phase_4_postprocessing/README.md` (180 lines) — data flow diagram
  for 4a→4b→4c→4d→4e
- `docs/PHASE4_ARCHITECTURE.md` (357 lines) — registry-as-catalog + group
  validation gate rules + per-question group requirements
- `docs/PHASE4B_REWRITE_ARCHITECTURE.md` (280 lines) — 4b.1/4b.2/4b.3/4b.4
  contracts including composite_flag / promotion_cap semantics
- `4c_group_validation/README.md` and the three integrated patches
  (01/02/03) in `patches/`
- Sampled sections of `STEP_C01f_hypothesis_tests.R` (2523 lines) — the
  integrated helpers (L286–561), main-loop call site (L2355–2425), verdict
  row construction, registry write path
- `STEP_C01i_d_seal.R` — write paths for q6_group_validation,
  q6_validation_promotion_cap, q1_composite_flag, family_linkage
- `schemas/internal_ancestry_composition.schema.json`,
  `schemas/frequency.v2.schema.json`, `schemas/internal_dynamics.schema.json`
  — for keys_extracted auditing
- `orchestrator/run_phase4.sh`, `run_phase4b.sh`, `LAUNCH_group_cheats_example.sh`
- `4d_group_dependent/` full contents (5 scripts + 2 aux), each preamble
- `4e_final_classification/compute_candidate_status.R` entire file
- `4e_final_classification/characterize_candidate.R` entire file
- `4c_group_validation/group_validation_gate.R` entire file
- `specs/INVERSION_REGISTRY_SPECIFICATION_v2.md` — the authoritative
  v2 key spec (352 keys)
- Past-chat search for authoritative key counts and architectural history

Total ≈ 5500 lines read + 6 docs.

## The mental model converged in this chat

Quentin confirmed the workflow intent:
- **Four orthogonal evidence streams** (Layer A PCA / Layer B SV / Layer C
  GHSL / Layer D Fisher OR) validate candidate existence from independent
  angles.
- **Recombinants are a SUB-PHENOMENON of inversions, not parallel
  arrangements.** They get registered for traceability but excluded from
  `comp` for T1/T2/T7/T8/T10 because their position-dependent class
  violates "one band per sample". Q2 internal_dynamics characterizes them;
  Q6 counts them separately from HOM_REF/HET/HOM_INV.
- **Jackknife is DESCRIPTIVE, not exclusionary.** `few_family_contributing`
  tells us how broadly an inversion distributes across families — not
  whether it's "real". Only `pca_family_confounded` (PCA groups track
  ancestry, not arrangement) is a SUSPECT demotion.
- **4e's job is to pigeonhole each candidate** with a 14-axis tag. The
  evolution analysis that uses those pigeonholes is downstream of 4e,
  not in it. 4e is the finish line.

## Fixes applied this session

### FIX 39 — SILENT — jackknife vocabulary mismatch

**Files:** `4c_group_validation/STEP_C01f_hypothesis_tests.R` L499, L546

T9 (`test_ancestry_jackknife`) emits `"few_family_contributing"` for the
2–3-family case (both Engine B branch L1366 and sim_mat fallback L1455).
`compute_group_validation` matched on `"few_family"` at L499 and L546.
Every `few_family_contributing` verdict fell through the whole
classifier: `family_linkage` stayed `"unknown"`, `internal_cap` stayed
`NA`, no `"restricted_family_spread"` quality flag. Same candidate also
missed the SUPPORTED promotion at L546 on T8 concordance + family
restriction.

T9's internal "_contributing" suffix is used symmetrically for both 2-3
family (`few_family_contributing`) and ≥3 family (`multi_family_contributing`).
Decision tree at L1812 already pairs them correctly. Patch 03's
documented table dropped the suffix on the `few_family` row by mistake.

**Fix:** L499 `"few_family"` → `"few_family_contributing"` (plus long
comment explaining), L546 `jk %in% c("single_family_fragile", "few_family")`
→ `jk %in% c("single_family_fragile", "few_family_contributing")`.

### FIX 40 — ARCHITECTURAL — composite_flag key consolidation (Option A + B)

**Files:** `schemas/internal_ancestry_composition.schema.json` L69 and
`4b_group_proposal/STEP_C01i_d_seal.R` L359–362

Two writers produced the same semantic value under two different key names:
- Schema `keys_extracted` directive: `q1_ancestry_composite_flag`
  (from nested_composition block's `composite_flag` field)
- Seal explicit add_evidence: `q1_composite_flag`

Both values trace to the same source (`nst$data$composite_flag`), so no
value conflict. But it was noise — key-name duplication, maintenance cost.
No consumer of either key exists in the current tarball, but 4e will read
exactly one of them when it's wired.

**Fix (Option A+B per Quentin):**
- Schema rename: `q1_ancestry_composite_flag` → `q1_composite_flag` (B)
- Seal explicit write removed; schema's keys_extracted is now the sole
  writer (A)
- Long comment in seal documenting where the value now comes from

Chat-9 considered keeping seal's write as belt-and-suspenders for "v9-era
registries that have add_evidence but not write_block", but Quentin
correctly pointed out that's theoretical — since only one composite_flag
exists in the system, the defensive backup is wasted complexity.

### FIX 41 — CRASH (active) — cheat30 backward-compat alias scope

**File:** `4d_group_dependent/cheat30_gds_by_genotype.R` L66 (before) → L100 (after)

Alias `compute_pairwise_ibs <- compute_pairwise_gds` was placed INSIDE
the function body of `compute_pairwise_gds` (L66), making it a local-frame
assignment that vanished when the function returned. Never reached
module scope. Downstream caller at L222 (`run_cheat30` → `ibs <-
compute_pairwise_ibs(dosage_matrix)`) would crash with "could not find
function 'compute_pairwise_ibs'" on every call.

**Fix:** Moved the alias out of the function body and placed it at
module scope immediately after `compute_pairwise_gds`'s definition.
Long comment documenting the original mis-scope.

### FIX 42 — ARCHITECTURAL (large) — 4e build_key_spec aligned to v2 spec

**File:** `4e_final_classification/compute_candidate_status.R` L40–233 (rewritten)

Current `build_key_spec()` listed 315 keys. Authoritative v2 spec
(`specs/INVERSION_REGISTRY_SPECIFICATION_v2.md`) specifies 352. Gap:
~37 keys across Q2, Q3, Q5, Q6, Q7 from schema-extracted outputs that
were never added to the spec list.

Also real bugs in the Q6 spec:
- `q6_n_HOM_STD` (v9 name) vs `q6_n_HOM_REF` (v10 canonical) — drift.
- `q6_hwe_chisq_p` (spec) vs `q6_hwe_p` (what schema writes) — drift.
- Missing `q6_freq_class`, `q6_genotype_balance`, `q6_hwe_verdict`,
  `q6_hwe_chi2`, `q6_polymorphism_class` — all in
  frequency.v2.schema.json's keys_extracted but not in the 4e spec.

**Fix:** full rewrite of `build_key_spec` + `compute_completion`:

1. **Expanded to 367 keys** (slightly over 352 due to intentional v9+v10
   dual aliases: both `HOM_REF` and `HOM_STD`, both `jackknife_verdict`
   and `jackknife_status`). Numeric target per Q:
   - Q1: 49 (includes `q1_composite_flag`, `q1_ancestry_dominant_structure`)
   - Q2: 55 (schema extractions from internal_dynamics + recombinant_map +
     internal_ancestry_composition now explicit)
   - Q3: 74 (close to v2 target 81; minor carrier keys added)
   - Q4: 47 (+ 9 tandem-repeat keys from cheat28)
   - Q5: 39
   - Q6: 28 (fixed naming drift + added schema-extracted keys)
   - Q7 total: 75 (includes Q7B audit + new `q7_t9_jackknife_verdict`)

2. **Added `*_aspir` per-Q lists** marking keys whose writers don't exist
   yet (writers live upstream in phase 2/3 or in modules whose registry
   wiring is pending).

3. **`compute_completion()` reworked** to respect aspirational marking:
   - `applicable = total - aspirational - not_assessed`
   - `pct = resolved / applicable * 100` (meaningful today)
   - `pct_of_spec = resolved / total * 100` (shows gap to v2)
   - Per-Q `n_aspirational` field so you can see the scope of pending
     wiring per question
   - Aspirational keys silently excluded from the `top_missing` list
     so it doesn't spam the report with known gaps

### FIX 43 — SILENT (impactful) — jackknife_status key never written, readers fall back

**Files:**
- `4c_group_validation/group_validation_gate.R` L64–79 (readers)
- `4e_final_classification/characterize_candidate.R` L448–459 (same pattern)

`q7_t9_jackknife_status` has FOUR readers (gate.R, characterize_q7,
assign_confidence_tier line 299, build_key_spec which lists it) but ZERO
writers anywhere in the codebase. C01f produces a `t9_jackknife_verdict`
TSV column but that never becomes a flat registry key. So `jk_status`
always falls through to `"unknown"`, `jk_robust` is always FALSE, and
the gate falls back to OR-only logic — capping any candidate without
Layer D at UNCERTAIN regardless of jackknife robustness.

The semantic info IS reliably available via `q6_family_linkage` — seal
initializes it to "unknown" and C01f overwrites it to
`multi_family`/`few_family`/`single_family`/`pca_family_confounded`.

**Fix:** consumers now prefer `q6_family_linkage` and fall back to the
aspirational `q7_t9_jackknife_verdict`/`q7_t9_jackknife_status` only
for v9-era registries.
- Gate: `family_linkage == "multi_family"` → treat as `jk_robust`
- Gate: `family_linkage == "pca_family_confounded"` → treat as `jk_fragile`
  (canonical fragile signal, not single_family_fragile)
- Gate `has_groups` check extended to accept `q6_n_HOM_REF` (v10),
  `q6_n_HOM_STD` (v9), or `q6_n_total`

### FIX 44 — CRASH (dormant) — characterize_candidate missing source()

**File:** `4e_final_classification/characterize_candidate.R` header (L1–48)

`characterize_candidate()` at L517 calls `assess_group_validation(keys)`
and at L537 calls `check_group_gate(q, group_val)`. Both functions are
defined in `4c_group_validation/group_validation_gate.R`. Neither
was sourced — script would crash with "could not find function
'assess_group_validation'" on every invocation.

**Fix:** Added `.gate_file_candidates` locator at top with four path
candidates:
1. relative to sourced file (sys.frame-based)
2. `"../4c_group_validation/group_validation_gate.R"`
3. `"4c_group_validation/group_validation_gate.R"` (cwd)
4. `GROUP_VALIDATION_GATE` env var override

Tries each in order. Emits warning (not error) if none found —
acceptable because a driver can source gate.R explicitly before
`characterize_candidate.R` without triggering the warning path.

## Findings flagged (not fixed this chat)

- **Finding T (chat 8 carry-forward):** README structure_type drift
  (4 listed, 7 produced). Doc pass.
- **Finding U (chat 8 carry-forward):** `q2_decomp_quality_flags` forward-
  declared; no consumer. Still stands — seal writes it; nobody reads it.
  Will become relevant when 4d quality-aware decisions land.
- **Finding V (chat 8 carry-forward):** cheat24 threshold divergence
  (50/200kb vs 100/500kb). Design call needed.
- **Finding W (new this chat):** `register_C01f_keys` external helper
  lives in `utils/registry_key_helpers.R` outside tarball. Needs
  verification it handles the 5 new verdict columns
  (`group_validation_after`, `quality_flags`, `family_linkage`,
  `promotion_cap_applied`, `comp_source`).
- **Finding Y (new this chat):** 4d production launcher missing.
  `run_phase4.sh` expects `${MODULE_DIR}/launchers/LAUNCH_group_cheats.sh`;
  only `LAUNCH_group_cheats_example.sh` exists with stale paths pointing
  at non-existent `cheats/`, `scripts/`, `burden/` subdirs.
- **Finding Z (new this chat):** 4d README stale. Says "not yet populated"
  but directory contains 7 scripts. Listed scripts (cheat20, Q5, Q6) don't
  match actual (cheat27/28/29/30 + audit + regenotype).
- **Finding AA (new this chat):** `population_regenotype.py` + its SLURM
  wrapper are misfiled in 4d. Its own header says it's SV-caller pipeline
  step (MODULE_4D/4E STEP_A03b) — belongs upstream in the DELLY/Manta
  module, not in phase 4d.
- **Finding AB (new this chat):** 4e key-count section comments in
  `build_key_spec` were off by a few (49→46, 40→36, etc). Partially
  fixed by FIX 42's full rewrite; residual cosmetic drift possible.
- **Finding AC (resolved by FIX 42):** 4e Q6 spec vs frequency schema drift.
- **Finding AD (resolved by FIX 42 via dual alias):** HOM_STD vs HOM_REF
  naming.
- **Finding AE (resolved by FIX 43):** `q7_t9_jackknife_status` readers
  with no writer.
- **Finding AF (new, pending next-chat work):** `run_characterize.R`
  driver doesn't exist. characterize_candidate.R is a function library
  with no invocation path. FIX 44 made the functions *findable* but
  nothing calls them end-to-end. The driver is the main deliverable
  for chat 10.
- **Finding AG (new, pending next-chat work):** `SPEC_VS_REALITY.md` doc
  needed to enumerate which of the 367 keys have writers vs which are
  aspirational. The `*_aspir` lists in `build_key_spec` are the
  machine-readable version of this doc, but a human-facing doc helps
  prioritise phase-2/3 writer wiring work.

## Parse-check backlog (cumulative)

Eleven R scripts + one Python script need LANTA `Rscript -e "parse(file=...)"`
+ smoke test before production:

Chat 7 (six files):
- `phase_3_refine/03_statistical_tests_and_seeds.py` (AST-parse OK; pipeline smoke OK)
- `4c_group_validation/STEP_C01f_hypothesis_tests.R`
- `4c_group_validation/group_validation_gate.R`
- `4e_final_classification/compute_candidate_status.R` (earlier edits by chat 7)
- `4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
- `4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`

Chat 8 (two files):
- `4b_group_proposal/STEP_C01i_b_multi_recomb.R`
- `4b_group_proposal/lib_decompose_helpers.R`

Chat 9 (five files + one schema):
- `4c_group_validation/STEP_C01f_hypothesis_tests.R` (FIX 39 — re-parse since chat 7's edits)
- `4c_group_validation/group_validation_gate.R` (FIX 43)
- `4d_group_dependent/cheat30_gds_by_genotype.R` (FIX 41)
- `4e_final_classification/compute_candidate_status.R` (FIX 42)
- `4e_final_classification/characterize_candidate.R` (FIX 43, FIX 44)
- `4b_group_proposal/STEP_C01i_d_seal.R` (FIX 40)
- `schemas/internal_ancestry_composition.schema.json` (FIX 40 — JSON
  parse verified in container, but schema-validator should re-check
  the v2 schema registration)

Total: eleven R + one schema + one Python needing LANTA parse-check
before next pipeline run.

## Files modified this chat (seven)

```
inversion_modules/phase_4_postprocessing/
  4b_group_proposal/STEP_C01i_d_seal.R          (FIX 40)
  4c_group_validation/STEP_C01f_hypothesis_tests.R  (FIX 39)
  4c_group_validation/group_validation_gate.R   (FIX 43)
  4d_group_dependent/cheat30_gds_by_genotype.R  (FIX 41)
  4e_final_classification/characterize_candidate.R (FIX 43, FIX 44)
  4e_final_classification/compute_candidate_status.R (FIX 42)
  schemas/internal_ancestry_composition.schema.json (FIX 40 schema rename)
```

Plus three doc files (this audit log + session summary append +
FIXES_APPLIED append).

## Recommendation for chat 10

**Priority 1: Finish 4e.**
- Write `run_characterize.R` driver (Finding AF). ~150 lines. Reuses
  `compute_candidate_status.R`'s registry-key loader pattern. Sources
  `group_validation_gate.R` + `characterize_candidate.R`. Iterates
  candidates, calls `characterize_candidate(keys)` per candidate, writes
  `characterization.tsv` and per-candidate `characterization.txt`.
- Write `SPEC_VS_REALITY.md` (Finding AG). Enumerates 367 spec keys,
  marks aspirational per Q, groups by producer module. The phase-2/3
  writer wiring TODO list.

**Priority 2: The orchestrator break.**
- `run_phase4.sh` calls `LAUNCH_C01i_decomposition.sh` as one job but
  4b is four scripts (decompose + multi_recomb + nested_comp + seal).
  Either update run_phase4.sh to call `run_phase4b.sh`, or update the
  launcher to do the 4-script DAG.
- `LAUNCH_group_cheats_example.sh` has stale paths. Either promote to
  production with corrected paths or create a proper
  `launchers/LAUNCH_group_cheats.sh`.

**Priority 3: The 4a audit.**
- C01d pass-1, C01e, C01g — not yet audited end-to-end. Would close
  the loop on phase 4 entirely.
- Check Layer C (GHSL) read contract — what does C01d do when Clair3
  is incomplete on some chromosomes?

**Priority 4: the phase 2/3 writer wiring.**
- The `*_aspir` lists in FIX 42 are the TODO. Each aspirational key is
  a script that computes a value but doesn't register it as a flat key.
  This is significant work — probably 3–5 chats of its own. Not urgent
  UNTIL the pipeline has run at least once end-to-end.

## Budget note

Six fixes applied this chat, considerably over the nominal three-fix
ceiling. Quentin explicitly lifted the ceiling mid-chat ("if there's
any code or thing missing, you code it") to authorize the 4e rework.
Each fix is bug-backed and individually understood. FIX 42 is the
largest by line count (~200 lines changed) but surgical in intent:
the spec rewrite is mechanical once the design choice ("expand to 352
with aspirational flag") is settled.

Remaining budget would have gone to writing the driver + docs; chat
was stopped here rather than start-half-finishing those items.
