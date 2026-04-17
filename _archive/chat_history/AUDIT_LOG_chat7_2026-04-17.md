# AUDIT LOG — chat 7 of the 4a audit series (2026-04-17) — PART D

Scope: apply the three chat 6 findings as patches, then begin the real
phase 4a audit. Framing: registry lib is not yet implemented; reader-
side orphan keys are forward-declarations pending the lib chat (per
chat 6 Findings 3/4).

## Fixes applied this session

### FIX 30 — CRASH — STEP03 seed writer KeyError

**File:** `inversion_modules/phase_3_refine/03_statistical_tests_and_seeds.py` L373–377

`r['total_score']` referenced a column the evidence TSV doesn't have.
Replaced with a defensive per-sample evidence score from the columns
that ARE in the TSV: `int(r.get('pe',0) or 0) + int(r.get('sr',0) or 0)`.
Extracted to a local `_ev_score(r)` helper with try/except fallback to
0 on malformed rows.

**Verified:** rebuilt the chat-6 scratch registry, reran STEP03 with
default thresholds. Strong-signal candidate qualifies (previously
crashed); seed file written cleanly: REF rows evidence_score=0,
HET rows =2 (pe=1+sr=1), INV rows =5 (pe=3+sr=2). Python AST-parse OK.

### FIX 31 — CRITICAL — C01f Layer D NA hard-code → registry read

**File:** `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R` L2316–2325

Replaced `layer_d_fisher_p = NA, layer_d_fisher_or = NA` in the
`compute_group_validation()` call with two `tryCatch` reads from
`reg$evidence$get_evidence(cid, key)`. Pattern mirrors the
`q6_group_validation` and `q6_validation_promotion_cap` reads two
lines above in the same file. NA fallback preserves behaviour when
STEP03 hasn't run for a candidate or registry is absent.

Also updated the patch template at
`patches/01_C01f_comp_from_registry.R` L276–304 to show the correct
reader instead of the NA hard-code, so future patch-reapplies don't
re-introduce the bug.

Brace balance verified (285 matched). Needs LANTA Rscript parse-check.

### FIX 32 — HIGH — Pathway A / validation gate missing OR > 5

**Files:**
- `4e_final_classification/compute_candidate_status.R` L272, L298, L309
- `4c_group_validation/group_validation_gate.R` L50–52, L89, L96

Added `fisher_or` reader (default 0 = conservative = blocks gate when
missing) and `&& fisher_or > 5` clause to both Pathway A gates in 4e
(4-layer convergence + 3-layer partial) and to `or_passed` in
`group_validation_gate.R`. Reason strings updated from `OR_p=...` to
`OR=..., p=...`. Now consistent with the Layer D schema description
and C01f's `compute_group_validation()` (after FIX 31).

Brace balance verified in both files. Needs LANTA Rscript parse-check.

## Phase 4a audit findings (partial — got through C01d only)

C01d (`STEP_C01d_candidate_scoring_wired_25_v934_registry.R`) —
1049 lines post-patches. Read order: header comment → CLI parsing →
per-candidate scoring loop → D11/D12 post-loop fill → final write.

### FIX 33 — SILENT — tier not recomputed after D11/D12 fill (CRITICAL)

**File:** C01d L819–832 (pre-patch); new block at L826–840 post-patch

The scoring loop at L348–376 computes `final_score`, `dim_positive`,
and `tier` inside the per-candidate loop. At that point `d11 = 0` and
`d12 = 0` (defaults assigned at L339 / L346 before the loop). After
the loop, the script reads C01g's boundary catalog and C01b's seeded
regions and fills `d11_boundary_concordance` / `d12_snake_concordance`
in place (L555–618, L467–546).

Then at L819–832 the script **recomputes `final_score` and
`dim_positive` with the real D11/D12**, but **never recomputes `tier`.**

Result: every candidate whose D11 or D12 crosses the 0.30 threshold
has `dim_positive` increased by 1 or 2 after the recompute — but
`tier` in `cand_dt` still reflects the value assigned in the loop
using `d11 = d12 = 0`. A candidate with 8 positive dimensions
(including both D11 and D12) could have been assigned Tier 2 in the
loop (6 positive dimensions at that point), and the master table
writes `dim_positive = 8, tier = 2` — internally inconsistent.
Downstream 4b/4c/4d/4e readers that use `tier` directly (without
re-checking `dim_positive`) inherit the deflated tier.

The peel-disappeared downgrade at L378 also runs with the pre-D11/D12
tier.

**Fix (applied):** after the `dim_positive` recompute at L827–832,
added a `tier` recomputation using the same ≥8/≥6/≥4 ladder as L373,
then re-applied the peel-disappeared downgrade guarded by column
existence. Both new blocks carry the same thresholds as the original
in-loop version so behaviour is identical to "what would have
happened if D11/D12 had been known in the loop."

### FIX 34 — cosmetic — stale dimension/tier header comment

**File:** C01d L32–51 header, L412 inline comment

Header described 10 dimensions with `Tier 1: ≥7/10` thresholds, but
live code is 12 dimensions with `≥8/≥6/≥4`. Updated header listing to
include D11 + D12 and corrected tier thresholds. Fixed the inline
comment at the output-table construction ("10 scoring dimensions" →
"12 scoring dimensions (D11+D12 filled after loop)").

### Finding 8 — NOTE (pending lib) — 4a writes no registry blocks

**Files:** C01d, C01e, C01g

Searched for `write_block`, `add_evidence`, `reg\$` patterns across
all three 4a scripts. Results:

- **C01d:** zero registry writes. Produces the flat
  `candidate_scores.tsv.gz` catalog; no block writes, no flat-key
  writes. README L11 claims "C01d writes Layer A block" — not accurate
  under current code.
- **C01g:** attempts registry writes at L1397–1426 via
  `register_C01g_boundary()` and `store_C01g_boundary()` helpers that
  are sourced from `utils/registry_key_helpers.R`. **That file does
  not exist anywhere in the working tree.** The `source()` call is
  tryCatch'd, the `exists(..., mode = "function")` guards silently
  return FALSE, the loop body never executes, and C01g prints
  `"[C01g] Registered 0 boundary evidence keys"` — which looks like a
  legitimate zero-result rather than an infrastructure miss.
- **C01e:** zero registry writes. Diagnostic figures only.

The 4a README's "Contract with phase 4b" section (L123–136) promises
`q1_layer_{a,b,c}_pass`, `q1_existence_score`, and four block files
(`existence_layer_{a,b,c}.json`, `boundary.json`). **None of the four
blocks is written by any 4a script today.** The `q1_layer_*` keys are
forward-declarations.

**This is consistent with the lib-not-yet-implemented framing.** When
the lib chat happens, the natural wiring is:
- C01d writes `existence_layer_a` block per candidate with data fields
  matching the existing schema's `keys_extracted` directives (which
  already extract `q1_d01..d12`, `q1_dim_positive`, `q1_shape_class`,
  `q1_composite_score`, `q7_tier`, etc. — see chat 6
  inventory).
- Layer B and Layer C blocks need writers in phase 2c (C00) and phase
  2e respectively, not in phase 4a — the 4a README's L11 table is the
  accurate one; the L123–136 "contract" section is overambitious.
- C01g writes `boundary` blocks (left + right per candidate) via the
  helper that needs to be created in `utils/registry_key_helpers.R`.

**Action this chat:** none. Flagged for the lib-design session. The
one change I would propose for the lib chat is making the C01g silent
failure mode noisy: when `.bridge_available` is true but the helper
file is missing, print a WARNING instead of letting the `exists()`
guard silently skip the body. That way a real dataset run surfaces
the missing helper immediately.

### Finding 9 — NOTE (pending lib) — silent skip on missing helper

Part of Finding 8 but worth calling out separately: the pattern
```r
for (.hf in c("utils/registry_key_helpers.R",
              "../utils/registry_key_helpers.R")) {
  if (file.exists(.hf)) { source(.hf); break }
}
if (exists("register_C01g_boundary", mode = "function")) { ... }
```
at C01g L1398–1417 silently does nothing when the helper is missing.
This is the same silent-standalone-fallback pattern chat 5 designed
for STEP03's optional registry — but in STEP03 the fallback is
explicit (`print("[registry] API not found ... skipped")`); in C01g
there's no fallback notice. When the lib chat happens, consider
adding a `WARNING` message before the `exists()` guard to make the
skip visible when `.bridge_available` is TRUE.

## Phase 4a audit — NOT yet covered

Ran out of budget after the C01d pass. Still to do:

- **C01g (1428 lines).** Verify: five source inputs per README
  (landscape, staircase, seeded regions, SV prior, blue-cross inner)
  are all read; fossil/Cheat 17 is truly archived (README L78–94
  claim); removed flags `--ref_fasta` and `--scores` are truly gone;
  `boundary_catalog_unified.tsv.gz` columns match what C01d's D11
  reader at L580–590 expects (`boundary_verdict`, `boundary_bp`,
  `matched_core_id`, `side`, `n_cheats_supporting`, `chrom`).
- **C01e (694 lines).** Lower priority — figures only, doesn't feed
  downstream. Verify no dead flashlight-style flags, and that it
  reads C01d's current column schema (d1..d12 + `dim_positive` +
  `tier` + `pattern`).
- **Cross-script contract end-to-end check.** The reader at C01d
  L580–590 assumes specific column names in C01g's output
  (`boundary_verdict`, `boundary_bp`, `matched_core_id`, `side`,
  `n_cheats_supporting`). Grep C01g for those column writes.
- **Validation that FIX 33's tier recompute doesn't break C01e
  figure generation** — if C01e keys figures off `tier`, the now-
  higher tier counts change the figure distribution. Not a bug, just
  worth noting.

## For chat 8

Start with the continuation list above: finish C01g (big, hour-ish),
scan C01e (fast), then end-to-end contract check. Close with an
honest assessment: can 4a actually produce the catalog in a state
phase 4b/4c/4d/4e can consume? Current answer under the lib framing
is "yes for flat `candidate_scores.tsv.gz`, no for registry blocks
until lib lands." That's fine — just make sure the catalog contract
is right.

## Parse-check status

- Python (STEP03 FIX 30): AST-parse OK, empirical smoke test
  confirmed.
- R (FIX 31, 32, 33, 34): brace + paren balance verified in all four
  edited files. Container has no R; LANTA
  `Rscript -e "parse(file=...)"` required before production.

Files edited this chat (six):
```
inversion_modules/phase_3_refine/03_statistical_tests_and_seeds.py
inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R
inversion_modules/phase_4_postprocessing/4c_group_validation/group_validation_gate.R
inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R
inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R
inversion_modules/phase_4_postprocessing/patches/01_C01f_comp_from_registry.R
```

Plus the three log/summary files (this audit log, FIXES_APPLIED
append, SESSION_SUMMARY append).

## Post-initial-writeup: C01g + C01e audit completed (same chat)

Finished the C01g + C01e audit in this same chat after clarifying
with Quentin that audit-only is not a rule — fix-as-we-go is the
posture unless a finding needs a design call.

### FIX 35 applied — SILENT — C01g dedup missing matched_core_id

**File:** C01g L465–481 (pre-patch) / now includes matched_core_id

The dedup aggregation at L465–481 produced `candidate_ids` as a
semicolon-joined string (all cluster candidates) but never produced a
single-value `matched_core_id` column. The registry-write block at
L1407 guards on `"matched_core_id" %in% names(confirmed)` — which was
always FALSE — so the per-boundary registry write loop would be
skipped entirely, even if the helper file from Finding 9 existed.

**Fix:** compute `matched_core_id` as the candidate_id of the
dedup's "best_idx" row (highest-priority boundary type in the
cluster), falling back to the first non-NA cluster candidate,
falling back to NA. `candidate_ids` (the joined list) is kept for
traceability.

Net effect: when the lib work from Finding 9 adds
`utils/registry_key_helpers.R`, the registry-write guard will now
succeed and per-boundary blocks will actually be written. Without
FIX 35, fixing Finding 9 alone wouldn't have been enough.

### Audit outcomes — C01g

- **Five source readers** (L212–438): all wired correctly; legacy
  filename fallbacks in place for pre-rename output dirs (snake1
  cores, sv_flashlight); `candidate_id` is populated from the
  source-specific ID column in each reader.
- **Dead flags** `--ref_fasta` and `--scores` (L106, L108): properly
  accepted-and-ignored per the chat-4 archival notes.
- **Cheat 17 archival** (L1197–1217): properly archived; `cheat17_*`
  columns emitted as NA for downstream schema stability;
  `boundary_activity` hard-coded to `"active"` with a clear comment.
- **Summary column `n_fossil`** (L1394): cosmetic leftover —
  always 0 because `boundary_activity` is always "active". Not
  worth a fix; the README archive note covers it.
- **Minor: `bp_pos == inv$bp1` float equality** (L400, L426):
  brittle for degenerate inversions (bp1 == bp2 → both positions
  labeled "left"). Extremely low probability; noting only.
- **Contract with C01d's D11 reader** (C01d L570–607 vs C01g
  output): all columns C01d reads (`boundary_bp`, `boundary_verdict`,
  `n_cheats_supporting`, `chrom`) are produced by C01g. Matched
  cleanly. The `matched_core_id` column C01d doesn't read is the
  one FIX 35 added — needed for the registry contract, not the D11
  scoring contract.

### Audit outcomes — C01e (694 lines)

- **Dead flags `--repeats` and `--het_dir`** (L82–83): properly
  accepted-and-ignored per chat-3 notes.
- **Triangle composition fallback** (L117–170): chat 3 already
  fixed the "everyone is HET" bug by replacing the retired
  `triangle_sample_composition.tsv.gz` reader with an on-demand
  k-means-on-PC1 compute_bands helper. Logic is clean.
- **C01d catalog consumption** (L107–113): reads only `interval_id`,
  `tier`, `final_score` — all still present in the C01d output after
  FIX 33/34. No breakage. FIX 33 means C01e now plots the correct
  tiered set (more candidates may land at Tier 1 than before);
  behavior change is correct.
- **No registry writes**: matches Finding 8. C01e is figures-only
  and doesn't need to write registry state.

### End-to-end contract check — phase 4a output consumption

Live readers of `candidate_scores.tsv.gz`:
- `4c/STEP_C01f_hypothesis_tests.R` L134 (fread, filters on tier)
- `4e/compute_candidate_status.R` L539 (fread from multiple
  fallback paths)
- `orchestrator/run_phase4b.sh` L37 (path export)

All three read `tier` and `final_score`. After FIX 33, these values
are correct (D11/D12-aware) for the first time.

No live reader of C01g's `boundary_catalog_unified.tsv.gz` exists
outside of C01d's D11 reader (verified: grep of `boundary_catalog_unified`
across the tree returns only C01d L557 + doc mentions). Under the
lib-not-yet-implemented framing this is expected — the registry
blocks that phase 4b/c/d/e would consume don't exist yet.

### Verdict on 4a readiness

- **Flat-TSV catalog**: C01d produces `candidate_scores.tsv.gz` with
  correct tier (post FIX 33), 12-dim scoring, and all downstream-
  expected columns. Phase 4b/c/d/e can consume it today.
- **Unified boundary catalog**: C01g produces
  `boundary_catalog_unified.tsv.gz` with all columns C01d's D11
  reader expects. D11 wiring works today.
- **Registry blocks from 4a**: not written by any 4a script. This
  is intentional (Finding 8 = lib-pending). When the lib chat
  happens, three things need to land together:
  1. `utils/registry_key_helpers.R` with
     `register_C01g_boundary` and `store_C01g_boundary` functions.
  2. Layer A block writer in C01d (the schema already exists at
     `registries/schemas/structured_block_schemas/existence_layer_a.schema.json`).
  3. Loader extension for `{from_derived: "block_exists"}` to
     synthesize the `q7_layer_{a,b,c,d}_{detected,tested}` flags
     from block existence (the 4e pathway gates depend on this —
     chat 6 Finding 3).

With FIX 35, when items (1) and (2) land, the C01g registry-write
block will now actually execute.

## Updated list of files edited (chat 7, final)

```
inversion_modules/phase_3_refine/03_statistical_tests_and_seeds.py      (FIX 30)
inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R   (FIX 31)
inversion_modules/phase_4_postprocessing/4c_group_validation/group_validation_gate.R        (FIX 32)
inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R (FIX 32)
inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R (FIX 33 + FIX 34)
inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R (FIX 35)
inversion_modules/phase_4_postprocessing/patches/01_C01f_comp_from_registry.R               (template update for FIX 31)
```

Six source files + one patch file + three doc files (audit log,
session summary, fixes applied).
