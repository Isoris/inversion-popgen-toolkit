# AUDIT LOG — chat 5 of the 4a audit series (2026-04-17) — PART B

Scope: phase 3 unpack + audit. Followed up on the Part A handoff's
structural-flatten + phase-3 audit remit. Additionally, resolved the
"phase 3 feeds nothing" finding that accumulated across chats 1–4 into
a proper wiring (phase 3 is the Layer D writer and the Layer B BND
rescue writer — it was never dead, just never wired through the
registry API).

## Summary

Nine findings this session: **six bug fixes** (FIX 23–28) applied to
phase 3 scripts, **one major architectural wiring** (FIX 29 v2) that
installs phase 3 as the writer for Layer D and supplementary Layer B
in the evidence registry, plus **one structural flatten** (MODULE_5A2
wrapper removed) and path-reference updates across eight documentation
files. All applied changes parse-check clean in-container (Python AST;
bash `-n`; JSON `json.load`). R files need an on-cluster `Rscript -e
"parse(file=...)"` — the container has no R installed.

One earlier mistake this session (a first-draft FIX 29 that added a
"D13" dimension to C01d) was explicitly **reverted** mid-session after
Quentin clarified the 4-layer architecture in which phase 3's OR test
is Layer D, not a sub-dimension of Layer A. The correct wiring — phase
3 writes registry blocks whose `keys_extracted` directives populate the
existing `q7_layer_d_*` / `q7b_bnd_*` flat keys that C01f and
`compute_candidate_status.R` already read — is FIX 29 v2.

## Structural change (done)

Flattened `inversion_modules/phase_3_refine/MODULE_5A2_breakpoint_validation/`
one level up. All 11 files now live directly in
`inversion_modules/phase_3_refine/`. The `MODULE_5A2_breakpoint_validation/`
wrapper directory was removed. No live script reference the old path;
eight documentation files were updated (listed in Part A handoff).
Seven retrospective mentions remain as intentional "was X, flattened
to Y" markers for historical context.

## Fixes applied

### FIX 23 — breakpoint_validator_standalone.py INV-vs-HET dead branch (DESIGN)

**File:** `phase_3_refine/breakpoint_validator_standalone.py` L888–896

**Problem:** The code
```python
or_het, p_het, ci_het = fisher_exact_2x2(a_h, b_h, c_h, d_h), None, None
if isinstance(or_het, tuple):
    or_het, p_het = or_het
    ci_het = fisher_exact_ci(a_h, b_h, c_h, d_h)
else:
    or_het_val = (a_h * d_h) / (b_h * c_h) if b_h * c_h > 0 else float('inf')
    _, p_het = fisher_exact_2x2(a_h, b_h, c_h, d_h)
    ci_het = fisher_exact_ci(a_h, b_h, c_h, d_h)
    or_het = or_het_val
```
tupled the 2-tuple return of `fisher_exact_2x2` into the first slot of
a 3-tuple, guaranteeing `isinstance(or_het, tuple)` was always True and
leaving the `else:` branch unreachable. Worse, the computed `or_het /
p_het / ci_het` values are **never consumed** — the `comparisons` list
at L898 contains only INV-vs-non-INV and INV-vs-REF; the docstring at
L880 promises three comparisons but only two are rendered.

**Fix:** simplified to straight-line code; added a TODO for wiring
INV-vs-HET into the `comparisons` list if the manuscript wants it.
Dead branch removed.

### FIX 24 — breakpoint_validator_standalone.py division by zero (CRASH)

**File:** `phase_3_refine/breakpoint_validator_standalone.py` L648

**Problem:** `f.write(f"  {g:<10} {y:>6} {n:>6} {t:>6} {y/t*100:>7.1f}%\n")`
crashes with `ZeroDivisionError` when any group (REF/HET/INV) is empty.
Every other division in the script is guarded with `if t > 0 else "–"`.

**Fix:** matched the guard pattern.

### FIX 25 — SLURM wrapper template filename mismatch (CRASH)

**File:** `phase_3_refine/breakpoint_validator_standalone.py` L1182

**Problem:** The auto-generated wrapper template references
`breakpoint_validator.py` (missing `_standalone` suffix). Users who run
the generated wrapper on the cluster hit "No such file or directory".

**Fix:** corrected filename to `breakpoint_validator_standalone.py`.

### FIX 26 — run_breakpoint_validation.sh STEP02/STEP05 wrong input filename (CRASH, critical)

**File:** `phase_3_refine/run_breakpoint_validation.sh` L45, L106

**Problem:** Both launcher calls pass
`'${BPV_EVIDENCE}/matched_delly_snake_candidates.tsv'` as `--candidates`,
but STEP01 writes the file as `matched_inv_candidates.tsv`
(01_extract_inv_candidates.sh L385). No symlink, no rename. **The
entire pipeline fails at STEP02 on first run.** Chat 4's "method is
architecturally complete, only plumbing remains" assessment missed this.

**Fix:** both launcher calls updated to `matched_inv_candidates.tsv`;
STEP02 and STEP05 `--help` strings aligned.

### FIX 27 — 03_statistical_tests_and_seeds.py seed_summary.tsv column drift (SILENT)

**File:** `phase_3_refine/03_statistical_tests_and_seeds.py` L222–228

**Problem:** `seed_summary.tsv` header is built from
`all_seeds[0].keys()`, but the `n_seeds_ref` key is conditionally added
only inside the `if qualifies:` branch. When qualification status varies
across candidates (some qualify, some don't), the TSV ends up with rows
of different widths under a header matched to whichever candidate went
first.

**Fix:** initialise `n_seeds_ref: 0` at dict creation alongside the
other seed-count defaults.

### FIX 28 — 05_delly_manta_concordance.py zero-Manta crash (CRASH)

**File:** `phase_3_refine/05_delly_manta_concordance.py` L138

**Problem:** Inside a `if n_delly_total > 0:` block, the code divides
by `n_manta_total` unconditionally. If DELLY finds calls but Manta
finds zero, the manuscript-sentence writer crashes.

**Fix:** split into four branches (both / delly-only / manta-only /
neither); updated `--help` string to reflect current STEP01 output name.

## FIX 29 v2 — PHASE 3 WIRED INTO PHASE 4 VIA EVIDENCE REGISTRY (architectural)

This is the heart of this session. Before-chat-5 state: phase 3 ran
end-to-end and produced flat TSVs (`all_candidates_tests.tsv`,
`seed_summary.tsv`, `bnd_inv_rescued_candidates.tsv`, etc.), but
**nothing in phase 4 read any of them**. The gate question raised in the
handoff ("does phase 3 feed 4a or 4d?") had a correct answer that made
everyone uncomfortable: phase 3 fed nothing; the `classify_inversions`
consumer didn't exist; C01i never touched the seed files; Layer B of
4a was sourced from C00 in phase 2c, not phase 3.

The correct resolution is **not** that phase 3 is dead — the user's
response-reply made this explicit. The resolution is that phase 3 owns
**Layer D** of the 4-layer evidence model entirely (the Fisher/Armitage
OR test linking PCA-INV carriers to physical breakpoint read evidence)
and supplements **Layer B** with BND rescue evidence (paired CT=3to3 +
CT=5to5 junctions that slipped past strict INV callers at low
coverage). Both Layer D and Layer B-rescue are already declared
consumers in phase 4 — C01f's `compute_group_validation()` gates
VALIDATED on `q7_layer_d_fisher_p < 0.05 AND q7_layer_d_fisher_or > 5`;
`compute_candidate_status.R` L188–190 expects the full
`q7_layer_d_*` key family. **The downstream was waiting for a writer.**

### The first-draft mistake (reverted)

My first-draft FIX 29 added a new "D13" scoring dimension to C01d that
read phase 3's flat TSV and combined concordance + inv_support_frac +
CA_pass into a 0–1 score, retuning composite weights and tier
thresholds. This confused Layer A (C01d's 12-dimension sim_mat-derived
composite) with Layer D (an independent evidence layer alongside A, B,
C). **Reverted in full.** C01d is restored to 12 dimensions with
original weights and `≥8 / ≥6 / ≥4` tier gates. The `--phase3_dir` CLI
flag was removed. A header comment documents the revert so nobody
retries the same mistake.

### The correct wiring (FIX 29 v2)

**Files changed:**

| File | Change |
|---|---|
| `phase_3_refine/03_statistical_tests_and_seeds.py` | Added `--registries_root` + `--candidate_map` CLI. Added `_load_registry()` helper that imports `registries/api/python/registry_loader.py` with graceful no-op if unavailable. Added `_load_candidate_map()` for inv_id → candidate_id translation. Per-candidate `reg.evidence.write_block(block_type="existence_layer_d", data={...})` emits a strict 2×2 HOM_INV-vs-HOM_REF contingency, Fisher OR + CI + p, Armitage Z + p, support counts, and sample lists per the schema's `required` fields. Closing summary prints block count. |
| `phase_3_refine/06_bnd_inversion_signal.py` | Added `--registries_root` + `--candidate_map` CLI. Added `_load_registry()` + `_load_candidate_map_bnd()` + `_match_cid()` (reciprocal-proximity lookup with 10kb tolerance per junction; falls back to synthetic `chrom:bp1-bp2:bnd_rescue` ID). Per paired junction, `reg.evidence.write_block(block_type="existence_layer_b_bnd_rescue", data={...})` emits orphan or catalog-matched pair info. Stats file now reports `registry_blocks_written`. |
| `phase_3_refine/00_breakpoint_validation_config.sh` | Added `REGISTRIES_ROOT` default (under `inversion_codebase_v8.5/registries`), and `CANDIDATE_MAP=""` stub. Both overridable from env. |
| `phase_3_refine/run_breakpoint_validation.sh` | Computes `REGISTRY_ARG` + `CANDMAP_ARG` before `sbatch` (the `--wrap=` heredoc is expanded by the submitting shell, not the worker). Passes both to STEP03 and STEP06. Echoes registry-enabled/standalone banner. |
| `registries/schemas/structured_block_schemas/existence_layer_b_bnd_rescue.schema.json` | New file. Co-exists with the primary `existence_layer_b` schema (owned by C00) so the BND-rescue contract doesn't mutate a C00-owned schema. 11 keys extracted: `q7b_bnd_rescued`, `q7b_bnd_rescue_source`, `q7b_bnd_pair_bp1/bp2/size_bp`, `q7b_bnd_match_type`, `q7b_bnd_matched_inv_id`, `q7b_bnd_left/right_pe`, `q7b_bnd_left/right_sr`. |
| `registries/schemas/structured_block_schemas/INDEX_remaining_blocks.json` | Added entry for `existence_layer_b_bnd_rescue`. |
| `phase_4_postprocessing/tests/test_registry_sanity.py` | Added fixtures for Layer D block (Fisher OR + Armitage + contingency) and BND-rescue block (cross-caller orphan pair). Asserts validation status on both. |

**Architectural invariants documented in code:**

- Schema-driven key extraction: the `keys_extracted` directives in each
  schema auto-materialise the flat keys. Writers do not hand-list keys;
  they hand the registry the block data and the machinery produces the
  flat `keys.tsv` append. This is existing registry behaviour, now
  actually exercised by phase 3.
- Silent standalone fallback: if `REGISTRIES_ROOT` is unset or points
  to a nonexistent directory, both writers print a one-line notice and
  keep running with flat-TSV-only output. Phase 3 remains useful
  outside the main pipeline.
- Contingency-table semantics: the Layer D schema requires a clean 2×2
  `HOM_INV` vs `HOM_REF` table, but phase 3's main Fisher test is
  INV-vs-non-INV (REF+HET pooled). The registry block emits the strict
  2×2 for the validation gate; the flat TSV keeps the INV-vs-non-INV
  Fisher as a sensitivity check. Armitage stays 3×2 (unambiguous
  trend) in both.

## Gate question — resolved

The handoff asked: *"Does phase 3 actually feed 4a, or only 4d? The
`phase_2_discovery/2c_precomp/README.md` claims `MODULE_5A2` feeds
C01g (phase 4a). Session notes from chat 1/2 said phase 3 feeds 4d via
`population_regenotype.py`."*

The resolved answer:

- **Neither 4a-direct nor 4d-via-regenotype was correct.** Phase 3
  feeds 4a and 4c (and 4e transitively) via the **evidence registry**,
  not via direct script reads of phase 3's flat TSVs. The registry is
  the indirection layer that decouples producers from consumers.
- 4d's `population_regenotype.py` is an **input** to phase 3 (produces
  dropout-rescued GT matrices that phase 3 can use for evidence
  extraction), not a consumer. The SLURM echo at
  `SLURM_A03b_population_regenotype.sh:156` ("use this as input for
  phase_3_refine") documents this direction.
- `classify_inversions` does not exist as a file; it's mentioned only
  in a single comment in a patches file. Document-only reference.
- `annotate_population_confidence.sh` living inside `phase_3_refine/`
  is a misfile — it sources `00_manta_config.sh`, writes into
  `MANTA_BASE/10_population_confidence/`, and is consumed by STEP_C00
  in phase 2c. Flagged; not moved this session (non-blocking cleanup).

## Found-but-deferred

- **C01i seeded init.** Phase 3's `04_deconvolution_seeds/{inv_id}_seeds.tsv`
  was designed as k-means initialization input for `STEP_C01i_decompose.R`
  (phase_3 README L42 plus the "How seeds feed into C01i" section).
  C01i currently does unsupervised k=3 on PC1. Wiring the seeded init
  is non-trivial — needs a `--seeds_dir` flag, per-candidate seed
  loading, and integration with whatever clustering call C01i makes —
  and is out of scope for the Part B audit. Flagged as future work.
- **`annotate_population_confidence.sh` relocation.** Should move out
  of `phase_3_refine/` to `MODULE_4H_ALL_Manta/` where it logically
  lives. Contract with STEP_C00 is already correct (column-mapping
  verified this session: `chrom, id, confidence_level, precision_class,
  pop_signal_score` all present in both writer and reader). The move
  is purely directory hygiene.
- **Layer B schema description cleanup.** The `existence_layer_b.schema.json`
  description still mentions "C00 build_flashlight" (old name;
  `build_sv_prior` is current). One-line edit; not touched this session
  because the schema is owned by 2c/C00 and was outside Part B scope.
- **D14 landscape classifier still orphaned** (FINDING 18 from chat 4).
  Unchanged.

## Parse-check

All modified Python files pass `python3 -c "import ast; ast.parse(...)"`.
Modified shell files pass `bash -n`. JSON files pass
`python3 -c "import json; json.load(...)"`. The modified R file
(`STEP_C01d_candidate_scoring_wired_25_v934_registry.R`) was restored
to its pre-FIX-29-draft state; the only net change is additive header
comment documenting the revert, which does not alter parse semantics.
Recommend `Rscript -e "parse(file=...)"` on LANTA before a production
run as a final sanity check.

## Revised understanding (for next chat's benefit)

Past chats oscillated between "phase 3 is the complete breakpoint
validation pipeline" (positive framing) and "phase 3 is dead, nobody
reads it" (negative framing). Both were partial. The truth is:

- Phase 3 **is** the complete breakpoint validation pipeline.
- Phase 3 **was** disconnected from phase 4 because nothing wrote its
  results into the registry.
- The registry schemas for Layer D and Layer B existed; the consumers
  in C01f and `compute_candidate_status.R` existed; the block type
  strings existed in the documented architecture. The only missing
  piece was the writer.
- Chat 5 wrote the writer. The 4-layer model is now actually 4 layers
  end-to-end.

If any phase 3 behaviour is questioned in future chats: first confirm
whether the question is about the flat-TSV side (STEP03 *did* write
those all along) or the registry side (STEP03 only writes those after
chat 5's FIX 29 v2). Both outputs now coexist.
