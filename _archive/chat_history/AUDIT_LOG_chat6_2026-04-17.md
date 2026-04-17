# AUDIT LOG — chat 6 of the 4a audit series (2026-04-17) — PART C

Scope: phase 4 consumer-side audit. Follow-up to chat 5's FIX 29 v2
(write-side wiring of Layer D and BND-rescue). This session traced the
read side: does phase 4 actually consume the keys that chat 5 wired?

## Framing note

The chat began with a mid-session clarification from Quentin: **the
registry library itself is not yet implemented.** The specs and the
plan were set in an earlier chat, but the decision was to complete the
pipeline code and review it first, then code the lib and wire it up.
So throughout the current codebase, scripts forward-declare the keys
they will consume ("I'll need `q7_layer_a_detected`") and the write
side declares the keys it will emit ("I'll write `q7b_bnd_rescued` when
the lib is there"). Reader/writer parity and schema-derived
`*_detected`/`*_tested` flags are lib-era concerns, not current-state
bugs.

This reframing matters for Findings 3 and 4 below — they are
intentional forward-declarations, not orphan bugs. The audit log
captures them as lib-design notes for whoever writes the lib work,
not as bugs to fix.

Findings 1, 2, and 5 are real bugs independent of lib status and
should be fixed.

## Targets run

- (A) Trace `q7_layer_d_*` and `q7b_bnd_*` consumption in
      `compute_candidate_status.R` — done (Finding 2)
- (B) Trace `compute_group_validation()` call site for Layer D keys in
      C01f — done (Finding 1)
- (C) Check pathway gate consistency with schema — done (Finding 2)
- (D) Empirical smoke test: run STEP03 on fake strong-signal evidence
      with a scratch registry, inspect `keys.tsv` — done (confirms
      Findings 3/4 are intentional forward-declarations; also surfaced
      Finding 5)

## Findings

### Finding 1 — CRITICAL — C01f hands Layer D hard-coded NA

**File:** `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R` L2316–2325

The live call to `compute_group_validation()` passes:
```r
  layer_d_fisher_p      = NA,   # written by STEP03/Layer D, not C01f
  layer_d_fisher_or     = NA,
```
The function's VALIDATED-promotion gate (L525–526) then requires
`is.finite(layer_d_fisher_p) && is.finite(layer_d_fisher_or)`, which
short-circuits to FALSE. **The VALIDATED tier is unreachable
regardless of Layer D evidence.**

The patch file `patches/01_C01f_comp_from_registry.R` L282–283 shows
the same NA hard-code in the template caller, with a comment: "not
computed in C01f v9.3.4" / "STEP03 writes Layer D separately." The
patch author intended the caller to read from the registry but never
supplied the reader code; the live script adopted the template
verbatim.

**Fix** (pattern already used two lines above at L2302–2314 for
`q6_group_validation` and `q6_validation_promotion_cap`):

```r
  layer_d_fisher_p <- tryCatch({
    if (exists("reg") && !is.null(reg$evidence) &&
        !is.null(reg$evidence$get_evidence)) {
      ev <- reg$evidence$get_evidence(cid, "q7_layer_d_fisher_p")
      if (!is.null(ev) && is.data.frame(ev) && nrow(ev) > 0) {
        as.numeric(ev$value[1])
      } else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)

  layer_d_fisher_or <- tryCatch({
    if (exists("reg") && !is.null(reg$evidence) &&
        !is.null(reg$evidence$get_evidence)) {
      ev <- reg$evidence$get_evidence(cid, "q7_layer_d_fisher_or")
      if (!is.null(ev) && is.data.frame(ev) && nrow(ev) > 0) {
        as.numeric(ev$value[1])
      } else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)
```

Then pass these into `compute_group_validation(...)` in place of `NA`.
The NA fallback preserves current behaviour when the registry is
absent or when STEP03 hasn't run for a candidate.

**Pre-existing the lib question.** This fix works today (using the
current R registry_loader at `registries/api/R/registry_loader.R`
which already implements `get_evidence`) and continues to work after
the lib work lands — the signature `reg$evidence$get_evidence(cid,
key)` is stable.

Not applied this session because deliberate-patch discipline from
earlier chats says audit-only sessions flag, they don't patch. Recommend
applying as FIX 31 in a short follow-up.

### Finding 2 — HIGH — Pathway A / validation gates missing OR > 5

**File A:** `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R`
L298 (Pathway A full 4-layer) and L309 (3-layer partial) gate on
`fisher_p < 0.05` only. The Layer D schema's description field says
"OR>5, p<0.05" and C01f's `compute_group_validation()` at L525–526
gates on both. Once Finding 1 is fixed, 4e and C01f will disagree on
candidates with significant-but-low-OR Fisher results (e.g., p=0.04
with OR=1.2) — C01f will withhold VALIDATED, 4e will award Pathway A.

**Fix:**
```r
  # add alongside fisher_p at L272:
  fisher_or <- as.numeric(keys[["q7_layer_d_fisher_or"]] %||% 0)
  # and in the Pathway A gates at L298, L309:
  if (layer_a && layer_b && layer_c_strong && layer_d &&
      fisher_p < 0.05 && fisher_or > 5) { ... }
```

**File B:** `inversion_modules/phase_4_postprocessing/4c_group_validation/group_validation_gate.R`
L50–52 reads `q7_layer_d_fisher_p` only, no OR. `or_passed` should be
`or_p < 0.05 && or_val > 5`. One-line read + one-line gate edit.

### Finding 3 — NOTE (not a bug) — Reader-declared keys pending registry lib

**Files:** `4e_final_classification/compute_candidate_status.R` and
`4e_final_classification/characterize_candidate.R`.

These readers declare the key family `q7_layer_{a,b,c}_detected` +
`q7_layer_d_tested` and use them as gate flags for the pathway engine.
These keys are not emitted by any writer in the current codebase and
not extracted by any schema's `keys_extracted` directive.

**This is intentional.** The registry lib is specified but not yet
implemented; the scripts forward-declare the keys they will consume
once the lib lands. When the lib is written, the expected path is
that each layer schema gains a `keys_extracted` directive of the form
`{ "key": "q7_layer_d_tested", "from_derived": "block_exists" }` (or
equivalent) so the existence of the block materialises the `*_tested`
flag automatically. The current loader's `_extract_keys_from_schema`
only walks dotted paths into the data dict; a `from_derived` branch
needs to be added.

The same holds for Layer A/B/C detail keys (`q7_layer_b_delly`,
`q7_layer_b_manta`, `q7_layer_c_ghsl_quality`, etc.): the schemas at
`registries/schemas/structured_block_schemas/existence_layer_{a,b,c}.schema.json`
currently extract a different family of keys (Layer A extracts
`q1_d*`/`q7_tier`/etc., Layer B extracts `q7_n_sv_calls`, Layer C
extracts `q7_ghsl_contrast_score`/etc.). When the lib lands, those
schemas will need additional `keys_extracted` entries to populate the
flat keys the 4e readers already expect.

**Interim consequence to be aware of.** Until the lib work is done,
running the 4e pathway engine end-to-end on real data will put every
candidate into the terminal `D_insufficient` tier (or `D_artifact` if
the verdict string matches the artifact pattern), because all the
layer-detected flags default to FALSE. This is a known state; no fix
is needed if the intent is to keep 4e dormant until the lib is wired.

**No action for this chat.** Flagged for the lib-design session.

### Finding 4 — NOTE (not a bug) — BND-rescue keys pending consumer wiring

**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R`
L212. Two BND-rescue keys are declared (`q7b_bnd_rescued`,
`q7b_bnd_rescue_concordance`) in the spec block but never read into
decision logic. The other 9 keys extracted by
`existence_layer_b_bnd_rescue.schema.json` (`q7b_bnd_rescue_source`,
`q7b_bnd_pair_bp1/bp2/size_bp`, `q7b_bnd_match_type`,
`q7b_bnd_matched_inv_id`, `q7b_bnd_{left,right}_{pe,sr}`) have no
consumer anywhere in phase 4.

**This is intentional.** Chat 5's FIX 29 v2 wired the write side as a
forward-declaration of the BND-rescue contract. Readers are part of
the lib work. The natural read-side wiring, when it happens, is:
Pathway E's `layer_b` flag becomes TRUE when either a catalog INV
exists or `q7b_bnd_rescued == TRUE`; an additional pathway variant
("E_bnd_rescued") distinguishes rescue-only from strict-caller
candidates. Concrete design is deferred to the lib chat.

**No action for this chat.** Flagged for the lib-design session.

### Finding 5 — CRITICAL — STEP03 seed writer crashes on qualifying candidates

**File:** `inversion_modules/phase_3_refine/03_statistical_tests_and_seeds.py` L373, L375, L377

Inside the `if qualifies:` branch that writes the per-candidate seed
file for C01i, the code references `r['total_score']`:

```python
f.write(f"{r['sample']}\tINV_confirmed\t{r['total_score']}\n")
```

The evidence TSV schema produced by STEP02 has columns `sample, group,
support, pe, sr, …` (per the STEP03 reader at L197–211 and the fake
data smoke test this session). There is no `total_score` column. When
a candidate qualifies for seeding, `r['total_score']` raises
`KeyError`, crashing the script.

**This is a new bug adjacent to chat 5's FIX 27.** FIX 27 fixed the
*header* drift in `seed_summary.tsv` when qualifying status varied
across candidates; it did not exercise the qualifying branch past the
header. Highest-quality candidates — those with strong enough Fisher
signal and concordance to qualify for seeded deconvolution, which are
the ones the seeds are most useful for — crash the script.

**Empirical confirmation this session.** Ran STEP03 on fake evidence
(20 REF / 20 HET / 20 INV with 10% / 40% / 90% support fractions);
Fisher OR=27.0, p=1.66e-06, Armitage p=3.96e-07. Would have qualified
for seeding; forced non-qualification by setting absurd thresholds
(`--seed_min_concordance 0.999 --seed_min_support_frac 0.999`) to get
past the crash and reach the end of main(). With default thresholds
the candidate qualifies and the crash fires.

**Fix requires a column choice from Quentin.** Best candidates:
- `pe + sr` (paired-end + split-read count; already in the TSV)
- `concordance` (scalar computed per-candidate at L248; not per-sample)
- Compute a per-sample score by summing PE+SR and cast to int

Cleanest one-liner:
```python
score = int(r.get('pe', 0)) + int(r.get('sr', 0))
f.write(f"{r['sample']}\tINV_confirmed\t{score}\n")
```
but Quentin should confirm what the seed file is meant to carry — if
C01i's seeded-init design (still deferred) uses the score for weighted
k-means, a richer score may be wanted.

**Not applied this session** — needs author input on the column
choice. Flag as FIX 30 for a follow-up.

## Empirical smoke test summary

Built a scratch registry at `/home/claude/chat6_work/scratch_registry/`
with the real `api/` and `schemas/` trees, built a fake strong-signal
evidence TSV at `/home/claude/chat6_work/smoke/evidence/fake_candidate_01_evidence.tsv`
(60 samples, clean signal), ran STEP03 with
`--registries_root` pointing at the scratch registry.

**Verified:**
1. Registry discovery works when `<REGISTRIES_ROOT>/api/python/registry_loader.py`
   exists. Silent standalone fallback when the path is wrong (as chat 5
   designed — prints "[registry] API not found … Layer D block write
   skipped" and continues).
2. STEP03 writes `existence_layer_d` block to
   `<root>/data/evidence_registry/per_candidate/fake_candidate_01/structured/existence_layer_d.json`
   with `validation_status=validated` and 8 extracted keys.
3. The resulting `keys.tsv` contains exactly these 8 keys:
   `q7_layer_d_fisher_or, _fisher_p, _fisher_ci_lower, _fisher_ci_upper,
   _armitage_z, _armitage_p, _n_inv_with_support, _n_inv_total`.
   `q7_layer_d_tested` is not present. (Consistent with the
   forward-declaration framing above.)
4. Surfaced the Finding 5 crash when thresholds allow qualification.

Block write pathway is clean; the chat 5 FIX 29 v2 architectural wiring
works as designed on real invocations.

## What was NOT done this session

- Finding 1, 2, 5 patches — audit discipline, flagged not applied.
  Recommend applying as FIX 30 (Finding 5) + FIX 31 (Finding 1) +
  FIX 32 (Finding 2) in a short follow-up.
- `--generate_test_data` end-to-end from
  `breakpoint_validator_standalone.py` into STEP03 — the standalone
  script is a parallel pipeline that doesn't produce the `{inv_id}_evidence.tsv`
  files STEP03 consumes. The two are deliberately separate (standalone
  for one-off candidate QC; modular for the full pipeline). Smoke test
  used a handwritten fake evidence TSV instead.
- STEP06 BND-rescue smoke test — requires real DELLY+Manta VCFs; not
  tractable in-container without a lot of plumbing. The write-side
  schema and the Python code are static-checked; no reason to think
  the write pathway misbehaves.
- C01i seeded-init wiring, `annotate_population_confidence.sh`
  relocation, D14 landscape classifier — deferred from chats 4/5;
  still deferred.
- AUDIT_LOG narrative compression — the chat 5 session proposed
  shrinking audit logs; this one is intentionally kept on the short
  side (this doc is ~200 lines, vs chat 5's ~270). Still room to
  trim retrospectively if a refactor chat happens.

## For the next chat

Short follow-up session to apply the three concrete fixes:
1. FIX 30: STEP03 `r['total_score']` → `int(r.get('pe',0)) + int(r.get('sr',0))`
   (or Quentin's preferred column) at L373/L375/L377.
2. FIX 31: C01f L2322–2323 — replace two `NA` hard-codes with
   registry reads (pattern in Finding 1 above).
3. FIX 32: 4e `compute_candidate_status.R` L272 add `fisher_or` read
   + L298/L309 add `&& fisher_or > 5`; `4c/group_validation_gate.R`
   L50 add OR read + L52 add OR gate.

All three fixes are surgical, parse-check locally (container has no R,
so R files need LANTA `Rscript -e "parse(file=...)"`). Should fit
comfortably in one chat.

Findings 3 and 4 stay deferred until the registry lib chat happens.
When it does, use this audit log as the consumer-side specification —
the 4e readers' declarations already tell the lib exactly which flat
keys the pathway engine expects.

## Parse-check

No code changes this session. Only filesystem artefacts are this
audit log and a session-summary update.
