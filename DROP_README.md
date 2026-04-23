# Cohort-identity fix — pass 3 (revised)

**Date:** 2026-04-24
**Scope:** correct 4 live docs + 1 code docstring that misidentified the 226-sample cohort as F1 hybrid when it's actually pure *Clarias gariepinus* hatchery broodstock. Additionally, rewrite the K=8 rationale in the design note because its original justification was scientifically wrong (not only because of the cohort conflation but also because the whole premise — that K matches cohort structure — was backwards).
**Risk:** zero on logic — only wording in comments/docs changes, no algorithmic change. Python syntax re-verified with `ast.parse()`.

**Supersedes the earlier pass-3 tarball** (which had an interim `[Quentin: verify...]` marker in the design note). Use this one instead. If you already dragged the earlier one, this one will cleanly overwrite it.

---

## Why this mattered

Your memory (and the architectural notes in `inversion_modules/_archive/chat_history/`) explicitly flag that three cohorts must never be conflated:

1. **F1 hybrid** (*C. gariepinus* × *C. macrocephalus*) — for the *genome assembly* paper
2. **226-sample pure *C. gariepinus* hatchery broodstock** — for `MS_Inversions_North_african_catfish` ← **the current project, this toolkit**
3. **Pure *C. macrocephalus* wild cohort** — for a future paper

5 places in the live tree incorrectly called cohort #2 an "F1 hybrid" cohort. That's the kind of error that would get caught in peer review of your inversion manuscript, or worse, propagate into someone else's citation chain.

A separate problem surfaced while fixing #3: the K=8 design-note rationale wasn't just using wrong cohort biology, it was using *wrong K-choice logic* — it claimed K=8 "matches the biology (~8 deep groupings)" as if K were a structure-match. Your actual logic is different: the hatchery has 20+ small families, so a structure-match K would be ~20, but K=8 is chosen as a **tool-choice for inversion-regime detection** (signal concentrates in 1-2 components, per-group sample size stays ≥25, composite-detection alphabet stays clean).

---

## What this tarball changes (overwrites the 5 files shown)

### 1. `inversion_modules/README.md` — the top of the pipeline-root README

**Before** (lines 3-5):
```
Consolidated deployment tree for the catfish F₁ hybrid inversion
pipeline (226 samples, *C. gariepinus × C. macrocephalus*). All paths
below are relative to this directory.
```

**After**: full paragraph describing the pure *C. gariepinus* hatchery cohort, mapping reference (the *gariepinus* haplotype from the F1 hybrid assembly, 963.9 Mb single haploid), and K clusters = hatchery broodline structure. Explicitly calls out that this is NOT the F1 hybrid cohort from the assembly paper.

### 2. `registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md`

Line 35: "226 F₁ hybrid *Clarias* fish" → "226 pure *Clarias gariepinus* hatchery broodstock fish".

### 3. `inversion_modules/phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`

**Full rewrite of Section 2 (lines 86-165)**, not just word replacement.

The old rationale said:
- "K should match the number of deep ancestry groups in the cohort"
- "~6 founder lines + two parent species → K ≈ 8 is biologically justified"
- "CONCLUSION: K=8 ... matches the biology (~8 deep groupings)"

Every one of these claims was wrong for this cohort. The new rationale says instead:
- This cohort is **NOT structure-matched by K**. The hatchery has 20+ small families, so a structure-match K would be ~20.
- K is used here as a **scan parameter for inversion-regime detection**, not a cohort-structure tool.
- The ancestry layer routinely sweeps K=2..12 (up to 20 for stability checks).
- K=8 is the default for per-candidate downstream steps because it balances three things:
  1. An inversion's dosage signal typically concentrates in 1-2 Q components at K=8 (clean, not fragmented)
  2. ~28 samples per component at K=8 = statistical tractability for Fst/jackknife (at K=20, 11/group is noise-dominated)
  3. 8-label alphabet for local `assigned_pop` keeps composite detection clean
- CONCLUSION explicitly reframes: "K=8 is the right default for inversion-regime detection — it is NOT chosen to match the hatchery's ancestry structure. Structure-matching would argue for K=20+ and would be the wrong choice for this pipeline's statistics."
- "WHEN TO DEVIATE" adds a new bullet for ancestry-structure questions (use the full K=2..12 sweep directly, don't rely on the K=8 default), and the "re-pick K by samples/group ≥25" guidance now correctly frames samples/group as the binding statistical constraint.

The IMPLEMENTATION paragraph at the bottom (make K a config variable via `ANCESTRY_K` env var) is preserved unchanged — that advice was always correct regardless of the cohort biology.

### 4. `inversion_modules/phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md`

Lines 170-175: "Require ... annotation on the hybrid genome ... may not be available for the F1 hybrid" → "Require ... annotation on the *C. gariepinus* reference (`fClaHyb_Gar_LG.fa`, the *gariepinus* haplotype extracted from the haplotype-resolved F1 hybrid assembly) ... coverage on this reference may be incomplete compared to long-established reference species."

The old wording muddled *reference* vs *cohort*. The live reference file is the 963.9 Mb single-haplotype *C. gariepinus* (not a merged 2-species genome); upstream `Modules/MODULE_2A_snp_discovery/docs/MODULE_2A_methods.md` already correctly describes it this way. This fix brings the phase-4 doc into line.

### 5. `inversion_modules/phase_4_postprocessing/4d_group_dependent/population_regenotype.py`

Docstring of `genotype_posteriors()`, lines 107-118. Two mentions of "F1 hybrids" changed to "hatchery broodstock with limited founders". The HWE-violation reasoning is unchanged (both cohort types violate HWE), just for the right biological reason now.

Python syntax re-verified with `python3 -c 'import ast; ast.parse(...)'` — parses cleanly.

---

## What this tarball does NOT change — and why

The `Modules/` tree (upstream read-prep, SNP discovery, DELLY SV calling) contains ~6 more "F1 hybrid" mentions. That tree may legitimately belong to the F1 hybrid *assembly* paper rather than the inversion paper. **Flag for pass 4:** if `Modules/MODULE_1_read_prep/helpers/MODULE_1_Methods_Revised.md` is the Methods section of your inversion manuscript, its cohort wording is much more visible and needs fixing too. If it's for the assembly paper, leave alone.

---

## Manual steps after drop

Just drag, commit, and push. No deletions needed — all 5 files are in-place overwrites.

---

## Commit message (suggested)

```
Fix cohort-identity conflation: 226-sample cohort is pure C. gariepinus

5 live files called the 226-sample inversion cohort "F1 hybrid" when
it's actually pure C. gariepinus hatchery broodstock. Additionally,
rewrite the K=8 design-note rationale whose logic was structurally
wrong (it framed K as a cohort-structure match; actually K is a
tuning parameter for inversion-regime detection, and the cohort has
20+ small families which would argue for K~20 if structure-matching
were the goal).

Files corrected:
- inversion_modules/README.md (pipeline-root description + reference
  genome clarification)
- registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md
- phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md
  (full rewrite of Section 2: K=8 framed as scan-parameter choice
  optimizing signal concentration, per-group sample size, and
  composite-detection alphabet — NOT a structure match)
- phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md
  (also corrects "hybrid genome" → "C. gariepinus reference" — the
  live reference file is 963.9 Mb single haplotype)
- phase_4_postprocessing/4d_group_dependent/population_regenotype.py
  (docstring: HWE violation attributed to hatchery-broodstock founder
  limits, not to F1 interspecies cross)

References to "F1 hybrid assembly" / "F1 hybrid cohort for the assembly
paper" are retained where they correctly disambiguate — the reference
genome really IS derived from an F1 hybrid assembly.

Upstream Modules/ tree (MODULE_1-4) not touched; its cohort-identity
wording may be for a different project (assembly paper) and needs
human judgment.

Python syntax verified with ast.parse; no logic change.
```
