# Pass 6 — manuscript-aligned cohort + reference cleanup

**Date:** 2026-04-24
**Scope:** use manuscript v19 as source of truth to align cohort-identity and reference-genome wording across the whole repo. Supersedes passes 3 and 5 cohort-identity fixes with improved, manuscript-consistent phrasing. Also fixes the `MODULE_1_Methods_Revised.md` §1.2 narrative (the legitimate Mash species-QC step was framed wrong as "F₁ verification" when it's actually "pure-species contamination screening") and the `MODULE_4D/WIKI.md` manuscript blockquote.
**Risk:** zero on logic — docs and one R / one Python / one bash comment only. Python AST and bash -n verified clean. Two apparent "F₁" mentions remaining are intentional and describe the *reference assembly* or the *hatchery business*, not the cohort samples.

---

## My misread of pass 5 — apology + correction

In pass 5's DROP_README I wrote that `MODULE_1_Methods_Revised.md` §1.2's Mash-distance check was "wouldn't make sense for a pure-species cohort." **That was wrong.** The Mash check is entirely correct practice — it's a read-level species-identity QC against farm mixups, mislabelling with pure *C. macrocephalus* broodstock, or contamination with morphologically cryptic *Clarias* (notably *C. batrachus*). A Thai hatchery that *produces F₁ seed* holds both parental broodstocks on site, so verifying which species a sample actually is becomes essential, not optional. I flagged a legitimate QC step as a potential bug.

What was actually wrong in the §1.2 text is the narrative framing around the check, not the check itself. The old prose said "all 226 samples are known F₁ interspecific hybrids" and "F₁ hybrids are expected to show approximately equal Mash distance to both parental references" — backwards for what was actually done. The new prose in this pass makes the framing match reality: pure *C. gariepinus* expected, near-equal or Mac-closer distances are QC failures to follow up on.

---

## Canonical vocabulary (from manuscript v19)

Using the manuscript's own phrasing as source of truth:

- **Cohort**: "226 *C. gariepinus* broodstock individuals sampled from a commercial Thai hatchery that maintains broodstock for F₁ hybrid catfish seed production"
- **Reference**: "*C. gariepinus* (Gar) subgenome reference (`fClaHyb_Gar_LG.fa`, 28 pseudochromosomes, ~964 Mb — the *gariepinus* haplotype extracted from the haplotype-resolved F₁ hybrid assembly)"
- **Manuscript structure note**: Section 1 = F₁ hybrid assembly (where the reference comes from); Section 2 = pop structure + inbreeding + variant catalogs on the 226 cohort; Section 3 = inversions. The *F₁ hybrid* refers only to the assembly subject in Section 1; the 226 samples are pure *C. gariepinus*.

---

## The 12 files this tarball overwrites

**Top priority — critical manuscript-aligned fix**

1. `Modules/MODULE_1_read_prep/helpers/MODULE_1_Methods_Revised.md`  
   – §1.1 cohort identity line: "226 *C. gariepinus* × *C. macrocephalus* F₁ hybrid catfish" → "226 *C. gariepinus* broodstock individuals sampled from a commercial Thai hatchery that maintains broodstock for F₁ hybrid catfish seed production"  
   – §1.2 full narrative rewrite: the Mash QC check is preserved exactly as it was implemented, but the prose around it is now pure-species-cohort framing (detect contamination with Mac broodstock, cryptic *C. batrachus*, or farm mixups) rather than "verify F₁ hybrid status"  
   – §1.3 reference-genome phrasing: "*C. gariepinus* haplotype-resolved reference genome" → "*C. gariepinus* (Gar) subgenome reference ... extracted from the haplotype-resolved F₁ hybrid assembly described in Section 1 of the manuscript"

2. `Modules/MODULE_4D_delly_inv/WIKI.md`  
   – Line 228 (Manuscript-ready prose snippet for Methods): "haplotype-resolved hybrid reference (fClaHyb_Gar_LG.fa, 28 chromosomes)" → "*C. gariepinus* (Gar) subgenome reference (`fClaHyb_Gar_LG.fa`, 28 pseudochromosomes, ~964 Mb — the *gariepinus* haplotype extracted from the haplotype-resolved F₁ hybrid assembly)". Also "226 samples" → "226 *C. gariepinus* samples".  
   – Lines 24 and 196 (three-snake prose): retired to 4-layer framework vocabulary (re-application of pass 5).

**Bugs missed in earlier passes**

3. `Modules/MODULE_4C_delly_dup/WIKI.md:307` — "Catfish hybrids (Mac × Gar) are diploid F1" → "This study cohort is diploid *C. gariepinus* (2n = 56, Maneechot et al. 2016)"

4. `Modules/MODULE_CONSERVATION_CORE/STEP_00_setup.sh:276` — shell echo "Genome: fClaHyb_Gar_LG (Clarias gariepinus × C. macrocephalus F1 hybrid, Gar haplotype)" → "Genome: fClaHyb_Gar_LG (Clarias gariepinus Gar subgenome, extracted from the haplotype-resolved F1 hybrid assembly)". `bash -n` clean.

5. `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R:243` — R comment "may not be run on the hybrid genome" → "may not be run on the Gar subgenome reference"

6. `inversion_modules/phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md:154` — "RepeatMasker to be run on the hybrid genome" → "RepeatMasker to be run on the Gar subgenome reference" (pass 3 only caught the similar line at 170; this one at 154 was missed)

**Pass 3 files re-applied with manuscript-aligned phrasing**

7. `inversion_modules/README.md` — pipeline-root description now uses manuscript's canonical cohort + reference phrasing, includes explicit Mash-QC provenance pointer, and notes "20+ small families" as the hatchery structure (per your message).

8. `inversion_modules/phase_4_postprocessing/docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md` — K=8 rationale updated: "hatchery has 20+ small families as the biological structure; the K that would actually match the ancestry structure is ~20" (replacing the prior interim "~8 broodlines" wording with your actual figure). The scientific substance is unchanged: K=8 remains a scan-parameter choice balancing signal concentration × per-group sample size × composite-detection alphabet.

9. `inversion_modules/phase_4_postprocessing/4d_group_dependent/population_regenotype.py` — two docstring fixes re-applied (F1 hybrids violate HWE → hatchery broodstock with limited founders violates HWE). Python AST clean.

10. `inversion_modules/phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md:170-175` — TE/GO aspirational-keys paragraph uses manuscript-aligned "Gar subgenome reference (`fClaHyb_Gar_LG.fa`, the *C. gariepinus* haplotype extracted from Section 1's haplotype-resolved F₁ hybrid assembly)".

11. `registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md` — "226 F₁ hybrid *Clarias* fish" → "226 *Clarias gariepinus* broodstock individuals from a commercial Thai hatchery (verified as pure *C. gariepinus* by read-level Mash screening...)".

**Pass 5 files re-applied with manuscript-aligned phrasing**

12. `Modules/MODULE_2A_snp_discovery/README.md` — line 7 (three-snake → 4-layer framework) and line 154 (F₁ hybrids violate HWE → hatchery broodstock) re-applied.  
    `Modules/MODULE_4B_delly_del/WIKI.md` and `Modules/MODULE_4A_clair3_snp_indel/WIKI.md` — "catfish hybrid reference" → "*C. gariepinus* Gar subgenome reference (extracted from the haplotype-resolved F1 hybrid assembly in Section 1 of the manuscript)".

---

## Intentional mentions of "F₁" that ARE correct and stay in

Three places correctly refer to F₁ in a non-cohort sense. These are deliberate and grep will show them:

- `MODULE_1_Methods_Revised.md:5` — "maintains broodstock for F₁ hybrid catfish seed production" (describes the *hatchery's business*, not the cohort)
- `MODULE_1_Methods_Revised.md:9` — "consistent with an F₁ hybrid" (describes the *expected QC failure pattern* the screen is designed to catch)
- `MODULE_4D_delly_inv/WIKI.md:228` and 4 other files — "extracted from the haplotype-resolved F₁ hybrid assembly" (describes the *assembly* the reference was derived from — Section 1 of the manuscript)

---

## ⚠ Flag for manuscript review (not a repo fix)

While cross-checking the repo against manuscript v19, I noticed an internal inconsistency in the manuscript's own Methods section that you should look at:

**The MAPQ=60 justification talks about "homeologous regions between subgenomes"** ("haplotype-resolved reference retains extensive homeologous regions where reads could map ambiguously between subgenomes") — but **the manuscript also states reads were aligned only to the Gar subgenome** (28 pseudochromosomes, ~964 Mb). If the mapping reference has only the Gar subgenome, then it doesn't "retain both subgenomes as separate sequences" and the homeologous-mismapping concern would not apply in that way.

Either (a) the MAPQ=60 justification should reference within-subgenome repeat/paralogue mapping (not cross-subgenome mismapping), or (b) the mapping reference includes both subgenomes and the paragraph describing "aligned to the Gar subgenome reference" needs updating. `MODULE_1_Methods_Revised.md` §1.6 inherits this same inconsistency — I did not touch it in this pass to avoid diverging from the manuscript, but whichever you settle on for the manuscript, the MODULE_1 §1.6 text should follow.

---

## Manual steps after drop

Drag-drop, commit, push. No deletions needed — all 12 files are in-place overwrites.

---

## Suggested commit message

```
Align cohort + reference wording to manuscript v19 (supersedes passes 3/5)

Use manuscript v19 as source of truth:
- Cohort: "226 C. gariepinus broodstock individuals sampled from a
  commercial Thai hatchery that maintains broodstock for F1 hybrid
  catfish seed production"
- Reference: "C. gariepinus (Gar) subgenome reference (fClaHyb_Gar_LG.fa,
  28 pseudochromosomes, ~964 Mb — the gariepinus haplotype extracted
  from the haplotype-resolved F1 hybrid assembly)"

Critical fix: MODULE_1_read_prep/helpers/MODULE_1_Methods_Revised.md §1.2
narrative was wrong. The Mash-distance species-identity QC step was
implemented correctly but the prose around it framed it as F1-hybrid
verification. It's pure-species contamination screening (detect mixups
with Mac broodstock on site, or cryptic C. batrachus). Rewritten to
match the actual experiment and manuscript v19 §Methods.

Fix MODULE_4D/WIKI.md:228 manuscript-ready Methods blockquote to match
the manuscript's own canonical "Gar subgenome reference" phrasing.

Pick up 4 cohort/reference-identity bugs that passes 3+5 missed:
- MODULE_4C/WIKI.md:307 "Catfish hybrids (Mac × Gar) are diploid F1"
- MODULE_CONSERVATION_CORE/STEP_00_setup.sh:276 (echo message)
- compute_candidate_status.R:243 ("hybrid genome" R comment)
- SPEC_VS_REALITY.md:154 ("hybrid genome" — only line 170 was fixed
  in pass 3)

Update DESIGN_NOTE K=8 rationale with actual broodline count: the
hatchery has 20+ small families, so structure-match K is ~20. K=8
remains chosen for inversion-regime detection (scan parameter, not
structure match).

Python AST and bash -n verified.

Flagged but not changed in this pass: manuscript v19 has an internal
inconsistency between the MAPQ=60 justification ("between subgenomes")
and the stated mapping reference ("Gar subgenome, 28 chromosomes").
Noted for manuscript review, not a repo fix.
```
