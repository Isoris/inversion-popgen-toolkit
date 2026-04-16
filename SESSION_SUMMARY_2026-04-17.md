# SESSION SUMMARY — 2026-04-17 phase 2 + phase 4a audit sweep

This session audited every script that feeds (or is fed by) phase
`4a_existence_layers/` end-to-end, and applied fixes to the real bugs
found. Output tarball: `inversion-popgen-toolkit_phase2_phase4a_fixes_2026-04-17.tar`.

## Files modified

| File | What changed |
|---|---|
| `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | Cheat 5 block rewritten (3 stacked bugs); `--flashlight_dir` silenced |
| `phase_4_postprocessing/4a_existence_layers/STEP_C01e_candidate_figures.R` | Panel F crash fixed; `compute_bands_for_candidate()` added for Panel B; `--repeats`/`--het_dir` silenced; metadata columns aligned with v9.3 C01d output; docstring corrected |
| `phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` | TE self-match bug fixed; `--ref_fasta` silenced; Cheat 17 fossil detection archived; `boundary_activity` simplified |
| `phase_4_postprocessing/4a_existence_layers/README.md` | Dead flags removed from C01d table; Cheat 17 archive section added |
| `phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R` | `cor()` short-intersect guard |
| `phase_2_discovery/2c_precomp/README.md` | Added 2026-04-17 CLI-flag cleanup notes + test_05 rename caveat |
| `phase_2_discovery/2d_candidate_detection/STEP_D05_ghsl_stability.R` | Header clarified — the unused file-based stub is now clearly labelled (sim_mat fallback is the real code path) |

## Files added

| Path | Purpose |
|---|---|
| `FIXES_APPLIED_2026-04-17.md` | Top-level fix log at repo root |
| `_archive_superseded/cheat17_fossil_detection/README.md` | Archive README explaining why Cheat 17 was removed + how to restore it if needed |

## The 14 fixes, grouped

### Crash-every-run (1)
- **BUG 5** — C01e Panel F: `gt_data$gt` referenced an undefined variable. `%||%` evaluates LHS first, so `gt_data$gt %||% 226` threw `object 'gt_data' not found` before the fallback could kick in. Fixed by deriving the sample count from `coseg_samples` and guarding `ref_list` with `exists()`.

### Silent data corruption (4)
- **BUG 9** — C01g TE fallback `grepl(chr, chr, fixed=TRUE)` always returned TRUE; on chromosome-name mismatches the fallback accepted the entire genome's TE annotations, inflating Cheat 21 verdicts. Fixed to match against `te_dt$chr`.
- **BUG 27** — C01d's Cheat 5 block had **three independent bugs stacked on the same ~20 lines**:
  (1) read `pc_obj$inv_likeness` but precomp stores the column inside `pc_obj$dt`;
  (2) read column `cheat5_family_fst_ratio` but C01a renamed it to `test05_family_fst_ratio`;
  (3) filtered on `mid_bp` which doesn't exist — only `start_bp`/`end_bp` do.
  Net effect: the family-vs-inversion F_ST ratio (a manuscript-relevant annotation) was always all-NA. Fixed by rewriting the whole block to read `pc_obj$dt`, accept both column names, and compute midpoints inline.
- **BUG 3** — C01e Panel B was plotting all 226 samples as grey. The `--triangles` flag read `triangle_sample_composition.tsv.gz`, a pre-v9.3 C01c artifact that no longer exists. Fixed by adding `compute_bands_for_candidate()` which runs k-means(3) on PC1 columns in precomp and returns REF/HET/INV assignments. Legacy file still accepted as override.
- **C01e metadata columns** — `candidate_metadata.tsv` was intersecting against pre-v9.3 dimension names (`d1_triangle..d10_ghsl`, `tube_stage_d/e`, `band1_n..band3_n`), so none of them survived the intersect and the per-candidate metadata file lost every dimension column silently. Fixed by aligning with current v9.3 C01d output: 12 dimensions, cheat25 viability, popgen annotations, morphology PA columns.

### Stale / dead flags (4)
- **BUG 1** — C01d `--flashlight_dir` parsed but never read. Silenced.
- **BUG 4** — C01e `--repeats` and `--het_dir` parsed but never read. Silenced.
- **BUG 6** — C01g `--ref_fasta` parsed but never read. Silenced.
- **ISSUE 7** — C01g `--scores` (Cheat 17 fossil detection) created a circular dependency with C01d that was never orchestrated as a two-pass run. Archived.

### Code quality / polish (5)
- **BUG 14** — PHASE_01C `cor()` guarded against < 3 samples shared. Prevents the R warning-and-NaN noise on short intersections.
- **C01e docstring** — claimed to produce `figure_composite.png`. It never did (composite assembly is external, via Inkscape or figrid). Docstring corrected to list actual outputs.
- **D05 header clarification** — earlier audit turn accused D05 of being a stub emitting all-NA partition stability. Re-reading showed the file contains TWO functions: the file-based `compute_ghsl_stability()` IS a stub, but `compute_ghsl_from_simmat()` is complete working code, and `run_all.R` calls the working one. Header now makes this unambiguous.
- **2c_precomp README** — noted the 2026-04-17 CLI-flag cleanup and the test_05 rename caveat (C01d output column is still `cheat5_family_fst_ratio` to avoid forcing a rename across 4b/4c/4e registry keys).
- **4a README** — removed the `--flashlight_dir` row from C01d's input contract table; added the Cheat 17 archive section.

### False alarm (1)
- **BUG 2** — claimed earlier that D05 emits all-NA partition_stability. Wrong. See the "code quality / polish" entry above.

## Verification

All five modified R files pass `Rscript -e "parse(...)"` — no syntax errors introduced. Full-repo parse-check finds only 6 errors, all in `phase_2_discovery/2c_precomp/patches/` (pre-existing legacy reference material, documented as such in its own README, never sourced at runtime).

Downstream column contract verified end-to-end:
- `phase_4_postprocessing/4e_final_classification/compute_candidate_status.R` reads `d11_boundary_concordance`, `d12_snake_concordance` via a name-to-key map at L581-600. C01d writes exactly these column names ✓
- `phase_4_postprocessing/tests/test_registry_sanity.py` references the same column names ✓
- 4b scripts (`STEP_C01i_decompose.R` et al.) read only `tier`, `pattern`, `final_score`, `start_mb`, `end_mb`, `chrom`, `interval_id` — all present ✓
- 4c `STEP_C01f_hypothesis_tests.R` reads `tier`, `pattern`, `final_score`, `start_mb`, `end_mb` — all present ✓

## What was NOT touched (by design)

- **Phase 3 (MODULE_5A2_breakpoint_validation)** — doesn't feed 4a directly. Feeds 4d via `population_regenotype.py`. The SV evidence that reaches C01g comes independently from C00 in 2c_precomp.
- **C01a internal math** — 1875 lines of PC1-band-geometry / inv_likeness / test_07 Beta / test_05 F_ST scan. I audited the output column contract (which is what everything else depends on); the internal maths is a science audit, not a plumbing audit.
- **Phase 2d deep-dive on D02/D07/D09/D13** — verified indirectly via the merged `scoring_table_<chr>.tsv` column contract being consumed by C01d. A direct audit of each D-module's output schema is worth doing in a future session but wasn't blocking.
- **The flashlight → sv_prior rename in the C01d output column name** — C01d writes `cheat5_family_fst_ratio` (legacy name) even though the underlying data comes from `test05_family_fst_ratio` in precomp. Renaming C01d's output would cascade to `compute_candidate_status.R` registry keys and `test_registry_sanity.py`; deferred to a future coordinated rename sweep.

## Remaining known issues (non-blocking)

- **BUG 19 / 20 / 21 / 22 / 23** from the C00 audit — robustness gaps in `STEP_C00_build_sv_prior.R` (dead `map_vcf_to_cga`, positional GT shortcut, INFO parser drops content after embedded `=`). None of these fire for the current catfish VCFs; left as-is.
- **Performance** — C00 is O(N²) on BND pair scan, C01g loads all 28 precomp RDS into memory. Fine for catfish scale; needs revisit if scaled to a larger species.
- **Cheat 17** — if the fossil vs active boundary distinction becomes manuscript-relevant, restore as a **post-hoc annotation script** that reads `boundary_catalog_unified.tsv.gz` + `candidate_scores.tsv.gz` after both exist. See `_archive_superseded/cheat17_fossil_detection/README.md` for how.

## For the next session

If running the pipeline:
1. The 14 fixes have been applied to the working copy. Do an actual LANTA test run on LG12 with the Tier-1 candidate at 52.4–54.6 Mb to confirm Panel B now plots three-band karyotypes and `cheat5_family_fst_ratio` populates for Tier 1/2 candidates.
2. Parse-check on LANTA with real R — a container-side parse-check is a necessary but not sufficient gate.

If continuing the audit:
1. **C01a internals** — worth a dedicated session, especially the test_05 F_ST scan at L1074-1248 (subsamples every 5th window and interpolates via `findInterval` — want to verify that's what's actually happening).
2. **Phase 2d D02/D09/D13** — verify each D-module's direct output schema against what Phase 8 in `run_all.R` actually reads.
3. **4c/4f hypothesis tests** — `STEP_C01f_hypothesis_tests.R` is 2189 lines and was flagged in the earlier session audit as needing the v10.1.1 `compute_group_validation()` return-type migration.
