# SESSION SUMMARY — 2026-04-17 phase 2 + phase 4a audit sweep

This session audited every script that feeds (or is fed by) phase
`4a_existence_layers/` end-to-end, and applied fixes to the real bugs
found. Output tarball: `inversion-popgen-toolkit_phase2_phase4a_fixes_2026-04-17.tar`.

**Update 2026-04-17 (chat 4 — Part A continuation):** Added FIX 15,
16, 17, 20 and FINDING 18, 19. See the corresponding section at the
bottom of this file and `AUDIT_LOG_chat4_2026-04-17.md` for details.

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

## Chat 4 additions (Part A — 2d/2e/C01f)

### Files modified (chat 4)

| File | What changed |
|---|---|
| `phase_2_discovery/2d_candidate_detection/run_all.R` | FIX 15 (drop D01 contrast before merge) + FIX 16 (reciprocal-overlap merge of NN sweep tree nn_birth onto scoring table) |
| `phase_2_discovery/2d_candidate_detection/STEP_D04_interior_cv.R` | FIX 17 (rename `patchiness` → `below_median_frac` to unblock D08's CV column) |
| `phase_2_discovery/2d_candidate_detection/STEP_D14_landscape_classifier.R` | FINDING 18 header documenting orphan status |
| `phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R` | FIX 20 (on-the-fly k-means fallback ported from C01e; three-tier fallback chain) |

### Files added (chat 4)

| Path | Purpose |
|---|---|
| `AUDIT_LOG_chat4_2026-04-17.md` | Detailed findings from chat 4 |
| `HANDOFF_PROMPT_next_chat_partB_2026-04-17.md` | Part B prompt with this chat's findings folded in |
| `phase_2_discovery/2e_ghsl/README.md` | Documents that 2e is standalone-figure-only (FINDING 19), not wired to C01d scoring |

### Chat 4 fixes, grouped by severity

**Silent data corruption (FIX 15, 17):**
- FIX 15 — D01/D08 `contrast` column collision. `merge()` produced
  `contrast.x`/`contrast.y`, C01d's `iv$contrast` was NULL,
  `d1_contrast = 0`, D1_block_strength silently collapsed on every
  candidate since v9.3. Reproduced in a fresh R session.
- FIX 17 — D04 `patchiness` (tail fraction) collided with D08
  `patchiness` (CV) by name. Merge produced NULL, C01d L222
  `1 - patchiness * 3` saturated patch_score at 1.0 (D5 always
  inflated). Renamed D04's column to `below_median_frac`.

**Design gap (FIX 16):**
- FIX 16 — NN sweep tree was built and saved but never merged onto
  the scoring table. C01d's D3 reads `iv$nn_birth` — always NA
  before this fix, so D3 capped at 0.6 instead of 1.0 on real
  inversions. Added reciprocal-overlap join.

**Silent candidate drop-out (FIX 20):**
- FIX 20 — C01f's `build_comp_with_fallback` fell back to an
  obsolete triangle file when the C01i registry was empty.
  Candidates with no registry entry were skipped silently. Ported
  C01e's `compute_bands_for_candidate()` as
  `comp_from_kmeans_fallback()`; new three-tier chain: registry →
  k-means → legacy file.

**Documented-not-fixed (FINDING 18, 19):**
- FINDING 18 — D14 landscape classifier is orphaned from run_all.R.
  C01d's pattern-classification columns (`category`,
  `landscape_category`, `n_children`) are absent from the merged
  table; C01d correctly falls through to its shape_class-based
  fallback. Not a bug, dead contract. Header comment added.
- FINDING 19 — 2e/C04 snake3 GHSL v5 is not wired to the C01d
  scoring path. Outputs consumed only by two phase_5_followup
  supplementary-figure scripts. Layer-C integration is a future
  project. `phase_2_discovery/2e_ghsl/README.md` written.

### Downstream impact of chat 4 fixes

These fixes affect the C01d-produced `candidate_scores.tsv.gz` and
the C01f hypothesis test verdicts. They do NOT affect:
- C01g boundary catalog (reads only `start_bp, end_bp, block_id`
  from `scoring_table_<chr>.tsv`)
- C00 sv_prior build
- C01a internals
- phase 3 / MODULE_5A2_breakpoint_validation
- 4b, 4c except C01f, 4d, 4e column contracts (all verified in
  chat 2)

### Parse-check (chat 4)

All modified R files pass `Rscript -e "parse(file=...)"` on
R 4.3.3 / data.table 1.14.10. A minimal data.table merge-collision
reproduction was run to confirm the FIX 15 failure mode.

## Chat 4 continuation additions — FIX 21, FIX 22

### FIX 21 — NN-smoothed sim_mats never produced (CRASH/DESIGN, critical)

**Scope:** `phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`,
`inversion_modules/00_inversion_config.sh`,
`phase_2_discovery/2d_candidate_detection/{run_all.R, 00_config.R,
LAUNCH_run_all.slurm}`

**The bug:** the live `STEP_C01a_precompute.R` (phase 2c) did NOT
produce `<chr>.sim_mat_nn<k>.rds` files. The NN-smoothing loop
existed only in the archived
`_archive/MODULE_5A2_Discovery_Core/snakes/STEP_C01a_snake1_precompute.R`
and was silently dropped when the live precomp was adopted.

**Downstream consequences, all silent:**
- `load_sim_mat()` in `run_all.R` looks for
  `<chr>.sim_mat_nn20/40/80.rds` and only found `nn=0` (via the
  precomp-RDS fallback).
- Phase 4 (D02 NN persistence, requires ≥2 scales) was SKIPPED on
  every run.
- Phase 5 (D09 NN sweep tree, requires ≥3 scales) was SKIPPED on
  every run.
- C01d's `iv$survives_nn40` and `iv$survives_nn80` were always NA
  → D3 score `0.30*nn40 + 0.30*nn80 + 0.40*nn_birth_score` = 0.
- C01d's `iv$nn_birth` was always NA → FIX 16 (tree→scoring merge)
  was a no-op because `tree` was always NULL upstream.
- The entire NN-persistence track — Layer-A's strongest structural
  evidence path — was dead.

**The fix (two parts):**

Part 1 — port the NN-smoothing loop from the archived snake1 into
the live precomp. Inserted right after the precomp RDS save in
`STEP_C01a_precompute.R`. Produces `<chr>.sim_mat_nn<k>.rds` under
`<precomp_dir>/sim_mats/` for each k in `NN_SIM_SCALES`.

Part 2 — calibrate the scale ladder to D09's classifier semantics.
D09 is a persistence-barcode: nn_birth = coarsest scale at which a
block first appears as its own unit. Roots are at the coarsest
scale; splits/disappears/novels track how blocks evolve under
increasing MDS-space smoothing. nn_birth ≥ 200 = INVERSION class.

New ladder: `20, 40, 80, 120, 160, 200, 240, 320` (was `20, 40, 80`)
— catches small (nn_birth=20-80) through multi-Mb (nn_birth=160-320)
inversions. Top scale 320 picks up ~16Mb of structure before
oversmoothing dissolves everything into chromosome-scale ambient.
Going beyond 320 not scientifically useful for this cohort
(largest suspected inversion ~10Mb = NN~200).

**Correct understanding of NN semantics (corrected twice during
this session):** NN is NOT sample-space merging (first wrong guess),
NOT genomic-adjacent smoothing (second wrong guess). It's
**MDS-space k-nearest-neighbor averaging**: for each window, find
its k most-similar (smallest dmat row) windows anywhere on the
chromosome, average their MDS coordinates, rebuild sim_mat. Windows
inside the same structural block share MDS neighbors; larger k
= smaller blocks dissolve into ambient structure. The tree tracks
when a block dies out under increasing k, giving a topological-
data-analysis persistence signature per block.

**What to run:** precomp MUST be rebuilt before any further 2d
detection runs. Updated launcher config is set. ~2× precomp
runtime vs the old 3-scale version (8 scales instead of 3, same
per-scale cost). Expected runtime for catfish cohort on LANTA:
a few minutes per chromosome per scale.

### FIX 22 — 2e GHSL v5 now wired into C01d scoring (DESIGN)

**Scope:** `phase_2_discovery/2d_candidate_detection/run_all.R`
(Phase 8 merge),
`phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
(D10 consumer), `phase_2_discovery/2e_ghsl/README.md` (updated).

**Before FIX 22:** 2e/C04 GHSL v5 was only consumed by two
phase_5_followup supplementary-figure scripts. C01d's D10 read
`partition_stability` from D05's sim_mat column-correlation
clustering. The Clair3-phased haplotype-divergence signal — the
most direct biological evidence for a karyotype difference — had
no path into the main scoring.

**Wiring:**

Part 1 (run_all.R Phase 8): reads
`<ghsl-dir>/annot/<chr>.ghsl_v5.annot.rds` if present. For each
block, aggregates C04's per-window annotations across
`global_window_id ∈ [block.start, block.end]`. Emits six new
columns onto the scoring table: `ghsl_v5_score_max`,
`ghsl_rank_stability_max`, `ghsl_div_contrast_z_max`,
`ghsl_div_bimodal_frac`, `ghsl_pass_frac`, `ghsl_n_scored_windows`.
Gracefully no-ops when `--ghsl-dir` is empty or the annot shard is
missing (Clair3 still running on that chromosome).

Part 2 (C01d D10): blends Clair3-phased and sim_mat signals. When
C04 has ≥3 scored windows in the block:
`d10 = 0.60 * d10_ghsl + 0.40 * d10_simmat`, where
`d10_ghsl = 0.50 * ghsl_v5_score_max + 0.25 * ghsl_div_bimodal_frac
+ 0.25 * ghsl_pass_frac`. When only sim_mat: `d10 = d10_simmat`
(preserves pre-FIX-22 behaviour). Emits `d10_source` column with
values `simmat_only` / `ghsl_and_simmat` for transparency.

**Status of the 4-layer framework after FIX 22:**

- Layer A (dosage/MDS): live via C01a precomp → D01-D08
- Layer B (SV calls): live via C00 sv_prior → C01g boundary catalog
- **Layer C (GHSL): now live via C04 → run_all Phase 8 → C01d D10**
- Layer D (genotype-breakpoint association): via hypothesis tests
  in C01f (phase 4c)

The 2c_precomp README's "4-layer independence framework" claim is
now mostly accurate, not aspirational.

### Files modified (chat 4 continuation)

| File | Changes |
|---|---|
| `phase_2_discovery/2c_precomp/STEP_C01a_precompute.R` | FIX 21 — NN smoothing loop ported from archive |
| `inversion_modules/00_inversion_config.sh` | FIX 21 — `NN_SIM_SCALES` expanded to `20,40,80,120,160,200,240,320` |
| `phase_2_discovery/2d_candidate_detection/run_all.R` | FIX 21 defaults + FIX 22 GHSL merge in Phase 8 |
| `phase_2_discovery/2d_candidate_detection/00_config.R` | FIX 21 CFG$NN_SCALES_DEFAULT updated |
| `phase_2_discovery/2d_candidate_detection/LAUNCH_run_all.slurm` | FIX 21 `--nn-scales` expanded |
| `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | FIX 22 — D10 dual-source scoring + d10_source column |
| `phase_2_discovery/2e_ghsl/README.md` | FINDING 19 superseded by FIX 22 |

### Still NOT fixed (documented for future)

- **NN ceiling vs chromosome size** — NN auto-clamps to
  `min(k, n_windows-1)` so small chromosomes just get fewer
  effective scales. No action needed unless a chromosome has <50
  windows (none expected in this genome).
- **C04 doesn't consume the C01i group registry.** It runs on
  phased Clair3 VCFs directly. Per-karyotype-class validation
  would require registry input. Deferred.
- **D14 landscape classifier still orphaned** (FINDING 18).
  Wiring it is a future project; C01d's shape_class fallback is
  sensible.

### Action required before running

1. **REBUILD PRECOMP** for all 28 chromosomes before any new 2d
   detection runs. FIX 21 changes the precomp output schema
   (adds `<precomp_dir>/sim_mats/<chr>.sim_mat_nn<k>.rds` files).
   Existing precomp RDS files from before FIX 21 will still load,
   but 2d's Phase 4/5 will remain broken until the nn-scaled
   sim_mats are produced.
2. After precomp rebuilds, 2d/run_all.R can be run per-chromosome
   with the new 8-scale ladder.
3. C04 (2e GHSL) can be run per-chromosome as Clair3 phasing
   completes. Phase 8 of run_all.R picks up whatever annot.rds
   shards exist at detection-run time.

## Chat 5 additions (Part B — phase 3 audit + registry wiring)

### Structural change

Flattened `inversion_modules/phase_3_refine/MODULE_5A2_breakpoint_validation/`
one level up. The wrapper directory was unnecessary nesting — phase 3 *is*
the breakpoint validation module. Eight documentation files had their
live path references updated; seven retrospective mentions remain as
intentional "was X, flattened to Y" markers.

### Bug fixes (FIX 23–28)

Six bugs found and fixed in the phase 3 scripts during the audit:

- **FIX 23 (DESIGN)** — `breakpoint_validator_standalone.py`: INV-vs-HET
  Fisher test had a dead `else:` branch (tupling a 2-tuple return into
  a 3-tuple's first slot made `isinstance(...)` always true), and the
  computed values were never consumed. Simplified and added TODO.
- **FIX 24 (CRASH)** — `breakpoint_validator_standalone.py` L648:
  unguarded `y/t*100` division in text report crashed on empty groups.
- **FIX 25 (CRASH)** — `breakpoint_validator_standalone.py` SLURM
  wrapper template referenced `breakpoint_validator.py` (missing
  `_standalone` suffix) → fails on cluster.
- **FIX 26 (CRASH, critical)** — `run_breakpoint_validation.sh` passed
  `matched_delly_snake_candidates.tsv` to STEP02/STEP05 but STEP01
  writes `matched_inv_candidates.tsv`. **The pipeline failed at STEP02
  on first real run.** Chat 4 missed this under the "plumbing only"
  framing.
- **FIX 27 (SILENT)** — STEP03 `seed_summary.tsv` had column drift when
  qualifying status varied across candidates (conditional key add).
- **FIX 28 (CRASH)** — STEP05 manuscript-sentence writer crashed when
  DELLY had calls but Manta had none. Split into four branches.

### FIX 29 v2 — phase 3 wired into phase 4 via evidence registry (architectural)

**Before chat 5:** phase 3 produced flat TSVs, but nothing in phase 4
read them. The architecture across chats described "Layer D OR test"
keys (`q7_layer_d_fisher_or/_fisher_p/...`) with documented consumers
in C01f's `compute_group_validation()` and `compute_candidate_status.R`
— but no writer existed. The only missing piece was the glue.

**After chat 5:**

- STEP03 writes `existence_layer_d` blocks per candidate via the Python
  registry API. The existing schema at
  `registries/schemas/structured_block_schemas/existence_layer_d.schema.json`
  auto-materialises the flat `q7_layer_d_*` keys that downstream code
  already reads. The VALIDATED promotion gate in C01f
  (`fisher_p < 0.05 AND fisher_or > 5`) is now reachable.
- STEP06 writes new `existence_layer_b_bnd_rescue` blocks per paired
  BND junction (both orphan rescues and catalog-matched pairs). This
  supplements the primary Layer B evidence from C00, so phase 4a/4e
  Layer B scoring can see BND-rescued inversions that strict INV
  callers missed.
- New schema file:
  `registries/schemas/structured_block_schemas/existence_layer_b_bnd_rescue.schema.json`.
  Extracts 11 keys (`q7b_bnd_rescued`, `q7b_bnd_rescue_source`,
  `q7b_bnd_pair_bp1/bp2/size_bp`, `q7b_bnd_match_type`, PE/SR per
  junction).
- Both writes are gated on `REGISTRIES_ROOT` being set (default
  `${BASE}/inversion_codebase_v8.5/registries`). Script remains runnable
  standalone without a registry. `CANDIDATE_MAP` TSV can optionally
  translate phase 3's `inv_id` strings into phase 4's canonical
  `candidate_id` naming.
- Earlier first-draft FIX 29 that added a "D13" scoring dimension to
  C01d was **reverted** mid-session — that confused Layer A's 12-way
  composite with Layer D's independent registry contract. C01d is back
  to its pre-chat-5 state with only a header comment documenting the
  revert.

**Files changed (FIX 29 v2 + 23–28):**

| File | Changes |
|---|---|
| `phase_3_refine/breakpoint_validator_standalone.py` | FIX 23/24/25 |
| `phase_3_refine/03_statistical_tests_and_seeds.py` | FIX 27 + FIX 29 v2 (Layer D writer) |
| `phase_3_refine/05_delly_manta_concordance.py` | FIX 28 + docstring alignment |
| `phase_3_refine/06_bnd_inversion_signal.py` | FIX 29 v2 (Layer B BND-rescue writer) |
| `phase_3_refine/02_extract_breakpoint_evidence.py` | docstring alignment (FIX 26 part) |
| `phase_3_refine/run_breakpoint_validation.sh` | FIX 26 + registry argument wiring |
| `phase_3_refine/00_breakpoint_validation_config.sh` | REGISTRIES_ROOT + CANDIDATE_MAP defaults |
| `phase_3_refine/README.md` | Full 4-layer architecture rewrite |
| `registries/schemas/structured_block_schemas/existence_layer_b_bnd_rescue.schema.json` | NEW |
| `registries/schemas/structured_block_schemas/INDEX_remaining_blocks.json` | New block indexed |
| `phase_4_postprocessing/tests/test_registry_sanity.py` | Layer D + BND-rescue fixtures |
| `phase_4_postprocessing/README.md` | Corrected phase-3 contract note |
| `phase_4_postprocessing/4a_existence_layers/README.md` | Updated Layer B/D source table |
| `phase_2_discovery/2c_precomp/README.md` | Corrected "two scripts" table for phase-3 outputs |
| `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | REVERTED first-draft D13 wiring back to 12-dim composite; header comment documenting the revert |

### Gate question — resolved

Phase 3 was never dead — it just wasn't wired through the registry.
The chat-4 handoff's observation that "phase 3 output filenames don't
appear anywhere in phase 4" was correct **at the flat-TSV level** and
misleading **at the architectural level**. Phase 3 contributes via the
registry-blocks indirection, which is the correct decoupling for a
pipeline this size. Chat 5 installed the writers; the readers have
been in place for weeks.

### Parse-check

All modified Python files pass `python3 -c "import ast; ast.parse(...)"`.
Modified shell files pass `bash -n`. JSON files pass `json.load`. The
modified R file (C01d) was reverted to its pre-chat-5 state plus an
additive header comment; recommend an on-cluster `Rscript -e "parse(...)"`
as final sanity before running.

### Still not fixed (deferred)

- **C01i seeded-init wiring.** The seed files
  `phase_3_refine/04_deconvolution_seeds/{inv_id}_seeds.tsv` are
  designed as k-means init input for `STEP_C01i_decompose.R`. Not yet
  wired; flagged as future work.
- **`annotate_population_confidence.sh` relocation.** Lives in
  `phase_3_refine/` but logically belongs at `MODULE_4H_ALL_Manta/`
  (sources `00_manta_config.sh`, consumed only by STEP_C00 in 2c).
  Directory-hygiene cleanup; contract with STEP_C00 verified clean.
- **`existence_layer_b` schema description cleanup.** Still says
  "C00 build_flashlight" (old name) — should say "build_sv_prior".
  Schema is 2c-owned, outside Part B scope.

## Chat 6 additions (Part C — phase 4 consumer-side audit)

### Framing

Mid-session clarification: the registry library itself is not yet
implemented. The specs and plan were set earlier; the decision was to
complete the pipeline code and review it first, then code the lib and
wire it up. So throughout the current codebase, readers forward-declare
the keys they will consume and writers declare the keys they will emit.
Reader/writer parity and schema-derived `*_detected`/`*_tested` flags
are lib-era concerns, not current-state bugs.

This reframes several "orphan key" observations as intentional
forward-declarations. Only findings that are real bugs independent of
lib status are tracked as fixes.

### Audit findings (chat 6)

Real bugs (three, all flagged not applied — audit-only session):

- **Finding 1 (CRITICAL, proposed FIX 31)** — C01f L2322–2323 hard-codes
  `layer_d_fisher_p = NA, layer_d_fisher_or = NA` in the call to
  `compute_group_validation()`. The VALIDATED-promotion gate at L525–526
  requires `is.finite()` on both, so it is unreachable. Chat 5's patch
  template (`patches/01_C01f_comp_from_registry.R` L282–283) has the
  same NA hard-code with a comment saying STEP03 writes Layer D
  separately — but never supplied the reader. The fix pattern already
  exists two lines above in the same file (used for
  `q6_group_validation` and `q6_validation_promotion_cap`): wrap
  `reg$evidence$get_evidence(cid, key)` in a tryCatch with NA fallback.
- **Finding 2 (HIGH, proposed FIX 32)** — 4e `compute_candidate_status.R`
  Pathway A gates (L298, L309) and `group_validation_gate.R` L50–52 all
  gate on `fisher_p < 0.05` without the `fisher_or > 5` clause. Schema
  description and C01f's `compute_group_validation()` require both.
  After Finding 1 is fixed, 4e and C01f disagree on significant-but-low-OR
  candidates. One-line fix per gate.
- **Finding 5 (CRITICAL, proposed FIX 30)** — `phase_3_refine/03_statistical_tests_and_seeds.py`
  L373/L375/L377 reference `r['total_score']` when writing seed files
  for qualifying candidates, but the evidence TSV schema has no such
  column. `KeyError` crash fires for every candidate that qualifies
  for seeding — the highest-quality candidates, which are the ones
  seeds are most useful for. Chat 5's FIX 27 fixed the adjacent header
  drift in `seed_summary.tsv` but didn't exercise the qualifying
  branch past the header. Fix requires a column choice from Quentin;
  best candidate is `int(r.get('pe',0)) + int(r.get('sr',0))`.

Lib-design notes (two, not bugs):

- **Finding 3** — 4e readers declare `q7_layer_{a,b,c}_detected` and
  `q7_layer_d_tested` as gate flags, and declare Layer A/B/C detail
  keys (e.g. `q7_layer_b_delly`, `q7_layer_c_ghsl_quality`). No writer
  exists, and the current schemas extract different key families.
  **Intentional pending lib.** When the lib lands, the natural design
  is schema-level `keys_extracted` entries of the form
  `{ "key": "q7_layer_d_tested", "from_derived": "block_exists" }` and
  the registry loader's `_extract_keys_from_schema` gets a
  `from_derived` branch. Until then, running 4e end-to-end puts every
  candidate into `D_insufficient`/`D_artifact` — expected state.
- **Finding 4** — 11 `q7b_bnd_*` keys written by STEP06 have no
  consumer in phase 4. **Intentional pending lib.** Natural read-side
  wiring: Pathway E's `layer_b` flag becomes TRUE when either a
  catalog INV exists or `q7b_bnd_rescued == TRUE`; add an "E_bnd_rescued"
  variant. Concrete design deferred to the lib chat.

### Empirical verification (chat 6)

Smoke test run end-to-end:

1. Built scratch registry at
   `/home/claude/chat6_work/scratch_registry/` pointing at the real
   `api/` and `schemas/` trees.
2. Built a fake strong-signal evidence TSV (60 samples, 10% REF /
   40% HET / 90% INV support fractions).
3. Ran `03_statistical_tests_and_seeds.py --registries_root ...` with
   thresholds tuned to avoid Finding 5's crash branch.
4. Confirmed: Fisher OR=27.0, p=1.66e-06, Armitage p=3.96e-07;
   `existence_layer_d` block written with `validation_status=validated`;
   8 keys extracted into `keys.tsv` (fisher_or, fisher_p, two CIs,
   armitage_z/p, n_inv_with_support, n_inv_total) — **no**
   `q7_layer_d_tested` key, as expected given current lib status.

Silent standalone fallback verified (wrong `--registries_root` →
prints skip message → continues with flat-TSV output).

### Files added (chat 6)

| File | Purpose |
|---|---|
| `AUDIT_LOG_chat6_2026-04-17.md` | Full audit narrative, all five findings, smoke test details, reframing after lib-status clarification |

No code files modified. This was an audit session by design; all
proposed fixes flagged for follow-up.

### Parse-check (chat 6)

No code changes. Only the audit log + this session summary
append were produced.

### For the next chat

Short fix-application session for FIX 30 + 31 + 32. All three
surgical, all three need at most three tryCatch-style edits total,
all three should parse-check in-container (Python + schema) or on
LANTA (R).

Findings 3 and 4 stay deferred for the registry-lib implementation
chat. Use chat-6 audit log as the consumer-side specification
document for that work: the 4e readers already tell the lib which
flat keys the pathway engine expects.

## Chat 7 additions (Part D — apply chat-6 fixes + begin 4a audit)

### Fixes applied (three from chat 6, one new from 4a audit)

- **FIX 30 (CRASH applied)** — STEP03 seed writer: `r['total_score']`
  → `int(r.get('pe',0) or 0) + int(r.get('sr',0) or 0)` with defensive
  cast. Empirical smoke test confirms: strong-signal candidate
  qualifies with default thresholds, no crash, seed file written
  (REF=0, HET=2, INV=5 evidence scores). Python AST OK.
- **FIX 31 (CRITICAL applied)** — C01f L2316–2325: replaced two `NA`
  hard-codes with `reg$evidence$get_evidence()` tryCatch reads
  mirroring the pattern used two lines above for `q6_group_validation`.
  Patch template at `patches/01_C01f_comp_from_registry.R` also
  updated so re-applying the patch doesn't re-introduce the bug.
- **FIX 32 (HIGH applied)** — 4e `compute_candidate_status.R` L272 +
  L298 + L309, and 4c `group_validation_gate.R` L50–59 + L89 + L96:
  added `fisher_or` reader and `fisher_or > 5` clause to the Pathway A
  gates and to `or_passed`. Reason strings updated. Now consistent
  with schema description and C01f after FIX 31.
- **FIX 33 (CRITICAL applied — newly found)** — C01d L826–840: added
  `tier` recomputation after the `dim_positive` recompute at L827–832.
  Prior behaviour: `tier` was assigned inside the per-candidate loop
  with `d11 = d12 = 0` defaults, then `dim_positive` got recomputed
  with the real D11/D12 but `tier` stayed frozen. Net effect: a
  candidate whose dim_positive was 8 after recompute could have
  `tier = 2` assigned during the loop, producing an internally
  inconsistent `(dim_positive = 8, tier = 2)` row in the catalog.
  Fix recomputes tier with same ≥8/≥6/≥4 ladder and re-applies the
  peel-disappeared downgrade guarded by column existence.
- **FIX 34 (COSMETIC applied)** — C01d L32–51 header + L412 inline:
  updated stale "10 scoring dimensions" + "Tier 1: ≥7/10" to "12
  dimensions" and "≥8/≥6/≥4" matching the code.

### Phase 4a audit — C01d pass complete

Verified C01d's 12-dimension scoring loop (D1–D12), weights sum to
exactly 1.00 (0.14+0.08+0.13+0.06+0.05+0.06+0.06+0.12+0.06+0.06+0.09+0.09
= 1.00 — checked). Confirmed D11 boundary-concordance reader at
L555–618 correctly matches C01g's column schema. Confirmed D12
seeded-region reader at L467–546 has legacy fallback for pre-rename
`snake1_core_regions_*.tsv.gz`. Confirmed chat-5 revert header is in
place (L99–113).

### Phase 4a audit — findings not fixed (pending lib work)

- **Finding 8 (NOTE)** — No 4a script writes any registry block.
  C01d: zero writes. C01e: zero writes (figures only). C01g: attempts
  writes via helpers in `utils/registry_key_helpers.R` but that file
  does not exist in the working tree. Silent skip (tryCatch +
  `exists(..., mode="function")` guards). Under the lib-not-yet-
  implemented framing this is expected; flagged for the lib chat.
- **Finding 9 (NOTE)** — C01g's silent skip when the helpers file is
  missing could cause a real-dataset run to silently produce no
  registry writes. Suggest adding a WARNING message when
  `.bridge_available` is TRUE but the helper file is missing. Flagged
  for the lib chat.

### Phase 4a audit — NOT yet done (deferred to chat 8)

- **C01g (1428 lines).** Biggest script in 4a; not yet audited.
  Needs: verify five source inputs are all read; Cheat 17 archival
  per README; removed flags are gone; output columns match C01d's
  reader expectations.
- **C01e (694 lines).** Low priority — figures only.
- **End-to-end contract check.** Verify C01g output columns match
  the names C01d's D11 reader expects: `boundary_verdict`,
  `boundary_bp`, `matched_core_id`, `side`, `n_cheats_supporting`,
  `chrom`.

### Files edited (chat 7)

| File | Changes |
|---|---|
| `phase_3_refine/03_statistical_tests_and_seeds.py` | FIX 30 |
| `phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R` | FIX 31 |
| `phase_4_postprocessing/4c_group_validation/group_validation_gate.R` | FIX 32 |
| `phase_4_postprocessing/4e_final_classification/compute_candidate_status.R` | FIX 32 |
| `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | FIX 33 + FIX 34 |
| `phase_4_postprocessing/patches/01_C01f_comp_from_registry.R` | template update for FIX 31 |

Plus `AUDIT_LOG_chat7_2026-04-17.md` added.

### Parse-check (chat 7)

Python (FIX 30) AST-parse OK. R files (FIX 31/32/33/34) have matched
brace and paren counts in all five edited files. No R parse checker
in-container; LANTA `Rscript -e "parse(file=...)"` required before
production run.

### For the next chat

Chat 8 = finish 4a audit. Start with C01g (largest, most inputs).
Then scan C01e. Then end-to-end contract check. Close with an
explicit verdict on whether 4a can produce a catalog phase 4b/c/d/e
can consume.

## Chat 7 completion — C01g + C01e audited, 4a verdict delivered

After the initial fix-pass writeup, Quentin clarified that audit-only
is not a rule — fix-as-we-go is the posture unless a finding needs a
design call. Completed C01g + C01e audit in the same chat.

### FIX 35 applied — SILENT — C01g dedup missing matched_core_id

C01g's dedup aggregated `candidate_ids` as semicolon-joined but never
produced a canonical single-value `matched_core_id`. The
registry-write block at L1407 guards on this column's existence, so
the per-boundary registry write loop was structurally unreachable —
separately from Finding 9's missing helper file. Fix derives
`matched_core_id` from the dedup's best_idx row with two fallback
levels.

### C01g audit — no other fixes needed

All five source readers wired correctly; dead flags properly
accepted-ignored; Cheat 17 archival clean; column contract with C01d's
D11 reader verified end-to-end.

### C01e audit — no fixes needed

Dead flags handled; chat-3's triangle-composition bug fix is in
place; reads only `interval_id`, `tier`, `final_score` from C01d — all
present after FIX 33/34.

### Phase 4a verdict

- **Flat-TSV path (candidate_scores + boundary_catalog_unified):**
  ready today. Phase 4b/c/d/e can consume C01d's catalog; C01d's D11
  reader consumes C01g's boundary catalog correctly.
- **Registry-block path (existence_layer_a/b/c blocks):** not
  produced by 4a today. Three items need to land together in the
  lib chat: (1) `utils/registry_key_helpers.R` with
  `register_C01g_boundary` + `store_C01g_boundary` (C01g is ready
  for this after FIX 35); (2) Layer A block writer in C01d; (3)
  registry loader extension for `{from_derived: "block_exists"}` so
  the `q7_layer_*_{detected,tested}` orphan family from chat 6
  Finding 3 gets populated.

### Final chat-7 file tally

Code edits: 6 source files + 1 patch template. All surgical, matched
brace/paren counts preserved. Needs LANTA `Rscript -e "parse(...)"`
for the five R files (container has no R).

---

## Chat 8 — phase 4b audit (Part E)

Scope: `4b_group_proposal/` — the four scripts (decompose,
multi_recomb, nested_composition, seal) that propose genotype groups
and register them in `sample_registry`. Handoff posture: fix-as-we-go
with a three-fix ceiling; flag lib-pending items.

### Files audited

- `README.md` — four-script architecture confirmed
- `STEP_C01i_decompose.R` (4b.1, 451 lines) — k-means on mean PC1
- `STEP_C01i_b_multi_recomb.R` (4b.2, 417 lines) — three-signal
  recombinant detector
- `STEP_C01i_d_seal.R` (4b.4, 417 lines) — synthesis + group
  registration
- `STEP_C01i_c_nested_composition.py` (4b.3, 446 lines) — Engine B
  composition analysis wrapper
- `nested_composition_core.py` (271 lines) — vendored classify_structure
  (spot-checked)
- `lib_decompose_helpers.R` (327 lines) — shared helpers, registry
  wrappers
- `engine_b_smoke_test.R` — skipped (diagnostic, not correctness-path)

Total ≈ 2280 lines read.

### Three fixes applied (ceiling reached)

**FIX 36 — COSMETIC — misleading decision-rule comment**

`STEP_C01i_b_multi_recomb.R` L316: inline comment described the
decision rule as `(S1 AND S2) OR S1 OR S3` (S1-alone is always
sufficient). Code actually implements `(S1 AND S2) OR S3 OR (S1 AND
mosaic>100kb)`, which matches the README. Comment rewritten to
describe the three-term rule and the noise-tolerance rationale. No
runtime change.

**FIX 37 — CRASH (dormant today) — undefined `het_carriers` in Signal 3**

`STEP_C01i_b_multi_recomb.R` L233 (function `flashlight_hemi_signal`):

```r
sample_dels <- dels[sample_id %in% unlist(het_carriers)]  # BUG
```

`het_carriers` is not defined anywhere in scope. This is a copy-paste
from the cheat2 pattern in decompose (L290–294) that never got
adapted for the per-sample Signal 3 context. Dormant today because
flashlight isn't universally available on LANTA; becomes a hard crash
("object 'het_carriers' not found") the first time `get_internal_dels()`
returns rows with flashlight loaded. Fix: per-row mask using the
scalar `sample_id` against each row's `het_carriers` list-column,
tryCatch'd to all-FALSE on failure.

**FIX 38 — SILENT — `%||%` precedence inflates mosaic lengths by 100 kb**

`lib_decompose_helpers.R` L88–89 (`extract_pc_loadings`):

```r
win_starts = wins$start_bp %||% wins$mid_bp - 50000L,
win_ends   = wins$end_bp   %||% wins$mid_bp + 50000L,
```

`%||%` (like all `%any%` operators in R) binds tighter than `+`/`-`,
so these parse as `(wins$start_bp %||% wins$mid_bp) - 50000L` —
subtracting 50 kb from every real window start when `start_bp` is
present. Every window interval inflated by 100 kb. Downstream in
multi_recomb, every `mosaic_length_bp` inflated by 100 kb.

**On the 50000 itself (Quentin's pushback was right):** the precomp
(C01a `STEP_C01a_precompute.R` L849) ALWAYS writes real `start_bp`
and `end_bp` columns, so the `%||%` fallback branch was never
supposed to fire. The 50000 was arbitrary — C01a's actual window
sizes come from the multiscale ladder (20/40/80/120/160/200/240/320),
not a fixed ±50 kb. **Fix v2:** strip the dead fallback entirely, use
the real columns directly.

**Downstream impact against real cheat24 thresholds** (50 kb short /
200 kb long, from `flashlight_v2/cheats/cheat24_recombinant_prior.R`):

1. Signal-1-alone threshold (multi_recomb L320): `mosaic_length_bp >
   100000`. A true-50kb mosaic reads as 150 kb, spuriously crossing
   the threshold → samples pulled into RECOMBINANT that should have
   stayed in HOM_REF/HET/HOM_INV.

2. cheat24 event classification: a real 10–40 kb gene conversion
   reads as 110–140 kb (no longer short) → misclassified as
   "ambiguous" instead of "gene_conversion". Most true GC events lost.

Both group composition and GC/DCO ratios biased.

### Three findings flagged (not fixed per handoff / budget)

- **Finding S — STEP03 seeds still unwired.** `04_deconvolution_seeds/
  {inv_id}_seeds.tsv` are written by STEP03 (post chat-7 FIX 30) but
  no 4b script reads them. Flashlight-seeded k-means works today when
  flashlight is available; STEP03 seeds would be a useful middle tier
  but wiring requires a design call (priority order when both
  available? disagreement resolution?). Chat 5 deferred this; chat 8
  continues to defer.
- **Finding T — README structure_type drift.** README lists 4
  structure types but the code produces 7. Doc-only drift; code
  contract between 4b.3 and 4b.4 intact. Also: decompose's header
  comment L16 promises a "per-window dosage track" output but the
  code actually produces a per-window class track (PC1 k-means
  labels), not dosage. Batch with Finding T doc pass.
- **Finding U — forward-declared flat keys** `q1_composite_flag`,
  `q2_decomp_quality_flags` have no live readers. Same pattern as
  chat 7 Finding 8; consistent with lib-not-yet-implemented framing.
- **Finding V — cheat24 threshold divergence (surfaced while tracing
  FIX 38).** Real cheat24 uses 50kb/200kb for short/long mosaic;
  inline fallback in multi_recomb uses 100kb/500kb. Same candidate
  gets different event class depending on flashlight availability.
  Needs Quentin to pick the authoritative thresholds before aligning.

### Contract verification summary (audit outcomes)

- **Group-registration contract (seal → comp_from_registry):** gids
  `inv_{cid}_{HOM_REF|HET|HOM_INV|RECOMBINANT}` match exactly. HOM_STD
  alias registered for v9 back-compat. ✅
- **Quality metric contract (decompose → seal):** field names
  (`silhouette_score`, `bic_gap_k3_vs_k2`, `phase_concordance`) match
  between writer and reader. ✅
- **Recombinant-detection contract (multi_recomb → seal):** verified
  post-FIX-37 and FIX-38. Signal combination rule correct (post-
  FIX-36 comment); mosaic lengths accurate (post-FIX-38); Signal 3
  no longer crashes (post-FIX-37). ✅
- **Tier flow:** all three 4b scripts filter `tier <= 3`. Post-FIX-33
  tier values accurate. Set of candidates entering 4b unchanged; only
  the internal tier distribution shifts. ✅
- **Promotion-cap contract (seal → C01f):** `q6_validation_promotion_cap`
  written by seal L368, read by C01f L2311 (chat-7 FIX 31). ✅

### Phase 4b verdict

Ready to feed phase 4c. FIX 37 and FIX 38 are correctness-meaningful
(crash-on-activation and silent-bias respectively); FIX 36 is purely
cosmetic. With the three fixes applied, phase 4b:
- Produces correct group compositions (post-FIX-38).
- Produces correct event-class distributions (post-FIX-38).
- Does not crash when flashlight is fully loaded (post-FIX-37).
- Has honest in-code documentation of its decision logic (post-FIX-36).

Deferred items (Findings S, T, U) are either design-gated (S) or
documentation/forward-declaration (T, U) and don't block the
flat-group-registration path that C01f consumes today.

### Parse-check backlog

Chat 8 adds two files to the existing chat 7 backlog:
- `4b_group_proposal/STEP_C01i_b_multi_recomb.R` (FIX 36 + 37)
- `4b_group_proposal/lib_decompose_helpers.R` (FIX 38)

Total after chat 8: nine R scripts + one Python script need LANTA
`Rscript -e "parse(file=...)"` + smoke tests before the next
pipeline run. Container has no R.

### For chat 9

**Recommended scope:** phase 4c full audit. `STEP_C01f_hypothesis_tests.R`
at 2512 lines is the biggest R script in the tree and the central
consumer of everything phase 4b produces. Worth a whole chat.

Alternative: combined 4d + 4e audit, which is fewer lines but more
cross-script contracts to verify.

After 4c (or 4d/4e) is done, the lib-design chat can start from a
fully-audited baseline.

### Final chat-8 file tally

Code edits: 2 source files (both in 4b_group_proposal/). Surgical,
brace/paren/bracket balance verified (0/0/0 across all four main 4b
files). Needs LANTA parse-check.

Documentation: this SESSION_SUMMARY append, AUDIT_LOG_chat8, FIXES_APPLIED
chat-8 section.

---

## Chat 9 addendum (2026-04-17) — phase 4 end-to-end audit + 4e rework

Scope: started as Option A (phase 4c audit, three-fix ceiling).
Quentin asked for full workflow understanding first, then authorized
scope expansion to "code anything missing" for phase 4e. Final scope:
phase 4 end-to-end connection audit + substantial 4e rework.

**Six fixes applied (39–44):**
- FIX 39: T9 jackknife vocabulary — `few_family_contributing` now classifies correctly in C01f's `compute_group_validation` (was silently falling through)
- FIX 40: composite_flag key consolidation — schema key renamed to `q1_composite_flag`, seal's redundant write removed (Option A+B)
- FIX 41: cheat30 `compute_pairwise_ibs` alias scope bug — was inside function body, now module-scope; would crash on every run_cheat30 call
- FIX 42: full rewrite of `build_key_spec` in 4e/compute_candidate_status to match v2 spec (367 keys) + per-Q aspirational-key marking + `compute_completion` rework (now has `pct` and `pct_of_spec` metrics)
- FIX 43: jackknife_status readers prefer `q6_family_linkage` (actually written) over `q7_t9_jackknife_status` (aspirational, no writer). gate.R + characterize_q7 + `has_groups` check extended to accept HOM_REF/HOM_STD aliases
- FIX 44: `characterize_candidate.R` now sources `group_validation_gate.R` via four-candidate path lookup; previously crashed on `assess_group_validation` call

**Seven files modified.** Plus this summary, audit log, and FIXES_APPLIED append.

**Not finished this chat:**
- `run_characterize.R` driver (the main pending 4e deliverable)
- `SPEC_VS_REALITY.md` writer-audit doc
- Orchestrator break (`run_phase4.sh` + 4d launchers)
- Phase 4a audit
- Phase 2/3 writer wiring (the aspirational-keys TODO)

Chat 10 handoff prompt: `HANDOFF_PROMPT_chat10_2026-04-17.md`

---

## Chat 10 addendum (2026-04-17) — Phase 4e finish line

Scope: close out Phase 4e per chat-9 handoff. Priorities 1, 2, and 3 all
delivered. Six fixes applied (FIX 45–50). Two of them (FIX 46, FIX 47)
were surfaced by smoke-testing and would have prevented every
`characterize_candidate()` call from ever succeeding — a parse-clean but
runtime-broken state chat 9 couldn't have detected without invocation.

**Six fixes applied (45–50):**
- FIX 45 (FEATURE): `run_characterize.R` — the 4e driver, 420 lines,
  sibling discovery + three-tier loader + crash-proof main loop + four
  outputs (characterization.tsv, per-candidate .txt, summary, merged
  full-report). Closes Finding AF.
- FIX 46 (CRASH): `characterize_candidate.R` %||% ordering bug — chat 9's
  FIX 44 gate-locator block used `%||%` before definition. Every
  `source()` of the file would crash with "could not find function %||%".
  Reordered helpers before gate-locator; also wrapped `sys.frame(1)$ofile`
  in tryCatch.
- FIX 47 (CRASH): `safe_num` zero-length failure — `is.finite(numeric(0))`
  returns `logical(0)`, and `if (logical(0))` throws. Every missing
  optional key triggered it. Fixed `safe_num`, `safe_bool`, `has_key` to
  guard length + NULL before the checks.
- FIX 48 (FEATURE): `run_phase4b.sh` emits `PHASE4B_SEAL_JID=<id>` trailer
  for parent orchestrators to parse deterministically.
- FIX 49 (ARCH): `run_phase4.sh` calls `run_phase4b.sh` sub-DAG instead
  of obsolete `LAUNCH_C01i_decomposition.sh` single-job launcher.
- FIX 50 (FEATURE): production `LAUNCH_group_cheats.sh` — corrected paths
  (`cheats/` → `4d_group_dependent/`), added cheat27/28/29 calls, left
  aspirational (Fst sub-block, burden regression) as commented-out gates.
  Closes Finding Y.

**SPEC_VS_REALITY.md delivered** (240 lines). 271/367 keys wired = 73.8%.
Breakdown per Q, aspirational keys grouped by 12 producer modules, wiring
roadmap with priority order for chat 11+. Closes Finding AG.

**Seven files modified; three new.** `run_characterize.R`,
`SPEC_VS_REALITY.md`, `LAUNCH_group_cheats.sh` are new.
`characterize_candidate.R`, `run_phase4.sh`, `run_phase4b.sh` modified.

**The finish line:** `compute_candidate_status.R` + `run_characterize.R`
now both runtime-work on a synthetic 3-candidate registry (healthy +
artifact + empty cases). Real-data run on LANTA is recommended as the
chat 11 priority-1 smoke test.

**Not finished this chat (preserved for chat 11):**
- Phase 4a audit (C01d / C01e / C01g end-to-end)
- Real-data LANTA smoke run (compute_candidate_status + run_characterize)
- Q7B breakpoint_evidence_audit flat-key extraction (25 keys, biggest
  aspirational block to close)
- Finding W verification (register_C01f_keys external helper)
- Phase 2/3 writer wiring for remaining 96 aspirational keys

Chat 10 audit log: `AUDIT_LOG_chat10_2026-04-17.md`

## Chat 11 — registry API extensions (buildout)

**22 new methods across R + Python bindings.** Sample and interval
registry APIs extended per chat-11 handoff. Evidence registry untouched
(mature; Phase 3 reference writer).

**Scope:** architectural library work. No upstream writer wiring (chats
12–17). No three-fix ceiling.

### Methods added

- **Interval registry (9):** `get_children`, `get_parent`, `get_ancestors`,
  `get_descendants`, `get_overlapping`, `get_nested_within`,
  `classify_relationship`, `update_candidate`, `bulk_add_candidates`.
  All support nested inversions, overlapping inversions, and the
  chat-13 seeded-region bulk writer.
- **Sample registry (10):** `get_sample_metadata`, `get_sample_groups`,
  `get_family`, `list_families`, `get_family_members`,
  `list_carriers_across_candidates`, `compute_group_overlap`,
  `list_recombinant_candidates`, `get_groups_for_candidate`,
  `find_co_segregating_groups`. Support family/pedigree queries,
  reverse lookup, cross-inversion sample tracking, and co-segregating
  repeat-call detection.
- **Query API composites (3, R only):** `nested_family_tree`,
  `overlap_clusters`, `sample_inversion_load`.

### Findings

- **AJ** (fixed): `get_master` was aliased to `full$list_groups` (wrong
  table). Rewired to `full$get_master`.
- **AK** (fixed): `utils/sample_registry.R::add_group()` column-class
  drift bug. `created` column auto-parsed as POSIXct on re-read but
  appended as character, fataling `rbind` on the second call in a
  session. One-line `colClasses` fix. Blocking for chat-13 bulk writers.
- **AL** (deferred to chat 12): `load_registry` name collision between
  `registries/api/R/registry_loader.R` and `utils/sample_registry.R`.
  Requires rename + cross-codebase call-site update.
- **AM, AN** (resolved in tarball): orphan duplicate files in
  `inversion_modules/utils/` (`sample_registry.R_v254`, `sample_map.R`)
  byte-identical to canonical `utils/` copies. Deleted.

### Files

- Modified: `registries/api/R/registry_loader.R` (+600 lines),
  `registries/api/python/registry_loader.py` (+409 lines),
  `utils/sample_registry.R` (+7 lines)
- Created: `registries/tests/test_interval_registry_extensions.R`,
  `registries/tests/test_sample_registry_extensions.R`,
  `registries/API_CHEATSHEET.md`, `AUDIT_LOG_chat11_2026-04-17.md`
- Deleted: `inversion_modules/utils/sample_registry.R_v254`,
  `inversion_modules/utils/sample_map.R`

### Not finished this chat (preserved for chat 12)

- Upstream writer wiring: precomp + SV prior + block detect
  (`STEP_C01a_precompute.R`, `STEP_C00_build_sv_prior.R`,
  `PHASE_01C_block_detect.R`) — these populate interval_registry windows
  and candidates, write `existence_layer_a` skeleton and
  `existence_layer_b` blocks
- Finding AL rename (`utils/sample_registry.R::load_registry` →
  `load_sample_groups_api`)

**Chat 11 end-state:** registry library is ready to receive upstream
writes. Every downstream consumer's query pattern (nested, overlapping,
recombinant, cross-inversion) is now a single-call method. C01f's
`comp_from_registry` function can collapse from 8 calls to 1 via
`get_groups_for_candidate`.

Chat 11 audit log: `AUDIT_LOG_chat11_2026-04-17.md`

## Chat 11.5 — decompose/multi_recomb redesign per chat-9 spec + C01j/k/l/m dispatch

**Context.** Review of chat-9 design document surfaced that the initial
draft of the decompose-v2 upgrade had drifted from the three-tier spec.
User also supplied four existing orphan phase-4 scripts
(STEP_C01j_regime_compatibility_engine.R,
STEP_C01l_local_structure_segments.R,
STEP_C01m_distance_concordance.R,
STEP_C01k_annotated_simmat.R) that were forgotten during the v10 rewrite
but are the proper implementations of compatibility-based regime
detection, per-segment ENA/Δ12 boundary sharpness, multi-scale sample
concordance, and the integrative per-chromosome figure, respectively.

**Decisions.**

1. **Dispatch the four orphan scripts into proper phase-4 homes**
   (verbatim, no logic changes; tuning deferred to chat 12):
   - `4a_existence_layers/STEP_C01j_regime_compatibility_engine.R`
   - `4a_existence_layers/STEP_C01l_local_structure_segments.R`
   - `4a_existence_layers/STEP_C01m_distance_concordance.R`
   - `4e_final_classification/STEP_C01k_annotated_simmat.R`
2. **Reframe chat-11.5 draft work.** The per-sample CUSUM that was
   originally a recombinant detector is NOT a recombinant detector —
   per-SNP CUSUM flags ≤500 bp gene-conversion tracts. Renamed to
   `gene_conversion_detector.R` with narrow scope: windowed binning
   (40 SNPs / 10-SNP step), max-tract-bins gate (default 5 bins).
   GC tracts become auxiliary annotations, never gate recombinant
   classification.
3. **Recombinant detection uses C01j.** The compatibility engine
   computes per-window sample × sample Hamming distance, Ward
   hierarchical clustering, adaptive k; samples are assigned to
   compatibility groups per window. Recombinants show regime changes
   within an interval. This is the Tier-2 per chat-9 design, in its
   correct form (not the drifted CUSUM).
4. **Boundary sharpness uses C01l.** Per-segment ENA/Δ12/entropy in
   5 segments (left_flank / inv_left_half / inv_core / inv_right_half /
   right_flank) replaces the ad-hoc staircase boundary_scan. The
   boundary_scan schema from chat 11 remains as a generic container;
   the contents come from C01l derivation.
5. **Family-vs-inversion signal uses C01m.** Multi-scale pair
   concordance (distances 80 / 160 / 320 / 640 windows) gives a
   sample-pair level signal distinguishing inversion-carrier LD
   (persistent across distance) from family LD (decaying).
6. **Per-chromosome integrative figure uses C01k.** Previously orphan;
   now the phase-4e reporting layer.

**Combination rule (per chat-9 §Combination):**
- R (regime change from C01j) AND G (GHSL SPLIT) → RECOMBINANT, HIGH
- R AND NOT G, GHSL insufficient → RECOMBINANT, MEDIUM (regime_only)
- R AND NOT G, GHSL sufficient → recomb_disputed
- NOT R AND G → recomb_ghsl_only (flag)
- NOT R AND NOT G → NOT_RECOMB (Tier 1 decompose class stands)
- Gene-conversion tracts (C) are orthogonal annotation, never gate.

**Four new Tier-2 evidence block schemas:** regime_segments,
local_structure_segments, distance_concordance, gene_conversion_tracts.
Four new query API methods on the registry: `reg$query$regime_segments`,
`local_structure_segments`, `distance_concordance`,
`gene_conversion_tracts`.

**Ported forward from earlier work:** `lib_step03_seed_loader.R`
(Tier-1 flashlight+STEP03 seed combiner, closes Finding S,
no-unsupervised-fallback), `lib_ghsl_confirmation.R` (Tier-3
confirmation).

**Pending for chat 12 (the next session):**
- Phase-4 end-to-end coherence audit (see `HANDOFF_PROMPT_chat12_*.md`).
- Rewrite `derive_R_from_regime` in `lib_recomb_combination.R` with
  DAG-based semantics: per-sample mini DAGs of regime transitions,
  cohort-level breadth of deviation, soft thresholds instead of rigid
  A→B→A matching. This is the primary code change for chat 12.
- Wire decompose to use `lib_step03_seed_loader.R` and remove
  unsupervised fallback.
- Wire multi_recomb to use `lib_recomb_combination.R`.
- Redirect C01j/k/l/m outputs from standalone TSV to
  `reg$evidence$write_block()` calls.
- Add per-sample DAG diagnostic plot.

Registry loader now 1328 lines (was 1166 after chat 11).

Chat 11.5 audit log: merged into this session summary (no separate
chat-11.5 audit; the work gets audited as part of chat-12's coherence
review).
