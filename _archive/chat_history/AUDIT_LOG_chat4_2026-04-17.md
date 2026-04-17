# AUDIT LOG — chat 4 of the 4a audit series (2026-04-17) — PART A

Scope: phase_2/2d detection core, phase_2/2e GHSL, and the C01f
triangle-fallback port (closing a loose end from chat 3). Part B
(phase 3 unpack + audit) deferred to a subsequent chat.

## Summary

Six findings added this session (FIX 15–17, 20 applied; FINDING 18–19
documented, not fixed). All applied fixes parse-check clean.

## Fixes applied

### FIX 15 — `contrast` column collision (SILENT, manuscript-relevant)

**File:** `inversion_modules/phase_2_discovery/2d_candidate_detection/run_all.R` Phase 8

**Problem:** D01's staircase detector writes a `contrast` column
(`block_height - bg_sim`) into the blocks table. D08's block-scoring
step also writes `contrast` (`inside_mean - flank_outside_mean`). In
run_all.R Phase 8 at line 412, `tab` (copied from candidates → blocks
→ has D01's `contrast`) is merged with `raw_scores[, .(block_id,
contrast, squareness, ...)]`. data.table's default `merge()` on a
shared non-join column produces `contrast.x` (D01) and `contrast.y`
(D08) — `tab$contrast` becomes NULL.

Downstream, C01d at L189 reads `iv$contrast` → NULL →
`safe_num(NULL, 0) = 0` → `d1_contrast = 0` → D1_block_strength
silently collapses to `0.35 * squareness + 0.25 * sharpness` on every
candidate, every chromosome. D1 has the highest weight in C01d's
final_score formula (0.14); the tier-count threshold (d1 ≥ 0.30)
becomes nearly impossible to hit without the contrast contribution.

**Verification:** reproduced the NULL behaviour in a fresh R 4.3.3
session with data.table 1.14.10:

```r
a <- data.table(block_id=1:3, contrast=c(0.1,0.2,0.3))
b <- data.table(block_id=1:3, contrast=c(0.15,0.25,0.35), shape="x")
m <- merge(a, b, by="block_id")
# m$contrast == NULL
# names(m) == c("block_id", "contrast.x", "contrast.y", "shape")
```

**Fix:** drop D01's `contrast` before the merge. D08's flank-based
variant-aware contrast is the one C01d was designed to read. D01's
`height` and `bg_sim` columns still survive for anyone wanting the
global-background contrast.

### FIX 16 — NN sweep tree never merged (DESIGN, manuscript-relevant)

**File:** `inversion_modules/phase_2_discovery/2d_candidate_detection/run_all.R` Phase 8

**Problem:** D09 builds an NN-sweep interval tree with the `nn_birth`
column (the coarsest NN at which each node appears separately —
high values indicate strong structure). Phase 5 saves it to
`nn_tree_<chr>.tsv` but Phase 8 never merges it onto `tab`. C01d
L210-212 reads `iv$nn_birth` for its D3 score:

```r
nn_birth <- safe_num(iv$nn_birth, 0)
nn_birth_score <- pmin(1, nn_birth / 200)
d3 <- 0.30 * nn40 + 0.30 * nn80 + 0.40 * nn_birth_score
```

Result: `iv$nn_birth` was always NA, `nn_birth_score = 0`, D3 capped
at 0.6 instead of 1.0 even for the strongest structural blocks.

**Fix:** reciprocal-overlap join. For each block, find the tree node
with the highest reciprocal overlap; take that node's `nn_birth` if
the overlap is ≥ 0.5. If no node overlaps, leave as NA
(`safe_num(NA, 0) = 0` in C01d handles it correctly). Diagnostic
message prints how many blocks got a nn_birth assignment.

### FIX 17 — `patchiness` column collision (SILENT, manuscript-relevant)

**File:** `inversion_modules/phase_2_discovery/2d_candidate_detection/STEP_D04_interior_cv.R`

**Problem:** D04 writes `patchiness = fraction of cells below
median/2` (a tail-fraction, range 0-1). D08 writes `patchiness =
sd/mean = coefficient of variation` (a spread, typical range
0-0.3). Same name, different semantics. D04 is merged first (via
`cv_ev`), then D08 at line 412, producing `patchiness.x`
(D04 tail fraction) and `patchiness.y` (D08 CV) — `tab$patchiness`
becomes NULL.

C01d L222 reads `iv$patchiness` and computes `patch_score = pmin(1,
pmax(0, 1 - patchiness * 3))`. The `1 - patchiness * 3` scaling
assumes the CV interpretation (CV of 0.3 gives patch_score of 0.1).
With `patchiness = NULL` → `safe_num(NULL, 0) = 0` → `patch_score =
1` — D5_interior_quality artificially saturated on every candidate.

**Fix:** rename D04's `patchiness` to `below_median_frac`. D08's CV
interpretation — which is what the L222 scaling was designed for —
now survives the merge. D04's tail fraction is preserved under a
semantically clearer name. No external consumer of D04's original
`patchiness` was found (grepped all phase 2/4/5 scripts except
archive); the only downstream readers of a column called `patchiness`
all use the D08/C01d CV interpretation.

### FIX 20 — C01f triangle fallback port (SILENT, closes chat 3 loose end)

**File:** `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`

**Problem:** `build_comp_with_fallback()` at L354-369 used a
two-tier fallback: C01i registry → legacy
`triangle_sample_composition.tsv.gz`. The legacy file isn't produced
by the current pipeline (same issue as BUG 3 in C01e, already fixed
in chat 2). When the registry was empty, the file-fallback closure
at L1983 (`comp_dt[chrom == chr & interval_id == iid]`) returned
empty; the candidate was silently skipped at L1985 with a misleading
"no composition (registry empty + no file entry)" message. Tier 1/2
candidates with no C01i registry entry were dropping out of
hypothesis testing entirely.

**Fix:** port `compute_bands_for_candidate()` from C01e as
`comp_from_kmeans_fallback()` — on-the-fly k-means(3) on PC_1_
average loadings within the candidate interval. New three-tier
fallback chain: registry → k-means → legacy file. Shape of
k-means output matches `comp_from_registry` exactly (columns:
`sample`, `band`, `source`, with `band1`=REF/low-PC1,
`band3`=INV/high-PC1, intersected with caller's `sample_names`).
Call site at L1980 updated to pass both `kmeans_fallback_fn` (new,
preferred) and `file_fallback_fn` (legacy, kept for back-compat).

## Findings documented, not fixed

### FINDING 18 — D14 landscape classifier is orphaned (DESIGN, non-blocking)

**File:** `inversion_modules/phase_2_discovery/2d_candidate_detection/STEP_D14_landscape_classifier.R`

D14 is not sourced by `run_all.R` and produces no files read by
Phase 8. C01d reads `iv$category` / `iv$landscape_category` /
`iv$n_children` at L333-355 for pattern-classification, but those
columns are absent from the merged scoring table. C01d correctly
falls through to its shape_class-based fallback
(`strong_inversion` / `diffuse_inversion` / `diagonal_band` /
`noise`), which is a sensible default. This is a dead contract, not
a silent corruption.

D13 plotting at L150 has a stale import of the old filename
`14_landscape_classifier.R` which will silently no-op.

**Action:** added a header comment to D14 documenting the orphan
status and the wire-in path for a future session.

### FINDING 19 — 2e GHSL is standalone-figure-only (DESIGN, non-blocking)

**Files:** `inversion_modules/phase_2_discovery/2e_ghsl/` (all)

2e/C04 snake3 GHSL v5 is NOT wired into the main C01d scoring
pipeline. Its outputs (`snake3v5_window_track.tsv.gz`,
`<chr>.ghsl_v5.annot.rds`, etc.) are consumed only by two scripts
in `phase_5_followup/codebase_8.5_experimental/current_followup/`
(STEP38, STEP40), both of which are supplementary-figure scripts.

C01d's D10 partition_stability comes from D05's
`compute_ghsl_from_simmat()` using sim_mat column-correlation
clustering — not from C04's Clair3-phased output. D05's
`compute_ghsl_stability()` is an explicit stub (header note
2026-04-17) that was designed to read C04 outputs but was never
finished.

The 2c_precomp README's "4-layer independence framework (A=dosage,
B=SV, C=GHSL, D=genotype-breakpoint)" describes Layer C as GHSL.
That's aspirational for the main scoring path — Layer C in the
current C01d pipeline is the sim_mat-based partition stability
(from D05 via Phase 7c), not C04's phased divergence.

C04's algorithm was skim-verified and is sound: rank stability
via Spearman correlation, bimodality via density-peak count,
z-score normalization against chromosome-wide baseline.

**Action:** created `phase_2_discovery/2e_ghsl/README.md`
documenting the standalone-figure status, the absence from the
C01d scoring path, and the Layer-C integration roadmap for a
future session.

## What was audited and is clean

Verified against their respective downstream consumers:

- **D01 staircase detector** — algorithm sound (vote accumulation
  → peak detection → two-pass block construction → parent/child
  nesting). One stylistic fragility at L457 where
  `parent_height <- blocks$height[r]` relies on `r` persisting
  from an inner `for` loop's broken-out state — functionally
  correct in R but ugly.

- **D07 matrix transforms** — six variants (raw, distcorr,
  localnorm, denoised, resid_bg, edge, support) are well-defined
  and sufficiently independent. `generate_all_variants` returns a
  named list consumed by D08.

- **D08 block scoring** — shape classifier thresholds defensible
  (strong_square / diffuse_square / diagonal_band / noise /
  ambiguous / unknown). Output schema matches C01d's reader at
  L189-200 except for the `contrast` / `patchiness` collisions
  (FIX 15, FIX 17).

- **D02 NN persistence** — emits
  `candidate_id, ref_*, nn<N>_match, nn<N>_status, max_survive_nn,
  survives_nn40, survives_nn80, nn_topology`; matches C01d's
  `iv$survives_nn40/nn80` at L206-209.

- **D03 flatness** — emits
  `candidate_id, flatness, far_near_ratio, monotonicity, ...`;
  matches C01d at L215-217.

- **D05 GHSL stability** — `compute_ghsl_from_simmat()` emits
  `candidate_id, partition_stability, partition_entropy,
  dominant_partition_frac, n_consistent_windows`; matches C01d's
  D10 at L283-287. The stub `compute_ghsl_stability()` is
  explicitly unused (chat 2 header clarification verified).

- **D06 SV overlap** — emits `candidate_id, best_sv_id, sv_type,
  sv_overlap_pct, n_sv_hits, ...`; matches C01d's D7 at
  L233-235.

- **D08b local PCA** — emits `candidate_id, pc1_silhouette,
  ica_kurtosis, pc1_trimodality, ...`; matches C01d's D9 at
  L278-280.

- **D09n peel diagnostic** — `effect_class` strings match C01d
  exactly (`stable` / `revealed_child` / `weakened` /
  `disappeared` / `ambiguous`). Long-format output (one row per
  block × peel_mode) is filtered in run_all.R Phase 9 on
  `peel_mode == "L1b_chrlocal_kin"` and `"L2_local_coseg"` and
  merged onto `tab` under names `l1b_effect` / `l2_effect`,
  which C01d reads at L239-240. The `sub` function-vs-matrix
  variable shadowing at L157, L406, L461, L512 is ugly but R's
  function-lookup rules make it work — calls to `sub(...)` as a
  function skip the matrix binding.

- **D09 NN sweep tree** — output schema (`node_id, start, end,
  width, start_mb, end_mb, nn_birth, nn_death, parent_node,
  children, height, topology, classification`) clean; now
  actually merged via FIX 16.

## Parse-check

All four modified R files pass `Rscript -e "parse(file=...)"` on
R 4.3.3 / data.table 1.14.10:

- `run_all.R`
- `STEP_D04_interior_cv.R`
- `STEP_D14_landscape_classifier.R`
- `STEP_C01f_hypothesis_tests.R`

## Non-scope for this chat (unchanged from chat 3)

- Phase 3 / MODULE_5A2_breakpoint_validation — deferred to Part B
- C01a internals (science audit, separate project)
- `flashlight → sv_prior` rename of the C01d output column
- Plotting scripts (D13, D15, D17) — skipped as per handoff
- D11a/D11b per-chr kinship one-shots — skipped as per handoff

## Meta-observation carried forward from chat 3

The 4:10 material-to-polish ratio from earlier chats has now
improved: this chat is 4:2 (four real fixes, two documented
findings). FIX 15 and FIX 17 were SILENT bugs that had been
shipping wrong D1 and D5 scores on every single run since the
v9.3 merge contract was established. FIX 16 capped D3 at 0.6 for
real inversions. FIX 20 was causing silent candidate drop-out in
the C01f hypothesis test. These are the kind of plumbing-is-science
bugs that only surface under a rigorous column-contract audit — the
detector algorithm papers over the individual dimension losses
because 12 dimensions average out, but the tier-count threshold
mechanism (each dimension ≥ threshold contributes a vote) is
directly sensitive to silent-NA dimensions failing their cutoff.

## For the next session (Part B)

The Part B handoff prompt
`HANDOFF_PROMPT_next_chat_partB_2026-04-17.md` is written at the
repo root. It summarizes this chat's findings so Part B has full
context without re-reading the codebase.

The 4 Part A fixes affect the scoring table and the C01f hypothesis
tests. None of them change the phase 3 / MODULE_5A2 breakpoint
validation logic, which is Part B's scope. Part B can proceed
independently.

---

## Chat 4 continuation (2026-04-17 late) — FIX 21 + FIX 22

### What prompted this continuation

Quentin pushed back on my NN interpretation twice. I had called it
"sample-space" (wrong) then "genomic-adjacent smoothing" (also
wrong). Correct reading: NN smoothing is MDS-space k-nearest-
neighbor averaging — for each window, find its k most-similar
windows anywhere on the chromosome (smallest dmat row), average
their MDS coordinates, rebuild sim_mat.

That prompted reading D09's `match_across_scales` in full, which
revealed: the tree is a **persistence barcode over NN scale**.
`nn_birth` is the coarsest scale at which each block first appears
as its own distinct unit. It's TDA-style persistence, not just a
classifier threshold. I had been reading the code too superficially.

Reading the archived snake1 precomp revealed the second finding:
**the NN-smoothing loop was lost** when the current live precomp
was adopted. The live
`phase_2_discovery/2c_precomp/STEP_C01a_precompute.R` writes only
nn=0. All multi-scale sim_mats that D02/D09 need have been absent
on every recent run. The NN-persistence track was dead — my own
FIX 16 was a no-op because `tree` was always NULL upstream.

### FIX 21 details

Ported the NN-smoothing block from
`_archive/MODULE_5A2_Discovery_Core/snakes/STEP_C01a_snake1_precompute.R`
(L1010-1042 in the archive) into
`phase_2_discovery/2c_precomp/STEP_C01a_precompute.R` right after
the precomp RDS save (after L1824).

The smoothing loop is O(n × k) per scale per chromosome where n is
window count and k is the NN value — for n=900 windows and 8
scales, the expensive `as.matrix(dist(smoothed))` dominates at
O(n²). Expected runtime per chromosome: single-digit minutes.
Across 28 chromosomes on a SLURM array: ≈1 hour total added.

Scale ladder calibrated to D09's persistence semantics, not to an
arbitrary "more is better" principle:

- 20, 40, 80 — fine scales; a small rare inversion spanning
  ~100kb-1Mb (2-20 windows) will have nn_birth in this range
- 120, 160 — middle; typical 1-5Mb inversions have nn_birth here
- 200, 240 — coarse; hits D09's `nn_birth >= 200 = INVERSION`
  classifier threshold, and reaches C01d's D3 saturation
- 320 — ceiling; confirms the largest candidate blocks persist
  through ~16Mb of MDS smoothing before dissolving

Going beyond 320 not useful for this cohort: largest suspected
real inversion ~10Mb (NN~200). Above 320, only chromosome-scale
family-LD survives smoothing, which isn't what D09 is classifying.

### FIX 22 details

Two-part wire for 2e/C04 into C01d:

Part 1 (`run_all.R` Phase 8): Added a block between the FIX 16
nn_birth merge and the FIX 15 block-scores merge. Reads
`<ghsl-dir>/annot/<chr>.ghsl_v5.annot.rds` shard. For each block,
aggregates C04's per-window annotations across `global_window_id ∈
[block.start, block.end]`:

- `ghsl_v5_score_max` — max over window scores (peak signal)
- `ghsl_rank_stability_max` — max over window rank-stability z
- `ghsl_div_contrast_z_max` — max over bimodality-contrast z
- `ghsl_div_bimodal_frac` — fraction of windows with ≥2 density
  peaks (bimodality is the signature of INV/INV vs INV/nonINV
  splitting)
- `ghsl_pass_frac` — fraction of windows with score > 0.65 (C04's
  internal PASS threshold at L478-482)
- `ghsl_n_scored_windows` — non-NA count

Part 2 (C01d D10): when `ghsl_n_scored_windows >= 3` and
`ghsl_v5_score_max` is finite, compute `d10_ghsl = 0.50 *
ghsl_v5_score_max + 0.25 * ghsl_div_bimodal_frac + 0.25 *
ghsl_pass_frac`, then blend `d10 = 0.60 * d10_ghsl + 0.40 *
d10_simmat`. Otherwise fall back to `d10 = d10_simmat` (pre-FIX-22
behaviour, preserves backward compat). Emits `d10_source` column
so downstream can see which path fired.

Blend weights reflect that Clair3-phased haplotype divergence is
direct biological evidence for a karyotype difference, while
sim_mat partition stability is an indirect structural proxy. When
both are available, the biological signal dominates (0.60) but the
structural signal still contributes (0.40) to handle cases where
Clair3 phase confidence is marginal.

### Revised understanding of the NN tree (for next chat's benefit)

`build_nn_sweep_tree(sim_mats_by_nn)` works on a dict of
`sim_mat_nn<k>.rds` files keyed by k. For each k it runs the full
D01 staircase detector, producing a block set. Then
`match_across_scales`:

1. Takes the coarsest k as roots — each block at the coarsest
   surviving scale is a root node, with `nn_birth = coarsest_nn`.
2. Walks scales from coarse to fine. At each k, for each new
   block, find the best reciprocal-overlap match among active
   parent leaves (`overlap_thresh = 0.70`).
3. Per parent, count matches. If 0 matches → parent "disappears"
   at this k (artifact of oversmoothing). If 1 match → parent
   "stable" with small boundary refinement. If >1 matches →
   parent "splits" into children (each with nn_birth = k).
4. Unmatched new blocks that fall inside an existing parent become
   "novel" children; otherwise become new roots with nn_birth = k.

So `nn_birth` is the coarsest k where that feature first
distinguishes itself as its own unit:

- A feature with `nn_birth = 320` is a root at the top scale —
  the block persists through the heaviest MDS-smoothing and is
  the strongest, most persistent structural unit.
- A feature with `nn_birth = 80` either (a) emerged as a novel
  block only at k=80 (didn't exist as a distinct unit at coarser
  scales — smaller / more local feature), or (b) is a child born
  from a split at k=80.
- A feature that "disappears_nn40" was visible from k>=40 but
  dissolves at k=20 (oversensitivity to noise at the finest
  scale).

C01d's scoring `nn_birth_score = pmin(1, nn_birth / 200)` gives
full credit to features persisting ≥200 scale deep, which is the
`INVERSION` class boundary in D09's own classifier. With
pre-FIX-21 scales capped at 80, nothing could score above 0.4 on
this component.

### Updated scope for Part B handoff

The Part B handoff prompt has been updated to include FIX 21/22
status and the precomp-rebuild dependency. Part B's phase 3 audit
is independent of these phase 2 changes.

