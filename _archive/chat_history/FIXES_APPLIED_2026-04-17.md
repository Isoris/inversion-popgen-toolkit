# FIXES_APPLIED_2026-04-17.md

This session's fixes applied to the working copy at /home/claude/work/fixed/.

Severity order (fixed first):
  CRASH  (BUG 5)  — C01e Panel F undefined gt_data
  SILENT (BUG 9)  — C01g TE self-match bug
  SILENT (BUG 27) — C01d test05 rename not propagated
  SILENT (BUG 3)  — C01e triangle_sample_composition obsolete read
  STALE  (BUG 1)  — C01d --flashlight_dir dead flag
  STALE  (BUG 4a) — C01e --repeats dead flag
  STALE  (BUG 4b) — C01e --het_dir dead flag  
  STALE  (BUG 6)  — C01g --ref_fasta dead flag
  META  (ISSUE 7) — C01g fossil detection circular dep → archive + remove
  META  (BUG 2)   — 2d D05 GHSL stub → wire to real GHSL or document NA
  POLISH (BUG 14) — PHASE_01C cor() guard for short intersect
  POLISH (C01e meta cols) — align metadata intersect with v9.3 columns

## Chat 4 additions (Part A — 2d/2e/C01f)

  SILENT (FIX 15) — run_all.R Phase 8: D01/D08 contrast collision → contrast.x/contrast.y → C01d D1 silently 0
  DESIGN (FIX 16) — run_all.R Phase 8: NN sweep tree nn_birth not merged → C01d D3 capped at 0.6
  SILENT (FIX 17) — STEP_D04_interior_cv.R: patchiness renamed to below_median_frac to stop collision with D08's CV patchiness
  SILENT (FIX 20) — STEP_C01f_hypothesis_tests.R: on-the-fly k-means fallback ported from C01e (closes chat 3 loose end)
  DOC    (FINDING 18) — STEP_D14_landscape_classifier.R: orphan status documented in header
  DOC    (FINDING 19) — phase_2_discovery/2e_ghsl/README.md: standalone-figure status documented

## Chat 4 continuation — FIX 21, FIX 22 (critical for tomorrow's precomp rebuild)

  CRASH  (FIX 21) — STEP_C01a_precompute.R: NN-smoothed sim_mats were never produced by the live precomp (the loop was only in the archived snake1 precomp). Entire D02/D09 NN-persistence track silently dead. Ported the MDS-space kNN smoothing loop; expanded scale ladder to 20,40,80,120,160,200,240,320 (from 20,40,80). ACTION REQUIRED: rebuild precomp for all 28 chromosomes before running 2d detection.
  DESIGN (FIX 22) — run_all.R Phase 8 + C01d D10: wired 2e/C04 GHSL v5 phased-Clair3 signal into the scoring table. Six GHSL columns now merged per block; C01d's D10 blends ghsl_v5_score_max with D05's sim_mat-based partition_stability (0.60 ghsl / 0.40 simmat when both available). d10_source column emitted for transparency. FINDING 19 superseded.

## Chat 5 additions (Part B — phase 3 audit + registry wiring)

  DESIGN (FIX 23)     — breakpoint_validator_standalone.py L888–896: INV-vs-HET Fisher block tupled 2-tuple return into first slot of a 3-tuple, making the else-branch unreachable. Computed values never consumed (docstring promises 3 comparisons, comparisons list has 2). Simplified to straight-line code + TODO.
  CRASH  (FIX 24)     — breakpoint_validator_standalone.py L648: unguarded y/t*100 in report writer crashes on empty groups. Added `if t > 0 else "–"` guard.
  CRASH  (FIX 25)     — breakpoint_validator_standalone.py L1182: SLURM wrapper template referenced breakpoint_validator.py (missing _standalone suffix). Users hit "No such file or directory" on cluster.
  CRASH  (FIX 26)     — run_breakpoint_validation.sh L45/L106: launcher passed matched_delly_snake_candidates.tsv as --candidates but STEP01 writes matched_inv_candidates.tsv. Pipeline failed at STEP02 on first run. Chat 4's "architecturally complete" assessment missed this. Fixed both launcher calls + script docs.
  SILENT (FIX 27)     — 03_statistical_tests_and_seeds.py L222: seed_summary.tsv column drift when qualifying status varies across candidates (n_seeds_ref conditionally added). Initialize all six seed keys upfront.
  CRASH  (FIX 28)     — 05_delly_manta_concordance.py L138: unguarded division by n_manta_total inside `if n_delly_total > 0` crashes when DELLY has calls but Manta has none. Split into four branches.
  DESIGN (FIX 29 v2)  — WIRING: phase 3 now writes into the evidence registry. STEP03 writes existence_layer_d blocks per candidate → q7_layer_d_fisher_or/_fisher_p/_armitage_z/_armitage_p/_n_inv_with_support/_n_inv_total flat keys that C01f's compute_group_validation() and compute_candidate_status.R already read for the VALIDATED promotion gate. STEP06 writes new existence_layer_b_bnd_rescue blocks per paired BND junction → q7b_bnd_rescued/_rescue_source/_pair_bp1/_pair_bp2/_match_type/etc flat keys that supplement Layer B evidence for candidates reconstructed from CT=3to3+CT=5to5 junction pairs. New schema: registries/schemas/structured_block_schemas/existence_layer_b_bnd_rescue.schema.json. Both STEP03 and STEP06 no-op gracefully if REGISTRIES_ROOT is unset. Launcher + config + INDEX_remaining_blocks updated. Earlier first-draft FIX 29 (a D13 dimension in C01d) was reverted mid-session — wrong abstraction; Layer D is not a sub-dimension of Layer A.
  STRUCT (chat 5)     — Flattened inversion_modules/phase_3_refine/MODULE_5A2_breakpoint_validation/ one level up. 11 files now in phase_3_refine/ directly; wrapper dir removed. Eight documentation files updated with corrected paths.
  DOC    (chat 5)     — Rewrote phase_3_refine/README.md with full 4-layer model, Layer D + BND-rescue contract tables, downstream consumer list, biological caveats. Corrected earlier chat-5 edits in phase_4/README.md, 4a_existence_layers/README.md, and 2c_precomp/README.md that had (prematurely) declared the phase-3→phase-4 contract inactive — it's now active via FIX 29 v2.
  TEST   (chat 5)     — test_registry_sanity.py: added Layer D and BND-rescue fixtures mirroring what STEP03 and STEP06 emit.

## Chat 6 additions (Part C — phase 4 consumer-side audit)

Note: audit-only session. All fixes below are **proposed, not applied** — flagged for a short follow-up session.

Also note: during chat 6, Quentin clarified that the registry library itself is not yet implemented. Findings that looked like reader/writer parity bugs (the `*_detected`/`*_tested` flag family; the BND-rescue key family) are therefore intentional forward-declarations, not bugs. Only bugs that are real independent of lib status are listed here.

  CRASH  (FIX 30 proposed) — 03_statistical_tests_and_seeds.py L373/L375/L377: seed writer references r['total_score'] but the evidence TSV has no such column. KeyError fires for every candidate that qualifies for seeding. Chat 5's FIX 27 fixed the adjacent header drift but didn't exercise the qualifying branch past the header. Best one-liner: `score = int(r.get('pe',0)) + int(r.get('sr',0))` — but Quentin should confirm the column choice since seed files feed C01i's (still-deferred) seeded-init.
  CRITICAL (FIX 31 proposed) — STEP_C01f_hypothesis_tests.R L2322–2323: live call to compute_group_validation() hard-codes layer_d_fisher_p=NA, layer_d_fisher_or=NA. The VALIDATED-promotion gate at L525–526 requires is.finite() on both → unreachable. Fix: two tryCatch blocks reading reg$evidence$get_evidence(cid, key) with NA fallback, pattern already in use two lines above for q6_group_validation. See AUDIT_LOG_chat6 Finding 1 for the exact snippet.
  HIGH   (FIX 32 proposed) — 4e compute_candidate_status.R L272/L298/L309 and 4c group_validation_gate.R L50/L52: Pathway A gates gate on `fisher_p < 0.05` only, missing the `fisher_or > 5` clause that the Layer D schema description and C01f's compute_group_validation() both require. Once FIX 31 lands, 4e and C01f will disagree on significant-but-low-OR candidates. One-line read of fisher_or + one-line gate update per site.

Empirical verification of chat 5 FIX 29 v2 (this session):
  Ran STEP03 on fake strong-signal evidence (20/20/20 REF/HET/INV with 10%/40%/90% support fractions) against a scratch registry. Block written cleanly, 8 keys extracted into keys.tsv, validation_status=validated. Silent standalone fallback verified (wrong --registries_root → skip message → continues). Chat 5's write-side wiring confirmed correct end-to-end.

## Chat 7 additions (Part D — apply chat-6 fixes + begin 4a audit)

  CRASH    (FIX 30)        — phase_3_refine/03_statistical_tests_and_seeds.py L373–377: seed writer referenced r['total_score'] but the evidence TSV has no such column. KeyError crashed every candidate that qualified for seeding. Replaced with defensive per-sample evidence score: int(r.get('pe',0) or 0) + int(r.get('sr',0) or 0). Empirical smoke test confirms strong-signal candidate now qualifies + writes seed file cleanly (REF=0, HET=2, INV=5 evidence scores). Python AST-parse OK.
  CRITICAL (FIX 31)        — 4c_group_validation/STEP_C01f_hypothesis_tests.R L2316–2325: replaced `layer_d_fisher_p = NA, layer_d_fisher_or = NA` hard-codes with `reg$evidence$get_evidence()` tryCatch reads. Pattern mirrors the `q6_group_validation` / `q6_validation_promotion_cap` reads two lines above. NA fallback preserves behaviour when STEP03 hasn't run. Also updated `patches/01_C01f_comp_from_registry.R` template so re-applying patches doesn't re-introduce the bug. Makes the VALIDATED promotion gate reachable for the first time.
  HIGH     (FIX 32)        — 4e_final_classification/compute_candidate_status.R L272 + L298 + L309 AND 4c_group_validation/group_validation_gate.R L50–59 + L89 + L96: added `fisher_or` reader (default 0 = conservative) and `fisher_or > 5` clause to both Pathway A gates in 4e (4-layer and 3-layer) and to `or_passed` in group_validation_gate. Reason strings updated from `OR_p=...` to `OR=..., p=...`. Now consistent with Layer D schema description and C01f VALIDATED gate (L525–526 of C01f).
  CRITICAL (FIX 33)        — 4a_existence_layers/STEP_C01d_candidate_scoring...R L826–840: added tier recomputation after the dim_positive recompute. Prior bug: tier was assigned inside the per-candidate loop using d11=d12=0 defaults, then dim_positive got recomputed with the real D11/D12 but tier stayed frozen. Every candidate with non-zero D11 or D12 had its tier systematically deflated; dim_positive and tier columns were internally inconsistent in the master catalog. Fix re-applies the ≥8/≥6/≥4 tier ladder AND the peel-disappeared downgrade after the dim_positive recompute. Guards the peel downgrade with column-existence check.
  COSMETIC (FIX 34)        — 4a_existence_layers/STEP_C01d_candidate_scoring...R L32–51 header + L412 inline: stale "10 scoring dimensions" / "Tier 1: ≥7/10" was describing pre-v9.3.2 C01d; live code is 12 dimensions with ≥8/≥6/≥4 threshold ladder. Updated both.

Phase 4a audit notes (not bugs, pending lib work):
  - Finding 8: No 4a script writes registry blocks (C01d + C01e: zero writes; C01g: attempts writes via `utils/registry_key_helpers.R` which doesn't exist yet). Expected under lib-not-yet-implemented framing; flagged for lib chat.
  - Finding 9: C01g silent-skip when helpers file is missing is invisible to operator at run time. Suggest WARNING message when .bridge_available is TRUE and helper is missing. Flagged for lib chat.

Phase 4a audit remaining (deferred to chat 8):
  - C01g (1428 lines) — largest script in 4a, not yet audited
  - C01e (694 lines) — figures only, low priority
  - End-to-end contract check — C01g output columns vs C01d D11 reader

## Chat 7 phase 4a audit continuation (post-initial writeup)

  SILENT (FIX 35) — 4a_existence_layers/STEP_C01g...R L465–481: dedup aggregation produced `candidate_ids` (semicolon-joined list) but never `matched_core_id` (single-value canonical id). Registry-write block at L1407 guards on `"matched_core_id" %in% names(confirmed)` → always FALSE → per-boundary registry writes never executed (second reason for "0 evidence keys" beyond Finding 9's missing helper). Fix computes `matched_core_id` from best_idx row's candidate_id with two fallback levels. Prepares C01g so when the lib work from Finding 9 lands, the registry writes will actually run.

## Phase 4a audit — verdict (chat 7)

  Flat-TSV catalog (candidate_scores.tsv.gz): ready for phase 4b/c/d/e consumption after FIX 33 corrected tier values.
  Unified boundary catalog (boundary_catalog_unified.tsv.gz): ready for C01d D11 consumption; column contract verified.
  Registry blocks from 4a: not yet produced (Finding 8 / Finding 9 = lib-pending). When the lib work lands, three items need to land together: (1) utils/registry_key_helpers.R with register_C01g_boundary + store_C01g_boundary, (2) Layer A block writer in C01d, (3) loader extension for {from_derived: "block_exists"} to synthesize q7_layer_*_{detected,tested} flags (chat 6 Finding 3). FIX 35 prepares item (1)'s landing.

## Chat 8 additions (Part E — phase 4b audit)

Posture: fix-as-we-go with a three-fix ceiling (per Quentin's handoff). Stop at three and flag the rest. Three applied below; one README-drift + one forward-declared-keys note flagged.

  COSMETIC (FIX 36) — 4b_group_proposal/STEP_C01i_b_multi_recomb.R L316: inline comment claimed decision rule was `(S1 AND S2) OR S1 OR S3` but code implements `(S1 AND S2) OR S3 OR (S1 AND mosaic>100kb)` (matches README §"Combined decision (4b.2)"). Code is correct; comment was misleading for auditors. Rewrote comment to describe the three-term rule and added the noise-tolerance rationale for why S1-alone is insufficient for short mosaics but sufficient for long ones. No runtime behaviour change.

  CRASH    (FIX 37) — 4b_group_proposal/STEP_C01i_b_multi_recomb.R L233 (function flashlight_hemi_signal): `sample_dels <- dels[sample_id %in% unlist(het_carriers)]` referenced an undefined object `het_carriers`. Copy-paste from the cheat2 pattern in STEP_C01i_decompose.R L290–294 (`bp_dels_left$het_carriers`) that was never adapted for per-sample filtering in the Signal 3 hemizygous detector. Function is called per-sample inside the main loop; `sample_id` is a scalar parameter, not a column. Bug is currently dormant because flashlight isn't universally available on LANTA, but activates as a hard crash ("object 'het_carriers' not found") the moment `get_internal_dels()` returns any rows with flashlight loaded. Fix: per-row mask using the scalar sample_id against each row's `het_carriers` list (same list-column convention as patches/patch_C01i_decomposition_flashlight.R L269), with tryCatch fallback to all-FALSE.

  SILENT   (FIX 38) — 4b_group_proposal/lib_decompose_helpers.R L88–89 (extract_pc_loadings): `%||%` binds tighter than `+`/`-` in R precedence, so `wins$start_bp %||% wins$mid_bp - 50000L` parsed as `(wins$start_bp %||% wins$mid_bp) - 50000L`, subtracting 50 kb from every window start even when `start_bp` was present. Same for `win_ends + 50000L`. Net effect: every window interval systematically inflated by 100 kb. Downstream, every `mosaic_length_bp` in multi_recomb was inflated by 100 kb, which (a) falsely promoted true-short mosaics past the 100 kb Signal-1-alone threshold → spurious RECOMBINANT calls, and (b) against the real cheat24 thresholds (50 kb / 200 kb from `flashlight_v2/cheats/cheat24_recombinant_prior.R`), most true gene-conversion events (tens of kb) were lost to "ambiguous" classification and some intermediate mosaics were mis-classified as double_crossover. **Verified by Quentin's pushback:** the precomp (C01a L849) always writes real `start_bp`/`end_bp` columns, so the `%||%` fallback branch was never supposed to fire and the 50000 magic number was arbitrary guesswork for dead code (unrelated to C01a's actual multiscale window ladder 20/40/80/120/160/200/240/320). **Fix (v2):** strip the fallback entirely and use the real columns directly (`win_starts <- wins$start_bp`). Grep confirmed no other instance of the same precedence pattern exists in the tree, and no caller in the tree produces a precomp without these columns.

Phase 4b audit notes (not bugs):

  - Finding S (flag only, per handoff): STEP03 seed writer at `04_deconvolution_seeds/{inv_id}_seeds.tsv` (fixed in FIX 30) produces seed files but no 4b script reads them. The only seeding path wired today is flashlight-based (STEP_C01i_decompose.R L266–301: `load_flashlight` → `get_sv_anchors` → `seeded_kmeans`). STEP03 seeds as designed k-means init for C01i remain deferred per chat 5. Needs a design call for the lib / seed-wiring chat. C01i's current flashlight-only approach works when flashlight is available; when it's not, falls back to unsupervised k-means — STEP03 seeds would be a useful middle path but requires deciding priority order (STEP03 > flashlight? flashlight > STEP03? combine?).

  - Finding T (README drift, lower priority than lib work): 4b_group_proposal/README.md §"Composite detection (4b.3)" lists four structure_type values (`single_block`, `two_block_composite`, `gradient`, `fragmented`) but nested_composition_core.py::classify_structure actually returns seven (`too_sparse`, `homogeneous`, `dominant_plus_secondary`, `two_block_composite`, `multi_block_fragmented`, `continuous_gradient`, `diffuse_mixed`). Seal's downstream logic references the real names (`two_block_composite` L139, `multi_block_fragmented` L132), so the contract between 4b.3 and 4b.4 is intact — the drift is README-only. Flag for future doc pass.

  - Finding U (forward-declared flat keys, parallels chat 7 Finding 8): seal writes `q1_composite_flag` (L359) and `q2_decomp_quality_flags` (L363) as registry evidence, but no live script reads these keys. They're forward-declarations for future 4c/4d/4e consumers. `q6_validation_promotion_cap` and the four `q6_*` / `q7_*` keys written by seal ARE consumed (by C01f post-FIX-31 and by compute_candidate_status.R). This is the same forward-declaration pattern as chat 7 Finding 8; consistent with lib-not-yet-implemented framing.

  - Finding V (cheat24 threshold divergence, surfaced while tracing FIX 38 impact): the real cheat24 at `flashlight_v2/cheats/cheat24_recombinant_prior.R` uses `MOSAIC_SHORT_BP=50000`, `MOSAIC_LONG_BP=200000` (gene_conversion <50kb, double_crossover ≥200kb). The inline fallback `classify_recombinant_event` in multi_recomb L99–121 uses 100000 / 500000. Same candidate gets different event_class depending on whether flashlight is loaded at runtime. Sensitivity/reproducibility issue. Not fixed this chat — needs Quentin to confirm which threshold pair is authoritative (likely the real cheat24) before aligning the fallback.

  - Doc drift in decompose header: `STEP_C01i_decompose.R` L16 header comment promises the output includes a "per-window dosage track". The code actually produces a per-window **class** track (HOM_REF/HET/HOM_INV assignments from PC1 k-means), not a dosage track. Batch with Finding T in the future doc pass.

## Phase 4b audit — verdict (chat 8)

  Group-registration contract: VERIFIED. Seal L254 registers groups as `inv_{cid}_{HOM_REF|HET|HOM_INV|RECOMBINANT}` matching exactly the gids that C01f::comp_from_registry() expects (L307–315, accepting HOM_REF or HOM_STD alias). HOM_STD backward-compat alias also registered at L270.

  Quality-flag contract: VERIFIED. Seal's determine_validation reads `silhouette_score`, `bic_gap_k3_vs_k2`, `phase_concordance` from the internal_dynamics block; decompose writes all three at L407–409 with matching names. Seal's quality_flags output (normal / composite_cap / low_silhouette / weak_k3_fit / low_phase / maybe_composite_flagged) is emitted but as noted in Finding U has no live reader yet.

  Structure-type contract: VERIFIED. nested_composition_core.py::classify_structure returns the seven values listed above; seal::resolve_final_classes matches against `two_block_composite` and `multi_block_fragmented` directly, other values pass through as pca_class.

  Recombinant-detection contract: VERIFIED post-FIX-37 + FIX-38. Signal 1 / Signal 2 / Signal 3 combination rule now matches the README (post-FIX-36 comment). Mosaic length calculation now accurate (post-FIX-38). Signal 3 no longer crashes when flashlight data is present (post-FIX-37). cheat24 event classification now receives accurate mosaic_length_bp inputs.

  Tier flow: VERIFIED. All three main 4b scripts filter on `tier <= 3` at their main loops (decompose L215, multi_recomb L264, seal L305). Post-FIX-33 tier values in candidate_scores.tsv.gz are now D11/D12-aware and internally consistent with dim_positive. All candidates passing the old (buggy) tier <= 3 filter also pass the new one (tier 3 is the bottom); the tier 1/2/3 distribution shifts but the set of candidates entering phase 4b is unchanged.

  STEP03 seed integration: STILL DEFERRED (Finding S). No behaviour change needed this chat.

  Overall: Phase 4b is ready to feed phase 4c cleanly. FIX 37 and FIX 38 are the two correctness-meaningful fixes this chat (FIX 36 cosmetic). Remaining items (README drift, STEP03 seed wiring, forward-declared key consumers) are deferrable.


## Chat 9 (2026-04-17)

### FIX 39 — SILENT — T9 jackknife vocabulary mismatch (C01f)
**File:** `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`
**Lines:** L499 (classifier branch), L546 (T8+family promotion)
T9 emits `few_family_contributing` for 2–3 family case; C01f classifier
matched on `few_family`. Every such verdict fell through whole classifier:
`family_linkage=unknown`, no cap, no quality flag. Fixed both call sites
to match `few_family_contributing`.

### FIX 40 — ARCHITECTURAL — composite_flag key consolidation (Option A+B)
**Files:**
- `inversion_modules/phase_4_postprocessing/schemas/internal_ancestry_composition.schema.json` L69
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R` L359–362
Schema extracted `q1_ancestry_composite_flag`; seal explicitly wrote
`q1_composite_flag`. Two names, same semantic value, same source
(`nst$data$composite_flag`). Fixed by renaming schema key to
`q1_composite_flag` (canonical) + removing seal's explicit write.
Single writer, single key name.

### FIX 41 — CRASH (active) — cheat30 alias scope bug (4d)
**File:** `inversion_modules/phase_4_postprocessing/4d_group_dependent/cheat30_gds_by_genotype.R`
**Lines:** L66 (before) → L100 (after)
`compute_pairwise_ibs <- compute_pairwise_gds` alias was inside function
body (local frame), not module scope. `run_cheat30` L222 calls
`compute_pairwise_ibs` — would crash "could not find function" on every
invocation. Fixed by moving alias outside the function definition.

### FIX 42 — ARCHITECTURAL (large) — 4e build_key_spec aligned to v2 spec
**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/compute_candidate_status.R`
**Lines:** L40–233 (`build_key_spec` rewritten), L560–630 (`compute_completion` reworked)
- Expanded from 315 to 367 keys. Intentional v9+v10 dual aliases
  (HOM_REF+HOM_STD, jackknife_verdict+jackknife_status) account for
  the overshoot versus the 352 v2 target.
- Added `*_aspir` per-Q lists marking keys whose writers don't yet
  produce them.
- `compute_completion()` now respects aspirational marking:
  `applicable = total - aspirational - not_assessed`; adds
  `pct_of_spec` as secondary metric.
- Fixed Q6 naming drift: schema writes `q6_hwe_p`/`q6_hwe_chi2`/
  `q6_hwe_verdict`/`q6_genotype_balance`/`q6_polymorphism_class`/
  `q6_freq_class`, none of which were in the old spec. Old spec had
  `q6_hwe_chisq_p` (wrong name) and `q6_n_HOM_STD` only (no HOM_REF).
  Both naming conventions now included.

### FIX 43 — SILENT (impactful) — jackknife_status readers → q6_family_linkage
**Files:**
- `inversion_modules/phase_4_postprocessing/4c_group_validation/group_validation_gate.R` L64–83
- `inversion_modules/phase_4_postprocessing/4e_final_classification/characterize_candidate.R` L448–459
`q7_t9_jackknife_status` has 4 readers, zero writers. Gate.R falls
through to OR-only logic → caps every candidate without Layer D at
UNCERTAIN regardless of jackknife robustness. Fix: prefer
`q6_family_linkage` (reliably written by seal + C01f). Canonical
fragile signal is now `family_linkage == "pca_family_confounded"`.
`has_groups` check also extended to accept `q6_n_HOM_REF` (v10 name)
alongside `q6_n_HOM_STD` (v9 alias).

### FIX 44 — CRASH (dormant) — characterize_candidate source()
**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/characterize_candidate.R`
**Lines:** L1–48 (header rewritten)
`characterize_candidate()` L517 calls `assess_group_validation` which
lives in `4c_group_validation/group_validation_gate.R` — never sourced
at module load. Would crash "could not find function". Fix: added
four-candidate path locator for the gate file at module load, emits
warning if not found, supports `GROUP_VALIDATION_GATE` env-var override.


### FIX 45 — FEATURE — run_characterize.R driver (closes Finding AF)
**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/run_characterize.R`
**Lines:** NEW, 420 lines
The missing 4e entry point. `characterize_candidate.R` is a function
library with no driver; past invocations were theoretical. This driver
does sibling discovery (four-path fallback for gate.R + cwd fallback
for characterize_candidate.R, env-var override, function-in-scope
assertion), three-tier loader mirroring `compute_candidate_status.R`
L727–827 (per-candidate *.keys.tsv → merged evidence_registry.tsv →
candidate_scores.tsv.gz), crash-proof main loop (empty keys → all-EMPTY
row; characterize throwing → tryCatch with error-logged EMPTY row),
and four outputs: `characterization.tsv` (26 columns), per-candidate
`characterization.txt` via `format_characterization()`,
`characterization_summary.txt` (per-Q status matrix + group distribution
+ top-25 contradicted + top-15 labels + error log), and conditional
`candidate_full_report.tsv` merging with `candidate_status.tsv` on
`candidate_id` when present. Parse OK; 3-candidate synthetic smoke
passes (healthy/artifact/empty cases all produce semantically correct
output).

### FIX 46 — CRASH (active) — characterize_candidate.R %||% ordering
**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/characterize_candidate.R`
**Lines:** L23–71 (reordered)
Chat 9's FIX 44 added a gate-locator block that used `%||%` at L28
before it was defined at L49. Every `source(characterize_candidate.R)`
crashed with "could not find function %||%" before any Q* function could
load. Also `sys.frame(1)$ofile` only resolves when sourced via
`source(file, chdir=TRUE)` — crashed otherwise. Fix: helpers
(`%||%`, `safe_num`, `safe_bool`, `has_key`) moved to top; gate-locator
now below them; `sys.frame(1)$ofile` read wrapped in tryCatch with null
fallback; the ofile-based candidate path only included when ofile is
non-empty. Parse-check alone couldn't find this; only runtime invocation.

### FIX 47 — CRASH (active) — safe_num zero-length failure
**File:** `inversion_modules/phase_4_postprocessing/4e_final_classification/characterize_candidate.R`
**Lines:** L24–40 (helper definitions)
`safe_num(x, d)` had `if (is.finite(x)) x else d` with no length guard.
For any missing key: `keys[[k]]` → NULL → `as.numeric(NULL)` →
`numeric(0)` → `is.finite(numeric(0))` → `logical(0)` → `if (logical(0))`
throws "argument is of length zero". Every characterize_qN function hit
this on its first check of an optional absent key. Q1/Q4 happened to
return EMPTY early and dodged it; Q2/Q3/Q5/Q6/Q7 crashed.
**This is why the 4e pipeline had never run end-to-end before chat 10.**
Fix: `safe_num`, `safe_bool`, `has_key` all now zero-length-safe
(null/empty check, scalarise input with `x[1]` before the is.finite /
is.logical branch). Verified by smoke test that all 7 Q* functions now
produce correct output.

### FIX 48 — FEATURE — run_phase4b.sh emits PHASE4B_SEAL_JID trailer
**File:** `inversion_modules/phase_4_postprocessing/orchestrator/run_phase4b.sh`
**Lines:** +6 lines at end
For parent `run_phase4.sh` (FIX 49) to gate 4c on the 4b.4 seal job, it
needs to know the seal jobid. Added
`echo "PHASE4B_SEAL_JID=${JID_SEAL}"` as the last stdout line.
Parent parses via `grep '^PHASE4B_SEAL_JID=' | tail -n1 | cut -d= -f2`.
Verified via fake-sbatch smoke — real runs produce one clean line.

### FIX 49 — ARCHITECTURAL — run_phase4.sh calls run_phase4b.sh sub-DAG
**File:** `inversion_modules/phase_4_postprocessing/orchestrator/run_phase4.sh`
**Lines:** L110–127 (replaced 3 lines with ~25)
Old orchestrator referenced `LAUNCH_C01i_decomposition.sh` as a single
4b job. 4b is now four scripts with its own DAG (decompose →
multi_recomb, nested_comp parallel, seal gates on all three) managed
by `run_phase4b.sh`. Fix: replaced the `LAUNCH_C01I` submit_job call
with `bash run_phase4b.sh`, captured stdout via `tee /dev/stderr`,
greps `PHASE4B_SEAL_JID` from the FIX-48 trailer. `JID_DECOMP` now
holds the seal jobid so the 4c `--dependency=afterok:${JID_DECOMP}`
gate still works. Dropped unused `LAUNCH_C01I` variable. Dry-run
forwarding via `${DRY_FLAG}`. Exits code 3 if grep empty (sub-orch
failed before emitting). Both scripts pass `bash -n`.

### FIX 50 — FEATURE — production LAUNCH_group_cheats.sh (closes Finding Y)
**File:** `inversion_modules/phase_4_postprocessing/orchestrator/LAUNCH_group_cheats.sh`
**Lines:** NEW file, renamed/corrected from `LAUNCH_group_cheats_example.sh`
`run_phase4.sh` L115 expects `launchers/LAUNCH_group_cheats.sh` (no
`_example` suffix). Example had stale paths (`cheats/` instead of
`4d_group_dependent/`) and missed cheat27/28/29. Fix: (1) corrected all
three existing cheat paths to `4d_group_dependent/`; (2) added cheat27
(SD/NAHR substrate), cheat28 (tandem-repeat context), cheat29 (junction
forensics) calls — all NONE-gated; (3) left aspirational scripts
(`compute_age_fst_subblock.R`, `STEP_C01f_c_burden_regression.R`)
as commented-out calls with SUPPORTED / VALIDATED gates still active
for diagnostic messages — uncommenting is a one-line change when
scripts arrive. Marked executable. Parse OK (`bash -n`).


## Chat 11 — registry API extensions

### FIX 51 — FEATURE — Interval registry extensions (9 methods)
**File:** `registries/api/R/registry_loader.R`
**Lines:** ~180 new in `load_intervals_api()` (L211–L434 region)

Added 9 methods supporting nested/overlapping/cross-chrom topology:
`get_children`, `get_parent`, `get_ancestors` (cycle-safe),
`get_descendants` (BFS, depth-capped), `get_overlapping` (with
`exclude_cid`), `get_nested_within` (coord-based, independent of
parent_id), `classify_relationship` (5-value enum),
`update_candidate` (partial row update, recomputes size_kb if bounds
changed), `bulk_add_candidates` (batch insert with dup-skip, returns
count). Read-fresh-per-call convention preserved. Internal helper
`read_cands()` factored out to avoid duplication across 10+ methods.

Parse-clean (`parse("registry_loader.R")` → 9 top-level exprs, same
count as pre-fix; append-only within the same closure).

### FIX 52 — FEATURE — Sample registry extensions (10 methods)
**File:** `registries/api/R/registry_loader.R`
**Lines:** ~200 new in `load_samples_api()` (L143–L434)

Refactored `load_samples_api()` to use a shared `build_extended()`
factory — real-path branch and shim branch both call it. Added 10
methods: `get_sample_metadata`, `get_sample_groups` (scans
`sample_groups.tsv` + group member files), `get_family`/`list_families`/
`get_family_members` (with warn-once guard when `family_id` column
absent), `list_carriers_across_candidates` (uses `subgroup` column when
present, else name-suffix fallback), `compute_group_overlap`,
`list_recombinant_candidates` (parses `inv_<cid>_RECOMBINANT[_GC|_DCO]?`),
`get_groups_for_candidate` (returns all 7 status groups in one call —
directly collapses `comp_from_registry` in C01f), `find_co_segregating_groups`
(preloads members to avoid O(N²) file reads, restricts to `inv_` prefix
by default).

Also fixed pre-existing Finding AJ in the same refactor: `get_master`
was aliased to `full$list_groups` (wrong table).

### FIX 53 — FEATURE — Query API composite extensions (3 methods)
**File:** `registries/api/R/registry_loader.R`
**Lines:** ~130 new in `load_query_api()`

`nested_family_tree(root_cid)` — recursive nested list with cycle
guard via `seen` set; returns NULL for unknown root.
`overlap_clusters(chrom)` — union-find via BFS over overlap edges;
enumerates chrom candidates via `get_overlapping(chrom, 1, MAX_INT)`
rather than requiring a new `list_candidates_on_chrom` method.
Returns data.table with `cluster_id` and `cluster_size`.
`sample_inversion_load(sid)` — parses groups from `get_sample_groups`,
breaks down by status (longer status prefixes first to avoid greedy
match), counts per-status plus total.

### FIX 54 — FEATURE — Python mirror (19 atomic methods)
**File:** `registries/api/python/registry_loader.py`
**Lines:** ~400 new (403 → 812)

Full interval API mirror (9 methods with `_has_parent` helper and
`_read_cands` internal). Full samples API mirror of 10 methods
(including warn-once `_fam_warned` class attr). Parity with R API
confirmed via smoke test against identical fixtures. Python intentionally
does NOT mirror the 3 composite query methods per handoff direction
(they live in R because they compose interval + sample API calls).

### FIX 55 — BUG — sample_registry.R column-class drift (Finding AK)
**File:** `utils/sample_registry.R`
**Lines:** L149 patched; +7 chars and +5 comment lines

`fread(groups_file)` without `colClasses` auto-detects `created`
column as POSIXct on re-read. `format(Sys.time(), "%Y-%m-%d %H:%M:%S")`
produces character. `rbind` of mixed-class columns fatals with
"Class attribute on column 8 of item 2 does not match". Fix:
`fread(groups_file, colClasses = c(created = "character"))` on the
one line that feeds into `rbind`. Other three `fread` call sites are
read-only queries (`list_groups`, `has_group`, `get_group`) and don't
rbind — left unchanged to minimize diff. Annotated in-file with
2026-04-17 Finding-AK explanation.

Blocking-severity for chat 11 smoke testing (fixture registers 6+
groups per session) and for chat 13 (seeded-region bulk writer will
hit this hard).

### FIX 56 — CLEANUP — delete orphan duplicate files (Findings AM, AN)
**Files deleted:**
- `inversion_modules/utils/sample_registry.R_v254`
- `inversion_modules/utils/sample_map.R`

Byte-identical duplicates of `utils/sample_registry.R` and
`utils/sample_map.R` respectively (md5 confirmed). The `_v254` suffix
suggests an informal archive marker. Grep across whole codebase
confirms no script sources these copies. Removed from chat-11 tarball
via `rsync --exclude`.

### FIX 57 — FEATURE — Test fixtures + cheatsheet
**Files created:**
- `registries/tests/test_interval_registry_extensions.R` (145 lines)
- `registries/tests/test_sample_registry_extensions.R` (175 lines)
- `registries/API_CHEATSHEET.md` (170 lines)

Each test fixture resolves the loader relative to its own script path
(script_dir → ../api/R/registry_loader.R, with cwd fallback), stands
up a temp registry, seeds a deterministic synthetic topology, exercises
every new method with assertions. Runs as `Rscript
registries/tests/test_*.R` from the repo root. Exits 0 on pass, stops
with `FAIL: <msg>` on first failure. Both fixtures verified passing
end-to-end against the staged tarball tree.

Cheatsheet indexes all methods with R/Py availability columns and
common usage patterns (C01f collapse, phase-4e characterization,
overlap reconciliation, QC flagging).

## Chat 11.5 — decompose redesign per chat-9 + orphan dispatch

### FIX 58 — DISPATCH — STEP_C01j regime compatibility engine (verbatim rehome)
**File moved:** `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01j_regime_compatibility_engine.R`
**Lines:** 927 (unchanged from source)
Previously orphan (missing from v10). Core of recombinant detection
per chat-9 Tier 2 spec: sliding-window sample × sample Hamming
distance → Ward hierarchical clustering with adaptive k →
compatibility groups → per-window sample memberships →
segment_regimes with 6 state labels (clean_inversion / structured_*
/ background_soup / weak_signal / transition). Parse-clean (60
top-level expressions).

Wiring to registry deferred to chat 13.

### FIX 59 — DISPATCH — STEP_C01l local structure segments (verbatim rehome)
**File moved:** `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01l_local_structure_segments.R`
**Lines:** 562 (unchanged from source)
5-segment per-candidate computation of delta12 / entropy / ENA from
local k-means membership proportions. Segments: left_flank,
inv_left_half, inv_core, inv_right_half, right_flank. Replaces the
ad-hoc staircase boundary_scan as the proper boundary-sharpness
signal. Parse-clean (53 top-level).

### FIX 60 — DISPATCH — STEP_C01m distance concordance (verbatim rehome)
**File moved:** `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01m_distance_concordance.R`
**Lines:** 567 (unchanged from source)
Multi-scale (d ∈ {80, 160, 320, 640} windows) sample-pair concordance.
Biologically: inversion-carrier pairs co-group at ALL tested distances;
family-LD pairs decay. Persistence profile distinguishes the two.
Parse-clean (21 top-level).

### FIX 61 — DISPATCH — STEP_C01k annotated simmat (verbatim rehome)
**File moved:** `inversion_modules/phase_4_postprocessing/4e_final_classification/STEP_C01k_annotated_simmat.R`
**Lines:** 317 (unchanged from source)
Per-chromosome integrative figure: raw simmat heatmap + triangle
outlines + system span bars + regime state strip (from C01j) +
tier/verdict annotations + family/inversion labels (from C01f) +
recombinant zones + bridge arcs. Parse-clean (37 top-level).

### FIX 62 — SCOPE — gene_conversion_detector.R (narrowed from recombinant)
**File:** `inversion_modules/phase_4_postprocessing/4b_group_proposal/gene_conversion_detector.R`
**Lines:** ~260 NEW
Renamed and narrowed. Originally drafted as a recombinant detector
using per-SNP CUSUM binary segmentation. Design review showed the
algorithm's sensitivity is wrong for that purpose (per-SNP scanning
flags ≤500 bp gene-conversion tracts as "changepoints" — false
positives for recombinants). Kept with narrowed scope:
- Windowed binning (40 SNPs/window, 10-SNP step by default)
- max_tract_bins gate (default 5) to EXCLUDE long tracts (those are
  recombinants, caught by C01j)
- min_tract_bins gate (default 1) plus magnitude threshold (default 0.5)
  to exclude single-window noise
- Direction classification: REF_in_INV_context, INV_in_REF_context,
  REF_in_HET_context, INV_in_HET_context
Output: per-sample short-tract catalog. Auxiliary annotation only —
does NOT gate recombinant calls per the combination rule.

### FIX 63 — PORT — lib_step03_seed_loader.R (Tier-1 enhancement)
**File:** `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_step03_seed_loader.R`
**Lines:** ~200 NEW
Closes chat-9 §Tier 1 (Finding S). Loads STEP03 deconvolution seeds
from `04_deconvolution_seeds/{inv_id}_seeds.tsv`. Combines with
flashlight seeds: agreement → union, disagreement → drop, priority
flashlight > STEP03 when forced to pick. `combine_tier1_seeds`
returns status ∈ {"ok", "flashlight_only", "step03_only",
"no_seeding"} with reason strings. NO unsupervised k-means fallback —
returns `no_seeding` status and caller skips the candidate.

### FIX 64 — PORT — lib_ghsl_confirmation.R (Tier-3 confirmation)
**File:** `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_ghsl_confirmation.R`
**Lines:** ~180 NEW
Closes chat-9 §Tier 3. Loads GHSL v5 karyotype runs from
`<chrom>.ghsl_v5.karyotypes.rds` (or equivalent). Per-sample status
in candidate interval ∈ {INV_INV, INV_nonINV, SPLIT, UNINFORMATIVE}.
SPLIT = sample has runs with BOTH call types in the interval =
recombinant confirmation. `ghsl_interval_resolution()` diagnostic
reports whether the interval has enough windows for GHSL to call runs
at all (default 10 windows — below that, Tier 3 is silent and Tier 2
stands alone at MEDIUM confidence).

### FIX 65 — FEATURE — lib_recomb_combination.R (chat-9 combination rule)
**File:** `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_recomb_combination.R`
**Lines:** ~280 NEW
Implements the chat-9 §Combination rule exactly. Consumes:
- R = recombinant signal (derived from C01j regime memberships)
- G = GHSL SPLIT (Tier 3)
- C = gene-conversion tracts (orthogonal annotation)
- B = boundary sharpness (derived from C01l segment summary; not in
  the gate path; feeds validation confidence in 4c)

Final per-sample status: RECOMBINANT/HIGH, RECOMBINANT/MEDIUM,
recomb_disputed, recomb_ghsl_only, NOT_RECOMB. Orthogonal annotation
distinguishes "RECOMB + GC tracts" from "GC-only".

KNOWN ISSUE: `derive_R_from_regime` uses naive consecutive-long-
segments-differ counting. Rigid A→B→A required to fire. User flagged
this as too rigid; correct semantics are DAG-based (see
HANDOFF_PROMPT_chat12). Fix scheduled for chat 12.

### FIX 66 — SCHEMA — 4 new Tier-2 evidence blocks
**Files NEW:** in `registries/schemas/structured_block_schemas/`:
- `regime_segments.schema.json` — C01j output
- `local_structure_segments.schema.json` — C01l output
- `distance_concordance.schema.json` — C01m output
- `gene_conversion_tracts.schema.json` — gene_conversion_detector output

Each has `keys_extracted` directives mapping sub-fields to flat
keys for the evidence registry (e.g. q2_regime_dominant_state,
q2_local_struct_core_delta12, q2_distance_conc_inv_vs_fam_score,
q2_gc_total_tracts).

### FIX 67 — FEATURE — 4 new query API methods on registry loader
**File:** `registries/api/R/registry_loader.R`
**Lines:** ~90 NEW in `load_query_api()`
- `reg$query$regime_segments(cid)` — returns segments + transitions +
  dominant_state + recomb flag
- `reg$query$local_structure_segments(cid)` — returns segment_summary
  + by_inversion_class + boundary_sharpness
- `reg$query$distance_concordance(cid)` — returns per_distance_summary
  + persistent_carriers + inversion_vs_family_score
- `reg$query$gene_conversion_tracts(cid)` — returns per_sample_summary
  + tracts + totals

Registry loader now 1328 lines (was 1166 after chat 11).
