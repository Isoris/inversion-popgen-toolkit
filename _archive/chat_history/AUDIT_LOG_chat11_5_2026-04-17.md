# AUDIT LOG — chat 11.5 (2026-04-17) — decompose redesign + orphan dispatch

## Scope

Bridge session between chat 11 (registry library buildout) and chat 12
(coherence audit + wiring). Covers:
- Recognition that the first draft of decompose-v2 had drifted from
  chat-9 three-tier design spec.
- Dispatch of four orphan phase-4 scripts (STEP_C01j / C01k / C01l / C01m)
  from user-provided uploads into their proper folders in the phase-4
  tree. These scripts were written BEFORE the v10 rewrite and forgotten;
  they are the correct implementations of the compatibility-based
  recombinant detector, per-segment ENA boundary sharpness,
  multi-scale sample concordance, and per-chromosome integrative figure.
- Reframing of the CUSUM changepoint detector as a gene-conversion
  detector (narrowed scope; not a recombinant detector).
- Port of STEP03 seed loader + GHSL confirmation libraries from
  chat-11.5 draft (the only pieces that were correctly aligned with
  chat-9 design).
- Combination rule implementation per chat-9 §Combination.
- Four new evidence block schemas + four new query API methods.
- Identification of the DAG-based recombinant semantics rewrite
  scheduled for chat 12.

## Reading done this chat

- `DECOMPOSE_UPGRADE_DESIGN_chat9.md` (421 lines) — the authoritative
  three-tier design specification.
- `STEP_C01j_regime_compatibility_engine.R` (927 lines) — full read,
  focus on `compute_compatibility_groups` (L136), `track_regimes`
  (L224), `segment_regimes` (L315), `analyze_membership_stability`
  (L424).
- `STEP_C01l_local_structure_segments.R` (562 lines) — full read,
  focus on `compute_local_memberships`, `compute_metrics`,
  `define_segments` (5-segment boundary definition).
- `STEP_C01m_distance_concordance.R` (567 lines) — scanning read,
  focus on `compute_concordance_at_distance` and
  `compute_concordance_fast`.
- `STEP_C01k_annotated_simmat.R` (317 lines) — scanning read of the
  figure composition layers.
- Existing phase-4 files for context: `STEP_C01i_decompose.R`,
  `STEP_C01i_b_multi_recomb.R`, `STEP_C01e_candidate_figures.R`,
  `STEP_C01i_c_nested_composition.py`,
  `docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md`.
- Chat-9's Snake 3 v5 to confirm GHSL karyotype output structure.
- Existing `utils/sample_registry.R` to ensure no new conflicts with
  Finding AL (load_registry shadowing).

Total ≈ 3500 lines read + cross-referenced across code.

## Critical design corrections (vs. first draft)

### Before (chat-11.5 first draft, WRONG)

- Proposed `classify_via_dosage` used `mean(dosage)` inside the
  interval as primary classifier. This is not a changepoint test;
  it's an interval-aggregate genotype call that doesn't detect
  recombinants.
- Kept unsupervised k-means fallback in `classify_via_pc1_kmeans`.
  Chat-9 design explicitly forbids unsupervised fallback.
- Did not wire STEP03 seeds (Finding S).
- Consumed GHSL in decompose instead of multi_recomb (wrong stage).
- Invented "S4 soft-boundary signal" instead of using the existing
  per-segment ENA from the C01l orphan.
- Per-SNP CUSUM binary segmentation as the recombinant detector.
  This algorithm's resolution is too fine — it flags 50–500bp gene
  conversion tracts as changepoints. User feedback: "what you have
  coded it's for gene conversion not for recombinants."

### After (chat-11.5 correct state)

Three-tier detector per chat-9 design, with proper separation of
concerns:

| Tier | Purpose | Implementation |
|---|---|---|
| 1 | Baseline HOM_REF/HET/HOM_INV classification | PC1 k-means with flashlight ∪ STEP03 seeds; no unsupervised fallback |
| 2 | Recombinant detection | C01j compatibility groups per sliding window (Hamming distance → Ward hierarchical clustering → adaptive k) |
| 3 | Independent confirmation | GHSL v5 karyotype runs (SPLIT = within-interval call change) |
| 4 | Diagnostic only | Existing PCA per-window class track, demoted |
| Orthogonal | Gene conversion | Windowed-binned excursion detector (`gene_conversion_detector.R`) — does NOT gate recomb |
| Boundary | Per-segment ENA/Δ12 | C01l — 5-segment characterization |
| Family-vs-inv | Sample concordance across distance | C01m |
| Reporting | Per-chrom figure | C01k |

## Dispatch map (orphan scripts to proper homes)

| Source | New home | Rationale |
|---|---|---|
| `STEP_C01j_regime_compatibility_engine.R` | `4a_existence_layers/` | Produces existence-layer evidence (regime state labels along candidate) |
| `STEP_C01l_local_structure_segments.R` | `4a_existence_layers/` | Per-candidate boundary/interior structure metrics; sits alongside C01g boundary catalog |
| `STEP_C01m_distance_concordance.R` | `4a_existence_layers/` | Per-chromosome existence-layer signal (family-vs-inversion distinguisher) |
| `STEP_C01k_annotated_simmat.R` | `4e_final_classification/` | Final integrative figure per chromosome, phase 4e reporting |

Filed verbatim (no logic changes). Tuning of parameters and wiring to
registry deferred to chat 12 (audit) and chat 13 (wiring).

## Findings surfaced this session

### Finding AO (chat 11.5) — per-SNP CUSUM is the wrong resolution for recombinants

Per-SNP changepoint scanning flags gene-conversion tracts (≤500 bp) as
double crossovers. This caused the initial draft to miscategorize GC
events as recombinants. **Resolved** by renaming the detector to
`gene_conversion_detector.R` with narrowed scope (max_tract_bins = 5,
windowed binning) and kept ONLY for GC detection.

### Finding AP (chat 11.5) — "mean dosage inside interval" is not changepoint detection

The `classify_via_dosage` function in the first draft used
`mean(dosage)` as the classifier. This is interval-aggregate, not
position-aware, and cannot detect recombinants. **Resolved** by
deleting `classify_via_dosage` entirely and using C01j's compatibility
groups as the primary recombinant signal.

### Finding AQ (chat 11.5) — unsupervised k-means fallback silently emits noisy groups

The first draft preserved the v9 unsupervised k-means fallback. At 9x
coverage with imbalanced class sizes this is genuinely unreliable.
Chat-9 design explicitly says to remove it. **Resolved** in
`lib_step03_seed_loader.R::combine_tier1_seeds` which returns
`status = "no_seeding"` when no seeds are available. Decompose caller
must check this status and skip the candidate.

### Finding AR (chat 11.5) — C01j's structure_score thresholds may not suit 9x cohort

Default thresholds in `segment_regimes`: `structure_score > 0.6` for
structured, `> 0.35` for weak_signal. These were tuned for the
original dataset (likely higher coverage). For 9x catfish cohort with
pooled-family LD, thresholds may need lowering. **Deferred** to chat 12
audit.

### Finding AS (chat 11.5) — C01l's flank_bp default 500kb may be too large

The 5-segment definition uses `flank_bp = 500000` as the flanking
region size around the candidate. For Quentin's ~Mb-scale inversions
this is ~50% of the interval on each side, which means the "flank"
statistics are computed on regions as large as the inversion interior.
Likely fine but worth verifying against specific cases. **Deferred** to
chat 12 audit.

### Finding AT (chat 11.5) — derive_R_from_regime logic is too rigid

Initial implementation counts consecutive-long-segments with different
values. User feedback: correct semantics are DAG-based. A sample's
regime track should be summarized as a mini DAG (nodes = distinct
regimes visited, edges = transitions), with soft-threshold
deviation_fraction as the R_fired gate, not rigid A→B→A match.
Cohort-level breadth (n_samples_with_deviation) feeds C01d scoring
separately. **Primary code change scheduled for chat 12.**

### Finding AU (chat 11.5) — C01e expects C01j output at standalone path

`STEP_C01e_candidate_figures.R` Panel H reads
`regime_windows.tsv.gz` from a standalone regime_dir argument. After
C01j moves into the registry-wired world, C01e needs to read from
`reg$query$regime_segments(cid)` instead. **Deferred** to chat 13 wiring.

### Finding AV (chat 11.5) — nested_composition and C01l both compute ENA

`STEP_C01i_c_nested_composition.py` already reads `delta12`, `entropy`,
`ena` from the local_Q artifact A per chat-4's design note, averaging
them across the interval. C01l computes per-segment ENA on PC1-derived
membership (different input). These are complementary but could cause
confusion. **Deferred** to chat 12 audit for a design decision: do
they merge into one module, or stay separate with clearly distinct
roles?

## Files created or modified this session

### Dispatched (verbatim, moved only)

- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01j_regime_compatibility_engine.R` (927 lines)
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01l_local_structure_segments.R` (562 lines)
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01m_distance_concordance.R` (567 lines)
- `inversion_modules/phase_4_postprocessing/4e_final_classification/STEP_C01k_annotated_simmat.R` (317 lines)

### Created

- `inversion_modules/phase_4_postprocessing/4b_group_proposal/gene_conversion_detector.R` (~260 lines)
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_step03_seed_loader.R` (~200 lines)
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_ghsl_confirmation.R` (~180 lines)
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/lib_recomb_combination.R` (~280 lines) — has known issue (Finding AT, DAG rewrite pending)
- `registries/schemas/structured_block_schemas/regime_segments.schema.json`
- `registries/schemas/structured_block_schemas/local_structure_segments.schema.json`
- `registries/schemas/structured_block_schemas/distance_concordance.schema.json`
- `registries/schemas/structured_block_schemas/gene_conversion_tracts.schema.json`
- `HANDOFF_PROMPT_chat12_2026-04-17.md` (phase-4 audit + DAG rewrite directive)

### Modified

- `registries/api/R/registry_loader.R` — added 4 query methods
  (regime_segments, local_structure_segments, distance_concordance,
  gene_conversion_tracts). Now 1328 lines (was 1166 after chat 11).

## Smoke tests

- `detect_cohort_changepoints` (repurposed as gene-conversion detector
  with windowing): 10-sample synthetic cohort, 1 clean single-crossover,
  1 double-crossover, 1 short GC tract. Detector correctly catalogs
  the short GC tract in CGA010, correctly does NOT fire for HOM
  samples. (Note: with the new framing, the CUSUM flagging of the
  SCO and DCO samples is not "correct" or "incorrect" — those are no
  longer in scope for this module. Recombinants are detected by
  C01j.)
- `classify_via_ghsl` / `ghsl_per_sample_in_interval`: synthetic GHSL
  karyo_dt with one SPLIT sample. Correctly classified.
- `combine_tier1_seeds`: agreement + conflict + sparse-source paths
  all tested.
- `derive_R_from_regime` with naive counter: found to incorrectly
  fire on short B burst (Finding AT). Counting logic patched to
  check "consecutive long segments with different values" but STILL
  wrong per DAG-based semantics; full fix deferred to chat 12.

## Parse-check backlog

All new files parse-clean:

- `STEP_C01j_regime_compatibility_engine.R` — 60 top-level
- `STEP_C01l_local_structure_segments.R` — 53 top-level
- `STEP_C01m_distance_concordance.R` — 21 top-level
- `STEP_C01k_annotated_simmat.R` — 37 top-level
- `gene_conversion_detector.R` — 8 top-level
- `lib_step03_seed_loader.R` — 4 top-level
- `lib_ghsl_confirmation.R` — 7 top-level
- `lib_recomb_combination.R` — 6 top-level
- `registry_loader.R` — 9 top-level (1328 lines)

## Cumulative findings carry-forward

Open from chat 11:
- AL (load_registry shadowing) — chat 12 to fix
- T (README structure_type drift) — doc pass
- V (cheat24 inline fallback threshold mismatch) — chat 12 to delete
  inline fallback per chat-9 §Deleted
- W (register_C01f_keys verification) — chat 15
- Z (4d README stale) — doc pass
- 8 (4a writes no registry blocks) — chat 13 wiring resolves
- 9 (C01g silent skip) — chat 15

New from chat 11.5:
- AO: per-SNP CUSUM wrong resolution for recombinants — resolved
- AP: mean-dosage is not changepoint detection — resolved
- AQ: unsupervised fallback — resolved (seed-loader returns
  no_seeding status)
- AR: C01j structure_score thresholds may need retuning for 9x cohort
  — deferred to chat 12
- AS: C01l flank_bp 500kb may be too large for ~Mb intervals —
  deferred to chat 12
- AT: derive_R_from_regime needs DAG-based rewrite — primary chat-12
  code change
- AU: C01e Panel H wiring to registry — deferred to chat 13
- AV: nested_composition and C01l both compute ENA — deferred to
  chat 12 design decision

## End-of-chat-11.5 state

Library: all 4 orphan scripts dispatched. 4 new schemas, 4 new
query methods, 4 new libraries (gc detector, step03 loader, ghsl
confirmation, combination rule). Decompose and multi_recomb NOT YET
modified to use the new libraries — that's chat 12's work after the
audit.

Chat 12 is the coherence audit. Chat 13 is registry wiring. Chat 14
is first HPC run. No mechanism changes after chat 13 until the HPC
run produces data to iterate on.
