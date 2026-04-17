# AUDIT_PHASE4_COHERENCE_2026-04-18.md

Chat 12 coherence audit of `inversion_modules/phase_4_postprocessing/`.

**Posture:** read before code, audit before wire. This is the pause before
chat 13 performs the registry wiring. Findings continue the A–AN series
from chat 11; new findings here are AO+.

**Corpus audited:** ~15,500 lines across 22 scripts in 4a/4b/4c/4e, plus
three chat-11.5 libraries (`lib_recomb_combination.R`,
`lib_step03_seed_loader.R`, `lib_ghsl_confirmation.R`) and the
`gene_conversion_detector.R` module. 4d (cheats) is a dispatch-level spot
check only — its per-cheat logic is out of scope for a coherence audit.

**Parallel code changes landed in chat 12** (alongside this document):

1. `lib_recomb_combination.R::derive_R_from_regime` rewritten with the
   DAG formulation (per_sample DAG list + cohort aggregates). Back-compat
   via `as_regime_dt()` unwrap and defensive unwrap in
   `combine_cohort_recomb`. Smoke-tested: 50/50 passing in
   `tests/test_dag_derive_R.R` including the AAABBA and AAABABBA patterns
   the handoff called out.
2. `gene_conversion_detector.R` rewritten as a per-SNP run-length
   detector. QC filter → diagnostic mask → per-sample flag → run-length
   with tolerance → length gate → geometric-prior confidence. Schema
   updated (`snp_qc` block, `run_length_flagged_snps`, `tolerance_snps`,
   `confidence`, `span_bp`). Smoke-tested: 33/33 passing in
   `tests/test_gc_detector.R`.
3. `regime_sample_dag.schema.json` (new). Carries the per-sample DAG
   block that feeds C01d Layer C in chat-13 wiring.
4. `STEP_C01i_decompose.R`: `--step03_seeds_dir` added, seed path routed
   through `combine_tier1_seeds()`, unsupervised k-means fallback
   **removed**, `decomp_status="no_seeding"` emitted on sparse-seed
   skips.
5. `STEP_C01i_b_multi_recomb.R`: full rewrite. Consumes
   `lib_recomb_combination.R`. R-and-G gate. Inline cheat24 fallback
   **deleted** (Finding V). 100kb mosaic threshold **removed**. Old
   Signal 1 demoted to Tier-4 diagnostic. Writes `regime_sample_dag`
   block alongside `recombinant_map`. Event-class split
   (RECOMBINANT / RECOMBINANT_GC / RECOMBINANT_DCO).
6. `STEP_C01f_hypothesis_tests.R::comp_from_registry` collapsed from
   4×`has_group`+4×`get_group` to one `get_groups_for_candidate(cid)`
   call. Legacy path preserved for pre-chat-11 registries.
7. `plot_sample_regime_dag.R` (new). Base-R per-sample DAG diagnostic.

What this document does NOT do: wire the new evidence blocks into C01d
Layer C scoring, seal registration, or phase-4e characterizer. That is
chat 13's job, as specified by the handoff DoD.

---

## 1. Per-script audit

### 1.1 `STEP_C01d_candidate_scoring_wired_25_v934_registry.R`  (1050 L)

**Purpose.** Layer-A existence scoring. Computes 12 dimensions (D1–D12),
composite `final_score`, tier (1–4) by threshold count, and a pattern
classification enum. Writes the candidate scoring table that everything
else in phase-4 reads.

**Function inventory** (key functions only):

- `safe_num(x, default=0)` — defensive numeric coercion.
- Main loop over `iv` rows computes D1..D10 inline from detector
  scoring_table columns.
- Post-loop passes fill D11 (boundary concordance from C01g catalog)
  and D12 (snake-staircase concordance from C01b core_regions).
- Pattern assignment at L382–405.

**What each D actually measures.**

| Dim | Source | Measures |
|-----|--------|----------|
| D1  | scoring_table `block_score` | Block contrast + shape strength |
| D2  | scoring_table `shape_class_score` | Shape-class match |
| D3  | scoring_table `nn_persist`  | Nearest-neighbour persistence |
| D4  | scoring_table `flat_decay`  | Inversion-like flat LD decay |
| D5  | scoring_table `interior`    | Interior homogeneity |
| D6  | scoring_table `multi_mat`   | Found in ≥2 matrix variants |
| D7  | flashlight `sv_concord`     | SV breakpoint concordance |
| D8  | peel diagnostic             | Peel stability / hypothesis |
| D9  | precomp PCA quality         | PCA clustering quality |
| D10 | bootstrap partition         | Partition stability |
| D11 | C01g boundary catalog       | Boundary concordance (post-loop) |
| D12 | C01b core_regions           | Snake-staircase overlap (post-loop) |

**Final-score formula (L354–357).**

    final_score = 0.14·d1 + 0.08·d2 + 0.13·d3 + 0.06·d4 + 0.05·d5
                + 0.06·d6 + 0.06·d7 + 0.12·d8 + 0.06·d9 + 0.06·d10
                + 0.09·d11 + 0.09·d12

Weights sum to 1.00 (checked). Tier assignment is count-based, not
score-based: `dim_positive >= {8,6,4}` maps to tiers {1,2,3}, else 4.
The count-and-score are not redundant — the score is continuous for
ranking, the count is discrete for promotion. OK.

**Drift from documentation.**

- Header claims "10 dimensions" in one place (v8.4 reference) and "12
  dimensions" elsewhere (v9.3.2+). Reality is 12. *Finding AO.*
- `n_children` at L383 is read from `iv$n_children %||% iv$n_nested`,
  but nothing upstream in the current phase-3 pipeline writes
  `n_children`. The fallback to `n_nested` works for the current
  scoring table, but the primary field is aspirational. Not a bug;
  mis-document. *Finding AO (continued).*

**Suspect thresholds.**

- Dim-positive cutoffs at L361–372 differ per dimension (0.20..0.60).
  D8 at 0.60 is tight; D7 at 0.20 is loose. D8 being strict is
  defensible because peel disappearance is the high-confidence confound
  flag. D7 at 0.20 accepts any flashlight hint; at 9x with small SV
  panels that's arguably correct. No change recommended, but worth
  recording the rationale.
- Downgrade rule at L380 is double-required: both peel levels must say
  disappeared AND tier currently ≤ 2. Tight. Leave.

**Pattern classification (L387–405).** The enum currently supports:
`complex_system`, `nested_fixed`, `nested_rare`, `family_ld`,
`strong_inversion`, `diffuse_inversion`, `diagonal_band`, `noise`,
`unclassified`. The Inv14.2-style "partial-overlap-neighbour" case
has no entry. Chat 13 should add `partial_overlap_neighbor` driven by
the new evidence: high `dominant_regime_support` with a non-trivial
`fraction_samples_deviating` AND a specific regime-transition pattern
in `regime_segments`. *Finding AP: enum missing a partial-overlap
bucket — will be reachable once C01j evidence wires into C01d.*

**Dead vs live outputs.**

- The table is written to the scoring path and consumed by C01e
  figures, C01f hypothesis tests, C01g boundary catalog, seal, and
  characterize. All downstream readers have been verified to find it
  (chat 10). No dead columns.
- What is MISSING from C01d's output that chat-13 wiring will add:
  `regime_dominance`, `boundary_sharpness`, `inversion_vs_family_score`
  from the new evidence blocks. These go into `existence_layer_c`,
  not `existence_layer_a`, so C01d itself stays the same — but the
  Layer-A scoring could optionally ingest them as D13–D15. Flag for
  chat 15 (not urgent).

**Silent-failure audit.** No obvious `tryCatch(..., error = function(e) NULL)`
swallowing real errors in the scoring loop. The post-loop D11/D12
fills do use `tryCatch`s to handle missing files, which is correct
behaviour (defaults to 0 with status="no_data").

**Verdict.** Stable. Wire-ready once C01j/C01l/C01m outputs feed Layer C.
No change in chat 12.

---

### 1.2 `STEP_C01e_candidate_figures.R`  (694 L)

**Purpose.** Per-candidate diagnostic figure composer (the 9-ish panel
plot). Panel H reads C01j regime output.

**Wiring check.** After the chat-11.5 move of `STEP_C01j_*.R` into
`4a_existence_layers/`, C01e's hardcoded paths need verification. I
did not fully trace Panel H in this audit; *Finding AQ — low
confidence that Panel H's C01j path resolution still works after
the directory move. Chat 13 should grep for hardcoded
`STEP_C01j` or `regime_memberships` path patterns and compare to
the new location.*

**Verdict.** Not touched in chat 12. Flag for chat 13.

---

### 1.3 `STEP_C01g_boundary_catalog_*.R`  (1446 L)

**Purpose.** Boundary catalog. Writes `boundary_left` / `boundary_right`
blocks via the registry. Currently single-point estimates.

**Audit question (handoff): should C01g consume C01l?** Answer: **YES,
eventually.** C01l produces per-segment ENA/Δ12 that let you compute
a "boundary fuzz width" — the span over which the segment metrics
transition from core to flank levels. The current C01g catalog has a
single `boundary_bp` point per side; wrapping a `fuzz_width_bp` field
around it would genuinely improve downstream D11 scoring and let the
manuscript's boundary-table be more informative.

However — this is a chat-13/14 wiring task, not a chat-12 audit fix.
The boundary schema is frozen and would need extension. The handoff
explicitly says "Only ADD new schemas … Never modify existing" — so
adding `fuzz_width_bp` requires a new `boundary_v2` schema with a
migration path for existing boundary blocks. *Finding AR: C01g +
C01l integration is real payoff but non-trivial; defer to chat 14
after first HPC run gives us real examples to calibrate against.*

**Suspect thresholds.** Not audited in depth — the boundary catalog
logic is from the pre-chat-9 era and worked in the v9.3.4 run. Leave.

**Verdict.** No change in chat 12. Plan the v2 schema for chat 14.

---

### 1.4 `STEP_C01j_regime_compatibility_engine.R`  (927 L)  **NEW in 4a**

**Purpose.** Sliding-window compatibility clustering over per-window
PC1 vectors. Emits `regime_memberships` (per-window per-sample group
ID) and derives regime segments (runs of consecutive windows with
similar structural properties).

**Function inventory.**

- `compute_compatibility_groups` (L136): per-window clustering of
  samples into compatibility groups from PC1 distances.
- `track_regimes` (L224): tracks group membership across windows,
  handling group-relabelling between windows.
- `segment_regimes` (L315): collapses consecutive similar windows into
  regime segments. Emits segment_id, state, start/end, mean metrics.
- `analyze_membership_stability` (L424): quantifies how stable a
  sample's group membership is across windows.

**Threshold audit.**

- `structure_score > 0.6` (L335, L338, L341) separates simple /
  moderate / complex flavours of clean-inversion call.
- `structure_score > 0.35` (L344) is the weak-signal cutoff; below
  0.35 is "background_soup" for thick windows, "weak_signal" otherwise.
- The 0.35 / 0.6 cuts are defensible for LD-style structure scores
  when signal is plentiful. At 9x with limited samples, I'm
  uncertain: the expected structure_score distribution over Quentin's
  cohort has not been characterised in the audit log. **Recommendation
  for chat 14:** run on LG12 (known strong inversion) and LG25 (the
  paralog-heavy control), plot the distribution of per-window
  structure_score, and tune if the 0.35/0.6 bands don't separate cleanly.
  *Finding AS.*

**Output routing.** Currently writes standalone TSVs to `--outdir`.
Chat-13 wiring task: redirect via `write_block("regime_segments", ...)`
for the segments and add a new evidence path for the raw
`regime_memberships` (which is too large to store in the registry as
JSON; it should probably stay as a sidecar TSV.gz with the registry
block pointing to it). **Don't rewrite the compute logic.** The logic
is trusted per the handoff.

**Verdict.** No compute-logic change in chat 12. Wiring task for chat 13.

---

### 1.5 `STEP_C01l_local_structure_segments.R`  (562 L)  **NEW in 4a**

**Purpose.** Five-segment local-structure characterization per
candidate: `left_flank / inv_left_half / inv_core / inv_right_half /
right_flank`. Computes per-segment mean delta12, entropy, ENA from
local k-means membership proportions.

**Functions.** `compute_local_memberships`, `compute_metrics`,
`define_segments`.

**Threshold audit.**

- `FLANK_BP = 500000` (L29, opt L65). For Quentin's ~Mb-scale
  inversions this is 50% of the average span — flanks are almost as
  big as the inversions themselves. This is fine for Quentin's cohort
  (the large inversions have plenty of outside-context window), but
  for smaller candidates (~300 kb span) it means flanks extend 1.67×
  the inversion length, which risks including unrelated structure.
  *Finding AT: `flank_bp` should scale with candidate span, e.g.
  `flank_bp = clamp(span_bp, 200kb, 500kb)`. Trivial chat-13 change.*

**Output routing.** Same as C01j — standalone TSV now; wire via
`write_block("local_structure_segments", ...)` in chat 13.

**Verdict.** No compute change. Threshold tweak (AT) is small and
optional; defer to chat 13.

---

### 1.6 `STEP_C01m_distance_concordance.R`  (567 L)  **NEW in 4a**

**Purpose.** Multi-scale pair concordance — for each pair of samples,
compute concordance at distances (80, 160, 320, 640) windows apart.
Used to disambiguate inversion signal (concordance decays with
distance within inversion carriers) from family structure
(concordance stays flat across distances).

**Threshold audit.**

- `DISTANCES = c(80, 160, 320, 640)` (L58). With 50-marker windows at
  Quentin's 1.6 kb/SNP density, 80 windows ≈ 128 kb, 640 windows ≈
  1.02 Mb. This covers the 100 kb – 1 Mb range, which matches the
  expected inversion-vs-family signal separation well. Defensible.

**Output routing.** Same — wire via
`write_block("distance_concordance", ...)` in chat 13.

**Verdict.** No change.

---

### 1.7 `STEP_C01i_decompose.R`  (452 L → 531 L after chat-12 patch)

**Purpose.** Per-candidate Tier-1 PC1 k-means decomposition into
HOM_REF / HET / HOM_INV. Writes `internal_dynamics` block.

**Changes landed in chat 12.**

1. Sources `lib_step03_seed_loader.R` alongside helpers.
2. Added CLI: `--step03_seeds_dir`, `--min_seeds_per_class`.
3. Routes seeds through `combine_tier1_seeds()` (flashlight ∪ STEP03,
   conflict drop).
4. **Removed unsupervised k-means fallback.** The chat-9 design was
   explicit: if no seeds, don't ship. Before this patch, a candidate
   with no SV anchors fell through to `kmeans_with_quality(mean_pc1,
   k=3, nstart=25)`, which could emit essentially random 3-way splits
   on pure noise intervals — and those splits got registered as
   groups. That was always going to pollute the downstream tests.
5. Emits minimal `internal_dynamics` block with
   `decomp_status="no_seeding"` + reason + seed_source + conflict
   count on skip. Downstream readers can key off `decomp_status` to
   distinguish "processed, produced result" from "skipped, no groups
   proposed".
6. Also catches the case where `seeded_kmeans()` internally fell back
   to unsupervised (`km$seeded = FALSE`): emit `no_seeding` there too.
   Without this second check, the internal fallback would silently
   leak past our removal of the outer fallback. *Finding AU: closed
   by chat 12.*

**Remaining suspect thresholds.**

- `silhouette_score` is written to the block but never used as a gate.
  Chat 9 planned a silhouette quality flag at 0.25. The audit log
  mentions 0.25 as too lenient at 9x — **recommendation:** set the
  flag (`decomp_quality = silhouette >= 0.40 ? "clean" : "noisy"`)
  as a non-gating annotation so downstream can filter. *Finding AV:
  silhouette cutoff 0.25 → 0.40 for 9x coverage. Small chat-13 item.*
- `bic_gap` cutoff 0.05 is written but not used as a gate. Same story.

**Verdict.** Patched.

---

### 1.8 `STEP_C01i_b_multi_recomb.R`  (436 L → 300 L after full rewrite)

**Purpose.** Per-candidate recombinant detection and classification.

**Changes landed in chat 12 — the full rewrite.**

1. Drops the three-signal combine (S1 ∧ S2) ∨ S3. The old S1 (per-
   window PCA class switch) was the noisy one; C01j regime
   memberships are the correct substitute.
2. Uses `derive_R_from_regime()` (DAG formulation) for the recombinant
   gate. R-and-G is now the primary decision.
3. **100 kb mosaic-length threshold removed.** Was a workaround for
   the old S1's false-positive rate at short mosaics. No longer needed
   because the DAG's `min_deviation_bp` gate serves the same purpose
   with a cleaner semantics.
4. **Inline cheat24 fallback deleted (Finding V closed).** Under
   investigation for months: the inline version disagreed with real
   cheat24 by up to 2× on the gene_conversion posterior for the same
   input. Best policy: if real cheat24 is absent, posterior = NA.
   Better no number than a wrong number.
5. Old Signal 1 demoted to Tier-4 diagnostic — recorded per sample for
   auditability (`s1_mosaic_start_bp` etc.) but not in the gate.
6. Writes the new `regime_sample_dag` block.
7. Event-class split: samples with `recomb_status="RECOMBINANT"` get
   further subdivided:
   - `RECOMBINANT_DCO` if `longest_deviation_bp >= min_dco_bp` (def 200kb)
   - `RECOMBINANT_GC` if GC tracts > 0
   - `RECOMBINANT` otherwise
   This drives the subgroup names seal will register in chat 13.

**Suspect thresholds in the rewrite.**

- `min_dev_frac = 0.15` (default) — 15% of interval deviating from
  dominant regime. Verified by the AAABABBA smoke test that the bp
  gate is the binding constraint for alternating-noise protection;
  the fraction gate is more lenient. Defensible.
- `min_dev_bp = 50000` (default) — 50 kb longest deviation. At 10 kb
  windows (typical C01j window span), this is 5 windows. The handoff
  said "a sample with 15%+ of windows in a non-dominant regime AND
  deviations long enough to not be pure noise". 5 windows ≈ "not pure
  noise" at Quentin's density.
- `min_dco_bp = 200000` (default) — 200 kb longest deviation = DCO.
  No principled justification for 200 kb specifically; chosen as
  "clearly longer than a GC event and longer than typical noise".
  *Finding AW: min_dco_bp = 200 kb is an educated guess; recalibrate
  after chat-14 HPC run on LG12.*

**Verdict.** Patched.

---

### 1.9 `STEP_C01i_c_nested_composition.py`

**Purpose.** Nested-composite flag using delta12/entropy/ena from
local_Q. Python, not R.

**Audit question (handoff): does this conflict with C01l?** Answer:
**they complement, they don't conflict.**

- Nested_composition operates on an **interval-averaged** ENA from
  local_Q. It emits a boolean `is_nested_composite`.
- C01l operates on **per-segment** ENA across 5 segments. It emits
  segment-level numbers and the boundary_sharpness classification.

The chat-13 design should **keep both**: C01l feeds the sharpness /
dominance signal for Layer C scoring; nested_composition's boolean
feeds the pattern enum choice. They answer different questions.

*Finding AX: no conflict — wire both independently.*

**Verdict.** No change.

---

### 1.10 `STEP_C01i_d_seal.R`  (422 L)

**Purpose.** Registers per-candidate groups to the sample registry
after decompose + multi_recomb produce their outputs.

**Function inventory.** `register_all_groups` (L221) is the main
entry. HOM_STD alias logic at around L270 (not re-read in this audit).

**Chat-13 wiring task.** Seal needs to read the new
`recomb_subgroup` field (RECOMBINANT / RECOMBINANT_GC /
RECOMBINANT_DCO) from the updated multi_recomb output and register
subgroup names accordingly. The `inv_<cid>_RECOMBINANT_GC` and
`inv_<cid>_RECOMBINANT_DCO` group IDs are already anticipated by the
chat-11 `get_groups_for_candidate` API (which returns these slots).

**Verdict.** No change in chat 12 (handoff explicitly told us not
to touch seal). Wiring in chat 13.

---

### 1.11 `gene_conversion_detector.R`  (261 L → 313 L after chat-12 rewrite)

**Full rewrite.** Per-SNP run-length detector replacing the windowed
binning. Previously, the 40-SNP bin floor was ~60 kb, 30× bigger
than the biology it was trying to detect. **This was the single
biggest correctness issue in the chat-11.5 delivery.**

Details in the per-file header and in `test_gc_detector.R`. 33 tests
passing.

**Why the old windowed approach was wrong.** It treated GC as a
low-frequency signal that could be smoothed. GC is actually a rare-
event, high-spatial-resolution signal — the right representation is
"one sample, a couple of SNPs, a few kb". The windowed approach's
floor scales with how many SNPs you average over; the right floor
scales with how many SNPs you need to reject genotyping noise, which
is 2–4 under a geometric-error prior.

**Suspect thresholds in the new detector.**

- `min_delta_af = 0.5` at step 2. A SNP is diagnostic if REF cohort
  and INV cohort have AF differing by ≥ 0.5. This is tight — at 5
  HOM_REF + 5 HOM_INV samples (smoke-test shape), sampling noise in AF
  is ~0.22 per group, so 0.5 is well above noise. For larger cohorts
  the threshold could go up to 0.7. Leave default at 0.5 for now.
- `max_het_fraction = 0.7` — flag SNPs where >70% samples look het.
  At 9x per-SNP genotyping, the HWE expectation at MAF=0.5 is het
  fraction 0.5; 0.7 is 2 SD above. Reasonable for paralog detection.
- `max_span_bp = 20000` — upper bound on GC span. Handoff specified
  5 kb typical, up to 5 kb outlier, so 20 kb is a 4× safety factor.
  Anything longer is in the gap between GC and recombinant (covered
  by C01j). Leave at 20 kb.
- `max_flagged_snps = 10` — upper bound on flag count. 10 flagged
  SNPs at 1.6 kb/SNP is ~16 kb, consistent with max_span_bp.

**Verdict.** Rewritten, tested, ready. Finding Z' closed (handoff
Finding Z was about the GC detector scope).

---

### 1.12 `lib_step03_seed_loader.R`  (243 L)

**Purpose.** Tier-1 seeding from STEP03 (Fisher/Armitage SV-based
assignments) combined with flashlight anchors. Conflict drop
rather than arbitrary pick. No unsupervised fallback.

**Function inventory.**

- `load_step03_seeds` (L59): reads `<dir>/<cid>_seeds.tsv`. Applies
  `p_value <= 0.05` and `n_supporting_snps >= 5` filters. Normalises
  class labels.
- `combine_tier1_seeds` (L124): union-agreed + drop-conflicts over
  flashlight and STEP03 seeds. Returns status ∈ {ok, flashlight_only,
  step03_only, no_seeding}. Sparse-seed cases collapse to no_seeding.

**Contract check vs chat-9 design.**

- ✔ "priority flashlight > STEP03 if they conflict." The code drops
  conflicting samples from both sources — it does NOT prefer
  flashlight's assignment. **This is actually MORE conservative than
  the chat-9 language; drop is safer than pick.** Arguably the
  documentation should be updated to reflect the stricter-than-specified
  behaviour. *Finding AY: documentation drift — drop-conflict is
  stricter than priority-flashlight. Update header comment; no code
  change.*
- ✔ "remove the unsupervised fallback." Yes — `no_seeding` is the
  explicit reject state.

**Suspect thresholds.**

- `p_value_thresh = 0.05`. Standard but debatable for per-candidate
  multi-testing over 200+ candidates. At alpha=0.05 you'd expect ~10
  false-positive seeds on noise. Probably OK because a false-positive
  seed gets dropped by conflict with flashlight on the same candidate.
  *Finding AZ: consider FDR-adjusted alpha, low priority.*
- `min_supporting_snps = 5`. 5 SNPs × 1.6 kb/SNP = 8 kb of support
  per sample. For ~1 Mb inversions, 8 kb of evidence is thin; for
  300 kb candidates it's reasonable. *Finding BA: consider scaling
  with candidate span. Low priority.*

**Verdict.** Good. Minor doc drift.

---

### 1.13 `lib_ghsl_confirmation.R`  (213 L)

**Purpose.** Tier-3 within-interval karyotype confirmation via the
ghsl (group-level haplotype similarity by locus) track. Checks that
inside a candidate interval, a sample has BOTH an INV_INV run and a
non-INV_INV run → consistent with a within-interval arrangement
change = recombinant.

**Audit target (handoff):** `resolve_karyo_bp` logic. The window-
index to bp mapping relies on `annot_dt` having one row per
`global_window_id` in chrom order.

I traced this. The function calls `annot_dt[global_window_id == X]`
and expects exactly one row. If `annot_dt` is multi-row per global_id
(which shouldn't happen at steady state, but COULD after a re-run
where two annotation passes both contributed rows), the function
silently picks the first and returns the wrong bp.

*Finding BB: add a `uniqueN(annot_dt$global_window_id) == nrow(annot_dt)`
assertion at the top of `resolve_karyo_bp`. Two-line change; defer
to chat 13.*

**Verdict.** Logic correct given contract; add the assertion.

---

### 1.14 `lib_recomb_combination.R`  (303 L → 458 L after chat-12 rewrite)

Already covered above. DAG rewrite of `derive_R_from_regime`. Back-
compat shim. Smoke-tested.

**Key insight from writing the fixture:** the two-gate (fraction AND
bp) structure is what separates the recombinant signal from
alternating-noise. The AAABABBA pattern has high deviation_fraction
(0.375) — without the bp gate, it would fire R. The bp gate (based
on the longest contiguous deviation) is the real discriminator. This
is why the handoff said the old counter over-counted: AAABBA is a
single recombinant event (one long B excursion); AAABABBA is
alternating noise (many short B excursions). Both have high
deviation_fraction; only the first has high longest_deviation_bp.

---

### 1.15 `STEP_C01f_hypothesis_tests.R`  (2523 L)

**Purpose.** The heavy consumer of the sample registry. Runs T1–T10
hypothesis tests per candidate using registered groups.

**Chat-12 change.** `comp_from_registry` (L300–348) collapsed from
4×`has_group` + 4×`get_group` to one `reg$get_groups_for_candidate(cid)`.
Legacy path preserved for pre-chat-11 registries. HOM_STD alias
priority unchanged.

**What I did NOT audit in depth.** The 2300+ lines of T1..T10
implementation. These were stable through chat 10 and the registry
collapse is localised — no reason to re-audit the tests themselves.

**Verdict.** Patched. Wire-ready.

---

### 1.16 `group_validation_gate.R`  (142 L)

**Purpose.** Reads the `q6_group_validation` key; emits NONE /
UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT verdict.

**Trace end-to-end.** Registry key `q6_group_validation` is written by
C01f (somewhere around L1900+ in the big file; not re-read). The gate
reads it, maps via a lookup table to the 5-level verdict, returns.
Straightforward.

**Verdict.** Not touched; no findings.

---

### 1.17 `4d_group_dependent/` (cheats)

Spot check of launcher dispatch only. Each cheat writes its own
block type. Missing schemas: none that I spotted. The cheats are
chat-13/14 wiring concerns; not a chat-12 audit target.

**Verdict.** No findings from coherence audit.

---

### 1.18 `compute_candidate_status.R`  (924 L)

**Purpose.** Final tier + verdict. Reads from all upstream evidence.

**Audit target:** `build_key_spec` L40–413, the 367-key spec.

I did not re-parse all 367 keys. The handoff's chat-10 figure of
73.8% wired is the starting point. Chat-12 adds keys from:

- `regime_segments` schema (4 `keys_extracted`)
- `local_structure_segments` schema (~6 `keys_extracted`)
- `distance_concordance` schema (~5 `keys_extracted`)
- `gene_conversion_tracts` schema (3 `keys_extracted`)
- `regime_sample_dag` schema (5 `keys_extracted`, new this chat)

Total ~23 new keys. If the 367-key spec already has slots for most
of these (Q2 has a lot of regime/structure keys), wiring brings us
to ~80% wired. If it's all-new, ~80% also but via expansion of the
spec. Chat 13 should compute the exact number.

*Finding BC: compute the wired-fraction after chat-13 wiring. Target
is 85%+ for the HPC run.*

**Verdict.** No change; fraction to be measured.

---

### 1.19 `characterize_candidate.R`  (695 L) — the Q1–Q7 characterizers

**Audit target (handoff):** confirm Q2 (structure) pulls from
regime_segments + local_structure_segments + distance_concordance.

**Reality:** Q2 is currently stub-populated for the new blocks. Chat
13 wiring will fill in:

- From regime_segments: `q2_regime_n_segments`,
  `q2_regime_dominant_state`, `q2_regime_has_recomb`,
  `q2_regime_n_samples_changed`.
- From local_structure_segments: `q2_boundary_sharpness_*`.
- From distance_concordance: `q2_distance_pattern`.
- From regime_sample_dag (new this chat):
  `q2_dag_fraction_samples_R_fired`,
  `q2_dag_fraction_samples_deviating`,
  `q2_dag_dominant_regime_cohort`, `q2_dag_dominant_regime_support`.

**This is the biggest payoff of the chat-12 work for phase-4e.** Q2
goes from a couple of scalar fields to a rich structural description.

**Verdict.** No change in chat 12. Priority wiring task for chat 13.

---

### 1.20 `run_characterize.R`  (557 L)

**Purpose.** Phase-4e driver. Smoke-tested by chat 10.

**Verdict.** No findings; three-tier loader still works with new keys
so long as the key-spec is expanded.

---

### 1.21 `STEP_C01k_annotated_simmat.R`  (317 L) — **NEW in 4e**

**Purpose.** Per-chromosome integrative figure.

**Audit question (handoff):** does it read the actual new evidence
blocks or stub paths?

I did not trace this end-to-end in chat 12. *Finding BD: C01k needs
its output paths verified against the registry block locations
produced by chat-13 wiring. Should be straightforward once wiring is
done.*

**Verdict.** Flag for chat 13.

---

## 2. Flow-of-evidence trace for one candidate

Tracing `inv_LG12_17` (the strong-inversion exemplar) from phase-2
discovery through phase-4e:

**Phase 2 (discovery, not in scope for this audit):** inv_detect
writes `scoring_table_LG12.tsv` row for interval 17.

**Phase 4a.**

- C01d reads scoring row, computes D1..D12, writes
  `existence_layer_a` block via `write_block`.
- C01j reads precomp PC1 matrix for LG12, writes `regime_memberships`
  (sidecar tsv.gz) and `regime_segments` (block, chat-13 wiring).
- C01l reads precomp, computes 5-segment metrics, writes
  `local_structure_segments` block (chat-13 wiring).
- C01m reads precomp, computes pair concordance at 4 distances,
  writes `distance_concordance` block (chat-13 wiring).
- C01g reads C01d output + C01b core_regions, writes boundary catalog.

**Phase 4b.**

- C01i_decompose reads candidates + flashlight + (chat-12) STEP03
  seeds, writes `internal_dynamics` block with `decomp_status`.
- C01i_b_multi_recomb reads:
  - C01j regime_memberships (as input to the DAG-based R signal)
  - GHSL per-sample (for G)
  - GC detector output (for C annotation)
  writes `recombinant_map` block AND (chat-12) `regime_sample_dag`
  block.
- C01i_c_nested_composition reads local_Q, writes nested-composite flag.
- C01i_d_seal reads decompose + multi_recomb + other phase-4b blocks,
  registers `inv_LG12_17_{HOM_REF,HET,HOM_INV,RECOMBINANT,
  RECOMBINANT_GC,RECOMBINANT_DCO}` groups (chat-13 wiring).

**Phase 4c.**

- C01f reads registry via `get_groups_for_candidate("LG12_17")`
  (chat-12 change). Runs T1..T10. Writes `q6_group_validation` +
  `hypothesis_verdict` block.
- group_validation_gate reads q6 key, emits 5-level verdict.

**Phase 4d (cheats).**

- Each cheat reads what it needs (typically the registered groups),
  writes its own block type. Not audited here.

**Phase 4e.**

- compute_candidate_status reads all the above via key_spec, emits
  final tier/verdict.
- characterize_candidate reads evidence blocks, produces Q1..Q7
  answers for the candidate. Q2 gets rich fields (chat-13 wiring).
- C01k composes the per-chromosome figure from multiple candidate
  blocks.

**Aspirational writes with no reader (dead outputs):** none that I
identified after chat-12 patches.

**Aspirational reads with no writer:** after chat 12:
- `n_children` in C01d (Finding AO).
- Everything C01d / characterize_candidate / C01k would read from
  `regime_segments` / `local_structure_segments` /
  `distance_concordance` / `regime_sample_dag` — these blocks don't
  yet get written by the registry path. **This is the chat-13 wiring
  work.**

---

## 3. Pattern-enum coverage

C01d's enum today: `strong_inversion`, `diffuse_inversion`,
`diagonal_band`, `nested_fixed`, `nested_rare`, `family_ld`,
`complex_system`, `noise`, `unclassified`.

**Reached in practice** (based on chat-10 run log):
`strong_inversion`, `nested_fixed`, `complex_system`, `unclassified`,
`noise`. The others are aspirational.

**Chat-12 evidence enables new buckets.**

- `partial_overlap_neighbor` (the Inv14.2 case). Triggers: high
  regime_dominance (>0.9 samples in dominant regime) but non-trivial
  deviation_fraction in a spatially-localised segment of the
  interval. Needs an additional C01j-derived field like
  `deviation_spatial_concentration`.
- `gene_conversion_heavy` (GC cohort signal): total_samples_with_gc
  / n_samples > 0.2 = this interval has above-baseline GC activity.
- `boundary_diffuse_but_core_strong`: C01l's
  `left_sharpness`/`right_sharpness` = "diffuse" but core_delta12 high.

*Finding BE: add 3 buckets after chat-13 wiring makes the triggering
fields available. Defer to chat 14.*

---

## 4. Threshold drift inventory

Flagged thresholds across the codebase, with source and verdict:

| Threshold | Where | Value | Verdict |
|-----------|-------|-------|---------|
| `mosaic_len > 100000` | old multi_recomb gate | 100 kb | **REMOVED chat 12** |
| cheat24 inline vs real | old multi_recomb fallback | 2× divergence | **REMOVED chat 12 (Finding V)** |
| `silhouette >= 0.25` | decompose quality flag | 0.25 | Too lenient at 9x → 0.40 (Finding AV) |
| `bic_gap > 0.05` | decompose quality flag | 0.05 | Not currently used as gate; leave |
| `phase_concordance` | decompose | 0.30 | Not gating; OK |
| `structure_score` cutoffs | C01j | 0.35 / 0.6 | Recalibrate after HPC run (Finding AS) |
| `flank_bp` | C01l | 500 kb | Should scale with candidate span (AT) |
| `DISTANCES` | C01m | 80/160/320/640 | OK |
| `min_dev_frac` | new DAG gate | 0.15 | OK |
| `min_dev_bp` | new DAG gate | 50 kb | OK |
| `min_dco_bp` | new DAG DCO split | 200 kb | Recalibrate (AW) |
| `min_delta_af` | GC detector | 0.5 | OK; could tighten with larger cohort |
| `max_het_fraction` | GC detector | 0.70 | OK |
| `max_span_bp` | GC detector | 20 kb | OK (4× safety over biology) |
| `p_value_thresh` | STEP03 seeds | 0.05 | Consider FDR (AZ, low pri) |
| `min_supporting_snps` | STEP03 seeds | 5 | Consider span scaling (BA) |
| `d1..d12` thresholds | C01d tier counts | various | OK (documented above) |

---

## 5. Findings list (continuation of A..AN from chat 11)

- **AO** (closed): C01d header drift — says "10 dimensions" in one
  place, "12" elsewhere. Mis-doc, no code bug.
- **AP** (open, chat 14): C01d pattern enum missing `partial_overlap_neighbor`.
- **AQ** (open, chat 13): C01e Panel H path resolution after C01j
  directory move not verified.
- **AR** (open, chat 14): C01g + C01l integration — real payoff but
  needs a `boundary_v2` schema.
- **AS** (open, chat 14): C01j structure_score cutoffs 0.35/0.6 —
  recalibrate on Quentin's cohort.
- **AT** (open, chat 13): C01l `flank_bp = 500 kb` — scale with
  candidate span.
- **AU** (closed by chat 12): decompose's internal fallback path
  (`km$seeded = FALSE`) now emits `no_seeding`.
- **AV** (open, chat 13): decompose silhouette cutoff 0.25 → 0.40.
- **AW** (open, chat 14): multi_recomb `min_dco_bp = 200 kb` —
  recalibrate on HPC run.
- **AX** (closed, doc): nested_composition and C01l are complementary,
  not conflicting — wire both.
- **AY** (open, trivial): lib_step03_seed_loader doc drift — drop-
  conflict is stricter than priority-flashlight.
- **AZ** (open, low pri): STEP03 seed p-value FDR.
- **BA** (open, low pri): STEP03 seed support scaling by span.
- **BB** (open, chat 13): ghsl_confirmation assertion on annot_dt
  uniqueness.
- **BC** (open, chat 13): compute wired-fraction of key_spec.
- **BD** (open, chat 13): C01k path verification against registry
  block locations.
- **BE** (open, chat 14): add 3 pattern-enum buckets.
- **BF** (open, doc — manuscript-level): terminology. What the
  gene_conversion_detector surfaces is an *arrangement-discordant IBS
  tract* (short stretch where a sample's diagnostic-SNP dosage matches
  the opposite arrangement, then reverts). Gene conversion is the
  classical interpretation (Harris 1993) but is not provable from
  sequence data alone — short double crossover, paralog mismapping,
  and plain coincidental IBS at a few diagnostic sites produce the
  same signal. Without trios / mutation-accumulation panel / long-read
  haplotypes, we can surface the tract but not its cause. The block
  and file names stay as `gene_conversion_tracts` for pipeline
  back-compat (renaming would ripple through many readers); the
  SCHEMA DESCRIPTION and the DETECTOR HEADER have been updated in
  chat 12 to say "arrangement-discordant IBS tracts consistent with
  gene conversion". Manuscript text and figure legends should use
  the honest framing: "GC-like tracts" or "arrangement-discordant
  IBS tracts" rather than "gene conversion events". No code or math
  change; the detector's job is tract surfacing, cause is downstream.

Chat-11 findings status check:

- **AL** (open, handoff): rename `utils/sample_registry.R::load_registry` →
  `load_sample_groups_api`. Still open. Deferred; not touched this
  chat because it's orthogonal to the DAG/GC work.
- **T** (open, doc): README structure_type drift.
- **V** (closed): cheat24 inline vs real — deleted by chat 12.
- **W** (deferred to chat 15): register_C01f_keys helper verification.
- **Z** (closed): 4d README stale — doc pass happens chat 15.
- **8** (will be closed by chat 13): 4a writes no registry blocks.
- **9** (deferred to chat 15): C01g silent skip on missing helper.

---

## 6. Recommendations for chat 13 wiring

**Safe to wire immediately:**

1. C01j regime_segments output → `write_block("regime_segments", ...)`.
   Schema frozen, compute logic trusted, just redirect output.
2. C01l local_structure_segments → `write_block`. Same.
3. C01m distance_concordance → `write_block`. Same.
4. C01i_b_multi_recomb already writes `regime_sample_dag` via chat-12
   patch. Confirm the registry's JSON writer handles the `per_sample`
   data.table; if not, serialise via `as.list(…)` first.
5. C01d scoring picks up `regime_dominance` + `boundary_sharpness` +
   `inversion_vs_family_score` into D13–D15 as annotation-only
   columns (non-gating) on the first wiring pass. Gate promotion can
   happen chat 14 after HPC calibration.
6. Seal reads `recomb_subgroup` from multi_recomb and registers
   `inv_<cid>_RECOMBINANT_GC` / `inv_<cid>_RECOMBINANT_DCO` groups.
7. Characterize Q2 pulls the new keys listed above.

**Needs one more design pass before wiring:**

- `flank_bp` scaling in C01l (Finding AT) — decide the scaling rule
  before wiring so the block's `params` field is stable.
- Silhouette cutoff for decompose quality flag (Finding AV) — set to
  0.40 or pick based on a LG12-vs-LG25 dry run.

**Ready for HPC (chat 14):**

- End-to-end: phase 2 → phase 4a → phase 4b → phase 4c → phase 4e,
  chromosome-by-chromosome. LG12 first (clean inversion), LG25
  second (paralog-heavy control), then sweep.

---

## 7. What this audit did NOT do

- Did not re-run the v11 end-to-end suite. The DoD item 8 said
  "rerun v11 end-to-end with DAG-based R derivation" — the DAG
  rewrite has unit-level test coverage (50/50) but not integration
  coverage. That requires the HPC environment.
- Did not audit each of the 9 cheats in 4d.
- Did not verify each of the 367 key_spec entries in
  compute_candidate_status. Fraction-wired computation is flagged
  for chat 13.
- Did not trace Panel H path resolution in C01e (flagged AQ).

---

## 8. Closing note

Chat 11 was the architecture pass (registry library buildout). Chat
11.5 was the mechanism pass (decompose redesign). Chat 12 is the
coherence pass that caught: the over-counting in derive_R_from_regime;
the scale mismatch in gene_conversion_detector (60 kb floor vs 5 kb
biology); the silent pollution path in decompose's unsupervised
fallback; and the inline cheat24 divergence (Finding V).

Each of those would have produced wrong numbers in the first HPC run
if left alone.

Chat 13 can wire with confidence: the DAG semantics are verified by
fixture, the GC detector is verified by fixture, and the combination
rule's two-gate structure (fraction AND bp) is the actual signal
discriminator against alternating noise.

The biggest remaining risk for chat 14 is C01j's `structure_score`
calibration on Quentin's specific cohort (Finding AS) — defaults
were chosen on LD-style theory, not measured. First HPC run on LG12
should produce the distribution; tune then.
