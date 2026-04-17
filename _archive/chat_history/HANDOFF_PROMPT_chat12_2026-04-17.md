# Handoff prompt for chat 12 — phase-4 coherence audit + DAG-based recomb + wiring

**SUPERSEDES** all prior chat-12 handoff drafts. This version reflects
Quentin's direction at end of chat 11 (two sessions condensed): after the
registry library work (chat 11) and the decompose/multi_recomb redesign
work (chat 11.5, the decompose-v2 / C01j dispatch session), the next chat
needs to STOP adding new mechanism and DO THREE THINGS in strict order:

1. **Read and understand every phase-4 script end-to-end, function by
   function.** Not skim — understand. Criticize where logic is weak,
   comment on what each function actually does vs. what it's documented
   to do. Identify every drift, every dead parameter, every
   almost-but-not-quite wiring. Produce a written audit.
2. **Fix the regime-change counting to match Quentin's DAG-based
   formulation** (details in "DAG-based recombinant semantics" below).
   This is a real logic change, not just a bug. Do it as part of the
   audit because it will interact with other corrections.
3. **Only AFTER the audit + fix**, wire everything into the registries.
   Do not start wiring before the audit is complete. Wiring the wrong
   mental model into the registry would be worse than not wiring at all.

After chat 12 finishes, chat 13 does the first HPC run.

---

## The reality of what's in this tarball

You get:

- **Everything from chat 11** (registry library extensions): sample +
  interval registry APIs with 22 new methods, test fixtures, API
  cheatsheet, orphan duplicates removed. All wired into
  `registries/api/R/registry_loader.R` (1328 lines now).
- **Four orphan phase-4 scripts dispatched into proper homes** (chat
  11.5): STEP_C01j (regime compatibility engine),
  STEP_C01l (local structure ENA/Δ12), STEP_C01m (distance concordance),
  STEP_C01k (annotated simmat integrative figure). These were written
  BEFORE the v10 rewrite and forgotten. They are the correct tools for
  recombinant detection, boundary characterization, family-vs-inversion
  disambiguation, and per-chromosome integrative reporting. The logic in
  them is trusted — tune and wire, don't rewrite.
- **Gene conversion detector** (chat 11.5 invention, kept with correct
  scope): `4b_group_proposal/gene_conversion_detector.R`. Detects short
  (1–5 bin) dosage excursion tracts — NOT recombinants. Windowed binning
  (40 SNPs / 10-SNP step default) avoids per-SNP noise. This module
  exists because the first draft of a recombinant detector was coded
  using per-SNP CUSUM — that design was wrong for recombinants but
  exactly right for gene conversion, so it was kept with narrowed scope.
- **Two partial libraries** ported forward: `lib_step03_seed_loader.R`
  (Tier-1 enhancement: flashlight ∪ STEP03 seeds, no unsupervised
  fallback), `lib_ghsl_confirmation.R` (Tier-3 confirmation).
- **A partial combination rule**: `lib_recomb_combination.R` drafted per
  chat-9 design, but its `derive_R_from_regime()` uses a naive
  "consecutive-long-segments-with-different-values" counter. Quentin
  explicitly pushed back: the correct semantics are DAG-based. That's
  the primary code change in chat 12.
- **Four new block schemas**: regime_segments, local_structure_segments,
  distance_concordance, gene_conversion_tracts. Contracts are frozen.
  Don't change.

---

## DAG-based recombinant semantics (PRIMARY CODE CHANGE)

### The scale spectrum: GC and CO/DCO are the same signal at different scales

This framing was made explicit at end-of-chat-11.5 and should guide the
chat-12 audit of the two detectors' parameter choices. CRITICAL
CORRECTION from end-of-chat: the first GC detector draft was wrong
about WHERE to do the detection. Re-specified below.

**Gene conversion (GC):**
- Same signal TYPE: sample's arrangement state switches within the interval
- BIOLOGICALLY VERY SHORT: 50 bp to ~1 kb typical, occasional up to 5 kb
- VISIBLE in a sample × SNP heatmap (samples on Y, SNPs on X) as: a
  SHORT HORIZONTAL STRIPE. One sample row shows the opposite-arrangement
  color across 1–3 consecutive SNP columns, then returns to its group's
  color. Horizontal because "one sample across a few SNPs" is one row
  spanning a few columns.
- VERTICAL stripes in such a heatmap are NOT GC events. They arise from:
  (a) non-diagnostic SNPs — columns where REF and INV groups happen to
      carry the same minor allele frequency; these columns don't
      discriminate arrangement at all and must be filtered before GC
      detection
  (b) paralog mismapping — reads from a paralogous locus mismap onto
      the reference at that position, causing many samples to have
      anomalous genotype calls at the same site; common in
      paralog-dense regions (e.g. immune gene clusters near candidate
      48 on LG25)
  (c) SNPs with excess heterozygosity (Hardy-Weinberg-violating sites)
      — genotyping artifacts that look like many samples being
      heterozygous at the same position
  The audit should flag that candidate 48 on LG25 (the diagnostic
  example candidate) sits near IGHM/IGH/IGHD/IGHA — a known
  paralog-heavy region. Most of the vertical stripes visible in its
  marker heatmap are likely paralog artifacts, NOT GC events.
- At Quentin's SNP density (~0.6 SNP/kb = ~1.6 kb per SNP), a 1 kb GC
  tract is AT MOST 1 SNP wide. A 2 kb tract is 1–2 SNPs. A 5 kb
  outlier tract is ~3 SNPs.
- Therefore GC MUST be detected at SNP resolution, with a run-length
  prior along the GENOMIC AXIS within a SINGLE SAMPLE ROW — NOT by
  windowing SNPs into bins of 40 and taking means.
- The first draft's "40-SNP bin" detector has a lower detection limit
  of ~60 kb (40 SNPs × 1.6 kb/SNP), which is larger than the biology.
  Real GC events are invisible to that detector. **DELETE the binned
  approach. REPLACE with per-SNP run-length detection (spec below).**

**Recombinant / crossover (CO) / double crossover (DCO):**
- Same signal TYPE: within-sample arrangement state change
- LARGE: tens of kb to hundreds of kb
- OFTEN MANY samples: segregating recombinant haplotype
- Detected by C01j compatibility groups at window scale (50-marker
  windows, step 10). C01j's min_run_windows gate naturally filters
  GC-scale events (too short to register as regime segments).

**Scale handoff:**
- GC detector: catches runs of 2–~10 flagged SNPs = ~3 kb to ~20 kb
  (wide range because SNP density varies across genome).
- C01j regime: catches regime changes persisting ≥ 3 windows ≈ ≥150
  markers ≈ ≥ 250 kb at Quentin's density.
- **Gap between ~20 kb and ~250 kb:** events in this gap are detected
  by NEITHER detector. Chat 12 should decide whether this gap matters
  biologically (is there a known class of intermediate-scale tracts
  in teleost inversions?) or if it's acceptable. If not, C01j's
  min_run_windows can be lowered to 2, shrinking the gap.

### Revised gene-conversion detector spec (REWRITE in chat 12)

DELETE the current `gene_conversion_detector.R` content (windowed
binning + max_tract_bins gate — wrong scale). REWRITE with the
per-SNP run-length algorithm below.

**Input per candidate:**
- Per-SNP genotype matrix (samples × SNPs) inside the candidate interval
- Tier-1 decompose class per sample: HOM_REF / HET / HOM_INV / AMBIGUOUS

**Step 1 — SNP QC pre-filter.** Before any per-sample analysis, filter
out SNPs that are likely to produce false positives. Exclude SNPs where:
- `n_samples_genotyped < 0.8 × total_cohort` (too much missingness)
- `observed_het_fraction > 0.7` (Hardy-Weinberg-violating; likely
  paralog-mismapping artifact)
- `depth_z_score > 3` or `< −3` (coverage anomaly — duplicated or
  deleted locus)
- `maf < 0.02` (rare alleles can't discriminate arrangement)

These filters remove the SOURCES of vertical stripes in the heatmap
that have nothing to do with GC. This step is essential for regions
like LG25:13.78–14.11 Mb near immune gene paralogs.

**Step 2 — Identify diagnostic SNPs.** A SNP passes the QC filter AND
has allele frequency differing between HOM_REF and HOM_INV groups by
≥ `min_delta_af` (default 0.5). Only diagnostic QC-clean SNPs carry
GC information. Write both the QC mask and the diagnostic mask to
the block — useful for downstream explanation of what's filtered and
why.

**Step 3 — Per-sample per-SNP flagging.** For each HOM_REF or HOM_INV
sample (skip HET — ambiguous per-haplotype assignment at 9×) at each
diagnostic QC-clean SNP:
- HOM_REF sample: consistent = dosage ≤ 0.5, inconsistent = dosage ≥ 1.5
- HOM_INV sample: consistent = dosage ≥ 1.5, inconsistent = dosage ≤ 0.5
- HET-looking intermediate dosages (0.5 < d < 1.5) count as neither
  flag nor consistent — they're tolerance states
- NA/missing = neither

**Step 4 — Run-length aggregation with tolerance.** Walk flagged SNPs
in genomic order. Define a candidate tract as a maximal run satisfying:
- Starts and ends with a flagged SNP
- Contains ≥ `min_run` consecutive flagged SNPs (default 2)
- Allows up to `max_tolerance` intervening non-flagged-non-consistent
  SNPs (default 1 — a tolerance for one HET-looking or missing call
  inside the tract)
- Run "ends" when encountering `max_tolerance + 1` consecutive
  non-flagged SNPs

**Step 5 — Length gate (upper bound, same spirit as before).**
- Tract of 2–~10 flagged SNPs with span ≤ `max_span_bp` (default
  20 kb): call it a GC tract candidate
- Tract > `max_span_bp` or > `max_flagged_snps` (default 10): NOT GC,
  discard (it's a recombinant; C01j handles it)

**Step 6 — Confidence from run length.** Geometric prior over
per-SNP genotyping error rate p ≈ 0.01 at 9×:
- P(k consecutive errors at diagnostic SNPs) ≈ p^k
- run_length = 2: LOW confidence (~0.01%, but k=2 is common so
  expected false positive count across the genome is non-trivial)
- run_length = 3: MEDIUM (~0.0001%)
- run_length ≥ 4: HIGH (<0.000001%)

Emit a `confidence` field per tract. Downstream can filter on
confidence. This is a more honest representation than a hard cutoff.

**Step 7 — Direction annotation.** Same as the current
`direction` field: `REF_in_INV_context`, `INV_in_REF_context`,
etc. Derived from sample's baseline class + tract dosage pattern.

**Step 8 — Per-cohort summary.** Fraction of samples with at least
one tract, total tract count, median tract length in SNPs and in bp.
Aggregates for the block's top-level fields.

This algorithm IS what's visible in the heatmap: runs of 2–3 SNP
columns where one sample is the opposite color. It matches the
biology (short tracts). It is simple to implement (~150 lines).

**Update the schema** `gene_conversion_tracts.schema.json` to carry
the new per-tract fields: `run_length_flagged_snps`, `tolerance_snps`,
`confidence` ∈ {LOW, MEDIUM, HIGH}, `span_bp`.

### What was wrong with the end-of-chat-11.5 draft

`derive_R_from_regime()` counted "number of transitions between consecutive
regime segments of length >= min_run_windows that have DIFFERENT values."
For a sample whose regime track is `A A A B B A A A A A A A` (min_run=3,
B-run length=2), the code correctly drops the B burst as noise — good.
But for `A A A B B B A A A A A A`, both A and the middle B are long
enough, so it counts one transition from A→B→A as 2 changes. That
over-counts, and the rigid "A→B→A" shape requirement is too strict.

### Quentin's DAG formulation (correct)

For each sample, build a minimal directed graph from the per-window
regime labels inside the candidate interval:

- **Nodes** = distinct regime labels the sample visits (after collapsing
  consecutive same-labels into runs; no min_run filter — keep short
  bursts as nodes but annotate them as short).
- **Edges** = transitions between consecutive distinct runs.
- **Node weights** = number of windows occupied by that run.
- **Edge order** = the position of the transition in the interval.

Per-sample summary:

- `dominant_regime` = argmax over run windows. "The full inversion is
  basically A" means dominant_regime = A for most samples, with some
  deviating.
- `n_distinct_nodes` = graph size (number of distinct labels visited).
- `n_edges` = number of transitions (for A-A-A-B-B-A: 2 transitions
  A→B, B→A).
- `terminal_matches_start` = TRUE if first run label == last run label.
  A recombinant's DAG might NOT start and end on the same label —
  that's actually a distinguishing feature vs. noise (noise tends to
  return to the dominant).
- `deviation_fraction` = 1 − (dominant_run_windows / total_windows).
  This is the key scalar: how much of this sample's interval deviates
  from its own dominant regime.
- `has_any_deviation` = deviation_fraction > deviation_threshold (e.g.
  0.10 — 10% of interval). Soft threshold, not a hard min_run cutoff.
- `longest_deviation_bp` = widest contiguous stretch of non-dominant
  regime.

### Cohort-level signal at candidate level

Per-sample `has_any_deviation` aggregates to:

- `n_samples_with_deviation` — how many samples show ANY regime
  deviation in the interval.
- `fraction_samples_deviating` = n_samples_with_deviation / n_samples.
- `candidate_regime_dominance` = single-regime call: "this interval is
  basically A with X% samples showing deviation."

The final inversion call uses BOTH layers:

- **Per-sample R signal** for the combination rule with G (GHSL):
  `R_fired = has_any_deviation AND deviation_fraction >= 0.15` (or
  similar — a sample with 15%+ of windows in a non-dominant regime
  AND deviations long enough to not be pure noise).
- **Candidate-level breadth** goes into the candidate's evidence block
  as a scalar feature feeding C01d Layer C scoring: an interval where
  90% of samples show some deviation is a more complex system than
  one where 5% do. This is `fraction_samples_deviating`.

### What the code change looks like

Rewrite `derive_R_from_regime()` to build these DAGs and produce per-
sample + candidate-level summaries. Approximately:

```r
derive_R_from_regime <- function(regime_memb, interval,
                                  min_deviation_fraction = 0.15,
                                  min_deviation_bp = 50000L) {
  # For each sample:
  per_sample <- regime_memb[..., by = sample_id, {
    # Build run-length encoding of group_id (no min-run filter)
    rle_g <- rle(group_id)
    run_lengths <- rle_g$lengths
    run_labels  <- rle_g$values

    # Dominant regime = most-windows label
    dominant <- run_labels[which.max(
      tapply(run_lengths, run_labels, sum)
    )]  # or similar

    n_dom_windows <- sum(run_lengths[run_labels == dominant])
    total_windows <- length(group_id)
    deviation_fraction <- 1 - n_dom_windows / total_windows

    # Longest deviation run
    non_dom_lengths <- run_lengths[run_labels != dominant]
    longest_dev_windows <- if (length(non_dom_lengths)) max(non_dom_lengths) else 0L
    # Map to bp using pos_mid_mb diffs or explicit window bp

    # DAG-based shape features
    n_distinct_nodes <- length(unique(run_labels))
    n_edges <- if (length(run_labels) <= 1) 0L else length(run_labels) - 1L
    terminal_matches_start <- run_labels[1] == run_labels[length(run_labels)]

    list(
      dominant_regime = dominant,
      deviation_fraction = deviation_fraction,
      longest_dev_windows = longest_dev_windows,
      n_distinct_nodes = n_distinct_nodes,
      n_edges = n_edges,
      terminal_matches_start = terminal_matches_start,
      R_fired = (deviation_fraction >= min_deviation_fraction)
                  & (longest_dev_windows * window_size_bp >= min_deviation_bp)
    )
  }]

  # Cohort-level
  cohort <- list(
    n_samples_total         = uniqueN(regime_memb$sample_id),
    n_samples_with_deviation= sum(per_sample$deviation_fraction > 0),
    fraction_samples_deviating = mean(per_sample$deviation_fraction > 0),
    dominant_regime_cohort  = ... mode-of-per-sample-dominants
  )

  list(per_sample = per_sample, cohort = cohort)
}
```

Both outputs flow into the combination rule AND into a new evidence
block (proposed schema: `regime_sample_dag`) that feeds C01d scoring.

### Tiny per-sample DAG plot

For diagnostics, write a per-sample tiny-graph plotter: the sample's
run-length-encoded regime track rendered as a horizontal strip with node
rectangles sized by run length, color-coded by regime label. Use
`ggplot2` + `ggraph` or plain base graphics (prefer plain R — fewer
dependencies). Output: one PDF panel per candidate, a grid of all
cohort samples' tiny DAGs. Useful for manual review of recombinant
calls.

---

## The phase-4 audit task

This is the biggest deliverable for chat 12. Read every file in
`inversion_modules/phase_4_postprocessing/`, understand what each
function does and why, criticize and comment, produce a written audit.
Not a surface skim — understand the math, the thresholds, the contract
every function implements vs. advertises.

### Reading order

**4a_existence_layers — the scoring and boundary layers:**

1. `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` (~1050 lines).
   12-dimension scoring (D1–D12). Pattern classification
   (nested_fixed / nested_rare / complex_system / strong_inversion /
   etc.). Tier assignment. This is the existence-layer-a writer.
   Questions to answer:
   - What does each of D1–D12 actually measure?
   - How are the Layer A / Layer B scores combined into `final_score`?
   - Where does `n_children` come from (it's not populated upstream yet)?
   - Does the pattern logic handle the Inv14.2-style partial-overlap
     case? (Probably no — flag as a needed extension.)
2. `STEP_C01e_candidate_figures.R`. Per-candidate diagnostic plot
   composer. Panel H reads C01j regime output — **confirm this wiring
   is still correct after C01j's move into 4a/** and update the hardcoded
   paths if not.
3. `STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`.
   Boundary catalog. Writes `boundary_left` and `boundary_right` blocks.
   Currently single-point estimates. **In light of C01l's per-segment
   ENA/Δ12, should C01g start consuming C01l output to compute fuzz
   width?** Answer this in the audit.
4. **NEW: `STEP_C01j_regime_compatibility_engine.R`**. 927 lines,
   sliding-window compatibility clustering. Read `compute_compatibility_groups`
   (L136), `track_regimes` (L224), `segment_regimes` (L315),
   `analyze_membership_stability` (L424). Verify the thresholds
   (structure_score cutoffs at 0.35 / 0.6) match what Quentin's cohort
   actually needs. Output paths need updating — currently writes to
   `--outdir` standalone; needs to write via the registry's
   `write_block("regime_segments", ...)` instead.
5. **NEW: `STEP_C01l_local_structure_segments.R`**. 562 lines, per-
   candidate 5-segment ENA/Δ12 computation. Read `compute_local_memberships`,
   `compute_metrics`, `define_segments`. Verify the segment definitions
   (left_flank uses `flank_bp=500kb` default — is this sensible for
   Quentin's ~Mb-scale inversions?). Same wiring note: standalone output
   now, needs registry `write_block("local_structure_segments", ...)`.
6. **NEW: `STEP_C01m_distance_concordance.R`**. 567 lines, multi-scale
   pair concordance. Read `compute_concordance_at_distance`,
   `compute_concordance_fast`. Verify: the `DISTANCES = 80,160,320,640`
   defaults match Quentin's window size. Same wiring note.

**4b_group_proposal — decomposition and recombinant calling:**

7. `STEP_C01i_decompose.R` (452 lines). Tier 1 k-means on mean_pc1.
   Read `seeded_kmeans` (L95), main loop L254+. Critical task for
   chat 12: **modify to use `lib_step03_seed_loader.R` instead of the
   inline seed logic.** Remove the unsupervised k-means fallback (per
   chat-9 design). When no seeds available, emit `decomp_status="no_seeding"`
   and skip. Don't ship noisy guesses.
8. `lib_decompose_helpers.R`. Shared utilities. Understand what's
   actually used by decompose + multi_recomb.
9. `STEP_C01i_b_multi_recomb.R` (436 lines). Read
   `find_switch_segment` (L143), `count_phase_switches` (L203),
   `flashlight_hemi_signal` (L227), and the old combination rule.
   Critical task for chat 12: **rewire to consume
   `lib_recomb_combination.R` with the new DAG-based
   `derive_R_from_regime`.** Remove the inline 100kb mosaic-length
   threshold. The old signal 1 (PCA per-window switch) becomes Tier 4
   diagnostic, no longer gating. New gate: R (from C01j regime) AND G
   (from GHSL).
10. `STEP_C01i_c_nested_composition.py`. Nested composite flag. Already
    uses delta12/entropy/ena from local_Q. Confirm it doesn't conflict
    with C01l (which does per-segment ENA). **Audit question: should
    C01l's per-segment ENA replace nested_composition's interval-averaged
    ENA, or do they complement?** Decision goes in the audit.
11. `STEP_C01i_d_seal.R` (422 lines). Read `register_all_groups`
    (L221). Understand the HOM_STD alias, subgroup registration, the
    promotion_cap logic. Does seal need to read the new combination
    rule's output for RECOMBINANT/RECOMBINANT_GC/RECOMBINANT_DCO
    subgroup names? Confirm the per-sample `recomb_status` from the
    new combination rule flows into seal's `class_groups` correctly.
12. `gene_conversion_detector.R`. Narrow-scope GC detector. Confirm
    the windowing parameters (40 SNPs default) are sensible for
    Quentin's density.
13. `lib_step03_seed_loader.R`. Read the conflict-drop rule. Confirm
    the STEP03 output path and format match what chat-9 §Data
    dependencies assumed.
14. `lib_ghsl_confirmation.R`. Read the resolve_karyo_bp logic. The
    window-index to bp mapping relies on `annot_dt` having one row per
    global_window_id in chrom order — verify this is actually what
    Snake 3 writes.
15. `lib_recomb_combination.R`. After the DAG-based rewrite of
    `derive_R_from_regime`, re-verify the combination rule still makes
    sense. Does R_fired from a DAG-based computation mean the same
    thing as R_fired from the old consecutive-long-segments counter?
    Update thresholds and wording if needed.

**4c_group_validation — hypothesis tests:**

16. `STEP_C01f_hypothesis_tests.R` (2523 lines). This is the biggest
    consumer of the sample registry. Read `comp_from_registry` (L300)
    — this is the function that should use the chat-11
    `get_groups_for_candidate` method (4 has_group + 4 get_group calls
    collapse to 1). Convert it in this audit, since it's a pure
    simplification.
17. `group_validation_gate.R` (143 lines). Gate logic: NONE /
    UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT. Where does the
    registry's `q6_group_validation` key get read? Trace it
    end-to-end.

**4d_group_dependent — cheats:**

18. Look at the cheat structure (~9 cheats). Don't audit each in
    depth; just confirm the launcher dispatches correctly and that
    each cheat writes its own block type. Are any missing schemas?

**4e_final_classification — characterization and final tier:**

19. `compute_candidate_status.R` (925 lines). Final tier + verdict
    computer. Read `build_key_spec` (L40–413) to confirm the 367-key
    spec still matches reality. The 73.8% wired figure from chat 10
    — is it higher now with chat-11 additions? Chat-12 additions
    will add regime/local_structure/distance_concordance keys — what
    fraction does that bring us to?
20. `characterize_candidate.R` (664 lines). Q1–Q7 characterizers.
    Confirm Q2 (structure) now pulls from the new regime_segments +
    local_structure_segments + distance_concordance evidence blocks.
    **This is the biggest payoff of the chat-12 audit** — Q2 gets
    richer fields it didn't have before.
21. `run_characterize.R` (420 lines). Phase-4e driver. Smoke-tested
    by chat 10. Confirm the three-tier loader still works with the
    new keys.
22. **NEW: `STEP_C01k_annotated_simmat.R`** (317 lines). Per-
    chromosome integrative figure. Read carefully. What does it
    read? Does it read the actual new evidence blocks (regime,
    local_struct, distance_conc) or stub paths? Wire it properly.

### The audit document

Write `AUDIT_PHASE4_COHERENCE_2026-04-18.md` as the primary deliverable.
Structure:

1. **Per-script section.** One per file, containing:
   - What the file is supposed to do (one paragraph).
   - Function inventory with each function's inputs/outputs/purpose.
   - What the file actually does (where does actual behavior drift from
     documentation?).
   - List of parameter defaults that look suspect (thresholds, cutoffs)
     with Quentin's-cohort-appropriateness judgment.
   - Outputs produced vs. outputs consumed downstream — any dead
     outputs? Any inputs the file ASSUMES but doesn't actually
     validate?
   - Any silent failures or tryCatch-swallowed errors that hide
     real problems?
2. **Flow-of-evidence section.** Trace a single candidate from
   phase-2 discovery through phase-4e final tier. At each hop, what
   gets written, by whom, and who reads it? Where are the writes
   aspirational (no reader)? Where are the reads aspirational (no
   writer)?
3. **Pattern classification coverage.** C01d's pattern enum
   (strong_inversion / nested_fixed / complex_system / etc.). Which
   ones actually get emitted in practice? Which are never reached?
   Does the chat-12 regime+local_struct+distance_conc evidence
   support additions to the enum (e.g. `partial_overlap_neighbor`
   for the Inv14.2 case)?
4. **Threshold drift list.** Every magic number in the codebase, with
   its source and whether it's defensible. Common drifters:
   - 100kb mosaic threshold (deleted in chat-9 design — verify it's
     really gone)
   - cheat24 thresholds for GC vs DCO vs suspicious
   - silhouette cutoff 0.25 in decompose quality flag
   - bic_gap cutoff 0.05
   - phase_concordance 0.30
   - structure_score cutoffs 0.35 / 0.6 in C01j
   - flank_bp 500kb in C01l
   - distances 80,160,320,640 in C01m
   - deviation_fraction threshold for R_fired (new, needs to be set)
5. **Findings list.** Continue the A–AN series from chat 11.
   Every issue found gets a new finding letter.
6. **Recommendations for chat 13 wiring.** What's safe to wire, what
   needs one more design pass, what's ready for HPC.

### Criticize, don't soft-pedal

When you find something weak — say so. "The silhouette threshold of
0.25 is too lenient for 9x coverage; it should probably be 0.4."
"The cheat24 inline fallback disagrees with the real cheat24 by 2×
on the gene_conversion posterior — same input, different output."
Quentin's project benefits from honest criticism of its own code.
Don't defend design choices that don't hold up.

---

## Files chat 12 SHOULD touch

- `STEP_C01i_decompose.R`: add STEP03 seed wiring, remove unsupervised
  fallback, emit `no_seeding` status. Minimal diff.
- `STEP_C01i_b_multi_recomb.R`: rewire to use `lib_recomb_combination.R`,
  remove 100kb threshold, remove old Signal 1 (PCA per-window),
  register RECOMBINANT/RECOMBINANT_GC/RECOMBINANT_DCO subgroup names
  per the new combination output.
- `STEP_C01f_hypothesis_tests.R` `comp_from_registry`: replace 4
  has_group + 4 get_group with 1 `get_groups_for_candidate` call.
  Small, localized.
- `STEP_C01j_regime_compatibility_engine.R`,
  `STEP_C01l_local_structure_segments.R`,
  `STEP_C01m_distance_concordance.R`,
  `STEP_C01k_annotated_simmat.R`: REHOME outputs via registry
  `write_block`. Current code writes standalone TSVs; next chat wires
  them to the evidence registry. BUT: don't rewrite the compute logic.
  Just redirect output.
- `lib_recomb_combination.R`: rewrite `derive_R_from_regime` with
  DAG-based semantics (the primary code change described above).
- Add: `plot_sample_regime_dag.R` (new) for the tiny per-sample DAG
  diagnostic plots.
- Add: `regime_sample_dag.schema.json` for the per-sample DAG evidence
  block.

## Files chat 12 should NOT touch

- `registries/api/*/registry_loader.*`: mature after chat 11. Don't
  touch the API.
- `registries/schemas/structured_block_schemas/*.schema.json`: frozen
  contracts. Don't change existing schemas (regime_segments,
  local_structure_segments, distance_concordance, gene_conversion_tracts,
  boundary_scan are already written). Only ADD new schemas
  (regime_sample_dag). Never modify existing.
- Any phase_2_discovery script except to read their outputs. That's
  chat 13/14.
- `characterize_candidate.R` and `compute_candidate_status.R`:
  mature from chat 10. Only change if the audit surfaces a bug that
  blocks wiring — otherwise leave for chat 13.

---

## Wiring work for the session AFTER chat 12

Chat 13 does the actual registry wiring:

- C01j / C01l / C01m / C01k output → `write_block()` calls
- C01d scoring picks up regime_dominance + boundary_sharpness +
  inversion_vs_family_score from the new evidence blocks
- multi_recomb writes recombinant_map block via
  `lib_recomb_combination.R` final output
- seal reads the RECOMBINANT group from the updated multi_recomb and
  registers groups via the new combination rule's per-sample records
- phase 4e characterize gains Q2 rich fields (regime dominance,
  boundary sharpness, family-vs-inversion contrast)

Chat 14: first real HPC run on LANTA. Chromosome-by-chromosome.
Evaluate, iterate one more pass.

---

## Posture

**Read before code. Audit before wire. Understand before change.**

Chat 11 was architectural (library buildout). Chat 11.5 was mechanism
(decompose redesign via C01j + C01l + C01m). Chat 12 is the pause
before wiring — the coherence audit that catches what's wrong with
the plan before it gets blessed into the registry.

If the audit reveals that a chat-11.5 decision was wrong, change it.
If it reveals an old chat-5 decision was wrong, flag it and propose
a fix. The goal of chat 12 is to produce a plan that chat 13 can
execute mechanically without having to second-guess the architecture.

**No three-fix ceiling.** Audit fixes are open-ended; each one is
defensible against "does this block wiring or produce wrong numbers?"

---

## Remaining findings from chat 11 still open

- **AL**: load_registry shadowing between `utils/sample_registry.R` and
  `registries/api/R/registry_loader.R`. Fix in chat 12: rename
  `utils/sample_registry.R::load_registry` to
  `load_sample_groups_api` (or rename file entirely to
  `sample_groups_api.R`). Update `load_bridge.R` step 4 (grep for
  `load_registry(BRIDGE_PATHS$REGISTRY_DIR)`).
- **T**: README structure_type drift — doc pass.
- **V**: cheat24 threshold divergence between inline fallback and real
  cheat24. Fix by deleting the inline fallback (chat-9 design item).
- **W**: `register_C01f_keys` helper verification — deferred to chat 15.
- **Z**: 4d README stale — doc pass.
- **8**: 4a writes no registry blocks — fixed by chat-13 wiring.
- **9**: C01g silent skip on missing helper — chat 15.

New findings expected from the chat-12 audit: AO, AP, AQ, ...

---

## Definition of done for chat 12

1. `AUDIT_PHASE4_COHERENCE_2026-04-18.md` — primary deliverable, every
   script audited.
2. `STEP_C01i_decompose.R` patched: seed loader wired, unsupervised
   fallback removed.
3. `STEP_C01i_b_multi_recomb.R` patched: uses `lib_recomb_combination.R`,
   old Signal 1 demoted to Tier 4 diagnostic, 100kb threshold
   removed, Finding V (inline fallback) resolved.
4. `STEP_C01f_hypothesis_tests.R` patched: `comp_from_registry` uses
   `get_groups_for_candidate` in one call.
5. `lib_recomb_combination.R::derive_R_from_regime` rewritten with
   DAG-based semantics.
6. `plot_sample_regime_dag.R` new — per-sample diagnostic plot.
7. `regime_sample_dag.schema.json` new — the DAG evidence block.
8. Smoke tests: rerun the v11 end-to-end with DAG-based R derivation.
   Add a new fixture for the DAG shape check (verify on
   A-A-A-B-B-A and A-A-A-B-A-B-B-A patterns).
9. `AUDIT_LOG_chat12_2026-04-18.md` — findings AO+.
10. Update `SESSION_SUMMARY_2026-04-18.md` and
    `FIXES_APPLIED_2026-04-18.md`.
11. Tarball as `inversion-popgen-toolkit_chat12_audit_dag_2026-04-18.tar`.

End state: chat 13 can start by reading the audit and mechanically
wiring without architectural doubts.
