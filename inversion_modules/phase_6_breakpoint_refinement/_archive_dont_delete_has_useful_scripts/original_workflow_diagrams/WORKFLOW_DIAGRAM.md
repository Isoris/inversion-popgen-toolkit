# Workflow Diagram — Dosage-Based Breakpoint Refinement

```
══════════════════════════════════════════════════════════════════════════════

                     INVERSION BREAKPOINT REFINEMENT
                        Dosage-based, per-carrier-voted

══════════════════════════════════════════════════════════════════════════════


    ┌──────────────────────────────────────────────────────────┐
    │                      UPSTREAM INPUTS                      │
    │                                                           │
    │  1. DOSAGE_DIR/<chr>.dosage.tsv.gz                        │
    │     Pre-computed sample × SNP dosage matrix (BEAGLE)      │
    │                                                           │
    │  2. <cand_dir>/candidate_pca_rotated.tsv                  │
    │     From STEP21: u, v coords + coarse_group_refined       │
    │     (HOMO_1 / HET / HOMO_2)                               │
    │                                                           │
    │  3. CANDIDATE_TABLE                                       │
    │     chrom, start_bp, end_bp — coarse shelf from Fst       │
    └────────────────────────┬─────────────────────────────────┘
                             │
          ╔══════════════════╪══════════════════╗
          ║                  ▼                  ║
          ║    ┌─────────────────────────────┐  ║
          ║    │    STEP_BP_01               │  ║
          ║    │    DOSAGE SIGNAL            │  ║
          ║    │─────────────────────────────│  ║
          ║    │                             │  ║
          ║    │  • Identify informative     │  ║
          ║    │    markers: SNPs where      │  ║
          ║    │    |h2_mean - h1_mean| is   │  ║
          ║    │    large                    │  ║
          ║    │                             │  ║
          ║    │  • Build core block:        │  ║
          ║    │    markers with high mutual │  ║
          ║    │    correlation across       │  ║
          ║    │    samples (C01i logic)     │  ║
          ║    │                             │  ║
          ║    │  • Extend block boundaries  │  ║
          ║    │    outward one SNP at a     │  ║
          ║    │    time; stop when          │  ║
          ║    │    correlation with core    │  ║
          ║    │    drops (C01i logic)       │  ║
          ║    │                             │  ║
          ║    └──────────────┬──────────────┘  ║
          ║                   │                 ║
          ║                   ▼                 ║
          ║      candidate_dosage_blocks.tsv    ║
          ║      candidate_dosage_informative   ║
          ║                  _markers.tsv       ║
          ║                   │                 ║
          ║                   ▼                 ║
          ║    ┌─────────────────────────────┐  ║
          ║    │    STEP_BP_02               │  ║
          ║    │    ANCESTRAL FRAGMENTS      │  ║
          ║    │─────────────────────────────│  ║
          ║    │                             │  ║
          ║    │  For each inversion carrier │  ║
          ║    │  (HET or HOMO_INV) sample:  │  ║
          ║    │                             │  ║
          ║    │  • Walk leftward from block │  ║
          ║    │    edge one SNP at a time   │  ║
          ║    │  • At each SNP, test        │  ║
          ║    │    whether this sample's    │  ║
          ║    │    dosage still correlates  │  ║
          ║    │    with the core consensus  │  ║
          ║    │  • Stop when correlation    │  ║
          ║    │    drops = this sample's    │  ║
          ║    │    personal boundary        │  ║
          ║    │                             │  ║
          ║    │  • Same walk rightward      │  ║
          ║    │                             │  ║
          ║    │  OUTPUT: per-sample         │  ║
          ║    │  (frag_start_bp,            │  ║
          ║    │   frag_end_bp)              │  ║
          ║    │                             │  ║
          ║    │  AGGREGATE: the             │  ║
          ║    │  distribution of            │  ║
          ║    │  fragment boundaries        │  ║
          ║    │  across all carriers        │  ║
          ║    │  gives:                     │  ║
          ║    │    • modal position =       │  ║
          ║    │      true breakpoint        │  ║
          ║    │    • distribution spread =  │  ║
          ║    │      CI around breakpoint   │  ║
          ║    │    • tail distribution =    │  ║
          ║    │      recombinant history    │  ║
          ║    │                             │  ║
          ║    └──────────────┬──────────────┘  ║
          ║                   │                 ║
          ║ candidate_ancestral_fragments.tsv   ║
          ║ candidate_fragment_distribution.tsv ║
          ╚═══════════════════╪═════════════════╝
                              │
          ┌───────────────────┼───────────────────┐
          │    OPTIONAL SIGNALS (used if available, skipped if not)
          │                                                  │
          │  • C01i marker_coseg_blocks.tsv.gz                │
          │    → additional block spans for multi-system      │
          │                                                   │
          │  • C01j regime_transitions.tsv                    │
          │    → enter_structured / exit_structured boundaries│
          │                                                   │
          │  • C01l segment_summary.tsv.gz                    │
          │    → Delta_12 drop at flanks                      │
          │                                                   │
          │  • STEP40 candidate_internal_breaks.tsv           │
          │    → window-level coherence breaks                │
          │                                                   │
          │  • STEP41 candidate_switching_events.tsv.gz       │
          │    → per-sample switch position pile-ups          │
          │                                                   │
          │  • STEP37 breakpoint_support_per_candidate.tsv    │
          │    (BOTTOM LAYER — weight 0.5, user instruction:  │
          │     "put sv on top if you have energy")           │
          │                                                   │
          └──────────────────┬────────────────────────────────┘
                             │
                             ▼
          ┌─────────────────────────────────────────┐
          │        STEP_BP_03                       │
          │        CONSENSUS MERGE                  │
          │─────────────────────────────────────────│
          │                                         │
          │  Collect breakpoint estimates from:     │
          │    • STEP_BP_02 fragment distribution   │
          │    • STEP_BP_01 block extension         │
          │    • C01i / C01j / C01l / STEP40 /      │
          │      STEP41 / STEP37 (if available)     │
          │                                         │
          │  Weighted merge (per side):             │
          │                                         │
          │    Method                        Weight │
          │    ─────────────────────────────────── │
          │    ancestral fragments (primary) 3.0   │
          │    block extension               2.0   │
          │    C01j regime transitions       2.0   │
          │    STEP40 coherence breaks       1.0   │
          │    STEP41 switch pile-ups        1.0   │
          │    C01l Delta_12 drops           1.0   │
          │    STEP37 SV clusters            0.5   │
          │                                         │
          │  Weighted median across methods →       │
          │    final_left_bp, final_right_bp        │
          │  Weighted MAD across methods →          │
          │    CI bounds                            │
          │  Count of methods within 2×MAD →        │
          │    n_methods_agreeing                   │
          │                                         │
          └────────────────┬────────────────────────┘
                           │
                           ▼
        candidate_breakpoints_consensus.tsv
        candidate_breakpoints_per_method.tsv
        breakpoints_catalog.tsv   (catalog-level)
                           │
                           ▼
          ┌─────────────────────────────────────────┐
          │        STEP_BP_04                       │
          │        DIAGNOSTIC FIGURE                │
          │─────────────────────────────────────────│
          │                                         │
          │  One multi-track PDF per candidate:     │
          │                                         │
          │  Track 1 (top):                         │
          │    Smoothed |h2_mean - h1_mean| along x │
          │                                         │
          │  Track 2:                               │
          │    Per-carrier ancestral fragment       │
          │    boundary RUG (left side)             │
          │                                         │
          │  Track 3:                               │
          │    Histogram inset of boundary          │
          │    distribution — shows the MODE        │
          │                                         │
          │  Track 4:                               │
          │    C01j regime state RIBBON             │
          │                                         │
          │  Track 5:                               │
          │    C01l Delta_12 profile                │
          │                                         │
          │  Track 6:                               │
          │    STEP41 switch pile-up rug            │
          │                                         │
          │  Track 7 (bottom, low alpha):           │
          │    STEP37 SV cluster points             │
          │                                         │
          │  Vertical red lines + CI ribbons:       │
          │    final consensus breakpoints          │
          │                                         │
          └────────────────┬────────────────────────┘
                           │
                           ▼
         candidate_breakpoint_diagnostic.pdf


══════════════════════════════════════════════════════════════════════════════

 TIMING / LOCATION

   STEP_BP_01   lives in followup pipeline
   STEP_BP_02   lives in followup pipeline
   STEP_BP_03   lives in followup pipeline
   STEP_BP_04   lives in followup pipeline

   STEP_Q10_breakpoint_refinement.sh  is a shell wrapper that bridges
   phase_qc_shelf into the followup pipeline for one chromosome at a time.
   It runs 01 → 02 → 03 → 04 in sequence.

══════════════════════════════════════════════════════════════════════════════

 INTEGRATION TARGETS

   Downstream: ancestry registry (schema TBD — user to confirm)
               → loads from:
                   candidate_breakpoints_consensus.tsv   (per-candidate row)
                   candidate_ancestral_fragments.tsv     (per-sample-per-cand)
                   candidate_breakpoints_per_method.tsv  (per-method per-cand)

   Cross-linked with previous session outputs:
     STEP36 (ordering comparison heatmaps) — consumer of breakpoint consensus
                                             for top-annotation vertical lines
     STEP38 (composite figure)             — can be upgraded to use
                                             consensus breakpoints instead of
                                             Fst shelf edges

══════════════════════════════════════════════════════════════════════════════
```

## Why this order (as opposed to putting ancestry first)

The ancestry-fragment scan (BP_02) depends on having a core block to scan outward from. That core block comes from BP_01. Hence: block first, fragments second. The consensus step (BP_03) then merges ALL signals including ancestry, so ancestry is *primary in weight* (weight = 3.0) but *second in computation order*.

## Why per-carrier voting is the right idea

Traditional breakpoint callers (DELLY, Manta, Fst) produce one population-level estimate with no native CI. Here, each of N_carriers samples produces its own breakpoint estimate through its own history of recombination. The distribution across carriers carries both:
- a robust central estimate (modal position), which outperforms mean because recombinant tails are asymmetric
- a natural CI (spread of the distribution) that reflects true biological uncertainty — samples whose ancestors recombined closer to the breakpoint vote with less certainty, and the distribution width captures that
- diagnostic information about recombinant history (tail length, bimodality)

No other method in the user's codebase produces all three of these together.
