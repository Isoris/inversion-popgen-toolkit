# Pipeline Diagram — Inversion Analysis (complete picture)

```
══════════════════════════════════════════════════════════════════════════════

         INVERSION ANALYSIS WORKFLOW — 226-sample C. gariepinus cohort

         Three modules · Shared polarity · Manuscript pipeline visible

══════════════════════════════════════════════════════════════════════════════


        ┌─────────────────────────────────────────────────────────┐
        │   COHORT & DATA                                          │
        │                                                          │
        │   226 samples  ·  ~9× short-read Illumina                │
        │   28 chromosomes  ·  BEAGLE-imputed dosage               │
        │   Clair3 VCFs per sample (where available)               │
        │   Target: LG28 inversion, 15.115 – 18.005 Mb (~2.89 Mb)  │
        │           p ≈ 0.5  ·  60 / 106 / 60 karyotypes           │
        └─────────────────────────────────────────────────────────┘


        ┌─────────────────────────────────────────────────────────┐
        │   PREREQUISITES (from earlier pipeline steps)            │
        │                                                          │
        │   DOSAGE_DIR/<chr>.dosage.tsv.gz                         │
        │   DOSAGE_DIR/<chr>.sites.tsv.gz                          │
        │   candidate_pca_rotated.tsv    (STEP21 karyotype)        │
        │   CANDIDATE_TABLE              (coarse shelf-Fst hits)   │
        │   STEP29 polarity output       (L1 / L2 per marker)      │
        └─────────────────────────────────────────────────────────┘


   ╔═══════════════════════════════════════════════════════════════════╗
   ║                                                                   ║
   ║   SHARED LAYER — L1 / L2 POLARITY (STEP29)                       ║
   ║                                                                   ║
   ║     L1 per marker:  sign(mean(HOMO_2) − mean(HOMO_1))             ║
   ║                     anchored on homozygotes only                  ║
   ║                                                                   ║
   ║     L2 per block:   smooth L1 to blocks ≥5 consecutive            ║
   ║                     same-sign markers                             ║
   ║                                                                   ║
   ║     Application:    dosage_polarized = (dosage − 1) × L2 + 1      ║
   ║                                                                   ║
   ║     Canonical scale after polarization:                           ║
   ║       HOMO_1 ≈ 0   ·   HET ≈ 1   ·   HOMO_2 ≈ 2                   ║
   ║                                                                   ║
   ╚═══════════════════════════════════════════════════════════════════╝
         │                       │                        │
         ▼                       ▼                        ▼
     Module 1              Module 2                 Module 3
     heatmaps              breakpoints              regime


════════════════════════════════════════════════════════════════════════════════

   MODULE 1 — MARKER DOSAGE HEATMAPS             [ scaffolded in _archive/ ]

════════════════════════════════════════════════════════════════════════════════

       ┌────────────────────────────────────────────────┐
       │   within_group_ordering.R  (shared helper)     │
       │                                                │
       │   Row ordering:                                │
       │     primary   = u (STEP21 rotated PCA)         │
       │     secondary = hclust / v / residual_pc1 /    │
       │                 error_profile / none           │
       │                                                │
       │   Column ordering:                             │
       │     genomic position, polarized via L2         │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   STEP36 (to rename → STEP42)                  │
       │   V1/V2/V3 ordering comparison heatmaps        │
       │   Side-by-side panel per candidate             │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   STEP41 (to rename → STEP44)                  │
       │   q08b_shelf_heatmap_architecture.R            │
       │   Publication-ready heatmap rendering          │
       └────────────────────────────────────────────────┘

   STATUS:  scaffolded in _archive/heatmap_ordering_upgrade/
            naming collisions with user's STEP36/40/41 — rename first


════════════════════════════════════════════════════════════════════════════════

   MODULE 2 — BREAKPOINT FINDING + VISUALIZATION      [ main output of session ]

════════════════════════════════════════════════════════════════════════════════

                      ╔═══════════════════════════╗
                      ║   CORE PIPELINE           ║
                      ║   (run in order)          ║
                      ╚═════════════╤═════════════╝
                                    ▼
       ┌────────────────────────────────────────────────┐
       │   01 — DOSAGE SIGNAL                           │
       │   delta_i per marker                           │
       │   |delta| ≥ 1.0 → informative                  │
       │   core block via |cor| ≥ 0.7                   │
       │   extend outward until cor drops below 0.6     │
       │                                                │
       │   Output:  candidate_dosage_blocks.tsv         │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   02 — ANCESTRAL FRAGMENTS                     │
       │   (the primary breakpoint signal)              │
       │                                                │
       │   Per carrier: walk outward SNP by SNP         │
       │   Stop when sample's dosage stops correlating  │
       │   with the core consensus.                     │
       │                                                │
       │   Aggregate across carriers:                   │
       │     modal position    = breakpoint             │
       │     spread            = CI (bootstrap)         │
       │     tail              = recombinant ancestry   │
       │                                                │
       │   Output:  candidate_ancestral_fragments.tsv   │
       │            candidate_fragment_distribution.tsv │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   03 — CONSENSUS MERGE                         │
       │                                                │
       │   Weighted merger across sources:              │
       │     ancestral_fragments    wt 3.0  ●           │
       │     block_extension        wt 2.0              │
       │     c01j_transitions       wt 2.0  ○           │
       │     step40_coherence       wt 1.0  ○           │
       │     step41_switches        wt 1.0  ○           │
       │     c01l_segments          wt 1.0  ○           │
       │     step37_sv              wt 0.5  ○           │
       │                                                │
       │     ● primary (always present)                 │
       │     ○ optional (read if available)             │
       │                                                │
       │   Output:  candidate_breakpoints_consensus.tsv │
       │              ⇦ MAIN RESULT                      │
       │            candidate_breakpoints_per_method.tsv│
       │            breakpoints_catalog.tsv             │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   04 — DIAGNOSTIC FIGURE                       │
       │   Multi-track PDF per candidate                │
       └────────────────────────────────────────────────┘


                    ╔═══════════════════════════╗
                    ║   OPTIONAL DIAGNOSTICS    ║
                    ║   (Module 2)              ║
                    ╚═════════════╤═════════════╝
                                  ▼
       ┌────────────────────────────────────────────────┐
       │   05 — FOUR-ENCODING ROBUSTNESS                │
       │                                                │
       │   Sample × sample distance under 4 encodings:  │
       │     minor · major · 012 · polarized (L2)       │
       │                                                │
       │   ARI between pair clusterings                 │
       │     All high → structure is robust             │
       │     Any low  → that encoding diverges          │
       │                                                │
       │   NOT in consensus. Enable via RUN_05=1.       │
       └────────────────────────────────────────────────┘

       ┌────────────────────────────────────────────────┐
       │   07 — BREAKPOINT EVIDENCE PILEUP              │
       │   (panel D replica)                            │
       │                                                │
       │   Python. Per-sample stacked reads at both BPs.│
       │   Tracks:                                      │
       │     - split reads  (red ticks)                 │
       │     - discordant pairs  (FF/RR arcs)           │
       │     - coverage  (mean + 10-90% ribbon)         │
       │     - assembly alignment  (schematic only)     │
       │                                                │
       │   Samples ordered INV → HET → REF              │
       │   Labels anonymizable as sample1, sample2, ... │
       │                                                │
       │   Modes:                                       │
       │     A) pre-extracted evidence TSV              │
       │     B) fresh pysam from BAMs                   │
       │                                                │
       │   Runtime-tested with synthetic data           │
       │   Example: docs/example_panel_D_output.png     │
       └────────────────────────────────────────────────┘


════════════════════════════════════════════════════════════════════════════════

   MODULE 3 — REGIME / ANCESTRY ANALYSIS               [ partly scaffolded ]

════════════════════════════════════════════════════════════════════════════════

       ┌────────────────────────────────────────────────┐
       │   C01j regime compatibility engine             │
       │   (user's existing codebase, not in bundle)    │
       │                                                │
       │   Produces:                                    │
       │     regime_memberships.tsv.gz                  │
       │     regime_windows.tsv.gz                      │
       │     regime_transitions.tsv                     │
       └───────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌────────────────────────────────────────────────┐
       │   06 — REGIME STREAM GRAPH                     │
       │   Stacked-area plot across the chromosome      │
       │   Bands = compatibility clusters               │
       │   Matching across windows via Hungarian-style  │
       │   Inversions → stable 3-band zones             │
       └────────────────────────────────────────────────┘

       ┌────────────────────────────────────────────────┐
       │   [DEFERRED] Rolling seed-and-extend analysis  │
       │                                                │
       │   Slide along chromosome, extend while regime  │
       │   stable, record block boundaries + shape:     │
       │     Sharp flip    → inversion breakpoint       │
       │     Gradual       → recombinant edge           │
       │     Bubble        → gene conversion            │
       │     Unstable      → family LD transition       │
       │                                                │
       │   Strictly more informative than breakpoints.  │
       │   NOT CODED — 5 design questions first.        │
       └────────────────────────────────────────────────┘

       ┌────────────────────────────────────────────────┐
       │   [TO WIRE] unified_ancestry integration       │
       │                                                │
       │   C01j consumes Q-matrix via dispatcher:       │
       │     get_region_stats(what = c("Q"))            │
       │                                                │
       │   instead of computing own distances.          │
       │   Engine B (instant_q.cpp) becomes Q source.   │
       └────────────────────────────────────────────────┘


════════════════════════════════════════════════════════════════════════════════

   MANUSCRIPT PIPELINE           [ draft chunks in docs/MANUSCRIPT_CHUNKS.md ]

════════════════════════════════════════════════════════════════════════════════

   ┌──────────────────────────────────────────────────────────────────┐
   │                                                                  │
   │   METHODS       Detection floor statement:                       │
   │                 size ≥ 50 kb  ·  MAF ≥ 15%                       │
   │                                                                  │
   │                                                                  │
   │   RESULTS       Sensitivity sweep figure                         │
   │                   (size × MAF, LG28 subsampling, 20 replicates)  │
   │                                                                  │
   │                 Breakpoint repeat-context enrichment figure      │
   │                   (SDs / LINEs / SINEs vs genome-wide null)      │
   │                                                                  │
   │                                                                  │
   │   DISCUSSION    Comparison to deer mice (Peichel et al.)         │
   │                   "our 50 kb size floor beats their 1 Mb;        │
   │                    our 15% MAF floor is set by cohort size"      │
   │                                                                  │
   │                                                                  │
   │   LIMITATIONS   (1) MAF < 15% undetectable                       │
   │                 (2) Size < 50 kb undetectable                    │
   │                 (3) Gene conversion collapsed                    │
   │                 (4) SD-rich regions inflate CI                   │
   │                 (5) No genome-wide split-read null (noise floor) │
   │                                                                  │
   └──────────────────────────────────────────────────────────────────┘

   Future scripts (NOT for next chat, recorded in MANUSCRIPT_CHUNKS.md):
     - sensitivity_sweep.R     ~150 lines
     - breakpoint_repeat_enrichment.R   ~100 lines


════════════════════════════════════════════════════════════════════════════════

   EXECUTION ORDER

════════════════════════════════════════════════════════════════════════════════

   Module 2 core (always):
       ./run_pipeline.sh  <chrom>  <config.R>  <cid>
       ⟹  01 → 02 → 03 → 04

   With Module 2 optional robustness check:
       RUN_05=1  ./run_pipeline.sh  <chrom>  <config.R>  <cid>

   Module 2 panel D (anytime, standalone):
       python  07_breakpoint_evidence_pileup.py  \\
           --bam-dir … --sample-list …  --out figure.pdf

   Module 3 stream graph (anytime, standalone):
       Rscript  06_regime_stream_graph.R  \\
           --memberships … --windows …  --out …

   Module 1 (no wrapper yet — apply helper to user's scripts):
       see  _archive/heatmap_ordering_upgrade/


════════════════════════════════════════════════════════════════════════════════

   MINIMAL SET TO RUN

════════════════════════════════════════════════════════════════════════════════

       01_dosage_signal.R
       02_ancestral_fragments.R
       03_consensus_merge.R
       04_diagnostic_figure.R
       run_pipeline.sh

   Everything else is optional:
       05, 06, 07  →  diagnostics and visualizations
       Module 1    →  _archive/heatmap_ordering_upgrade/
       Module 3    →  awaiting unified_ancestry wiring


════════════════════════════════════════════════════════════════════════════════

   WHAT THE NEXT CHAT DOES           [ user mandate — WIRING, not new code ]

════════════════════════════════════════════════════════════════════════════════

   1. Registry adapter for each of the 7 scripts
      (read inputs by key, write outputs with provenance + sha256)

   2. Genome-location collector / driver
      (walks all candidates × all chromosomes, aggregates results)

   3. Cross-module API
      (scripts ask the registry for inputs, no hardcoded paths)

   4. Module 3 plugs into unified_ancestry via dispatcher
      (stop computing own Q-matrices)


       TRAP TO AVOID
       "let's improve X while we wire." NO.
       Wiring is mechanical. It must finish in ONE session.
       Improvements go in later sessions.


   FILES USER MUST UPLOAD AT START OF WIRING CHAT

       - DATABASE_DESIGN.md              (registry schema)
       - results_registry schema         (current state)
       - region_stats_dispatcher.R
       - instant_q.R                     (or the cpp + wrapper)
       - current latest  C01j
       - one example  C01f STEP02 output TSV
       - one example  config_inversion_followup.R


════════════════════════════════════════════════════════════════════════════════

   DOCUMENTATION LAYOUT

════════════════════════════════════════════════════════════════════════════════

   breakpoint_pipeline/
     README.md                          ⟵ read first, plain-English
     PIPELINE_DIAGRAM.md                ⟵ this file
     run_pipeline.sh                    ⟵ the wrapper
     01 – 04                            ⟵ Module 2 core
     05                                 ⟵ Module 2 robustness
     06                                 ⟵ Module 3 visualization
     07 + test_07                       ⟵ Module 2 panel D + smoke test

     docs/
       HANDOFF.md                       ⟵ where you stopped, next-chat plan
       METHODOLOGY.md                   ⟵ math + provenance
       SESSION_AUDIT.md                 ⟵ decisions D1–D20 + corrections
       MANUSCRIPT_CHUNKS.md             ⟵ draft paragraphs + figure sketches
       example_panel_D_output.png       ⟵ reference render

     _archive/
       README.md
       heatmap_ordering_upgrade/        ⟵ Module 1 scaffolding
       architecture_classifier/         ⟵ abandoned (do not use)
       original_workflow_diagrams/      ⟵ superseded


   LOAD ORDER FOR NEXT CHAT

       HANDOFF.md  →  METHODOLOGY.md  →  SESSION_AUDIT.md  →  this diagram


════════════════════════════════════════════════════════════════════════════════
```
