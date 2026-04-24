# Workflow Diagram — 4-Encoding Sample Comparison Module

```
══════════════════════════════════════════════════════════════════════════════

              STEP_BP_05  —  FOUR-ENCODING SAMPLE COMPARISON
                        (Diagnostic, not breakpoint caller)

══════════════════════════════════════════════════════════════════════════════

   ┌──────────────────────────────────────────────────────────┐
   │                      INPUTS                               │
   │                                                           │
   │  Required:                                                │
   │    • dosage matrix  (markers × samples, 0-2)              │
   │                                                           │
   │  Optional:                                                │
   │    • sites TSV      (marker, pos — for genomic ordering   │
   │                       and region filtering)                │
   │    • polarity TSV   (marker, l2_sign from STEP29 —        │
   │                       enables the 4th encoding)            │
   │    • groups TSV     (sample, group — from STEP21 u/v      │
   │                       or invgt output)                     │
   └────────────────────────┬─────────────────────────────────┘
                            │
                            ▼
          ╔═════════════════════════════════════════╗
          ║   BUILD FOUR ENCODED MATRICES            ║
          ║   (one preprocessing step)               ║
          ╠═════════════════════════════════════════╣
          ║                                          ║
          ║   1.  MINOR    : X  (identity)           ║
          ║                                          ║
          ║   2.  MAJOR    : 2 - X                    ║
          ║                                          ║
          ║   3.  012      : threshold at GT_THRESH  ║
          ║                  (0 if x<0.5,             ║
          ║                   1 if 0.5≤x<1.5,        ║
          ║                   2 if x≥1.5)             ║
          ║                                          ║
          ║   4.  POLARIZED: flip markers where      ║
          ║                  L2 sign = -1            ║
          ║                  (skipped if no polarity ║
          ║                  file provided)          ║
          ║                                          ║
          ╚═══════════════════╤═════════════════════╝
                              │
        ┌─────────────────────┼────────────────────┬────────────────────┐
        ▼                     ▼                    ▼                    ▼
  ┌───────────┐         ┌───────────┐       ┌───────────┐       ┌───────────┐
  │  MINOR    │         │  MAJOR    │       │   012     │       │ POLARIZED │
  │           │         │           │       │           │       │           │
  │ manhattan │         │ manhattan │       │ manhattan │       │ manhattan │
  │ sample×   │         │ sample×   │       │ sample×   │       │ sample×   │
  │ sample    │         │ sample    │       │ sample    │       │ sample    │
  │ distance  │         │ distance  │       │ distance  │       │ distance  │
  └─────┬─────┘         └─────┬─────┘       └─────┬─────┘       └─────┬─────┘
        │                     │                    │                    │
        ▼                     ▼                    ▼                    ▼
    hclust                hclust               hclust                hclust
    (ward.D2)            (ward.D2)            (ward.D2)            (ward.D2)
    cutree(k=3)          cutree(k=3)          cutree(k=3)          cutree(k=3)
        │                     │                    │                    │
        ▼                     ▼                    ▼                    ▼
   cluster_minor       cluster_major          cluster_012          cluster_polarized
        │                     │                    │                    │
        └─────────┬───────────┴──────────┬─────────┴──────┬────────────┘
                  │                      │                │
                  ▼                      ▼                ▼
           ┌─────────────────────────────────────────────┐
           │   INTER-ENCODING AGREEMENT                   │
           │                                              │
           │   ARI(minor, major)                          │
           │   ARI(minor, 012)                            │
           │   ARI(minor, polarized)                      │
           │   ARI(major, 012)                            │
           │   ARI(major, polarized)                      │
           │   ARI(012,   polarized)                      │
           │                                              │
           │   High ARI everywhere  → structure is        │
           │                          robust to encoding  │
           │                          choice. Signal is   │
           │                          trustworthy.        │
           │                                              │
           │   Low ARI for one pair  → that encoding      │
           │                           diverges. Read     │
           │                           it to see what     │
           │                           polarity matters.  │
           └─────────────────────┬───────────────────────┘
                                 │
                                 ▼
           ┌─────────────────────────────────────────────┐
           │   MDS PROJECTION PER ENCODING                │
           │   (for visual inspection of how samples      │
           │    cluster in 2D under each encoding)        │
           └─────────────────────┬───────────────────────┘
                                 │
                                 ▼
           ┌─────────────────────────────────────────────┐
           │                 OUTPUTS                      │
           │                                              │
           │  {label}_distance_minor.tsv.gz               │
           │  {label}_distance_major.tsv.gz               │
           │  {label}_distance_012.tsv.gz                 │
           │  {label}_distance_polarized.tsv.gz           │
           │                                              │
           │  {label}_clusters.tsv                        │
           │     sample × cluster assignment per encoding │
           │                                              │
           │  {label}_agreement.tsv                       │
           │     ARI for all encoding pairs               │
           │                                              │
           │  {label}_heatmaps_panel.pdf + .png           │
           │     4 heatmaps side by side                  │
           │                                              │
           │  {label}_mds_panel.pdf + .png                │
           │     4 MDS scatter plots side by side         │
           │                                              │
           │  {label}_dendrograms.pdf                     │
           │     4 dendrograms stacked                    │
           └─────────────────────────────────────────────┘


══════════════════════════════════════════════════════════════════════════════

 WHY THIS EXISTS

 Committing to one allele-polarity convention (minor vs major, raw vs
 discretized, unpolarized vs STEP29-polarized) is an arbitrary choice.
 This script refuses to commit. It computes sample distances under all
 four conventions and reports whether the clustering depends on the
 choice.

 Interpretation of the ARI matrix:
   - If all four encodings give ARI ≈ 1 → the biological structure is
     robust, signal is strong, commit to any one encoding going forward.
   - If two encodings agree and two disagree → the disagreeing pair
     captures something the agreeing pair misses (or vice versa).
     Read the MDS panels to see what.
   - If all four disagree → the signal is too weak to trust any one
     clustering. Treat the candidate as composite/ambiguous.

══════════════════════════════════════════════════════════════════════════════

 WHERE THIS SITS IN THE BREAKPOINT WORKFLOW

 BP_05 is a DIAGNOSTIC LAYER, not a breakpoint caller.

   BP_01 → BP_02 → BP_03 → BP_04         (core breakpoint workflow)
                             │
                             └─── BP_05 (optional diagnostic; does NOT
                                          feed into consensus)

 Invoke from the wrapper with environment variable RUN_BP05=1:
   RUN_BP05=1 STEP_Q10_breakpoint_refinement.sh C_gar_LG28 config.R 1

 BP_05 output does not modify the breakpoint call. Use it to decide:
   - Whether the candidate is simple (all 4 encodings agree) or
     composite (encodings disagree).
   - Whether STEP29 polarity harmonization matters here.
   - Whether sub-structure claims (BB_1 vs BB_2) survive across encodings.

══════════════════════════════════════════════════════════════════════════════

 HONEST CAVEAT

 Manhattan distance on dosage is PARTIALLY polarity-invariant by
 construction: a flipped marker contributes |x_i - x_j| to sample-sample
 distance either way. This means for strong 3-band signals (HOMO_1/HET/
 HOMO_2 cleanly separated), all four encodings will typically give ARI ≈ 1
 because the signal is strong enough to survive any polarity scheme.

 The 4-encoding comparison earns its value at SUBTLE structure:
   - Sub-variants within a karyotype (BB_1 vs BB_2)
   - Noisy low-coverage sites where dosage is smeared
   - Recombinant boundaries with mid-sample transitions

 Synthetic validation (performed during development): even with 50% of
 markers polarity-flipped, Manhattan distance on the 3-band signal
 recovered ground truth with ARI = 1.0 across all four encodings. The
 tool is therefore MOST USEFUL for confirming that claimed sub-structure
 is not an encoding artifact — which is exactly the question you care
 about for composite candidates.

══════════════════════════════════════════════════════════════════════════════
```
