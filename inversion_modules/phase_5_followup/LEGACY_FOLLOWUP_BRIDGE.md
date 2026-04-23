# =============================================================================
# LEGACY FOLLOWUP BRIDGE: STEP20-41 → v8.5 INTEGRATION MAP
# =============================================================================
#
# The MODULE_5B followup scripts (STEP20-41) contain mature implementations
# that should be integrated into the v8.5 pipeline. This document maps
# each legacy script to its v8.5 equivalent and integration priority.
#
# Scripts live in: MODULE_5B_Inversion_Followup/current_followup/
# They are NOT deleted — they serve as reference implementations.
#
# =============================================================================
# IMMEDIATE USE (publication figures — use as-is)
# =============================================================================
#
# STEP28  Marker heatmap (ComplexHeatmap, polarity, quality tier)
#   → USE for Panel D genotype heatmap in candidate figures
#   → Input: dosage matrix from C01i core package
#
# STEP35  Harmonized heatmaps (polarity-corrected, multiple modes)
#   → THE publication genotype heatmap with proper polarity correction
#   → Requires: STEP29 polarity output
#
# STEP37  Breakpoint support from SV (724 lines, DELLY2+Manta, Fisher test)
#   → THE breakpoint evidence module
#   → Input: DELLY/Manta VCFs + sample group assignments from C01i
#
# STEP38  Composite figure (537 lines, multi-panel publication layout)
#   → THE manuscript figure generator
#   → Input: all per-candidate outputs
#
# STEP39  Breakpoint visualization (genome-wide + focal views)
#   → Supplementary figure for breakpoint density
#
# =============================================================================
# INTEGRATE INTO v8.5 SCRIPTS
# =============================================================================
#
# STEP19  FST scan (direct from dosage, no vcftools)
#   → Replace C01e's vcftools Fst commands with STEP19's direct computation
#   → Target: STEP_C01e_candidate_figures.R Panel D Fst section
#
# STEP20  Window profile engine (012 strings, Hamming, consensus)
#   → Core functions already reused in C01i/C01j
#   → Remaining: encode_string(), consensus_012() utilities → helpers/
#
# STEP21  GMM state assignment (better than k-means)
#   → Upgrade C01c triangle composition from crude k-means to GMM+rotation
#   → Target: STEP_C01c_triangle_regimes.R sample_composition function
#   → Priority: MEDIUM (k-means works, GMM is cleaner)
#
# STEP22  Subclustering (DBSCAN within stripes, shape classification)
#   → Feeds into C01g tube graph stripe geometry
#   → Target: STEP_C01g_layer4_tube_graph.R node characterization
#
# STEP23  Ancestry association (Cramér's V, chi-square)
#   → Replace C01f_b ancestry_cor with formal statistical tests
#   → Target: STEP_C01f_b_negative_controls.R ancestry instability section
#
# STEP24  Het + marker support (direct from dosage)
#   → Replace C01e's vcftools het commands
#   → Target: STEP_C01e_candidate_figures.R Panel D het section
#
# STEP29  Coherence + polarity
#   → Feed polarity into C01i core package extraction
#   → Polarity = which allele is REF vs INV at each marker
#   → Critical for STEP35 harmonized heatmaps
#   → Target: STEP_C01i_multi_inversion_decomposition.R core_package function
#
# STEP30  Window trajectories (per-sample per-window state tracking)
#   → Complements C01j at sample level (C01j is group level)
#   → Target: addon to STEP_C01j_regime_compatibility_engine.R
#
# STEP32  Multiscale persistence (250/500/1000 SNP windows)
#   → Add as D11 scoring dimension in C01d
#   → Tests if signal survives scale change → more confident
#   → Target: STEP_C01d_candidate_scoring.R new dimension
#
# STEP33  Within-stripe analysis (PCA on HET only, geometry metrics)
#   → Replace C01f T7 carrier substructure test
#   → Target: STEP_C01f_hypothesis_tests.R T7 section
#
# STEP34  Anchor resolution (clean homo as anchor, HET decomposition)
#   → Improve C01g anchor selection
#   → Target: STEP_C01g_layer4_tube_graph.R anchor tracking
#
# STEP36  Clair3 local signatures (phase blocks, weak indels, markers)
#   → Enrich C01i core package + C01h recombinant scanner
#   → Target: STEP_C01i + STEP_C01h
#
# STEP40  Internal coherence (Snake 2+3 across all windows)
#   → Feeds split decisions in C01g tube graph
#   → Target: STEP_C01g_layer4_tube_graph.R split decision logic
#
# STEP41  Membership trajectories (per-sample PC2 sub-clustering)
#   → Detailed complement to C01j stability analysis
#   → Target: addon to STEP_C01j_regime_compatibility_engine.R
#
# =============================================================================
# INTEGRATION PRIORITY
# =============================================================================
#
# P1 (now):  STEP28/35/37/38 — publication figures, use as-is
# P2 (next): STEP29 polarity → C01i, STEP19 Fst → C01e
# P3 (later): STEP21 GMM → C01c, STEP32 multiscale → C01d
# P4 (refine): STEP22/23/33/34/36/40/41 → various upgrades
#
# =============================================================================
# INTEGRATION TABLE (one-row summary)
# =============================================================================
#
# Legacy   Lines  v8.5 Target              Priority  Action
# ------   -----  ----------------------  --------  ------
# STEP19    178   C01e Fst section         P2        Replace vcftools commands
# STEP20    587   C01i/C01j (done)         done      Primitives already shared
# STEP21    265   C01c composition         P3        GMM replaces k-means
# STEP22    475   C01g stripe geometry     P4        Shape classification
# STEP23    309   C01f_b ancestry          P4        Formal stats
# STEP24    238   C01e het section         P2        Direct het from dosage
# STEP25B   246   C01d scoring             merged    Already in candidate_scores
# STEP26    462   C01e figures             P1        More figure types
# STEP27    284   C01e faceted PCA         P3        Unique diagnostic
# STEP28    112   C01e Panel D             P1        USE DIRECTLY
# STEP29    315   C01i polarity            P2        Feed core package
# STEP30    221   C01j trajectories        P3        Per-sample complement
# STEP31    219   C01e diagnostics         merged    Merged with C01e
# STEP32    253   C01d D11 dimension       P3        Multiscale persistence
# STEP33    303   C01f T7 replacement      P4        Stripe geometry
# STEP34    256   C01g anchor              P4        Better anchor selection
# STEP35    270   C01e polarity heatmap    P1        USE DIRECTLY
# STEP36    358   C01i + C01h              P4        Clair3 enrichment
# STEP37    724   Breakpoint module        P1        USE DIRECTLY
# STEP38    537   Composite figure         P1        USE DIRECTLY
# STEP39    603   Breakpoint viz           P1        Supplementary figure
# STEP40    580   C01g coherence           P4        Split decisions
# STEP41    457   C01j stability           P4        Detailed tracking
