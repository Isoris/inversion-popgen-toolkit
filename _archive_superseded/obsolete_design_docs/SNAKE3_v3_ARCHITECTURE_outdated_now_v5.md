# =============================================================================
# SNAKE 3 v3 — GHSL HAPLOTYPE CONTRAST ARCHITECTURE
# =============================================================================
# Draft: 2026-04-07
# Status: DESIGN (not yet coded)
#
# This document defines the complete redesign of Snake 3 from a per-sample
# heterozygosity counter (v2) to a true group-contrast haplotype sharing
# layer that can:
#   1. Confirm inversions independently from eigenvalue/sim_mat methods
#   2. Resolve the het (middle) band that eigenvalues can't separate from noise
#   3. Produce GHSL heatmaps showing inversion block structure
#   4. Compare divergent karyotypes at sample×sample resolution
#   5. Integrate with Snake 1 candidates for scoring
#
# =============================================================================
# TABLE OF CONTENTS
# =============================================================================
#
#   1. CORE CONCEPT
#   2. DATA FLOW
#   3. INPUT FORMAT
#   4. PROCESSING STEPS (5 stages)
#   5. OUTPUT FORMAT
#   6. DIAGNOSTIC FIGURES
#   7. INTEGRATION WITH PIPELINE
#   8. IMPLEMENTATION PLAN
#
# =============================================================================
# 1. CORE CONCEPT
# =============================================================================
#
# WHAT GHSL MEASURES:
#
#   For each local-PCA window (100 SNPs), we ask:
#   "Do samples in the same PC1 group share the same phased haplotype
#    configuration, while samples in different groups have different ones?"
#
#   If YES → strong GHSL contrast → real inversion (groups carry different
#            haplotype backgrounds across the region)
#   If NO  → low GHSL contrast → family noise or random clustering
#
# WHY THIS WORKS IN HATCHERY:
#
#   Ne~20 means founder haplotypes are long (low recombination over few
#   generations). If an inversion captured a specific haplotype background,
#   every carrier shares that exact haplotype across the entire inverted
#   region. Family structure HELPS here because it preserves the original
#   inversion-linked haplotype intact.
#
# THE PHASING LIMITATION:
#
#   At ~9x Illumina, WhatsHap produces phase blocks of ~30-800 bp (short).
#   But our windows are 100 SNPs = ~20-100 kb. So a window may contain
#   10-50 independent phase blocks per sample. Within each block, phase
#   is reliable. Across blocks, strand assignment is arbitrary.
#
# THE SOLUTION — PHASE-BLOCK-AWARE CONTRAST:
#
#   We don't stitch phase blocks. Instead:
#   (a) Within each phase block that spans ≥2 het SNPs in the window,
#       extract the haplotype pattern (sequence of 0|1 vs 1|0 calls)
#   (b) Compare patterns WITHIN phase blocks across samples
#   (c) Two samples sharing an inversion allele will have matching (or
#       perfectly flipped) patterns within overlapping phase blocks
#   (d) Two samples with different alleles will have uncorrelated patterns
#
#   The flip ambiguity (hapA/hapB labeling) is handled by trying both
#   orientations and taking the better concordance.
#
# HOMOZYGOUS SAMPLES:
#
#   Hom-REF and hom-ALT samples are fully informative for GHSL:
#   - They have NO het sites (no phase ambiguity)
#   - Their genotype is unambiguous: 0/0 or 1/1
#   - They contribute the clearest signal to their group
#   - For haplotype comparison: hom samples provide the "anchor"
#     against which het samples are evaluated
#
#   In practice: hom samples define what the inversion haplotype looks
#   like. Het samples should show a mix. This is the key to resolving
#   the middle band.
#
# =============================================================================
# 2. DATA FLOW
# =============================================================================
#
#   INPUTS:
#     ┌─────────────────────────┐    ┌──────────────────────────────────┐
#     │ precomp/<chr>.precomp.rds │    │ postprocess_results/<chr>/<sample>│
#     │  • window coordinates    │    │  • all_variants_with_phase.tsv   │
#     │  • PC1 loadings per sample│    │  • phase_blocks.tsv              │
#     │  • k=3 band assignments  │    │                                  │
#     └──────────┬──────────────┘    └──────────────┬───────────────────┘
#                │                                   │
#                ▼                                   ▼
#     ┌──────────────────────────────────────────────────────────┐
#     │  STAGE 1: EXTRACT — per-sample × per-window haplotypes  │
#     └──────────────────────┬───────────────────────────────────┘
#                            │
#                            ▼
#     ┌──────────────────────────────────────────────────────────┐
#     │  STAGE 2: COMPARE — pairwise haplotype concordance      │
#     │           within overlapping phase blocks                │
#     └──────────────────────┬───────────────────────────────────┘
#                            │
#                            ▼
#     ┌──────────────────────────────────────────────────────────┐
#     │  STAGE 3: CONTRAST — within-group vs between-group      │
#     │           using PC1 k=3 band assignments                 │
#     └──────────────────────┬───────────────────────────────────┘
#                            │
#                            ▼
#     ┌──────────────────────────────────────────────────────────┐
#     │  STAGE 4: SCORE — per-window GHSL contrast score        │
#     │           + per-sample group confidence                  │
#     └──────────────────────┬───────────────────────────────────┘
#                            │
#                            ▼
#     ┌──────────────────────────────────────────────────────────┐
#     │  STAGE 5: VISUALIZE — heatmaps, strip plots, figures    │
#     └─────────────────────────────────────────────────────────┘
#
# =============================================================================
# 3. INPUT FORMAT
# =============================================================================
#
# FROM PRECOMP (per chromosome):
#   precomp/<chr>.precomp.rds contains:
#     $dt — data.table with columns:
#       global_window_id, start_bp, end_bp, chrom
#       PC_1_<sample> — PC1 loading per sample per window
#     $sim_mat — window×window similarity matrix
#     $mds_mat — MDS coordinates
#
#   The PC1 columns give us the group assignments via k=3 clustering
#   (same as used in Snake 1 cores and het_contrast).
#
# FROM CLAIR3 POSTPROCESS (per sample × per chromosome):
#   all_variants_with_phase.tsv — key columns:
#     CHROM, POS, REF, ALT1
#     GT_CLASS        — "het", "hom_ref", "hom_alt"
#     IS_SIMPLE_BIALLELIC — TRUE/FALSE
#     IS_PHASED       — TRUE/FALSE
#     PHASE_GT        — "0|1", "1|0" (for phased hets)
#     PS              — phase set start position
#     PHASE_BLOCK_ID  — block identifier
#     PHASE_TIER      — "TIER_1_WHATSHAP" or "TIER_2_READPAIR"
#     IS_SNP          — TRUE/FALSE
#     QUAL, GQ, DP    — quality metrics
#
# =============================================================================
# 4. PROCESSING STEPS
# =============================================================================
#
# ─── STAGE 1: EXTRACT ───────────────────────────────────────────────
#
# For each chromosome, for each window, for each sample:
#
#   1a. Load all_variants_with_phase.tsv for the sample
#   1b. Filter to: IS_SIMPLE_BIALLELIC=TRUE, IS_SNP=TRUE, QUAL>=20, GQ>=10
#   1c. For variants in [window_start, window_end]:
#       - Separate into HET (phased) and HOM sites
#       - For HET+phased: group by PHASE_BLOCK_ID
#       - For each phase block with ≥2 het SNPs in the window:
#           Extract the haplotype pattern as a binary vector
#           (0 = "0|1", 1 = "1|0")
#       - For HOM sites: record genotype (0=hom_ref, 2=hom_alt)
#
#   Output: per-sample-window structure:
#     {
#       sample_id, window_id,
#       n_het_phased,        # phased het SNPs in window
#       n_hom,               # hom SNPs in window
#       n_phase_blocks,      # phase blocks covering the window
#       phase_block_patterns: [  # list of patterns per block
#         { block_id, positions: [pos1, pos2, ...],
#           hap_pattern: [0,1,1,0,...] }
#       ],
#       hom_genotypes: { pos1: 0, pos2: 2, ... }
#     }
#
# PERFORMANCE NOTE:
#   Loading 226 per-sample TSV files per chromosome is expensive.
#   Strategy: precompute a merged per-chr file during prep stage
#   (one-time cost), then the R script reads a single file per chr.
#
# ─── STAGE 2: COMPARE ──────────────────────────────────────────────
#
# For each window, compute pairwise haplotype concordance between
# all sample pairs. Two comparison modes:
#
# MODE A — Phase-block concordance (for het×het pairs):
#   For samples i and j, find phase blocks that overlap in both samples.
#   Within each overlapping block region:
#     - Extract the haplotype patterns for both samples
#     - At positions where BOTH samples have phased hets:
#       Compare: same phase (0|1 vs 0|1) or flipped (0|1 vs 1|0)
#     - Try both orientations (flip one sample's pattern)
#     - Concordance = max(fraction_same, fraction_flipped)
#       If concordance > 0.5 → they share the same haplotype background
#       (with known orientation)
#     - Average across all overlapping blocks
#
# MODE B — Genotype concordance (for hom×hom and hom×het pairs):
#   At positions where both samples have called genotypes:
#     - IBS (identity by state): fraction of matching alleles
#     - This works for hom×hom (clear), hom×het (informative),
#       and unphased het×het (less informative but still usable)
#
# COMBINED CONCORDANCE:
#   concordance(i,j) = w_phase * phase_concordance + w_geno * geno_concordance
#   where:
#     w_phase = n_shared_phased_positions / n_total_shared_positions
#     w_geno  = 1 - w_phase
#   This naturally upweights phase information when available and
#   falls back to genotype IBS when not.
#
#   Output: per-window concordance matrix (226×226, symmetric)
#
# SPARSE OPTIMIZATION:
#   Don't compute all 226×226 = 25,538 pairs per window.
#   Only compute pairs where both samples have ≥ MIN_SHARED_SITES
#   shared positions (e.g., ≥5). Store as sparse triplets.
#
# ─── STAGE 3: CONTRAST ─────────────────────────────────────────────
#
# For each window, use the PC1-based k=3 group assignments:
#   band1 (low PC1), band2 (mid PC1), band3 (high PC1)
#
# Compute:
#   within_sim  = mean concordance for pairs WITHIN the same band
#   between_sim = mean concordance for pairs in DIFFERENT bands
#
#   ghsl_contrast = within_sim - between_sim
#     Range: [-1, +1]
#     Positive = groups have distinct haplotypes (inversion signal)
#     Zero = no difference (noise)
#     Negative = shouldn't happen (means grouping is wrong)
#
# ADDITIONAL METRICS:
#   ghsl_ratio = within_sim / (between_sim + epsilon)
#     Ratio form, >1 = inversion signal
#
#   band_separation[1v2] = mean_concordance(band1×band1) - mean_concordance(band1×band2)
#   band_separation[1v3] = same for band1 vs band3
#   band_separation[2v3] = same for band2 vs band3
#     These tell you WHICH pair of bands is most distinct
#
#   het_band_identity:
#     For each sample in band2 (the het/middle band):
#       mean concordance with band1 members
#       mean concordance with band3 members
#     If het sample has ~equal concordance with both →
#       it's truly heterozygous for the inversion
#     If het sample matches one band strongly →
#       misassigned, should be in that band
#
# ─── STAGE 4: SCORE ────────────────────────────────────────────────
#
# Per-window GHSL score:
#   ghsl_score = ghsl_contrast (from stage 3)
#   ghsl_quality = based on:
#     - n_pairs_computed (coverage of the concordance matrix)
#     - n_samples_with_phasing (how many samples contributed)
#     - mean_shared_sites (how many positions per comparison)
#   ghsl_status = PASS / WEAK / FAIL based on adaptive thresholds
#
# Per-sample confidence:
#   For each sample in each window:
#     within_mean  = mean concordance with own-band members
#     between_mean = mean concordance with other-band members
#     sample_ghsl  = within_mean - between_mean
#     assignment_confidence = STRONG / MODERATE / AMBIGUOUS
#       STRONG: sample_ghsl > 2*median_sample_ghsl
#       MODERATE: sample_ghsl > 0
#       AMBIGUOUS: sample_ghsl ≤ 0 (might be misassigned)
#
# Per-window band verification:
#   Do the GHSL-derived groups agree with PC1-derived groups?
#   ghsl_band_agreement: fraction of samples where GHSL confirms PC1 band
#   If low → the eigenvalue grouping might be wrong at this window
#
# ─── STAGE 5: VISUALIZE ────────────────────────────────────────────
#
# See section 6 for full figure catalog.
#
# =============================================================================
# 5. OUTPUT FORMAT
# =============================================================================
#
# CORE OUTPUTS (per chromosome, TSV.gz):
#
#   snake3v3_window_track_<chr>.tsv.gz
#     Per-window summary:
#       chrom, global_window_id, start_bp, end_bp, pos_mb
#       ghsl_contrast          — within minus between similarity
#       ghsl_ratio             — within / between
#       band_sep_1v2           — band1 vs band2 contrast
#       band_sep_1v3           — band1 vs band3 contrast
#       band_sep_2v3           — band2 vs band3 contrast
#       n_pairs_total          — total pairs computed
#       n_pairs_within         — within-group pairs
#       n_pairs_between        — between-group pairs
#       n_samples_with_data    — samples with ≥MIN_SHARED phasing
#       mean_shared_sites      — mean shared positions per pair
#       mean_within_sim        — mean within-group concordance
#       mean_between_sim       — mean between-group concordance
#       ghsl_quality           — GOOD / MODERATE / LOW
#       ghsl_status            — PASS / WEAK / FAIL
#       phase_coverage         — fraction of samples with phasing
#       snake3_pa              — 1 if PASS, 0 otherwise (for PA matrix)
#
#   snake3v3_sample_scores_<chr>.tsv.gz
#     Per-sample × per-window detail:
#       sample_id, global_window_id
#       pc1_band               — PC1-derived band (1/2/3)
#       n_het_phased           — phased het SNPs in window
#       n_hom                  — hom SNPs in window
#       n_phase_blocks         — phase blocks in window
#       within_mean            — mean concordance with own band
#       between_mean           — mean concordance with other bands
#       sample_ghsl            — within_mean - between_mean
#       assignment_confidence  — STRONG / MODERATE / AMBIGUOUS
#       best_match_band        — which band this sample matches best
#       ghsl_confirms_pc1      — TRUE if best_match == pc1_band
#
#   snake3v3_concordance_<chr>.tsv.gz (OPTIONAL — large)
#     Full concordance matrix per window (sparse format):
#       global_window_id, sample_i, sample_j, concordance,
#       n_shared_sites, n_phase_shared, n_geno_shared
#     Only stored for windows with ghsl_status == PASS
#     This enables downstream heatmap generation without recomputation.
#
#   snake3v3_het_resolution_<chr>.tsv.gz
#     Middle-band analysis:
#       sample_id, global_window_id
#       concordance_with_band1 — mean concordance with band1 members
#       concordance_with_band3 — mean concordance with band3 members
#       het_score              — |conc_band1 - conc_band3| (low = true het)
#       inferred_genotype      — INV_HOM_REF / INV_HET / INV_HOM_ALT
#       confidence             — from het_score magnitude
#
# SUMMARY OUTPUTS:
#
#   snake3v3_summary.tsv
#     Per-chromosome: n_windows, n_pass, n_weak, n_fail,
#     mean_contrast, max_contrast, phase_coverage
#
#   snake3v3_window_pa.tsv.gz
#     Roary-style PA matrix for integration with C01d scoring:
#       chrom, global_window_id, snake3_pa (0/1), ghsl_contrast
#
# =============================================================================
# 6. DIAGNOSTIC FIGURES
# =============================================================================
#
# ─── FIGURE G1: GHSL Contrast Profile ──────────────────────────────
#   X-axis: genomic position (Mb)
#   Y-axis: ghsl_contrast score
#   Color: PASS (green), WEAK (orange), FAIL (grey)
#   Overlay: Snake 1 candidate regions as shaded rectangles
#   Purpose: see where GHSL signal aligns with eigenvalue candidates
#
# ─── FIGURE G2: GHSL Heatmap (the key figure) ─────────────────────
#   Full sample×sample concordance heatmap for a candidate region.
#   Rows/columns: 226 samples, ordered by PC1 band assignment
#   Fill: concordance (0=blue, 1=red)
#   Side annotation: band (1=blue, 2=green, 3=red), ancestry (K=8 colors)
#   Block structure: if inversion is real, you see 3×3 block pattern:
#     ┌─────────┬─────────┬─────────┐
#     │  HIGH   │   MED   │   LOW   │  ← band1×band1 = high concordance
#     ├─────────┼─────────┼─────────┤
#     │   MED   │  MIXED  │   MED   │  ← band2 = intermediate (het)
#     ├─────────┼─────────┼─────────┤
#     │   LOW   │   MED   │  HIGH   │  ← band3×band3 = high concordance
#     └─────────┴─────────┴─────────┘
#   The off-diagonal blocks (band1×band3) should be LOW (different haplotypes).
#   The band2 row/column should show INTERMEDIATE values with both groups.
#
#   This is the publication-quality figure that proves an inversion exists.
#
#   Multiple versions:
#     G2a: average concordance across the entire candidate region
#     G2b: per-window heatmaps (small multiples, 4-6 windows from the region)
#     G2c: delta heatmap (region concordance minus background concordance)
#
# ─── FIGURE G3: Het Band Resolution ───────────────────────────────
#   X-axis: concordance with band1
#   Y-axis: concordance with band3
#   Points: one per sample (colored by PC1 band assignment)
#   Expected pattern for real inversion:
#     - Band1 samples cluster at (high, low) — match band1, differ from band3
#     - Band3 samples cluster at (low, high) — match band3, differ from band1
#     - Band2 (het) samples cluster at (medium, medium) — partial match to both
#   For family noise:
#     - All samples scatter without clear separation
#
# ─── FIGURE G4: Phase Coverage Heatmap ─────────────────────────────
#   X-axis: window position (Mb)
#   Y-axis: 226 samples
#   Fill: fraction of window with phased data (0=white, 1=dark blue)
#   Purpose: shows where Snake 3 has enough data to be reliable
#   Also shows if certain samples consistently lack phasing (low coverage)
#
# ─── FIGURE G5: GHSL vs Eigenvalue Comparison ──────────────────────
#   4-panel figure per chromosome:
#     Panel A: inv-likeness profile (from C01a)
#     Panel B: ghsl_contrast profile (from Snake 3)
#     Panel C: het_contrast profile (from C01a v8.5)
#     Panel D: concordance (correlation between A, B, C tracks)
#   Purpose: see where the three independent signals agree/disagree
#
# ─── FIGURE G6: Pairwise Karyotype Divergence ─────────────────────
#   For pairs of samples with known different karyotypes (from PC1 bands):
#   X-axis: genomic position along chromosome
#   Y-axis: GHSL concordance between the pair
#   Multiple pairs overlaid (thin lines) with mean (thick line)
#   Expected: concordance drops sharply at inversion boundaries
#   This gives you breakpoint resolution from haplotype data
#
# ─── FIGURE G7: Per-Candidate Composite ───────────────────────────
#   Multi-panel publication figure for each confirmed candidate:
#     Top: eigenvalue spectrum + inv-likeness
#     Middle: GHSL heatmap (G2a style)
#     Bottom-left: het resolution scatter (G3)
#     Bottom-right: pairwise divergence profile (G6)
#   This is THE manuscript figure for each inversion.
#
# ─── FIGURE G8: Genome-wide GHSL Summary Strip ────────────────────
#   All 28 chromosomes, ideogram-style strip:
#   Color: ghsl_contrast (grey=low, red=high)
#   Overlay: candidate regions from Snake 1
#   Purpose: genome-wide view of where GHSL confirms vs contradicts
#
# =============================================================================
# 7. INTEGRATION WITH PIPELINE
# =============================================================================
#
# WHERE SNAKE 3 v3 FITS:
#
#   PHASE_01: LANDSCAPE
#     01A precompute ← provides PC1 bands + window grid
#
#   PHASE_02: DETECTION
#     02A cores ← independent
#     02B merge ← independent
#     02C triangles ← independent
#     02D Snake 3 v3 ← NEW, runs in parallel with 02A-02C
#                       only needs precomp + Clair3 postprocess
#
#   PHASE_03: SCORING
#     03A scoring ← Snake 3 PA matrix added as new dimension
#       D11_ghsl_contrast: ghsl_contrast score (0-1 rescaled)
#       pa_snake3_frac: fraction of windows with PASS GHSL
#
#   PHASE_04: CHARACTERIZATION
#     04A decomposition ← Snake 3 concordance used for group validation
#     04G het resolution ← NEW, uses Snake 3 het_resolution output
#       Replaces crude PC1 middle-band assignment with GHSL-verified
#       inversion genotypes (INV_HOM_REF / INV_HET / INV_HOM_ALT)
#
#   PHASE_05: FIGURES
#     05B candidate figures ← G7 composite added per candidate
#     05D all-layers strip ← Snake 3 track added as layer
#
# SCORING INTEGRATION (C01d):
#   Snake 3 provides two things to the 10D scoring:
#     1. ghsl_contrast as dimension D11 (11th scoring dimension)
#     2. PA matrix (snake3_pa per window) for cross-layer agreement
#   A candidate with pa_layers_agree = 5 (core + merge + triangle +
#   regime + GHSL all agree) is virtually certain to be real.
#
# HET RESOLUTION INTEGRATION (C01i/C01j):
#   Snake 3's het_resolution output provides GHSL-verified genotype
#   assignments that replace the crude PC1-based band assignments.
#   This is critical for:
#     - Accurate carrier frequency estimation
#     - Correct Fst computation (need to know who is het)
#     - Recombinant detection (need to identify het breakpoints)
#     - Founder pack analysis (MODULE_6 needs clean genotypes)
#
# =============================================================================
# 8. IMPLEMENTATION PLAN
# =============================================================================
#
# SCRIPTS TO CREATE:
#
#   snakes/STEP_C04_snake3_ghsl_v3.R          (main engine, ~800-1000 lines)
#     Stages 1-4: extract, compare, contrast, score
#     Reads: precomp RDS + Clair3 postprocess TSV files
#     Writes: all core outputs (track, sample_scores, het_resolution, PA)
#
#   snakes/STEP_C04b_snake3_ghsl_figures.R    (figures, ~500-700 lines)
#     Stage 5: all diagnostic figures G1-G8
#     Reads: Snake 3 outputs + precomp + candidate regions
#
#   utils/prepare_ghsl_merged_variants.sh     (prep, ~100 lines)
#     One-time: merge per-sample postprocess TSVs into a single
#     per-chromosome file for fast loading in R
#     Columns: sample_id, pos, gt_class, phase_gt, ps, phase_block_id,
#              phase_tier, is_phased, ref, alt
#     Filter: IS_SIMPLE_BIALLELIC=TRUE, IS_SNP=TRUE
#     Output: ghsl_prep/<chr>.merged_phased_snps.tsv.gz
#
#   launchers/LAUNCH_C04_snake3_ghsl.slurm    (SLURM launcher)
#
# PREP SCRIPT DETAILS:
#   The prep merges 226 × all_variants_with_phase.tsv into one file per chr.
#   This is a one-time cost (~10 min per chr) that avoids loading 226
#   separate files in R. The merged file is ~500MB per chr compressed.
#
#   Alternatively: prep as SLURM array (1 task per chr, 28 tasks).
#
# MAIN ENGINE PSEUDOCODE:
#
#   for each chromosome:
#     load precomp RDS → get window grid + PC1 columns
#     compute k=3 bands from PC1 (same as het_contrast)
#     load merged_phased_snps.tsv.gz for this chr
#
#     for each window:
#       # Stage 1: Extract
#       subset variants to [window_start, window_end]
#       for each sample:
#         extract phased het patterns per phase block
#         extract hom genotypes
#
#       # Stage 2: Compare
#       for each sample pair (i,j) with ≥MIN_SHARED sites:
#         compute phase_concordance (within shared phase blocks)
#         compute geno_concordance (IBS at shared sites)
#         combined = weighted average
#
#       # Stage 3: Contrast
#       split pairs into within-band and between-band
#       ghsl_contrast = mean(within) - mean(between)
#       compute band-specific separations
#
#       # Stage 4: Score
#       per-sample: within_mean, between_mean, confidence
#       het resolution: band2 members → concordance with band1 vs band3
#
# COMPUTATIONAL COMPLEXITY:
#   Per window: O(N^2) pairwise comparisons where N=226
#   = 25,425 pairs per window
#   Per chromosome: ~5000-10000 windows × 25k pairs = ~125M-250M comparisons
#   This is expensive. Mitigations:
#     1. Sparse: only compute pairs with ≥5 shared sites (~50-80% of pairs)
#     2. Vectorized: data.table operations, not nested loops
#     3. Parallel: mclapply across windows within a chromosome
#     4. Pre-filter: skip windows where <20 samples have phasing data
#
#   Estimated runtime: 30-60 min per chromosome on 28 cores.
#   Total: ~14-28 hours for all 28 chromosomes (SLURM array recommended).
#
# MEMORY:
#   The merged phased SNP file + precomp for one chr: ~2-4 GB
#   The concordance matrix per window (sparse): negligible
#   Peak memory: ~16-32 GB per chromosome.
#
# =============================================================================
# 9. KEY DESIGN DECISIONS
# =============================================================================
#
# Q: Should we use the same k=3 bands as Snake 1, or re-cluster?
# A: Use the same PC1 k=3 bands from precomp. This ensures consistency
#    across all layers. GHSL's job is to VALIDATE the grouping, not
#    create a new one.
#
# Q: What if a window has <3 groups (e.g., no hom-alt carriers)?
# A: Fall back to k=2. If the smallest group has <5% of samples,
#    use k=2 (split at PC1 median). Report which k was used.
#    This handles rare inversions where only 2 groups exist.
#
# Q: Should we store the full concordance matrix?
# A: Only for PASS windows in candidate regions. The full matrix
#    is needed for figure G2 (heatmap) but too large to store
#    genome-wide. For FAIL windows, only store the summary stats.
#
# Q: How to handle TIER_2 (read-pair) phasing vs TIER_1 (WhatsHap)?
# A: Compute both: ghsl_t1_only (conservative) and ghsl_combined
#    (includes TIER_2). Report both. TIER_2 adds coverage but may
#    be noisier. Let the user decide which to trust.
#
# Q: Should Snake 3 run per-chr or genome-wide?
# A: Per-chr (consistent with all other snakes). Each chromosome
#    gets its own adaptive thresholds based on its own phase coverage
#    and baseline concordance distribution.
#
# Q: What about the existing Snake 3 v2 code?
# A: Keep as reference (STEP10h_snake3_ghsl_haplotype_contrast_v2.R).
#    v3 is a clean rewrite with different logic (group contrast vs
#    individual het rate). The prep script can be shared.
#
# =============================================================================
# 10. PREP FILE SPECIFICATION
# =============================================================================
#
# The prep step produces one file per chromosome:
#   ghsl_prep/<chr>.merged_phased_snps.tsv.gz
#
# Columns:
#   sample_id          — e.g., CGA407
#   pos                — variant position (1-based)
#   ref                — reference allele
#   alt                — alternate allele
#   gt_class           — het / hom_ref / hom_alt
#   is_phased          — TRUE / FALSE
#   phase_gt           — "0|1" / "1|0" / NA (for hom or unphased)
#   ps                 — phase set position (for phased hets)
#   phase_block_id     — block ID (for phased hets)
#   phase_tier         — TIER_1_WHATSHAP / TIER_2_READPAIR / NA
#   qual               — variant quality
#   gq                 — genotype quality
#   dp                 — read depth
#
# Filter applied during prep:
#   IS_SIMPLE_BIALLELIC = TRUE
#   IS_SNP = TRUE (biallelic SNPs only for v3; indels later)
#   FILTER = "PASS" or STATUS includes pass/rescue
#   Dedup: one row per sample × position
#
# Sort: by pos, then sample_id (enables efficient windowed access)
#
# Expected size: ~200-500 MB per chromosome (compressed)
# Expected rows: ~50M-150M per chromosome (226 samples × 200k-700k SNPs)
#
# =============================================================================
