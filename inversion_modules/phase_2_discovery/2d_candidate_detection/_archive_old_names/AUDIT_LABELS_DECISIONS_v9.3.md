# ============================================================================
# FULL AUDIT — inv_detect v9.3: Labels, Decisions, Data Flow
# ============================================================================
# 25 files | 6,907 lines | 9 detection phases + 4 peel levels +
# 6 matrix cheats + 1 landscape classifier + 3 plotters
# ============================================================================


# ============================================================================
# 1. STAIRCASE DETECTOR (01) — Vote-based boundary detection
# ============================================================================
#
# INPUT:  sim_mat (N×N)
# OUTPUT: blocks, step_table, rank_table, vote_profile, boundaries
#
# ALGORITHM:
#   For each of N windows, read sim_mat row in both directions.
#   Find positions where similarity drops by ≥ min_drop AND stays
#   dropped for ≥ persist_n consecutive bins. Record positions + heights.
#   Count votes: position X gets +1 for each window that sees a step there.
#   Boundaries = peaks in vote profile above threshold.
#   Blocks = intervals between boundaries with elevated interior sim.
#
# THRESHOLDS (from 00_config.R):
#   STAIR_MIN_BLOCK_WIDTH = 5      # bins (= 250 kb at 50kb/window)
#   STAIR_MIN_DROP        = 0.03   # similarity units
#   STAIR_PERSIST_N       = 3      # bins sustained
#   STAIR_SMOOTH_SPAN     = 3      # running median span
#   vote_blur             = 2      # ±bins for jitter tolerance
#   vote_threshold        = max(min_block_width, 3)  # = 5 votes minimum
#   bg_sim filter         = median(upper_tri) + min_drop
#
# PER-WINDOW OUTPUTS:
#   step_positions[]   — absolute positions of each step
#   step_heights[]     — similarity level BEFORE each step
#   step_drops[]       — similarity level AFTER each step
#   background         — final level after all steps (far-distance sim)
#
# LABELS / FLAGS:
#   family_ld_flag     — TRUE if bg_level > chr_median + 2*MAD
#                        (window sees elevated similarity to distant windows)
#   inside_block       — TRUE if step_position falls inside a detected block
#   containing_block_id — which block contains this internal step
#
# VOTE PROFILE COLUMNS:
#   vote_right    — votes from right-direction steps
#   vote_left     — votes from left-direction steps
#   vote_combined — pmin(vote_right, vote_left)  [boundary from both sides]
#   vote_smooth   — running median smoothed
#   bg_level      — per-window background similarity
#   family_ld     — per-window flag
#
# RANK TABLE:
#   signature     — "540_620_820" = sorted step positions as string
#                   windows with same signature = same structural context
#
# BLOCK TABLE COLUMNS:
#   block_id, start, end, width, start_bp, end_bp, start_mb, end_mb,
#   height (median interior sim), bg_sim, contrast (height - bg_sim),
#   step_down_left, step_down_right, parent_id, n_artifacts
#
# DECISION RULES:
#   step detected       ← v < current_level - min_drop for persist_n bins
#   boundary position   ← vote peak: local max ≥ vote_threshold, min_gap apart
#   block kept          ← interior median > bg_sim + min_drop
#   parent/child        ← block A contains block B by coordinate containment
#
# ISSUES FOUND IN AUDIT:
#   - current_level tracks with 0.95/0.05 EMA which can drift on long
#     gradual declines (family LD). May call false steps on smooth slopes.
#     MITIGATION: persist_n=3 + vote_threshold filters most of these.
#   - vote_blur=2 means a step at position 540 votes for 538-542.
#     Two real boundaries 4 bins apart could merge. Unlikely for real
#     inversions (boundaries are >50 bins apart) but possible for
#     dense micro-inversions.
#   - The 0.95/0.05 EMA constant is not in config. Should be exposed.


# ============================================================================
# 2. NN PERSISTENCE (02) — Independent staircase per NN scale
# ============================================================================
#
# INPUT:  sim_mats at nn0/20/40/80 (named list)
# OUTPUT: per-block NN survival table
#
# LABELS (per-block per-NN):
#   nn<N>_status:
#     "stable"     — block matched at this NN (reciprocal overlap ≥ 0.35)
#     "disappears" — no match at this NN
#     "splits"     — multiple matches (block fragments)
#     "absent"     — no sim_mat available for this NN
#
# SUMMARY LABELS:
#   survives_nn40  — TRUE/FALSE (max_survive_nn ≥ 40)
#   survives_nn80  — TRUE/FALSE
#   nn_topology:
#     "stable"              — all NN scales agree
#     "disappears_nn<N>"    — first disappearance at nnN
#     "splits"              — splits at some scale
#     "complex"             — mixed behavior
#
# THRESHOLD:
#   NN_OVERLAP_THRESH = 0.70 (reciprocal overlap for matching)
#   uses 0.50 * threshold for rough matching = 0.35 effective


# ============================================================================
# 3. FLATNESS (03) — Decay profile shape
# ============================================================================
#
# INPUT:  candidates + sim_mat
# OUTPUT: per-block decay shape metrics
#
# COLUMNS:
#   flatness        — 1 - range/mean (high = flat plateau = inversion)
#   far_near_ratio  — median(far 20%) / median(near 20%) = squareness proxy
#   monotonicity    — fraction of adjacent pairs that are non-increasing
#   sim_q10..q90    — similarity at 5 quantile distances
#   slope_near/mid/far — Theil-Sen slopes in 3 thirds
#   n_sign_flips    — total slope direction changes
#   n_artifact_dips — neg→pos flips (dip then recovery = artifact)
#   has_artifact    — TRUE if any artifact dips
#
# NO CLASSIFICATION LABELS — raw metrics only.


# ============================================================================
# 4. INTERIOR CV (04) — Block homogeneity
# ============================================================================
#
# INPUT:  candidates + sim_mat
# OUTPUT: per-block uniformity
#
# COLUMNS:
#   interior_cv   — sd/mean of upper-triangle values
#   homogeneity   — 1 - cv
#   stripe_count  — local minima in column means below 70% of block mean
#   patchiness    — fraction of values below 50% of block median


# ============================================================================
# 5. GHSL (05) — Sample partition stability
# ============================================================================
#
# INPUT:  candidates + sim_mat
# OUTPUT: per-block partition consistency
#
# COLUMNS:
#   partition_stability     — fraction of adjacent window pairs in same cluster
#   n_consistent_windows    — longest run of same-cluster windows
#   dominant_partition_frac — fraction in largest cluster
#   partition_entropy       — Shannon entropy of partition (bits)


# ============================================================================
# 6. SV OVERLAP (06) — Breakpoint matching
# ============================================================================
#
# INPUT:  candidates + SV file + chr
# OUTPUT: per-block best SV match
#
# COLUMNS:
#   best_sv_id, best_sv_caller, sv_type,
#   sv_overlap_pct, sv_left_dist_kb, sv_right_dist_kb,
#   sv_span_ratio, n_sv_hits
#
# THRESHOLD:
#   SV_MAX_DIST_KB = 500


# ============================================================================
# 7. MATRIX CHEATS (07) — 6 downstream transforms
# ============================================================================
#
# INPUT:  raw sim_mat
# OUTPUT: named list of matrix variants
#
# VARIANT NAMES:
#   raw        — original sim_mat (always present)
#   distcorr   — A1: subtract median sim at each off-diagonal distance
#   localnorm  — A2: local z-score (median/MAD in ±50 bin window)
#   denoised   — A3: iterated 2D median filter (fused lasso proxy)
#   resid_bg   — A4: sim - SVD rank-3 reconstruction (local residual)
#   edge       — A5: Sobel gradient magnitude (boundaries light up)
#   support    — A6: soft threshold at quantile 0.75
#
# THRESHOLDS:
#   CHEAT_LOCALNORM_WINDOW = 50      # ±bins
#   CHEAT_DENOISED_LAMBDA  = 1.0     # proxy iterations ∝ lambda
#   CHEAT_RESIDBG_RANK     = 3       # SVD components for background
#   CHEAT_SUPPORT_QUANTILE = 0.75


# ============================================================================
# 8. BLOC SCORING (08) — Per-block quality metrics
# ============================================================================
#
# INPUT:  blocks + sim_mat (per variant)
# OUTPUT: per-block per-variant scores
#
# COLUMNS:
#   inside_mean, inside_med, outside_mean, contrast,
#   squareness (far_med / near_med), near_med, far_med,
#   sharpness (edge_inside - edge_outside), occupancy,
#   patchiness (cv), shape_class, bg_median, bg_mad
#
# SHAPE CLASSES (from classify_shape):
#   strong_square   — sq>0.7, contrast>0.05, occ>0.5, patch<0.3
#   diffuse_square  — sq>0.5, contrast>0.03, occ>0.3
#   diagonal_band   — sq<0.4, contrast>0.02
#   noise           — contrast<0.02 or occ<0.2
#   ambiguous       — none of the above
#   unknown         — squareness or contrast is NA
#
# THRESHOLDS:
#   BLOC_NEAR_FRAC  = 0.2     # 20% of block width = "near"
#   BLOC_FAR_FRAC   = 0.2     # 20% of block width = "far"
#   BLOC_EDGE_DEPTH = 3       # bins for edge contrast


# ============================================================================
# 9. NN SWEEP TREE (09) — Interval tree with nn_birth
# ============================================================================
#
# INPUT:  sim_mats at multiple NN scales
# OUTPUT: tree of blocks with birth/death NN scales
#
# NODE COLUMNS:
#   node_id, start, end, width, start_mb, end_mb,
#   nn_birth, nn_death, parent_node, children, height, topology
#
# NN_BIRTH CLASSIFICATION:
#   nn_birth ≥ 200   → "INVERSION"
#   nn_birth 80-199  → "CANDIDATE"
#   nn_birth 40-79   → "WEAK_CANDIDATE"
#   nn_birth < 40    → "FAMILY_LD"
#
# TOPOLOGY LABELS:
#   "root"        — first node at coarsest scale
#   "child"       — split from parent at finer scale
#   "novel"       — appeared at finer scale with no parent match
#   "splits"      — this node splits into children at finer scale
#   "disappears"  — this node vanishes at finer scale


# ============================================================================
# 10. CONSENSUS (10) — Cross-variant matching
# ============================================================================
#
# INPUT:  blocks from each matrix variant
# OUTPUT: consensus table
#
# CONFIDENCE LABELS:
#   "HIGH"     — block in ≥ 3 variants (CONSENSUS_MIN_VARIANTS)
#   "MODERATE" — block in 2 variants
#   "LOW"      — block in 1 variant only
#
# RELATIONSHIP LABELS:
#   "universal"         — found in ALL variants
#   "near_universal"    — found in all but 1
#   "hidden_under_haze" — not in raw, found after distcorr/resid_bg
#   "too_patchy_for_raw"— not in raw, found after denoising
#   "raw_only"          — only in raw, not in treated variants (noise?)
#   "partial"           — found in some variants
#
# THRESHOLD:
#   CONSENSUS_MIN_VARIANTS = 3
#   CONSENSUS_OVERLAP      = 0.50


# ============================================================================
# 11. PEELING DIAGNOSTIC (C01n) — Sample removal test
# ============================================================================
#
# INPUT:  precomp RDS (dt with PC_1_<sample>), blocks, optional pruned list
# OUTPUT: before/after comparison per block per peel level
#
# PEEL LEVELS:
#   L1_pruned_list    — remove samples from ngsRelate pruned list
#                       (genome-wide OR per-chromosome table auto-detected)
#   L1b_chrlocal_kin  — remove chr-local PC1 relatives (cor > 0.7)
#                       (computed from precomp, no external file needed)
#   L2_local_coseg    — cluster samples by PC1 trajectory within block,
#                       remove dominant local cluster
#   L3_block_drivers  — remove high PC1 leverage + C01i HOM_INV samples
#   L3_hom_inv_only   — remove only C01i HOM_INV samples
#
# EFFECT CLASSES (from classify_peel_effect):
#   "stable"          — contrast ratio ≥ 0.80
#   "weakened"        — contrast ratio 0.15-0.80
#   "disappeared"     — contrast ratio < 0.15 or after_contrast < 0.01
#   "revealed_child"  — squareness IMPROVED by > 0.10 after peeling
#   "ambiguous"       — before or after contrast is NA
#
# MATH: v_peeled = v_original * mask (mask[fish]=0 for peeled, 1 otherwise)
#        angular distance recomputed from masked loading vectors
#
# THRESHOLDS:
#   max_peel_frac      = 0.3   # never remove >30% of samples
#   local_kin_thresh   = 0.7   # PC1 correlation for L1b
#   local_k            = 5     # k for L2 clustering
#   kin_threshold      = 0.05  # ngsRelate rab threshold for L1


# ============================================================================
# 12. BRIDGE TO CODEBASE (12) — Format adapter
# ============================================================================
#
# INPUT:  detector blocks + scores + NN + peel
# OUTPUT: triangle_intervals.tsv.gz in C01c format
#
# SHAPE → INTERVAL_TYPE MAPPING:
#   strong_square  → "strong_triangle"
#   diffuse_square → "moderate_triangle"
#   diagonal_band  → "sharp_but_not_square"
#   noise          → "weak_zone"
#   ambiguous      → "diffuse_zone"
#   (default)      → "diffuse_zone"
#
# EXTRA COLUMNS PASSED THROUGH (C01d ignores unknown columns):
#   squareness, shape_class, survives_nn40, l1b_peel_effect, source


# ============================================================================
# 13. ANNOTATED PLOT (13) — Chromosome-wide sim_mat with outlines
# ============================================================================
#
# COLOR SYSTEM:
#   Color = inversion SYSTEM (root block). All children inherit root's color.
#   15-color palette: black, blue, red, green, brown, purple, cyan, amber,
#   indigo, rose, teal, yellow-brown, violet, sky, pink.
#   Each independent root block gets the next color.
#
# THICKNESS = DEPTH:
#   depth 0 → 1.3pt solid    depth 3 → 0.55pt dashed
#   depth 1 → 0.95pt solid   depth 4 → 0.40pt dotted
#   depth 2 → 0.70pt dashed  depth 5+ → 0.35pt dotted
#
# LABEL SIZE: 3.2pt at depth 0, decreasing to 1.2pt at depth 6+
#
# MINI-BLOCK DENSITY ZONES:
#   Sliding 2 Mb window. If ≥3 small blocks (<50 bins) in one window →
#   orange dotted rectangle with italic label.
#
# INTERNAL FEATURES: Red × marks on diagonal (shape=4, size=0.8)


# ============================================================================
# 14. LANDSCAPE CLASSIFIER (14) — Full chromosome categorization
# ============================================================================
#
# BLOCK CATEGORIES:
#   strong_inversion    — sq>0.65, cn>0.05, nn40 survives
#   diffuse_inversion   — sq>0.40, cn>0.03
#   complex_system      — n_children ≥ 2
#   nested_fixed        — depth≥1, height>0.85
#   nested_rare         — depth≥1, occupancy<0.3
#   nested_sub_block    — depth≥1, other
#   family_ld_band      — sq<0.35, patch>0.25, peel_l1b="disappeared"
#   diffuse_diagonal    — sq<0.35, patch>0.25, peel_l1b≠"disappeared"
#   diffuse_high_density— n_artifacts>3, sq<0.5
#   weak_signal         — contrast<0.02
#   unclassified        — none of the above
#
# GAP CATEGORIES (regions between blocks):
#   clean_background     — gap_elevation < bg_mad
#   family_ld_band       — far_median > bg_median + bg_mad
#   extended_suppression — gap_elevation > 3 * bg_mad
#   transition_zone      — everything else
#
# CONFIDENCE SCORE (0-100%):
#   Weighted sum of:
#     squareness          × weight 3  (scaled 0-1 direct)
#     contrast            × weight 2  (scaled by /0.15)
#     occupancy           × weight 1  (direct)
#     survives_nn40       × weight 3  (binary 0/1)
#     survives_nn80       × weight 2  (binary 0/1)
#     peel_l1b stability  × weight 3  (stable=1, weakened=0.5, disappeared=0)
#     peel_l2 stability   × weight 2  (stable=1, weakened=0.6, disappeared=0.1)
#     step sharpness      × weight 1  (scaled by /0.1)
#     vote strength       × weight 2  (scaled by /50)
#     low patchiness      × weight 1  (1 - 3*patchiness)
#   Total weights: up to 20 (if all metrics available)
#
# RICH LABEL FORMAT:
#   "B7 Outer Strong Inv [87%]"
#   "24.0-28.0 Mb (4.0 Mb)"
#   "nn40+ peel:OK"


# ============================================================================
# 15. ZOOMED PLOTS (15) — Saturation-driven sub-plots
# ============================================================================
#
# SATURATION CURVE:
#   Order blocks by start position. Accumulate unique covered windows.
#   Plateau = dense packing. Steep rise = new territory.
#
# ZOOM REGION DECISION:
#   Cut when a block adds > gap_frac × N new windows at once.
#   Default gap_frac = 0.05 (5% of chromosome).
#   Each segment between cuts = one zoom plot.
#
# NO FIXED THRESHOLDS for "big vs small" — the saturation shape decides.


# ============================================================================
# DATA FLOW SUMMARY
# ============================================================================
#
# sim_mat (from C01a precomp)
#   │
#   ├─ 01_staircase → blocks + step_table + vote_profile + rank_table
#   │
#   ├─ 07_cheats → 6 treated variants
#   │     │
#   │     └─ 08_bloc_scoring → shape_class per variant
#   │
#   ├─ 02_nn_persistence → survives_nn40, nn_topology
#   │
#   ├─ 09_nn_sweep_tree → nn_birth, classification (INVERSION/CANDIDATE/...)
#   │
#   ├─ 10_consensus → confidence (HIGH/MODERATE/LOW), relationship
#   │
#   ├─ 03-06 evidence → flatness, cv, ghsl, sv columns
#   │
#   ├─ C01n peeling → effect_class per peel level
#   │
#   └─ 14_landscape → final category + confidence_pct + rich_label
#         │
#         ├─ 13_plot → chromosome-wide annotated heatmap
#         │
#         ├─ 15_plot_zoom → saturation curve + per-region zoomed plots
#         │
#         └─ 12_bridge → triangle_intervals.tsv.gz → C01d scoring
#
#
# ============================================================================
# KNOWN GAPS / ISSUES
# ============================================================================
#
# 1. EMA CONSTANT (0.95/0.05) in staircase traversal is hardcoded.
#    Should be in config. On long gradual family LD declines, the EMA
#    tracks downward and may call false steps.
#
# 2. LANDSCAPE CLASSIFIER uses fixed thresholds (sq>0.65 for strong,
#    sq>0.40 for diffuse, etc.). These were set without calibration
#    against real data. First run on LG01 will reveal if they're right.
#
# 3. CONFIDENCE SCORE scales contrast by /0.15 and votes by /50.
#    These divisors are arbitrary. Need calibration: what's the typical
#    contrast of a known inversion? What's the typical vote count?
#
# 4. PEEL DIAGNOSTIC uses local_kin_thresh=0.7 for L1b.
#    Too high = misses real local relatives. Too low = removes too many.
#    The per-chr ngsRelate (11a/11b) provides a more rigorous alternative.
#
# 5. VOTE BLUR = 2 bins. Two real boundaries 4 bins apart would merge.
#    For most inversions (>50 bins wide) this is fine. For micro-inversions
#    (<10 bins) it could be a problem.
#
# 6. MATRIX CHEATS are SLOW on 9000×9000 matrices. A2 (local contrast)
#    is O(N² × window²). Use --skip-cheats for initial testing.
#
# 7. NO UNIT TESTS. The staircase, voting, and classification logic has
#    no automated test suite. Should build synthetic sim_mats with known
#    blocks and verify detection.
#
# 8. STEP TABLE can be HUGE (N windows × max_steps × 2 directions).
#    For N=9000 and max_steps=8: up to 144,000 rows. Saved as .tsv.gz
#    but still large. Consider summarizing for production.
#
# 9. SATURATION-BASED ZOOM (15) uses gap_frac=0.05. On a chromosome
#    where 90% is covered by blocks, this might produce only 1 zoom
#    covering the whole chromosome (useless). Need a fallback that
#    splits large zoom regions at the highest internal boundaries.
#
# 10. BRIDGE (12) maps shape_class → interval_type but C01d also expects
#     sample_composition, bridges, offdiag_linkage files. These are
#     written as empty data.tables. C01d handles empty gracefully but
#     some scoring dimensions (D3-D5, D7-D8) will be NA.
#
# ============================================================================
