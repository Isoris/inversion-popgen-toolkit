// shared/band_tracking/index.js
// =====================================================================
// Public API for the band-tracking layer.
//
// Pipeline (top-down):
//
//   per-window K-means labels  (existing — produced by shared/per_l2_cluster.js
//                                or upstream R-side)
//          ↓
//   single_band_track_from_seed  →  one-band trajectory across windows
//   single_band_score_continuity →  per-step retained/lost/gained + jaccard
//          ↓
//   het_detect_candidate_band    →  intermediate-PC1 bands per window
//   het_track_skeleton           →  forward+backward stitched het track
//   het_define_interval          →  bp coordinates from a skeleton
//          ↓
//   hom_anchor_to_het            →  HOM_A / HOM_B sample sets for the interval
//          ↓
//   iv_call_samples_from_skeleton →  per-sample karyotype calls
//          ↓
//   iv_merge_het_tracks          →  optional: merge nearby intervals with
//                                   shared sample cores (post-call cleanup
//                                   for noisy stretches)
//
// All math reuses the existing shared/contingency.js primitives —
// nothing here rewrites the L3 contingency core.
// =====================================================================

export {
  bandMembers, bandJaccard,
  single_band_score_continuity,
  single_band_track_from_seed,
} from './single_band.js';

export {
  meanPc1PerBand,
  het_detect_candidate_band,
  het_track_skeleton,
  het_define_interval,
  iv_merge_het_tracks,
} from './het.js';

export {
  hom_anchor_in_window,
  hom_anchor_to_het,
} from './hom.js';

export {
  iv_call_samples_from_skeleton,
} from './iv.js';

// ---------------------------------------------------------------------
// LAYER 1 — Trajectory similarity (fast macro-grouping aid)
// ---------------------------------------------------------------------
export {
  pickPc1OrientationReferenceSamples,
  computePc1SignAnchors,
  band_compute_pc1_trajectory,
  band_pairwise_trajectory_correlation,
  band_group_by_trajectory_similarity,
} from './trajectory.js';

// ---------------------------------------------------------------------
// LAYER 2 — BandSet Projection (set-based authority)
// ---------------------------------------------------------------------
export {
  PATTERN_CLASS,
  bp_compute_projection_vector,
  bp_detect_visited_excluded_bands,
  bp_classify_projection_pattern,
  bp_project_bandset_to_target_bands,
} from './projection.js';

// ---------------------------------------------------------------------
// COMBINER — Layer 1 + Layer 2 → macro-band groups → karyotype model
// ---------------------------------------------------------------------
export {
  kt_combine_trajectory_and_projection_evidence,
  kt_infer_macro_band_groups,
  kt_resolve_karyotype_model,
} from './karyotype_model.js';

// ---------------------------------------------------------------------
// VOTE EVIDENCE — raw vote extraction + co-association matrix
// (shared foundation for both band-centric and locus-centric views)
// ---------------------------------------------------------------------
export {
  extract_votes,
  build_coassociation_matrix,
  build_per_band_vote_index,
  voteRecords_from_projections,
} from './vote_evidence.js';

// ---------------------------------------------------------------------
// LAYER A — band-centric: "standing on the band, looking at who votes for us"
// ---------------------------------------------------------------------
export {
  collect_voters_for_band,
  compute_partner_affinities,
  derive_partner_sets,
  compute_voter_consensus,
  compute_overlap_conflict,
  build_per_band_view,
} from './band_voters.js';

// ---------------------------------------------------------------------
// LAYER B math — bruteforce partition enumeration + scoring
// ---------------------------------------------------------------------
export {
  enumerate_partitions_as_blocks,
  score_partition_against_coassoc,
  enumerate_and_score_all_partitions,
  select_top_partitions_adaptive,
  partitions_are_compatible,
  bellNumber,
} from './partition_enumerate.js';

// ---------------------------------------------------------------------
// LAYER C — orchestrator: full consensus pipeline + 6 classes + QC fields
// ---------------------------------------------------------------------
export {
  consensus_partition,
  CONSENSUS_CLASS,
  RESOLVING_POWER_CLASS,
} from './partition_consensus.js';
