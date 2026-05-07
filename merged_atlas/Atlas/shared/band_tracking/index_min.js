// shared/band_tracking/index_min.js  (DRIVER FOR OPTION B ONLY)
// =====================================================================
// Minimal re-export surface for running consensus_partition on real
// LG28 data. Skips the upstream modules (single_band, het, hom, iv,
// trajectory, karyotype_model) that aren't needed once voteRecords
// already exist.
//
// Use this from run_LG28.mjs only. The production index.js stays the
// canonical surface for atlas-side code.
// =====================================================================

export { PATTERN_CLASS } from './projection.js';

export {
  extract_votes,
  build_coassociation_matrix,
  build_per_band_vote_index,
  voteRecords_from_projections,
} from './vote_evidence.js';

export {
  collect_voters_for_band,
  compute_partner_affinities,
  derive_partner_sets,
  compute_voter_consensus,
  compute_overlap_conflict,
  build_per_band_view,
} from './band_voters.js';

export {
  enumerate_partitions_as_blocks,
  score_partition_against_coassoc,
  enumerate_and_score_all_partitions,
  select_top_partitions_adaptive,
  partitions_are_compatible,
  bellNumber,
} from './partition_enumerate.js';

export {
  consensus_partition,
  CONSENSUS_CLASS,
  RESOLVING_POWER_CLASS,
} from './partition_consensus.js';
