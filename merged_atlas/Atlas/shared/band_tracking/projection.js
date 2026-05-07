// shared/band_tracking/projection.js  (MINIMAL STUB FOR OPTION B DRIVER)
// =====================================================================
// This stub exists only to satisfy the `import { PATTERN_CLASS } from
// './projection.js'` in vote_evidence.js while running consensus_partition
// on real LG28 data, WITHOUT pulling in the full Layer-2 projection
// machinery (bp_project_bandset_to_target_bands, etc.) from the earlier
// tarballs.
//
// In production this file is the full Layer-2 projection module shipped
// by `Atlas_band_layers_2026-05-05.tar.gz` (tarball #3). Do not commit
// this stub on top of that — let the production file win.
// =====================================================================

export const PATTERN_CLASS = Object.freeze({
  SINGLE:       'SINGLE',
  SUBSET:       'SUBSET',
  SUBSET_SPLIT: 'SUBSET_SPLIT',
  SPLIT_TWO:    'SPLIT_TWO',
  FAN:          'FAN',
  SCATTER:      'SCATTER',
  EMPTY:        'EMPTY',
});
