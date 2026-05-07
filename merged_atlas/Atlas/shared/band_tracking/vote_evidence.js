// shared/band_tracking/vote_evidence.js
// =====================================================================
// Raw vote extraction — shared foundation for both the band-centric
// view (band_voters.js) and the locus-centric view (partition_consensus.js).
//
// A "vote" is a single piece of evidence from one (source_band × source_window)
// projection onto the locus's target bands. Per the agreed semantics:
//
//   The voter (a source band's members at one source window) says:
//     "Your locus bands {visited} are like my samples.
//      Your locus bands {excluded} are not like my samples."
//
// That's directional information — visited and excluded are not symmetric
// (visited carries positive same-state evidence; excluded carries
// other-state evidence). When we ask later "do bands b_i and b_j belong
// together?", we count BOTH "both visited" and "both excluded" as evidence
// of co-grouping, but tagged so we can audit which evidence type fired.
//
// Vote quality: only votes whose source projection has pattern_class
// SUBSET, SUBSET_SPLIT, SPLIT_TWO, or SINGLE are kept. SCATTER/FAN/EMPTY
// votes don't express a partition opinion (per Q1c) and would dilute
// signal — they're recorded but flagged uninformative.
//
// This module produces the data structure both views consume:
//   - per-vote records (the directional evidence)
//   - the co-association matrix (locus-level summary)
//   - per-band incoming-vote indices (band-level summary)
// =====================================================================

import { PATTERN_CLASS } from './projection.js';

// ---------------------------------------------------------------------
// Vote weighting by pattern class quality
// ---------------------------------------------------------------------

// Pattern classes ranked by how confidently they express a partition.
// Used as the default vote weight; can be overridden via opts.
const DEFAULT_VOTE_WEIGHTS = Object.freeze({
  [PATTERN_CLASS.SINGLE]:       1.0,   // strongest — single dominant target
  [PATTERN_CLASS.SUBSET]:       1.0,   // strong — clean visited / excluded split
  [PATTERN_CLASS.SUBSET_SPLIT]: 0.7,   // strong subset, but split internally
  [PATTERN_CLASS.SPLIT_TWO]:    0.5,   // borderline — two-way split, partition-ish
  [PATTERN_CLASS.FAN]:          0.0,   // uninformative
  [PATTERN_CLASS.SCATTER]:      0.0,   // uninformative
  [PATTERN_CLASS.EMPTY]:        0.0,   // empty source — no vote
});

const INFORMATIVE_PATTERNS = new Set([
  PATTERN_CLASS.SINGLE,
  PATTERN_CLASS.SUBSET,
  PATTERN_CLASS.SUBSET_SPLIT,
  PATTERN_CLASS.SPLIT_TWO,
]);

// ---------------------------------------------------------------------
// (1) extract_votes
// ---------------------------------------------------------------------

/**
 * Extract raw votes from the per-source-band projection results.
 *
 * Input: a "voter source" — a list of (source_window, source_band, projection)
 * tuples. Each projection has visited_bands / excluded_bands referring to
 * the LOCUS's bands (the targets). Each is a single vote.
 *
 * In our pipeline, the upstream `kt_combine_*_evidence` produces per-band
 * `projections[k]` where each `projections[k].per_window[i]` is the result
 * of source-band-k projecting onto target-window-i. To get votes ABOUT
 * the locus, we need the inverse: long-range source bands projecting onto
 * the LOCUS. This function operates on whatever shape the caller assembles.
 * Pass an array of per-vote records.
 *
 * @param {object[]} voteRecords  array of:
 *   {
 *     source_w: int,            // source window index (the voter's address)
 *     source_k: int,            // source band id within that window
 *     visited_bands: number[],  // locus bands the voter says are "like me"
 *     excluded_bands: number[], // locus bands the voter says are "not like me"
 *     pattern_class: string,    // from PATTERN_CLASS
 *     n_source: int,            // size of the source band-set
 *     purity_vector?: Float32Array,  // optional; for tie-breaking
 *   }
 * @param {object} [opts]
 * @param {boolean} [opts.dropUninformative=true]   filter out FAN/SCATTER/EMPTY
 * @param {object}  [opts.weights]                  override DEFAULT_VOTE_WEIGHTS
 * @returns {{
 *   votes:         Vote[],
 *   K_locus:       int,           // inferred from max band id seen
 *   n_voters:      int,
 *   n_dropped_uninformative: int,
 * }}
 *
 * Vote = {
 *   id: int,                      // index into votes[]
 *   source_w, source_k: int,
 *   visited_bands: number[],
 *   excluded_bands: number[],
 *   pattern_class: string,
 *   weight: number,               // ∈ [0, 1]
 *   informative: bool,
 * }
 */
export function extract_votes(voteRecords, opts) {
  opts = opts || {};
  const dropUninformative = opts.dropUninformative !== false;
  const weights = Object.assign({}, DEFAULT_VOTE_WEIGHTS, opts.weights || {});

  const votes = [];
  let K_locus = 0;
  let n_dropped = 0;

  for (let i = 0; i < voteRecords.length; i++) {
    const r = voteRecords[i];
    const w = weights[r.pattern_class] != null ? +weights[r.pattern_class] : 0;
    const informative = INFORMATIVE_PATTERNS.has(r.pattern_class);
    if (dropUninformative && !informative) { n_dropped++; continue; }
    // Track maximum band id mentioned for K inference
    for (const b of r.visited_bands  || []) if (b + 1 > K_locus) K_locus = b + 1;
    for (const b of r.excluded_bands || []) if (b + 1 > K_locus) K_locus = b + 1;
    votes.push({
      id: votes.length,
      source_w: r.source_w,
      source_k: r.source_k,
      visited_bands:  Array.from(r.visited_bands  || []).slice().sort((a,b)=>a-b),
      excluded_bands: Array.from(r.excluded_bands || []).slice().sort((a,b)=>a-b),
      pattern_class: r.pattern_class,
      weight: w,
      informative,
      n_source: r.n_source != null ? r.n_source : null,
    });
  }

  // Optional: band_groups[b] = group id for band b. Two bands in the SAME
  // group cannot both be visited by any single vote (they're alternative
  // outcomes within one observation unit, e.g. two clusters of one
  // window/envelope). Used by build_coassociation_matrix to suppress the
  // spurious "both excluded → together" co-grouping for same-group pairs.
  // See FINDING_2026-05-06_both_excluded_bug.md.
  let band_groups = null;
  if (opts.band_groups) {
    if (opts.band_groups.length < K_locus) {
      throw new Error(
        `extract_votes: opts.band_groups length ${opts.band_groups.length} ` +
        `< inferred K_locus ${K_locus}`);
    }
    band_groups = Int32Array.from(opts.band_groups);
  }

  return {
    votes, K_locus,
    n_voters: votes.length,
    n_dropped_uninformative: n_dropped,
    band_groups,
  };
}

// ---------------------------------------------------------------------
// (2) build_coassociation_matrix
//
// For each pair (b_i, b_j), the matrix entry is the SIGNED, WEIGHTED
// co-association count. We track three numbers per pair:
//
//   together[i][j] = sum of weights of votes where b_i and b_j BOTH visited
//                    (or BOTH excluded — same-side evidence — UNLESS the
//                     two bands belong to the same `band_group`, in which
//                     case both-excluded is treated as `apart`; see
//                     FINDING_2026-05-06_both_excluded_bug.md)
//   apart[i][j]    = sum of weights of votes where one is visited and
//                    the other is excluded (different-side evidence)
//   total[i][j]    = sum of weights of votes that mention BOTH bands
//                    (= together + apart; for normalization)
//
// Then: coassoc[i][j] = together / total ∈ [0, 1]
//   = "given a vote that mentions both, how often does it group them"
//
// "Mentions" here means b is in the vote's visited OR excluded set. Bands
// not mentioned by a vote (because the source's projection had ambiguous
// purity for that band) don't contribute to that pair's evidence.
//
// The `band_groups` correction: when two bands are in the same group
// (e.g. two clusters of the same window/envelope — alternative
// arrangements observed at one locus position), they cannot both be
// visited by any single vote (the source has a single arrangement at
// that position). A vote that excludes BOTH of them is therefore not
// evidence that they group together; it's just evidence that the
// source's arrangement differs from both. Without the correction this
// drives spurious co-association inflation when K_arrangements > 2.
// ---------------------------------------------------------------------

/**
 * @param {object} voteResult  output of extract_votes
 * @returns {{
 *   K_locus:       int,
 *   together:      Float64Array(K*K),   // weighted co-grouping count
 *   apart:         Float64Array(K*K),   // weighted opposing-group count
 *   total:         Float64Array(K*K),   // weighted "either side mentioned"
 *   coassoc:       Float64Array(K*K),   // ∈ [0,1] when total > 0; 0 otherwise
 *   n_mentions:    Int32Array(K),       // count of votes mentioning each band
 *   total_weight:  number,              // sum of vote weights
 * }}
 */
export function build_coassociation_matrix(voteResult) {
  const { votes, K_locus, band_groups } = voteResult;
  const K = K_locus;
  const together = new Float64Array(K * K);
  const apart    = new Float64Array(K * K);
  const total    = new Float64Array(K * K);
  const n_mentions = new Int32Array(K);
  let total_weight = 0;
  const useGroups = !!band_groups;

  for (const v of votes) {
    if (v.weight <= 0) continue;
    total_weight += v.weight;
    const vis = new Set(v.visited_bands);
    const exc = new Set(v.excluded_bands);
    for (const b of vis) n_mentions[b]++;
    for (const b of exc) n_mentions[b]++;

    // For every pair (i, j) where i < j, decide the relationship:
    // - both visited      → together
    // - both excluded     → together IF different group, ELSE apart
    //                       (same-group pairs can't be co-grouped: they
    //                        are alternative arrangements at one locus
    //                        position; see FINDING)
    // - one visited, one excluded → apart
    // - either unmentioned → no contribution
    for (let i = 0; i < K; i++) {
      const iVis = vis.has(i), iExc = exc.has(i);
      if (!iVis && !iExc) continue;
      for (let j = i + 1; j < K; j++) {
        const jVis = vis.has(j), jExc = exc.has(j);
        if (!jVis && !jExc) continue;
        const idx = i * K + j;
        const idxR = j * K + i;
        total[idx] += v.weight; total[idxR] += v.weight;

        let isTogether;
        if (iVis && jVis) {
          isTogether = true;
        } else if (iExc && jExc) {
          // Suppress same-group both-excluded co-grouping
          isTogether = useGroups ? (band_groups[i] !== band_groups[j]) : true;
        } else {
          isTogether = false;
        }

        if (isTogether) {
          together[idx] += v.weight; together[idxR] += v.weight;
        } else {
          apart[idx] += v.weight; apart[idxR] += v.weight;
        }
      }
    }
  }

  const coassoc = new Float64Array(K * K);
  for (let i = 0; i < K; i++) {
    coassoc[i * K + i] = 1;   // diagonal — a band always groups with itself
    for (let j = i + 1; j < K; j++) {
      const idx = i * K + j;
      const t = total[idx];
      const c = t > 0 ? together[idx] / t : 0;
      coassoc[idx] = c;
      coassoc[j * K + i] = c;
    }
  }
  return { K_locus: K, together, apart, total, coassoc, n_mentions, total_weight,
           band_groups: band_groups || null };
}

// ---------------------------------------------------------------------
// (3) build_per_band_vote_index
//
// For the band-centric view: index every vote by which target bands it
// mentions. Lets band_voters.js answer "show me everyone who voted ABOUT
// b" in O(votes_for_b) time.
// ---------------------------------------------------------------------

/**
 * @param {object} voteResult  output of extract_votes
 * @returns {{
 *   votes_visiting:   Map<band_id, vote_id[]>,   // votes where this band IS visited
 *   votes_excluding:  Map<band_id, vote_id[]>,   // votes where this band IS excluded
 *   votes_mentioning: Map<band_id, vote_id[]>,   // either of the above
 *   votes_silent_on:  Map<band_id, vote_id[]>,   // band not mentioned by this vote
 * }}
 */
export function build_per_band_vote_index(voteResult) {
  const { votes, K_locus } = voteResult;
  const visiting   = new Map();
  const excluding  = new Map();
  const mentioning = new Map();
  const silent     = new Map();
  for (let b = 0; b < K_locus; b++) {
    visiting.set(b, []); excluding.set(b, []);
    mentioning.set(b, []); silent.set(b, []);
  }
  for (const v of votes) {
    const visSet = new Set(v.visited_bands);
    const excSet = new Set(v.excluded_bands);
    for (let b = 0; b < K_locus; b++) {
      if (visSet.has(b)) {
        visiting.get(b).push(v.id);
        mentioning.get(b).push(v.id);
      } else if (excSet.has(b)) {
        excluding.get(b).push(v.id);
        mentioning.get(b).push(v.id);
      } else {
        silent.get(b).push(v.id);
      }
    }
  }
  return {
    votes_visiting: visiting,
    votes_excluding: excluding,
    votes_mentioning: mentioning,
    votes_silent_on: silent,
  };
}

// ---------------------------------------------------------------------
// Convenience: build voteRecords array from kt_combine_*_evidence's
// projections. The combiner's `projections[k]` is "source band k at the
// LOCUS projecting onto target windows" — that's the OPPOSITE direction
// from what voting needs (we want "long-range source bands voting ABOUT
// the locus's bands"). So this convenience function is provided for
// tests/synthetic data only; real-data assembly will need its own
// reverse-direction projection runner.
//
// For the symmetric case where the SAME projection serves both directions
// (when source and target are interchangeable cohort partitions), this
// helper gives you a clean record list.
// ---------------------------------------------------------------------

/**
 * Reshape per-source-band projections into a flat voteRecords array.
 * Useful for tests; in production you'll typically build voteRecords
 * directly from a long-range → locus projection sweep.
 */
export function voteRecords_from_projections(projectionsBySourceK) {
  const out = [];
  for (let k = 0; k < projectionsBySourceK.length; k++) {
    for (const pw of projectionsBySourceK[k].per_window) {
      out.push({
        source_w: pw.w,
        source_k: k,
        visited_bands:  pw.ve.visited_bands,
        excluded_bands: pw.ve.excluded_bands,
        pattern_class:  pw.cls.pattern_class,
        n_source:       pw.pv.n_source,
        purity_vector:  pw.pv.purity_vector,
      });
    }
  }
  return out;
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._extract_votes               = extract_votes;
  window._build_coassociation_matrix  = build_coassociation_matrix;
  window._build_per_band_vote_index   = build_per_band_vote_index;
  window._voteRecords_from_projections = voteRecords_from_projections;
  window._DEFAULT_VOTE_WEIGHTS         = DEFAULT_VOTE_WEIGHTS;
}
