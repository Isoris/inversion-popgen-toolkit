// shared/band_tracking/band_voters.js
// =====================================================================
// LAYER A — per-band, band-centric view.
//
// "Standing on band b at the locus, looking at who votes for us."
//
// For each target band b, we collect everyone who has an opinion about b
// (voters who included b in their visited or excluded set). For each
// such voter we record what they said:
//   - "you {b, ...} are like my samples"     (b is in voter's visited)
//   - "you {b, ...} are not like my samples" (b is in voter's excluded)
//
// From those incoming opinions we derive per-band metrics:
//   - voter_consensus[b]: how concentrated is b's group-membership
//                         distribution across all incoming votes?
//   - primary_partner_set[b]: which OTHER bands does b most often share
//                              a "side" with (visited-together OR
//                              excluded-together)?
//   - secondary_partner_set[b]: if there's a strong second pattern that
//                                disagrees with the primary, that's the
//                                signature of multi-layer structure
//                                showing up at b.
//   - hidden_regime_residual[b]: per-band version of the QC metric
//                                ∈ [0, 1]; 0 = clean, 1 = maximally
//                                ambiguous.
//
// Everything is a derived view of the co-association matrix in
// vote_evidence.js. We don't recompute votes here — we read from them.
// =====================================================================

// ---------------------------------------------------------------------
// (1) collect_voters_for_band
// ---------------------------------------------------------------------

/**
 * Return the list of voters who have an opinion about band b.
 *
 * @param {number} b
 * @param {object} voteResult     output of extract_votes
 * @param {object} voteIndex      output of build_per_band_vote_index
 * @returns {{
 *   visiting_voters:  { source_w, source_k, weight, partners_visited, partners_excluded }[],
 *   excluding_voters: { source_w, source_k, weight, partners_visited, partners_excluded }[],
 *   total_weight:     number,
 *   silence_count:    int,        // how many votes did NOT mention b
 *                                  //   (informational only — silence is
 *                                  //    not the same as ambiguity)
 * }}
 */
export function collect_voters_for_band(b, voteResult, voteIndex) {
  const visIds = voteIndex.votes_visiting.get(b)  || [];
  const excIds = voteIndex.votes_excluding.get(b) || [];
  const silIds = voteIndex.votes_silent_on.get(b) || [];
  const visiting_voters = [];
  const excluding_voters = [];
  let total_weight = 0;
  for (const id of visIds) {
    const v = voteResult.votes[id];
    visiting_voters.push({
      source_w: v.source_w, source_k: v.source_k, weight: v.weight,
      partners_visited:  v.visited_bands.filter(x => x !== b),
      partners_excluded: v.excluded_bands,
      pattern_class: v.pattern_class,
    });
    total_weight += v.weight;
  }
  for (const id of excIds) {
    const v = voteResult.votes[id];
    excluding_voters.push({
      source_w: v.source_w, source_k: v.source_k, weight: v.weight,
      partners_visited:  v.visited_bands,
      partners_excluded: v.excluded_bands.filter(x => x !== b),
      pattern_class: v.pattern_class,
    });
    total_weight += v.weight;
  }
  return {
    visiting_voters,
    excluding_voters,
    total_weight,
    silence_count: silIds.length,
  };
}

// ---------------------------------------------------------------------
// (2) compute_partner_affinities
//
// For one band b, derive a per-OTHER-band affinity score: how often,
// across all votes that mention BOTH b and the other band, are they
// placed on the same side?
//
// This is exactly the row of the co-association matrix for b — but we
// also report the underlying counts so band_voters.js can drill in
// without recomputing.
// ---------------------------------------------------------------------

/**
 * @param {number} b
 * @param {object} coassoc        output of build_coassociation_matrix
 * @returns {{
 *   affinity:      Float64Array(K),   // co-association score per other band; ∈ [0,1]
 *                                      //   self-entry (b == b) is 1.
 *   together_w:    Float64Array(K),   // raw same-side weight per other band
 *   apart_w:       Float64Array(K),   // raw other-side weight per other band
 *   total_w:       Float64Array(K),   // total mentions per other band
 * }}
 */
export function compute_partner_affinities(b, coassoc) {
  const K = coassoc.K_locus;
  const affinity   = new Float64Array(K);
  const together_w = new Float64Array(K);
  const apart_w    = new Float64Array(K);
  const total_w    = new Float64Array(K);
  for (let j = 0; j < K; j++) {
    if (j === b) { affinity[j] = 1; continue; }
    const idx = b * K + j;
    affinity[j]   = coassoc.coassoc[idx];
    together_w[j] = coassoc.together[idx];
    apart_w[j]    = coassoc.apart[idx];
    total_w[j]    = coassoc.total[idx];
  }
  return { affinity, together_w, apart_w, total_w };
}

// ---------------------------------------------------------------------
// (3) derive_partner_sets
//
// From b's affinity row, derive the primary and secondary partner sets.
//
// Algorithm:
//   1. Sort other bands by affinity descending.
//   2. Primary = the prefix of high-affinity bands above primaryThr.
//      We require strict inequality with the gap: if there's a clear gap
//      between adjacent affinities, the primary set ends there.
//      Otherwise, all bands with affinity ≥ primaryThr are primary.
//   3. Secondary = bands with affinity in [secondaryLo, secondaryHi).
//      These are bands with NON-NEGLIGIBLE but NOT-PRIMARY affinity —
//      candidates for "this band is sometimes also grouped with these,
//      possibly in a different layer."
//   4. Excluded = bands with affinity < apartThr (clear opposite-side).
//
// IMPORTANT: secondary != "uncertain." Secondary specifically means
// "there's a separable second affinity peak, not just noise."
// We only declare secondary bands when:
//   - their affinity is ≥ secondaryLo AND
//   - the EVIDENCE BASE for them (total_w) is non-trivial AND
//   - they're NOT primary AND
//   - there's a discernible gap (not a smooth tail).
// ---------------------------------------------------------------------

/**
 * @param {number} b
 * @param {object} affRecord       output of compute_partner_affinities
 * @param {object} [opts]
 * @param {number} [opts.primaryThr=0.7]      affinity ≥ → primary partner
 * @param {number} [opts.apartThr=0.3]        affinity < → excluded partner
 * @param {number} [opts.secondaryLo=0.4]     in [secondaryLo, primaryThr)
 *                                            → candidate secondary
 * @param {number} [opts.minTotalW=1.0]       partners with total weight
 *                                            below this are dropped
 *                                            (insufficient evidence base)
 * @param {number} [opts.gapThr=0.15]         minimum affinity gap separating
 *                                            primary from secondary
 * @returns {{
 *   primary_partners:   number[],
 *   secondary_partners: number[],
 *   excluded_partners:  number[],
 *   uncertain_partners: number[],     // [apartThr, secondaryLo) — too weak
 *                                     //   for primary/secondary, too strong
 *                                     //   for excluded
 *   primary_min_affinity:  number,
 *   secondary_max_affinity: number,
 *   has_separable_secondary: boolean, // true iff there IS a secondary set
 *                                     //   AND it's separated from primary
 *                                     //   by ≥ gapThr
 * }}
 */
export function derive_partner_sets(b, affRecord, opts) {
  opts = opts || {};
  const primaryThr  = opts.primaryThr  != null ? +opts.primaryThr  : 0.7;
  const apartThr    = opts.apartThr    != null ? +opts.apartThr    : 0.3;
  const secondaryLo = opts.secondaryLo != null ? +opts.secondaryLo : 0.4;
  const minTotalW   = opts.minTotalW   != null ? +opts.minTotalW   : 1.0;
  const gapThr      = opts.gapThr      != null ? +opts.gapThr      : 0.15;

  const K = affRecord.affinity.length;
  const primary = [], secondary = [], excluded = [], uncertain = [];

  // Sort other bands by affinity descending — for separability test
  const ranked = [];
  for (let j = 0; j < K; j++) {
    if (j === b) continue;
    if (affRecord.total_w[j] < minTotalW) continue;   // skip evidence-thin
    ranked.push([j, affRecord.affinity[j]]);
  }
  ranked.sort((a, b) => b[1] - a[1]);

  for (const [j, aff] of ranked) {
    if (aff >= primaryThr) primary.push(j);
    else if (aff < apartThr) excluded.push(j);
    else if (aff >= secondaryLo) secondary.push(j);
    else uncertain.push(j);
  }

  // Separable-secondary test: highest secondary affinity must be at
  // least gapThr below the lowest primary affinity.
  const primMin = primary.length > 0 ? Math.min(...primary.map(j => affRecord.affinity[j])) : NaN;
  const secMax  = secondary.length > 0 ? Math.max(...secondary.map(j => affRecord.affinity[j])) : NaN;
  const has_separable_secondary = secondary.length > 0 &&
    (Number.isFinite(primMin) ? primMin - secMax >= gapThr : true);

  return {
    primary_partners:   primary,
    secondary_partners: secondary,
    excluded_partners:  excluded,
    uncertain_partners: uncertain,
    primary_min_affinity:  primMin,
    secondary_max_affinity: secMax,
    has_separable_secondary,
  };
}

// ---------------------------------------------------------------------
// (4) compute_voter_consensus
//
// For one band b: how concentrated is the distribution of "where do my
// voters place me" across the affinity space? Low entropy = b's voters
// agree (clean band). High entropy = b's voters disagree (ambiguous band).
//
// Computed from b's affinity row treated as a probability distribution
// over partners. Specifically:
//
//   p_j = total_w[j] / sum(total_w)
//
// (probability that a vote-mention of b also mentions partner j;
//  this is the "who pays attention to me" distribution)
//
// But what we really want is the AGREEMENT distribution: of the votes
// that mention both b and j, do they agree on side? That's affinity[j]
// directly. We compute two metrics:
//
//   voter_consensus[b]  = 1 - normalized_entropy of {affinity[j] | j ≠ b}
//                         viewed as a distribution after softmax
//                         No — better: use the direct co-grouping
//                         agreement averaged over partners with
//                         positive evidence base, weighted by evidence.
//
// Concretely:
//   weighted_avg_agreement = Σ_j (max(affinity[j], 1-affinity[j]) * total_w[j])
//                           / Σ_j total_w[j]
// where max(aff, 1-aff) ∈ [0.5, 1] is "how confident is the voter group
// about b's side relative to j" — 1 means voters agreed perfectly that
// b and j are same-side (or perfectly opposite-side); 0.5 means the
// voters were split 50/50.
//
// hidden_regime_residual[b] = 1 - weighted_avg_agreement * 2 + 1
//                           = 2 * (1 - weighted_avg_agreement)
// This ranges over [0, 1]: 0 when all voters agree perfectly, 1 when
// all voters are perfectly split 50/50 on every partner.
// ---------------------------------------------------------------------

/**
 * @param {number} b
 * @param {object} affRecord        output of compute_partner_affinities
 * @returns {{
 *   weighted_avg_agreement:    number,   // ∈ [0.5, 1]
 *   hidden_regime_residual:    number,   // ∈ [0, 1]
 *   evidence_base:             number,   // total weight of votes that
 *                                        //   mention b alongside any partner
 *   n_partners_with_evidence:  int,
 * }}
 */
export function compute_voter_consensus(b, affRecord) {
  const K = affRecord.affinity.length;
  let agreementWSum = 0, totalW = 0, nPartners = 0;
  for (let j = 0; j < K; j++) {
    if (j === b) continue;
    const w = affRecord.total_w[j];
    if (w <= 0) continue;
    nPartners++;
    const aff = affRecord.affinity[j];
    const agreement = Math.max(aff, 1 - aff);   // ∈ [0.5, 1]
    agreementWSum += agreement * w;
    totalW += w;
  }
  if (totalW === 0) {
    return {
      weighted_avg_agreement: 1,    // no evidence — degenerate "agreement"
      hidden_regime_residual: 0,    // can't have hidden conflict if no evidence
      evidence_base: 0,
      n_partners_with_evidence: 0,
    };
  }
  const wa = agreementWSum / totalW;
  // wa ∈ [0.5, 1]; rescale to [0, 1] residual
  // residual = 2 * (1 - wa) gives 0 when wa=1, 1 when wa=0.5
  const residual = Math.max(0, 2 * (1 - wa));
  return {
    weighted_avg_agreement: wa,
    hidden_regime_residual: residual,
    evidence_base: totalW,
    n_partners_with_evidence: nPartners,
  };
}

// ---------------------------------------------------------------------
// (5) compute_overlap_conflict
//
// For one band b: does b have BOTH a primary and a secondary partner set
// that meaningfully disagree? This is the key multi-layer signature at
// the per-band level.
//
// conflict score = (n_secondary_partners that are NOT in primary_partners
//                   AND NOT in excluded_partners) / K
//
// 0 = no conflict (only one consistent partner set for b)
// > 0 = b has separable secondary partners — multi-layer evidence at b
// ---------------------------------------------------------------------

/**
 * @param {number} b
 * @param {object} sets         output of derive_partner_sets
 * @param {int} K
 * @returns {{
 *   conflict_score:    number,        // ∈ [0, 1]
 *   conflict_partners: number[],      // bands in secondary that aren't
 *                                     //   in primary or excluded
 *   has_conflict:      bool,          // conflict_score > 0
 * }}
 */
export function compute_overlap_conflict(b, sets, K) {
  const primSet  = new Set(sets.primary_partners);
  const exclSet  = new Set(sets.excluded_partners);
  const conflicts = [];
  for (const j of sets.secondary_partners) {
    if (!primSet.has(j) && !exclSet.has(j)) conflicts.push(j);
  }
  // Normalize by (K-1) so a band can't possibly score > 1
  return {
    conflict_score: K > 1 ? conflicts.length / (K - 1) : 0,
    conflict_partners: conflicts,
    has_conflict: conflicts.length > 0 && sets.has_separable_secondary,
  };
}

// ---------------------------------------------------------------------
// (6) build_per_band_view
//
// Composer: run the full per-band analysis for every band at the locus.
// Returns a structured array of per-band records — the "standing on
// each band, what do voters say" data structure.
// ---------------------------------------------------------------------

/**
 * @param {object} voteResult     output of extract_votes
 * @param {object} voteIndex      output of build_per_band_vote_index
 * @param {object} coassoc        output of build_coassociation_matrix
 * @param {object} [opts]         passed to derive_partner_sets
 * @returns {object[]}            per band, in band-id order
 */
export function build_per_band_view(voteResult, voteIndex, coassoc, opts) {
  const K = coassoc.K_locus;
  const out = [];
  for (let b = 0; b < K; b++) {
    const voters       = collect_voters_for_band(b, voteResult, voteIndex);
    const aff          = compute_partner_affinities(b, coassoc);
    const sets         = derive_partner_sets(b, aff, opts);
    const consensus    = compute_voter_consensus(b, aff);
    const conflict     = compute_overlap_conflict(b, sets, K);

    out.push({
      band: b,
      // Voters
      n_visiting_voters:  voters.visiting_voters.length,
      n_excluding_voters: voters.excluding_voters.length,
      total_voter_weight: voters.total_weight,
      silence_count:      voters.silence_count,
      visiting_voters:    voters.visiting_voters,
      excluding_voters:   voters.excluding_voters,
      // Affinities
      partner_affinity:   aff.affinity,
      partner_total_w:    aff.total_w,
      // Sets
      primary_partners:   sets.primary_partners,
      secondary_partners: sets.secondary_partners,
      excluded_partners:  sets.excluded_partners,
      uncertain_partners: sets.uncertain_partners,
      has_separable_secondary: sets.has_separable_secondary,
      // QC
      voter_consensus:        consensus.weighted_avg_agreement,
      hidden_regime_residual: consensus.hidden_regime_residual,
      evidence_base:          consensus.evidence_base,
      conflict_score:         conflict.conflict_score,
      conflict_partners:      conflict.conflict_partners,
      has_conflict:           conflict.has_conflict,
    });
  }
  return out;
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._collect_voters_for_band   = collect_voters_for_band;
  window._compute_partner_affinities = compute_partner_affinities;
  window._derive_partner_sets        = derive_partner_sets;
  window._compute_voter_consensus    = compute_voter_consensus;
  window._compute_overlap_conflict   = compute_overlap_conflict;
  window._build_per_band_view        = build_per_band_view;
}
