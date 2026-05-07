// shared/band_tracking/partition_consensus.js
// =====================================================================
// LAYER C — orchestrator + final classification + QC fields.
//
// Ties together:
//   - vote_evidence.js   (raw votes + co-association matrix)
//   - band_voters.js     (per-band view; primary/secondary partner sets;
//                          per-band hidden_regime_residual; conflict score)
//   - partition_enumerate.js  (bruteforce locus-level partitions)
//
// Outputs both views simultaneously (per the "do both" instruction):
//   - locus-centric:   top-N partitions, classification class,
//                      five locus-level QC fields
//   - band-centric:    per-band records with their own QC and
//                      partner-set evidence
//
// Final classification classes (six of them, per spec):
//   CLEAN_PARTITION
//   SOFT_PARTITION
//   AMBIGUOUS_BAND
//   OVERLAPPING_VOTES
//   MULTI_LAYER_STRUCTURE
//   NO_CLEAN_CONSENSUS
//
// Resolving-power class (separate axis):
//   CLEAN_RESOLUTION
//   COMPLEX_BUT_RESOLVABLE
//   LOW_RESOLUTION_REVIEW
//   INVALID_HIDDEN_REGIME_CONFLICT
// =====================================================================

import {
  extract_votes,
  build_coassociation_matrix,
  build_per_band_vote_index,
} from './vote_evidence.js';

import { build_per_band_view } from './band_voters.js';

import {
  enumerate_and_score_all_partitions,
  select_top_partitions_adaptive,
  partitions_are_compatible,
  bellNumber,
} from './partition_enumerate.js';

// ---------------------------------------------------------------------
// Classification classes
// ---------------------------------------------------------------------

export const CONSENSUS_CLASS = Object.freeze({
  CLEAN_PARTITION:       'CLEAN_PARTITION',
  SOFT_PARTITION:        'SOFT_PARTITION',
  AMBIGUOUS_BAND:        'AMBIGUOUS_BAND',
  OVERLAPPING_VOTES:     'OVERLAPPING_VOTES',
  MULTI_LAYER_STRUCTURE: 'MULTI_LAYER_STRUCTURE',
  NO_CLEAN_CONSENSUS:    'NO_CLEAN_CONSENSUS',
});

export const RESOLVING_POWER_CLASS = Object.freeze({
  CLEAN_RESOLUTION:               'CLEAN_RESOLUTION',
  COMPLEX_BUT_RESOLVABLE:         'COMPLEX_BUT_RESOLVABLE',
  LOW_RESOLUTION_REVIEW:          'LOW_RESOLUTION_REVIEW',
  INVALID_HIDDEN_REGIME_CONFLICT: 'INVALID_HIDDEN_REGIME_CONFLICT',
});

// ---------------------------------------------------------------------
// Helper: partition entropy
//
// Given top-N partitions with scores s_1, s_2, ..., s_N, normalize
// to a distribution p_i = s_i / Σ s_j and compute Shannon entropy in
// nats. Normalize to [0, 1] by dividing by ln(N).
// N = 1 → 0 (no spread).
// ---------------------------------------------------------------------

function partition_entropy(picked) {
  if (picked.length <= 1) return 0;
  let sum = 0;
  for (const p of picked) {
    if (Number.isFinite(p.score) && p.score > 0) sum += p.score;
  }
  if (sum === 0) return 0;
  let H = 0;
  for (const p of picked) {
    if (!Number.isFinite(p.score) || p.score <= 0) continue;
    const q = p.score / sum;
    if (q > 0) H -= q * Math.log(q);
  }
  const Hmax = Math.log(picked.length);
  return Hmax > 0 ? H / Hmax : 0;
}

// ---------------------------------------------------------------------
// Decision tree for classification
//
// Given:
//   - top picked partitions (1+; sorted by score desc)
//   - per-band view (each band has hidden_regime_residual, has_conflict, etc.)
// Decide which of the 6 classes fires.
//
// Rules (in order):
//
// 1. NO_CLEAN_CONSENSUS:
//    best partition score < strongScoreThr (default 0.50)
//    OR total_pair_weight is too thin (< minEvidence — default 1.0)
//
// 2. MULTI_LAYER_STRUCTURE:
//    ≥ 2 picked partitions AND any pair among them is INCOMPATIBLE
//    (per partitions_are_compatible) AND each has score ≥
//    multiLayerMinScore (default 0.55).
//
// 3. AMBIGUOUS_BAND:
//    one or more bands with hidden_regime_residual ≥ ambiguousBandThr
//    (default 0.5) AND a clear winner partition exists
//    (best score >= softScoreThr, default 0.65).
//
// 4. OVERLAPPING_VOTES:
//    no MULTI_LAYER detected (no incompatible high-scoring partitions),
//    no ambiguous band, BUT best partition score is in the SOFT range
//    (softScoreThr ≤ best < strongScoreThr) AND co-association matrix
//    has more than fragmentationThr fraction of pairs in the [apartThr,
//    primaryThr] grey zone — i.e., raw votes don't cleanly support the
//    answer.
//
// 5. CLEAN_PARTITION:
//    best score ≥ cleanScoreThr (default 0.85), top picks size = 1.
//
// 6. SOFT_PARTITION:
//    softScoreThr ≤ best < cleanScoreThr, single dominant partition.
// ---------------------------------------------------------------------

function classify_consensus(picked, perBand, coassoc, opts) {
  opts = opts || {};
  const strongScoreThr     = opts.strongScoreThr     != null ? +opts.strongScoreThr     : 0.50;
  const softScoreThr       = opts.softScoreThr       != null ? +opts.softScoreThr       : 0.65;
  const cleanScoreThr      = opts.cleanScoreThr      != null ? +opts.cleanScoreThr      : 0.85;
  const ambiguousBandThr   = opts.ambiguousBandThr   != null ? +opts.ambiguousBandThr   : 0.5;
  const multiLayerMinScore = opts.multiLayerMinScore != null ? +opts.multiLayerMinScore : 0.55;
  const fragmentationThr   = opts.fragmentationThr   != null ? +opts.fragmentationThr   : 0.3;
  const minEvidence        = opts.minEvidence        != null ? +opts.minEvidence        : 1.0;
  const reasoning = [];

  if (picked.length === 0) {
    reasoning.push('NO_CLEAN_CONSENSUS: no scored partitions');
    return { consensus_class: CONSENSUS_CLASS.NO_CLEAN_CONSENSUS, reasoning };
  }

  const best = picked[0];

  // Rule 1 — insufficient evidence or low best score
  if (!Number.isFinite(best.score) || best.score < strongScoreThr ||
      best.total_pair_weight < minEvidence) {
    reasoning.push(`NO_CLEAN_CONSENSUS: best score=${(best.score ?? NaN).toFixed(3)} ` +
                   `< ${strongScoreThr}, or total_pair_weight=` +
                   `${best.total_pair_weight.toFixed(2)} < ${minEvidence}`);
    return { consensus_class: CONSENSUS_CLASS.NO_CLEAN_CONSENSUS, reasoning };
  }

  // Rule 2 — multi-layer vs ambiguous-band.
  // When ≥ 2 incompatible high-scoring partitions exist, we have to
  // decide: is this a TRUE cross-product layer split (different bands
  // disagree) or just ONE ambiguous band dragging the answer (rest of
  // bands clean)?
  //
  // The discriminator: count bands with hidden_regime_residual ≥
  // ambiguousBandThr. If exactly 1 (or a small minority of all bands),
  // the disagreement is concentrated → AMBIGUOUS_BAND. If many bands
  // are ambiguous → MULTI_LAYER_STRUCTURE.
  const ambBands = perBand.filter(b => b.hidden_regime_residual >= ambiguousBandThr);

  if (picked.length >= 2) {
    const eligible = picked.filter(p => p.score >= multiLayerMinScore);
    // Multi-layer detection: a partition counts as a "genuine alternative
    // layer" only if it is incompatible with the BEST partition (i.e., it
    // is neither a refinement of the best nor an aggregation of the best).
    // Pairs of siblings that are mutually incompatible but both refinements
    // of the best are NOT multi-layer evidence — they're noise neighbors of
    // a single truth (e.g. the K=3 partition with one band peeled into a
    // singleton, in nine different ways, produces 9 mutually-incompatible
    // siblings, none of which constitutes a layer).
    let foundIncompatibleAgainstBest = false;
    let firstAlternative = null;
    for (let i = 1; i < eligible.length && !foundIncompatibleAgainstBest; i++) {
      const compat = partitions_are_compatible(best.blocks, eligible[i].blocks);
      if (!compat.compatible) {
        foundIncompatibleAgainstBest = true;
        firstAlternative = eligible[i];
      }
    }
    if (foundIncompatibleAgainstBest) {
      // Concentration check: how many bands have high residual?
      // If only a small minority, this is AMBIGUOUS_BAND, not MULTI_LAYER.
      const K = perBand.length;
      const concentrationThreshold = Math.max(1, Math.ceil(K * 0.34));
      // ↑ default: ≤ 1/3 of bands ambiguous → AMBIGUOUS_BAND;
      //   ≥ 1/3 ambiguous → MULTI_LAYER (genuine cross-product)
      if (ambBands.length > 0 && ambBands.length < concentrationThreshold &&
          best.score >= softScoreThr) {
        reasoning.push(
          `AMBIGUOUS_BAND: ${ambBands.length} of ${K} band(s) ambiguous ` +
          `(below concentration threshold of ${concentrationThreshold}); ` +
          `the multi-layer-looking pattern is dragged by these specific bands: ` +
          ambBands.map(b => `b${b.band}=${b.hidden_regime_residual.toFixed(3)}`).join(', ')
        );
        return {
          consensus_class: CONSENSUS_CLASS.AMBIGUOUS_BAND, reasoning,
          ambiguous_band_ids: ambBands.map(b => b.band),
        };
      }
      // Otherwise: genuine multi-layer
      reasoning.push(
        `MULTI_LAYER_STRUCTURE: best partition (score ${best.score.toFixed(3)}) and ` +
        `alternative (score ${firstAlternative.score.toFixed(3)}) are INCOMPATIBLE ` +
        `(neither refines the other), and ${ambBands.length} band(s) are ambiguous ` +
        `(≥ ${concentrationThreshold} threshold for genuine cross-product)`
      );
      reasoning.push(`  best:        [${best.blocks.map(b => '{' + b.join(',') + '}').join(' | ')}]`);
      reasoning.push(`  alternative: [${firstAlternative.blocks.map(b => '{' + b.join(',') + '}').join(' | ')}]`);
      return { consensus_class: CONSENSUS_CLASS.MULTI_LAYER_STRUCTURE, reasoning };
    }
  }

  // Rule 3 — ambiguous band (without multi-layer geometry):
  // any band with high residual AND clear winner partition exists
  if (ambBands.length > 0 && best.score >= softScoreThr) {
    reasoning.push(
      `AMBIGUOUS_BAND: best partition score=${best.score.toFixed(3)} ≥ ${softScoreThr} ` +
      `but ${ambBands.length} band(s) have hidden_regime_residual ≥ ${ambiguousBandThr}: ` +
      ambBands.map(b => `b${b.band}=${b.hidden_regime_residual.toFixed(3)}`).join(', ')
    );
    return {
      consensus_class: CONSENSUS_CLASS.AMBIGUOUS_BAND, reasoning,
      ambiguous_band_ids: ambBands.map(b => b.band),
    };
  }

  // Rule 4 — overlapping votes: best in the soft band AND coassoc has
  // grey-zone fraction above threshold
  if (best.score >= softScoreThr && best.score < strongScoreThr) {
    // Already excluded by Rule 1; this branch never fires. Logical safety.
    reasoning.push('OVERLAPPING_VOTES: best score in soft band');
    return { consensus_class: CONSENSUS_CLASS.OVERLAPPING_VOTES, reasoning };
  }
  if (best.score < cleanScoreThr) {
    // Inspect grey zone in coassoc — fraction of pairs in [0.3, 0.7]
    const K = coassoc.K_locus;
    let nGrey = 0, nWithEvidence = 0;
    for (let i = 0; i < K; i++) {
      for (let j = i + 1; j < K; j++) {
        if (coassoc.total[i * K + j] > 0) {
          nWithEvidence++;
          const c = coassoc.coassoc[i * K + j];
          if (c >= 0.3 && c <= 0.7) nGrey++;
        }
      }
    }
    const greyFrac = nWithEvidence > 0 ? nGrey / nWithEvidence : 0;
    if (greyFrac >= fragmentationThr) {
      reasoning.push(
        `OVERLAPPING_VOTES: best score=${best.score.toFixed(3)} in soft band, ` +
        `${nGrey}/${nWithEvidence} pairs in grey zone [0.3, 0.7] ` +
        `(fraction=${greyFrac.toFixed(3)} ≥ ${fragmentationThr})`
      );
      return { consensus_class: CONSENSUS_CLASS.OVERLAPPING_VOTES, reasoning };
    }
    reasoning.push(
      `SOFT_PARTITION: best score=${best.score.toFixed(3)} in [${softScoreThr}, ${cleanScoreThr}), ` +
      `grey-zone fraction=${greyFrac.toFixed(3)} < ${fragmentationThr}`
    );
    return { consensus_class: CONSENSUS_CLASS.SOFT_PARTITION, reasoning };
  }

  reasoning.push(`CLEAN_PARTITION: best score=${best.score.toFixed(3)} ≥ ${cleanScoreThr}`);
  return { consensus_class: CONSENSUS_CLASS.CLEAN_PARTITION, reasoning };
}

// ---------------------------------------------------------------------
// Resolving-power class — separate axis based on hidden_regime_residual
// ---------------------------------------------------------------------

function classify_resolving_power(hiddenRegimeResidual) {
  if (!Number.isFinite(hiddenRegimeResidual)) return RESOLVING_POWER_CLASS.LOW_RESOLUTION_REVIEW;
  if (hiddenRegimeResidual < 0.20) return RESOLVING_POWER_CLASS.CLEAN_RESOLUTION;
  if (hiddenRegimeResidual < 0.40) return RESOLVING_POWER_CLASS.COMPLEX_BUT_RESOLVABLE;
  if (hiddenRegimeResidual < 0.60) return RESOLVING_POWER_CLASS.LOW_RESOLUTION_REVIEW;
  return RESOLVING_POWER_CLASS.INVALID_HIDDEN_REGIME_CONFLICT;
}

// ---------------------------------------------------------------------
// (1) consensus_partition — top-level orchestrator
// ---------------------------------------------------------------------

/**
 * Run the full pipeline:
 *   raw vote records → vote evidence
 *                    → per-band view
 *                    → bruteforce partitions + scoring
 *                    → top-N adaptive selection
 *                    → consensus class + QC fields
 *
 * @param {object[]} voteRecords    array of raw vote records
 *                                  (see vote_evidence.js extract_votes)
 * @param {object}   [opts]
 * @param {object}   [opts.partnerSetOpts]   forwarded to derive_partner_sets
 * @param {object}   [opts.classifyOpts]     forwarded to classify_consensus
 * @param {object}   [opts.enumOpts]         forwarded to enumerate_and_score_all_partitions
 * @param {object}   [opts.selectOpts]       forwarded to select_top_partitions_adaptive
 * @returns {{
 *   K_locus:                    int,
 *   votes:                      ReturnType<extract_votes>,
 *   coassoc:                    ReturnType<build_coassociation_matrix>,
 *   per_band:                   ReturnType<build_per_band_view>,
 *   top_partitions:             scoredPartition[],
 *   consensus_class:            CONSENSUS_CLASS,
 *   ambiguous_band_ids?:        number[],
 *   reasoning:                  string[],
 *   pca_vote_consensus_score:   number,
 *   pca_hidden_regime_residual: number,
 *   pca_partition_entropy:      number,
 *   pca_overlap_conflict_score: number,
 *   pca_resolving_power_class:  RESOLVING_POWER_CLASS,
 *   bruteforce_used:            bool,        // false if K exceeded cap
 *   n_partitions_total:         int,
 * }}
 */
export function consensus_partition(voteRecords, opts) {
  opts = opts || {};
  const reasoning = [];

  // Stage 1 — vote evidence
  const voteResult = extract_votes(voteRecords, opts);
  reasoning.push(`Stage 1: extracted ${voteResult.n_voters} informative votes ` +
                 `(${voteResult.n_dropped_uninformative} dropped as uninformative)`);
  if (voteResult.band_groups) {
    const nGroups = new Set(Array.from(voteResult.band_groups)).size;
    reasoning.push(`Stage 1: band_groups correction active (${nGroups} groups across ${voteResult.K_locus} bands)`);
  }
  if (voteResult.K_locus === 0) {
    return {
      K_locus: 0, votes: voteResult, coassoc: null, per_band: [],
      top_partitions: [], consensus_class: CONSENSUS_CLASS.NO_CLEAN_CONSENSUS,
      reasoning: [...reasoning, 'No bands mentioned in any vote — empty locus'],
      pca_vote_consensus_score: 0, pca_hidden_regime_residual: 1,
      pca_partition_entropy: 0, pca_overlap_conflict_score: 0,
      pca_resolving_power_class: RESOLVING_POWER_CLASS.INVALID_HIDDEN_REGIME_CONFLICT,
      bruteforce_used: false, n_partitions_total: 0,
    };
  }

  const coassoc = build_coassociation_matrix(voteResult);
  const voteIndex = build_per_band_vote_index(voteResult);

  // Stage 2 — per-band view (band-centric)
  const perBand = build_per_band_view(voteResult, voteIndex, coassoc, opts.partnerSetOpts);
  reasoning.push(`Stage 2: per-band view built for K_locus=${voteResult.K_locus}`);

  // Stage 3 — locus-centric bruteforce
  let scoredResult = null;
  let bruteforce_used = true;
  try {
    scoredResult = enumerate_and_score_all_partitions(coassoc, opts.enumOpts);
  } catch (e) {
    reasoning.push(`Stage 3: bruteforce FAILED — ${e.message}`);
    bruteforce_used = false;
  }
  if (!scoredResult) {
    reasoning.push(`Stage 3: K=${voteResult.K_locus} exceeded K_cap; ` +
                   `bruteforce skipped. Caller should fall back to greedy ` +
                   `(kt_infer_macro_band_groups).`);
    bruteforce_used = false;
  } else {
    reasoning.push(`Stage 3: enumerated ${scoredResult.n_partitions_total} partitions ` +
                   `(Bell number for K=${voteResult.K_locus}); ` +
                   `kept top ${scoredResult.scored.length} by score`);
  }

  const picked = scoredResult
    ? select_top_partitions_adaptive(scoredResult.scored, opts.selectOpts)
    : [];
  if (picked.length > 0) {
    reasoning.push(`Stage 3: top-${picked.length} partitions within δ of best ` +
                   `(scores: ${picked.map(p => p.score.toFixed(3)).join(', ')})`);
  }

  // Stage 4 — classification
  const cls = classify_consensus(picked, perBand, coassoc, opts.classifyOpts);
  for (const r of cls.reasoning) reasoning.push(`Stage 4: ${r}`);

  // Stage 5 — QC fields
  const best = picked[0];
  const pca_vote_consensus_score = best && Number.isFinite(best.score) ? best.score : 0;
  // Per-locus hidden_regime_residual is NOT just 1 - vote consensus; it
  // also folds in the per-band ambiguity. We define it as the MAX of:
  //   - 1 - pca_vote_consensus_score  (locus-level)
  //   - mean of per_band[*].hidden_regime_residual  (band-level)
  // The locus is "as broken as its weakest reading."
  let perBandResidualSum = 0, perBandResidualN = 0;
  for (const b of perBand) {
    if (b.evidence_base > 0) {
      perBandResidualSum += b.hidden_regime_residual;
      perBandResidualN++;
    }
  }
  const meanBandResidual = perBandResidualN > 0 ? perBandResidualSum / perBandResidualN : 0;
  const locusResidual    = 1 - pca_vote_consensus_score;
  const pca_hidden_regime_residual = Math.max(locusResidual, meanBandResidual);

  const pca_partition_entropy = partition_entropy(picked);

  // pca_overlap_conflict_score: per-band conflict score, averaged over
  // bands with non-trivial evidence
  let conflictSum = 0, conflictN = 0;
  for (const b of perBand) {
    if (b.evidence_base > 0) {
      conflictSum += b.conflict_score;
      conflictN++;
    }
  }
  const pca_overlap_conflict_score = conflictN > 0 ? conflictSum / conflictN : 0;

  const pca_resolving_power_class = classify_resolving_power(pca_hidden_regime_residual);

  reasoning.push(`Stage 5: pca_vote_consensus_score=${pca_vote_consensus_score.toFixed(3)}`);
  reasoning.push(`Stage 5: pca_hidden_regime_residual=${pca_hidden_regime_residual.toFixed(3)} ` +
                 `(locus=${locusResidual.toFixed(3)}, ` +
                 `band_mean=${meanBandResidual.toFixed(3)})`);
  reasoning.push(`Stage 5: pca_partition_entropy=${pca_partition_entropy.toFixed(3)}`);
  reasoning.push(`Stage 5: pca_overlap_conflict_score=${pca_overlap_conflict_score.toFixed(3)}`);
  reasoning.push(`Stage 5: pca_resolving_power_class=${pca_resolving_power_class}`);

  return {
    K_locus: voteResult.K_locus,
    votes: voteResult,
    coassoc,
    per_band: perBand,
    top_partitions: picked,
    consensus_class: cls.consensus_class,
    ambiguous_band_ids: cls.ambiguous_band_ids,
    reasoning,
    pca_vote_consensus_score,
    pca_hidden_regime_residual,
    pca_partition_entropy,
    pca_overlap_conflict_score,
    pca_resolving_power_class,
    bruteforce_used,
    n_partitions_total: scoredResult ? scoredResult.n_partitions_total : 0,
  };
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._consensus_partition         = consensus_partition;
  window._classify_consensus           = classify_consensus;
  window._classify_resolving_power     = classify_resolving_power;
  window._partition_entropy            = partition_entropy;
  window._CONSENSUS_CLASS              = CONSENSUS_CLASS;
  window._RESOLVING_POWER_CLASS        = RESOLVING_POWER_CLASS;
}
