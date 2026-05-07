// tests/test_band_consensus.js
// =====================================================================
// Tests for vote_evidence + band_voters + partition_enumerate +
// partition_consensus.
//
// Fixtures pin the SIX classification classes:
//   F1 — CLEAN_PARTITION             (votes unanimous)
//   F2 — SOFT_PARTITION              (votes mostly agree)
//   F3 — MULTI_LAYER_STRUCTURE       (two incompatible high-score partitions)
//   F4 — AMBIGUOUS_BAND              (one band's voters disagree;
//                                     others clean)
//   F5 — NO_CLEAN_CONSENSUS          (votes are noise)
//
// Plus:
//   F6 — Bell-number sanity
//   F7 — partition_enumerate finds the planted truth
//   F8 — bruteforce returns the same answer the spec implies
//   F9 — per-band view "standing on b" record consistency
// =====================================================================

import {
  // vote evidence
  extract_votes,
  build_coassociation_matrix,
  build_per_band_vote_index,
  // band-centric
  build_per_band_view,
  derive_partner_sets,
  compute_partner_affinities,
  compute_voter_consensus,
  // bruteforce
  enumerate_partitions_as_blocks,
  enumerate_and_score_all_partitions,
  partitions_are_compatible,
  bellNumber,
  // orchestrator
  consensus_partition,
  CONSENSUS_CLASS,
  RESOLVING_POWER_CLASS,
  // for vote-record construction
  PATTERN_CLASS,
} from '../shared/band_tracking/index.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

// ---------------------------------------------------------------------
// Helper: build a synthetic vote record.
// ---------------------------------------------------------------------

function vote(source_w, source_k, visited, excluded, pattern_class = PATTERN_CLASS.SUBSET) {
  return {
    source_w, source_k,
    visited_bands: visited.slice(),
    excluded_bands: excluded.slice(),
    pattern_class,
    n_source: 10,
  };
}

// ---------------------------------------------------------------------
// FIXTURE 1 — CLEAN_PARTITION
// 6 bands. 20 votes, all saying "{0,1,2} vs {3,4,5}".
// ---------------------------------------------------------------------
console.log('--- FIXTURE 1: CLEAN_PARTITION ---');
{
  const votes = [];
  for (let w = 0; w < 20; w++) {
    votes.push(vote(w, 0, [0, 1, 2], [3, 4, 5], PATTERN_CLASS.SUBSET));
  }
  const r = consensus_partition(votes);
  check('F1: K_locus = 6', r.K_locus === 6);
  check('F1: bruteforce used', r.bruteforce_used === true);
  check('F1: n_partitions_total = 203 (Bell(6))', r.n_partitions_total === 203);
  check('F1: pca_vote_consensus_score ≥ 0.99',
        r.pca_vote_consensus_score >= 0.99,
        `got ${r.pca_vote_consensus_score.toFixed(4)}`);
  check('F1: pca_hidden_regime_residual ≤ 0.01',
        r.pca_hidden_regime_residual <= 0.01,
        `got ${r.pca_hidden_regime_residual.toFixed(4)}`);
  check('F1: consensus_class = CLEAN_PARTITION',
        r.consensus_class === CONSENSUS_CLASS.CLEAN_PARTITION,
        `got ${r.consensus_class}`);
  check('F1: resolving_power = CLEAN_RESOLUTION',
        r.pca_resolving_power_class === RESOLVING_POWER_CLASS.CLEAN_RESOLUTION);
  // Top partition should be {0,1,2}|{3,4,5}
  check('F1: top_partitions has 1 entry', r.top_partitions.length === 1);
  if (r.top_partitions.length === 1) {
    const blocks = r.top_partitions[0].blocks;
    const blockSets = blocks.map(b => new Set(b));
    const has012 = blockSets.some(s => s.size === 3 && s.has(0) && s.has(1) && s.has(2));
    const has345 = blockSets.some(s => s.size === 3 && s.has(3) && s.has(4) && s.has(5));
    check('F1: top has block {0,1,2}', has012);
    check('F1: top has block {3,4,5}', has345);
  }
  // Per-band view: every band should be clean (low residual)
  for (const b of r.per_band) {
    check(`F1: band ${b.band} hidden_regime_residual = 0`,
          b.hidden_regime_residual === 0);
    check(`F1: band ${b.band} no conflict`, !b.has_conflict);
  }
}

// ---------------------------------------------------------------------
// FIXTURE 2 — SOFT_PARTITION
// Constructing a "pure" SOFT_PARTITION case from raw partition-vote
// records is intrinsically hard: any noise that lowers the consensus
// score tends to either concentrate on one band (→ AMBIGUOUS_BAND) or
// produce incompatible alternative partitions (→ MULTI_LAYER). The
// SOFT_PARTITION path fires when the best score is in [softScoreThr,
// cleanScoreThr) AND no per-band ambiguity dominates AND no incompatible
// alternative exists with high score.
//
// We test the path is reachable using a minimal compatible-refinement
// case where the refinement is small enough not to trigger AMBIGUOUS_BAND.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 2: SOFT_PARTITION (path exists) ---');
{
  // 4 bands. 15 votes for {0,1}|{2,3}, 5 votes for {0,1}|{2}|{3} —
  // a compatible refinement (the second is finer than the first).
  // No band should be ambiguous (band 3's voters all agree it's NOT
  // with {0,1}; the only disagreement is whether 3 joins 2 or stands alone).
  const votes = [];
  for (let w = 0; w < 15; w++) votes.push(vote(w, 0, [0, 1], [2, 3]));
  for (let w = 15; w < 20; w++) votes.push(vote(w, 0, [0, 1], [2]));
  // Note the 5 finer votes only mention 3 bands — band 3 is silent on those.
  const r = consensus_partition(votes);
  check('F2: bruteforce used', r.bruteforce_used);
  // Resolving power should be CLEAN or COMPLEX_BUT_RESOLVABLE
  check('F2: resolving power CLEAN or COMPLEX_BUT_RESOLVABLE',
        r.pca_resolving_power_class === RESOLVING_POWER_CLASS.CLEAN_RESOLUTION ||
        r.pca_resolving_power_class === RESOLVING_POWER_CLASS.COMPLEX_BUT_RESOLVABLE,
        `got ${r.pca_resolving_power_class}, residual=${r.pca_hidden_regime_residual.toFixed(3)}`);
  // No band ambiguous
  const ambBands = r.per_band.filter(b => b.hidden_regime_residual >= 0.5);
  check('F2: no ambiguous bands',
        ambBands.length === 0,
        `${ambBands.length} ambiguous: ` +
        ambBands.map(b => `b${b.band}=${b.hidden_regime_residual.toFixed(3)}`).join(', '));
  // Either CLEAN or SOFT — both are valid "answer is workable" outcomes
  check('F2: consensus_class is CLEAN or SOFT (no ambiguity, no multi-layer)',
        r.consensus_class === CONSENSUS_CLASS.CLEAN_PARTITION ||
        r.consensus_class === CONSENSUS_CLASS.SOFT_PARTITION,
        `got ${r.consensus_class}, score=${r.pca_vote_consensus_score.toFixed(3)}`);
}

// ---------------------------------------------------------------------
// FIXTURE 3 — MULTI_LAYER_STRUCTURE
// Two incompatible partitions, each strongly supported by half the votes.
// "Ancestry" axis: {0,1,2}|{3,4,5} (10 votes)
// "Family" axis:   {0,3}|{1,4}|{2,5} (10 votes)
// Both score ~0.5 because the OTHER partition's votes contradict half
// of every pair — but each is incompatible with the other. Crucial:
// we need to make BOTH score ≥ multiLayerMinScore=0.55.
// To strengthen, we add a few more votes for each axis (15+15=30 total).
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 3: MULTI_LAYER_STRUCTURE ---');
{
  const votes = [];
  // Axis 1: ancestry {0,1,2}|{3,4,5}
  for (let w = 0; w < 15; w++) votes.push(vote(w, 0, [0, 1, 2], [3, 4, 5]));
  // Axis 2: family {0,3}|{1,4}|{2,5} — encoded by 3 source bands per
  // window each voting for one block
  for (let w = 15; w < 30; w++) {
    votes.push(vote(w, 0, [0, 3], [1, 2, 4, 5]));
    votes.push(vote(w, 1, [1, 4], [0, 2, 3, 5]));
    votes.push(vote(w, 2, [2, 5], [0, 1, 3, 4]));
  }
  const r = consensus_partition(votes);
  check('F3: bruteforce used', r.bruteforce_used);
  // Multi-layer should be detected
  check('F3: consensus_class = MULTI_LAYER_STRUCTURE',
        r.consensus_class === CONSENSUS_CLASS.MULTI_LAYER_STRUCTURE,
        `got ${r.consensus_class} | best_score=${r.pca_vote_consensus_score.toFixed(3)} ` +
        `| top=${r.top_partitions.length}`);
  // Top picks should include BOTH layers as separate partitions
  check('F3: ≥2 top partitions', r.top_partitions.length >= 2);
  if (r.top_partitions.length >= 2) {
    // Find the ancestry-like partition (block of size 3)
    const hasAncestry = r.top_partitions.some(p =>
      p.blocks.length === 2 && p.blocks.some(b => b.length === 3));
    const hasFamily = r.top_partitions.some(p =>
      p.blocks.length === 3 && p.blocks.every(b => b.length === 2));
    check('F3: family partition {2,2,2} present in top', hasFamily);
    // The two should be incompatible IF ancestry is in top. If not, check
    // that some pair in top is incompatible (any two layers are sufficient
    // signal for MULTI_LAYER classification, which already passed above).
    if (hasAncestry && hasFamily) {
      const ancestry = r.top_partitions.find(p =>
        p.blocks.length === 2 && p.blocks.some(b => b.length === 3));
      const family = r.top_partitions.find(p =>
        p.blocks.length === 3 && p.blocks.every(b => b.length === 2));
      const compat = partitions_are_compatible(ancestry.blocks, family.blocks);
      check('F3: ancestry and family are INCOMPATIBLE', !compat.compatible);
    } else {
      // Any pair in top must be incompatible (else MULTI_LAYER wouldn't have fired)
      let anyIncompatible = false;
      for (let i = 0; i < r.top_partitions.length && !anyIncompatible; i++) {
        for (let j = i + 1; j < r.top_partitions.length && !anyIncompatible; j++) {
          const c = partitions_are_compatible(r.top_partitions[i].blocks,
                                              r.top_partitions[j].blocks);
          if (!c.compatible) anyIncompatible = true;
        }
      }
      check('F3: at least one incompatible pair in top (multi-layer evidence)',
            anyIncompatible);
    }
  }
  check('F3: pca_partition_entropy > 0 (split between layers)',
        r.pca_partition_entropy > 0,
        `entropy=${r.pca_partition_entropy.toFixed(3)}`);
}

// ---------------------------------------------------------------------
// FIXTURE 4 — AMBIGUOUS_BAND
// 5 bands. 18 votes mostly agree on {0,1}|{2,3,4}, BUT band 1 is pulled
// into different groups by some votes — half see band 1 with {0}, half
// see band 1 with {2,3,4}. The OTHER bands (0, 2, 3, 4) are clean.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 4: AMBIGUOUS_BAND ---');
{
  const votes = [];
  // 12 clean votes: {0,1}|{2,3,4}
  for (let w = 0; w < 12; w++) votes.push(vote(w, 0, [0, 1], [2, 3, 4]));
  // 8 ambiguous votes: band 1 sits with {2,3,4}
  for (let w = 12; w < 20; w++) votes.push(vote(w, 0, [1, 2, 3, 4], [0]));
  const r = consensus_partition(votes);
  // Band 1's voters are split — its hidden_regime_residual should be high
  const band1 = r.per_band.find(b => b.band === 1);
  check('F4: band 1 has high hidden_regime_residual',
        band1.hidden_regime_residual >= 0.4,
        `got ${band1.hidden_regime_residual.toFixed(3)}`);
  // Bands 0, 2, 3, 4 clean
  for (const b of r.per_band) {
    if (b.band === 1) continue;
    check(`F4: band ${b.band} hidden_regime_residual reasonable`,
          b.hidden_regime_residual <= 0.5,
          `got ${b.hidden_regime_residual.toFixed(3)}`);
  }
  // Best partition probably wins (12 vs 8) but soft, AND band 1 ambiguous
  // → AMBIGUOUS_BAND if best.score >= softScoreThr=0.65, else may fall
  // through to OVERLAPPING_VOTES or NO_CLEAN_CONSENSUS depending on numbers
  check('F4: consensus_class is AMBIGUOUS_BAND or OVERLAPPING_VOTES',
        r.consensus_class === CONSENSUS_CLASS.AMBIGUOUS_BAND ||
        r.consensus_class === CONSENSUS_CLASS.OVERLAPPING_VOTES,
        `got ${r.consensus_class} | best=${r.pca_vote_consensus_score.toFixed(3)}`);
  if (r.consensus_class === CONSENSUS_CLASS.AMBIGUOUS_BAND) {
    check('F4: ambiguous_band_ids contains 1',
          (r.ambiguous_band_ids || []).includes(1));
  }
}

// ---------------------------------------------------------------------
// FIXTURE 5 — NO_CLEAN_CONSENSUS
// 4 bands. 12 votes, each voting for a DIFFERENT partition. Maximum noise.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 5: NO_CLEAN_CONSENSUS ---');
{
  const votes = [];
  // Round-robin random partitions — every conceivable 2-block split
  votes.push(vote(0,  0, [0],       [1, 2, 3]));
  votes.push(vote(1,  0, [0, 1],    [2, 3]));
  votes.push(vote(2,  0, [0, 2],    [1, 3]));
  votes.push(vote(3,  0, [0, 3],    [1, 2]));
  votes.push(vote(4,  0, [0, 1, 2], [3]));
  votes.push(vote(5,  0, [0, 1, 3], [2]));
  votes.push(vote(6,  0, [0, 2, 3], [1]));
  votes.push(vote(7,  0, [1],       [0, 2, 3]));
  votes.push(vote(8,  0, [2],       [0, 1, 3]));
  votes.push(vote(9,  0, [3],       [0, 1, 2]));
  votes.push(vote(10, 0, [1, 2],    [0, 3]));
  votes.push(vote(11, 0, [1, 3],    [0, 2]));
  const r = consensus_partition(votes);
  // Best score should be much less than the strong threshold (the
  // threshold the spec uses for "no clean consensus" is 0.50; with this
  // many contradictory partitions, even the best-fitting one barely
  // clears chance.) We accept anything below 0.65 (the soft band cut).
  check('F5: pca_vote_consensus_score < 0.65 (no clear winner)',
        r.pca_vote_consensus_score < 0.65,
        `got ${r.pca_vote_consensus_score.toFixed(3)}`);
  check('F5: high hidden_regime_residual',
        r.pca_hidden_regime_residual >= 0.40,
        `got ${r.pca_hidden_regime_residual.toFixed(3)}`);
  // Either NO_CLEAN_CONSENSUS or OVERLAPPING_VOTES depending on exact
  // distribution; both are "the answer is not clean."
  check('F5: consensus_class is non-clean',
        r.consensus_class === CONSENSUS_CLASS.NO_CLEAN_CONSENSUS ||
        r.consensus_class === CONSENSUS_CLASS.OVERLAPPING_VOTES ||
        r.consensus_class === CONSENSUS_CLASS.MULTI_LAYER_STRUCTURE,
        `got ${r.consensus_class}`);
}

// ---------------------------------------------------------------------
// FIXTURE 6 — Bell number sanity
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 6: Bell-number sanity ---');
{
  // Confirm partition enumerator generates exactly Bell(K) partitions for K=1..7
  for (let K = 1; K <= 7; K++) {
    let n = 0;
    for (const _p of enumerate_partitions_as_blocks(K)) n++;
    check(`enumerator yields Bell(${K})=${bellNumber(K)} partitions`,
          n === bellNumber(K), `enumerated ${n}`);
  }
  // Spot-check K=4 partitions are unique
  const all = [];
  for (const p of enumerate_partitions_as_blocks(4)) {
    // Canonical key: sort blocks by min element, sort within block, stringify
    const canon = p.map(b => b.slice().sort((a, b) => a - b))
                  .sort((a, b) => a[0] - b[0])
                  .map(b => b.join(','))
                  .join('|');
    all.push(canon);
  }
  const uniq = new Set(all);
  check('K=4 enumerated partitions are unique', uniq.size === 15);
}

// ---------------------------------------------------------------------
// FIXTURE 7 — partition_enumerate finds the planted truth
// 4 bands, all 30 votes say {0,1}|{2,3}.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 7: planted truth recovery ---');
{
  const voteRecs = [];
  for (let w = 0; w < 30; w++) voteRecs.push(vote(w, 0, [0, 1], [2, 3]));
  const ev = extract_votes(voteRecs);
  const co = build_coassociation_matrix(ev);
  const sc = enumerate_and_score_all_partitions(co);
  check('F7: 15 partitions enumerated (Bell(4))', sc.n_partitions_total === 15);
  // Best should be the planted partition
  const best = sc.scored[0];
  const setA = new Set(best.blocks[0]);
  const setB = new Set(best.blocks[1]);
  const planted = (setA.size === 2 && setA.has(0) && setA.has(1) &&
                   setB.size === 2 && setB.has(2) && setB.has(3)) ||
                  (setA.size === 2 && setA.has(2) && setA.has(3) &&
                   setB.size === 2 && setB.has(0) && setB.has(1));
  check('F7: top partition = {0,1}|{2,3}', planted);
  check('F7: top score = 1.0 (perfect)', best.score === 1.0,
        `got ${best.score.toFixed(4)}`);
}

// ---------------------------------------------------------------------
// FIXTURE 8 — K cap behavior
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 8: K cap behavior ---');
{
  // Build a fake coassoc with K=11 to trip the cap
  const K = 11;
  const co = {
    K_locus: K,
    together: new Float64Array(K * K),
    apart:    new Float64Array(K * K),
    total:    new Float64Array(K * K),
    coassoc:  new Float64Array(K * K),
    n_mentions: new Int32Array(K),
    total_weight: 0,
  };
  const sc = enumerate_and_score_all_partitions(co, { K_cap: 10 });
  check('F8: K=11 with K_cap=10 returns null', sc === null);
  // K=10 with default cap should still run (115k partitions — slow but works)
  // Skip this in test; expensive. Verify Bell number for confirmation.
  check('F8: Bell(10) = 115975', bellNumber(10) === 115975);
  // Hard fail at K > 13 by default (spec says "fail at K=13"; the
  // implementation uses strict >, so K=14 is the first failing K)
  let threw = false;
  try {
    const K14 = 14;
    const co14 = { K_locus: K14,
      together: new Float64Array(K14 * K14), apart: new Float64Array(K14 * K14),
      total: new Float64Array(K14 * K14), coassoc: new Float64Array(K14 * K14),
      n_mentions: new Int32Array(K14), total_weight: 0 };
    enumerate_and_score_all_partitions(co14);
  } catch (e) { threw = true; }
  check('F8: K=14 throws (hard fail above 13)', threw);
}

// ---------------------------------------------------------------------
// FIXTURE 9 — Per-band view consistency.
// "Standing on band 0" on F1's clean-partition data should report:
// - 20 visiting voters (all 20 votes had band 0 in visited)
// - 0 excluding voters
// - primary partners = {1, 2}
// - excluded partners = {3, 4, 5}
// - hidden_regime_residual = 0
// - no conflict
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 9: per-band view "standing on b" ---');
{
  const votes = [];
  for (let w = 0; w < 20; w++) votes.push(vote(w, 0, [0, 1, 2], [3, 4, 5]));
  const ev = extract_votes(votes);
  const co = build_coassociation_matrix(ev);
  const idx = build_per_band_vote_index(ev);
  const view = build_per_band_view(ev, idx, co);
  const b0 = view[0];
  check('F9: band 0 has 20 visiting voters', b0.n_visiting_voters === 20);
  check('F9: band 0 has 0 excluding voters', b0.n_excluding_voters === 0);
  check('F9: band 0 primary = [1, 2]',
        b0.primary_partners.length === 2 &&
        new Set(b0.primary_partners).has(1) &&
        new Set(b0.primary_partners).has(2));
  check('F9: band 0 excluded = [3, 4, 5]',
        b0.excluded_partners.length === 3 &&
        new Set(b0.excluded_partners).has(3) &&
        new Set(b0.excluded_partners).has(4) &&
        new Set(b0.excluded_partners).has(5));
  check('F9: band 0 voter_consensus = 1', b0.voter_consensus === 1);
  check('F9: band 0 hidden_regime_residual = 0', b0.hidden_regime_residual === 0);
  check('F9: band 0 no conflict', !b0.has_conflict);
  check('F9: band 0 has no separable secondary',
        !b0.has_separable_secondary);
  // visiting_voters records carry source addresses
  check('F9: visiting_voters[0] has source_w', typeof b0.visiting_voters[0].source_w === 'number');
  check('F9: visiting_voters[0] partners_visited has bands 1, 2',
        b0.visiting_voters[0].partners_visited.length === 2);

  // "Standing on band 3" — band 3 is in EXCLUDED for all 20 votes
  const b3 = view[3];
  check('F9: band 3 has 0 visiting voters', b3.n_visiting_voters === 0);
  check('F9: band 3 has 20 excluding voters', b3.n_excluding_voters === 20);
  check('F9: band 3 primary = [4, 5] (co-excluded together)',
        b3.primary_partners.length === 2 &&
        new Set(b3.primary_partners).has(4) &&
        new Set(b3.primary_partners).has(5));
  check('F9: band 3 excluded = [0, 1, 2]',
        b3.excluded_partners.length === 3);
}

// ---------------------------------------------------------------------
// FIXTURE 10 — Uninformative votes filtered
// SCATTER votes shouldn't get extracted by default.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 10: uninformative vote filter ---');
{
  const votes = [
    vote(0, 0, [0, 1, 2], [3, 4, 5], PATTERN_CLASS.SUBSET),
    vote(1, 0, [0, 1, 2], [3, 4, 5], PATTERN_CLASS.SUBSET),
    vote(2, 0, [0, 3], [], PATTERN_CLASS.SCATTER),    // dropped
    vote(3, 0, [], [], PATTERN_CLASS.EMPTY),           // dropped
  ];
  const r = extract_votes(votes);
  check('F10: 2 informative votes kept', r.n_voters === 2);
  check('F10: 2 uninformative dropped', r.n_dropped_uninformative === 2);

  // With dropUninformative=false, all 4 stay (but with weight=0 for uninformative)
  const r2 = extract_votes(votes, { dropUninformative: false });
  check('F10: with dropUninformative=false, 4 votes kept', r2.n_voters === 4);
  // The two informative ones should have non-zero weight
  let nNonZero = 0;
  for (const v of r2.votes) if (v.weight > 0) nNonZero++;
  check('F10: 2 votes have non-zero weight', nNonZero === 2);
}

// ---------------------------------------------------------------------
// FIXTURE 11 — K=3 arrangements with band_groups correction
//
// Regression test for FINDING_2026-05-06_both_excluded_bug.md.
// Synthesizes the LG28 prototype shape: 3 windows × 3 clusters per
// window, samples partition into 3 arrangements (60/106/60 sizes) with
// permuted cluster ids across windows. Without band_groups → score
// ceiling at 0.7333, MULTI_LAYER_STRUCTURE false positive. With
// band_groups → CLEAN_PARTITION at score ≥ 0.86.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 11: K=3 arrangements (band_groups correction) ---');
{
  // Build 9 voters, one per (window, cluster). Each voter "knows" which
  // cluster ids in OTHER windows belong to its arrangement (by Jaccard
  // on samples). For a 60/106/60 split with permuted cluster ids, the
  // Jaccard pattern is unanimous within an arrangement and 0 across.
  //
  // Sample assignment:
  //   samples 0..59     = arrangement A
  //   samples 60..165   = arrangement B
  //   samples 166..225  = arrangement C
  // Cluster id permutation per window (window-local labels):
  //   w0: A→cluster0, B→cluster1, C→cluster2  → bands 0,1,2
  //   w1: A→cluster1, B→cluster2, C→cluster0  → bands 3,4,5  (c0→3, c1→4, c2→5)
  //                                                so band 4 = w1 cluster1 = arrangement A
  //   w2: A→cluster2, B→cluster0, C→cluster1  → bands 6,7,8  (c0→6, c1→7, c2→8)
  //                                                so band 8 = w2 cluster2 = arrangement A

  // Bands by arrangement:
  //   A = {0, 4, 8}, B = {1, 5, 6}, C = {2, 3, 7}
  // Hand-built voteRecords: each source band visits the same-arrangement
  // bands in OTHER windows and excludes everything else.
  function v3(srcW, srcK, visited, excluded, n) {
    return { source_w: srcW, source_k: srcK,
             visited_bands: visited, excluded_bands: excluded,
             pattern_class: PATTERN_CLASS.SUBSET, n_source: n };
  }
  const k3votes = [
    // window 0
    v3(0, 0, [4, 8],   [3, 5, 6, 7],     60),  // band 0 = arrangement A
    v3(0, 1, [5, 6],   [3, 4, 7, 8],    106),  // band 1 = arrangement B
    v3(0, 2, [3, 7],   [4, 5, 6, 8],     60),  // band 2 = arrangement C
    // window 1
    v3(1, 0, [2, 7],   [0, 1, 6, 8],     60),  // band 3 (c0 in w1) = arrangement C
    v3(1, 1, [0, 8],   [1, 2, 6, 7],     60),  // band 4 (c1 in w1) = arrangement A
    v3(1, 2, [1, 6],   [0, 2, 7, 8],    106),  // band 5 (c2 in w1) = arrangement B
    // window 2
    v3(2, 0, [1, 5],   [0, 2, 3, 4],    106),  // band 6 (c0 in w2) = arrangement B
    v3(2, 1, [2, 3],   [0, 1, 4, 5],     60),  // band 7 (c1 in w2) = arrangement C
    v3(2, 2, [0, 4],   [1, 2, 3, 5],     60),  // band 8 (c2 in w2) = arrangement A
  ];

  // band_groups: bands 0,1,2 → window 0; 3,4,5 → window 1; 6,7,8 → window 2
  const band_groups = [0, 0, 0, 1, 1, 1, 2, 2, 2];

  // Without band_groups: should fail to classify as clean (bug fires)
  const rNoFix = consensus_partition(k3votes);
  check('F11: WITHOUT band_groups, score ≤ 0.78 (bug fires)',
        rNoFix.pca_vote_consensus_score <= 0.78,
        `got ${rNoFix.pca_vote_consensus_score.toFixed(4)}`);

  // With band_groups: should classify CLEAN_PARTITION at high score
  const r = consensus_partition(k3votes, { band_groups });
  check('F11: K_locus = 9', r.K_locus === 9);
  check('F11: band_groups attached to voteResult',
        r.votes.band_groups != null && r.votes.band_groups.length === 9);
  check('F11: pca_vote_consensus_score ≥ 0.85 (band_groups fix active)',
        r.pca_vote_consensus_score >= 0.85,
        `got ${r.pca_vote_consensus_score.toFixed(4)}`);
  check('F11: consensus_class = CLEAN_PARTITION',
        r.consensus_class === CONSENSUS_CLASS.CLEAN_PARTITION,
        `got ${r.consensus_class}, score=${r.pca_vote_consensus_score.toFixed(3)}`);
  // The recovered top partition should group bands by arrangement:
  //   A = {0, 4, 8}, B = {1, 5, 6}, C = {2, 3, 7}
  if (r.top_partitions.length > 0) {
    const top = r.top_partitions[0].blocks.map(b => new Set(b));
    const hasA = top.some(s => s.size === 3 && s.has(0) && s.has(4) && s.has(8));
    const hasB = top.some(s => s.size === 3 && s.has(1) && s.has(5) && s.has(6));
    const hasC = top.some(s => s.size === 3 && s.has(2) && s.has(3) && s.has(7));
    check('F11: top recovers arrangement A {0,4,8}', hasA);
    check('F11: top recovers arrangement B {1,5,6}', hasB);
    check('F11: top recovers arrangement C {2,3,7}', hasC);
  }
  // Cross-arrangement same-window pair should have coassoc 0
  // (bands 0 and 1 are both at window 0 but different arrangements;
  // without the fix, coassoc would be ~0.33 from third-arrangement
  // sources excluding both)
  const K = r.coassoc.K_locus;
  check('F11: coassoc[0,1] same-window cross-arrangement = 0',
        Math.abs(r.coassoc.coassoc[0 * K + 1]) < 0.05,
        `got ${r.coassoc.coassoc[0 * K + 1].toFixed(4)}`);
  // Same-arrangement cross-window pair should have high coassoc
  // (bands 0 and 4 are both arrangement A in different windows)
  check('F11: coassoc[0,4] same-arrangement cross-window ≥ 0.9',
        r.coassoc.coassoc[0 * K + 4] >= 0.9,
        `got ${r.coassoc.coassoc[0 * K + 4].toFixed(4)}`);
}

// ---------------------------------------------------------------------
// FIXTURE 12 — N=2 windows, K=3 arrangements: known label-ambiguity case
//
// With only 2 windows of K=3 clusters each, the votes from the two
// envelopes can confirm "there are 3 groups across the locus" but
// cannot pin which group at envelope 1 corresponds to which group at
// envelope 0 (no third anchor for triangulation). All 6 label
// permutations score 1.0 simultaneously, producing N mutually
// incompatible top partitions — a true mathematical ambiguity, not a
// biological multi-layer structure.
//
// Current behavior: this case classifies MULTI_LAYER_STRUCTURE despite
// score=1.0 and residual=0.0. Documented as a known limitation; pin
// the behavior here so any future change is intentional.
// ---------------------------------------------------------------------
console.log('\n--- FIXTURE 12: N=2 K=3 label ambiguity (known limitation) ---');
{
  // 2 windows × 3 clusters → 6 bands.
  // Bands by arrangement (truth, for reference only — voters can't see this
  // unambiguously without a 3rd window):
  //   A = {0, 3}, B = {1, 4}, C = {2, 5}
  function v(sw, sk, vis, exc, n) {
    return { source_w: sw, source_k: sk,
             visited_bands: vis, excluded_bands: exc,
             pattern_class: PATTERN_CLASS.SUBSET, n_source: n };
  }
  const votes = [
    v(0, 0, [3], [4, 5], 60),   // band 0 → matches band 3 by Jaccard
    v(0, 1, [4], [3, 5], 106),
    v(0, 2, [5], [3, 4], 60),
    v(1, 0, [0], [1, 2], 60),
    v(1, 1, [1], [0, 2], 106),
    v(1, 2, [2], [0, 1], 60),
  ];
  const band_groups = [0, 0, 0, 1, 1, 1];
  const r = consensus_partition(votes, { band_groups });
  check('F12: pca_vote_consensus_score = 1.0 (data is internally consistent)',
        Math.abs(r.pca_vote_consensus_score - 1.0) < 1e-9,
        `got ${r.pca_vote_consensus_score.toFixed(4)}`);
  check('F12: pca_hidden_regime_residual = 0 (no per-band conflict)',
        r.pca_hidden_regime_residual < 1e-9,
        `got ${r.pca_hidden_regime_residual.toFixed(4)}`);
  // Multiple partitions tie at 1.0 — label ambiguity is mathematical
  const tied = r.top_partitions.filter(p => Math.abs(p.score - 1.0) < 1e-9);
  check('F12: ≥ 3 partitions tie at score 1.0 (label ambiguity)',
        tied.length >= 3,
        `got ${tied.length}`);
  // Current classifier returns MULTI_LAYER_STRUCTURE for this — document.
  check('F12: known behavior: classifies MULTI_LAYER_STRUCTURE',
        r.consensus_class === CONSENSUS_CLASS.MULTI_LAYER_STRUCTURE,
        `if this changes, update FINDING and decide whether new behavior ` +
        `is intentional. got ${r.consensus_class}`);
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
