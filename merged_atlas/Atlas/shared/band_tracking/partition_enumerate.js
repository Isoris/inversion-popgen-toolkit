// shared/band_tracking/partition_enumerate.js
// =====================================================================
// LAYER B math — bruteforce enumeration and scoring of all set-partitions
// of K target bands.
//
// Bell numbers (number of set partitions of [K]):
//   K=1: 1     K=2: 2     K=3: 5      K=4: 15
//   K=5: 52    K=6: 203   K=7: 877    K=8: 4140
//   K=9: 21147 K=10: 115975  K=11: 678570  K=12: 4213597
//
// Default cap: K=10 (115k partitions, sub-second). Hard fail at K=13.
// Above the cap, callers should fall back to greedy (kt_infer_*).
//
// Each candidate partition P is scored against the co-association matrix
// using pair-based agreement: for every pair (b_i, b_j) in [0..K),
//   if P co-groups them (same block):
//      score gains coassoc[i,j] * total[i,j]
//   else (P puts them in different blocks):
//      score gains (1 - coassoc[i,j]) * total[i,j]
// summed over all pairs, normalized by total possible pair-weight.
//
// This is exactly the "fraction of pair-weighted votes that agree with
// this partition" — pca_vote_consensus_score for partition P.
// =====================================================================

// ---------------------------------------------------------------------
// (1) enumerate_partitions
//
// Generate all set-partitions of [0..K-1] using restricted growth strings.
// A restricted growth string a[0..K-1] satisfies a[0]=0 and
// a[i] <= 1 + max(a[0..i-1]). Each RGS bijects to a partition: the
// block of band i is a[i], with blocks in order of first appearance.
// ---------------------------------------------------------------------

/**
 * @param {int} K
 * @returns {Generator<number[]>}  yields each RGS as an Int32Array-like array;
 *                                  reuse-safe (yielded array is owned by caller
 *                                  if they keep a reference, but we always
 *                                  return a fresh slice to avoid aliasing bugs)
 */
function* enumerate_rgs(K) {
  if (K <= 0) { yield []; return; }
  const a = new Int32Array(K);
  // Track running max for each position so we know the bound at i.
  const m = new Int32Array(K);   // m[i] = max(a[0..i])
  let i = K - 1;
  // Initial = [0,0,0,...,0]
  yield Array.from(a);
  while (true) {
    // Increment the rightmost position whose a[i] < m[i-1] + 1
    i = K - 1;
    while (i > 0 && a[i] === m[i - 1] + 1) i--;
    if (i === 0) return;
    a[i]++;
    m[i] = Math.max(m[i - 1], a[i]);
    for (let j = i + 1; j < K; j++) {
      a[j] = 0;
      m[j] = m[j - 1];
    }
    yield Array.from(a);
  }
}

/**
 * Convert an RGS into the explicit partition (array of blocks, each
 * block an array of band indices).
 *
 * @param {number[]} rgs
 * @returns {number[][]}    block lists, sorted by first-appearance band
 */
function rgs_to_blocks(rgs) {
  const K = rgs.length;
  let nBlocks = 0;
  for (const v of rgs) if (v + 1 > nBlocks) nBlocks = v + 1;
  const blocks = [];
  for (let g = 0; g < nBlocks; g++) blocks.push([]);
  for (let i = 0; i < K; i++) blocks[rgs[i]].push(i);
  return blocks;
}

/**
 * Public: enumerate partitions as block arrays.
 *
 * @param {int} K
 * @yields {number[][]}    each partition is an array of blocks
 */
export function* enumerate_partitions_as_blocks(K) {
  for (const rgs of enumerate_rgs(K)) yield rgs_to_blocks(rgs);
}

// ---------------------------------------------------------------------
// (2) score_partition_against_coassoc
//
// Pair-based scoring per the spec:
//
//   total_pair_weight  = Σ_{i<j} total[i,j]
//   agreement          = Σ_{i<j} ((same_block(i,j) ? coassoc[i,j]
//                                                  : 1 - coassoc[i,j])
//                                 * total[i,j])
//   pca_vote_consensus_score = agreement / total_pair_weight
//
// If total_pair_weight == 0, every partition gets 0; return NaN to
// flag insufficient evidence.
// ---------------------------------------------------------------------

/**
 * @param {number[][]} blocks       a partition of [0..K-1]
 * @param {object}     coassoc      output of build_coassociation_matrix
 * @returns {{
 *   score:        number,           // pca_vote_consensus_score for THIS partition
 *   agreements:   number,           // numerator (weighted)
 *   disagreements: number,
 *   total_pair_weight: number,      // denominator
 * }}
 */
export function score_partition_against_coassoc(blocks, coassoc) {
  const K = coassoc.K_locus;
  // Map each band to its block id
  const blockOf = new Int32Array(K);
  for (let g = 0; g < blocks.length; g++) {
    for (const b of blocks[g]) blockOf[b] = g;
  }
  let agree = 0, disagree = 0, totalPair = 0;
  for (let i = 0; i < K; i++) {
    for (let j = i + 1; j < K; j++) {
      const idx = i * K + j;
      const t = coassoc.total[idx];
      if (t === 0) continue;
      totalPair += t;
      const sameBlock = blockOf[i] === blockOf[j];
      const c = coassoc.coassoc[idx];
      // P agrees with the votes about pair (i,j) iff
      //   sameBlock AND votes co-grouped them (high coassoc), or
      //   !sameBlock AND votes split them (low coassoc).
      // Per-pair contribution to the agreement = (sameBlock ? c : 1-c) * t.
      const contrib = (sameBlock ? c : 1 - c) * t;
      agree += contrib;
      disagree += t - contrib;
    }
  }
  return {
    score: totalPair > 0 ? agree / totalPair : NaN,
    agreements: agree,
    disagreements: disagree,
    total_pair_weight: totalPair,
  };
}

// ---------------------------------------------------------------------
// (3) enumerate_and_score_all_partitions
// ---------------------------------------------------------------------

const DEFAULT_K_CAP = 10;
const HARD_K_FAIL   = 13;

/**
 * Enumerate every partition of K bands and score against coassoc.
 *
 * @param {object} coassoc
 * @param {object} [opts]
 * @param {int}    [opts.K_cap=10]      warn (return null with reason)
 *                                      if K exceeds this
 * @param {int}    [opts.hard_fail=13]  throw if K exceeds this
 * @param {int}    [opts.top_n_cap=200] keep at most this many partitions
 *                                      in the returned `scored` array
 *                                      (top by score; rest discarded)
 * @returns {{
 *   scored:     { blocks: number[][], score: number,
 *                 agreements: number, disagreements: number,
 *                 total_pair_weight: number, rgs: number[] }[],
 *               // sorted descending by score
 *   K:          int,
 *   n_partitions_total: int,
 *   K_cap_exceeded: bool,
 *   total_pair_weight: number,    // shared across all partitions
 * } | null}                       returns null if K_cap exceeded
 *                                  (caller should fall back to greedy)
 */
export function enumerate_and_score_all_partitions(coassoc, opts) {
  opts = opts || {};
  const K_cap     = opts.K_cap     != null ? +opts.K_cap     : DEFAULT_K_CAP;
  const hard_fail = opts.hard_fail != null ? +opts.hard_fail : HARD_K_FAIL;
  const top_n_cap = opts.top_n_cap != null ? +opts.top_n_cap : 200;
  const K = coassoc.K_locus;

  if (K === 0) {
    return { scored: [], K: 0, n_partitions_total: 0, K_cap_exceeded: false,
             total_pair_weight: 0 };
  }
  if (K > hard_fail) {
    throw new Error(`enumerate_and_score_all_partitions: K=${K} exceeds hard_fail=${hard_fail}`);
  }
  if (K > K_cap) {
    return null;   // signal caller to fall back
  }

  // Min-heap keyed by score (we'd want a max-heap of size top_n_cap, but
  // simpler: keep a sorted array of length top_n_cap, insert with binary
  // search; for top_n_cap = 200 this is fine).
  // Actually simpler still: collect all (Bell(K) ≤ 115975 for K=10),
  // sort once at end. Memory cost is acceptable up to K=10.
  const scored = [];
  let totalPairWeight = 0;
  let n_total = 0;
  // Pre-fetch totalPairWeight from a single partition (it's invariant).
  // Compute on the fly for the first partition.
  for (const rgs of enumerate_rgs(K)) {
    n_total++;
    const blocks = rgs_to_blocks(rgs);
    const sc = score_partition_against_coassoc(blocks, coassoc);
    if (n_total === 1) totalPairWeight = sc.total_pair_weight;
    scored.push({ blocks, ...sc, rgs: Array.from(rgs) });
  }
  scored.sort((a, b) => {
    if (Number.isNaN(a.score) && Number.isNaN(b.score)) return 0;
    if (Number.isNaN(a.score)) return 1;
    if (Number.isNaN(b.score)) return -1;
    return b.score - a.score;
  });
  return {
    scored: scored.slice(0, top_n_cap),
    K,
    n_partitions_total: n_total,
    K_cap_exceeded: false,
    total_pair_weight: totalPairWeight,
  };
}

// ---------------------------------------------------------------------
// (4) select_top_partitions_adaptive
//
// Per Q5: return all partitions whose score is within δ of the best.
// If δ=0.10 and best score = 0.85, every partition with score ≥ 0.75 is
// returned. If only one partition is within δ of the best, only one is
// returned (clean single-partition answer).
// ---------------------------------------------------------------------

/**
 * @param {object[]} scored      output[`scored`] from enumerate_and_score_all_partitions
 * @param {object} [opts]
 * @param {number} [opts.delta=0.10]    score window below the best
 * @param {int}    [opts.n_min=1]       always return at least this many
 * @param {int}    [opts.n_max=10]      cap regardless of delta
 * @returns {object[]}
 */
export function select_top_partitions_adaptive(scored, opts) {
  opts = opts || {};
  const delta = opts.delta != null ? +opts.delta : 0.10;
  const n_min = opts.n_min != null ? +opts.n_min : 1;
  const n_max = opts.n_max != null ? +opts.n_max : 10;
  if (scored.length === 0) return [];
  const best = scored[0].score;
  const cutoff = Number.isFinite(best) ? best - delta : -Infinity;
  const picked = [];
  for (const s of scored) {
    if (picked.length >= n_max) break;
    if (s.score >= cutoff || picked.length < n_min) picked.push(s);
    else break;
  }
  return picked;
}

// ---------------------------------------------------------------------
// (5) partitions_are_compatible
//
// Two partitions are "compatible" if one is a refinement of the other,
// i.e., every block of the finer partition is a subset of some block of
// the coarser partition. Compatible partitions describe the SAME layer
// at different granularities (e.g., {0,1,2}|{3,4,5} and
// {0,1,2}|{3,4}|{5}).
//
// INCOMPATIBLE partitions describe DIFFERENT layers (the cross-product
// case: {0,1,2}|{3,4,5} vs {0,3}|{1,4}|{2,5}).
// ---------------------------------------------------------------------

/**
 * @param {number[][]} pA
 * @param {number[][]} pB
 * @returns {{
 *   compatible: bool,
 *   refinement_direction: 'A_finer' | 'B_finer' | 'equal' | 'incompatible',
 * }}
 */
export function partitions_are_compatible(pA, pB) {
  // K is implicit; every band must appear exactly once in each partition.
  const blockOfA = new Map(), blockOfB = new Map();
  for (let g = 0; g < pA.length; g++) for (const b of pA[g]) blockOfA.set(b, g);
  for (let g = 0; g < pB.length; g++) for (const b of pB[g]) blockOfB.set(b, g);
  // Check A finer than B: every A block is contained in one B block.
  let aFinerThanB = true;
  for (const block of pA) {
    if (block.length === 0) continue;
    const bId = blockOfB.get(block[0]);
    for (const x of block) {
      if (blockOfB.get(x) !== bId) { aFinerThanB = false; break; }
    }
    if (!aFinerThanB) break;
  }
  let bFinerThanA = true;
  for (const block of pB) {
    if (block.length === 0) continue;
    const aId = blockOfA.get(block[0]);
    for (const x of block) {
      if (blockOfA.get(x) !== aId) { bFinerThanA = false; break; }
    }
    if (!bFinerThanA) break;
  }
  if (aFinerThanB && bFinerThanA) {
    return { compatible: true, refinement_direction: 'equal' };
  }
  if (aFinerThanB) return { compatible: true, refinement_direction: 'A_finer' };
  if (bFinerThanA) return { compatible: true, refinement_direction: 'B_finer' };
  return { compatible: false, refinement_direction: 'incompatible' };
}

// ---------------------------------------------------------------------
// (6) Bell numbers — exposed for testing and for caller-side
// "should I bother" decisions.
// ---------------------------------------------------------------------

const BELL_NUMBERS = [
  1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597,
];

export function bellNumber(K) {
  if (K < 0) return 0;
  if (K < BELL_NUMBERS.length) return BELL_NUMBERS[K];
  // For K beyond the table, use Bell triangle recurrence
  let row = [1];
  for (let n = 1; n <= K; n++) {
    const next = [row[row.length - 1]];
    for (let k = 0; k < row.length; k++) next.push(next[k] + row[k]);
    row = next;
  }
  return row[0];
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._enumerate_partitions_as_blocks   = enumerate_partitions_as_blocks;
  window._score_partition_against_coassoc  = score_partition_against_coassoc;
  window._enumerate_and_score_all_partitions = enumerate_and_score_all_partitions;
  window._select_top_partitions_adaptive   = select_top_partitions_adaptive;
  window._partitions_are_compatible        = partitions_are_compatible;
  window._bellNumber                       = bellNumber;
}
