// shared/hungarian.js
// =====================================================================
// Hungarian-style label alignment + chain projection + concordance.
//
// Extracted from legacy Inversion_atlas.html (turn 165 close binary).
// All functions are pure.
//
// Source line refs:
//   alignLabels                    line 30843
//   permutations                   line 30879
//   _hungarianChainProjection      line 39120
//   _finalizeChain                 line 39190
//   _concordanceMatrix             line 39212
//   _LINEAGE_CHAIN_BREAK_AGREEMENT line 39090
//
// Note: the Hungarian alignment here uses brute-force K! permutation
// enumeration (correct for K up to ~8; the atlas uses K=3 and K=6 so
// 720 perms max — well within budget). For larger K the canonical
// O(K^3) Hungarian assignment would be needed, but it isn't today.
// =====================================================================

/**
 * Threshold below which the chain projection breaks and a fresh chain
 * is started anchored at the new L2's own labels. Default = 0.50;
 * legacy code at line 39090 used this exact value.
 */
export const LINEAGE_CHAIN_BREAK_AGREEMENT = 0.50;

/**
 * Align two K-clustering label arrays by enumerating K! permutations
 * and picking the one that maximises the diagonal of the contingency
 * table.
 *
 * Returns { table, raw_table, perm, invPerm, aligned, concord, n }:
 *   - table:     K×K contingency after applying the best permutation
 *   - raw_table: K×K contingency before any permutation
 *   - perm:      bestPerm[r] = the g2-label that aligns to g1-label r
 *   - invPerm:   inverse of perm
 *   - aligned:   Int8Array of g2 with permuted labels
 *   - concord:   diagonal-sum / n  (∈ [0, 1])
 *
 * @param {number[]|Int8Array} g1
 * @param {number[]|Int8Array} g2
 * @param {number} K
 */
export function alignLabels(g1, g2, K) {
  const n = g1.length;
  const T = Array.from({ length: K }, () => new Int32Array(K));
  for (let i = 0; i < n; i++) T[g1[i]][g2[i]]++;
  const perms = permutations(K);
  let bestPerm = perms[0], bestSum = -1;
  for (const p of perms) {
    let s = 0;
    for (let r = 0; r < K; r++) s += T[r][p[r]];
    if (s > bestSum) { bestSum = s; bestPerm = p; }
  }
  const invPerm = new Array(K);
  for (let r = 0; r < K; r++) invPerm[bestPerm[r]] = r;
  const aligned = new Int8Array(n);
  for (let i = 0; i < n; i++) aligned[i] = invPerm[g2[i]];
  const Ta = Array.from({ length: K }, () => new Int32Array(K));
  for (let i = 0; i < n; i++) Ta[g1[i]][aligned[i]]++;
  let diag = 0;
  for (let r = 0; r < K; r++) diag += Ta[r][r];
  return {
    table: Ta,
    raw_table: T,
    perm: bestPerm,
    invPerm,
    aligned,
    concord: n > 0 ? diag / n : 0,
    n,
  };
}

/**
 * Generate all K! permutations of [0..K-1]. Stable Heap-style recursion.
 *
 * @param {number} K
 * @returns {number[][]}
 */
export function permutations(K) {
  const out = [];
  const arr = Array.from({ length: K }, (_, i) => i);
  function rec(start) {
    if (start === K) { out.push(arr.slice()); return; }
    for (let i = start; i < K; i++) {
      [arr[start], arr[i]] = [arr[i], arr[start]];
      rec(start + 1);
      [arr[start], arr[i]] = [arr[i], arr[start]];
    }
  }
  rec(0);
  return out;
}

/**
 * Walk a chain of L2 indices, aligning each to its predecessor via
 * `alignLabels`, breaking the chain when concordance drops below
 * `LINEAGE_CHAIN_BREAK_AGREEMENT`. Returns:
 *
 *   { chains, n_samples, n_chains, n_total_L2 }
 *
 * Each chain is { l2_indices, projected (Int8Array, row-major), start_idx, n_L2 }.
 *
 * @param {number[]} l2_indices         Ordered L2 indices to project.
 * @param {number} K                    K used for clustering at every L2.
 * @param {(idx:number) => Int8Array|null} getLabelsForL2
 *        Caller-supplied: returns the per-sample label array for the
 *        given L2 index (or null if unavailable). The legacy default
 *        was a `getL2Cluster(idx).fixedKLabels` lookup; that requires
 *        the discovery sub-atlas's `getL2Cluster` to be in scope, so
 *        the caller passes it explicitly here.
 */
export function hungarianChainProjection(l2_indices, K, getLabelsForL2) {
  if (typeof getLabelsForL2 !== 'function') {
    throw new TypeError('hungarianChainProjection: getLabelsForL2 callback is required');
  }
  const n_L2 = l2_indices.length;
  if (n_L2 === 0) return { chains: [], n_samples: 0, n_chains: 0, n_total_L2: 0 };

  let n_samples = 0;
  let firstUsableIdx = -1;
  let firstLabels = null;
  for (let i = 0; i < n_L2; i++) {
    const lab = getLabelsForL2(l2_indices[i]);
    if (lab && lab.length > 0) {
      firstUsableIdx = i; firstLabels = lab; n_samples = lab.length; break;
    }
  }
  if (firstUsableIdx < 0) {
    return { chains: [], n_samples: 0, n_chains: 0, n_total_L2: 0 };
  }

  const chains = [];
  let curChain = {
    l2_indices: [l2_indices[firstUsableIdx]],
    projectedRows: [firstLabels],
    start_idx: firstUsableIdx,
  };
  let refLabels = firstLabels;

  for (let i = firstUsableIdx + 1; i < n_L2; i++) {
    const lab = getLabelsForL2(l2_indices[i]);
    if (!lab || lab.length !== n_samples) continue;
    let aligned = lab;
    let agreement = 1;
    try {
      const r = alignLabels(refLabels, lab, K);
      if (r && r.aligned && typeof r.concord === 'number') {
        aligned = r.aligned;
        agreement = r.concord;
      }
    } catch (_) { /* fall back to raw labels */ }
    if (agreement < LINEAGE_CHAIN_BREAK_AGREEMENT) {
      chains.push(finalizeChain(curChain, n_samples));
      curChain = {
        l2_indices: [l2_indices[i]],
        projectedRows: [lab],
        start_idx: i,
      };
      refLabels = lab;
    } else {
      curChain.l2_indices.push(l2_indices[i]);
      curChain.projectedRows.push(aligned);
      refLabels = aligned;
    }
  }
  chains.push(finalizeChain(curChain, n_samples));

  let n_total_L2 = 0;
  for (const ch of chains) n_total_L2 += ch.l2_indices.length;

  return { chains, n_samples, n_chains: chains.length, n_total_L2 };
}

function finalizeChain(chain, n_samples) {
  const n = chain.l2_indices.length;
  const flat = new Int8Array(n * n_samples);
  for (let row = 0; row < n; row++) {
    const src = chain.projectedRows[row];
    for (let s = 0; s < n_samples; s++) flat[row * n_samples + s] = src[s];
  }
  return {
    l2_indices: chain.l2_indices,
    projected:  flat,
    start_idx:  chain.start_idx,
    n_L2:       n,
  };
}

/**
 * Per-fish-pair concordance matrix from a chain projection.
 *
 *   concordance(i, j) = (# L2s where projected[i] == projected[j]) /
 *                       (# L2s where both labels are >= 0)
 *
 * Returns Float32Array of length n_samples^2 (full square; symmetric;
 * diagonal = 1). Storage is square (not packed) so callers can index
 * without translation.
 *
 * @param {{chains, n_samples, n_total_L2}} projection  Output of hungarianChainProjection
 * @returns {Float32Array}
 */
export function concordanceMatrix(projection) {
  const nS = projection.n_samples;
  const out = new Float32Array(nS * nS);
  for (let i = 0; i < nS; i++) out[i * nS + i] = 1;
  if (nS === 0 || projection.n_total_L2 === 0) return out;

  const matchCount = new Uint32Array(nS * nS);
  const validCount = new Uint32Array(nS * nS);
  for (const ch of projection.chains) {
    const proj = ch.projected;
    const nL2 = ch.n_L2;
    for (let row = 0; row < nL2; row++) {
      const base = row * nS;
      for (let i = 0; i < nS; i++) {
        const li = proj[base + i];
        if (li < 0) continue;
        for (let j = i + 1; j < nS; j++) {
          const lj = proj[base + j];
          if (lj < 0) continue;
          const key = i * nS + j;
          validCount[key]++;
          if (li === lj) matchCount[key]++;
        }
      }
    }
  }
  for (let i = 0; i < nS; i++) {
    for (let j = i + 1; j < nS; j++) {
      const key = i * nS + j;
      const v = validCount[key];
      const c = (v > 0) ? (matchCount[key] / v) : 0;
      out[key] = c;
      out[j * nS + i] = c;
    }
  }
  return out;
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window.alignLabels                 = alignLabels;
  window.permutations                = permutations;
  window._hungarianChainProjection   = hungarianChainProjection;
  window._concordanceMatrix          = concordanceMatrix;
  window._LINEAGE_CHAIN_BREAK_AGREEMENT = LINEAGE_CHAIN_BREAK_AGREEMENT;
}
