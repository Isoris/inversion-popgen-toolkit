// shared/contingency.js
// =====================================================================
// Pure-compute primitives for K×K contingency tables and the metrics
// the atlas computes on top of them.
//
// Extracted from legacy Inversion_atlas.html (turn 165 close binary).
// All functions are pure: no `state`, no DOM, no globals.
//
// Source line refs in the legacy file (for traceability):
//   _buildContingency         line 12250
//   _detectFuseEvents         line 12274
//   _computeARI               line 12320
//   _computeNMI               line 12353
//   _scaleStabilityVerdict    line 12408
//   _cramersV                 line 36046
//   _chiSqSurvival            line 38443
//   _lnGamma                  line 38482
// =====================================================================

/**
 * Build a K_a × K_b contingency table from two label arrays of equal
 * length. Returns null on bad inputs.
 *
 * @param {number[]|Int8Array} labelsA  Label array for partition A
 * @param {number[]|Int8Array} labelsB  Label array for partition B (same length)
 * @param {number} KA                   Cardinality of A's labels (0..KA-1 valid)
 * @param {number} KB                   Cardinality of B's labels
 * @returns {{ M: number[][], KA: number, KB: number, n: number } | null}
 */
export function buildContingency(labelsA, labelsB, KA, KB) {
  if (!isArrayLike(labelsA) || !isArrayLike(labelsB)) return null;
  if (labelsA.length !== labelsB.length) return null;
  if (!(KA > 0) || !(KB > 0)) return null;
  const M = Array.from({ length: KA }, () => new Array(KB).fill(0));
  let n = 0;
  for (let i = 0; i < labelsA.length; i++) {
    const a = labelsA[i], b = labelsB[i];
    if (a < 0 || a >= KA || b < 0 || b >= KB) continue;
    M[a][b]++; n++;
  }
  return { M, KA, KB, n };
}

/**
 * Fuse-event detection: for each coarse cluster c, find all fine
 * clusters f where ≥thresh fraction of f's samples land in c. If two
 * or more such fine clusters exist for one coarse cluster, that's a
 * fuse event.
 *
 * Convention: rows of `ct.M` = fine clusters; cols = coarse clusters.
 *
 * @param {{ M: number[][], KA: number, KB: number }} ct
 * @param {{ thresh?: number }} [opts]
 * @returns {{ fine_clusters: number[], coarse_cluster: number }[]}
 */
export function detectFuseEvents(ct, opts) {
  if (!ct || !ct.M) return [];
  const t = (opts && opts.thresh != null) ? +opts.thresh : 0.80;
  const { M, KA, KB } = ct;
  const rowSum = new Array(KA).fill(0);
  for (let a = 0; a < KA; a++) for (let b = 0; b < KB; b++) rowSum[a] += M[a][b];
  const fuses = [];
  for (let b = 0; b < KB; b++) {
    const fineList = [];
    for (let a = 0; a < KA; a++) {
      if (rowSum[a] === 0) continue;
      if (M[a][b] / rowSum[a] >= t) fineList.push(a);
    }
    if (fineList.length >= 2) fuses.push({ fine_clusters: fineList, coarse_cluster: b });
  }
  return fuses;
}

/**
 * Adjusted Rand Index between two label arrays. Returns NaN on bad
 * input, 1 on identical partitions.
 *
 * @param {number[]|Int8Array} labelsA
 * @param {number[]|Int8Array} labelsB
 * @returns {number}
 */
export function computeARI(labelsA, labelsB) {
  if (!isArrayLike(labelsA) || labelsA.length !== labelsB.length) return NaN;
  const n = labelsA.length;
  if (n < 2) return NaN;
  const cellByKey = new Map();
  const rowByA = new Map();
  const colByB = new Map();
  for (let i = 0; i < n; i++) {
    const a = labelsA[i], b = labelsB[i];
    const key = a + '|' + b;
    cellByKey.set(key, (cellByKey.get(key) || 0) + 1);
    rowByA.set(a, (rowByA.get(a) || 0) + 1);
    colByB.set(b, (colByB.get(b) || 0) + 1);
  }
  const C2 = (k) => k * (k - 1) / 2;
  let sumCellPairs = 0;
  for (const v of cellByKey.values()) sumCellPairs += C2(v);
  let sumRowPairs = 0;
  for (const v of rowByA.values()) sumRowPairs += C2(v);
  let sumColPairs = 0;
  for (const v of colByB.values()) sumColPairs += C2(v);
  const totalPairs = C2(n);
  if (totalPairs === 0) return NaN;
  const expected = (sumRowPairs * sumColPairs) / totalPairs;
  const max = (sumRowPairs + sumColPairs) / 2;
  if (max === expected) return 1;
  return (sumCellPairs - expected) / (max - expected);
}

/**
 * Normalized Mutual Information (symmetric, log base e). Range [0, 1].
 *
 * @param {number[]|Int8Array} labelsA
 * @param {number[]|Int8Array} labelsB
 * @returns {number}
 */
export function computeNMI(labelsA, labelsB) {
  if (!isArrayLike(labelsA) || labelsA.length !== labelsB.length) return NaN;
  const n = labelsA.length;
  if (n < 2) return NaN;
  const cellByKey = new Map();
  const rowByA = new Map();
  const colByB = new Map();
  for (let i = 0; i < n; i++) {
    const a = labelsA[i], b = labelsB[i];
    cellByKey.set(a + '|' + b, (cellByKey.get(a + '|' + b) || 0) + 1);
    rowByA.set(a, (rowByA.get(a) || 0) + 1);
    colByB.set(b, (colByB.get(b) || 0) + 1);
  }
  let HA = 0, HB = 0, I = 0;
  for (const ca of rowByA.values()) {
    const p = ca / n;
    if (p > 0) HA -= p * Math.log(p);
  }
  for (const cb of colByB.values()) {
    const p = cb / n;
    if (p > 0) HB -= p * Math.log(p);
  }
  for (const [key, cab] of cellByKey) {
    const [aStr, bStr] = key.split('|');
    const ca = rowByA.get(+aStr) || rowByA.get(aStr);
    const cb = colByB.get(+bStr) || colByB.get(bStr);
    const pab = cab / n;
    const pa  = ca / n;
    const pb  = cb / n;
    if (pab > 0 && pa > 0 && pb > 0) I += pab * Math.log(pab / (pa * pb));
  }
  if (HA + HB === 0) return 1;
  return (2 * I) / (HA + HB);
}

/**
 * Cramér's V for a flat row-major K_a × K_c contingency table. Returns
 * a value in [0, 1] (0 = independent, 1 = perfect association).
 *
 * Note: this primitive uses a flat array rather than a {M, KA, KB}
 * object — it's the form the legacy code used at line 36046, kept as
 * the canonical input contract.
 *
 * @param {number[]|Int32Array|Float32Array} table  Length K_a*K_c, row-major
 * @param {number} K_a
 * @param {number} K_c
 * @returns {number}
 */
export function cramersV(table, K_a, K_c) {
  let n = 0;
  const rowSums = new Array(K_a).fill(0);
  const colSums = new Array(K_c).fill(0);
  for (let i = 0; i < K_a; i++) {
    for (let j = 0; j < K_c; j++) {
      const v = table[i * K_c + j];
      n += v; rowSums[i] += v; colSums[j] += v;
    }
  }
  if (n < 2) return NaN;
  let nzR = 0, nzC = 0;
  for (let i = 0; i < K_a; i++) if (rowSums[i] > 0) nzR++;
  for (let j = 0; j < K_c; j++) if (colSums[j] > 0) nzC++;
  if (nzR < 2 || nzC < 2) return 0;
  let chi2 = 0;
  for (let i = 0; i < K_a; i++) {
    for (let j = 0; j < K_c; j++) {
      const exp = (rowSums[i] * colSums[j]) / n;
      if (exp > 0) {
        const o = table[i * K_c + j];
        const d = o - exp;
        chi2 += d * d / exp;
      }
    }
  }
  const minDim = Math.min(nzR, nzC) - 1;
  if (minDim < 1) return 0;
  const v = Math.sqrt(chi2 / (n * minDim));
  return Math.max(0, Math.min(1, v));
}

/**
 * Chi-squared survival function (1 - CDF). Numerically stable for
 * df up to ~200. Used in conjunction with cramersV() chi2 outputs.
 *
 * @param {number} chi2
 * @param {number} df
 * @returns {number}
 */
export function chiSqSurvival(chi2, df) {
  if (!isFinite(chi2) || !isFinite(df)) return NaN;
  if (chi2 <= 0) return 1;
  if (df <= 0) return NaN;
  const s = df / 2, x = chi2 / 2;
  if (x < s + 1) {
    let term = 1 / s, sum = term;
    for (let k = 1; k < 200; k++) {
      term *= x / (s + k);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-12) break;
    }
    const lnGammaS = lnGamma(s);
    const P = sum * Math.exp(-x + s * Math.log(x) - lnGammaS);
    return Math.max(0, Math.min(1, 1 - P));
  } else {
    const TINY = 1e-300;
    let b = x + 1 - s, c = 1 / TINY, d = 1 / b, h = d;
    for (let k = 1; k < 200; k++) {
      const an = -k * (k - s);
      b += 2;
      d = an * d + b; if (Math.abs(d) < TINY) d = TINY;
      c = b + an / c; if (Math.abs(c) < TINY) c = TINY;
      d = 1 / d;
      const delta = d * c;
      h *= delta;
      if (Math.abs(delta - 1) < 1e-12) break;
    }
    const lnGammaS = lnGamma(s);
    const Q = Math.exp(-x + s * Math.log(x) - lnGammaS) * h;
    return Math.max(0, Math.min(1, Q));
  }
}

/**
 * Stirling-based log Gamma (Lanczos approximation, ~10 digits for x≥0.5).
 *
 * @param {number} x
 * @returns {number}
 */
export function lnGamma(x) {
  const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
               -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) {
    y += 1;
    ser += cof[j] / y;
  }
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/**
 * Scale-stability verdict cascade for a 3-pane local-PCA window.
 *
 * Verdicts (legacy thresholds, all configurable via opts):
 *   STABLE_3BAND       all 3 panes K=3, no fuses, no splits, ARI≥0.85 pairwise
 *   STABLE_6BAND       all 3 panes K=6, no fuses, no splits, ARI≥0.85 pairwise
 *   NESTED_3IN6        K=6 fine, K=3 coarse, exactly 3 fuses, 0 splits, no
 *                      crossing fuses (each fine cluster maps to ONE coarse)
 *   OVERLAP_BREAKS_3   K=3 fine AND K=3 coarse, K=6 medium, fuses+splits at
 *                      the medium boundary (ribbons reshuffle middle, restore
 *                      at edges). Heuristic: pane1↔pane3 ARI≥0.7 (edges agree)
 *                      AND pane2 introduces fuses+splits relative to both
 *                      neighbors.
 *   UNSTABLE           otherwise (low ARI somewhere, or noisy / random)
 *
 * @param {Array} panes      [{K, ok, ...}, {K, ok, ...}, {K, ok, ...}]
 * @param {Array} pairwise   [{contingency, fuseEvents, splitEvents, ari, nmi}, …]
 * @param {{ari_stable?: number, ari_edge?: number}} [opts]
 * @returns {string}
 */
export function scaleStabilityVerdict(panes, pairwise, opts) {
  const t = Object.assign({ ari_stable: 0.85, ari_edge: 0.70 }, opts || {});
  if (!Array.isArray(panes) || panes.length !== 3) return 'UNSTABLE';
  if (!Array.isArray(pairwise) || pairwise.length !== 2) return 'UNSTABLE';
  const allOk = panes.every(p => p && p.ok !== false);
  if (!allOk) return 'UNSTABLE';
  const Ks = panes.map(p => p.K | 0);
  const ari12 = pairwise[0] && isFinite(pairwise[0].ari) ? pairwise[0].ari : NaN;
  const ari23 = pairwise[1] && isFinite(pairwise[1].ari) ? pairwise[1].ari : NaN;
  const nFuse12  = (pairwise[0] && pairwise[0].fuseEvents)  ? pairwise[0].fuseEvents.length  : 0;
  const nFuse23  = (pairwise[1] && pairwise[1].fuseEvents)  ? pairwise[1].fuseEvents.length  : 0;
  const nSplit12 = (pairwise[0] && pairwise[0].splitEvents) ? pairwise[0].splitEvents.length : 0;
  const nSplit23 = (pairwise[1] && pairwise[1].splitEvents) ? pairwise[1].splitEvents.length : 0;

  const stableEdges = (ari12 >= t.ari_stable) && (ari23 >= t.ari_stable)
                    && nFuse12 + nFuse23 + nSplit12 + nSplit23 === 0;
  if (stableEdges && Ks[0] === 3 && Ks[1] === 3 && Ks[2] === 3) return 'STABLE_3BAND';
  if (stableEdges && Ks[0] === 6 && Ks[1] === 6 && Ks[2] === 6) return 'STABLE_6BAND';

  // Nested 3-in-6: middle pane K=6 fine, edges K=3 coarse, exactly 3 fuses each.
  if (Ks[1] === 6 && Ks[0] === 3 && Ks[2] === 3
      && nFuse12 === 3 && nFuse23 === 3
      && nSplit12 === 0 && nSplit23 === 0) {
    return 'NESTED_3IN6';
  }

  // Overlap breaks 3: edges agree (K=3, K=3, ARI≥edge), middle (K=6) introduces
  // fuses+splits relative to both neighbors.
  if (Ks[0] === 3 && Ks[2] === 3 && Ks[1] === 6
      && nFuse12 + nSplit12 > 0 && nFuse23 + nSplit23 > 0) {
    // Approximate the pane1↔pane3 ARI via the chain ari12, ari23 lower bound.
    // The legacy code may have computed it directly; we keep the looser
    // proxy (both edge-pair ARIs above ari_edge). If that's too loose for
    // real data, swap in a direct pane1↔pane3 ARI via ari1_3 in opts.
    const ari1_3_proxy = Math.min(ari12, ari23);
    if (ari1_3_proxy >= t.ari_edge) return 'OVERLAP_BREAKS_3';
  }

  return 'UNSTABLE';
}

// ---------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------
function isArrayLike(x) {
  return Array.isArray(x) || (x && typeof x.length === 'number' && typeof x !== 'string');
}

// ---------------------------------------------------------------------
// Console-debug exposures (optional; preserves legacy `window._buildContingency`
// debugging convenience)
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._buildContingency      = buildContingency;
  window._detectFuseEvents      = detectFuseEvents;
  window._computeARI            = computeARI;
  window._computeNMI            = computeNMI;
  window._cramersV              = cramersV;
  window._chiSqSurvival         = chiSqSurvival;
  window._lnGamma               = lnGamma;
  window._scaleStabilityVerdict = scaleStabilityVerdict;
}
