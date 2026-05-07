// shared/per_l2_cluster.js
// =====================================================================
// Per-L2 PC aggregation + clustering, refactored to take an explicit
// context object instead of reading a global `state`.
//
// Source line refs (legacy Inversion_atlas.html, turn 165 close):
//   getPC                line 9951   (wrapped here as ctx.getPC)
//   aggregateL2          line 10249
//   sampleSpreadL2       line 10294  (-> sampleSpreadRange via env._s0/._e0)
//   sampleSpreadRange    line 10304
//   sigmaProfileL2       line 10342
//   clusterL2            line 10503
//   clusterL2AtK         line 10694
//   ensureL2Cache        line 10665
//   getL2Cluster         line 10673
//   getL2ClusterAt       line 10770
//
// Design: the legacy code reads `state.data`, `state.k`, `state.aggMethod`,
// `state.kMode`, `state.kRange`, `state.silThreshold`, `state.minNGroup`,
// `state.minNWin`, `state.flipPC1`, `state.pc1Sign` directly. The new
// API takes a `ctx` object that bundles those reads as fields/callbacks.
// This is the pattern from `hungarianChainProjection`: pure function +
// caller-supplied accessors. No global state, no caches.
//
// Caller manages caching (typically a Map keyed by l2idx + ctx.cacheKey).
// =====================================================================

import { kmeans1D, kmeans2D, adaptiveK1D } from './kmeans.js';

/**
 * Build a clustering context object from a state-like input. Helper
 * for callers to keep the call site short.
 *
 * @param {object} state         the discovery sub-atlas's state
 * @returns {ClusterContext}
 *
 * @typedef {object} ClusterContext
 * @property {object} data           state.data
 * @property {number} k              state.k
 * @property {string} aggMethod      'mean_pc1' | 'mean_pc12' | 'median_pc1'
 * @property {string} kMode          'fixed' | 'adaptive'
 * @property {number[]} kRange       [kMin, kMax] for adaptive
 * @property {number} silThreshold
 * @property {number} minNGroup
 * @property {number} minNWin
 * @property {(winIdx:number) => {pc1:Array, pc2:Array, sign:number}} getPC
 *           callback returning per-window PC arrays + sign correction
 *           (caller wraps state.data.windows[w] + state.pc1Sign + state.flipPC1)
 */
export function contextFromState(state) {
  return {
    data: state.data,
    k: state.k,
    aggMethod: state.aggMethod,
    kMode: state.kMode,
    kRange: state.kRange,
    silThreshold: state.silThreshold,
    minNGroup: state.minNGroup,
    minNWin: state.minNWin,
    getPC: (winIdx) => {
      const w = state.data.windows[winIdx];
      const s = state.flipPC1 && state.pc1Sign ? state.pc1Sign[winIdx] : 1;
      return { pc1: w.pc1, pc2: w.pc2, sign: s };
    },
  };
}

// =====================================================================
// 1. Per-L2 aggregation (mean / median PC1 across windows in the L2)
// =====================================================================

/**
 * Aggregate a sample's PC1 (and optionally PC2) across all windows in
 * an L2 envelope. Returns per-sample mean/median values to be fed
 * into K-means.
 *
 * @param {ClusterContext} ctx
 * @param {number} l2idx          L2 envelope index in ctx.data.l2_envelopes
 * @returns {{xs: Float64Array, ys: Float64Array|null, nW: number} | null}
 */
export function aggregateL2(ctx, l2idx) {
  const d = ctx.data;
  const env = d.l2_envelopes[l2idx];
  const s = env._s0, e = env._e0;
  const nW = e - s + 1;
  const nS = d.n_samples;
  if (nW <= 0) return null;

  const xs = new Float64Array(nS);
  const ys = ctx.aggMethod === 'mean_pc12' ? new Float64Array(nS) : null;

  if (ctx.aggMethod === 'median_pc1') {
    const tmp = new Float64Array(nW);
    for (let si = 0; si < nS; si++) {
      for (let w = 0; w < nW; w++) {
        const { pc1, sign } = ctx.getPC(s + w);
        tmp[w] = pc1[si] * sign;
      }
      const sorted = Array.from(tmp).sort((a, b) => a - b);
      xs[si] = sorted.length % 2 === 1
        ? sorted[(sorted.length - 1) >> 1]
        : 0.5 * (sorted[sorted.length / 2 - 1] + sorted[sorted.length / 2]);
    }
  } else {
    for (let w = 0; w < nW; w++) {
      const { pc1, pc2, sign } = ctx.getPC(s + w);
      for (let si = 0; si < nS; si++) {
        xs[si] += pc1[si] * sign;
        if (ys) ys[si] += pc2[si];
      }
    }
    for (let si = 0; si < nS; si++) {
      xs[si] /= nW;
      if (ys) ys[si] /= nW;
    }
  }
  return { xs, ys, nW };
}

// =====================================================================
// 2. Per-sample spread (std-dev of PC1 across the L2's windows)
// =====================================================================

/**
 * Per-sample standard deviation of sign-aligned PC1 across an arbitrary
 * window range [s, e] (both inclusive, 0-based). Used by sampleSpreadL2
 * (full L2 range) and by the candidate-page σ profile (across multiple L2s).
 *
 * @param {ClusterContext} ctx
 * @param {number} s            first window (inclusive)
 * @param {number} e            last window (inclusive)
 * @returns {Float64Array | null}  per-sample σ; null on input <2 windows
 */
export function sampleSpreadRange(ctx, s, e) {
  const d = ctx.data;
  if (!d) return null;
  const nW = e - s + 1;
  const nS = d.n_samples;
  if (nW < 2) return null;
  const mean = new Float64Array(nS);
  const sumSq = new Float64Array(nS);
  for (let w = 0; w < nW; w++) {
    const { pc1, sign } = ctx.getPC(s + w);
    for (let si = 0; si < nS; si++) mean[si] += pc1[si] * sign;
  }
  for (let si = 0; si < nS; si++) mean[si] /= nW;
  for (let w = 0; w < nW; w++) {
    const { pc1, sign } = ctx.getPC(s + w);
    for (let si = 0; si < nS; si++) {
      const v = pc1[si] * sign - mean[si];
      sumSq[si] += v * v;
    }
  }
  const sd = new Float64Array(nS);
  for (let si = 0; si < nS; si++) sd[si] = Math.sqrt(sumSq[si] / (nW - 1));
  return sd;
}

/**
 * Per-sample σ across an L2 envelope's full window range.
 *
 * @param {ClusterContext} ctx
 * @param {number} l2idx
 * @returns {Float64Array | null}
 */
export function sampleSpreadL2(ctx, l2idx) {
  const d = ctx.data;
  const env = d.l2_envelopes[l2idx];
  if (!env) return null;
  return sampleSpreadRange(ctx, env._s0, env._e0);
}

// =====================================================================
// 3. σ profile + bimodality verdict (legacy `sigmaProfileL2`)
// =====================================================================

/**
 * Test the "two inversions vs double crossovers" hypothesis for an L2.
 *
 * Returns null on insufficient data. Otherwise:
 *   { q50, q90, q95, n_high, ratio_high, is_bimodal,
 *     bimodality_coef, verdict, reason }
 *
 * @param {ClusterContext} ctx
 * @param {number} l2idx
 * @param {number} usedK   the K actually used by clusterL2 (may differ from ctx.k)
 */
export function sigmaProfileL2(ctx, l2idx, usedK) {
  const sd = sampleSpreadL2(ctx, l2idx);
  if (!sd) return null;
  const vals = Array.from(sd).filter(v => isFinite(v));
  if (vals.length < 10) return null;
  vals.sort((a, b) => a - b);
  const q = (p) => vals[Math.min(vals.length - 1, Math.floor(p * vals.length))];
  const q50 = q(0.50), q90 = q(0.90), q95 = q(0.95);

  // High-σ tail: samples whose σ exceeds 2× q50 (heuristic from D16)
  const high_thr = q50 * 2;
  let n_high = 0;
  for (const v of vals) if (v > high_thr) n_high++;
  const ratio_high = n_high / vals.length;

  // Sarle's bimodality coefficient: (skew^2 + 1) / kurt
  const n = vals.length;
  let mean = 0;
  for (const v of vals) mean += v;
  mean /= n;
  let m2 = 0, m3 = 0, m4 = 0;
  for (const v of vals) {
    const dv = v - mean;
    m2 += dv * dv; m3 += dv * dv * dv; m4 += dv * dv * dv * dv;
  }
  m2 /= n; m3 /= n; m4 /= n;
  const sd_pop = Math.sqrt(m2);
  const skew = sd_pop > 1e-12 ? m3 / (sd_pop * sd_pop * sd_pop) : 0;
  const kurt = m2 > 1e-12 ? m4 / (m2 * m2) : 0;
  const bimodality_coef = kurt > 1e-12 ? (skew * skew + 1) / kurt : 0;
  const is_bimodal = bimodality_coef > (5 / 9) && ratio_high < 0.25 && ratio_high > 0.02;

  let verdict, reason;
  if (is_bimodal && usedK > 3) {
    verdict = 'DOUBLE_CROSSOVER_LIKELY';
    reason = `bimodal σ + small high-σ tail + K=${usedK}>3 — high-σ samples may be het samples appearing as hom inside converted segments`;
  } else if (q50 < 0.05 && usedK > 3) {
    verdict = 'STACKED_INVERSIONS';
    reason = `σ low for all samples + K=${usedK}>3 — joint karyotypes appear stable`;
  } else if (q90 > 0.5) {
    verdict = 'NOISY';
    reason = `q90=${q90.toFixed(3)} too high — region is noisy or low power`;
  } else {
    verdict = 'NORMAL';
    reason = `single inversion or low complexity (q50=${q50.toFixed(3)}, K=${usedK})`;
  }
  return { q50, q90, q95, n_high, ratio_high, is_bimodal, bimodality_coef, verdict, reason };
}

// =====================================================================
// 4. clusterL2 — the main per-L2 K-means pipeline
// =====================================================================

/**
 * Cluster an L2 envelope: aggregate PC1 (mean/median or 2D), run K-means
 * (fixed or adaptive), compute family purity + halves coherence.
 *
 * @param {ClusterContext} ctx
 * @param {number} l2idx
 * @returns {object}  result with shape:
 *     { ok, reason, labels, n_per_group, centers, centers_y, nW,
 *       n_below_threshold, fam_purity, fam_per_cluster,
 *       coherence, incoherent, usedK, silhouette, fixedKLabels }
 */
export function clusterL2(ctx, l2idx) {
  const d = ctx.data;
  const env = d.l2_envelopes[l2idx];
  const nW = env._e0 - env._s0 + 1;
  if (nW < 1) return { ok: false, reason: 'NO_WINDOWS' };

  const agg = aggregateL2(ctx, l2idx);
  if (!agg) return { ok: false, reason: 'NO_WINDOWS' };

  let result;
  let usedK = ctx.k;
  let silhouette = null;
  let fixedKLabels = null;
  let fixedKResult = null;
  if (ctx.kMode === 'adaptive' && ctx.aggMethod !== 'mean_pc12') {
    const ak = adaptiveK1D(agg.xs, ctx.kRange[0], ctx.kRange[1], ctx.silThreshold, ctx.minNGroup);
    if (ak != null) {
      result = { labels: ak.labels, centers: ak.centers, n_per_group: ak.n_per_group };
      usedK = ak.k;
      silhouette = ak.silhouette;
    } else {
      result = kmeans1D(agg.xs, ctx.k);
      usedK = ctx.k;
    }
    fixedKResult = kmeans1D(agg.xs, ctx.k);
    fixedKLabels = fixedKResult.labels;
  } else if (ctx.aggMethod === 'mean_pc12') {
    result = kmeans2D(agg.xs, agg.ys, ctx.k);
    usedK = ctx.k;
    fixedKLabels = result.labels;
  } else {
    result = kmeans1D(agg.xs, ctx.k);
    usedK = ctx.k;
    fixedKLabels = result.labels;
  }
  const ok = result.n_per_group.every(c => c >= ctx.minNGroup);
  const reason = ok ? null
    : (nW < ctx.minNWin ? 'LOW_NWIN' : 'LOW_GROUP_N');

  // ---- Family purity (legacy STEP_D16 v3.2) ----
  let fam_purity = NaN;
  let fam_per_cluster = null;
  if (d.family_source && d.family_source !== 'none' && d.samples) {
    const labels = result.labels;
    const purity_vec = [];
    const size_vec = [];
    fam_per_cluster = [];
    for (let k = 0; k < ctx.k; k++) {
      const counts = new Map();
      let total = 0;
      for (let si = 0; si < labels.length; si++) {
        if (labels[si] !== k) continue;
        const f = d.samples[si].family_id;
        if (f == null || f === -1) continue;
        counts.set(f, (counts.get(f) || 0) + 1);
        total++;
      }
      if (total === 0) {
        fam_per_cluster.push({ purity: NaN, total: 0, top_family: null, top_count: 0 });
        continue;
      }
      let top_family = null, top_count = 0;
      for (const [f, c] of counts) {
        if (c > top_count) { top_count = c; top_family = f; }
      }
      const purity = top_count / total;
      purity_vec.push(purity);
      size_vec.push(total);
      fam_per_cluster.push({ purity, total, top_family, top_count });
    }
    if (purity_vec.length > 0) {
      let num = 0, den = 0;
      for (let i = 0; i < purity_vec.length; i++) {
        num += purity_vec[i] * size_vec[i];
        den += size_vec[i];
      }
      if (den > 0) fam_purity = num / den;
    }
  }

  // ---- Halves coherence (legacy STEP_D16 v3.2) ----
  let coherence = NaN;
  const s0 = env._s0, e0 = env._e0;
  if (e0 - s0 + 1 >= 6) {
    const half = Math.floor((e0 - s0 + 1) / 2);
    const rangeA = [s0, s0 + half - 1];
    const rangeB = [e0 - half + 1, e0];
    const nS = d.n_samples;
    const meanA = new Float64Array(nS);
    const meanB = new Float64Array(nS);
    let nA = 0, nB = 0;
    for (let w = rangeA[0]; w <= rangeA[1]; w++) {
      const { pc1, sign } = ctx.getPC(w);
      for (let si = 0; si < nS; si++) meanA[si] += pc1[si] * sign;
      nA++;
    }
    for (let w = rangeB[0]; w <= rangeB[1]; w++) {
      const { pc1, sign } = ctx.getPC(w);
      for (let si = 0; si < nS; si++) meanB[si] += pc1[si] * sign;
      nB++;
    }
    if (nA > 0 && nB > 0) {
      for (let si = 0; si < nS; si++) { meanA[si] /= nA; meanB[si] /= nB; }
      let mAx = 0, mBx = 0, n_finite = 0;
      for (let si = 0; si < nS; si++) {
        if (isFinite(meanA[si]) && isFinite(meanB[si])) {
          mAx += meanA[si]; mBx += meanB[si]; n_finite++;
        }
      }
      if (n_finite >= 5) {
        mAx /= n_finite; mBx /= n_finite;
        let cov = 0, vA = 0, vB = 0;
        for (let si = 0; si < nS; si++) {
          if (!isFinite(meanA[si]) || !isFinite(meanB[si])) continue;
          const da = meanA[si] - mAx, db = meanB[si] - mBx;
          cov += da * db; vA += da * da; vB += db * db;
        }
        if (vA > 1e-12 && vB > 1e-12) coherence = cov / Math.sqrt(vA * vB);
      }
    }
  }
  const incoherent = isFinite(coherence) && coherence < 0.60;

  return {
    labels: result.labels,
    n_per_group: result.n_per_group,
    centers: result.centers || result.cx,
    centers_y: result.cy || null,
    nW,
    n_below_threshold: nW < ctx.minNWin,
    ok,
    reason,
    fam_purity,
    fam_per_cluster,
    coherence,
    incoherent,
    usedK,
    silhouette,
    fixedKLabels,
  };
}

/**
 * Cluster an L2 at a specific K, ignoring `ctx.kMode` (always fixed-K).
 * Used by the L3 panel when `l3KMode === 'k6'` or `'both'`.
 *
 * @param {ClusterContext} ctx
 * @param {number} l2idx
 * @param {number} K
 * @returns {object}  same shape as clusterL2 minus adaptive fields
 */
export function clusterL2AtK(ctx, l2idx, K) {
  const d = ctx.data;
  const env = d.l2_envelopes[l2idx];
  const nW = env._e0 - env._s0 + 1;
  if (nW < 1) return { ok: false, reason: 'NO_WINDOWS' };
  const agg = aggregateL2(ctx, l2idx);
  if (!agg) return { ok: false, reason: 'NO_WINDOWS' };
  let result;
  if (ctx.aggMethod === 'mean_pc12') {
    result = kmeans2D(agg.xs, agg.ys, K);
  } else {
    result = kmeans1D(agg.xs, K);
  }
  const ok = result.n_per_group.every(c => c >= ctx.minNGroup);
  const reason = ok ? null : (nW < ctx.minNWin ? 'LOW_NWIN' : 'LOW_GROUP_N');
  return {
    labels: result.labels,
    n_per_group: result.n_per_group,
    centers: result.centers || result.cx,
    centers_y: result.cy || null,
    nW,
    n_below_threshold: nW < ctx.minNWin,
    ok,
    reason,
    usedK: K,
    fixedKLabels: result.labels,   // at this K, raw labels ARE the fixed-K labels
  };
}

// =====================================================================
// 5. Cache key helpers — for callers that want to memoize cluster results
// =====================================================================

/**
 * Derive a deterministic cache key from a context. Caller compares
 * this against the previous key to know when to invalidate.
 *
 * @param {ClusterContext} ctx
 * @returns {string}
 */
export function clusterCacheKey(ctx) {
  return [
    ctx.k, ctx.kMode, ctx.aggMethod,
    (ctx.data && (ctx.data.flipPC1 ?? 1)),    // included for parity with legacy ensureL2Cache
    ctx.minNGroup, ctx.minNWin, ctx.silThreshold,
    Array.isArray(ctx.kRange) ? ctx.kRange.join(',') : '',
  ].join('|');
}

/**
 * A simple Map-backed cache. Caller picks when to invalidate.
 * Mirrors the legacy `ensureL2Cache` + `state.l2GroupCache` pattern but
 * without the global state.
 *
 * Usage:
 *   const cache = new ClusterCache();
 *   const cl = cache.getOrCompute(ctx, l2idx);
 *   // ... later ...
 *   cache.invalidateIfChanged(ctx);
 */
export class ClusterCache {
  constructor() {
    this._map = new Map();
    this._key = null;
  }
  invalidateIfChanged(ctx) {
    const k = clusterCacheKey(ctx);
    if (this._key !== k) {
      this._map.clear();
      this._key = k;
    }
  }
  getOrCompute(ctx, l2idx) {
    this.invalidateIfChanged(ctx);
    if (this._map.has(l2idx)) return this._map.get(l2idx);
    const r = clusterL2(ctx, l2idx);
    this._map.set(l2idx, r);
    return r;
  }
  clear() { this._map.clear(); this._key = null; }
  get size() { return this._map.size; }
}

// ---------------------------------------------------------------------
// Console-debug exposures (legacy `getL2Cluster` had `window.getL2Cluster`)
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._aggregateL2       = aggregateL2;
  window._sampleSpreadL2    = sampleSpreadL2;
  window._sampleSpreadRange = sampleSpreadRange;
  window._sigmaProfileL2    = sigmaProfileL2;
  window._clusterL2         = clusterL2;
  window._clusterL2AtK      = clusterL2AtK;
  window._clusterCacheKey   = clusterCacheKey;
  window._ClusterCache      = ClusterCache;
}
