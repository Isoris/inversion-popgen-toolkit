// shared/kmeans.js
// =====================================================================
// Pure K-means primitives. No `state`, no DOM. Drop-in replacements
// for the legacy `kmeans1D`, `kmeans2D`, `silhouette1D`, `adaptiveK1D`.
//
// Source line refs (legacy Inversion_atlas.html, turn 165 close):
//   kmeans1D       line 10092
//   kmeans2D       line 10133
//   silhouette1D   line 10178
//   adaptiveK1D    line 10224
//
// Convention preserved from legacy:
//   - Centers always sorted ascending after fit. Label 0 = lowest center.
//   - Labels are Int8Array (so K up to 127 is fine; atlas uses K=3, K=6).
//   - Empty input returns shape-compatible result with zero counts.
// =====================================================================

/**
 * 1-D K-means with quantile init + 12 Lloyd iterations.
 *
 * @param {ArrayLike<number>} values  input scalars (length n)
 * @param {number} k                  cluster count (k >= 1)
 * @returns {{labels: Int8Array, centers: Float64Array, n_per_group: number[]}}
 */
export function kmeans1D(values, k) {
  const n = values.length;
  const labels = new Int8Array(n);
  if (n === 0) return { labels, centers: new Float64Array(k), n_per_group: new Array(k).fill(0) };
  const sorted = Float64Array.from(values).sort();
  const centers = new Float64Array(k);
  for (let i = 0; i < k; i++) {
    const q = (i + 0.5) / k;
    centers[i] = sorted[Math.min(n - 1, Math.max(0, Math.floor(q * n)))];
  }
  for (let it = 0; it < 12; it++) {
    let changed = false;
    for (let i = 0; i < n; i++) {
      let best = 0, bd = Infinity;
      for (let j = 0; j < k; j++) {
        const dd = Math.abs(values[i] - centers[j]);
        if (dd < bd) { bd = dd; best = j; }
      }
      if (labels[i] !== best) { labels[i] = best; changed = true; }
    }
    const sums = new Float64Array(k);
    const cnts = new Int32Array(k);
    for (let i = 0; i < n; i++) { sums[labels[i]] += values[i]; cnts[labels[i]]++; }
    for (let j = 0; j < k; j++) if (cnts[j] > 0) centers[j] = sums[j] / cnts[j];
    if (!changed) break;
  }
  // Sort centers ascending so label 0 = lowest, label k-1 = highest.
  const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => centers[a] - centers[b]);
  const remap = new Int8Array(k);
  for (let i = 0; i < k; i++) remap[order[i]] = i;
  const newCenters = new Float64Array(k);
  for (let i = 0; i < k; i++) newCenters[i] = centers[order[i]];
  for (let i = 0; i < n; i++) labels[i] = remap[labels[i]];
  const npg = new Array(k).fill(0);
  for (let i = 0; i < n; i++) npg[labels[i]]++;
  return { labels, centers: newCenters, n_per_group: npg };
}

/**
 * 2-D K-means on (PC1, PC2) with PC1-quantile init.
 *
 * @param {ArrayLike<number>} xs  PC1 values
 * @param {ArrayLike<number>} ys  PC2 values (same length)
 * @param {number} k              cluster count
 * @returns {{labels: Int8Array, cx: Float64Array, cy: Float64Array, n_per_group: number[]}}
 */
export function kmeans2D(xs, ys, k) {
  const n = xs.length;
  const labels = new Int8Array(n);
  if (n === 0) return { labels, cx: new Float64Array(k), cy: new Float64Array(k), n_per_group: new Array(k).fill(0) };
  // Init: sort by PC1, take quantile points
  const idxSort = Array.from({ length: n }, (_, i) => i).sort((a, b) => xs[a] - xs[b]);
  const cx = new Float64Array(k), cy = new Float64Array(k);
  for (let i = 0; i < k; i++) {
    const q = (i + 0.5) / k;
    const ii = idxSort[Math.min(n - 1, Math.max(0, Math.floor(q * n)))];
    cx[i] = xs[ii]; cy[i] = ys[ii];
  }
  for (let it = 0; it < 12; it++) {
    let changed = false;
    for (let i = 0; i < n; i++) {
      let best = 0, bd = Infinity;
      for (let j = 0; j < k; j++) {
        const dx = xs[i] - cx[j], dy = ys[i] - cy[j];
        const dd = dx * dx + dy * dy;
        if (dd < bd) { bd = dd; best = j; }
      }
      if (labels[i] !== best) { labels[i] = best; changed = true; }
    }
    const sx = new Float64Array(k), sy = new Float64Array(k);
    const cn = new Int32Array(k);
    for (let i = 0; i < n; i++) { sx[labels[i]] += xs[i]; sy[labels[i]] += ys[i]; cn[labels[i]]++; }
    for (let j = 0; j < k; j++) if (cn[j] > 0) { cx[j] = sx[j] / cn[j]; cy[j] = sy[j] / cn[j]; }
    if (!changed) break;
  }
  // Sort centers by cx so label 0 = lowest PC1
  const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => cx[a] - cx[b]);
  const remap = new Int8Array(k);
  for (let i = 0; i < k; i++) remap[order[i]] = i;
  const ncx = new Float64Array(k), ncy = new Float64Array(k);
  for (let i = 0; i < k; i++) { ncx[i] = cx[order[i]]; ncy[i] = cy[order[i]]; }
  for (let i = 0; i < n; i++) labels[i] = remap[labels[i]];
  const npg = new Array(k).fill(0);
  for (let i = 0; i < n; i++) npg[labels[i]]++;
  return { labels, cx: ncx, cy: ncy, n_per_group: npg };
}

/**
 * Silhouette score for a 1-D K-means clustering. Returns mean silhouette
 * in [-1, +1]. Higher = better-separated clusters.
 *
 * Formal: for each sample i, a(i) = mean abs distance to other samples
 * in same cluster, b(i) = min over other clusters of mean distance.
 * silhouette(i) = (b - a) / max(a, b). Mean over all i.
 *
 * Returns NaN when n < 4 or any cluster has fewer than 2 samples.
 *
 * @param {ArrayLike<number>} values
 * @param {ArrayLike<number>} labels   integer labels in [0, k)
 * @param {number} k
 * @returns {number}
 */
export function silhouette1D(values, labels, k) {
  const n = values.length;
  if (n < 4 || k < 2) return NaN;
  const indicesByK = Array.from({ length: k }, () => []);
  for (let i = 0; i < n; i++) indicesByK[labels[i]].push(i);
  for (let kk = 0; kk < k; kk++) if (indicesByK[kk].length < 2) return NaN;
  let total = 0, ncount = 0;
  for (let i = 0; i < n; i++) {
    const my = labels[i];
    let aSum = 0, aCnt = 0;
    for (const j of indicesByK[my]) {
      if (j === i) continue;
      aSum += Math.abs(values[i] - values[j]); aCnt++;
    }
    const a = aCnt > 0 ? aSum / aCnt : 0;
    let bMin = Infinity;
    for (let other = 0; other < k; other++) {
      if (other === my) continue;
      let bSum = 0, bCnt = 0;
      for (const j of indicesByK[other]) {
        bSum += Math.abs(values[i] - values[j]); bCnt++;
      }
      if (bCnt > 0) {
        const b = bSum / bCnt;
        if (b < bMin) bMin = b;
      }
    }
    if (!isFinite(bMin)) continue;
    const s = (bMin - a) / Math.max(a, bMin, 1e-12);
    total += s; ncount++;
  }
  return ncount > 0 ? total / ncount : NaN;
}

/**
 * Adaptive K selection (matches D16 v3's adaptive_k logic):
 * for each k in [k_min..k_max] run kmeans1D, score by silhouette,
 * pick the k with highest silhouette ABOVE silThreshold. If none clear
 * the threshold, fall back to the smallest k tested.
 *
 * Returns null if values.length < kMin * minNGroup (insufficient data).
 *
 * @param {ArrayLike<number>} values
 * @param {number} kMin
 * @param {number} kMax
 * @param {number} silThreshold   lower bound for "good" silhouette
 * @param {number} minNGroup      minimum members in each cluster
 * @returns {{k: number, silhouette: number, labels: Int8Array, centers: Float64Array, n_per_group: number[]} | null}
 */
export function adaptiveK1D(values, kMin, kMax, silThreshold, minNGroup) {
  if (values.length < kMin * minNGroup) return null;
  let bestK = kMin, bestSil = -Infinity, bestResult = null;
  for (let k = kMin; k <= kMax; k++) {
    if (values.length < k * minNGroup) break;
    const r = kmeans1D(values, k);
    if (r.n_per_group.some(c => c < minNGroup)) continue;
    const sil = silhouette1D(values, r.labels, k);
    if (!isFinite(sil)) continue;
    if (sil > bestSil) { bestSil = sil; bestK = k; bestResult = r; }
  }
  if (bestResult == null) {
    bestResult = kmeans1D(values, kMin);
    bestK = kMin;
    bestSil = silhouette1D(values, bestResult.labels, kMin);
  }
  return { k: bestK, silhouette: bestSil, ...bestResult };
}

// ---------------------------------------------------------------------
// Console-debug exposures (preserves legacy `window.kmeans1D`)
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window.kmeans1D     = kmeans1D;
  window.kmeans2D     = kmeans2D;
  window.silhouette1D = silhouette1D;
  window.adaptiveK1D  = adaptiveK1D;
}
