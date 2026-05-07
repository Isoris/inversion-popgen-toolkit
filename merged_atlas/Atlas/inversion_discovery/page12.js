// Atlas/inversion_discovery/page12.js
//
// Same six-panel layout as page1 but driven by θπ.
// Empty-state until R-pipeline ships theta_pi layers.
//
// Extracted LITERALLY from legacy/Inversion_atlas.html.
// Per HANDOFF_BATCH_1: state is now passed as first arg (was a global in legacy).
// All other unresolved references are flagged TODO_MISSING — the merge chat
// will resolve them once all 5 batches land.
//
// Entry points (in extraction order):
  // _refreshThetaPiLayerStatus() — legacy lines 53045-53061
  // _refreshThetaPiPanelVisibility() — legacy lines 53097-53133
  // _drawThCusumHero() — legacy lines 53174-53399
  // _drawThLinesPanel() — legacy lines 53432-53633
  // _drawThSimMatPanel() — legacy lines 53660-53776
  // _drawThZPanel() — legacy lines 53793-53915
  // _drawThAnchorStripPanel() — legacy lines 53939-54022
  // _drawThPcaPanel() — legacy lines 54045-54168
//

import { contextFromState, clusterL2, ClusterCache } from '../shared/per_l2_cluster.js';
import { hetRateColor } from '../shared/het_rate.js';
import { alignLabels, hungarianChainProjection, concordanceMatrix } from '../shared/hungarian.js';
import { buildContingency, computeARI, computeNMI, cramersV } from '../shared/contingency.js';
import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

// ---------------------------------------------------------------------------
// TODO_MISSING — functions referenced but not defined in this module or shared/*.
// (These are called by the extracted bodies below. The merge chat will either
//  pull them from another batch, or hoist them into shared/ if they are reused.)
// ---------------------------------------------------------------------------
// TODO_MISSING(colorFor)
// TODO_MISSING(fillFor)
// TODO_MISSING(has)
// TODO_MISSING(kColor)
// TODO_MISSING(palette)
// TODO_MISSING(q)
// TODO_MISSING(showHide)
// TODO_MISSING(toX)
// TODO_MISSING(toY)
// TODO_MISSING(xAt)
// TODO_MISSING(xToPx)
// TODO_MISSING(yAt)
// TODO_MISSING(yToPx)

// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT — state.<slot> references that are NOT in
// shared/state.js SLOT_REGISTRY. Some may be ad-hoc geometry caches the
// legacy added imperatively (state._simGeom, state._zGeom, etc.); some may
// be slot names we need to register. Merge chat decides per-slot.
// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT(state._simGeom)
// TODO_MISSING_SLOT(state._thSimGeom)

// ---------------------------------------------------------------------------
// Extracted bodies
// ---------------------------------------------------------------------------

// --- _refreshThetaPiLayerStatus() — legacy lines 53045-53061 ---
export function _refreshThetaPiLayerStatus(state) {
  if (typeof document === 'undefined') return;
  const indicators = document.querySelectorAll('[data-th-layer]');
  if (!indicators || indicators.length === 0) return;
  indicators.forEach(node => {
    const layerName = node.dataset && node.dataset.thLayer;
    if (!layerName) return;
    const present = !!(state.layersPresent && state.layersPresent.has(layerName));
    if (present) {
      node.textContent = '🟢 loaded';
      node.style.color = 'var(--good)';
    } else {
      node.textContent = '⚪ not loaded';
      node.style.color = 'var(--ink-dimmer)';
    }
  });
}

// --- _refreshThetaPiPanelVisibility() — legacy lines 53097-53133 ---
export function _refreshThetaPiPanelVisibility(state) {
  if (typeof document === 'undefined') return;
  const has = (name) => !!(state.layersPresent && state.layersPresent.has(name));

  const hasPerWindow  = has('theta_pi_per_window');
  const hasLocalPCA   = has('theta_pi_local_pca');
  const hasEnvelopes  = has('theta_pi_envelopes');
  const hasCusumTheta = has('cusum_theta');
  const hasAny = hasPerWindow || hasLocalPCA || hasEnvelopes || hasCusumTheta;

  // Per-panel visibility. The visible-display value is 'flex' for the
  // ctrl bar (it lays out its children horizontally) and 'block' for
  // the .panel divs (default block layout).
  const showHide = (id, visible, displayWhenVisible) => {
    const el = document.getElementById(id);
    if (!el) return;
    el.style.display = visible ? (displayWhenVisible || 'block') : 'none';
  };

  showHide('thCtrlBar',                     hasAny,         'flex');
  showHide('thCusumHeroPanel',              hasCusumTheta,  'block');
  showHide('thSimPanel',                    hasLocalPCA,    'block');
  showHide('thZPanel',                      hasLocalPCA,    'block');
  showHide('thLinesPanel',                  hasPerWindow,   'block');
  showHide('thAnchorStripPanel',            hasEnvelopes,   'block');
  showHide('thPcaPanel',                    hasLocalPCA,    'block');
  showHide('thTrackedSamplesPanelCompact',  hasPerWindow,   'block');
  showHide('thL3Panel',                     hasLocalPCA,    'block');

  // Empty-state block hides as soon as ANY θπ layer arrives. The
  // empty-state still serves as the "what this page is, what data
  // it needs" documentation when no θπ data has loaded.
  const empty = document.getElementById('thetaPiEmpty');
  if (empty) {
    empty.style.display = hasAny ? 'none' : 'block';
  }
}

// --- _drawThCusumHero() — legacy lines 53174-53399 ---
export function _drawThCusumHero(state) {
  if (typeof document === 'undefined') return;
  const data = state && state.data && state.data.cusum_theta;
  if (!data || !Array.isArray(data.persample) || data.persample.length === 0) {
    return;  // nothing to draw; visibility wiring keeps panel hidden
  }

  // Determine genomic range. Prefer explicit range_bp; fall back to
  // min/max of cp_bp.
  let xMin, xMax;
  if (data.range_bp && Number.isFinite(data.range_bp.start) &&
      Number.isFinite(data.range_bp.end) && data.range_bp.end > data.range_bp.start) {
    xMin = data.range_bp.start;
    xMax = data.range_bp.end;
  } else {
    let lo = Infinity, hi = -Infinity;
    for (const r of data.persample) {
      if (Number.isFinite(r.cp_bp)) {
        if (r.cp_bp < lo) lo = r.cp_bp;
        if (r.cp_bp > hi) hi = r.cp_bp;
      }
    }
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) return;
    // Pad 5% on each side so end ticks aren't flush against the canvas edge.
    const pad = (hi - lo) * 0.05;
    xMin = lo - pad;
    xMax = hi + pad;
  }
  const xToPx = (bp, w) => ((bp - xMin) / (xMax - xMin)) * w;

  // Karyotype color palette — matches atlas convention elsewhere.
  // Ungrouped (no karyotype field) gets a neutral grey.
  const kColor = (k) => {
    if (k === 'HOM_REF') return '#3a7dde';   // blue
    if (k === 'HET')     return '#d97842';   // orange
    if (k === 'HOM_INV') return '#7c4ad9';   // purple
    return '#9aa3ad';                         // grey: no karyotype
  };

  // ── Lane 1: per-carrier CUSUM strip ──
  // One short horizontal track per carrier, sorted by cp_bp ascending.
  // Tick at the changepoint position. Tick height encodes strength
  // (clamped to [0.1, 1.0] of row height).
  const stripCanvas = document.getElementById('thCusumStripCanvas');
  if (stripCanvas && typeof stripCanvas.getContext === 'function') {
    // Sort a copy by cp_bp so we don't mutate input.
    const sorted = data.persample
      .filter(r => Number.isFinite(r.cp_bp))
      .slice()
      .sort((a, b) => a.cp_bp - b.cp_bp);

    // Compute strength normalization (max strength → full tick height).
    let maxStrength = 0;
    for (const r of sorted) {
      const s = Number.isFinite(r.strength) ? Math.abs(r.strength) : 0;
      if (s > maxStrength) maxStrength = s;
    }
    if (maxStrength <= 0) maxStrength = 1;

    // Use clientWidth/clientHeight if non-zero (real DOM), else fall back
    // to 800 × 110 (vm sandbox / pre-layout). Set canvas pixel dims to
    // match logical CSS dims so coords are 1:1.
    const w = (stripCanvas.clientWidth  || stripCanvas.width  || 800);
    const h = (stripCanvas.clientHeight || stripCanvas.height || 110);
    stripCanvas.width = w;
    stripCanvas.height = h;
    const ctx = stripCanvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, w, h);
      // Light vertical grid every ~Mb for readability.
      ctx.strokeStyle = 'rgba(120,120,120,0.15)';
      ctx.lineWidth = 1;
      const mbStart = Math.ceil(xMin / 1e6);
      const mbEnd = Math.floor(xMax / 1e6);
      for (let mb = mbStart; mb <= mbEnd; mb++) {
        const px = xToPx(mb * 1e6, w);
        ctx.beginPath();
        ctx.moveTo(px, 0);
        ctx.lineTo(px, h);
        ctx.stroke();
      }

      // Per-carrier rows. Pack rows into available height: row height =
      // h / n, capped at 4px so ticks stay visible.
      const n = sorted.length;
      const rowH = Math.max(1, Math.min(4, h / Math.max(n, 1)));
      const tickPad = 1;  // 1 px gap between rows
      for (let i = 0; i < n; i++) {
        const r = sorted[i];
        const yTop = i * rowH;
        const yBot = yTop + Math.max(1, rowH - tickPad);
        const px = xToPx(r.cp_bp, w);
        const sNorm = Math.max(0.1, Math.min(1.0, Math.abs(r.strength || 0) / maxStrength));
        const tickH = (yBot - yTop) * sNorm;
        ctx.fillStyle = kColor(r.karyotype);
        // Asymmetry: rising changepoints (sign +1) drawn from yTop down,
        // falling (-1) drawn from yBot up. Default rising if unspecified.
        const yStart = (r.asymmetry === -1) ? (yBot - tickH) : yTop;
        ctx.fillRect(px, yStart, 1, tickH);
      }
    }
  }

  // ── Lane 2: cohort changepoint histogram ──
  // Simple binned histogram of cp_bp across all carriers. Bin count
  // chosen via Freedman–Diaconis-ish heuristic (clamped to [20, 100]).
  const histCanvas = document.getElementById('thCusumHistCanvas');
  if (histCanvas && typeof histCanvas.getContext === 'function') {
    const w = (histCanvas.clientWidth  || histCanvas.width  || 800);
    const h = (histCanvas.clientHeight || histCanvas.height || 70);
    histCanvas.width = w;
    histCanvas.height = h;
    const ctx = histCanvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, w, h);

      const cps = data.persample
        .map(r => r.cp_bp)
        .filter(v => Number.isFinite(v));
      if (cps.length > 0) {
        // Bin count: ~ sqrt(n), clamped to [20, 100]
        const nBins = Math.max(20, Math.min(100, Math.ceil(Math.sqrt(cps.length) * 4)));
        const binW = (xMax - xMin) / nBins;
        const counts = new Array(nBins).fill(0);
        for (const v of cps) {
          let bin = Math.floor((v - xMin) / binW);
          if (bin < 0) bin = 0;
          if (bin >= nBins) bin = nBins - 1;
          counts[bin]++;
        }
        let maxCount = 0;
        for (const c of counts) if (c > maxCount) maxCount = c;
        if (maxCount > 0) {
          // Median + IQR shading first (so bars draw on top).
          const sorted = cps.slice().sort((a, b) => a - b);
          const q = (p) => {
            const idx = (sorted.length - 1) * p;
            const lo = Math.floor(idx), hi = Math.ceil(idx);
            return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
          };
          const q25 = q(0.25), q50 = q(0.50), q75 = q(0.75);
          const px25 = xToPx(q25, w);
          const px50 = xToPx(q50, w);
          const px75 = xToPx(q75, w);
          // IQR rectangle
          ctx.fillStyle = 'rgba(245,196,58,0.12)';
          ctx.fillRect(Math.min(px25, px75), 0, Math.abs(px75 - px25), h);
          // Median line
          ctx.strokeStyle = 'rgba(245,196,58,0.65)';
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(px50, 0);
          ctx.lineTo(px50, h);
          ctx.stroke();
          // Bars
          ctx.fillStyle = 'rgba(140,170,210,0.85)';
          for (let i = 0; i < nBins; i++) {
            if (counts[i] === 0) continue;
            const x0 = xToPx(xMin + i * binW, w);
            const x1 = xToPx(xMin + (i + 1) * binW, w);
            const barH = (counts[i] / maxCount) * (h - 2);
            ctx.fillRect(x0, h - barH, Math.max(1, x1 - x0 - 0.5), barH);
          }
        }
      }
    }
  }

  // ── Spread badge: IQR-based descriptive classifier ──
  // Thresholds from chat d0ac5b6e manuscript framing — display only,
  // not statistical claims.
  const spreadBadge = document.getElementById('thCusumSpreadBadge');
  if (spreadBadge) {
    const cps = data.persample
      .map(r => r.cp_bp)
      .filter(v => Number.isFinite(v))
      .sort((a, b) => a - b);
    if (cps.length < 4) {
      spreadBadge.textContent = 'spread —';
      spreadBadge.title = 'Not enough carriers to compute IQR (need ≥4).';
    } else {
      const q = (p) => {
        const idx = (cps.length - 1) * p;
        const lo = Math.floor(idx), hi = Math.ceil(idx);
        return cps[lo] + (cps[hi] - cps[lo]) * (idx - lo);
      };
      const iqr = q(0.75) - q(0.25);
      let label, color;
      if (iqr < 100e3) {
        label = 'tight'; color = '#3cc08a';
      } else if (iqr < 500e3) {
        label = 'intermediate'; color = '#f5a524';
      } else {
        label = 'ragged'; color = '#e07a7a';
      }
      const iqrKb = (iqr / 1e3).toFixed(0);
      spreadBadge.textContent = 'spread ' + label + ' · IQR ' + iqrKb + ' kb';
      spreadBadge.style.color = color;
      spreadBadge.title = 'Spread class — descriptive, IQR of carrier ' +
        'changepoint positions. tight < 100 kb, intermediate 100–500 kb, ' +
        'ragged > 500 kb. Not a statistical test of distributional shape.';
    }
  }

  // ── Concord badge: cross-stream concordance with GHSL ──
  // Reads cusum_concordance layer if present; otherwise placeholder.
  const concordBadge = document.getElementById('thCusumConcordBadge');
  if (concordBadge) {
    const cc = state && state.data && state.data.cusum_concordance;
    // Expected shape (future-spec, STEP_DC06):
    //   { theta_ghsl_within_50kb: <int>, n_carriers_compared: <int>, ... }
    if (cc && Number.isFinite(cc.theta_ghsl_within_50kb) &&
        Number.isFinite(cc.n_carriers_compared) && cc.n_carriers_compared > 0) {
      const frac = cc.theta_ghsl_within_50kb / cc.n_carriers_compared;
      const pct = (frac * 100).toFixed(0);
      concordBadge.textContent = 'concord ' + cc.theta_ghsl_within_50kb +
        '/' + cc.n_carriers_compared + ' · ' + pct + '%';
      concordBadge.title = 'Cross-stream concord — fraction of carriers ' +
        'whose θπ + GHSL changepoints land within 50 kb of each other.';
    } else {
      concordBadge.textContent = 'concord —';
      concordBadge.title = 'Cross-stream concord requires the ' +
        'cusum_concordance layer (STEP_DC06_cusum_concordance.R) — not loaded yet.';
    }
  }
}

// --- _drawThLinesPanel() — legacy lines 53432-53633 ---
export function _drawThLinesPanel(state) {
  if (typeof document === 'undefined') return;
  const data = state && state.data && state.data.theta_pi_per_window;

  // Two valid shapes — try Shape B first, fall back to Shape A.
  let samples = null;        // [{ sample_id, theta_pi: [...] }]
  let windowMb = null;       // [<n_windows>] of Mb positions
  let nWin = 0, nSamp = 0;

  if (data && Array.isArray(data.samples) && Array.isArray(data.windows) &&
      data.samples.length > 0 && data.windows.length > 0) {
    // Shape B
    samples = data.samples;
    windowMb = data.windows.map(w => {
      if (Number.isFinite(w.center_mb)) return w.center_mb;
      if (Number.isFinite(w.start_bp))  return w.start_bp / 1e6;
      return null;
    });
    nWin = windowMb.length;
    nSamp = samples.length;
  } else if (state.data && Array.isArray(state.data.windows) &&
             state.data.windows.length > 0 &&
             Array.isArray(state.data.windows[0].theta) &&
             state.data.windows[0].theta.length > 0) {
    // Shape A — recover samples[] from per-window theta arrays.
    // The atlas's top-level state.data.samples carries sample_ids in
    // matching order with windows[i].theta values.
    const winArr = state.data.windows;
    const sampleIds = (state.data.samples && Array.isArray(state.data.samples))
      ? state.data.samples.map((s, i) => (s && s.sample_id) || ('S' + i))
      : winArr[0].theta.map((_, i) => 'S' + i);
    nSamp = winArr[0].theta.length;
    nWin = winArr.length;
    // Build samples[] shape on the fly. Don't mutate state.
    samples = sampleIds.map((sid, si) => {
      const tp = new Array(nWin);
      for (let w = 0; w < nWin; w++) {
        const row = winArr[w].theta;
        tp[w] = (row && Number.isFinite(row[si])) ? row[si] : NaN;
      }
      return { sample_id: sid, theta_pi: tp };
    });
    windowMb = winArr.map(w => {
      if (Number.isFinite(w.center_mb)) return w.center_mb;
      if (Number.isFinite(w.start_bp))  return w.start_bp / 1e6;
      return null;
    });
  } else {
    return;  // no data; visibility wiring keeps panel hidden
  }

  if (nWin < 2 || nSamp === 0) return;

  // Find or create the single canvas inside the container.
  const container = document.getElementById('thLinesCanvasContainer');
  if (!container) return;
  let cv = container.querySelector('canvas');
  if (!cv && typeof document.createElement === 'function') {
    cv = document.createElement('canvas');
    cv.id = 'thLinesCanvas';
    cv.style.display = 'block';
    cv.style.width = '100%';
    cv.style.height = '180px';
    if (typeof container.appendChild === 'function') container.appendChild(cv);
  }
  if (!cv || typeof cv.getContext !== 'function') return;

  const w = (cv.clientWidth  || cv.width  || 800);
  const h = (cv.clientHeight || cv.height || 180);
  cv.width = w;
  cv.height = h;
  const ctx = cv.getContext('2d');
  if (!ctx) return;

  // Compute X bounds from windowMb (skip nulls).
  let mbMin = Infinity, mbMax = -Infinity;
  for (const m of windowMb) {
    if (Number.isFinite(m)) {
      if (m < mbMin) mbMin = m;
      if (m > mbMax) mbMax = m;
    }
  }
  if (!Number.isFinite(mbMin) || !Number.isFinite(mbMax) || mbMax <= mbMin) return;

  // Compute Y bounds across all sample θπ values, ignoring NaN.
  let yMin = Infinity, yMax = -Infinity;
  for (let s = 0; s < nSamp; s++) {
    const tp = samples[s].theta_pi;
    if (!Array.isArray(tp)) continue;
    for (let w0 = 0; w0 < nWin; w0++) {
      const v = tp[w0];
      if (Number.isFinite(v)) {
        if (v < yMin) yMin = v;
        if (v > yMax) yMax = v;
      }
    }
  }
  if (!Number.isFinite(yMin) || !Number.isFinite(yMax) || yMax <= yMin) return;
  // Pad Y by 5% on each side
  const yPad = (yMax - yMin) * 0.05;
  yMin -= yPad;
  yMax += yPad;

  // Layout pads (mirrors page 1 lines panel pad: l 44 / r 16 / t 6 / b 8)
  const pad = { l: 44, r: 16, t: 6, b: 8 };
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b;
  if (plotW <= 0 || plotH <= 0) return;

  const xToPx = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;
  const yToPx = (v)  => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

  ctx.clearRect(0, 0, w, h);

  // Light Mb gridlines (matches the hero panel's gridline style)
  ctx.strokeStyle = 'rgba(120,120,120,0.15)';
  ctx.lineWidth = 1;
  const mbStart = Math.ceil(mbMin);
  const mbEnd = Math.floor(mbMax);
  for (let mb = mbStart; mb <= mbEnd; mb++) {
    const px = xToPx(mb);
    ctx.beginPath();
    ctx.moveTo(px, pad.t);
    ctx.lineTo(px, pad.t + plotH);
    ctx.stroke();
  }

  // Karyotype coloring source — try cluster_labels_theta first, then
  // active candidate's locked_labels, else null (neutral grey).
  // Returns an array of length nSamp, values in {-1, 0, 1, 2} where
  // -1 = no label (grey), 0 = HOM_REF/REF band, 1 = HET, 2 = HOM_INV.
  // For dosage candidates locked_labels uses 0..K-1 indices that map
  // to band order; we use a simple 3-color palette indexed by label %3.
  const labels = (() => {
    const out = new Array(nSamp).fill(-1);
    const clt = state.data && state.data.cluster_labels_theta;
    if (clt && Array.isArray(clt.labels) && clt.labels.length === nSamp) {
      for (let i = 0; i < nSamp; i++) {
        const v = clt.labels[i];
        if (Number.isInteger(v) && v >= 0 && v <= 2) out[i] = v;
      }
      return out;
    }
    const cand = state.candidate;
    if (cand && cand.locked_labels && cand.locked_labels.length === nSamp) {
      for (let i = 0; i < nSamp; i++) {
        const v = cand.locked_labels[i];
        if (Number.isInteger(v) && v >= 0) out[i] = v % 3;
      }
      return out;
    }
    return out;
  })();

  // Color palette (matches hero ticks). Ungrouped = neutral grey, low alpha.
  const colorFor = (lab) => {
    if (lab === 0) return 'rgba(58,125,222,0.55)';   // HOM_REF blue
    if (lab === 1) return 'rgba(217,120,66,0.55)';   // HET orange
    if (lab === 2) return 'rgba(124,74,217,0.55)';   // HOM_INV purple
    return 'rgba(154,163,173,0.30)';                 // ungrouped grey, faint
  };

  // Draw polylines. Two-pass: ungrouped first (so grouped lines paint on top).
  ctx.lineWidth = 0.8;
  for (let pass = 0; pass < 2; pass++) {
    for (let s = 0; s < nSamp; s++) {
      const lab = labels[s];
      if (pass === 0 && lab !== -1) continue;
      if (pass === 1 && lab === -1) continue;
      const tp = samples[s].theta_pi;
      if (!Array.isArray(tp)) continue;

      ctx.strokeStyle = colorFor(lab);
      ctx.beginPath();
      let started = false;
      for (let wi = 0; wi < nWin; wi++) {
        const mb = windowMb[wi];
        const v = tp[wi];
        if (!Number.isFinite(mb) || !Number.isFinite(v)) {
          // Break the polyline at gaps so NaN segments don't bridge.
          started = false;
          continue;
        }
        const px = xToPx(mb);
        const py = yToPx(v);
        if (!started) {
          ctx.moveTo(px, py);
          started = true;
        } else {
          ctx.lineTo(px, py);
        }
      }
      ctx.stroke();
    }
  }

  // Y-axis label (single tick at midpoint to keep this scaffold light).
  ctx.fillStyle = 'rgba(154,163,173,0.85)';
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'right';
  ctx.fillText('θπ', pad.l - 6, pad.t + 10);
}

// --- _drawThSimMatPanel() — legacy lines 53660-53776 ---
export function _drawThSimMatPanel(state) {
  if (typeof document === 'undefined') return;
  const tp = state && state.data && state.data.theta_pi_local_pca;
  if (!tp) return;
  const sim = tp.sim_mat;
  if (!sim || (!Array.isArray(sim) && !ArrayBuffer.isView(sim))) return;

  // Determine matrix side length
  let N = 0;
  if (Number.isFinite(tp.n_windows_thumb) && tp.n_windows_thumb > 0) {
    N = tp.n_windows_thumb | 0;
  } else {
    const s = Math.round(Math.sqrt(sim.length));
    if (s * s === sim.length) N = s;
  }
  if (N < 2) return;

  const canvas = document.getElementById('thSimCanvas');
  if (!canvas || typeof canvas.getContext !== 'function') return;
  const w = (canvas.clientWidth  || canvas.width  || 800);
  const h = (canvas.clientHeight || canvas.height || 600);
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.clearRect(0, 0, w, h);

  // Centered square geometry (mirrors page 1's drawSim layout, simpler
  // — no legend strip in this slice).
  const padTop = 18, padBottom = 18, padLeft = 36, padRight = 36;
  const availW = Math.max(50, w - padLeft - padRight);
  const availH = Math.max(50, h - padTop - padBottom);
  const side = Math.max(50, Math.min(availW, availH));
  const x0 = padLeft + Math.floor((availW - side) / 2);
  const y0 = padTop  + Math.floor((availH - side) / 2);

  // Compute color range from sim_mat min/max (skip diagonal which we'll
  // paint yellow as a visual orientation aid).
  let mn = Infinity, mx = -Infinity;
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      if (i === j) continue;
      const v = sim[i * N + j];
      if (Number.isFinite(v)) {
        if (v < mn) mn = v;
        if (v > mx) mx = v;
      }
    }
  }
  if (!Number.isFinite(mn) || !Number.isFinite(mx) || mx <= mn) return;
  const rng = mx - mn;

  // Build off-screen ImageData at native resolution (N×N), then blit to
  // canvas with smoothing if upscaling significantly.
  let img;
  if (typeof ctx.createImageData === 'function') {
    img = ctx.createImageData(N, N);
  }
  if (img && img.data) {
    // Use page-1's simColor palette when it's globally available so the two
    // pages are visually consistent.
    const palette = (typeof simColor === 'function') ? simColor : null;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const k = i * N + j;
        let r, g, b;
        if (i === j) {
          r = 232; g = 197; b = 71;  // diagonal yellow (matches page 1)
        } else {
          const v0 = sim[k];
          const v = Number.isFinite(v0) ? Math.max(0, Math.min(1, (v0 - mn) / rng)) : 0;
          if (palette) {
            const c = palette(v);
            r = c[0]; g = c[1]; b = c[2];
          } else {
            // Local fallback ramp: dark-blue → bright-yellow
            r = Math.round(30 + v * 220);
            g = Math.round(50 + v * 200);
            b = Math.round(120 - v * 90);
          }
        }
        img.data[k * 4]     = r;
        img.data[k * 4 + 1] = g;
        img.data[k * 4 + 2] = b;
        img.data[k * 4 + 3] = 255;
      }
    }
    // Blit. Prefer createElement('canvas') + putImageData + drawImage
    // (lets us upscale with smoothing). If createElement isn't available
    // (vm sandbox), fall back to putImageData directly at top-left.
    const off = (typeof document.createElement === 'function')
      ? document.createElement('canvas') : null;
    if (off && typeof off.getContext === 'function') {
      off.width = N; off.height = N;
      const offCtx = off.getContext('2d');
      if (offCtx && typeof offCtx.putImageData === 'function') {
        offCtx.putImageData(img, 0, 0);
        if (side / N > 1.5) ctx.imageSmoothingEnabled = true;
        if (typeof ctx.drawImage === 'function') {
          ctx.drawImage(off, x0, y0, side, side);
        }
      }
    } else if (typeof ctx.putImageData === 'function') {
      ctx.putImageData(img, x0, y0);
    }
  }

  // Frame
  ctx.strokeStyle = 'rgba(140,140,140,0.4)';
  ctx.lineWidth = 1;
  ctx.strokeRect(x0 + 0.5, y0 + 0.5, side, side);

  // Persist geometry for any future click-to-jump handler (matches page-1
  // state._simGeom convention so handlers port cleanly).
  state._thSimGeom = { x0, y0, x1: x0 + side, y1: y0 + side, side };
}

// --- _drawThZPanel() — legacy lines 53793-53915 ---
export function _drawThZPanel(state) {
  if (typeof document === 'undefined') return;
  const tp = state && state.data && state.data.theta_pi_local_pca;
  if (!tp) return;
  const z = tp.z;
  if (!z || (!Array.isArray(z) && !ArrayBuffer.isView(z))) return;
  const nWin = z.length;
  if (nWin < 2) return;

  const canvas = document.getElementById('thZCanvas');
  if (!canvas || typeof canvas.getContext !== 'function') return;
  const w = (canvas.clientWidth  || canvas.width  || 800);
  const h = (canvas.clientHeight || canvas.height || 140);
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.clearRect(0, 0, w, h);

  // Layout pad mirrors page 1's drawZ uncollapsed pad.
  const pad = { l: 44, r: 16, t: 6, b: 16 };
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b;
  if (plotW <= 0 || plotH <= 0) return;

  // X-axis: window index → Mb if window position info is available, else
  // raw window index. Prefer state.data.windows[].center_mb (the precomp
  // grid the dosage scrubber uses) so X stays aligned with #thLinesPanel
  // and the hero. Fall back to window indices (0..nWin-1) when absent.
  const winArr = state && state.data && Array.isArray(state.data.windows)
    ? state.data.windows : null;
  let xAt;
  let xMin, xMax;
  if (winArr && winArr.length === nWin) {
    const mbs = new Array(nWin);
    for (let i = 0; i < nWin; i++) {
      const wn = winArr[i];
      mbs[i] = Number.isFinite(wn.center_mb) ? wn.center_mb
             : Number.isFinite(wn.start_bp)  ? wn.start_bp / 1e6
             : null;
    }
    let mbMin = Infinity, mbMax = -Infinity;
    for (const m of mbs) {
      if (Number.isFinite(m)) {
        if (m < mbMin) mbMin = m;
        if (m > mbMax) mbMax = m;
      }
    }
    if (Number.isFinite(mbMin) && Number.isFinite(mbMax) && mbMax > mbMin) {
      xMin = mbMin; xMax = mbMax;
      xAt = (i) => Number.isFinite(mbs[i])
        ? pad.l + ((mbs[i] - xMin) / (xMax - xMin)) * plotW
        : null;
    }
  }
  if (!xAt) {
    xMin = 0; xMax = nWin - 1;
    xAt = (i) => pad.l + (i / (nWin - 1)) * plotW;
  }

  // Y-axis: from 0 (|Z| is always ≥ 0) to max(|z|) padded 10%.
  let yMax = 0;
  for (let i = 0; i < nWin; i++) {
    const v = Math.abs(z[i]);
    if (Number.isFinite(v) && v > yMax) yMax = v;
  }
  if (yMax <= 0) return;
  yMax *= 1.10;

  const yAt = (v) => pad.t + plotH - (v / yMax) * plotH;

  // Mb gridlines (matches hero + lines panels)
  ctx.strokeStyle = 'rgba(120,120,120,0.15)';
  ctx.lineWidth = 1;
  if (Number.isFinite(xMin) && Number.isFinite(xMax)) {
    const mbStart = Math.ceil(xMin);
    const mbEnd = Math.floor(xMax);
    for (let mb = mbStart; mb <= mbEnd; mb++) {
      const px = pad.l + ((mb - xMin) / (xMax - xMin)) * plotW;
      ctx.beginPath();
      ctx.moveTo(px, pad.t);
      ctx.lineTo(px, pad.t + plotH);
      ctx.stroke();
    }
  }

  // Threshold-suggesting horizontal line at |Z|=2 if it falls in range —
  // matches the conventional "noteworthy" threshold without making it a
  // hard claim. Dimmed.
  if (yMax > 2) {
    ctx.strokeStyle = 'rgba(120,120,120,0.25)';
    ctx.beginPath();
    const py = yAt(2);
    ctx.moveTo(pad.l, py);
    ctx.lineTo(pad.l + plotW, py);
    ctx.stroke();
  }

  // Waveform — single polyline, breaks at NaN.
  ctx.strokeStyle = 'rgba(140,170,210,0.85)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  let started = false;
  for (let i = 0; i < nWin; i++) {
    const v = Math.abs(z[i]);
    const px = xAt(i);
    if (!Number.isFinite(v) || px === null) {
      started = false;
      continue;
    }
    const py = yAt(v);
    if (!started) { ctx.moveTo(px, py); started = true; }
    else          { ctx.lineTo(px, py); }
  }
  ctx.stroke();

  // Y-axis label — small "θπ |Z|" at top-left
  ctx.fillStyle = 'rgba(154,163,173,0.85)';
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'right';
  ctx.fillText('|Z|', pad.l - 6, pad.t + 10);
}

// --- _drawThAnchorStripPanel() — legacy lines 53939-54022 ---
export function _drawThAnchorStripPanel(state) {
  if (typeof document === 'undefined') return;
  const env = state && state.data && state.data.theta_pi_envelopes;
  if (!env) return;
  const l1 = Array.isArray(env.l1) ? env.l1 : [];
  const l2 = Array.isArray(env.l2) ? env.l2 : [];
  if (l1.length === 0 && l2.length === 0) return;

  const canvas = document.getElementById('thAnchorStripCanvas');
  if (!canvas || typeof canvas.getContext !== 'function') return;
  const w = (canvas.clientWidth  || canvas.width  || 800);
  const h = (canvas.clientHeight || canvas.height || 24);
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.clearRect(0, 0, w, h);

  // X-axis range: prefer state.data.windows Mb range so we stay aligned
  // with #thLinesPanel and the |Z| panel; fall back to envelope extents.
  let mbMin = Infinity, mbMax = -Infinity;
  const winArr = state.data && Array.isArray(state.data.windows)
    ? state.data.windows : null;
  if (winArr && winArr.length > 0) {
    for (const wn of winArr) {
      const m = Number.isFinite(wn.center_mb) ? wn.center_mb
              : Number.isFinite(wn.start_bp)  ? wn.start_bp / 1e6
              : null;
      if (Number.isFinite(m)) {
        if (m < mbMin) mbMin = m;
        if (m > mbMax) mbMax = m;
      }
    }
  }
  if (!Number.isFinite(mbMin) || !Number.isFinite(mbMax) || mbMax <= mbMin) {
    // Fall back to envelope extents
    for (const arr of [l1, l2]) {
      for (const e of arr) {
        if (Number.isFinite(e.start_bp)) {
          const m = e.start_bp / 1e6;
          if (m < mbMin) mbMin = m;
        }
        if (Number.isFinite(e.end_bp)) {
          const m = e.end_bp / 1e6;
          if (m > mbMax) mbMax = m;
        }
      }
    }
  }
  if (!Number.isFinite(mbMin) || !Number.isFinite(mbMax) || mbMax <= mbMin) return;

  const padL = 44, padR = 16;
  const plotW = Math.max(10, w - padL - padR);
  const xToPx = (mb) => padL + ((mb - mbMin) / (mbMax - mbMin)) * plotW;

  // Top half = L1 (blue), bottom half = L2 (green). Matches the panel
  // description and page 1's #zPanel L1/L2 zone bar convention.
  const halfH = Math.max(2, Math.floor(h / 2) - 1);
  const yL1 = 0;
  const yL2 = h - halfH;

  ctx.fillStyle = 'rgba(48,116,200,0.55)';
  for (const e of l1) {
    if (!Number.isFinite(e.start_bp) || !Number.isFinite(e.end_bp)) continue;
    const x0 = xToPx(e.start_bp / 1e6);
    const x1 = xToPx(e.end_bp / 1e6);
    ctx.fillRect(x0, yL1, Math.max(1, x1 - x0), halfH);
  }
  ctx.fillStyle = 'rgba(82,168,124,0.55)';
  for (const e of l2) {
    if (!Number.isFinite(e.start_bp) || !Number.isFinite(e.end_bp)) continue;
    const x0 = xToPx(e.start_bp / 1e6);
    const x1 = xToPx(e.end_bp / 1e6);
    ctx.fillRect(x0, yL2, Math.max(1, x1 - x0), halfH);
  }

  // Side labels (L1 / L2) — small, dimmed
  ctx.fillStyle = 'rgba(154,163,173,0.85)';
  ctx.font = '9px ui-monospace, monospace';
  ctx.textAlign = 'right';
  ctx.fillText('L1', padL - 4, yL1 + halfH * 0.7);
  ctx.fillText('L2', padL - 4, yL2 + halfH * 0.7);
}

// --- _drawThPcaPanel() — legacy lines 54045-54168 ---
export function _drawThPcaPanel(state) {
  if (typeof document === 'undefined') return;
  const tp = state && state.data && state.data.theta_pi_local_pca;
  if (!tp) return;

  // Resolve pc1 + pc2 vectors. Prefer per-window slice when state.cur is
  // valid and per-window arrays are present.
  let pc1 = null, pc2 = null, nSamp = 0;
  const cur = state.cur;
  if (Number.isFinite(cur) && cur >= 0 &&
      Array.isArray(tp.pc1_by_window) && Array.isArray(tp.pc2_by_window) &&
      Array.isArray(tp.pc1_by_window[cur]) && Array.isArray(tp.pc2_by_window[cur])) {
    pc1 = tp.pc1_by_window[cur];
    pc2 = tp.pc2_by_window[cur];
    nSamp = pc1.length;
  } else if ((Array.isArray(tp.pc1) || ArrayBuffer.isView(tp.pc1)) &&
             (Array.isArray(tp.pc2) || ArrayBuffer.isView(tp.pc2))) {
    pc1 = tp.pc1;
    pc2 = tp.pc2;
    nSamp = pc1.length;
  } else {
    return;  // no usable PC data
  }
  if (nSamp < 2 || pc2.length !== nSamp) return;

  const canvas = document.getElementById('thPcaCanvas');
  if (!canvas || typeof canvas.getContext !== 'function') return;
  const w = (canvas.clientWidth  || canvas.width  || 600);
  const h = (canvas.clientHeight || canvas.height || 400);
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.clearRect(0, 0, w, h);

  // Compute bounds, ignoring NaN
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (let i = 0; i < nSamp; i++) {
    const x = pc1[i], y = pc2[i];
    if (Number.isFinite(x)) {
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
    }
    if (Number.isFinite(y)) {
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    }
  }
  if (!Number.isFinite(xMin) || !Number.isFinite(xMax) || xMax <= xMin) return;
  if (!Number.isFinite(yMin) || !Number.isFinite(yMax) || yMax <= yMin) return;

  // 8% pad on each axis (matches page 1)
  const xPad = (xMax - xMin) * 0.08 || 0.01;
  const yPad = (yMax - yMin) * 0.08 || 0.01;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

  const pad = { l: 38, r: 16, t: 12, b: 24 };
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b;
  if (plotW <= 0 || plotH <= 0) return;

  const toX = (v) => pad.l + ((v - xMin) / (xMax - xMin)) * plotW;
  const toY = (v) => pad.t + (1 - (v - yMin) / (yMax - yMin)) * plotH;

  // Frame
  ctx.strokeStyle = 'rgba(140,140,140,0.4)';
  ctx.lineWidth = 1;
  ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, plotW, plotH);

  // Axis label "PC1" / "PC2"
  ctx.fillStyle = 'rgba(154,163,173,0.85)';
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'center';
  ctx.fillText('θπ PC1', pad.l + plotW / 2, h - 6);
  ctx.save();
  ctx.translate(10, pad.t + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('θπ PC2', 0, 0);
  ctx.restore();

  // Resolve per-sample labels (same cascade as lines panel)
  const labels = (() => {
    const out = new Array(nSamp).fill(-1);
    const clt = state.data && state.data.cluster_labels_theta;
    if (clt && Array.isArray(clt.labels) && clt.labels.length === nSamp) {
      for (let i = 0; i < nSamp; i++) {
        const v = clt.labels[i];
        if (Number.isInteger(v) && v >= 0 && v <= 2) out[i] = v;
      }
      return out;
    }
    const cand = state.candidate;
    if (cand && cand.locked_labels && cand.locked_labels.length === nSamp) {
      for (let i = 0; i < nSamp; i++) {
        const v = cand.locked_labels[i];
        if (Number.isInteger(v) && v >= 0) out[i] = v % 3;
      }
      return out;
    }
    return out;
  })();

  const fillFor = (lab) => {
    if (lab === 0) return 'rgba(58,125,222,0.80)';   // HOM_REF blue
    if (lab === 1) return 'rgba(217,120,66,0.80)';   // HET orange
    if (lab === 2) return 'rgba(124,74,217,0.80)';   // HOM_INV purple
    return 'rgba(154,163,173,0.55)';                 // ungrouped grey
  };

  // Two-pass: ungrouped first, grouped on top
  const r = 2.5;
  for (let pass = 0; pass < 2; pass++) {
    for (let i = 0; i < nSamp; i++) {
      const lab = labels[i];
      if (pass === 0 && lab !== -1) continue;
      if (pass === 1 && lab === -1) continue;
      const x = pc1[i], y = pc2[i];
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      ctx.fillStyle = fillFor(lab);
      ctx.fillRect(toX(x) - r, toY(y) - r, r * 2, r * 2);
    }
  }
}
