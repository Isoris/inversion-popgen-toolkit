// =============================================================================
// atlas_renderers_turn4.js — Turn 4 of chat A
// =============================================================================
//
// Adds three missing canvas renderers to the popstats stack and wires up a
// static/live/split mode chip on the toolbar. Also ships a small response
// adapter that turns popstats_server payloads into the {mb, series[]} shape
// the multi-line renderer consumes.
//
// Modules:
//   drawPopstatsMultiline   — per-group lines / fill / heatmap
//   drawPopstatsBars        — vertical bars (count-y values)
//   drawPopstatsCategoricalStrip — colored cells for categorical maxQ_label
//   adaptPopstatsResponse   — server payload → per-track {mb, series} map
//   adaptHobsResponse       — hobs payload → per-track {mb, series} map
//   modeChip                — UI element for the popstats toolbar
//   popstatsTrackOverrides  — getData() resolvers that respect mode chip
//
// This file is structured as a registration function so the atlas's existing
// renderer dispatch can call us. After loading, the atlas patches its
// renderer dispatch to include 'multiline', 'bars', and 'categorical_strip'.
//
// API (window.popgenRenderers):
//   .drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, trackDef)
//   .drawPopstatsBars(ctx, toX, pad, plotW, plotH, data, trackDef)
//   .drawPopstatsCategoricalStrip(ctx, toX, pad, plotW, plotH, data, trackDef)
//   .adaptPopstatsResponse(payload) → { theta_invgt, fst_hom1_hom2, ... }
//   .adaptHobsResponse(payload)     → { hobs_hexp }
//   .getMultilineMode(trackId)      → 'curves' | 'fill' | 'heatmap'
//   .setMultilineMode(trackId, m)
//   .getPopstatsMode()              → 'static' | 'live' | 'split'
//   .setPopstatsMode(m)
//   .makeModeChip()                 → DOM element for the toolbar
//   .makeMultilineModeChip(trackId) → DOM element for the track header
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenRenderers = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Color palettes
  // ===========================================================================

  // Default per-series palette for multi-line tracks. Q04 PDF uses muted
  // greys for θπ-by-class; we lean similar but slightly more distinguishable.
  // Track-def can override by passing series with explicit colors.
  const DEFAULT_SERIES_PALETTE = [
    '#1f4e79',  // deep blue   — H1/H1 default
    '#f5a524',  // amber       — H1/H2 default
    '#3cc08a',  // teal-green  — H2/H2 default
    '#e0555c',  // coral       — H1/H3 default
    '#b07cf7',  // violet      — H2/H3 default
    '#7ad3db',  // cyan        — H3/H3 default
    '#ff8c6e',  // salmon
    '#9bbf7c',  // olive
    '#d28cb8',  // pink
    '#5a8fb3',  // slate blue
  ];

  // Categorical strip palette (matches atlas's ANC_PALETTE_BASE for
  // consistency with page 7 ancestry rendering).
  const CATEGORICAL_PALETTE = [
    '#1f4e79', '#f5a524', '#3cc08a', '#e0555c', '#b07cf7',
    '#f0d56a', '#7ad3db', '#ff8c6e', '#9bbf7c', '#d28cb8',
    '#5a8fb3', '#ce9b66',
  ];

  // Stable hash → palette index for unknown labels so they get a consistent
  // color across renders.
  function _hashIdx(label, mod) {
    let h = 0;
    const s = String(label);
    for (let i = 0; i < s.length; i++) h = ((h << 5) - h + s.charCodeAt(i)) | 0;
    return Math.abs(h) % mod;
  }

  function colorForSeries(name, idx) {
    if (typeof idx === 'number' && idx >= 0 && idx < DEFAULT_SERIES_PALETTE.length) {
      return DEFAULT_SERIES_PALETTE[idx];
    }
    return DEFAULT_SERIES_PALETTE[_hashIdx(name, DEFAULT_SERIES_PALETTE.length)];
  }

  function colorForCategory(label) {
    return CATEGORICAL_PALETTE[_hashIdx(label, CATEGORICAL_PALETTE.length)];
  }

  // Theme color helper — atlas exposes themeColor(name); fall back to neutral
  // colors when running outside the atlas (e.g. during unit tests).
  function _theme(name) {
    if (typeof themeColor === 'function') {
      try { return themeColor(name); } catch (_) {}
    }
    const fallback = {
      'ink':  '#0e1116',
      'dim':  '#555e69',
      'rule': '#c8cdd2',
      'bg':   '#ffffff',
    };
    return fallback[name] || '#999';
  }

  // ===========================================================================
  // drawPopstatsMultiline
  // ===========================================================================
  // data shape:
  //   { mb: number[],
  //     series: [
  //       { name: string, color?: string, values: (number|null)[],
  //         dashed?: boolean, alpha?: number }, ...
  //     ],
  //     refLine?: number,
  //     yMin?, yMax?,
  //     mode_override?: 'curves'|'fill'|'heatmap'  // bypass per-track localStorage
  //   }
  //
  // Render decisions:
  //   - mode = trackDef.multilineMode || getMultilineMode(trackDef.id) ||
  //            (series.length === 2 ? 'fill' : 'curves')
  //   - mode='fill' valid only for series.length === 2; falls back to 'curves'
  //   - mode='heatmap' valid for series.length >= 1; renders one row per series

  function drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, trackDef) {
    if (!data || !Array.isArray(data.mb) || !Array.isArray(data.series)) return;
    const series = data.series.filter(s => s && Array.isArray(s.values) &&
                                           s.values.length === data.mb.length);
    if (series.length === 0) return;

    const mode = data.mode_override
              || (trackDef && trackDef.multilineMode)
              || getMultilineMode(trackDef && trackDef.id)
              || (series.length === 2 ? 'fill' : 'curves');

    if (mode === 'heatmap') {
      _drawMultilineHeatmap(ctx, toX, pad, plotW, plotH, data, series, trackDef);
      return;
    }
    if (mode === 'fill' && series.length === 2) {
      _drawMultilineFill(ctx, toX, pad, plotW, plotH, data, series, trackDef);
      return;
    }
    _drawMultilineCurves(ctx, toX, pad, plotW, plotH, data, series, trackDef);
  }

  // -- y-range helper: shared across modes -------------------------------------
  function _computeYRange(data, series) {
    let yMin = Infinity, yMax = -Infinity;
    if (isFinite(data.yMin) && isFinite(data.yMax)) {
      return [data.yMin, data.yMax];
    }
    for (const s of series) {
      for (const v of s.values) {
        if (v == null || !isFinite(v)) continue;
        if (v < yMin) yMin = v;
        if (v > yMax) yMax = v;
      }
    }
    if (!isFinite(yMin) || !isFinite(yMax)) return [null, null];
    if (yMin === yMax) { yMin -= 0.5; yMax += 0.5; }
    return [yMin, yMax];
  }

  // -- Curves mode: N overlaid lines with halo+core stroke --------------------
  function _drawMultilineCurves(ctx, toX, pad, plotW, plotH, data, series, trackDef) {
    const [yMin, yMax] = _computeYRange(data, series);
    if (yMin == null) return;
    const toY = (v) => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

    // Reference line (e.g., Fst = 0)
    if (isFinite(data.refLine)) {
      ctx.save();
      ctx.strokeStyle = '#c0504d';
      ctx.lineWidth = 0.6;
      ctx.setLineDash([4, 3]);
      const yr = toY(data.refLine);
      ctx.beginPath();
      ctx.moveTo(pad.l, yr); ctx.lineTo(pad.l + plotW, yr); ctx.stroke();
      ctx.restore();
    }

    // For each series: white halo (1.4px) under colored core (0.5px).
    // Same aesthetic as drawPopstatsLine. Multi-line uses slightly slimmer
    // strokes so overlapping lines don't blur into each other.
    const mb = data.mb;
    for (let s_i = 0; s_i < series.length; s_i++) {
      const s = series[s_i];
      const color = s.color || colorForSeries(s.name, s_i);
      const alpha = (typeof s.alpha === 'number') ? s.alpha : 1.0;

      // Halo
      ctx.save();
      ctx.globalAlpha = 0.65 * alpha;
      ctx.lineWidth = 1.4;
      ctx.strokeStyle = 'rgba(255,255,255,0.85)';
      _strokeSeries(ctx, mb, s.values, toX, toY, false);
      ctx.restore();

      // Core
      ctx.save();
      ctx.globalAlpha = alpha;
      ctx.lineWidth = 0.7;
      ctx.strokeStyle = color;
      if (s.dashed) ctx.setLineDash([4, 3]);
      _strokeSeries(ctx, mb, s.values, toX, toY, false);
      ctx.restore();
    }

    _drawSeriesLegend(ctx, pad, plotW, series);
    _drawYTicks(ctx, pad, plotH, yMin, yMax);
  }

  // -- Fill mode: K=2 with one as foreground line + the other as fill area ----
  function _drawMultilineFill(ctx, toX, pad, plotW, plotH, data, series, trackDef) {
    const [yMin, yMax] = _computeYRange(data, series);
    if (yMin == null) return;
    const toY = (v) => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

    const fgSeries = series[0];
    const bgSeries = series[1];
    const bgColor  = bgSeries.color || colorForSeries(bgSeries.name, 1);
    const fgColor  = fgSeries.color || colorForSeries(fgSeries.name, 0);

    // Background fill
    ctx.save();
    ctx.globalAlpha = 0.30;
    ctx.fillStyle = bgColor;
    _fillSeries(ctx, data.mb, bgSeries.values, toX, toY, plotH, pad);
    ctx.restore();

    // Background line on top of fill
    ctx.save();
    ctx.globalAlpha = 0.85;
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = bgColor;
    _strokeSeries(ctx, data.mb, bgSeries.values, toX, toY, false);
    ctx.restore();

    // Foreground line with halo
    ctx.save();
    ctx.lineWidth = 1.4;
    ctx.strokeStyle = 'rgba(255,255,255,0.85)';
    _strokeSeries(ctx, data.mb, fgSeries.values, toX, toY, false);
    ctx.lineWidth = 0.7;
    ctx.strokeStyle = fgColor;
    _strokeSeries(ctx, data.mb, fgSeries.values, toX, toY, false);
    ctx.restore();

    _drawSeriesLegend(ctx, pad, plotW, series);
    _drawYTicks(ctx, pad, plotH, yMin, yMax);
  }

  // -- Heatmap mode: K rows, value → viridis-ish color ramp -------------------
  function _drawMultilineHeatmap(ctx, toX, pad, plotW, plotH, data, series, trackDef) {
    const [yMin, yMax] = _computeYRange(data, series);
    if (yMin == null) return;

    const N = data.mb.length;
    const rowH = Math.max(2, Math.floor(plotH / series.length));
    const cellW = Math.max(1, plotW / N);

    for (let s_i = 0; s_i < series.length; s_i++) {
      const s = series[s_i];
      const y0 = pad.t + s_i * rowH;
      // Per-row label on the left
      ctx.save();
      ctx.fillStyle = _theme('dim');
      ctx.font = '9px ui-monospace, monospace';
      ctx.textAlign = 'right';
      ctx.fillText(s.name, pad.l - 6, y0 + rowH / 2 + 3);
      ctx.restore();

      for (let i = 0; i < N; i++) {
        const v = s.values[i];
        if (v == null || !isFinite(v)) continue;
        const t = (v - yMin) / (yMax - yMin);
        ctx.fillStyle = _viridis(t);
        const x = toX(data.mb[i]) - cellW / 2;
        ctx.fillRect(x, y0, cellW + 0.5, rowH - 1);
      }
    }
    _drawSeriesLegend(ctx, pad, plotW, series);
  }

  // -- Helpers ----------------------------------------------------------------

  function _strokeSeries(ctx, mb, values, toX, toY, fill) {
    ctx.beginPath();
    let first = true;
    for (let i = 0; i < mb.length; i++) {
      const v = values[i];
      if (v == null || !isFinite(v)) { first = true; continue; }
      const x = toX(mb[i]), y = toY(v);
      if (first) { ctx.moveTo(x, y); first = false; }
      else { ctx.lineTo(x, y); }
    }
    if (fill) ctx.fill();
    else ctx.stroke();
  }

  function _fillSeries(ctx, mb, values, toX, toY, plotH, pad) {
    let started = false;
    let lastX = null;
    ctx.beginPath();
    for (let i = 0; i < mb.length; i++) {
      const v = values[i];
      const x = toX(mb[i]);
      if (v == null || !isFinite(v)) {
        if (started) {
          // Close current segment to the baseline
          ctx.lineTo(lastX, pad.t + plotH);
          ctx.closePath();
          ctx.fill();
          ctx.beginPath();
          started = false;
        }
        continue;
      }
      const y = toY(v);
      if (!started) {
        ctx.moveTo(x, pad.t + plotH);
        ctx.lineTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
      lastX = x;
    }
    if (started) {
      ctx.lineTo(lastX, pad.t + plotH);
      ctx.closePath();
      ctx.fill();
    }
  }

  function _drawSeriesLegend(ctx, pad, plotW, series) {
    if (!series || series.length === 0) return;
    ctx.save();
    ctx.font = '9px ui-monospace, monospace';
    ctx.textBaseline = 'middle';
    let x = pad.l + 6;
    const y = pad.t + 8;
    for (let s_i = 0; s_i < series.length; s_i++) {
      const s = series[s_i];
      const color = s.color || colorForSeries(s.name, s_i);
      // swatch
      ctx.fillStyle = color;
      ctx.fillRect(x, y - 4, 8, 8);
      ctx.strokeStyle = 'rgba(0,0,0,0.2)';
      ctx.lineWidth = 0.5;
      ctx.strokeRect(x + 0.5, y - 3.5, 7, 7);
      // label
      x += 11;
      ctx.fillStyle = _theme('ink');
      ctx.textAlign = 'left';
      ctx.fillText(s.name, x, y);
      x += ctx.measureText(s.name).width + 10;
      // wrap if needed
      if (x > pad.l + plotW - 30) break;
    }
    ctx.restore();
  }

  function _drawYTicks(ctx, pad, plotH, yMin, yMax) {
    ctx.save();
    ctx.fillStyle = _theme('dim');
    ctx.font = '8px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(yMax.toFixed(2), pad.l - 4, pad.t + 8);
    ctx.fillText(yMin.toFixed(2), pad.l - 4, pad.t + plotH);
    ctx.restore();
  }

  // -- Viridis-ish ramp (5-stop linear interp) --------------------------------
  // 0 → dark purple, 0.25 → deep blue, 0.5 → green-teal, 0.75 → yellow-green,
  // 1 → bright yellow. Cheaper than full viridis; visually similar enough for
  // popstats heatmaps.
  const _RAMP = [
    [0.0, [68, 1, 84]],
    [0.25, [59, 82, 139]],
    [0.5, [33, 145, 140]],
    [0.75, [94, 201, 98]],
    [1.0, [253, 231, 37]],
  ];
  function _viridis(t) {
    if (!isFinite(t)) return '#555';
    if (t <= 0) return 'rgb(68,1,84)';
    if (t >= 1) return 'rgb(253,231,37)';
    for (let i = 1; i < _RAMP.length; i++) {
      const [t1, c1] = _RAMP[i];
      const [t0, c0] = _RAMP[i - 1];
      if (t <= t1) {
        const f = (t - t0) / (t1 - t0);
        const r = Math.round(c0[0] + f * (c1[0] - c0[0]));
        const g = Math.round(c0[1] + f * (c1[1] - c0[1]));
        const b = Math.round(c0[2] + f * (c1[2] - c0[2]));
        return 'rgb(' + r + ',' + g + ',' + b + ')';
      }
    }
    return '#555';
  }

  // ===========================================================================
  // drawPopstatsBars — vertical bars per window
  // ===========================================================================
  // data: { mb: number[], values: (number|null)[], color?: string }
  // Used by tracks declared with renderer: 'bars' (e.g., low_cov_count).

  function drawPopstatsBars(ctx, toX, pad, plotW, plotH, data, trackDef) {
    if (!data || !Array.isArray(data.mb) || !Array.isArray(data.values)) return;
    const mb = data.mb;
    const values = data.values;
    if (mb.length !== values.length) return;
    let yMin = 0, yMax = -Infinity;
    for (const v of values) if (isFinite(v) && v > yMax) yMax = v;
    if (!isFinite(yMax) || yMax <= 0) yMax = 1;
    if (isFinite(data.yMax)) yMax = data.yMax;
    if (isFinite(data.yMin)) yMin = data.yMin;
    const toY = (v) => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;
    const cellW = Math.max(1, plotW / mb.length * 0.8);
    const color = (trackDef && trackDef.color) || data.color || '#c0504d';
    ctx.fillStyle = color;
    for (let i = 0; i < mb.length; i++) {
      const v = values[i];
      if (v == null || !isFinite(v) || v <= yMin) continue;
      const x = toX(mb[i]) - cellW / 2;
      const y = toY(v);
      const h = pad.t + plotH - y;
      if (h <= 0) continue;
      ctx.fillRect(x, y, cellW, h);
    }
    _drawYTicks(ctx, pad, plotH, yMin, yMax);
  }

  // ===========================================================================
  // drawPopstatsCategoricalStrip — colored cells per window
  // ===========================================================================
  // data: { mb: number[], labels: string[], colorMap?: { label: '#hex' } }
  // The Q04 "dom. Q" strip above Δ12 is exactly this. Atlas's
  // state.data.ancestry_window.maxQ_label[] is the natural data source.

  function drawPopstatsCategoricalStrip(ctx, toX, pad, plotW, plotH, data, trackDef) {
    if (!data || !Array.isArray(data.mb) || !Array.isArray(data.labels)) return;
    if (data.mb.length !== data.labels.length) return;
    const cellW = Math.max(1, plotW / data.mb.length);
    const colorMap = data.colorMap || {};
    const stripH = Math.min(plotH - 4, 28);
    const y = pad.t + (plotH - stripH) / 2;
    for (let i = 0; i < data.mb.length; i++) {
      const lbl = data.labels[i];
      if (lbl == null) continue;
      const color = colorMap[lbl] || colorForCategory(String(lbl));
      const x = toX(data.mb[i]) - cellW / 2;
      ctx.fillStyle = color;
      ctx.fillRect(x, y, cellW + 0.5, stripH);
    }
    // Frame
    ctx.strokeStyle = _theme('rule');
    ctx.lineWidth = 0.5;
    ctx.strokeRect(pad.l + 0.5, y + 0.5, plotW, stripH);
  }

  // ===========================================================================
  // Response adapters — server payload → per-track {mb, series} maps
  // ===========================================================================
  // The popstats server returns one envelope with all metrics + groups in
  // long-form columns. The atlas wants per-track shaped data. These adapters
  // do the projection.
  //
  // Input: payload from popstatsGroupwise (a windows[] array with column-style
  // numeric fields like theta_pi_HOM1, Fst_HOM1_HOM2, etc.) + the parsed
  // columns metadata.
  //
  // Output: { trackId → {mb, series, refLine?} }

  function adaptPopstatsResponse(payload) {
    if (!payload || !Array.isArray(payload.windows) || payload.windows.length === 0) {
      return {};
    }
    const wins = payload.windows;
    const cols = payload.columns || Object.keys(wins[0]);

    // mb axis: prefer center_mb if computed; else compute from start+end.
    const mb = wins.map(w => {
      if (typeof w.center_mb === 'number') return w.center_mb;
      if (typeof w.start === 'number' && typeof w.end === 'number') {
        return (w.start + w.end) / 2 / 1e6;
      }
      return null;
    });

    const out = {};

    // Per-group θπ → multiline track 'theta_invgt'
    const thetaCols = cols.filter(c => c.startsWith('theta_pi_'));
    if (thetaCols.length > 0) {
      out.theta_invgt = {
        mb,
        series: thetaCols.map((c, i) => ({
          name: c.slice('theta_pi_'.length),
          color: colorForSeries(c.slice('theta_pi_'.length), i),
          values: wins.map(w => _numOrNull(w[c])),
        })),
      };
    }

    // Per-pair Fst → if exactly one pair, single 'fst_hom1_hom2' style line.
    // If multiple pairs, a multiline track 'fst_pairs' with one series per pair.
    const fstCols = cols.filter(c => c.startsWith('Fst_'));
    if (fstCols.length === 1) {
      // Canonical placeholder track id 'fst_hom1_hom2' uses single line.
      const c = fstCols[0];
      out.fst_hom1_hom2 = {
        mb,
        // single-series multiline: fall through to curves with one line, but
        // also expose .values for the existing line renderer in case the track
        // def is renderer:'line'
        series: [{
          name: c.slice('Fst_'.length),
          color: '#7b3294',
          values: wins.map(w => _numOrNull(w[c])),
        }],
        values: wins.map(w => _numOrNull(w[c])),
        refLine: 0,
      };
    } else if (fstCols.length > 1) {
      out.fst_pairs = {
        mb,
        series: fstCols.map((c, i) => ({
          name: c.slice('Fst_'.length).replace('_', ' vs '),
          color: colorForSeries(c, i),
          values: wins.map(w => _numOrNull(w[c])),
        })),
        refLine: 0,
      };
    }

    // dXY same shape as Fst
    const dxyCols = cols.filter(c => c.startsWith('dXY_'));
    if (dxyCols.length > 0) {
      out.dxy_pairs = {
        mb,
        series: dxyCols.map((c, i) => ({
          name: c.slice('dXY_'.length).replace('_', ' vs '),
          color: colorForSeries(c, i),
          values: wins.map(w => _numOrNull(w[c])),
        })),
      };
    }

    // dA, MI, MInorm — same pattern, made available for the track gallery
    for (const prefix of ['dA_', 'MI_', 'MInorm_']) {
      const matched = cols.filter(c => c.startsWith(prefix));
      if (matched.length === 0) continue;
      const trackId = prefix.toLowerCase().replace(/_$/, '') + '_pairs';
      out[trackId] = {
        mb,
        series: matched.map((c, i) => ({
          name: c.slice(prefix.length).replace('_', ' vs '),
          color: colorForSeries(c, i),
          values: wins.map(w => _numOrNull(w[c])),
        })),
      };
    }

    // Cohort-wide θπ + Tajima's D + S — single-line tracks, useful as
    // "always-on" reference lines in the gallery
    if (cols.includes('theta_pi')) {
      out.theta_pi_cohort = {
        mb, series: [{ name: 'cohort θπ', color: '#7a8398',
                       values: wins.map(w => _numOrNull(w.theta_pi)) }],
        values: wins.map(w => _numOrNull(w.theta_pi)),
      };
    }
    if (cols.includes('tajD')) {
      out.tajima_d = {
        mb, series: [{ name: "Tajima's D", color: '#9bbf7c',
                       values: wins.map(w => _numOrNull(w.tajD)) }],
        values: wins.map(w => _numOrNull(w.tajD)),
        refLine: 0,
      };
    }

    return out;
  }

  function _numOrNull(v) {
    if (v == null) return null;
    const n = Number(v);
    return isFinite(n) ? n : null;
  }

  // -- Hobs / HoverE response adapter ----------------------------------------
  // Hobs payload has groups[gName].scales[scale_label].mean_Hobs / mean_Hexp /
  // HoverE arrays. Adapter extracts the primary scale (default 10kb) and
  // produces a multiline track per metric.
  function adaptHobsResponse(payload, opts) {
    opts = opts || {};
    const scale = opts.scale || '10kb';
    if (!payload || !payload.groups) return {};
    const out = {};
    const groupNames = Object.keys(payload.groups);
    if (groupNames.length === 0) return out;

    // Probe the first group's scale to pull out a center_bp axis (shared
    // across groups — the server emits identical window grids per group)
    const firstScale = payload.groups[groupNames[0]].scales &&
                       payload.groups[groupNames[0]].scales[scale];
    if (!firstScale || !Array.isArray(firstScale.center_bp)) return out;
    const mb = firstScale.center_bp.map(bp => bp / 1e6);

    for (const metric of ['mean_Hobs', 'mean_Hexp', 'HoverE']) {
      const series = groupNames.map((gn, i) => {
        const sc = payload.groups[gn].scales[scale];
        const values = sc && Array.isArray(sc[metric]) ? sc[metric] : [];
        return {
          name: gn,
          color: colorForSeries(gn, i),
          values,
        };
      });
      const trackId = (metric === 'HoverE') ? 'hobs_hexp'
                    : (metric === 'mean_Hobs') ? 'hobs_mean_per_group'
                    :                           'hexp_mean_per_group';
      out[trackId] = {
        mb,
        series,
        refLine: (metric === 'HoverE') ? 1.0 : null,
      };
    }
    return out;
  }

  // -- Ancestry-Q response adapter --------------------------------------------
  // Per-group mean Q-vector → for each Q-component, a multi-line series across
  // groups. So 'ancestry_q1_per_group' has G series (one per group), values
  // are the K1 component mean per window. K times that = many tracks.
  function adaptAncestryResponse(payload) {
    if (!payload || !payload.groups || !Array.isArray(payload.window_mid_bp)) {
      return {};
    }
    const groupNames = Object.keys(payload.groups);
    if (groupNames.length === 0) return {};
    const mb = payload.window_mid_bp.map(bp => bp / 1e6);
    const K = payload.K || (
      payload.groups[groupNames[0]].Q_mean &&
      payload.groups[groupNames[0]].Q_mean[0] ?
      payload.groups[groupNames[0]].Q_mean[0].length : 0);
    if (!K) return {};
    const out = {};
    for (let qi = 0; qi < K; qi++) {
      const trackId = 'ancestry_q' + (qi + 1) + '_per_group';
      out[trackId] = {
        mb,
        series: groupNames.map((gn, gi) => ({
          name: gn,
          color: colorForSeries(gn, gi),
          values: (payload.groups[gn].Q_mean || []).map(row => {
            return (Array.isArray(row) && row[qi] != null) ? Number(row[qi]) : null;
          }),
        })),
      };
    }
    return out;
  }

  // ===========================================================================
  // Mode chip state — static / live / split (popstats page)
  // ===========================================================================

  const POPSTATS_MODE_LS = 'inversion_atlas.popstatsMode';
  const VALID_MODES = ['static', 'live', 'split'];

  function getPopstatsMode() {
    if (typeof localStorage === 'undefined') return 'static';
    try {
      const v = localStorage.getItem(POPSTATS_MODE_LS);
      if (VALID_MODES.indexOf(v) >= 0) return v;
    } catch (_) {}
    return 'static';
  }

  function setPopstatsMode(m) {
    if (VALID_MODES.indexOf(m) < 0) m = 'static';
    if (typeof localStorage !== 'undefined') {
      try { localStorage.setItem(POPSTATS_MODE_LS, m); } catch (_) {}
    }
    // Notify listeners
    if (typeof window !== 'undefined' && window.dispatchEvent) {
      window.dispatchEvent(new CustomEvent('popgen:mode-changed', {
        detail: { mode: m },
      }));
    }
  }

  // -- DOM element factory for the toolbar -----------------------------------
  function makeModeChip(opts) {
    if (typeof document === 'undefined') return null;
    opts = opts || {};
    const wrap = document.createElement('span');
    wrap.className = 'ps-mode-chip-group';
    wrap.style.cssText =
      'display: inline-flex; align-items: center; gap: 4px; ' +
      'margin-left: 12px; padding: 2px 4px; ' +
      'background: var(--panel-2, #f0f1f3); border: 1px solid var(--rule, #c8cdd2); ' +
      'border-radius: 14px; font-family: var(--mono, monospace); font-size: 10px;';
    const label = document.createElement('span');
    label.textContent = 'mode:';
    label.style.cssText = 'color: var(--ink-dim, #555); padding: 0 4px;';
    wrap.appendChild(label);

    const buttons = {};
    for (const m of VALID_MODES) {
      const b = document.createElement('button');
      b.textContent = m;
      b.dataset.mode = m;
      b.style.cssText =
        'padding: 2px 8px; border: 0; background: transparent; ' +
        'color: var(--ink, #0e1116); cursor: pointer; border-radius: 10px; ' +
        'font-family: inherit; font-size: 10px;';
      b.title = ({
        'static': 'Static — read precomputed JSONs (manuscript-reproducible).',
        'live':   'Live — recompute via server using current group composition.',
        'split':  'Split — overlay both static (solid) and live (dashed) for diff.',
      })[m];
      b.addEventListener('click', () => {
        setPopstatsMode(m);
        for (const x of VALID_MODES) {
          buttons[x].style.background = (x === m)
            ? 'var(--accent, #f5a524)'
            : 'transparent';
          buttons[x].style.color = (x === m)
            ? '#0e1116'
            : 'var(--ink, #0e1116)';
          buttons[x].style.fontWeight = (x === m) ? '600' : 'normal';
        }
        if (typeof opts.onChange === 'function') opts.onChange(m);
      });
      buttons[m] = b;
      wrap.appendChild(b);
    }
    // Initialize selection from persisted state
    const cur = getPopstatsMode();
    if (buttons[cur]) {
      buttons[cur].style.background = 'var(--accent, #f5a524)';
      buttons[cur].style.color = '#0e1116';
      buttons[cur].style.fontWeight = '600';
    }
    return wrap;
  }

  // ===========================================================================
  // Per-track multiline mode (curves / fill / heatmap)
  // ===========================================================================

  const MULTILINE_MODE_LS = 'inversion_atlas.multilineMode';

  function _readMultilineModeMap() {
    if (typeof localStorage === 'undefined') return {};
    try {
      const raw = localStorage.getItem(MULTILINE_MODE_LS);
      if (raw) return JSON.parse(raw) || {};
    } catch (_) {}
    return {};
  }

  function _writeMultilineModeMap(map) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(MULTILINE_MODE_LS, JSON.stringify(map)); } catch (_) {}
  }

  function getMultilineMode(trackId) {
    if (!trackId) return null;
    const map = _readMultilineModeMap();
    return map[trackId] || null;
  }

  function setMultilineMode(trackId, mode) {
    if (!trackId) return;
    const valid = ['curves', 'fill', 'heatmap'];
    if (valid.indexOf(mode) < 0) return;
    const map = _readMultilineModeMap();
    map[trackId] = mode;
    _writeMultilineModeMap(map);
    if (typeof window !== 'undefined' && window.dispatchEvent) {
      window.dispatchEvent(new CustomEvent('popgen:multiline-mode-changed', {
        detail: { trackId, mode },
      }));
    }
  }

  function makeMultilineModeChip(trackId) {
    if (typeof document === 'undefined') return null;
    const wrap = document.createElement('span');
    wrap.className = 'ps-multiline-mode-chip';
    wrap.style.cssText =
      'display: inline-flex; gap: 2px; ' +
      'background: rgba(0,0,0,0.05); border-radius: 8px; ' +
      'font-family: var(--mono, monospace); font-size: 9px; padding: 1px;';
    const cur = getMultilineMode(trackId) || 'curves';
    const valid = ['curves', 'fill', 'heatmap'];
    for (const m of valid) {
      const b = document.createElement('button');
      b.textContent = m;
      b.style.cssText =
        'padding: 1px 5px; border: 0; cursor: pointer; ' +
        'background: ' + (m === cur ? 'rgba(245,165,36,0.7)' : 'transparent') + '; ' +
        'border-radius: 6px; color: var(--ink, #0e1116); font-family: inherit; ' +
        'font-size: inherit;';
      b.title = ({
        'curves':  'N overlaid lines',
        'fill':    'foreground line + background fill (K=2 only)',
        'heatmap': 'one row per group, value → color ramp',
      })[m];
      b.addEventListener('click', () => {
        setMultilineMode(trackId, m);
        for (const child of wrap.children) {
          child.style.background = (child.textContent === m)
            ? 'rgba(245,165,36,0.7)'
            : 'transparent';
        }
      });
      wrap.appendChild(b);
    }
    return wrap;
  }

  // ===========================================================================
  // Renderer dispatch patch — call this after the atlas's existing
  // drawPopstatsTracks loads so our new branches are reachable.
  // ===========================================================================
  // The atlas has a hard-coded if/else chain at line 52261 that handles
  // 'ideogram', 'sim_collapse', and (anything-with-getData → drawPopstatsLine).
  // We can't surgically modify that without editing the atlas HTML; the
  // intended integration is to replace that 4-line block with a single
  // dispatch call to popgenRenderers.dispatch(). The patch comment below
  // describes the exact splice.

  function dispatch(t, ctx, toX, pad, plotW, plotH, mbMin, mbMax, bps, wins, state) {
    // This is the unified dispatch. The atlas's existing renderer block
    // becomes:
    //
    //   if (popgenRenderers.dispatch(t, ctx, toX, pad, plotW, plotH,
    //                                 mbMin, mbMax, bps, wins, state)) {
    //     // handled by popgenRenderers
    //   } else if (t.renderer === 'ideogram') {
    //     drawPopstatsIdeogram(...);
    //   } else if (t.renderer === 'sim_collapse') {
    //     drawPopstatsSimCollapse(...);
    //   } else if (typeof t.getData === 'function') {
    //     const data = t.getData(state);
    //     if (data) drawPopstatsLine(ctx, toX, pad, plotW, plotH, data, t);
    //   }
    //
    // Returns true if it handled the track, false if the atlas should fall
    // through to its existing dispatch.
    if (!t || typeof t.renderer !== 'string') return false;
    if (t.renderer === 'multiline') {
      const data = (typeof t.getData === 'function') ? t.getData(state) : null;
      if (!data) return true;  // claimed; nothing to draw
      drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, t);
      return true;
    }
    if (t.renderer === 'bars') {
      const data = (typeof t.getData === 'function') ? t.getData(state) : null;
      if (!data) return true;
      drawPopstatsBars(ctx, toX, pad, plotW, plotH, data, t);
      return true;
    }
    if (t.renderer === 'categorical_strip') {
      const data = (typeof t.getData === 'function') ? t.getData(state) : null;
      if (!data) return true;
      drawPopstatsCategoricalStrip(ctx, toX, pad, plotW, plotH, data, t);
      return true;
    }
    return false;
  }

  // ===========================================================================
  // Public exports
  // ===========================================================================

  const api = {
    // Renderers
    drawPopstatsMultiline,
    drawPopstatsBars,
    drawPopstatsCategoricalStrip,

    // Adapters
    adaptPopstatsResponse,
    adaptHobsResponse,
    adaptAncestryResponse,

    // Mode chips + state
    getPopstatsMode, setPopstatsMode, makeModeChip,
    getMultilineMode, setMultilineMode, makeMultilineModeChip,

    // Dispatch helper
    dispatch,

    // Color helpers (also used by the gallery in turn 6.5)
    colorForSeries, colorForCategory,

    // Constants exposed for tests
    DEFAULT_SERIES_PALETTE, CATEGORICAL_PALETTE, VALID_MODES,
  };

  return api;
}));
