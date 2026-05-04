// =============================================================================
// test_atlas_renderers_turn4.js
// =============================================================================
// Tests for the turn-4 renderers, adapters, and mode chip state.
// Mocks CanvasRenderingContext2D to capture draw calls + DOM/localStorage so
// we can run in plain Node.
//
// Run:
//   node test_atlas_renderers_turn4.js
// =============================================================================

'use strict';

// ----- localStorage shim -----------------------------------------------------
const _ls = new Map();
globalThis.localStorage = {
  getItem(k) { return _ls.has(k) ? _ls.get(k) : null; },
  setItem(k, v) { _ls.set(k, String(v)); },
  removeItem(k) { _ls.delete(k); },
  clear() { _ls.clear(); },
};

// ----- Minimal DOM shim for the chip-factory tests --------------------------
function makeFakeElement(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    children: [],
    listeners: new Map(),
    style: { cssText: '' },
    dataset: {},
    attributes: {},
    _parent: null,
    appendChild(child) { child._parent = this; this.children.push(child); return child; },
    addEventListener(type, fn) {
      if (!this.listeners.has(type)) this.listeners.set(type, []);
      this.listeners.get(type).push(fn);
    },
    dispatchEvent(ev) {
      const list = this.listeners.get(ev.type) || [];
      for (const fn of list) fn(ev);
    },
    setAttribute(k, v) { this.attributes[k] = v; },
    set className(v) { this.attributes.class = v; },
    get className() { return this.attributes.class || ''; },
    set textContent(v) { this._text = String(v); },
    get textContent() { return this._text || ''; },
  };
  return el;
}
globalThis.document = { createElement: makeFakeElement };
const _events = [];
globalThis.window = {
  dispatchEvent: (ev) => { _events.push(ev); },
  CustomEvent: function (type, init) {
    return { type, detail: (init && init.detail) || {} };
  },
};
globalThis.CustomEvent = globalThis.window.CustomEvent;

// ----- Mock canvas context — records draw operations -----------------------
function makeMockCtx() {
  const ops = [];
  const ctx = {
    ops,
    fillStyle: '#000',
    strokeStyle: '#000',
    lineWidth: 1,
    globalAlpha: 1,
    font: '10px sans-serif',
    textAlign: 'left',
    textBaseline: 'alphabetic',

    save() { ops.push({ op: 'save' }); },
    restore() { ops.push({ op: 'restore' }); },
    beginPath() { ops.push({ op: 'beginPath' }); },
    closePath() { ops.push({ op: 'closePath' }); },
    moveTo(x, y) { ops.push({ op: 'moveTo', x, y }); },
    lineTo(x, y) { ops.push({ op: 'lineTo', x, y }); },
    stroke() {
      ops.push({ op: 'stroke', strokeStyle: this.strokeStyle,
                 lineWidth: this.lineWidth, alpha: this.globalAlpha });
    },
    fill() {
      ops.push({ op: 'fill', fillStyle: this.fillStyle, alpha: this.globalAlpha });
    },
    fillRect(x, y, w, h) {
      ops.push({ op: 'fillRect', x, y, w, h, fillStyle: this.fillStyle });
    },
    strokeRect(x, y, w, h) {
      ops.push({ op: 'strokeRect', x, y, w, h, strokeStyle: this.strokeStyle });
    },
    clearRect(x, y, w, h) {
      ops.push({ op: 'clearRect', x, y, w, h });
    },
    setLineDash(arr) { ops.push({ op: 'setLineDash', arr }); },
    fillText(text, x, y) {
      ops.push({ op: 'fillText', text, x, y, fillStyle: this.fillStyle,
                 font: this.font });
    },
    measureText(text) { return { width: text.length * 6 }; },
  };
  return ctx;
}

// ----- import the module ----------------------------------------------------
const R = require('./atlas_renderers_turn4.js');

// ----- runner ---------------------------------------------------------------
let _pass = 0, _fail = 0;
function ok(label, cond) {
  if (cond) { console.log('  ok  ' + label); _pass++; }
  else      { console.log('FAIL ' + label); _fail++; }
}
function eq(label, a, b) {
  const A = JSON.stringify(a), B = JSON.stringify(b);
  if (A === B) { console.log('  ok  ' + label); _pass++; }
  else         { console.log('FAIL ' + label + ': ' + A + ' != ' + B); _fail++; }
}
function group(name, fn) { console.log(name); fn(); }

// ----- common stub args -----------------------------------------------------
function commonArgs() {
  const ctx = makeMockCtx();
  const pad = { l: 80, r: 60, t: 14, b: 6 };
  const plotW = 600, plotH = 80;
  const mbMin = 0, mbMax = 100;
  const toX = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;
  return { ctx, pad, plotW, plotH, toX };
}

// =============================================================================
// Tests
// =============================================================================

group('drawPopstatsMultiline — curves mode', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const data = {
    mb: [10, 20, 30, 40, 50, 60, 70, 80],
    series: [
      { name: 'H1/H1', values: [0.1, 0.12, 0.11, null, 0.13, 0.14, 0.12, 0.11] },
      { name: 'H1/H2', values: [0.18, 0.19, 0.21, 0.20, 0.20, 0.19, 0.21, 0.18] },
      { name: 'H2/H2', values: [0.05, 0.06, 0.04, 0.05, 0.07, 0.05, 0.06, 0.05] },
    ],
    yMin: 0, yMax: 0.25,
    refLine: 0.1,
  };
  R.drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, {
    id: 'theta_invgt', renderer: 'multiline',
  });
  // Should have at least one stroke per series + halo + ref line
  const strokes = ctx.ops.filter(o => o.op === 'stroke').length;
  ok('multiple strokes (halo + core × 3 series + refLine)', strokes >= 7);

  // moveTo/lineTo span all three series
  const moveTos = ctx.ops.filter(o => o.op === 'moveTo').length;
  ok('several moveTo for separate series + null breaks', moveTos >= 4);

  // Y ticks rendered
  const yTicks = ctx.ops.filter(o => o.op === 'fillText' &&
                                       /^0\.\d+$/.test(String(o.text)));
  ok('y-tick labels rendered', yTicks.length >= 2);

  // Legend with 3 swatches
  const swatchFills = ctx.ops.filter(o => o.op === 'fillRect' && o.w === 8 && o.h === 8);
  eq('3 legend swatches', swatchFills.length, 3);
});

group('drawPopstatsMultiline — fill mode (K=2)', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const data = {
    mb: [10, 20, 30, 40, 50],
    series: [
      { name: 'foreground', values: [0.5, 0.6, 0.55, 0.7, 0.65] },
      { name: 'background', values: [0.3, 0.35, 0.32, 0.4, 0.38] },
    ],
    yMin: 0, yMax: 1,
    mode_override: 'fill',
  };
  R.drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, {
    id: 'fill_demo', renderer: 'multiline',
  });
  // fill mode renders an actual fill operation
  const fills = ctx.ops.filter(o => o.op === 'fill').length;
  ok('fill operation present', fills >= 1);
});

group('drawPopstatsMultiline — heatmap mode', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const N = 20;
  const data = {
    mb: Array.from({length: N}, (_, i) => 5 + i * 4),
    series: [
      { name: 'g0', values: Array.from({length: N}, (_, i) => Math.sin(i / 3)) },
      { name: 'g1', values: Array.from({length: N}, (_, i) => Math.cos(i / 3)) },
      { name: 'g2', values: Array.from({length: N}, (_, i) => i / N) },
    ],
    mode_override: 'heatmap',
  };
  R.drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, {
    id: 'heat_demo', renderer: 'multiline',
  });
  // 3 series × 20 windows = 60 cells minimum
  const cells = ctx.ops.filter(o => o.op === 'fillRect' &&
                                     /^rgb/.test(String(o.fillStyle)));
  ok('60+ heatmap cells drawn', cells.length >= 60);
  // Per-row labels (3) PLUS legend labels at top (3) = 6 fillText calls with
  // those exact strings. Row labels are right-aligned at left margin; legend
  // labels are at top with swatches. We verify both groups appear.
  const rowAndLegendLabels = ctx.ops.filter(o => o.op === 'fillText' &&
                                         ['g0','g1','g2'].includes(String(o.text)));
  eq('3 row labels + 3 legend labels = 6', rowAndLegendLabels.length, 6);
});

group('drawPopstatsMultiline — auto-fallback fill→curves for K≠2', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  // Set per-track mode = 'fill' but provide K=3 series → should fall back
  R.setMultilineMode('threehoz', 'fill');
  const data = {
    mb: [10, 20, 30],
    series: [
      { name: 'a', values: [1, 2, 3] },
      { name: 'b', values: [2, 3, 4] },
      { name: 'c', values: [3, 4, 5] },
    ],
    yMin: 0, yMax: 5,
  };
  R.drawPopstatsMultiline(ctx, toX, pad, plotW, plotH, data, {
    id: 'threehoz', renderer: 'multiline',
  });
  // No `fill` op (only stroke), since K=3 falls back to curves
  const fills = ctx.ops.filter(o => o.op === 'fill').length;
  eq('no fill op (fell back to curves)', fills, 0);
  const strokes = ctx.ops.filter(o => o.op === 'stroke').length;
  ok('strokes present', strokes >= 6);
});

group('drawPopstatsBars', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const data = {
    mb: [10, 30, 50, 70, 90],
    values: [3, 0, 5, 2, null],
    color: '#c0504d',
  };
  R.drawPopstatsBars(ctx, toX, pad, plotW, plotH, data, { id: 'low_cov_count' });
  // Three non-zero bars expected (skipping 0 and null)
  const bars = ctx.ops.filter(o => o.op === 'fillRect' && o.fillStyle === '#c0504d');
  eq('3 non-zero non-null bars', bars.length, 3);
});

group('drawPopstatsCategoricalStrip', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const data = {
    mb: [10, 20, 30, 40, 50],
    labels: ['K1', 'K1', 'K2', 'K2', 'K3'],
    colorMap: { K1: '#1f4e79', K2: '#f5a524', K3: '#3cc08a' },
  };
  R.drawPopstatsCategoricalStrip(ctx, toX, pad, plotW, plotH, data, {
    id: 'dom_q_strip',
  });
  // 5 cells + 1 frame
  const cells = ctx.ops.filter(o => o.op === 'fillRect');
  eq('5 categorical cells', cells.length, 5);
  // Cell 0 should be K1 color
  eq('cell 0 color = K1', cells[0].fillStyle, '#1f4e79');
});

// =============================================================================
// Adapter tests — pure data transforms
// =============================================================================

group('adaptPopstatsResponse — single Fst pair', () => {
  const payload = {
    columns: ['window_id', 'chrom', 'start', 'end', 'center_mb',
              'theta_pi', 'tajD', 'theta_pi_HOM1', 'theta_pi_HET',
              'theta_pi_HOM2', 'Fst_HOM1_HOM2'],
    windows: [
      { window_id: 'w1', chrom: 'C_gar_LG28', start: 0,        end: 50000,
        center_mb: 0.025, theta_pi: 0.005, tajD: -0.3,
        theta_pi_HOM1: 0.004, theta_pi_HET: 0.008, theta_pi_HOM2: 0.003,
        Fst_HOM1_HOM2: 0.05 },
      { window_id: 'w2', chrom: 'C_gar_LG28', start: 50000,    end: 100000,
        center_mb: 0.075, theta_pi: 0.006, tajD: -0.2,
        theta_pi_HOM1: 0.005, theta_pi_HET: 0.009, theta_pi_HOM2: 0.004,
        Fst_HOM1_HOM2: 0.07 },
    ],
  };
  const tracks = R.adaptPopstatsResponse(payload);
  ok('theta_invgt produced',     tracks.theta_invgt != null);
  ok('fst_hom1_hom2 produced',   tracks.fst_hom1_hom2 != null);
  ok('theta_pi_cohort produced', tracks.theta_pi_cohort != null);
  ok('tajima_d produced',        tracks.tajima_d != null);

  eq('theta_invgt has 3 series', tracks.theta_invgt.series.length, 3);
  eq('theta_invgt mb len',       tracks.theta_invgt.mb.length, 2);
  eq('theta_invgt HOM1 values',  tracks.theta_invgt.series[0].values, [0.004, 0.005]);
  eq('Fst single line',          tracks.fst_hom1_hom2.values, [0.05, 0.07]);
  eq('Tajima ref line at 0',     tracks.tajima_d.refLine, 0);
});

group('adaptPopstatsResponse — multiple Fst pairs', () => {
  const payload = {
    columns: ['center_mb', 'Fst_A_B', 'Fst_A_C', 'Fst_B_C'],
    windows: [
      { center_mb: 1.0, Fst_A_B: 0.10, Fst_A_C: 0.20, Fst_B_C: 0.05 },
      { center_mb: 2.0, Fst_A_B: 0.12, Fst_A_C: 0.18, Fst_B_C: 0.06 },
    ],
  };
  const tracks = R.adaptPopstatsResponse(payload);
  ok('fst_pairs (multiple) produced', tracks.fst_pairs != null);
  ok('no fst_hom1_hom2 single-line track', tracks.fst_hom1_hom2 == null);
  eq('3 series', tracks.fst_pairs.series.length, 3);
  eq('first series name', tracks.fst_pairs.series[0].name, 'A vs B');
});

group('adaptHobsResponse', () => {
  const payload = {
    groups: {
      HOM1: { scales: { '10kb': {
        center_bp: [25000, 75000, 125000],
        mean_Hobs: [0.05, 0.06, 0.04],
        mean_Hexp: [0.10, 0.11, 0.09],
        HoverE:    [0.5,  0.55, 0.44],
      }}},
      HOM2: { scales: { '10kb': {
        center_bp: [25000, 75000, 125000],
        mean_Hobs: [0.04, 0.05, 0.03],
        mean_Hexp: [0.10, 0.10, 0.09],
        HoverE:    [0.4,  0.5,  0.33],
      }}},
    },
  };
  const tracks = R.adaptHobsResponse(payload);
  ok('hobs_hexp track', tracks.hobs_hexp != null);
  ok('hobs_mean_per_group', tracks.hobs_mean_per_group != null);
  ok('hexp_mean_per_group', tracks.hexp_mean_per_group != null);
  eq('hobs_hexp mb',     tracks.hobs_hexp.mb, [0.025, 0.075, 0.125]);
  eq('hobs_hexp K series', tracks.hobs_hexp.series.length, 2);
  eq('refLine at 1.0',   tracks.hobs_hexp.refLine, 1.0);
});

group('adaptAncestryResponse', () => {
  const payload = {
    K: 3,
    window_mid_bp: [50_000, 150_000, 250_000],
    groups: {
      HOM1: { Q_mean: [
        [0.7, 0.2, 0.1],
        [0.65, 0.25, 0.10],
        [0.60, 0.30, 0.10],
      ]},
      HOM2: { Q_mean: [
        [0.3, 0.5, 0.2],
        [0.35, 0.45, 0.20],
        [0.40, 0.40, 0.20],
      ]},
    },
  };
  const tracks = R.adaptAncestryResponse(payload);
  ok('Q1 track', tracks.ancestry_q1_per_group != null);
  ok('Q2 track', tracks.ancestry_q2_per_group != null);
  ok('Q3 track', tracks.ancestry_q3_per_group != null);
  eq('Q1 series count', tracks.ancestry_q1_per_group.series.length, 2);
  eq('Q1 HOM1 values',  tracks.ancestry_q1_per_group.series[0].values, [0.7, 0.65, 0.6]);
});

// =============================================================================
// Mode chip tests
// =============================================================================

group('popstats mode chip state', () => {
  R.setPopstatsMode('live');
  eq('mode persisted',   R.getPopstatsMode(), 'live');
  R.setPopstatsMode('static');
  eq('mode revert',      R.getPopstatsMode(), 'static');
  R.setPopstatsMode('badvalue');
  eq('invalid → static', R.getPopstatsMode(), 'static');
});

group('popstats mode chip event', () => {
  _events.length = 0;
  R.setPopstatsMode('split');
  const modeEvents = _events.filter(e => e.type === 'popgen:mode-changed');
  ok('event dispatched',         modeEvents.length === 1);
  eq('event detail',             modeEvents[0].detail.mode, 'split');
});

group('multiline mode chip state', () => {
  R.setMultilineMode('theta_invgt', 'curves');
  R.setMultilineMode('hobs_hexp', 'heatmap');
  eq('theta_invgt = curves',  R.getMultilineMode('theta_invgt'), 'curves');
  eq('hobs_hexp   = heatmap', R.getMultilineMode('hobs_hexp'), 'heatmap');
  // Invalid value rejected
  R.setMultilineMode('theta_invgt', 'invalid');
  eq('invalid kept curves',   R.getMultilineMode('theta_invgt'), 'curves');
});

group('makeModeChip DOM construction', () => {
  const chip = R.makeModeChip({ onChange: () => {} });
  ok('chip element',          chip != null);
  // Three buttons + one label
  eq('children = label + 3 buttons', chip.children.length, 4);
  eq('button text 0 = mode:',         chip.children[0].textContent, 'mode:');
  eq('button text 1 = static',        chip.children[1].textContent, 'static');
  eq('button text 2 = live',          chip.children[2].textContent, 'live');
  eq('button text 3 = split',         chip.children[3].textContent, 'split');
});

group('makeMultilineModeChip DOM construction', () => {
  const chip = R.makeMultilineModeChip('theta_invgt');
  ok('chip element',     chip != null);
  eq('3 buttons',        chip.children.length, 3);
  eq('button 1 = curves',  chip.children[0].textContent, 'curves');
  eq('button 2 = fill',    chip.children[1].textContent, 'fill');
  eq('button 3 = heatmap', chip.children[2].textContent, 'heatmap');
});

// =============================================================================
// Dispatch tests
// =============================================================================

group('dispatch — multiline track is handled', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const tDef = {
    id: 'theta_invgt',
    renderer: 'multiline',
    getData: () => ({
      mb: [1, 2, 3],
      series: [
        { name: 'g0', values: [0.1, 0.2, 0.15] },
        { name: 'g1', values: [0.2, 0.3, 0.25] },
      ],
      yMin: 0, yMax: 0.5,
    }),
  };
  const handled = R.dispatch(tDef, ctx, toX, pad, plotW, plotH, 0, 100, [], [], {});
  ok('returned true (handled)', handled === true);
  ok('strokes were drawn',      ctx.ops.some(o => o.op === 'stroke'));
});

group('dispatch — bars track is handled', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const tDef = {
    id: 'low_cov_count',
    renderer: 'bars',
    getData: () => ({ mb: [10, 20], values: [3, 5] }),
  };
  const handled = R.dispatch(tDef, ctx, toX, pad, plotW, plotH, 0, 100, [], [], {});
  ok('returned true (handled)', handled === true);
  ok('fillRect drawn',          ctx.ops.some(o => o.op === 'fillRect'));
});

group('dispatch — non-claimed renderers fall through', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  // 'line' is claimed by atlas existing dispatch, not our patch
  const tDef = { id: 'whatever', renderer: 'line', getData: () => ({ mb: [], values: [] }) };
  const handled = R.dispatch(tDef, ctx, toX, pad, plotW, plotH, 0, 100, [], [], {});
  ok('returned false (fall through)', handled === false);
  ok('no draw operations on ctx',    ctx.ops.length === 0);
});

group('dispatch — multiline with null getData claims slot but draws nothing', () => {
  const { ctx, pad, plotW, plotH, toX } = commonArgs();
  const tDef = { id: 'no_data', renderer: 'multiline', getData: () => null };
  const handled = R.dispatch(tDef, ctx, toX, pad, plotW, plotH, 0, 100, [], [], {});
  ok('claimed (true)', handled === true);
  ok('no draw ops',    ctx.ops.length === 0);
});

// ----- summary --------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);
