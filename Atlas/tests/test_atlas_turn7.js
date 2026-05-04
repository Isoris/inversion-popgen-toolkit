// =============================================================================
// test_atlas_turn7.js
// =============================================================================
// Tests for Q09b shelf-LD test + page-7 stratified Δ12 resolver.
//
// Run: node test_atlas_turn7.js
// =============================================================================

'use strict';

// ----- localStorage --------------------------------------------------------
const _ls = new Map();
globalThis.localStorage = {
  getItem(k) { return _ls.has(k) ? _ls.get(k) : null; },
  setItem(k, v) { _ls.set(k, String(v)); },
  removeItem(k) { _ls.delete(k); },
  clear() { _ls.clear(); },
};

// ----- DOM shim (slim — only need createElement for makeShelfLDPanel) ------
function makeEl(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    nodeType: 1, children: [], parentNode: null,
    style: null, dataset: {}, attributes: {}, listeners: new Map(),
    width: 360, height: 360,
    appendChild(c) { c.parentNode = this; this.children.push(c); return c; },
    removeChild(c) {
      const i = this.children.indexOf(c);
      if (i >= 0) this.children.splice(i, 1);
      c.parentNode = null; return c;
    },
    addEventListener(type, fn) {
      if (!this.listeners.has(type)) this.listeners.set(type, []);
      this.listeners.get(type).push(fn);
    },
    dispatchEvent(ev) {
      const list = this.listeners.get(ev.type) || [];
      for (const fn of list) fn(ev);
    },
    setAttribute(k, v) {
      this.attributes[k] = String(v);
      if (k.startsWith('data-')) this.dataset[k.slice(5)] = String(v);
    },
    getAttribute(k) { return this.attributes[k] !== undefined ? this.attributes[k] : null; },
    set className(v) { this.attributes.class = v; },
    get className()  { return this.attributes.class || ''; },
    set textContent(v) { this._text = String(v); this.children.length = 0; },
    get textContent()  { return this._text || ''; },
    set id(v) { this.attributes.id = v; },
    get id()  { return this.attributes.id || ''; },
    set title(v) { this.attributes.title = v; },
    getContext() {
      // Stub canvas context that records ops
      return {
        ops: [],
        clearRect() {}, fillRect() {}, strokeRect() {},
        save() {}, restore() {},
        beginPath() {}, moveTo() {}, lineTo() {}, stroke() {}, fill() {},
        translate() {}, rotate() {},
        fillText() {},
        set fillStyle(v) {}, get fillStyle() { return ''; },
        set strokeStyle(v) {}, get strokeStyle() { return ''; },
        set lineWidth(v) {}, get lineWidth() { return 1; },
        set font(v) {}, get font() { return '10px monospace'; },
        set textAlign(v) {}, get textAlign() { return 'left'; },
      };
    },
  };
  const props = {};
  el.style = new Proxy({}, {
    get(_, p) { return props[String(p)] !== undefined ? props[String(p)] : ''; },
    set(_, p, v) { props[String(p)] = String(v); return true; },
  });
  return el;
}
globalThis.document = {
  body: makeEl('body'), head: makeEl('head'), readyState: 'complete',
  createElement: makeEl,
  getElementById(id) {
    function walk(n) {
      if (n.attributes && n.attributes.id === id) return n;
      for (const c of n.children || []) { const f = walk(c); if (f) return f; }
      return null;
    }
    return walk(document.head) || walk(document.body);
  },
  addEventListener() {}, removeEventListener() {},
};
globalThis.window = {
  innerWidth: 1280, innerHeight: 800,
  addEventListener() {}, removeEventListener() {},
  dispatchEvent() {},
  popgenRenderers: null,
  state: null,
};

// ----- runner -------------------------------------------------------------
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
function approx(label, a, b, tol) {
  tol = tol || 0.01;
  if (Math.abs(a - b) <= tol) { console.log('  ok  ' + label); _pass++; }
  else         { console.log('FAIL ' + label + ': ' + a + ' ≉ ' + b); _fail++; }
}
function group(name, fn) { console.log(name); return fn(); }

// ----- import the module under test --------------------------------------
const mod = require('./atlas_turn7.js');

// =============================================================================
// Tests
// =============================================================================

(async () => {

group('pearsonCorrelationMatrix — perfect correlation', () => {
  // Two columns: y = 2x + 5 across 5 rows
  const M = [
    [1, 7],
    [2, 9],
    [3, 11],
    [4, 13],
    [5, 15],
  ];
  const c = mod.pearsonCorrelationMatrix(M);
  approx('r(0,1) ≈ 1', c[0][1], 1.0);
  approx('r(1,0) ≈ 1', c[1][0], 1.0);
  eq('diagonal = 1', c[0][0], 1.0);
});

group('pearsonCorrelationMatrix — anti-correlation', () => {
  const M = [
    [1, 5],
    [2, 4],
    [3, 3],
    [4, 2],
    [5, 1],
  ];
  const c = mod.pearsonCorrelationMatrix(M);
  approx('r(0,1) ≈ -1', c[0][1], -1.0);
});

group('pearsonCorrelationMatrix — no correlation', () => {
  // x = [1,2,3,4,5,6,7,8], y = [3, 1, 4, 1, 5, 9, 2, 6]  (≈ uncorrelated)
  const M = [
    [1, 3], [2, 1], [3, 4], [4, 1], [5, 5], [6, 9], [7, 2], [8, 6],
  ];
  const c = mod.pearsonCorrelationMatrix(M);
  ok('|r(0,1)| < 0.5', Math.abs(c[0][1]) < 0.5);
});

group('pearsonCorrelationMatrix — pairwise NA handling', () => {
  // Column 0 has NaN at row 1, column 1 has NaN at row 3
  // After dropping pairs with NA: rows 0, 2, 4 remain
  const M = [
    [1, 2],
    [NaN, 5],
    [3, 6],
    [4, NaN],
    [5, 10],
  ];
  const c = mod.pearsonCorrelationMatrix(M);
  // y = 2x exactly on rows 0,2,4 → r = 1
  approx('r(0,1) ≈ 1 from 3 pairs', c[0][1], 1.0);
});

group('pearsonCorrelationMatrix — too few paired obs returns NaN', () => {
  const M = [
    [1, NaN],
    [NaN, 2],
    [3, NaN],
  ];
  const c = mod.pearsonCorrelationMatrix(M);
  ok('NaN when no pairs', isNaN(c[0][1]));
});

group('verdictFromCorrMatrix — uniform high → SINGLE_INVERSION', () => {
  const N = 10;
  const c = [];
  for (let i = 0; i < N; i++) {
    const row = [];
    for (let j = 0; j < N; j++) row.push(i === j ? 1.0 : 0.85);
    c.push(row);
  }
  const v = mod.verdictFromCorrMatrix(c, N);
  eq('verdict = SINGLE_INVERSION', v.verdict, 'SINGLE_INVERSION');
  approx('between ≈ 0.85', v.between_halves_corr, 0.85);
});

group('verdictFromCorrMatrix — block structure → LIKELY_MULTIPLE_ARRANGEMENTS', () => {
  const N = 10;
  const c = [];
  for (let i = 0; i < N; i++) {
    const row = [];
    for (let j = 0; j < N; j++) {
      if (i === j) row.push(1.0);
      else if ((i < N/2) === (j < N/2)) row.push(0.85);   // same block
      else row.push(0.05);                                // cross block
    }
    c.push(row);
  }
  const v = mod.verdictFromCorrMatrix(c, N);
  eq('verdict = LIKELY_MULTIPLE_ARRANGEMENTS', v.verdict, 'LIKELY_MULTIPLE_ARRANGEMENTS');
  approx('within-A ≈ 0.85', v.within_first_half_corr, 0.85);
  approx('between ≈ 0.05', v.between_halves_corr, 0.05);
});

group('verdictFromCorrMatrix — mid value → AMBIGUOUS', () => {
  const N = 6;
  const c = [];
  for (let i = 0; i < N; i++) {
    const row = [];
    for (let j = 0; j < N; j++) row.push(i === j ? 1.0 : 0.5);
    c.push(row);
  }
  const v = mod.verdictFromCorrMatrix(c, N);
  eq('verdict = AMBIGUOUS', v.verdict, 'AMBIGUOUS_OR_PARTIAL_LD');
});

// =============================================================================
// runShelfLDTest — synthetic data
// =============================================================================
//
// Build a chunk that simulates a clean single-inversion shelf:
// - 60 samples: 20 Hom1, 20 Het, 20 Hom2
// - 200 SNPs across [10..18 Mb] with one polarity (Hom2 dosage > Hom1)
// - Per-sample dosages drawn deterministically to match the karyotype

function makeSingleInversionChunk(N_HOM = 20, N_HET = 20, N_HOM2 = 20, N_SNPS = 200) {
  const samples = [];
  for (let i = 0; i < N_HOM;  i++) samples.push('h1_' + i);
  for (let i = 0; i < N_HET;  i++) samples.push('he_' + i);
  for (let i = 0; i < N_HOM2; i++) samples.push('h2_' + i);
  const N = samples.length;

  const markers = [];
  const dosage  = [];
  for (let m = 0; m < N_SNPS; m++) {
    const bp = Math.round(10_000_000 + (m / N_SNPS) * 8_000_000);
    markers.push({ marker_id: 'snp_' + m, pos_bp: bp, missingness: 0 });
    const row = [];
    // Hom1: dosage 0; Het: dosage 1; Hom2: dosage 2 (polarity = +1)
    for (let i = 0; i < N_HOM;  i++) row.push(0);
    for (let i = 0; i < N_HET;  i++) row.push(1);
    for (let i = 0; i < N_HOM2; i++) row.push(2);
    dosage.push(row);
  }
  return { samples, markers, dosage };
}

function makeInvgtMap(chunk, n_hom1, n_het, n_hom2) {
  const m = new Map();
  let i = 0;
  for (let k = 0; k < n_hom1; k++) m.set(chunk.samples[i++], 'Hom1');
  for (let k = 0; k < n_het;  k++) m.set(chunk.samples[i++], 'Het');
  for (let k = 0; k < n_hom2; k++) m.set(chunk.samples[i++], 'Hom2');
  return m;
}

group('runShelfLDTest — clean single inversion', () => {
  const chunk = makeSingleInversionChunk();
  const invgt = makeInvgtMap(chunk, 20, 20, 20);
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    n_bins: 20,
  });
  ok('ok',                r.ok);
  eq('200 shelf SNPs',    r.n_shelf_snps, 200);
  // Every SNP is 100% diagnostic in this construction
  eq('200 strong diag',   r.n_diag_strong, 200);
  // Verdict should be SINGLE_INVERSION
  eq('verdict = SINGLE_INVERSION', r.summary.verdict, 'SINGLE_INVERSION');
  ok('between > 0.7',     r.summary.between_halves_corr > 0.7);
});

// Now build a TWO-arrangement chunk:
// - Same 60 samples, but the LEFT HALF of the shelf and the RIGHT HALF
//   independently partition samples into different "Hom1" / "Hom2" groups.
// - Specifically: in left half, Hom1 = first 20, Hom2 = last 20.
//   In right half, Hom1 = first 30, Hom2 = last 30 (with a different cut).
// - Real fish_calls labels follow LEFT half (so the diagnostic SNPs are
//   defined off the left), and RIGHT half SNPs will be high-diagnostic too
//   but their dosages don't predict left-half assignments.

function makeTwoArrangementChunk() {
  const N_HOM = 20, N_HET = 20, N_HOM2 = 20;
  const samples = [];
  for (let i = 0; i < N_HOM;  i++) samples.push('a_' + i);
  for (let i = 0; i < N_HET;  i++) samples.push('b_' + i);
  for (let i = 0; i < N_HOM2; i++) samples.push('c_' + i);
  const N = samples.length;
  const N_SNPS = 200;

  const markers = [];
  const dosage  = [];
  for (let m = 0; m < N_SNPS; m++) {
    const bp = Math.round(10_000_000 + (m / N_SNPS) * 8_000_000);
    markers.push({ marker_id: 'snp_' + m, pos_bp: bp, missingness: 0 });
    const row = [];
    const inLeftHalf = m < N_SNPS / 2;
    if (inLeftHalf) {
      // Karyotype pattern matches the labels: a=Hom1 (dosage 0), b=Het (1),
      // c=Hom2 (2)
      for (let i = 0; i < N_HOM;  i++) row.push(0);
      for (let i = 0; i < N_HET;  i++) row.push(1);
      for (let i = 0; i < N_HOM2; i++) row.push(2);
    } else {
      // RIGHT half: scrambled — dosage doesn't track left labels.
      // a samples: 0..9 → 0, 10..19 → 2 (split half-and-half)
      // b samples: alternating 1
      // c samples: 0..9 → 2, 10..19 → 0
      // Hom1 mean dosage on RIGHT = (0*10 + 2*10) / 20 = 1
      // Hom2 mean dosage on RIGHT = (2*10 + 0*10) / 20 = 1
      // → diagnostic = |1 - 1| = 0, so SNPs in right half won't even be
      // included as "diagnostic" → bin scores in the right half are NaN
      // (insufficient diagnostic SNPs). Need another scrambling that keeps
      // RIGHT diagnostic but uncorrelated with LEFT.
      //
      // Better: in RIGHT half, define a fresh partition based on EVEN/ODD
      // sample index. Even samples → 2, odd → 0. Then RIGHT diagnostic
      // wrt the left labels:
      //   Hom1 (a_0..a_19): mean dosage = (10*2 + 10*0)/20 = 1 again
      // Same trap. To get high RIGHT-half diagnostic relative to LEFT
      // labels, we want a NEW partition that doesn't co-segregate with
      // (a, c). Solution: shift the partition by one sample position.
      // Specifically: a_0..a_9 + b_0..b_9 → 0; a_10..a_19 + b_10..b_19 +
      // c_0..c_19 → 2.
      // → On Hom1=a samples: half are 0, half are 2 → mean 1
      // → On Hom2=c samples: all 2 → mean 2
      // → diagnostic = |1 - 2| = 1.0  ✓
      // But the per-sample arrangement score for left-half labels' Hom1
      // (a) on the RIGHT bin will be: half of them are 0, half are 2
      // → mean within-Hom1 group = 1, but per-sample value is split.
      // So RIGHT-bin per-sample scores form a different cluster than
      // LEFT-bin per-sample scores → low cor between halves.  ✓

      // Polarity setup:
      // a_0..a_9 → 0; a_10..a_19 → 2; all b → 1; all c → 2
      for (let i = 0; i < N_HOM;  i++) row.push(i < 10 ? 0 : 2);
      for (let i = 0; i < N_HET;  i++) row.push(1);
      for (let i = 0; i < N_HOM2; i++) row.push(2);
    }
    dosage.push(row);
  }
  return { samples, markers, dosage };
}

group('runShelfLDTest — two arrangements (independent left/right structure)', () => {
  const chunk = makeTwoArrangementChunk();
  const invgt = makeInvgtMap(chunk, 20, 20, 20);
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    n_bins: 20,
  });
  ok('ok',  r.ok);
  // LEFT-half SNPs: highly diagnostic with the labels (perfect)
  // RIGHT-half SNPs: also diagnostic (mean a = 1, mean c = 2 → |diag| = 1.0)
  // The KEY question: does the per-sample arrangement score in left bins
  // correlate with right bins? It should NOT.
  ok('verdict = LIKELY_MULTIPLE_ARRANGEMENTS or AMBIGUOUS',
    r.summary.verdict === 'LIKELY_MULTIPLE_ARRANGEMENTS' ||
    r.summary.verdict === 'AMBIGUOUS_OR_PARTIAL_LD');
  ok('between < 0.7 (not single)', r.summary.between_halves_corr < 0.7);
  // within-A should be high (LEFT bins all agree)
  ok('within-A high', r.summary.within_first_half_corr > 0.7);
  // within-B should be high (RIGHT bins all agree among themselves)
  ok('within-B high', r.summary.within_second_half_corr > 0.7);
});

group('runShelfLDTest — error: no shelf SNPs', () => {
  const chunk = makeSingleInversionChunk();
  const invgt = makeInvgtMap(chunk, 20, 20, 20);
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 50_000_000, shelf_end_bp: 60_000_000,   // outside chunk
    n_bins: 20,
  });
  ok('not ok', !r.ok);
  ok('error mentions shelf', /shelf|SNP/i.test(r.error));
});

group('runShelfLDTest — error: too few Hom1', () => {
  const chunk = makeSingleInversionChunk();
  const invgt = new Map();
  invgt.set(chunk.samples[0], 'Hom1');   // only 1 Hom1
  invgt.set(chunk.samples[20], 'Het');
  invgt.set(chunk.samples[40], 'Hom2');
  invgt.set(chunk.samples[41], 'Hom2');
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    n_bins: 20,
  });
  ok('not ok', !r.ok);
  ok('error mentions group sizes', /Hom1|Hom2/i.test(r.error));
});

group('runShelfLDTest — handles H-system labels (H1/H1 alias)', () => {
  const chunk = makeSingleInversionChunk();
  const invgt = new Map();
  let i = 0;
  for (let k = 0; k < 20; k++) invgt.set(chunk.samples[i++], 'H1/H1');
  for (let k = 0; k < 20; k++) invgt.set(chunk.samples[i++], 'H1/H2');
  for (let k = 0; k < 20; k++) invgt.set(chunk.samples[i++], 'H2/H2');
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    n_bins: 20,
  });
  ok('ok with H-system labels', r.ok);
  eq('verdict = SINGLE_INVERSION', r.summary.verdict, 'SINGLE_INVERSION');
});

group('runShelfLDTest — NA dosage handling (-1 = NA)', () => {
  const chunk = makeSingleInversionChunk();
  // Set 5% of dosages to -1 (missing)
  for (let m = 0; m < chunk.markers.length; m++) {
    if (m % 20 === 0) {
      chunk.dosage[m][0] = -1;
      chunk.dosage[m][30] = -1;
    }
  }
  const invgt = makeInvgtMap(chunk, 20, 20, 20);
  const r = mod.runShelfLDTest({
    chunk, sample_invgt: invgt,
    shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    n_bins: 20,
  });
  ok('ok with NA',  r.ok);
  // Should still classify as SINGLE_INVERSION
  eq('still SINGLE_INVERSION', r.summary.verdict, 'SINGLE_INVERSION');
});

// =============================================================================
// makeShelfLDPanel — DOM construction smoke test
// =============================================================================

group('makeShelfLDPanel constructs DOM', () => {
  const panel = mod.makeShelfLDPanel({
    inputProvider: () => ({
      chunk: makeSingleInversionChunk(),
      sample_invgt: makeInvgtMap(makeSingleInversionChunk(), 20, 20, 20),
      shelf_start_bp: 10_000_000,
      shelf_end_bp:   18_000_000,
    }),
  });
  ok('panel created', !!panel);
  eq('data-popgen-shelf-ld root',
    panel.attributes['data-popgen-shelf-ld'], 'root');
  // Children: header, status, canvas, verdict
  ok('has 4 children',  panel.children.length === 4);
  // Run button is inside header (children[0])
  const header = panel.children[0];
  const runBtn = header.children.find(c =>
    c.attributes && c.attributes['data-popgen-shelf-ld'] === 'run');
  ok('run button found', !!runBtn);
});

await group('makeShelfLDPanel run button triggers compute', async () => {
  const provided = makeSingleInversionChunk();
  const invgt = makeInvgtMap(provided, 20, 20, 20);
  const panel = mod.makeShelfLDPanel({
    inputProvider: () => ({
      chunk: provided, sample_invgt: invgt,
      shelf_start_bp: 10_000_000, shelf_end_bp: 18_000_000,
    }),
  });
  const header = panel.children[0];
  const runBtn = header.children.find(c =>
    c.attributes && c.attributes['data-popgen-shelf-ld'] === 'run');
  // Simulate click (handler is async; resolves on next tick)
  runBtn.dispatchEvent({
    type: 'click', target: runBtn,
    preventDefault() {}, stopPropagation() {},
  });
  // Yield enough to let the async handler finish
  await new Promise(r => setTimeout(r, 20));
  await new Promise(r => setTimeout(r, 0));
  // Status should now reflect the result
  const status = panel.children[1];
  ok('status updated', /200 shelf SNPs/.test(status._text || status.textContent));
  const verdict = panel.children[3];
  ok('verdict set', verdict.attributes['data-verdict'] === 'SINGLE_INVERSION');
});

// =============================================================================
// Page 7 stratifier
// =============================================================================

group('makePage7AncestryStratifier — falls back to cohort delta12', () => {
  const stratifier = mod.makePage7AncestryStratifier();
  const s = {
    data: {
      windows: [
        { center_mb: 0.5 }, { center_mb: 1.5 }, { center_mb: 2.5 },
      ],
      ancestry_window: { delta12: [0.55, 0.62, 0.5] },
    },
  };
  const result = stratifier(s);
  ok('returns data',    !!result);
  eq('series length 1', result.series.length, 1);
  eq('values',          result.series[0].values, [0.55, 0.62, 0.5]);
});

group('makePage7AncestryStratifier — returns null when nothing available', () => {
  const stratifier = mod.makePage7AncestryStratifier();
  const r = stratifier({ data: null });
  ok('null', r == null);
});

group('makePage7AncestryStratifier — uses live ancestry response when present', () => {
  // Fake renderers stub
  const renderers = {
    adaptAncestryResponse: () => ({
      ancestry_q1_per_group: {
        mb: [0.5, 1.5],
        series: [
          { name: 'HOM1', values: [0.7, 0.6] },
          { name: 'HOM2', values: [0.3, 0.4] },
        ],
      },
    }),
  };
  globalThis.window.popgenRenderers = renderers;
  const stratifier = mod.makePage7AncestryStratifier();
  const s = {
    data: { windows: [], ancestry_window: null },
    popstatsLive: { lastAncestryResponse: { K: 2 } },
  };
  const result = stratifier(s);
  ok('returns data', !!result);
  eq('two series',   result.series.length, 2);
});

// ----- summary -------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);

})().catch(e => { console.error('UNCAUGHT', e); process.exit(1); });
