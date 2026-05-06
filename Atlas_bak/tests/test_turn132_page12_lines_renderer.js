// =============================================================================
// turn 132 — Slice 5: per-sample θπ lines renderer for #thLinesPanel.
//
// _drawThLinesPanel() reads state.data.theta_pi_per_window and paints
// 226 polylines into a single canvas inside #thLinesCanvasContainer.
//
// Supports BOTH canonical shapes from the existing detector:
//   Shape A (legacy): data.windows[w].theta = [<n_samples values>]
//   Shape B (current): data.theta_pi_per_window.{samples, windows, values}
//
// Coloring:
//   - state.data.cluster_labels_theta.labels (preferred), or
//   - state.candidate.locked_labels (fallback), or
//   - neutral grey (no labels available)
//
// This slice is the minimum viable lines mirror — single source, no lasso,
// no tracked-sample thickening, no K=3/K=6 K-means coloring, no caching.
// Those are page-1-specific control state and ship in later slices if
// needed.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Source-level checks
// ============================================================================
console.log('\n=== 1. Source-level checks ===');

ok('_drawThLinesPanel function defined',
   /function _drawThLinesPanel\(\)/.test(html));

ok('renderer reads from state.data.theta_pi_per_window',
   /_drawThLinesPanel[\s\S]{0,1500}state\.data\.theta_pi_per_window/.test(html));

ok('renderer is vm-safe (early-returns when document is undefined)',
   /_drawThLinesPanel[\s\S]{0,300}typeof document === 'undefined'/.test(html));

ok('renderer fires after _drawThCusumHero in onDataLoad path',
   /_drawThCusumHero[\s\S]{0,500}_drawThLinesPanel/.test(html));

ok('Shape A support documented (legacy windows[w].theta arrays)',
   /Shape A[\s\S]{0,200}embedded per-window/.test(html));

ok('Shape B support documented (top-level theta_pi_per_window)',
   /Shape B[\s\S]{0,200}top-level layer/.test(html));

ok('renderer falls back to cluster_labels_theta for coloring',
   /_drawThLinesPanel[\s\S]{0,8000}cluster_labels_theta/.test(html));

ok('renderer falls back to candidate.locked_labels for coloring',
   /_drawThLinesPanel[\s\S]{0,8000}state\.candidate[\s\S]{0,200}locked_labels/.test(html));

ok('lasso/tracked/K-means deferral comments present',
   /Deferred to later slices[\s\S]{0,400}Lasso[\s\S]{0,200}Tracked-sample[\s\S]{0,200}K-means/.test(html));

// ============================================================================
// 2. Behavioral — sandbox exec with synthetic data
// ============================================================================
console.log('\n=== 2. Behavioral — sandbox exec ===');

function makeCanvasMock(initialW, initialH) {
  const calls = [];
  const ctx = {
    clearRect: (...a) => calls.push(['clearRect', ...a]),
    fillRect:  (...a) => calls.push(['fillRect',  ...a]),
    fillText:  (...a) => calls.push(['fillText',  ...a]),
    beginPath: ()     => calls.push(['beginPath']),
    moveTo:    (...a) => calls.push(['moveTo',    ...a]),
    lineTo:    (...a) => calls.push(['lineTo',    ...a]),
    stroke:    ()     => calls.push(['stroke']),
    set fillStyle(v)  { calls.push(['fillStyle', v]); },
    get fillStyle()   { return ''; },
    set strokeStyle(v){ calls.push(['strokeStyle', v]); },
    get strokeStyle() { return ''; },
    set lineWidth(v)  { calls.push(['lineWidth', v]); },
    get lineWidth()   { return 1; },
    set font(v)       { calls.push(['font', v]); },
    get font()        { return ''; },
    set textAlign(v)  { calls.push(['textAlign', v]); },
    get textAlign()   { return ''; },
  };
  return {
    clientWidth:  initialW,
    clientHeight: initialH,
    width:        initialW,
    height:       initialH,
    getContext: () => ctx,
    style: {},
    _calls: calls,
  };
}

function makeContainer(canvas) {
  return {
    children: canvas ? [canvas] : [],
    querySelector(sel) {
      if (sel !== 'canvas') return null;
      return this.children.find(c => c && c._isCanvas) || null;
    },
    appendChild(c) { c._isCanvas = true; this.children.push(c); },
  };
}

function buildSandbox(theta_pi_per_window, options) {
  options = options || {};
  const canvas = options.preExisting ? makeCanvasMock(800, 180) : null;
  if (canvas) canvas._isCanvas = true;
  const container = makeContainer(canvas);

  const elements = { thLinesCanvasContainer: container };

  return {
    document: {
      elements,
      createElement(tag) {
        if (tag !== 'canvas') return null;
        const c = makeCanvasMock(800, 180);
        return c;
      },
      getElementById(id) { return this.elements[id] || null; },
    },
    state: {
      data: Object.assign(
        { theta_pi_per_window },
        options.cluster_labels_theta ? { cluster_labels_theta: options.cluster_labels_theta } : {}
      ),
      candidate: options.candidate || null,
    },
    _container: container,
    _preCanvas: canvas,
  };
}

function runFnIn(sandbox) {
  const fnSrc = html.match(/function _drawThLinesPanel\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document,
    state:    sandbox.state,
    Number:   Number,
    Array:    Array,
    Math:     Math,
    Infinity: Infinity,
    NaN:      NaN,
    isNaN:    isNaN,
  });
  vm.runInContext(fnSrc + '\n_drawThLinesPanel();', ctx);
}

function buildShapeBData(nWin, nSamp) {
  const samples = [];
  const windows = [];
  for (let w = 0; w < nWin; w++) {
    windows.push({ start_bp: w * 50000, center_mb: (w * 50000) / 1e6 });
  }
  for (let s = 0; s < nSamp; s++) {
    const tp = new Array(nWin);
    for (let w = 0; w < nWin; w++) {
      // synthetic: sample s, window w → 0.001 + s*0.0001 + sin(w/5)*0.0005
      tp[w] = 0.001 + s * 0.0001 + Math.sin(w / 5) * 0.0005;
    }
    samples.push({ sample_id: 'S' + s, theta_pi: tp });
  }
  return { samples, windows };
}

// 2a. No theta_pi_per_window data → renderer no-ops
{
  const sandbox = buildSandbox(undefined);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2a: no data does not crash', !crashed);
  ok('2a: no data → no canvas created',
     sandbox._container.children.length === 0);
}

// 2b. Shape B with samples + windows → renderer creates canvas, paints
{
  const data = buildShapeBData(60, 10);
  const sandbox = buildSandbox(data);
  runFnIn(sandbox);

  ok('2b: Shape B data → canvas created in container',
     sandbox._container.children.length === 1);

  const canvas = sandbox._container.children[0];
  ok('2b: created canvas has correct id',
     canvas.id === 'thLinesCanvas');

  ok('2b: canvas was cleared',
     canvas._calls.some(c => c[0] === 'clearRect'));

  // 10 samples × stroke per sample (in pass 0 since no labels = all ungrouped)
  const strokeCalls = canvas._calls.filter(c => c[0] === 'stroke').length;
  ok('2b: at least 10 polyline strokes (one per sample)', strokeCalls >= 10,
     'got ' + strokeCalls + ' stroke calls');
}

// 2c. Pre-existing canvas reused, not re-created
{
  const data = buildShapeBData(30, 5);
  const sandbox = buildSandbox(data, { preExisting: true });
  runFnIn(sandbox);

  ok('2c: pre-existing canvas reused (no second canvas)',
     sandbox._container.children.length === 1);
  ok('2c: pre-existing canvas painted (cleared)',
     sandbox._preCanvas._calls.some(c => c[0] === 'clearRect'));
}

// 2d. Empty samples / empty windows → no-op
{
  const sandbox = buildSandbox({ samples: [], windows: [] });
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2d: empty data does not crash', !crashed);
  ok('2d: empty data → no canvas created',
     sandbox._container.children.length === 0);
}

// 2e. Single window (n<2) → no-op
{
  const sandbox = buildSandbox(buildShapeBData(1, 5));
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2e: single window does not crash', !crashed);
  ok('2e: single window → no canvas created',
     sandbox._container.children.length === 0);
}

// 2f. cluster_labels_theta coloring — uses HOM_REF/HET/HOM_INV palette
{
  const data = buildShapeBData(30, 5);
  const labels = [0, 1, 2, 0, 1];   // 5 samples: REF, HET, INV, REF, HET
  const sandbox = buildSandbox(data, {
    cluster_labels_theta: { labels },
  });
  runFnIn(sandbox);

  const canvas = sandbox._container.children[0];
  const fillStyles = canvas._calls
    .filter(c => c[0] === 'strokeStyle')
    .map(c => c[1]);

  ok('2f: HOM_REF (label 0) blue color used',
     fillStyles.some(s => s.indexOf('58,125,222') !== -1));
  ok('2f: HET (label 1) orange color used',
     fillStyles.some(s => s.indexOf('217,120,66') !== -1));
  ok('2f: HOM_INV (label 2) purple color used',
     fillStyles.some(s => s.indexOf('124,74,217') !== -1));
}

// 2g. candidate.locked_labels fallback when cluster_labels_theta absent
{
  const data = buildShapeBData(30, 4);
  const candidate = {
    locked_labels: [0, 1, 2, 1],
  };
  const sandbox = buildSandbox(data, { candidate });
  runFnIn(sandbox);

  const canvas = sandbox._container.children[0];
  const fillStyles = canvas._calls
    .filter(c => c[0] === 'strokeStyle')
    .map(c => c[1]);

  ok('2g: locked_labels fallback uses karyotype palette',
     fillStyles.some(s => s.indexOf('58,125,222') !== -1) &&
     fillStyles.some(s => s.indexOf('217,120,66') !== -1) &&
     fillStyles.some(s => s.indexOf('124,74,217') !== -1));
}

// 2h. No labels at all → all samples drawn in neutral grey
{
  const data = buildShapeBData(30, 4);
  const sandbox = buildSandbox(data);
  runFnIn(sandbox);

  const canvas = sandbox._container.children[0];
  const fillStyles = canvas._calls
    .filter(c => c[0] === 'strokeStyle')
    .map(c => c[1]);

  ok('2h: no labels → grey palette used (154,163,173)',
     fillStyles.some(s => s.indexOf('154,163,173') !== -1));
  ok('2h: no labels → no karyotype palette colors used',
     !fillStyles.some(s =>
       s.indexOf('58,125,222') !== -1 ||
       s.indexOf('217,120,66') !== -1 ||
       s.indexOf('124,74,217') !== -1));
}

// 2i. NaN values in theta_pi break polylines, don't crash
{
  const data = buildShapeBData(30, 3);
  // Inject NaN into sample 0 at windows 5..10
  for (let w = 5; w <= 10; w++) data.samples[0].theta_pi[w] = NaN;
  // And one entire missing array for sample 1
  data.samples[1].theta_pi = null;

  const sandbox = buildSandbox(data);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2i: NaN values + null arrays do not crash', !crashed);
}

// 2j. Idempotency: calling twice yields same final state (no canvas leaks)
{
  const data = buildShapeBData(30, 5);
  const sandbox = buildSandbox(data);
  runFnIn(sandbox);
  const childCountAfter1 = sandbox._container.children.length;
  runFnIn(sandbox);
  const childCountAfter2 = sandbox._container.children.length;

  ok('2j: idempotent — second call does not create extra canvas',
     childCountAfter1 === 1 && childCountAfter2 === 1);
}

// 2k. Shape A (legacy): data.windows[w].theta arrays
{
  const nWin = 30, nSamp = 5;
  const windows = [];
  for (let w = 0; w < nWin; w++) {
    const theta = new Array(nSamp);
    for (let s = 0; s < nSamp; s++) {
      theta[s] = 0.001 + s * 0.0001 + Math.sin(w / 5) * 0.0005;
    }
    windows.push({
      start_bp: w * 50000,
      center_mb: (w * 50000) / 1e6,
      theta,
    });
  }

  // Shape A: state.data.windows[*].theta WITHOUT theta_pi_per_window key.
  const canvas = null;
  const container = makeContainer(canvas);
  const sandbox = {
    document: {
      elements: { thLinesCanvasContainer: container },
      createElement(tag) {
        if (tag !== 'canvas') return null;
        return makeCanvasMock(800, 180);
      },
      getElementById(id) { return this.elements[id] || null; },
    },
    state: {
      data: {
        // No theta_pi_per_window key; only top-level windows[]
        windows,
        samples: Array.from({ length: nSamp }, (_, i) => ({ sample_id: 'S' + i })),
      },
      candidate: null,
    },
    _container: container,
  };

  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2k: Shape A (legacy windows[w].theta) does not crash', !crashed);
  ok('2k: Shape A → canvas created', sandbox._container.children.length === 1);
  if (sandbox._container.children.length === 1) {
    const cv = sandbox._container.children[0];
    ok('2k: Shape A → canvas painted',
       cv._calls.some(c => c[0] === 'clearRect'));
    ok('2k: Shape A → at least 5 polyline strokes',
       cv._calls.filter(c => c[0] === 'stroke').length >= 5);
  }
}

// 2l. Defensive: invalid label values (out of range) → fall through to grey
{
  const data = buildShapeBData(30, 4);
  const sandbox = buildSandbox(data, {
    cluster_labels_theta: { labels: [99, -5, 'bogus', null] },
  });
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2l: invalid label values do not crash', !crashed);
}

// ============================================================================
// 3. Slice 1-3 still intact
// ============================================================================
console.log('\n=== 3. Earlier slices still intact ===');

ok('Slice 1: #thLinesPanel still in DOM',
   /id="thLinesPanel"/.test(html));
ok('Slice 1: #thLinesCanvasContainer still in DOM',
   /id="thLinesCanvasContainer"/.test(html));
ok('Slice 2: visibility wiring for #thLinesPanel still in source',
   /showHide\('thLinesPanel',\s*hasPerWindow/.test(html));
ok('Slice 3: _drawThCusumHero still defined',
   /function _drawThCusumHero\(\)/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
