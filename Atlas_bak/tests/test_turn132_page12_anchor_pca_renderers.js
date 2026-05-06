// =============================================================================
// turn 132 — Slice 7: θπ envelope anchor strip + PC1×PC2 scatter renderers
// for #thAnchorStripPanel and #thPcaPanel.
//
// _drawThAnchorStripPanel reads state.data.theta_pi_envelopes.{l1, l2}
// and paints L1 (blue, top half) + L2 (green, bottom half) rectangles
// along the genomic Mb axis. Naming "anchor strip" mirrors page 1 for
// visual parity but the SEMANTICS differ: page 1 paints anchor-concord
// values, page 12 paints envelope geometry. Comment block in source
// documents the deliberate naming clash.
//
// _drawThPcaPanel reads state.data.theta_pi_local_pca.{pc1, pc2} (flat
// per-sample arrays — collapsed view) and paints a K-means colored
// scatter. Per-window override via pc1_by_window/pc2_by_window if
// state.cur is valid. Coloring cascade: cluster_labels_theta →
// candidate.locked_labels → grey (matches lines panel + hero strip).
//
// Deferred to later slices: window-trail, lasso, sign-flip rules,
// per-cluster ellipses, hit-testing.
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

ok('_drawThAnchorStripPanel defined',
   /function _drawThAnchorStripPanel\(\)/.test(html));
ok('_drawThPcaPanel defined',
   /function _drawThPcaPanel\(\)/.test(html));

ok('anchor strip reads theta_pi_envelopes.l1 / l2',
   /_drawThAnchorStripPanel[\s\S]{0,1500}theta_pi_envelopes[\s\S]{0,500}\.l1[\s\S]{0,500}\.l2/.test(html));
ok('PCA reads theta_pi_local_pca.pc1 / pc2',
   /_drawThPcaPanel[\s\S]{0,2000}theta_pi_local_pca[\s\S]{0,500}\.pc1[\s\S]{0,500}\.pc2/.test(html));

ok('PCA supports per-window override via pc1_by_window / pc2_by_window',
   /_drawThPcaPanel[\s\S]{0,3000}pc1_by_window[\s\S]{0,200}pc2_by_window/.test(html));

ok('PCA falls back to cluster_labels_theta for coloring',
   /_drawThPcaPanel[\s\S]{0,8000}cluster_labels_theta/.test(html));
ok('PCA falls back to candidate.locked_labels for coloring',
   /_drawThPcaPanel[\s\S]{0,8000}state\.candidate[\s\S]{0,200}locked_labels/.test(html));

ok('both renderers fire after sim_mat + Z renderers',
   /_drawThZPanel[\s\S]{0,500}_drawThAnchorStripPanel[\s\S]{0,300}_drawThPcaPanel/.test(html));

ok('anchor strip is vm-safe',
   /_drawThAnchorStripPanel[\s\S]{0,300}typeof document === 'undefined'/.test(html));
ok('PCA is vm-safe',
   /_drawThPcaPanel[\s\S]{0,300}typeof document === 'undefined'/.test(html));

ok('anchor strip documents naming-clash deliberation',
   /NAMING NOTE[\s\S]{0,800}page 12[\s\S]{0,300}envelope geometry/.test(html));
ok('PCA documents deferred features (window-trail, lasso, sign-flip)',
   /Window-trail[\s\S]{0,300}Sign-flip[\s\S]{0,300}Lasso[\s\S]{0,500}function _drawThPcaPanel\(\)/.test(html));

// ============================================================================
// 2. Behavioral — sandbox exec
// ============================================================================
console.log('\n=== 2. Behavioral — sandbox exec ===');

function makeCanvasMock(initialW, initialH) {
  const calls = [];
  const ctx = {
    clearRect:  (...a) => calls.push(['clearRect', ...a]),
    fillRect:   (...a) => calls.push(['fillRect', ...a]),
    strokeRect: (...a) => calls.push(['strokeRect', ...a]),
    fillText:   (...a) => calls.push(['fillText', ...a]),
    beginPath:  () => calls.push(['beginPath']),
    moveTo:     (...a) => calls.push(['moveTo', ...a]),
    lineTo:     (...a) => calls.push(['lineTo', ...a]),
    stroke:     () => calls.push(['stroke']),
    save:       () => calls.push(['save']),
    restore:    () => calls.push(['restore']),
    translate:  (...a) => calls.push(['translate', ...a]),
    rotate:     (...a) => calls.push(['rotate', ...a]),
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
    clientWidth: initialW, clientHeight: initialH,
    width: initialW, height: initialH,
    getContext: () => ctx,
    _calls: calls,
  };
}

function buildAnchorSandbox(theta_pi_envelopes, windows) {
  const main = makeCanvasMock(800, 24);
  return {
    document: {
      getElementById(id) {
        if (id === 'thAnchorStripCanvas') return main;
        return null;
      },
    },
    state: {
      data: Object.assign(
        { theta_pi_envelopes },
        windows ? { windows } : {}
      ),
    },
    _main: main,
  };
}

function buildPcaSandbox(theta_pi_local_pca, options) {
  options = options || {};
  const main = makeCanvasMock(600, 400);
  return {
    document: {
      getElementById(id) {
        if (id === 'thPcaCanvas') return main;
        return null;
      },
    },
    state: {
      data: Object.assign(
        { theta_pi_local_pca },
        options.cluster_labels_theta ? { cluster_labels_theta: options.cluster_labels_theta } : {}
      ),
      candidate: options.candidate || null,
      cur: options.cur,
    },
    _main: main,
  };
}

function runAnchorIn(sandbox) {
  const fnSrc = html.match(/function _drawThAnchorStripPanel\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document, state: sandbox.state,
    Number, Array, Math, Infinity, ArrayBuffer,
  });
  vm.runInContext(fnSrc + '\n_drawThAnchorStripPanel();', ctx);
}

function runPcaIn(sandbox) {
  const fnSrc = html.match(/function _drawThPcaPanel\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document, state: sandbox.state,
    Number, Array, Math, Infinity, ArrayBuffer,
  });
  vm.runInContext(fnSrc + '\n_drawThPcaPanel();', ctx);
}

// ─── anchor strip tests ─────────────────────────────────────────────────

// 2a. No envelopes → no-op
{
  const sandbox = buildAnchorSandbox(undefined);
  let crashed = false;
  try { runAnchorIn(sandbox); } catch (e) { crashed = true; }
  ok('2a: no envelopes does not crash', !crashed);
  ok('2a: no envelopes → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2b. Both l1 and l2 empty → no-op
{
  const sandbox = buildAnchorSandbox({ l1: [], l2: [] });
  let crashed = false;
  try { runAnchorIn(sandbox); } catch (e) { crashed = true; }
  ok('2b: empty envelopes does not crash', !crashed);
  ok('2b: empty envelopes → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2c. L1 + L2 envelopes with windows[] for axis range → paints rects
{
  const windows = [];
  for (let i = 0; i < 100; i++) {
    windows.push({ start_bp: i * 50000, center_mb: (i * 50000) / 1e6 });
  }
  const sandbox = buildAnchorSandbox({
    l1: [
      { start_bp: 1.0e6, end_bp: 1.5e6 },
      { start_bp: 3.0e6, end_bp: 3.8e6 },
    ],
    l2: [
      { start_bp: 1.1e6, end_bp: 1.4e6 },
      { start_bp: 3.1e6, end_bp: 3.7e6 },
      { start_bp: 4.0e6, end_bp: 4.2e6 },
    ],
  }, windows);
  runAnchorIn(sandbox);

  ok('2c: canvas cleared',
     sandbox._main._calls.some(c => c[0] === 'clearRect'));

  const fillRects = sandbox._main._calls.filter(c => c[0] === 'fillRect');
  ok('2c: 5 envelope rectangles painted (2 L1 + 3 L2)',
     fillRects.length === 5,
     'got ' + fillRects.length);

  const fillStyles = sandbox._main._calls
    .filter(c => c[0] === 'fillStyle')
    .map(c => c[1]);
  ok('2c: L1 blue color used',
     fillStyles.some(s => s.indexOf('48,116,200') !== -1));
  ok('2c: L2 green color used',
     fillStyles.some(s => s.indexOf('82,168,124') !== -1));

  ok('2c: L1/L2 side labels drawn',
     sandbox._main._calls.filter(c => c[0] === 'fillText').length >= 2);
}

// 2d. L1 only (no L2) → still paints
{
  const sandbox = buildAnchorSandbox({
    l1: [{ start_bp: 1e6, end_bp: 2e6 }],
    l2: [],
  });
  let crashed = false;
  try { runAnchorIn(sandbox); } catch (e) { crashed = true; }
  ok('2d: L1-only does not crash', !crashed);
  ok('2d: L1-only → at least one fillRect for envelope',
     sandbox._main._calls.filter(c => c[0] === 'fillRect').length >= 1);
}

// 2e. Falls back to envelope extents when windows[] absent
{
  const sandbox = buildAnchorSandbox({
    l1: [{ start_bp: 5e6, end_bp: 6e6 }],
    l2: [{ start_bp: 7e6, end_bp: 8e6 }],
  });
  let crashed = false;
  try { runAnchorIn(sandbox); } catch (e) { crashed = true; }
  ok('2e: no windows → falls back to envelope extents (no crash)', !crashed);
  ok('2e: still paints rectangles',
     sandbox._main._calls.filter(c => c[0] === 'fillRect').length >= 2);
}

// 2f. Idempotent
{
  const sandbox = buildAnchorSandbox({
    l1: [{ start_bp: 1e6, end_bp: 2e6 }],
    l2: [{ start_bp: 1.2e6, end_bp: 1.8e6 }],
  });
  runAnchorIn(sandbox);
  const cnt1 = sandbox._main._calls.filter(c => c[0] === 'clearRect').length;
  runAnchorIn(sandbox);
  const cnt2 = sandbox._main._calls.filter(c => c[0] === 'clearRect').length;
  ok('2f: idempotent — second draw cleared canvas again',
     cnt2 === cnt1 + 1);
}

// 2g. Envelope rows with missing start_bp / end_bp are skipped
{
  const sandbox = buildAnchorSandbox({
    l1: [
      { start_bp: 1e6, end_bp: 2e6 },
      { start_bp: 1e6 },                  // missing end_bp
      { end_bp: 2e6 },                    // missing start_bp
      { start_bp: 'foo', end_bp: 2e6 },   // non-numeric
    ],
    l2: [],
  });
  let crashed = false;
  try { runAnchorIn(sandbox); } catch (e) { crashed = true; }
  ok('2g: malformed envelope rows do not crash', !crashed);
  // Only 1 valid row should produce a fillRect for envelopes
  const fillRects = sandbox._main._calls.filter(c => c[0] === 'fillRect');
  ok('2g: only valid rows painted (1 of 4 L1 entries)',
     fillRects.length === 1);
}

// ─── PCA scatter tests ──────────────────────────────────────────────────

// 2h. No theta_pi_local_pca → no-op
{
  const sandbox = buildPcaSandbox(undefined);
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2h: no PCA data does not crash', !crashed);
  ok('2h: no PCA data → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2i. Missing pc1/pc2 → no-op
{
  const sandbox = buildPcaSandbox({ z: [1, 2, 3] });
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2i: missing pc1/pc2 does not crash', !crashed);
  ok('2i: missing pc1/pc2 → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2j. Valid pc1 + pc2 → scatter painted
{
  const nSamp = 60;
  const pc1 = new Array(nSamp);
  const pc2 = new Array(nSamp);
  for (let i = 0; i < nSamp; i++) {
    pc1[i] = (i % 3 === 0 ? -1 : 1) + (Math.random() - 0.5) * 0.1;
    pc2[i] = (i % 2 === 0 ? -1 : 1) + (Math.random() - 0.5) * 0.1;
  }
  const sandbox = buildPcaSandbox({ pc1, pc2 });
  runPcaIn(sandbox);

  ok('2j: PCA canvas cleared',
     sandbox._main._calls.some(c => c[0] === 'clearRect'));
  ok('2j: frame stroked',
     sandbox._main._calls.some(c => c[0] === 'strokeRect'));
  ok('2j: nSamp dots painted (fillRect calls)',
     sandbox._main._calls.filter(c => c[0] === 'fillRect').length === nSamp);
  ok('2j: axis labels painted',
     sandbox._main._calls.filter(c => c[0] === 'fillText').length >= 2);
}

// 2k. cluster_labels_theta coloring uses karyotype palette
{
  const nSamp = 6;
  const pc1 = [-1, -1, 0, 0, 1, 1];
  const pc2 = [-1, 1, -1, 1, -1, 1];
  const labels = [0, 1, 2, 0, 1, 2];
  const sandbox = buildPcaSandbox(
    { pc1, pc2 },
    { cluster_labels_theta: { labels } }
  );
  runPcaIn(sandbox);

  const fillStyles = sandbox._main._calls
    .filter(c => c[0] === 'fillStyle')
    .map(c => c[1]);

  ok('2k: HOM_REF blue used',
     fillStyles.some(s => s.indexOf('58,125,222') !== -1));
  ok('2k: HET orange used',
     fillStyles.some(s => s.indexOf('217,120,66') !== -1));
  ok('2k: HOM_INV purple used',
     fillStyles.some(s => s.indexOf('124,74,217') !== -1));
}

// 2l. candidate.locked_labels fallback
{
  const nSamp = 4;
  const sandbox = buildPcaSandbox(
    { pc1: [-1, 0, 1, 0.5], pc2: [-1, 1, -1, 0.2] },
    { candidate: { locked_labels: [0, 1, 2, 1] } }
  );
  runPcaIn(sandbox);

  const fillStyles = sandbox._main._calls
    .filter(c => c[0] === 'fillStyle')
    .map(c => c[1]);

  ok('2l: locked_labels fallback uses karyotype palette',
     fillStyles.some(s => s.indexOf('58,125,222') !== -1) &&
     fillStyles.some(s => s.indexOf('217,120,66') !== -1) &&
     fillStyles.some(s => s.indexOf('124,74,217') !== -1));
}

// 2m. No labels → all dots grey
{
  const nSamp = 4;
  const sandbox = buildPcaSandbox(
    { pc1: [-1, 0, 1, 0.5], pc2: [-1, 1, -1, 0.2] }
  );
  runPcaIn(sandbox);

  const fillStyles = sandbox._main._calls
    .filter(c => c[0] === 'fillStyle')
    .map(c => c[1]);

  ok('2m: no labels → grey dots used (154,163,173)',
     fillStyles.some(s => s.indexOf('154,163,173') !== -1));
  ok('2m: no labels → no karyotype palette used',
     !fillStyles.some(s =>
       s.indexOf('58,125,222') !== -1 ||
       s.indexOf('217,120,66') !== -1 ||
       s.indexOf('124,74,217') !== -1));
}

// 2n. Per-window override via pc1_by_window / pc2_by_window when state.cur valid
{
  const nSamp = 4;
  const pc1_by_window = [
    [-1, -1, 1, 1],
    [-2, -2, 2, 2],   // current window
    [-3, -3, 3, 3],
  ];
  const pc2_by_window = [
    [-1, 1, -1, 1],
    [-2, 2, -2, 2],   // current window
    [-3, 3, -3, 3],
  ];
  const sandbox = buildPcaSandbox(
    {
      pc1: [-100, -100, 100, 100],   // collapsed (should NOT be used)
      pc2: [-100, 100, -100, 100],
      pc1_by_window,
      pc2_by_window,
    },
    { cur: 1 }
  );
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2n: per-window override works', !crashed);
  ok('2n: per-window override → 4 dots painted',
     sandbox._main._calls.filter(c => c[0] === 'fillRect').length === nSamp);
}

// 2o. Invalid state.cur → falls back to collapsed pc1 / pc2
{
  const sandbox = buildPcaSandbox(
    {
      pc1: [-1, 0, 1, 0.5],
      pc2: [-1, 1, -1, 0.2],
      pc1_by_window: [[1, 2, 3, 4]],   // only 1 window; cur=999 invalid
      pc2_by_window: [[1, 2, 3, 4]],
    },
    { cur: 999 }
  );
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2o: invalid state.cur falls back to collapsed view',
     !crashed && sandbox._main._calls.filter(c => c[0] === 'fillRect').length === 4);
}

// 2p. NaN values in pc1 / pc2 are filtered (no crash)
{
  const sandbox = buildPcaSandbox({
    pc1: [-1, NaN, 1, 0.5],
    pc2: [-1, 1, NaN, 0.2],
  });
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2p: NaN values do not crash', !crashed);
}

// 2q. Float32Array pc1 / pc2 supported
{
  const nSamp = 6;
  const pc1 = new Float32Array(nSamp);
  const pc2 = new Float32Array(nSamp);
  for (let i = 0; i < nSamp; i++) {
    pc1[i] = (i - 3) * 0.5;
    pc2[i] = (i % 2) - 0.5;
  }
  const sandbox = buildPcaSandbox({ pc1, pc2 });
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2q: Float32Array pc1/pc2 works', !crashed);
  ok('2q: Float32Array pc1/pc2 → dots painted',
     sandbox._main._calls.filter(c => c[0] === 'fillRect').length === nSamp);
}

// 2r. nSamp < 2 → no-op
{
  const sandbox = buildPcaSandbox({ pc1: [0.5], pc2: [0.5] });
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2r: nSamp=1 does not crash', !crashed);
  ok('2r: nSamp=1 → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2s. Mismatched pc1/pc2 lengths → no-op
{
  const sandbox = buildPcaSandbox({ pc1: [0.1, 0.2, 0.3], pc2: [0.5, 0.6] });
  let crashed = false;
  try { runPcaIn(sandbox); } catch (e) { crashed = true; }
  ok('2s: mismatched lengths does not crash', !crashed);
  ok('2s: mismatched lengths → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// ============================================================================
// 3. Earlier slices still intact
// ============================================================================
console.log('\n=== 3. Earlier slices still intact ===');

ok('Slice 1: #thAnchorStripPanel still in DOM', /id="thAnchorStripPanel"/.test(html));
ok('Slice 1: #thAnchorStripCanvas still in DOM', /id="thAnchorStripCanvas"/.test(html));
ok('Slice 1: #thPcaPanel still in DOM', /id="thPcaPanel"/.test(html));
ok('Slice 1: #thPcaCanvas still in DOM', /id="thPcaCanvas"/.test(html));
ok('Slice 2: visibility wires #thAnchorStripPanel to hasEnvelopes',
   /showHide\('thAnchorStripPanel',\s*hasEnvelopes/.test(html));
ok('Slice 2: visibility wires #thPcaPanel to hasLocalPCA',
   /showHide\('thPcaPanel',\s*hasLocalPCA/.test(html));
ok('Slice 6: _drawThSimMatPanel still defined',
   /function _drawThSimMatPanel\(\)/.test(html));
ok('Slice 6: _drawThZPanel still defined',
   /function _drawThZPanel\(\)/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
