// =============================================================================
// turn 132 — Slice 6: θπ sim_mat + |Z| renderers for #thSimPanel and #thZPanel.
//
// Both consume state.data.theta_pi_local_pca:
//   - _drawThSimMatPanel reads .sim_mat (flat row-major n×n) and .n_windows_thumb
//   - _drawThZPanel reads .z (per-window |Z| vector)
//
// Minimum-viable mirrors of page 1's drawSim and drawZ:
//   - Centered square sim_mat heatmap, palette matches simColor when global,
//     diagonal painted yellow for orientation
//   - Single-line |Z| waveform, dim |Z|=2 reference line, Mb gridlines
//
// Deferred (NOT in this slice): L1/L2 overlays, click-to-jump,
// PDF-style upper/lower split, candidate strip, settings popover.
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

ok('_drawThSimMatPanel defined',
   /function _drawThSimMatPanel\(\)/.test(html));
ok('_drawThZPanel defined',
   /function _drawThZPanel\(\)/.test(html));

ok('sim renderer reads theta_pi_local_pca.sim_mat',
   /_drawThSimMatPanel[\s\S]{0,1500}theta_pi_local_pca[\s\S]{0,200}sim_mat/.test(html));
ok('Z renderer reads theta_pi_local_pca.z',
   /_drawThZPanel[\s\S]{0,1500}theta_pi_local_pca[\s\S]{0,200}\.z\b/.test(html));

ok('sim renderer reuses page-1 simColor palette when available',
   /_drawThSimMatPanel[\s\S]{0,3500}typeof simColor === 'function'/.test(html));
ok('Z renderer plots |Z| (Math.abs)',
   /_drawThZPanel[\s\S]{0,3500}Math\.abs\(z\[i\]\)/.test(html));

ok('both renderers fire after _drawThLinesPanel',
   /_drawThLinesPanel[\s\S]{0,500}_drawThSimMatPanel[\s\S]{0,300}_drawThZPanel/.test(html));

ok('sim renderer is vm-safe (early-returns when document is undefined)',
   /_drawThSimMatPanel[\s\S]{0,300}typeof document === 'undefined'/.test(html));
ok('Z renderer is vm-safe',
   /_drawThZPanel[\s\S]{0,300}typeof document === 'undefined'/.test(html));

ok('sim renderer documents L1/L2 overlay deferral',
   /Deferred to later slices[\s\S]{0,500}L1[\s\S]{0,80}L2[\s\S]{0,200}envelope/.test(html));
ok('Z renderer documents L1/L2 zone bar deferral',
   /L1[^\n]{0,80}\/[^\n]{0,80}L2[^\n]{0,80}zone bar[\s\S]{0,500}function _drawThZPanel\(\)/.test(html));

// ============================================================================
// 2. Behavioral — sandbox exec
// ============================================================================
console.log('\n=== 2. Behavioral — sandbox exec ===');

function makeCanvasMock(initialW, initialH, withCreateImageData) {
  const calls = [];
  const ctx = {
    clearRect:    (...a) => calls.push(['clearRect', ...a]),
    fillRect:     (...a) => calls.push(['fillRect', ...a]),
    strokeRect:   (...a) => calls.push(['strokeRect', ...a]),
    fillText:     (...a) => calls.push(['fillText', ...a]),
    beginPath:    () => calls.push(['beginPath']),
    moveTo:       (...a) => calls.push(['moveTo', ...a]),
    lineTo:       (...a) => calls.push(['lineTo', ...a]),
    stroke:       () => calls.push(['stroke']),
    drawImage:    (...a) => calls.push(['drawImage', a[0] && a[0].width, a[1], a[2], a[3], a[4]]),
    putImageData: (...a) => calls.push(['putImageData']),
    createImageData: withCreateImageData
      ? (n, m) => ({ data: new Uint8ClampedArray(n * m * 4), width: n, height: m })
      : undefined,
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
    set imageSmoothingEnabled(v) { calls.push(['imageSmoothingEnabled', v]); },
    get imageSmoothingEnabled()  { return false; },
  };
  return {
    clientWidth: initialW, clientHeight: initialH,
    width: initialW, height: initialH,
    getContext: () => ctx,
    _calls: calls,
  };
}

function makeOffCanvasMock() {
  const calls = [];
  const ctx = {
    putImageData: (...a) => calls.push(['putImageData']),
  };
  return {
    width: 0, height: 0,
    getContext: () => ctx,
    _calls: calls,
  };
}

function buildSimSandbox(theta_pi_local_pca) {
  const main = makeCanvasMock(800, 600, /*withCreateImageData*/ true);
  return {
    document: {
      getElementById(id) {
        if (id === 'thSimCanvas') return main;
        return null;
      },
      createElement(tag) {
        if (tag !== 'canvas') return null;
        return makeOffCanvasMock();
      },
    },
    state: { data: { theta_pi_local_pca } },
    _main: main,
  };
}

function buildZSandbox(theta_pi_local_pca, windows) {
  const main = makeCanvasMock(800, 140, false);
  return {
    document: {
      getElementById(id) {
        if (id === 'thZCanvas') return main;
        return null;
      },
    },
    state: {
      data: Object.assign(
        { theta_pi_local_pca },
        windows ? { windows } : {}
      ),
    },
    _main: main,
  };
}

function runSimIn(sandbox) {
  const fnSrc = html.match(/function _drawThSimMatPanel\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document, state: sandbox.state,
    Number, Array, Math, Infinity, ArrayBuffer,
  });
  vm.runInContext(fnSrc + '\n_drawThSimMatPanel();', ctx);
}

function runZIn(sandbox) {
  const fnSrc = html.match(/function _drawThZPanel\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document, state: sandbox.state,
    Number, Array, Math, Infinity, ArrayBuffer,
  });
  vm.runInContext(fnSrc + '\n_drawThZPanel();', ctx);
}

// ─── sim_mat tests ───────────────────────────────────────────────────────

// 2a. No theta_pi_local_pca → no-op
{
  const sandbox = buildSimSandbox(undefined);
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2a: no theta_pi_local_pca does not crash', !crashed);
  ok('2a: no data → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2b. Missing sim_mat → no-op
{
  const sandbox = buildSimSandbox({ z: [1, 2, 3] });  // has z but no sim_mat
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; }
  ok('2b: missing sim_mat does not crash', !crashed);
  ok('2b: missing sim_mat → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2c. Valid sim_mat (10×10) → canvas cleared + drawImage called
{
  const N = 10;
  const sim = new Array(N * N);
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      sim[i * N + j] = (i === j) ? 1.0 : Math.exp(-Math.abs(i - j) / 3);
    }
  }
  const sandbox = buildSimSandbox({
    sim_mat: sim,
    n_windows_thumb: N,
    z: new Array(N).fill(0),
  });
  runSimIn(sandbox);

  ok('2c: canvas cleared',
     sandbox._main._calls.some(c => c[0] === 'clearRect'));
  ok('2c: heatmap blitted via drawImage',
     sandbox._main._calls.some(c => c[0] === 'drawImage'));
  ok('2c: frame stroked (strokeRect call)',
     sandbox._main._calls.some(c => c[0] === 'strokeRect'));
  ok('2c: state._thSimGeom persisted',
     sandbox.state._thSimGeom &&
     Number.isFinite(sandbox.state._thSimGeom.x0) &&
     Number.isFinite(sandbox.state._thSimGeom.side));
}

// 2d. n_windows_thumb absent → infers N from sqrt(sim.length)
{
  const N = 8;
  const sim = new Array(N * N);
  for (let k = 0; k < sim.length; k++) sim[k] = Math.random();
  const sandbox = buildSimSandbox({ sim_mat: sim, z: new Array(N).fill(0) });
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; }
  ok('2d: missing n_windows_thumb → infers N from sqrt(length)', !crashed);
  ok('2d: heatmap painted', sandbox._main._calls.some(c => c[0] === 'drawImage'));
}

// 2e. Non-square sim.length (sqrt isn't integer) and no n_windows_thumb → no-op
{
  const sandbox = buildSimSandbox({ sim_mat: new Array(50), z: [] });
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; }
  ok('2e: non-square sim.length → no crash', !crashed);
  ok('2e: non-square sim.length → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2f. Idempotent — second call yields same drawImage count
{
  const N = 6;
  const sim = new Array(N * N);
  for (let k = 0; k < sim.length; k++) sim[k] = Math.random();
  const sandbox = buildSimSandbox({
    sim_mat: sim, n_windows_thumb: N, z: new Array(N).fill(0),
  });
  runSimIn(sandbox);
  const cnt1 = sandbox._main._calls.filter(c => c[0] === 'drawImage').length;
  runSimIn(sandbox);
  const cnt2 = sandbox._main._calls.filter(c => c[0] === 'drawImage').length;
  ok('2f: idempotent — second draw adds another drawImage (no leaks)',
     cnt2 === cnt1 + 1);
}

// 2g. NaN values in sim_mat don't crash
{
  const N = 5;
  const sim = new Array(N * N);
  for (let k = 0; k < sim.length; k++) sim[k] = (k % 7 === 0) ? NaN : Math.random();
  const sandbox = buildSimSandbox({
    sim_mat: sim, n_windows_thumb: N, z: new Array(N).fill(0),
  });
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; }
  ok('2g: NaN in sim_mat does not crash', !crashed);
}

// 2h. Constant sim_mat (no range) → renderer no-ops gracefully
{
  const N = 5;
  const sim = new Array(N * N).fill(1.0);
  const sandbox = buildSimSandbox({
    sim_mat: sim, n_windows_thumb: N, z: new Array(N).fill(0),
  });
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; }
  ok('2h: constant sim_mat (zero range) does not crash', !crashed);
}

// ─── Z waveform tests ────────────────────────────────────────────────────

// 2i. No data → no-op
{
  const sandbox = buildZSandbox(undefined);
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; }
  ok('2i: Z renderer with no data does not crash', !crashed);
  ok('2i: no data → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2j. Empty z → no-op (length < 2)
{
  const sandbox = buildZSandbox({ z: [], sim_mat: [] });
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; }
  ok('2j: empty z does not crash', !crashed);
  ok('2j: empty z → canvas not painted',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2k. Valid z (50 windows) → canvas painted + waveform stroked
{
  const z = new Array(50);
  for (let i = 0; i < 50; i++) {
    z[i] = Math.abs(Math.sin(i / 5)) * 3;
  }
  const sandbox = buildZSandbox({ z, sim_mat: [] });
  runZIn(sandbox);

  ok('2k: Z canvas cleared',
     sandbox._main._calls.some(c => c[0] === 'clearRect'));
  ok('2k: waveform stroke calls present',
     sandbox._main._calls.some(c => c[0] === 'stroke'));
  ok('2k: |Z|=2 reference line drawn (multiple stroke calls beyond gridlines)',
     sandbox._main._calls.filter(c => c[0] === 'stroke').length >= 2);
}

// 2l. Z with windows[].center_mb available → uses Mb axis
{
  const z = new Array(20);
  for (let i = 0; i < 20; i++) z[i] = Math.abs(Math.sin(i / 3));
  const windows = [];
  for (let i = 0; i < 20; i++) {
    windows.push({ start_bp: i * 100000, center_mb: (i * 100000) / 1e6 });
  }
  const sandbox = buildZSandbox({ z, sim_mat: [] }, windows);
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2l: Z with Mb-bearing windows does not crash', !crashed);
  ok('2l: Z waveform painted',
     sandbox._main._calls.some(c => c[0] === 'stroke'));
}

// 2m. NaN in z breaks the polyline gracefully
{
  const z = new Array(30);
  for (let i = 0; i < 30; i++) z[i] = (i === 15) ? NaN : Math.abs(Math.sin(i / 3));
  const sandbox = buildZSandbox({ z, sim_mat: [] });
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; }
  ok('2m: NaN in z does not crash', !crashed);
}

// 2n. All-zero z → renderer detects degenerate range, no-ops gracefully
{
  const sandbox = buildZSandbox({ z: new Array(20).fill(0), sim_mat: [] });
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; }
  ok('2n: all-zero z does not crash', !crashed);
}

// 2o. Idempotent — second call reuses the same canvas
{
  const z = new Array(30).fill(0).map((_, i) => i / 10);
  const sandbox = buildZSandbox({ z, sim_mat: [] });
  runZIn(sandbox);
  runZIn(sandbox);
  // Two clearRect calls confirm both passes ran
  ok('2o: Z renderer idempotent — both calls cleared canvas',
     sandbox._main._calls.filter(c => c[0] === 'clearRect').length === 2);
}

// 2p. Typed array (Float32Array) for z works
{
  const z = new Float32Array(30);
  for (let i = 0; i < 30; i++) z[i] = Math.abs(Math.sin(i / 3));
  const sandbox = buildZSandbox({ z, sim_mat: [] });
  let crashed = false;
  try { runZIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2p: Float32Array z works', !crashed);
  ok('2p: Float32Array z → waveform painted',
     sandbox._main._calls.some(c => c[0] === 'stroke'));
}

// 2q. Typed array (Float32Array) for sim_mat works
{
  const N = 6;
  const sim = new Float32Array(N * N);
  for (let k = 0; k < sim.length; k++) sim[k] = Math.random();
  const sandbox = buildSimSandbox({
    sim_mat: sim, n_windows_thumb: N, z: new Float32Array(N),
  });
  let crashed = false;
  try { runSimIn(sandbox); } catch (e) { crashed = true; console.error(e); }
  ok('2q: Float32Array sim_mat works', !crashed);
  ok('2q: Float32Array sim_mat → drawImage called',
     sandbox._main._calls.some(c => c[0] === 'drawImage'));
}

// ============================================================================
// 3. Earlier slices still intact
// ============================================================================
console.log('\n=== 3. Earlier slices still intact ===');

ok('Slice 1: #thSimPanel still in DOM', /id="thSimPanel"/.test(html));
ok('Slice 1: #thSimCanvas still in DOM', /id="thSimCanvas"/.test(html));
ok('Slice 1: #thZPanel still in DOM', /id="thZPanel"/.test(html));
ok('Slice 1: #thZCanvas still in DOM', /id="thZCanvas"/.test(html));
ok('Slice 2: visibility wires #thSimPanel to hasLocalPCA',
   /showHide\('thSimPanel',\s*hasLocalPCA/.test(html));
ok('Slice 2: visibility wires #thZPanel to hasLocalPCA',
   /showHide\('thZPanel',\s*hasLocalPCA/.test(html));
ok('Slice 3: _drawThCusumHero still defined',
   /function _drawThCusumHero\(\)/.test(html));
ok('Slice 5: _drawThLinesPanel still defined',
   /function _drawThLinesPanel\(\)/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
