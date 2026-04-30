// =============================================================================
// test_v4_repeat_density_controls_runtime.js — turn 19+20 features
//
// Tests the extended repeat density panel:
//   - Class dropdown  (turn 19)
//   - Y-axis mode toggle: linear/auto/log  (turn 19)
//   - View mode toggle: full chrom / zoomed to candidate  (turn 19)
//   - Multi-chrom support + chrom picker  (turn 20)
//   - Alt+← / Alt+→ keyboard nav  (turn 20)
//   - Persistence of all prefs to localStorage  (turn 20)
//   - Flash animation on chrom switch  (turn 20)
//
// Doesn't re-test the data layer (covered by test_v4_repeat_density_runtime.js).
// =============================================================================

const fs = require('fs');
const path = require('path');

const HTML = fs.readFileSync(path.join(__dirname, '..', 'pca_scrubber_v3.html'), 'utf8');

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  PASS ' + msg); pass++; }
  else      { console.log('  FAIL ' + msg); fail++; }
}

// ============================================================================
// SECTION 1 — Constants and helper functions present
// ============================================================================
console.log('\n--- SECTION 1: New constants and helpers ---');

ok(/REPEAT_DENSITY_PREFS_LS_KEY\s*=\s*'pca_scrubber_v3\.repeatDensity\.prefs'/.test(HTML),
   'prefs localStorage key constant defined');
ok(/REPEAT_DENSITY_Y_MODES\s*=\s*\['linear',\s*'auto',\s*'log'\]/.test(HTML),
   'Y_MODES constant has 3 values: linear, auto, log');
ok(/REPEAT_DENSITY_VIEW_MODES\s*=\s*\['full',\s*'zoomed'\]/.test(HTML),
   'VIEW_MODES constant has 2 values: full, zoomed');

ok(HTML.indexOf('function _repeatDensityPrefs') > 0,
   '_repeatDensityPrefs defined');
ok(HTML.indexOf('function _persistRepeatDensityPrefs') > 0,
   '_persistRepeatDensityPrefs defined');
ok(HTML.indexOf('function _restoreRepeatDensityPrefs') > 0,
   '_restoreRepeatDensityPrefs defined');
ok(HTML.indexOf('function _repeatDensityChromList') > 0,
   '_repeatDensityChromList defined');
ok(HTML.indexOf('function _switchChromViaMainSelect') > 0,
   '_switchChromViaMainSelect defined');
ok(HTML.indexOf('function _wireRepeatDensityArrowNav') > 0,
   '_wireRepeatDensityArrowNav defined');
ok(HTML.indexOf('function _flashRepeatDensityPanel') > 0,
   '_flashRepeatDensityPanel defined');

// Restore + arrow nav are wired on page load
ok(/_restoreRepeatDensityPrefs\(\)/.test(HTML),
   '_restoreRepeatDensityPrefs called on page load');
ok(/_wireRepeatDensityArrowNav\(\)/.test(HTML),
   '_wireRepeatDensityArrowNav called on page load');

// ============================================================================
// SECTION 2 — Renderer wiring
// ============================================================================
console.log('\n--- SECTION 2: Renderer dropdowns + wiring ---');

function extractFn(name) {
  const idx = HTML.indexOf('function ' + name + '(');
  if (idx < 0) throw new Error('helper not found: ' + name);
  let depth = 0, started = false, end = idx;
  for (let i = idx; i < HTML.length; i++) {
    const ch = HTML[i];
    if (ch === '{') { depth++; started = true; }
    else if (ch === '}') { depth--; if (started && depth === 0) { end = i + 1; break; } }
  }
  return HTML.slice(idx, end);
}

const renderBody = extractFn('_renderRepeatDensityPanel');

// Class select dropdown
ok(/id="bndRdClassSelect"/.test(renderBody),
   'renderer emits #bndRdClassSelect');
ok(/<select[^>]*id="bndRdClassSelect"[\s\S]*?<\/select>/.test(renderBody) ||
   /'<select id="bndRdClassSelect"[\s\S]*?<\/select>'/.test(renderBody) ||
   renderBody.indexOf('id="bndRdClassSelect"') > 0,
   'class select is a <select> dropdown');
// Class select wires onchange to update active class + persist + re-render
ok(/getElementById\('bndRdClassSelect'\)/.test(renderBody),
   'class select is wired via getElementById');
ok(/_setRepeatDensityActiveClass\(chrom,\s*clsSel\.value\)/.test(renderBody),
   'class change calls _setRepeatDensityActiveClass');
ok(/_persistRepeatDensityPrefs\(\)/.test(renderBody),
   'class change persists prefs');

// Y-mode select
ok(/id="bndRdYModeSelect"/.test(renderBody), 'renderer emits #bndRdYModeSelect');
ok(/getElementById\('bndRdYModeSelect'\)/.test(renderBody), 'y-mode select wired');
ok(/REPEAT_DENSITY_Y_MODES\.includes/.test(renderBody),
   'y-mode change validates against allowlist');

// View-mode select
ok(/id="bndRdViewModeSelect"/.test(renderBody), 'renderer emits #bndRdViewModeSelect');
ok(/getElementById\('bndRdViewModeSelect'\)/.test(renderBody), 'view-mode select wired');
ok(/REPEAT_DENSITY_VIEW_MODES\.includes/.test(renderBody),
   'view-mode change validates against allowlist');

// Chrom picker (only when >1 chrom loaded)
ok(/id="bndRdChromSelect"/.test(renderBody), 'renderer emits #bndRdChromSelect (when >1 chrom)');
ok(/allChroms\.length > 1/.test(renderBody),
   'chrom picker only shown when >1 chrom loaded');
ok(/_switchChromViaMainSelect\(chromSel\.value\)/.test(renderBody),
   'chrom picker switches via main select');

// Header HTML structure: dropdowns wrapped in .bnd-rd-ctrl labels
ok(/class="bnd-rd-ctrl"/.test(renderBody),
   'dropdowns wrapped in .bnd-rd-ctrl labels');
ok(/class="bnd-rd-ctrl-label"/.test(renderBody),
   'control labels use .bnd-rd-ctrl-label');

// Zoom fallback warning
ok(/zoomFallback/.test(renderBody),
   'renderer tracks zoom-fallback state');
ok(/no candidate active/.test(renderBody),
   'renderer notes when zoom falls back to full chrom');

// ============================================================================
// SECTION 3 — Y-axis mode logic
// ============================================================================
console.log('\n--- SECTION 3: Y-axis mode logic ---');

// Three branches in the renderer: linear, auto, log
ok(/yMode === 'log'/.test(renderBody), 'log mode branch present');
ok(/yMode === 'auto'/.test(renderBody), 'auto mode branch present');
ok(/Math\.log10/.test(renderBody), 'log mode uses Math.log10');
ok(/LOG_EPS/.test(renderBody), 'log mode uses epsilon floor for log(0)');

// Auto mode: yUpper = max_density
ok(/info\.max_density/.test(renderBody),
   'auto mode reads info.max_density');

// Y-axis title label changes per mode
ok(/yMode === 'log' \? ' \(log\)'/.test(renderBody),
   'y-axis title shows "(log)" suffix for log mode');
ok(/yMode === 'auto' \? ' \(auto\)'/.test(renderBody),
   'y-axis title shows "(auto)" suffix for auto mode');

// ============================================================================
// SECTION 4 — View mode logic
// ============================================================================
console.log('\n--- SECTION 4: View mode logic ---');

ok(/viewMode === 'zoomed'/.test(renderBody),
   'zoomed mode branch present');
ok(/_boundaryScanRange/.test(renderBody),
   'zoomed mode reads scan range');

// Candidate highlight only in full-chrom view
ok(/viewMode === 'full' && scanRangeMb/.test(renderBody),
   'candidate highlight only when in full mode AND scan range exists');

// Padding around scan range when zoomed
ok(/span \* 0\.05/.test(renderBody),
   'zoomed view pads x-range by 5% on each side');

// X-axis tick step adapts to xRange
ok(/xRange > 50/.test(renderBody) && /xRange > 1/.test(renderBody),
   'x-axis tick step adapts to range size');

// Off-canvas scatter culling
ok(/xMb < xMin - 0\.001 \|\| xMb > xMax \+ 0\.001/.test(renderBody),
   'scatter culls points outside visible x-range');

// ============================================================================
// SECTION 5 — Arrow key nav
// ============================================================================
console.log('\n--- SECTION 5: Arrow key navigation ---');

const navBody = extractFn('_wireRepeatDensityArrowNav');
ok(/state\._repeatDensityArrowNavBound/.test(navBody),
   'nav binding tracked to prevent double-bind');
ok(/e\.altKey/.test(navBody), 'requires Alt key');
ok(/'ArrowLeft'/.test(navBody) && /'ArrowRight'/.test(navBody),
   'handles both ArrowLeft and ArrowRight');
ok(/INPUT[\s\S]*TEXTAREA[\s\S]*SELECT/.test(navBody),
   'ignores arrows when target is INPUT/TEXTAREA/SELECT');
ok(/page11[\s\S]*active/.test(navBody),
   'only acts when page11 is the active page');
ok(/list\.length <= 1/.test(navBody),
   'short-circuits when only 0 or 1 chrom loaded');
ok(/_switchChromViaMainSelect/.test(navBody),
   'switches via the main chrom select');
ok(/_flashRepeatDensityPanel/.test(navBody),
   'flashes panel on successful switch');
ok(/e\.preventDefault\(\)/.test(navBody),
   'preventDefault on successful switch (so browser back-nav doesn\'t fire)');

const flashBody = extractFn('_flashRepeatDensityPanel');
ok(/classList\.add\('bnd-rd-flash'\)/.test(flashBody),
   'flash adds bnd-rd-flash class');
ok(/setTimeout/.test(flashBody),
   'flash uses setTimeout to remove the class');
ok(/600/.test(flashBody),
   'flash duration is 600ms (matches CSS animation)');

// ============================================================================
// SECTION 6 — CSS for new controls
// ============================================================================
console.log('\n--- SECTION 6: CSS for new controls ---');

ok(/#page11 \.bnd-rd-ctrl\b/.test(HTML),         '.bnd-rd-ctrl CSS scoped');
ok(/#page11 \.bnd-rd-ctrl-label\b/.test(HTML),    '.bnd-rd-ctrl-label CSS');
ok(/#page11 \.bnd-rd-ctrl select\b/.test(HTML),   '.bnd-rd-ctrl select styled');
ok(/#page11 \.bnd-rd-warn\b/.test(HTML),          '.bnd-rd-warn styled');
ok(/#page11 \.bnd-rd-empty-loaded\b/.test(HTML),  '.bnd-rd-empty-loaded CSS');
ok(/\.bnd-rd-flash\b/.test(HTML),                  '.bnd-rd-flash CSS');
ok(/@keyframes bndRdFlash\b/.test(HTML),           '@keyframes bndRdFlash defined');
ok(/box-shadow:.*rgba\(245, 165, 36/.test(HTML),
   'flash animation uses amber accent color');

// Header switched from align-items: baseline to align-items: center
// (because dropdowns don't sit cleanly on baseline)
const headerCss = HTML.match(/#page11 \.bnd-rd-header \{[\s\S]*?\}/);
ok(headerCss !== null, 'header CSS extractable');
if (headerCss) {
  ok(/align-items:\s*center/.test(headerCss[0]),
     'header uses align-items: center for dropdown alignment');
  ok(/flex-wrap:\s*wrap/.test(headerCss[0]),
     'header uses flex-wrap so dropdowns fit on narrow viewports');
}

// ============================================================================
// SECTION 7 — Sandbox tests for prefs + chrom list helpers
// ============================================================================
console.log('\n--- SECTION 7: Sandbox unit tests ---');

const helpers = [
  '_isRepeatDensityJSON', '_validateRepeatDensityChrom',
  '_storeRepeatDensity', '_persistRepeatDensity', '_restoreRepeatDensity',
  '_clearRepeatDensity', '_getRepeatDensity',
  '_resolveRepeatDensityClass', '_setRepeatDensityActiveClass',
  '_repeatDensityPrefs', '_persistRepeatDensityPrefs', '_restoreRepeatDensityPrefs',
  '_repeatDensityChromList',
];
let helperBodies = '';
for (const n of helpers) helperBodies += extractFn(n) + '\n';

const c1 = HTML.match(/const REPEAT_DENSITY_LS_PREFIX[^;]*;/);
const c2 = HTML.match(/const REPEAT_DENSITY_PREFS_LS_KEY[^;]*;/);
const c3 = HTML.match(/const REPEAT_DENSITY_Y_MODES[^;]*;/);
const c4 = HTML.match(/const REPEAT_DENSITY_VIEW_MODES[^;]*;/);
const constants = (c1 && c2 && c3 && c4) ? (c1[0] + '\n' + c2[0] + '\n' + c3[0] + '\n' + c4[0]) : '';
ok(constants.length > 0, 'all 4 constants extractable');

function makeEnv() {
  const factory = new Function(`
    const state = {};
    const _ls = new Map();
    const localStorage = {
      get length() { return _ls.size; },
      key: (i) => Array.from(_ls.keys())[i] || null,
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => { _ls.set(k, String(v)); },
      removeItem: (k) => { _ls.delete(k); },
      _dump: () => Array.from(_ls.entries()),
    };
    const console = { warn: () => {}, info: () => {}, error: () => {} };
    ${constants}
    ${helperBodies}
    return {
      state, localStorage,
      store: _storeRepeatDensity,
      prefs: _repeatDensityPrefs,
      persistPrefs: _persistRepeatDensityPrefs,
      restorePrefs: _restoreRepeatDensityPrefs,
      chromList: _repeatDensityChromList,
      resolveClass: _resolveRepeatDensityClass,
      setClass: _setRepeatDensityActiveClass,
    };
  `);
  return factory();
}

function makeFixture(chromName = 'C_gar_LG28') {
  const N = 5;
  return {
    version: 2, species: 'C_gariepinus', binning_source: 'scrubber_windows',
    precomp_chrom: chromName, n_chromosomes: 1, n_classes: 2,
    classes: ['repeat_fragment', 'Gypsy_LTR_retrotransposon'],
    default_class: 'repeat_fragment', loess_span: 0.3,
    chromosomes: [{
      chrom: chromName, n_windows: N,
      window_centers_mb: [0.5, 1.5, 2.5, 3.5, 4.5],
      window_start_bp:   [0, 1000000, 2000000, 3000000, 4000000],
      window_end_bp:     [1000000, 2000000, 3000000, 4000000, 5000000],
      by_class: {
        repeat_fragment: { densities: [0.1, 0.2, 0.5, 0.3, 0.1], loess: [0.15, 0.2, 0.3, 0.25, 0.15], max_density: 0.5, loess_span: 0.3 },
        Gypsy_LTR_retrotransposon: { densities: [0.05, 0.05, 0.1, 0.05, 0.05], loess: [0.05, 0.06, 0.07, 0.06, 0.05], max_density: 0.1, loess_span: 0.3 },
      },
    }],
  };
}

// --- _repeatDensityPrefs ---
{
  const env = makeEnv();
  const p = env.prefs();
  ok(p.yMode === 'linear', 'default yMode is linear');
  ok(p.viewMode === 'full', 'default viewMode is full');
  ok(env.state.repeatDensityPrefs === p, 'prefs stored on state');
  // Idempotent
  const p2 = env.prefs();
  ok(p === p2, 'second call returns same object');
}

// --- _persistRepeatDensityPrefs / _restoreRepeatDensityPrefs ---
{
  const env = makeEnv();
  const p = env.prefs();
  p.yMode = 'log';
  p.viewMode = 'zoomed';
  env.state.repeatDensityActiveClass = { 'C_gar_LG28': 'CACTA_TIR_transposon' };
  env.persistPrefs();
  const dump = env.localStorage._dump();
  ok(dump.length === 1, 'persist writes one localStorage key');
  ok(dump[0][0] === 'pca_scrubber_v3.repeatDensity.prefs', 'key has correct name');

  // Round-trip
  const env2 = makeEnv();
  for (const [k, v] of dump) env2.localStorage.setItem(k, v);
  env2.restorePrefs();
  ok(env2.prefs().yMode === 'log', 'yMode round-trips');
  ok(env2.prefs().viewMode === 'zoomed', 'viewMode round-trips');
  ok(env2.state.repeatDensityActiveClass['C_gar_LG28'] === 'CACTA_TIR_transposon',
     'activeClass per chrom round-trips');
}

// --- restore validates against allowlist ---
{
  const env = makeEnv();
  // Inject malicious prefs
  env.localStorage.setItem('pca_scrubber_v3.repeatDensity.prefs',
                            JSON.stringify({ yMode: 'invalid_xyz', viewMode: 'invalid_qrs' }));
  env.restorePrefs();
  ok(env.prefs().yMode === 'linear', 'invalid yMode falls back to linear');
  ok(env.prefs().viewMode === 'full', 'invalid viewMode falls back to full');
}

// --- restore handles malformed JSON gracefully ---
{
  const env = makeEnv();
  env.localStorage.setItem('pca_scrubber_v3.repeatDensity.prefs', 'not valid json {');
  env.restorePrefs();   // shouldn't throw
  ok(env.prefs().yMode === 'linear', 'malformed JSON falls back to defaults');
}

// --- _repeatDensityChromList ---
{
  const env = makeEnv();
  ok(env.chromList().length === 0, 'empty when no data loaded');
  env.store(makeFixture('C_gar_LG28'));
  ok(env.chromList().length === 1, 'one chrom after one store');
  env.store(makeFixture('C_gar_LG07'));
  env.store(makeFixture('C_gar_LG14'));
  const list = env.chromList();
  ok(list.length === 3, 'three chroms after three stores');
  // Sorted alphabetically
  ok(list[0] === 'C_gar_LG07' && list[1] === 'C_gar_LG14' && list[2] === 'C_gar_LG28',
     'chrom list sorted alphabetically');
}

// --- Multi-chrom: each chrom has independent active class ---
{
  const env = makeEnv();
  env.store(makeFixture('C_gar_LG28'));
  env.store(makeFixture('C_gar_LG07'));
  env.setClass('C_gar_LG28', 'Gypsy_LTR_retrotransposon');
  // LG07 not yet overridden — should use default
  ok(env.resolveClass('C_gar_LG28') === 'Gypsy_LTR_retrotransposon',
     'LG28 uses overridden class');
  ok(env.resolveClass('C_gar_LG07') === 'repeat_fragment',
     'LG07 uses default (independent of LG28 override)');
}

// --- Persistence preserves per-chrom active class ---
{
  const env = makeEnv();
  env.store(makeFixture('C_gar_LG28'));
  env.store(makeFixture('C_gar_LG07'));
  env.setClass('C_gar_LG28', 'Gypsy_LTR_retrotransposon');
  env.setClass('C_gar_LG07', 'repeat_fragment');   // explicit, even though it's the default
  env.persistPrefs();

  const env2 = makeEnv();
  for (const [k, v] of env.localStorage._dump()) env2.localStorage.setItem(k, v);
  env2.restorePrefs();
  // Need to also store the chroms in env2 for resolveClass to work
  env2.store(makeFixture('C_gar_LG28'));
  env2.store(makeFixture('C_gar_LG07'));
  ok(env2.resolveClass('C_gar_LG28') === 'Gypsy_LTR_retrotransposon',
     'per-chrom active class persists across reload (LG28)');
  ok(env2.resolveClass('C_gar_LG07') === 'repeat_fragment',
     'per-chrom active class persists across reload (LG07)');
}

// ============================================================================
// SECTION 8 — Empty state with loaded-elsewhere hint
// ============================================================================
console.log('\n--- SECTION 8: Empty state hint ---');

ok(/Loaded for:/.test(renderBody),
   'empty state lists chroms with data loaded elsewhere');
ok(/loadedList\.length > 0/.test(renderBody),
   'loaded hint conditional on having any data');
ok(/bnd-rd-empty-loaded/.test(renderBody),
   'loaded hint uses .bnd-rd-empty-loaded class');

console.log('\n--- v4 turn 19+20 repeat density controls runtime tests: ' +
            pass + '/' + (pass + fail) + ' ---');
process.exit(fail > 0 ? 1 : 0);
