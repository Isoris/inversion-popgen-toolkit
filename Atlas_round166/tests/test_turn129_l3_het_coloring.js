// =============================================================================
// turn 129 test — L3 het-rate dot coloring (Slice 1 of
// SPEC_l3_het_dosage_coloring.md)
//
// Quentin's request (turn 128d):
//   "Now that we have dosage, can we get Het (dosage) coloring implemented or
//   all panels of l3 contingency? if activated? Add a button for it. but
//   always active for tracked sample panel unless toggled off?"
//
// Slice 1 ships:
//   - state.l3HetColoring boolean slot (default false)
//   - L3 toolbar [ ] het checkbox + change handler
//   - localStorage persistence (pca_scrubber_v3.l3HetColoring)
//   - _computeHetRateForL2(l2idx) → Float32Array[n_samples] of rates ∈ [0,1] | NaN
//   - _hetRateColor(rate) → CSS color (cold / neutral / warm) | dim for NaN
//   - drawMiniPCA hook overriding non-tracked dot fill when toggle on
//   - Color-bar legend on focal mini-PCA (focal pane only)
//   - Cache invalidation when dosage_chunks layer (re)loads
//   - Toggle disabled when dosage_chunks layer is missing
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

// =============================================================================
// Source-level checks
// =============================================================================
console.log('\n=== Source-level checks ===');

ok('state slot l3HetColoring exists with default false',
   /l3HetColoring:\s*false/.test(html));

ok('L3 toolbar checkbox markup present (id=l3HetToggle)',
   /<input[^>]*id="l3HetToggle"[^>]*type="checkbox"|<input[^>]*type="checkbox"[^>]*id="l3HetToggle"/.test(html));

ok('L3 toolbar het label container present (id=l3HetToggleLabel)',
   /id="l3HetToggleLabel"/.test(html));

ok('_hetRateColor function defined',
   /function _hetRateColor\(rate\)/.test(html));

ok('_computeHetRateForL2 function defined',
   /function _computeHetRateForL2\(l2idx\)/.test(html));

ok('_invalidateHetRateCache function defined',
   /function _invalidateHetRateCache\(\)/.test(html));

ok('_getHetRateCache function defined',
   /function _getHetRateCache\(\)/.test(html));

ok('_hetRateColor exposed on window',
   /window\._hetRateColor\s*=\s*_hetRateColor/.test(html));

ok('_computeHetRateForL2 exposed on window',
   /window\._computeHetRateForL2\s*=\s*_computeHetRateForL2/.test(html));

ok('_invalidateHetRateCache exposed on window',
   /window\._invalidateHetRateCache\s*=\s*_invalidateHetRateCache/.test(html));

ok('drawMiniPCA references state.l3HetColoring (the het branch is wired)',
   /function drawMiniPCA\([\s\S]{0,30000}state\.l3HetColoring/.test(html));

ok('drawMiniPCA references _hetRateColor (color override active)',
   /function drawMiniPCA\([\s\S]{0,30000}_hetRateColor/.test(html));

ok('drawMiniPCA calls _computeHetRateForL2 (per-L2 rate fetch)',
   /function drawMiniPCA\([\s\S]{0,30000}_computeHetRateForL2\(l2idx\)/.test(html));

ok('color-bar legend gated on focal pane (paneOffset === 0)',
   /paneOffset === 0 && state\.l3HetColoring/.test(html));

ok('localStorage key uses pca_scrubber_v3.l3HetColoring',
   /pca_scrubber_v3\.l3HetColoring/.test(html));

ok('change handler restored from localStorage',
   /localStorage\.getItem\('pca_scrubber_v3\.l3HetColoring'\)/.test(html));

ok('change handler writes to localStorage on toggle',
   /localStorage\.setItem\('pca_scrubber_v3\.l3HetColoring'/.test(html));

ok('dosage_chunks layer-load invalidates het cache',
   /added\.indexOf\('dosage_chunks'\) >= 0[\s\S]{0,400}_invalidateHetRateCache/.test(html));

ok('toggle availability refreshes after layer load',
   /added\.indexOf\('dosage_chunks'\) >= 0[\s\S]{0,500}_refreshL3HetToggleAvailability/.test(html));

ok('_refreshL3HetToggleAvailability function defined',
   /function _refreshL3HetToggleAvailability\(\)/.test(html));

ok('toggle force-off when dosage layer disappears mid-session',
   /!hasDosage && state\.l3HetColoring[\s\S]{0,200}state\.l3HetColoring = false/.test(html));

ok('withAlpha extended to handle rgb() format from _hetRateColor',
   /col\.indexOf\('rgb\('\) === 0/.test(html));

// =============================================================================
// Behavioural tests (sandboxed)
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        i++;
      }
    } else if (ch === '/' && src[i+1] === '/') {
      while (i < src.length && src[i] !== '\n') i++;
    } else if (ch === '/' && src[i+1] === '*') {
      i += 2;
      while (i < src.length - 1 && !(src[i] === '*' && src[i+1] === '/')) i++;
      i++;
    }
    i++;
  }
  return src.substring(start, i);
}

const fnHetColor = pullFunction(html, '_hetRateColor');
const fnComputeHet = pullFunction(html, '_computeHetRateForL2');
// turn 151: _computeHetRateForL2 was refactored into a thin wrapper that
// delegates to _computeHetRateForRange. Tests that exercise the L2 helper
// must now also load the range core into the sandbox.
const fnComputeHetRange = pullFunction(html, '_computeHetRateForRange');
const fnGetCache = pullFunction(html, '_getHetRateCache');
const fnInvalidateCache = pullFunction(html, '_invalidateHetRateCache');
const fnWithAlpha = pullFunction(html, 'withAlpha');

ok('_hetRateColor extractable',         !!fnHetColor);
ok('_computeHetRateForL2 extractable',  !!fnComputeHet);
ok('_computeHetRateForRange extractable (turn 151 split)', !!fnComputeHetRange);
ok('_getHetRateCache extractable',      !!fnGetCache);
ok('_invalidateHetRateCache extractable', !!fnInvalidateCache);
ok('withAlpha extractable',             !!fnWithAlpha);

// Pull the _HET_RAMP constant block too — the color helper needs it.
const rampMatch = html.match(/const _HET_RAMP = \{[\s\S]*?\};/);
ok('_HET_RAMP constant extractable', !!rampMatch);

function makeColorSandbox() {
  const sandbox = { console };
  vm.createContext(sandbox);
  vm.runInContext(rampMatch[0], sandbox);
  vm.runInContext(fnHetColor, sandbox);
  return sandbox;
}

// --- Test A: _hetRateColor returns dim for NaN / null
console.log('\nTest A: _hetRateColor handles NaN / null');
{
  const sandbox = makeColorSandbox();
  vm.runInContext('var c1 = _hetRateColor(NaN);',  sandbox);
  vm.runInContext('var c2 = _hetRateColor(null);', sandbox);
  vm.runInContext('var c3 = _hetRateColor(undefined);', sandbox);
  ok('NaN → var(--ink-dimmer)',  sandbox.c1 === 'var(--ink-dimmer)');
  ok('null → var(--ink-dimmer)', sandbox.c2 === 'var(--ink-dimmer)');
  ok('undefined → var(--ink-dimmer)', sandbox.c3 === 'var(--ink-dimmer)');
}

// --- Test B: _hetRateColor endpoints + midpoint
console.log('\nTest B: _hetRateColor ramp endpoints');
{
  const sandbox = makeColorSandbox();
  vm.runInContext('var c0 = _hetRateColor(0);',   sandbox);
  vm.runInContext('var c5 = _hetRateColor(0.5);', sandbox);
  vm.runInContext('var c1 = _hetRateColor(1);',   sandbox);
  ok('rate=0 is the cold-blue (#2166AC) endpoint',
     sandbox.c0 === 'rgb(33,102,172)',
     'got: ' + sandbox.c0);
  ok('rate=0.5 is the neutral (#F7F7F7) midpoint',
     sandbox.c5 === 'rgb(247,247,247)',
     'got: ' + sandbox.c5);
  ok('rate=1 is the warm-red (#B2182B) endpoint',
     sandbox.c1 === 'rgb(178,24,43)',
     'got: ' + sandbox.c1);
}

// --- Test C: _hetRateColor interpolates smoothly between endpoints
console.log('\nTest C: _hetRateColor smooth interpolation');
{
  const sandbox = makeColorSandbox();
  vm.runInContext('var c25 = _hetRateColor(0.25);', sandbox);
  vm.runInContext('var c75 = _hetRateColor(0.75);', sandbox);
  // 0.25 should be midway between cold and neutral
  // R: 33 → (33+247)/2 = 140 ; G: 102 → (102+247)/2 = 174.5 ; B: 172 → (172+247)/2 = 209.5
  ok('rate=0.25 is between cold and neutral',
     sandbox.c25 === 'rgb(140,175,210)' || sandbox.c25 === 'rgb(140,174,210)' ||
     sandbox.c25 === 'rgb(140,175,209)' || sandbox.c25 === 'rgb(140,174,209)',
     'got: ' + sandbox.c25);
  // 0.75 should be midway between neutral and warm
  ok('rate=0.75 is between neutral and warm',
     /^rgb\(\d+,\d+,\d+\)$/.test(sandbox.c75) &&
     // R: 247 → (247+178)/2 = 212.5 → 213 (or 212 by rounding)
     (sandbox.c75.startsWith('rgb(212,') || sandbox.c75.startsWith('rgb(213,')),
     'got: ' + sandbox.c75);
}

// --- Test D: clamping out-of-range
console.log('\nTest D: clamping out-of-range');
{
  const sandbox = makeColorSandbox();
  vm.runInContext('var c_neg = _hetRateColor(-0.5);', sandbox);
  vm.runInContext('var c_big = _hetRateColor(2);', sandbox);
  ok('rate < 0 clamps to cold-blue',
     sandbox.c_neg === 'rgb(33,102,172)');
  ok('rate > 1 clamps to warm-red',
     sandbox.c_big === 'rgb(178,24,43)');
}

// --- Test E: _computeHetRateForL2 on a synthetic chunk
console.log('\nTest E: _computeHetRateForL2 on a synthetic chunk');
{
  // Build a 4-window L2 (bp 1000–1003) with 3 samples.
  // Sample S0: dosage [0, 1, 1, 0] → het=2, valid=4 → 0.50
  // Sample S1: dosage [0, 0, 0, 0] → het=0, valid=4 → 0.00
  // Sample S2: dosage [-1, -1, -1, -1] → all NA → NaN
  const chunk = {
    samples: ['S0', 'S1', 'S2'],
    markers: [
      { pos_bp: 1000 },
      { pos_bp: 1001 },
      { pos_bp: 1002 },
      { pos_bp: 1003 },
    ],
    dosage: [
      [0, 0, -1],   // marker 0 across 3 samples
      [1, 0, -1],   // marker 1
      [1, 0, -1],   // marker 2
      [0, 0, -1],   // marker 3
    ],
  };
  const stateMock = {
    data: {
      n_samples: 3,
      samples: [{ id: 'S0' }, { id: 'S1' }, { id: 'S2' }],
      l2_envelopes: [{ start_bp: 1000, end_bp: 1003 }],
    },
  };
  const popgenMock = {
    getCachedChunk(start_bp, end_bp) {
      return chunk;
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: popgenMock },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnInvalidateCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var rates = _computeHetRateForL2(0);', sandbox);
  ok('returns Float32Array',
     sandbox.rates && sandbox.rates.constructor && sandbox.rates.constructor.name === 'Float32Array',
     'ctor: ' + (sandbox.rates && sandbox.rates.constructor && sandbox.rates.constructor.name));
  ok('returns array of length n_samples (3)', sandbox.rates.length === 3);
  ok('S0 het rate is 0.5 (2 het / 4 valid)',
     Math.abs(sandbox.rates[0] - 0.5) < 1e-6,
     'got: ' + sandbox.rates[0]);
  ok('S1 het rate is 0.0 (0 het / 4 valid)',
     Math.abs(sandbox.rates[1] - 0.0) < 1e-6,
     'got: ' + sandbox.rates[1]);
  ok('S2 het rate is NaN (all NA)',
     Number.isNaN(sandbox.rates[2]),
     'got: ' + sandbox.rates[2]);
}

// --- Test F: cache hit returns same Float32Array reference
console.log('\nTest F: per-L2 cache hits avoid recomputation');
{
  const chunk = {
    samples: ['S0'],
    markers: [{ pos_bp: 1000 }, { pos_bp: 1001 }],
    dosage: [[1], [0]],
  };
  const stateMock = {
    data: {
      n_samples: 1,
      samples: [{ id: 'S0' }],
      l2_envelopes: [{ start_bp: 1000, end_bp: 1001 }],
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: { getCachedChunk: () => chunk } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnInvalidateCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var r1 = _computeHetRateForL2(0); var r2 = _computeHetRateForL2(0);', sandbox);
  ok('first call returns Float32Array',
     sandbox.r1 && sandbox.r1.constructor && sandbox.r1.constructor.name === 'Float32Array');
  ok('second call returns SAME reference (cached)', sandbox.r1 === sandbox.r2);

  // After invalidation, fresh array is returned.
  vm.runInContext('_invalidateHetRateCache(); var r3 = _computeHetRateForL2(0);', sandbox);
  ok('after invalidation, new array is returned', sandbox.r3 !== sandbox.r1);
  ok('after invalidation, contents still correct',
     Math.abs(sandbox.r3[0] - 0.5) < 1e-6);
}

// --- Test G: no covering chunk → all-NaN
console.log('\nTest G: no covering chunk → all-NaN');
{
  const stateMock = {
    data: {
      n_samples: 3,
      samples: [{ id: 'S0' }, { id: 'S1' }, { id: 'S2' }],
      l2_envelopes: [{ start_bp: 1000, end_bp: 2000 }],
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: { getCachedChunk: () => null } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var rates = _computeHetRateForL2(0);', sandbox);
  ok('still returns Float32Array of length n_samples',
     sandbox.rates && sandbox.rates.constructor && sandbox.rates.constructor.name === 'Float32Array' && sandbox.rates.length === 3);
  ok('all entries are NaN',
     Number.isNaN(sandbox.rates[0]) &&
     Number.isNaN(sandbox.rates[1]) &&
     Number.isNaN(sandbox.rates[2]));
}

// --- Test H: missing L2 envelope → all-NaN, no crash
console.log('\nTest H: missing L2 envelope → all-NaN, no crash');
{
  const stateMock = {
    data: {
      n_samples: 2,
      samples: [{ id: 'S0' }, { id: 'S1' }],
      l2_envelopes: [],
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: { getCachedChunk: () => null } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var rates = _computeHetRateForL2(99);', sandbox);
  ok('returns Float32Array even for unknown l2idx',
     sandbox.rates && sandbox.rates.constructor && sandbox.rates.constructor.name === 'Float32Array');
  ok('all entries are NaN', Number.isNaN(sandbox.rates[0]) && Number.isNaN(sandbox.rates[1]));
}

// --- Test I: chunk samples not present in cohort → those entries stay NaN
console.log('\nTest I: chunk-only samples are projected, missing ones stay NaN');
{
  const chunk = {
    samples: ['S0', 'X_NOT_IN_COHORT'],   // S0 in cohort, X not
    markers: [{ pos_bp: 1000 }],
    dosage: [[1, 1]],   // both het
  };
  const stateMock = {
    data: {
      n_samples: 2,
      samples: [{ id: 'S0' }, { id: 'S1' }],   // S1 isn't in chunk
      l2_envelopes: [{ start_bp: 1000, end_bp: 1000 }],
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: { getCachedChunk: () => chunk } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var rates = _computeHetRateForL2(0);', sandbox);
  ok('S0 (in chunk + cohort) gets the het rate',
     Math.abs(sandbox.rates[0] - 1.0) < 1e-6,
     'got: ' + sandbox.rates[0]);
  ok('S1 (in cohort, not in chunk) stays NaN',
     Number.isNaN(sandbox.rates[1]),
     'got: ' + sandbox.rates[1]);
}

// --- Test J: withAlpha now handles rgb() format
console.log('\nTest J: withAlpha extension for rgb() inputs');
{
  const sandbox = { console };
  vm.createContext(sandbox);
  vm.runInContext(fnWithAlpha, sandbox);
  vm.runInContext('var a = withAlpha("rgb(100,200,50)", 0.5);', sandbox);
  vm.runInContext('var b = withAlpha("#ABCDEF", 0.3);', sandbox);
  vm.runInContext('var c = withAlpha("rgba(1,2,3,1)", 0.4);', sandbox);
  vm.runInContext('var d = withAlpha("var(--ink)", 0.4);', sandbox);
  ok('rgb() input → rgba() with alpha',
     sandbox.a === 'rgba(100,200,50,0.5)',
     'got: ' + sandbox.a);
  ok('# input still produces rgba() (regression check)',
     sandbox.b === 'rgba(171,205,239,0.3)',
     'got: ' + sandbox.b);
  ok('rgba() input passes through (assumed already alpha-set)',
     sandbox.c === 'rgba(1,2,3,1)');
  ok('CSS var input passes through',
     sandbox.d === 'var(--ink)');
}

// --- Test K: marker bp filtering excludes markers outside L2
console.log('\nTest K: marker bp filtering scopes het rate to L2');
{
  const chunk = {
    samples: ['S0'],
    markers: [
      { pos_bp: 999 },    // BEFORE L2
      { pos_bp: 1000 },   // in L2
      { pos_bp: 2000 },   // in L2
      { pos_bp: 2001 },   // AFTER L2
    ],
    dosage: [
      [1],   // het, but outside L2 → ignored
      [1],   // het, in L2 → counted
      [0],   // homo, in L2 → counted (denominator only)
      [1],   // het, but outside L2 → ignored
    ],
  };
  const stateMock = {
    data: {
      n_samples: 1,
      samples: [{ id: 'S0' }],
      l2_envelopes: [{ start_bp: 1000, end_bp: 2000 }],
    },
  };
  const sandbox = {
    state: stateMock,
    window: { popgenDosage: { getCachedChunk: () => chunk } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnGetCache, sandbox);
  vm.runInContext(fnComputeHetRange, sandbox);
  vm.runInContext(fnComputeHet, sandbox);
  vm.runInContext('var rates = _computeHetRateForL2(0);', sandbox);
  // Only markers 1 and 2 are in-L2: 1 het + 1 homo = 0.5
  ok('only in-L2 markers contribute (rate = 0.5, not 0.75 or 1.0)',
     Math.abs(sandbox.rates[0] - 0.5) < 1e-6,
     'got: ' + sandbox.rates[0]);
}

// =============================================================================
// Final tally
// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
