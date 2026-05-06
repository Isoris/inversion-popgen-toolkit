// =============================================================================
// turn 130 test — Lineage UI surfacing (Slice 2 of
// SPEC_distant_band_concordance_fish_trajectory.md)
//
// Slice 1 shipped the compute (runLineageCompute, _concordanceMatrix,
// _hungarianChainProjection, _lineageClustering). Slice 2 surfaces the
// result in the UI:
//   - state.linesColorMode === 'lineage' added to the picker
//   - _lineageColor(si) resolver reading state.lineageResult
//   - Both _resolveSampleColorByMode and _resolveSampleScopeColor
//     dispatch to lineage mode
//   - _drawLineageStrip drawer above the per-sample-lines panel
//   - state.linesLineageStripOn slot (default true) + setter +
//     localStorage persistence
//   - applyData clears state.lineageResult on chrom switch
//   - Toggle checkbox in the lines header
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

ok('lineage mode registered in _LINES_COLOR_MODES',
   /_LINES_COLOR_MODES = \[[\s\S]{0,1500}id:\s*'lineage',\s*layer:\s*null/.test(html));

ok('lineage option present in picker dropdown',
   /<option value="lineage"[\s\S]{0,400}lineage<\/option>/.test(html));

ok('_lineageColor function defined',
   /function _lineageColor\(si\)/.test(html));

ok('_lineageColor uses golden-angle hue rotation',
   /goldenAngle = 137\.508/.test(html));

ok('_lineageColor exposed on window',
   /window\._lineageColor\s*=\s*_lineageColor/.test(html));

ok('_resolveSampleColorByMode dispatches lineage',
   /function _resolveSampleColorByMode\([\s\S]{0,800}mode === 'lineage'[\s\S]{0,80}_lineageColor\(si\)/.test(html));

ok('_resolveSampleScopeColor dispatches lineage',
   /function _resolveSampleScopeColor\([\s\S]{0,400}mode === 'lineage'[\s\S]{0,80}_lineageColor\(si\)/.test(html));

ok('_lineageColor auto-triggers compute when result missing',
   /function _lineageColor\([\s\S]{0,1500}runLineageCompute === 'function'[\s\S]{0,400}_lineageComputeScheduled/.test(html));

ok('_drawLineageStrip function defined',
   /function _drawLineageStrip\(ctx, pad, plotW, plotH, mbMin, mbMax\)/.test(html));

ok('_drawLineageStrip exposed on window',
   /window\._drawLineageStrip\s*=\s*_drawLineageStrip/.test(html));

ok('_drawLineageStrip respects linesLineageStripOn === false (off path)',
   /_state\.linesLineageStripOn === false/.test(html));

ok('lineage strip uses same golden-angle palette as _lineageColor',
   /_drawLineageStrip[\s\S]{0,4000}goldenAngle = 137\.508/.test(html));

ok('lineage strip wired into drawLinesPanel render path',
   /isPC1 && typeof _drawLineageStrip === 'function'[\s\S]{0,200}_drawLineageStrip\(ctx, pad, plotW, plotH, mbMin, mbMax\)/.test(html));

ok('state slot linesLineageStripOn defined (default true)',
   /linesLineageStripOn:\s*true/.test(html));

ok('setLinesLineageStripOn function defined',
   /function setLinesLineageStripOn\(b\)/.test(html));

ok('setLinesLineageStripOn writes to localStorage',
   /function setLinesLineageStripOn\([\s\S]{0,400}localStorage\.setItem\(_LINES_LINEAGE_STRIP_KEY/.test(html));

ok('localStorage key for lineage strip is namespaced',
   /_LINES_LINEAGE_STRIP_KEY = 'inversion_atlas\.linesLineageStripOn'/.test(html));

ok('lineage strip toggle wireup reads localStorage on init',
   /linesLineageStripToggle[\s\S]{0,600}localStorage\.getItem\(_LINES_LINEAGE_STRIP_KEY\)/.test(html));

ok('lineage strip toggle defaults checked in markup (default ON contract)',
   /<input[^>]*id="linesLineageStripToggle"[^>]*checked|<input[^>]*checked[^>]*id="linesLineageStripToggle"/.test(html));

ok('applyData invalidates lineage cache on chrom switch',
   /state\.data = data;[\s\S]{0,400}invalidateLineageCache/.test(html));

ok('applyData resets _lineageComputeScheduled on chrom switch',
   /state\._lineageComputeScheduled = false/.test(html));

ok('lineage strip falls back gracefully when state.lineageResult is null',
   /_drawLineageStrip[\s\S]{0,1500}!result \|\| !result\.lineage_id_per_sample/.test(html));

ok('lineage strip respects Hungarian chain breaks (separate visual treatment)',
   /_drawLineageStrip[\s\S]{0,3500}!inChain\.has\(l2idx\)/.test(html));

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

const fnLineageColor = pullFunction(html, '_lineageColor');
const fnResolveByMode = pullFunction(html, '_resolveSampleColorByMode');
const fnResolveScope = pullFunction(html, '_resolveSampleScopeColor');

ok('_lineageColor extractable',                !!fnLineageColor);
ok('_resolveSampleColorByMode extractable',    !!fnResolveByMode);
ok('_resolveSampleScopeColor extractable',     !!fnResolveScope);

function makeSandbox(stateOverride) {
  const sandbox = {
    state: stateOverride || {},
    console,
    // Stub these out so resolver bodies don't choke
    familyColor: function() { return '#888888'; },
    runLineageCompute: function() { /* noop in tests */ },
    requestIdleCallback: undefined,
    setTimeout: setTimeout,
  };
  // Register on window to mirror production globalness.
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fnLineageColor,  sandbox);
  vm.runInContext(fnResolveByMode, sandbox);
  vm.runInContext(fnResolveScope,  sandbox);
  return sandbox;
}

// --- Test 1: _lineageColor returns null when state.lineageResult missing
console.log('\nTest 1: _lineageColor returns null when result missing');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },     // 3 L2s so auto-trigger is allowed
  });
  vm.runInContext('var c = _lineageColor(0);', sandbox);
  ok('returns null', sandbox.c === null);
  ok('schedules auto-trigger on state',
     sandbox.state._lineageComputeScheduled === true);
}

// --- Test 2: _lineageColor uses golden-angle palette and returns hsl()
console.log('\nTest 2: _lineageColor returns hsl() with golden-angle hue');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },
    lineageResult: {
      lineage_id_per_sample: new Int32Array([0, 1, 2, 3]),
      n_samples: 4,
      n_lineages: 4,
    },
  });
  vm.runInContext(`
    var c0 = _lineageColor(0);
    var c1 = _lineageColor(1);
    var c2 = _lineageColor(2);
    var c3 = _lineageColor(3);
  `, sandbox);
  // Expected hues: 210 + 0*137.508 = 210
  //                210 + 1*137.508 = 347.508
  //                210 + 2*137.508 = 485.016 % 360 = 125.016
  //                210 + 3*137.508 = 622.524 % 360 = 262.524
  ok('lineage 0 → hue 210', /^hsl\(210\.0,/.test(sandbox.c0), 'got: ' + sandbox.c0);
  ok('lineage 1 → hue 347.5', /^hsl\(347\.5,/.test(sandbox.c1), 'got: ' + sandbox.c1);
  ok('lineage 2 → hue 125.0', /^hsl\(125\.0,/.test(sandbox.c2), 'got: ' + sandbox.c2);
  ok('lineage 3 → hue 262.5', /^hsl\(262\.5,/.test(sandbox.c3), 'got: ' + sandbox.c3);
  ok('every color is HSL with 70% saturation',
     [sandbox.c0, sandbox.c1, sandbox.c2, sandbox.c3].every(c => /, 70%,/.test(c)));
  ok('every color is HSL with 55% lightness',
     [sandbox.c0, sandbox.c1, sandbox.c2, sandbox.c3].every(c => /, 55%\)$/.test(c)));
}

// --- Test 3: _lineageColor returns null for invalid sample index
console.log('\nTest 3: _lineageColor handles bad inputs');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },
    lineageResult: {
      lineage_id_per_sample: new Int32Array([0, 1]),
      n_samples: 2,
      n_lineages: 2,
    },
  });
  vm.runInContext(`
    var oob = _lineageColor(99);
    var neg = _lineageColor(-1);
  `, sandbox);
  ok('out-of-range si returns null',  sandbox.oob === null);
  ok('negative si returns null',       sandbox.neg === null);
}

// --- Test 4: _lineageColor returns null for invalid lineage_id (-1)
console.log('\nTest 4: _lineageColor handles invalid lineage_id (no-call)');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },
    lineageResult: {
      lineage_id_per_sample: new Int32Array([0, -1]),
      n_samples: 2,
      n_lineages: 1,
    },
  });
  vm.runInContext(`
    var c0 = _lineageColor(0);
    var c1 = _lineageColor(1);
  `, sandbox);
  ok('valid lineage 0 → color',  /^hsl\(/.test(sandbox.c0));
  ok('lineage_id -1 → null',     sandbox.c1 === null);
}

// --- Test 5: _resolveSampleColorByMode dispatches 'lineage' mode
console.log('\nTest 5: _resolveSampleColorByMode dispatches lineage');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },
    lineageResult: {
      lineage_id_per_sample: new Int32Array([0]),
      n_samples: 1,
      n_lineages: 1,
    },
  });
  vm.runInContext(`
    var c = _resolveSampleColorByMode(0, null, 'lineage');
    var k = _resolveSampleColorByMode(0, null, 'kmeans');
  `, sandbox);
  ok('mode=lineage produces hsl color', /^hsl\(/.test(sandbox.c));
  ok('mode=kmeans returns null (falls through to existing logic)',
     sandbox.k === null);
}

// --- Test 6: _resolveSampleScopeColor dispatches 'lineage' mode
console.log('\nTest 6: _resolveSampleScopeColor dispatches lineage');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}, {}] },
    lineageResult: {
      lineage_id_per_sample: new Int32Array([0]),
      n_samples: 1,
      n_lineages: 1,
    },
  });
  vm.runInContext(`
    var c = _resolveSampleScopeColor(0, 'lineage');
    var f = _resolveSampleScopeColor(0, 'family');
  `, sandbox);
  ok('mode=lineage produces hsl color', /^hsl\(/.test(sandbox.c));
  // Family color path is stubbed to '#888888' in our sandbox
  ok('mode=family delegates to familyColor',  sandbox.f === '#888888');
}

// --- Test 7: _lineageColor doesn't auto-trigger when fewer than 3 L2s
console.log('\nTest 7: auto-trigger guards on _LINEAGE_MIN_L2_FOR_COMPUTE');
{
  const sandbox = makeSandbox({
    data: { l2_envelopes: [{}, {}] },     // only 2 L2s
  });
  vm.runInContext('var c = _lineageColor(0);', sandbox);
  ok('returns null',  sandbox.c === null);
  ok('does NOT schedule auto-trigger',
     !sandbox.state._lineageComputeScheduled);
}

// =============================================================================
// Final tally
// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
