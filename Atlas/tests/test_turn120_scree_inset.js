// =============================================================================
// turn 120 integration test — PCA scree-plot inset on the tracked-sample plane
//
// Tests:
//   1. Source-level: helpers, CSS classes, toggle wired to drawPCA.
//   2. Renderer behavior under three precomp scenarios:
//        a. Full lam_top_k array of 5 eigenvalues (schema 2.16) → 5 bars.
//        b. Legacy precomp with only lam1/lam2 → 2 bars + fallback hint.
//        c. Neither field → empty render (returns '').
//   3. Empty states: scree disabled, no data loaded, no window selected.
//   4. Defensive sort: lam_top_k that's not already sorted descending
//      gets sorted before rendering.
//   5. Cap at 7 bars max (don't blow up the inset width on a freak large k).
//   6. Bar fill colors match the documented palette for at least PC1 and PC2.
//   7. λ₁/λ₂ ratio line: shown with 1 decimal precision.
//   8. Fallback hint appears only when isFallback path is taken.
// =============================================================================

const fs = require('fs');
const vm = require('vm');

const html = fs.readFileSync('/home/claude/work/build/Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// --- Source-level checks
console.log('\n=== Source-level checks ===');
ok('_renderScreeInsetHTML defined',  /function _renderScreeInsetHTML\(/.test(html));
ok('_refreshScreeInset defined',     /function _refreshScreeInset\(/.test(html));
ok('_SCREE_BAR_COLORS defined',      /const _SCREE_BAR_COLORS = \[/.test(html));
ok('state.screePlotEnabled init',    /state\.screePlotEnabled\s*=/.test(html));
ok('localStorage persist key',
   html.indexOf("'inversion_atlas.screePlotEnabled'") >= 0);
ok('toggle handler wired to checkbox',
   /screeToggle.*addEventListener.*'change'/s.test(html));
ok('_refreshScreeInset hooked into drawPCA',
   /try \{ _refreshScreeInset\(\); \}/.test(html));
ok('screeInset div in pcaCanvasWrap',
   /<div id="screeInset" class="scree-inset"/.test(html));
ok('scree toggle checkbox in pcaTrackedAside',
   /<input type="checkbox" id="screeToggle"/.test(html));

// --- CSS contract
console.log('\n=== CSS contract ===');
const cssClasses = [
  'scree-inset', 'scree-inset-label', 'scree-inset-svg',
  'scree-inset-bar', 'scree-inset-ratio', 'scree-inset-fallback-hint',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined',
     html.indexOf('.' + cls) >= 0);
}

// =============================================================================
// Behavioural test: pull the renderer into a sandbox + exercise scenarios.
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\$&') + '\\s*\\(', 'm');
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
        if (quote === '`' && src[i] === '$' && src[i+1] === '{') {
          i += 2; let d = 1;
          while (i < src.length && d > 0) {
            if (src[i] === '{') d++;
            else if (src[i] === '}') d--;
            i++;
          }
          continue;
        }
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

// Pull the renderer
const fn = pullFunction(html, '_renderScreeInsetHTML');
if (!fn) {
  fail++; console.log('  FAIL could not extract _renderScreeInsetHTML');
} else {
  // Pull the palette constant
  const constMatch = html.match(/^const _SCREE_BAR_COLORS = \[[\s\S]*?\];/m);
  const constSrc = constMatch ? constMatch[0] : '';

  function makeSandbox(stateOverlay) {
    const state = Object.assign({
      data: null, cur: -1, screePlotEnabled: false,
    }, stateOverlay);
    const sandbox = {
      state, console,
      Number, Array, Math, JSON, Object, String, isNaN, isFinite,
    };
    vm.createContext(sandbox);
    if (constSrc) vm.runInContext(constSrc, sandbox);
    vm.runInContext(fn, sandbox);
    return sandbox;
  }

  // --- Test 1: scree disabled
  console.log('\nTest 1: scree disabled → empty');
  const s1 = makeSandbox({ screePlotEnabled: false,
    data: { windows: [{ lam1: 0.8, lam2: 0.3, lam_top_k: [0.8, 0.3, 0.2] }] }, cur: 0 });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s1);
  ok('disabled toggle → empty string', s1.__out === '');

  // --- Test 2: no data
  console.log('\nTest 2: no data loaded → empty');
  const s2 = makeSandbox({ screePlotEnabled: true, data: null });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s2);
  ok('no data → empty string', s2.__out === '');

  // --- Test 3: no window selected
  console.log('\nTest 3: cur < 0 → empty');
  const s3 = makeSandbox({ screePlotEnabled: true,
    data: { windows: [{ lam1: 0.8, lam2: 0.3 }] }, cur: -1 });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s3);
  ok('cur=-1 → empty string', s3.__out === '');

  // --- Test 4: full lam_top_k (5 eigenvalues, schema 2.16)
  console.log('\nTest 4: full lam_top_k → 5 bars');
  const s4 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.82, lam2: 0.31,
                        lam_top_k: [0.82, 0.31, 0.18, 0.12, 0.08] }] },
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s4);
  ok('full spectrum: returns string', typeof s4.__out === 'string' && s4.__out.length > 0);
  const nRects4 = (s4.__out.match(/<rect/g) || []).length;
  ok('full spectrum: 5 <rect> elements', nRects4 === 5, 'got ' + nRects4);
  ok('full spectrum: NO fallback hint',
     s4.__out.indexOf('scree-inset-fallback-hint') < 0);
  ok('full spectrum: ratio λ₁/λ₂ shown',
     /λ₁\/λ₂ = 2\.6/.test(s4.__out),
     'expected "λ₁/λ₂ = 2.6" (0.82/0.31 = 2.645)');

  // --- Test 5: legacy precomp (lam1/lam2 only)
  console.log('\nTest 5: legacy [lam1, lam2] → 2 bars + hint');
  const s5 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.65, lam2: 0.42 }] },   // no lam_top_k
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s5);
  const nRects5 = (s5.__out.match(/<rect/g) || []).length;
  ok('legacy: 2 <rect> elements', nRects5 === 2, 'got ' + nRects5);
  ok('legacy: fallback hint present',
     s5.__out.indexOf('scree-inset-fallback-hint') >= 0);
  ok('legacy: hint mentions schema 2.16',
     /precomp ≥2\.16/.test(s5.__out));

  // --- Test 6: neither field → empty
  console.log('\nTest 6: no eigenvalues at all → empty');
  const s6 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ z: 1.5 }] },   // window object with no eigenvalues
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s6);
  ok('no eigenvalues → empty string', s6.__out === '');

  // --- Test 7: defensive sort — unsorted lam_top_k gets sorted
  console.log('\nTest 7: defensive descending sort');
  const s7 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.30, lam2: 0.80,
                        lam_top_k: [0.30, 0.80, 0.50, 0.10, 0.25] }] },
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s7);
  // After defensive sort, descending = [0.80, 0.50, 0.30, 0.25, 0.10]
  // The first bar (PC1, color #7ad3db) should be the tallest.
  // We can verify by checking the inset shows λ₁/λ₂ = 0.80/0.50 = 1.6
  ok('unsorted input: ratio uses sorted values',
     /λ₁\/λ₂ = 1\.6/.test(s7.__out),
     'expected "1.6" from sorted [0.80, 0.50, ...] not from raw [0.30, 0.80, ...]');

  // --- Test 8: cap at 7 bars max
  console.log('\nTest 8: cap at 7 bars max');
  const s8 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 1.0, lam2: 0.5,
                        lam_top_k: [1.0, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.08, 0.05, 0.02] }] },
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s8);
  const nRects8 = (s8.__out.match(/<rect/g) || []).length;
  ok('cap: 10-element lam_top_k → 7 bars', nRects8 === 7, 'got ' + nRects8);

  // --- Test 9: NaN/zero filter
  console.log('\nTest 9: NaN/zero filtering');
  const s9 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.8, lam2: 0.3,
                        lam_top_k: [0.8, 0.3, NaN, 0, -0.001, 0.1] }] },
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s9);
  // Filtered: [0.8, 0.3, 0.1] (NaN, 0, -0.001 dropped) → 3 bars
  const nRects9 = (s9.__out.match(/<rect/g) || []).length;
  ok('NaN/zero filtered: 3 bars from valid subset', nRects9 === 3, 'got ' + nRects9);

  // --- Test 10: PC1 bar color matches palette
  console.log('\nTest 10: PC1 bar uses #7ad3db');
  ok('PC1 color in output', s4.__out.indexOf('#7ad3db') >= 0);

  // --- Test 11: window index in label
  console.log('\nTest 11: label includes "w<index>"');
  ok('label has window index',
     /PC eigenvalues · w0/.test(s4.__out));

  // --- Test 12: out-of-range cur
  console.log('\nTest 12: cur out of range → empty');
  const s12 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.8, lam2: 0.3 }] },
    cur: 5,   // only 1 window in array
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s12);
  ok('cur > windows.length-1 → empty string', s12.__out === '');

  // --- Test 13: only one positive eigenvalue → empty (need at least 2)
  console.log('\nTest 13: < 2 valid eigenvalues → empty');
  const s13 = makeSandbox({
    screePlotEnabled: true,
    data: { windows: [{ lam1: 0.8, lam2: NaN }] },
    cur: 0,
  });
  vm.runInContext('var __out = _renderScreeInsetHTML();', s13);
  ok('only 1 valid eigenvalue → empty string', s13.__out === '');
}

// =============================================================================
console.log('\n=========================================');
console.log('passed: ' + pass + ' / failed: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
