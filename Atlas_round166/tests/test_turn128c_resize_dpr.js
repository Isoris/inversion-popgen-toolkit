// =============================================================================
// turn 128c integration test — resize + DPR-change canvas refit
//
// Quentin's bug report:
//   "When we unzoom with CTRL - the pca panels of local PCA of L3
//    contingency panel are distorted and not rerendered?"
//
// Two failure modes addressed:
//   (a) Chrome/Edge sometimes doesn't fire `resize` on CTRL+- because window
//       dimensions in CSS pixels don't change — only devicePixelRatio does.
//       Fix: also subscribe to matchMedia('(resolution: <DPR>dppx)') and
//       trigger redraw on DPR change. Re-arm with the new DPR on each fire.
//   (b) Even when `resize` does fire, calling fitCanvas synchronously inside
//       the handler reads stale bounding rects during smooth-zoom animations.
//       Fix: requestAnimationFrame to let layout settle before re-drawing.
//
// Tests:
//   1. Source: _redrawAllOnResize is defined and registered on window resize.
//   2. Source: _installDprListener is defined and called at module init.
//   3. Source: redraw is wrapped in requestAnimationFrame.
//   4. Source: redraw is debounced by a flag (no double-schedule per frame).
//   5. Source: matchMedia uses the (resolution: Xdppx) media query.
//   6. Source: matchMedia listener re-installs itself on each fire (DPR
//      changes, so the previous mq becomes stale).
//   7. Source: matchMedia listener wrapped in try/catch for ancient browsers.
//   8. Source: handles both modern .addEventListener('change') and legacy
//      .addListener API.
//   9. Source: redraw still calls drawSim/drawZ/drawPCA/renderL3Panel.
//  10. Behavioural: with rAF stub, calling _redrawAllOnResize twice in the
//      same tick only schedules one rAF.
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

ok('_redrawAllOnResize defined',
   /function\s+_redrawAllOnResize\s*\(/.test(html));

ok('_installDprListener defined',
   /function\s+_installDprListener\s*\(/.test(html));

ok('window.addEventListener("resize", _redrawAllOnResize) wired',
   /window\.addEventListener\(['"]resize['"], _redrawAllOnResize\)/.test(html));

ok('_installDprListener invoked at module init',
   /^_installDprListener\(\);/m.test(html),
   'must be called once on script load to subscribe to the initial DPR');

ok('redraw wraps body in requestAnimationFrame',
   /_redrawAllOnResize[\s\S]{0,500}requestAnimationFrame\(/.test(html));

ok('redraw uses a debounce flag (_resizeRafScheduled)',
   /let\s+_resizeRafScheduled\s*=\s*false/.test(html));

ok('matchMedia uses (resolution: <DPR>dppx) query',
   /window\.matchMedia\('\(resolution: ' \+ dpr \+ 'dppx\)'\)/.test(html));

ok('matchMedia listener re-installs itself on each fire (recursion)',
   /_installDprListener[\s\S]{0,800}_installDprListener\(\);/.test(html),
   'the function must call itself inside its onChange callback to re-arm with the new DPR');

ok('matchMedia install wrapped in try/catch (degrade for old browsers)',
   /function\s+_installDprListener\s*\(\)\s*\{\s*try\s*\{/.test(html));

ok('matchMedia handles both modern (addEventListener) and legacy (addListener) APIs',
   /addEventListener\(['"]change['"]/.test(html) &&
   /addListener\(/.test(html));

ok('redraw calls drawSim/drawZ/drawPCA/renderL3Panel',
   /_redrawAllOnResize[\s\S]{0,800}drawSim\(\);[\s\S]{0,200}drawZ\(\);[\s\S]{0,200}drawPCA\(\);[\s\S]{0,200}renderL3Panel\(\);/.test(html));

ok('redraw early-exits when state.data is null',
   /_redrawAllOnResize[\s\S]{0,300}if \(!state\.data\) return;/.test(html));

// Old one-line resize handler (no rAF, no DPR) must be gone.
const oldHandler = /window\.addEventListener\('resize',\s*\(\)\s*=>\s*\{\s*if \(state\.data\) \{ drawSim\(\); drawZ\(\); drawPCA\(\); renderL3Panel\(\); \}\s*\}\);/;
ok('old synchronous resize handler is gone',
   !oldHandler.test(html));

// =============================================================================
// Behavioural test — debounce works
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
      const q = ch; i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === q) break;
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

const fnRedraw = pullFunction(html, '_redrawAllOnResize');
ok('_redrawAllOnResize extractable', !!fnRedraw);

// --- Test: rAF debounce — two synchronous calls = one rAF
console.log('\nTest: rAF debounce coalesces same-tick calls');
{
  let rafScheduled = 0;
  const rafQueue = [];
  const sandbox = {
    state: { data: { foo: 1 } },
    drawSim:       () => {},
    drawZ:         () => {},
    drawPCA:       () => {},
    renderL3Panel: () => {},
    requestAnimationFrame: (cb) => { rafScheduled++; rafQueue.push(cb); return rafScheduled; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext('let _resizeRafScheduled = false;', sandbox);
  vm.runInContext(fnRedraw, sandbox);

  vm.runInContext('_redrawAllOnResize(); _redrawAllOnResize(); _redrawAllOnResize();', sandbox);
  ok('three calls in same tick → exactly 1 rAF scheduled', rafScheduled === 1,
     'got ' + rafScheduled);

  // Now flush the rAF; the flag should reset and new calls trigger again
  rafQueue[0]();
  vm.runInContext('_redrawAllOnResize();', sandbox);
  ok('after rAF flush, next call schedules a new rAF', rafScheduled === 2);
}

// --- Test: redraw early-exits when state.data is null (no rAF scheduled)
console.log('\nTest: no state.data → no work scheduled');
{
  let rafScheduled = 0;
  const sandbox = {
    state: { data: null },
    drawSim: () => {}, drawZ: () => {}, drawPCA: () => {}, renderL3Panel: () => {},
    requestAnimationFrame: () => { rafScheduled++; return 1; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext('let _resizeRafScheduled = false;', sandbox);
  vm.runInContext(fnRedraw, sandbox);
  vm.runInContext('_redrawAllOnResize();', sandbox);
  ok('no rAF when state.data is null', rafScheduled === 0);
}

// --- Test: rAF callback resets debounce flag and calls all four draw fns
console.log('\nTest: rAF callback runs all four draw functions');
{
  let rafScheduled = 0;
  const rafQueue = [];
  const calls = [];
  const sandbox = {
    state: { data: {} },
    drawSim:       () => calls.push('drawSim'),
    drawZ:         () => calls.push('drawZ'),
    drawPCA:       () => calls.push('drawPCA'),
    renderL3Panel: () => calls.push('renderL3Panel'),
    requestAnimationFrame: (cb) => { rafScheduled++; rafQueue.push(cb); return rafScheduled; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext('let _resizeRafScheduled = false;', sandbox);
  vm.runInContext(fnRedraw, sandbox);
  vm.runInContext('_redrawAllOnResize();', sandbox);
  rafQueue[0]();   // flush
  ok('all four draw functions called',
     calls.length === 4 &&
     calls.indexOf('drawSim') >= 0 &&
     calls.indexOf('drawZ') >= 0 &&
     calls.indexOf('drawPCA') >= 0 &&
     calls.indexOf('renderL3Panel') >= 0,
     'calls: ' + JSON.stringify(calls));
  // The debounce flag should be cleared so subsequent calls re-arm
  vm.runInContext('_redrawAllOnResize();', sandbox);
  ok('after rAF flush, flag cleared (new rAF scheduled)', rafScheduled === 2);
}

// --- Test: each draw call's throw doesn't break the others (try/catch each)
console.log('\nTest: drawSim throw doesnt break the chain');
{
  const rafQueue = [];
  const calls = [];
  const sandbox = {
    state: { data: {} },
    drawSim:       () => { throw new Error('boom drawSim'); },
    drawZ:         () => calls.push('drawZ'),
    drawPCA:       () => calls.push('drawPCA'),
    renderL3Panel: () => calls.push('renderL3Panel'),
    requestAnimationFrame: (cb) => { rafQueue.push(cb); return 1; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext('let _resizeRafScheduled = false;', sandbox);
  vm.runInContext(fnRedraw, sandbox);
  vm.runInContext('_redrawAllOnResize();', sandbox);
  let threw = false;
  try { rafQueue[0](); } catch (_) { threw = true; }
  ok('rAF callback does not propagate drawSim throw', !threw);
  ok('subsequent draws still happened despite drawSim throw',
     calls.indexOf('drawZ') >= 0 &&
     calls.indexOf('drawPCA') >= 0 &&
     calls.indexOf('renderL3Panel') >= 0);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
