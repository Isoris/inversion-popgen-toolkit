// =============================================================================
// turn 128 integration test — window-mode cycler "lit" state + sidebar sync
//
// Bug from the latest user message:
//   "When we are in window mode for cursor navigation in page 1 we have 3
//   choices window (1) window (5) and window (10) and it works great but
//   their button is not lighting up when its active and selected"
//
// Two layers of fix:
//   (B1) Cycler #jumpToWindowsBtn now toggles class .is-active-step when
//        state.stepMode is one of win1 / win5 / win10. New CSS rule lights
//        it up (accent border + accent text) without making it look like a
//        zoom button (which would conflict with Quentin's earlier "keep them
//        visually distinct" intent).
//   (B2) _mirrorStepModeFromState was previously REFERENCED but never
//        DEFINED — meaning every header-cycler click silently failed to
//        sync the sidebar #stepModeBar buttons. Now defined; sidebar tracks
//        every cycler click. Conversely, sidebar-button clicks now also
//        refresh the cycler so both surfaces stay reflective.
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

ok('CSS rule for #jumpToWindowsBtn.is-active-step exists',
   /#viewModeBar\s+#jumpToWindowsBtn\.is-active-step\s*\{/.test(html));

ok('_mirrorStepModeFromState is now DEFINED (not just referenced)',
   /function\s+_mirrorStepModeFromState\s*\(/.test(html));

ok('_refreshStepSizeBtn toggles is-active-step class',
   /_refreshStepSizeBtn[\s\S]{0,2500}classList\.toggle\(['"]is-active-step['"]/.test(html));

ok('sidebar #stepModeBar click handler now also calls _refreshStepSizeBtn',
   /querySelectorAll\(['"]#stepModeBar button['"]\)\.forEach\(btn[\s\S]{0,1200}_refreshStepSizeBtn\(\)/.test(html),
   'sidebar click should sync the header cycler');

// turn 147b updated _refreshStepSizeBtn to handle l2/winN modes
// distinctly. The label is built via a `label` variable (not inline
// template literal) — internal restructure, same output format.
ok('cycler still emits "📊 Windows (N)" format for cycler-active modes',
   /label\s*=\s*`📊 Windows \(\$\{n\}\)`/.test(html));

// =============================================================================
// Behavioural tests
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

// Pull _STEP_CYCLE_ORDER
const orderMatch = html.match(/^const _STEP_CYCLE_ORDER = \[[\s\S]*?\];/m);
const orderSrc = orderMatch ? orderMatch[0] : '';
const stepSizeFn = pullFunction(html, '_stepSizeForMode');
const refreshFn  = pullFunction(html, '_refreshStepSizeBtn');
const mirrorFn   = pullFunction(html, '_mirrorStepModeFromState');

ok('_STEP_CYCLE_ORDER constant extractable', orderSrc !== '');
ok('_stepSizeForMode extractable',           !!stepSizeFn);
ok('_refreshStepSizeBtn extractable',        !!refreshFn);
ok('_mirrorStepModeFromState extractable',   !!mirrorFn);

// Build a sandbox with a fake DOM. We need:
//   - getElementById('jumpToWindowsBtn') returning a fake btn with classList,
//     textContent, etc.
//   - querySelectorAll('#stepModeBar button') returning fake sidebar btns
//   - getElementById('stepModeInfo') returning a fake info span
function makeBadge() {
  const classes = new Set();
  return {
    textContent: '',
    classList: {
      add:    (c) => classes.add(c),
      remove: (c) => classes.delete(c),
      toggle: (c, force) => {
        if (typeof force === 'boolean') {
          if (force) classes.add(c); else classes.delete(c);
        } else {
          if (classes.has(c)) classes.delete(c); else classes.add(c);
        }
      },
      contains: (c) => classes.has(c),
    },
    _classes: classes,
  };
}

function makeSidebarBtn(stepMode) {
  const classes = new Set();
  return {
    dataset: { step: stepMode },
    classList: {
      add:    (c) => classes.add(c),
      remove: (c) => classes.delete(c),
      toggle: (c, force) => {
        if (typeof force === 'boolean') {
          if (force) classes.add(c); else classes.delete(c);
        } else {
          if (classes.has(c)) classes.delete(c); else classes.add(c);
        }
      },
      contains: (c) => classes.has(c),
    },
    _classes: classes,
  };
}

function makeSandbox(stepMode) {
  const cyclerBtn = makeBadge();
  const sidebarBtns = ['l2', 'win1', 'win5', 'win10', 'winN'].map(makeSidebarBtn);
  const infoSpan = { textContent: '' };

  const fakeDoc = {
    getElementById: function (id) {
      if (id === 'jumpToWindowsBtn') return cyclerBtn;
      if (id === 'stepModeInfo')      return infoSpan;
      return null;
    },
    querySelectorAll: function (sel) {
      if (sel === '#stepModeBar button') return sidebarBtns;
      return [];
    },
  };

  const sandbox = {
    state:    { stepMode: stepMode },
    document: fakeDoc,
    console:  console,
    _stepModeLabel: function (m) { return 'label-for-' + m; },
    Set, Map, Object, Array, String, Number,
  };
  vm.createContext(sandbox);
  vm.runInContext(orderSrc,  sandbox);
  vm.runInContext(stepSizeFn, sandbox);
  vm.runInContext(refreshFn, sandbox);
  vm.runInContext(mirrorFn,  sandbox);
  return { sandbox, cyclerBtn, sidebarBtns, infoSpan };
}

// --- Test: cycler lights up for win1
console.log('\nTest: cycler lit for win1');
{
  const { sandbox, cyclerBtn } = makeSandbox('win1');
  vm.runInContext('_refreshStepSizeBtn();', sandbox);
  ok('text reads "📊 Windows (1)"',
     cyclerBtn.textContent === '📊 Windows (1)',
     'got "' + cyclerBtn.textContent + '"');
  ok('class is-active-step is set for win1',
     cyclerBtn.classList.contains('is-active-step'));
}

// --- Test: cycler lights up for win5
console.log('\nTest: cycler lit for win5');
{
  const { sandbox, cyclerBtn } = makeSandbox('win5');
  vm.runInContext('_refreshStepSizeBtn();', sandbox);
  ok('text reads "📊 Windows (5)"',
     cyclerBtn.textContent === '📊 Windows (5)');
  ok('class is-active-step is set for win5',
     cyclerBtn.classList.contains('is-active-step'));
}

// --- Test: cycler lights up for win10
console.log('\nTest: cycler lit for win10');
{
  const { sandbox, cyclerBtn } = makeSandbox('win10');
  vm.runInContext('_refreshStepSizeBtn();', sandbox);
  ok('text reads "📊 Windows (10)"',
     cyclerBtn.textContent === '📊 Windows (10)');
  ok('class is-active-step is set for win10',
     cyclerBtn.classList.contains('is-active-step'));
}

// --- Test: cycler dims for L2 mode (sidebar selected L2)
// turn 147b — label switched from "📊 Windows (1)" to "📊 Windows (—)" so
// the user can tell at a glance that the cycler is NOT the active step
// setting (arrows jump L2 envelopes, not windows). Pre-147b the "(1)"
// label was misleading — looked like step=1 was active when it wasn't.
console.log('\nTest: cycler dim when stepMode=l2 (sidebar took over)');
{
  const { sandbox, cyclerBtn } = makeSandbox('l2');
  vm.runInContext('_refreshStepSizeBtn();', sandbox);
  ok('text shows "📊 Windows (—)" (cycler not the live setting)',
     cyclerBtn.textContent === '📊 Windows (—)');
  ok('class is-active-step is NOT set for l2',
     !cyclerBtn.classList.contains('is-active-step'),
     'l2 isnt in the cycler order so cycler should be dim');
}

// --- Test: cycler engages for winN mode (sidebar's custom N)
// turn 147b — Quentin's bug fix. Pre-147b the cycler stayed dim for
// winN, hiding the fact that arrows ARE stepping in window units. Now
// the cycler lights up AND shows the custom N value (e.g. "📊 Windows
// (15)") so the user knows their custom step size is live.
console.log('\nTest: cycler ENGAGES when stepMode=winN (turn 147b fix)');
{
  const { sandbox, cyclerBtn } = makeSandbox('winN');
  // Set a custom N value
  sandbox.state.stepModeN = 15;
  vm.runInContext('_refreshStepSizeBtn();', sandbox);
  ok('class is-active-step IS set for winN (turn 147b fix)',
     cyclerBtn.classList.contains('is-active-step'));
  ok('label reflects the custom N value',
     cyclerBtn.textContent === '📊 Windows (15)');
}

// --- Test: _mirrorStepModeFromState toggles sidebar buttons correctly
console.log('\nTest: _mirrorStepModeFromState mirrors state to sidebar');
{
  const { sandbox, sidebarBtns, infoSpan } = makeSandbox('win5');
  vm.runInContext('_mirrorStepModeFromState();', sandbox);
  // Only the win5 button should have .active
  const activeBtns = sidebarBtns.filter(b => b.classList.contains('active'));
  ok('exactly one sidebar button is active', activeBtns.length === 1,
     'got ' + activeBtns.length + ' active');
  ok('the active sidebar button is the win5 one',
     activeBtns[0] && activeBtns[0].dataset.step === 'win5',
     activeBtns[0] ? 'got ' + activeBtns[0].dataset.step : 'no active btn');
  ok('stepModeInfo text reflects current mode',
     infoSpan.textContent === 'label-for-win5',
     'got "' + infoSpan.textContent + '"');
}

// --- Test: mirror works for L2 and winN too
console.log('\nTest: mirror works for L2 / winN');
{
  for (const mode of ['l2', 'win1', 'win10', 'winN']) {
    const { sandbox, sidebarBtns } = makeSandbox(mode);
    vm.runInContext('_mirrorStepModeFromState();', sandbox);
    const activeBtns = sidebarBtns.filter(b => b.classList.contains('active'));
    ok('mode=' + mode + ': exactly the matching sidebar button is active',
       activeBtns.length === 1 && activeBtns[0].dataset.step === mode);
  }
}

// --- Test: graceful when btn missing (no-op, doesn't throw)
console.log('\nTest: graceful when DOM elements absent');
{
  const sandbox = {
    state: { stepMode: 'win1' },
    document: {
      getElementById: () => null,
      querySelectorAll: () => [],
    },
    console: console,
    _stepModeLabel: function (m) { return ''; },
  };
  vm.createContext(sandbox);
  vm.runInContext(orderSrc,  sandbox);
  vm.runInContext(stepSizeFn, sandbox);
  vm.runInContext(refreshFn, sandbox);
  vm.runInContext(mirrorFn,  sandbox);
  let threw = false;
  try {
    vm.runInContext('_refreshStepSizeBtn(); _mirrorStepModeFromState();', sandbox);
  } catch (e) { threw = true; }
  ok('no throw when buttons absent', !threw);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
