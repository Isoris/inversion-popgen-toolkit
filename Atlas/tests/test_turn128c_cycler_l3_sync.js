// =============================================================================
// turn 128c integration test — cycler → L3 compareUnit sync
//
// Quentin's bug report:
//   "When we select the window 5 or 10 in the multi button, it also
//    switches to 1 or 5 or 10 windows unit in the L3 contingency tables panel."
//
// (Quentin is describing the EXPECTED behaviour after the fix — clicking
//  the cycler should propagate to L3. Currently it doesn't, because
//  cycleStepSize was calling a never-defined `_mirrorStepModeIntoCompare`
//  helper instead of the real `_syncStepModeToCompareUnit`.)
//
// Tests:
//   1. cycleStepSize source no longer references the dead
//      _mirrorStepModeIntoCompare name.
//   2. cycleStepSize now calls _syncStepModeToCompareUnit.
//   3. _syncStepModeToCompareUnit, given state.stepModeSync=true and
//      state.stepMode=win5, sets state.compareUnit='win5' AND toggles the
//      .active class on the matching #l3CompareUnit button AND persists
//      to localStorage AND calls renderL3Panel.
//   4. _syncStepModeToCompareUnit is a no-op when state.stepModeSync=false.
//   5. The cycler's stepMode value is also written to
//      'pca_scrubber_v3.stepmode' (not just _STEP_CYCLE_KEY), so a reload
//      restores the same state both surfaces use.
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

ok('cycleStepSize no longer calls the dead _mirrorStepModeIntoCompare',
   !/cycleStepSize[\s\S]{0,800}_mirrorStepModeIntoCompare/.test(html),
   'this name was referenced but never defined; should be replaced');

ok('cycleStepSize calls _syncStepModeToCompareUnit (the real sync helper)',
   /cycleStepSize[\s\S]{0,2000}_syncStepModeToCompareUnit\(\)/.test(html));

ok('cycleStepSize persists stepMode to pca_scrubber_v3.stepmode (not just cycle key)',
   /cycleStepSize[\s\S]{0,2000}localStorage\.setItem\(['"]pca_scrubber_v3\.stepmode['"]/.test(html));

ok('cycleStepSize wraps _syncStepModeToCompareUnit in try/catch (defensive)',
   /typeof _syncStepModeToCompareUnit === 'function'\) \{\s*try \{ _syncStepModeToCompareUnit\(\); \}/.test(html));

// =============================================================================
// Behavioural test: extract _syncStepModeToCompareUnit + cycleStepSize
// into a sandbox and verify state propagation.
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

const fnSync = pullFunction(html, '_syncStepModeToCompareUnit');
const fnCycle = pullFunction(html, 'cycleStepSize');
const orderMatch = html.match(/^const _STEP_CYCLE_ORDER = \[[\s\S]*?\];/m);
const cycleKeyMatch = html.match(/^const _STEP_CYCLE_KEY = .*?;/m);

ok('_syncStepModeToCompareUnit extractable', !!fnSync);
ok('cycleStepSize extractable',               !!fnCycle);
ok('_STEP_CYCLE_ORDER extractable',           !!orderMatch);
ok('_STEP_CYCLE_KEY extractable',             !!cycleKeyMatch);

function makeBtn(unit) {
  const classes = new Set();
  return {
    dataset: { l3unit: unit },
    classList: {
      toggle: (c, force) => { force ? classes.add(c) : classes.delete(c); },
      contains: (c) => classes.has(c),
    },
    _classes: classes,
  };
}

function makeStepModeBtn(step) {
  const classes = new Set();
  return {
    dataset: { step: step },
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

function makeStorage() {
  const map = new Map();
  return {
    setItem: (k, v) => map.set(k, String(v)),
    getItem: (k) => map.has(k) ? map.get(k) : null,
    removeItem: (k) => map.delete(k),
    _map: map,
  };
}

function makeBtn0(id) {
  const classes = new Set();
  return {
    id, textContent: '',
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

function makeSandbox(initialState) {
  const l3Btns = ['L2', 'win1', 'win5', 'win10', 'winN'].map(makeBtn);
  const sidebarBtns = ['l2', 'win1', 'win5', 'win10', 'winN'].map(makeStepModeBtn);
  const cyclerBtn = makeBtn0('jumpToWindowsBtn');
  const stepModeInfo = { textContent: '' };
  const cuInput = { value: '' };
  const renderCalls = [];

  const sandbox = {
    state: Object.assign({
      stepMode: 'l2',
      stepModeN: 5,
      stepModeSync: true,
      compareUnit: 'L2',
      compareUnitN: 5,
    }, initialState || {}),
    document: {
      getElementById: function (id) {
        if (id === 'jumpToWindowsBtn') return cyclerBtn;
        if (id === 'l3CompareUnitN')   return cuInput;
        if (id === 'stepModeInfo')     return stepModeInfo;
        return null;
      },
      querySelectorAll: function (sel) {
        if (sel === '#l3CompareUnit button') return l3Btns;
        if (sel === '#stepModeBar button')   return sidebarBtns;
        return [];
      },
    },
    localStorage: makeStorage(),
    renderL3Panel: function () { renderCalls.push('renderL3Panel'); },
    _stepModeLabel: function (m) { return 'lbl-' + m; },
    console,
    Set, Map, Array, Object, String, Number, Math,
  };
  vm.createContext(sandbox);
  // Define a stub for _mirrorStepModeFromState so cycleStepSize works
  vm.runInContext(`function _mirrorStepModeFromState() {}`, sandbox);
  vm.runInContext(`function _refreshStepSizeBtn() {}`, sandbox);
  vm.runInContext(orderMatch[0],     sandbox);
  vm.runInContext(cycleKeyMatch[0],  sandbox);
  vm.runInContext(fnSync,            sandbox);
  vm.runInContext(fnCycle,           sandbox);
  return { sandbox, l3Btns, sidebarBtns, renderCalls };
}

// --- Test: cycler from win1 → win5 propagates to L3 compareUnit
console.log('\nTest: cycler propagates win1 → win5 to L3');
{
  const { sandbox, l3Btns, renderCalls } = makeSandbox({ stepMode: 'win1', compareUnit: 'win1' });
  // Fix: need a starting compareUnit (sync uses it as the "from" baseline)
  vm.runInContext('cycleStepSize();', sandbox);
  ok('state.stepMode advanced to win5',
     sandbox.state.stepMode === 'win5',
     'got ' + sandbox.state.stepMode);
  ok('state.compareUnit advanced to win5',
     sandbox.state.compareUnit === 'win5',
     'got ' + sandbox.state.compareUnit);
  // L3 button for win5 should be active
  const win5Btn = l3Btns.find(b => b.dataset.l3unit === 'win5');
  ok('L3 #l3CompareUnit win5 button has .active',
     win5Btn && win5Btn.classList.contains('active'));
  // L3 button for win1 should NOT be active
  const win1Btn = l3Btns.find(b => b.dataset.l3unit === 'win1');
  ok('L3 #l3CompareUnit win1 button NO LONGER active',
     win1Btn && !win1Btn.classList.contains('active'));
  // localStorage persists pca_scrubber_v3.stepmode
  ok('pca_scrubber_v3.stepmode persisted',
     sandbox.localStorage.getItem('pca_scrubber_v3.stepmode') === 'win5');
  // localStorage persists pca_scrubber_v3.compareunit (via sync helper)
  ok('pca_scrubber_v3.compareunit persisted',
     sandbox.localStorage.getItem('pca_scrubber_v3.compareunit') === 'win5');
  // renderL3Panel was called by the sync helper
  ok('renderL3Panel called once',
     renderCalls.indexOf('renderL3Panel') >= 0);
}

// --- Test: cycler from win5 → win10 propagates
console.log('\nTest: cycler win5 → win10 propagates');
{
  const { sandbox } = makeSandbox({ stepMode: 'win5', compareUnit: 'win5' });
  vm.runInContext('cycleStepSize();', sandbox);
  ok('cycler advanced to win10', sandbox.state.stepMode === 'win10');
  ok('compareUnit synced to win10', sandbox.state.compareUnit === 'win10');
}

// --- Test: cycler from win10 → win1 (wrap-around) propagates
console.log('\nTest: cycler win10 → win1 (wrap)');
{
  const { sandbox } = makeSandbox({ stepMode: 'win10', compareUnit: 'win10' });
  vm.runInContext('cycleStepSize();', sandbox);
  ok('cycler wrapped to win1', sandbox.state.stepMode === 'win1');
  ok('compareUnit synced to win1', sandbox.state.compareUnit === 'win1');
}

// --- Test: cycler from L2 mode jumps to win1 + propagates
console.log('\nTest: cycler from L2 starts → win1 propagates');
{
  const { sandbox } = makeSandbox({ stepMode: 'l2', compareUnit: 'L2' });
  vm.runInContext('cycleStepSize();', sandbox);
  ok('cycler entered cycle at win1', sandbox.state.stepMode === 'win1');
  ok('compareUnit synced to win1', sandbox.state.compareUnit === 'win1');
}

// --- Test: when stepModeSync is off, cycler still updates state.stepMode
//          but does NOT propagate to L3
console.log('\nTest: sync OFF — cycler does not push to L3');
{
  const { sandbox, renderCalls } = makeSandbox({
    stepMode: 'win1', compareUnit: 'L2', stepModeSync: false,
  });
  vm.runInContext('cycleStepSize();', sandbox);
  ok('cycler still advanced state.stepMode to win5',
     sandbox.state.stepMode === 'win5');
  ok('compareUnit unchanged (sync off)',
     sandbox.state.compareUnit === 'L2');
  // renderL3Panel should NOT have been called from the sync (it only runs
  // when sync is on)
  ok('renderL3Panel NOT called when sync off',
     renderCalls.indexOf('renderL3Panel') < 0);
}

// --- Test: also verify the bare _syncStepModeToCompareUnit independently
console.log('\nTest: _syncStepModeToCompareUnit standalone');
{
  const { sandbox, l3Btns } = makeSandbox({
    stepMode: 'win10', compareUnit: 'L2', stepModeSync: true,
  });
  vm.runInContext('_syncStepModeToCompareUnit();', sandbox);
  ok('compareUnit pushed to win10', sandbox.state.compareUnit === 'win10');
  const win10Btn = l3Btns.find(b => b.dataset.l3unit === 'win10');
  ok('L3 win10 button active', win10Btn.classList.contains('active'));
}

// --- Test: winN with a custom N also propagates
console.log('\nTest: winN with custom N propagates compareUnitN');
{
  const { sandbox } = makeSandbox({
    stepMode: 'winN', stepModeN: 17, compareUnit: 'win1', compareUnitN: 5,
    stepModeSync: true,
  });
  vm.runInContext('_syncStepModeToCompareUnit();', sandbox);
  ok('compareUnit became winN', sandbox.state.compareUnit === 'winN');
  ok('compareUnitN became 17',  sandbox.state.compareUnitN === 17);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
