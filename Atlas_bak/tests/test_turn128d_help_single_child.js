// =============================================================================
// turn 128d integration test — single-child stage pill collapse
//
// Quentin's request:
//   "The help category has just the help button so why not just remove the
//   help button and leave it as the category button since it just expands
//   itself. Like Help > Help. Just having Help is enough."
//
// Implementation: generic. Any stage pill whose data-stage matches exactly
// one page-tab gets data-single-child="1". CSS hides the redundant child
// tab + the pill's expand arrow + its child-count chip. Clicking the pill
// programmatically clicks the (hidden) child tab so the page switch still
// happens through the existing handler.
//
// Tests:
//   1. CSS: .pill-arrow + .pill-count are display:none on data-single-child="1"
//   2. CSS: button[data-page][data-single-child="1"] is display:none
//   3. JS: _initTabGroupPills detects single-child stages by counting tabs
//      with the same data-stage and tags them.
//   4. JS: clicking a single-child pill calls .click() on the child tab
//      (page switch goes through existing tab handler).
//   5. JS: clicking a multi-child pill still calls _setActiveStage (the
//      legacy expand behaviour is preserved for non-single stages).
//   6. JS: tab-button click sets data-as-page="1" on the matching single-
//      child parent pill, clears it on others.
//   7. The 'help' stage in the actual atlas markup is the canonical test
//      target — it has exactly one child (page5).
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
// CSS contract
// =============================================================================
console.log('\n=== CSS contract ===');

ok('CSS: pill-arrow hidden on data-single-child="1"',
   /\.tab-stage-pill\[data-single-child="1"\] \.pill-arrow\b[\s\S]*?display:\s*none/.test(html));

ok('CSS: pill-count hidden on data-single-child="1"',
   /\.tab-stage-pill\[data-single-child="1"\] \.pill-count\b[\s\S]*?display:\s*none/.test(html));

ok('CSS: button[data-page][data-single-child="1"] hidden permanently',
   /#tabBar button\[data-page\]\[data-single-child="1"\]\s*\{\s*display:\s*none/.test(html));

ok('CSS: pill[data-as-page="1"] gets active-tab styling',
   /\.tab-stage-pill\[data-single-child="1"\]\[data-as-page="1"\][\s\S]*?font-weight:\s*600/.test(html));

// =============================================================================
// Source: _initTabGroupPills changes
// =============================================================================
console.log('\n=== _initTabGroupPills updates ===');

ok('_initTabGroupPills queries tabs per stage (single-child detection)',
   /_initTabGroupPills[\s\S]{0,1500}#tabBar button\[data-page\]\[data-stage="['"\s\S]*?singleChild = '1'/.test(html));

ok('_initTabGroupPills marks both pill and child with data-single-child="1"',
   /pill\.dataset\.singleChild = '1';\s*\n\s*tabs\[0\]\.dataset\.singleChild = '1'/.test(html));

ok('Pill click routes single-child stages to child.click()',
   /const child = singleChildMap\.get\(stage\);\s*\n\s*if \(child\) \{\s*\n\s*child\.click\(\);\s*\n\s*return;/.test(html));

ok('Pill click for multi-child stages still calls _setActiveStage (preserve legacy)',
   /const child = singleChildMap\.get\(stage\);[\s\S]{0,200}\}\s*\n\s*_setActiveStage\(stage\);/.test(html));

// =============================================================================
// Source: tab-button click sets data-as-page on single-child pill
// =============================================================================
console.log('\n=== tab-button click → data-as-page mirror ===');

ok('Tab click handler iterates all stage pills and toggles data-as-page',
   /document\.querySelectorAll\('#tabBar \.tab-stage-pill'\)\.forEach\(p => \{\s*\n\s*const isActiveSingle = \(p\.dataset\.singleChild === '1' && p\.dataset\.stage === btn\.dataset\.stage\)/.test(html));

ok('Tab click handler sets data-as-page="1" on matching single-child pill',
   /if \(isActiveSingle\) p\.dataset\.asPage = '1';/.test(html));

ok('Tab click handler clears data-as-page on non-matching pills',
   /else delete p\.dataset\.asPage/.test(html));

// =============================================================================
// Behavioural: simulate _initTabGroupPills against a fake DOM
// =============================================================================
console.log('\n=== Behavioural tests (fake DOM) ===');

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

const fnInit = pullFunction(html, '_initTabGroupPills');
ok('_initTabGroupPills extractable', !!fnInit);

// Fake DOM helpers
function makePill(stage) {
  const dataset = { stage };
  const handlers = [];
  return {
    dataset,
    addEventListener: (evt, cb) => { if (evt === 'click') handlers.push(cb); },
    click() { handlers.forEach(h => h()); },
  };
}

function makeTab(stage, pageId) {
  const dataset = { stage, page: pageId };
  let clickedTimes = 0;
  return {
    dataset,
    click() { clickedTimes++; },
    get _clicks() { return clickedTimes; },
  };
}

function makeFakeDoc(pills, tabs) {
  return {
    querySelector: () => null,
    querySelectorAll: function (sel) {
      if (sel === '#tabBar .tab-stage-pill') return pills;
      // Match #tabBar button[data-page][data-stage="<X>"]
      const m = sel.match(/data-stage="([^"]+)"/);
      if (m) {
        const wanted = m[1];
        return tabs.filter(t => t.dataset.stage === wanted);
      }
      return [];
    },
  };
}

// --- Test: pills with single child get data-single-child="1"
console.log('\nTest: single-child stages tagged correctly');
{
  const pills = [
    makePill('discovery'), makePill('refinement'), makePill('help'),
  ];
  const tabs = [
    makeTab('discovery', 'page1'), makeTab('discovery', 'page12'),  // 2
    makeTab('refinement', 'page11'),                                // 1 (single)
    makeTab('help', 'page5'),                                       // 1 (single)
  ];
  const sandbox = {
    document: makeFakeDoc(pills, tabs),
    localStorage: { getItem: () => null, setItem: () => {} },
    console,
    Map, Set, Array,
  };
  vm.createContext(sandbox);
  vm.runInContext('function _setActiveStage() {}', sandbox);
  vm.runInContext(fnInit, sandbox);
  vm.runInContext('_initTabGroupPills();', sandbox);

  ok('discovery pill NOT tagged (2 children)',
     pills[0].dataset.singleChild !== '1');
  ok('refinement pill tagged (1 child)',
     pills[1].dataset.singleChild === '1');
  ok('help pill tagged (1 child)',
     pills[2].dataset.singleChild === '1');
  ok('refinement child tab tagged',
     tabs[2].dataset.singleChild === '1');
  ok('help child tab tagged',
     tabs[3].dataset.singleChild === '1');
  ok('discovery children NOT tagged',
     tabs[0].dataset.singleChild !== '1' && tabs[1].dataset.singleChild !== '1');
}

// --- Test: clicking a single-child pill routes to child.click()
console.log('\nTest: single-child pill click → child tab click');
{
  const pills = [makePill('help')];
  const tabs = [makeTab('help', 'page5')];
  const setActiveCalls = [];
  const sandbox = {
    document: makeFakeDoc(pills, tabs),
    localStorage: { getItem: () => null, setItem: () => {} },
    console,
    Map, Set, Array,
  };
  vm.createContext(sandbox);
  vm.runInContext('function _setActiveStage(stage) { __setActiveCalls.push(stage); }', sandbox);
  vm.runInContext('var __setActiveCalls = [];', sandbox);
  vm.runInContext(fnInit, sandbox);
  vm.runInContext('_initTabGroupPills();', sandbox);
  pills[0].click();
  ok('child tab .click() fired exactly once',
     tabs[0]._clicks === 1, 'got ' + tabs[0]._clicks);
  ok('_setActiveStage NOT called (single-child path skips it)',
     sandbox.__setActiveCalls.length === 0,
     'got ' + JSON.stringify(sandbox.__setActiveCalls));
}

// --- Test: clicking a multi-child pill still calls _setActiveStage
console.log('\nTest: multi-child pill click → _setActiveStage');
{
  const pills = [makePill('discovery')];
  const tabs = [
    makeTab('discovery', 'page1'),
    makeTab('discovery', 'page12'),
  ];
  const sandbox = {
    document: makeFakeDoc(pills, tabs),
    localStorage: { getItem: () => null, setItem: () => {} },
    console,
    Map, Set, Array,
  };
  vm.createContext(sandbox);
  vm.runInContext('var __setActiveCalls = []; function _setActiveStage(stage) { __setActiveCalls.push(stage); }', sandbox);
  vm.runInContext(fnInit, sandbox);
  vm.runInContext('_initTabGroupPills();', sandbox);
  pills[0].click();
  ok('_setActiveStage called with "discovery"',
     sandbox.__setActiveCalls.length === 1 && sandbox.__setActiveCalls[0] === 'discovery');
  ok('child tabs were NOT auto-clicked (multi-child stage)',
     tabs[0]._clicks === 0 && tabs[1]._clicks === 0);
}

// =============================================================================
// Confirm against the actual atlas markup that 'help' is exactly single-child
// =============================================================================
console.log('\n=== Atlas markup sanity ===');

// Count <button data-page=...data-stage="help"> in the actual HTML
const helpPageBtnRegex = /<button[^>]*data-page="page\d+"[^>]*data-stage="help"/g;
const helpMatches = html.match(helpPageBtnRegex) || [];
ok('exactly 1 page button with data-stage="help" exists in markup',
   helpMatches.length === 1, 'got ' + helpMatches.length);

const helpPillRegex = /<button[^>]*class="tab-stage-pill"[^>]*data-stage="help"/g;
const helpPillMatches = html.match(helpPillRegex) || [];
ok('exactly 1 tab-stage-pill with data-stage="help" exists',
   helpPillMatches.length === 1);

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
