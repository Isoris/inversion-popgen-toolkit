// =============================================================================
// turn 147 — L3 toolbar disclosure (consolidate 3-row sprawl → 1 row)
// =============================================================================
// Quentin (turn 144 screenshot 1): "L3 contingency-table toolbar sprawls
// across 3 rows on page 1, ~120px wasted vertical."
//
// The fix:
//   - 14 secondary controls get class `l3-more-item`
//   - #l3Panel ships with class `l3-more-collapsed` so they're hidden by default
//   - A new ⋯ more / ⋯ less button (#l3MoreToggleBtn) toggles the class
//   - State persists to localStorage 'pca_scrubber_v3.l3MoreExpanded'
//   - Functionality preserved — hidden items keep their listeners, just
//     not visible until the user opts in
//   - The redundant manual line break (.l3-head-break) was removed since
//     flex-wrap: wrap on the parent handles narrow viewports naturally
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
// 1. CSS rules
// ============================================================================
console.log('\n=== 1. CSS rules ===');
ok('CSS rule for .l3-more-collapsed .l3-more-item present',
   /#l3Panel\.l3-more-collapsed\s+\.l3-more-item\s*\{[\s\S]{0,80}?display:\s*none\s*!important/.test(html));
ok('CSS rule uses !important to beat inline display style',
   /#l3Panel\.l3-more-collapsed\s+\.l3-more-item[\s\S]{0,80}?display:\s*none\s*!important/.test(html));
ok('Toggle button has its own CSS rule',
   /#l3MoreToggleBtn\s*\{[\s\S]{0,300}?cursor:\s*pointer/.test(html));
ok('Toggle gets accent style when expanded',
   /#l3Panel:not\(\.l3-more-collapsed\)\s+#l3MoreToggleBtn[\s\S]{0,160}?background:\s*var\(--accent\)/.test(html));

// ============================================================================
// 2. HTML markers
// ============================================================================
console.log('\n=== 2. HTML markers ===');
ok('#l3Panel ships with l3-more-collapsed class by default',
   /<div class="panel l3-more-collapsed" id="l3Panel">/.test(html));
ok('#l3MoreToggleBtn declared in markup',
   /id="l3MoreToggleBtn"/.test(html));
ok('toggle button text starts as "⋯ more"',
   /<button id="l3MoreToggleBtn"[\s\S]{0,400}?>⋯ more<\/button>/.test(html));
ok('toggle has descriptive title',
   /id="l3MoreToggleBtn"[\s\S]{0,400}?title="Show \/ hide secondary L3 controls/.test(html));

// ============================================================================
// 3. The 14 secondary controls each carry l3-more-item
// ============================================================================
console.log('\n=== 3. Secondary controls marked ===');
const secondaryIds = [
  'l3ColorMode',
  'l3KMode',
  'l3HetToggleLabel',
  'l2SweepToggleLabel',
  'l2SweepInspectBtn',
  'gPanelOpenBtn',
  'l3BcScope',
  'l3ScaleMode',
  'l3ScalePanes',
  'l3DetailedBtn',
  'l3SpotlightTrackedBtn',
  'l3SpotlightClearBtn',
  'pinL2Btn',
  'dhCursorToggleBtn',
];
for (const id of secondaryIds) {
  // Match the element's class attribute and confirm l3-more-item is in it.
  // Tolerant of attribute order: search id="…" within ~400 chars before/after
  // a class= containing l3-more-item, OR the reverse.
  const re = new RegExp(
    'id="' + id + '"[\\s\\S]{0,1200}?class="[^"]*\\bl3-more-item\\b' +
    '|' +
    'class="[^"]*\\bl3-more-item\\b[^"]*"[\\s\\S]{0,1200}?id="' + id + '"');
  ok('#' + id + ' has l3-more-item class', re.test(html));
}

// ============================================================================
// 4. Workhorses are NOT marked (stay always-visible)
// ============================================================================
console.log('\n=== 4. Workhorse controls stay visible ===');
const workhorseIds = [
  'l3Layout',           // Focal/+L/+R/+L/R/±2/Dual
  'l3CompareUnit',      // L2/1w/5w/10w/Nw + N input  (Quentin's screenshot 3 priority)
  'l3ReclusterMode',    // recluster + draft action row + edit row
  'l3PromoteCandBtn',   // ★ promote
  'l3CollapseBtn',      // ▼ panel collapse
];
for (const id of workhorseIds) {
  // Match the element's class attribute. Class must NOT contain l3-more-item.
  // We pull a tight window around the id="…" attribute and check the same tag's class.
  const tagRe = new RegExp('<[^>]+id="' + id + '"[^>]*>');
  const m = html.match(tagRe);
  if (!m) { ok('#' + id + ' tag found',  false, 'tag not found'); continue; }
  ok('#' + id + ' has NO l3-more-item class', !/class="[^"]*\bl3-more-item\b/.test(m[0]));
}

// ============================================================================
// 5. Manual line break removed (let flex-wrap do its job)
// ============================================================================
console.log('\n=== 5. Redundant line-break removed ===');
ok('.l3-head-break div is gone from DOM',
   !/<div class="l3-head-break"/.test(html));
ok('explanatory comment about removal is present',
   /turn 147[\s\S]{0,200}?manual line-break removed/.test(html));

// ============================================================================
// 6. JS wiring — _applyL3MoreState + click handler + persistence
// ============================================================================
console.log('\n=== 6. JS wiring ===');
ok('_applyL3MoreState function defined',
   /function _applyL3MoreState\(expanded\)/.test(html));
ok('toggles the l3-more-collapsed class',
   /classList\.toggle\('l3-more-collapsed',\s*!expanded\)/.test(html));
ok('button text flips between "⋯ more" / "⋯ less"',
   /textContent\s*=\s*expanded\s*\?\s*'⋯ less'\s*:\s*'⋯ more'/.test(html));
ok('reads localStorage on boot',
   /localStorage\.getItem\('pca_scrubber_v3\.l3MoreExpanded'\)/.test(html));
ok('persists choice on click',
   /localStorage\.setItem\('pca_scrubber_v3\.l3MoreExpanded'/.test(html));
ok('applies state at startup',
   /_applyL3MoreState\(_l3MoreExpanded\)/.test(html));
ok('click handler bound to the toggle',
   /_moreBtn\.addEventListener\('click'/.test(html));

// ============================================================================
// 7. Sandboxed JS behavior — round-trip through localStorage
// ============================================================================
console.log('\n=== 7. Sandboxed behavior ===');
{
  // Extract just the IIFE body that wires the toggle.
  // We need _applyL3MoreState + the read-from-storage + click wiring.
  // Easier: extract the entire IIFE that wraps both the resize handle
  // and our toggle wiring (pre-existing function, plus our additions).
  // We then drive it with mocked DOM.

  // Build minimal mock DOM
  const _ls = new Map();
  const panelEl = {
    id: 'l3Panel',
    _classes: new Set(['panel', 'l3-more-collapsed']),
    classList: {
      toggle(name, force) {
        if (force === true) panelEl._classes.add(name);
        else if (force === false) panelEl._classes.delete(name);
        else if (panelEl._classes.has(name)) panelEl._classes.delete(name);
        else panelEl._classes.add(name);
      },
      contains(name) { return panelEl._classes.has(name); },
      add(name) { panelEl._classes.add(name); },
      remove(name) { panelEl._classes.delete(name); },
    },
  };
  const moreBtn = {
    id: 'l3MoreToggleBtn',
    textContent: '',
    title: '',
    _listeners: { click: [] },
    addEventListener(type, fn) {
      if (!this._listeners[type]) this._listeners[type] = [];
      this._listeners[type].push(fn);
    },
    _click() { (this._listeners.click || []).forEach(fn => fn({})); },
  };

  // Stub out everything else _applyL3MoreState wiring touches
  const elements = new Map([
    ['l3Panel', panelEl],
    ['l3MoreToggleBtn', moreBtn],
    ['l3Resize', null],
    ['l3CollapseBtn', null],
    ['l3TitleClickable', null],
  ]);

  // Reconstruct the toggle wiring exactly as it appears in the IIFE.
  // The block lives between the "turn 147" comment and the "})();" IIFE
  // close. Anchor to that close so the brace count is balanced.
  const re2 = /\/\/ turn 147 — ⋯ more disclosure for the L3 toolbar[\s\S]*?_applyL3MoreState\(_l3MoreExpanded\);\s*\n\s*\}\);\s*\n\s*\}/;
  const m = html.match(re2);
  if (!m) { ok('extract turn-147 wiring', false, 'not found'); }
  else {
    const sandbox = {
      document: { getElementById: (id) => elements.get(id) || null },
      localStorage: {
        getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
        setItem: (k, v) => _ls.set(k, String(v)),
      },
      console,
    };
    const ctx = vm.createContext(sandbox);
    vm.runInContext(m[0], ctx);

    ok('bootstrapped collapsed (no localStorage entry)',
       panelEl._classes.has('l3-more-collapsed'));
    ok('button text starts at "⋯ more"',
       moreBtn.textContent === '⋯ more');

    // Click → expand
    moreBtn._click();
    ok('after click: l3-more-collapsed removed',
       !panelEl._classes.has('l3-more-collapsed'));
    ok('after click: button text → "⋯ less"',
       moreBtn.textContent === '⋯ less');
    ok('after click: localStorage saved "1"',
       _ls.get('pca_scrubber_v3.l3MoreExpanded') === '1');

    // Click again → collapse
    moreBtn._click();
    ok('second click: l3-more-collapsed re-added',
       panelEl._classes.has('l3-more-collapsed'));
    ok('second click: button text → "⋯ more"',
       moreBtn.textContent === '⋯ more');
    ok('second click: localStorage saved "0"',
       _ls.get('pca_scrubber_v3.l3MoreExpanded') === '0');

    // Bootstrap with a previously-expanded value: simulate by re-running with
    // localStorage already containing "1"
    const _ls2 = new Map([['pca_scrubber_v3.l3MoreExpanded', '1']]);
    const panelEl2 = {
      id: 'l3Panel',
      _classes: new Set(['panel', 'l3-more-collapsed']),  // ships collapsed
      classList: {
        toggle(name, force) {
          if (force === true) panelEl2._classes.add(name);
          else if (force === false) panelEl2._classes.delete(name);
          else if (panelEl2._classes.has(name)) panelEl2._classes.delete(name);
          else panelEl2._classes.add(name);
        },
        contains(name) { return panelEl2._classes.has(name); },
        add(name) { panelEl2._classes.add(name); },
        remove(name) { panelEl2._classes.delete(name); },
      },
    };
    const moreBtn2 = {
      id: 'l3MoreToggleBtn', textContent: '', title: '',
      _listeners: { click: [] },
      addEventListener(t, fn) { (this._listeners[t] = this._listeners[t] || []).push(fn); },
    };
    const elements2 = new Map([
      ['l3Panel', panelEl2],
      ['l3MoreToggleBtn', moreBtn2],
      ['l3Resize', null], ['l3CollapseBtn', null], ['l3TitleClickable', null],
    ]);
    const sandbox2 = {
      document: { getElementById: (id) => elements2.get(id) || null },
      localStorage: { getItem: (k) => _ls2.has(k) ? _ls2.get(k) : null,
                       setItem: (k, v) => _ls2.set(k, String(v)) },
      console,
    };
    const ctx2 = vm.createContext(sandbox2);
    vm.runInContext(m[0], ctx2);
    ok('boot with saved="1": panel becomes expanded',
       !panelEl2._classes.has('l3-more-collapsed'));
    ok('boot with saved="1": button text starts "⋯ less"',
       moreBtn2.textContent === '⋯ less');
  }
}

// ============================================================================
// 8. Functional preservation — hidden items keep their IDs and listeners
// ============================================================================
console.log('\n=== 8. Hidden items still functional ===');
// All 14 element IDs still appear exactly once each in the markup
for (const id of secondaryIds) {
  const matches = html.match(new RegExp('id="' + id + '"', 'g')) || [];
  ok('#' + id + ' present exactly once in DOM', matches.length === 1);
}
// And the addEventListener wiring for each is preserved (they're outside the
// HTML toolbar block, in the JS section that wires existing buttons by ID)
ok('l3DetailedBtn click wiring still present',
   /document\.getElementById\('l3DetailedBtn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));
ok('l3SpotlightTrackedBtn click wiring still present',
   /document\.getElementById\('l3SpotlightTrackedBtn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));
ok('l3SpotlightClearBtn click wiring still present',
   /document\.getElementById\('l3SpotlightClearBtn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));
ok('pinL2Btn click wiring still present',
   /document\.getElementById\('pinL2Btn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));
ok('dhCursorToggleBtn click wiring still present',
   /document\.getElementById\('dhCursorToggleBtn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));
ok('gPanelOpenBtn click wiring still present',
   /document\.getElementById\('gPanelOpenBtn'\)[\s\S]{0,400}?addEventListener\('click'/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
