// =============================================================================
// turn 128d integration test — compact tracked-samples panel parity with fixed
//
// Quentin's bug report: "The fixed layout for the settings of the tracked
// samples hasn't yet been applied to the compact layout tracked samples,
// which still has the old settings style."
//
// Fixed-mode #pcaTrackedAside has (since v4 turn 2 ask 2/3/4):
//   - Row A: trails | sign-align PC1 | lasso | scree
//   - Lasso confirm/clear bar (pcaLassoBar)
//   - Row B: trail back slider
//   - Row C: N tracked slider
//   - Row D: K(N) cycle button + "all" + 6 K-colored numeric band buttons
//   - Row E: Auto-pick + Clear (dark)
//
// Compact-mode #trackedSamplesPanelCompact previously had only:
//   - trails + sign-align PC1
//   - trail back / N tracked sliders
//   - 4-button band grid (all + 3 named bands) — NO K=6 support
//   - Auto-pick + Clear
//
// turn 128d: ported the fixed-mode layout to compact. New IDs use Compact
// suffix; mirrors are wired:
//   - pcaLassoToggleCompact / screeToggleCompact via _wireCompactCheckboxMirror
//   - kCycleBtnCompact via .click() delegation to kCycleBtnAside
//   - data-band-compact + data-k-band on the new 6 numeric buttons (auto-
//     coloured/enabled by _syncTrackedCompactUI's existing data-k-band loop)
//   - pcaLassoBarCompact mirrored from _updatePcaLassoUI
//   - pcaLassoConfirmBtnCompact / pcaLassoClearBtnCompact via .click() delegation
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
// Markup checks: new compact-mode controls exist in the right shape
// =============================================================================
console.log('\n=== Compact panel markup ===');

// Lasso checkbox + label
ok('pcaLassoToggleCompact <input type="checkbox"> exists',
   /<input[^>]*type="checkbox"[^>]*id="pcaLassoToggleCompact"/.test(html));
ok('pcaLassoLabelCompact wraps the checkbox',
   /<label[^>]*id="pcaLassoLabelCompact"[\s\S]{0,500}id="pcaLassoToggleCompact"[\s\S]{0,300}<\/label>/.test(html));

// Scree checkbox + label
ok('screeToggleCompact <input type="checkbox"> exists',
   /<input[^>]*type="checkbox"[^>]*id="screeToggleCompact"/.test(html));
ok('screeToggleLabelCompact wraps the checkbox',
   /<label[^>]*id="screeToggleLabelCompact"[\s\S]{0,500}id="screeToggleCompact"[\s\S]{0,300}<\/label>/.test(html));

// Lasso bar + badge + buttons (compact)
ok('pcaLassoBarCompact <div> exists with display:none initially',
   /<div[^>]*id="pcaLassoBarCompact"[^>]*display:\s*none/.test(html));
ok('pcaLassoBadgeCompact span exists',
   /<span[^>]*id="pcaLassoBadgeCompact"/.test(html));
ok('pcaLassoConfirmBtnCompact button exists',
   /<button[^>]*id="pcaLassoConfirmBtnCompact"/.test(html));
ok('pcaLassoClearBtnCompact button exists',
   /<button[^>]*id="pcaLassoClearBtnCompact"/.test(html));

// K-cycle button + 6 K-band buttons + "all"
ok('kCycleBtnCompact button exists with K=3 default label',
   /<button[^>]*id="kCycleBtnCompact"[\s\S]{0,1500}>K=3</.test(html));
ok('compact has data-band-compact="all" button (the new style)',
   /<button[^>]*data-band-compact="all"[^>]*class="active"[\s\S]{0,500}>all</.test(html));

// Six numeric buttons with both data-band-compact AND data-k-band attributes
for (let k = 0; k < 6; k++) {
  ok(`compact band button data-band-compact="${k}" data-k-band="${k}" present`,
     new RegExp(`<button[^>]*data-band-compact="${k}"[^>]*data-k-band="${k}"`).test(html));
}

// The OLD compact band buttons (band 1 (lo), band 2 (mid), band 3 (hi))
// must be GONE from the panel. They were inside #trackedSamplesPanelCompactBody.
// The auxiliary side-column #pcaBandPickColCompact still has them — that's
// a separate UI (and out of scope here).
const compactBodyBlock = (() => {
  const start = html.indexOf('id="trackedSamplesPanelCompactBody"');
  if (start < 0) return '';
  const end = html.indexOf('</div><!-- /#trackedSamplesPanelCompactBody -->', start);
  // Fall back to a generous slice if the close-comment isn't there
  const cap = html.substring(start, start + 20000);
  // Find the matching </div> by depth-counting? Easier: just slice 12000 chars
  // and trim — for the regex check below we just need to confirm that the
  // OLD labels aren't present. Use a 14000-char window which should cover
  // the body without leaking into the next panel.
  return cap.substring(0, 14000);
})();
ok('OLD compact body labels "band 1 (lo)" / "band 2 (mid)" / "band 3 (hi)" are gone',
   !/band 1 \(lo\)/.test(compactBodyBlock) &&
   !/band 2 \(mid\)/.test(compactBodyBlock) &&
   !/band 3 \(hi\)/.test(compactBodyBlock));

// =============================================================================
// JS wiring: mirror handlers
// =============================================================================
console.log('\n=== JS wiring ===');

ok('_wireCompactCheckboxMirror helper defined',
   /function\s+_wireCompactCheckboxMirror\s*\(/.test(html));

ok('_wireCompactCheckboxMirror called for pcaLassoToggleCompact ↔ pcaLassoToggle',
   /_wireCompactCheckboxMirror\('pcaLassoToggleCompact', 'pcaLassoToggle'\)/.test(html));

ok('_wireCompactCheckboxMirror called for screeToggleCompact ↔ screeToggle',
   /_wireCompactCheckboxMirror\('screeToggleCompact',\s*'screeToggle'\)/.test(html));

ok('kCycleBtnCompact click delegates to kCycleBtnAside',
   /kCycleBtnCompact[\s\S]{0,500}getElementById\('kCycleBtnAside'\)[\s\S]{0,200}\.click\(\)/.test(html));

ok('pcaLassoConfirmBtnCompact delegates to pcaLassoConfirmBtn',
   /_pcaLassoConfirmBtnCompact[\s\S]{0,300}getElementById\('pcaLassoConfirmBtn'\)[\s\S]{0,100}\.click\(\)/.test(html));

ok('pcaLassoClearBtnCompact delegates to pcaLassoClearBtn',
   /_pcaLassoClearBtnCompact[\s\S]{0,300}getElementById\('pcaLassoClearBtn'\)[\s\S]{0,100}\.click\(\)/.test(html));

ok('_updatePcaLassoUI mirrors onto pcaLassoToggleCompact',
   /_updatePcaLassoUI[\s\S]{0,1500}getElementById\('pcaLassoToggleCompact'\)/.test(html));

ok('_updatePcaLassoUI mirrors onto pcaLassoBarCompact',
   /_updatePcaLassoUI[\s\S]{0,1500}getElementById\('pcaLassoBarCompact'\)/.test(html));

ok('_syncTrackedCompactUI updates kCycleBtnCompact label',
   /_syncTrackedCompactUI[\s\S]{0,4500}getElementById\('kCycleBtnCompact'\)/.test(html));

// =============================================================================
// Behavioural tests — _wireCompactCheckboxMirror
// =============================================================================
console.log('\n=== _wireCompactCheckboxMirror behaviour ===');

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

const fnMirror = pullFunction(html, '_wireCompactCheckboxMirror');
ok('_wireCompactCheckboxMirror extractable', !!fnMirror);

function makeCheckbox() {
  const handlers = { change: [] };
  return {
    checked: false,
    addEventListener: (evt, cb) => { handlers[evt] = handlers[evt] || []; handlers[evt].push(cb); },
    dispatchEvent: (evt) => {
      const list = handlers[evt.type] || [];
      list.forEach(cb => cb(evt));
    },
    _handlers: handlers,
  };
}

function makeFakeDoc(map) {
  return {
    getElementById: (id) => map[id] || null,
  };
}

// Case 1: compact change → mirrors to fixed
console.log('\nTest: compact → fixed propagates');
{
  const compact = makeCheckbox();
  const fixed   = makeCheckbox();
  // Simulate fixed's external "change" listener mutating state
  let stateMut = false;
  fixed.addEventListener('change', () => { stateMut = true; });
  const sandbox = {
    document: makeFakeDoc({ pcaLassoToggleCompact: compact, pcaLassoToggle: fixed }),
    Event: function (type, opts) { return { type, bubbles: !!(opts && opts.bubbles) }; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnMirror, sandbox);
  vm.runInContext("_wireCompactCheckboxMirror('pcaLassoToggleCompact', 'pcaLassoToggle');", sandbox);
  // User toggles compact ON
  compact.checked = true;
  compact.dispatchEvent({ type: 'change' });
  ok('fixed.checked propagated to true',  fixed.checked === true);
  ok('fixed change handler fired (state mutation)', stateMut === true);
}

// Case 2: fixed change → mirrors to compact (reverse direction)
console.log('\nTest: fixed → compact propagates (reverse)');
{
  const compact = makeCheckbox();
  const fixed   = makeCheckbox();
  const sandbox = {
    document: makeFakeDoc({ pcaLassoToggleCompact: compact, pcaLassoToggle: fixed }),
    Event: function (type, opts) { return { type, bubbles: !!(opts && opts.bubbles) }; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnMirror, sandbox);
  vm.runInContext("_wireCompactCheckboxMirror('pcaLassoToggleCompact', 'pcaLassoToggle');", sandbox);
  // External code (sidebar) toggles fixed ON
  fixed.checked = true;
  fixed.dispatchEvent({ type: 'change' });
  ok('compact.checked propagated to true', compact.checked === true);
}

// Case 3: no infinite recursion when both fire
console.log('\nTest: no infinite recursion when checked is already in sync');
{
  const compact = makeCheckbox();
  const fixed   = makeCheckbox();
  let dispatchCount = 0;
  const _origDispatch = compact.dispatchEvent;
  compact.dispatchEvent = function (evt) {
    dispatchCount++;
    if (dispatchCount > 100) throw new Error('infinite recursion');
    return _origDispatch.call(this, evt);
  };
  const sandbox = {
    document: makeFakeDoc({ pcaLassoToggleCompact: compact, pcaLassoToggle: fixed }),
    Event: function (type, opts) { return { type, bubbles: !!(opts && opts.bubbles) }; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnMirror, sandbox);
  vm.runInContext("_wireCompactCheckboxMirror('pcaLassoToggleCompact', 'pcaLassoToggle');", sandbox);
  // Trigger change on compact
  compact.checked = true;
  let threw = false;
  try { compact.dispatchEvent({ type: 'change' }); } catch (_) { threw = true; }
  ok('no infinite recursion on compact change', !threw);
}

// Case 4: graceful when one side is missing
console.log('\nTest: graceful when one side is missing');
{
  const sandbox = {
    document: makeFakeDoc({ pcaLassoToggleCompact: null, pcaLassoToggle: null }),
    Event: function (type, opts) { return { type, bubbles: !!(opts && opts.bubbles) }; },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnMirror, sandbox);
  let threw = false;
  try {
    vm.runInContext("_wireCompactCheckboxMirror('pcaLassoToggleCompact', 'pcaLassoToggle');", sandbox);
  } catch (_) { threw = true; }
  ok('no throw when both sides missing', !threw);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
