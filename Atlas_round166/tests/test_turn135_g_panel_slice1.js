// =============================================================================
// turn 135 — Slice 1: G-panel unified groups popup (SPEC_g_panel_unified_groups.md)
// =============================================================================
// Implements:
//   - _renderGPanelModal()                — top-level modal renderer
//   - _gPanelClose() / _gPanelToggle()    — open/close/toggle helpers
//   - _gPanelRenderTabKaryotype()         — Slice 2 placeholder body
//   - _gPanelRenderTabInheritance()       — Slice 3 placeholder body
//   - _gPanelRenderTabManual()            — re-host of manual groups list
//   - _GPANEL_TABS                        — tab metadata array
//   - state.gPanelOpen / state.gPanelTab  — popup state
//   - DOM: #gPanelOverlay, #gPanelOpenBtn, #gPanelClose, #gPanelBody,
//          #manualGroupsListPopup (rendered into manual tab)
//   - Hotkey 'g' (lowercase, no modifiers, guards against input fields)
//   - renderManualGroupsList extended to ALSO populate
//     #manualGroupsListPopup when present
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
// 1. Source-level: function definitions + constants
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_renderGPanelModal defined',          /function _renderGPanelModal\b/.test(html));
ok('_gPanelClose defined',                /function _gPanelClose\b/.test(html));
ok('_gPanelToggle defined',               /function _gPanelToggle\b/.test(html));
ok('_gPanelRenderTabKaryotype defined',   /function _gPanelRenderTabKaryotype\b/.test(html));
ok('_gPanelRenderTabInheritance defined', /function _gPanelRenderTabInheritance\b/.test(html));
ok('_gPanelRenderTabManual defined',      /function _gPanelRenderTabManual\b/.test(html));

ok('_GPANEL_TABS array defined',          /const _GPANEL_TABS\s*=\s*\[/.test(html));
ok("_GPANEL_TABS includes 'karyotype'",   /key:\s*'karyotype'/.test(html));
ok("_GPANEL_TABS includes 'inheritance'", /key:\s*'inheritance'/.test(html));
ok("_GPANEL_TABS includes 'manual'",      /key:\s*'manual'/.test(html));

ok('window._renderGPanelModal exported',
   /window\._renderGPanelModal\b/.test(html));
ok('window._gPanelToggle exported',
   /window\._gPanelToggle\b/.test(html));

// ============================================================================
// 2. DOM elements: trigger button + overlay
// ============================================================================
console.log('\n=== 2. DOM elements ===');

ok('#gPanelOverlay div present',
   /id="gPanelOverlay"/.test(html));
ok('#gPanelOverlay starts display:none',
   /id="gPanelOverlay"[\s\S]{0,300}display:\s*none/.test(html));
ok('#gPanelOverlay sibling to other modal overlays (z-index 1000)',
   /id="gPanelOverlay"[\s\S]{0,400}z-index:\s*1000/.test(html));
ok('#gPanelOpenBtn trigger button present',
   /id="gPanelOpenBtn"/.test(html));
ok('trigger button labeled "G ▾"',
   /id="gPanelOpenBtn"[\s\S]{0,800}>G\s*▾</.test(html));
ok('trigger button sibling to L2-sweep inspect button',
   /l2SweepInspectBtn[\s\S]{0,2500}gPanelOpenBtn/.test(html));

// ============================================================================
// 3. state init: gPanelOpen / gPanelTab
// ============================================================================
console.log('\n=== 3. state init slots ===');

ok('state.gPanelOpen: false in init block',
   /gPanelOpen:\s*false/.test(html));
ok("state.gPanelTab: 'manual' in init block (default to working tab)",
   /gPanelTab:\s*'manual'/.test(html));

// ============================================================================
// 4. renderManualGroupsList extended to populate #manualGroupsListPopup
// ============================================================================
console.log('\n=== 4. renderManualGroupsList extension ===');

ok('renderManualGroupsList looks up #manualGroupsListPopup',
   /renderManualGroupsList[\s\S]{0,1500}getElementById\('manualGroupsListPopup'\)/.test(html));
ok('the popup container is added to the containers array',
   /containers\.push\(popup\)/.test(html));

// ============================================================================
// 5. Hotkey wiring — 'g' (lowercase, no modifiers, input-field guards)
// ============================================================================
console.log('\n=== 5. Hotkey wiring ===');

ok("_wireGPanelHotkey IIFE present",
   /_wireGPanelHotkey/.test(html));
ok('hotkey listens for key === "g"',
   /e\.key !== 'g'/.test(html));

// Pull the actual IIFE body (not the state-init comment that mentions the
// IIFE name) by anchoring on the function declaration form
const hotkeyMatch = html.match(/function _wireGPanelHotkey\b[\s\S]{0,2000}/);
ok('hotkey IIFE body extractable',
   hotkeyMatch !== null && hotkeyMatch[0].length > 200);

const hotkeyBody = hotkeyMatch ? hotkeyMatch[0] : '';
ok('hotkey rejects modifier keys',
   /altKey \|\| e\.ctrlKey \|\| e\.metaKey \|\| e\.shiftKey/.test(hotkeyBody));
ok('hotkey guards against INPUT/TEXTAREA/SELECT',
   /tag === 'INPUT' \|\| tag === 'TEXTAREA' \|\| tag === 'SELECT'/.test(hotkeyBody));
ok('hotkey guards against contentEditable',
   /isContentEditable/.test(hotkeyBody));
ok('hotkey calls preventDefault before _gPanelToggle',
   /e\.preventDefault\(\);[\s\S]{0,200}_gPanelToggle/.test(hotkeyBody));

// ============================================================================
// 6. Behavioral — extract & run G-panel block in sandbox
// ============================================================================
console.log('\n=== 6. Sandbox: tab body renderers ===');

function extractGPanelBlock() {
  const start = html.indexOf('// G-PANEL UNIFIED GROUPS MODAL');
  if (start < 0) return null;
  const endMarker = 'window._gPanelRenderTabManual     = _gPanelRenderTabManual;';
  const endIdx = html.indexOf(endMarker, start);
  if (endIdx < 0) return null;
  const after = html.indexOf('}', endIdx);
  if (after < 0) return null;
  return html.substring(start, after + 1);
}
const gpanelBlock = extractGPanelBlock();
ok('G-panel block extracts cleanly', gpanelBlock !== null && gpanelBlock.length > 1000);

function makeGPanelSandbox(scenario) {
  const stateStub = scenario.state || {
    candidate: null,
    candidateList: [],
    candidates: {},
    k: 3,
    activeMode: 'default',
    gPanelOpen: false,
    gPanelTab: 'manual',
  };
  // Fake DOM with overlay + body + buttons
  const elements = {};
  function makeEl(id, opts) {
    const el = {
      id,
      tagName: (opts && opts.tag) || 'DIV',
      innerHTML: '',
      style: {},
      dataset: {},
      classList: { add() {}, remove() {}, contains() { return false; } },
      _listeners: {},
      addEventListener(ev, fn) {
        if (!this._listeners[ev]) this._listeners[ev] = [];
        this._listeners[ev].push(fn);
      },
      removeEventListener() {},
      querySelectorAll() { return []; },
      querySelector() { return null; },
      click() {
        if (this._listeners['click']) {
          for (const fn of this._listeners['click']) fn({ target: this });
        }
      },
      getAttribute(name) { return this[name] || (this.dataset && this.dataset[name.replace('data-', '')]) || null; },
    };
    if (opts && opts.disabled !== undefined) el.disabled = opts.disabled;
    return el;
  }
  elements.gPanelOverlay = makeEl('gPanelOverlay');
  elements.gPanelOpenBtn = makeEl('gPanelOpenBtn');
  const docStub = {
    getElementById(id) { return elements[id] || null; },
    addEventListener() {},
    removeEventListener() {},
    body: { appendChild() {} },
    querySelector() { return null; },
    readyState: 'complete',
  };
  const ctx = {
    state: stateStub,
    window: undefined,
    document: docStub,
    console: { warn: () => {}, log: () => {} },
    renderManualGroupsList: scenario.renderManualGroupsList || function () {},
    Number, Array, Object, Math, JSON, Set, Date,
    isFinite, parseInt, parseFloat,
    setTimeout: () => {},
  };
  ctx.window = ctx;
  const context = vm.createContext(ctx);
  vm.runInContext(gpanelBlock, context);
  return { ctx: context, state: stateStub, elements };
}

// 6a. Karyotype tab body — empty-state (no candidate focused)
//
// turn 136: Slice 2 shipped. The empty-state still says "No candidate
// currently focused" but no longer mentions "Slice 2 pending".
{
  console.log('\n--- 6a. Karyotype tab body (no candidate focused) ---');
  const sb = makeGPanelSandbox({});
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('returns string', typeof body === 'string');
  ok('Slice 2 shipped — no "Slice 2 pending" placeholder',
     !/Slice 2 pending/.test(body));
  ok('shows no-candidate empty-state', /No candidate currently focused/.test(body));
  ok('mentions catalogue or page-1 pin', /catalogue|page-1 pin/.test(body));
}

// 6b. Karyotype tab body — focused candidate
//
// turn 136: Slice 2 shipped, replacing the "Slice 2 pending" placeholder.
// The candidate-info header is now a one-line breadcrumb, and the body
// shows real per-band rows (or an empty-state if locked_labels is missing).
// This sandbox test passes a candidate WITHOUT locked_labels so it sees
// the "no locked_labels" branch — that's the cleanest source-level
// guarantee that Slice 2 is in place. Full Slice 2 rendering is exercised
// by test_turn136_g_panel_karyotype_slice2.js.
{
  console.log('\n--- 6b. Karyotype tab body (with focused candidate, no labels) ---');
  const sb = makeGPanelSandbox({
    state: {
      candidate: { id: 'cand_42', start_bp: 14e6, end_bp: 16.5e6, K: 3 },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('mentions the candidate id', /cand_42/.test(body));
  ok('formats bp range as Mb', /14\.00–16\.50 Mb/.test(body));
  ok('shows K=3', /K=3/.test(body));
  ok('Slice 2 shipped — no "Slice 2 pending" placeholder',
     !/Slice 2 pending/.test(body));
  ok('shows no-labels empty-state for candidate without locked_labels',
     /no\s*<code>locked_labels<\/code>|locked_labels.*length 0|missing/.test(body));
}

// 6c. Inheritance tab body — turn 152 Slice 3 inverted these assertions:
// the Slice 1 placeholder is REMOVED (replaced by the real implementation)
// and the body now wires _gatherActiveCandidatesForInheritance + threshold
// slider + recompute / matrix / mkmg buttons. See test_turn152 for the
// new assertions.
{
  console.log('\n--- 6c. Inheritance tab body (Slice 3 — placeholder removed by turn 152) ---');
  // Sandbox the inheritance tab. Stub the helpers it now depends on.
  let sb = makeGPanelSandbox({});
  // Inject the IGC default constant (turn 152's _gpInhEnsureThreshold reads
  // it as a default-fallback when localStorage is absent).
  vm.runInContext('const _IGC_DEFAULT_COSINE_DIST_THRESHOLD = 0.15;', sb.ctx);
  // Stub _gatherActiveCandidatesForInheritance + _gpInhEnsureThreshold for
  // the no-candidates case. The real helpers live elsewhere; we only need
  // the tab body to render without ReferenceErrors.
  vm.runInContext(
    'var _gatherActiveCandidatesForInheritance = () => [];\n' +
    'var _gpInhEnsureThreshold = () => 0.15;\n',
    sb.ctx
  );
  let body = sb.ctx._gPanelRenderTabInheritance();
  ok('Slice 1 "pending" copy is GONE', !/Slice 3 pending/.test(body));
  ok('Slice 1 "confirmed candidates available" copy is GONE',
     !/Confirmed candidates available for inheritance compute/.test(body));
  ok('renders the new "Inheritance groups" header',
     /Inheritance groups/.test(body));
  ok('renders the empty-state when items.length < 2',
     /Need ≥2 candidates with locked labels/.test(body));
  ok('reports eligible candidate count (the new metric)',
     /eligible candidate/.test(body));

  // ≥2 eligible: should render the slider + buttons instead of the empty state
  sb = makeGPanelSandbox({
    state: {
      candidate: null, candidates: {}, k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'inheritance',
      candidateList: [
        { id: 'a', confirmed: true },
        { id: 'b', confirmed: true },
        { id: 'c', confirmed: true },
        { id: 'd', confirmed: false },
      ],
    },
  });
  vm.runInContext('const _IGC_DEFAULT_COSINE_DIST_THRESHOLD = 0.15;', sb.ctx);
  vm.runInContext(
    'var _gatherActiveCandidatesForInheritance = () => [' +
    '  { id: "a", labels: new Int8Array([0,1,2]), K: 3, start_bp: 0, end_bp: 1e6, seq_num: 1 },' +
    '  { id: "b", labels: new Int8Array([0,1,2]), K: 3, start_bp: 1e6, end_bp: 2e6, seq_num: 2 }' +
    '];\n' +
    'var _gpInhEnsureThreshold = () => 0.15;\n',
    sb.ctx
  );
  body = sb.ctx._gPanelRenderTabInheritance();
  ok('with ≥2 candidates, renders the threshold slider',
     /id="gpInhThresholdSlider"/.test(body));
  ok('with ≥2 candidates, renders the recompute button',
     /id="gpInhComputeBtn"/.test(body));
  ok('with ≥2 candidates, renders the matrix-view shortcut',
     /id="gpInhMatrixBtn"/.test(body));
}

// 6d. Manual tab body — embeds the popup container + all action buttons
{
  console.log('\n--- 6d. Manual tab body structure ---');
  const sb = makeGPanelSandbox({});
  const body = sb.ctx._gPanelRenderTabManual();
  ok('contains #manualGroupsListPopup container',
     /id="manualGroupsListPopup"/.test(body));
  ok('contains #gpMgAddBtn',          /id="gpMgAddBtn"/.test(body));
  ok('contains #gpMgExportBtn',       /id="gpMgExportBtn"/.test(body));
  ok('contains #gpMgImportBtn',       /id="gpMgImportBtn"/.test(body));
  ok('contains #gpMgBandPickBar',     /id="gpMgBandPickBar"/.test(body));
  ok('contains 6 K-band picker buttons',
     (body.match(/data-gpmgband="\d"/g) || []).length === 6);
  ok('K-band buttons range k0..k5',
     /data-gpmgband="0"/.test(body) && /data-gpmgband="5"/.test(body));
  ok('mentions "single source of truth" / sidebar parity',
     /sidebar list/.test(body));
}

// ============================================================================
// 7. Sandbox: _gPanelToggle + _gPanelClose state machine
// ============================================================================
console.log('\n=== 7. Toggle / close state machine ===');

{
  const sb = makeGPanelSandbox({});
  ok('state starts gPanelOpen=false', sb.state.gPanelOpen === false);
  // Toggle once → opens (renderModal sets gPanelOpen=true)
  sb.ctx._gPanelToggle();
  ok('first toggle opens the popup', sb.state.gPanelOpen === true);
  ok('overlay display set to flex',
     sb.elements.gPanelOverlay.style.display === 'flex');
  ok('overlay innerHTML populated (>1000 chars)',
     sb.elements.gPanelOverlay.innerHTML.length > 1000);
  // Toggle again → closes
  sb.ctx._gPanelToggle();
  ok('second toggle closes the popup', sb.state.gPanelOpen === false);
  ok('overlay display set to none after close',
     sb.elements.gPanelOverlay.style.display === 'none');
}

// ============================================================================
// 8. Sandbox: tab switching changes state.gPanelTab and re-renders
// ============================================================================
console.log('\n=== 8. Tab switching ===');

{
  // Render once with default tab (manual)
  const sb = makeGPanelSandbox({});
  sb.ctx._renderGPanelModal();
  ok("default state.gPanelTab = 'manual'", sb.state.gPanelTab === 'manual');
  ok("manual rendered (contains popup container)",
     /manualGroupsListPopup/.test(sb.elements.gPanelOverlay.innerHTML));

  // Switch to karyotype directly via state mutation + re-render
  sb.state.gPanelTab = 'karyotype';
  sb.ctx._renderGPanelModal();
  // turn 136: Slice 2 shipped — no more "Slice 2 pending" placeholder.
  // The karyotype body now starts with "Karyotype groups" header, then
  // (since this sandbox has no state.candidate) shows the no-candidate
  // empty-state.
  ok('karyotype rendered after state change (Slice 2 body)',
     /Karyotype groups/.test(sb.elements.gPanelOverlay.innerHTML) &&
     /No candidate currently focused/.test(sb.elements.gPanelOverlay.innerHTML));

  // Invalid tab values fall back to 'manual'
  sb.state.gPanelTab = 'bogus';
  sb.ctx._renderGPanelModal();
  ok("invalid tab falls back to 'manual'", sb.state.gPanelTab === 'manual');
}

// ============================================================================
// 9. Sandbox: rendered HTML contains expected structural elements
// ============================================================================
console.log('\n=== 9. Rendered HTML structure ===');

{
  const sb = makeGPanelSandbox({});
  sb.ctx._renderGPanelModal();
  const innerHTML = sb.elements.gPanelOverlay.innerHTML;
  ok('header contains "G · unified groups"', /G · unified groups/.test(innerHTML));
  ok('header has close button #gPanelClose',
     /id="gPanelClose"/.test(innerHTML));
  ok('header now mentions all-three-flavours-wired framing (turn 152 inverted)',
     /all three flavours wired/.test(innerHTML) &&
     !/Slice 1 ships manual tab/.test(innerHTML));
  ok('tab strip has 3 _gpTabBtn buttons',
     (innerHTML.match(/class="_gpTabBtn"/g) || []).length === 3);
  ok('tab buttons have data-gptab attributes',
     /data-gptab="karyotype"/.test(innerHTML) &&
     /data-gptab="inheritance"/.test(innerHTML) &&
     /data-gptab="manual"/.test(innerHTML));
  ok('body container has id="gPanelBody"',
     /id="gPanelBody"/.test(innerHTML));
}

// ============================================================================
// 10. renderManualGroupsList extension exercised in sandbox
// ============================================================================
console.log('\n=== 10. renderManualGroupsList populates popup container ===');

// Extract just renderManualGroupsList to validate its container collection
{
  const re = /function renderManualGroupsList\(\)[\s\S]*?\n\}/m;
  const match = html.match(re);
  ok('renderManualGroupsList source extractable', match !== null);
  if (match) {
    // Build a minimal sandbox with three containers and call the function
    const containers = {
      manualGroupsList:        { id: 'manualGroupsList', innerHTML: '' },
      manualGroupsListCompact: { id: 'manualGroupsListCompact', innerHTML: '' },
      manualGroupsListPopup:   { id: 'manualGroupsListPopup', innerHTML: '' },
    };
    const ctx = {
      state: { manualGroups: [
        { id: 'g1', name: 'wild_type', members: ['CGA001', 'CGA002'], color: '#ff0000', scope: 'chrom' },
      ] },
      document: {
        getElementById(id) { return containers[id] || null; },
      },
      escapeHtml: function (s) { return String(s); },
      window: undefined,
    };
    ctx.window = ctx;
    const context = vm.createContext(ctx);
    vm.runInContext(match[0], context);
    context.renderManualGroupsList();
    ok('sidebar container populated',
       containers.manualGroupsList.innerHTML.indexOf('wild_type') !== -1);
    ok('compact container populated',
       containers.manualGroupsListCompact.innerHTML.indexOf('wild_type') !== -1);
    ok('popup container populated (NEW this turn)',
       containers.manualGroupsListPopup.innerHTML.indexOf('wild_type') !== -1);

    // Without the popup container, only the original two are populated
    delete containers.manualGroupsListPopup;
    const c2 = { ...containers, manualGroupsListPopup: undefined };
    const ctx2 = {
      state: { manualGroups: [
        { id: 'g2', name: 'inversion', members: ['CGA003'], color: '#0000ff', scope: 'chrom' },
      ] },
      document: { getElementById(id) { return c2[id] || null; } },
      escapeHtml: function (s) { return String(s); },
      window: undefined,
    };
    ctx2.window = ctx2;
    const context2 = vm.createContext(ctx2);
    vm.runInContext(match[0], context2);
    context2.renderManualGroupsList();
    // The first two should now contain 'inversion' but no popup container exists
    ok('renderManualGroupsList no-ops on missing popup container (graceful)',
       c2.manualGroupsList.innerHTML.indexOf('inversion') !== -1);
  }
}

// ============================================================================
// 11. Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
