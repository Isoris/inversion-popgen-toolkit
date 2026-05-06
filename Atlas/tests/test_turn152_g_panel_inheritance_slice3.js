// =============================================================================
// turn 152 — G-panel Slice 3: inheritance tab
// =============================================================================
// Replaces the Slice 1 placeholder body of the inheritance tab with a real
// renderer + interactive controls. Compute is unchanged (turn 115's
// runInheritanceCompute / inheritanceGroupClustering); this slice wires
// the result into the G-panel as one of three group-flavor surfaces.
//
// What this turn adds:
//   - state.gPanelInheritanceThreshold (default 0.15, persists to localStorage)
//   - _gpInhEnsureThreshold()  — load + clamp + write back
//   - _gpInhSetThreshold(t)    — clamp + persist + invalidate cache
//   - _gpInhSiToCga(si)        — sample-index → CGA helper
//   - _gpInhMembersForGroup(result, gid) — Int32Array of fish indices
//   - _gpInhMakeManualGroup(result, gid) — wraps addToManualGroup
//   - _gpInhItemLabel(result, idx)       — "I1" / "I2" formatter
//   - _gpInhRenderEmpty(reason, items)   — empty-state body
//   - _gpInhRenderResult(result, items, threshold)
//   - _gPanelRenderTabInheritance() — top-level tab body, replaces placeholder
//   - Modal interaction wiring: slider (debounced), recompute, matrix, mkmg
//
// Header subtitle + tab strip "(slice N)" tags removed since all three
// slices now ship.
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
// 1. Source-pattern: state slot + helpers declared
// ============================================================================
console.log('\n=== 1. Helpers declared ===');

ok('state.gPanelInheritanceThreshold default = 0.15',
   /gPanelInheritanceThreshold:\s*0\.15/.test(html));

ok('_GP_INH_THRESHOLD_KEY constant declared',
   /const\s+_GP_INH_THRESHOLD_KEY\s*=\s*'inversion_atlas\.gPanelInheritanceThreshold'/.test(html));

ok('_GP_INH_THRESHOLD_MIN = 0.05',
   /const\s+_GP_INH_THRESHOLD_MIN\s*=\s*0\.05/.test(html));

ok('_GP_INH_THRESHOLD_MAX = 0.40',
   /const\s+_GP_INH_THRESHOLD_MAX\s*=\s*0\.40/.test(html));

ok('_gpInhEnsureThreshold() function declared',
   /function\s+_gpInhEnsureThreshold\s*\(\s*\)/.test(html));

ok('_gpInhSetThreshold(t) function declared',
   /function\s+_gpInhSetThreshold\s*\(\s*t\s*\)/.test(html));

ok('_gpInhSiToCga(si) function declared',
   /function\s+_gpInhSiToCga\s*\(\s*si\s*\)/.test(html));

ok('_gpInhMembersForGroup(result, groupId) function declared',
   /function\s+_gpInhMembersForGroup\s*\(\s*result\s*,\s*groupId\s*\)/.test(html));

ok('_gpInhMakeManualGroup(result, groupId) function declared',
   /function\s+_gpInhMakeManualGroup\s*\(\s*result\s*,\s*groupId\s*\)/.test(html));

ok('_gpInhItemLabel(result, itemIdx) function declared',
   /function\s+_gpInhItemLabel\s*\(\s*result\s*,\s*itemIdx\s*\)/.test(html));

ok('_gpInhRenderEmpty(reason, items) function declared',
   /function\s+_gpInhRenderEmpty\s*\(\s*reason\s*,\s*items\s*\)/.test(html));

ok('_gpInhRenderResult(result, items, threshold) function declared',
   /function\s+_gpInhRenderResult\s*\(\s*result\s*,\s*items\s*,\s*threshold\s*\)/.test(html));

ok('all 6 inheritance helpers exposed on window',
   /window\._gpInhEnsureThreshold\b/.test(html) &&
   /window\._gpInhSetThreshold\b/.test(html) &&
   /window\._gpInhMembersForGroup\b/.test(html) &&
   /window\._gpInhMakeManualGroup\b/.test(html) &&
   /window\._gpInhItemLabel\b/.test(html) &&
   /window\._gpInhRenderEmpty\b/.test(html));

// ============================================================================
// 2. Source-pattern: tab body uses live compute, not placeholder text
// ============================================================================
console.log('\n=== 2. Inheritance tab body — placeholder removed, real body shipped ===');

const inhTabRe = /function\s+_gPanelRenderTabInheritance\s*\(\s*\)\s*\{([\s\S]*?)\nfunction\s/;
const inhTabMatch = html.match(inhTabRe);
ok('_gPanelRenderTabInheritance body extracted', !!inhTabMatch);
const inhTabBody = inhTabMatch ? inhTabMatch[1] : '';

ok('placeholder "Slice 3 pending" text REMOVED',
   !/Slice 3 pending/.test(inhTabBody),
   'placeholder copy should no longer be present');

ok('placeholder "the inheritance compute would produce a result" REMOVED',
   !/the inheritance compute would produce a result/.test(inhTabBody));

ok('tab body calls _gatherActiveCandidatesForInheritance',
   /_gatherActiveCandidatesForInheritance\s*\(/.test(inhTabBody));

ok('tab body calls _gpInhEnsureThreshold',
   /_gpInhEnsureThreshold\s*\(/.test(inhTabBody));

ok('tab body renders threshold slider with correct attributes',
   /id="gpInhThresholdSlider"/.test(inhTabBody) &&
   /min="'\s*\+\s*_GP_INH_THRESHOLD_MIN/.test(inhTabBody) &&
   /max="'\s*\+\s*_GP_INH_THRESHOLD_MAX/.test(inhTabBody));

ok('tab body renders recompute button (id=gpInhComputeBtn)',
   /id="gpInhComputeBtn"/.test(inhTabBody));

ok('tab body renders matrix-view shortcut (id=gpInhMatrixBtn)',
   /id="gpInhMatrixBtn"/.test(inhTabBody));

ok('tab body falls through to _gpInhRenderResult when state.inheritanceResult exists',
   /_gpInhRenderResult\s*\(\s*result\s*,\s*items\s*,\s*threshold\s*\)/.test(inhTabBody));

ok('tab body falls through to _gpInhRenderEmpty when items.length < 2',
   /_gpInhRenderEmpty\s*\(/.test(inhTabBody));

// ============================================================================
// 3. Source-pattern: result renderer surfaces per-group cards
// ============================================================================
console.log('\n=== 3. Result renderer — per-group cards + actions ===');

const inhResRe = /function\s+_gpInhRenderResult\s*\([^)]*\)\s*\{([\s\S]*?)\nfunction\s/;
const inhResMatch = html.match(inhResRe);
ok('_gpInhRenderResult body extracted', !!inhResMatch);
const inhResBody = inhResMatch ? inhResMatch[1] : '';

ok('renders per-candidate group-counts summary',
   /per_item_n_groups/.test(inhResBody));

ok('renders per-group make-manual button with data-gpinh-mkmg-gid',
   /data-gpinh-mkmg-gid/.test(inhResBody));

ok('uses _gpInhMembersForGroup to compute n_fish per group',
   /_gpInhMembersForGroup\s*\(\s*result\s*,\s*gid\s*\)/.test(inhResBody));

ok('renders n_fish + n_bands per group card',
   /n_fish=/.test(inhResBody) && /n_bands=/.test(inhResBody));

ok('handles empty group list with helpful message',
   /Compute returned 0 inheritance groups/.test(inhResBody));

// ============================================================================
// 4. Source-pattern: modal wiring
// ============================================================================
console.log('\n=== 4. Modal wiring — slider, recompute, matrix, mkmg ===');

const modalRe = /function\s+_renderGPanelModal\s*\(\s*\)\s*\{([\s\S]*?)\nfunction\s+_gPanelClose/;
const modalMatch = html.match(modalRe);
ok('_renderGPanelModal body extracted', !!modalMatch);
const modalBody = modalMatch ? modalMatch[1] : '';

ok('modal has "if (activeTab === \'inheritance\')" branch',
   /if\s*\(\s*activeTab\s*===\s*'inheritance'\s*\)/.test(modalBody));

ok('slider input handler debounces with setTimeout',
   /'gpInhThresholdSlider'/.test(modalBody) &&
   /setTimeout/.test(modalBody));

ok('slider live-updates value label before debounce',
   /gpInhThresholdValue/.test(modalBody) &&
   /textContent\s*=\s*v\.toFixed/.test(modalBody));

ok('slider debounce calls _gpInhSetThreshold + runInheritanceCompute',
   /_gpInhSetThreshold\s*\(\s*v\s*\)/.test(modalBody) &&
   /runInheritanceCompute\s*\(\s*\{[^}]*force:\s*true[^}]*threshold:\s*v[^}]*\}\s*\)/.test(modalBody));

ok('recompute button (gpInhComputeBtn) wired',
   /'gpInhComputeBtn'/.test(modalBody));

ok('matrix button (gpInhMatrixBtn) calls openInheritanceMatrix',
   /'gpInhMatrixBtn'/.test(modalBody) &&
   /openInheritanceMatrix\s*\(\s*\)/.test(modalBody));

ok('per-group "+ make manual" buttons wired via [data-gpinh-mkmg-gid]',
   /querySelectorAll\s*\(\s*'\[data-gpinh-mkmg-gid\]'\s*\)/.test(modalBody));

ok('mkmg handler calls _gpInhMakeManualGroup with the parsed gid',
   /_gpInhMakeManualGroup\s*\(\s*result\s*,\s*gid\s*\)/.test(modalBody));

// ============================================================================
// 5. Source-pattern: stale Slice-N labels removed from header / tab strip
// ============================================================================
console.log('\n=== 5. Header + tab strip cleaned of stale Slice-N copy ===');

ok('header subtitle no longer says "Slice 1 ships manual tab"',
   !/Slice 1 ships manual tab; karyotype \(Slice 2\) and inheritance \(Slice 3\) are placeholders/.test(html));

ok('header subtitle now says "all three flavours wired"',
   /all three flavours wired/.test(html));

ok('tab strip no longer renders "(slice N)" tags',
   !/'\s*\+\s*sliceTag/.test(html),
   'sliceTag concat should be removed; the for-loop just emits t.label');

// The const _GPANEL_TABS still has slice numbers (unchanged contract for
// any caller that reads it), but they're no longer rendered.
const tabsConstMatch = html.match(/_GPANEL_TABS\s*=\s*\[[^\]]*\]/);
const tabsConstSrc = tabsConstMatch ? tabsConstMatch[0] : '';
ok('_GPANEL_TABS const still carries slice metadata for downstream readers',
   /slice:\s*1/.test(tabsConstSrc) &&
   /slice:\s*2/.test(tabsConstSrc) &&
   /slice:\s*3/.test(tabsConstSrc));

// ============================================================================
// 6. Sandboxed: helper unit tests
// ============================================================================
console.log('\n=== 6. Helper unit tests (sandboxed) ===');

function pullFunction(src, name) {
  const decl = `function ${name}`;
  const i = src.indexOf(decl);
  if (i < 0) return null;
  // Find first '{' after the opening paren
  let p = src.indexOf('(', i);
  if (p < 0) return null;
  let braceStart = src.indexOf('{', p);
  if (braceStart < 0) return null;
  let depth = 1, j = braceStart + 1;
  while (j < src.length && depth > 0) {
    const ch = src[j];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    j++;
  }
  return src.substring(i, j);
}

function pullConst(src, name) {
  const re = new RegExp('const\\s+' + name + '\\s*=[^;]*;', 'm');
  const m = src.match(re);
  return m ? m[0] : null;
}

function makeSandbox(stateOverrides, opts) {
  opts = opts || {};
  const ctx = {
    state: Object.assign({
      data: { n_samples: 6, samples: [
        { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
        { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
      ]},
      gPanelInheritanceThreshold: 0.15,
      manualGroups: [],
      inheritanceResult: null,
    }, stateOverrides || {}),
    Float32Array, Int32Array, Uint8Array, Map, Set, Number, console,
    window: opts.window || {},
    parseFloat, parseInt, isFinite, Math, Array, Object, JSON, String,
  };
  ctx.window.state = ctx.state;
  if (opts.localStorage) {
    ctx.localStorage = opts.localStorage;
  }
  vm.createContext(ctx);
  // Inject IGC default + minimum constants the helpers reference.
  vm.runInContext(
    'const _IGC_DEFAULT_COSINE_DIST_THRESHOLD = 0.15;\n' +
    'const _IGC_MIN_BANDS_FOR_CLUSTERING = 2;\n',
    ctx
  );
  // Pull threshold const block + helper functions.
  const consts = [
    pullConst(html, '_GP_INH_THRESHOLD_KEY'),
    pullConst(html, '_GP_INH_THRESHOLD_MIN'),
    pullConst(html, '_GP_INH_THRESHOLD_MAX'),
  ].filter(Boolean).join('\n');
  vm.runInContext(consts, ctx);

  const fns = [
    '_gpInhEnsureThreshold',
    '_gpInhSetThreshold',
    '_gpInhSiToCga',
    '_gpInhMembersForGroup',
    '_gpInhItemLabel',
  ];
  for (const n of fns) {
    const code = pullFunction(html, n);
    if (code) vm.runInContext(code, ctx);
  }
  return ctx;
}

// 6a. _gpInhEnsureThreshold defaults
{
  const sb = makeSandbox();
  const v = sb._gpInhEnsureThreshold();
  ok('default threshold = 0.15 when no localStorage',
     Math.abs(v - 0.15) < 1e-9);
}

// 6b. _gpInhEnsureThreshold reads localStorage
{
  const ls = {
    _store: { 'inversion_atlas.gPanelInheritanceThreshold': '0.20' },
    getItem(k) { return this._store[k]; },
    setItem(k, v) { this._store[k] = v; },
  };
  const sb = makeSandbox({ gPanelInheritanceThreshold: 0.15 }, { localStorage: ls });
  const v = sb._gpInhEnsureThreshold();
  ok('threshold loaded from localStorage = 0.20',
     Math.abs(v - 0.20) < 1e-9);
  ok('state slot updated to 0.20',
     Math.abs(sb.state.gPanelInheritanceThreshold - 0.20) < 1e-9);
}

// 6c. _gpInhEnsureThreshold ignores invalid localStorage values
{
  const ls = {
    _store: { 'inversion_atlas.gPanelInheritanceThreshold': '99' },
    getItem(k) { return this._store[k]; },
    setItem(k, v) { this._store[k] = v; },
  };
  const sb = makeSandbox({ gPanelInheritanceThreshold: 0.15 }, { localStorage: ls });
  const v = sb._gpInhEnsureThreshold();
  ok('out-of-range localStorage value falls back to state default',
     Math.abs(v - 0.15) < 1e-9);
}

// 6d. _gpInhSetThreshold clamps + persists
{
  const ls = {
    _store: {},
    getItem(k) { return this._store[k]; },
    setItem(k, v) { this._store[k] = v; },
  };
  const sb = makeSandbox({}, { localStorage: ls });
  // Stub invalidateInheritanceCache so the helper's optional invalidate path doesn't blow up
  vm.runInContext('var invalidateInheritanceCache = function(){};', sb);
  sb._gpInhSetThreshold(0.25);
  ok('set 0.25 → state slot is 0.25',
     Math.abs(sb.state.gPanelInheritanceThreshold - 0.25) < 1e-9);
  ok('set 0.25 → localStorage persists',
     ls._store['inversion_atlas.gPanelInheritanceThreshold'] === '0.25');
  // Above max → clamp
  sb._gpInhSetThreshold(0.99);
  ok('set 0.99 clamps to MAX (0.40)',
     Math.abs(sb.state.gPanelInheritanceThreshold - 0.40) < 1e-9);
  // Below min → clamp
  sb._gpInhSetThreshold(0.001);
  ok('set 0.001 clamps to MIN (0.05)',
     Math.abs(sb.state.gPanelInheritanceThreshold - 0.05) < 1e-9);
  // NaN → falls back to default
  sb._gpInhSetThreshold(NaN);
  ok('NaN falls back to default (0.15)',
     Math.abs(sb.state.gPanelInheritanceThreshold - 0.15) < 1e-9);
}

// 6e. _gpInhSetThreshold invalidates inheritance cache when available
{
  const sb = makeSandbox();
  let invalidatedTimes = 0;
  vm.runInContext(
    'var invalidateInheritanceCache = function() { _invCalls++; };',
    sb
  );
  sb._invCalls = 0;
  sb._gpInhSetThreshold(0.20);
  ok('invalidateInheritanceCache called on threshold change',
     sb._invCalls === 1);
}

// 6f. _gpInhSiToCga
{
  const sb = makeSandbox();
  ok('_gpInhSiToCga(0) returns "CGA001"',
     sb._gpInhSiToCga(0) === 'CGA001');
  ok('_gpInhSiToCga(99) returns null (out of range)',
     sb._gpInhSiToCga(99) === null);
}

// 6g. _gpInhItemLabel formats with seq_num
{
  const sb = makeSandbox();
  const result = {
    items_meta: [
      { id: 'A', seq_num: 1 },
      { id: 'B', seq_num: 2 },
      { id: 'C' },                       // no seq_num
    ],
  };
  ok('item 0 → "I1"', sb._gpInhItemLabel(result, 0) === 'I1');
  ok('item 1 → "I2"', sb._gpInhItemLabel(result, 1) === 'I2');
  ok('item 2 (no seq_num) → "C" (id fallback)',
     sb._gpInhItemLabel(result, 2) === 'C');
  ok('out-of-range item → "?"',
     sb._gpInhItemLabel(result, 99) === '?');
}

// 6h. _gpInhMembersForGroup unions fish across bands of one group
{
  const sb = makeSandbox();
  // Build a synthetic result: 2 items × 2 bands = 4 bands. Group g0 has
  // bands 0 and 2 (fish from item-0-band-0 and item-1-band-0). Group g1
  // has bands 1 and 3 (fish from item-0-band-1 and item-1-band-1).
  // 6 fish total. Item 0: fish 0,1,2 in band 0 ; fish 3,4,5 in band 1.
  // Item 1: fish 0,2,4 in band 0 ; fish 1,3,5 in band 1.
  // Group g0 union = {0,1,2,4} ; Group g1 union = {1,3,4,5}.
  const fishMasks = [
    new Uint8Array([1,1,1,0,0,0]),   // band 0 (item 0, band 0) → fish 0,1,2
    new Uint8Array([0,0,0,1,1,1]),   // band 1 (item 0, band 1) → fish 3,4,5
    new Uint8Array([1,0,1,0,1,0]),   // band 2 (item 1, band 0) → fish 0,2,4
    new Uint8Array([0,1,0,1,0,1]),   // band 3 (item 1, band 1) → fish 1,3,5
  ];
  const result = {
    n_bands_total: 4,
    fish_masks: fishMasks,
    cut: { group_id_per_band: new Int32Array([0, 1, 0, 1]) },
    rtab: {
      per_group: {
        0: { 0: [0], 1: [0] },   // any non-null is enough
        1: { 0: [1], 1: [1] },
      },
    },
  };
  const g0 = sb._gpInhMembersForGroup(result, 0);
  const g1 = sb._gpInhMembersForGroup(result, 1);
  ok('group 0 union = {0,1,2,4}',
     g0.length === 4 && g0[0] === 0 && g0[1] === 1 && g0[2] === 2 && g0[3] === 4);
  ok('group 1 union = {1,3,4,5}',
     g1.length === 4 && g1[0] === 1 && g1[1] === 3 && g1[2] === 4 && g1[3] === 5);
}

// 6i. _gpInhMembersForGroup handles missing/empty result
{
  const sb = makeSandbox();
  const empty = sb._gpInhMembersForGroup(null, 0);
  ok('null result → empty Int32Array',
     empty instanceof sb.Int32Array && empty.length === 0);
  const noBands = sb._gpInhMembersForGroup({ rtab: {}, n_bands_total: 0 }, 0);
  ok('empty result → empty Int32Array',
     noBands instanceof sb.Int32Array && noBands.length === 0);
}

// 6j. _gpInhMembersForGroup returns sorted indices
{
  const sb = makeSandbox();
  const fishMasks = [
    new Uint8Array([0,0,1,0,0,1]),   // fish 2, 5
    new Uint8Array([1,0,0,0,1,0]),   // fish 0, 4
  ];
  const result = {
    n_bands_total: 2,
    fish_masks: fishMasks,
    cut: { group_id_per_band: new Int32Array([7, 7]) },
    rtab: { per_group: { 7: { 0: [0], 1: [0] } } },
  };
  const g7 = sb._gpInhMembersForGroup(result, 7);
  ok('union {2,5,0,4} → sorted [0,2,4,5]',
     g7.length === 4 && g7[0] === 0 && g7[1] === 2 && g7[2] === 4 && g7[3] === 5);
}

// ============================================================================
// 7. Source-pattern: window exports
// ============================================================================
console.log('\n=== 7. Window exports for tests + downstream ===');

const exportRe = /window\._gpInh[A-Za-z]+\s*=\s*_gpInh[A-Za-z]+/g;
const exportMatches = html.match(exportRe) || [];
ok('at least 6 _gpInh* helpers on window',
   exportMatches.length >= 6,
   'found ' + exportMatches.length);

ok('threshold constants exposed for tests',
   /window\._GP_INH_THRESHOLD_KEY\b/.test(html) &&
   /window\._GP_INH_THRESHOLD_MIN\b/.test(html) &&
   /window\._GP_INH_THRESHOLD_MAX\b/.test(html));

// ============================================================================
// 8. Negative / regression: existing flow preserved
// ============================================================================
console.log('\n=== 8. Existing flow preserved ===');

ok('runInheritanceCompute still exists (unchanged)',
   /function\s+runInheritanceCompute\s*\(\s*opts\s*\)/.test(html));

ok('inheritanceGroupClustering still exists (unchanged)',
   /function\s+inheritanceGroupClustering\s*\(\s*items\s*,\s*opts\s*\)/.test(html));

ok('openInheritanceMatrix still exists (unchanged)',
   /function\s+openInheritanceMatrix\s*\(\s*opts\s*\)/.test(html));

ok('manual tab body still exists (Slice 1 unchanged)',
   /function\s+_gPanelRenderTabManual\s*\(\s*\)/.test(html));

ok('karyotype tab body still exists (Slice 2 unchanged)',
   /function\s+_gPanelRenderTabKaryotype\s*\(\s*\)/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
