// =============================================================================
// turn 153 — Auto-register inheritance on candidate-list mutations
// =============================================================================
// Pre-turn-153 the inheritance compute was kept in sync with the candidate
// set by *draw-time* stale-cache guards inside _drawInheritanceLabelsStrip
// (and a sibling guard in the cockpit purity path). The G-panel inheritance
// tab and the matrix popup did NOT have these guards — they read
// state.inheritanceResult directly and rendered whatever was there.
//
// This turn moves invalidation upstream to the candidate-mutation funnel
// itself (persistCandidateList), so all consumers see a consistent picture
// without each having to implement its own staleness check.
//
// What this turn adds:
//   - _autoRegisterInheritanceOnCandidateChange()
//       Compares cached cache-key vs expected; if mismatch:
//         (a) invalidates state.inheritanceResult + state.inheritanceCacheKey
//         (b) schedules _notifyInheritanceConsumers via queueMicrotask
//             (debounced through _inheritanceNotifyScheduled)
//   - _notifyInheritanceConsumers()
//       Calls window.requestRepaint() (covers lines-strip),
//       re-renders G-panel modal IF gPanelOpen && gPanelTab === 'inheritance',
//       re-draws matrix popup IF its modal is in the DOM and visible.
//   - persistCandidateList() now calls _autoRegister... after _rebuild.
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
// 1. Source-pattern: helpers declared + exposed
// ============================================================================
console.log('\n=== 1. Helpers declared ===');

ok('_autoRegisterInheritanceOnCandidateChange function declared',
   /function\s+_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)/.test(html));

ok('_notifyInheritanceConsumers function declared',
   /function\s+_notifyInheritanceConsumers\s*\(\s*\)/.test(html));

ok('window._autoRegisterInheritanceOnCandidateChange exposed',
   /window\._autoRegisterInheritanceOnCandidateChange\s*=\s*_autoRegisterInheritanceOnCandidateChange/.test(html));

ok('window._notifyInheritanceConsumers exposed',
   /window\._notifyInheritanceConsumers\s*=\s*_notifyInheritanceConsumers/.test(html));

// ============================================================================
// 2. Source-pattern: persistCandidateList wires the hook
// ============================================================================
console.log('\n=== 2. persistCandidateList → auto-register hook ===');

const pclRe = /function\s+persistCandidateList\s*\(\s*\)\s*\{([\s\S]*?)\nfunction\s+loadCandidateList/;
const pclMatch = html.match(pclRe);
ok('persistCandidateList body extracted', !!pclMatch);
const pclBody = pclMatch ? pclMatch[1] : '';

ok('persistCandidateList still calls _rebuildCandidateRegistries (turn 129 contract preserved)',
   /_rebuildCandidateRegistries\s*\(\s*\)/.test(pclBody));

ok('persistCandidateList NOW calls _autoRegisterInheritanceOnCandidateChange',
   /_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)/.test(pclBody));

ok('auto-register hook is positioned AFTER _rebuildCandidateRegistries (registries fresh first)',
   pclBody.indexOf('_rebuildCandidateRegistries') < pclBody.indexOf('_autoRegisterInheritanceOnCandidateChange'));

ok('auto-register hook is positioned BEFORE localStorage write (deterministic order)',
   pclBody.indexOf('_autoRegisterInheritanceOnCandidateChange') <
   pclBody.indexOf('localStorage.setItem'));

ok('auto-register hook is guarded by typeof check',
   /typeof\s+_autoRegisterInheritanceOnCandidateChange\s*===\s*'function'/.test(pclBody));

ok('auto-register hook is wrapped in try/catch (non-blocking)',
   /try\s*\{\s*_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)\s*;?\s*\}\s*catch/.test(pclBody));

// ============================================================================
// 3. Source-pattern: _autoRegister body
// ============================================================================
console.log('\n=== 3. _autoRegister body ===');

const autoRe = /function\s+_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)\s*\{([\s\S]*?)\nfunction\s+_notifyInheritanceConsumers/;
const autoMatch = html.match(autoRe);
ok('_autoRegister body extracted', !!autoMatch);
const autoBody = autoMatch ? autoMatch[1] : '';

ok('reads candidates via _gatherActiveCandidatesForInheritance',
   /_gatherActiveCandidatesForInheritance\s*\(\s*\)/.test(autoBody));

ok('computes expected cache key via _inheritanceCacheKey',
   /_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*\)/.test(autoBody));

ok('compares against state.inheritanceCacheKey',
   /_state\.inheritanceCacheKey/.test(autoBody));

ok('no-ops on cache-key match (idempotent)',
   /cachedKey\s*===\s*expectedKey/.test(autoBody));

ok('calls invalidateInheritanceCache on mismatch',
   /invalidateInheritanceCache\s*\(\s*\)/.test(autoBody));

ok('debounces notify via _inheritanceNotifyScheduled flag',
   /_inheritanceNotifyScheduled/.test(autoBody));

ok('uses queueMicrotask when available (fastest debounce)',
   /queueMicrotask/.test(autoBody));

ok('falls back to setTimeout when queueMicrotask unavailable',
   /setTimeout\s*\(\s*fire\s*,\s*0\s*\)/.test(autoBody));

ok('guards against missing helpers (sandbox-safety)',
   /typeof\s+_gatherActiveCandidatesForInheritance\s*!==\s*'function'/.test(autoBody) &&
   /typeof\s+_inheritanceCacheKey\s*!==\s*'function'/.test(autoBody));

// ============================================================================
// 4. Source-pattern: _notifyInheritanceConsumers body
// ============================================================================
console.log('\n=== 4. _notifyInheritanceConsumers body ===');

const notifyRe = /function\s+_notifyInheritanceConsumers\s*\(\s*\)\s*\{([\s\S]*?)\n(?:window\.|function\s+)/;
const notifyMatch = html.match(notifyRe);
ok('_notifyInheritanceConsumers body extracted', !!notifyMatch);
const notifyBody = notifyMatch ? notifyMatch[1] : '';

ok('calls window.requestRepaint() when available (lines-strip path)',
   /window\.requestRepaint\s*\(\s*\)/.test(notifyBody));

ok('re-renders G-panel modal only when open AND on inheritance tab',
   /_state\.gPanelOpen/.test(notifyBody) &&
   /_state\.gPanelTab\s*===\s*'inheritance'/.test(notifyBody) &&
   /_renderGPanelModal\s*\(\s*\)/.test(notifyBody));

ok('re-draws matrix popup only when its modal is visible',
   /getElementById\s*\(\s*'inheritanceMatrixModal'\s*\)/.test(notifyBody) &&
   /matrixModal\.style\.display\s*!==\s*'none'/.test(notifyBody) &&
   /_redrawInheritanceMatrix\s*\(\s*\)/.test(notifyBody));

ok('each consumer is independently guarded (typeof + null checks)',
   (notifyBody.match(/typeof\s+\w[\w\.]*\s*===\s*'function'/g) || []).length >= 2 &&
   /matrixModal\s*&&/.test(notifyBody),
   '2+ typeof guards (requestRepaint, _renderGPanelModal, _redrawInheritanceMatrix) ' +
   'plus one null-check on the matrix modal');

ok('each consumer call is wrapped in try/catch (partial notify is OK)',
   (notifyBody.match(/try\s*\{[^}]*\}\s*catch/g) || []).length >= 3);

// ============================================================================
// 5. Sandboxed: _autoRegister behavioural tests
// ============================================================================
console.log('\n=== 5. _autoRegister sandboxed behaviour ===');

function pullFunction(src, name) {
  const decl = `function ${name}`;
  const i = src.indexOf(decl);
  if (i < 0) return null;
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

function makeSandbox(opts) {
  opts = opts || {};
  const ctx = {
    state: Object.assign({
      activeMode: 'default',
      inheritanceResult: opts.initialResult || null,
      inheritanceCacheKey: opts.initialKey || null,
      gPanelOpen: false,
      gPanelTab: 'manual',
    }, opts.state || {}),
    Float32Array, Int32Array, Uint8Array, Map, Set, Number, console,
    Object, Array, JSON, String, Math,
  };
  ctx.window = opts.window || {};
  ctx.window.state = ctx.state;
  ctx.queueMicrotask = opts.queueMicrotask || queueMicrotask;
  if (opts.document !== undefined) ctx.document = opts.document;
  vm.createContext(ctx);
  // Stub helpers the auto-register references.
  let gatherCallCount = 0;
  let cacheKeyCallCount = 0;
  ctx._gatherCallCount = () => gatherCallCount;
  ctx._cacheKeyCallCount = () => cacheKeyCallCount;
  ctx._mockGatherItems = opts.gatherItems || [];
  ctx._mockCacheKey = opts.expectedKey || 'k1';
  vm.runInContext(
    'var _gatherActiveCandidatesForInheritance = function() { ' +
    '  _gatherCallCount(); return _mockGatherItems; };\n' +
    'var _inheritanceCacheKey = function(items, mode) { ' +
    '  _cacheKeyCallCount(); return _mockCacheKey; };\n' +
    'var invalidateInheritanceCache = function() { ' +
    '  state.inheritanceResult = null; state.inheritanceCacheKey = null; };\n' +
    'var _renderGPanelModal = function() { ' +
    '  state._renderModalCalls = (state._renderModalCalls||0) + 1; };\n' +
    'var _redrawInheritanceMatrix = function() { ' +
    '  state._redrawMatrixCalls = (state._redrawMatrixCalls||0) + 1; };\n',
    ctx
  );
  // Inject the two real functions.
  const fnAuto = pullFunction(html, '_autoRegisterInheritanceOnCandidateChange');
  const fnNotify = pullFunction(html, '_notifyInheritanceConsumers');
  if (fnAuto) vm.runInContext(fnAuto, ctx);
  if (fnNotify) vm.runInContext(fnNotify, ctx);
  return ctx;
}

// 5a. No-op when cache key matches
{
  const sb = makeSandbox({
    initialKey: 'k1',
    initialResult: { foo: 'bar' },
    expectedKey: 'k1',
  });
  sb._autoRegisterInheritanceOnCandidateChange();
  ok('cache-key match → result preserved (no-op)',
     sb.state.inheritanceResult && sb.state.inheritanceResult.foo === 'bar');
  ok('cache-key match → cacheKey preserved',
     sb.state.inheritanceCacheKey === 'k1');
  ok('cache-key match → notify NOT scheduled',
     !sb.state._inheritanceNotifyScheduled);
}

// 5b. Invalidates + schedules notify on cache-key mismatch
{
  const sb = makeSandbox({
    initialKey: 'k_old',
    initialResult: { foo: 'bar' },
    expectedKey: 'k_new',
  });
  sb._autoRegisterInheritanceOnCandidateChange();
  ok('cache-key mismatch → result invalidated',
     sb.state.inheritanceResult === null);
  ok('cache-key mismatch → cacheKey invalidated',
     sb.state.inheritanceCacheKey === null);
  ok('cache-key mismatch → notify scheduled',
     sb.state._inheritanceNotifyScheduled === true);
}

// 5c. Debounce: multiple rapid calls coalesce to ONE notify
{
  let microtaskFires = 0;
  const queued = [];
  const sb = makeSandbox({
    initialKey: 'k_old',
    expectedKey: 'k_new',
    queueMicrotask: (fn) => { queued.push(fn); microtaskFires++; },
  });
  sb._autoRegisterInheritanceOnCandidateChange();
  sb._autoRegisterInheritanceOnCandidateChange();
  sb._autoRegisterInheritanceOnCandidateChange();
  ok('three rapid calls → only ONE microtask scheduled',
     microtaskFires === 1);
  // Drain the queued microtask to verify the flag resets
  for (const fn of queued) fn();
  ok('after microtask fires, _inheritanceNotifyScheduled resets to false',
     sb.state._inheritanceNotifyScheduled === false);
}

// 5d. After microtask resets, a new mismatch CAN re-schedule
{
  let microtaskFires = 0;
  const queued = [];
  const sb = makeSandbox({
    initialKey: 'k_old',
    expectedKey: 'k_new',
    queueMicrotask: (fn) => { queued.push(fn); microtaskFires++; },
  });
  sb._autoRegisterInheritanceOnCandidateChange();
  // Drain
  for (const fn of queued) fn();
  queued.length = 0;
  // Now the cache key is re-set inside auto-register? No — we never re-set it
  // in the function. Re-update mock key so next call invalidates again.
  sb.state.inheritanceCacheKey = 'k_new';   // simulate someone wrote new key
  sb._mockCacheKey = 'k_newer';
  sb._autoRegisterInheritanceOnCandidateChange();
  ok('after drain, second mismatch reschedules notify',
     microtaskFires === 2);
}

// 5e. Missing helpers → fail-soft (no throw, no state mutation)
{
  // Build a sandbox WITHOUT injecting the helper stubs.
  const ctx = {
    state: { inheritanceCacheKey: 'k1', inheritanceResult: { foo: 1 } },
    Object, Array, JSON, console, Math,
  };
  ctx.window = { state: ctx.state };
  vm.createContext(ctx);
  // Inject ONLY the auto-register, not the helpers it references.
  const fnAuto = pullFunction(html, '_autoRegisterInheritanceOnCandidateChange');
  vm.runInContext('var invalidateInheritanceCache = function() {};', ctx);
  vm.runInContext(fnAuto, ctx);
  let threw = null;
  try { ctx._autoRegisterInheritanceOnCandidateChange(); }
  catch (e) { threw = e; }
  ok('missing _gatherActiveCandidatesForInheritance → no throw',
     threw === null);
  ok('missing helpers → no state mutation',
     ctx.state.inheritanceResult && ctx.state.inheritanceResult.foo === 1);
}

// ============================================================================
// 6. Sandboxed: _notifyInheritanceConsumers behavioural tests
// ============================================================================
console.log('\n=== 6. _notifyInheritanceConsumers sandboxed behaviour ===');

// 6a. requestRepaint always called when present
{
  let repaintCalls = 0;
  const sb = makeSandbox({
    window: {
      requestRepaint: () => { repaintCalls++; },
    },
    // No DOM — matrix popup branch should fail-soft
    document: undefined,
  });
  sb._notifyInheritanceConsumers();
  ok('requestRepaint called once', repaintCalls === 1);
}

// 6b. G-panel modal re-rendered ONLY when open AND on inheritance tab
{
  // Closed
  let sb = makeSandbox({
    state: { gPanelOpen: false, gPanelTab: 'inheritance' },
  });
  sb._notifyInheritanceConsumers();
  ok('panel closed → modal NOT re-rendered',
     !sb.state._renderModalCalls);

  // Open but on manual tab
  sb = makeSandbox({
    state: { gPanelOpen: true, gPanelTab: 'manual' },
  });
  sb._notifyInheritanceConsumers();
  ok('panel open on manual tab → modal NOT re-rendered',
     !sb.state._renderModalCalls);

  // Open on inheritance tab
  sb = makeSandbox({
    state: { gPanelOpen: true, gPanelTab: 'inheritance' },
  });
  sb._notifyInheritanceConsumers();
  ok('panel open on inheritance tab → modal re-rendered exactly once',
     sb.state._renderModalCalls === 1);
}

// 6c. Matrix popup re-drawn only when visible
{
  // No matrix modal in DOM
  let sb = makeSandbox({
    document: { getElementById: () => null },
  });
  sb._notifyInheritanceConsumers();
  ok('no matrix modal in DOM → no redraw',
     !sb.state._redrawMatrixCalls);

  // Matrix modal in DOM but display:none
  sb = makeSandbox({
    document: {
      getElementById: (id) => id === 'inheritanceMatrixModal'
                              ? { style: { display: 'none' } } : null,
    },
  });
  sb._notifyInheritanceConsumers();
  ok('matrix modal hidden (display:none) → no redraw',
     !sb.state._redrawMatrixCalls);

  // Matrix modal in DOM and visible
  sb = makeSandbox({
    document: {
      getElementById: (id) => id === 'inheritanceMatrixModal'
                              ? { style: { display: 'flex' } } : null,
    },
  });
  sb._notifyInheritanceConsumers();
  ok('matrix modal visible → redraw called exactly once',
     sb.state._redrawMatrixCalls === 1);
}

// 6d. All three consumers fire when all three are present + active
{
  let repaintCalls = 0;
  const sb = makeSandbox({
    state: { gPanelOpen: true, gPanelTab: 'inheritance' },
    window: { requestRepaint: () => { repaintCalls++; } },
    document: {
      getElementById: (id) => id === 'inheritanceMatrixModal'
                              ? { style: { display: 'flex' } } : null,
    },
  });
  sb._notifyInheritanceConsumers();
  ok('all three consumers fire when all are active',
     repaintCalls === 1 &&
     sb.state._renderModalCalls === 1 &&
     sb.state._redrawMatrixCalls === 1);
}

// 6e. Partial notify: failure in one consumer does not block others
{
  let repaintCalls = 0;
  const sb = makeSandbox({
    state: { gPanelOpen: true, gPanelTab: 'inheritance' },
    window: { requestRepaint: () => { repaintCalls++; throw new Error('paint blew up'); } },
    document: {
      getElementById: (id) => id === 'inheritanceMatrixModal'
                              ? { style: { display: 'flex' } } : null,
    },
  });
  let threw = null;
  try { sb._notifyInheritanceConsumers(); } catch (e) { threw = e; }
  ok('thrown error in requestRepaint → no propagation',
     threw === null);
  ok('thrown error in requestRepaint → G-panel still re-rendered',
     sb.state._renderModalCalls === 1);
  ok('thrown error in requestRepaint → matrix still redrawn',
     sb.state._redrawMatrixCalls === 1);
}

// ============================================================================
// 7. Negative / regression: existing flow preserved
// ============================================================================
console.log('\n=== 7. Existing flow preserved ===');

ok('runInheritanceCompute still exists',
   /function\s+runInheritanceCompute\s*\(\s*opts\s*\)/.test(html));

ok('invalidateInheritanceCache still exists (just expanded around it)',
   /function\s+invalidateInheritanceCache\s*\(\s*\)/.test(html));

// Turn 155 extended _inheritanceCacheKey from (items, mode) to
// (items, mode, threshold) with threshold optional (defaults to
// state.gPanelInheritanceThreshold). The function still exists; the
// signature is now backwards-compatible with old (items, mode) callers.
ok('_inheritanceCacheKey still exists (turn 155: signature extended with optional threshold)',
   /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(html));

ok('_drawInheritanceLabelsStrip still has its draw-time stale-cache guard',
   /_state\.inheritanceCacheKey\s*!==\s*expectedKey/.test(html));

ok('_rebuildCandidateRegistries still exists (turn 129 contract)',
   /function\s+_rebuildCandidateRegistries\s*\(\s*\)/.test(html));

ok('persistCandidateList still writes localStorage',
   /localStorage\.setItem\s*\(\s*key\s*,\s*JSON\.stringify\s*\(\s*arr\s*\)/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
