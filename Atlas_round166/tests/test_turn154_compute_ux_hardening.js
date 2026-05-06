// =============================================================================
// turn 154 — G-panel inheritance compute UX hardening
// =============================================================================
// Pre-turn-154 the user clicked ↻ recompute and saw nothing change in the UI
// when compute returned null (insufficient items, all-K=0, exception). The
// "No compute result yet" placeholder rendered identically to "you haven't
// clicked yet". And the tab did NOT auto-compute on first open — every fresh
// open required a manual click even when items were stable.
//
// What this turn adds:
//   1. state.inheritanceLastStatus written by every runInheritanceCompute
//      call (success / cached / insufficient_items / null_result / exception)
//   2. _gPanelRenderTabInheritance auto-computes on mount when stale
//   3. Status pill above the threshold row surfaces last-compute outcome
//   4. Distinct empty-result message ("succeeded but produced 0 groups")
//      vs failed-compute message ("see status above")
//   5. ↻ recompute button flashes "computing…" and disables during click
//   6. invalidateInheritanceCache also clears inheritanceLastStatus
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
// 1. Source-pattern: runInheritanceCompute writes status on every path
// ============================================================================
console.log('\n=== 1. runInheritanceCompute status tracking ===');

const ricRe = /function\s+runInheritanceCompute\s*\(\s*opts\s*\)\s*\{([\s\S]*?)\nfunction\s+invalidateInheritanceCache/;
const ricMatch = html.match(ricRe);
ok('runInheritanceCompute body extracted', !!ricMatch);
const ricBody = ricMatch ? ricMatch[1] : '';

ok('declares _setStatus helper', /_setStatus\s*=\s*\(s\)\s*=>/.test(ricBody));

ok('writes status on insufficient_items path',
   /reason:\s*'insufficient_items'/.test(ricBody));

ok('writes status on cached-hit path',
   /reason:\s*'cached'/.test(ricBody));

ok('writes status on exception path',
   /reason:\s*'exception'/.test(ricBody));

ok('writes status on null_result path',
   /reason:\s*'null_result'/.test(ricBody));

ok('writes status on successful compute path',
   /reason:\s*'computed'/.test(ricBody));

ok('wraps inheritanceGroupClustering in try/catch (turn 154)',
   /try\s*\{[\s\S]*?result\s*=\s*inheritanceGroupClustering/.test(ricBody));

ok('status includes ok / message / n_items fields',
   /ok:\s*true/.test(ricBody) && /message:/.test(ricBody) && /n_items:/.test(ricBody));

ok('status timestamps via Date.now()',
   /Date\.now\s*\(\s*\)/.test(ricBody));

// ============================================================================
// 2. Source-pattern: invalidateInheritanceCache clears status
// ============================================================================
console.log('\n=== 2. invalidateInheritanceCache clears status ===');

const invRe = /function\s+invalidateInheritanceCache\s*\(\s*\)\s*\{([\s\S]*?)\n\}/;
const invMatch = html.match(invRe);
const invBody = invMatch ? invMatch[1] : '';

ok('invalidateInheritanceCache also nulls inheritanceLastStatus',
   /_state\.inheritanceLastStatus\s*=\s*null/.test(invBody));

// ============================================================================
// 3. Source-pattern: tab body auto-computes on mount + renders status pill
// ============================================================================
console.log('\n=== 3. Tab body auto-compute + status pill ===');

const tabRe = /function\s+_gPanelRenderTabInheritance\s*\(\s*\)\s*\{([\s\S]*?)\n\/\/\s*=+/;
const tabMatch = html.match(tabRe);
ok('tab body extracted', !!tabMatch);
const tabBody = tabMatch ? tabMatch[1] : '';

ok('tab body auto-computes when items.length >= 2 AND cache stale',
   /items\.length\s*>=\s*2/.test(tabBody) &&
   /cacheStale/.test(tabBody) &&
   /runInheritanceCompute\s*\(\s*\{\s*threshold:\s*threshold\s*\}\s*\)/.test(tabBody));

ok('auto-compute uses _inheritanceCacheKey to detect staleness',
   /_inheritanceCacheKey\s*\(\s*items\s*,/.test(tabBody));

ok('auto-compute does NOT pass force: true (lets cache hit)',
   !/force:\s*true[\s\S]*items\.length\s*>=\s*2/.test(tabBody.split('Threshold slider')[0] || tabBody));

ok('renders status pill (id=gpInhStatusPill)',
   /id="gpInhStatusPill"/.test(tabBody));

ok('status pill colors green when ok=true',
   /var\(--good\)/.test(tabBody) && /status\.ok/.test(tabBody));

ok('status pill colors amber on exception/null_result',
   /var\(--accent[^)]*\)/.test(tabBody));

ok('status pill renders last-compute message',
   /status\.message/.test(tabBody));

ok('status pill renders different icon per outcome',
   /●/.test(tabBody) && /⚠/.test(tabBody));

// ============================================================================
// 4. Source-pattern: tab body distinguishes 0-groups from null-result
// ============================================================================
console.log('\n=== 4. Tab body — distinct messages for 0-groups vs null-result ===');

ok('tab body has dedicated branch for n_groups === 0',
   /n_groups\s*===\s*0/.test(tabBody));

ok('0-groups branch suggests lowering threshold',
   /Try lowering the threshold/.test(tabBody) ||
   /lower.*threshold/i.test(tabBody));

ok('null-result branch points user to status pill',
   /No groups to display.*see status above/.test(tabBody) ||
   /see status above/.test(tabBody));

// ============================================================================
// 5. Source-pattern: recompute button flash
// ============================================================================
console.log('\n=== 5. ↻ recompute click feedback ===');

// Extract the recompute click handler — it lives inside _renderGPanelModal
const modalRe = /function\s+_renderGPanelModal\s*\(\s*\)\s*\{([\s\S]*?)\nfunction\s+_gPanelClose/;
const modalBody = (html.match(modalRe) || ['', ''])[1];

ok('recompute click sets button text to "computing…"',
   /textContent\s*=\s*'computing…'/.test(modalBody));

ok('recompute click disables button during compute',
   /recomputeBtn\.disabled\s*=\s*true/.test(modalBody));

ok('on re-render fail, recompute button restored',
   /textContent\s*=\s*orig/.test(modalBody) &&
   /disabled\s*=\s*false/.test(modalBody));

ok('on exception, recompute button still restored',
   (modalBody.match(/recomputeBtn\.textContent\s*=\s*orig/g) || []).length >= 2);

// ============================================================================
// 6. Sandboxed: runInheritanceCompute status writes
// ============================================================================
console.log('\n=== 6. Sandboxed status-write behaviour ===');

function pullFunction(src, name, terminator) {
  const decl = `function ${name}`;
  const i = src.indexOf(decl);
  if (i < 0) return null;
  const p = src.indexOf('(', i);
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
      inheritanceResult: null,
      inheritanceCacheKey: null,
      inheritanceLastStatus: null,
    }, opts.state || {}),
    Date, Object, Array, JSON, Math, console, String, Number,
    Float32Array, Int32Array, Uint8Array,
  };
  ctx.window = { state: ctx.state };
  vm.createContext(ctx);
  // Constants the helper reads
  vm.runInContext('const _IGC_DEFAULT_COSINE_DIST_THRESHOLD = 0.15;', ctx);
  // Mockable upstream helpers
  ctx._mockGatherItems = opts.gatherItems != null ? opts.gatherItems : [];
  ctx._mockCacheKey = opts.cacheKey || 'k1';
  ctx._mockClusterReturn = opts.clusterReturn !== undefined ? opts.clusterReturn : null;
  ctx._mockClusterThrows = opts.clusterThrows || null;
  vm.runInContext(
    'var _gatherActiveCandidatesForInheritance = function() { ' +
    '  return _mockGatherItems; };\n' +
    'var _inheritanceCacheKey = function(items, mode) { ' +
    '  return _mockCacheKey; };\n' +
    'var inheritanceGroupClustering = function(items, opts) { ' +
    '  if (_mockClusterThrows) throw new Error(_mockClusterThrows); ' +
    '  return _mockClusterReturn; };\n',
    ctx
  );
  // Inject the real runInheritanceCompute
  const fn = pullFunction(html, 'runInheritanceCompute');
  vm.runInContext(fn, ctx);
  return ctx;
}

// 6a. insufficient_items status
{
  const sb = makeSandbox({ gatherItems: [] });
  const out = sb.runInheritanceCompute({ force: true });
  ok('insufficient_items: returns null', out === null);
  const s = sb.state.inheritanceLastStatus;
  ok('insufficient_items: status.ok=false', s && s.ok === false);
  ok('insufficient_items: reason set', s && s.reason === 'insufficient_items');
  ok('insufficient_items: n_items=0', s && s.n_items === 0);
  ok('insufficient_items: message includes "Need ≥2"', s && /Need\s+≥2/.test(s.message));
}

// 6b. cached-hit status
{
  const cachedResult = { rtab: { group_ids: [0, 1, 2] } };
  const sb = makeSandbox({
    gatherItems: [{id:'a'},{id:'b'}],
    state: {
      activeMode: 'default',
      inheritanceCacheKey: 'k1',
      inheritanceResult: cachedResult,
    },
  });
  const out = sb.runInheritanceCompute();   // no force
  ok('cached-hit: returns the cached result', out === cachedResult);
  const s = sb.state.inheritanceLastStatus;
  ok('cached-hit: status.ok=true', s && s.ok === true);
  ok('cached-hit: reason="cached"', s && s.reason === 'cached');
  ok('cached-hit: n_groups read from cached result', s && s.n_groups === 3);
}

// 6c. exception status
{
  const sb = makeSandbox({
    gatherItems: [{id:'a',K:3,labels:[0,1,2]},{id:'b',K:3,labels:[0,1,2]}],
    clusterThrows: 'mock cluster failure',
  });
  const out = sb.runInheritanceCompute({ force: true });
  ok('exception: returns null', out === null);
  const s = sb.state.inheritanceLastStatus;
  ok('exception: status.ok=false', s && s.ok === false);
  ok('exception: reason="exception"', s && s.reason === 'exception');
  ok('exception: error string captured',
     s && /mock cluster failure/.test(s.error || ''));
  ok('exception: state.inheritanceResult cleared',
     sb.state.inheritanceResult === null);
}

// 6d. null_result status
{
  const sb = makeSandbox({
    gatherItems: [{id:'a',K:3,labels:[0,1,2]},{id:'b',K:3,labels:[0,1,2]}],
    clusterReturn: null,
  });
  const out = sb.runInheritanceCompute({ force: true });
  ok('null_result: returns null', out === null);
  const s = sb.state.inheritanceLastStatus;
  ok('null_result: status.ok=false', s && s.ok === false);
  ok('null_result: reason="null_result"', s && s.reason === 'null_result');
  ok('null_result: message hints at K=0 issue',
     s && /K\s*values|usable bands/.test(s.message));
}

// 6e. computed status
{
  const sb = makeSandbox({
    gatherItems: [
      {id:'a',K:3,labels:[0,1,2],start_bp:0,end_bp:1e6,seq_num:1},
      {id:'b',K:3,labels:[0,1,2],start_bp:1e6,end_bp:2e6,seq_num:2},
    ],
    clusterReturn: { rtab: { group_ids: [10, 20] } },
    cacheKey: 'fresh_key',
  });
  const out = sb.runInheritanceCompute({ force: true, threshold: 0.20 });
  ok('computed: returns the new result', out !== null);
  ok('computed: result has items_meta attached',
     out.items_meta && out.items_meta.length === 2);
  const s = sb.state.inheritanceLastStatus;
  ok('computed: status.ok=true', s && s.ok === true);
  ok('computed: reason="computed"', s && s.reason === 'computed');
  ok('computed: n_groups=2', s && s.n_groups === 2);
  ok('computed: threshold=0.20', s && Math.abs(s.threshold - 0.20) < 1e-9);
  ok('computed: state.inheritanceCacheKey updated',
     sb.state.inheritanceCacheKey === 'fresh_key');
}

// 6f. status timestamp
{
  const sb = makeSandbox({ gatherItems: [] });
  const beforeT = Date.now();
  sb.runInheritanceCompute({ force: true });
  const afterT = Date.now();
  const s = sb.state.inheritanceLastStatus;
  ok('status.at is a recent timestamp',
     s && Number.isFinite(s.at) && s.at >= beforeT && s.at <= afterT);
}

// ============================================================================
// 7. Sandboxed: invalidateInheritanceCache clears status
// ============================================================================
console.log('\n=== 7. invalidateInheritanceCache clears status ===');

{
  const ctx = {
    state: {
      inheritanceResult: { foo: 1 },
      inheritanceCacheKey: 'k1',
      inheritanceLastStatus: { ok: true, message: 'test' },
    },
    console,
  };
  ctx.window = { state: ctx.state };
  vm.createContext(ctx);
  const fn = pullFunction(html, 'invalidateInheritanceCache');
  vm.runInContext(fn, ctx);
  ctx.invalidateInheritanceCache();
  ok('invalidate clears inheritanceResult', ctx.state.inheritanceResult === null);
  ok('invalidate clears inheritanceCacheKey', ctx.state.inheritanceCacheKey === null);
  ok('invalidate clears inheritanceLastStatus', ctx.state.inheritanceLastStatus === null);
}

// ============================================================================
// 8. Existing flow preserved
// ============================================================================
console.log('\n=== 8. Existing flow preserved ===');

ok('_gPanelRenderTabInheritance still exists',
   /function\s+_gPanelRenderTabInheritance\s*\(\s*\)/.test(html));

ok('runInheritanceCompute still augments result with items_meta',
   /result\.items_meta\s*=/.test(html));

ok('_gpInhRenderResult still exists',
   /function\s+_gpInhRenderResult\s*\(\s*result\s*,\s*items\s*,\s*threshold\s*\)/.test(html));

ok('cache invalidation on threshold change still works (turn 152 contract)',
   /invalidateInheritanceCache\s*\(\s*\)/.test(html));

ok('turn 153 auto-register hook still in persistCandidateList',
   /_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
