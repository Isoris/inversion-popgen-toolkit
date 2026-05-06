// =============================================================================
// turn 155 — Threshold folded into _inheritanceCacheKey
// =============================================================================
// Pre-turn-155 the cache key only captured (mode, items). Two computes at
// different thresholds would BOTH look like the same cached entry, causing
// runInheritanceCompute's cache-hit branch to return stale clusters when the
// threshold changed. The workaround in _gpInhSetThreshold (calling
// invalidateInheritanceCache from the slider) was paper-thin: any consumer
// that built an expected key without first reading the threshold could
// mis-detect staleness, and any caller that passed { threshold: x } without
// updating the slider state would write a wrongly-keyed cache entry.
//
// What this turn changes:
//   1. _inheritanceCacheKey signature: (items, mode) → (items, mode, threshold).
//      threshold is OPTIONAL — defaults to state.gPanelInheritanceThreshold,
//      then to _IGC_DEFAULT_COSINE_DIST_THRESHOLD. Backwards compatible: old
//      (items, mode) callers all keep working with the new default.
//   2. runInheritanceCompute resolves threshold BEFORE building the cache key
//      so the key folds in the correct value. Default threshold also now
//      reads state.gPanelInheritanceThreshold so all consumers agree.
//   3. inheritanceGroupClustering call is given the resolved threshold
//      explicitly (was: forwarded original opts; could disagree with the
//      cache key's threshold).
//   4. _gpInhSetThreshold's invalidateInheritanceCache call kept as
//      defense-in-depth (comment updated).
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

// ============================================================================
// 1. Source-pattern: signature extended
// ============================================================================
console.log('\n=== 1. Signature extended with threshold ===');

ok('_inheritanceCacheKey takes 3 params (items, mode, threshold)',
   /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(html));

const ckRe = /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)\s*\{([\s\S]*?)\nfunction\s/;
const ckMatch = html.match(ckRe);
ok('_inheritanceCacheKey body extracted', !!ckMatch);
const ckBody = ckMatch ? ckMatch[1] : '';

ok('threshold defaults to state.gPanelInheritanceThreshold when null',
   /_state\.gPanelInheritanceThreshold/.test(ckBody) &&
   /Number\.isFinite\s*\(\s*_state\.gPanelInheritanceThreshold\s*\)/.test(ckBody));

ok('falls back to _IGC_DEFAULT_COSINE_DIST_THRESHOLD when state has no threshold',
   /_IGC_DEFAULT_COSINE_DIST_THRESHOLD/.test(ckBody));

ok('rounds to 4 decimals to avoid floating-point key fragmentation',
   /Math\.round\s*\(\s*\(?\+?t/.test(ckBody) &&
   /1e4/.test(ckBody) &&
   /toFixed\s*\(\s*4\s*\)/.test(ckBody));

ok('key string includes the threshold ("@t" prefix)',
   /@t\$\{tStr\}/.test(ckBody) || /\$\{mode\}@t/.test(ckBody));

// ============================================================================
// 2. Source-pattern: runInheritanceCompute uses resolved threshold
// ============================================================================
console.log('\n=== 2. runInheritanceCompute uses resolved threshold ===');

const ricRe = /function\s+runInheritanceCompute\s*\(\s*opts\s*\)\s*\{([\s\S]*?)\nfunction\s+invalidateInheritanceCache/;
const ricMatch = html.match(ricRe);
ok('runInheritanceCompute body extracted', !!ricMatch);
const ricBody = ricMatch ? ricMatch[1] : '';

ok('threshold resolved BEFORE cache-key build',
   ricBody.indexOf('let threshold') < ricBody.indexOf('_inheritanceCacheKey(items, mode, threshold)'));

ok('cache key is built with the resolved threshold',
   /_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(ricBody));

ok('threshold falls back through state.gPanelInheritanceThreshold',
   /_state\.gPanelInheritanceThreshold/.test(ricBody));

ok('inheritanceGroupClustering receives the resolved threshold explicitly',
   /inheritanceGroupClustering\s*\(\s*items\s*,\s*Object\.assign\s*\([^)]*threshold:\s*threshold/.test(ricBody));

// ============================================================================
// 3. Source-pattern: _gpInhSetThreshold workaround retained as defense-in-depth
// ============================================================================
console.log('\n=== 3. _gpInhSetThreshold defense-in-depth ===');

const setRe = /function\s+_gpInhSetThreshold\s*\(\s*t\s*\)\s*\{([\s\S]*?)\n(?:\}\s*\n\s*\/\/|\}\s*\nfunction)/;
const setMatch = html.match(setRe);
ok('_gpInhSetThreshold body extracted', !!setMatch);
const setBody = setMatch ? setMatch[1] : '';

ok('still calls invalidateInheritanceCache (defense-in-depth, retained)',
   /invalidateInheritanceCache\s*\(\s*\)/.test(setBody));

ok('comment explains the cache-key fix makes invalidate redundant for runInheritanceCompute',
   /defense-in-depth/.test(setBody) || /turn 155/.test(setBody));

// ============================================================================
// 4. Sandboxed: keys differ when threshold differs
// ============================================================================
console.log('\n=== 4. Keys differ when threshold differs (sandboxed) ===');

function makeSandbox(opts) {
  opts = opts || {};
  const ctx = {
    state: opts.state || {},
    Float32Array, Int32Array, Uint8Array, Map, Set, Number, console,
    Object, Array, JSON, String, Math,
  };
  ctx.window = opts.window || {};
  ctx.window.state = ctx.state;
  vm.createContext(ctx);
  // Inject the constants and helpers _inheritanceCacheKey references.
  vm.runInContext(
    'const _IGC_DEFAULT_COSINE_DIST_THRESHOLD = 0.15;\n' +
    'function _hashLockedLabels(labels) {\n' +
    '  if (!labels || !labels.length) return 0;\n' +
    '  let h = 0x811c9dc5 | 0;\n' +
    '  for (let i = 0; i < labels.length; i++) h = (h * 31 + (labels[i] + 2)) | 0;\n' +
    '  return h >>> 0;\n' +
    '}\n',
    ctx
  );
  const fnCk = pullFunction(html, '_inheritanceCacheKey');
  if (fnCk) vm.runInContext(fnCk, ctx);
  return ctx;
}

// 4a. Same items + same mode + DIFFERENT thresholds → different keys
{
  const sb = makeSandbox();
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  const k015 = sb._inheritanceCacheKey(items, 'default', 0.15);
  const k020 = sb._inheritanceCacheKey(items, 'default', 0.20);
  const k025 = sb._inheritanceCacheKey(items, 'default', 0.25);
  ok('threshold 0.15 vs 0.20 → different keys', k015 !== k020);
  ok('threshold 0.20 vs 0.25 → different keys', k020 !== k025);
  ok('threshold 0.15 vs 0.25 → different keys', k015 !== k025);
}

// 4b. Same items + same mode + same threshold → same key (idempotent)
{
  const sb = makeSandbox();
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  ok('same inputs → identical keys',
     sb._inheritanceCacheKey(items, 'default', 0.15) ===
     sb._inheritanceCacheKey(items, 'default', 0.15));
}

// 4c. Tiny floating-point jitter coalesces to same key
{
  const sb = makeSandbox();
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  // 0.15 vs 0.150000000001 — slider math jitter
  const a = sb._inheritanceCacheKey(items, 'default', 0.15);
  const b = sb._inheritanceCacheKey(items, 'default', 0.150000000001);
  ok('floating-point jitter (0.15 vs 0.150000000001) coalesces to same key',
     a === b);
}

// 4d. Different items at same threshold → different keys
{
  const sb = makeSandbox();
  const itemsA = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  const itemsB = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'C', K: 3, start_bp: 2000, end_bp: 3000, labels: new Int8Array([0,1,2]) },   // different item
  ];
  ok('changing items → different keys (existing turn 129 contract preserved)',
     sb._inheritanceCacheKey(itemsA, 'default', 0.15) !==
     sb._inheritanceCacheKey(itemsB, 'default', 0.15));
}

// 4e. Different mode at same threshold + items → different keys
{
  const sb = makeSandbox();
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  ok('changing mode → different keys (existing turn 129 contract preserved)',
     sb._inheritanceCacheKey(items, 'default', 0.15) !==
     sb._inheritanceCacheKey(items, 'detailed', 0.15));
}

// ============================================================================
// 5. Sandboxed: optional threshold defaults from state
// ============================================================================
console.log('\n=== 5. Optional threshold defaults to state ===');

// 5a. Omitting threshold pulls from state.gPanelInheritanceThreshold
{
  const sb = makeSandbox({
    state: { gPanelInheritanceThreshold: 0.20 },
  });
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  const explicitKey = sb._inheritanceCacheKey(items, 'default', 0.20);
  const defaultKey = sb._inheritanceCacheKey(items, 'default');   // omit threshold
  ok('omitted threshold uses state.gPanelInheritanceThreshold = 0.20',
     explicitKey === defaultKey);
}

// 5b. State threshold absent → falls back to _IGC_DEFAULT (0.15)
{
  const sb = makeSandbox({ state: {} });
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  const explicitKey = sb._inheritanceCacheKey(items, 'default', 0.15);
  const defaultKey = sb._inheritanceCacheKey(items, 'default');   // omit threshold
  ok('omitted threshold + no state → falls back to _IGC_DEFAULT (0.15)',
     explicitKey === defaultKey);
}

// 5c. Backwards compat — old (items, mode) callers still work
{
  const sb = makeSandbox({ state: { gPanelInheritanceThreshold: 0.15 } });
  const items = [
    { id: 'A', K: 3, start_bp: 0, end_bp: 1000, labels: new Int8Array([0,1,2]) },
    { id: 'B', K: 3, start_bp: 1000, end_bp: 2000, labels: new Int8Array([0,1,2]) },
  ];
  let key = null;
  let threw = null;
  try {
    key = sb._inheritanceCacheKey(items, 'default');   // 2-arg call
  } catch (e) { threw = e; }
  ok('2-arg call (items, mode) still works (backwards compat)',
     threw === null && key !== null && typeof key === 'string');
}

// ============================================================================
// 6. Existing test contracts preserved
// ============================================================================
console.log('\n=== 6. Existing contracts preserved ===');

ok('runInheritanceCompute still exists',
   /function\s+runInheritanceCompute\s*\(\s*opts\s*\)/.test(html));

ok('inheritanceGroupClustering still exists',
   /function\s+inheritanceGroupClustering\s*\(\s*items\s*,\s*opts\s*\)/.test(html));

ok('invalidateInheritanceCache still exists',
   /function\s+invalidateInheritanceCache\s*\(\s*\)/.test(html));

ok('turn 153 _autoRegisterInheritanceOnCandidateChange still exists',
   /function\s+_autoRegisterInheritanceOnCandidateChange\s*\(\s*\)/.test(html));

ok('turn 154 inheritanceLastStatus still tracked',
   /inheritanceLastStatus/.test(html));

ok('cache-hit fast path still exists in runInheritanceCompute',
   /_state\.inheritanceCacheKey\s*===\s*cacheKey/.test(ricBody));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
