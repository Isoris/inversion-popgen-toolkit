// =============================================================================
// turn 161 — Band-trace UI layer
// =============================================================================
// SPEC_distant_band_concordance_fish_trajectory.md Slice 4, UI half.
// Compute layer shipped turn 160.
//
// What this turn ships:
//   - state slots: bandTraceFishSet, bandTraceOn, bandTraceCache,
//                  bandTraceCacheKey
//   - setBandTraceFishSet(arr) — persists to localStorage, invalidates cache
//   - setBandTraceOn(bool) — persists to localStorage
//   - _bandTraceFromFocalCandidate(opts) — fills fish-set from the focal
//     candidate's largest band (or specified bandIdx)
//   - _bandTraceCacheKey(chrom, fishSet, K, n_envelopes) — deterministic key
//   - _bandTraceGetOrCompute() — cache-aware accessor
//   - _drawBandTraceStrip(ctx, pad, plotW, plotH, mbMin, mbMax) — Canvas drawer
//   - _bandTraceClearCache() — cache invalidation hook
//   - lines header: checkbox toggle (default OFF) + "🔍 trace" button
//   - drawLinesPanel paint loop: invokes _drawBandTraceStrip when toggle on
//   - applyData chrom-switch: clears cache + fish-set
//
// Tests:
//   1. Source-pattern checks (state slots, exports, button/toggle DOM ids)
//   2. _bandTraceCacheKey is deterministic + sensitive to all inputs
//   3. setBandTraceFishSet round-trips, dedupes, and invalidates cache
//   4. setBandTraceOn round-trips
//   5. _bandTraceFromFocalCandidate picks largest band by default
//   6. _bandTraceFromFocalCandidate respects bandIdx override
//   7. _bandTraceFromFocalCandidate returns null when no candidate
//   8. _drawBandTraceStrip skips silently when toggle off
//   9. _drawBandTraceStrip skips silently when no fish-set
//  10. _drawBandTraceStrip paints expected primitives when active
//  11. Cache key invalidates on chrom switch (different chrom → different key)
//  12. Regression: turn 130 lineage strip + turn 159 propose-all still wired
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
// 1. Source-pattern checks
// ============================================================================
console.log('\n=== 1. Source-pattern checks ===');

// State slots
ok('state.bandTraceFishSet slot in init',
   /bandTraceFishSet:\s*null/.test(html));
ok('state.bandTraceOn slot in init',
   /bandTraceOn:\s*false/.test(html));
ok('state.bandTraceCache slot in init',
   /bandTraceCache:\s*null/.test(html));
ok('state.bandTraceCacheKey slot in init',
   /bandTraceCacheKey:\s*null/.test(html));

// Functions
ok('_drawBandTraceStrip defined',
   /function\s+_drawBandTraceStrip\s*\(/.test(html));
ok('_bandTraceCacheKey defined',
   /function\s+_bandTraceCacheKey\s*\(/.test(html));
ok('_bandTraceGetOrCompute defined',
   /function\s+_bandTraceGetOrCompute\s*\(/.test(html));
ok('_bandTraceFromFocalCandidate defined',
   /function\s+_bandTraceFromFocalCandidate\s*\(/.test(html));
ok('_bandTraceClearCache defined',
   /function\s+_bandTraceClearCache\s*\(/.test(html));
ok('setBandTraceFishSet defined',
   /function\s+setBandTraceFishSet\s*\(/.test(html));
ok('setBandTraceOn defined',
   /function\s+setBandTraceOn\s*\(/.test(html));

// Window exports
ok('window._drawBandTraceStrip exported',
   /window\._drawBandTraceStrip\s*=\s*_drawBandTraceStrip/.test(html));
ok('window.setBandTraceFishSet exported',
   /window\.setBandTraceFishSet\s*=\s*setBandTraceFishSet/.test(html));
ok('window.setBandTraceOn exported',
   /window\.setBandTraceOn\s*=\s*setBandTraceOn/.test(html));
ok('window._bandTraceFromFocalCandidate exported',
   /window\._bandTraceFromFocalCandidate\s*=\s*_bandTraceFromFocalCandidate/.test(html));

// DOM
ok('linesBandTraceToggle checkbox in lines header',
   html.indexOf('id="linesBandTraceToggle"') > -1);
ok('linesBandTraceFromCandBtn button in lines header',
   html.indexOf('id="linesBandTraceFromCandBtn"') > -1);
ok('default-OFF on the checkbox (no `checked` attr)',
   /id="linesBandTraceToggle"[^>]*\/>/.test(html) &&
   !/id="linesBandTraceToggle"[^>]*checked[^>]*\/>/.test(html));

// Wiring
ok('drawLinesPanel calls _drawBandTraceStrip',
   /_drawBandTraceStrip\(ctx,\s*pad,\s*plotW,\s*plotH,\s*mbMin,\s*mbMax\)/.test(html));
ok('applyData clears bandTraceFishSet on chrom switch',
   /state\.bandTraceFishSet\s*=\s*null/.test(html));
ok('applyData calls _bandTraceClearCache',
   /_bandTraceClearCache\(\)/.test(html));

// localStorage keys
ok('_BTRACE_ON_LS_KEY constant declared',
   /_BTRACE_ON_LS_KEY\s*=\s*['"]inversion_atlas\.bandTraceOn['"]/.test(html));
ok('_BTRACE_FISH_SET_LS_KEY constant declared',
   /_BTRACE_FISH_SET_LS_KEY\s*=\s*['"]inversion_atlas\.bandTraceFishSet['"]/.test(html));

// Regime palette
ok('_BTRACE_REGIME_COLOR has co_seg key',
   /co_seg:\s*['"]#7ad394['"]/.test(html));
ok('_BTRACE_REGIME_COLOR has fanned key',
   /fanned:\s*['"]#9aa3ad['"]/.test(html));

// Build sandbox with the helpers + a fake state + a fake localStorage
// + a stub _bandTraceForFishSet so we can isolate UI-layer testing
// from compute-layer (already covered by turn 160).
const sbx = {
  console,
  state: {
    bandTraceFishSet: null,
    bandTraceOn: false,
    bandTraceCache: null,
    bandTraceCacheKey: null,
    k: 3,
    candidate: null,
    data: null,
  },
  window: {},
  _drawCalls: 0,
};
sbx.window.state = sbx.state;
sbx.localStorage = (function() {
  const _store = {};
  return {
    setItem: (k, v) => { _store[k] = String(v); },
    getItem: (k) => (k in _store) ? _store[k] : null,
    removeItem: (k) => { delete _store[k]; },
    _store,
  };
})();
vm.createContext(sbx);

vm.runInContext('var localStorage = this.localStorage;', sbx);
vm.runInContext('function drawLinesPanel(){ _drawCalls++; }', sbx);
// Stub the compute fn — turn 160 already tests the real one.
vm.runInContext(`
  function _bandTraceForFishSet(fishSet, opts) {
    if (!fishSet || !fishSet.length) return null;
    const K = (opts && opts.K) || 3;
    const l2 = (opts && opts.l2_indices) || [];
    return {
      n_fish_selected: fishSet.length,
      n_chains: 1, n_total_L2: l2.length, K,
      per_l2: l2.map((idx, i) => ({
        l2_idx: idx, chain_idx: 0, chain_position: i,
        n_valid: fishSet.length,
        band_counts: new Int32Array(K),
        band_fractions: (function(){
          const f = new Float32Array(K);
          // Fake: i % K gets all the mass
          f[idx % K] = 1;
          return f;
        })(),
        entropy: 0, dominant_band: idx % K, dominant_fraction: 1, regime: 'co_seg',
      })),
    };
  }
`, sbx);

// Inject the real helpers from the atlas
const fnCacheKey = pullFunction(html, '_bandTraceCacheKey');
const fnSetFish  = pullFunction(html, 'setBandTraceFishSet');
const fnSetOn    = pullFunction(html, 'setBandTraceOn');
const fnFromCand = pullFunction(html, '_bandTraceFromFocalCandidate');
const fnGoC      = pullFunction(html, '_bandTraceGetOrCompute');
const fnDraw     = pullFunction(html, '_drawBandTraceStrip');
const fnClear    = pullFunction(html, '_bandTraceClearCache');

ok('pulled _bandTraceCacheKey',          !!fnCacheKey);
ok('pulled setBandTraceFishSet',         !!fnSetFish);
ok('pulled setBandTraceOn',              !!fnSetOn);
ok('pulled _bandTraceFromFocalCandidate', !!fnFromCand);
ok('pulled _bandTraceGetOrCompute',      !!fnGoC);
ok('pulled _drawBandTraceStrip',         !!fnDraw);
ok('pulled _bandTraceClearCache',        !!fnClear);

vm.runInContext('const _BTRACE_ON_LS_KEY = "inversion_atlas.bandTraceOn";', sbx);
vm.runInContext('const _BTRACE_FISH_SET_LS_KEY = "inversion_atlas.bandTraceFishSet";', sbx);
vm.runInContext('const _BTRACE_REGIME_COLOR = {co_seg:"#7ad394",partial:"#f5a524",fanned:"#9aa3ad",sparse:"#5a6068",no_valid:"#3a3e44"};', sbx);
vm.runInContext(fnCacheKey, sbx);
vm.runInContext(fnGoC,      sbx);
vm.runInContext(fnSetFish,  sbx);
vm.runInContext(fnSetOn,    sbx);
vm.runInContext(fnFromCand, sbx);
vm.runInContext(fnClear,    sbx);
vm.runInContext(fnDraw,     sbx);

// ============================================================================
// 2. _bandTraceCacheKey determinism
// ============================================================================
console.log('\n=== 2. Cache key determinism ===');

const k1 = vm.runInContext('_bandTraceCacheKey("LG28", [3, 1, 2], 3, 30)', sbx);
const k2 = vm.runInContext('_bandTraceCacheKey("LG28", [1, 2, 3], 3, 30)', sbx);
ok('order-insensitive in fish-set', k1 === k2);

const k3 = vm.runInContext('_bandTraceCacheKey("LG28", [1, 2, 3], 3, 30)', sbx);
const k4 = vm.runInContext('_bandTraceCacheKey("LG28", [1, 2, 3], 3, 31)', sbx);
ok('sensitive to n_envelopes', k3 !== k4);

const k5 = vm.runInContext('_bandTraceCacheKey("LG14", [1, 2, 3], 3, 30)', sbx);
ok('sensitive to chrom', k1 !== k5);

const k6 = vm.runInContext('_bandTraceCacheKey("LG28", [1, 2, 3], 6, 30)', sbx);
ok('sensitive to K', k1 !== k6);

const k7 = vm.runInContext('_bandTraceCacheKey("LG28", null, 3, 30)', sbx);
ok('null fishSet → null key', k7 === null);

// ============================================================================
// 3. setBandTraceFishSet round-trip + dedupe
// ============================================================================
console.log('\n=== 3. setBandTraceFishSet ===');

vm.runInContext('setBandTraceFishSet([5, 3, 5, 1, 3, 2]);', sbx);
ok('fish-set deduped', sbx.state.bandTraceFishSet.length === 4);
ok('fish-set contents correct',
   JSON.stringify(sbx.state.bandTraceFishSet.slice().sort((a,b)=>a-b)) === JSON.stringify([1,2,3,5]));
ok('fish-set persisted to localStorage',
   sbx.localStorage.getItem('inversion_atlas.bandTraceFishSet') === JSON.stringify([5,3,1,2]));
ok('drawLinesPanel called after set', sbx._drawCalls > 0);

vm.runInContext('setBandTraceFishSet([]);', sbx);
ok('empty array clears fish-set', sbx.state.bandTraceFishSet === null);
ok('empty array clears localStorage',
   sbx.localStorage.getItem('inversion_atlas.bandTraceFishSet') === null);

vm.runInContext('setBandTraceFishSet([10, -1, 20, -5, 30]);', sbx);
ok('negative indices filtered out', sbx.state.bandTraceFishSet.length === 3);

// Cache invalidation
sbx.state.bandTraceCache = { _stub: true };
sbx.state.bandTraceCacheKey = 'old';
vm.runInContext('setBandTraceFishSet([1,2,3]);', sbx);
ok('cache cleared on set', sbx.state.bandTraceCache === null);
ok('cache key cleared on set', sbx.state.bandTraceCacheKey === null);

// ============================================================================
// 4. setBandTraceOn round-trip
// ============================================================================
console.log('\n=== 4. setBandTraceOn ===');

vm.runInContext('setBandTraceOn(true);', sbx);
ok('state.bandTraceOn = true', sbx.state.bandTraceOn === true);
ok('persisted to LS as 1', sbx.localStorage.getItem('inversion_atlas.bandTraceOn') === '1');

vm.runInContext('setBandTraceOn(false);', sbx);
ok('state.bandTraceOn = false', sbx.state.bandTraceOn === false);
ok('persisted to LS as 0', sbx.localStorage.getItem('inversion_atlas.bandTraceOn') === '0');

// Boolean coercion
vm.runInContext('setBandTraceOn("yes");', sbx);
ok('"yes" coerced to true', sbx.state.bandTraceOn === true);
vm.runInContext('setBandTraceOn(0);', sbx);
ok('0 coerced to false', sbx.state.bandTraceOn === false);

// ============================================================================
// 5. _bandTraceFromFocalCandidate — largest-band default
// ============================================================================
console.log('\n=== 5. _bandTraceFromFocalCandidate default (largest band) ===');

// 6 samples, K=3, labels = [0, 0, 1, 1, 1, 2] → band 1 is largest (3 members)
sbx.state.candidate = {
  id: 'I1', K: 3, locked_labels: [0, 0, 1, 1, 1, 2],
};
const r5 = vm.runInContext('_bandTraceFromFocalCandidate()', sbx);
ok('returns array of size 3', Array.isArray(r5) && r5.length === 3);
ok('picked band 1 members [2,3,4]',
   JSON.stringify(r5.slice().sort((a,b)=>a-b)) === JSON.stringify([2,3,4]));
ok('side effect: state.bandTraceFishSet updated',
   sbx.state.bandTraceFishSet && sbx.state.bandTraceFishSet.length === 3);

// ============================================================================
// 6. _bandTraceFromFocalCandidate — bandIdx override
// ============================================================================
console.log('\n=== 6. bandIdx override ===');

const r6 = vm.runInContext('_bandTraceFromFocalCandidate({bandIdx: 0})', sbx);
ok('bandIdx=0 picks band 0 members [0,1]',
   JSON.stringify(r6.slice().sort((a,b)=>a-b)) === JSON.stringify([0,1]));

const r6b = vm.runInContext('_bandTraceFromFocalCandidate({bandIdx: 2})', sbx);
ok('bandIdx=2 picks band 2 members [5]',
   JSON.stringify(r6b) === JSON.stringify([5]));

// Out-of-range bandIdx falls back to largest
const r6c = vm.runInContext('_bandTraceFromFocalCandidate({bandIdx: 99})', sbx);
ok('out-of-range bandIdx falls back to largest band',
   r6c.length === 3);

// ============================================================================
// 7. _bandTraceFromFocalCandidate — null when no candidate
// ============================================================================
console.log('\n=== 7. No-candidate / no-labels handling ===');

sbx.state.candidate = null;
const r7a = vm.runInContext('_bandTraceFromFocalCandidate()', sbx);
ok('no candidate → null', r7a === null);

sbx.state.candidate = { id: 'X', K: 3, locked_labels: null };
const r7b = vm.runInContext('_bandTraceFromFocalCandidate()', sbx);
ok('candidate without locked_labels → null', r7b === null);

sbx.state.candidate = { id: 'X', K: 3, locked_labels: [-1, -1, -1] };
const r7c = vm.runInContext('_bandTraceFromFocalCandidate()', sbx);
ok('all-invalid labels → null', r7c === null);

// ============================================================================
// 8 & 9. _drawBandTraceStrip skip conditions
// ============================================================================
console.log('\n=== 8 & 9. Drawer skip conditions ===');

// Build a fake ctx that records calls
function makeFakeCtx() {
  const calls = [];
  return {
    save: () => calls.push(['save']),
    restore: () => calls.push(['restore']),
    fillRect: (...a) => calls.push(['fillRect', ...a]),
    strokeRect: (...a) => calls.push(['strokeRect', ...a]),
    set fillStyle(v) { this._fill = v; calls.push(['fillStyle', v]); },
    get fillStyle() { return this._fill; },
    set strokeStyle(v) { this._stroke = v; calls.push(['strokeStyle', v]); },
    get strokeStyle() { return this._stroke; },
    set lineWidth(v) { calls.push(['lineWidth', v]); },
    set globalAlpha(v) { calls.push(['globalAlpha', v]); },
    _calls: calls,
  };
}

// Setup data so the drawer would have something to paint
sbx.state.data = {
  chrom: 'LG28',
  l2_envelopes: [
    { start_bp: 1_000_000, end_bp: 2_000_000 },
    { start_bp: 2_500_000, end_bp: 3_000_000 },
    { start_bp: 3_500_000, end_bp: 4_000_000 },
  ],
};

// Toggle off: drawer should produce no calls beyond an early return
sbx.state.bandTraceOn = false;
sbx.state.bandTraceFishSet = [1, 2, 3];
sbx._fakeCtx = makeFakeCtx();
vm.runInContext('_drawBandTraceStrip(_fakeCtx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx);
ok('toggle off → drawer is no-op', sbx._fakeCtx._calls.length === 0);

// Toggle on but no fish-set: same
sbx.state.bandTraceOn = true;
sbx.state.bandTraceFishSet = null;
sbx._fakeCtx = makeFakeCtx();
vm.runInContext('_drawBandTraceStrip(_fakeCtx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx);
ok('no fish-set → drawer is no-op', sbx._fakeCtx._calls.length === 0);

// Toggle on + fish-set + data: should paint
sbx.state.bandTraceOn = true;
sbx.state.bandTraceFishSet = [1, 2, 3];
sbx.state.bandTraceCache = null;
sbx.state.bandTraceCacheKey = null;
sbx._fakeCtx = makeFakeCtx();
vm.runInContext('_drawBandTraceStrip(_fakeCtx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx);
ok('drawer painted backdrop fillRect',
   sbx._fakeCtx._calls.some(c => c[0] === 'fillRect'));
ok('drawer painted at least one regime stripe (3 L2s × top stripe + stacked = ≥4 fillRects)',
   sbx._fakeCtx._calls.filter(c => c[0] === 'fillRect').length >= 4);
ok('drawer drew frame strokeRect',
   sbx._fakeCtx._calls.some(c => c[0] === 'strokeRect'));
ok('drawer cache populated after first paint',
   sbx.state.bandTraceCache !== null && sbx.state.bandTraceCacheKey !== null);

// Second paint should hit cache (no recompute)
const cacheBefore = sbx.state.bandTraceCache;
vm.runInContext('_drawBandTraceStrip(_fakeCtx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx);
ok('second paint reuses cached trace (same object identity)',
   sbx.state.bandTraceCache === cacheBefore);

// ============================================================================
// 10. Cache key invalidates on chrom switch
// ============================================================================
console.log('\n=== 10. Cache key invalidates on chrom switch ===');

const oldKey = sbx.state.bandTraceCacheKey;
sbx.state.data.chrom = 'LG14';   // simulate chrom switch
sbx._fakeCtx = makeFakeCtx();
vm.runInContext('_drawBandTraceStrip(_fakeCtx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx);
ok('cache key changed on chrom switch', sbx.state.bandTraceCacheKey !== oldKey);

// ============================================================================
// 11. _bandTraceClearCache
// ============================================================================
console.log('\n=== 11. _bandTraceClearCache ===');

sbx.state.bandTraceCache = { _stub: true };
sbx.state.bandTraceCacheKey = 'something';
vm.runInContext('_bandTraceClearCache();', sbx);
ok('cache cleared', sbx.state.bandTraceCache === null);
ok('key cleared', sbx.state.bandTraceCacheKey === null);

// ============================================================================
// 12. Regression
// ============================================================================
console.log('\n=== 12. Regression checks ===');

ok('turn 130 _drawLineageStrip still defined',
   /function\s+_drawLineageStrip\s*\(/.test(html));
ok('turn 130 lineage toggle still in place',
   html.indexOf('id="linesLineageStripToggle"') > -1);
ok('turn 159 _gpKaryoProposeAll still defined',
   /function\s+_gpKaryoProposeAll\s*\(/.test(html));
ok('turn 160 _bandTraceForFishSet still defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('turn 160 _bandTraceRegimeRuns still defined',
   /function\s+_bandTraceRegimeRuns\s*\(/.test(html));
ok('turn 156 V-shape diagnostic still present',
   /function\s+openVShapePlot\s*\(/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
