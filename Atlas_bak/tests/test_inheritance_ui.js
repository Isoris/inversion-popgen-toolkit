// Tests for turn 2c — inheritance UI integration.
// Tests the JavaScript logic (gather, cache, format, compute trigger).
// Canvas drawing isn't directly verifiable in JSDOM but we exercise the
// drawing code path to catch crashes.

const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
});
const w = dom.window;

function run() {
  const failures = [];
  let testNum = 0;
  function t(name, fn) {
    testNum++;
    try { fn(); console.log(`  PASS [${testNum}] ${name}`); }
    catch (e) { failures.push({ name, err: e.message }); console.log(`  FAIL [${testNum}] ${name}: ${e.message}`); }
  }
  function eq(a, b, m) { if (a !== b) throw new Error(`${m||''} expected ${JSON.stringify(b)}, got ${JSON.stringify(a)}`); }

  if (typeof w.runInheritanceCompute !== 'function'
      || typeof w.invalidateInheritanceCache !== 'function'
      || typeof w._formatInheritanceLabel !== 'function'
      || typeof w._gatherActiveCandidatesForInheritance !== 'function'
      || typeof w._drawInheritanceLabelsStrip !== 'function') {
    console.log('  FAIL: turn 2c functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2c functions exposed on window');

  // Helper: synthesize a state with N candidates carrying labels arrays
  function setupState(candConfigs, mode) {
    // candConfigs: [{id, K, start_bp, end_bp, labels}]
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = mode || 'default';
    _state.k = (mode === 'detailed') ? 6 : 3;
    const target = (_state.activeMode === 'detailed') ? 'candidates_detailed' : 'candidates';
    _state[target] = {};
    for (const c of candConfigs) {
      _state[target][c.id] = {
        id: c.id,
        K: c.K,
        K_used: c.K,
        locked_labels: c.labels,
        start_bp: c.start_bp,
        end_bp: c.end_bp,
        source: 'L2',
      };
    }
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
    _state._inheritanceComputeScheduled = false;
    _state.linesInheritanceLabelsOn = true;
  }

  // ============================================================
  // _formatInheritanceLabel
  // ============================================================
  t('format label: single candidate, 3 groups -> "I1·3g"', () => {
    eq(w._formatInheritanceLabel(1, 1, 3), 'I1·3g');
  });
  t('format label: range 3-5, 6 groups -> "I3-5·6g"', () => {
    eq(w._formatInheritanceLabel(3, 5, 6), 'I3-5·6g');
  });
  t('format label: single, 1 group -> "I1·1g"', () => {
    eq(w._formatInheritanceLabel(1, 1, 1), 'I1·1g');
  });

  // ============================================================
  // _gatherActiveCandidatesForInheritance
  // ============================================================
  t('gather: returns sorted candidates with seq_num', () => {
    // Note: out of order start_bp to test sorting
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels },
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    eq(items.length, 3);
    eq(items[0].id, 'A');
    eq(items[1].id, 'B');
    eq(items[2].id, 'C');
    eq(items[0].seq_num, 1);
    eq(items[1].seq_num, 2);
    eq(items[2].seq_num, 3);
  });

  t('gather: skips candidates with no locked_labels', () => {
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: [] }]);
    const items = w._gatherActiveCandidatesForInheritance();
    eq(items.length, 0);
  });

  t('gather: respects activeMode', () => {
    const labels = [];
    for (let k = 0; k < 6; k++) for (let i = 0; i < 20; i++) labels.push(k);
    setupState([
      { id: 'X', K: 6, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], 'detailed');
    const items = w._gatherActiveCandidatesForInheritance();
    eq(items.length, 1);
    eq(items[0].K, 6);
  });

  t('gather: empty state returns empty array', () => {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.candidates = {};
    _state.candidates_detailed = {};
    eq(w._gatherActiveCandidatesForInheritance().length, 0);
  });

  // ============================================================
  // runInheritanceCompute
  // ============================================================
  t('runInheritanceCompute: returns null with <2 candidates', () => {
    const labels = [0,1,2,0,1,2];
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    eq(w.runInheritanceCompute(), null);
  });

  t('runInheritanceCompute: produces result with 2 perfectly linked candidates', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    const r = w.runInheritanceCompute();
    if (!r) throw new Error('expected result, got null');
    eq(r.cut.n_groups, 3);
    eq(r.items_meta.length, 2);
    eq(r.items_meta[0].seq_num, 1);
    eq(r.items_meta[1].seq_num, 2);
  });

  t('runInheritanceCompute: caches result, returns same on repeat call', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    const r1 = w.runInheritanceCompute();
    const r2 = w.runInheritanceCompute();
    eq(r1, r2);   // same object identity = cache hit
  });

  t('runInheritanceCompute: force=true re-runs', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    const r1 = w.runInheritanceCompute();
    const r2 = w.runInheritanceCompute({ force: true });
    if (r1 === r2) throw new Error('force=true should re-run, got same identity');
  });

  // ============================================================
  // invalidateInheritanceCache
  // ============================================================
  t('invalidate: clears cache', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    w.runInheritanceCompute();
    if (!w.state.inheritanceResult) throw new Error('cache should be populated');
    w.invalidateInheritanceCache();
    if (w.state.inheritanceResult !== null) throw new Error('cache should be cleared');
  });

  // ============================================================
  // Cache key changes when candidate set changes
  // ============================================================
  t('cache key: changes when candidate added', () => {
    const labA = [], labB = [], labC = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) {
      labA.push(k); labB.push(k); labC.push(k);
    }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    const items1 = w._gatherActiveCandidatesForInheritance();
    const key1 = w._inheritanceCacheKey(items1, 'default');
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels: labC },
    ]);
    const items2 = w._gatherActiveCandidatesForInheritance();
    const key2 = w._inheritanceCacheKey(items2, 'default');
    if (key1 === key2) throw new Error('keys should differ');
  });

  t('cache key: differs between default and detailed mode', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const items = w._gatherActiveCandidatesForInheritance();
    const k1 = w._inheritanceCacheKey(items, 'default');
    const k2 = w._inheritanceCacheKey(items, 'detailed');
    if (k1 === k2) throw new Error('mode should be in cache key');
  });

  // ============================================================
  // _drawInheritanceLabelsStrip — exercise without crashing
  // ============================================================
  // Stub canvas context — JSDOM doesn't ship canvas backend
  function makeStubCtx() {
    return {
      save() {}, restore() {},
      fillRect() {}, strokeRect() {}, fillText() {},
      beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
      measureText(s) { return { width: s.length * 6 }; },
      _set: function(k, v) { this[k] = v; },
      get fillStyle() { return this._fillStyle; },
      set fillStyle(v) { this._fillStyle = v; },
      get strokeStyle() { return this._strokeStyle; },
      set strokeStyle(v) { this._strokeStyle = v; },
      get font() { return this._font; },
      set font(v) { this._font = v; },
      get textAlign() { return this._textAlign; },
      set textAlign(v) { this._textAlign = v; },
      get textBaseline() { return this._textBaseline; },
      set textBaseline(v) { this._textBaseline = v; },
      get lineWidth() { return this._lineWidth; },
      set lineWidth(v) { this._lineWidth = v; },
    };
  }

  t('draw: exercises without throwing when there is data', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    w.runInheritanceCompute();   // populate cache
    const ctx = makeStubCtx();
    const pad = { l: 40, r: 20, t: 20, b: 30 };
    w._drawInheritanceLabelsStrip(ctx, pad, 740, 150, 0, 10);
  });

  t('draw: no-op when toggle is off', () => {
    setupState([]);
    w.state.linesInheritanceLabelsOn = false;
    const ctx = makeStubCtx();
    const pad = { l: 40, r: 20, t: 20, b: 30 };
    w._drawInheritanceLabelsStrip(ctx, pad, 740, 150, 0, 10);
  });

  t('draw: no-op when fewer than 2 candidates', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const ctx = makeStubCtx();
    const pad = { l: 40, r: 20, t: 20, b: 30 };
    w._drawInheritanceLabelsStrip(ctx, pad, 740, 150, 0, 10);
  });

  t('draw: schedules deferred compute when cache empty', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    eq(w.state.inheritanceResult, null);
    const ctx = makeStubCtx();
    const pad = { l: 40, r: 20, t: 20, b: 30 };
    w._drawInheritanceLabelsStrip(ctx, pad, 740, 150, 0, 10);
    if (!w.state._inheritanceComputeScheduled) {
      throw new Error('expected _inheritanceComputeScheduled to be true');
    }
  });

  // ============================================================
  // setLinesInheritanceLabelsOn
  // ============================================================
  t('toggle: setLinesInheritanceLabelsOn updates state', () => {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    w.setLinesInheritanceLabelsOn(false);
    eq(w.state.linesInheritanceLabelsOn, false);
    w.setLinesInheritanceLabelsOn(true);
    eq(w.state.linesInheritanceLabelsOn, true);
  });

  // ============================================================
  // Realistic 3-candidate scenario
  // ============================================================
  t('realistic: 3 candidates with mixed linkage produce expected labels', () => {
    const N = 226;
    const labA = [], labB = [], labC = [];
    const groups = [40, 38, 38, 38, 36, 36];
    for (let g = 0; g < 6; g++) {
      for (let i = 0; i < groups[g]; i++) {
        labA.push(g);
        labB.push(g);                  // B perfectly linked to A
        labC.push(Math.floor(g / 2));  // C lower-resolution version
      }
    }
    setupState([
      { id: 'cand1', K: 6, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'cand2', K: 6, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
      { id: 'cand3', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels: labC },
    ], 'default');
    const r = w.runInheritanceCompute();
    if (!r) throw new Error('expected result');
    // Cand1 and Cand2 each should have 6 inheritance groups; cand3 has 3
    eq(r.rtab.per_item_n_groups[0], 6);
    eq(r.rtab.per_item_n_groups[1], 6);
    eq(r.rtab.per_item_n_groups[2], 3);
    // seq_num assignments
    eq(r.items_meta[0].seq_num, 1);
    eq(r.items_meta[1].seq_num, 2);
    eq(r.items_meta[2].seq_num, 3);
    // Labels would render as I1·6g, I2·6g, I3·3g
    eq(w._formatInheritanceLabel(r.items_meta[0].seq_num, r.items_meta[0].seq_num, r.rtab.per_item_n_groups[0]), 'I1·6g');
    eq(w._formatInheritanceLabel(r.items_meta[2].seq_num, r.items_meta[2].seq_num, r.rtab.per_item_n_groups[2]), 'I3·3g');
  });

  return failures;
}

setTimeout(() => {
  const failures = run();
  if (failures.length > 0) {
    console.log(`\n${failures.length} test(s) failed`);
    failures.forEach(f => console.log(`  - ${f.name}: ${f.err}`));
    process.exit(1);
  } else {
    console.log('\nAll tests passed');
    process.exit(0);
  }
}, 500);
