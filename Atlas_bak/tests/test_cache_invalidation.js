// Tests for turn 2n — auto cache invalidation on candidate edits
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
  url: 'http://localhost/',
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

  if (typeof w._hashLockedLabels !== 'function'
      || typeof w._inheritanceCacheKey !== 'function'
      || typeof w.runInheritanceCompute !== 'function') {
    console.log('  FAIL: turn 2n functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2n functions exposed');

  function setupCandidates(candConfigs) {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.k = 3;
    _state.candidates = {};
    _state.candidates_detailed = {};
    let nSamples = 0;
    for (const c of candConfigs) {
      _state.candidates[c.id] = {
        id: c.id, chrom: 'TEST',
        K: c.K, K_used: c.K,
        locked_labels: c.labels,
        start_w: 0, end_w: 30,
        start_bp: c.start_bp,
        end_bp: c.end_bp,
        source: 'test',
      };
      if (c.labels && c.labels.length > nSamples) nSamples = c.labels.length;
    }
    _state.data = {
      chrom: 'TEST', n_windows: 30,
      samples: Array.from({ length: nSamples }, (_, i) => ({ id: `S${i}` })),
    };
    _state.tracked = [];
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
    _state._inheritanceComputeScheduled = false;
  }

  // ============================================================
  // _hashLockedLabels
  // ============================================================
  t('hash: empty/null labels returns 0', () => {
    eq(w._hashLockedLabels(null), 0);
    eq(w._hashLockedLabels([]), 0);
  });

  t('hash: same labels return same hash', () => {
    const a = [0, 1, 2, 0, 1, 2, 0, 1, 2];
    const b = [0, 1, 2, 0, 1, 2, 0, 1, 2];
    eq(w._hashLockedLabels(a), w._hashLockedLabels(b));
  });

  t('hash: different labels return different hashes', () => {
    const a = [0, 1, 2, 0, 1, 2];
    const b = [0, 1, 2, 0, 1, 1];   // last element changed
    if (w._hashLockedLabels(a) === w._hashLockedLabels(b)) {
      throw new Error('hashes should differ');
    }
  });

  t('hash: handles -1 (unassigned) correctly', () => {
    const a = [0, 1, 2, -1, 0];
    const b = [0, 1, 2, 0, 0];
    if (w._hashLockedLabels(a) === w._hashLockedLabels(b)) {
      throw new Error('-1 should hash differently from 0');
    }
  });

  t('hash: order-sensitive', () => {
    const a = [0, 1, 2, 3];
    const b = [3, 2, 1, 0];
    if (w._hashLockedLabels(a) === w._hashLockedLabels(b)) {
      throw new Error('order should matter');
    }
  });

  t('hash: works on Int8Array (production type)', () => {
    const a = new Int8Array([0, 1, 2, 0, 1, 2]);
    const b = new Int8Array([0, 1, 2, 0, 1, 2]);
    eq(w._hashLockedLabels(a), w._hashLockedLabels(b));
  });

  // ============================================================
  // Cache key changes when locked_labels change (the BUG fixed in 2n)
  // ============================================================
  t('cache key: changes when locked_labels change with same K', () => {
    const labels1 = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels1.push(k);
    const labels2 = labels1.slice();
    labels2[0] = 1;   // flip one
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: labels1 }]);
    const items1 = w._gatherActiveCandidatesForInheritance();
    const k1 = w._inheritanceCacheKey(items1, 'default');
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: labels2 }]);
    const items2 = w._gatherActiveCandidatesForInheritance();
    const k2 = w._inheritanceCacheKey(items2, 'default');
    if (k1 === k2) throw new Error(`cache keys should differ: ${k1} vs ${k2}`);
  });

  t('cache key: changes when start_bp shifts (boundary refinement)', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([{ id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels }]);
    const items1 = w._gatherActiveCandidatesForInheritance();
    const k1 = w._inheritanceCacheKey(items1, 'default');
    setupCandidates([{ id: 'A', K: 3, start_bp: 1_500_000, end_bp: 2_000_000, labels }]);
    const items2 = w._gatherActiveCandidatesForInheritance();
    const k2 = w._inheritanceCacheKey(items2, 'default');
    if (k1 === k2) throw new Error('cache key should change when start_bp shifts');
  });

  t('cache key: changes when end_bp shifts', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([{ id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels }]);
    const items1 = w._gatherActiveCandidatesForInheritance();
    const k1 = w._inheritanceCacheKey(items1, 'default');
    setupCandidates([{ id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_500_000, labels }]);
    const items2 = w._gatherActiveCandidatesForInheritance();
    const k2 = w._inheritanceCacheKey(items2, 'default');
    if (k1 === k2) throw new Error('cache key should change when end_bp shifts');
  });

  t('cache key: stays same when nothing meaningful changes', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const items1 = w._gatherActiveCandidatesForInheritance();
    const k1 = w._inheritanceCacheKey(items1, 'default');
    // Same setup again
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: labels.slice() }]);
    const items2 = w._gatherActiveCandidatesForInheritance();
    const k2 = w._inheritanceCacheKey(items2, 'default');
    eq(k1, k2);
  });

  // ============================================================
  // runInheritanceCompute respects new cache invalidation
  // ============================================================
  t('compute: re-runs when locked_labels changed but ids/Ks identical', () => {
    const labels1 = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels1.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: labels1 },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: labels1.slice() },
    ]);
    const r1 = w.runInheritanceCompute();
    if (!r1) throw new Error('expected initial result');
    // Mutate B's labels in place (simulate re-cluster with same K)
    const newLabels = labels1.slice();
    newLabels[0] = 1; newLabels[1] = 1; newLabels[30] = 0; newLabels[31] = 0;
    w.state.candidates['B'].locked_labels = newLabels;
    const r2 = w.runInheritanceCompute();
    if (r1 === r2) throw new Error('compute should have re-run after labels mutation');
  });

  t('compute: cache hit when nothing changed', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: labels.slice() },
    ]);
    const r1 = w.runInheritanceCompute();
    const r2 = w.runInheritanceCompute();   // no changes
    eq(r1, r2);   // same identity = cache hit
  });

  // ============================================================
  // computeTrackedLinkageProjection drops inh when stale
  // ============================================================
  t('linkage: stale inheritance result is treated as null', () => {
    const labels1 = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels1.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labels1 },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labels1.slice() },
    ]);
    w.runInheritanceCompute();
    // Inheritance is now cached. Verify a projection includes inh_group_id.
    const fishIdx = [];
    for (let i = 0; i < 10; i++) fishIdx.push(i);
    const proj1 = w.computeTrackedLinkageProjection(fishIdx);
    if (proj1.per_candidate[0].inh_group_id == null) {
      throw new Error('expected inh_group_id in fresh state');
    }
    // Now mutate A's labels (cache becomes stale but we don't recompute yet)
    const newLabels = labels1.slice();
    newLabels[0] = 2; newLabels[1] = 2;
    w.state.candidates['A'].locked_labels = newLabels;
    // The cached inheritanceCacheKey now mismatches the live items
    const proj2 = w.computeTrackedLinkageProjection(fishIdx);
    // Stale: inh_group_id should be null (degraded gracefully)
    if (proj2.per_candidate[0].inh_group_id != null) {
      throw new Error('expected inh_group_id null after stale cache');
    }
    // But purity is still computed
    if (proj2.per_candidate[0].purity == null) {
      throw new Error('purity should still compute');
    }
  });

  t('linkage: stale-cache schedules a deferred recompute', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: labels.slice() },
    ]);
    w.runInheritanceCompute();
    // Mutate
    const newLabels = labels.slice();
    newLabels[0] = 2;
    w.state.candidates['A'].locked_labels = newLabels;
    // Sanity: scheduled flag is false before
    eq(w.state._inheritanceComputeScheduled, false);
    // Trigger via the linkage projection
    const fishIdx = [0, 1, 2, 3, 4, 5];
    w.computeTrackedLinkageProjection(fishIdx);
    // Flag should now be set (deferred recompute scheduled)
    if (!w.state._inheritanceComputeScheduled) {
      throw new Error('expected _inheritanceComputeScheduled = true');
    }
  });

  t('linkage: fresh cache does not schedule recompute', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: labels.slice() },
    ]);
    w.runInheritanceCompute();
    eq(w.state._inheritanceComputeScheduled, false);
    const fishIdx = [0, 1, 2, 3, 4, 5];
    w.computeTrackedLinkageProjection(fishIdx);
    // Cache is fresh; should NOT have scheduled
    eq(w.state._inheritanceComputeScheduled, false);
  });

  // ============================================================
  // Realistic end-to-end: refinement → cache invalidation → fresh result
  // ============================================================
  t('end-to-end: refining a candidate boundary invalidates inheritance', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labels.slice() },
    ]);
    const r1 = w.runInheritanceCompute();
    if (!r1) throw new Error('initial compute failed');
    // Refine A's boundary (start_bp shifts left)
    w.state.candidates['A'].start_bp = 800_000;
    const r2 = w.runInheritanceCompute();
    if (r1 === r2) throw new Error('compute should have re-run after boundary shift');
    // The new result should reflect the new bp range
    eq(r2.items_meta[0].start_bp, 800_000);
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
