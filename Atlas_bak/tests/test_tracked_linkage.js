// Tests for turn 2i — tracked-linkage projection.
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

function makeStubCtx() {
  const ctx = {
    save() {}, restore() {},
    fillRect() {}, strokeRect() {}, fillText() {},
    beginPath() {}, moveTo() {}, lineTo() {}, stroke() {},
    measureText(s) { return { width: s.length * 6 }; },
  };
  ['fillStyle','strokeStyle','font','textAlign','textBaseline','lineWidth'].forEach(k => {
    let v = ''; Object.defineProperty(ctx, k, { get: () => v, set: nv => { v = nv; } });
  });
  return ctx;
}

function run() {
  const failures = [];
  let testNum = 0;
  function t(name, fn) {
    testNum++;
    try { fn(); console.log(`  PASS [${testNum}] ${name}`); }
    catch (e) { failures.push({ name, err: e.message }); console.log(`  FAIL [${testNum}] ${name}: ${e.message}`); }
  }
  function eq(a, b, m) { if (a !== b) throw new Error(`${m||''} expected ${JSON.stringify(b)}, got ${JSON.stringify(a)}`); }
  function near(a, b, eps, m) { if (Math.abs(a - b) > eps) throw new Error(`${m||''} expected ${b} ±${eps}, got ${a}`); }

  if (typeof w.computeTrackedLinkageProjection !== 'function'
      || typeof w._drawTrackedLinkageStrip !== 'function'
      || typeof w._tlpInhGroupColor !== 'function') {
    console.log('  FAIL: turn 2i functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2i functions exposed');

  // Setup helper — same as test_inheritance_matrix
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
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
    _state.tracked = [];
  }

  // ============================================================
  // _tlpInhGroupColor
  // ============================================================
  t('color: null returns grey fallback', () => {
    const c = w._tlpInhGroupColor(null);
    eq(c, '#7a8398');
  });
  t('color: group 0 returns blue', () => {
    const c = w._tlpInhGroupColor(0);
    eq(c, '#4fa3ff');
  });
  t('color: group 1 returns orange', () => {
    const c = w._tlpInhGroupColor(1);
    eq(c, '#f5a524');
  });
  t('color: cycles for high group ids', () => {
    const c0 = w._tlpInhGroupColor(0);
    const c10 = w._tlpInhGroupColor(10);   // 10 % 10 = 0
    eq(c0, c10);
  });

  // ============================================================
  // computeTrackedLinkageProjection — basics
  // ============================================================
  t('projection: empty fishIdx returns n_total=0 + empty', () => {
    setupCandidates([]);
    const proj = w.computeTrackedLinkageProjection([]);
    eq(proj.n_total, 0);
    eq(proj.per_candidate.length, 0);
  });

  t('projection: < MIN_FISH returns empty per_candidate', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    const proj = w.computeTrackedLinkageProjection([0, 1]);   // only 2
    eq(proj.n_total, 2);
    eq(proj.per_candidate.length, 0);   // bailed out
  });

  t('projection: 100% pure selection returns purity=1.0', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    // Indices 0..29 are all in band 0
    const fishIdx = Array.from({ length: 30 }, (_, i) => i);
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.n_total, 30);
    eq(proj.per_candidate.length, 1);
    eq(proj.per_candidate[0].dominant_band, 0);
    eq(proj.per_candidate[0].purity, 1.0);
    eq(proj.per_candidate[0].dominant_count, 30);
  });

  t('projection: mixed selection reports dominant + purity correctly', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    // 25 fish from band 0 + 5 from band 1 → purity 25/30 ≈ 0.833
    const fishIdx = [];
    for (let i = 0; i < 25; i++) fishIdx.push(i);            // band 0
    for (let i = 30; i < 35; i++) fishIdx.push(i);            // band 1
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.per_candidate[0].dominant_band, 0);
    near(proj.per_candidate[0].purity, 25 / 30, 0.001);
    eq(proj.per_candidate[0].per_band_counts[0], 25);
    eq(proj.per_candidate[0].per_band_counts[1], 5);
    eq(proj.per_candidate[0].per_band_counts[2], 0);
  });

  t('projection: walks all candidates', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels: [...labels] },
    ]);
    const fishIdx = Array.from({ length: 10 }, (_, i) => i);   // 10 fish in band 0
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.per_candidate.length, 3);
    // Sorted by start_bp via _gatherActiveCandidatesForInheritance
    eq(proj.per_candidate[0].seq_num, 1);
    eq(proj.per_candidate[1].seq_num, 2);
    eq(proj.per_candidate[2].seq_num, 3);
    // All three should report band 0 dominant with purity 1.0 (since all candidates have identical labels)
    for (const r of proj.per_candidate) {
      eq(r.dominant_band, 0);
      eq(r.purity, 1.0);
    }
  });

  t('projection: includes inh_group_id when inheritance result exists', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    // Run actual inheritance compute
    w.runInheritanceCompute();
    const fishIdx = [];
    for (let i = 0; i < 10; i++) fishIdx.push(i);
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.per_candidate.length, 2);
    if (proj.per_candidate[0].inh_group_id == null) throw new Error('expected inh_group_id');
    if (typeof proj.per_candidate[0].inh_group_color !== 'string') throw new Error('expected color');
  });

  t('projection: handles fish indices outside any band gracefully', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    // Add a -1 (unassigned) and an out-of-range to the fishIdx
    const fishIdx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    // All in band 0, purity 1.0
    eq(proj.per_candidate[0].purity, 1.0);
  });

  t('projection: skips candidates where no lassoed fish hit any band', () => {
    // Build state with TWO candidates with disjoint label arrays — second
    // candidate's locked_labels has -1 for all our lassoed fish positions.
    const labelsA = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labelsA.push(k);
    const labelsB = new Array(90).fill(-1);   // all unassigned
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labelsA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labelsB },
    ]);
    const fishIdx = [0, 1, 2, 3, 4];
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.per_candidate.length, 1);   // only A; B had no hits
    eq(proj.per_candidate[0].candidate_id, 'A');
  });

  // ============================================================
  // _drawTrackedLinkageStrip — renderer integration
  // ============================================================
  t('draw: no-op when no tracked', () => {
    setupCandidates([]);
    w.state.tracked = [];
    const ctx = makeStubCtx();
    w._drawTrackedLinkageStrip(ctx, { l: 40, r: 20, t: 20, b: 30 }, 600, 200, 0, 100);
  });

  t('draw: no-op when linesTrackedLinkageOn is explicitly false', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    w.state.tracked = [0, 1, 2, 3, 4];
    w.state.linesTrackedLinkageOn = false;
    const ctx = makeStubCtx();
    w._drawTrackedLinkageStrip(ctx, { l: 40, r: 20, t: 20, b: 30 }, 600, 200, 0, 100);
    // No throw = pass
    w.state.linesTrackedLinkageOn = true;   // restore for later tests
  });

  t('draw: exercises with valid tracked + candidates', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    w.state.tracked = [0, 1, 2, 3, 4, 5];
    const ctx = makeStubCtx();
    w._drawTrackedLinkageStrip(ctx, { l: 40, r: 20, t: 20, b: 30 }, 600, 200, 0, 10);
  });

  // ============================================================
  // printTrackedLinkageTable — console output (smoke test)
  // ============================================================
  t('printTable: handles empty tracked gracefully', () => {
    setupCandidates([]);
    w.state.tracked = [];
    const r = w.printTrackedLinkageTable();
    eq(r, null);
  });

  t('printTable: returns projection object on data', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    w.state.tracked = [0, 1, 2, 3, 4];
    const r = w.printTrackedLinkageTable();
    if (!r || !r.per_candidate) throw new Error('expected projection object');
  });

  // ============================================================
  // Realistic LG28-style scenario
  // ============================================================
  t('end-to-end: lasso of 40 fish with 95% purity reports correctly', () => {
    const N = 226;
    const labels = [];
    const groups = [40, 38, 38, 38, 36, 36];
    for (let g = 0; g < 6; g++) for (let i = 0; i < groups[g]; i++) labels.push(g);
    // Cohort: fish 0-39 are in band 0, 40-77 in band 1, etc.
    setupCandidates([
      { id: 'LG28_cand_01', K: 6, start_bp: 15_000_000, end_bp: 18_000_000, labels },
    ]);
    // Replace the auto-generated samples with N=226
    w.state.data.samples = Array.from({ length: N }, (_, i) => ({ id: `S${i}` }));
    // Lasso 40 fish, 38 from band 0 + 2 from band 1 (mistake) → 38/40 = 95%
    const fishIdx = [];
    for (let i = 0; i < 38; i++) fishIdx.push(i);                    // band 0
    fishIdx.push(40); fishIdx.push(41);                                // band 1 (lassoed by mistake)
    const proj = w.computeTrackedLinkageProjection(fishIdx);
    eq(proj.n_total, 40);
    eq(proj.per_candidate[0].dominant_band, 0);
    near(proj.per_candidate[0].purity, 0.95, 0.001);
    eq(proj.per_candidate[0].per_band_counts[0], 38);
    eq(proj.per_candidate[0].per_band_counts[1], 2);
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
