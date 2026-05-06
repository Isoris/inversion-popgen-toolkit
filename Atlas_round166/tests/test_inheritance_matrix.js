// Tests for turn 2h — heatmap matrix UI (Cramér's V cell grid + detail panel).
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
    translate() {}, rotate() {},
    measureText(s) { return { width: s.length * 6 }; },
    createLinearGradient() {
      return { addColorStop() {} };
    },
  };
  // Property setters/getters
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

  if (typeof w.renderInheritanceMatrix !== 'function'
      || typeof w.openInheritanceMatrix !== 'function'
      || typeof w.closeInheritanceMatrix !== 'function'
      || typeof w._renderInheritanceCellDetail !== 'function'
      || typeof w._imxHeatColor !== 'function') {
    console.log('  FAIL: turn 2h functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2h functions exposed');

  // Stub canvas getContext for JSDOM
  function stubCanvas(canvas) {
    canvas.getContext = () => makeStubCtx();
    return canvas;
  }

  function clearLS() { try { w.localStorage.clear(); } catch (_) {} }

  function setupCandidates(candConfigs, opts) {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = (opts && opts.mode) || 'default';
    _state.k = (_state.activeMode === 'detailed') ? 6 : 3;
    const target = (_state.activeMode === 'detailed') ? 'candidates_detailed' : 'candidates';
    _state[target] = {};
    let nSamples = 0;
    for (const c of candConfigs) {
      _state[target][c.id] = {
        id: c.id,
        chrom: c.chrom || 'TEST',
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
      chrom: 'TEST',
      n_windows: 30,
      samples: Array.from({ length: nSamples }, (_, i) => ({ id: `S${i}` })),
    };
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
  }

  // ============================================================
  // _imxHeatColor
  // ============================================================
  t('heat color: V=0 returns valid hsl', () => {
    const c = w._imxHeatColor(0);
    if (typeof c !== 'string' || !c.startsWith('hsl(')) throw new Error(`bad color: ${c}`);
  });
  t('heat color: V=1 returns valid hsl', () => {
    const c = w._imxHeatColor(1);
    if (typeof c !== 'string' || !c.startsWith('hsl(')) throw new Error(`bad color: ${c}`);
  });
  t('heat color: V=0.5 returns valid hsl', () => {
    const c = w._imxHeatColor(0.5);
    if (typeof c !== 'string' || !c.startsWith('hsl(')) throw new Error(`bad color: ${c}`);
  });
  t('heat color: out-of-range clamped', () => {
    const c1 = w._imxHeatColor(-0.5);
    const c2 = w._imxHeatColor(1.5);
    if (typeof c1 !== 'string' || typeof c2 !== 'string') throw new Error('clamp failed');
  });

  // ============================================================
  // renderInheritanceMatrix — empty / single / multi candidate cases
  // ============================================================
  t('render: 0 candidates -> empty placeholder', () => {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.candidates = {};
    _state.candidates_detailed = {};
    const canvas = w.document.createElement('canvas');
    canvas.width = 640; canvas.height = 480;
    stubCanvas(canvas);
    const r = w.renderInheritanceMatrix(canvas);
    if (!r || r.items.length !== 0) throw new Error('expected 0 items');
  });

  t('render: 1 candidate -> minimum-2-needed placeholder', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    const canvas = w.document.createElement('canvas');
    canvas.width = 640; canvas.height = 480;
    stubCanvas(canvas);
    const r = w.renderInheritanceMatrix(canvas);
    eq(r.items.length, 1);
    eq(r.matrix, null);   // returned null when N<2
  });

  t('render: 2 perfectly linked candidates -> matrix has cramer_v populated', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    const canvas = w.document.createElement('canvas');
    canvas.width = 640; canvas.height = 480;
    stubCanvas(canvas);
    const r = w.renderInheritanceMatrix(canvas);
    eq(r.items.length, 2);
    if (!r.matrix) throw new Error('matrix missing');
    if (!r.matrix.cramer_v) throw new Error('cramer_v missing');
    // Off-diagonal V should be 1.0 (perfect linkage)
    const v = r.matrix.cramer_v[0 * 2 + 1];
    if (!(v > 0.95)) throw new Error(`expected V>0.95 for perfect linkage, got ${v}`);
    if (!r.hitTest) throw new Error('hitTest missing');
  });

  t('render: hitTest resolves canvas (x,y) to (i,j)', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
      { id: 'C', K: 3, start_bp: 9, end_bp: 10, labels: [...labels] },
    ]);
    const canvas = w.document.createElement('canvas');
    canvas.width = 640; canvas.height = 480;
    stubCanvas(canvas);
    const r = w.renderInheritanceMatrix(canvas);
    if (!r.hitTest) throw new Error('hitTest missing');
    // Click in upper-left first cell
    const ox = r.origin.ox, oy = r.origin.oy, cellPx = r.cellPx;
    const cell0 = r.hitTest(ox + cellPx / 2, oy + cellPx / 2);
    eq(cell0.i, 0); eq(cell0.j, 0);
    // Click in (1, 2) cell
    const cell12 = r.hitTest(ox + 2 * cellPx + 5, oy + 1 * cellPx + 5);
    eq(cell12.i, 1); eq(cell12.j, 2);
    // Click outside the matrix
    const out = r.hitTest(5, 5);
    eq(out, null);
  });

  // ============================================================
  // _renderInheritanceCellDetail
  // ============================================================
  t('detail: diagonal cell shows "no information" message', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    const panel = w.document.createElement('div');
    w._renderInheritanceCellDetail(panel, items, null, 0, 0);
    if (panel.innerHTML.indexOf('Diagonal') < 0) throw new Error('expected Diagonal text');
  });

  t('detail: off-diagonal cell renders contingency table + interpretation', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    const mat = w.crossCandidateMatrix(items);
    const panel = w.document.createElement('div');
    w._renderInheritanceCellDetail(panel, items, mat, 0, 1);
    // Should contain Cramér's V text
    if (panel.innerHTML.indexOf("Cramér's V") < 0) throw new Error('missing V text');
    // Should contain contingency table
    if (panel.innerHTML.indexOf('<table') < 0) throw new Error('missing table');
    // Strong linkage interpretation should appear (V close to 1)
    if (panel.innerHTML.indexOf('strong linkage') < 0) throw new Error('missing strong linkage');
  });

  t('detail: handles null panel gracefully', () => {
    w._renderInheritanceCellDetail(null, [], null, 0, 1);
  });

  t('detail: weak linkage (independent candidates) shows correct interpretation', () => {
    // Make two completely independent random partitions — V will be near 0
    const labelsA = [];
    const labelsB = [];
    for (let i = 0; i < 90; i++) {
      labelsA.push(i % 3);
      labelsB.push(Math.floor(i / 30));   // perfectly anti-correlated with A's pattern
    }
    // Shuffle B so it becomes independent of A
    let seed = 7;
    function rng() { seed = (seed * 9301 + 49297) % 233280; return seed / 233280; }
    for (let i = labelsB.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      [labelsB[i], labelsB[j]] = [labelsB[j], labelsB[i]];
    }
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels: labelsA },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: labelsB },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    const mat = w.crossCandidateMatrix(items);
    const panel = w.document.createElement('div');
    w._renderInheritanceCellDetail(panel, items, mat, 0, 1);
    // The interpretation should not be "strong linkage" for shuffled labels
    if (panel.innerHTML.indexOf('strong linkage') >= 0) {
      throw new Error('shuffled labels should not show strong linkage');
    }
  });

  // ============================================================
  // openInheritanceMatrix / closeInheritanceMatrix
  // ============================================================
  t('open: creates modal in DOM if missing', () => {
    // Ensure no modal exists
    const existing = w.document.getElementById('inheritanceMatrixModal');
    if (existing) existing.remove();

    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    // Stub canvas getContext on document.createElement
    const origCreateElement = w.document.createElement.bind(w.document);
    w.document.createElement = function(tag) {
      const el = origCreateElement(tag);
      if (tag === 'canvas') stubCanvas(el);
      return el;
    };

    w.openInheritanceMatrix();
    const modal = w.document.getElementById('inheritanceMatrixModal');
    if (!modal) throw new Error('modal not created');
    if (modal.style.display !== 'flex') throw new Error('modal not visible');

    w.document.createElement = origCreateElement;
  });

  t('close: hides modal', () => {
    w.closeInheritanceMatrix();
    const modal = w.document.getElementById('inheritanceMatrixModal');
    if (modal && modal.style.display !== 'none') throw new Error('modal not hidden');
  });

  t('open + close: idempotent', () => {
    const origCreateElement = w.document.createElement.bind(w.document);
    w.document.createElement = function(tag) {
      const el = origCreateElement(tag);
      if (tag === 'canvas') stubCanvas(el);
      return el;
    };
    w.openInheritanceMatrix();
    w.closeInheritanceMatrix();
    w.openInheritanceMatrix();
    const modal = w.document.getElementById('inheritanceMatrixModal');
    if (modal.style.display !== 'flex') throw new Error('reopen failed');
    w.closeInheritanceMatrix();
    w.document.createElement = origCreateElement;
  });

  // ============================================================
  // End-to-end: realistic 4-candidate scenario
  // ============================================================
  t('end-to-end: 4 candidates with mixed linkage produces matrix + detail', () => {
    const N = 90;
    // 2 linked candidates + 2 independent ones
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    // C and D are independent of A/B — shuffle
    const labC = labA.slice();
    const labD = labA.slice();
    let seed = 42;
    function rng() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
    for (let i = labC.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      [labC[i], labC[j]] = [labC[j], labC[i]];
      [labD[i], labD[j]] = [labD[j], labD[i]];
    }
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels: labC },
      { id: 'D', K: 3, start_bp: 13_000_000, end_bp: 14_000_000, labels: labD },
    ]);
    const canvas = w.document.createElement('canvas');
    canvas.width = 640; canvas.height = 480;
    stubCanvas(canvas);
    const r = w.renderInheritanceMatrix(canvas);
    eq(r.items.length, 4);
    if (!r.matrix) throw new Error('matrix null');
    // A-B linkage should be ~1.0 (V[0,1])
    const v_AB = r.matrix.cramer_v[0 * 4 + 1];
    if (!(v_AB > 0.95)) throw new Error(`expected V_AB>0.95, got ${v_AB}`);
    // A-C linkage should be much lower
    const v_AC = r.matrix.cramer_v[0 * 4 + 2];
    if (!(v_AC < v_AB)) throw new Error(`expected V_AC<V_AB, got V_AB=${v_AB} V_AC=${v_AC}`);
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
