// Tests for turn 2k — annotation cockpit page (cursor mechanics + render)
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
    setTransform() {}, translate() {}, rotate() {},
    measureText(s) { return { width: s.length * 6 }; },
  };
  ['fillStyle','strokeStyle','font','textAlign','textBaseline','lineWidth','globalAlpha'].forEach(k => {
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

  if (typeof w.refreshAnnotationCockpit !== 'function'
      || typeof w._annoCockpitOnKey !== 'function'
      || typeof w._annoCockpitOnClick !== 'function'
      || typeof w._annoCockpitCandidateAtCursor !== 'function'
      || typeof w._annoCockpitChromExtent !== 'function') {
    console.log('  FAIL: turn 2k functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2k functions exposed');

  function setupChrom(candConfigs, opts) {
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
      chrom: 'TEST',
      n_windows: 30,
      chrom_len_bp: (opts && opts.chrom_len_bp) || 50_000_000,
      samples: Array.from({ length: nSamples }, (_, i) => ({ id: `S${i}` })),
      windows: (opts && opts.windows) || null,
    };
    _state.cockpitCursor = { mb: null, candidate_id: null };
    _state.tracked = [];
    _state.inheritanceResult = null;
  }

  // Stub canvas for testing
  function stubPage() {
    // Make sure page21 exists in DOM (it's in the original HTML)
    const canvas = w.document.getElementById('annoCockpitCanvas');
    if (canvas) canvas.getContext = () => makeStubCtx();
  }

  // ============================================================
  // _annoCockpitChromExtent
  // ============================================================
  t('extent: uses chrom_len_bp when present', () => {
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: [0,1,2] },
    ], { chrom_len_bp: 50_000_000 });
    const items = w._gatherActiveCandidatesForInheritance();
    const ext = w._annoCockpitChromExtent(items);
    eq(ext.mbMin, 0);
    eq(ext.mbMax, 50);
  });

  t('extent: falls back to candidate range when no chrom_len_bp', () => {
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: [0,1,2] },
      { id: 'B', K: 3, start_bp: 30_000_000, end_bp: 32_000_000, labels: [0,1,2] },
    ]);
    delete w.state.data.chrom_len_bp;
    const items = w._gatherActiveCandidatesForInheritance();
    const ext = w._annoCockpitChromExtent(items);
    eq(ext.mbMin, 0);
    if (ext.mbMax < 32) throw new Error(`expected mbMax >= 32, got ${ext.mbMax}`);
  });

  t('extent: returns null with no candidates and no data', () => {
    w.state = w.state || {};
    w.state.data = {};
    w.state.candidates = {};
    w.state.candidates_detailed = {};
    const items = [];
    const ext = w._annoCockpitChromExtent(items);
    eq(ext, null);
  });

  // ============================================================
  // _annoCockpitCandidateAtCursor
  // ============================================================
  t('cursor: returns candidate when cursor is inside its range', () => {
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: [0,1,2] },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [0,1,2] },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    const at = w._annoCockpitCandidateAtCursor(1.5, items);
    eq(at.id, 'A');
    const at2 = w._annoCockpitCandidateAtCursor(5.5, items);
    eq(at2.id, 'B');
  });

  t('cursor: returns null when cursor is outside any candidate range', () => {
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: [0,1,2] },
    ]);
    const items = w._gatherActiveCandidatesForInheritance();
    eq(w._annoCockpitCandidateAtCursor(10, items), null);
    eq(w._annoCockpitCandidateAtCursor(0.5, items), null);
  });

  t('cursor: returns null with null mb', () => {
    eq(w._annoCockpitCandidateAtCursor(null, []), null);
  });

  // ============================================================
  // refreshAnnotationCockpit (DOM presence)
  // ============================================================
  t('refresh: 0 candidates shows empty state', () => {
    setupChrom([]);
    stubPage();
    w.refreshAnnotationCockpit();
    const empty = w.document.getElementById('annotationCockpitEmpty');
    const body = w.document.getElementById('annotationCockpitBody');
    eq(empty.style.display, 'block');
    eq(body.style.display, 'none');
  });

  t('refresh: with candidates hides empty state, shows body', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    const empty = w.document.getElementById('annotationCockpitEmpty');
    const body = w.document.getElementById('annotationCockpitBody');
    eq(empty.style.display, 'none');
    eq(body.style.display, 'block');
  });

  // ============================================================
  // _annoCockpitOnKey: cursor movement
  // ============================================================
  t('key: ArrowRight first time initializes cursor to chrom midpoint', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    eq(w.state.cockpitCursor.mb, null);
    const ev = { key: 'ArrowRight', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    near(w.state.cockpitCursor.mb, 15, 1);   // chrom midpoint = 15 Mb
  });

  t('key: ArrowRight steps cursor forward', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { chrom_len_bp: 100_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 5;
    const ev = { key: 'ArrowRight', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    if (w.state.cockpitCursor.mb <= 5) throw new Error('cursor should advance');
  });

  t('key: ArrowLeft steps cursor backward', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { chrom_len_bp: 100_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 50;
    const ev = { key: 'ArrowLeft', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    if (w.state.cockpitCursor.mb >= 50) throw new Error('cursor should retreat');
  });

  t('key: Shift+ArrowRight jumps to next candidate boundary', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
      { id: 'B', K: 3, start_bp: 15_000_000, end_bp: 17_000_000, labels: [...labels] },
      { id: 'C', K: 3, start_bp: 25_000_000, end_bp: 27_000_000, labels: [...labels] },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 0;
    const ev = { key: 'ArrowRight', shiftKey: true, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    // Should jump to first boundary which is start_bp of A = 5
    near(w.state.cockpitCursor.mb, 5, 0.01);
    w._annoCockpitOnKey(ev);   // jump to next (end of A = 7)
    near(w.state.cockpitCursor.mb, 7, 0.01);
    w._annoCockpitOnKey(ev);   // jump to start of B = 15
    near(w.state.cockpitCursor.mb, 15, 0.01);
  });

  t('key: Shift+ArrowLeft jumps to previous boundary', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
      { id: 'B', K: 3, start_bp: 15_000_000, end_bp: 17_000_000, labels: [...labels] },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 20;
    const ev = { key: 'ArrowLeft', shiftKey: true, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    // Should jump back to end of B = 17
    near(w.state.cockpitCursor.mb, 17, 0.01);
  });

  t('key: Escape clears cursor', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { chrom_len_bp: 50_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 25;
    const ev = { key: 'Escape', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    eq(w.state.cockpitCursor.mb, null);
  });

  t('key: ArrowKey is ignored when no candidates exist', () => {
    setupChrom([]);
    stubPage();
    const ev = { key: 'ArrowRight', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);   // should not throw
    if (w.state.cockpitCursor.mb !== null) throw new Error('cursor should remain null');
  });

  // ============================================================
  // Cursor stays clamped to chrom extent
  // ============================================================
  t('key: ArrowRight clamps to mbMax', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 30;   // already at max
    const ev = { key: 'ArrowRight', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    if (w.state.cockpitCursor.mb > 30) throw new Error('should be clamped to mbMax');
  });

  t('key: ArrowLeft clamps to mbMin', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 0;
    const ev = { key: 'ArrowLeft', shiftKey: false, preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    if (w.state.cockpitCursor.mb < 0) throw new Error('should be clamped to mbMin');
  });

  // ============================================================
  // Readouts
  // ============================================================
  t('readouts: cohort count appears', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    const cohortReadout = w.document.getElementById('annoCohortReadout');
    if (cohortReadout.textContent.indexOf('2 candidate') < 0) {
      throw new Error(`expected '2 candidate' in '${cohortReadout.textContent}'`);
    }
  });

  t('readouts: cursor position appears with candidate info', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitUpdateReadouts();
    const cursorReadout = w.document.getElementById('annoCursorReadout');
    if (cursorReadout.innerHTML.indexOf('I1') < 0) throw new Error('expected I1 in readout');
    if (cursorReadout.innerHTML.indexOf('6.00 Mb') < 0) throw new Error('expected mb position in readout');
  });

  t('readouts: between-candidates state shows italic placeholder', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 10_000_000, end_bp: 11_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 5;   // between A and B
    w._annoCockpitUpdateReadouts();
    const cursorReadout = w.document.getElementById('annoCursorReadout');
    if (cursorReadout.innerHTML.indexOf('between candidates') < 0) {
      throw new Error('expected between-candidates label');
    }
  });

  t('readouts: candidate info pane shows band sizes (turn 2l: as clickable buttons)', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    // Force re-render of the hap panel by clearing the lastHapCandId
    w.state._cockpitLastHapCandId = null;
    w._annoCockpitUpdateReadouts();
    const info = w.document.getElementById('annoCandidateInfo');
    // Turn 2l renders band sizes as buttons "0: 30", "1: 30", "2: 30"
    if (info.innerHTML.indexOf('0: 30') < 0) throw new Error('expected 0:30');
    if (info.innerHTML.indexOf('1: 30') < 0) throw new Error('expected 1:30');
    if (info.innerHTML.indexOf('2: 30') < 0) throw new Error('expected 2:30');
  });

  // ============================================================
  // _annoCockpitOnClick — sets cursor at the clicked mb
  // ============================================================
  t('click: sets cursor at the clicked x position', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    const canvas = w.document.getElementById('annoCockpitCanvas');
    // Mock getBoundingClientRect
    canvas.getBoundingClientRect = () => ({ left: 0, top: 0, right: 1200, bottom: 500, width: 1200, height: 500 });
    canvas.focus = () => {};
    // Click at x=600 (middle of canvas, so middle of plot area)
    const ev = {
      currentTarget: canvas,
      clientX: 600, clientY: 250,
    };
    w._annoCockpitOnClick(ev);
    // The cursor should be roughly at the chromosome midpoint (15 Mb)
    if (w.state.cockpitCursor.mb == null) throw new Error('cursor not set');
    near(w.state.cockpitCursor.mb, 15, 2);
  });

  t('click: outside plot area does not set cursor', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ], { chrom_len_bp: 30_000_000 });
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 5;   // pre-set
    const canvas = w.document.getElementById('annoCockpitCanvas');
    canvas.getBoundingClientRect = () => ({ left: 0, top: 0, right: 1200, bottom: 500, width: 1200, height: 500 });
    canvas.focus = () => {};
    const ev = {
      currentTarget: canvas,
      clientX: 10, clientY: 250,   // x=10 is in the left padding (60px)
    };
    w._annoCockpitOnClick(ev);
    eq(w.state.cockpitCursor.mb, 5);   // unchanged
  });

  // ============================================================
  // End-to-end: tab activation populates the cockpit
  // ============================================================
  t('end-to-end: tab21 button exists and activation works', () => {
    const btn = w.document.querySelector('#tabBar button[data-page="page21"]');
    if (!btn) throw new Error('page21 button missing from tab bar');
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    stubPage();
    // Simulate tab click via direct refresh
    w.refreshAnnotationCockpit();
    const body = w.document.getElementById('annotationCockpitBody');
    eq(body.style.display, 'block');
    const cohortReadout = w.document.getElementById('annoCohortReadout');
    if (cohortReadout.textContent.indexOf('2 candidate') < 0) {
      throw new Error('cohort readout did not update');
    }
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
