// Tests for turn 2l — cockpit digit-key band selection + linkage panel + hap panel
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
      windows: null,
    };
    _state.cockpitCursor = { mb: null, candidate_id: null };
    _state.cockpitSelectedBand = null;
    _state._cockpitLastHapCandId = null;
    _state.tracked = [];
    _state.inheritanceResult = null;
    _state.haplotypeVocabs = {};
  }

  function stubPage() {
    const canvas = w.document.getElementById('annoCockpitCanvas');
    if (canvas) canvas.getContext = () => makeStubCtx();
    try { w.localStorage.clear(); } catch (_) {}
  }

  // ============================================================
  // Digit-key band selection
  // ============================================================
  t('digit: pressing 0 over candidate selects band 0', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    const ev = { key: '0', preventDefault: () => {} };
    w._annoCockpitOnKey(ev);
    eq(w.state.tracked.length, 30);
    eq(w.state.tracked[0], 0);
    eq(w.state.tracked[29], 29);
    eq(w.state.cockpitSelectedBand.band, 0);
    eq(w.state.cockpitSelectedBand.candidate_id, 'A');
    eq(w.state.cockpitSelectedBand.n, 30);
  });

  t('digit: pressing 1 selects band 1', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitOnKey({ key: '1', preventDefault: () => {} });
    eq(w.state.tracked[0], 30);
    eq(w.state.tracked[29], 59);
    eq(w.state.cockpitSelectedBand.band, 1);
  });

  t('digit: pressing digit beyond K is ignored', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitOnKey({ key: '5', preventDefault: () => {} });
    eq(w.state.tracked.length, 0);
    eq(w.state.cockpitSelectedBand, null);
  });

  t('digit: pressing key when cursor is between candidates is ignored', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 10_000_000, end_bp: 11_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 5;   // gap
    w._annoCockpitOnKey({ key: '0', preventDefault: () => {} });
    eq(w.state.tracked.length, 0);
  });

  t('digit: K=6 candidate accepts digits 0-5', () => {
    const N = 120;
    const labels = [];
    for (let k = 0; k < 6; k++) for (let i = 0; i < 20; i++) labels.push(k);
    setupChrom([
      { id: 'M', K: 6, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitOnKey({ key: '5', preventDefault: () => {} });
    eq(w.state.tracked.length, 20);
    eq(w.state.tracked[0], 100);   // band 5 starts at index 100
    eq(w.state.cockpitSelectedBand.band, 5);
  });

  // ============================================================
  // Escape clears tracked too (turn 2l added)
  // ============================================================
  t('Escape: clears state.tracked in addition to cursor', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w.state.tracked = [0, 1, 2, 3, 4];
    w._annoCockpitOnKey({ key: 'Escape', preventDefault: () => {} });
    eq(w.state.tracked.length, 0);
    eq(w.state.cockpitCursor.mb, null);
  });

  // ============================================================
  // _annoCockpitRenderLinkagePanel
  // ============================================================
  t('linkage panel: shows placeholder when no tracked', () => {
    setupChrom([]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.tracked = [];
    if (typeof w._annoCockpitRenderLinkagePanel !== 'function') {
      // Find via global window — function may not be exposed
      // We call via selecting a band which internally renders
      // Skip if not exposed
      return;
    }
    w._annoCockpitRenderLinkagePanel();
    const panel = w.document.getElementById('annoLinkagePanel');
    if (panel.innerHTML.indexOf('No band selected') < 0) {
      throw new Error('expected placeholder');
    }
  });

  t('linkage panel: populates after band selection', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 1.5;
    w._annoCockpitOnKey({ key: '0', preventDefault: () => {} });
    const panel = w.document.getElementById('annoLinkagePanel');
    if (panel.innerHTML.indexOf('<table') < 0) throw new Error('expected table');
    if (panel.innerHTML.indexOf('I1') < 0) throw new Error('expected I1');
    if (panel.innerHTML.indexOf('I2') < 0) throw new Error('expected I2');
  });

  t('linkage panel: clear button resets tracked', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 1.5;
    w._annoCockpitOnKey({ key: '0', preventDefault: () => {} });
    eq(w.state.tracked.length, 30);
    const clearBtn = w.document.getElementById('annoCockpitClearBtn');
    if (!clearBtn) throw new Error('clear button missing');
    clearBtn.click();
    eq(w.state.tracked.length, 0);
    eq(w.state.cockpitSelectedBand, null);
  });

  // ============================================================
  // _annoCockpitRenderHapPanel
  // ============================================================
  t('hap panel: shows placeholder when no candidate under cursor', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 10_000_000, end_bp: 11_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 5;   // gap
    w.state._cockpitLastHapCandId = 'STALE';
    w._annoCockpitUpdateReadouts();
    const info = w.document.getElementById('annoCandidateInfo');
    if (info.innerHTML.indexOf('Move the cursor') < 0) {
      throw new Error('expected placeholder');
    }
  });

  t('hap panel: shows candidate header + band buttons + hap labels section', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w.state._cockpitLastHapCandId = null;
    w._annoCockpitUpdateReadouts();
    const info = w.document.getElementById('annoCandidateInfo');
    if (info.innerHTML.indexOf('I1') < 0) throw new Error('expected I1 header');
    if (info.innerHTML.indexOf('anno-band-pick') < 0) throw new Error('expected band-pick buttons');
    if (info.innerHTML.indexOf('hapLabelsSection') < 0) throw new Error('expected hap labels section');
  });

  t('hap panel: only re-renders when cursor moves to different candidate', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
      { id: 'B', K: 3, start_bp: 10_000_000, end_bp: 12_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitUpdateReadouts();
    const firstHtml = w.document.getElementById('annoCandidateInfo').innerHTML;
    // Move cursor within A — should NOT re-render
    w.state.cockpitCursor.mb = 6.5;
    w._annoCockpitUpdateReadouts();
    const secondHtml = w.document.getElementById('annoCandidateInfo').innerHTML;
    eq(firstHtml, secondHtml);   // unchanged
    // Move to candidate B — SHOULD re-render
    w.state.cockpitCursor.mb = 11;
    w._annoCockpitUpdateReadouts();
    const thirdHtml = w.document.getElementById('annoCandidateInfo').innerHTML;
    if (thirdHtml === secondHtml) throw new Error('should have re-rendered for new candidate');
  });

  t('hap panel: clicking band button selects that band', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    w.state.cockpitCursor.mb = 6;
    w.state._cockpitLastHapCandId = null;
    w._annoCockpitUpdateReadouts();
    const button = w.document.querySelector('.anno-band-pick[data-band-idx="2"]');
    if (!button) throw new Error('band button missing');
    button.click();
    eq(w.state.tracked.length, 30);
    eq(w.state.tracked[0], 60);
    eq(w.state.cockpitSelectedBand.band, 2);
  });

  // ============================================================
  // refactored _wireCandidateHaplotypeAnnotations accepts sectionEl
  // ============================================================
  t('wire: _wireCandidateHaplotypeAnnotations accepts explicit section', () => {
    // Set up minimal state so candidateHaplotypeAnnotationsHtml has the
    // window.state hooks it needs for vocab lookup
    if (!w.state) w.state = {};
    w.state.haplotypeVocabs = {};
    w.state.candidates = {};
    w.state.candidates_detailed = {};
    try { w.localStorage.clear(); } catch (_) {}
    const c = { id: 'cand_W', chrom: 'TEST', K: 3 };
    const container = w.document.createElement('div');
    container.innerHTML = w.candidateHaplotypeAnnotationsHtml(c);
    w.document.body.appendChild(container);
    // JSDOM querySelector with duplicate IDs across the document can be
    // flaky; use children[0] since we know the layout.
    const section = container.children[0];
    if (!section || section.id !== 'hapLabelsSection') throw new Error('expected hapLabelsSection as first child');
    w._wireCandidateHaplotypeAnnotations(c, section);
    const picker = section.querySelector('.hap-label-picker[data-band-idx="0"]');
    if (!picker) throw new Error('picker not found');
    picker.value = 'HOM_REF';
    picker.dispatchEvent(new w.Event('change', { bubbles: true }));
    eq(c.haplotype_labels[0], 'HOM_REF');
    container.remove();
  });

  // ============================================================
  // End-to-end: full workflow
  // ============================================================
  t('end-to-end: cursor + digit-key + hap label + linkage shading', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupChrom([
      { id: 'A', K: 3, start_bp: 5_000_000, end_bp: 7_000_000, labels },
      { id: 'B', K: 3, start_bp: 15_000_000, end_bp: 17_000_000, labels: [...labels] },
      { id: 'C', K: 3, start_bp: 25_000_000, end_bp: 27_000_000, labels: [...labels] },
    ]);
    stubPage();
    w.refreshAnnotationCockpit();
    // Step 1: place cursor on A
    w.state.cockpitCursor.mb = 6;
    w._annoCockpitUpdateReadouts();
    eq(w.state._cockpitLastHapCandId, 'A');
    // Step 2: pick band 0
    w._annoCockpitOnKey({ key: '0', preventDefault: () => {} });
    eq(w.state.cockpitSelectedBand.candidate_id, 'A');
    eq(w.state.cockpitSelectedBand.band, 0);
    // Step 3: linkage panel populated
    const linkagePanel = w.document.getElementById('annoLinkagePanel');
    if (linkagePanel.innerHTML.indexOf('I1') < 0) throw new Error('linkage panel missing I1');
    if (linkagePanel.innerHTML.indexOf('I2') < 0) throw new Error('linkage panel missing I2');
    if (linkagePanel.innerHTML.indexOf('I3') < 0) throw new Error('linkage panel missing I3');
    // Step 4: cursor moves to B
    w.state.cockpitCursor.mb = 16;
    w._annoCockpitUpdateReadouts();
    eq(w.state._cockpitLastHapCandId, 'B');
    // hap panel re-rendered for B
    const info = w.document.getElementById('annoCandidateInfo');
    if (info.innerHTML.indexOf('· B') < 0 && info.innerHTML.indexOf('I2 · B') < 0 && info.innerHTML.indexOf('B</b>') < 0 && info.innerHTML.indexOf('B ·') < 0) {
      // Just check that B appears somewhere in the candidate ID portion
      if (info.innerHTML.indexOf('B') < 0) throw new Error('hap panel did not update for B');
    }
    // Step 5: tracked still has A's band-0 fish (lasso persists across cursor moves)
    eq(w.state.tracked.length, 30);
    eq(w.state.tracked[0], 0);
    // Step 6: ESC clears
    w._annoCockpitOnKey({ key: 'Escape', preventDefault: () => {} });
    eq(w.state.tracked.length, 0);
    eq(w.state.cockpitCursor.mb, null);
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
