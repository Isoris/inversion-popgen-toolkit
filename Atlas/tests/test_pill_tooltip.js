// Tests for turn 2p — inheritance pill tooltip
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

  if (typeof w._inhTooltipEnsureEl !== 'function'
      || typeof w._inhTooltipBuildHtml !== 'function'
      || typeof w._inhTooltipShow !== 'function'
      || typeof w._inhTooltipHide !== 'function'
      || typeof w._wireInheritancePillTooltip !== 'function') {
    console.log('  FAIL: turn 2p functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2p functions exposed');

  function setupTwoCands() {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.k = 3;
    _state.candidates = {};
    _state.candidates_detailed = {};
    _state.linesInheritanceLabelsOn = true;
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    _state.candidates['A'] = {
      id: 'A', chrom: 'TEST', K: 3, K_used: 3,
      locked_labels: labels, start_w: 0, end_w: 30,
      start_bp: 1_000_000, end_bp: 2_000_000, source: 'test',
    };
    _state.candidates['B'] = {
      id: 'B', chrom: 'TEST', K: 3, K_used: 3,
      locked_labels: labels.slice(), start_w: 0, end_w: 30,
      start_bp: 5_000_000, end_bp: 6_000_000, source: 'test',
    };
    _state.data = {
      chrom: 'TEST', n_windows: 30,
      samples: Array.from({ length: 90 }, (_, i) => ({ id: `S${i}` })),
    };
    _state._inhPillHitRegions = [];
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
    return _state;
  }

  // ============================================================
  // Tooltip element creation
  // ============================================================
  t('ensureEl: creates tooltip div on first call', () => {
    // Remove any previous
    const existing = w.document.getElementById('inhPillTooltip');
    if (existing) existing.remove();
    const tip = w._inhTooltipEnsureEl();
    if (!tip) throw new Error('no tooltip element returned');
    eq(tip.id, 'inhPillTooltip');
    if (!w.document.getElementById('inhPillTooltip')) throw new Error('tip not in DOM');
  });

  t('ensureEl: idempotent — second call returns same element', () => {
    const tip1 = w._inhTooltipEnsureEl();
    const tip2 = w._inhTooltipEnsureEl();
    if (tip1 !== tip2) throw new Error('second call should return same element');
    const all = w.document.querySelectorAll('#inhPillTooltip');
    eq(all.length, 1);
  });

  // ============================================================
  // Hit region recording during draw
  // ============================================================
  t('draw: records hit regions for visible pills', () => {
    const _state = setupTwoCands();
    w.runInheritanceCompute();
    // Simulate a draw — call the strip drawer directly
    const ctx = makeStubCtx();
    const pad = { l: 60, r: 30, t: 30, b: 50 };
    w._drawInheritanceLabelsStrip(ctx, pad, 600, 400, 0, 10);
    // Hit regions should have been recorded for both candidates
    if (!_state._inhPillHitRegions || _state._inhPillHitRegions.length < 2) {
      throw new Error(`expected ≥2 hit regions, got ${_state._inhPillHitRegions ? _state._inhPillHitRegions.length : 0}`);
    }
    const r0 = _state._inhPillHitRegions[0];
    if (!r0.candidate_id) throw new Error('hit region missing candidate_id');
    if (r0.seq_num !== 1 && r0.seq_num !== 2) throw new Error('seq_num not set');
    if (typeof r0.x !== 'number' || typeof r0.y !== 'number') throw new Error('xy not numeric');
    if (typeof r0.w !== 'number' || typeof r0.h !== 'number') throw new Error('wh not numeric');
  });

  t('draw: resets hit regions on each call', () => {
    const _state = setupTwoCands();
    w.runInheritanceCompute();
    const ctx = makeStubCtx();
    const pad = { l: 60, r: 30, t: 30, b: 50 };
    w._drawInheritanceLabelsStrip(ctx, pad, 600, 400, 0, 10);
    const firstCount = _state._inhPillHitRegions.length;
    // Pollute with junk
    _state._inhPillHitRegions.push({ x: 999, y: 999, w: 1, h: 1, candidate_id: 'JUNK' });
    // Redraw — should reset
    w._drawInheritanceLabelsStrip(ctx, pad, 600, 400, 0, 10);
    eq(_state._inhPillHitRegions.length, firstCount);
    for (const r of _state._inhPillHitRegions) {
      if (r.candidate_id === 'JUNK') throw new Error('stale hit region survived');
    }
  });

  // ============================================================
  // Tooltip HTML content
  // ============================================================
  t('buildHtml: shows candidate id and seq_num', () => {
    setupTwoCands();
    w.runInheritanceCompute();
    const html = w._inhTooltipBuildHtml({
      candidate_id: 'A', seq_num: 1, n_groups: 2,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
    });
    if (html.indexOf('I1') < 0) throw new Error('expected I1');
    if (html.indexOf('A') < 0) throw new Error('expected candidate id A');
    if (html.indexOf('1.00') < 0) throw new Error('expected start Mb');
    if (html.indexOf('2.00') < 0) throw new Error('expected end Mb');
  });

  t('buildHtml: shows per-band fish counts', () => {
    setupTwoCands();
    w.runInheritanceCompute();
    const html = w._inhTooltipBuildHtml({
      candidate_id: 'A', seq_num: 1, n_groups: 2,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
    });
    // 30 fish per band
    if (html.indexOf('30 fish') < 0) throw new Error('expected band size in tooltip');
  });

  t('buildHtml: shows group IDs', () => {
    setupTwoCands();
    w.runInheritanceCompute();
    const html = w._inhTooltipBuildHtml({
      candidate_id: 'A', seq_num: 1, n_groups: 2,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
    });
    if (html.indexOf('g0') < 0 && html.indexOf('g1') < 0) {
      throw new Error('expected group IDs (g0/g1)');
    }
  });

  t('buildHtml: handles missing inheritance result gracefully', () => {
    if (!w.state) w.state = {};
    w.state.inheritanceResult = null;
    const html = w._inhTooltipBuildHtml({
      candidate_id: 'X', seq_num: 5, n_groups: 1,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
    });
    if (html.indexOf('unavailable') < 0 && html.indexOf('I5') < 0) {
      throw new Error('expected fallback rendering');
    }
  });

  // ============================================================
  // Show/hide
  // ============================================================
  t('show: makes tooltip visible with content', () => {
    setupTwoCands();
    w.runInheritanceCompute();
    const hit = {
      candidate_id: 'A', seq_num: 1, n_groups: 2,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
      x: 100, y: 50, w: 30, h: 14,
    };
    w._inhTooltipShow(hit, 200, 100);
    const tip = w.document.getElementById('inhPillTooltip');
    eq(tip.style.display, 'block');
    if (tip.innerHTML.indexOf('I1') < 0) throw new Error('tip content missing');
  });

  t('hide: makes tooltip invisible', () => {
    w._inhTooltipShow({
      candidate_id: 'A', seq_num: 1, n_groups: 1,
      start_bp: 1, end_bp: 2, item_idx: 0,
    }, 100, 100);
    w._inhTooltipHide();
    const tip = w.document.getElementById('inhPillTooltip');
    eq(tip.style.display, 'none');
  });

  // ============================================================
  // Mousemove handler
  // ============================================================
  t('wire: idempotent — second call does not re-attach', () => {
    const cv = w.document.createElement('canvas');
    cv.width = 800; cv.height = 400;
    let attaches = 0;
    const origAdd = cv.addEventListener.bind(cv);
    cv.addEventListener = function(type, fn) { if (type === 'mousemove') attaches++; return origAdd(type, fn); };
    w._wireInheritancePillTooltip(cv);
    w._wireInheritancePillTooltip(cv);
    eq(attaches, 1);
  });

  t('mousemove: shows tooltip when over a pill', () => {
    const _state = setupTwoCands();
    w.runInheritanceCompute();
    // Manually set hit regions
    _state._inhPillHitRegions = [{
      x: 100, y: 20, w: 40, h: 14,
      candidate_id: 'A', seq_num: 1, n_groups: 2,
      start_bp: 1_000_000, end_bp: 2_000_000, item_idx: 0,
    }];
    const cv = w.document.createElement('canvas');
    cv.width = 800; cv.height = 400;
    cv.getBoundingClientRect = () => ({ left: 0, top: 0, right: 800, bottom: 400, width: 800, height: 400 });
    w.document.body.appendChild(cv);
    w._wireInheritancePillTooltip(cv);
    // Simulate mousemove inside the pill region
    const ev = new w.MouseEvent('mousemove', { clientX: 120, clientY: 27 });
    cv.dispatchEvent(ev);
    const tip = w.document.getElementById('inhPillTooltip');
    eq(tip.style.display, 'block');
    cv.remove();
  });

  t('mousemove: hides tooltip when off all pills', () => {
    const _state = setupTwoCands();
    _state._inhPillHitRegions = [{
      x: 100, y: 20, w: 40, h: 14,
      candidate_id: 'A', seq_num: 1, n_groups: 1,
      start_bp: 1, end_bp: 2, item_idx: 0,
    }];
    const cv = w.document.createElement('canvas');
    cv.width = 800; cv.height = 400;
    cv.getBoundingClientRect = () => ({ left: 0, top: 0, right: 800, bottom: 400, width: 800, height: 400 });
    w.document.body.appendChild(cv);
    w._wireInheritancePillTooltip(cv);
    // First go inside, then outside
    cv.dispatchEvent(new w.MouseEvent('mousemove', { clientX: 120, clientY: 27 }));
    cv.dispatchEvent(new w.MouseEvent('mousemove', { clientX: 500, clientY: 200 }));
    const tip = w.document.getElementById('inhPillTooltip');
    eq(tip.style.display, 'none');
    cv.remove();
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
