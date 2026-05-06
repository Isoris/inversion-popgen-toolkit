// =============================================================================
// turn 132 — Slice 3: cusum_theta hero renderer for page 12.
//
// Builds on Slices 1 (DOM scaffold) and 2 (visibility wiring). Adds
// _drawThCusumHero() which reads state.data.cusum_theta and paints
// #thCusumStripCanvas + #thCusumHistCanvas + updates the spread/concord
// badges.
//
// JSON contract (state.data.cusum_theta):
//   { schema_version, layer, chrom, source_script, range_bp, persample[] }
// Where persample[i] = { sample_id, cp_bp, strength, asymmetry, karyotype }.
//
// Per the walked-back lib design (chat 487c7f04), JSON carries only
// per-carrier observations. Atlas computes display-only spread badge
// descriptively. Concord badge reads cusum_concordance if present.
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
// 1. Source-level checks — function defined, JSON contract documented, wired
// ============================================================================
console.log('\n=== 1. Source-level checks ===');

ok('_drawThCusumHero function defined',
   /function _drawThCusumHero\(\)/.test(html));

ok('JSON contract documented (persample[] shape)',
   /JSON contract \(state\.data\.cusum_theta\)/.test(html));

ok('contract documents persample[].cp_bp / strength / asymmetry / karyotype',
   /persample[\s\S]{0,500}cp_bp[\s\S]{0,200}strength[\s\S]{0,200}asymmetry[\s\S]{0,200}karyotype/.test(html));

ok('contract documents range_bp { start, end }',
   /range_bp:\s*\{\s*start:[\s\S]{0,80}end:/.test(html));

ok('contract documents source_script: STEP_T05_theta_cusum.R',
   /STEP_T05_theta_cusum\.R/.test(html));

ok('renderer fires after _refreshThetaPiPanelVisibility in onDataLoad path',
   /_refreshThetaPiPanelVisibility[\s\S]{0,500}_drawThCusumHero/.test(html));

ok('renderer reads from state.data.cusum_theta',
   /_drawThCusumHero[\s\S]{0,1500}state\.data\.cusum_theta/.test(html));

ok('renderer reads from state.data.cusum_concordance for concord badge',
   /_drawThCusumHero[\s\S]{0,20000}state\.data\.cusum_concordance/.test(html));

ok('renderer is vm-safe (early-returns when document is undefined)',
   /_drawThCusumHero[\s\S]{0,400}typeof document === 'undefined'/.test(html));

ok('walked-back design comment present (observe first, no aggregated modes)',
   /walked-back lib contract[\s\S]{0,400}per-carrier observations/.test(html));

// ============================================================================
// 2. Behavioral — sandbox exec with synthetic cusum_theta data
// ============================================================================
console.log('\n=== 2. Behavioral — sandbox exec ===');

// Minimal canvas mock that records the calls we care about.
function makeCanvasMock(initialW, initialH) {
  const calls = [];
  const ctx = {
    clearRect: (...a) => calls.push(['clearRect', ...a]),
    fillRect:  (...a) => calls.push(['fillRect',  ...a]),
    beginPath: ()     => calls.push(['beginPath']),
    moveTo:    (...a) => calls.push(['moveTo',    ...a]),
    lineTo:    (...a) => calls.push(['lineTo',    ...a]),
    stroke:    ()     => calls.push(['stroke']),
    set fillStyle(v)  { calls.push(['fillStyle', v]); },
    get fillStyle()   { return ''; },
    set strokeStyle(v){ calls.push(['strokeStyle', v]); },
    get strokeStyle() { return ''; },
    set lineWidth(v)  { calls.push(['lineWidth', v]); },
    get lineWidth()   { return 1; },
  };
  return {
    clientWidth:  initialW,
    clientHeight: initialH,
    width:        initialW,
    height:       initialH,
    getContext: () => ctx,
    _calls: calls,
  };
}

function makeBadgeMock(id) {
  return {
    id,
    textContent: '',
    title: '',
    style: {},
  };
}

function buildSandbox(cusumTheta, cusumConcordance) {
  const elements = {
    thCusumStripCanvas: makeCanvasMock(800, 110),
    thCusumHistCanvas:  makeCanvasMock(800,  70),
    thCusumSpreadBadge: makeBadgeMock('thCusumSpreadBadge'),
    thCusumConcordBadge: makeBadgeMock('thCusumConcordBadge'),
  };
  return {
    document: {
      elements,
      getElementById(id) { return this.elements[id] || null; },
    },
    state: {
      data: {
        cusum_theta: cusumTheta,
        cusum_concordance: cusumConcordance,
      },
    },
  };
}

function runFnIn(sandbox) {
  const fnSrc = html.match(/function _drawThCusumHero\(\)[\s\S]*?\n\}\n/)[0];
  const ctx = vm.createContext({
    document: sandbox.document,
    state:    sandbox.state,
    Number:   Number,
    Array:    Array,
    Math:     Math,
    Infinity: Infinity,
  });
  vm.runInContext(fnSrc + '\n_drawThCusumHero();', ctx);
}

// 2a. No cusum_theta data → renderer no-ops (no crash, badges untouched)
{
  const sandbox = buildSandbox(undefined, undefined);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2a: no cusum_theta data does not crash', !crashed);
  ok('2a: no data → strip canvas not painted (no clearRect calls)',
     sandbox.document.elements.thCusumStripCanvas._calls.filter(c => c[0] === 'clearRect').length === 0);
  ok('2a: no data → spread badge unchanged',
     sandbox.document.elements.thCusumSpreadBadge.textContent === '');
}

// 2b. Empty persample[] array → no-op
{
  const sandbox = buildSandbox({
    schema_version: 1,
    layer: 'cusum_theta',
    chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample: [],
  }, undefined);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2b: empty persample does not crash', !crashed);
  ok('2b: empty persample → strip canvas not painted',
     sandbox.document.elements.thCusumStripCanvas._calls.filter(c => c[0] === 'clearRect').length === 0);
}

// 2c. Tight cluster (IQR < 100kb) → spread = "tight"
{
  const persample = [];
  // 60 carriers tightly clustered around 14.0 Mb (±20 kb)
  for (let i = 0; i < 60; i++) {
    persample.push({
      sample_id: 'S' + i,
      cp_bp: 14000000 + Math.round((i - 30) * 600),  // span ~ 36 kb total
      strength: 1.0 + (i % 5) * 0.1,
      asymmetry: i % 2 === 0 ? 1 : -1,
      karyotype: 'HET',
    });
  }
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample,
  }, undefined);
  runFnIn(sandbox);

  ok('2c: tight cluster → both canvases painted (clearRect called)',
     sandbox.document.elements.thCusumStripCanvas._calls.some(c => c[0] === 'clearRect') &&
     sandbox.document.elements.thCusumHistCanvas._calls.some(c => c[0] === 'clearRect'));

  ok('2c: tight cluster → strip painted with fillRect calls (per-carrier ticks)',
     sandbox.document.elements.thCusumStripCanvas._calls.filter(c => c[0] === 'fillRect').length >= 60);

  ok('2c: tight cluster → spread badge says "tight"',
     /^spread tight/.test(sandbox.document.elements.thCusumSpreadBadge.textContent),
     'got: ' + sandbox.document.elements.thCusumSpreadBadge.textContent);

  ok('2c: spread badge includes IQR in kb',
     /IQR \d+ kb/.test(sandbox.document.elements.thCusumSpreadBadge.textContent));
}

// 2d. Intermediate spread (IQR 100k–500k) → spread = "intermediate"
{
  const persample = [];
  // 60 carriers spread evenly over 14.0–14.3 Mb (IQR ~150 kb)
  for (let i = 0; i < 60; i++) {
    persample.push({
      sample_id: 'S' + i,
      cp_bp: 14000000 + i * 5000,  // 0..295 kb span; IQR ~ 150 kb
      strength: 1.0,
      asymmetry: 1,
      karyotype: 'HET',
    });
  }
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample,
  }, undefined);
  runFnIn(sandbox);

  ok('2d: intermediate spread → spread badge says "intermediate"',
     /^spread intermediate/.test(sandbox.document.elements.thCusumSpreadBadge.textContent),
     'got: ' + sandbox.document.elements.thCusumSpreadBadge.textContent);
}

// 2e. Ragged spread (IQR > 500k) → spread = "ragged"
{
  const persample = [];
  // 60 carriers spread over 1.5 Mb → IQR ~ 750 kb
  for (let i = 0; i < 60; i++) {
    persample.push({
      sample_id: 'S' + i,
      cp_bp: 14000000 + i * 25000,  // 0..1.475 Mb span
      strength: 1.0,
      asymmetry: 1,
      karyotype: 'HET',
    });
  }
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample,
  }, undefined);
  runFnIn(sandbox);

  ok('2e: ragged spread → spread badge says "ragged"',
     /^spread ragged/.test(sandbox.document.elements.thCusumSpreadBadge.textContent),
     'got: ' + sandbox.document.elements.thCusumSpreadBadge.textContent);
}

// 2f. <4 carriers → spread badge says "—" (cannot compute IQR)
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample: [
      { sample_id: 'A', cp_bp: 14e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'B', cp_bp: 14.1e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.2e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
    ],
  }, undefined);
  runFnIn(sandbox);

  ok('2f: <4 carriers → spread badge text is "spread —"',
     sandbox.document.elements.thCusumSpreadBadge.textContent === 'spread —',
     'got: ' + sandbox.document.elements.thCusumSpreadBadge.textContent);

  ok('2f: <4 carriers → spread badge title explains why',
     /Not enough carriers/.test(sandbox.document.elements.thCusumSpreadBadge.title));
}

// 2g. Concord badge: no cusum_concordance → "concord —"
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample: [
      { sample_id: 'A', cp_bp: 14.00e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'B', cp_bp: 14.01e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.02e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'D', cp_bp: 14.03e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
    ],
  }, undefined);
  runFnIn(sandbox);

  ok('2g: no concordance layer → concord badge "concord —"',
     sandbox.document.elements.thCusumConcordBadge.textContent === 'concord —');

  ok('2g: concord badge title mentions STEP_DC06',
     /STEP_DC06/.test(sandbox.document.elements.thCusumConcordBadge.title));
}

// 2h. Concord badge: with cusum_concordance → "concord N/M · X%"
{
  const persample = [];
  for (let i = 0; i < 50; i++) {
    persample.push({ sample_id: 'S' + i, cp_bp: 14e6 + i * 1000,
                     strength: 1.0, asymmetry: 1, karyotype: 'HET' });
  }
  const sandbox = buildSandbox(
    { schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
      range_bp: { start: 0, end: 30e6 }, persample },
    { theta_ghsl_within_50kb: 38, n_carriers_compared: 54 }
  );
  runFnIn(sandbox);

  ok('2h: with concordance → concord badge "concord 38/54 · 70%"',
     sandbox.document.elements.thCusumConcordBadge.textContent === 'concord 38/54 · 70%',
     'got: ' + sandbox.document.elements.thCusumConcordBadge.textContent);
}

// 2i. Karyotype-aware coloring: ticks use HOM_REF/HET/HOM_INV palette
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 13e6, end: 16e6 },
    persample: [
      { sample_id: 'A', cp_bp: 14.0e6, strength: 1.0, asymmetry: 1, karyotype: 'HOM_REF' },
      { sample_id: 'B', cp_bp: 14.1e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.2e6, strength: 1.0, asymmetry: 1, karyotype: 'HOM_INV' },
      { sample_id: 'D', cp_bp: 14.3e6, strength: 1.0, asymmetry: 1, karyotype: null },
    ],
  }, undefined);
  runFnIn(sandbox);

  // Inspect fillStyle assignments — should include HOM_REF blue, HET orange,
  // HOM_INV purple, and ungrouped grey.
  const fillStyles = sandbox.document.elements.thCusumStripCanvas._calls
    .filter(c => c[0] === 'fillStyle')
    .map(c => c[1]);

  ok('2i: HOM_REF carrier rendered with blue (#3a7dde)',
     fillStyles.includes('#3a7dde'));
  ok('2i: HET carrier rendered with orange (#d97842)',
     fillStyles.includes('#d97842'));
  ok('2i: HOM_INV carrier rendered with purple (#7c4ad9)',
     fillStyles.includes('#7c4ad9'));
  ok('2i: null karyotype rendered with grey (#9aa3ad)',
     fillStyles.includes('#9aa3ad'));
}

// 2j. range_bp absent → renderer infers from min/max(cp_bp) + 5% pad
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    // no range_bp
    persample: [
      { sample_id: 'A', cp_bp: 14.0e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'B', cp_bp: 14.5e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.2e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'D', cp_bp: 14.3e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
    ],
  }, undefined);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2j: missing range_bp → no crash (infers from data)', !crashed);
  ok('2j: missing range_bp → strip still painted',
     sandbox.document.elements.thCusumStripCanvas._calls.some(c => c[0] === 'clearRect'));
}

// 2k. Idempotency — calling renderer twice yields same final state
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample: [
      { sample_id: 'A', cp_bp: 14.0e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'B', cp_bp: 14.1e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.2e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'D', cp_bp: 14.3e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
    ],
  }, undefined);
  runFnIn(sandbox);
  const t1 = sandbox.document.elements.thCusumSpreadBadge.textContent;
  runFnIn(sandbox);
  const t2 = sandbox.document.elements.thCusumSpreadBadge.textContent;
  ok('2k: idempotent — second draw yields same spread badge',
     t1 === t2, 't1=' + t1 + ' t2=' + t2);
}

// 2l. Defensive: persample with NaN cp_bp values are filtered out
{
  const sandbox = buildSandbox({
    schema_version: 1, layer: 'cusum_theta', chrom: 'LG28',
    range_bp: { start: 0, end: 30e6 },
    persample: [
      { sample_id: 'A', cp_bp: 14.0e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'B', cp_bp: NaN, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'C', cp_bp: 14.2e6, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
      { sample_id: 'D', cp_bp: undefined, strength: 1.0, asymmetry: 1, karyotype: 'HET' },
    ],
  }, undefined);
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('2l: NaN/undefined cp_bp values do not crash', !crashed);
}

// ============================================================================
// 3. Slice 1 + 2 still intact (no regressions in the panel scaffold)
// ============================================================================
console.log('\n=== 3. Slices 1 & 2 still intact ===');

ok('Slice 1: #thCusumHeroPanel still in DOM',
   /id="thCusumHeroPanel"/.test(html));
ok('Slice 1: #thCusumStripCanvas still in DOM',
   /id="thCusumStripCanvas"/.test(html));
ok('Slice 1: #thCusumHistCanvas still in DOM',
   /id="thCusumHistCanvas"/.test(html));
ok('Slice 2: _refreshThetaPiPanelVisibility still defined',
   /function _refreshThetaPiPanelVisibility\(\)/.test(html));
ok('Slice 2: cusum_theta layer still gates #thCusumHeroPanel visibility',
   /showHide\('thCusumHeroPanel',\s*hasCusumTheta/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
