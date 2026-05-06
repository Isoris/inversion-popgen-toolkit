// =============================================================================
// turn 156 — V-shape (u, agreement_fraction) per-candidate diagnostic
// =============================================================================
// Port of the FIG_C30 V-shape signature from STEP29's coherence pipeline.
// Per-sample agreement_fraction (from computeStripeQuality) plotted against
// u-axis position (from _getOrComputeUVRotation). Clean candidates show a V
// dip at u≈het_u with high homo plateaus.
//
// What this turn ships:
//   - _buildVShapeData(candidate, chunk, sqRows) — unifies UV rotation and
//     stripe-quality outputs into a (u, agreement, group) triples array.
//   - _vShapeColor(group) — group → hex color (reuses _BR_KARYO_COLORS).
//   - _drawVShapePlot(canvas, data, opts) — canvas renderer with axes,
//     gridlines, reference lines at het_u/hom1_u/hom2_u and agreement=0.5,
//     legend, group-colored dots.
//   - openVShapePlot(candidate) / closeVShapePlot() — modal popup with
//     graceful empty/missing-chunk/missing-sq states.
//   - "🔍 V-shape diagnostic" button in the G-panel karyotype tab footer.
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
// 1. Source-pattern: helpers + popup declared and exposed
// ============================================================================
console.log('\n=== 1. Helpers + popup declared ===');

ok('_buildVShapeData declared',
   /function\s+_buildVShapeData\s*\(\s*candidate\s*,\s*chunk\s*,\s*sqRows\s*\)/.test(html));

ok('_vShapeColor declared',
   /function\s+_vShapeColor\s*\(\s*group\s*\)/.test(html));

ok('_drawVShapePlot declared',
   /function\s+_drawVShapePlot\s*\(\s*canvas\s*,\s*data\s*,\s*opts\s*\)/.test(html));

ok('openVShapePlot declared',
   /function\s+openVShapePlot\s*\(\s*candidate\s*\)/.test(html));

ok('closeVShapePlot declared',
   /function\s+closeVShapePlot\s*\(\s*\)/.test(html));

ok('all 5 functions exposed on window',
   /window\._buildVShapeData\b/.test(html) &&
   /window\._vShapeColor\b/.test(html) &&
   /window\._drawVShapePlot\b/.test(html) &&
   /window\.openVShapePlot\b/.test(html) &&
   /window\.closeVShapePlot\b/.test(html));

// ============================================================================
// 2. Source-pattern: _buildVShapeData body
// ============================================================================
console.log('\n=== 2. _buildVShapeData body ===');

const buildBody = pullFunction(html, '_buildVShapeData');
ok('_buildVShapeData body extracted', !!buildBody);

ok('checks candidate.ref_l2 (returns MISSING_REF_L2 if absent)',
   /reason:\s*'MISSING_REF_L2'/.test(buildBody));

ok('calls _getOrComputeUVRotation(candidate.ref_l2)',
   /_getOrComputeUVRotation\s*\(\s*candidate\.ref_l2\s*\)/.test(buildBody));

ok('checks chunk presence (returns MISSING_CHUNK if absent)',
   /reason:\s*'MISSING_CHUNK'/.test(buildBody));

ok('checks sqRows is non-empty array (returns MISSING_SQ if absent)',
   /reason:\s*'MISSING_SQ'/.test(buildBody));

ok('builds idToCohort map by walking state.data.samples',
   /idToCohort/.test(buildBody) && /state\.data\.samples/.test(buildBody));

ok('returns reason=NO_POINTS when no valid points',
   /reason:\s*'NO_POINTS'/.test(buildBody));

ok('point shape: { u, agreement, group, sample_idx, cga, coherence_class, stripe_quality }',
   /u:\s*u/.test(buildBody) &&
   /agreement:\s*agr/.test(buildBody) &&
   /group:/.test(buildBody) &&
   /sample_idx:/.test(buildBody) &&
   /cga:/.test(buildBody));

ok('output carries het_u, hom1_u, hom2_u from rotation',
   /het_u:\s*rot\.het_u/.test(buildBody) &&
   /hom1_u:\s*rot\.hom1_u/.test(buildBody) &&
   /hom2_u:\s*rot\.hom2_u/.test(buildBody));

// ============================================================================
// 3. Source-pattern: _drawVShapePlot body
// ============================================================================
console.log('\n=== 3. _drawVShapePlot body ===');

const drawBody = pullFunction(html, '_drawVShapePlot');
ok('_drawVShapePlot body extracted', !!drawBody);

ok('returns reason from data when !data.ok',
   /!data\.ok/.test(drawBody) && /reason:/.test(drawBody));

ok('handles missing canvas context (returns NO_CTX)',
   /reason:\s*'NO_CTX'/.test(drawBody));

ok('respects devicePixelRatio for crisp rendering',
   /devicePixelRatio/.test(drawBody) && /setTransform/.test(drawBody));

ok('default canvas size 560 × 360',
   /560/.test(drawBody) && /360/.test(drawBody));

ok('draws horizontal gridlines at agreement = 0/0.25/0.5/0.75/1',
   /a\s*\+=\s*0\.25/.test(drawBody));

ok('draws vertical reference line at het_u',
   /Number\.isFinite\s*\(\s*data\.het_u\s*\)/.test(drawBody) &&
   /het_u/.test(drawBody));

ok('draws vertical reference lines at hom1_u + hom2_u',
   /hom1_u/.test(drawBody) && /hom2_u/.test(drawBody));

ok('draws horizontal reference at agreement = 0.5 (random-guess baseline)',
   /yOf\s*\(\s*0\.5\s*\)/.test(drawBody));

ok('draws points in HOMO_1, HOMO_2, HET order so HET is on top',
   /'HOMO_1',\s*'HOMO_2',\s*'HET'/.test(drawBody));

ok('renders legend with three groups',
   /legendItems/.test(drawBody) &&
   /'HOMO_1'/.test(drawBody) &&
   /'HET'/.test(drawBody) &&
   /'HOMO_2'/.test(drawBody));

ok('returns hits[] for future tooltip wiring',
   /hits:\s*hits/.test(drawBody));

// ============================================================================
// 4. Source-pattern: openVShapePlot popup
// ============================================================================
console.log('\n=== 4. openVShapePlot popup ===');

const openBody = pullFunction(html, 'openVShapePlot');
ok('openVShapePlot body extracted', !!openBody);

ok('lazy-creates #vShapePlotModal',
   /id\s*=\s*['"]vShapePlotModal['"]/.test(openBody) &&
   /createElement\s*\(\s*'div'\s*\)/.test(openBody));

ok('binds Escape key to close',
   /key\s*===\s*'Escape'/.test(openBody));

ok('binds click-outside (modal === ev.target) to close',
   /ev\.target\s*===\s*modal/.test(openBody));

ok('uses _ensureDosageHmState for chunk + last_sq lookup',
   /_ensureDosageHmState/.test(openBody));

ok('checks last_sq.candidate_id matches candidate.id (no cross-candidate stale data)',
   /last_sq\.candidate_id\s*===\s*cand\.id/.test(openBody));

ok('handles missing chunk gracefully (renders helpful message)',
   /load dosage chunk first/.test(openBody));

ok('handles missing stripe-quality gracefully',
   /run stripe-quality compute first/.test(openBody) ||
   /stripe.quality/.test(openBody));

ok('handles no-candidate gracefully',
   /no candidate focused/.test(openBody));

ok('status line reports per-group n + mean_agreement',
   /mean_agreement/.test(openBody) &&
   /HOMO_1/.test(openBody) &&
   /HET/.test(openBody) &&
   /HOMO_2/.test(openBody));

// ============================================================================
// 5. Source-pattern: G-panel karyotype tab button + wiring
// ============================================================================
console.log('\n=== 5. G-panel karyotype tab integration ===');

ok('karyotype tab footer contains gpKaryoVShapeBtn',
   /id="gpKaryoVShapeBtn"/.test(html));

ok('button click handler wired to openVShapePlot',
   /document\.getElementById\s*\(\s*'gpKaryoVShapeBtn'\s*\)/.test(html) &&
   /openVShapePlot\s*\(\s*cand\s*\)/.test(html));

// ============================================================================
// 6. Sandboxed: _buildVShapeData behaviour
// ============================================================================
console.log('\n=== 6. _buildVShapeData sandboxed behaviour ===');

function makeSandbox(opts) {
  opts = opts || {};
  const ctx = {
    state: opts.state || {},
    Float32Array, Float64Array, Int8Array, Int32Array, Uint8Array,
    Map, Set, Number, console, Object, Array, JSON, String, Math,
  };
  ctx.window = opts.window || {};
  ctx.window.state = ctx.state;
  vm.createContext(ctx);
  // Stub _getOrComputeUVRotation. Test injects per-call.
  vm.runInContext('var _getOrComputeUVRotation = null;', ctx);
  // Inject the function under test
  const fn = pullFunction(html, '_buildVShapeData');
  vm.runInContext(fn, ctx);
  return ctx;
}

// 6a. Missing ref_l2 → MISSING_REF_L2
{
  const sb = makeSandbox();
  const out = sb._buildVShapeData({ id: 'X' }, {}, [{}]);
  ok('candidate without ref_l2 → MISSING_REF_L2',
     out.ok === false && out.reason === 'MISSING_REF_L2' && Array.isArray(out.points) && out.points.length === 0);
}

// 6b. UV rotation fail → UV_FAIL
{
  const sb = makeSandbox();
  vm.runInContext('_getOrComputeUVRotation = () => ({ ok: false, reason: "NO_ENV" });', sb);
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, { samples: ['s1'] }, [{}]);
  ok('UV rotation fails → UV_FAIL',
     out.ok === false && out.reason === 'UV_FAIL');
}

// 6c. Missing chunk → MISSING_CHUNK
{
  const sb = makeSandbox();
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, us: new Float64Array([0,1,2]), het_u: 1, hom1_u: 0, hom2_u: 2 });',
    sb
  );
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, null, [{}]);
  ok('missing chunk → MISSING_CHUNK',
     out.ok === false && out.reason === 'MISSING_CHUNK');
}

// 6d. Missing sqRows → MISSING_SQ
{
  const sb = makeSandbox();
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, us: new Float64Array([0,1,2]), het_u: 1, hom1_u: 0, hom2_u: 2 });',
    sb
  );
  const out1 = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, { samples: ['s1'] }, null);
  ok('missing sqRows → MISSING_SQ',
     out1.ok === false && out1.reason === 'MISSING_SQ');
  const out2 = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, { samples: ['s1'] }, []);
  ok('empty sqRows → MISSING_SQ',
     out2.ok === false && out2.reason === 'MISSING_SQ');
}

// 6e. Happy path: full V-shape data
{
  const sb = makeSandbox({
    state: {
      data: {
        samples: [
          { id: 'CGA001' }, { id: 'CGA002' }, { id: 'CGA003' },
          { id: 'CGA004' }, { id: 'CGA005' }, { id: 'CGA006' },
        ],
      },
    },
  });
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, ' +
    '  us: new Float64Array([-2.0, -1.8, 0.1, -0.05, 1.9, 2.1]), ' +
    '  het_u: 0.05, hom1_u: -1.9, hom2_u: 2.0 ' +
    '});',
    sb
  );
  const chunk = { samples: ['CGA001', 'CGA002', 'CGA003', 'CGA004', 'CGA005', 'CGA006'] };
  const sqRows = [
    { sample: 'CGA001', coarse_group: 'HOMO_1', agreement_fraction: 0.95, coherence_class: 'coherent', stripe_quality: 'core' },
    { sample: 'CGA002', coarse_group: 'HOMO_1', agreement_fraction: 0.92, coherence_class: 'coherent', stripe_quality: 'core' },
    { sample: 'CGA003', coarse_group: 'HET',    agreement_fraction: 0.55, coherence_class: 'coherent', stripe_quality: 'core' },
    { sample: 'CGA004', coarse_group: 'HET',    agreement_fraction: 0.62, coherence_class: 'coherent', stripe_quality: 'core' },
    { sample: 'CGA005', coarse_group: 'HOMO_2', agreement_fraction: 0.88, coherence_class: 'coherent', stripe_quality: 'core' },
    { sample: 'CGA006', coarse_group: 'HOMO_2', agreement_fraction: 0.91, coherence_class: 'coherent', stripe_quality: 'core' },
  ];
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, chunk, sqRows);
  ok('happy path → ok=true', out.ok === true);
  ok('happy path → 6 points',
     out.points && out.points.length === 6);
  ok('happy path → het_u/hom1_u/hom2_u carried through',
     Math.abs(out.het_u - 0.05) < 1e-9 &&
     Math.abs(out.hom1_u - -1.9) < 1e-9 &&
     Math.abs(out.hom2_u - 2.0) < 1e-9);
  ok('happy path → first point has cohort u + agreement + group',
     out.points[0].u === -2.0 &&
     out.points[0].agreement === 0.95 &&
     out.points[0].group === 'HOMO_1' &&
     out.points[0].cga === 'CGA001');
  // V-shape signal: HET points should have agreement around 0.5 < 0.7,
  // homo points should have agreement > 0.85.
  const homo1Agreements = out.points.filter(p => p.group === 'HOMO_1').map(p => p.agreement);
  const hetAgreements   = out.points.filter(p => p.group === 'HET').map(p => p.agreement);
  const homo2Agreements = out.points.filter(p => p.group === 'HOMO_2').map(p => p.agreement);
  const meanH1  = homo1Agreements.reduce((a, b) => a + b, 0) / homo1Agreements.length;
  const meanHet = hetAgreements.reduce((a, b) => a + b, 0) / hetAgreements.length;
  const meanH2  = homo2Agreements.reduce((a, b) => a + b, 0) / homo2Agreements.length;
  ok('V-shape signal: mean(HOMO_1) > mean(HET)',
     meanH1 > meanHet);
  ok('V-shape signal: mean(HOMO_2) > mean(HET)',
     meanH2 > meanHet);
}

// 6f. Sample id not in cohort → skipped (not crashed)
{
  const sb = makeSandbox({
    state: { data: { samples: [{ id: 'A' }, { id: 'B' }] } },
  });
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, us: new Float64Array([0, 1]), het_u: 0.5, hom1_u: 0, hom2_u: 1 });',
    sb
  );
  const chunk = { samples: ['A', 'GHOST', 'B'] };
  const sqRows = [
    { sample: 'A',     coarse_group: 'HOMO_1', agreement_fraction: 0.9 },
    { sample: 'GHOST', coarse_group: 'HET',    agreement_fraction: 0.5 },
    { sample: 'B',     coarse_group: 'HOMO_2', agreement_fraction: 0.85 },
  ];
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, chunk, sqRows);
  ok('unknown sample skipped silently → 2 points (not 3, not crash)',
     out.ok === true && out.points.length === 2);
}

// 6g. Non-finite u → skipped
{
  const sb = makeSandbox({
    state: { data: { samples: [{ id: 'A' }, { id: 'B' }] } },
  });
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, us: new Float64Array([NaN, 1]), het_u: 0.5, hom1_u: 0, hom2_u: 1 });',
    sb
  );
  const chunk = { samples: ['A', 'B'] };
  const sqRows = [
    { sample: 'A', coarse_group: 'HOMO_1', agreement_fraction: 0.9 },
    { sample: 'B', coarse_group: 'HOMO_2', agreement_fraction: 0.85 },
  ];
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, chunk, sqRows);
  ok('NaN u → that sample is dropped (1 point survives)',
     out.ok === true && out.points.length === 1 && out.points[0].cga === 'B');
}

// 6h. Non-finite agreement → skipped
{
  const sb = makeSandbox({
    state: { data: { samples: [{ id: 'A' }, { id: 'B' }] } },
  });
  vm.runInContext(
    '_getOrComputeUVRotation = () => ({ ok: true, us: new Float64Array([0, 1]), het_u: 0.5, hom1_u: 0, hom2_u: 1 });',
    sb
  );
  const chunk = { samples: ['A', 'B'] };
  const sqRows = [
    { sample: 'A', coarse_group: 'HOMO_1', agreement_fraction: NaN },
    { sample: 'B', coarse_group: 'HOMO_2', agreement_fraction: 0.85 },
  ];
  const out = sb._buildVShapeData({ id: 'X', ref_l2: 5 }, chunk, sqRows);
  ok('NaN agreement → that sample dropped',
     out.ok === true && out.points.length === 1 && out.points[0].cga === 'B');
}

// ============================================================================
// 7. Sandboxed: _vShapeColor
// ============================================================================
console.log('\n=== 7. _vShapeColor ===');

{
  const ctx = { Object, Array, console };
  vm.createContext(ctx);
  const fn = pullFunction(html, '_vShapeColor');
  vm.runInContext(fn, ctx);
  ok('HOMO_1 → blue (#3b6fb6)', ctx._vShapeColor('HOMO_1') === '#3b6fb6');
  ok('HET → green (#3cc08a)',   ctx._vShapeColor('HET') === '#3cc08a');
  ok('HOMO_2 → orange (#d97a2c)', ctx._vShapeColor('HOMO_2') === '#d97a2c');
  ok('unknown group → grey #888888', ctx._vShapeColor('whatever') === '#888888');
  ok('null/undefined → grey #888888',
     ctx._vShapeColor(null) === '#888888' && ctx._vShapeColor(undefined) === '#888888');
}

// ============================================================================
// 8. Existing flow preserved
// ============================================================================
console.log('\n=== 8. Existing flow preserved ===');

ok('computeStripeQuality still exists',
   /function\s+computeStripeQuality\s*\(\s*chunk\s*,\s*sampleGroup\s*,\s*samplePC1\s*\)/.test(html));

ok('_getOrComputeUVRotation still exists',
   /function\s+_getOrComputeUVRotation\s*\(\s*l2idx\s*\)/.test(html));

ok('_BR_KARYO_COLORS still defined',
   /const\s+_BR_KARYO_COLORS\s*=\s*\{/.test(html));

ok('inheritance matrix popup (turn 122) still exists',
   /function\s+openInheritanceMatrix\s*\(\s*opts\s*\)/.test(html));

ok('G-panel karyotype tab body still has color/export/calibration buttons',
   /id="gpKaryoColorByBtn"/.test(html) &&
   /id="gpKaryoExportBtn"/.test(html) &&
   /id="gpHlabelCalibrationBtn"/.test(html));

ok('turn 155 _inheritanceCacheKey signature still extended (no regression)',
   /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
