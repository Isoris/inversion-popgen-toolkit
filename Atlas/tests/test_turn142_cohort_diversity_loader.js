// =============================================================================
// turn 142 — cohort_diversity_v1 loader (foundation for breeding-readiness card)
// =============================================================================
// SPEC: specs_new_turn131/SPEC_per_candidate_breeding_readiness_card.md
// (Slice 1, Turn A — data plumbing only; computation lands in Turn B,
//  render in Turn C, catalogue export in Turn D.)
//
// What this turn ships:
//   - state.cohortDiversity   — the loaded layer + byCGA index
//   - _isCohortDiversityJSON  — accepts both wrapped + raw-array shapes
//   - _normalizeCohortDiversityRow — coerces NaN / wrong types to null
//   - _storeCohortDiversity   — validates, normalizes, indexes
//   - _persistCohortDiversity / _restoreCohortDiversity — localStorage
//   - _clearCohortDiversity
//   - _diversityForSampleIdx(si) — chrom-sample → diversity row resolver
//   - _diversityForCGA(cga)      — direct CGA → diversity row resolver
//   - _cohortDiversityCoverageOnCurrentChrom() — diagnostic
//   - File-picker dispatch in loadJSON
//   - Startup restore wired in
//
// The fixture json/cohort_diversity_v1.json (226 real samples derived
// from Diversity_atlas dt_S1) is also part of this turn. We test against
// it to catch any drift from the Diversity-atlas-paste convenience path.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const FIXTURE_PATH = path.resolve(__dirname, '..', 'json', 'cohort_diversity_v1.json');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// Helper: build a fresh sandbox with the loader region loaded.
// ============================================================================
function buildSandbox() {
  const re = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?window\.COHORT_DIVERSITY_LS_KEY\s*=\s*COHORT_DIVERSITY_LS_KEY;\s*\n}/;
  const m = html.match(re);
  if (!m) throw new Error('cohort_diversity loader region not found');
  const _ls = new Map();
  const ctx = vm.createContext({
    state: { data: null, cohortDiversity: null },
    window: {},
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
    console,
  });
  vm.runInContext(m[0], ctx);
  // Expose the localStorage map for tests that need to inspect it
  ctx.__ls = _ls;
  return ctx;
}

// ============================================================================
// 1. Source-level definitions present
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('COHORT_DIVERSITY_TOOL constant',
   /const COHORT_DIVERSITY_TOOL\s*=\s*['"]cohort_diversity_v1['"]/.test(html));
ok('COHORT_DIVERSITY_LS_KEY constant',
   /const COHORT_DIVERSITY_LS_KEY\s*=\s*['"]inversion_atlas\.cohort_diversity['"]/.test(html));
ok('_COHORT_DIVERSITY_NUMERIC_COLS constant',
   /const _COHORT_DIVERSITY_NUMERIC_COLS\s*=\s*\[/.test(html));
ok('numeric col list includes f_roh',
   /['"]f_roh['"]/.test(html.match(/_COHORT_DIVERSITY_NUMERIC_COLS\s*=\s*\[[^\]]*\]/)[0]));
ok('numeric col list includes h, f_hom, callable_bp',
   /['"]h['"]/.test(html) && /['"]f_hom['"]/.test(html) && /['"]callable_bp['"]/.test(html));

ok('_isCohortDiversityJSON defined',         /function _isCohortDiversityJSON\(/.test(html));
ok('_normalizeCohortDiversityRow defined',    /function _normalizeCohortDiversityRow\(/.test(html));
ok('_storeCohortDiversity defined',           /function _storeCohortDiversity\(/.test(html));
ok('_persistCohortDiversity defined',         /function _persistCohortDiversity\(/.test(html));
ok('_restoreCohortDiversity defined',         /function _restoreCohortDiversity\(/.test(html));
ok('_clearCohortDiversity defined',           /function _clearCohortDiversity\(/.test(html));
ok('_diversityForSampleIdx defined',          /function _diversityForSampleIdx\(/.test(html));
ok('_diversityForCGA defined',                /function _diversityForCGA\(/.test(html));
ok('_cohortDiversityCoverageOnCurrentChrom defined',
   /function _cohortDiversityCoverageOnCurrentChrom\(/.test(html));

// ============================================================================
// 2. Window exports
// ============================================================================
console.log('\n=== 2. Window exports ===');

ok('window._isCohortDiversityJSON',
   /window\._isCohortDiversityJSON\s*=\s*_isCohortDiversityJSON/.test(html));
ok('window._storeCohortDiversity',
   /window\._storeCohortDiversity\s*=\s*_storeCohortDiversity/.test(html));
ok('window._diversityForSampleIdx',
   /window\._diversityForSampleIdx\s*=\s*_diversityForSampleIdx/.test(html));
ok('window._diversityForCGA',
   /window\._diversityForCGA\s*=\s*_diversityForCGA/.test(html));
ok('window._cohortDiversityCoverageOnCurrentChrom',
   /window\._cohortDiversityCoverageOnCurrentChrom\s*=\s*_cohortDiversityCoverageOnCurrentChrom/.test(html));

// ============================================================================
// 3. File-picker dispatch wired
// ============================================================================
console.log('\n=== 3. File-picker dispatch wired ===');

ok('classifier returns cohort_diversity kind',
   /return\s+['"]cohort_diversity['"]/.test(html));
ok('classifier checks _isCohortDiversityJSON before chromosome',
   // Quick proxy: the classifier text mentions cohort_diversity ABOVE 'chromosome'
   (() => {
     const cIdx = html.indexOf("return 'cohort_diversity'");
     const chIdx = html.indexOf("return 'chromosome'");
     return cIdx > 0 && chIdx > 0 && cIdx < chIdx;
   })());
ok('loadJSON dispatch calls _storeCohortDiversity',
   /_isCohortDiversityJSON\(data\)\)\s*\{[\s\S]{0,500}_storeCohortDiversity\(data\)/.test(html));
ok('loadJSON dispatch calls _persistCohortDiversity',
   /_storeCohortDiversity\(data\)\)\s*\{[\s\S]{0,500}_persistCohortDiversity/.test(html));
ok('loadJSON cohort_diversity branch returns to next()',
   // Roughly: after the storeCohortDiversity branch, there's a next();return; block
   /_storeCohortDiversity[\s\S]{0,800}next\(\);\s*\n\s*return;/.test(html));

// ============================================================================
// 4. Startup restore wired
// ============================================================================
console.log('\n=== 4. Startup restore wired ===');

ok('_restoreCohortDiversity called at startup',
   /_restoreCohortDiversity === 'function'\) _restoreCohortDiversity\(\)/.test(html));

// ============================================================================
// 5. Behavioral: detection
// ============================================================================
console.log('\n=== 5. Detection behavior ===');

const ctx = buildSandbox();
const W = ctx.window;

ok('null → false',           W._isCohortDiversityJSON(null) === false);
ok('undefined → false',      W._isCohortDiversityJSON(undefined) === false);
ok('"" → false',              W._isCohortDiversityJSON('') === false);
ok('empty array → false',    W._isCohortDiversityJSON([]) === false);
ok('empty object → false',   W._isCohortDiversityJSON({}) === false);
ok('chromosome-shape → false',
   W._isCohortDiversityJSON({ chrom: 'LG28', n_windows: 100, samples: [{ ind: 'I1', cga: 'CGA001' }] }) === false);
ok('cs_breakpoints-shape → false',
   W._isCohortDiversityJSON({ tool: 'cross_species_breakpoints_v1', schema_version: 1, breakpoints: [] }) === false);
ok('wrapped form → true',
   W._isCohortDiversityJSON({
     tool: 'cohort_diversity_v1', schema_version: 1, samples: [{ sample_id: 'CGA001', f_roh: 0.2 }],
   }) === true);
ok('wrapped without schema_version → false',
   W._isCohortDiversityJSON({ tool: 'cohort_diversity_v1', samples: [] }) === false);
ok('wrapped without samples[] → false',
   W._isCohortDiversityJSON({ tool: 'cohort_diversity_v1', schema_version: 1 }) === false);
ok('raw-array with sample + f_roh → true',
   W._isCohortDiversityJSON([{ sample: 'CGA001', f_roh: 0.2 }]) === true);
ok('raw-array with sample_id + f_roh → true',
   W._isCohortDiversityJSON([{ sample_id: 'CGA001', f_roh: 0.2 }]) === true);
ok('raw-array without f_roh → false',
   W._isCohortDiversityJSON([{ sample: 'CGA001', h: 0.005 }]) === false);
ok('raw-array without sample id → false',
   W._isCohortDiversityJSON([{ k8: 'K1', f_roh: 0.2 }]) === false);

// ============================================================================
// 6. Behavioral: row normalization
// ============================================================================
console.log('\n=== 6. Row normalization ===');

const goodRow = {
  sample_id: 'CGA009', k8: 'K1', pruned81: true,
  h: 0.0046, f_hom: -0.032, f_roh: 0.254,
  roh_total_bp: 1.4e8, roh_n: 3190, roh_longest_bp: 3221911, roh_mean_bp: 45959,
  th_in: 0.001, th_out: 0.005, th_ratio: 0.24, callable_bp: 5e8,
};
const norm = W._normalizeCohortDiversityRow(goodRow);
ok('good row normalizes', norm !== null);
ok('  sample_id preserved', norm.sample_id === 'CGA009');
ok('  k8 preserved as string', norm.k8 === 'K1');
ok('  pruned81 preserved as bool', norm.pruned81 === true);
ok('  f_roh preserved', norm.f_roh === 0.254);
ok('  h preserved', norm.h === 0.0046);
ok('  callable_bp preserved', norm.callable_bp === 5e8);

ok('null row → null',                    W._normalizeCohortDiversityRow(null) === null);
ok('row missing sample_id → null',       W._normalizeCohortDiversityRow({ f_roh: 0.2 }) === null);
ok('row with empty sample_id → null',    W._normalizeCohortDiversityRow({ sample_id: '', f_roh: 0.2 }) === null);

const altRow = W._normalizeCohortDiversityRow({ sample: 'CGA042', f_roh: 0.3 }); // dt_S1 shape
ok('alt-key (sample) accepted',          altRow !== null && altRow.sample_id === 'CGA042');

const trimRow = W._normalizeCohortDiversityRow({ sample_id: '  CGA005  ', f_roh: 0.2 });
ok('sample_id trimmed',                  trimRow.sample_id === 'CGA005');

const k8intRow = W._normalizeCohortDiversityRow({ sample_id: 'CGA001', k8: 7, f_roh: 0.2 });
ok('k8 int → string',                    k8intRow.k8 === '7');

const pruned81str = W._normalizeCohortDiversityRow({ sample_id: 'CGA001', pruned81: 'true', f_roh: 0.2 });
ok('pruned81 "true" → true',             pruned81str.pruned81 === true);
const pruned81num = W._normalizeCohortDiversityRow({ sample_id: 'CGA001', pruned81: 1, f_roh: 0.2 });
ok('pruned81 1 → true',                  pruned81num.pruned81 === true);
const pruned81bad = W._normalizeCohortDiversityRow({ sample_id: 'CGA001', pruned81: 'maybe', f_roh: 0.2 });
ok('pruned81 garbage → null',            pruned81bad.pruned81 === null);

const nanRow = W._normalizeCohortDiversityRow({ sample_id: 'CGA001', h: NaN, f_roh: 'not numeric', roh_n: 'NA' });
ok('NaN h → null',                       nanRow.h === null);
ok('non-numeric f_roh → null',           nanRow.f_roh === null);
ok('non-numeric roh_n → null',           nanRow.roh_n === null);

// ============================================================================
// 7. Behavioral: store + index
// ============================================================================
console.log('\n=== 7. Store + index ===');

const ctx2 = buildSandbox();
const W2 = ctx2.window;

const wrapped = {
  tool: 'cohort_diversity_v1', schema_version: 1,
  generated_at: '2026-05-05T00:00:00Z',
  cohort: { n_samples: 3 },
  samples: [
    { sample_id: 'CGA009', k8: 'K1', pruned81: true,  h: 0.0046, f_roh: 0.254, callable_bp: 5e8 },
    { sample_id: 'CGA010', k8: 'K3', pruned81: true,  h: 0.0045, f_roh: 0.270, callable_bp: 5e8 },
    { sample_id: 'CGA021', k8: 'K2', pruned81: false, h: 0.0045, f_roh: 0.323, callable_bp: 5e8 },
  ],
};

ok('store wrapped → true',                  W2._storeCohortDiversity(wrapped) === true);
ok('state.cohortDiversity populated',       ctx2.state.cohortDiversity !== null);
ok('  schema_version preserved',            ctx2.state.cohortDiversity.schema_version === 1);
ok('  tool stamped',                         ctx2.state.cohortDiversity.tool === 'cohort_diversity_v1');
ok('  generated_at preserved',               ctx2.state.cohortDiversity.generated_at === '2026-05-05T00:00:00Z');
ok('  cohort metadata deep-copied',          ctx2.state.cohortDiversity.cohort.n_samples === 3);
ok('  samples[] length matches',             ctx2.state.cohortDiversity.samples.length === 3);
ok('  byCGA is a Map',                       (() => {
  const m = ctx2.state.cohortDiversity.byCGA;
  return m && typeof m.get === 'function' && typeof m.set === 'function' && typeof m.has === 'function';
})());
ok('  byCGA size matches',                   ctx2.state.cohortDiversity.byCGA.size === 3);
ok('  diagnostics.n_loaded',                 ctx2.state.cohortDiversity.diagnostics.n_loaded === 3);
ok('  diagnostics.n_dropped_no_id',          ctx2.state.cohortDiversity.diagnostics.n_dropped_no_id === 0);
ok('  loaded_at present',                    typeof ctx2.state.cohortDiversity.loaded_at === 'string');

// Duplicate id detection
const dupCtx = buildSandbox();
dupCtx.window._storeCohortDiversity([
  { sample: 'CGA001', f_roh: 0.1 },
  { sample: 'CGA001', f_roh: 0.5 },  // duplicate, last writer wins
  { sample: 'CGA002', f_roh: 0.2 },
]);
ok('duplicate id counted',                   dupCtx.state.cohortDiversity.diagnostics.n_duplicate_id === 1);
ok('duplicate id last writer wins',          dupCtx.window._diversityForCGA('CGA001').f_roh === 0.5);

// Drop rows with no sample id
const dropCtx = buildSandbox();
dropCtx.window._storeCohortDiversity({
  tool: 'cohort_diversity_v1', schema_version: 1, samples: [
    { sample_id: 'CGA001', f_roh: 0.1 },
    { f_roh: 0.5 },                   // no id, dropped
    { sample_id: 'CGA002', f_roh: 0.2 },
  ],
});
ok('dropped no-id row counted',              dropCtx.state.cohortDiversity.diagnostics.n_dropped_no_id === 1);
ok('valid rows survive drop',                dropCtx.state.cohortDiversity.samples.length === 2);

// All-bad → false
const allbadCtx = buildSandbox();
ok('all-bad-rows wrapped → false',
   allbadCtx.window._storeCohortDiversity({
     tool: 'cohort_diversity_v1', schema_version: 1, samples: [{}, { f_roh: 0.1 }],
   }) === false);

// Non-detection input
ok('non-detection input rejected',
   W2._storeCohortDiversity({ chrom: 'LG28', samples: [] }) === false);

// ============================================================================
// 8. Behavioral: byCGA lookup is case-insensitive
// ============================================================================
console.log('\n=== 8. byCGA lookup ===');

ok('CGA009 finds row',                       W2._diversityForCGA('CGA009').f_roh === 0.254);
ok('cga009 (lowercase) finds row',           W2._diversityForCGA('cga009').f_roh === 0.254);
ok('Cga009 (mixed) finds row',               W2._diversityForCGA('Cga009').f_roh === 0.254);
ok('CGA999 (missing) → null',                W2._diversityForCGA('CGA999') === null);
ok('null arg → null',                         W2._diversityForCGA(null) === null);
ok('non-string arg → null',                  W2._diversityForCGA(42) === null);
ok('"" → null',                               W2._diversityForCGA('') === null);

// ============================================================================
// 9. Behavioral: _diversityForSampleIdx (chrom-sample → row)
// ============================================================================
console.log('\n=== 9. _diversityForSampleIdx ===');

ctx2.state.data = {
  samples: [
    { ind: 'Ind0', cga: 'CGA009' },
    { ind: 'Ind1', cga: 'CGA010' },
    { ind: 'Ind2', cga: 'CGA999' },     // not in cohort diversity
    { ind: 'Ind3' },                      // no cga at all
  ],
};

ok('si=0 (CGA009) resolves',                W2._diversityForSampleIdx(0).f_roh === 0.254);
ok('si=1 (CGA010) resolves',                W2._diversityForSampleIdx(1).f_roh === 0.270);
ok('si=2 (CGA999, missing) → null',         W2._diversityForSampleIdx(2) === null);
ok('si=3 (no cga) → null',                  W2._diversityForSampleIdx(3) === null);
ok('si=99 (out of range) → null',           W2._diversityForSampleIdx(99) === null);
ok('si=-1 → null',                            W2._diversityForSampleIdx(-1) === null);
ok('si=non-int → null',                      W2._diversityForSampleIdx(1.5) === null);
ok('si=string → null',                       W2._diversityForSampleIdx('0') === null);
ok('si=null → null',                          W2._diversityForSampleIdx(null) === null);

// State without data
const noDataCtx = buildSandbox();
noDataCtx.window._storeCohortDiversity(wrapped);
ok('si lookup with no state.data → null',
   noDataCtx.window._diversityForSampleIdx(0) === null);

// State with data but no cohort
const noCohortCtx = buildSandbox();
noCohortCtx.state.data = { samples: [{ cga: 'CGA009' }] };
ok('si lookup with no cohort → null',
   noCohortCtx.window._diversityForSampleIdx(0) === null);

// ============================================================================
// 10. Coverage helper
// ============================================================================
console.log('\n=== 10. Coverage diagnostic ===');

const cov = W2._cohortDiversityCoverageOnCurrentChrom();
ok('coverage.n_total',                       cov.n_total === 4);
ok('coverage.n_resolved',                    cov.n_resolved === 2);
ok('coverage.n_unresolved',                  cov.n_unresolved === 2);
ok('coverage lists unresolved cgas',         Array.isArray(cov.unresolved_cgas));
ok('coverage caps unresolved at 10', (() => {
  // Build a sandbox with 15 unresolved samples
  const c = buildSandbox();
  c.state.data = {
    samples: Array.from({ length: 15 }, (_, i) => ({ ind: 'I' + i, cga: 'CGA9' + (1000 + i) })),
  };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1, samples: [],
  });
  // store fails because samples is empty — reload with one matching row
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [{ sample_id: 'CGA9999', f_roh: 0.1 }],
  });
  const cv = c.window._cohortDiversityCoverageOnCurrentChrom();
  return cv.n_unresolved === 15 && cv.unresolved_cgas.length === 10;
})());

// Coverage when no cohort loaded
const noCohortCov = noCohortCtx.window._cohortDiversityCoverageOnCurrentChrom();
ok('coverage with no cohort: n_resolved=0',  noCohortCov.n_resolved === 0);
ok('coverage with no cohort: n_unresolved=n_total', noCohortCov.n_unresolved === noCohortCov.n_total);

// Coverage when no chrom loaded
const noChromCov = noDataCtx.state ? null : null;
const cleanCtx = buildSandbox();
const noChromCov2 = cleanCtx.window._cohortDiversityCoverageOnCurrentChrom();
ok('coverage with no chrom: n_total=0',      noChromCov2.n_total === 0);

// ============================================================================
// 11. Persist + restore round-trip
// ============================================================================
console.log('\n=== 11. localStorage round-trip ===');

const psCtx = buildSandbox();
psCtx.window._storeCohortDiversity(wrapped);
psCtx.window._persistCohortDiversity();
ok('persist writes localStorage',
   psCtx.__ls.has('inversion_atlas.cohort_diversity'));
const rawLs = psCtx.__ls.get('inversion_atlas.cohort_diversity');
ok('persisted value is JSON',                (() => { try { JSON.parse(rawLs); return true; } catch (_) { return false; } })());
const persistedParsed = JSON.parse(rawLs);
ok('persisted does NOT include byCGA Map',   !('byCGA' in persistedParsed));
ok('persisted preserves samples',             persistedParsed.samples.length === 3);

// Clear in-memory, restore from localStorage
psCtx.state.cohortDiversity = null;
ok('restore returns true',                   psCtx.window._restoreCohortDiversity() === true);
ok('restored state has 3 samples',           psCtx.state.cohortDiversity.samples.length === 3);
ok('restored state rebuilds byCGA Map',       (() => {
  const m = psCtx.state.cohortDiversity.byCGA;
  return m && typeof m.get === 'function' && typeof m.set === 'function' && typeof m.has === 'function';
})());
ok('restored byCGA finds CGA009',             psCtx.window._diversityForCGA('CGA009').f_roh === 0.254);

// Restore with no key returns false
const emptyCtx = buildSandbox();
ok('restore with empty localStorage → false', emptyCtx.window._restoreCohortDiversity() === false);

// Restore with corrupted key returns false (does not throw)
const corruptCtx = buildSandbox();
corruptCtx.__ls.set('inversion_atlas.cohort_diversity', 'NOT JSON');
ok('restore with corrupt localStorage → false (no throw)',
   corruptCtx.window._restoreCohortDiversity() === false);

// ============================================================================
// 12. Clear
// ============================================================================
console.log('\n=== 12. Clear ===');

const clrCtx = buildSandbox();
clrCtx.window._storeCohortDiversity(wrapped);
clrCtx.window._persistCohortDiversity();
ok('before clear: state populated',          clrCtx.state.cohortDiversity !== null);
ok('before clear: localStorage populated',   clrCtx.__ls.has('inversion_atlas.cohort_diversity'));

clrCtx.window._clearCohortDiversity();
ok('after clear: state.cohortDiversity null', clrCtx.state.cohortDiversity === null);
ok('after clear: localStorage removed',      !clrCtx.__ls.has('inversion_atlas.cohort_diversity'));

// ============================================================================
// 13. Real fixture (json/cohort_diversity_v1.json) end-to-end
// ============================================================================
console.log('\n=== 13. Real fixture ===');

let fixture = null;
try { fixture = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf8')); } catch (_) {}

ok('fixture file exists + parses',           fixture !== null);
if (fixture) {
  ok('fixture is wrapped form',              fixture.tool === 'cohort_diversity_v1');
  ok('fixture has 226 samples',              Array.isArray(fixture.samples) && fixture.samples.length === 226);

  const fxCtx = buildSandbox();
  ok('fixture passes detection',             fxCtx.window._isCohortDiversityJSON(fixture) === true);
  ok('fixture stores cleanly',               fxCtx.window._storeCohortDiversity(fixture) === true);
  ok('all 226 samples loaded',               fxCtx.state.cohortDiversity.samples.length === 226);
  ok('zero drops',                            fxCtx.state.cohortDiversity.diagnostics.n_dropped_no_id === 0);
  ok('zero duplicates',                       fxCtx.state.cohortDiversity.diagnostics.n_duplicate_id === 0);

  // Spot-check a few known rows from the Diversity-atlas data
  const r9 = fxCtx.window._diversityForCGA('CGA009');
  ok('CGA009 f_roh ≈ 0.254',                  r9 && Math.abs(r9.f_roh - 0.25413458) < 1e-6, r9 && r9.f_roh);
  ok('CGA009 k8 = K1',                         r9 && r9.k8 === 'K1');

  const r322 = fxCtx.window._diversityForCGA('CGA322');
  ok('CGA322 has very low f_roh (outlier)',   r322 && r322.f_roh < 0.05, r322 && r322.f_roh);
}

// ============================================================================
// 14. Raw-array form (Diversity-atlas dt_S1 paste convenience)
// ============================================================================
console.log('\n=== 14. Raw-array form (dt_S1 paste) ===');

// Simulate what a user would do: open Diversity_atlas, copy the JSON
// blob inside <script id="dt_S1">, save as a .json file, drop into
// the file picker. Detection + store should both work.
const diversityHtml = fs.readFileSync(
  path.resolve(__dirname, '..', 'Diversity_atlas.html'), 'utf8');
const dtMatch = diversityHtml.match(
  /<script type="application\/json" id="dt_S1">([\s\S]*?)<\/script>/);

if (dtMatch) {
  const dtArr = JSON.parse(dtMatch[1]);
  const dtCtx = buildSandbox();
  ok('dt_S1 paste: 226 rows',                dtArr.length === 226);
  ok('dt_S1 paste: detection succeeds',      dtCtx.window._isCohortDiversityJSON(dtArr) === true);
  ok('dt_S1 paste: store succeeds',          dtCtx.window._storeCohortDiversity(dtArr) === true);
  ok('dt_S1 paste: 226 samples loaded',      dtCtx.state.cohortDiversity.samples.length === 226);
  ok('dt_S1 paste: byCGA built',              dtCtx.state.cohortDiversity.byCGA.size === 226);

  // Sample ID was 'sample' in the source; check it's now 'sample_id'
  const someRow = dtCtx.state.cohortDiversity.samples[0];
  ok('dt_S1 paste: sample → sample_id',       typeof someRow.sample_id === 'string' &&
                                              someRow.sample_id.startsWith('CGA'));
}

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
