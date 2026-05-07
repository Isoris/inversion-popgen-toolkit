// Atlas/tests/test_catalogue_page18.js
// Smoke test — page18 (marker readiness panel) module imports + key
// surface + AF scoring purity check.
import * as page18 from '../inversion_catalogue/page18.js';

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  ✓ ' + msg); pass++; }
  else      { console.log('  ✗ ' + msg); fail++; }
}

console.log('page18 (marker readiness panel) smoke test');
ok(typeof page18 === 'object',                              'module imports as object');
ok(page18.__MODULE_ID__ === 'inversion_catalogue/page18',   'correct __MODULE_ID__');
ok(typeof page18.renderMarkerPanelPage === 'function',      'renderMarkerPanelPage exported');

// Key helpers
ok(typeof page18._mpScoreVariantAf === 'function',          '_mpScoreVariantAf exported');
ok(typeof page18._mpAnnotateGelVisibility === 'function',   '_mpAnnotateGelVisibility exported');
ok(typeof page18._mpDeriveAutoPanel === 'function',         '_mpDeriveAutoPanel exported');
ok(typeof page18._mpIsValidPanelJson === 'function',        '_mpIsValidPanelJson exported');
ok(typeof page18._mpIsValidVariantAfsJson === 'function',   '_mpIsValidVariantAfsJson exported');
ok(typeof page18._mpParseTsv === 'function',                '_mpParseTsv exported');

// _mpScoreVariantAf is a pure function — exercise on a Tier-1-clean variant
//   AF_STD <= 0.02, AF_HET in [0.25, 0.75], AF_INV >= 0.80
const t1clean = page18._mpScoreVariantAf({ af_std: 0.01, af_het: 0.50, af_inv: 0.95 });
ok(t1clean && typeof t1clean === 'object',                  'score returns object');
ok(t1clean.tier_from_af === 1,                              'Tier-1-clean variant scores tier=1');
ok(t1clean.private_score > 0.9,                             'private_score > 0.9 for clean variant');
ok(Math.abs(t1clean.dosage_score - 1.0) < 1e-9,             'dosage_score = 1.0 at AF_HET=0.5');

// Validation helpers should accept invalid inputs gracefully.
// (Legacy returns falsy, not strict false, on invalid input.)
ok(!page18._mpIsValidPanelJson(null),                       '_mpIsValidPanelJson(null) is falsy');
ok(!page18._mpIsValidPanelJson({}),                         '_mpIsValidPanelJson({}) is falsy');
ok(!!page18._mpIsValidPanelJson({ markers: [] }),           '_mpIsValidPanelJson({markers:[]}) is truthy');
ok(!page18._mpIsValidVariantAfsJson(null),                  '_mpIsValidVariantAfsJson(null) is falsy');
ok(!page18._mpIsValidVariantAfsJson({}),                    '_mpIsValidVariantAfsJson({}) is falsy');
ok(!!page18._mpIsValidVariantAfsJson({ variants_by_inversion: {} }), '_mpIsValidVariantAfsJson(valid shape) is truthy');

console.log('pass: ' + pass + '   fail: ' + fail);
process.exit(fail > 0 ? 1 : 0);
