// Atlas/tests/test_catalogue_page17.js
// Smoke test — page17 (stats profile) module imports + key surface checks.
import * as page17 from '../inversion_catalogue/page17.js';

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  ✓ ' + msg); pass++; }
  else      { console.log('  ✗ ' + msg); fail++; }
}

console.log('page17 (stats profile) smoke test');
ok(typeof page17 === 'object',                                  'module imports as object');
ok(page17.__MODULE_ID__ === 'inversion_catalogue/page17',       'correct __MODULE_ID__');
ok(typeof page17.renderStatsProfilePage === 'function',         'renderStatsProfilePage exported');

// Constants
ok(page17.SP_ALPHA === 0.05,                                    'SP_ALPHA = 0.05');
ok(typeof page17.SP_EFFECT_COLORS === 'object',                 'SP_EFFECT_COLORS exported');
ok('enriched'   in page17.SP_EFFECT_COLORS,                     'SP_EFFECT_COLORS.enriched present');
ok('depleted'   in page17.SP_EFFECT_COLORS,                     'SP_EFFECT_COLORS.depleted present');
ok(Array.isArray(page17.SP_DEFAULT_ROWS),                       'SP_DEFAULT_ROWS is an array');
ok(page17.SP_DEFAULT_ROWS.length > 0,                           'SP_DEFAULT_ROWS non-empty');

// Helpers
ok(typeof page17._spEnsureState === 'function',                 '_spEnsureState exported');
ok(typeof page17._spDeriveAllRows === 'function',               '_spDeriveAllRows exported');
ok(typeof page17._spIsValidProfileJson === 'function',          '_spIsValidProfileJson exported');
ok(typeof page17._spParseTsv === 'function',                    '_spParseTsv exported');

// Validation helpers should be pure — exercise them with simple inputs.
// (Legacy returns falsy, not strict false, on invalid input.)
ok(!page17._spIsValidProfileJson(null),                         '_spIsValidProfileJson(null) is falsy');
ok(!page17._spIsValidProfileJson({}),                           '_spIsValidProfileJson({}) is falsy');
ok(!page17._spIsValidProfileJson({ metadata: {} }),             '_spIsValidProfileJson(no summary_rows) is falsy');
ok(!!page17._spIsValidProfileJson({ metadata: {}, summary_rows: [] }), '_spIsValidProfileJson(valid shape) is truthy');

console.log('pass: ' + pass + '   fail: ' + fail);
process.exit(fail > 0 ? 1 : 0);
