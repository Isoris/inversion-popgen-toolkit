// Atlas/tests/test_catalogue_page3.js
// Smoke test — page3 module imports and exposes its public entry points.
import * as page3 from '../inversion_catalogue/page3.js';

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  ✓ ' + msg); pass++; }
  else      { console.log('  ✗ ' + msg); fail++; }
}

console.log('page3 (catalogue) smoke test');
ok(typeof page3 === 'object',                          'module imports as object');
ok(page3.__MODULE_ID__ === 'inversion_catalogue/page3','correct __MODULE_ID__');
ok(typeof page3.renderCataloguePage === 'function',    'renderCataloguePage exported');
ok(typeof page3.initCataloguePage === 'function',      'initCataloguePage exported');
// renderCataloguePage with no document (Node) should not throw
try { page3.renderCataloguePage(); ok(true, 'renderCataloguePage runs in Node (no DOM) without throwing'); }
catch (e) { ok(false, 'renderCataloguePage threw: ' + e.message); }

console.log('pass: ' + pass + '   fail: ' + fail);
process.exit(fail > 0 ? 1 : 0);
