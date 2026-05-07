// Atlas/tests/test_catalogue_page9.js
// Smoke test — page9 module imports and exposes its public entry points.
import * as page9 from '../inversion_catalogue/page9.js';

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  ✓ ' + msg); pass++; }
  else      { console.log('  ✗ ' + msg); fail++; }
}

console.log('page9 (confirmed carousel) smoke test');
ok(typeof page9 === 'object',                          'module imports as object');
ok(page9.__MODULE_ID__ === 'inversion_catalogue/page9','correct __MODULE_ID__');
ok(typeof page9.refreshConfirmedCarousel === 'function','refreshConfirmedCarousel exported');
ok(typeof page9.initConfirmedCarousel === 'function',   'initConfirmedCarousel exported');
try { page9.refreshConfirmedCarousel(); ok(true, 'refreshConfirmedCarousel runs in Node without throwing'); }
catch (e) { ok(false, 'refreshConfirmedCarousel threw: ' + e.message); }
try { page9.initConfirmedCarousel(); ok(true, 'initConfirmedCarousel runs in Node without throwing'); }
catch (e) { ok(false, 'initConfirmedCarousel threw: ' + e.message); }

console.log('pass: ' + pass + '   fail: ' + fail);
process.exit(fail > 0 ? 1 : 0);
