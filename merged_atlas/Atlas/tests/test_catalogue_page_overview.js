// tests/test_catalogue_page_overview.js
// Smoke tests for inversion_catalogue/page_overview.js.
// Legacy page_overview is an empty stub (legacy line 9322), so the
// surface here is intentionally tiny — we just confirm the module
// imports and wires without throwing.

import * as pageOv from '../inversion_catalogue/page_overview.js';

let pass = 0, fail = 0;
function check(name, ok) {
  if (ok) { console.log('  ✓', name); pass++; }
  else    { console.log('  ✗', name); fail++; }
}

check('exports wirePageOverview', typeof pageOv.wirePageOverview === 'function');
check('exports default', typeof pageOv.default === 'function');

const handle = pageOv.wirePageOverview({ data: {} });
check('wirePageOverview returns object', handle && typeof handle === 'object');
check('handle has renderPageOverview', typeof handle.renderPageOverview === 'function');

let threw = false;
try { handle.renderPageOverview(); } catch (e) { threw = true; }
check('renderPageOverview is a no-op (does not throw)', !threw);

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
