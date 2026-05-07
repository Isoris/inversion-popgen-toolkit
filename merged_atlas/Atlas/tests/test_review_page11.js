// tests/test_review_page11.js
// Smoke test: page11.js imports without throwing.

import * as page11 from '../inversion_review/page11.js';

let pass = 0;
let fail = 0;

function check(name, cond) {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}`); fail++; }
}

console.log('test_review_page11 — page11 module surface');
check('page11 imports as object', typeof page11 === 'object');
check('renderBoundariesPage exported',  typeof page11.renderBoundariesPage  === 'function');
check('_bndKeyHandler exported',        typeof page11._bndKeyHandler        === 'function');
check('_bndAttachHotkeys exported',     typeof page11._bndAttachHotkeys     === 'function');
check('_bndDetachHotkeys exported',     typeof page11._bndDetachHotkeys     === 'function');

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
