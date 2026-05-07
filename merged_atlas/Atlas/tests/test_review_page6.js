// tests/test_review_page6.js
// Smoke test: page6.js imports without throwing.

import * as page6 from '../inversion_review/page6.js';

let pass = 0;
let fail = 0;

function check(name, cond) {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}`); fail++; }
}

console.log('test_review_page6 — popstats module surface');
check('page6 imports as object', typeof page6 === 'object');
check('showPopstatsPage exported',    typeof page6.showPopstatsPage    === 'function');
check('refreshPopstatsPage exported', typeof page6.refreshPopstatsPage === 'function');

// Calling without window should be a graceful no-op.
let threw = false;
try { page6.showPopstatsPage({ candidate: null }); } catch (_) { threw = true; }
check('showPopstatsPage tolerates missing renderer (no throw)', !threw);

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
