// tests/test_review_page7.js
// Smoke test: page7.js imports without throwing.

import * as page7 from '../inversion_review/page7.js';

let pass = 0;
let fail = 0;

function check(name, cond) {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}`); fail++; }
}

console.log('test_review_page7 — ancestry module surface');
check('page7 imports as object', typeof page7 === 'object');
check('showAncestryPage exported',    typeof page7.showAncestryPage    === 'function');
check('refreshAncestryPage exported', typeof page7.refreshAncestryPage === 'function');

// Calling without window should be a graceful no-op.
let threw = false;
try { page7.showAncestryPage({ candidate: null }); } catch (_) { threw = true; }
check('showAncestryPage tolerates missing renderer (no throw)', !threw);

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
