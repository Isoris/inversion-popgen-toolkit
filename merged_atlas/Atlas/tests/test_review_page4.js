// tests/test_review_page4.js
// Smoke test: page4.js imports without throwing.

import * as page4 from '../inversion_review/page4.js';

let pass = 0;
let fail = 0;

function check(name, cond) {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}`); fail++; }
}

console.log('test_review_page4 — karyotype/tier module surface');
check('page4 imports as object', typeof page4 === 'object');
check('renderCandidateKaryotype exported',  typeof page4.renderCandidateKaryotype  === 'function');
check('renderCandidateTier exported',       typeof page4.renderCandidateTier       === 'function');
check('_refreshSubviewButtonStyles exported', typeof page4._refreshSubviewButtonStyles === 'function');
check('setKaryoSubview exported',           typeof page4.setKaryoSubview           === 'function');
check('karyoState exported with subview',
      typeof page4.karyoState === 'object' &&
      (page4.karyoState.subview === 'karyotype' || page4.karyoState.subview === 'tier'));

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
