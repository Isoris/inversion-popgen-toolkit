// tests/test_review_page_sv_evidence.js
// Smoke test: page_sv_evidence.js imports without throwing.

import * as svEv from '../inversion_review/page_sv_evidence.js';

let pass = 0;
let fail = 0;

function check(name, cond) {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}`); fail++; }
}

console.log('test_review_page_sv_evidence — page_sv_evidence module surface');
check('page_sv_evidence imports as object', typeof svEv === 'object');
check('showSvEvidencePage exported', typeof svEv.showSvEvidencePage === 'function');
check('hideSvEvidencePage exported', typeof svEv.hideSvEvidencePage === 'function');

// Calling without a global window or AtlasSVEvidence should be a no-op.
let threw = false;
try { svEv.showSvEvidencePage({ candidate: null }); } catch (_) { threw = true; }
check('showSvEvidencePage tolerates missing AtlasSVEvidence (no throw)', !threw);

let threw2 = false;
try { svEv.hideSvEvidencePage(); } catch (_) { threw2 = true; }
check('hideSvEvidencePage tolerates missing AtlasSVEvidence (no throw)', !threw2);

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
