// v3.98 tests: NULL_CONCORD K=6 fix + v3.80 UI text rewrite
//
// Two independent items shipped this increment:
//   (A) NULL_CONCORD K=6 hardcoded 0.30 → 0.24 (Monte Carlo at N=226,
//       5000 trials each, balanced random labels, Hungarian alignment).
//       Also updated the fallback default from 0.30 → 0.24.
//   (B) v3.80 UI text rewrite — enforce the v3.80 model in user-facing
//       text: intervals are PROMOTED to candidates; fish are ASSIGNED
//       to regimes; "draft" is an interval, not a candidate.
//
// TESTS (8)
//   1. SCRUBBER_VERSION is v3.98+
//   2. NULL_CONCORD has K=6 entry with value 0.24; fallback is 0.24
//   3. NULL_CONCORD comment cites Monte Carlo / N=226 / 5000 trials
//      (so future-me sees where the value came from)
//   4. Page 2 empty-state mentions "promote", "candidate list",
//      "provisional", "page 4 catalogue" (corrected tab number)
//   5. Page 5 (karyotype) empty-state mentions "promoted candidate",
//      "promote" via page 1, and "page 4 catalogue" (corrected tab)
//   6. Catalogue button title says "Promote ... provisional candidates"
//   7. Help-page Vocabulary section exists with rows for L1/L2/L3,
//      Interval, Candidate, Regime, Boundary zone
//   8. Candidate-mode tooltip uses "draft interval" (not "draft candidate")

const fs = require('fs');
const html = fs.readFileSync('/home/claude/v3/pca_scrubber_v3.html', 'utf8');

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  PASS ' + msg); pass++; }
  else      { console.log('  FAIL ' + msg); fail++; }
}

// =============================================================================
// TEST 1: SCRUBBER_VERSION v3.98+
// =============================================================================
console.log('--- TEST 1: SCRUBBER_VERSION v3.98+ ---');
const verMatch = html.match(/const SCRUBBER_VERSION = 'v3\.(\d+)/);
const verNum = verMatch ? parseInt(verMatch[1], 10) : 0;
ok(verNum >= 98, `version is v3.98+ (got v3.${verNum})`);

// =============================================================================
// TEST 2: NULL_CONCORD K=6 entry + fallback both = 0.24
// =============================================================================
console.log('\n--- TEST 2: NULL_CONCORD K=6 fix ---');
const ncTable = html.match(/const NULL_CONCORD = \{([^}]+)\}/);
ok(ncTable !== null, 'NULL_CONCORD table found');
if (ncTable) {
  const tableStr = ncTable[1];
  ok(/2:\s*0\.55/.test(tableStr), 'K=2: 0.55');
  ok(/3:\s*0\.39/.test(tableStr), 'K=3: 0.39');
  ok(/4:\s*0\.30/.test(tableStr), 'K=4: 0.30');
  ok(/5:\s*0\.25/.test(tableStr), 'K=5: 0.25');
  ok(/6:\s*0\.24/.test(tableStr), 'K=6: 0.24 (was missing pre-v3.98, fell through to 0.30)');
}
const fallbackMatch = html.match(/const nullConcord = NULL_CONCORD\[K\] \|\| ([\d.]+)/);
ok(fallbackMatch !== null && fallbackMatch[1] === '0.24',
   `fallback = 0.24 (got ${fallbackMatch ? fallbackMatch[1] : 'NOT FOUND'})`);

// =============================================================================
// TEST 3: NULL_CONCORD comment cites empirical Monte Carlo
// =============================================================================
console.log('\n--- TEST 3: NULL_CONCORD provenance comment ---');
// Locate the comment block immediately above the NULL_CONCORD declaration
const ncIdx = html.indexOf('const NULL_CONCORD = {');
const commentBlock = ncIdx > 0 ? html.slice(Math.max(0, ncIdx - 1500), ncIdx) : '';
ok(/Monte Carlo/i.test(commentBlock),
   'comment references "Monte Carlo"');
ok(/N=?226|N\s*=\s*226/.test(commentBlock),
   'comment references N=226 (LG28 cohort)');
ok(/5000\s*trials/.test(commentBlock),
   'comment references 5000 trials');
ok(/N-dependent|N\s*dependent/i.test(commentBlock),
   'comment notes the value is N-dependent');

// =============================================================================
// TEST 4: Page 2 empty-state vocabulary
// =============================================================================
console.log('\n--- TEST 4: Page 2 empty-state v3.80 vocabulary ---');
const empty2Match = html.match(/<div id="candidateEmpty"[^>]*>([\s\S]+?)<\/div>\s*<\/div>\s*<div id="page3"/);
ok(empty2Match !== null, 'page 2 empty-state block found');
if (empty2Match) {
  const txt = empty2Match[1];
  ok(/Promote an interval/i.test(txt) || /promote.*candidate list/i.test(txt),
     'mentions promoting an interval to the candidate list');
  ok(/page 4 catalogue/i.test(txt),
     'uses corrected tab number "page 4 catalogue"');
  ok(/page 1 diagnostic/i.test(txt),
     'mentions "page 1 diagnostic"');
  ok(/provisional/i.test(txt),
     'mentions "provisional"');
  ok(/mark confirmed/i.test(txt),
     'mentions "mark confirmed"');
  ok(!/page 3 catalogue/i.test(txt),
     'no longer says "page 3 catalogue" (was wrong tab number after v3.97)');
}

// =============================================================================
// TEST 5: Page 5 karyotype empty-state
// =============================================================================
console.log('\n--- TEST 5: Page 5 karyotype empty-state ---');
const empty5Match = html.match(/<div id="candKaryoEmpty"[\s\S]+?<\/div>\s*<\/div>\s*<div id="candKaryoContent"/);
ok(empty5Match !== null, 'karyotype empty-state block found');
if (empty5Match) {
  const txt = empty5Match[0];
  ok(/promote/i.test(txt) && /candidate/i.test(txt),
     'uses "promote" + "candidate" vocabulary');
  ok(/page 4 catalogue/i.test(txt),
     'uses corrected tab number "page 4 catalogue"');
  ok(!/page 3 catalogue/i.test(txt),
     'no longer says wrong tab "page 3 catalogue"');
}

// =============================================================================
// TEST 6: Catalogue button title
// =============================================================================
console.log('\n--- TEST 6: Catalogue button title ---');
const catBtnMatch = html.match(/id="catViewAsCandidate"[^>]+title="([^"]+)"/);
ok(catBtnMatch !== null, 'catalogue button found');
if (catBtnMatch) {
  const title = catBtnMatch[1];
  ok(/Promote/i.test(title),
     `title says "Promote..." (got: ${title.slice(0, 60)}...)`);
  ok(/provisional/i.test(title),
     'title says "provisional"');
}

// =============================================================================
// TEST 7: Help-page Vocabulary section
// =============================================================================
console.log('\n--- TEST 7: Help-page Vocabulary section ---');
const vocabMatch = html.match(/<h4>Vocabulary<\/h4>([\s\S]+?)<h4>/);
ok(vocabMatch !== null, 'Vocabulary section exists in help page');
if (vocabMatch) {
  const txt = vocabMatch[1];
  ok(/L1\s*\/\s*L2\s*\/\s*L3/.test(txt), 'has L1/L2/L3 row');
  ok(/<b>Interval<\/b>/.test(txt), 'has Interval row');
  ok(/<b>Candidate<\/b>/.test(txt), 'has Candidate row');
  ok(/<b>Regime<\/b>/.test(txt), 'has Regime row');
  ok(/<b>Boundary zone<\/b>/.test(txt), 'has Boundary zone row');
  // Vocabulary primer enforces the v3.80 lifecycle
  ok(/promoted/i.test(txt) && /candidate list/i.test(txt),
     'mentions promotion to the candidate list');
  ok(/g0\s*\/\s*g1\s*\/\s*g2/.test(txt) || /homo1.*het.*homo2/i.test(txt),
     'regime row defines g0/g1/g2 = homo1/het/homo2');
  ok(/fish_regime_calls/i.test(txt) || /fish.regime/i.test(txt),
     'regime row references fish_regime_calls.tsv');
}

// =============================================================================
// TEST 8: Candidate-mode tooltip
// =============================================================================
console.log('\n--- TEST 8: Candidate-mode tooltip ---');
const candModeMatch = html.match(/id="candidateModeBtn"[^>]+title="([^"]+)"/);
ok(candModeMatch !== null, 'candidate-mode button found');
if (candModeMatch) {
  const title = candModeMatch[1];
  // The fix: "draft candidate" → "draft interval"
  ok(/draft interval/i.test(title),
     `says "draft interval" (got: ${title.slice(0, 80)}...)`);
  ok(!/draft candidate/i.test(title),
     'no longer says "draft candidate"');
  ok(/promote/i.test(title),
     'mentions promote');
}

// =============================================================================
// SUMMARY
// =============================================================================
console.log(`\n--- v3.98 tests: ${pass}/${pass + fail} ---`);
if (fail > 0) process.exit(1);
