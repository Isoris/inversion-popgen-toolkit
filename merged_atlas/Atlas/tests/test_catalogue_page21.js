// Atlas/tests/test_catalogue_page21.js
// Smoke test — page21 (annotation cockpit) module imports + key surface +
// pure helper sanity checks (cursor lookup, chrom-extent computation,
// band-color palette).
import * as page21 from '../inversion_catalogue/page21.js';

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('  ✓ ' + msg); pass++; }
  else      { console.log('  ✗ ' + msg); fail++; }
}

console.log('page21 (annotation cockpit) smoke test');
ok(typeof page21 === 'object',                                  'module imports as object');
ok(page21.__MODULE_ID__ === 'inversion_catalogue/page21',       'correct __MODULE_ID__');
ok(typeof page21.refreshAnnotationCockpit === 'function',       'refreshAnnotationCockpit exported');

// Internal helpers
ok(typeof page21._annoCockpitCandidateAtCursor === 'function',  '_annoCockpitCandidateAtCursor exported');
ok(typeof page21._annoCockpitChromExtent === 'function',        '_annoCockpitChromExtent exported');
ok(typeof page21._ackBandColor === 'function',                  '_ackBandColor exported');
ok(typeof page21._ackEnsureState === 'function',                '_ackEnsureState exported');

// Constants
ok(page21._ACK_PAD && typeof page21._ACK_PAD === 'object',      '_ACK_PAD exported');
ok(Array.isArray(page21._ACK_BAND_PALETTE),                     '_ACK_BAND_PALETTE is array');
ok(page21._ACK_BAND_PALETTE.length === 6,                       '_ACK_BAND_PALETTE has 6 colors');

// _ackBandColor is pure — null/negative → grey, valid index → palette entry
ok(page21._ackBandColor(null) === '#666',                       '_ackBandColor(null) → grey');
ok(page21._ackBandColor(-1) === '#666',                         '_ackBandColor(-1) → grey');
ok(page21._ackBandColor(0) === page21._ACK_BAND_PALETTE[0],     '_ackBandColor(0) → palette[0]');
ok(page21._ackBandColor(7) === page21._ACK_BAND_PALETTE[1],     '_ackBandColor(7) wraps modulo palette length');

// _annoCockpitCandidateAtCursor — pure: returns null on null/empty input
ok(page21._annoCockpitCandidateAtCursor(null, []) === null,     'cursor lookup with null mb → null');
ok(page21._annoCockpitCandidateAtCursor(5.0, null) === null,    'cursor lookup with null items → null');

// Build a tiny synthetic candidate set and look up by mb
const items = [
  { start_bp: 1_000_000, end_bp: 2_000_000, candidate_id: 'A' },
  { start_bp: 5_000_000, end_bp: 6_000_000, candidate_id: 'B' },
];
const hit = page21._annoCockpitCandidateAtCursor(1.5, items);  // 1.5 Mb -> A
ok(hit && hit.candidate_id === 'A',                             'cursor 1.5 Mb hits candidate A');
const miss = page21._annoCockpitCandidateAtCursor(3.0, items); // 3 Mb -> none
ok(miss === null,                                               'cursor 3 Mb (gap) → null');

console.log('pass: ' + pass + '   fail: ' + fail);
process.exit(fail > 0 ? 1 : 0);
