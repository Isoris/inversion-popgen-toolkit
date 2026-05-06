// Tests for v4 turn 88 — Parallel candidate registry (data layer).
// Coverage:
//   1. Helpers exposed
//   2. Default activeMode is 'default' on first run
//   3. setActiveMode persists in localStorage; rejects bogus values
//   4. ensureState idempotency: parallel slots created once
//   5. getActiveCandidate routes default vs detailed correctly
//   6. setActiveCandidate writes to the correct slot, never crosses
//   7. Deep clone preserves Int8Array, nested objects, primitives
//   8. initDetailedFromDefault duplicates 1-to-1 with deep clones
//   9. Detailed clone is independent: mutating clone doesn't affect original
//  10. assertCandidateMode catches cross-mode contamination
//  11. clearDetailedState resets only detailed slots
//  12. Mode switch preserves data in both slots

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[pcr] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[pcr] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['_pcrEnsureState', 'getActiveMode', 'setActiveMode',
                    'getActiveCandidate', 'setActiveCandidate',
                    'getActiveCandidateList', 'getActiveCandidatesMap',
                    'initDetailedFromDefault', 'assertCandidateMode',
                    'clearDetailedState', '_pcrDeepCloneCandidate']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed on window');
}
if (!Array.isArray(win._PCR_VALID_MODES)) fail('_PCR_VALID_MODES not exposed');
ok('all 11 parallel-registry helpers + 1 constant exposed on window');

// ---- 2. Default activeMode on first run ----
// Defensively clear any localStorage carryover from previous tests
try { win.localStorage.removeItem('inversion_atlas.activeMode'); } catch (_) {}
delete win.state.activeMode;
const m0 = win.getActiveMode();
if (m0 !== 'default') fail('default mode should be "default", got "' + m0 + '"');
ok('default activeMode on first run: "default"');

// ---- 3. setActiveMode persistence + validation ----
const setRes = win.setActiveMode('detailed');
if (!setRes) fail('setActiveMode("detailed") should return true');
if (win.getActiveMode() !== 'detailed') fail('mode should be detailed');
let stored = null;
try { stored = win.localStorage.getItem('inversion_atlas.activeMode'); } catch (_) {}
if (stored !== 'detailed') fail('detailed mode should persist in localStorage');
ok('setActiveMode("detailed") persists to localStorage');

// Bogus mode rejected
const bogusRes = win.setActiveMode('foobar');
if (bogusRes !== false) fail('bogus mode should return false');
if (win.getActiveMode() !== 'detailed') fail('bogus mode should not change state');
ok('setActiveMode rejects bogus values');

// Round-trip back
win.setActiveMode('default');

// ---- 4. ensureState idempotency ----
const s1 = win._pcrEnsureState();
const s2 = win._pcrEnsureState();
if (s1 !== s2) fail('ensureState should return same state object');
if (s1.candidate_detailed !== null) fail('candidate_detailed should init to null');
if (!Array.isArray(s1.candidateList_detailed)) fail('candidateList_detailed should be array');
if (typeof s1.candidates_detailed !== 'object') fail('candidates_detailed should be object');
ok('ensureState: idempotent, parallel slots created (candidate_detailed=null, list=[], map={})');

// ---- 5. getActiveCandidate routes default vs detailed ----
// Set a default candidate
win.state.candidate = { id: 'cand_default_1', K: 3, locked_labels: new Int8Array([0, 1, 2]) };
win.setActiveMode('default');
const ad1 = win.getActiveCandidate();
if (!ad1 || ad1.id !== 'cand_default_1') fail('default mode should return cand_default_1');
ok('getActiveCandidate (default mode): returns state.candidate');

win.setActiveMode('detailed');
const ad2 = win.getActiveCandidate();
if (ad2 !== null) fail('detailed mode with no detailed candidate should return null, got ' + JSON.stringify(ad2));
ok('getActiveCandidate (detailed mode, none set): returns null');

// ---- 6. setActiveCandidate writes to correct slot ----
win.setActiveCandidate({ id: 'cand_detailed_1', K: 3, _system: 'detailed' });
if (win.state.candidate_detailed.id !== 'cand_detailed_1') fail('setActiveCandidate(detailed) should set candidate_detailed');
if (win.state.candidate.id !== 'cand_default_1') fail('detailed write should NOT affect default candidate');
ok('setActiveCandidate (detailed mode): writes to candidate_detailed only');

win.setActiveMode('default');
win.setActiveCandidate({ id: 'cand_default_2', K: 4 });
if (win.state.candidate.id !== 'cand_default_2') fail('default write should set state.candidate');
if (win.state.candidate_detailed.id !== 'cand_detailed_1') fail('default write should NOT affect detailed');
ok('setActiveCandidate (default mode): writes to candidate only, detailed unchanged');

// ---- 7. Deep clone preserves types ----
const original = {
  id: 'test',
  K: 3,
  locked_labels: new Int8Array([1, 2, 3]),
  pc_scores: new Float32Array([0.1, 0.2, 0.3]),
  l2_indices: [10, 20, 30],
  notes: { author: 'qa', tags: ['x', 'y'] },
  primitive_string: 'hello',
  primitive_number: 42,
  primitive_bool: true,
  primitive_null: null,
};
const clone = win._pcrDeepCloneCandidate(original);
if (clone === original) fail('clone should be a different reference');
if (!(clone.locked_labels instanceof win.Int8Array || clone.locked_labels instanceof Int8Array)) {
  fail('locked_labels should remain Int8Array, got ' + clone.locked_labels.constructor.name);
}
if (clone.locked_labels === original.locked_labels) fail('Int8Array should be cloned, not shared');
if (clone.locked_labels[0] !== 1) fail('Int8Array values should be preserved');
if (!(clone.pc_scores instanceof win.Float32Array || clone.pc_scores instanceof Float32Array)) {
  fail('pc_scores should remain Float32Array');
}
if (clone.l2_indices === original.l2_indices) fail('l2_indices array should be cloned');
if (clone.l2_indices[0] !== 10) fail('l2_indices values should be preserved');
if (clone.notes === original.notes) fail('nested object should be cloned');
if (clone.notes.author !== 'qa') fail('nested string should be preserved');
if (clone.notes.tags === original.notes.tags) fail('nested array should be cloned');
if (clone.primitive_null !== null) fail('null should pass through');
if (clone.primitive_number !== 42) fail('number should pass through');
if (clone.primitive_bool !== true) fail('bool should pass through');
ok('_pcrDeepCloneCandidate: preserves Int8Array, Float32Array, nested objects, primitives');

// Mutate clone, verify original unchanged
clone.locked_labels[0] = 99;
if (original.locked_labels[0] !== 1) fail('mutating clone Int8Array should not affect original');
clone.notes.author = 'changed';
if (original.notes.author !== 'qa') fail('mutating clone nested object should not affect original');
ok('clone independence: mutations on clone do not affect original');

// ---- 8. initDetailedFromDefault ----
// Clear and set up a known default state
win.state.candidate = { id: 'def_a', K: 3, locked_labels: new Int8Array([0, 0, 1]) };
win.state.candidates = {
  'def_a': { id: 'def_a', K: 3, locked_labels: new Int8Array([0, 0, 1]) },
  'def_b': { id: 'def_b', K: 4, locked_labels: new Int8Array([0, 1, 2, 3]) },
};
win.state.candidateList = [
  { id: 'def_a', K: 3 },
  { id: 'def_b', K: 4 },
];
win.clearDetailedState();
const dupCount = win.initDetailedFromDefault();
if (dupCount !== 2) fail('expected 2 candidates duplicated, got ' + dupCount);
if (!win.state.candidates_detailed['def_a']) fail('def_a should be in detailed map');
if (!win.state.candidates_detailed['def_b']) fail('def_b should be in detailed map');
if (win.state.candidateList_detailed.length !== 2) fail('detailed list should have 2 entries');
if (win.state.candidate_detailed.id !== 'def_a') fail('detailed active candidate should match default');
ok('initDetailedFromDefault: 2 candidates duplicated, list + map + active populated');

// All detailed candidates should have _system='detailed' tag
for (const id of Object.keys(win.state.candidates_detailed)) {
  const c = win.state.candidates_detailed[id];
  if (c._system !== 'detailed') fail('detailed candidate ' + id + ' missing _system tag');
}
if (win.state.candidate_detailed._system !== 'detailed') fail('active detailed candidate missing _system tag');
ok('initDetailedFromDefault: all duplicated candidates tagged with _system="detailed"');

// ---- 9. Detailed clone independence ----
win.state.candidates_detailed['def_a'].locked_labels[0] = 99;
if (win.state.candidates['def_a'].locked_labels[0] !== 0) {
  fail('mutating detailed.locked_labels should not affect default.locked_labels');
}
ok('detailed candidate independence: mutating detailed does not affect default');

// ---- 10. assertCandidateMode ----
win.setActiveMode('default');
const defaultCand = { id: 'x', _system: 'default' };
if (!win.assertCandidateMode(defaultCand)) fail('default cand in default mode should pass');
const detailedCand = { id: 'y', _system: 'detailed' };
if (win.assertCandidateMode(detailedCand)) fail('detailed cand in default mode should fail');
ok('assertCandidateMode (default): default cand passes, detailed cand fails');

win.setActiveMode('detailed');
if (win.assertCandidateMode(defaultCand)) fail('default cand in detailed mode should fail');
if (!win.assertCandidateMode(detailedCand)) fail('detailed cand in detailed mode should pass');
ok('assertCandidateMode (detailed): detailed cand passes, default cand fails');

// Pre-turn-88 candidates with no _system tag should be treated as default
const legacyCand = { id: 'legacy' };
win.setActiveMode('default');
if (!win.assertCandidateMode(legacyCand)) fail('untagged cand should default to "default"');
win.setActiveMode('detailed');
if (win.assertCandidateMode(legacyCand)) fail('untagged cand should fail in detailed mode');
ok('assertCandidateMode: untagged candidates treated as "default" mode');

// Null candidate is always fine
if (!win.assertCandidateMode(null)) fail('null candidate should always pass');
ok('assertCandidateMode: null candidate always passes');

// Override mode argument works
win.setActiveMode('default');
if (win.assertCandidateMode(detailedCand, 'default')) fail('detailed cand should fail in default override');
if (!win.assertCandidateMode(detailedCand, 'detailed')) fail('detailed cand should pass in detailed override');
ok('assertCandidateMode: modeOverride argument respected');

// ---- 11. clearDetailedState resets only detailed slots ----
win.setActiveMode('default');
win.state.candidate = { id: 'preserve_me', K: 3 };
win.state.candidates = { 'preserve_me': { id: 'preserve_me', K: 3 } };
win.state.candidateList = [{ id: 'preserve_me' }];
win.state.candidate_detailed = { id: 'detailed_x', _system: 'detailed' };
win.state.candidates_detailed = { 'detailed_x': { id: 'detailed_x' } };
win.state.candidateList_detailed = [{ id: 'detailed_x' }];

win.clearDetailedState();

if (win.state.candidate_detailed !== null) fail('candidate_detailed should be null after clear');
if (Object.keys(win.state.candidates_detailed).length !== 0) fail('candidates_detailed should be empty');
if (win.state.candidateList_detailed.length !== 0) fail('candidateList_detailed should be empty');
if (!win.state.candidate || win.state.candidate.id !== 'preserve_me') fail('default candidate should be preserved');
if (!win.state.candidates['preserve_me']) fail('default candidates map should be preserved');
if (win.state.candidateList.length !== 1) fail('default candidateList should be preserved');
ok('clearDetailedState: clears only detailed slots, preserves default');

// ---- 12. Mode switch preserves data in both slots ----
win.state.candidate = { id: 'def_active' };
win.state.candidate_detailed = { id: 'detailed_active', _system: 'detailed' };

win.setActiveMode('default');
const a1 = win.getActiveCandidate();
if (a1.id !== 'def_active') fail('default mode should return def_active');

win.setActiveMode('detailed');
const a2 = win.getActiveCandidate();
if (a2.id !== 'detailed_active') fail('detailed mode should return detailed_active');

win.setActiveMode('default');
const a3 = win.getActiveCandidate();
if (a3.id !== 'def_active') fail('switching back should return def_active');
ok('mode switch: both candidate slots preserved across switches');

// ---- 13. getActiveCandidateList + getActiveCandidatesMap routing ----
win.state.candidateList = [{ id: 'a' }, { id: 'b' }];
win.state.candidateList_detailed = [{ id: 'd1' }];
win.setActiveMode('default');
if (win.getActiveCandidateList().length !== 2) fail('default list should be length 2');
win.setActiveMode('detailed');
if (win.getActiveCandidateList().length !== 1) fail('detailed list should be length 1');
ok('getActiveCandidateList: routes correctly between default (n=2) and detailed (n=1)');

win.state.candidates = { 'a': {}, 'b': {} };
win.state.candidates_detailed = { 'd1': {} };
win.setActiveMode('default');
if (Object.keys(win.getActiveCandidatesMap()).length !== 2) fail('default map should have 2 entries');
win.setActiveMode('detailed');
if (Object.keys(win.getActiveCandidatesMap()).length !== 1) fail('detailed map should have 1 entry');
ok('getActiveCandidatesMap: routes correctly between default and detailed');

// ---- 14. initDetailedFromDefault is idempotent ----
win.setActiveMode('default');
win.state.candidate = { id: 'idem_test', K: 3 };
win.state.candidates = { 'idem_test': { id: 'idem_test', K: 3 } };
win.state.candidateList = [{ id: 'idem_test' }];
win.clearDetailedState();
win.initDetailedFromDefault();
const list1Length = win.state.candidateList_detailed.length;
win.initDetailedFromDefault();
const list2Length = win.state.candidateList_detailed.length;
if (list1Length !== list2Length) {
  fail('idempotent init: expected list size unchanged, got ' + list1Length + ' then ' + list2Length);
}
ok('initDetailedFromDefault: idempotent (no growth on repeated calls — overwrites cleanly)');

console.log('\n[pcr] ALL CHECKS PASSED');
