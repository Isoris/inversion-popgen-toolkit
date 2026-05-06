// Tests for v4 turn 89 — Merge isolation + contingency-table separation.
// Coverage:
//   1. Helpers exposed
//   2. assertSameMode: matching tags pass
//   3. assertSameMode: mismatched tags fail with warning
//   4. assertSameMode: legacy untagged candidates pass (default + default)
//   5. assertSameMode: null candidates always pass (no-op)
//   6. buildContingencyForCandidates: same-mode, valid candidates → contingency
//   7. buildContingencyForCandidates: cross-mode → null
//   8. buildContingencyForCandidates: missing locked_labels → null
//   9. mergeIsolationAudit: clean state passes
//  10. mergeIsolationAudit: cross-contamination caught
//  11. Audit catches default-tag in detailed slot
//  12. Audit catches detailed-tag in default slot

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[mi] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[mi] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// Suppress console.warn output during tests
let warnCount = 0;
const realWarn = win.console.warn;
win.console.warn = function(...args) { warnCount++; };

// ---- 1. Helpers exposed ----
for (const name of ['assertSameMode', 'buildContingencyForCandidates',
                    'mergeIsolationAudit']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed on window');
}
ok('all 3 merge-isolation helpers exposed on window');

// ---- 2. Same mode passes ----
const candA_def = { id: 'A', _system: 'default', K: 3, locked_labels: [0, 1, 2] };
const candB_def = { id: 'B', _system: 'default', K: 3, locked_labels: [0, 0, 1] };
warnCount = 0;
if (!win.assertSameMode(candA_def, candB_def)) fail('two default candidates should match');
if (warnCount !== 0) fail('matching default candidates should not warn');
ok('assertSameMode: two default candidates pass, no warning');

const candA_det = { id: 'A', _system: 'detailed', K: 3, locked_labels: [0, 1, 2] };
const candB_det = { id: 'B', _system: 'detailed', K: 3, locked_labels: [0, 0, 1] };
warnCount = 0;
if (!win.assertSameMode(candA_det, candB_det)) fail('two detailed candidates should match');
if (warnCount !== 0) fail('matching detailed candidates should not warn');
ok('assertSameMode: two detailed candidates pass, no warning');

// ---- 3. Mismatched modes fail with warning ----
warnCount = 0;
if (win.assertSameMode(candA_def, candA_det)) fail('cross-mode should fail');
if (warnCount === 0) fail('cross-mode should emit a warning');
ok('assertSameMode: cross-mode (default + detailed) rejected with warning');

warnCount = 0;
if (win.assertSameMode(candA_det, candA_def)) fail('reverse cross-mode should also fail');
if (warnCount === 0) fail('reverse cross-mode should also warn');
ok('assertSameMode: reverse cross-mode also rejected with warning');

// ---- 4. Untagged candidates treated as default ----
const legacyA = { id: 'legacyA', K: 3, locked_labels: [0, 1, 2] };
const legacyB = { id: 'legacyB', K: 3, locked_labels: [1, 0, 2] };
warnCount = 0;
if (!win.assertSameMode(legacyA, legacyB)) fail('two untagged candidates should match (default+default)');
if (warnCount !== 0) fail('untagged + untagged should not warn');
ok('assertSameMode: two untagged candidates treated as default+default, pass');

// Untagged + tagged-default also passes
if (!win.assertSameMode(legacyA, candA_def)) fail('untagged + default-tagged should match');
ok('assertSameMode: untagged + tagged-default pass');

// Untagged + detailed-tagged fails (untagged is default-by-default)
warnCount = 0;
if (win.assertSameMode(legacyA, candA_det)) fail('untagged + detailed should fail');
if (warnCount === 0) fail('untagged vs detailed should warn');
ok('assertSameMode: untagged + tagged-detailed rejected (untagged defaults to "default")');

// ---- 5. Null candidates always pass ----
warnCount = 0;
if (!win.assertSameMode(null, candA_def)) fail('null + default should pass');
if (!win.assertSameMode(candA_det, null)) fail('detailed + null should pass');
if (!win.assertSameMode(null, null)) fail('null + null should pass');
if (warnCount !== 0) fail('null arguments should not warn');
ok('assertSameMode: null arguments always pass without warning (no-op)');

// ---- 6. buildContingencyForCandidates: same-mode → contingency ----
const ct1 = win.buildContingencyForCandidates(candA_def, candB_def);
if (!ct1) fail('same-mode contingency should not be null');
if (!Array.isArray(ct1.M)) fail('contingency M should be array');
if (ct1.KA !== 3) fail('KA should be 3');
if (ct1.KB !== 3) fail('KB should be 3');
if (ct1.n !== 3) fail('n should be 3 (3 samples in both)');
ok('buildContingencyForCandidates: default+default → KA=3, KB=3, n=3 contingency');

// Verify the contingency math
//   labels A = [0, 1, 2], B = [0, 0, 1]
//   M[a][b]:  M[0][0]=1 (sample 0: A=0, B=0)
//             M[1][0]=1 (sample 1: A=1, B=0)
//             M[2][1]=1 (sample 2: A=2, B=1)
if (ct1.M[0][0] !== 1) fail('M[0][0] should be 1');
if (ct1.M[1][0] !== 1) fail('M[1][0] should be 1');
if (ct1.M[2][1] !== 1) fail('M[2][1] should be 1');
ok('contingency cell counts correct: M[0][0]=1, M[1][0]=1, M[2][1]=1');

// ---- 7. Cross-mode contingency → null ----
warnCount = 0;
const ctCross = win.buildContingencyForCandidates(candA_def, candA_det);
if (ctCross !== null) fail('cross-mode contingency should be null');
if (warnCount === 0) fail('cross-mode contingency should warn');
ok('buildContingencyForCandidates: cross-mode → null + warning');

// ---- 8. Missing locked_labels → null ----
const candNoLabels = { id: 'X', _system: 'default', K: 3 };
const ctNoLabels = win.buildContingencyForCandidates(candA_def, candNoLabels);
if (ctNoLabels !== null) fail('missing locked_labels should give null');
ok('buildContingencyForCandidates: missing locked_labels → null');

// Both null
if (win.buildContingencyForCandidates(null, null) !== null) fail('null inputs should give null');
ok('buildContingencyForCandidates: null inputs → null');

// ---- 9. Audit clean state passes ----
// Clean state: clear everything first
win.state.candidate = null;
win.state.candidates = {};
win.state.candidateList = [];
win.clearDetailedState();
const audit1 = win.mergeIsolationAudit();
if (!audit1.ok) fail('clean state audit should pass, got violations: ' + JSON.stringify(audit1.violations));
if (audit1.violations.length !== 0) fail('clean state should have 0 violations');
ok('mergeIsolationAudit: clean state passes (0 violations)');

// ---- 10. Populated clean state passes ----
win.state.candidate = { id: 'def_active', _system: 'default' };
win.state.candidates = {
  'def_a': { id: 'def_a', _system: 'default' },
  'def_b': { id: 'def_b' },   // legacy untagged - should be fine
};
win.state.candidateList = [{ id: 'def_a', _system: 'default' }];
win.state.candidate_detailed = { id: 'det_active', _system: 'detailed' };
win.state.candidates_detailed = {
  'det_a': { id: 'det_a', _system: 'detailed' },
};
win.state.candidateList_detailed = [{ id: 'det_a', _system: 'detailed' }];

const audit2 = win.mergeIsolationAudit();
if (!audit2.ok) {
  fail('populated clean state should pass, got violations: ' + JSON.stringify(audit2.violations));
}
ok('mergeIsolationAudit: populated clean state passes (no cross-contamination)');

// ---- 11. Detailed-tag in default slot caught ----
win.state.candidates['contaminant'] = { id: 'contaminant', _system: 'detailed' };
const audit3 = win.mergeIsolationAudit();
if (audit3.ok) fail('detailed candidate in default slot should be caught');
const v = audit3.violations.find(x => x.id === 'contaminant');
if (!v) fail('violation for contaminant should be reported');
if (v.expected !== 'default') fail('expected should be "default"');
if (v.actual !== 'detailed') fail('actual should be "detailed"');
ok('mergeIsolationAudit: detailed-tagged candidate in candidates[] caught');

// Also catches in candidateList
win.state.candidateList.push({ id: 'list_contam', _system: 'detailed' });
const audit4 = win.mergeIsolationAudit();
const lv = audit4.violations.find(x => x.id === 'list_contam');
if (!lv) fail('violation for list_contam should be reported');
ok('mergeIsolationAudit: detailed-tagged candidate in candidateList caught');

// ---- 12. Default-tag in detailed slot caught ----
// Clean and re-populate
win.state.candidate = null;
win.state.candidates = {};
win.state.candidateList = [];
win.clearDetailedState();
win.state.candidates_detailed['def_in_detailed'] = { id: 'def_in_detailed', _system: 'default' };
const audit5 = win.mergeIsolationAudit();
if (audit5.ok) fail('default candidate in detailed slot should be caught');
const dv = audit5.violations.find(x => x.id === 'def_in_detailed');
if (!dv) fail('violation for def_in_detailed should be reported');
if (dv.expected !== 'detailed') fail('expected should be "detailed"');
if (dv.actual !== 'default') fail('actual should be "default"');
ok('mergeIsolationAudit: default-tagged candidate in candidates_detailed caught');

// Untagged candidate in detailed slot is also caught (treated as default)
win.clearDetailedState();
win.state.candidates_detailed['untagged_in_det'] = { id: 'untagged_in_det' };
const audit6 = win.mergeIsolationAudit();
if (audit6.ok) fail('untagged candidate in detailed slot should be caught');
const uv = audit6.violations.find(x => x.id === 'untagged_in_det');
if (!uv) fail('untagged-in-detailed violation should be reported');
ok('mergeIsolationAudit: untagged candidate in detailed slot caught');

// ---- 13. Audit returns clean after fix ----
win.clearDetailedState();
win.state.candidate = null;
win.state.candidates = {};
win.state.candidateList = [];
const auditFinal = win.mergeIsolationAudit();
if (!auditFinal.ok) fail('after cleanup, audit should pass');
ok('mergeIsolationAudit: clean state after cleanup');

// ---- 14. Round-trip: init from default produces clean state ----
win.state.candidate = { id: 'init_test_a', _system: 'default', K: 3, locked_labels: new Int8Array([0, 1, 2]) };
win.state.candidates = {
  'init_test_a': { id: 'init_test_a', _system: 'default', K: 3, locked_labels: new Int8Array([0, 1, 2]) },
};
win.state.candidateList = [{ id: 'init_test_a', _system: 'default' }];
win.clearDetailedState();
win.initDetailedFromDefault();
const auditAfterInit = win.mergeIsolationAudit();
if (!auditAfterInit.ok) {
  fail('after initDetailedFromDefault, audit should pass — violations: ' + JSON.stringify(auditAfterInit.violations));
}
ok('round-trip: initDetailedFromDefault produces audit-clean state (all detailed candidates tagged correctly)');

// Restore real warn
win.console.warn = realWarn;

console.log('\n[mi] ALL CHECKS PASSED');
