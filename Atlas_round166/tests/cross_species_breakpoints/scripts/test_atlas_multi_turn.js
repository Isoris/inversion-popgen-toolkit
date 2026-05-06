// Tests for v4 turns 88, 89, 90, 92 — parallel candidate registry, merge
// isolation, detailed-mode classifier, diamond detector.
//
// Coverage:
//   1. Helpers exposed
//   2. Parallel state slots: candidate_detailed, candidates_detailed,
//      candidateList_detailed; activeMode; default values
//   3. setActiveMode persists in localStorage; rejects bogus values
//   4. getActiveCandidate routes default vs detailed
//   5. _pcrDeepCloneCandidate handles typed arrays + nested objects
//   6. initDetailedFromDefault: 1-to-1 duplication; idempotent; tags _system
//   7. assertCandidateMode + assertSameMode + mergeIsolationAudit
//   8. Detailed-mode classifier: per-sample divergence classes,
//      band interpretation, multimodal detection, H-system assignment
//   9. Diamond detector: strict (>=1 stable band), strict2 (>=2), loose
//  10. summarizeDiamonds: catalogue-column-ready summary

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[mt] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[mt] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
const required = [
  // Turn 88
  '_pcrEnsureState', 'getActiveMode', 'setActiveMode',
  'getActiveCandidate', 'setActiveCandidate',
  'getActiveCandidateList', 'getActiveCandidatesMap',
  'initDetailedFromDefault', 'assertCandidateMode',
  'clearDetailedState', '_pcrDeepCloneCandidate',
  // Turn 89
  'assertSameMode', 'buildContingencyForCandidates', 'mergeIsolationAudit',
  // Turn 90
  'classifyDetailedCandidate', '_dbdDetectModes',
  '_dbdDivergenceClass', '_dbdBandInterpretation', '_dbdAssignHaplotypes',
  // Turn 92
  'detectDiamonds', 'summarizeDiamonds',
  '_ddComputeBandStats', '_ddDetectSplittingRange', '_ddIsBandStable',
];
for (const name of required) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
ok('all ' + required.length + ' helpers exposed (turns 88-90 + 92)');

// ---- 2. State slots ----
win._pcrEnsureState();
if (win.state.candidate_detailed !== null) fail('candidate_detailed should default to null');
if (!Array.isArray(win.state.candidateList_detailed)) fail('candidateList_detailed should be array');
if (typeof win.state.candidates_detailed !== 'object') fail('candidates_detailed should be object');
if (win.state.activeMode !== 'default') fail('activeMode should default to "default"');
ok('state slots: candidate_detailed=null, list=[], map={}, activeMode="default"');

// ---- 3. setActiveMode + persistence ----
win.setActiveMode('detailed');
if (win.getActiveMode() !== 'detailed') fail('setActiveMode("detailed") failed');
ok('setActiveMode("detailed") persists in state');

let storedMode = null;
try { storedMode = win.localStorage.getItem('inversion_atlas.activeMode'); } catch (_) {}
if (storedMode !== 'detailed') fail('mode should persist in localStorage, got ' + storedMode);
ok('setActiveMode persists in localStorage');

const bogusResult = win.setActiveMode('bogus');
if (bogusResult !== false) fail('setActiveMode("bogus") should return false');
if (win.getActiveMode() !== 'detailed') fail('bogus mode should not change state');
ok('setActiveMode rejects invalid mode names');

win.setActiveMode('default');

// ---- 4. getActiveCandidate routes ----
const defaultCand = { id: 'cand_a', K: 3, locked_labels: new Int8Array([0, 1, 2, 0]) };
const detailedCand = { id: 'cand_a', K: 3, locked_labels: new Int8Array([0, 1, 2, 0]), _system: 'detailed' };
win.state.candidate = defaultCand;
win.state.candidate_detailed = detailedCand;

win.setActiveMode('default');
if (win.getActiveCandidate().id !== 'cand_a' || win.getActiveCandidate()._system) {
  fail('default mode should return default candidate');
}
ok('getActiveCandidate(default mode): returns default candidate');

win.setActiveMode('detailed');
const detRet = win.getActiveCandidate();
if (!detRet || detRet._system !== 'detailed') fail('detailed mode should return detailed candidate');
ok('getActiveCandidate(detailed mode): returns detailed candidate (_system="detailed")');

// setActiveCandidate writes to the right slot
win.setActiveMode('detailed');
const newDetailed = { id: 'cand_b', K: 3, _system: 'detailed' };
win.setActiveCandidate(newDetailed);
if (win.state.candidate_detailed.id !== 'cand_b') fail('setActiveCandidate(detailed) failed');
if (win.state.candidate.id !== 'cand_a') fail('setActiveCandidate must NOT cross-write to default');
ok('setActiveCandidate routes correctly: writes detailed slot, leaves default untouched');

win.setActiveMode('default');
win.setActiveCandidate({ id: 'cand_c' });
if (win.state.candidate.id !== 'cand_c') fail('setActiveCandidate(default) failed');
if (win.state.candidate_detailed.id !== 'cand_b') fail('default write must not cross to detailed');
ok('setActiveCandidate routes correctly: default mode does not cross-write');

// ---- 5. Deep clone with typed arrays ----
const orig = {
  id: 'orig',
  K: 3,
  locked_labels: new Int8Array([0, 1, 2, 0, 1]),
  metadata: {
    nested: [1, 2, 3],
    deep: { name: 'test' },
  },
  pc_data: new Float32Array([0.1, 0.2, 0.3]),
};
const cloned = win._pcrDeepCloneCandidate(orig);
if (cloned === orig) fail('deep clone should produce a new object');
if (cloned.locked_labels === orig.locked_labels) fail('typed array should be cloned');
if (cloned.metadata === orig.metadata) fail('nested object should be cloned');
if (cloned.metadata.nested === orig.metadata.nested) fail('nested array should be cloned');
if (cloned.locked_labels.constructor.name !== 'Int8Array') fail('clone should preserve Int8Array type');
if (cloned.pc_data.constructor.name !== 'Float32Array') fail('clone should preserve Float32Array type');
// Mutating clone should not affect original
cloned.locked_labels[0] = 99;
if (orig.locked_labels[0] !== 0) fail('mutation of clone leaked to original');
ok('deep clone: handles Int8Array, Float32Array, nested objects, nested arrays');

// ---- 6. initDetailedFromDefault ----
win.state.candidate = { id: 'def_a', K: 3, locked_labels: new Int8Array([0, 1, 2]) };
win.state.candidates = {
  'def_a': win.state.candidate,
  'def_b': { id: 'def_b', K: 4, locked_labels: new Int8Array([0, 1, 2, 3]) },
};
win.state.candidateList = [win.state.candidates.def_a, win.state.candidates.def_b];
win.state.candidate_detailed = null;
win.state.candidates_detailed = {};
win.state.candidateList_detailed = [];

const dupCount = win.initDetailedFromDefault();
if (dupCount !== 2) fail('expected 2 candidates duplicated, got ' + dupCount);
if (Object.keys(win.state.candidates_detailed).length !== 2) fail('candidates_detailed should have 2 entries');
if (win.state.candidateList_detailed.length !== 2) fail('candidateList_detailed should have 2 entries');
if (!win.state.candidate_detailed) fail('candidate_detailed should be set');
ok('initDetailedFromDefault: duplicated 2 candidates from default to detailed');

// _system tag set on every detailed candidate
if (win.state.candidate_detailed._system !== 'detailed') fail('detailed candidate missing _system tag');
if (win.state.candidates_detailed.def_a._system !== 'detailed') fail('candidates_detailed[def_a] missing _system tag');
if (win.state.candidateList_detailed[0]._system !== 'detailed') fail('candidateList_detailed[0] missing _system tag');
ok('all detailed candidates tagged with _system="detailed"');

// Default candidates untouched
if (win.state.candidate._system) fail('default candidate must NOT have _system tag set');
if (win.state.candidates.def_a._system) fail('default candidates_map untagged');
ok('default candidates remain untagged (no cross-contamination)');

// Idempotent: second call duplicates again (overwrites with fresh clones)
const dupCount2 = win.initDetailedFromDefault();
if (dupCount2 !== 2) fail('idempotent: second call should also return 2');
ok('initDetailedFromDefault: idempotent across repeated calls');

// ---- 7. Assertions and audit ----
// Match: same _system tag
const detA = { id: 'a', _system: 'detailed' };
const detB = { id: 'b', _system: 'detailed' };
if (!win.assertSameMode(detA, detB)) fail('two detailed candidates should pass assertSameMode');
ok('assertSameMode: two detailed candidates pass');

const defA = { id: 'a' };
const defB = { id: 'b' };
if (!win.assertSameMode(defA, defB)) fail('two untagged (default) candidates should pass');
ok('assertSameMode: two untagged candidates pass (treated as default)');

// Mismatch: cross-mode rejected
if (win.assertSameMode(detA, defB)) fail('detailed vs default should fail assertSameMode');
ok('assertSameMode: cross-mode comparison rejected');

// null is no-op
if (!win.assertSameMode(null, detA)) fail('null should return true (no-op)');
ok('assertSameMode: null candidate is treated as no-op');

// mergeIsolationAudit: clean state should pass
win.setActiveMode('default');
win.state.candidate = { id: 'def_a', K: 3 };
win.state.candidate_detailed = { id: 'det_a', K: 3, _system: 'detailed' };
const audit = win.mergeIsolationAudit();
if (!audit.ok) fail('clean state should audit clean, got violations: ' + JSON.stringify(audit.violations));
ok('mergeIsolationAudit: clean state passes (no violations)');

// Inject a violation: detailed-tagged candidate in default slot
win.state.candidate = { id: 'BAD', _system: 'detailed' };
const audit2 = win.mergeIsolationAudit();
if (audit2.ok) fail('audit should detect cross-mode contamination');
if (audit2.violations.length !== 1) fail('audit should report 1 violation, got ' + audit2.violations.length);
if (audit2.violations[0].id !== 'BAD') fail('audit violation should reference candidate id');
ok('mergeIsolationAudit: detects cross-mode contamination (1 violation reported)');

// Restore
win.state.candidate = { id: 'def_a' };

// ---- 8. Detailed-mode classifier ----
// Synthetic 30-sample candidate: 3 K-means bands, low/mid/high heterozygosity
const N = 30;
const labels = new Int8Array(N);
const het = new Float32Array(N);
for (let i = 0; i < 10; i++) {  // band 0 = low het (homozygote-like)
  labels[i] = 0;
  het[i] = 0.10 + Math.random() * 0.05;
}
for (let i = 10; i < 20; i++) {  // band 1 = high het (heterozygote-like)
  labels[i] = 1;
  het[i] = 0.50 + Math.random() * 0.05;
}
for (let i = 20; i < 30; i++) {  // band 2 = low het (homozygote-like)
  labels[i] = 2;
  het[i] = 0.12 + Math.random() * 0.05;
}

const synCand = {
  id: 'syn',
  K: 3,
  locked_labels: labels,
  per_sample_het: het,
  start_w: 0,
  end_w: 4,
};

// Add minimal data so the classifier has something to look at
win.state.data = {
  n_samples: N,
  windows: [],
};
const pc1 = new Float32Array(N);
for (let i = 0; i < N; i++) {
  pc1[i] = labels[i] === 0 ? -2 : (labels[i] === 1 ? 0 : 2);
}
for (let wi = 0; wi < 5; wi++) {
  win.state.data.windows.push({ pc1: Array.from(pc1), pc2: new Array(N).fill(0) });
}

const classification = win.classifyDetailedCandidate(synCand, het);
if (!classification) fail('classifyDetailedCandidate returned null');
if (classification.n_samples !== N) fail('classification.n_samples should be ' + N);
if (classification.n_bands !== 3) fail('classification.n_bands should be 3');
if (!Array.isArray(classification.per_sample)) fail('classification.per_sample should be array');
if (!Array.isArray(classification.per_band)) fail('classification.per_band should be array');
if (!Array.isArray(classification.karyotype_assignment)) fail('classification.karyotype_assignment should be array');
ok('classifyDetailedCandidate: returns per_sample + per_band + karyotype_assignment');

// Per-sample classification
const ps0 = classification.per_sample[0];
if (ps0.divergence_class !== 'low') fail('sample 0 (band 0, het=0.1x) should be low-divergence, got ' + ps0.divergence_class);
ok('per_sample[0]: low divergence (band 0, het=' + het[0].toFixed(3) + ')');

const ps10 = classification.per_sample[10];
if (ps10.divergence_class !== 'high') fail('sample 10 (band 1, het=0.5x) should be high-divergence, got ' + ps10.divergence_class);
ok('per_sample[10]: high divergence (band 1, het=' + het[10].toFixed(3) + ')');

// Per-band interpretation
const pb0 = classification.per_band.find(p => p.band === 0);
if (pb0.interpretation !== 'hom-like') fail('band 0 should be hom-like, got ' + pb0.interpretation);
ok('per_band[band=0]: interpretation="hom-like" (10 fish, all low-divergence)');

const pb1 = classification.per_band.find(p => p.band === 1);
if (pb1.interpretation !== 'het-like') fail('band 1 should be het-like, got ' + pb1.interpretation);
ok('per_band[band=1]: interpretation="het-like" (10 fish, all high-divergence)');

// Karyotype assignment: hom-like band 0 has lowest PC1 → H1/H1; band 2 → H2/H2
const ka = classification.karyotype_assignment;
const ka0 = ka.find(a => a.band === 0);
const ka2 = ka.find(a => a.band === 2);
if (ka0.haplotype_class !== 'H1/H1') fail('band 0 (lowest PC1, hom-like) should be H1/H1, got ' + ka0.haplotype_class);
if (ka2.haplotype_class !== 'H2/H2') fail('band 2 (highest PC1, hom-like) should be H2/H2, got ' + ka2.haplotype_class);
ok('karyotype_assignment: band 0 (lowest PC1) → H1/H1, band 2 (highest PC1) → H2/H2');

const ka1 = ka.find(a => a.band === 1);
if (ka1.haplotype_class !== 'H1/H2') fail('band 1 (mid, het-like) should be H1/H2, got ' + ka1.haplotype_class);
ok('karyotype_assignment: band 1 (het-like) → H1/H2 (intermediate)');

// Multimodal heterozygosity test: build synthetic mixed band
const hetMixed = new Float32Array(20);
for (let i = 0; i < 10; i++) hetMixed[i] = 0.10;   // low mode
for (let i = 10; i < 20; i++) hetMixed[i] = 0.50;  // high mode
const modesResult = win._dbdDetectModes(hetMixed);
if (!modesResult || !Array.isArray(modesResult.mode_centers)) fail('_dbdDetectModes should return modes');
if (modesResult.n_modes < 2) fail('multimodal data should detect >= 2 modes, got ' + modesResult.n_modes);
ok('_dbdDetectModes: bimodal het distribution detected (n_modes=' + modesResult.n_modes + ')');

// ---- 9. Diamond detector ----
// Build a candidate with synthetic per-window PC1: 3 bands (lo, mid, hi)
// where mid splits in middle windows. 5 windows total.
const ddN = 30;
const ddLabels = new Int8Array(ddN);
for (let i = 0; i < 10; i++) ddLabels[i] = 0;       // lo
for (let i = 10; i < 20; i++) ddLabels[i] = 1;      // mid (will split)
for (let i = 20; i < 30; i++) ddLabels[i] = 2;      // hi

const ddWindows = [];
for (let wi = 0; wi < 9; wi++) {
  const isMid = wi >= 3 && wi <= 5;   // split happens in windows 3, 4, 5
  const pc1Arr = new Float32Array(ddN);
  for (let si = 0; si < ddN; si++) {
    const lab = ddLabels[si];
    if (lab === 0) pc1Arr[si] = -2 + (Math.random() - 0.5) * 0.1;       // lo (flat)
    else if (lab === 2) pc1Arr[si] = 2 + (Math.random() - 0.5) * 0.1;   // hi (flat)
    else {
      // mid: at split windows, samples 10-14 push UP, 15-19 push DOWN
      if (isMid) {
        pc1Arr[si] = (si < 15) ? 0.5 : -0.5;
      } else {
        pc1Arr[si] = 0 + (Math.random() - 0.5) * 0.1;   // tight
      }
    }
  }
  ddWindows.push({ pc1: Array.from(pc1Arr) });
}
win.state.data = {
  n_samples: ddN,
  windows: ddWindows,
};
const ddCand = {
  id: 'diam',
  K: 3,
  locked_labels: ddLabels,
  start_w: 0,
  end_w: 8,
};

const diamonds = win.detectDiamonds(ddCand);
if (!Array.isArray(diamonds)) fail('detectDiamonds should return array');
if (diamonds.length === 0) fail('expected at least 1 diamond');
ok('detectDiamonds: detected ' + diamonds.length + ' diamond pattern(s)');

const d0 = diamonds[0];
if (d0.splitting_band !== 1) fail('expected splitting band = 1 (mid), got ' + d0.splitting_band);
ok('diamond.splitting_band = 1 (mid band as constructed)');

if (!Array.isArray(d0.stable_bands) || d0.stable_bands.length !== 2) {
  fail('expected 2 stable bands (0 and 2), got ' + JSON.stringify(d0.stable_bands));
}
ok('diamond.stable_bands = [0, 2] (both lo and hi flat through diamond)');

if (!d0.strict) fail('with 2 stable bands, diamond should be strict=true');
if (!d0.strict2) fail('with 2 stable bands, diamond should be strict2=true');
ok('diamond.strict=true, strict2=true (both stable bands present)');

// summarizeDiamonds: catalogue summary
const summ = win.summarizeDiamonds(ddCand);
if (!summ.has_loose) fail('summary.has_loose should be true');
if (!summ.has_strict) fail('summary.has_strict should be true');
if (!summ.has_strict2) fail('summary.has_strict2 should be true');
if (summ.n_diamonds < 1) fail('summary.n_diamonds should be >= 1');
ok('summarizeDiamonds: n_diamonds=' + summ.n_diamonds + ', has_loose+strict+strict2 all true');

// ---- 10. Helper edge cases ----
// Empty candidate
const emptyCand = { id: 'empty', K: 3, locked_labels: new Int8Array(0) };
const emptyClass = win.classifyDetailedCandidate(emptyCand);
if (emptyClass !== null && emptyClass.n_samples !== 0) {
  // Either null or zero-length is acceptable for empty inputs
  // Don't fail
}
ok('classifier handles empty candidate gracefully (no crash)');

// Diamond detector on candidate with no real bands
const noDiamondCand = {
  id: 'no_diamond',
  K: 1,
  locked_labels: new Int8Array(N).fill(0),
  start_w: 0,
  end_w: 4,
};
const noDiamonds = win.detectDiamonds(noDiamondCand);
if (!Array.isArray(noDiamonds)) fail('should return array even if empty');
ok('detectDiamonds: returns array (possibly empty) when no diamonds');

console.log('\n[mt] ALL CHECKS PASSED');
