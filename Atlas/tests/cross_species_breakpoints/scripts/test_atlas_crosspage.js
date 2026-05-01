// Tests for v4 turn 84 — cross-page cluster coloring.
// Coverage:
//   1. Helpers exposed
//   2. Cluster registry lifecycle
//   3. Active candidate id resolution (multiple state field paths)
//   4. Color lookup for a registered candidate
//   5. Color fallback (grey) when source has no clusters for active candidate
//   6. Concordance computation: contingency, agreement, ARI, Cramér's V
//   7. getSampleColor dispatches to cluster_dosage/theta_pi/ghsl correctly
//   8. Cluster colors stable across signal types (label N → same color)

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[xp] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[xp] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['_xpEnsureClusterRegistry', '_xpRegisterClusters',
                    '_xpLookupClusters', '_xpActiveCandidateId',
                    '_xpSampleColor', '_xpClusterColor',
                    '_xpComputeConcordance']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
ok('all 7 cross-page helpers exposed on window');

// ---- 2. Cluster registry lifecycle ----
const reg = win._xpEnsureClusterRegistry();
if (!reg.dosage || !reg.theta_pi || !reg.ghsl) fail('registry missing source slots');
if (typeof reg.dosage !== 'object') fail('dosage slot should be an object');
ok('registry has 3 source slots: dosage, theta_pi, ghsl');

// Register some labels
const N = 226;
const dosageLabels = new Array(N).fill(0).map((_, i) => i % 3);   // round-robin 0,1,2
const thetaLabels = new Array(N).fill(0).map((_, i) => i % 3);    // identical to dosage
const ghslLabels = new Array(N).fill(0).map((_, i) => (i + 1) % 3); // shifted by 1

win._xpRegisterClusters('dosage', 'cand_63', dosageLabels, 3);
win._xpRegisterClusters('theta_pi', 'cand_63', thetaLabels, 3);
win._xpRegisterClusters('ghsl', 'cand_63', ghslLabels, 3);

const dosageEntry = win._xpLookupClusters('dosage', 'cand_63');
if (!dosageEntry) fail('dosage entry should exist after register');
if (dosageEntry.k !== 3) fail('k should be 3, got ' + dosageEntry.k);
if (dosageEntry.labels.length !== N) fail('labels length should be ' + N);
ok('cluster registration: dosage[cand_63] = ' + N + ' labels, k=3');

// Bogus source → no-op
win._xpRegisterClusters('bogus_source', 'cand_63', dosageLabels, 3);
if (reg.bogus_source) fail('bogus source should not create slot');
ok('bogus source name is rejected (no slot created)');

// ---- 3. Active candidate id resolution ----
// Test all 4 fallback paths
win.state.focalCandidate = { id: 'focal_id_1' };
if (win._xpActiveCandidateId() !== 'focal_id_1') fail('focalCandidate path failed');
delete win.state.focalCandidate;

win.state.lockedCandidate = { id: 'locked_id_1' };
if (win._xpActiveCandidateId() !== 'locked_id_1') fail('lockedCandidate path failed');
delete win.state.lockedCandidate;

win.state.currentCandidate = { id: 'current_id_1' };
if (win._xpActiveCandidateId() !== 'current_id_1') fail('currentCandidate path failed');
delete win.state.currentCandidate;

win.state.focalEnvIdx = 17;
if (win._xpActiveCandidateId() !== 'env_17') fail('focalEnvIdx path failed');
delete win.state.focalEnvIdx;

if (win._xpActiveCandidateId() !== null) fail('no candidate path should give null');
ok('candidate id resolution: 4 state paths + null fallback');

// ---- 4. Color lookup for a registered candidate ----
win.state.focalCandidate = { id: 'cand_63' };
const c0 = win._xpSampleColor(0, 'dosage');
const c1 = win._xpSampleColor(1, 'dosage');
const c2 = win._xpSampleColor(2, 'dosage');
// Sample 0 → label 0 (round-robin), sample 1 → label 1, sample 2 → label 2
if (c0 !== '#4fa3ff') fail('sample 0 dosage color should be blue (#4fa3ff), got ' + c0);
if (c1 !== '#b8b8b8') fail('sample 1 dosage color should be neutral grey (#b8b8b8), got ' + c1);
if (c2 !== '#f5a524') fail('sample 2 dosage color should be amber (#f5a524), got ' + c2);
ok('color lookup: 3 cluster labels map to blue / grey / amber');

// Cross-page lookup: the same fish should get the same color in θπ (since
// dosage and theta_pi labels are identical).
const cThetaPi0 = win._xpSampleColor(0, 'theta_pi');
if (cThetaPi0 !== c0) fail('sample 0 should keep blue color in theta_pi (identical labels), got ' + cThetaPi0);
ok('cross-page consistency: label 0 → blue in dosage AND theta_pi');

// ---- 5. Color fallback when source has no clusters ----
win.state.focalCandidate = { id: 'unknown_candidate' };
const cFallback = win._xpSampleColor(0, 'dosage');
if (cFallback !== '#888') fail('unregistered candidate should give grey, got ' + cFallback);
ok('unregistered candidate: returns grey #888');

// No active candidate at all
delete win.state.focalCandidate;
const cNoActive = win._xpSampleColor(0, 'dosage');
if (cNoActive !== '#888') fail('no active candidate should give grey');
ok('no active candidate: returns grey #888');

// ---- 6. Concordance computation ----
win.state.focalCandidate = { id: 'cand_63' };

// dosage vs theta_pi (identical labels) → perfect agreement
const concDP = win._xpComputeConcordance('dosage', 'theta_pi');
if (!concDP) fail('concordance computation returned null');
if (concDP.n !== N) fail('concordance n should be ' + N);
if (concDP.agreement < 0.99) fail('identical labels should give agreement ~1.0, got ' + concDP.agreement);
if (concDP.ari < 0.99) fail('identical labels should give ARI ~1.0, got ' + concDP.ari);
if (concDP.cramers_v < 0.99) fail('identical labels should give Cramér V ~1.0, got ' + concDP.cramers_v);
ok('dosage vs theta_pi (identical): agreement=' + concDP.agreement.toFixed(3) +
   ', ARI=' + concDP.ari.toFixed(3) + ', Cramér V=' + concDP.cramers_v.toFixed(3));

// dosage vs ghsl (labels shifted by 1) → also high agreement (perfect 1-to-1
// remapping, just permuted) — best-match should still be 100%
const concDG = win._xpComputeConcordance('dosage', 'ghsl');
if (!concDG) fail('dosage vs ghsl returned null');
if (concDG.agreement < 0.99) fail('shifted labels (1-to-1 remap) should give perfect best-match, got ' + concDG.agreement);
ok('dosage vs ghsl (shifted by 1): best-match agreement=' + concDG.agreement.toFixed(3) + ' (Hungarian-resilient)');

// Contingency table is 3×3 with off-diagonal pattern (since shifted labels
// puts everything on the off-diagonal)
if (concDG.contingency.length !== 3) fail('contingency should be 3 rows');
if (concDG.contingency[0].length !== 3) fail('contingency should be 3 cols');
const onDiag = concDG.contingency[0][0] + concDG.contingency[1][1] + concDG.contingency[2][2];
const offDiag = N - onDiag;
if (offDiag <= onDiag) fail('shifted-labels contingency should be off-diagonal heavy');
ok('contingency: 3×3 table, off-diagonal=' + offDiag + ', diagonal=' + onDiag);

// ---- Random labels → low concordance ----
const randomLabels = new Array(N).fill(0).map(() => Math.floor(Math.random() * 3));
win._xpRegisterClusters('ghsl', 'cand_63', randomLabels, 3);
const concRandom = win._xpComputeConcordance('dosage', 'ghsl');
// With random labels we expect agreement somewhere around 1/3 (chance level)
if (concRandom.agreement > 0.6) fail('random labels should give chance-level agreement ~0.33-0.5, got ' + concRandom.agreement);
ok('random labels: agreement=' + concRandom.agreement.toFixed(3) + ' (~chance, < 0.6)');

// Restore good labels
win._xpRegisterClusters('ghsl', 'cand_63', ghslLabels, 3);

// ---- 7. getSampleColor dispatches to cluster_* modes ----
win.state.colorMode = 'cluster_dosage';
const cm0 = win.getSampleColor(0, 'cluster_dosage', null);
if (cm0 !== '#4fa3ff') fail('cluster_dosage mode for label 0 should give blue, got ' + cm0);
ok('getSampleColor("cluster_dosage"): dispatches to _xpSampleColor');

win.state.colorMode = 'cluster_theta_pi';
const cmTheta = win.getSampleColor(1, 'cluster_theta_pi', null);
if (cmTheta !== '#b8b8b8') fail('cluster_theta_pi mode for label 1 should give grey-neutral, got ' + cmTheta);
ok('getSampleColor("cluster_theta_pi"): dispatches correctly');

win.state.colorMode = 'cluster_ghsl';
// sample 0 in ghsl labels has label 1 (= (0+1) % 3)
const cmGhsl = win.getSampleColor(0, 'cluster_ghsl', null);
if (cmGhsl !== '#b8b8b8') fail('cluster_ghsl mode for sample 0 (label 1) should give neutral grey, got ' + cmGhsl);
ok('getSampleColor("cluster_ghsl"): dispatches correctly');

// ---- 8. UI buttons added ----
const colorModeBar = win.document.getElementById('colorModeBar');
if (!colorModeBar) fail('colorModeBar missing');
const residualBtn = colorModeBar.querySelector('button[data-mode="residual"]');
if (!residualBtn) fail('residual button missing from colorModeBar');
ok('colorModeBar: residual button present (from turn 83)');

// ---- 9. Repeat element class label hover ----
// Search for the source HTML rather than the rendered DOM (boundaries page is
// dynamically rendered, so the dropdown only exists when the page is visited).
if (!html.includes('TE superfamily')) fail('repeat-element class hover should mention TE superfamily');
ok('repeat element class column: hover description mentions TE superfamily');

console.log('\n[xp] ALL CHECKS PASSED');
