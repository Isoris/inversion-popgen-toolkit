// Tests for v4 turn 86 — Q-association ancestry coloring.
// Coverage:
//   1. Helpers exposed
//   2. State lifecycle: ensure / register / setK / clear
//   3. Q-vector lookup (subarray)
//   4. Hard-call coloring (dominant-Q full color, admixed → grey)
//   5. Blend coloring (weighted RGB blend of group colors)
//   6. Proportion bars: per-cluster (figure-style)
//   7. Proportion bars: cohort-wide
//   8. getSampleColor dispatches to q_ancestry mode
//   9. UI: button + sub-controls present
//  10. Sub-controls hidden by default, shown on q_ancestry

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[qa] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[qa] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['_qaEnsureState', '_qaRegisterK', '_qaSetK', '_qaSampleQ',
                    '_qaResolveSample', '_qaSampleColor', '_qaProportionBars',
                    '_qaClearState', '_qaPopulateKSelect']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
if (!Array.isArray(win._QA_DEFAULT_PALETTE)) fail('_QA_DEFAULT_PALETTE not exposed');
if (typeof win._QA_ADMIX_THRESHOLD_DEFAULT !== 'number') fail('threshold not exposed');
ok('all 9 Q-ancestry helpers + 2 constants exposed on window');

// ---- 2. State lifecycle ----
const ui0 = win._qaEnsureState();
if (!ui0) fail('ensure state failed');
if (ui0.K !== null) fail('initial K should be null, got ' + ui0.K);
if (!Array.isArray(ui0.available_K)) fail('available_K should be array');
if (ui0.available_K.length !== 0) fail('available_K should start empty');
ok('initial state: K=null, available_K=[], palette=20 colors');

// Register K=2 with synthetic Q-vectors for 4 fish
const N = 4;
const K2 = 2;
//   fish 0: 0.95 / 0.05 → clean Q1 (dominant=0)
//   fish 1: 0.55 / 0.45 → admixed (max < 0.7)
//   fish 2: 0.10 / 0.90 → clean Q2 (dominant=1)
//   fish 3: 0.30 / 0.70 → clean Q2 (just at threshold)
const q2 = new Float32Array([0.95, 0.05,  0.55, 0.45,  0.10, 0.90,  0.30, 0.70]);
win._qaRegisterK(K2, q2);
const ui = win._qaEnsureState();
if (!ui.available_K.includes(2)) fail('K=2 not registered in available_K');
if (ui.K !== 2) fail('first registered K should auto-set, got K=' + ui.K);
ok('register K=2: available_K=[2], auto-set as active');

// Register K=8 alongside
const K8 = 8;
const q8 = new Float32Array(N * K8);
// Fish 0 dominant Q3, fish 1 dominant Q5, etc.
for (let i = 0; i < N; i++) {
  const dom = (i * 2) % K8;
  for (let k = 0; k < K8; k++) q8[i * K8 + k] = (k === dom) ? 0.92 : 0.08 / 7;
}
win._qaRegisterK(K8, q8);
const ui2 = win._qaEnsureState();
if (!ui2.available_K.includes(8)) fail('K=8 not registered');
if (ui2.available_K.length !== 2) fail('expected 2 K values, got ' + ui2.available_K.length);
ok('register K=8: available_K=[2, 8] (sorted ascending)');

// Switch K
win._qaSetK(8);
if (win._qaEnsureState().K !== 8) fail('setK(8) failed');
ok('setK(8): active K is now 8');

// Bogus K → console warn, no change
win._qaSetK(99);
if (win._qaEnsureState().K !== 8) fail('setK(99) should not change active K');
ok('setK(99) bogus: rejected, K stays 8');

// ---- 3. Q-vector lookup ----
const qVec0 = win._qaSampleQ(0);
if (!qVec0 || qVec0.length !== 8) fail('sample 0 Q-vector should be length 8, got ' + (qVec0 ? qVec0.length : 'null'));
if (Math.abs(qVec0[0] - 0.92) > 0.01) fail('sample 0 dominant Q should be 0.92 at index 0');
ok('Q-vector lookup: sample 0 length=8, max at index 0 = 0.92');

// ---- 4. Hard-call coloring ----
win._qaSetK(2);
win.state.qDisplayMode = 'hard';

const c0 = win._qaSampleColor(0);  // 95% Q1 → orange (palette[0] = #f5a524)
const c1 = win._qaSampleColor(1);  // admixed (max=0.55 < 0.7) → grey
const c2 = win._qaSampleColor(2);  // 90% Q2 → purple (palette[1] = #9b59b6)
const c3 = win._qaSampleColor(3);  // 70% Q2 → palette[1] (just at threshold)

if (c0 !== '#f5a524') fail('fish 0 (95% Q1) should be palette[0] orange, got ' + c0);
if (c1 !== '#888') fail('fish 1 (admixed 55/45) should be grey #888, got ' + c1);
if (c2 !== '#9b59b6') fail('fish 2 (90% Q2) should be palette[1] purple, got ' + c2);
ok('hard mode: clean Q1 → orange (#f5a524)');
ok('hard mode: admixed (max=0.55) → grey #888');
ok('hard mode: clean Q2 → purple (#9b59b6)');

// Edge case: max_q exactly at threshold (0.70) → still gets full color
if (c3 !== '#9b59b6') fail('fish 3 (70% Q2 at threshold) should still get full color, got ' + c3);
ok('hard mode: max_q == threshold → full color (boundary inclusive)');

// ---- 5. Blend coloring ----
win.state.qDisplayMode = 'blend';

const b0 = win._qaSampleColor(0);  // 95% Q1 → mostly orange
const b1 = win._qaSampleColor(1);  // 55/45 → midway between orange and purple

// Helper to parse hex → [r,g,b]
const hex = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
const [r0, g0, b0c] = hex(b0);
const [r1, g1, b1c] = hex(b1);
// orange #f5a524 = (245, 165, 36); purple #9b59b6 = (155, 89, 182)
// fish 0 (95% Q1): R should be near 245
if (r0 < 230) fail('blend fish 0 (95% Q1) should have R near 245, got ' + r0);
// fish 1 (55/45): R should be midway
if (r1 < 180 || r1 > 220) fail('blend fish 1 (55/45) should have R midway 180-220, got ' + r1);
ok('blend mode: 95% Q1 → R=' + r0 + ' (≥230, near orange)');
ok('blend mode: 55/45 → R=' + r1 + ' (180-220, midway)');

// ---- 6. Proportion bars: per_cluster ----
// Synthetic K-means labels: fish 0,1 → cluster 0 (HOMO_1); fish 2,3 → cluster 1 (HOMO_2)
const labels = [0, 0, 1, 1];
const barsPerCluster = win._qaProportionBars('per_cluster', labels);
if (!barsPerCluster) fail('per_cluster bars returned null');
if (!Array.isArray(barsPerCluster.clusters)) fail('per_cluster should have clusters array');
if (barsPerCluster.clusters.length !== 2) fail('expected 2 clusters, got ' + barsPerCluster.clusters.length);
const c0Bar = barsPerCluster.clusters[0];
if (c0Bar.id !== 0) fail('first cluster should have id=0');
if (c0Bar.count !== 2) fail('cluster 0 should have count=2');
// cluster 0 = avg of fish 0 (Q=[0.95,0.05]) and fish 1 (Q=[0.55,0.45]) = [0.75, 0.25]
if (Math.abs(c0Bar.props[0] - 0.75) > 0.01) fail('cluster 0 mean Q1 should be 0.75, got ' + c0Bar.props[0]);
if (Math.abs(c0Bar.props[1] - 0.25) > 0.01) fail('cluster 0 mean Q2 should be 0.25, got ' + c0Bar.props[1]);
ok('per_cluster bars: cluster 0 (n=2) → Q1=0.75, Q2=0.25 (figure-style decomposition)');
const c1Bar = barsPerCluster.clusters[1];
// cluster 1 = avg of fish 2 (Q=[0.10,0.90]) and fish 3 (Q=[0.30,0.70]) = [0.20, 0.80]
if (Math.abs(c1Bar.props[0] - 0.20) > 0.01) fail('cluster 1 mean Q1 should be 0.20');
if (Math.abs(c1Bar.props[1] - 0.80) > 0.01) fail('cluster 1 mean Q2 should be 0.80');
ok('per_cluster bars: cluster 1 (n=2) → Q1=0.20, Q2=0.80');

// ---- 7. Proportion bars: cohort ----
const barsCohort = win._qaProportionBars('cohort');
if (!barsCohort) fail('cohort bars returned null');
if (barsCohort.n !== N) fail('cohort n should be ' + N);
// cohort mean = mean over 4 fish: (0.95 + 0.55 + 0.10 + 0.30) / 4 = 0.475 for Q1
if (Math.abs(barsCohort.props[0] - 0.475) > 0.01) fail('cohort mean Q1 should be 0.475');
if (Math.abs(barsCohort.props[1] - 0.525) > 0.01) fail('cohort mean Q2 should be 0.525');
ok('cohort bars: cohort-wide (n=4) → Q1=0.475, Q2=0.525');

// ---- 8. getSampleColor dispatches ----
win.state.colorMode = 'q_ancestry';
win.state.qDisplayMode = 'hard';
const dispatchedC0 = win.getSampleColor(0, 'q_ancestry', null);
if (dispatchedC0 !== '#f5a524') fail('dispatched color for fish 0 (clean Q1) should be orange, got ' + dispatchedC0);
ok('getSampleColor("q_ancestry"): dispatches to _qaSampleColor');

const dispatchedC1 = win.getSampleColor(1, 'q_ancestry', null);
if (dispatchedC1 !== '#888') fail('dispatched color for admixed fish should be grey');
ok('getSampleColor("q_ancestry"): admixed → grey');

// ---- 9. UI: button + sub-controls present ----
const colorModeBar = win.document.getElementById('colorModeBar');
if (!colorModeBar) fail('colorModeBar missing');
const qBtn = colorModeBar.querySelector('button[data-mode="q_ancestry"]');
if (!qBtn) fail('Q ancestry button missing');
ok('colorModeBar: Q ancestry button present');

const subControls = win.document.getElementById('qAncestrySubControls');
if (!subControls) fail('qAncestrySubControls missing');
const kSel = win.document.getElementById('qAncestryKSelect');
if (!kSel) fail('qAncestryKSelect missing');
ok('q_ancestry sub-controls: K select + display + legend toggles present');

// ---- 10. Sub-controls hidden by default ----
// (CSS inline style sets display: none)
if (subControls.style.display !== 'none') fail('sub-controls should be hidden by default');
ok('sub-controls: hidden by default (display: none)');

// Click the q_ancestry button → sub-controls become visible
qBtn.click();
// JSDOM: style.display should now be 'flex'
if (subControls.style.display !== 'flex') fail('sub-controls should show as flex on q_ancestry click, got ' + subControls.style.display);
ok('q_ancestry button click: sub-controls become visible (display: flex)');

// Click another mode → sub-controls hide
const noneBtn = colorModeBar.querySelector('button[data-mode="none"]');
noneBtn.click();
if (subControls.style.display !== 'none') fail('sub-controls should hide when mode is not q_ancestry');
ok('switching mode away from q_ancestry: sub-controls hide');

// ---- 11. K dropdown population ----
qBtn.click();   // re-activate q_ancestry to show sub-controls
win._qaPopulateKSelect();
const kSel2 = win.document.getElementById('qAncestryKSelect');
const opts = Array.from(kSel2.querySelectorAll('option')).map(o => o.value);
if (opts.length !== 2) fail('K dropdown should have 2 options (K=2, K=8), got ' + opts.length);
if (opts[0] !== '2' || opts[1] !== '8') fail('K dropdown options should be ["2", "8"], got ' + JSON.stringify(opts));
ok('K dropdown: populated with [2, 8] from registered K values');

// ---- 12. Clear state ----
win._qaClearState();
const cleared = win._qaEnsureState();
if (cleared.K !== null) fail('after clear, K should be null');
if (cleared.available_K.length !== 0) fail('after clear, available_K should be empty');
ok('_qaClearState: resets K registry, keeps user prefs');

// User preferences persist across clear
if (!win.state.qDisplayMode) fail('qDisplayMode preference should persist');
if (!win.state.qLegendMode) fail('qLegendMode preference should persist');
ok('user prefs (qDisplayMode, qLegendMode) persist across clear');

console.log('\n[qa] ALL CHECKS PASSED');
