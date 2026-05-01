// Tests for v4 turn 103 — band-reach + regime-breadth.
//
// Coverage:
//   1. Helpers exposed
//   2. Empty state: no L2 envelopes → empty result
//   3. Single L2 (all fish in one band each) → all reach=1
//   4. 3-L2 chain with stable lanes → narrow regime breadth
//   5. 3-L2 chain with chaotic crossings → wide regime breadth
//   6. Multi-band populated: bands_populated counts bands with >=5 fish
//   7. computeWindowedBandReachPerL2 returns one entry per L2

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[breach] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[breach] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
const required = ['computeBandReachAcrossL2s', 'computeWindowedBandReachPerL2'];
for (const name of required) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
ok('all 2 band-reach helpers exposed');

// ---- 2. Empty state ----
win.state.data = null;
const r1 = win.computeBandReachAcrossL2s([0]);
if (r1 !== null) fail('empty state should return null, got ' + JSON.stringify(r1));
ok('empty state: returns null');

const r2 = win.computeWindowedBandReachPerL2();
if (!Array.isArray(r2) || r2.length !== 0) fail('empty state should return []');
ok('computeWindowedBandReachPerL2: returns [] on empty state');

// ---- 3. Mock state.data + getL2Cluster + alignLabels ----
const N = 30;

// Build label sets:
//   L2[0]: 10 in band 0, 10 in band 1, 10 in band 2 (clean lanes)
//   L2[1]: same as L2[0] (each fish stays in its lane → narrow)
//   L2[2]: same as L2[0] (still narrow)
//   L2[3]: random shuffling (50% chance each fish in any band → wide)
const labels0 = new Int8Array(N);
const labels1 = new Int8Array(N);
const labels2 = new Int8Array(N);
for (let i = 0; i < 10; i++) { labels0[i] = 0; labels1[i] = 0; labels2[i] = 0; }
for (let i = 10; i < 20; i++) { labels0[i] = 1; labels1[i] = 1; labels2[i] = 1; }
for (let i = 20; i < 30; i++) { labels0[i] = 2; labels1[i] = 2; labels2[i] = 2; }

// L2[3]: every fish gets band that's NOT its labels2 band
//   For labels2=0: shift to band 1
//   For labels2=1: shift to band 2
//   For labels2=2: shift to band 0
const labels3 = new Int8Array(N);
for (let i = 0; i < N; i++) labels3[i] = (labels2[i] + 1) % 3;

// L2[4]: every fish gets the remaining band (the one not in labels2 or labels3)
//   labels2 + labels3 + labels4 = {0, 1, 2}
const labels4 = new Int8Array(N);
for (let i = 0; i < N; i++) labels4[i] = (labels2[i] + 2) % 3;
// Result: every fish visits all 3 bands across [labels2, labels3, labels4]

const fakeClusters = [
  { fixedKLabels: labels0, labels: labels0, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { fixedKLabels: labels1, labels: labels1, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { fixedKLabels: labels2, labels: labels2, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { fixedKLabels: labels3, labels: labels3, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { fixedKLabels: labels4, labels: labels4, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
];
win.getL2Cluster = (l2idx) => fakeClusters[l2idx] || null;
// Mock alignLabels to be identity (perm = [0, 1, 2]) — labels are already aligned
win.alignLabels = (a, b, K) => ({ perm: [0, 1, 2], concord: 1.0 });

win.state.data = {
  n_samples: N,
  l2_envelopes: [
    { start_bp: 1000000,  end_bp: 2000000 },
    { start_bp: 2000001,  end_bp: 3000000 },
    { start_bp: 3000001,  end_bp: 4000000 },
    { start_bp: 4000001,  end_bp: 5000000 },
    { start_bp: 5000001,  end_bp: 6000000 },
  ],
};
win.state.k = 3;

// ---- 3. Single L2: all fish reach=1 ----
const r3 = win.computeBandReachAcrossL2s([0]);
if (!r3) fail('single L2 should return result');
if (r3.K !== 3) fail('K should be 3');
let allReachOne = true;
for (let s = 0; s < N; s++) {
  if (r3.per_sample_band_reach[s] !== 1) { allReachOne = false; break; }
}
if (!allReachOne) fail('single L2: every fish should have band_reach=1');
ok('single L2: every fish reach=1, narrow_fraction=' + r3.narrow_fraction.toFixed(2));
if (r3.regime_breadth !== 'narrow') fail('single L2 with 3 bands populated should be narrow regime, got ' + r3.regime_breadth);
ok('single L2 with 3 bands populated: regime_breadth = "narrow"');

// ---- 4. 3-L2 chain with stable lanes ----
const r4 = win.computeBandReachAcrossL2s([0, 1, 2]);
if (!r4) fail('3-L2 chain returned null');
let stillNarrow = true;
for (let s = 0; s < N; s++) {
  if (r4.per_sample_band_reach[s] !== 1) { stillNarrow = false; break; }
}
if (!stillNarrow) fail('3 L2s with identical labels: every fish should still have reach=1');
ok('3-L2 chain (identical labels): every fish reach=1');
if (r4.regime_breadth !== 'narrow') fail('expected narrow, got ' + r4.regime_breadth);
ok('3-L2 chain (stable lanes): regime_breadth = "narrow"');

// ---- 5. 3-L2 chain with chaotic crossings (L2[2], L2[3], L2[4] all shuffle) ----
// labels2 = original lanes
// labels3 = (i+1) % 3 — every fish moved one band
// labels4 = (i+2) % 3 — every fish moved another band
// So every fish visits 3 distinct bands across this chain
const r5 = win.computeBandReachAcrossL2s([2, 3, 4]);
if (!r5) fail('chaotic chain returned null');
let avgReach = 0;
for (let s = 0; s < N; s++) avgReach += r5.per_sample_band_reach[s];
avgReach /= N;
if (avgReach < 2.9) fail('chaotic chain: avg reach should be ~3.0 (every fish visits all bands), got ' + avgReach.toFixed(2));
ok('3-L2 chain (chaotic): avg reach=' + avgReach.toFixed(2) + ' (every fish visits all 3 bands), narrow_frac=' + r5.narrow_fraction.toFixed(2));
if (r5.regime_breadth !== 'wide') fail('expected wide, got ' + r5.regime_breadth);
ok('3-L2 chain (chaotic crossings): regime_breadth = "wide"');

// ---- 6. bands_populated count ----
if (r5.bands_populated !== 3) fail('bands_populated should be 3 (all 3 bands have >=5 fish), got ' + r5.bands_populated);
ok('bands_populated: 3 (all bands have ≥5 fish visiting)');

// Per-band visit counts: with 30 fish each visiting 3 bands across this chain,
// each band should have 30 visitors (everyone visits every band)
let allPopulated = true;
for (let k = 0; k < 3; k++) {
  if (r5.per_band_visit_count[k] !== 30) {
    allPopulated = false;
    fail('band ' + k + ' visit count should be 30, got ' + r5.per_band_visit_count[k]);
  }
}
ok('per_band_visit_count: each of 3 bands has 30 visiting fish (chaotic chain)');

// ---- 7. computeWindowedBandReachPerL2 ----
const windowed = win.computeWindowedBandReachPerL2();
if (!Array.isArray(windowed) || windowed.length !== 5) fail('expected 5 entries (one per L2), got ' + windowed.length);
ok('computeWindowedBandReachPerL2: returns 5 entries (one per L2 envelope)');
for (const w of windowed) {
  if (typeof w.center_l2 !== 'number') fail('windowed entry missing center_l2');
  if (!isFinite(w.center_mb)) fail('windowed entry missing center_mb');
}
ok('windowed entries: each has center_l2 + center_mb populated');

// L2[0] is at the boundary (only L2[0] and L2[1] in window) — should still
// be narrow because labels0 = labels1
const w0 = windowed[0];
if (w0.regime_breadth !== 'narrow') fail('L2[0] should be narrow (stable lanes), got ' + w0.regime_breadth);
ok('windowed L2[0]: narrow (stable lanes with neighbor)');

// L2[3] is in the chaotic zone (window covers L2[2], L2[3], L2[4])
const w3 = windowed[3];
if (w3.regime_breadth !== 'wide') fail('L2[3] should be wide, got ' + w3.regime_breadth);
ok('windowed L2[3]: wide (chaotic neighbors)');

console.log('\n[breach] ALL CHECKS PASSED');
