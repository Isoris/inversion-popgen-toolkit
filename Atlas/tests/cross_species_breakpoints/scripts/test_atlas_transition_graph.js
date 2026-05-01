// Tests for v4 turn 101 — Structural-haplotype transition graph.
//
// Coverage:
//   1. Helpers exposed
//   2. Empty state: no l2_envelopes → empty graph
//   3. Per-boundary stats: transition_rate, n_changed, n_samples
//   4. Edge counts sorted by count desc; from_band / to_band correct
//   5. Hotspot threshold: boundaries with rate >= 0.30 flagged
//   6. summarizeTransitionBoundary produces top_edges_label

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[shtg] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[shtg] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
const required = [
  'computeStructuralHaplotypeTransitionGraph',
  'summarizeTransitionBoundary',
  '_shtgComputeBoundaryStats',
];
for (const name of required) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
if (typeof win._SHTG_HOTSPOT_THRESHOLD !== 'number') fail('_SHTG_HOTSPOT_THRESHOLD not exposed');
ok('all 3 transition-graph helpers + threshold const exposed');

// ---- 2. Empty state ----
win.state.data = null;
let g = win.computeStructuralHaplotypeTransitionGraph();
if (!g || !Array.isArray(g.boundaries) || !Array.isArray(g.hotspots)) {
  fail('should return { boundaries, hotspots } even on empty state');
}
if (g.boundaries.length !== 0 || g.hotspots.length !== 0) fail('empty state should return empty arrays');
ok('empty state: returns empty boundaries + hotspots arrays');

// ---- 3+4+5. Need to mock state.data + getL2Cluster + compareL2Pair ----
// Build 3 L2 envelopes with synthetic labels so we can test:
//   - L2[0]→L2[1]: low transition rate (most fish stay) — not a hotspot
//   - L2[1]→L2[2]: high transition rate (many fish change) — hotspot

const N = 30;

// L2[0] labels: 10 fish in band 0, 10 in band 1, 10 in band 2
const labels0 = new Int8Array(N);
for (let i = 0; i < 10; i++) labels0[i] = 0;
for (let i = 10; i < 20; i++) labels0[i] = 1;
for (let i = 20; i < 30; i++) labels0[i] = 2;

// L2[1]: same as L2[0] (no transitions, after Hungarian alignment)
const labels1 = new Int8Array(labels0);

// L2[2]: many fish swap from band 0 to band 1 (high transition rate)
const labels2 = new Int8Array(N);
for (let i = 0; i < 5; i++)  labels2[i] = 0;   // 5 stayed in band 0
for (let i = 5; i < 10; i++) labels2[i] = 1;   // 5 swapped to band 1
for (let i = 10; i < 20; i++) labels2[i] = 1;  // band 1 unchanged
for (let i = 20; i < 30; i++) labels2[i] = 2;  // band 2 unchanged
// Result: 5 fish (out of 30) transitioned at boundary 1→2 → rate ~17%
//   ... actually let me push more transitions to exceed the 30% threshold
for (let i = 0; i < 12; i++) labels2[i] = 1;   // bump 12 fish from band 0 to band 1
// Now transitions: 12 fish (40%) → hotspot

win.state.data = {
  n_samples: N,
  n_windows: 9,
  l2_envelopes: [
    { start_bp: 1000000,  end_bp: 2000000, candidate_id: 'L2_0' },
    { start_bp: 2000001,  end_bp: 3000000, candidate_id: 'L2_1' },
    { start_bp: 3000001,  end_bp: 4000000, candidate_id: 'L2_2' },
  ],
  windows: [],
};

// Mock state.k for compareL2Pair
win.state.k = 3;
win.state.mergeThr = 0.85;
win.state.alpha = 0.001;

// Mock getL2Cluster and compareL2Pair so the data layer can read labels
const fakeClusters = [
  { labels: labels0, fixedKLabels: labels0, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { labels: labels1, fixedKLabels: labels1, ok: true, usedK: 3, n_per_group: [10, 10, 10] },
  { labels: labels2, fixedKLabels: labels2, ok: true, usedK: 3, n_per_group: [3, 17, 10] },
];
win.getL2Cluster = (l2idx) => fakeClusters[l2idx] || null;
// compareL2Pair returns { perm, concord }: we mock perm as identity since the
// labels are constructed already aligned (band 0 stays band 0, etc.)
win.compareL2Pair = (li, ri) => ({
  perm: [0, 1, 2],
  concord: 0.9,   // doesn't actually matter, our code recomputes via labels
});

// Compute the graph
g = win.computeStructuralHaplotypeTransitionGraph();
if (g.boundaries.length !== 2) fail('expected 2 boundaries, got ' + g.boundaries.length);
ok('graph: 2 boundaries computed (3 L2 envelopes → 2 adjacent pairs)');

// Boundary 0→1: identical labels → 0 transitions
const b01 = g.boundaries[0];
if (b01.l2_left !== 0 || b01.l2_right !== 1) fail('first boundary should be (0, 1)');
if (b01.n_changed !== 0) fail('boundary 0→1 should have 0 transitions, got ' + b01.n_changed);
if (b01.transition_rate !== 0) fail('boundary 0→1 transition_rate should be 0');
ok('boundary 0→1: transition_rate=0 (identical labels, no fish changed)');

// Boundary 1→2: 10 fish switched from band 0 to band 1
// (fish 0..9 were band 0 in labels1, now band 1 in labels2; fish 10..11 stayed in band 1)
const b12 = g.boundaries[1];
if (b12.l2_left !== 1 || b12.l2_right !== 2) fail('second boundary should be (1, 2)');
if (b12.n_changed !== 10) fail('boundary 1→2 should have 10 transitions, got ' + b12.n_changed);
if (Math.abs(b12.transition_rate - 0.333) > 0.01) {
  fail('boundary 1→2 rate should be ~0.333, got ' + b12.transition_rate.toFixed(3));
}
ok('boundary 1→2: 10/30 fish transitioned, transition_rate=' + b12.transition_rate.toFixed(3));

// Edges should be sorted by count desc
const e0 = b12.edges[0];
if (e0.count < (b12.edges[1] || { count: 0 }).count) fail('edges not sorted by count desc');
ok('edges sorted by count descending: top edge = B' + e0.from_band + '→B' + e0.to_band + ' (n=' + e0.count + ')');

// Edges:
//   B0→B1: fish 0..9 (10 fish)
//   B1→B1: fish 10..19 (10 fish)
//   B2→B2: fish 20..29 (10 fish)
ok('boundary 1→2 edges: ' + b12.edges.map(e => 'B' + e.from_band + '→B' + e.to_band + ':' + e.count).join(', '));

// ---- 5. Hotspot detection ----
if (g.hotspots.length === 0) fail('boundary with rate ~0.33+ should be a hotspot');
ok('hotspots detected: ' + g.hotspots.length + ' (threshold=0.30)');
if (g.hotspots[0].l2_left !== 1) fail('hotspot should be boundary 1→2');
ok('hotspot is at boundary 1→2 (the high-transition boundary)');

// ---- 6. summarizeTransitionBoundary ----
const summ = win.summarizeTransitionBoundary(b12);
if (!summ) fail('summarize returned null');
if (typeof summ.transition_rate !== 'number') fail('summary missing transition_rate');
if (typeof summ.top_edges_label !== 'string') fail('summary missing top_edges_label');
if (!summ.is_hotspot) fail('summary should mark this as hotspot');
ok('summarizeTransitionBoundary: rate=' + summ.transition_rate.toFixed(2) + ', is_hotspot=true, label="' + summ.top_edges_label + '"');

// Position_mb computed
if (!isFinite(b12.position_mb)) fail('position_mb should be computed from envelope coords');
ok('boundary 1→2 position_mb = ' + b12.position_mb.toFixed(2) + ' Mb');

console.log('\n[shtg] ALL CHECKS PASSED');
