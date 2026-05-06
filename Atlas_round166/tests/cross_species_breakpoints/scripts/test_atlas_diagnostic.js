// Tests for v4 turn 83 diagnostic mode (residual coloring + suspicion scoring).
// Coverage:
//   1. Helpers exposed
//   2. Residual color ramp: blue at z=0, amber mid, red at z>=3
//   3. Per-window residual computation (atlas-side from synthetic PC1+PC2)
//   4. Per-window residual computation (precomp-provided path)
//   5. Per-candidate suspicion summary (mean, max, consistency, suspicious flag)
//   6. getSampleColor dispatches to residual mode when state.colorMode === 'residual'
//   7. Cache invalidation

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[diag] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[diag] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['_diagComputeWindowResiduals', '_diagComputeCandidateSuspicion',
                    '_diagResidualColor', '_diagSampleColor', '_diagClearCache']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
if (typeof win._DIAG_RESIDUAL_SUSPICIOUS_Z !== 'number') fail('_DIAG_RESIDUAL_SUSPICIOUS_Z not exposed');
if (typeof win._DIAG_BAND_CONSISTENCY_THRESH !== 'number') fail('_DIAG_BAND_CONSISTENCY_THRESH not exposed');
ok('all 5 diagnostic helpers + 2 thresholds exposed on window');

// ---- 2. Residual color ramp ----
const colZ0 = win._diagResidualColor(0);
const colZmid = win._diagResidualColor(2);
const colZhigh = win._diagResidualColor(4);
// At z=0 we should be in blue territory: B > R, B > G
const hex = h => [parseInt(h.slice(1,3),16), parseInt(h.slice(3,5),16), parseInt(h.slice(5,7),16)];
const [r0, g0, b0] = hex(colZ0);
const [r4, g4, b4] = hex(colZhigh);
if (b0 < r0) fail('z=0 should be blue-dominant, got R='+r0+' B='+b0);
if (r4 < b4) fail('z>=4 should be red-dominant, got R='+r4+' B='+b4);
ok('residual color ramp: z=0 → blue (' + colZ0 + '), z=4 → red (' + colZhigh + ')');

// Null/invalid → grey fallback
if (win._diagResidualColor(null) !== '#888') fail('null should give grey');
if (win._diagResidualColor(NaN) !== '#888') fail('NaN should give grey');
ok('residual color: null/NaN → #888 (grey fallback)');

// ---- 3. Per-window residual computation, atlas-side ----
// Construct synthetic state.data: 2 windows with 3 K-means bands of ~10 fish each.
// Band 1 fish at PC1≈-2, Band 2 at 0, Band 3 at +2. Add ONE outlier at +5 in band 3.
const N = 30;
const pc1 = new Array(N).fill(0).map((_, i) => {
  if (i < 10) return -2 + (Math.random() - 0.5) * 0.3;       // band 1
  if (i < 20) return  0 + (Math.random() - 0.5) * 0.3;       // band 2
  if (i === 25) return 8;                                    // outlier (band 3)
  return 2 + (Math.random() - 0.5) * 0.3;                    // band 3
});
const pc2 = new Array(N).fill(0).map(() => (Math.random() - 0.5) * 0.5);

win.state.data = {
  n_samples: N,
  windows: [
    { pc1: pc1.slice(), pc2: pc2.slice() },
    { pc1: pc1.slice(), pc2: pc2.slice() },
  ],
};
win.state.cur = 0;
win.state.k = 3;
win._diagClearCache();

const resW0 = win._diagComputeWindowResiduals(0);
if (!resW0) fail('atlas-computed residual returned null');
if (resW0.source !== 'atlas') fail('expected source=atlas, got ' + resW0.source);
if (resW0.residuals.length !== N) fail('residuals length mismatch');
ok('atlas-side residual computation: window 0, n=' + N + ', source=atlas');

// Outlier fish should have a high residual; deep-cluster fish should have low
const outlierZ = resW0.residuals[25];
const insiderZ = resW0.residuals[0];   // first band-1 fish, should be deep inside
if (outlierZ < 2) fail('outlier (PC1=5 against band-3 centroid ~2) should have high residual_z, got ' + outlierZ);
if (insiderZ > 2) fail('insider should have low residual_z, got ' + insiderZ);
ok('outlier residual_z=' + outlierZ.toFixed(2) + ' (high), insider=' + insiderZ.toFixed(2) + ' (low)');

// ---- 4. Per-window residual: precomp-provided path ----
// If state.data.windows[i].band_residual_z is present, use it directly
const precompResid = new Array(N).fill(0).map((_, i) => i === 5 ? 4.5 : 0.5);
const precompBand = new Array(N).fill(0).map((_, i) => i < 10 ? 0 : (i < 20 ? 1 : 2));
win.state.data.windows[0].band_residual_z = precompResid;
win.state.data.windows[0].band = precompBand;
win._diagClearCache();
const rPre = win._diagComputeWindowResiduals(0);
if (rPre.source !== 'precomp') fail('expected source=precomp when band_residual_z present');
if (Math.abs(rPre.residuals[5] - 4.5) > 0.01) fail('precomp residual not propagated, got ' + rPre.residuals[5]);
ok('precomp-provided path: source=precomp, residual values pass through unchanged');

// Clean up for next tests
delete win.state.data.windows[0].band_residual_z;
delete win.state.data.windows[0].band;
win._diagClearCache();

// ---- 5. Per-candidate suspicion summary ----
const summary = win._diagComputeCandidateSuspicion(0, 1, 3);
if (!summary) fail('candidate suspicion summary returned null');
if (summary.meanZ.length !== N) fail('meanZ length wrong');
if (summary.maxZ.length !== N) fail('maxZ length wrong');
if (summary.consistency.length !== N) fail('consistency length wrong');
if (summary.suspicious.length !== N) fail('suspicious length wrong');
if (summary.n_windows !== 2) fail('expected n_windows=2, got ' + summary.n_windows);
ok('candidate suspicion summary: meanZ + maxZ + consistency + suspicious flag, n_windows=2');

// Outlier (idx 25) should be flagged suspicious
if (!summary.suspicious[25]) fail('outlier (idx 25) should be flagged suspicious');
ok('outlier idx 25: suspicious=true (maxZ=' + summary.maxZ[25].toFixed(2) + ')');

// Insider (idx 0) should NOT be flagged
if (summary.suspicious[0]) fail('insider (idx 0) should NOT be flagged suspicious');
ok('insider idx 0: suspicious=false (maxZ=' + summary.maxZ[0].toFixed(2) + ')');

// ---- 6. getSampleColor dispatches to residual mode ----
win.state.colorMode = 'cluster';
win.state.cur = 0;
win._diagClearCache();
// In cluster mode without groupLabels argument, fallback grey
const clusterCol = win.getSampleColor(25, 'cluster', null);
if (clusterCol !== '#888') fail('cluster mode without labels should give grey, got ' + clusterCol);
// In residual mode, the outlier should get a red-ish color
const residualCol = win.getSampleColor(25, 'residual', null);
if (residualCol === '#888') fail('residual mode should give residual-based color, not grey');
const [rr, gg, bb] = hex(residualCol);
if (rr <= bb) fail('outlier residual color should be red-dominant, got ' + residualCol);
ok('getSampleColor("residual"): outlier dot colored red (' + residualCol + ')');

// Insider in residual mode should be blue-ish
const insiderCol = win.getSampleColor(0, 'residual', null);
const [ir, ig, ib] = hex(insiderCol);
if (ib <= ir) fail('insider residual color should be blue-dominant, got ' + insiderCol);
ok('getSampleColor("residual"): insider dot colored blue (' + insiderCol + ')');

// ---- 7. Cache invalidation ----
win._diagComputeWindowResiduals(0);
if (Object.keys(win.state._diagCache).length === 0) fail('cache should populate after compute');
win._diagClearCache();
if (Object.keys(win.state._diagCache).length !== 0) fail('clearCache should empty cache');
ok('_diagClearCache: empties the residual cache');

// ---- 8. Edge cases ----
// Empty windows array
win.state.data = { n_samples: 0, windows: [] };
win._diagClearCache();
const empty = win._diagComputeWindowResiduals(0);
if (empty !== null) fail('empty window data should return null');
ok('empty data: returns null gracefully');

// Missing pc1/pc2 in window
win.state.data = { n_samples: 5, windows: [{}] };
win._diagClearCache();
const noPc = win._diagComputeWindowResiduals(0);
if (noPc !== null) fail('missing pc1/pc2 should return null');
ok('missing pc1/pc2: returns null gracefully');

console.log('\n[diag] ALL CHECKS PASSED');
