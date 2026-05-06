// Tests for v4 turn 92 — Diamond detector algorithm.
// Coverage:
//   1. Helpers exposed
//   2. No diamond on flat 3-band data → empty array
//   3. Loose diamond detected when one band splits
//   4. Strict diamond requires ≥1 stable band
//   5. Strict2 diamond requires ≥2 stable bands
//   6. Slanting bands tracked separately from stable bands
//   7. Multiple bands splitting → multiple diamonds
//   8. summarizeDiamonds counts loose / strict / strict2 correctly

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[dd] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[dd] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['detectDiamonds', 'summarizeDiamonds',
                    '_ddComputeBandStats', '_ddDetectSplittingRange',
                    '_ddIsBandStable']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed on window');
}
ok('all 5 diamond-detector helpers exposed on window');

// Helper: build synthetic state.data with N windows, each containing per-sample
// PC1 values according to a generator function `pc1AtWindow(wi, si)`.
function buildSyntheticData(nSamples, nWindows, pc1AtWindow) {
  const wins = [];
  for (let wi = 0; wi < nWindows; wi++) {
    const pc1 = new Array(nSamples);
    for (let si = 0; si < nSamples; si++) pc1[si] = pc1AtWindow(wi, si);
    wins.push({ pc1, pc2: new Array(nSamples).fill(0) });
  }
  return { n_samples: nSamples, windows: wins };
}

// ---- 2. No diamond on flat 3-band data ----
// 30 samples in 3 flat bands at PC1=-2, 0, +2 across 20 windows
const nS = 30;
const flatLabels = new Int8Array(nS);
for (let si = 0; si < 10; si++) flatLabels[si] = 0;
for (let si = 10; si < 20; si++) flatLabels[si] = 1;
for (let si = 20; si < 30; si++) flatLabels[si] = 2;

win.state.data = buildSyntheticData(nS, 20, (wi, si) => {
  const band = flatLabels[si];
  const center = [-2, 0, 2][band];
  return center + (((si * 7919 + wi * 31) % 100) / 100 - 0.5) * 0.2;   // small jitter
});

const candFlat = {
  K: 3,
  locked_labels: flatLabels,
  start_w: 0,
  end_w: 19,
};
const flatDiamonds = win.detectDiamonds(candFlat);
if (flatDiamonds.length !== 0) fail('flat 3-band data should produce no diamonds, got ' + flatDiamonds.length);
ok('flat 3-band data: no diamonds detected');

const flatSummary = win.summarizeDiamonds(candFlat);
if (flatSummary.has_loose) fail('flat data summary should not have loose diamond');
if (flatSummary.n_diamonds !== 0) fail('flat data n_diamonds should be 0');
ok('summarizeDiamonds (flat): n_diamonds=0, has_loose=false');

// ---- 3. Loose diamond detected when one band splits ----
// Band 1 (middle) splits into two sub-bands at windows 6-13.
// Bands 0 and 2 stay flat throughout.
const splitLabels = new Int8Array(nS);
for (let si = 0; si < 10; si++) splitLabels[si] = 0;
for (let si = 10; si < 20; si++) splitLabels[si] = 1;
for (let si = 20; si < 30; si++) splitLabels[si] = 2;

win.state.data = buildSyntheticData(nS, 20, (wi, si) => {
  const band = splitLabels[si];
  if (band === 0) return -2 + (((si * 31) % 100) / 100 - 0.5) * 0.1;   // flat
  if (band === 2) return  2 + (((si * 31) % 100) / 100 - 0.5) * 0.1;   // flat
  // Band 1: in windows 6-13, split into ±0.6 from center
  if (band === 1) {
    const subBand = (si - 10) < 5 ? -1 : 1;   // first 5 sub-band low, last 5 high
    if (wi >= 6 && wi <= 13) return 0 + subBand * 0.6;
    return 0 + (((si * 31) % 100) / 100 - 0.5) * 0.05;   // flat outside diamond
  }
  return 0;
});

const candSplit = {
  K: 3,
  locked_labels: splitLabels,
  start_w: 0,
  end_w: 19,
};
const splitDiamonds = win.detectDiamonds(candSplit);
if (splitDiamonds.length === 0) fail('splitting band should produce at least one diamond');
const d0 = splitDiamonds[0];
if (d0.splitting_band !== 1) fail('splitting band should be band 1, got ' + d0.splitting_band);
if (d0.diamond_start_w < 5 || d0.diamond_start_w > 7) fail('diamond start should be ~window 6, got ' + d0.diamond_start_w);
if (d0.diamond_end_w < 12 || d0.diamond_end_w > 14) fail('diamond end should be ~window 13, got ' + d0.diamond_end_w);
ok('loose diamond detected: splitting_band=1, windows ' + d0.diamond_start_w + '-' + d0.diamond_end_w);

// ---- 4. Strict diamond: ≥1 stable band ----
// Bands 0 and 2 are flat, so both should be stable
if (!d0.strict) fail('diamond should be strict (bands 0 and 2 are flat)');
if (d0.stable_bands.length < 1) fail('strict diamond needs ≥1 stable band, got ' + d0.stable_bands.length);
ok('strict diamond: ' + d0.stable_bands.length + ' stable bands [' + d0.stable_bands.join(', ') + ']');

// ---- 5. Strict2 diamond: ≥2 stable bands ----
if (!d0.strict2) fail('with 2 flat bands, diamond should be strict2');
ok('strict2 diamond: ≥2 stable bands present');

// ---- 6. Slanting bands tracked separately ----
// Build a case where band 0 slants and band 2 stays flat
win.state.data = buildSyntheticData(nS, 20, (wi, si) => {
  const band = splitLabels[si];
  if (band === 0) {
    // Slant: band 0 starts at -2.5 and rises to -1.5 across the diamond windows
    const slope = (wi >= 6 && wi <= 13) ? (wi - 6) * 0.1 : (wi < 6 ? 0 : 0.7);
    return -2.5 + slope;
  }
  if (band === 2) return 2 + (((si * 31) % 100) / 100 - 0.5) * 0.05;  // flat
  // Band 1: same diamond as before
  if (band === 1) {
    const subBand = (si - 10) < 5 ? -1 : 1;
    if (wi >= 6 && wi <= 13) return 0 + subBand * 0.6;
    return 0 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
  }
  return 0;
});

const slantDiamonds = win.detectDiamonds(candSplit);
if (slantDiamonds.length === 0) fail('slanting case should still produce a diamond');
const ds = slantDiamonds[0];
// Band 0 should now be slanting (drift > 5%), band 2 still stable
if (!ds.slanting_bands.includes(0)) fail('band 0 should be in slanting_bands, got ' + JSON.stringify(ds.slanting_bands));
if (!ds.stable_bands.includes(2)) fail('band 2 should be stable, got ' + JSON.stringify(ds.stable_bands));
if (ds.strict2) fail('with only 1 stable band, should not be strict2');
ok('slanting bands tracked: stable=[' + ds.stable_bands.join(',') + '], slanting=[' + ds.slanting_bands.join(',') + ']');

// ---- 7. Multiple bands splitting → multiple diamonds ----
// Build a case where band 0 AND band 1 both split (different positions).
// 4 bands total: 0 splits at windows 4-7, 1 splits at windows 12-16.
const nS4 = 40;
const fourLabels = new Int8Array(nS4);
for (let si = 0; si < 10; si++) fourLabels[si] = 0;
for (let si = 10; si < 20; si++) fourLabels[si] = 1;
for (let si = 20; si < 30; si++) fourLabels[si] = 2;
for (let si = 30; si < 40; si++) fourLabels[si] = 3;

win.state.data = buildSyntheticData(nS4, 20, (wi, si) => {
  const band = fourLabels[si];
  if (band === 0) {
    const subBand = (si - 0) < 5 ? -1 : 1;
    if (wi >= 4 && wi <= 7) return -3 + subBand * 0.6;
    return -3 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
  }
  if (band === 1) {
    const subBand = (si - 10) < 5 ? -1 : 1;
    if (wi >= 12 && wi <= 16) return -1 + subBand * 0.6;
    return -1 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
  }
  if (band === 2) return  1 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
  if (band === 3) return  3 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
  return 0;
});

const candFour = {
  K: 4,
  locked_labels: fourLabels,
  start_w: 0,
  end_w: 19,
};
const multiDiamonds = win.detectDiamonds(candFour);
if (multiDiamonds.length !== 2) fail('two splitting bands should produce 2 diamonds, got ' + multiDiamonds.length);
const d_b0 = multiDiamonds.find(d => d.splitting_band === 0);
const d_b1 = multiDiamonds.find(d => d.splitting_band === 1);
if (!d_b0 || !d_b1) fail('expected diamonds for bands 0 and 1');
if (d_b0.diamond_start_w > 5) fail('band 0 diamond should start near window 4');
if (d_b1.diamond_start_w < 11) fail('band 1 diamond should start near window 12');
ok('multiple diamonds: band 0 (windows ' + d_b0.diamond_start_w + '-' + d_b0.diamond_end_w +
   '), band 1 (windows ' + d_b1.diamond_start_w + '-' + d_b1.diamond_end_w + ')');

// ---- 8. summarizeDiamonds ----
const sumMulti = win.summarizeDiamonds(candFour);
if (sumMulti.n_diamonds !== 2) fail('n_diamonds should be 2');
if (!sumMulti.has_loose) fail('has_loose should be true');
if (!sumMulti.has_strict) fail('has_strict should be true (bands 2 and 3 are flat)');
if (!sumMulti.has_strict2) fail('has_strict2 should be true (≥2 stable bands)');
if (sumMulti.n_strict2 !== 2) fail('n_strict2 should be 2 (both diamonds have ≥2 stable bands)');
ok('summarizeDiamonds (multi): n=2, n_strict=2, n_strict2=2');

// ---- 9. Edge case: candidate without locked_labels ----
const empty = win.detectDiamonds({ K: 3 });
if (empty.length !== 0) fail('candidate without labels should give empty diamond array');
ok('candidate without locked_labels: empty diamond array');

// ---- 10. Edge case: out-of-range window indices ----
const badRange = win.detectDiamonds({ K: 3, locked_labels: flatLabels, start_w: -5, end_w: -1 });
// start_w=-5 still gets passed to _ddDetectSplittingRange but loop just yields no windows
if (badRange.length !== 0) fail('out-of-range window indices should give empty diamond array');
ok('out-of-range window indices: graceful empty result');

// ---- 11. Edge case: candidate without start_w/end_w ----
const noRange = win.detectDiamonds({ K: 3, locked_labels: flatLabels });
if (noRange.length !== 0) fail('no window range should give empty diamond array');
ok('no start_w/end_w: empty diamond array');

console.log('\n[dd] ALL CHECKS PASSED');
