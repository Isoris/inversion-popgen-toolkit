// Tests for v4 turn 90 — Detailed-mode band derivation classifier.
// Coverage:
//   1. Helpers exposed
//   2. Mode detection: unimodal het distribution → 1 mode
//   3. Mode detection: clearly bimodal → 2 modes
//   4. Mode detection: small-N edge cases → 1 mode (insufficient samples)
//   5. Divergence class thresholds
//   6. Band interpretation: hom-like, het-like, mixed, ambiguous
//   7. Haplotype assignment for 3-band biallelic system → H1/H1, H1/H2, H2/H2
//   8. Haplotype assignment for 6-band 3-haplotype system
//   9. Mixed band sub-band propagation
//  10. Atlas-proxy fallback when no per_sample_het is provided
//  11. Source tagging: 'precomp_het' vs 'atlas_proxy'
//  12. Edge case: empty band

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[dbd] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[dbd] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['classifyDetailedCandidate', '_dbdDetectModes',
                    '_dbdDivergenceClass', '_dbdBandInterpretation',
                    '_dbdAssignHaplotypes', '_dbdHetProxyFromPCResidual']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed on window');
}
if (typeof win._DBD_LOW_HET_THRESHOLD !== 'number') fail('_DBD_LOW_HET_THRESHOLD missing');
if (typeof win._DBD_HIGH_HET_THRESHOLD !== 'number') fail('_DBD_HIGH_HET_THRESHOLD missing');
ok('all 6 classifier helpers + 2 thresholds exposed on window');

// ---- 2. Unimodal mode detection ----
const unimodal = new Float32Array([0.2, 0.21, 0.22, 0.19, 0.20, 0.21, 0.18, 0.22, 0.20, 0.19]);
const m1 = win._dbdDetectModes(unimodal);
if (m1.n_modes !== 1) fail('unimodal het distribution should give 1 mode, got ' + m1.n_modes);
if (m1.mode_centers.length !== 1) fail('unimodal should give 1 center');
ok('mode detection: unimodal het distribution → 1 mode (center ' + m1.mode_centers[0].toFixed(3) + ')');

// ---- 3. Bimodal mode detection ----
// Two clear clusters: low-het ~0.2 and high-het ~0.6
const bimodal = new Float32Array([
  0.18, 0.20, 0.22, 0.19, 0.21, 0.20,    // low cluster
  0.58, 0.60, 0.62, 0.59, 0.61, 0.60,    // high cluster
]);
const m2 = win._dbdDetectModes(bimodal);
if (m2.n_modes !== 2) fail('bimodal het distribution should give 2 modes, got ' + m2.n_modes);
if (m2.mode_centers.length !== 2) fail('bimodal should give 2 centers');
if (m2.mode_centers[0] >= m2.mode_centers[1]) fail('first center should be lower');
const lowC = m2.mode_centers[0];
const highC = m2.mode_centers[1];
if (Math.abs(lowC - 0.20) > 0.05) fail('low center should be ~0.20, got ' + lowC.toFixed(3));
if (Math.abs(highC - 0.60) > 0.05) fail('high center should be ~0.60, got ' + highC.toFixed(3));
ok('mode detection: bimodal (0.2 vs 0.6) → 2 modes, centers at ' + lowC.toFixed(2) + ' / ' + highC.toFixed(2));

// ---- 4. Small-N edge case ----
const tooSmall = new Float32Array([0.2, 0.6, 0.3]);
const mSmall = win._dbdDetectModes(tooSmall);
if (mSmall.n_modes !== 1) fail('N=3 should yield 1 mode (below MIN_MODE_SAMPLES per cluster)');
ok('mode detection: small-N (N=3) → 1 mode (insufficient samples per cluster)');

// Empty input
const mEmpty = win._dbdDetectModes(new Float32Array(0));
if (mEmpty.n_modes !== 1) fail('empty input should give 1 mode');
ok('mode detection: empty input → 1 mode (graceful fallback)');

// ---- 5. Divergence class thresholds ----
if (win._dbdDivergenceClass(0.1) !== 'low') fail('het=0.1 should be low');
if (win._dbdDivergenceClass(0.45) !== 'mid') fail('het=0.45 should be mid');
if (win._dbdDivergenceClass(0.9) !== 'high') fail('het=0.9 should be high');
ok('divergence class: 0.1→low, 0.45→mid, 0.9→high');

// Boundary cases
if (win._dbdDivergenceClass(0.39) !== 'low') fail('het=0.39 should still be low (just below 0.40)');
if (win._dbdDivergenceClass(0.41) !== 'mid') fail('het=0.41 should be mid');
if (win._dbdDivergenceClass(0.49) !== 'mid') fail('het=0.49 should still be mid');
if (win._dbdDivergenceClass(0.50) !== 'high') fail('het=0.50 should be high (>=)');
ok('divergence boundary: thresholds at 0.40 and 0.50 (low<0.40, mid in [0.40, 0.50), high>=0.50)');

// ---- 6. Band interpretation ----
if (win._dbdBandInterpretation(0.85, 0.10, 0.05, 1) !== 'hom-like') fail('80%+ low → hom-like');
if (win._dbdBandInterpretation(0.05, 0.10, 0.85, 1) !== 'het-like') fail('80%+ high → het-like');
if (win._dbdBandInterpretation(0.40, 0.40, 0.20, 2) !== 'mixed') fail('bimodal → mixed');
if (win._dbdBandInterpretation(0.10, 0.70, 0.20, 1) !== 'ambiguous') fail('mostly mid → ambiguous');
ok('band interpretation: hom-like (low>=0.7) / het-like (high>=0.7) / mixed (multimodal) / ambiguous (mostly mid)');

// ---- 7. 3-band biallelic system → H1/H1, H1/H2, H2/H2 ----
// Construct a simple candidate: 3 bands, hom-like-low + het-like-mid + hom-like-high
// All hom-like + het-like classifications driven by per-sample het
const candidate3 = {
  K: 3,
  locked_labels: new Int8Array([0,0,0,0,0,  1,1,1,1,1,  2,2,2,2,2]),  // 5 in each band
  start_w: 10,
  end_w: 20,
  // Provide samplePC1 directly so classifier doesn't need to read state.data
  samplePC1: [
    -2,-2.1,-1.9,-2.0,-2.1,    // band 0 (low PC1)
    0,0.1,-0.1,0.05,-0.05,     // band 1 (mid PC1)
    2,1.9,2.1,2.0,2.05,        // band 2 (high PC1)
  ],
};
// Per-sample het: bands 0 and 2 = low (0.2), band 1 = high (0.6)
const het3 = new Float32Array([
  0.20,0.18,0.22,0.21,0.19,    // band 0: low het
  0.60,0.62,0.58,0.61,0.59,    // band 1: high het
  0.20,0.21,0.19,0.22,0.18,    // band 2: low het
]);
const result3 = win.classifyDetailedCandidate(candidate3, het3);
if (!result3) fail('classify returned null');
if (result3.n_samples !== 15) fail('n_samples should be 15');
if (result3.n_bands !== 3) fail('n_bands should be 3');
if (result3.source !== 'precomp_het') fail('source should be precomp_het, got ' + result3.source);
ok('3-band biallelic: classified, source=precomp_het, n=15, K=3');

// Per-band interpretation: band 0 = hom-like, band 1 = het-like, band 2 = hom-like
const band0 = result3.per_band.find(b => b.band === 0);
const band1 = result3.per_band.find(b => b.band === 1);
const band2 = result3.per_band.find(b => b.band === 2);
if (band0.interpretation !== 'hom-like') fail('band 0 should be hom-like, got ' + band0.interpretation);
if (band1.interpretation !== 'het-like') fail('band 1 should be het-like, got ' + band1.interpretation);
if (band2.interpretation !== 'hom-like') fail('band 2 should be hom-like, got ' + band2.interpretation);
ok('3-band: band 0=hom-like, band 1=het-like, band 2=hom-like (correct)');

// Karyotype assignment: hom-like bands ordered by median PC1 → H1/H1 (low), H2/H2 (high)
// het-like band 1 is between them → H1/H2
const ka0 = result3.karyotype_assignment.find(a => a.band === 0);
const ka1 = result3.karyotype_assignment.find(a => a.band === 1);
const ka2 = result3.karyotype_assignment.find(a => a.band === 2);
if (ka0.haplotype_class !== 'H1/H1') fail('band 0 (low PC1) should be H1/H1, got ' + ka0.haplotype_class);
if (ka1.haplotype_class !== 'H1/H2') fail('band 1 (het, between) should be H1/H2, got ' + ka1.haplotype_class);
if (ka2.haplotype_class !== 'H2/H2') fail('band 2 (high PC1) should be H2/H2, got ' + ka2.haplotype_class);
ok('3-band haplotype assignment: H1/H1 (low PC1) / H1/H2 (het-like) / H2/H2 (high PC1)');

// ---- 8. 6-band 3-haplotype system ----
// 3 hom-like (low PC1, mid-low PC1, high PC1) + 3 het-like
const candidate6 = {
  K: 6,
  locked_labels: new Int8Array([
    0,0,0,0,    // band 0: hom-like, very low PC1 → H1/H1
    1,1,1,1,    // band 1: het-like, between H1 and H2 → H1/H2
    2,2,2,2,    // band 2: hom-like, mid PC1 → H2/H2
    3,3,3,3,    // band 3: het-like, between H1 and H3 → H1/H3
    4,4,4,4,    // band 4: het-like, between H2 and H3 → H2/H3
    5,5,5,5,    // band 5: hom-like, high PC1 → H3/H3
  ]),
  start_w: 0,
  end_w: 5,
  samplePC1: [
    -3,-3,-3,-3,           // band 0: lowest
    -1.5,-1.5,-1.5,-1.5,   // band 1: between band 0 and 2
    0,0,0,0,               // band 2: mid (will be H2/H2)
    1,1,1,1,               // band 3: between band 0 and band 5
    1.5,1.5,1.5,1.5,       // band 4: between band 2 and 5
    3,3,3,3,               // band 5: highest
  ],
};
// het: bands 0, 2, 5 = low; bands 1, 3, 4 = high
const het6 = new Float32Array([
  0.2,0.2,0.2,0.2,    // band 0
  0.6,0.6,0.6,0.6,    // band 1
  0.2,0.2,0.2,0.2,    // band 2
  0.7,0.7,0.7,0.7,    // band 3
  0.7,0.7,0.7,0.7,    // band 4
  0.2,0.2,0.2,0.2,    // band 5
]);
const result6 = win.classifyDetailedCandidate(candidate6, het6);
if (!result6) fail('K=6 classify returned null');
if (result6.n_bands !== 6) fail('K=6 n_bands should be 6');

// Find each band's assignment
const ka_b0 = result6.karyotype_assignment.find(a => a.band === 0);
const ka_b2 = result6.karyotype_assignment.find(a => a.band === 2);
const ka_b5 = result6.karyotype_assignment.find(a => a.band === 5);
if (ka_b0.haplotype_class !== 'H1/H1') fail('band 0 (lowest PC1) should be H1/H1, got ' + ka_b0.haplotype_class);
if (ka_b2.haplotype_class !== 'H2/H2') fail('band 2 (mid PC1) should be H2/H2, got ' + ka_b2.haplotype_class);
if (ka_b5.haplotype_class !== 'H3/H3') fail('band 5 (highest PC1) should be H3/H3, got ' + ka_b5.haplotype_class);
ok('K=6 haplotype assignment: H1/H1 (lowest), H2/H2 (mid), H3/H3 (highest)');

// Het-like bands should have intermediate-pair labels
const ka_b1 = result6.karyotype_assignment.find(a => a.band === 1);
const ka_b3 = result6.karyotype_assignment.find(a => a.band === 3);
const ka_b4 = result6.karyotype_assignment.find(a => a.band === 4);
// Each is between two homs; the closest pair determines its label
if (!ka_b1.haplotype_class || !ka_b1.haplotype_class.startsWith('H1/H')) fail('band 1 should be a het-like H1/H? class, got ' + ka_b1.haplotype_class);
if (!ka_b3.haplotype_class || !ka_b3.haplotype_class.includes('H')) fail('band 3 should have an H-class, got ' + ka_b3.haplotype_class);
ok('K=6 het-like bands: band 1 → ' + ka_b1.haplotype_class + ', band 3 → ' + ka_b3.haplotype_class + ', band 4 → ' + ka_b4.haplotype_class);

// ---- 9. Mixed band sub-band propagation ----
const candidateMixed = {
  K: 2,
  locked_labels: new Int8Array([0,0,0,0,0,0,0,0,0,0,  1,1,1,1,1,1,1,1,1,1]),
  samplePC1: new Array(20).fill(0),
  start_w: 0, end_w: 0,
};
// Band 0: bimodal het — half low, half high
// Band 1: unimodal low
const hetMixed = new Float32Array([
  0.2,0.2,0.2,0.2,0.22,    // band 0 low half
  0.6,0.6,0.6,0.6,0.62,    // band 0 high half
  0.2,0.21,0.2,0.2,0.19,   // band 1 unimodal low
  0.21,0.2,0.2,0.21,0.2,   // band 1 unimodal low
]);
const resultMixed = win.classifyDetailedCandidate(candidateMixed, hetMixed);
const bandMixed = resultMixed.per_band.find(b => b.band === 0);
if (!bandMixed.is_mixed) fail('band 0 should be flagged is_mixed');
if (bandMixed.n_modes !== 2) fail('band 0 should have n_modes=2, got ' + bandMixed.n_modes);
if (bandMixed.interpretation !== 'mixed') fail('band 0 interpretation should be "mixed"');
ok('mixed band: is_mixed=true, n_modes=2, interpretation="mixed"');

// Sub-band propagated to per_sample
const sample0 = resultMixed.per_sample[0];   // first low fish in band 0
const sample5 = resultMixed.per_sample[5];   // first high fish in band 0
const sample10 = resultMixed.per_sample[10]; // first band-1 fish
if (sample0.sub_band == null) fail('mixed-band fish should have sub_band assigned');
if (sample5.sub_band == null) fail('mixed-band fish should have sub_band assigned');
if (sample0.sub_band === sample5.sub_band) fail('low and high fish in mixed band should have different sub_bands');
if (sample10.sub_band !== null) fail('unimixed-band fish should have sub_band=null, got ' + sample10.sub_band);
ok('sub-band propagation: low+high fish in mixed band get distinct sub_bands; unimodal-band fish get null');

// ---- 10. Atlas-proxy fallback ----
// Without per_sample_het, classifier should fall back to PC-residual proxy
win.state.data = {
  n_samples: 15,
  windows: candidate3.samplePC1 ? (() => {
    // Construct windows array with pc1/pc2 so _diagComputeCandidateSuspicion works
    const wins = [];
    for (let wi = 10; wi <= 20; wi++) {
      wins.push({ pc1: candidate3.samplePC1.slice(), pc2: new Array(15).fill(0) });
    }
    return wins;
  })() : []
};
// Pad windows array up to index 20
while (win.state.data.windows.length < 21) win.state.data.windows.unshift({});
win.state.k = 3;
if (typeof win._diagClearCache === 'function') win._diagClearCache();

const proxy = win._dbdHetProxyFromPCResidual(candidate3, 15);
if (!proxy) fail('PC-residual proxy returned null');
if (proxy.length !== 15) fail('proxy length should be 15');
ok('atlas-proxy fallback: returns Float32Array length 15');

// Use the proxy in classifier (no per_sample_het)
const resultProxy = win.classifyDetailedCandidate(candidate3, null);
if (!resultProxy) fail('classifier with no het should still return result via proxy');
if (resultProxy.source !== 'atlas_proxy') fail('source should be atlas_proxy, got ' + resultProxy.source);
ok('classifier: no per_sample_het → falls back to atlas_proxy (source tagged)');

// ---- 11. Source tagging ----
const r1 = win.classifyDetailedCandidate(candidate3, het3);
const r2 = win.classifyDetailedCandidate(candidate3, null);
if (r1.source !== 'precomp_het') fail('with het → precomp_het');
if (r2.source !== 'atlas_proxy') fail('without het → atlas_proxy');
ok('source tagging correct: precomp_het when het provided, atlas_proxy when not');

// ---- 12. Edge cases ----
// Empty band — band index that has no samples
const candidateSparse = {
  K: 4,
  locked_labels: new Int8Array([0,0,2,2]),    // band 1 and band 3 are empty
  samplePC1: [-1, -1, 1, 1],
  start_w: 10, end_w: 20,
};
const hetSparse = new Float32Array([0.2, 0.2, 0.2, 0.2]);
const resultSparse = win.classifyDetailedCandidate(candidateSparse, hetSparse);
const emptyBand = resultSparse.per_band.find(b => b.band === 1);
if (emptyBand.n !== 0) fail('empty band should have n=0');
if (emptyBand.interpretation !== 'ambiguous') fail('empty band interpretation should be ambiguous');
ok('empty band: n=0, interpretation=ambiguous (no crash)');

// Null candidate
if (win.classifyDetailedCandidate(null, null) !== null) fail('null candidate should give null');
ok('null candidate: returns null gracefully');

// Candidate without locked_labels
const noLabels = { K: 3 };
if (win.classifyDetailedCandidate(noLabels, null) !== null) fail('no locked_labels should give null');
ok('candidate without locked_labels: returns null gracefully');

console.log('\n[dbd] ALL CHECKS PASSED');
