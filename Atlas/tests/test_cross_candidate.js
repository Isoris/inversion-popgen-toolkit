// Tests for cross-candidate inheritance compute (turn 2a).
// We don't load full atlas state; we exercise the pure compute functions
// with controlled label arrays.

const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');

const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
});
const w = dom.window;

function run() {
  const failures = [];
  let testNum = 0;
  function t(name, fn) {
    testNum++;
    try {
      fn();
      console.log(`  PASS [${testNum}] ${name}`);
    } catch (e) {
      failures.push({ name, err: e.message });
      console.log(`  FAIL [${testNum}] ${name}: ${e.message}`);
    }
  }
  function eq(actual, expected, msg) {
    if (actual !== expected) {
      throw new Error(`${msg || ''} expected ${JSON.stringify(expected)}, got ${JSON.stringify(actual)}`);
    }
  }
  function approx(actual, expected, tol, msg) {
    if (!isFinite(actual) || Math.abs(actual - expected) > tol) {
      throw new Error(`${msg || ''} expected ~${expected} (tol ${tol}), got ${actual}`);
    }
  }

  // Sanity: functions exposed
  if (typeof w.crossCandidateContingency !== 'function'
      || typeof w.crossCandidateMatrix !== 'function'
      || typeof w._bandFingerprintFor !== 'function'
      || typeof w._chiSqSurvival !== 'function') {
    failures.push({ name: 'functions exposed', err: 'one or more functions missing' });
    console.log('  FAIL: functions not exposed on window');
    return failures;
  }
  console.log('  All compute functions exposed on window');

  // ============================================================
  // _lnGamma sanity
  // ============================================================
  t('lnGamma(1) = 0', () => approx(w._lnGamma(1), 0, 1e-9));
  t('lnGamma(2) = 0', () => approx(w._lnGamma(2), 0, 1e-9));
  t('lnGamma(5) = ln(24) ~ 3.1781', () => approx(w._lnGamma(5), Math.log(24), 1e-6));
  t('lnGamma(10) = ln(362880) ~ 12.8018', () => approx(w._lnGamma(10), Math.log(362880), 1e-6));

  // ============================================================
  // _chiSqSurvival sanity (compare to known critical values)
  // ============================================================
  // df=1: chi2=3.841 -> p=0.05; chi2=6.635 -> p=0.01
  t('chi2 survival df=1, x=3.841 ~ 0.05', () => approx(w._chiSqSurvival(3.841, 1), 0.05, 0.005));
  t('chi2 survival df=1, x=6.635 ~ 0.01', () => approx(w._chiSqSurvival(6.635, 1), 0.01, 0.002));
  // df=4: chi2=9.488 -> p=0.05
  t('chi2 survival df=4, x=9.488 ~ 0.05', () => approx(w._chiSqSurvival(9.488, 4), 0.05, 0.005));
  // df=25: chi2=37.65 -> p=0.05
  t('chi2 survival df=25, x=37.65 ~ 0.05', () => approx(w._chiSqSurvival(37.65, 25), 0.05, 0.005));
  // Boundaries
  t('chi2 survival x=0 returns 1', () => eq(w._chiSqSurvival(0, 5), 1));
  t('chi2 survival df<=0 returns NaN', () => {
    const r = w._chiSqSurvival(5, 0);
    if (!isNaN(r)) throw new Error(`expected NaN, got ${r}`);
  });

  // ============================================================
  // _bandFingerprintFor: one band's distribution at target
  // ============================================================
  t('bandFingerprintFor: perfect linkage (all band 0 of A go to band 0 of B)', () => {
    // 50 fish: 25 in (A=0, B=0), 25 in (A=1, B=1)
    const labA = [], labB = [];
    for (let i = 0; i < 25; i++) { labA.push(0); labB.push(0); }
    for (let i = 0; i < 25; i++) { labA.push(1); labB.push(1); }
    const fp = w._bandFingerprintFor(labA, labB, 2, 2, 0, 'A', 'B');
    eq(fp.n, 25);
    eq(fp.mode_band, 0);
    approx(fp.mode_share, 1.0, 1e-9);
    approx(fp.entropy, 0, 1e-9);
    approx(fp.concentration, 1, 1e-9);
  });

  t('bandFingerprintFor: uniform distribution -> concentration ~ 0', () => {
    // 60 fish in band 0 of A spread evenly across 6 bands of B
    const labA = [], labB = [];
    for (let b = 0; b < 6; b++) {
      for (let i = 0; i < 10; i++) { labA.push(0); labB.push(b); }
    }
    const fp = w._bandFingerprintFor(labA, labB, 1, 6, 0, 'A', 'B');
    eq(fp.n, 60);
    approx(fp.concentration, 0, 1e-6);
    approx(fp.entropy, Math.log(6), 1e-6);
    approx(fp.mode_share, 1/6, 1e-6);
  });

  t('bandFingerprintFor: n=0 (no fish in source band) returns concentration NaN', () => {
    const labA = [0, 0, 0], labB = [0, 1, 2];
    const fp = w._bandFingerprintFor(labA, labB, 2, 3, 1, 'A', 'B');  // band 1 of A is empty
    eq(fp.n, 0);
    if (!isNaN(fp.concentration)) throw new Error(`expected NaN concentration, got ${fp.concentration}`);
  });

  t('bandFingerprintFor: invalid band index returns null', () => {
    const fp = w._bandFingerprintFor([0, 1], [0, 1], 2, 2, 5, 'A', 'B');
    if (fp !== null) throw new Error('expected null for out-of-range band');
  });

  // ============================================================
  // crossCandidateContingency: pairwise full output
  // ============================================================
  t('CCC: perfect 2x2 linkage (Cramér V=1, p~0)', () => {
    const labA = [], labB = [];
    for (let i = 0; i < 50; i++) { labA.push(0); labB.push(0); }
    for (let i = 0; i < 50; i++) { labA.push(1); labB.push(1); }
    const r = w.crossCandidateContingency(labA, labB, 2, 2);
    approx(r.cramer_v, 1.0, 1e-9);
    eq(r.df, 1);
    if (r.p_value > 1e-10) throw new Error(`expected p~0, got ${r.p_value}`);
    eq(r.n, 100);
  });

  t('CCC: independence (Cramér V ~ 0, p large)', () => {
    // 200 fish, A and B labels independent uniform random over 2 bands.
    // Build deterministically by interleaving so independence is exact.
    const labA = [], labB = [];
    for (let a = 0; a < 2; a++) {
      for (let b = 0; b < 2; b++) {
        for (let i = 0; i < 50; i++) { labA.push(a); labB.push(b); }
      }
    }
    // Now table is exactly 50/50/50/50. Cramer V = 0.
    const r = w.crossCandidateContingency(labA, labB, 2, 2);
    approx(r.cramer_v, 0, 1e-9);
    approx(r.chi2, 0, 1e-9);
    if (r.p_value < 0.99) throw new Error(`expected p~1 for exact independence, got ${r.p_value}`);
  });

  t('CCC: partial linkage (Cramér V intermediate)', () => {
    // 200 fish: 80 in (0,0), 20 in (0,1), 20 in (1,0), 80 in (1,1)
    // Strong but imperfect linkage.
    const labA = [], labB = [];
    for (let i = 0; i < 80; i++) { labA.push(0); labB.push(0); }
    for (let i = 0; i < 20; i++) { labA.push(0); labB.push(1); }
    for (let i = 0; i < 20; i++) { labA.push(1); labB.push(0); }
    for (let i = 0; i < 80; i++) { labA.push(1); labB.push(1); }
    const r = w.crossCandidateContingency(labA, labB, 2, 2);
    // chi2 = sum (O-E)^2/E ; E_ij = 100*100/200 = 50 for all cells
    // chi2 = 4 * (80-50)^2/50 mixed... actually 2*(80-50)^2/50 + 2*(20-50)^2/50
    //       = 2*30^2/50 + 2*30^2/50 = 4*900/50 = 72
    // V = sqrt(72/(200*1)) = sqrt(0.36) = 0.6
    approx(r.chi2, 72, 1e-6);
    approx(r.cramer_v, 0.6, 1e-6);
    eq(r.df, 1);
    if (r.p_value > 1e-15) throw new Error(`expected p very small, got ${r.p_value}`);
  });

  t('CCC: 3x3 table with K=3 candidates', () => {
    // Symmetric strong linkage at K=3
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
      for (let i = 0; i < 5; i++)  { labA.push(k); labB.push((k+1)%3); }
    }
    const r = w.crossCandidateContingency(labA, labB, 3, 3);
    eq(r.KA, 3); eq(r.KB, 3); eq(r.n, 105);
    eq(r.df, 4);
    if (r.cramer_v < 0.7) throw new Error(`expected high V, got ${r.cramer_v}`);
    if (r.cramer_v > 1) throw new Error(`V > 1 invalid: ${r.cramer_v}`);
  });

  t('CCC: 6x6 K=6 multi-haplotype scenario', () => {
    // 6 H-classes, 30 fish each, all on diagonal (perfect linkage)
    const labA = [], labB = [];
    for (let k = 0; k < 6; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    }
    const r = w.crossCandidateContingency(labA, labB, 6, 6);
    eq(r.n, 180);
    approx(r.cramer_v, 1.0, 1e-9);
    eq(r.dominant_patterns.length, 6);
    // Each fingerprint should be perfectly concentrated
    for (let i = 0; i < 6; i++) {
      approx(r.band_fingerprints[i].concentration, 1, 1e-9);
      eq(r.band_fingerprints[i].mode_band, i);
    }
  });

  t('CCC: handles -1 labels (invalid samples) by excluding them', () => {
    const labA = [0, 0, 1, 1, -1, -1, 0, 1];
    const labB = [0, 0, 1, 1, 0, 1, -1, -1];
    const r = w.crossCandidateContingency(labA, labB, 2, 2, { min_n: 4 });
    eq(r.n, 4);  // only 4 fish have both labels valid
  });

  t('CCC: low-n returns warning result', () => {
    const labA = [0, 1, 0];
    const labB = [0, 1, 0];
    const r = w.crossCandidateContingency(labA, labB, 2, 2, { min_n: 10 });
    eq(r.n, 3);
    if (!r.low_count_warning) throw new Error('expected low_count_warning=true');
    if (r.dominant_patterns.length > 0) throw new Error('expected no dominant_patterns at low n');
  });

  t('CCC: shape mismatch returns null', () => {
    const r = w.crossCandidateContingency([0, 1], [0, 1, 0], 2, 2);
    if (r !== null) throw new Error(`expected null, got ${JSON.stringify(r)}`);
  });

  t('CCC: dominant_patterns sorted by count desc, capped at 6', () => {
    // 7 distinct hot cells of varying sizes
    const labA = [], labB = [];
    const cells = [[0,0,40], [0,1,35], [1,0,30], [1,1,25], [2,0,20], [2,1,15], [3,0,10]];
    for (const [a, b, n] of cells) {
      for (let i = 0; i < n; i++) { labA.push(a); labB.push(b); }
    }
    const r = w.crossCandidateContingency(labA, labB, 4, 2);
    eq(r.dominant_patterns.length, 6);
    eq(r.dominant_patterns[0].count, 40);
    eq(r.dominant_patterns[1].count, 35);
    eq(r.dominant_patterns[5].count, 15);
    // The 7th cell (count 10) should be excluded
    for (const p of r.dominant_patterns) {
      if (p.count < 15) throw new Error('dominant_patterns has cell below threshold');
    }
  });

  t('CCC: low_count_warning on sparse table', () => {
    // 3x3 table with most cells ~1
    const labA = [], labB = [];
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        labA.push(a); labB.push(b);  // 1 each
      }
    }
    // Beef up one cell to clear min_n
    for (let i = 0; i < 20; i++) { labA.push(0); labB.push(0); }
    const r = w.crossCandidateContingency(labA, labB, 3, 3);
    if (!r.low_count_warning) throw new Error('expected low_count_warning on sparse 3x3');
  });

  t('CCC: mismatched K (KA=3, KB=6) works correctly', () => {
    // Each band of A maps to 2 bands of B (1:2 mapping)
    const labA = [], labB = [];
    for (let a = 0; a < 3; a++) {
      for (let i = 0; i < 25; i++) { labA.push(a); labB.push(2*a); }
      for (let i = 0; i < 25; i++) { labA.push(a); labB.push(2*a + 1); }
    }
    const r = w.crossCandidateContingency(labA, labB, 3, 6);
    eq(r.KA, 3); eq(r.KB, 6); eq(r.n, 150);
    eq(r.df, (3-1)*(6-1));  // 10
    // V = sqrt(chi2 / (n * min(KA-1,KB-1))) = sqrt(chi2 / (150*2))
    // Each row has 25+25 in two cells, 0 in others. Strong linkage.
    if (r.cramer_v < 0.9) throw new Error(`expected V > 0.9, got ${r.cramer_v}`);
  });

  t('CCC: band_fingerprints have correct length (one per row of A)', () => {
    const labA = [], labB = [];
    for (let i = 0; i < 100; i++) { labA.push(i % 4); labB.push(i % 6); }
    const r = w.crossCandidateContingency(labA, labB, 4, 6);
    eq(r.band_fingerprints.length, 4);  // one per band of A
    for (const fp of r.band_fingerprints) {
      eq(fp.dist.length, 6);  // KB-long
      eq(fp.target_K, 6);
    }
  });

  // ============================================================
  // crossCandidateMatrix
  // ============================================================
  t('CCM: 3 candidates, all linked', () => {
    // Build 3 candidates with perfectly correlated labels
    const labels = [];
    for (let i = 0; i < 90; i++) labels.push(i % 3);  // 0,1,2,0,1,2,...
    const items = [
      { id: 'A', labels, K: 3 },
      { id: 'B', labels, K: 3 },  // same labels as A
      { id: 'C', labels, K: 3 },
    ];
    const m = w.crossCandidateMatrix(items);
    eq(m.n_items, 3);
    eq(m.ids[0], 'A'); eq(m.ids[1], 'B'); eq(m.ids[2], 'C');
    // Diagonal = 1
    approx(m.cramer_v[0*3 + 0], 1, 1e-9);
    approx(m.cramer_v[1*3 + 1], 1, 1e-9);
    approx(m.cramer_v[2*3 + 2], 1, 1e-9);
    // Off-diagonal also = 1 (perfect linkage)
    approx(m.cramer_v[0*3 + 1], 1, 1e-9);
    approx(m.cramer_v[0*3 + 2], 1, 1e-9);
    approx(m.cramer_v[1*3 + 2], 1, 1e-9);
    // Symmetric
    approx(m.cramer_v[1*3 + 0], m.cramer_v[0*3 + 1], 1e-9);
  });

  t('CCM: 3 candidates, 2 linked + 1 independent', () => {
    const labA = [], labB = [], labC = [];
    // A and B perfectly correlated, C independent. To generate independence
    // exactly we tile the joint distribution P(A,C) = P(A) * P(C) uniformly:
    // for each (a, c) pair in {0,1,2} x {0,1,2}, write 22 fish.
    // That gives 198 fish, every joint cell equal -> Cramér V = 0.
    for (let a = 0; a < 3; a++) {
      for (let c = 0; c < 3; c++) {
        for (let i = 0; i < 22; i++) {
          labA.push(a);
          labB.push(a);  // B = A always
          labC.push(c);
        }
      }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
      { id: 'C', labels: labC, K: 3 },
    ];
    const m = w.crossCandidateMatrix(items);
    approx(m.cramer_v[0*3 + 1], 1, 1e-9, 'A-B should be perfect');
    approx(m.cramer_v[0*3 + 2], 0, 1e-9, 'A-C should be exactly 0');
    approx(m.cramer_v[1*3 + 2], 0, 1e-9, 'B-C should be exactly 0');
  });

  t('CCM: full_details=true populates pair_details cache', () => {
    const labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
    const items = [
      { id: 'X', labels, K: 2 },
      { id: 'Y', labels, K: 2 },
    ];
    const m = w.crossCandidateMatrix(items, { full_details: true, min_n: 5 });
    if (!m.pair_details.has('0:1')) throw new Error('pair_details cache empty');
    const det = m.pair_details.get('0:1');
    eq(det.KA, 2);
  });

  t('CCM: empty / invalid input returns null', () => {
    eq(w.crossCandidateMatrix([]), null);
    eq(w.crossCandidateMatrix(null), null);
  });

  t('CCM: handles items with mismatched K (K=3 vs K=6)', () => {
    const labA = [];
    const labB = [];
    for (let i = 0; i < 100; i++) {
      labA.push(i % 3);     // K=3
      labB.push(i % 6);     // K=6
    }
    const items = [
      { id: 'k3', labels: labA, K: 3 },
      { id: 'k6', labels: labB, K: 6 },
    ];
    const m = w.crossCandidateMatrix(items);
    eq(m.K_per_item[0], 3);
    eq(m.K_per_item[1], 6);
    if (!isFinite(m.cramer_v[0*2 + 1])) throw new Error('cross-K should yield finite V');
  });

  t('CCM: skips items missing labels gracefully', () => {
    const items = [
      { id: 'A', labels: [0, 1, 0, 1], K: 2 },
      { id: 'B', labels: null, K: 2 },
      { id: 'C', labels: [0, 0, 1, 1], K: 2 },
    ];
    const m = w.crossCandidateMatrix(items, { min_n: 2 });
    // diagonals still 1
    approx(m.cramer_v[0*3 + 0], 1, 1e-9);
    approx(m.cramer_v[1*3 + 1], 1, 1e-9);
    // pair involving B should be NaN (or 0, depending; we set to NaN)
    if (isFinite(m.cramer_v[0*3 + 1])) throw new Error('expected NaN for B pair');
  });

  // ============================================================
  // Realistic scenario: 226 fish, K=6, two linked candidates + one independent
  // ============================================================
  t('Realistic scenario: 226 fish, K=6, 2 linked candidates + 1 independent', () => {
    const N = 226;
    const labA = [], labB = [], labC = [];
    // Distribute 226 fish across 6 bands with realistic counts: 50/40/40/40/30/26
    const groups = [50, 40, 40, 40, 30, 26];
    let idx = 0;
    for (let g = 0; g < 6; g++) {
      for (let i = 0; i < groups[g]; i++) {
        labA.push(g);
        labB.push(g);  // perfect linkage A-B
        // C: independent — use modular hash to break correlation
        labC.push((g * 7 + i * 13) % 6);
        idx++;
      }
    }
    eq(idx, N);
    const items = [
      { id: 'cand_inv_LG28_A', labels: labA, K: 6 },
      { id: 'cand_inv_LG28_B', labels: labB, K: 6 },
      { id: 'cand_inv_LG12_C', labels: labC, K: 6 },
    ];
    const m = w.crossCandidateMatrix(items);
    // A-B perfect
    approx(m.cramer_v[0*3 + 1], 1, 1e-6);
    // A-C and B-C should be much smaller
    if (m.cramer_v[0*3 + 2] > 0.5) throw new Error(`A-C V=${m.cramer_v[0*3 + 2]} should be < 0.5`);
    if (m.cramer_v[1*3 + 2] > 0.5) throw new Error(`B-C V=${m.cramer_v[1*3 + 2]} should be < 0.5`);
  });

  return failures;
}

setTimeout(() => {
  const failures = run();
  if (failures.length > 0) {
    console.log(`\n${failures.length} test(s) failed`);
    failures.forEach(f => console.log(`  - ${f.name}: ${f.err}`));
    process.exit(1);
  } else {
    console.log('\nAll tests passed');
    process.exit(0);
  }
}, 500);
