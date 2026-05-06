// Test bandReachBimodality directly. We don't load full atlas state; we
// stub computeBandReachAcrossL2s to return synthetic per_sample_band_reach
// distributions and verify the verdict logic.

const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');

// Run atlas in JSDOM but bypass canvas-init crashes by stripping out
// inline scripts that reference DOM-canvas APIs at parse time.
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),  // silence
});

const w = dom.window;

// Wait for atlas script to evaluate. The function we need is exposed
// late but our insertion added it after computeBandReachAcrossL2s.
// JSDOM evaluates scripts synchronously in pretendToBeVisual mode.

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
    if (Math.abs(actual - expected) > tol) {
      throw new Error(`${msg || ''} expected ~${expected} (tol ${tol}), got ${actual}`);
    }
  }

  // Verify the function exists
  if (typeof w.bandReachBimodality !== 'function') {
    failures.push({ name: 'function exists', err: 'window.bandReachBimodality not defined' });
    console.log('  FAIL: window.bandReachBimodality not exposed');
    return failures;
  }
  console.log('  function bandReachBimodality is defined');

  // Stub computeBandReachAcrossL2s to return controlled inputs.
  // Each test case = a per-sample reach distribution.
  function makeStubReach(reachArray, K) {
    return function() {
      const Int8 = w.Int8Array || Int8Array;
      const Int32 = w.Int32Array || Int32Array;
      return {
        l2_indices: [0, 1, 2],
        n_samples: reachArray.length,
        n_valid: reachArray.filter(r => r > 0).length,
        K,
        per_sample_band_reach: new Int8(reachArray),
        per_band_visit_count: new Int32(K),
        narrow_fraction: 0,   // not used by the bimodality function
        bands_populated: K,
        regime_breadth: 'medium',
      };
    };
  }

  // Build a reach array with a given count per reach value
  // counts[r] = how many fish have reach == r. counts[0] = invalid fish.
  function buildReach(counts) {
    const arr = [];
    for (let r = 0; r < counts.length; r++) {
      for (let i = 0; i < counts[r]; i++) arr.push(r);
    }
    return arr;
  }

  // === Test cases ===

  t('K=3 returns NA verdict', () => {
    w.computeBandReachAcrossL2s = makeStubReach(buildReach([0, 100, 80, 40]), 3);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'NA');
  });

  t('Tiny n returns NA', () => {
    w.computeBandReachAcrossL2s = makeStubReach(buildReach([0, 3, 2, 1, 1]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'NA');
  });

  t('Pure narrow distribution -> UNIMODAL_NARROW', () => {
    // 90 narrow (reach 1 or 2), 5 mid, 5 wide
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 60, 30, 5, 3, 2]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'UNIMODAL_NARROW',
       `narrow=${r.narrow_fraction}, wide=${r.wide_fraction}, reason=${r.reason}`);
    approx(r.narrow_fraction, 0.90, 0.01);
  });

  t('Pure wide distribution -> UNIMODAL_WIDE', () => {
    // Most fish wandering, K=5
    // Use 10 narrow, 10 mid, 80 wide. wide_fraction should be 0.80.
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 5, 5, 10, 40, 40]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'UNIMODAL_WIDE',
       `narrow=${r.narrow_fraction}, wide=${r.wide_fraction}, reason=${r.reason}`);
  });

  t('Bimodal (Quentin scenario: half narrow, half wide) -> BIMODAL_STACKED', () => {
    // 113 fish narrow at reach=1-2, 113 fish wide at reach=4-5, only ~5 at reach=3
    // narrow=113, wide=108 (reach 4)+5(reach 5)=113... use 90/10 to distinguish
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 70, 40, 5, 60, 50]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'BIMODAL_STACKED',
       `narrow=${r.narrow_fraction}, wide=${r.wide_fraction}, trough=${r.trough_count}, reason=${r.reason}`);
    if (!r.trough_below_both) throw new Error('trough_below_both should be true');
  });

  t('Bimodal but trough not low enough -> UNDETERMINED', () => {
    // 60 narrow, 50 mid (NOT a trough), 60 wide.
    // narrow=0.35, wide=0.35, trough=0.30 — fails the dip test.
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 40, 20, 50, 40, 20]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'UNDETERMINED',
       `narrow=${r.narrow_fraction}, wide=${r.wide_fraction}, trough=${r.trough_count}, reason=${r.reason}`);
  });

  t('Equal narrow and wide with deep trough -> BIMODAL_STACKED', () => {
    // Symmetric: 100 narrow, 0 trough, 100 wide.
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 50, 50, 0, 50, 50]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.verdict, 'BIMODAL_STACKED');
    eq(r.trough_count, 0);
    if (!r.trough_below_both) throw new Error('trough_below_both should be true with trough=0');
  });

  t('Histogram counts match input', () => {
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 10, 20, 5, 30, 15]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    eq(r.reach_histogram[1], 10);
    eq(r.reach_histogram[2], 20);
    eq(r.reach_histogram[3], 5);
    eq(r.reach_histogram[4], 30);
    eq(r.reach_histogram[5], 15);
    eq(r.n_valid, 80);
  });

  t('K=6 mode (multi-haplotype) works', () => {
    // 226-fish realistic scenario at K=6 with clear bimodal signature
    // 80 narrow (reach 1) + 50 narrow (reach 2) = 130 narrow
    // 16 mid (reach 3)
    // 30 + 30 + 20 = 80 wide
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 80, 50, 16, 30, 30, 20]), 6);
    const r = w.bandReachBimodality([0,1,2]);
    // narrow=130/226=0.575, wide=80/226=0.354, trough=16/226=0.071.
    // Both exceed 0.25 mins. Smaller peak=80, dip thr=0.75*80=60. trough=16<60 OK.
    eq(r.K, 6);
    eq(r.verdict, 'BIMODAL_STACKED',
       `narrow=${r.narrow_fraction}, wide=${r.wide_fraction}, reason=${r.reason}`);
  });

  t('Empty / null L2 indices returns null gracefully', () => {
    w.computeBandReachAcrossL2s = function() { return null; };
    const r = w.bandReachBimodality([]);
    if (r !== null) throw new Error(`expected null, got ${JSON.stringify(r)}`);
  });

  t('Verdict reason string is informative', () => {
    w.computeBandReachAcrossL2s = makeStubReach(
      buildReach([0, 50, 50, 0, 50, 50]), 5);
    const r = w.bandReachBimodality([0,1,2]);
    if (typeof r.reason !== 'string' || r.reason.length < 10) {
      throw new Error(`reason too short: "${r.reason}"`);
    }
  });

  t('Constants are exposed on window', () => {
    eq(typeof w._BR_BIMODAL_NARROW_MIN, 'number');
    eq(typeof w._BR_BIMODAL_WIDE_MIN, 'number');
    eq(typeof w._BR_BIMODAL_DIP_FRAC, 'number');
    eq(typeof w._BR_NARROW_REACH_MAX, 'number');
    eq(typeof w._BR_WIDE_REACH_MIN, 'number');
    eq(typeof w._BR_TROUGH_REACH, 'number');
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
