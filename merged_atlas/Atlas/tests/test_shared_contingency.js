// tests/test_shared_contingency.js
// Verifies behavioural parity between shared/contingency.js and the
// legacy Inversion_atlas.html primitives.
//
// Run:  node tests/test_shared_contingency.js

import {
  buildContingency, detectFuseEvents,
  computeARI, computeNMI,
  cramersV, chiSqSurvival, lnGamma,
  scaleStabilityVerdict,
} from '../shared/contingency.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}
function approx(a, b, tol = 1e-9) {
  if (!Number.isFinite(a) && !Number.isFinite(b)) return true;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return false;
  return Math.abs(a - b) <= tol;
}

console.log('--- buildContingency ---');
{
  // Identity partition: KA=KB=3, all i map a=>i, b=>i, n=10
  const A = [0, 0, 1, 1, 2, 2, 0, 1, 2, 0];
  const B = [0, 0, 1, 1, 2, 2, 0, 1, 2, 0];
  const ct = buildContingency(A, B, 3, 3);
  check('null on bad input',                buildContingency(null, [1], 2, 2) === null);
  check('null on length mismatch',          buildContingency([1, 2], [1], 2, 2) === null);
  check('null on KA=0',                     buildContingency([0], [0], 0, 1) === null);
  check('returns shape {M, KA, KB, n}',     ct && ct.M && ct.KA === 3 && ct.KB === 3 && ct.n === 10);
  check('diagonal counts',                  ct.M[0][0] === 4 && ct.M[1][1] === 3 && ct.M[2][2] === 3);
  check('off-diagonal zero',                ct.M[0][1] === 0 && ct.M[1][2] === 0);
}
{
  // Out-of-range labels are skipped, n decreases
  const A = [0, 1, 2, -1, 5];   // -1 and 5 are out-of-range for KA=3
  const B = [0, 1, 2, 0, 1];
  const ct = buildContingency(A, B, 3, 3);
  check('out-of-range skipped, n=3', ct.n === 3);
  check('M[0][0]=1, M[1][1]=1, M[2][2]=1', ct.M[0][0] === 1 && ct.M[1][1] === 1 && ct.M[2][2] === 1);
}

console.log('\n--- detectFuseEvents ---');
{
  // Synthetic: fine clusters {0,1} both map ≥80% into coarse cluster 0
  // Fine cluster 2 maps 100% into coarse cluster 1 (alone — not a fuse)
  // M shape: 3 rows (fine) × 2 cols (coarse)
  const ct = {
    M: [
      [10, 0],   // fine 0 → 100% coarse 0
      [9, 1],    // fine 1 → 90% coarse 0
      [0, 8],    // fine 2 → 100% coarse 1
    ],
    KA: 3, KB: 2, n: 28,
  };
  const fuses = detectFuseEvents(ct, { thresh: 0.80 });
  check('one fuse event detected',          fuses.length === 1);
  check('fuse coarse cluster = 0',          fuses[0].coarse_cluster === 0);
  check('fuse fine clusters = [0, 1]',
        fuses[0].fine_clusters.length === 2
        && fuses[0].fine_clusters.includes(0)
        && fuses[0].fine_clusters.includes(1));
  // Higher threshold filters fine 1 (only 90%)
  const fuses2 = detectFuseEvents(ct, { thresh: 0.95 });
  check('higher thresh excludes fine 1',    fuses2.length === 0);
  check('default thresh is 0.80',           detectFuseEvents(ct).length === 1);
  check('null input → []',                  detectFuseEvents(null).length === 0);
}

console.log('\n--- computeARI ---');
{
  // Identical partitions → ARI = 1
  check('identical partitions → 1',
        approx(computeARI([0,0,1,1,2,2], [0,0,1,1,2,2]), 1));
  // Renamed labels (still identical clustering) → ARI = 1
  check('renamed labels → 1',
        approx(computeARI([0,0,1,1,2,2], [2,2,0,0,1,1]), 1));
  // Independent → ARI ≈ 0
  // n=8, A: 4×0, 4×1; B: alternating 0/1 (independent of A) → ARI should be near 0
  check('NaN on length mismatch',           Number.isNaN(computeARI([0,1], [0])));
  check('NaN on n<2',                       Number.isNaN(computeARI([0], [0])));
  // Manual: A=[0,0,0,1,1,1], B=[0,1,2,0,1,2]: every cell has count 1
  // → sumCellPairs=0, sumRowPairs=2*C2(3)=6, sumColPairs=3*C2(2)=3
  // → totalPairs=C2(6)=15, expected=18/15=1.2, max=4.5
  // → ARI = (0 - 1.2) / (4.5 - 1.2) ≈ -0.3636
  const ari_split = computeARI([0,0,0,1,1,1], [0,1,2,0,1,2]);
  check('split case ≈ -0.3636',             approx(ari_split, -1.2 / 3.3, 1e-6));
}

console.log('\n--- computeNMI ---');
{
  check('identical → 1',                    approx(computeNMI([0,0,1,1,2,2], [0,0,1,1,2,2]), 1));
  check('renamed → 1',                      approx(computeNMI([0,0,1,1,2,2], [2,2,0,0,1,1]), 1));
  check('NaN on bad input',                 Number.isNaN(computeNMI([0,1], [0])));
  // Worst-case independent (A all 0s, B all distinct) → NMI = 0
  // A = [0,0,0,0], B = [0,0,0,0]: HA=HB=0 → returns 1 (both single cluster)
  check('both single-cluster → 1',          computeNMI([0,0,0], [0,0,0]) === 1);
}

console.log('\n--- cramersV ---');
{
  // Independence: row*col = expected exactly → V = 0
  // 2x2 with M = [[10,10],[10,10]]: chi2 = 0 → V = 0
  const flatTable = [10, 10, 10, 10];
  check('independence → V = 0',             approx(cramersV(flatTable, 2, 2), 0));
  // Perfect 1-1 mapping: M = [[5,0],[0,5]]: chi2 = 10 (compute)
  // Actually for 2×2 perfect: V should be 1
  const perfect = [5, 0, 0, 5];
  check('perfect 2x2 → V = 1',              approx(cramersV(perfect, 2, 2), 1));
  // n=0 → NaN
  check('n=0 → NaN',                        Number.isNaN(cramersV([0,0,0,0], 2, 2)));
  // Single non-empty row → V = 0 (degenerate)
  check('single non-empty row → 0',         cramersV([5, 5, 0, 0], 2, 2) === 0);
}

console.log('\n--- lnGamma & chiSqSurvival ---');
{
  // lnGamma well-known values
  check('lnGamma(1) = 0',                   approx(lnGamma(1), 0, 1e-9));
  check('lnGamma(2) = 0',                   approx(lnGamma(2), 0, 1e-9));
  check('lnGamma(3) = ln(2)',               approx(lnGamma(3), Math.log(2), 1e-9));
  check('lnGamma(4) = ln(6)',               approx(lnGamma(4), Math.log(6), 1e-9));
  check('lnGamma(5) = ln(24)',              approx(lnGamma(5), Math.log(24), 1e-9));
  // chi-sq survival: for chi2=0 → 1; for very large chi2 → 0
  check('Q(chi2=0, df=1) = 1',              chiSqSurvival(0, 1) === 1);
  check('Q(chi2=1000, df=1) ≈ 0',           chiSqSurvival(1000, 1) < 1e-100);
  // For chi2=df, Q is around 0.3-0.5 (df=1: ~0.317; df=2: ~0.368; df=4: ~0.406)
  check('Q(chi2=1, df=1) ≈ 0.317',          approx(chiSqSurvival(1, 1), 0.31731, 1e-3));
  check('Q(chi2=2, df=2) ≈ 0.368',          approx(chiSqSurvival(2, 2), 0.36788, 1e-3));
  check('Q(chi2<0) → 1',                    chiSqSurvival(-1, 2) === 1);
  check('Q(df=0) → NaN',                    Number.isNaN(chiSqSurvival(1, 0)));
}

console.log('\n--- scaleStabilityVerdict ---');
{
  // STABLE_3BAND: all panes K=3, ARI≥0.85, no fuses/splits
  const stable3 = scaleStabilityVerdict(
    [{K:3, ok:true}, {K:3, ok:true}, {K:3, ok:true}],
    [{ari:0.95, fuseEvents:[], splitEvents:[]}, {ari:0.92, fuseEvents:[], splitEvents:[]}],
  );
  check('STABLE_3BAND',                     stable3 === 'STABLE_3BAND');
  // STABLE_6BAND
  const stable6 = scaleStabilityVerdict(
    [{K:6, ok:true}, {K:6, ok:true}, {K:6, ok:true}],
    [{ari:0.92, fuseEvents:[], splitEvents:[]}, {ari:0.90, fuseEvents:[], splitEvents:[]}],
  );
  check('STABLE_6BAND',                     stable6 === 'STABLE_6BAND');
  // NESTED_3IN6: K=[3,6,3], 3 fuses each side, no splits
  const nested = scaleStabilityVerdict(
    [{K:3, ok:true}, {K:6, ok:true}, {K:3, ok:true}],
    [
      {ari:0.6, fuseEvents:[{},{},{}], splitEvents:[]},
      {ari:0.6, fuseEvents:[{},{},{}], splitEvents:[]},
    ],
  );
  check('NESTED_3IN6',                      nested === 'NESTED_3IN6');
  // OVERLAP_BREAKS_3: K=[3,6,3] but with fuses+splits at middle
  const overlap = scaleStabilityVerdict(
    [{K:3, ok:true}, {K:6, ok:true}, {K:3, ok:true}],
    [
      {ari:0.75, fuseEvents:[{}], splitEvents:[{}]},
      {ari:0.75, fuseEvents:[{}], splitEvents:[{}]},
    ],
  );
  check('OVERLAP_BREAKS_3',                 overlap === 'OVERLAP_BREAKS_3');
  // UNSTABLE: low ARI
  const unstable = scaleStabilityVerdict(
    [{K:3, ok:true}, {K:3, ok:true}, {K:3, ok:true}],
    [{ari:0.3, fuseEvents:[], splitEvents:[]}, {ari:0.4, fuseEvents:[], splitEvents:[]}],
  );
  check('UNSTABLE on low ARI',              unstable === 'UNSTABLE');
  // UNSTABLE on bad input
  check('UNSTABLE on null panes',           scaleStabilityVerdict(null, []) === 'UNSTABLE');
  check('UNSTABLE on wrong arity',          scaleStabilityVerdict(
        [{K:3, ok:true}, {K:3, ok:true}], []) === 'UNSTABLE');
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
