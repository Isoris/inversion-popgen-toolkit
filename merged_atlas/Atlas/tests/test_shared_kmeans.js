// tests/test_shared_kmeans.js

import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}
function approx(a, b, tol = 1e-9) {
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  return Math.abs(a - b) <= tol;
}

console.log('--- kmeans1D ---');
{
  // Empty input
  const r0 = kmeans1D([], 3);
  check('empty input → labels.length=0',  r0.labels.length === 0);
  check('empty input → centers.length=3', r0.centers.length === 3);
  check('empty input → n_per_group=[0,0,0]',
        r0.n_per_group.length === 3 && r0.n_per_group.every(c => c === 0));

  // 3 well-separated clusters
  const vals = [0, 0.1, 0.2,  5, 5.1, 5.2,  10, 10.1, 10.2];
  const r = kmeans1D(vals, 3);
  check('labels Int8Array',                 r.labels instanceof Int8Array);
  check('centers Float64Array',             r.centers instanceof Float64Array);
  check('3 centers',                        r.centers.length === 3);
  check('centers ascending',                r.centers[0] < r.centers[1] && r.centers[1] < r.centers[2]);
  check('center 0 ≈ 0.1',                   approx(r.centers[0], 0.1, 1e-9));
  check('center 1 ≈ 5.1',                   approx(r.centers[1], 5.1, 1e-9));
  check('center 2 ≈ 10.1',                  approx(r.centers[2], 10.1, 1e-9));
  check('cluster 0 has 3 members',           r.n_per_group[0] === 3);
  check('cluster 1 has 3 members',           r.n_per_group[1] === 3);
  check('cluster 2 has 3 members',           r.n_per_group[2] === 3);
  check('label 0 = lowest value',            r.labels[0] === 0);
  check('label 1 = middle value',            r.labels[3] === 1);
  check('label 2 = highest value',           r.labels[6] === 2);
}

{
  // Identical values → all in same cluster (after sort, label 0)
  const r = kmeans1D([5, 5, 5, 5, 5], 1);
  check('K=1 identical → all label 0',
        Array.from(r.labels).every(l => l === 0));
  check('K=1 center = 5',                    approx(r.centers[0], 5));
}

console.log('\n--- kmeans2D ---');
{
  // 3 clusters in 2D
  const xs = [0, 0.1, 0.2,  5, 5.1, 5.2,  10, 10.1, 10.2];
  const ys = [0, 0.1, -0.1, 5, 4.9, 5.1, 10, 9.9,  10.1];
  const r = kmeans2D(xs, ys, 3);
  check('labels length n=9',                 r.labels.length === 9);
  check('cx ascending',                       r.cx[0] < r.cx[1] && r.cx[1] < r.cx[2]);
  check('3 members per cluster',              r.n_per_group.every(c => c === 3));
  check('cx[0] ≈ 0.1',                        approx(r.cx[0], 0.1, 1e-9));
  check('cy[1] ≈ 5.0',                        approx(r.cy[1], 5.0, 0.1));
}

{
  // Empty input
  const r = kmeans2D([], [], 3);
  check('empty 2D → labels empty',           r.labels.length === 0);
  check('empty 2D → cx length = 3',          r.cx.length === 3);
}

console.log('\n--- silhouette1D ---');
{
  // 3 well-separated clusters → silhouette near 1
  const vals = [0, 0.1, 0.2, 5, 5.1, 5.2, 10, 10.1, 10.2];
  const r = kmeans1D(vals, 3);
  const s = silhouette1D(vals, r.labels, 3);
  check('clean 3-clusters silhouette > 0.9',  s > 0.9, `got ${s.toFixed(4)}`);

  // Random labels → silhouette near 0 or negative
  const badLabels = new Int8Array([0, 1, 2, 0, 1, 2, 0, 1, 2]);
  const sBad = silhouette1D(vals, badLabels, 3);
  check('shuffled labels silhouette < 0.5',   sBad < 0.5, `got ${sBad.toFixed(4)}`);
}

{
  // n < 4 → NaN
  check('n<4 → NaN',                Number.isNaN(silhouette1D([0,1,2], new Int8Array([0,0,0]), 1)));
  // k < 2 → NaN
  check('k<2 → NaN',                Number.isNaN(silhouette1D([0,1,2,3], new Int8Array([0,0,0,0]), 1)));
  // Singleton cluster → NaN
  check('singleton cluster → NaN',
        Number.isNaN(silhouette1D([0,1,2,3,100], new Int8Array([0,0,0,0,1]), 2)));
}

console.log('\n--- adaptiveK1D ---');
{
  // 3 clean clusters, search [2..5]: should pick K=3
  const vals = [0, 0.1, 0.2, 0.15,  5, 5.1, 5.2, 5.15,  10, 10.1, 10.2, 10.15];
  const r = adaptiveK1D(vals, 2, 5, 0.45, 3);
  check('adaptive picks K=3',          r && r.k === 3, r ? `got K=${r.k}` : 'null');
  check('best silhouette > 0.7',       r && r.silhouette > 0.7);
  check('result has labels',           r && r.labels && r.labels.length === 12);

  // Insufficient data
  check('too few values → null',       adaptiveK1D([1, 2, 3], 2, 5, 0.45, 5) === null);

  // Falls back to kMin if no K passes silhouette threshold
  const noisy = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
  const r2 = adaptiveK1D(noisy, 2, 5, 0.99, 2);  // silThreshold so high nothing passes
  check('falls back to kMin',          r2 && r2.k === 2);
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
