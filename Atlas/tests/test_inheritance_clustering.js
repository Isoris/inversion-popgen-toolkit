// Tests for inheritance group clustering (turn 2b, Jaccard-based).
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
    try { fn(); console.log(`  PASS [${testNum}] ${name}`); }
    catch (e) { failures.push({ name, err: e.message }); console.log(`  FAIL [${testNum}] ${name}: ${e.message}`); }
  }
  function eq(a, b, m) { if (a !== b) throw new Error(`${m||''} expected ${JSON.stringify(b)}, got ${JSON.stringify(a)}`); }
  function approx(a, b, t, m) { if (!isFinite(a) || Math.abs(a-b) > t) throw new Error(`${m||''} expected ~${b} (tol ${t}), got ${a}`); }

  if (typeof w.inheritanceGroupClustering !== 'function'
      || typeof w._buildBandFishMask !== 'function'
      || typeof w._jaccardDistance !== 'function') {
    console.log('  FAIL: required functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  inheritanceGroupClustering + helpers exposed');

  t('jaccard: identical sets distance = 0', () => {
    const a = new Uint8Array([1,1,0,1,0]);
    const b = new Uint8Array([1,1,0,1,0]);
    approx(w._jaccardDistance(a, b), 0, 1e-9);
  });

  t('jaccard: disjoint sets distance = 1', () => {
    const a = new Uint8Array([1,1,0,0,0]);
    const b = new Uint8Array([0,0,1,1,1]);
    approx(w._jaccardDistance(a, b), 1, 1e-9);
  });

  t('jaccard: half overlap distance = 2/3', () => {
    const a = new Uint8Array([1,1,0,0]);
    const b = new Uint8Array([0,1,1,0]);
    approx(w._jaccardDistance(a, b), 2/3, 1e-9);
  });

  t('jaccard: empty masks distance = 1', () => {
    const a = new Uint8Array([0,0,0]);
    const b = new Uint8Array([0,0,0]);
    eq(w._jaccardDistance(a, b), 1);
  });

  t('jaccard: shape mismatch returns NaN', () => {
    const a = new Uint8Array([1,0]);
    const b = new Uint8Array([1,0,1]);
    if (!isNaN(w._jaccardDistance(a, b))) throw new Error('expected NaN');
  });

  t('fish mask: extracts correct fish for band', () => {
    const items = [{ id: 'A', labels: [0, 1, 2, 0, 1, 2, 0, 0], K: 3 }];
    const mask = w._buildBandFishMask(items, [3], 0, 0);
    eq(mask.length, 8);
    eq(mask[0], 1); eq(mask[3], 1); eq(mask[6], 1); eq(mask[7], 1);
    eq(mask[1], 0); eq(mask[2], 0); eq(mask[4], 0); eq(mask[5], 0);
    eq(mask._count, 4);
  });

  t('fish mask: invalid band returns null', () => {
    const items = [{ id: 'A', labels: [0,1,2], K: 3 }];
    const mask = w._buildBandFishMask(items, [3], 0, 5);
    if (mask !== null) throw new Error('expected null');
  });

  t('agglomerative: 4 points cluster correctly', () => {
    const N = 4;
    const dm = new Float32Array(16);
    function set(i, j, d) { dm[i*N+j] = d; dm[j*N+i] = d; }
    set(0, 1, 0.05); set(2, 3, 0.05);
    set(0, 2, 0.9); set(0, 3, 0.9); set(1, 2, 0.9); set(1, 3, 0.9);
    const dend = w._agglomerativeAverageLinkage(dm, N);
    eq(dend.length, 3);
    const firstTwoDists = [dend[0].dist, dend[1].dist].sort();
    approx(firstTwoDists[0], 0.05, 1e-6);
    approx(firstTwoDists[1], 0.05, 1e-6);
    approx(dend[2].dist, 0.9, 1e-6);
  });

  t('cut: medium threshold gives 2 groups', () => {
    const N = 4;
    const dend = [
      { left: 0, right: 1, dist: 0.1, size: 2, members: [0, 1] },
      { left: 2, right: 3, dist: 0.2, size: 2, members: [2, 3] },
      { left: 4, right: 5, dist: 0.8, size: 4, members: [0, 1, 2, 3] },
    ];
    const c = w._cutDendrogram(dend, N, 0.5);
    eq(c.n_groups, 2);
    eq(c.group_id_per_band[0], c.group_id_per_band[1]);
    eq(c.group_id_per_band[2], c.group_id_per_band[3]);
  });

  t('IGC: 2 perfectly linked candidates -> 3 inheritance groups', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.n_bands_total, 6);
    eq(r.cut.n_groups, 3);
    eq(r.rtab.per_item_n_groups[0], 3);
    eq(r.rtab.per_item_n_groups[1], 3);
    for (const gid of r.rtab.group_ids) {
      eq(Object.keys(r.rtab.per_group[gid]).length, 2);
    }
  });

  t('IGC: 2 independent candidates -> 6 groups (no merging)', () => {
    const labA = [], labB = [];
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        for (let i = 0; i < 22; i++) { labA.push(a); labB.push(b); }
      }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.cut.n_groups, 6);
  });

  t('IGC: 3 candidates, 2 linked + 1 independent', () => {
    const labA = [], labB = [], labC = [];
    for (let a = 0; a < 3; a++) {
      for (let c = 0; c < 3; c++) {
        for (let i = 0; i < 25; i++) {
          labA.push(a); labB.push(a); labC.push(c);
        }
      }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
      { id: 'C', labels: labC, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.cut.n_groups, 6);
    eq(r.rtab.per_item_n_groups[0], 3);
    eq(r.rtab.per_item_n_groups[1], 3);
    eq(r.rtab.per_item_n_groups[2], 3);
    for (let b = 0; b < 3; b++) {
      const idxA = b, idxB = 3 + b;
      eq(r.cut.group_id_per_band[idxA], r.cut.group_id_per_band[idxB]);
    }
    for (let b = 0; b < 3; b++) {
      const idxA = b, idxC = 6 + b;
      if (r.cut.group_id_per_band[idxA] === r.cut.group_id_per_band[idxC]) {
        throw new Error(`A${b} and C${b} should not share group`);
      }
    }
  });

  t('IGC: K=6 multi-haplotype scenario, 2 linked candidates -> 6 groups', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 6; k++) {
      for (let i = 0; i < 20; i++) { labA.push(k); labB.push(k); }
    }
    const items = [
      { id: 'A', labels: labA, K: 6 },
      { id: 'B', labels: labB, K: 6 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.n_bands_total, 12);
    eq(r.cut.n_groups, 6);
    for (const gid of r.rtab.group_ids) {
      eq(Object.keys(r.rtab.per_group[gid]).length, 2);
    }
  });

  t('IGC: noisy linkage tolerated by looser threshold', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 27; i++) { labA.push(k); labB.push(k); }
      for (let i = 0; i < 3; i++)  { labA.push(k); labB.push((k+1)%3); }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const rTight = w.inheritanceGroupClustering(items, { threshold: 0.10 });
    if (rTight.cut.n_groups < 6) throw new Error(`tight threshold should not merge, got ${rTight.cut.n_groups}`);
    const rLoose = w.inheritanceGroupClustering(items, { threshold: 0.25 });
    eq(rLoose.cut.n_groups, 3);
  });

  t('IGC: threshold parameter works', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const r1 = w.inheritanceGroupClustering(items, { threshold: 0.0001 });
    eq(r1.cut.n_groups, 3);
    // To force everything to merge, threshold must be > 1 (Jaccard dist
    // between disjoint within-candidate bands is exactly 1).
    const r2 = w.inheritanceGroupClustering(items, { threshold: 1.01 });
    eq(r2.cut.n_groups, 1);
  });

  t('IGC: too few items returns null', () => {
    eq(w.inheritanceGroupClustering([]), null);
    eq(w.inheritanceGroupClustering([{id: 'A', labels: [0,1,2], K: 3}]), null);
  });

  t('IGC: dendrogram has N-1 merges', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.dendrogram.length, r.n_bands_total - 1);
  });

  t('IGC: realistic 226-fish 3-candidate scenario (lower-res C)', () => {
    const N = 226;
    const labA = [], labB = [], labC = [];
    const groups = [40, 38, 38, 38, 36, 36];
    let idx = 0;
    for (let g = 0; g < 6; g++) {
      for (let i = 0; i < groups[g]; i++) {
        labA.push(g);
        labB.push(g);
        labC.push(Math.floor(g / 2));
        idx++;
      }
    }
    eq(idx, N);
    const items = [
      { id: 'A', labels: labA, K: 6 },
      { id: 'B', labels: labB, K: 6 },
      { id: 'C', labels: labC, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    // A and B merge perfectly → 6 groups. C bands are supersets of A/B bands;
    // Jaccard distance ≈ 0.49, > default threshold 0.15, so C stays separate
    // → +3 more groups = 9 total.
    eq(r.cut.n_groups, 9);
    eq(r.rtab.per_item_n_groups[0], 6);
    eq(r.rtab.per_item_n_groups[1], 6);
    eq(r.rtab.per_item_n_groups[2], 3);
  });

  t('IGC: Rtab structure is correct', () => {
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) {
      for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    }
    const items = [
      { id: 'A', labels: labA, K: 3 },
      { id: 'B', labels: labB, K: 3 },
    ];
    const r = w.inheritanceGroupClustering(items);
    eq(r.rtab.group_ids.length, 3);
    for (const gid of r.rtab.group_ids) {
      const g = r.rtab.per_group[gid];
      eq(Object.keys(g).length, 2);
      eq(g[0].length, 1);
      eq(g[1].length, 1);
    }
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
