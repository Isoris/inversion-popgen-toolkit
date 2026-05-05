// =============================================================================
// turn 150 — Slab-aware U/V cluster modes
// =============================================================================
// Closes the U/V gap turn 149 visibly flagged with the "↩ slab fallback"
// notice. The five U/V modes (uv-rotated, uv-denoise, uv-dbscan,
// uv-dist-rank, uv-dist-fuzzy, plus the distance-uv alias) now run their
// real algorithms against slab ranges instead of falling back to kmeans-K3.
//
// Architecture:
//   - _aggregateWindowRangeForUV(s, e) — extracted phase 1 (per-sample mean
//     PC1*sign / PC2 across window range). Shared between L2 and slab.
//   - _computeUVRotationCore(xs, ys, nS) — extracted phase 2 (K-means3 +
//     principal-axis rotation → (u, v) coords). Shared.
//   - _getOrComputeUVRotation(l2idx) — refactored: aggregation + core,
//     L2 cache (state.l2UVRotationCache).
//   - _getOrComputeUVRotationSlab(s, e) — NEW: aggregation + core,
//     slab cache (state.slabUVRotationCache).
//   - _clusterFromRotation_UV{Rotated,Denoise,DBSCAN,DistRank,DistFuzzy}(rot)
//     — extracted phase 3 (post-rotation clustering). Shared.
//   - clusterL2_UV* / clusterSlab_UV* — thin wrappers calling phase 3 with
//     L2 / slab rotations respectively.
//   - getSlabClusterByMode(s, e, mode) — NEW: dispatcher mirroring
//     getL2ClusterByMode.
//   - compareSlabPair_byMode — rewritten: U/V modes route through
//     getSlabClusterByMode and build their own contingency. fellBack:false
//     for all known modes; fellBack:true only on genuine cluster failure.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Phase-1 / phase-2 extraction — shared aggregation + core
// ============================================================================
console.log('\n=== 1. Aggregation + rotation core extracted ===');

ok('_aggregateWindowRangeForUV(s, e) function declared',
   /function\s+_aggregateWindowRangeForUV\s*\(\s*s\s*,\s*e\s*\)/.test(html));

ok('_computeUVRotationCore(xs, ys, nS) function declared',
   /function\s+_computeUVRotationCore\s*\(\s*xs\s*,\s*ys\s*,\s*nS\s*\)/.test(html));

const aggFnRe = /function\s+_aggregateWindowRangeForUV\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const aggMatch = html.match(aggFnRe);
ok('_aggregateWindowRangeForUV body extracted', !!aggMatch);
const aggBody = aggMatch ? aggMatch[1] : '';

ok('aggregation always uses MEAN (no aggMethod branching)',
   !/state\.aggMethod/.test(aggBody),
   'rotation pipeline uses mean unconditionally per L2 convention');

ok('aggregation iterates pc1[si] * sign and accumulates pc2[si]',
   /pc1\[si\]\s*\*\s*sign[\s\S]{0,200}?pc2\[si\]/.test(aggBody));

ok('aggregation divides by nW for mean',
   /xs\[si\]\s*\/=\s*nW[\s\S]{0,80}?ys\[si\]\s*\/=\s*nW/.test(aggBody));

const coreFnRe = /function\s+_computeUVRotationCore\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const coreMatch = html.match(coreFnRe);
ok('_computeUVRotationCore body extracted', !!coreMatch);
const coreBody = coreMatch ? coreMatch[1] : '';

ok('core returns ok:true with us, vs, angle, centroids',
   /return\s*\{[\s\S]{0,800}?ok:\s*true,[\s\S]{0,200}?us,\s*vs,[\s\S]{0,200}?angle:/.test(coreBody));

ok('core returns ok:false on KMEANS_FAILED',
   /reason:\s*'KMEANS_FAILED'/.test(coreBody));

ok('core uses dx/dy from cx[2]−cx[0] and cy[2]−cy[0] (Hom1↔Hom2 axis)',
   /dx\s*=\s*k3\.cx\[2\]\s*-\s*k3\.cx\[0\][\s\S]{0,80}?dy\s*=\s*k3\.cy\[2\]\s*-\s*k3\.cy\[0\]/.test(coreBody));

// ============================================================================
// 2. _getOrComputeUVRotation refactor — uses extracted helpers
// ============================================================================
console.log('\n=== 2. L2 rotation refactored to use shared helpers ===');

const l2RotFnRe = /function\s+_getOrComputeUVRotation\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const l2RotMatch = html.match(l2RotFnRe);
ok('_getOrComputeUVRotation body extracted', !!l2RotMatch);
const l2RotBody = l2RotMatch ? l2RotMatch[1] : '';

ok('L2 rotation calls _aggregateWindowRangeForUV',
   /_aggregateWindowRangeForUV\(env\._s0,\s*env\._e0\)/.test(l2RotBody));

ok('L2 rotation calls _computeUVRotationCore',
   /_computeUVRotationCore\(agg\.xs,\s*agg\.ys,\s*d\.n_samples\)/.test(l2RotBody));

ok('L2 rotation still uses state.l2UVRotationCache',
   /state\.l2UVRotationCache\.set\(l2idx,\s*result\)/.test(l2RotBody));

// The old inline aggregation loop should be GONE from L2 rotation
ok('L2 rotation no longer has inline pc1/pc2 aggregation loop (moved to helper)',
   !/for \(let w = 0; w < nW; w\+\+\) \{[\s\S]{0,200}?pc1\[si\]/.test(l2RotBody),
   'aggregation loop should be in _aggregateWindowRangeForUV, not inlined here');

// ============================================================================
// 3. _getOrComputeUVRotationSlab — NEW slab rotation cache
// ============================================================================
console.log('\n=== 3. Slab rotation cache ===');

ok('_getOrComputeUVRotationSlab function declared',
   /function\s+_getOrComputeUVRotationSlab\s*\(\s*s\s*,\s*e\s*\)/.test(html));

const slabRotFnRe = /function\s+_getOrComputeUVRotationSlab\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const slabRotMatch = html.match(slabRotFnRe);
ok('_getOrComputeUVRotationSlab body extracted', !!slabRotMatch);
const slabRotBody = slabRotMatch ? slabRotMatch[1] : '';

ok('slab rotation uses state.slabUVRotationCache (NOT state.l2UVRotationCache)',
   /state\.slabUVRotationCache/.test(slabRotBody) &&
   !/state\.l2UVRotationCache/.test(slabRotBody));

ok('slab rotation cache key includes both s and e',
   /cacheKey\s*=\s*`\$\{s\}_\$\{e\}`/.test(slabRotBody));

ok('slab rotation calls _aggregateWindowRangeForUV with s, e',
   /_aggregateWindowRangeForUV\(s,\s*e\)/.test(slabRotBody));

ok('slab rotation calls _computeUVRotationCore',
   /_computeUVRotationCore\(agg\.xs,\s*agg\.ys,\s*d\.n_samples\)/.test(slabRotBody));

ok('slab rotation invalidates on data identity change',
   /_slabUVRotationCacheDataKey/.test(slabRotBody));

ok('slab rotation defends against bad ranges (s < 0, e >= n_windows, s > e)',
   /BAD_RANGE/.test(slabRotBody));

// ============================================================================
// 4. clusterSlab_UV* — five new slab clusterers
// ============================================================================
console.log('\n=== 4. Slab U/V clusterers (5 modes) ===');

for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const fn = new RegExp(`function\\s+clusterSlab_${mode}\\s*\\(\\s*s\\s*,\\s*e\\s*\\)`);
  ok(`clusterSlab_${mode}(s, e) declared`, fn.test(html));
}

// Each slab clusterer should be a thin wrapper: get rotation, call
// _clusterFromRotation_UV* core. Verify by extracting bodies.
for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const re = new RegExp(`function\\s+clusterSlab_${mode}\\s*\\([^)]*\\)\\s*\\{([\\s\\S]*?)\\n\\}`);
  const m = html.match(re);
  if (m) {
    ok(`clusterSlab_${mode} body uses _getOrComputeUVRotationSlab`,
       /_getOrComputeUVRotationSlab\(s,\s*e\)/.test(m[1]));
    ok(`clusterSlab_${mode} body calls _clusterFromRotation_${mode}`,
       new RegExp(`_clusterFromRotation_${mode}\\(rot\\)`).test(m[1]));
  }
}

// ============================================================================
// 5. _clusterFromRotation_UV* shared cores
// ============================================================================
console.log('\n=== 5. Shared post-rotation cores ===');

for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const fn = new RegExp(`function\\s+_clusterFromRotation_${mode}\\s*\\(\\s*rot\\s*\\)`);
  ok(`_clusterFromRotation_${mode}(rot) declared`, fn.test(html));
}

// Each core should call _wrapKmeansResultAsCluster with null l2idx
// turn 150: count using non-greedy multiline match since the first arg can
// be { labels, n_per_group: npg } (object literal) which `[^,]*` can't span.
const wrapNullMatches = html.match(/_wrapKmeansResultAsCluster\([\s\S]*?,\s*null,\s*3/g) || [];
ok('all _wrapKmeansResultAsCluster calls in cores pass null for l2idx (>= 5 expected)',
   wrapNullMatches.length >= 5,
   `found ${wrapNullMatches.length} matches`);

// ============================================================================
// 6. clusterL2_UV* still exist (thin wrappers now)
// ============================================================================
console.log('\n=== 6. L2 clusterers preserved as thin wrappers ===');

for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const fn = new RegExp(`function\\s+clusterL2_${mode}\\s*\\(\\s*l2idx\\s*\\)`);
  ok(`clusterL2_${mode}(l2idx) preserved`, fn.test(html));
}

// Each L2 wrapper now should be just two lines: get rotation, call core
for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const re = new RegExp(`function\\s+clusterL2_${mode}\\s*\\([^)]*\\)\\s*\\{([\\s\\S]*?)\\n\\}`);
  const m = html.match(re);
  if (m) {
    ok(`clusterL2_${mode} now calls shared core via _getOrComputeUVRotation`,
       /_getOrComputeUVRotation\(l2idx\)/.test(m[1]) &&
       new RegExp(`_clusterFromRotation_${mode}\\(rot\\)`).test(m[1]));
  }
}

// ============================================================================
// 7. getSlabClusterByMode — dispatcher
// ============================================================================
console.log('\n=== 7. getSlabClusterByMode dispatcher ===');

ok('getSlabClusterByMode(s, e, mode) declared',
   /function\s+getSlabClusterByMode\s*\(\s*s\s*,\s*e\s*,\s*mode\s*\)/.test(html));

const dispatchFnRe = /function\s+getSlabClusterByMode\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const dispatchMatch = html.match(dispatchFnRe);
ok('getSlabClusterByMode body extracted', !!dispatchMatch);
const dispatchBody = dispatchMatch ? dispatchMatch[1] : '';

ok('kmeans-K3 default → getSlabClusterAt(s, e, state.k)',
   /'kmeans-K3'[\s\S]{0,100}?getSlabClusterAt\(s,\s*e,\s*state\.k\)/.test(dispatchBody));

ok('kmeans-K6 → getSlabClusterAt(s, e, 6)',
   /'kmeans-K6'[\s\S]{0,100}?getSlabClusterAt\(s,\s*e,\s*6\)/.test(dispatchBody));

ok('uses state.slabGroupCacheByMode for U/V modes',
   /state\.slabGroupCacheByMode/.test(dispatchBody));

ok('cache key for U/V modes is "${s}_${e}_${mode}"',
   /cacheKey\s*=\s*`\$\{s\}_\$\{e\}_\$\{mode\}`/.test(dispatchBody));

for (const [optionVal, fnSuffix] of [
  ['uv-rotated', 'UVRotated'],
  ['distance-uv', 'UVRotated'],          // alias
  ['uv-denoise', 'UVDenoise'],
  ['uv-dbscan', 'UVDBSCAN'],
  ['uv-dist-rank', 'UVDistRank'],
  ['uv-dist-fuzzy', 'UVDistFuzzy'],
]) {
  const re = new RegExp(`case '${optionVal.replace(/-/g, '\\-')}':[\\s\\S]{0,300}?clusterSlab_${fnSuffix}\\(s,\\s*e\\)`);
  ok(`${optionVal} → clusterSlab_${fnSuffix}`, re.test(dispatchBody));
}

ok('default case falls back to getSlabClusterAt at state.k (not crash)',
   /default:[\s\S]{0,500}?getSlabClusterAt\(s,\s*e,\s*state\.k\)/.test(dispatchBody));

// ============================================================================
// 8. compareSlabPair_byMode rewritten — U/V actually runs (not fallback)
// ============================================================================
console.log('\n=== 8. compareSlabPair_byMode now runs U/V for real ===');

const cmpFnRe = /function\s+compareSlabPair_byMode\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const cmpMatch = html.match(cmpFnRe);
const cmpBody = cmpMatch ? cmpMatch[1] : '';

ok('U/V branch calls getSlabClusterByMode for both ranges',
   /getSlabClusterByMode\(leftRange\[0\],\s+leftRange\[1\],\s+m\)[\s\S]{0,200}?getSlabClusterByMode\(rightRange\[0\],\s*rightRange\[1\],\s*m\)/.test(cmpBody));

ok('U/V branch uses alignLabels to build contingency',
   /alignLabels\(llab,\s*rlab,\s*K\)/.test(cmpBody));

ok('U/V branch sets verdict via mergeThr / LOW_POWER (parity with compareSlabPair)',
   /verdict\s*=\s*'LOW_POWER'[\s\S]{0,200}?verdict\s*=\s*'MERGE'[\s\S]{0,200}?verdict\s*=\s*'SEPARATE'/.test(cmpBody));

ok('U/V success path returns isSlabPair: true',
   /isSlabPair:\s*true/.test(cmpBody));

ok('U/V success path returns reclusterMode: m, fellBack',
   /reclusterMode:\s*m,[\s\S]{0,200}?fellBack/.test(cmpBody));

ok('U/V failure path keeps fellBack:true for genuine cluster failures',
   /if \(!cl\s*\|\|\s*!cr[\s\S]{0,400}?fellBack\s*=\s*true/.test(cmpBody));

// ============================================================================
// 9. Sandboxed run of the new slab U/V pipeline end-to-end
// ============================================================================
console.log('\n=== 9. Sandboxed end-to-end: rotation core + dispatcher ===');

{
  // Set up a synthetic 3-band cohort: 9 samples, K=3 well-separated clusters
  // in PC1xPC2 space. Test that _computeUVRotationCore correctly identifies
  // the principal axis and rotates Hom1↔Hom2 to horizontal.
  const sandbox = {
    state: { k: 3, mergeThr: 0.85 },
    console,
  };
  const ctx = vm.createContext(sandbox);

  // We need kmeans2D + the rotation core. Extract them.
  const k2dMatch = html.match(/(function\s+kmeans2D\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  ok('extracted kmeans2D', !!k2dMatch);

  const coreMatch2 = html.match(/(function\s+_computeUVRotationCore\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  ok('extracted _computeUVRotationCore', !!coreMatch2);

  if (k2dMatch && coreMatch2) {
    vm.runInContext(k2dMatch[1], ctx);
    vm.runInContext(coreMatch2[1], ctx);

    // Synthetic data: 3 clusters at (-1, 0), (0, 0), (1, 0) — already aligned
    // on the x-axis. Rotation should be ~0° and (u, v) should ≈ (xs, ys).
    const setup = `
      const xs = new Float64Array(9);
      const ys = new Float64Array(9);
      // Cluster 0 (Hom1) at (-1, 0)
      for (let i = 0; i < 3; i++) { xs[i] = -1 + i*0.01; ys[i] = i*0.01; }
      // Cluster 1 (Het) at (0, 0)
      for (let i = 3; i < 6; i++) { xs[i] = (i-3)*0.01; ys[i] = (i-3)*0.01; }
      // Cluster 2 (Hom2) at (1, 0)
      for (let i = 6; i < 9; i++) { xs[i] = 1 + (i-6)*0.01; ys[i] = (i-6)*0.01; }
      const rot = _computeUVRotationCore(xs, ys, 9);
    `;
    vm.runInContext(setup, ctx);

    ok('rotation succeeds for clean 3-band data',
       vm.runInContext(`rot.ok === true`, ctx));
    ok('rotation produces us, vs typed arrays of length 9',
       vm.runInContext(`rot.us instanceof Float64Array && rot.vs instanceof Float64Array && rot.us.length === 9`, ctx));
    ok('axis is approximately horizontal (angle ≈ 0)',
       vm.runInContext(`Math.abs(rot.angle) < 5 || Math.abs(rot.angle - 180) < 5`, ctx),
       'angle should be near 0° or 180° (sign convention)');
    ok('Hom1 / Hom2 centroids span u-axis (|hom1_u - hom2_u| > 0.5)',
       vm.runInContext(`Math.abs(rot.hom1_u - rot.hom2_u) > 0.5`, ctx));
    ok('baseLabels populated with K=3 cluster assignments',
       vm.runInContext(`rot.baseLabels && rot.baseLabels.length === 9`, ctx));

    // Test rotation of an off-axis data: Hom1 at (-1, -1), Hom2 at (1, 1)
    const setup2 = `
      const xs2 = new Float64Array(9);
      const ys2 = new Float64Array(9);
      for (let i = 0; i < 3; i++) { xs2[i] = -1 + i*0.01; ys2[i] = -1 + i*0.01; }
      for (let i = 3; i < 6; i++) { xs2[i] = (i-3)*0.01; ys2[i] = (i-3)*0.01; }
      for (let i = 6; i < 9; i++) { xs2[i] = 1 + (i-6)*0.01; ys2[i] = 1 + (i-6)*0.01; }
      const rot2 = _computeUVRotationCore(xs2, ys2, 9);
    `;
    vm.runInContext(setup2, ctx);

    ok('rotation succeeds for 45° off-axis data',
       vm.runInContext(`rot2.ok === true`, ctx));
    ok('rotation angle is ~45° for diagonal data',
       vm.runInContext(`Math.abs(Math.abs(rot2.angle) - 45) < 5`, ctx));
    ok('post-rotation, Hom1 and Hom2 are still at extreme u (axis aligned)',
       vm.runInContext(`Math.abs(rot2.hom1_u - rot2.hom2_u) > 1.5`, ctx));
    ok('post-rotation, vs values are small (rotation collapsed perpendicular)',
       vm.runInContext(`Math.max(...Array.from(rot2.vs).map(Math.abs)) < 0.2`, ctx));
  }
}

// ============================================================================
// 10. Sandboxed UV cluster cores — they actually produce clusters
// ============================================================================
console.log('\n=== 10. Sandboxed _clusterFromRotation_UVRotated ===');

{
  const sandbox = { state: { k: 3, minNGroup: 2, data: { samples: [], n_samples: 9 } }, console };
  const ctx = vm.createContext(sandbox);

  const k2dMatch = html.match(/(function\s+kmeans2D\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  const wrapMatch = html.match(/(function\s+_wrapKmeansResultAsCluster\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  const fromRotMatch = html.match(/(function\s+_clusterFromRotation_UVRotated\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);

  if (k2dMatch && wrapMatch && fromRotMatch) {
    vm.runInContext(k2dMatch[1], ctx);
    vm.runInContext(wrapMatch[1], ctx);
    vm.runInContext(fromRotMatch[1], ctx);

    // Synthetic rot result with us/vs already in 3 clean clusters
    const setup = `
      const us = new Float64Array(9);
      const vs = new Float64Array(9);
      for (let i = 0; i < 3; i++) { us[i] = -1 + i*0.01; vs[i] = i*0.01; }
      for (let i = 3; i < 6; i++) { us[i] = (i-3)*0.01; vs[i] = (i-3)*0.01; }
      for (let i = 6; i < 9; i++) { us[i] = 1 + (i-6)*0.01; vs[i] = (i-6)*0.01; }
      const rot = {
        ok: true, us, vs,
        hom1_u: -1, hom1_v: 0, het_u: 0, het_v: 0, hom2_u: 1, hom2_v: 0,
        baseLabels: new Int8Array([0,0,0,1,1,1,2,2,2]),
        baseN: [3, 3, 3], degenerate: false,
      };
      const result = _clusterFromRotation_UVRotated(rot);
    `;
    vm.runInContext(setup, ctx);

    ok('UVRotated produces ok:true cluster',
       vm.runInContext(`result.ok === true`, ctx));
    ok('UVRotated produces K=3 cluster labels',
       vm.runInContext(`result.usedK === 3`, ctx));
    ok('UVRotated produces 3 groups via n_per_group',
       vm.runInContext(`result.n_per_group.length === 3`, ctx));
    ok('UVRotated handles degenerate: true (falls back to baseLabels)',
       vm.runInContext(`
         (function() {
           const degenRot = { ok: true, degenerate: true, baseLabels: new Int8Array([0,0,1,1,2,2,2,2,2]), baseN: [2,2,5] };
           const r = _clusterFromRotation_UVRotated(degenRot);
           return r.ok === true && r.reason === 'UVRotated-degenerate-fallback';
         })()
       `, ctx));
    ok('UVRotated returns ok:false on rot.ok:false',
       vm.runInContext(`
         (function() {
           const failRot = { ok: false, reason: 'KMEANS_FAILED' };
           const r = _clusterFromRotation_UVRotated(failRot);
           return r.ok === false && r.reason === 'KMEANS_FAILED';
         })()
       `, ctx));
  }
}

// ============================================================================
// 11. Sandboxed getSlabClusterByMode — routing test
// ============================================================================
console.log('\n=== 11. Sandboxed getSlabClusterByMode dispatch ===');

{
  const _routes = [];
  const sandbox = {
    state: { k: 3 },
    getSlabClusterAt: (s, e, K) => { _routes.push({fn:'getSlabClusterAt', s, e, K}); return {ok:true, labels:[]}; },
    clusterSlab_UVRotated: (s, e) => { _routes.push({fn:'UVRotated', s, e}); return {ok:true, labels:[]}; },
    clusterSlab_UVDenoise: (s, e) => { _routes.push({fn:'UVDenoise', s, e}); return {ok:true, labels:[]}; },
    clusterSlab_UVDBSCAN: (s, e) => { _routes.push({fn:'UVDBSCAN', s, e}); return {ok:true, labels:[]}; },
    clusterSlab_UVDistRank: (s, e) => { _routes.push({fn:'UVDistRank', s, e}); return {ok:true, labels:[]}; },
    clusterSlab_UVDistFuzzy: (s, e) => { _routes.push({fn:'UVDistFuzzy', s, e}); return {ok:true, labels:[]}; },
    console,
  };
  const ctx = vm.createContext(sandbox);

  const dispatchSrc = html.match(/(function\s+getSlabClusterByMode\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  if (dispatchSrc) {
    vm.runInContext(dispatchSrc[1], ctx);

    for (const [mode, expectedFn] of [
      ['kmeans-K3', 'getSlabClusterAt'],
      ['kmeans-K6', 'getSlabClusterAt'],
      ['uv-rotated', 'UVRotated'],
      ['distance-uv', 'UVRotated'],            // alias
      ['uv-denoise', 'UVDenoise'],
      ['uv-dbscan', 'UVDBSCAN'],
      ['uv-dist-rank', 'UVDistRank'],
      ['uv-dist-fuzzy', 'UVDistFuzzy'],
    ]) {
      _routes.length = 0;
      vm.runInContext(`getSlabClusterByMode(0, 4, '${mode}')`, ctx);
      ok(`mode '${mode}' routes to ${expectedFn}`,
         _routes.length >= 1 && _routes[0].fn === expectedFn);
    }

    // Cache test: second call to same (s, e, mode) should hit cache, not re-call.
    // Use a fresh range NOT tested in the loop above so the cache starts cold.
    _routes.length = 0;
    vm.runInContext(`getSlabClusterByMode(20, 24, 'uv-rotated'); getSlabClusterByMode(20, 24, 'uv-rotated');`, ctx);
    ok('repeat call hits cache (clusterSlab_UVRotated called only once)',
       _routes.length === 1);

    // Different range → fresh call (not from cache)
    _routes.length = 0;
    vm.runInContext(`getSlabClusterByMode(30, 34, 'uv-rotated'); getSlabClusterByMode(35, 39, 'uv-rotated');`, ctx);
    ok('different (s, e) → fresh call',
       _routes.length === 2);
  }
}

// ============================================================================
// 12. compareSlabPair_byMode no longer always returns fellBack:true for U/V
// ============================================================================
console.log('\n=== 12. fellBack:false for U/V when slab cluster succeeds ===');

// turn 149 had: U/V → fellBack:true unconditionally. After turn 150,
// fellBack should be false when getSlabClusterByMode returns valid clusters
// (the normal case).
ok('compareSlabPair_byMode no longer hardcodes fellBack:true in U/V else-branch',
   !/m === 'kmeans-K3'\s*\|\|\s*m === 'kmeans-K6'\)\s*\{[^}]*\}\s*else\s*\{[^}]*?fellBack\s*=\s*true/.test(cmpBody),
   'old turn 149 unconditional fallback structure should be gone');

// ============================================================================
// 13. Slab renderer — fallback notice now rare (only on cluster failure)
// ============================================================================
console.log('\n=== 13. Slab renderer fallback notice still wired ===');

const slabFnRe2 = /function\s+renderL3PanelSlab\s*\(\s*\)\s*\{([\s\S]*?)\n\}\s*\n\s*\/\/[^\n]*Compact focal-content HTML for a slab/;
const slabFnMatch2 = html.match(slabFnRe2);
const slabBody2 = slabFnMatch2 ? slabFnMatch2[1] : '';

// The notice code is unchanged from turn 149 — it reads cmp.fellBack which
// is still produced by compareSlabPair_byMode (just much rarer now).
ok('slab renderer still reads cmp.fellBack to surface notice',
   /if \(cmp\.fellBack/.test(slabBody2));

ok('slab renderer still reads cmp.requestedMode for notice',
   /cmp\.requestedMode/.test(slabBody2));

// ============================================================================
// 14. L2 path unchanged
// ============================================================================
console.log('\n=== 14. L2 mode preserved ===');

ok('getL2ClusterByMode still exists',
   /function\s+getL2ClusterByMode\s*\(\s*l2idx\s*,\s*mode\s*\)/.test(html));

ok('compareL2Pair_byMode still exists',
   /function\s+compareL2Pair_byMode\s*\(\s*leftIdx\s*,\s*rightIdx\s*,\s*mode\s*\)/.test(html));

// L2 wrappers still take l2idx
for (const mode of ['UVRotated', 'UVDenoise', 'UVDBSCAN', 'UVDistRank', 'UVDistFuzzy']) {
  const fn = new RegExp(`function\\s+clusterL2_${mode}\\s*\\(\\s*l2idx\\s*\\)`);
  ok(`clusterL2_${mode} still takes l2idx (L2 path unchanged)`, fn.test(html));
}

// ============================================================================
// Summary
// ============================================================================
console.log(`\n=========================================`);
console.log(`  PASS: ${pass}   FAIL: ${fail}`);
console.log(`=========================================\n`);
process.exit(fail > 0 ? 1 : 0);
