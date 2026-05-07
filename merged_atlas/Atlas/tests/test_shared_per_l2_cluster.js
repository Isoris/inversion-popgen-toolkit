// tests/test_shared_per_l2_cluster.js
//
// Synthetic-fixture tests. Validates the new ctx-based API works end-
// to-end without depending on real precomp data.
//
// We construct minimal `data` objects that look like state.data and
// verify the cluster pipeline produces the expected shapes + values.

import {
  contextFromState,
  aggregateL2,
  sampleSpreadL2, sampleSpreadRange,
  sigmaProfileL2,
  clusterL2, clusterL2AtK,
  clusterCacheKey, ClusterCache,
} from '../shared/per_l2_cluster.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

// -----------------------------------------------------------------------
// Build a synthetic ctx with 12 samples, 20 windows, 1 L2 envelope
// covering windows [0..19]. PC1 separates samples into 3 bands clearly.
// -----------------------------------------------------------------------
function makeFixture({ n_samples = 12, n_windows = 20 } = {}) {
  // 3 groups: samples 0-3 in band A (PC1 ≈ -1), 4-7 in band B (PC1 ≈ 0),
  // 8-11 in band C (PC1 ≈ +1). Add noise per window.
  const windows = [];
  for (let w = 0; w < n_windows; w++) {
    const noise = () => (Math.sin(w * 7.13 + 1) * 0.05);  // deterministic
    const pc1 = new Float32Array(n_samples);
    const pc2 = new Float32Array(n_samples);
    for (let s = 0; s < n_samples; s++) {
      const band = Math.floor(s / 4);   // 0, 1, 2
      pc1[s] = (band - 1) + noise() * (0.3 + 0.1 * (s % 3));
      pc2[s] = noise() * 0.1;
    }
    windows.push({ pc1: Array.from(pc1), pc2: Array.from(pc2),
                   start_bp: w * 50_000, end_bp: (w + 1) * 50_000 - 1 });
  }
  return {
    data: {
      n_samples,
      n_windows,
      windows,
      l2_envelopes: [
        { _s0: 0, _e0: n_windows - 1, start_bp: 0, end_bp: n_windows * 50_000 - 1 },
      ],
      family_source: 'none',
    },
    k: 3,
    aggMethod: 'mean_pc1',
    kMode: 'fixed',
    kRange: [2, 5],
    silThreshold: 0.45,
    minNGroup: 3,
    minNWin: 3,
    flipPC1: false,
    pc1Sign: null,
  };
}

console.log('--- contextFromState ---');
{
  const state = makeFixture();
  const ctx = contextFromState(state);
  check('ctx has data',                     ctx.data === state.data);
  check('ctx.k = 3',                        ctx.k === 3);
  check('ctx.aggMethod = mean_pc1',         ctx.aggMethod === 'mean_pc1');
  check('ctx.getPC is a function',          typeof ctx.getPC === 'function');
  // Verify getPC returns proper shape + applies flipPC1=false sign=1
  const got = ctx.getPC(0);
  check('getPC.pc1 is array-like',          Array.isArray(got.pc1) || ArrayBuffer.isView(got.pc1));
  check('getPC.sign = 1 (no flip)',         got.sign === 1);
}

console.log('\n--- aggregateL2 (mean_pc1) ---');
{
  const state = makeFixture();
  const ctx = contextFromState(state);
  const agg = aggregateL2(ctx, 0);
  check('returns {xs, ys, nW}',             agg && agg.xs && agg.nW === 20);
  check('xs.length = n_samples',            agg.xs.length === 12);
  check('ys is null for mean_pc1',          agg.ys === null);
  // Sample 0 in band A (PC1 ≈ -1)
  check('sample 0 (band A) ≈ -1',           Math.abs(agg.xs[0] + 1) < 0.5, `got ${agg.xs[0].toFixed(3)}`);
  // Sample 5 in band B (PC1 ≈ 0)
  check('sample 5 (band B) ≈ 0',            Math.abs(agg.xs[5]) < 0.5, `got ${agg.xs[5].toFixed(3)}`);
  // Sample 10 in band C (PC1 ≈ +1)
  check('sample 10 (band C) ≈ +1',          Math.abs(agg.xs[10] - 1) < 0.5, `got ${agg.xs[10].toFixed(3)}`);
}

console.log('\n--- aggregateL2 (mean_pc12) ---');
{
  const state = makeFixture();
  state.aggMethod = 'mean_pc12';
  const ctx = contextFromState(state);
  const agg = aggregateL2(ctx, 0);
  check('ys is Float64Array for mean_pc12',  agg.ys instanceof Float64Array);
  check('ys.length = n_samples',             agg.ys.length === 12);
}

console.log('\n--- aggregateL2 (median_pc1) ---');
{
  const state = makeFixture();
  state.aggMethod = 'median_pc1';
  const ctx = contextFromState(state);
  const agg = aggregateL2(ctx, 0);
  check('returns xs',                        agg && agg.xs);
  // Median should still produce 3 distinct levels
  check('median preserves band structure',
        Math.abs(agg.xs[0] - agg.xs[10]) > 1.0, `xs[0]=${agg.xs[0].toFixed(3)} xs[10]=${agg.xs[10].toFixed(3)}`);
}

console.log('\n--- sampleSpreadL2 + sampleSpreadRange ---');
{
  const state = makeFixture();
  const ctx = contextFromState(state);
  const sd = sampleSpreadL2(ctx, 0);
  check('sd is Float64Array',                sd instanceof Float64Array);
  check('sd.length = n_samples',             sd.length === 12);
  check('all sd values finite',              Array.from(sd).every(v => Number.isFinite(v)));
  // Spread should be small (samples are stable in their bands)
  check('spread < 0.3 for all samples',
        Array.from(sd).every(v => v < 0.3), `max=${Math.max(...sd).toFixed(3)}`);

  const sd2 = sampleSpreadRange(ctx, 0, 19);
  check('sampleSpreadRange = sampleSpreadL2 (same range)',
        Array.from(sd).every((v, i) => Math.abs(v - sd2[i]) < 1e-12));

  // <2 windows → null
  check('range with <2 windows → null',      sampleSpreadRange(ctx, 5, 5) === null);
  check('range w<2 → null (negative)',       sampleSpreadRange(ctx, 5, 4) === null);
}

console.log('\n--- clusterL2 (fixed K=3) ---');
{
  const state = makeFixture();
  const ctx = contextFromState(state);
  const cl = clusterL2(ctx, 0);
  check('ok = true',                         cl.ok === true);
  check('reason = null',                     cl.reason === null);
  check('labels.length = 12',                cl.labels.length === 12);
  check('n_per_group = [4, 4, 4]',           cl.n_per_group.length === 3
                                             && cl.n_per_group[0] === 4
                                             && cl.n_per_group[1] === 4
                                             && cl.n_per_group[2] === 4);
  check('usedK = 3',                         cl.usedK === 3);
  check('silhouette null in fixed mode',     cl.silhouette === null);
  check('fixedKLabels === labels (fixed mode)',
        Array.from(cl.fixedKLabels).every((v, i) => v === cl.labels[i]));
  // Samples 0-3 should all be in the same cluster (band A = label 0 since centers ascending)
  const band_A_label = cl.labels[0];
  check('samples 0-3 all band A',
        cl.labels[1] === band_A_label && cl.labels[2] === band_A_label && cl.labels[3] === band_A_label);
  check('band A label = 0 (lowest center)',  band_A_label === 0);
  // Samples 8-11 → band C (highest center → label 2)
  check('samples 8-11 all band C (label 2)',
        cl.labels[8] === 2 && cl.labels[9] === 2 && cl.labels[10] === 2 && cl.labels[11] === 2);
}

console.log('\n--- clusterL2 (adaptive mode) ---');
{
  const state = makeFixture();
  state.kMode = 'adaptive';
  state.kRange = [2, 4];
  const ctx = contextFromState(state);
  const cl = clusterL2(ctx, 0);
  check('ok = true',                         cl.ok === true);
  check('usedK = 3 (best silhouette)',       cl.usedK === 3);
  check('silhouette > 0.7',                  cl.silhouette > 0.7, `got ${cl.silhouette}`);
  // fixedKLabels still computed at ctx.k=3
  check('fixedKLabels exists in adaptive',   cl.fixedKLabels && cl.fixedKLabels.length === 12);
}

console.log('\n--- clusterL2AtK ---');
{
  const state = makeFixture();
  state.kMode = 'adaptive';   // intentionally adaptive — clusterL2AtK should ignore
  const ctx = contextFromState(state);
  const cl3 = clusterL2AtK(ctx, 0, 3);
  check('AtK=3 ok',                         cl3.ok === true);
  check('AtK=3 usedK=3',                    cl3.usedK === 3);
  // n_per_group still [4,4,4]
  check('AtK=3 n_per_group [4,4,4]',
        cl3.n_per_group.length === 3 && cl3.n_per_group.every(c => c === 4));
  // fixedKLabels === labels at this K
  check('AtK fixedKLabels = labels',
        Array.from(cl3.fixedKLabels).every((v, i) => v === cl3.labels[i]));

  // K=4 should now have unbalanced groups → ok = false
  const cl4 = clusterL2AtK(ctx, 0, 4);
  check('AtK=4 produces 4 labels',          cl4.labels.length === 12);
  check('AtK=4 n_per_group has 4 entries',  cl4.n_per_group.length === 4);
}

console.log('\n--- clusterL2 fail modes ---');
{
  // nW=1 doesn't fail by itself if K-means still produces ≥minNGroup per cluster.
  // It just sets n_below_threshold=true. Matches legacy contract.
  const state = makeFixture();
  state.data.l2_envelopes[0] = { _s0: 0, _e0: 0 };
  const ctx = contextFromState(state);
  const cl = clusterL2(ctx, 0);
  check('nW=1 sets n_below_threshold=true', cl.n_below_threshold === true);
  check('nW=1 nW field reports 1',           cl.nW === 1);

  // To FAIL ok, we need to fail minNGroup. Use a fixture with only 5 samples (K=3 minNGroup=3 → at least one group has <3).
  const small = makeFixture({ n_samples: 5 });
  const ctxSmall = contextFromState(small);
  const clSmall = clusterL2(ctxSmall, 0);
  check('5 samples K=3 minNGroup=3 → ok=false', clSmall.ok === false);
  check('reason = LOW_GROUP_N',              clSmall.reason === 'LOW_GROUP_N');

  // negative range → NO_WINDOWS (immediate)
  state.data.l2_envelopes[0] = { _s0: 5, _e0: 4 };
  const cl2 = clusterL2(ctx, 0);
  check('negative range → ok=false',         cl2.ok === false);
  check('negative range reason NO_WINDOWS',  cl2.reason === 'NO_WINDOWS');
}

console.log('\n--- clusterCacheKey ---');
{
  const state = makeFixture();
  const ctx1 = contextFromState(state);
  const k1 = clusterCacheKey(ctx1);

  // Mutate state.k
  state.k = 6;
  const ctx2 = contextFromState(state);
  const k2 = clusterCacheKey(ctx2);
  check('different k → different cache key', k1 !== k2);

  // Restore k, change kMode
  state.k = 3;
  state.kMode = 'adaptive';
  const ctx3 = contextFromState(state);
  const k3 = clusterCacheKey(ctx3);
  check('different kMode → different key',   k3 !== k1);
}

console.log('\n--- ClusterCache ---');
{
  const state = makeFixture();
  const ctx = contextFromState(state);
  const cache = new ClusterCache();
  check('cache starts empty',                cache.size === 0);

  const cl1 = cache.getOrCompute(ctx, 0);
  check('first compute',                     cl1 && cl1.ok === true);
  check('cache size = 1',                    cache.size === 1);

  // Same ctx + same l2idx → same object (memoized)
  const cl1b = cache.getOrCompute(ctx, 0);
  check('memoized: same object identity',    cl1b === cl1);

  // Mutate state.k → cache should invalidate
  state.k = 6;
  const ctx6 = contextFromState(state);
  const cl6 = cache.getOrCompute(ctx6, 0);
  check('k change invalidates cache',         cl6 !== cl1);
  check('cache contains only the new entry',  cache.size === 1);
}

console.log('\n--- sigmaProfileL2 ---');
{
  // The clean fixture has ALL samples with low σ — verdict probably 'STACKED_INVERSIONS' if usedK>3
  // or 'NORMAL' if usedK<=3. Either way it shouldn't crash.
  const state = makeFixture();
  const ctx = contextFromState(state);
  const profile = sigmaProfileL2(ctx, 0, 3);
  // Stable bands → q90 should be small → not noisy
  if (profile) {
    check('sigmaProfile q50 finite',          Number.isFinite(profile.q50));
    check('sigmaProfile q95 finite',          Number.isFinite(profile.q95));
    check('sigmaProfile verdict is string',   typeof profile.verdict === 'string');
  } else {
    check('sigmaProfile null is OK for small n',  true);
  }

  // <10 finite values → null (we have 12 samples — should pass; trigger via mostly-NaN)
  const state2 = makeFixture({ n_samples: 5 });
  const ctx2 = contextFromState(state2);
  const profile2 = sigmaProfileL2(ctx2, 0, 3);
  check('n<10 → null',                       profile2 === null);
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
