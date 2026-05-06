// =============================================================================
// turn 130 test — Fish-trajectory lineage compute (Slice 1 of
// SPEC_distant_band_concordance_fish_trajectory.md)
//
// Recovered from old chat (f74cf5d4) where it was proposed but never built.
// Implements Approach 2: pairwise fish concordance over Hungarian-projected
// L2-band trajectories → hierarchical clustering → per-fish lineage labels.
//
// Tests:
//   1.  All exports present + extractable.
//   2.  Hungarian chain projection: simple identity case (no Hungarian
//       needed — labels already aligned).
//   3.  Hungarian chain projection: Hungarian needed (labels swapped).
//   4.  Hungarian chain projection: chain break when concord drops.
//   5.  Concordance matrix: diagonal = 1, symmetric, range [0, 1].
//   6.  Concordance matrix: simple all-same trajectory case → all 1.
//   7.  Concordance matrix: two-lineage planted case → block structure.
//   8.  Concordance matrix: -1 no-call labels excluded from valid count.
//   9.  Lineage clustering: planted 2-lineage → 2 lineages recovered.
//   10. Lineage clustering: all-identical trajectories → 1 lineage.
//   11. Lineage clustering: every-fish-unique → n_samples lineages.
//   12. runLineageCompute: full pipeline + cache hit returns same result.
//   13. runLineageCompute: bails out when fewer than 3 L2s.
//   14. invalidateLineageCache: clears state slot.
//   15. Cache key: includes chrom + mode + threshold.
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

// =============================================================================
// Source-level checks
// =============================================================================
console.log('\n=== Source-level checks ===');

ok('runLineageCompute defined',
   /function runLineageCompute\(opts\)/.test(html));
ok('_hungarianChainProjection defined',
   /function _hungarianChainProjection\(l2_indices, K, getLabelsForL2\)/.test(html));
ok('_concordanceMatrix defined',
   /function _concordanceMatrix\(projection\)/.test(html));
ok('_lineageClustering defined',
   /function _lineageClustering\(concordanceMatrix, n_samples, threshold\)/.test(html));
ok('_lineageCacheKey defined',
   /function _lineageCacheKey\(l2_indices, K, threshold, mode, chrom\)/.test(html));
ok('invalidateLineageCache defined',
   /function invalidateLineageCache\(\)/.test(html));

ok('runLineageCompute exposed on window',
   /window\.runLineageCompute\s*=\s*runLineageCompute/.test(html));
ok('_concordanceMatrix exposed on window',
   /window\._concordanceMatrix\s*=\s*_concordanceMatrix/.test(html));
ok('_LINEAGE_DEFAULT_THRESHOLD constant defined',
   /const _LINEAGE_DEFAULT_THRESHOLD = 0\.50/.test(html));
ok('_LINEAGE_CHAIN_BREAK_AGREEMENT constant defined',
   /const _LINEAGE_CHAIN_BREAK_AGREEMENT = 0\.50/.test(html));

ok('lineage compute reuses _agglomerativeAverageLinkage',
   /_lineageClustering[\s\S]{0,800}_agglomerativeAverageLinkage/.test(html));
ok('lineage compute reuses _cutDendrogram',
   /_lineageClustering[\s\S]{0,800}_cutDendrogram/.test(html));

// =============================================================================
// Behavioural tests (sandboxed)
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        i++;
      }
    } else if (ch === '/' && src[i+1] === '/') {
      while (i < src.length && src[i] !== '\n') i++;
    } else if (ch === '/' && src[i+1] === '*') {
      i += 2;
      while (i < src.length - 1 && !(src[i] === '*' && src[i+1] === '/')) i++;
      i++;
    }
    i++;
  }
  return src.substring(start, i);
}

// Pull all the helpers we need, plus the existing alignLabels + agglomerative
// + cut so the sandbox can run end-to-end.
const sources = {
  alignLabels:                pullFunction(html, 'alignLabels'),
  permutations:               pullFunction(html, 'permutations'),
  _hungarianChainProjection:  pullFunction(html, '_hungarianChainProjection'),
  _finalizeChain:             pullFunction(html, '_finalizeChain'),
  _concordanceMatrix:         pullFunction(html, '_concordanceMatrix'),
  _agglomerativeAverageLinkage: pullFunction(html, '_agglomerativeAverageLinkage'),
  _cutDendrogram:             pullFunction(html, '_cutDendrogram'),
  _lineageClustering:         pullFunction(html, '_lineageClustering'),
  _lineageCacheKey:           pullFunction(html, '_lineageCacheKey'),
  runLineageCompute:          pullFunction(html, 'runLineageCompute'),
  invalidateLineageCache:     pullFunction(html, 'invalidateLineageCache'),
};
for (const k of Object.keys(sources)) {
  ok(`${k} extractable`, !!sources[k]);
}

// Constants — pull as raw text since we need them in the sandbox.
const constsBlock = (html.match(/const _LINEAGE_DEFAULT_THRESHOLD = 0\.50;[\s\S]*?const _LINEAGE_MIN_FISH_PER_LINEAGE = 1;/) || [''])[0];
ok('lineage constants block extractable', !!constsBlock);

function makeSandbox(stateOverride) {
  const sandbox = {
    state: stateOverride || {},
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(constsBlock, sandbox);
  vm.runInContext(sources.permutations, sandbox);
  vm.runInContext(sources.alignLabels, sandbox);
  vm.runInContext(sources._finalizeChain, sandbox);
  vm.runInContext(sources._hungarianChainProjection, sandbox);
  vm.runInContext(sources._concordanceMatrix, sandbox);
  vm.runInContext(sources._agglomerativeAverageLinkage, sandbox);
  vm.runInContext(sources._cutDendrogram, sandbox);
  vm.runInContext(sources._lineageClustering, sandbox);
  vm.runInContext(sources._lineageCacheKey, sandbox);
  vm.runInContext(sources.runLineageCompute, sandbox);
  vm.runInContext(sources.invalidateLineageCache, sandbox);
  return sandbox;
}

// --- Test 2: Hungarian chain projection — identity case
console.log('\nTest 2: Hungarian chain projection (identity)');
{
  // 3 L2s, 4 fish, K=3. Labels already aligned.
  // Fish 0: always band 0
  // Fish 1: always band 0
  // Fish 2: always band 1
  // Fish 3: always band 2
  const labelsPerL2 = [
    new Int8Array([0, 0, 1, 2]),
    new Int8Array([0, 0, 1, 2]),
    new Int8Array([0, 0, 1, 2]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var getLabels = function(i) { return labelsPerL2[i]; };
    var proj = _hungarianChainProjection([0, 1, 2], 3, getLabels);
  `, sandbox);
  ok('1 chain when all L2s align',     sandbox.proj.n_chains === 1);
  ok('all 3 L2s in single chain',      sandbox.proj.chains[0].n_L2 === 3);
  ok('n_samples = 4',                  sandbox.proj.n_samples === 4);
  ok('n_total_L2 = 3',                 sandbox.proj.n_total_L2 === 3);
}

// --- Test 3: Hungarian chain projection — labels swapped (Hungarian needed)
console.log('\nTest 3: Hungarian chain projection (Hungarian needed)');
{
  // L2 #0: fish [A=0, B=0, C=1, D=2]   (band 0 is the AB fish, band 1 is C, band 2 is D)
  // L2 #1: same fish but K-means labelled them as 2/2/0/1 (cyclic permutation).
  //   After Hungarian alignment to L2 #0, this should map back to 0/0/1/2.
  const labelsPerL2 = [
    new Int8Array([0, 0, 1, 2]),
    new Int8Array([2, 2, 0, 1]),   // perm (0→2, 1→0, 2→1) of L2 #0
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var getLabels = function(i) { return labelsPerL2[i]; };
    var proj = _hungarianChainProjection([0, 1], 3, getLabels);
  `, sandbox);
  ok('1 chain (Hungarian aligns both L2s)',
     sandbox.proj.n_chains === 1, 'got ' + sandbox.proj.n_chains);
  // After alignment, projected[1] should equal projected[0] (labels match perfectly).
  ok('projected L2 #1 row matches L2 #0 row after alignment',
     sandbox.proj.chains[0].projected[4] === 0 &&    // row 1 fish 0
     sandbox.proj.chains[0].projected[5] === 0 &&    // row 1 fish 1
     sandbox.proj.chains[0].projected[6] === 1 &&    // row 1 fish 2
     sandbox.proj.chains[0].projected[7] === 2);     // row 1 fish 3
}

// --- Test 4: Hungarian chain projection — chain break
console.log('\nTest 4: Hungarian chain projection (chain break)');
{
  // 3 L2s. L2 #0 and L2 #1 align well. L2 #2 has labels totally
  // uncorrelated → max Hungarian agreement < 0.5 → chain breaks.
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2]),
    new Int8Array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2]),
    // Random shuffle that doesn't have a good permutation match.
    new Int8Array([0, 1, 2, 0, 1, 2, 0, 1, 2, 0]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var getLabels = function(i) { return labelsPerL2[i]; };
    var proj = _hungarianChainProjection([0, 1, 2], 3, getLabels);
  `, sandbox);
  ok('2 chains formed (L2 #2 breaks)', sandbox.proj.n_chains === 2);
  ok('first chain has 2 L2s',  sandbox.proj.chains[0].n_L2 === 2);
  ok('second chain has 1 L2',  sandbox.proj.chains[1].n_L2 === 1);
  ok('n_total_L2 still equals 3', sandbox.proj.n_total_L2 === 3);
}

// --- Test 5+6: concordance matrix on simple cases
console.log('\nTest 5: concordance matrix (all-same trajectory)');
{
  // 4 fish, 3 L2s, all identical labels.
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 0]),
    new Int8Array([0, 0, 0, 0]),
    new Int8Array([0, 0, 0, 0]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0, 1, 2], 3, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
  `, sandbox);
  ok('matrix is Float32Array',
     sandbox.c && sandbox.c.constructor.name === 'Float32Array');
  ok('matrix has length n_samples^2 = 16', sandbox.c.length === 16);
  ok('all entries == 1 (everyone always together)',
     Array.from(sandbox.c).every(v => Math.abs(v - 1) < 1e-6));
}

console.log('\nTest 6: concordance matrix (two-lineage planted)');
{
  // 6 fish, 4 L2s. Fish 0,1,2 always go to band 0; Fish 3,4,5 always go to band 1.
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0,1,2,3], 2, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
  `, sandbox);
  // Expect: c(0,1) = c(0,2) = c(1,2) = 1; c(3,4) = c(3,5) = c(4,5) = 1.
  // Cross-block c(0,3) = c(0,4) = c(0,5) = 0 (never in same band).
  ok('within-block concordance == 1 (fish 0,1)',
     Math.abs(sandbox.c[0 * 6 + 1] - 1) < 1e-6,
     'got: ' + sandbox.c[0*6+1]);
  ok('within-block concordance == 1 (fish 3,4)',
     Math.abs(sandbox.c[3 * 6 + 4] - 1) < 1e-6);
  ok('cross-block concordance == 0 (fish 0,3)',
     Math.abs(sandbox.c[0 * 6 + 3] - 0) < 1e-6);
  ok('matrix is symmetric (c[0,3] == c[3,0])',
     Math.abs(sandbox.c[0 * 6 + 3] - sandbox.c[3 * 6 + 0]) < 1e-6);
  ok('diagonal == 1', Math.abs(sandbox.c[0] - 1) < 1e-6);
}

// --- Test 7: -1 no-call labels excluded
console.log('\nTest 7: -1 no-call labels excluded from validCount');
{
  // 3 fish, 4 L2s. Fish 0 and Fish 1 always together; Fish 2 has all -1s.
  const labelsPerL2 = [
    new Int8Array([0, 0, -1]),
    new Int8Array([1, 1, -1]),
    new Int8Array([0, 0, -1]),
    new Int8Array([1, 1, -1]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0,1,2,3], 2, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
  `, sandbox);
  ok('fish 0 and 1 concordance is 1 (always together when both called)',
     Math.abs(sandbox.c[0 * 3 + 1] - 1) < 1e-6);
  ok('fish 0 and 2 concordance is 0 (fish 2 always no-call → no valid pairs)',
     Math.abs(sandbox.c[0 * 3 + 2] - 0) < 1e-6);
}

// --- Test 8: lineage clustering (planted two-lineage)
console.log('\nTest 8: lineage clustering recovers planted lineages');
{
  // Same as test 6 but follow through to lineage labels.
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0,1,2,3], 2, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
    var clu = _lineageClustering(c, 6, 0.5);
  `, sandbox);
  ok('exactly 2 lineages found',
     sandbox.clu.n_lineages === 2,
     'got: ' + sandbox.clu.n_lineages);
  // Fish 0,1,2 in one lineage; Fish 3,4,5 in another.
  const labels = Array.from(sandbox.clu.lineage_id_per_sample);
  ok('fish 0,1,2 share a lineage',
     labels[0] === labels[1] && labels[1] === labels[2]);
  ok('fish 3,4,5 share a lineage',
     labels[3] === labels[4] && labels[4] === labels[5]);
  ok('the two groups are different lineages',
     labels[0] !== labels[3]);
}

// --- Test 9: all-identical → 1 lineage
console.log('\nTest 9: all-identical trajectories → 1 lineage');
{
  const labelsPerL2 = [
    new Int8Array([0, 0, 0]),
    new Int8Array([0, 0, 0]),
    new Int8Array([0, 0, 0]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0,1,2], 2, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
    var clu = _lineageClustering(c, 3, 0.5);
  `, sandbox);
  ok('1 lineage when everyone identical',
     sandbox.clu.n_lineages === 1);
}

// --- Test 10: every fish unique → n_samples lineages (with strict threshold)
console.log('\nTest 10: every fish unique → many lineages (strict cut)');
{
  // 3 fish, 3 L2s, every fish in a different band at every L2.
  // Strict threshold 0.0 → only exact matches merge → no merges → 3 lineages.
  const labelsPerL2 = [
    new Int8Array([0, 1, 2]),
    new Int8Array([0, 1, 2]),
    new Int8Array([0, 1, 2]),
  ];
  const sandbox = makeSandbox({});
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var proj = _hungarianChainProjection([0,1,2], 3, function(i){ return labelsPerL2[i]; });
    var c = _concordanceMatrix(proj);
    var clu = _lineageClustering(c, 3, 0.001);
  `, sandbox);
  ok('3 lineages with strict threshold', sandbox.clu.n_lineages === 3);
}

// --- Test 11: runLineageCompute end-to-end + cache hit
console.log('\nTest 11: runLineageCompute end-to-end + cache hit');
{
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
  ];
  const sandbox = makeSandbox({
    data: {
      chrom: 'LG28',
      l2_envelopes: [{}, {}, {}],
      n_samples: 6,
    },
    k: 2,
    activeMode: 'default',
  });
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    var getLabels = function(i) { return labelsPerL2[i]; };
    var r1 = runLineageCompute({ getLabelsForL2: getLabels });
    var r2 = runLineageCompute({ getLabelsForL2: getLabels });
  `, sandbox);
  ok('first run returns a result', !!sandbox.r1);
  ok('result has n_samples = 6', sandbox.r1.n_samples === 6);
  ok('result has 2 lineages', sandbox.r1.n_lineages === 2);
  ok('result has chrom = LG28', sandbox.r1.chrom === 'LG28');
  ok('cache hit returns same reference (r1 === r2)', sandbox.r1 === sandbox.r2);
  ok('state.lineageResult is set', sandbox.state.lineageResult === sandbox.r1);
  ok('state.lineageCacheKey is set', typeof sandbox.state.lineageCacheKey === 'string');
}

// --- Test 12: bails out when fewer than 3 L2s
console.log('\nTest 12: bails out when fewer than _LINEAGE_MIN_L2_FOR_COMPUTE L2s');
{
  const sandbox = makeSandbox({
    data: {
      chrom: 'LG01',
      l2_envelopes: [{}, {}],   // only 2
      n_samples: 6,
    },
    k: 2,
    activeMode: 'default',
  });
  vm.runInContext(`
    var r = runLineageCompute({ getLabelsForL2: function(){ return new Int8Array([0,0,0,1,1,1]); } });
  `, sandbox);
  ok('returns null with < 3 L2s', sandbox.r === null);
  ok('state.lineageResult is null', sandbox.state.lineageResult === null);
}

// --- Test 13: invalidateLineageCache clears the slot
console.log('\nTest 13: invalidateLineageCache');
{
  const labelsPerL2 = [
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
    new Int8Array([0, 0, 0, 1, 1, 1]),
  ];
  const sandbox = makeSandbox({
    data: {
      chrom: 'LG28',
      l2_envelopes: [{}, {}, {}],
      n_samples: 6,
    },
    k: 2,
    activeMode: 'default',
  });
  sandbox.labelsPerL2 = labelsPerL2;
  vm.runInContext(`
    runLineageCompute({ getLabelsForL2: function(i){ return labelsPerL2[i]; } });
    invalidateLineageCache();
  `, sandbox);
  ok('lineageResult cleared',  sandbox.state.lineageResult === null);
  ok('lineageCacheKey cleared', sandbox.state.lineageCacheKey === null);
}

// --- Test 14: cache key includes chrom + mode + threshold
console.log('\nTest 14: cache key includes chrom + mode + threshold');
{
  const sandbox = makeSandbox({});
  vm.runInContext(`
    var k1 = _lineageCacheKey([0,1,2], 3, 0.5, 'default', 'LG28');
    var k2 = _lineageCacheKey([0,1,2], 3, 0.5, 'default', 'LG12');
    var k3 = _lineageCacheKey([0,1,2], 3, 0.5, 'detailed', 'LG28');
    var k4 = _lineageCacheKey([0,1,2], 3, 0.7, 'default', 'LG28');
    var k5 = _lineageCacheKey([0,1,2], 3, 0.5, 'default', 'LG28');
  `, sandbox);
  ok('different chrom → different key', sandbox.k1 !== sandbox.k2);
  ok('different mode → different key',  sandbox.k1 !== sandbox.k3);
  ok('different threshold → different key', sandbox.k1 !== sandbox.k4);
  ok('same args → same key', sandbox.k1 === sandbox.k5);
}

// =============================================================================
// Final tally
// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
