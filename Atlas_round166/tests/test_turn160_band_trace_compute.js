// =============================================================================
// turn 160 — Band-trace inheritance map (compute layer)
// =============================================================================
// SPEC_distant_band_concordance_fish_trajectory.md Slice 4, compute half.
//
// What this turn ships:
//   - _bandTraceShannonEntropy(fractions, K) — normalized Shannon entropy
//     in [0, 1]. 0 = degenerate (all in one band), 1 = uniform across K.
//   - _bandTraceForFishSet(fishSet, opts) — per-L2 trace of which bands a
//     fish-set occupies, with regime classification (co_seg / partial /
//     fanned / sparse / no_valid).
//   - _bandTraceRegimeRuns(trace, opts) — runs of consecutive co-seg L2s
//     with stable dominant band; runs are inversion-footprint candidates.
//
// What this turn does NOT ship (deferred to next turn):
//   - UI: trace strip, band-set picker, combinatorial enumeration. The
//     compute layer is callable from the console / programmatically.
//
// Tests:
//   1. _bandTraceShannonEntropy basics (degenerate, uniform, partial)
//   2. _bandTraceForFishSet on a synthetic two-inversion planted layout
//   3. Regime classification thresholds
//   4. Empty / null inputs handled
//   5. _bandTraceRegimeRuns extracts contiguous co-seg stretches
//   6. Runs respect chain breaks (no run spans a chain break)
//   7. min_run_length filter drops single-L2 blips
//   8. allow_partial flag controls inclusion of partial regimes
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
function near(a, b, eps) { return Math.abs(a - b) < (eps || 1e-6); }

function pullFunction(src, name) {
  const decl = `function ${name}`;
  const i = src.indexOf(decl);
  if (i < 0) return null;
  let p = src.indexOf('(', i);
  if (p < 0) return null;
  let braceStart = src.indexOf('{', p);
  if (braceStart < 0) return null;
  let depth = 1, j = braceStart + 1;
  while (j < src.length && depth > 0) {
    const ch = src[j];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    j++;
  }
  return src.substring(i, j);
}

// ============================================================================
// 1. Source-pattern checks
// ============================================================================
console.log('\n=== 1. Source-pattern checks ===');

ok('_bandTraceShannonEntropy defined',
   /function\s+_bandTraceShannonEntropy\s*\(/.test(html));
ok('_bandTraceForFishSet defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('_bandTraceRegimeRuns defined',
   /function\s+_bandTraceRegimeRuns\s*\(/.test(html));
ok('window._bandTraceForFishSet exported',
   /window\._bandTraceForFishSet\s*=\s*_bandTraceForFishSet/.test(html));
ok('window._bandTraceRegimeRuns exported',
   /window\._bandTraceRegimeRuns\s*=\s*_bandTraceRegimeRuns/.test(html));
ok('window._bandTraceShannonEntropy exported',
   /window\._bandTraceShannonEntropy\s*=\s*_bandTraceShannonEntropy/.test(html));
ok('_BTRACE_COSEG_ENTROPY_MAX constant declared',
   /_BTRACE_COSEG_ENTROPY_MAX\s*=\s*0\.40/.test(html));
ok('_BTRACE_FANNED_ENTROPY_MIN constant declared',
   /_BTRACE_FANNED_ENTROPY_MIN\s*=\s*0\.85/.test(html));

// Build sandbox with helpers + their dependencies
const sbx = { console, state: {}, window: {} };
sbx.window.state = sbx.state;
vm.createContext(sbx);

// Inject Hungarian helpers — the trace needs _hungarianChainProjection,
// which itself needs alignLabels. We pull the real functions for high-fidelity
// behaviour, then stub _LINEAGE_CHAIN_BREAK_AGREEMENT explicitly.
const fnAlignLabels = pullFunction(html, 'alignLabels');
const fnHungarian   = pullFunction(html, '_hungarianChainProjection');
const fnFinalize    = pullFunction(html, '_finalizeChain');
const fnEntropy     = pullFunction(html, '_bandTraceShannonEntropy');
const fnTrace       = pullFunction(html, '_bandTraceForFishSet');
const fnRuns        = pullFunction(html, '_bandTraceRegimeRuns');

ok('pulled alignLabels',                !!fnAlignLabels);
ok('pulled _hungarianChainProjection',  !!fnHungarian);
ok('pulled _finalizeChain',             !!fnFinalize);
ok('pulled _bandTraceShannonEntropy',   !!fnEntropy);
ok('pulled _bandTraceForFishSet',       !!fnTrace);
ok('pulled _bandTraceRegimeRuns',       !!fnRuns);

vm.runInContext('const _LINEAGE_CHAIN_BREAK_AGREEMENT = 0.50;', sbx);
vm.runInContext('const _BTRACE_COSEG_ENTROPY_MAX  = 0.40;', sbx);
vm.runInContext('const _BTRACE_FANNED_ENTROPY_MIN = 0.85;', sbx);
vm.runInContext('const _BTRACE_MIN_VALID_FISH     = 3;',    sbx);
vm.runInContext('const _BTRACE_MIN_RUN_LENGTH     = 2;',    sbx);
vm.runInContext(fnAlignLabels, sbx);
vm.runInContext(fnFinalize,    sbx);
vm.runInContext(fnHungarian,   sbx);
vm.runInContext(fnEntropy,     sbx);
vm.runInContext(fnTrace,       sbx);
vm.runInContext(fnRuns,        sbx);

// ============================================================================
// 2. _bandTraceShannonEntropy basics
// ============================================================================
console.log('\n=== 2. _bandTraceShannonEntropy basics ===');

const e_zero = vm.runInContext('_bandTraceShannonEntropy([1, 0, 0], 3)', sbx);
ok('all-in-one-band → entropy 0', near(e_zero, 0));

const e_uniform = vm.runInContext('_bandTraceShannonEntropy([1/3, 1/3, 1/3], 3)', sbx);
ok('uniform K=3 → entropy 1', near(e_uniform, 1, 1e-5));

const e_two = vm.runInContext('_bandTraceShannonEntropy([0.5, 0.5, 0], 3)', sbx);
ok('half-half over 2 of 3 → 0 < entropy < 1',
   e_two > 0.5 && e_two < 0.95);

const e_empty = vm.runInContext('_bandTraceShannonEntropy([], 3)', sbx);
ok('empty fractions → entropy 0', e_empty === 0);

const e_K1 = vm.runInContext('_bandTraceShannonEntropy([1], 1)', sbx);
ok('K=1 → entropy 0 (no normalization possible)', e_K1 === 0);

// Unnormalized counts also work (function divides by total)
const e_counts = vm.runInContext('_bandTraceShannonEntropy([10, 0, 0], 3)', sbx);
ok('unnormalized counts work too', near(e_counts, 0));

// ============================================================================
// 3. _bandTraceForFishSet — synthetic two-inversion planted layout
// ============================================================================
console.log('\n=== 3. _bandTraceForFishSet behaviour ===');

// Setup: 30 samples, 10 L2s, K=3.
//   Inversion A spans L2s 2–4: fish 0–9 in band 0, 10–19 in band 1, 20–29 in band 2
//   Outside A: fish bands are random-but-balanced (use a seeded pattern that's
//              fanned for the fish-set we're interested in)
//   Inversion B spans L2s 6–8: same fish 0–9 land in band 1 (NOT band 0 — different
//              inversion, but same lineage so same fish co-segregate)
const n_samples = 30;
const n_L2 = 10;
const K = 3;
function buildLabels(l2_idx) {
  const labels = new Int8Array(n_samples);
  if (l2_idx >= 2 && l2_idx <= 4) {
    // Inversion A footprint: clean partition 0/1/2
    for (let s = 0; s < n_samples; s++) labels[s] = Math.floor(s / 10);
  } else if (l2_idx >= 6 && l2_idx <= 8) {
    // Inversion B footprint: fish 0-9 → band 1, fish 10-19 → band 2, fish 20-29 → band 0
    for (let s = 0; s < n_samples; s++) {
      const grp = Math.floor(s / 10);
      labels[s] = (grp + 1) % 3;
    }
  } else {
    // Fanned / random: use s % K so fish 0-9 spread across all 3 bands
    for (let s = 0; s < n_samples; s++) labels[s] = s % K;
  }
  return labels;
}

// Inject the labels-getter into the sandbox
sbx._labelsByL2 = {};
for (let i = 0; i < n_L2; i++) sbx._labelsByL2[i] = buildLabels(i);
vm.runInContext('function _testGetLabelsForL2(idx) { return _labelsByL2[idx]; }', sbx);

// Track fish 0–9 (band 0 at L2 #2). Should co-seg in inv A range, fan elsewhere.
sbx._fishSet = [];
for (let s = 0; s < 10; s++) sbx._fishSet.push(s);

const traceJSON = vm.runInContext(`
  (function(){
    const tr = _bandTraceForFishSet(_fishSet, {
      l2_indices: [0,1,2,3,4,5,6,7,8,9],
      K: 3,
      getLabelsForL2: _testGetLabelsForL2,
    });
    return JSON.stringify({
      n_fish_selected: tr.n_fish_selected,
      n_chains: tr.n_chains,
      n_total_L2: tr.n_total_L2,
      K: tr.K,
      per_l2: tr.per_l2.map(e => ({
        l2_idx: e.l2_idx, chain_idx: e.chain_idx, n_valid: e.n_valid,
        band_counts: Array.from(e.band_counts),
        band_fractions: Array.from(e.band_fractions),
        entropy: e.entropy, dominant_band: e.dominant_band,
        dominant_fraction: e.dominant_fraction, regime: e.regime,
      })),
    });
  })();
`, sbx);
const trace = JSON.parse(traceJSON);

ok('trace returned', !!trace);
ok('n_fish_selected = 10', trace && trace.n_fish_selected === 10);
ok('K = 3', trace && trace.K === 3);
ok('per_l2 has 10 entries', trace && trace.per_l2.length === 10);

// At L2 #2: all 10 fish in band 0 → entropy 0 → co_seg
const e2 = trace.per_l2.find(e => e.l2_idx === 2);
ok('L2 #2: dominant_band = 0',  e2 && e2.dominant_band === 0);
ok('L2 #2: dominant_fraction = 1.0', e2 && near(e2.dominant_fraction, 1));
ok('L2 #2: regime = co_seg',    e2 && e2.regime === 'co_seg');
ok('L2 #2: entropy ≈ 0',        e2 && near(e2.entropy, 0));

// At L2 #7: same fish concentrate in band 1 → still co_seg, different band
const e7 = trace.per_l2.find(e => e.l2_idx === 7);
ok('L2 #7: dominant_band = 1',  e7 && e7.dominant_band === 1);
ok('L2 #7: regime = co_seg',    e7 && e7.regime === 'co_seg');

// At L2 #5 (between inversions): fanned — fish 0-9 spread by s%K, ~10/3 in each
const e5 = trace.per_l2.find(e => e.l2_idx === 5);
ok('L2 #5: regime is fanned',
   e5 && e5.regime === 'fanned',
   'got regime=' + (e5 ? e5.regime : 'NULL') + ' entropy=' + (e5 ? e5.entropy : 'NULL'));
ok('L2 #5: entropy > 0.85',     e5 && e5.entropy > 0.85);

// ============================================================================
// 4. Regime classification thresholds — synthetic distributions
// ============================================================================
console.log('\n=== 4. Regime thresholds ===');

// Force a "partial" by giving 70/30 split → entropy ≈ 0.61
sbx._labelsByL2_partial = {};
sbx._labelsByL2_partial[0] = new Int8Array(10);
for (let s = 0; s < 7; s++) sbx._labelsByL2_partial[0][s] = 0;
for (let s = 7; s < 10; s++) sbx._labelsByL2_partial[0][s] = 1;
vm.runInContext('function _getPartialLabels(idx){return _labelsByL2_partial[idx];}', sbx);
vm.runInContext('_partialFishSet = [0,1,2,3,4,5,6,7,8,9];', sbx);

const partialTraceJSON = vm.runInContext(`
  (function(){
    const tr = _bandTraceForFishSet(_partialFishSet, {
      l2_indices: [0], K: 3, getLabelsForL2: _getPartialLabels,
    });
    return JSON.stringify({entropy: tr.per_l2[0].entropy, regime: tr.per_l2[0].regime});
  })();
`, sbx);
const pt = JSON.parse(partialTraceJSON);
ok('70/30 split → regime "partial"', pt.regime === 'partial');
ok('70/30 split → 0.40 < entropy < 0.85',
   pt.entropy > 0.40 && pt.entropy < 0.85,
   'got ' + pt.entropy);

// Test sparse: only 2 valid fish → below min_valid=3
sbx._labelsByL2_sparse = {};
sbx._labelsByL2_sparse[0] = new Int8Array(10);
for (let s = 0; s < 10; s++) sbx._labelsByL2_sparse[0][s] = -1;   // most invalid
sbx._labelsByL2_sparse[0][3] = 0;
sbx._labelsByL2_sparse[0][4] = 0;
vm.runInContext('function _getSparseLabels(idx){return _labelsByL2_sparse[idx];}', sbx);
const sparseJSON = vm.runInContext(`
  (function(){
    const tr = _bandTraceForFishSet([0,1,2,3,4,5,6,7,8,9], {
      l2_indices: [0], K: 3, getLabelsForL2: _getSparseLabels,
    });
    return JSON.stringify({n_valid: tr.per_l2[0].n_valid, regime: tr.per_l2[0].regime});
  })();
`, sbx);
const sp = JSON.parse(sparseJSON);
ok('only 2 valid → regime "sparse"', sp.regime === 'sparse' && sp.n_valid === 2);

// no_valid case: all fish have label -1
sbx._labelsByL2_zero = {};
sbx._labelsByL2_zero[0] = new Int8Array(10);
for (let s = 0; s < 10; s++) sbx._labelsByL2_zero[0][s] = -1;
vm.runInContext('function _getZeroLabels(idx){return _labelsByL2_zero[idx];}', sbx);
const zeroJSON = vm.runInContext(`
  (function(){
    const tr = _bandTraceForFishSet([0,1,2,3,4,5,6,7,8,9], {
      l2_indices: [0], K: 3, getLabelsForL2: _getZeroLabels,
    });
    return JSON.stringify({n_valid: tr.per_l2[0].n_valid, regime: tr.per_l2[0].regime});
  })();
`, sbx);
const zr = JSON.parse(zeroJSON);
ok('all invalid → regime "no_valid"', zr.regime === 'no_valid' && zr.n_valid === 0);

// ============================================================================
// 5. Empty / null inputs
// ============================================================================
console.log('\n=== 5. Empty / null inputs ===');

const nullTrace = vm.runInContext('_bandTraceForFishSet(null, {})', sbx);
ok('null fishSet → null', nullTrace === null);

const emptyTrace = vm.runInContext('_bandTraceForFishSet([], {})', sbx);
ok('empty fishSet → null', emptyTrace === null);

const noL2 = vm.runInContext(`
  _bandTraceForFishSet([0,1,2], {
    l2_indices: [], K: 3, getLabelsForL2: _testGetLabelsForL2,
  });
`, sbx);
ok('empty l2_indices → null', noL2 === null);

// ============================================================================
// 6. _bandTraceRegimeRuns — contiguous co-seg stretches
// ============================================================================
console.log('\n=== 6. _bandTraceRegimeRuns extraction ===');

const runsJSON = vm.runInContext(`
  (function(){
    const tr = _bandTraceForFishSet(_fishSet, {
      l2_indices: [0,1,2,3,4,5,6,7,8,9], K: 3,
      getLabelsForL2: _testGetLabelsForL2,
    });
    return JSON.stringify(_bandTraceRegimeRuns(tr));
  })();
`, sbx);
const runs = JSON.parse(runsJSON);

ok('two runs detected (inv A + inv B)', runs.length === 2);
ok('run[0] starts at L2 #2', runs[0] && runs[0].start_l2_idx === 2);
ok('run[0] ends at L2 #4',   runs[0] && runs[0].end_l2_idx === 4);
ok('run[0] dominant_band = 0', runs[0] && runs[0].dominant_band === 0);
ok('run[0] n_L2 = 3',          runs[0] && runs[0].n_L2 === 3);

ok('run[1] starts at L2 #6', runs[1] && runs[1].start_l2_idx === 6);
ok('run[1] ends at L2 #8',   runs[1] && runs[1].end_l2_idx === 8);
ok('run[1] dominant_band = 1', runs[1] && runs[1].dominant_band === 1);

// ============================================================================
// 7. Runs respect chain breaks
// ============================================================================
console.log('\n=== 7. Chain breaks split runs ===');

// Build a trace with a chain break in the middle of a co-seg stretch.
// Easiest approach: hand-craft a per_l2 array and feed it to _bandTraceRegimeRuns
// directly (the function takes any object with .per_l2).
const handcrafted = {
  K: 3,
  per_l2: [
    { l2_idx: 0, chain_idx: 0, chain_position: 0, n_valid: 10, entropy: 0,
      dominant_band: 0, dominant_fraction: 1, regime: 'co_seg' },
    { l2_idx: 1, chain_idx: 0, chain_position: 1, n_valid: 10, entropy: 0,
      dominant_band: 0, dominant_fraction: 1, regime: 'co_seg' },
    // Chain break between L2 #1 and L2 #2 — even though regime stays co_seg,
    // the chain_idx changes so the run must split.
    { l2_idx: 2, chain_idx: 1, chain_position: 0, n_valid: 10, entropy: 0,
      dominant_band: 0, dominant_fraction: 1, regime: 'co_seg' },
    { l2_idx: 3, chain_idx: 1, chain_position: 1, n_valid: 10, entropy: 0,
      dominant_band: 0, dominant_fraction: 1, regime: 'co_seg' },
  ],
};
sbx._handcraftedTrace = handcrafted;
const splitRunsJSON = vm.runInContext('JSON.stringify(_bandTraceRegimeRuns(_handcraftedTrace))', sbx);
const splitRuns = JSON.parse(splitRunsJSON);
ok('chain break splits run into 2', splitRuns.length === 2);
ok('split run[0] all in chain 0', splitRuns[0] && splitRuns[0].chain_idx === 0);
ok('split run[1] all in chain 1', splitRuns[1] && splitRuns[1].chain_idx === 1);

// ============================================================================
// 8. min_run_length filter + allow_partial flag
// ============================================================================
console.log('\n=== 8. Filters ===');

// Single co-seg L2 surrounded by fanned should be filtered out at min_run_length=2
const blip = {
  K: 3,
  per_l2: [
    { l2_idx: 0, chain_idx: 0, chain_position: 0, n_valid: 10, entropy: 0.9,
      dominant_band: 0, dominant_fraction: 0.4, regime: 'fanned' },
    { l2_idx: 1, chain_idx: 0, chain_position: 1, n_valid: 10, entropy: 0,
      dominant_band: 0, dominant_fraction: 1, regime: 'co_seg' },
    { l2_idx: 2, chain_idx: 0, chain_position: 2, n_valid: 10, entropy: 0.9,
      dominant_band: 1, dominant_fraction: 0.4, regime: 'fanned' },
  ],
};
sbx._blip = blip;
const blipDefaultJSON = vm.runInContext('JSON.stringify(_bandTraceRegimeRuns(_blip))', sbx);
const blipDefault = JSON.parse(blipDefaultJSON);
ok('single-L2 blip filtered at default min_run_length', blipDefault.length === 0);

const blipKeepJSON = vm.runInContext(
  'JSON.stringify(_bandTraceRegimeRuns(_blip, {min_run_length: 1}))', sbx);
const blipKeep = JSON.parse(blipKeepJSON);
ok('single-L2 blip kept at min_run_length=1', blipKeep.length === 1);

// allow_partial=false drops partial-only runs
const partialOnly = {
  K: 3,
  per_l2: [
    { l2_idx: 0, chain_idx: 0, chain_position: 0, n_valid: 10, entropy: 0.6,
      dominant_band: 0, dominant_fraction: 0.7, regime: 'partial' },
    { l2_idx: 1, chain_idx: 0, chain_position: 1, n_valid: 10, entropy: 0.6,
      dominant_band: 0, dominant_fraction: 0.7, regime: 'partial' },
  ],
};
sbx._partialOnly = partialOnly;
const partialIncludedJSON = vm.runInContext(
  'JSON.stringify(_bandTraceRegimeRuns(_partialOnly))', sbx);
const partialIncluded = JSON.parse(partialIncludedJSON);
ok('default (allow_partial=true): partial run included', partialIncluded.length === 1);

const partialExcludedJSON = vm.runInContext(
  'JSON.stringify(_bandTraceRegimeRuns(_partialOnly, {allow_partial: false}))', sbx);
const partialExcluded = JSON.parse(partialExcludedJSON);
ok('allow_partial=false: partial run excluded', partialExcluded.length === 0);

// ============================================================================
// 9. Regression — turn 159 helpers still wired
// ============================================================================
console.log('\n=== 9. Regression checks ===');

ok('turn 159 _gpKaryoProposeAll still defined',
   /function\s+_gpKaryoProposeAll\s*\(/.test(html));
ok('turn 159 window export still in place',
   /window\._gpKaryoProposeAll\s*=\s*_gpKaryoProposeAll/.test(html));
ok('turn 156 V-shape diagnostic still present',
   /function\s+openVShapePlot\s*\(/.test(html));
ok('turn 130 lineage compute still present',
   /function\s+runLineageCompute\s*\(/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
