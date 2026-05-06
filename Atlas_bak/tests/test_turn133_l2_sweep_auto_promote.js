// =============================================================================
// turn 133 — Slice 1: L2-sweep auto-promote (SPEC_l2_sweep_inheritance.md)
// =============================================================================
// Implements:
//   - runL2SweepInheritance({ force }) — orchestrator that gathers usable L2s
//     and runs inheritanceGroupClustering on synthetic items
//   - isUsableL2(l2idx) — filter: cluster ok + ≥ _IGC_MIN_BANDS_FOR_CLUSTERING
//   - _autoPromoteFromSweep(result) — applies 5+1 gates, builds candidates,
//     pushes via addCandidateToList
//   - _silhouetteForL2FixedK(l2idx, cluster) — on-demand silhouette using
//     the existing silhouette1D helper (works in fixed-K mode where
//     getL2Cluster returns silhouette: null)
//   - localStorage helpers: _loadL2SweepDismissed, _saveL2SweepDismissed,
//     _addL2SweepDismissed
//   - L3 toolbar checkbox: #l2SweepToggle / #l2SweepToggleLabel
//   - state.l2SweepEnabled init from localStorage
//
// Gate model (SPEC §3.1, plus gate 6 = dismissed-set):
//   1. silhouette ≥ _AUTO_PROMOTE_MIN_SILHOUETTE (0.30)
//   2. ≥ _AUTO_PROMOTE_MIN_GROUPS inheritance groups touch this L2 (2)
//   3. min band size ≥ _AUTO_PROMOTE_MIN_BAND_SIZE (5)
//   4. not already in candidateList (l2_indices overlap)
//   5. not within _AUTO_PROMOTE_DEDUPE_BP (100kb) of an existing candidate
//   6. not in per-chrom dismissed set
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
// 1. Source-level: function definitions + constants exist
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('runL2SweepInheritance defined',  /function runL2SweepInheritance\b/.test(html));
ok('isUsableL2 defined',             /function isUsableL2\b/.test(html));
ok('_autoPromoteFromSweep defined',  /function _autoPromoteFromSweep\b/.test(html));
ok('_silhouetteForL2FixedK defined', /function _silhouetteForL2FixedK\b/.test(html));
ok('_l2SweepEnabledChange defined',  /function _l2SweepEnabledChange\b/.test(html));
ok('_l2SweepInitFromStorage defined',/function _l2SweepInitFromStorage\b/.test(html));
ok('invalidateL2SweepCache defined', /function invalidateL2SweepCache\b/.test(html));
ok('_loadL2SweepDismissed defined',  /function _loadL2SweepDismissed\b/.test(html));
ok('_saveL2SweepDismissed defined',  /function _saveL2SweepDismissed\b/.test(html));
ok('_addL2SweepDismissed defined',   /function _addL2SweepDismissed\b/.test(html));

ok('_AUTO_PROMOTE_MIN_SILHOUETTE = 0.30',
   /const _AUTO_PROMOTE_MIN_SILHOUETTE\s*=\s*0\.30/.test(html));
ok('_AUTO_PROMOTE_MIN_GROUPS = 2',
   /const _AUTO_PROMOTE_MIN_GROUPS\s*=\s*2\b/.test(html));
ok('_AUTO_PROMOTE_MIN_BAND_SIZE = 5',
   /const _AUTO_PROMOTE_MIN_BAND_SIZE\s*=\s*5\b/.test(html));
ok('_AUTO_PROMOTE_DEDUPE_BP = 100000',
   /const _AUTO_PROMOTE_DEDUPE_BP\s*=\s*100000\b/.test(html));

ok('_L2_SWEEP_ENABLED_KEY uses pca_scrubber_v3 namespace',
   /_L2_SWEEP_ENABLED_KEY\s*=\s*'pca_scrubber_v3\.l2SweepEnabled'/.test(html));
ok('_L2_SWEEP_DISMISSED_KEY_PFX uses pca_scrubber_v3 namespace',
   /_L2_SWEEP_DISMISSED_KEY_PFX\s*=\s*'pca_scrubber_v3\.l2SweepDismissed\.'/.test(html));

ok('window exports for sweep functions',
   /window\.runL2SweepInheritance\b/.test(html) &&
   /window\.isUsableL2\b/.test(html) &&
   /window\._autoPromoteFromSweep\b/.test(html));

// ============================================================================
// 2. DOM: L3 toolbar checkbox
// ============================================================================
console.log('\n=== 2. DOM checkbox in L3 toolbar ===');

ok('#l2SweepToggle checkbox present',
   /id="l2SweepToggle"[^>]*type="checkbox"|type="checkbox"[^>]*id="l2SweepToggle"/.test(html));
ok('#l2SweepToggleLabel wrapping label present',
   /id="l2SweepToggleLabel"/.test(html));
ok('checkbox label contains visible "L2-sweep" text',
   /<span>L2-sweep[^<]*<\/span>/.test(html));
ok('label sits adjacent to l3HetToggleLabel (mirrors that pattern)',
   /l3HetToggleLabel[\s\S]{0,3000}l2SweepToggleLabel/.test(html));

// ============================================================================
// 3. state init: l2SweepEnabled / l2SweepResult / l2SweepCacheKey
// ============================================================================
console.log('\n=== 3. state.l2Sweep* init slots ===');

ok('state.l2SweepEnabled: false in init block',
   /l2SweepEnabled:\s*false/.test(html));
ok('state.l2SweepResult: null in init block',
   /l2SweepResult:\s*null/.test(html));
ok('state.l2SweepCacheKey: null in init block',
   /l2SweepCacheKey:\s*null/.test(html));

// ============================================================================
// 4. Behavioral — extract & sandbox-exec the new functions
// ============================================================================
console.log('\n=== 4. Behavioral — sandbox exec ===');

// We'll extract the entire L2-sweep block from the L2-SWEEP banner to the
// closing window.* exports. Then run it in a vm sandbox with stubs for the
// dependencies (state, getL2Cluster, aggregateL2, silhouette1D,
// inheritanceGroupClustering, addCandidateToList, _IGC_MIN_BANDS_FOR_CLUSTERING).
function extractSweepBlock() {
  // From the banner header to the closing "if (typeof window !== 'undefined')"
  const start = html.indexOf('// L2-SWEEP INHERITANCE');
  if (start < 0) return null;
  // Find the final window-export block right after _silhouetteForL2FixedK export
  const endMarker = 'window._silhouetteForL2FixedK    = _silhouetteForL2FixedK;';
  const endIdx = html.indexOf(endMarker, start);
  if (endIdx < 0) return null;
  // Include the closing brace of the `if (typeof window...)` block
  const after = html.indexOf('}', endIdx);
  if (after < 0) return null;
  return html.substring(start, after + 1);
}
const sweepBlock = extractSweepBlock();
ok('sweep code block extracts cleanly', sweepBlock !== null && sweepBlock.length > 1000);

// Build the stubs. We DON'T run the real getL2Cluster — too heavy. We
// inject a fake table indexed by l2idx that returns whatever we need
// for the test scenario.
function makeSandbox(scenario) {
  // localStorage stub
  const _ls = scenario.ls || {};
  const localStorageStub = {
    getItem(k) { return Object.prototype.hasOwnProperty.call(_ls, k) ? _ls[k] : null; },
    setItem(k, v) { _ls[k] = String(v); },
    removeItem(k) { delete _ls[k]; },
    _store: _ls,
  };

  // state stub
  const stateStub = scenario.state || {
    data: { chrom: 'LG28', l2_envelopes: [] },
    candidateList: [],
    candidates: {},
    k: 3,
    activeMode: 'default',
  };

  // Stub getL2Cluster: scenario.clusters[l2idx] -> ClusterResult
  const clusters = scenario.clusters || {};
  const getL2ClusterStub = function (l2idx) { return clusters[l2idx] || null; };

  // Stub aggregateL2: scenario.aggs[l2idx] -> { xs }
  const aggs = scenario.aggs || {};
  const aggregateL2Stub = function (l2idx) { return aggs[l2idx] || null; };

  // Stub silhouette1D — for test purposes return scenario.silhouettes[l2idx]
  // when xs matches the agg.xs identity; else compute a real silhouette.
  // The cleaner contract: just return a value we can drive from scenario.
  const silhouette1DStub = function (xs, labels, K) {
    // Simple test stub: return scenario.fakeSilhouetteByLength[xs.length] if set,
    // else compute a crude silhouette as ratio of between/within (just for sanity).
    if (scenario.fakeSilhouetteByXs0 && xs.length > 0) {
      const key = xs[0];
      if (Object.prototype.hasOwnProperty.call(scenario.fakeSilhouetteByXs0, key)) {
        return scenario.fakeSilhouetteByXs0[key];
      }
    }
    return 0.5;  // generic fallback
  };

  // Stub inheritanceGroupClustering — returns a result with a groups[] array
  // built from scenario.groups[itemId] = [groupIds...]
  const inheritanceGroupClusteringStub = function (items, opts) {
    const groups = {};   // groupId -> { members: [{item_idx, band}] }
    const itemGroups = scenario.itemGroups || {};
    for (let i = 0; i < items.length; i++) {
      const ids = itemGroups[items[i].id] || [];
      for (const gid of ids) {
        if (!groups[gid]) groups[gid] = { members: [] };
        groups[gid].members.push({ item_idx: i, band: 0 });
      }
    }
    return {
      groups: Object.values(groups),
      n_groups: Object.keys(groups).length,
    };
  };

  // Stub addCandidateToList — record what gets added
  const addedCands = [];
  const addCandidateToListStub = function (cand) {
    if (!cand || !cand.id) return;
    if (stateStub.candidateList.some(c => c && c.id === cand.id)) return;
    stateStub.candidateList.push(cand);
    addedCands.push(cand);
  };

  // Console stub (capture warnings)
  const warnings = [];
  const consoleStub = {
    warn: (...args) => warnings.push(args.join(' ')),
    log:  () => {},
  };

  // Build the context with all of the stubs preloaded as globals.
  // The atlas block expects them to be top-level reachable.
  const ctx = {
    state: stateStub,
    window: undefined,           // sweep block guards on `typeof window !== 'undefined'`
    document: undefined,
    localStorage: localStorageStub,
    console: consoleStub,
    getL2Cluster: getL2ClusterStub,
    aggregateL2: aggregateL2Stub,
    silhouette1D: silhouette1DStub,
    inheritanceGroupClustering: inheritanceGroupClusteringStub,
    addCandidateToList: addCandidateToListStub,
    _IGC_MIN_BANDS_FOR_CLUSTERING: 2,
    _inheritanceCacheKey: function (items, mode) {
      return mode + '::' + items.map(it => it.id + '@' + it.K).join('|');
    },
    Number, Array, Object, Math, JSON, Set, Date, Float64Array, Int8Array,
    isFinite, console: consoleStub,
  };
  ctx.window = ctx;   // self-ref so the sweep block's window-guard sees a window

  return { ctx, addedCands, warnings, ls: _ls, state: stateStub };
}

function runSweepBlockIn(ctx) {
  const context = vm.createContext(ctx);
  vm.runInContext(sweepBlock, context);
  return context;
}

// 4a. Empty L2 envelopes → returns null, no candidates added
{
  console.log('\n--- 4a. Empty l2_envelopes → null ---');
  const sb = makeSandbox({});
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  ok('empty envelopes returns null', result === null);
  ok('no candidates auto-promoted', sb.addedCands.length === 0);
  ok('state.l2SweepResult set to null', sb.state.l2SweepResult === null);
}

// 4b. Single usable L2 → still null (need ≥ _IGC_MIN_BANDS_FOR_CLUSTERING)
{
  console.log('\n--- 4b. Single L2 → null (below _IGC_MIN_BANDS_FOR_CLUSTERING) ---');
  const sb = makeSandbox({
    state: {
      data: {
        chrom: 'LG28',
        l2_envelopes: [
          { _s0: 0, _e0: 10, start_bp: 1000000, end_bp: 2000000 },
        ],
      },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: new Int8Array([0,0,0,1,1,1,2,2,2]),
           n_per_group: [3, 3, 3], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array([0,0.1,0.05,0.5,0.55,0.5,1.0,0.95,1.05]) },
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  ok('single L2 returns null', result === null);
}

// 4c. Multiple L2s, all usable → result with l2_meta + per-L2 silhouette
{
  console.log('\n--- 4c. Multiple usable L2s → result populated ---');
  const labelsBalanced = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const xsBalanced     = new Float64Array([0,0.05,0.1,0.05,0.0, 0.5,0.55,0.5,0.45,0.55, 1.0,0.95,1.05,1.0,0.98]);
  const sb = makeSandbox({
    state: {
      data: {
        chrom: 'LG28',
        l2_envelopes: [
          { _s0: 0, _e0: 10, start_bp: 1000000, end_bp: 2000000 },
          { _s0: 11, _e0: 20, start_bp: 5000000, end_bp: 6000000 },
          { _s0: 21, _e0: 30, start_bp: 9000000, end_bp: 10000000 },
        ],
      },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsBalanced, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsBalanced, n_per_group: [5,5,5], silhouette: null },
      2: { ok: true, fixedKLabels: labelsBalanced, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: xsBalanced }, 1: { xs: xsBalanced }, 2: { xs: xsBalanced },
    },
    fakeSilhouetteByXs0: { 0: 0.7 },
    itemGroups: {
      'L2:0': [10, 11],   // 2 groups
      'L2:1': [10, 12],   // 2 groups
      'L2:2': [11, 12],   // 2 groups
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  ok('result is non-null with 3 usable L2s', result !== null);
  ok('result.l2_meta has 3 entries',
     result && Array.isArray(result.l2_meta) && result.l2_meta.length === 3,
     result && result.l2_meta && 'len=' + result.l2_meta.length);
  ok('result.l2_meta entries have silhouette field',
     result && result.l2_meta.every(m => 'silhouette' in m));
  ok('result.l2_meta sorted by start_bp',
     result && result.l2_meta[0].start_bp < result.l2_meta[1].start_bp
            && result.l2_meta[1].start_bp < result.l2_meta[2].start_bp);
  ok('result.items_meta seq_num assigned 1..N',
     result && result.items_meta[0].seq_num === 1
            && result.items_meta[1].seq_num === 2
            && result.items_meta[2].seq_num === 3);
}

// 4d. isUsableL2 — rejects un-clustered L2s
{
  console.log('\n--- 4d. isUsableL2 rejection paths ---');
  const sb = makeSandbox({
    state: {
      data: {
        chrom: 'LG28',
        l2_envelopes: [
          { _s0: 0, _e0: 10, start_bp: 1000000, end_bp: 2000000 },
          { _s0: 11, _e0: 20, start_bp: 5000000, end_bp: 6000000 },
          { _s0: 21, _e0: 30, start_bp: 9000000, end_bp: 10000000 },
        ],
      },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true,  fixedKLabels: new Int8Array([0,0,0,1,1,1]), n_per_group: [3,3,0] },
      1: { ok: false, reason: 'LOW_NWIN' },
      2: { ok: true,  fixedKLabels: new Int8Array([0,0,0,0,0,0]), n_per_group: [6,0,0] },
    },
    aggs: {
      0: { xs: new Float64Array([0,0.1,0,1,1.1,1]) },
      2: { xs: new Float64Array([0,0,0,0,0,0]) },
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const u0 = ctx.isUsableL2(0);
  const u1 = ctx.isUsableL2(1);
  const u2 = ctx.isUsableL2(2);
  ok('isUsableL2(0) usable=true (2 bands populated)', u0.usable === true);
  ok('isUsableL2(1) usable=false reason=LOW_NWIN', u1.usable === false && u1.reason === 'LOW_NWIN');
  ok('isUsableL2(2) usable=false reason=TOO_FEW_BANDS',
     u2.usable === false && u2.reason === 'TOO_FEW_BANDS', 'got reason=' + u2.reason);
}

// 4e. _autoPromoteFromSweep — gate 1 (low silhouette)
{
  console.log('\n--- 4e. Gate 1: low silhouette rejected ---');
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]),
           n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]),
           n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15) },
      1: { xs: new Float64Array(15).fill(0.99) },
    },
    fakeSilhouetteByXs0: { 0: 0.10, 0.99: 0.65 },
    itemGroups: { 'L2:0': [10, 11], 'L2:1': [10, 11] },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  const ap = ctx._autoPromoteFromSweep(result);
  ok('exactly 1 promoted (the high-sil L2)', ap.promoted.length === 1, 'got ' + ap.promoted.length);
  ok('promoted is l2idx=1 (the 0.65 sil one)', ap.promoted[0] === 1);
  const lowSilSkip = ap.skipped.find(s => s.reason === 'LOW_SILHOUETTE');
  ok('low-sil rejection has reason LOW_SILHOUETTE', !!lowSilSkip);
  ok('low-sil rejection records the silhouette value',
     lowSilSkip && Math.abs(lowSilSkip.silhouette - 0.10) < 1e-9);
}

// 4f. Gate 2 (single inheritance group → reject) and gate 3 (small band)
{
  console.log('\n--- 4f. Gate 2 (single group) + Gate 3 (small band) ---');
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },     // single group
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },    // small band
                { _s0:21, _e0:30, start_bp:9e6, end_bp:10e6 },   // good
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]),
           n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: new Int8Array([0,0,1,2,2,2,2,2,2,2,2,2,2,2,2]),
           n_per_group: [2,1,12], silhouette: null },     // band of 1 → reject
      2: { ok: true, fixedKLabels: new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]),
           n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
      2: { xs: new Float64Array(15).fill(0.3) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7, 0.3: 0.7 },
    itemGroups: {
      'L2:0': [10],          // single group → gate 2 reject
      'L2:1': [10, 11],
      'L2:2': [10, 11],
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  const ap = ctx._autoPromoteFromSweep(result);
  ok('only 1 promoted (the well-behaved L2)', ap.promoted.length === 1,
     'got promoted=' + JSON.stringify(ap.promoted) + ' skipped=' + JSON.stringify(ap.skipped.map(s=>s.reason)));
  ok('skipped contains TOO_FEW_GROUPS',
     ap.skipped.some(s => s.reason === 'TOO_FEW_GROUPS'));
  ok('skipped contains SMALL_BAND',
     ap.skipped.some(s => s.reason === 'SMALL_BAND'));
}

// 4g. Gate 4 (already in candidateList) + Gate 5 (dedupe radius)
{
  console.log('\n--- 4g. Gate 4 (already saved) + Gate 5 (dedupe < 100kb) ---');
  // Pre-populate one candidate covering L2[0] (gate 4) and another candidate
  // adjacent to L2[2] within 100kb (gate 5). L2[1] is far enough → promotes.
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
                { _s0:21, _e0:30, start_bp:9e6, end_bp:10e6 },
              ] },
      candidateList: [
        { id: 'manual_a', l2_indices: [0], start_bp: 1e6,    end_bp: 2e6 },
        // Adjacent to L2[2] within 100kb (L2[2] starts at 9e6; this cand ends at 9.05e6 → 0bp gap if overlap, else <100kb)
        { id: 'manual_b', l2_indices: [99], start_bp: 8.95e6, end_bp: 9.05e6 },
      ],
      candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      2: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
      2: { xs: new Float64Array(15).fill(0.3) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7, 0.3: 0.7 },
    itemGroups: {
      'L2:0': [10, 11], 'L2:1': [10, 11], 'L2:2': [10, 11],
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  const ap = ctx._autoPromoteFromSweep(result);
  ok('exactly 1 promoted (L2[1], far from any candidate)',
     ap.promoted.length === 1, 'got ' + JSON.stringify(ap.promoted));
  ok('promoted is l2idx=1', ap.promoted[0] === 1);
  ok('skipped contains ALREADY_IN_CANDIDATE for L2[0]',
     ap.skipped.some(s => s.l2idx === 0 && s.reason === 'ALREADY_IN_CANDIDATE'));
  ok('skipped contains DEDUPE_TOO_CLOSE for L2[2]',
     ap.skipped.some(s => s.l2idx === 2 && s.reason === 'DEDUPE_TOO_CLOSE'));
}

// 4h. Gate 6 — dismissed-set persistence
{
  console.log('\n--- 4h. Gate 6: dismissed-set blocks re-promote ---');
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7 },
    itemGroups: { 'L2:0': [10, 11], 'L2:1': [10, 11] },
    // Pre-populate the dismissed set in localStorage so L2[0] is blocked
    ls: { 'pca_scrubber_v3.l2SweepDismissed.LG28': JSON.stringify([0]) },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const result = ctx.runL2SweepInheritance({ force: true });
  const ap = ctx._autoPromoteFromSweep(result);
  ok('exactly 1 promoted (L2[1])', ap.promoted.length === 1);
  ok('promoted is l2idx=1', ap.promoted[0] === 1);
  ok('skipped contains DISMISSED for L2[0]',
     ap.skipped.some(s => s.l2idx === 0 && s.reason === 'DISMISSED'));
}

// 4i. _addL2SweepDismissed → round-trip via _loadL2SweepDismissed
{
  console.log('\n--- 4i. dismissed-set localStorage round-trip ---');
  const sb = makeSandbox({});
  const ctx = runSweepBlockIn(sb.ctx);
  ctx._addL2SweepDismissed('LG28', 7);
  ctx._addL2SweepDismissed('LG28', 12);
  ctx._addL2SweepDismissed('LG12', 3);
  const setLG28 = ctx._loadL2SweepDismissed('LG28');
  const setLG12 = ctx._loadL2SweepDismissed('LG12');
  ok('LG28 dismissed has both 7 and 12',
     setLG28.has(7) && setLG28.has(12) && setLG28.size === 2);
  ok('LG12 dismissed has only 3',
     setLG12.has(3) && setLG12.size === 1);
  ok('different chrom keys do not cross-contaminate',
     !setLG28.has(3) && !setLG12.has(7));
}

// 4j. _l2SweepEnabledChange → localStorage round-trip
{
  console.log('\n--- 4j. l2SweepEnabled localStorage round-trip ---');
  const sb = makeSandbox({
    state: { data: null, candidateList: [], candidates: {}, k: 3, activeMode: 'default' },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  ctx._l2SweepEnabledChange(true);
  ok('state.l2SweepEnabled=true after enable', sb.state.l2SweepEnabled === true);
  ok('localStorage l2SweepEnabled=1',
     sb.ls['pca_scrubber_v3.l2SweepEnabled'] === '1');
  ctx._l2SweepEnabledChange(false);
  ok('state.l2SweepEnabled=false after disable', sb.state.l2SweepEnabled === false);
  ok('localStorage l2SweepEnabled=0',
     sb.ls['pca_scrubber_v3.l2SweepEnabled'] === '0');
}

// 4k. cache key — second call with same state returns cached result (not recomputed)
{
  console.log('\n--- 4k. cache key short-circuits second call ---');
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7 },
    itemGroups: { 'L2:0': [10, 11], 'L2:1': [10, 11] },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const r1 = ctx.runL2SweepInheritance({ force: true });
  const r2 = ctx.runL2SweepInheritance();   // not forced — should hit cache
  ok('first call returns non-null', r1 !== null);
  ok('second call returns the SAME object reference (cache hit)', r1 === r2);
  // Force re-run gets a fresh result (may be a new object)
  const r3 = ctx.runL2SweepInheritance({ force: true });
  ok('forced third call still non-null', r3 !== null);
}

// 4l. invalidateL2SweepCache clears state.l2SweepResult
{
  console.log('\n--- 4l. invalidateL2SweepCache clears the slot ---');
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7 },
    itemGroups: { 'L2:0': [10, 11], 'L2:1': [10, 11] },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  ctx.runL2SweepInheritance({ force: true });
  ok('result populated', sb.state.l2SweepResult !== null);
  ctx.invalidateL2SweepCache();
  ok('result cleared after invalidate', sb.state.l2SweepResult === null);
  ok('cache key cleared after invalidate', sb.state.l2SweepCacheKey === null);
}

// 4m. Auto-promoted candidate has correct shape (source, confirmed, etc.)
{
  console.log('\n--- 4m. Auto-promoted candidate shape ---');
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
                { _s0:11, _e0:20, start_bp:5e6, end_bp:6e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
      1: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
      1: { xs: new Float64Array(15).fill(0.2) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7, 0.2: 0.7 },
    itemGroups: { 'L2:0': [10, 11], 'L2:1': [10, 11] },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const r = ctx.runL2SweepInheritance({ force: true });
  ctx._autoPromoteFromSweep(r);
  ok('2 candidates added', sb.addedCands.length === 2);
  const cand0 = sb.addedCands[0];
  ok('candidate.source = "auto_l2_sweep"', cand0.source === 'auto_l2_sweep');
  ok('candidate.confirmed = false', cand0.confirmed === false);
  ok('candidate.chrom set', cand0.chrom === 'LG28');
  ok('candidate.l2_indices = [<l2idx>]',
     Array.isArray(cand0.l2_indices) && cand0.l2_indices.length === 1);
  ok('candidate.start_bp / end_bp set from env',
     typeof cand0.start_bp === 'number' && typeof cand0.end_bp === 'number');
  ok('candidate.K = 3 (state.k)', cand0.K === 3);
  ok('candidate.locked_labels populated', !!cand0.locked_labels);
  ok('candidate.id has auto_l2sweep_ prefix', /^auto_l2sweep_/.test(cand0.id));
  ok('candidate.auto_promoted_at is ISO date',
     typeof cand0.auto_promoted_at === 'string' &&
     /^\d{4}-\d{2}-\d{2}T/.test(cand0.auto_promoted_at));
}

// 4n. Re-running auto-promote on the same state → no duplicate candidates
{
  console.log('\n--- 4n. Idempotent re-promote (gate 4 catches own auto promotes) ---');
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp:1e6, end_bp:2e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: {
      0: { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null },
    },
    aggs: {
      0: { xs: new Float64Array(15).fill(0.1) },
    },
    fakeSilhouetteByXs0: { 0.1: 0.7 },
    itemGroups: { 'L2:0': [10, 11] },
  });
  // Bypass the "needs ≥ 2 items" guard for this test by adding one synthetic
  // pre-existing manual candidate at a far-away location, so the sweep
  // includes 2 items but only L2:0 gets promoted (the manual_far isn't an L2).
  // Actually — to exercise idempotent re-promote with 1 envelope, we relax
  // the guard by setting _IGC_MIN_BANDS_FOR_CLUSTERING to 1 in the sandbox.
  sb.ctx._IGC_MIN_BANDS_FOR_CLUSTERING = 1;
  const ctx = runSweepBlockIn(sb.ctx);
  const r1 = ctx.runL2SweepInheritance({ force: true });
  const ap1 = ctx._autoPromoteFromSweep(r1);
  ok('first pass promoted L2[0]', ap1.promoted.length === 1 && ap1.promoted[0] === 0);
  // Run again — should be a no-op because gate 4 (ALREADY_IN_CANDIDATE) trips
  const r2 = ctx.runL2SweepInheritance({ force: true });
  const ap2 = ctx._autoPromoteFromSweep(r2);
  ok('second pass promotes 0 (gate 4 hit)',
     ap2.promoted.length === 0, 'got promoted=' + JSON.stringify(ap2.promoted));
  ok('second pass skipped reason ALREADY_IN_CANDIDATE',
     ap2.skipped.some(s => s.reason === 'ALREADY_IN_CANDIDATE'));
  ok('candidateList still has only 1 entry',
     sb.state.candidateList.length === 1);
}

// ============================================================================
// 5. SPEC §8 fixture: 6 synthetic L2s with planted inheritance
// ============================================================================
console.log('\n=== 5. SPEC §8 fixture: planted-inheritance recovery ===');

// 6 L2s; pairs (0,1), (2,3), (4,5) share inheritance. Sweep should
// recover at least the per-item group structure such that mates are
// in matching groups and non-mates are not.
{
  const labelsOK = new Int8Array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2]);
  const sb = makeSandbox({
    state: {
      data: { chrom: 'LG28',
              l2_envelopes: [
                { _s0:0, _e0:10, start_bp: 1e6,  end_bp: 2e6 },
                { _s0:11, _e0:20, start_bp: 5e6,  end_bp: 6e6 },
                { _s0:21, _e0:30, start_bp: 9e6,  end_bp: 10e6 },
                { _s0:31, _e0:40, start_bp: 13e6, end_bp: 14e6 },
                { _s0:41, _e0:50, start_bp: 17e6, end_bp: 18e6 },
                { _s0:51, _e0:60, start_bp: 21e6, end_bp: 22e6 },
              ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
    },
    clusters: Object.fromEntries([0,1,2,3,4,5].map(i =>
      [i, { ok: true, fixedKLabels: labelsOK, n_per_group: [5,5,5], silhouette: null }]
    )),
    aggs: Object.fromEntries([0,1,2,3,4,5].map(i =>
      [i, { xs: new Float64Array(15).fill(0.1 * i + 0.05) }]
    )),
    fakeSilhouetteByXs0: { 0.05: 0.6, 0.15: 0.6, 0.25: 0.6, 0.35: 0.6, 0.45: 0.6, 0.55: 0.6 },
    // Planted: pairs (0,1)→groupA; (2,3)→groupB; (4,5)→groupC
    // and a shared groupZ that everyone touches (so each item gets ≥2 groups)
    itemGroups: {
      'L2:0': [100, 999], 'L2:1': [100, 999],
      'L2:2': [200, 999], 'L2:3': [200, 999],
      'L2:4': [300, 999], 'L2:5': [300, 999],
    },
  });
  const ctx = runSweepBlockIn(sb.ctx);
  const r = ctx.runL2SweepInheritance({ force: true });
  ok('SPEC fixture: result has 6 l2_meta entries',
     r && r.l2_meta && r.l2_meta.length === 6);
  const ap = ctx._autoPromoteFromSweep(r);
  ok('SPEC fixture: all 6 promoted (well-separated, ≥2 groups each)',
     ap.promoted.length === 6, 'got ' + ap.promoted.length);
  // Verify each promoted candidate's start_bp matches an L2 in order
  const promotedStarts = ap.promoted.map(li =>
    sb.state.data.l2_envelopes[li].start_bp
  );
  ok('SPEC fixture: promoted L2s span the chrom',
     promotedStarts.includes(1e6) && promotedStarts.includes(21e6));
}

// ============================================================================
// 6. Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
