// =============================================================================
// turn 163 — Band-picker dropdown + runs TSV export + chain-break ticks
// =============================================================================
// Three small follow-ups from the turn 162 §6 NEXT list (items 4, 5, 6).
// All observation-only — no interpretation overlays, no ranked output.
//
// What this turn ships:
//   - linesBandTracePickSelect dropdown in the lines header (default
//     "largest" preserves turn 161 behaviour; per-band options expose
//     the focal candidate's bands by index + fish count)
//   - _updateBandTracePickOptions() — rebuild dropdown options when
//     the focal candidate changes; preserves selection when valid
//   - _readBandTracePickValue() — bridge that the 🔍 trace handler uses
//     to read the dropdown's current value; returns int|null
//   - 🔍 trace click handler now reads the dropdown and passes
//     {bandIdx: N} to _bandTraceFromFocalCandidate (or no opts for
//     "largest")
//   - linesBandTraceExportRunsBtn (📊 runs) in the lines header
//   - _bandTraceRunsToTSV(runs, opts) — pure serializer, one row per run
//   - _bandTraceDownloadRunsTSV() — Blob/anchor download path with
//     same headless-safe contract as the per-L2 TSV
//   - _drawBandTraceStrip: thin red vertical tick at every L2 where
//     chain_idx changes from the previous painted L2
//   - _BTRACE_CHAIN_BREAK_COLOR constant
//
// Tests:
//   1. Source-pattern checks (function defs, exports, DOM ids, constant)
//   2. _updateBandTracePickOptions — populates from candidate, preserves
//      selection, resets on candidate change, disables when no candidate
//   3. _readBandTracePickValue — null for "largest" / missing / invalid;
//      integer for valid band index
//   4. 🔍 trace handler reads dropdown via _readBandTracePickValue
//   5. _bandTraceRunsToTSV — header columns, one row per run, NaN-safe
//   6. _bandTraceRunsToTSV — empty runs → header only; null runs → null
//   7. _bandTraceRunsToTSV — start_bp / end_bp from envelopes
//   8. _bandTraceDownloadRunsTSV — filename format, returns null with no fish-set
//   9. 📊 runs button DOM + click handler
//  10. Chain-break tick — drawn when chain_idx changes
//  11. Chain-break tick — NOT drawn when chain_idx stable
//  12. Chain-break tick — NOT drawn at first painted L2
//  13. _BTRACE_CHAIN_BREAK_COLOR constant + window export
//  14. Regression — all turn 160/161/162 surfaces still wired
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

ok('_updateBandTracePickOptions defined',
   /function\s+_updateBandTracePickOptions\s*\(/.test(html));
ok('_readBandTracePickValue defined',
   /function\s+_readBandTracePickValue\s*\(/.test(html));
ok('_bandTraceRunsToTSV defined',
   /function\s+_bandTraceRunsToTSV\s*\(/.test(html));
ok('_bandTraceDownloadRunsTSV defined',
   /function\s+_bandTraceDownloadRunsTSV\s*\(/.test(html));

ok('window._updateBandTracePickOptions exported',
   /window\._updateBandTracePickOptions\s*=\s*_updateBandTracePickOptions/.test(html));
ok('window._readBandTracePickValue exported',
   /window\._readBandTracePickValue\s*=\s*_readBandTracePickValue/.test(html));
ok('window._bandTraceRunsToTSV exported',
   /window\._bandTraceRunsToTSV\s*=\s*_bandTraceRunsToTSV/.test(html));
ok('window._bandTraceDownloadRunsTSV exported',
   /window\._bandTraceDownloadRunsTSV\s*=\s*_bandTraceDownloadRunsTSV/.test(html));

ok('linesBandTracePickSelect dropdown in lines header',
   html.indexOf('id="linesBandTracePickSelect"') > -1);
ok('linesBandTraceExportRunsBtn button in lines header',
   html.indexOf('id="linesBandTraceExportRunsBtn"') > -1);
ok('dropdown has default "largest" option',
   /id="linesBandTracePickSelect"[\s\S]{0,1500}<option\s+value="largest">largest<\/option>/.test(html));
ok('runs export button has 📊 runs label',
   /id="linesBandTraceExportRunsBtn"[\s\S]{0,1500}📊 runs/.test(html));

ok('_BTRACE_CHAIN_BREAK_COLOR constant declared',
   /_BTRACE_CHAIN_BREAK_COLOR\s*=\s*['"]#[0-9a-fA-F]{6}['"]/.test(html));
ok('_BTRACE_CHAIN_BREAK_COLOR exported on window',
   /window\._BTRACE_CHAIN_BREAK_COLOR\s*=\s*_BTRACE_CHAIN_BREAK_COLOR/.test(html));

const fnDraw = pullFunction(html, '_drawBandTraceStrip');
ok('_drawBandTraceStrip references chain-break color',
   fnDraw && /_BTRACE_CHAIN_BREAK_COLOR/.test(fnDraw));
ok('_drawBandTraceStrip tracks prevChainIdx',
   fnDraw && /prevChainIdx/.test(fnDraw));
ok('_drawBandTraceStrip prevChainIdx initialised to -1',
   fnDraw && /let\s+prevChainIdx\s*=\s*-1/.test(fnDraw));

// Trace button reads dropdown
ok('🔍 trace click handler calls _readBandTracePickValue',
   /_readBandTracePickValue\(\)/.test(html));
// Runs export button click handler
ok('runs export button click handler calls _bandTraceDownloadRunsTSV',
   /linesBandTraceExportRunsBtn[\s\S]{0,800}_bandTraceDownloadRunsTSV\(\)/.test(html));

// ============================================================================
// 2. _updateBandTracePickOptions
// ============================================================================
console.log('\n=== 2. _updateBandTracePickOptions ===');

// Build a sandbox with a fake document that mirrors the bits the
// helper uses (getElementById, createElement, sel.options, sel.remove).
const sbx = {
  console,
  state: { candidate: null, k: 3 },
  window: {},
};
sbx.window.state = sbx.state;

vm.createContext(sbx);
vm.runInContext(`
  function _makeFakeSelect() {
    const sel = {
      _opts: [],
      value: '',
      disabled: false,
      get options() { return this._opts; },
      remove(idx) { this._opts.splice(idx, 1); },
      appendChild(opt) { this._opts.push(opt); },
    };
    return sel;
  }
  function _makeFakeOption() {
    return { value: '', textContent: '' };
  }
  var _selStore = null;
  var document = {
    getElementById(id) {
      if (id === 'linesBandTracePickSelect') return _selStore;
      return null;
    },
    createElement(tag) {
      if (tag === 'option') return _makeFakeOption();
      return null;
    },
  };
  // Seed the select with the default "largest" option as the source HTML does.
  _selStore = _makeFakeSelect();
  _selStore._opts.push({ value: 'largest', textContent: 'largest' });
  _selStore.value = 'largest';
`, sbx);

const fnUpdate = pullFunction(html, '_updateBandTracePickOptions');
ok('pulled _updateBandTracePickOptions', !!fnUpdate);
vm.runInContext(fnUpdate, sbx);

// No candidate → only "largest" option, disabled
sbx.state.candidate = null;
vm.runInContext('_updateBandTracePickOptions()', sbx);
ok('no candidate: only "largest" option remains',
   sbx.window && (sbx.eval || true) && vm.runInContext('_selStore._opts.length', sbx) === 1);
ok('no candidate: dropdown disabled',
   vm.runInContext('_selStore.disabled', sbx) === true);
ok('no candidate: value reset to "largest"',
   vm.runInContext('_selStore.value', sbx) === 'largest');

// With a candidate K=3 (labels = [0,0,0,1,1,2]) → 4 options total
sbx.state.candidate = {
  id: 'I1', K: 3,
  locked_labels: [0, 0, 0, 1, 1, 2],
};
vm.runInContext('_updateBandTracePickOptions()', sbx);
ok('candidate K=3: 4 options total (largest + b0/b1/b2)',
   vm.runInContext('_selStore._opts.length', sbx) === 4);
ok('candidate K=3: dropdown enabled',
   vm.runInContext('_selStore.disabled', sbx) === false);
ok('candidate K=3: option b0 (n=3)',
   vm.runInContext('_selStore._opts[1].textContent', sbx) === 'b0 (n=3)');
ok('candidate K=3: option b1 (n=2)',
   vm.runInContext('_selStore._opts[2].textContent', sbx) === 'b1 (n=2)');
ok('candidate K=3: option b2 (n=1)',
   vm.runInContext('_selStore._opts[3].textContent', sbx) === 'b2 (n=1)');
ok('candidate K=3: option values are stringified band indices',
   vm.runInContext('_selStore._opts[1].value', sbx) === '0' &&
   vm.runInContext('_selStore._opts[2].value', sbx) === '1' &&
   vm.runInContext('_selStore._opts[3].value', sbx) === '2');

// User picks b1, then candidate changes to a smaller K=3 with b1 still valid:
// selection should be preserved.
vm.runInContext('_selStore.value = "1"', sbx);
sbx.state.candidate = {
  id: 'I2', K: 3,
  locked_labels: [0, 1, 1, 1, 2, 2],
};
vm.runInContext('_updateBandTracePickOptions()', sbx);
ok('preserves selection when band index still valid (b1 → b1)',
   vm.runInContext('_selStore.value', sbx) === '1');

// User picks b2, then candidate changes to K=2 (labels in {0,1}):
// b2 is no longer valid, selection should reset to "largest".
vm.runInContext('_selStore.value = "2"', sbx);
sbx.state.candidate = {
  id: 'I3', K: 2,
  locked_labels: [0, 0, 1, 1, 1],
};
vm.runInContext('_updateBandTracePickOptions()', sbx);
ok('resets to "largest" when previous band index invalid',
   vm.runInContext('_selStore.value', sbx) === 'largest');
ok('K=2 candidate: 3 options total (largest + b0/b1)',
   vm.runInContext('_selStore._opts.length', sbx) === 3);

// Repopulating doesn't keep stale options around.
sbx.state.candidate = {
  id: 'I4', K: 3,
  locked_labels: [0, 0, 1, 1, 1, 2],
};
vm.runInContext('_updateBandTracePickOptions()', sbx);
vm.runInContext('_updateBandTracePickOptions()', sbx);  // idempotent
ok('idempotent: still 4 options after second call',
   vm.runInContext('_selStore._opts.length', sbx) === 4);

// Candidate with empty locked_labels → treated as no candidate
sbx.state.candidate = { id: 'I5', K: 3, locked_labels: [] };
vm.runInContext('_updateBandTracePickOptions()', sbx);
ok('empty locked_labels: only "largest" option', vm.runInContext('_selStore._opts.length', sbx) === 1);
ok('empty locked_labels: dropdown disabled', vm.runInContext('_selStore.disabled', sbx) === true);

// ============================================================================
// 3. _readBandTracePickValue
// ============================================================================
console.log('\n=== 3. _readBandTracePickValue ===');

const fnRead = pullFunction(html, '_readBandTracePickValue');
ok('pulled _readBandTracePickValue', !!fnRead);
vm.runInContext(fnRead, sbx);

vm.runInContext('_selStore.value = "largest"', sbx);
ok('"largest" → null', vm.runInContext('_readBandTracePickValue()', sbx) === null);

vm.runInContext('_selStore.value = "0"', sbx);
ok('"0" → 0', vm.runInContext('_readBandTracePickValue()', sbx) === 0);

vm.runInContext('_selStore.value = "2"', sbx);
ok('"2" → 2', vm.runInContext('_readBandTracePickValue()', sbx) === 2);

vm.runInContext('_selStore.value = ""', sbx);
ok('"" → null', vm.runInContext('_readBandTracePickValue()', sbx) === null);

vm.runInContext('_selStore.value = "garbage"', sbx);
ok('non-numeric value → null', vm.runInContext('_readBandTracePickValue()', sbx) === null);

vm.runInContext('_selStore.value = "-1"', sbx);
ok('negative value → null', vm.runInContext('_readBandTracePickValue()', sbx) === null);

// Missing dropdown → null
vm.runInContext('document.getElementById = function(){ return null; }', sbx);
ok('missing dropdown → null', vm.runInContext('_readBandTracePickValue()', sbx) === null);

// ============================================================================
// 5, 6, 7. _bandTraceRunsToTSV
// ============================================================================
console.log('\n=== 5-7. _bandTraceRunsToTSV ===');

const sbx2 = { console };
vm.createContext(sbx2);
const fnRunsTSV = pullFunction(html, '_bandTraceRunsToTSV');
ok('pulled _bandTraceRunsToTSV', !!fnRunsTSV);
vm.runInContext(fnRunsTSV, sbx2);

// Build synthetic runs (output shape of _bandTraceRegimeRuns)
vm.runInContext(`
  var _runs = [
    { start_l2_idx: 5, end_l2_idx: 8,
      start_chain_position: 5, end_chain_position: 8,
      chain_idx: 0, n_L2: 4, n_co_seg: 3, n_partial: 1,
      dominant_band: 0, mean_dominant_fraction: 0.825, mean_entropy: 0.235 },
    { start_l2_idx: 14, end_l2_idx: 17,
      start_chain_position: 14, end_chain_position: 17,
      chain_idx: 0, n_L2: 4, n_co_seg: 4, n_partial: 0,
      dominant_band: 2, mean_dominant_fraction: 0.91, mean_entropy: 0.13 },
    { start_l2_idx: 22, end_l2_idx: 23,
      start_chain_position: 0, end_chain_position: 1,
      chain_idx: 1, n_L2: 2, n_co_seg: 2, n_partial: 0,
      dominant_band: 1, mean_dominant_fraction: 0.95, mean_entropy: 0.075 },
  ];
  var _envs = [];
  for (var i = 0; i < 30; i++) {
    _envs.push({ start_bp: i * 1_000_000, end_bp: (i + 1) * 1_000_000 });
  }
`, sbx2);

const tsv = vm.runInContext('_bandTraceRunsToTSV(_runs, {chrom: "LG28", envelopes: _envs})', sbx2);
const tsvLines = tsv.split('\n').filter(l => l.length > 0);

ok('runs TSV non-empty', tsvLines.length > 0);
ok('runs TSV header begins with chrom', tsvLines[0].split('\t')[0] === 'chrom');
ok('runs TSV header has core columns',
   /chrom\trun_idx\tchain_idx\tstart_l2_idx\tend_l2_idx\tstart_chain_position\tend_chain_position\tstart_bp\tend_bp\tn_L2\tn_co_seg\tn_partial\tdominant_band\tmean_dominant_fraction\tmean_entropy/
     .test(tsvLines[0]));
ok('runs TSV has 3 data rows', tsvLines.length === 4);

const r0 = tsvLines[1].split('\t');
ok('row 0 chrom = LG28', r0[0] === 'LG28');
ok('row 0 run_idx = 0', r0[1] === '0');
ok('row 0 chain_idx = 0', r0[2] === '0');
ok('row 0 start_l2_idx = 5', r0[3] === '5');
ok('row 0 end_l2_idx = 8', r0[4] === '8');
ok('row 0 start_bp = 5000000 (envelope 5)', r0[7] === '5000000');
ok('row 0 end_bp = 9000000 (envelope 8)', r0[8] === '9000000');
ok('row 0 n_L2 = 4', r0[9] === '4');
ok('row 0 n_co_seg = 3', r0[10] === '3');
ok('row 0 n_partial = 1', r0[11] === '1');
ok('row 0 dominant_band = 0', r0[12] === '0');
ok('row 0 mean_dominant_fraction formatted to 6 dp',
   /^0\.825\d{3}$/.test(r0[13]));
ok('row 0 mean_entropy formatted to 6 dp',
   /^0\.235\d{3}$/.test(r0[14]));

const r2 = tsvLines[3].split('\t');
ok('row 2 chain_idx = 1 (chain break)', r2[2] === '1');
ok('row 2 dominant_band = 1', r2[12] === '1');
ok('row 2 run_idx = 2', r2[1] === '2');

// Empty runs → header only
const tsv_empty = vm.runInContext('_bandTraceRunsToTSV([], {chrom: "LG14"})', sbx2);
const lines_empty = tsv_empty.split('\n').filter(l => l.length > 0);
ok('empty runs → header only', lines_empty.length === 1);

// Null runs → null
ok('null runs → null', vm.runInContext('_bandTraceRunsToTSV(null, {})', sbx2) === null);

// Missing envelopes → blank bp columns
const tsv_no_env = vm.runInContext('_bandTraceRunsToTSV(_runs, {chrom: "LG14"})', sbx2);
const lines_no_env = tsv_no_env.split('\n').filter(l => l.length > 0);
const r0_no_env = lines_no_env[1].split('\t');
ok('no envelopes: start_bp blank', r0_no_env[7] === '');
ok('no envelopes: end_bp blank', r0_no_env[8] === '');
ok('no envelopes: rest of row still populated', r0_no_env[2] === '0' && r0_no_env[9] === '4');

// Run with NaN means → empty in TSV
vm.runInContext(`
  var _badRuns = [
    { start_l2_idx: 0, end_l2_idx: 1, chain_idx: 0,
      start_chain_position: 0, end_chain_position: 1,
      n_L2: 2, n_co_seg: 2, n_partial: 0,
      dominant_band: 0, mean_dominant_fraction: NaN, mean_entropy: undefined },
  ];
`, sbx2);
const tsv_bad = vm.runInContext('_bandTraceRunsToTSV(_badRuns, {chrom: "LG28"})', sbx2);
const lines_bad = tsv_bad.split('\n').filter(l => l.length > 0);
const r_bad = lines_bad[1].split('\t');
ok('NaN mean_dominant_fraction → empty', r_bad[13] === '');
ok('undefined mean_entropy → empty', r_bad[14] === '');

// ============================================================================
// 8. _bandTraceDownloadRunsTSV — filename + null path
// ============================================================================
console.log('\n=== 8. _bandTraceDownloadRunsTSV ===');

// Build a sandbox that includes the compute layer + cache + clearCache.
const sbx3 = {
  console,
  state: {
    bandTraceFishSet: null,
    bandTraceOn: false,
    bandTraceCache: null,
    bandTraceCacheKey: null,
    _btraceHits: null,
    k: 3,
    candidate: null,
    data: null,
  },
  window: {},
};
sbx3.window.state = sbx3.state;
sbx3.localStorage = (function() {
  const _s = {};
  return { setItem: (k,v)=>{_s[k]=String(v);}, getItem: k=>(k in _s?_s[k]:null), removeItem: k=>{delete _s[k];}, _store:_s };
})();
vm.createContext(sbx3);
vm.runInContext('var localStorage = this.localStorage;', sbx3);

vm.runInContext(`
  function _bandTraceForFishSet(fishSet, opts) {
    if (!fishSet || !fishSet.length) return null;
    const K = (opts && opts.K) || 3;
    const l2 = (opts && opts.l2_indices) || [];
    return {
      n_fish_selected: fishSet.length,
      n_chains: 1, n_total_L2: l2.length, K,
      per_l2: l2.map((idx, i) => ({
        l2_idx: idx, chain_idx: 0, chain_position: i,
        n_valid: fishSet.length,
        band_counts: new Int32Array(K),
        band_fractions: (function(){ const f = new Float32Array(K); f[idx % K] = 1; return f; })(),
        entropy: 0.05, dominant_band: idx % K, dominant_fraction: 1, regime: 'co_seg',
      })),
    };
  }
  function _bandTraceRegimeRuns(trace, opts) {
    // Stub: return a single run covering the whole trace.
    if (!trace || !trace.per_l2 || trace.per_l2.length === 0) return [];
    const first = trace.per_l2[0];
    const last  = trace.per_l2[trace.per_l2.length - 1];
    return [{
      start_l2_idx: first.l2_idx, end_l2_idx: last.l2_idx,
      start_chain_position: 0, end_chain_position: trace.per_l2.length - 1,
      chain_idx: 0, n_L2: trace.per_l2.length,
      n_co_seg: trace.per_l2.length, n_partial: 0,
      dominant_band: 0, mean_dominant_fraction: 1, mean_entropy: 0,
    }];
  }
  const _BTRACE_ON_LS_KEY        = 'inversion_atlas.bandTraceOn';
  const _BTRACE_FISH_SET_LS_KEY  = 'inversion_atlas.bandTraceFishSet';
`, sbx3);

const fnCacheKey = pullFunction(html, '_bandTraceCacheKey');
const fnGoC      = pullFunction(html, '_bandTraceGetOrCompute');
const fnRunsTSV2 = pullFunction(html, '_bandTraceRunsToTSV');
const fnDl       = pullFunction(html, '_bandTraceDownloadRunsTSV');
vm.runInContext(fnCacheKey, sbx3);
vm.runInContext(fnGoC, sbx3);
vm.runInContext(fnRunsTSV2, sbx3);
vm.runInContext(fnDl, sbx3);

// No fish-set → null
ok('no fish-set → returns null',
   vm.runInContext('_bandTraceDownloadRunsTSV()', sbx3) === null);

// With cached trace + fish-set + envelopes — headless: returns filename
sbx3.state.data = {
  chrom: 'LG28',
  l2_envelopes: [
    { start_bp: 1_000_000, end_bp: 2_000_000 },
    { start_bp: 2_500_000, end_bp: 3_000_000 },
    { start_bp: 3_500_000, end_bp: 4_000_000 },
  ],
};
sbx3.state.bandTraceFishSet = [1, 2, 3, 4];
sbx3.state.bandTraceCache = null;
sbx3.state.bandTraceCacheKey = null;
const fname = vm.runInContext('_bandTraceDownloadRunsTSV()', sbx3);
ok('returns a filename string', typeof fname === 'string' && fname.length > 0);
ok('filename starts with band_trace_runs_LG28', fname.indexOf('band_trace_runs_LG28') === 0);
ok('filename embeds n_fish (n4)', /_n4_/.test(fname));
ok('filename embeds K (K3)', /_K3\.tsv$/.test(fname));

// Empty l2_envelopes → null
sbx3.state.data = { chrom: 'LG28', l2_envelopes: [] };
sbx3.state.bandTraceCache = null;
sbx3.state.bandTraceCacheKey = null;
ok('empty l2_envelopes → null',
   vm.runInContext('_bandTraceDownloadRunsTSV()', sbx3) === null);

// ============================================================================
// 10, 11, 12. Chain-break tick rendering
// ============================================================================
console.log('\n=== 10-12. Chain-break tick ===');

// Build a fresh sandbox to exercise the strip drawer with a multi-chain trace.
const sbx4 = {
  console,
  state: {
    bandTraceFishSet: [1, 2, 3],
    bandTraceOn: true,
    bandTraceCache: null,
    bandTraceCacheKey: null,
    _btraceHits: null,
    k: 3,
    data: {
      chrom: 'LG28',
      l2_envelopes: [
        { start_bp: 1_000_000, end_bp: 2_000_000 },
        { start_bp: 2_500_000, end_bp: 3_000_000 },
        { start_bp: 3_500_000, end_bp: 4_000_000 },
        { start_bp: 4_500_000, end_bp: 5_000_000 },
      ],
    },
  },
  window: {},
};
sbx4.window.state = sbx4.state;
sbx4.localStorage = (function() {
  const _s = {};
  return { setItem: (k,v)=>{_s[k]=String(v);}, getItem: k=>(k in _s?_s[k]:null), removeItem: k=>{delete _s[k];}, _store:_s };
})();
vm.createContext(sbx4);
vm.runInContext('var localStorage = this.localStorage;', sbx4);

// Inject a stub _bandTraceForFishSet that produces a chain break between
// L2 1 and L2 2 (L2s 0,1 in chain 0; L2s 2,3 in chain 1).
vm.runInContext(`
  function _bandTraceForFishSet(fishSet, opts) {
    if (!fishSet || !fishSet.length) return null;
    const K = (opts && opts.K) || 3;
    const l2 = (opts && opts.l2_indices) || [];
    return {
      n_fish_selected: fishSet.length,
      n_chains: 2, n_total_L2: l2.length, K,
      per_l2: l2.map((idx, i) => ({
        l2_idx: idx, chain_idx: i < 2 ? 0 : 1, chain_position: i < 2 ? i : i - 2,
        n_valid: fishSet.length,
        band_counts: new Int32Array(K),
        band_fractions: (function(){ const f = new Float32Array(K); f[idx % K] = 1; return f; })(),
        entropy: 0.05, dominant_band: idx % K, dominant_fraction: 1, regime: 'co_seg',
      })),
    };
  }
  const _BTRACE_REGIME_COLOR = {co_seg:'#7ad394',partial:'#f5a524',fanned:'#9aa3ad',sparse:'#5a6068',no_valid:'#3a3e44'};
  const _BTRACE_CHAIN_BREAK_COLOR = '#e85a5a';
  const _BTRACE_ON_LS_KEY        = 'inversion_atlas.bandTraceOn';
  const _BTRACE_FISH_SET_LS_KEY  = 'inversion_atlas.bandTraceFishSet';
  function drawLinesPanel() {}
`, sbx4);

vm.runInContext(pullFunction(html, '_bandTraceCacheKey'), sbx4);
vm.runInContext(pullFunction(html, '_bandTraceGetOrCompute'), sbx4);
vm.runInContext(pullFunction(html, '_drawBandTraceStrip'), sbx4);

// Fake ctx that records every fillRect with its fillStyle at call time.
vm.runInContext(`
  function makeCtx() {
    const calls = [];
    let _curFill = null;
    return {
      save: () => calls.push(['save']),
      restore: () => calls.push(['restore']),
      fillRect: (...a) => calls.push(['fillRect', _curFill, ...a]),
      strokeRect: (...a) => calls.push(['strokeRect', ...a]),
      set fillStyle(v) { _curFill = v; },
      get fillStyle() { return _curFill; },
      set strokeStyle(v) {},
      set lineWidth(v) {},
      set globalAlpha(v) {},
      _calls: calls,
    };
  }
  var _ctx = makeCtx();
  _drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 6);
`, sbx4);

// Count fillRect calls with the chain-break color
const tickCalls = vm.runInContext(`
  _ctx._calls.filter(c => c[0] === 'fillRect' && c[1] === '#e85a5a')
`, sbx4);
ok('chain-break tick: exactly one tick painted (between L2 1 and L2 2)',
   Array.isArray(tickCalls) && tickCalls.length === 1);
ok('chain-break tick: width = 1 (single pixel)',
   tickCalls.length > 0 && tickCalls[0][4] === 1);
ok('chain-break tick: spans full strip height (h=7)',
   tickCalls.length > 0 && tickCalls[0][5] === 7);

// Now exercise a SINGLE-chain trace — no ticks should be drawn.
vm.runInContext(`
  function _bandTraceForFishSet(fishSet, opts) {
    if (!fishSet || !fishSet.length) return null;
    const K = (opts && opts.K) || 3;
    const l2 = (opts && opts.l2_indices) || [];
    return {
      n_fish_selected: fishSet.length,
      n_chains: 1, n_total_L2: l2.length, K,
      per_l2: l2.map((idx, i) => ({
        l2_idx: idx, chain_idx: 0, chain_position: i,
        n_valid: fishSet.length,
        band_counts: new Int32Array(K),
        band_fractions: (function(){ const f = new Float32Array(K); f[idx % K] = 1; return f; })(),
        entropy: 0.05, dominant_band: idx % K, dominant_fraction: 1, regime: 'co_seg',
      })),
    };
  }
`, sbx4);

// Reset cache + replay
sbx4.state.bandTraceCache = null;
sbx4.state.bandTraceCacheKey = null;
sbx4.state._btraceHits = null;
vm.runInContext(`
  _ctx = makeCtx();
  _drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 6);
`, sbx4);
const tickCalls2 = vm.runInContext(`
  _ctx._calls.filter(c => c[0] === 'fillRect' && c[1] === '#e85a5a')
`, sbx4);
ok('single chain: no chain-break ticks',
   Array.isArray(tickCalls2) && tickCalls2.length === 0);

// First-painted-L2 never gets a tick — even if chain_idx !== 0
vm.runInContext(`
  function _bandTraceForFishSet(fishSet, opts) {
    if (!fishSet || !fishSet.length) return null;
    const K = (opts && opts.K) || 3;
    const l2 = (opts && opts.l2_indices) || [];
    return {
      n_fish_selected: fishSet.length,
      n_chains: 1, n_total_L2: l2.length, K,
      per_l2: l2.map((idx, i) => ({
        l2_idx: idx, chain_idx: 7,   // non-zero from the start
        chain_position: i,
        n_valid: fishSet.length,
        band_counts: new Int32Array(K),
        band_fractions: (function(){ const f = new Float32Array(K); f[idx % K] = 1; return f; })(),
        entropy: 0.05, dominant_band: idx % K, dominant_fraction: 1, regime: 'co_seg',
      })),
    };
  }
`, sbx4);
sbx4.state.bandTraceCache = null;
sbx4.state.bandTraceCacheKey = null;
sbx4.state._btraceHits = null;
vm.runInContext(`
  _ctx = makeCtx();
  _drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 6);
`, sbx4);
const tickCalls3 = vm.runInContext(`
  _ctx._calls.filter(c => c[0] === 'fillRect' && c[1] === '#e85a5a')
`, sbx4);
ok('first painted L2 with chain_idx=7 still no tick (no previous to compare)',
   Array.isArray(tickCalls3) && tickCalls3.length === 0);

// ============================================================================
// 13. Constant + window export
// ============================================================================
console.log('\n=== 13. Constants ===');

ok('_BTRACE_CHAIN_BREAK_COLOR is a hex color',
   /_BTRACE_CHAIN_BREAK_COLOR\s*=\s*['"]#[0-9a-fA-F]{6}['"]/.test(html));

// ============================================================================
// 14. Regression
// ============================================================================
console.log('\n=== 14. Regression ===');

// Turn 162
ok('turn 162 _bandTraceToTSV still defined',
   /function\s+_bandTraceToTSV\s*\(/.test(html));
ok('turn 162 _bandTraceDownloadTSV still defined',
   /function\s+_bandTraceDownloadTSV\s*\(/.test(html));
ok('turn 162 _wireBandTraceTooltip still defined',
   /function\s+_wireBandTraceTooltip\s*\(/.test(html));
ok('turn 162 linesBandTraceExportBtn still in lines header',
   html.indexOf('id="linesBandTraceExportBtn"') > -1);

// Turn 161
ok('turn 161 _drawBandTraceStrip still defined',
   /function\s+_drawBandTraceStrip\s*\(/.test(html));
ok('turn 161 setBandTraceFishSet still defined',
   /function\s+setBandTraceFishSet\s*\(/.test(html));
ok('turn 161 setBandTraceOn still defined',
   /function\s+setBandTraceOn\s*\(/.test(html));
ok('turn 161 _bandTraceFromFocalCandidate still defined',
   /function\s+_bandTraceFromFocalCandidate\s*\(/.test(html));
ok('turn 161 linesBandTraceToggle still in lines header',
   html.indexOf('id="linesBandTraceToggle"') > -1);
ok('turn 161 linesBandTraceFromCandBtn still in lines header',
   html.indexOf('id="linesBandTraceFromCandBtn"') > -1);

// Turn 160
ok('turn 160 _bandTraceForFishSet still defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('turn 160 _bandTraceRegimeRuns still defined',
   /function\s+_bandTraceRegimeRuns\s*\(/.test(html));

// Turn 130
ok('turn 130 _drawLineageStrip still defined',
   /function\s+_drawLineageStrip\s*\(/.test(html));
ok('turn 130 lineage toggle still in lines header',
   html.indexOf('id="linesLineageStripToggle"') > -1);

// Turn 157A / 122
ok('turn 157A V-shape tooltip still defined',
   /function\s+_wireVShapeTooltip\s*\(/.test(html));
ok('turn 122 inheritance pill tooltip still defined',
   /function\s+_wireInheritancePillTooltip\s*\(/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
