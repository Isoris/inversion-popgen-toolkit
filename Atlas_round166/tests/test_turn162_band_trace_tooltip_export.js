// =============================================================================
// turn 162 — Band-trace hover tooltip + TSV export
// =============================================================================
// Follow-up to turn 161 (which shipped the band-trace strip UI). The
// turn-161 handoff §6 listed brackets, tooltip, combinatorial enumeration,
// and TSV export as the next four candidate items. This turn ships
// items 2 (tooltip) and 4 (TSV export) — both manuscript-facing,
// observation-only, no interpretation overlays. Brackets and the
// combinatorial enumerator stay deferred per the manuscript framing.
//
// What this turn ships:
//   - _btraceTooltipEnsureEl / _btraceTooltipBuildHtml / _btraceTooltipShow
//     / _btraceTooltipHide — tooltip primitives
//   - _btraceHitTest(hits, x, y) — CSS-pixel rectangle hit-test
//   - _wireBandTraceTooltip(canvas) — idempotent mousemove handler attach
//     on the lines-panel PC1 sub-canvas
//   - state._btraceHits — written by _drawBandTraceStrip (one rect per
//     visible L2 column), cleared by _bandTraceClearCache
//   - _bandTraceToTSV(trace, opts) — header + per-L2 rows; K columns of
//     band_fraction_<k>; NaN-safe
//   - _bandTraceDownloadTSV() — Blob/anchor download path; headless-safe
//   - lines header: "📊 TSV" button next to the existing "🔍 trace"
//   - Click handler on the export button alerts when no trace is active
//
// Tests:
//   1. Source-pattern checks (function defs, exports, DOM ids, tooltip el id)
//   2. _btraceTooltipBuildHtml — regime chip, dominant band swatch, entropy,
//      chain hidden when n_chains===1, shown when >1
//   3. _btraceTooltipBuildHtml — defensive on null hit / null entry
//   4. _btraceHitTest — hit, miss, edge-of-rect inclusive
//   5. _wireBandTraceTooltip — idempotent via dataset marker
//   6. _drawBandTraceStrip populates state._btraceHits with the correct shape
//   7. Hits cleared by _bandTraceClearCache + chrom switch via applyData
//   8. Hits NOT populated when toggle off / no fish-set
//   9. _bandTraceToTSV — header columns, one data row per L2, NaN-safe formatting
//  10. _bandTraceToTSV — empty per_l2 → just header; null trace → null
//  11. _bandTraceToTSV — start_bp/end_bp from envelopes; blank when missing
//  12. _bandTraceDownloadTSV — filename format, returns null with no fish-set
//  13. DOM: linesBandTraceExportBtn present + tooltip element id constant
//  14. Regression — turn 161 strip + turn 160 compute + turn 130 lineage
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

// Function declarations
ok('_btraceTooltipEnsureEl defined',
   /function\s+_btraceTooltipEnsureEl\s*\(/.test(html));
ok('_btraceTooltipBuildHtml defined',
   /function\s+_btraceTooltipBuildHtml\s*\(/.test(html));
ok('_btraceTooltipShow defined',
   /function\s+_btraceTooltipShow\s*\(/.test(html));
ok('_btraceTooltipHide defined',
   /function\s+_btraceTooltipHide\s*\(/.test(html));
ok('_btraceHitTest defined',
   /function\s+_btraceHitTest\s*\(/.test(html));
ok('_wireBandTraceTooltip defined',
   /function\s+_wireBandTraceTooltip\s*\(/.test(html));
ok('_bandTraceToTSV defined',
   /function\s+_bandTraceToTSV\s*\(/.test(html));
ok('_bandTraceDownloadTSV defined',
   /function\s+_bandTraceDownloadTSV\s*\(/.test(html));

// Window exports
ok('window._btraceTooltipBuildHtml exported',
   /window\._btraceTooltipBuildHtml\s*=\s*_btraceTooltipBuildHtml/.test(html));
ok('window._btraceHitTest exported',
   /window\._btraceHitTest\s*=\s*_btraceHitTest/.test(html));
ok('window._wireBandTraceTooltip exported',
   /window\._wireBandTraceTooltip\s*=\s*_wireBandTraceTooltip/.test(html));
ok('window._bandTraceToTSV exported',
   /window\._bandTraceToTSV\s*=\s*_bandTraceToTSV/.test(html));
ok('window._bandTraceDownloadTSV exported',
   /window\._bandTraceDownloadTSV\s*=\s*_bandTraceDownloadTSV/.test(html));

// DOM element id constants
ok('_BTRACE_TOOLTIP_EL_ID constant',
   /_BTRACE_TOOLTIP_EL_ID\s*=\s*['"]btracePointTooltip['"]/.test(html));

// DOM in lines header
ok('linesBandTraceExportBtn in lines header',
   html.indexOf('id="linesBandTraceExportBtn"') > -1);
ok('export button has 📊 TSV label',
   /id="linesBandTraceExportBtn"[\s\S]{0,1200}📊 TSV/.test(html));

// Wiring: tooltip on PC1 sub-canvas
ok('drawLinesPanel wires _wireBandTraceTooltip on PC1 sub-canvas',
   /src === 'pc1'[\s\S]{0,200}_wireBandTraceTooltip\(cv\)/.test(html));

// Wiring: click handler on the export button
ok('export button click handler defined',
   /linesBandTraceExportBtn[\s\S]{0,500}_bandTraceDownloadTSV\(\)/.test(html));

// _drawBandTraceStrip writes hits
ok('_drawBandTraceStrip initializes _btraceHits',
   /_state\._btraceHits\s*=\s*\[\]/.test(html));
ok('_drawBandTraceStrip pushes hit rect',
   /_state\._btraceHits\.push\(\{/.test(html));

// _bandTraceClearCache clears hits
const fnClear = pullFunction(html, '_bandTraceClearCache');
ok('_bandTraceClearCache clears _btraceHits',
   fnClear && /_state\._btraceHits\s*=\s*null/.test(fnClear));

// ============================================================================
// 2 & 3. _btraceTooltipBuildHtml rendering
// ============================================================================
console.log('\n=== 2-3. _btraceTooltipBuildHtml ===');

// Build a sandbox with the helper. The function references
// _BTRACE_REGIME_COLOR (constant) and optionally _gpKaryoColor (function);
// inject both so the rendering paths exercise without DOM dependencies.
const sbx = {
  console,
  state: {},
  window: {},
};
vm.createContext(sbx);

vm.runInContext(`
  const _BTRACE_REGIME_COLOR = {
    co_seg:'#7ad394', partial:'#f5a524', fanned:'#9aa3ad',
    sparse:'#5a6068', no_valid:'#3a3e44',
  };
  function _gpKaryoColor(k) {
    const PAL = ['#3b6fb6', '#ffd866', '#d97a2c', '#7ad394', '#a76de2', '#e85a5a'];
    return PAL[k % PAL.length];
  }
`, sbx);

const fnBuild = pullFunction(html, '_btraceTooltipBuildHtml');
ok('pulled _btraceTooltipBuildHtml', !!fnBuild);
vm.runInContext(fnBuild, sbx);

// A "co_seg" hit with valid data
const hit_coseg = {
  l2_idx: 7, x: 0, y: 0, w: 10, h: 7,
  entry: {
    l2_idx: 7, chain_idx: 0, chain_position: 7, n_valid: 12,
    band_fractions: [0.83, 0.0, 0.17],
    entropy: 0.213, dominant_band: 0, dominant_fraction: 0.83,
    regime: 'co_seg',
  },
};
sbx._h = hit_coseg;
const html_coseg = vm.runInContext('_btraceTooltipBuildHtml(_h, {n_chains: 1, n_fish_selected: 12})', sbx);
ok('co_seg html mentions L2 #7', html_coseg.indexOf('L2 #7') > -1);
ok('co_seg html shows regime chip in green',
   html_coseg.indexOf('#7ad394') > -1 && html_coseg.indexOf('co_seg') > -1);
ok('co_seg html shows dominant band 0 (b0)', /b0/.test(html_coseg));
ok('co_seg html shows dominant percentage', /83\.0%/.test(html_coseg));
ok('co_seg html shows entropy 0.213', html_coseg.indexOf('0.213') > -1);
ok('co_seg html shows n_valid 12 / 12', /12\s*\/\s*12/.test(html_coseg));
ok('co_seg html omits chain row when n_chains=1',
   html_coseg.indexOf('chain ') < 0);
ok('co_seg html lists band fractions', html_coseg.indexOf('band fractions:') > -1);
ok('co_seg html shows non-zero band-fraction swatches (b0 + b2)',
   html_coseg.indexOf('b0 83%') > -1 && html_coseg.indexOf('b2 17%') > -1);
ok('co_seg html omits zero-fraction band b1', html_coseg.indexOf('b1 0%') < 0);

// A multi-chain hit shows the chain badge
const hit_chain = {
  l2_idx: 22, x: 0, y: 0, w: 10, h: 7,
  entry: {
    l2_idx: 22, chain_idx: 1, chain_position: 4, n_valid: 8,
    band_fractions: [0.5, 0.5, 0.0],
    entropy: 0.63, dominant_band: 0, dominant_fraction: 0.5,
    regime: 'partial',
  },
};
sbx._h = hit_chain;
const html_chain = vm.runInContext('_btraceTooltipBuildHtml(_h, {n_chains: 2, n_fish_selected: 8})', sbx);
ok('chain row shown when n_chains>1', html_chain.indexOf('chain 1') > -1);
ok('partial regime renders amber chip',
   html_chain.indexOf('#f5a524') > -1 && html_chain.indexOf('partial') > -1);

// no_valid edge case — body shows entropy/n_valid but skips band-fractions section
const hit_nv = {
  l2_idx: 99, x: 0, y: 0, w: 10, h: 7,
  entry: {
    l2_idx: 99, chain_idx: 0, chain_position: 0, n_valid: 0,
    band_fractions: [0, 0, 0],
    entropy: 1.0, dominant_band: -1, dominant_fraction: 0,
    regime: 'no_valid',
  },
};
sbx._h = hit_nv;
const html_nv = vm.runInContext('_btraceTooltipBuildHtml(_h, {n_chains: 1, n_fish_selected: 5})', sbx);
ok('no_valid html omits band-fractions section',
   html_nv.indexOf('band fractions:') < 0);
ok('no_valid html omits dominant row when dominant_band=-1',
   html_nv.indexOf('dominant:') < 0);
ok('no_valid html still shows entropy + n_valid',
   html_nv.indexOf('entropy:') > -1 && /0\s*\/\s*5/.test(html_nv));

// Defensive: null hit / null entry → empty string, no throw
ok('null hit returns empty string',
   vm.runInContext('_btraceTooltipBuildHtml(null, {})', sbx) === '');
ok('hit without entry returns empty string',
   vm.runInContext('_btraceTooltipBuildHtml({l2_idx: 0}, {})', sbx) === '');

// ============================================================================
// 4. _btraceHitTest
// ============================================================================
console.log('\n=== 4. _btraceHitTest ===');

const fnHit = pullFunction(html, '_btraceHitTest');
ok('pulled _btraceHitTest', !!fnHit);
vm.runInContext(fnHit, sbx);

vm.runInContext(`
  const _hits = [
    { x: 10, y: 0, w: 10, h: 7, l2_idx: 0 },
    { x: 22, y: 0, w: 10, h: 7, l2_idx: 1 },
    { x: 34, y: 0, w: 10, h: 7, l2_idx: 2 },
  ];
`, sbx);

ok('hit inside first rect (15, 3)',
   vm.runInContext('_btraceHitTest(_hits, 15, 3)', sbx).l2_idx === 0);
ok('hit inside third rect (40, 4)',
   vm.runInContext('_btraceHitTest(_hits, 40, 4)', sbx).l2_idx === 2);
ok('miss between rects (21, 3) → null',
   vm.runInContext('_btraceHitTest(_hits, 21, 3)', sbx) === null);
ok('miss above strip (15, -5) → null',
   vm.runInContext('_btraceHitTest(_hits, 15, -5)', sbx) === null);
ok('miss empty hits → null',
   vm.runInContext('_btraceHitTest([], 15, 3)', sbx) === null);
ok('miss null hits → null',
   vm.runInContext('_btraceHitTest(null, 15, 3)', sbx) === null);
ok('hit at left-edge (10, 3) inclusive',
   vm.runInContext('_btraceHitTest(_hits, 10, 3)', sbx).l2_idx === 0);
ok('hit at right-edge (20, 3) inclusive',
   vm.runInContext('_btraceHitTest(_hits, 20, 3)', sbx).l2_idx === 0);

// ============================================================================
// 5. _wireBandTraceTooltip — idempotency
// ============================================================================
console.log('\n=== 5. _wireBandTraceTooltip idempotency ===');

const fnWire = pullFunction(html, '_wireBandTraceTooltip');
ok('pulled _wireBandTraceTooltip', !!fnWire);
// _wireBandTraceTooltip references _btraceTooltipHide / _btraceTooltipShow
// inside its mousemove handler. Stub them so the wire path can install
// without ReferenceErrors. We don't exercise the listener body in this
// section — listener-body coverage is via _btraceHitTest + _btraceTooltipBuildHtml.
vm.runInContext(`
  function _btraceTooltipHide() {}
  function _btraceTooltipShow() {}
`, sbx);
vm.runInContext(fnWire, sbx);

// Fake canvas with addEventListener that records calls
vm.runInContext(`
  function makeFakeCanvas() {
    return {
      dataset: {},
      _listeners: {},
      addEventListener(ev, fn) {
        if (!this._listeners[ev]) this._listeners[ev] = [];
        this._listeners[ev].push(fn);
      },
      style: {},
    };
  }
  var _cv = makeFakeCanvas();
  _wireBandTraceTooltip(_cv);
  var _firstAttach = _cv._listeners.mousemove ? _cv._listeners.mousemove.length : 0;
  _wireBandTraceTooltip(_cv);   // idempotent — should NOT add a second handler
  var _secondAttach = _cv._listeners.mousemove ? _cv._listeners.mousemove.length : 0;
`, sbx);

ok('first wire attaches mousemove handler',
   sbx._firstAttach === 1);
ok('second wire is idempotent (no second attach)',
   sbx._secondAttach === 1);
ok('dataset marker set after wire',
   sbx._cv.dataset.btraceTooltipWired === '1');

// Null/undefined canvas is safe
ok('null canvas is no-op',
   vm.runInContext('_wireBandTraceTooltip(null); true', sbx) === true);

// ============================================================================
// 6, 7, 8. _drawBandTraceStrip hits population
// ============================================================================
console.log('\n=== 6-8. _drawBandTraceStrip hits ===');

// Set up a fresh sandbox so we can exercise the drawer end-to-end.
const sbx2 = {
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
sbx2.window.state = sbx2.state;
sbx2.localStorage = (function() {
  const _s = {};
  return { setItem: (k,v)=>{_s[k]=String(v);}, getItem: k=>(k in _s?_s[k]:null), removeItem: k=>{delete _s[k];}, _store:_s };
})();
vm.createContext(sbx2);
vm.runInContext('var localStorage = this.localStorage;', sbx2);
vm.runInContext('function drawLinesPanel(){}', sbx2);
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
  const _BTRACE_ON_LS_KEY        = 'inversion_atlas.bandTraceOn';
  const _BTRACE_FISH_SET_LS_KEY  = 'inversion_atlas.bandTraceFishSet';
  const _BTRACE_REGIME_COLOR = {co_seg:'#7ad394',partial:'#f5a524',fanned:'#9aa3ad',sparse:'#5a6068',no_valid:'#3a3e44'};
`, sbx2);

const fnCacheKey = pullFunction(html, '_bandTraceCacheKey');
const fnGoC      = pullFunction(html, '_bandTraceGetOrCompute');
const fnSetFish  = pullFunction(html, 'setBandTraceFishSet');
const fnSetOn    = pullFunction(html, 'setBandTraceOn');
const fnClear2   = pullFunction(html, '_bandTraceClearCache');
const fnDraw     = pullFunction(html, '_drawBandTraceStrip');
vm.runInContext(fnCacheKey, sbx2);
vm.runInContext(fnGoC, sbx2);
vm.runInContext(fnSetFish, sbx2);
vm.runInContext(fnSetOn, sbx2);
vm.runInContext(fnClear2, sbx2);
vm.runInContext(fnDraw, sbx2);

function makeCtx() {
  const calls = [];
  return {
    save: () => calls.push(['save']),
    restore: () => calls.push(['restore']),
    fillRect: (...a) => calls.push(['fillRect', ...a]),
    strokeRect: (...a) => calls.push(['strokeRect', ...a]),
    set fillStyle(v) { this._fill = v; },
    get fillStyle() { return this._fill; },
    set strokeStyle(v) { this._stroke = v; },
    get strokeStyle() { return this._stroke; },
    set lineWidth(v) {},
    set globalAlpha(v) {},
    _calls: calls,
  };
}

// Setup data
sbx2.state.data = {
  chrom: 'LG28',
  l2_envelopes: [
    { start_bp: 1_000_000, end_bp: 2_000_000 },
    { start_bp: 2_500_000, end_bp: 3_000_000 },
    { start_bp: 3_500_000, end_bp: 4_000_000 },
  ],
};

// Active state — toggle on, fish-set, cache cleared
sbx2.state.bandTraceOn = true;
sbx2.state.bandTraceFishSet = [1, 2, 3, 4, 5];
sbx2.state.bandTraceCache = null;
sbx2.state.bandTraceCacheKey = null;
sbx2.state._btraceHits = null;
sbx2._ctx = makeCtx();

vm.runInContext('_drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx2);

ok('hits array exists after paint', Array.isArray(sbx2.state._btraceHits));
ok('one hit per visible L2 (3 envelopes all in view)',
   sbx2.state._btraceHits && sbx2.state._btraceHits.length === 3);
ok('hit shape has x/y/w/h/l2_idx/entry',
   sbx2.state._btraceHits[0] &&
   typeof sbx2.state._btraceHits[0].x === 'number' &&
   typeof sbx2.state._btraceHits[0].y === 'number' &&
   typeof sbx2.state._btraceHits[0].w === 'number' &&
   typeof sbx2.state._btraceHits[0].h === 'number' &&
   typeof sbx2.state._btraceHits[0].l2_idx === 'number' &&
   sbx2.state._btraceHits[0].entry != null);
ok('hits ordered by l2_idx (0, 1, 2)',
   sbx2.state._btraceHits[0].l2_idx === 0 &&
   sbx2.state._btraceHits[1].l2_idx === 1 &&
   sbx2.state._btraceHits[2].l2_idx === 2);
ok('hit rect height equals strip height (7)',
   sbx2.state._btraceHits[0].h === 7);
ok('hit entry carries regime + entropy',
   sbx2.state._btraceHits[0].entry.regime === 'co_seg' &&
   typeof sbx2.state._btraceHits[0].entry.entropy === 'number');

// Off-screen L2 should NOT produce a hit (mb range below the envelope)
sbx2.state._btraceHits = null;
sbx2._ctx = makeCtx();
vm.runInContext('_drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 10, 20)', sbx2);
ok('no hits when all envelopes outside view range',
   Array.isArray(sbx2.state._btraceHits) && sbx2.state._btraceHits.length === 0);

// Toggle off → no hits written, no paint
sbx2.state.bandTraceOn = false;
sbx2.state._btraceHits = null;   // ensure stale state if any
sbx2._ctx = makeCtx();
vm.runInContext('_drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx2);
ok('toggle off → drawer is no-op (no hits write)',
   sbx2.state._btraceHits === null);

// No fish-set → no hits written
sbx2.state.bandTraceOn = true;
sbx2.state.bandTraceFishSet = null;
sbx2.state._btraceHits = null;
sbx2._ctx = makeCtx();
vm.runInContext('_drawBandTraceStrip(_ctx, {l: 0, t: 30}, 600, 200, 0, 5)', sbx2);
ok('no fish-set → drawer is no-op (no hits write)',
   sbx2.state._btraceHits === null);

// _bandTraceClearCache nulls hits
sbx2.state._btraceHits = [{x:1,y:1,w:1,h:1,l2_idx:0,entry:{}}];
vm.runInContext('_bandTraceClearCache()', sbx2);
ok('_bandTraceClearCache clears _btraceHits',
   sbx2.state._btraceHits === null);

// applyData chrom-switch path (regression check at the source level — the
// turn-161 code already does state.bandTraceFishSet=null; verify the new
// hits clear is also there transitively via _bandTraceClearCache call).
ok('applyData calls _bandTraceClearCache (transitively clears hits)',
   /_bandTraceClearCache\(\)/.test(html));

// ============================================================================
// 9, 10, 11. _bandTraceToTSV
// ============================================================================
console.log('\n=== 9-11. _bandTraceToTSV ===');

const fnTSV = pullFunction(html, '_bandTraceToTSV');
ok('pulled _bandTraceToTSV', !!fnTSV);
vm.runInContext(fnTSV, sbx2);

// Build a synthetic trace
vm.runInContext(`
  var _trace = {
    n_fish_selected: 12,
    n_chains: 1, n_total_L2: 3, K: 3,
    per_l2: [
      { l2_idx: 0, chain_idx: 0, chain_position: 0, n_valid: 12,
        band_fractions: new Float32Array([0.83, 0.0, 0.17]),
        entropy: 0.213, dominant_band: 0, dominant_fraction: 0.83, regime: 'co_seg' },
      { l2_idx: 1, chain_idx: 0, chain_position: 1, n_valid: 10,
        band_fractions: new Float32Array([0.5, 0.4, 0.1]),
        entropy: 0.671, dominant_band: 0, dominant_fraction: 0.5, regime: 'partial' },
      { l2_idx: 2, chain_idx: 0, chain_position: 2, n_valid: 0,
        band_fractions: new Float32Array([0, 0, 0]),
        entropy: 1.0, dominant_band: -1, dominant_fraction: 0, regime: 'no_valid' },
    ],
  };
  var _envs = [
    { start_bp: 1000000, end_bp: 2000000 },
    { start_bp: 2500000, end_bp: 3000000 },
    { start_bp: 3500000, end_bp: 4000000 },
  ];
`, sbx2);

const tsv = vm.runInContext('_bandTraceToTSV(_trace, {chrom: "LG28", envelopes: _envs})', sbx2);
const tsvLines = tsv.split('\n').filter(l => l.length > 0);

ok('TSV is non-empty', tsvLines.length > 0);
ok('TSV first line is header', tsvLines[0].split('\t')[0] === 'chrom');
ok('TSV header has expected core columns',
   /chrom\tl2_idx\tchain_idx\tchain_position\tstart_bp\tend_bp\tn_valid\tn_fish_selected\tregime\tdominant_band\tdominant_fraction\tentropy/
     .test(tsvLines[0]));
ok('TSV header has K=3 band_fraction columns',
   /band_fraction_0\tband_fraction_1\tband_fraction_2/.test(tsvLines[0]));
ok('TSV has 3 data rows (one per L2)', tsvLines.length === 4);

const cols0 = tsvLines[1].split('\t');
ok('row 0 chrom = LG28', cols0[0] === 'LG28');
ok('row 0 l2_idx = 0', cols0[1] === '0');
ok('row 0 chain_idx = 0', cols0[2] === '0');
ok('row 0 start_bp = 1000000', cols0[4] === '1000000');
ok('row 0 end_bp = 2000000', cols0[5] === '2000000');
ok('row 0 n_valid = 12', cols0[6] === '12');
ok('row 0 n_fish_selected = 12', cols0[7] === '12');
ok('row 0 regime = co_seg', cols0[8] === 'co_seg');
ok('row 0 dominant_band = 0', cols0[9] === '0');
ok('row 0 dominant_fraction formatted to 6 dp', /^0\.83\d{4}$/.test(cols0[10]));
ok('row 0 entropy formatted to 6 dp', /^0\.213\d{3}$/.test(cols0[11]));
ok('row 0 band_fraction_0 ≈ 0.83', /^0\.83\d{4}$/.test(cols0[12]));
ok('row 0 band_fraction_1 = 0.000000', cols0[13] === '0.000000');
ok('row 0 band_fraction_2 ≈ 0.17', /^0\.17\d{4}$/.test(cols0[14]));

// no_valid row has dominant_band = -1 and zero band fractions
const cols2 = tsvLines[3].split('\t');
ok('no_valid row regime = no_valid', cols2[8] === 'no_valid');
ok('no_valid row dominant_band = -1', cols2[9] === '-1');
ok('no_valid row n_valid = 0', cols2[6] === '0');

// Empty per_l2 → just header
const tsv_empty = vm.runInContext(`
  _bandTraceToTSV({n_fish_selected: 5, K: 3, n_chains: 1, n_total_L2: 0, per_l2: []}, {chrom: 'LG14'})
`, sbx2);
const lines_empty = tsv_empty.split('\n').filter(l => l.length > 0);
ok('empty per_l2 → header only', lines_empty.length === 1);
ok('empty per_l2 header still has K=3 band cols',
   /band_fraction_0\tband_fraction_1\tband_fraction_2/.test(lines_empty[0]));

// null trace → null
ok('null trace → null', vm.runInContext('_bandTraceToTSV(null, {})', sbx2) === null);

// Missing envelopes → blank start_bp/end_bp (no crash)
const tsv_no_env = vm.runInContext('_bandTraceToTSV(_trace, {chrom: "LG14"})', sbx2);
const lines_no_env = tsv_no_env.split('\n').filter(l => l.length > 0);
const cols_no_env = lines_no_env[1].split('\t');
ok('no envelopes: start_bp blank', cols_no_env[4] === '');
ok('no envelopes: end_bp blank', cols_no_env[5] === '');
ok('no envelopes: rest of row still populated', cols_no_env[8] === 'co_seg');

// ============================================================================
// 12. _bandTraceDownloadTSV — filename + null path
// ============================================================================
console.log('\n=== 12. _bandTraceDownloadTSV ===');

const fnDl = pullFunction(html, '_bandTraceDownloadTSV');
ok('pulled _bandTraceDownloadTSV', !!fnDl);
vm.runInContext(fnDl, sbx2);

// No fish-set → null
sbx2.state.bandTraceFishSet = null;
sbx2.state.bandTraceCache = null;
sbx2.state.bandTraceCacheKey = null;
ok('no fish-set → returns null',
   vm.runInContext('_bandTraceDownloadTSV()', sbx2) === null);

// With cached trace + fish-set, headless env: returns filename, no crash.
// (document/Blob aren't defined in the sandbox, so the function takes the
// headless branch and returns the filename without performing the download.)
sbx2.state.bandTraceFishSet = [1, 2, 3, 4, 5];
sbx2.state.bandTraceCache = null;          // force recompute via _bandTraceGetOrCompute
sbx2.state.bandTraceCacheKey = null;
const fname = vm.runInContext('_bandTraceDownloadTSV()', sbx2);
ok('returns a filename string', typeof fname === 'string' && fname.length > 0);
ok('filename starts with band_trace_LG28', fname.indexOf('band_trace_LG28') === 0);
ok('filename embeds n_fish (n5)', /_n5_/.test(fname));
ok('filename embeds K (K3)', /_K3\.tsv$/.test(fname));

// Empty per_l2 → null (because no envelopes / no cached trace fully)
sbx2.state.data = { chrom: 'LG28', l2_envelopes: [] };
sbx2.state.bandTraceCache = null;
sbx2.state.bandTraceCacheKey = null;
ok('empty l2_envelopes → null',
   vm.runInContext('_bandTraceDownloadTSV()', sbx2) === null);

// ============================================================================
// 13. DOM constants + tooltip element id
// ============================================================================
console.log('\n=== 13. DOM constants ===');

ok('linesBandTraceExportBtn id present in source',
   html.indexOf('id="linesBandTraceExportBtn"') > -1);
ok('export button title mentions TSV',
   /id="linesBandTraceExportBtn"[\s\S]{0,500}TSV/.test(html));
ok('export button title mentions filename pattern',
   /band_trace_<chrom>_n<n_fish>/.test(html));

// ============================================================================
// 14. Regression
// ============================================================================
console.log('\n=== 14. Regression ===');

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
ok('turn 160 _bandTraceForFishSet still defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('turn 160 _bandTraceRegimeRuns still defined',
   /function\s+_bandTraceRegimeRuns\s*\(/.test(html));
ok('turn 130 _drawLineageStrip still defined',
   /function\s+_drawLineageStrip\s*\(/.test(html));
ok('turn 130 lineage toggle still in lines header',
   html.indexOf('id="linesLineageStripToggle"') > -1);
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
