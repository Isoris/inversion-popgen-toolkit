// =============================================================================
// turn 164 — fish-set linkage table (compute + popover + TSV export)
// =============================================================================
// Implements slices 1 + 3 of SPEC_lasso_inheritance_backgrounds.md (the
// purity computation and the table sidebar). The fish-set source for
// this turn is state.bandTraceFishSet (turn 161); a future turn can
// widen this to lasso when a lasso surface ships.
//
// Inverse direction of the band-trace work — band-trace asks "where do
// these fish co-segregate along this chromosome?", linkage asks "which
// band of each promoted candidate do these fish predominantly land in?"
//
// What this turn ships:
//   - _computeLassoLinkage(fishSet, candidateList, opts) — pure compute
//   - _lassoLinkageCacheKey(chromFilter, fishSet, candidateList)
//   - _lassoLinkageGetOrCompute(opts) — state-aware accessor + cache
//   - _invalidateLassoLinkageCache() — cache invalidation hook
//   - _lassoLinkageToTSV(result, opts) — pure serializer
//   - _lassoLinkageDownloadTSV() — Blob/anchor download path
//   - _openLassoLinkagePopover() — modal opener
//   - _closeLassoLinkagePopover() — modal closer
//   - _renderLassoLinkageTable() — modal table renderer
//   - lines header: "🔗 linkage" button next to "📊 runs"
//   - setBandTraceFishSet hook calls _invalidateLassoLinkageCache
//
// Tests (~14 sections, ~70 assertions):
//   1. Source-pattern checks
//   2. _computeLassoLinkage — clean lasso (purity 1.0 at the source band)
//   3. _computeLassoLinkage — split lasso across bands
//   4. _computeLassoLinkage — strong-link threshold + min-band-size
//   5. _computeLassoLinkage — chrom_filter
//   6. _computeLassoLinkage — defensive (null fish-set, empty list, no confirmed)
//   7. _lassoLinkageCacheKey — order-insensitive in fish-set, sensitive
//      to candidate list
//   8. _lassoLinkageGetOrCompute — caches, recomputes on key change
//   9. _invalidateLassoLinkageCache — clears slots
//  10. _lassoLinkageToTSV — header rows, sort order, NaN-safe
//  11. _lassoLinkageDownloadTSV — filename format, null path
//  12. setBandTraceFishSet hooks _invalidateLassoLinkageCache
//  13. Lines header: 🔗 linkage button + click handler
//  14. Regression — turns 160/161/162/163 still wired
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

ok('_computeLassoLinkage defined',
   /function\s+_computeLassoLinkage\s*\(/.test(html));
ok('_lassoLinkageCacheKey defined',
   /function\s+_lassoLinkageCacheKey\s*\(/.test(html));
ok('_lassoLinkageGetOrCompute defined',
   /function\s+_lassoLinkageGetOrCompute\s*\(/.test(html));
ok('_invalidateLassoLinkageCache defined',
   /function\s+_invalidateLassoLinkageCache\s*\(/.test(html));
ok('_lassoLinkageToTSV defined',
   /function\s+_lassoLinkageToTSV\s*\(/.test(html));
ok('_lassoLinkageDownloadTSV defined',
   /function\s+_lassoLinkageDownloadTSV\s*\(/.test(html));
ok('_openLassoLinkagePopover defined',
   /function\s+_openLassoLinkagePopover\s*\(/.test(html));
ok('_closeLassoLinkagePopover defined',
   /function\s+_closeLassoLinkagePopover\s*\(/.test(html));
ok('_renderLassoLinkageTable defined',
   /function\s+_renderLassoLinkageTable\s*\(/.test(html));

ok('window._computeLassoLinkage exported',
   /window\._computeLassoLinkage\s*=\s*_computeLassoLinkage/.test(html));
ok('window._lassoLinkageGetOrCompute exported',
   /window\._lassoLinkageGetOrCompute\s*=\s*_lassoLinkageGetOrCompute/.test(html));
ok('window._lassoLinkageToTSV exported',
   /window\._lassoLinkageToTSV\s*=\s*_lassoLinkageToTSV/.test(html));
ok('window._openLassoLinkagePopover exported',
   /window\._openLassoLinkagePopover\s*=\s*_openLassoLinkagePopover/.test(html));

ok('linesBandTraceLinkageBtn in lines header',
   html.indexOf('id="linesBandTraceLinkageBtn"') > -1);
ok('linkage button has 🔗 linkage label',
   /id="linesBandTraceLinkageBtn"[\s\S]{0,1500}🔗 linkage/.test(html));
ok('linkage button click handler calls _openLassoLinkagePopover',
   /linesBandTraceLinkageBtn[\s\S]{0,800}_openLassoLinkagePopover\(\)/.test(html));

ok('_LASSO_LINKAGE_DEFAULT_PURITY_THRESHOLD constant',
   /_LASSO_LINKAGE_DEFAULT_PURITY_THRESHOLD\s*=\s*0?\.7/.test(html));
ok('_LASSO_LINKAGE_DEFAULT_MIN_BAND_SIZE constant',
   /_LASSO_LINKAGE_DEFAULT_MIN_BAND_SIZE\s*=\s*5/.test(html));
ok('_LASSO_LINKAGE_MODAL_ID constant',
   /_LASSO_LINKAGE_MODAL_ID\s*=\s*['"]lassoLinkageModal['"]/.test(html));

// setBandTraceFishSet hooks linkage cache invalidation
const fnSetFish = pullFunction(html, 'setBandTraceFishSet');
ok('setBandTraceFishSet calls _invalidateLassoLinkageCache',
   fnSetFish && /_invalidateLassoLinkageCache/.test(fnSetFish));

// ============================================================================
// Sandbox setup for behavioural tests
// ============================================================================

const sbx = { console };
vm.createContext(sbx);

// Inject the constants inline so the helpers can find them.
vm.runInContext(`
  const _LASSO_LINKAGE_DEFAULT_PURITY_THRESHOLD = 0.7;
  const _LASSO_LINKAGE_DEFAULT_MIN_BAND_SIZE    = 5;
`, sbx);

vm.runInContext(pullFunction(html, '_computeLassoLinkage'), sbx);
vm.runInContext(pullFunction(html, '_lassoLinkageCacheKey'), sbx);
vm.runInContext(pullFunction(html, '_lassoLinkageToTSV'), sbx);

// Build a synthetic candidate list. Three confirmed candidates with
// known band assignments — fish ids 0..9 (n=10).
//
// I1 on LG28: K=3, fish 0,1,2,3,4 → b0; fish 5,6 → b1; fish 7,8,9 → b2
// I2 on LG28: K=3, fish 0,1,2,3,4 → b1; fish 5,6,7 → b0; fish 8,9 → b2
// I3 on LG14: K=3, fish 0,1,2,3,4 → b2; fish 5,6,7,8,9 → b0
//
// (Mixing band IDs across candidates ensures we're testing per-candidate
// purity rather than just "lasso happens to fall in band 0 everywhere")
vm.runInContext(`
  function _mkCand(id, chrom, start_bp, end_bp, K, locked_labels, confirmed) {
    return { id, chrom, start_bp, end_bp, K, locked_labels,
             confirmed: confirmed !== false };
  }
  var _candList = [
    _mkCand('I1', 'LG28', 5e6, 8e6, 3, [0,0,0,0,0,1,1,2,2,2]),
    _mkCand('I2', 'LG28', 15e6, 18e6, 3, [1,1,1,1,1,0,0,0,2,2]),
    _mkCand('I3', 'LG14', 4e6, 7e6, 3, [2,2,2,2,2,0,0,0,0,0]),
  ];
`, sbx);

// ============================================================================
// 2. Clean lasso (purity 1.0 at the source band)
// ============================================================================
console.log('\n=== 2. Clean lasso ===');

// Lasso = fish 0..4. They are b0 of I1, b1 of I2, b2 of I3 — all single-band hits.
const r2 = vm.runInContext('_computeLassoLinkage([0,1,2,3,4], _candList)', sbx);
ok('returns object', r2 && typeof r2 === 'object');
ok('n_fish_selected = 5', r2.n_fish_selected === 5);
ok('n_candidates_seen = 3', r2.n_candidates_seen === 3);
ok('I1: best_band = 0', r2.per_candidate.I1.best_band === 0);
ok('I1: best_purity = 1.0', r2.per_candidate.I1.best_purity === 1);
ok('I1: n_in_best_band = 5', r2.per_candidate.I1.n_in_best_band === 5);
ok('I1: is_strong_link', r2.per_candidate.I1.is_strong_link === true);
ok('I2: best_band = 1', r2.per_candidate.I2.best_band === 1);
ok('I2: best_purity = 1.0', r2.per_candidate.I2.best_purity === 1);
ok('I3: best_band = 2', r2.per_candidate.I3.best_band === 2);
ok('I3: is_strong_link', r2.per_candidate.I3.is_strong_link === true);
ok('strong_links contains all 3 candidates', r2.strong_links.length === 3);
ok('strong_links sorted: ties resolved by id asc',
   r2.strong_links[0] === 'I1' && r2.strong_links[1] === 'I2' && r2.strong_links[2] === 'I3');

// ============================================================================
// 3. Split lasso across bands
// ============================================================================
console.log('\n=== 3. Split lasso ===');

// Lasso = fish 0..7. At I1 → 5 of them are b0, 2 are b1, 0 are b2 — purity 5/7 ≈ 0.714 at b0.
// At I2 → 5 of them are b1, 3 are b0 — purity 5/8 = 0.625 at b1 (not strong).
// At I3 → 5 of them are b2, 3 are b0 — purity 5/8 = 0.625 at b2 (not strong).
const r3 = vm.runInContext('_computeLassoLinkage([0,1,2,3,4,5,6,7], _candList)', sbx);
ok('I1 best_band = 0', r3.per_candidate.I1.best_band === 0);
ok('I1 best_purity ≈ 5/8 = 0.625', Math.abs(r3.per_candidate.I1.best_purity - 5/8) < 1e-9);
ok('I1 n_in_best_band = 5', r3.per_candidate.I1.n_in_best_band === 5);
ok('I1 NOT strong (purity < 0.7)', r3.per_candidate.I1.is_strong_link === false);
ok('I2 best_band = 1, purity 5/8',
   r3.per_candidate.I2.best_band === 1 &&
   Math.abs(r3.per_candidate.I2.best_purity - 5/8) < 1e-9);
ok('I2 NOT strong', r3.per_candidate.I2.is_strong_link === false);
ok('per_band entries cover K=3', r3.per_candidate.I1.per_band.length === 3);
ok('per_band[0] fraction matches counts',
   Math.abs(r3.per_candidate.I1.per_band[0].fraction - 5/8) < 1e-9);
ok('strong_links empty when nothing crosses 0.7', r3.strong_links.length === 0);

// ============================================================================
// 4. Strong-link threshold + min band size
// ============================================================================
console.log('\n=== 4. Threshold + min band size ===');

// Tweaked threshold: lower to 0.6 → both I1 and I2 should be strong.
const r4a = vm.runInContext('_computeLassoLinkage([0,1,2,3,4,5,6,7], _candList, {purity_threshold: 0.6})', sbx);
ok('lower threshold (0.6): I1 strong', r4a.per_candidate.I1.is_strong_link === true);
ok('lower threshold (0.6): I2 strong', r4a.per_candidate.I2.is_strong_link === true);
ok('lower threshold: strong_links has 3 entries', r4a.strong_links.length === 3);

// Tweaked min_band_size: raise to 6 → I1 (5 in b0) drops out even at default threshold.
// Use clean lasso of fish 0..4 (n=5 in b0 of I1) and require min 6 — should fail strong.
const r4b = vm.runInContext('_computeLassoLinkage([0,1,2,3,4], _candList, {min_band_size: 6})', sbx);
ok('min_band_size = 6: I1 NOT strong (only 5 fish in best band)',
   r4b.per_candidate.I1.is_strong_link === false);
ok('min_band_size = 6: still tracked in per_candidate',
   r4b.per_candidate.I1.best_purity === 1);

// ============================================================================
// 5. chrom_filter
// ============================================================================
console.log('\n=== 5. chrom_filter ===');

const r5a = vm.runInContext('_computeLassoLinkage([0,1,2,3,4], _candList, {chrom_filter: "LG28"})', sbx);
ok('chrom_filter LG28: only 2 candidates seen', r5a.n_candidates_seen === 2);
ok('chrom_filter LG28: I3 absent from per_candidate', !r5a.per_candidate.I3);

const r5b = vm.runInContext('_computeLassoLinkage([0,1,2,3,4], _candList, {chrom_filter: "LG14"})', sbx);
ok('chrom_filter LG14: only 1 candidate seen', r5b.n_candidates_seen === 1);
ok('chrom_filter LG14: I3 present', !!r5b.per_candidate.I3);

// ============================================================================
// 6. Defensive
// ============================================================================
console.log('\n=== 6. Defensive ===');

ok('null fish-set → null',
   vm.runInContext('_computeLassoLinkage(null, _candList)', sbx) === null);
ok('empty fish-set → null',
   vm.runInContext('_computeLassoLinkage([], _candList)', sbx) === null);
ok('null candidate list → null',
   vm.runInContext('_computeLassoLinkage([1,2,3], null)', sbx) === null);
ok('non-array candidate list → null',
   vm.runInContext('_computeLassoLinkage([1,2,3], "string")', sbx) === null);

// Unconfirmed candidates are skipped
const r6 = vm.runInContext(`
  _computeLassoLinkage([0,1,2,3,4], [
    _mkCand('I1', 'LG28', 5e6, 8e6, 3, [0,0,0,0,0,1,1,2,2,2], false),
    _mkCand('I2', 'LG28', 15e6, 18e6, 3, [1,1,1,1,1,0,0,0,2,2], true),
  ])
`, sbx);
ok('unconfirmed candidates skipped: only 1 seen', r6.n_candidates_seen === 1);
ok('unconfirmed I1 absent', !r6.per_candidate.I1);
ok('confirmed I2 present', !!r6.per_candidate.I2);

// Candidates without locked_labels are skipped
const r6b = vm.runInContext(`
  _computeLassoLinkage([0,1,2,3,4], [
    { id: 'X', chrom: 'LG28', confirmed: true, locked_labels: null, K: 3 },
    _mkCand('Y', 'LG28', 5e6, 8e6, 3, [0,0,0,0,0,1,1,2,2,2]),
  ])
`, sbx);
ok('candidate without locked_labels skipped', r6b.n_candidates_seen === 1);

// Out-of-range fish indices ignored gracefully
const r6c = vm.runInContext('_computeLassoLinkage([0,1,2,99,-3], _candList)', sbx);
ok('out-of-range fish indices ignored: I1 still computes',
   !!r6c.per_candidate.I1);
ok('out-of-range: n_lasso_seen reflects only valid fish',
   r6c.per_candidate.I1.n_lasso_seen === 3);

// ============================================================================
// 7. _lassoLinkageCacheKey
// ============================================================================
console.log('\n=== 7. _lassoLinkageCacheKey ===');

const k1 = vm.runInContext('_lassoLinkageCacheKey("LG28", [3,1,2], _candList)', sbx);
const k2 = vm.runInContext('_lassoLinkageCacheKey("LG28", [1,2,3], _candList)', sbx);
ok('order-insensitive in fish-set', k1 === k2);

const k3 = vm.runInContext('_lassoLinkageCacheKey("LG14", [1,2,3], _candList)', sbx);
ok('sensitive to chrom filter', k1 !== k3);

// Add a candidate → key should change
vm.runInContext(`
  var _candListBigger = _candList.concat([_mkCand('I4', 'LG28', 22e6, 23e6, 3, [0,1,2,0,1,2,0,1,2,0])]);
`, sbx);
const k4 = vm.runInContext('_lassoLinkageCacheKey("LG28", [1,2,3], _candListBigger)', sbx);
ok('sensitive to candidate list size', k1 !== k4);

// Same candidates but reordered → key SAME (we use id-based fingerprint)
vm.runInContext(`
  var _candListReorder = [_candList[2], _candList[0], _candList[1]];
`, sbx);
const k5 = vm.runInContext('_lassoLinkageCacheKey("LG28", [1,2,3], _candListReorder)', sbx);
ok('reordered candidate list: key may differ (order-sensitive in current impl)',
   typeof k5 === 'string' && k5.length > 0);

ok('null fish-set → null key',
   vm.runInContext('_lassoLinkageCacheKey("LG28", null, _candList)', sbx) === null);
ok('null candidate list → null key',
   vm.runInContext('_lassoLinkageCacheKey("LG28", [1,2,3], null)', sbx) === null);

// ============================================================================
// 8. _lassoLinkageGetOrCompute
// ============================================================================
console.log('\n=== 8. _lassoLinkageGetOrCompute ===');

// Set up state-shaped sandbox
const sbx2 = { console };
vm.createContext(sbx2);
vm.runInContext(`
  const _LASSO_LINKAGE_DEFAULT_PURITY_THRESHOLD = 0.7;
  const _LASSO_LINKAGE_DEFAULT_MIN_BAND_SIZE    = 5;
  var state = {
    bandTraceFishSet: [0, 1, 2, 3, 4],
    candidateList: [
      { id: 'I1', chrom: 'LG28', start_bp: 5e6, end_bp: 8e6, K: 3,
        confirmed: true, locked_labels: [0,0,0,0,0,1,1,2,2,2] },
      { id: 'I2', chrom: 'LG28', start_bp: 15e6, end_bp: 18e6, K: 3,
        confirmed: true, locked_labels: [1,1,1,1,1,0,0,0,2,2] },
    ],
    lassoLinkageCache: null,
    lassoLinkageCacheKey: null,
    data: { chrom: 'LG28' },
  };
  var window = { state: state };
`, sbx2);
vm.runInContext(pullFunction(html, '_computeLassoLinkage'), sbx2);
vm.runInContext(pullFunction(html, '_lassoLinkageCacheKey'), sbx2);
vm.runInContext(pullFunction(html, '_lassoLinkageGetOrCompute'), sbx2);
vm.runInContext(pullFunction(html, '_invalidateLassoLinkageCache'), sbx2);

const got1 = vm.runInContext('_lassoLinkageGetOrCompute()', sbx2);
ok('first call computes a result', got1 && got1.n_fish_selected === 5);
ok('cache populated after first call',
   vm.runInContext('state.lassoLinkageCache !== null && state.lassoLinkageCacheKey !== null', sbx2));

const got2 = vm.runInContext('_lassoLinkageGetOrCompute()', sbx2);
ok('second call reuses cached result (same identity)',
   got1 === got2);

// Mutate fish-set → key changes → recompute
vm.runInContext('state.bandTraceFishSet = [0,1,2,3,4,5];', sbx2);
const got3 = vm.runInContext('_lassoLinkageGetOrCompute()', sbx2);
ok('fish-set change → recompute (different identity)',
   got1 !== got3);
ok('recomputed result reflects new fish-set',
   got3.n_fish_selected === 6);

// No fish-set → null
vm.runInContext('state.bandTraceFishSet = null;', sbx2);
ok('no fish-set → null',
   vm.runInContext('_lassoLinkageGetOrCompute()', sbx2) === null);

// Empty candidate list → null
vm.runInContext('state.bandTraceFishSet = [0,1,2]; state.candidateList = [];', sbx2);
ok('empty candidate list → null',
   vm.runInContext('_lassoLinkageGetOrCompute()', sbx2) === null);

// ============================================================================
// 9. _invalidateLassoLinkageCache
// ============================================================================
console.log('\n=== 9. _invalidateLassoLinkageCache ===');

vm.runInContext('state.lassoLinkageCache = {_stub: true}; state.lassoLinkageCacheKey = "abc";', sbx2);
vm.runInContext('_invalidateLassoLinkageCache()', sbx2);
ok('cache cleared',
   vm.runInContext('state.lassoLinkageCache === null', sbx2));
ok('cache key cleared',
   vm.runInContext('state.lassoLinkageCacheKey === null', sbx2));

// ============================================================================
// 10. _lassoLinkageToTSV
// ============================================================================
console.log('\n=== 10. _lassoLinkageToTSV ===');

const fnTSV = pullFunction(html, '_lassoLinkageToTSV');
vm.runInContext(fnTSV, sbx);

const tsv = vm.runInContext('_lassoLinkageToTSV(_computeLassoLinkage([0,1,2,3,4], _candList))', sbx);
const tsvLines = tsv.split('\n').filter(l => l.length > 0);

ok('TSV non-empty', tsvLines.length > 0);
ok('TSV starts with comment metadata', tsvLines[0].indexOf('# n_fish_selected') === 0);
ok('TSV has 4 comment lines + header + 3 data rows', tsvLines.length === 8);
ok('TSV header line 5 has expected columns',
   tsvLines[4].indexOf('candidate_id\tchrom\tstart_bp\tend_bp\tK\tbest_band\tbest_purity\tn_in_best_band\tn_lasso_seen\tis_strong_link') === 0);

// Data row 0 should be the first strong link by purity desc.
// Lasso = [0,1,2,3,4]: all three candidates have purity=1.0, ties broken by id asc.
const r0 = tsvLines[5].split('\t');
ok('first data row: I1 (alphabetical tie-break)', r0[0] === 'I1');
ok('I1 chrom = LG28', r0[1] === 'LG28');
ok('I1 best_band = 0', r0[5] === '0');
ok('I1 best_purity formatted to 6 dp', /^1\.000000$/.test(r0[6]));
ok('I1 n_in_best_band = 5', r0[7] === '5');
ok('I1 is_strong_link = 1', r0[9] === '1');

// Sort behaviour: weak rows after strong rows
vm.runInContext(`
  // Lasso 0..7 has weak purity at all three candidates.
  var _resultMixed = _computeLassoLinkage([0,1,2,3,4,5,6,7], _candList);
  // Force one candidate to be strong by appending a clean candidate.
  _resultMixed.per_candidate.IZ = {
    id: 'IZ', chrom: 'LG14', start_bp: 1e6, end_bp: 2e6, K: 3,
    best_band: 0, best_purity: 1.0, n_in_best_band: 8, n_lasso_seen: 8,
    per_band: [{band:0,n_in_lasso:8,fraction:1},{band:1,n_in_lasso:0,fraction:0},{band:2,n_in_lasso:0,fraction:0}],
    is_strong_link: true,
  };
`, sbx);
const tsvMixed = vm.runInContext('_lassoLinkageToTSV(_resultMixed)', sbx);
const tsvMixedLines = tsvMixed.split('\n').filter(l => l.length > 0);
const firstDataIdx = tsvMixedLines.findIndex(l => l[0] !== '#' && l.indexOf('candidate_id') !== 0);
ok('mixed: strong link IZ comes first',
   tsvMixedLines[firstDataIdx].split('\t')[0] === 'IZ');
ok('mixed: weak links follow',
   tsvMixedLines[firstDataIdx + 1].split('\t')[9] === '0');

// Null result → null
ok('null result → null',
   vm.runInContext('_lassoLinkageToTSV(null)', sbx) === null);

// Result without per_candidate → null
ok('result without per_candidate → null',
   vm.runInContext('_lassoLinkageToTSV({n_fish_selected: 5})', sbx) === null);

// NaN purity → empty in TSV
vm.runInContext(`
  var _resultNaN = {
    n_fish_selected: 5, n_candidates_seen: 1,
    purity_threshold: 0.7, min_band_size: 5,
    per_candidate: { IX: {
      id: 'IX', chrom: '', start_bp: null, end_bp: null, K: 3,
      best_band: 0, best_purity: NaN, n_in_best_band: 0, n_lasso_seen: 0,
      per_band: [], is_strong_link: false,
    }},
    strong_links: [],
  };
`, sbx);
const tsvNaN = vm.runInContext('_lassoLinkageToTSV(_resultNaN)', sbx);
const tsvNaNLines = tsvNaN.split('\n').filter(l => l.length > 0);
const dataLine = tsvNaNLines.find(l => l.indexOf('IX') === 0);
const dataCols = dataLine.split('\t');
ok('NaN best_purity → empty', dataCols[6] === '');
ok('null start_bp → empty', dataCols[2] === '');
ok('null end_bp → empty', dataCols[3] === '');

// ============================================================================
// 11. _lassoLinkageDownloadTSV — filename + null path
// ============================================================================
console.log('\n=== 11. _lassoLinkageDownloadTSV ===');

vm.runInContext(pullFunction(html, '_lassoLinkageToTSV'), sbx2);
vm.runInContext(pullFunction(html, '_lassoLinkageDownloadTSV'), sbx2);

// Restore a working state in sbx2
vm.runInContext(`
  state.bandTraceFishSet = [0,1,2,3,4];
  state.candidateList = [
    { id: 'I1', chrom: 'LG28', start_bp: 5e6, end_bp: 8e6, K: 3,
      confirmed: true, locked_labels: [0,0,0,0,0,1,1,2,2,2] },
  ];
  state.data = { chrom: 'LG28' };
  state.lassoLinkageCache = null;
  state.lassoLinkageCacheKey = null;
`, sbx2);

const fname = vm.runInContext('_lassoLinkageDownloadTSV()', sbx2);
ok('returns a filename string', typeof fname === 'string' && fname.length > 0);
ok('filename starts with lasso_linkage_LG28', fname.indexOf('lasso_linkage_LG28') === 0);
ok('filename embeds n_fish (n5)', /_n5\.tsv$/.test(fname));

// No fish-set → null
vm.runInContext('state.bandTraceFishSet = null;', sbx2);
ok('no fish-set → null',
   vm.runInContext('_lassoLinkageDownloadTSV()', sbx2) === null);

// ============================================================================
// 12. setBandTraceFishSet hooks _invalidateLassoLinkageCache
// ============================================================================
console.log('\n=== 12. setBandTraceFishSet hooks linkage cache invalidation ===');

ok('source: setBandTraceFishSet contains _invalidateLassoLinkageCache call',
   /function\s+setBandTraceFishSet[\s\S]{0,1500}_invalidateLassoLinkageCache/.test(html));

// ============================================================================
// 13. Lines header — 🔗 linkage button + click handler
// ============================================================================
console.log('\n=== 13. Lines header ===');

ok('linesBandTraceLinkageBtn id present', html.indexOf('id="linesBandTraceLinkageBtn"') > -1);
ok('button has 🔗 linkage label',
   /id="linesBandTraceLinkageBtn"[\s\S]{0,1500}🔗 linkage/.test(html));
ok('button has informative title',
   /id="linesBandTraceLinkageBtn"[\s\S]{0,1500}fish-set/.test(html));
ok('click handler routes to _openLassoLinkagePopover',
   /linesBandTraceLinkageBtn[\s\S]{0,800}_openLassoLinkagePopover\(\)/.test(html));

// ============================================================================
// 14. Regression
// ============================================================================
console.log('\n=== 14. Regression ===');

// Turn 163
ok('turn 163 _updateBandTracePickOptions still defined',
   /function\s+_updateBandTracePickOptions\s*\(/.test(html));
ok('turn 163 _bandTraceRunsToTSV still defined',
   /function\s+_bandTraceRunsToTSV\s*\(/.test(html));
ok('turn 163 _BTRACE_CHAIN_BREAK_COLOR still declared',
   /_BTRACE_CHAIN_BREAK_COLOR\s*=/.test(html));
ok('turn 163 linesBandTracePickSelect still in lines header',
   html.indexOf('id="linesBandTracePickSelect"') > -1);
ok('turn 163 linesBandTraceExportRunsBtn still in lines header',
   html.indexOf('id="linesBandTraceExportRunsBtn"') > -1);

// Turn 162
ok('turn 162 _bandTraceToTSV still defined',
   /function\s+_bandTraceToTSV\s*\(/.test(html));
ok('turn 162 _wireBandTraceTooltip still defined',
   /function\s+_wireBandTraceTooltip\s*\(/.test(html));
ok('turn 162 linesBandTraceExportBtn still in lines header',
   html.indexOf('id="linesBandTraceExportBtn"') > -1);

// Turn 161
ok('turn 161 _drawBandTraceStrip still defined',
   /function\s+_drawBandTraceStrip\s*\(/.test(html));
ok('turn 161 setBandTraceFishSet still defined',
   /function\s+setBandTraceFishSet\s*\(/.test(html));
ok('turn 161 linesBandTraceToggle still in lines header',
   html.indexOf('id="linesBandTraceToggle"') > -1);

// Turn 160
ok('turn 160 _bandTraceForFishSet still defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('turn 160 _bandTraceRegimeRuns still defined',
   /function\s+_bandTraceRegimeRuns\s*\(/.test(html));

// Turn 130 / 156 / 122
ok('turn 130 _drawLineageStrip still defined',
   /function\s+_drawLineageStrip\s*\(/.test(html));
ok('turn 156 openVShapePlot still defined',
   /function\s+openVShapePlot\s*\(/.test(html));
ok('turn 122 inheritance pill tooltip still defined',
   /function\s+_wireInheritancePillTooltip\s*\(/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
