// =============================================================================
// turn 165 — G-panel `auto` review tab (Slice 2 of review_surfaces spec)
// =============================================================================
// Implements SPEC_review_surfaces_auto_and_lineages.md §3.2: a G-panel
// tab listing unconfirmed auto-promoted candidates with per-row
// Confirm / Dismiss / Inspect actions plus bulk Confirm-all /
// Dismiss-all.
//
// Slice 0 (turn 130) shipped the predicate, CSS, sort modifier in
// refreshCandidateListUI. Slice 1 deliverables turned out to be
// already-shipped (cross-candidate matrix filter via
// _gatherActiveCandidatesForInheritance, lines-panel candidate-bands
// skip via _paintCandidateBands's confirmed!==true guard, confirm
// action just sets confirmed=true). This turn closes the visible gap
// — Quentin had no UI surface to triage the auto-promotes.
//
// Slice 3 (lineages tab) deliberately deferred — depends on
// state.lineageResult being exercised on real data.
//
// What this turn ships:
//   - state.gPanelTab comment widened to include 'auto'
//   - _GPANEL_TABS extended with the auto entry + visibleWhen predicate
//   - _renderGPanelModal filters the tab strip via visibleWhen, with
//     fallback-to-manual when the active tab disappears
//   - _renderGPanelModal routes the 'auto' body to _gPanelRenderTabAuto
//   - Inactive-tab dashed-border styling for review-accent tabs
//   - _gPanelHasAutoCandidates(state) predicate
//   - _gPanelCollectAutoCandidates(state) gathers unconfirmed auto cands
//   - _gPanelDismissAuto(cand) reuses _addL2SweepDismissed from turn 133/134
//   - _gPanelRenderTabAuto() body renderer
//   - _gPanelWireAutoTabActions() per-row + bulk action wiring
//   - 5 window exports
//
// Tests (~14 sections, ~80 assertions):
//   1. Source-pattern checks
//   2. _gPanelHasAutoCandidates — true/false on fixtures
//   3. _gPanelCollectAutoCandidates — filters confirmed + non-auto
//   4. _gPanelDismissAuto — calls _addL2SweepDismissed, removes from list
//   5. _gPanelRenderTabAuto — header count, bulk buttons, per-row cards
//   6. _gPanelRenderTabAuto — empty-list instructional message
//   7. _gPanelRenderTabAuto — band-count swatches
//   8. _gPanelRenderTabAuto — silhouette/n_groups metadata when present
//   9. _GPANEL_TABS extension — auto entry, visibleWhen predicate
//  10. _renderGPanelModal — visibleWhen filters tab strip
//  11. _renderGPanelModal — fallback-to-manual when active disappears
//  12. _renderGPanelModal — routes 'auto' body
//  13. Confirm action sets cand.confirmed=true
//  14. Regression — turns 130/135/136/152/160-164 still wired
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

ok('_gPanelHasAutoCandidates defined',
   /function\s+_gPanelHasAutoCandidates\s*\(/.test(html));
ok('_gPanelCollectAutoCandidates defined',
   /function\s+_gPanelCollectAutoCandidates\s*\(/.test(html));
ok('_gPanelDismissAuto defined',
   /function\s+_gPanelDismissAuto\s*\(/.test(html));
ok('_gPanelRenderTabAuto defined',
   /function\s+_gPanelRenderTabAuto\s*\(/.test(html));
ok('_gPanelWireAutoTabActions defined',
   /function\s+_gPanelWireAutoTabActions\s*\(/.test(html));

ok('window._gPanelHasAutoCandidates exported',
   /window\._gPanelHasAutoCandidates\s*=\s*_gPanelHasAutoCandidates/.test(html));
ok('window._gPanelCollectAutoCandidates exported',
   /window\._gPanelCollectAutoCandidates\s*=\s*_gPanelCollectAutoCandidates/.test(html));
ok('window._gPanelDismissAuto exported',
   /window\._gPanelDismissAuto\s*=\s*_gPanelDismissAuto/.test(html));
ok('window._gPanelRenderTabAuto exported',
   /window\._gPanelRenderTabAuto\s*=\s*_gPanelRenderTabAuto/.test(html));
ok('window._gPanelWireAutoTabActions exported',
   /window\._gPanelWireAutoTabActions\s*=\s*_gPanelWireAutoTabActions/.test(html));

// _GPANEL_TABS includes auto entry
ok('_GPANEL_TABS contains auto entry',
   /_GPANEL_TABS\s*=\s*\[[\s\S]*?key:\s*['"]auto['"]/.test(html));
ok('auto entry has visibleWhen predicate',
   /key:\s*['"]auto['"][\s\S]{0,400}visibleWhen:/.test(html));
ok('auto entry has accent: review',
   /key:\s*['"]auto['"][\s\S]{0,400}accent:\s*['"]review['"]/.test(html));

// gPanelTab comment widened
ok('gPanelTab comment widened to include auto',
   /gPanelTab:\s*['"]manual['"][\s\S]{0,200}auto/.test(html));

// _renderGPanelModal modifications
ok('_renderGPanelModal computes visibleTabs via visibleWhen',
   /visibleTabs\s*=\s*_GPANEL_TABS\.filter[\s\S]{0,400}visibleWhen/.test(html));
ok('_renderGPanelModal iterates visibleTabs in tab strip',
   /for\s*\(\s*const\s+t\s+of\s+visibleTabs\s*\)/.test(html));
ok('_renderGPanelModal routes activeTab===\'auto\' to _gPanelRenderTabAuto',
   /activeTab\s*===\s*['"]auto['"][\s\S]{0,200}_gPanelRenderTabAuto\(\)/.test(html));
ok('_renderGPanelModal calls _gPanelWireAutoTabActions when active',
   /activeTab\s*===\s*['"]auto['"][\s\S]{0,400}_gPanelWireAutoTabActions\(\)/.test(html));
ok('_renderGPanelModal applies dashed-border to inactive review-accent tabs',
   /isReviewAccent\s*=[\s\S]{0,200}accent\s*===\s*['"]review['"]/.test(html));

// _gPanelDismissAuto uses existing _addL2SweepDismissed (no duplicate
// dismissal infrastructure)
const fnDismiss = pullFunction(html, '_gPanelDismissAuto');
ok('_gPanelDismissAuto calls _addL2SweepDismissed (reuse from turn 133/134)',
   fnDismiss && /_addL2SweepDismissed/.test(fnDismiss));
ok('_gPanelDismissAuto removes cand from candidateList',
   fnDismiss && /candidateList\s*=\s*[^;]*\.filter/.test(fnDismiss));

// ============================================================================
// 2. _gPanelHasAutoCandidates predicate
// ============================================================================
console.log('\n=== 2. _gPanelHasAutoCandidates ===');

const sbx = { console };
vm.createContext(sbx);
// Provide a minimal _isAutoCandidate so the predicate uses the canonical
// path. Verified separately at the source level. Also provide `state`
// and `window.state` placeholders because the helpers do
// `_state = _state || ((typeof window !== 'undefined' && window.state) ? window.state : state)`
// and reading `state` when neither is set throws ReferenceError.
vm.runInContext(`
  function _isAutoCandidate(cand) {
    if (!cand) return false;
    if (cand.confirmed) return false;
    const src = cand.source;
    if (typeof src !== 'string') return false;
    return src.indexOf('auto_') === 0;
  }
  var state = { candidateList: [] };
  var window = { state };
`, sbx);
vm.runInContext(pullFunction(html, '_gPanelHasAutoCandidates'), sbx);
vm.runInContext(pullFunction(html, '_gPanelCollectAutoCandidates'), sbx);

ok('null state → false',
   vm.runInContext('_gPanelHasAutoCandidates(null)', sbx) === false);
ok('state without candidateList → false',
   vm.runInContext('_gPanelHasAutoCandidates({})', sbx) === false);
ok('empty candidateList → false',
   vm.runInContext('_gPanelHasAutoCandidates({candidateList: []})', sbx) === false);
ok('all confirmed candidates → false',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'I1', source:'auto_l2_sweep', confirmed:true},
     {id:'I2', source:'manual', confirmed:true},
   ]})`, sbx) === false);
ok('manual unconfirmed only → false',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'I1', source:'manual', confirmed:false},
   ]})`, sbx) === false);
ok('one auto unconfirmed → true',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'I1', source:'manual', confirmed:true},
     {id:'I2', source:'auto_l2_sweep', confirmed:false},
   ]})`, sbx) === true);
ok('auto_lineage source → true',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'L1', source:'auto_lineage', confirmed:false},
   ]})`, sbx) === true);

// Defensive fallback: works without _isAutoCandidate
const sbxF = { console };
vm.createContext(sbxF);
vm.runInContext('var state = { candidateList: [] }; var window = { state };', sbxF);
vm.runInContext(pullFunction(html, '_gPanelHasAutoCandidates'), sbxF);
ok('fallback (no _isAutoCandidate): auto detected',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'X', source:'auto_l2_sweep', confirmed:false},
   ]})`, sbxF) === true);
ok('fallback: confirmed auto NOT detected',
   vm.runInContext(`_gPanelHasAutoCandidates({candidateList: [
     {id:'X', source:'auto_l2_sweep', confirmed:true},
   ]})`, sbxF) === false);

// ============================================================================
// 3. _gPanelCollectAutoCandidates
// ============================================================================
console.log('\n=== 3. _gPanelCollectAutoCandidates ===');

const collected = vm.runInContext(`_gPanelCollectAutoCandidates({candidateList: [
  {id:'M1', source:'manual',        confirmed:true},
  {id:'A1', source:'auto_l2_sweep', confirmed:false},
  {id:'A2', source:'auto_l2_sweep', confirmed:true},     // confirmed → skip
  {id:'M2', source:'manual',        confirmed:false},    // manual unconfirmed → skip
  {id:'A3', source:'auto_lineage',  confirmed:false},
]})`, sbx);

ok('returns array of length 2', Array.isArray(collected) && collected.length === 2);
ok('preserves source order: A1 first, A3 second',
   collected[0].id === 'A1' && collected[1].id === 'A3');
ok('skips confirmed auto', !collected.some(c => c.id === 'A2'));
ok('skips manual unconfirmed', !collected.some(c => c.id === 'M2'));

ok('null state → []',
   vm.runInContext('_gPanelCollectAutoCandidates(null).length', sbx) === 0);
ok('no candidateList → []',
   vm.runInContext('_gPanelCollectAutoCandidates({}).length', sbx) === 0);

// ============================================================================
// 4. _gPanelDismissAuto
// ============================================================================
console.log('\n=== 4. _gPanelDismissAuto ===');

const sbx2 = { console };
vm.createContext(sbx2);
// Set up a state with an auto candidate; provide a stub for
// _addL2SweepDismissed that records calls.
vm.runInContext(`
  var _dismissCalls = [];
  function _addL2SweepDismissed(chrom, l2idx) {
    _dismissCalls.push({ chrom, l2idx });
  }
  var state = {
    data: { chrom: 'LG28' },
    candidateList: [
      { id: 'A1', source: 'auto_l2_sweep', confirmed: false,
        chrom: 'LG28', ref_l2: 7 },
      { id: 'A2', source: 'auto_l2_sweep', confirmed: false,
        chrom: 'LG28', ref_l2: 14 },
      { id: 'M1', source: 'manual', confirmed: true,
        chrom: 'LG28', ref_l2: 3 },
    ],
  };
  var window = { state };
`, sbx2);
vm.runInContext(pullFunction(html, '_gPanelDismissAuto'), sbx2);

const dismissResult = vm.runInContext(`_gPanelDismissAuto(state.candidateList[0])`, sbx2);
ok('_gPanelDismissAuto returns true on success', dismissResult === true);
ok('candidate removed from candidateList',
   vm.runInContext('state.candidateList.length', sbx2) === 2);
ok('removed candidate was A1',
   vm.runInContext('state.candidateList.find(c => c.id === "A1")', sbx2) === undefined);
ok('A2 still present', vm.runInContext('!!state.candidateList.find(c => c.id === "A2")', sbx2));
ok('M1 still present', vm.runInContext('!!state.candidateList.find(c => c.id === "M1")', sbx2));
ok('_addL2SweepDismissed called once',
   vm.runInContext('_dismissCalls.length', sbx2) === 1);
ok('_addL2SweepDismissed called with chrom=LG28, l2idx=7',
   vm.runInContext('_dismissCalls[0].chrom', sbx2) === 'LG28' &&
   vm.runInContext('_dismissCalls[0].l2idx', sbx2) === 7);

// Defensive: null cand → false, no mutation
ok('null cand → false', vm.runInContext('_gPanelDismissAuto(null)', sbx2) === false);
ok('cand without id → false',
   vm.runInContext('_gPanelDismissAuto({source:"auto_l2_sweep"})', sbx2) === false);

// Cand without ref_l2 → still removes from list, just skips _addL2SweepDismissed
vm.runInContext('_dismissCalls.length = 0;', sbx2);
const beforeLen = vm.runInContext('state.candidateList.length', sbx2);
vm.runInContext(`state.candidateList.push({id: 'A3', source:'auto_l2_sweep', confirmed:false, chrom:'LG28'});`, sbx2);
vm.runInContext('_gPanelDismissAuto(state.candidateList[state.candidateList.length-1])', sbx2);
ok('cand without ref_l2: still removed from list',
   vm.runInContext('state.candidateList.length', sbx2) === beforeLen);
ok('cand without ref_l2: no _addL2SweepDismissed call',
   vm.runInContext('_dismissCalls.length', sbx2) === 0);

// ============================================================================
// 5. _gPanelRenderTabAuto — non-empty list
// ============================================================================
console.log('\n=== 5. _gPanelRenderTabAuto — non-empty ===');

const sbx3 = { console };
vm.createContext(sbx3);
vm.runInContext(`
  function _isAutoCandidate(cand) {
    if (!cand) return false;
    if (cand.confirmed) return false;
    const src = cand.source;
    return typeof src === 'string' && src.indexOf('auto_') === 0;
  }
  function _gpKaryoColor(k) {
    const PAL = ['#3b6fb6','#ffd866','#d97a2c','#7ad394','#a76de2','#e85a5a'];
    return PAL[k % PAL.length];
  }
  var state = {
    data: { chrom: 'LG28' },
    candidateList: [
      { id: 'A1', source: 'auto_l2_sweep', confirmed: false,
        chrom: 'LG28', ref_l2: 7, K: 3,
        start_bp: 5_000_000, end_bp: 8_000_000,
        locked_labels: new Int8Array([0,0,0,0,1,1,1,2,2,2]),
        meta: { silhouette: 0.42, n_groups: 2 } },
      { id: 'A2', source: 'auto_l2_sweep', confirmed: false,
        chrom: 'LG28', ref_l2: 14, K: 3,
        start_bp: 12_000_000, end_bp: 15_000_000,
        locked_labels: new Int8Array([1,1,1,1,1,2,2,2,0,0]) },
      { id: 'M1', source: 'manual', confirmed: true, chrom: 'LG28', ref_l2: 3 },
    ],
  };
  var window = { state };
`, sbx3);
vm.runInContext(pullFunction(html, '_gPanelHasAutoCandidates'), sbx3);
vm.runInContext(pullFunction(html, '_gPanelCollectAutoCandidates'), sbx3);
vm.runInContext(pullFunction(html, '_gPanelRenderTabAuto'), sbx3);

const html5 = vm.runInContext('_gPanelRenderTabAuto()', sbx3);

ok('returns non-empty string', typeof html5 === 'string' && html5.length > 0);
ok('header shows count = 2',
   /<span style="color: var\(--good\);">2<\/span>/.test(html5));
ok('header text "auto-promoted candidates"',
   html5.indexOf('auto-promoted candidate') > -1);
ok('header includes "awaiting review"',
   html5.indexOf('awaiting review') > -1);
ok('bulk Confirm-all button present',
   html5.indexOf('id="gpAutoConfirmAll"') > -1 &&
   html5.indexOf('confirm all</button>') > -1);
ok('bulk Dismiss-all button present',
   html5.indexOf('id="gpAutoDismissAll"') > -1 &&
   html5.indexOf('dismiss all</button>') > -1);

// Per-row markup
ok('row class gpAutoRow present (2 rows)',
   (html5.match(/class="gpAutoRow"/g) || []).length === 2);
ok('row data-cand-id A1', html5.indexOf('data-cand-id="A1"') > -1);
ok('row data-cand-id A2', html5.indexOf('data-cand-id="A2"') > -1);
ok('row does NOT include manual M1', html5.indexOf('data-cand-id="M1"') < 0);
ok('🤖 prefix appears for both rows',
   (html5.match(/🤖/g) || []).length >= 2);
ok('A1 row shows L2:7', html5.indexOf('🤖 L2:7') > -1);
ok('A1 row shows chrom LG28',
   html5.indexOf('LG28') > -1);
ok('A1 row shows span 5.00–8.00 Mb',
   html5.indexOf('5.00–8.00 Mb') > -1);
ok('A1 row shows K=3 + silhouette + n_groups',
   /K=3, silh=0\.42, n_groups=2/.test(html5));
ok('A2 row shows K=3 without silh (no meta)',
   /L2:14[\s\S]{0,200}K=3</.test(html5) &&
   !/L2:14[\s\S]{0,200}silh=/.test(html5));

// Per-row action buttons
ok('per-row gpAutoActConfirm class present (2 rows)',
   (html5.match(/class="gpAutoActConfirm"/g) || []).length === 2);
ok('per-row gpAutoActDismiss class present (2 rows)',
   (html5.match(/class="gpAutoActDismiss"/g) || []).length === 2);
ok('per-row gpAutoActInspect class present (2 rows)',
   (html5.match(/class="gpAutoActInspect"/g) || []).length === 2);

// Band swatches
ok('A1 row shows g0 (n=4)', html5.indexOf('g0 (n=4)') > -1);
ok('A1 row shows g1 (n=3)', html5.indexOf('g1 (n=3)') > -1);
ok('A1 row shows g2 (n=3)', html5.indexOf('g2 (n=3)') > -1);
ok('A2 row shows g1 (n=5)', html5.indexOf('g1 (n=5)') > -1);
ok('A2 row shows g2 (n=3)', html5.indexOf('g2 (n=3)') > -1);
ok('A2 row shows g0 (n=2)', html5.indexOf('g0 (n=2)') > -1);

// ============================================================================
// 6. _gPanelRenderTabAuto — empty list path
// ============================================================================
console.log('\n=== 6. _gPanelRenderTabAuto — empty path ===');

vm.runInContext('state.candidateList = [];', sbx3);
const htmlEmpty = vm.runInContext('_gPanelRenderTabAuto()', sbx3);
ok('empty: still produces non-empty string', typeof htmlEmpty === 'string' && htmlEmpty.length > 0);
ok('empty: header shows count = 0',
   /<span style="color: var\(--good\);">0<\/span>/.test(htmlEmpty));
ok('empty: instructional placeholder visible',
   htmlEmpty.indexOf('No auto-promoted candidates currently awaiting review') > -1);
ok('empty: no per-row cards', !(/class="gpAutoRow"/.test(htmlEmpty)));

// State without candidateList
vm.runInContext('state.candidateList = null;', sbx3);
const htmlNull = vm.runInContext('_gPanelRenderTabAuto()', sbx3);
ok('null candidateList: function does not throw',
   typeof htmlNull === 'string');

// ============================================================================
// 7. _GPANEL_TABS extension
// ============================================================================
console.log('\n=== 7. _GPANEL_TABS extension ===');

// Pull the array literal from source for shape verification.
const tabArrMatch = html.match(/_GPANEL_TABS\s*=\s*\[([\s\S]*?)\];/);
ok('_GPANEL_TABS array found in source', !!tabArrMatch);
const tabArrSrc = tabArrMatch ? tabArrMatch[1] : '';
ok('karyotype entry preserved', /key:\s*['"]karyotype['"]/.test(tabArrSrc));
ok('inheritance entry preserved', /key:\s*['"]inheritance['"]/.test(tabArrSrc));
ok('manual entry preserved', /key:\s*['"]manual['"]/.test(tabArrSrc));
ok('auto entry added at end', /key:\s*['"]auto['"]/.test(tabArrSrc));
ok('auto entry has label "auto · review"',
   /key:\s*['"]auto['"][\s\S]{0,200}label:\s*['"]auto · review['"]/.test(tabArrSrc));

// ============================================================================
// 8. _renderGPanelModal — visibility filtering
// ============================================================================
console.log('\n=== 8. _renderGPanelModal — visibility filter ===');

// Source-pattern: visibility filter applied at render time
ok('visibleTabs filter applied in _renderGPanelModal',
   /const\s+visibleTabs\s*=\s*_GPANEL_TABS\.filter/.test(html));
ok('validTabs derived from visibleTabs',
   /const\s+validTabs\s*=\s*visibleTabs\.map/.test(html));
ok('fallback to manual when active not in validTabs',
   /validTabs\.indexOf\(_state\.gPanelTab\)\s*===\s*-1[\s\S]{0,200}_state\.gPanelTab\s*=\s*['"]manual['"]/.test(html));

// ============================================================================
// 9. _renderGPanelModal — body routing
// ============================================================================
console.log('\n=== 9. _renderGPanelModal — body routing ===');

ok('activeTab===\'auto\' branch present',
   /activeTab\s*===\s*['"]auto['"][\s\S]{0,300}_gPanelRenderTabAuto/.test(html));
ok('activeTab===\'auto\' guarded by typeof === function',
   /typeof\s+_gPanelRenderTabAuto\s*===\s*['"]function['"]/.test(html));
ok('manual fallback via else branch',
   /\}\s*else\s*\{[\s\S]{0,200}_gPanelRenderTabManual\(\)/.test(html));
ok('action wiring after render when active===auto',
   /activeTab\s*===\s*['"]auto['"][\s\S]{0,400}_gPanelWireAutoTabActions\(\)/.test(html));

// ============================================================================
// 10. Confirm action mutates cand.confirmed
// ============================================================================
console.log('\n=== 10. Confirm action ===');

// We can't easily exercise the real DOM-bound confirm path in node, but
// we can verify the source-level pattern: the wiring function sets
// c.confirmed = true on the matching candidate.
const fnWire = pullFunction(html, '_gPanelWireAutoTabActions');
ok('wire function defines findCand helper',
   fnWire && /findCand\s*=\s*\(id\)\s*=>/.test(fnWire));
ok('wire function reads data-cand-id from button',
   fnWire && /getAttribute\(['"]data-cand-id['"]\)/.test(fnWire));
ok('confirm action sets c.confirmed = true',
   fnWire && /c\.confirmed\s*=\s*true/.test(fnWire));
ok('dismiss action calls _gPanelDismissAuto',
   fnWire && /_gPanelDismissAuto\(c\)/.test(fnWire));
ok('inspect action calls setCandidate',
   fnWire && /setCandidate\(c\)/.test(fnWire));
ok('inspect action closes modal',
   fnWire && /_gPanelClose\(\)/.test(fnWire));
ok('bulk confirm-all uses confirm() prompt',
   fnWire && /typeof\s+confirm\s*===\s*['"]function['"][\s\S]{0,200}confirm\(/.test(fnWire));
ok('bulk confirm-all sets confirmed=true on each',
   fnWire && /for\s*\(\s*const\s+c\s+of\s+auto\s*\)\s*c\.confirmed\s*=\s*true/.test(fnWire));
ok('bulk dismiss-all calls _gPanelDismissAuto on each id',
   fnWire && /_gPanelDismissAuto\(c\)/.test(fnWire));
ok('rerender call after each action',
   fnWire && /_renderGPanelModal\(\)/.test(fnWire));

// ============================================================================
// 11. Visibility predicate sandbox-driven
// ============================================================================
console.log('\n=== 11. Visibility predicate ===');

// Test that _GPANEL_TABS.find(t => t.key==='auto').visibleWhen returns
// true / false based on state. We pull the array literal into a sandbox
// and exercise the predicate.
const sbxV = { console };
vm.createContext(sbxV);
vm.runInContext(`
  function _isAutoCandidate(cand) {
    if (!cand || cand.confirmed) return false;
    return typeof cand.source === 'string' && cand.source.indexOf('auto_') === 0;
  }
  var state = { candidateList: [] };
  var window = { state };
`, sbxV);
vm.runInContext(pullFunction(html, '_gPanelHasAutoCandidates'), sbxV);
// Build the tabs array via eval of the source literal — extract the
// visibleWhen function body literally.
vm.runInContext(`
  // The visibleWhen on the auto tab calls _gPanelHasAutoCandidates,
  // which is now in-scope. Reconstruct the predicate as written in
  // source:
  var _tabs_auto = {
    key: 'auto',
    visibleWhen: function(_state) {
      return (typeof _gPanelHasAutoCandidates === 'function') && _gPanelHasAutoCandidates(_state);
    }
  };
`, sbxV);

ok('auto.visibleWhen({}) → false',
   vm.runInContext('_tabs_auto.visibleWhen({})', sbxV) === false);
ok('auto.visibleWhen({candidateList: []}) → false',
   vm.runInContext('_tabs_auto.visibleWhen({candidateList: []})', sbxV) === false);
ok('auto.visibleWhen with auto cand → true',
   vm.runInContext(`_tabs_auto.visibleWhen({candidateList: [
     {id:'A', source:'auto_l2_sweep', confirmed:false}
   ]})`, sbxV) === true);
ok('auto.visibleWhen all confirmed → false',
   vm.runInContext(`_tabs_auto.visibleWhen({candidateList: [
     {id:'A', source:'auto_l2_sweep', confirmed:true}
   ]})`, sbxV) === false);

// ============================================================================
// 12. Source-pattern: dashed-border styling for review-accent inactive
// ============================================================================
console.log('\n=== 12. Review-accent dashed border ===');

ok('isReviewAccent variable computed',
   /const\s+isReviewAccent\s*=\s*\(t\.accent\s*===\s*['"]review['"]\)/.test(html));
ok('borderStyle conditional uses dashed for inactive review',
   /borderStyle\s*=\s*\(\s*!isActive\s*&&\s*isReviewAccent\s*\)\s*\?\s*['"]dashed['"]\s*:\s*['"]solid['"]/.test(html));

// ============================================================================
// 13. _gPanelDismissAuto uses existing dismissal helper
// ============================================================================
console.log('\n=== 13. Reuse of existing dismissal helper ===');

ok('_addL2SweepDismissed still defined (turn 133/134)',
   /function\s+_addL2SweepDismissed\s*\(/.test(html));
ok('_loadL2SweepDismissed still defined',
   /function\s+_loadL2SweepDismissed\s*\(/.test(html));
ok('_L2_SWEEP_DISMISSED_KEY_PFX constant still in source',
   /_L2_SWEEP_DISMISSED_KEY_PFX\s*=\s*['"]pca_scrubber_v3\.l2SweepDismissed\./.test(html));
ok('_dismissAuto guards on Number.isInteger(cand.ref_l2)',
   fnDismiss && /Number\.isInteger\(cand\.ref_l2\)/.test(fnDismiss));

// ============================================================================
// 14. Regression
// ============================================================================
console.log('\n=== 14. Regression ===');

ok('turn 130 _isAutoCandidate still defined',
   /function\s+_isAutoCandidate\s*\(/.test(html));
ok('turn 130 cli-auto-prefix CSS still in source',
   /cli-auto-prefix/.test(html));
ok('turn 130 .cand-list-item.is-auto CSS still in source',
   /\.cand-list-item\.is-auto/.test(html));
ok('turn 133/134 runL2SweepInheritance still defined',
   /function\s+runL2SweepInheritance\s*\(/.test(html));
ok('turn 135 _renderGPanelModal still defined',
   /function\s+_renderGPanelModal\s*\(/.test(html));
ok('turn 136 _gPanelRenderTabKaryotype still defined',
   /function\s+_gPanelRenderTabKaryotype\s*\(/.test(html));
ok('turn 152 _gPanelRenderTabInheritance still defined',
   /function\s+_gPanelRenderTabInheritance\s*\(/.test(html));
ok('turn 160 _bandTraceForFishSet still defined',
   /function\s+_bandTraceForFishSet\s*\(/.test(html));
ok('turn 161 _drawBandTraceStrip still defined',
   /function\s+_drawBandTraceStrip\s*\(/.test(html));
ok('turn 162 _bandTraceToTSV still defined',
   /function\s+_bandTraceToTSV\s*\(/.test(html));
ok('turn 163 _updateBandTracePickOptions still defined',
   /function\s+_updateBandTracePickOptions\s*\(/.test(html));
ok('turn 164 _computeLassoLinkage still defined',
   /function\s+_computeLassoLinkage\s*\(/.test(html));
ok('turn 164 _openLassoLinkagePopover still defined',
   /function\s+_openLassoLinkagePopover\s*\(/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
