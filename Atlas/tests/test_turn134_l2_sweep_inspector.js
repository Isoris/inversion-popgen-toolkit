// =============================================================================
// turn 134 — Slice 2: L2-sweep inspector modal (SPEC_l2_sweep_inheritance.md §4.2)
// =============================================================================
// Builds on turn 133 Slice 1. Adds:
//   - _renderL2SweepInspectorModal() — DOM renderer for the modal
//   - _l2siBuildRows() — derives displayable rows from state.l2SweepResult
//   - _l2siRowMatchesFilter(meta, status, filter) — filter predicate
//   - _l2siGroupsTouching(result, item_idx) — count groups touching a row
//   - _l2siStatusFor(l2idx, dismissed, candidateList) — promoted/dismissed/pending
//   - _l2siPromoteRow(l2idx) — explicit per-row promote (gate bypass)
//   - _l2siDismissRow(l2idx) — explicit per-row dismiss + drop auto candidate
//   - _l2siFmtBpRange(s, e) — bp range to "1.00–2.00 Mb (1.00 Mb)" string
//   - _l2siFilterState() — lazy filter state on state.l2SweepInspectorFilter
//   - _refreshL2SweepInspectBtnAvailability() — disable trigger when no sweep
//
// Trigger button: #l2SweepInspectBtn in L3 toolbar, sibling to L2-sweep
// checkbox. Modal overlay: #l2SweepInspectorOverlay, sibling to schema +
// JS-scripts modals. Same close pattern (✕ button, click outside, Esc).
//
// applyData hook ALSO added in Slice 1 follow-up — chrom switch with
// toggle on re-runs the sweep. Tests verify the hook block presence.
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
// 1. Source-level: function definitions + DOM elements + chrom-load hook
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_renderL2SweepInspectorModal defined',
   /function _renderL2SweepInspectorModal\b/.test(html));
ok('_l2siBuildRows defined',         /function _l2siBuildRows\b/.test(html));
ok('_l2siRowMatchesFilter defined',  /function _l2siRowMatchesFilter\b/.test(html));
ok('_l2siGroupsTouching defined',    /function _l2siGroupsTouching\b/.test(html));
ok('_l2siStatusFor defined',         /function _l2siStatusFor\b/.test(html));
ok('_l2siPromoteRow defined',        /function _l2siPromoteRow\b/.test(html));
ok('_l2siDismissRow defined',        /function _l2siDismissRow\b/.test(html));
ok('_l2siFmtBpRange defined',        /function _l2siFmtBpRange\b/.test(html));
ok('_l2siFilterState defined',       /function _l2siFilterState\b/.test(html));
ok('_l2siStatusChip defined',        /function _l2siStatusChip\b/.test(html));
ok('_refreshL2SweepInspectBtnAvailability defined',
   /function _refreshL2SweepInspectBtnAvailability\b/.test(html));

// Window exports
ok('window._renderL2SweepInspectorModal exported',
   /window\._renderL2SweepInspectorModal\b/.test(html));
ok('window._l2siBuildRows exported',
   /window\._l2siBuildRows\b/.test(html));
ok('window._l2siPromoteRow / _l2siDismissRow exported',
   /window\._l2siPromoteRow\b/.test(html) && /window\._l2siDismissRow\b/.test(html));

// Filter defaults
ok('_L2SI_DEFAULT_MIN_SIL = 0', /const _L2SI_DEFAULT_MIN_SIL\s*=\s*0\b/.test(html));
ok('_L2SI_DEFAULT_MIN_GRPS = 0', /const _L2SI_DEFAULT_MIN_GRPS\s*=\s*0\b/.test(html));
ok("_L2SI_DEFAULT_STATUS = 'all'",
   /const _L2SI_DEFAULT_STATUS\s*=\s*'all'/.test(html));

// ============================================================================
// 2. DOM: trigger button + overlay div
// ============================================================================
console.log('\n=== 2. DOM elements ===');

ok('#l2SweepInspectorOverlay div present',
   /id="l2SweepInspectorOverlay"/.test(html));
ok('#l2SweepInspectorOverlay sibling to other overlays (z-index 1000)',
   /id="l2SweepInspectorOverlay"[\s\S]{0,400}z-index:\s*1000/.test(html));
ok('#l2SweepInspectorOverlay starts display:none',
   /id="l2SweepInspectorOverlay"[\s\S]{0,200}display:\s*none/.test(html));
ok('#l2SweepInspectBtn trigger button present',
   /id="l2SweepInspectBtn"/.test(html));
ok('#l2SweepInspectBtn starts disabled',
   /id="l2SweepInspectBtn"[\s\S]{0,800}disabled\s*>/.test(html));
ok('inspector button shows 🔍 inspect',
   /id="l2SweepInspectBtn"[\s\S]{0,600}🔍\s*inspect/.test(html));
ok('inspector button sibling to L2-sweep checkbox',
   /l2SweepToggleLabel[\s\S]{0,2500}l2SweepInspectBtn/.test(html));

// ============================================================================
// 3. applyData chrom-load hook for L2-sweep
// ============================================================================
console.log('\n=== 3. applyData chrom-load hook ===');

ok('applyData calls invalidateL2SweepCache on chrom change',
   /applyData[\s\S]{0,8000}invalidateL2SweepCache/.test(html));
ok('applyData re-runs sweep when state.l2SweepEnabled is true',
   /state\.l2SweepEnabled[\s\S]{0,500}runL2SweepInheritance\(\s*\{\s*force:\s*true\s*\}\s*\)/.test(html));
ok('applyData calls _autoPromoteFromSweep after sweep',
   /sweepRes\s*=\s*runL2SweepInheritance[\s\S]{0,200}_autoPromoteFromSweep\(sweepRes\)/.test(html));
ok('applyData refreshes inspect-btn availability after sweep',
   /\[l2sweep\]\s*applyData hook failed[\s\S]{0,800}_refreshL2SweepInspectBtnAvailability/.test(html));

// ============================================================================
// 4. Behavioral — filter predicate (pure function, easy to sandbox)
// ============================================================================
console.log('\n=== 4. _l2siRowMatchesFilter behavior ===');

function extractFnSrc(name) {
  const re = new RegExp('function ' + name + '\\([\\s\\S]*?\\n\\}', 'm');
  const m = html.match(re);
  return m ? m[0] : null;
}

const filterFnSrc = extractFnSrc('_l2siRowMatchesFilter');
ok('extracted _l2siRowMatchesFilter source', filterFnSrc !== null);

const filterCtx = vm.createContext({
  Number, isFinite, Array,
});
vm.runInContext(filterFnSrc, filterCtx);

// 4a. No filter → always passes
ok('null filter passes all',
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.1, n_per_group: [1,1,1] }, 'pending', null) === true);

// 4b. minSil filters out low-sil
ok('minSil 0.5 rejects sil 0.3',
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.3, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0.5, minGroups: 0, status: 'all' }) === false);
ok('minSil 0.5 accepts sil 0.7',
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0.5, minGroups: 0, status: 'all' }) === true);
ok('minSil > 0 rejects NaN silhouette',
   filterCtx._l2siRowMatchesFilter({ silhouette: NaN, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0.1, minGroups: 0, status: 'all' }) === false);
ok('minSil = 0 accepts NaN silhouette',
   filterCtx._l2siRowMatchesFilter({ silhouette: NaN, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0, minGroups: 0, status: 'all' }) === true);

// 4c. minGroups filters by populated bands
ok('minGroups 3 rejects [5,5,0] (only 2 populated)',
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,0] }, 'pending',
                                    { minSil: 0, minGroups: 3, status: 'all' }) === false);
ok('minGroups 3 accepts [5,5,5]',
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0, minGroups: 3, status: 'all' }) === true);

// 4d. status filter
ok("status='promoted' rejects pending row",
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,5] }, 'pending',
                                    { minSil: 0, minGroups: 0, status: 'promoted' }) === false);
ok("status='promoted' accepts promoted row",
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,5] }, 'promoted',
                                    { minSil: 0, minGroups: 0, status: 'promoted' }) === true);
ok("status='all' accepts every status",
   filterCtx._l2siRowMatchesFilter({ silhouette: 0.7, n_per_group: [5,5,5] }, 'dismissed',
                                    { minSil: 0, minGroups: 0, status: 'all' }) === true);

// ============================================================================
// 5. Behavioral — _l2siGroupsTouching (count groups per item_idx)
// ============================================================================
console.log('\n=== 5. _l2siGroupsTouching behavior ===');

const groupsFnSrc = extractFnSrc('_l2siGroupsTouching');
const groupsCtx = vm.createContext({ Number, Array, Set });
vm.runInContext(groupsFnSrc, groupsCtx);

const fakeResult = {
  groups: [
    { members: [{ item_idx: 0, band: 0 }, { item_idx: 1, band: 0 }] },
    { members: [{ item_idx: 0, band: 1 }, { item_idx: 2, band: 0 }] },
    { members: [{ item_idx: 1, band: 1 }, { item_idx: 2, band: 1 }] },
  ],
};
ok('item 0 touches 2 groups (#0 and #1)',
   groupsCtx._l2siGroupsTouching(fakeResult, 0) === 2);
ok('item 1 touches 2 groups (#0 and #2)',
   groupsCtx._l2siGroupsTouching(fakeResult, 1) === 2);
ok('item 2 touches 2 groups (#1 and #2)',
   groupsCtx._l2siGroupsTouching(fakeResult, 2) === 2);
ok('item 99 touches 0 groups',
   groupsCtx._l2siGroupsTouching(fakeResult, 99) === 0);

// Alternate shape (cells instead of members) — should also work
const altResult = {
  groups: [
    { cells: [{ itemIdx: 5, band: 0 }] },
    { cells: [{ itemIdx: 5, band: 1 }] },
  ],
};
ok('alternate "cells" shape with itemIdx works',
   groupsCtx._l2siGroupsTouching(altResult, 5) === 2);

// Empty groups → 0
ok('null result returns 0', groupsCtx._l2siGroupsTouching(null, 0) === 0);
ok('result without groups returns 0',
   groupsCtx._l2siGroupsTouching({}, 0) === 0);

// ============================================================================
// 6. Behavioral — _l2siStatusFor (uses dismissed Set + candidateList)
// ============================================================================
console.log('\n=== 6. _l2siStatusFor behavior ===');

const statusFnSrc = extractFnSrc('_l2siStatusFor');
const statusCtx = vm.createContext({ Array });
vm.runInContext(statusFnSrc, statusCtx);

const dismissedSet = new Set([7, 12]);
const candList = [
  { id: 'a', l2_indices: [3, 5] },
  { id: 'b', l2_indices: [9] },
];
ok('dismissed L2 returns "dismissed"',
   statusCtx._l2siStatusFor(7, dismissedSet, candList) === 'dismissed');
ok('L2 in candidate returns "promoted"',
   statusCtx._l2siStatusFor(3, dismissedSet, candList) === 'promoted');
ok('L2 in candidate (different position) returns "promoted"',
   statusCtx._l2siStatusFor(9, dismissedSet, candList) === 'promoted');
ok('L2 not in either set returns "pending"',
   statusCtx._l2siStatusFor(99, dismissedSet, candList) === 'pending');
ok('null dismissed set still works (no crash)',
   statusCtx._l2siStatusFor(99, null, candList) === 'pending');
ok('null candidateList still works (no crash)',
   statusCtx._l2siStatusFor(99, dismissedSet, null) === 'pending');
// Edge: dismissed AND in candidateList → spec says dismissed wins (the
// dismissed set is checked first; that's the user's most recent statement).
ok('dismissed wins over promoted (user just dismissed)',
   statusCtx._l2siStatusFor(7, dismissedSet, [{ id: 'x', l2_indices: [7] }]) === 'dismissed');

// ============================================================================
// 7. Behavioral — _l2siFmtBpRange formatting
// ============================================================================
console.log('\n=== 7. _l2siFmtBpRange behavior ===');

const fmtFnSrc = extractFnSrc('_l2siFmtBpRange');
const fmtCtx = vm.createContext({});
vm.runInContext(fmtFnSrc, fmtCtx);

ok('1Mb-2Mb formats correctly',
   fmtCtx._l2siFmtBpRange(1e6, 2e6) === '1.00–2.00 Mb (1.00 Mb)');
ok('14Mb-16.5Mb formats correctly',
   fmtCtx._l2siFmtBpRange(14e6, 16.5e6) === '14.00–16.50 Mb (2.50 Mb)');
ok('null start returns —',
   fmtCtx._l2siFmtBpRange(null, 1e6) === '—');
ok('non-number end returns —',
   fmtCtx._l2siFmtBpRange(1e6, 'foo') === '—');

// ============================================================================
// 8. Behavioral — full sandbox: build rows + filter + actions
// ============================================================================
console.log('\n=== 8. Full sandbox integration ===');

function extractInspectorBlock() {
  const start = html.indexOf('// L2-SWEEP INSPECTOR MODAL');
  if (start < 0) return null;
  const endMarker = 'window._l2siFilterState                   = _l2siFilterState;';
  const endIdx = html.indexOf(endMarker, start);
  if (endIdx < 0) return null;
  const after = html.indexOf('}', endIdx);
  if (after < 0) return null;
  return html.substring(start, after + 1);
}
const inspectorBlock = extractInspectorBlock();
ok('inspector code block extracts cleanly',
   inspectorBlock !== null && inspectorBlock.length > 1000);

function makeInspectorSandbox(scenario) {
  const _ls = scenario.ls || {};
  const localStorageStub = {
    getItem(k) { return Object.prototype.hasOwnProperty.call(_ls, k) ? _ls[k] : null; },
    setItem(k, v) { _ls[k] = String(v); },
    removeItem(k) { delete _ls[k]; },
  };
  const stateStub = scenario.state || {
    data: { chrom: 'LG28', l2_envelopes: [] },
    candidateList: [],
    candidates: {},
    k: 3,
    activeMode: 'default',
    l2SweepResult: null,
  };
  // Stubs for promote/dismiss path
  const removedIds = [];
  const removeCandidateFullyStub = function (id) {
    stateStub.candidateList = stateStub.candidateList.filter(c => c.id !== id);
    removedIds.push(id);
    return true;
  };
  const addedCands = [];
  const addCandidateToListStub = function (cand) {
    if (!cand || !cand.id) return;
    if (stateStub.candidateList.some(c => c && c.id === cand.id)) return;
    stateStub.candidateList.push(cand);
    addedCands.push(cand);
  };
  const getL2ClusterStub = function (l2idx) {
    return scenario.clusters ? scenario.clusters[l2idx] : null;
  };
  // Also need _addL2SweepDismissed-style helpers — paste minimal versions
  // into the sandbox so the inspector block can call them.
  const dismissedHelpers = `
    const _L2_SWEEP_DISMISSED_KEY_PFX = 'pca_scrubber_v3.l2SweepDismissed.';
    function _loadL2SweepDismissed(chrom) {
      if (!chrom) return new Set();
      try {
        const raw = localStorage.getItem(_L2_SWEEP_DISMISSED_KEY_PFX + chrom);
        if (!raw) return new Set();
        const arr = JSON.parse(raw);
        if (!Array.isArray(arr)) return new Set();
        const out = new Set();
        for (const v of arr) {
          const n = Number(v);
          if (Number.isInteger(n) && n >= 0) out.add(n);
        }
        return out;
      } catch (_) { return new Set(); }
    }
    function _saveL2SweepDismissed(chrom, set) {
      if (!chrom || !(set instanceof Set)) return;
      try {
        localStorage.setItem(_L2_SWEEP_DISMISSED_KEY_PFX + chrom,
                             JSON.stringify(Array.from(set)));
      } catch (_) {}
    }
    function _addL2SweepDismissed(chrom, l2idx) {
      if (!chrom || !Number.isInteger(l2idx)) return;
      const s = _loadL2SweepDismissed(chrom);
      s.add(l2idx);
      _saveL2SweepDismissed(chrom, s);
    }
  `;
  const ctx = {
    state: stateStub,
    window: undefined,
    document: undefined,                          // disable DOM-ready wireup
    localStorage: localStorageStub,
    console: { warn: () => {}, log: () => {} },
    addCandidateToList: addCandidateToListStub,
    removeCandidateFully: removeCandidateFullyStub,
    removeCandidateFromList: removeCandidateFullyStub,
    getL2Cluster: getL2ClusterStub,
    runL2SweepInheritance: () => stateStub.l2SweepResult,
    Number, Array, Object, Math, JSON, Set, Date, Float64Array, Int8Array,
    isFinite, parseInt, parseFloat,
    setTimeout: () => {},
  };
  ctx.window = ctx;
  // Run the dismissed-set helpers first so the inspector block's calls
  // resolve.
  const context = vm.createContext(ctx);
  vm.runInContext(dismissedHelpers, context);
  vm.runInContext(inspectorBlock, context);
  return { ctx: context, addedCands, removedIds, ls: _ls, state: stateStub };
}

// 8a. _l2siBuildRows reads from state.l2SweepResult.l2_meta + candidateList
{
  console.log('\n--- 8a. _l2siBuildRows + status derivation ---');
  const sb = makeInspectorSandbox({
    state: {
      data: { chrom: 'LG28', l2_envelopes: [
        { _s0: 0, _e0: 10, start_bp: 1e6, end_bp: 2e6 },
        { _s0: 11, _e0: 20, start_bp: 5e6, end_bp: 6e6 },
        { _s0: 21, _e0: 30, start_bp: 9e6, end_bp: 10e6 },
      ] },
      candidateList: [
        { id: 'manual_a', l2_indices: [1] },
      ],
      candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: {
        l2_meta: [
          { l2idx: 0, item_idx: 0, seq_num: 1, start_bp: 1e6, end_bp: 2e6,
            silhouette: 0.6, n_per_group: [5,5,5], n_bands: 3 },
          { l2idx: 1, item_idx: 1, seq_num: 2, start_bp: 5e6, end_bp: 6e6,
            silhouette: 0.5, n_per_group: [5,5,5], n_bands: 3 },
          { l2idx: 2, item_idx: 2, seq_num: 3, start_bp: 9e6, end_bp: 10e6,
            silhouette: 0.7, n_per_group: [5,5,5], n_bands: 3 },
        ],
        groups: [
          { members: [{ item_idx: 0 }, { item_idx: 1 }] },
          { members: [{ item_idx: 0 }, { item_idx: 2 }] },
          { members: [{ item_idx: 1 }, { item_idx: 2 }] },
        ],
      },
    },
    ls: { 'pca_scrubber_v3.l2SweepDismissed.LG28': JSON.stringify([2]) },
  });
  const rows = sb.ctx._l2siBuildRows();
  ok('builds 3 rows', rows.length === 3);
  ok('L2:0 status = pending', rows[0].status === 'pending');
  ok('L2:1 status = promoted (in candidateList)', rows[1].status === 'promoted');
  ok('L2:2 status = dismissed (in localStorage set)', rows[2].status === 'dismissed');
  ok('L2:0 n_groups = 2 (touches groups #0 and #1)',
     rows[0].n_groups === 2, 'got ' + rows[0].n_groups);
  ok('L2:0 silhouette preserved', Math.abs(rows[0].silhouette - 0.6) < 1e-9);
  ok('L2:0 n_per_group preserved as array',
     Array.isArray(rows[0].n_per_group) && rows[0].n_per_group.length === 3);
}

// 8b. _l2siPromoteRow adds to candidateList with auto_l2_sweep tag
{
  console.log('\n--- 8b. _l2siPromoteRow gate-bypass behavior ---');
  const sb = makeInspectorSandbox({
    state: {
      data: { chrom: 'LG28', l2_envelopes: [
        { _s0: 0, _e0: 10, start_bp: 1e6, end_bp: 2e6 },
      ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: null,
    },
    clusters: {
      0: { fixedKLabels: new Int8Array([0,0,1,1,2,2]) },
    },
  });
  const r1 = sb.ctx._l2siPromoteRow(0);
  ok('first promote returns true', r1 === true);
  ok('candidate added to list', sb.state.candidateList.length === 1);
  ok('candidate has source=auto_l2_sweep',
     sb.addedCands[0].source === 'auto_l2_sweep');
  ok('candidate has confirmed=false',
     sb.addedCands[0].confirmed === false);
  ok('candidate l2_indices = [0]',
     Array.isArray(sb.addedCands[0].l2_indices) &&
     sb.addedCands[0].l2_indices[0] === 0);
  ok('candidate notes mentions inspector',
     /inspector/i.test(sb.addedCands[0].notes));
  // Re-promote → false (already in list)
  const r2 = sb.ctx._l2siPromoteRow(0);
  ok('re-promote returns false', r2 === false);
  ok('candidateList still 1', sb.state.candidateList.length === 1);
}

// 8c. _l2siDismissRow adds to dismissed set + drops auto candidate
{
  console.log('\n--- 8c. _l2siDismissRow behavior ---');
  const sb = makeInspectorSandbox({
    state: {
      data: { chrom: 'LG28', l2_envelopes: [
        { _s0: 0, _e0: 10, start_bp: 1e6, end_bp: 2e6 },
        { _s0: 11, _e0: 20, start_bp: 5e6, end_bp: 6e6 },
      ] },
      candidateList: [
        { id: 'manual_a', source: null,            confirmed: true,  l2_indices: [0] },
        { id: 'auto_a',   source: 'auto_l2_sweep', confirmed: false, l2_indices: [1] },
      ],
      candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: null,
    },
  });
  // Dismiss L2:1 (auto-promoted)
  const r1 = sb.ctx._l2siDismissRow(1);
  ok('dismiss L2:1 returns true', r1 === true);
  ok('auto candidate removed', sb.removedIds.includes('auto_a'));
  ok('manual candidate still present',
     sb.state.candidateList.some(c => c.id === 'manual_a'));
  // Dismissed set persisted
  const stored = JSON.parse(sb.ls['pca_scrubber_v3.l2SweepDismissed.LG28']);
  ok('dismissed set persisted with l2idx=1',
     Array.isArray(stored) && stored.indexOf(1) !== -1);

  // Dismissing L2:0 (manual candidate) — should add to dismissed set
  // but NOT touch the manual candidate.
  const r2 = sb.ctx._l2siDismissRow(0);
  ok('dismiss L2:0 returns true', r2 === true);
  ok('manual candidate STILL present (not dropped by dismiss)',
     sb.state.candidateList.some(c => c.id === 'manual_a'));
  const stored2 = JSON.parse(sb.ls['pca_scrubber_v3.l2SweepDismissed.LG28']);
  ok('dismissed set now has both 0 and 1',
     stored2.indexOf(0) !== -1 && stored2.indexOf(1) !== -1);
}

// 8d. _l2siFilterState lazy-init on state.l2SweepInspectorFilter
{
  console.log('\n--- 8d. _l2siFilterState lazy init ---');
  const sb = makeInspectorSandbox({
    state: {
      data: null, candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: null,
    },
  });
  ok('filter slot starts undefined',
     sb.state.l2SweepInspectorFilter === undefined);
  const f1 = sb.ctx._l2siFilterState();
  ok('first call creates the slot', !!sb.state.l2SweepInspectorFilter);
  ok('default minSil = 0', f1.minSil === 0);
  ok("default status = 'all'", f1.status === 'all');
  // Mutating the returned object reflects on state (object ref)
  f1.minSil = 0.4;
  const f2 = sb.ctx._l2siFilterState();
  ok('second call returns same ref (mutation persists)',
     f2.minSil === 0.4);
}

// 8e. _refreshL2SweepInspectBtnAvailability — guard on state.l2SweepResult
{
  console.log('\n--- 8e. inspect-btn availability gating ---');
  // Need a fake document with the button
  const fakeBtn = {
    id: 'l2SweepInspectBtn',
    disabled: true,
    style: {},
    dataset: {},
    addEventListener() {},
  };
  const fakeDoc = {
    getElementById(id) { return id === 'l2SweepInspectBtn' ? fakeBtn : null; },
  };
  const sb = makeInspectorSandbox({
    state: {
      data: { chrom: 'LG28', l2_envelopes: [] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: null,
    },
  });
  // Re-context with document populated and re-eval the inspector block so
  // _refreshL2SweepInspectBtnAvailability picks up the document.
  const ctx2 = {
    state: sb.state,
    window: undefined,
    document: fakeDoc,
    localStorage: { getItem: () => null, setItem: () => {} },
    console: { warn: () => {}, log: () => {} },
    addCandidateToList: () => {},
    removeCandidateFully: () => {},
    removeCandidateFromList: () => {},
    getL2Cluster: () => null,
    runL2SweepInheritance: () => null,
    Number, Array, Object, Math, JSON, Set, Date, Float64Array, Int8Array,
    isFinite, parseInt, parseFloat,
    setTimeout: () => {},
  };
  ctx2.window = ctx2;
  const c2 = vm.createContext(ctx2);
  vm.runInContext(`
    function _loadL2SweepDismissed() { return new Set(); }
    function _saveL2SweepDismissed() {}
    function _addL2SweepDismissed() {}
  `, c2);
  vm.runInContext(inspectorBlock, c2);
  // No sweep result → button disabled
  c2._refreshL2SweepInspectBtnAvailability();
  ok('button disabled when no sweep result', fakeBtn.disabled === true);
  ok('button cursor not-allowed when disabled',
     fakeBtn.style.cursor === 'not-allowed');
  // Populate sweep result → button enabled
  sb.state.l2SweepResult = { l2_meta: [] };
  c2._refreshL2SweepInspectBtnAvailability();
  ok('button enabled when sweep result present', fakeBtn.disabled === false);
  ok('button cursor pointer when enabled',
     fakeBtn.style.cursor === 'pointer');
}

// ============================================================================
// 9. Modal renderer — DOM rendering test (uses fake document)
// ============================================================================
console.log('\n=== 9. Modal renderer DOM output ===');

function makeFullDom(overlayInitDisplay) {
  const overlay = {
    innerHTML: '',
    style: { display: overlayInitDisplay || 'none' },
    addEventListener() {},
    querySelectorAll() { return []; },
  };
  const elements = { l2SweepInspectorOverlay: overlay };
  return {
    overlay,
    document: {
      getElementById(id) { return elements[id] || null; },
      createElement() { return { addEventListener() {} }; },
      addEventListener() {},
      removeEventListener() {},
      body: { appendChild() {} },
    },
  };
}

// 9a. Empty sweep result → renders empty-state message
{
  console.log('\n--- 9a. Empty sweep result message ---');
  const dom = makeFullDom('none');
  const ctx2 = {
    state: {
      data: { chrom: 'LG28', l2_envelopes: [] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: { l2_meta: [], groups: [] },
    },
    window: undefined,
    document: dom.document,
    localStorage: { getItem: () => null, setItem: () => {} },
    console: { warn: () => {}, log: () => {} },
    addCandidateToList: () => {}, removeCandidateFully: () => {},
    removeCandidateFromList: () => {}, getL2Cluster: () => null,
    runL2SweepInheritance: () => null,
    Number, Array, Object, Math, JSON, Set, Date, Float64Array, Int8Array,
    isFinite, parseInt, parseFloat, setTimeout: () => {},
  };
  ctx2.window = ctx2;
  const c2 = vm.createContext(ctx2);
  vm.runInContext(`
    function _loadL2SweepDismissed() { return new Set(); }
    function _saveL2SweepDismissed() {}
    function _addL2SweepDismissed() {}
  `, c2);
  vm.runInContext(inspectorBlock, c2);
  c2._renderL2SweepInspectorModal();
  ok('overlay set to display:flex',
     dom.overlay.style.display === 'flex');
  ok('rendered HTML contains chrom name',
     dom.overlay.innerHTML.indexOf('LG28') !== -1);
  ok('rendered HTML contains "0 usable L2s" header',
     /0 usable L2s/.test(dom.overlay.innerHTML));
  ok('rendered HTML contains empty-state message',
     /No usable L2s in the current sweep result/.test(dom.overlay.innerHTML));
  ok('rendered HTML has a close button',
     /id="l2SweepInspectClose"/.test(dom.overlay.innerHTML));
  ok('rendered HTML has filter inputs',
     /id="l2siFilterMinSil"/.test(dom.overlay.innerHTML) &&
     /id="l2siFilterMinGrps"/.test(dom.overlay.innerHTML) &&
     /id="l2siFilterStatus"/.test(dom.overlay.innerHTML));
  ok('bulk-action buttons present',
     /id="l2siBulkPromote"/.test(dom.overlay.innerHTML) &&
     /id="l2siBulkDismiss"/.test(dom.overlay.innerHTML) &&
     /id="l2siResweep"/.test(dom.overlay.innerHTML));
}

// 9b. Populated sweep → renders one row per L2_meta entry
{
  console.log('\n--- 9b. Populated sweep → row rendering ---');
  const dom = makeFullDom('none');
  const ctx2 = {
    state: {
      data: { chrom: 'LG28', l2_envelopes: [
        { _s0: 0, _e0: 10, start_bp: 1e6, end_bp: 2e6 },
        { _s0: 11, _e0: 20, start_bp: 5e6, end_bp: 6e6 },
      ] },
      candidateList: [], candidates: {}, k: 3, activeMode: 'default',
      l2SweepResult: {
        l2_meta: [
          { l2idx: 0, item_idx: 0, seq_num: 1, start_bp: 1e6, end_bp: 2e6,
            silhouette: 0.6, n_per_group: [5,5,5], n_bands: 3 },
          { l2idx: 1, item_idx: 1, seq_num: 2, start_bp: 5e6, end_bp: 6e6,
            silhouette: 0.4, n_per_group: [5,5,5], n_bands: 3 },
        ],
        groups: [
          { members: [{ item_idx: 0 }, { item_idx: 1 }] },
          { members: [{ item_idx: 0 }] },
          { members: [{ item_idx: 1 }] },
        ],
      },
    },
    window: undefined,
    document: dom.document,
    localStorage: { getItem: () => null, setItem: () => {} },
    console: { warn: () => {}, log: () => {} },
    addCandidateToList: () => {}, removeCandidateFully: () => {},
    removeCandidateFromList: () => {}, getL2Cluster: () => null,
    runL2SweepInheritance: () => null,
    Number, Array, Object, Math, JSON, Set, Date, Float64Array, Int8Array,
    isFinite, parseInt, parseFloat, setTimeout: () => {},
  };
  ctx2.window = ctx2;
  const c2 = vm.createContext(ctx2);
  vm.runInContext(`
    function _loadL2SweepDismissed() { return new Set(); }
    function _saveL2SweepDismissed() {}
    function _addL2SweepDismissed() {}
  `, c2);
  vm.runInContext(inspectorBlock, c2);
  c2._renderL2SweepInspectorModal();
  const innerHTML = dom.overlay.innerHTML;
  ok('contains "2 usable L2s" header', /2 usable L2s/.test(innerHTML));
  ok('contains range "1.00–2.00 Mb"', /1\.00–2\.00 Mb/.test(innerHTML));
  ok('contains range "5.00–6.00 Mb"', /5\.00–6\.00 Mb/.test(innerHTML));
  ok('contains promote button per row', (innerHTML.match(/_l2siPromoteBtn/g) || []).length === 2);
  ok('contains dismiss button per row', (innerHTML.match(/_l2siDismissBtn/g) || []).length === 2);
  ok('contains seq numbers I1 and I2',
     /I1\s*<\/div>/.test(innerHTML) && /I2\s*<\/div>/.test(innerHTML));
  ok('contains silhouette 0.60 and 0.40',
     innerHTML.indexOf('0.60') !== -1 && innerHTML.indexOf('0.40') !== -1);
  ok('contains "⚪ pending" status chips for both rows',
     (innerHTML.match(/⚪ pending/g) || []).length === 2);
}

// ============================================================================
// 10. Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
