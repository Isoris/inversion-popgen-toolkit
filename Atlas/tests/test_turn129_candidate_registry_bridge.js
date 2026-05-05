// =============================================================================
// turn 129 test — _rebuildCandidateRegistries() bridge
//
// Quentin's bug report (turn 128d):
//   "The grouping panel is bugged — it doesn't auto select or compute when
//    we click it for the candidate. It doesn't find the 'groups' also
//    shouldn't it be related or linked to the groups in candidate page as
//    well? Or its a different grouping?"
//
// Diagnosis (handoff §3): two divergent registries.
//   - state.candidateList (array) ← every user action mutates this
//   - state.candidates    (dict)  ← only Save-Session import populates this
//   - _gatherActiveCandidatesForInheritance() reads state.candidates → empty
//   - runInheritanceCompute() bails with items.length < 2 → I·g pills empty
//
// Fix (option c, "auto-sync"): on every mutation of state.candidateList,
// rebuild state.candidates as a derived index. Single helper hooked into:
//   - persistCandidateList()  (~5 mutation entry points)
//   - loadCandidateList()     (F5 / cold start)
//   - commitL3Draft()         (draft → candidateList)
//
// Tests verify:
//   1.  _rebuildCandidateRegistries is defined and exposed on window.
//   2.  Bridge populates state.candidates (dict) from state.candidateList (array).
//   3.  Bridge populates state.candidates_detailed from state.candidateList_detailed.
//   4.  Bridge replaces the dict (doesn't merge): removed candidates drop out.
//   5.  Bridge handles empty / undefined arrays gracefully.
//   6.  After bridge, _gatherActiveCandidatesForInheritance returns the array's
//       contents (unblocks the I·g compute).
//   7.  Single-candidate gather → inheritance returns null (need ≥ 2; verified
//       by reading the gatherer's contract).
//   8.  persistCandidateList wires the bridge.
//   9.  loadCandidateList wires the bridge.
//   10. commitL3Draft now persists (was broken pre-129) which transitively
//       triggers the bridge.
//   11. Session-load merges payload.candidates into state.candidateList so
//       the bridge can't subsequently wipe the loaded data.
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

ok('_rebuildCandidateRegistries function is defined',
   /function _rebuildCandidateRegistries\(\)/.test(html));

ok('_rebuildCandidateRegistries exposed on window',
   /window\._rebuildCandidateRegistries\s*=\s*_rebuildCandidateRegistries/.test(html));

ok('persistCandidateList calls _rebuildCandidateRegistries',
   /function persistCandidateList\(\)[\s\S]{0,800}_rebuildCandidateRegistries\(\)/.test(html));

ok('loadCandidateList calls _rebuildCandidateRegistries',
   /function loadCandidateList\(\)[\s\S]{0,3000}_rebuildCandidateRegistries\(\)/.test(html));

ok('commitL3Draft calls persistCandidateList',
   /function commitL3Draft\(\)[\s\S]{0,2500}persistCandidateList\(\)/.test(html));

ok('session-load handler folds payload.candidates into candidateList',
   /payload\.candidates[\s\S]{0,600}state\.candidateList\.push\(c\)/.test(html),
   'Save-Session import should now populate the array, not just the dict');

ok('session-load also calls _rebuildCandidateRegistries',
   /payload\.candidates[\s\S]{0,800}_rebuildCandidateRegistries\(\)/.test(html));

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

const fnRebuild  = pullFunction(html, '_rebuildCandidateRegistries');
const fnGather   = pullFunction(html, '_gatherActiveCandidatesForInheritance');
ok('_rebuildCandidateRegistries extractable', !!fnRebuild);
ok('_gatherActiveCandidatesForInheritance extractable', !!fnGather);

function makeSandbox(initialState) {
  const sandbox = { state: initialState, console };
  vm.createContext(sandbox);
  vm.runInContext(fnRebuild, sandbox);
  vm.runInContext(fnGather, sandbox);
  return sandbox;
}

function makeCandidate(id, start_bp, end_bp, K) {
  return {
    id: id,
    start_bp: start_bp,
    end_bp: end_bp,
    K: K || 3,
    locked_labels: new Int8Array([0, 1, 2, 0, 1]),
  };
}

// --- Test 1: bridge populates dict from array
console.log('\nTest 1: bridge populates dict from array');
{
  const sandbox = makeSandbox({
    candidateList: [
      makeCandidate('c1', 1000, 2000),
      makeCandidate('c2', 5000, 6000),
      makeCandidate('c3', 9000, 10000),
    ],
  });
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('rebuild returned default_n=3', sandbox.__r.default_n === 3,
     'got: ' + JSON.stringify(sandbox.__r));
  ok('state.candidates is a dict with 3 entries',
     typeof sandbox.state.candidates === 'object' &&
     Object.keys(sandbox.state.candidates).length === 3);
  ok('state.candidates["c1"] is the c1 candidate',
     sandbox.state.candidates.c1 && sandbox.state.candidates.c1.id === 'c1');
  ok('state.candidates["c2"] is the c2 candidate',
     sandbox.state.candidates.c2 && sandbox.state.candidates.c2.id === 'c2');
  ok('state.candidates["c3"] is the c3 candidate',
     sandbox.state.candidates.c3 && sandbox.state.candidates.c3.id === 'c3');
}

// --- Test 2: bridge replaces (doesn't merge) — removed candidates drop out
console.log('\nTest 2: bridge replaces — removed candidates drop out');
{
  const sandbox = makeSandbox({
    candidateList: [makeCandidate('c1', 1000, 2000), makeCandidate('c2', 5000, 6000)],
    // Pre-existing stale entry that should be dropped (not in array)
    candidates: { c_stale: makeCandidate('c_stale', 7000, 8000) },
  });
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  ok('stale c_stale is dropped',
     sandbox.state.candidates.c_stale === undefined);
  ok('current c1 is present',
     sandbox.state.candidates.c1 && sandbox.state.candidates.c1.id === 'c1');
  ok('current c2 is present',
     sandbox.state.candidates.c2 && sandbox.state.candidates.c2.id === 'c2');
  ok('dict has exactly 2 entries (matches array)',
     Object.keys(sandbox.state.candidates).length === 2);
}

// --- Test 3: bridge handles detailed-mode array
console.log('\nTest 3: bridge handles detailed-mode array');
{
  const sandbox = makeSandbox({
    candidateList: [makeCandidate('cd', 100, 200)],
    candidateList_detailed: [makeCandidate('cd', 100, 200), makeCandidate('cd2', 300, 400)],
  });
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('rebuild returned default_n=1',  sandbox.__r.default_n  === 1);
  ok('rebuild returned detailed_n=2', sandbox.__r.detailed_n === 2);
  ok('candidates_detailed has 2 entries',
     Object.keys(sandbox.state.candidates_detailed).length === 2);
  ok('candidates (default) has 1 entry',
     Object.keys(sandbox.state.candidates).length === 1);
}

// --- Test 4: bridge handles empty / undefined arrays
console.log('\nTest 4: bridge handles empty / undefined arrays');
{
  const sandbox = makeSandbox({});  // no candidateList at all
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('no crash on undefined arrays', sandbox.__r.default_n === 0);
  ok('state.candidates is empty object', typeof sandbox.state.candidates === 'object' &&
     Object.keys(sandbox.state.candidates).length === 0);
  ok('state.candidates_detailed is empty object', typeof sandbox.state.candidates_detailed === 'object' &&
     Object.keys(sandbox.state.candidates_detailed).length === 0);
}
{
  const sandbox = makeSandbox({ candidateList: [] });  // empty array
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('no crash on empty array', sandbox.__r.default_n === 0);
}
{
  const sandbox = makeSandbox({ candidateList: [null, undefined, { id: '' }, makeCandidate('c1', 1, 2)] });
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('skips falsy items (null/undefined)', sandbox.__r.default_n === 1);
  ok('skips empty-id items', !sandbox.state.candidates['']);
  ok('keeps valid items', !!sandbox.state.candidates.c1);
}

// --- Test 5: bridge unblocks the gatherer (the actual bug-fix proof)
console.log('\nTest 5: bridge unblocks _gatherActiveCandidatesForInheritance');
{
  const sandbox = makeSandbox({
    activeMode: 'default',
    candidateList: [
      makeCandidate('c1', 1000, 2000),
      makeCandidate('c2', 5000, 6000),
      makeCandidate('c3', 9000, 10000),
    ],
    k: 3,
  });

  // BEFORE bridge: state.candidates is undefined → gather returns []
  vm.runInContext('var __preBridge = _gatherActiveCandidatesForInheritance();', sandbox);
  ok('pre-bridge gather returns empty (the bug)',
     Array.isArray(sandbox.__preBridge) && sandbox.__preBridge.length === 0,
     'this is the broken state from turn 128d');

  // AFTER bridge: state.candidates is populated → gather returns 3
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  vm.runInContext('var __postBridge = _gatherActiveCandidatesForInheritance();', sandbox);
  ok('post-bridge gather returns 3 candidates (the fix)',
     Array.isArray(sandbox.__postBridge) && sandbox.__postBridge.length === 3);
  ok('post-bridge candidates are sorted by start_bp',
     sandbox.__postBridge[0].id === 'c1' &&
     sandbox.__postBridge[1].id === 'c2' &&
     sandbox.__postBridge[2].id === 'c3');
  ok('post-bridge candidates have seq_num assigned (I1/I2/I3)',
     sandbox.__postBridge[0].seq_num === 1 &&
     sandbox.__postBridge[1].seq_num === 2 &&
     sandbox.__postBridge[2].seq_num === 3);
}

// --- Test 6: gatherer respects activeMode after bridge
console.log('\nTest 6: bridge respects activeMode (default vs detailed)');
{
  const sandbox = makeSandbox({
    activeMode: 'detailed',
    candidateList: [makeCandidate('c_def', 1000, 2000)],
    candidateList_detailed: [makeCandidate('c_det1', 100, 200), makeCandidate('c_det2', 300, 400)],
  });
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  vm.runInContext('var __out = _gatherActiveCandidatesForInheritance();', sandbox);
  ok('detailed-mode gather returns detailed candidates',
     sandbox.__out.length === 2 &&
     sandbox.__out[0].id === 'c_det1' &&
     sandbox.__out[1].id === 'c_det2');
}

// --- Test 7: single-candidate path → still null (gatherer's contract: < 2 ⇒ null)
console.log('\nTest 7: single-candidate inheritance path returns 1 (caller bails)');
{
  const sandbox = makeSandbox({
    activeMode: 'default',
    candidateList: [makeCandidate('only', 1000, 2000)],
  });
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  vm.runInContext('var __out = _gatherActiveCandidatesForInheritance();', sandbox);
  ok('1 candidate gathered (runInheritanceCompute will bail downstream)',
     sandbox.__out.length === 1,
     'gatherer returns the 1; caller checks items.length < 2 and returns null');
}

// --- Test 8: idempotency — calling rebuild twice is a no-op the second time
console.log('\nTest 8: idempotency');
{
  const sandbox = makeSandbox({
    candidateList: [makeCandidate('c1', 1000, 2000), makeCandidate('c2', 5000, 6000)],
  });
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  const firstSnapshot = JSON.parse(JSON.stringify(
    Object.keys(sandbox.state.candidates).map(k => sandbox.state.candidates[k].id)
  ));
  vm.runInContext('_rebuildCandidateRegistries();', sandbox);
  const secondSnapshot = JSON.parse(JSON.stringify(
    Object.keys(sandbox.state.candidates).map(k => sandbox.state.candidates[k].id)
  ));
  ok('first and second rebuild produce same dict',
     JSON.stringify(firstSnapshot) === JSON.stringify(secondSnapshot));
}

// --- Test 9: window.state vs module state — bridge prefers window.state
console.log('\nTest 9: bridge uses window.state when present');
{
  const sandbox = {
    state: { candidateList: [makeCandidate('module_c', 1, 2)] },
    window: { state: { candidateList: [makeCandidate('window_c', 100, 200)] } },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnRebuild, sandbox);
  vm.runInContext('var __r = _rebuildCandidateRegistries();', sandbox);
  ok('bridge wrote to window.state.candidates',
     sandbox.window.state.candidates && sandbox.window.state.candidates.window_c);
  ok('bridge did NOT write to module-state.candidates',
     !sandbox.state.candidates ||
     !sandbox.state.candidates.module_c);
}

// =============================================================================
// Final tally
// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
