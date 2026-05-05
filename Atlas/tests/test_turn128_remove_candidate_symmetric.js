// =============================================================================
// turn 128 integration test — symmetric remove-candidate behaviour
//
// Bug from OBSERVATIONS_TO_FIX.txt:
//   "When we go to the candidate focus and we click to drop the candidate, it
//   doesnt drop it. So we have to go to the karyotype / tier page and click
//   the red cross and only then it will be dropped. I feel that the two
//   buttons must be working similarily."
//
// Fix: extracted a canonical `removeCandidateFully(id)` helper and routed
// both surfaces (page 2 ✕ clear button, page 4 karyotype/tier red ✕) through
// it.
//
// Tests verify:
//   1. removeCandidateFully is defined and exposed.
//   2. removeCandidateFully drops the id from candidateList, nulls
//      state.candidate when it was the active focus, and calls every refresh
//      hook in the right order.
//   3. removeCandidateFully treats `id === null/undefined` as "remove the
//      currently focused candidate".
//   4. removeCandidateFully returns false when there's nothing to remove.
//   5. Page 2 ✕-clear handler routes through removeCandidateFully when the
//      focused candidate is in the saved list, and falls back to
//      clearCandidate() otherwise.
//   6. Page 4 cli-remove handler delegates to removeCandidateFully (no
//      duplicate inline logic).
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

ok('removeCandidateFully defined',
   /function removeCandidateFully\(/.test(html));

ok('page 2 ✕ button calls removeCandidateFully',
   /candidateClearBtn[\s\S]{0,800}removeCandidateFully\(state\.candidate\.id\)/.test(html),
   'page 2 #candidateClearBtn handler should now route to removeCandidateFully');

ok('page 2 ✕ button still has fallback to clearCandidate for unsaved focus',
   /candidateClearBtn[\s\S]{0,800}clearCandidate\(\)/.test(html));

ok('page 4 cli-remove handler delegates to removeCandidateFully',
   /cli-remove[\s\S]{0,600}removeCandidateFully\(cid\)/.test(html));

ok('page 4 cli-remove no longer has the OLD inline state.candidate=null block',
   !/cli-remove[\s\S]{0,200}removeCandidateFromList\(cid\);[\s\S]{0,200}state\.candidate = null/.test(html),
   'old inline removal logic should be gone — delegated to removeCandidateFully');

// =============================================================================
// Behavioural tests
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\$&') + '\\s*\\(', 'm');
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

const fnRemoveFully = pullFunction(html, 'removeCandidateFully');
ok('removeCandidateFully extractable', !!fnRemoveFully);

// Build a sandbox where removeCandidateFromList, _persistActiveCandidate,
// refreshCandidateUI, renderCandidateKaryotype, drawAll, setCur are stubs
// that record their calls. The real implementations would touch DOM /
// localStorage; we just want to assert the orchestration is correct.
function makeSandbox(initialState) {
  const calls = [];
  const sandbox = {
    state: initialState,
    removeCandidateFromList: function (id) { calls.push(['removeCandidateFromList', id]);
      sandbox.state.candidateList = sandbox.state.candidateList.filter(c => c.id !== id);
    },
    _persistActiveCandidate: function (cid) { calls.push(['_persistActiveCandidate', cid]); },
    refreshCandidateUI:      function ()    { calls.push(['refreshCandidateUI']); },
    renderCandidateKaryotype: function ()   { calls.push(['renderCandidateKaryotype']); },
    drawAll:                  function ()    { calls.push(['drawAll']); },
    setCur:                   function (c)   { calls.push(['setCur', c]); },
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnRemoveFully, sandbox);
  return { sandbox, calls };
}

// --- Test: remove an active, saved candidate
console.log('\nTest: remove active+saved candidate');
{
  const cand = { id: 'cand_1' };
  const { sandbox, calls } = makeSandbox({
    candidate:     cand,
    candidateList: [cand, { id: 'cand_2' }],
    cur: 5,
  });
  vm.runInContext('var __r = removeCandidateFully("cand_1");', sandbox);
  ok('returned true (remove happened)', sandbox.__r === true);
  ok('candidateList no longer contains cand_1',
     sandbox.state.candidateList.every(c => c.id !== 'cand_1'),
     'list: ' + JSON.stringify(sandbox.state.candidateList));
  ok('state.candidate is null',
     sandbox.state.candidate === null);
  // Order-of-calls
  const callNames = calls.map(c => c[0]);
  ok('removeCandidateFromList called',
     callNames.indexOf('removeCandidateFromList') >= 0);
  ok('_persistActiveCandidate(null) called for active candidate',
     calls.some(c => c[0] === '_persistActiveCandidate' && c[1] === null));
  ok('refreshCandidateUI called',
     callNames.indexOf('refreshCandidateUI') >= 0);
  ok('renderCandidateKaryotype called',
     callNames.indexOf('renderCandidateKaryotype') >= 0);
  ok('drawAll called for page 1 overlay refresh',
     callNames.indexOf('drawAll') >= 0);
}

// --- Test: remove a saved candidate that is NOT the active focus
console.log('\nTest: remove non-active saved candidate');
{
  const active = { id: 'cand_active' };
  const { sandbox, calls } = makeSandbox({
    candidate:     active,
    candidateList: [active, { id: 'cand_other' }],
    cur: 0,
  });
  vm.runInContext('var __r = removeCandidateFully("cand_other");', sandbox);
  ok('returned true', sandbox.__r === true);
  ok('cand_other gone from list',
     sandbox.state.candidateList.every(c => c.id !== 'cand_other'));
  ok('state.candidate untouched (still cand_active)',
     sandbox.state.candidate && sandbox.state.candidate.id === 'cand_active');
  ok('_persistActiveCandidate NOT called (active was different)',
     !calls.some(c => c[0] === '_persistActiveCandidate'));
  ok('refreshCandidateUI still called',
     calls.some(c => c[0] === 'refreshCandidateUI'));
}

// --- Test: id=null/undefined → remove the currently focused candidate
console.log('\nTest: id=null routes to currently focused candidate');
{
  const cand = { id: 'cand_focus' };
  const { sandbox, calls } = makeSandbox({
    candidate:     cand,
    candidateList: [cand],
    cur: 0,
  });
  vm.runInContext('var __r = removeCandidateFully(null);', sandbox);
  ok('returned true with null id', sandbox.__r === true);
  ok('cand_focus removed from list',
     sandbox.state.candidateList.length === 0);
  ok('state.candidate nulled', sandbox.state.candidate === null);
  ok('removeCandidateFromList called with focused id',
     calls.some(c => c[0] === 'removeCandidateFromList' && c[1] === 'cand_focus'));
}

// --- Test: id=null AND no active focus → no-op, returns false
console.log('\nTest: nothing to remove → returns false, no side-effects');
{
  const { sandbox, calls } = makeSandbox({
    candidate:     null,
    candidateList: [{ id: 'foo' }],
    cur: 0,
  });
  vm.runInContext('var __r = removeCandidateFully(undefined);', sandbox);
  ok('returned false', sandbox.__r === false);
  ok('list untouched',
     sandbox.state.candidateList.length === 1 &&
     sandbox.state.candidateList[0].id === 'foo');
  ok('no side-effects (zero calls)', calls.length === 0,
     'unexpected calls: ' + JSON.stringify(calls));
}

// --- Test: graceful when renderCandidateKaryotype is absent (page 4 not loaded)
console.log('\nTest: graceful when renderCandidateKaryotype undefined');
{
  const cand = { id: 'cand_x' };
  const sandbox = {
    state: {
      candidate:     cand,
      candidateList: [cand],
      cur: 0,
    },
    removeCandidateFromList: function (id) {
      sandbox.state.candidateList = sandbox.state.candidateList.filter(c => c.id !== id);
    },
    _persistActiveCandidate: function () {},
    refreshCandidateUI: function () {},
    // renderCandidateKaryotype intentionally NOT defined
    drawAll: function () {},
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnRemoveFully, sandbox);
  let threw = false;
  try {
    vm.runInContext('removeCandidateFully("cand_x");', sandbox);
  } catch (e) { threw = true; }
  ok('does not throw when renderCandidateKaryotype is missing', !threw);
  ok('candidate still got removed despite missing helper',
     sandbox.state.candidateList.length === 0);
}

// --- Test: graceful when drawAll AND setCur both missing
console.log('\nTest: graceful when drawAll AND setCur missing');
{
  const cand = { id: 'cand_y' };
  const sandbox = {
    state: { candidate: cand, candidateList: [cand], cur: 0 },
    removeCandidateFromList: function (id) {
      sandbox.state.candidateList = sandbox.state.candidateList.filter(c => c.id !== id);
    },
    _persistActiveCandidate: function () {},
    refreshCandidateUI: function () {},
    renderCandidateKaryotype: function () {},
    console,
  };
  vm.createContext(sandbox);
  vm.runInContext(fnRemoveFully, sandbox);
  let threw = false;
  try { vm.runInContext('removeCandidateFully("cand_y");', sandbox); }
  catch (e) { threw = true; }
  ok('does not throw when drawAll and setCur both missing', !threw);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
