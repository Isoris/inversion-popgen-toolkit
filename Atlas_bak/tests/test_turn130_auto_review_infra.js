// =============================================================================
// turn 130 follow-up — Auto-review infrastructure (Slice 0 of
// SPEC_review_surfaces_auto_and_lineages.md)
//
// Slice 0 lays the data-flow plumbing for auto-promoted candidates so that
// when L2-sweep ships its auto-promote mechanic, the review surfaces work
// out-of-the-box. This slice ships:
//   - _isAutoCandidate(cand) predicate as the central source of truth
//   - _gatherActiveCandidatesForInheritance skips auto candidates (so I·g
//     pills + cross-candidate matrix don't include unreviewed proposals)
//   - .is-auto CSS class on candidate list rows + 🤖 prefix
//   - Sort modifier in refreshCandidateListUI sorts auto candidates last
//
// No producers exist yet (L2-sweep ships next). Tests use synthetic
// candidates with cand.source = 'auto_l2_sweep' to validate the plumbing.
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

ok('_isAutoCandidate function defined',
   /function _isAutoCandidate\(cand\)/.test(html));

ok('_isAutoCandidate exposed on window',
   /window\._isAutoCandidate\s*=\s*_isAutoCandidate/.test(html));

ok('_isAutoCandidate documents auto_l2_sweep source tag',
   /'auto_l2_sweep'/.test(html));

ok('_isAutoCandidate documents auto_lineage source tag',
   /'auto_lineage'/.test(html));

ok('_gatherActiveCandidatesForInheritance skips auto candidates',
   /_gatherActiveCandidatesForInheritance[\s\S]{0,2500}_isAutoCandidate\(c\)\)\s*continue/.test(html));

ok('_gatherActiveCandidatesForInheritance has defensive fallback for sandboxes',
   /_gatherActiveCandidatesForInheritance[\s\S]{0,3000}typeof _isAutoCandidate === 'function'[\s\S]{0,400}c\.source\.indexOf\('auto_'\) === 0/.test(html));

ok('refreshCandidateListUI sorts auto candidates to the bottom',
   /refreshCandidateListUI[\s\S]{0,4000}_isAutoCandidate[\s\S]{0,200}aAuto !== bAuto/.test(html));

ok('refreshCandidateListUI emits .is-auto class on auto rows',
   /refreshCandidateListUI[\s\S]{0,5500}isAuto \? ' is-auto' : ''/.test(html));

ok('refreshCandidateListUI emits 🤖 prefix on auto rows',
   /refreshCandidateListUI[\s\S]{0,5500}cli-auto-prefix[\s\S]{0,200}🤖/.test(html));

ok('CSS rule for .cand-list-item.is-auto present',
   /#candListContainer \.cand-list-item\.is-auto[\s\S]{0,300}border-style: dashed/.test(html));

ok('CSS rule for .cli-auto-prefix present',
   /\.cli-auto-prefix[\s\S]{0,200}font-size: 11px/.test(html));

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

const fnIsAuto       = pullFunction(html, '_isAutoCandidate');
const fnGather       = pullFunction(html, '_gatherActiveCandidatesForInheritance');

ok('_isAutoCandidate extractable',                          !!fnIsAuto);
ok('_gatherActiveCandidatesForInheritance extractable',     !!fnGather);

function makeSandbox(stateOverride) {
  const sandbox = { state: stateOverride || {}, console };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fnIsAuto,    sandbox);
  vm.runInContext(fnGather,    sandbox);
  return sandbox;
}

// --- Test 1: _isAutoCandidate basic predicate
console.log('\nTest 1: _isAutoCandidate predicate');
{
  const sandbox = makeSandbox({});
  vm.runInContext(`
    var r1 = _isAutoCandidate(null);
    var r2 = _isAutoCandidate(undefined);
    var r3 = _isAutoCandidate({});
    var r4 = _isAutoCandidate({ source: 'manual' });
    var r5 = _isAutoCandidate({ source: 'auto_l2_sweep' });
    var r6 = _isAutoCandidate({ source: 'auto_lineage' });
    var r7 = _isAutoCandidate({ source: 'auto_l2_sweep', confirmed: true });
    var r8 = _isAutoCandidate({ source: '' });
    var r9 = _isAutoCandidate({ source: 42 });
  `, sandbox);
  ok('null → false',                              sandbox.r1 === false);
  ok('undefined → false',                         sandbox.r2 === false);
  ok('empty object → false',                      sandbox.r3 === false);
  ok('source=manual → false',                     sandbox.r4 === false);
  ok('source=auto_l2_sweep → true',               sandbox.r5 === true);
  ok('source=auto_lineage → true',                sandbox.r6 === true);
  ok('confirmed auto candidate → false',          sandbox.r7 === false);
  ok('empty source → false',                      sandbox.r8 === false);
  ok('non-string source → false',                 sandbox.r9 === false);
}

// --- Test 2: gather() skips auto candidates
console.log('\nTest 2: _gatherActiveCandidatesForInheritance skips auto');
{
  const sandbox = makeSandbox({
    activeMode: 'default',
    k: 3,
    candidates: {
      'cand_manual_1': {
        source: 'manual',
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 1000000,
        end_bp:   2000000,
      },
      'cand_auto_1': {
        source: 'auto_l2_sweep',
        confirmed: false,
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 3000000,
        end_bp:   4000000,
      },
      'cand_manual_2': {
        source: 'manual',
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 5000000,
        end_bp:   6000000,
      },
      'cand_auto_confirmed': {
        // Confirmed auto candidate — should be included (user accepted it)
        source: 'auto_l2_sweep',
        confirmed: true,
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 7000000,
        end_bp:   8000000,
      },
    },
  });
  vm.runInContext('var items = _gatherActiveCandidatesForInheritance();', sandbox);
  ok('returns 3 items (skipped one unconfirmed auto)',
     sandbox.items.length === 3, 'got: ' + sandbox.items.length);
  const ids = sandbox.items.map(it => it.id).sort();
  ok('manual_1 included',           ids.includes('cand_manual_1'));
  ok('manual_2 included',           ids.includes('cand_manual_2'));
  ok('auto_confirmed included',     ids.includes('cand_auto_confirmed'));
  ok('unconfirmed auto skipped',    !ids.includes('cand_auto_1'));
}

// --- Test 3: defensive fallback works without _isAutoCandidate
console.log('\nTest 3: defensive fallback runs without _isAutoCandidate');
{
  // Build a sandbox containing ONLY the gather function — emulates the
  // turn-129 registry bridge test that sandboxes just this function.
  const sandbox = { state: {
    activeMode: 'default',
    k: 3,
    candidates: {
      'manual_a': {
        source: 'manual',
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 1000000,
        end_bp:   2000000,
      },
      'auto_a': {
        source: 'auto_l2_sweep',
        confirmed: false,
        locked_labels: new Int8Array([0, 1, 2]),
        start_bp: 3000000,
        end_bp:   4000000,
      },
    },
  }, console };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(fnGather, sandbox);
  vm.runInContext('var items = _gatherActiveCandidatesForInheritance();', sandbox);
  // Defensive fallback: even though _isAutoCandidate isn't defined here,
  // the inline source check should still skip the auto candidate.
  ok('defensive path skips auto candidate without _isAutoCandidate',
     sandbox.items.length === 1 && sandbox.items[0].id === 'manual_a');
}

// --- Test 4: gather sequence numbering
console.log('\nTest 4: gather seq_num assignment');
{
  const sandbox = makeSandbox({
    activeMode: 'default',
    k: 3,
    candidates: {
      'cand_b': {
        source: 'manual',
        locked_labels: new Int8Array([0, 1]),
        start_bp: 5000000, end_bp: 6000000,
      },
      'cand_a': {
        source: 'manual',
        locked_labels: new Int8Array([0, 1]),
        start_bp: 1000000, end_bp: 2000000,
      },
    },
  });
  vm.runInContext('var items = _gatherActiveCandidatesForInheritance();', sandbox);
  // Sorted by start_bp ASC; seq_num assigned 1..n
  ok('items sorted by start_bp ASC',
     sandbox.items[0].id === 'cand_a' && sandbox.items[1].id === 'cand_b');
  ok('seq_num 1, 2 assigned in order',
     sandbox.items[0].seq_num === 1 && sandbox.items[1].seq_num === 2);
}

// --- Test 5: sort modifier (sandboxed sort logic)
console.log('\nTest 5: refreshCandidateListUI sort puts auto last');
{
  // Replicate the sort comparator that ships in refreshCandidateListUI
  // and verify the ordering it produces. We can't call the full function
  // (it touches DOM); we just verify the comparator semantics.
  const sandbox = makeSandbox({});
  vm.runInContext(`
    var list = [
      { id: 'old_manual',  source: 'manual',         created_at: 100 },
      { id: 'auto_new',    source: 'auto_l2_sweep',  created_at: 999, confirmed: false },
      { id: 'new_manual',  source: 'manual',         created_at: 200 },
      { id: 'auto_old',    source: 'auto_l2_sweep',  created_at: 1,   confirmed: false },
      { id: 'auto_done',   source: 'auto_l2_sweep',  created_at: 50,  confirmed: true },
    ];
    var sorted = list.slice().sort(function(a, b) {
      var aAuto = _isAutoCandidate(a);
      var bAuto = _isAutoCandidate(b);
      if (aAuto !== bAuto) return aAuto ? 1 : -1;
      return b.created_at - a.created_at;
    });
    var orderedIds = sorted.map(function(x) { return x.id; });
  `, sandbox);
  // Expected order:
  //   non-auto first by created_at desc:
  //     new_manual (200) → old_manual (100) → auto_done (50, but confirmed)
  //   then auto unconfirmed by created_at desc:
  //     auto_new (999) → auto_old (1)
  ok('non-auto candidates first (manual + confirmed auto)',
     sandbox.orderedIds[0] === 'new_manual' &&
     sandbox.orderedIds[1] === 'old_manual' &&
     sandbox.orderedIds[2] === 'auto_done');
  ok('auto unconfirmed last, sorted by created_at desc',
     sandbox.orderedIds[3] === 'auto_new' &&
     sandbox.orderedIds[4] === 'auto_old');
}

// =============================================================================
// Final tally
// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
