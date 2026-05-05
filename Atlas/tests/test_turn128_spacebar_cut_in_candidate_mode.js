// =============================================================================
// turn 128 integration test — spacebar contextual binding
//
// Quentin's request: "I would like to have spacebar to cut the interval when
// in candidate mode we can cut L2 intervals to make more L2. Its easier to
// merge."
//
// Implementation (minimum viable interpretation):
//   When state.candidateMode is ON and toggleL3DraftCutAtCursor places/removes
//   a cut, spacebar consumes the event. Otherwise (no candidate mode, or no
//   valid draft to cut), spacebar falls through to togglePlay() — the global
//   play/pause role is preserved.
//
// The bigger feature ("literally split L2 envelopes into two new L2s on the
// registry") is deferred — it's a major architectural mutation that would
// affect catalogue persistence, K-cluster caching, page-4 sample regimes,
// page-8 windows table, and Save Session round-trip. See the spec stub
// added under specs_todo/.
//
// Tests are source-level: verify that the keyboard handler now calls
// toggleL3DraftCutAtCursor for spacebar when candidateMode is on, and falls
// back to togglePlay otherwise. We also assert that the C-key path is
// untouched (still calls toggleL3DraftCutAtCursor unconditionally — i.e.,
// the function itself gates on candidateMode + draft state).
// =============================================================================

const fs = require('fs');
const path = require('path');

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

// The spacebar branch must check state.candidateMode AND call
// toggleL3DraftCutAtCursor. After turn 128c, there's no togglePlay fallback —
// when not in candidate mode (or no valid cut), the handler quietly returns
// without preventDefault, so the browser's native spacebar behaviour stands.
ok('spacebar branch checks state.candidateMode',
   /if \(e\.key === ' '\) \{[\s\S]{0,500}state\.candidateMode/.test(html),
   'spacebar handler should gate on candidateMode');

ok('spacebar branch calls toggleL3DraftCutAtCursor when gated true',
   /if \(e\.key === ' '\) \{[\s\S]{0,500}toggleL3DraftCutAtCursor\(\)/.test(html));

ok('spacebar consumed via e.preventDefault when cut placed',
   /if \(toggleL3DraftCutAtCursor\(\)\) \{ e\.preventDefault\(\); return; \}/.test(html));

ok('spacebar wraps cut call in try/catch for fail-soft fallback',
   /if \(state\.candidateMode && typeof toggleL3DraftCutAtCursor === 'function'\) \{\s*try \{/.test(html));

// turn 128c: spacebar no longer falls through to togglePlay. Quentin asked
// to unbind the global play/pause shortcut: "You can unbind play button
// space bar since it's not useful really." So the spacebar branch must NOT
// reach togglePlay anywhere in its body. The Play button itself (#playBtn)
// keeps its mouse-click handler — only the global key shortcut is gone.
//
// We sniff the spacebar block (between `if (e.key === ' ') {` and its
// matching `}`) and assert togglePlay is NOT inside it.
function _extractSpacebarBlock(src) {
  const m = src.match(/if \(e\.key === ' '\) \{/);
  if (!m) return null;
  let i = m.index + m[0].length;
  let depth = 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const q = ch; i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === q) break;
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
  return src.substring(m.index, i);
}

const spacebarBlock = _extractSpacebarBlock(html);
ok('spacebar block extractable', !!spacebarBlock);
if (spacebarBlock) {
  ok('spacebar block does NOT call togglePlay anywhere',
     spacebarBlock.indexOf('togglePlay') < 0,
     'turn 128c: global play/pause shortcut unbound');
  // The block should still preventDefault the cut path — but NOT preventDefault
  // when cut didn't happen (so the browser can do native spacebar = page scroll).
  ok('spacebar block has at least one preventDefault (the cut-success path)',
     /e\.preventDefault\(\)/.test(spacebarBlock));
  // The block should NOT do an UNCONDITIONAL preventDefault — i.e., the LAST
  // statement in the block should be a bare `return;`, not preventDefault.
  ok('spacebar block ends with bare return (no unconditional preventDefault)',
     /\/\* no global play\/pause anymore[\s\S]{0,200}return;\s*\}\s*$/.test(spacebarBlock));
}

// Old spacebar binding (unconditional togglePlay only) must be gone.
const oldOneLiner = /if \(e\.key === ' '\)\s*\{\s*e\.preventDefault\(\); togglePlay\(\); \}\s*else if \(e\.key === 'ArrowRight'/;
ok('old single-purpose spacebar handler is gone',
   !oldOneLiner.test(html),
   'spacebar should now route through candidate-mode cut shortcut only');

// C key handler is untouched — still calls toggleL3DraftCutAtCursor directly
ok('C-key handler still wired to toggleL3DraftCutAtCursor',
   /if \(e\.key === 'c' \|\| e\.key === 'C'\) \{\s*if \(typeof toggleL3DraftCutAtCursor === 'function' && toggleL3DraftCutAtCursor\(\)\)/.test(html));

// =============================================================================
// Logic test: simulate the handler's branching with stubs
// =============================================================================
console.log('\n=== Branching logic ===');

function simulateSpacebarBranch(opts) {
  const state = { candidateMode: opts.candidateMode };
  const calls = [];
  function toggleL3DraftCutAtCursor() {
    calls.push('cut');
    return opts.cutReturns;     // true if it placed/removed a cut
  }
  function togglePlay() {
    calls.push('play');          // we now expect this to be NEVER called
  }
  // Inline reproduction of the handler logic for the spacebar branch
  // (matching the turn 128c version — no togglePlay fallback).
  let preventedDefault = false;
  function preventDefault() { preventedDefault = true; }
  if (true) {   // simulating "key === ' '"
    if (state.candidateMode && typeof toggleL3DraftCutAtCursor === 'function') {
      try {
        if (toggleL3DraftCutAtCursor()) { preventDefault(); return { calls, preventedDefault, returned: true }; }
      } catch (_) {}
    }
    /* no global play/pause anymore — return without preventDefault */
    return { calls, preventedDefault, returned: false };
  }
  return { calls, preventedDefault, returned: false };
}

// Case 1: candidate mode ON, cut succeeded → cut called, preventDefault, exit early
{
  const r = simulateSpacebarBranch({ candidateMode: true, cutReturns: true });
  ok('candidate mode ON + cut succeeds → only "cut" called',
     r.calls.length === 1 && r.calls[0] === 'cut');
  ok('candidate mode ON + cut succeeds → preventDefault called',
     r.preventedDefault === true);
  ok('candidate mode ON + cut succeeds → handler exits early',
     r.returned === true);
  ok('candidate mode ON + cut succeeds → togglePlay NOT called',
     r.calls.indexOf('play') < 0);
}

// Case 2: candidate mode ON but cut returned false (cursor at edge etc.)
//         → "cut" attempted, no "play" since togglePlay is unbound, NO preventDefault
//           so the browser does its native spacebar behaviour (page scroll).
{
  const r = simulateSpacebarBranch({ candidateMode: true, cutReturns: false });
  ok('candidate mode ON + cut returns false → only "cut" attempted',
     r.calls.length === 1 && r.calls[0] === 'cut');
  ok('candidate mode ON + cut returns false → togglePlay NOT called',
     r.calls.indexOf('play') < 0);
  ok('candidate mode ON + cut returns false → preventDefault NOT called (browser handles space natively)',
     r.preventedDefault === false);
}

// Case 3: candidate mode OFF — spacebar does NOTHING in our handler.
//         Browser's native spacebar handling stands.
{
  const r = simulateSpacebarBranch({ candidateMode: false, cutReturns: true });
  ok('candidate mode OFF → "cut" never called',
     r.calls.indexOf('cut') < 0);
  ok('candidate mode OFF → "play" never called (turn 128c unbind)',
     r.calls.indexOf('play') < 0);
  ok('candidate mode OFF → preventDefault NOT called (browser handles space)',
     r.preventedDefault === false);
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
