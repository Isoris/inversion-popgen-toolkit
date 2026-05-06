// =============================================================================
// turn 147b — Windows-(N) button engagement state
// =============================================================================
// Quentin (turn 144 screenshot 2): "📊 Windows (N) button at top-right shows
// no engagement state when N≠1".
//
// Pre-fix behavior:
//   - state.stepMode='l2'    → label "📊 Windows (1)", dim, MISLEADING (no
//                              cycler step is active; arrows jump L2)
//   - state.stepMode='winN'  → label "📊 Windows (1)", dim, WRONG (the cycler
//                              isn't live, but the user IS stepping N
//                              windows — the button hides this)
//
// Post-fix behavior:
//   - state.stepMode='l2'    → label "📊 Windows (—)", DIM (correctly
//                              indicates "I'm not the live step setting")
//   - state.stepMode='win1'  → label "📊 Windows (1)", LIT (cycler live)
//   - state.stepMode='win5'  → label "📊 Windows (5)", LIT (cycler live)
//   - state.stepMode='win10' → label "📊 Windows (10)", LIT
//   - state.stepMode='winN'  → label "📊 Windows (N)" with actual N from
//                              state.stepModeN, LIT (so user knows arrows
//                              ARE stepping in window units)
//
// Plus: typing a new N value into the sidebar number input now refreshes
// the cycler label immediately (was silent before).
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
// 1. Source-level changes
// ============================================================================
console.log('\n=== 1. Source-level changes ===');
ok('_refreshStepSizeBtn updated with l2 case',
   /m === 'l2' \|\| m == null/.test(html));
ok('l2 mode emits "Windows \\(—\\)" label',
   /label\s*=\s*'📊 Windows \(—\)'/.test(html));
ok('winN mode reads state.stepModeN',
   /m === 'winN'[\s\S]{0,100}?state\.stepModeN/.test(html));
ok('winN mode emits the actual N in the label',
   /m === 'winN'[\s\S]{0,250}?label\s*=\s*`📊 Windows \(\$\{n\}\)`/.test(html));
ok('winN mode lights up (lit=true)',
   /m === 'winN'[\s\S]{0,300}?lit\s*=\s*true/.test(html));
ok('button title is set with explanatory tooltip',
   /btn\.title\s*=\s*title/.test(html));
ok('stepModeNInput change triggers cycler refresh',
   /_stepModeNInput[\s\S]{0,800}?turn 147[\s\S]{0,200}?_refreshStepSizeBtn/.test(html));

// ============================================================================
// 2. Sandboxed _refreshStepSizeBtn behavior
// ============================================================================
console.log('\n=== 2. Sandboxed cycler-refresh behavior ===');
{
  // Extract _stepSizeForMode + _refreshStepSizeBtn together with the
  // _STEP_CYCLE_ORDER const.
  const helperRe = /const _STEP_CYCLE_KEY[\s\S]*?function _refreshStepSizeBtn\(\)\s*\{[\s\S]*?\n\}/;
  const m = html.match(helperRe);
  if (!m) { ok('extract refresh-cycler helpers', false, 'not found'); }
  else {
    // Build a tiny mock DOM
    function makeBtn() {
      return {
        textContent: '',
        title: '',
        _classes: new Set(),
        classList: {
          toggle(name, force) {
            if (force === true)  this._owner._classes.add(name);
            else if (force === false) this._owner._classes.delete(name);
            else if (this._owner._classes.has(name)) this._owner._classes.delete(name);
            else this._owner._classes.add(name);
          },
          contains(name) { return this._owner._classes.has(name); },
        },
        _setOwner() { this.classList._owner = this; return this; },
      };
    }

    const btn = makeBtn()._setOwner();
    const sandbox = {
      state: { stepMode: 'l2', stepModeN: 1 },
      document: { getElementById: (id) => (id === 'jumpToWindowsBtn') ? btn : null },
      console,
    };
    const ctx = vm.createContext(sandbox);
    vm.runInContext(m[0], ctx);

    // Test l2 mode
    sandbox.state.stepMode = 'l2';
    ctx._refreshStepSizeBtn();
    ok('stepMode=l2 → label "📊 Windows (—)"',  btn.textContent === '📊 Windows (—)');
    ok('stepMode=l2 → cycler dim (no is-active-step)',
       !btn._classes.has('is-active-step'));
    ok('stepMode=l2 → tooltip mentions L2 envelopes',
       /L2 envelopes/.test(btn.title));

    // Test win1
    sandbox.state.stepMode = 'win1';
    ctx._refreshStepSizeBtn();
    ok('stepMode=win1 → label "📊 Windows (1)"',  btn.textContent === '📊 Windows (1)');
    ok('stepMode=win1 → cycler lit',              btn._classes.has('is-active-step'));

    // Test win5
    sandbox.state.stepMode = 'win5';
    ctx._refreshStepSizeBtn();
    ok('stepMode=win5 → label "📊 Windows (5)"',  btn.textContent === '📊 Windows (5)');
    ok('stepMode=win5 → cycler lit',              btn._classes.has('is-active-step'));

    // Test win10
    sandbox.state.stepMode = 'win10';
    ctx._refreshStepSizeBtn();
    ok('stepMode=win10 → label "📊 Windows (10)"', btn.textContent === '📊 Windows (10)');
    ok('stepMode=win10 → cycler lit',              btn._classes.has('is-active-step'));

    // Test winN with custom value 15 (the SCREENSHOT 2 case)
    sandbox.state.stepMode  = 'winN';
    sandbox.state.stepModeN = 15;
    ctx._refreshStepSizeBtn();
    ok('stepMode=winN, N=15 → label "📊 Windows (15)"',
       btn.textContent === '📊 Windows (15)');
    ok('stepMode=winN → cycler LIT (the fix!)',
       btn._classes.has('is-active-step'));
    ok('stepMode=winN → tooltip mentions custom value',
       /custom value/.test(btn.title));

    // Test winN with default 1
    sandbox.state.stepMode  = 'winN';
    sandbox.state.stepModeN = 1;
    ctx._refreshStepSizeBtn();
    ok('stepMode=winN, N=1 → label "📊 Windows (1)"',
       btn.textContent === '📊 Windows (1)');
    ok('stepMode=winN, N=1 → cycler still lit',
       btn._classes.has('is-active-step'));

    // Test winN with very large value
    sandbox.state.stepMode  = 'winN';
    sandbox.state.stepModeN = 200;
    ctx._refreshStepSizeBtn();
    ok('stepMode=winN, N=200 → label "📊 Windows (200)"',
       btn.textContent === '📊 Windows (200)');

    // Test stepMode=null defensively
    sandbox.state.stepMode = null;
    ctx._refreshStepSizeBtn();
    ok('stepMode=null → label "📊 Windows (—)" (defensive)',
       btn.textContent === '📊 Windows (—)');

    // Test winN with bad stepModeN value (defensive)
    sandbox.state.stepMode  = 'winN';
    sandbox.state.stepModeN = -5;
    ctx._refreshStepSizeBtn();
    ok('stepMode=winN, bad N → falls back to N=1',
       btn.textContent === '📊 Windows (1)');
  }
}

// ============================================================================
// 3. Existing CSS rule still matches the cycler-lit state
// ============================================================================
console.log('\n=== 3. CSS rule still wired ===');
ok('CSS rule for #jumpToWindowsBtn.is-active-step still present',
   /#viewModeBar\s+#jumpToWindowsBtn\.is-active-step\s*\{[\s\S]{0,200}?color:\s*var\(--accent\)/.test(html));
ok('CSS rule applies accent border + accent text',
   /#viewModeBar\s+#jumpToWindowsBtn\.is-active-step[\s\S]{0,200}?border-color:\s*var\(--accent\)/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
