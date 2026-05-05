// =============================================================================
// tests/sv_evidence/test_step4_5_select.js
// =============================================================================
// Step 4.5 regression test for atlas_sv_evidence.js. Verifies:
//
//   - public API additions: setSelectMode, setSelection
//   - _setSelectMode toggles state.selectMode + clears in-progress drag
//   - state.selection filter integrates into _applyFilters (only SVs whose
//     position_bp falls inside [startBp, endBp] survive)
//   - _commitSelectDrag commits a drag span ≥ 1 kb as a sticky selection
//     and treats shorter drags as a click (no selection, just cursor)
//   - selection co-exists with other filters (composes correctly)
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step4_5_select.js
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

// In-memory localStorage shim
const _ls = {};
global.window = global.window || {};
global.window.localStorage = {
  getItem: (k) => (k in _ls) ? _ls[k] : null,
  setItem: (k, v) => { _ls[k] = String(v); },
  removeItem: (k) => { delete _ls[k]; },
};

const MODULE_PATH  = path.resolve(__dirname, '../../js/atlas_sv_evidence.js');
const FIXTURE_PATH = path.resolve(__dirname, 'fixture_sv_genotype_counts_v1.json');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  ✓ ' + msg); }
  else      { fail++; console.error('  ✗ ' + msg); }
}
function section(name) { console.log('\n— ' + name + ' —'); }

const mod = require(MODULE_PATH);
const C   = mod._const;
const h   = mod._internals;
const layer = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf-8'));

// ---------------------------------------------------------------------------
section('public API');

['setSelectMode', 'setSelection'].forEach(k => {
  assert(typeof mod[k] === 'function', 'public API has .' + k);
});

// ---------------------------------------------------------------------------
section('_setSelectMode');

// Reset state cleanly
mod._state.selectMode = 'zoom';
mod._state.selection = null;
mod._state._selectDrag = null;
mod._state.rootEl = null;  // disable DOM side-effects

h._setSelectMode('select');
assert(mod._state.selectMode === 'select', 'selectMode toggled to "select"');

h._setSelectMode('zoom');
assert(mod._state.selectMode === 'zoom', 'selectMode back to "zoom"');

// Invalid value rejected
h._setSelectMode('garbage');
assert(mod._state.selectMode === 'zoom', 'invalid mode rejected (state unchanged)');

// In-progress drag is cleared on mode switch
mod._state._selectDrag = { startBp: 100, currentBp: 200 };
h._setSelectMode('select');
assert(mod._state._selectDrag == null, 'mode switch clears in-progress drag');

// Sticky selection is preserved on mode switch
mod._state.selection = { startBp: 1000, endBp: 2000 };
h._setSelectMode('zoom');
assert(mod._state.selection != null, 'mode switch preserves sticky selection');

// Cleanup
mod._state.selection = null;
mod._state._selectDrag = null;
mod._state.selectMode = 'zoom';

// ---------------------------------------------------------------------------
section('selection filter integrates into _applyFilters');

const calls = layer.sv_calls;
const noFilters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };

// Without any selection, all calls survive
mod._state.selection = null;
let filtered = h._applyFilters(calls, noFilters);
assert(filtered.length === calls.length,
  'no selection → all ' + calls.length + ' calls pass');

// Selection covering whole range — same result
mod._state.selection = { startBp: 0, endBp: 1e10 };
filtered = h._applyFilters(calls, noFilters);
assert(filtered.length === calls.length,
  'selection [0, 1e10] → all calls pass');

// Tight selection around SV001 (15.142380 Mb) and SV002 (15.138900 Mb)
mod._state.selection = { startBp: 15100000, endBp: 15200000 };
filtered = h._applyFilters(calls, noFilters);
const ids = filtered.map(c => c.sv_id).sort();
assert(ids.includes('SV001') && ids.includes('SV002'),
  'tight selection 15.10–15.20 Mb keeps SV001 + SV002');
assert(!ids.includes('SV004') && !ids.includes('SV005'),
  'tight selection excludes right-boundary SV004 + SV005');

// Selection at right boundary
mod._state.selection = { startBp: 18100000, endBp: 18200000 };
filtered = h._applyFilters(calls, noFilters);
const idsRight = filtered.map(c => c.sv_id).sort();
assert(idsRight.includes('SV004') && idsRight.includes('SV005'),
  'right-boundary selection keeps SV004 + SV005');
assert(!idsRight.includes('SV001'),
  'right-boundary selection excludes SV001');

// Selection composes with other filters: tight zone + sv_type=BND
mod._state.selection = { startBp: 15100000, endBp: 15200000 };
filtered = h._applyFilters(calls, { ...noFilters, sv_type: 'BND' });
const idsCompose = filtered.map(c => c.sv_id);
assert(idsCompose.includes('SV001'),
  'compose: selection ∩ sv_type=BND keeps SV001 (BND)');
assert(!idsCompose.includes('SV002'),
  'compose: selection ∩ sv_type=BND drops SV002 (INV)');

// Selection that excludes everything
mod._state.selection = { startBp: 99999000, endBp: 99999999 };
filtered = h._applyFilters(calls, noFilters);
assert(filtered.length === 0,
  'selection outside all SV positions → 0 calls');

// Cleanup
mod._state.selection = null;

// ---------------------------------------------------------------------------
section('public setSelection wrapper');

mod.setSelection(15100000, 15200000);
assert(mod._state.selection != null,
  'setSelection populates _state.selection');
assert(mod._state.selection.startBp === 15100000 &&
       mod._state.selection.endBp === 15200000,
  'setSelection records start + end');

// setSelection auto-orders endpoints
mod.setSelection(15200000, 15100000);
assert(mod._state.selection.startBp === 15100000 &&
       mod._state.selection.endBp === 15200000,
  'setSelection auto-orders reversed args');

// setSelection(null, null) clears
mod.setSelection(null, null);
assert(mod._state.selection == null, 'setSelection(null, null) clears');

// ---------------------------------------------------------------------------
section('_commitSelectDrag — short drag = click, long drag = selection');

// Short drag (< 1 kb) — clears any selection, sets cursor at midpoint
mod._state.selection = { startBp: 5000, endBp: 9000 }; // pre-existing
mod._state._selectDrag = { startBp: 14500000, currentBp: 14500500 }; // 500 bp drag
h._commitSelectDrag();
assert(mod._state.selection == null,
  'short drag (<1 kb) clears any pre-existing selection');
assert(Math.abs(mod._state.cursorBp - 14500250) < 1,
  'short drag places cursor at midpoint of drag');

// Long drag (≥ 1 kb) — commits as selection
mod._state._selectDrag = { startBp: 14500000, currentBp: 14600000 }; // 100 kb
h._commitSelectDrag();
assert(mod._state.selection != null,
  'long drag (≥1 kb) commits as selection');
assert(mod._state.selection.startBp === 14500000 &&
       mod._state.selection.endBp === 14600000,
  'commit records correct start + end');
assert(mod._state._selectDrag == null,
  'commit clears in-progress drag');

// Long drag with reversed direction — start > current
mod._state.selection = null;
mod._state._selectDrag = { startBp: 16000000, currentBp: 15000000 }; // dragged left
h._commitSelectDrag();
assert(mod._state.selection.startBp === 15000000 &&
       mod._state.selection.endBp === 16000000,
  'reversed drag direction → selection low/high orderered correctly');

// _commitSelectDrag with no in-flight drag is a no-op
mod._state._selectDrag = null;
mod._state.selection = { startBp: 1, endBp: 2 };
h._commitSelectDrag();
assert(mod._state.selection.startBp === 1,
  'commit with no drag in flight is a no-op');

// Cleanup
mod._state.selection = null;
mod._state._selectDrag = null;
mod._state.cursorBp = null;

// ---------------------------------------------------------------------------
section('selection survives non-applicable filter changes');

mod._state.selection = { startBp: 14640000, endBp: 15640000 }; // left-boundary zone
filtered = h._applyFilters(calls, { ...noFilters, zone: 'boundary_500kb' });
const ids2 = filtered.map(c => c.sv_id);
assert(ids2.includes('SV001') && ids2.includes('SV002'),
  'selection + zone=boundary_500kb: SV001 + SV002 survive (left boundary)');
assert(!ids2.includes('SV007'),
  'selection + zone=boundary_500kb: SV007 (body) filtered out');

mod._state.selection = null;

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
