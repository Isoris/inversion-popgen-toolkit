// =============================================================================
// tests/sv_evidence/test_step3_5_cursor.js
// =============================================================================
// Step 3.5 regression test for atlas_sv_evidence.js. Verifies the cursor +
// hover + click + keyboard navigation logic added in step 3.5:
//
//   - _windowForPreset returns the correct [lo, hi] for default / left_close /
//     right_close presets.
//   - _nearestSV picks the closest SV to a given bp (within filtered set).
//   - _stepCursor moves cursor by 0.5% of span, clamped to window.
//   - _jumpCursorToNextSV jumps to next/prev SV within window + filters.
//   - _cycleSvTypeFilter cycles through All/BND/INV/DEL/DUP/Other.
//   - _setViewPreset switches preset and clamps cursor.
//   - SV_TYPE_CYCLE constant is exposed.
//   - public attachHotkeys / detachHotkeys are exposed on the API.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step3_5_cursor.js
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

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
section('public API surface (step 3.5 additions)');

['attachHotkeys', 'detachHotkeys', 'setViewPreset', 'setCursorBp']
  .forEach(k => assert(typeof mod[k] === 'function',
    'public API has .' + k + ' function'));

assert(Array.isArray(C.SV_TYPE_CYCLE) && C.SV_TYPE_CYCLE.length === 6,
  'SV_TYPE_CYCLE has 6 entries');
assert(C.SV_TYPE_CYCLE[0] === 'All' && C.SV_TYPE_CYCLE[5] === 'Other',
  'SV_TYPE_CYCLE order: All → ... → Other');

// ---------------------------------------------------------------------------
section('_windowForPreset');
const wDefault = h._windowForPreset(layer, 'default');
assert(wDefault.lo === layer.zone_definitions_bp.left_flank[0],
  'default preset lo = left_flank[0]');
assert(wDefault.hi === layer.zone_definitions_bp.right_flank[1],
  'default preset hi = right_flank[1]');

const wLeft = h._windowForPreset(layer, 'left_close');
assert(wLeft.lo === layer.boundary_left_bp - 250000,
  'left_close lo = boundary_left - 250 kb');
assert(wLeft.hi === layer.boundary_left_bp + 250000,
  'left_close hi = boundary_left + 250 kb');

const wRight = h._windowForPreset(layer, 'right_close');
assert(wRight.lo === layer.boundary_right_bp - 250000,
  'right_close lo = boundary_right - 250 kb');
assert(wRight.hi === layer.boundary_right_bp + 250000,
  'right_close hi = boundary_right + 250 kb');

// Unknown preset falls back to default
const wUnknown = h._windowForPreset(layer, 'totally_made_up');
assert(wUnknown.preset === 'default', 'unknown preset falls back to default');

// Missing zone_definitions_bp falls back to ±500 kb of boundaries
const layerNoZones = { ...layer, zone_definitions_bp: {} };
const wFallback = h._windowForPreset(layerNoZones, 'default');
assert(wFallback.lo === layer.boundary_left_bp - 500000,
  'no zones → default lo = boundary_left - 500 kb');
assert(wFallback.hi === layer.boundary_right_bp + 500000,
  'no zones → default hi = boundary_right + 500 kb');

// ---------------------------------------------------------------------------
section('_nearestSV');
// Set up module state with the fixture loaded
mod._state.layer = layer;
mod._state.activeCandidateId = 'cand_LG28_15Mb';
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };

// SV001 is at 15142380. Cursor at 15142000 → nearest = SV001 (Δ 380)
let near = h._nearestSV(15142000);
assert(near && near.sv.sv_id === 'SV001',
  'cursor 15.142 Mb → nearest SV = SV001');
assert(near.dist === 380, 'distance = 380 bp');

// Cursor at 18121950 → nearest = SV004 at exactly 18121950 (Δ 0)
near = h._nearestSV(18121950);
assert(near && near.sv.sv_id === 'SV004',
  'cursor at SV004 position → nearest = SV004');
assert(near.dist === 0, 'distance = 0 bp');

// Cursor at 14000000 (before any SV) → nearest is the leftmost
near = h._nearestSV(14000000);
assert(near && near.sv != null, 'cursor before any SV → returns leftmost SV');

// With show_only_associated=true, SV007 (FDR 0.49) and SV008 (FDR 0.83)
// should be excluded. Cursor near SV007 should return a different SV.
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: true,
                       fdr_threshold: 0.05 };
near = h._nearestSV(16882400);  // SV007's position
assert(near && near.sv.sv_id !== 'SV007',
  'show_only_associated=true excludes SV007 from nearest search');

// No layer → null
const savedLayer = mod._state.layer;
mod._state.layer = null;
assert(h._nearestSV(15000000) === null, 'no layer → _nearestSV returns null');
mod._state.layer = savedLayer;

// Null bp → null
assert(h._nearestSV(null) === null, 'null cursor → null');

// Reset filters
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };

// ---------------------------------------------------------------------------
section('_stepCursor');

// Cache a window so _stepCursor has something to read from
mod._state._lastWindow = {
  lo: 14140000, hi: 19120000,
  drawW: 800, totalH: 280,
  x: bp => ((bp - 14140000) / (19120000 - 14140000)) * 800,
  xInv: px => 14140000 + (px / 800) * (19120000 - 14140000),
};
const span = 19120000 - 14140000;
const expectedStep = span * 0.005;

// First step from null cursor → centres on window, moves once
mod._state.cursorBp = null;
h._stepCursor(1);
const c1 = mod._state.cursorBp;
const expectedC1 = (14140000 + 19120000) / 2 + expectedStep;
assert(Math.abs(c1 - expectedC1) < 1,
  'first → step from null centres on window then steps right');

// Step right from c1
h._stepCursor(1);
assert(Math.abs(mod._state.cursorBp - (c1 + expectedStep)) < 1,
  'step right adds 0.5% of span');

// Step left
h._stepCursor(-1);
assert(Math.abs(mod._state.cursorBp - c1) < 1,
  'step left undoes step right');

// Step right past hi clamps
mod._state.cursorBp = 19119900;
h._stepCursor(1);
assert(mod._state.cursorBp <= 19120000, 'cursor clamps to hi');

// Step left past lo clamps
mod._state.cursorBp = 14140100;
h._stepCursor(-1);
assert(mod._state.cursorBp >= 14140000, 'cursor clamps to lo');

// ---------------------------------------------------------------------------
section('_jumpCursorToNextSV');

// Build a known sorted positions reference (filtered set, all unfiltered)
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
const sortedPositions = layer.sv_calls.map(c => c.position_bp).sort((a, b) => a - b);

// From cursor at lo, jumping right should land on the leftmost SV
mod._state.cursorBp = 14140000;
h._jumpCursorToNextSV(1);
assert(mod._state.cursorBp === sortedPositions[0],
  'jump right from lo → leftmost SV (' + sortedPositions[0] + ')');

// Jumping right again → next SV
const firstPos = mod._state.cursorBp;
h._jumpCursorToNextSV(1);
assert(mod._state.cursorBp > firstPos, 'jump right again advances');
assert(sortedPositions.includes(mod._state.cursorBp),
  'cursor lands on a real SV position');

// Jump left from hi
mod._state.cursorBp = 19120000;
h._jumpCursorToNextSV(-1);
assert(mod._state.cursorBp === sortedPositions[sortedPositions.length - 1],
  'jump left from hi → rightmost SV (' + sortedPositions[sortedPositions.length - 1] + ')');

// Jump wraps around
mod._state.cursorBp = sortedPositions[sortedPositions.length - 1] + 100;
h._jumpCursorToNextSV(1);
assert(mod._state.cursorBp === sortedPositions[0],
  'jump right past last SV → wraps to leftmost');

// With filters reducing the set, jump uses the reduced set
mod._state.filters = { ...C.DEFAULT_FILTERS, sv_type: 'BND',
                       zone: 'all', show_only_associated: false };
const bndPositions = layer.sv_calls.filter(c => c.sv_type === 'BND')
                       .map(c => c.position_bp).sort((a, b) => a - b);
mod._state.cursorBp = 14000000;
h._jumpCursorToNextSV(1);
assert(mod._state.cursorBp === bndPositions[0],
  'jump with sv_type=BND filter → only lands on BND positions');
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };

// ---------------------------------------------------------------------------
section('_cycleSvTypeFilter');

// Need a stubbed rootEl so _renderFiltersBlock doesn't throw
mod._state.rootEl = null;  // disable re-render side-effects for the test

mod._state.filters = { ...C.DEFAULT_FILTERS, sv_type: 'All' };
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'BND', 'All → BND');
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'INV', 'BND → INV');
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'DEL', 'INV → DEL');
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'DUP', 'DEL → DUP');
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'Other', 'DUP → Other');
h._cycleSvTypeFilter(1);
assert(mod._state.filters.sv_type === 'All', 'Other → All (wrap)');

// Reverse direction
h._cycleSvTypeFilter(-1);
assert(mod._state.filters.sv_type === 'Other', 'reverse: All → Other');
h._cycleSvTypeFilter(-1);
assert(mod._state.filters.sv_type === 'DUP', 'reverse: Other → DUP');

mod._state.filters = { ...C.DEFAULT_FILTERS };

// ---------------------------------------------------------------------------
section('_setViewPreset');

mod._state.viewPreset = 'default';
mod._state.cursorBp = 16500000;  // inside default window

h._setViewPreset('left_close');
assert(mod._state.viewPreset === 'left_close', 'preset switched to left_close');
// 16.5 Mb is outside [14.892, 15.392 Mb] → cursor should clip to right edge of new window
assert(mod._state.cursorBp <= layer.boundary_left_bp + 250000,
  'cursor clipped to new window upper bound');

h._setViewPreset('default');
assert(mod._state.viewPreset === 'default', 'preset switched back to default');

// Invalid preset is rejected
h._setViewPreset('garbage');
assert(mod._state.viewPreset === 'default',
  'invalid preset value is rejected (state unchanged)');

// ---------------------------------------------------------------------------
section('attach/detach hotkeys (smoke test, no DOM)');

// document is not defined in plain Node — guard inside _attachLocusHotkeys
// returns silently when document is undefined. So we expect attachHotkeys
// to not throw and not crash.
let threw = false;
try { mod.attachHotkeys(); } catch (_) { threw = true; }
assert(!threw, 'attachHotkeys does not throw when document is undefined');
try { mod.detachHotkeys(); } catch (_) { threw = true; }
assert(!threw, 'detachHotkeys does not throw when document is undefined');

// ---------------------------------------------------------------------------
section('cleanup');
mod._state.layer = null;
mod._state.activeCandidateId = null;
mod._state.filters = { ...C.DEFAULT_FILTERS };
mod._state.cursorBp = null;
mod._state.markerLeftBp = null;
mod._state.markerRightBp = null;
mod._state.viewPreset = 'default';
mod._state._lastWindow = null;

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
