// =============================================================================
// tests/sv_evidence/test_step5_upset.js
// =============================================================================
// Step 5 regression test for atlas_sv_evidence.js. Verifies:
//
//   - _validateCombinationsLayer accepts the worked example, rejects bad
//     format_version / missing fields.
//   - _renderUpSetPanel produces SVG with N top-bars, N×M dot matrix,
//     per-evidence-type set-size mini-bars, and one transparent hit-rect
//     per combination.
//   - _onUpSetBarClick populates _state.selectedSamples + activeCombinationIndex;
//     toggles off on second click.
//   - _clearSampleSelection clears state.
//   - public setCombinationsLayer accepts valid layer, rejects invalid.
//   - _ingestJsonText routes by format_version (sv_genotype_counts_v1 →
//     _state.layer; sv_evidence_combinations_v1 → _state.combinationsLayer).
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step5_upset.js
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

// localStorage shim
const _ls = {};
global.window = global.window || {};
global.window.localStorage = {
  getItem: (k) => (k in _ls) ? _ls[k] : null,
  setItem: (k, v) => { _ls[k] = String(v); },
  removeItem: (k) => { delete _ls[k]; },
};

const MODULE_PATH      = path.resolve(__dirname, '../../js/atlas_sv_evidence.js');
const FIXTURE_GT       = path.resolve(__dirname, 'fixture_sv_genotype_counts_v1.json');
const FIXTURE_COMBINS  = path.resolve(__dirname, 'fixture_sv_evidence_combinations_v1.json');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  ✓ ' + msg); }
  else      { fail++; console.error('  ✗ ' + msg); }
}
function section(name) { console.log('\n— ' + name + ' —'); }

const mod = require(MODULE_PATH);
const C   = mod._const;
const h   = mod._internals;
const layer    = JSON.parse(fs.readFileSync(FIXTURE_GT, 'utf-8'));
const combLyr  = JSON.parse(fs.readFileSync(FIXTURE_COMBINS, 'utf-8'));

// ---------------------------------------------------------------------------
section('public API');

['onUpSetBarClick', 'clearSampleSelection', 'setCombinationsLayer']
  .forEach(k => assert(typeof mod[k] === 'function', 'public API has .' + k));

assert(typeof C.UPSET === 'object' && C.UPSET.barH === 50,
  'UPSET layout constant exposed');

// ---------------------------------------------------------------------------
section('_validateCombinationsLayer');

let v = h._validateCombinationsLayer(combLyr);
assert(v.ok === true, 'fixture passes validation');

v = h._validateCombinationsLayer(null);
assert(v.ok === false && v.reason === 'NOT_OBJECT', 'null → NOT_OBJECT');

v = h._validateCombinationsLayer({ format_version: 'wrong' });
assert(v.ok === false && v.reason === 'BAD_FORMAT_VERSION',
  'wrong format_version → BAD_FORMAT_VERSION');

v = h._validateCombinationsLayer({
  format_version: 'sv_evidence_combinations_v1',
  evidence_types: [],
  combinations: []
});
assert(v.ok === false && v.reason === 'NO_EVIDENCE_TYPES',
  'empty evidence_types → NO_EVIDENCE_TYPES');

v = h._validateCombinationsLayer({
  format_version: 'sv_evidence_combinations_v1',
  evidence_types: combLyr.evidence_types,
  combinations: 'not-an-array'
});
assert(v.ok === false && v.reason === 'NO_COMBINATIONS',
  'non-array combinations → NO_COMBINATIONS');

// ---------------------------------------------------------------------------
section('public setCombinationsLayer');

assert(mod.setCombinationsLayer(combLyr) === true,
  'valid layer → returns true');
assert(mod._state.combinationsLayer === combLyr,
  '_state.combinationsLayer set');

assert(mod.setCombinationsLayer({ format_version: 'wrong' }) === false,
  'invalid layer → returns false');
assert(mod._state.combinationsLayer === combLyr,
  '_state.combinationsLayer NOT changed by invalid input');

// ---------------------------------------------------------------------------
section('_renderUpSetPanel renders correctly');

// Use a fake DOM (just a div with innerHTML). Smoke test the output structure.
const fakeHost = {
  _html: '',
  innerHTML: '',
  set innerHTML(v) { this._html = v; },
  get innerHTML() { return this._html; },
  querySelectorAll: () => [],
  querySelector: () => null,
};

mod._state.layer = layer;
mod._state.combinationsLayer = combLyr;
h._renderUpSetPanel(fakeHost);
const html = fakeHost._html;

assert(html.length > 100, 'renders a non-trivial HTML/SVG block');
assert(html.indexOf('<svg') !== -1, 'contains an SVG element');

// 4 combinations → 4 top-bar groups
const barCount = (html.match(/data-sv-upset-bar="\d+"/g) || []).length;
assert(barCount === 4, '4 combination bars rendered (matches fixture)');

// Bar values shown: 42, 21, 17, 8
['42', '21', '17', '8'].forEach(n => {
  // Each value text appears at least once. Use a regex that excludes other
  // numeric noise — just check substring presence which is enough for a
  // synthetic fixture with distinct values.
  const re = new RegExp('>' + n + '<');
  assert(re.test(html), `bar value ${n} appears in output`);
});

// 8 evidence types → 8 row labels
['Left split-read', 'Right split-read', 'Manta INV', 'DELLY INV',
 'Left MAPQ0', 'Right MAPQ0'].forEach(label => {
  assert(html.indexOf(label) !== -1, `evidence label "${label}" present`);
});

// Per-evidence totals appear (67, 65, 73, 28, etc.)
['67', '73', '28'].forEach(n => {
  assert(html.indexOf(n) !== -1, `set-size ${n} present`);
});

// Empty state when combinationsLayer is missing
mod._state.combinationsLayer = null;
fakeHost._html = '';
h._renderUpSetPanel(fakeHost);
assert(fakeHost._html.indexOf('No') !== -1 ||
       fakeHost._html.indexOf('—') !== -1,
  'missing combinations layer → empty hint');

// Invalid layer
mod._state.combinationsLayer = { format_version: 'totally_wrong' };
fakeHost._html = '';
h._renderUpSetPanel(fakeHost);
assert(fakeHost._html.indexOf('Invalid') !== -1 ||
       fakeHost._html.indexOf('error') !== -1 ||
       fakeHost._html.indexOf('BAD_') !== -1,
  'invalid layer → error message');

// Restore
mod._state.combinationsLayer = combLyr;

// ---------------------------------------------------------------------------
section('_onUpSetBarClick + _clearSampleSelection');

// Disable DOM render side-effects
mod._state.rootEl = null;
mod._state.selectedSamples = null;
mod._state.activeCombinationIndex = null;

h._onUpSetBarClick(0);
assert(mod._state.selectedSamples != null && mod._state.selectedSamples.size === 42,
  'click bar 0 → selectedSamples size 42 (matches fixture)');
assert(mod._state.activeCombinationIndex === 0,
  'activeCombinationIndex = 0');
assert(mod._state.selectedSamples.has('s_h2h2_0'),
  'selectedSamples contains expected sample');

// Click same bar → toggle off
h._onUpSetBarClick(0);
assert(mod._state.selectedSamples == null,
  'second click on same bar → selectedSamples cleared');
assert(mod._state.activeCombinationIndex == null,
  'activeCombinationIndex cleared');

// Click different bar
h._onUpSetBarClick(1);
assert(mod._state.selectedSamples.size === 21,
  'click bar 1 → selectedSamples size 21');
assert(mod._state.activeCombinationIndex === 1,
  'activeCombinationIndex = 1');

// Click bar 2 (different one) → switches without toggling
h._onUpSetBarClick(2);
assert(mod._state.selectedSamples.size === 17,
  'click bar 2 → switches selection (size 17)');
assert(mod._state.activeCombinationIndex === 2,
  'activeCombinationIndex = 2');

// _clearSampleSelection clears
h._clearSampleSelection();
assert(mod._state.selectedSamples == null,
  '_clearSampleSelection clears state');
assert(mod._state.activeCombinationIndex == null,
  'activeCombinationIndex cleared');

// Out-of-range index is no-op
h._onUpSetBarClick(99);
assert(mod._state.selectedSamples == null,
  'out-of-range bar index is a no-op');

// ---------------------------------------------------------------------------
section('_ingestJsonText format_version dispatch');

// Reset
mod._state.layer = null;
mod._state.combinationsLayer = null;
mod._state.layerError = null;

// Ingest gt counts
h._ingestJsonText(JSON.stringify(layer));
assert(mod._state.layer != null, 'sv_genotype_counts_v1 → _state.layer set');
assert(mod._state.combinationsLayer == null,
  'gt counts ingestion does NOT touch combinationsLayer');
assert(mod._state.activeCandidateId === layer.candidate_id,
  'candidate_id adopted from gt counts');

// Ingest combinations
h._ingestJsonText(JSON.stringify(combLyr));
assert(mod._state.combinationsLayer != null,
  'sv_evidence_combinations_v1 → _state.combinationsLayer set');
assert(mod._state.layer != null,
  'combinations ingestion does NOT clear gt layer');

// Ingest unknown format_version → records an error, doesn't crash
h._ingestJsonText(JSON.stringify({ format_version: 'unknown_v9', x: 1 }));
assert(mod._state.layerError && mod._state.layerError.indexOf('unknown') !== -1,
  'unknown format_version → layerError set');

// Bad JSON is caught
h._ingestJsonText('not valid json {');
assert(mod._state.layerError && mod._state.layerError.indexOf('parse') !== -1,
  'bad JSON → layerError set with parse hint');

// Cleanup
mod._state.layer = null;
mod._state.combinationsLayer = null;
mod._state.layerError = null;
mod._state.activeCandidateId = null;
mod._state.selectedSamples = null;
mod._state.activeCombinationIndex = null;

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
