// =============================================================================
// tests/sv_evidence/test_step6_heatmap.js
// =============================================================================
// Step 6 regression test for atlas_sv_evidence.js (Sample × SV dosage heatmap).
// Verifies:
//
//   - _validateSupportLayer accepts the worked example, rejects bad
//     format_version, missing samples / sv_ids, dosage row-count mismatch.
//   - _supportCell returns the right char for known cells; null for OOB.
//   - _renderHeatmapPanel produces SVG with the expected backdrop rects,
//     overlay cells for HET / HOM-ALT, group separator lines, group labels,
//     SV row labels, and the inline legend.
//   - _renderHeatmapPanel routes to compact vs large depending on
//     _state.heatmapView.
//   - _renderHeatmapPanel emits an empty-state hint when supportLayer is null.
//   - _renderHeatmapPanel emits an "Invalid support layer" error block when
//     the layer is structurally bad.
//   - _hasNoCarrierInSelection correctly identifies SVs with zero carriers
//     in the active sample selection (the locus-glyph dim hook).
//   - _ingestJsonText routes sv_support_by_sample_v1 → _state.supportLayer
//     without disturbing _state.layer or _state.combinationsLayer.
//   - public setSupportLayer accepts valid, rejects invalid.
//   - public setHeatmapView toggles _state.heatmapView and removes overlay
//     when collapsing.
//   - HEATMAP layout constants (compact + large + colours) exposed.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step6_heatmap.js
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

const MODULE_PATH    = path.resolve(__dirname, '../../js/atlas_sv_evidence.js');
const FIXTURE_GT     = path.resolve(__dirname, 'fixture_sv_genotype_counts_v1.json');
const FIXTURE_COMB   = path.resolve(__dirname, 'fixture_sv_evidence_combinations_v1.json');
const FIXTURE_SUPP   = path.resolve(__dirname, 'fixture_sv_support_by_sample_v1.json');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  \u2713 ' + msg); }
  else      { fail++; console.error('  \u2717 ' + msg); }
}
function section(name) { console.log('\n\u2014 ' + name + ' \u2014'); }

const mod   = require(MODULE_PATH);
const C     = mod._const;
const h     = mod._internals;
const layer = JSON.parse(fs.readFileSync(FIXTURE_GT,   'utf-8'));
const cmLyr = JSON.parse(fs.readFileSync(FIXTURE_COMB, 'utf-8'));
const supp  = JSON.parse(fs.readFileSync(FIXTURE_SUPP, 'utf-8'));

// Helper to deep-clone a JSON-safe object so we can mutate without leaking.
const clone = o => JSON.parse(JSON.stringify(o));

// ---------------------------------------------------------------------------
section('public API surface');

['setSupportLayer', 'setHeatmapView']
  .forEach(k => assert(typeof mod[k] === 'function', 'public API has .' + k));

// ---------------------------------------------------------------------------
section('HEATMAP layout constants');

assert(typeof C.HEATMAP === 'object', 'HEATMAP constant exposed');
assert(C.HEATMAP.compact && typeof C.HEATMAP.compact.cellW === 'number',
       'HEATMAP.compact.cellW is a number');
assert(C.HEATMAP.compact.cellH > 0 && C.HEATMAP.compact.cellH < 100,
       'HEATMAP.compact.cellH in plausible range');
assert(C.HEATMAP.large && typeof C.HEATMAP.large.cellW === 'number',
       'HEATMAP.large.cellW is a number');
assert(C.HEATMAP.large.cellH > C.HEATMAP.compact.cellH,
       'large cells taller than compact cells');
assert(C.HEATMAP.colours && typeof C.HEATMAP.colours['0'] === 'string',
       'HEATMAP.colours has AA entry');
assert(C.HEATMAP.colours['1'] && C.HEATMAP.colours['2'],
       'HEATMAP.colours has HET and HOM-ALT entries');
assert(C.HEATMAP.colours['.'] === 'transparent',
       'HEATMAP.colours: missing maps to transparent');
assert(typeof C.HEATMAP.compact.maxRows === 'number' &&
       C.HEATMAP.compact.maxRows >= 4,
       'HEATMAP.compact.maxRows is a sensible cap');

// ---------------------------------------------------------------------------
section('_validateSupportLayer');

let v = h._validateSupportLayer(supp);
assert(v.ok === true, 'fixture passes validation');

v = h._validateSupportLayer(null);
assert(v.ok === false && v.reason === 'NOT_OBJECT',
       'null \u2192 NOT_OBJECT');

v = h._validateSupportLayer(undefined);
assert(v.ok === false && v.reason === 'NOT_OBJECT',
       'undefined \u2192 NOT_OBJECT');

v = h._validateSupportLayer({ format_version: 'wrong_v9' });
assert(v.ok === false && v.reason === 'BAD_FORMAT_VERSION',
       'wrong format_version \u2192 BAD_FORMAT_VERSION');

v = h._validateSupportLayer({
  format_version: 'sv_support_by_sample_v1',
  samples: [],
  sv_ids: ['SV001'],
  dosage_compact: [],
});
assert(v.ok === false && v.reason === 'NO_SAMPLES',
       'empty samples \u2192 NO_SAMPLES');

v = h._validateSupportLayer({
  format_version: 'sv_support_by_sample_v1',
  samples: ['s1'],
  sv_ids: [],
  dosage_compact: ['0'],
});
assert(v.ok === false && v.reason === 'NO_SV_IDS',
       'empty sv_ids \u2192 NO_SV_IDS');

v = h._validateSupportLayer({
  format_version: 'sv_support_by_sample_v1',
  samples: ['s1', 's2'],
  sv_ids: ['SV001'],
  dosage_compact: ['0'],   // length mismatch
});
assert(v.ok === false && v.reason === 'DOSAGE_ROW_COUNT_MISMATCH',
       'dosage rows != samples \u2192 DOSAGE_ROW_COUNT_MISMATCH');

v = h._validateSupportLayer({
  format_version: 'sv_support_by_sample_v1',
  samples: ['s1'],
  sv_ids: ['SV001'],
  dosage_compact: 'not-an-array',
});
assert(v.ok === false && v.reason === 'DOSAGE_ROW_COUNT_MISMATCH',
       'non-array dosage_compact \u2192 DOSAGE_ROW_COUNT_MISMATCH');

// ---------------------------------------------------------------------------
section('_supportCell');

// Fixture row 0 is "00...0..." (8 chars, AA / miss / miss / miss / AA / miss / miss / miss)
const ch0_0 = h._supportCell(supp, 0, 0);
assert(ch0_0 === '0', 'row 0, col 0 \u2192 \'0\' (AA)');
const ch0_1 = h._supportCell(supp, 0, 1);
assert(ch0_1 === '0', 'row 0, col 1 \u2192 \'0\' (AA)');
const ch0_2 = h._supportCell(supp, 0, 2);
assert(ch0_2 === '.', 'row 0, col 2 \u2192 \'.\' (miss)');

// HET block — first row of het is "10100100"
const hetStartIdx = supp.row_groups['H1/H2'][0];
const chHet0 = h._supportCell(supp, hetStartIdx, 0);
assert(chHet0 === '1', 'first HET row, col 0 \u2192 \'1\' (HET dosage)');
const chHet2 = h._supportCell(supp, hetStartIdx, 2);
assert(chHet2 === '1', 'first HET row, col 2 \u2192 \'1\'');

// H2/H2 first row is "22222200"
const h22StartIdx = supp.row_groups['H2/H2'][0];
const chH22_0 = h._supportCell(supp, h22StartIdx, 0);
assert(chH22_0 === '2', 'first H2/H2 row, col 0 \u2192 \'2\' (HOM-ALT)');
const chH22_6 = h._supportCell(supp, h22StartIdx, 6);
assert(chH22_6 === '0', 'first H2/H2 row, col 6 \u2192 \'0\' (no SV here)');

// OOB
assert(h._supportCell(supp, -1, 0) == null || h._supportCell(supp, -1, 0) === '',
       'negative rowIdx \u2192 nullish');
assert(h._supportCell(supp, 99999, 0) == null,
       'OOB rowIdx \u2192 null');
const chOob = h._supportCell(supp, 0, 99);
assert(chOob == null || chOob === '',
       'OOB colIdx \u2192 null/empty');

// Null layer
assert(h._supportCell(null, 0, 0) == null,
       'null layer \u2192 null');

// ---------------------------------------------------------------------------
section('_renderHeatmapPanel — empty / error / valid (compact)');

const fakeHost = {
  _html: '',
  innerHTML: '',
  set innerHTML(v) { this._html = v; },
  get innerHTML() { return this._html; },
  querySelectorAll: () => [],
  querySelector: () => null,
  classList: { contains: () => false },
};

// Reset state and try empty
mod._state.layer        = null;
mod._state.supportLayer = null;
mod._state.heatmapView  = 'compact';

h._renderHeatmapPanel(fakeHost);
assert(fakeHost._html.length > 0, 'renders something for empty state');
assert(fakeHost._html.indexOf('No') !== -1 || fakeHost._html.indexOf('not loaded') !== -1
       || fakeHost._html.indexOf('Drop') !== -1,
       'empty state mentions missing / drop folder');

// Invalid layer
mod._state.supportLayer = { format_version: 'totally_wrong' };
fakeHost._html = '';
h._renderHeatmapPanel(fakeHost);
assert(fakeHost._html.indexOf('Invalid') !== -1 ||
       fakeHost._html.indexOf('BAD_') !== -1,
       'invalid layer \u2192 error block');

// Valid layer
mod._state.layer        = layer;
mod._state.supportLayer = supp;
fakeHost._html = '';
h._renderHeatmapPanel(fakeHost);
const html = fakeHost._html;
assert(html.length > 200, 'renders a non-trivial HTML/SVG block (compact)');
assert(html.indexOf('<svg') !== -1, 'contains an SVG element');
assert(html.indexOf('data-sv-hm-svg') !== -1, 'SVG carries data-sv-hm-svg attribute');

// AA backdrop: one rect per visible row spanning the full sample width
const backdropMatches = html.match(/data-sv-row="\d+"/g) || [];
assert(backdropMatches.length >= 1,
       'at least one AA backdrop rect rendered (got ' + backdropMatches.length + ')');

// Group separators — 2 vertical lines between H1/H1 ↔ H1/H2 ↔ H2/H2
const sepCount = (html.match(/<line[^>]*stroke-width="0\.7"/g) || []).length;
assert(sepCount >= 2,
       'at least 2 group-separator lines (got ' + sepCount + ')');

// Group labels: H1/H1 (n=61), H1/H2 (n=103), H2/H2 (n=62) — substring search
// (slashes are not HTML-encoded; the plain string lives in the SVG text node).
['H1/H1 (n=61)', 'H1/H2 (n=103)', 'H2/H2 (n=62)'].forEach(lbl => {
  assert(html.indexOf(lbl) !== -1,
    'group label "' + lbl + '" appears in output');
});

// SV row labels — under the default filter (boundary_500kb + show_only_associated)
// only the 5 boundary-zoned associated SVs render; SV006/007/008 are in the
// inversion body / flanks and are filtered out. That filtering is part of the
// designed contract: the heatmap stays in sync with the table.
['SV001', 'SV002', 'SV003'].forEach(svid => {
  assert(html.indexOf(svid) !== -1, 'SV row label "' + svid + '" present');
});
const labelCount = (html.match(/data-sv-hm-svid="SV00\d"/g) || []).length;
assert(labelCount >= 4 && labelCount <= 8,
  'between 4 and 8 SV row labels render under default filters (got ' + labelCount + ')');

// data-sv-hm-svid attribute on SV row labels (for click handlers)
const svidLabels = (html.match(/data-sv-hm-svid="SV00\d"/g) || []).length;
assert(svidLabels >= 4, 'SV row labels carry data-sv-hm-svid (got ' + svidLabels + ')');

// With every filter relaxed, all 8 SVs in the support fixture should render
// (n_svs = 8, well below compact maxRows = 12)
const savedFilters = mod._state.filters;
mod._state.filters = {
  sv_type: 'All', quality: 'All', zone: 'All',
  min_samples: 0, show_only_associated: false, fdr_threshold: 1.0
};
const allHtml = h._renderHeatmapHTML(supp, 'compact');
const allLabels = (allHtml.match(/data-sv-hm-svid="SV00\d"/g) || []).length;
assert(allLabels === 8,
  'with relaxed filters, all 8 SVs render in the heatmap (got ' + allLabels + ')');
mod._state.filters = savedFilters;

// Inline legend
assert(html.indexOf('AA') !== -1 && html.indexOf('AB') !== -1 &&
       html.indexOf('BB') !== -1, 'inline legend present (AA / AB / BB)');

// Expand button (compact mode)
assert(html.indexOf('data-sv-hm-action="expand"') !== -1,
       'compact mode shows expand button');

// HET / HOM-ALT overlay rects — at least one of each colour fills should appear
const hetCol = C.HEATMAP.colours['1'];
const homCol = C.HEATMAP.colours['2'];
assert(html.indexOf('fill="' + hetCol + '"') !== -1,
       'HET-coloured rect rendered');
assert(html.indexOf('fill="' + homCol + '"') !== -1,
       'HOM-ALT coloured rect rendered');

// ---------------------------------------------------------------------------
section('_renderHeatmapPanel — large mode placeholder');

mod._state.heatmapView = 'large';
fakeHost._html = '';
h._renderHeatmapPanel(fakeHost);
assert(fakeHost._html.indexOf('main area') !== -1 ||
       fakeHost._html.indexOf('close it') !== -1,
       'compact host shows a placeholder when large is open');

// Reset
mod._state.heatmapView = 'compact';

// ---------------------------------------------------------------------------
section('_renderHeatmapHTML — truncation cap');

// Synthesize a layer with 50 SVs — must truncate to compact maxRows (12)
const bigSupp = clone(supp);
bigSupp.sv_ids = [];
for (let i = 0; i < 50; i++) bigSupp.sv_ids.push('SVB' + i.toString().padStart(3, '0'));
// Pad each row out to 50 chars (all AA so the test is simple)
bigSupp.dosage_compact = bigSupp.dosage_compact.map(r => '0'.repeat(50));

mod._state.layer = null;        // bypass filter; just rely on the cap
mod._state.supportLayer = bigSupp;
const truncHtml = h._renderHeatmapHTML(bigSupp, 'compact');
assert(truncHtml.indexOf('showing top') !== -1,
       'truncation note shown when n_svs > maxRows.compact');
assert(truncHtml.indexOf('SVB000') !== -1,
       'first SV rendered when truncated');
assert(truncHtml.indexOf('SVB049') === -1,
       'last SV NOT rendered when truncated (cap kicks in)');

// Restore
mod._state.layer        = layer;
mod._state.supportLayer = supp;

// ---------------------------------------------------------------------------
section('_renderHeatmapHTML — filter integration');

// Apply a sv_type filter that only matches 1-2 SVs in the fixture
mod._state.filters = { ...C.DEFAULT_FILTERS, sv_type: 'BND' };
const filtHtml = h._renderHeatmapHTML(supp, 'compact');
// At least one SV should be rendered (the BND ones); not all 8
const svidsInFiltered = (filtHtml.match(/data-sv-hm-svid="SV00\d"/g) || []).length;
assert(svidsInFiltered >= 1,
       'filtered heatmap still shows >= 1 row (BNDs from fixture)');
assert(svidsInFiltered < supp.sv_ids.length,
       'filtered heatmap shows fewer rows than full fixture');

// Reset filters
mod._state.filters = { ...C.DEFAULT_FILTERS };

// ---------------------------------------------------------------------------
section('_hasNoCarrierInSelection');

// No selection \u2192 returns false (don't dim anything)
mod._state.selectedSamples = null;
assert(h._hasNoCarrierInSelection('SV001') === false,
       'no selection \u2192 _hasNoCarrierInSelection returns false');

// No support layer \u2192 returns false
mod._state.selectedSamples = new Set(['s_h1h1_0']);
const tmpSupp = mod._state.supportLayer;
mod._state.supportLayer = null;
assert(h._hasNoCarrierInSelection('SV001') === false,
       'no support layer \u2192 returns false');
mod._state.supportLayer = tmpSupp;

// Selection contains only H1/H1 samples (rows 0-60), all AA for SV001 \u2192 dim
mod._state.selectedSamples = new Set(supp.samples.slice(0, 5));
assert(h._hasNoCarrierInSelection('SV001') === true,
       'selection of H1/H1 samples (all AA) for SV001 \u2192 should dim (no carriers)');

// Selection contains some HET samples (rows 61-163), HET for SV001 \u2192 don't dim
mod._state.selectedSamples = new Set(supp.samples.slice(70, 80));
assert(h._hasNoCarrierInSelection('SV001') === false,
       'selection containing HET carriers \u2192 should NOT dim');

// Unknown SV id \u2192 returns false (can't determine)
mod._state.selectedSamples = new Set(supp.samples.slice(70, 80));
assert(h._hasNoCarrierInSelection('SV_UNKNOWN') === false,
       'SV not in support matrix \u2192 returns false');

// Cleanup
mod._state.selectedSamples = null;

// ---------------------------------------------------------------------------
section('_ingestJsonText routes sv_support_by_sample_v1');

mod._state.layer            = null;
mod._state.combinationsLayer = null;
mod._state.supportLayer     = null;
mod._state.layerError       = null;

// Ingest gt counts first
h._ingestJsonText(JSON.stringify(layer));
assert(mod._state.layer != null, 'gt counts \u2192 _state.layer set');
assert(mod._state.supportLayer == null,
       'gt counts ingestion does NOT touch supportLayer');

// Ingest support layer
h._ingestJsonText(JSON.stringify(supp));
assert(mod._state.supportLayer != null,
       'support layer \u2192 _state.supportLayer set');
assert(mod._state.layer != null,
       'support ingestion does NOT clear gt layer');
assert(mod._state.combinationsLayer == null,
       'support ingestion does NOT touch combinationsLayer');

// Ingest combinations to be sure all three coexist
h._ingestJsonText(JSON.stringify(cmLyr));
assert(mod._state.combinationsLayer != null,
       'combinations layer \u2192 _state.combinationsLayer set');
assert(mod._state.supportLayer != null,
       'combinations ingestion does NOT clear supportLayer');
assert(mod._state.layer != null,
       'combinations ingestion does NOT clear layer');

// Bad support layer (right format_version, missing samples) \u2192 layerError
mod._state.supportLayer = null;
mod._state.layerError = null;
h._ingestJsonText(JSON.stringify({
  format_version: 'sv_support_by_sample_v1',
  samples: [], sv_ids: ['x'], dosage_compact: [],
}));
assert(mod._state.layerError != null && mod._state.layerError !== '',
       'invalid support layer \u2192 layerError set');
assert(mod._state.supportLayer == null,
       'invalid support layer \u2192 supportLayer NOT set');

// ---------------------------------------------------------------------------
section('public setSupportLayer');

mod._state.supportLayer = null;
mod._state.layerError   = null;
mod._state.rootEl       = null;     // disable DOM render

assert(mod.setSupportLayer(supp) === true,
       'valid layer \u2192 setSupportLayer returns true');
assert(mod._state.supportLayer === supp,
       '_state.supportLayer set after setSupportLayer');

assert(mod.setSupportLayer({ format_version: 'wrong' }) === false,
       'invalid layer \u2192 setSupportLayer returns false');
assert(mod._state.supportLayer === supp,
       '_state.supportLayer NOT changed by invalid input');

// ---------------------------------------------------------------------------
section('public setHeatmapView');

mod._state.heatmapView = 'compact';
mod.setHeatmapView('large');
assert(mod._state.heatmapView === 'large',
       'setHeatmapView("large") flips state');

mod.setHeatmapView('compact');
assert(mod._state.heatmapView === 'compact',
       'setHeatmapView("compact") flips back');

// Invalid value is no-op
mod.setHeatmapView('bogus');
assert(mod._state.heatmapView === 'compact',
       'invalid view value is a no-op');

// ---------------------------------------------------------------------------
section('_recomputeGtCountsView — selection-scoped recount (deferred-feature unlock)');

// No selection or no support \u2192 view is null
mod._state.layer            = layer;
mod._state.supportLayer     = supp;
mod._state.selectedSamples  = null;
mod._state.gtCountsView     = null;
let view = h._recomputeGtCountsView();
assert(view === null, 'no selection \u2192 view is null');
assert(mod._state.gtCountsView === null, 'state.gtCountsView is null');

mod._state.selectedSamples = new Set(['s_h1h1_0']);
mod._state.supportLayer    = null;
view = h._recomputeGtCountsView();
assert(view === null, 'no support layer \u2192 view is null');

// Restore support; pick a selection that's well-defined.
mod._state.supportLayer = supp;
// Select 5 H1/H1 samples (all-AA rows) + 5 H2/H2 samples (mostly-BB rows).
const selH11 = supp.samples.slice(0, 5);
const selH22 = supp.samples.slice(supp.row_groups['H2/H2'][0],
                                  supp.row_groups['H2/H2'][0] + 5);
mod._state.selectedSamples = new Set([...selH11, ...selH22]);
view = h._recomputeGtCountsView();
assert(view !== null, 'with selection + support \u2192 view computed');
assert(typeof view === 'object', 'view is an object');

// Per-group N counts within the selection
assert(view._groupNs['H1/H1'] === 5, 'H1/H1 in-selection count = 5');
assert(view._groupNs['H1/H2'] === 0, 'H1/H2 in-selection count = 0');
assert(view._groupNs['H2/H2'] === 5, 'H2/H2 in-selection count = 5');

// SV001 in fixture: H1/H1 rows are all "0" \u2192 5 AA. H2/H2 rows are "2..." \u2192 5 BB.
const v001 = view['SV001'];
assert(v001 != null, 'view has SV001 entry');
assert(v001['H1/H1'].AA === 5,  'SV001 H1/H1.AA = 5 (all-AA)');
assert(v001['H1/H1'].AB === 0,  'SV001 H1/H1.AB = 0');
assert(v001['H1/H1'].BB === 0,  'SV001 H1/H1.BB = 0');
assert(v001['H1/H1'].miss === 0, 'SV001 H1/H1.miss = 0');
assert(v001['H1/H2'].AA === 0,  'SV001 H1/H2.AA = 0 (no H1/H2 in sel)');
assert(v001['H2/H2'].BB === 5,  'SV001 H2/H2.BB = 5 (all-HOM-ALT block)');
assert(v001['H2/H2'].AA === 0,  'SV001 H2/H2.AA = 0');

// Now check the table accessor (gc) reads the view, not call.genotype_counts
const gtTotal = layer.sv_calls.find(c => c.sv_id === 'SV001').genotype_counts['H1/H1'].AA;
assert(gtTotal !== 5,
       'sanity: full-cohort SV001 H1/H1 AA != 5 (different denominator)');
const sv001Call = layer.sv_calls.find(c => c.sv_id === 'SV001');
const aaCol = mod._const.TABLE_COLUMNS.find(col => col.id === 'h1h1_AA');
assert(typeof aaCol === 'object', 'table has h1h1_AA column');
const cellVal = aaCol.get(sv001Call);
assert(cellVal === 5,
       `gc() reads from gtCountsView when active (got ${cellVal}, expected 5)`);

// Clear selection \u2192 falls back to full-cohort counts
mod._state.gtCountsView = null;
mod._state.selectedSamples = null;
const fallbackVal = aaCol.get(sv001Call);
assert(fallbackVal === gtTotal,
       'gc() falls back to call.genotype_counts when view is null');

// SV not in support matrix \u2192 gc() falls back transparently for that SV
mod._state.gtCountsView = view;        // restore active view
const fakeCall = { sv_id: 'SV_GHOST',
                   genotype_counts: { 'H1/H1': { AA: 99, AB: 0, BB: 0, miss: 0 },
                                      'H1/H2': {}, 'H2/H2': {} } };
const ghostVal = aaCol.get(fakeCall);
assert(ghostVal === 99,
       'SV not in view \u2192 falls back to call.genotype_counts (got ' + ghostVal + ')');

// ---------------------------------------------------------------------------
section('_onUpSetBarClick triggers recount when supportLayer is loaded');

// Set up: load all 3 layers
mod._state.layer            = layer;
mod._state.combinationsLayer = cmLyr;
mod._state.supportLayer     = supp;
mod._state.selectedSamples  = null;
mod._state.activeCombinationIndex = null;
mod._state.gtCountsView     = null;
mod._state.rootEl           = null;     // disable DOM render

// Click bar 0 (top combination) \u2192 should populate gtCountsView
h._onUpSetBarClick(0);
assert(mod._state.selectedSamples != null && mod._state.selectedSamples.size > 0,
       'UpSet bar 0 click \u2192 selectedSamples populated');
assert(mod._state.gtCountsView != null,
       'UpSet bar 0 click \u2192 gtCountsView computed');
assert(typeof mod._state.gtCountsView._groupNs === 'object',
       'gtCountsView has _groupNs');

// Click again (toggle off) \u2192 view cleared
h._onUpSetBarClick(0);
assert(mod._state.selectedSamples == null,
       'second click \u2192 selectedSamples cleared');
assert(mod._state.gtCountsView == null,
       'second click \u2192 gtCountsView cleared');

// Click bar without support layer \u2192 view stays null
mod._state.supportLayer = null;
mod._state.gtCountsView = null;
h._onUpSetBarClick(0);
assert(mod._state.selectedSamples != null,
       'click without support \u2192 selectedSamples still set');
assert(mod._state.gtCountsView == null,
       'click without support \u2192 gtCountsView stays null (graceful no-op)');

// Restore for fallthrough
mod._state.supportLayer = supp;
h._clearSampleSelection();
assert(mod._state.gtCountsView == null,
       '_clearSampleSelection clears gtCountsView');

// Set selection BEFORE support layer (race condition); then load support
mod._state.combinationsLayer = cmLyr;
mod._state.supportLayer = null;
mod._state.gtCountsView = null;
h._onUpSetBarClick(2);
assert(mod._state.selectedSamples != null && mod._state.gtCountsView == null,
       'selection set without support \u2192 view null');
// Now set support via public API \u2192 should backfill the recount
mod.setSupportLayer(supp);
assert(mod._state.gtCountsView != null,
       'setSupportLayer with active selection \u2192 backfills gtCountsView');

// Cleanup
h._clearSampleSelection();
mod._state.supportLayer     = null;
mod._state.combinationsLayer = null;

// ---------------------------------------------------------------------------
section('Cleanup');

mod._state.layer            = null;
mod._state.combinationsLayer = null;
mod._state.supportLayer     = null;
mod._state.layerError       = null;
mod._state.activeCandidateId = null;
mod._state.selectedSamples  = null;
mod._state.activeCombinationIndex = null;
mod._state.heatmapView      = 'compact';
mod._state.heatmapHoverCell = null;
assert(true, 'state reset');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
