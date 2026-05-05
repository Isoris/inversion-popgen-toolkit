// =============================================================================
// tests/sv_evidence/test_step4_rightrail.js
// =============================================================================
// Step 4 regression test for atlas_sv_evidence.js. Verifies:
//
//   - _renderBoundarySummaryHTML produces 2 stacked tables (left + right)
//     with one row per SV type + a totals row, and click-targets that
//     carry both side and type.
//   - _renderLegendHTML produces 2 columns (FDR + pattern labels).
//   - _readAnnotations / _writeAnnotations / _cycleAnnotation: round-trip,
//     null → L → R → ★ → null cycle, persisted to localStorage.
//   - _stepHighlightedRow: navigates through filtered+sorted rows in the
//     correct order, wraps, paginates if needed.
//   - TSV export now includes a 'User annotation' column at the end.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step4_rightrail.js
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

// In-memory localStorage shim (Node has no window.localStorage)
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
section('boundary summary HTML');

const bsHtml = h._renderBoundarySummaryHTML(layer);
assert(typeof bsHtml === 'string', 'returns string');

// Two blocks, one per side
const leftMatches  = (bsHtml.match(/sv-bsblock-left\b/g)  || []).length;
const rightMatches = (bsHtml.match(/sv-bsblock-right\b/g) || []).length;
assert(leftMatches  === 1, 'one Left boundary block');
assert(rightMatches === 1, 'one Right boundary block');

// Each block has 5 SV-type rows + 1 totals row → 6 <tr>s per block, 12 total
const trMatches = (bsHtml.match(/<tr/g) || []).length;
assert(trMatches >= 14, 'has at least 14 <tr> (header + 5 type + totals × 2 sides)');

// Each row carries side + type data attributes for click handler
['BND','INV','DEL','DUP','Other'].forEach(t => {
  assert(new RegExp(`data-sv-bs-side="left"\\s+data-sv-bs-type="${t}"`).test(bsHtml),
    'left.' + t + ' row has data attrs');
  assert(new RegExp(`data-sv-bs-side="right"\\s+data-sv-bs-type="${t}"`).test(bsHtml),
    'right.' + t + ' row has data attrs');
});

// Total row reflects fixture: left has 24+11+38+22+7 = 102 total, 9+5+7+4+0 = 25 assoc
assert(/<tr class="sv-bs-total">[\s\S]*?>Total<[\s\S]*?>102</.test(bsHtml),
  'left total = 102 SVs');
assert(/<tr class="sv-bs-total">[\s\S]*?>Total<[\s\S]*?>25</.test(bsHtml),
  'left total assoc = 25');

// Empty / partial: missing boundary_summary section
const layerNoBs = { ...layer, boundary_summary: { left: null, right: null } };
const bsHtmlEmpty = h._renderBoundarySummaryHTML(layerNoBs);
assert(bsHtmlEmpty.indexOf('no summary in layer') !== -1,
  'missing summary → "no summary in layer" hint per side');

// ---------------------------------------------------------------------------
section('legend HTML');

const legHtml = h._renderLegendHTML();
assert(legHtml.indexOf('Association (FDR)') !== -1, 'has FDR column title');
assert(legHtml.indexOf('Pattern labels')    !== -1, 'has Pattern labels column title');

// 5 FDR tiers
const fdrItems = (legHtml.match(/data-sv-legend-fdr=/g) || []).length;
assert(fdrItems === 5, '5 FDR swatches (one per tier)');

// 6 pattern labels
const patternItems = (legHtml.match(/data-sv-legend-pattern=/g) || []).length;
assert(patternItems === 6, '6 pattern-label swatches');

// Each pattern colour appears
Object.values(C.PATTERN_LABELS).forEach(pl => {
  assert(legHtml.indexOf(pl.color) !== -1, 'pattern colour ' + pl.color + ' present');
});

// ---------------------------------------------------------------------------
section('row annotations');

// Initial state — no annotations stored
mod._state.activeCandidateId = 'cand_LG28_15Mb';
mod._state.rowAnnotations = h._readAnnotations('cand_LG28_15Mb');
assert(typeof mod._state.rowAnnotations === 'object' && mod._state.rowAnnotations !== null,
  'fresh annotations is an object');
assert(Object.keys(mod._state.rowAnnotations).length === 0,
  'no annotations on first read');

// Cycle: null → L → R → star → null
h._cycleAnnotation('SV001');
assert(mod._state.rowAnnotations.SV001 === 'L', 'cycle 1: null → L');
h._cycleAnnotation('SV001');
assert(mod._state.rowAnnotations.SV001 === 'R', 'cycle 2: L → R');
h._cycleAnnotation('SV001');
assert(mod._state.rowAnnotations.SV001 === 'star', 'cycle 3: R → star');
h._cycleAnnotation('SV001');
assert(!('SV001' in mod._state.rowAnnotations), 'cycle 4: star → null (key removed)');

// Mix two SVs
h._cycleAnnotation('SV004');                     // L
h._cycleAnnotation('SV005'); h._cycleAnnotation('SV005');  // R
assert(mod._state.rowAnnotations.SV004 === 'L', 'SV004 → L independent');
assert(mod._state.rowAnnotations.SV005 === 'R', 'SV005 → R independent');

// Persist + reload
const reloaded = h._readAnnotations('cand_LG28_15Mb');
assert(reloaded.SV004 === 'L' && reloaded.SV005 === 'R',
  'annotations persisted to localStorage and round-trip');

// Different candidate → different namespace
const otherCid = h._readAnnotations('cand_LG02_99Mb');
assert(Object.keys(otherCid).length === 0,
  'different candidate has its own annotation set');

// _renderAnnoChip output
assert(h._renderAnnoChip('L').indexOf('sv-anno-L') !== -1,
  'chip for "L" has sv-anno-L class');
assert(h._renderAnnoChip('R').indexOf('sv-anno-R') !== -1,
  'chip for "R" has sv-anno-R class');
assert(h._renderAnnoChip('star').indexOf('sv-anno-star') !== -1,
  'chip for "star" has sv-anno-star class');
assert(h._renderAnnoChip(null).indexOf('sv-anno-empty') !== -1,
  'chip for null has sv-anno-empty class');

// Cleanup
mod._state.rowAnnotations = {};
mod._state.activeCandidateId = null;

// ---------------------------------------------------------------------------
section('_stepHighlightedRow');

// Set up state with the fixture layer
mod._state.layer = layer;
mod._state.activeCandidateId = 'cand_LG28_15Mb';
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
mod._state.sortColumnId = 'fdr';
mod._state.sortDirection = 'asc';
mod._state.pageSize = 25;
mod._state.pageIndex = 0;
mod._state.highlightedSvId = null;

// We can't call _stepHighlightedRow directly because it triggers DOM
// rendering. Replicate its core decision logic here to verify the helper's
// effect on _state.highlightedSvId without needing a DOM. Use the public
// internal directly — the renderer guards on rootEl being null.
mod._state.rootEl = null; // disable render side effects

// Sorted-by-FDR-asc order should put SV004 first
h._stepHighlightedRow(1);
assert(mod._state.highlightedSvId === 'SV004',
  'first ↓ press → highlight SV004 (lowest FDR)');

h._stepHighlightedRow(1);
assert(mod._state.highlightedSvId === 'SV001',
  'second ↓ press → SV001 (next-lowest FDR, 2.1e-9)');

h._stepHighlightedRow(-1);
assert(mod._state.highlightedSvId === 'SV004',
  '↑ press → back to SV004');

// Wrap from first → last
h._stepHighlightedRow(-1);
// SV008 has highest FDR (0.83) so when sorted ascending, it's last
assert(mod._state.highlightedSvId === 'SV008',
  '↑ from first row wraps to last (SV008, highest FDR)');

// Now test pagination: set pageSize=2, ensure ↓ flips pages
mod._state.pageSize = 2;
mod._state.pageIndex = 0;
mod._state.highlightedSvId = 'SV004';  // page 0 (rank 0)
h._stepHighlightedRow(1);  // SV001 (rank 1, still page 0)
assert(mod._state.pageIndex === 0, 'page index stays 0 for rank 1 with pageSize 2');
h._stepHighlightedRow(1);  // rank 2 → page 1
assert(mod._state.pageIndex === 1, 'page index flips to 1 when crossing page boundary');

// Cleanup
mod._state.pageSize = 25;
mod._state.pageIndex = 0;
mod._state.highlightedSvId = null;
mod._state.layer = null;
mod._state.activeCandidateId = null;
mod._state.filters = { ...C.DEFAULT_FILTERS };

// ---------------------------------------------------------------------------
section('TSV export now includes annotation column');

mod._state.layer = layer;
mod._state.activeCandidateId = 'cand_LG28_15Mb';
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
mod._state.sortColumnId = 'fdr';
mod._state.sortDirection = 'asc';
mod._state.rowAnnotations = { SV001: 'L', SV004: 'star' };

const tsv = h._buildFilteredTSV();
const lines = tsv.split('\n');
const headerLine = lines.find(l => !l.startsWith('#') && l.length > 0);
assert(headerLine.split('\t').slice(-1)[0] === 'User annotation',
  'TSV header last column is "User annotation"');

// Find SV001 row → annotation column should read "L"
const sv001Line = lines.find(l => l.startsWith('SV001\t'));
assert(sv001Line, 'SV001 row exists in TSV');
if (sv001Line) {
  const cells = sv001Line.split('\t');
  assert(cells[cells.length - 1] === 'L',
    'SV001 row annotation column is "L"');
}

// Find SV004 row → annotation column should read "star"
const sv004Line = lines.find(l => l.startsWith('SV004\t'));
if (sv004Line) {
  const cells = sv004Line.split('\t');
  assert(cells[cells.length - 1] === 'star',
    'SV004 row annotation column is "star"');
}

// Find a row with no annotation → empty cell
const sv006Line = lines.find(l => l.startsWith('SV006\t'));
if (sv006Line) {
  const cells = sv006Line.split('\t');
  assert(cells[cells.length - 1] === '',
    'unannotated row (SV006) has empty annotation cell');
}

// Provenance now includes n_annotations
const annoMetaLine = lines.find(l => l.startsWith('# n_annotations:'));
assert(annoMetaLine && annoMetaLine.includes('2'),
  'TSV provenance includes n_annotations: 2');

// Cleanup
mod._state.rowAnnotations = {};
mod._state.layer = null;
mod._state.activeCandidateId = null;
mod._state.filters = { ...C.DEFAULT_FILTERS };

// ---------------------------------------------------------------------------
section('public API additions');
// All step-4 internals exposed
['_readAnnotations','_writeAnnotations','_cycleAnnotation',
 '_renderAnnoChip','_renderBoundarySummaryHTML','_renderLegendHTML',
 '_stepHighlightedRow'].forEach(k => {
  assert(typeof h[k] === 'function', '_internals.' + k + ' is exposed');
});

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
