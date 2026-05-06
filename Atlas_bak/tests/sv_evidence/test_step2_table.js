// =============================================================================
// tests/sv_evidence/test_step2_table.js
// =============================================================================
// Step 2 regression test for atlas_sv_evidence.js. Verifies:
//   - TABLE_COLUMNS registry + TABLE_GROUP_HEADERS exposed.
//   - Column registry covers all spec §3.5 columns (sv_id … notes).
//   - Each column has the required functions (get / fmt / sortKey / tsv).
//   - _applyFilters: SV type / quality / zone / min_samples /
//     show_only_associated / fdr_threshold each behave per spec §3.3.
//   - _sortRows: ascending + descending with null-handling and stable
//     tiebreak on sv_id.
//   - _paginateRows: page index clamping, page-size, startRank/endRank.
//   - _computeVisibleRows composes filter → sort → paginate correctly.
//   - _buildFilteredTSV emits a header + provenance comment block + N rows.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step2_table.js
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
const fix = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf-8'));
const calls = fix.sv_calls;

// ---------------------------------------------------------------------------
section('column registry');

assert(Array.isArray(C.TABLE_COLUMNS) && C.TABLE_COLUMNS.length >= 22,
       'TABLE_COLUMNS has ≥ 22 columns');
assert(Array.isArray(C.TABLE_GROUP_HEADERS) && C.TABLE_GROUP_HEADERS.length === 7,
       'TABLE_GROUP_HEADERS has 7 group spans');

const colIds = C.TABLE_COLUMNS.map(c => c.id);
const expectedIds = [
  'sv_id', 'sv_type', 'zone',
  'position_bp', 'distance_to_edge_bp', 'n_samples',
  'h1h1_AA','h1h1_AB','h1h1_BB','h1h1_miss',
  'h1h2_AA','h1h2_AB','h1h2_BB','h1h2_miss',
  'h2h2_AA','h2h2_AB','h2h2_BB','h2h2_miss',
  'odds_ratio','p_value','fdr',
  'pattern','notes',
];
expectedIds.forEach(id => assert(colIds.indexOf(id) !== -1, 'has column id: ' + id));

C.TABLE_COLUMNS.forEach(col => {
  assert(typeof col.get === 'function',     col.id + '.get is a function');
  assert(typeof col.fmt === 'function',     col.id + '.fmt is a function');
  assert(typeof col.sortKey === 'function', col.id + '.sortKey is a function');
  assert(typeof col.tsv === 'function',     col.id + '.tsv is a function');
});

// Group-header column-count sums to total cols
const totalCols = C.TABLE_GROUP_HEADERS.reduce((a, g) => a + g.cols, 0);
assert(totalCols === C.TABLE_COLUMNS.length,
       'group-header colspans sum to ' + C.TABLE_COLUMNS.length + ' (got ' + totalCols + ')');

// ---------------------------------------------------------------------------
section('column getters / formatters on real fixture row SV001');
const sv001 = calls.find(c => c.sv_id === 'SV001');
const colById = id => C.TABLE_COLUMNS.find(c => c.id === id);

assert(colById('sv_id').get(sv001) === 'SV001',          'sv_id getter');
assert(colById('sv_type').get(sv001) === 'BND',          'sv_type getter');
assert(colById('zone').get(sv001) === 'left_boundary',   'zone getter');
assert(colById('zone').fmt(colById('zone').get(sv001)) === 'Left boundary',
       'zone formatter → "Left boundary"');
assert(colById('position_bp').get(sv001) === 15142380,   'position_bp getter');
assert(colById('position_bp').fmt(15142380) === (15142380).toLocaleString(),
       'position_bp formatter uses toLocaleString');
assert(colById('h1h1_AA').get(sv001) === 55,             'h1h1_AA getter');
assert(colById('h1h2_AB').get(sv001) === 81,             'h1h2_AB getter');
assert(colById('h2h2_BB').get(sv001) === 55,             'h2h2_BB getter');
assert(colById('odds_ratio').get(sv001) === 42.1,        'odds_ratio getter');
assert(colById('odds_ratio').fmt(42.1) === '42.1',       'odds_ratio formatter (10–100)');
assert(colById('odds_ratio').fmt(2.13) === '2.13',       'odds_ratio formatter (<10)');
assert(colById('odds_ratio').fmt(125)  === '125',        'odds_ratio formatter (≥100)');
assert(colById('fdr').get(sv001) === 2.1e-9,             'fdr getter');
assert(colById('pattern').get(sv001) === 'canonical_breakpoint_marker',
       'pattern getter');

// ---------------------------------------------------------------------------
section('_fmtSci edge cases');
assert(h._fmtSci(2.1e-9) === '2.1e-9',  '_fmtSci(2.1e-9)');
assert(h._fmtSci(1.2e-12) === '1.2e-12','_fmtSci(1.2e-12)');
assert(h._fmtSci(0.0041) === '0.0041',  '_fmtSci(0.0041)');
assert(h._fmtSci(0.43)   === '0.43',    '_fmtSci(0.43)');
assert(h._fmtSci(0.49)   === '0.49',    '_fmtSci(0.49)');
assert(h._fmtSci(0.18)   === '0.18',    '_fmtSci(0.18)');
assert(h._fmtSci(1)      === '1.00',    '_fmtSci(1)');
assert(h._fmtSci(0)      === '0',       '_fmtSci(0)');
assert(h._fmtSci(null)   === '—',       '_fmtSci(null)');

// ---------------------------------------------------------------------------
section('_applyFilters — SV type');
let filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, sv_type: 'All',
  show_only_associated: false, zone: 'all' });
assert(filtered.length === calls.length, 'sv_type=All returns all calls');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, sv_type: 'BND',
  show_only_associated: false, zone: 'all' });
assert(filtered.every(c => c.sv_type === 'BND'), 'sv_type=BND returns only BND');
assert(filtered.length === 3, '3 BND calls in fixture');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, sv_type: 'INV',
  show_only_associated: false, zone: 'all' });
assert(filtered.length === 2, '2 INV calls in fixture');

section('_applyFilters — Quality');
filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, quality: 'PASS',
  show_only_associated: false, zone: 'all' });
assert(filtered.length === calls.length, 'all fixture calls are PASS');

section('_applyFilters — Zone');
filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, zone: 'boundary_500kb',
  show_only_associated: false });
assert(filtered.every(c => c.zone === 'left_boundary' || c.zone === 'right_boundary'),
       'boundary_500kb keeps only left/right boundary');
assert(filtered.length === 5, 'fixture has 5 boundary-zone calls');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, zone: 'inversion_body',
  show_only_associated: false });
assert(filtered.every(c => c.zone === 'inversion_body'), 'inversion_body filter');
assert(filtered.length === 2, '2 inversion-body calls');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, zone: 'flanks',
  show_only_associated: false });
assert(filtered.every(c => c.zone === 'left_flank' || c.zone === 'right_flank'),
       'flanks filter');
assert(filtered.length === 1, '1 flank call');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS, zone: 'all',
  show_only_associated: false });
assert(filtered.length === calls.length, 'zone=all returns all');

section('_applyFilters — show_only_associated');
filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS,
  zone: 'all', show_only_associated: true, fdr_threshold: 0.05 });
assert(filtered.every(c => c.fisher && c.fisher.fdr_bh < 0.05),
       'show_only_associated drops FDR ≥ 0.05');
assert(filtered.length === 5, '5 associated calls at FDR < 0.05');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS,
  zone: 'all', show_only_associated: true, fdr_threshold: 1e-6 });
assert(filtered.every(c => c.fisher.fdr_bh < 1e-6), 'fdr<1e-6 filter');
assert(filtered.length === 4, '4 calls at FDR < 1e-6 (SV001/2/4/5)');

section('_applyFilters — min_samples');
filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS,
  zone: 'all', show_only_associated: false, min_samples: 226 });
assert(filtered.length === 2,
  'min_samples=226 keeps only the two calls with n=226 (SV001, SV002)');

filtered = h._applyFilters(calls, { ...C.DEFAULT_FILTERS,
  zone: 'all', show_only_associated: false, min_samples: 1 });
assert(filtered.length === calls.length, 'min_samples=1 keeps all');

// ---------------------------------------------------------------------------
section('_sortRows');
let sorted = h._sortRows(calls, 'fdr', 'asc');
const fdrs = sorted.map(c => c.fisher.fdr_bh);
let isAsc = true;
for (let i = 1; i < fdrs.length; i++) if (fdrs[i] < fdrs[i-1]) { isAsc = false; break; }
assert(isAsc, 'sort by fdr ascending: lowest FDR first');
assert(sorted[0].sv_id === 'SV004', 'first row by fdr asc is SV004 (6.7e-11)');

sorted = h._sortRows(calls, 'fdr', 'desc');
assert(sorted[0].sv_id === 'SV008', 'first row by fdr desc is SV008 (0.83)');

sorted = h._sortRows(calls, 'odds_ratio', 'desc');
assert(sorted[0].sv_id === 'SV004', 'first row by OR desc is SV004 (56.9)');
assert(sorted[1].sv_id === 'SV001', 'second row by OR desc is SV001 (42.1)');

sorted = h._sortRows(calls, 'sv_id', 'asc');
assert(sorted.map(c => c.sv_id).join(',') === 'SV001,SV002,SV003,SV004,SV005,SV006,SV007,SV008',
       'sort by sv_id asc');

// Null-handling: distance_to_edge_bp is null on SV006 / SV007
sorted = h._sortRows(calls, 'distance_to_edge_bp', 'asc');
assert(sorted[sorted.length - 1].distance_to_edge_bp == null ||
       sorted[sorted.length - 2].distance_to_edge_bp == null,
       'null distances sort to the end (asc)');

// ---------------------------------------------------------------------------
section('_paginateRows');
const rows = Array.from({ length: 25 }, (_, i) => ({ sv_id: 'X' + i }));

let p = h._paginateRows(rows, 10, 0);
assert(p.rows.length === 10,                'page 0 of 25 size 10 → 10 rows');
assert(p.startRank === 1 && p.endRank === 10, 'startRank=1 endRank=10');
assert(p.lastPage === 2,                    'lastPage = 2 (3 pages 0,1,2)');

p = h._paginateRows(rows, 10, 2);
assert(p.rows.length === 5,                 'page 2 of 25 size 10 → 5 rows');
assert(p.startRank === 21 && p.endRank === 25, 'startRank=21 endRank=25');

p = h._paginateRows(rows, 10, 999);
assert(p.pageIndex === 2,                   'page index clamped to lastPage');
assert(p.rows.length === 5,                 'clamped page yields 5 rows');

p = h._paginateRows([], 10, 0);
assert(p.rows.length === 0 && p.total === 0 && p.startRank === 0 && p.endRank === 0,
       'empty input → empty paged output');

// ---------------------------------------------------------------------------
section('_computeVisibleRows (compose)');
const v = h._computeVisibleRows(fix,
  { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false },
  'fdr', 'asc', 25, 0);
assert(v.total === calls.length,        'compose: total = layer.sv_calls.length');
assert(v.rows.length === calls.length,  'compose: page contains all rows (one page)');
assert(v.rows[0].sv_id === 'SV004',     'compose: top row is SV004 (lowest FDR)');
assert(Array.isArray(v.filteredAll),    'compose: filteredAll exposed for TSV');
assert(v.filteredAll.length === calls.length, 'compose: filteredAll matches total');

// With show_only_associated default ON (DEFAULT_FILTERS), only 5 rows pass
const vAssoc = h._computeVisibleRows(fix,
  { ...C.DEFAULT_FILTERS, zone: 'all' /* keep show_only_associated default */ },
  'fdr', 'asc', 25, 0);
assert(vAssoc.total === 5, 'compose with show_only_associated=true → 5 rows');

// ---------------------------------------------------------------------------
section('_buildFilteredTSV');

// Stub the global atlas state to hand the karyotype counts to the TSV builder
global.window = global.window || {};
global.window.state = {
  candidateState: {
    cand_LG28_15Mb: {
      locked_labels: (function () {
        const out = [];
        for (let i = 0; i < 61; i++)  out.push({ sample_id: 'a' + i, label: 'HOMO_1' });
        for (let i = 0; i < 103; i++) out.push({ sample_id: 'b' + i, label: 'HET'    });
        for (let i = 0; i < 62; i++)  out.push({ sample_id: 'c' + i, label: 'HOMO_2' });
        return out;
      }()),
    },
  },
};
mod._state.layer = fix;
mod._state.activeCandidateId = 'cand_LG28_15Mb';
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
mod._state.sortColumnId = 'fdr';
mod._state.sortDirection = 'asc';

const tsv = h._buildFilteredTSV();
assert(typeof tsv === 'string' && tsv.length > 0, 'TSV output is non-empty string');

const tsvLines = tsv.split('\n');
const commentLines = tsvLines.filter(l => l.startsWith('#'));
assert(commentLines.length >= 8, 'TSV has ≥ 8 provenance comment lines');
assert(commentLines.some(l => l.includes('candidate_id: cand_LG28_15Mb')),
       'TSV provenance includes candidate_id');
assert(commentLines.some(l => l.includes('chrom: C_gar_LG28')),
       'TSV provenance includes chrom');
assert(commentLines.some(l => l.includes('n_rows: ' + calls.length)),
       'TSV provenance includes n_rows');
assert(commentLines.some(l => l.includes('sort: fdr asc')),
       'TSV provenance includes sort directive');

// First non-comment line is the header
const firstData = tsvLines.findIndex(l => !l.startsWith('#') && l.length > 0);
const headerCells = tsvLines[firstData].split('\t');
// Step 2 had 23 columns. Step 4 appended a "User annotation" column at the
// end (per Quentin's row-annotation feature), so the count is ≥ 23.
assert(headerCells.length >= C.TABLE_COLUMNS.length,
       'TSV header has ≥ ' + C.TABLE_COLUMNS.length + ' columns (got ' + headerCells.length + ')');
assert(headerCells.includes('sv_id'),  'TSV header has sv_id');
assert(headerCells.some(h => h.startsWith('H1/H1 (n=61)')),
       'TSV header annotates H1/H1 group with n=61');
assert(headerCells.some(h => h.startsWith('H1/H2 (n=103)')),
       'TSV header annotates H1/H2 group with n=103');
assert(headerCells.some(h => h.startsWith('H2/H2 (n=62)')),
       'TSV header annotates H2/H2 group with n=62');

// Body rows
const bodyRows = tsvLines.slice(firstData + 1).filter(l => l.length > 0);
assert(bodyRows.length === calls.length, 'TSV body has ' + calls.length + ' rows');
const firstRow = bodyRows[0].split('\t');
assert(firstRow[0] === 'SV004', 'first TSV row is SV004 (lowest FDR)');
assert(firstRow.length === headerCells.length,
       'first TSV row column count matches header (' + headerCells.length + ')');

// Verify the filtered-only TSV (with show_only_associated=true)
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all' /* show_only_assoc=true default */ };
const tsv2 = h._buildFilteredTSV();
const lines2 = tsv2.split('\n');
const body2 = lines2.slice(lines2.findIndex(l => !l.startsWith('#') && l.length > 0) + 1)
                    .filter(l => l.length > 0);
assert(body2.length === 5, 'filtered TSV (show_only_associated) has 5 rows');

// Cleanup
mod._state.layer = null;
mod._state.activeCandidateId = null;
mod._state.filters = { ...C.DEFAULT_FILTERS };
delete global.window.state;

// ---------------------------------------------------------------------------
section('exportFilteredTSV public API');
const empty = mod.exportFilteredTSV();
assert(typeof empty === 'string', 'exportFilteredTSV returns string when no layer');
// The empty path should not throw and should still return a header-shape line
const eHead = empty.split('\n').find(l => !l.startsWith('#') && l.length > 0);
assert(eHead && eHead.split('\t').length >= 5,
       'empty exportFilteredTSV header has ≥ 5 columns');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
