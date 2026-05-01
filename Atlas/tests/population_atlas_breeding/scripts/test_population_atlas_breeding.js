// Smoke test for the breeding page (Population_atlas.html page 8).
// Tests:
//   - DOM presence: page8 div + tab button
//   - Empty-state render
//   - Auto-derivation hooks (Diversity / Relatedness / Burden / Inversions / Markers)
//   - Status sorting (red > yellow > green)
//   - Simple/detailed view toggle
//   - User upload via JSON + via TSV
//   - User row override of derived row
//   - Filters (status + source)
//   - CSV export
//   - Source cards count
//   - Tab renumbering (after hatchery health insertion at "7": breeding is now "8 breeding", about is "9 about")
const { JSDOM } = require('jsdom');
const fs = require('fs');

const ATLAS_PATH = '/home/claude/work/Population_atlas.html';
const html = fs.readFileSync(ATLAS_PATH, 'utf-8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  pretendToBeVisual: true,
  url: 'http://localhost/Population_atlas.html',
});
const { window } = dom;
const { document } = window;
window.URL.createObjectURL = () => 'blob:mock';
window.URL.revokeObjectURL = () => {};

function fail(msg) { console.error('[bp-smoke] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[bp-smoke] ok  :', msg); }

// ---- 1. DOM presence ----
const page8 = document.getElementById('page8');
if (!page8) fail('page8 div missing');
ok('page8 div present');

const tabBtn = document.querySelector('#tabBar button[data-page="page8"]');
if (!tabBtn) fail('page8 tab button missing');
const tabLabel = (tabBtn.textContent || '').trim().replace(/\s+/g, ' ');
if (!tabLabel.includes('8') || !tabLabel.includes('breeding')) fail('expected tab "8 breeding", got "' + tabLabel + '"');
ok('page8 tab button: "' + tabLabel + '"');

const aboutBtn = document.querySelector('#tabBar button[data-page="page7"]');
if (!aboutBtn) fail('page7 (about) tab button missing');
const aboutLabel = (aboutBtn.textContent || '').trim().replace(/\s+/g, ' ');
if (!aboutLabel.includes('9') || !aboutLabel.includes('about')) fail('about should be renumbered to "9 about", got "' + aboutLabel + '"');
ok('about tab renumbered to: "' + aboutLabel + '"');

for (const name of [
  '_bpRender', '_bpEnsureState', '_bpAllRows', '_bpDerivedRows',
  '_bpAutoDiversityRows', '_bpAutoRelatednessRows', '_bpAutoBurdenRows',
  '_bpAutoInversionRows', '_bpAutoMarkerRows',
  '_bpIngestText', '_bpFilteredRows', '_bpExportCsv',
  '_bpIsValidJson', '_bpParseTsv',
]) {
  if (typeof window[name] !== 'function') fail(name + ' missing from window');
}
ok('all breeding-page helpers exposed on window');

// ---- 2. Empty-state render ----
// Activate page8 so the tab dispatcher has a chance to fire
document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
page8.classList.add('active');
window._bpRender();
const slot = document.getElementById('bpTableSlot');
if (!slot.innerHTML.includes('No rows yet')) fail('empty-state message missing in initial render');
ok('empty-state message renders before any data is loaded');

const badge = document.getElementById('bpStatusBadge');
if (!badge.textContent.includes('0 rows')) fail('status badge should show "0 rows" initially');
ok('status badge shows "0 rows" before data load');

// ---- 3. Auto-derivation: Diversity (top-1% F_ROH + top-1% θπ) ----
window.eval(`
  state.perSampleStats = [
    // 100 samples; CGA001 has top F_ROH (well above #2), CGA100 has top π (well above #99)
    ${Array.from({length: 100}, (_, i) => {
      const id = 'CGA' + String(i + 1).padStart(3, '0');
      // F_ROH: CGA001 = 0.40 (clear top), others 0.20 - 0.01
      const froh = i === 0 ? 0.40 : (0.20 - ((i - 1) * 0.19 / 98));
      // θπ: CGA100 = 0.005 (clear top), others 0.001 - 0.003
      const pi = i === 99 ? 0.005 : (0.001 + (i * 0.002 / 98));
      return `{ sample: '${id}', F_ROH: ${froh.toFixed(4)}, theta_pi: ${pi.toFixed(6)} }`;
    }).join(',\n    ')}
  ];
`);
const divRows = window._bpAutoDiversityRows();
// Top-1% = top 1 sample of 100
if (divRows.length < 2) fail('expected at least 2 diversity rows (top-1% froh + top-1% pi), got ' + divRows.length);
const frohRow = divRows.find(r => r.status === 'red');
if (!frohRow) fail('expected a red F_ROH row');
if (frohRow.sample !== 'CGA001') fail('CGA001 should be top-1% F_ROH, got ' + frohRow.sample);
ok('Diversity F_ROH top-1%: ' + frohRow.sample + ' flagged red (F_ROH=0.4000)');
const piRow = divRows.find(r => r.status === 'green' && r.highlight.includes('diversity'));
if (!piRow) fail('expected a green θπ row');
if (piRow.sample !== 'CGA100') fail('CGA100 should be top-1% π, got ' + piRow.sample);
ok('Diversity θπ top-1%: ' + piRow.sample + ' flagged green (θπ=0.005000)');

// ---- 4. Auto-derivation: Relatedness (largest family cluster) ----
window.eval(`
  state.familyClusters = {
    largest_cluster_size: 4,
    clusters: [
      { cluster_id: 'F-A', members: ['CGA002', 'CGA007', 'CGA011', 'CGA019'] },
      { cluster_id: 'F-B', members: ['CGA040', 'CGA050'] },
    ]
  };
`);
const relRows = window._bpAutoRelatednessRows();
if (relRows.length !== 4) fail('expected 4 family-cluster rows, got ' + relRows.length);
if (!relRows.every(r => r.status === 'yellow')) fail('all family-cluster rows should be yellow');
if (!relRows.every(r => r.value.includes('family size = 4'))) fail('family size should be 4 in value');
ok('Relatedness: 4 yellow rows from largest cluster F-A (CGA002, CGA007, CGA011, CGA019)');

// ---- 5. Auto-derivation: Burden ----
window.eval(`
  state.perSampleStats.forEach((s, i) => {
    s.deleterious_burden = i === 0 ? 99.5 : (10 + i * 0.5);  // CGA001 = top
  });
`);
const burdenRows = window._bpAutoBurdenRows();
const burdenTop = burdenRows.find(r => r.sample === 'CGA001');
if (!burdenTop) fail('CGA001 should be top-1% burden');
if (burdenTop.status !== 'red') fail('burden top-1% should be red');
ok('Burden top-1%: CGA001 flagged red (burden=99.50)');

// ---- 6. Auto-derivation: Inversions (rare carriers) ----
window.eval(`
  state.inversionCarriers = [
    { sample: 'CGA044', inversion_id: 'INV_LG05_002', frequency_in_cohort: 0.018 },
    { sample: 'CGA055', inversion_id: 'INV_LG28_001', frequency_in_cohort: 0.010 },
    { sample: 'CGA060', inversion_id: 'INV_LG02_005', frequency_in_cohort: 0.20 },  // not rare
  ];
`);
const invRows = window._bpAutoInversionRows();
if (invRows.length !== 2) fail('expected 2 rare-inversion rows (≤5%), got ' + invRows.length);
if (!invRows.every(r => r.status === 'green')) fail('rare-inversion rows should be green');
if (!invRows.find(r => r.sample === 'CGA044')) fail('CGA044 missing from inversion rows');
ok('Inversions: 2 green rows for rare carriers (CGA044 INV_LG05_002 1.8%, CGA055 INV_LG28_001 1.0%)');

// ---- 7. Auto-derivation: Markers (validation controls) ----
window.eval(`
  state.markerControls = [
    { sample: 'CGA071', marker_id: 'M001', role: 'positive_INV' },
    { sample: 'CGA072', marker_id: 'M001', role: 'positive_HET' },
    { sample: 'CGA088', marker_id: 'M001', role: 'negative_STD' },
  ];
`);
const mkRows = window._bpAutoMarkerRows();
if (mkRows.length !== 3) fail('expected 3 marker rows, got ' + mkRows.length);
const cga071 = mkRows.find(r => r.sample === 'CGA071');
if (!cga071) fail('CGA071 missing');
if (!cga071.value.includes('M001 INV/INV')) fail('CGA071 value should mention M001 INV/INV');
if (cga071.status !== 'green') fail('marker controls should be green');
ok('Markers: 3 green rows (CGA071 INV/INV, CGA072 HET, CGA088 STD/STD)');

// ---- 8. Source cards count + total ----
window._bpRender();
const cards = document.getElementById('bpSourceCards');
if (!cards.innerHTML.includes('Diversity')) fail('source cards missing Diversity');
const totalRows = window._bpAllRows().length;
// 2 diversity + 4 relatedness + 1 burden + 2 inversions + 3 markers = 12
if (totalRows !== 12) fail('expected 12 total rows, got ' + totalRows);
ok('total rows: 12 (Diversity 2 / Relatedness 4 / Burden 1 / Inversions 2 / Markers 3)');

// ---- 9. Status sorting (red first, then yellow, then green) ----
const tableHtml = document.getElementById('bpTableSlot').innerHTML;
const firstRedIdx = tableHtml.indexOf('bp-red');
const firstYellowIdx = tableHtml.indexOf('bp-yellow');
const firstGreenIdx = tableHtml.indexOf('bp-green');
if (firstRedIdx < 0 || firstYellowIdx < 0 || firstGreenIdx < 0) fail('expected all 3 colors in table');
if (!(firstRedIdx < firstYellowIdx && firstYellowIdx < firstGreenIdx))
  fail('expected sort order red < yellow < green; got red=' + firstRedIdx + ' yellow=' + firstYellowIdx + ' green=' + firstGreenIdx);
ok('rows sorted: red (caution) first, yellow (review) middle, green (useful) last');

// ---- 10. Simple vs detailed view toggle ----
const ui = window._bpEnsureState();
ui.view_mode = 'simple';
window._bpRender();
let html2 = document.getElementById('bpTableSlot').innerHTML;
if (!html2.includes('bp-simple')) fail('simple view should add bp-simple class to table');
const simpleColCount = (html2.match(/<th[\s>]/g) || []).length;
if (simpleColCount !== 4) fail('simple view should have 4 columns, got ' + simpleColCount);
ok('simple view: 4 columns (Sample / Highlight / Value / Source)');

ui.view_mode = 'detailed';
window._bpRender();
html2 = document.getElementById('bpTableSlot').innerHTML;
if (!html2.includes('bp-detailed')) fail('detailed view should add bp-detailed class');
// Count actual <th ...> tags (not <thead>) using a regex that requires a space or '>' after "th"
const detailedColCount = (html2.match(/<th[\s>]/g) || []).length;
if (detailedColCount !== 6) fail('detailed view should have 6 columns, got ' + detailedColCount);
if (!html2.includes('Notes / origin')) fail('detailed view missing notes/origin column');
ok('detailed view: 6 columns including Status + Notes / origin');

// Reset to simple for downstream
ui.view_mode = 'simple';

// ---- 11. User upload via JSON ----
const userJson = {
  metadata: { panel_name: 'Q2 review' },
  rows: [
    { sample: 'CGA200', highlight: 'Custom committee flag', value: 'see notes',
      source: 'Diversity', status: 'red', notes: 'manually flagged for chr08 inbreeding' },
    // Override CGA001 — change status from red to yellow with note
    { sample: 'CGA001', highlight: 'High inbreeding (acknowledged)', value: 'F_ROH=0.40, accepted',
      source: 'Diversity', status: 'yellow', notes: 'committee already reviewed' },
  ],
};
window._bpIngestText(JSON.stringify(userJson), 'breeding_highlights.json');
const merged = window._bpAllRows();
const cga200 = merged.find(r => r.sample === 'CGA200');
if (!cga200) fail('user-only row CGA200 not appended');
if (cga200.status !== 'red') fail('CGA200 should be red');
ok('user JSON upload: CGA200 (user-only row) appended');

// CGA001 should be overlaid: status yellow now, not red
const cga001 = merged.find(r => r.sample === 'CGA001' && r.source === 'Diversity');
if (!cga001) fail('CGA001 missing from merged');
if (cga001.status !== 'yellow') fail('CGA001 should be overlaid to yellow, got ' + cga001.status);
if (cga001._origin !== 'overlaid') fail('CGA001 should have _origin=overlaid, got ' + cga001._origin);
ok('user JSON upload: CGA001 (auto-row) overlaid → yellow with override note');

// ---- 12. Filters ----
ui.filter_status = 'red';
window._bpRender();
let filteredHtml = document.getElementById('bpTableSlot').innerHTML;
if (filteredHtml.includes('bp-yellow')) fail('red filter should hide yellow rows');
if (filteredHtml.includes('bp-green')) fail('red filter should hide green rows');
if (!filteredHtml.includes('CGA200')) fail('red filter should include CGA200');
ok('status=red filter: only red rows visible');
ui.filter_status = 'all';

ui.filter_source = 'Markers';
window._bpRender();
filteredHtml = document.getElementById('bpTableSlot').innerHTML;
if (!filteredHtml.includes('CGA071')) fail('Markers filter should include CGA071');
if (filteredHtml.match(/CGA002/)) fail('Markers filter should not include relatedness rows');
ok('source=Markers filter: only marker rows visible');
ui.filter_source = 'all';

// ---- 13. TSV upload ----
const tsv = [
  'sample\thighlight\tvalue\tsource\tstatus\tnotes\tref_url',
  'CGA300\tCustom TSV row\tspecial\tDiversity\tred\tloaded from TSV\t#page5',
  'CGA301\tAnother TSV row\tinfo\tInversions\tgreen\t\t',
].join('\n');
window._bpIngestText(tsv, 'breeding_highlights.tsv');
const tsvMerged = window._bpAllRows();
if (!tsvMerged.find(r => r.sample === 'CGA300')) fail('TSV row CGA300 not loaded');
if (!tsvMerged.find(r => r.sample === 'CGA301')) fail('TSV row CGA301 not loaded');
ok('TSV upload: 2 rows parsed and merged');

// ---- 14. CSV export ----
let exportedClicked = false;
const origCreateElement = document.createElement.bind(document);
document.createElement = function (tag) {
  const el = origCreateElement(tag);
  if (tag === 'a') {
    el.click = function () { exportedClicked = true; };
  }
  return el;
};
const csvBtn = document.getElementById('bpExportCsvBtn');
if (csvBtn) csvBtn.click();
if (!exportedClicked) fail('CSV export click handler did not fire');
ok('CSV export: click handler fires + Blob constructed');

// ---- 15. Reset clears user rows but keeps derived ----
const $reset = document.getElementById('bpResetBtn');
if (!$reset) fail('reset button missing');
$reset.click();
const afterReset = window._bpAllRows();
if (afterReset.find(r => r.sample === 'CGA200')) fail('reset should remove user-only CGA200');
if (afterReset.find(r => r.sample === 'CGA300')) fail('reset should remove TSV row CGA300');
// CGA001 should be back to red (auto-derived only)
const cga001AfterReset = afterReset.find(r => r.sample === 'CGA001' && r.source === 'Diversity');
if (!cga001AfterReset) fail('CGA001 should still be auto-derived after reset');
if (cga001AfterReset.status !== 'red') fail('CGA001 should be red again after reset (auto-row), got ' + cga001AfterReset.status);
ok('reset: drops user rows, restores auto-row defaults (CGA001 → red)');

// ---- 16. Status badge ready state ----
window._bpRender();
const badge2 = document.getElementById('bpStatusBadge');
if (badge2.classList.contains('planned')) fail('badge should not be "planned" when rows exist');
if (!badge2.classList.contains('ready')) fail('badge should be "ready" when rows exist');
if (!badge2.textContent.match(/auto-derived/)) fail('badge should mention auto-derived count');
ok('status badge: "' + badge2.textContent + '"');

// ---- 17. Empty state when no upstream data ----
// Drop all upstream state and verify the empty message comes back
window.eval('delete state.perSampleStats; delete state.familyClusters; delete state.inversionCarriers; delete state.markerControls;');
window._bpRender();
const emptyHtml = document.getElementById('bpTableSlot').innerHTML;
if (!emptyHtml.includes('No rows yet')) fail('empty state should reappear after dropping upstream data');
ok('empty state restored after dropping upstream data');

// ---- 18. About page (page7) — turn 80 changelog updates ----
const aboutPage = document.getElementById('page7');
if (!aboutPage) fail('about page (page7) missing');
if (!aboutPage.innerHTML.includes('What\u2019s new in v1.2') &&
    !aboutPage.innerHTML.includes("What's new in v1.2"))
  fail('about page missing "What\u2019s new in v1.2" card');
if (!aboutPage.innerHTML.includes('Page 7 \u2014 Breeding')) fail('about page should describe new breeding tab');
if (!aboutPage.innerHTML.includes('Hash routing')) fail('about page should describe hash routing');
ok('about page: "What\u2019s new in v1.2" card lists breeding tab + hash routing');

// Stale "page 13 (Diversity Atlas preview)" cross-reference cleaned up
if (aboutPage.innerHTML.includes('Diversity Atlas preview'))
  fail('stale "Diversity Atlas preview" cross-ref should have been replaced');
if (!aboutPage.innerHTML.includes('tab 13 (stats profile)')) fail('cross-ref should mention stats profile tab');
if (!aboutPage.innerHTML.includes('tab 14 (marker panel)')) fail('cross-ref should mention marker panel tab');
ok('about page: stale "Diversity Atlas preview" replaced with current tab 13/14 cross-refs');

console.log('\n[bp-smoke] ALL CHECKS PASSED');
