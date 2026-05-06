// Smoke test for the F_ROH × H plane / hatchery health page (Population
// Atlas page 9, visible label "7 hatchery health"). Tests:
//   - DOM presence (page9 div, "7 hatchery health" tab, renumbering)
//   - Quadrant classification (4 corners)
//   - H_ref resolution (explicit, species-peer fallback, global fallback)
//   - Focal-cohort detection (explicit flag, highest-F_ROH fallback)
//   - Empty state
//   - JSON upload + plot render + table render + summary cards
//   - Verdict text uses correct colored CSS class
//   - Breeding page auto-row appears with colored highlight_html + correct status
//   - CSV export
//   - TSV upload
//   - Reset

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

function fail(msg) { console.error('[fhh-smoke] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[fhh-smoke] ok  :', msg); }

// ---- 1. DOM presence + renumbering ----
const page9 = document.getElementById('page9');
if (!page9) fail('page9 (hatchery health) div missing');
ok('page9 hatchery health div present');

const fhhTab = document.querySelector('#tabBar button[data-page="page9"]');
if (!fhhTab) fail('hatchery health tab missing');
const fhhLabel = (fhhTab.textContent || '').trim().replace(/\s+/g, ' ');
if (!fhhLabel.includes('7') || !fhhLabel.toLowerCase().includes('hatchery')) fail('expected "7 hatchery health", got "' + fhhLabel + '"');
ok('hatchery health tab: "' + fhhLabel + '"');

// Breeding should now be "8"
const bpTab = document.querySelector('#tabBar button[data-page="page8"]');
if (!bpTab) fail('breeding tab missing');
const bpLabel = (bpTab.textContent || '').trim().replace(/\s+/g, ' ');
if (!bpLabel.startsWith('8') || !bpLabel.toLowerCase().includes('breeding')) fail('expected "8 breeding", got "' + bpLabel + '"');
ok('breeding tab renumbered to: "' + bpLabel + '"');

// About should now be "9"
const aboutTab = document.querySelector('#tabBar button[data-page="page7"]');
if (!aboutTab) fail('about tab missing');
const aboutLabel = (aboutTab.textContent || '').trim().replace(/\s+/g, ' ');
if (!aboutLabel.startsWith('9') || !aboutLabel.toLowerCase().includes('about')) fail('expected "9 about", got "' + aboutLabel + '"');
ok('about tab renumbered to: "' + aboutLabel + '"');

// ---- 2. Helpers exposed ----
for (const name of ['_fhhRender', '_fhhEnsureState', '_fhhClassify',
                    '_fhhFocalCohort', '_fhhResolveHRef', '_fhhIngestText',
                    '_fhhExportCsv', '_fhhParseTsv', '_fhhIsValidJson',
                    '_bpAutoHatcheryHealthRow']) {
  if (typeof window[name] !== 'function') fail(name + ' not exposed');
}
if (!window.FHH_VERDICTS) fail('FHH_VERDICTS not exposed');
ok('all hatchery health helpers exposed on window');

// ---- 3. Quadrant classification (4 corners) ----
// Use a 4-cohort test set so medians are well-defined: F_ROH median = 0.15,
// H/H_ref median = 0.6 (using H_ref = 0.01 for everyone for clean math).
const testCohorts4 = [
  { cohort: 'A_low_high',  F_ROH: 0.05, H_per_site: 0.009, H_ref: 0.01 },  // healthy
  { cohort: 'B_high_high', F_ROH: 0.30, H_per_site: 0.008, H_ref: 0.01 },  // recent inbreeding only
  { cohort: 'C_low_low',   F_ROH: 0.05, H_per_site: 0.003, H_ref: 0.01 },  // erosion only
  { cohort: 'D_high_low',  F_ROH: 0.30, H_per_site: 0.002, H_ref: 0.01 },  // caution
];

const ra = window._fhhClassify(testCohorts4[0], testCohorts4);
if (ra.verdict !== 'good') fail('A (low F_ROH, high H/H_ref) should be "good", got ' + ra.verdict);
ok('quadrant A (low F_ROH 0.05, high H/H_ref 0.9) → good (Healthy)');

const rb = window._fhhClassify(testCohorts4[1], testCohorts4);
if (rb.verdict !== 'review_inbreeding') fail('B should be "review_inbreeding", got ' + rb.verdict);
ok('quadrant B (high F_ROH 0.30, high H/H_ref 0.8) → review_inbreeding');

const rc = window._fhhClassify(testCohorts4[2], testCohorts4);
if (rc.verdict !== 'review_erosion') fail('C should be "review_erosion", got ' + rc.verdict);
ok('quadrant C (low F_ROH 0.05, low H/H_ref 0.3) → review_erosion');

const rd = window._fhhClassify(testCohorts4[3], testCohorts4);
if (rd.verdict !== 'caution') fail('D should be "caution", got ' + rd.verdict);
ok('quadrant D (high F_ROH 0.30, low H/H_ref 0.2) → caution');

// ---- 4. H_ref resolution ----
// Explicit H_ref wins
const hrefExplicit = window._fhhResolveHRef(
  { species: 'X', H_per_site: 0.005, H_ref: 0.012 }, []);
if (hrefExplicit !== 0.012) fail('explicit H_ref should win, got ' + hrefExplicit);
ok('H_ref resolution: explicit per-cohort H_ref wins');

// Species-peer fallback
const peerCohorts = [
  { species: 'X', H_per_site: 0.020, H_ref: undefined },
  { species: 'X', H_per_site: 0.005 },
];
const hrefPeer = window._fhhResolveHRef({ species: 'X', H_per_site: 0.005 }, peerCohorts);
if (hrefPeer !== 0.020) fail('peer fallback should give 0.020, got ' + hrefPeer);
ok('H_ref resolution: species-peer fallback uses max H of same species');

// Global fallback (no species match)
const hrefGlobal = window._fhhResolveHRef(
  { species: 'Y', H_per_site: 0.005 }, peerCohorts);
if (hrefGlobal !== 0.020) fail('global fallback should give 0.020, got ' + hrefGlobal);
ok('H_ref resolution: falls back to global max when no species match');

// ---- 5. Focal-cohort detection ----
window._fhhEnsureState().cohorts = testCohorts4.slice();
// No is_focal flagged → falls back to highest F_ROH
let focal = window._fhhFocalCohort();
if (focal.cohort !== 'B_high_high' && focal.cohort !== 'D_high_low') fail('focal fallback should pick a high-F_ROH cohort, got ' + focal.cohort);
ok('focal fallback: picks highest-F_ROH cohort (' + focal.cohort + ')');

// Explicit is_focal flag wins
window._fhhEnsureState().cohorts = testCohorts4.map((c, i) =>
  Object.assign({}, c, { is_focal: i === 0 }));
focal = window._fhhFocalCohort();
if (focal.cohort !== 'A_low_high') fail('explicit is_focal should win, got ' + focal.cohort);
ok('focal: explicit is_focal flag wins (' + focal.cohort + ')');

// ---- 6. Empty state ----
window._fhhEnsureState().cohorts = [];
window._fhhRender();
const plotSlot = document.getElementById('fhhPlotSlot');
if (!plotSlot.innerHTML.includes('No cohorts loaded')) fail('empty plot state missing');
const tableSlot = document.getElementById('fhhTableSlot');
if (!tableSlot.innerHTML.includes('No cohorts loaded')) fail('empty table state missing');
const badge = document.getElementById('fhhStatusBadge');
if (!badge.classList.contains('planned')) fail('badge should be planned in empty state');
ok('empty state: plot + table show "No cohorts loaded", badge stays planned');

// ---- 7. JSON upload + render ----
const jsonPayload = {
  metadata: { panel_name: 'F_ROH × H plane', date: '2026-04-30' },
  cohorts: [
    { cohort: 'this study (226 broodstock)',
      species: 'C. gariepinus', is_focal: true,
      F_ROH: 0.277, H_per_site: 0.00455, H_ref: 0.0156,
      notes: 'Thai hatchery, 9x WGS, ngsF-HMM ≥1 Mb' },
    { cohort: 'Nigerian wild C. gariepinus',
      species: 'C. gariepinus',
      F_ROH: 0.05, H_per_site: 0.0156, H_ref: 0.0156,
      citation: 'Sanda 2026' },
    { cohort: 'rainbow trout (D Ambrosio)',
      species: 'O. mykiss',
      F_ROH: 0.20, H_per_site: 0.0030, H_ref: 0.005,
      citation: "D'Ambrosio 2019" },
    { cohort: 'alligator (Glenn et al)',
      species: 'A. mississippiensis',
      F_ROH: 0.08, H_per_site: 0.0040, H_ref: 0.008 },
  ],
};
window._fhhIngestText(JSON.stringify(jsonPayload), 'hatchery_health.json');

const fhhState = window._fhhEnsureState();
if (fhhState.cohorts.length !== 4) fail('expected 4 cohorts, got ' + fhhState.cohorts.length);
ok('JSON upload: 4 cohorts loaded');

// ---- 8. Plot SVG renders with correct colored points ----
const plotHtml = document.getElementById('fhhPlotSlot').innerHTML;
if (!plotHtml.includes('<svg')) fail('SVG plot missing after upload');
if (!plotHtml.includes('class="fhh-plot-point fhh-caution')) fail('SVG should include caution-class point (focal cohort)');
ok('SVG plot rendered with verdict-colored points');

// Focal cohort gets the bigger radius + bold label
if (!plotHtml.includes('fhh-focal')) fail('focal cohort should have fhh-focal class');
if (!plotHtml.includes('fhh-plot-label-focal')) fail('focal label should have bold class');
ok('focal cohort styled distinctly (larger point + bold label)');

// ---- 9. Table renders with colored verdict text ----
const tableHtml = document.getElementById('fhhTableSlot').innerHTML;
if (!tableHtml.includes('<table class="fhh-table">')) fail('table not rendered');
if (!tableHtml.includes('fhh-caution')) fail('caution-class verdict missing in table');
if (!tableHtml.includes('Caution \u2014 historical erosion')) fail('caution verdict label missing');
if (!tableHtml.includes('FOCAL')) fail('FOCAL marker missing in table');
ok('table rendered with FOCAL marker + colored verdict text');

// Cohort name + numeric values stay in normal ink (NOT in colored span)
const focalRow = tableHtml.match(/<tr>[^<]*<td class="fhh-cohort">[^<]*<span class="fhh-focal-marker">FOCAL<\/span>[^<]*this study[\s\S]*?<\/tr>/);
if (!focalRow) fail('could not isolate focal row from table');
// Confirm the cohort name itself is NOT wrapped in fhh-caution
const cohortCell = focalRow[0].match(/<td class="fhh-cohort">[\s\S]*?<\/td>/)[0];
if (cohortCell.includes('fhh-caution') || cohortCell.includes('fhh-good') || cohortCell.includes('fhh-review')) {
  fail('cohort name cell should NOT have verdict color');
}
ok('cohort name in normal ink (not colored); only verdict text is colored');

// ---- 10. Summary cards populated ----
const cardsHtml = document.getElementById('fhhSummaryCards').innerHTML;
if (!cardsHtml.includes('fhh-good')) fail('cards missing fhh-good');
if (!cardsHtml.includes('fhh-caution')) fail('cards missing fhh-caution');
if (!cardsHtml.includes('Cohorts loaded')) fail('cards missing total count');
ok('summary cards rendered with quadrant counts');

// Status badge updated
const badge2 = document.getElementById('fhhStatusBadge');
if (badge2.classList.contains('planned')) fail('badge should not be planned after load');
if (!badge2.innerHTML.includes('focal:')) fail('badge should mention focal cohort verdict');
if (!badge2.innerHTML.includes('fhh-')) fail('badge should include colored verdict span');
ok('status badge shows colored focal verdict');

// ---- 11. Breeding page auto-row appears + carries colored highlight ----
const bpRows = window._bpAutoHatcheryHealthRow();
if (bpRows.length !== 1) fail('expected 1 hatchery health row, got ' + bpRows.length);
const hRow = bpRows[0];
if (hRow.source !== 'Diversity') fail('breeding row source should be "Diversity"');
if (hRow.status !== 'red') fail('focal=this study (caution) should map to status=red, got ' + hRow.status);
if (!hRow.highlight_html || !hRow.highlight_html.includes('fhh-caution')) {
  fail('highlight_html should carry fhh-caution class');
}
if (!hRow.value.includes('F_ROH 0.277')) fail('breeding row value should include F_ROH 0.277');
if (!hRow.ref_url || !hRow.ref_url.includes('page9')) fail('ref_url should point to page9');
ok('breeding auto-row: source=Diversity, status=red, highlight_html carries fhh-caution');

// ---- 12. Breeding page renders the auto-row with colored verdict ----
window._bpRender();
const bpSlot = document.getElementById('bpTableSlot');
if (!bpSlot.innerHTML.includes('this study (226 broodstock)')) fail('breeding table missing focal cohort name');
if (!bpSlot.innerHTML.includes('class="fhh-caution"')) fail('breeding table should preserve fhh-caution span from highlight_html');
ok('breeding page table includes hatchery health row with colored verdict text');

// ---- 13. CSV export ----
let exportClicked = false;
const origCreateElement = document.createElement.bind(document);
document.createElement = function (tag) {
  const el = origCreateElement(tag);
  if (tag === 'a') el.click = function () { exportClicked = true; };
  return el;
};
const csvBtn = document.getElementById('fhhExportCsvBtn');
if (csvBtn) csvBtn.click();
if (!exportClicked) fail('CSV export click did not fire');
ok('CSV export click handler fires');

// ---- 14. TSV upload ----
const tsv = [
  'cohort\tspecies\tF_ROH\tH_per_site\tH_ref\tcitation\tnotes\tis_focal',
  'cohort_X\tC. gariepinus\t0.18\t0.006\t0.012\tTSV test\tnotes here\tfalse',
].join('\n');
const tsvParsed = window._fhhParseTsv(tsv);
if (!tsvParsed || tsvParsed.cohorts.length !== 1) fail('TSV parse failed');
if (tsvParsed.cohorts[0].F_ROH !== 0.18) fail('TSV F_ROH wrong');
ok('TSV parser: 1 cohort with F_ROH=0.18 parsed correctly');

// ---- 15. Reset ----
const resetBtn = document.getElementById('fhhResetBtn');
if (resetBtn) resetBtn.click();
const afterReset = window._fhhEnsureState();
if (afterReset.cohorts.length !== 0) fail('reset should clear cohorts');
ok('reset: cohorts cleared');

// After reset, the breeding-page hatchery health row should also disappear
window._bpRender();
const bpAfterReset = document.getElementById('bpTableSlot').innerHTML;
if (bpAfterReset.includes('class="fhh-caution"')) fail('breeding row should be removed after FHH reset');
ok('breeding page: hatchery health row removed after FHH reset');

// ---- 16. CSS classes for colored verdicts present ----
const styleText = Array.from(document.querySelectorAll('style')).map(s => s.textContent).join('\n');
for (const cls of ['fhh-good', 'fhh-review', 'fhh-caution']) {
  if (!styleText.includes('.' + cls + ' ')) fail('CSS missing rule for .' + cls);
}
ok('CSS rules present for fhh-good / fhh-review / fhh-caution');

// ---- 17. Hash routing — #page9 lands on hatchery health ----
// (This is testable because the existing hash-routing infrastructure runs
// any matching button click on hashchange.)
window.location.hash = '#page9';
window.dispatchEvent(new window.Event('hashchange'));
// Give it a moment
setTimeout(() => {
  const active = document.querySelector('.page.active');
  if (!active || active.id !== 'page9') {
    console.error('[fhh-smoke] FAIL: #page9 hash should activate page9, got #' + (active ? active.id : 'none'));
    process.exit(1);
  }
  ok('#page9 hash routes to hatchery health page');

  console.log('\n[fhh-smoke] ALL CHECKS PASSED');
}, 30);
