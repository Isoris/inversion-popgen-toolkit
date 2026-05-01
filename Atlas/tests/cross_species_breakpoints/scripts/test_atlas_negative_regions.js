// Tests for v4 turn 82:
//   1. Caution banner on Population Atlas hatchery health page
//   2. Caution banner on Inversion Atlas negative regions page
//   3. Negative regions page (page19, visible "6 negative regions")
//      - Tab insertion + stage band
//      - Renumbering of downstream tabs (karyotype 7, popstats 8, ..., help 16)
//      - 5 region_status values + colored pills
//      - Evidence layer mini-pills (pass / limited / fail)
//      - JSON upload + table render + summary cards
//      - TSV upload
//      - CSV export
//      - Reset
//      - Methods note + caution banner present

const { JSDOM } = require('jsdom');
const fs = require('fs');

function fail(msg) { console.error('[caution-nr-smoke] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[caution-nr-smoke] ok  :', msg); }

// ============================================================================
// Suite 1 — Caution banner on Population Atlas hatchery health page
// ============================================================================
console.log('\n=== Suite 1: Population Atlas hatchery health caution banner ===');
{
  const html = fs.readFileSync('/home/claude/work/Population_atlas.html', 'utf-8');
  const dom = new JSDOM(html, {
    runScripts: 'dangerously', pretendToBeVisual: true,
    url: 'http://localhost/Population_atlas.html',
  });
  const { window: win } = dom;
  win.URL.createObjectURL = () => 'blob:mock';
  win.URL.revokeObjectURL = () => {};

  const page9 = win.document.getElementById('page9');
  if (!page9) fail('hatchery health page (page9) missing');

  const banner = page9.querySelector('.caution-banner');
  if (!banner) fail('caution banner missing on hatchery health page');
  ok('caution banner present on hatchery health page');

  const bannerText = banner.textContent;
  if (!bannerText.includes('New feature')) fail('banner should mention "New feature"');
  if (!bannerText.includes('no published comparator')) fail('banner should explain no comparator baseline');
  if (!bannerText.includes('novel framework')) fail('banner should call F_ROH | H novel framework');
  ok('banner text mentions: new feature, novel framework, no published comparator');

  // Banner uses the warning icon (⚠ U+26A0)
  if (!banner.querySelector('.caution-banner-icon')) fail('banner icon missing');
  const iconText = banner.querySelector('.caution-banner-icon').textContent;
  if (!iconText.includes('\u26a0')) fail('banner icon should be ⚠ (U+26A0)');
  ok('banner icon: ⚠ (U+26A0)');

  // CSS class for amber/accent theme
  const styleText = Array.from(win.document.querySelectorAll('style')).map(s => s.textContent).join('\n');
  if (!styleText.includes('.caution-banner ')) fail('CSS rule for .caution-banner missing');
  if (!styleText.includes('var(--accent)')) fail('caution banner should use --accent token');
  ok('CSS rule for .caution-banner present (uses --accent token)');

  // The v1.3 What's new card should have its escape codes properly rendered now
  const aboutPage = win.document.getElementById('page7');
  const aboutText = aboutPage ? aboutPage.textContent : '';
  if (aboutText.includes('\\u2014') || aboutText.includes('\\u00d7')) {
    fail('v1.3 What\u2019s new card should not contain literal escape sequences (rendered text should use real unicode)');
  }
  ok('v1.3 What\u2019s new card: literal \\u escape codes converted to real unicode');
}

// ============================================================================
// Suite 2 — Negative regions page on Inversion Atlas
// ============================================================================
console.log('\n=== Suite 2: Inversion Atlas negative regions page ===');
{
  const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
  const dom = new JSDOM(html, {
    runScripts: 'dangerously', pretendToBeVisual: true,
    url: 'http://localhost/Inversion_atlas.html',
  });
  const { window: win } = dom;
  win.URL.createObjectURL = () => 'blob:mock';
  win.URL.revokeObjectURL = () => {};

  // ---- DOM presence ----
  const page19 = win.document.getElementById('page19');
  if (!page19) fail('negative regions page (page19) div missing');
  ok('page19 negative regions div present');

  const nrTab = win.document.querySelector('#tabBar button[data-page="page19"]');
  if (!nrTab) fail('negative regions tab missing');
  const tabLabel = (nrTab.textContent || '').trim().replace(/\s+/g, ' ');
  if (!tabLabel.startsWith('6') || !tabLabel.toLowerCase().includes('negative')) {
    fail('expected tab "6 negative regions", got "' + tabLabel + '"');
  }
  ok('tab label: "' + tabLabel + '"');

  // Stage band
  if (nrTab.dataset.stage !== 'refinement') {
    fail('negative regions tab should have data-stage="refinement", got "' + nrTab.dataset.stage + '"');
  }
  ok('negative regions tab has data-stage="refinement" (cyan band)');

  // ---- Renumbering of downstream tabs ----
  const expectedRenumbering = [
    { page: 'page4', expectedStart: '7', label: 'karyotype' },
    { page: 'page6', expectedStart: '8', label: 'popstats' },
    { page: 'page7', expectedStart: '9', label: 'ancestry' },
    { page: 'page8', expectedStart: '10', label: 'windows' },
    { page: 'page9', expectedStart: '11', label: 'confirmed' },
    { page: 'page10', expectedStart: '12', label: 'markers' },
    { page: 'page16', expectedStart: '13', label: 'cross-species' },
    { page: 'page17', expectedStart: '14', label: 'stats profile' },
    { page: 'page18', expectedStart: '15', label: 'marker panel' },
    { page: 'page5', expectedStart: '16', label: 'help' },
  ];
  for (const e of expectedRenumbering) {
    const btn = win.document.querySelector('#tabBar button[data-page="' + e.page + '"]');
    if (!btn) fail('tab ' + e.page + ' missing');
    const numSpan = btn.querySelector('.num');
    if (!numSpan || numSpan.textContent.trim() !== e.expectedStart) {
      fail(e.page + ' (' + e.label + ') should be renumbered to "' + e.expectedStart + '", got "' +
        (numSpan ? numSpan.textContent : '<no .num>') + '"');
    }
  }
  ok('all 10 downstream tabs renumbered correctly (karyotype=7, ..., help=16)');

  // ---- Caution banner on negative regions page ----
  const banner = page19.querySelector('.caution-banner');
  if (!banner) fail('caution banner missing on negative regions page');
  const bannerText = banner.textContent;
  if (!bannerText.includes('never prove')) fail('banner should explain "never prove" 100% absence');
  if (!bannerText.includes('region_status')) fail('banner should reference region_status field');
  ok('caution banner: "never prove no inversions with 100% certainty" framing present');

  // ---- 5 region_status values defined ----
  if (typeof win.NR_STATUS_DEFS !== 'object') fail('NR_STATUS_DEFS not exposed');
  const expectedStatus = [
    'inversion_candidate',
    'complex_candidate',
    'no_detectable_inversion_high_confidence',
    'no_detectable_inversion_low_power',
    'not_testable',
  ];
  for (const s of expectedStatus) {
    if (!win.NR_STATUS_DEFS[s]) fail('NR_STATUS_DEFS missing key ' + s);
  }
  ok('all 5 region_status values defined (inversion_candidate, complex_candidate, high-conf neg, low-power neg, not_testable)');

  // ---- Vocabulary card lists all 5 statuses ----
  const vocabHtml = page19.innerHTML;
  for (const s of expectedStatus) {
    if (!vocabHtml.includes(s)) fail('vocabulary card missing status: ' + s);
  }
  ok('vocabulary card lists all 5 region_status values with code highlighting');

  // ---- 8 evidence layers defined ----
  if (!Array.isArray(win.NR_EVIDENCE_LAYERS)) fail('NR_EVIDENCE_LAYERS not exposed');
  if (win.NR_EVIDENCE_LAYERS.length !== 8) {
    fail('expected 8 evidence layers, got ' + win.NR_EVIDENCE_LAYERS.length);
  }
  const layerKeys = win.NR_EVIDENCE_LAYERS.map(L => L.key);
  for (const k of ['local_pca', 'ghsl', 'sv_callers', 'heterozygosity', 'ld', 'synteny', 'read_depth', 'callable_mask']) {
    if (!layerKeys.includes(k)) fail('missing evidence layer: ' + k);
  }
  ok('all 8 evidence layers defined: ' + layerKeys.join(', '));

  // Evidence layer card lists them all
  for (const L of win.NR_EVIDENCE_LAYERS) {
    if (!vocabHtml.includes(L.label)) fail('evidence layer card missing label: ' + L.label);
  }
  ok('evidence layers card lists all 8 layers with expected patterns');

  // ---- Helpers exposed ----
  for (const name of ['_nrRender', '_nrEnsureState', '_nrIngestText', '_nrParseTsv', '_nrExportCsv']) {
    if (typeof win[name] !== 'function') fail(name + ' not exposed');
  }
  ok('all 5 helpers exposed: _nrRender, _nrEnsureState, _nrIngestText, _nrParseTsv, _nrExportCsv');

  // ---- Empty state ----
  win._nrEnsureState().regions = [];
  win._nrRender();
  const tableSlot = win.document.getElementById('nrTableSlot');
  if (!tableSlot.innerHTML.includes('No regions loaded')) fail('empty state missing');
  ok('empty state: "No regions loaded" message shown');

  // ---- JSON upload ----
  const jsonPayload = {
    metadata: { panel_name: 'LG28 sweep', date: '2026-04-30' },
    regions: [
      { region_id: 'LG28:0-1.5Mb', chr: 'LG28', start_bp: 0, end_bp: 1500000,
        region_status: 'no_detectable_inversion_high_confidence',
        evidence: { local_pca: 'pass', ghsl: 'pass', sv_callers: 'pass', heterozygosity: 'pass', callable_mask: 'pass' },
        snp_density: 0.012, callable_fraction: 0.94, n_samples: 226,
        notes: 'stable, high-coverage region' },
      { region_id: 'LG28:5-7.5Mb', chr: 'LG28', start_bp: 5000000, end_bp: 7500000,
        region_status: 'no_detectable_inversion_low_power',
        evidence: { local_pca: 'pass', ghsl: 'limited', sv_callers: 'pass', heterozygosity: 'pass', callable_mask: 'limited' },
        snp_density: 0.0035, callable_fraction: 0.62,
        notes: 'repetitive, low SNP density' },
      { region_id: 'LG28:10-11Mb', chr: 'LG28', start_bp: 10000000, end_bp: 11000000,
        region_status: 'inversion_candidate',
        evidence: { local_pca: 'pass', ghsl: 'pass', sv_callers: 'pass' },
        notes: 'cross-listed positive call' },
      { region_id: 'LG28:14-15Mb', chr: 'LG28', start_bp: 14000000, end_bp: 15000000,
        region_status: 'complex_candidate',
        evidence: { local_pca: 'pass', ghsl: 'pass' },
        notes: 'multi-arrangement region' },
      { region_id: 'LG28:18-19Mb', chr: 'LG28', start_bp: 18000000, end_bp: 19000000,
        region_status: 'not_testable',
        evidence: {},
        notes: 'assembly gap' },
    ],
  };
  win._nrIngestText(JSON.stringify(jsonPayload), 'negative_regions.json');
  if (win._nrEnsureState().regions.length !== 5) fail('expected 5 regions, got ' + win._nrEnsureState().regions.length);
  ok('JSON upload: 5 regions loaded covering all 5 region_status values');

  // ---- Table renders all 5 status pills ----
  const tableHtml = win.document.getElementById('nrTableSlot').innerHTML;
  if (!tableHtml.includes('<table class="nr-table">')) fail('table not rendered');
  for (const s of expectedStatus) {
    const def = win.NR_STATUS_DEFS[s];
    if (!tableHtml.includes(def.klass)) fail('table missing pill class for: ' + s + ' (' + def.klass + ')');
  }
  ok('table renders all 5 status pills with correct CSS classes (positive/complex/good/low-power/untested)');

  // ---- Evidence mini-pills ----
  if (!tableHtml.includes('class="nr-evidence-mini pass"')) fail('table missing pass-class evidence pill');
  if (!tableHtml.includes('class="nr-evidence-mini limited"')) fail('table missing limited-class evidence pill');
  ok('table renders evidence mini-pills with pass/limited classes');

  // ---- Summary cards ----
  const cardsHtml = win.document.getElementById('nrSummaryCards').innerHTML;
  if (!cardsHtml.includes('Total regions')) fail('summary cards missing total');
  if (!cardsHtml.includes('High-confidence negative')) fail('summary cards missing high-conf neg label');
  if (!cardsHtml.includes('Negative \u2014 limited power')) fail('summary cards missing low-power label');
  ok('summary cards: Total regions + per-status counts');

  // Status badge updated
  const tableBadge = win.document.getElementById('nrTableBadge');
  if (!tableBadge.textContent.includes('5 regions loaded')) fail('table badge should show "5 regions loaded"');
  ok('table badge shows "5 regions loaded"');

  // ---- TSV upload ----
  const tsv = [
    'region_id\tchr\tstart_bp\tend_bp\tregion_status\tsnp_density\tcallable_fraction\tn_samples\tnotes',
    'TSV:1-2Mb\tLG01\t1000000\t2000000\tno_detectable_inversion_high_confidence\t0.014\t0.96\t226\tTSV test',
  ].join('\n');
  const tsvParsed = win._nrParseTsv(tsv);
  if (!tsvParsed || tsvParsed.regions.length !== 1) fail('TSV parse failed');
  if (tsvParsed.regions[0].region_status !== 'no_detectable_inversion_high_confidence') fail('TSV region_status wrong');
  ok('TSV parser: 1 region with high-conf neg status');

  // ---- CSV export ----
  let exportClicked = false;
  const origCreateElement = win.document.createElement.bind(win.document);
  win.document.createElement = function (tag) {
    const el = origCreateElement(tag);
    if (tag === 'a') el.click = function () { exportClicked = true; };
    return el;
  };
  const csvBtn = win.document.getElementById('nrExportCsvBtn');
  if (csvBtn) csvBtn.click();
  if (!exportClicked) fail('CSV export click did not fire');
  ok('CSV export click handler fires');

  // ---- Reset ----
  const resetBtn = win.document.getElementById('nrResetBtn');
  if (resetBtn) resetBtn.click();
  if (win._nrEnsureState().regions.length !== 0) fail('reset should clear regions');
  ok('reset: regions cleared');

  // ---- Methods note + manuscript wording ----
  // textContent strips HTML and joins runs; check on the flat text not innerHTML
  const page19Text = page19.textContent.replace(/\s+/g, ' ');
  if (!page19Text.includes('inversion-negative')) fail('page should use the term "inversion-negative"');
  if (!page19Text.includes('not interpreted as definitively')) fail('methods note should include "not interpreted as definitively inversion-free" wording');
  if (!page19Text.includes('Absence calls indicate')) fail('Discussion sentence ("Absence calls indicate...") should be present');
  ok('methods note includes both Methods + Discussion manuscript wordings');

  // ---- CSS rules present ----
  const invStyleText = Array.from(win.document.querySelectorAll('style')).map(s => s.textContent).join('\n');
  for (const cls of ['nr-status-positive', 'nr-status-complex', 'nr-status-good', 'nr-status-low-power', 'nr-status-untested']) {
    if (!invStyleText.includes('.' + cls + ' ')) fail('CSS missing rule for .' + cls);
  }
  ok('CSS rules for all 5 status pill classes present');

  // ---- Caution banner CSS reused ----
  if (!invStyleText.includes('.caution-banner ')) fail('caution banner CSS missing on inversion atlas');
  ok('caution-banner CSS rule present in inversion atlas (consistent with population atlas)');
}

console.log('\n[caution-nr-smoke] ALL CHECKS PASSED');
