// =============================================================================
// turn 132 — Slice 1: page 12 (θπ scrubber) frontend panel scaffold mirrors
// page 1, with one extra hero panel (#thCusumHeroPanel) for the per-carrier
// CUSUM boundary-distribution visualization.
//
// This slice ships the DOM scaffold only:
//   - #thCtrlBar (mirror of #ctrlBar)
//   - #thCusumHeroPanel (slightly-different-from-page-1 hero)
//   - #thSimPanel (mirror of #simPanel)
//   - #thZPanel (mirror of #zPanel)
//   - #thLinesPanel (mirror of #linesPanel)
//   - #thAnchorStripPanel (mirror of #anchorStripPanel)
//   - #thPcaPanel (mirror of #pcaPanel)
//   - #thTrackedSamplesPanelCompact (mirror of #trackedSamplesPanelCompact)
//   - #thL3Panel (mirror of #l3Panel)
//
// Plus: cusum_theta layer added to #thetaPiLayerStatus grid (R-side
// STEP_T05). Plus: panel-layout description grid extended to list the
// hero + anchor-strip + tracked-samples panels.
//
// All panels start display: none. Visibility wiring + renderers ship in
// later slices when the corresponding theta_pi_* / cusum_theta data lands.
// =============================================================================

const fs = require('fs');
const path = require('path');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Empty-state block still present (we add panels alongside, never remove)
// ============================================================================
console.log('\n=== 1. Empty-state preservation ===');

ok('#thetaPiEmpty empty-state block still exists',
   /id="thetaPiEmpty"/.test(html));

ok('#thetaPiEmpty contains "local PCA θπ" header',
   /id="thetaPiEmpty"[\s\S]{0,800}local PCA θπ — chromosome-wide diversity scanner/.test(html));

ok('#thetaPiLayerStatus grid still present',
   /id="thetaPiLayerStatus"/.test(html));

// ============================================================================
// 2. cusum_theta added to layer-status grid (4th row)
// ============================================================================
console.log('\n=== 2. cusum_theta layer row added ===');

ok('#thetaPiLayerStatus has cusum_theta row',
   /data-th-layer="cusum_theta"/.test(html));

ok('cusum_theta row references STEP_T05',
   /<code>cusum_theta<\/code>[\s\S]{0,300}<code>STEP_T05<\/code>/.test(html));

ok('cusum_theta row documents hero panel feed',
   /data-th-layer="cusum_theta"[\s\S]{0,400}hero panel · per-carrier CUSUM strip \+ boundary histogram/.test(html));

ok('cusum_theta row starts as ⚪ not loaded',
   /data-th-layer="cusum_theta"[\s\S]{0,150}⚪ not loaded/.test(html));

// ============================================================================
// 3. All 9 panels present with the documented `th`-prefixed ids
// ============================================================================
console.log('\n=== 3. All 9 mirrored panels present ===');

const expectedPanelIds = [
  'thCtrlBar',
  'thCusumHeroPanel',
  'thSimPanel',
  'thZPanel',
  'thLinesPanel',
  'thAnchorStripPanel',
  'thPcaPanel',
  'thTrackedSamplesPanelCompact',
  'thL3Panel',
];

for (const panelId of expectedPanelIds) {
  const re = new RegExp('id="' + panelId + '"');
  ok('panel #' + panelId + ' exists', re.test(html));
}

// ============================================================================
// 4. All 9 panels in correct top-to-bottom order
// ============================================================================
console.log('\n=== 4. Panel ordering ===');

const positions = expectedPanelIds.map(id => ({
  id,
  pos: html.indexOf('id="' + id + '"'),
}));

ok('all panels found (no -1 positions)',
   positions.every(p => p.pos > -1),
   positions.filter(p => p.pos === -1).map(p => p.id).join(','));

// Check strict ascending order
let inOrder = true;
let orderError = null;
for (let i = 1; i < positions.length; i++) {
  if (positions[i].pos <= positions[i - 1].pos) {
    inOrder = false;
    orderError = positions[i - 1].id + ' should precede ' + positions[i].id;
    break;
  }
}
ok('panels appear in documented top-to-bottom order', inOrder, orderError);

// ============================================================================
// 5. All 9 panels nested inside #page12 (not leaked into #page15 or sibling)
// ============================================================================
console.log('\n=== 5. Panels nested inside #page12 ===');

const page12Start = html.indexOf('<div id="page12"');
const page15Start = html.indexOf('<div id="page15"');
ok('#page12 found', page12Start > -1);
ok('#page15 found (sibling)', page15Start > -1);
ok('#page12 precedes #page15', page12Start < page15Start);

for (const panelId of expectedPanelIds) {
  const panelPos = html.indexOf('id="' + panelId + '"');
  ok('#' + panelId + ' is between #page12 and #page15',
     panelPos > page12Start && panelPos < page15Start,
     'pos=' + panelPos + ' page12=' + page12Start + ' page15=' + page15Start);
}

// ============================================================================
// 6. All 9 panels start display: none
// ============================================================================
console.log('\n=== 6. All panels start display: none ===');

for (const panelId of expectedPanelIds) {
  // Match the opening tag for the panel id and inspect its style attribute.
  const re = new RegExp('id="' + panelId + '"[^>]*style="[^"]*display:\\s*none', 'i');
  ok('#' + panelId + ' starts display: none', re.test(html));
}

// ============================================================================
// 7. Hero panel (#thCusumHeroPanel) has the documented inner structure
// ============================================================================
console.log('\n=== 7. Hero panel inner structure ===');

ok('#thCusumHeroPanel contains #thCusumStripCanvas (lane 1: per-carrier ticks)',
   /id="thCusumHeroPanel"[\s\S]{0,1500}id="thCusumStripCanvas"/.test(html));

ok('#thCusumHeroPanel contains #thCusumHistCanvas (lane 2: KDE histogram)',
   /id="thCusumHeroPanel"[\s\S]{0,1500}id="thCusumHistCanvas"/.test(html));

ok('#thCusumHeroPanel contains #thCusumSpreadBadge (tight/intermediate/ragged)',
   /id="thCusumHeroPanel"[\s\S]{0,1800}id="thCusumSpreadBadge"/.test(html));

ok('#thCusumHeroPanel contains #thCusumConcordBadge (cross-stream concord)',
   /id="thCusumHeroPanel"[\s\S]{0,1800}id="thCusumConcordBadge"/.test(html));

ok('hero panel documents page 1 has no equivalent',
   /<code>#thCusumHeroPanel<\/code>[\s\S]{0,200}page 1 has none/.test(html));

// ============================================================================
// 8. Required canvas/container children for non-hero panels
// ============================================================================
console.log('\n=== 8. Inner elements for mirrored panels ===');

const innerElements = [
  ['thSimPanel', 'thSimCanvas'],
  ['thZPanel', 'thZCanvas'],
  ['thLinesPanel', 'thLinesCanvasContainer'],
  ['thAnchorStripPanel', 'thAnchorStripCanvas'],
  ['thPcaPanel', 'thPcaCanvas'],
  ['thTrackedSamplesPanelCompact', 'thTrackedSamplesContainer'],
  ['thL3Panel', 'thL3Container'],
];

for (const [panelId, innerId] of innerElements) {
  const re = new RegExp(
    'id="' + panelId + '"[\\s\\S]{0,1200}id="' + innerId + '"'
  );
  ok('#' + panelId + ' contains #' + innerId, re.test(html));
}

// ============================================================================
// 9. Ctrl bar mirrors page 1 controls
// ============================================================================
console.log('\n=== 9. #thCtrlBar mirrors #ctrlBar controls ===');

ok('#thCtrlBar contains #thPlayBtn',
   /id="thCtrlBar"[\s\S]{0,1500}id="thPlayBtn"/.test(html));

ok('#thCtrlBar contains #thScrubber',
   /id="thCtrlBar"[\s\S]{0,1500}id="thScrubber"/.test(html));

ok('#thCtrlBar contains #thViewModeBar with genome/L1/L2 buttons',
   /id="thViewModeBar"[\s\S]{0,800}data-thviewmode="genome"[\s\S]{0,400}data-thviewmode="l1"[\s\S]{0,400}data-thviewmode="l2"/.test(html));

ok('#thCtrlBar contains #thWinIdx label',
   /id="thCtrlBar"[\s\S]{0,1500}id="thWinIdx"/.test(html));

// ============================================================================
// 10. Panel-layout description grid was extended to include all 9 panels
// ============================================================================
console.log('\n=== 10. Panel-layout description grid extended ===');

const descPanels = [
  '#thCtrlBar', '#thCusumHeroPanel', '#thSimPanel', '#thZPanel',
  '#thLinesPanel', '#thAnchorStripPanel', '#thPcaPanel',
  '#thTrackedSamplesPanelCompact', '#thL3Panel',
];

for (const dp of descPanels) {
  const re = new RegExp('<code>' + dp.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '<\\/code>');
  ok('panel-layout grid documents ' + dp, re.test(html));
}

// ============================================================================
// 11. Panels do NOT leak into #page15 (GHSL page)
// ============================================================================
console.log('\n=== 11. No leakage into #page15 ===');

const page15ToEnd = html.slice(page15Start);
for (const panelId of expectedPanelIds) {
  ok('#' + panelId + ' does NOT appear inside #page15',
     !new RegExp('id="' + panelId + '"').test(page15ToEnd));
}

// ============================================================================
// 12. No accidental JS render functions leaked in (pure scaffold slice)
// ============================================================================
console.log('\n=== 12. Slice is scaffold-only — no renderers yet ===');

// We deliberately do NOT define drawing functions in this slice. If any
// of these names appear, that means a future slice's code crept in early.
const tooEarlyFunctions = [
  'function drawThSimMat',
  'function drawThZ',
  'function drawThCusumHero',
  'function drawThPca',
  'function _wireThetaPiPanels',
];

for (const fn of tooEarlyFunctions) {
  ok('no premature ' + fn,
     html.indexOf(fn) === -1,
     'unexpected: ' + fn + ' should ship in later slice');
}

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
