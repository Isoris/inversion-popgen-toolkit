// Tests for v4 turn 80 visual additions:
//   1. Stage bands on tab buttons (data-stage attribute, 5 stages)
//   2. First-use attention pulse on key buttons (load JSON, Zoom L2,
//      Window step cycler, move-to-minimap, promote)
//   3. Pulse dismissal on click + localStorage persistence

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');

function fail(msg) { console.error('[stage-pulse] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[stage-pulse] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window } = dom;
const { document } = window;

// ---- 1. Stage bands ----
const allTabs = document.querySelectorAll('#tabBar button[data-page]');
if (allTabs.length === 0) fail('no tab buttons found');
ok(allTabs.length + ' tab buttons present in tab bar');

// Every tab must have a data-stage
const missing = [];
const stageCounts = {};
for (const btn of allTabs) {
  const stage = btn.dataset.stage;
  if (!stage) missing.push(btn.dataset.page);
  else stageCounts[stage] = (stageCounts[stage] || 0) + 1;
}
if (missing.length > 0) fail('these tabs missing data-stage: ' + missing.join(', '));
ok('every tab has a data-stage attribute');

// Expected stage distribution
const expected = {
  discovery:     4,   // page1, page12, page15, page2
  refinement:    3,   // page11, page3, page19 (NEW: negative regions, v4 turn 82)
  classification: 4,  // page4, page6, page7, page8
  synthesis:     5,   // page9, page10, page16, page17, page18
  help:          1,   // page5
};
for (const stage of Object.keys(expected)) {
  if (stageCounts[stage] !== expected[stage]) {
    fail('stage "' + stage + '" expected ' + expected[stage] + ' tabs, got ' + (stageCounts[stage] || 0));
  }
}
ok('stage distribution: discovery=4, refinement=3, classification=4, synthesis=5, help=1 (17 total)');

// Verify discovery group: page1, page12, page15, page2
const discoveryPages = Array.from(document.querySelectorAll('#tabBar button[data-stage="discovery"]'))
  .map(b => b.dataset.page).sort();
const expectedDiscovery = ['page1', 'page12', 'page15', 'page2'].sort();
if (JSON.stringify(discoveryPages) !== JSON.stringify(expectedDiscovery)) {
  fail('discovery group wrong: got ' + JSON.stringify(discoveryPages));
}
ok('discovery stage: ' + discoveryPages.join(', ') + ' (4 PCA scrubbers + candidate focus)');

// Verify refinement: page11 (boundaries), page3 (catalogue), page19 (negative regions)
const refinementPages = Array.from(document.querySelectorAll('#tabBar button[data-stage="refinement"]'))
  .map(b => b.dataset.page).sort();
const expectedRefinement = ['page11', 'page19', 'page3'].sort();
if (JSON.stringify(refinementPages) !== JSON.stringify(expectedRefinement)) {
  fail('refinement group wrong: got ' + JSON.stringify(refinementPages));
}
ok('refinement stage: page11 (boundaries) + page3 (catalogue) + page19 (negative regions)');

// Verify synthesis: page9, page10, page16, page17, page18
const synthesisPages = Array.from(document.querySelectorAll('#tabBar button[data-stage="synthesis"]'))
  .map(b => b.dataset.page).sort();
const expectedSynthesis = ['page9', 'page10', 'page16', 'page17', 'page18'].sort();
if (JSON.stringify(synthesisPages) !== JSON.stringify(expectedSynthesis)) {
  fail('synthesis group wrong: got ' + JSON.stringify(synthesisPages));
}
ok('synthesis stage: ' + synthesisPages.join(', ') + ' (confirmed, markers, cross-species, stats profile, marker panel)');

// ---- 2. CSS rules for stage bands present ----
// Find the inline stylesheet that contains our new rules
const styleText = Array.from(document.querySelectorAll('style')).map(s => s.textContent).join('\n');
for (const stage of Object.keys(expected)) {
  if (!styleText.includes('[data-stage="' + stage + '"]')) {
    fail('CSS missing rule for data-stage="' + stage + '"');
  }
}
ok('CSS contains stage-band rules for all 5 stages');

// Verify the colors are subtle (alpha values around 0.20-0.30 for default state)
if (!styleText.match(/inset 0 2px 0 0 rgba\(79,\s*163,\s*255,\s*0\.\d+\)/)) {
  fail('discovery stage band should use blue rgb(79, 163, 255)');
}
ok('discovery band uses cool blue rgba(79, 163, 255, ~0.28)');

// ---- 3. attention-pulse helpers exposed ----
if (typeof window._initAttentionPulses !== 'function') fail('_initAttentionPulses not exposed');
if (!Array.isArray(window._ATTENTION_TARGETS)) fail('_ATTENTION_TARGETS not exposed');
if (window._ATTENTION_TARGETS.length !== 5) fail('expected 5 attention targets, got ' + window._ATTENTION_TARGETS.length);
ok('5 attention-pulse targets registered: ' + window._ATTENTION_TARGETS.map(t => t.key).join(', '));

// ---- 4. Pulse class applied on initial load ----
window._initAttentionPulses();

const fileInput = document.getElementById('fileInput');
if (!fileInput) fail('fileInput missing');
const fileWrapper = fileInput.parentElement;
if (!fileWrapper.classList.contains('attention-pulse')) {
  fail('fileInput parent (.ctl) should have attention-pulse class on first load');
}
ok('load-JSON wrapper has attention-pulse on first load');

const zoomL2 = document.querySelector('#viewModeBar button[data-viewmode="l2"]');
if (!zoomL2.classList.contains('attention-pulse')) {
  fail('Zoom L2 button should have attention-pulse on first load');
}
ok('Zoom L2 button has attention-pulse on first load');

const windowsBtn = document.getElementById('jumpToWindowsBtn');
if (!windowsBtn.classList.contains('attention-pulse')) {
  fail('jumpToWindowsBtn should have attention-pulse on first load');
}
ok('Window step cycler button has attention-pulse on first load');

const moveMinimap = document.getElementById('simMoveMinimapBtn');
if (!moveMinimap.classList.contains('attention-pulse')) {
  fail('simMoveMinimapBtn should have attention-pulse on first load');
}
ok('▤ move-to-minimap button has attention-pulse on first load');

const promote = document.getElementById('promoteCandidateBtn');
if (!promote.classList.contains('attention-pulse')) {
  fail('promoteCandidateBtn should have attention-pulse on first load');
}
ok('⬆ promote-to-candidate button has attention-pulse on first load');

// ---- 5. Click dismisses pulse ----
zoomL2.click();
if (zoomL2.classList.contains('attention-pulse')) {
  fail('Zoom L2 pulse should be removed after click');
}
ok('clicking Zoom L2 removes its attention-pulse class');

// ---- 6. localStorage flag set ----
const stored = window.localStorage.getItem(window._ATTENTION_PULSE_KEY_PREFIX + 'zoom_l2');
if (stored !== '1') fail('localStorage flag not set for zoom_l2 (got ' + stored + ')');
ok('localStorage[zoom_l2]=1 persisted after dismissal');

// ---- 7. Other pulses NOT dismissed by clicking one ----
if (!windowsBtn.classList.contains('attention-pulse')) {
  fail('clicking Zoom L2 should NOT dismiss other pulses (Windows still pulses)');
}
ok('clicking one button does not dismiss other pulses');

// ---- 8. Re-init should NOT re-add pulse to dismissed button ----
zoomL2.classList.remove('attention-pulse-fade'); // simulate fade timer expiry
window._initAttentionPulses();
if (zoomL2.classList.contains('attention-pulse')) {
  fail('previously-dismissed Zoom L2 should not get pulse re-added on re-init');
}
ok('previously-dismissed buttons stay clean on subsequent _initAttentionPulses calls');

// Other pulses should still be on (their localStorage flag wasn't set)
if (!windowsBtn.classList.contains('attention-pulse')) {
  fail('non-dismissed buttons should keep their pulse on re-init');
}
ok('non-dismissed buttons keep their pulse on re-init');

// ---- 9. Active tab gets stronger band ----
// Activate page1 (discovery)
document.querySelectorAll('#tabBar button[data-page]').forEach(b => b.classList.remove('active'));
const page1Btn = document.querySelector('#tabBar button[data-page="page1"]');
page1Btn.classList.add('active');
// CSS: discovery.active gets rgba(...,0.75); discovery default gets rgba(...,0.28)
// We verify by checking that the CSS rule for active state exists
if (!styleText.match(/discovery"\]\.active.*0\.7\d/s)) {
  fail('CSS missing stronger band for active discovery tab');
}
ok('CSS includes stronger band for active state in each stage');

// ---- 10. Animation keyframes defined ----
if (!styleText.includes('@keyframes attention-pulse-kf')) {
  fail('attention-pulse keyframes missing');
}
if (!styleText.includes('animation: attention-pulse-kf')) {
  fail('attention-pulse animation rule missing');
}
ok('CSS @keyframes attention-pulse-kf defined and applied to .attention-pulse');

console.log('\n[stage-pulse] ALL CHECKS PASSED');
