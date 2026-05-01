// Tests for v4 turns 91, 93, 94, 95 — UI hooks for detailed mode toggle,
// catalogue diamond column, diamond zone overlay, SNP density strip.
//
// Coverage:
//   1. Helpers exposed
//   2. Turn 91: activeMode toggle UI present + responds to clicks
//   3. Turn 91: detailed mode banner toggles visibility correctly
//   4. Turn 93: catalogue Diamond column wired with strictness toggle
//   5. Turn 93: 3-button strictness toolbar present + clicks switch mode
//   6. Turn 94: _drawDiamondOverlay function present + safe under empty state
//   7. Turn 95: SNP density toggle UI + state.linesSnpDensityOn persists
//   8. Turn 95: _snpDensityForWindow falls back through priorities

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[ui-mt] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[ui-mt] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
const required = [
  // Turn 91 helpers
  '_activeModeRefreshBanner', '_activeModeApplyToButtons',
  // Turn 93 helpers
  '_catComputeRowDiamondCount', '_catComputeRowDiamondSummary',
  '_catBuildTransientCandidateFromL2', 'setCatDiamondMode',
  // Turn 94
  '_drawDiamondOverlay',
  // Turn 95
  '_drawSnpDensityStrip', '_snpDensityForWindow', '_snpDensitySource',
  'setLinesSnpDensityOn', 'setLinesSnpDensityMode',
];
for (const name of required) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed (or wrong type)');
}
ok('all ' + required.length + ' helpers exposed (turns 91 + 93 + 94 + 95)');

// ---- 2. Turn 91: activeMode toggle UI ----
const modeBar = win.document.getElementById('activeModeBar');
if (!modeBar) fail('activeModeBar div missing');
const defaultBtn = modeBar.querySelector('button[data-mode="default"]');
const detailedBtn = modeBar.querySelector('button[data-mode="detailed"]');
if (!defaultBtn) fail('default mode button missing');
if (!detailedBtn) fail('detailed mode button missing');
ok('activeMode toggle: bar present with default + detailed buttons');

// Initially default should be active
if (!defaultBtn.classList.contains('active')) {
  fail('default button should be active on load (got class list: ' + defaultBtn.className + ')');
}
ok('initial state: default button is active');

// Click detailed → should switch
detailedBtn.click();
if (win.getActiveMode() !== 'detailed') fail('clicking detailed should switch mode to detailed');
if (!detailedBtn.classList.contains('active')) fail('detailed button should now be active');
if (defaultBtn.classList.contains('active')) fail('default button should no longer be active');
ok('clicking detailed: state.activeMode → "detailed", buttons toggle correctly');

// Banner becomes visible in detailed mode
const banner = win.document.getElementById('detailedModeBanner');
if (!banner) fail('detailedModeBanner missing');
if (banner.style.display !== 'block') fail('banner should be visible in detailed mode (display=block), got ' + banner.style.display);
ok('detailedModeBanner: visible in detailed mode');

// Click default → banner hides
defaultBtn.click();
if (win.getActiveMode() !== 'default') fail('clicking default should switch back');
if (banner.style.display !== 'none') fail('banner should hide in default mode');
ok('clicking default: banner hides');

// ---- 3. Turn 93: catalogue diamond column ----
// CAT_COLS_BASE should include a 'diamond' column
const html2 = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
if (!html2.includes("key: 'diamond'")) fail("CAT_COLS_BASE should have key: 'diamond'");
if (!html2.includes('Diamond strictness')) fail("toolbar should have diamond strictness controls");
ok('catalogue Diamond column registered in CAT_COLS_BASE');

// 3-button strictness toolbar present
const dLoose = win.document.getElementById('catDiamondLoose');
const dStrict = win.document.getElementById('catDiamondStrict');
const dStrict2 = win.document.getElementById('catDiamondStrict2');
if (!dLoose) fail('catDiamondLoose button missing');
if (!dStrict) fail('catDiamondStrict button missing');
if (!dStrict2) fail('catDiamondStrict2 button missing');
ok('3-button diamond strictness toolbar: loose + strict + strict-2 all present');

// Initial active = loose
if (!dLoose.classList.contains('active')) fail('loose should be initial active');
ok('initial diamond mode: loose (active class on loose button)');

// Click strict → switches
dStrict.click();
if (!dStrict.classList.contains('active')) fail('strict should be active after click');
if (dLoose.classList.contains('active')) fail('loose should no longer be active');
ok('clicking strict button: switches active class correctly');

// Click strict-2
dStrict2.click();
if (!dStrict2.classList.contains('active')) fail('strict-2 should be active');
if (dStrict.classList.contains('active')) fail('strict should no longer be active');
ok('clicking strict-2 button: active class moves correctly');

// localStorage persistence
let storedMode = null;
try { storedMode = win.localStorage.getItem('pca_scrubber_v3.catDiamondMode'); } catch (_) {}
if (storedMode !== 'strict2') fail('expected strict2 in localStorage, got ' + storedMode);
ok('catDiamondMode persists in localStorage: strict2');

// _catBuildTransientCandidateFromL2 returns proper candidate shape
const env = { start_w: 0, end_w: 5, start_bp: 1000, end_bp: 5000 };
const cl  = { labels: new Int8Array([0, 1, 2, 0, 1, 2, 0, 1]), usedK: 3 };
const transient = win._catBuildTransientCandidateFromL2(0, env, cl);
if (!transient) fail('_catBuildTransientCandidateFromL2 returned null');
if (transient.K !== 3) fail('transient.K should be 3, got ' + transient.K);
if (transient.start_w !== 0 || transient.end_w !== 5) fail('start_w/end_w not propagated');
if (!transient.locked_labels) fail('locked_labels not propagated');
ok('_catBuildTransientCandidateFromL2: builds candidate with K, start_w, end_w, locked_labels');

// Reset to loose for downstream tests
dLoose.click();

// ---- 4. Turn 94: diamond overlay function ----
// Should be safe to call with no candidate — early returns
win.state.candidate = null;
let didThrow = false;
try {
  // Build a fake ctx that records calls
  const fakeCtx = {
    save: () => {},
    restore: () => {},
    fillRect: () => {},
    strokeRect: () => {},
    beginPath: () => {},
    moveTo: () => {},
    lineTo: () => {},
    stroke: () => {},
    fillText: () => {},
    setLineDash: () => {},
    fillStyle: '', strokeStyle: '', lineWidth: 1, font: '', textAlign: '', textBaseline: '',
  };
  const pad = { l: 44, r: 16, t: 6, b: 8 };
  win._drawDiamondOverlay(fakeCtx, pad, 600, 200, 0, 100, 660, 220);
} catch (e) {
  didThrow = true;
}
if (didThrow) fail('_drawDiamondOverlay should not throw on empty state');
ok('_drawDiamondOverlay: safe to call with no active candidate (no throw)');

// ---- 5. Turn 95+99: SNP density 3-button toggle ----
const snpBar = win.document.getElementById('linesSnpDensityBar');
if (!snpBar) fail('linesSnpDensityBar missing');
const snpOff   = snpBar.querySelector('button[data-snpdens-mode="off"]');
const snpStrip = snpBar.querySelector('button[data-snpdens-mode="strip"]');
const snpShade = snpBar.querySelector('button[data-snpdens-mode="shade"]');
if (!snpOff || !snpStrip || !snpShade) fail('SNP density 3-button toggle missing');
ok('SNP density 3-button toggle: off + strip + shade all present in DOM');

// State default
if (win.state.linesSnpDensityMode && win.state.linesSnpDensityMode !== 'off') {
  fail('state.linesSnpDensityMode should default to "off", got ' + win.state.linesSnpDensityMode);
}
ok('state.linesSnpDensityMode: defaults to off');

// Click strip
snpStrip.click();
if (win.state.linesSnpDensityMode !== 'strip') fail('clicking strip should set mode=strip');
if (!snpStrip.classList.contains('active')) fail('strip button should be active');
let storedSnpMode = null;
try { storedSnpMode = win.localStorage.getItem('inversion_atlas.linesSnpDensityMode'); } catch (_) {}
if (storedSnpMode !== 'strip') fail('strip mode should persist, got ' + storedSnpMode);
// Legacy boolean state synced
if (!win.state.linesSnpDensityOn) fail('legacy linesSnpDensityOn should be true when mode=strip');
ok('clicking strip: state.linesSnpDensityMode → "strip", localStorage persists, legacy bool synced');

// Click shade
snpShade.click();
if (win.state.linesSnpDensityMode !== 'shade') fail('clicking shade should set mode=shade');
if (!snpShade.classList.contains('active')) fail('shade button should be active');
if (snpStrip.classList.contains('active')) fail('strip should no longer be active');
ok('clicking shade: state.linesSnpDensityMode → "shade", active class moves');

// Click off
snpOff.click();
if (win.state.linesSnpDensityMode !== 'off') fail('clicking off should set mode=off');
if (win.state.linesSnpDensityOn) fail('legacy linesSnpDensityOn should be false when mode=off');
ok('clicking off: state.linesSnpDensityMode → "off", legacy bool synced false');

// _drawSnpDensityShade present + safe
if (typeof win._drawSnpDensityShade !== 'function') fail('_drawSnpDensityShade not exposed');
let didThrowShade = false;
try {
  const fakeCtx = {
    save: () => {}, restore: () => {}, fillRect: () => {},
    fillStyle: '',
  };
  const pad = { l: 44, r: 16, t: 6, b: 8 };
  win._drawSnpDensityShade(fakeCtx, pad, 600, 200, 0, 10);
} catch (_) { didThrowShade = true; }
if (didThrowShade) fail('_drawSnpDensityShade should be safe when off');
ok('_drawSnpDensityShade: safe when off (early-return)');

// ---- Turn 102: trans-rate strip toggle UI ----
const transCb = win.document.getElementById('linesTransRateToggle');
if (!transCb) fail('linesTransRateToggle checkbox missing');
ok('trans-rate toggle: checkbox present in DOM');

if (typeof win.setLinesTransRateOn !== 'function') fail('setLinesTransRateOn not exposed');
if (typeof win._drawTransitionRateStrip !== 'function') fail('_drawTransitionRateStrip not exposed');
ok('setLinesTransRateOn + _drawTransitionRateStrip exposed');

// Default off
if (win.state.linesTransRateOn) fail('state.linesTransRateOn should default to false');
ok('state.linesTransRateOn defaults to false');

// Toggle on
win.setLinesTransRateOn(true);
if (!win.state.linesTransRateOn) fail('setLinesTransRateOn(true) failed');
let storedTrans = null;
try { storedTrans = win.localStorage.getItem('inversion_atlas.linesTransRateOn'); } catch (_) {}
if (storedTrans !== '1') fail('trans-rate state should persist');
ok('setLinesTransRateOn(true): state + localStorage persist');

// Drawer safe when off
win.setLinesTransRateOn(false);
let didThrowTrans = false;
try {
  const fakeCtx = { save: () => {}, restore: () => {}, fillRect: () => {}, strokeRect: () => {}, beginPath: () => {}, moveTo: () => {}, lineTo: () => {}, stroke: () => {}, fillStyle: '', strokeStyle: '', lineWidth: 0 };
  win._drawTransitionRateStrip(fakeCtx, { l: 44, r: 16, t: 6, b: 8 }, 600, 200, 0, 10);
} catch (_) { didThrowTrans = true; }
if (didThrowTrans) fail('_drawTransitionRateStrip should be safe when off');
ok('_drawTransitionRateStrip: safe when off');

// ---- 6. Turn 95: _snpDensityForWindow priority ----
// Set up a fake window with each priority field
win.state.data = {
  windows: [
    { center_mb: 1.0, n_snps: 42 },                      // priority 1
    { center_mb: 2.0, snp_count: 35 },                    // priority 1 (alt name)
    { center_mb: 3.0, lambda1: 12.5 },                    // priority 2
    { center_mb: 4.0, variance_pc1: 8.7 },                // priority 2 (alt)
    { center_mb: 5.0, pc1: [0.1, 0.5, -0.2, 0.3, -0.1] }, // priority 3 (proxy)
    { center_mb: 6.0 },                                    // none
  ],
};

const v0 = win._snpDensityForWindow(0);
if (v0 !== 42) fail('window 0 should return n_snps=42, got ' + v0);
ok('_snpDensityForWindow[0]: n_snps=42 (priority 1: precomp:n_snps)');

const v1 = win._snpDensityForWindow(1);
if (v1 !== 35) fail('window 1 should return snp_count=35, got ' + v1);
ok('_snpDensityForWindow[1]: snp_count=35 (priority 1 alt)');

const v2 = win._snpDensityForWindow(2);
if (v2 !== 12.5) fail('window 2 should return lambda1=12.5, got ' + v2);
ok('_snpDensityForWindow[2]: lambda1=12.5 (priority 2: precomp lambda1)');

const v3 = win._snpDensityForWindow(3);
if (v3 !== 8.7) fail('window 3 should return variance_pc1=8.7, got ' + v3);
ok('_snpDensityForWindow[3]: variance_pc1=8.7 (priority 2 alt)');

const v4 = win._snpDensityForWindow(4);
if (v4 == null || !isFinite(v4)) fail('window 4 should return pc1-variance proxy, got ' + v4);
ok('_snpDensityForWindow[4]: pc1 variance proxy = ' + v4.toFixed(4));

const v5 = win._snpDensityForWindow(5);
if (v5 !== null) fail('window 5 should return null (no source available)');
ok('_snpDensityForWindow[5]: null (no source available)');

// Source labels match priority
if (win._snpDensitySource(0) !== 'precomp:n_snps') fail('source[0] mismatch');
if (win._snpDensitySource(2) !== 'precomp:lambda1') fail('source[2] mismatch');
if (win._snpDensitySource(4) !== 'proxy:pc1_variance') fail('source[4] mismatch');
if (win._snpDensitySource(5) !== 'none') fail('source[5] should be "none"');
ok('_snpDensitySource: priority labels (precomp:n_snps / precomp:lambda1 / proxy:pc1_variance / none) all correct');

// _drawSnpDensityStrip safe to call when off (turn 99: now mode-based)
let didThrowStrip = false;
try {
  win.state.linesSnpDensityMode = 'off';
  win.state.linesSnpDensityOn = false;
  const fakeCtx = {
    save: () => {}, restore: () => {}, fillRect: () => {}, strokeRect: () => {},
    fillStyle: '', strokeStyle: '', lineWidth: 0,
  };
  const pad = { l: 44, r: 16, t: 6, b: 8 };
  win._drawSnpDensityStrip(fakeCtx, pad, 600, 200, 0, 10);
} catch (_) { didThrowStrip = true; }
if (didThrowStrip) fail('_drawSnpDensityStrip should be safe when mode=off');
ok('_drawSnpDensityStrip: safe when mode=off (early-return)');

console.log('\n[ui-mt] ALL CHECKS PASSED');
