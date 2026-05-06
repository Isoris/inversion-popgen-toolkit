// Tests for v4 turn 93 — Catalogue diamond column.
// Coverage:
//   1. Helpers exposed
//   2. Diamond column added to CAT_COLS_BASE
//   3. Toolbar button present + cycles loose → strict → strict2
//   4. Diamond mode persists in localStorage
//   5. _catComputeRowDiamondCount returns 0 when no detector input
//   6. _catComputeRowDiamondCount returns count from summary
//   7. _catBuildTransientCandidateFromL2 creates valid candidate shape
//   8. catState.diamondMode default and switching

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[catd] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[catd] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['_catComputeRowDiamondCount', '_catComputeRowDiamondSummary',
                    '_catBuildTransientCandidateFromL2', 'setCatDiamondMode']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed on window');
}
if (!Array.isArray(win._CAT_DIAMOND_MODES)) fail('_CAT_DIAMOND_MODES not exposed');
if (win._CAT_DIAMOND_MODES.length !== 3) fail('expected 3 diamond modes');
ok('all 4 catalogue-diamond helpers + 1 modes constant exposed on window');

// ---- 2. Diamond column added ----
const cols = win.getCatCols();
const diamondCol = cols.find(c => c.key === 'diamond');
if (!diamondCol) fail('diamond column missing from catalogue columns');
if (diamondCol.type !== 'num') fail('diamond column type should be "num"');
if (typeof diamondCol.fmtHtml !== 'function') fail('diamond column should have fmtHtml');
if (!diamondCol.tooltip) fail('diamond column should have tooltip');
ok('diamond column present: type=num, has fmtHtml + tooltip');

// fmt for empty
if (diamondCol.fmt(0) !== '—') fail('count=0 should format as "—"');
if (diamondCol.fmt(NaN) !== '—') fail('NaN should format as "—"');
if (diamondCol.fmt(3) !== '3') fail('count=3 should format as "3"');
ok('diamond column fmt: 0/NaN → "—", positive → number string');

// fmtHtml for non-empty contains diamond glyph
const html3 = diamondCol.fmtHtml(3);
if (!html3.includes('◆')) fail('fmtHtml should include diamond glyph');
if (!html3.includes('3')) fail('fmtHtml should include count');
ok('diamond column fmtHtml: glyph ◆ + count rendered');

// ---- 3. Toolbar: 3 buttons present ----
const btnLoose = win.document.getElementById('catDiamondLoose');
const btnStrict = win.document.getElementById('catDiamondStrict');
const btnStrict2 = win.document.getElementById('catDiamondStrict2');
if (!btnLoose) fail('catDiamondLoose button missing from DOM');
if (!btnStrict) fail('catDiamondStrict button missing');
if (!btnStrict2) fail('catDiamondStrict2 button missing');
ok('catalogue toolbar: 3 diamond buttons (loose / strict / strict-2) present');

// Initial active class on loose
if (!btnLoose.classList.contains('active')) fail('initial: loose button should have active class');
ok('initial active class on "loose" button');

// Click strict button
btnStrict.click();
if (win.catState.diamondMode !== 'strict') fail('after click strict, mode should be "strict", got "' + win.catState.diamondMode + '"');
if (!btnStrict.classList.contains('active')) fail('strict button should have active class');
if (btnLoose.classList.contains('active')) fail('loose button should no longer be active');
ok('click strict button → mode="strict", active class moved');

// Click strict2 button
btnStrict2.click();
if (win.catState.diamondMode !== 'strict2') fail('after click strict2, mode should be "strict2"');
if (!btnStrict2.classList.contains('active')) fail('strict2 button should have active class');
ok('click strict-2 button → mode="strict2"');

// Click back to loose
btnLoose.click();
if (win.catState.diamondMode !== 'loose') fail('after click loose, mode should be "loose"');
if (!btnLoose.classList.contains('active')) fail('loose button should be active again');
ok('click loose button → mode="loose" (full cycle)');

// ---- 4. Mode persistence ----
win.setCatDiamondMode('strict2');
let stored = null;
try { stored = win.localStorage.getItem('pca_scrubber_v3.catDiamondMode'); } catch (_) {}
if (stored !== 'strict2') fail('strict2 mode should persist in localStorage');
ok('mode persistence: strict2 stored in localStorage');

// Bogus mode falls back to loose
win.setCatDiamondMode('bogus');
if (win.catState.diamondMode !== 'loose') fail('bogus mode should fall back to "loose"');
ok('bogus mode value: falls back to "loose"');

// ---- 5. _catComputeRowDiamondCount with no data ----
// Reset state
win.state.data = null;
const noDataCount = win._catComputeRowDiamondCount(0, {}, { labels: null });
if (noDataCount !== 0) fail('no detector input should give 0, got ' + noDataCount);
ok('count with no detector input: 0 (graceful fallback)');

// ---- 6. _catComputeRowDiamondCount with synthetic diamond data ----
// Build a synthetic state.data that triggers a strict2 diamond
const nS = 30;
const labels = new Int8Array(nS);
for (let si = 0; si < 10; si++) labels[si] = 0;
for (let si = 10; si < 20; si++) labels[si] = 1;
for (let si = 20; si < 30; si++) labels[si] = 2;

const wins = [];
for (let wi = 0; wi < 20; wi++) {
  const pc1 = new Array(nS);
  for (let si = 0; si < nS; si++) {
    const band = labels[si];
    if (band === 0) pc1[si] = -2 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
    else if (band === 2) pc1[si] = 2 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
    else {
      // Band 1 splits in middle windows
      const subBand = (si - 10) < 5 ? -1 : 1;
      if (wi >= 6 && wi <= 13) pc1[si] = 0 + subBand * 0.6;
      else pc1[si] = 0 + (((si * 31) % 100) / 100 - 0.5) * 0.05;
    }
  }
  wins.push({ pc1, pc2: new Array(nS).fill(0) });
}
win.state.data = { n_samples: nS, windows: wins };

const env = { start_w: 0, end_w: 19, start_bp: 0, end_bp: 1000000, candidate_id: 'test_l2' };
const cl = { labels, usedK: 3 };

const summary = win._catComputeRowDiamondSummary(0, env, cl);
if (!summary) fail('summary should not be null with synthetic diamond data');
if (summary.n_loose < 1) fail('synthetic diamond should give loose count >= 1, got ' + summary.n_loose);
if (summary.n_strict2 < 1) fail('synthetic diamond should give strict2 count >= 1');
ok('summary: loose=' + summary.n_loose + ', strict=' + summary.n_strict + ', strict2=' + summary.n_strict2);

// Set mode and verify count flips
win.setCatDiamondMode('loose');
const cLoose = win._catComputeRowDiamondCount(0, env, cl);
if (cLoose !== summary.n_loose) fail('loose count should match summary.n_loose');
ok('catDiamondMode=loose: count = ' + cLoose);

win.setCatDiamondMode('strict');
const cStrict = win._catComputeRowDiamondCount(0, env, cl);
if (cStrict !== summary.n_strict) fail('strict count should match summary.n_strict');
ok('catDiamondMode=strict: count = ' + cStrict);

win.setCatDiamondMode('strict2');
const cStrict2 = win._catComputeRowDiamondCount(0, env, cl);
if (cStrict2 !== summary.n_strict2) fail('strict2 count should match summary.n_strict2');
ok('catDiamondMode=strict2: count = ' + cStrict2);

// ---- 7. _catBuildTransientCandidateFromL2 ----
const cand = win._catBuildTransientCandidateFromL2(0, env, cl);
if (!cand) fail('transient candidate should be built');
if (cand.K !== 3) fail('candidate K should be 3');
if (!cand.locked_labels) fail('candidate should have locked_labels');
if (cand.start_w !== 0) fail('candidate start_w should be 0');
if (cand.end_w !== 19) fail('candidate end_w should be 19');
ok('_catBuildTransientCandidateFromL2: valid shape (K=3, start_w=0, end_w=19)');

// Without labels → null
const noLabels = win._catBuildTransientCandidateFromL2(0, env, { labels: null });
if (noLabels !== null) fail('no labels should give null');
ok('_catBuildTransientCandidateFromL2: null when no labels');

// ---- 8. catState.diamondMode field exists with correct default ----
// (After our cycling above it's now strict2; reset and check default behavior)
delete win.catState.diamondMode;
try { win.localStorage.removeItem('pca_scrubber_v3.catDiamondMode'); } catch (_) {}
// Re-evaluate the init — since the script already ran, just check the field
// is present after setMode
win.setCatDiamondMode('loose');
if (win.catState.diamondMode !== 'loose') fail('catState.diamondMode should be "loose" after init');
ok('catState.diamondMode field: present, default-able to "loose"');

console.log('\n[catd] ALL CHECKS PASSED');
