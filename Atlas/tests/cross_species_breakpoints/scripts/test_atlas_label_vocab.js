// Tests for v4 turn 87 — karyotype label vocabulary toggle.
// Coverage:
//   1. Helpers exposed
//   2. Default vocab is 'legacy' on first run
//   3. Legacy labels for K=3 (band 1 (lo) / band 2 (mid) / band 3 (hi))
//   4. Legacy labels for K!=3 (fallback to "band N" 1-indexed)
//   5. Detailed labels for K=3 (H1/H1, H1/H2, H2/H2)
//   6. Detailed labels for K=4..6 (full H-system)
//   7. Detailed labels for K>6 (graceful fallback to "band N")
//   8. Out-of-range bandIdx returns "?"
//   9. Caveat tooltip in detailed mode
//  10. setVocab persists in localStorage
//  11. setVocab rejects bogus values
//  12. UI: button bar present, default highlighted
//  13. UI: clicking button switches vocab + active class
//  14. Edge case: bandIdx = 0 boundary

const { JSDOM } = require('jsdom');
const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');
function fail(msg) { console.error('[lv] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[lv] ok  :', msg); }

const dom = new JSDOM(html, {
  runScripts: 'dangerously', pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window: win } = dom;

// ---- 1. Helpers exposed ----
for (const name of ['getKaryotypeLabel', 'getKaryotypeLabelCaveat',
                    '_karyoLabelEnsureVocab', '_karyoLabelSetVocab']) {
  if (typeof win[name] !== 'function') fail(name + ' not exposed');
}
if (typeof win._KARYO_DETAILED_LABELS !== 'object') fail('_KARYO_DETAILED_LABELS not exposed');
if (!Array.isArray(win._KARYO_LEGACY_LABELS_K3)) fail('_KARYO_LEGACY_LABELS_K3 not exposed');
ok('all 4 vocabulary helpers + 2 constants exposed on window');

// ---- 2. Default vocab is 'legacy' ----
// Make sure no localStorage carryover from previous tests
try { win.localStorage.removeItem('inversion_atlas.labelVocab'); } catch (_) {}
delete win.state.labelVocab;
const v0 = win._karyoLabelEnsureVocab();
if (v0 !== 'legacy') fail('default vocab should be legacy, got ' + v0);
ok('default vocab on first run: legacy');

// ---- 3. Legacy labels for K=3 ----
const lL0 = win.getKaryotypeLabel(0, 3, 'legacy');
const lL1 = win.getKaryotypeLabel(1, 3, 'legacy');
const lL2 = win.getKaryotypeLabel(2, 3, 'legacy');
if (lL0 !== 'band 1 (lo)')  fail('K=3 band 0 legacy should be "band 1 (lo)", got ' + lL0);
if (lL1 !== 'band 2 (mid)') fail('K=3 band 1 legacy should be "band 2 (mid)", got ' + lL1);
if (lL2 !== 'band 3 (hi)')  fail('K=3 band 2 legacy should be "band 3 (hi)", got ' + lL2);
ok('legacy K=3: band 1 (lo) / band 2 (mid) / band 3 (hi)');

// ---- 4. Legacy labels for K!=3 (fallback) ----
const lL_K4_0 = win.getKaryotypeLabel(0, 4, 'legacy');
const lL_K4_3 = win.getKaryotypeLabel(3, 4, 'legacy');
const lL_K6_5 = win.getKaryotypeLabel(5, 6, 'legacy');
if (lL_K4_0 !== 'band 1') fail('K=4 band 0 legacy should fall back to "band 1", got ' + lL_K4_0);
if (lL_K4_3 !== 'band 4') fail('K=4 band 3 legacy should be "band 4", got ' + lL_K4_3);
if (lL_K6_5 !== 'band 6') fail('K=6 band 5 legacy should be "band 6", got ' + lL_K6_5);
ok('legacy K!=3: 1-indexed "band N" fallback (K=4: band 1..4; K=6 band 5: band 6)');

// ---- 5. Detailed labels for K=3 ----
const lD0 = win.getKaryotypeLabel(0, 3, 'detailed');
const lD1 = win.getKaryotypeLabel(1, 3, 'detailed');
const lD2 = win.getKaryotypeLabel(2, 3, 'detailed');
if (lD0 !== 'H1/H1') fail('K=3 band 0 detailed should be H1/H1, got ' + lD0);
if (lD1 !== 'H1/H2') fail('K=3 band 1 detailed should be H1/H2, got ' + lD1);
if (lD2 !== 'H2/H2') fail('K=3 band 2 detailed should be H2/H2, got ' + lD2);
ok('detailed K=3: H1/H1 / H1/H2 / H2/H2 (biallelic 3-class system)');

// ---- 6. Detailed labels for K=4..6 ----
const lD_K4 = [0,1,2,3].map(i => win.getKaryotypeLabel(i, 4, 'detailed'));
const expectedK4 = ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3'];
for (let i = 0; i < 4; i++) {
  if (lD_K4[i] !== expectedK4[i]) fail('K=4 detailed band ' + i + ' should be ' + expectedK4[i] + ', got ' + lD_K4[i]);
}
ok('detailed K=4: ' + lD_K4.join(' / '));

const lD_K5 = [0,1,2,3,4].map(i => win.getKaryotypeLabel(i, 5, 'detailed'));
const expectedK5 = ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3', 'H3/H3'];
for (let i = 0; i < 5; i++) {
  if (lD_K5[i] !== expectedK5[i]) fail('K=5 detailed band ' + i + ' should be ' + expectedK5[i] + ', got ' + lD_K5[i]);
}
ok('detailed K=5: ' + lD_K5.join(' / '));

const lD_K6 = [0,1,2,3,4,5].map(i => win.getKaryotypeLabel(i, 6, 'detailed'));
const expectedK6 = ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3', 'H2/H3', 'H3/H3'];
for (let i = 0; i < 6; i++) {
  if (lD_K6[i] !== expectedK6[i]) fail('K=6 detailed band ' + i + ' should be ' + expectedK6[i] + ', got ' + lD_K6[i]);
}
ok('detailed K=6: ' + lD_K6.join(' / ') + ' (full 3-haplotype 6-class system)');

// ---- 7. Detailed K>6 fallback ----
const lD_K7 = win.getKaryotypeLabel(0, 7, 'detailed');
if (lD_K7 !== 'band 1') fail('K=7 detailed should fall back to "band 1", got ' + lD_K7);
const lD_K8_3 = win.getKaryotypeLabel(3, 8, 'detailed');
if (lD_K8_3 !== 'band 4') fail('K=8 band 3 detailed should fall back to "band 4", got ' + lD_K8_3);
ok('detailed K>6: graceful fallback to "band N" 1-indexed');

// ---- 8. Out-of-range bandIdx ----
if (win.getKaryotypeLabel(null, 3, 'detailed') !== '?') fail('null bandIdx should return "?"');
if (win.getKaryotypeLabel(-1, 3, 'detailed') !== '?') fail('negative bandIdx should return "?"');
ok('out-of-range bandIdx (null, -1): returns "?"');

// ---- 9. Caveat tooltip ----
win._karyoLabelSetVocab('legacy');
if (win.getKaryotypeLabelCaveat() !== null) fail('caveat in legacy mode should be null');
ok('caveat in legacy mode: null');

win._karyoLabelSetVocab('detailed');
const caveat = win.getKaryotypeLabelCaveat();
if (typeof caveat !== 'string' || caveat.length < 50) fail('caveat in detailed mode should be a substantial string, got ' + caveat);
if (!caveat.toLowerCase().includes('operational')) fail('caveat should mention "operational"');
if (!caveat.toLowerCase().includes('heterozygosity') && !caveat.toLowerCase().includes('haplotype')) {
  fail('caveat should reference heterozygosity / haplotype evidence');
}
ok('caveat in detailed mode: ' + caveat.length + ' chars, mentions "operational" + heterozygosity');

// ---- 10. setVocab persists in localStorage ----
win._karyoLabelSetVocab('detailed');
const stored = win.localStorage.getItem('inversion_atlas.labelVocab');
if (stored !== 'detailed') fail('setVocab(detailed) should persist in localStorage, got ' + stored);
ok('setVocab persists in localStorage: inversion_atlas.labelVocab = detailed');

// ---- 11. setVocab rejects bogus values ----
const ok1 = win._karyoLabelSetVocab('bogus');
if (ok1 !== false) fail('bogus vocab should be rejected, got ' + ok1);
const stillDetailed = win.localStorage.getItem('inversion_atlas.labelVocab');
if (stillDetailed !== 'detailed') fail('rejected setVocab should not change stored value');
ok('setVocab rejects bogus values, leaves stored vocab unchanged');

// ---- 12. UI: button bar present, default highlighted ----
const bar = win.document.getElementById('labelVocabBar');
if (!bar) fail('labelVocabBar missing from DOM');
const btnLegacy = bar.querySelector('button[data-vocab="legacy"]');
const btnDetailed = bar.querySelector('button[data-vocab="detailed"]');
if (!btnLegacy || !btnDetailed) fail('legacy/detailed buttons missing');
ok('UI: labelVocabBar present with legacy + detailed buttons');

// Re-init UI to reflect the stored 'detailed' from earlier
// (the init function ran on JSDOM load before we changed localStorage)
win.document.querySelectorAll('#labelVocabBar button[data-vocab]').forEach(b => {
  b.classList.toggle('active', b.dataset.vocab === win._karyoLabelEnsureVocab());
});
if (!btnDetailed.classList.contains('active')) fail('detailed button should be active after setVocab(detailed)');
ok('UI: detailed button highlighted as active when state matches');

// ---- 13. UI: clicking switches vocab + active class ----
btnLegacy.click();
if (win._karyoLabelEnsureVocab() !== 'legacy') fail('click on legacy button should set vocab to legacy');
if (!btnLegacy.classList.contains('active')) fail('legacy button should now be active');
if (btnDetailed.classList.contains('active')) fail('detailed button should no longer be active');
ok('UI: clicking legacy button → vocab=legacy, active class moves');

btnDetailed.click();
if (win._karyoLabelEnsureVocab() !== 'detailed') fail('click on detailed button should set vocab to detailed');
if (!btnDetailed.classList.contains('active')) fail('detailed button should now be active');
if (btnLegacy.classList.contains('active')) fail('legacy button should no longer be active');
ok('UI: clicking detailed button → vocab=detailed, active class moves');

// ---- 14. State default applied via state.labelVocab ----
delete win.state.labelVocab;
try { win.localStorage.removeItem('inversion_atlas.labelVocab'); } catch (_) {}
const v = win._karyoLabelEnsureVocab();
if (v !== 'legacy') fail('after clearing localStorage, default should be legacy, got ' + v);
ok('clean state: defaults back to legacy');

// ---- 15. Help-page H-system section present (HTML check, not DOM since the help page is dynamic) ----
if (!html.includes('Operational structural-haplotype framework')) {
  fail('help-page H-system section header missing');
}
if (!html.includes('H1/H1') || !html.includes('H1/H2')) {
  fail('help-page should reference H-system labels in body text');
}
if (!html.includes('manuscript-safe')) {
  fail('help-page should reference manuscript-safe wording');
}
ok('help-page: H-system framework section present, references H-labels and manuscript wording');

// ---- 16. K=3 button tagging (data-vocab-band-k3) ----
const taggedBtns = win.document.querySelectorAll('[data-vocab-band-k3]');
if (taggedBtns.length < 9) {
  fail('expected at least 9 tagged K=3 buttons (3 sites x 3 bands), got ' + taggedBtns.length);
}
ok('K=3 button tagging: ' + taggedBtns.length + ' buttons marked with data-vocab-band-k3 (>=9 expected)');

// Each tag must be 0/1/2
const tagSet = new Set();
for (const b of taggedBtns) {
  const v = parseInt(b.dataset.vocabBandK3, 10);
  if (!Number.isInteger(v) || v < 0 || v > 2) {
    fail('bad data-vocab-band-k3 value: ' + b.dataset.vocabBandK3);
  }
  tagSet.add(v);
}
if (tagSet.size !== 3) fail('tagged buttons should cover bands 0/1/2; got ' + Array.from(tagSet));
ok('all tagged buttons have valid band index (0/1/2)');

// ---- 17. _karyoLabelRefreshButtons rewrites tagged labels ----
if (typeof win._karyoLabelRefreshButtons !== 'function') {
  fail('_karyoLabelRefreshButtons not exposed on window');
}
ok('_karyoLabelRefreshButtons exposed on window');

// Switch to detailed and refresh
win._karyoLabelSetVocab('detailed');
win._karyoLabelRefreshButtons();
const firstBtn = win.document.querySelector('[data-vocab-band-k3="0"]');
if (firstBtn.textContent !== 'H1/H1') {
  fail('after detailed refresh, band-0 button should say "H1/H1", got "' + firstBtn.textContent + '"');
}
ok('refreshButtons (detailed): band-0 button -> "H1/H1"');

const secondBtn = win.document.querySelector('[data-vocab-band-k3="1"]');
if (secondBtn.textContent !== 'H1/H2') {
  fail('band-1 button should say "H1/H2", got "' + secondBtn.textContent + '"');
}
const thirdBtn = win.document.querySelector('[data-vocab-band-k3="2"]');
if (thirdBtn.textContent !== 'H2/H2') {
  fail('band-2 button should say "H2/H2", got "' + thirdBtn.textContent + '"');
}
ok('refreshButtons (detailed): band-1 -> "H1/H2", band-2 -> "H2/H2"');

// Round-trip back to legacy
win._karyoLabelSetVocab('legacy');
win._karyoLabelRefreshButtons();
const firstBtnAgain = win.document.querySelector('[data-vocab-band-k3="0"]');
if (firstBtnAgain.textContent !== 'band 1 (lo)') {
  fail('after legacy refresh, band-0 button should say "band 1 (lo)", got "' + firstBtnAgain.textContent + '"');
}
ok('refreshButtons (legacy round-trip): band-0 -> "band 1 (lo)"');

// ---- 18. Idempotency: refreshButtons can be called repeatedly ----
win._karyoLabelRefreshButtons();
win._karyoLabelRefreshButtons();
const stillCorrect = win.document.querySelector('[data-vocab-band-k3="0"]');
if (stillCorrect.textContent !== 'band 1 (lo)') {
  fail('refreshButtons should be idempotent');
}
ok('refreshButtons: idempotent across repeated calls');

console.log('\n[lv] ALL CHECKS PASSED');
