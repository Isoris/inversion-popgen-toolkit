// =============================================================================
// turn 127 integration test — header reorganization
//
// Verifies that:
//   1. Metadata strip is cleaned (no L1=/L2=/L2 peaks=); active count not duplicated.
//   2. Three colour-coded folders exist in the header (blue/green/yellow).
//   3. Each folder has the right children (with original IDs preserved):
//      - blue (session): saveSessionBtn, loadSessionBtn, layoutModeBtn, resetLayoutBtn
//      - green (mode): candidateModeBtn, atlasToolsAutofill, atlasToolsLabelsToggle
//      - yellow (data): atlasToolsExport, atlasToolsMatrix, activeSamplesBadge,
//                       and the dynamically-appended server status badge target
//   4. Old #atlasToolsGroup div in <nav> tab-bar is gone.
//   5. CSS for header-folder + header-folder-panel + colour variants exists.
//   6. Hover-reveal pattern uses the same :hover/:focus-within idiom as
//      #atlasModeIndicator.
//   7. atlasToolsGroup ID still exists exactly once (now inside the yellow folder)
//      so the dynamic server-badge code can still find it.
//   8. No duplicate IDs anywhere in the header.
//   9. Regression: prior turns 117-126 all still wired and findable in source.
// =============================================================================

const fs = require('fs');

const html = fs.readFileSync('/home/claude/work/build/Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// =============================================================================
// 1. Metadata strip cleanup
// =============================================================================
console.log('\n=== Metadata strip cleanup ===');

// The ONLY headerMeta innerHTML build site (line 41635) should now be the
// trimmed version without L1/L2/peaks counts.
const metaBuilds = html.match(/headerMeta'\)\.innerHTML\s*=\s*`[^`]+`/g) || [];
ok('exactly one headerMeta build site', metaBuilds.length === 1);
if (metaBuilds.length === 1) {
  const tpl = metaBuilds[0];
  ok('headerMeta no longer contains L1=', !/L1=\$\{l1c\}/.test(tpl));
  ok('headerMeta no longer contains L2=', !/L2=\$\{l2c\}/.test(tpl));
  ok('headerMeta no longer contains L2 peaks=', !/L2 peaks=\$\{bc\}/.test(tpl));
  ok('headerMeta keeps chrom', /\$\{data\.chrom\}/.test(tpl));
  ok('headerMeta keeps n_windows', /\$\{data\.n_windows\}/.test(tpl));
  ok('headerMeta keeps n_samples', /\$\{data\.n_samples\}/.test(tpl));
}

// =============================================================================
// 2. Three folders exist with the right colours
// =============================================================================
console.log('\n=== Three colour-coded folders in header ===');

ok('blue folder (data-color="blue") present',
   /class="header-folder" data-color="blue"/.test(html));
ok('green folder (data-color="green") present',
   /class="header-folder" data-color="green"/.test(html));
ok('yellow folder (data-color="yellow") present',
   /class="header-folder" data-color="yellow"/.test(html));

// =============================================================================
// 3. Each folder has the right children with preserved IDs
// =============================================================================
console.log('\n=== Folder children (original IDs preserved) ===');

// Find each folder block in order. They appear as
// <div class="header-folder" data-color="X" ...> ... </div>
function extractFolderBlock(src, color) {
  const re = new RegExp(
    '<div class="header-folder" data-color="' + color + '"[\\s\\S]*?</div>\\s*</div>',
    'm');
  const m = src.match(re);
  return m ? m[0] : null;
}

const blueBlk   = extractFolderBlock(html, 'blue');
const greenBlk  = extractFolderBlock(html, 'green');
const yellowBlk = extractFolderBlock(html, 'yellow');

ok('blue folder block extractable',   !!blueBlk);
ok('green folder block extractable',  !!greenBlk);
ok('yellow folder block extractable', !!yellowBlk);

// blue (session) children
ok('blue folder contains saveSessionBtn',
   blueBlk && /id="saveSessionBtn"/.test(blueBlk));
ok('blue folder contains loadSessionBtn',
   blueBlk && /id="loadSessionBtn"/.test(blueBlk));
ok('blue folder contains layoutModeBtn',
   blueBlk && /id="layoutModeBtn"/.test(blueBlk));
ok('blue folder contains resetLayoutBtn',
   blueBlk && /id="resetLayoutBtn"/.test(blueBlk));
ok('blue folder contains hidden file input loadSessionInput',
   blueBlk && /id="loadSessionInput"/.test(blueBlk));

// green (mode) children
ok('green folder contains candidateModeBtn',
   greenBlk && /id="candidateModeBtn"/.test(greenBlk));
ok('green folder contains atlasToolsAutofill',
   greenBlk && /id="atlasToolsAutofill"/.test(greenBlk));
ok('green folder contains atlasToolsLabelsToggle',
   greenBlk && /id="atlasToolsLabelsToggle"/.test(greenBlk));

// yellow (data) children — note: this is also the atlasToolsGroup
ok('yellow folder contains atlasToolsGroup ID (server-badge append target)',
   yellowBlk && /id="atlasToolsGroup"/.test(yellowBlk));
ok('yellow folder contains atlasToolsExport',
   yellowBlk && /id="atlasToolsExport"/.test(yellowBlk));
ok('yellow folder contains atlasToolsMatrix',
   yellowBlk && /id="atlasToolsMatrix"/.test(yellowBlk));
ok('yellow folder contains activeSamplesBadge',
   yellowBlk && /id="activeSamplesBadge"/.test(yellowBlk));

// =============================================================================
// 4. Old atlasToolsGroup is gone from <nav>
// =============================================================================
console.log('\n=== Old atlasToolsGroup removed from <nav> ===');

// There should be exactly ONE id="atlasToolsGroup" in the entire HTML, and
// it should sit inside the yellow folder, not inside <nav>.
const atlasToolsGroupCount = (html.match(/id="atlasToolsGroup"/g) || []).length;
ok('exactly one #atlasToolsGroup in HTML', atlasToolsGroupCount === 1);

// The legacy decl had `margin-left:auto;display:flex;align-items:center;gap:4px`
// inline style — should no longer be present.
ok('legacy inline-styled #atlasToolsGroup div is gone',
   !/id="atlasToolsGroup" style="margin-left:auto/.test(html));

// And specifically: the #atlasToolsGroup must NOT appear inside <nav id="tabBar">.
// (Anchor past the CSS comment at the top that mentions the literal tag.)
const navOpenStr = '<nav id="tabBar">';
// Find the LAST occurrence of the opening tag — the real one, not the comment.
const navIdx = html.lastIndexOf(navOpenStr);
const navEndIdx = html.indexOf('</nav>', navIdx);
const navBlk = html.substring(navIdx, navEndIdx);
ok('atlasToolsGroup not in <nav id="tabBar">',
   navBlk.indexOf('atlasToolsGroup') < 0);

// =============================================================================
// 5. CSS for folder buttons + popover panels
// =============================================================================
console.log('\n=== CSS for folders + popovers ===');

const cssChecks = [
  ['header .header-folder',                'folder root selector'],
  ['header .header-folder-btn',            'folder button'],
  ['header .header-folder-panel',          'folder popover'],
  ['header .header-folder[data-color="blue"]',   'blue variant'],
  ['header .header-folder[data-color="green"]',  'green variant'],
  ['header .header-folder[data-color="yellow"]', 'yellow variant'],
];
for (const [sel, label] of cssChecks) {
  ok('CSS selector ' + label + ' exists', html.indexOf(sel) >= 0);
}

// Hover/focus reveal pattern (mirrors atlasModeIndicator)
ok('hover-reveal pattern uses :hover .header-folder-panel',
   /\.header-folder:hover .header-folder-panel/.test(html));
ok('focus-within reveal too',
   /\.header-folder:focus-within .header-folder-panel/.test(html));

// Override of inline-styled atlas-tool-btn so the in-folder buttons size correctly
ok('CSS override for .atlas-tool-btn inside folder panel exists',
   /header \.header-folder-panel \.atlas-tool-btn/.test(html));

// Active-mode highlight on green folder when candidateMode active
ok('active-state ring rule for green folder exists',
   /header \.header-folder\[data-has-active="1"\]\[data-color="green"\]/.test(html));

// =============================================================================
// 6. Spacer pushes the right-side cluster right
// =============================================================================
console.log('\n=== Right-aligned right-side cluster via flex spacer ===');

// One <div style="flex: 1;"></div> spacer should sit before the first folder.
// Anchor to <body><header> via the h1 — there's a CSS comment containing
// the literal '<header>' string up in the stylesheet that we need to skip.
const h1Anchor = html.indexOf('<h1>Inversion Atlas</h1>');
const headerStart = html.lastIndexOf('<header>', h1Anchor);
const headerEnd = html.indexOf('</header>', headerStart);
const headerBlk = html.substring(headerStart, headerEnd);
ok('spacer div exists in header',
   /<div style="flex: 1;">/.test(headerBlk));

// The spacer must precede the blue folder (which is the leftmost of the cluster)
const spacerIdx = headerBlk.indexOf('<div style="flex: 1;">');
const blueIdx = headerBlk.indexOf('data-color="blue"');
ok('spacer precedes blue folder',
   spacerIdx > 0 && blueIdx > spacerIdx);

// And blue precedes green precedes yellow
const greenIdx = headerBlk.indexOf('data-color="green"');
const yellowIdx = headerBlk.indexOf('data-color="yellow"');
ok('folder order is blue → green → yellow',
   blueIdx > 0 && greenIdx > blueIdx && yellowIdx > greenIdx);

// And yellow precedes the atlas-mode indicator + theme toggle
const atlasModeIdx = headerBlk.indexOf('id="atlasModeIndicator"');
const themeIdx = headerBlk.indexOf('id="themeToggleBtn"');
ok('yellow folder precedes atlas mode indicator',
   yellowIdx > 0 && atlasModeIdx > yellowIdx);
ok('atlas mode indicator precedes theme toggle',
   atlasModeIdx > 0 && themeIdx > atlasModeIdx);

// And the legacy margin-left:auto rules are gone
ok('legacy themeToggleBtn margin-left:auto rule removed',
   !/header #themeToggleBtn \{ margin-left: auto/.test(html));
ok('legacy resetLayoutBtn margin-left:auto rule removed',
   !/header #resetLayoutBtn \{ margin-left: auto/.test(html));

// =============================================================================
// 7. No duplicate IDs across the affected section
// =============================================================================
console.log('\n=== No duplicate IDs ===');
const idsToCheck = [
  'saveSessionBtn', 'loadSessionBtn', 'layoutModeBtn', 'resetLayoutBtn',
  'candidateModeBtn', 'atlasToolsExport', 'atlasToolsMatrix',
  'atlasToolsAutofill', 'atlasToolsLabelsToggle', 'activeSamplesBadge',
  'atlasToolsGroup', 'loadSessionInput',
];
for (const id of idsToCheck) {
  const re = new RegExp('id="' + id + '"', 'g');
  const cnt = (html.match(re) || []).length;
  ok('ID #' + id + ' appears exactly once', cnt === 1,
     'count = ' + cnt);
}

// =============================================================================
// 8. Server-badge dynamic append target still works
// =============================================================================
console.log('\n=== Server-badge dynamic append still works ===');
ok('server badge code still calls getElementById(\'atlasToolsGroup\')',
   /document\.getElementById\('atlasToolsGroup'\)/.test(html));
ok('server badge appendChild call preserved',
   /group\.appendChild\(badge\)/.test(html));

// =============================================================================
// 9. Regression: prior turn fingerprints
// =============================================================================
console.log('\n=== Regression: prior turns still wired ===');
const fingerprints = [
  ['turn 117 _renderCrossSpeciesFocalVsBg', 'function _renderCrossSpeciesFocalVsBg('],
  ['turn 121 page16b tab', 'data-page="page16b"'],
  ['turn 122 _msAutoSuggestAgeModel', 'function _msAutoSuggestAgeModel('],
  ['turn 123 _isCompTEFragilityJSON', 'function _isCompTEFragilityJSON('],
  ['turn 124 _isKaryotypeLineageJSON', 'function _isKaryotypeLineageJSON('],
  ['turn 125 focal_lineage_fission', "verdict: 'focal_lineage_fission'"],
  ['turn 126 _msSummarizeRefinement', 'function _msSummarizeRefinement('],
  ['turn 126 _msPolarizeKaryotypeEventWithRefinement',
                  'function _msPolarizeKaryotypeEventWithRefinement('],
];
for (const [label, needle] of fingerprints) {
  ok(label + ' still in source', html.indexOf(needle) >= 0);
}

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
