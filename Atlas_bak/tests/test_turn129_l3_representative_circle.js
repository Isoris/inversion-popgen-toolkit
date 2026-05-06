// =============================================================================
// turn 129 test — Representative individual nearest each group center
//
// Quentin's request (from screenshot of an R figure):
//   "You see this circle in the middle that shows the het group. I feel that
//   its useful can we have it?"
//
// Implementation: for each band g0..gK-1 with samples in the L3 mini-PCA,
// find the sample whose (PC1, PC2) is closest to the band's centroid, and
// circle it with a thin black open ring. Mirrors the R "circled =
// representative individual nearest each group center" overlay.
//
// Always-on when labels exist + colorMode === 'cluster'. Skipped on
// q_ancestry / family / none modes where the band concept doesn't apply
// the same way. Drawn between non-tracked and tracked passes so a tracked
// sample that's also a representative gets both rings.
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

console.log('\n=== Source-level checks ===');

ok('representative-circle pass added',
   /Representative individual nearest each group center/.test(html));

ok('rep pass gated on labels + cluster mode',
   /labels && state\.colorMode === 'cluster'/.test(html));

ok('rep pass computes centroid via Float64Array sumX / sumY',
   /sumX = new Float64Array\(K\)[\s\S]{0,200}sumY = new Float64Array\(K\)/.test(html));

ok('rep pass uses Euclidean distance argmin',
   /bestD = Infinity[\s\S]{0,500}d < bestD/.test(html));

ok('rep pass draws a thin black ring (themeColor ink)',
   /strokeStyle = themeColor\('ink'\);[\s\S]{0,200}arc\(rx, ry, 4\.0/.test(html));

ok('rep pass sits between non-tracked dot fill and tracked-sample loop',
   /arc\(x, y, 1\.8, 0, Math\.PI \* 2\); ctx\.fill\(\);[\s\S]{0,200}Representative individual nearest each group center[\s\S]{0,3500}Tracked samples on top/.test(html),
   'rep block must be after the non-tracked dot fill and before the tracked-sample loop');

ok('rep pass guards against finite PC values',
   /Representative individual nearest each group center[\s\S]{0,4000}!Number\.isFinite\(xv\) \|\| !Number\.isFinite\(yv\)/.test(html));

ok('rep pass guards against empty band (cnt[k] === 0)',
   /Representative individual nearest each group center[\s\S]{0,4000}cnt\[k\] === 0/.test(html));

ok('rep pass skips no-call labels (l < 0)',
   /Representative individual nearest each group center[\s\S]{0,4000}l < 0/.test(html));

ok('rep pass derives K from labels max + 1 (no hard-coded K)',
   /kMax = -1[\s\S]{0,500}kMax \+ 1/.test(html));

console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
