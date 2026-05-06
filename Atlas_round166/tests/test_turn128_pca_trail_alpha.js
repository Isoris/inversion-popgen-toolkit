// =============================================================================
// turn 128 integration test — PCA trail alpha attenuation on non-PC1 axes
//
// Bug from the latest user message:
//   "When PCA has 2 PCs the tracked samples 'line that tracks their location
//   when we move, must be more alpha because the lines are not horizontal
//   anymore when using 2 PCs"
//
// Interpretation: in the default pc1×pc2 view, trail lines stay relatively
// short on the X axis (PC1 changes slowly across windows for any given
// sample). When the user picks two non-PC1 axes (e.g., pc2×pc3), both axes
// vary a lot per window so the trails wander in 2D and visually dominate the
// scatter. Drop the trail line alpha in that case so the wandering is faded.
//
// This is a source-level test (no canvas evaluation) — verifies that drawPCA
// computes a lower trail alpha when neither axisX nor axisY is 'pc1', and
// keeps the previous 0.85 alpha otherwise.
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

// =============================================================================
// Source-level checks
// =============================================================================
console.log('\n=== Source-level checks ===');

// The alpha-attenuation predicate must check both axes against 'pc1'.
ok('predicate _trailNeitherIsPc1 defined as (axisX !== "pc1" && axisY !== "pc1")',
   /_trailNeitherIsPc1\s*=\s*\(axisX !== ['"]pc1['"]\s*&&\s*axisY !== ['"]pc1['"]\)/.test(html));

// The alpha must be conditionally chosen between attenuated and 0.85 default.
ok('_trailLineAlpha conditional uses _trailNeitherIsPc1 ? <attenuated> : 0.85',
   /_trailLineAlpha\s*=\s*_trailNeitherIsPc1\s*\?\s*[\d.]+\s*:\s*0\.85/.test(html));

// Specifically the attenuated value should be in (0, 0.5) — visibly lower
// than the original 0.85 but still readable.
const m = html.match(/_trailLineAlpha\s*=\s*_trailNeitherIsPc1\s*\?\s*([\d.]+)\s*:/);
ok('attenuated alpha extractable', !!m);
if (m) {
  const v = parseFloat(m[1]);
  ok('attenuated alpha is a valid number', isFinite(v));
  ok('attenuated alpha is meaningfully lower than 0.85',
     v > 0 && v < 0.5,
     'got ' + m[1]);
}

// The trail line stroke must use _trailLineAlpha (not the previous hardcoded 0.85).
ok('trail line stroke uses globalAlpha = _trailLineAlpha',
   /ctx\.globalAlpha\s*=\s*_trailLineAlpha/.test(html));

// The OLD `ctx.globalAlpha = 0.85` literal in drawPCA's trail block must be gone.
// It should appear nowhere in the file as a bare constant assigned to globalAlpha
// inside the trail rendering. We check for the specific old pattern that used
// to exist in drawPCA right before the trail line beginPath.
const oldPattern = /for \(const si of state\.tracked\) \{\s*const col = trackedColor\(si\);\s*ctx\.strokeStyle = col; ctx\.lineWidth = 1\.3;\s*ctx\.globalAlpha = 0\.85;/;
ok('old hardcoded ctx.globalAlpha = 0.85 in trail block is gone',
   !oldPattern.test(html),
   'the literal 0.85 alpha line in the trail block should now route through _trailLineAlpha');

// The dot-fade alpha (0.15 + 0.5 * t) is preserved — it implied temporal
// direction and that's still useful regardless of which axes are showing.
ok('dot-fade alpha 0.15 + 0.5 * t is preserved',
   /ctx\.globalAlpha = 0\.15 \+ 0\.5 \* t/.test(html));

// =============================================================================
// Logic test — pure boolean check of the predicate against synthetic axis combos
// =============================================================================
console.log('\n=== Predicate logic ===');

function neitherIsPc1(axisX, axisY) {
  return (axisX !== 'pc1' && axisY !== 'pc1');
}

// Default view: pc1 x pc2 → trail stays at 0.85
ok('pc1×pc2 → predicate false (full alpha)',  neitherIsPc1('pc1', 'pc2') === false);
ok('pc1×pc3 → predicate false (full alpha)',  neitherIsPc1('pc1', 'pc3') === false);
ok('pc2×pc1 → predicate false (full alpha)',  neitherIsPc1('pc2', 'pc1') === false);
ok('pc4×pc1 → predicate false (full alpha)',  neitherIsPc1('pc4', 'pc1') === false);

// Both non-PC1 — trail attenuates
ok('pc2×pc3 → predicate true (attenuated)',   neitherIsPc1('pc2', 'pc3') === true);
ok('pc3×pc2 → predicate true (attenuated)',   neitherIsPc1('pc3', 'pc2') === true);
ok('pc2×pc4 → predicate true (attenuated)',   neitherIsPc1('pc2', 'pc4') === true);
ok('pc3×pc4 → predicate true (attenuated)',   neitherIsPc1('pc3', 'pc4') === true);

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
