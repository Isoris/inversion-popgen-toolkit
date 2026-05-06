// =============================================================================
// turn 128d test — L3 default layout = leftright (was: focal)
//
// Quentin's bug report:
//   "The L3 contingency table should default as +1L +1R view and not single
//   tab panel."
//
// Fix: state.l3Layout default = 'leftright'; markup .active class moved
// from data-layout="focal" to data-layout="leftright".
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

console.log('\n=== L3 default layout ===');

ok('state.l3Layout default is "leftright"',
   /l3Layout:\s*'leftright'/.test(html),
   'should default to the +L/R two-neighbor view');

ok('state.l3Layout default is NOT "focal" anymore',
   !/l3Layout:\s*'focal'\s*,/.test(html));

ok('markup: data-layout="leftright" carries class="active"',
   /<button[^>]*data-layout="leftright"[^>]*class="active"|<button[^>]*class="active"[^>]*data-layout="leftright"/.test(html));

ok('markup: data-layout="focal" no longer has class="active"',
   !/<button[^>]*data-layout="focal"[^>]*class="active"/.test(html) &&
   !/<button[^>]*class="active"[^>]*data-layout="focal"/.test(html));

// Other layouts must still exist (we didn't drop them)
for (const layout of ['focal', 'left', 'right', 'leftright', 'four']) {
  ok('layout button data-layout="' + layout + '" still present',
     new RegExp('<button[^>]*data-layout="' + layout + '"').test(html));
}

// Exactly one button has the active class in the l3Layout group
const layoutBlock = (() => {
  const start = html.indexOf('id="l3Layout"');
  if (start < 0) return '';
  const end = html.indexOf('</div>', start);
  return html.substring(start, end + 6);
})();
const activeMatches = layoutBlock.match(/class="active"/g) || [];
ok('exactly one l3Layout button is active in markup',
   activeMatches.length === 1, 'got ' + activeMatches.length + ' active');

console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
