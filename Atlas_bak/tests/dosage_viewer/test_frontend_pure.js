// =============================================================================
// tests/dosage_viewer/test_frontend_pure.js
// =============================================================================
// Tests for the pure (DOM-free) functions in dosage_viewer/static/app.js.
//
// We can't run the full app (it needs window/document/canvas), but we
// CAN extract the colour-ramp and the LRU-pushRecent logic and run them
// in a Node sandbox. That covers:
//
//   - dosageRGB colour interpolation: 0 -> blue, 1 -> white, 2 -> red,
//     intermediates are linearly interpolated, NA returns null.
//   - LRU semantics: cap of 16 enforced; reinsertion moves to the front;
//     keys are derived from chrom:start-end:mode.
//
// Run:
//   cd Atlas && node tests/dosage_viewer/test_frontend_pure.js
// =============================================================================

'use strict';

const fs = require('fs');
const path = require('path');
const vm = require('vm');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  \u2713 ' + msg); }
  else      { fail++; console.error('  \u2717 ' + msg); }
}
function section(name) { console.log('\n\u2014 ' + name + ' \u2014'); }

const APP = fs.readFileSync(
  path.resolve(__dirname, '../../dosage_viewer/static/app.js'),
  'utf-8'
);

// ---------------------------------------------------------------------------
// Extract pure functions from app.js by string slicing. This is brittle
// in principle but the functions in question are short and stable.
// ---------------------------------------------------------------------------
function extractFn(name) {
  const re = new RegExp(
    '\\nfunction ' + name + '\\([^)]*\\)\\s*\\{[\\s\\S]*?\\n\\}\\n',
    'm'
  );
  const m = APP.match(re);
  if (!m) return null;
  return m[0];
}

const _hexToRgb_src = extractFn('_hexToRgb');
const dosageRGB_src = extractFn('dosageRGB');

assert(_hexToRgb_src, '_hexToRgb extracted');
assert(dosageRGB_src, 'dosageRGB extracted');

// Re-define RGB0/RGB1/RGB2 in the sandbox using the same hex values
const sandbox = {};
vm.createContext(sandbox);
vm.runInContext(_hexToRgb_src, sandbox);
vm.runInContext(`
  const RAMP = ['#2166AC', '#F7F7F7', '#B2182B'];
  const RGB0 = _hexToRgb(RAMP[0]);
  const RGB1 = _hexToRgb(RAMP[1]);
  const RGB2 = _hexToRgb(RAMP[2]);
`, sandbox);
vm.runInContext(dosageRGB_src, sandbox);

// ---------------------------------------------------------------------------
section('dosageRGB colour ramp');

// d == 0 should be exactly the first ramp colour
const c0 = sandbox.dosageRGB(0);
assert(c0 && c0[0] === 0x21 && c0[1] === 0x66 && c0[2] === 0xAC,
       `d=0 -> #2166AC (got rgb(${c0}))`);

// d == 1 should be exactly the middle ramp colour
const c1 = sandbox.dosageRGB(1);
assert(c1 && c1[0] === 0xF7 && c1[1] === 0xF7 && c1[2] === 0xF7,
       `d=1 -> #F7F7F7 (got rgb(${c1}))`);

// d == 2 should be exactly the last ramp colour
const c2 = sandbox.dosageRGB(2);
assert(c2 && c2[0] === 0xB2 && c2[1] === 0x18 && c2[2] === 0x2B,
       `d=2 -> #B2182B (got rgb(${c2}))`);

// d == 0.5 should be midway between #2166AC and #F7F7F7
const c05 = sandbox.dosageRGB(0.5);
assert(Math.abs(c05[0] - (0x21 + 0xF7) / 2) < 1.0,
       `d=0.5 R is midway (got ${c05[0]})`);
assert(Math.abs(c05[1] - (0x66 + 0xF7) / 2) < 1.0,
       `d=0.5 G is midway (got ${c05[1]})`);

// d == 1.5 should be midway between #F7F7F7 and #B2182B
const c15 = sandbox.dosageRGB(1.5);
assert(Math.abs(c15[0] - (0xF7 + 0xB2) / 2) < 1.0,
       `d=1.5 R is midway (got ${c15[0]})`);

// NA-style inputs return null
assert(sandbox.dosageRGB(-1) === null, 'd=-1 -> null (NA)');
assert(sandbox.dosageRGB(null) === null, 'd=null -> null');
assert(sandbox.dosageRGB(undefined) === null, 'd=undefined -> null');
assert(sandbox.dosageRGB(NaN) === null, 'd=NaN -> null');

// Out-of-range clamps (not NA — we still pick the nearest endpoint)
const cBig = sandbox.dosageRGB(10);
assert(cBig && cBig[0] === 0xB2 && cBig[2] === 0x2B,
       'd=10 clamps to ramp end (red)');

// Aggregate-mode floats (non-integer) interpolate smoothly
const cFrac = sandbox.dosageRGB(0.25);
const cFracNeighbor = sandbox.dosageRGB(0.26);
const dist = Math.abs(cFrac[0] - cFracNeighbor[0]) + Math.abs(cFrac[1] - cFracNeighbor[1]) + Math.abs(cFrac[2] - cFracNeighbor[2]);
assert(dist <= 5, `d=0.25 vs 0.26 close (got Manhattan ${dist})`);

// ---------------------------------------------------------------------------
section('LRU semantics — pushRecent');

// Re-implement pushRecent's logic in the sandbox without DOM dependencies.
// We extract the function and stub renderRecent to be a no-op.
const pushRecent_src = extractFn('pushRecent');
assert(pushRecent_src, 'pushRecent extracted');

const lruSandbox = {
  state: {
    recent: new Map(),
  },
  Date: Date,
  // No-op stub for renderRecent (it touches the DOM).
  renderRecent: function () {},
};
vm.createContext(lruSandbox);
vm.runInContext(`
  const RECENT_LRU_CAP = 16;
` + pushRecent_src, lruSandbox);

// Push 5 unique entries
for (let i = 0; i < 5; i++) {
  lruSandbox.pushRecent(
    { chrom: 'C_gar_LG28', start: i * 1000, end: (i + 1) * 1000, mode: 'raw',
      max_sites: 1000, seed: 1 },
    { n_sites_total: 100, n_sites_returned: 100 }
  );
}
assert(lruSandbox.state.recent.size === 5, `5 entries after 5 unique pushes (got ${lruSandbox.state.recent.size})`);

// Re-push the first one — should still have 5, but now it's most-recent
lruSandbox.pushRecent(
  { chrom: 'C_gar_LG28', start: 0, end: 1000, mode: 'raw',
    max_sites: 1000, seed: 1 },
  { n_sites_total: 100, n_sites_returned: 100 }
);
assert(lruSandbox.state.recent.size === 5,
  `re-push doesn't grow size (got ${lruSandbox.state.recent.size})`);

const keys = [...lruSandbox.state.recent.keys()];
assert(keys[keys.length - 1].startsWith('C_gar_LG28:0-1000'),
  'reinserted entry is now most-recent');

// Push 16 more to test the cap
for (let i = 0; i < 16; i++) {
  lruSandbox.pushRecent(
    { chrom: 'C_gar_LG12', start: i * 100, end: (i + 1) * 100, mode: 'even',
      max_sites: 500, seed: i },
    { n_sites_total: 50, n_sites_returned: 50 }
  );
}
assert(lruSandbox.state.recent.size === 16,
  `LRU capped at 16 (got ${lruSandbox.state.recent.size})`);

// Earliest entries should have been evicted
const finalKeys = [...lruSandbox.state.recent.keys()];
const oldestIsLG12 = finalKeys.every(k => k.startsWith('C_gar_LG12:'));
assert(oldestIsLG12, 'oldest entries (LG28 ones) all evicted');

// Different mode = different cache key
lruSandbox.pushRecent(
  { chrom: 'C_gar_LG12', start: 0, end: 100, mode: 'random',
    max_sites: 500, seed: 0 },
  { n_sites_total: 50, n_sites_returned: 50 }
);
const matchingMode = [...lruSandbox.state.recent.keys()].filter(k => k.includes('[random]'));
assert(matchingMode.length === 1, 'mode is part of the cache key (1 entry with [random])');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
