// =============================================================================
// turn 151 — Slab het coloring parity with L2 mode
// =============================================================================
// Closes the het-coloring gap that turn 129 (SPEC_l3_het_dosage_coloring
// Slice 1) left for slab mode. The L2 path already had _computeHetRateForL2
// + _hetRateColor + drawMiniPCA override + state.l3HetColoring toggle. Slab
// mode (introduced in turn 148-150) drew uniform K-cluster-colored dots and
// ignored the toggle.
//
// Architecture:
//   - _computeHetRateForRange(start_bp, end_bp, cacheKey) — extracted core,
//     range-keyed (was inlined inside _computeHetRateForL2). Cache key is
//     arbitrary string|number; lookup short-circuits when the cache has it.
//   - _computeHetRateForL2(l2idx) — refactored: thin wrapper, pulls bp span
//     from env.start_bp/env.end_bp, delegates to range core with l2idx as key.
//   - _computeHetRateForSlab(start_w, end_w) — NEW: translates window indices
//     to bp via state.data.windows.start_bp/end_bp and delegates to range core
//     with cache key 'slab:<s>:<e>'.
//   - drawSlabMiniPCA — adds the same _hetMode / _hetRates intercept as
//     drawMiniPCA's lines 47120–47144, plus the same gradient swatch legend
//     at the bottom-left of the plot.
//
// Shared cache: state.__hetRateCache is a Map keyed by either l2idx or the
// 'slab:s:e' string. _invalidateHetRateCache clears both kinds in one drop.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Source-pattern: function declarations land
// ============================================================================
console.log('\n=== 1. Helpers declared ===');

ok('_computeHetRateForRange(start_bp, end_bp, cacheKey) declared',
   /function\s+_computeHetRateForRange\s*\(\s*start_bp\s*,\s*end_bp\s*,\s*cacheKey\s*\)/.test(html));

ok('_computeHetRateForL2(l2idx) still declared (refactored)',
   /function\s+_computeHetRateForL2\s*\(\s*l2idx\s*\)/.test(html));

ok('_computeHetRateForSlab(start_w, end_w) declared',
   /function\s+_computeHetRateForSlab\s*\(\s*start_w\s*,\s*end_w\s*\)/.test(html));

ok('_computeHetRateForRange exposed on window for tests',
   /window\._computeHetRateForRange\s*=\s*_computeHetRateForRange/.test(html));

ok('_computeHetRateForSlab exposed on window for tests',
   /window\._computeHetRateForSlab\s*=\s*_computeHetRateForSlab/.test(html));

// ============================================================================
// 2. Source-pattern: L2 path delegates to range core (no behaviour drift)
// ============================================================================
console.log('\n=== 2. L2 path now delegates to range core ===');

const l2FnRe = /function\s+_computeHetRateForL2\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const l2Match = html.match(l2FnRe);
ok('_computeHetRateForL2 body extracted', !!l2Match);
const l2Body = l2Match ? l2Match[1] : '';

ok('L2 wrapper calls _computeHetRateForRange',
   /_computeHetRateForRange\s*\(\s*env\.start_bp\s*,\s*env\.end_bp\s*,\s*l2idx\s*\)/.test(l2Body));

ok('L2 wrapper no longer contains inlined chunk fetch',
   !/popgenDosage\.getCachedChunk/.test(l2Body),
   'inlined fetch should now live in _computeHetRateForRange only');

ok('L2 wrapper preserves all-NaN fallback for missing env',
   /Float32Array\s*\(\s*nS\s*\)/.test(l2Body) &&
   /out\[si\]\s*=\s*NaN/.test(l2Body));

// ============================================================================
// 3. Source-pattern: slab path translates windows → bp
// ============================================================================
console.log('\n=== 3. Slab helper: window-index → bp translation ===');

const slabFnRe = /function\s+_computeHetRateForSlab\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const slabMatch = html.match(slabFnRe);
ok('_computeHetRateForSlab body extracted', !!slabMatch);
const slabBody = slabMatch ? slabMatch[1] : '';

ok('slab helper reads windows.start_bp / windows.end_bp',
   /w\.start_bp\[start_w\]/.test(slabBody) &&
   /w\.end_bp\[end_w\]/.test(slabBody));

ok('slab cache key is "slab:<s>:<e>"',
   /'slab:'\s*\+\s*start_w\s*\+\s*':'\s*\+\s*end_w/.test(slabBody));

ok('slab helper validates window-index range bounds',
   /start_w\s*<\s*0/.test(slabBody) &&
   /end_w\s*>=\s*w\.start_bp\.length/.test(slabBody));

ok('slab helper returns all-NaN buffer when windows unavailable',
   /!w\.start_bp\b/.test(slabBody) &&
   /Float32Array\s*\(\s*nS\s*\)/.test(slabBody));

// ============================================================================
// 4. Source-pattern: drawSlabMiniPCA wires het mode
// ============================================================================
console.log('\n=== 4. drawSlabMiniPCA references het mode ===');

const slabPcaRe = /function\s+drawSlabMiniPCA\s*\([^)]*\)\s*\{([\s\S]*?)\nfunction\s/;
const slabPcaMatch = html.match(slabPcaRe);
ok('drawSlabMiniPCA body extracted', !!slabPcaMatch);
const slabPcaBody = slabPcaMatch ? slabPcaMatch[1] : '';

ok('drawSlabMiniPCA reads state.l3HetColoring',
   /state\.l3HetColoring/.test(slabPcaBody));

ok('drawSlabMiniPCA calls _computeHetRateForSlab',
   /_computeHetRateForSlab\s*\(\s*range\[0\]\s*,\s*range\[1\]\s*\)/.test(slabPcaBody));

ok('drawSlabMiniPCA uses _hetRateColor for fill override',
   /_hetRateColor\s*\(\s*_hetRates\[si\]\s*\)/.test(slabPcaBody));

ok('drawSlabMiniPCA falls back to groupColor when het mode is off',
   /groupColor\s*\(\s*k\s*\)/.test(slabPcaBody));

ok('drawSlabMiniPCA renders gradient swatch when het mode is on',
   /createLinearGradient/.test(slabPcaBody) &&
   /addColorStop\s*\(\s*0\.50/.test(slabPcaBody));

ok('drawSlabMiniPCA swatch labels: 0 / 0.5 / 1',
   /fillText\s*\(\s*'0'/.test(slabPcaBody) &&
   /fillText\s*\(\s*'0\.5'/.test(slabPcaBody) &&
   /fillText\s*\(\s*'1'/.test(slabPcaBody));

ok('drawSlabMiniPCA swatch caption is "het"',
   /fillText\s*\(\s*'het'/.test(slabPcaBody));

// ============================================================================
// 5. Sandboxed: _computeHetRateForRange unit tests
// ============================================================================
console.log('\n=== 5. _computeHetRateForRange unit tests ===');

// Build a minimal sandbox with the helper extracted.
function makeSandbox(opts) {
  opts = opts || {};
  const ctx = {
    state: opts.state || {
      data: {
        n_samples: 4,
        samples: [
          { id: 'CGA001' }, { id: 'CGA002' },
          { id: 'CGA003' }, { id: 'CGA004' },
        ],
      },
    },
    Float32Array, Int32Array, Map, Number, console,
    window: opts.window || {},
  };
  // Need _getHetRateCache too.
  ctx.state.__hetRateCache = ctx.state.__hetRateCache || new Map();
  vm.createContext(ctx);
  // Extract just the three helpers + _hetRateColor + cache helpers.
  const fnNames = [
    '_HET_RAMP', '_hetRateColor',
    '_getHetRateCache', '_invalidateHetRateCache',
    '_computeHetRateForRange', '_computeHetRateForL2', '_computeHetRateForSlab',
  ];
  // Pick out each by its body. We search for the function/const declaration
  // and then find the matching close brace by counting.
  const extracted = [];
  for (const name of fnNames) {
    let re;
    if (name === '_HET_RAMP') {
      re = new RegExp('const\\s+' + name + '\\s*=\\s*\\{[\\s\\S]*?\\n\\};', 'm');
    } else {
      re = new RegExp('function\\s+' + name + '\\s*\\([^)]*\\)\\s*\\{[\\s\\S]*?\\n\\}', 'm');
    }
    const m = html.match(re);
    if (m) extracted.push(m[0]);
  }
  const code = extracted.join('\n\n');
  vm.runInContext(code, ctx);
  return ctx;
}

const sb = makeSandbox();
ok('sandbox: _computeHetRateForRange exists',
   typeof sb._computeHetRateForRange === 'function');

// 5a. No chunk available → all-NaN buffer
{
  const out = sb._computeHetRateForRange(1000, 2000, 'test_no_chunk');
  ok('no chunk → returns Float32Array of length n_samples',
     out instanceof Float32Array && out.length === 4);
  ok('no chunk → all entries are NaN',
     out.every(v => Number.isNaN(v)));
}

// 5b. Cache hit short-circuits
{
  const sb2 = makeSandbox();
  const first = sb2._computeHetRateForRange(0, 1, 'kkey');
  // Mutate to detect re-compute
  first[0] = 0.42;
  const second = sb2._computeHetRateForRange(0, 1, 'kkey');
  ok('cache hit returns same buffer (mutation visible)',
     Math.abs(second[0] - 0.42) < 1e-6);
}

// 5c. Invalid range → all-NaN buffer
{
  const out = sb._computeHetRateForRange(NaN, 1000, 'bad1');
  ok('NaN start_bp → all-NaN buffer',
     out instanceof Float32Array && out.every(v => Number.isNaN(v)));
}
{
  const out = sb._computeHetRateForRange(1000, 500, 'bad2');
  ok('end_bp < start_bp → all-NaN buffer',
     out instanceof Float32Array && out.every(v => Number.isNaN(v)));
}

// 5d. Real chunk path: synthetic chunk at bp 100..200, 4 samples, 5 markers
{
  const sb3 = makeSandbox();
  // Synthetic chunk: 5 markers at pos 110,120,130,140,150. 4 samples.
  // Sample 0: all hom-ref (0). Sample 1: all het (1). Sample 2: half het.
  // Sample 3: missing (NaN).
  const chunk = {
    samples: ['CGA001', 'CGA002', 'CGA003', 'CGA004'],
    markers: [
      { pos_bp: 110 }, { pos_bp: 120 }, { pos_bp: 130 },
      { pos_bp: 140 }, { pos_bp: 150 },
    ],
    dosage: [
      [0, 1, 1, NaN],
      [0, 1, 0, NaN],
      [0, 1, 1, NaN],
      [0, 1, 0, NaN],
      [0, 1, 1, NaN],
    ],
  };
  sb3.window.popgenDosage = {
    getCachedChunk: () => chunk,
  };
  const out = sb3._computeHetRateForRange(100, 200, 'real');
  ok('real chunk → sample 0 all-hom → het rate 0',
     out[0] === 0);
  ok('real chunk → sample 1 all-het → het rate 1',
     out[1] === 1);
  ok('real chunk → sample 2 mixed → het rate 0.6 (3/5)',
     Math.abs(out[2] - 0.6) < 1e-6);
  ok('real chunk → sample 3 all-NaN → het rate NaN',
     Number.isNaN(out[3]));
}

// 5e. Marker filter excludes out-of-range markers
{
  const sb4 = makeSandbox();
  const chunk = {
    samples: ['CGA001', 'CGA002', 'CGA003', 'CGA004'],
    markers: [
      { pos_bp: 50 },   // out of range (below)
      { pos_bp: 110 },  // in range
      { pos_bp: 120 },  // in range
      { pos_bp: 250 },  // out of range (above)
    ],
    dosage: [
      // The outside-range markers would push sample 0's rate up if included.
      [1, 0, 0, 1],   // pos 50: ignored
      [0, 0, 0, 0],   // pos 110: in
      [0, 0, 0, 0],   // pos 120: in
      [1, 1, 1, 1],   // pos 250: ignored
    ],
  };
  sb4.window.popgenDosage = { getCachedChunk: () => chunk };
  const out = sb4._computeHetRateForRange(100, 200, 'filt');
  ok('marker filter: out-of-range markers excluded (sample 0 rate = 0)',
     out[0] === 0);
}

// ============================================================================
// 6. Sandboxed: _computeHetRateForSlab unit tests
// ============================================================================
console.log('\n=== 6. _computeHetRateForSlab unit tests ===');

{
  const sb = makeSandbox({
    state: {
      data: {
        n_samples: 4,
        samples: [
          { id: 'CGA001' }, { id: 'CGA002' },
          { id: 'CGA003' }, { id: 'CGA004' },
        ],
        windows: {
          start_bp: [100, 200, 300, 400, 500],
          end_bp:   [199, 299, 399, 499, 599],
        },
      },
    },
  });
  // 6a. No chunk → all-NaN
  const out = sb._computeHetRateForSlab(1, 3);
  ok('slab no-chunk → Float32Array of length n_samples',
     out instanceof Float32Array && out.length === 4);
  ok('slab no-chunk → all NaN',
     out.every(v => Number.isNaN(v)));

  // 6b. Out-of-bounds window indices
  const oob = sb._computeHetRateForSlab(-1, 3);
  ok('slab negative start_w → all-NaN buffer',
     oob instanceof Float32Array && oob.every(v => Number.isNaN(v)));

  const oob2 = sb._computeHetRateForSlab(1, 99);
  ok('slab end_w >= n_windows → all-NaN buffer',
     oob2 instanceof Float32Array && oob2.every(v => Number.isNaN(v)));

  const oob3 = sb._computeHetRateForSlab(3, 1);
  ok('slab end_w < start_w → all-NaN buffer',
     oob3 instanceof Float32Array && oob3.every(v => Number.isNaN(v)));
}

// 6c. Slab path translates window indices to bp range via the windows arrays
{
  const sb = makeSandbox({
    state: {
      data: {
        n_samples: 4,
        samples: [
          { id: 'CGA001' }, { id: 'CGA002' },
          { id: 'CGA003' }, { id: 'CGA004' },
        ],
        windows: {
          start_bp: [100, 200, 300, 400, 500],
          end_bp:   [199, 299, 399, 499, 599],
        },
      },
    },
  });
  // Slab covers windows 1..3 → bp 200..499. Place markers in/out of range.
  const chunk = {
    samples: ['CGA001', 'CGA002', 'CGA003', 'CGA004'],
    markers: [
      { pos_bp: 150 },  // in window 0 (out of slab)
      { pos_bp: 250 },  // in window 1 (in slab)
      { pos_bp: 350 },  // in window 2 (in slab)
      { pos_bp: 450 },  // in window 3 (in slab)
      { pos_bp: 550 },  // in window 4 (out of slab)
    ],
    dosage: [
      [1, 1, 1, 1],   // out: would all read 1 if included
      [0, 1, 1, 0],   // in
      [0, 1, 1, 0],   // in
      [0, 1, 1, 0],   // in
      [1, 1, 1, 1],   // out
    ],
  };
  sb.window.popgenDosage = { getCachedChunk: () => chunk };
  const out = sb._computeHetRateForSlab(1, 3);
  // Each sample sees 3 in-range markers.
  ok('slab marker filter: sample 0 → 0 het / 3 = 0',
     out[0] === 0);
  ok('slab marker filter: sample 1 → 3 het / 3 = 1',
     out[1] === 1);
  ok('slab marker filter: sample 2 → 3 het / 3 = 1',
     out[2] === 1);
  ok('slab marker filter: sample 3 → 0 het / 3 = 0',
     out[3] === 0);

  // Slab cache key uses 'slab:1:3'
  ok('slab path populates cache under "slab:s:e" key',
     sb.state.__hetRateCache.has('slab:1:3'));
}

// 6d. L2 and slab paths are independent (different cache keys)
{
  const sb = makeSandbox({
    state: {
      data: {
        n_samples: 4,
        samples: [
          { id: 'CGA001' }, { id: 'CGA002' },
          { id: 'CGA003' }, { id: 'CGA004' },
        ],
        windows: {
          start_bp: [100, 200, 300],
          end_bp:   [199, 299, 399],
        },
        l2_envelopes: [
          { start_bp: 100, end_bp: 399 },   // l2idx=0 same bp range
        ],
      },
    },
  });
  const chunk = {
    samples: ['CGA001', 'CGA002', 'CGA003', 'CGA004'],
    markers: [{ pos_bp: 250 }],
    dosage: [[1, 0, 1, 0]],
  };
  sb.window.popgenDosage = { getCachedChunk: () => chunk };
  const l2Out = sb._computeHetRateForL2(0);
  const slabOut = sb._computeHetRateForSlab(0, 2);
  ok('L2 path produces same het rates as slab over identical bp range',
     l2Out[0] === slabOut[0] && l2Out[1] === slabOut[1] &&
     l2Out[2] === slabOut[2] && l2Out[3] === slabOut[3]);
  ok('but the two paths use different cache keys',
     sb.state.__hetRateCache.has(0) &&
     sb.state.__hetRateCache.has('slab:0:2'));

  // Cache invalidation drops both kinds in one call.
  sb._invalidateHetRateCache();
  ok('_invalidateHetRateCache clears both L2 and slab keys',
     !sb.state.__hetRateCache.has(0) &&
     !sb.state.__hetRateCache.has('slab:0:2'));
}

// ============================================================================
// 7. Negative / regression: L2 path preserved
// ============================================================================
console.log('\n=== 7. L2 path preserved (no regression) ===');

{
  const sb = makeSandbox({
    state: {
      data: {
        n_samples: 3,
        samples: [{ id: 'A' }, { id: 'B' }, { id: 'C' }],
        l2_envelopes: [
          { start_bp: 1000, end_bp: 2000 },
        ],
      },
    },
  });
  // 7a. Missing env → all-NaN
  const out = sb._computeHetRateForL2(99);
  ok('L2 missing env → all-NaN buffer',
     out instanceof Float32Array && out.length === 3 &&
     out.every(v => Number.isNaN(v)));

  // 7b. Valid l2idx but no chunk → all-NaN, cached
  const out2 = sb._computeHetRateForL2(0);
  ok('L2 no-chunk → all-NaN buffer',
     out2 instanceof Float32Array && out2.every(v => Number.isNaN(v)));
  ok('L2 no-chunk → cached under l2idx key (number, not string)',
     sb.state.__hetRateCache.has(0));
}

// ============================================================================
// 8. Source-pattern: legend swatch sized + positioned consistently
// ============================================================================
console.log('\n=== 8. Slab swatch dimensions match L2 swatch ===');

ok('slab swatch width = 36 (matches L2)',
   /swW\s*=\s*36/.test(slabPcaBody));
ok('slab swatch height = 6 (matches L2)',
   /swH\s*=\s*6/.test(slabPcaBody));
ok('slab swatch backdrop uses themeColor("bg") at 0.92 alpha',
   /withAlpha\s*\(\s*themeColor\s*\(\s*'bg'\s*\)\s*,\s*0\.92\s*\)/.test(slabPcaBody));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
