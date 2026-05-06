// =============================================================================
// turn 157 — V-shape tooltip on hover
// =============================================================================
// Wires the hits[] array returned by _drawVShapePlot (turn 156) into a
// hover tooltip showing per-sample CGA, group, u, agreement_fraction,
// coherence_class, stripe_quality.
//
// What this turn ships:
//   - _vsTooltipEnsureEl(): lazy-creates #vsPointTooltip div on body
//   - _vsTooltipBuildHtml(hit): renders sample-detail HTML
//   - _vsTooltipShow(hit, x, y): position with edge-clamping
//   - _vsTooltipHide(): hide
//   - _vsHitTest(hits, mx, my): closest-point hit test (handles overlapping
//     points by picking smallest-distance hit within radius)
//   - _wireVShapeTooltip(canvas): mousemove + mouseleave handlers, idempotent
//   - _vsLastRender module-level snapshot, populated inside openVShapePlot
//     after _drawVShapePlot returns
//   - openVShapePlot wires _vsLastRender + _wireVShapeTooltip after draw
//   - closeVShapePlot also calls _vsTooltipHide
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

function pullFunction(src, name) {
  const decl = `function ${name}`;
  const i = src.indexOf(decl);
  if (i < 0) return null;
  let p = src.indexOf('(', i);
  if (p < 0) return null;
  let braceStart = src.indexOf('{', p);
  if (braceStart < 0) return null;
  let depth = 1, j = braceStart + 1;
  while (j < src.length && depth > 0) {
    const ch = src[j];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    j++;
  }
  return src.substring(i, j);
}

// ============================================================================
// 1. Source-pattern: helpers declared and exposed
// ============================================================================
console.log('\n=== 1. Tooltip helpers declared ===');

ok('_vsTooltipEnsureEl declared',
   /function\s+_vsTooltipEnsureEl\s*\(\s*\)/.test(html));

ok('_vsTooltipBuildHtml declared',
   /function\s+_vsTooltipBuildHtml\s*\(\s*hit\s*\)/.test(html));

ok('_vsTooltipShow declared',
   /function\s+_vsTooltipShow\s*\(\s*hit\s*,\s*clientX\s*,\s*clientY\s*\)/.test(html));

ok('_vsTooltipHide declared',
   /function\s+_vsTooltipHide\s*\(\s*\)/.test(html));

ok('_vsHitTest declared',
   /function\s+_vsHitTest\s*\(\s*hits\s*,\s*mouseX\s*,\s*mouseY\s*\)/.test(html));

ok('_wireVShapeTooltip declared',
   /function\s+_wireVShapeTooltip\s*\(\s*canvas\s*\)/.test(html));

ok('_vsLastRender declared as module-level let',
   /let\s+_vsLastRender\s*=\s*null/.test(html));

ok('all 6 tooltip functions exposed on window',
   /window\._vsTooltipEnsureEl\b/.test(html) &&
   /window\._vsTooltipBuildHtml\b/.test(html) &&
   /window\._vsTooltipShow\b/.test(html) &&
   /window\._vsTooltipHide\b/.test(html) &&
   /window\._vsHitTest\b/.test(html) &&
   /window\._wireVShapeTooltip\b/.test(html));

// ============================================================================
// 2. Source-pattern: _vsTooltipEnsureEl behaviour
// ============================================================================
console.log('\n=== 2. _vsTooltipEnsureEl ===');

const ensureBody = pullFunction(html, '_vsTooltipEnsureEl');
ok('uses unique id #vsPointTooltip (no collision with #inhPillTooltip)',
   /'vsPointTooltip'/.test(ensureBody) && !/inhPillTooltip/.test(ensureBody));

ok('returns null when document is undefined (sandbox-safe)',
   /typeof\s+document\s*===\s*'undefined'/.test(ensureBody));

ok('reuses existing element on repeat call (idempotent)',
   /if\s*\(\s*tip\s*\)\s*return\s+tip/.test(ensureBody));

ok('appends to document.body',
   /document\.body\.appendChild/.test(ensureBody));

ok('positions absolutely, pointer-events:none, high z-index',
   /position\s*:\s*fixed/.test(ensureBody) &&
   /pointer-events\s*:\s*none/.test(ensureBody) &&
   /z-index\s*:\s*10000/.test(ensureBody));

// ============================================================================
// 3. Source-pattern: _vsTooltipBuildHtml content
// ============================================================================
console.log('\n=== 3. _vsTooltipBuildHtml content ===');

const buildBody = pullFunction(html, '_vsTooltipBuildHtml');

ok('returns empty string for missing hit / hit.point',
   /!hit\s*\|\|\s*!hit\.point/.test(buildBody));

ok('renders cga as header (with fallback "?")',
   /p\.cga\s*\|\|\s*'\?'/.test(buildBody));

ok('renders group with color chip + text',
   /_vShapeColor\s*\(\s*p\.group\s*\)/.test(buildBody));

ok('renders u (formatted to 4 decimals)',
   /p\.u\.toFixed\s*\(\s*4\s*\)/.test(buildBody));

ok('renders agreement (formatted to 4 decimals)',
   /p\.agreement\.toFixed\s*\(\s*4\s*\)/.test(buildBody));

ok('renders coherence_class with semantic color',
   /p\.coherence_class/.test(buildBody) &&
   /'coherent'/.test(buildBody) &&
   /'discordant'/.test(buildBody));

ok('renders stripe_quality tier with semantic color',
   /p\.stripe_quality/.test(buildBody) &&
   /'core'/.test(buildBody) &&
   /'junk'/.test(buildBody));

// ============================================================================
// 4. Source-pattern: _vsTooltipShow positioning
// ============================================================================
console.log('\n=== 4. _vsTooltipShow positioning ===');

const showBody = pullFunction(html, '_vsTooltipShow');

ok('default position: cursor + (margin, margin)',
   /clientX\s*\+\s*margin/.test(showBody) &&
   /clientY\s*\+\s*margin/.test(showBody));

ok('flips to cursor - rect.width when right edge would clip',
   /clientX\s*-\s*rect\.width\s*-\s*margin/.test(showBody));

ok('flips above when bottom edge would clip',
   /clientY\s*-\s*rect\.height\s*-\s*margin/.test(showBody));

ok('clamps to min 4px from viewport edge',
   /Math\.max\s*\(\s*4\s*,/.test(showBody));

// ============================================================================
// 5. Source-pattern: _wireVShapeTooltip
// ============================================================================
console.log('\n=== 5. _wireVShapeTooltip ===');

const wireBody = pullFunction(html, '_wireVShapeTooltip');

ok('idempotent — bails when canvas.dataset.vsTooltipWired already set',
   /dataset\.vsTooltipWired\s*===\s*'1'/.test(wireBody));

ok('sets the dataset marker on first call',
   /dataset\.vsTooltipWired\s*=\s*'1'/.test(wireBody));

ok('reads from _vsLastRender (not from canvas attributes)',
   /_vsLastRender/.test(wireBody));

ok('uses CSS-pixel coords directly (NO DPR scaling vs inheritance pill pattern)',
   !/canvas\.width\s*\/\s*rect\.width/.test(wireBody) &&
   /ev\.clientX\s*-\s*rect\.left/.test(wireBody));

ok('sets cursor:pointer on hover, default off',
   /cursor\s*=\s*'pointer'/.test(wireBody) &&
   /cursor\s*=\s*'default'/.test(wireBody));

ok('listens to mouseleave to hide tooltip',
   /'mouseleave'/.test(wireBody));

// ============================================================================
// 6. Source-pattern: openVShapePlot wires _vsLastRender + tooltip
// ============================================================================
console.log('\n=== 6. openVShapePlot integration ===');

const openBody = pullFunction(html, 'openVShapePlot');

ok('openVShapePlot stores _vsLastRender after _drawVShapePlot returns',
   /_vsLastRender\s*=\s*\{/.test(openBody));

ok('snapshot includes hits[] and data',
   /hits\s*:[\s\S]*?result\.hits/.test(openBody) &&
   /data\s*:\s*data/.test(openBody));

ok('snapshot includes candidate_id (defensive — survives re-open)',
   /candidate_id\s*:\s*cand\.id/.test(openBody));

ok('openVShapePlot calls _wireVShapeTooltip(canvas) after draw',
   /_wireVShapeTooltip\s*\(\s*canvas\s*\)/.test(openBody));

ok('_wireVShapeTooltip call is in try/catch (non-blocking)',
   /try\s*\{\s*_wireVShapeTooltip\s*\(\s*canvas\s*\)/.test(openBody));

const closeBody = pullFunction(html, 'closeVShapePlot');
ok('closeVShapePlot also hides tooltip (no orphan tooltip on close)',
   /_vsTooltipHide\s*\(\s*\)/.test(closeBody));

// ============================================================================
// 7. Sandboxed: _vsHitTest behaviour
// ============================================================================
console.log('\n=== 7. _vsHitTest sandboxed ===');

function makeHitTestSandbox() {
  const ctx = { Math, Number, Array, Object, console, Infinity };
  vm.createContext(ctx);
  const fn = pullFunction(html, '_vsHitTest');
  vm.runInContext(fn, ctx);
  return ctx;
}

// 7a. Empty / null hits → null
{
  const sb = makeHitTestSandbox();
  ok('null hits → null', sb._vsHitTest(null, 100, 100) === null);
  ok('empty hits → null', sb._vsHitTest([], 100, 100) === null);
  ok('non-array hits → null', sb._vsHitTest('foo', 100, 100) === null);
}

// 7b. Direct hit on point center
{
  const sb = makeHitTestSandbox();
  const hits = [
    { x: 100, y: 100, r: 3.2, point: { cga: 'A' } },
    { x: 200, y: 200, r: 3.2, point: { cga: 'B' } },
  ];
  const hit = sb._vsHitTest(hits, 100, 100);
  ok('direct hit on point A center returns A', hit && hit.point.cga === 'A');
}

// 7c. Within fudge radius
{
  const sb = makeHitTestSandbox();
  const hits = [{ x: 100, y: 100, r: 3.2, point: { cga: 'A' } }];
  // r=3.2 + fudge 1.5 = 4.7. 4.0 away should hit.
  const hit1 = sb._vsHitTest(hits, 104, 100);
  ok('cursor 4px from center → hit', hit1 && hit1.point.cga === 'A');
  // 6.0 away should not hit.
  const hit2 = sb._vsHitTest(hits, 106, 100);
  ok('cursor 6px from center → miss', hit2 === null);
}

// 7d. Closest-point picking when overlapping
{
  const sb = makeHitTestSandbox();
  // Two points right next to each other, both within radius of cursor.
  const hits = [
    { x: 99,  y: 100, r: 3.2, point: { cga: 'FAR' } },   // dist^2 = 1
    { x: 100, y: 100, r: 3.2, point: { cga: 'NEAR' } },   // dist^2 = 0
  ];
  const hit = sb._vsHitTest(hits, 100, 100);
  ok('overlapping points → picks closest (NEAR)', hit && hit.point.cga === 'NEAR');
}

// 7e. Default radius when r unset
{
  const sb = makeHitTestSandbox();
  const hits = [{ x: 100, y: 100, point: { cga: 'A' } }];   // no r
  const hit = sb._vsHitTest(hits, 102, 100);   // 2px away
  ok('missing r uses 3.2 default', hit && hit.point.cga === 'A');
}

// ============================================================================
// 8. Sandboxed: _vsTooltipBuildHtml
// ============================================================================
console.log('\n=== 8. _vsTooltipBuildHtml sandboxed ===');

function makeBuildSandbox() {
  const ctx = { Number, Object, Array, console };
  vm.createContext(ctx);
  // _vShapeColor is referenced; stub it.
  vm.runInContext('var _vShapeColor = (g) => g === "HET" ? "#3cc08a" : "#888";', ctx);
  const fn = pullFunction(html, '_vsTooltipBuildHtml');
  vm.runInContext(fn, ctx);
  return ctx;
}

// 8a. Missing hit → empty string
{
  const sb = makeBuildSandbox();
  ok('null hit → empty string', sb._vsTooltipBuildHtml(null) === '');
  ok('hit without point → empty string', sb._vsTooltipBuildHtml({}) === '');
}

// 8b. Full point renders all fields
{
  const sb = makeBuildSandbox();
  const hit = {
    point: {
      cga: 'CGA042',
      group: 'HET',
      u: 0.12345,
      agreement: 0.6789,
      coherence_class: 'coherent',
      stripe_quality: 'core',
    },
  };
  const html2 = sb._vsTooltipBuildHtml(hit);
  ok('rendered HTML includes CGA', html2.indexOf('CGA042') >= 0);
  ok('rendered HTML includes group',  html2.indexOf('HET') >= 0);
  ok('rendered HTML includes u to 4 decimals', html2.indexOf('0.1235') >= 0 || html2.indexOf('0.1234') >= 0);
  ok('rendered HTML includes agreement to 4 decimals', html2.indexOf('0.6789') >= 0);
  ok('rendered HTML includes coherence_class', html2.indexOf('coherent') >= 0);
  ok('rendered HTML includes stripe_quality',  html2.indexOf('core') >= 0);
}

// 8c. Missing optional fields → still renders header
{
  const sb = makeBuildSandbox();
  const hit = {
    point: { cga: 'X', group: 'HOMO_1', u: 1.0, agreement: 0.9 },
  };
  const html3 = sb._vsTooltipBuildHtml(hit);
  ok('without coherence_class/stripe_quality → still includes cga + u + agreement',
     html3.indexOf('X') >= 0 && html3.indexOf('1.0000') >= 0 && html3.indexOf('0.9000') >= 0);
  ok('without optional fields → no "core" or "coherent" leaked',
     html3.indexOf('core') < 0 && html3.indexOf('coherent') < 0);
}

// 8d. Non-finite u/agreement → "NA" instead of crashing
{
  const sb = makeBuildSandbox();
  const hit = {
    point: { cga: 'Y', group: 'HET', u: NaN, agreement: undefined },
  };
  const html4 = sb._vsTooltipBuildHtml(hit);
  ok('NaN u + undefined agreement → "NA" placeholders',
     (html4.match(/NA/g) || []).length >= 2);
}

// ============================================================================
// 9. Existing flow preserved
// ============================================================================
console.log('\n=== 9. Existing flow preserved ===');

ok('turn 156 _drawVShapePlot still returns hits[]',
   /hits\s*:\s*hits/.test(html));

ok('turn 156 _buildVShapeData still exists',
   /function\s+_buildVShapeData\s*\(\s*candidate\s*,\s*chunk\s*,\s*sqRows\s*\)/.test(html));

ok('turn 156 _vShapeColor still exists',
   /function\s+_vShapeColor\s*\(\s*group\s*\)/.test(html));

ok('inheritance pill tooltip pattern (turn 122) still exists (no name collision)',
   /function\s+_inhTooltipShow\s*\(/.test(html) &&
   /function\s+_wireInheritancePillTooltip\s*\(/.test(html));

ok('turn 152 _gPanelRenderTabInheritance still exists',
   /function\s+_gPanelRenderTabInheritance\s*\(\s*\)/.test(html));

ok('turn 155 _inheritanceCacheKey still has 3-arg signature',
   /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
