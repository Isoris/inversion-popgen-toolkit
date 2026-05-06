// =============================================================================
// turn 158 — popstats track categorization (QC vs popstats) + off-by-default
// =============================================================================
// Closes Quentin's items #5 and #6 from the turn 157 handoff:
//   #5 — separate QC vs popstats sections in the right-hand selector
//   #6 — turn off popstats AND QC tracks by default
//
// What this turn adds:
//   - `category` field on each POPSTATS_TRACKS entry: 'always' | 'qc' | 'popstats'
//   - _popstatsCategoryOf(track) — defaults to 'other' for auto-discovered
//     tracks without explicit category (legacy behaviour preserved)
//   - loadPopstatsView returns { hidden: Set, shown: Set } instead of a bare
//     hidden Set. The shown set is new and stores explicit user opt-ins for
//     qc/popstats categories.
//   - savePopstatsView accepts the same {hidden, shown} shape (or a bare
//     Set for back-compat with any caller that hadn't migrated).
//   - renderPopstatsPage's chip strip groups chips by category with section
//     labels ("QC:", "popstats:"). Visibility check is per-category:
//       * always  → visible if hasData
//       * qc      → visible iff explicitly shown
//       * popstats → visible iff explicitly shown
//       * other (auto-discovered) → visible if hasData (legacy)
//   - Chip click handler dispatches based on category: qc/popstats toggle
//     the shown set; always/other toggle the legacy hidden set.
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
// 1. Source-pattern: every POPSTATS_TRACKS entry has category
// ============================================================================
console.log('\n=== 1. Track categorization ===');

// Pull the POPSTATS_TRACKS array body
const trackArrRe = /const\s+POPSTATS_TRACKS\s*=\s*\[([\s\S]*?)\n\];/;
const trackArrMatch = html.match(trackArrRe);
ok('POPSTATS_TRACKS array body extracted', !!trackArrMatch);
const trackArrBody = trackArrMatch ? trackArrMatch[1] : '';

// Each entry has the shape `{ id: 'X', label: ..., category: '...' }`. Count.
const entryCount = (trackArrBody.match(/\{\s*id:\s*'[^']+'/g) || []).length;
ok('POPSTATS_TRACKS has 12 entries (no entries dropped)',
   entryCount === 12, 'got ' + entryCount);

const catCount = (trackArrBody.match(/category:\s*'(?:always|qc|popstats)'/g) || []).length;
ok('every track entry has a category field',
   catCount === entryCount, catCount + '/' + entryCount + ' have category');

// Specific category assignments
ok('ideogram is category=always',
   /id:\s*'ideogram'[\s\S]{0,150}?category:\s*'always'/.test(trackArrBody));

ok('sim_collapse is category=always',
   /id:\s*'sim_collapse'[\s\S]{0,150}?category:\s*'always'/.test(trackArrBody));

ok('z is category=always',
   /id:\s*'z'[\s\S]{0,200}?category:\s*'always'/.test(trackArrBody));

for (const id of ['snp_density', 'beagle_unc', 'coverage', 'low_cov_count']) {
  ok(id + ' is category=qc',
     new RegExp("id:\\s*'" + id + "'[\\s\\S]{0,200}?category:\\s*'qc'").test(trackArrBody));
}

for (const id of ['theta_invgt', 'fst_hom1_hom2', 'hobs_hexp', 'delta12', 'delta12_multi']) {
  ok(id + ' is category=popstats',
     new RegExp("id:\\s*'" + id + "'[\\s\\S]{0,200}?category:\\s*'popstats'").test(trackArrBody));
}

// ============================================================================
// 2. Source-pattern: _popstatsCategoryOf helper
// ============================================================================
console.log('\n=== 2. _popstatsCategoryOf helper ===');

ok('_popstatsCategoryOf declared',
   /function\s+_popstatsCategoryOf\s*\(\s*track\s*\)/.test(html));

ok('_popstatsCategoryOf falls back to "other" for missing/non-string category',
   /return\s*\([\s\S]{0,200}?'other'/.test(pullFunction(html, '_popstatsCategoryOf') || ''));

ok('_popstatsCategoryOf exposed on window',
   /window\._popstatsCategoryOf\s*=\s*_popstatsCategoryOf/.test(html));

// ============================================================================
// 3. Source-pattern: loadPopstatsView/savePopstatsView new shape
// ============================================================================
console.log('\n=== 3. loadPopstatsView/savePopstatsView shape ===');

const loadBody = pullFunction(html, 'loadPopstatsView');
ok('loadPopstatsView body extracted', !!loadBody);
ok('loadPopstatsView returns { hidden, shown } (both Sets)',
   /hidden:\s*new Set\(\)[\s\S]{0,200}?shown:\s*new Set\(\)/.test(loadBody));
ok('loadPopstatsView reads hiddenChips array from storage',
   /Array\.isArray\(\s*v\.hiddenChips\s*\)/.test(loadBody));
ok('loadPopstatsView reads shownChips array from storage',
   /Array\.isArray\(\s*v\.shownChips\s*\)/.test(loadBody));

const saveBody = pullFunction(html, 'savePopstatsView');
ok('savePopstatsView body extracted', !!saveBody);
ok('savePopstatsView accepts {hidden, shown} object',
   /view\.hidden\s+instanceof\s+Set/.test(saveBody) &&
   /view\.shown\s+instanceof\s+Set/.test(saveBody));
ok('savePopstatsView accepts bare Set for back-compat',
   /view\s+instanceof\s+Set/.test(saveBody));
ok('savePopstatsView persists both hiddenChips + shownChips',
   /hiddenChips:\s*Array\.from\(hidden\)/.test(saveBody) &&
   /shownChips:\s*Array\.from\(shown\)/.test(saveBody));

// ============================================================================
// 4. Source-pattern: renderPopstatsPage visibility logic + grouped chips
// ============================================================================
console.log('\n=== 4. renderPopstatsPage visibility + chip grouping ===');

const rpBody = pullFunction(html, 'renderPopstatsPage');
ok('renderPopstatsPage body extracted', !!rpBody);

ok('visibility uses view.hidden + view.shown sets',
   /view\.hidden\.has\s*\(\s*t\.id\s*\)/.test(rpBody) &&
   /view\.shown\.has\s*\(\s*t\.id\s*\)/.test(rpBody));

ok('alwaysOn skips all checks',
   /if\s*\(\s*t\.alwaysOn\s*\)\s*return\s+true/.test(rpBody));

ok('always category respects hasData',
   /cat\s*===\s*'always'/.test(rpBody) && /t\.hasData/.test(rpBody));

ok('qc category off by default (visible iff in shown)',
   /cat\s*===\s*'qc'[\s\S]{0,200}?view\.shown\.has/.test(rpBody) ||
   /'qc'[\s\S]{0,100}?'popstats'[\s\S]{0,200}?view\.shown\.has/.test(rpBody));

ok('popstats category off by default',
   /cat\s*===\s*'popstats'/.test(rpBody) || /'qc'\s*\|\|\s*cat\s*===\s*'popstats'/.test(rpBody));

ok('legacy auto-discovered tracks default-on when hasData',
   /return\s+!!t\.hasData/.test(rpBody));

ok('chips sorted by category order',
   /_CAT_ORDER/.test(rpBody) &&
   /sortedTracks/.test(rpBody));

ok('chip strip inserts section labels for qc and popstats',
   /_CAT_LABELS/.test(rpBody) &&
   /'QC'/.test(rpBody) &&
   /'popstats'/.test(rpBody));

ok('chip carries data-category attribute',
   /chip\.dataset\.category\s*=\s*cat/.test(rpBody));

ok('chip click toggles shown set for qc/popstats',
   /trackCat\s*===\s*'qc'\s*\|\|\s*trackCat\s*===\s*'popstats'/.test(rpBody) &&
   /cur\.shown\.add/.test(rpBody) &&
   /cur\.shown\.delete/.test(rpBody));

ok('chip click toggles hidden set for legacy/always categories',
   /cur\.hidden\.add/.test(rpBody) &&
   /cur\.hidden\.delete/.test(rpBody));

ok('chip tooltip surfaces "off by default" hint for qc/popstats',
   /off by default/.test(rpBody));

// ============================================================================
// 5. Sandboxed: loadPopstatsView / savePopstatsView round-trip
// ============================================================================
console.log('\n=== 5. Storage round-trip ===');

function makeStorageSandbox() {
  const store = {};
  const ctx = {
    localStorage: {
      getItem(k) { return Object.prototype.hasOwnProperty.call(store, k) ? store[k] : null; },
      setItem(k, v) { store[k] = String(v); },
    },
    JSON, Set, Object, Array, console,
  };
  vm.createContext(ctx);
  // Inject the constant + helpers
  vm.runInContext("const POPSTATS_STORAGE_KEY = 'scrubber_v3_popstats';", ctx);
  vm.runInContext(pullFunction(html, 'loadPopstatsView'), ctx);
  vm.runInContext(pullFunction(html, 'savePopstatsView'), ctx);
  return { ctx, store };
}

// 5a. Empty storage → both sets empty
{
  const { ctx } = makeStorageSandbox();
  const view = ctx.loadPopstatsView();
  ok('empty storage → view.hidden empty', view.hidden instanceof ctx.Set && view.hidden.size === 0);
  ok('empty storage → view.shown empty',  view.shown  instanceof ctx.Set && view.shown.size === 0);
}

// 5b. Round-trip {hidden, shown}
{
  const { ctx, store } = makeStorageSandbox();
  const v = { hidden: new ctx.Set(['z']), shown: new ctx.Set(['fst_hom1_hom2', 'hobs_hexp']) };
  ctx.savePopstatsView(v);
  const persisted = JSON.parse(store['scrubber_v3_popstats']);
  ok('persisted JSON has hiddenChips array', Array.isArray(persisted.hiddenChips));
  ok('persisted JSON has shownChips array', Array.isArray(persisted.shownChips));
  ok('persisted hiddenChips contains z',
     persisted.hiddenChips.includes('z') && persisted.hiddenChips.length === 1);
  ok('persisted shownChips contains both popstats ids',
     persisted.shownChips.includes('fst_hom1_hom2') &&
     persisted.shownChips.includes('hobs_hexp') &&
     persisted.shownChips.length === 2);
  const v2 = ctx.loadPopstatsView();
  ok('round-trip preserves hidden set',
     v2.hidden.has('z') && v2.hidden.size === 1);
  ok('round-trip preserves shown set',
     v2.shown.has('fst_hom1_hom2') && v2.shown.has('hobs_hexp') && v2.shown.size === 2);
}

// 5c. Back-compat: bare Set passed to savePopstatsView treated as hidden
{
  const { ctx, store } = makeStorageSandbox();
  ctx.savePopstatsView(new ctx.Set(['legacy_track']));
  const persisted = JSON.parse(store['scrubber_v3_popstats']);
  ok('bare Set → hidden persisted',
     persisted.hiddenChips && persisted.hiddenChips.includes('legacy_track'));
  ok('bare Set → shown empty',
     Array.isArray(persisted.shownChips) && persisted.shownChips.length === 0);
}

// 5d. Back-compat: legacy storage with only hiddenChips loads correctly
{
  const { ctx, store } = makeStorageSandbox();
  store['scrubber_v3_popstats'] = JSON.stringify({ hiddenChips: ['old_track_a', 'old_track_b'] });
  const view = ctx.loadPopstatsView();
  ok('legacy storage → hidden populated',
     view.hidden.has('old_track_a') && view.hidden.has('old_track_b'));
  ok('legacy storage → shown empty',
     view.shown.size === 0);
}

// 5e. Corrupted JSON in storage → empty sets, no throw
{
  const { ctx, store } = makeStorageSandbox();
  store['scrubber_v3_popstats'] = 'not-json {{{';
  let threw = null;
  let view = null;
  try { view = ctx.loadPopstatsView(); } catch (e) { threw = e; }
  ok('corrupted JSON → no throw', threw === null);
  ok('corrupted JSON → empty hidden + shown sets',
     view.hidden.size === 0 && view.shown.size === 0);
}

// ============================================================================
// 6. Sandboxed: _popstatsCategoryOf
// ============================================================================
console.log('\n=== 6. _popstatsCategoryOf ===');

{
  const ctx = { console };
  vm.createContext(ctx);
  vm.runInContext(pullFunction(html, '_popstatsCategoryOf'), ctx);
  ok('explicit always category → "always"',
     ctx._popstatsCategoryOf({ category: 'always' }) === 'always');
  ok('explicit qc category → "qc"',
     ctx._popstatsCategoryOf({ category: 'qc' }) === 'qc');
  ok('explicit popstats category → "popstats"',
     ctx._popstatsCategoryOf({ category: 'popstats' }) === 'popstats');
  ok('missing category → "other"',
     ctx._popstatsCategoryOf({ id: 'auto_discovered' }) === 'other');
  ok('null track → "other"',
     ctx._popstatsCategoryOf(null) === 'other');
  ok('non-string category → "other"',
     ctx._popstatsCategoryOf({ category: 42 }) === 'other');
}

// ============================================================================
// 7. Existing flow preserved
// ============================================================================
console.log('\n=== 7. Existing flow preserved ===');

ok('renderPopstatsPage still defined',
   /^function\s+renderPopstatsPage\s*\(\s*\)/m.test(html));

ok('window.renderPopstatsPage export still in place',
   /window\.renderPopstatsPage\s*=\s*renderPopstatsPage/.test(html));

ok('collectPopstatsTracks still defined',
   /^function\s+collectPopstatsTracks\s*\(\s*\)/m.test(html));

ok('alwaysOn flag still on ideogram',
   /id:\s*'ideogram'[\s\S]{0,200}?alwaysOn:\s*true/.test(html));

ok('turn 156 V-shape diagnostic still present',
   /function\s+_buildVShapeData\s*\(\s*candidate\s*,\s*chunk\s*,\s*sqRows\s*\)/.test(html));

ok('turn 157 dosage bridge auto-install still wired',
   /_wireDosageBridgeAutoInstall/.test(html));

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
