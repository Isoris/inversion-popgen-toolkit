// =============================================================================
// turn 123 integration test — comparative TE breakpoint fragility (page16b)
//
// Covers:
//   1. Source-level: detector + store + persist + restore + clear,
//      classifier + dispatcher + restore-on-load wiring, CSS classes.
//   2. JSON detection: comparative_te_breakpoint_fragility_v1.
//   3. Behavioural: _msGetTEFragilityForBreakpoint match-by-bp_id and
//      match-by-chrom-position.
//   4. _msBuildTEFragilityStripHtml empty-states + populated table.
//   5. Integration with center column rendering.
//   6. Regression: turn 121 + 122 still wired.
// =============================================================================

const fs = require('fs');
const vm = require('vm');

const html = fs.readFileSync('/home/claude/work/build/Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// =============================================================================
// Source-level
// =============================================================================
console.log('\n=== Source-level: layer wiring ===');

ok('_isCompTEFragilityJSON defined',     /function _isCompTEFragilityJSON\(/.test(html));
ok('_storeCompTEFragility defined',       /function _storeCompTEFragility\(/.test(html));
ok('_persistCompTEFragility defined',     /function _persistCompTEFragility\(/.test(html));
ok('_restoreCompTEFragility defined',     /function _restoreCompTEFragility\(/.test(html));
ok('_clearCompTEFragility defined',       /function _clearCompTEFragility\(/.test(html));
ok('_msGetTEFragilityForBreakpoint defined',
   /function _msGetTEFragilityForBreakpoint\(/.test(html));
ok('_msBuildTEFragilityStripHtml defined',
   /function _msBuildTEFragilityStripHtml\(/.test(html));

ok('classifier routes comp_te_fragility',
   /_isCompTEFragilityJSON\(data\)\)\s*return\s+'comp_te_fragility'/.test(html));
ok('loadMultipleJSONs dispatches comp_te_fragility',
   /_isCompTEFragilityJSON[\s\S]{0,400}_storeCompTEFragility/.test(html));
ok('restore wired on page load',
   /_restoreCompTEFragility/.test(html) && /_restoreCompTEFragility\(\)/.test(html));
ok('localStorage key for teFragility',
   html.indexOf("'inversion_atlas.teFragility.v1'") >= 0);

ok('window exposes _isCompTEFragilityJSON',
   /window\._isCompTEFragilityJSON\s*=/.test(html));
ok('window exposes _msGetTEFragilityForBreakpoint',
   /window\._msGetTEFragilityForBreakpoint\s*=/.test(html));
ok('window exposes _msBuildTEFragilityStripHtml',
   /window\._msBuildTEFragilityStripHtml\s*=/.test(html));

console.log('\n=== Source-level: center column integration ===');
ok('center column calls _msBuildTEFragilityStripHtml (lineage-loaded path)',
   html.indexOf('_msBuildTEFragilityStripHtml(bp)') >= 0);
// Strip should appear in BOTH branches: when lineageDist exists AND when it doesn't.
const stripCalls = (html.match(/_msBuildTEFragilityStripHtml\(bp\)/g) || []).length;
ok('strip rendered in both lineage-loaded and lineage-empty branches (≥2 calls)',
   stripCalls >= 2, 'found ' + stripCalls + ' call(s)');

console.log('\n=== CSS contract ===');
const cssClasses = [
  'ms-tef-table', 'ms-tef-empty',
  'ms-tef-fold-chip',
  'ms-tef-fold-depleted', 'ms-tef-fold-neutral',
  'ms-tef-fold-elevated', 'ms-tef-fold-strong',
  'ms-tef-legend',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined', html.indexOf('.' + cls) >= 0);
}

console.log('\n=== Regression: prior turns still wired ===');
ok('turn 121 _isSyntenyMultispeciesJSON still defined',
   /function _isSyntenyMultispeciesJSON\(/.test(html));
ok('turn 121 _renderMultiSpeciesPage still defined',
   /function _renderMultiSpeciesPage\(/.test(html));
ok('turn 122 _msAutoSuggestAgeModel still defined',
   /function _msAutoSuggestAgeModel\(/.test(html));
ok('turn 122 _msBuildClassificationTSV still defined',
   /function _msBuildClassificationTSV\(/.test(html));
ok('turn 121 page16b tab button still in tabBar',
   /data-page="page16b" data-stage="synthesis"/.test(html));

// =============================================================================
// Behavioural — pull functions into a sandbox
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '/' && i + 1 < src.length && src[i+1] === '/') {
      i = src.indexOf('\n', i);
      if (i < 0) i = src.length;
      continue;
    }
    else if (ch === '/' && i + 1 < src.length && src[i+1] === '*') {
      i = src.indexOf('*/', i + 2);
      if (i < 0) i = src.length;
      else i += 2;
      continue;
    }
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        // Bail on newlines for quote/apostrophe — broken string literal
        if ((quote === '"' || quote === "'") && src[i] === '\n') break;
        i++;
      }
    }
    i++;
  }
  return src.substring(start, i);
}
function pullConst(src, name) {
  const re = new RegExp('^const\\s+' + name + '\\s*=\\s*(\\[[\\s\\S]*?\\]);', 'm');
  const m = src.match(re);
  return m ? m[1] : null;
}

const fnIs       = pullFunction(html, '_isCompTEFragilityJSON');
const fnStore    = pullFunction(html, '_storeCompTEFragility');
const fnLookup   = pullFunction(html, '_msGetTEFragilityForBreakpoint');
const fnBuild    = pullFunction(html, '_msBuildTEFragilityStripHtml');
const fnGetSp    = pullFunction(html, '_msGetEffectiveSpeciesList');
const defaultSp  = pullConst(html, '_MS_DEFAULT_SPECIES');

ok('extracted _isCompTEFragilityJSON source',     !!fnIs);
ok('extracted _storeCompTEFragility source',      !!fnStore);
ok('extracted _msGetTEFragilityForBreakpoint',    !!fnLookup);
ok('extracted _msBuildTEFragilityStripHtml',      !!fnBuild);
ok('extracted _msGetEffectiveSpeciesList',        !!fnGetSp);
ok('extracted _MS_DEFAULT_SPECIES',               !!defaultSp);

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  const state = { teFragility: null, syntenyMultispecies: null };
  const _MS_DEFAULT_SPECIES = ${defaultSp || '[]'};
  ${fnIs}
  ${fnStore}
  ${fnGetSp}
  ${fnLookup}
  ${fnBuild}
  module.exports = {
    _isCompTEFragilityJSON,
    _storeCompTEFragility,
    _msGetEffectiveSpeciesList,
    _msGetTEFragilityForBreakpoint,
    _msBuildTEFragilityStripHtml,
    state,
  };
`;
const ctx = { module: { exports: {} } };
vm.createContext(ctx);
try {
  vm.runInContext(sandboxSrc, ctx);
} catch (e) {
  console.log('  FATAL sandbox load failed:', e.message);
  process.exit(1);
}
const m = ctx.module.exports;

// === JSON detection
console.log('\n=== JSON detection ===');
const validV1 = {
  tool: 'comparative_te_breakpoint_fragility_v1',
  schema_version: 1,
  generated_at: '2026-05-03T00:00:00Z',
  per_breakpoint_per_species: [
    {
      bp_id: 'csbp_lg28_15_1Mb',
      species: 'Cmac',
      gar_chr: 'C_gar_LG28',
      gar_pos_bp: 15115000,
      homologous_chrom: 'C_mac_LG01',
      focal_lo_bp: 14000000, focal_hi_bp: 16200000,
      focal_te_density: 0.612,
      bg_te_density_chrom: 0.318,
      fold_enrichment: 1.92,
      percentile: 96
    },
  ],
};
ok('valid v1 detected',                      m._isCompTEFragilityJSON(validV1) === true);
ok('legacy alias te_fragility_v1 detected',
   m._isCompTEFragilityJSON({ ...validV1, tool: 'te_fragility_v1' }) === true);
ok('wrong tool not detected',
   m._isCompTEFragilityJSON({ ...validV1, tool: 'something_else' }) === false);
ok('missing per_breakpoint_per_species not detected',
   m._isCompTEFragilityJSON({ tool: 'comparative_te_breakpoint_fragility_v1', schema_version: 1 }) === false);
ok('null not detected',                      m._isCompTEFragilityJSON(null) === false);
ok('non-object not detected',                m._isCompTEFragilityJSON('hi') === false);

// === store + lookup
console.log('\n=== Store + lookup ===');
const stored = m._storeCompTEFragility(validV1);
ok('store returns true on valid JSON',  stored === true);
ok('state.teFragility populated',        m.state.teFragility != null);
ok('store deep-clones (mutation safety)',
   m.state.teFragility.per_breakpoint_per_species[0] !== validV1.per_breakpoint_per_species[0]);

// Lookup by bp_id
const bp1 = { id: 'csbp_lg28_15_1Mb', gar_chr: 'C_gar_LG28', gar_pos_mb: 15.115 };
const lookupCmac = m._msGetTEFragilityForBreakpoint(bp1, 'Cmac');
ok('lookup by bp_id returns entry',     lookupCmac != null);
ok('lookup correct species',            lookupCmac && lookupCmac.species === 'Cmac');
ok('lookup carries fold_enrichment',    lookupCmac && lookupCmac.fold_enrichment === 1.92);

// Lookup by chrom + position (within 100 kb)
const bp2 = { id: 'OTHER_ID', gar_chr: 'C_gar_LG28', gar_pos_mb: 15.118 };
const lookupByPos = m._msGetTEFragilityForBreakpoint(bp2, 'Cmac');
ok('lookup by chrom+pos within 100kb',  lookupByPos != null);

// Lookup misses
const bp3 = { id: 'X', gar_chr: 'C_gar_LG14', gar_pos_mb: 5.0 };
ok('wrong chrom → null',                 m._msGetTEFragilityForBreakpoint(bp3, 'Cmac') === null);
const bp4 = { id: 'X', gar_chr: 'C_gar_LG28', gar_pos_mb: 20.0 };  // 5 Mb away
ok('chrom match but >100kb away → null', m._msGetTEFragilityForBreakpoint(bp4, 'Cmac') === null);

// Lookup wrong species
ok('wrong species → null', m._msGetTEFragilityForBreakpoint(bp1, 'Cfus') === null);

// Lookup with no data
m.state.teFragility = null;
ok('no teFragility → null',  m._msGetTEFragilityForBreakpoint(bp1, 'Cmac') === null);
ok('no bp → null',           m._msGetTEFragilityForBreakpoint(null, 'Cmac') === null);
ok('no species → null',      m._msGetTEFragilityForBreakpoint(bp1, null) === null);

// === Build strip HTML
console.log('\n=== Build strip HTML ===');
m.state.teFragility = null;
const html_empty = m._msBuildTEFragilityStripHtml(bp1);
ok('no JSON → empty-state hint with filename',
   /comparative_te_breakpoint_fragility_v1\.json/.test(html_empty));
ok('no JSON → ms-cls-block wrapper',  /class="ms-cls-block"/.test(html_empty));
ok('no JSON → ms-cls-empty span',     /class="ms-cls-empty"/.test(html_empty));
ok('no bp → empty string',            m._msBuildTEFragilityStripHtml(null) === '');

// Now reload and render with data
m._storeCompTEFragility(validV1);
const html_full = m._msBuildTEFragilityStripHtml(bp1);
ok('with data → renders ms-tef-table', /class="ms-tef-table"/.test(html_full));
ok('with data → contains Cmac fold value', /1\.92.*\u00D7|1\.92×/.test(html_full));
ok('with data → contains percentile',  /96.{0,2}%ile|96%ile|96\u00B0/.test(html_full));
ok('with data → contains density %',   /61\.2%/.test(html_full));
ok('with data → ms-tef-fold-elevated chip used (fold ≥ 1.5)',
   /ms-tef-fold-elevated/.test(html_full));
ok('with data → "no data" rows for missing species',
   /class="ms-tef-empty"/.test(html_full));
ok('with data → focal stars on Cgar/Cmac',
   /\u2605/.test(html_full));
ok('with data → contains legend section',
   /class="ms-tef-legend"/.test(html_full));
ok('with data → legend mentions fragility proxy disclaimer',
   /[Ff]ragility proxy/.test(html_full));

// === fold chip threshold scoring
console.log('\n=== Fold chip thresholds ===');
const validBands = {
  tool: 'comparative_te_breakpoint_fragility_v1',
  schema_version: 1,
  per_breakpoint_per_species: [
    { species: 'Cmac',  bp_id: 'b1', fold_enrichment: 0.7,  focal_te_density: 0.20, percentile: 30 },
    { species: 'Cfus',  bp_id: 'b1', fold_enrichment: 1.05, focal_te_density: 0.32, percentile: 55 },
    { species: 'Capus', bp_id: 'b1', fold_enrichment: 1.80, focal_te_density: 0.45, percentile: 88 },
    { species: 'Phyp',  bp_id: 'b1', fold_enrichment: 3.10, focal_te_density: 0.71, percentile: 99 },
  ],
};
m._storeCompTEFragility(validBands);
const html_bands = m._msBuildTEFragilityStripHtml({ id: 'b1', gar_chr: 'C_gar_LG01', gar_pos_mb: 5.0 });
ok('fold 0.7 → depleted band',  /ms-tef-fold-depleted">[^<]*0\.70/.test(html_bands));
ok('fold 1.05 → neutral band',  /ms-tef-fold-neutral">[^<]*1\.05/.test(html_bands));
ok('fold 1.80 → elevated band', /ms-tef-fold-elevated">[^<]*1\.80/.test(html_bands));
ok('fold 3.10 → strong band',   /ms-tef-fold-strong">[^<]*3\.10/.test(html_bands));

// no-matching-bp case
m._storeCompTEFragility(validBands);
const html_nomatch = m._msBuildTEFragilityStripHtml({ id: 'unknown_bp', gar_chr: 'C_gar_LG99', gar_pos_mb: 99 });
ok('JSON loaded but no bp match → friendly empty-state',
   /no entry matching this breakpoint/.test(html_nomatch));
ok('no-match empty-state cites breakpoint id', /unknown_bp/.test(html_nomatch));

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
