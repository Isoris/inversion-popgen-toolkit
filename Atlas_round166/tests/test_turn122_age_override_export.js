// =============================================================================
// turn 122 integration test — age-model auto-suggest, dxy_per_inversion layer,
// manual classification override, manuscript-export TSV, TE-fragility panel.
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
// Source-level checks
// =============================================================================
console.log('\n=== Source: dxy_per_inversion layer ===');
ok('_isDxyPerInversionJSON defined',  /function _isDxyPerInversionJSON\(/.test(html));
ok('_storeDxyPerInversion defined',   /function _storeDxyPerInversion\(/.test(html));
ok('_persistDxyPerInversion defined', /function _persistDxyPerInversion\(/.test(html));
ok('_restoreDxyPerInversion defined', /function _restoreDxyPerInversion\(/.test(html));
ok('_clearDxyPerInversion defined',   /function _clearDxyPerInversion\(/.test(html));
ok('_msGetDxyForBreakpoint defined',  /function _msGetDxyForBreakpoint\(/.test(html));
ok('classifier routes dxy_per_inversion',
   /_isDxyPerInversionJSON\(data\)\)\s*return\s+'dxy_per_inversion'/.test(html));
ok('loadMultipleJSONs dispatches dxy_per_inversion',
   /_isDxyPerInversionJSON[\s\S]{0,400}_storeDxyPerInversion/.test(html));
ok('restore wired on load (dxy)', /_restoreDxyPerInversion/.test(html));
ok('localStorage key for dxyPerInversion',
   html.indexOf("'inversion_atlas.dxyPerInversion.v1'") >= 0);

console.log('\n=== Source: age-model auto-suggest ===');
ok('_msAutoSuggestAgeModel defined',  /function _msAutoSuggestAgeModel\(/.test(html));
ok('YOUNG-POP rule emitted in source',  /age_model:\s*'YOUNG-POP'/.test(html));
ok('OLD-POLY rule emitted in source',   /age_model:\s*'OLD-POLY'/.test(html));
ok('OLD-BP-YOUNG-INV rule in source',   /age_model:\s*'OLD-BP-YOUNG-INV'/.test(html));
ok('LINEAGE-KARYO rule in source',      /age_model:\s*'LINEAGE-KARYO'/.test(html));
ok('MULTI-AGE-HOTSPOT rule in source',  /age_model:\s*'MULTI-AGE-HOTSPOT'/.test(html));

console.log('\n=== Source: manual override editor ===');
ok('_msInitClassifications defined', /function _msInitClassifications\(/.test(html));
ok('_msSetClassification defined',   /function _msSetClassification\(/.test(html));
ok('_msGetClassification defined',   /function _msGetClassification\(/.test(html));
ok('_msClearClassification defined', /function _msClearClassification\(/.test(html));
ok('_msPersistClassifications defined', /function _msPersistClassifications\(/.test(html));
ok('localStorage key for classifications',
   html.indexOf("'inversion_atlas.classifications.v1'") >= 0);
ok('override-save button wired',  /id="msOverrideSaveBtn"/.test(html));
ok('override-clear button wired', /id="msOverrideClearBtn"/.test(html));
ok('override-arch select wired',  /id="msOverrideArch"/.test(html));
ok('override-age select wired',   /id="msOverrideAge"/.test(html));
ok('override-conf select wired',  /id="msOverrideConf"/.test(html));
ok('override-notes textarea wired',/id="msOverrideNotes"/.test(html));

console.log('\n=== Source: TSV export ===');
ok('_msBuildClassificationTSV defined',     /function _msBuildClassificationTSV\(/.test(html));
ok('_msDownloadClassificationTSV defined',  /function _msDownloadClassificationTSV\(/.test(html));
ok('export button id wired',                 /id="msExportTSVBtn"/.test(html));
ok('TSV header includes architecture_class', /architecture_class/.test(html));
ok('TSV header includes age_model',          /'age_model'/.test(html));
ok('TSV header includes architecture_source',/'architecture_source'/.test(html));

console.log('\n=== Source: focal-vs-bg in detail column ===');
ok('msDetailFocalVsBg DOM slot',          /id="msDetailFocalVsBg"/.test(html));
ok('focal-vs-bg call site in detail',     /makeFocalVsBgPanel[\s\S]{0,500}'page16b'/.test(html));
ok('focal-vs-bg shares state._focalVsBg', /state\._focalVsBg\.csRadiusBp/.test(html));

console.log('\n=== Source: CSS for new pieces ===');
const cssClasses = [
  'ms-age-YOUNG_POP', 'ms-age-OLD_POLY', 'ms-age-OLD_BP_YOUNG_INV',
  'ms-age-LINEAGE_KARYO', 'ms-age-MULTI_AGE_HOTSPOT',
  'ms-cls-override', 'ms-cls-override-grid', 'ms-cls-override-buttons',
  'ms-cls-select', 'ms-cls-textarea',
  'ms-btn', 'ms-btn-primary',
  'ms-cls-export', 'ms-cls-export-hint', 'ms-cls-dxy-note',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls, html.indexOf('.' + cls) >= 0);
}

console.log('\n=== Source: window exports ===');
const windowExports = [
  '_isDxyPerInversionJSON', '_storeDxyPerInversion',
  '_msAutoSuggestAgeModel',
  '_msSetClassification', '_msGetClassification', '_msClearClassification',
  '_msBuildClassificationTSV', '_msDownloadClassificationTSV',
];
for (const name of windowExports) {
  ok('window.' + name, new RegExp('window\\.' + name + '\\s*=').test(html));
}

// =============================================================================
// Behavioural — sandbox _msAutoSuggestAgeModel
// =============================================================================
console.log('\n=== Behavioural — age-model auto-suggest ===');

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
    if (ch === '/' && src[i+1] === '/') {
      const nl = src.indexOf('\n', i);
      i = nl < 0 ? src.length : nl + 1;
      continue;
    }
    if (ch === '/' && src[i+1] === '*') {
      const close = src.indexOf('*/', i + 2);
      i = close < 0 ? src.length : close + 2;
      continue;
    }
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        if (src[i] === '\n' && (quote === '"' || quote === "'")) break;
        i++;
      }
    }
    i++;
  }
  return src.substring(start, i);
}

const fnAge = pullFunction(html, '_msAutoSuggestAgeModel');
ok('extracted _msAutoSuggestAgeModel source', !!fnAge);

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  ${fnAge}
  module.exports = { _msAutoSuggestAgeModel };
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

// MULTI-AGE-HOTSPOT — strongest case, mixed types, 3+ species
const ml_multiAge = m._msAutoSuggestAgeModel(
  { id: 'bp1', event_type: 'inversion' },
  { Cgar: 'boundary_present', Cmac: 'boundary_present',
    Phyp: 'fission_in_target', Tfulv: 'fusion_in_target',
    Cfus: 'unknown' },
  null
);
ok('mixed types + 4 events → MULTI-AGE-HOTSPOT',
   ml_multiAge.age_model === 'MULTI-AGE-HOTSPOT');
ok('MULTI-AGE-HOTSPOT confidence is medium',
   ml_multiAge.confidence === 'medium');

// LINEAGE-KARYO — fission event_type + sister-species fission
const ml_lineage = m._msAutoSuggestAgeModel(
  { id: 'bp2', event_type: 'fission' },
  { Cgar: 'boundary_present', Phyp: 'fission_in_target' },
  null
);
ok('fission event + comparative fission → LINEAGE-KARYO',
   ml_lineage.age_model === 'LINEAGE-KARYO');

// OLD-BP-YOUNG-INV — 3+ species share boundary, dxy not elevated
const ml_oldbpYoung = m._msAutoSuggestAgeModel(
  { id: 'bp3', event_type: 'inversion' },
  { Cgar: 'boundary_present', Cmac: 'boundary_present',
    Cfus: 'boundary_present', Capus: 'boundary_present',
    Phyp: 'unknown' },
  { fold_elevation_inside_vs_flank: 1.1, dxy_within_inversion_ref_vs_inv: 0.0015 }
);
ok('4 present + low fold (1.1x) → OLD-BP-YOUNG-INV',
   ml_oldbpYoung.age_model === 'OLD-BP-YOUNG-INV');
ok('OLD-BP-YOUNG-INV with dxy data has medium confidence',
   ml_oldbpYoung.confidence === 'medium');

// OLD-POLY — high fold elevation + sharing
const ml_oldPoly = m._msAutoSuggestAgeModel(
  { id: 'bp4', event_type: 'inversion' },
  { Cgar: 'boundary_present', Cmac: 'boundary_present' },
  { fold_elevation_inside_vs_flank: 2.1, dxy_within_inversion_ref_vs_inv: 0.0042 }
);
ok('2 present + high fold (2.1x) → OLD-POLY',
   ml_oldPoly.age_model === 'OLD-POLY');

// YOUNG-POP — low fold, only focal species
const ml_young = m._msAutoSuggestAgeModel(
  { id: 'bp5', event_type: 'inversion' },
  { Cgar: 'boundary_present', Cmac: 'boundary_absent_internal_to_block',
    Cfus: 'boundary_absent_internal_to_block' },
  { fold_elevation_inside_vs_flank: 1.05, dxy_within_inversion_ref_vs_inv: 0.0008 }
);
ok('1 present + low fold (1.05x) → YOUNG-POP',
   ml_young.age_model === 'YOUNG-POP');

// Uncertain — no dxy, no lineage
const ml_uncertain = m._msAutoSuggestAgeModel(
  { id: 'bp6', event_type: 'inversion' },
  null, null
);
ok('no dxy + no lineage → null age_model',
   ml_uncertain.age_model === null);
ok('uncertain rationale mentions dxy_per_inversion',
   /dxy_per_inversion_v1/.test(ml_uncertain.rationale));

// Uncertain — has lineage but no dxy
const ml_noDxy = m._msAutoSuggestAgeModel(
  { id: 'bp7', event_type: 'inversion' },
  { Cgar: 'boundary_present' },
  null
);
ok('lineage but no dxy → uncertain or rule-fallback',
   ml_noDxy.age_model === null || ml_noDxy.age_model === 'YOUNG-POP' ||
   ml_noDxy.age_model === 'OLD-BP-YOUNG-INV');

// No bp at all
const ml_noBp = m._msAutoSuggestAgeModel(null, null, null);
ok('no bp → null age_model + unknown confidence',
   ml_noBp.age_model === null && ml_noBp.confidence === 'unknown');

// Signals object always populated
ok('signals object always returned',
   ml_multiAge.signals && typeof ml_multiAge.signals.n_present === 'number');
ok('signals.n_rearrangement_types reflects mixed types',
   ml_multiAge.signals.n_rearrangement_types >= 2);

// =============================================================================
// JSON detector tests
// =============================================================================
console.log('\n=== JSON detectors ===');

const fnIsDxy = pullFunction(html, '_isDxyPerInversionJSON');
ok('extracted _isDxyPerInversionJSON', !!fnIsDxy);
const ctx2 = { module: { exports: {} } };
vm.createContext(ctx2);
vm.runInContext(fnIsDxy + '\nmodule.exports = { _isDxyPerInversionJSON };', ctx2);
const m2 = ctx2.module.exports;

const dxyValid = {
  tool: 'dxy_per_inversion_v1',
  schema_version: 1,
  per_inversion: [{ candidate_id: 'cand_lg28', dxy_within_inversion_ref_vs_inv: 0.004 }],
};
ok('valid dxy JSON detected', m2._isDxyPerInversionJSON(dxyValid) === true);
ok('wrong tool not detected',
   m2._isDxyPerInversionJSON({ ...dxyValid, tool: 'something_else' }) === false);
ok('missing per_inversion not detected',
   m2._isDxyPerInversionJSON({ tool: 'dxy_per_inversion_v1', schema_version: 1 }) === false);
ok('null not detected', m2._isDxyPerInversionJSON(null) === false);

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
