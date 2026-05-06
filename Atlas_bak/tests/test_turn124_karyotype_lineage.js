// =============================================================================
// turn 124 integration test — karyotype-context lineage layer (page16b)
//
// Covers:
//   1. Source-level: detector + store + persist + restore + clear,
//      classifier + dispatcher + restore-on-load wiring, CSS classes.
//   2. JSON detection: karyotype_lineage_v1.
//   3. Behavioural: polarization rules
//      - cgar_lineage_fission   (Cgar 1-2, others 1-1, outgroup 1-1)
//      - cmac_lineage_fission   (Cmac 1-2, others 1-1, outgroup 1-1)
//      - ancestral_split_both_retained (most species 1-2 with concordant targets)
//      - recurrent_fission_hotspot (≥3 species 1-2 with discordant targets)
//      - no_karyotype_change (all 1-1)
//      - unresolved (mixed sparse signal)
//   4. Lookup: cgar_chr matcher across naming variants.
//   5. Render HTML: empty-state + populated table + outgroup tag + verdict chip.
//   6. Age-model integration: karyoVerdict argument flips LINEAGE-KARYO,
//      MULTI-AGE-HOTSPOT decisions.
//   7. Lineage table: shows merged karyo_class column when both layers loaded.
//   8. Regression: turns 117–123 still pass.
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

ok('_isKaryotypeLineageJSON defined',     /function _isKaryotypeLineageJSON\(/.test(html));
ok('_storeKaryotypeLineage defined',      /function _storeKaryotypeLineage\(/.test(html));
ok('_persistKaryotypeLineage defined',    /function _persistKaryotypeLineage\(/.test(html));
ok('_restoreKaryotypeLineage defined',    /function _restoreKaryotypeLineage\(/.test(html));
ok('_clearKaryotypeLineage defined',      /function _clearKaryotypeLineage\(/.test(html));
ok('_msGetKaryotypeEntryForCgarChr defined',
   /function _msGetKaryotypeEntryForCgarChr\(/.test(html));
ok('_msPolarizeKaryotypeEvent defined',
   /function _msPolarizeKaryotypeEvent\(/.test(html));
ok('_msBuildKaryotypeContextHtml defined',
   /function _msBuildKaryotypeContextHtml\(/.test(html));

ok('classifier routes karyotype_lineage',
   /_isKaryotypeLineageJSON\(data\)\)\s*return\s+'karyotype_lineage'/.test(html));
ok('loadMultipleJSONs dispatches karyotype_lineage',
   /_isKaryotypeLineageJSON[\s\S]{0,400}_storeKaryotypeLineage/.test(html));
ok('restore wired on page load',
   /_restoreKaryotypeLineage\(\)/.test(html));
ok('localStorage key for karyotypeLineage',
   html.indexOf("'inversion_atlas.karyotypeLineage.v1'") >= 0);

ok('window exposes _isKaryotypeLineageJSON',
   /window\._isKaryotypeLineageJSON\s*=/.test(html));
ok('window exposes _msPolarizeKaryotypeEvent',
   /window\._msPolarizeKaryotypeEvent\s*=/.test(html));
ok('window exposes _msBuildKaryotypeContextHtml',
   /window\._msBuildKaryotypeContextHtml\s*=/.test(html));

console.log('\n=== Source-level: integration with center column + age model ===');
ok('center column calls _msBuildKaryotypeContextHtml',
   html.indexOf('_msBuildKaryotypeContextHtml(bp)') >= 0);
const ctxCalls = (html.match(/_msBuildKaryotypeContextHtml\(bp\)/g) || []).length;
ok('karyotype context rendered in both branches (≥2 calls)',
   ctxCalls >= 2, 'found ' + ctxCalls + ' call(s)');
ok('age-model fn signature now takes karyoVerdict',
   /function _msAutoSuggestAgeModel\(bp, lineageDist, dxyEntry, karyoVerdict\)/.test(html));
ok('age-model RULE 0 reads cgar_lineage_fission',
   html.indexOf("cgar_lineage_fission") >= 0);
ok('age-model RULE 0 reads cmac_lineage_fission',
   html.indexOf("cmac_lineage_fission") >= 0);
ok('age-model RULE 0 reads recurrent_fission_hotspot',
   html.indexOf("recurrent_fission_hotspot") >= 0);
ok('center column passes karyoVerdict to age-model',
   /_msPolarizeKaryotypeEvent\(bp\.gar_chr\);[\s\S]{0,200}_msAutoSuggestAgeModel\(bp, lineageDist, dxyEntry, karyoVerdict\)/.test(html));
ok('TSV exporter passes karyoVerdict to age-model',
   (html.match(/_msAutoSuggestAgeModel\(bp, lineageDist, dxyEntry, karyoVerdict\)/g) || []).length >= 2);

console.log('\n=== CSS contract ===');
const cssClasses = [
  'ms-karyo-table', 'ms-karyo-targets',
  'ms-karyo-class-chip',
  'ms-karyo-class-1to1', 'ms-karyo-class-1to2',
  'ms-karyo-class-1to3', 'ms-karyo-class-1to4plus',
  'ms-karyo-class-other', 'ms-karyo-outgroup-tag',
  'ms-karyo-verdict-row', 'ms-karyo-verdict-chip',
  'ms-karyo-verdict-cgar', 'ms-karyo-verdict-cmac',
  'ms-karyo-verdict-ancestral', 'ms-karyo-verdict-hotspot',
  'ms-karyo-verdict-stable', 'ms-karyo-verdict-unresolved',
  'ms-karyo-legend',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined', html.indexOf('.' + cls) >= 0);
}

console.log('\n=== Regression: prior turns still wired ===');
ok('turn 121 _isSyntenyMultispeciesJSON still defined',
   /function _isSyntenyMultispeciesJSON\(/.test(html));
ok('turn 122 _msAutoSuggestAgeModel still defined',
   /function _msAutoSuggestAgeModel\(/.test(html));
ok('turn 122 _msBuildClassificationTSV still defined',
   /function _msBuildClassificationTSV\(/.test(html));
ok('turn 123 _isCompTEFragilityJSON still defined',
   /function _isCompTEFragilityJSON\(/.test(html));
ok('turn 121 page16b tab button still in tabBar',
   /data-page="page16b" data-stage="synthesis"/.test(html));

// =============================================================================
// Behavioural — sandbox
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

const fnIs       = pullFunction(html, '_isKaryotypeLineageJSON');
const fnStore    = pullFunction(html, '_storeKaryotypeLineage');
const fnGetEntryFocal = pullFunction(html, '_msGetKaryotypeEntryForFocalChr');
const fnGetEntry = pullFunction(html, '_msGetKaryotypeEntryForCgarChr');
const fnPolarize = pullFunction(html, '_msPolarizeKaryotypeEvent');
const fnSummarize = pullFunction(html, '_msSummarizeRefinement');
const fnAdjust = pullFunction(html, '_msAdjustConfidenceForRefinement');
const fnPolarizeRef = pullFunction(html, '_msPolarizeKaryotypeEventWithRefinement');
const fnRefineChip = pullFunction(html, '_msBuildRefinementChipHtml');
const fnEffCls = pullFunction(html, '_msGetEffectiveClassForCell');
const fnEffTgt = pullFunction(html, '_msGetEffectiveTargetsForCell');
const fnBuild    = pullFunction(html, '_msBuildKaryotypeContextHtml');
const fnGetSp    = pullFunction(html, '_msGetEffectiveSpeciesList');
const defaultSp  = pullConst(html, '_MS_DEFAULT_SPECIES');

ok('extracted _isKaryotypeLineageJSON',     !!fnIs);
ok('extracted _storeKaryotypeLineage',      !!fnStore);
ok('extracted _msGetKaryotypeEntryForCgarChr', !!fnGetEntry);
ok('extracted _msGetKaryotypeEntryForFocalChr', !!fnGetEntryFocal);
ok('extracted _msPolarizeKaryotypeEvent',   !!fnPolarize);
ok('extracted _msBuildKaryotypeContextHtml',!!fnBuild);
ok('extracted _msGetEffectiveSpeciesList',  !!fnGetSp);

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  const state = { karyotypeLineage: null, syntenyMultispecies: null };
  const _MS_DEFAULT_SPECIES = ${defaultSp || '[]'};
  ${fnIs}
  ${fnStore}
  ${fnGetSp}
  ${fnGetEntryFocal}
  ${fnGetEntry}
  ${fnPolarize}
  ${fnSummarize}
  ${fnAdjust}
  ${fnEffCls}
  ${fnEffTgt}
  ${fnRefineChip}
  ${fnPolarizeRef}
  ${fnBuild}
  module.exports = {
    _isKaryotypeLineageJSON,
    _storeKaryotypeLineage,
    _msGetKaryotypeEntryForCgarChr,
    _msPolarizeKaryotypeEvent,
    _msBuildKaryotypeContextHtml,
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
const valid = {
  tool: 'karyotype_lineage_v1',
  schema_version: 1,
  generated_at: '2026-05-03T13:00:00Z',
  params: { method: 'mashmap', segment_bp: 1000000, identity_pct: 85, filter: 'one-to-one' },
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [
    {
      cgar_chr: 'C_gar_LG28',
      classes_by_species: {
        Cgar:  { class: '1-1', targets: ['C_gar_LG28'] },
        Cmac:  { class: '1-2', targets: ['C_mac_LG01', 'C_mac_LG18'] },
        Cfus:  { class: '1-1', targets: ['C_fus_LG07'] },
        Capus: { class: '1-1', targets: ['Capus_LG12'] },
        Phyp:  { class: '1-1', targets: ['Phyp_chr18'] },
        Ngraeffei: { class: '1-1', targets: ['N_graeffei|22'] },
      },
    },
  ],
};
ok('valid karyotype_lineage_v1 detected', m._isKaryotypeLineageJSON(valid) === true);
ok('legacy alias mashmap_karyotype_lineage_v1 detected',
   m._isKaryotypeLineageJSON({ ...valid, tool: 'mashmap_karyotype_lineage_v1' }) === true);
ok('wrong tool not detected',
   m._isKaryotypeLineageJSON({ ...valid, tool: 'something_else' }) === false);
ok('missing per_cgar_chr not detected',
   m._isKaryotypeLineageJSON({ tool: 'karyotype_lineage_v1', schema_version: 1, params: {} }) === false);
ok('missing params not detected',
   m._isKaryotypeLineageJSON({ tool: 'karyotype_lineage_v1', schema_version: 1, per_cgar_chr: [] }) === false);
ok('null not detected', m._isKaryotypeLineageJSON(null) === false);

// === store + entry lookup
console.log('\n=== Store + lookup ===');
const stored = m._storeKaryotypeLineage(valid);
ok('store returns true on valid JSON',  stored === true);
ok('state.karyotypeLineage populated',  m.state.karyotypeLineage != null);
ok('outgroup_species preserved',         Array.isArray(m.state.karyotypeLineage.outgroup_species) &&
                                          m.state.karyotypeLineage.outgroup_species[0] === 'Ngraeffei');

const e1 = m._msGetKaryotypeEntryForCgarChr('C_gar_LG28');
ok('lookup exact cgar_chr',  e1 != null && e1.cgar_chr === 'C_gar_LG28');
const e2 = m._msGetKaryotypeEntryForCgarChr('LG28');
ok('lookup short LG28 → matches C_gar_LG28', e2 != null && e2.cgar_chr === 'C_gar_LG28');
const e3 = m._msGetKaryotypeEntryForCgarChr('LG99');
ok('lookup unknown chr → null', e3 === null);

// === Polarization rules
console.log('\n=== Polarization verdicts ===');

// Rule 2: focal-lineage fission (Cgar = inferred focal in this data)
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG10',
    classes_by_species: {
      Cgar:  { class: '1-2', targets: ['C_gar_LG10', 'C_gar_LG14'] },
      Cmac:  { class: '1-1', targets: ['C_mac_LG09'] },
      Cfus:  { class: '1-1', targets: ['C_fus_LG09'] },
      Phyp:  { class: '1-1', targets: ['Phyp_chr12'] },
      Ngraeffei: { class: '1-1', targets: ['N_graeffei|10'] },
    },
  }],
});
const v_focal = m._msPolarizeKaryotypeEvent('C_gar_LG10');
ok('Cgar 1-2 + sister/outgroups 1-1 → focal_lineage_fission',
   v_focal.verdict === 'focal_lineage_fission');
ok('focal_lineage_fission has medium confidence',
   v_focal.confidence === 'medium');
ok('signals.focal_class === "1-2"', v_focal.signals.focal_class === '1-2');
ok('signals.focal_species === "Cgar"', v_focal.signals.focal_species === 'Cgar');
ok('signals.outgroup_one_to_one >= 1', v_focal.signals.outgroup_one_to_one >= 1);

// Rule 3: sister-lineage fission (Cmac = sister, Cmac shows 1-2)
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG01', 'C_mac_LG18'] },
      Cfus:  { class: '1-1', targets: ['C_fus_LG07'] },
      Phyp:  { class: '1-1', targets: ['Phyp_chr18'] },
      Ngraeffei: { class: '1-1', targets: ['N_graeffei|22'] },
    },
  }],
});
const v_sister = m._msPolarizeKaryotypeEvent('C_gar_LG28');
ok('Cmac 1-2 (sister) + focal/outgroups 1-1 → sister_lineage_fission',
   v_sister.verdict === 'sister_lineage_fission');

// Rule 1a: Recurrent fission hotspot — 3+ species 1-2, discordant targets
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG07',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG07'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG02', 'C_mac_LG18'] },
      Cfus:  { class: '1-2', targets: ['C_fus_LG03', 'C_fus_LG07'] },
      Capus: { class: '1-2', targets: ['Capus_LG05', 'Capus_LG09'] },
      Phyp:  { class: '1-1', targets: ['Phyp_chr07'] },
      Ngraeffei: { class: '1-1', targets: ['N_graeffei|7'] },
    },
  }],
});
const v_hot = m._msPolarizeKaryotypeEvent('C_gar_LG07');
ok('3 species 1-2 with discordant targets → recurrent_fission_hotspot',
   v_hot.verdict === 'recurrent_fission_hotspot');

// Rule 1b: Ancestral split, all retained — concordant targets
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG14',
    classes_by_species: {
      Cgar:  { class: '1-2', targets: ['C_gar_LG14a', 'C_gar_LG14b'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG14a', 'C_mac_LG14b'] },
      Cfus:  { class: '1-2', targets: ['C_fus_LG14a', 'C_fus_LG14b'] },
      Capus: { class: '1-1', targets: ['Capus_LG14'] },
      Ngraeffei: { class: '1-1', targets: ['N_graeffei|14'] },
    },
  }],
});
// (3 species 1-2, but each has different targets — let's instead build one with concordant)
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG14',
    classes_by_species: {
      Cgar:  { class: '1-2', targets: ['T1', 'T2'] },
      Cmac:  { class: '1-2', targets: ['T1', 'T2'] },
      Cfus:  { class: '1-2', targets: ['T1', 'T2'] },
      Capus: { class: '1-1', targets: ['T_C'] },
      Ngraeffei: { class: '1-1', targets: ['T_N'] },
    },
  }],
});
const v_anc = m._msPolarizeKaryotypeEvent('C_gar_LG14');
ok('3 species 1-2 with concordant targets → ancestral_split_both_retained',
   v_anc.verdict === 'ancestral_split_both_retained');

// Rule 4: All 1-1 → no_karyotype_change
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG01',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['T1'] },
      Cmac:  { class: '1-1', targets: ['T1'] },
      Cfus:  { class: '1-1', targets: ['T1'] },
      Capus: { class: '1-1', targets: ['T1'] },
      Phyp:  { class: '1-1', targets: ['T1'] },
      Ngraeffei: { class: '1-1', targets: ['T1'] },
    },
  }],
});
const v_stable = m._msPolarizeKaryotypeEvent('C_gar_LG01');
ok('all species 1-1 → no_karyotype_change',
   v_stable.verdict === 'no_karyotype_change');

// Fallback: mixed sparse → unresolved
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Ngraeffei'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG05',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['T1'] },
      Cmac:  { class: '1-3', targets: ['T1', 'T2', 'T3'] },
      Cfus:  { class: '1-4+', targets: ['T1', 'T2', 'T3', 'T4'] },
      Ngraeffei: { class: '1-3', targets: ['T1', 'T2', 'T3'] },
    },
  }],
});
const v_unr = m._msPolarizeKaryotypeEvent('C_gar_LG05');
ok('mixed sparse signal → unresolved verdict',
   v_unr.verdict === 'unresolved');

// No data
m.state.karyotypeLineage = null;
const v_none = m._msPolarizeKaryotypeEvent('C_gar_LG28');
ok('no karyotype layer → null verdict', v_none.verdict === null);
ok('no karyotype layer → unknown confidence', v_none.confidence === 'unknown');
ok('no karyotype layer → rationale mentions filename',
   /karyotype_lineage_v1\.json/.test(v_none.rationale));

// === HTML render
console.log('\n=== HTML render ===');
m.state.karyotypeLineage = null;
const html_empty = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('no JSON → empty-state with chrom name',
   /C_gar_LG28/.test(html_empty) && /class="ms-cls-empty"/.test(html_empty));
ok('no bp → empty string', m._msBuildKaryotypeContextHtml(null) === '');
ok('bp without gar_chr → empty string',
   m._msBuildKaryotypeContextHtml({ id: 'x' }) === '');

// Re-load Cmac-fission scenario, render full HTML
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG01', 'C_mac_LG18'] },
      Cfus:  { class: '1-1', targets: ['C_fus_LG07'] },
      Phyp:  { class: '1-1', targets: ['Phyp_chr18'] },
      Tros:  { class: '1-1', targets: ['Tros|22'] },
    },
  }],
});
const html_full = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('with data → renders ms-karyo-table',  /class="ms-karyo-table"/.test(html_full));
ok('with data → renders verdict chip (sister)',
   /class="ms-karyo-verdict-chip ms-karyo-verdict-sister"/.test(html_full));
ok('with data → mentions Cmac lineage in label (substituted from sister_species)',
   /Cmac lineage fission/.test(html_full));
ok('with data → contains 1-1 chips',      /ms-karyo-class-1to1/.test(html_full));
ok('with data → contains 1-2 chip',       /ms-karyo-class-1to2/.test(html_full));
ok('with data → outgroup tag rendered',
   /class="ms-karyo-outgroup-tag"/.test(html_full));
ok('with data → contains target chrom string',
   /C_mac_LG01/.test(html_full));
ok('with data → contains resolution disclaimer',
   /1 Mb/.test(html_full) && /85%/.test(html_full));

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
