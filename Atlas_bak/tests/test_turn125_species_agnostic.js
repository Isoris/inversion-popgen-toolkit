// =============================================================================
// turn 125 integration test — species-agnostic karyotype layer
//
// Verifies that the karyotype lineage layer can be driven by a JSON declaring
// arbitrary focal/sister/outgroup species. The atlas should not care whether
// focal is "Cgar" or "Cmac" or some entirely different species.
//
// Covers:
//   1. Source-level: new/renamed schema fields (focal_species, sister_species,
//      per_focal_chr) accepted; legacy (per_cgar_chr) accepted as alias.
//   2. Polarization is species-agnostic:
//      a. Cgar-focal scenario → focal_lineage_fission when Cgar 1-2 + others 1-1
//      b. Cmac-focal scenario (different paper!) → focal_lineage_fission when
//         Cmac 1-2 + others 1-1
//      c. Custom non-catfish-named species → still works
//   3. Sister-species detection: when sister_species is declared, sister 1-2
//      with focal/outgroup 1-1 → sister_lineage_fission.
//   4. Verdict label substitution: "focal_lineage_fission" → "Cgar lineage fission"
//      uses the actual species name from the JSON.
//   5. focal_species inference: when not declared, atlas infers from data.
//   6. Backward compat: legacy verdicts (cgar_/cmac_) still trigger LINEAGE-KARYO
//      in the age-model RULE 0 (so old demo files keep working).
//   7. CSS classes for new verdict variants exist.
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
console.log('\n=== Source-level: schema + helpers ===');

ok('detector accepts per_focal_chr',
   /Array\.isArray\(data\.per_focal_chr\)/.test(html));
ok('detector keeps legacy per_cgar_chr alias',
   /Array\.isArray\(data\.per_cgar_chr\)/.test(html));
ok('store normalizes cgar_chr → focal_chr',
   /copy\.focal_chr\s*=\s*copy\.cgar_chr/.test(html));
ok('store reads focal_species from JSON',
   /parsed\.focal_species/.test(html));
ok('store reads sister_species from JSON',
   /parsed\.sister_species/.test(html));
ok('store has focal-species inference fallback',
   /e\.targets\[0\]\s*===\s*first\.focal_chr/.test(html));
ok('atlas state.karyotypeLineage.per_focal_chr (canonical)',
   /per_focal_chr:\s*normalizedEntries/.test(html));
ok('atlas state.karyotypeLineage.focal_species (canonical)',
   /focal_species:\s*focalSpecies/.test(html));
ok('new lookup _msGetKaryotypeEntryForFocalChr defined',
   /function _msGetKaryotypeEntryForFocalChr\(/.test(html));
ok('legacy lookup _msGetKaryotypeEntryForCgarChr is wrapper',
   /function _msGetKaryotypeEntryForCgarChr\(focalChr\)\s*\{[\s\S]{0,80}return\s+_msGetKaryotypeEntryForFocalChr/.test(html));

ok('window exposes _msGetKaryotypeEntryForFocalChr',
   /window\._msGetKaryotypeEntryForFocalChr\s*=/.test(html));

console.log('\n=== Source-level: rules + verdicts ===');
ok('focal_lineage_fission verdict emitted',
   /verdict:\s*'focal_lineage_fission'/.test(html));
ok('sister_lineage_fission verdict emitted',
   /verdict:\s*'sister_lineage_fission'/.test(html));
ok('hardcoded cgar_class signal removed (replaced with focal_class)',
   !/signals\.cgar_class\s*===\s*'1-2'/.test(html));
ok('hardcoded cmac_class signal removed (replaced with sister_one_to_two)',
   !/signals\.cmac_class\s*===\s*'1-2'/.test(html));
ok('signals.focal_species exposed',
   /focal_species:\s*focalSp/.test(html));
ok('signals.sister_classes dict',
   /sister_classes:\s*\{\}/.test(html));

console.log('\n=== Source-level: age-model RULE 0 backward compat ===');
ok('age-model recognizes focal_lineage_fission',
   /lineageVerdicts\s*=\s*\[[^\]]*'focal_lineage_fission'/.test(html));
ok('age-model recognizes sister_lineage_fission',
   /lineageVerdicts\s*=\s*\[[^\]]*'sister_lineage_fission'/.test(html));
ok('age-model still recognizes cgar_lineage_fission (backward compat)',
   /lineageVerdicts\s*=\s*\[[^\]]*'cgar_lineage_fission'/.test(html));
ok('age-model still recognizes cmac_lineage_fission (backward compat)',
   /lineageVerdicts\s*=\s*\[[^\]]*'cmac_lineage_fission'/.test(html));

console.log('\n=== Source-level: CSS for new verdict classes ===');
ok('CSS .ms-karyo-verdict-focal defined',
   html.indexOf('.ms-karyo-verdict-focal') >= 0);
ok('CSS .ms-karyo-verdict-sister defined',
   html.indexOf('.ms-karyo-verdict-sister') >= 0);
ok('CSS .ms-karyo-verdict-cgar still defined (backward compat)',
   html.indexOf('.ms-karyo-verdict-cgar') >= 0);
ok('CSS .ms-karyo-verdict-cmac still defined (backward compat)',
   html.indexOf('.ms-karyo-verdict-cmac') >= 0);

console.log('\n=== Source-level: HTML renderer uses focal_species for substitution ===');
ok('HTML renderer substitutes focal species name in verdict label',
   /focalSp[\s\S]{0,80}'focal_lineage_fission'[\s\S]{0,200}lineage fission/.test(html));
ok('HTML renderer substitutes sister species name in verdict label',
   /sisterCarriers[\s\S]{0,200}lineage fission/.test(html));

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
const fnGetFocal = pullFunction(html, '_msGetKaryotypeEntryForFocalChr');
const fnGetCgar  = pullFunction(html, '_msGetKaryotypeEntryForCgarChr');
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

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  const state = { karyotypeLineage: null };
  const _MS_DEFAULT_SPECIES = ${defaultSp || '[]'};
  ${fnIs}
  ${fnStore}
  ${fnGetSp}
  ${fnGetFocal}
  ${fnGetCgar}
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
    _msGetKaryotypeEntryForFocalChr,
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

// === Schema acceptance — new species-agnostic shape
console.log('\n=== JSON schema: species-agnostic shape ===');
const dataNew = {
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG10',
    classes_by_species: {
      Cgar:  { class: '1-2', targets: ['C_gar_LG10', 'C_gar_LG14'] },
      Cmac:  { class: '1-1', targets: ['C_mac_LG09'] },
      Tros:  { class: '1-1', targets: ['Tros_chr14'] },
    },
  }],
};
ok('new shape (per_focal_chr) accepted', m._isKaryotypeLineageJSON(dataNew) === true);
ok('new shape stores',                    m._storeKaryotypeLineage(dataNew) === true);
ok('focal_species preserved',             m.state.karyotypeLineage.focal_species === 'Cgar');
ok('sister_species preserved',
   m.state.karyotypeLineage.sister_species[0] === 'Cmac');
ok('per_focal_chr available canonically',
   m.state.karyotypeLineage.per_focal_chr[0].focal_chr === 'C_gar_LG10');

// === Schema legacy shape with per_cgar_chr still works
console.log('\n=== JSON schema: legacy per_cgar_chr accepted ===');
const dataLegacy = {
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  outgroup_species: ['Tros'],
  per_cgar_chr: [{ cgar_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG01', 'C_mac_LG18'] },
      Tros:  { class: '1-1', targets: ['Tros|22'] },
    },
  }],
};
ok('legacy shape accepted', m._isKaryotypeLineageJSON(dataLegacy) === true);
ok('legacy shape stores',   m._storeKaryotypeLineage(dataLegacy) === true);
ok('legacy normalized to per_focal_chr',
   Array.isArray(m.state.karyotypeLineage.per_focal_chr) &&
   m.state.karyotypeLineage.per_focal_chr[0].focal_chr === 'C_gar_LG28');
ok('legacy: focal_species inferred from data (Cgar maps 1-1 to itself)',
   m.state.karyotypeLineage.focal_species === 'Cgar');
ok('legacy: sister_species defaults to []',
   Array.isArray(m.state.karyotypeLineage.sister_species) &&
   m.state.karyotypeLineage.sister_species.length === 0);

// === Polarization — Cgar-focal scenario
console.log('\n=== Polarization with Cgar as focal ===');
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG10',
    classes_by_species: {
      Cgar:  { class: '1-2', targets: ['C_gar_LG10', 'C_gar_LG14'] },
      Cmac:  { class: '1-1', targets: ['C_mac_LG09'] },
      Cfus:  { class: '1-1', targets: ['C_fus_LG09'] },
      Tros:  { class: '1-1', targets: ['Tros_chr14'] },
    },
  }],
});
const v_focal_cgar = m._msPolarizeKaryotypeEvent('C_gar_LG10');
ok('Cgar focal + Cgar 1-2 → focal_lineage_fission',
   v_focal_cgar.verdict === 'focal_lineage_fission');
ok('signals.focal_species === "Cgar"',
   v_focal_cgar.signals.focal_species === 'Cgar');
ok('signals.focal_class === "1-2"',
   v_focal_cgar.signals.focal_class === '1-2');

// === Polarization — Cmac-focal scenario (different paper!)
console.log('\n=== Polarization with Cmac as focal (different paper) ===');
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cmac',
  sister_species: ['Cgar'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_mac_LG14',
    classes_by_species: {
      Cmac: { class: '1-2', targets: ['C_mac_LG14a', 'C_mac_LG14b'] },
      Cgar: { class: '1-1', targets: ['C_gar_LG12'] },
      Cfus: { class: '1-1', targets: ['C_fus_LG12'] },
      Tros: { class: '1-1', targets: ['Tros_chr19'] },
    },
  }],
});
const v_focal_cmac = m._msPolarizeKaryotypeEvent('C_mac_LG14');
ok('Cmac focal + Cmac 1-2 → focal_lineage_fission (NOT cmac_lineage_fission)',
   v_focal_cmac.verdict === 'focal_lineage_fission');
ok('signals.focal_species === "Cmac"',
   v_focal_cmac.signals.focal_species === 'Cmac');

// === Polarization — sister-species detection
console.log('\n=== Sister-species detection ===');
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'] },
      Cmac:  { class: '1-2', targets: ['C_mac_LG01', 'C_mac_LG18'] },
      Cfus:  { class: '1-1', targets: ['C_fus_LG07'] },
      Tros:  { class: '1-1', targets: ['Tros|22'] },
    },
  }],
});
const v_sister_cmac = m._msPolarizeKaryotypeEvent('C_gar_LG28');
ok('Cmac sister 1-2 + Cgar focal 1-1 → sister_lineage_fission',
   v_sister_cmac.verdict === 'sister_lineage_fission');
ok('signals.sister_one_to_two === 1',
   v_sister_cmac.signals.sister_one_to_two === 1);

// === Custom non-catfish species names — atlas should not care
console.log('\n=== Generic species names (atlas is species-agnostic) ===');
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Salmosalar',
  sister_species: ['Othrutta'],
  outgroup_species: ['Esoxniger'],
  per_focal_chr: [{ focal_chr: 'Ssal_chr05',
    classes_by_species: {
      Salmosalar: { class: '1-2', targets: ['Ssal_chr05', 'Ssal_chr19'] },
      Othrutta:   { class: '1-1', targets: ['Otru_chr08'] },
      Esoxniger:  { class: '1-1', targets: ['Enig_chr11'] },
    },
  }],
});
const v_salmon = m._msPolarizeKaryotypeEvent('Ssal_chr05');
ok('non-catfish focal species → focal_lineage_fission',
   v_salmon.verdict === 'focal_lineage_fission');
ok('non-catfish: signals.focal_species correctly read',
   v_salmon.signals.focal_species === 'Salmosalar');

// === HTML render — focal species name substituted in label
console.log('\n=== Verdict label substitution ===');
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG10',
    classes_by_species: {
      Cgar: { class: '1-2', targets: ['T1', 'T2'] },
      Cmac: { class: '1-1', targets: ['M1'] },
      Cfus: { class: '1-1', targets: ['F1'] },
      Tros: { class: '1-1', targets: ['R1'] },
    },
  }],
});
const html_focal = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG10' });
ok('verdict label uses "Cgar lineage fission" (substituted)',
   /Cgar lineage fission/.test(html_focal));
ok('verdict CSS class is .ms-karyo-verdict-focal',
   /class="ms-karyo-verdict-chip ms-karyo-verdict-focal"/.test(html_focal));

m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 1, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac', 'Cfus'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar: { class: '1-1', targets: ['T0'] },
      Cmac: { class: '1-2', targets: ['T1', 'T2'] },
      Cfus: { class: '1-1', targets: ['T3'] },
      Tros: { class: '1-1', targets: ['T4'] },
    },
  }],
});
const html_sister = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('verdict label uses "Cmac lineage fission" (sister carrier)',
   /Cmac lineage fission/.test(html_sister));
ok('verdict CSS class is .ms-karyo-verdict-sister',
   /class="ms-karyo-verdict-chip ms-karyo-verdict-sister"/.test(html_sister));

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
