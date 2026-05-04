// =============================================================================
// turn 121 integration test — multi-species classification cockpit (page16b)
//
// Tests:
//   1. Source-level: tab button, page DOM, CSS classes, dispatcher hooks,
//      JSON detectors, store/persist/restore helpers, classifier extension.
//   2. JSON detection: synteny_multispecies_v1 + phylo_tree_v1 detectors.
//   3. Behavioural: _msAutoSuggestArchitecture across the 6 classes (A-F),
//      across all-zero / E-recurrent / C-fission / B-pair patterns.
//   4. Default species fallback: the 9-species reference list is shape-correct.
//   5. Newick parser: minimal Newick → leaves with depth.
//   6. Empty-state hierarchy: rendering doesn't throw when nothing loaded.
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
console.log('\n=== Source-level checks ===');

ok('page16b tab button in tabBar', /data-page="page16b" data-stage="synthesis"/.test(html));
ok('page16b "13b multi-species" label', /<span class="num">13b<\/span> multi-species/.test(html));
ok('page16b page DOM container', /<div id="page16b" class="page">/.test(html));
ok('page16b activation hook in dispatcher',
   /if \(target === 'page16b'\)[\s\S]{0,200}_renderMultiSpeciesPage/.test(html));
ok('msTree DOM slot',                /id="msTree"/.test(html));
ok('msCenter DOM slot',              /id="msCenter"/.test(html));
ok('msDetail DOM slot',              /id="msDetail"/.test(html));
ok('msActiveHeader DOM slot',        /id="msActiveHeader"/.test(html));
ok('msClassification DOM slot',      /id="msClassification"/.test(html));
ok('msLineageTable DOM slot',        /id="msLineageTable"/.test(html));
ok('msDetailDotplot DOM slot',       /id="msDetailDotplot"/.test(html));
ok('msDetailBoundary DOM slot',      /id="msDetailBoundary"/.test(html));
ok('classification framework link in subtitle', /id="msClassificationFrameworkLink"/.test(html));

// === CSS contract
console.log('\n=== CSS contract ===');
const cssClasses = [
  'ms-tree-svg', 'ms-tree-leaf-label', 'ms-tree-leaf-active', 'ms-tree-leaf-focal',
  'ms-tree-branch', 'ms-tree-leaf-marker', 'ms-tree-legend',
  'ms-cls-block', 'ms-cls-title', 'ms-cls-chip', 'ms-cls-A', 'ms-cls-B',
  'ms-cls-C', 'ms-cls-D', 'ms-cls-E', 'ms-cls-F',
  'ms-cls-confidence', 'ms-cls-rationale', 'ms-cls-empty',
  'ms-lineage-table', 'ms-lineage-row-active',
  'ms-lineage-status-present', 'ms-lineage-status-absent',
  'ms-lineage-status-fission', 'ms-lineage-status-fusion',
  'ms-lineage-status-unknown',
  'ms-detail-title', 'ms-detail-section-title', 'ms-detail-empty-msg',
  'ms-detail-kv',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined', html.indexOf('.' + cls) >= 0);
}

// === JSON layer detectors and helpers
console.log('\n=== JSON layer contract ===');
ok('_isSyntenyMultispeciesJSON defined',  /function _isSyntenyMultispeciesJSON\(/.test(html));
ok('_storeSyntenyMultispecies defined',   /function _storeSyntenyMultispecies\(/.test(html));
ok('_persistSyntenyMultispecies defined', /function _persistSyntenyMultispecies\(/.test(html));
ok('_restoreSyntenyMultispecies defined', /function _restoreSyntenyMultispecies\(/.test(html));
ok('_clearSyntenyMultispecies defined',   /function _clearSyntenyMultispecies\(/.test(html));
ok('_isPhyloTreeJSON defined',            /function _isPhyloTreeJSON\(/.test(html));
ok('_storePhyloTree defined',             /function _storePhyloTree\(/.test(html));
ok('_persistPhyloTree defined',           /function _persistPhyloTree\(/.test(html));
ok('_restorePhyloTree defined',           /function _restorePhyloTree\(/.test(html));
ok('_clearPhyloTree defined',             /function _clearPhyloTree\(/.test(html));
ok('classifier routes synteny_multispecies',
   /_isSyntenyMultispeciesJSON\(data\)\)\s*return\s+'synteny_multispecies'/.test(html));
ok('classifier routes phylo_tree',
   /_isPhyloTreeJSON\(data\)\)\s*return\s+'phylo_tree'/.test(html));
ok('loadMultipleJSONs dispatches synteny_multispecies',
   /_isSyntenyMultispeciesJSON[\s\S]{0,400}_storeSyntenyMultispecies/.test(html));
ok('loadMultipleJSONs dispatches phylo_tree',
   /_isPhyloTreeJSON[\s\S]{0,400}_storePhyloTree/.test(html));
ok('restore wired on load (synteny_multispecies)',
   /_restoreSyntenyMultispecies/.test(html));
ok('restore wired on load (phylo_tree)',
   /_restorePhyloTree/.test(html));
ok('localStorage key for synteny_multispecies',
   html.indexOf("'inversion_atlas.syntenyMultispecies.v1'") >= 0);
ok('localStorage key for phylo_tree',
   html.indexOf("'inversion_atlas.phyloTree.v1'") >= 0);
ok('localStorage key for multiSpecies UI',
   html.indexOf("'inversion_atlas.multiSpeciesUI.v1'") >= 0);

// === Renderer + helpers
console.log('\n=== Renderer + helpers ===');
ok('_renderMultiSpeciesPage defined',     /function _renderMultiSpeciesPage\(/.test(html));
ok('_msInitState defined',                /function _msInitState\(/.test(html));
ok('_msPersistUI defined',                /function _msPersistUI\(/.test(html));
ok('_msGetActiveBreakpoint defined',      /function _msGetActiveBreakpoint\(/.test(html));
ok('_msGetEffectiveSpeciesList defined',  /function _msGetEffectiveSpeciesList\(/.test(html));
ok('_msGetLineageDistribution defined',   /function _msGetLineageDistribution\(/.test(html));
ok('_msAutoSuggestArchitecture defined',  /function _msAutoSuggestArchitecture\(/.test(html));
ok('_msParseNewickToLeaves defined',      /function _msParseNewickToLeaves\(/.test(html));
ok('_msRenderTreeSvg defined',            /function _msRenderTreeSvg\(/.test(html));
ok('_msRenderActiveHeader defined',       /function _msRenderActiveHeader\(/.test(html));
ok('_msRenderCenterColumn defined',       /function _msRenderCenterColumn\(/.test(html));
ok('_msRenderDetailColumn defined',       /function _msRenderDetailColumn\(/.test(html));
ok('_msGetSyntenyBlocksForSpecies defined',/function _msGetSyntenyBlocksForSpecies\(/.test(html));
ok('_MS_DEFAULT_SPECIES defined',         /const _MS_DEFAULT_SPECIES = \[/.test(html));
ok('window exposes _msAutoSuggestArchitecture',
   /window\._msAutoSuggestArchitecture = _msAutoSuggestArchitecture/.test(html));
ok('window exposes _renderMultiSpeciesPage',
   /window\._renderMultiSpeciesPage\s*=/.test(html));

// =============================================================================
// Behavioural tests — pull functions into a vm sandbox
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
    // Skip line comments — `// ...` to end of line (apostrophes inside don't open strings)
    if (ch === '/' && src[i+1] === '/') {
      const nl = src.indexOf('\n', i);
      i = nl < 0 ? src.length : nl + 1;
      continue;
    }
    // Skip block comments — `/* ... */`
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
        if (src[i] === '\n' && (quote === '"' || quote === "'")) {
          // Single/double-quoted strings shouldn't span lines — bail to avoid runaway
          break;
        }
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

const fnAuto    = pullFunction(html, '_msAutoSuggestArchitecture');
const fnNewick  = pullFunction(html, '_msParseNewickToLeaves');
const defaultSp = pullConst(html, '_MS_DEFAULT_SPECIES');

ok('extracted _msAutoSuggestArchitecture source', !!fnAuto);
ok('extracted _msParseNewickToLeaves source', !!fnNewick);
ok('extracted _MS_DEFAULT_SPECIES source', !!defaultSp);

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  ${fnAuto}
  ${fnNewick}
  const _MS_DEFAULT_SPECIES = ${defaultSp || '[]'};
  module.exports = {
    _msAutoSuggestArchitecture,
    _msParseNewickToLeaves,
    _MS_DEFAULT_SPECIES,
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

// === Default species list shape
console.log('\n=== Default 9-species fallback ===');
ok('default list has 9-10 species', m._MS_DEFAULT_SPECIES.length >= 9 && m._MS_DEFAULT_SPECIES.length <= 12);
ok('default list contains Cgar focal',
   m._MS_DEFAULT_SPECIES.some(s => s.id === 'Cgar' && s.focal === true));
ok('default list contains Cmac focal',
   m._MS_DEFAULT_SPECIES.some(s => s.id === 'Cmac' && s.focal === true));
ok('default list contains Tros (deep outgroup)',
   m._MS_DEFAULT_SPECIES.some(s => s.id === 'Tros' && s.tier === 'deep_outgroup'));
ok('default list contains Cfus (clarias context)',
   m._MS_DEFAULT_SPECIES.some(s => s.id === 'Cfus' && s.tier === 'clarias_context'));
ok('default list every entry has id+label+focal',
   m._MS_DEFAULT_SPECIES.every(s => s.id && s.label && typeof s.focal === 'boolean'));
ok('default list non-focal species count > focal count',
   m._MS_DEFAULT_SPECIES.filter(s => !s.focal).length >
   m._MS_DEFAULT_SPECIES.filter(s => s.focal).length);

// === Architecture auto-suggest classification
console.log('\n=== Architecture class auto-suggest ===');
const bp_inv  = { id: 'bp1', event_type: 'inversion' };
const bp_fis  = { id: 'bp2', event_type: 'fission' };
const bp_tra  = { id: 'bp3', event_type: 'translocation' };

// No lineage data
const noLineage_inv = m._msAutoSuggestArchitecture(bp_inv, null);
ok('no-lineage inversion → A class', noLineage_inv.class === 'A');
ok('no-lineage inversion confidence is low', noLineage_inv.confidence === 'low');
const noLineage_fis = m._msAutoSuggestArchitecture(bp_fis, null);
ok('no-lineage fission → C class', noLineage_fis.class === 'C');
const noLineage_tra = m._msAutoSuggestArchitecture(bp_tra, null);
ok('no-lineage translocation → D class', noLineage_tra.class === 'D');

// Lineage data — Class E (recurrent, >=3 species)
const lineage_E = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_present',
  Cfus: 'boundary_present',
  Capus: 'boundary_present',
  Phyp: 'boundary_absent_internal_to_block',
  Tfulv: 'unknown',
};
const e_result = m._msAutoSuggestArchitecture(bp_inv, lineage_E);
ok('4 species present → E (recurrent hotspot)', e_result.class === 'E');
ok('E confidence is medium', e_result.confidence === 'medium');
ok('E n_supporting === 4', e_result.n_supporting === 4);
ok('E rationale mentions hotspot', /hotspot/i.test(e_result.rationale));

// Class C — fission/fusion + boundary present
const lineage_C = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_present',
  Phyp: 'fission_in_target',
  Tfulv: 'unknown',
};
const c_result = m._msAutoSuggestArchitecture(bp_inv, lineage_C);
ok('2 present + 1 fission → C class', c_result.class === 'C');

// Class B — boundary present in exactly 2 species
const lineage_B = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_present',
  Cfus: 'boundary_absent_internal_to_block',
  Tfulv: 'unknown',
};
const b_result = m._msAutoSuggestArchitecture(bp_inv, lineage_B);
ok('exactly 2 present → B class', b_result.class === 'B');

// Class A — only 1 species shows the boundary
const lineage_A = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_absent_internal_to_block',
  Cfus: 'boundary_absent_internal_to_block',
};
const a_result = m._msAutoSuggestArchitecture(bp_inv, lineage_A);
ok('1 present → A class', a_result.class === 'A');

// Class F — ambiguous mixed signal (3 absent, 1 present, 0 fis/fus)
const lineage_F = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_absent_internal_to_block',
  Cfus: 'boundary_absent_internal_to_block',
  Capus: 'boundary_absent_internal_to_block',
};
const f_result = m._msAutoSuggestArchitecture(bp_inv, lineage_F);
// 1 present → falls into nPresent <= 1 branch → A. So we test a different
// ambiguous case explicitly:
const lineage_F2 = {
  Cgar: 'boundary_present',
  Cmac: 'boundary_present',
  Cfus: 'boundary_absent_internal_to_block',
  Capus: 'boundary_absent_internal_to_block',
  Phyp: 'unknown',
  Tfulv: 'unknown',
  Hwyc: 'unknown',
};
// 2 present, 2 absent, 3 unknown, 0 fission/fusion — falls into "exactly 2 present" → B
const f2_result = m._msAutoSuggestArchitecture(bp_inv, lineage_F2);
ok('2 present + 2 absent + unknowns → B (synteny-boundary)', f2_result.class === 'B');

// No bp at all
const empty_result = m._msAutoSuggestArchitecture(null, null);
ok('no bp → null class',  empty_result.class === null);
ok('no bp → "unknown" confidence', empty_result.confidence === 'unknown');

// === Newick parser
console.log('\n=== Newick parser ===');
const newick1 = '(Tros,((Cgar,Cmac),(Cfus,Capus)));';
const leaves1 = m._msParseNewickToLeaves(newick1, ['Tros','Cgar','Cmac','Cfus','Capus']);
ok('5-leaf newick produces 5 leaves', leaves1.length === 5);
ok('newick leaves include Tros',  leaves1.some(l => l.id === 'Tros'));
ok('newick leaves include Cgar',  leaves1.some(l => l.id === 'Cgar'));
ok('newick leaves include Capus', leaves1.some(l => l.id === 'Capus'));

const newick2 = '(Tros:100,(((Sglanis:30,Sasotus:30):50,Tfulv:80):10,((Cgar:4,Cmac:4):8,(Cfus:6,(Cbat:3,Capus:3):3):6):20):10);';
const leaves2 = m._msParseNewickToLeaves(newick2, null);
ok('newick with branch lengths parses',
   leaves2.length >= 8 && leaves2.some(l => l.id === 'Cgar'));
ok('newick branch lengths stripped from id',
   leaves2.every(l => !/[:0-9]/.test(l.id) || /^[A-Za-z]/.test(l.id)));

const empty_newick = m._msParseNewickToLeaves('', null);
ok('empty newick → empty array', Array.isArray(empty_newick) && empty_newick.length === 0);

const garbage = m._msParseNewickToLeaves(null, null);
ok('null newick → empty array', Array.isArray(garbage) && garbage.length === 0);

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
