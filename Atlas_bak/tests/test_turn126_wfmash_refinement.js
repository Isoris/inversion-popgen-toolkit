// =============================================================================
// turn 126 integration test — wfmash refinement of mashmap karyotype calls
//
// Covers:
//   1. Source-level: refinement helpers exist, exposed on window, render
//      uses refinement-aware polarization.
//   2. _msSummarizeRefinement counts states correctly across cells.
//   3. _msAdjustConfidenceForRefinement boosts/drops confidence per rules.
//   4. _msGetEffectiveClassForCell + Targets respect refuted/refined overrides.
//   5. _msBuildRefinementChipHtml produces appropriate chip for each state.
//   6. End-to-end: HTML render shows refinement column, summary line,
//      strikethrough on refuted classes, and adjusted verdict confidence.
//   7. Empty-state: no refinement fields → no refinement column, no summary line,
//      verdict unchanged.
//   8. Regression: turns 121-125 still pass.
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
console.log('\n=== Source-level: refinement helpers ===');

ok('_msSummarizeRefinement defined',
   /function _msSummarizeRefinement\(/.test(html));
ok('_msAdjustConfidenceForRefinement defined',
   /function _msAdjustConfidenceForRefinement\(/.test(html));
ok('_msPolarizeKaryotypeEventWithRefinement defined',
   /function _msPolarizeKaryotypeEventWithRefinement\(/.test(html));
ok('_msGetEffectiveClassForCell defined',
   /function _msGetEffectiveClassForCell\(/.test(html));
ok('_msGetEffectiveTargetsForCell defined',
   /function _msGetEffectiveTargetsForCell\(/.test(html));
ok('_msBuildRefinementChipHtml defined',
   /function _msBuildRefinementChipHtml\(/.test(html));

ok('window exports _msSummarizeRefinement',
   /window\._msSummarizeRefinement\s*=/.test(html));
ok('window exports _msAdjustConfidenceForRefinement',
   /window\._msAdjustConfidenceForRefinement\s*=/.test(html));
ok('window exports _msPolarizeKaryotypeEventWithRefinement',
   /window\._msPolarizeKaryotypeEventWithRefinement\s*=/.test(html));
ok('window exports _msGetEffectiveClassForCell',
   /window\._msGetEffectiveClassForCell\s*=/.test(html));
ok('window exports _msBuildRefinementChipHtml',
   /window\._msBuildRefinementChipHtml\s*=/.test(html));

console.log('\n=== Source-level: render uses refinement-aware polarization ===');
ok('render calls _msPolarizeKaryotypeEventWithRefinement',
   /_msPolarizeKaryotypeEventWithRefinement\(focalChr\)/.test(html));
ok('render builds refinement summary line',
   /ms-karyo-refine-summary/.test(html));
ok('render uses _msBuildRefinementChipHtml in row',
   /_msBuildRefinementChipHtml\(c\)/.test(html));
ok('render adds strikethrough class on refuted',
   /ms-karyo-class-refuted/.test(html));
ok('render conditional refinement column',
   /hasAnyRefinement/.test(html));

console.log('\n=== Source-level: CSS classes ===');
const cssClasses = [
  'ms-refine-chip', 'ms-refine-chip-confirmed', 'ms-refine-chip-refuted',
  'ms-refine-chip-refined', 'ms-refine-chip-failed', 'ms-refine-chip-pending',
  'ms-karyo-class-refuted', 'ms-karyo-refine-col',
  'ms-karyo-refine-summary', 'ms-karyo-refine-summary-label',
  'ms-refine-note', 'ms-refine-note-warn',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined', html.indexOf('.' + cls) >= 0);
}

console.log('\n=== Regression: prior turns still wired ===');
ok('turn 121 _isSyntenyMultispeciesJSON still defined',
   /function _isSyntenyMultispeciesJSON\(/.test(html));
ok('turn 122 _msAutoSuggestAgeModel still defined',
   /function _msAutoSuggestAgeModel\(/.test(html));
ok('turn 123 _isCompTEFragilityJSON still defined',
   /function _isCompTEFragilityJSON\(/.test(html));
ok('turn 124 _isKaryotypeLineageJSON still defined',
   /function _isKaryotypeLineageJSON\(/.test(html));
ok('turn 125 species-agnostic verdicts still emitted',
   /verdict:\s*'focal_lineage_fission'/.test(html));

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

const fns = {};
for (const fnName of [
  '_isKaryotypeLineageJSON', '_storeKaryotypeLineage',
  '_msGetKaryotypeEntryForFocalChr', '_msGetKaryotypeEntryForCgarChr',
  '_msGetEffectiveSpeciesList',
  '_msPolarizeKaryotypeEvent',
  '_msSummarizeRefinement', '_msAdjustConfidenceForRefinement',
  '_msPolarizeKaryotypeEventWithRefinement',
  '_msGetEffectiveClassForCell', '_msGetEffectiveTargetsForCell',
  '_msBuildRefinementChipHtml',
  '_msBuildKaryotypeContextHtml',
]) {
  fns[fnName] = pullFunction(html, fnName);
  if (!fns[fnName]) {
    console.log('  FAIL pulling ' + fnName);
    process.exit(1);
  }
}
const defaultSp = pullConst(html, '_MS_DEFAULT_SPECIES');

const sandboxSrc = `
  function _esc(s){ return String(s == null ? '' : s)
    .replace(/&/g,'&amp;').replace(/</g,'&lt;')
    .replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  const state = { karyotypeLineage: null };
  const _MS_DEFAULT_SPECIES = ${defaultSp || '[]'};
  ${Object.values(fns).join('\n')}
  module.exports = { state, ${Object.keys(fns).map(k => k.replace(/^_/, '_')).join(', ')} };
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

// === _msGetEffectiveClassForCell
console.log('\n=== Effective class/targets respect refinement ===');
ok('null cell → null',
   m._msGetEffectiveClassForCell(null) === null);
ok('no refinement → mashmap class',
   m._msGetEffectiveClassForCell({ class: '1-2', targets: ['T1','T2'] }) === '1-2');
ok('confirmed → mashmap class',
   m._msGetEffectiveClassForCell({ class: '1-2', refined_by_wfmash: 'confirmed' }) === '1-2');
ok('refuted with wfmash_class → wfmash class',
   m._msGetEffectiveClassForCell({ class: '1-2', refined_by_wfmash: 'refuted', wfmash_class: '1-1' }) === '1-1');
ok('refined with wfmash_class → wfmash class',
   m._msGetEffectiveClassForCell({ class: '1-2', refined_by_wfmash: 'refined', wfmash_class: '1-3' }) === '1-3');
ok('refuted but no wfmash_class → still mashmap class (defensive)',
   m._msGetEffectiveClassForCell({ class: '1-2', refined_by_wfmash: 'refuted' }) === '1-2');

ok('targets fallback when not refuted',
   JSON.stringify(m._msGetEffectiveTargetsForCell({ class: '1-2', targets: ['A','B'] }))
   === JSON.stringify(['A','B']));
ok('refuted → wfmash_targets override',
   JSON.stringify(m._msGetEffectiveTargetsForCell({
     class: '1-2', targets: ['A','B'],
     refined_by_wfmash: 'refuted', wfmash_class: '1-1',
     wfmash_targets: [{ chrom: 'X', start_bp: 0, end_bp: 1000 }],
   })) === JSON.stringify(['X']));

// === _msSummarizeRefinement
console.log('\n=== Refinement summary counts ===');
const entry1 = {
  focal_chr: 'C_gar_LG28',
  classes_by_species: {
    Cgar:  { class: '1-1', refined_by_wfmash: 'confirmed' },
    Cmac:  { class: '1-2', refined_by_wfmash: 'confirmed' },
    Cfus:  { class: '1-1', refined_by_wfmash: 'confirmed' },
    Capus: { class: '1-2', refined_by_wfmash: 'refuted', wfmash_class: '1-1' },
    Phyp:  { class: '1-1' /* not_attempted */ },
    Tros:  { class: '1-1', refined_by_wfmash: 'failed' },
  },
};
const s1 = m._msSummarizeRefinement(entry1);
ok('confirmed count',     s1.confirmed === 3);
ok('refuted count',       s1.refuted === 1);
ok('refined count',       s1.refined === 0);
ok('failed count',        s1.failed === 1);
ok('not_attempted count', s1.not_attempted === 1);
ok('total count',         s1.total === 6);

const empty_summary = m._msSummarizeRefinement(null);
ok('null entry → all zeros', empty_summary.total === 0);
ok('empty classes_by_species → all zeros',
   m._msSummarizeRefinement({ focal_chr: 'X', classes_by_species: {} }).total === 0);

// === _msAdjustConfidenceForRefinement
console.log('\n=== Confidence adjustment ===');
const baseVerdict = {
  verdict: 'focal_lineage_fission',
  confidence: 'medium',
  rationale: 'baseline rationale',
  signals: {},
};

// All confirmed → boost medium → high
const adjBoost = m._msAdjustConfidenceForRefinement(baseVerdict, {
  confirmed: 5, refuted: 0, refined: 0, failed: 0, not_attempted: 0, total: 5,
});
ok('all confirmed → confidence boosted to high', adjBoost.confidence === 'high');
ok('all confirmed → rationale gets boost note', /Boosted by wfmash/.test(adjBoost.rationale));

// 80% confirmed, 0 refuted → still boost
const adj80 = m._msAdjustConfidenceForRefinement(baseVerdict, {
  confirmed: 4, refuted: 0, refined: 0, failed: 1, not_attempted: 0, total: 5,
});
ok('80% confirmed → still boosted', adj80.confidence === 'high');

// 70% confirmed, 0 refuted → not enough to boost
const adj70 = m._msAdjustConfidenceForRefinement(baseVerdict, {
  confirmed: 3, refuted: 0, refined: 0, failed: 1, not_attempted: 1, total: 5,
});
ok('70% confirmed → no boost (still medium)', adj70.confidence === 'medium');

// ≥30% refuted → drop tier
const adjDrop = m._msAdjustConfidenceForRefinement(baseVerdict, {
  confirmed: 1, refuted: 2, refined: 0, failed: 0, not_attempted: 2, total: 5,
});
ok('40% refuted → confidence dropped to low', adjDrop.confidence === 'low');
ok('drop case → rationale gets caution note', /Caution.*refuted/.test(adjDrop.rationale));

// <50% touched → no adjustment
const adjSparse = m._msAdjustConfidenceForRefinement(baseVerdict, {
  confirmed: 1, refuted: 0, refined: 0, failed: 0, not_attempted: 4, total: 5,
});
ok('only 20% touched → no adjustment', adjSparse.confidence === 'medium');

// Unresolved verdict — never adjust
const unresolvedVerdict = { verdict: 'unresolved', confidence: 'low', rationale: '', signals: {} };
const adjUnres = m._msAdjustConfidenceForRefinement(unresolvedVerdict, {
  confirmed: 5, refuted: 0, refined: 0, failed: 0, not_attempted: 0, total: 5,
});
ok('unresolved verdict not adjusted', adjUnres.confidence === 'low');

// Null verdict
const adjNull = m._msAdjustConfidenceForRefinement(null, { total: 5 });
ok('null verdict → returned as-is', adjNull === null);

// === _msBuildRefinementChipHtml
console.log('\n=== Refinement chip HTML ===');
ok('null cell → empty string',
   m._msBuildRefinementChipHtml(null) === '');
ok('no refinement → pending chip',
   /ms-refine-chip-pending/.test(m._msBuildRefinementChipHtml({ class: '1-2' })));
ok('confirmed → confirmed chip',
   /ms-refine-chip-confirmed/.test(m._msBuildRefinementChipHtml({ class: '1-2', refined_by_wfmash: 'confirmed' })));
ok('confirmed chip has checkmark',
   /\u2713/.test(m._msBuildRefinementChipHtml({ class: '1-2', refined_by_wfmash: 'confirmed' })));
ok('refuted chip with wfmash_class shows arrow',
   /\u2192 1-1/.test(m._msBuildRefinementChipHtml({ class: '1-2', refined_by_wfmash: 'refuted', wfmash_class: '1-1' })));
ok('refined chip shows new class',
   /1-3/.test(m._msBuildRefinementChipHtml({ class: '1-2', refined_by_wfmash: 'refined', wfmash_class: '1-3' })));
ok('failed chip → failed class',
   /ms-refine-chip-failed/.test(m._msBuildRefinementChipHtml({ class: '1-2', refined_by_wfmash: 'failed' })));

// === End-to-end render
console.log('\n=== End-to-end render with refinement ===');

m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 2, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'], refined_by_wfmash: 'confirmed' },
      Cmac:  { class: '1-2', targets: ['T1','T2'],   refined_by_wfmash: 'confirmed' },
      Cfus:  { class: '1-1', targets: ['F1'],         refined_by_wfmash: 'confirmed' },
      Capus: { class: '1-1', targets: ['A1'],         refined_by_wfmash: 'confirmed' },
      Tros:  { class: '1-1', targets: ['R1'],         refined_by_wfmash: 'confirmed' },
    },
  }],
});
const htmlConfirmed = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('all-confirmed render shows refinement summary line',
   /ms-karyo-refine-summary/.test(htmlConfirmed));
ok('all-confirmed render contains "5</b> confirmed"',
   /<b>5<\/b> confirmed/.test(htmlConfirmed));
ok('all-confirmed render contains refinement column header',
   /<th>wfmash refinement<\/th>/.test(htmlConfirmed));
ok('all-confirmed render contains confirmed chips in rows',
   (htmlConfirmed.match(/ms-refine-chip-confirmed/g) || []).length >= 5);
ok('all-confirmed → confidence boosted to high',
   /\(high confidence\)/.test(htmlConfirmed));

// Now refute the sister cell
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 2, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar:  { class: '1-1', targets: ['C_gar_LG28'], refined_by_wfmash: 'confirmed' },
      Cmac:  { class: '1-2', targets: ['T1','T2'],   refined_by_wfmash: 'refuted', wfmash_class: '1-1' },
      Cfus:  { class: '1-1', targets: ['F1'],         refined_by_wfmash: 'confirmed' },
      Tros:  { class: '1-1', targets: ['R1'],         refined_by_wfmash: 'confirmed' },
    },
  }],
});
const htmlRefuted = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('refuted render contains refuted chip',
   /ms-refine-chip-refuted/.test(htmlRefuted));
ok('refuted render contains strikethrough class',
   /ms-karyo-class-refuted/.test(htmlRefuted));
ok('refuted render contains refuted count in summary',
   /<b>1<\/b> refuted/.test(htmlRefuted));

// No refinement at all → no extra column / summary
m._storeKaryotypeLineage({
  tool: 'karyotype_lineage_v1', schema_version: 2, params: {},
  focal_species: 'Cgar',
  sister_species: ['Cmac'],
  outgroup_species: ['Tros'],
  per_focal_chr: [{ focal_chr: 'C_gar_LG28',
    classes_by_species: {
      Cgar: { class: '1-1', targets: ['C_gar_LG28'] },
      Cmac: { class: '1-2', targets: ['T1','T2'] },
      Tros: { class: '1-1', targets: ['R1'] },
    },
  }],
});
const htmlNoRef = m._msBuildKaryotypeContextHtml({ id: 'b', gar_chr: 'C_gar_LG28' });
ok('no-refinement render hides summary line',
   !/ms-karyo-refine-summary/.test(htmlNoRef));
ok('no-refinement render hides refinement column header',
   !/<th>wfmash refinement<\/th>/.test(htmlNoRef));
ok('no-refinement render verdict still works',
   /sister_lineage_fission|sister/i.test(htmlNoRef));

// =============================================================================
console.log('\n=============================================================');
console.log('PASSED: ' + pass + ' / ' + (pass + fail));
if (fail > 0) {
  console.log('FAILED: ' + fail);
  process.exit(1);
} else {
  console.log('ALL CHECKS PASSED');
}
