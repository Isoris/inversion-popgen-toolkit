// =============================================================================
// turn 159 — Propose-all karyotype groups workflow
// =============================================================================
// Adds the "G hotkey → karyotype tab → 🎯 propose all → Enter accepts → auto-jump
// to inheritance" UX layer on top of the existing turn-136 karyotype tab.
//
// What this turn ships:
//   - _gpKaryoProposeAll(candidate) — pure: returns proposal bundle with
//     per-band {name, bandIdx, bandLabel, classification, color,
//     sampleIdxList, n}. Skips empty bands. Embeds H-label classification
//     in name when classifier disagrees with band's nominal label.
//   - _gpKaryoAcceptProposals(bundle) — calls addToManualGroup ×N, returns
//     {accepted, skipped, statuses}. After accept, switches gPanelTab to
//     'inheritance' and re-renders so the user sees the cross-candidate view.
//   - UI: 🎯 propose all groups button in the karyotype-tab footer
//   - UI: proposal strip rendered when state._gpKaryoProposals is set,
//     showing per-band names + classification chips + accept/cancel buttons
//   - Enter accepts when the strip is showing; Esc cancels
//   - Both helpers exported to window
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
// 1. Source-pattern checks — both helpers and UI elements present
// ============================================================================
console.log('\n=== 1. Source-pattern checks ===');

ok('_gpKaryoProposeAll defined',
   /function\s+_gpKaryoProposeAll\s*\(/.test(html));
ok('_gpKaryoAcceptProposals defined',
   /function\s+_gpKaryoAcceptProposals\s*\(/.test(html));
ok('window._gpKaryoProposeAll exported',
   /window\._gpKaryoProposeAll\s*=\s*_gpKaryoProposeAll/.test(html));
ok('window._gpKaryoAcceptProposals exported',
   /window\._gpKaryoAcceptProposals\s*=\s*_gpKaryoAcceptProposals/.test(html));

ok('propose-all button id present in karyo tab',
   html.indexOf('id="gpKaryoProposeAllBtn"') > -1);
ok('proposal strip id present',
   html.indexOf('id="_gpKaryoProposalStrip"') > -1);
ok('accept button id present',
   html.indexOf('id="gpKaryoProposalAcceptBtn"') > -1);
ok('cancel button id present',
   html.indexOf('id="gpKaryoProposalCancelBtn"') > -1);

ok('Enter/Esc handler IIFE wired',
   /_wireGPanelKaryoProposalKeys/.test(html));
ok('Enter/Esc handler binds only when proposals active',
   /if \(!stateNow\._gpKaryoProposals\) return/.test(html));

ok('Existing G hotkey still wired (regression check)',
   /_wireGPanelHotkey/.test(html));
ok('Existing _gpKaryoMakeManualGroup still wired (regression check)',
   /function _gpKaryoMakeManualGroup/.test(html));

// ============================================================================
// 2. Behavioural — _gpKaryoProposeAll on a synthetic candidate
// ============================================================================
console.log('\n=== 2. _gpKaryoProposeAll behaviour ===');

const fnExtract  = pullFunction(html, '_gpKaryoExtractBands');
const fnGetLabel = pullFunction(html, 'getKaryotypeLabel');
const fnLabelEnsure = pullFunction(html, '_karyoLabelEnsureVocab');
const fnColor    = pullFunction(html, '_gpKaryoColor');
const fnPropose  = pullFunction(html, '_gpKaryoProposeAll');
const fnAccept   = pullFunction(html, '_gpKaryoAcceptProposals');

ok('pulled _gpKaryoExtractBands', !!fnExtract);
ok('pulled _gpKaryoColor',        !!fnColor);
ok('pulled _gpKaryoProposeAll',   !!fnPropose);
ok('pulled _gpKaryoAcceptProposals', !!fnAccept);
ok('pulled getKaryotypeLabel',    !!fnGetLabel);

// Build a sandbox with all needed helpers + a fake state.data.samples and
// a fake addToManualGroup that captures calls.
const sbx = {
  console,
  state: {
    gPanelOpen: false,
    gPanelTab: 'karyotype',
    data: {
      samples: [
        { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
        { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
      ],
    },
    k: 3,
    manualGroups: [],
    _gpKaryoProposals: null,
  },
  window: {},
  _addedGroups: [],
  // Constants the helpers reference
  _GP_PALETTE: ['#3b6fb6', '#ffd866', '#d97a2c', '#7ad394', '#a76de2'],
};
sbx.window.state = sbx.state;
sbx.addToManualGroup = function(name, samples, opts) {
  // Mimic the real implementation's exclusivity + return-shape contract
  const cgas = samples.map(si => sbx.state.data.samples[si] && sbx.state.data.samples[si].cga).filter(Boolean);
  let g = sbx.state.manualGroups.find(x => x.name === name);
  if (!g) {
    g = { id: 'mg_' + sbx.state.manualGroups.length, name,
          color: (opts && opts.color) || '#888', scope: 'chrom', members: [] };
    sbx.state.manualGroups.push(g);
  }
  for (const c of cgas) if (!g.members.includes(c)) g.members.push(c);
  sbx.window.addToManualGroup = sbx.addToManualGroup;
  return g;
};
sbx.window.addToManualGroup = sbx.addToManualGroup;
// vm needs addToManualGroup as a global since the function uses bare name
vm.createContext(sbx);

// Inject helpers — use stubs for getKaryotypeLabel/_karyoLabelEnsureVocab
// so we don't have to drag in their constant dependencies. The body
// of _gpKaryoExtractBands and _gpKaryoProposeAll only call these as
// "function ? function() : fallback" so a stub is fine.
vm.runInContext('function _karyoLabelEnsureVocab(){return "legacy";}', sbx);
vm.runInContext('function getKaryotypeLabel(k,K){const lab = ["HOM_REF","HET","HOM_INV"]; return (K===3 && k<3) ? lab[k] : "band " + (k+1);}', sbx);
vm.runInContext('function _gpKaryoColor(k){return _GP_PALETTE[k%_GP_PALETTE.length];}', sbx);
vm.runInContext(fnExtract, sbx);
vm.runInContext(fnPropose, sbx);
vm.runInContext(fnAccept,  sbx);

// Synthetic candidate: 6 samples, K=3, labels = [0,0,1,1,2,-1]
// → band 0: 2 members, band 1: 2 members, band 2: 1 member, 1 ungrouped
const candA = {
  id: 'I3',
  K: 3,
  start_bp: 15000000,
  end_bp:   16000000,
  locked_labels: [0, 0, 1, 1, 2, -1],
};

const bundle = vm.runInContext('_gpKaryoProposeAll(' + JSON.stringify(candA) + ')', sbx);
ok('bundle returned', !!bundle);
ok('bundle.candidateId === "I3"', bundle && bundle.candidateId === 'I3');
ok('bundle.K === 3', bundle && bundle.K === 3);
ok('bundle.n_proposals === 3 (one per non-empty band)',
   bundle && bundle.n_proposals === 3);
ok('bundle.n_empty_bands === 0', bundle && bundle.n_empty_bands === 0);
ok('bundle.n_ungrouped === 1 (the -1 sample)',
   bundle && bundle.n_ungrouped === 1);
ok('bundle.classifier_ran === false (no _classifyHLabelBands in sandbox)',
   bundle && bundle.classifier_ran === false);

const props = bundle ? bundle.proposals : [];
ok('proposal[0] has expected shape',
   props[0] && typeof props[0].name === 'string'
            && Array.isArray(props[0].sampleIdxList)
            && Number.isInteger(props[0].bandIdx));
ok('proposal[0].name starts with "I3_"',
   props[0] && /^I3_/.test(props[0].name));
ok('proposal[0].n === 2', props[0] && props[0].n === 2);
ok('proposal[2].n === 1', props[2] && props[2].n === 1);
ok('proposal sampleIdxList sums to 5 (excludes the -1 ungrouped)',
   props.reduce((s, p) => s + p.sampleIdxList.length, 0) === 5);

// ============================================================================
// 3. Empty bands skipped
// ============================================================================
console.log('\n=== 3. Empty bands skipped ===');

// labels = [0,0,0,0,0,0] — only band 0 has members, bands 1+2 empty
const candAllOne = {
  id: 'I4', K: 3, start_bp: 1, end_bp: 2,
  locked_labels: [0, 0, 0, 0, 0, 0],
};
const bundle2 = vm.runInContext(
  '_gpKaryoProposeAll(' + JSON.stringify(candAllOne) + ')', sbx);
ok('all-one bundle: only 1 proposal',
   bundle2 && bundle2.n_proposals === 1);
ok('all-one bundle: 2 empty bands counted',
   bundle2 && bundle2.n_empty_bands === 2);

// ============================================================================
// 4. _gpKaryoAcceptProposals — calls addToManualGroup ×N
// ============================================================================
console.log('\n=== 4. _gpKaryoAcceptProposals ===');

// Reset manual groups
sbx.state.manualGroups = [];
const acceptResult = vm.runInContext(
  '_gpKaryoAcceptProposals(' + JSON.stringify(bundle) + ')', sbx);

ok('accept result.accepted === 3', acceptResult && acceptResult.accepted === 3);
ok('accept result.skipped === 0', acceptResult && acceptResult.skipped === 0);
ok('accept result.statuses array length 3',
   acceptResult && Array.isArray(acceptResult.statuses) && acceptResult.statuses.length === 3);
ok('all 3 statuses ok',
   acceptResult && acceptResult.statuses.every(s => s.ok === true));
ok('state.manualGroups now has 3 entries',
   sbx.state.manualGroups.length === 3);
ok('first manual group has 2 members (band 0)',
   sbx.state.manualGroups[0] && sbx.state.manualGroups[0].members.length === 2);

// ============================================================================
// 5. Tab auto-jump to inheritance
// ============================================================================
console.log('\n=== 5. Tab auto-jump to inheritance ===');

// gPanelOpen=false → no jump. gPanelOpen=true → jump.
sbx.state.gPanelOpen = false;
sbx.state.gPanelTab  = 'karyotype';
vm.runInContext('_gpKaryoAcceptProposals(' + JSON.stringify(bundle) + ')', sbx);
ok('panel closed → tab unchanged',
   sbx.state.gPanelTab === 'karyotype');

sbx.state.gPanelOpen = true;
sbx.state.gPanelTab  = 'karyotype';
sbx.state.manualGroups = [];
vm.runInContext('_gpKaryoAcceptProposals(' + JSON.stringify(bundle) + ')', sbx);
ok('panel open → tab switched to inheritance',
   sbx.state.gPanelTab === 'inheritance');

// ============================================================================
// 6. Empty/null inputs handled gracefully
// ============================================================================
console.log('\n=== 6. Empty/null inputs ===');

const emptyAccept = vm.runInContext('_gpKaryoAcceptProposals(null)', sbx);
ok('null bundle → accepted=0, skipped=0',
   emptyAccept && emptyAccept.accepted === 0 && emptyAccept.skipped === 0);

const emptyAccept2 = vm.runInContext(
  '_gpKaryoAcceptProposals({proposals: []})', sbx);
ok('empty proposals → accepted=0',
   emptyAccept2 && emptyAccept2.accepted === 0);

// Candidate with no locked_labels
const candNoLabels = vm.runInContext(
  '_gpKaryoProposeAll({id: "X", K: 3})', sbx);
ok('candidate without locked_labels → null bundle', candNoLabels === null);

// ============================================================================
// 7. Naming — H-label suffix logic (without classifier in sandbox)
// ============================================================================
console.log('\n=== 7. Naming sanity ===');

ok('proposal names have <candId>_ prefix',
   props.every(p => p.name.startsWith('I3_')));
ok('proposal names are filesystem-safe (no spaces in candTag)',
   props.every(p => /^[A-Za-z0-9_/-]+$/.test(p.name.split('_')[0])));

// ============================================================================
// 8. Re-render of modal called on tab jump
// ============================================================================
console.log('\n=== 8. Re-render hook fires ===');

// Use an object on the sandbox so the function in vm-context can mutate it
sbx._renderCalls = { count: 0 };
vm.runInContext('function _renderGPanelModal() { _renderCalls.count++; }', sbx);

sbx.state.gPanelOpen = true;
sbx.state.gPanelTab  = 'karyotype';
sbx.state.manualGroups = [];
const renderBefore = sbx._renderCalls.count;
vm.runInContext('_gpKaryoAcceptProposals(' + JSON.stringify(bundle) + ')', sbx);
ok('_renderGPanelModal called after accept (with panel open)',
   sbx._renderCalls.count > renderBefore);

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail > 0 ? 1 : 0);
