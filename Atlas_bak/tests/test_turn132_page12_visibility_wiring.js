// =============================================================================
// turn 132 — Slice 2: page 12 (θπ scrubber) panel-visibility wiring.
//
// Builds on Slice 1's DOM scaffold. Adds _refreshThetaPiPanelVisibility()
// which reads state.layersPresent and toggles each panel's display
// based on which θπ layer(s) are present:
//
//   theta_pi_per_window   → #thLinesPanel, #thTrackedSamplesPanelCompact
//   theta_pi_local_pca    → #thSimPanel, #thZPanel, #thPcaPanel, #thL3Panel
//   theta_pi_envelopes    → #thAnchorStripPanel
//   cusum_theta           → #thCusumHeroPanel
//   ANY θπ layer          → #thCtrlBar shows, #thetaPiEmpty hides
//
// Renderers (canvas drawing) ship in later slices.
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
// 1. Source-level checks — function defined and wired
// ============================================================================
console.log('\n=== 1. Source-level checks ===');

ok('_refreshThetaPiPanelVisibility function defined',
   /function _refreshThetaPiPanelVisibility\(\)/.test(html));

ok('function lives next to _refreshGhslLayerStatus / hasLayer family',
   /_refreshGhslLayerStatus[\s\S]{0,2200}function _refreshThetaPiPanelVisibility/.test(html));

ok('function fires in the onDataLoad refresh path',
   /_refreshThetaPiLayerStatus[\s\S]{0,500}_refreshThetaPiPanelVisibility\(\)/.test(html));

ok('function uses the same hasLayer-style API (state.layersPresent.has)',
   /_refreshThetaPiPanelVisibility[\s\S]{0,1200}state\.layersPresent && state\.layersPresent\.has/.test(html));

ok('function early-returns when document is undefined (vm-safe)',
   /_refreshThetaPiPanelVisibility[\s\S]{0,400}typeof document === 'undefined'/.test(html));

// ============================================================================
// 2. Source-level checks — panel-to-layer mapping
// ============================================================================
console.log('\n=== 2. Panel-to-layer mapping documented in source ===');

const fnBody = (() => {
  const m = html.match(/function _refreshThetaPiPanelVisibility\([\s\S]*?\n\}/);
  return m ? m[0] : '';
})();

ok('function body extracted', fnBody.length > 0);

ok("source maps theta_pi_per_window → #thLinesPanel + #thTrackedSamplesPanelCompact",
   /has\('theta_pi_per_window'\)/.test(fnBody) &&
   /showHide\('thLinesPanel',\s*hasPerWindow/.test(fnBody) &&
   /showHide\('thTrackedSamplesPanelCompact',\s*hasPerWindow/.test(fnBody));

ok("source maps theta_pi_local_pca → #thSimPanel/#thZPanel/#thPcaPanel/#thL3Panel",
   /has\('theta_pi_local_pca'\)/.test(fnBody) &&
   /showHide\('thSimPanel',\s*hasLocalPCA/.test(fnBody) &&
   /showHide\('thZPanel',\s*hasLocalPCA/.test(fnBody) &&
   /showHide\('thPcaPanel',\s*hasLocalPCA/.test(fnBody) &&
   /showHide\('thL3Panel',\s*hasLocalPCA/.test(fnBody));

ok("source maps theta_pi_envelopes → #thAnchorStripPanel",
   /has\('theta_pi_envelopes'\)/.test(fnBody) &&
   /showHide\('thAnchorStripPanel',\s*hasEnvelopes/.test(fnBody));

ok("source maps cusum_theta → #thCusumHeroPanel",
   /has\('cusum_theta'\)/.test(fnBody) &&
   /showHide\('thCusumHeroPanel',\s*hasCusumTheta/.test(fnBody));

ok("source uses 'flex' for #thCtrlBar (it's a flex row)",
   /showHide\('thCtrlBar',\s*hasAny,\s*'flex'\)/.test(fnBody));

ok("source uses 'block' for the .panel divs",
   /showHide\('thCusumHeroPanel',\s*hasCusumTheta,\s*'block'\)/.test(fnBody));

ok("source toggles #thetaPiEmpty inversely to hasAny",
   /getElementById\('thetaPiEmpty'\)[\s\S]{0,300}hasAny \? 'none' : 'block'/.test(fnBody));

// ============================================================================
// 3. Behavioral test — exec the function in a JSDOM-like sandbox
// ============================================================================
console.log('\n=== 3. Behavioral — exec in sandbox ===');

// Build a minimal DOM for the page 12 panels. We don't need real layout —
// just elements with the right ids and a "style" object that records the
// last assigned `display` value.
function makeStyleObj() {
  const s = { display: 'none' };
  return s;
}
function makeElement(id, initialDisplay) {
  return {
    id,
    style: { display: initialDisplay !== undefined ? initialDisplay : 'none' },
  };
}

const panelIds = [
  'thCtrlBar',
  'thCusumHeroPanel',
  'thSimPanel',
  'thZPanel',
  'thLinesPanel',
  'thAnchorStripPanel',
  'thPcaPanel',
  'thTrackedSamplesPanelCompact',
  'thL3Panel',
  'thetaPiEmpty',
];

function makeFakeDoc() {
  const elements = {};
  for (const id of panelIds) {
    elements[id] = makeElement(id, id === 'thetaPiEmpty' ? 'block' : 'none');
  }
  return {
    elements,
    getElementById(id) { return this.elements[id] || null; },
  };
}

function runFnIn(sandboxLike) {
  // Extract function source and run it in a sandbox where document and
  // state are provided.
  const fnSrc = html.match(/function _refreshThetaPiPanelVisibility\([\s\S]*?\n\}/)[0];
  const ctx = vm.createContext({
    document: sandboxLike.document,
    state: sandboxLike.state,
  });
  vm.runInContext(fnSrc + '\n_refreshThetaPiPanelVisibility();', ctx);
}

// 3a. No layers present → ALL panels hidden, empty-state visible
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set() },
  };
  runFnIn(sandbox);

  ok('3a: empty layers → #thCtrlBar hidden',
     doc.elements.thCtrlBar.style.display === 'none');
  ok('3a: empty layers → all .panel mirrors hidden',
     ['thCusumHeroPanel', 'thSimPanel', 'thZPanel', 'thLinesPanel',
      'thAnchorStripPanel', 'thPcaPanel', 'thTrackedSamplesPanelCompact',
      'thL3Panel'].every(id => doc.elements[id].style.display === 'none'));
  ok('3a: empty layers → #thetaPiEmpty visible',
     doc.elements.thetaPiEmpty.style.display === 'block');
}

// 3b. Only theta_pi_per_window → #thLinesPanel + #thTrackedSamples reveal
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set(['theta_pi_per_window']) },
  };
  runFnIn(sandbox);

  ok('3b: per_window only → #thCtrlBar visible (any-layer)',
     doc.elements.thCtrlBar.style.display === 'flex');
  ok('3b: per_window only → #thLinesPanel visible',
     doc.elements.thLinesPanel.style.display === 'block');
  ok('3b: per_window only → #thTrackedSamplesPanelCompact visible',
     doc.elements.thTrackedSamplesPanelCompact.style.display === 'block');
  ok('3b: per_window only → #thSimPanel still hidden (needs local_pca)',
     doc.elements.thSimPanel.style.display === 'none');
  ok('3b: per_window only → #thCusumHeroPanel still hidden (needs cusum_theta)',
     doc.elements.thCusumHeroPanel.style.display === 'none');
  ok('3b: per_window only → #thAnchorStripPanel still hidden (needs envelopes)',
     doc.elements.thAnchorStripPanel.style.display === 'none');
  ok('3b: per_window only → #thetaPiEmpty hidden',
     doc.elements.thetaPiEmpty.style.display === 'none');
}

// 3c. Only theta_pi_local_pca → 4 panels reveal
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set(['theta_pi_local_pca']) },
  };
  runFnIn(sandbox);

  ok('3c: local_pca only → #thSimPanel, #thZPanel, #thPcaPanel, #thL3Panel visible',
     ['thSimPanel', 'thZPanel', 'thPcaPanel', 'thL3Panel']
       .every(id => doc.elements[id].style.display === 'block'));
  ok('3c: local_pca only → #thLinesPanel hidden (needs per_window)',
     doc.elements.thLinesPanel.style.display === 'none');
  ok('3c: local_pca only → #thAnchorStripPanel hidden',
     doc.elements.thAnchorStripPanel.style.display === 'none');
  ok('3c: local_pca only → #thetaPiEmpty hidden',
     doc.elements.thetaPiEmpty.style.display === 'none');
}

// 3d. Only theta_pi_envelopes → only anchor strip + ctrl bar
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set(['theta_pi_envelopes']) },
  };
  runFnIn(sandbox);

  ok('3d: envelopes only → #thAnchorStripPanel visible',
     doc.elements.thAnchorStripPanel.style.display === 'block');
  ok('3d: envelopes only → #thCtrlBar visible',
     doc.elements.thCtrlBar.style.display === 'flex');
  ok('3d: envelopes only → #thSimPanel hidden (needs local_pca)',
     doc.elements.thSimPanel.style.display === 'none');
  ok('3d: envelopes only → #thLinesPanel hidden (needs per_window)',
     doc.elements.thLinesPanel.style.display === 'none');
  ok('3d: envelopes only → #thCusumHeroPanel hidden (needs cusum_theta)',
     doc.elements.thCusumHeroPanel.style.display === 'none');
}

// 3e. Only cusum_theta → only hero + ctrl bar
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set(['cusum_theta']) },
  };
  runFnIn(sandbox);

  ok('3e: cusum_theta only → #thCusumHeroPanel visible',
     doc.elements.thCusumHeroPanel.style.display === 'block');
  ok('3e: cusum_theta only → #thCtrlBar visible',
     doc.elements.thCtrlBar.style.display === 'flex');
  ok('3e: cusum_theta only → other 7 panels hidden',
     ['thSimPanel', 'thZPanel', 'thLinesPanel', 'thAnchorStripPanel',
      'thPcaPanel', 'thTrackedSamplesPanelCompact', 'thL3Panel']
       .every(id => doc.elements[id].style.display === 'none'));
  ok('3e: cusum_theta only → #thetaPiEmpty hidden',
     doc.elements.thetaPiEmpty.style.display === 'none');
}

// 3f. ALL four θπ layers → ALL 9 panels visible, empty hidden
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set([
      'theta_pi_per_window', 'theta_pi_local_pca',
      'theta_pi_envelopes', 'cusum_theta',
    ]) },
  };
  runFnIn(sandbox);

  ok('3f: all layers → #thCtrlBar visible (flex)',
     doc.elements.thCtrlBar.style.display === 'flex');
  ok('3f: all layers → all 8 .panel mirrors visible (block)',
     ['thCusumHeroPanel', 'thSimPanel', 'thZPanel', 'thLinesPanel',
      'thAnchorStripPanel', 'thPcaPanel', 'thTrackedSamplesPanelCompact',
      'thL3Panel'].every(id => doc.elements[id].style.display === 'block'));
  ok('3f: all layers → #thetaPiEmpty hidden',
     doc.elements.thetaPiEmpty.style.display === 'none');
}

// 3g. Idempotency — calling twice with same state gives same result
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: { layersPresent: new Set(['theta_pi_local_pca']) },
  };
  runFnIn(sandbox);
  const snapshot1 = panelIds.map(id => doc.elements[id].style.display).join(',');
  runFnIn(sandbox);
  const snapshot2 = panelIds.map(id => doc.elements[id].style.display).join(',');
  ok('3g: idempotent — second call yields same display states',
     snapshot1 === snapshot2,
     'snap1=' + snapshot1 + ' snap2=' + snapshot2);
}

// 3h. Layer removal — going from many layers to none re-hides panels
{
  const doc = makeFakeDoc();
  // First call with all layers
  let sandbox = {
    document: doc,
    state: { layersPresent: new Set([
      'theta_pi_per_window', 'theta_pi_local_pca',
      'theta_pi_envelopes', 'cusum_theta',
    ]) },
  };
  runFnIn(sandbox);
  // Now layers go away (e.g. user loaded a different chrom with no θπ data)
  sandbox = {
    document: doc,
    state: { layersPresent: new Set() },
  };
  runFnIn(sandbox);

  ok('3h: layers removed → all panels hidden again',
     ['thCtrlBar', 'thCusumHeroPanel', 'thSimPanel', 'thZPanel',
      'thLinesPanel', 'thAnchorStripPanel', 'thPcaPanel',
      'thTrackedSamplesPanelCompact', 'thL3Panel']
       .every(id => doc.elements[id].style.display === 'none'));
  ok('3h: layers removed → #thetaPiEmpty visible again',
     doc.elements.thetaPiEmpty.style.display === 'block');
}

// 3i. Defensive: state without layersPresent doesn't crash
{
  const doc = makeFakeDoc();
  const sandbox = {
    document: doc,
    state: {},  // no layersPresent at all
  };
  let crashed = false;
  try { runFnIn(sandbox); } catch (e) { crashed = true; }
  ok('3i: state without layersPresent does not crash', !crashed);
  ok('3i: with no layersPresent → #thetaPiEmpty visible',
     doc.elements.thetaPiEmpty.style.display === 'block');
  ok('3i: with no layersPresent → all panels hidden',
     ['thCtrlBar', 'thCusumHeroPanel', 'thSimPanel']
       .every(id => doc.elements[id].style.display === 'none'));
}

// 3j. Defensive: missing DOM elements don't crash
{
  const sparse = {
    elements: { thSimPanel: makeElement('thSimPanel', 'none') },  // only one element exists
    getElementById(id) { return this.elements[id] || null; },
  };
  const sandbox = {
    document: sparse,
    state: { layersPresent: new Set(['theta_pi_local_pca']) },
  };
  let crashed = false;
  try {
    const fnSrc = html.match(/function _refreshThetaPiPanelVisibility\([\s\S]*?\n\}/)[0];
    const ctx = vm.createContext({ document: sandbox.document, state: sandbox.state });
    vm.runInContext(fnSrc + '\n_refreshThetaPiPanelVisibility();', ctx);
  } catch (e) { crashed = true; }
  ok('3j: missing DOM elements do not crash', !crashed);
  ok('3j: existing element still updated (#thSimPanel → block)',
     sparse.elements.thSimPanel.style.display === 'block');
}

// ============================================================================
// 4. Existing layer-status indicators still work (no regression)
// ============================================================================
console.log('\n=== 4. Existing _refreshThetaPiLayerStatus untouched ===');

ok('_refreshThetaPiLayerStatus still defined',
   /function _refreshThetaPiLayerStatus\(\)/.test(html));

ok('_refreshThetaPiLayerStatus still iterates [data-th-layer]',
   /_refreshThetaPiLayerStatus[\s\S]{0,600}querySelectorAll\('\[data-th-layer\]'\)/.test(html));

ok('cusum_theta indicator row still in #thetaPiLayerStatus grid',
   /id="thetaPiLayerStatus"[\s\S]{0,3500}data-th-layer="cusum_theta"/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
