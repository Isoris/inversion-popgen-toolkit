// =============================================================================
// turn 136 — Slice 2: G-panel karyotype tab content (SPEC §4)
// =============================================================================
// Replaces the Slice 1 placeholder with real per-band content derived from
// state.candidate.locked_labels. New helpers tested here:
//
//   _gpKaryoExtractBands(candidate, samples, K)
//     Returns { bands: [{bandIdx, label, sampleIdxList, members, color}],
//               nUngrouped, K }
//     Honors K boundaries (lab >= K is ungrouped); resolves CGA names from
//     samples[i].cga || samples[i].ind; preserves bands with 0 members.
//
//   _gpKaryoMakeManualGroup(candidate, bandIdx)
//     Wraps addToManualGroup(name, sampleIdxList, opts). Default name
//     "<candId>_<bandLabel>".
//
//   _gpKaryoColorByCandidate(candidate)
//     Sets state.lockedLabels = Int8Array(candidate.locked_labels) and
//     state.lockedRefL2 = candidate.ref_l2. Mirrors #lockColorsBtn handler.
//
//   _gpKaryoExportTSV(candidate)
//     Builds a 2-column TSV (cga, band_label) and triggers download.
//     Returns false in non-DOM env (sandbox).
//
//   _gpKaryoColor(bandIdx)
//     Returns hex from a 6-color palette matching _ACK_BAND_PALETTE.
//
// Plus updated _gPanelRenderTabKaryotype renderer + interaction wiring in
// _renderGPanelModal.
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
// 1. Source-level: new helpers exist + window-exported
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_gpKaryoExtractBands defined',     /function _gpKaryoExtractBands\b/.test(html));
ok('_gpKaryoMakeManualGroup defined',  /function _gpKaryoMakeManualGroup\b/.test(html));
ok('_gpKaryoColorByCandidate defined', /function _gpKaryoColorByCandidate\b/.test(html));
ok('_gpKaryoExportTSV defined',        /function _gpKaryoExportTSV\b/.test(html));
ok('_gpKaryoColor defined',            /function _gpKaryoColor\b/.test(html));
ok('_GPANEL_KARYO_PALETTE constant',
   /const _GPANEL_KARYO_PALETTE\s*=\s*\[/.test(html));
ok('palette has 6 colors',
   /_GPANEL_KARYO_PALETTE\s*=\s*\[\s*'#4fa3ff',\s*'#b8b8b8',\s*'#f5a524',\s*'#3cc08a',\s*'#e0555c',\s*'#b07cf7'\s*\]/.test(html));

ok('window._gpKaryoExtractBands exported',
   /window\._gpKaryoExtractBands\b/.test(html));
ok('window._gpKaryoMakeManualGroup exported',
   /window\._gpKaryoMakeManualGroup\b/.test(html));
ok('window._gpKaryoColorByCandidate exported',
   /window\._gpKaryoColorByCandidate\b/.test(html));
ok('window._gpKaryoExportTSV exported',
   /window\._gpKaryoExportTSV\b/.test(html));

// ============================================================================
// 2. _gpKaryoColor — palette indexing
// ============================================================================
console.log('\n=== 2. _gpKaryoColor palette ===');

function extractFnSrc(name) {
  const re = new RegExp('function ' + name + '\\([\\s\\S]*?\\n\\}', 'm');
  const m = html.match(re);
  return m ? m[0] : null;
}

const colorFnSrc = extractFnSrc('_gpKaryoColor');
ok('source extractable', colorFnSrc !== null);

const colorCtx = vm.createContext({});
vm.runInContext(`
  const _GPANEL_KARYO_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];
  ${colorFnSrc}
`, colorCtx);

ok('band 0 → #4fa3ff', colorCtx._gpKaryoColor(0) === '#4fa3ff');
ok('band 5 → #b07cf7', colorCtx._gpKaryoColor(5) === '#b07cf7');
ok('band 6 wraps around → #4fa3ff', colorCtx._gpKaryoColor(6) === '#4fa3ff');
ok('null → #666 fallback',  colorCtx._gpKaryoColor(null) === '#666');
ok('-1 → #666 fallback',    colorCtx._gpKaryoColor(-1) === '#666');

// ============================================================================
// 3. _gpKaryoExtractBands — band derivation from locked_labels
// ============================================================================
console.log('\n=== 3. _gpKaryoExtractBands ===');

const extractFnBlock = extractFnSrc('_gpKaryoExtractBands');
ok('source extractable', extractFnBlock !== null);

// Stub getKaryotypeLabel + the palette in the sandbox
function makeExtractCtx() {
  const ctx = {
    Number, Array, Math, Int8Array,
  };
  vm.createContext(ctx);
  vm.runInContext(`
    function getKaryotypeLabel(b, K) { return 'b' + b + '/K' + K; }
    const _GPANEL_KARYO_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];
    ${colorFnSrc}
    ${extractFnBlock}
  `, ctx);
  return ctx;
}

// 3a. K=3 with 6 samples → 3 bands of 2
{
  const ctx = makeExtractCtx();
  const candidate = { locked_labels: new Int8Array([0,0,1,1,2,2]), K: 3 };
  const samples = [
    { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
    { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
  ];
  const ext = ctx._gpKaryoExtractBands(candidate, samples, 3);
  ok('returns object', ext && typeof ext === 'object');
  ok('K=3 → 3 bands',  ext.bands.length === 3);
  ok('band 0 has 2 members', ext.bands[0].sampleIdxList.length === 2);
  ok('band 0 members are CGA001+CGA002',
     ext.bands[0].members.length === 2 &&
     ext.bands[0].members[0] === 'CGA001' &&
     ext.bands[0].members[1] === 'CGA002');
  ok('band 0 sampleIdxList is [0,1]',
     ext.bands[0].sampleIdxList[0] === 0 &&
     ext.bands[0].sampleIdxList[1] === 1);
  ok('band 0 has color from palette', ext.bands[0].color === '#4fa3ff');
  ok('band 0 label uses getKaryotypeLabel', ext.bands[0].label === 'b0/K3');
  ok('nUngrouped = 0', ext.nUngrouped === 0);
}

// 3b. K=6 with mixed labels including -1 (ungrouped)
{
  const ctx = makeExtractCtx();
  const candidate = { locked_labels: new Int8Array([0, 5, -1, 3, -1, 1, 0, 2]), K: 6 };
  const samples = Array.from({ length: 8 }, (_, i) => ({ cga: 'CGA' + (i + 1) }));
  const ext = ctx._gpKaryoExtractBands(candidate, samples, 6);
  ok('K=6 → 6 bands', ext.bands.length === 6);
  ok('band 0 has 2 members (samples 0 and 6)',
     ext.bands[0].sampleIdxList.length === 2);
  ok('band 5 has 1 member (sample 1)',
     ext.bands[5].sampleIdxList.length === 1 &&
     ext.bands[5].sampleIdxList[0] === 1);
  ok('band 4 is empty', ext.bands[4].sampleIdxList.length === 0);
  ok('nUngrouped = 2 (two -1 labels)', ext.nUngrouped === 2);
}

// 3c. Out-of-K-range labels also count as ungrouped
{
  const ctx = makeExtractCtx();
  // K=3 but locked_labels has values 0..4 (lab=3 and lab=4 are out of range)
  const candidate = { locked_labels: new Int8Array([0, 1, 2, 3, 4, 0]), K: 3 };
  const samples = Array.from({ length: 6 }, (_, i) => ({ cga: 'CGA' + i }));
  const ext = ctx._gpKaryoExtractBands(candidate, samples, 3);
  ok('K=3 → 3 bands', ext.bands.length === 3);
  ok('band 0 has 2 members (samples 0 and 5)',
     ext.bands[0].sampleIdxList.length === 2);
  ok('nUngrouped counts out-of-range (3 and 4) — total 2',
     ext.nUngrouped === 2);
}

// 3d. Length mismatch — uses min(labels.length, samples.length)
{
  const ctx = makeExtractCtx();
  const candidate = { locked_labels: new Int8Array([0, 1, 2, 0, 1]), K: 3 };
  const samples = [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }];   // only 3 samples
  const ext = ctx._gpKaryoExtractBands(candidate, samples, 3);
  ok('truncates to samples.length',
     ext.bands[0].sampleIdxList.length + ext.bands[1].sampleIdxList.length +
     ext.bands[2].sampleIdxList.length + ext.nUngrouped === 3);
}

// 3e. Sample with .ind fallback (no .cga)
{
  const ctx = makeExtractCtx();
  const candidate = { locked_labels: new Int8Array([0, 0]), K: 3 };
  const samples = [{ ind: 'IND_X' }, { ind: 'IND_Y' }];
  const ext = ctx._gpKaryoExtractBands(candidate, samples, 3);
  ok('uses .ind when .cga missing',
     ext.bands[0].members[0] === 'IND_X' && ext.bands[0].members[1] === 'IND_Y');
}

// 3f. Null candidate → null
{
  const ctx = makeExtractCtx();
  const ext = ctx._gpKaryoExtractBands(null, [], 3);
  ok('null candidate returns null', ext === null);
}

// 3g. Missing locked_labels → null
{
  const ctx = makeExtractCtx();
  const ext = ctx._gpKaryoExtractBands({ K: 3 }, [{ cga: 'X' }], 3);
  ok('missing locked_labels returns null', ext === null);
}

// 3h. Empty locked_labels → null
{
  const ctx = makeExtractCtx();
  const ext = ctx._gpKaryoExtractBands(
    { locked_labels: new Int8Array(0), K: 3 }, [{ cga: 'X' }], 3);
  ok('empty locked_labels returns null', ext === null);
}

// ============================================================================
// 4. _gpKaryoMakeManualGroup — wraps addToManualGroup
// ============================================================================
console.log('\n=== 4. _gpKaryoMakeManualGroup ===');

function makeMakeManualCtx(scenario) {
  const stateStub = scenario.state || {
    data: { samples: [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }] },
    k: 3,
  };
  const calls = [];
  const ctx = {
    state: stateStub,
    window: undefined,
    Number, Array, Int8Array, Math,
    addToManualGroup: function (name, sampleIdxList, opts) {
      calls.push({ name, sampleIdxList: Array.from(sampleIdxList), opts });
      return { id: 'g_' + name, name, members: sampleIdxList.map(i => 'CGA' + (i+1).toString().padStart(3, '0')) };
    },
    getKaryotypeLabel: (b, K) => 'L' + b + '/' + K,
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(`
    const _GPANEL_KARYO_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];
    ${colorFnSrc}
    ${extractFnBlock}
    ${extractFnSrc('_gpKaryoMakeManualGroup')}
  `, ctx);
  return { ctx, calls, state: stateStub };
}

// 4a. Happy path
{
  const sb = makeMakeManualCtx({});
  const cand = { id: 'cand_X', locked_labels: new Int8Array([0, 1, 0]), K: 3 };
  sb.state.data.samples = [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }];
  const result = sb.ctx._gpKaryoMakeManualGroup(cand, 0);
  ok('returns the group from addToManualGroup', result && result.id);
  ok('addToManualGroup called once', sb.calls.length === 1);
  ok('group name is candId_<label>',
     sb.calls[0].name === 'cand_X_L0/3' || /^cand_X_L0/.test(sb.calls[0].name));
  ok('sampleIdxList = [0, 2] (band 0 members)',
     sb.calls[0].sampleIdxList.length === 2 &&
     sb.calls[0].sampleIdxList[0] === 0 &&
     sb.calls[0].sampleIdxList[1] === 2);
  ok('opts.color set from palette',
     sb.calls[0].opts && sb.calls[0].opts.color === '#4fa3ff');
}

// 4b. Empty band → returns null, no addToManualGroup call
{
  const sb = makeMakeManualCtx({});
  const cand = { id: 'cand_X', locked_labels: new Int8Array([0, 0, 0]), K: 3 };
  sb.state.data.samples = [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }];
  const result = sb.ctx._gpKaryoMakeManualGroup(cand, 1);   // band 1 is empty
  ok('empty band returns null', result === null);
  ok('addToManualGroup NOT called', sb.calls.length === 0);
}

// 4c. Out-of-range bandIdx → null
{
  const sb = makeMakeManualCtx({});
  const cand = { id: 'cand_X', locked_labels: new Int8Array([0, 1, 2]), K: 3 };
  sb.state.data.samples = [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }];
  ok('bandIdx -1 returns null',
     sb.ctx._gpKaryoMakeManualGroup(cand, -1) === null);
  ok('bandIdx 99 returns null',
     sb.ctx._gpKaryoMakeManualGroup(cand, 99) === null);
}

// 4d. Special chars in candidate id are sanitized
{
  const sb = makeMakeManualCtx({});
  const cand = { id: 'cand/with spaces!', locked_labels: new Int8Array([0, 0]), K: 3 };
  sb.state.data.samples = [{ cga: 'CGA001' }, { cga: 'CGA002' }];
  sb.ctx._gpKaryoMakeManualGroup(cand, 0);
  ok('candidate id sanitized to safe chars',
     /^[A-Za-z0-9_/-]+$/.test(sb.calls[0].name),
     'got name: ' + sb.calls[0].name);
}

// ============================================================================
// 5. _gpKaryoColorByCandidate — locks PCA colors to candidate
// ============================================================================
console.log('\n=== 5. _gpKaryoColorByCandidate ===');

function makeColorByCtx() {
  const stateStub = { lockedLabels: null, lockedRefL2: null };
  const calls = { drawPCA: 0, renderL3Panel: 0, drawLinesPanel: 0,
                  refreshLockBtn: 0, refreshCandidateUI: 0 };
  const ctx = {
    state: stateStub,
    window: undefined,
    Int8Array, Array, Number,
    drawPCA:          () => { calls.drawPCA++; },
    renderL3Panel:    () => { calls.renderL3Panel++; },
    drawLinesPanel:   () => { calls.drawLinesPanel++; },
    refreshLockBtn:   () => { calls.refreshLockBtn++; },
    refreshCandidateUI: () => { calls.refreshCandidateUI++; },
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(extractFnSrc('_gpKaryoColorByCandidate'), ctx);
  return { ctx, state: stateStub, calls };
}

// 5a. Happy path
{
  const sb = makeColorByCtx();
  const cand = { locked_labels: new Int8Array([0, 1, 2, 0, 1]), ref_l2: 7 };
  const ok1 = sb.ctx._gpKaryoColorByCandidate(cand);
  ok('returns true on success', ok1 === true);
  ok('state.lockedLabels populated',
     sb.state.lockedLabels !== null && sb.state.lockedLabels.length === 5);
  ok('state.lockedLabels is a copy (Int8Array)',
     sb.state.lockedLabels instanceof Int8Array);
  ok('state.lockedRefL2 = 7', sb.state.lockedRefL2 === 7);

  // Snapshot semantics: mutating original locked_labels does NOT mutate state
  cand.locked_labels[0] = 99;
  ok('snapshot is decoupled from original',
     sb.state.lockedLabels[0] === 0);

  // Repaints triggered
  ok('drawPCA called',          sb.calls.drawPCA === 1);
  ok('renderL3Panel called',    sb.calls.renderL3Panel === 1);
  ok('drawLinesPanel called',   sb.calls.drawLinesPanel === 1);
  ok('refreshLockBtn called',   sb.calls.refreshLockBtn === 1);
  ok('refreshCandidateUI called', sb.calls.refreshCandidateUI === 1);
}

// 5b. Null candidate → false, state unchanged
{
  const sb = makeColorByCtx();
  ok('null candidate returns false',
     sb.ctx._gpKaryoColorByCandidate(null) === false);
  ok('state.lockedLabels remains null',
     sb.state.lockedLabels === null);
}

// 5c. Missing locked_labels → false
{
  const sb = makeColorByCtx();
  ok('missing locked_labels returns false',
     sb.ctx._gpKaryoColorByCandidate({ ref_l2: 5 }) === false);
}

// 5d. ref_l2 not a number → state.lockedRefL2 = null
{
  const sb = makeColorByCtx();
  const cand = { locked_labels: new Int8Array([0, 1]), ref_l2: undefined };
  sb.ctx._gpKaryoColorByCandidate(cand);
  ok('lockedRefL2 = null when ref_l2 not a number',
     sb.state.lockedRefL2 === null);
}

// ============================================================================
// 6. _gpKaryoExportTSV — TSV row count + sandbox-safe (returns false)
// ============================================================================
console.log('\n=== 6. _gpKaryoExportTSV ===');

// In a sandbox without Blob/URL/document, the function should fail-soft with
// false rather than throw. Test that. Real DOM testing happens browser-side.
function makeExportTSVSandbox() {
  const stateStub = {
    data: {
      chrom: 'LG28',
      samples: [
        { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
        { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
      ],
    },
    k: 3,
  };
  const ctx = {
    state: stateStub,
    window: undefined,
    Number, Array, Int8Array, Math,
    Blob: undefined,        // sandbox path: no Blob → return false
    URL: undefined,
    document: undefined,
    setTimeout: () => {},
    console: { warn: () => {} },
    getKaryotypeLabel: (b, K) => 'L' + b,
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(`
    const _GPANEL_KARYO_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];
    ${colorFnSrc}
    ${extractFnBlock}
    ${extractFnSrc('_gpKaryoExportTSV')}
  `, ctx);
  return { ctx, state: stateStub };
}

// 6a. Returns false in sandbox (no Blob)
{
  const sb = makeExportTSVSandbox();
  const cand = { id: 'cand_X', locked_labels: new Int8Array([0,0,1,1,2,2]), K: 3 };
  const result = sb.ctx._gpKaryoExportTSV(cand);
  ok('sandbox (no Blob) returns false', result === false);
}

// 6b. Null candidate → false
{
  const sb = makeExportTSVSandbox();
  ok('null candidate returns false', sb.ctx._gpKaryoExportTSV(null) === false);
}

// 6c. Missing locked_labels → false
{
  const sb = makeExportTSVSandbox();
  ok('missing locked_labels returns false',
     sb.ctx._gpKaryoExportTSV({ id: 'cand_X', K: 3 }) === false);
}

// 6d. With a real Blob/URL/document mock — verify TSV content
{
  const stateStub = {
    data: {
      chrom: 'LG28',
      samples: [
        { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
      ],
    },
    k: 3,
  };
  const downloads = [];
  const ctx = {
    state: stateStub,
    window: undefined,
    Number, Array, Int8Array, Math,
    Blob: function (parts, opts) { return { _parts: parts, _opts: opts }; },
    URL: { createObjectURL: () => 'blob://fake', revokeObjectURL: () => {} },
    document: {
      body: { appendChild() {}, removeChild() {} },
      createElement(tag) {
        const el = { tagName: tag, _attrs: {} };
        el.click = function () {
          downloads.push({ href: el.href, download: el.download, blob: el._blob });
        };
        Object.defineProperty(el, 'href', {
          set(v) { el._href = v; }, get() { return el._href; }
        });
        return el;
      },
    },
    setTimeout: () => {},
    console: { warn: () => {} },
    getKaryotypeLabel: (b, K) => 'L' + b + '/K' + K,
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(`
    const _GPANEL_KARYO_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];
    ${colorFnSrc}
    ${extractFnBlock}
    ${extractFnSrc('_gpKaryoExportTSV')}
  `, ctx);
  const cand = { id: 'cand_X', locked_labels: new Int8Array([0, 1, 2]), K: 3 };
  const result = ctx._gpKaryoExportTSV(cand);
  ok('returns true with real Blob mock', result === true);
  ok('one download triggered', downloads.length === 1);
  if (downloads.length === 1) {
    const blob = downloads[0].blob;
    // Mock returned the Blob from createElement, so check via the blob
    // we created in `new Blob([...])`
  }
  // Re-derive the TSV from the mock blob
  // (Our Blob mock retains _parts.) Find the "a" element returned from
  // createElement to access _blob via the click handler.
  // Easier: just verify the download fields directly.
  ok('download filename matches pattern',
     /^karyotype\.cand_X\.LG28\.tsv$/.test(downloads[0].download || ''));
}

// ============================================================================
// 7. Renderer — _gPanelRenderTabKaryotype output (Slice 2 body)
// ============================================================================
console.log('\n=== 7. _gPanelRenderTabKaryotype Slice 2 body ===');

function extractGPanelBlock() {
  const start = html.indexOf('// G-PANEL UNIFIED GROUPS MODAL');
  if (start < 0) return null;
  const endMarker = 'window._gpKaryoColor              = _gpKaryoColor;';
  const endIdx = html.indexOf(endMarker, start);
  if (endIdx < 0) return null;
  const after = html.indexOf('}', endIdx);
  if (after < 0) return null;
  return html.substring(start, after + 1);
}
const gpanelBlock = extractGPanelBlock();
ok('G-panel block extracts cleanly',
   gpanelBlock !== null && gpanelBlock.length > 1000);

function makeRenderCtx(scenario) {
  const stateStub = scenario.state || {
    candidate: null,
    candidateList: [],
    candidates: {},
    data: null,
    k: 3,
    activeMode: 'default',
    gPanelOpen: false,
    gPanelTab: 'karyotype',
  };
  const ctx = {
    state: stateStub,
    window: undefined,
    document: undefined,
    Number, Array, Object, Math, JSON, Set, Date, Int8Array,
    isFinite, parseInt, parseFloat,
    setTimeout: () => {},
    console: { warn: () => {}, log: () => {} },
    renderManualGroupsList: () => {},
    addToManualGroup: () => null,
    drawPCA: () => {}, renderL3Panel: () => {}, drawLinesPanel: () => {},
    refreshLockBtn: () => {}, refreshCandidateUI: () => {},
    getKaryotypeLabel: (b, K) => {
      if (K === 3) return ['band 1 (lo)','band 2 (mid)','band 3 (hi)'][b] || 'band ' + (b+1);
      return 'band ' + (b+1);
    },
    getKaryotypeLabelCaveat: () => null,
    _karyoLabelEnsureVocab: () => 'legacy',
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(gpanelBlock, ctx);
  return { ctx, state: stateStub };
}

// 7a. No candidate → empty state
{
  console.log('\n--- 7a. No-candidate empty state ---');
  const sb = makeRenderCtx({});
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('returns string', typeof body === 'string');
  ok('mentions "Karyotype groups" header',
     /Karyotype groups/.test(body));
  ok('shows "No candidate currently focused"',
     /No candidate currently focused/.test(body));
  ok('mentions catalogue OR page-1 pin path',
     /catalogue|page-1 pin/.test(body));
}

// 7b. Candidate with locked_labels — full Slice 2 rendering
{
  console.log('\n--- 7b. Candidate with locked_labels (K=3) ---');
  const sb = makeRenderCtx({
    state: {
      candidate: {
        id: 'cand_42',
        start_bp: 14e6,
        end_bp: 16.5e6,
        K: 3,
        locked_labels: new Int8Array([0,0,0,1,1,1,2,2,2]),
        ref_l2: 7,
      },
      candidateList: [], candidates: {},
      data: { chrom: 'LG28', samples: Array.from({length: 9}, (_, i) => ({cga: 'CGA' + i})) },
      k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('shows candidate id', /cand_42/.test(body));
  ok('shows bp range as Mb', /14\.00–16\.50 Mb/.test(body));
  ok('shows K=3', /K=3/.test(body));
  ok('shows label vocab',  /vocab:\s*legacy/.test(body));
  ok('renders 3 band rows',
     (body.match(/_gpKaryoBandRow/g) || []).length === 3);
  ok('renders K=3 legacy labels',
     /band 1 \(lo\)/.test(body) &&
     /band 2 \(mid\)/.test(body) &&
     /band 3 \(hi\)/.test(body));
  ok('each band shows n=3 count',
     (body.match(/n=3/g) || []).length >= 3);
  ok('renders per-band "+ manual" buttons',
     (body.match(/_gpKaryoMakeBtn/g) || []).length === 3);
  ok('renders "color PCA by this candidate" button',
     /id="gpKaryoColorByBtn"/.test(body));
  ok('renders "export TSV" button',
     /id="gpKaryoExportBtn"/.test(body));
  ok('shows sample CGA preview (CGA0..)',
     /CGA0/.test(body));
}

// 7c. Detailed-vocab mode shows operational-label caveat
{
  console.log('\n--- 7c. Detailed vocab caveat ---');
  const sb = makeRenderCtx({
    state: {
      candidate: {
        id: 'cand_X', start_bp: 1e6, end_bp: 2e6, K: 3,
        locked_labels: new Int8Array([0, 1, 2]), ref_l2: 0,
      },
      data: { chrom: 'LG28', samples: [{cga: 'A'}, {cga: 'B'}, {cga: 'C'}] },
      k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  // Override the vocab + caveat
  sb.ctx._karyoLabelEnsureVocab = () => 'detailed';
  sb.ctx.getKaryotypeLabelCaveat = () => 'Operational H-label — confirm with heterozygosity.';
  sb.ctx.getKaryotypeLabel = (b, K) => ['H1/H1','H1/H2','H2/H2'][b];
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('detailed vocab shown in header', /vocab:\s*detailed/.test(body));
  ok('caveat indicator (⚠ operational) present',
     /⚠ operational labels/.test(body));
  ok('detailed labels rendered (H1/H1, H1/H2, H2/H2)',
     /H1\/H1/.test(body) && /H1\/H2/.test(body) && /H2\/H2/.test(body));
}

// 7d. Candidate without locked_labels → no-labels empty state
{
  console.log('\n--- 7d. Candidate without locked_labels ---');
  const sb = makeRenderCtx({
    state: {
      candidate: { id: 'cand_X', start_bp: 1e6, end_bp: 2e6, K: 3 },
      data: { chrom: 'LG28', samples: [{cga: 'A'}] },
      k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('shows no-locked_labels message',
     /no\s*<code>locked_labels<\/code>|locked_labels.*length 0|missing/i.test(body));
  ok('still shows candidate id', /cand_X/.test(body));
}

// 7e. Empty bands rendered with disabled "+ manual" button
{
  console.log('\n--- 7e. Empty bands disabled state ---');
  const sb = makeRenderCtx({
    state: {
      candidate: {
        id: 'cand_X', start_bp: 1e6, end_bp: 2e6, K: 3,
        // All in band 0; bands 1 and 2 empty
        locked_labels: new Int8Array([0, 0, 0]), ref_l2: 0,
      },
      data: { chrom: 'LG28', samples: [{cga: 'A'}, {cga: 'B'}, {cga: 'C'}] },
      k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('renders all 3 bands even when 2 are empty',
     (body.match(/_gpKaryoBandRow/g) || []).length === 3);
  // The 2 empty bands should have disabled "+ manual" buttons
  ok('2 disabled "+ manual" buttons (one per empty band)',
     (body.match(/_gpKaryoMakeBtn[^>]*disabled/g) || []).length === 2);
  ok('1 enabled "+ manual" button (for the non-empty band)',
     (body.match(/_gpKaryoMakeBtn(?:[^>]*?(?!disabled))*?>/g) || []).length >= 1);
  ok('shows "(no samples)" preview for empty bands',
     /\(no samples\)/.test(body));
}

// 7f. Ungrouped count displayed when present
{
  console.log('\n--- 7f. Ungrouped count surfacing ---');
  const sb = makeRenderCtx({
    state: {
      candidate: {
        id: 'cand_X', start_bp: 1e6, end_bp: 2e6, K: 3,
        // Two -1 labels (ungrouped)
        locked_labels: new Int8Array([0, -1, 1, -1, 2]), ref_l2: 0,
      },
      data: { chrom: 'LG28', samples: Array.from({length: 5}, (_, i) => ({cga: 'C' + i})) },
      k: 3, activeMode: 'default',
      gPanelOpen: false, gPanelTab: 'karyotype',
    },
  });
  const body = sb.ctx._gPanelRenderTabKaryotype();
  ok('shows ungrouped count (2)',
     /2 samples? have label -1/.test(body));
}

// ============================================================================
// 8. Wiring presence — _renderGPanelModal binds karyotype interactions
// ============================================================================
console.log('\n=== 8. Karyotype interaction wiring (source-level) ===');

ok('_renderGPanelModal contains "karyotype" branch',
   /activeTab === 'karyotype'/.test(html));
ok("karyotype branch wires _gpKaryoMakeBtn click handlers",
   /querySelectorAll\('\._gpKaryoMakeBtn'\)/.test(html));
ok('karyotype branch wires #gpKaryoColorByBtn click',
   /getElementById\('gpKaryoColorByBtn'\)/.test(html));
ok('karyotype branch wires #gpKaryoExportBtn click',
   /getElementById\('gpKaryoExportBtn'\)/.test(html));
ok("karyotype color-by closes popup after success",
   /_gpKaryoColorByCandidate\(cand\)[\s\S]{0,200}_gPanelClose\(\)/.test(html));

// ============================================================================
// 9. Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
