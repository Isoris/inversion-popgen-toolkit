// =============================================================================
// test_atlas_dock_tabs.js
// =============================================================================
// Tests for the three tab UIs from turn 5b. Reuses the proxy-style DOM mock
// from the dock test, plus mocks for the engine + live layer.
//
// Run: node test_atlas_dock_tabs.js
// =============================================================================

'use strict';

// ----- localStorage shim -----------------------------------------------------
const _ls = new Map();
globalThis.localStorage = {
  getItem(k) { return _ls.has(k) ? _ls.get(k) : null; },
  setItem(k, v) { _ls.set(k, String(v)); },
  removeItem(k) { _ls.delete(k); },
  clear() { _ls.clear(); },
};

// ----- DOM shim with Proxy-based style -------------------------------------
function makeEl(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    nodeType: 1,
    children: [],
    parentNode: null,
    style: null,
    dataset: {},
    attributes: {},
    listeners: new Map(),
    isContentEditable: false,
    _text: '',
    _value: '',
    offsetWidth: 100,
    offsetHeight: 30,
    appendChild(child) { child.parentNode = this; this.children.push(child); return child; },
    removeChild(child) {
      const i = this.children.indexOf(child);
      if (i >= 0) this.children.splice(i, 1);
      child.parentNode = null;
      return child;
    },
    addEventListener(type, fn) {
      if (!this.listeners.has(type)) this.listeners.set(type, []);
      this.listeners.get(type).push(fn);
    },
    removeEventListener(type, fn) {
      const list = this.listeners.get(type) || [];
      const i = list.indexOf(fn);
      if (i >= 0) list.splice(i, 1);
    },
    dispatchEvent(ev) {
      ev.target = ev.target || this;
      const list = this.listeners.get(ev.type) || [];
      for (const fn of list) fn(ev);
    },
    setAttribute(k, v) {
      this.attributes[k] = String(v);
      if (k.startsWith('data-')) this.dataset[k.slice(5)] = String(v);
    },
    getAttribute(k) { return this.attributes[k] !== undefined ? this.attributes[k] : null; },
    set className(v) { this.attributes.class = v; },
    get className()  { return this.attributes.class || ''; },
    set textContent(v) { this._text = String(v); this.children.length = 0; },
    get textContent()  {
      if (this._text) return this._text;
      return this.children.map(c => c.textContent || '').join('');
    },
    set id(v)        { this.attributes.id = v; },
    get id()         { return this.attributes.id || ''; },
    set value(v)     { this._value = v; },
    get value()      { return this._value; },
    set disabled(v)  { this.attributes.disabled = v ? 'true' : 'false'; },
    get disabled()   { return this.attributes.disabled === 'true'; },
    get firstChild() { return this.children[0]; },
    querySelector(sel) {
      const m = /^\[([\w-]+)="([^"]+)"\]$/.exec(sel);
      if (!m) return null;
      const [_, attr, val] = m;
      function walk(node) {
        if (node.attributes && node.attributes[attr] === val) return node;
        for (const c of node.children || []) { const f = walk(c); if (f) return f; }
        return null;
      }
      return walk(this);
    },
  };
  const props = {};
  el.style = new Proxy({ _props: props }, {
    get(_t, p) {
      if (p === '_props') return props;
      if (p === 'cssText') return Object.entries(props).map(([k,v]) => `${k}: ${v}`).join('; ');
      return props[String(p)] !== undefined ? props[String(p)] : '';
    },
    set(_t, p, v) {
      if (p === '_props') return true;
      if (p === 'cssText') {
        for (const k of Object.keys(props)) delete props[k];
        if (!v) return true;
        for (const part of String(v).split(';')) {
          const [k, val] = part.split(':');
          if (k && val) props[k.trim()] = val.trim();
        }
        return true;
      }
      props[String(p)] = String(v);
      return true;
    },
  });
  return el;
}

const _docKeyListeners = [];
const _alerts = [], _prompts = [];
let _promptResponses = [];

globalThis.document = {
  body: makeEl('body'),
  head: makeEl('head'),
  readyState: 'complete',
  activeElement: null,
  createElement: makeEl,
  createTextNode(t) { return { nodeType: 3, textContent: String(t) }; },
  getElementById(id) {
    function walk(node) {
      if (node.attributes && node.attributes.id === id) return node;
      for (const c of node.children || []) { const f = walk(c); if (f) return f; }
      return null;
    }
    return walk(document.head) || walk(document.body);
  },
  addEventListener(type, fn) {
    if (type === 'keydown') _docKeyListeners.push(fn);
  },
  removeEventListener() {},
};

globalThis.window = {
  innerWidth: 1280, innerHeight: 800,
  popgen: null,        // set below from real engine module
  popgenLive: null,    // set below
  popgenDock: null,    // set below
  popgenDockTabs: null,
  state: null,
};
globalThis.alert  = (msg) => { _alerts.push(msg); };
globalThis.prompt = (msg, def) => {
  _prompts.push(msg);
  if (_promptResponses.length > 0) return _promptResponses.shift();
  return def || '';
};

// ----- runner ---------------------------------------------------------------
let _pass = 0, _fail = 0;
function ok(label, cond) {
  if (cond) { console.log('  ok  ' + label); _pass++; }
  else      { console.log('FAIL ' + label); _fail++; }
}
function eq(label, a, b) {
  const A = JSON.stringify(a), B = JSON.stringify(b);
  if (A === B) { console.log('  ok  ' + label); _pass++; }
  else         { console.log('FAIL ' + label + ': ' + A + ' != ' + B); _fail++; }
}
function group(name, fn) { console.log(name); return fn(); }

// Async wrapper for sequencing — call each group's promise serially
async function runTests(testsFn) { await testsFn(); }

// ----- minimal atlas state mock for dim resolution --------------------------
function makeAtlasState() {
  const N = 12;   // bumped from 6 → 12 so each regime has 4 samples (passes min_n=3)
  const samples = [];
  for (let i = 0; i < N; i++) samples.push('CGA_' + String(i + 1).padStart(3,'0'));
  const fish_calls = [];
  for (let si = 0; si < N; si++) {
    let regime;
    if (si < 4) regime = 0;
    else if (si < 8) regime = 1;
    else regime = 2;
    fish_calls.push({
      sample_idx: si, sample_id: samples[si], regime,
      confidence: 0.95, n_supporting: 3, n_intervals: 3,
      ambiguous: false, votes: [regime,regime,regime],
    });
  }
  return {
    data: {
      chrom: 'C_gar_LG28', cohort_key: 'tab_test', n_samples: N,
      samples,
      windows: [{idx:0,center_mb:0.5,start_bp:0,end_bp:1e6,z:0.1}],
      candidates: [{
        id: 'cand_LG28_15Mb', start_mb: 15, end_mb: 18,
        start_bp: 15e6, end_bp: 18e6, k: 3, _system: 'detailed', fish_calls,
      }],
      family: [1,1,1,2,2,2,null,null,null,null,null,null],
      relatedness_hub: [null,null,null,null,30,null,null,null,null,null,null,null],
      f_roh: [0.01,0.02,0.05,0.03,0.08,0.12,0.04,0.05,0.02,0.01,0.03,0.04],
      f_roh_threshold: 0.10,
    },
    cur: 0, tracked: [0,4,8],
    candidate: null,
    candidates: null,
    activeSampleSet: null, k: 3, labelVocab: 'detailed',
  };
}
const atlasState = makeAtlasState();
atlasState.candidate = atlasState.data.candidates[0];
atlasState.candidates = atlasState.data.candidates;
globalThis.state = atlasState;
globalThis.window.state = atlasState;

// ----- import the engine + dock tabs ----------------------------------------
const engine = require('/home/claude/atlas_groups/atlas_group_engine.js');
engine.invalidateDimensions();
engine.loadFromLocalStorage();
globalThis.popgen = engine;
globalThis.window.popgen = engine;

// Stub the dock with bodies that setTabContent can write to
const dockBodies = { groups: makeEl('div'), lasso: makeEl('div'), snapshots: makeEl('div') };
const dockMock = {
  setTabContent(id, el) {
    const body = dockBodies[id];
    while (body.firstChild) body.removeChild(body.firstChild);
    body.appendChild(el);
  },
  getTabBody(id) { return dockBodies[id]; },
  onSaveSnapshot() {},
  onTabChange()    {},
};
globalThis.popgenDock = dockMock;
globalThis.window.popgenDock = dockMock;

// Stub live layer with deterministic responses
let _liveCalls = 0, _liveLastReq = null;
const liveMock = {
  popstatsGroupwise: async (req, opts) => {
    _liveCalls++; _liveLastReq = req;
    return {
      ok: true, kind: 'popstats_groupwise', cache_key: 'fakekey',
      cacheState: 'miss', payload: { windows: [], n_windows: 0 },
      request_ms: 42, fetched_ms: 30, ts: Date.now(),
    };
  },
};
globalThis.popgenLive = liveMock;
globalThis.window.popgenLive = liveMock;

// import tab module
const tabs = require('./atlas_dock_tabs.js');
globalThis.popgenDockTabs = tabs;
globalThis.window.popgenDockTabs = tabs;

// =============================================================================
// Tests
// =============================================================================

(async () => {

group('makeGroupsTab structure', () => {
  const t = tabs.makeGroupsTab();
  ok('returns element',   !!t);
  eq('data-popgen-tab',   t.attributes['data-popgen-tab'], 'groups');
  // Should contain four major sections: dim, draft, slots, compute strip
  ok('has 4+ children',   t.children.length >= 4);
});

group('makeLassoTab structure', () => {
  const t = tabs.makeLassoTab();
  ok('returns element',   !!t);
  eq('data-popgen-tab',   t.attributes['data-popgen-tab'], 'lasso');
});

group('makeSnapshotsTab structure', () => {
  const t = tabs.makeSnapshotsTab();
  ok('returns element',   !!t);
  eq('data-popgen-tab',   t.attributes['data-popgen-tab'], 'snapshots');
  // should have a "save snapshot" button at top
  ok('first child is button',  t.children[0].tagName === 'BUTTON');
});

group('install + dock integration', () => {
  tabs.install();
  ok('groups tab installed in dock',
    dockBodies.groups.children.length === 1);
  ok('lasso tab installed',
    dockBodies.lasso.children.length === 1);
  ok('snapshots tab installed',
    dockBodies.snapshots.children.length === 1);
});

group('draft slot — predicates round-trip', () => {
  tabs._internals._clearDraft();
  const draft = tabs.draftSlot();
  eq('starts empty',           draft.predicates, []);
  tabs._internals._addDraftPredicate({
    dim: 'diploid_class@cand_LG28_15Mb', op: 'eq', value: 'H1/H1',
  });
  tabs._internals._addDraftPredicate({
    dim: 'family', op: 'eq', value: 'family_1',
  });
  eq('two predicates', tabs.draftSlot().predicates.length, 2);
});

group('draft slot — converts to expression AST', () => {
  // From state above, two eq predicates joined by AND
  const expr = tabs._internals._draftToExpr();
  eq('top-level type',    expr.type, 'and');
  eq('two children',      expr.children.length, 2);
  eq('first eq dim',      expr.children[0].dim, 'diploid_class@cand_LG28_15Mb');
  eq('first eq value',    expr.children[0].value, 'H1/H1');
});

group('draft slot — single predicate is bare expr (no and wrapper)', () => {
  tabs._internals._clearDraft();
  tabs._internals._addDraftPredicate({
    dim: 'roh_carrier', op: 'truthy',
  });
  const expr = tabs._internals._draftToExpr();
  eq('bare type', expr.type, 'truthy');
  eq('bare dim',  expr.dim, 'roh_carrier');
});

await group('compute fires live request with right shape', async () => {
  // Set up a real slot so compute has something to request
  engine.setSlot(0, {
    name: 'HOM1', source: 'inline',
    ref: { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H1/H1' },
    color: '#1f4e79',
  });
  engine.setSlot(1, {
    name: 'HOM2', source: 'inline',
    ref: { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2' },
    color: '#3cc08a',
  });
  _liveCalls = 0; _liveLastReq = null;
  const env = await tabs.compute({ scope: 'candidate' });
  ok('one live call',          _liveCalls === 1);
  ok('request has chrom',      _liveLastReq && _liveLastReq.chrom === 'C_gar_LG28');
  ok('request has region',     _liveLastReq && _liveLastReq.region != null);
  eq('request region start',   _liveLastReq.region.start_bp, 15_000_000);
  eq('request region end',     _liveLastReq.region.end_bp,   18_000_000);
  eq('two groups',             Object.keys(_liveLastReq.groups).sort(), ['HOM1','HOM2']);
  ok('returns ok envelope',    env.ok === true);
  // Should have stashed the response into atlas state
  ok('stashed in state.popstatsLive',
    atlasState.popstatsLive && atlasState.popstatsLive.lastResponse != null);
});

await group('compute with no slots returns error envelope', async () => {
  engine.clearSlot(0); engine.clearSlot(1);
  // Clear all slots from earlier tests
  for (let i = 0; i < engine.MAX_SLOTS; i++) engine.clearSlot(i);
  _liveCalls = 0;
  const env = await tabs.compute({ scope: 'candidate' });
  ok('returns error',     env.ok === false);
  ok('error mentions groups', /groups/i.test(env.error));
  ok('no live call',      _liveCalls === 0);
});

group('lasso row save-as-set creates a saved set', () => {
  // Record a fake lasso so the row exists
  engine.recordSelection({
    sample_ids: ['CGA_001', 'CGA_002'],
    source_page: 'page1', kind: 'rect_lasso',
  });
  // Set up the prompt response BEFORE clicking
  _promptResponses = ['my_test_lasso'];
  // Refresh lasso tab so the new row appears
  tabs.refresh('lasso');
  // The lasso tab's first child is the list div; first row inside has our entry
  const lassoTab = tabs._internals.lassoTab;
  const list = lassoTab.children[0];
  ok('one lasso row',       list.children.length === 1);
  // Find the "save as set" button (text '⤴')
  const row = list.children[0];
  // Buttons in the row: pin, save-as-set, drop. Find by text.
  const saveBtn = row.children.find(c => c.tagName === 'BUTTON' &&
                                          c.textContent === '⤴');
  ok('save-as-set button found',  !!saveBtn);
  // Click → should trigger lassoToSet
  saveBtn.dispatchEvent({ type: 'click', target: saveBtn,
                          preventDefault() {}, stopPropagation() {} });
  ok('saved set was created',
    engine.listSets().includes('my_test_lasso'));
  eq('saved set members',
    engine.getSet('my_test_lasso'), ['CGA_001','CGA_002']);
});

group('snapshot save → restore round-trip', () => {
  // Set a slot, take snapshot, mutate, restore
  engine.setSlot(0, {
    name: 'TESTSLOT', source: 'inline',
    ref: { type: 'eq', dim: 'family', value: 'family_1' },
  });
  const snapId = engine.saveSnapshot({ note: 'testsnap' });
  ok('snapshot saved', !!snapId);
  // Mutate
  engine.clearSlot(0);
  ok('slot cleared', engine.getSlots()[0] == null);
  // Restore
  engine.restoreSnapshot(snapId);
  const slotsAfter = engine.getSlots();
  ok('slot restored',  slotsAfter[0] != null);
  eq('slot name kept', slotsAfter[0].name, 'TESTSLOT');
});

group('refresh propagates after restore', () => {
  // Refresh all tabs — should not throw
  let threw = false;
  try { tabs.refreshAll(); }
  catch (e) { threw = true; console.error(e); }
  ok('refreshAll runs without error', !threw);
});

group('scope chip persistence', () => {
  // Direct localStorage probe
  localStorage.setItem('inversion_atlas.popgen.scopeMode', '5w');
  // Re-create groups tab to pick up the new value
  const t = tabs.makeGroupsTab();
  // Find the scope chip (last section with class pg-scope-chip)
  function findByClass(node, cls) {
    if (node.attributes && (node.attributes.class || '').includes(cls)) return node;
    for (const c of node.children || []) {
      const f = findByClass(c, cls); if (f) return f;
    }
    return null;
  }
  const chip = findByClass(t, 'pg-scope-chip');
  ok('scope chip found',  !!chip);
  // The active button should be '5w'
  const active = (chip.children || []).find(c => c.attributes['data-on'] === 'true');
  ok('5w marked active',  active && active.textContent === '5w');
});

// ----- summary -------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);

})().catch(e => { console.error('UNCAUGHT', e); process.exit(1); });
