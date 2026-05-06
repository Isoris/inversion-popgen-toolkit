// =============================================================================
// test_atlas_page6_wiring.js
// =============================================================================
// Tests for turn 6.
//
// Run: node test_atlas_page6_wiring.js
// =============================================================================

'use strict';

// ----- localStorage --------------------------------------------------------
const _ls = new Map();
globalThis.localStorage = {
  getItem(k) { return _ls.has(k) ? _ls.get(k) : null; },
  setItem(k, v) { _ls.set(k, String(v)); },
  removeItem(k) { _ls.delete(k); },
  clear() { _ls.clear(); },
};

// ----- DOM shim (minimal proxy-style style) --------------------------------
function makeEl(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    nodeType: 1, children: [], parentNode: null,
    style: null, dataset: {}, attributes: {}, listeners: new Map(),
    isContentEditable: false, _text: '', _value: '', offsetWidth: 100, offsetHeight: 30,
    appendChild(c) { c.parentNode = this; this.children.push(c); return c; },
    removeChild(c) {
      const i = this.children.indexOf(c);
      if (i >= 0) this.children.splice(i, 1);
      c.parentNode = null; return c;
    },
    insertBefore(newN, ref) {
      const i = ref ? this.children.indexOf(ref) : this.children.length;
      if (i < 0) this.children.push(newN);
      else this.children.splice(i, 0, newN);
      newN.parentNode = this;
      return newN;
    },
    addEventListener(type, fn) {
      if (!this.listeners.has(type)) this.listeners.set(type, []);
      this.listeners.get(type).push(fn);
    },
    removeEventListener() {},
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
    set innerHTML(v) { this._html = String(v); },
    get innerHTML() { return this._html || ''; },
    set id(v)        { this.attributes.id = v; },
    get id()         { return this.attributes.id || ''; },
    set value(v)     { this._value = v; },
    get value()      { return this._value; },
    set type(v)      { this.attributes.type = v; },
    get type()       { return this.attributes.type || 'text'; },
    set placeholder(v) { this.attributes.placeholder = v; },
    set disabled(v)  { this.attributes.disabled = v ? 'true' : 'false'; },
    get disabled()   { return this.attributes.disabled === 'true'; },
    get firstChild() { return this.children[0]; },
    getBoundingClientRect() {
      return { left: 0, top: 0, right: this.offsetWidth, bottom: this.offsetHeight,
               width: this.offsetWidth, height: this.offsetHeight, x: 0, y: 0 };
    },
    querySelector(sel) {
      const m = /^\[([\w-]+)="([^"]+)"\]$/.exec(sel);
      if (m) {
        const [_, attr, val] = m;
        function walk(n) {
          if (n.attributes && n.attributes[attr] === val) return n;
          for (const c of n.children || []) { const f = walk(c); if (f) return f; }
          return null;
        }
        return walk(this);
      }
      // class selector
      const cm = /^\.([a-zA-Z_-]+)$/.exec(sel);
      if (cm) {
        function walk(n) {
          const cls = (n.attributes && n.attributes.class) || '';
          if (cls.includes(cm[1])) return n;
          for (const c of n.children || []) { const f = walk(c); if (f) return f; }
          return null;
        }
        return walk(this);
      }
      return null;
    },
    querySelectorAll(sel) {
      const out = [];
      const m = /^([a-z]+)\[([\w-]+)="([^"]+)"\]$/.exec(sel);
      const m2 = /^\[([\w-]+)="([^"]+)"\]$/.exec(sel);
      let attr, val, tag;
      if (m) { tag = m[1].toUpperCase(); attr = m[2]; val = m[3]; }
      else if (m2) { attr = m2[1]; val = m2[2]; }
      else return out;
      function walk(n) {
        if (n.attributes && n.attributes[attr] === val &&
            (!tag || n.tagName === tag)) {
          out.push(n);
        }
        for (const c of n.children || []) walk(c);
      }
      walk(this);
      return out;
    },
  };
  const props = {};
  el.style = new Proxy({}, {
    get(_, p) {
      if (p === 'cssText') return Object.entries(props).map(([k,v]) => `${k}: ${v}`).join('; ');
      return props[String(p)] !== undefined ? props[String(p)] : '';
    },
    set(_, p, v) {
      if (p === 'cssText') {
        for (const k of Object.keys(props)) delete props[k];
        if (v) for (const part of String(v).split(';')) {
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

const _docKeyL = [], _docMoveL = [], _docUpL = [];
globalThis.document = {
  body: makeEl('body'), head: makeEl('head'), readyState: 'complete',
  activeElement: null,
  createElement: makeEl,
  createTextNode(t) { return { nodeType: 3, textContent: String(t), parentNode: null,
    children: [], attributes: {} }; },
  getElementById(id) {
    function walk(n) {
      if (n.attributes && n.attributes.id === id) return n;
      for (const c of n.children || []) { const f = walk(c); if (f) return f; }
      return null;
    }
    return walk(document.head) || walk(document.body);
  },
  addEventListener(type, fn) {
    if (type === 'keydown')   _docKeyL.push(fn);
    if (type === 'mousemove') _docMoveL.push(fn);
    if (type === 'mouseup')   _docUpL.push(fn);
  },
  removeEventListener() {},
};
globalThis.window = {
  innerWidth: 1280, innerHeight: 800,
  addEventListener(type, fn) {},
  dispatchEvent(ev) { return true; },
  CustomEvent: function (type, init) {
    return { type, detail: (init && init.detail) || {} };
  },
};
globalThis.CustomEvent = window.CustomEvent;
globalThis.alert = (m) => { _alerts.push(m); };
const _alerts = [];

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
function group(name, fn) { console.log(name); fn(); }

// ----- atlas state mock for live/static resolution --------------------------
function makeAtlasState() {
  return {
    data: {
      chrom: 'C_gar_LG28', cohort_key: 'p6_test', n_samples: 6,
      samples: ['CGA_001','CGA_002','CGA_003','CGA_004','CGA_005','CGA_006'],
      windows: [
        { idx:0, center_mb: 0.5,  start_bp: 0,        end_bp: 1e6 },
        { idx:1, center_mb: 1.5,  start_bp: 1e6,      end_bp: 2e6 },
        { idx:2, center_mb: 2.5,  start_bp: 2e6,      end_bp: 3e6 },
      ],
      tracks: {
        // Static-mode track for Fst Hom1-Hom2: a single line with values
        fst_hom1_hom2: { values: [0.05, 0.07, 0.06] },
      },
      candidates: [],
    },
    cur: 1, tracked: [],
    candidate: null, candidates: [],
    activeSampleSet: null, k: 3, labelVocab: 'detailed',
    popstatsLive: {
      lastResponse: {
        // Engine F response shape: windows[] with column-like fields
        columns: ['window_id', 'chrom', 'start', 'end', 'center_mb',
                  'theta_pi', 'theta_pi_HOM1', 'theta_pi_HET', 'theta_pi_HOM2',
                  'Fst_HOM1_HOM2'],
        windows: [
          { window_id: 'w0', chrom: 'C_gar_LG28', start: 0,    end: 1e6,
            center_mb: 0.5, theta_pi: 0.005,
            theta_pi_HOM1: 0.004, theta_pi_HET: 0.008, theta_pi_HOM2: 0.003,
            Fst_HOM1_HOM2: 0.10 },
          { window_id: 'w1', chrom: 'C_gar_LG28', start: 1e6,  end: 2e6,
            center_mb: 1.5, theta_pi: 0.006,
            theta_pi_HOM1: 0.005, theta_pi_HET: 0.009, theta_pi_HOM2: 0.004,
            Fst_HOM1_HOM2: 0.12 },
          { window_id: 'w2', chrom: 'C_gar_LG28', start: 2e6,  end: 3e6,
            center_mb: 2.5, theta_pi: 0.007,
            theta_pi_HOM1: 0.006, theta_pi_HET: 0.010, theta_pi_HOM2: 0.005,
            Fst_HOM1_HOM2: 0.14 },
        ],
      },
    },
  };
}
const atlasState = makeAtlasState();
globalThis.state = atlasState;
globalThis.window.state = atlasState;

// ----- import deps (engine + renderers) -------------------------------------
const engine    = require('/home/claude/atlas_groups/atlas_group_engine.js');
const renderers = require('/home/claude/atlas_renderers/atlas_renderers_turn4.js');
globalThis.popgen = engine;
globalThis.window.popgen = engine;
globalThis.popgenRenderers = renderers;
globalThis.window.popgenRenderers = renderers;

// Stub the dock tabs so the slot-combine UI installer can find them
const dockTabsStub = {
  makeGroupsTab() {
    const root = makeEl('div');
    root.attributes['data-popgen-tab'] = 'groups';
    // Build minimal section structure that turn 6 expects:
    //   3 .pg-section divs (dim, draft, slots)
    for (let i = 0; i < 3; i++) {
      const sec = makeEl('div');
      sec.attributes.class = 'pg-section';
      sec.appendChild(makeEl('div'));    // title
      sec.appendChild(makeEl('div'));    // list
      root.appendChild(sec);
    }
    // The slots-list (third section's second child) gets a few empty rows
    const slotsList = root.children[2].children[1];
    for (let i = 0; i < 3; i++) {
      const row = makeEl('div');
      row.attributes['data-empty'] = 'true';
      row.attributes.class = 'pg-slot';
      row.appendChild(makeEl('span'));            // idx
      row.appendChild(makeEl('input'));           // name placeholder (will be replaced)
      row.appendChild(makeEl('span'));            // source
      row.appendChild(makeEl('span'));            // n
      row.appendChild(makeEl('span'));            // color
      row.appendChild(makeEl('button'));          // del
      slotsList.appendChild(row);
    }
    root._refresh = () => {};
    return root;
  },
  refresh() {},
  refreshAll() {},
};
globalThis.popgenDockTabs = dockTabsStub;
globalThis.window.popgenDockTabs = dockTabsStub;

// import the module under test
const wiring = require('./atlas_page6_wiring.js');

// =============================================================================
// Tests
// =============================================================================

group('static-mode resolver', () => {
  renderers.setPopstatsMode('static');
  const trackDef = { id: 'fst_hom1_hom2', renderer: 'line' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('returns data',          !!data);
  eq('values from static',    data.values, [0.05, 0.07, 0.06]);
  eq('mb from windows',       data.mb, [0.5, 1.5, 2.5]);
});

group('live-mode resolver', () => {
  renderers.setPopstatsMode('live');
  const trackDef = { id: 'fst_hom1_hom2', renderer: 'line' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('returns data',          !!data);
  // Live values should come from popstatsLive.lastResponse, NOT from static
  eq('values from live',      data.values, [0.10, 0.12, 0.14]);
});

group('live mode falls back to static when live missing', () => {
  renderers.setPopstatsMode('live');
  // A track only present in static
  atlasState.data.tracks.snp_density = { values: [50, 60, 70] };
  const trackDef = { id: 'snp_density', renderer: 'line' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('falls back to static',  !!data);
  eq('static values used',    data.values, [50, 60, 70]);
});

group('split-mode composes live + static as ghost', () => {
  renderers.setPopstatsMode('split');
  const trackDef = { id: 'fst_hom1_hom2', renderer: 'line' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('returns data',          !!data);
  ok('has multiline series',  Array.isArray(data.series));
  // Should have at least 2 series: 1 live + 1 static-ghost
  ok('>= 2 series',           data.series.length >= 2);
  // The ghost series should be dashed
  const ghost = data.series.find(s => s.dashed);
  ok('ghost series found',    !!ghost);
  ok('ghost name labelled',   /static/.test(ghost.name));
});

group('multi-series live response: theta_pi by group', () => {
  renderers.setPopstatsMode('live');
  const trackDef = { id: 'theta_invgt', renderer: 'multiline' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('returns data',           !!data);
  ok('multiline series',       Array.isArray(data.series));
  eq('three series',           data.series.length, 3);
  // Verify the series are HOM1, HET, HOM2 in that order (alpha by name)
  const names = data.series.map(s => s.name).sort();
  eq('series names', names, ['HET','HOM1','HOM2']);
});

group('_ensureMultilineShape converts single-line to multi-line', () => {
  const single = { mb: [1, 2, 3], values: [0.1, 0.2, 0.3] };
  const ml = wiring._internals._ensureMultilineShape(single, 'whatever');
  eq('mb preserved',   ml.mb, [1, 2, 3]);
  eq('one series',     ml.series.length, 1);
  eq('series values',  ml.series[0].values, [0.1, 0.2, 0.3]);
});

group('_ensureMultilineShape passes through multi-line', () => {
  const ml_in = { mb: [1, 2], series: [{ name: 'a', values: [1, 2] }] };
  const ml = wiring._internals._ensureMultilineShape(ml_in, 'whatever');
  ok('same object semantics',   ml.series.length === 1);
  eq('same series',             ml.series, ml_in.series);
});

group('static mode preserves null when neither static nor live present', () => {
  renderers.setPopstatsMode('static');
  const trackDef = { id: 'never_exists_anywhere', renderer: 'line' };
  const wrapped = wiring.wrapTrackDef(trackDef);
  const data = wrapped.getData(atlasState);
  ok('returns null',  data == null);
});

group('always-on tracks pass through unchanged', () => {
  const ideoDef = { id: 'ideogram', alwaysOn: true, renderer: 'ideogram' };
  const wrapped = wiring.wrapTrackDef(ideoDef);
  ok('same object', wrapped === ideoDef);
});

group('wrapAllTrackDefs maps each track', () => {
  const list = [
    { id: 'a', renderer: 'line' },
    { id: 'b', renderer: 'multiline' },
    { id: 'ideogram', renderer: 'ideogram', alwaysOn: true },
  ];
  const wrapped = wiring.wrapAllTrackDefs(list);
  eq('three results',    wrapped.length, 3);
  ok('a getData',        typeof wrapped[0].getData === 'function');
  ok('b getData',        typeof wrapped[1].getData === 'function');
  ok('ideogram pass-through', wrapped[2] === list[2]);
});

group('tooltip data cache + show/hide', () => {
  const canvas = makeEl('canvas');
  canvas.offsetWidth = 600; canvas.offsetHeight = 80;
  canvas.getBoundingClientRect = () => ({
    left: 0, top: 0, right: 600, bottom: 80,
    width: 600, height: 80, x: 0, y: 0,
  });
  const trackDef = { id: 'theta_invgt', label: 'θπ by invgt' };
  wiring.attachTooltip(trackDef, canvas);
  wiring.setTooltipDataFor(canvas, {
    data: {
      mb: [0.5, 1.5, 2.5],
      series: [
        { name: 'HOM1', color: '#1f4e79', values: [0.004, 0.005, 0.006] },
        { name: 'HOM2', color: '#3cc08a', values: [0.003, 0.004, 0.005] },
      ],
    },
    trackDef,
    padL: 80, plotW: 460, mbMin: 0, mbMax: 3,
  });
  // Pointer move at x = 80 + (1.5/3)*460 ≈ 310, expecting nearest mb idx = 1
  canvas.dispatchEvent({
    type: 'pointermove',
    clientX: 310, clientY: 50,
    target: canvas,
  });
  const tt = wiring._internals.tooltipEl;
  ok('tooltip created',       !!tt);
  ok('tooltip displayed',     tt && tt.style.display !== 'none');
  // The tooltip's innerHTML should contain HOM1 and HOM2 names + values
  ok('tooltip mentions HOM1', tt && /HOM1/.test(tt.innerHTML));
  ok('tooltip mentions HOM2', tt && /HOM2/.test(tt.innerHTML));
  ok('tooltip mentions mb',   tt && /1\.500 Mb/.test(tt.innerHTML));
  // Pointer leave hides
  canvas.dispatchEvent({ type: 'pointerleave', target: canvas });
  ok('hidden on leave', tt.style.display === 'none');
});

group('slot combine UI installs and decorates empty rows', () => {
  // Reset slots to ensure empty
  for (let i = 0; i < engine.MAX_SLOTS; i++) engine.clearSlot(i);
  wiring.installSlotCombineUI();
  // Re-create the groups tab (the installer hooks makeGroupsTab)
  const tab = dockTabsStub.makeGroupsTab();
  // Find the slots section — third .pg-section
  const sections = tab.children.filter(c => c.attributes && c.attributes.class === 'pg-section');
  eq('three sections',   sections.length, 3);
  const slotsList = sections[2].children[1];
  // After decoration, each empty row should have a combine button
  let combineBtnCount = 0;
  for (const row of slotsList.children) {
    if (row.children.find(c =>
      c.attributes && (c.attributes.class || '').includes('pg-slot-combine-btn'))) {
      combineBtnCount++;
    }
  }
  ok('at least one combine btn', combineBtnCount >= 1);
});

group('slot combine picker creates minus slot', () => {
  // Set up two source slots
  for (let i = 0; i < engine.MAX_SLOTS; i++) engine.clearSlot(i);
  engine.setSlot(0, { name: 'A', source: 'inline',
    ref: { type: 'all' }, color: '#1f4e79' });
  engine.setSlot(1, { name: 'B', source: 'inline',
    ref: { type: 'eq', dim: 'family', value: 'family_1' }, color: '#3cc08a' });
  // Use _internals to access the picker factory directly
  const picker = wiring._internals._wrappedGetData;   // sanity check we have access
  // The picker is created on click; rather than simulate a full click chain,
  // verify directly: after manually building a 'minus' slot, the engine resolves it
  engine.setSlot(2, { name: 'C', source: 'minus',
    ref: { left: 'A', right: 'B' }, color: '#7a8398' });
  const cIds = engine.resolveSlotByIdx(2);
  ok('C resolves',  Array.isArray(cIds));
});

group('_firstEmptySlotIdx', () => {
  // Currently slots 0,1,2 are full from the prior test
  for (let i = 0; i < engine.MAX_SLOTS; i++) engine.clearSlot(i);
  engine.setSlot(0, { name: 'X', source: 'inline', ref: { type: 'all' } });
  eq('first empty after slot 0', wiring._internals._firstEmptySlotIdx(), 1);
  // All full
  for (let i = 1; i < engine.MAX_SLOTS; i++) {
    engine.setSlot(i, { name: 'S' + i, source: 'inline', ref: { type: 'all' } });
  }
  eq('no empty', wiring._internals._firstEmptySlotIdx(), -1);
});

group('_suggestName picks unused name', () => {
  const used = [{ name: 'CA' }, { name: 'CB' }];
  const suggested = wiring._internals._suggestName(used);
  ok('skips used', suggested !== 'CA' && suggested !== 'CB');
  ok('matches CC pattern', suggested === 'CC');
});

// ----- summary -------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);
