// =============================================================================
// test_atlas_track_gallery.js
// =============================================================================
// Run: node test_atlas_track_gallery.js
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

// ----- DOM shim (Proxy style for tag queries + attr setting) ---------------
function makeEl(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    nodeType: 1, children: [], parentNode: null,
    style: null, dataset: {}, attributes: {}, listeners: new Map(),
    isContentEditable: false, _text: '', _value: '',
    offsetWidth: 240, offsetHeight: 600,
    appendChild(c) { c.parentNode = this; this.children.push(c); return c; },
    removeChild(c) {
      const i = this.children.indexOf(c);
      if (i >= 0) this.children.splice(i, 1);
      c.parentNode = null; return c;
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
    set id(v) { this.attributes.id = v; },
    get id()  { return this.attributes.id || ''; },
    set value(v) { this._value = v; },
    get value()  { return this._value; },
    set placeholder(v) { this.attributes.placeholder = v; },
    set title(v)       { this.attributes.title = v; },
    set disabled(v) { this.attributes.disabled = v ? 'true' : 'false'; },
    get disabled()  { return this.attributes.disabled === 'true'; },
    get firstChild() { return this.children[0]; },
    querySelector(sel) {
      const m = /^\[([\w-]+)="([^"]+)"\]$/.exec(sel);
      if (!m) return null;
      const [_, attr, val] = m;
      function walk(n) {
        if (n.attributes && n.attributes[attr] === val) return n;
        for (const c of n.children || []) { const f = walk(c); if (f) return f; }
        return null;
      }
      return walk(this);
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

globalThis.document = {
  body: makeEl('body'), head: makeEl('head'), readyState: 'complete',
  createElement: makeEl,
  createTextNode(t) { return { nodeType: 3, textContent: String(t),
                                 parentNode: null, children: [], attributes: {} }; },
  getElementById(id) {
    function walk(n) {
      if (n.attributes && n.attributes.id === id) return n;
      for (const c of n.children || []) { const f = walk(c); if (f) return f; }
      return null;
    }
    return walk(document.head) || walk(document.body);
  },
  addEventListener() {}, removeEventListener() {},
};
globalThis.window = {
  innerWidth: 1280, innerHeight: 800,
  addEventListener() {}, removeEventListener() {},
  dispatchEvent() {},
  CustomEvent: function(type, init) { return { type, detail: (init && init.detail) || {} }; },
};
const _alerts = [], _prompts = [];
let _promptResponses = [];
globalThis.alert = (m) => { _alerts.push(m); };
globalThis.prompt = (m, def) => {
  _prompts.push(m);
  if (_promptResponses.length > 0) return _promptResponses.shift();
  return def || '';
};
globalThis.CustomEvent = window.CustomEvent;

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

// ----- atlas state mock -----------------------------------------------------
function makePayload() {
  return {
    columns: ['window_id', 'chrom', 'start', 'end', 'center_mb',
              'n_sites', 'n_sites_used', 'S',
              'theta_pi', 'theta_w', 'tajD', 'het',
              'theta_pi_HOM1', 'theta_pi_HET', 'theta_pi_HOM2',
              'Fst_HOM1_HOM2', 'Fst_HOM1_HET', 'Fst_HET_HOM2',
              'dXY_HOM1_HOM2', 'dA_HOM1_HOM2',
              'MI_HOM1_HOM2', 'MInorm_HOM1_HOM2'],
    windows: [
      { window_id: 'w0', chrom: 'C_gar_LG28', start: 0, end: 1e6,
        center_mb: 0.5, n_sites: 100, n_sites_used: 90, S: 12,
        theta_pi: 0.005, theta_w: 0.006, tajD: -0.3, het: 0.1,
        theta_pi_HOM1: 0.004, theta_pi_HET: 0.008, theta_pi_HOM2: 0.003,
        Fst_HOM1_HOM2: 0.10, Fst_HOM1_HET: 0.05, Fst_HET_HOM2: 0.06,
        dXY_HOM1_HOM2: 0.012, dA_HOM1_HOM2: 0.008,
        MI_HOM1_HOM2: 0.4, MInorm_HOM1_HOM2: 0.7 },
    ],
  };
}
const atlasState = {
  data: { chrom: 'C_gar_LG28', cohort_key: 'gallery_test',
          n_samples: 12, samples: [], windows: [],
          tracks: { snp_density: { values: [50, 60, 70] } } },
  popstatsLive: { lastResponse: makePayload() },
  cur: 0, tracked: [], k: 3,
};
globalThis.state = atlasState;
globalThis.window.state = atlasState;

// engine + renderers (small stubs for cohort key + adapt)
globalThis.popgen = {
  getCohortKey: () => 'gallery_test',
};
globalThis.window.popgen = globalThis.popgen;
globalThis.popgenRenderers = require('/home/claude/atlas_renderers/atlas_renderers_turn4.js');
globalThis.window.popgenRenderers = globalThis.popgenRenderers;

const gallery = require('./atlas_track_gallery.js');

// =============================================================================
// Tests
// =============================================================================

group('column → track-def derivation', () => {
  const inn = gallery._internals;
  // Cohort θπ
  let d = inn._trackDefForColumn('theta_pi', null);
  ok('theta_pi → cohort track', d && d.id === 'theta_pi_cohort');
  eq('theta_pi family',  d.family, 'cohort');
  // Tajima D
  d = inn._trackDefForColumn('tajD', null);
  ok('tajD → tajima_d',  d && d.id === 'tajima_d');
  // Per-group θπ
  d = inn._trackDefForColumn('theta_pi_HOM1', null);
  ok('theta_pi_HOM1 → theta_pi_HOM1', d && d.id === 'theta_pi_HOM1');
  eq('per-group family',           d.family, 'theta_pi');
  ok('per-group has color',         typeof d.color === 'string');
  // Pairwise Fst
  d = inn._trackDefForColumn('Fst_HOM1_HOM2', null);
  ok('Fst_HOM1_HOM2 → fst_hom1_hom2', d && d.id === 'fst_hom1_hom2');
  eq('label format',  d.label, 'Fst HOM1–HOM2');
  // dXY
  d = inn._trackDefForColumn('dXY_HOM1_HOM2', null);
  ok('dXY_HOM1_HOM2 → dxy_*',  d && d.id === 'dxy_hom1_hom2');
  // Skipped columns
  d = inn._trackDefForColumn('window_id', null);
  ok('window_id skipped',  d == null);
  d = inn._trackDefForColumn('start', null);
  ok('start skipped',      d == null);
  d = inn._trackDefForColumn('n_sites', null);
  ok('n_sites skipped',    d == null);
});

group('discoverTracks finds canonical + per-column tracks', () => {
  const tracks = gallery.discoverTracks(makePayload());
  // Always-on canonical
  ok('ideogram present',         tracks.find(t => t.id === 'ideogram') != null);
  ok('sim_collapse present',     tracks.find(t => t.id === 'sim_collapse') != null);
  ok('z present',                tracks.find(t => t.id === 'z') != null);
  // Q04 canonical placeholders
  ok('theta_invgt canonical',    tracks.find(t => t.id === 'theta_invgt') != null);
  ok('fst_hom1_hom2 canonical',  tracks.find(t => t.id === 'fst_hom1_hom2') != null);
  // Per-group θπ from columns
  ok('theta_pi_HOM1 derived',    tracks.find(t => t.id === 'theta_pi_HOM1') != null);
  ok('theta_pi_HET derived',     tracks.find(t => t.id === 'theta_pi_HET') != null);
  ok('theta_pi_HOM2 derived',    tracks.find(t => t.id === 'theta_pi_HOM2') != null);
  // Pairwise Fst (3 pairs)
  ok('fst_hom1_het derived',     tracks.find(t => t.id === 'fst_hom1_het') != null);
  ok('fst_het_hom2 derived',     tracks.find(t => t.id === 'fst_het_hom2') != null);
  // dXY / dA / MI / MInorm
  ok('dxy_hom1_hom2 derived',    tracks.find(t => t.id === 'dxy_hom1_hom2') != null);
  ok('da_hom1_hom2 derived',     tracks.find(t => t.id === 'da_hom1_hom2') != null);
  ok('mi_hom1_hom2 derived',     tracks.find(t => t.id === 'mi_hom1_hom2') != null);
  ok('minorm_hom1_hom2 derived', tracks.find(t => t.id === 'minorm_hom1_hom2') != null);
  // Cohort θπ
  ok('theta_pi_cohort derived',  tracks.find(t => t.id === 'theta_pi_cohort') != null);
  // Tajima
  ok('tajima_d derived',         tracks.find(t => t.id === 'tajima_d') != null);
});

group('discoverTracks works with no payload', () => {
  const tracks = gallery.discoverTracks(null);
  // Only canonical defs
  ok('still has ideogram',  tracks.find(t => t.id === 'ideogram') != null);
  ok('still has theta_invgt', tracks.find(t => t.id === 'theta_invgt') != null);
  // No per-group derived
  ok('no theta_pi_HOM1',    tracks.find(t => t.id === 'theta_pi_HOM1') == null);
});

group('pattern matcher', () => {
  const m = gallery._internals._matchPattern;
  ok('exact match',         m('theta_pi_HOM1', 'theta_pi_HOM1'));
  ok('prefix glob match',   m('theta_pi_HOM1', 'theta_pi_*'));
  ok('prefix glob negative', !m('Fst_HOM1_HOM2', 'theta_pi_*'));
  ok('star wildcard',       m('anything', '*'));
  ok('null pattern',        !m('x', null));
});

group('preset resolution: q04_stack returns Q04_STACK', () => {
  const tracks = gallery.discoverTracks(makePayload());
  const ids = gallery._internals._resolvePreset('q04_stack', tracks);
  eq('exact Q04_STACK', ids, gallery.Q04_STACK);
});

group('preset resolution: theta_deep matches per-group θπ + cohort + tajD', () => {
  const tracks = gallery.discoverTracks(makePayload());
  const ids = gallery._internals._resolvePreset('theta_deep', tracks);
  ok('contains theta_pi_HOM1',    ids.includes('theta_pi_HOM1'));
  ok('contains theta_pi_HOM2',    ids.includes('theta_pi_HOM2'));
  ok('contains theta_pi_cohort',  ids.includes('theta_pi_cohort'));
  ok('contains tajima_d',         ids.includes('tajima_d'));
  ok('does NOT contain Fst',     !ids.find(x => x.startsWith('fst_')));
  ok('ideogram pinned at top',    ids[0] === 'ideogram');
  ok('sim_collapse 2nd',          ids[1] === 'sim_collapse');
});

group('preset resolution: fst_all_pairs', () => {
  const tracks = gallery.discoverTracks(makePayload());
  const ids = gallery._internals._resolvePreset('fst_all_pairs', tracks);
  ok('fst_hom1_hom2',  ids.includes('fst_hom1_hom2'));
  ok('fst_hom1_het',   ids.includes('fst_hom1_het'));
  ok('fst_het_hom2',   ids.includes('fst_het_hom2'));
  ok('dxy_hom1_hom2',  ids.includes('dxy_hom1_hom2'));
});

group('active track ids: default = Q04_STACK', () => {
  const ids = gallery.getActiveTrackIds();
  eq('default = Q04_STACK', ids, gallery.Q04_STACK);
});

group('active track ids: toggle', () => {
  // Start fresh
  gallery.setActiveTrackIds([]);
  gallery.toggleTrack('theta_pi_HOM1');
  ok('added',  gallery.isTrackActive('theta_pi_HOM1'));
  gallery.toggleTrack('theta_pi_HOM1');
  ok('removed', !gallery.isTrackActive('theta_pi_HOM1'));
});

group('applyPreset switches active set', () => {
  gallery.applyPreset('fst_all_pairs');
  const ids = gallery.getActiveTrackIds();
  ok('fst_hom1_hom2 active',  ids.includes('fst_hom1_hom2'));
  ok('theta_invgt NOT active', !ids.includes('theta_invgt'));
});

group('user presets — save / list / delete', () => {
  // Save the current active set
  gallery.applyPreset('q04_stack');
  gallery.savePreset('my_q04_lite', ['ideogram', 'theta_invgt', 'fst_hom1_hom2']);
  const presets = gallery.listPresets();
  ok('listed',  presets.find(p => p.name === 'my_q04_lite') != null);
  // Apply user preset
  gallery.applyPreset('my_q04_lite');
  const ids = gallery.getActiveTrackIds();
  eq('user preset applied', ids.sort(),
    ['fst_hom1_hom2', 'ideogram', 'theta_invgt']);
  // Delete
  gallery.deletePreset('my_q04_lite');
  const after = gallery.listPresets();
  ok('removed from list', after.find(p => p.name === 'my_q04_lite') == null);
});

group('user preset — name validation', () => {
  let threw = false;
  try { gallery.savePreset('bad name', ['ideogram']); }
  catch (_) { threw = true; }
  ok('rejects spaces', threw);
  threw = false;
  try { gallery.savePreset('q04_stack', ['ideogram']); }
  catch (_) { threw = true; }
  ok('rejects builtin name', threw);
});

group('search filter', () => {
  const tracks = gallery.discoverTracks(makePayload());
  const inn = gallery._internals;
  const matched = tracks.filter(t => inn._matchSearch(t, 'fst'));
  ok('Fst filter', matched.length >= 3);
  const no = tracks.filter(t => inn._matchSearch(t, 'nonexistent_xyz'));
  eq('no matches', no.length, 0);
});

group('_groupTracksByFamily: Q04 composites + per-group families', () => {
  const tracks = gallery.discoverTracks(makePayload());
  const grouped = gallery._internals._groupTracksByFamily(tracks);
  // First group should be 'composite'
  eq('first family', grouped[0].family.key, 'composite');
  ok('composite has ideogram',
    grouped[0].tracks.find(t => t.id === 'ideogram') != null);
  // Should have a theta_pi family with at least HOM1/HET/HOM2
  const thetaFam = grouped.find(g => g.family.key === 'theta_pi');
  ok('theta_pi family present',  !!thetaFam);
  ok('contains HOM1 in theta family',
    thetaFam.tracks.find(t => t.id === 'theta_pi_HOM1') != null);
  // Should have an fst family
  const fstFam = grouped.find(g => g.family.key === 'fst');
  ok('fst family present',  !!fstFam);
  ok('fst family has 3+ pairs', fstFam.tracks.length >= 3);
});

group('makeGallery DOM construction', () => {
  const root = gallery.makeGallery({ onTrackToggle: () => {} });
  ok('root element',  !!root);
  eq('data-popgen-gallery=root',
    root.attributes['data-popgen-gallery'], 'root');
  // Header has title + preset select + close (no close because no onClose)
  const header = root.children.find(c =>
    c.attributes && c.attributes['data-popgen-gallery'] === 'header');
  ok('header present',   !!header);
  // Search input
  const search = root.children.find(c =>
    c.attributes && c.attributes['data-popgen-gallery'] === 'search');
  ok('search input present',  !!search);
  // Body
  const body = root.children.find(c =>
    c.attributes && c.attributes['data-popgen-gallery'] === 'body');
  ok('body present',     !!body);
  // Body has at least one category div with chips
  ok('body has cat',     body.children.length >= 1);
  // Active chips count: each Q04_STACK id present in discovered tracks
  // should render with data-on=true
  function findChips(node, found) {
    if (node.attributes && node.attributes['data-popgen-gallery'] === 'chip') {
      found.push(node);
    }
    for (const c of node.children || []) findChips(c, found);
  }
  const chips = []; findChips(root, chips);
  ok('chips rendered',  chips.length >= 5);
  const onChips = chips.filter(c => c.attributes['data-on'] === 'true');
  ok('some chips active by default',  onChips.length >= 1);
});

group('makeGallery toggle chip flips active state', () => {
  // Reset to Q04 stack
  gallery.applyPreset('q04_stack');
  const root = gallery.makeGallery({ onTrackToggle: () => {} });
  // Find a chip for a non-active track (theta_pi_HOM1)
  function findChips(node, found) {
    if (node.attributes && node.attributes['data-popgen-gallery'] === 'chip') {
      found.push(node);
    }
    for (const c of node.children || []) findChips(c, found);
  }
  const chips = []; findChips(root, chips);
  const targetChip = chips.find(c => c.attributes.title === 'theta_pi_HOM1');
  ok('target chip found',  !!targetChip);
  eq('initially off',  targetChip.attributes['data-on'], 'false');
  // Click
  targetChip.dispatchEvent({ type: 'click', target: targetChip,
                              preventDefault() {}, stopPropagation() {} });
  eq('now on',  targetChip.attributes['data-on'], 'true');
  ok('engine state updated',
    gallery.getActiveTrackIds().includes('theta_pi_HOM1'));
});

// ----- summary -------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);
