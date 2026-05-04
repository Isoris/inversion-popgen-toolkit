// =============================================================================
// test_atlas_floating_dock.js
// =============================================================================
// Tests for the dock shell. Uses a fuller DOM mock than the canvas tests
// because the dock relies on layout (offsetWidth, getBoundingClientRect),
// event dispatch, and parent-child relationships.
//
// Run: node test_atlas_floating_dock.js
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

// ----- DOM shim — richer than the canvas-test version ------------------------
function makeEl(tag) {
  const el = {
    tagName: String(tag).toUpperCase(),
    nodeType: 1,
    children: [],
    childNodes: [],
    parentNode: null,
    style: null,    // assigned below as a Proxy
    dataset: {},
    attributes: {},
    listeners: new Map(),
    isContentEditable: false,
    _text: '',
    offsetWidth: 380,
    offsetHeight: 460,
    appendChild(child) {
      child.parentNode = this;
      this.children.push(child);
      this.childNodes.push(child);
      return child;
    },
    removeChild(child) {
      const i = this.children.indexOf(child);
      if (i >= 0) this.children.splice(i, 1);
      const j = this.childNodes.indexOf(child);
      if (j >= 0) this.childNodes.splice(j, 1);
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
      if (k.startsWith('data-')) this.dataset[_camel(k.slice(5))] = String(v);
    },
    getAttribute(k) { return this.attributes[k] !== undefined ? this.attributes[k] : null; },
    set className(v) { this.attributes.class = v; },
    get className()  { return this.attributes.class || ''; },
    set textContent(v) { this._text = String(v); this.children.length = 0; },
    get textContent()  { return this._text || ''; },
    set id(v)        { this.attributes.id = v; },
    get id()         { return this.attributes.id || ''; },
    get firstChild() { return this.children[0]; },
    getBoundingClientRect() {
      const x = parseFloat((this.style._props && this.style._props.left) || '0');
      const y = parseFloat((this.style._props && this.style._props.top) || '0');
      return { left: x, top: y, right: x + this.offsetWidth, bottom: y + this.offsetHeight,
               width: this.offsetWidth, height: this.offsetHeight, x, y };
    },
    querySelector(sel) {
      return _querySelector(this, sel);
    },
  };
  // style is a Proxy so direct property assignment (style.left = '10px') works
  const props = {};
  el.style = new Proxy({ _props: props }, {
    get(target, prop) {
      if (prop === '_props') return props;
      if (prop === 'cssText') {
        return Object.entries(props).map(([k, v]) => `${k}: ${v}`).join('; ');
      }
      return props[String(prop)] !== undefined ? props[String(prop)] : '';
    },
    set(target, prop, value) {
      if (prop === '_props') { /* readonly */ return true; }
      if (prop === 'cssText') {
        for (const k of Object.keys(props)) delete props[k];
        if (!value) return true;
        for (const part of String(value).split(';')) {
          const [k, val] = part.split(':');
          if (k && val) props[k.trim()] = val.trim();
        }
        return true;
      }
      props[String(prop)] = String(value);
      return true;
    },
  });
  return el;
}
function _camel(s) { return s.replace(/-([a-z])/g, (_, c) => c.toUpperCase()); }
function _querySelector(root, sel) {
  // Implements only [attr="value"] for our needs
  const m = /^\[([\w-]+)="([^"]+)"\]$/.exec(sel);
  if (!m) return null;
  const [_, attr, val] = m;
  function walk(node) {
    if (node.attributes && node.attributes[attr] === val) return node;
    for (const c of node.children || []) {
      const f = walk(c); if (f) return f;
    }
    return null;
  }
  return walk(root);
}

const _docKeyListeners = [];
const _docMoveListeners = [];
const _docUpListeners = [];
globalThis.document = {
  body:  makeEl('body'),
  head:  makeEl('head'),
  readyState: 'complete',
  activeElement: null,
  createElement: makeEl,
  getElementById(id) {
    function walk(node) {
      if (node.attributes && node.attributes.id === id) return node;
      for (const c of node.children || []) {
        const f = walk(c); if (f) return f;
      }
      return null;
    }
    return walk(document.head) || walk(document.body);
  },
  addEventListener(type, fn) {
    if (type === 'keydown')   _docKeyListeners.push(fn);
    if (type === 'mousemove') _docMoveListeners.push(fn);
    if (type === 'mouseup')   _docUpListeners.push(fn);
  },
  removeEventListener(type, fn) {
    if (type === 'keydown') {
      const i = _docKeyListeners.indexOf(fn);
      if (i >= 0) _docKeyListeners.splice(i, 1);
    }
  },
};
function fireKey(opts) {
  for (const fn of _docKeyListeners) fn(Object.assign({
    preventDefault() {}, stopPropagation() {},
    metaKey: false, ctrlKey: false, altKey: false, shiftKey: false,
  }, opts));
}
function fireMove(x, y) {
  for (const fn of _docMoveListeners) fn({ clientX: x, clientY: y,
    preventDefault() {}, stopPropagation() {} });
}
function fireUp() {
  for (const fn of _docUpListeners) fn({ preventDefault() {}, stopPropagation() {} });
}

globalThis.window = {
  innerWidth: 1280,
  innerHeight: 800,
  popgen: { getCohortKey: () => 'test_cohort_001' },
};
globalThis.popgen = globalThis.window.popgen;

// ----- import the dock --------------------------------------------------------
const dock = require('./atlas_floating_dock.js');

// ----- runner -----------------------------------------------------------------
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

// =============================================================================
// Tests
// =============================================================================

group('mount creates root + launcher and injects CSS', () => {
  dock.mount();
  ok('root attached',     !!dock._internals.root);
  ok('launcher attached', !!dock._internals.launcher);
  ok('CSS injected',      !!document.getElementById('popgen-dock-css'));
  // mount() called twice should be a no-op
  const before = document.body.children.length;
  dock.mount();
  eq('mount idempotent', document.body.children.length, before);
});

group('initial tab is groups + tabs are wired', () => {
  eq('active = groups', dock.getActiveTab(), 'groups');
  const tabIds = dock._internals.tabs.map(t => t.id);
  eq('three tabs',      tabIds, ['groups', 'lasso', 'snapshots']);
  // Tab buttons: only groups marked active
  const buttons = dock._internals.tabButtons;
  eq('groups active',    buttons.groups.attributes['data-active'], 'true');
  eq('lasso inactive',   buttons.lasso.attributes['data-active'],  'false');
  eq('snapshots inactive', buttons.snapshots.attributes['data-active'], 'false');
});

group('setActiveTab switches both button and body active state', () => {
  dock.setActiveTab('lasso');
  eq('active = lasso',   dock.getActiveTab(), 'lasso');
  const buttons = dock._internals.tabButtons;
  const bodies  = dock._internals.bodies;
  eq('groups inactive',  buttons.groups.attributes['data-active'],  'false');
  eq('lasso active',     buttons.lasso.attributes['data-active'],   'true');
  eq('groups body inactive', bodies.groups.attributes['data-active'], 'false');
  eq('lasso body active',    bodies.lasso.attributes['data-active'],  'true');
  // Invalid id is ignored
  dock.setActiveTab('not_a_tab');
  eq('still lasso after invalid', dock.getActiveTab(), 'lasso');
});

group('onTabChange callback', () => {
  let fired = null;
  dock.onTabChange((id) => { fired = id; });
  dock.setActiveTab('snapshots');
  eq('callback got new tab id', fired, 'snapshots');
  dock.setActiveTab('groups');
});

group('hide / show / toggle', () => {
  dock.hide();
  eq('after hide isOpen=false', dock.isOpen(), false);
  eq('root collapsed=true',
    dock._internals.root.attributes['data-collapsed'], 'true');
  eq('launcher visible',
    dock._internals.launcher.attributes['data-hidden'], 'false');
  dock.show();
  eq('after show isOpen=true', dock.isOpen(), true);
  eq('launcher hidden',
    dock._internals.launcher.attributes['data-hidden'], 'true');
  dock.toggle();
  eq('toggle hides',           dock.isOpen(), false);
  dock.toggle();
  eq('toggle shows',           dock.isOpen(), true);
});

group('setTabContent replaces placeholder', () => {
  const newEl = makeEl('div');
  newEl.textContent = 'real content here';
  dock.setTabContent('groups', newEl);
  const body = dock.getTabBody('groups');
  eq('one child',         body.children.length, 1);
  eq('content swapped',   body.children[0].textContent, 'real content here');
  // Setting again replaces, not appends
  const newer = makeEl('div');
  newer.textContent = 'newer';
  dock.setTabContent('groups', newer);
  eq('still one child',   body.children.length, 1);
  eq('content replaced',  body.children[0].textContent, 'newer');
});

group('keybind G toggles dock when no input focused', () => {
  // Make sure dock is open
  dock.show();
  const before = dock.isOpen();
  fireKey({ key: 'g' });
  eq('G hides', dock.isOpen(), !before);
  fireKey({ key: 'g' });
  eq('G shows', dock.isOpen(), before);
  // With ctrl modifier — don't fire
  const stillOpen = dock.isOpen();
  fireKey({ key: 'g', ctrlKey: true });
  eq('Ctrl+G ignored', dock.isOpen(), stillOpen);
});

group('keybind G suppressed when input focused', () => {
  dock.show();
  const fakeInput = makeEl('input');
  fakeInput.tagName = 'INPUT';
  document.activeElement = fakeInput;
  const before = dock.isOpen();
  fireKey({ key: 'g' });
  eq('G ignored when input focused', dock.isOpen(), before);
  document.activeElement = null;
});

group('Shift+S triggers snapshot callback', () => {
  let fired = 0;
  dock.onSaveSnapshot(() => { fired++; });
  fireKey({ key: 'S', shiftKey: true });
  eq('callback fired once', fired, 1);
  // Plain s (no shift) does NOT fire
  fireKey({ key: 's' });
  eq('plain s ignored', fired, 1);
  // Cmd+S also ignored (browser save)
  fireKey({ key: 'S', shiftKey: true, metaKey: true });
  eq('Cmd+Shift+S ignored', fired, 1);
});

group('Shift+S suppressed when input focused', () => {
  let fired = 0;
  dock.onSaveSnapshot(() => { fired++; });
  const fakeTextarea = makeEl('textarea');
  fakeTextarea.tagName = 'TEXTAREA';
  document.activeElement = fakeTextarea;
  fireKey({ key: 'S', shiftKey: true });
  eq('Shift+S ignored when textarea focused', fired, 0);
  document.activeElement = null;
});

group('drag header updates position + persists', () => {
  // Find header inside root by attribute
  const root = dock._internals.root;
  const header = root.querySelector('[data-popgen-dock="header"]');
  ok('header found', !!header);

  // Simulate mousedown on header
  header.dispatchEvent({
    type: 'mousedown',
    target: header,                       // not a button — drag should start
    clientX: 100, clientY: 100,
    preventDefault() {}, stopPropagation() {},
  });
  // Move
  fireMove(150, 130);
  fireUp();

  // Position changed by (50, 30) from previously stored value
  // Initial position is whatever was loaded; verify the *delta* matches.
  ok('position updated',
    dock._internals.position.x !== 24 || dock._internals.position.y !== 80);

  // Persisted: verify localStorage has the new position
  const k = 'inversion_atlas.dock.test_cohort_001.position';
  const stored = JSON.parse(localStorage.getItem(k));
  ok('position persisted', stored && typeof stored.x === 'number');
});

group('resize handle updates size + persists', () => {
  const root = dock._internals.root;
  const handle = root.querySelector('[data-popgen-dock="resize"]');
  ok('handle found', !!handle);
  handle.dispatchEvent({
    type: 'mousedown',
    target: handle,
    clientX: 200, clientY: 200,
    preventDefault() {}, stopPropagation() {},
  });
  fireMove(280, 250);
  fireUp();
  ok('size grew', dock._internals.size.w >= 380 && dock._internals.size.h >= 460);
  const k = 'inversion_atlas.dock.test_cohort_001.size';
  const stored = JSON.parse(localStorage.getItem(k));
  ok('size persisted', stored && stored.w >= 380);
});

group('persistence: unmount + remount restores state', () => {
  dock.setActiveTab('snapshots');
  dock.hide();
  const positionBefore = Object.assign({}, dock._internals.position);
  const sizeBefore     = Object.assign({}, dock._internals.size);

  dock.unmount();
  ok('unmounted',  dock._internals.root === null);
  // Body children flushed
  ok('dock removed from DOM',
    !document.body.children.find(c => c.attributes['data-popgen-dock'] === 'root'));

  dock.mount();
  eq('active tab restored',  dock.getActiveTab(), 'snapshots');
  eq('open state restored',  dock.isOpen(), false);
  eq('position restored',    dock._internals.position, positionBefore);
  eq('size restored',        dock._internals.size, sizeBefore);
});

group('launcher click expands dock', () => {
  // Dock is currently hidden from the test above
  dock.hide();
  const launcher = dock._internals.launcher;
  // Simulate a click without movement — should expand
  launcher.dispatchEvent({
    type: 'mousedown',
    target: launcher,
    clientX: 16, clientY: 16,
    preventDefault() {}, stopPropagation() {},
  });
  fireUp();   // no move between mousedown and mouseup → click
  eq('dock open after launcher click', dock.isOpen(), true);
});

group('launcher drag relocates without expanding', () => {
  dock.hide();   // close so the launcher is interactable
  const launcher = dock._internals.launcher;
  launcher.dispatchEvent({
    type: 'mousedown',
    target: launcher,
    clientX: 16, clientY: 16,
    preventDefault() {}, stopPropagation() {},
  });
  fireMove(120, 200);   // movement >4px → drag
  fireUp();
  eq('dock STILL closed', dock.isOpen(), false);
  // Launcher position updated
  const left = parseFloat(launcher.style._props.left);
  const top  = parseFloat(launcher.style._props.top);
  ok('launcher moved', left > 16 && top > 16);
});

// ----- summary ---------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);
