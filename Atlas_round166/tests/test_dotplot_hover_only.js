// =============================================================================
// Smoke test for atlas_dotplot.js after pin removal.
//
// Verifies:
//  1. The module loads without syntax errors and exposes makeDotplotPanel
//  2. makeDotplotPanel called with a tiny fixture returns an HTMLElement-like
//     wrapper without throwing
//  3. The popup state object no longer has a `pinned` field
//  4. Calling closePopup repeatedly is idempotent (no document.removeEventListener
//     errors from the now-deleted onDocClick reference)
// =============================================================================

const fs = require('fs');
const path = require('path');

// Minimal DOM shim
function makeEl(tag) {
  const el = {
    tagName: tag.toUpperCase(),
    children: [],
    parentNode: null,
    _innerHTML: '',
    style: {},
    attributes: {},
    classList: {
      _cls: new Set(),
      add(c) { this._cls.add(c); },
      remove(c) { this._cls.delete(c); },
      contains(c) { return this._cls.has(c); },
    },
    _listeners: {},
    appendChild(c) { this.children.push(c); c.parentNode = this; return c; },
    removeChild(c) {
      const i = this.children.indexOf(c);
      if (i >= 0) { this.children.splice(i, 1); c.parentNode = null; }
      return c;
    },
    insertAdjacentHTML(_pos, html) { this._innerHTML += html; },
    addEventListener(ev, fn) { (this._listeners[ev] = this._listeners[ev] || []).push(fn); },
    removeEventListener(ev, fn) {
      const arr = this._listeners[ev] || [];
      const i = arr.indexOf(fn);
      if (i >= 0) arr.splice(i, 1);
    },
    setAttribute(k, v) { this.attributes[k] = String(v); },
    getAttribute(k) { return this.attributes[k] != null ? this.attributes[k] : null; },
    querySelector() { return null; },
    querySelectorAll() { return []; },
    getBoundingClientRect() { return { left: 0, top: 0, width: 300, height: 300 }; },
    contains(other) { return other === this; },
    get innerHTML() { return this._innerHTML; },
    set innerHTML(v) { this._innerHTML = v; this.children = []; },
    get offsetWidth() { return 200; },
    get offsetHeight() { return 40; },
    dispatchEvent() {},
  };
  return el;
}
const document = {
  body: makeEl('body'),
  _listeners: {},
  createElement(tag) { return makeEl(tag); },
  addEventListener(ev, fn) {
    (this._listeners[ev] = this._listeners[ev] || []).push(fn);
  },
  removeEventListener(ev, fn) {
    const arr = this._listeners[ev] || [];
    const i = arr.indexOf(fn);
    if (i >= 0) arr.splice(i, 1);
  },
};
const window = {};
const localStorage = {
  _store: {},
  getItem(k) { return this._store[k] != null ? this._store[k] : null; },
  setItem(k, v) { this._store[k] = String(v); },
  removeItem(k) { delete this._store[k]; },
};

// Load module
const modSrc = fs.readFileSync(
  path.join(__dirname, 'atlas_dotplot.js'), 'utf8');
new Function('window', 'document', 'localStorage', 'self',
             modSrc + '\nreturn window.popgenDotplot;')(
             window, document, localStorage, window);

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// --- Test 1: module loaded
ok('window.popgenDotplot defined', !!window.popgenDotplot);
ok('makeDotplotPanel exposed',
   typeof window.popgenDotplot.makeDotplotPanel === 'function');

// --- Test 2: tiny fixture renders without throwing
const fixture = {
  species_query:  { name: 'Cgar' },
  species_target: { name: 'Cmac' },
  chrom_lengths_query:  { 'LG28': 22_000_000 },
  chrom_lengths_target: { 'LG28': 21_500_000, 'LG12': 18_000_000 },
  wfmash_blocks: [
    { q_chr: 'LG28', q_start: 1_000_000, q_end: 5_000_000,
      t_chr: 'LG28', t_start: 1_500_000, t_end: 5_500_000,
      strand: '+', pi: 0.92 },
    { q_chr: 'LG28', q_start: 8_000_000, q_end: 12_000_000,
      t_chr: 'LG28', t_start: 9_000_000, t_end: 13_000_000,
      strand: '+', pi: 0.88 },
  ],
};
let panel = null, threw = null;
try {
  panel = window.popgenDotplot.makeDotplotPanel({
    data: fixture,
    defaultResolution: 'wfmash',
    miniSize: 200,
    enlargedSize: 500,
  });
} catch (e) { threw = e; }
ok('makeDotplotPanel does not throw', threw === null,
   threw && threw.message);
ok('returns a wrapper element', panel && typeof panel === 'object');

// --- Test 3: confirm pin removal at the source level by re-reading the file
const modText = fs.readFileSync(
  path.join(__dirname, 'atlas_dotplot.js'), 'utf8');
ok('no `pinned` field references anywhere in module',
   modText.indexOf('pinned') < 0,
   'still present at offset ' + modText.indexOf('pinned'));
ok('no `onDocClick` references',
   modText.indexOf('onDocClick') < 0);
ok('no "click to pin" hint text',
   modText.indexOf('click to pin') < 0);
ok('hover-only hint text present',
   modText.indexOf("'hover to enlarge'") >= 0);
ok('PINNED state diagram comment removed',
   modText.indexOf('PINNED ──×─or─outside-click') < 0);

// --- Test 4: simulate hover open + leave close (the new pure-hover flow)
// Find the miniHolder by its cssText (DOM shim stores raw cssText, not
// parsed properties).
const miniHolder = panel && panel.children.find(c =>
  c && c.style && typeof c.style.cssText === 'string' &&
  c.style.cssText.indexOf('cursor:pointer') >= 0);
ok('miniHolder found in wrap.children (cursor:pointer)',
   miniHolder != null,
   'children: ' + (panel ? panel.children.length : 0));

// Trigger mouseenter — should call renderEnlarged → opens popup
let openOK = false;
try {
  const enterFns = miniHolder._listeners['mouseenter'] || [];
  for (const fn of enterFns) fn();
  // The popup is appended to document.body. Check that body got a child.
  openOK = document.body.children.length > 0;
} catch (e) { threw = e; }
ok('mouseenter opens popup (body has new child)', openOK,
   'body.children.length=' + document.body.children.length);

// Trigger mouseleave — schedules a 80ms close timer
let leaveOK = false;
try {
  const leaveFns = miniHolder._listeners['mouseleave'] || [];
  for (const fn of leaveFns) fn();
  leaveOK = true;
} catch (e) { threw = e; }
ok('mouseleave does not throw', leaveOK, threw && threw.message);

// Wait for the 80ms close timer + a small margin, then check the popup is gone
setTimeout(() => {
  // Trigger mouseover outside both mini and popup to flush the close path
  const overFns = (document._listeners['mouseover'] || []);
  // Body is "outside" both mini and popup
  const fakeEvent = { target: document.body };
  for (const fn of overFns) fn(fakeEvent);
  // Now wait another 100ms for the timer
  setTimeout(() => {
    ok('popup eventually closes after mouseleave',
       document.body.children.length === 0,
       'body.children.length=' + document.body.children.length);

    // --- Test 5: idempotent close — calling close handlers multiple times
    // shouldn't crash now that onDocClick is gone (it would have been a
    // ReferenceError if any code path still referenced it).
    let multiCloseOK = true;
    try {
      // Re-open and re-close several times
      for (let i = 0; i < 3; i++) {
        const enterFns = miniHolder._listeners['mouseenter'] || [];
        for (const fn of enterFns) fn();
        const leaveFns = miniHolder._listeners['mouseleave'] || [];
        for (const fn of leaveFns) fn();
      }
    } catch (e) {
      multiCloseOK = false;
      console.log('  multi-close caught:', e.message);
    }
    ok('open/close cycles do not throw', multiCloseOK);

    // Print results
    console.log('\n=========================================');
    console.log('passed: ' + pass + ' / failed: ' + fail);
    console.log('=========================================');
    process.exit(fail > 0 ? 1 : 0);
  }, 100);
}, 200);
