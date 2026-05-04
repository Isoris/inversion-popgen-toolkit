// =============================================================================
// atlas_floating_dock.js — Turn 5a of chat A
// =============================================================================
//
// The floating dock — a small, always-accessible panel above any atlas page.
// Three tabs:
//   1. Groups          — group composition editor (turn 5b)
//   2. Lasso history   — selection log with pin/restore (turn 5b)
//   3. Snapshots       — full-state captures with restore (turn 5b)
//
// This file delivers the SHELL only:
//   - mount/unmount lifecycle
//   - draggable header
//   - tab-bar with three tabs whose bodies are placeholders for turn 5b
//   - keyboard shortcuts: G (toggle), Shift+S (save snapshot when 5b lands)
//   - position + size + open/closed + active-tab persistence (localStorage)
//   - resize handle (corner)
//   - close button (×) collapses the dock to a tiny "G" launcher
//   - launcher is itself draggable independently of the dock
//   - z-index 9000 — below atlas's 9999/10000 modals
//
// API (window.popgenDock):
//   .mount()                — attach to document.body, restore persisted state
//   .unmount()              — detach (rare; dev/test only)
//   .show()                 — open dock, hide launcher
//   .hide()                 — close dock, show launcher
//   .toggle()               — flip
//   .setActiveTab(id)       — 'groups' | 'lasso' | 'snapshots'
//   .getActiveTab()
//   .setTabContent(id, el)  — used by turn 5b to inject tab content
//   .getTabBody(id)         — return the DOM element for a tab body
//   .onSaveSnapshot(fn)     — register callback for Shift+S
//   .onTabChange(fn)        — register callback for tab change
//   .isOpen()
//
// Persistence keys (under inversion_atlas.dock.<cohort>):
//   .position  — { x, y }
//   .size      — { w, h }
//   .open      — bool
//   .activeTab — id
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenDock = factory();
    // Auto-mount in browser if document is already loaded; otherwise wait.
    if (typeof document !== 'undefined') {
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', () => root.popgenDock.mount());
      } else {
        root.popgenDock.mount();
      }
    }
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const LS_PREFIX = 'inversion_atlas.dock.';

  const DEFAULT_POSITION = { x: 24, y: 80 };
  const DEFAULT_SIZE     = { w: 380, h: 460 };
  const MIN_SIZE         = { w: 280, h: 240 };
  const LAUNCHER_DEFAULT_POS = { x: 16, y: 16 };

  const TABS = [
    { id: 'groups',     label: 'Groups',  hint: 'Compose groups + Compute (G)' },
    { id: 'lasso',      label: 'Lasso',   hint: 'Selection history' },
    { id: 'snapshots',  label: 'Snaps',   hint: 'Save/restore state (Shift+S)' },
  ];

  const Z_INDEX_DOCK     = 9000;
  const Z_INDEX_LAUNCHER = 9001;  // launcher slightly above so it's clickable
                                  // even if the dock is at min size

  // ===========================================================================
  // State
  // ===========================================================================

  let _root = null;          // dock root element
  let _launcher = null;      // collapsed-state launcher element
  let _bodies = {};          // tab id → body element
  let _tabButtons = {};      // tab id → button element
  let _active = 'groups';
  let _open = true;
  let _position = Object.assign({}, DEFAULT_POSITION);
  let _size = Object.assign({}, DEFAULT_SIZE);
  let _launcherPos = Object.assign({}, LAUNCHER_DEFAULT_POS);
  let _mounted = false;

  // Listener registries
  const _saveSnapshotCallbacks = [];
  const _tabChangeCallbacks    = [];

  // ===========================================================================
  // Cohort-scoped persistence
  // ===========================================================================
  // Each cohort gets its own dock-state record so a user's preferred dock
  // position survives chrom switches but is independent across cohorts. Falls
  // back to a single shared key if the cohort isn't known.

  function _getCohortKey() {
    if (typeof window !== 'undefined' && window.popgen &&
        typeof window.popgen.getCohortKey === 'function') {
      try { return window.popgen.getCohortKey(); } catch (_) {}
    }
    if (typeof globalThis !== 'undefined' && globalThis.popgen &&
        typeof globalThis.popgen.getCohortKey === 'function') {
      try { return globalThis.popgen.getCohortKey(); } catch (_) {}
    }
    return 'default';
  }

  function _lsKey(suffix) {
    return LS_PREFIX + _getCohortKey() + '.' + suffix;
  }

  function _lsGet(suffix, fallback) {
    if (typeof localStorage === 'undefined') return fallback;
    try {
      const raw = localStorage.getItem(_lsKey(suffix));
      if (raw == null) return fallback;
      return JSON.parse(raw);
    } catch (_) { return fallback; }
  }

  function _lsSet(suffix, value) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(_lsKey(suffix), JSON.stringify(value)); } catch (_) {}
  }

  function _persistAll() {
    _lsSet('position',     _position);
    _lsSet('size',         _size);
    _lsSet('open',         _open);
    _lsSet('activeTab',    _active);
    _lsSet('launcherPos',  _launcherPos);
  }

  function _loadPersisted() {
    _position    = _lsGet('position',    DEFAULT_POSITION);
    _size        = _lsGet('size',        DEFAULT_SIZE);
    _open        = _lsGet('open',        true);
    _active      = _lsGet('activeTab',   'groups');
    _launcherPos = _lsGet('launcherPos', LAUNCHER_DEFAULT_POS);
    // Sanity: clamp size
    if (!_size || typeof _size.w !== 'number') _size = Object.assign({}, DEFAULT_SIZE);
    if (_size.w < MIN_SIZE.w) _size.w = MIN_SIZE.w;
    if (_size.h < MIN_SIZE.h) _size.h = MIN_SIZE.h;
    if (TABS.findIndex(t => t.id === _active) < 0) _active = 'groups';
  }

  // ===========================================================================
  // CSS injection (one stylesheet, scoped under [data-popgen-dock])
  // ===========================================================================

  const DOCK_CSS = `
    [data-popgen-dock="root"] {
      position: fixed; z-index: ${Z_INDEX_DOCK};
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      box-shadow: 0 6px 24px rgba(0,0,0,0.18);
      border-radius: 8px;
      font-family: var(--mono, ui-monospace, "SF Mono", monospace);
      font-size: 11px;
      color: var(--ink, #0e1116);
      display: flex; flex-direction: column;
      overflow: hidden;
      user-select: none;
    }
    [data-popgen-dock="root"][data-collapsed="true"] {
      display: none;
    }
    [data-popgen-dock="header"] {
      display: flex; align-items: center; gap: 6px;
      padding: 6px 10px; height: 26px;
      background: var(--panel-2, #f0f1f3);
      border-bottom: 1px solid var(--rule, #c8cdd2);
      cursor: move;
      flex-shrink: 0;
    }
    [data-popgen-dock="title"] {
      font-weight: 600;
      letter-spacing: 0.04em; text-transform: uppercase;
      font-size: 10px;
      color: var(--ink-dim, #555e69);
      flex: 1;
      pointer-events: none;
    }
    [data-popgen-dock="kbd"] {
      font-size: 9px;
      color: var(--ink-dimmer, #8a909a);
      padding: 0 4px;
      pointer-events: none;
    }
    [data-popgen-dock="close"] {
      width: 18px; height: 18px;
      border: 0; background: transparent;
      color: var(--ink-dim, #555e69);
      cursor: pointer; font-size: 14px; line-height: 1;
      border-radius: 4px;
      padding: 0;
    }
    [data-popgen-dock="close"]:hover {
      background: var(--bad, #e0555c); color: white;
    }
    [data-popgen-dock="tabbar"] {
      display: flex;
      background: var(--panel-2, #f0f1f3);
      border-bottom: 1px solid var(--rule, #c8cdd2);
      flex-shrink: 0;
    }
    [data-popgen-dock="tab"] {
      flex: 1;
      padding: 6px 8px;
      border: 0; background: transparent;
      color: var(--ink-dim, #555e69);
      cursor: pointer;
      font-family: inherit; font-size: 10px;
      font-weight: 600;
      letter-spacing: 0.05em; text-transform: uppercase;
      border-right: 1px solid var(--rule, #c8cdd2);
      border-bottom: 2px solid transparent;
    }
    [data-popgen-dock="tab"]:last-child { border-right: 0; }
    [data-popgen-dock="tab"]:hover {
      background: var(--panel-3, #e6e8ec);
      color: var(--ink, #0e1116);
    }
    [data-popgen-dock="tab"][data-active="true"] {
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-bottom-color: var(--accent, #f5a524);
    }
    [data-popgen-dock="bodies"] {
      flex: 1;
      overflow: hidden;
      position: relative;
    }
    [data-popgen-dock="body"] {
      position: absolute; inset: 0;
      overflow: auto;
      padding: 10px 12px;
      display: none;
    }
    [data-popgen-dock="body"][data-active="true"] {
      display: block;
    }
    [data-popgen-dock="body"] .placeholder {
      color: var(--ink-dim, #555e69);
      font-style: italic;
      padding: 18px 4px;
      line-height: 1.5;
    }
    [data-popgen-dock="resize"] {
      position: absolute;
      bottom: 0; right: 0;
      width: 14px; height: 14px;
      cursor: nwse-resize;
      background:
        linear-gradient(135deg,
          transparent 0%, transparent 50%,
          var(--ink-dim, #555e69) 50%, var(--ink-dim, #555e69) 60%,
          transparent 60%, transparent 70%,
          var(--ink-dim, #555e69) 70%, var(--ink-dim, #555e69) 80%,
          transparent 80%);
      opacity: 0.4;
    }
    [data-popgen-dock="resize"]:hover { opacity: 0.8; }

    [data-popgen-launcher="root"] {
      position: fixed; z-index: ${Z_INDEX_LAUNCHER};
      width: 32px; height: 32px;
      border-radius: 16px;
      background: var(--accent, #f5a524);
      color: #0e1116;
      display: flex; align-items: center; justify-content: center;
      cursor: move;
      font-family: var(--mono, ui-monospace, monospace);
      font-size: 14px; font-weight: 700;
      box-shadow: 0 3px 10px rgba(0,0,0,0.3);
      border: 0;
      user-select: none;
    }
    [data-popgen-launcher="root"][data-hidden="true"] { display: none; }
    [data-popgen-launcher="root"]:hover {
      box-shadow: 0 4px 14px rgba(0,0,0,0.4);
      transform: scale(1.05);
    }
  `;

  function _injectCSS() {
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-dock-css')) return;
    const style = document.createElement('style');
    style.id = 'popgen-dock-css';
    style.textContent = DOCK_CSS;
    document.head.appendChild(style);
  }

  // ===========================================================================
  // Drag + resize helpers
  // ===========================================================================

  function _makeDraggable(handle, root, onMove) {
    if (!handle || !root) return;
    let startX = 0, startY = 0, baseX = 0, baseY = 0, dragging = false;
    handle.addEventListener('mousedown', (e) => {
      // Don't start drag from interactive children (buttons, inputs)
      if (e.target.tagName === 'BUTTON' || e.target.tagName === 'INPUT') return;
      dragging = true;
      startX = e.clientX; startY = e.clientY;
      const rect = root.getBoundingClientRect();
      baseX = rect.left; baseY = rect.top;
      e.preventDefault();
    });
    document.addEventListener('mousemove', (e) => {
      if (!dragging) return;
      let nx = baseX + (e.clientX - startX);
      let ny = baseY + (e.clientY - startY);
      // Clamp to viewport so dock can't be dragged offscreen entirely
      const margin = 24;
      const w = root.offsetWidth, h = root.offsetHeight;
      const vw = window.innerWidth, vh = window.innerHeight;
      nx = Math.max(-w + margin, Math.min(vw - margin, nx));
      ny = Math.max(0, Math.min(vh - margin, ny));
      root.style.left = nx + 'px';
      root.style.top  = ny + 'px';
      if (onMove) onMove(nx, ny);
    });
    document.addEventListener('mouseup', () => {
      if (dragging) { dragging = false; }
    });
  }

  function _makeResizable(handle, root, onResize) {
    if (!handle || !root) return;
    let startX = 0, startY = 0, baseW = 0, baseH = 0, resizing = false;
    handle.addEventListener('mousedown', (e) => {
      resizing = true;
      startX = e.clientX; startY = e.clientY;
      baseW = root.offsetWidth; baseH = root.offsetHeight;
      e.preventDefault(); e.stopPropagation();
    });
    document.addEventListener('mousemove', (e) => {
      if (!resizing) return;
      const nw = Math.max(MIN_SIZE.w, baseW + (e.clientX - startX));
      const nh = Math.max(MIN_SIZE.h, baseH + (e.clientY - startY));
      root.style.width  = nw + 'px';
      root.style.height = nh + 'px';
      if (onResize) onResize(nw, nh);
    });
    document.addEventListener('mouseup', () => {
      if (resizing) { resizing = false; }
    });
  }

  // ===========================================================================
  // DOM construction
  // ===========================================================================

  function _buildDock() {
    const root = document.createElement('div');
    root.setAttribute('data-popgen-dock', 'root');
    root.style.left   = _position.x + 'px';
    root.style.top    = _position.y + 'px';
    root.style.width  = _size.w + 'px';
    root.style.height = _size.h + 'px';
    root.setAttribute('data-collapsed', _open ? 'false' : 'true');

    // Header
    const header = document.createElement('div');
    header.setAttribute('data-popgen-dock', 'header');
    const title = document.createElement('span');
    title.setAttribute('data-popgen-dock', 'title');
    title.textContent = 'popgen';
    const kbd = document.createElement('span');
    kbd.setAttribute('data-popgen-dock', 'kbd');
    kbd.textContent = 'G · Shift+S';
    const close = document.createElement('button');
    close.setAttribute('data-popgen-dock', 'close');
    close.textContent = '×';
    close.title = 'Collapse to launcher (G)';
    close.addEventListener('click', () => hide());
    header.appendChild(title);
    header.appendChild(kbd);
    header.appendChild(close);

    // Tab bar
    const tabbar = document.createElement('div');
    tabbar.setAttribute('data-popgen-dock', 'tabbar');
    for (const t of TABS) {
      const btn = document.createElement('button');
      btn.setAttribute('data-popgen-dock', 'tab');
      btn.setAttribute('data-tab-id', t.id);
      btn.setAttribute('data-active', String(t.id === _active));
      btn.title = t.hint;
      btn.textContent = t.label;
      btn.addEventListener('click', () => setActiveTab(t.id));
      _tabButtons[t.id] = btn;
      tabbar.appendChild(btn);
    }

    // Bodies
    const bodies = document.createElement('div');
    bodies.setAttribute('data-popgen-dock', 'bodies');
    for (const t of TABS) {
      const body = document.createElement('div');
      body.setAttribute('data-popgen-dock', 'body');
      body.setAttribute('data-tab-id', t.id);
      body.setAttribute('data-active', String(t.id === _active));
      // Placeholder until turn 5b injects content
      const ph = document.createElement('div');
      ph.className = 'placeholder';
      ph.textContent = '— ' + t.label + ' tab — content lands in turn 5b. ' +
                       'Engine state is reachable now via window.popgen.';
      body.appendChild(ph);
      _bodies[t.id] = body;
      bodies.appendChild(body);
    }

    // Resize handle
    const resize = document.createElement('div');
    resize.setAttribute('data-popgen-dock', 'resize');
    resize.title = 'Drag to resize';

    root.appendChild(header);
    root.appendChild(tabbar);
    root.appendChild(bodies);
    root.appendChild(resize);

    _makeDraggable(header, root, (x, y) => {
      _position.x = x; _position.y = y;
      _persistAll();
    });
    _makeResizable(resize, root, (w, h) => {
      _size.w = w; _size.h = h;
      _persistAll();
    });

    return root;
  }

  function _buildLauncher() {
    const launcher = document.createElement('button');
    launcher.setAttribute('data-popgen-launcher', 'root');
    launcher.setAttribute('data-hidden', _open ? 'true' : 'false');
    launcher.style.left = _launcherPos.x + 'px';
    launcher.style.top  = _launcherPos.y + 'px';
    launcher.textContent = 'G';
    launcher.title = 'Open popgen dock (G)';
    // Click expands; drag relocates. Distinguish via a small drag-threshold.
    let downX = 0, downY = 0, dragged = false, dragging = false;
    let startLX = 0, startLY = 0;
    launcher.addEventListener('mousedown', (e) => {
      downX = e.clientX; downY = e.clientY;
      dragged = false; dragging = true;
      const rect = launcher.getBoundingClientRect();
      startLX = rect.left; startLY = rect.top;
      e.preventDefault();
    });
    document.addEventListener('mousemove', (e) => {
      if (!dragging) return;
      const dx = e.clientX - downX, dy = e.clientY - downY;
      if (!dragged && Math.hypot(dx, dy) > 4) dragged = true;
      if (dragged) {
        let nx = startLX + dx, ny = startLY + dy;
        const m = 4, w = launcher.offsetWidth, h = launcher.offsetHeight;
        nx = Math.max(m, Math.min(window.innerWidth  - w - m, nx));
        ny = Math.max(m, Math.min(window.innerHeight - h - m, ny));
        launcher.style.left = nx + 'px';
        launcher.style.top  = ny + 'px';
        _launcherPos.x = nx; _launcherPos.y = ny;
      }
    });
    document.addEventListener('mouseup', () => {
      if (!dragging) return;
      dragging = false;
      if (dragged) { _persistAll(); }
      else         { show(); }
    });
    return launcher;
  }

  // ===========================================================================
  // Public lifecycle
  // ===========================================================================

  function mount() {
    if (_mounted) return;
    if (typeof document === 'undefined') return;
    _injectCSS();
    _loadPersisted();
    _root     = _buildDock();
    _launcher = _buildLauncher();
    document.body.appendChild(_root);
    document.body.appendChild(_launcher);
    _bindGlobalKeys();
    _mounted = true;
  }

  function unmount() {
    if (!_mounted) return;
    _unbindGlobalKeys();
    if (_root && _root.parentNode)         _root.parentNode.removeChild(_root);
    if (_launcher && _launcher.parentNode) _launcher.parentNode.removeChild(_launcher);
    _root = null; _launcher = null;
    _bodies = {}; _tabButtons = {};
    _mounted = false;
  }

  function show() {
    _open = true;
    if (_root)     _root.setAttribute('data-collapsed', 'false');
    if (_launcher) _launcher.setAttribute('data-hidden', 'true');
    _persistAll();
  }

  function hide() {
    _open = false;
    if (_root)     _root.setAttribute('data-collapsed', 'true');
    if (_launcher) _launcher.setAttribute('data-hidden', 'false');
    _persistAll();
  }

  function toggle() { if (_open) hide(); else show(); }

  function isOpen() { return _open; }

  function setActiveTab(id) {
    if (TABS.findIndex(t => t.id === id) < 0) return;
    _active = id;
    for (const tid of Object.keys(_tabButtons)) {
      _tabButtons[tid].setAttribute('data-active', String(tid === id));
    }
    for (const tid of Object.keys(_bodies)) {
      _bodies[tid].setAttribute('data-active', String(tid === id));
    }
    _persistAll();
    for (const fn of _tabChangeCallbacks) {
      try { fn(id); } catch (_) {}
    }
  }

  function getActiveTab() { return _active; }

  function setTabContent(id, el) {
    if (!_bodies[id] || !el) return;
    // Replace placeholder/whatever with the new element
    while (_bodies[id].firstChild) _bodies[id].removeChild(_bodies[id].firstChild);
    _bodies[id].appendChild(el);
  }

  function getTabBody(id) { return _bodies[id] || null; }

  function onSaveSnapshot(fn) {
    if (typeof fn === 'function') _saveSnapshotCallbacks.push(fn);
  }

  function onTabChange(fn) {
    if (typeof fn === 'function') _tabChangeCallbacks.push(fn);
  }

  // ===========================================================================
  // Keyboard shortcuts
  // ===========================================================================

  function _isInputFocused() {
    if (typeof document === 'undefined') return false;
    const a = document.activeElement;
    if (!a) return false;
    const tag = (a.tagName || '').toUpperCase();
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return true;
    if (a.isContentEditable) return true;
    return false;
  }

  let _keyHandler = null;
  function _bindGlobalKeys() {
    if (typeof document === 'undefined') return;
    _keyHandler = (e) => {
      if (_isInputFocused()) return;
      // G — toggle dock
      if (e.key === 'g' || e.key === 'G') {
        if (e.metaKey || e.ctrlKey || e.altKey) return;
        if (e.shiftKey) return;
        toggle();
        e.preventDefault();
        return;
      }
      // Shift+S — save snapshot
      if ((e.key === 'S' || e.key === 's') && e.shiftKey &&
          !e.metaKey && !e.ctrlKey && !e.altKey) {
        for (const fn of _saveSnapshotCallbacks) {
          try { fn(); } catch (_) {}
        }
        e.preventDefault();
        return;
      }
    };
    document.addEventListener('keydown', _keyHandler);
  }

  function _unbindGlobalKeys() {
    if (typeof document === 'undefined' || !_keyHandler) return;
    document.removeEventListener('keydown', _keyHandler);
    _keyHandler = null;
  }

  // ===========================================================================
  // Public API
  // ===========================================================================

  return {
    mount, unmount, show, hide, toggle, isOpen,
    setActiveTab, getActiveTab,
    setTabContent, getTabBody,
    onSaveSnapshot, onTabChange,

    // Internals exposed for tests
    _internals: {
      get root()      { return _root; },
      get launcher()  { return _launcher; },
      get bodies()    { return _bodies; },
      get tabButtons(){ return _tabButtons; },
      get position()  { return _position; },
      get size()      { return _size; },
      get tabs()      { return TABS.slice(); },
      _persistAll,
      _loadPersisted,
    },
  };
}));
