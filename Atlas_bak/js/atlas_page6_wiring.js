// =============================================================================
// atlas_page6_wiring.js — Turn 6 of chat A
// =============================================================================
//
// Wires page 6 (popstats stack) to the live-server response stashed in
// state.popstatsLive.lastResponse. Adds:
//
//   1. wrapTrackDef(trackDef) — returns a clone of the placeholder track def
//      with a getData() that resolves static / live / split based on the
//      popstats mode chip from turn 4.
//   2. wrapAllTrackDefs(trackList) — convenience over an array.
//   3. attachTooltip(trackDef, canvas) — install pointer-tracking on a single
//      multi-line track canvas. Hover shows per-series values at x.
//   4. installSlotCombineUI() — adds a "combine" button to each empty slot
//      row in the dock's Groups tab. Opens a small picker for minus / union /
//      intersect of existing slots.
//
// API (window.popgenPage6):
//   .wrapTrackDef(t)
//   .wrapAllTrackDefs(list)
//   .attachTooltip(t, canvas)
//   .installSlotCombineUI()
//   .invalidatePage6() — call this after Compute lands a new response so the
//                        atlas re-renders the popstats stack
//
// SPLICE: see SPLICE_POINTS.md. The atlas's collectPopstatsTracks() should
// call wrapAllTrackDefs() on its output before passing to drawPopstatsTracks().
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenPage6 = factory();
    if (typeof document !== 'undefined') {
      // Auto-listen for live computes — when the dock fires a Compute, we
      // need to ask the atlas to re-render the popstats stack.
      const tryAutoWire = () => {
        if (root.popgenLive && root.popgenRenderers) {
          // Patch popgenLive.popstatsGroupwise so we re-render after each call
          _patchLiveLayer(root.popgenLive);
          // Listen for mode changes so static<->live re-renders
          if (typeof window !== 'undefined' && window.addEventListener) {
            window.addEventListener('popgen:mode-changed', () => _triggerRerender());
            window.addEventListener('popgen:multiline-mode-changed',
              () => _triggerRerender());
          }
          // Install the slot-combine UI once the dock tabs exist
          const tryInstallSlotCombine = () => {
            if (root.popgenDockTabs) {
              root.popgenPage6.installSlotCombineUI();
            } else {
              setTimeout(tryInstallSlotCombine, 100);
            }
          };
          tryInstallSlotCombine();
        } else {
          setTimeout(tryAutoWire, 100);
        }
      };
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', tryAutoWire);
      } else { tryAutoWire(); }
    }
  }

  function _patchLiveLayer(live) {
    if (live.__page6Patched) return;
    const orig = live.popstatsGroupwise;
    live.popstatsGroupwise = async function (...args) {
      const env = await orig.apply(this, args);
      if (env && env.ok) _triggerRerender();
      return env;
    };
    // Same wrapping for hobs and ancestry
    if (live.hobsGroupwise) {
      const oh = live.hobsGroupwise;
      live.hobsGroupwise = async function (...args) {
        const env = await oh.apply(this, args);
        if (env && env.ok) _triggerRerender();
        return env;
      };
    }
    if (live.ancestryGroupwiseQ) {
      const oa = live.ancestryGroupwiseQ;
      live.ancestryGroupwiseQ = async function (...args) {
        const env = await oa.apply(this, args);
        if (env && env.ok) _triggerRerender();
        return env;
      };
    }
    live.__page6Patched = true;
  }

  function _triggerRerender() {
    // Defer to next tick so the live envelope is fully stashed in state
    if (typeof setTimeout === 'undefined') return;
    setTimeout(() => {
      if (typeof window === 'undefined') return;
      // Best-effort: call atlas's renderPopstatsPage if exposed
      if (typeof window.renderPopstatsPage === 'function') {
        try { window.renderPopstatsPage(); } catch (_) {}
      } else if (window.dispatchEvent) {
        window.dispatchEvent(new CustomEvent('popgen:page6-rerender-requested'));
      }
    }, 0);
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Atlas / engine / live / renderer accessors
  // ===========================================================================

  function _renderers() {
    if (typeof window !== 'undefined' && window.popgenRenderers) return window.popgenRenderers;
    if (typeof globalThis !== 'undefined' && globalThis.popgenRenderers) return globalThis.popgenRenderers;
    return null;
  }
  function _engine() {
    if (typeof window !== 'undefined' && window.popgen) return window.popgen;
    if (typeof globalThis !== 'undefined' && globalThis.popgen) return globalThis.popgen;
    return null;
  }
  function _atlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }
  function _dockTabs() {
    if (typeof window !== 'undefined' && window.popgenDockTabs) return window.popgenDockTabs;
    if (typeof globalThis !== 'undefined' && globalThis.popgenDockTabs) return globalThis.popgenDockTabs;
    return null;
  }

  // ===========================================================================
  // Static-mode resolvers — read from existing state.data.tracks / etc.
  // ===========================================================================
  // For each known track id, a fallback resolver that produces data from
  // precomp JSON. These match the atlas's existing static-mode shapes; the
  // wrapper falls back to these when mode='static' or live data is absent.

  function _staticResolver(trackId) {
    return (s) => {
      if (!s || !s.data) return null;
      // Tracks that mirror values straight from state.data.tracks dict
      if (s.data.tracks && s.data.tracks[trackId]) {
        const t = s.data.tracks[trackId];
        if (!Array.isArray(t.values)) return null;
        let mb;
        if (Array.isArray(t.pos_bp) && t.pos_bp.length === t.values.length) {
          mb = t.pos_bp.map(bp => bp / 1e6);
        } else if (s.data.windows && s.data.windows.length === t.values.length) {
          mb = s.data.windows.map(w => w.center_mb);
        } else { return null; }
        return { mb, values: t.values, min: t.min, max: t.max };
      }
      return null;
    };
  }

  // ===========================================================================
  // Live-mode resolvers — read from state.popstatsLive.lastResponse
  // ===========================================================================

  function _liveResolver(trackId) {
    return (s) => {
      if (!s || !s.popstatsLive) return null;
      const resp = s.popstatsLive.lastResponse;
      if (!resp) return null;
      // Run the response through the popstats adapter to get per-track shape
      const renderers = _renderers();
      if (!renderers) return null;
      const adapted = renderers.adaptPopstatsResponse(resp);
      if (adapted && adapted[trackId]) return adapted[trackId];
      // Hobs response is stashed under a different key in state if present
      if (s.popstatsLive.lastHobsResponse) {
        const hobsAdapted = renderers.adaptHobsResponse(s.popstatsLive.lastHobsResponse);
        if (hobsAdapted && hobsAdapted[trackId]) return hobsAdapted[trackId];
      }
      return null;
    };
  }

  // ===========================================================================
  // Wrapped getData — dispatches on mode chip
  // ===========================================================================

  function _wrappedGetData(trackId, originalGetData) {
    return (s) => {
      const renderers = _renderers();
      const mode = renderers ? renderers.getPopstatsMode() : 'static';
      const live = _liveResolver(trackId)(s);
      const stat = (typeof originalGetData === 'function') ? originalGetData(s) :
                    _staticResolver(trackId)(s);
      if (mode === 'live') {
        // Prefer live; if missing, fall back to static so the track still renders
        return live || stat;
      }
      if (mode === 'split') {
        // Compose: render both, with live as primary and static as ghost
        if (!live) return stat;
        if (!stat) return live;
        // Multi-line shape: append static as a dashed ghost series.
        // Single-line shape: convert to multi-line if both present.
        const liveAsML = _ensureMultilineShape(live, trackId);
        const statAsML = _ensureMultilineShape(stat, trackId);
        const ghostSeries = (statAsML.series || []).map(s2 => ({
          name: 'static · ' + s2.name,
          color: s2.color, dashed: true, alpha: 0.55,
          values: s2.values,
        }));
        return {
          mb: liveAsML.mb,
          series: liveAsML.series.concat(ghostSeries),
          refLine: liveAsML.refLine,
          yMin: liveAsML.yMin, yMax: liveAsML.yMax,
        };
      }
      // mode === 'static' (default)
      return stat;
    };
  }

  function _ensureMultilineShape(data, trackId) {
    if (!data) return { mb: [], series: [] };
    if (Array.isArray(data.series)) return data;
    // Single-line: wrap as one-series multi-line
    if (Array.isArray(data.values)) {
      return {
        mb: data.mb || [],
        series: [{ name: trackId, color: '#7b3294', values: data.values }],
        refLine: data.refLine,
        yMin: data.yMin, yMax: data.yMax,
      };
    }
    return { mb: [], series: [] };
  }

  // ===========================================================================
  // Public wrapping API
  // ===========================================================================

  function wrapTrackDef(t) {
    if (!t || typeof t !== 'object') return t;
    // Always-on tracks (ideogram, sim_collapse) keep their renderer; they
    // don't have getData, just custom renderers
    if (t.alwaysOn) return t;
    // Tracks already with getData get the wrapper applied OVER their static
    const original = (typeof t.getData === 'function') ? t.getData : null;
    return Object.assign({}, t, {
      getData: _wrappedGetData(t.id, original),
      // For multiline tracks, the renderer chip in the header gets a chance
      // (turn 4 already shipped this; nothing to add here)
      hasData: true,   // turn 6 onward: assume data may exist via live
    });
  }

  function wrapAllTrackDefs(list) {
    if (!Array.isArray(list)) return [];
    return list.map(wrapTrackDef);
  }

  function invalidatePage6() {
    if (typeof setTimeout === 'undefined') return;
    setTimeout(() => {
      if (typeof window === 'undefined') return;
      if (typeof window.renderPopstatsPage === 'function') {
        try { window.renderPopstatsPage(); } catch (_) {}
      }
    }, 0);
  }

  // ===========================================================================
  // Tooltip — pointer tracking over multi-line track canvases
  // ===========================================================================
  // The tooltip is a single floating div appended to <body>. attachTooltip(t,
  // canvas) installs pointermove on that canvas. On move we compute which
  // window index the cursor is over, then read the per-series values from
  // the cached data and render them.

  let _tooltipEl = null;
  let _tooltipDataCache = new WeakMap();   // canvas → {data, trackDef, padL, plotW, mbMin, mbMax}

  function _ensureTooltipEl() {
    if (typeof document === 'undefined') return null;
    if (_tooltipEl) return _tooltipEl;
    const el = document.createElement('div');
    el.id = 'popgen-popstats-tooltip';
    el.style.cssText = 'position:fixed; pointer-events:none; z-index:9500; ' +
      'background: var(--panel, #fbfcfd); border: 1px solid var(--rule, #c8cdd2); ' +
      'border-radius: 4px; padding: 4px 6px; font-family: var(--mono, monospace); ' +
      'font-size: 10px; color: var(--ink, #0e1116); display: none; ' +
      'box-shadow: 0 2px 8px rgba(0,0,0,0.15); max-width: 220px;';
    document.body.appendChild(el);
    _tooltipEl = el;
    return el;
  }

  function attachTooltip(trackDef, canvas, opts) {
    if (!canvas || !trackDef) return;
    opts = opts || {};
    const padL  = (typeof opts.padL  === 'number') ? opts.padL  : 80;
    const padR  = (typeof opts.padR  === 'number') ? opts.padR  : 60;
    if (canvas.__popgenTooltipAttached) return;
    canvas.__popgenTooltipAttached = true;

    canvas.addEventListener('pointermove', (e) => {
      const cache = _tooltipDataCache.get(canvas);
      if (!cache || !cache.data || !Array.isArray(cache.data.mb) ||
          cache.data.mb.length === 0) return;
      const rect = canvas.getBoundingClientRect();
      const x    = e.clientX - rect.left;
      const plotW = rect.width - padL - padR;
      if (x < padL || x > padL + plotW) { _hideTooltip(); return; }
      const { mbMin, mbMax, data, trackDef } = cache;
      const mbAtX = mbMin + ((x - padL) / plotW) * (mbMax - mbMin);
      // Find nearest mb index
      let nearest = 0, nearestD = Infinity;
      for (let i = 0; i < data.mb.length; i++) {
        const d = Math.abs(data.mb[i] - mbAtX);
        if (d < nearestD) { nearestD = d; nearest = i; }
      }
      _showTooltip(e.clientX, e.clientY, trackDef, data, nearest);
    });
    canvas.addEventListener('pointerleave', _hideTooltip);
  }

  function setTooltipDataFor(canvas, payload) {
    if (!canvas) return;
    _tooltipDataCache.set(canvas, payload);
  }

  function _showTooltip(x, y, trackDef, data, idx) {
    const el = _ensureTooltipEl();
    if (!el) return;
    const lines = [];
    const mb = data.mb[idx];
    lines.push(_fmtMb(mb));
    if (Array.isArray(data.series)) {
      for (const s of data.series) {
        const v = s.values && s.values[idx];
        const valStr = (v == null || isNaN(v)) ? '—' : _fmtNum(v);
        const dot = '<span style="display:inline-block;width:8px;height:8px;' +
                    'background:' + (s.color || '#999') + ';' +
                    'border-radius:1px;margin-right:4px;"></span>';
        lines.push(dot + s.name + ': <b>' + valStr + '</b>');
      }
    } else if (Array.isArray(data.values)) {
      const v = data.values[idx];
      lines.push((trackDef.label || trackDef.id || 'value') + ': <b>' +
                 (v == null || isNaN(v) ? '—' : _fmtNum(v)) + '</b>');
    }
    el.innerHTML = lines.join('<br>');
    el.style.display = 'block';
    // Position near cursor; clamp to viewport
    const tw = el.offsetWidth || 180;
    const th = el.offsetHeight || 60;
    const vw = (typeof window !== 'undefined') ? window.innerWidth  : 1280;
    const vh = (typeof window !== 'undefined') ? window.innerHeight : 800;
    let lx = x + 12;
    let ly = y + 12;
    if (lx + tw > vw - 8) lx = x - tw - 12;
    if (ly + th > vh - 8) ly = y - th - 12;
    el.style.left = lx + 'px';
    el.style.top  = ly + 'px';
  }

  function _hideTooltip() {
    if (_tooltipEl) _tooltipEl.style.display = 'none';
  }

  function _fmtMb(mb) {
    if (!isFinite(mb)) return '?';
    return mb.toFixed(3) + ' Mb';
  }
  function _fmtNum(n) {
    if (Math.abs(n) >= 100)    return n.toFixed(1);
    if (Math.abs(n) >= 1)      return n.toFixed(3);
    if (Math.abs(n) >= 0.001)  return n.toFixed(4);
    return n.toExponential(2);
  }

  // ===========================================================================
  // Slot-combine UI
  // ===========================================================================
  // Adds a small "+ combine" affordance to each empty slot in the dock's
  // Groups tab. Clicking opens a tiny inline picker:
  //   - operator: minus | union | intersect
  //   - operands: two slot-name dropdowns (or multi-select for union/intersect)
  //   - name field
  //   - confirm / cancel
  // The result calls engine.setSlot() with the right shape.

  let _slotComboInstalled = false;

  function installSlotCombineUI() {
    if (_slotComboInstalled) return;
    if (typeof document === 'undefined') return;
    const dockTabs = _dockTabs();
    if (!dockTabs) return;
    // Inject extra CSS for the combine picker
    if (!document.getElementById('popgen-slot-combine-css')) {
      const style = document.createElement('style');
      style.id = 'popgen-slot-combine-css';
      style.textContent = `
        .pg-slot-combine-btn {
          font-size: 10px; padding: 1px 6px; border-radius: 8px;
          background: var(--panel-2, #f0f1f3);
          border: 1px solid var(--rule, #c8cdd2);
          color: var(--ink-dim, #555e69); cursor: pointer;
          font-family: var(--mono, monospace);
        }
        .pg-slot-combine-btn:hover {
          background: var(--accent, #f5a524); color: #0e1116;
          border-color: var(--accent, #f5a524);
        }
        .pg-combine-picker {
          padding: 8px; background: var(--panel, #fbfcfd);
          border: 1px solid var(--accent, #f5a524); border-radius: 4px;
          display: flex; flex-direction: column; gap: 4px;
          font-family: var(--mono, monospace); font-size: 10px;
        }
        .pg-combine-picker select, .pg-combine-picker input {
          font-family: inherit; font-size: 10px;
          padding: 2px 4px;
          border: 1px solid var(--rule, #c8cdd2);
          background: var(--panel, #fbfcfd);
          color: var(--ink, #0e1116);
          border-radius: 3px;
        }
        .pg-combine-row {
          display: flex; gap: 4px; align-items: center;
        }
        .pg-combine-row label { min-width: 40px; color: var(--ink-dim, #555e69); }
      `;
      document.head.appendChild(style);
    }
    // Re-decorate Groups tab on every refresh
    const origMakeGroupsTab = dockTabs.makeGroupsTab;
    dockTabs.makeGroupsTab = function () {
      const tab = origMakeGroupsTab.apply(this, arguments);
      _decorateSlotsWithCombineButton(tab);
      // Wrap the original _refresh to also re-decorate
      const origRefresh = tab._refresh;
      tab._refresh = function () {
        if (typeof origRefresh === 'function') origRefresh();
        _decorateSlotsWithCombineButton(tab);
      };
      return tab;
    };
    // Re-install if a tab is already mounted
    dockTabs.refresh && dockTabs.refresh('groups');
    _slotComboInstalled = true;
  }

  function _decorateSlotsWithCombineButton(groupsTab) {
    if (!groupsTab) return;
    // Find the slots section — the third .pg-section
    const sections = (groupsTab.children || []).filter(c =>
      c.attributes && (c.attributes.class || '').includes('pg-section'));
    if (sections.length < 3) return;
    const slotSection = sections[2];   // dim, draft, slots, computeStrip — slots is index 2
    if (!slotSection) return;
    const slotList = slotSection.children[1];   // title, list
    if (!slotList) return;
    // Each empty slot row gets a "+ combine" button after the idx cell if not present
    for (const row of slotList.children) {
      const isEmpty = row.attributes && row.attributes['data-empty'] === 'true';
      if (!isEmpty) continue;
      // Already decorated?
      const hasBtn = row.children && row.children.find(c =>
        c.attributes && (c.attributes.class || '').includes('pg-slot-combine-btn'));
      if (hasBtn) continue;
      const btn = _makeCombineButton(row, groupsTab);
      // Replace the disabled placeholder name input with the combine button
      // (the row has 6 children from turn 5b: idx, name, source, n, color, del)
      if (row.children && row.children.length >= 2) {
        // Remove the disabled name input and replace with button container
        row.removeChild(row.children[1]);
        // Insert button at index 1
        if (typeof row.insertBefore === 'function' && row.children.length > 1) {
          // The row's appendChild is the only mutation method we have;
          // just append (will sit at end). Acceptable visual.
          row.appendChild(btn);
        } else {
          row.appendChild(btn);
        }
      } else {
        row.appendChild(btn);
      }
    }
  }

  function _makeCombineButton(slotRow, groupsTab) {
    const btn = document.createElement('button');
    btn.className = 'pg-slot-combine-btn';
    btn.textContent = '+ combine';
    btn.title = 'Combine existing slots with minus/union/intersect';
    btn.addEventListener('click', () => {
      const picker = _makeCombinePicker(slotRow, groupsTab);
      // Mount picker BELOW the slot row by appending to its parent at the
      // right position. The picker takes over until confirmed/cancelled.
      const parent = slotRow.parentNode;
      if (!parent) return;
      // Insert picker right after this row's parent index
      parent.appendChild(picker);
    });
    return btn;
  }

  function _makeCombinePicker(slotRow, groupsTab) {
    const eng = _engine();
    const wrap = document.createElement('div');
    wrap.className = 'pg-combine-picker';

    const slots = eng ? eng.getSlots().filter(s => s) : [];

    const opRow = document.createElement('div');
    opRow.className = 'pg-combine-row';
    opRow.appendChild(_label('op'));
    const opSel = document.createElement('select');
    for (const v of ['minus', 'union', 'intersect']) {
      const o = document.createElement('option'); o.value = v; o.textContent = v;
      opSel.appendChild(o);
    }
    opRow.appendChild(opSel);

    const memRow = document.createElement('div');
    memRow.className = 'pg-combine-row';
    memRow.appendChild(_label('args'));
    // For minus: two single-selects (left, right)
    // For union/intersect: a small list of dropdowns with "+ add"
    const memContainer = document.createElement('div');
    memContainer.style.cssText = 'display:flex; gap:3px; flex-wrap:wrap;';
    memRow.appendChild(memContainer);

    function renderMembers() {
      while (memContainer.firstChild) memContainer.removeChild(memContainer.firstChild);
      const op = opSel.value;
      if (op === 'minus') {
        memContainer.appendChild(_slotSelect(slots, 'left'));
        memContainer.appendChild(document.createTextNode(' − '));
        memContainer.appendChild(_slotSelect(slots, 'right'));
      } else {
        // Start with two members
        memContainer.appendChild(_slotSelect(slots, 'm0'));
        memContainer.appendChild(_slotSelect(slots, 'm1'));
        const addBtn = document.createElement('button');
        addBtn.className = 'pg-slot-combine-btn';
        addBtn.textContent = '+';
        addBtn.title = 'Add another operand';
        addBtn.addEventListener('click', () => {
          const idx = memContainer.querySelectorAll('select[data-slot-arg]').length;
          memContainer.insertBefore(_slotSelect(slots, 'm' + idx), addBtn);
        });
        memContainer.appendChild(addBtn);
      }
    }
    opSel.addEventListener('change', renderMembers);
    renderMembers();

    const nameRow = document.createElement('div');
    nameRow.className = 'pg-combine-row';
    nameRow.appendChild(_label('name'));
    const nameInp = document.createElement('input');
    nameInp.type = 'text';
    nameInp.placeholder = 'slot name';
    nameInp.value = _suggestName(slots);
    nameRow.appendChild(nameInp);

    const actionRow = document.createElement('div');
    actionRow.className = 'pg-combine-row';
    const okBtn = document.createElement('button');
    okBtn.className = 'pg-slot-combine-btn';
    okBtn.textContent = '✓ confirm';
    okBtn.style.background = 'var(--accent, #f5a524)';
    okBtn.style.color = '#0e1116';
    okBtn.addEventListener('click', () => {
      const op = opSel.value;
      const name = (nameInp.value || '').trim();
      if (!name || !/^[A-Za-z0-9_]+$/.test(name)) {
        alert('Slot name: [A-Za-z0-9_]+');
        return;
      }
      let ref;
      if (op === 'minus') {
        const left  = memContainer.querySelector('select[data-slot-arg="left"]').value;
        const right = memContainer.querySelector('select[data-slot-arg="right"]').value;
        if (!left || !right) { alert('Pick both sides for minus.'); return; }
        ref = { left, right };
      } else {
        const sels = memContainer.querySelectorAll('select[data-slot-arg]');
        const members = [];
        for (const s of sels) if (s.value) members.push(s.value);
        if (members.length < 2) {
          alert('Pick at least two operands.');
          return;
        }
        ref = { members };
      }
      const targetIdx = _firstEmptySlotIdx();
      if (targetIdx < 0) {
        alert('No empty slot available; delete one first.');
        return;
      }
      eng.setSlot(targetIdx, { name, source: op, ref, color: '#7a8398' });
      _close();
      const dt = _dockTabs();
      if (dt && dt.refresh) dt.refresh('groups');
    });
    const cancelBtn = document.createElement('button');
    cancelBtn.className = 'pg-slot-combine-btn';
    cancelBtn.textContent = '✕ cancel';
    cancelBtn.addEventListener('click', _close);
    actionRow.appendChild(okBtn);
    actionRow.appendChild(cancelBtn);

    wrap.appendChild(opRow);
    wrap.appendChild(memRow);
    wrap.appendChild(nameRow);
    wrap.appendChild(actionRow);

    function _close() {
      if (wrap.parentNode) wrap.parentNode.removeChild(wrap);
    }
    return wrap;
  }

  function _label(text) {
    const l = document.createElement('label');
    l.textContent = text;
    return l;
  }

  function _slotSelect(slots, role) {
    const sel = document.createElement('select');
    sel.setAttribute('data-slot-arg', role);
    const blank = document.createElement('option');
    blank.value = ''; blank.textContent = '(pick)';
    sel.appendChild(blank);
    for (const s of slots) {
      const o = document.createElement('option');
      o.value = s.name; o.textContent = s.name;
      sel.appendChild(o);
    }
    return sel;
  }

  function _firstEmptySlotIdx() {
    const eng = _engine();
    if (!eng) return -1;
    const slots = eng.getSlots();
    for (let i = 0; i < eng.MAX_SLOTS; i++) {
      if (!slots[i]) return i;
    }
    return -1;
  }

  function _suggestName(existing) {
    const used = new Set(existing.map(s => s.name));
    for (let i = 0; i < 26; i++) {
      const n = 'C' + String.fromCharCode(65 + i);   // CA, CB, CC...
      if (!used.has(n)) return n;
    }
    return 'combined_' + Date.now().toString(36).slice(-4);
  }

  // ===========================================================================
  // Public exports
  // ===========================================================================

  return {
    wrapTrackDef, wrapAllTrackDefs,
    invalidatePage6,
    attachTooltip, setTooltipDataFor,
    installSlotCombineUI,
    _internals: {
      _staticResolver, _liveResolver, _wrappedGetData,
      _ensureMultilineShape,
      _firstEmptySlotIdx, _suggestName,
      get tooltipEl() { return _tooltipEl; },
      get tooltipDataCache() { return _tooltipDataCache; },
    },
  };
}));
