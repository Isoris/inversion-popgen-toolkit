// =============================================================================
// atlas_dock_tabs.js — Turn 5b of chat A
// =============================================================================
//
// Builds the three tab UIs for the floating dock from turn 5a:
//
//   - Groups tab: dimension picker + draft-slot scratchpad + slot grid +
//                 Compute button with cursor-scope chip
//   - Lasso tab:  selection history table with pin / drop / load-as-set
//   - Snaps tab:  snapshot list with restore / pin / drop / inline note edit
//
// Wires to:
//   window.popgen      — engine (turn 2)
//   window.popgenLive  — request layer (turn 3)
//   window.popgenDock  — dock shell (turn 5a)
//
// API (window.popgenDockTabs):
//   .install()          — build all three tabs and inject into dock
//   .refreshAll()       — re-render every tab from current engine state
//   .refresh(tabId)     — refresh just one tab
//   .makeGroupsTab()    — DOM element for Groups tab (testable in isolation)
//   .makeLassoTab()
//   .makeSnapshotsTab()
//   .draftSlot()        — return the current draft-slot record
//   .compute()          — fire popstatsGroupwise for current slots+scope
//
// Persistence: draft slot is in localStorage under the popgen prefix.
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenDockTabs = factory();
    if (typeof document !== 'undefined') {
      const tryInstall = () => {
        if (root.popgenDock && root.popgen) root.popgenDockTabs.install();
        else setTimeout(tryInstall, 50);
      };
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', tryInstall);
      } else { tryInstall(); }
    }
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const LS_DRAFT = 'inversion_atlas.popgen.draftSlot';

  // ===========================================================================
  // Atlas plumbing accessors
  // ===========================================================================

  function _engine() {
    if (typeof window !== 'undefined' && window.popgen) return window.popgen;
    if (typeof globalThis !== 'undefined' && globalThis.popgen) return globalThis.popgen;
    return null;
  }
  function _live() {
    if (typeof window !== 'undefined' && window.popgenLive) return window.popgenLive;
    if (typeof globalThis !== 'undefined' && globalThis.popgenLive) return globalThis.popgenLive;
    return null;
  }
  function _dock() {
    if (typeof window !== 'undefined' && window.popgenDock) return window.popgenDock;
    if (typeof globalThis !== 'undefined' && globalThis.popgenDock) return globalThis.popgenDock;
    return null;
  }
  function _atlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }

  // ===========================================================================
  // CSS injection (extends the dock's stylesheet — same scoping pattern)
  // ===========================================================================

  const TABS_CSS = `
    [data-popgen-tab] { font-family: var(--mono, ui-monospace, monospace); font-size: 11px; }
    [data-popgen-tab="groups"] { display: flex; flex-direction: column; gap: 8px; }

    .pg-section { display: flex; flex-direction: column; gap: 4px; }
    .pg-section-title {
      font-size: 9px; font-weight: 600; letter-spacing: 0.06em;
      text-transform: uppercase; color: var(--ink-dim, #555e69);
      padding: 2px 0;
    }
    .pg-empty { color: var(--ink-dim, #555e69); font-style: italic;
                font-size: 10px; padding: 4px 0; }

    /* Dimension picker */
    .pg-dim-cat { background: var(--panel-2, #f0f1f3);
                  border: 1px solid var(--rule, #c8cdd2); border-radius: 4px;
                  margin-bottom: 4px; }
    .pg-dim-cat-header {
      display: flex; align-items: center; gap: 4px;
      padding: 3px 6px; cursor: pointer; user-select: none;
      font-size: 10px; color: var(--ink, #0e1116); }
    .pg-dim-cat-header:hover { background: var(--panel-3, #e6e8ec); }
    .pg-dim-cat-toggle { width: 10px; display: inline-block; color: var(--ink-dim, #555e69); }
    .pg-dim-cat-count { color: var(--ink-dim, #555e69); margin-left: auto;
                         font-size: 9px; }
    .pg-dim-cat-body { padding: 3px 6px; display: none; }
    .pg-dim-cat[data-open="true"] .pg-dim-cat-body { display: block; }
    .pg-dim-cat[data-open="true"] .pg-dim-cat-toggle::before { content: '▾'; }
    .pg-dim-cat[data-open="false"] .pg-dim-cat-toggle::before { content: '▸'; }

    .pg-dim-row { display: flex; align-items: center; gap: 4px;
                   padding: 2px 0; flex-wrap: wrap; }
    .pg-dim-name { font-size: 9px; color: var(--ink-dim, #555e69);
                    width: 60px; flex-shrink: 0;
                    overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .pg-dim-values { display: flex; flex-wrap: wrap; gap: 2px; flex: 1; }
    .pg-dim-chip {
      padding: 1px 5px; border-radius: 8px;
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      font-size: 9px; cursor: pointer; user-select: none;
      color: var(--ink, #0e1116);
    }
    .pg-dim-chip:hover { background: var(--panel-3, #e6e8ec); border-color: var(--accent, #f5a524); }
    .pg-dim-chip[data-on="true"] {
      background: var(--accent, #f5a524); color: #0e1116;
      border-color: var(--accent, #f5a524); font-weight: 600; }
    .pg-dim-bool {
      padding: 1px 5px; border-radius: 8px;
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      font-size: 9px; cursor: pointer; user-select: none;
    }
    .pg-dim-bool:hover { border-color: var(--accent, #f5a524); }

    /* Draft slot */
    .pg-draft {
      background: rgba(245, 165, 36, 0.08);
      border: 1px dashed var(--accent, #f5a524);
      border-radius: 4px;
      padding: 6px 8px;
      display: flex; flex-direction: column; gap: 4px;
    }
    .pg-draft-pills { display: flex; flex-wrap: wrap; gap: 3px; min-height: 18px; }
    .pg-draft-pill {
      display: inline-flex; align-items: center; gap: 3px;
      padding: 2px 6px;
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--accent, #f5a524);
      border-radius: 8px; font-size: 9px;
    }
    .pg-draft-pill button {
      border: 0; background: transparent; cursor: pointer;
      color: var(--ink-dim, #555e69);
      padding: 0 2px; font-size: 11px; line-height: 1;
    }
    .pg-draft-pill button:hover { color: var(--bad, #e0555c); }
    .pg-draft-actions { display: flex; gap: 4px; flex-wrap: wrap; align-items: center; }
    .pg-draft-count { font-size: 10px; color: var(--ink, #0e1116); flex: 1;
                       font-weight: 600; }
    .pg-btn {
      padding: 3px 8px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      border-radius: 4px; cursor: pointer;
      font-family: inherit; font-size: 10px;
      color: var(--ink, #0e1116);
    }
    .pg-btn:hover { background: var(--panel-3, #e6e8ec); }
    .pg-btn-primary {
      background: var(--accent, #f5a524); color: #0e1116;
      border-color: var(--accent, #f5a524); font-weight: 600;
    }
    .pg-btn-primary:hover { background: #d99119; }
    .pg-btn-danger:hover { background: var(--bad, #e0555c); color: white; }

    /* Slot grid */
    .pg-slot {
      display: grid;
      grid-template-columns: 18px 60px 1fr 30px 14px 14px;
      gap: 3px; align-items: center;
      padding: 3px 4px;
      border-bottom: 1px solid var(--rule, #c8cdd2);
      font-size: 10px;
    }
    .pg-slot[data-empty="true"] { opacity: 0.4; font-style: italic; }
    .pg-slot-idx { font-size: 9px; color: var(--ink-dim, #555e69); text-align: center; }
    .pg-slot-name { font-family: inherit; font-size: 10px;
                     border: 1px solid transparent; background: transparent;
                     padding: 1px 3px; color: var(--ink, #0e1116);
                     border-radius: 3px; min-width: 0; }
    .pg-slot-name:hover, .pg-slot-name:focus {
      border-color: var(--rule, #c8cdd2); background: var(--panel, #fbfcfd); }
    .pg-slot-source {
      font-size: 9px; color: var(--ink-dim, #555e69);
      overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }
    .pg-slot-n { font-size: 9px; color: var(--ink, #0e1116); text-align: right;
                  font-weight: 600; }
    .pg-slot-color {
      width: 12px; height: 12px; border-radius: 6px;
      border: 1px solid rgba(0,0,0,0.2); cursor: pointer;
    }
    .pg-slot-del { font-size: 12px; color: var(--ink-dim, #555e69);
                   background: transparent; border: 0; cursor: pointer; padding: 0; }
    .pg-slot-del:hover { color: var(--bad, #e0555c); }

    /* Compute strip */
    .pg-compute-strip {
      display: flex; gap: 6px; align-items: center;
      padding: 6px 0; border-top: 1px solid var(--rule, #c8cdd2);
    }
    .pg-scope-chip {
      display: inline-flex; padding: 1px;
      background: var(--panel-2, #f0f1f3);
      border-radius: 8px; font-size: 9px;
    }
    .pg-scope-chip button {
      padding: 1px 5px; border: 0; cursor: pointer;
      background: transparent; border-radius: 6px;
      font-family: inherit; font-size: inherit;
      color: var(--ink, #0e1116);
    }
    .pg-scope-chip button[data-on="true"] {
      background: var(--accent, #f5a524); font-weight: 600;
    }

    .pg-result { font-size: 9px; color: var(--ink-dim, #555e69); }
    .pg-result[data-state="error"] { color: var(--bad, #e0555c); }
    .pg-result[data-state="ok"]    { color: var(--good, #3cc08a); }

    /* Lasso + snapshots tables */
    [data-popgen-tab="lasso"], [data-popgen-tab="snapshots"] {
      display: flex; flex-direction: column; gap: 4px;
    }
    .pg-row {
      display: grid;
      grid-template-columns: 1fr 30px 14px 14px;
      gap: 3px; align-items: center;
      padding: 3px 4px;
      border-bottom: 1px solid var(--rule, #c8cdd2);
      font-size: 10px;
    }
    .pg-row[data-pinned="true"] { background: rgba(245,165,36,0.05); }
    .pg-row-meta { color: var(--ink, #0e1116);
                    overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .pg-row-meta .pg-sub { color: var(--ink-dim, #555e69); font-size: 9px; }
    .pg-row-n { font-size: 9px; text-align: right; color: var(--ink, #0e1116);
                font-weight: 600; }
    .pg-row-act {
      background: transparent; border: 0; cursor: pointer; padding: 0;
      font-size: 12px; color: var(--ink-dim, #555e69);
    }
    .pg-row-act:hover { color: var(--accent, #f5a524); }
  `;

  function _injectCSS() {
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-tabs-css')) return;
    const el = document.createElement('style');
    el.id = 'popgen-tabs-css';
    el.textContent = TABS_CSS;
    document.head.appendChild(el);
  }

  // ===========================================================================
  // Tiny DOM helpers
  // ===========================================================================

  function el(tag, props, children) {
    const e = document.createElement(tag);
    if (props) {
      for (const k of Object.keys(props)) {
        if (k === 'class')          e.className = props[k];
        else if (k === 'text')      e.textContent = props[k];
        else if (k === 'style')     e.style.cssText = props[k];
        else if (k === 'on') {
          for (const ev of Object.keys(props.on)) {
            e.addEventListener(ev, props.on[ev]);
          }
        }
        else if (k.startsWith('data-') || k === 'title' || k === 'placeholder' ||
                 k === 'value' || k === 'type')        e.setAttribute(k, props[k]);
        else                                            e[k] = props[k];
      }
    }
    if (children) for (const c of children) {
      if (c == null) continue;
      if (typeof c === 'string') e.appendChild(document.createTextNode(c));
      else                       e.appendChild(c);
    }
    return e;
  }

  function clearChildren(node) {
    while (node && node.firstChild) node.removeChild(node.firstChild);
  }

  // ===========================================================================
  // Draft slot — the scratchpad expression the user is composing right now
  // ===========================================================================
  // Stored as { predicates: [{dim, op, value(s)}, ...], joinOp: 'and'|'or' }.
  // The "Save / send to slot" actions translate this to an inline AND-of-eq
  // expression that goes into a slot record.

  function _readDraft() {
    if (typeof localStorage === 'undefined') return _newDraft();
    try {
      const raw = localStorage.getItem(LS_DRAFT);
      if (raw) {
        const d = JSON.parse(raw);
        if (d && Array.isArray(d.predicates)) return d;
      }
    } catch (_) {}
    return _newDraft();
  }
  function _writeDraft(d) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(LS_DRAFT, JSON.stringify(d)); } catch (_) {}
  }
  function _newDraft() { return { predicates: [], joinOp: 'and' }; }
  function _clearDraft() { _draftSlot = _newDraft(); _writeDraft(_draftSlot); }

  let _draftSlot = _readDraft();
  function draftSlot() { return _draftSlot; }

  // Convert the draft's predicates into an expression AST.
  function _draftToExpr() {
    const preds = _draftSlot.predicates;
    if (preds.length === 0) return { type: 'all' };
    const nodes = preds.map(p => {
      if (p.op === 'eq')      return { type: 'eq',      dim: p.dim, value: p.value };
      if (p.op === 'in')      return { type: 'in',      dim: p.dim, values: p.values };
      if (p.op === 'truthy')  return { type: 'truthy',  dim: p.dim };
      if (p.op === 'falsy')   return { type: 'falsy',   dim: p.dim };
      if (p.op === 'not_in')  return { type: 'not_in',  dim: p.dim, values: p.values };
      return { type: 'all' };
    });
    if (nodes.length === 1) return nodes[0];
    return { type: _draftSlot.joinOp || 'and', children: nodes };
  }

  function _addDraftPredicate(pred) {
    _draftSlot.predicates.push(pred);
    _writeDraft(_draftSlot);
  }
  function _removeDraftPredicateAt(i) {
    _draftSlot.predicates.splice(i, 1);
    _writeDraft(_draftSlot);
  }

  // ===========================================================================
  // Groups tab
  // ===========================================================================

  function makeGroupsTab() {
    const root = el('div', { 'data-popgen-tab': 'groups' });
    const dimSection    = el('div', { class: 'pg-section' });
    const dimTitle      = el('div', { class: 'pg-section-title', text: 'dimensions' });
    const dimList       = el('div');
    dimSection.appendChild(dimTitle);
    dimSection.appendChild(dimList);

    const draftSection  = el('div', { class: 'pg-section' });
    const draftTitle    = el('div', { class: 'pg-section-title', text: 'draft' });
    const draftBox      = el('div', { class: 'pg-draft' });
    const draftPills    = el('div', { class: 'pg-draft-pills' });
    const draftActions  = el('div', { class: 'pg-draft-actions' });
    draftBox.appendChild(draftPills);
    draftBox.appendChild(draftActions);
    draftSection.appendChild(draftTitle);
    draftSection.appendChild(draftBox);

    const slotSection   = el('div', { class: 'pg-section' });
    const slotTitle     = el('div', { class: 'pg-section-title', text: 'slots' });
    const slotList      = el('div');
    slotSection.appendChild(slotTitle);
    slotSection.appendChild(slotList);

    const computeStrip  = el('div', { class: 'pg-compute-strip' });
    const scopeChip     = el('span', { class: 'pg-scope-chip' });
    const computeBtn    = el('button', { class: 'pg-btn pg-btn-primary', text: 'Compute' });
    const resultLine    = el('span', { class: 'pg-result' });
    computeStrip.appendChild(scopeChip);
    computeStrip.appendChild(computeBtn);
    computeStrip.appendChild(resultLine);

    root.appendChild(dimSection);
    root.appendChild(draftSection);
    root.appendChild(slotSection);
    root.appendChild(computeStrip);

    // Refresh sub-functions ----------------------------------------------------

    function refreshDims() {
      const eng = _engine();
      clearChildren(dimList);
      if (!eng) {
        dimList.appendChild(el('div', { class: 'pg-empty', text: 'engine not loaded' }));
        return;
      }
      const dims = eng.listDimensions();
      const cats = _groupDimsByCategory(dims);
      for (const catName of Object.keys(cats)) {
        const cat = cats[catName];
        const wrap = el('div', { class: 'pg-dim-cat', 'data-open': cat.open ? 'true' : 'false' });
        const header = el('div', { class: 'pg-dim-cat-header' }, [
          el('span', { class: 'pg-dim-cat-toggle' }),
          el('span', { text: catName }),
          el('span', { class: 'pg-dim-cat-count', text: '(' + cat.dims.length + ')' }),
        ]);
        header.addEventListener('click', () => {
          const open = wrap.getAttribute('data-open') !== 'true';
          wrap.setAttribute('data-open', open ? 'true' : 'false');
          cat.open = open;
        });
        const body = el('div', { class: 'pg-dim-cat-body' });
        for (const d of cat.dims) {
          body.appendChild(_dimRow(d, refreshDraft));
        }
        wrap.appendChild(header);
        wrap.appendChild(body);
        dimList.appendChild(wrap);
      }
    }

    function refreshDraft() {
      const eng = _engine();
      if (!eng) return;
      clearChildren(draftPills);
      const preds = _draftSlot.predicates;
      if (preds.length === 0) {
        draftPills.appendChild(el('span', { class: 'pg-empty',
          text: 'click dimension values above to start' }));
      } else {
        for (let i = 0; i < preds.length; i++) {
          const p = preds[i];
          const label = _predicateLabel(p);
          const pill = el('span', { class: 'pg-draft-pill' }, [
            el('span', { text: label }),
            el('button', { text: '×', title: 'remove',
                            on: { click: () => { _removeDraftPredicateAt(i); refreshDraft(); } } }),
          ]);
          draftPills.appendChild(pill);
        }
      }
      // Update count and actions
      clearChildren(draftActions);
      let count = 0;
      try {
        const expr = _draftToExpr();
        count = eng.countMembers(expr);
      } catch (e) { count = 0; }
      draftActions.appendChild(el('span', { class: 'pg-draft-count',
        text: 'n = ' + count }));
      // Join op chip (and/or)
      if (preds.length >= 2) {
        const joinChip = el('span', { class: 'pg-scope-chip' });
        for (const op of ['and', 'or']) {
          const b = el('button', {
            text: op, 'data-on': (op === (_draftSlot.joinOp || 'and')) ? 'true' : 'false',
            on: { click: () => {
              _draftSlot.joinOp = op; _writeDraft(_draftSlot); refreshDraft();
            }},
          });
          joinChip.appendChild(b);
        }
        draftActions.appendChild(joinChip);
      }
      const sendBtn = el('button', { class: 'pg-btn',
        text: '→ slot',
        title: 'Save current draft into the next empty slot',
        on: { click: () => _draftToNextSlot(refreshSlots, refreshDraft) },
      });
      const saveExprBtn = el('button', { class: 'pg-btn',
        text: 'save expr',
        title: 'Save as a named expression (re-evaluates live)',
        on: { click: () => _saveDraftAsExpression(refreshDims) },
      });
      const saveSetBtn = el('button', { class: 'pg-btn',
        text: 'save set',
        title: 'Save as a frozen set of sample IDs',
        on: { click: () => _saveDraftAsSet(refreshDims) },
      });
      const clearBtn = el('button', { class: 'pg-btn pg-btn-danger',
        text: 'clear', on: { click: () => { _clearDraft(); refreshDraft(); } },
      });
      draftActions.appendChild(sendBtn);
      draftActions.appendChild(saveExprBtn);
      draftActions.appendChild(saveSetBtn);
      draftActions.appendChild(clearBtn);
    }

    function refreshSlots() {
      const eng = _engine();
      clearChildren(slotList);
      if (!eng) return;
      const slots = eng.getSlots();
      for (let i = 0; i < eng.MAX_SLOTS; i++) {
        slotList.appendChild(_slotRow(i, slots[i] || null, refreshSlots));
      }
    }

    function refreshComputeStrip() {
      clearChildren(scopeChip);
      const cur = _scopeMode();
      const scopes = ['1w', '5w', '10w', 'L2', 'candidate', 'chrom'];
      for (const s of scopes) {
        const b = el('button', {
          text: s, 'data-on': (s === cur) ? 'true' : 'false',
          on: { click: () => { _setScopeMode(s); refreshComputeStrip(); } },
        });
        scopeChip.appendChild(b);
      }
    }

    computeBtn.addEventListener('click', async () => {
      resultLine.setAttribute('data-state', '');
      resultLine.textContent = 'computing…';
      const r = await compute({ scope: _scopeMode() });
      if (r && r.ok) {
        resultLine.setAttribute('data-state', 'ok');
        resultLine.textContent = 'ok ' + r.cacheState +
          ' (' + (r.request_ms || 0) + 'ms)';
      } else {
        resultLine.setAttribute('data-state', 'error');
        resultLine.textContent = 'err: ' + (r && r.error ? r.error.slice(0, 40) : 'no live layer');
      }
    });

    // initial paint
    refreshDims();
    refreshDraft();
    refreshSlots();
    refreshComputeStrip();

    root._refresh = () => {
      refreshDims(); refreshDraft(); refreshSlots(); refreshComputeStrip();
    };
    return root;
  }

  // -- Dimension rendering ----------------------------------------------------

  function _groupDimsByCategory(dims) {
    const cats = {};
    function ensureCat(name, openByDefault) {
      if (!cats[name]) cats[name] = { dims: [], open: !!openByDefault };
      return cats[name];
    }
    for (const d of dims) {
      if (d.candidate) ensureCat('per candidate · ' + d.candidate, true).dims.push(d);
      else if (d.user)  ensureCat('saved sets', true).dims.push(d);
      else if (d.lasso) ensureCat('lasso history', false).dims.push(d);
      else              ensureCat('global', true).dims.push(d);
    }
    return cats;
  }

  function _dimRow(dim, onChange) {
    const eng = _engine();
    const row = el('div', { class: 'pg-dim-row' });
    row.appendChild(el('span', { class: 'pg-dim-name', text: _dimDisplayName(dim.name),
                                  title: dim.name }));
    const values = el('div', { class: 'pg-dim-values' });
    if (dim.kind === 'categorical') {
      // Enumerate distinct values from the dim's map
      const dm = eng._internals.lookupDimension(dim.name);
      const seen = new Set();
      for (const [, v] of dm) if (v != null) seen.add(String(v));
      const sorted = Array.from(seen).sort();
      for (const v of sorted) {
        const chip = el('span', { class: 'pg-dim-chip', text: v,
          on: { click: () => {
            _addDraftPredicate({ dim: dim.name, op: 'eq', value: v });
            if (onChange) onChange();
          }}
        });
        values.appendChild(chip);
      }
      if (sorted.length === 0) values.appendChild(el('span', { class: 'pg-empty', text: '(empty)' }));
    } else if (dim.kind === 'boolean') {
      const tBtn = el('span', { class: 'pg-dim-bool', text: '= true',
        on: { click: () => {
          _addDraftPredicate({ dim: dim.name, op: 'truthy' });
          if (onChange) onChange();
        }}
      });
      const fBtn = el('span', { class: 'pg-dim-bool', text: '= false',
        on: { click: () => {
          _addDraftPredicate({ dim: dim.name, op: 'falsy' });
          if (onChange) onChange();
        }}
      });
      values.appendChild(tBtn);
      values.appendChild(fBtn);
    } else if (dim.kind === 'integer' || dim.kind === 'numeric') {
      values.appendChild(el('span', { class: 'pg-empty', text: '(numeric — UI lands turn 6)' }));
    }
    row.appendChild(values);
    return row;
  }

  function _dimDisplayName(name) {
    const at = name.indexOf('@');
    if (at > 0) return name.slice(0, at);
    if (name.startsWith('saved_')) return name.slice(6);
    if (name.startsWith('lasso_')) return name.slice(6, 16) + '…';
    return name;
  }

  function _predicateLabel(p) {
    const dn = _dimDisplayName(p.dim);
    if (p.op === 'eq')     return dn + ' = ' + p.value;
    if (p.op === 'neq')    return dn + ' ≠ ' + p.value;
    if (p.op === 'in')     return dn + ' ∈ {' + (p.values || []).join(',') + '}';
    if (p.op === 'not_in') return dn + ' ∉ {' + (p.values || []).join(',') + '}';
    if (p.op === 'truthy') return dn;
    if (p.op === 'falsy')  return '¬' + dn;
    return p.op + ' ' + dn;
  }

  // -- Slot rendering ---------------------------------------------------------

  function _slotRow(idx, slot, onChange) {
    const eng = _engine();
    const empty = !slot;
    const row = el('div', { class: 'pg-slot', 'data-empty': empty ? 'true' : 'false' });
    row.appendChild(el('span', { class: 'pg-slot-idx', text: String(idx + 1) }));
    if (empty) {
      row.appendChild(el('input', { class: 'pg-slot-name', placeholder: '— empty —', disabled: true }));
      row.appendChild(el('span', { class: 'pg-slot-source', text: '' }));
      row.appendChild(el('span', { class: 'pg-slot-n', text: '' }));
      row.appendChild(el('span', { class: 'pg-slot-color', style: 'visibility:hidden' }));
      row.appendChild(el('button', { class: 'pg-slot-del', text: '', disabled: true,
                                     style: 'visibility:hidden' }));
      return row;
    }
    const nameInput = el('input', {
      class: 'pg-slot-name', value: slot.name, title: slot.name,
      on: {
        change: (e) => {
          const newName = e.target.value.trim();
          if (!/^[A-Za-z0-9_]+$/.test(newName)) {
            e.target.value = slot.name;
            return;
          }
          slot.name = newName;
          eng.setSlot(idx, slot);
          if (onChange) onChange();
        },
      },
    });
    row.appendChild(nameInput);
    row.appendChild(el('span', { class: 'pg-slot-source', text: _slotSourceLabel(slot),
                                  title: _slotSourceLabel(slot) }));
    let n = 0;
    try { n = eng.resolveSlotByIdx(idx).length; } catch (_) {}
    row.appendChild(el('span', { class: 'pg-slot-n', text: String(n) }));
    const swatch = el('span', { class: 'pg-slot-color',
                                  style: 'background:' + (slot.color || '#999') });
    swatch.addEventListener('click', () => {
      const c = prompt('hex color for slot ' + slot.name + ':', slot.color || '#');
      if (c) { slot.color = c; eng.setSlot(idx, slot); if (onChange) onChange(); }
    });
    row.appendChild(swatch);
    const delBtn = el('button', { class: 'pg-slot-del', text: '×', title: 'delete slot',
      on: { click: () => { eng.clearSlot(idx); if (onChange) onChange(); } } });
    row.appendChild(delBtn);
    return row;
  }

  function _slotSourceLabel(slot) {
    if (!slot) return '';
    if (slot.source === 'expression')  return 'expr: ' + slot.ref;
    if (slot.source === 'set')         return 'set: '  + slot.ref;
    if (slot.source === 'inline')      return 'inline: ' + _briefExpr(slot.ref);
    if (slot.source === 'complement_of') return '¬ ' + slot.ref;
    if (slot.source === 'minus' && slot.ref) return slot.ref.left + ' − ' + slot.ref.right;
    if (slot.source === 'union' && slot.ref) return '∪{' + (slot.ref.members || []).join(',') + '}';
    if (slot.source === 'intersect' && slot.ref) return '∩{' + (slot.ref.members || []).join(',') + '}';
    return slot.source || '';
  }

  function _briefExpr(expr) {
    if (!expr || typeof expr !== 'object') return '?';
    if (expr.type === 'all')  return '*';
    if (expr.type === 'none') return '∅';
    if (expr.type === 'eq')   return _dimDisplayName(expr.dim) + '=' + expr.value;
    if (expr.type === 'and' || expr.type === 'or')
      return '(' + (expr.children || []).map(_briefExpr).join(expr.type === 'and' ? '∧' : '∨') + ')';
    if (expr.type === 'not')   return '¬' + _briefExpr(expr.child);
    if (expr.type === 'minus') return _briefExpr(expr.left) + '−' + _briefExpr(expr.right);
    return expr.type;
  }

  // -- Draft actions ----------------------------------------------------------

  function _draftToNextSlot(refreshSlots, refreshDraft) {
    const eng = _engine();
    if (!eng) return;
    const slots = eng.getSlots();
    let target = -1;
    for (let i = 0; i < eng.MAX_SLOTS; i++) {
      if (!slots[i]) { target = i; break; }
    }
    if (target < 0) {
      alert('All ' + eng.MAX_SLOTS + ' slots are full. Delete one first.');
      return;
    }
    const expr = _draftToExpr();
    const proposedName = String.fromCharCode(65 + target);   // A, B, C...
    const name = (prompt('Slot name (A-Z, _, 0-9):', proposedName) || '').trim();
    if (!name) return;
    if (!/^[A-Za-z0-9_]+$/.test(name)) {
      alert('Slot names: [A-Za-z0-9_]+ only.');
      return;
    }
    eng.setSlot(target, { name, source: 'inline', ref: expr,
                          color: _autoColor(target) });
    _clearDraft();
    if (refreshDraft) refreshDraft();
    if (refreshSlots) refreshSlots();
  }

  function _saveDraftAsExpression(refreshDims) {
    const eng = _engine();
    if (!eng) return;
    const name = (prompt('Save expression as:') || '').trim();
    if (!name) return;
    if (!/^[A-Za-z0-9_]+$/.test(name)) {
      alert('Expression names: [A-Za-z0-9_]+ only.');
      return;
    }
    eng.saveExpression(name, _draftToExpr());
    if (refreshDims) refreshDims();
  }

  function _saveDraftAsSet(refreshDims) {
    const eng = _engine();
    if (!eng) return;
    const name = (prompt('Save frozen set as:') || '').trim();
    if (!name) return;
    if (!/^[A-Za-z0-9_]+$/.test(name)) {
      alert('Set names: [A-Za-z0-9_]+ only.');
      return;
    }
    eng.freezeSetFromExpression(name, _draftToExpr());
    if (refreshDims) refreshDims();
  }

  // -- Per-slot palette (deterministic by slot index) -------------------------
  const SLOT_PALETTE = [
    '#1f4e79', '#f5a524', '#3cc08a', '#e0555c', '#b07cf7',
    '#7ad3db', '#ff8c6e', '#9bbf7c', '#d28cb8', '#5a8fb3',
  ];
  function _autoColor(idx) { return SLOT_PALETTE[idx % SLOT_PALETTE.length]; }

  // -- Scope mode -------------------------------------------------------------
  const LS_SCOPE = 'inversion_atlas.popgen.scopeMode';
  function _scopeMode() {
    if (typeof localStorage === 'undefined') return 'candidate';
    try { return localStorage.getItem(LS_SCOPE) || 'candidate'; } catch (_) { return 'candidate'; }
  }
  function _setScopeMode(s) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(LS_SCOPE, s); } catch (_) {}
  }

  // -- Compute (fire popstats request) ----------------------------------------
  async function compute(opts) {
    opts = opts || {};
    const eng = _engine();
    const live = _live();
    if (!eng || !live) return { ok: false, error: 'engine or live layer missing' };
    const req = eng.buildPopstatsRequest({ scope: opts.scope || 'candidate' });
    if (!req.groups || Object.keys(req.groups).length === 0) {
      return { ok: false, error: 'no groups defined (fill at least one slot)' };
    }
    const env = await live.popstatsGroupwise(req, { debounce_ms: 0 });
    if (env.ok) {
      // Stash result for the renderers (turn 6 will hook getData() to this)
      const aS = _atlasState();
      if (aS) {
        aS.popstatsLive = aS.popstatsLive || {};
        aS.popstatsLive.lastResponse = env.payload;
        aS.popstatsLive.lastRequest  = req;
        aS.popstatsLive.lastEnvelope = env;
      }
    }
    return env;
  }

  // ===========================================================================
  // Lasso tab
  // ===========================================================================

  function makeLassoTab() {
    const root = el('div', { 'data-popgen-tab': 'lasso' });
    const list = el('div');
    root.appendChild(list);

    function refresh() {
      const eng = _engine();
      clearChildren(list);
      if (!eng) {
        list.appendChild(el('div', { class: 'pg-empty', text: 'engine not loaded' }));
        return;
      }
      const entries = eng.listSelections();
      if (entries.length === 0) {
        list.appendChild(el('div', { class: 'pg-empty',
          text: 'No lassos yet — drag-select on page 1 to record one.' }));
        return;
      }
      for (const e of entries) {
        list.appendChild(_lassoRow(e, refresh));
      }
    }

    refresh();
    root._refresh = refresh;
    return root;
  }

  function _lassoRow(e, onChange) {
    const eng = _engine();
    const row = el('div', { class: 'pg-row', 'data-pinned': e.pinned ? 'true' : 'false' });
    const meta = el('div', { class: 'pg-row-meta' });
    const subParts = [];
    if (e.source_page) subParts.push(e.source_page);
    if (e.cursor_chrom) subParts.push(e.cursor_chrom);
    if (e.cursor_candidate) subParts.push(e.cursor_candidate);
    meta.appendChild(el('span', { text: _shortTime(e.ts) + ' · ' + (e.kind || 'lasso') }));
    if (subParts.length > 0) {
      meta.appendChild(el('br'));
      meta.appendChild(el('span', { class: 'pg-sub', text: subParts.join(' · ') }));
    }
    row.appendChild(meta);
    row.appendChild(el('span', { class: 'pg-row-n', text: String(e.sample_ids.length) }));
    const pinBtn = el('button', { class: 'pg-row-act', text: e.pinned ? '★' : '☆',
      title: 'pin / unpin',
      on: { click: () => { eng.pinSelection(e.id, !e.pinned); if (onChange) onChange(); } },
    });
    row.appendChild(pinBtn);
    const loadBtn = el('button', { class: 'pg-row-act', text: '⤴', title: 'save as named set',
      on: { click: () => {
        const name = (prompt('Save lasso as set name:') || '').trim();
        if (!name) return;
        if (!/^[A-Za-z0-9_]+$/.test(name)) {
          alert('Set names: [A-Za-z0-9_]+ only.'); return;
        }
        try { eng.lassoToSet(e.id, name); }
        catch (err) { alert(String(err.message || err)); }
        // Re-render every tab so the new set shows up under "Saved" in groups
        const tabs = (typeof window !== 'undefined' && window.popgenDockTabs) ?
                     window.popgenDockTabs : null;
        if (tabs) tabs.refreshAll();
      }}
    });
    row.appendChild(loadBtn);
    const delBtn = el('button', { class: 'pg-row-act', text: '×', title: 'drop',
      on: { click: () => { eng.dropSelection(e.id); if (onChange) onChange(); }}});
    row.appendChild(delBtn);
    return row;
  }

  // ===========================================================================
  // Snapshots tab
  // ===========================================================================

  function makeSnapshotsTab() {
    const root = el('div', { 'data-popgen-tab': 'snapshots' });
    const list = el('div');
    const newBtn = el('button', { class: 'pg-btn', text: 'save snapshot now (Shift+S)',
      style: 'margin-bottom:6px;',
      on: { click: () => {
        const eng = _engine();
        if (!eng) return;
        const note = prompt('Snapshot note (optional):') || '';
        eng.saveSnapshot({ note });
        refresh();
      }}
    });
    root.appendChild(newBtn);
    root.appendChild(list);

    function refresh() {
      const eng = _engine();
      clearChildren(list);
      if (!eng) {
        list.appendChild(el('div', { class: 'pg-empty', text: 'engine not loaded' }));
        return;
      }
      const snaps = eng.listSnapshots();
      if (snaps.length === 0) {
        list.appendChild(el('div', { class: 'pg-empty',
          text: 'No snapshots yet — Shift+S to save.' }));
        return;
      }
      for (const s of snaps) list.appendChild(_snapshotRow(s, refresh));
    }

    refresh();
    root._refresh = refresh;
    return root;
  }

  function _snapshotRow(s, onChange) {
    const eng = _engine();
    const row = el('div', { class: 'pg-row', 'data-pinned': s.pinned ? 'true' : 'false' });
    const meta = el('div', { class: 'pg-row-meta' });
    meta.appendChild(el('span', { text: _shortTime(s.ts) +
                                  (s.note ? ' · ' + s.note : '') }));
    const subParts = [];
    if (s.chrom)             subParts.push(s.chrom);
    if (s.candidate && s.candidate.id) subParts.push(s.candidate.id);
    if (Array.isArray(s.slots)) {
      const nSlots = s.slots.filter(x => x).length;
      if (nSlots > 0) subParts.push(nSlots + ' slots');
    }
    if (subParts.length > 0) {
      meta.appendChild(el('br'));
      meta.appendChild(el('span', { class: 'pg-sub', text: subParts.join(' · ') }));
    }
    row.appendChild(meta);
    const slotsCount = Array.isArray(s.slots) ? s.slots.filter(x => x).length : 0;
    row.appendChild(el('span', { class: 'pg-row-n', text: String(slotsCount) }));
    const pinBtn = el('button', { class: 'pg-row-act', text: s.pinned ? '★' : '☆',
      title: 'pin / unpin',
      on: { click: () => { eng.pinSnapshot(s.id, !s.pinned); if (onChange) onChange(); } } });
    row.appendChild(pinBtn);
    const restoreBtn = el('button', { class: 'pg-row-act', text: '↻', title: 'restore',
      on: { click: () => {
        try { eng.restoreSnapshot(s.id); }
        catch (err) { alert(String(err.message || err)); return; }
        // Re-render every tab — the restored slots/expressions need to show
        const tabs = (typeof window !== 'undefined' && window.popgenDockTabs) ?
                     window.popgenDockTabs : null;
        if (tabs) tabs.refreshAll();
      }}
    });
    row.appendChild(restoreBtn);
    const delBtn = el('button', { class: 'pg-row-act', text: '×', title: 'drop',
      on: { click: () => { eng.dropSnapshot(s.id); if (onChange) onChange(); } } });
    row.appendChild(delBtn);
    return row;
  }

  function _shortTime(ts) {
    if (!ts) return '?';
    const d = new Date(ts);
    const ago = Date.now() - ts;
    if (ago < 60_000)       return 'just now';
    if (ago < 3_600_000)    return Math.floor(ago / 60_000) + ' min ago';
    if (ago < 86_400_000)   return Math.floor(ago / 3_600_000) + ' h ago';
    return d.toLocaleDateString() + ' ' + d.toLocaleTimeString().slice(0,5);
  }

  // ===========================================================================
  // Install — mount tabs into the dock + wire global callbacks
  // ===========================================================================

  let _installed = false;
  let _groupsTab = null, _lassoTab = null, _snapsTab = null;

  function install() {
    if (_installed) return;
    if (typeof document === 'undefined') return;
    const dock = _dock();
    if (!dock) return;
    _injectCSS();
    _groupsTab = makeGroupsTab();
    _lassoTab  = makeLassoTab();
    _snapsTab  = makeSnapshotsTab();
    dock.setTabContent('groups',    _groupsTab);
    dock.setTabContent('lasso',     _lassoTab);
    dock.setTabContent('snapshots', _snapsTab);
    // Wire Shift+S
    dock.onSaveSnapshot(() => {
      const eng = _engine();
      if (!eng) return;
      eng.saveSnapshot({ note: '' });
      refresh('snapshots');
    });
    // Also refresh on tab change
    dock.onTabChange((id) => refresh(id));
    _installed = true;
  }

  function refresh(tabId) {
    const tab = (tabId === 'groups') ? _groupsTab
              : (tabId === 'lasso')  ? _lassoTab
              : (tabId === 'snapshots') ? _snapsTab
              : null;
    if (tab && typeof tab._refresh === 'function') tab._refresh();
  }
  function refreshAll() {
    refresh('groups'); refresh('lasso'); refresh('snapshots');
  }

  return {
    install, refresh, refreshAll,
    makeGroupsTab, makeLassoTab, makeSnapshotsTab,
    draftSlot, compute,
    _internals: {
      _draftToExpr,
      get groupsTab() { return _groupsTab; },
      get lassoTab()  { return _lassoTab; },
      get snapsTab()  { return _snapsTab; },
      _addDraftPredicate, _clearDraft,
    },
  };
}));
