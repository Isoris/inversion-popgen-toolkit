// =============================================================================
// atlas_overview.js — Turn 9 of chat A
// =============================================================================
//
// Spreadsheet-style overview of all candidate inversions. Rows = candidates,
// columns = identity + quick-glance summary stats + manuscript flags.
//
// Click a row → focal candidate switches. Sort by any column. Filter by
// chrom, k, tier, confirmed. Multi-select for comparison view.
//
// Design choices:
//
//   - Data source: a `getCandidates()` callback (default: state.data.candidates
//     for the current chrom, fallback to state.atlas_catalogue if present).
//     Decoupled from atlas state so the overview works in narrow (single-chrom)
//     or wide (cross-chrom catalogue) modes.
//
//   - Live stats are LAZY: rows show "—" until row is selected or "Compute all"
//     is pressed. Compute reuses the turn-3 popgenLive layer; cold cache fills
//     in batches with progress indicator. Warm cache: instant.
//
//   - Per-row stats come from per-candidate compute (firing popstatsGroupwise
//     against the candidate's bp range). The path to a server-side batch
//     endpoint is documented in SPLICE_POINTS.md but not implemented here.
//
//   - Manuscript flags (confirmed, exclude, supplement) persist in a separate
//     localStorage namespace, so the overview can be used as a triage tool
//     during manuscript prep.
//
// API (window.popgenOverview):
//   .makeOverview(opts)              — DOM element with sortable, filterable
//                                       table. opts.getCandidates is the data
//                                       source (defaults to atlas state).
//   .computeRowStats(candidate)      — fire popstatsGroupwise for one row,
//                                       returns {fst_hom_hom, dxy_hom_hom,
//                                       het_hom1, het_het, het_hom2,
//                                       theta_pi_depression, q09b_verdict?}
//   .computeAllRows(candidates, opts) — batched compute over a list, with
//                                       progress callback
//   .setRowFlag(candId, flag, value) — persist a manuscript flag
//   .getRowFlag(candId, flag)
//   .exportTSV(candidates)           — synth a TSV string of the visible rows
//                                       (manuscript supplement)
//   .makeComparisonView(rows)        — small inline diff view for 2+ selected
//                                       rows, side-by-side
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenOverview = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const LS_PREFIX     = 'inversion_atlas.overview.';
  const LS_FLAGS      = LS_PREFIX + 'flags';      // candId → {confirmed, supplement, exclude, note}
  const LS_VIEW       = LS_PREFIX + 'view';       // sort/filter/visible-cols
  const LS_ROW_STATS  = LS_PREFIX + 'rowStats';   // candId → computed stats (cached)

  // Column registry. Each column has: id, label, kind, source, formatter,
  // sortable. kind in {'identity','live','flag'}. source ∈ {'cand','live','flag'}.
  const COLUMNS = [
    // Identity (always available, no compute)
    { id: 'chrom',       label: 'chrom',     kind: 'identity', source: 'cand', sortable: true,
      get: c => c.chrom || (c._chrom_inferred || ''),
      fmt: v => String(v || '—') },
    { id: 'start_mb',    label: 'start',     kind: 'identity', source: 'cand', sortable: true,
      get: c => _bpToMb(c.start_bp != null ? c.start_bp : c.start_mb * 1e6),
      fmt: v => v == null ? '—' : v.toFixed(2) },
    { id: 'end_mb',      label: 'end',       kind: 'identity', source: 'cand', sortable: true,
      get: c => _bpToMb(c.end_bp != null ? c.end_bp : c.end_mb * 1e6),
      fmt: v => v == null ? '—' : v.toFixed(2) },
    { id: 'width_mb',    label: 'width',     kind: 'identity', source: 'cand', sortable: true,
      get: c => {
        const a = (c.start_bp != null) ? c.start_bp : c.start_mb * 1e6;
        const b = (c.end_bp   != null) ? c.end_bp   : c.end_mb   * 1e6;
        return _bpToMb(b - a);
      },
      fmt: v => v == null ? '—' : v.toFixed(2) },
    { id: 'k',           label: 'k',         kind: 'identity', source: 'cand', sortable: true,
      get: c => c.k || c.n_clusters || null,
      fmt: v => v == null ? '—' : String(v) },
    { id: 'n_HOM1',      label: 'n HOM1',    kind: 'identity', source: 'cand', sortable: true,
      get: c => _countByRegime(c, 0),
      fmt: v => v == null ? '—' : String(v) },
    { id: 'n_HET',       label: 'n HET',     kind: 'identity', source: 'cand', sortable: true,
      get: c => _countByRegime(c, 1),
      fmt: v => v == null ? '—' : String(v) },
    { id: 'n_HOM2',      label: 'n HOM2',    kind: 'identity', source: 'cand', sortable: true,
      get: c => _countByRegime(c, 2),
      fmt: v => v == null ? '—' : String(v) },
    { id: 'tier',        label: 'tier',      kind: 'identity', source: 'cand', sortable: true,
      get: c => c.tier || (c.scoring && c.scoring.tier) || null,
      fmt: v => v == null ? '—' : String(v) },
    // Live (lazy compute)
    { id: 'fst_hh',      label: 'Fst H1·H2', kind: 'live', source: 'live', sortable: true,
      get: c => _liveStat(c, 'fst_hom_hom'),
      fmt: v => v == null ? '—' : v.toFixed(3) },
    { id: 'dxy_hh',      label: 'dXY H1·H2', kind: 'live', source: 'live', sortable: true,
      get: c => _liveStat(c, 'dxy_hom_hom'),
      fmt: v => v == null ? '—' : v.toFixed(4) },
    { id: 'theta_drop',  label: 'θπ drop',   kind: 'live', source: 'live', sortable: true,
      get: c => _liveStat(c, 'theta_pi_depression'),
      fmt: v => v == null ? '—' : v.toFixed(2) },
    { id: 'het_het',     label: 'het@HET',   kind: 'live', source: 'live', sortable: true,
      get: c => _liveStat(c, 'het_het'),
      fmt: v => v == null ? '—' : v.toFixed(3) },
    { id: 'q09b',        label: 'Q09b',      kind: 'live', source: 'live', sortable: true,
      get: c => _liveStat(c, 'q09b_verdict_short'),
      fmt: v => v == null ? '—' : v },
    // Flags (manual triage)
    { id: 'confirmed',   label: '✓',         kind: 'flag', source: 'flag', sortable: true,
      get: c => _flagFor(c, 'confirmed'),
      fmt: v => v ? '✓' : '' },
    { id: 'supplement',  label: 'sup',       kind: 'flag', source: 'flag', sortable: true,
      get: c => _flagFor(c, 'supplement'),
      fmt: v => v ? 'S' : '' },
    { id: 'exclude',     label: 'excl',      kind: 'flag', source: 'flag', sortable: true,
      get: c => _flagFor(c, 'exclude'),
      fmt: v => v ? '×' : '' },
    { id: 'note',        label: 'note',      kind: 'flag', source: 'flag', sortable: false,
      get: c => _flagFor(c, 'note'),
      fmt: v => v ? String(v).slice(0, 24) : '' },
  ];

  function _bpToMb(bp) {
    if (bp == null || !isFinite(bp)) return null;
    return bp / 1e6;
  }

  function _countByRegime(cand, regime) {
    if (!cand || !Array.isArray(cand.fish_calls)) return null;
    let n = 0;
    for (const fc of cand.fish_calls) {
      if (fc && fc.regime === regime) n++;
    }
    return n;
  }

  // ===========================================================================
  // Engine accessors
  // ===========================================================================

  function _atlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }
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
  function _renderers() {
    if (typeof window !== 'undefined' && window.popgenRenderers) return window.popgenRenderers;
    if (typeof globalThis !== 'undefined' && globalThis.popgenRenderers) return globalThis.popgenRenderers;
    return null;
  }
  function _turn7() {
    if (typeof window !== 'undefined' && window.popgenTurn7) return window.popgenTurn7;
    if (typeof globalThis !== 'undefined' && globalThis.popgenTurn7) return globalThis.popgenTurn7;
    return null;
  }

  // ===========================================================================
  // Row-stats cache (per cohort, per candidate)
  // ===========================================================================

  function _cohortKey() {
    const eng = _engine();
    if (eng && typeof eng.getCohortKey === 'function') {
      try { return eng.getCohortKey(); } catch (_) {}
    }
    return 'default';
  }

  function _readCache() {
    if (typeof localStorage === 'undefined') return {};
    try {
      const raw = localStorage.getItem(LS_ROW_STATS + '.' + _cohortKey());
      return raw ? JSON.parse(raw) : {};
    } catch (_) { return {}; }
  }
  function _writeCache(c) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(LS_ROW_STATS + '.' + _cohortKey(), JSON.stringify(c)); }
    catch (_) {}
  }

  function _liveStat(cand, key) {
    const c = _readCache();
    const row = c[cand.id || cand._id];
    if (!row) return null;
    return row[key];
  }

  function _setLiveStats(candId, stats) {
    const c = _readCache();
    c[candId] = Object.assign({}, c[candId], stats, { _ts: Date.now() });
    _writeCache(c);
  }

  function clearRowStatsCache() {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.removeItem(LS_ROW_STATS + '.' + _cohortKey()); } catch (_) {}
  }

  // ===========================================================================
  // Manuscript flags
  // ===========================================================================

  function _readFlags() {
    if (typeof localStorage === 'undefined') return {};
    try {
      const raw = localStorage.getItem(LS_FLAGS + '.' + _cohortKey());
      return raw ? JSON.parse(raw) : {};
    } catch (_) { return {}; }
  }
  function _writeFlags(f) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(LS_FLAGS + '.' + _cohortKey(), JSON.stringify(f)); }
    catch (_) {}
  }

  function _flagFor(cand, flagName) {
    const f = _readFlags();
    const row = f[cand.id || cand._id];
    if (!row) return null;
    return row[flagName];
  }

  function setRowFlag(candId, flag, value) {
    const f = _readFlags();
    f[candId] = Object.assign({}, f[candId], { [flag]: value });
    _writeFlags(f);
  }

  function getRowFlag(candId, flag) {
    const f = _readFlags();
    if (!f[candId]) return null;
    const v = f[candId][flag];
    return (v === undefined) ? null : v;
  }

  // ===========================================================================
  // Per-row compute — fires popstatsGroupwise for one candidate
  // ===========================================================================
  // Composes HOM1/HET/HOM2 inline groups from cand.fish_calls, fires Compute,
  // extracts shelf-summary stats from the response.

  async function computeRowStats(cand, opts) {
    opts = opts || {};
    const live = _live();
    if (!live) return { ok: false, error: 'popgenLive not available' };
    if (!cand || !Array.isArray(cand.fish_calls)) {
      return { ok: false, error: 'candidate has no fish_calls' };
    }

    // Build group → sample IDs from fish_calls
    const samples = (opts.atlasState && opts.atlasState.data && opts.atlasState.data.samples)
                    || (_atlasState() && _atlasState().data && _atlasState().data.samples)
                    || [];
    const groups = { HOM1: [], HET: [], HOM2: [] };
    for (const fc of cand.fish_calls) {
      if (!fc || fc.regime == null || fc.regime < 0) continue;
      const id = fc.sample_id || samples[fc.sample_idx] || ('idx_' + fc.sample_idx);
      if (fc.regime === 0) groups.HOM1.push(id);
      else if (fc.regime === 1) groups.HET.push(id);
      else if (fc.regime === 2) groups.HOM2.push(id);
    }
    // Drop groups smaller than min_n=3 (matches engine F default)
    if (groups.HOM1.length < 3 || groups.HOM2.length < 3) {
      return { ok: false, error: 'too few HOM1/HOM2 samples (<3 each)' };
    }
    if (groups.HET.length < 3) {
      delete groups.HET;
    }

    const start_bp = (cand.start_bp != null) ? cand.start_bp : Math.round(cand.start_mb * 1e6);
    const end_bp   = (cand.end_bp   != null) ? cand.end_bp   : Math.round(cand.end_mb   * 1e6);
    const chrom    = cand.chrom    || cand._chrom_inferred ||
                     (_atlasState() && _atlasState().data && _atlasState().data.chrom);

    if (!chrom || start_bp == null || end_bp == null) {
      return { ok: false, error: 'cannot resolve chrom/start/end' };
    }

    const req = {
      chrom,
      region: { start_bp, end_bp },
      groups,
      // Larger windows for the overview — we only want shelf summary stats,
      // not per-window resolution
      fixed_win: opts.fixed_win || [200_000, 100_000],
    };

    const env = await live.popstatsGroupwise(req, { debounce_ms: 0 });
    if (!env || !env.ok) {
      return { ok: false, error: (env && env.error) || 'compute failed' };
    }

    const stats = _summarizeShelfStats(env.payload, cand);
    _setLiveStats(cand.id, stats);
    return { ok: true, stats, cacheState: env.cacheState, request_ms: env.request_ms };
  }

  // Summarize a popstats response over the shelf into single-number per-row
  // stats. Means over windows in the response.
  function _summarizeShelfStats(payload, cand) {
    if (!payload || !Array.isArray(payload.windows) || payload.windows.length === 0) {
      return {};
    }
    const out = {};
    function meanOf(col) {
      let sum = 0, n = 0;
      for (const w of payload.windows) {
        const v = w[col];
        if (v == null || isNaN(v)) continue;
        sum += v; n++;
      }
      return n > 0 ? sum / n : null;
    }
    // Fst HOM1-HOM2 — column name from engine F is Fst_HOM1_HOM2
    out.fst_hom_hom = meanOf('Fst_HOM1_HOM2');
    // dXY likewise
    out.dxy_hom_hom = meanOf('dXY_HOM1_HOM2');
    // Per-class θπ — depression = mean(HOM1+HOM2)/2 vs HET
    const theta_h1 = meanOf('theta_pi_HOM1');
    const theta_he = meanOf('theta_pi_HET');
    const theta_h2 = meanOf('theta_pi_HOM2');
    if (theta_h1 != null && theta_h2 != null && theta_he != null && theta_he > 0) {
      const homMean = (theta_h1 + theta_h2) / 2;
      out.theta_pi_depression = homMean / theta_he;
    }
    // het@HET — proxy for HoverE (when hobs response not separately fired)
    const het_h1 = meanOf('het_HOM1');
    const het_he = meanOf('het_HET');
    const het_h2 = meanOf('het_HOM2');
    if (het_h1 != null) out.het_hom1 = het_h1;
    if (het_he != null) out.het_het  = het_he;
    if (het_h2 != null) out.het_hom2 = het_h2;
    // Q09b verdict left null — must be triggered separately because it
    // requires a dosage chunk fetch
    return out;
  }

  // ===========================================================================
  // Batch compute over many rows
  // ===========================================================================
  // Fires per-row computes serially with optional concurrency cap. Reports
  // progress via onProgress callback.

  async function computeAllRows(candidates, opts) {
    opts = opts || {};
    const concurrency = Math.max(1, Math.min(8, opts.concurrency || 3));
    const onProgress = opts.onProgress || (() => {});
    const onRowDone  = opts.onRowDone  || (() => {});
    const total = candidates.length;
    let done = 0, errors = 0;

    // Simple worker-pool pattern over a queue
    const queue = candidates.slice();
    async function worker() {
      while (queue.length > 0) {
        const cand = queue.shift();
        if (!cand) break;
        try {
          const result = await computeRowStats(cand, opts);
          if (!result.ok) errors++;
          onRowDone(cand, result);
        } catch (e) {
          errors++;
          onRowDone(cand, { ok: false, error: String(e) });
        }
        done++;
        onProgress({ done, total, errors });
        if (opts.cancelled && opts.cancelled()) break;
      }
    }
    const workers = [];
    for (let i = 0; i < concurrency; i++) workers.push(worker());
    await Promise.all(workers);
    return { done, errors, total };
  }

  // ===========================================================================
  // TSV export
  // ===========================================================================

  function exportTSV(candidates, opts) {
    opts = opts || {};
    const cols = (opts.columns && opts.columns.length > 0)
                  ? COLUMNS.filter(c => opts.columns.includes(c.id))
                  : COLUMNS.filter(c => c.id !== 'note');   // drop note from default export
    const lines = [];
    // Header
    lines.push(['candidate_id'].concat(cols.map(c => c.id)).join('\t'));
    for (const cand of candidates) {
      const row = [cand.id || ''];
      for (const c of cols) {
        const v = c.get(cand);
        if (v == null) row.push('');
        else if (typeof v === 'number') row.push(_fmtTSVNum(v));
        else row.push(String(v).replace(/\t/g, ' '));
      }
      lines.push(row.join('\t'));
    }
    return lines.join('\n');
  }
  function _fmtTSVNum(n) {
    if (Math.abs(n) >= 100)    return n.toFixed(2);
    if (Math.abs(n) >= 1)      return n.toFixed(4);
    if (Math.abs(n) >= 0.001)  return n.toFixed(5);
    return n.toExponential(3);
  }

  // ===========================================================================
  // Sorting + filtering
  // ===========================================================================

  function _sortRows(rows, sortKey, sortDir) {
    if (!sortKey) return rows;
    const col = COLUMNS.find(c => c.id === sortKey);
    if (!col) return rows;
    const dir = (sortDir === 'desc') ? -1 : 1;
    return rows.slice().sort((a, b) => {
      const va = col.get(a);
      const vb = col.get(b);
      // null/undefined sorts to end always
      if (va == null && vb == null) return 0;
      if (va == null) return 1;
      if (vb == null) return -1;
      if (typeof va === 'number' && typeof vb === 'number') {
        return (va - vb) * dir;
      }
      return String(va).localeCompare(String(vb)) * dir;
    });
  }

  function _filterRows(rows, filters) {
    if (!filters) return rows;
    return rows.filter(c => {
      if (filters.chrom && c.chrom !== filters.chrom) return false;
      if (filters.k != null && (c.k || 0) !== filters.k) return false;
      if (filters.confirmed === true  && !_flagFor(c, 'confirmed')) return false;
      if (filters.confirmed === false &&  _flagFor(c, 'confirmed')) return false;
      if (filters.exclude_hidden &&     _flagFor(c, 'exclude'))    return false;
      if (filters.tier && c.tier !== filters.tier) return false;
      if (filters.text) {
        const t = String(filters.text).toLowerCase();
        const blob = (c.id + ' ' + (c.chrom || '')).toLowerCase();
        if (!blob.includes(t)) return false;
      }
      return true;
    });
  }

  // ===========================================================================
  // DOM builder
  // ===========================================================================

  let _injectedCSS = false;

  function _injectCSS() {
    if (_injectedCSS) return;
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-overview-css')) {
      _injectedCSS = true; return;
    }
    const style = document.createElement('style');
    style.id = 'popgen-overview-css';
    style.textContent = OVERVIEW_CSS;
    document.head.appendChild(style);
    _injectedCSS = true;
  }

  function makeOverview(opts) {
    opts = opts || {};
    if (typeof document === 'undefined') return null;
    _injectCSS();

    const root = document.createElement('div');
    root.setAttribute('data-popgen-overview', 'root');

    // Toolbar
    const toolbar = document.createElement('div');
    toolbar.setAttribute('data-popgen-overview', 'toolbar');
    const title = document.createElement('span');
    title.setAttribute('data-popgen-overview', 'title');
    title.textContent = 'CANDIDATE OVERVIEW';
    toolbar.appendChild(title);

    const filterInp = document.createElement('input');
    filterInp.setAttribute('data-popgen-overview', 'filter');
    filterInp.placeholder = 'filter candidates…';
    filterInp.type = 'text';
    toolbar.appendChild(filterInp);

    // Chrom dropdown
    const chromSel = document.createElement('select');
    chromSel.setAttribute('data-popgen-overview', 'chrom-filter');
    toolbar.appendChild(chromSel);

    // Hide-excluded toggle
    const hideExcl = document.createElement('button');
    hideExcl.setAttribute('data-popgen-overview', 'hide-excluded');
    hideExcl.setAttribute('data-on', 'true');
    hideExcl.title = 'Hide rows flagged as excluded';
    hideExcl.textContent = 'hide ×';
    toolbar.appendChild(hideExcl);

    // Compute all
    const computeAllBtn = document.createElement('button');
    computeAllBtn.setAttribute('data-popgen-overview', 'compute-all');
    computeAllBtn.textContent = '▶ compute all visible';
    computeAllBtn.title = 'Run popstats compute for every visible row that doesn\'t already have stats';
    toolbar.appendChild(computeAllBtn);

    // Export TSV
    const exportBtn = document.createElement('button');
    exportBtn.setAttribute('data-popgen-overview', 'export');
    exportBtn.textContent = '↓ TSV';
    exportBtn.title = 'Download visible rows as TSV (ready for manuscript Supp Table)';
    toolbar.appendChild(exportBtn);

    // Progress
    const progress = document.createElement('span');
    progress.setAttribute('data-popgen-overview', 'progress');
    toolbar.appendChild(progress);

    root.appendChild(toolbar);

    // Table container
    const tableWrap = document.createElement('div');
    tableWrap.setAttribute('data-popgen-overview', 'table-wrap');
    const table = document.createElement('table');
    table.setAttribute('data-popgen-overview', 'table');
    const thead = document.createElement('thead');
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    table.appendChild(tbody);
    tableWrap.appendChild(table);
    root.appendChild(tableWrap);

    // Status line at bottom
    const status = document.createElement('div');
    status.setAttribute('data-popgen-overview', 'status');
    root.appendChild(status);

    // ---- State ----
    const view = _readView();
    let _candidates = [];
    let _selectedId = null;
    let _multiSelected = new Set();

    function _readView() {
      if (typeof localStorage === 'undefined') return _defaultView();
      try {
        const raw = localStorage.getItem(LS_VIEW + '.' + _cohortKey());
        const v = raw ? JSON.parse(raw) : {};
        return Object.assign(_defaultView(), v);
      } catch (_) { return _defaultView(); }
    }
    function _defaultView() {
      return { sortKey: 'chrom', sortDir: 'asc',
               filters: { exclude_hidden: true, text: '', chrom: '' },
               visibleCols: COLUMNS.map(c => c.id) };
    }
    function _saveView() {
      if (typeof localStorage === 'undefined') return;
      try { localStorage.setItem(LS_VIEW + '.' + _cohortKey(), JSON.stringify(view)); }
      catch (_) {}
    }

    // ---- Render functions ----

    function refreshHeader() {
      while (thead.firstChild) thead.removeChild(thead.firstChild);
      const tr = document.createElement('tr');
      // First column = row marker / select indicator
      const thMark = document.createElement('th');
      thMark.setAttribute('data-popgen-overview', 'mark');
      thMark.textContent = '·';
      tr.appendChild(thMark);
      // Identity column for the candidate id
      const thId = document.createElement('th');
      thId.setAttribute('data-popgen-overview', 'th');
      thId.setAttribute('data-col', 'id');
      thId.textContent = 'candidate';
      thId.addEventListener('click', () => _setSort('id'));
      tr.appendChild(thId);
      // Each registered column
      for (const c of COLUMNS) {
        if (!view.visibleCols.includes(c.id)) continue;
        const th = document.createElement('th');
        th.setAttribute('data-popgen-overview', 'th');
        th.setAttribute('data-col', c.id);
        th.setAttribute('data-kind', c.kind);
        th.textContent = c.label;
        th.title = c.label + ' (' + c.kind + ')';
        if (view.sortKey === c.id) {
          th.setAttribute('data-sort', view.sortDir);
        }
        if (c.sortable) {
          th.style.cursor = 'pointer';
          th.addEventListener('click', () => _setSort(c.id));
        }
        tr.appendChild(th);
      }
      thead.appendChild(tr);
    }

    function _setSort(key) {
      if (key === 'id') {
        if (view.sortKey === 'id') view.sortDir = (view.sortDir === 'asc') ? 'desc' : 'asc';
        else { view.sortKey = 'id'; view.sortDir = 'asc'; }
      } else {
        const col = COLUMNS.find(c => c.id === key);
        if (!col || !col.sortable) return;
        if (view.sortKey === key) view.sortDir = (view.sortDir === 'asc') ? 'desc' : 'asc';
        else { view.sortKey = key; view.sortDir = 'asc'; }
      }
      _saveView();
      refreshAll();
    }

    function refreshBody() {
      while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
      let rows = _filterRows(_candidates, view.filters);
      if (view.sortKey === 'id') {
        rows = rows.slice().sort((a, b) => {
          const dir = (view.sortDir === 'desc') ? -1 : 1;
          return String(a.id || '').localeCompare(String(b.id || '')) * dir;
        });
      } else {
        rows = _sortRows(rows, view.sortKey, view.sortDir);
      }
      for (const cand of rows) {
        tbody.appendChild(_makeRow(cand));
      }
      status.textContent = rows.length + ' visible / ' + _candidates.length + ' total';
    }

    function _makeRow(cand) {
      const tr = document.createElement('tr');
      tr.setAttribute('data-popgen-overview', 'row');
      tr.setAttribute('data-cand-id', cand.id || '');
      const isFocal = (_atlasState() && _atlasState().candidate &&
                        _atlasState().candidate.id === cand.id);
      const isSelected = (_selectedId === cand.id);
      const isMulti    = _multiSelected.has(cand.id);
      if (isFocal)    tr.setAttribute('data-focal', 'true');
      if (isSelected) tr.setAttribute('data-selected', 'true');
      if (isMulti)    tr.setAttribute('data-multi', 'true');
      if (_flagFor(cand, 'exclude'))   tr.setAttribute('data-excluded', 'true');
      if (_flagFor(cand, 'confirmed')) tr.setAttribute('data-confirmed', 'true');

      // Mark cell — clickable, toggles multi-select on shift-click
      const tdMark = document.createElement('td');
      tdMark.setAttribute('data-popgen-overview', 'mark');
      tdMark.textContent = isFocal ? '◉' : (isMulti ? '⊕' : (isSelected ? '◎' : '·'));
      tdMark.addEventListener('click', (e) => {
        if (e.shiftKey) {
          if (_multiSelected.has(cand.id)) _multiSelected.delete(cand.id);
          else                              _multiSelected.add(cand.id);
        } else {
          _selectedId = cand.id;
          if (opts.onRowSelect) opts.onRowSelect(cand);
        }
        refreshBody();
      });
      tr.appendChild(tdMark);

      // Candidate id column
      const tdId = document.createElement('td');
      tdId.setAttribute('data-popgen-overview', 'cand-id');
      tdId.textContent = cand.id || '?';
      tdId.title = cand.id;
      tdId.addEventListener('click', () => {
        if (opts.onRowFocus) opts.onRowFocus(cand);
      });
      tdId.style.cursor = 'pointer';
      tr.appendChild(tdId);

      for (const c of COLUMNS) {
        if (!view.visibleCols.includes(c.id)) continue;
        const td = document.createElement('td');
        td.setAttribute('data-popgen-overview', 'td');
        td.setAttribute('data-col', c.id);
        td.setAttribute('data-kind', c.kind);
        const v = c.get(cand);
        td.textContent = c.fmt(v);
        if (c.kind === 'flag') {
          td.style.cursor = 'pointer';
          td.addEventListener('click', () => {
            if (c.id === 'note') {
              const cur = _flagFor(cand, 'note') || '';
              const next = (typeof prompt === 'function')
                ? prompt('note for ' + cand.id, cur) : cur;
              if (next !== null && next !== undefined) {
                setRowFlag(cand.id, 'note', next.trim());
                refreshBody();
              }
            } else {
              setRowFlag(cand.id, c.id, !_flagFor(cand, c.id));
              refreshBody();
            }
          });
        }
        if (c.kind === 'live' && v == null) {
          td.style.opacity = '0.5';
        }
        tr.appendChild(td);
      }
      return tr;
    }

    function refreshChromFilter() {
      while (chromSel.firstChild) chromSel.removeChild(chromSel.firstChild);
      const allOpt = document.createElement('option');
      allOpt.value = ''; allOpt.textContent = '— all chroms —';
      chromSel.appendChild(allOpt);
      const seen = new Set();
      for (const c of _candidates) {
        const ch = c.chrom || c._chrom_inferred;
        if (ch && !seen.has(ch)) { seen.add(ch); }
      }
      const sorted = Array.from(seen).sort();
      for (const ch of sorted) {
        const o = document.createElement('option');
        o.value = ch; o.textContent = ch;
        chromSel.appendChild(o);
      }
      chromSel.value = view.filters.chrom || '';
    }

    function refreshAll() {
      refreshHeader();
      refreshChromFilter();
      refreshBody();
    }

    // Public refresh hook so the atlas can poke when candidates change
    root._refresh = refreshAll;

    // ---- Toolbar wiring ----

    filterInp.value = view.filters.text || '';
    let filterDebounce = null;
    filterInp.addEventListener('input', () => {
      clearTimeout(filterDebounce);
      filterDebounce = setTimeout(() => {
        view.filters.text = filterInp.value;
        _saveView();
        refreshBody();
      }, 150);
    });

    chromSel.addEventListener('change', () => {
      view.filters.chrom = chromSel.value;
      _saveView();
      refreshBody();
    });

    hideExcl.addEventListener('click', () => {
      const on = hideExcl.getAttribute('data-on') === 'true';
      hideExcl.setAttribute('data-on', on ? 'false' : 'true');
      view.filters.exclude_hidden = !on;
      _saveView();
      refreshBody();
    });

    computeAllBtn.addEventListener('click', async () => {
      const visibleRows = _filterRows(_candidates, view.filters);
      const todo = visibleRows.filter(c => _liveStat(c, 'fst_hom_hom') == null);
      if (todo.length === 0) {
        status.textContent = 'all visible rows already computed';
        return;
      }
      progress.textContent = '0 / ' + todo.length;
      computeAllBtn.disabled = true;
      const result = await computeAllRows(todo, {
        concurrency: 3,
        onProgress: ({ done, total, errors }) => {
          progress.textContent = done + ' / ' + total +
            (errors > 0 ? ' (' + errors + ' err)' : '');
        },
        onRowDone: () => refreshBody(),
      });
      computeAllBtn.disabled = false;
      progress.textContent = 'done: ' + result.done + ' / ' + result.total +
                             (result.errors > 0 ? ' · ' + result.errors + ' errors' : '');
      refreshBody();
    });

    exportBtn.addEventListener('click', () => {
      const rows = _filterRows(_candidates, view.filters);
      const tsv = exportTSV(rows);
      _downloadString(tsv,
        'overview_' + (_atlasState() && _atlasState().data && _atlasState().data.cohort_key
                        ? _atlasState().data.cohort_key : 'cohort') + '.tsv',
        'text/tab-separated-values');
    });

    // ---- Initial data load ----

    function reloadCandidates() {
      const fn = opts.getCandidates || _defaultGetCandidates;
      const list = fn() || [];
      _candidates = list.map(c => Object.assign({}, c, {
        // Inject _chrom_inferred so rows from the current chrom always have a chrom value
        _chrom_inferred: c.chrom || (
          _atlasState() && _atlasState().data && _atlasState().data.chrom
        ),
      }));
      refreshAll();
    }
    root._reload = reloadCandidates;

    reloadCandidates();
    return root;
  }

  function _defaultGetCandidates() {
    const s = _atlasState();
    if (!s || !s.data) return [];
    if (Array.isArray(s.atlas_catalogue)) return s.atlas_catalogue;
    if (Array.isArray(s.data.candidates)) return s.data.candidates;
    return [];
  }

  function _downloadString(content, filename, mime) {
    if (typeof document === 'undefined') return;
    if (typeof Blob === 'undefined' || typeof URL === 'undefined') return;
    try {
      const blob = new Blob([content], { type: mime || 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = filename;
      document.body.appendChild(a);
      a.click();
      setTimeout(() => {
        URL.revokeObjectURL(url);
        if (a.parentNode) a.parentNode.removeChild(a);
      }, 200);
    } catch (e) {
      console.warn('download failed', e);
    }
  }

  // ===========================================================================
  // Comparison view — small inline diff for 2+ selected rows
  // ===========================================================================
  // Used by the floating dock or page 6: when the user shift-clicks 2+ rows
  // in the overview, the comparison view shows them side-by-side with a
  // per-stat delta column.

  function makeComparisonView(rows) {
    if (!Array.isArray(rows) || rows.length === 0) return null;
    if (typeof document === 'undefined') return null;
    _injectCSS();
    const wrap = document.createElement('div');
    wrap.setAttribute('data-popgen-overview', 'compare');
    const title = document.createElement('div');
    title.setAttribute('data-popgen-overview', 'compare-title');
    title.textContent = 'comparing ' + rows.length + ' candidates';
    wrap.appendChild(title);
    const tbl = document.createElement('table');
    tbl.setAttribute('data-popgen-overview', 'compare-table');
    // Header: stat name | each row's value | delta (only if 2 rows)
    const thead = document.createElement('thead');
    const trh = document.createElement('tr');
    trh.appendChild(_th('stat'));
    for (const r of rows) {
      trh.appendChild(_th(r.id || '?'));
    }
    if (rows.length === 2) trh.appendChild(_th('Δ'));
    thead.appendChild(trh);
    tbl.appendChild(thead);
    const tbody = document.createElement('tbody');
    for (const c of COLUMNS) {
      if (c.kind === 'flag') continue;
      const tr = document.createElement('tr');
      tr.appendChild(_td(c.label));
      const vals = rows.map(r => c.get(r));
      for (const v of vals) tr.appendChild(_td(c.fmt(v)));
      if (rows.length === 2 && typeof vals[0] === 'number' && typeof vals[1] === 'number') {
        const d = vals[1] - vals[0];
        const dt = _td(d.toFixed(4));
        dt.setAttribute('data-popgen-overview', 'delta');
        if (d > 0) dt.setAttribute('data-sign', 'pos');
        if (d < 0) dt.setAttribute('data-sign', 'neg');
        tr.appendChild(dt);
      } else if (rows.length === 2) {
        tr.appendChild(_td('—'));
      }
      tbody.appendChild(tr);
    }
    tbl.appendChild(tbody);
    wrap.appendChild(tbl);
    return wrap;
  }
  function _th(t) {
    const e = document.createElement('th');
    e.textContent = t;
    return e;
  }
  function _td(t) {
    const e = document.createElement('td');
    e.textContent = t;
    return e;
  }

  // ===========================================================================
  // CSS
  // ===========================================================================

  const OVERVIEW_CSS = `
    [data-popgen-overview="root"] {
      display: flex; flex-direction: column;
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      font-family: var(--mono, ui-monospace, monospace);
      font-size: 11px; color: var(--ink, #0e1116);
      max-height: 80vh;
    }
    [data-popgen-overview="toolbar"] {
      display: flex; gap: 8px; align-items: center;
      padding: 6px 10px;
      background: var(--panel-2, #f0f1f3);
      border-bottom: 1px solid var(--rule, #c8cdd2);
      flex-shrink: 0;
    }
    [data-popgen-overview="title"] {
      font-weight: 700; letter-spacing: 0.06em;
      color: var(--ink-dim, #555e69); font-size: 10px;
    }
    [data-popgen-overview="filter"] {
      flex: 1; max-width: 220px;
      padding: 3px 6px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-radius: 3px;
      font-family: inherit; font-size: 11px;
    }
    [data-popgen-overview="chrom-filter"] {
      padding: 3px 6px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-radius: 3px;
      font-family: inherit; font-size: 11px;
    }
    [data-popgen-overview="hide-excluded"],
    [data-popgen-overview="compute-all"],
    [data-popgen-overview="export"] {
      padding: 3px 8px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-radius: 3px; cursor: pointer;
      font-family: inherit; font-size: 11px;
    }
    [data-popgen-overview="hide-excluded"][data-on="true"] {
      background: var(--accent, #f5a524); color: #0e1116;
      border-color: var(--accent, #f5a524);
    }
    [data-popgen-overview="compute-all"]:hover,
    [data-popgen-overview="export"]:hover {
      background: var(--panel-3, #e6e8ec);
    }
    [data-popgen-overview="compute-all"]:disabled {
      opacity: 0.5; cursor: not-allowed;
    }
    [data-popgen-overview="progress"] {
      font-size: 10px; color: var(--ink-dim, #555e69);
      margin-left: auto;
    }
    [data-popgen-overview="table-wrap"] {
      flex: 1; overflow: auto;
    }
    [data-popgen-overview="table"] {
      width: 100%; border-collapse: collapse;
      font-size: 10px; font-family: var(--mono, monospace);
    }
    [data-popgen-overview="table"] thead {
      position: sticky; top: 0;
      background: var(--panel-2, #f0f1f3);
      border-bottom: 1px solid var(--rule, #c8cdd2);
    }
    [data-popgen-overview="th"] {
      padding: 5px 6px; text-align: left;
      font-weight: 600; user-select: none;
      border-right: 1px solid var(--rule, #c8cdd2);
      white-space: nowrap;
    }
    [data-popgen-overview="th"][data-kind="live"] { color: var(--accent, #f5a524); }
    [data-popgen-overview="th"][data-kind="flag"] { color: var(--ink-dim, #555e69); }
    [data-popgen-overview="th"][data-sort="asc"]::after {
      content: ' ▲'; font-size: 8px;
    }
    [data-popgen-overview="th"][data-sort="desc"]::after {
      content: ' ▼'; font-size: 8px;
    }
    [data-popgen-overview="th"][data-popgen-overview="mark"] { width: 16px; }
    [data-popgen-overview="row"] {
      cursor: pointer;
      border-bottom: 1px solid var(--rule, #c8cdd2);
    }
    [data-popgen-overview="row"]:hover {
      background: rgba(245, 165, 36, 0.08);
    }
    [data-popgen-overview="row"][data-focal="true"] {
      background: rgba(245, 165, 36, 0.18);
      font-weight: 600;
    }
    [data-popgen-overview="row"][data-multi="true"] {
      background: rgba(91, 126, 79, 0.10);
    }
    [data-popgen-overview="row"][data-excluded="true"] {
      opacity: 0.4; text-decoration: line-through;
    }
    [data-popgen-overview="row"][data-confirmed="true"] {
      box-shadow: inset 3px 0 0 var(--good, #3cc08a);
    }
    [data-popgen-overview="td"], [data-popgen-overview="cand-id"] {
      padding: 3px 6px;
      border-right: 1px solid var(--rule, #c8cdd2);
      white-space: nowrap;
      font-variant-numeric: tabular-nums;
    }
    [data-popgen-overview="cand-id"] {
      font-weight: 500;
      max-width: 180px; overflow: hidden; text-overflow: ellipsis;
    }
    [data-popgen-overview="td"][data-kind="live"] { color: var(--accent, #f5a524); }
    [data-popgen-overview="td"][data-kind="flag"] { text-align: center; }
    [data-popgen-overview="mark"] {
      text-align: center; cursor: pointer; user-select: none;
      width: 16px; padding: 3px 0;
    }
    [data-popgen-overview="status"] {
      padding: 4px 10px;
      background: var(--panel-2, #f0f1f3);
      border-top: 1px solid var(--rule, #c8cdd2);
      font-size: 10px; color: var(--ink-dim, #555e69);
      flex-shrink: 0;
    }

    /* Comparison view */
    [data-popgen-overview="compare"] {
      padding: 8px;
      background: var(--panel-2, #f0f1f3);
      border-radius: 3px;
      margin-top: 8px;
    }
    [data-popgen-overview="compare-title"] {
      font-size: 10px; font-weight: 600;
      letter-spacing: 0.06em;
      color: var(--ink-dim, #555e69);
      margin-bottom: 6px;
    }
    [data-popgen-overview="compare-table"] {
      width: 100%; border-collapse: collapse;
      font-family: var(--mono, monospace); font-size: 10px;
    }
    [data-popgen-overview="compare-table"] th,
    [data-popgen-overview="compare-table"] td {
      padding: 3px 6px; text-align: left;
      border-bottom: 1px solid var(--rule, #c8cdd2);
    }
    [data-popgen-overview="delta"][data-sign="pos"] { color: var(--good, #3cc08a); }
    [data-popgen-overview="delta"][data-sign="neg"] { color: var(--bad,  #e0555c); }
  `;

  // ===========================================================================
  // Public exports
  // ===========================================================================

  return {
    makeOverview,
    computeRowStats,
    computeAllRows,
    setRowFlag, getRowFlag,
    clearRowStatsCache,
    exportTSV,
    makeComparisonView,
    COLUMNS,
    _internals: {
      _summarizeShelfStats, _sortRows, _filterRows,
      _readCache, _writeCache, _readFlags, _writeFlags,
      _setLiveStats, _liveStat, _flagFor, _bpToMb, _countByRegime,
    },
  };
}));
