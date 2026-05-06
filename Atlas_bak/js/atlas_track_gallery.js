// =============================================================================
// atlas_track_gallery.js — Turn 6.5 of chat A
// =============================================================================
//
// Discovers all addressable tracks from a popstats server response and
// exposes them as a gallery sidebar on page 6. Provides:
//
//   1. discoverTracks(payload) — scan response columns, return a flat list of
//      track defs (one per metric-group combination), grouped by metric family
//   2. makeGallery() — DOM element for the gallery sidebar (collapsible tray)
//   3. View presets — named curated subsets of track ids:
//        - 'q04_stack'      (default): the canonical 10-track Q04 layout
//        - 'theta_deep'             : every per-group θπ + cohort θπ + Tajima D
//        - 'fst_all_pairs'          : every Fst_<a>_<b> + dXY_<a>_<b>
//        - 'het_signature'          : HoverE + Hobs + Hexp per group + Δ12
//        - 'all_metrics'            : everything discovered
//        - user-saved presets via savePreset(name, trackIds)
//
// API (window.popgenGallery):
//   .discoverTracks(payload, hobsPayload?)
//   .makeGallery({ onTrackToggle, onPresetApply })
//   .getActiveTrackIds()
//   .setActiveTrackIds(ids)
//   .listPresets()
//   .applyPreset(name)
//   .savePreset(name, ids)
//   .deletePreset(name)
//
// Persistence:
//   - Active track id list, scoped per cohort, in localStorage
//   - User-saved presets, scoped per cohort
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenGallery = factory();
    if (typeof document !== 'undefined') {
      const tryAutoWire = () => {
        if (root.popgenLive && root.popgen) {
          // After every successful Compute, refresh discovered tracks
          if (root.popgenLive.popstatsGroupwise && !root.popgenLive.__galleryPatched) {
            const orig = root.popgenLive.popstatsGroupwise;
            root.popgenLive.popstatsGroupwise = async function () {
              const env = await orig.apply(this, arguments);
              if (env && env.ok && env.payload) {
                root.popgenGallery.refreshFromState();
              }
              return env;
            };
            root.popgenLive.__galleryPatched = true;
          }
        } else { setTimeout(tryAutoWire, 100); }
      };
      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', tryAutoWire);
      } else { tryAutoWire(); }
    }
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const LS_PREFIX        = 'inversion_atlas.gallery.';
  const LS_ACTIVE_TRACKS = LS_PREFIX + 'activeTracks';
  const LS_PRESETS       = LS_PREFIX + 'userPresets';

  // The Q04 default — exactly the 10 tracks the static PDF carries.
  // These are CANONICAL track ids; the discovery layer builds matching defs
  // when the underlying data is present.
  const Q04_STACK = [
    'ideogram', 'sim_collapse', 'z',
    'snp_density', 'beagle_unc', 'coverage', 'low_cov_count',
    'theta_invgt', 'fst_hom1_hom2', 'hobs_hexp',
    'delta12', 'delta12_multi',
  ];

  // Built-in preset definitions (resolve dynamically against discovered tracks)
  const BUILTIN_PRESETS = {
    q04_stack: {
      label: 'Q04 stack (default)',
      // Always-on items + Q04 placeholders. Used as the initial active set.
      track_ids: Q04_STACK,
    },
    theta_deep: {
      label: 'θπ deep dive',
      // Filter pattern: 'theta_pi_*', 'theta_pi_cohort', 'tajima_d', 'theta_w'
      pattern: ['theta_pi_*', 'theta_pi_cohort', 'tajima_d'],
    },
    fst_all_pairs: {
      label: 'Fst all-pairs',
      pattern: ['fst_*', 'dxy_*', 'da_*'],
    },
    het_signature: {
      label: 'HET signature',
      pattern: ['hobs_hexp', 'hobs_mean_per_group', 'hexp_mean_per_group', 'delta12*'],
    },
    all_metrics: {
      label: 'All metrics',
      pattern: ['*'],
    },
  };

  // Metric family classification for grouping in the sidebar
  const METRIC_FAMILIES = [
    { key: 'composite', label: 'Composite (canonical)',
      ids: ['ideogram', 'sim_collapse', 'z',
            'snp_density', 'beagle_unc', 'coverage', 'low_cov_count',
            'theta_invgt', 'fst_hom1_hom2', 'hobs_hexp',
            'delta12', 'delta12_multi'] },
    { key: 'theta_pi', label: 'θπ (per group + cohort)',
      pattern: 'theta_pi_*' },
    { key: 'fst',      label: 'Fst (pairwise)',
      pattern: 'fst_*' },
    { key: 'dxy',      label: 'dXY (pairwise)',
      pattern: 'dxy_*' },
    { key: 'da',       label: 'dA (pairwise)',
      pattern: 'da_*' },
    { key: 'mi',       label: 'MI (pairwise)',
      pattern: 'mi_*' },
    { key: 'minorm',   label: 'MInorm (pairwise)',
      pattern: 'minorm_*' },
    { key: 'het',      label: 'HoverE / Hobs / Hexp',
      pattern: ['hobs_*', 'hexp_*'] },
    { key: 'cohort',   label: 'Cohort stats',
      ids: ['theta_pi_cohort', 'tajima_d', 'theta_w_cohort'] },
    { key: 'ancestry', label: 'Ancestry',
      pattern: ['ancestry_q*', 'delta*'] },
    { key: 'tracksdict', label: 'Static tracks dict (precomp)',
      pattern: 'tracksdict_*' },
  ];

  // ===========================================================================
  // Engine accessors
  // ===========================================================================

  function _atlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }
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
  function _cohortKey() {
    const eng = _engine();
    if (eng && typeof eng.getCohortKey === 'function') {
      try { return eng.getCohortKey(); } catch (_) {}
    }
    return 'default';
  }

  // ===========================================================================
  // Persistence
  // ===========================================================================

  function _lsKey(suffix) { return suffix + '.' + _cohortKey(); }
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

  // ===========================================================================
  // Track discovery — scan a popstats response, generate track defs
  // ===========================================================================

  // Each discovered track has the same shape the renderer dispatch expects:
  //   { id, label, height, renderer, color?, yLabel?, edgeTop?, edgeBot?,
  //     family, fromColumn, getData(s) }
  //
  // family classifies for the sidebar grouping. fromColumn is the source
  // column id from the response (e.g., 'theta_pi_HOM1'). getData reads from
  // state.popstatsLive at render time.

  function discoverTracks(payload, hobsPayload) {
    const out = [];

    // 1. Composite/canonical tracks — always available; data presence checked
    //    by the existing static path. We expose the full Q04 stack as
    //    "addressable" so the user can re-add them if they removed.
    const canonicalDefs = _canonicalTrackDefs();
    for (const def of canonicalDefs) out.push(def);

    // 2. Discover from popstats response columns
    if (payload && Array.isArray(payload.windows) && payload.windows.length > 0) {
      const cols = payload.columns || Object.keys(payload.windows[0]);
      const seen = new Set(out.map(t => t.id));
      for (const col of cols) {
        const def = _trackDefForColumn(col, payload);
        if (!def) continue;
        if (seen.has(def.id)) continue;
        seen.add(def.id);
        out.push(def);
      }
    }

    // 3. Discover from hobs response (per-group HoverE/Hobs/Hexp at scales)
    if (hobsPayload && hobsPayload.groups) {
      const renderers = _renderers();
      if (renderers) {
        const adapted = renderers.adaptHobsResponse(hobsPayload);
        for (const trackId of Object.keys(adapted)) {
          if (out.find(t => t.id === trackId)) continue;
          out.push({
            id:       trackId,
            label:    _prettyTrackName(trackId),
            height:   90,
            renderer: 'multiline',
            yLabel:   trackId,
            family:   'het',
            fromHobs: true,
          });
        }
      }
    }

    return out;
  }

  function _canonicalTrackDefs() {
    return [
      { id: 'ideogram',     label: 'ideogram',  height: 50, renderer: 'ideogram',
        alwaysOn: true, family: 'composite' },
      { id: 'sim_collapse', label: 'sim_mat',   height: 50, renderer: 'sim_collapse',
        family: 'composite' },
      { id: 'z',            label: 'Z',         height: 110, renderer: 'line',
        color: '#1f4e79', yLabel: 'Robust Z', edgeTop: 'outlier', edgeBot: 'typical',
        family: 'composite' },
      { id: 'snp_density',  label: 'SNP density', height: 90, renderer: 'line',
        color: '#2c7a39', yLabel: 'SNPs/10kb',
        edgeTop: 'dense', edgeBot: 'sparse', family: 'composite' },
      { id: 'beagle_unc',   label: 'BEAGLE uncert', height: 90, renderer: 'line',
        color: '#7b3294', yLabel: 'uncert frac',
        edgeTop: 'low-conf', edgeBot: 'confident', family: 'composite' },
      { id: 'coverage',     label: 'coverage', height: 90, renderer: 'line',
        color: '#c0504d', yLabel: 'mean cov', family: 'composite' },
      { id: 'low_cov_count', label: '# low cov', height: 60, renderer: 'bars',
        color: '#c0504d', yLabel: '# low cov',
        edgeTop: 'many low', edgeBot: 'all OK', family: 'composite' },
      { id: 'theta_invgt',  label: 'θπ by invgt', height: 90, renderer: 'multiline',
        yLabel: 'θπ invgt', edgeTop: 'diverse', edgeBot: 'low π',
        family: 'composite' },
      { id: 'fst_hom1_hom2', label: 'Fst Hom1-Hom2', height: 90, renderer: 'line',
        color: '#7b3294', yLabel: 'Fst',
        edgeTop: 'differentiated', edgeBot: 'panmictic', family: 'composite' },
      { id: 'hobs_hexp',    label: 'Hobs/Hexp',  height: 90, renderer: 'multiline',
        yLabel: 'Hobs/Hexp',
        edgeTop: 'het excess (~2)', edgeBot: 'hom deficit (~0)',
        family: 'composite' },
      { id: 'delta12',      label: 'ancestry Δ12', height: 90, renderer: 'line',
        color: '#2c7a39', yLabel: 'Δ12',
        edgeTop: 'clear', edgeBot: 'ambiguous', family: 'composite' },
      { id: 'delta12_multi', label: 'Δ12 multi-scale', height: 90, renderer: 'multiline',
        yLabel: 'Δ12 (1×/5×/10×)',
        edgeTop: 'scale-stable', edgeBot: 'scale-dep', family: 'composite' },
    ];
  }

  // Column → track def. Some columns are atomic (a single line per window),
  // some belong to canonical composites (e.g., theta_pi_HOM1 also feeds the
  // theta_invgt aggregate, but stays addressable as its own track).
  function _trackDefForColumn(col, payload) {
    if (!col || typeof col !== 'string') return null;
    // Skip non-data columns
    if (['window_id', 'chrom', 'start', 'end', 'center_mb',
         'n_sites', 'n_sites_used', 'S'].indexOf(col) >= 0) return null;

    if (col === 'theta_pi') {
      return {
        id: 'theta_pi_cohort', label: 'θπ (cohort)', height: 80,
        renderer: 'line', color: '#7a8398',
        yLabel: 'θπ', family: 'cohort', fromColumn: col,
      };
    }
    if (col === 'theta_w') {
      return {
        id: 'theta_w_cohort', label: 'Watterson θ', height: 80,
        renderer: 'line', color: '#7a8398',
        yLabel: 'θ_w', family: 'cohort', fromColumn: col,
      };
    }
    if (col === 'tajD') {
      return {
        id: 'tajima_d', label: "Tajima's D", height: 80,
        renderer: 'line', color: '#9bbf7c',
        yLabel: "Tajima's D", family: 'cohort', fromColumn: col,
      };
    }
    if (col === 'het') {
      return {
        id: 'het_cohort', label: 'het (cohort)', height: 80,
        renderer: 'line', color: '#5a8fb3',
        yLabel: 'het', family: 'het', fromColumn: col,
      };
    }
    // Per-group θπ: theta_pi_<groupName>
    let m = /^theta_pi_(.+)$/.exec(col);
    if (m) {
      const g = m[1];
      return {
        id: 'theta_pi_' + g, label: 'θπ ' + g, height: 80,
        renderer: 'line', color: _colorForGroup(g),
        yLabel: 'θπ ' + g, family: 'theta_pi', fromColumn: col,
      };
    }
    // Pairwise Fst: Fst_<a>_<b>
    m = /^Fst_(.+)_(.+)$/.exec(col);
    if (m) {
      const a = m[1], b = m[2];
      return {
        id: 'fst_' + a.toLowerCase() + '_' + b.toLowerCase(),
        label: 'Fst ' + a + '–' + b, height: 80,
        renderer: 'line', color: '#7b3294',
        yLabel: 'Fst', refLine: 0,
        family: 'fst', fromColumn: col,
      };
    }
    m = /^dXY_(.+)_(.+)$/.exec(col);
    if (m) {
      const a = m[1], b = m[2];
      return {
        id: 'dxy_' + a.toLowerCase() + '_' + b.toLowerCase(),
        label: 'dXY ' + a + '–' + b, height: 80,
        renderer: 'line', color: '#3b6c8c',
        yLabel: 'dXY',
        family: 'dxy', fromColumn: col,
      };
    }
    m = /^dA_(.+)_(.+)$/.exec(col);
    if (m) {
      const a = m[1], b = m[2];
      return {
        id: 'da_' + a.toLowerCase() + '_' + b.toLowerCase(),
        label: 'dA ' + a + '–' + b, height: 80,
        renderer: 'line', color: '#3cc08a',
        yLabel: 'dA',
        family: 'da', fromColumn: col,
      };
    }
    m = /^MI_(.+)_(.+)$/.exec(col);
    if (m) {
      const a = m[1], b = m[2];
      return {
        id: 'mi_' + a.toLowerCase() + '_' + b.toLowerCase(),
        label: 'MI ' + a + '–' + b, height: 80,
        renderer: 'line', color: '#b07cf7',
        yLabel: 'MI',
        family: 'mi', fromColumn: col,
      };
    }
    m = /^MInorm_(.+)_(.+)$/.exec(col);
    if (m) {
      const a = m[1], b = m[2];
      return {
        id: 'minorm_' + a.toLowerCase() + '_' + b.toLowerCase(),
        label: 'MInorm ' + a + '–' + b, height: 80,
        renderer: 'line', color: '#d28cb8',
        yLabel: 'MInorm',
        family: 'minorm', fromColumn: col,
      };
    }
    return null;
  }

  // Deterministic color per group name
  const _GROUP_PALETTE = [
    '#1f4e79', '#f5a524', '#3cc08a', '#e0555c', '#b07cf7',
    '#7ad3db', '#ff8c6e', '#9bbf7c', '#d28cb8', '#5a8fb3',
  ];
  function _colorForGroup(name) {
    let h = 0;
    for (let i = 0; i < name.length; i++) h = ((h << 5) - h + name.charCodeAt(i)) | 0;
    return _GROUP_PALETTE[Math.abs(h) % _GROUP_PALETTE.length];
  }

  function _prettyTrackName(id) {
    return id.replace(/_/g, ' ');
  }

  // ===========================================================================
  // Active tracks state
  // ===========================================================================

  function getActiveTrackIds() {
    const stored = _lsGet('activeTracks', null);
    if (Array.isArray(stored) && stored.length > 0) return stored.slice();
    return Q04_STACK.slice();   // initial default
  }

  function setActiveTrackIds(ids) {
    if (!Array.isArray(ids)) return;
    const clean = ids.filter(x => typeof x === 'string');
    _lsSet('activeTracks', clean);
    _emit('popgen:gallery-changed', { activeTrackIds: clean });
  }

  function isTrackActive(id) {
    const set = new Set(getActiveTrackIds());
    return set.has(id);
  }

  function toggleTrack(id) {
    const cur = getActiveTrackIds();
    const i = cur.indexOf(id);
    if (i >= 0) cur.splice(i, 1);
    else        cur.push(id);
    setActiveTrackIds(cur);
  }

  function _emit(name, detail) {
    if (typeof window !== 'undefined' && window.dispatchEvent) {
      window.dispatchEvent(new CustomEvent(name, { detail }));
    }
  }

  // ===========================================================================
  // Presets
  // ===========================================================================

  function listPresets() {
    const userPresets = _lsGet('userPresets', {}) || {};
    const out = [];
    for (const k of Object.keys(BUILTIN_PRESETS)) {
      out.push({ name: k, label: BUILTIN_PRESETS[k].label, builtin: true });
    }
    for (const k of Object.keys(userPresets)) {
      out.push({ name: k, label: '★ ' + k, builtin: false });
    }
    return out;
  }

  // Resolve a preset to a concrete track-id list given the currently
  // discovered universe of tracks.
  function _resolvePreset(name, discoveredTracks) {
    const userPresets = _lsGet('userPresets', {}) || {};
    if (BUILTIN_PRESETS[name]) {
      const p = BUILTIN_PRESETS[name];
      if (Array.isArray(p.track_ids)) return p.track_ids.slice();
      if (Array.isArray(p.pattern) || typeof p.pattern === 'string') {
        const patterns = Array.isArray(p.pattern) ? p.pattern : [p.pattern];
        const ids = [];
        for (const t of discoveredTracks) {
          for (const pat of patterns) {
            if (_matchPattern(t.id, pat)) { ids.push(t.id); break; }
          }
        }
        // Always keep ideogram + sim_collapse pinned at the top
        if (!ids.includes('ideogram'))     ids.unshift('ideogram');
        if (!ids.includes('sim_collapse')) ids.splice(1, 0, 'sim_collapse');
        return ids;
      }
    }
    if (userPresets[name] && Array.isArray(userPresets[name])) {
      return userPresets[name].slice();
    }
    return null;
  }

  function _matchPattern(id, pattern) {
    if (!pattern) return false;
    if (pattern === '*') return true;
    if (pattern === id)  return true;
    // Glob-style: 'theta_pi_*' matches 'theta_pi_HOM1'
    if (pattern.endsWith('*')) {
      return id.startsWith(pattern.slice(0, -1));
    }
    return false;
  }

  function applyPreset(name) {
    const a = _atlasState();
    const payload = (a && a.popstatsLive) ? a.popstatsLive.lastResponse : null;
    const hobsPayload = (a && a.popstatsLive) ? a.popstatsLive.lastHobsResponse : null;
    const tracks = discoverTracks(payload, hobsPayload);
    const ids = _resolvePreset(name, tracks);
    if (!ids) return false;
    setActiveTrackIds(ids);
    return true;
  }

  function savePreset(name, ids) {
    if (!/^[A-Za-z0-9_]+$/.test(name)) {
      throw new Error('preset name: [A-Za-z0-9_]+');
    }
    if (BUILTIN_PRESETS[name]) {
      throw new Error('cannot overwrite built-in preset: ' + name);
    }
    if (!Array.isArray(ids)) throw new Error('ids must be an array');
    const userPresets = _lsGet('userPresets', {}) || {};
    userPresets[name] = ids.slice();
    _lsSet('userPresets', userPresets);
  }

  function deletePreset(name) {
    const userPresets = _lsGet('userPresets', {}) || {};
    if (BUILTIN_PRESETS[name]) {
      throw new Error('cannot delete built-in preset: ' + name);
    }
    delete userPresets[name];
    _lsSet('userPresets', userPresets);
  }

  // ===========================================================================
  // Group discovered tracks by family for sidebar rendering
  // ===========================================================================

  function _groupTracksByFamily(tracks) {
    const out = [];
    for (const fam of METRIC_FAMILIES) {
      const matched = [];
      for (const t of tracks) {
        if (fam.ids && fam.ids.includes(t.id)) { matched.push(t); continue; }
        if (fam.pattern) {
          const pats = Array.isArray(fam.pattern) ? fam.pattern : [fam.pattern];
          for (const p of pats) {
            if (_matchPattern(t.id, p)) { matched.push(t); break; }
          }
        } else if (t.family === fam.key) {
          matched.push(t);
        }
      }
      if (matched.length > 0) out.push({ family: fam, tracks: matched });
    }
    // Anything not matched goes into "Other"
    const matchedIds = new Set();
    for (const g of out) for (const t of g.tracks) matchedIds.add(t.id);
    const orphan = tracks.filter(t => !matchedIds.has(t.id));
    if (orphan.length > 0) {
      out.push({ family: { key: 'other', label: 'Other' }, tracks: orphan });
    }
    return out;
  }

  // ===========================================================================
  // Gallery DOM
  // ===========================================================================

  let _galleryRoot = null;
  let _searchQuery = '';
  let _expandedFamilies = new Set(['composite', 'theta_pi', 'fst', 'het']);

  function makeGallery(opts) {
    opts = opts || {};
    _injectCSS();
    const root = document.createElement('div');
    root.setAttribute('data-popgen-gallery', 'root');

    // Header: title + preset dropdown + close
    const header = document.createElement('div');
    header.setAttribute('data-popgen-gallery', 'header');
    const title = document.createElement('span');
    title.textContent = 'tracks';
    title.setAttribute('data-popgen-gallery', 'title');
    header.appendChild(title);

    const presetSel = document.createElement('select');
    presetSel.setAttribute('data-popgen-gallery', 'preset');
    presetSel.title = 'Apply a preset view';
    _populatePresetSelect(presetSel);
    presetSel.addEventListener('change', () => {
      if (!presetSel.value) return;
      if (presetSel.value === '__save__') {
        const name = (prompt('Save current track selection as preset:') || '').trim();
        if (!name) { presetSel.value = ''; return; }
        try {
          savePreset(name, getActiveTrackIds());
          _populatePresetSelect(presetSel);
          presetSel.value = name;
        } catch (e) {
          alert('Save failed: ' + (e.message || e));
          presetSel.value = '';
        }
      } else {
        applyPreset(presetSel.value);
        _refreshContent(root);
        if (opts.onPresetApply) opts.onPresetApply(presetSel.value);
      }
    });
    header.appendChild(presetSel);

    if (opts.onClose) {
      const closeBtn = document.createElement('button');
      closeBtn.setAttribute('data-popgen-gallery', 'close');
      closeBtn.textContent = '×';
      closeBtn.title = 'Hide gallery';
      closeBtn.addEventListener('click', opts.onClose);
      header.appendChild(closeBtn);
    }

    root.appendChild(header);

    // Search box
    const searchBox = document.createElement('input');
    searchBox.setAttribute('data-popgen-gallery', 'search');
    searchBox.placeholder = 'filter tracks…';
    searchBox.value = _searchQuery;
    searchBox.addEventListener('input', () => {
      _searchQuery = searchBox.value;
      _refreshContent(root);
    });
    root.appendChild(searchBox);

    // Content body
    const body = document.createElement('div');
    body.setAttribute('data-popgen-gallery', 'body');
    root.appendChild(body);

    _galleryRoot = root;
    _refreshContent(root, opts.onTrackToggle);
    return root;
  }

  function _populatePresetSelect(sel) {
    while (sel.firstChild) sel.removeChild(sel.firstChild);
    const blank = document.createElement('option');
    blank.value = ''; blank.textContent = '— preset —';
    sel.appendChild(blank);
    for (const p of listPresets()) {
      const o = document.createElement('option');
      o.value = p.name; o.textContent = p.label;
      sel.appendChild(o);
    }
    const saveOpt = document.createElement('option');
    saveOpt.value = '__save__'; saveOpt.textContent = '+ save current as preset…';
    sel.appendChild(saveOpt);
  }

  function _refreshContent(root, onToggle) {
    if (!root) return;
    const body = root.querySelector('[data-popgen-gallery="body"]');
    if (!body) return;
    while (body.firstChild) body.removeChild(body.firstChild);
    const a = _atlasState();
    const payload = (a && a.popstatsLive) ? a.popstatsLive.lastResponse : null;
    const hobsPayload = (a && a.popstatsLive) ? a.popstatsLive.lastHobsResponse : null;
    const tracks = discoverTracks(payload, hobsPayload);
    const filtered = _searchQuery
      ? tracks.filter(t => _matchSearch(t, _searchQuery))
      : tracks;
    const grouped = _groupTracksByFamily(filtered);
    const activeSet = new Set(getActiveTrackIds());

    if (grouped.length === 0) {
      const empty = document.createElement('div');
      empty.setAttribute('data-popgen-gallery', 'empty');
      empty.textContent = _searchQuery
        ? '— no matches —'
        : '— no tracks discovered yet. Compute to populate. —';
      body.appendChild(empty);
      return;
    }

    for (const grp of grouped) {
      const wrap = document.createElement('div');
      wrap.setAttribute('data-popgen-gallery', 'cat');
      const open = _expandedFamilies.has(grp.family.key) || !!_searchQuery;
      wrap.setAttribute('data-open', open ? 'true' : 'false');
      const head = document.createElement('div');
      head.setAttribute('data-popgen-gallery', 'cat-head');
      const tog = document.createElement('span');
      tog.setAttribute('data-popgen-gallery', 'cat-toggle');
      head.appendChild(tog);
      const labelSpan = document.createElement('span');
      labelSpan.textContent = grp.family.label;
      head.appendChild(labelSpan);
      const cnt = document.createElement('span');
      cnt.setAttribute('data-popgen-gallery', 'cat-count');
      cnt.textContent = '(' + grp.tracks.length + ')';
      head.appendChild(cnt);
      head.addEventListener('click', () => {
        const isOpen = wrap.getAttribute('data-open') === 'true';
        wrap.setAttribute('data-open', isOpen ? 'false' : 'true');
        if (isOpen) _expandedFamilies.delete(grp.family.key);
        else        _expandedFamilies.add(grp.family.key);
      });
      wrap.appendChild(head);

      const chips = document.createElement('div');
      chips.setAttribute('data-popgen-gallery', 'chips');
      for (const t of grp.tracks) {
        const chip = document.createElement('span');
        chip.setAttribute('data-popgen-gallery', 'chip');
        chip.setAttribute('data-on', activeSet.has(t.id) ? 'true' : 'false');
        chip.title = t.id;
        chip.textContent = t.label || t.id;
        chip.addEventListener('click', () => {
          toggleTrack(t.id);
          chip.setAttribute('data-on',
            getActiveTrackIds().indexOf(t.id) >= 0 ? 'true' : 'false');
          if (onToggle) onToggle(t.id);
        });
        chips.appendChild(chip);
      }
      wrap.appendChild(chips);
      body.appendChild(wrap);
    }
  }

  function _matchSearch(track, q) {
    q = q.toLowerCase();
    if ((track.id || '').toLowerCase().includes(q)) return true;
    if ((track.label || '').toLowerCase().includes(q)) return true;
    if ((track.family || '').toLowerCase().includes(q)) return true;
    return false;
  }

  function refreshFromState() {
    if (!_galleryRoot) return;
    _refreshContent(_galleryRoot);
  }

  // ===========================================================================
  // CSS injection
  // ===========================================================================

  const GALLERY_CSS = `
    [data-popgen-gallery="root"] {
      display: flex; flex-direction: column;
      width: 240px; height: 100%;
      background: var(--panel, #fbfcfd);
      border-left: 1px solid var(--rule, #c8cdd2);
      font-family: var(--mono, ui-monospace, monospace);
      font-size: 11px; color: var(--ink, #0e1116);
      overflow: hidden;
    }
    [data-popgen-gallery="header"] {
      display: flex; align-items: center; gap: 6px;
      padding: 6px 8px;
      background: var(--panel-2, #f0f1f3);
      border-bottom: 1px solid var(--rule, #c8cdd2);
      flex-shrink: 0;
    }
    [data-popgen-gallery="title"] {
      font-weight: 600; letter-spacing: 0.04em; text-transform: uppercase;
      font-size: 10px; color: var(--ink-dim, #555e69); flex: 1;
    }
    [data-popgen-gallery="preset"] {
      font-family: inherit; font-size: 10px;
      padding: 2px 4px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-radius: 3px;
      max-width: 130px;
    }
    [data-popgen-gallery="close"] {
      width: 20px; height: 20px; padding: 0; border: 0;
      background: transparent; cursor: pointer;
      font-size: 14px; line-height: 1; color: var(--ink-dim, #555e69);
      border-radius: 4px;
    }
    [data-popgen-gallery="close"]:hover { background: var(--bad, #e0555c); color: white; }
    [data-popgen-gallery="search"] {
      padding: 4px 8px; margin: 6px 8px;
      border: 1px solid var(--rule, #c8cdd2);
      background: var(--panel, #fbfcfd);
      color: var(--ink, #0e1116);
      border-radius: 3px;
      font-family: inherit; font-size: 10px;
    }
    [data-popgen-gallery="body"] { flex: 1; overflow: auto; padding: 6px 8px; }
    [data-popgen-gallery="empty"] {
      color: var(--ink-dim, #555e69); font-style: italic;
      padding: 12px 4px; font-size: 10px;
    }
    [data-popgen-gallery="cat"] {
      margin-bottom: 4px;
      border: 1px solid var(--rule, #c8cdd2); border-radius: 4px;
      background: var(--panel-2, #f0f1f3);
    }
    [data-popgen-gallery="cat-head"] {
      display: flex; align-items: center; gap: 4px;
      padding: 4px 6px; cursor: pointer; user-select: none;
      font-size: 10px;
    }
    [data-popgen-gallery="cat-head"]:hover { background: var(--panel-3, #e6e8ec); }
    [data-popgen-gallery="cat-toggle"] {
      width: 10px; display: inline-block; color: var(--ink-dim, #555e69);
    }
    [data-popgen-gallery="cat"][data-open="true"] [data-popgen-gallery="cat-toggle"]::before { content: '▾'; }
    [data-popgen-gallery="cat"][data-open="false"] [data-popgen-gallery="cat-toggle"]::before { content: '▸'; }
    [data-popgen-gallery="cat-count"] {
      margin-left: auto; color: var(--ink-dim, #555e69);
      font-size: 9px;
    }
    [data-popgen-gallery="chips"] {
      display: none; flex-wrap: wrap; gap: 3px;
      padding: 4px 6px;
      background: var(--panel, #fbfcfd);
      border-top: 1px solid var(--rule, #c8cdd2);
    }
    [data-popgen-gallery="cat"][data-open="true"] [data-popgen-gallery="chips"] {
      display: flex;
    }
    [data-popgen-gallery="chip"] {
      padding: 2px 6px; border-radius: 8px;
      background: var(--panel-2, #f0f1f3);
      border: 1px solid var(--rule, #c8cdd2);
      font-size: 9px; cursor: pointer;
      user-select: none;
      color: var(--ink, #0e1116);
    }
    [data-popgen-gallery="chip"]:hover {
      border-color: var(--accent, #f5a524);
    }
    [data-popgen-gallery="chip"][data-on="true"] {
      background: var(--accent, #f5a524); color: #0e1116;
      border-color: var(--accent, #f5a524); font-weight: 600;
    }
    [data-popgen-gallery="chip"][data-on="true"]::before { content: '● '; }
    [data-popgen-gallery="chip"][data-on="false"]::before { content: '○ '; }
  `;

  function _injectCSS() {
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-gallery-css')) return;
    const style = document.createElement('style');
    style.id = 'popgen-gallery-css';
    style.textContent = GALLERY_CSS;
    document.head.appendChild(style);
  }

  // ===========================================================================
  // Public exports
  // ===========================================================================

  return {
    discoverTracks,
    makeGallery,
    refreshFromState,
    getActiveTrackIds, setActiveTrackIds, isTrackActive, toggleTrack,
    listPresets, applyPreset, savePreset, deletePreset,
    Q04_STACK,
    BUILTIN_PRESETS,
    METRIC_FAMILIES,
    _internals: {
      _trackDefForColumn, _resolvePreset, _matchPattern, _groupTracksByFamily,
      _canonicalTrackDefs, _matchSearch,
    },
  };
}));
