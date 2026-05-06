// =============================================================================
// dosage_viewer/static/app.js
// =============================================================================
// Standalone dosage heatmap viewer. Talks to 02_run_server.py over /api/*.
//
// Architecture
// ------------
// State shape (all client-side, no localStorage by default):
//
//   state = {
//     manifest:        result of /api/manifest (chroms, samples, candidates_loaded)
//     activeRegion:    last successfully-rendered envelope (used for redraw + tooltips)
//     recent:          LRU<key, { region, ts, result_n_returned }>  (cap 16)
//     pending:         current AbortController (so a Render-during-fetch cancels)
//     theme:           'light' | 'dark'
//     candidates:      [{ id, chrom, start, end, ... }]  (built from manifest)
//   }
//
// Render flow (Render button click):
//   1. Read top-bar inputs -> RegionRequest
//   2. Push into LRU; abort any in-flight fetch
//   3. fetch /api/region with the request
//   4. On 200: stash the envelope as state.activeRegion, draw the canvas
//      (sized to n_sites_returned x n_samples), update meta strip
//   5. On non-200: surface error in the meta strip (no redraw)
//
// Spec invariants honoured:
//   - No auto-rerender on slider drag (Render button only)
//   - Recent-region LRU capped at 16
//   - FIG_C08-style colour scheme
//   - Per-site modes: matrix is int with -1 = NA
//   - Aggregate mode: matrix is float-or-null; null cells render as NA pattern

'use strict';

// ---------------------------------------------------------------------------
// Config & constants
// ---------------------------------------------------------------------------

const RECENT_LRU_CAP = 16;
const CANVAS_WIDTH_PX = 1024;          // fixed per spec; sample rows resize-fit
const CANVAS_MIN_HEIGHT = 60;          // never shrink below this
const CANVAS_MAX_HEIGHT = 800;         // never grow above
const RAMP = ['#2166AC', '#F7F7F7', '#B2182B'];   // ColorBrewer RdBu

const API_BASE = '';   // same-origin
const ENDPOINT_MANIFEST   = '/api/manifest';
const ENDPOINT_REGION     = '/api/region';
const ENDPOINT_CANDIDATE  = '/api/candidate/';
const ENDPOINT_BREAKPOINT = '/api/breakpoint/';
const ENDPOINT_TILE       = '/api/tile';

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
const state = {
  manifest: null,
  activeRegion: null,
  recent: new Map(),     // insertion-order = LRU order
  pending: null,
  theme: 'light',
  candidates: [],
};

// ---------------------------------------------------------------------------
// DOM helpers
// ---------------------------------------------------------------------------
function $(id) { return document.getElementById(id); }
function setMeta(id, text, cls) {
  const el = $(id);
  if (!el) return;
  el.textContent = text;
  el.className = 'meta-pill' + (cls ? (' ' + cls) : '');
}
function setBadge(id, text) { const el = $(id); if (el) el.textContent = text; }

// ---------------------------------------------------------------------------
// Colour ramp — interpolate dosage 0..2 to the FIG_C08 RdBu scheme
// ---------------------------------------------------------------------------
function _hexToRgb(hex) {
  const s = hex.replace('#', '');
  return [
    parseInt(s.slice(0, 2), 16),
    parseInt(s.slice(2, 4), 16),
    parseInt(s.slice(4, 6), 16),
  ];
}
const RGB0 = _hexToRgb(RAMP[0]);
const RGB1 = _hexToRgb(RAMP[1]);
const RGB2 = _hexToRgb(RAMP[2]);

function dosageRGB(d) {
  // d in [0, 2] (per-site modes use ints 0/1/2 with -1=NA; aggregate uses
  // float means).
  if (d == null || !Number.isFinite(d) || d < 0) return null;   // NA
  const x = Math.max(0, Math.min(2, d));
  const a = x < 1 ? RGB0 : RGB1;
  const b = x < 1 ? RGB1 : RGB2;
  const t = x < 1 ? x : x - 1;
  return [
    Math.round(a[0] + (b[0] - a[0]) * t),
    Math.round(a[1] + (b[1] - a[1]) * t),
    Math.round(a[2] + (b[2] - a[2]) * t),
  ];
}

// NA pattern — diagonal stripes baked into a small image data buffer
let _naPatternRGBA = null;
function naPatternRGBA(theme) {
  if (_naPatternRGBA) return _naPatternRGBA;
  const a = theme === 'dark' ? [74, 77, 84] : [204, 204, 204];
  const b = theme === 'dark' ? [31, 33, 39] : [255, 255, 255];
  _naPatternRGBA = { a: a, b: b };
  return _naPatternRGBA;
}

// ---------------------------------------------------------------------------
// Initial bootstrap — fetch manifest, populate selects
// ---------------------------------------------------------------------------
async function bootstrap() {
  setBadge('storeBadge', 'connecting…');
  try {
    const r = await fetch(API_BASE + ENDPOINT_MANIFEST);
    if (!r.ok) throw new Error('HTTP ' + r.status);
    state.manifest = await r.json();
  } catch (e) {
    setBadge('storeBadge', 'offline (' + e.message + ')');
    setMeta('mp_warning', 'cannot reach /api/manifest', 'error');
    return;
  }
  setBadge('storeBadge',
    'connected · ' + state.manifest.n_samples + ' samples · ' +
    state.manifest.chroms.length + ' chroms');

  // Populate chrom selector
  const chromSel = $('chromSelect');
  chromSel.innerHTML = '';
  for (const c of state.manifest.chroms) {
    const opt = document.createElement('option');
    opt.value = c.name;
    opt.textContent = c.name + ' (' + Math.round(c.length_bp / 1e6) + ' Mb · ' +
                       c.n_sites.toLocaleString() + ' sites)';
    chromSel.appendChild(opt);
  }
  // Pick first chrom + set sensible default range
  const first = state.manifest.chroms[0];
  if (first) {
    chromSel.value = first.name;
    $('startInput').value = 1;
    $('endInput').value = Math.min(first.length_bp, 2_000_000);
  }

  // Candidate selector — populated from /api/candidate has no list endpoint
  // in P4 spec, so we display whatever the manifest's
  // candidate_regions_loaded flag says. (For full population we'd add a
  // /api/candidates list endpoint; deferred to next slice.)
  if (state.manifest.candidate_regions_loaded) {
    setBadge('lastFetch', 'candidate config loaded');
  } else {
    setBadge('lastFetch', 'no candidate_regions.tsv configured');
    // Disable the candidate row controls
    $('loadCandBtn').disabled = true;
    $('loadLeftBtn').disabled = true;
    $('loadRightBtn').disabled = true;
  }
}

// ---------------------------------------------------------------------------
// LRU recent-region tracking
// ---------------------------------------------------------------------------
function pushRecent(region, env) {
  const key = `${region.chrom}:${region.start}-${region.end} [${region.mode}]`;
  if (state.recent.has(key)) state.recent.delete(key);
  state.recent.set(key, {
    region,
    ts: Date.now(),
    n_total: env.n_sites_total,
    n_returned: env.n_sites_returned,
  });
  while (state.recent.size > RECENT_LRU_CAP) {
    const k0 = state.recent.keys().next().value;
    state.recent.delete(k0);
  }
  renderRecent();
}

function renderRecent() {
  const ul = $('recentList');
  ul.innerHTML = '';
  // Newest-first
  const entries = [...state.recent.entries()].reverse();
  for (const [key, item] of entries) {
    const li = document.createElement('li');
    const tag = document.createElement('span');
    tag.className = 'recent-tag';
    tag.textContent = item.region.mode;
    li.appendChild(tag);
    const txt = document.createElement('span');
    const ago = Math.round((Date.now() - item.ts) / 1000);
    txt.textContent = `${item.region.chrom} ${item.region.start.toLocaleString()}–${item.region.end.toLocaleString()} bp · ${item.n_returned}/${item.n_total} sites · ${ago}s ago`;
    li.appendChild(txt);
    li.addEventListener('click', () => {
      $('chromSelect').value = item.region.chrom;
      $('startInput').value = item.region.start;
      $('endInput').value = item.region.end;
      $('modeSelect').value = item.region.mode;
      $('maxSitesInput').value = item.region.max_sites;
      $('seedInput').value = item.region.seed;
      doRender();
    });
    ul.appendChild(li);
  }
}

// ---------------------------------------------------------------------------
// Read controls -> RegionRequest
// ---------------------------------------------------------------------------
function readRegionRequest() {
  return {
    chrom: $('chromSelect').value,
    start: parseInt($('startInput').value, 10),
    end:   parseInt($('endInput').value, 10),
    max_sites: parseInt($('maxSitesInput').value, 10),
    mode: $('modeSelect').value,
    seed: parseInt($('seedInput').value, 10) || 1,
  };
}

// ---------------------------------------------------------------------------
// Render flow — fetch /api/region and draw
// ---------------------------------------------------------------------------
async function doRender() {
  const req = readRegionRequest();
  if (!req.chrom) {
    setMeta('mp_warning', 'no chrom selected', 'error');
    return;
  }
  if (!Number.isFinite(req.start) || !Number.isFinite(req.end) ||
      req.end < req.start) {
    setMeta('mp_warning', 'bad start/end', 'error');
    return;
  }

  // Cancel any pending request
  if (state.pending) {
    try { state.pending.abort(); } catch (_) {}
  }
  const ctrl = new AbortController();
  state.pending = ctrl;
  setMeta('mp_warning', 'fetching…', '');

  const params = new URLSearchParams({
    chrom: req.chrom,
    start: String(req.start),
    end:   String(req.end),
    max_sites: String(req.max_sites),
    mode: req.mode,
    seed: String(req.seed),
  });

  const t0 = performance.now();
  let env;
  try {
    const r = await fetch(API_BASE + ENDPOINT_REGION + '?' + params,
                          { signal: ctrl.signal });
    if (!r.ok) {
      const text = await r.text();
      let msg = 'HTTP ' + r.status;
      try {
        const j = JSON.parse(text);
        if (j.detail && j.detail.error) msg = j.detail.error;
        if (j.detail && j.detail.message) msg += ': ' + j.detail.message;
      } catch (_) {}
      setMeta('mp_warning', msg, 'error');
      state.pending = null;
      return;
    }
    env = await r.json();
  } catch (e) {
    if (e.name === 'AbortError') return;
    setMeta('mp_warning', 'fetch failed: ' + (e.message || e), 'error');
    state.pending = null;
    return;
  }
  state.pending = null;
  const ms = Math.round(performance.now() - t0);
  setBadge('lastFetch', `last fetch ${ms} ms`);

  state.activeRegion = env;
  pushRecent(req, env);
  drawHeatmap(env);
  updateMeta(env);
}

// ---------------------------------------------------------------------------
// Canvas rendering
// ---------------------------------------------------------------------------
function drawHeatmap(env) {
  const canvas = $('heatmapCanvas');
  const ctx = canvas.getContext('2d');

  const nSites = env.n_sites_returned || 0;
  const nSamples = (env.samples || []).length;
  if (nSites === 0 || nSamples === 0) {
    canvas.width = CANVAS_WIDTH_PX;
    canvas.height = CANVAS_MIN_HEIGHT;
    ctx.fillStyle = getCSSVar('--bg-panel') || '#fff';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = getCSSVar('--text-dim') || '#888';
    ctx.font = '12px sans-serif';
    ctx.fillText('No sites in region.', 10, 20);
    return;
  }

  const isAggregate = !!(env._meta && env._meta.is_aggregate);
  const tileWPx = Math.max(1, Math.floor(CANVAS_WIDTH_PX / nSites));
  const tileHPx = Math.max(2, Math.min(8, Math.floor(720 / nSamples)));
  const heightPx = Math.max(CANVAS_MIN_HEIGHT,
                             Math.min(CANVAS_MAX_HEIGHT, nSamples * tileHPx));
  // Width may shrink if (tileWPx * nSites) < 1024; that's intentional —
  // we'd rather see crisp 1+ pixel tiles than blurry sub-pixel ones.
  const widthPx = Math.max(CANVAS_WIDTH_PX, tileWPx * nSites);

  canvas.width = widthPx;
  canvas.height = heightPx;
  // Clear canvas to bg-panel
  ctx.fillStyle = getCSSVar('--bg-panel') || '#fff';
  ctx.fillRect(0, 0, widthPx, heightPx);

  // Build an ImageData buffer for performance — single putImageData beats
  // n_sites × n_samples fillRect calls by a lot.
  const id = ctx.createImageData(widthPx, heightPx);
  const buf = id.data;

  const naColors = naPatternRGBA(state.theme);

  for (let i = 0; i < nSites; i++) {
    const row = env.matrix[i];
    if (!row) continue;
    const xStart = i * tileWPx;
    for (let j = 0; j < nSamples; j++) {
      const v = row[j];
      const isNA = (v == null || (typeof v === 'number' && (v < 0 || !Number.isFinite(v))));
      const yStart = j * tileHPx;
      const yEnd = Math.min(heightPx, yStart + tileHPx);
      const xEnd = Math.min(widthPx, xStart + tileWPx);

      let r, g, b;
      if (isNA) {
        // For aggregate, null = empty bin (all-NA); render with hatched NA pattern.
        // For per-site -1 cells, same treatment.
        // We use a simple checker between two greys; pixels with (px+py) even
        // get colour A, odd get colour B.
        for (let py = yStart; py < yEnd; py++) {
          for (let px = xStart; px < xEnd; px++) {
            const odd = (((px - xStart) + (py - yStart)) & 1) === 1;
            const c = odd ? naColors.a : naColors.b;
            const idx = (py * widthPx + px) * 4;
            buf[idx]     = c[0];
            buf[idx + 1] = c[1];
            buf[idx + 2] = c[2];
            buf[idx + 3] = 255;
          }
        }
        continue;
      }
      const rgb = dosageRGB(v);
      if (!rgb) continue;
      r = rgb[0]; g = rgb[1]; b = rgb[2];
      for (let py = yStart; py < yEnd; py++) {
        for (let px = xStart; px < xEnd; px++) {
          const idx = (py * widthPx + px) * 4;
          buf[idx]     = r;
          buf[idx + 1] = g;
          buf[idx + 2] = b;
          buf[idx + 3] = 255;
        }
      }
    }
  }
  ctx.putImageData(id, 0, 0);

  // Stash tile dims for hover lookup
  canvas.dataset.tileW = String(tileWPx);
  canvas.dataset.tileH = String(tileHPx);
  canvas.dataset.nSites = String(nSites);
  canvas.dataset.nSamples = String(nSamples);
}

function getCSSVar(name) {
  if (typeof getComputedStyle !== 'function') return '';
  return getComputedStyle(document.body).getPropertyValue(name).trim();
}

// ---------------------------------------------------------------------------
// Meta-strip update
// ---------------------------------------------------------------------------
function updateMeta(env) {
  setMeta('mp_chrom', env.chrom);
  setMeta('mp_range',
    env.start_bp.toLocaleString() + '–' + env.end_bp.toLocaleString() + ' bp');
  setMeta('mp_sites',
    env.n_sites_returned + ' / ' + env.n_sites_total + ' sites' +
    (env.downsampled ? ' (downsampled)' : ''));
  setMeta('mp_mode', 'mode: ' + env.sampling_mode);
  setMeta('mp_n_samples', env.samples.length + ' samples');
  setMeta('mp_warning', env.warning || '', env.warning ? 'warn' : '');
  // Legend note for aggregate
  const note = $('legendNote');
  if (note) {
    note.textContent = (env._meta && env._meta.is_aggregate)
      ? '(aggregate: cells = mean dosage per bin per sample)'
      : '';
  }
}

// ---------------------------------------------------------------------------
// Hover tooltip
// ---------------------------------------------------------------------------
function attachHoverHandlers() {
  const canvas = $('heatmapCanvas');
  const tip = $('hoverInfo');
  canvas.addEventListener('mousemove', (ev) => {
    if (!state.activeRegion) return;
    const rect = canvas.getBoundingClientRect();
    const px = ev.clientX - rect.left;
    const py = ev.clientY - rect.top;
    const tw = parseInt(canvas.dataset.tileW || '1', 10);
    const th = parseInt(canvas.dataset.tileH || '1', 10);
    const i = Math.floor(px / tw);
    const j = Math.floor(py / th);
    if (i < 0 || j < 0) { tip.hidden = true; return; }
    const env = state.activeRegion;
    if (i >= env.n_sites_returned || j >= env.samples.length) {
      tip.hidden = true; return;
    }
    const pos = env.positions[i];
    const sample = env.samples[j];
    const v = env.matrix[i] && env.matrix[i][j];
    const isAgg = !!(env._meta && env._meta.is_aggregate);
    const valueStr = (v == null || (typeof v === 'number' && v < 0))
      ? 'NA'
      : (isAgg ? Number(v).toFixed(2) : String(v));
    const sid = env.site_ids && env.site_ids[i] ? env.site_ids[i] : '';
    tip.textContent =
      `${env.chrom}:${pos != null ? pos.toLocaleString() : '?'}  · ` +
      `${sample}  ·  dosage = ${valueStr}` +
      (sid ? `  ·  ${sid}` : '');
    tip.style.left = (px + 12) + 'px';
    tip.style.top = (py + 12) + 'px';
    tip.hidden = false;
  });
  canvas.addEventListener('mouseleave', () => { tip.hidden = true; });
}

// ---------------------------------------------------------------------------
// Candidate / breakpoint / tile shortcuts
// ---------------------------------------------------------------------------
async function loadCandidate() {
  const cid = $('candidateSelect').value;
  if (!cid) return;
  const params = new URLSearchParams({
    max_sites: $('maxSitesInput').value,
    mode: $('modeSelect').value,
    seed: $('seedInput').value || '1',
  });
  await fetchAndRender(API_BASE + ENDPOINT_CANDIDATE + encodeURIComponent(cid) + '?' + params);
}

async function loadBreakpoint(side) {
  const cid = $('candidateSelect').value;
  if (!cid) return;
  const window_kbp = parseInt($('bpWindowInput').value, 10) || 500;
  const params = new URLSearchParams({
    window: String(window_kbp * 1000),
    max_sites: $('maxSitesInput').value,
    mode: $('modeSelect').value,
  });
  await fetchAndRender(API_BASE + ENDPOINT_BREAKPOINT +
                        encodeURIComponent(cid) + '/' + side + '?' + params);
}

async function loadTile() {
  const chrom = $('chromSelect').value;
  if (!chrom) return;
  const width = Math.min(2000, parseInt($('maxSitesInput').value, 10) || 500);
  const params = new URLSearchParams({ chrom, width: String(width) });
  await fetchAndRender(API_BASE + ENDPOINT_TILE + '?' + params);
}

async function fetchAndRender(url) {
  if (state.pending) { try { state.pending.abort(); } catch (_) {} }
  const ctrl = new AbortController();
  state.pending = ctrl;
  setMeta('mp_warning', 'fetching…', '');
  try {
    const r = await fetch(url, { signal: ctrl.signal });
    if (!r.ok) {
      const text = await r.text();
      let msg = 'HTTP ' + r.status;
      try {
        const j = JSON.parse(text);
        if (j.detail && j.detail.error) msg = j.detail.error;
      } catch (_) {}
      setMeta('mp_warning', msg, 'error');
      return;
    }
    const env = await r.json();
    state.activeRegion = env;
    // Reflect the new region in the inputs
    $('chromSelect').value = env.chrom;
    $('startInput').value = env.start_bp;
    $('endInput').value = env.end_bp;
    pushRecent({
      chrom: env.chrom, start: env.start_bp, end: env.end_bp,
      mode: env.sampling_mode, max_sites: env.n_sites_returned, seed: 1,
    }, env);
    drawHeatmap(env);
    updateMeta(env);
  } finally {
    state.pending = null;
  }
}

// ---------------------------------------------------------------------------
// Wire events
// ---------------------------------------------------------------------------
function wire() {
  $('renderBtn').addEventListener('click', doRender);
  $('themeToggle').addEventListener('click', () => {
    state.theme = (state.theme === 'light') ? 'dark' : 'light';
    document.body.dataset.theme = state.theme;
    _naPatternRGBA = null;   // recompute NA palette
    if (state.activeRegion) drawHeatmap(state.activeRegion);
  });
  $('loadCandBtn').addEventListener('click', loadCandidate);
  $('loadLeftBtn').addEventListener('click', () => loadBreakpoint('left'));
  $('loadRightBtn').addEventListener('click', () => loadBreakpoint('right'));
  $('tileBtn').addEventListener('click', loadTile);

  // Pressing Enter inside any region-input should not auto-render — per spec.
  // Escape clears any pending fetch.
  document.addEventListener('keydown', (ev) => {
    if (ev.key === 'Escape' && state.pending) {
      try { state.pending.abort(); } catch (_) {}
    }
  });
  attachHoverHandlers();
}

// ---------------------------------------------------------------------------
// Entry
// ---------------------------------------------------------------------------
(async function init() {
  wire();
  await bootstrap();
})();
