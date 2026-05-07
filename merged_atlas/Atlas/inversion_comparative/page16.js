// page16.js — Cross-species breakpoints page (cross-species comparative dashboard)
// Extracted verbatim from legacy/Inversion_atlas.html for batch 5 of the parallel migration.
//
// Source line ranges in legacy:
//   20971-21114 — Section header + cross-species init/IO + constants (CROSS_SPECIES_*, CS_EVENT_DEF)
//   23717-26025 — Cross-species runtime: filter/sort, render*, ideograms, flank charts, synteny, dotplot, focal-vs-bg
//   28367-28419 — _csBuildPermResultHtml (permutation test result HTML)
//
// All function bodies are unmodified. No imports yet — see TODO_MISSING markers below.
// At top level the module is plain JS; the merge chat will wrap into ES module exports as needed.

//
// =============================================================================
// TODO_MISSING markers (frequency × symbol — likely legacy line)
// =============================================================================
// All references below are unresolved within this module; the merge chat
// resolves them by either promoting to shared/ (multi-page consumers) or by
// copying locally (single-page consumers).
//
//  TODO_MISSING(_esc)                  ×61   legacy line 13925   HTML escape helper, used everywhere
//  TODO_MISSING(_getRepeatDensity)     × 3   legacy line 14477   accessor for state.repeatDensity by chrom
//  TODO_MISSING(setCur)                × 1   legacy line 51747   scrubber: set current window index (jump)
//  TODO_MISSING(drawZ)                 × 1   legacy line 31839   scrubber: |Z| panel draw
//  TODO_MISSING(drawSim)               × 1   legacy line 31331   scrubber: similarity-matrix panel draw
//  TODO_MISSING(drawLinesPanel)        × 1   legacy line 34894   scrubber: lines/strip panel draw
//  TODO_MISSING(drawWinSumStrip)       × 1   not defined in legacy — guarded by typeof check, optional hook
//
// External library globals (vendor):
//  TODO_MISSING(window.popgenDotplot)        external — popgen dotplot library (csDotplot wiring)
//  TODO_MISSING(window.popgenFocalVsBg)      external — popgen focal-vs-background library
//
// Top-level state slots this module reads (handled by shared/state.js):
//   state.crossSpecies, state.repeatDensity, state.candidateList, state.dotplotMashmap,
//   state.cur, state.data, state._crossSpeciesUI, state._csSyntenyCache,
//   state._csSyntenyEdgesCache, state._csInversionContextCache,
//   state._csOverlayIndex, state._focalVsBg, state._csDotplotPanel*,
//   state._csHoverActive, state._csHoverRaf, state._csHeaderScrollWired,
//   state._crossSpeciesKeysBound
// =============================================================================
//

// =============================================================================
// BLOCK A — section header + constants + IO/state functions (legacy 20971-21114)
// =============================================================================
// =============================================================================
// v4 turn 74 — Cross-species breakpoints page (page16)
// -----------------------------------------------------------------------------
// Catalogue of chromosome-scale rearrangements between Cgar and Cmac, derived
// from a wfmash 1-to-1 alignment of the two haplotypes (cs_breakpoints_v1
// JSON, output of STEP_CS01_extract_breakpoints.py). Each breakpoint is
// rendered with both species' coordinates, a syntenic-block linking line,
// and the flanking repeat-element density on both species (drawing from
// the boundaries-page TEfull JSONs already loaded into state.repeatDensity).
// Spalax-style enrichment of all_TE at breakpoints is the manuscript hook.
//
// Schema (cs_breakpoints_v1):
//   {tool, schema_version, generated_at, species_query, species_target,
//    input_paf, params, n_breakpoints, n_by_event_type, breakpoints: [...]}
//
// Each breakpoint:
//   {id, event_type, event_type_refined?, gar_chr, gar_pos_start,
//    gar_pos_end, gar_pos_mb, n_member_breakpoints?,
//    prev_block: {mac_chr, mac_start_bp, mac_end_bp, strand, block_size_bp, mapping_quality},
//    next_block: {...},
//    flanking_repeat_density_gar: {<class>: {mean, max, n_windows}, ...},
//    flanking_repeat_density_mac: {prev: {mac_chr, anchor_bp, by_class}, next: {...}},
//    candidate_overlap?: [<candidate_id>, ...],
//    manuscript_note?: <string>}
// =============================================================================

const CROSS_SPECIES_LS_KEY  = 'pca_scrubber_v3.crossSpecies.v1';
const CROSS_SPECIES_TOOL    = 'cross_species_breakpoints_v1';
const CROSS_SPECIES_FLANK_DEFAULT_BP = 100000;
// Event-type label/colour map. Keys match both event_type and event_type_refined.
const CS_EVENT_DEF = {
  inversion:               { label: 'inversion',        cls: 'cs-evt-inv',   icon: '\u25C6' },
  translocation:           { label: 'translocation',    cls: 'cs-evt-trans', icon: '\u21A6' },
  fission_or_fusion:       { label: 'fission/fusion',   cls: 'cs-evt-fis',   icon: '\u2702' },
  translocation_or_fission:{ label: 'transloc.|fission', cls: 'cs-evt-trans', icon: '\u21A6' },
  mixed:                   { label: 'mixed',            cls: 'cs-evt-mixed', icon: '\u25E7' },
};

// Detector — runs in the file picker route to recognize cs_breakpoints_v1 JSONs.
function _isCrossSpeciesJSON(data) {
  return !!(
    data &&
    data.tool === CROSS_SPECIES_TOOL &&
    typeof data.schema_version === 'number' &&
    Array.isArray(data.breakpoints)
  );
}

// State shape:
//   state.crossSpecies = {
//     loaded_at, species_query, species_target, input_paf, params,
//     n_breakpoints, n_by_event_type, breakpoints[],
//     // UI state — not persisted across reloads
//     active_id, filter: { events: Set<string>, search: string },
//     sort: 'gar_chr' | 'block_size' | 'gar_all_te' | 'gar_pos',
//     flank_bp: <number>,
//   }
function _ensureCrossSpeciesState() {
  if (!state.crossSpecies || typeof state.crossSpecies !== 'object') {
    state.crossSpecies = null;
  }
  if (!state._crossSpeciesUI) {
    state._crossSpeciesUI = {
      active_id: null,
      filter: { events: null, search: '' },   // events=null means "all"
      sort: 'gar_pos',
      flank_bp: CROSS_SPECIES_FLANK_DEFAULT_BP,
    };
  }
  return state._crossSpeciesUI;
}

function _storeCrossSpecies(parsed) {
  if (!_isCrossSpeciesJSON(parsed)) return false;
  // Defensive copy of breakpoints — guards against the JSON being mutated
  // by shared-reference paths. Each breakpoint already has serializable
  // shape so a JSON round-trip is the safest deep clone.
  const cleaned = {
    schema_version:  parsed.schema_version,
    tool:            parsed.tool,
    generated_at:    parsed.generated_at || null,
    species_query:   parsed.species_query  || { name: 'unknown', haplotype: '?' },
    species_target:  parsed.species_target || { name: 'unknown', haplotype: '?' },
    input_paf:       parsed.input_paf || null,
    params:          parsed.params || {},
    n_breakpoints:   typeof parsed.n_breakpoints === 'number' ? parsed.n_breakpoints
                     : parsed.breakpoints.length,
    n_by_event_type: parsed.n_by_event_type || {},
    breakpoints:     JSON.parse(JSON.stringify(parsed.breakpoints)),
    loaded_at:       new Date().toISOString(),
    // v4 turn 75: schema v2 fields (synteny_blocks + chrom lengths). Stored
    // as-is when present; absent on v1 JSONs (atlas falls back to "no
    // synteny analysis available" message).
    synteny_blocks:        Array.isArray(parsed.synteny_blocks)
                            ? JSON.parse(JSON.stringify(parsed.synteny_blocks)) : null,
    n_synteny_blocks:      typeof parsed.n_synteny_blocks === 'number'
                            ? parsed.n_synteny_blocks
                            : (Array.isArray(parsed.synteny_blocks) ? parsed.synteny_blocks.length : 0),
    chrom_lengths_query:   parsed.chrom_lengths_query  || {},
    chrom_lengths_target:  parsed.chrom_lengths_target || {},
  };
  state.crossSpecies = cleaned;
  // Invalidate derived synteny caches since the data has changed.
  state._csSyntenyCache = null;
  state._csSyntenyEdgesCache = null;
  state._csInversionContextCache = null;
  // v4 turn 114a: invalidate the per-window overlay index built lazily by
  // _ensureCsOverlayIndex. Cross-page rendering (lines panel, |z| track,
  // sim_mat, future boundaries page) reads from this index. Recomputed
  // on next access.
  state._csOverlayIndex = null;
  return true;
}

function _persistCrossSpecies() {
  try {
    if (state.crossSpecies) {
      localStorage.setItem(CROSS_SPECIES_LS_KEY, JSON.stringify(state.crossSpecies));
    } else {
      localStorage.removeItem(CROSS_SPECIES_LS_KEY);
    }
  } catch (_) { /* fail-soft: quota exceeded just means restore won't work */ }
}

function _restoreCrossSpecies() {
  try {
    const raw = localStorage.getItem(CROSS_SPECIES_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeCrossSpecies(parsed);
  } catch (_) { return false; }
}

function _clearCrossSpecies() {
  state.crossSpecies = null;
  if (state._crossSpeciesUI) state._crossSpeciesUI.active_id = null;
  // Invalidate derived synteny caches (v4 turn 75)
  state._csSyntenyCache = null;
  state._csSyntenyEdgesCache = null;
  state._csInversionContextCache = null;
  // v4 turn 114a: invalidate per-window overlay index too
  state._csOverlayIndex = null;
  _persistCrossSpecies();
}

// =============================================================================
// BLOCK B — cross-species runtime: filter/sort, render*, ideograms, flank charts,
// synteny, dotplot, focal-vs-bg (legacy 23717-26025)
// =============================================================================
// -----------------------------------------------------------------------------
// Filter / sort helpers
// -----------------------------------------------------------------------------

function _csEventTypeOf(bp) {
  // Prefer the refined classification when present
  return bp.event_type_refined || bp.event_type || 'unknown';
}

// v4 turn 114a: cs-breakpoint cross-page overlay index.
// Builds (and caches) a structure that maps each precomp window index to
// the list of cs-breakpoints whose Gar genomic position falls inside that
// window. Used by drawLinesPanel, drawZ, drawSim, and future cross-page
// renderers. Invalidated when state.crossSpecies changes (see _storeCrossSpecies).
//
// Output shape:
//   {
//     chrom: state.data.chrom,
//     // Only breakpoints whose gar_chr matches the currently-loaded chrom
//     bps: [{ id, win, mb, event_type, ... }, ...],
//     byWindow: Map<windowIdx, [bp_ref, ...]>,
//   }
//
// Returns null if no crossSpecies data, no precomp data, or chrom mismatch.
function _ensureCsOverlayIndex() {
  if (state._csOverlayIndex) return state._csOverlayIndex;
  if (!state.data || !state.crossSpecies) return null;
  const cs = state.crossSpecies;
  if (!Array.isArray(cs.breakpoints) || cs.breakpoints.length === 0) return null;
  const chrom = state.data.chrom;
  if (!chrom) return null;
  // Match cs-breakpoints to the loaded chrom. We tolerate exact match or
  // a case-insensitive match on the gar_chr field (the precomp chrom is
  // always Gar reference per project convention).
  const norm = (s) => (typeof s === 'string') ? s.trim() : '';
  const want = norm(chrom).toLowerCase();
  const wins = state.data.windows;
  if (!Array.isArray(wins) || wins.length === 0) return null;

  // Pre-compute a sorted ascending array of window start_bp for binary search
  // (windows are already in genomic order in standard precomp output, but we
  // don't assume — we re-derive). Each window is [start_bp, end_bp); some
  // older precomps use start_bp / end_bp keys, others use start / end.
  const winStart = new Float64Array(wins.length);
  const winEnd   = new Float64Array(wins.length);
  for (let i = 0; i < wins.length; i++) {
    const w = wins[i];
    winStart[i] = (typeof w.start_bp === 'number') ? w.start_bp
                : (typeof w.start === 'number') ? w.start : 0;
    winEnd[i]   = (typeof w.end_bp === 'number') ? w.end_bp
                : (typeof w.end === 'number') ? w.end
                : (winStart[i] + 50000);  // fall-back 50 kb step if absent
  }
  // Binary search: returns the window index that contains pos, or -1 if pos
  // is outside the chrom range. Searches for the largest i with winStart[i]
  // <= pos AND winEnd[i] > pos. If pos beyond last window's end_bp, returns -1.
  const _bpToWindow = (pos) => {
    if (!Number.isFinite(pos)) return -1;
    if (pos < winStart[0]) return -1;
    if (pos >= winEnd[winEnd.length - 1]) return -1;
    let lo = 0, hi = winStart.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >>> 1;
      if (pos < winStart[mid]) {
        hi = mid - 1;
      } else if (pos >= winEnd[mid]) {
        lo = mid + 1;
      } else {
        return mid;
      }
    }
    // Edge case: position fell into a gap between windows. Return the
    // closest preceding window so the overlay still shows roughly there.
    return Math.max(0, Math.min(winStart.length - 1, lo - 1));
  };

  const out = { chrom, bps: [], byWindow: new Map(), bpToWindow: _bpToWindow };
  for (const bp of cs.breakpoints) {
    if (norm(bp.gar_chr).toLowerCase() !== want) continue;
    // Use gar_pos_start as primary anchor (start of breakpoint interval).
    // gar_pos_mb is also available; fall back to it × 1e6 if needed.
    let pos = (typeof bp.gar_pos_start === 'number') ? bp.gar_pos_start
            : (typeof bp.gar_pos_bp === 'number') ? bp.gar_pos_bp
            : (typeof bp.gar_pos_mb === 'number') ? bp.gar_pos_mb * 1e6
            : NaN;
    const wi = _bpToWindow(pos);
    if (wi < 0) continue;
    const entry = {
      id: bp.id,
      win: wi,
      pos_bp: pos,
      mb: pos / 1e6,
      event_type: _csEventTypeOf(bp),
      raw: bp,
    };
    out.bps.push(entry);
    if (!out.byWindow.has(wi)) out.byWindow.set(wi, []);
    out.byWindow.get(wi).push(entry);
  }
  state._csOverlayIndex = out;
  return out;
}

// Quick accessor (returns array, possibly empty) for a given window index.
function _csOverlayBpsForWindow(wi) {
  const idx = _ensureCsOverlayIndex();
  if (!idx) return [];
  return idx.byWindow.get(wi) || [];
}

// v4 turn 114c (remaining): cross-page click-to-jump helpers for cs-bp lines.
// Used by zCanvas / sim_mat / lines panel / winSumStrip / boundaries TE
// density. Caller is responsible for the canvas-specific projection from
// each bp to a target screen coordinate; this just walks the bps array,
// computes squared distance, and returns the closest hit within tolerance.
//
// Hit tolerance is intentionally generous (4 px ≈ the dashed line's visual
// "thickness" with the 1.2 px stroke + 4-on-3-off dash pattern) so users
// don't have to land pixel-perfect on a thin dashed line.
const _CS_BP_HIT_TOL_PX = 4;

// Linear hit-test for a 1-D vertical-line layout. Caller passes a function
// `bpToX(bp) -> number | null` that maps each bp to the same x position the
// renderer drew at; bps outside the visible range can return null and are
// skipped. Returns the closest bp within tolerance, or null.
function _csBpHitTestX(clickX, tolerancePx) {
  const idx = _ensureCsOverlayIndex();
  if (!idx || idx.bps.length === 0) return null;
  const tol = (tolerancePx != null) ? tolerancePx : _CS_BP_HIT_TOL_PX;
  return _csBpHitTestXFromList(idx.bps, clickX, tolerancePx, arguments[2]);
}
// Internal — split out so per-canvas callers can pass their own bpToX
// projection (different canvases use different mb→x mappings) and an
// already-filtered subset of bps if they want.
function _csBpHitTestXFromList(bps, clickX, tolerancePx, bpToX) {
  if (!Array.isArray(bps) || bps.length === 0) return null;
  if (typeof bpToX !== 'function') return null;
  const tol = (tolerancePx != null) ? tolerancePx : _CS_BP_HIT_TOL_PX;
  let bestBp = null, bestDx = tol + 0.001;
  for (const bp of bps) {
    const xx = bpToX(bp);
    if (!Number.isFinite(xx)) continue;
    const dx = Math.abs(xx - clickX);
    if (dx < bestDx) { bestDx = dx; bestBp = bp; }
  }
  return bestBp;
}

// 2-D hit-test for the sim_mat red cross overlay. Each cs-bp renders at
// (toPx(bp.win), toPy(bp.win)) on the diagonal; this returns the closest
// bp within tolerance of the click in 2-D.
function _csBpHitTest2D(bps, clickX, clickY, bpToXY, tolerancePx) {
  if (!Array.isArray(bps) || bps.length === 0) return null;
  if (typeof bpToXY !== 'function') return null;
  const tol = (tolerancePx != null) ? tolerancePx : _CS_BP_HIT_TOL_PX;
  const tol2 = tol * tol;
  let bestBp = null, bestD2 = tol2 + 0.001;
  for (const bp of bps) {
    const xy = bpToXY(bp);
    if (!xy || !Number.isFinite(xy.x) || !Number.isFinite(xy.y)) continue;
    const dx = xy.x - clickX, dy = xy.y - clickY;
    const d2 = dx * dx + dy * dy;
    if (d2 < bestD2) { bestD2 = d2; bestBp = bp; }
  }
  return bestBp;
}

// Common jump action — center cursor on the breakpoint's window and trigger
// the standard cross-panel redraw cascade. setCur already handles the panel
// invalidation; we just wrap it for symmetric behavior across canvases.
function _csBpJumpToWindow(bp) {
  if (!bp || typeof bp.win !== 'number' || bp.win < 0) return false;
  if (typeof setCur === 'function') {
    setCur(bp.win);
    return true;
  }
  // Fallback: mutate state directly + fire a minimal redraw.
  state.cur = bp.win;
  if (typeof drawZ === 'function') { try { drawZ(); } catch (_) {} }
  if (typeof drawSim === 'function') { try { drawSim(); } catch (_) {} }
  if (typeof drawLinesPanel === 'function') { try { drawLinesPanel(); } catch (_) {} }
  if (typeof drawWinSumStrip === 'function') { try { drawWinSumStrip(); } catch (_) {} }
  return true;
}
// Expose for tests
if (typeof window !== 'undefined') {
  window._csBpHitTestXFromList = _csBpHitTestXFromList;
  window._csBpHitTest2D = _csBpHitTest2D;
  window._csBpJumpToWindow = _csBpJumpToWindow;
}

// =============================================================================
// v4 turn 114d — cross-panel hover glow on cs-breakpoint proximity.
// =============================================================================
// When the cursor moves close to a cs-breakpoint line on any source panel,
// the PCA scatter (#pcaPanel) and the L3 contingency focal column
// (#l3FocalCol) light up with a soft blue glow via the .cs-bp-hover-glow
// CSS class. Pure class-toggle, zero canvas redraws on the target panels.
//
// Source panels:
//   - zCanvas, simCanvas, winSumStripCanvas (single canvases)
//   - lines-subpanel canvases (multiple, wired during drawLinesPanel's loop)
//   - SVG hit-zone rects on the boundaries TE-density panel (use native
//     mouseenter/mouseleave on the data-cs-bp-id rects)
//
// Hit detection on canvas panels uses _csBpHitTestXFromList /
// _csBpHitTest2D — the same helpers as the click-to-jump path. mousemove
// is throttled via requestAnimationFrame so we run hit-testing at most
// once per frame regardless of mousemove fire rate.
//
// State:
//   state._csHoverActive  — currently-hovered bp.id, or null
//   state._csHoverRaf     — pending rAF handle, used for throttling
//
// We expose _csBpHoverEnter / _csBpHoverLeave / _wireCsBpHoverOnCanvas on
// window for tests + for drawLinesPanel's per-subpanel wiring.

const _CS_BP_HOVER_TARGETS = ['pcaPanel', 'l3FocalCol'];
const _CS_BP_HOVER_TOL_PX = 5;   // 1 px more generous than the click tol

function _csBpHoverEnter(bp) {
  if (!bp) return;
  const bpId = (bp.id != null) ? String(bp.id) : ('w' + bp.win);
  if (state._csHoverActive === bpId) return;   // already in this state
  state._csHoverActive = bpId;
  for (const id of _CS_BP_HOVER_TARGETS) {
    const el = document.getElementById(id);
    if (el) el.classList.add('cs-bp-hover-glow');
  }
}

function _csBpHoverLeave() {
  if (state._csHoverActive == null) return;   // already in this state
  state._csHoverActive = null;
  for (const id of _CS_BP_HOVER_TARGETS) {
    const el = document.getElementById(id);
    if (el) el.classList.remove('cs-bp-hover-glow');
  }
}

// Helper used by drawLinesPanel's subpanel loop and by the one-shot wiring
// for #zCanvas / #simCanvas / #winSumStripCanvas. The mode argument picks
// between linear (1-D x-axis) and diagonal (sim_mat 2-D) hit-testing.
//
// projection: function called with (canvas, evt) → { hitTest(bps) → bp|null }
// where hitTest does the per-canvas geometry lookup (mb→x, win→{x,y}, etc.).
//
// Returns nothing; attaches mousemove + pointerleave listeners to the canvas.
// Idempotent via canvas.__csBpHoverWired flag.
function _wireCsBpHoverOnCanvas(canvas, projection) {
  if (!canvas || canvas.__csBpHoverWired) return;
  canvas.__csBpHoverWired = true;
  canvas.addEventListener('mousemove', (evt) => {
    // Throttle via rAF — at most one hit-test per frame.
    if (state._csHoverRaf) return;
    state._csHoverRaf = requestAnimationFrame(() => {
      state._csHoverRaf = null;
      try {
        const csIdx = (typeof _ensureCsOverlayIndex === 'function')
          ? _ensureCsOverlayIndex() : null;
        if (!csIdx || csIdx.bps.length === 0) {
          if (state._csHoverActive != null) _csBpHoverLeave();
          return;
        }
        const hit = projection(canvas, evt, csIdx);
        if (hit) _csBpHoverEnter(hit);
        else     _csBpHoverLeave();
      } catch (_) { /* fail-soft */ }
    });
  });
  // Leaving the canvas should always clear, even if rAF would normally fire
  // a final hit-test. Run it synchronously here for snappy fade-out.
  canvas.addEventListener('pointerleave', () => {
    if (state._csHoverRaf) {
      cancelAnimationFrame(state._csHoverRaf);
      state._csHoverRaf = null;
    }
    _csBpHoverLeave();
  });
}

// Expose for tests + for drawLinesPanel to call per subpanel
if (typeof window !== 'undefined') {
  window._csBpHoverEnter = _csBpHoverEnter;
  window._csBpHoverLeave = _csBpHoverLeave;
  window._wireCsBpHoverOnCanvas = _wireCsBpHoverOnCanvas;
  window._CS_BP_HOVER_TOL_PX = _CS_BP_HOVER_TOL_PX;
}

function _csGarAllTeFlankMean(bp) {
  return ((bp.flanking_repeat_density_gar || {}).all_TE || {}).mean;
}

function _csIsSpalaxFlagged(bp) {
  // Either an explicit manuscript_note OR a derived check using the loaded
  // chrom-wide all_TE mean (when the boundaries-page TEfull layer has it).
  if (bp.manuscript_note) return true;
  const flank = _csGarAllTeFlankMean(bp);
  if (typeof flank !== 'number') return false;
  const chrom = bp.gar_chr;
  const entry = (typeof _getRepeatDensity === 'function') ? _getRepeatDensity(chrom) : null;
  if (!entry || !entry.chrom_block || !entry.chrom_block.by_class) return false;
  const dens = (entry.chrom_block.by_class.all_TE || {}).densities || [];
  if (dens.length === 0) return false;
  const sum = dens.reduce((a, d) => a + (Number.isFinite(d) ? d : 0), 0);
  const cnt = dens.reduce((a, d) => a + (Number.isFinite(d) ? 1 : 0), 0);
  if (cnt === 0) return false;
  const chromMean = sum / cnt;
  return flank >= 1.5 * chromMean;
}

function _csFilteredSorted() {
  if (!state.crossSpecies) return [];
  const ui = _ensureCrossSpeciesState();
  let bps = state.crossSpecies.breakpoints.slice();
  // Event filter
  if (ui.filter.events && ui.filter.events.size > 0) {
    bps = bps.filter(b => ui.filter.events.has(_csEventTypeOf(b)));
  }
  // Search (case-insensitive substring match against gar_chr + mac_chrs + id)
  if (ui.filter.search) {
    const q = ui.filter.search.trim().toLowerCase();
    if (q) {
      bps = bps.filter(b => {
        const hay = [
          b.id, b.gar_chr,
          (b.prev_block || {}).mac_chr,
          (b.next_block || {}).mac_chr,
        ].filter(Boolean).join(' ').toLowerCase();
        return hay.includes(q);
      });
    }
  }
  // Sort
  const cmp = (a, b) => {
    if (ui.sort === 'gar_pos') {
      const c = String(a.gar_chr).localeCompare(String(b.gar_chr));
      if (c) return c;
      return (a.gar_pos_start || 0) - (b.gar_pos_start || 0);
    }
    if (ui.sort === 'block_size') {
      const sa = ((a.prev_block || {}).block_size_bp || 0) + ((a.next_block || {}).block_size_bp || 0);
      const sb = ((b.prev_block || {}).block_size_bp || 0) + ((b.next_block || {}).block_size_bp || 0);
      return sb - sa;
    }
    if (ui.sort === 'gar_all_te') {
      const fa = _csGarAllTeFlankMean(a);
      const fb = _csGarAllTeFlankMean(b);
      // null-pushed-to-bottom
      if (typeof fa !== 'number' && typeof fb !== 'number') return 0;
      if (typeof fa !== 'number') return 1;
      if (typeof fb !== 'number') return -1;
      return fb - fa;   // descending
    }
    if (ui.sort === 'gar_chr') {
      return String(a.gar_chr).localeCompare(String(b.gar_chr));
    }
    return 0;
  };
  bps.sort(cmp);
  return bps;
}

// -----------------------------------------------------------------------------
// Renderer: toolbar
// -----------------------------------------------------------------------------

function _renderCrossSpeciesToolbar() {
  const slot = document.getElementById('csToolbar');
  if (!slot) return;
  const ui = _ensureCrossSpeciesState();
  const cs = state.crossSpecies;
  if (!cs) {
    slot.innerHTML = '<span class="cs-toolbar-meta" style="color: var(--ink-dimmer);">' +
                     'Drop a <code style="color:var(--ink-dim);">cs_breakpoints_v1.json</code> file to populate.' +
                     '</span>';
    return;
  }
  // Per-event-type counts (post-filter for the active ones, post-search-only for chip totals)
  const counts = {};
  for (const b of cs.breakpoints) {
    const et = _csEventTypeOf(b);
    counts[et] = (counts[et] || 0) + 1;
  }
  // Chip order: known events first (in CS_EVENT_DEF order), then any others
  const known = Object.keys(CS_EVENT_DEF);
  const others = Object.keys(counts).filter(k => !known.includes(k));
  const order = known.filter(k => counts[k] > 0).concat(others);
  const chips = order.map(et => {
    const def = CS_EVENT_DEF[et] || { label: et, cls: 'cs-evt-mixed', icon: '\u25CF' };
    const on  = !!(ui.filter.events && ui.filter.events.has(et));
    return '<span class="cs-chip ' + (on ? 'cs-chip-on' : '') +
           '" data-cs-evt="' + _esc(et) + '" title="Filter to ' + _esc(def.label) +
           ' breakpoints (click to toggle)">' +
           '<span style="opacity:0.85;">' + def.icon + '</span> ' +
           _esc(def.label) +
           '<span class="cs-chip-count">' + counts[et] + '</span>' +
           '</span>';
  }).join('');
  // Sort dropdown
  const sortOpts = [
    ['gar_pos',     'Gar pos'],
    ['gar_chr',     'Gar chrom'],
    ['block_size',  'block size'],
    ['gar_all_te',  'gar all_TE flank'],
  ].map(([v, lbl]) =>
    '<option value="' + v + '"' + (ui.sort === v ? ' selected' : '') + '>' + _esc(lbl) + '</option>'
  ).join('');
  // Visible count
  const visible = _csFilteredSorted().length;
  const total   = cs.breakpoints.length;
  // Assemble
  const sq = (cs.species_query || {}).name || '?';
  const st = (cs.species_target || {}).name || '?';
  slot.innerHTML =
    '<span class="cs-toolbar-meta"><b>' + _esc(sq) + '</b> \u2194 <b>' + _esc(st) + '</b></span>' +
    chips +
    '<input type="search" id="csSearch" class="cs-search" placeholder="search chrom or id"' +
    ' value="' + _esc(ui.filter.search || '') + '" />' +
    '<label class="cs-toolbar-meta">sort ' +
      '<select id="csSortSelect" class="cs-sort">' + sortOpts + '</select>' +
    '</label>' +
    '<span class="cs-toolbar-spacer"></span>' +
    '<span class="cs-toolbar-meta"><b>' + visible + '</b> / ' + total + ' breakpoints</span>' +
    (cs.generated_at
      ? '<span class="cs-toolbar-meta" style="font-size:10px;opacity:0.7;">generated ' + _esc(cs.generated_at) + '</span>'
      : '') +
    '<button id="csClearBtn" class="cs-chip" title="Forget the loaded cross-species data (does not delete the file on disk)">' +
      'forget</button>';
  // Wire interactions
  slot.querySelectorAll('[data-cs-evt]').forEach(el => {
    el.addEventListener('click', () => {
      const et = el.getAttribute('data-cs-evt');
      if (!ui.filter.events) ui.filter.events = new Set();
      if (ui.filter.events.has(et)) ui.filter.events.delete(et);
      else ui.filter.events.add(et);
      // empty set behaves like "all"
      if (ui.filter.events.size === 0) ui.filter.events = null;
      _renderCrossSpeciesToolbar();
      _renderCrossSpeciesCatalogue();
    });
  });
  const search = slot.querySelector('#csSearch');
  if (search) {
    search.addEventListener('input', () => {
      ui.filter.search = search.value || '';
      _renderCrossSpeciesCatalogue();
    });
  }
  const sortSel = slot.querySelector('#csSortSelect');
  if (sortSel) {
    sortSel.addEventListener('change', () => {
      ui.sort = sortSel.value;
      _renderCrossSpeciesCatalogue();
    });
  }
  const clearBtn = slot.querySelector('#csClearBtn');
  if (clearBtn) {
    clearBtn.addEventListener('click', () => {
      if (!confirm('Forget the loaded cross-species breakpoint catalogue? ' +
                   'This only clears the atlas state — the JSON file on disk is untouched.')) return;
      _clearCrossSpecies();
      _renderCrossSpeciesPage();
    });
  }
}

// -----------------------------------------------------------------------------
// Renderer: catalogue list
// -----------------------------------------------------------------------------

function _renderCrossSpeciesCatalogue() {
  const slot = document.getElementById('csCatalogue');
  if (!slot) return;
  const ui = _ensureCrossSpeciesState();
  const bps = _csFilteredSorted();
  if (bps.length === 0) {
    slot.innerHTML = '<div class="cs-cat-empty">' +
      (state.crossSpecies
        ? 'No breakpoints match the current filter.'
        : 'No data loaded yet.') +
      '</div>';
    return;
  }
  let html = '';
  for (const bp of bps) {
    const et = _csEventTypeOf(bp);
    const def = CS_EVENT_DEF[et] || { label: et, cls: 'cs-evt-mixed', icon: '\u25CF' };
    const isActive = (ui.active_id === bp.id);
    const flagSpalax = _csIsSpalaxFlagged(bp);
    const hasCand    = Array.isArray(bp.candidate_overlap) && bp.candidate_overlap.length > 0;
    const flanksGar  = _csGarAllTeFlankMean(bp);
    const flanksGarStr = (typeof flanksGar === 'number')
      ? ' \u00B7 gar all_TE ' + flanksGar.toFixed(2)
      : '';
    const macChrs = [];
    if (bp.prev_block && bp.prev_block.mac_chr) macChrs.push(bp.prev_block.mac_chr);
    if (bp.next_block && bp.next_block.mac_chr && bp.next_block.mac_chr !== bp.prev_block?.mac_chr) {
      macChrs.push(bp.next_block.mac_chr);
    }
    const macStr = macChrs.join(' \u2192 ') || '?';
    const flagStr =
      (flagSpalax ? '<span class="cs-row-flag cs-flag-spalax" title="Spalax-style all_TE enrichment at gar flank (\u22651.5x chrom mean)">\u25CF</span>' : '') +
      (hasCand    ? '<span class="cs-row-flag cs-flag-cand" title="Overlaps polymorphic candidate(s): ' +
                    _esc(bp.candidate_overlap.join(', ')) + '">\u272A</span>' : '');
    html +=
      '<div class="cs-row ' + (isActive ? 'cs-row-active' : '') +
      '" data-cs-bp-id="' + _esc(bp.id) + '" title="' + _esc(bp.id) + '">' +
        '<span class="cs-row-icon ' + def.cls + '">' + def.icon + '</span>' +
        '<span class="cs-row-main">' +
          '<div class="cs-row-main-line1">' +
            '<code>' + _esc(bp.gar_chr) + '</code> ' +
            '<span class="cs-arrow">\u2192</span> ' +
            '<code>' + _esc(macStr) + '</code>' +
          '</div>' +
          '<div class="cs-row-main-line2">' +
            _esc(def.label) + ' \u00B7 ' + (bp.gar_pos_mb || 0).toFixed(2) + ' Mb' +
            flanksGarStr +
          '</div>' +
        '</span>' +
        '<span>' + flagStr + '</span>' +
      '</div>';
  }
  slot.innerHTML = html;
  // Wire row-click
  slot.querySelectorAll('[data-cs-bp-id]').forEach(el => {
    el.addEventListener('click', () => {
      const id = el.getAttribute('data-cs-bp-id');
      ui.active_id = id;
      _renderCrossSpeciesCatalogue();
      _renderCrossSpeciesFocus();
    });
  });
}

// -----------------------------------------------------------------------------
// Renderer: focus panel (header + ideogram + flank charts)
// -----------------------------------------------------------------------------

function _csFindActive() {
  const ui = _ensureCrossSpeciesState();
  if (!state.crossSpecies || !ui.active_id) return null;
  return state.crossSpecies.breakpoints.find(b => b.id === ui.active_id) || null;
}

function _renderCrossSpeciesFocus() {
  const empty   = document.getElementById('csEmpty');
  const header  = document.getElementById('csFocusHeader');
  const ideos   = document.getElementById('csIdeograms');
  const flanks  = document.getElementById('csFlankCharts');
  if (!empty || !header || !ideos || !flanks) return;
  const bp = _csFindActive();
  if (!bp) {
    // Show empty state when no breakpoint is selected.
    empty.style.display  = state.crossSpecies ? 'none' : 'block';
    header.style.display = 'none';
    ideos.style.display  = 'none';
    flanks.style.display = 'none';
    if (state.crossSpecies) {
      // Data is loaded but nothing selected yet — show a hint
      header.style.display = 'block';
      header.innerHTML =
        '<div class="cs-focus-title" style="color:var(--ink-dim);">' +
          'Select a breakpoint from the catalogue on the left to inspect.' +
        '</div>';
    }
    return;
  }
  empty.style.display = 'none';
  header.style.display = 'block';
  ideos.style.display  = 'block';
  flanks.style.display = 'block';
  // Header
  const et = _csEventTypeOf(bp);
  const def = CS_EVENT_DEF[et] || { label: et, cls: 'cs-evt-mixed', icon: '\u25CF' };
  const sq = (state.crossSpecies.species_query || {}).name || 'query';
  const st = (state.crossSpecies.species_target || {}).name || 'target';
  const note = bp.manuscript_note
    ? '<div class="cs-focus-note">\u25CF ' + _esc(bp.manuscript_note) + '</div>'
    : '';
  const candStr = (Array.isArray(bp.candidate_overlap) && bp.candidate_overlap.length > 0)
    ? '<div class="cs-focus-note-cand">\u272A coincides with polymorphic candidate(s): <code>' +
      bp.candidate_overlap.map(c => _esc(c)).join('</code>, <code>') + '</code>' +
      ' \u2014 the same locus is recurrently rearrangement-prone across species AND segregating in the 226-sample cohort.</div>'
    : '';
  const memberStr = bp.n_member_breakpoints && bp.n_member_breakpoints > 1
    ? ' \u00B7 ' + bp.n_member_breakpoints + ' member breakpoints' : '';
  const prevStrand = (bp.prev_block || {}).strand || '?';
  const nextStrand = (bp.next_block || {}).strand || '?';
  header.innerHTML =
    '<div class="cs-focus-title">' +
      '<code>' + _esc(bp.id) + '</code>' +
      '<span class="cs-evt-pill ' + def.cls + '">' + _esc(def.label) + '</span>' +
      memberStr +
    '</div>' +
    '<div class="cs-focus-coords">' +
      '<i>' + _esc(sq) + '</i> <code>' + _esc(bp.gar_chr) + '</code> @ ' +
      '<code>' + (bp.gar_pos_mb || 0).toFixed(3) + ' Mb</code>' +
      ' \u2014 prev block: <code>' + _esc((bp.prev_block || {}).mac_chr || '?') + '</code> ' + prevStrand +
      ', next block: <code>' + _esc((bp.next_block || {}).mac_chr || '?') + '</code> ' + nextStrand +
    '</div>' +
    note + candStr;
  // Ideograms
  ideos.innerHTML = '<div class="cs-ideograms">' + _csBuildIdeogramsSvg(bp) + '</div>';
  // Flank charts: per-species TE-density tracks + bar-table
  flanks.innerHTML = _csBuildFlankCharts(bp);
  // turn 117: focal-vs-background widget. Wired here (rather than in
  // _renderCrossSpeciesPage) because the catalogue row-click + ↑/↓ step
  // handlers both call _renderCrossSpeciesFocus directly, not the page
  // renderer — so this is the single hook that catches all active_id
  // changes. Fail-soft.
  if (typeof _renderCrossSpeciesFocalVsBg === 'function') {
    try { _renderCrossSpeciesFocalVsBg(); } catch (e) {
      console.warn('[crossSpecies] focal-vs-bg render failed:',
                   e && e.message ? e.message : e);
    }
  }
}

// -----------------------------------------------------------------------------
// Ideogram SVG: two parallel chrom strips (Gar above, Mac below) with the
// breakpoint marker on Gar and the prev/next blocks' Mac coordinates marked
// on the Mac strip (one or two strips depending on whether prev.mac_chr ===
// next.mac_chr). A connecting line shows which Mac region maps to the Gar
// breakpoint flank.
// -----------------------------------------------------------------------------

function _csIdeogramChromLength(species, chrom) {
  // Try to get the chrom length from the loaded TEfull JSON's window_end_bp[-1]
  // — that's a good approximation. species is 'gar' or 'mac'; chrom is the
  // chrom name string. Returns a number (bp) or null.
  if (!chrom) return null;
  const entry = (typeof _getRepeatDensity === 'function') ? _getRepeatDensity(chrom) : null;
  if (entry && entry.chrom_block && Array.isArray(entry.chrom_block.window_end_bp)) {
    const arr = entry.chrom_block.window_end_bp;
    if (arr.length > 0) return arr[arr.length - 1];
  }
  // Fall back: if Gar, use the active state.data n_bp when chrom matches
  if (species === 'gar' && state.data && state.data.chrom === chrom && state.data.n_bp) {
    return state.data.n_bp;
  }
  return null;
}

function _csBuildIdeogramsSvg(bp) {
  const W = 760;            // overall width (the panel pads internally)
  const padL = 92, padR = 24;
  const padT = 22, padB = 28;
  const rowH = 18;
  const rowGap = 60;
  const innerW = W - padL - padR;
  const garLen = _csIdeogramChromLength('gar', bp.gar_chr);
  const macPrevChr = (bp.prev_block || {}).mac_chr;
  const macNextChr = (bp.next_block || {}).mac_chr;
  const sameMac = macPrevChr && macNextChr && macPrevChr === macNextChr;
  const macLen1 = _csIdeogramChromLength('mac', macPrevChr);
  const macLen2 = sameMac ? null : _csIdeogramChromLength('mac', macNextChr);
  const H = padT + rowH + rowGap + rowH + (sameMac ? 0 : (rowGap + rowH)) + padB;

  // Helpers
  const xOnRow = (pos, len) => {
    if (!len || !Number.isFinite(len) || len <= 0) return null;
    return padL + (pos / len) * innerW;
  };
  const formatBp = (bp_) => {
    if (bp_ == null) return '?';
    if (bp_ >= 1e6) return (bp_ / 1e6).toFixed(2) + ' Mb';
    if (bp_ >= 1e3) return (bp_ / 1e3).toFixed(0) + ' kb';
    return bp_ + ' bp';
  };
  const drawRow = (yTop, label, lenBp, color) => {
    if (!lenBp) {
      // Fallback: dashed unknown-length strip
      return '<rect x="' + padL + '" y="' + yTop + '" width="' + innerW + '" height="' + rowH +
             '" fill="none" stroke="' + color + '" stroke-width="1" stroke-dasharray="3 3" rx="3" />' +
             '<text x="' + (padL - 8) + '" y="' + (yTop + rowH * 0.7) + '" font-size="10.5" ' +
             'fill="' + color + '" text-anchor="end" font-family="ui-monospace, monospace">' +
             _esc(label) + '</text>' +
             '<text x="' + (padL + innerW + 4) + '" y="' + (yTop + rowH * 0.7) + '" ' +
             'font-size="9.5" fill="var(--ink-dimmer)" font-family="ui-monospace, monospace">no length</text>';
    }
    return '<rect x="' + padL + '" y="' + yTop + '" width="' + innerW + '" height="' + rowH +
           '" fill="' + color + '" fill-opacity="0.18" stroke="' + color + '" stroke-width="1" rx="3" />' +
           '<text x="' + (padL - 8) + '" y="' + (yTop + rowH * 0.7) + '" font-size="10.5" ' +
           'fill="' + color + '" text-anchor="end" font-family="ui-monospace, monospace">' +
           _esc(label) + '</text>' +
           '<text x="' + (padL + innerW + 4) + '" y="' + (yTop + rowH * 0.7) + '" ' +
           'font-size="9.5" fill="var(--ink-dimmer)" font-family="ui-monospace, monospace">' +
           formatBp(lenBp) + '</text>';
  };

  // Row 1: Gar
  const garY = padT;
  let svg = '';
  svg += drawRow(garY, bp.gar_chr, garLen, '#4fa3ff');
  // Breakpoint marker on Gar (vertical line + diamond)
  if (garLen) {
    const gx = xOnRow((bp.gar_pos_start + bp.gar_pos_end) / 2, garLen);
    const gxs = xOnRow(bp.gar_pos_start, garLen);
    const gxe = xOnRow(bp.gar_pos_end, garLen);
    if (gxs != null && gxe != null) {
      svg +=
        '<rect x="' + Math.max(padL, gxs - 0.5) + '" y="' + (garY - 2) +
        '" width="' + Math.max(2, gxe - gxs + 1) + '" height="' + (rowH + 4) +
        '" fill="var(--accent)" fill-opacity="0.45" />';
    }
    if (gx != null) {
      svg +=
        '<line x1="' + gx.toFixed(1) + '" x2="' + gx.toFixed(1) +
        '" y1="' + (garY - 6) + '" y2="' + (garY + rowH + 6) + '" ' +
        'stroke="var(--accent)" stroke-width="1.4" />' +
        '<polygon points="' + gx.toFixed(1) + ',' + (garY - 9) + ' ' +
                              (gx + 4).toFixed(1) + ',' + (garY - 5) + ' ' +
                              gx.toFixed(1) + ',' + (garY - 1) + ' ' +
                              (gx - 4).toFixed(1) + ',' + (garY - 5) + '" ' +
        'fill="var(--accent)" />' +
        '<text x="' + gx.toFixed(1) + '" y="' + (garY - 12) + '" ' +
        'font-size="10" fill="var(--accent)" text-anchor="middle" font-family="ui-monospace, monospace">' +
        (bp.gar_pos_mb || 0).toFixed(2) + ' Mb</text>';
    }
  }

  // Row 2: Mac (prev block) — same row as 'next block' if sameMac, else two rows
  const mac1Y = garY + rowH + rowGap;
  if (sameMac) {
    svg += drawRow(mac1Y, macPrevChr, macLen1, '#3cc08a');
  } else {
    svg += drawRow(mac1Y, macPrevChr || '?', macLen1, '#3cc08a');
  }

  // Markers + connecting lines for the prev block on Mac row 1
  const drawBlockMarker = (yTop, lenBp, startBp, endBp, strand, color) => {
    if (!lenBp) return '';
    const x1 = xOnRow(startBp, lenBp);
    const x2 = xOnRow(endBp,   lenBp);
    if (x1 == null || x2 == null) return '';
    const xa = Math.min(x1, x2);
    const xb = Math.max(x1, x2);
    let m =
      '<rect x="' + xa.toFixed(1) + '" y="' + (yTop - 2) + '" width="' +
      Math.max(2, xb - xa).toFixed(1) + '" height="' + (rowH + 4) +
      '" fill="' + color + '" fill-opacity="0.55" stroke="' + color + '" stroke-width="0.8" />';
    // Strand arrow
    const cx = (xa + xb) / 2;
    const arrowChar = (strand === '-') ? '\u2190' : '\u2192';
    m +=
      '<text x="' + cx.toFixed(1) + '" y="' + (yTop + rowH * 0.7) + '" ' +
      'font-size="11" fill="white" text-anchor="middle" font-family="ui-monospace, monospace" font-weight="600">' +
      arrowChar + '</text>';
    return m;
  };
  // prev block marker on Mac row 1
  svg += drawBlockMarker(
    mac1Y, macLen1,
    (bp.prev_block || {}).mac_start_bp || 0,
    (bp.prev_block || {}).mac_end_bp || 0,
    (bp.prev_block || {}).strand,
    '#3cc08a',
  );
  // next block: either on same Mac row (sameMac) or new row
  let mac2Y = null;
  if (sameMac) {
    svg += drawBlockMarker(
      mac1Y, macLen1,
      (bp.next_block || {}).mac_start_bp || 0,
      (bp.next_block || {}).mac_end_bp || 0,
      (bp.next_block || {}).strand,
      '#9b59b6',   // distinct color for the next block on same chrom
    );
  } else {
    mac2Y = mac1Y + rowH + rowGap;
    svg += drawRow(mac2Y, macNextChr || '?', macLen2, '#9b59b6');
    svg += drawBlockMarker(
      mac2Y, macLen2,
      (bp.next_block || {}).mac_start_bp || 0,
      (bp.next_block || {}).mac_end_bp || 0,
      (bp.next_block || {}).strand,
      '#9b59b6',
    );
  }

  // Connecting lines: from Gar breakpoint anchor to the inner edges of the
  // prev/next blocks on the Mac rows (i.e. the block-boundary closest to
  // the gar breakpoint). Curved bezier for visual separation.
  if (garLen) {
    const gx = xOnRow((bp.gar_pos_start + bp.gar_pos_end) / 2, garLen);
    if (gx != null) {
      // prev block: anchor at mac_end (block ends before the bp)
      if (macLen1) {
        const mPrevX = xOnRow((bp.prev_block || {}).mac_end_bp || 0, macLen1);
        if (mPrevX != null) {
          const y1 = garY + rowH;
          const y2 = mac1Y;
          const cy = (y1 + y2) / 2;
          svg +=
            '<path d="M ' + gx.toFixed(1) + ' ' + (y1 + 2) +
            ' C ' + gx.toFixed(1) + ' ' + cy + ', ' + mPrevX.toFixed(1) + ' ' + cy + ', ' +
            mPrevX.toFixed(1) + ' ' + (y2 - 2) + '" ' +
            'fill="none" stroke="#3cc08a" stroke-width="1" stroke-opacity="0.55" stroke-dasharray="2 2" />';
        }
      }
      // next block: anchor at mac_start (block starts after the bp)
      const macLenForNext = sameMac ? macLen1 : macLen2;
      const y2_next = sameMac ? mac1Y : mac2Y;
      if (macLenForNext && y2_next != null) {
        const mNextX = xOnRow((bp.next_block || {}).mac_start_bp || 0, macLenForNext);
        if (mNextX != null) {
          const y1 = garY + rowH;
          const y2 = y2_next;
          const cy = (y1 + y2) / 2;
          svg +=
            '<path d="M ' + gx.toFixed(1) + ' ' + (y1 + 2) +
            ' C ' + gx.toFixed(1) + ' ' + cy + ', ' + mNextX.toFixed(1) + ' ' + cy + ', ' +
            mNextX.toFixed(1) + ' ' + (y2 - 2) + '" ' +
            'fill="none" stroke="#9b59b6" stroke-width="1" stroke-opacity="0.55" stroke-dasharray="2 2" />';
        }
      }
    }
  }

  return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H +
         '" width="100%" height="' + H + '" preserveAspectRatio="xMinYMin meet" ' +
         'style="max-width:100%; display:block;">' + svg + '</svg>';
}

// -----------------------------------------------------------------------------
// Flank charts: two repeat-density panels (gar flank + mac flanks) plus a
// flanking-density bar table summarising mean/max for each annotated class.
// -----------------------------------------------------------------------------

function _csBuildFlankCharts(bp) {
  // Per-species panels: full-flank density curve from the loaded TEfull JSONs.
  // If the TEfull layer isn't loaded for the chrom, show an empty-state with
  // a "load it on the boundaries page first" hint.
  const garPanel = _csBuildOneFlankPanel(
    'gar', bp.gar_chr, (bp.gar_pos_start + bp.gar_pos_end) / 2,
    'all_TE',
  );
  // Mac panels: prev block end on prev mac_chr, next block start on next mac_chr.
  // When same mac_chr, draw one panel centred on the inversion centre.
  const macPrev = bp.prev_block || {};
  const macNext = bp.next_block || {};
  let macPanels = '';
  if (macPrev.mac_chr === macNext.mac_chr && macPrev.mac_chr) {
    // Inversion case — one mac chrom, anchor at the centre between block ends
    const center = (macPrev.mac_end_bp + macNext.mac_start_bp) / 2;
    macPanels += _csBuildOneFlankPanel(
      'mac', macPrev.mac_chr, center, 'all_TE',
      'mac flank ' + _esc(macPrev.mac_chr),
    );
  } else {
    if (macPrev.mac_chr) {
      macPanels += _csBuildOneFlankPanel(
        'mac', macPrev.mac_chr, macPrev.mac_end_bp, 'all_TE',
        'mac flank ' + _esc(macPrev.mac_chr) + ' (prev)',
      );
    }
    if (macNext.mac_chr) {
      macPanels += _csBuildOneFlankPanel(
        'mac', macNext.mac_chr, macNext.mac_start_bp, 'all_TE',
        'mac flank ' + _esc(macNext.mac_chr) + ' (next)',
      );
    }
  }

  // Flank bar table (the at-a-glance per-class summary)
  const flankTable = _csBuildFlankTable(bp);

  return garPanel + macPanels + flankTable;
}

function _csBuildOneFlankPanel(species, chrom, anchorBp, defaultClass, titleOverride) {
  const ui = _ensureCrossSpeciesState();
  const flank = ui.flank_bp || CROSS_SPECIES_FLANK_DEFAULT_BP;
  const speciesLabel = (species === 'gar')
    ? ((state.crossSpecies && (state.crossSpecies.species_query || {}).name) || 'gar')
    : ((state.crossSpecies && (state.crossSpecies.species_target || {}).name) || 'mac');
  const titleHtml =
    '<div class="cs-rd-panel-title">' +
      (titleOverride || ('flank \u00B7 <i>' + _esc(speciesLabel) + '</i> ' +
                          '<code>' + _esc(chrom) + '</code> ' +
                          'around <code>' + (anchorBp / 1e6).toFixed(3) + ' Mb</code>')) +
      ' <span class="cs-rd-source">\u00B7 ' + _esc(defaultClass) + ' \u00B1' +
      (flank / 1000) + ' kb</span>' +
    '</div>';
  // Source: TEfull JSON for chrom, if loaded
  const entry = (typeof _getRepeatDensity === 'function') ? _getRepeatDensity(chrom) : null;
  if (!entry || !entry.chrom_block || !entry.chrom_block.by_class || !entry.chrom_block.by_class[defaultClass]) {
    return '<div class="cs-rd-panel">' + titleHtml +
      '<div class="cs-rd-panel-empty">' +
        'No <code>' + _esc(defaultClass) + '</code> density loaded for <code>' + _esc(chrom) + '</code>. ' +
        'Drop the <code>' + _esc(chrom) + '_repeat_density_TEfull.json</code> on the boundaries page first.' +
      '</div></div>';
  }
  const block = entry.chrom_block;
  const dens = block.by_class[defaultClass].densities || [];
  const wcMb = block.window_centers_mb || [];
  // Filter to the flank window
  const lo = anchorBp - flank;
  const hi = anchorBp + flank;
  const xs = [];
  const ys = [];
  let chromMeanSum = 0, chromMeanCnt = 0;
  for (let i = 0; i < dens.length; i++) {
    const d = dens[i];
    if (Number.isFinite(d)) { chromMeanSum += d; chromMeanCnt += 1; }
    const cMb = wcMb[i];
    if (cMb == null) continue;
    const cBp = cMb * 1e6;
    if (cBp < lo || cBp > hi) continue;
    if (!Number.isFinite(d)) continue;
    xs.push(cBp);
    ys.push(d);
  }
  const chromMean = chromMeanCnt > 0 ? chromMeanSum / chromMeanCnt : null;
  if (xs.length === 0) {
    return '<div class="cs-rd-panel">' + titleHtml +
      '<div class="cs-rd-panel-empty">No windows in the \u00B1' + (flank / 1000) + ' kb flank around ' +
      (anchorBp / 1e6).toFixed(3) + ' Mb on <code>' + _esc(chrom) + '</code>.</div>' +
      '</div>';
  }
  // SVG plot — Spalax-figure idiom: scatter dots + simple smoothed line.
  const W = 760, H = 130;
  const margins = { top: 14, right: 18, bottom: 24, left: 42 };
  const innerW = W - margins.left - margins.right;
  const innerH = H - margins.top - margins.bottom;
  const xMin = lo, xMax = hi, xRange = xMax - xMin;
  const yMin = 0, yMax = 1;
  const toX = (bp_) => margins.left + ((bp_ - xMin) / xRange) * innerW;
  const toY = (d) => margins.top + (1 - (d - yMin) / (yMax - yMin)) * innerH;

  // Border + grid
  let svg =
    '<rect x="' + margins.left + '" y="' + margins.top + '" width="' + innerW + '" height="' + innerH +
    '" fill="none" stroke="var(--rule)" stroke-width="0.8" />';
  // Y-axis grid lines + labels
  for (const yv of [0, 0.25, 0.5, 0.75, 1.0]) {
    const y = toY(yv);
    svg +=
      '<line x1="' + margins.left + '" x2="' + (margins.left + innerW) +
      '" y1="' + y.toFixed(1) + '" y2="' + y.toFixed(1) +
      '" stroke="var(--rule)" stroke-width="0.5" stroke-opacity="0.5" />' +
      '<text x="' + (margins.left - 6) + '" y="' + (y + 3).toFixed(1) +
      '" font-size="9.5" fill="var(--ink-dim)" text-anchor="end" font-family="ui-monospace, monospace">' +
      yv.toFixed(2) + '</text>';
  }
  // X-axis ticks (5 evenly spaced)
  for (let i = 0; i <= 4; i++) {
    const xv = xMin + (i / 4) * xRange;
    const x = toX(xv);
    svg +=
      '<line x1="' + x.toFixed(1) + '" x2="' + x.toFixed(1) +
      '" y1="' + (margins.top + innerH) + '" y2="' + (margins.top + innerH + 3) +
      '" stroke="var(--ink-dim)" stroke-width="0.5" />' +
      '<text x="' + x.toFixed(1) + '" y="' + (margins.top + innerH + 14) +
      '" font-size="9.5" fill="var(--ink-dim)" text-anchor="middle" font-family="ui-monospace, monospace">' +
      (xv / 1e6).toFixed(2) + '</text>';
  }
  // Anchor vertical line
  const ax = toX(anchorBp);
  svg +=
    '<line x1="' + ax.toFixed(1) + '" x2="' + ax.toFixed(1) +
    '" y1="' + margins.top + '" y2="' + (margins.top + innerH) +
    '" stroke="var(--accent)" stroke-width="1.2" stroke-dasharray="3 2" />';
  // Chrom-mean horizontal reference line (helps the eye judge enrichment)
  if (chromMean != null) {
    const my = toY(chromMean);
    svg +=
      '<line x1="' + margins.left + '" x2="' + (margins.left + innerW) +
      '" y1="' + my.toFixed(1) + '" y2="' + my.toFixed(1) +
      '" stroke="var(--ink-dimmer)" stroke-width="0.8" stroke-dasharray="4 3" />' +
      '<text x="' + (margins.left + innerW - 4) + '" y="' + (my - 3).toFixed(1) +
      '" font-size="9" fill="var(--ink-dimmer)" text-anchor="end" font-family="ui-monospace, monospace">' +
      'chrom mean ' + chromMean.toFixed(2) + '</text>';
  }
  // Scatter
  for (let i = 0; i < xs.length; i++) {
    const x = toX(xs[i]);
    const y = toY(ys[i]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    const fillColor = (species === 'gar') ? '#4fa3ff' : '#3cc08a';
    svg += '<circle cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) +
           '" r="2" fill="' + fillColor + '" fill-opacity="0.5" />';
  }
  // Simple LOESS-ish smoothed line (moving average over 5 points, light-weight)
  if (xs.length >= 4) {
    const win = Math.max(3, Math.min(11, Math.round(xs.length / 6)));
    const half = Math.floor(win / 2);
    const path = [];
    for (let i = 0; i < xs.length; i++) {
      let s = 0, c = 0;
      for (let j = Math.max(0, i - half); j <= Math.min(xs.length - 1, i + half); j++) {
        s += ys[j]; c += 1;
      }
      const yv = s / c;
      const x = toX(xs[i]);
      const y = toY(yv);
      if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      path.push(x.toFixed(1) + ',' + y.toFixed(1));
    }
    if (path.length >= 2) {
      const lineColor = (species === 'gar') ? '#2c5d8c' : '#1f7a52';
      svg += '<polyline points="' + path.join(' ') + '" fill="none" stroke="' +
             lineColor + '" stroke-width="1.5" stroke-opacity="0.85" />';
    }
  }
  // Y-axis title
  svg += '<text x="13" y="' + (margins.top + innerH / 2) + '" font-size="10" fill="var(--ink-dim)" ' +
         'text-anchor="middle" font-family="ui-monospace, monospace" ' +
         'transform="rotate(-90, 13, ' + (margins.top + innerH / 2) + ')">' +
         _esc(defaultClass) + '</text>';
  // X-axis title
  svg += '<text x="' + (margins.left + innerW / 2) + '" y="' + (H - 4) +
         '" font-size="10" fill="var(--ink-dim)" text-anchor="middle" font-family="ui-monospace, monospace">' +
         'position on <code>' + _esc(chrom) + '</code> (Mb)</text>';
  return '<div class="cs-rd-panel">' + titleHtml +
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H +
    '" width="100%" height="' + H + '" preserveAspectRatio="none" ' +
    'style="max-width:100%;">' + svg + '</svg>' +
    '</div>';
}

function _csBuildFlankTable(bp) {
  // Build a per-class table: class | gar mean (bar) | mac prev mean | mac next mean
  // Sort: priority classes (all_TE first), then alphabetical.
  const garMap = bp.flanking_repeat_density_gar || {};
  const macWrap = bp.flanking_repeat_density_mac || {};
  const macPrev = (macWrap.prev || {}).by_class || {};
  const macNext = (macWrap.next || {}).by_class || {};
  const allClasses = new Set();
  for (const c of Object.keys(garMap))  allClasses.add(c);
  for (const c of Object.keys(macPrev)) allClasses.add(c);
  for (const c of Object.keys(macNext)) allClasses.add(c);
  if (allClasses.size === 0) {
    return '<div class="cs-rd-panel">' +
      '<div class="cs-rd-panel-title">flanking density by class</div>' +
      '<div class="cs-rd-panel-empty">No flanking-density data on this breakpoint. ' +
      'Re-run <code>STEP_CS01_extract_breakpoints.py</code> with <code>--te-dir</code> set.</div>' +
      '</div>';
  }
  const PRIO = (typeof PRIORITY_CLASSES !== 'undefined') ? PRIORITY_CLASSES : ['all_TE'];
  const TSDS = (typeof TSD_CLASSES !== 'undefined') ? TSD_CLASSES : [];
  const order = Array.from(allClasses).sort((a, b) => {
    const pa = PRIO.indexOf(a); const pb = PRIO.indexOf(b);
    if (pa !== -1 || pb !== -1) {
      if (pa === -1) return 1;
      if (pb === -1) return -1;
      return pa - pb;
    }
    const ta = TSDS.includes(a), tb = TSDS.includes(b);
    if (ta && !tb) return -1;
    if (!ta && tb) return 1;
    return a.localeCompare(b);
  });
  const fmt = (v) => (typeof v === 'number') ? v.toFixed(3) : '\u2014';
  const bar = (v, side) => {
    if (typeof v !== 'number') return '';
    const w = Math.min(100, Math.max(0, v * 100));
    return '<span class="cs-flank-bar"><span class="cs-flank-bar-fill cs-flank-bar-' + side +
           '" style="width:' + w.toFixed(1) + '%;"></span></span>';
  };
  const showMacPrev = (macWrap.prev || {}).mac_chr;
  const showMacNext = (macWrap.next || {}).mac_chr;
  const sameMac = showMacPrev && showMacNext && showMacPrev === showMacNext;
  let html = '<div class="cs-rd-panel">' +
    '<div class="cs-rd-panel-title">flanking density by class' +
      ' <span class="cs-rd-source">\u00B7 mean over windows in \u00B1' +
      (((bp.flanking_repeat_density_gar || {}).all_TE || {}).n_windows != null
        ? ('the recorded flank window')
        : 'flank') + '</span>' +
    '</div>' +
    '<table class="cs-flank-table">' +
    '<thead><tr>' +
      '<th>class</th>' +
      '<th colspan="2">gar flank mean</th>' +
      (sameMac ? '<th colspan="2">mac flank mean</th>'
               : '<th colspan="2">mac prev mean</th><th colspan="2">mac next mean</th>') +
    '</tr></thead><tbody>';
  for (const c of order) {
    const g = (garMap[c] || {}).mean;
    const mp = (macPrev[c] || {}).mean;
    const mn = (macNext[c] || {}).mean;
    html +=
      '<tr>' +
        '<td class="cs-flank-class"><code>' + _esc(c) + '</code></td>' +
        '<td class="cs-flank-num">' + fmt(g) + '</td>' +
        '<td>' + bar(g, 'gar') + '</td>';
    if (sameMac) {
      // Use prev mean as the "same mac" representative; show next in a secondary row would over-clutter
      const merged = (typeof mp === 'number' && typeof mn === 'number') ? (mp + mn) / 2
                   : (typeof mp === 'number' ? mp : (typeof mn === 'number' ? mn : null));
      html +=
        '<td class="cs-flank-num">' + fmt(merged) + '</td>' +
        '<td>' + bar(merged, 'mac') + '</td>';
    } else {
      html +=
        '<td class="cs-flank-num">' + fmt(mp) + '</td>' +
        '<td>' + bar(mp, 'mac') + '</td>' +
        '<td class="cs-flank-num">' + fmt(mn) + '</td>' +
        '<td>' + bar(mn, 'mac') + '</td>';
    }
    html += '</tr>';
  }
  html += '</tbody></table></div>';
  return html;
}

// -----------------------------------------------------------------------------
// Top-level page renderer + nav helpers
// -----------------------------------------------------------------------------

function _renderCrossSpeciesPage() {
  _renderCrossSpeciesToolbar();
  _renderCrossSpeciesCatalogue();
  _renderCrossSpeciesFocus();
  // v4 turn 75: macro-synteny section. Renders only when cs_breakpoints v2
  // (with synteny_blocks) is loaded. Fail-soft so a partial v1 atlas still
  // works for users whose cs_breakpoints JSON is from before the schema bump.
  if (typeof _renderCrossSpeciesSynteny === 'function') {
    try { _renderCrossSpeciesSynteny(); } catch (e) {
      console.warn('[crossSpecies] synteny render failed:', e && e.message ? e.message : e);
    }
  }
  // turn 115: dot plot panel. Renders when wfmash synteny_blocks (cs v2) or
  // mashmap multi-resolution (state.dotplotMashmap) is loaded. Fail-soft.
  if (typeof _renderCrossSpeciesDotplot === 'function') {
    try { _renderCrossSpeciesDotplot(); } catch (e) {
      console.warn('[crossSpecies] dotplot render failed:', e && e.message ? e.message : e);
    }
  }
  // Sticky-disappearing header on scroll. When user scrolls down inside
  // either the catalogue list or the focus panel, hide the page-16 header
  // so the breakpoint detail can use the full vertical space. Show it again
  // when scrolling up. Wired once per session via state flag; idempotent
  // across re-renders.
  if (!state._csHeaderScrollWired) {
    state._csHeaderScrollWired = true;
    const header = document.getElementById('page16Header');
    if (header) {
      const SCROLL_THRESHOLD = 24;            // px before we react (anti-jitter)
      const HEADER_HIDE_AFTER = 80;           // px scrolled before allowing hide
      const lastY = { focus: 0, cat: 0 };
      const apply = (delta, depth) => {
        if (depth < HEADER_HIDE_AFTER) {
          // Near the top → always show
          header.style.transition = 'opacity 0.18s ease, max-height 0.18s ease, padding 0.18s ease, margin 0.18s ease';
          header.style.opacity = '1';
          header.style.maxHeight = '';
          header.style.overflow = '';
          header.style.paddingTop = '';
          header.style.paddingBottom = '';
          header.style.marginBottom = '';
          return;
        }
        if (delta > SCROLL_THRESHOLD) {
          // Scrolling down past threshold → hide
          header.style.transition = 'opacity 0.18s ease, max-height 0.18s ease, padding 0.18s ease, margin 0.18s ease';
          header.style.opacity = '0';
          header.style.maxHeight = '0px';
          header.style.overflow = 'hidden';
          header.style.paddingTop = '0';
          header.style.paddingBottom = '0';
          header.style.marginBottom = '0';
        } else if (delta < -SCROLL_THRESHOLD) {
          // Scrolling up past threshold → show
          header.style.transition = 'opacity 0.18s ease, max-height 0.18s ease, padding 0.18s ease, margin 0.18s ease';
          header.style.opacity = '1';
          header.style.maxHeight = '';
          header.style.overflow = '';
          header.style.paddingTop = '';
          header.style.paddingBottom = '';
          header.style.marginBottom = '';
        }
      };
      const onFocus = (e) => {
        const y = e.target.scrollTop;
        apply(y - lastY.focus, y);
        lastY.focus = y;
      };
      const onCat = (e) => {
        const y = e.target.scrollTop;
        apply(y - lastY.cat, y);
        lastY.cat = y;
      };
      // Wire to both scrollable children. Re-renders replace innerHTML of
      // those panels' contents, but the wrapper div elements (#csFocus and
      // #csCatalogue) stay, so listeners attached to them survive.
      const focusEl = document.getElementById('csFocus');
      const catEl = document.getElementById('csCatalogue');
      if (focusEl) focusEl.addEventListener('scroll', onFocus, { passive: true });
      if (catEl) catEl.addEventListener('scroll', onCat, { passive: true });
    }
  }
}

// Step the active selection up/down within the currently-filtered list.
function _crossSpeciesStep(delta) {
  const ui = _ensureCrossSpeciesState();
  const bps = _csFilteredSorted();
  if (bps.length === 0) return;
  let idx = bps.findIndex(b => b.id === ui.active_id);
  if (idx < 0) idx = (delta > 0) ? -1 : bps.length;
  let nextIdx = idx + delta;
  if (nextIdx < 0) nextIdx = bps.length - 1;
  if (nextIdx >= bps.length) nextIdx = 0;
  ui.active_id = bps[nextIdx].id;
  _renderCrossSpeciesCatalogue();
  _renderCrossSpeciesFocus();
  // Scroll the active row into view
  const slot = document.getElementById('csCatalogue');
  if (slot) {
    const el = slot.querySelector('.cs-row-active');
    if (el && typeof el.scrollIntoView === 'function') {
      el.scrollIntoView({ block: 'nearest' });
    }
  }
}

function _crossSpeciesToggleEventFilter(et) {
  const ui = _ensureCrossSpeciesState();
  if (!ui.filter.events) ui.filter.events = new Set();
  if (ui.filter.events.has(et)) ui.filter.events.delete(et);
  else ui.filter.events.add(et);
  if (ui.filter.events.size === 0) ui.filter.events = null;
  _renderCrossSpeciesPage();
}

function _crossSpeciesClearFilters() {
  const ui = _ensureCrossSpeciesState();
  ui.filter.events = null;
  ui.filter.search = '';
  _renderCrossSpeciesPage();
}

// Keyboard wiring — only acts when page16 is the active page.
function _wireCrossSpeciesKeys() {
  if (state._crossSpeciesKeysBound) return;
  document.addEventListener('keydown', (e) => {
    const page = document.getElementById('page16');
    if (!page || !page.classList.contains('active')) return;
    const tag = (e.target && e.target.tagName) || '';
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
    if (e.target && e.target.isContentEditable) return;
    if (e.altKey || e.ctrlKey || e.metaKey) return;
    if (e.key === 'ArrowDown' && !e.shiftKey) { e.preventDefault(); _crossSpeciesStep(+1); return; }
    if (e.key === 'ArrowUp'   && !e.shiftKey) { e.preventDefault(); _crossSpeciesStep(-1); return; }
    if (e.key === 'i' || e.key === 'I') { e.preventDefault(); _crossSpeciesToggleEventFilter('inversion'); return; }
    if (e.key === 't' || e.key === 'T') { e.preventDefault(); _crossSpeciesToggleEventFilter('translocation'); return; }
    if (e.key === 'f' || e.key === 'F') { e.preventDefault(); _crossSpeciesToggleEventFilter('fission_or_fusion'); return; }
    if (e.key === 'Escape') { e.preventDefault(); _crossSpeciesClearFilters(); return; }
  });
  state._crossSpeciesKeysBound = true;
}

// =============================================================================
// v4 turn 75 — Cross-species synteny module (page16, sub-section)
// -----------------------------------------------------------------------------
// Macro-synteny interpretation layer on top of the cs_breakpoints_v1 v2 data.
// Renders inside the focus column of page16, BELOW the per-breakpoint focus
// view. Components:
//   1. Sankey diagram: Cgar LGs (left) ↔ classification tier (middle) ↔
//      Cmac LGs (right). Flow widths = aligned bp. Three classification
//      classes: one_to_one, fusion_in_target, fission_in_target, many_to_many.
//   2. Per-Cgar-chromosome relationship table.
//   3. Inversion-candidate overlay: each promoted candidate from
//      state.candidateList is annotated against the synteny blocks on its
//      Gar chrom — distance to nearest synteny edge, inside-vs-spanning,
//      proximity to fusion/fission boundary, evolutionary interpretation.
//   4. Permutation test: 10,000 random-position repositionings of each
//      candidate on its chrom; reports observed-vs-expected proximity to
//      nearest synteny edge.
//
// All synteny analysis is derived on-the-fly from
// state.crossSpecies.synteny_blocks. No new file format. Cached via
// state._csSyntenyCache (invalidated on cs_breakpoints reload or chrom
// list change).
// =============================================================================

const CS_SYNTENY_NOISE_FLOOR_FRAC = 0.05;   // <5% of chrom aligned bp = ignore
const CS_SYNTENY_PERM_N           = 10000;
const CS_SYNTENY_PROXIMITY_QUANTILE = 0.25; // "near edge" = bottom 25% of distances
// Classification colors. Match the chip palette from the breakpoints layer.
const CS_REL_DEF = {
  one_to_one:        { label: '1:1 conserved',       color: '#3cc08a' },
  fusion_in_target:  { label: 'fusion in target',    color: '#f0a35e' },
  fission_in_target: { label: 'fission in target',   color: '#7ad3db' },
  many_to_many:      { label: 'many-to-many',        color: '#b07cf7' },
};

// -----------------------------------------------------------------------------
// Synteny analysis (atlas-side)
// -----------------------------------------------------------------------------

// Returns the synteny blocks if the v2 data layer is loaded, else null.
function _csGetSyntenyBlocks() {
  const cs = state.crossSpecies;
  if (!cs || !Array.isArray(cs.synteny_blocks)) return null;
  return cs.synteny_blocks;
}

// Compute per-chromosome aligned-bp coverage maps and classify the macro
// synteny relationship per Cgar chrom and per Cmac chrom.
//
// Returns:
//   {
//     gar_partners: { <gar_chr>: [ {mac_chr, aligned_bp, frac, dominant: bool}, ... ] }
//     mac_partners: { <mac_chr>: [ {gar_chr, aligned_bp, frac, dominant: bool}, ... ] }
//     gar_class:    { <gar_chr>: <one_of CS_REL_DEF keys> }
//     mac_class:    { <mac_chr>: <one_of CS_REL_DEF keys> }
//     flows:        [ {gar_chr, mac_chr, aligned_bp, gar_frac, mac_frac, class}, ... ]
//   }
//
// Logic:
//   - For each gar chrom, count aligned_bp per mac chrom partner.
//   - Drop partners with frac < CS_SYNTENY_NOISE_FLOOR_FRAC (5%) — those are
//     short syntenic islands, not chromosome-level relationships.
//   - The "dominant" partners that survive determine the gar chrom's
//     relationship class: 1 partner → check the partner's gar list; if
//     that mac chrom has multiple dominant gar partners then the gar chrom
//     participates in a "fusion in target" (multiple gar → 1 mac).
//   - Multiple gar partners → "fission in target" (1 gar → multiple mac).
//   - Combined → "many_to_many".
//
// Caching: results stored in state._csSyntenyCache, invalidated on data
// reload (clear in _storeCrossSpecies / _clearCrossSpecies).
function _csComputeSynteny() {
  if (state._csSyntenyCache) return state._csSyntenyCache;
  const blocks = _csGetSyntenyBlocks();
  if (!blocks) return null;
  // Per-(gar,mac) aligned bp tally
  const tally = new Map();   // key="gar||mac" -> {gar_chr, mac_chr, aligned_bp}
  const garTotal = new Map();
  const macTotal = new Map();
  for (const b of blocks) {
    const key = b.gar_chr + '||' + b.mac_chr;
    const cur = tally.get(key);
    const sz = b.block_size_bp || (b.gar_end - b.gar_start) || 0;
    if (cur) cur.aligned_bp += sz;
    else tally.set(key, { gar_chr: b.gar_chr, mac_chr: b.mac_chr, aligned_bp: sz });
    garTotal.set(b.gar_chr, (garTotal.get(b.gar_chr) || 0) + sz);
    macTotal.set(b.mac_chr, (macTotal.get(b.mac_chr) || 0) + sz);
  }
  // Build partner lists
  const gar_partners = {};
  const mac_partners = {};
  for (const t of tally.values()) {
    const garFrac = t.aligned_bp / (garTotal.get(t.gar_chr) || 1);
    const macFrac = t.aligned_bp / (macTotal.get(t.mac_chr) || 1);
    if (!gar_partners[t.gar_chr]) gar_partners[t.gar_chr] = [];
    if (!mac_partners[t.mac_chr]) mac_partners[t.mac_chr] = [];
    gar_partners[t.gar_chr].push({
      mac_chr: t.mac_chr, aligned_bp: t.aligned_bp,
      frac: garFrac, dominant: garFrac >= CS_SYNTENY_NOISE_FLOOR_FRAC,
    });
    mac_partners[t.mac_chr].push({
      gar_chr: t.gar_chr, aligned_bp: t.aligned_bp,
      frac: macFrac, dominant: macFrac >= CS_SYNTENY_NOISE_FLOOR_FRAC,
    });
  }
  // Sort partner lists by aligned_bp desc
  for (const c in gar_partners) gar_partners[c].sort((a, b) => b.aligned_bp - a.aligned_bp);
  for (const c in mac_partners) mac_partners[c].sort((a, b) => b.aligned_bp - a.aligned_bp);

  // Classify per gar chrom
  const gar_class = {};
  for (const c in gar_partners) {
    const dom = gar_partners[c].filter(p => p.dominant);
    if (dom.length === 0) {
      gar_class[c] = 'many_to_many';   // shouldn't happen but safe default
      continue;
    }
    if (dom.length === 1) {
      // Only one mac partner. Does that mac have multiple gar partners?
      const macC = dom[0].mac_chr;
      const macDoms = (mac_partners[macC] || []).filter(p => p.dominant);
      if (macDoms.length === 1) gar_class[c] = 'one_to_one';
      else                       gar_class[c] = 'fusion_in_target';
    } else {
      // Multiple mac partners — fission in target (gar → multiple mac)
      // Distinguish from many_to_many: do ANY of those mac partners ALSO have
      // multiple gar partners? If yes, it's many_to_many.
      const anyMacWithMultGar = dom.some(p =>
        ((mac_partners[p.mac_chr] || []).filter(x => x.dominant).length > 1));
      gar_class[c] = anyMacWithMultGar ? 'many_to_many' : 'fission_in_target';
    }
  }
  // Symmetric classification from the Mac perspective. We flip the labels:
  //   - mac with 1 gar partner whose gar has multiple mac partners → mac is
  //     one of the pieces of a "gar→multiple_mac" pattern, which from mac's
  //     POV means it shares its gar source with siblings → fusion_in_target
  //     when target=gar (i.e. mac merges into gar) — but our convention
  //     keeps target=mac, so this is "fission_in_target" being seen from the
  //     other side. To avoid mislabelling, mac_class stores the SAME class
  //     LABEL as the corresponding gar chroms — easier to interpret.
  const mac_class = {};
  for (const c in mac_partners) {
    const dom = mac_partners[c].filter(p => p.dominant);
    if (dom.length === 0) { mac_class[c] = 'many_to_many'; continue; }
    if (dom.length === 1) {
      const garC = dom[0].gar_chr;
      const garDoms = (gar_partners[garC] || []).filter(p => p.dominant);
      // mac has 1 gar partner. If gar has multiple mac partners, this mac
      // is one piece of a fission (gar→multi_mac) — label matches gar_class
      // on garC, which is 'fission_in_target'.
      if (garDoms.length === 1) mac_class[c] = 'one_to_one';
      else                       mac_class[c] = 'fission_in_target';
    } else {
      // Multiple gar partners → fusion (multi_gar → 1_mac)
      const anyGarWithMultMac = dom.some(p =>
        ((gar_partners[p.gar_chr] || []).filter(x => x.dominant).length > 1));
      mac_class[c] = anyGarWithMultMac ? 'many_to_many' : 'fusion_in_target';
    }
  }

  // Flat flow list (for Sankey drawing)
  const flows = [];
  for (const t of tally.values()) {
    const garFrac = t.aligned_bp / (garTotal.get(t.gar_chr) || 1);
    const macFrac = t.aligned_bp / (macTotal.get(t.mac_chr) || 1);
    flows.push({
      gar_chr: t.gar_chr,
      mac_chr: t.mac_chr,
      aligned_bp: t.aligned_bp,
      gar_frac: garFrac,
      mac_frac: macFrac,
      class: gar_class[t.gar_chr] || 'many_to_many',
    });
  }
  flows.sort((a, b) => b.aligned_bp - a.aligned_bp);

  const result = { gar_partners, mac_partners, gar_class, mac_class, flows };
  state._csSyntenyCache = result;
  return result;
}

// Compute per-gar-chrom synteny block boundaries — the "edges" inversion
// breakpoints will be measured against. An edge is the START of a block
// or the END of a block ON A GAR CHROM. Returns:
//   { <gar_chr>: [<edge_bp>, ...] sorted ascending }
// Edges include the chromosome ends (0 and chrom_length).
function _csSyntenyEdgesByChrom() {
  if (state._csSyntenyEdgesCache) return state._csSyntenyEdgesCache;
  const blocks = _csGetSyntenyBlocks();
  if (!blocks) return null;
  const edges = {};
  for (const b of blocks) {
    if (!edges[b.gar_chr]) edges[b.gar_chr] = new Set();
    edges[b.gar_chr].add(b.gar_start);
    edges[b.gar_chr].add(b.gar_end);
  }
  // Add chrom ends
  const lens = (state.crossSpecies && state.crossSpecies.chrom_lengths_query) || {};
  for (const c in edges) {
    edges[c].add(0);
    if (lens[c]) edges[c].add(lens[c]);
  }
  const out = {};
  for (const c in edges) out[c] = Array.from(edges[c]).sort((a, b) => a - b);
  state._csSyntenyEdgesCache = out;
  return out;
}

// Distance from `pos` (bp) to the nearest synteny edge on `chrom`.
// Returns null if chrom has no edges loaded.
function _csDistanceToNearestEdge(chrom, pos) {
  const edges = _csSyntenyEdgesByChrom();
  if (!edges || !edges[chrom] || edges[chrom].length === 0) return null;
  const arr = edges[chrom];
  // Binary search
  let lo = 0, hi = arr.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] < pos) lo = mid + 1;
    else hi = mid;
  }
  // Now arr[lo] is the smallest edge >= pos. Compare against arr[lo-1] too.
  let best = Math.abs(arr[lo] - pos);
  if (lo > 0) best = Math.min(best, Math.abs(arr[lo - 1] - pos));
  return best;
}

// For each promoted inversion candidate, compute its synteny context.
// `state.candidateList` is the canonical source — same as on page2/page3.
// Returns a list of:
//   { candidate_id, chrom, start_bp, end_bp, len_bp, center_bp,
//     dist_left_edge_bp, dist_right_edge_bp, dist_min_edge_bp,
//     inside_single_block: bool, spanning_n_blocks: int,
//     chrom_class, partner_mac_chr, interpretation: <string> }
//
// `dist_min_edge_bp` is min(dist from start to nearest edge, dist from end
// to nearest edge) — this is the "is at least one breakpoint near a synteny
// edge?" metric the manuscript narrative cares about. The permutation test
// uses this as the test statistic.
function _csInversionContexts() {
  if (state._csInversionContextCache) return state._csInversionContextCache;
  const synteny = _csComputeSynteny();
  if (!synteny) return null;
  const blocks = _csGetSyntenyBlocks();
  if (!blocks) return null;
  // Group blocks by gar_chr for efficient inside-block lookup
  const byGar = {};
  for (const b of blocks) {
    if (!byGar[b.gar_chr]) byGar[b.gar_chr] = [];
    byGar[b.gar_chr].push(b);
  }
  // Promoted candidates only — match what the rest of the atlas considers
  // canonical inversion intervals.
  const cands = (state.candidateList || []).filter(c =>
    c && c.chrom && (Number.isFinite(c.start_bp) || Number.isFinite(c.startBp)));
  const out = [];
  for (const c of cands) {
    const chrom = c.chrom;
    const start = c.start_bp != null ? c.start_bp : c.startBp;
    const end   = c.end_bp   != null ? c.end_bp   : c.endBp;
    if (!Number.isFinite(start) || !Number.isFinite(end)) continue;
    const center = (start + end) / 2;
    const distLeft  = _csDistanceToNearestEdge(chrom, start);
    const distRight = _csDistanceToNearestEdge(chrom, end);
    const distMin = (distLeft != null && distRight != null)
      ? Math.min(distLeft, distRight)
      : (distLeft != null ? distLeft : distRight);
    // Count synteny blocks that overlap [start, end] on this chrom
    const chromBlocks = byGar[chrom] || [];
    let overlapping = 0;
    for (const b of chromBlocks) {
      if (b.gar_end >= start && b.gar_start <= end) overlapping++;
    }
    const insideSingle = (overlapping === 1);
    const chromClass = synteny.gar_class[chrom] || 'unknown';
    // Pick the dominant mac partner for this gar chrom for the table
    const partners = (synteny.gar_partners[chrom] || []).filter(p => p.dominant);
    const partnerMac = partners.length === 1
      ? partners[0].mac_chr
      : (partners.length > 1
         ? partners.map(p => p.mac_chr).join('+')
         : '?');
    const interpretation = _csInversionInterpretation({
      chromClass, insideSingle, distMin,
      spanning_n_blocks: overlapping,
    });
    out.push({
      candidate_id: c.id || c.candidate_id || ('CAND_' + (out.length + 1).toString().padStart(3, '0')),
      chrom, start_bp: start, end_bp: end, len_bp: end - start, center_bp: center,
      dist_left_edge_bp:  distLeft,
      dist_right_edge_bp: distRight,
      dist_min_edge_bp:   distMin,
      inside_single_block: insideSingle,
      spanning_n_blocks: overlapping,
      chrom_class: chromClass,
      partner_mac_chr: partnerMac,
      interpretation,
    });
  }
  state._csInversionContextCache = out;
  return out;
}

function _csInversionInterpretation(ctx) {
  // Heuristic phrase. The permutation test (separate function) returns a
  // refined p-value-aware interpretation that overrides this default once
  // it's been run.
  const cls = ctx.chromClass;
  if (ctx.insideSingle) {
    if (cls === 'one_to_one') {
      return 'inside conserved 1:1 syntenic block — likely recent within-lineage polymorphism';
    }
    if (cls === 'fusion_in_target' || cls === 'fission_in_target') {
      return 'inside synteny block on a chromosome that participated in fusion/fission — recent inversion in a previously-rearranged chrom';
    }
    return 'inside a single synteny block on a complex chromosome';
  }
  // Spanning multiple blocks → at least one breakpoint is at a synteny edge
  if (ctx.spanning_n_blocks >= 2) {
    if (cls === 'fusion_in_target' || cls === 'fission_in_target') {
      return 'breakpoints near synteny edges on a fusion/fission chrom — possible reuse of ancient rearrangement boundary';
    }
    if (cls === 'many_to_many') {
      return 'breakpoints near synteny edges on a complex chrom — proceed with caution, may be assembly artifact';
    }
    return 'spans \u22652 synteny blocks — at least one breakpoint coincides with a synteny edge';
  }
  return 'context unclear';
}

// Permutation test: are inversion candidates' breakpoints closer to synteny
// edges than expected for random matched intervals on the same chromosome?
//
// Test statistic per candidate: dist_min_edge_bp = min(dist from start to
// nearest edge, dist from end to nearest edge). The null distribution is
// constructed by repeatedly drawing a random START position uniformly on
// the chrom (preserving the inversion length) and computing the same
// statistic for the random interval. This keeps inversion length matched —
// long inversions naturally have a higher chance of one endpoint being
// near an edge by virtue of spanning more chromosome.
//
// Args: opts = { n_perm = 10000, quantile = 0.25 }
// Returns: { per_candidate: [...], global: {...} }
function _csPermutationTest(opts) {
  opts = opts || {};
  const n_perm = opts.n_perm || CS_SYNTENY_PERM_N;
  const q = opts.quantile != null ? opts.quantile : CS_SYNTENY_PROXIMITY_QUANTILE;
  const ctxs = _csInversionContexts();
  if (!ctxs) return null;
  const lens = (state.crossSpecies && state.crossSpecies.chrom_lengths_query) || {};

  // Per-candidate: draw n_perm random start positions, compute matching
  // dist_min_edge for each, sort. Length-matched null preserves the
  // "longer inversions have shorter expected min-distance" effect that
  // a position-only null would miss.
  const per_candidate = [];
  let nearEdgeObs = 0, nearEdgeExp = 0, nValid = 0;
  for (const ctx of ctxs) {
    const len = lens[ctx.chrom];
    if (!len || len <= 0 || ctx.dist_min_edge_bp == null) {
      per_candidate.push({ ...ctx, percentile: null, p_value: null,
                           expected_dist_median: null,
                           is_near_edge_observed: null });
      continue;
    }
    const invLen = ctx.len_bp;
    if (invLen <= 0 || invLen >= len) {
      per_candidate.push({ ...ctx, percentile: null, p_value: null,
                           expected_dist_median: null,
                           is_near_edge_observed: null });
      continue;
    }
    const seed = _csSeedFromString(ctx.chrom + '|' + ctx.candidate_id);
    const rand = _csMulberry32(seed);
    const dists = new Float64Array(n_perm);
    for (let i = 0; i < n_perm; i++) {
      const startRand = rand() * (len - invLen);
      const endRand   = startRand + invLen;
      const dl = _csDistanceToNearestEdge(ctx.chrom, startRand);
      const dr = _csDistanceToNearestEdge(ctx.chrom, endRand);
      const dm = (dl != null && dr != null) ? Math.min(dl, dr)
              : (dl != null ? dl : (dr != null ? dr : NaN));
      dists[i] = dm;
    }
    // Sort the random distribution
    const sorted = Array.from(dists).filter(v => Number.isFinite(v)).sort((a, b) => a - b);
    if (sorted.length === 0) {
      per_candidate.push({ ...ctx, percentile: null, p_value: null,
                           expected_dist_median: null,
                           is_near_edge_observed: null });
      continue;
    }
    const obs = ctx.dist_min_edge_bp;
    // Percentile of observed distance in the random distribution (lower
    // percentile = closer to an edge than random)
    let lo = 0, hi = sorted.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (sorted[mid] < obs) lo = mid + 1; else hi = mid;
    }
    const pct = lo / sorted.length;
    const p_value = (lo + 1) / (sorted.length + 1);   // smoothed
    const medianRandom = sorted[Math.floor(sorted.length / 2)];
    // "Near edge" threshold: the lower q-quantile of the random distribution
    const threshIdx = Math.max(0, Math.floor(sorted.length * q) - 1);
    const isNearEdge = obs <= sorted[threshIdx];
    per_candidate.push({
      ...ctx,
      percentile: pct, p_value,
      expected_dist_median: medianRandom,
      is_near_edge_observed: isNearEdge,
    });
    nValid += 1;
    if (isNearEdge) nearEdgeObs += 1;
    nearEdgeExp += q;   // expected fraction "near edge" under null = q
  }
  const global = {
    n_candidates: ctxs.length,
    n_valid: nValid,
    n_perm,
    quantile: q,
    near_edge_observed: nearEdgeObs,
    near_edge_expected: nearEdgeExp,
    enrichment_ratio: nearEdgeExp > 0 ? (nearEdgeObs / nearEdgeExp) : null,
    binomial_p: _csBinomialTailP(nValid, q, nearEdgeObs),
  };
  return { per_candidate, global };
}

// Mulberry32 PRNG (deterministic, fast). Accepts a 32-bit seed.
function _csMulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function _csSeedFromString(s) {
  // Simple FNV-1a-ish hash → 32-bit seed
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return h >>> 0;
}

// One-sided upper-tail binomial p-value: P(X >= k | n, p).
// Computed via the regularized incomplete beta or, for our small sizes,
// by direct summation up to n. We use direct summation since n_candidates
// is typically <100.
function _csBinomialTailP(n, p, k) {
  if (n <= 0) return null;
  if (k <= 0) return 1.0;
  if (k > n) return 0.0;
  // log-space pmf for stability
  function logFact(n) {
    if (n < 2) return 0;
    // Stirling for large n; exact for small
    if (n < 200) {
      let s = 0;
      for (let i = 2; i <= n; i++) s += Math.log(i);
      return s;
    }
    // Stirling
    return n * Math.log(n) - n + 0.5 * Math.log(2 * Math.PI * n);
  }
  let logSum = -Infinity;
  for (let x = k; x <= n; x++) {
    const logPmf = logFact(n) - logFact(x) - logFact(n - x)
                 + x * Math.log(p) + (n - x) * Math.log(1 - p);
    // log-sum-exp
    if (logSum === -Infinity) logSum = logPmf;
    else logSum = Math.max(logSum, logPmf) + Math.log1p(Math.exp(-Math.abs(logSum - logPmf)));
  }
  return Math.exp(logSum);
}

// -----------------------------------------------------------------------------
// Sankey rendering
// -----------------------------------------------------------------------------
// 3-tier Sankey: gar (left) → classification (middle) → mac (right)
// Width/height proportional to aligned_bp. Flows colored by classification.
//
// Layout strategy:
//   - Tier 1 (gar): vertical column on the LEFT. Stacked rectangles, height
//     proportional to chrom aligned bp (sum of partner aligned_bp), in
//     natural-sort chrom order.
//   - Tier 2 (class): a thin band in the MIDDLE, four sub-bands stacked
//     by class. Each flow passes THROUGH its class sub-band.
//   - Tier 3 (mac): same as gar, on the RIGHT.
//   - Flows: smooth Bezier ribbons from gar to (class anchor) to mac, each
//     ribbon's vertical extent on the gar side == that flow's aligned bp.
//
// We don't try to minimise crossings — natural sort ordering is fine for
// visual interpretation. This is for chromosome-scale relationships, not
// gene-level synteny where crossing minimisation matters.

function _csNaturalSortChroms(arr) {
  // Sort LG01, LG02, LG10, LG28 etc. — natural numeric sort on the suffix
  return arr.slice().sort((a, b) => {
    const na = String(a).match(/(\d+)/);
    const nb = String(b).match(/(\d+)/);
    if (na && nb) {
      const ia = parseInt(na[1], 10), ib = parseInt(nb[1], 10);
      if (ia !== ib) return ia - ib;
    }
    return String(a).localeCompare(String(b));
  });
}

function _csBuildSankeySvg() {
  const synteny = _csComputeSynteny();
  if (!synteny) return '';
  const flows = synteny.flows;
  if (flows.length === 0) return '';

  const W = 760, H = 460;
  const padL = 70, padR = 70, padT = 26, padB = 24;
  const innerW = W - padL - padR;
  const innerH = H - padT - padB;

  // Aggregate per-gar / per-mac totals from the flows (so the Sankey doesn't
  // double-count any block that the noise floor dropped at the edges)
  const garTotals = new Map(), macTotals = new Map();
  for (const f of flows) {
    garTotals.set(f.gar_chr, (garTotals.get(f.gar_chr) || 0) + f.aligned_bp);
    macTotals.set(f.mac_chr, (macTotals.get(f.mac_chr) || 0) + f.aligned_bp);
  }
  const totalAligned = Array.from(garTotals.values()).reduce((s, v) => s + v, 0);
  if (totalAligned <= 0) return '';
  // Layout: stack chroms vertically with a small gap between each
  const garChroms = _csNaturalSortChroms(Array.from(garTotals.keys()));
  const macChroms = _csNaturalSortChroms(Array.from(macTotals.keys()));
  const chromGap = 3;
  // Compute available height after subtracting gaps
  const garUsable = innerH - chromGap * (garChroms.length - 1);
  const macUsable = innerH - chromGap * (macChroms.length - 1);
  // Per-chrom y-positions (top edge) and heights
  const garLayout = {};
  let y = padT;
  for (const c of garChroms) {
    const h = (garTotals.get(c) / totalAligned) * garUsable;
    garLayout[c] = { y, h };
    y += h + chromGap;
  }
  const macLayout = {};
  y = padT;
  for (const c of macChroms) {
    const h = (macTotals.get(c) / totalAligned) * macUsable;
    macLayout[c] = { y, h };
    y += h + chromGap;
  }
  // Class tier in the middle: 4 horizontal bands, height = sum of class flows
  const classBands = ['one_to_one', 'fusion_in_target', 'fission_in_target', 'many_to_many'];
  const classTotal = {};
  for (const k of classBands) classTotal[k] = 0;
  for (const f of flows) classTotal[f.class] = (classTotal[f.class] || 0) + f.aligned_bp;
  const classLayout = {};
  const midX = padL + innerW / 2;
  const midBandW = 80;
  const midBandX = midX - midBandW / 2;
  let cy = padT;
  for (const k of classBands) {
    const h = (classTotal[k] / totalAligned) * innerH;
    classLayout[k] = { y: cy, h, label: CS_REL_DEF[k].label, color: CS_REL_DEF[k].color };
    cy += h;
  }

  // SVG header
  let svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' + W + ' ' + H +
            '" width="100%" height="' + H + '" preserveAspectRatio="xMidYMid meet" style="max-width:100%;">';

  // Tier headers
  svg += '<text x="' + (padL + 6) + '" y="' + (padT - 10) +
         '" font-size="11" fill="var(--ink)" font-family="ui-monospace, monospace" font-weight="600">Cgar LGs</text>';
  svg += '<text x="' + midX + '" y="' + (padT - 10) +
         '" font-size="11" fill="var(--ink)" font-family="ui-monospace, monospace" font-weight="600" text-anchor="middle">classification</text>';
  svg += '<text x="' + (W - padR - 6) + '" y="' + (padT - 10) +
         '" font-size="11" fill="var(--ink)" font-family="ui-monospace, monospace" font-weight="600" text-anchor="end">Cmac LGs</text>';

  // Flow ribbons. For each flow, we need:
  //   - the gar slice (a sub-rectangle of gar_chrom's column at the left)
  //   - the class slice (a sub-rectangle of the class band in the middle)
  //   - the mac slice (a sub-rectangle of mac_chrom's column at the right)
  // We track a running offset within each of these so multiple flows from
  // the same gar chrom stack vertically (and similarly for mac and class).
  // NOTE: the class band needs flows from all gar chroms to stack within it
  // — easiest is to track running offset PER (gar, class) and PER (class, mac).
  // For a clean visualization we process flows in (class, gar_chrom) order.
  flows.sort((a, b) => {
    const ka = classBands.indexOf(a.class);
    const kb = classBands.indexOf(b.class);
    if (ka !== kb) return ka - kb;
    const ga = garChroms.indexOf(a.gar_chr);
    const gb = garChroms.indexOf(b.gar_chr);
    if (ga !== gb) return ga - gb;
    return a.mac_chr.localeCompare(b.mac_chr);
  });
  // Running offsets
  const garOffset = {}; for (const c of garChroms) garOffset[c] = 0;
  const macOffset = {}; for (const c of macChroms) macOffset[c] = 0;
  const classOffset = {}; for (const k of classBands) classOffset[k] = 0;

  // Draw flows BEFORE drawing the chrom rects so the rects sit on top
  for (const f of flows) {
    const garRect = garLayout[f.gar_chr];
    const macRect = macLayout[f.mac_chr];
    const clsRect = classLayout[f.class];
    if (!garRect || !macRect || !clsRect) continue;
    const flowFrac = f.aligned_bp / totalAligned;
    const flowH = flowFrac * innerH;   // approximate ribbon height
    const garY = garRect.y + garOffset[f.gar_chr];
    const macY = macRect.y + macOffset[f.mac_chr];
    const clsY = clsRect.y + classOffset[f.class];
    // Width of each slice — use the ribbon's flow ratio scaled to chrom height
    const gh = flowFrac / (garTotals.get(f.gar_chr) / totalAligned) * garRect.h;
    const mh = flowFrac / (macTotals.get(f.mac_chr) / totalAligned) * macRect.h;
    const ch = flowFrac / (classTotal[f.class] / totalAligned) * clsRect.h;
    // Anchor x's
    const x0 = padL + 16;            // right edge of gar column
    const x1 = midBandX;             // left edge of class band
    const x2 = midBandX + midBandW;  // right edge of class band
    const x3 = W - padR - 16;        // left edge of mac column
    const c0x = (x0 + x1) / 2;
    const c1x = (x2 + x3) / 2;
    // Bezier ribbons
    const color = clsRect.color;
    // Path: gar (top) → class (top) at x1, then class band straight, then class → mac
    const path =
      'M ' + x0 + ' ' + garY +
      ' C ' + c0x + ' ' + garY + ', ' + c0x + ' ' + clsY + ', ' + x1 + ' ' + clsY +
      ' L ' + x2 + ' ' + clsY +
      ' C ' + c1x + ' ' + clsY + ', ' + c1x + ' ' + macY + ', ' + x3 + ' ' + macY +
      ' L ' + x3 + ' ' + (macY + mh) +
      ' C ' + c1x + ' ' + (macY + mh) + ', ' + c1x + ' ' + (clsY + ch) + ', ' + x2 + ' ' + (clsY + ch) +
      ' L ' + x1 + ' ' + (clsY + ch) +
      ' C ' + c0x + ' ' + (clsY + ch) + ', ' + c0x + ' ' + (garY + gh) + ', ' + x0 + ' ' + (garY + gh) +
      ' Z';
    svg += '<path d="' + path + '" fill="' + color + '" fill-opacity="0.32" stroke="' + color +
           '" stroke-width="0.4" stroke-opacity="0.6">' +
           '<title>' + _esc(f.gar_chr) + ' \u2192 ' + _esc(f.mac_chr) +
           '\n' + (f.aligned_bp / 1e6).toFixed(2) + ' Mb aligned' +
           '\n' + (f.gar_frac * 100).toFixed(1) + '% of ' + _esc(f.gar_chr) +
           ', ' + (f.mac_frac * 100).toFixed(1) + '% of ' + _esc(f.mac_chr) +
           '\nclass: ' + _esc(CS_REL_DEF[f.class].label) +
           '</title></path>';
    garOffset[f.gar_chr]   += gh;
    macOffset[f.mac_chr]   += mh;
    classOffset[f.class]   += ch;
  }

  // Tier 1: gar chrom rects
  for (const c of garChroms) {
    const r = garLayout[c];
    svg += '<rect x="' + padL + '" y="' + r.y + '" width="16" height="' + r.h +
           '" fill="#4fa3ff" fill-opacity="0.85" stroke="var(--rule)" stroke-width="0.5">' +
           '<title>Cgar ' + _esc(c) + '\n' + (garTotals.get(c) / 1e6).toFixed(2) +
           ' Mb aligned</title></rect>';
    // Label to the LEFT of the rect
    svg += '<text x="' + (padL - 4) + '" y="' + (r.y + r.h / 2 + 3) + '" font-size="10" ' +
           'fill="var(--ink)" text-anchor="end" font-family="ui-monospace, monospace">' +
           _esc(c) + '</text>';
  }
  // Tier 3: mac chrom rects
  for (const c of macChroms) {
    const r = macLayout[c];
    svg += '<rect x="' + (W - padR - 16) + '" y="' + r.y + '" width="16" height="' + r.h +
           '" fill="#3cc08a" fill-opacity="0.85" stroke="var(--rule)" stroke-width="0.5">' +
           '<title>Cmac ' + _esc(c) + '\n' + (macTotals.get(c) / 1e6).toFixed(2) +
           ' Mb aligned</title></rect>';
    svg += '<text x="' + (W - padR + 4) + '" y="' + (r.y + r.h / 2 + 3) + '" font-size="10" ' +
           'fill="var(--ink)" text-anchor="start" font-family="ui-monospace, monospace">' +
           _esc(c) + '</text>';
  }
  // Tier 2: class band rects + labels
  for (const k of classBands) {
    const r = classLayout[k];
    if (r.h <= 0) continue;
    svg += '<rect x="' + midBandX + '" y="' + r.y + '" width="' + midBandW + '" height="' + r.h +
           '" fill="' + r.color + '" fill-opacity="0.20" stroke="' + r.color +
           '" stroke-width="0.6" stroke-opacity="0.9" />';
    if (r.h >= 14) {
      svg += '<text x="' + (midX) + '" y="' + (r.y + r.h / 2 + 3) +
             '" font-size="10" fill="' + r.color + '" font-weight="600" text-anchor="middle" ' +
             'font-family="ui-monospace, monospace">' + _esc(r.label) + '</text>';
    }
  }
  svg += '</svg>';
  return svg;
}

// -----------------------------------------------------------------------------
// Renderer: synteny section (sits below the focus panel)
// -----------------------------------------------------------------------------

function _csBuildSyntenySection() {
  const blocks = _csGetSyntenyBlocks();
  if (!blocks) {
    return '<div class="cs-synteny-empty">' +
      'Synteny analysis requires <code>cs_breakpoints_v1</code> schema_version \u22652. ' +
      'Re-run <code>STEP_CS01_extract_breakpoints.py</code> to refresh the JSON.' +
      '</div>';
  }
  const synteny = _csComputeSynteny();
  if (!synteny) return '';
  // Sankey
  const sankey = _csBuildSankeySvg();
  // Per-chrom relationship table
  const tableRows = _csNaturalSortChroms(Object.keys(synteny.gar_partners)).map(c => {
    const partners = synteny.gar_partners[c].filter(p => p.dominant);
    const partnerStr = partners.map(p => _esc(p.mac_chr) + ' (' + (p.frac * 100).toFixed(0) + '%)').join(', ');
    const cls = synteny.gar_class[c] || 'unknown';
    const def = CS_REL_DEF[cls] || { label: cls, color: 'var(--ink-dim)' };
    return '<tr>' +
      '<td><code>' + _esc(c) + '</code></td>' +
      '<td>' + (partnerStr || '\u2014') + '</td>' +
      '<td><span class="cs-rel-pill" style="color:' + def.color + ';border-color:' + def.color + ';">' +
      _esc(def.label) + '</span></td>' +
      '</tr>';
  }).join('');
  // Inversion contexts
  const ctxs = _csInversionContexts();
  let inversionTable = '';
  if (ctxs && ctxs.length > 0) {
    const fmtDist = (d) => {
      if (d == null) return '\u2014';
      if (d < 1000) return d + ' bp';
      if (d < 1e6)  return (d / 1000).toFixed(1) + ' kb';
      return (d / 1e6).toFixed(2) + ' Mb';
    };
    const rows = ctxs.map(ctx => {
      const cls = ctx.chrom_class;
      const def = CS_REL_DEF[cls] || { label: cls, color: 'var(--ink-dim)' };
      const isInside = ctx.inside_single_block;
      const insideStr = isInside ? 'inside' : 'spans ' + ctx.spanning_n_blocks;
      // Show left + right distances; bold the closer one. The permutation
      // test on dist_min_edge_bp (= min(left, right)) is what the global
      // enrichment test uses.
      const dl = ctx.dist_left_edge_bp,  dr = ctx.dist_right_edge_bp;
      const lTag = (dl != null && dr != null && dl <= dr) ? 'b' : 'span';
      const rTag = (dl != null && dr != null && dr <  dl) ? 'b' : 'span';
      return '<tr>' +
        '<td><code>' + _esc(ctx.candidate_id) + '</code></td>' +
        '<td><code>' + _esc(ctx.chrom) + '</code></td>' +
        '<td>' + (ctx.start_bp / 1e6).toFixed(2) + '\u2013' + (ctx.end_bp / 1e6).toFixed(2) + ' Mb</td>' +
        '<td><' + lTag + '>' + fmtDist(dl) + '</' + lTag + '> / ' +
                  '<' + rTag + '>' + fmtDist(dr) + '</' + rTag + '></td>' +
        '<td>' + insideStr + '</td>' +
        '<td><span class="cs-rel-pill" style="color:' + def.color + ';border-color:' + def.color + ';">' +
        _esc(def.label) + '</span></td>' +
        '<td class="cs-interp">' + _esc(ctx.interpretation) + '</td>' +
        '</tr>';
    }).join('');
    inversionTable =
      '<div class="cs-synteny-subtitle">Inversion candidates in synteny context (' + ctxs.length + ')</div>' +
      '<table class="cs-synteny-table">' +
      '<thead><tr>' +
        '<th>candidate</th><th>chrom</th><th>position</th>' +
        '<th>L / R bp distance to edge</th><th>blocks</th>' +
        '<th>chrom rel.</th><th>interpretation</th>' +
      '</tr></thead><tbody>' + rows + '</tbody></table>' +
      '<button class="cs-perm-btn" id="csPermRunBtn" title="Run a permutation test (' +
      CS_SYNTENY_PERM_N + ' draws per candidate, length-matched random intervals on the same chrom) ' +
      'to assess whether inversion breakpoints are closer to synteny edges than expected by chance.">' +
      'Run permutation test (' + CS_SYNTENY_PERM_N + ' draws \u00B7 length-matched)</button>' +
      '<div id="csPermResult" class="cs-perm-result" style="display:none;"></div>';
  } else {
    inversionTable = '<div class="cs-synteny-subtitle">Inversion candidates in synteny context</div>' +
      '<div class="cs-synteny-empty-small">No promoted inversion candidates in <code>state.candidateList</code> ' +
      'on chromosomes covered by the synteny blocks. Promote candidates on page 3 (catalogue) ' +
      'or page 1 (local PCA |z|) to populate this section.</div>';
  }

  // Legend
  const legend = '<div class="cs-synteny-legend">' +
    Object.keys(CS_REL_DEF).map(k => {
      const d = CS_REL_DEF[k];
      return '<span class="cs-rel-pill" style="color:' + d.color + ';border-color:' + d.color + ';">' +
        _esc(d.label) + '</span>';
    }).join(' ') + '</div>';

  return '<div class="cs-synteny">' +
    '<div class="cs-synteny-header">' +
      '<h3 class="cs-synteny-title">Macro-synteny (Cgar \u2194 Cmac)</h3>' +
      '<div class="cs-synteny-subhint">Derived from <b>' + blocks.length + '</b> wfmash synteny blocks ' +
      '(\u22655% chrom-aligned-bp threshold). Hover Sankey ribbons for per-flow stats.</div>' +
    '</div>' +
    legend +
    '<div class="cs-synteny-sankey">' + sankey + '</div>' +
    '<div class="cs-synteny-subtitle">Per-chromosome relationships</div>' +
    '<table class="cs-synteny-table">' +
    '<thead><tr><th>Cgar chrom</th><th>Cmac partners (\u22655% aligned)</th><th>relationship</th></tr></thead>' +
    '<tbody>' + tableRows + '</tbody></table>' +
    inversionTable +
    '</div>';
}

// Render the synteny section into its slot. Wires up the permutation
// button after the DOM is in place.
function _renderCrossSpeciesSynteny() {
  const slot = document.getElementById('csSynteny');
  if (!slot) return;
  // Render only when we have v2 data
  const blocks = _csGetSyntenyBlocks();
  if (!blocks || !state.crossSpecies) {
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  slot.style.display = 'block';
  slot.innerHTML = _csBuildSyntenySection();
  // Wire permutation button
  const btn = slot.querySelector('#csPermRunBtn');
  if (btn) {
    btn.addEventListener('click', () => {
      btn.disabled = true;
      btn.textContent = 'running...';
      // Defer to next tick so the UI updates before the synchronous compute
      setTimeout(() => {
        let result;
        try {
          result = _csPermutationTest();
        } catch (e) {
          console.warn('[crossSpecies] perm test failed:', e);
        }
        const out = slot.querySelector('#csPermResult');
        if (out && result) {
          out.style.display = 'block';
          out.innerHTML = _csBuildPermResultHtml(result);
        }
        btn.textContent = 'Re-run permutation test';
        btn.disabled = false;
      }, 0);
    });
  }
}

// =============================================================================
// turn 115: Cross-species dot plot wiring.
// -----------------------------------------------------------------------------
// Calls window.popgenDotplot.makeDotplotPanel(opts) to render the dot plot
// panel into #csDotplot. Two data sources are supported, both optional:
//
//   1. wfmash synteny_blocks — already inside cs_breakpoints_v1 schema_v2;
//      reads from state.crossSpecies. No new file format. The default
//      resolution when only this is loaded is 'wfmash'.
//
//   2. mashmap multi-resolution — a separate dotplot_mashmap_v1.json layer
//      stored on state.dotplotMashmap. Pre-computed once on the cluster
//      across multiple segment-length × percent-identity tiers and dropped
//      into the atlas via the file picker. Recognized by the standard
//      detector chain (_isDotplotMashmapJSON below).
//
// At least ONE of those two sources must be present for the panel to render.
// If neither is loaded, the #csDotplot wrapper stays hidden — same fail-soft
// pattern as #csSynteny.
//
// Cached panel: state._csDotplotPanel holds the live HTMLElement. Re-render
// only replaces the panel when the underlying data sources change (new
// cs_breakpoints load, new mashmap load), so the user's "pinned" popup state
// survives unrelated re-renders of page 16. Cache key is a small fingerprint
// of the data shape.
// =============================================================================
function _renderCrossSpeciesDotplot() {
  const slot = document.getElementById('csDotplot');
  if (!slot) return;
  if (typeof window === 'undefined' || !window.popgenDotplot ||
      typeof window.popgenDotplot.makeDotplotPanel !== 'function') {
    // Module not loaded yet (or load failed). Hide gracefully.
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  const cs = state.crossSpecies;
  const mm = state.dotplotMashmap;
  const haveWfmash = !!(cs && Array.isArray(cs.synteny_blocks) && cs.synteny_blocks.length > 0);
  const haveMashmap = !!(mm && Array.isArray(mm.resolutions) && mm.resolutions.length > 0);
  if (!haveWfmash && !haveMashmap) {
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  // Build a small fingerprint to know when to actually re-render.
  const fp = [
    haveWfmash ? ('wf:' + cs.n_synteny_blocks + ':' + (cs.loaded_at || '')) : 'wf:0',
    haveMashmap ? ('mm:' + mm.resolutions.length + ':' + (mm.generated_at || '')) : 'mm:0',
  ].join('|');
  if (state._csDotplotPanelFp === fp && state._csDotplotPanel &&
      slot.firstChild === state._csDotplotPanel) {
    slot.style.display = 'block';
    return;   // already mounted, data unchanged → leave the panel alone
  }
  // Build a fresh data object for the dotplot module — adapt cs_breakpoints
  // synteny_blocks shape if present, and pass through mashmap resolutions
  // as-is (they were emitted in the unified shape).
  const dataObj = {
    species_query:        (cs && cs.species_query)  || (mm && mm.species_query)  || { name: 'query' },
    species_target:       (cs && cs.species_target) || (mm && mm.species_target) || { name: 'target' },
    chrom_lengths_query:  (cs && cs.chrom_lengths_query)  || (mm && mm.chrom_lengths_query)  || {},
    chrom_lengths_target: (cs && cs.chrom_lengths_target) || (mm && mm.chrom_lengths_target) || {},
    wfmash_blocks:        haveWfmash ? cs.synteny_blocks : null,
    mashmap_resolutions:  haveMashmap ? mm.resolutions : null,
  };
  let panel;
  try {
    panel = window.popgenDotplot.makeDotplotPanel({
      data: dataObj,
      defaultResolution: haveWfmash ? 'wfmash'
                                    : (mm.resolutions[0] && mm.resolutions[0].label),
      miniSize: 300,
      enlargedSize: 720,
    });
  } catch (e) {
    console.warn('[crossSpecies] dotplot render failed:',
                 e && e.message ? e.message : e);
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  slot.innerHTML = '';
  slot.appendChild(panel);
  slot.style.display = 'block';
  state._csDotplotPanel = panel;
  state._csDotplotPanelFp = fp;
}

// -----------------------------------------------------------------------------
// turn 117: cross-species focal-vs-background widget renderer.
//
// Mounts popgenFocalVsBg.makeFocalVsBgPanel into #csFocalVsBg with a focal
// interval centered on the active cs_breakpoint's Gar coordinates ±
// state._focalVsBg.csRadiusBp (default 2 Mb, Spalax-matching).
//
// The widget reads per-window metric values from `state.data` (the 226-fish
// cohort precomp). This means the widget is biologically meaningful only
// when the active breakpoint's gar_chr matches the cohort's currently-loaded
// chrom — otherwise we'd be filtering LG12 windows by LG28 Mb coords. So we
// hide the panel when:
//   - the popgenFocalVsBg module hasn't loaded
//   - no breakpoint is active
//   - bp.gar_chr does not match state.data.chrom (cohort loaded for a
//     different chromosome)
// Fail-soft on render errors (try/catch).
// -----------------------------------------------------------------------------
function _renderCrossSpeciesFocalVsBg() {
  const slot = document.getElementById('csFocalVsBg');
  if (!slot) return;
  // Initialise state slot if missing (handoff: do this first thing in the
  // renderer so any code path that reads state._focalVsBg sees a valid object).
  state._focalVsBg = state._focalVsBg || { csRadiusBp: 2000000 };
  // Module not loaded → hide and bail
  if (typeof window === 'undefined' || !window.popgenFocalVsBg ||
      typeof window.popgenFocalVsBg.makeFocalVsBgPanel !== 'function') {
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  // No active breakpoint → hide
  const bp = (typeof _csFindActive === 'function') ? _csFindActive() : null;
  if (!bp) {
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  // Cohort chrom mismatch — widget reads per-window metrics from state.data
  // (the 226-fish cohort, single haplotype). If the user is browsing a
  // breakpoint on LG28 but has LG12 cohort data loaded, the widget can't
  // sensibly compute z / fst / etc. for the focal interval. Show a brief
  // hint instead of bogus stats.
  const cohortChrom = (state.data && state.data.chrom) || null;
  if (!cohortChrom || cohortChrom !== bp.gar_chr) {
    slot.style.display = 'block';
    slot.innerHTML =
      '<div style="background:var(--panel-2);border:1px solid var(--rule);' +
      'border-radius:4px;padding:14px 16px;font-family:var(--mono);' +
      'font-size:11px;color:var(--ink-dim);line-height:1.6;">' +
      'focal vs background \u00B7 ' +
      '<span style="color:var(--ink);">' + _esc(bp.gar_chr || '?') + '</span> ' +
      'breakpoint @ ' + (bp.gar_pos_mb || 0).toFixed(2) + ' Mb' +
      '<div style="margin-top:5px;color:var(--ink-dimmer);font-size:10px;">' +
      'cohort chrom in scrubber: <code>' + _esc(cohortChrom || '(none loaded)') +
      '</code> \u2014 load a precomp for <code>' + _esc(bp.gar_chr || '?') +
      '</code> on page 1 to compute focal-vs-background statistics here.' +
      '</div></div>';
    return;
  }
  // Build focal range from breakpoint Gar coords ± csRadiusBp. The breakpoint
  // span [gar_pos_start, gar_pos_end] is typically narrow (synteny edge); use
  // its midpoint as the anchor and expand symmetrically.
  const radius = (state._focalVsBg && Number.isFinite(state._focalVsBg.csRadiusBp)
                  && state._focalVsBg.csRadiusBp > 0) ? state._focalVsBg.csRadiusBp : 2000000;
  const anchorBp = ((bp.gar_pos_start || 0) + (bp.gar_pos_end || bp.gar_pos_start || 0)) / 2;
  const chromLen = (state.data && state.data.n_bp) ||
                   (typeof _csIdeogramChromLength === 'function'
                     ? _csIdeogramChromLength('gar', bp.gar_chr) : null);
  let focal_lo_bp = Math.max(0, anchorBp - radius);
  let focal_hi_bp = anchorBp + radius;
  if (chromLen != null && Number.isFinite(chromLen) && chromLen > 0) {
    focal_hi_bp = Math.min(chromLen, focal_hi_bp);
  }
  // Anchor marks: dashed lines at the start/end of the breakpoint span (or
  // just the midpoint if start==end). Drawn by the widget's scatter renderer.
  const anchor_marks_mb = [];
  if (Number.isFinite(bp.gar_pos_start)) anchor_marks_mb.push(bp.gar_pos_start / 1e6);
  if (Number.isFinite(bp.gar_pos_end) && bp.gar_pos_end !== bp.gar_pos_start) {
    anchor_marks_mb.push(bp.gar_pos_end / 1e6);
  }
  let panel;
  try {
    panel = window.popgenFocalVsBg.makeFocalVsBgPanel({
      page_id:           'page16',
      anchor_id:         bp.id,
      anchor_type:       'cs_breakpoint',
      chrom:             bp.gar_chr,
      focal_lo_bp:       focal_lo_bp,
      focal_hi_bp:       focal_hi_bp,
      anchor_marks_mb:   anchor_marks_mb,
      radius_source:     'self',
      radius_default_bp: radius,
      get_radius_bp:     () => state._focalVsBg.csRadiusBp,
      set_radius_bp:     (r) => { state._focalVsBg.csRadiusBp = r; },
      on_radius_change:  () => {
        // Radius changed inside the widget: re-render this same renderer so
        // focal_lo_bp / focal_hi_bp recompute from the new radius. Other
        // page-16 panels don't care about the focal-vs-bg radius.
        try { _renderCrossSpeciesFocalVsBg(); } catch (_) {}
      },
      cohort_selector:   null,           // day 1: single cohort, no picker
    });
  } catch (e) {
    console.warn('[crossSpecies] focal-vs-bg render failed:',
                 e && e.message ? e.message : e);
    slot.style.display = 'none';
    slot.innerHTML = '';
    return;
  }
  slot.innerHTML = '';
  slot.appendChild(panel);
  slot.style.display = 'block';
}

// NOTE: dotplot_mashmap_v1 constants and IO loader (legacy 26023-26124) live in
// page16b.js — page16 only reads state.dotplotMashmap as an optional overlay.


// =============================================================================
// BLOCK C — _csBuildPermResultHtml (legacy 28367-28419)
// Note: positioned in legacy after the multi-species block by accident of patch
// order. Belongs logically with the cross-species permutation test.
// =============================================================================
function _csBuildPermResultHtml(result) {
  const g = result.global;
  const enr = g.enrichment_ratio != null ? g.enrichment_ratio.toFixed(2) : '\u2014';
  const p   = g.binomial_p != null ? (g.binomial_p < 0.001 ? g.binomial_p.toExponential(2)
                                                            : g.binomial_p.toFixed(3)) : '\u2014';
  const interp = (g.binomial_p != null && g.binomial_p < 0.05 && g.enrichment_ratio > 1)
    ? 'Inversion breakpoints are SIGNIFICANTLY closer to synteny edges than expected by chance. ' +
      'This supports the "ancient/fragile architecture reuse" interpretation.'
    : (g.binomial_p != null && g.binomial_p < 0.05 && g.enrichment_ratio < 1)
    ? 'Inversion breakpoints are SIGNIFICANTLY FARTHER from synteny edges than expected. ' +
      'Most candidates are inside conserved blocks — consistent with recent within-lineage polymorphism.'
    : 'No significant deviation from random expectation. Mixed signal — ' +
      'individual breakpoints may still be informative; inspect per-candidate columns.';
  // Per-candidate quick-view: top 10 by ascending p
  const topByP = result.per_candidate.slice()
    .filter(r => r.p_value != null)
    .sort((a, b) => a.p_value - b.p_value)
    .slice(0, 10);
  const rows = topByP.map(r => {
    const sig = (r.p_value != null && r.p_value < 0.05) ? '*' : '';
    return '<tr>' +
      '<td><code>' + _esc(r.candidate_id) + '</code></td>' +
      '<td><code>' + _esc(r.chrom) + '</code></td>' +
      '<td>' + (r.dist_min_edge_bp != null
        ? (r.dist_min_edge_bp / 1e6).toFixed(3) + ' Mb' : '\u2014') + '</td>' +
      '<td>' + (r.expected_dist_median != null
        ? (r.expected_dist_median / 1e6).toFixed(3) + ' Mb' : '\u2014') + '</td>' +
      '<td>' + (r.percentile != null ? (r.percentile * 100).toFixed(1) + '%' : '\u2014') + '</td>' +
      '<td>' + (r.p_value != null
        ? (r.p_value < 0.001 ? r.p_value.toExponential(2) : r.p_value.toFixed(3)) : '\u2014') +
      '<span class="cs-perm-sig">' + sig + '</span></td>' +
      '</tr>';
  }).join('');
  return '<div class="cs-perm-summary">' +
    '<div><b>Global enrichment test</b> (n=' + g.n_valid +
    ' candidates, ' + g.n_perm.toLocaleString() + ' draws each, quantile=' +
    (g.quantile * 100).toFixed(0) + '%):</div>' +
    '<div class="cs-perm-stats">' +
      '<span>observed near-edge: <b>' + g.near_edge_observed + '</b> / ' + g.n_valid + '</span>' +
      ' \u00B7 <span>expected: ' + g.near_edge_expected.toFixed(1) + '</span>' +
      ' \u00B7 <span>enrichment ratio: <b>' + enr + '</b></span>' +
      ' \u00B7 <span>binomial p: <b>' + p + '</b></span>' +
    '</div>' +
    '<div class="cs-perm-interp">' + _esc(interp) + '</div>' +
    '</div>' +
    '<div class="cs-synteny-subtitle">Top 10 candidates by permutation p-value</div>' +
    '<table class="cs-synteny-table">' +
    '<thead><tr><th>candidate</th><th>chrom</th><th>obs dist</th>' +
    '<th>median random dist</th><th>percentile</th><th>p-value</th></tr></thead>' +
    '<tbody>' + rows + '</tbody></table>';
}


