// =============================================================================
// inversion_review/page11.js — "4 boundaries" tab
// =============================================================================
// Stage:        refinement
// Legacy DOM:   <div id="page11"> (legacy lines 7906–7921)
// Renderer:     renderBoundariesPage()  (legacy lines 30175–30314)
// Hotkey hub:   _bndKeyHandler / _bndAttachHotkeys / _bndDetachHotkeys
//                (legacy lines 30318–30345)
//
// What this page does
// -------------------
// Refines a promoted candidate's [start_bp, end_bp] into approximate left/right
// boundary zones using multiple evidence tracks. Auto-propose runs the
// BOUNDARY_TRACK_WEIGHTS-weighted algorithm; manual override via E/F at the
// scrubber cursor. Save (B) persists onto the candidate; Reset (R) clears.
//
// Terminology contract (kept from legacy comments):
//   "boundary zone"   = default verdict; what auto-propose ever produces.
//   "exact breakpoint" = reserved for junction-level evidence only.
//
// Extraction notes (Batch 2)
// --------------------------
// Only the page-level entry points are extracted here. The ~25 boundary
// helpers (`_bnd*`, `_ensureBoundariesState`, `_findSVAnchorsInZone`,
// `_buildBoundaryTrackScores`, `_computeBoundaryEdges`, `_buildBoundaryRecord`,
// `_renderCandidateNavInline`, the four `_wireRepeatDensity*` wirers, etc.)
// are marked TODO_MISSING below — the merge chat decides whether to promote
// them to shared/ or copy them in here.
// =============================================================================

import { SLOT_REGISTRY } from '../shared/state.js';

// -----------------------------------------------------------------------------
// Unresolved external deps — to be filled by the merge chat.
// -----------------------------------------------------------------------------
//
// TODO_MISSING(_ensureBoundariesState)         — legacy line 17763
// TODO_MISSING(_boundaryScanRange)             — legacy line 17813
// TODO_MISSING(_findSVAnchorsInZone)           — legacy line 17895
// TODO_MISSING(_buildBoundaryTrackScores)      — legacy line 17942
// TODO_MISSING(_computeBoundaryEdges)          — legacy line 18127
// TODO_MISSING(_buildBoundaryRecord)           — legacy line 18256
// TODO_MISSING(_bndFindCandidate)              — legacy line 18297
// TODO_MISSING(_bndFmtBp)                      — legacy line 18308
// TODO_MISSING(_bndCloneRecord)                — legacy line 18317
// TODO_MISSING(_bndStageFromCandidate)         — legacy line 18337
// TODO_MISSING(_bndPopulateCandidateSelect)    — legacy line 18351
// TODO_MISSING(_bndUpdateRadiusButtons)        — legacy line 18383
// TODO_MISSING(_bndUpdateSaveButton)           — legacy line 18394
// TODO_MISSING(_bndUpdateStatusSelect)         — legacy line 18404
// TODO_MISSING(_bndSelectCandidate)            — legacy line 18415
// TODO_MISSING(_bndAutoPropose)                — legacy line 18424
// TODO_MISSING(_bndManualOverride)             — legacy line 18479
// TODO_MISSING(_bndOverrideLeft)               — legacy line 18622
// TODO_MISSING(_bndOverrideRight)              — legacy line 18623
// TODO_MISSING(_bndReset)                      — legacy line 18626
// TODO_MISSING(_bndSave)                       — legacy line 18640
// TODO_MISSING(_bndRefreshSummary)             — legacy line 18658
// TODO_MISSING(_bndDrawTracks)                 — legacy line 18726
// TODO_MISSING(_bndRefreshUI)                  — legacy line 18911
// TODO_MISSING(_renderBndFocalVsBg)            — legacy line 20507
// TODO_MISSING(_wireRepeatDensityEscapeReset)  — legacy line 20583
// TODO_MISSING(_wireRepeatDensityClassCycle)   — legacy line 20856
// TODO_MISSING(_wireRepeatDensityAllTeToggle)  — legacy line 20895
// TODO_MISSING(_wireRepeatDensityArrowNav)     — legacy line 20941
// TODO_MISSING(_renderCandidateNavInline)      — shared candidate-nav helper
//                                                (used by pages 4/6/7/11 too)
//
// All of the `_bnd*` helpers operate on a single closure-scoped `state`
// object (the legacy global). When porting, they should accept `state` as
// an explicit first argument. For now this module presumes the merge chat
// promotes a `state` symbol via the wrapper convention used in main.js.
// -----------------------------------------------------------------------------


// =============================================================================
// renderBoundariesPage — public entry, called by tab dispatcher
// =============================================================================
// Legacy: lines 30175–30314. Verbatim except for one change:
//   - all bare `_bnd*`, `_renderCandidateNavInline`, `_ensureBoundariesState`
//     references stay as global lookups; the merge chat patches them.
// =============================================================================
export function renderBoundariesPage() {
  const slot = document.getElementById('page11Content');
  if (!slot) return;
  const bs = _ensureBoundariesState();
  // v3.99 turn 13 ask 1: candidate-nav bar above the existing toolbar. Always
  // re-rendered so prev/next button state stays in sync. Boundaries page also
  // has its own #bndCandSelect dropdown, but the prev/next buttons are
  // faster for stepping through; both stay in sync via _navigateToCandidate.
  const page11 = document.getElementById('page11');
  if (page11) {
    const oldNav = page11.querySelector('.cand-nav-inline');
    if (oldNav) oldNav.remove();
    const navBar = _renderCandidateNavInline({ idPrefix: 'bnd' });
    navBar.style.margin = '12px 32px 0';
    // Insert AFTER the page11Header so the nav sits between title and content
    const header = document.getElementById('page11Header');
    if (header && header.nextSibling) {
      page11.insertBefore(navBar, header.nextSibling);
    } else {
      page11.appendChild(navBar);
    }
  }
  // Build toolbar + structure once; populate dynamic parts each time
  if (!slot.__bndBuilt) {
    slot.innerHTML = `
      <div class="bnd-toolbar">
        <span class="bnd-lbl">candidate:</span>
        <select id="bndCandSelect"></select>
        <span class="bnd-lbl">scan radius:</span>
        <div class="bnd-radius-group">
          <button class="bnd-radius-btn" data-radius="1000"
                  title="±1 kb scan — for SV-supported breakpoints with split-read evidence; most PCA tracks won't render meaningfully at this scale">1 kb</button>
          <button class="bnd-radius-btn" data-radius="5000"
                  title="±5 kb scan — fine breakpoint inspection; matches MODULE_3 win5000.step1000 θπ scale">5 kb</button>
          <button class="bnd-radius-btn" data-radius="10000"
                  title="±10 kb scan — boundary zone refinement; matches MODULE_3 win10000.step2000 θπ scale">10 kb</button>
          <button class="bnd-radius-btn" data-radius="25000"
                  title="±25 kb scan — typical breakpoint working scale for boundary_zone_only verdicts">25 kb</button>
          <button class="bnd-radius-btn" data-radius="100000"
                  title="±100 kb scan — broader boundary context; matches the hobs_profile / fst_profile distances_kb default range in SCHEMA §21">100 kb</button>
          <span class="bnd-radius-sep"
                style="display: inline-block; width: 1px; height: 16px;
                       background: var(--rule); margin: 0 4px;
                       vertical-align: middle;"></span>
          <button class="bnd-radius-btn" data-radius="1000000"
                  title="±1 Mb scan — whole-region context, useful when SV calls and PCA boundaries disagree by hundreds of kb">1 Mb</button>
          <button class="bnd-radius-btn" data-radius="1500000"
                  title="±1.5 Mb scan">1.5 Mb</button>
          <button class="bnd-radius-btn" data-radius="2000000"
                  title="±2 Mb scan">2 Mb</button>
          <button class="bnd-radius-btn" data-radius="5000000"
                  title="±5 Mb scan — broadest context, useful for very large or hugely-asymmetric candidates (kicks in automatically at &gt;CANDIDATE_HUGE_BP via _boundaryScanRange)">5 Mb</button>
        </div>
        <button id="bndAutoProposeBtn" class="primary" title="Run the auto-propose algorithm (or hotkey A)">auto-propose</button>
        <button id="bndOverrideLBtn" title="Override LEFT boundary at the current focal window (hotkey E)">override left (E)</button>
        <button id="bndOverrideRBtn" title="Override RIGHT boundary at the current focal window (hotkey F)">override right (F)</button>
        <button id="bndResetBtn" title="Clear both boundaries (hotkey R)">reset (R)</button>
        <button id="bndSaveBtn" title="Persist boundary annotation onto the candidate (hotkey B)">save (B)</button>
        <span class="bnd-status-wrap">
          <span class="bnd-lbl">status:</span>
          <select id="bndStatusSel" title="Strongest available evidence — auto only ever sets boundary_zone_only; promote manually after reviewing SV/junction evidence.">
            <option value="boundary_zone_only">boundary_zone_only</option>
            <option value="SV_supported">SV_supported</option>
            <option value="junction_supported">junction_supported</option>
          </select>
        </span>
      </div>
      <div class="bnd-info" id="bndInfo">No candidate selected.</div>
      <div class="bnd-tracks" id="bndTracks"></div>
      <div class="bnd-summary" id="bndSummary"></div>
      <!-- v4 turn 25: class summary auto-scan -->
      <div class="bnd-class-summary" id="bndClassSummary"></div>
      <!-- v4 turn 18: repeat density panel -->
      <div class="bnd-repeat-density" id="bndRepeatDensity"></div>
      <!-- turn 116: ncRNA density panel (rRNA / tRNA / ncRNA-other) -->
      <div class="bnd-ncrna-density" id="bndNcRNADensity"></div>
      <!-- turn 117: focal-vs-background widget -->
      <div class="bnd-focal-vs-bg" id="bndFocalVsBg"></div>
    `;
    // Wire events
    const candSel = document.getElementById('bndCandSelect');
    if (candSel) {
      candSel.addEventListener('change', () => {
        const v = candSel.value;
        const id = (v === '' || v == null) ? null : (Number.isFinite(parseInt(v, 10)) ? parseInt(v, 10) : v);
        _bndSelectCandidate(id);
      });
    }
    const radiusBtns = document.querySelectorAll('#page11 .bnd-radius-btn');
    radiusBtns.forEach(b => {
      b.addEventListener('click', () => {
        const r = parseInt(b.getAttribute('data-radius'), 10);
        if (Number.isFinite(r) && r > 0) {
          const bs2 = _ensureBoundariesState();
          bs2.scan_radius_bp = r;
          if (bs2.active_cand_id != null) bs2.cache.delete(bs2.active_cand_id);
          _bndRefreshUI();
        }
      });
    });
    const autoBtn = document.getElementById('bndAutoProposeBtn');
    if (autoBtn) autoBtn.addEventListener('click', _bndAutoPropose);
    const lBtn = document.getElementById('bndOverrideLBtn');
    if (lBtn) lBtn.addEventListener('click', _bndOverrideLeft);
    const rBtn = document.getElementById('bndOverrideRBtn');
    if (rBtn) rBtn.addEventListener('click', _bndOverrideRight);
    const resBtn = document.getElementById('bndResetBtn');
    if (resBtn) resBtn.addEventListener('click', _bndReset);
    const saveBtn = document.getElementById('bndSaveBtn');
    if (saveBtn) saveBtn.addEventListener('click', _bndSave);
    const statSel = document.getElementById('bndStatusSel');
    if (statSel) {
      statSel.addEventListener('change', () => {
        const bs2 = _ensureBoundariesState();
        bs2.staging.breakpoint_status = statSel.value;
        bs2.staging.dirty = true;
        _bndUpdateSaveButton();
      });
    }
    slot.__bndBuilt = true;
  }
  _bndPopulateCandidateSelect();
  _bndRefreshUI();
}


// =============================================================================
// Hotkey wiring — E (override left) / F (override right) / B (save)
//                  R (reset) / A (auto-propose). Only fires when page11 active.
// Legacy: lines 30318–30345.
// =============================================================================
let _bndKeyHandlerAttached = false;

export function _bndKeyHandler(e) {
  const visible = document.getElementById('page11');
  if (!visible || !visible.classList || !visible.classList.contains('active')) return;
  const ae = document.activeElement;
  if (ae && /^(INPUT|TEXTAREA|SELECT)$/.test(ae.tagName)) return;
  // Don't trigger on Ctrl/Meta/Alt combos (those are reserved for browser shortcuts)
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  const k = (e.key || '').toLowerCase();
  if      (k === 'e') { e.preventDefault(); _bndOverrideLeft(); }
  else if (k === 'f') { e.preventDefault(); _bndOverrideRight(); }
  else if (k === 'b') { e.preventDefault(); _bndSave(); }
  else if (k === 'r') { e.preventDefault(); _bndReset(); }
  else if (k === 'a') { e.preventDefault(); _bndAutoPropose(); }
}

export function _bndAttachHotkeys() {
  if (_bndKeyHandlerAttached) return;
  if (typeof document === 'undefined') return;
  document.addEventListener('keydown', _bndKeyHandler);
  _bndKeyHandlerAttached = true;
}

export function _bndDetachHotkeys() {
  if (!_bndKeyHandlerAttached) return;
  if (typeof document === 'undefined') return;
  document.removeEventListener('keydown', _bndKeyHandler);
  _bndKeyHandlerAttached = false;
}


// =============================================================================
// state references not in SLOT_REGISTRY (caught at extraction time)
// =============================================================================
//
// state.repeatDensity   — per-chrom transposable-element density layer; used
//                         by `_wireRepeatDensity*` wirers and the bnd repeat
//                         density panel. Lazily created by data loader.
// state.ncRNADensity    — per-chrom ncRNA (rRNA/tRNA) density layer; same
//                         pattern as repeatDensity.
// state.candidate.*     — IS in SLOT_REGISTRY, but the `.locked_labels`,
//                         `.K`, `.tracks`, `.ref_l2`, `.start_bp`, `.end_bp`,
//                         `.id` sub-fields are intra-candidate schema (see
//                         shared/state_io.js KNOWN_LAYERS for canonical shape).
// state.candidateList   — IS in SLOT_REGISTRY (cross_atlas array); used by
//                         _bndPopulateCandidateSelect.
// state.data            — IS in SLOT_REGISTRY (transient); used as
//                         state.data.final_classification consumer (page4 too).
//
// (No genuine SLOT_REGISTRY-missing slots in renderBoundariesPage itself —
// all the cache/anchor business is keyed off the closure-scoped `bs` object
// from _ensureBoundariesState, not state.)
// =============================================================================
