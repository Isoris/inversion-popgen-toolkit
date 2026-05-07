// Atlas/inversion_catalogue/page9.js
// =============================================================================
// page9 — Confirmed candidates carousel
// (`<div id="page9">` contains the prev/next nav bar and meta slot)
//
// Source: legacy/Inversion_atlas.html lines 7782–7812 (HTML shell only)
//
// IMPORTANT: A grep of the legacy file reveals the page9 carousel JS does
// NOT EXIST in legacy/Inversion_atlas.html. The HTML shell is wired with
// IDs (#confirmedNavBar, #confirmedNavPrev, #confirmedNavNext,
// #confirmedNavInfo, #confirmedCandidateMeta, #confirmedEmpty) but no
// JavaScript file in the legacy drop populates or wires them. The page is
// effectively a stub even in the legacy build — pressing the page9 tab
// shows the empty-state message at #confirmedEmpty.
//
// Confirmed by:
//   $ grep -n 'confirmedNav' legacy/Inversion_atlas.html
//   7783: <div id="confirmedNavBar" ...
//   7787: <button id="confirmedNavPrev" ...
//   7793: <button id="confirmedNavNext" ...
//   (all hits are HTML only — no JS handlers)
//
// What the page is supposed to do (per the page-tab tooltip at legacy
// line 5076 and the empty-state copy at legacy lines 7802–7811):
//   "Carousel walk-through of all candidates marked confirmed on page 2."
//   - state.candidateList.filter(c => c.confirmed === true)
//   - prev/next buttons cycle through the confirmed set
//   - reuses the page2 candidate-focus rendering for the current candidate
//
// External dependencies (when wired up):
//   TODO_MISSING(_renderConfirmedCarousel)
//     — full carousel renderer (does not exist yet in legacy)
//   TODO_MISSING(_wireConfirmedCarouselNav)
//     — keydown ←/→ + button click handlers (does not exist yet)
//   TODO_MISSING(renderCandidateFocus)  — owned by Batch 1 (page2)
//     — the page9 carousel reuses page2's candidate-focus renderer for
//       the currently-displayed confirmed candidate
//   global `state`                — reads state.candidateList,
//                                   state.confirmedCarouselIndex (new)
//
// Decision for this batch: ship a no-op shell with the public entry,
// matching the legacy behaviour (showing the empty-state message). The
// merge chat — or a follow-up batch — implements the real carousel.
// =============================================================================

const state = (typeof window !== 'undefined' && window.state) ? window.state : {};

/**
 * Public entry: refresh the confirmed-candidates carousel.
 *
 * Current behaviour (matches legacy stub): show #confirmedEmpty, hide
 * #confirmedNavBar and #confirmedCandidateMeta. This is a placeholder; the
 * full carousel is TODO_MISSING.
 */
export function refreshConfirmedCarousel() {
  if (typeof document === 'undefined') return;
  const navBar = document.getElementById('confirmedNavBar');
  const meta   = document.getElementById('confirmedCandidateMeta');
  const empty  = document.getElementById('confirmedEmpty');

  // TODO_MISSING(_renderConfirmedCarousel) — full carousel implementation.
  // For now: count confirmed candidates. If zero, keep the empty message.
  // If non-zero, still keep the empty message because the carousel renderer
  // is not yet implemented; the merge chat / a follow-up batch ships it.
  const confirmedCount = (Array.isArray(state.candidateList))
    ? state.candidateList.filter(c => c && c.confirmed === true).length
    : 0;

  if (confirmedCount === 0) {
    if (navBar) navBar.style.display = 'none';
    if (meta)   meta.style.display = 'none';
    if (empty)  empty.style.display = 'block';
    return;
  }

  // Non-zero confirmed candidates exist — but the carousel logic is not
  // ported yet. Show a placeholder telling the user where to look.
  if (empty) {
    empty.style.display = 'block';
    empty.innerHTML =
      '<div style="font-size:14px; margin-bottom:14px;">' +
        confirmedCount + ' confirmed candidate' +
        (confirmedCount === 1 ? '' : 's') + '.' +
      '</div>' +
      '<div style="font-size:12px; line-height:1.7; max-width:540px; margin:0 auto;">' +
        'Carousel rendering is not yet wired in the modular build. ' +
        'View confirmed candidates on <b>page 2 candidate focus</b>.' +
      '</div>';
  }
  if (navBar) navBar.style.display = 'none';
  if (meta)   meta.style.display = 'none';
}

/**
 * Public entry: wire the prev/next/keydown handlers. Idempotent.
 *
 * TODO_MISSING(_wireConfirmedCarouselNav) — full handler wiring.
 */
export function initConfirmedCarousel() {
  // TODO_MISSING(_wireConfirmedCarouselNav)
  // Stub: no-op until the carousel renderer is implemented.
  return;
}

export const __MODULE_ID__ = 'inversion_catalogue/page9';
