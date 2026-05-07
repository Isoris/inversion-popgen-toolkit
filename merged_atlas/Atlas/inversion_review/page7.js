// =============================================================================
// inversion_review/page7.js — "9 ancestry" tab
// =============================================================================
// Stage:        classification
// Legacy DOM:   <div id="page7"> (legacy lines 7660–7670)
// Renderer:     renderAncestryPage()  — referenced via typeof guard at legacy
//                                       lines 59629–59630, 59732–59733.
//                                       NOT defined inline in legacy HTML —
//                                       lives in an external script.
//
// What this page does
// -------------------
// Per-window ancestry view for the active chromosome. Reuses the popstats
// stack architecture (canvas tracks gated by chip toggles in #ancViewChips).
//
// Three view chips (per legacy line 8974 description):
//   - K-cluster label   (per-window argmax-Q assignment)
//   - Q-value heatmap   (top-1 ancestry component intensity)
//   - delta12           (top-1 minus top-2 Q — confidence)
//
// Driven by a `<chrom>_phase4_ancestry.json` layer that the user drag-drops
// into the page. Renders sample × window heatmaps for whichever chips are
// active. The empty-state in #ancNoChrom prompts the user to load the
// precomp first, then the ancestry layer.
//
// DOM contract:
//   #ancToolbar      — sticky chip bar (<span id="ancViewChips">) + meta
//   #ancStack        — vertical canvas-section stack
//   #ancNoChrom      — empty state with drag-drop instructions
//
// Extraction notes (Batch 2)
// --------------------------
// `renderAncestryPage` is called by the tab dispatcher but never defined
// inside Inversion_atlas.html. It must be in a sibling external script
// alongside the popstats wiring (page-6 and page-7 share the canvas-track
// architecture per the page-6 description: "Reuses the popstats stack
// architecture", legacy line 8974). The merge chat may wish to extract the
// canvas-track helpers shared between page6 and page7 into a sibling
// module under inversion_review/.
//
// Decision (BATCH_2_NOTES): same lifecycle-wrapper pattern as page6.js.
// =============================================================================


// -----------------------------------------------------------------------------
// External-file deps
// -----------------------------------------------------------------------------
//
// TODO_MISSING(renderAncestryPage)         — external (sibling of atlas_page6_wiring)
//                                              (dispatched at legacy 59629 + 59732)
//
// Layer schema:
// TODO_MISSING(ancestry layer loader)      — drag-drop loader for
//                                            <chrom>_phase4_ancestry.json.
//                                            Likely lives in main bootstrap.
// -----------------------------------------------------------------------------


/**
 * Show the ancestry page. The renderer is defined externally (sibling of
 * atlas_page6_wiring.js); this wrapper exists so the tab dispatcher has a
 * single ES-module entry to call.
 *
 * @param {object} state  shared state (not used directly here — the legacy
 *                        renderer reads window.state)
 * @returns {void}
 */
export function showAncestryPage(state) {
  if (typeof window === 'undefined') return;
  const fn = /** @type {any} */ (window).renderAncestryPage;
  if (typeof fn === 'function') {
    try {
      // TODO_MISSING(renderAncestryPage)
      fn();
    } catch (err) {
      console.warn('[page7] renderAncestryPage threw:', err);
    }
    return;
  }

  // Fallback empty state — surfaces the missing-dep clearly to the user.
  const stack = document.getElementById('ancStack');
  const noChrom = document.getElementById('ancNoChrom');
  if (noChrom) {
    noChrom.style.display = 'block';
    noChrom.innerHTML =
      'Ancestry renderer not loaded. The renderer is the sibling of ' +
      '<code>atlas_page6_wiring.js</code>; drop the ancestry script ' +
      'alongside this HTML and reload.';
  }
  if (stack) stack.innerHTML = '';
}


/**
 * Re-render the ancestry page after state changes (e.g. ancestry layer
 * loaded, chrom switch). Idempotent.
 *
 * @returns {void}
 */
export function refreshAncestryPage() {
  showAncestryPage();
}


// =============================================================================
// state references not in SLOT_REGISTRY
// =============================================================================
//
// state.data                 — IS in SLOT_REGISTRY (transient). The ancestry
//                              renderer reads state.data.ancestry (an
//                              optional sub-layer) for per-window K labels,
//                              Q-values, and delta12.
// state.ancestryViewChips    — Set of currently-active view chips (e.g.
//                              'label', 'qvalue', 'delta12'). Page-7
//                              private; not in SLOT_REGISTRY. Recommendation:
//                              keep page-private (UI state, transient).
// state.candidate            — IS in SLOT_REGISTRY. Used for the "selected
//                              candidate" cursor band on each ancestry section.
//
// No new SLOT_REGISTRY entries are introduced by this page.
// =============================================================================
