// =============================================================================
// inversion_review/page6.js — "8 popstats" tab
// =============================================================================
// Stage:        classification
// Legacy DOM:   <div id="page6"> (legacy lines 7647–7658)
// Renderer:     renderPopstatsPage()  — defined in external js/atlas_page6_wiring.js
// Tab dispatch: legacy lines 59626–59627, 59729–59730
//
// What this page does
// -------------------
// Stack of population-genetic tracks aligned to the chromosome:
//   - |Z| / score
//   - SNP density
//   - BEAGLE imputation uncertainty
//   - depth / coverage
//   - θπ (Tajima's per-window pi)
//   - F_ST (between karyotype groups)
//   - Hobs / Hexp (observed vs expected heterozygosity)
//   - ancestry Δ12 (top-1 minus top-2 Q)
//
// The track inventory lives in the request layer (atlas_request_layer.js)
// which talks to a popstats live server (POST /api/popstats/*); the page-6
// wiring (atlas_page6_wiring.js) wraps that with track-def slot-combine UI
// and tooltips. Click chips in #psChips to toggle visibility.
//
// DOM contract:
//   #psToolbar    — sticky chip bar (<span id="psChips">)
//   #psStack      — vertical canvas-track stack
//   #psNoChrom    — empty state when no precomp loaded
//   #psGalleryTray — collapsible track-discovery sidebar (turn 6.5)
//
// Extraction notes (Batch 2)
// --------------------------
// The popstats stack is driven entirely by external JS files (per the script
// inventory at legacy line 54963–54972):
//   - js/atlas_request_layer.js  → window.popgenLive
//   - js/atlas_page6_wiring.js   → window.popgenPage6
//   - js/atlas_track_gallery.js  → window.popgenGallery
//
// Inside Inversion_atlas.html itself there is NO `function renderPopstatsPage`
// — the legacy code only references it via `if (typeof renderPopstatsPage
// === 'function')` guards. So the page module here is a thin lifecycle
// wrapper, like page_sv_evidence.js.
//
// Decision (BATCH_2_NOTES): kept the popstats stack as an external-script
// dep rather than promoting to shared. The whole thing depends on a
// popstats live server and IndexedDB caching that doesn't fit the
// "pure-helper" shape of the shared/ modules.
// =============================================================================


// -----------------------------------------------------------------------------
// External-file deps
// -----------------------------------------------------------------------------
//
// TODO_MISSING(renderPopstatsPage)        — js/atlas_page6_wiring.js
//                                            (dispatched at legacy 59626 + 59729)
// TODO_MISSING(popgenLive)                — js/atlas_request_layer.js
//                                            (POST /api/popstats/* wrappers)
// TODO_MISSING(popgenPage6)               — js/atlas_page6_wiring.js
// TODO_MISSING(popgenGallery)             — js/atlas_track_gallery.js
// -----------------------------------------------------------------------------


/**
 * Show the popstats page. The renderer is defined in atlas_page6_wiring.js;
 * this wrapper exists so the tab dispatcher has a single ES-module entry
 * to call.
 *
 * @param {object} state  shared state (not used directly here — the legacy
 *                        renderer reads window.state)
 * @returns {void}
 */
export function showPopstatsPage(state) {
  if (typeof window === 'undefined') return;
  const fn = /** @type {any} */ (window).renderPopstatsPage;
  if (typeof fn === 'function') {
    try {
      // TODO_MISSING(renderPopstatsPage)
      fn();
    } catch (err) {
      console.warn('[page6] renderPopstatsPage threw:', err);
    }
    return;
  }

  // Fallback empty state — surfaces the missing-dep clearly to the user.
  const stack = document.getElementById('psStack');
  const noChrom = document.getElementById('psNoChrom');
  if (noChrom) {
    noChrom.style.display = 'block';
    noChrom.innerHTML =
      'Popstats wiring (<code>atlas_page6_wiring.js</code>) not loaded. ' +
      'Drop the page-6 script bundle alongside this HTML and reload.';
  }
  if (stack) stack.innerHTML = '';
}


/**
 * Re-render the popstats page after state changes (e.g. chrom switch,
 * candidate change). Idempotent — safe to call repeatedly.
 *
 * @returns {void}
 */
export function refreshPopstatsPage() {
  // Same dispatch as show — the legacy renderer is itself idempotent.
  showPopstatsPage();
}


// =============================================================================
// state references not in SLOT_REGISTRY
// =============================================================================
//
// state.data                 — IS in SLOT_REGISTRY (transient). The popstats
//                              renderer reads its layers (theta_pi, fst,
//                              hobs/hexp, etc.) from inside state.data.
// state.candidate            — IS in SLOT_REGISTRY (cross_atlas).
// state.popstatsLive         — page-6 cache, NOT in SLOT_REGISTRY. Holds
//                              IndexedDB-backed responses keyed by
//                              { chrom, group_set_id, metric }. Owned by
//                              atlas_request_layer.js. Recommendation:
//                              keep page-private; do NOT promote.
// state.popstatsTracksOn     — Set of currently-active chip IDs. Page-6
//                              private. Recommendation: keep page-private.
// state.popstatsGalleryOpen  — gallery tray collapsed/expanded flag. Page-6
//                              private.
//
// All three popstats-* state slots are page-private and intentionally
// stay outside SLOT_REGISTRY (they're caches/UI flags, not data).
// =============================================================================
