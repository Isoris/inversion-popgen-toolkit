// Atlas/inversion_catalogue/page3.js
// =============================================================================
// page3 — Catalogue (sortable, filterable table of all L2 envelopes)
// (`<div id="page3">` contains the toolbar and `<table id="catTable">`)
//
// Source: legacy/Inversion_atlas.html lines 7261–7368 (HTML shell)
//
// IMPORTANT — legacy state of page3 rendering:
//
//   The `renderCatalogue` function is REFERENCED 9 times in the legacy file
//   (lines 13750, 13963, 13974, 14055, 14084, 55642, 55853, 59551, 64293)
//   but is NEVER DEFINED. Every call site uses the
//     if (typeof renderCatalogue === 'function') renderCatalogue();
//   guard, meaning the legacy build expects an EXTERNAL renderer to install
//   `window.renderCatalogue` at runtime. None of the legacy JS in
//   Inversion_atlas.html does so. Confirmation:
//
//     $ grep -nE 'function renderCatalogue|var renderCatalogue|let renderCatalogue|const renderCatalogue|renderCatalogue\s*=|window\.renderCatalogue' legacy/Inversion_atlas.html
//     (only the call-site `typeof` guards above match)
//
//   The catalogue table (#catTable, #catHead, #catBody) and toolbar
//   (#catFilter, #catVerdictFilter, #catSelectAll, #catClearSel,
//   #catExportTSV, #catExportMD, #catExportJSON, the gallery + breeding
//   exports, view-mode buttons, display-mode buttons, diamond-mode
//   buttons) all live in the HTML at lines 7261–7368, but the JavaScript
//   that populates them is missing from this legacy drop.
//
//   The ONE piece of catalogue-toolbar JS that DOES live in legacy is the
//   breeding-export wiring at lines 23668–23715
//   (`_wireCatalogueBreedingExportBtns`). It wires #catBreedingTierSel,
//   #catExportBreedingHTML, #catExportBreedingJSON to dispatchers in the
//   Turn 146 breeding-card section (legacy ~23236–23715). The dispatchers
//   themselves (`_exportBreedingCardsHTML`, `_exportBreedingCardsJSON`,
//   `_buildBreedingCardsCombinedHTML`, `_buildBreedingCardsJSONBundle`,
//   `_filterCandsForBreedingExport`, `_BREEDING_EXPORT_TIER_MODES`) are
//   all in that ~480-line block — too large to extract as part of page3
//   without overstepping this batch's scope. Marked TODO_MISSING.
//
// Decision for this batch: ship page3 as a SHELL — public entry
// `renderCataloguePage()` does the minimum useful thing (shows the
// "Load a JSON to populate the catalogue." empty-state message), and
// all real catalogue logic is left as TODO_MISSING for the merge chat.
//
// External dependencies (TODO_MISSING):
//
//   TODO_MISSING(_buildCatalogueRows)
//     — the row-builder that turns state.l2Envelopes / state.candidateList
//       / state.regimeRegistry into table rows. Does not exist in legacy.
//   TODO_MISSING(_filterCatalogueRows)
//     — applies #catFilter, #catVerdictFilter, view-mode (L2 / L1 /
//       favorites / L3), simple/detailed display mode. Does not exist.
//   TODO_MISSING(_sortCatalogueRows)
//     — header-click sort. Does not exist.
//   TODO_MISSING(_paintCatalogueRow)
//     — per-row HTML builder including the verdict pill, diamond pill,
//       cand-confirmed CSS class, etc. Does not exist.
//   TODO_MISSING(_exportCatalogueTSV)
//     — wired to #catExportTSV. Does not exist.
//   TODO_MISSING(_exportCatalogueMarkdown)
//     — wired to #catExportMD. Does not exist.
//   TODO_MISSING(_exportCatalogueJSON)
//     — wired to #catExportJSON. Does not exist.
//   TODO_MISSING(_exportCatalogueGallerySVG / PNG / PDF)
//     — wired to the 📷 gallery export buttons. Does not exist.
//   TODO_MISSING(_wireCatalogueBreedingExportBtns)
//     — DOES exist in legacy at lines 23668–23715. Depends on the Turn
//       146 export dispatchers in legacy lines ~23236–23715.
//   TODO_MISSING(_exportBreedingCardsHTML)
//     — legacy ~23260+ (Turn 146 dispatcher).
//   TODO_MISSING(_exportBreedingCardsJSON)
//     — legacy ~23260+ (Turn 146 dispatcher).
//   TODO_MISSING(_filterCandsForBreedingExport)
//     — legacy Turn 146.
//   TODO_MISSING(_buildBreedingCardsCombinedHTML)
//     — legacy Turn 146.
//   TODO_MISSING(_buildBreedingCardsJSONBundle)
//     — legacy Turn 146.
//   TODO_MISSING(_BREEDING_EXPORT_TIER_MODES)
//     — legacy Turn 146 constant table.
//   TODO_MISSING(_openRegimeRegistryDialog)
//     — wired to #catRegimeRegistry.
//   TODO_MISSING(_assignSelectedToRegime)
//     — wired to #catRegimeAssignSel.
//   TODO_MISSING(_promoteSelectedAsCandidate)
//     — wired to #catViewAsCandidate.
//   TODO_MISSING(makeShelfLDPanel)
//     — populates #page3_q09b_slot. Comment in the HTML says
//       "atlas_turn7.js's makeShelfLDPanel()" — separate JS file.
//   TODO_MISSING(makeLDSplitPanel)
//     — populates #page3_ld_slot. Comment in the HTML says
//       "atlas_ld.js's makeLDSplitPanel()" — separate JS file.
//
//   global `state`               — when the renderer ships, will read:
//                                    state.l2Envelopes
//                                    state.candidateList
//                                    state.regimeRegistry
//                                    state.catalogueViewMode
//                                    state.catalogueDisplayMode
//                                    state.diamondMode
//                                    state.catalogueSort
//                                    state.catalogueSelected
//                                    state.catalogueFavorites
//
// =============================================================================

const state = (typeof window !== 'undefined' && window.state) ? window.state : {};

/**
 * Public entry: render the catalogue page.
 *
 * Current implementation (matches legacy reality where `renderCatalogue` is
 * not defined): show the "Load a JSON to populate the catalogue." empty
 * message in #catEmpty, leave the table empty.
 *
 * The merge chat (or a follow-up batch) implements:
 *   _buildCatalogueRows -> _filterCatalogueRows -> _sortCatalogueRows
 *   -> per-row paint -> wire toolbar handlers.
 */
export function renderCataloguePage() {
  if (typeof document === 'undefined') return;
  const head  = document.getElementById('catHead');
  const body  = document.getElementById('catBody');
  const empty = document.getElementById('catEmpty');
  const selInfo = document.getElementById('catSelInfo');

  // TODO_MISSING(_buildCatalogueRows) — fill #catHead and #catBody.
  if (head) head.innerHTML = '';
  if (body) body.innerHTML = '';
  if (empty) {
    empty.style.display = 'block';
    empty.textContent = 'Load a JSON to populate the catalogue.';
  }
  if (selInfo) selInfo.textContent = '0 selected of 0';
}

/**
 * Public entry: idempotent toolbar wiring.
 *
 * Currently a no-op. The merge chat wires:
 *   - view-mode buttons (#catViewFav, #catViewL3, #catViewL2, #catViewL1)
 *   - display-mode buttons (#catDispSimple, #catDispDetailed)
 *   - diamond-mode buttons (#catDiamondLoose, #catDiamondStrict, #catDiamondStrict2)
 *   - filter inputs (#catFilter, #catVerdictFilter)
 *   - selection buttons (#catSelectAll, #catClearSel)
 *   - regime buttons (#catRegimeRegistry, #catRegimeAssignSel)
 *   - candidate-promotion (#catViewAsCandidate)
 *   - export buttons (#catExportTSV, #catExportMD, #catExportJSON,
 *                     #catExportGallerySVG/PNG/PDF,
 *                     #catBreedingTierSel, #catExportBreedingHTML/JSON)
 */
export function initCataloguePage() {
  // TODO_MISSING(_wireCatalogueToolbar)
  // TODO_MISSING(_wireCatalogueBreedingExportBtns) — legacy lines 23668–23715
  return;
}

export const __MODULE_ID__ = 'inversion_catalogue/page3';
