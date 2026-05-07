// inversion_catalogue/page_overview.js
// =====================================================================
// Page "Overview" — synthesis-stage tab, declared but empty in legacy.
//
// Legacy state (Inversion_atlas.html line 9322):
//     <div id="page_overview" class="page"></div>
//
// The tab button exists at legacy line 5138 (data-page="page_overview"
// data-stage="synthesis") but no body and no render function were ever
// shipped. There are zero references to `renderOverview`, `renderPage_overview`,
// or `page_overview` in any JS scope of the legacy file — verified by
// `grep -niE "(renderOverview|page_overview|renderPageOverview)"`,
// which only returns the tab button and the empty <div>.
//
// This module exists so the page registry has a non-throwing entry for
// `page_overview`. If/when the synthesis overview gets designed, replace
// the body of renderPageOverview() with the real render logic and update
// page_overview.html with the panel skeleton.
// =====================================================================

export function wirePageOverview(state) {
  function renderPageOverview() {
    // TODO_MISSING(synthesis_overview_design) — legacy ships an empty stub.
    // The synthesis-stage overview was scoped but never implemented in
    // Inversion_atlas.html. Decide whether to:
    //   (a) drop the tab from the new build, or
    //   (b) populate it with a high-level workflow summary
    //       (counts of candidates per stage, layer-presence checklist, etc.)
    // No-op for now — keeps the empty <div> exactly as the legacy did.
  }

  return { renderPageOverview };
}

export default wirePageOverview;
