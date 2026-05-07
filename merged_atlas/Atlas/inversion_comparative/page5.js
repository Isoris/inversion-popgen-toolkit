// page5.js — Quick-reference help page
//
// IMPORTANT — handoff doc / actual content discrepancy:
//   HANDOFF_BATCH_5.md describes page5 as 'Multi-species comparison page (Pairwise
//   dotplots, synteny graph, breakpoint annotations)' at legacy lines 8159-9316.
//   In legacy/Inversion_atlas.html, lines 8159-9316 are actually the static
//   QUICK-REFERENCE / HELP page (data-page='page5', data-stage='help').
//   See:  legacy line 5142  →  <button data-page="page5" data-stage="help">
//                              <span class="num">16</span> help</button>
//   The actual multi-species cockpit is page16b (lines 8033-8100; renderer
//   _renderMultiSpeciesPage at legacy line 27497). page16b lives in page16b.js.
//
// Because page5 is purely declarative HTML (the help/vocabulary/hotkeys/pipeline
// reference content), it has NO dedicated render function in the legacy file.
// CSS for #page5 lives at legacy lines 2054-2071. The only JS reference to page5
// is the tabBar router at legacy line 5142 and a single dispatch at line 27509
// (which is part of _renderMultiSpeciesPage's internal tabBar refresh logic).
//
// This module is therefore a stub. It exists so the page-loader can dispatch to
// it the same way it dispatches to the other pages (uniform routing).

export function renderPage5(/* state */) {
  // Help page is static HTML — nothing to do at render time. The HTML content
  // lives in page5.html and is mounted by the shell when the help tab activates.
  return;
}

export const PAGE5_META = {
  id: 'page5',
  stage: 'help',
  label: 'help',
  num: 16,  // tab number from legacy
  static: true,
};
