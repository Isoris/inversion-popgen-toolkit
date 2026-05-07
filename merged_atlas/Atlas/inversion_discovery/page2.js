// Atlas/inversion_discovery/page2.js
//
// Multi-panel deep-dive on a single promoted candidate.
// renderCandidateMetadata composes ~15 sub-panels (header, σ profile, bands,
//   notes, etc) by calling helpers like candidateNavHtml(), candidateBandsHtml(),
//   etc — all flagged TODO_MISSING below.
//
// Extracted LITERALLY from legacy/Inversion_atlas.html.
// Per HANDOFF_BATCH_1: state is now passed as first arg (was a global in legacy).
// All other unresolved references are flagged TODO_MISSING — the merge chat
// will resolve them once all 5 batches land.
//
// Entry points (in extraction order):
  // renderCandidateMetadata() — legacy lines 58421-58481
  // wireCandidateNav() — legacy lines 59491-59595
//

import { contextFromState, clusterL2, ClusterCache } from '../shared/per_l2_cluster.js';
import { hetRateColor } from '../shared/het_rate.js';
import { alignLabels, hungarianChainProjection, concordanceMatrix } from '../shared/hungarian.js';
import { buildContingency, computeARI, computeNMI, cramersV } from '../shared/contingency.js';
import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

// ---------------------------------------------------------------------------
// TODO_MISSING — functions referenced but not defined in this module or shared/*.
// (These are called by the extracted bodies below. The merge chat will either
//  pull them from another batch, or hoist them into shared/ if they are reused.)
// ---------------------------------------------------------------------------
// TODO_MISSING(_navigateToCandidate)
// TODO_MISSING(_wireCandidateBandClicks)
// TODO_MISSING(_wireCandidateBlockChips)
// TODO_MISSING(_wireCandidateDosageHeatmap)
// TODO_MISSING(_wireCandidateHaplotypeAnnotations)
// TODO_MISSING(_wireCandidateRegimeRow)
// TODO_MISSING(addCandidateToList)
// TODO_MISSING(candidateAgeOriginHtml)
// TODO_MISSING(candidateAncestryConfoundHtml)
// TODO_MISSING(candidateBandComposition)
// TODO_MISSING(candidateBandsHtml)
// TODO_MISSING(candidateBlockChipsHtml)
// TODO_MISSING(candidateDosageHeatmapHtml)
// TODO_MISSING(candidateFromJSON)
// TODO_MISSING(candidateHaplotypeAnnotationsHtml)
// TODO_MISSING(candidateHeaderHtml)
// TODO_MISSING(candidateHetShapeHtml)
// TODO_MISSING(candidateListClosestIndex)
// TODO_MISSING(candidateListIndexOf)
// TODO_MISSING(candidateListSortedByPos)
// TODO_MISSING(candidateNavHtml)
// TODO_MISSING(candidateNotesHtml)
// TODO_MISSING(candidateProfileHtml)
// TODO_MISSING(candidateRegimeRowHtml)
// TODO_MISSING(candidateRichCardHtml)
// TODO_MISSING(candidateSigmaChartHtml)
// TODO_MISSING(candidateSubbandHtml)
// TODO_MISSING(candidateSummaryHtml)
// TODO_MISSING(candidateToJSON)
// TODO_MISSING(drawCandGHSLPerBand)
// TODO_MISSING(drawCandLinesPanel)
// TODO_MISSING(drawCandLocalPCA)
// TODO_MISSING(drawCandidateLocationStrip)
// TODO_MISSING(drawCandidateSigmaChart)
// TODO_MISSING(persistCandidateList)
// TODO_MISSING(refreshCandidateListUI)
// TODO_MISSING(refreshCandidateUI)
// TODO_MISSING(renderCatalogue)
// TODO_MISSING(sigmaProfileCandidate)
// TODO_MISSING(wireCandidateAncestryConfound)
// TODO_MISSING(wireCandidateButtons)

// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT — state.<slot> references that are NOT in
// shared/state.js SLOT_REGISTRY. Some may be ad-hoc geometry caches the
// legacy added imperatively (state._simGeom, state._zGeom, etc.); some may
// be slot names we need to register. Merge chat decides per-slot.
// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT(state.candidatePageMode)

// ---------------------------------------------------------------------------
// Extracted bodies
// ---------------------------------------------------------------------------

// --- renderCandidateMetadata() — legacy lines 58421-58481 ---
export function renderCandidateMetadata(state) {
  // Full deep-dive view (Turn B). Composes:
  //   - candidate header (chrom, span, source, K, action buttons)
  //   - σ across candidate verdict + chart + drifters
  //   - per-band composition cards
  //   - notes textarea
  const slot = document.getElementById('candidateMeta');
  const empty = document.getElementById('candidateEmpty');
  if (!slot) return;
  if (!state.candidate) {
    slot.innerHTML = '';
    slot.style.display = 'none';
    if (empty) empty.style.display = 'block';
    return;
  }
  const c = state.candidate;
  if (empty) empty.style.display = 'none';
  slot.style.display = 'block';

  const profile = sigmaProfileCandidate(c);
  const bands = candidateBandComposition(c);

  slot.innerHTML =
    candidateNavHtml(c) +                  // v3.47: prev / next + confirmed toggle
    candidateHeaderHtml(c) +
    candidateBlockChipsHtml(c) +           // v3.99 t14e+ continue: 18-block chip row + inspector mount (SCHEMA §21)
    candidateRichCardHtml(c) +             // v3.47: location strip + multi-panel grid
    candidateSummaryHtml(c, profile, bands) +
    candidateSubbandHtml(c) +              // v3.82: K=6 nesting verdict + per-group purity
    candidateHetShapeHtml(c) +             // v3.93: FIG_C07-style ridgeline (placement B)
    candidateDosageHeatmapHtml(c) +        // v3.94: FIG_C08-style dosage heatmap (static view)
    candidateProfileHtml(c, profile) +
    candidateSigmaChartHtml(c, profile) +
    candidateBandsHtml(c, bands) +
    candidateHaplotypeAnnotationsHtml(c) +     // turn 118: per-band haplotype labels
    candidateAncestryConfoundHtml(c) +     // v4 turn 9: ancestry-confound panel
    candidateRegimeRowHtml(c) +            // v4 turn 15b: regime registry strip
    candidateAgeOriginHtml(c) +            // turn 119: cheat30 age & origin panel
    candidateNotesHtml(c);

  // Wire interactions after DOM insertion
  wireCandidateButtons(c, profile);
  wireCandidateNav(c);                     // v3.47
  if (typeof _wireCandidateBlockChips === 'function') {
    try { _wireCandidateBlockChips(); } catch (e) { console.warn('_wireCandidateBlockChips:', e.message); }
  }
  try { wireCandidateAncestryConfound(c); } catch (e) { console.warn('wireCandidateAncestryConfound:', e.message); }
  try { _wireCandidateHaplotypeAnnotations(c); } catch (e) { console.warn('haplotype annotations:', e.message); }
  try { _wireCandidateBandClicks(c, bands); } catch (e) { console.warn('band clicks:', e.message); }
  try { _wireCandidateRegimeRow(c); } catch (e) { console.warn('_wireCandidateRegimeRow:', e.message); }
  drawCandidateLocationStrip(c);           // v3.47: simdat / L1 / karyogram minis
  // v3.48: ready-data analysis panels
  requestAnimationFrame(() => {
    try { drawCandLocalPCA(c); } catch (e) { console.warn('drawCandLocalPCA:', e.message); }
    try { drawCandLinesPanel(c); } catch (e) { console.warn('drawCandLinesPanel:', e.message); }
    try { drawCandGHSLPerBand(c); } catch (e) { console.warn('drawCandGHSLPerBand:', e.message); }
    // v3.94: dosage heatmap (static view) — wired after DOM mount
    try { _wireCandidateDosageHeatmap(c); } catch (e) { console.warn('dosage heatmap:', e.message); }
  });
  if (profile && profile.sd) drawCandidateSigmaChart(c, profile);
}

// --- wireCandidateNav() — legacy lines 59491-59595 ---
export function wireCandidateNav(state, c) {
  const prevBtn = document.getElementById('candNavPrev');
  const nextBtn = document.getElementById('candNavNext');
  const confirmBtn = document.getElementById('candConfirmedToggle');
  const sorted = candidateListSortedByPos();
  if (prevBtn) {
    prevBtn.addEventListener('click', () => {
      if (sorted.length === 0) return;
      let idx = candidateListIndexOf(c);
      if (idx < 0) {
        // Active candidate not in list; jump to nearest neighbor on the LEFT
        // of cur position, falling back to closest if no left neighbor.
        const nearest = candidateListClosestIndex(c);
        idx = (nearest > 0 && sorted[nearest].start_bp > c.start_bp) ? nearest - 1 : nearest;
      } else {
        idx = Math.max(0, idx - 1);
      }
      const target = sorted[idx];
      if (target) _navigateToCandidate(target);
    });
  }
  if (nextBtn) {
    nextBtn.addEventListener('click', () => {
      if (sorted.length === 0) return;
      let idx = candidateListIndexOf(c);
      if (idx < 0) {
        const nearest = candidateListClosestIndex(c);
        idx = (nearest < sorted.length - 1 && sorted[nearest].start_bp < c.start_bp)
              ? nearest + 1 : nearest;
      } else {
        idx = Math.min(sorted.length - 1, idx + 1);
      }
      const target = sorted[idx];
      if (target) _navigateToCandidate(target);
    });
  }
  if (confirmBtn) {
    confirmBtn.addEventListener('click', () => {
      // Mutate active candidate's confirmed flag
      c.confirmed = !c.confirmed;
      // v3.99 turn 12 fix: if the candidate isn't in the saved list yet (e.g.
      // it was opened via catalogue "view as candidate" without explicit save),
      // auto-add it now when the user confirms. Without this, "mark confirmed"
      // had no effect on the saved list — the page-9 confirmed walkthrough
      // would never see the candidate. We only auto-add on confirm=true; on
      // un-confirm of an unsaved candidate we keep the list untouched (the
      // user clearly never wanted this in the list).
      let inList = state.candidateList.find(x => x.id === c.id);
      if (!inList && c.confirmed) {
        // Snapshot mirrors the page-2 list-toggle path
        const snapshot = candidateFromJSON(candidateToJSON(c));
        snapshot.confirmed = true;
        addCandidateToList(snapshot);
        inList = state.candidateList.find(x => x.id === c.id);
      }
      // Mirror the flag onto the saved copy
      if (inList) inList.confirmed = c.confirmed;
      persistCandidateList();
      // v3.99 turn 12: refresh the catalogue so the new ✓ green tint + ID
      // chip on confirmed rows reflects the change immediately. Cheap re-render.
      if (typeof renderCatalogue === 'function') {
        try { renderCatalogue(); } catch (_) {}
      }
      // v3.56: if we're in confirmed-only mode (page 8) and the user just
      // UN-confirmed the current candidate, it disappears from the navigation.
      // Advance to the next confirmed candidate, or show the empty state if
      // none remain.
      if (state.candidatePageMode === 'confirmed' && !c.confirmed) {
        const remaining = state.candidateList.filter(x => x && x.confirmed);
        if (remaining.length > 0) {
          // Pick the closest by position to the candidate we just unconfirmed
          const cMid = (c.start_bp + c.end_bp) / 2;
          let bestI = 0, bestD = Infinity;
          for (let i = 0; i < remaining.length; i++) {
            const mid = (remaining[i].start_bp + remaining[i].end_bp) / 2;
            const d = Math.abs(mid - cMid);
            if (d < bestD) { bestD = d; bestI = i; }
          }
          _navigateToCandidate(remaining[bestI]);
          refreshCandidateListUI();
          return;
        } else {
          // No more confirmed candidates — show the page-8 empty state
          state.candidate = null;
          const slot = document.getElementById('candidateMeta');
          const empty = document.getElementById('candidateEmpty');
          if (slot) { slot.innerHTML = ''; slot.style.display = 'none'; }
          if (empty) {
            empty.style.display = 'block';
            empty.innerHTML =
              '<div style="font-size:14px; margin-bottom:14px;">' +
              'No more confirmed candidates.</div>' +
              '<div style="font-size:12px; line-height:1.7; max-width:540px; margin:0 auto;">' +
              'Confirm candidates from <b>page 2 candidate focus</b> to populate this view.' +
              '</div>';
          }
          refreshCandidateListUI();
          return;
        }
      }
      refreshCandidateUI();
      refreshCandidateListUI();
    });
  }
}
