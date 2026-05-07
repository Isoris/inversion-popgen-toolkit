// =============================================================================
// inversion_review/page4.js — "7 karyotype / tier" tab
// =============================================================================
// Stage:        refinement
// Legacy DOM:   <div id="page4"> (legacy lines 7573–7645)
// Renderers:    renderCandidateKaryotype()       (legacy lines 62800–62836)
//               _renderCandidateKaryotypeBody()  (legacy lines 63091–63289)
//               renderCandidateTier()            (legacy lines 63411–63492)
//               _refreshSubviewButtonStyles()    (legacy lines 62840–62861)
//               _wireCandSubviewButtons()        (legacy line 63547)
//               karyoState (page-local UI state) (legacy lines 62508–62526)
//
// What this page does
// -------------------
// Two sub-views of the active candidate, toggleable with the buttons in
// #candKaryoSubviewBar:
//
//   (1) Karyotype — per-sample regime breakdown. Each row is one fish
//       showing its locked K=3 label (HOMO_1 / HET / HOMO_2 in the
//       operational H-system, ordered by median PC1). Two-track candidates
//       get the extra track-membership pills (active_bands per track).
//       Sortable, filterable, band-filterable.
//
//   (2) Tier — 14-axis classification view. Renders an empty-state card
//       until the cluster-side R-pipeline ships `final_classification.json`
//       (cross-references SCHEMA_V2.md §19; sources mentioned: §13 evidence
//       framework, §9 cluster-emit `classification` layer, candidate
//       schema §13 completion + characterization blocks).
//
// Left side of the page is the candidate list pane (#candListPane) with
// import/export/registry/manuscript-bundle/clear actions. The candidate-list
// rendering itself is `refreshCandidateListUI()` (legacy line 62528) — used
// by other tabs too, so it belongs in the merge chat's promotion list.
//
// Extraction notes (Batch 2)
// --------------------------
// Three functions are extracted verbatim:
//   - renderCandidateKaryotype  (dispatcher: routes to body or tier)
//   - renderCandidateTier       (14-axis grid renderer + empty state)
//   - _refreshSubviewButtonStyles (idempotent button-style updater)
//
// `_renderCandidateKaryotypeBody` is ~200 lines of table-building logic
// that touches ~15 helpers (sigmaProfileCandidate, _isKaryoTwoTrack,
// getKaryotypeLabel, getKaryotypeLabelCaveat, groupColor, buildKaryotypeRows,
// ...). It's left as TODO_MISSING — it's a plain copy-port and the merge
// chat can pull it in en bloc.
// =============================================================================

import { SLOT_REGISTRY } from '../shared/state.js';

// -----------------------------------------------------------------------------
// Page-local UI state (not in SLOT_REGISTRY — page-private)
// Legacy: lines 62508–62526
// -----------------------------------------------------------------------------
const _CAND_SUBVIEW_KEY = 'pca_scrubber_v3.candSubview';

export const karyoState = {
  sortKey: 'k_label',     // 'cga' | 'k_label' | 'sigma' | 'family_id'
  sortAsc: true,
  filter: '',
  bandFilter: '',         // '' = all bands; 'k0', 'k1', ...
  // Page-4 sub-view toggle. 'karyotype' (default) = existing per-sample
  // regime table. 'tier' = 14-axis classification view (cluster-side
  // R-pipeline output, planned). When set to 'tier' and no
  // final_classification layer is loaded, renders an empty state.
  subview: 'karyotype',
};

// Restore persisted subview choice
try {
  const saved = (typeof localStorage !== 'undefined') ? localStorage.getItem(_CAND_SUBVIEW_KEY) : null;
  if (saved === 'karyotype' || saved === 'tier') karyoState.subview = saved;
} catch (_) { /* no-op */ }


// -----------------------------------------------------------------------------
// Unresolved external deps — to be filled by the merge chat.
// -----------------------------------------------------------------------------
//
// Karyotype body + supporting helpers:
// TODO_MISSING(_renderCandidateKaryotypeBody)        — legacy line 63091
// TODO_MISSING(buildKaryotypeRows)                   — legacy line 62703
// TODO_MISSING(_isKaryoTwoTrack)                     — likely sigma helpers
// TODO_MISSING(sigmaProfileCandidate)                — sigma-profile helper
// TODO_MISSING(getKaryotypeLabel)                    — legacy line 37022
// TODO_MISSING(getKaryotypeLabelCaveat)              — legacy line 37040
// TODO_MISSING(groupColor)                           — palette helper
// TODO_MISSING(refreshCandidateListUI)               — legacy line 62528
//                                                       (used by every page
//                                                        with a cand-list pane)
//
// Tier-specific:
// TODO_MISSING(_renderTierAxesGrid)                  — legacy ~63495 (right
//                                                       after renderCandidateTier)
// -----------------------------------------------------------------------------


// =============================================================================
// renderCandidateKaryotype — public entry, dispatcher
// =============================================================================
// Legacy: lines 62800–62836. Verbatim except for the import-time karyoState.
// =============================================================================
export function renderCandidateKaryotype() {
  const empty = document.getElementById('candKaryoEmpty');
  const content = document.getElementById('candKaryoContent');
  const tierContent = document.getElementById('candTierContent');
  const subviewBar = document.getElementById('candKaryoSubviewBar');
  if (!empty || !content) return;
  if (!state.candidate) {
    empty.style.display = 'block';
    content.style.display = 'none';
    if (tierContent) { tierContent.style.display = 'none'; tierContent.innerHTML = ''; }
    if (subviewBar) subviewBar.style.display = 'none';
    content.innerHTML = '';
    return;
  }
  empty.style.display = 'none';
  if (subviewBar) subviewBar.style.display = 'flex';
  // Ensure the subview-toggle buttons are wired (idempotent guard inside)
  if (typeof _wireCandSubviewButtons === 'function') {
    try { _wireCandSubviewButtons(); } catch (_) { /* swallow */ }
  }

  // Subview routing. Show karyotype or tier panel based on karyoState.subview.
  const subview = (karyoState.subview === 'tier') ? 'tier' : 'karyotype';
  if (subview === 'tier') {
    content.style.display = 'none';
    if (tierContent) tierContent.style.display = 'block';
    _refreshSubviewButtonStyles(subview);
    renderCandidateTier();
    return;
  }
  // Karyotype (default): existing render path
  content.style.display = 'block';
  if (tierContent) tierContent.style.display = 'none';
  _refreshSubviewButtonStyles(subview);
  _renderCandidateKaryotypeBody();
}


// =============================================================================
// _refreshSubviewButtonStyles — idempotent button-state updater
// =============================================================================
// Legacy: lines 62840–62861. Verbatim.
// =============================================================================
export function _refreshSubviewButtonStyles(active) {
  const kBtn = document.getElementById('candSubviewKaryoBtn');
  const tBtn = document.getElementById('candSubviewTierBtn');
  const setActive = (btn, isActive) => {
    if (!btn || !btn.style) return;
    if (isActive) {
      btn.classList && btn.classList.add('active');
      btn.style.background = 'rgba(245,165,36,0.15)';
      btn.style.borderColor = 'var(--accent, #f5a524)';
      btn.style.color = 'var(--ink)';
      btn.style.fontWeight = '600';
    } else {
      btn.classList && btn.classList.remove('active');
      btn.style.background = 'var(--panel-2)';
      btn.style.borderColor = 'var(--rule)';
      btn.style.color = 'var(--ink-dim)';
      btn.style.fontWeight = '400';
    }
  };
  setActive(kBtn, active === 'karyotype');
  setActive(tBtn, active === 'tier');
}


// =============================================================================
// renderCandidateTier — 14-axis classification view + empty state
// =============================================================================
// Legacy: lines 63411–63492. Verbatim.
//
// Reads:    state.candidate, state.data.final_classification[<cid>]
// Renders:  empty-state card if final_classification not loaded, else the
//           14-axis grid via _renderTierAxesGrid (TODO_MISSING).
// =============================================================================
export function renderCandidateTier() {
  const tierContent = document.getElementById('candTierContent');
  const layerStatus = document.getElementById('candTierLayerStatus');
  if (!tierContent) return;
  if (!state.candidate) {
    tierContent.innerHTML = '';
    return;
  }
  const c = state.candidate;
  const layer = (state.data && state.data.final_classification) || null;
  const layerLoaded = !!layer && typeof layer === 'object';
  // Each candidate's classification keyed by candidate id (preferred) or
  // fallback to ref_l2. The R-side spec uses `candidate_id` as the key.
  const candKey = c.id || (c.ref_l2 != null ? `ref_l2_${c.ref_l2}` : null);
  const candEntry = (layerLoaded && candKey) ? (layer[candKey] || layer[c.id] || null) : null;

  // Update the small "layer status" text in the subview bar
  if (layerStatus) {
    if (layerLoaded && candEntry) {
      layerStatus.textContent = `final_classification: loaded · ${Object.keys(candEntry).length} axes for this candidate`;
      layerStatus.style.color = 'var(--good, #3cc08a)';
    } else if (layerLoaded) {
      layerStatus.textContent = 'final_classification: loaded · no entry for this candidate';
      layerStatus.style.color = 'var(--ink-dim)';
    } else {
      layerStatus.textContent = 'final_classification: NOT LOADED — render empty';
      layerStatus.style.color = 'var(--ink-dim)';
    }
  }

  if (!layerLoaded || !candEntry) {
    // EMPTY STATE — explain the contract
    tierContent.innerHTML = `
      <div style="padding: 32px 28px; max-width: 880px; margin: 0 auto;">
        <h3 style="margin: 0 0 6px; font-family: var(--serif); font-weight: 500;">
          14-axis classification for this candidate
        </h3>
        <div style="color: var(--ink-dim); font-size: 12px; margin-bottom: 16px;">
          candidate ${c.id ? c.id.replace(/^cand_/, '') : '?'} · ${c.chrom} ${(c.start_bp/1e6).toFixed(2)}–${(c.end_bp/1e6).toFixed(2)} Mb
        </div>
        <div style="background: var(--panel); border: 1px solid var(--rule); border-radius: 6px;
                    padding: 16px 20px; font-size: 12px; line-height: 1.6; color: var(--ink-dim);">
          <div style="font-weight: 600; color: var(--ink); margin-bottom: 8px;">
            ⚠ The <code>final_classification</code> JSON layer isn't loaded yet.
          </div>
          <p style="margin: 0 0 8px;">The Tier view consumes a per-candidate classification produced by the cluster-side R-pipeline (<code>characterize_candidate.R</code> + <code>classify_inversions.R</code>, planned). When that pipeline emits its output, the scrubber reads it as <code>state.data.final_classification</code> (an optional top-level JSON layer) and renders the 14 axes here.</p>
          <p style="margin: 0 0 8px;">Until then, this view shows the axis schema only — <b>no values</b>. See SCHEMA_V2.md §19 for the full contract.</p>
          <p style="margin: 0;">In the meantime, partial information lives in <code>state.data.classification</code> (§9 cluster-emit layer) and the candidate's own <code>completion</code> + <code>characterization</code> blocks (schema 2.12 §13). The Tier view will eventually unify all three sources into one read-only display.</p>
        </div>

        <h4 style="margin: 22px 0 8px; font-size: 12px; text-transform: uppercase;
                   letter-spacing: 0.08em; color: var(--ink-dim);
                   font-family: var(--mono); font-weight: 600;">
          14-axis schema (read-only preview)
        </h4>
        ${_renderTierAxesGrid(null)}
      </div>
    `;
    return;
  }

  // VALUE STATE — render the 14 axes with their values
  const span_mb = (c.end_bp - c.start_bp) / 1e6;
  tierContent.innerHTML = `
    <div style="padding: 32px 28px; max-width: 1080px; margin: 0 auto;">
      <h3 style="margin: 0 0 6px; font-family: var(--serif); font-weight: 500;">
        14-axis classification
      </h3>
      <div style="color: var(--ink-dim); font-size: 12px; margin-bottom: 16px;">
        candidate ${c.id ? c.id.replace(/^cand_/, '') : '?'} · ${c.chrom}
        ${(c.start_bp/1e6).toFixed(2)}–${(c.end_bp/1e6).toFixed(2)} Mb · ${span_mb.toFixed(2)} Mb · K=${c.K}
      </div>
      ${_renderTierAxesGrid(candEntry)}
      <div style="margin-top: 20px; padding: 12px 16px; background: var(--panel);
                  border: 1px solid var(--rule); border-radius: 4px;
                  font-size: 11px; color: var(--ink-dim); line-height: 1.5;">
        <div><b>Source</b>: <code>state.data.final_classification[${candKey || '?'}]</code></div>
        <div><b>Schema</b>: SCHEMA_V2.md §19 (planned schema 2.16 — cluster-side emit) · cross-references §13 (three-axis evidence framework).</div>
      </div>
    </div>
  `;
}


// -----------------------------------------------------------------------------
// Persist subview-toggle helper. Called from the wire-buttons hook (legacy
// 63547) when the user clicks the karyotype/tier chip. Exposed here so the
// merge chat can stitch it into the dispatcher.
// -----------------------------------------------------------------------------
export function setKaryoSubview(next) {
  if (next !== 'karyotype' && next !== 'tier') return;
  karyoState.subview = next;
  try {
    if (typeof localStorage !== 'undefined') {
      localStorage.setItem(_CAND_SUBVIEW_KEY, next);
    }
  } catch (_) { /* no-op */ }
  renderCandidateKaryotype();
}


// =============================================================================
// state references not in SLOT_REGISTRY
// =============================================================================
//
// state.candidate            — IS in SLOT_REGISTRY (cross_atlas).
// state.data                 — IS in SLOT_REGISTRY (transient).
// state.data.classification  — sub-field, §9 cluster-emit layer
//                               (loaded into state.data by data IO; not a
//                                top-level slot).
// state.data.final_classification
//                            — sub-field of state.data; the R-pipeline
//                              output keyed by candidate_id. Read-only.
//                              NOT a SLOT_REGISTRY slot; it's a layer
//                              attached to the precomp JSON.
// state.activeMode           — referenced inside _renderCandidateKaryotypeBody
//                              (filtered out as TODO_MISSING).
//                              Legacy two-mode toggle ("simple" vs "detailed").
//                              NOT in current SLOT_REGISTRY.
//
// Recommendation for merge chat: state.activeMode is a candidate for being
// promoted to SLOT_REGISTRY as a `persisted` slot (key `inversion_atlas.activeMode`).
// =============================================================================
