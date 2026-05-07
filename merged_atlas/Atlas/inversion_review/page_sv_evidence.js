// =============================================================================
// inversion_review/page_sv_evidence.js — "5b SV evidence" tab
// =============================================================================
// Stage:        refinement
// Legacy DOM:   <div id="page_sv_evidence"> with inner #sv_evidence_root
//               (legacy lines 9329–9331)
// Renderer:     window.AtlasSVEvidence — loaded from js/atlas_sv_evidence.js
//               (external file, NOT inlined in legacy Inversion_atlas.html)
// Spec:         specs_todo/SPEC_sv_evidence_page.md
//
// What this page does
// -------------------
// Read-only view of SV calls clustered around a candidate's boundaries,
// scored against the karyotype groups. Loads `json/sv_genotype_counts/<cid>.json`
// per candidate. Renders three things stacked top-to-bottom:
//   - SV table (columns: caller, type, chrom, pos, len, support, gt-counts)
//   - UpSet plot of caller intersections
//   - dosage heatmap (step 6 spec)
//
// Extraction notes (Batch 2)
// --------------------------
// Unlike the other review pages, page_sv_evidence is NOT defined inline in
// the legacy HTML. The whole renderer lives in `js/atlas_sv_evidence.js`,
// loaded as an external script and exposed as `window.AtlasSVEvidence`.
//
// Inside legacy Inversion_atlas.html the only references are:
//   - the mount slot DOM (lines 9329–9331)
//   - the script-tag inventory entry (line 54981) listing
//     `{ file: 'js/atlas_sv_evidence.js', global: 'AtlasSVEvidence' }`
//
// So this module:
//   1. Provides the lifecycle adapter the merge chat will wire into the
//      tab dispatcher (init when the tab first activates, loadCandidate
//      whenever state.candidate changes).
//   2. Treats `window.AtlasSVEvidence` as the contract surface — the
//      external JS file is shipped alongside the atlas as-is.
// =============================================================================


// -----------------------------------------------------------------------------
// External-file dep — the page renderer lives in atlas_sv_evidence.js
// -----------------------------------------------------------------------------
//
// TODO_MISSING(AtlasSVEvidence.init)            — js/atlas_sv_evidence.js
// TODO_MISSING(AtlasSVEvidence.loadCandidate)   — js/atlas_sv_evidence.js
// TODO_MISSING(AtlasSVEvidence.destroy)         — js/atlas_sv_evidence.js
//
// These are not in the legacy HTML at all — the merge chat needs to either
// (a) bundle js/atlas_sv_evidence.js as a sibling module, or
// (b) port its contents into shared/ if we go fully ES-module.
//
// Decision (BATCH_2_NOTES): kept as external-file dep for now. Promoting
// the SV evidence renderer to ES modules is out of scope for batch 2 — it
// touches IndexedDB caching and the UpSet plot lib, which deserves its own
// extraction pass.
// -----------------------------------------------------------------------------


/**
 * Show the SV-evidence page. Initialises the AtlasSVEvidence module on
 * first activation and dispatches a candidate-load if one is selected.
 *
 * Called by the tab dispatcher when the user clicks the "5b SV evidence"
 * chip. Idempotent: subsequent activations only re-render if the active
 * candidate has changed.
 *
 * @param {object} state   shared state object (not used directly here;
 *                         AtlasSVEvidence reads from window.state in legacy)
 * @returns {void}
 */
export function showSvEvidencePage(state) {
  if (typeof window === 'undefined') return;
  const mod = /** @type {any} */ (window).AtlasSVEvidence;
  if (!mod) {
    // External script not loaded — leave a friendly empty-state in the slot
    const root = document.getElementById('sv_evidence_root');
    if (root && !root.__svInitFailed) {
      root.innerHTML =
        '<div style="padding:32px 28px;color:var(--ink-dim);font-family:var(--mono);font-size:11.5px;line-height:1.6;">' +
        '<div style="font-weight:600;color:var(--ink);margin-bottom:8px;">SV evidence module not loaded</div>' +
        '<p style="margin:0;max-width:720px;">' +
        'The <code>AtlasSVEvidence</code> renderer ships in <code>js/atlas_sv_evidence.js</code>, ' +
        'an external script that wasn\'t included in this build. Drop the file alongside ' +
        '<code>inversion_review.html</code> and reload to enable this page.' +
        '</p></div>';
      root.__svInitFailed = true;
    }
    return;
  }

  // First-activation init. The legacy module owns its own init guard, but
  // we double-guard here so the dispatcher never has to care.
  if (!mod.__pageInitDone) {
    try {
      // TODO_MISSING(AtlasSVEvidence.init)
      mod.init({ rootSelector: '#sv_evidence_root' });
      mod.__pageInitDone = true;
    } catch (err) {
      console.warn('[page_sv_evidence] AtlasSVEvidence.init failed:', err);
      return;
    }
  }

  // Re-load whenever the active candidate changed since last activation.
  const cand = state && state.candidate;
  const cid = cand ? (cand.id || null) : null;
  if (mod.__lastCid !== cid) {
    try {
      // TODO_MISSING(AtlasSVEvidence.loadCandidate)
      mod.loadCandidate(cid);
      mod.__lastCid = cid;
    } catch (err) {
      console.warn('[page_sv_evidence] AtlasSVEvidence.loadCandidate failed:', err);
    }
  }
}


/**
 * Hide / tear down the SV-evidence page when the user navigates away.
 * Optional — most legacy code just leaves the DOM in place. Provided so the
 * dispatcher can call it symmetrically alongside show().
 *
 * @returns {void}
 */
export function hideSvEvidencePage() {
  if (typeof window === 'undefined') return;
  const mod = /** @type {any} */ (window).AtlasSVEvidence;
  if (mod && typeof mod.destroy === 'function') {
    try {
      // TODO_MISSING(AtlasSVEvidence.destroy)
      mod.destroy();
    } catch (_) { /* swallow */ }
  }
}


// =============================================================================
// state references not in SLOT_REGISTRY
// =============================================================================
//
// state.candidate    — IS in SLOT_REGISTRY (cross_atlas). The renderer reads
//                      `.id` to fetch json/sv_genotype_counts/<cid>.json.
//
// (AtlasSVEvidence also reads its own internal state via window.state in the
//  legacy build — that's the closure-scoped legacy global, not anything new.
//  No new SLOT_REGISTRY entries are introduced by this page.)
// =============================================================================
