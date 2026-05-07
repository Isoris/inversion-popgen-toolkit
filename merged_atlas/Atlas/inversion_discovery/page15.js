// Atlas/inversion_discovery/page15.js
//
// Same layout as page1/page12 but driven by GHSL haplotype divergence.
// Currently almost-pure empty-state — only the layer-status badge has a handler.
//
// Extracted LITERALLY from legacy/Inversion_atlas.html.
// Per HANDOFF_BATCH_1: state is now passed as first arg (was a global in legacy).
// All other unresolved references are flagged TODO_MISSING — the merge chat
// will resolve them once all 5 batches land.
//
// Entry points (in extraction order):
  // _refreshGhslLayerStatus() — legacy lines 53065-53081
//

import { contextFromState, clusterL2, ClusterCache } from '../shared/per_l2_cluster.js';
import { hetRateColor } from '../shared/het_rate.js';
import { alignLabels, hungarianChainProjection, concordanceMatrix } from '../shared/hungarian.js';
import { buildContingency, computeARI, computeNMI, cramersV } from '../shared/contingency.js';
import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

// ---------------------------------------------------------------------------
// Extracted bodies
// ---------------------------------------------------------------------------

// --- _refreshGhslLayerStatus() — legacy lines 53065-53081 ---
export function _refreshGhslLayerStatus(state) {
  if (typeof document === 'undefined') return;
  const indicators = document.querySelectorAll('[data-gh-layer]');
  if (!indicators || indicators.length === 0) return;
  indicators.forEach(node => {
    const layerName = node.dataset && node.dataset.ghLayer;
    if (!layerName) return;
    const present = !!(state.layersPresent && state.layersPresent.has(layerName));
    if (present) {
      node.textContent = '🟢 loaded';
      node.style.color = 'var(--good)';
    } else {
      node.textContent = '⚪ not loaded';
      node.style.color = 'var(--ink-dimmer)';
    }
  });
}
