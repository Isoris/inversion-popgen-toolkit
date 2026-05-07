// Atlas/inversion_discovery/page19.js
//
// Inversion-negative regions catalogue — complement of page3.
// NO JS handlers exist in legacy yet — page is pure HTML scaffold (see page19.html).
// Loaders (nrLoadBtn, nrExportCsvBtn, nrTableSlot, nrSummaryCards) will be authored fresh.
//
// Extracted LITERALLY from legacy/Inversion_atlas.html.
// Per HANDOFF_BATCH_1: state is now passed as first arg (was a global in legacy).
// All other unresolved references are flagged TODO_MISSING — the merge chat
// will resolve them once all 5 batches land.
//
// Entry points (in extraction order):
//

import { contextFromState, clusterL2, ClusterCache } from '../shared/per_l2_cluster.js';
import { hetRateColor } from '../shared/het_rate.js';
import { alignLabels, hungarianChainProjection, concordanceMatrix } from '../shared/hungarian.js';
import { buildContingency, computeARI, computeNMI, cramersV } from '../shared/contingency.js';
import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

// ---------------------------------------------------------------------------
// Extracted bodies
// ---------------------------------------------------------------------------
