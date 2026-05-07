// shared/state.js
// =====================================================================
// Canonical state contract for sub-atlases.
//
// The legacy Inversion_atlas.html uses one big global `state` object
// with ~110 declared slots and ~120 more runtime-created cache slots
// (~232 distinct state.X references in the source). This module does
// NOT try to redeclare all of them — that's intra-atlas business. What
// it DOES do:
//
//   1. makeState(opts)            → factory for a fresh state object,
//                                    with the sub-atlas's own initial slots
//                                    plus the canonical defaults from §1.
//   2. CROSS_ATLAS_SLOTS           → the small subset that round-trips
//                                    between sub-atlases via review/<workflow>/sessions/.
//   3. serializeState(state)       → JSON-friendly snapshot of cross-atlas slots
//   4. deserializeState(snapshot)  → hydrate Set/Map/typed-array fields back
//   5. mergeStateFromSession(state, snapshot) → apply a hydrated snapshot to a state
//
// The contract: each sub-atlas owns its full state in memory. Cross-
// atlas slots get written to review JSON when changed (so other atlases
// can pick them up); that's the only state plumbing across HTMLs.
//
// Slot inventory taken from legacy Inversion_atlas.html (line 9344..9789).
// =====================================================================

// ---------------------------------------------------------------------
// 1. CANONICAL SLOT DEFAULTS
//
// These are the slots from the legacy state declaration. We don't use
// them all in every sub-atlas; we list them here so:
//
//   - serialization knows the type of each slot (for typed-array reconstruction)
//   - sub-atlases can opt into a curated subset by passing `slots:` to makeState
//   - the auto-promote/review pages don't accidentally reinvent slot names
//
// Where the legacy code would do `state.k = 3`, the canonical default
// is `3`. Where the legacy code creates a slot lazily (e.g.
// `state.__hetRateCache = new Map()`), the slot is initialised to its
// empty-state default here so first reads don't return undefined.
//
// Slot tags drive serialization behaviour:
//   - 'persisted'   — written to localStorage; not sent to other atlases
//   - 'cross_atlas' — written to data/review/<workflow>/sessions/ AND
//                     read back on next session in any sub-atlas
//   - 'transient'   — never serialized (UI state, caches, animation flags)
//   - 'derived'     — computed from other slots; never serialized
// ---------------------------------------------------------------------

/**
 * Slot inventory. Each entry: { default, tag, type? }.
 *   - default:  initial value at sub-atlas startup
 *   - tag:      'persisted' | 'cross_atlas' | 'transient' | 'derived'
 *   - type:     'int' | 'float' | 'bool' | 'string' | 'array' | 'object'
 *               | 'Set' | 'Map' | 'Float32Array' | 'Int8Array' | 'Int32Array' | 'Uint32Array'
 *               | 'function' | 'null'
 *   - persistKey?: localStorage key for 'persisted' slots
 *
 * Inventory is incomplete by design — we add slots as sub-atlases are
 * extracted. Adding a slot here without a sub-atlas using it is a
 * waste; missing one is a bug. Use the legacy file as ground truth.
 */
export const SLOT_REGISTRY = Object.freeze({

  // -------------------------------------------------------------------
  // Cross-atlas: the small set that round-trips between HTMLs.
  // -------------------------------------------------------------------

  /** Active candidate at the top of the scrubber. May be null. */
  candidate:                { default: null, tag: 'cross_atlas',
                              type: 'object',
                              note: 'object with { id, chrom, start_bp, end_bp, ... }' },

  /** Full candidate list for the current chromosome. */
  candidateList:            { default: [], tag: 'cross_atlas',
                              type: 'array',
                              note: 'array of candidate objects; updated by auto-promote and manual flag' },

  /** Detailed candidate variants (legacy parallel structure). */
  candidate_detailed:       { default: null, tag: 'cross_atlas',
                              type: 'object' },
  candidateList_detailed:   { default: [], tag: 'cross_atlas',
                              type: 'array' },

  /** Per-sample karyotype lock. When non-null, K-means labels frozen to a reference L2. */
  lockedLabels:             { default: null, tag: 'cross_atlas',
                              type: 'Int8Array',
                              note: 'Int8Array[n_samples] | null' },
  lockedRefL2:              { default: null, tag: 'cross_atlas',
                              type: 'int' },

  /** User-defined manual groups (from lasso etc.). */
  manualGroups:             { default: null, tag: 'cross_atlas',
                              type: 'object',
                              note: 'object whose values are sample-index arrays' },

  /** Tracked samples (CGA tags pinned by the user). */
  tracked:                  { default: [], tag: 'cross_atlas',
                              type: 'array',
                              note: 'array of sample indices (ints)' },

  /** Sweep auto-promote results (turn 133-134). */
  l2SweepResult:            { default: null, tag: 'cross_atlas',
                              type: 'object' },

  /** Active-samples mask (turn 128 AS1). */
  activeSampleSet:          { default: null, tag: 'cross_atlas',
                              type: 'Set',
                              note: 'Set<sample_idx> | null (null = all active)' },
  activeSampleReasons:      { default: null, tag: 'cross_atlas',
                              type: 'Map',
                              note: 'Map<sample_idx, string>' },
  activeSampleRules:        { default: [], tag: 'cross_atlas',
                              type: 'array' },

  /** Inversion-review tab decisions from the review atlas. */
  candidate_review_decisions: { default: {}, tag: 'cross_atlas',
                                type: 'object',
                                note: 'object keyed by candidate_id → {decision, reason, by, at}' },

  /** Locked karyotype groups (manual locks; H1/H1, H1/H2, H2/H2). */
  locked_karyotype_groups:  { default: {}, tag: 'cross_atlas',
                              type: 'object' },

  // -------------------------------------------------------------------
  // Persisted: in localStorage on this sub-atlas. Not cross-atlas.
  // -------------------------------------------------------------------

  /** Active K (default 3; can be 6). */
  k:                        { default: 3, tag: 'persisted',
                              type: 'int',
                              persistKey: 'pca_scrubber_v3.k' },
  kMode:                    { default: 'fixed', tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.kMode' },
  kRange:                   { default: [2, 5], tag: 'persisted', type: 'array',
                              persistKey: 'pca_scrubber_v3.kRange' },

  /** Heat-map / signal-detection thresholds from the L3 controls. */
  silThreshold:             { default: 0.45, tag: 'persisted', type: 'float',
                              persistKey: 'pca_scrubber_v3.silThreshold' },
  mergeThr:                 { default: 0.85, tag: 'persisted', type: 'float',
                              persistKey: 'pca_scrubber_v3.mergeThr' },
  alpha:                    { default: 0.05, tag: 'persisted', type: 'float',
                              persistKey: 'pca_scrubber_v3.alpha' },
  minNGroup:                { default: 5,    tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.minNGroup' },
  minNWin:                  { default: 3,    tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.minNWin' },
  aggMethod:                { default: 'mean_pc1', tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.aggMethod' },

  /** L3 panel layout / coloring controls. */
  l3Layout:                 { default: 'leftright', tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.l3Layout' },
  l3ColorMode:              { default: 'shared',    tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.l3ColorMode' },
  l3KMode:                  { default: 'k3',        tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.l3KMode' },
  l3HetColoring:            { default: false,       tag: 'persisted', type: 'bool',
                              persistKey: 'pca_scrubber_v3.l3HetColoring' },
  l3BcScope:                { default: 'carousel',  tag: 'persisted', type: 'string',
                              persistKey: 'pca_scrubber_v3.l3BcScope' },
  l3ShowDetails:            { default: false,       tag: 'persisted', type: 'bool',
                              persistKey: 'pca_scrubber_v3.l3ShowDetails' },
  l2SweepEnabled:           { default: false,       tag: 'persisted', type: 'bool',
                              persistKey: 'pca_scrubber_v3.l2SweepEnabled' },

  /** Per-sample-lines panel toggles. */
  linesColorMode:           { default: 'kmeans',    tag: 'persisted', type: 'string',
                              persistKey: 'inversion_atlas.linesColorMode' },
  linesSnpDensityMode:      { default: 'off',       tag: 'persisted', type: 'string',
                              persistKey: 'inversion_atlas.linesSnpDensityMode' },
  linesTransRateOn:         { default: false,       tag: 'persisted', type: 'bool',
                              persistKey: 'inversion_atlas.linesTransRateOn' },
  linesRegimeBreadthOn:     { default: false,       tag: 'persisted', type: 'bool',
                              persistKey: 'inversion_atlas.linesRegimeBreadthOn' },
  linesLineageStripOn:      { default: true,        tag: 'persisted', type: 'bool',
                              persistKey: 'inversion_atlas.linesLineageStripOn' },
  linesPanelCandidateBands: { default: true,        tag: 'persisted', type: 'bool',
                              persistKey: 'inversion_atlas.linesPanelCandidateBands' },

  /** G-panel popup state. */
  gPanelTab:                { default: 'manual', tag: 'persisted', type: 'string',
                              persistKey: 'inversion_atlas.gPanelTab' },
  gPanelInheritanceThreshold: { default: 0.15, tag: 'persisted', type: 'float',
                              persistKey: 'inversion_atlas.gPanelInheritanceThreshold' },

  /** Band-trace toggles (turn 161). */
  bandTraceOn:              { default: false, tag: 'persisted', type: 'bool',
                              persistKey: 'inversion_atlas.bandTraceOn' },
  bandTraceFishSet:         { default: null,  tag: 'cross_atlas', type: 'Set',
                              note: 'Set<sample_idx> | null — single-cohort tracker' },

  /** Panel resize heights. */
  simPanelH:                { default: 520, tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.simPanelH' },
  zPanelH:                  { default: 100, tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.zPanelH' },
  linesPanelH:              { default: 180, tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.linesPanelH' },
  pcaPanelH:                { default: 280, tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.pcaPanelH' },

  // -------------------------------------------------------------------
  // Transient: not persisted at all. Reset on reload.
  // -------------------------------------------------------------------

  data:                     { default: null,  tag: 'transient', type: 'object',
                              note: 'precomp main JSON for the active chromosome' },
  cur:                      { default: 0,     tag: 'transient', type: 'int',
                              note: 'scrubber position (window index)' },
  playing:                  { default: false, tag: 'transient', type: 'bool' },
  playTimer:                { default: null,  tag: 'transient', type: 'function' },
  pc1Sign:                  { default: null,  tag: 'derived',   type: 'Float32Array' },
  layersPresent:            { default: null,  tag: 'derived',   type: 'Set' },
  trailN:                   { default: 15,    tag: 'transient', type: 'int' },
  trailOn:                  { default: true,  tag: 'transient', type: 'bool' },
  flipPC1:                  { default: true,  tag: 'transient', type: 'bool' },

  // Sim_mat display
  simScale:                 { default: null,  tag: 'transient', type: 'string' },
  pdfStyle:                 { default: true,  tag: 'transient', type: 'bool' },
  zCollapsed:               { default: false, tag: 'transient', type: 'bool' },
  pcaCollapsed:             { default: false, tag: 'transient', type: 'bool' },
  l3Collapsed:              { default: false, tag: 'transient', type: 'bool' },

  // Lasso state — UI-only, transient
  pcaLassoActive:           { default: false, tag: 'transient', type: 'bool' },
  linesLassoActive:         { default: false, tag: 'transient', type: 'bool' },
  linesLassoRect:           { default: null,  tag: 'transient', type: 'object' },
  linesLassoCommitted:      { default: null,  tag: 'transient', type: 'object' },
  linesLassoSelected:       { default: [],    tag: 'transient', type: 'array' },

  spotlight:                { default: null,  tag: 'transient', type: 'int' },
  spotlightTrackedAll:      { default: false, tag: 'transient', type: 'bool' },

  // G-panel UI flag
  gPanelOpen:               { default: false, tag: 'transient', type: 'bool' },

  // Anchor (for color-locking against a reference L2)
  trackingAnchor:           { default: null,  tag: 'transient', type: 'object' },
  anchorConcord:            { default: null,  tag: 'derived',   type: 'Float32Array' },
  anchorStripH:             { default: 18,    tag: 'persisted', type: 'int',
                              persistKey: 'pca_scrubber_v3.anchorStripH' },

  // Lineage compute results (turn 130)
  l2SweepCacheKey:          { default: null,  tag: 'derived',   type: 'string' },

  // Band-trace runtime
  bandTraceCache:           { default: null,  tag: 'transient', type: 'object' },
  bandTraceCacheKey:        { default: null,  tag: 'derived',   type: 'string' },

  // -------------------------------------------------------------------
  // Derived caches — never serialized, computed on first access.
  // The double-underscore prefix is the legacy convention.
  // -------------------------------------------------------------------

  __hetRateCache:           { default: null,  tag: 'transient', type: 'Map',
                              note: 'lazy: created as new Map() on first access' },
  __dosageChunkCache:       { default: null,  tag: 'transient', type: 'Map' },
  __dosageHm:               { default: null,  tag: 'transient', type: 'object' },
  __linesCache:             { default: null,  tag: 'transient', type: 'object' },
  __linesGeom:              { default: null,  tag: 'transient', type: 'object' },
  __pcaPlotRect:            { default: null,  tag: 'transient', type: 'object' },
  __pcaScreenXY:            { default: null,  tag: 'transient', type: 'Float32Array' },
  __boundaries:             { default: null,  tag: 'transient', type: 'object' },
  __bdByL2:                 { default: null,  tag: 'transient', type: 'Map' },
  __candJumperMarkers:      { default: null,  tag: 'transient', type: 'array' },

  // -------------------------------------------------------------------
  // viewControls (composite — handled specially)
  // -------------------------------------------------------------------
  viewControls:             { default: () => ({
                                pcaXY: ['pc1', 'pc2'],
                                linesYsources: ['pc1'],
                                linked: true,
                              }), tag: 'persisted', type: 'object',
                              persistKey: 'inversion_atlas.viewControls' },

  // The total list is incomplete; missing slots fall through to undefined
  // (which is fine — the legacy code does the same on lazy slots).
});

// ---------------------------------------------------------------------
// 2. CROSS_ATLAS_SLOTS (derived)
// The names of the slots tagged 'cross_atlas' — the set written to
// review/<workflow>/sessions/. Used by serializeState/deserializeState.
// ---------------------------------------------------------------------
export const CROSS_ATLAS_SLOTS = Object.freeze(
  Object.entries(SLOT_REGISTRY)
    .filter(([_, def]) => def.tag === 'cross_atlas')
    .map(([name]) => name)
);

// ---------------------------------------------------------------------
// 3. PERSISTED_SLOTS (derived)
// Slots written to localStorage. The persistKey is the localStorage key.
// ---------------------------------------------------------------------
export const PERSISTED_SLOTS = Object.freeze(
  Object.entries(SLOT_REGISTRY)
    .filter(([_, def]) => def.tag === 'persisted' && def.persistKey)
    .map(([name, def]) => ({ name, key: def.persistKey, type: def.type }))
);

// ---------------------------------------------------------------------
// 4. makeState — factory
// ---------------------------------------------------------------------

/**
 * Build a fresh state object for a sub-atlas. By default, populates ALL
 * slots in SLOT_REGISTRY with their defaults. Pass `slots: ['a','b']`
 * to only populate a curated subset (the rest stay undefined; legacy
 * behaviour for un-init slots).
 *
 * @param {object} [opts]
 * @param {string[]} [opts.slots]  Subset of slot names; default = all
 * @returns {object} a fresh state object
 */
export function makeState({ slots } = {}) {
  const out = {};
  const names = slots || Object.keys(SLOT_REGISTRY);
  for (const name of names) {
    const def = SLOT_REGISTRY[name];
    if (!def) {
      // Unknown slot — treat as null (legacy uninitialised semantics)
      out[name] = null;
      continue;
    }
    out[name] = (typeof def.default === 'function')
      ? def.default()
      : cloneDefault(def.default);
  }
  return out;
}

// Deep-clone a default so two calls to makeState() get fresh array/object refs
function cloneDefault(v) {
  if (v === null || v === undefined) return v;
  if (Array.isArray(v)) return v.slice();
  if (v instanceof Set) return new Set(v);
  if (v instanceof Map) return new Map(v);
  if (typeof v === 'object') return { ...v };
  return v;   // primitives
}

// ---------------------------------------------------------------------
// 5. Serializer — convert cross-atlas slots to JSON-friendly form
// ---------------------------------------------------------------------

/**
 * Convert state's cross-atlas slots into a JSON-serializable snapshot.
 * Handles Set, Map, Int8Array/Float32Array/etc. Returns an object
 * suitable for JSON.stringify().
 *
 *   { schema_version: 1,
 *     workflow: 'inversion',
 *     created_at: '...',
 *     slots: { candidate: {...}, candidateList: [...], lockedLabels: [...] }
 *   }
 *
 * @param {object} state
 * @param {object} [opts]
 * @param {string} [opts.workflow]
 * @returns {object}
 */
export function serializeState(state, { workflow = 'inversion' } = {}) {
  const slots = {};
  for (const name of CROSS_ATLAS_SLOTS) {
    if (!(name in state)) continue;
    const v = state[name];
    if (v === null || v === undefined) {
      slots[name] = null;
      continue;
    }
    slots[name] = encodeValue(v);
  }
  return {
    schema_version: 1,
    workflow,
    created_at: new Date().toISOString(),
    slots,
  };
}

/**
 * Inverse of serializeState. Reconstruct typed-array / Set / Map fields
 * back into their runtime forms.
 *
 * @param {object} snapshot  output of serializeState
 * @returns {object}         partial state object (only cross-atlas slots)
 */
export function deserializeState(snapshot) {
  if (!snapshot || typeof snapshot !== 'object') return {};
  if (snapshot.schema_version !== 1) {
    console.warn(`[state] unknown snapshot schema_version ${snapshot.schema_version}; treating as v1`);
  }
  const out = {};
  const slots = snapshot.slots || {};
  for (const name of CROSS_ATLAS_SLOTS) {
    if (!(name in slots)) continue;
    const def = SLOT_REGISTRY[name];
    out[name] = decodeValue(slots[name], def && def.type);
  }
  return out;
}

/**
 * Apply a snapshot's slots onto an existing state object (in place).
 * Use after the sub-atlas has done its initial load and we're picking
 * up cross-atlas state from a sibling.
 *
 * @param {object} state
 * @param {object} snapshot
 */
export function mergeStateFromSession(state, snapshot) {
  const partial = deserializeState(snapshot);
  for (const name of Object.keys(partial)) {
    state[name] = partial[name];
  }
  return state;
}

// ---------------------------------------------------------------------
// 6. Persistence helpers (localStorage)
// ---------------------------------------------------------------------

/**
 * Read all 'persisted' slot values from localStorage, applying type
 * coercion. Used at sub-atlas startup to hydrate persisted state.
 *
 * @param {Storage} [storage]  defaults to localStorage; pass mock for testing
 * @returns {object}            partial state with persisted slot values
 */
export function readPersistedSlots(storage) {
  const ls = storage || (typeof localStorage !== 'undefined' ? localStorage : null);
  if (!ls) return {};
  const out = {};
  for (const { name, key, type } of PERSISTED_SLOTS) {
    let raw;
    try { raw = ls.getItem(key); } catch { raw = null; }
    if (raw === null || raw === undefined) continue;
    out[name] = coerceFromStorage(raw, type);
  }
  return out;
}

/**
 * Write one slot's value to localStorage if it's persisted.
 * Caller invokes this whenever a persisted slot changes.
 */
export function writePersistedSlot(name, value, storage) {
  const ls = storage || (typeof localStorage !== 'undefined' ? localStorage : null);
  if (!ls) return false;
  const def = SLOT_REGISTRY[name];
  if (!def || def.tag !== 'persisted' || !def.persistKey) return false;
  try {
    ls.setItem(def.persistKey, encodeForStorage(value, def.type));
    return true;
  } catch { return false; }
}

// ---------------------------------------------------------------------
// Type coercion helpers — encode/decode values for JSON & localStorage
// ---------------------------------------------------------------------

function encodeValue(v) {
  if (v === null || v === undefined) return null;
  if (v instanceof Set)  return { __type: 'Set',  items: Array.from(v) };
  if (v instanceof Map)  return { __type: 'Map',  items: Array.from(v.entries()) };
  if (v instanceof Int8Array)    return { __type: 'Int8Array',    items: Array.from(v) };
  if (v instanceof Int32Array)   return { __type: 'Int32Array',   items: Array.from(v) };
  if (v instanceof Uint32Array)  return { __type: 'Uint32Array',  items: Array.from(v) };
  if (v instanceof Float32Array) return { __type: 'Float32Array', items: Array.from(v) };
  if (v instanceof Float64Array) return { __type: 'Float64Array', items: Array.from(v) };
  if (Array.isArray(v))  return v.map(encodeValue);
  if (typeof v === 'object') {
    const out = {};
    for (const k of Object.keys(v)) out[k] = encodeValue(v[k]);
    return out;
  }
  return v;
}

function decodeValue(v, hintedType) {
  if (v === null || v === undefined) return v;
  if (Array.isArray(v)) return v.map(x => decodeValue(x));
  if (typeof v === 'object' && v.__type) {
    switch (v.__type) {
      case 'Set':            return new Set(v.items);
      case 'Map':            return new Map(v.items);
      case 'Int8Array':      return Int8Array.from(v.items);
      case 'Int32Array':     return Int32Array.from(v.items);
      case 'Uint32Array':    return Uint32Array.from(v.items);
      case 'Float32Array':   return Float32Array.from(v.items);
      case 'Float64Array':   return Float64Array.from(v.items);
      default:               return v;   // unknown __type: return as-is
    }
  }
  if (typeof v === 'object') {
    const out = {};
    for (const k of Object.keys(v)) out[k] = decodeValue(v[k]);
    return out;
  }
  return v;
}

function coerceFromStorage(raw, type) {
  if (raw === '' || raw === null || raw === undefined) return null;
  switch (type) {
    case 'int':    { const n = parseInt(raw, 10);    return Number.isFinite(n) ? n : null; }
    case 'float':  { const n = parseFloat(raw);      return Number.isFinite(n) ? n : null; }
    case 'bool':   return raw === 'true' || raw === '1';
    case 'string': return String(raw);
    case 'array':
    case 'object':
      try { return JSON.parse(raw); } catch { return null; }
    default:       return raw;
  }
}

function encodeForStorage(v, type) {
  switch (type) {
    case 'int':
    case 'float':  return String(v);
    case 'bool':   return v ? 'true' : 'false';
    case 'string': return String(v ?? '');
    case 'array':
    case 'object': return JSON.stringify(v);
    default:       return String(v ?? '');
  }
}

// ---------------------------------------------------------------------
// Console-debug exposures
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._SLOT_REGISTRY      = SLOT_REGISTRY;
  window._CROSS_ATLAS_SLOTS  = CROSS_ATLAS_SLOTS;
  window._serializeState     = serializeState;
  window._deserializeState   = deserializeState;
}
