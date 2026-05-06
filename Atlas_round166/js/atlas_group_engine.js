// =============================================================================
// atlas_group_engine.js — Turn 2 of chat A
// =============================================================================
//
// Group composition data layer for the atlas. No UI yet. Provides:
//
//   1. window.popgen — global namespace with all the public methods
//   2. Group expression AST + evaluator
//   3. Dimension extractors that read existing atlas state
//      (state.candidate.fish_calls, state.data.ancestry_sample, family palette,
//      tracked samples, lasso history, etc.)
//   4. Slot system — named groups (A, B, C, ..., up to 10 — engine F's MAX_GROUPS)
//      that compile to {name: [sample_ids]} ready for the popstats server
//   5. Cursor-region resolver — '1w' | '5w' | '10w' | 'L2' | 'candidate' |
//      'chrom' | null → {start_bp, end_bp} or null for whole chromosome
//   6. Lasso/selection history log — every lasso, rectangle drag, manual confirm
//      pushes an entry with timestamp + page-of-origin + cursor-context
//   7. Snapshot system — Shift+S captures (chrom, cursor, expressions, last
//      result, note, ts); Shift+R opens the snapshot list
//   8. localStorage persistence — keyed by cohort, scoped per-chrom for
//      cursor-derived state, global for saved expressions and snapshots
//
// This file is intentionally UI-free. Run window.popgen.demo() in the browser
// console to exercise every code path against the current atlas state.
// =============================================================================

(function (root, factory) {
  // UMD-ish: works as an inline <script> in the atlas HTML and also as a
  // Node module loaded by jsdom for unit tests.
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgen = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const POPGEN_LS_PREFIX        = 'inversion_atlas.popgen.';
  const POPGEN_LS_EXPRESSIONS   = POPGEN_LS_PREFIX + 'expressions';
  const POPGEN_LS_SAVED_SETS    = POPGEN_LS_PREFIX + 'savedSets';
  const POPGEN_LS_SNAPSHOTS     = POPGEN_LS_PREFIX + 'snapshots';
  const POPGEN_LS_LASSO_HIST    = POPGEN_LS_PREFIX + 'lassoHistory';
  const POPGEN_LS_SLOTS         = POPGEN_LS_PREFIX + 'slots';

  // Engine F's MAX_GROUPS. Hard cap on simultaneous slots in one popstats
  // request. The atlas can have more SAVED expressions than this; only the
  // ACTIVE ones (mapped to slot positions) get sent in a single Compute call.
  const MAX_SLOTS               = 10;

  // Lasso history is unbounded by intent (every action is recorded), but we
  // ring-buffer to localStorage to bound storage. Pinned entries persist.
  const LASSO_HISTORY_MAX       = 200;

  // Snapshot retention. Pinned snapshots never auto-evict.
  const SNAPSHOT_MAX            = 100;

  // Cursor scopes that the user can pick. The resolver returns the bp window
  // covered by the current scope.
  const CURSOR_SCOPES = ['1w', '5w', '10w', 'L2', 'candidate', 'chrom', null];

  // ===========================================================================
  // Atlas state accessor
  // ===========================================================================
  // The engine reads the existing atlas global `state` object. Because we may
  // load before the atlas defines it (or in a Node test environment with no
  // atlas at all), every read goes through getAtlasState() which returns null
  // if the atlas isn't there. All extractors handle null gracefully.

  function getAtlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }

  // The cohort key is what scopes localStorage entries. Reuses the atlas's own
  // cohort identifier if present (set by AS1 in turn 128). Falls back to
  // chrom-bound when the cohort key isn't loaded yet.
  function getCohortKey() {
    const s = getAtlasState();
    if (!s || !s.data) return 'no_data';
    if (s.data.cohort_key) return String(s.data.cohort_key);
    if (s.data.dataset)   return String(s.data.dataset);
    if (s.data.chrom)     return 'chrom_' + String(s.data.chrom);
    return 'unknown';
  }

  // Map sample_idx (atlas's internal) ↔ sample_id (the canonical string used by
  // the server). Both directions. Falls back to numeric strings if the atlas
  // hasn't loaded a sample_meta yet.
  function sampleIdxToId(idx) {
    const s = getAtlasState();
    if (!s || !s.data) return String(idx);
    if (Array.isArray(s.data.samples) && idx >= 0 && idx < s.data.samples.length) {
      return String(s.data.samples[idx]);
    }
    if (s.data.sample_meta && Array.isArray(s.data.sample_meta.cga) &&
        idx >= 0 && idx < s.data.sample_meta.cga.length) {
      return String(s.data.sample_meta.cga[idx]);
    }
    return String(idx);
  }

  function sampleIdToIdx(id) {
    const s = getAtlasState();
    if (!s || !s.data) return -1;
    if (Array.isArray(s.data.samples)) {
      for (let i = 0; i < s.data.samples.length; i++) {
        if (String(s.data.samples[i]) === String(id)) return i;
      }
    }
    if (s.data.sample_meta && Array.isArray(s.data.sample_meta.cga)) {
      for (let i = 0; i < s.data.sample_meta.cga.length; i++) {
        if (String(s.data.sample_meta.cga[i]) === String(id)) return i;
      }
    }
    return -1;
  }

  function allSampleIds() {
    const s = getAtlasState();
    if (!s || !s.data) return [];
    if (Array.isArray(s.data.samples)) return s.data.samples.map(String);
    if (s.data.sample_meta && Array.isArray(s.data.sample_meta.cga)) {
      return s.data.sample_meta.cga.map(String);
    }
    if (typeof s.data.n_samples === 'number') {
      const out = [];
      for (let i = 0; i < s.data.n_samples; i++) out.push(String(i));
      return out;
    }
    return [];
  }

  // ===========================================================================
  // Detailed-mode H-label tables (mirrors atlas line 27892)
  // ===========================================================================
  // Operational diploid labels for K=3..6, ordered by median PC1 ascending.
  // These are NOT confirmed biological assignments — heterozygosity evidence
  // is required to claim e.g. that "H1/H2" is actually heterozygote-like.

  const BAND_LABELS_DETAILED = {
    3: ['H1/H1', 'H1/H2', 'H2/H2'],
    4: ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3'],
    5: ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3', 'H3/H3'],
    6: ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3', 'H2/H3', 'H3/H3'],
  };

  const BAND_LABELS_LEGACY = {
    3: ['g0', 'g1', 'g2'],
    4: ['g0', 'g1', 'g2', 'g3'],
    5: ['g0', 'g1', 'g2', 'g3', 'g4'],
    6: ['g0', 'g1', 'g2', 'g3', 'g4', 'g5'],
  };

  // Given a diploid label like "H1/H2", return [hap_a, hap_b] sorted.
  // Returns null for legacy g* labels (no haplotype meaning).
  function parseHaplotypePair(label) {
    if (typeof label !== 'string') return null;
    const m = /^H(\d+)\/H(\d+)$/.exec(label);
    if (!m) return null;
    const a = 'H' + m[1], b = 'H' + m[2];
    return a < b ? [a, b] : [b, a];
  }

  // Carrier predicates derived from a diploid label.
  function carriesHaplotype(label, hap) {
    const pair = parseHaplotypePair(label);
    if (!pair) return false;
    return pair[0] === hap || pair[1] === hap;
  }

  function isHomozygousFor(label, hap) {
    const pair = parseHaplotypePair(label);
    if (!pair) return false;
    return pair[0] === hap && pair[1] === hap;
  }

  // ===========================================================================
  // Dimension registry
  // ===========================================================================
  // A dimension produces, for every sample in the cohort, either:
  //   - a categorical label (string) — for `eq` / `in` / `not_in` predicates
  //   - a boolean — for `tracked` / direct toggle predicates
  //   - a number — for threshold predicates (future)
  //
  // Dimensions are CACHED. The cache key is (dimension_name, atlas_state_token)
  // where atlas_state_token is a coarse identifier for "what state version".
  // For per-candidate dimensions, the token incorporates the candidate id so
  // switching candidates invalidates the per-candidate dims but not the
  // global ones (family, ROH, etc.).
  //
  // Dimension naming convention:
  //   Global dims:           'family'  'relatedness_hub'  'roh_carrier'
  //                          'ancestry_dominantQ_K8'  'tracked'  'active_samples'
  //   Per-candidate dims:    'diploid_class@<candidate_id>'
  //                          'carries_H1@<candidate_id>'
  //                          'homozygous_H2@<candidate_id>'
  //                          'subband_stable@<candidate_id>'
  //                          'votes_consensus@<candidate_id>'
  //                          'fish_regime@<candidate_id>'  (raw regime int)
  //                          'fish_ambiguous@<candidate_id>' (boolean)
  //   Lasso/saved:           'lasso_<lasso_id>'  'saved_<saved_set_name>'

  const _dimCache = new Map();
  let _dimCacheToken = 0;

  function _bumpDimCache() { _dimCacheToken++; _dimCache.clear(); }

  // Atlas should call popgen.invalidateDimensions() whenever upstream data
  // changes (focal candidate switched, K changed, family palette reloaded,
  // etc.). The dimension extractors don't bother subscribing; the consumer
  // is responsible for refresh. This keeps the engine cheap — extractors only
  // run when an expression evaluator pulls them.
  function invalidateDimensions() { _bumpDimCache(); }

  function _cacheDim(name, fn) {
    const key = name + '@token=' + _dimCacheToken;
    if (_dimCache.has(key)) return _dimCache.get(key);
    const value = fn();
    _dimCache.set(key, value);
    return value;
  }

  // ---------------------------------------------------------------------------
  // Extractor: per-candidate diploid class (the H-system labels)
  // ---------------------------------------------------------------------------
  // Reads from cand.fish_calls (cand = state.candidate or any cand object).
  // Returns: Map<sample_id_string, label_string | null>
  // Labels are detailed-mode strings ('H1/H1', 'H1/H2', ...) when atlas is in
  // detailed mode, legacy ('g0', 'g1', ...) otherwise.

  function _extractDiploidClass(candObj, useDetailed) {
    if (!candObj || !Array.isArray(candObj.fish_calls)) return new Map();
    const labelsForK = useDetailed ? BAND_LABELS_DETAILED : BAND_LABELS_LEGACY;
    // Determine K: prefer candObj.k, fall back to atlas state.k, else 3.
    let K = candObj.k;
    const s = getAtlasState();
    if (!K && s) K = s.k;
    if (!K) K = 3;
    const labels = labelsForK[K] || labelsForK[3];
    const out = new Map();
    for (const fc of candObj.fish_calls) {
      if (!fc) continue;
      const id = fc.sample_id ? String(fc.sample_id) : sampleIdxToId(fc.sample_idx);
      const r = (typeof fc.regime === 'number') ? fc.regime : -1;
      if (r < 0 || r >= labels.length) {
        out.set(id, null);  // unassigned (fish_calls allows -1)
      } else {
        out.set(id, labels[r]);
      }
    }
    return out;
  }

  // Register a per-candidate dimension family for an arbitrary cand object.
  // Returns the diploidClass map so other extractors can derive carrier preds.
  function getDiploidClassMap(candId) {
    return _cacheDim('diploid_class@' + candId, () => {
      const cand = _resolveCandidateById(candId);
      if (!cand) return new Map();
      // Detect detailed-mode by looking at cand._system if set (turn 88+) else
      // fall back to atlas's current display mode (label vocabulary toggle).
      const useDetailed = (cand._system === 'detailed') || _isAtlasInDetailedMode();
      return _extractDiploidClass(cand, useDetailed);
    });
  }

  function _resolveCandidateById(candId) {
    const s = getAtlasState();
    if (!s) return null;
    if (s.candidate && s.candidate.id === candId) return s.candidate;
    if (Array.isArray(s.candidates)) {
      for (const c of s.candidates) if (c && c.id === candId) return c;
    }
    if (s.data && Array.isArray(s.data.candidates)) {
      for (const c of s.data.candidates) if (c && c.id === candId) return c;
    }
    return null;
  }

  function _isAtlasInDetailedMode() {
    const s = getAtlasState();
    if (!s) return false;
    // Atlas's label-vocabulary toggle. Best-effort detection across multiple
    // possible state field names (the toggle has gone through several names).
    if (s.labelVocab === 'detailed') return true;
    if (s.bandLabels === 'detailed') return true;
    if (s.detailedMode === true) return true;
    return false;
  }

  // ---------------------------------------------------------------------------
  // Extractor: per-candidate carrier / homozygosity predicates
  // ---------------------------------------------------------------------------

  function getCarriesHapMap(candId, hap) {
    return _cacheDim('carries_' + hap + '@' + candId, () => {
      const dc = getDiploidClassMap(candId);
      const out = new Map();
      for (const [id, label] of dc) out.set(id, carriesHaplotype(label, hap));
      return out;
    });
  }

  function getHomozygousHapMap(candId, hap) {
    return _cacheDim('homozygous_' + hap + '@' + candId, () => {
      const dc = getDiploidClassMap(candId);
      const out = new Map();
      for (const [id, label] of dc) out.set(id, isHomozygousFor(label, hap));
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: per-candidate fish-call quality dims
  // ---------------------------------------------------------------------------

  function getFishAmbiguousMap(candId) {
    return _cacheDim('fish_ambiguous@' + candId, () => {
      const cand = _resolveCandidateById(candId);
      if (!cand || !Array.isArray(cand.fish_calls)) return new Map();
      const out = new Map();
      for (const fc of cand.fish_calls) {
        if (!fc) continue;
        const id = fc.sample_id ? String(fc.sample_id) : sampleIdxToId(fc.sample_idx);
        out.set(id, !!fc.ambiguous);
      }
      return out;
    });
  }

  function getFishConfidenceMap(candId) {
    return _cacheDim('fish_confidence@' + candId, () => {
      const cand = _resolveCandidateById(candId);
      if (!cand || !Array.isArray(cand.fish_calls)) return new Map();
      const out = new Map();
      for (const fc of cand.fish_calls) {
        if (!fc) continue;
        const id = fc.sample_id ? String(fc.sample_id) : sampleIdxToId(fc.sample_idx);
        out.set(id, (typeof fc.confidence === 'number') ? fc.confidence : 0);
      }
      return out;
    });
  }

  // Subband-stability boolean: did the fish vote consistently across all L2s
  // inside the candidate (i.e. no parent crossover within the candidate span)?
  function getSubbandStableMap(candId) {
    return _cacheDim('subband_stable@' + candId, () => {
      const cand = _resolveCandidateById(candId);
      if (!cand || !Array.isArray(cand.fish_calls)) return new Map();
      const out = new Map();
      for (const fc of cand.fish_calls) {
        if (!fc) continue;
        const id = fc.sample_id ? String(fc.sample_id) : sampleIdxToId(fc.sample_idx);
        // Stable = n_supporting === n_intervals AND not ambiguous
        const stable = !fc.ambiguous &&
                       (fc.n_supporting === fc.n_intervals) &&
                       (fc.n_intervals > 0);
        out.set(id, !!stable);
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: family palette (global, precomp-derived)
  // ---------------------------------------------------------------------------

  function getFamilyMap() {
    return _cacheDim('family', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s || !s.data) return out;
      // Atlas family palette is loaded as an array indexed by sample_idx with
      // family ids as values. Multiple possible state locations.
      let arr = null;
      if (Array.isArray(s.data.family))            arr = s.data.family;
      else if (Array.isArray(s.data.family_id))    arr = s.data.family_id;
      else if (s.data.sample_meta && Array.isArray(s.data.sample_meta.family))
                                                   arr = s.data.sample_meta.family;
      if (!arr) return out;
      for (let i = 0; i < arr.length; i++) {
        if (arr[i] == null) continue;
        out.set(sampleIdxToId(i), 'family_' + String(arr[i]));
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: relatedness hub
  // ---------------------------------------------------------------------------

  function getRelatednessHubMap() {
    return _cacheDim('relatedness_hub', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s || !s.data) return out;
      let arr = null;
      if (Array.isArray(s.data.relatedness_hub))   arr = s.data.relatedness_hub;
      else if (Array.isArray(s.data.hub_id))       arr = s.data.hub_id;
      else if (s.data.sample_meta && Array.isArray(s.data.sample_meta.hub))
                                                   arr = s.data.sample_meta.hub;
      if (!arr) return out;
      for (let i = 0; i < arr.length; i++) {
        if (arr[i] == null || arr[i] === '' || arr[i] === -1) continue;
        out.set(sampleIdxToId(i), 'hub_' + String(arr[i]));
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: ROH carrier (boolean)
  // ---------------------------------------------------------------------------

  function getRohCarrierMap() {
    return _cacheDim('roh_carrier', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s || !s.data) return out;
      // F_ROH per-sample threshold. Atlas uses different field names; probe.
      let arr = null, threshold = 0.05;  // default
      if (Array.isArray(s.data.f_roh))             arr = s.data.f_roh;
      else if (Array.isArray(s.data.F_ROH))        arr = s.data.F_ROH;
      else if (s.data.sample_meta && Array.isArray(s.data.sample_meta.f_roh))
                                                   arr = s.data.sample_meta.f_roh;
      if (s.data.f_roh_threshold != null) threshold = s.data.f_roh_threshold;
      if (!arr) return out;
      for (let i = 0; i < arr.length; i++) {
        if (arr[i] == null || isNaN(arr[i])) continue;
        out.set(sampleIdxToId(i), arr[i] > threshold);
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: ancestry dominant Q (per-window OR per-sample-globally)
  // ---------------------------------------------------------------------------
  // For chat A we expose the GLOBAL dominant-Q (mean across all windows) as a
  // simple categorical dim. Per-window slicing is a cursor-region feature
  // handled later — for now, a sample's global dominant Q is sufficient for
  // most groupings ("compare K1-dominant to K2-dominant samples").

  function getAncestryDominantQMap() {
    return _cacheDim('ancestry_dominantQ', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s || !s.data || !s.data.ancestry_sample) return out;
      const as = s.data.ancestry_sample;
      // Long format reader. Wide format would need a sniffer like the server's.
      if (Array.isArray(as.samples) && Array.isArray(as.maxQ_label)) {
        // maxQ_label[w][s] — we want mode across windows per sample.
        const n_w = as.maxQ_label.length;
        const n_s = as.samples.length;
        for (let si = 0; si < n_s; si++) {
          const counts = Object.create(null);
          for (let wi = 0; wi < n_w; wi++) {
            const v = as.maxQ_label[wi] && as.maxQ_label[wi][si];
            if (v == null) continue;
            counts[v] = (counts[v] || 0) + 1;
          }
          let bestK = null, bestC = -1;
          for (const k in counts) if (counts[k] > bestC) { bestC = counts[k]; bestK = k; }
          if (bestK != null) out.set(String(as.samples[si]), 'K' + String(bestK));
        }
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: tracked samples (boolean)
  // ---------------------------------------------------------------------------

  function getTrackedMap() {
    return _cacheDim('tracked', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s) return out;
      const allIds = allSampleIds();
      // Default: false for all known samples
      for (const id of allIds) out.set(id, false);
      if (Array.isArray(s.tracked)) {
        for (const idx of s.tracked) out.set(sampleIdxToId(idx), true);
      }
      return out;
    });
  }

  // ---------------------------------------------------------------------------
  // Extractor: active sample set (boolean) — AS1 mask from turn 128
  // ---------------------------------------------------------------------------

  function getActiveSampleMap() {
    return _cacheDim('active_samples', () => {
      const s = getAtlasState();
      const out = new Map();
      if (!s) return out;
      const allIds = allSampleIds();
      if (s.activeSampleSet == null) {
        // null means "all active"
        for (const id of allIds) out.set(id, true);
      } else {
        for (const id of allIds) out.set(id, false);
        for (const idx of s.activeSampleSet) out.set(sampleIdxToId(idx), true);
      }
      return out;
    });
  }

  // ===========================================================================
  // Dimension dispatch — central lookup
  // ===========================================================================
  // Maps a dimension name (string) to its Map<sample_id, value>. Handles both
  // global names ('family') and per-candidate names ('diploid_class@<id>').

  function _lookupDimension(name) {
    if (typeof name !== 'string') return new Map();
    // Per-candidate dims have @ in the name
    const at = name.indexOf('@');
    if (at > 0) {
      const dimType = name.slice(0, at);
      const candId  = name.slice(at + 1);
      switch (dimType) {
        case 'diploid_class':       return getDiploidClassMap(candId);
        case 'carries_H1':          return getCarriesHapMap(candId, 'H1');
        case 'carries_H2':          return getCarriesHapMap(candId, 'H2');
        case 'carries_H3':          return getCarriesHapMap(candId, 'H3');
        case 'homozygous_H1':       return getHomozygousHapMap(candId, 'H1');
        case 'homozygous_H2':       return getHomozygousHapMap(candId, 'H2');
        case 'homozygous_H3':       return getHomozygousHapMap(candId, 'H3');
        case 'subband_stable':      return getSubbandStableMap(candId);
        case 'fish_ambiguous':      return getFishAmbiguousMap(candId);
        case 'fish_confidence':     return getFishConfidenceMap(candId);
        case 'fish_regime':         {
          // raw regime int (0..K-1 by PC1 rank, -1 = unassigned)
          return _cacheDim('fish_regime@' + candId, () => {
            const cand = _resolveCandidateById(candId);
            if (!cand || !Array.isArray(cand.fish_calls)) return new Map();
            const out = new Map();
            for (const fc of cand.fish_calls) {
              if (!fc) continue;
              const id = fc.sample_id ? String(fc.sample_id) : sampleIdxToId(fc.sample_idx);
              out.set(id, (typeof fc.regime === 'number') ? fc.regime : -1);
            }
            return out;
          });
        }
        default: return new Map();
      }
    }
    // Saved sets and lasso refs — pull from in-memory tables
    if (name.startsWith('saved_')) {
      const setName = name.slice(6);
      return _savedSetToMap(setName);
    }
    if (name.startsWith('lasso_')) {
      const lid = name.slice(6);
      return _lassoEntryToMap(lid);
    }
    // Global dims
    switch (name) {
      case 'family':                return getFamilyMap();
      case 'relatedness_hub':       return getRelatednessHubMap();
      case 'roh_carrier':           return getRohCarrierMap();
      case 'ancestry_dominantQ':    return getAncestryDominantQMap();
      case 'tracked':               return getTrackedMap();
      case 'active_samples':        return getActiveSampleMap();
      default:                      return new Map();
    }
  }

  // List every dimension name the engine currently knows about, so a UI layer
  // (turn 3) can populate dropdowns. Dynamic — derives from atlas state.
  function listDimensions() {
    const dims = [];
    // Globals (always present, even if data is empty)
    dims.push({ name: 'family',                kind: 'categorical' });
    dims.push({ name: 'relatedness_hub',       kind: 'categorical' });
    dims.push({ name: 'roh_carrier',           kind: 'boolean' });
    dims.push({ name: 'ancestry_dominantQ',    kind: 'categorical' });
    dims.push({ name: 'tracked',               kind: 'boolean' });
    dims.push({ name: 'active_samples',        kind: 'boolean' });
    // Per-candidate: enumerate candidates from atlas state
    const s = getAtlasState();
    const cands = [];
    if (s && s.candidate)             cands.push(s.candidate);
    if (s && Array.isArray(s.candidates)) for (const c of s.candidates) if (c && cands.indexOf(c) < 0) cands.push(c);
    if (s && s.data && Array.isArray(s.data.candidates))
      for (const c of s.data.candidates) if (c && cands.indexOf(c) < 0) cands.push(c);
    for (const cand of cands) {
      if (!cand || !cand.id) continue;
      const id = cand.id;
      dims.push({ name: 'diploid_class@'    + id, kind: 'categorical', candidate: id });
      dims.push({ name: 'fish_regime@'      + id, kind: 'integer',     candidate: id });
      dims.push({ name: 'fish_ambiguous@'   + id, kind: 'boolean',     candidate: id });
      dims.push({ name: 'fish_confidence@'  + id, kind: 'numeric',     candidate: id });
      dims.push({ name: 'subband_stable@'   + id, kind: 'boolean',     candidate: id });
      for (const h of ['H1', 'H2', 'H3']) {
        dims.push({ name: 'carries_'    + h + '@' + id, kind: 'boolean', candidate: id });
        dims.push({ name: 'homozygous_' + h + '@' + id, kind: 'boolean', candidate: id });
      }
    }
    // Saved sets
    for (const name of Object.keys(_savedSets)) {
      dims.push({ name: 'saved_' + name, kind: 'boolean', user: true });
    }
    // Lasso entries (recent, unpinned trimmed elsewhere)
    for (const entry of _lassoHistory) {
      dims.push({
        name: 'lasso_' + entry.id,
        kind: 'boolean',
        lasso: true,
        ts: entry.ts,
        n: entry.sample_ids.length,
      });
    }
    return dims;
  }

  // ===========================================================================
  // Expression AST + evaluator
  // ===========================================================================
  // AST nodes are plain JSON-compatible objects so they serialize trivially to
  // localStorage. Evaluator returns a Set of sample IDs (strings).

  // Validate an expression node before evaluation. Returns null if valid,
  // or a string error message if not. UI should call this before storing
  // user-typed expressions.
  function validateExpression(expr, depth) {
    depth = depth || 0;
    if (depth > 20) return 'expression too deep (max 20 levels)';
    if (!expr || typeof expr !== 'object') return 'expression must be an object';
    switch (expr.type) {
      case 'all':
      case 'none':
        return null;
      case 'eq':
      case 'neq':
        if (typeof expr.dim !== 'string') return 'eq/neq missing dim';
        if (expr.value == null) return 'eq/neq missing value';
        return null;
      case 'in':
      case 'not_in':
        if (typeof expr.dim !== 'string') return 'in/not_in missing dim';
        if (!Array.isArray(expr.values)) return 'in/not_in missing values array';
        return null;
      case 'tracked':
        if (typeof expr.set !== 'string') return 'tracked missing saved-set name';
        return null;
      case 'lasso':
        if (typeof expr.id !== 'string') return 'lasso missing id';
        return null;
      case 'gt':
      case 'gte':
      case 'lt':
      case 'lte':
        if (typeof expr.dim !== 'string') return 'comparison missing dim';
        if (typeof expr.value !== 'number') return 'comparison value must be number';
        return null;
      case 'truthy':
      case 'falsy':
        if (typeof expr.dim !== 'string') return 'truthy/falsy missing dim';
        return null;
      case 'and':
      case 'or':
        if (!Array.isArray(expr.children) || expr.children.length === 0)
          return 'and/or needs non-empty children';
        for (const c of expr.children) {
          const e = validateExpression(c, depth + 1);
          if (e) return e;
        }
        return null;
      case 'not':
        if (!expr.child) return 'not needs child';
        return validateExpression(expr.child, depth + 1);
      case 'minus':
        if (!expr.left || !expr.right) return 'minus needs left + right';
        const el = validateExpression(expr.left, depth + 1);
        if (el) return el;
        const er = validateExpression(expr.right, depth + 1);
        if (er) return er;
        return null;
      case 'ref':
        if (typeof expr.name !== 'string') return 'ref missing name';
        return null;
      default:
        return 'unknown expression type: ' + String(expr.type);
    }
  }

  // Evaluate to Set<sample_id>. _refStack guards against cycles.
  function _evaluate(expr, _refStack) {
    _refStack = _refStack || new Set();
    if (!expr || typeof expr !== 'object') return new Set();
    const allIds = allSampleIds();

    switch (expr.type) {
      case 'all':  return new Set(allIds);
      case 'none': return new Set();

      case 'eq': {
        const dm = _lookupDimension(expr.dim);
        const out = new Set();
        for (const [id, v] of dm) if (v === expr.value) out.add(id);
        return out;
      }
      case 'neq': {
        const dm = _lookupDimension(expr.dim);
        const out = new Set();
        // neq applies only to samples present in the dimension; samples
        // missing from the dim are NOT included (avoid silent cohort-wide
        // hits when the dim is empty).
        for (const [id, v] of dm) if (v !== expr.value) out.add(id);
        return out;
      }
      case 'in': {
        const dm = _lookupDimension(expr.dim);
        const valSet = new Set(expr.values);
        const out = new Set();
        for (const [id, v] of dm) if (valSet.has(v)) out.add(id);
        return out;
      }
      case 'not_in': {
        const dm = _lookupDimension(expr.dim);
        const valSet = new Set(expr.values);
        const out = new Set();
        for (const [id, v] of dm) if (!valSet.has(v)) out.add(id);
        return out;
      }
      case 'gt': case 'gte': case 'lt': case 'lte': {
        const dm = _lookupDimension(expr.dim);
        const out = new Set();
        for (const [id, v] of dm) {
          if (typeof v !== 'number' || isNaN(v)) continue;
          let pass = false;
          if (expr.type === 'gt')  pass = v >  expr.value;
          if (expr.type === 'gte') pass = v >= expr.value;
          if (expr.type === 'lt')  pass = v <  expr.value;
          if (expr.type === 'lte') pass = v <= expr.value;
          if (pass) out.add(id);
        }
        return out;
      }
      case 'truthy': {
        const dm = _lookupDimension(expr.dim);
        const out = new Set();
        for (const [id, v] of dm) if (v) out.add(id);
        return out;
      }
      case 'falsy': {
        const dm = _lookupDimension(expr.dim);
        const out = new Set();
        for (const [id, v] of dm) if (!v) out.add(id);
        return out;
      }
      case 'tracked': {
        return _savedSetToSet(expr.set);
      }
      case 'lasso': {
        const e = _findLassoById(expr.id);
        if (!e) return new Set();
        return new Set(e.sample_ids);
      }
      case 'and': {
        if (expr.children.length === 0) return new Set();
        let acc = _evaluate(expr.children[0], _refStack);
        for (let i = 1; i < expr.children.length; i++) {
          const next = _evaluate(expr.children[i], _refStack);
          const intersected = new Set();
          for (const id of acc) if (next.has(id)) intersected.add(id);
          acc = intersected;
        }
        return acc;
      }
      case 'or': {
        const out = new Set();
        for (const c of expr.children) {
          for (const id of _evaluate(c, _refStack)) out.add(id);
        }
        return out;
      }
      case 'not': {
        const inner = _evaluate(expr.child, _refStack);
        const out = new Set();
        for (const id of allIds) if (!inner.has(id)) out.add(id);
        return out;
      }
      case 'minus': {
        const l = _evaluate(expr.left, _refStack);
        const r = _evaluate(expr.right, _refStack);
        const out = new Set();
        for (const id of l) if (!r.has(id)) out.add(id);
        return out;
      }
      case 'ref': {
        if (_refStack.has(expr.name)) {
          // Cycle guard
          return new Set();
        }
        const e = _expressions[expr.name];
        if (!e) return new Set();
        const newStack = new Set(_refStack);
        newStack.add(expr.name);
        return _evaluate(e, newStack);
      }
      default:
        return new Set();
    }
  }

  // Evaluate against a default cohort scope. Returns an array of sorted IDs
  // (sorted for stable cache keys server-side).
  function evaluate(expr) {
    const err = validateExpression(expr);
    if (err) throw new Error('invalid expression: ' + err);
    const set = _evaluate(expr);
    return Array.from(set).sort();
  }

  // Count without materializing the full ID list — useful for live UI
  // counters next to expression slots.
  function countMembers(expr) {
    const err = validateExpression(expr);
    if (err) return 0;
    return _evaluate(expr).size;
  }

  // ===========================================================================
  // Saved expressions and sets
  // ===========================================================================
  // Two stores:
  //   _expressions[name] = AST node (composable, re-evaluates on dimension
  //                        invalidation)
  //   _savedSets[name]   = Set<sample_id_string> (frozen at save time, does
  //                        not re-evaluate)
  //
  // The user picks: "save as expression" (live) vs "save as set" (snapshot).

  let _expressions = {};
  let _savedSets   = {};

  function saveExpression(name, expr) {
    if (typeof name !== 'string' || !name) throw new Error('name required');
    const err = validateExpression(expr);
    if (err) throw new Error('cannot save invalid expression: ' + err);
    _expressions[name] = expr;
    _persistExpressions();
  }

  function deleteExpression(name) {
    delete _expressions[name];
    _persistExpressions();
  }

  function getExpression(name) { return _expressions[name] || null; }
  function listExpressions()   { return Object.keys(_expressions).slice(); }

  function saveSet(name, ids) {
    if (typeof name !== 'string' || !name) throw new Error('name required');
    if (!Array.isArray(ids)) throw new Error('ids must be an array');
    _savedSets[name] = ids.slice().map(String);
    _persistSavedSets();
  }

  function freezeSetFromExpression(name, expr) {
    saveSet(name, evaluate(expr));
  }

  function deleteSet(name) {
    delete _savedSets[name];
    _persistSavedSets();
  }

  function getSet(name) {
    return _savedSets[name] ? _savedSets[name].slice() : null;
  }

  function listSets() { return Object.keys(_savedSets).slice(); }

  function _savedSetToMap(name) {
    const out = new Map();
    const set = _savedSets[name];
    if (!set) return out;
    const setIds = new Set(set);
    for (const id of allSampleIds()) out.set(id, setIds.has(id));
    return out;
  }

  function _savedSetToSet(name) {
    const arr = _savedSets[name] || [];
    return new Set(arr);
  }

  // ===========================================================================
  // Slots — the "live" group set sent to the popstats server
  // ===========================================================================
  // Slots are positional named groups. Each slot points to an expression
  // (by-ref) or a saved set (by-ref). At Compute time the engine resolves
  // slots to ID lists.
  //
  // Slot record:
  //   { name: 'A', source: 'expression'|'set'|'inline'|'complement_of'
  //                       |'minus'|'union'|'intersect',
  //     ref: 'name' | inline_expr | { left, right } | { members: [name,...] },
  //     color: '#...' }
  //
  // Slot relations let the UI compose without nested trees. Slot C can be
  //   { source: 'minus', ref: { left: 'A', right: 'B' } }
  // and Slot E can be
  //   { source: 'union', ref: { members: ['C', 'D'] } }
  // Cycle-guarded by _seenSlots Set.

  let _slots = [];   // ordered, positional

  function setSlot(idx, slot) {
    if (idx < 0 || idx >= MAX_SLOTS) throw new Error('slot idx out of range');
    if (slot && typeof slot === 'object') {
      if (typeof slot.name !== 'string' || !/^[A-Za-z0-9_]+$/.test(slot.name)) {
        throw new Error('slot.name must match [A-Za-z0-9_]+ (engine F constraint)');
      }
    }
    _slots[idx] = slot;
    _persistSlots();
  }

  function clearSlot(idx) { _slots[idx] = null; _persistSlots(); }

  function getSlots() { return _slots.slice(0, MAX_SLOTS); }

  function _resolveSlotToIds(slot, _seenSlots) {
    _seenSlots = _seenSlots || new Set();
    if (!slot) return [];
    if (slot.name && _seenSlots.has(slot.name)) return [];   // global cycle guard
    const childSeen = slot.name ? new Set([..._seenSlots, slot.name]) : _seenSlots;
    if (slot.source === 'expression') {
      const expr = _expressions[slot.ref];
      if (!expr) return [];
      return evaluate(expr);
    }
    if (slot.source === 'set') {
      return getSet(slot.ref) || [];
    }
    if (slot.source === 'inline') {
      return evaluate(slot.ref);
    }
    if (slot.source === 'complement_of') {
      const target = _slots.find(s => s && s.name === slot.ref);
      if (!target) return [];
      const targetIds = new Set(_resolveSlotToIds(target, childSeen));
      const out = [];
      for (const id of allSampleIds()) if (!targetIds.has(id)) out.push(id);
      return out.sort();
    }
    if (slot.source === 'minus') {
      // ref = { left: slotName, right: slotName }
      const refObj = slot.ref || {};
      const leftSlot  = _slots.find(s => s && s.name === refObj.left);
      const rightSlot = _slots.find(s => s && s.name === refObj.right);
      if (!leftSlot)  return [];
      const leftIds  = _resolveSlotToIds(leftSlot, childSeen);
      const rightIds = rightSlot ? new Set(_resolveSlotToIds(rightSlot, childSeen)) : new Set();
      return leftIds.filter(id => !rightIds.has(id));
    }
    if (slot.source === 'union' || slot.source === 'intersect') {
      const refObj = slot.ref || {};
      const members = Array.isArray(refObj.members) ? refObj.members : [];
      if (members.length === 0) return [];
      const memberSlots = members.map(n => _slots.find(s => s && s.name === n)).filter(Boolean);
      if (memberSlots.length === 0) return [];
      const memberIdSets = memberSlots.map(s => new Set(_resolveSlotToIds(s, childSeen)));
      if (slot.source === 'union') {
        const out = new Set();
        for (const ms of memberIdSets) for (const id of ms) out.add(id);
        return Array.from(out).sort();
      } else {
        // intersect
        if (memberIdSets.length === 0) return [];
        let acc = memberIdSets[0];
        for (let i = 1; i < memberIdSets.length; i++) {
          const next = memberIdSets[i];
          const intersected = new Set();
          for (const id of acc) if (next.has(id)) intersected.add(id);
          acc = intersected;
        }
        return Array.from(acc).sort();
      }
    }
    return [];
  }

  // Compile the active slots into the {name: [ids]} shape ready for POST to
  // /api/popstats/groupwise. Drops empty slots and slots resolving to fewer
  // than `min_n` members (server-side filter is the authoritative gate, but
  // we reject obviously-too-small partitions client-side too).
  function compileGroupsForRequest(opts) {
    opts = opts || {};
    const min_n = (typeof opts.min_n === 'number') ? opts.min_n : 1;
    const groups = {};
    const warnings = [];
    for (const slot of _slots) {
      if (!slot) continue;
      const ids = _resolveSlotToIds(slot);
      if (ids.length < min_n) {
        warnings.push(
          'slot "' + slot.name + '" resolved to ' + ids.length +
          ' members (< min_n=' + min_n + '), dropped');
        continue;
      }
      if (groups[slot.name]) {
        warnings.push('duplicate slot name "' + slot.name + '" — using first');
        continue;
      }
      groups[slot.name] = ids;
    }
    return { groups, warnings };
  }

  // Resolve a single slot to its sample IDs. Used by the dock UI's slot grid
  // to display the live count next to each slot. Returns sorted IDs.
  function resolveSlotByIdx(idx) {
    if (idx < 0 || idx >= MAX_SLOTS) return [];
    return _resolveSlotToIds(_slots[idx]);
  }

  // Resolve any slot record (not necessarily one currently in _slots) — used
  // by the draft-slot scratchpad in the UI to show "if I save this, I'll get N".
  function resolveSlotRecord(slot) { return _resolveSlotToIds(slot); }

  // ===========================================================================
  // Cursor scope resolver
  // ===========================================================================
  // Returns {start_bp, end_bp} for a scope, or null for whole-chromosome.
  // Scopes:
  //   '1w'        ± 0 windows around state.cur (single window)
  //   '5w'        ± 2 windows
  //   '10w'       ± 5 windows
  //   'L2'        the L2 envelope containing state.cur (if any)
  //   'candidate' state.candidate's bp span
  //   'chrom'     null (whole chromosome)
  //   null        null (whole chromosome)

  function cursorRegion(scope) {
    const s = getAtlasState();
    if (!s || !s.data) return null;
    if (scope == null || scope === 'chrom') return null;
    if (scope === 'candidate') {
      const c = s.candidate;
      if (!c) return null;
      // start/end may be in mb or bp depending on candidate source.
      let start_bp = null, end_bp = null;
      if (typeof c.start_bp === 'number' && typeof c.end_bp === 'number') {
        start_bp = c.start_bp; end_bp = c.end_bp;
      } else if (typeof c.start_mb === 'number' && typeof c.end_mb === 'number') {
        start_bp = Math.round(c.start_mb * 1e6); end_bp = Math.round(c.end_mb * 1e6);
      }
      if (start_bp == null || end_bp == null) return null;
      return { start_bp, end_bp };
    }
    if (scope === 'L2') {
      if (!Array.isArray(s.data.l2_envelopes) || !s.windowToL2) return null;
      const li = s.windowToL2[s.cur];
      if (li < 0) return null;
      const env = s.data.l2_envelopes[li];
      if (!env) return null;
      const wins = s.data.windows;
      if (!wins) return null;
      const sw = wins[env._s0], ew = wins[env._e0];
      if (!sw || !ew) return null;
      const sBp = (typeof sw.start_bp === 'number') ? sw.start_bp
                : (typeof sw.center_mb === 'number') ? Math.round(sw.center_mb * 1e6) : null;
      const eBp = (typeof ew.end_bp === 'number') ? ew.end_bp
                : (typeof ew.center_mb === 'number') ? Math.round(ew.center_mb * 1e6) : null;
      if (sBp == null || eBp == null) return null;
      return { start_bp: sBp, end_bp: eBp };
    }
    // Window-count scopes: '1w', '5w', '10w'
    const m = /^(\d+)w$/.exec(scope);
    if (m) {
      const n = parseInt(m[1], 10);
      const half = Math.floor((n - 1) / 2);
      const wins = s.data.windows;
      if (!Array.isArray(wins) || wins.length === 0) return null;
      const cur = (typeof s.cur === 'number') ? s.cur : 0;
      const lo = Math.max(0, cur - half);
      const hi = Math.min(wins.length - 1, cur + half);
      const sw = wins[lo], ew = wins[hi];
      if (!sw || !ew) return null;
      const sBp = (typeof sw.start_bp === 'number') ? sw.start_bp
                : (typeof sw.center_mb === 'number') ? Math.round(sw.center_mb * 1e6) : null;
      const eBp = (typeof ew.end_bp === 'number') ? ew.end_bp
                : (typeof ew.center_mb === 'number') ? Math.round(ew.center_mb * 1e6) : null;
      if (sBp == null || eBp == null) return null;
      return { start_bp: sBp, end_bp: eBp };
    }
    return null;
  }

  function getCurrentChrom() {
    const s = getAtlasState();
    if (!s || !s.data) return null;
    return s.data.chrom || null;
  }

  // Build a fully-formed POST body for /api/popstats/groupwise. Slots become
  // groups; cursor scope becomes region. The atlas's request layer (turn 3)
  // takes this and ships it.
  function buildPopstatsRequest(opts) {
    opts = opts || {};
    const cohort_min_n = (typeof opts.min_n === 'number') ? opts.min_n : 3;
    const compiled = compileGroupsForRequest({ min_n: cohort_min_n });
    return {
      chrom:    opts.chrom    || getCurrentChrom(),
      region:   opts.region != null ? opts.region : cursorRegion(opts.scope || 'candidate'),
      groups:   compiled.groups,
      metrics:  opts.metrics  || ['fst', 'dxy', 'theta_pi'],
      win_bp:   opts.win_bp   || 50000,
      step_bp:  opts.step_bp  || 10000,
      win_type: (typeof opts.win_type === 'number') ? opts.win_type : 2,
      downsample: (typeof opts.downsample === 'number') ? opts.downsample : 1,
      ncores:   opts.ncores   || 4,
      _warnings: compiled.warnings,
    };
  }

  // ===========================================================================
  // Lasso / selection history
  // ===========================================================================

  let _lassoHistory = [];

  function recordSelection(rec) {
    if (!rec || !Array.isArray(rec.sample_ids)) {
      throw new Error('recordSelection needs {sample_ids: []}');
    }
    const e = {
      id:               rec.id || ('lasso_' + Date.now().toString(36) + '_' +
                                   Math.random().toString(36).slice(2, 6)),
      ts:               rec.ts || Date.now(),
      source_page:      rec.source_page || 'unknown',
      kind:             rec.kind || 'lasso',
      sample_ids:       rec.sample_ids.map(String),
      cursor_chrom:     rec.cursor_chrom || getCurrentChrom(),
      cursor_window:    (typeof rec.cursor_window === 'number') ? rec.cursor_window
                          : (function () { const s = getAtlasState(); return s ? s.cur : null; })(),
      cursor_candidate: rec.cursor_candidate ||
                          (function () { const s = getAtlasState();
                                         return (s && s.candidate) ? s.candidate.id : null; })(),
      note:             rec.note || '',
      pinned:           !!rec.pinned,
    };
    _lassoHistory.unshift(e);
    // Trim to max, keeping pinned
    const pinned = _lassoHistory.filter(x => x.pinned);
    const unpinned = _lassoHistory.filter(x => !x.pinned);
    const trimmed = unpinned.slice(0, Math.max(0, LASSO_HISTORY_MAX - pinned.length));
    _lassoHistory = pinned.concat(trimmed).sort((a, b) => b.ts - a.ts);
    _persistLassoHistory();
    return e.id;
  }

  function pinSelection(id, pinned) {
    const e = _lassoHistory.find(x => x.id === id);
    if (!e) return false;
    e.pinned = (pinned !== false);
    _persistLassoHistory();
    return true;
  }

  function dropSelection(id) {
    const before = _lassoHistory.length;
    _lassoHistory = _lassoHistory.filter(x => x.id !== id);
    if (_lassoHistory.length !== before) _persistLassoHistory();
    return _lassoHistory.length !== before;
  }

  function listSelections() { return _lassoHistory.slice(); }

  function _findLassoById(id) { return _lassoHistory.find(x => x.id === id) || null; }

  function _lassoEntryToMap(id) {
    const e = _findLassoById(id);
    const out = new Map();
    if (!e) return out;
    const setIds = new Set(e.sample_ids);
    for (const sid of allSampleIds()) out.set(sid, setIds.has(sid));
    return out;
  }

  // Materialize a lasso into a saved set with a user-given name.
  function lassoToSet(lassoId, savedName) {
    const e = _findLassoById(lassoId);
    if (!e) throw new Error('no lasso with id ' + lassoId);
    saveSet(savedName, e.sample_ids);
    return savedName;
  }

  // ===========================================================================
  // Snapshots — capture full state of interest for "come back later"
  // ===========================================================================

  let _snapshots = [];

  function saveSnapshot(opts) {
    opts = opts || {};
    const s = getAtlasState();
    const snap = {
      id:           opts.id || ('snap_' + Date.now().toString(36) + '_' +
                                Math.random().toString(36).slice(2, 6)),
      ts:           Date.now(),
      note:         opts.note || '',
      cohort:       getCohortKey(),
      chrom:        getCurrentChrom(),
      cursor:       (s && typeof s.cur === 'number') ? s.cur : null,
      cursor_scope: opts.cursor_scope || null,
      candidate:    (s && s.candidate) ? {
        id: s.candidate.id,
        start_mb: s.candidate.start_mb,
        end_mb:   s.candidate.end_mb,
      } : null,
      slots:        getSlots(),
      expressions:  Object.assign({}, _expressions),
      saved_sets:   Object.assign({}, _savedSets),
      last_result:  opts.last_result || null,
      pinned:       !!opts.pinned,
    };
    _snapshots.unshift(snap);
    const pinned   = _snapshots.filter(x => x.pinned);
    const unpinned = _snapshots.filter(x => !x.pinned);
    const trimmed  = unpinned.slice(0, Math.max(0, SNAPSHOT_MAX - pinned.length));
    _snapshots = pinned.concat(trimmed).sort((a, b) => b.ts - a.ts);
    _persistSnapshots();
    return snap.id;
  }

  function pinSnapshot(id, pinned) {
    const e = _snapshots.find(x => x.id === id);
    if (!e) return false;
    e.pinned = (pinned !== false);
    _persistSnapshots();
    return true;
  }

  function dropSnapshot(id) {
    const before = _snapshots.length;
    _snapshots = _snapshots.filter(x => x.id !== id);
    if (_snapshots.length !== before) _persistSnapshots();
    return _snapshots.length !== before;
  }

  function listSnapshots() { return _snapshots.slice(); }

  // restoreSnapshot returns the full record. It does NOT mutate atlas state
  // (jumping to a chrom/cursor/candidate is the atlas-side UI's job in turn 3).
  // It DOES restore expressions, saved_sets, and slots so that the Compute
  // button works against the same partition the snapshot used.
  function restoreSnapshot(id) {
    const snap = _snapshots.find(x => x.id === id);
    if (!snap) throw new Error('no snapshot with id ' + id);
    if (snap.expressions) _expressions = Object.assign({}, snap.expressions);
    if (snap.saved_sets)  _savedSets   = Object.assign({}, snap.saved_sets);
    if (Array.isArray(snap.slots)) _slots = snap.slots.slice();
    _persistExpressions();
    _persistSavedSets();
    _persistSlots();
    _bumpDimCache();
    return snap;
  }

  // ===========================================================================
  // Persistence (localStorage) — key per cohort, scoped per-chrom where useful
  // ===========================================================================

  function _safeGet(key) {
    if (typeof localStorage === 'undefined') return null;
    try { return localStorage.getItem(key); } catch (_) { return null; }
  }
  function _safeSet(key, val) {
    if (typeof localStorage === 'undefined') return;
    try { localStorage.setItem(key, val); } catch (_) {}
  }

  function _persistExpressions() {
    _safeSet(POPGEN_LS_EXPRESSIONS + '.' + getCohortKey(),
             JSON.stringify(_expressions));
  }
  function _persistSavedSets() {
    _safeSet(POPGEN_LS_SAVED_SETS + '.' + getCohortKey(),
             JSON.stringify(_savedSets));
  }
  function _persistSnapshots() {
    _safeSet(POPGEN_LS_SNAPSHOTS + '.' + getCohortKey(),
             JSON.stringify(_snapshots));
  }
  function _persistLassoHistory() {
    _safeSet(POPGEN_LS_LASSO_HIST + '.' + getCohortKey(),
             JSON.stringify(_lassoHistory));
  }
  function _persistSlots() {
    _safeSet(POPGEN_LS_SLOTS + '.' + getCohortKey(),
             JSON.stringify(_slots));
  }

  function loadFromLocalStorage() {
    const ck = getCohortKey();
    try {
      const r1 = _safeGet(POPGEN_LS_EXPRESSIONS + '.' + ck);
      if (r1) _expressions = JSON.parse(r1) || {};
    } catch (_) { _expressions = {}; }
    try {
      const r2 = _safeGet(POPGEN_LS_SAVED_SETS + '.' + ck);
      if (r2) _savedSets = JSON.parse(r2) || {};
    } catch (_) { _savedSets = {}; }
    try {
      const r3 = _safeGet(POPGEN_LS_SNAPSHOTS + '.' + ck);
      if (r3) _snapshots = JSON.parse(r3) || [];
    } catch (_) { _snapshots = []; }
    try {
      const r4 = _safeGet(POPGEN_LS_LASSO_HIST + '.' + ck);
      if (r4) _lassoHistory = JSON.parse(r4) || [];
    } catch (_) { _lassoHistory = []; }
    try {
      const r5 = _safeGet(POPGEN_LS_SLOTS + '.' + ck);
      if (r5) _slots = JSON.parse(r5) || [];
    } catch (_) { _slots = []; }
  }

  // ===========================================================================
  // Demo / self-test
  // ===========================================================================
  // Console-runnable demo. Builds the Quentin-example expression
  // (g3 ∪ g2) MINUS hub_30 ∪ g6_no_crossover, prints the compiled POST body.

  function demo() {
    const lines = [];
    lines.push('=== popgen demo ===');
    const s = getAtlasState();
    if (!s || !s.data) {
      lines.push('NO ATLAS STATE — running structural smoke tests only.');
    } else {
      lines.push('chrom: ' + getCurrentChrom() +
                 '   n_samples: ' + (s.data.n_samples || allSampleIds().length));
      const cId = (s.candidate && s.candidate.id) ? s.candidate.id : '<none>';
      lines.push('focal candidate: ' + cId);
      lines.push('detailed mode: ' + _isAtlasInDetailedMode());
    }
    lines.push('');
    lines.push('--- known dimensions ---');
    const dims = listDimensions();
    for (const d of dims) lines.push('  ' + d.name + ' (' + d.kind + ')');
    lines.push('');
    lines.push('--- example expression: H1/H2 ∩ NOT family_3 ---');
    const expr = {
      type: 'and',
      children: [
        // requires a focal candidate — use the first one we know about
        (function () {
          const c = s && s.candidate;
          if (!c || !c.id) return { type: 'all' };
          return { type: 'eq', dim: 'diploid_class@' + c.id, value: 'H1/H2' };
        })(),
        { type: 'not', child: {
          type: 'eq', dim: 'family', value: 'family_3'
        } },
      ],
    };
    const err = validateExpression(expr);
    lines.push('  valid: ' + (err ? 'NO — ' + err : 'yes'));
    if (!err) {
      const ids = evaluate(expr);
      lines.push('  resolved n: ' + ids.length);
      lines.push('  first 5: ' + ids.slice(0, 5).join(', '));
    }
    lines.push('');
    lines.push('--- compileGroupsForRequest with empty slots ---');
    const compiled = compileGroupsForRequest({ min_n: 1 });
    lines.push('  groups: ' + JSON.stringify(Object.keys(compiled.groups)));
    lines.push('  warnings: ' + compiled.warnings.length);
    lines.push('');
    lines.push('--- buildPopstatsRequest with cursor scope=candidate ---');
    const req = buildPopstatsRequest({ scope: 'candidate' });
    lines.push('  chrom: ' + req.chrom);
    lines.push('  region: ' + JSON.stringify(req.region));
    lines.push('  win_bp: ' + req.win_bp + '/' + req.step_bp);
    lines.push('');
    lines.push('=== demo done ===');
    const text = lines.join('\n');
    if (typeof console !== 'undefined') console.log(text);
    return text;
  }

  // ===========================================================================
  // Exposed API
  // ===========================================================================

  const api = {
    // Constants
    MAX_SLOTS,
    CURSOR_SCOPES,
    BAND_LABELS_DETAILED,
    BAND_LABELS_LEGACY,

    // Cohort + sample helpers
    getCohortKey, allSampleIds, sampleIdxToId, sampleIdToIdx,

    // Dimensions
    listDimensions, invalidateDimensions,
    parseHaplotypePair, carriesHaplotype, isHomozygousFor,

    // Expressions
    validateExpression, evaluate, countMembers,
    saveExpression, deleteExpression, getExpression, listExpressions,

    // Saved sets
    saveSet, freezeSetFromExpression, deleteSet, getSet, listSets,

    // Slots
    setSlot, clearSlot, getSlots, compileGroupsForRequest,
    resolveSlotByIdx, resolveSlotRecord,

    // Cursor + request
    cursorRegion, getCurrentChrom, buildPopstatsRequest,

    // Lasso history
    recordSelection, pinSelection, dropSelection, listSelections, lassoToSet,

    // Snapshots
    saveSnapshot, pinSnapshot, dropSnapshot, listSnapshots, restoreSnapshot,

    // Persistence
    loadFromLocalStorage,

    // Demo
    demo,

    // Internal exposure for testing
    _internals: {
      lookupDimension: _lookupDimension,
      evaluate: _evaluate,
      bumpDimCache: _bumpDimCache,
      get expressions() { return _expressions; },
      get savedSets()   { return _savedSets; },
      get lassoHistory() { return _lassoHistory; },
      get snapshots()   { return _snapshots; },
      get slots()       { return _slots; },
    },
  };

  // Auto-load persisted state on first import (no-op if no localStorage).
  loadFromLocalStorage();

  return api;
}));
