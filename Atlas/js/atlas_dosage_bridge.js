// =============================================================================
// atlas_dosage_bridge.js — atlas-side installer for the dosage chunk bridge
// =============================================================================
//
// Connects the live `renderDosageHeatmap` path in Inversion_atlas.html to the
// popstats server's GET /api/dosage/chunk endpoint. The renderer reads chunks
// from `state.data.dosage_chunks.chunks[i].url`; this module synthesises that
// index from a probe to /api/dosage/manifest, with templated URLs that get
// substituted to the requested region at fetch time.
//
// Architecture summary
// --------------------
//   1. atlasInstallServerDosageBridge() probes /api/dosage/manifest via
//      atlasServer.url. If unreachable or empty, no-op.
//   2. For each chrom in the manifest, synthesise one chunk with bounds
//      [1, length_bp] (full chromosome) and URL pattern:
//         <atlasServer.url>/api/dosage/chunk?chrom=__CHROM__&start=__START__&end=__END__&cap=__CAP__
//      The placeholders are substituted in `_resolveChunkUrl` (already added
//      to Inversion_atlas.html as the second surgical edit for this slice).
//   3. Replace state.data.dosage_chunks with `{ chrom, chunks: [...], source: 'atlas_server' }`
//      where chrom is the active chrom (mirrors the pre-existing schema).
//   4. The renderer is none the wiser — it iterates chunks via _findCoveringChunk
//      and fetches whatever URL comes back, with substitution happening
//      transparently in _fetchAndCacheChunk.
//
// Spec: specs_todo/from_turn129/S6_dosage_heatmap_streaming_viewer.md (P4.1).
// Sister tests: tests/dosage_bridge/test_dosage_bridge.py (Python side),
//               tests/dosage_bridge/test_atlas_dosage_bridge.js (this file).
//
// Public API:
//   AtlasDosageBridge.install({ baseUrl, activeChrom, cap, force })
//   AtlasDosageBridge.uninstall()
//   AtlasDosageBridge.buildSyntheticIndex(manifest, { baseUrl, activeChrom, cap })
//   AtlasDosageBridge.getStatus()        -> { installed: bool, chroms: int, baseUrl, ts }
//
// Idempotent: install() can be called repeatedly; subsequent calls update
// the active chrom + refresh the manifest. uninstall() restores whatever was
// in `state.data.dosage_chunks` before install (or null if there was nothing).
// =============================================================================

(function (root, factory) {
  'use strict';
  if (typeof module === 'object' && module.exports) {
    module.exports = factory();
  } else {
    root.AtlasDosageBridge = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Config
  // ---------------------------------------------------------------------------
  const DEFAULT_BASE_URL    = 'http://localhost:8765';
  const DEFAULT_CAP_DEFAULT = 1000;
  const DEFAULT_TIMEOUT_MS  = 1500;
  const MANIFEST_PATH       = '/api/dosage/manifest';
  const CHUNK_PATH          = '/api/dosage/chunk';
  const SCHEMA_VERSION      = 'dosage_manifest_v1';

  // Saved state for uninstall + getStatus
  let _installed   = false;
  let _previousIdx = null;
  let _baseUrl     = DEFAULT_BASE_URL;
  let _ts          = 0;
  let _chromsCount = 0;
  let _activeChrom = null;

  // ---------------------------------------------------------------------------
  // Manifest probe — fetch /api/dosage/manifest with a short timeout
  // ---------------------------------------------------------------------------
  function _fetchManifest(baseUrl, timeoutMs) {
    if (typeof fetch === 'undefined') {
      return Promise.resolve(null);
    }
    const url = (baseUrl || DEFAULT_BASE_URL).replace(/\/+$/, '') + MANIFEST_PATH;
    let timer = null;
    const ctrl = (typeof AbortController === 'function') ? new AbortController() : null;
    if (ctrl) timer = setTimeout(() => ctrl.abort(), timeoutMs || DEFAULT_TIMEOUT_MS);
    return fetch(url, ctrl ? { signal: ctrl.signal, cache: 'no-store' } : { cache: 'no-store' })
      .then(r => {
        if (timer) clearTimeout(timer);
        if (!r || !r.ok) return null;
        return r.json();
      })
      .catch(() => null);
  }

  // ---------------------------------------------------------------------------
  // buildSyntheticIndex — pure transform from a manifest payload to the
  // chunk-index shape that `state.data.dosage_chunks` expects.
  //
  // The renderer's _findCoveringChunk picks the first chunk where
  // start_bp <= region.start_bp && end_bp >= region.end_bp. We give each
  // chrom one chunk spanning [1, length_bp], so any in-bounds region
  // resolves to it.
  //
  // The URL template uses placeholder tokens that _resolveChunkUrl in
  // Inversion_atlas.html substitutes at fetch time:
  //   __CHROM__ → chunkRef.chrom (URL-encoded)
  //   __START__ → region.start_bp
  //   __END__   → region.end_bp
  //   __CAP__   → chunkRef.cap or state.data.dosage_chunks.cap_default
  // ---------------------------------------------------------------------------
  function buildSyntheticIndex(manifest, opts) {
    opts = opts || {};
    if (!manifest || !Array.isArray(manifest.chroms) || manifest.chroms.length === 0) {
      return null;
    }
    const baseUrl  = (opts.baseUrl || DEFAULT_BASE_URL).replace(/\/+$/, '');
    const cap      = (typeof opts.cap === 'number' && opts.cap > 0) ? opts.cap : DEFAULT_CAP_DEFAULT;
    const activeCh = opts.activeChrom || (manifest.chroms[0] && manifest.chroms[0].name) || '';
    const urlTpl   = baseUrl + CHUNK_PATH +
      '?chrom=__CHROM__&start=__START__&end=__END__&cap=__CAP__';

    const chunks = manifest.chroms.map(c => {
      // Defensive: if length_bp is 0 (couldn't read), give a generous upper
      // bound so requests aren't silently rejected by the covering check.
      const len = (typeof c.length_bp === 'number' && c.length_bp > 0)
                    ? c.length_bp : 1e9;
      return {
        chrom:    c.name,
        start_bp: 1,
        end_bp:   len,
        url:      urlTpl,
        cap:      cap,
        // Mark as synthetic so callers can tell static and live chunks apart
        _synthetic: true,
      };
    });

    return {
      chrom:        activeCh,
      chunks:       chunks,
      source:       'atlas_server',
      base_url:     baseUrl,
      cap_default:  cap,
      schema_version: SCHEMA_VERSION,
      manifest_n_samples: manifest.n_samples || 0,
      manifest_limits: manifest.limits || null,
      ts:           Date.now(),
    };
  }

  // ---------------------------------------------------------------------------
  // install — fetch the manifest, build the index, write it onto state.data.
  // ---------------------------------------------------------------------------
  // Options:
  //   baseUrl     base URL of the popstats server (default 'http://localhost:8765')
  //   activeChrom which chrom the atlas is currently focused on; sets
  //               state.data.dosage_chunks.chrom (the renderer reads it
  //               for badging / status text). Defaults to the first chrom
  //               in the manifest if not given.
  //   cap         default cap to put on synthetic chunks. The renderer's
  //               selectTopMarkers does its own capping, but the cap value
  //               flows through __CAP__ when the renderer doesn't override.
  //   force       skip the cache-lifetime check; always re-probe.
  //   stateRef    optional: the global state object. Defaults to window.state.
  //   atlasServer optional: the global atlasServer adapter (for url + status).
  //   manifest    optional: skip the network probe and use this object.
  //               Used by tests to inject a fake manifest.
  //
  // Returns Promise<{installed: bool, reason?: string, chroms?: number}>.
  // ---------------------------------------------------------------------------
  async function install(opts) {
    opts = opts || {};
    const stateRef = opts.stateRef ||
      (typeof window !== 'undefined' && window.state) || null;
    if (!stateRef) {
      return { installed: false, reason: 'no_state' };
    }
    if (!stateRef.data) stateRef.data = {};

    let baseUrl = opts.baseUrl;
    if (!baseUrl) {
      const adapter = opts.atlasServer ||
        (typeof window !== 'undefined' && window.atlasServer) || null;
      baseUrl = (adapter && adapter.url) || DEFAULT_BASE_URL;
    }

    // Either use the provided manifest or fetch one
    const manifest = opts.manifest || await _fetchManifest(baseUrl, opts.timeoutMs);
    if (!manifest) {
      return { installed: false, reason: 'manifest_unavailable', baseUrl: baseUrl };
    }
    if (!Array.isArray(manifest.chroms) || manifest.chroms.length === 0) {
      return { installed: false, reason: 'manifest_empty', baseUrl: baseUrl };
    }
    if (manifest.schema_version &&
        manifest.schema_version !== SCHEMA_VERSION) {
      return { installed: false, reason: 'schema_mismatch',
               got: manifest.schema_version, expected: SCHEMA_VERSION };
    }

    const idx = buildSyntheticIndex(manifest, {
      baseUrl: baseUrl,
      activeChrom: opts.activeChrom ||
        (stateRef.data && (stateRef.data.chrom || stateRef.data.activeChrom)),
      cap: opts.cap,
    });
    if (!idx) {
      return { installed: false, reason: 'index_build_failed' };
    }

    // Save whatever was there before so uninstall can restore it. We
    // capture the prior index whenever the existing dosage_chunks isn't
    // one we made — i.e. its `source` is not 'atlas_server'. That way:
    //   - First install over a static index: snapshot it.
    //   - Re-install while already synthetic: don't overwrite (the
    //     intermediate synthetic isn't "original").
    //   - Install over no prior index: snapshot null.
    // This handles the chrom-switch flow (multiple installs in a row)
    // without losing the user's original static index.
    const existing = stateRef.data.dosage_chunks;
    const existingIsSynthetic = !!(existing && existing.source === 'atlas_server');
    if (!existingIsSynthetic) {
      _previousIdx = existing || null;
    }
    stateRef.data.dosage_chunks = idx;

    _installed   = true;
    _baseUrl     = baseUrl;
    _ts          = idx.ts;
    _chromsCount = idx.chunks.length;
    _activeChrom = idx.chrom;

    return {
      installed: true,
      chroms: _chromsCount,
      baseUrl: _baseUrl,
      activeChrom: _activeChrom,
    };
  }

  // ---------------------------------------------------------------------------
  // uninstall — restore the pre-install dosage_chunks (or null).
  // ---------------------------------------------------------------------------
  function uninstall(opts) {
    opts = opts || {};
    const stateRef = opts.stateRef ||
      (typeof window !== 'undefined' && window.state) || null;
    if (!stateRef || !stateRef.data) {
      _installed = false;
      _previousIdx = null;
      return { installed: false };
    }
    if (_installed) {
      stateRef.data.dosage_chunks = _previousIdx;
    }
    _installed   = false;
    _previousIdx = null;
    _ts          = 0;
    _chromsCount = 0;
    _activeChrom = null;
    return { installed: false };
  }

  // ---------------------------------------------------------------------------
  // setActiveChrom — update state.data.dosage_chunks.chrom on a chrom switch
  // without re-probing the manifest. Cheap; called when the atlas's "active
  // chromosome" changes (e.g. user navigates from LG28 to LG12).
  // ---------------------------------------------------------------------------
  function setActiveChrom(chrom, opts) {
    opts = opts || {};
    const stateRef = opts.stateRef ||
      (typeof window !== 'undefined' && window.state) || null;
    if (!stateRef || !stateRef.data || !stateRef.data.dosage_chunks) {
      return false;
    }
    const idx = stateRef.data.dosage_chunks;
    if (idx.source !== 'atlas_server') return false;   // don't touch static indices
    idx.chrom = chrom;
    _activeChrom = chrom;
    return true;
  }

  // ---------------------------------------------------------------------------
  // getStatus — diagnostics for the toolbar badge / debug logging
  // ---------------------------------------------------------------------------
  function getStatus() {
    return {
      installed:   _installed,
      chroms:      _chromsCount,
      baseUrl:     _baseUrl,
      activeChrom: _activeChrom,
      ts:          _ts,
    };
  }

  // ---------------------------------------------------------------------------
  return {
    install: install,
    uninstall: uninstall,
    setActiveChrom: setActiveChrom,
    getStatus: getStatus,
    buildSyntheticIndex: buildSyntheticIndex,
    _const: {
      DEFAULT_BASE_URL,
      DEFAULT_CAP_DEFAULT,
      MANIFEST_PATH,
      CHUNK_PATH,
      SCHEMA_VERSION,
    },
  };
}));
