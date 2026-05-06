// =============================================================================
// atlas_request_layer.js — Turn 3 of chat A
// =============================================================================
//
// Atlas-side request layer for live popstats. Sits on top of the existing
// `atlasServer` adapter (line 43471 of Inversion_atlas.html) and adds:
//
//   1. atlasServer.popstatsGroupwise(req)   — wrap POST /api/popstats/groupwise
//   2. atlasServer.hobsGroupwise(req)       — wrap POST /api/popstats/hobs_groupwise
//   3. atlasServer.ancestryGroupwiseQ(req)  — wrap POST /api/ancestry/groupwise_q
//   4. atlasServer.popstatsHealth()         — refresh of the existing health
//                                              check + capture of engine hashes
//
//   5. window.popgenLive — top-level request orchestrator. Methods:
//        .request(kind, body, opts)
//          → debounce + cancel-stale + IndexedDB cache + returns Promise
//             with a result envelope { cacheState, payload, request_ms, ... }
//        .cancel(token)
//        .clearCache(prefix?)
//        .cacheStats()
//
// All public methods follow the same envelope shape so consumers don't have
// to special-case kind. The atlas's existing static path (drag-drop precomp
// JSONs) is untouched; live mode is additive.
//
// This file does NOT modify the atlasServer global directly when loaded in
// Node tests — instead it returns an installer fn that mutates the passed
// adapter. In-browser, popgenLive auto-installs against window.atlasServer.
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenLive = factory()(root.atlasServer || _makeStubServer());
  }

  function _makeStubServer() {
    // When the atlas adapter isn't loaded yet (script ordering), provide a
    // minimal stub so popgenLive doesn't blow up. The real adapter takes
    // over once it loads — popgenLive re-attaches via .attach(server).
    return {
      url: 'http://localhost:8765',
      status: 'unknown',
      isAvailable: async () => false,
      compute: async () => ({ ok: false, status: 0, error: 'no adapter' }),
    };
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const ENDPOINTS = {
    popstats_groupwise:    '/api/popstats/groupwise',
    hobs_groupwise:        '/api/popstats/hobs_groupwise',
    ancestry_groupwise_q:  '/api/ancestry/groupwise_q',
    health:                '/api/health',
    cache_keys:            '/api/cache/keys',
  };

  // Default debounce window (ms). User-configurable via popgenLive.config.
  const DEFAULT_DEBOUNCE_MS = 250;

  // Per-request timeout. Q07b's first run can take 30+ s; we set a long ceiling.
  const DEFAULT_TIMEOUT_MS  = 120_000;

  // IndexedDB config
  const IDB_NAME       = 'inversion_atlas_popgen';
  const IDB_VERSION    = 1;
  const IDB_STORE      = 'live_results';

  // LRU bound on the in-memory cache index. The IDB store holds the actual
  // payloads; this just bounds how many entries we track for fast probes.
  const MEM_INDEX_MAX  = 1000;

  // ===========================================================================
  // Cache key derivation — must match the server's cache_key logic byte-exact
  // so the atlas's "cached?" probe matches the server's "we have this" probe
  // and IndexedDB hits skip the network entirely.
  // ===========================================================================

  function _canonicalGroups(groups) {
    if (!groups) return [];
    return Object.keys(groups).sort().map(name => {
      const members = (groups[name] || []).slice().map(String).sort();
      return [name, members];
    });
  }

  // Browser-safe sha256 → 32 hex chars (matches server's truncated key length)
  async function _sha256_32(payload) {
    const enc = new TextEncoder().encode(payload);
    if (typeof crypto !== 'undefined' && crypto.subtle && crypto.subtle.digest) {
      const buf = await crypto.subtle.digest('SHA-256', enc);
      const hex = Array.from(new Uint8Array(buf))
        .map(b => b.toString(16).padStart(2, '0')).join('');
      return hex.slice(0, 32);
    }
    // Node fallback for tests
    if (typeof require !== 'undefined') {
      const { createHash } = require('crypto');
      return createHash('sha256').update(payload).digest('hex').slice(0, 32);
    }
    throw new Error('no crypto.subtle and no node crypto module');
  }

  function _stableStringify(obj) {
    // Deterministic JSON: keys sorted at every nesting level. Matches server's
    // json.dumps(..., sort_keys=True, separators=(',', ':')).
    if (obj === null || typeof obj !== 'object') return JSON.stringify(obj);
    if (Array.isArray(obj)) {
      return '[' + obj.map(_stableStringify).join(',') + ']';
    }
    const keys = Object.keys(obj).sort();
    const parts = keys.map(k => JSON.stringify(k) + ':' + _stableStringify(obj[k]));
    return '{' + parts.join(',') + '}';
  }

  async function popstatsCacheKey(req, engineHash) {
    const payload = _stableStringify({
      kind: 'popstats_groupwise.v1',
      chrom: req.chrom,
      region: req.region || null,
      groups: _canonicalGroups(req.groups),
      metrics: (req.metrics || ['fst','dxy','theta_pi']).slice().sort(),
      win_bp: req.win_bp || 50000,
      step_bp: req.step_bp || 10000,
      type: (req.win_type != null ? req.win_type : 2),
      downsample: (req.downsample != null ? req.downsample : 1),
      engine_region_popstats: engineHash || '?',
    });
    return _sha256_32(payload);
  }

  async function hobsCacheKey(req, engineHashHobs, engineHashAngsd) {
    const payload = _stableStringify({
      kind: 'hobs_groupwise.v1',
      chrom: req.chrom,
      region: req.region || null,
      groups: _canonicalGroups(req.groups),
      scales: (req.scales || ['10kb']).slice().sort(),
      engine_hobs_windower: engineHashHobs || '?',
      engine_angsd: engineHashAngsd || '?',
    });
    return _sha256_32(payload);
  }

  async function ancestryCacheKey(req) {
    const payload = _stableStringify({
      kind: 'ancestry_groupwise_q.v1',
      chrom: req.chrom,
      region: req.region || null,
      groups: _canonicalGroups(req.groups),
      K: req.K || 8,
      scale: req.scale || 'dense',
    });
    return _sha256_32(payload);
  }

  // ===========================================================================
  // IndexedDB wrapper — small, dependency-free, handles upgrade/version
  // ===========================================================================

  let _idb = null;
  let _idbOpenPromise = null;
  let _idbAvailable = (typeof indexedDB !== 'undefined');

  function _openIDB() {
    if (!_idbAvailable) return Promise.resolve(null);
    if (_idb) return Promise.resolve(_idb);
    if (_idbOpenPromise) return _idbOpenPromise;
    _idbOpenPromise = new Promise((resolve) => {
      const req = indexedDB.open(IDB_NAME, IDB_VERSION);
      req.onupgradeneeded = (ev) => {
        const db = ev.target.result;
        if (!db.objectStoreNames.contains(IDB_STORE)) {
          const store = db.createObjectStore(IDB_STORE, { keyPath: 'cache_key' });
          store.createIndex('kind', 'kind', { unique: false });
          store.createIndex('ts', 'ts', { unique: false });
        }
      };
      req.onsuccess = () => { _idb = req.result; resolve(_idb); };
      req.onerror = () => {
        // Disable IDB silently if access blocked (private browsing, quota)
        _idbAvailable = false;
        resolve(null);
      };
    });
    return _idbOpenPromise;
  }

  async function _idbGet(cacheKey) {
    const db = await _openIDB();
    if (!db) return null;
    return new Promise((resolve) => {
      let req;
      try {
        const tx = db.transaction(IDB_STORE, 'readonly');
        req = tx.objectStore(IDB_STORE).get(cacheKey);
      } catch (_) { return resolve(null); }
      req.onsuccess = () => resolve(req.result || null);
      req.onerror   = () => resolve(null);
    });
  }

  async function _idbPut(record) {
    const db = await _openIDB();
    if (!db) return false;
    return new Promise((resolve) => {
      try {
        const tx = db.transaction(IDB_STORE, 'readwrite');
        const req = tx.objectStore(IDB_STORE).put(record);
        req.onsuccess = () => resolve(true);
        req.onerror   = () => resolve(false);
      } catch (_) { resolve(false); }
    });
  }

  async function _idbDelete(cacheKey) {
    const db = await _openIDB();
    if (!db) return false;
    return new Promise((resolve) => {
      try {
        const tx = db.transaction(IDB_STORE, 'readwrite');
        const req = tx.objectStore(IDB_STORE).delete(cacheKey);
        req.onsuccess = () => resolve(true);
        req.onerror   = () => resolve(false);
      } catch (_) { resolve(false); }
    });
  }

  async function _idbList(prefix) {
    const db = await _openIDB();
    if (!db) return [];
    return new Promise((resolve) => {
      const out = [];
      try {
        const tx = db.transaction(IDB_STORE, 'readonly');
        const store = tx.objectStore(IDB_STORE);
        const req = store.openCursor();
        req.onsuccess = (ev) => {
          const cursor = ev.target.result;
          if (!cursor) return resolve(out);
          if (!prefix || String(cursor.key).startsWith(prefix)) {
            const v = cursor.value;
            out.push({
              cache_key: v.cache_key,
              kind: v.kind,
              ts: v.ts,
              chrom: v.chrom || null,
              size: (v.payload ? JSON.stringify(v.payload).length : 0),
            });
          }
          cursor.continue();
        };
        req.onerror = () => resolve(out);
      } catch (_) { resolve(out); }
    });
  }

  async function _idbClearByPrefix(prefix) {
    const db = await _openIDB();
    if (!db) return 0;
    return new Promise((resolve) => {
      let n = 0;
      try {
        const tx = db.transaction(IDB_STORE, 'readwrite');
        const store = tx.objectStore(IDB_STORE);
        const req = store.openCursor();
        req.onsuccess = (ev) => {
          const cursor = ev.target.result;
          if (!cursor) return resolve(n);
          if (!prefix || String(cursor.key).startsWith(prefix)) {
            cursor.delete(); n++;
          }
          cursor.continue();
        };
        req.onerror = () => resolve(n);
      } catch (_) { resolve(n); }
    });
  }

  // ===========================================================================
  // In-memory request registry — for cancel + de-dup of in-flight identical
  // requests. Two simultaneous Compute clicks with identical body produce
  // a single network roundtrip with both Promises resolving to the shared
  // result.
  // ===========================================================================

  const _inflight = new Map();   // cache_key -> { promise, controller, ts }
  const _memIndex = new Map();   // cache_key -> { kind, ts } (LRU touch on hit)

  function _bumpMemIndex(cacheKey, kind) {
    if (_memIndex.has(cacheKey)) _memIndex.delete(cacheKey);
    _memIndex.set(cacheKey, { kind, ts: Date.now() });
    while (_memIndex.size > MEM_INDEX_MAX) {
      const oldest = _memIndex.keys().next().value;
      _memIndex.delete(oldest);
    }
  }

  // ===========================================================================
  // Debouncer — per-channel. Each kind ('popstats', 'hobs', etc.) has its
  // own debounce queue so a fast hobs request doesn't get blocked by a slow
  // popstats request.
  // ===========================================================================

  const _debouncers = new Map();  // channel -> { timer, latest_call }

  function _debounce(channel, ms, fn) {
    return new Promise((resolve, reject) => {
      const existing = _debouncers.get(channel);
      if (existing && existing.timer) {
        clearTimeout(existing.timer);
        // Reject the previous waiter — that call's debounce window was preempted
        if (existing.reject) existing.reject({ cancelled: 'preempted' });
      }
      const rec = {
        timer: null,
        resolve, reject,
      };
      _debouncers.set(channel, rec);
      rec.timer = setTimeout(async () => {
        try {
          const result = await fn();
          if (rec === _debouncers.get(channel)) {
            _debouncers.delete(channel);
            rec.resolve(result);
          } else {
            // someone preempted us mid-flight; result still resolves, but to
            // a no-op since the rejection already fired
          }
        } catch (e) {
          if (rec === _debouncers.get(channel)) {
            _debouncers.delete(channel);
            rec.reject(e);
          }
        }
      }, ms);
    });
  }

  // ===========================================================================
  // Engine hash registry — cached from /api/health, used to mint cache keys
  // that are byte-equivalent to the server's keys.
  // ===========================================================================

  let _engineHashes = {
    region_popstats: null,
    hobs_windower: null,
    angsd_patched: null,
    instant_q: null,
    last_refresh: 0,
  };
  const HEALTH_REFRESH_MS = 30_000;

  async function _refreshEngineHashes(server) {
    if (!server || !server.url) return;
    const stale = (Date.now() - _engineHashes.last_refresh) > HEALTH_REFRESH_MS;
    if (!stale && _engineHashes.region_popstats) return;
    try {
      const resp = await fetch(server.url + ENDPOINTS.health, { cache: 'no-store' });
      if (!resp || !resp.ok) return;
      const data = await resp.json();
      if (data && data.engines) {
        _engineHashes.region_popstats = data.engines.region_popstats || null;
        _engineHashes.hobs_windower   = data.engines.hobs_windower   || null;
        _engineHashes.angsd_patched   = data.engines.angsd_patched   || null;
        _engineHashes.instant_q       = data.engines.instant_q       || null;
        _engineHashes.last_refresh    = Date.now();
      }
    } catch (_) { /* server unreachable — leave hashes stale */ }
  }

  // ===========================================================================
  // Result envelope
  // ===========================================================================
  // Every public request returns one of these shapes so consumers can be
  // uniform across kinds.
  //
  //   { ok: true,
  //     cacheState: 'hit'|'miss'|'computing',  // 'computing' means in-flight
  //     cache_key: '...',
  //     kind: 'popstats_groupwise',
  //     payload: {...server response...},
  //     request_ms: 47,                        // total round-trip incl. cache
  //     fetched_ms: 30,                        // network only (0 on cache hit)
  //     ts: <epoch>,
  //   }
  //
  //   { ok: false, error: '...', http_status: 503 }

  function _envelopeOk(kind, cacheKey, payload, cacheState, fetched_ms, request_ms) {
    return {
      ok: true,
      kind,
      cache_key: cacheKey,
      cacheState,
      payload,
      request_ms,
      fetched_ms,
      ts: Date.now(),
    };
  }

  function _envelopeErr(kind, msg, status, cacheKey) {
    return {
      ok: false,
      kind,
      cache_key: cacheKey || null,
      error: String(msg),
      http_status: status || 0,
      ts: Date.now(),
    };
  }

  // ===========================================================================
  // Network POST helper with timeout + AbortController
  // ===========================================================================

  async function _postJSON(url, body, controller, timeoutMs) {
    const ctrl = controller || (typeof AbortController !== 'undefined'
                                ? new AbortController() : null);
    const timer = (typeof setTimeout !== 'undefined') ? setTimeout(() => {
      if (ctrl && typeof ctrl.abort === 'function') ctrl.abort();
    }, timeoutMs || DEFAULT_TIMEOUT_MS) : null;
    let resp;
    try {
      resp = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body),
        signal: ctrl ? ctrl.signal : undefined,
      });
    } finally {
      if (timer) clearTimeout(timer);
    }
    let bodyJson = null;
    let bodyText = null;
    try { bodyJson = await resp.json(); }
    catch (_) {
      try { bodyText = await resp.text(); } catch (__) {}
    }
    return { resp, bodyJson, bodyText };
  }

  // ===========================================================================
  // Public request orchestrator
  // ===========================================================================

  function makePopgenLive(server) {

    if (!server) {
      throw new Error('popgenLive requires a server adapter');
    }

    // Derive cache key for a given (kind, body). Awaits engine-hash refresh
    // because keys depend on engine hashes for popstats + hobs.
    async function deriveCacheKey(kind, body) {
      await _refreshEngineHashes(server);
      switch (kind) {
        case 'popstats_groupwise':
          return popstatsCacheKey(body, _engineHashes.region_popstats);
        case 'hobs_groupwise':
          return hobsCacheKey(body, _engineHashes.hobs_windower,
                                    _engineHashes.angsd_patched);
        case 'ancestry_groupwise_q':
          return ancestryCacheKey(body);
        default:
          throw new Error('unknown kind: ' + kind);
      }
    }

    // Probe the local cache without firing a network request. Returns the
    // cached envelope or null.
    async function probeCache(kind, body) {
      const cacheKey = await deriveCacheKey(kind, body);
      const rec = await _idbGet(cacheKey);
      if (!rec) return null;
      _bumpMemIndex(cacheKey, kind);
      return _envelopeOk(kind, cacheKey, rec.payload, 'hit', 0,
                         /*request_ms*/ 0);
    }

    // Fire the request. Internal — wrapped by .request() with debounce/cancel.
    async function _fireRequest(kind, body, opts) {
      const t0 = (typeof performance !== 'undefined') ? performance.now() : Date.now();
      const cacheKey = await deriveCacheKey(kind, body);

      // Cache hit?
      if (!opts || !opts.no_cache) {
        const cached = await _idbGet(cacheKey);
        if (cached) {
          _bumpMemIndex(cacheKey, kind);
          const t1 = (typeof performance !== 'undefined') ? performance.now() : Date.now();
          return _envelopeOk(kind, cacheKey, cached.payload, 'hit',
                             /*fetched_ms*/ 0, t1 - t0);
        }
      }

      // De-dup in-flight: if an identical request is already live, await it
      const existing = _inflight.get(cacheKey);
      if (existing && existing.promise) {
        return existing.promise;
      }

      const ctrl = (typeof AbortController !== 'undefined') ? new AbortController() : null;
      const url = server.url + ENDPOINTS[kind];
      const fetchedT0 = (typeof performance !== 'undefined') ? performance.now() : Date.now();

      const promise = (async () => {
        try {
          const { resp, bodyJson, bodyText } =
            await _postJSON(url, body, ctrl, opts && opts.timeout_ms);
          const fetchedT1 = (typeof performance !== 'undefined') ? performance.now() : Date.now();
          const fetched_ms = Math.round(fetchedT1 - fetchedT0);
          if (!resp.ok) {
            const detail = (bodyJson && bodyJson.detail) || bodyText || resp.statusText;
            return _envelopeErr(kind, detail, resp.status, cacheKey);
          }
          if (!bodyJson) {
            return _envelopeErr(kind, 'empty or non-JSON response', resp.status, cacheKey);
          }
          // Persist to IDB
          await _idbPut({
            cache_key: cacheKey,
            kind,
            chrom: body.chrom || null,
            ts: Date.now(),
            payload: bodyJson,
          });
          _bumpMemIndex(cacheKey, kind);
          const t1 = (typeof performance !== 'undefined') ? performance.now() : Date.now();
          return _envelopeOk(kind, cacheKey, bodyJson, 'miss', fetched_ms, Math.round(t1 - t0));
        } catch (e) {
          if (e && e.name === 'AbortError') {
            return _envelopeErr(kind, 'request aborted', 0, cacheKey);
          }
          return _envelopeErr(kind, (e && e.message) || String(e), 0, cacheKey);
        } finally {
          _inflight.delete(cacheKey);
        }
      })();

      _inflight.set(cacheKey, { promise, controller: ctrl, ts: Date.now(), kind });
      return promise;
    }

    // Cancel an in-flight request by cache key (or all of them).
    function cancel(cacheKey) {
      if (cacheKey == null) {
        // Cancel everything
        for (const [, rec] of _inflight) {
          if (rec.controller && typeof rec.controller.abort === 'function') {
            try { rec.controller.abort(); } catch (_) {}
          }
        }
        _inflight.clear();
        return;
      }
      const rec = _inflight.get(cacheKey);
      if (rec && rec.controller && typeof rec.controller.abort === 'function') {
        try { rec.controller.abort(); } catch (_) {}
        _inflight.delete(cacheKey);
      }
    }

    // Public request — debounce + cancel-stale + cache + dedup
    function request(kind, body, opts) {
      opts = opts || {};
      const debounceMs = (opts.debounce_ms != null) ? opts.debounce_ms
                                                    : DEFAULT_DEBOUNCE_MS;
      // Channel = kind + chrom; a popstats LG28 request and a popstats LG14
      // request belong to different channels (no preemption).
      const channel = kind + ':' + (body.chrom || '?');
      if (debounceMs > 0) {
        return _debounce(channel, debounceMs, () => _fireRequest(kind, body, opts));
      }
      return _fireRequest(kind, body, opts);
    }

    // Convenience wrappers — three named methods so callers don't pass kind
    function popstatsGroupwise(body, opts)   { return request('popstats_groupwise',   body, opts); }
    function hobsGroupwise(body, opts)       { return request('hobs_groupwise',       body, opts); }
    function ancestryGroupwiseQ(body, opts)  { return request('ancestry_groupwise_q', body, opts); }

    // Health refresh (also done lazily on first request)
    async function refreshHealth() {
      await _refreshEngineHashes(server);
      return Object.assign({}, _engineHashes);
    }

    // Cache management
    async function clearCache(prefix) {
      _memIndex.clear();
      return _idbClearByPrefix(prefix || '');
    }

    async function cacheStats() {
      const all = await _idbList('');
      const byKind = {};
      let totalBytes = 0;
      for (const e of all) {
        byKind[e.kind] = (byKind[e.kind] || 0) + 1;
        totalBytes += e.size;
      }
      return {
        n_entries: all.length,
        by_kind: byKind,
        bytes: totalBytes,
        in_flight: _inflight.size,
        engine_hashes: Object.assign({}, _engineHashes),
      };
    }

    // Server adapter swap (atlas's atlasServer adapter loads after this script
    // in some builds). Call attach() with the real adapter once it loads.
    function attach(realServer) { server = realServer; }

    // For consumers that want to install named methods onto the existing
    // atlasServer adapter directly (so other code can call e.g.
    // atlasServer.popstatsGroupwise(req)), provide installer:
    function installOnAdapter(adapter) {
      if (!adapter || typeof adapter !== 'object') return;
      adapter.popstatsGroupwise   = popstatsGroupwise;
      adapter.hobsGroupwise       = hobsGroupwise;
      adapter.ancestryGroupwiseQ  = ancestryGroupwiseQ;
      adapter.popstatsHealth      = refreshHealth;
    }

    return {
      // Top-level requests
      request,
      popstatsGroupwise,
      hobsGroupwise,
      ancestryGroupwiseQ,

      // Cache + introspection
      probeCache,
      deriveCacheKey,
      clearCache,
      cacheStats,
      cancel,

      // Server lifecycle
      attach,
      refreshHealth,
      installOnAdapter,

      // Test hooks
      _internals: {
        _idbGet, _idbPut, _idbDelete, _idbList,
        _engineHashes,
        get inflight() { return _inflight; },
        get memIndex() { return _memIndex; },
      },
    };
  }

  // The factory returns a function — caller passes in the server adapter.
  // In the browser the IIFE already invokes it with window.atlasServer.
  return function (server) { return makePopgenLive(server); };
}));
