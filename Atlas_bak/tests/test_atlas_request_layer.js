// =============================================================================
// test_atlas_request_layer.js
// =============================================================================
// Tests for the request layer. Mocks fetch + indexedDB so the file is
// browser-shape-equivalent but runs in plain Node.
//
// Run:
//   node test_atlas_request_layer.js
// =============================================================================

'use strict';

// ----- Performance shim -------------------------------------------------------
if (typeof globalThis.performance === 'undefined') {
  globalThis.performance = { now: () => Date.now() };
}

// ----- IndexedDB shim (in-memory) --------------------------------------------
// Just enough of the IDB API for the request layer's needs: open(), object
// store with keyPath: 'cache_key', get/put/delete, openCursor.

function makeFakeIDB() {
  const store = new Map();   // key -> record

  function makeReq() {
    return {
      _resolve(v) { this.result = v; if (this.onsuccess) this.onsuccess({ target: this }); },
      _reject(e)  { this.error = e;  if (this.onerror)   this.onerror({ target: this }); },
    };
  }

  return {
    open(name, version) {
      const req = makeReq();
      const db = {
        objectStoreNames: { contains: (n) => n === 'live_results' },
        transaction(_storeName, _mode) {
          return {
            objectStore(_n) {
              return {
                get(key) {
                  const r = makeReq();
                  setTimeout(() => r._resolve(store.has(key) ? store.get(key) : undefined), 0);
                  return r;
                },
                put(record) {
                  const r = makeReq();
                  setTimeout(() => { store.set(record.cache_key, record); r._resolve(); }, 0);
                  return r;
                },
                delete(key) {
                  const r = makeReq();
                  setTimeout(() => { store.delete(key); r._resolve(); }, 0);
                  return r;
                },
                openCursor() {
                  const r = makeReq();
                  setTimeout(() => {
                    const entries = [...store.entries()];
                    let i = 0;
                    function step() {
                      if (i >= entries.length) { r._resolve(null); return; }
                      const [k, v] = entries[i++];
                      const cursor = {
                        key: k,
                        value: v,
                        delete() { store.delete(k); },
                        continue() { setTimeout(step, 0); },
                      };
                      r._resolve(cursor);
                    }
                    step();
                  }, 0);
                  return r;
                },
                createIndex(_n, _kp, _opts) {},
              };
            },
          };
        },
      };
      setTimeout(() => req._resolve(db), 0);
      return req;
    },
  };
}

globalThis.indexedDB = makeFakeIDB();

// ----- fetch mock with deterministic responses --------------------------------
// Records all calls; returns whatever the caller queues. Supports Abort.

const fetchCalls = [];
let fetchQueue = [];   // { delay_ms, ok, status, body } per call

function queueFetch(response) {
  fetchQueue.push(response);
}

globalThis.fetch = function (url, opts) {
  const call = { url, opts: opts || {}, ts: Date.now() };
  fetchCalls.push(call);
  const queued = fetchQueue.shift() || { ok: true, status: 200, body: { ok: true } };
  const delay = queued.delay_ms || 0;
  return new Promise((resolve, reject) => {
    let aborted = false;
    if (opts && opts.signal) {
      const sig = opts.signal;
      if (sig.aborted) { aborted = true; }
      else if (typeof sig.addEventListener === 'function') {
        sig.addEventListener('abort', () => { aborted = true; });
      }
    }
    setTimeout(() => {
      if (aborted) {
        const err = new Error('aborted'); err.name = 'AbortError';
        return reject(err);
      }
      resolve({
        ok: queued.ok,
        status: queued.status,
        statusText: queued.statusText || '',
        json: () => Promise.resolve(queued.body),
        text: () => Promise.resolve(typeof queued.body === 'string' ? queued.body
                                  : JSON.stringify(queued.body)),
        headers: { get() { return null; } },
      });
    }, delay);
  });
};

// ----- crypto.subtle shim using Node's crypto --------------------------------
if (typeof globalThis.crypto === 'undefined' || !globalThis.crypto.subtle) {
  const { createHash } = require('crypto');
  globalThis.crypto = {
    subtle: {
      digest(alg, data) {
        return Promise.resolve(
          createHash(alg.replace('-', '').toLowerCase())
            .update(Buffer.from(data))
            .digest()
            .buffer
        );
      },
    },
  };
}

// ----- Import + boot ---------------------------------------------------------
const factory = require('./atlas_request_layer.js');
// Build a stub server adapter
const server = {
  url: 'http://localhost:8765',
  status: 'unknown',
  isAvailable: async () => true,
};
const live = factory(server);

// ----- Test runner ----------------------------------------------------------
let _pass = 0, _fail = 0;
function ok(label, cond) {
  if (cond) { console.log('  ok  ' + label); _pass++; }
  else      { console.log('FAIL ' + label); _fail++; }
}
function eq(label, a, b) {
  const A = JSON.stringify(a), B = JSON.stringify(b);
  if (A === B) { console.log('  ok  ' + label); _pass++; }
  else         { console.log('FAIL ' + label + ': ' + A + ' != ' + B); _fail++; }
}
async function group(name, fn) { console.log(name); await fn(); }

const sleep = (ms) => new Promise(r => setTimeout(r, ms));

// =============================================================================
// Tests
// =============================================================================

(async () => {

  // ---------------------------------------------------------------------------
  await group('cache-key derivation parity', async () => {
    // Same body, group reordering and member reordering: same key.
    queueFetch({ ok: true, status: 200, body: {
      ok: true, engines: { region_popstats: 'enghashABC', hobs_windower: 'h', angsd_patched: 'a' } }
    });
    queueFetch({ ok: true, status: 200, body: {
      ok: true, engines: { region_popstats: 'enghashABC', hobs_windower: 'h', angsd_patched: 'a' } }
    });

    const reqA = {
      chrom: 'C_gar_LG28',
      region: { start_bp: 15_000_000, end_bp: 18_000_000 },
      groups: { HOM1: ['s1','s2','s3'], HOM2: ['s8','s9','s10'] },
      metrics: ['fst','dxy','theta_pi'],
      win_bp: 50000, step_bp: 10000, win_type: 2, downsample: 1,
    };
    const reqB = {
      chrom: 'C_gar_LG28',
      region: { start_bp: 15_000_000, end_bp: 18_000_000 },
      groups: { HOM2: ['s10','s9','s8'], HOM1: ['s3','s2','s1'] },
      metrics: ['theta_pi','dxy','fst'],
      win_bp: 50000, step_bp: 10000, win_type: 2, downsample: 1,
    };
    const k1 = await live.deriveCacheKey('popstats_groupwise', reqA);
    const k2 = await live.deriveCacheKey('popstats_groupwise', reqB);
    eq('same content, reordered → same key', k1, k2);
    ok('key is 32 hex chars', /^[0-9a-f]{32}$/.test(k1));

    // Different region → different key
    const reqC = Object.assign({}, reqA, {
      region: { start_bp: 16_000_000, end_bp: 19_000_000 },
    });
    const k3 = await live.deriveCacheKey('popstats_groupwise', reqC);
    ok('different region → different key', k3 !== k1);

    // Cross-language parity: hash this fixed payload and compare to a
    // pre-computed server-side hash. If atlas-side and server-side ever
    // diverge in canonicalization, this test catches it.
    //
    // Fixed expected key was obtained by running:
    //   popstats_server.popstats_cache_key(
    //     chrom='C_gar_LG28',
    //     region={'start_bp':15000000,'end_bp':18000000},
    //     groups={'HOM1':['s1','s2','s3'],'HOM2':['s8','s9','s10']},
    //     metrics=['fst','dxy','theta_pi'],
    //     win_bp=50000, step_bp=10000, win_type=2, downsample=1,
    //     engines={'region_popstats':'enghashABC'},
    //   )
    // == 'c2046458dd73d703995acb1f002b709a'
    eq('cross-language parity with server', k1, 'c2046458dd73d703995acb1f002b709a');
  });

  // ---------------------------------------------------------------------------
  await group('cache miss → IDB persist → cache hit', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    // The /api/health probe (refresh) — and we already refreshed once above
    // but the freshness check makes us hit it again only after HEALTH_REFRESH_MS.
    // We skip queueing health because the cached engine_hashes are still fresh.
    queueFetch({ ok: true, status: 200, delay_ms: 5, body: {
      kind: 'popstats_groupwise.v1',
      windows: [{ window_id: 'w1', chrom: 'LG28', theta_pi: 0.01 }],
      n_windows: 1,
    }});

    const req = {
      chrom: 'C_gar_LG28', region: null,
      groups: { HOM1: ['CGA_001','CGA_002','CGA_003'] },
      metrics: ['theta_pi'],
      win_bp: 50000, step_bp: 10000,
    };
    const r1 = await live.popstatsGroupwise(req, { debounce_ms: 0 });
    eq('first call cacheState=miss', r1.cacheState, 'miss');
    ok('payload returned',     r1.payload && r1.payload.n_windows === 1);
    eq('one fetch occurred',   fetchCalls.length, 1);

    // Second identical call → should hit IDB, no network
    fetchCalls.length = 0;
    const r2 = await live.popstatsGroupwise(req, { debounce_ms: 0 });
    eq('second call cacheState=hit', r2.cacheState, 'hit');
    eq('no fetch on hit',            fetchCalls.length, 0);
    eq('cache key matches',          r1.cache_key, r2.cache_key);
  });

  // ---------------------------------------------------------------------------
  await group('error envelope on 4xx', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    queueFetch({ ok: false, status: 400, body: {
      detail: "group 'HOM1' has n=1 < min_group_n=3"
    }});
    const r = await live.popstatsGroupwise({
      chrom: 'C_gar_LG28', region: null,
      groups: { HOM1: ['CGA_001'] },
    }, { debounce_ms: 0 });
    ok('envelope.ok = false', r.ok === false);
    eq('http_status forwarded', r.http_status, 400);
    ok('detail in error message', /min_group_n/.test(r.error));
  });

  // ---------------------------------------------------------------------------
  await group('in-flight de-dup', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    // Slow response so two concurrent requests overlap
    queueFetch({ ok: true, status: 200, delay_ms: 100, body: {
      kind: 'popstats_groupwise.v1',
      windows: [{ window_id: 'wDD', theta_pi: 0.02 }],
      n_windows: 1,
    }});

    const req = {
      chrom: 'C_gar_LG28', region: null,
      groups: { GROUPDD: ['CGA_004','CGA_005','CGA_006'] },
      metrics: ['theta_pi'],
      win_bp: 50000, step_bp: 10000,
    };
    // Fire two identical requests in parallel; the second should de-dup
    // onto the first's promise.
    const [a, b] = await Promise.all([
      live.popstatsGroupwise(req, { debounce_ms: 0 }),
      live.popstatsGroupwise(req, { debounce_ms: 0 }),
    ]);
    eq('only one fetch despite two concurrent calls', fetchCalls.length, 1);
    eq('both got the same cache_key', a.cache_key, b.cache_key);
    eq('both cacheState=miss',        a.cacheState, b.cacheState);
  });

  // ---------------------------------------------------------------------------
  await group('debounce: rapid calls collapse to one', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    queueFetch({ ok: true, status: 200, body: {
      kind: 'popstats_groupwise.v1',
      windows: [{ window_id: 'wDB', theta_pi: 0.03 }],
      n_windows: 1,
    }});

    const baseReq = {
      chrom: 'C_gar_LG28', region: null,
      groups: { GROUPX: ['CGA_007','CGA_008','CGA_009'] },
      metrics: ['theta_pi'],
      win_bp: 50000, step_bp: 10000,
    };
    // Fire 5 calls with same body in quick succession, 50ms debounce
    const promises = [];
    for (let i = 0; i < 5; i++) {
      // First 4 should be preempted by the 5th; await each individually
      promises.push(
        live.popstatsGroupwise(baseReq, { debounce_ms: 30 })
            .then(r => ({ ok: true, r }))
            .catch(e => ({ ok: false, e }))
      );
      await sleep(5);  // tiny gap so they actually queue separately
    }
    await sleep(80);  // let the final debounce fire
    const results = await Promise.all(promises);
    const fired = results.filter(x => x.ok && x.r && x.r.payload).length;
    const preempted = results.filter(x => !x.ok || (x.r && x.r.cacheState === undefined)).length;
    ok('exactly one fetch fired', fetchCalls.length === 1);
    ok('one promise resolved with payload', fired === 1);
    ok('four were preempted',     preempted >= 4 || results.length - fired === 4);
  });

  // ---------------------------------------------------------------------------
  await group('different chroms = different debounce channels (no preempt)', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    queueFetch({ ok: true, status: 200, body: {
      kind: 'popstats_groupwise.v1',
      windows: [], n_windows: 0,
    }});
    queueFetch({ ok: true, status: 200, body: {
      kind: 'popstats_groupwise.v1',
      windows: [], n_windows: 0,
    }});

    const reqLG28 = {
      chrom: 'C_gar_LG28', region: null,
      groups: { G: ['CGA_010','CGA_011','CGA_012'] },
      metrics: ['theta_pi'], win_bp: 50000, step_bp: 10000,
    };
    const reqLG14 = Object.assign({}, reqLG28, { chrom: 'C_gar_LG14' });
    const [a, b] = await Promise.all([
      live.popstatsGroupwise(reqLG28, { debounce_ms: 30 }),
      live.popstatsGroupwise(reqLG14, { debounce_ms: 30 }),
    ]);
    ok('two fetches (different channels)',  fetchCalls.length === 2);
    ok('both succeeded',                    a.ok && b.ok);
    ok('different cache keys',              a.cache_key !== b.cache_key);
  });

  // ---------------------------------------------------------------------------
  await group('cancel cancels in-flight', async () => {
    fetchCalls.length = 0; fetchQueue.length = 0;
    queueFetch({ ok: true, status: 200, delay_ms: 200, body: {
      kind: 'popstats_groupwise.v1', windows: [], n_windows: 0,
    }});
    const req = {
      chrom: 'C_gar_LG28', region: null,
      groups: { GG: ['CGA_001','CGA_002','CGA_003'] },
      metrics: ['theta_pi'], win_bp: 50000, step_bp: 10000,
    };
    const p = live.popstatsGroupwise(req, { debounce_ms: 0 });
    await sleep(20);  // let the request fire
    live.cancel();    // cancel everything
    const r = await p;
    ok('result is error envelope', !r.ok);
    ok('error mentions abort or cancel', /abort|cancel/i.test(r.error));
  });

  // ---------------------------------------------------------------------------
  await group('cacheStats shape', async () => {
    const stats = await live.cacheStats();
    ok('n_entries number', typeof stats.n_entries === 'number');
    ok('by_kind has popstats', stats.by_kind.popstats_groupwise > 0);
    ok('engine_hashes carried', !!stats.engine_hashes);
    ok('in_flight integer',     typeof stats.in_flight === 'number');
  });

  // ---------------------------------------------------------------------------
  await group('clearCache by prefix', async () => {
    const before = await live.cacheStats();
    const cleared = await live.clearCache('');
    const after = await live.cacheStats();
    ok('clearCache returned a number', typeof cleared === 'number');
    eq('after clear, n_entries=0', after.n_entries, 0);
    ok('cleared >= before',       cleared >= before.n_entries);
  });

  // ---------------------------------------------------------------------------
  console.log('');
  console.log('=========================================');
  console.log(_pass + ' passed, ' + _fail + ' failed');
  console.log('=========================================');
  if (_fail > 0) process.exit(1);

})().catch(e => {
  console.error('UNCAUGHT', e);
  process.exit(1);
});
