// Tests for turn 2o — atlas server adapter (atlasServer.* methods)
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
  url: 'http://localhost/',
});
const w = dom.window;

function run() {
  const failures = [];
  let testNum = 0;
  function t(name, fn) {
    testNum++;
    return Promise.resolve()
      .then(() => fn())
      .then(() => console.log(`  PASS [${testNum}] ${name}`))
      .catch(e => { failures.push({ name, err: e.message }); console.log(`  FAIL [${testNum}] ${name}: ${e.message}`); });
  }
  function eq(a, b, m) { if (a !== b) throw new Error(`${m||''} expected ${JSON.stringify(b)}, got ${JSON.stringify(a)}`); }

  // ============================================================
  // Adapter exists
  // ============================================================
  if (!w.atlasServer) {
    console.log('  FAIL: atlasServer not exposed');
    return Promise.resolve([{ name: 'fns exposed', err: 'missing' }]);
  }
  console.log('  atlasServer exposed');

  function makeFetchStub(handlers) {
    // handlers: map of urlPattern → function returning {ok, status, body, contentType}
    return async function stubFetch(url, opts) {
      for (const pattern of Object.keys(handlers)) {
        if (url.indexOf(pattern) >= 0) {
          const r = await handlers[pattern](url, opts);
          return {
            ok: r.ok !== false,
            status: r.status || (r.ok === false ? 500 : 200),
            headers: {
              get: (k) => k.toLowerCase() === 'content-type' ? (r.contentType || 'text/plain') : null,
            },
            text: async () => r.body || '',
            json: async () => (typeof r.body === 'string' ? JSON.parse(r.body) : r.body),
          };
        }
      }
      throw new Error(`no handler for ${url}`);
    };
  }

  function failingFetch() {
    return async () => { throw new Error('connection refused'); };
  }

  function pchain() {
    let p = Promise.resolve();
    return {
      add(name, fn) { p = p.then(() => t(name, fn)); return this; },
      done() { return p.then(() => failures); }
    };
  }

  return pchain()
    // ============================================================
    // Defaults + URL persistence
    // ============================================================
    .add('default url is localhost:8765', () => {
      eq(w.atlasServer.url, 'http://localhost:8765');
    })
    .add('initial status is unknown', () => {
      // Reset (other tests may have run before)
      w.atlasServer.status = 'unknown';
      w.atlasServer.lastChecked = 0;
      eq(w.atlasServer.status, 'unknown');
      eq(w.atlasServer.mode, 'view');   // unknown → view
    })
    .add('setUrl persists in localStorage and resets status', () => {
      w.atlasServer.status = 'available';
      w.atlasServer.setUrl('http://localhost:9999');
      eq(w.atlasServer.url, 'http://localhost:9999');
      eq(w.atlasServer.status, 'unknown');
      const saved = w.localStorage.getItem('inversion_atlas.serverUrl');
      eq(saved, 'http://localhost:9999');
      // Restore default
      w.atlasServer.setUrl('http://localhost:8765');
    })

    // ============================================================
    // isAvailable: caching
    // ============================================================
    .add('isAvailable: returns true when fetch succeeds', async () => {
      w.atlasServer.status = 'unknown';
      w.atlasServer.lastChecked = 0;
      let calls = 0;
      w.fetch = makeFetchStub({
        '/health': () => { calls++; return { ok: true, body: '{"ok":true}', contentType: 'application/json' }; },
      });
      const result = await w.atlasServer.isAvailable();
      eq(result, true);
      eq(w.atlasServer.status, 'available');
      eq(w.atlasServer.mode, 'server');
      eq(calls, 1);
    })
    .add('isAvailable: cached for 30s, no second fetch', async () => {
      // First call should hit fetch; immediate second call should skip it
      w.atlasServer.status = 'unknown';
      w.atlasServer.lastChecked = 0;
      let calls = 0;
      w.fetch = makeFetchStub({
        '/health': () => { calls++; return { ok: true, body: '{"ok":true}' }; },
      });
      await w.atlasServer.isAvailable();
      eq(calls, 1, 'first call should fetch');
      // Second call within cache window — should NOT increment
      await w.atlasServer.isAvailable();
      eq(calls, 1, 'second call should hit cache');
    })
    .add('isAvailable: forceRefresh skips cache', async () => {
      let calls = 0;
      w.fetch = makeFetchStub({
        '/health': () => { calls++; return { ok: true, body: '{"ok":true}' }; },
      });
      await w.atlasServer.isAvailable();
      const startCalls = calls;
      await w.atlasServer.isAvailable(true);
      if (calls === startCalls) throw new Error('forceRefresh did not refetch');
    })
    .add('isAvailable: returns false on failed fetch', async () => {
      w.atlasServer.status = 'unknown';
      w.atlasServer.lastChecked = 0;
      w.fetch = failingFetch();
      const result = await w.atlasServer.isAvailable();
      eq(result, false);
      eq(w.atlasServer.status, 'unavailable');
      eq(w.atlasServer.mode, 'view');
    })
    .add('isAvailable: returns false on non-200 response', async () => {
      w.atlasServer.status = 'unknown';
      w.atlasServer.lastChecked = 0;
      w.fetch = makeFetchStub({
        '/health': () => ({ ok: false, status: 500, body: 'broken' }),
      });
      const result = await w.atlasServer.isAvailable();
      eq(result, false);
    })

    // ============================================================
    // read: GET /file/<path>
    // ============================================================
    .add('read: success returns text', async () => {
      w.fetch = makeFetchStub({
        '/file/foo.txt': () => ({ ok: true, body: 'hello world', contentType: 'text/plain' }),
      });
      const r = await w.atlasServer.read('foo.txt');
      eq(r.ok, true);
      eq(r.text, 'hello world');
    })
    .add('read: JSON content-type returns parsed json', async () => {
      w.fetch = makeFetchStub({
        '/file/data.json': () => ({ ok: true, body: '{"x":42}', contentType: 'application/json' }),
      });
      const r = await w.atlasServer.read('data.json');
      eq(r.ok, true);
      eq(r.json.x, 42);
    })
    .add('read: opts.as=json forces JSON parsing', async () => {
      w.fetch = makeFetchStub({
        '/file/data': () => ({ ok: true, body: '{"y":7}', contentType: 'text/plain' }),
      });
      const r = await w.atlasServer.read('data', { as: 'json' });
      eq(r.ok, true);
      eq(r.json.y, 7);
    })
    .add('read: 404 returns ok=false', async () => {
      w.fetch = makeFetchStub({
        '/file/nope': () => ({ ok: false, status: 404, body: '' }),
      });
      const r = await w.atlasServer.read('nope');
      eq(r.ok, false);
      eq(r.status, 404);
    })
    .add('read: empty path returns error', async () => {
      const r = await w.atlasServer.read('');
      eq(r.ok, false);
    })
    .add('read: connection error returns ok=false with error', async () => {
      w.fetch = failingFetch();
      const r = await w.atlasServer.read('foo.txt');
      eq(r.ok, false);
      if (!r.error) throw new Error('expected error message');
    })

    // ============================================================
    // write: POST /file/<path>
    // ============================================================
    .add('write: serializes object as JSON', async () => {
      let captured = null;
      w.fetch = makeFetchStub({
        '/file/out.json': (url, opts) => { captured = opts; return { ok: true, body: '{"ok":true}' }; },
      });
      const r = await w.atlasServer.write('out.json', { hello: 'world' });
      eq(r.ok, true);
      eq(captured.method, 'POST');
      eq(captured.headers['Content-Type'], 'application/json');
      eq(captured.body, '{"hello":"world"}');
    })
    .add('write: passes string body verbatim', async () => {
      let captured = null;
      w.fetch = makeFetchStub({
        '/file/out.txt': (url, opts) => { captured = opts; return { ok: true, body: '{}' }; },
      });
      const r = await w.atlasServer.write('out.txt', 'plain text');
      eq(r.ok, true);
      eq(captured.body, 'plain text');
      eq(captured.headers['Content-Type'], 'text/plain');
    })
    .add('write: empty path returns error', async () => {
      const r = await w.atlasServer.write('', {});
      eq(r.ok, false);
    })

    // ============================================================
    // compute: POST /compute/<name>
    // ============================================================
    .add('compute: sends args as JSON, parses JSON response', async () => {
      let captured = null;
      w.fetch = makeFetchStub({
        '/compute/echo': (url, opts) => {
          captured = opts;
          return { ok: true, body: '{"ok":true,"result":{"echoed":"yes"}}', contentType: 'application/json' };
        },
      });
      const r = await w.atlasServer.compute('echo', { x: 1 });
      eq(r.ok, true);
      eq(JSON.parse(captured.body).x, 1);
      eq(r.json.ok, true);
      eq(r.json.result.echoed, 'yes');
    })
    .add('compute: missing name returns error', async () => {
      const r = await w.atlasServer.compute('', {});
      eq(r.ok, false);
    })
    .add('compute: server error returns ok=false', async () => {
      w.fetch = makeFetchStub({
        '/compute/broken': () => ({ ok: false, status: 500, body: 'oops' }),
      });
      const r = await w.atlasServer.compute('broken', {});
      eq(r.ok, false);
      eq(r.status, 500);
    })

    // ============================================================
    // Mode badge
    // ============================================================
    .add('badge: exists in DOM after init', () => {
      // The init function runs on DOMContentLoaded or via setTimeout
      // Force a re-init so we can check DOM
      try { w._atlasServerInitBadge(); } catch (_) {}
      const badge = w.document.getElementById('atlasServerBadge');
      if (!badge) throw new Error('badge not in DOM');
    })
    .add('badge: clickable element', () => {
      const badge = w.document.getElementById('atlasServerBadge');
      eq(badge.tagName.toLowerCase(), 'button');
    })
    .add('badge: only one badge created (idempotent init)', () => {
      try { w._atlasServerInitBadge(); } catch (_) {}
      try { w._atlasServerInitBadge(); } catch (_) {}
      const badges = w.document.querySelectorAll('#atlasServerBadge');
      eq(badges.length, 1);
    })

    .done();
}

setTimeout(() => {
  run().then(failures => {
    if (failures.length > 0) {
      console.log(`\n${failures.length} test(s) failed`);
      failures.forEach(f => console.log(`  - ${f.name}: ${f.err}`));
      process.exit(1);
    } else {
      console.log('\nAll tests passed');
      process.exit(0);
    }
  });
}, 500);
