// =============================================================================
// tests/dosage_bridge/test_html_surgical_edits.js
// =============================================================================
// Verifies the two atlas-side surgical edits made for the S6.P4.1 bridge:
//
//   1. _resolveChunkUrl(chunkRef, region) substitutes __CHROM__/__START__/
//      __END__/__CAP__ placeholders in templated URLs, AND returns static
//      URLs unchanged.
//
//   2. _fetchAndCacheChunk(chunkRef, region) builds a cache key that
//      INCLUDES the requested region for templated chunks (so different
//      regions don't collide on the same cache entry), but stays
//      backward-compatible (no region tag) for static chunks.
//
// We extract the new function bodies out of Inversion_atlas.html via regex
// so this test is faithful to what's in the file (no risk of drift between
// test and reality). Then we eval them in a sandbox with mocks for fetch/
// state/AbortController and exercise both code paths.
//
// Run:
//   cd Atlas && node tests/dosage_bridge/test_html_surgical_edits.js
// =============================================================================

'use strict';

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const HTML_PATH = path.resolve(__dirname, '../../Inversion_atlas.html');
const html = fs.readFileSync(HTML_PATH, 'utf-8');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  \u2713 ' + msg); }
  else      { fail++; console.error('  \u2717 ' + msg); }
}
function section(name) { console.log('\n\u2014 ' + name + ' \u2014'); }

// ---------------------------------------------------------------------------
section('Extract function bodies from Inversion_atlas.html');

function extractFn(name) {
  // Match the function definition through to its closing brace at column 0.
  // The functions we care about all start at indent level 0 in the inline
  // <script> block (i.e. the top of an outer function or IIFE).
  // Use a non-greedy match up through the next "^}\n" line.
  const re = new RegExp(
    '\\nfunction ' + name + '\\([^)]*\\)\\s*\\{[\\s\\S]*?\\n\\}\\n',
    'm'
  );
  const m = html.match(re);
  if (!m) return null;
  return m[0];
}

const findCoveringSrc   = extractFn('_findCoveringChunk');
const fetchAndCacheSrc  = extractFn('_fetchAndCacheChunk');
const resolveUrlSrc     = extractFn('_resolveChunkUrl');

assert(findCoveringSrc, '_findCoveringChunk source extracted');
assert(fetchAndCacheSrc, '_fetchAndCacheChunk source extracted');
assert(resolveUrlSrc, '_resolveChunkUrl source extracted');

assert(fetchAndCacheSrc.indexOf('region') !== -1,
       '_fetchAndCacheChunk now takes a region parameter');
assert(fetchAndCacheSrc.indexOf('_resolveChunkUrl') !== -1,
       '_fetchAndCacheChunk calls _resolveChunkUrl');
assert(findCoveringSrc.indexOf('_fetchAndCacheChunk(') !== -1 &&
       findCoveringSrc.indexOf(', region)') !== -1,
       '_findCoveringChunk passes region to _fetchAndCacheChunk');

// ---------------------------------------------------------------------------
section('_resolveChunkUrl in isolation');

// Build a sandbox with just the helpers we need.
function buildSandbox(opts) {
  opts = opts || {};
  const fetchCalls = [];
  const cache = new Map();
  const sandbox = {
    state: opts.state || { data: { chrom: 'C_gar_LG28' } },
    console: console,
    fetch: opts.fetch || function (url, _opts) {
      fetchCalls.push(url);
      return Promise.resolve({
        ok: true,
        json: () => Promise.resolve({ tag: 'fake', url: url }),
      });
    },
    AbortController: function () { this.abort = function () {}; this.signal = null; },
    setTimeout: setTimeout,
    clearTimeout: clearTimeout,
    Promise: Promise,
    // Helpers the source expects to find in scope
    _getDosageChunkCache: function () { return cache; },
    _ensureDosageHmState: function () { return { pending_fetch: null }; },
    _regionKey: function (region) {
      return `${region.start_bp}:${region.end_bp}`;
    },
  };
  vm.createContext(sandbox);
  // Inject the three functions
  vm.runInContext(resolveUrlSrc, sandbox);
  vm.runInContext(fetchAndCacheSrc, sandbox);
  vm.runInContext(findCoveringSrc, sandbox);
  return { sandbox, fetchCalls, cache };
}

const { sandbox } = buildSandbox();

// Case 1: static URL (no placeholders) returns unchanged
const staticUrl = sandbox._resolveChunkUrl(
  { url: 'static_chunk_LG28.json' },
  { start_bp: 17000000, end_bp: 17500000 }
);
assert(staticUrl === 'static_chunk_LG28.json',
       'static URL passes through unchanged');

// Case 2: templated URL with all 4 placeholders substituted
const tpl = 'http://h:8765/api/dosage/chunk?chrom=__CHROM__&start=__START__&end=__END__&cap=__CAP__';
const resolved = sandbox._resolveChunkUrl(
  { url: tpl, chrom: 'C_gar_LG28', cap: 1500 },
  { start_bp: 17000000, end_bp: 17500000 }
);
assert(resolved.indexOf('chrom=C_gar_LG28') !== -1,
       'template: __CHROM__ substituted from chunkRef.chrom');
assert(resolved.indexOf('start=17000000') !== -1,
       'template: __START__ substituted from region.start_bp');
assert(resolved.indexOf('end=17500000') !== -1,
       'template: __END__ substituted from region.end_bp');
assert(resolved.indexOf('cap=1500') !== -1,
       'template: __CAP__ substituted from chunkRef.cap');
assert(resolved.indexOf('__') === -1,
       'no unresolved placeholders left');

// Case 3: chrom from state.data when chunkRef.chrom is missing
const fromState = sandbox._resolveChunkUrl(
  { url: tpl, cap: 500 },
  { start_bp: 1, end_bp: 100 }
);
assert(fromState.indexOf('chrom=C_gar_LG28') !== -1,
       'template: missing chunkRef.chrom -> falls back to state.data.chrom');

// Case 4: nullish url returns nullish
assert(sandbox._resolveChunkUrl({ url: null }, {}) == null,
       'null url -> null');
assert(sandbox._resolveChunkUrl({}, {}) == null,
       'no url -> nullish');

// Case 5: chrom URL-encoded for safety
const encoded = sandbox._resolveChunkUrl(
  { url: tpl, chrom: 'C/Bad_LG12', cap: 100 },
  { start_bp: 1, end_bp: 100 }
);
assert(encoded.indexOf('chrom=C%2FBad_LG12') !== -1,
       'chrom value URL-encoded in resolved URL');

// ---------------------------------------------------------------------------
section('_fetchAndCacheChunk — cache key includes region for templated URLs');

const env2 = buildSandbox();
const tplChunk = {
  url: tpl, chrom: 'C_gar_LG28', start_bp: 1, end_bp: 41203400, cap: 1000,
};
// Two different regions on the same templated chunk should produce two
// distinct cache entries (otherwise a region-A query would return a cached
// region-B payload — wrong content).
const r1 = { start_bp: 17000000, end_bp: 17500000 };
const r2 = { start_bp: 18000000, end_bp: 18500000 };

return env2.sandbox._fetchAndCacheChunk(tplChunk, r1).then(a => {
  assert(a && a.tag === 'fake', 'first fetch returns a chunk');
  assert(env2.cache.size === 1, 'cache has 1 entry after first fetch');
  return env2.sandbox._fetchAndCacheChunk(tplChunk, r2);
}).then(b => {
  assert(b && b.tag === 'fake', 'second fetch returns a chunk');
  assert(env2.cache.size === 2,
         'cache has 2 entries after fetching a different region '
         + '(templated chunks key on region — no collision)');
  // Same region again -> cache hit, no new entry
  return env2.sandbox._fetchAndCacheChunk(tplChunk, r1);
}).then(c => {
  assert(c && c.tag === 'fake', 'third fetch (same region as first) returns chunk');
  assert(env2.cache.size === 2,
         'cache size unchanged on repeat region (cache hit)');
  // Confirm the URL actually substituted in fetch
  assert(env2.fetchCalls.length >= 2,
         '2 fetches issued (the third was a cache hit)');
  assert(env2.fetchCalls[0].indexOf('start=17000000') !== -1,
         'first fetch URL resolved with region 1');
  assert(env2.fetchCalls[1].indexOf('start=18000000') !== -1,
         'second fetch URL resolved with region 2');
}).then(() => {
  // ------------------------------------------------------------------------
  section('_fetchAndCacheChunk — static URL backward compat');

  const env3 = buildSandbox();
  const staticChunk = {
    url: 'static_LG28.json', chrom: 'C_gar_LG28', start_bp: 1, end_bp: 41203400,
  };
  return env3.sandbox._fetchAndCacheChunk(staticChunk,
                                            { start_bp: 1000, end_bp: 2000 })
    .then(a => {
      assert(a && a.tag === 'fake', 'static URL fetch returns chunk');
      // Static URLs must NOT include the region in the cache key, matching
      // pre-existing behavior. For backward compat the file URL is fetched
      // verbatim.
      assert(env3.fetchCalls[0] === 'static_LG28.json',
             'static URL fetched verbatim (no substitution)');
      // Repeat with a different region — for static URLs this must be a
      // cache hit because the URL is identical and the chunk content
      // doesn't depend on region.
      return env3.sandbox._fetchAndCacheChunk(staticChunk,
                                              { start_bp: 9999, end_bp: 99999 });
    }).then(b => {
      assert(env3.fetchCalls.length === 1,
             'static URL: second fetch with different region is a cache hit');
      assert(env3.cache.size === 1, 'static URL: cache still has 1 entry');
    });
}).then(() => {
  // ------------------------------------------------------------------------
  console.log('\n=========================================');
  console.log('  PASS: ' + pass + '   FAIL: ' + fail);
  console.log('=========================================');
  process.exit(fail === 0 ? 0 : 1);
}).catch(err => {
  console.error('uncaught:', err);
  process.exit(2);
});
