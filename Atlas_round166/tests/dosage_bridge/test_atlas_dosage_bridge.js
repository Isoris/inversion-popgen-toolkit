// =============================================================================
// tests/dosage_bridge/test_atlas_dosage_bridge.js
// =============================================================================
// Tests for js/atlas_dosage_bridge.js — the atlas-side installer that
// synthesises state.data.dosage_chunks from a popstats /api/dosage/manifest
// payload. Covers:
//
//   - buildSyntheticIndex(manifest, opts) shape, URL template, cap defaults
//   - install({ manifest, stateRef }) writes state.data.dosage_chunks
//   - install respects pre-existing dosage_chunks (uninstall restores it)
//   - install handles missing manifest gracefully
//   - install handles empty manifest gracefully
//   - install rejects schema_version mismatches
//   - setActiveChrom updates the chrom field without re-probing
//   - setActiveChrom is a no-op for static (file-based) indices
//   - uninstall restores the previous dosage_chunks
//   - getStatus reports current state
//   - URL template carries placeholders __CHROM__/__START__/__END__/__CAP__
//   - Integration with the renderer's _resolveChunkUrl substitution
//     (mock the substitution; verify the URL we build is acceptable)
//
// Run:
//   cd Atlas && node tests/dosage_bridge/test_atlas_dosage_bridge.js
// =============================================================================

'use strict';

const path = require('path');

const MODULE = path.resolve(__dirname, '../../js/atlas_dosage_bridge.js');
const Bridge = require(MODULE);

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  \u2713 ' + msg); }
  else      { fail++; console.error('  \u2717 ' + msg); }
}
function section(name) { console.log('\n\u2014 ' + name + ' \u2014'); }

// ---------------------------------------------------------------------------
section('module shape');

assert(typeof Bridge.install === 'function', 'install is a function');
assert(typeof Bridge.uninstall === 'function', 'uninstall is a function');
assert(typeof Bridge.setActiveChrom === 'function', 'setActiveChrom is a function');
assert(typeof Bridge.getStatus === 'function', 'getStatus is a function');
assert(typeof Bridge.buildSyntheticIndex === 'function', 'buildSyntheticIndex is a function');
assert(Bridge._const.SCHEMA_VERSION === 'dosage_manifest_v1',
       'SCHEMA_VERSION constant exposed');
assert(Bridge._const.MANIFEST_PATH === '/api/dosage/manifest',
       'MANIFEST_PATH = /api/dosage/manifest');
assert(Bridge._const.CHUNK_PATH === '/api/dosage/chunk',
       'CHUNK_PATH = /api/dosage/chunk');

// ---------------------------------------------------------------------------
section('buildSyntheticIndex — basic');

const manifestOk = {
  ok: true,
  schema_version: 'dosage_manifest_v1',
  n_samples: 226,
  chroms: [
    { name: 'C_gar_LG01', length_bp: 41203400 },
    { name: 'C_gar_LG12', length_bp: 25700000 },
    { name: 'C_gar_LG28', length_bp: 28100000 },
  ],
  limits: { default_max_sites: 1000, hard_cap: 20000, max_region_bp: 50000000 },
};

const idx = Bridge.buildSyntheticIndex(manifestOk, {
  baseUrl: 'http://localhost:8765',
  cap: 1000,
});
assert(idx != null, 'index built');
assert(idx.source === 'atlas_server', 'index.source = atlas_server');
assert(idx.schema_version === 'dosage_manifest_v1', 'schema_version preserved');
assert(idx.cap_default === 1000, 'cap_default = 1000');
assert(Array.isArray(idx.chunks) && idx.chunks.length === 3,
       '3 synthetic chunks (one per chrom)');
assert(idx.manifest_n_samples === 226, 'n_samples passes through');
assert(idx.chrom === 'C_gar_LG01',
       'default activeChrom = first chrom in manifest');

// One chunk per chrom: bounds [1, length_bp], URL contains placeholders
const c0 = idx.chunks[0];
assert(c0.chrom === 'C_gar_LG01', 'chunk 0 chrom');
assert(c0.start_bp === 1, 'chunk start = 1');
assert(c0.end_bp === 41203400, 'chunk end = length_bp');
assert(c0.url.indexOf('__CHROM__') !== -1, 'URL has __CHROM__ placeholder');
assert(c0.url.indexOf('__START__') !== -1, 'URL has __START__ placeholder');
assert(c0.url.indexOf('__END__')   !== -1, 'URL has __END__ placeholder');
assert(c0.url.indexOf('__CAP__')   !== -1, 'URL has __CAP__ placeholder');
assert(c0.url.indexOf('http://localhost:8765') === 0,
       'URL prefixed with the base URL');
assert(c0._synthetic === true, 'chunk marked _synthetic');
assert(c0.cap === 1000, 'chunk-local cap = 1000');

// ---------------------------------------------------------------------------
section('buildSyntheticIndex — edge cases');

assert(Bridge.buildSyntheticIndex(null) === null,
       'null manifest returns null');
assert(Bridge.buildSyntheticIndex({ chroms: [] }) === null,
       'empty chroms returns null');
assert(Bridge.buildSyntheticIndex({}) === null,
       'no chroms field returns null');

// length_bp == 0 -> defensive upper bound
const idxZero = Bridge.buildSyntheticIndex(
  { chroms: [{ name: 'C_gar_LG99', length_bp: 0 }] },
  { baseUrl: 'http://h:1' }
);
assert(idxZero.chunks[0].end_bp >= 1e9,
       'zero length_bp -> defensive upper bound (>= 1e9)');

// activeChrom override
const idxAC = Bridge.buildSyntheticIndex(manifestOk, {
  baseUrl: 'http://localhost:8765',
  activeChrom: 'C_gar_LG28',
});
assert(idxAC.chrom === 'C_gar_LG28', 'activeChrom override applied');

// Custom cap
const idxCap = Bridge.buildSyntheticIndex(manifestOk, {
  baseUrl: 'http://localhost:8765',
  cap: 5000,
});
assert(idxCap.cap_default === 5000, 'custom cap_default applied');
assert(idxCap.chunks[0].cap === 5000, 'chunk-local cap matches');

// Trailing-slash baseUrl is normalised
const idxSlash = Bridge.buildSyntheticIndex(manifestOk, {
  baseUrl: 'http://localhost:8765/',
});
assert(idxSlash.chunks[0].url.indexOf('http://localhost:8765/api/dosage/chunk') === 0,
       'trailing slash on baseUrl normalised before /api/dosage/chunk');

// ---------------------------------------------------------------------------
section('install — manifest provided directly');

(async function main() {

// Fresh state with an unrelated dosage_chunks already present
let state = { data: { dosage_chunks: { source: 'static', chunks: [] } } };
let res = await Bridge.install({ manifest: manifestOk, stateRef: state });
assert(res.installed === true, 'install returns installed:true');
assert(res.chroms === 3, 'install reports 3 chroms');
assert(state.data.dosage_chunks.source === 'atlas_server',
       'state.data.dosage_chunks.source = atlas_server');
assert(state.data.dosage_chunks.chunks.length === 3,
       '3 chunks in state.data.dosage_chunks');

const status1 = Bridge.getStatus();
assert(status1.installed === true, 'getStatus reports installed');
assert(status1.chroms === 3, 'getStatus reports 3 chroms');

// uninstall restores the previous index
const restoredRes = Bridge.uninstall({ stateRef: state });
assert(restoredRes.installed === false, 'uninstall returns installed:false');
assert(state.data.dosage_chunks &&
       state.data.dosage_chunks.source === 'static',
       'uninstall restored the previous static index');

// ---------------------------------------------------------------------------
section('install — no prior dosage_chunks');

state = { data: {} };
res = await Bridge.install({ manifest: manifestOk, stateRef: state });
assert(res.installed === true, 'install on empty state succeeds');
assert(state.data.dosage_chunks.source === 'atlas_server',
       'index installed');
Bridge.uninstall({ stateRef: state });
assert(state.data.dosage_chunks == null,
       'uninstall on empty-prior state restores null');

// ---------------------------------------------------------------------------
section('install — refuses bad manifest');

state = { data: {} };
res = await Bridge.install({ manifest: null, stateRef: state });
assert(res.installed === false && res.reason === 'manifest_unavailable',
       'null manifest -> manifest_unavailable');

res = await Bridge.install({ manifest: { chroms: [] }, stateRef: state });
assert(res.installed === false && res.reason === 'manifest_empty',
       'empty manifest -> manifest_empty');

res = await Bridge.install({
  manifest: { schema_version: 'wrong_v9', chroms: [{ name: 'X' }] },
  stateRef: state,
});
assert(res.installed === false && res.reason === 'schema_mismatch',
       'wrong schema_version -> schema_mismatch');
assert(state.data.dosage_chunks == null,
       'state untouched on rejection');

// ---------------------------------------------------------------------------
section('setActiveChrom — switches active chrom on synthetic index');

state = { data: {} };
await Bridge.install({ manifest: manifestOk, stateRef: state });
assert(state.data.dosage_chunks.chrom === 'C_gar_LG01',
       'initial active chrom = C_gar_LG01');
const switched = Bridge.setActiveChrom('C_gar_LG28', { stateRef: state });
assert(switched === true, 'setActiveChrom returns true on success');
assert(state.data.dosage_chunks.chrom === 'C_gar_LG28',
       'active chrom switched to LG28');
assert(Bridge.getStatus().activeChrom === 'C_gar_LG28',
       'getStatus reflects new active chrom');

// ---------------------------------------------------------------------------
section('setActiveChrom — no-op on static index');

state = { data: { dosage_chunks: { source: 'static', chunks: [], chrom: 'X' } } };
const noop = Bridge.setActiveChrom('Y', { stateRef: state });
assert(noop === false, 'setActiveChrom returns false for static index');
assert(state.data.dosage_chunks.chrom === 'X',
       'static index chrom NOT overwritten');

// ---------------------------------------------------------------------------
section('URL template — substitution acceptance');

// Simulate the renderer's _resolveChunkUrl behavior to confirm our URL
// pattern is sane. (Real substitution lives in Inversion_atlas.html.)
function _resolveChunkUrl(chunkRef, region) {
  const u = chunkRef && chunkRef.url;
  if (!u || u.indexOf('__') < 0) return u;
  const r = region || {};
  return u
    .replace('__CHROM__', encodeURIComponent(chunkRef.chrom || ''))
    .replace('__START__', String(r.start_bp != null ? r.start_bp : chunkRef.start_bp))
    .replace('__END__',   String(r.end_bp   != null ? r.end_bp   : chunkRef.end_bp))
    .replace('__CAP__',   String(chunkRef.cap || 1000));
}

state = { data: {} };
await Bridge.install({ manifest: manifestOk, stateRef: state, baseUrl: 'http://h:9' });
const lg28 = state.data.dosage_chunks.chunks.find(c => c.chrom === 'C_gar_LG28');
const resolved = _resolveChunkUrl(lg28, { start_bp: 17000000, end_bp: 17500000 });
assert(resolved.indexOf('chrom=C_gar_LG28') !== -1, 'resolved URL has chrom');
assert(resolved.indexOf('start=17000000') !== -1, 'resolved URL has start');
assert(resolved.indexOf('end=17500000') !== -1, 'resolved URL has end');
assert(resolved.indexOf('cap=1000') !== -1, 'resolved URL has cap');
assert(resolved.indexOf('__') === -1,
       'no unresolved placeholders in final URL');

// ---------------------------------------------------------------------------
section('install — re-install updates without losing the original prior');

state = { data: { dosage_chunks: { source: 'static', tag: 'original' } } };
await Bridge.install({ manifest: manifestOk, stateRef: state });
assert(state.data.dosage_chunks.source === 'atlas_server', 'first install OK');
// Re-install (e.g. after a chrom switch); previous-snapshot should NOT be
// overwritten with the synthetic one.
await Bridge.install({
  manifest: manifestOk,
  stateRef: state,
  activeChrom: 'C_gar_LG12',
});
assert(state.data.dosage_chunks.chrom === 'C_gar_LG12',
       'second install switches active chrom');
Bridge.uninstall({ stateRef: state });
assert(state.data.dosage_chunks &&
       state.data.dosage_chunks.tag === 'original',
       'uninstall restores the ORIGINAL static prior, not the intermediate synthetic one');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);

})().catch(err => {
  console.error('uncaught:', err);
  process.exit(2);
});
