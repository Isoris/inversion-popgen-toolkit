// inversion_review/main.js
// Smoke-test entry. Real review-page logic ships next turn.

import {
  State, KNOWN_LAYERS, parseTsv, downloadReviewJson,
  assertSharedVersion, SHARED_VERSION,
} from '../shared/state_io.js';

// Workflow packs assert the shared/ version they were built against.
// If shared/ has been replaced with an older pack, this fails loudly.
assertSharedVersion(1, 'inversion_review');

const state = new State({ workflow: 'inversion', baseUrl: './data' });

// One smoke test: parseTsv round-trip
function smokeParseTsv() {
  const sample = '# K = 3\n# JACCARD_MIN = 0.65\nfoo\tbar\tbaz\n1\t2\t3\n4\t5\t6\n';
  const out = parseTsv(sample);
  const ok = out.params.K === '3'
          && out.header.length === 3
          && out.rows.length === 2
          && out.rows[0].foo === '1'
          && out.rows[1].baz === '6';
  return { name: 'parseTsv round-trip', pass: ok, detail: out };
}

// One smoke test: known layer registry has the SPEC BLOCK 1 R-module outputs
function smokeKnownLayers() {
  const stems = KNOWN_LAYERS.precompFiles.map(f => f.stem);
  const required = [
    'band_nodes', 'band_edges', 'transition_events', 'band_trajectories',
    'het_band_backbones', 'candidate_track_proposals', 'candidate_tracks',
    'manual_review_queue',
  ];
  const missing = required.filter(r => !stems.includes(r));
  return {
    name: 'KNOWN_LAYERS.precompFiles has all SPEC BLOCK 1 stems',
    pass: missing.length === 0,
    detail: missing.length ? `missing: ${missing.join(', ')}` : 'all 8 present',
  };
}

// One smoke test: downloadReviewJson headless mode
function smokeDownloadHeadless() {
  const obj = { foo: 1, bar: [2, 3] };
  const r = downloadReviewJson('inversion', 'manual_overrides.json', obj, { headless: true });
  const ok = r.filename === 'inversion__manual_overrides.json'
          && JSON.parse(r.text).bar[1] === 3;
  return { name: 'downloadReviewJson headless', pass: ok, detail: r.filename };
}

function runAndDisplay() {
  const root = document.getElementById('root') || document.body;
  const tests = [smokeParseTsv(), smokeKnownLayers(), smokeDownloadHeadless()];
  const lines = tests.map(t => `${t.pass ? '✓' : '✗'} ${t.name} — ${typeof t.detail === 'string' ? t.detail : JSON.stringify(t.detail)}`);
  const allPass = tests.every(t => t.pass);
  root.innerHTML =
    `<h1 style="font:14px ui-monospace,Menlo,monospace">inversion_review.html — modular pattern smoke test</h1>` +
    `<p style="font:12px sans-serif;color:${allPass ? 'green' : 'red'}">${allPass ? 'all pass' : 'FAIL'}</p>` +
    `<pre style="font:12px ui-monospace,Menlo,monospace;background:#fafafa;padding:8px;border:1px solid #ddd">${lines.join('\n')}</pre>` +
    `<p style="font:12px sans-serif;color:#666">workflow: <code>${state.workflow}</code> · baseUrl: <code>${state.baseUrl}</code></p>`;
  // Also surface results to console so headless harness can read them
  if (typeof window !== 'undefined') {
    window.__SMOKE_RESULTS__ = tests;
  }
}

if (typeof document !== 'undefined') {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', runAndDisplay);
  } else {
    runAndDisplay();
  }
}

// Export for headless test runner (Node-side smoke test)
export { smokeParseTsv, smokeKnownLayers, smokeDownloadHeadless };
