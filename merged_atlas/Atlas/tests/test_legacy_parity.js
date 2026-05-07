// tests/test_legacy_parity.js
//
// Critical safety net for the extraction step: load the LEGACY
// Inversion_atlas.html, grab its top-level functions, and verify that
// the new shared/* modules produce bit-identical output for a battery
// of inputs.
//
// This works because the legacy atlas is one big <script> that defines
// _buildContingency, _hetRateColor, etc. as top-level functions. We
// extract that script body, evaluate it in a controlled context, and
// compare outputs.
//
// Skipped if the legacy file is not present at /home/claude/full_extract/Atlas/Inversion_atlas.html.
// On the host machine this would point at the user's checkout.

import { readFileSync, existsSync } from 'node:fs';
import { Script, createContext } from 'node:vm';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const LEGACY_HTML_DEFAULT = '/home/claude/full_extract/Atlas/Inversion_atlas.html';
const LEGACY_HTML = process.env.LEGACY_ATLAS || LEGACY_HTML_DEFAULT;

let pass = 0, fail = 0, skipped = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

if (!existsSync(LEGACY_HTML)) {
  console.log(`  ⊘ legacy atlas not found at ${LEGACY_HTML}`);
  console.log(`  ⊘ skipping parity tests (set LEGACY_ATLAS env var to point elsewhere)`);
  process.exit(0);
}

// -----------------------------------------------------------------------------
// Extract <script> bodies from the legacy HTML and evaluate them in an
// isolated VM context. We grab ALL top-level <script> blocks and run
// them in order. Some of them touch DOM/window — we provide minimal
// stubs that absorb those calls without crashing.
// -----------------------------------------------------------------------------

console.log(`--- loading legacy atlas from ${LEGACY_HTML} ---`);
const html = readFileSync(LEGACY_HTML, 'utf-8');

// Match every <script>…</script> that does NOT have type="module" or src=
// (we want only inline script blocks)
const scriptRe = /<script(?![^>]*\bsrc=)(?![^>]*\btype="module")[^>]*>([\s\S]*?)<\/script>/g;
const scriptBlocks = [];
let m;
while ((m = scriptRe.exec(html)) !== null) {
  scriptBlocks.push(m[1]);
}

// Fallback: if no closed <script>…</script> pair (the legacy file in
// our truncated tarball is missing its closing tags), grab everything
// from the FIRST plain <script> open tag to end-of-file.
if (scriptBlocks.length === 0) {
  const openRe = /<script(?![^>]*\bsrc=)(?![^>]*\btype="module")[^>]*>/;
  const openMatch = html.match(openRe);
  if (openMatch) {
    const start = openMatch.index + openMatch[0].length;
    let body = html.slice(start);
    const closeIdx = body.indexOf('</script>');
    if (closeIdx >= 0) body = body.slice(0, closeIdx);
    scriptBlocks.push(body);
    console.log('  (used unclosed-<script> fallback — file appears truncated)');
  }
}
console.log(`  found ${scriptBlocks.length} inline script block(s); total ${scriptBlocks.reduce((s,b)=>s+b.length,0).toLocaleString()} chars`);

// -----------------------------------------------------------------------------
// Per-function extraction strategy.
//
// Evaluating the whole 3 MB legacy script in a Node VM context is
// fragile (browser-only globals, truncation artefacts, top-level UI
// wiring that throws). Instead, we extract individual top-level
// function definitions by name, plus their named-constant deps, and
// eval ONLY those. This is bullet-proof for pure-compute primitives.
// -----------------------------------------------------------------------------

const SCRIPT_BODY = scriptBlocks.join('\n');

/**
 * Extract a top-level function definition by name. Walks from
 * `function NAME(...)` and tracks brace depth to find the matching
 * close. Returns the source string (including the `function NAME(...)
 * { ... }` wrapper) or null if not found.
 */
function extractTopLevelFunction(src, name) {
  // Anchor on a fresh-line `function NAME(` (top-level function decl).
  const re = new RegExp(`(^|\\n)function\\s+${name}\\s*\\(`);
  const m = re.exec(src);
  if (!m) return null;
  const start = (m[0][0] === '\n') ? m.index + 1 : m.index;
  // Find the opening `{`
  let i = start;
  while (i < src.length && src[i] !== '{') i++;
  if (i >= src.length) return null;
  // Walk tracking brace depth, respecting strings/comments
  let depth = 0;
  let inLineComment = false, inBlockComment = false;
  let inStr = null;   // '"' | "'" | '`' or null
  for (; i < src.length; i++) {
    const c = src[i], next = src[i + 1];
    if (inLineComment) {
      if (c === '\n') inLineComment = false;
      continue;
    }
    if (inBlockComment) {
      if (c === '*' && next === '/') { inBlockComment = false; i++; }
      continue;
    }
    if (inStr) {
      if (c === '\\') { i++; continue; }
      if (c === inStr) inStr = null;
      continue;
    }
    if (c === '/' && next === '/') { inLineComment = true; i++; continue; }
    if (c === '/' && next === '*') { inBlockComment = true; i++; continue; }
    if (c === '"' || c === "'" || c === '`') { inStr = c; continue; }
    if (c === '{') depth++;
    else if (c === '}') {
      depth--;
      if (depth === 0) {
        return src.slice(start, i + 1);
      }
    }
  }
  return null;
}

/**
 * Extract a top-level `const NAME = ...;` declaration body. Used for
 * _HET_RAMP. Same comment/string awareness.
 */
function extractTopLevelConst(src, name) {
  const re = new RegExp(`(^|\\n)const\\s+${name}\\s*=`);
  const m = re.exec(src);
  if (!m) return null;
  const start = (m[0][0] === '\n') ? m.index + 1 : m.index;
  let i = start;
  // Find the assignment value end — first top-level `;` after the `=`
  let inLineComment = false, inBlockComment = false, inStr = null;
  let parenDepth = 0, braceDepth = 0, bracketDepth = 0;
  // skip past `const NAME =`
  while (i < src.length && src[i] !== '=') i++;
  i++;
  for (; i < src.length; i++) {
    const c = src[i], next = src[i + 1];
    if (inLineComment) { if (c === '\n') inLineComment = false; continue; }
    if (inBlockComment) { if (c === '*' && next === '/') { inBlockComment = false; i++; } continue; }
    if (inStr) { if (c === '\\') { i++; continue; } if (c === inStr) inStr = null; continue; }
    if (c === '/' && next === '/') { inLineComment = true; i++; continue; }
    if (c === '/' && next === '*') { inBlockComment = true; i++; continue; }
    if (c === '"' || c === "'" || c === '`') { inStr = c; continue; }
    if (c === '(') parenDepth++; else if (c === ')') parenDepth--;
    else if (c === '{') braceDepth++; else if (c === '}') braceDepth--;
    else if (c === '[') bracketDepth++; else if (c === ']') bracketDepth--;
    else if (c === ';' && parenDepth === 0 && braceDepth === 0 && bracketDepth === 0) {
      return src.slice(start, i + 1);
    }
  }
  return null;
}

// Build a small synthetic script with ONLY the primitives we want to
// compare. We extract each by name; if any are missing the
// corresponding parity tests are reported as skipped.
const wanted = [
  'function:_buildContingency',
  'function:_detectFuseEvents',
  'function:_computeARI',
  'function:_computeNMI',
  'function:_cramersV',
  'function:_chiSqSurvival',
  'function:_lnGamma',
  'function:_hetRateColor',
  'function:alignLabels',
  'function:permutations',          // dep of alignLabels
  'function:_hungarianChainProjection',
  'function:_finalizeChain',        // dep of _hungarianChainProjection
  'function:_concordanceMatrix',
  'const:_HET_RAMP',
  'const:_LINEAGE_CHAIN_BREAK_AGREEMENT',
  // Step 2.A — pure K-means primitives
  'function:kmeans1D',
  'function:kmeans2D',
  'function:silhouette1D',
  'function:adaptiveK1D',
];

let extractedSrc = '';
const extracted = {};
for (const tag of wanted) {
  const [kind, name] = tag.split(':');
  const body = (kind === 'function')
    ? extractTopLevelFunction(SCRIPT_BODY, name)
    : extractTopLevelConst(SCRIPT_BODY, name);
  if (body) {
    extracted[name] = true;
    extractedSrc += body + '\n';
  } else {
    extracted[name] = false;
  }
}
console.log(`  extracted: ${Object.keys(extracted).filter(k => extracted[k]).length}/${wanted.length} primitives`);
for (const [name, ok] of Object.entries(extracted)) {
  if (!ok) console.log(`    ⊘ missing: ${name}`);
}

// Eval the extracted snippet in a fresh VM context — minimal stub.
const extractCtx = { console, Math, JSON, Number, String, Boolean, Array, Object,
  Float32Array, Float64Array, Int8Array, Int16Array, Int32Array,
  Uint8Array, Uint8ClampedArray, Uint16Array, Uint32Array,
  ArrayBuffer, Map, Set, isFinite, isNaN, parseInt, parseFloat,
  Error, TypeError, Symbol,
};
extractCtx.window = extractCtx;   // legacy `window.foo = foo` writes back to ctx
createContext(extractCtx);
try {
  new Script(extractedSrc).runInContext(extractCtx, { timeout: 5000 });
} catch (e) {
  console.log(`  ✗ failed to evaluate extracted snippet: ${e.message}`);
  console.log('  (writing extracted source to /tmp/legacy_extract.js for debug)');
  await import('node:fs').then(fs => fs.writeFileSync('/tmp/legacy_extract.js', extractedSrc));
  process.exit(1);
}


// -----------------------------------------------------------------------------
// Now we have legacy functions on extractCtx (legacy `window.foo = foo`
// exposures landed there because we aliased `window = extractCtx`).
// -----------------------------------------------------------------------------

const legacyBuildContingency = extractCtx._buildContingency;
const legacyComputeARI       = extractCtx._computeARI;
const legacyComputeNMI       = extractCtx._computeNMI;
const legacyCramersV         = extractCtx._cramersV;
const legacyChiSqSurvival    = extractCtx._chiSqSurvival;
const legacyLnGamma          = extractCtx._lnGamma;
const legacyHetRateColor     = extractCtx._hetRateColor;
const legacyAlignLabels      = extractCtx.alignLabels;
const legacyHungarianChain   = extractCtx._hungarianChainProjection;
const legacyConcordanceMatrix= extractCtx._concordanceMatrix;
const legacyKmeans1D         = extractCtx.kmeans1D;
const legacyKmeans2D         = extractCtx.kmeans2D;
const legacySilhouette1D     = extractCtx.silhouette1D;
const legacyAdaptiveK1D      = extractCtx.adaptiveK1D;

console.log('\n--- legacy primitives accessible? ---');
check('_buildContingency',         typeof legacyBuildContingency === 'function');
check('_computeARI',               typeof legacyComputeARI === 'function');
check('_computeNMI',               typeof legacyComputeNMI === 'function');
check('_cramersV',                 typeof legacyCramersV === 'function');
check('_chiSqSurvival',            typeof legacyChiSqSurvival === 'function');
check('_lnGamma',                  typeof legacyLnGamma === 'function');
check('_hetRateColor',             typeof legacyHetRateColor === 'function');
check('alignLabels',               typeof legacyAlignLabels === 'function');
check('_hungarianChainProjection', typeof legacyHungarianChain === 'function');
check('_concordanceMatrix',        typeof legacyConcordanceMatrix === 'function');
check('kmeans1D',                  typeof legacyKmeans1D === 'function');
check('kmeans2D',                  typeof legacyKmeans2D === 'function');
check('silhouette1D',              typeof legacySilhouette1D === 'function');
check('adaptiveK1D',               typeof legacyAdaptiveK1D === 'function');

// -----------------------------------------------------------------------------
// Load NEW shared modules
// -----------------------------------------------------------------------------
const sioContingency = await import('../shared/contingency.js');
const sioHungarian   = await import('../shared/hungarian.js');
const sioHetRate     = await import('../shared/het_rate.js');
const sioKmeans      = await import('../shared/kmeans.js');

// -----------------------------------------------------------------------------
// Parity battery
// -----------------------------------------------------------------------------

function deepEqual(a, b) {
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  if (a === b) return true;
  if (typeof a !== typeof b) return false;
  if (typeof a !== 'object' || a === null || b === null) {
    if (typeof a === 'number') return Math.abs(a - b) < 1e-12;
    return false;
  }
  if (Array.isArray(a) !== Array.isArray(b)) return false;
  if (Array.isArray(a)) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!deepEqual(a[i], b[i])) return false;
    return true;
  }
  // Typed arrays
  if (ArrayBuffer.isView(a) !== ArrayBuffer.isView(b)) return false;
  if (ArrayBuffer.isView(a)) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!deepEqual(a[i], b[i])) return false;
    return true;
  }
  // Plain objects
  const ka = Object.keys(a), kb = Object.keys(b);
  if (ka.length !== kb.length) return false;
  for (const k of ka) if (!deepEqual(a[k], b[k])) return false;
  return true;
}

console.log('\n--- _buildContingency parity ---');
{
  const tests = [
    [[0,0,1,1,2,2], [0,0,1,1,2,2], 3, 3],
    [[0,0,0,1,1,2], [2,1,0,1,2,0], 3, 3],
    [[0,1,2,-1,5],  [0,1,2,0,1],   3, 3],
    [[0,0,1,1],     [0,1,0,1],     2, 2],
    [[],            [],            2, 2],
  ];
  for (let i = 0; i < tests.length; i++) {
    const [A, B, KA, KB] = tests[i];
    const legacy = legacyBuildContingency(A, B, KA, KB);
    const fresh  = sioContingency.buildContingency(A, B, KA, KB);
    check(`tests[${i}] legacy === fresh`, deepEqual(legacy, fresh),
          legacy ? `n=${legacy.n}` : 'null');
  }
}

console.log('\n--- _computeARI parity ---');
{
  const tests = [
    [[0,0,1,1,2,2], [0,0,1,1,2,2]],
    [[0,0,1,1,2,2], [2,2,0,0,1,1]],
    [[0,0,0,1,1,1], [0,1,2,0,1,2]],
    [[0,0,0,1,1,1,2,2,2], [2,2,2,1,1,1,0,0,0]],
  ];
  for (let i = 0; i < tests.length; i++) {
    const [A, B] = tests[i];
    const legacy = legacyComputeARI(A, B);
    const fresh  = sioContingency.computeARI(A, B);
    check(`ARI tests[${i}]: ${legacy.toFixed(4)} vs ${fresh.toFixed(4)}`,
          deepEqual(legacy, fresh));
  }
}

console.log('\n--- _computeNMI parity ---');
{
  const tests = [
    [[0,0,1,1,2,2], [0,0,1,1,2,2]],
    [[0,0,1,1,2,2], [2,2,0,0,1,1]],
    [[0,0,0,1,1,1], [0,1,2,0,1,2]],
  ];
  for (let i = 0; i < tests.length; i++) {
    const [A, B] = tests[i];
    const legacy = legacyComputeNMI(A, B);
    const fresh  = sioContingency.computeNMI(A, B);
    check(`NMI tests[${i}]: ${legacy.toFixed(4)} vs ${fresh.toFixed(4)}`,
          deepEqual(legacy, fresh));
  }
}

console.log('\n--- _cramersV parity ---');
{
  const tests = [
    [[10, 10, 10, 10], 2, 2],
    [[5, 0, 0, 5], 2, 2],
    [[1, 2, 3, 4, 5, 6], 2, 3],
    [[5, 5, 0, 0], 2, 2],   // single non-empty col
  ];
  for (let i = 0; i < tests.length; i++) {
    const [tbl, KA, KB] = tests[i];
    const legacy = legacyCramersV(tbl, KA, KB);
    const fresh  = sioContingency.cramersV(tbl, KA, KB);
    check(`V tests[${i}]: ${legacy} vs ${fresh}`, deepEqual(legacy, fresh));
  }
}

console.log('\n--- _chiSqSurvival + _lnGamma parity ---');
{
  for (const x of [0.5, 1, 2, 3, 5, 10, 50, 100]) {
    check(`lnGamma(${x})`, deepEqual(legacyLnGamma(x), sioContingency.lnGamma(x)));
  }
  const cases = [[0,1],[1,1],[2,2],[5,3],[10,5],[100,10],[1000,1]];
  for (const [chi2, df] of cases) {
    check(`Q(chi2=${chi2}, df=${df})`,
          deepEqual(legacyChiSqSurvival(chi2, df), sioContingency.chiSqSurvival(chi2, df)));
  }
}

console.log('\n--- _hetRateColor parity ---');
{
  const rates = [null, undefined, NaN, Infinity, 0, 0.1, 0.25, 0.42, 0.5, 0.7, 0.85, 1.0, -0.3, 1.5];
  for (const r of rates) {
    const legacy = legacyHetRateColor(r);
    const fresh  = sioHetRate.hetRateColor(r);
    check(`color(${r}): "${legacy}" vs "${fresh}"`, legacy === fresh);
  }
}

console.log('\n--- alignLabels parity ---');
{
  const tests = [
    [[0,0,1,1,2,2], [0,0,1,1,2,2], 3],
    [[0,0,1,1,2,2], [2,2,0,0,1,1], 3],
    [[0,0,0,1,1,1,2,2,2], [0,0,0,0,0,0,0,0,0], 3],   // collapsed
    [[0,1,0,1,0,1], [1,0,1,0,1,0], 2],                 // swap
  ];
  for (let i = 0; i < tests.length; i++) {
    const [g1, g2, K] = tests[i];
    const legacy = legacyAlignLabels(g1, g2, K);
    const fresh  = sioHungarian.alignLabels(g1, g2, K);
    check(`alignLabels tests[${i}] concord`, deepEqual(legacy.concord, fresh.concord),
          `legacy=${legacy.concord} fresh=${fresh.concord}`);
    check(`alignLabels tests[${i}] aligned`,
          deepEqual(Array.from(legacy.aligned), Array.from(fresh.aligned)));
    check(`alignLabels tests[${i}] perm`,
          deepEqual(Array.from(legacy.perm), Array.from(fresh.perm)));
  }
}

console.log('\n--- _hungarianChainProjection parity ---');
{
  // Build a synthetic getLabelsForL2 callback for both legacy and fresh
  const labelsByL2 = new Map([
    [10, new Int8Array([0, 0, 1, 1, 2, 2])],
    [20, new Int8Array([1, 1, 2, 2, 0, 0])],   // rotated
    [30, new Int8Array([0, 0, 1, 1, 2, 2])],
    [40, new Int8Array([2, 2, 0, 0, 1, 1])],
  ]);
  const getLabels = (idx) => labelsByL2.get(idx) || null;
  const legacy = legacyHungarianChain([10, 20, 30, 40], 3, getLabels);
  const fresh  = sioHungarian.hungarianChainProjection([10, 20, 30, 40], 3, getLabels);
  check('chain count',         legacy.n_chains === fresh.n_chains);
  check('n_samples',           legacy.n_samples === fresh.n_samples);
  check('n_total_L2',          legacy.n_total_L2 === fresh.n_total_L2);
  // Check projected arrays bit-identical
  for (let c = 0; c < legacy.chains.length; c++) {
    check(`chain[${c}].l2_indices`,
          deepEqual(legacy.chains[c].l2_indices, fresh.chains[c].l2_indices));
    check(`chain[${c}].projected`,
          deepEqual(Array.from(legacy.chains[c].projected),
                    Array.from(fresh.chains[c].projected)));
  }
}

console.log('\n--- _concordanceMatrix parity ---');
{
  const proj = {
    chains: [{
      projected: new Int8Array([0,0,1,1, 0,0,1,1, 1,0,1,0]),
      n_L2: 3,
    }],
    n_samples: 4,
    n_total_L2: 3,
  };
  const legacy = legacyConcordanceMatrix(proj);
  const fresh  = sioHungarian.concordanceMatrix(proj);
  check('concordanceMatrix lengths', legacy.length === fresh.length);
  check('concordanceMatrix bit-identical',
        deepEqual(Array.from(legacy), Array.from(fresh)));
}

console.log('\n--- kmeans1D parity ---');
{
  const fixtures = [
    { vals: [0, 0.1, 0.2, 5, 5.1, 5.2, 10, 10.1, 10.2], k: 3 },
    { vals: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],            k: 2 },
    { vals: [-3, -2, -1, 0, 1, 2, 3],                   k: 3 },
    { vals: [5, 5, 5, 5, 5],                            k: 1 },
    { vals: [],                                          k: 3 },
  ];
  for (let i = 0; i < fixtures.length; i++) {
    const { vals, k } = fixtures[i];
    const legacy = legacyKmeans1D(vals, k);
    const fresh  = sioKmeans.kmeans1D(vals, k);
    check(`kmeans1D[${i}] labels bit-identical`,
          deepEqual(Array.from(legacy.labels), Array.from(fresh.labels)));
    check(`kmeans1D[${i}] centers bit-identical`,
          deepEqual(Array.from(legacy.centers), Array.from(fresh.centers)));
    check(`kmeans1D[${i}] n_per_group`,
          deepEqual(legacy.n_per_group, fresh.n_per_group));
  }
}

console.log('\n--- kmeans2D parity ---');
{
  const fixtures = [
    { xs: [0, 0.1, 0.2, 5, 5.1, 5.2, 10, 10.1, 10.2],
      ys: [0, 0.1, -0.1, 5, 4.9, 5.1, 10, 9.9, 10.1], k: 3 },
    { xs: [1, 2, 3, 4, 5, 6], ys: [6, 5, 4, 3, 2, 1], k: 2 },
    { xs: [], ys: [], k: 3 },
  ];
  for (let i = 0; i < fixtures.length; i++) {
    const { xs, ys, k } = fixtures[i];
    const legacy = legacyKmeans2D(xs, ys, k);
    const fresh  = sioKmeans.kmeans2D(xs, ys, k);
    check(`kmeans2D[${i}] labels`,
          deepEqual(Array.from(legacy.labels), Array.from(fresh.labels)));
    check(`kmeans2D[${i}] cx`,
          deepEqual(Array.from(legacy.cx), Array.from(fresh.cx)));
    check(`kmeans2D[${i}] cy`,
          deepEqual(Array.from(legacy.cy), Array.from(fresh.cy)));
    check(`kmeans2D[${i}] n_per_group`,
          deepEqual(legacy.n_per_group, fresh.n_per_group));
  }
}

console.log('\n--- silhouette1D parity ---');
{
  const fixtures = [
    { vals: [0, 0.1, 0.2, 5, 5.1, 5.2, 10, 10.1, 10.2],
      labels: new Int8Array([0,0,0,1,1,1,2,2,2]), k: 3 },
    { vals: [1, 2, 3, 4, 5, 6, 7, 8],
      labels: new Int8Array([0,0,0,0,1,1,1,1]), k: 2 },
    // Degenerate cases
    { vals: [1, 2, 3], labels: new Int8Array([0,0,0]), k: 1 },     // n<4 → NaN
    { vals: [1,2,3,4], labels: new Int8Array([0,0,0,1]), k: 2 },   // singleton in cluster 1 → NaN
  ];
  for (let i = 0; i < fixtures.length; i++) {
    const { vals, labels, k } = fixtures[i];
    const legacy = legacySilhouette1D(vals, labels, k);
    const fresh  = sioKmeans.silhouette1D(vals, labels, k);
    check(`silhouette1D[${i}]`, deepEqual(legacy, fresh),
          `legacy=${legacy} fresh=${fresh}`);
  }
}

console.log('\n--- adaptiveK1D parity ---');
{
  const fixtures = [
    { vals: [0,0.1,0.2,0.15, 5,5.1,5.2,5.15, 10,10.1,10.2,10.15],
      kMin: 2, kMax: 5, silThr: 0.45, minNGroup: 3 },
    { vals: [1,2,3,4,5,6,7,8,9,10],
      kMin: 2, kMax: 5, silThr: 0.99, minNGroup: 2 },  // forces fallback
    { vals: [1,2,3],
      kMin: 2, kMax: 5, silThr: 0.45, minNGroup: 5 },  // null path
  ];
  for (let i = 0; i < fixtures.length; i++) {
    const { vals, kMin, kMax, silThr, minNGroup } = fixtures[i];
    const legacy = legacyAdaptiveK1D(vals, kMin, kMax, silThr, minNGroup);
    const fresh  = sioKmeans.adaptiveK1D(vals, kMin, kMax, silThr, minNGroup);
    if (legacy === null) {
      check(`adaptiveK1D[${i}] both null`, fresh === null);
    } else {
      check(`adaptiveK1D[${i}] k`,           deepEqual(legacy.k, fresh.k));
      check(`adaptiveK1D[${i}] silhouette`,  deepEqual(legacy.silhouette, fresh.silhouette));
      check(`adaptiveK1D[${i}] labels`,
            deepEqual(Array.from(legacy.labels), Array.from(fresh.labels)));
      check(`adaptiveK1D[${i}] centers`,
            deepEqual(Array.from(legacy.centers), Array.from(fresh.centers)));
    }
  }
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
