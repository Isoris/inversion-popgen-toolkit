// tests/test_modular_smoke.js
// Run with: node --experimental-vm-modules tests/test_modular_smoke.js
//   (or just: node tests/test_modular_smoke.js  on Node 22+)
//
// Verifies:
//   1. shared/state_io.js loads cleanly as an ES module in Node
//   2. parseTsv works
//   3. KNOWN_LAYERS contains all SPEC BLOCK 1 stems
//   4. downloadReviewJson(headless: true) round-trips correctly
//   5. State class instantiates with right defaults
//   6. flatten.py output runs through node --check without parse errors

import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '..');

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) {
    console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`);
    pass++;
  } else {
    console.log(`  ✗ ${name}  ${detail}`);
    fail++;
  }
}

console.log('--- imports ---');
const sio = await import(resolve(ROOT, 'shared/state_io.js'));
check('shared/state_io.js loads', !!sio);
check('exports KNOWN_LAYERS',     typeof sio.KNOWN_LAYERS === 'object');
check('exports parseTsv',         typeof sio.parseTsv === 'function');
check('exports State',            typeof sio.State === 'function');
check('exports downloadReviewJson', typeof sio.downloadReviewJson === 'function');
check('exports pathPrecompMain',  typeof sio.pathPrecompMain === 'function');
check('exports SHARED_VERSION',   typeof sio.SHARED_VERSION === 'number');
check('exports assertSharedVersion', typeof sio.assertSharedVersion === 'function');

console.log('\n--- SHARED_VERSION assertion ---');
let threw1 = false;
try { sio.assertSharedVersion(sio.SHARED_VERSION, 'test'); } catch { threw1 = true; }
check('assert at current version succeeds', !threw1);
let threw2 = false;
try { sio.assertSharedVersion(sio.SHARED_VERSION + 1, 'test'); } catch { threw2 = true; }
check('assert at higher version throws',    threw2);
let threw3 = false;
try { sio.assertSharedVersion(0, 'test'); } catch { threw3 = true; }
check('assert at lower version succeeds',   !threw3);

console.log('\n--- workflow tagging on layers ---');
const inversionPrecomp = sio.KNOWN_LAYERS.precompFiles.filter(f => f.workflow === 'inversion');
check('all 8 SPEC BLOCK 1 stems tagged inversion', inversionPrecomp.length === 8);
const corePrecomp = sio.KNOWN_LAYERS.precompFiles.filter(f => f.workflow === 'core');
check('repeat_density tagged core',
      corePrecomp.some(f => f.stem === 'repeat_density.scrubber_windows'));
check('all cohort entries tagged core',
      sio.KNOWN_LAYERS.cohort.every(f => f.workflow === 'core'));
const candidatesByWorkflow = {};
for (const c of sio.KNOWN_LAYERS.candidate) {
  candidatesByWorkflow[c.workflow] = (candidatesByWorkflow[c.workflow] || 0) + 1;
}
check('candidate entries tagged inversion (currently)',
      candidatesByWorkflow.inversion === sio.KNOWN_LAYERS.candidate.length);

console.log('\n--- parseTsv ---');
const sample = '# K = 3\n# JACCARD_MIN = 0.65\nfoo\tbar\tbaz\n1\t2\t3\n4\t5\t6\n';
const t = sio.parseTsv(sample);
check('parses 2 params',      Object.keys(t.params).length === 2);
check('K param',              t.params.K === '3');
check('header has 3 cols',    t.header.length === 3);
check('rows is length 2',     t.rows.length === 2);
check('row 0 col foo = "1"',  t.rows[0].foo === '1');
check('row 1 col baz = "6"',  t.rows[1].baz === '6');
check('empty input ok',       sio.parseTsv('').rows.length === 0);
check('null input ok',        sio.parseTsv(null).rows.length === 0);

console.log('\n--- KNOWN_LAYERS ---');
const stems = sio.KNOWN_LAYERS.precompFiles.map(f => f.stem);
const required = [
  'band_nodes', 'band_edges', 'transition_events', 'band_trajectories',
  'het_band_backbones', 'candidate_track_proposals', 'candidate_tracks',
  'manual_review_queue',
];
for (const r of required) {
  check(`precompFiles includes ${r}`, stems.includes(r));
}
check('cohort layers include relatedness',
      sio.KNOWN_LAYERS.cohort.some(c => c.name === 'relatedness'));
check('reviewInversion includes manual_overrides.json',
      sio.KNOWN_LAYERS.reviewInversion.includes('manual_overrides.json'));

console.log('\n--- path helpers ---');
check('pathPrecompMain', sio.pathPrecompMain('./data', 'LG28') === './data/precomp/LG28/LG28.json');
check('pathPrecompLayer json',
      sio.pathPrecompLayer('./data', 'LG28', 'band_tracks', 'json')
        === './data/precomp/LG28/LG28.band_tracks.json');
check('pathPrecompLayer tsv',
      sio.pathPrecompLayer('./data', 'LG28', 'band_edges', 'tsv')
        === './data/precomp/LG28/LG28.band_edges.tsv');
check('pathCohort',     sio.pathCohort('./data', 'relatedness.json') === './data/cohort/relatedness.json');
check('pathCandidate',  sio.pathCandidate('./data', 'cand_LG28_15Mb', 'sv_genotype_counts.json')
                          === './data/candidates/cand_LG28_15Mb/sv_genotype_counts.json');
check('pathReview',     sio.pathReview('./data', 'inversion', 'manual_overrides.json')
                          === './data/review/inversion/manual_overrides.json');
check('trailing slash on baseUrl tolerated',
      sio.pathPrecompMain('./data/', 'LG28') === './data/precomp/LG28/LG28.json');

console.log('\n--- downloadReviewJson (headless) ---');
const r = sio.downloadReviewJson('inversion', 'manual_overrides.json',
                                  { foo: 1, bar: [2, 3] }, { headless: true });
check('filename has workflow prefix',  r.filename === 'inversion__manual_overrides.json');
check('text is valid JSON',            JSON.parse(r.text).foo === 1);
check('text round-trips arrays',       JSON.parse(r.text).bar[1] === 3);

console.log('\n--- downloadReviewTsv (headless) ---');
const tsv = sio.downloadReviewTsv('inversion', 'queue.tsv', {
  params: { K: 3, mode: 'het_anchored' },
  header: ['cid', 'class', 'n_carriers'],
  rows: [
    { cid: 'cand_LG28_15Mb', class: 'clean_skeleton', n_carriers: 47 },
    { cid: 'cand_LG28_22Mb', class: 'fragmented',     n_carriers: 12 },
  ],
}, { headless: true });
check('tsv filename', tsv.filename === 'inversion__queue.tsv');
const tsvLines = tsv.text.split('\n');
check('first line is # K param',     tsvLines[0] === '# K = 3');
check('second line is # mode param', tsvLines[1] === '# mode = het_anchored');
check('third line is header',        tsvLines[2] === 'cid\tclass\tn_carriers');
check('fourth line is row 1',        tsvLines[3] === 'cand_LG28_15Mb\tclean_skeleton\t47');

console.log('\n--- State class ---');
const state = new sio.State({ workflow: 'inversion' });
check('default baseUrl is ./data',   state.baseUrl === './data');
check('workflow set',                state.workflow === 'inversion');
check('cohort initially null',       state.cohort === null);
check('precomp initially empty',     Object.keys(state.precomp).length === 0);
check('hasLayer returns false for unknown chrom',
      state.hasLayer('LG99', 'windows') === false);
check('layerNames returns [] for unknown chrom',
      Array.isArray(state.layerNames('LG99')) && state.layerNames('LG99').length === 0);

let threw = false;
try { new sio.State({}); } catch (e) { threw = true; }
check('State without workflow throws', threw);

console.log('\n--- flatten.py output node-checks (optional — run `python3 build/flatten.py inversion_review.html` first) ---');
// Verify the flattened HTML's inlined JS is at minimum syntactically
// valid by extracting the <script> body and running node --check.
const flatPath = resolve(ROOT, 'dist/inversion_review_flat.html');
let flatExists = false;
try { readFileSync(flatPath, 'utf-8'); flatExists = true; } catch {}

if (!flatExists) {
  console.log(`  ⊘ flat HTML not built yet — skipping (run: python3 build/flatten.py inversion_review.html)`);
} else {
  check('flatten.py output exists', flatExists, flatPath);
  const flat = readFileSync(flatPath, 'utf-8');
  const m = flat.match(/<script type="module">([\s\S]*?)<\/script>/);
  check('flat HTML has one inlined <script type="module">', !!m);
  if (m) {
    // Write the script body to a temp file and run `node --check`
    const tmp = '/tmp/flat_module_check.mjs';
    const { writeFileSync } = await import('node:fs');
    writeFileSync(tmp, m[1], 'utf-8');
    const r = spawnSync('node', ['--check', tmp]);
    check('node --check passes on flattened bundle',
          r.status === 0,
          r.status === 0 ? 'syntax ok' : (r.stderr.toString().split('\n')[0] || 'syntax error'));
  }
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
