#!/usr/bin/env node
// =============================================================================
// cross_check_M06.js
//
// Cross-check that JSONs produced by STEP_M06_emit_boundary_evidence.R
// satisfy the scrubber's actual layer-detection rules and that
// `_buildBoundaryTrackScores` lights up the four boundary_evidence-derived
// tracks (fst_edge, theta_pi_step, discordant_pile, sv_anchor).
//
// Mirrors the v3.95 pattern (cross_check_scrubber_layers.js +
// end_to_end_check.js for STEP_M04/STEP_M05).
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

// =============================================================================
// Find the M06 smoke-test work dir (stable location set by smoke_test_M06.R)
// =============================================================================
const work = '/tmp/m06_smoke_latest';
if (!fs.existsSync(work)) {
  console.error('FATAL: M06 smoke-test work dir not found at:', work);
  console.error('  Run smoke_test_M06.R first.');
  process.exit(2);
}
console.log('using work dir:', work);

const CHROM = 'C_gar_LG28';
const beJson = path.join(work, 'scrubber_data', CHROM, `${CHROM}_boundary_evidence.json`);
if (!fs.existsSync(beJson)) {
  console.error('FATAL: boundary_evidence JSON missing:', beJson);
  process.exit(2);
}
const be = JSON.parse(fs.readFileSync(beJson, 'utf8'));
console.log('layer keys:', Object.keys(be).filter(k => !k.startsWith('_')).join(', '));
console.log('schema_version:', be.schema_version);
console.log('candidates in layer:', be.boundary_evidence.length);

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { console.log('ok:   ' + msg); pass++; }
  else      { console.log('FAIL: ' + msg); fail++; }
}

// =============================================================================
// Part 1: layer detection — scrubber's rule for boundary_evidence
//   Array.isArray(data.boundary_evidence) AND first row has candidate_id + tracks
// =============================================================================
console.log('\n--- Part 1: layer detection ---');
ok(Array.isArray(be.boundary_evidence),
   'boundary_evidence is an array');
const r0 = be.boundary_evidence[0];
ok(r0 && Number.isInteger(r0.candidate_id),
   'first row has integer candidate_id');
ok(r0 && r0.tracks && typeof r0.tracks === 'object',
   'first row has tracks{} object');
ok(Number.isInteger(r0.scan_start_bp) && Number.isInteger(r0.scan_end_bp),
   'scan_start_bp + scan_end_bp are integers');
ok(Number.isInteger(r0.scan_window_bp) && r0.scan_window_bp > 0,
   'scan_window_bp is positive integer');

// =============================================================================
// Part 2: track shape contract
// =============================================================================
console.log('\n--- Part 2: track shape contract ---');
const expectedLen = Math.floor((r0.scan_end_bp - r0.scan_start_bp) / r0.scan_window_bp);
const t = r0.tracks;

// Each numerical track is array of length n_windows
for (const trackName of ['fst', 'theta_pi_homo1', 'theta_pi_het', 'theta_pi_homo2',
                          'discordant_pair_pileup']) {
  if (Array.isArray(t[trackName])) {
    ok(t[trackName].length === expectedLen,
       `${trackName} length = ${expectedLen} (got ${t[trackName].length})`);
  } else {
    console.log(`  (${trackName} absent — optional)`);
  }
}

// sv_anchors: array of {kind, pos_bp, [ct], [qual]}
if (Array.isArray(t.sv_anchors)) {
  ok(t.sv_anchors.length > 0, `sv_anchors has ${t.sv_anchors.length} entries`);
  const a0 = t.sv_anchors[0];
  ok(typeof a0.kind === 'string' && a0.kind.length > 0,
     `first sv_anchor.kind is string: "${a0.kind}"`);
  ok(Number.isInteger(a0.pos_bp),
     `first sv_anchor.pos_bp is integer: ${a0.pos_bp}`);
}

// =============================================================================
// Part 3: extract REAL _buildBoundaryTrackScores from scrubber HTML, run it
// =============================================================================
console.log('\n--- Part 3: REAL _buildBoundaryTrackScores ---');

const html = fs.readFileSync('/home/claude/v3/pca_scrubber_v3.html', 'utf8');

// Helper: extract a top-level function by brace-balance
function extractFn(src, name) {
  const re = new RegExp(`function ${name}\\s*\\(`);
  const m = src.match(re);
  if (!m) throw new Error(`function ${name} not found`);
  const i = m.index;
  const open = src.indexOf('{', i);
  let depth = 0, j = open;
  for (; j < src.length; j++) {
    if (src[j] === '{') depth++;
    else if (src[j] === '}') { depth--; if (depth === 0) { j++; break; } }
  }
  return src.slice(i, j);
}

// Build a sandbox with state.data containing the boundary_evidence layer,
// plus a synthetic state.data.windows that covers a scan range matching the
// candidate's scan range. This way _buildBoundaryTrackScores can resample
// from M06's 5-kb grid into the L1 window grid.
const ctx = {};
vm.createContext(ctx);

// Build state.data.windows with proper start_bp/end_bp arrays
const N_WIN = 200;   // 200 L1 windows of 25 kb each = 5 Mb
const WIN_BP = 25000;
const wStartArr = new Array(N_WIN);
const wEndArr   = new Array(N_WIN);
const scanStart = r0.scan_start_bp;
for (let i = 0; i < N_WIN; i++) {
  wStartArr[i] = scanStart + i * WIN_BP;
  wEndArr[i]   = scanStart + (i + 1) * WIN_BP - 1;
}

ctx.state = {
  data: {
    boundary_evidence: be.boundary_evidence,
    windows: { start_bp: wStartArr, end_bp: wEndArr },
  },
  __boundaries: null,
};

// Inject the function (and its dependencies — none, _resample is local)
vm.runInContext(extractFn(html, '_buildBoundaryTrackScores'), ctx);

// Run it on a fake candidate matching M06's candidate_id=1
const cand = { id: 1, chrom: CHROM, start_bp: 12000000, end_bp: 14000000 };
const scanRange = { win_lo: 0, win_hi: N_WIN - 1 };
const result = ctx._buildBoundaryTrackScores(cand, scanRange);

ok(result && result.tracks,
   'result has tracks object');
ok(result.len === N_WIN,
   `result.len = ${N_WIN} (got ${result.len})`);

// =============================================================================
// Part 4: check each boundary_evidence-derived track is populated
// =============================================================================
console.log('\n--- Part 4: boundary_evidence-derived tracks ---');
const out = result.tracks;

// Duck-type check for Float64Array — VM contexts have different constructors,
// so `instanceof Float64Array` fails across VM boundary. The `.BYTES_PER_ELEMENT`
// property is set to 8 on Float64Array typed-array instances.
const isFloat64Array = (a) => a != null
  && a.BYTES_PER_ELEMENT === 8
  && typeof a.length === 'number';

// fst_edge — from beRow.tracks.fst
ok(isFloat64Array(out.fst_edge),
   'fst_edge is Float64Array');
ok(out.fst_edge.length === N_WIN,
   `fst_edge length = ${N_WIN} (got ${out.fst_edge.length})`);
const fstFinite = Array.from(out.fst_edge).filter(Number.isFinite);
ok(fstFinite.length > 0,
   `fst_edge has ${fstFinite.length} finite values`);
const fstMean = fstFinite.reduce((a, b) => a + b, 0) / fstFinite.length;
console.log(`  fst_edge: ${fstFinite.length}/${N_WIN} finite, mean ${fstMean.toFixed(4)}`);
ok(fstMean > 0,
   'fst_edge mean is positive (signal present)');

// theta_pi_step — averaged across the three regime-conditional θπ tracks
ok(isFloat64Array(out.theta_pi_step),
   'theta_pi_step is Float64Array');
const tpFinite = Array.from(out.theta_pi_step).filter(Number.isFinite);
ok(tpFinite.length > 0,
   `theta_pi_step has ${tpFinite.length} finite values`);

// discordant_pile — from beRow.tracks.discordant_pair_pileup
ok(isFloat64Array(out.discordant_pile),
   'discordant_pile is Float64Array');
const discFinite = Array.from(out.discordant_pile).filter(Number.isFinite);
ok(discFinite.length > 0,
   `discordant_pile has ${discFinite.length} finite values`);
const discMax = Math.max(...discFinite);
console.log(`  discordant_pile: max value ${discMax}`);
ok(discMax > 0,
   'discordant_pile has positive values (pileup spikes propagated)');

// sv_anchor — 1 per window if any anchor falls in window range, 0 otherwise
ok(isFloat64Array(out.sv_anchor),
   'sv_anchor is Float64Array');
const svPositive = Array.from(out.sv_anchor).filter(v => v > 0);
console.log(`  sv_anchor: ${svPositive.length} windows with anchors (out of ${N_WIN})`);
ok(svPositive.length > 0,
   'sv_anchor has at least one window with anchor present');
ok(svPositive.length <= 5,
   'sv_anchor windows count is sensible (<=5; we put 3 in)');

// =============================================================================
// Part 5: signal landing in the right place
// =============================================================================
console.log('\n--- Part 5: signal landing in right windows ---');

// FST signal should be elevated near 11.95 Mb (left edge) and 14.05 Mb (right edge)
// — find the windows covering those positions
function windowAt(bp) {
  for (let i = 0; i < N_WIN; i++) {
    if (wStartArr[i] <= bp && bp <= wEndArr[i]) return i;
  }
  return -1;
}
const wLeft  = windowAt(11950000);
const wRight = windowAt(14025000);
console.log(`  left edge window = ${wLeft}, right edge window = ${wRight}`);

if (wLeft >= 0 && wRight >= 0) {
  // FST should peak near these windows
  const fstAtLeft  = out.fst_edge[wLeft];
  const fstAtRight = out.fst_edge[wRight];
  console.log(`  FST at left edge: ${fstAtLeft.toFixed(4)}`);
  console.log(`  FST at right edge: ${fstAtRight.toFixed(4)}`);
  ok(Number.isFinite(fstAtLeft) && fstAtLeft > 0.3,
     `FST elevated at left edge (>0.3; got ${fstAtLeft.toFixed(3)})`);
  ok(Number.isFinite(fstAtRight) && fstAtRight > 0.3,
     `FST elevated at right edge (>0.3)`);

  // sv_anchor should fire at exactly the windows containing our 3 SV calls
  ok(out.sv_anchor[wLeft]  > 0,
     'sv_anchor fires at left-edge window (DELLY INV at 11.95 Mb)');
  ok(out.sv_anchor[wRight] > 0,
     'sv_anchor fires at right-edge window (DELLY INV at 14.025 Mb)');
}

// =============================================================================
// Part 6: empty boundary_evidence — should not crash
// =============================================================================
console.log('\n--- Part 6: empty boundary_evidence handled gracefully ---');
const ctx2 = {};
vm.createContext(ctx2);
ctx2.state = {
  data: {
    boundary_evidence: [],   // empty
    windows: { start_bp: wStartArr, end_bp: wEndArr },
  },
  __boundaries: null,
};
vm.runInContext(extractFn(html, '_buildBoundaryTrackScores'), ctx2);
const empty = ctx2._buildBoundaryTrackScores(cand, scanRange);
ok(empty && empty.tracks,
   'empty boundary_evidence returns valid result');
ok(!empty.tracks.fst_edge,
   'empty boundary_evidence does not produce fst_edge');
ok(!empty.tracks.sv_anchor,
   'empty boundary_evidence does not produce sv_anchor');

// =============================================================================
// Part 7: missing candidate — should not crash, no boundary_evidence tracks
// =============================================================================
console.log('\n--- Part 7: missing candidate_id handled gracefully ---');
const cand999 = { id: 999, chrom: CHROM, start_bp: 12000000, end_bp: 14000000 };
const result999 = ctx._buildBoundaryTrackScores(cand999, scanRange);
ok(result999 && result999.tracks,
   'missing candidate_id returns valid result');
ok(!result999.tracks.fst_edge,
   'missing candidate_id does not produce fst_edge (no row matched)');

// =============================================================================
// Summary
// =============================================================================
console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail === 0 ? 0 : 1);
