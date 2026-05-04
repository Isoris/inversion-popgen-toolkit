// =============================================================================
// turn 118 integration test — enriched manuscript bundle export.
//
// Tests the new _bundle*Block helpers and the enriched _candidateTSV under
// several state-loading scenarios:
//   1. No enrichment layers loaded → bundle gracefully degrades to the
//      pre-118 output shape (back-compat)
//   2. Each layer loaded individually → the corresponding block appears
//      and only that block
//   3. All layers loaded → full enrichment renders correctly
//   4. Edge case: candidate not present in any enrichment layer → still
//      renders the base block, no enrichment sub-blocks
// =============================================================================

const fs = require('fs');
const vm = require('vm');

// Extract just the helper functions we want to test from the atlas inline JS
// (the simpler path is to read the html, pull inline scripts, eval in a
// shim context, and call exported functions). But the atlas uses many
// state-dependent helpers (sigmaProfileCandidate, candidateBandComposition,
// _mpEnsureState, etc.) so we'll re-implement just the enrichment helpers
// here using the SAME logic and test them in isolation.
//
// This is a re-implementation test (not a runtime test against the live
// atlas) — the canonical assertion is that the logic in the atlas module
// matches what we've reimplemented here. We verify by reading the source
// substring and confirming the patterns are present.

// First, syntactic / source-level checks — confirm the helpers exist in
// the atlas with the expected names and call patterns
const html = fs.readFileSync('/home/claude/work/build/Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// --- Source-level checks
console.log('\n=== Source-level checks ===');
ok('_bundleBoundaryBlock defined',
   /function _bundleBoundaryBlock\(c\)/.test(html));
ok('_bundleTierBlock defined',
   /function _bundleTierBlock\(c\)/.test(html));
ok('_bundlePopstatsAggregateBlock defined',
   /function _bundlePopstatsAggregateBlock\(c\)/.test(html));
ok('_bundleCsOverlapBlock defined',
   /function _bundleCsOverlapBlock\(c\)/.test(html));
ok('_bundleMarkerReadinessBlock defined',
   /function _bundleMarkerReadinessBlock\(c\)/.test(html));
ok('_bundleFocalVsBgBlock defined',
   /function _bundleFocalVsBgBlock\(c\)/.test(html));
ok('_bundleEnrichmentForCandidate dispatcher defined',
   /function _bundleEnrichmentForCandidate\(c\)/.test(html));
ok('_bundleEnrichmentScalarsForCandidate defined',
   /function _bundleEnrichmentScalarsForCandidate\(c\)/.test(html));
ok('_bundleEnrichmentForCandidate is called from _perCandidateDetailMarkdown',
   /enrichBlocks = _bundleEnrichmentForCandidate\(c\)/.test(html));
ok('_bundleEnrichmentScalarsForCandidate is called from _candidateTSV',
   /_bundleEnrichmentScalarsForCandidate\(c\)/.test(html));

// --- TSV column header checks (the canonical contract for downstream R)
console.log('\n=== TSV column contract ===');
const expectedNewCols = [
  'boundary_left_zone_start_bp', 'boundary_left_zone_end_bp',
  'boundary_right_zone_start_bp', 'boundary_right_zone_end_bp',
  'breakpoint_status', 'confidence_tier',
  'popstats_n_windows_inside', 'popstats_peak_z',
  'popstats_mean_fst_hom1_hom2', 'popstats_mean_hobs_hexp',
  'popstats_mean_theta_pi_log2_ratio',
  'cs_overlap_n', 'cs_overlap_ids',
  'marker_best_tier', 'marker_n_total',
  'marker_n_t1', 'marker_n_t2', 'marker_n_t3', 'marker_n_t4',
  'fvb_p_z', 'fvb_p_fst_hom1_hom2',
  'fvb_p_hobs_hexp_ratio', 'fvb_p_theta_pi_log2_ratio',
];
for (const col of expectedNewCols) {
  ok('TSV column "' + col + '" present',
     html.indexOf("'" + col + "'") >= 0);
}

// --- Enrichment-layers status block in source-data summary
console.log('\n=== Enrichment status block ===');
ok('Source-data summary has the enrichment layers block',
   html.indexOf('Enrichment layers loaded at export time') >= 0);

// --- Caveat additions
console.log('\n=== Caveat additions ===');
ok('Popstats aggregates caveat added',
   /Popstats aggregates.*peak Z, mean F_ST/.test(html));
ok('Tier classification caveat added',
   /atlas does NOT compute tier itself/.test(html));
ok('Focal-vs-bg multiple-testing caveat added',
   /Multiple-testing across many candidates is NOT corrected/.test(html));
ok('Cross-species coincidence framing caveat added',
   /recurrent fragility, not shared inheritance/.test(html));

// =============================================================================
// Behavioural test: build a tiny vm context with the helpers + a fake state,
// and verify each block renders correctly under different layer-loaded
// scenarios.
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

// Pull just the enrichment helpers + dependencies out of the atlas source
// using a regex match. We need: _bundleBoundaryBlock, _bundleTierBlock,
// _bundlePopstatsAggregateBlock, _bundleCsOverlapBlock,
// _bundleMarkerReadinessBlock, _bundleFocalVsBgBlock,
// _bundleEnrichmentForCandidate, _bundleEnrichmentScalarsForCandidate.
function pullFunction(src, fnName) {
  // Match `function fnName(...) { ... }` with brace-balancing
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  // Find the opening brace after the args
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      // skip a string
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) { break; }
        if (quote === '`' && src[i] === '$' && src[i+1] === '{') {
          // skip template literal interpolation by tracking braces
          i += 2; let d = 1;
          while (i < src.length && d > 0) {
            if (src[i] === '{') d++;
            else if (src[i] === '}') d--;
            i++;
          }
          continue;
        }
        i++;
      }
    } else if (ch === '/' && src[i+1] === '/') {
      // line comment
      while (i < src.length && src[i] !== '\n') i++;
    } else if (ch === '/' && src[i+1] === '*') {
      i += 2;
      while (i < src.length - 1 && !(src[i] === '*' && src[i+1] === '/')) i++;
      i++;
    }
    i++;
  }
  return src.substring(start, i);
}

const fnNames = [
  '_bundleBoundaryBlock',
  '_bundleTierBlock',
  '_bundlePopstatsAggregateBlock',
  '_bundleCsOverlapBlock',
  '_bundleMarkerReadinessBlock',
  '_bundleFocalVsBgBlock',
  '_bundleEnrichmentForCandidate',
  '_bundleEnrichmentScalarsForCandidate',
];
const fnSrcs = [];
for (const fn of fnNames) {
  const s = pullFunction(html, fn);
  if (!s) {
    fail++;
    console.log('  FAIL could not extract ' + fn + ' from atlas source');
  } else {
    fnSrcs.push(s);
  }
}

if (fnSrcs.length === fnNames.length) {
  // Build a sandbox with a fake state + window + minimal _mpEnsureState
  function makeSandbox(stateOverlay) {
    const state = Object.assign({
      data: null, popstatsLive: null, crossSpecies: null, candidateList: [],
    }, stateOverlay);
    const sandbox = {
      state,
      console,
      window: {},   // no popgenFocalVsBg by default
      Number, Array, Math, JSON, Object, String,
      _mpEnsureState: () => null,   // default: no marker panel state
    };
    // Allow override of window.popgenFocalVsBg
    if (stateOverlay && stateOverlay.__window_popgenFocalVsBg) {
      sandbox.window.popgenFocalVsBg = stateOverlay.__window_popgenFocalVsBg;
    }
    if (stateOverlay && stateOverlay.__mpPanel != null) {
      sandbox.window.__mpPanel = stateOverlay.__mpPanel;
      sandbox._mpEnsureState = () => ({ panel: stateOverlay.__mpPanel });
    }
    vm.createContext(sandbox);
    for (const s of fnSrcs) {
      vm.runInContext(s, sandbox);
    }
    return sandbox;
  }

  // --- Test 1: empty state, no candidate enrichment → all blocks return null
  console.log('\nTest 1: empty state → no enrichment blocks emitted');
  const s1 = makeSandbox({});
  const cand1 = { id: 'c1', candidate_id: 'c1', chrom: 'LG28',
                  start_bp: 15000000, end_bp: 18000000, K: 3 };
  const blocks1 = vm.runInContext('_bundleEnrichmentForCandidate(state.__cand1)',
    Object.assign(s1, { state: Object.assign(s1.state, { __cand1: cand1 }) }));
  // Re-run properly with cand1 visible inside the helper closure
  s1.cand1 = cand1;
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' +
    JSON.stringify(cand1) + ');', s1);
  ok('empty state → 0 blocks',
     s1.__out && s1.__out.length === 0,
     'got ' + (s1.__out && s1.__out.length) + ' blocks');

  // --- Test 2: only boundary annotation present → exactly 1 block
  console.log('\nTest 2: boundary annotation only → 1 block (boundary)');
  const cand2 = {
    id: 'c2', candidate_id: 'c2', chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3,
    boundary_left:  { zone_start_bp: 14990000, zone_end_bp: 15010000, score: 0.85,
                      source: 'sv_anchor', support_class: 'STRONG',
                      sv_anchors_in_zone: [{type:'DEL'}], notes: '' },
    boundary_right: { zone_start_bp: 17990000, zone_end_bp: 18010000, score: 0.72,
                      source: 'flanking_pc1', support_class: 'MODERATE',
                      sv_anchors_in_zone: [], notes: '' },
    breakpoint_status: 'sharp', boundary_notes: 'cleaner left edge than right',
  };
  const s2 = makeSandbox({});
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand2) + ');', s2);
  ok('boundary only → 1 block', s2.__out.length === 1, 'got ' + s2.__out.length);
  ok('block mentions zone_start_bp', s2.__out[0].indexOf('zone_start_bp') >= 0);
  ok('block mentions sharp status', s2.__out[0].indexOf('sharp') >= 0);

  // --- Test 3: tier classification only
  console.log('\nTest 3: tier only → 1 block (tier)');
  const cand3 = { id: 'c3', candidate_id: 'c3', chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3 };
  const s3 = makeSandbox({
    data: { final_classification: { c3: { confidence_tier: 1, existence_a: 'VALIDATED',
            boundary_quality: 'sharp', mechanism: 'TE-mediated' } } },
  });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand3) + ');', s3);
  ok('tier only → 1 block', s3.__out.length === 1, 'got ' + s3.__out.length);
  ok('block mentions confidence_tier', s3.__out[0].indexOf('confidence_tier') >= 0);
  ok('block mentions VALIDATED', s3.__out[0].indexOf('VALIDATED') >= 0);

  // --- Test 4: popstats aggregates
  console.log('\nTest 4: popstats only → 1 block (popstats aggregate)');
  const cand4 = { id: 'c4', candidate_id: 'c4', chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3 };
  const wins4 = [];
  // Build 50 windows: 25 inside the inversion, 25 outside
  for (let i = 0; i < 50; i++) {
    const mb = i * 0.5 + 5;   // 5..29.5 Mb in 0.5 increments
    wins4.push({
      center_mb: mb, z: i % 7, fst_hom1_hom2: 0.1 + (i % 5) / 10,
      hobs: 0.3, hexp: 0.4, theta_pi_hom_inv: 0.005, theta_pi_hom_ref: 0.003,
    });
  }
  const s4 = makeSandbox({ popstatsLive: { lastResponse: { windows: wins4 } } });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand4) + ');', s4);
  ok('popstats only → 1 block', s4.__out.length === 1, 'got ' + s4.__out.length);
  ok('block mentions peak |Z|', s4.__out[0].indexOf('peak |Z|') >= 0);
  ok('block reports n_windows_inside', s4.__out[0].indexOf('n_windows_inside') >= 0);

  // Verify the scalar export
  vm.runInContext('var __scalars = _bundleEnrichmentScalarsForCandidate(' +
    JSON.stringify(cand4) + ');', s4);
  ok('popstats scalar peak_z is finite',
     s4.__scalars.popstats_peak_z != null && Number.isFinite(s4.__scalars.popstats_peak_z));
  ok('popstats scalar n_windows_inside > 0',
     s4.__scalars.popstats_n_windows_inside > 0,
     'got ' + s4.__scalars.popstats_n_windows_inside);
  ok('popstats theta ratio is finite (constant data)',
     s4.__scalars.popstats_mean_theta_pi_log2_ratio != null &&
     Number.isFinite(s4.__scalars.popstats_mean_theta_pi_log2_ratio));

  // --- Test 5: cs overlap
  console.log('\nTest 5: cs overlap only → 1 block (cs)');
  const cand5 = { id: 'LG28_cand_01', candidate_id: 'LG28_cand_01',
                  chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3 };
  const s5 = makeSandbox({
    crossSpecies: {
      breakpoints: [
        { id: 'cs_bp_42', gar_chr: 'LG28', gar_pos_mb: 16.5, event_type: 'INV',
          candidate_overlap: ['LG28_cand_01'], manuscript_note: 'recurrent fragility' },
        { id: 'cs_bp_43', gar_chr: 'LG28', gar_pos_mb: 22.1, event_type: 'TRA',
          candidate_overlap: ['some_other_cand'] },
      ],
    },
  });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand5) + ');', s5);
  ok('cs only → 1 block', s5.__out.length === 1, 'got ' + s5.__out.length);
  ok('block mentions cs_bp_42', s5.__out[0].indexOf('cs_bp_42') >= 0);
  ok('block does NOT mention cs_bp_43', s5.__out[0].indexOf('cs_bp_43') < 0);
  ok('block has the manuscript framing',
     s5.__out[0].indexOf('recurrent fragility hotspot') >= 0);

  // --- Test 6: marker readiness
  console.log('\nTest 6: marker panel only → 1 block (markers)');
  const cand6 = { id: 'c6', candidate_id: 'c6', chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3 };
  const s6 = makeSandbox({
    __mpPanel: [
      { candidate_id: 'c6', tier: 1, validation_status: 'pilot_tested',
        controls: { positive_controls_INV: ['s1'], negative_controls_STD: ['s2'] } },
      { candidate_id: 'c6', tier: 2, validation_status: 'untested', controls: {} },
      { candidate_id: 'c6', tier: 1, validation_status: 'untested', controls: {} },
      { candidate_id: 'OTHER', tier: 1, validation_status: 'pilot_tested', controls: {} },
    ],
  });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand6) + ');', s6);
  ok('markers only → 1 block', s6.__out.length === 1, 'got ' + s6.__out.length);
  ok('best tier = 1', s6.__out[0].indexOf('Tier 1') >= 0);
  ok('T1=2, T2=1', s6.__out[0].indexOf('T1=2') >= 0 && s6.__out[0].indexOf('T2=1') >= 0);

  // --- Test 7: focal-vs-bg
  console.log('\nTest 7: focal-vs-bg only → 1 block (fvb)');
  const cand7 = { id: 'c7', candidate_id: 'c7', chrom: 'LG28', start_bp: 15000000, end_bp: 18000000, K: 3 };
  // Build a fake popgenFocalVsBg that returns a single metric with clear separation
  const s7 = makeSandbox({
    __window_popgenFocalVsBg: {
      availableMetrics: () => [{
        id: 'z', label: 'robust |Z|',
        getData: () => ({
          mb:     [ 5,  8, 10, 14, 16, 17, 18, 22, 28],   // 5 outside, 4 inside (actually 3 inside: 16, 17, 18 only barely; let me say 14..18 inside, so 4 inside)
          values: [ 1,  1,  1,  1,  4,  4,  4,  1,  1],   // focal mean ~3.25, bg mean ~1
        }),
      }],
      computeFocalBgStats: (opts) => ({
        focal_arr: [4, 4, 4], bg_arr: [1, 1, 1, 1, 1, 1],
        n_focal: 3, n_bg: 6, focal_mean: 4, bg_mean: 1,
        wilcoxon: { p_two_sided: 0.012, z_stat: 2.51, n_focal: 3, n_bg: 6 },
        scatter_xy: [],
      }),
    },
  });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand7) + ');', s7);
  ok('fvb only → 1 block', s7.__out.length === 1, 'got ' + s7.__out.length);
  ok('block mentions Wilcoxon', s7.__out[0].indexOf('Wilcoxon') >= 0);
  ok('block has metric label', s7.__out[0].indexOf('robust |Z|') >= 0);
  ok('block has P value', s7.__out[0].indexOf('0.012') >= 0);

  // --- Test 8: all layers loaded → 6 blocks
  console.log('\nTest 8: all layers loaded → 6 blocks');
  const cand8 = Object.assign({}, cand2, {
    id: 'LG28_cand_01', candidate_id: 'LG28_cand_01',
  });
  const s8 = makeSandbox({
    data: { final_classification: { LG28_cand_01: { confidence_tier: 1, mechanism: 'TE-mediated' } } },
    popstatsLive: { lastResponse: { windows: wins4 } },
    crossSpecies: {
      breakpoints: [{
        id: 'cs_bp_42', gar_chr: 'LG28', gar_pos_mb: 16.5, event_type: 'INV',
        candidate_overlap: ['LG28_cand_01'], manuscript_note: 'recurrent fragility',
      }],
    },
    __mpPanel: [{ candidate_id: 'LG28_cand_01', tier: 1, validation_status: 'pilot_tested', controls: {} }],
    __window_popgenFocalVsBg: {
      availableMetrics: () => [{
        id: 'z', label: 'robust |Z|',
        getData: () => ({ mb: [16, 17], values: [4, 4] }),
      }],
      computeFocalBgStats: () => ({
        focal_arr: [4, 4], bg_arr: [1, 1], n_focal: 2, n_bg: 2,
        focal_mean: 4, bg_mean: 1,
        wilcoxon: { p_two_sided: 0.05, z_stat: 1.96, n_focal: 2, n_bg: 2 },
        scatter_xy: [],
      }),
    },
  });
  vm.runInContext('var __out = _bundleEnrichmentForCandidate(' + JSON.stringify(cand8) + ');', s8);
  ok('all layers → 6 blocks', s8.__out.length === 6, 'got ' + s8.__out.length);
  ok('boundary block is FIRST',         /Boundary refinement/.test(s8.__out[0]));
  ok('tier block is SECOND',            /14-axis tier classification/.test(s8.__out[1]));
  ok('popstats block is THIRD',         /Popstats aggregates/.test(s8.__out[2]));
  ok('cs block is FOURTH',              /Cross-species coincidence/.test(s8.__out[3]));
  ok('marker block is FIFTH',           /Marker readiness/.test(s8.__out[4]));
  ok('fvb block is SIXTH',              /Focal-vs-background Wilcoxon/.test(s8.__out[5]));

  // --- Test 9: scalars round-trip — compare scalar output across scenarios
  console.log('\nTest 9: scalar round-trip across scenarios');
  vm.runInContext('var __scalars = _bundleEnrichmentScalarsForCandidate(' +
    JSON.stringify(cand8) + ');', s8);
  ok('scalar boundary_left_zone_start_bp set',
     s8.__scalars.boundary_left_zone_start_bp === 14990000);
  ok('scalar confidence_tier = 1',
     s8.__scalars.confidence_tier === '1');
  ok('scalar cs_overlap_n = 1', s8.__scalars.cs_overlap_n === 1);
  ok('scalar cs_overlap_ids = "cs_bp_42"',
     s8.__scalars.cs_overlap_ids === 'cs_bp_42');
  ok('scalar marker_best_tier = 1', s8.__scalars.marker_best_tier === 1);
  ok('scalar marker_n_total = 1', s8.__scalars.marker_n_total === 1);
  ok('scalar fvb_p_z = 0.05', Math.abs(s8.__scalars.fvb_p_z - 0.05) < 1e-9);

  // --- Test 10: empty scalars when nothing loaded
  console.log('\nTest 10: scalars all-null with empty state');
  vm.runInContext('var __scalars = _bundleEnrichmentScalarsForCandidate(' +
    JSON.stringify({id:'cE', candidate_id:'cE', start_bp:0, end_bp:1e6, K:3}) + ');', s1);
  ok('boundary_left_zone_start_bp null', s1.__scalars.boundary_left_zone_start_bp == null);
  ok('confidence_tier null', s1.__scalars.confidence_tier == null);
  ok('popstats_peak_z null', s1.__scalars.popstats_peak_z == null);
  ok('cs_overlap_n null', s1.__scalars.cs_overlap_n == null);
  ok('marker_best_tier null', s1.__scalars.marker_best_tier == null);
  ok('fvb_p_z null', s1.__scalars.fvb_p_z == null);
}

// =============================================================================
console.log('\n=========================================');
console.log('passed: ' + pass + ' / failed: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
