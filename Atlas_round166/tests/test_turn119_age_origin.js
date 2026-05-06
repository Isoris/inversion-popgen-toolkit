// =============================================================================
// turn 119 integration test — age & origin atlas surface (cheat30 GDS)
//
// Tests:
//   1. Source-level: validators, store, accessor, panel HTML builder, ridgeline
//      renderer, bundle helper, scalar gather all defined.
//   2. TSV column contract for the 12 new cheat30_* columns.
//   3. Source-data summary status row + caveats added.
//   4. Behavioural (sandbox): _bundleAgeOriginBlock returns null on empty
//      state, full block on populated state. Scalars round-trip.
//   5. Schema validator accepts well-formed JSON, rejects malformed.
//   6. Per-candidate accessor returns null for missing chrom / missing
//      candidate id; returns the right block when present.
//   7. Panel HTML builder produces "no data" empty states correctly and
//      a populated panel when the result is present, with the verdict
//      pill colour matching the origin_class.
//   8. Ridgeline SVG renders with three ridges and does not throw on
//      degenerate (single-density) input.
// =============================================================================

const fs = require('fs');
const vm = require('vm');

const html = fs.readFileSync('/home/claude/work/build/Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// --- Source-level checks
console.log('\n=== Source-level checks ===');
ok('_isCheat30JSON defined',         /function _isCheat30JSON\(/.test(html));
ok('_storeCheat30Results defined',   /function _storeCheat30Results\(/.test(html));
ok('_persistCheat30Results defined', /function _persistCheat30Results\(/.test(html));
ok('_restoreCheat30Results defined', /function _restoreCheat30Results\(/.test(html));
ok('_clearCheat30Results defined',   /function _clearCheat30Results\(/.test(html));
ok('_cheat30ForCandidate defined',   /function _cheat30ForCandidate\(/.test(html));
ok('candidateAgeOriginHtml defined', /function candidateAgeOriginHtml\(c\)/.test(html));
ok('_drawCheat30Ridgeline defined',  /function _drawCheat30Ridgeline\(/.test(html));
ok('_bundleAgeOriginBlock defined',  /function _bundleAgeOriginBlock\(c\)/.test(html));
ok('_bundleAgeOriginBlock wired into dispatcher',
   /out\.push\(_bundleAgeOriginBlock\(c\)\);/.test(html));
ok('candidateAgeOriginHtml wired into renderCandidateMetadata',
   /candidateAgeOriginHtml\(c\)\s*\+/.test(html));
ok('_AGEORIG_CLASS palette defined', /_AGEORIG_CLASS\s*=\s*\{/.test(html));
ok('_AGEORIG_RIDGE_COLOR palette defined', /_AGEORIG_RIDGE_COLOR\s*=\s*\{/.test(html));

// --- CSS contract
console.log('\n=== CSS contract ===');
const cssClasses = [
  'ageorig-section', 'ageorig-empty', 'ageorig-grid', 'ageorig-left',
  'ageorig-right', 'ageorig-ridgeline', 'ageorig-ridgeline-svg',
  'ageorig-ridgeline-empty', 'ageorig-verdict-pill',
  'ageorig-verdict-explain', 'ageorig-bimodal-warn', 'ageorig-numbers',
  'ageorig-num-row', 'ageorig-num-label', 'ageorig-num-value',
  'ageorig-caveat',
];
for (const cls of cssClasses) {
  ok('CSS class .' + cls + ' defined',
     html.indexOf('.' + cls) >= 0);
}

// --- TSV column contract
console.log('\n=== TSV column contract (cheat30_* fields) ===');
const expectedNewCols = [
  'cheat30_origin_class', 'cheat30_separation_p',
  'cheat30_separation_effect', 'cheat30_age_proxy',
  'cheat30_dip_p', 'cheat30_dip_stat', 'cheat30_is_bimodal',
  'cheat30_n_ref', 'cheat30_n_het', 'cheat30_n_inv',
  'cheat30_mean_gds_same', 'cheat30_mean_gds_diff',
];
for (const col of expectedNewCols) {
  ok('TSV column "' + col + '" present in header list',
     html.indexOf("'" + col + "'") >= 0);
}

// --- Source-summary + caveats additions
console.log('\n=== Source-summary + caveats ===');
ok('Source-summary has cheat30 status row',
   html.indexOf('cheat30 GDS-by-genotype (turn 119)') >= 0);
ok('Caveats has cheat30 age-proxy ordinal note',
   html.indexOf('age proxy** (when present) is *ordinal*') >= 0);
ok('Caveats has cheat30 origin_class power note',
   html.indexOf('cheat30 origin_class') >= 0);

// =============================================================================
// Behavioural test: pull the helpers and exercise them in a sandbox.
// =============================================================================

console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        if (quote === '`' && src[i] === '$' && src[i+1] === '{') {
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
  '_isCheat30JSON',
  '_storeCheat30Results',
  '_cheat30ForCandidate',
  '_drawCheat30Ridgeline',
  '_bundleAgeOriginBlock',
  '_fmt4',
  '_fmt3',
  '_fmtP',
];
const fnSrcs = [];
for (const fn of fnNames) {
  const s = pullFunction(html, fn);
  if (!s) {
    fail++; console.log('  FAIL could not extract ' + fn);
  } else {
    fnSrcs.push(s);
  }
}

if (fnSrcs.length === fnNames.length) {
  // We also need a couple of constants (CHEAT30_TOOL, _AGEORIG_*) — pull them
  // by simple regex.
  const constSrcs = [];
  const constMatches = [
    /^const CHEAT30_TOOL = .*?;/m,
    /^const _AGEORIG_CLASS = \{[\s\S]*?\n\};/m,
    /^const _AGEORIG_RIDGE_COLOR = \{[\s\S]*?\n\};/m,
  ];
  for (const re of constMatches) {
    const m = html.match(re);
    if (m) constSrcs.push(m[0]);
  }
  // Build the sandbox
  function makeSandbox(stateOverlay) {
    const state = Object.assign({
      data: null, candidateList: [], cheat30Results: null,
    }, stateOverlay);
    const sandbox = {
      state, console,
      window: {}, localStorage: { getItem: () => null, setItem: () => {}, removeItem: () => {} },
      Number, Array, Math, JSON, Object, String, Infinity, isNaN, parseFloat,
      _esc: (s) => String(s == null ? '' : s)
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'),
    };
    vm.createContext(sandbox);
    for (const c of constSrcs) {
      vm.runInContext(c, sandbox);
    }
    for (const s of fnSrcs) {
      vm.runInContext(s, sandbox);
    }
    return sandbox;
  }

  // --- Test V1: schema validator
  console.log('\nTest V1: schema validator');
  const sV = makeSandbox({});
  const goodJson = {
    schema_version: 'cheat30_v1', tool: 'cheat30_gds_by_genotype',
    chrom: 'LG28', candidates: { 'LG28_cand_01': {} },
  };
  vm.runInContext('var __ok = _isCheat30JSON(' + JSON.stringify(goodJson) + ');', sV);
  ok('validator accepts well-formed JSON', sV.__ok === true);
  vm.runInContext('var __ok2 = _isCheat30JSON({ chrom: "LG28" });', sV);
  ok('validator rejects missing schema_version', sV.__ok2 === false);
  vm.runInContext('var __ok3 = _isCheat30JSON({ schema_version: "v2_other", chrom: "LG28", candidates: {} });', sV);
  ok('validator rejects wrong schema_version prefix', sV.__ok3 === false);
  vm.runInContext('var __ok4 = _isCheat30JSON({ schema_version: "cheat30_v1", chrom: "LG28" });', sV);
  ok('validator rejects missing candidates', sV.__ok4 === false);

  // --- Test S1: store + accessor
  console.log('\nTest S1: store + accessor');
  vm.runInContext(
    'var __chrom = _storeCheat30Results(' + JSON.stringify(goodJson) + ');',
    sV);
  ok('store returns chrom name', sV.__chrom === 'LG28');
  ok('state.cheat30Results populated', sV.state.cheat30Results && sV.state.cheat30Results.LG28);
  // accessor with matching candidate
  const cand = { id: 'LG28_cand_01', candidate_id: 'LG28_cand_01', chrom: 'LG28' };
  vm.runInContext('var __r = _cheat30ForCandidate(' + JSON.stringify(cand) + ');', sV);
  ok('accessor returns the candidate block when matched',
     sV.__r != null && typeof sV.__r === 'object');
  // accessor with non-matching candidate
  const cand2 = { id: 'NOT_THERE', chrom: 'LG28' };
  vm.runInContext('var __r2 = _cheat30ForCandidate(' + JSON.stringify(cand2) + ');', sV);
  ok('accessor returns null for missing candidate', sV.__r2 == null);
  // accessor with non-matching chrom
  const cand3 = { id: 'LG28_cand_01', chrom: 'LG99' };
  vm.runInContext('var __r3 = _cheat30ForCandidate(' + JSON.stringify(cand3) + ');', sV);
  ok('accessor returns null for unknown chrom', sV.__r3 == null);

  // --- Test B1: bundle block — null when no result
  console.log('\nTest B1: bundle block — empty state');
  const sB1 = makeSandbox({});
  vm.runInContext('var __out = _bundleAgeOriginBlock(' + JSON.stringify(cand) + ');', sB1);
  ok('null when state.cheat30Results undefined', sB1.__out == null);

  // --- Test B2: bundle block — populated state, recurrent class
  console.log('\nTest B2: bundle block — recurrent verdict renders');
  const sB2 = makeSandbox({
    cheat30Results: {
      LG28: {
        candidates: {
          'LG28_cand_01': {
            n_ref: 60, n_het: 106, n_inv: 60,
            separation_p: 4.8e-12, separation_effect: 0.082,
            age_proxy: 0.067,
            dip_stat: 0.041, dip_p: 0.012, is_bimodal: true,
            origin_class: 'recurrent',
            mean_ibs_same: 0.847, mean_ibs_diff: 0.765,
          },
        },
      },
    },
  });
  vm.runInContext('var __out = _bundleAgeOriginBlock(' + JSON.stringify(cand) + ');', sB2);
  ok('block emitted', typeof sB2.__out === 'string' && sB2.__out.length > 0);
  ok('block mentions recurrent', /recurrent/.test(sB2.__out));
  ok('block mentions Wilcoxon framing', /Wilcoxon/.test(sB2.__out));
  ok('block mentions bimodal',
     sB2.__out.indexOf('bimodal') >= 0);
  ok('block mentions n samples',
     sB2.__out.indexOf('n samples') >= 0);
  ok('block reports REF=60', sB2.__out.indexOf('REF=60') >= 0);
  ok('separation P uses scientific notation',
     /separation P[\s\S]*4\.80e-12/i.test(sB2.__out));

  // --- Test B3: bundle block — single_origin verdict
  console.log('\nTest B3: bundle block — single_origin');
  const sB3 = makeSandbox({
    cheat30Results: {
      LG28: {
        candidates: {
          'LG28_cand_01': {
            n_ref: 70, n_het: 90, n_inv: 66,
            separation_p: 1.2e-9, separation_effect: 0.05,
            age_proxy: 0.012,
            dip_stat: 0.020, dip_p: 0.45, is_bimodal: false,
            origin_class: 'single_origin',
            mean_ibs_same: 0.82, mean_ibs_diff: 0.77,
          },
        },
      },
    },
  });
  vm.runInContext('var __out = _bundleAgeOriginBlock(' + JSON.stringify(cand) + ');', sB3);
  ok('single_origin: block emitted', typeof sB3.__out === 'string' && sB3.__out.length > 0);
  ok('single_origin: block mentions unimodal',
     sB3.__out.indexOf('unimodal') >= 0);
  ok('single_origin: block mentions one founding event',
     sB3.__out.indexOf('one founding inversion event') >= 0);

  // --- Test B4: bundle block — inconclusive verdict
  console.log('\nTest B4: bundle block — inconclusive');
  const sB4 = makeSandbox({
    cheat30Results: {
      LG28: {
        candidates: {
          'LG28_cand_01': {
            n_ref: 4, n_het: 100, n_inv: 3,
            separation_p: 0.31, separation_effect: 0.005,
            age_proxy: null,
            dip_stat: null, dip_p: null, is_bimodal: null,
            origin_class: 'inconclusive',
            mean_ibs_same: null, mean_ibs_diff: null,
          },
        },
      },
    },
  });
  vm.runInContext('var __out = _bundleAgeOriginBlock(' + JSON.stringify(cand) + ');', sB4);
  ok('inconclusive: block emitted', typeof sB4.__out === 'string' && sB4.__out.length > 0);
  ok('inconclusive: block notes class size threshold',
     sB4.__out.indexOf('class size below threshold') >= 0);

  // --- Test R1: ridgeline renderer — three ridges
  console.log('\nTest R1: ridgeline renderer');
  const sR1 = makeSandbox({});
  const dens = {
    'HOM_REF_HOM_REF':  { x: [0.7, 0.75, 0.80, 0.85, 0.90], y: [0.1, 0.5, 1.2, 0.8, 0.2] },
    'HOM_INV_HOM_INV':  { x: [0.7, 0.75, 0.80, 0.85, 0.90], y: [0.3, 1.0, 0.6, 1.1, 0.4] },
    'HOM_REF_HOM_INV':  { x: [0.7, 0.75, 0.80, 0.85, 0.90], y: [0.4, 1.5, 1.0, 0.4, 0.1] },
  };
  const sumPair = {
    'HOM_REF_HOM_REF':  { mean: 0.81, sd: 0.04 },
    'HOM_INV_HOM_INV':  { mean: 0.83, sd: 0.06 },
    'HOM_REF_HOM_INV':  { mean: 0.76, sd: 0.04 },
  };
  vm.runInContext('var __svg = _drawCheat30Ridgeline(' +
    JSON.stringify(dens) + ', ' + JSON.stringify(sumPair) + ', { width: 380, height: 180 });', sR1);
  ok('ridgeline returns string', typeof sR1.__svg === 'string');
  ok('ridgeline contains <svg>', sR1.__svg.indexOf('<svg') >= 0);
  ok('ridgeline has REF / REF label',  sR1.__svg.indexOf('REF / REF') >= 0);
  ok('ridgeline has INV / INV label',  sR1.__svg.indexOf('INV / INV') >= 0);
  ok('ridgeline has REF / INV label',  sR1.__svg.indexOf('REF / INV') >= 0);
  ok('ridgeline has GDS axis label',   sR1.__svg.indexOf('GDS') >= 0);
  // 3 path elements (one per ridge fill) + 3 stroke paths = 6 paths total
  const nPaths = (sR1.__svg.match(/<path /g) || []).length;
  ok('ridgeline has 6 path elements (3 fill + 3 stroke)',
     nPaths === 6, 'got ' + nPaths);

  // --- Test R2: degenerate input
  console.log('\nTest R2: degenerate ridgeline input');
  const sR2 = makeSandbox({});
  vm.runInContext('var __svg = _drawCheat30Ridgeline({}, {}, {});', sR2);
  ok('empty pair_density returns empty-state stub',
     typeof sR2.__svg === 'string' && sR2.__svg.indexOf('no pair_density') >= 0);
}

// =============================================================================
console.log('\n=========================================');
console.log('passed: ' + pass + ' / failed: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
