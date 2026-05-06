// =============================================================================
// tests/sv_evidence/test_step1_skeleton.js
// =============================================================================
// Step 1 regression test for atlas_sv_evidence.js. Verifies:
//   - Module loads cleanly under Node (UMD path).
//   - Public API surface matches the spec §4.3 contract.
//   - Constants (DEFAULT_FILTERS, FDR_COLOURS, PATTERN_LABELS, SV_TYPE_STYLE,
//     KARYOTYPE_REMAP, GROUP_COLOURS) are present and shaped correctly.
//   - Internal helpers (_fdrColour, _fmtBp, _fmtMb, _resolveKaryotypeGroups)
//     produce expected outputs.
//   - The synthetic fixture JSON validates against the spec §2.2 schema.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step1_skeleton.js
//
// This test does NOT exercise DOM rendering — that requires a browser
// environment. Step 2's test will introduce a jsdom-style harness for the
// table renderer.
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

const MODULE_PATH  = path.resolve(__dirname, '../../js/atlas_sv_evidence.js');
const FIXTURE_PATH = path.resolve(__dirname, 'fixture_sv_genotype_counts_v1.json');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  ✓ ' + msg); }
  else      { fail++; console.error('  ✗ ' + msg); }
}
function section(name) { console.log('\n— ' + name + ' —'); }

// ---------------------------------------------------------------------------
section('module load');
const mod = require(MODULE_PATH);
assert(mod && typeof mod === 'object', 'module exports an object');

// ---------------------------------------------------------------------------
section('public API surface (spec §4.3)');
['init', 'loadCandidate', 'refresh', 'setFilters', 'exportFilteredTSV', '_state']
  .forEach(k => assert(k in mod, '.' + k + ' is exported'));
assert(typeof mod.init                === 'function', '.init is a function');
assert(typeof mod.loadCandidate       === 'function', '.loadCandidate is a function');
assert(typeof mod.refresh             === 'function', '.refresh is a function');
assert(typeof mod.setFilters          === 'function', '.setFilters is a function');
assert(typeof mod.exportFilteredTSV   === 'function', '.exportFilteredTSV is a function');
assert(typeof mod._state              === 'object',   '._state is an object');

// ---------------------------------------------------------------------------
section('constants (spec §2.3, §3.3, §5.2)');
const C = mod._const;
assert(C, '_const is exposed');

// DEFAULT_FILTERS — spec §3.3
const f = C.DEFAULT_FILTERS;
assert(f.sv_type === 'All',                    'DEFAULT_FILTERS.sv_type defaults to All');
assert(f.quality === 'PASS',                   'DEFAULT_FILTERS.quality defaults to PASS');
assert(f.zone === 'boundary_500kb',            'DEFAULT_FILTERS.zone defaults to boundary_500kb');
assert(f.min_samples === 5,                    'DEFAULT_FILTERS.min_samples defaults to 5');
assert(f.show_only_associated === true,        'DEFAULT_FILTERS.show_only_associated defaults to true');
assert(f.fdr_threshold === 0.05,               'DEFAULT_FILTERS.fdr_threshold defaults to 0.05');

// FDR_COLOURS — spec §5.2 (5 tiers, order matters)
assert(Array.isArray(C.FDR_COLOURS) && C.FDR_COLOURS.length === 5,
       'FDR_COLOURS has 5 tiers');
assert(C.FDR_COLOURS[0].color === '#bf616a', 'FDR tier 1 (< 1e-6) is red #bf616a');
assert(C.FDR_COLOURS[1].color === '#d08770', 'FDR tier 2 (1e-6 – 1e-4) is orange #d08770');
assert(C.FDR_COLOURS[2].color === '#ebcb8b', 'FDR tier 3 (1e-4 – 0.01) is yellow #ebcb8b');
assert(C.FDR_COLOURS[3].color === '#a3be8c', 'FDR tier 4 (0.01 – 0.05) is green #a3be8c');
assert(C.FDR_COLOURS[4].color === '#888888', 'FDR tier 5 (> 0.05) is grey #888888');

// PATTERN_LABELS — spec §2.3 (6 labels with exact colours)
const pl = C.PATTERN_LABELS;
assert(pl.canonical_breakpoint_marker.color === '#7d4cdb', 'canonical_breakpoint_marker is violet #7d4cdb');
assert(pl.dominant_presence_marker.color    === '#d08770', 'dominant_presence_marker is amber #d08770');
assert(pl.het_specific_marker.color         === '#5e81ac', 'het_specific_marker is teal #5e81ac');
assert(pl.sub_haplotype_marker.color        === '#88c0d0', 'sub_haplotype_marker is blue #88c0d0');
assert(pl.internal_linked_marker.color      === '#ebcb8b', 'internal_linked_marker is yellow #ebcb8b');
assert(pl.uninformative.color               === '#888888', 'uninformative is grey #888888');

// SV_TYPE_STYLE — 5 entries
['BND','INV','DEL','DUP','Other'].forEach(t =>
  assert(t in C.SV_TYPE_STYLE, 'SV_TYPE_STYLE has ' + t));

// KARYOTYPE_REMAP — spec §5.1
assert(C.KARYOTYPE_REMAP.HOMO_1 === 'H1/H1', 'remap HOMO_1 → H1/H1');
assert(C.KARYOTYPE_REMAP.HET    === 'H1/H2', 'remap HET → H1/H2');
assert(C.KARYOTYPE_REMAP.HOMO_2 === 'H2/H2', 'remap HOMO_2 → H2/H2');

// ---------------------------------------------------------------------------
section('helpers');
const h = mod._internals;
assert(h._fdrColour(2.1e-9) === '#bf616a',  '_fdrColour(2.1e-9) → red');
assert(h._fdrColour(5e-5)   === '#d08770',  '_fdrColour(5e-5) → orange');
assert(h._fdrColour(0.005)  === '#ebcb8b',  '_fdrColour(0.005) → yellow');
assert(h._fdrColour(0.03)   === '#a3be8c',  '_fdrColour(0.03) → green');
assert(h._fdrColour(0.5)    === '#888888',  '_fdrColour(0.5) → grey');
assert(h._fdrColour(null)   === '#888888',  '_fdrColour(null) → grey');
assert(h._fdrColour(undefined) === '#888888', '_fdrColour(undefined) → grey');

assert(h._fmtMb(15142000) === '15.14 Mb', '_fmtMb formats Mb to 2 dp');
assert(h._fmtMb(null)     === '—',        '_fmtMb(null) → em dash');
assert(h._fmtBp(380)      === '380 bp',   '_fmtBp formats bp');
assert(h._fmtBp(15142)    === '15.1 kb',  '_fmtBp formats kb');
assert(h._fmtBp(15142000) === '15.14 Mb', '_fmtBp formats Mb');

// _resolveKaryotypeGroups under empty/missing global state
const noState = h._resolveKaryotypeGroups('cand_LG28_15Mb');
assert(noState.locked === false, '_resolveKaryotypeGroups → locked:false when no state');
assert(noState.groupOfSample instanceof Map, '_resolveKaryotypeGroups → returns a Map');

// _resolveKaryotypeGroups under a stubbed window.state
global.window = global.window || {};
global.window.state = {
  candidateState: {
    cand_LG28_15Mb: {
      locked_labels: [
        { sample_id: 's1', label: 'HOMO_1' },
        { sample_id: 's2', label: 'HET'    },
        { sample_id: 's3', label: 'HOMO_2' },
        { sample_id: 's4', label: 'HOMO_1' },
      ],
    },
  },
};
const stubbed = h._resolveKaryotypeGroups('cand_LG28_15Mb');
assert(stubbed.locked === true,                   '_resolveKaryotypeGroups → locked:true with stubbed state');
assert(stubbed.counts['H1/H1'] === 2,             '  H1/H1 count = 2');
assert(stubbed.counts['H1/H2'] === 1,             '  H1/H2 count = 1');
assert(stubbed.counts['H2/H2'] === 1,             '  H2/H2 count = 1');
assert(stubbed.groupOfSample.get('s1') === 'H1/H1', '  s1 → H1/H1');
assert(stubbed.groupOfSample.get('s2') === 'H1/H2', '  s2 → H1/H2');
delete global.window.state;

// ---------------------------------------------------------------------------
section('fixture validity (spec §2.2 schema)');
const raw = fs.readFileSync(FIXTURE_PATH, 'utf-8');
const fix = JSON.parse(raw);

assert(fix.format_version === 'sv_genotype_counts_v1', 'format_version is sv_genotype_counts_v1');
assert(typeof fix.candidate_id === 'string',           'candidate_id present');
assert(typeof fix.chrom === 'string',                  'chrom present');
assert(Number.isFinite(fix.boundary_left_bp),          'boundary_left_bp is numeric');
assert(Number.isFinite(fix.boundary_right_bp),         'boundary_right_bp is numeric');
assert(fix.boundary_right_bp > fix.boundary_left_bp,   'right > left');

// Zone definitions: 5 zones, all 2-tuples in ascending bp order
const Z = fix.zone_definitions_bp;
['left_flank','left_boundary','inversion_body','right_boundary','right_flank']
  .forEach(z => {
    assert(Array.isArray(Z[z]) && Z[z].length === 2, 'zone ' + z + ' is a 2-tuple');
    assert(Z[z][0] < Z[z][1], 'zone ' + z + ' is ascending');
  });

// Groups
['H1/H1','H1/H2','H2/H2'].forEach(g => {
  assert(g in fix.groups_used, 'groups_used has ' + g);
  assert(Number.isFinite(fix.groups_used[g].n), '  ' + g + '.n is numeric');
});

// SV calls — at least one of each type, all required fields present
assert(Array.isArray(fix.sv_calls) && fix.sv_calls.length >= 5, 'sv_calls has ≥ 5 entries');
const required = ['sv_id','sv_type','chrom','position_bp','zone','quality',
                  'genotype_counts','fisher','pattern_label'];
fix.sv_calls.forEach(call => {
  required.forEach(k => assert(k in call,
    'sv_call ' + call.sv_id + ' has ' + k));
  // Every group in genotype_counts must have AA/AB/BB/miss
  ['H1/H1','H1/H2','H2/H2'].forEach(g => {
    const gc = call.genotype_counts[g] || {};
    ['AA','AB','BB','miss'].forEach(c => assert(typeof gc[c] === 'number',
      '  ' + call.sv_id + '.' + g + '.' + c + ' is numeric'));
  });
  // pattern_label is one of the six known values
  assert(call.pattern_label in C.PATTERN_LABELS,
    '  ' + call.sv_id + ' has known pattern_label: ' + call.pattern_label);
});

// Boundary summary — left + right blocks
assert(fix.boundary_summary && fix.boundary_summary.left && fix.boundary_summary.right,
       'boundary_summary has left + right');
['BND','INV','DEL','DUP','Other'].forEach(t => {
  assert(t in fix.boundary_summary.left.by_sv_type,
    'boundary_summary.left.by_sv_type has ' + t);
});

// UpSet top combinations
assert(Array.isArray(fix.upset_top_combinations) && fix.upset_top_combinations.length > 0,
       'upset_top_combinations is populated');
fix.upset_top_combinations.forEach((c, i) => {
  assert(Array.isArray(c.members) && c.members.length >= 1,
    'upset[' + i + '] has members[]');
  assert(Number.isFinite(c.intersection_size),
    'upset[' + i + '].intersection_size is numeric');
});

// ---------------------------------------------------------------------------
section('exportFilteredTSV stub');
const tsv = mod.exportFilteredTSV();
assert(typeof tsv === 'string',                   'exportFilteredTSV returns a string');
assert(tsv.split('\n')[0].split('\t').length >= 5, 'exportFilteredTSV header has ≥ 5 columns');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
