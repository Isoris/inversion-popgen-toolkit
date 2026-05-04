// =============================================================================
// test_atlas_group_engine.js
// =============================================================================
// Tests in plain Node. Mocks the atlas's `state` global and `localStorage`,
// then exercises every code path in atlas_group_engine.js.
//
// Run:
//   node test_atlas_group_engine.js
// =============================================================================

'use strict';

// ----- localStorage shim ------------------------------------------------------
const _ls = new Map();
globalThis.localStorage = {
  getItem(k) { return _ls.has(k) ? _ls.get(k) : null; },
  setItem(k, v) { _ls.set(k, String(v)); },
  removeItem(k) { _ls.delete(k); },
  clear() { _ls.clear(); },
};

// ----- atlas state mock -------------------------------------------------------
// Build a small but realistic state: 12 samples, 1 chromosome with 100
// windows (1 Mb each → 100 Mb chrom), 1 focal candidate at 15-18 Mb with
// fish_calls populated, families assigned to half the samples.
function makeState() {
  const N = 12;
  const samples = [];
  for (let i = 0; i < N; i++) samples.push('CGA_' + String(i + 1).padStart(3, '0'));
  const wins = [];
  for (let w = 0; w < 100; w++) {
    wins.push({
      idx: w,
      center_mb: w + 0.5,
      start_bp: w * 1_000_000,
      end_bp: (w + 1) * 1_000_000,
      z: Math.abs(((w * 17) % 50) - 25) / 10,  // fake Z
    });
  }
  // K=3, regime int 0..2 for the 12 samples — simulate a real LG28-15Mb
  // candidate where samples 0..3 = H1/H1 (regime 0), 4..7 = H1/H2 (regime 1),
  // 8..11 = H2/H2 (regime 2).
  const fish_calls = [];
  for (let si = 0; si < N; si++) {
    let regime;
    if (si < 4) regime = 0;
    else if (si < 8) regime = 1;
    else regime = 2;
    fish_calls.push({
      sample_idx: si,
      sample_id:  samples[si],
      regime,
      confidence: 0.95,
      n_supporting: 3,
      n_intervals: 3,
      ambiguous: (si === 5),  // sample 5 marked ambiguous
      votes: [regime, regime, regime],
    });
  }
  // Sample 7 gets a different regime at L2 #1 — subband not stable
  fish_calls[7].n_supporting = 2;
  fish_calls[7].votes = [1, 1, 0];

  const candidate = {
    id: 'cand_LG28_15Mb',
    start_mb: 15.0,
    end_mb: 18.0,
    start_bp: 15_000_000,
    end_bp: 18_000_000,
    k: 3,
    fish_calls,
    l2_indices: [14, 16, 17],
    _system: 'detailed',  // detailed mode → H1/H1, H1/H2, H2/H2
  };

  // Second candidate on a different LG, for cross-LG predicates
  const candidate2 = {
    id: 'cand_LG14_22Mb',
    start_mb: 22.0,
    end_mb: 24.0,
    start_bp: 22_000_000,
    end_bp: 24_000_000,
    k: 3,
    fish_calls: fish_calls.map((fc, i) => ({
      ...fc,
      // At LG14 the partition is shuffled — only some samples are H2/H2 at both.
      regime: (i % 3),
    })),
    _system: 'detailed',
  };

  return {
    data: {
      chrom: 'C_gar_LG28',
      cohort_key: 'all_226_test',
      n_samples: N,
      samples,
      windows: wins,
      candidates: [candidate, candidate2],
      // family palette: samples 0..2 = family_1, 3..5 = family_2, 6..11 = no family
      family: [1, 1, 1, 2, 2, 2, null, null, null, null, null, null],
      // relatedness hubs: sample 5 in hub 30
      relatedness_hub: [null, null, null, null, null, 30, null, null, null, null, null, null],
      // F_ROH: sample 8 above threshold
      f_roh: [0.01, 0.02, 0.01, 0.0, 0.05, 0.03, 0.04, 0.01, 0.12, 0.05, 0.04, 0.02],
      f_roh_threshold: 0.10,
      // ancestry_sample minimal: maxQ_label[w][s] — give samples 0..5 mode K1, 6..11 mode K2
      ancestry_sample: {
        samples,
        maxQ_label: (function () {
          const m = [];
          for (let w = 0; w < 100; w++) {
            const row = [];
            for (let s = 0; s < N; s++) {
              row.push(s < 6 ? '1' : '2');
            }
            m.push(row);
          }
          return m;
        })(),
      },
    },
    cur: 50,
    tracked: [0, 4, 8],            // tracked = first sample of each regime
    activeSampleSet: null,
    candidate,                    // focal
    candidates: [candidate, candidate2],
    k: 3,
    labelVocab: 'detailed',
  };
}

globalThis.state = makeState();

// ----- import the engine -----------------------------------------------------
const engine = require('./atlas_group_engine.js');
// Force a re-load of localStorage now that we have it shimmed
engine.loadFromLocalStorage();
engine.invalidateDimensions();

// ----- test runner ------------------------------------------------------------
let _pass = 0, _fail = 0;
function ok(label, cond) {
  if (cond) { console.log('  ok  ' + label); _pass++; }
  else      { console.log('FAIL ' + label); _fail++; }
}
function eq(label, a, b) {
  const A = JSON.stringify(a), B = JSON.stringify(b);
  if (A === B) { console.log('  ok  ' + label); _pass++; }
  else         { console.log('FAIL ' + label + ': ' + A + ' != ' + B); _fail++; }
}
function group(name, fn) { console.log(name); fn(); }

// =============================================================================
// Tests
// =============================================================================

group('haplotype label parsing', () => {
  eq('parse H1/H2', engine.parseHaplotypePair('H1/H2'), ['H1', 'H2']);
  eq('parse H2/H2', engine.parseHaplotypePair('H2/H2'), ['H2', 'H2']);
  eq('parse H3/H2 sorts',  engine.parseHaplotypePair('H3/H2'), ['H2', 'H3']);
  eq('parse legacy g0',  engine.parseHaplotypePair('g0'), null);
  ok('carries H1 in H1/H2', engine.carriesHaplotype('H1/H2', 'H1'));
  ok('does NOT carry H3 in H1/H2', !engine.carriesHaplotype('H1/H2', 'H3'));
  ok('homozygous H2 in H2/H2', engine.isHomozygousFor('H2/H2', 'H2'));
  ok('NOT homozygous H1 in H1/H2', !engine.isHomozygousFor('H1/H2', 'H1'));
});

group('cohort + sample helpers', () => {
  eq('cohort key', engine.getCohortKey(), 'all_226_test');
  eq('n samples', engine.allSampleIds().length, 12);
  eq('idx 0 -> CGA_001', engine.sampleIdxToId(0), 'CGA_001');
  eq('CGA_005 -> idx 4', engine.sampleIdToIdx('CGA_005'), 4);
});

group('per-candidate diploid_class extractor (detailed mode)', () => {
  const dims = engine.listDimensions();
  const names = dims.map(d => d.name);
  ok('lists diploid_class@cand_LG28_15Mb',
    names.includes('diploid_class@cand_LG28_15Mb'));
  ok('lists carries_H1@cand_LG28_15Mb',
    names.includes('carries_H1@cand_LG28_15Mb'));
  ok('lists homozygous_H2@cand_LG28_15Mb',
    names.includes('homozygous_H2@cand_LG28_15Mb'));

  // Direct evaluator probe
  const ids_H1H1 = engine.evaluate({
    type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H1/H1'
  });
  eq('H1/H1 = first 4 samples', ids_H1H1, ['CGA_001','CGA_002','CGA_003','CGA_004']);

  const ids_H2H2 = engine.evaluate({
    type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2'
  });
  eq('H2/H2 = last 4 samples', ids_H2H2, ['CGA_009','CGA_010','CGA_011','CGA_012']);

  const carriesH1 = engine.evaluate({
    type: 'truthy', dim: 'carries_H1@cand_LG28_15Mb'
  });
  eq('carries H1 = samples 1..8 (H1/H1 ∪ H1/H2)', carriesH1.length, 8);

  const homH2 = engine.evaluate({
    type: 'truthy', dim: 'homozygous_H2@cand_LG28_15Mb'
  });
  eq('homozygous H2 = last 4', homH2, ['CGA_009','CGA_010','CGA_011','CGA_012']);

  const ambiguous = engine.evaluate({
    type: 'truthy', dim: 'fish_ambiguous@cand_LG28_15Mb'
  });
  eq('ambiguous = sample 6', ambiguous, ['CGA_006']);

  const stable = engine.evaluate({
    type: 'truthy', dim: 'subband_stable@cand_LG28_15Mb'
  });
  ok('subband_stable excludes sample 6 (ambig) and 8 (n_supp != n_int)',
    !stable.includes('CGA_006') && !stable.includes('CGA_008'));
});

group('global dimensions', () => {
  const fam1 = engine.evaluate({ type: 'eq', dim: 'family', value: 'family_1' });
  eq('family_1 = first 3', fam1, ['CGA_001','CGA_002','CGA_003']);

  const hub30 = engine.evaluate({ type: 'eq', dim: 'relatedness_hub', value: 'hub_30' });
  eq('hub_30 = sample 6', hub30, ['CGA_006']);

  const rohC = engine.evaluate({ type: 'truthy', dim: 'roh_carrier' });
  eq('roh carrier = sample 9 (f_roh=0.12 > 0.10)', rohC, ['CGA_009']);

  const tracked = engine.evaluate({ type: 'truthy', dim: 'tracked' });
  eq('tracked = samples 1, 5, 9 (idx 0, 4, 8)', tracked,
    ['CGA_001','CGA_005','CGA_009']);

  const k1 = engine.evaluate({ type: 'eq', dim: 'ancestry_dominantQ', value: 'K1' });
  eq('dominantQ K1 = first 6', k1.length, 6);
});

group('compound expressions', () => {
  // The Quentin example variant: (H1/H2 ∪ H2/H2) MINUS family_2
  const expr = {
    type: 'minus',
    left: {
      type: 'or',
      children: [
        { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H1/H2' },
        { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2' },
      ],
    },
    right: { type: 'eq', dim: 'family', value: 'family_2' },
  };
  ok('valid', engine.validateExpression(expr) === null);
  const ids = engine.evaluate(expr);
  // (samples 5..12) minus (samples 4,5,6) = samples 5..12 with 5,6 removed.
  // Sample 5 is in H1/H2 and family_2 → drop. Sample 6 is in H1/H2 and family_2 → drop.
  // Sample 4 is in H1/H1 (not in left), so MINUS doesn't drop it.
  // Result: 7,8,9,10,11,12 (= idx 6..11)
  eq('compound result', ids,
    ['CGA_007','CGA_008','CGA_009','CGA_010','CGA_011','CGA_012']);
});

group('cross-LG expression (the Result-3-row-4 question)', () => {
  // "Are samples that are H2/H2 at LG28 also H2/H2 at LG14?"
  const candIds = engine.evaluate({
    type: 'and',
    children: [
      { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2' },
      { type: 'eq', dim: 'diploid_class@cand_LG14_22Mb', value: 'H2/H2' },
    ],
  });
  // LG28 H2/H2 = samples 9..12 (regime 2)
  // LG14: regime = i % 3, so regime 2 → samples where i % 3 === 2 → idx 2,5,8,11
  // Intersection: idx 8,11 → CGA_009, CGA_012
  eq('cross-LG H2/H2 intersection', candIds, ['CGA_009','CGA_012']);
});

group('saved expressions + ref + cycle guard', () => {
  engine.saveExpression('hom2_lg28', {
    type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2'
  });
  engine.saveExpression('hom2_minus_fam2', {
    type: 'minus',
    left:  { type: 'ref', name: 'hom2_lg28' },
    right: { type: 'eq', dim: 'family', value: 'family_2' },
  });
  const ids = engine.evaluate({ type: 'ref', name: 'hom2_minus_fam2' });
  eq('ref expression resolves', ids,
    ['CGA_009','CGA_010','CGA_011','CGA_012']);
  // Cycle: A -> B -> A
  engine.saveExpression('cycA', { type: 'ref', name: 'cycB' });
  engine.saveExpression('cycB', { type: 'ref', name: 'cycA' });
  const cycRes = engine.evaluate({ type: 'ref', name: 'cycA' });
  eq('cycle returns empty', cycRes, []);
});

group('saved sets', () => {
  engine.saveSet('three_picked', ['CGA_002', 'CGA_005', 'CGA_011']);
  eq('list', engine.listSets(), ['three_picked']);
  const ids = engine.evaluate({ type: 'truthy', dim: 'saved_three_picked' });
  eq('saved set evaluates', ids, ['CGA_002','CGA_005','CGA_011']);

  // Freeze an expression to a saved set
  engine.freezeSetFromExpression('frozen_h2h2', {
    type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H2/H2'
  });
  eq('frozen set members', engine.getSet('frozen_h2h2'),
    ['CGA_009','CGA_010','CGA_011','CGA_012']);
});

group('slots + compileGroupsForRequest', () => {
  // Use the saved expressions and sets to build a request body.
  engine.setSlot(0, { name: 'HOM1', source: 'expression', ref: 'hom2_lg28' });
  engine.setSlot(1, { name: 'HET',  source: 'inline', ref: {
    type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H1/H2'
  }});
  engine.setSlot(2, { name: 'BG',   source: 'complement_of', ref: 'HOM1' });
  const compiled = engine.compileGroupsForRequest({ min_n: 1 });
  eq('three slots compiled', Object.keys(compiled.groups).sort(), ['BG','HET','HOM1']);
  eq('HOM1 size', compiled.groups.HOM1.length, 4);
  eq('HET size',  compiled.groups.HET.length, 4);
  eq('BG size = N - 4', compiled.groups.BG.length, 8);

  // Slot name validation
  let threw = false;
  try { engine.setSlot(3, { name: 'bad name', source: 'inline', ref: { type: 'all' } }); }
  catch (_) { threw = true; }
  ok('rejects bad slot name', threw);
});

group('cursor scope resolver', () => {
  eq('chrom scope = null', engine.cursorRegion('chrom'), null);
  const r1 = engine.cursorRegion('1w');
  eq('1w at cur=50 → window 50', r1, { start_bp: 50_000_000, end_bp: 51_000_000 });
  const r5 = engine.cursorRegion('5w');
  // 5 windows centered at cur=50, half = (5-1)/2 = 2 → windows 48..52
  // Window 48 starts at 48 Mb; window 52 ends at 53 Mb
  eq('5w → ±2 windows around cur=50', r5, { start_bp: 48_000_000, end_bp: 53_000_000 });
  const rL2 = engine.cursorRegion('L2');
  eq('L2 (no envelopes loaded) = null', rL2, null);
  const rC = engine.cursorRegion('candidate');
  eq('candidate scope', rC, { start_bp: 15_000_000, end_bp: 18_000_000 });
});

group('buildPopstatsRequest', () => {
  // Slots already set in earlier test.
  const req = engine.buildPopstatsRequest({ scope: 'candidate' });
  eq('request chrom',  req.chrom,  'C_gar_LG28');
  eq('request region', req.region, { start_bp: 15_000_000, end_bp: 18_000_000 });
  eq('request groups', Object.keys(req.groups).sort(), ['BG','HET','HOM1']);
  eq('request win_bp', req.win_bp, 50000);
  eq('request step_bp', req.step_bp, 10000);
  eq('warnings', req._warnings.length, 0);
});

group('lasso history + snapshots', () => {
  const lid1 = engine.recordSelection({
    sample_ids: ['CGA_001','CGA_002'], source_page: 'page1', kind: 'rect_lasso'
  });
  const lid2 = engine.recordSelection({
    sample_ids: ['CGA_005','CGA_006'], source_page: 'page1', kind: 'pca_lasso'
  });
  eq('history length', engine.listSelections().length, 2);
  // Lasso evaluator
  const r = engine.evaluate({ type: 'lasso', id: lid1 });
  eq('lasso evaluates', r, ['CGA_001','CGA_002']);
  // Pin + drop
  ok('pin', engine.pinSelection(lid1, true));
  ok('drop', engine.dropSelection(lid2));
  eq('after drop', engine.listSelections().length, 1);

  // Materialize lasso into saved set
  engine.lassoToSet(lid1, 'lassoed_pair');
  eq('materialized', engine.getSet('lassoed_pair'), ['CGA_001','CGA_002']);

  // Snapshot
  const sid = engine.saveSnapshot({ note: 'looking at LG28 hom2' });
  ok('snapshot listed', engine.listSnapshots().length === 1);
  // Mutate state and restore
  engine.deleteExpression('hom2_lg28');
  engine.restoreSnapshot(sid);
  ok('restored expression', engine.getExpression('hom2_lg28') !== null);
});

group('persistence round-trip', () => {
  // Verify localStorage entries exist and round-trip
  const keys = [];
  for (const k of _ls.keys()) if (k.startsWith('inversion_atlas.popgen.')) keys.push(k);
  ok('at least 3 popgen ls keys', keys.length >= 3);
  // Reload and verify expressions survive
  engine.loadFromLocalStorage();
  ok('hom2_minus_fam2 still present after reload',
    engine.getExpression('hom2_minus_fam2') !== null);
});

group('expression validation rejects malformed AST', () => {
  ok('null',     engine.validateExpression(null) !== null);
  ok('no type',  engine.validateExpression({}) !== null);
  ok('eq no dim', engine.validateExpression({ type: 'eq', value: 'x' }) !== null);
  ok('and empty', engine.validateExpression({ type: 'and', children: [] }) !== null);
  ok('unknown op', engine.validateExpression({ type: 'banana' }) !== null);
  // Depth bomb
  let deep = { type: 'all' };
  for (let i = 0; i < 25; i++) deep = { type: 'not', child: deep };
  ok('rejects too-deep', engine.validateExpression(deep) !== null);
});

group('slot relations: minus, union, intersect', () => {
  // Set up: slot A = HOM1 (idx 0..3), slot B = family_1 (idx 0..2)
  engine.setSlot(5, { name: 'A_hom1', source: 'inline',
    ref: { type: 'eq', dim: 'diploid_class@cand_LG28_15Mb', value: 'H1/H1' }});
  engine.setSlot(6, { name: 'B_fam1', source: 'inline',
    ref: { type: 'eq', dim: 'family', value: 'family_1' }});
  // C = A minus B → samples in HOM1 but NOT family_1
  engine.setSlot(7, { name: 'C_amb',  source: 'minus',
    ref: { left: 'A_hom1', right: 'B_fam1' }});
  const cIds = engine.resolveSlotByIdx(7);
  // HOM1 = idx 0..3 (CGA_001..004); family_1 = idx 0..2 (CGA_001..003)
  // C = CGA_004 only
  eq('A minus B = CGA_004', cIds, ['CGA_004']);

  // D = A union B  → HOM1 ∪ family_1 = CGA_001..004
  engine.setSlot(8, { name: 'D_un', source: 'union',
    ref: { members: ['A_hom1', 'B_fam1'] }});
  eq('A union B', engine.resolveSlotByIdx(8),
    ['CGA_001','CGA_002','CGA_003','CGA_004']);

  // E = A intersect B → HOM1 ∩ family_1 = CGA_001..003
  engine.setSlot(9, { name: 'E_in', source: 'intersect',
    ref: { members: ['A_hom1', 'B_fam1'] }});
  eq('A intersect B', engine.resolveSlotByIdx(9),
    ['CGA_001','CGA_002','CGA_003']);
});

group('slot cycle guard (minus pointing back)', () => {
  // F = G minus F (cycle) → should not infinite-loop
  engine.setSlot(0, { name: 'F', source: 'minus',
    ref: { left: 'G', right: 'F' }});
  engine.setSlot(1, { name: 'G', source: 'minus',
    ref: { left: 'F', right: 'F' }});
  // No throw, returns []
  ok('cycle resolves to empty without error',
    Array.isArray(engine.resolveSlotByIdx(0)));
});

// ----- summary ----------------------------------------------------------------
console.log('');
console.log('=========================================');
console.log(_pass + ' passed, ' + _fail + ' failed');
console.log('=========================================');
if (_fail > 0) process.exit(1);
