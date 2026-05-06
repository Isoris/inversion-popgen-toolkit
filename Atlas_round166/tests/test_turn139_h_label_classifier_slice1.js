// =============================================================================
// turn 139 — H-label classifier Slice 1 (SPEC_observable_allele_h_label_system)
// =============================================================================
// Implements:
//   _classifyHLabelBands(candidate, opts?)
//   _invalidateHLabelClassificationCache()
//   _hlabelBuildRegimeHistogramTSV(opts?)
//   _hlabelExportRegimeHistogramTSV()
//   _hlabelDosageStateFromRate(rate)
//   _hlabelImpliedHFromHetCount(nHet)
//   _hlabelRegimeConsistency(nHom, impliedH)
//
// Constants (window-exposed):
//   _HLABEL_HET_FRACTION_HIGH = 0.70
//   _HLABEL_HET_FRACTION_LOW  = 0.30
//   _HLABEL_MIN_DOSAGE_N      = 5
//   _HLABEL_SAMPLE_HET_RATE_LO = 0.30
//   _HLABEL_SAMPLE_HET_RATE_HI = 0.70
//
// Slice 1 ships CLASSIFICATION + REGIME INFERENCE only — no H-labels yet
// (those come in Slice 2 after Quentin runs the calibration TSV on real
// data). Karyotype tab additions: per-band classification chip, regime
// summary line, calibration export button.
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Source-level — function definitions + constants
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_classifyHLabelBands defined',
   /function _classifyHLabelBands\b/.test(html));
ok('_invalidateHLabelClassificationCache defined',
   /function _invalidateHLabelClassificationCache\b/.test(html));
ok('_hlabelBuildRegimeHistogramTSV defined',
   /function _hlabelBuildRegimeHistogramTSV\b/.test(html));
ok('_hlabelExportRegimeHistogramTSV defined',
   /function _hlabelExportRegimeHistogramTSV\b/.test(html));
ok('_hlabelDosageStateFromRate defined',
   /function _hlabelDosageStateFromRate\b/.test(html));
ok('_hlabelImpliedHFromHetCount defined',
   /function _hlabelImpliedHFromHetCount\b/.test(html));
ok('_hlabelRegimeConsistency defined',
   /function _hlabelRegimeConsistency\b/.test(html));

ok('_HLABEL_HET_FRACTION_HIGH constant 0.70',
   /const _HLABEL_HET_FRACTION_HIGH\s*=\s*0\.70/.test(html));
ok('_HLABEL_HET_FRACTION_LOW constant 0.30',
   /const _HLABEL_HET_FRACTION_LOW\s*=\s*0\.30/.test(html));
ok('_HLABEL_MIN_DOSAGE_N constant 5',
   /const _HLABEL_MIN_DOSAGE_N\s*=\s*5/.test(html));
ok('_HLABEL_SAMPLE_HET_RATE_LO constant 0.30',
   /const _HLABEL_SAMPLE_HET_RATE_LO\s*=\s*0\.30/.test(html));
ok('_HLABEL_SAMPLE_HET_RATE_HI constant 0.70',
   /const _HLABEL_SAMPLE_HET_RATE_HI\s*=\s*0\.70/.test(html));

ok('window._classifyHLabelBands exported',
   /window\._classifyHLabelBands\b/.test(html));
ok('window._invalidateHLabelClassificationCache exported',
   /window\._invalidateHLabelClassificationCache\b/.test(html));
ok('window._hlabelBuildRegimeHistogramTSV exported',
   /window\._hlabelBuildRegimeHistogramTSV\b/.test(html));
ok('window._hlabelExportRegimeHistogramTSV exported',
   /window\._hlabelExportRegimeHistogramTSV\b/.test(html));

// ============================================================================
// 2. Helper: extract a function's source for sandboxing
// ============================================================================
function extractFnSrc(name) {
  const re = new RegExp('function ' + name + '\\([\\s\\S]*?\\n\\}', 'm');
  const m = html.match(re);
  return m ? m[0] : null;
}

// ============================================================================
// 3. _hlabelDosageStateFromRate — per-sample rate → state
// ============================================================================
console.log('\n=== 3. _hlabelDosageStateFromRate ===');

{
  const ctx = vm.createContext({});
  vm.runInContext(`
    const _HLABEL_SAMPLE_HET_RATE_LO = 0.30;
    const _HLABEL_SAMPLE_HET_RATE_HI = 0.70;
    ${extractFnSrc('_hlabelDosageStateFromRate')}
  `, ctx);
  ok('rate=NaN → NA',                ctx._hlabelDosageStateFromRate(NaN) === 'NA');
  ok('rate=null → NA',               ctx._hlabelDosageStateFromRate(null) === 'NA');
  ok('rate=undefined → NA',          ctx._hlabelDosageStateFromRate(undefined) === 'NA');
  ok('rate=Infinity → NA',           ctx._hlabelDosageStateFromRate(Infinity) === 'NA');
  ok('rate=0.0 → HOM',               ctx._hlabelDosageStateFromRate(0.0) === 'HOM');
  ok('rate=0.10 → HOM',              ctx._hlabelDosageStateFromRate(0.10) === 'HOM');
  ok('rate=0.25 → HOM (just below LO)', ctx._hlabelDosageStateFromRate(0.25) === 'HOM');
  ok('rate=0.30 → HET (boundary inclusive)', ctx._hlabelDosageStateFromRate(0.30) === 'HET');
  ok('rate=0.50 → HET',              ctx._hlabelDosageStateFromRate(0.50) === 'HET');
  ok('rate=0.70 → HET (boundary inclusive)', ctx._hlabelDosageStateFromRate(0.70) === 'HET');
  ok('rate=0.85 → HOM (above HI)',   ctx._hlabelDosageStateFromRate(0.85) === 'HOM');
  ok('rate=1.0 → HOM',               ctx._hlabelDosageStateFromRate(1.0) === 'HOM');
}

// ============================================================================
// 4. _hlabelImpliedHFromHetCount — n_het → implied H
// ============================================================================
console.log('\n=== 4. _hlabelImpliedHFromHetCount ===');

{
  const ctx = vm.createContext({ Number });
  vm.runInContext(extractFnSrc('_hlabelImpliedHFromHetCount'), ctx);
  ok('n_het=0 → H=1', ctx._hlabelImpliedHFromHetCount(0) === 1);
  ok('n_het=1 → H=2', ctx._hlabelImpliedHFromHetCount(1) === 2);
  ok('n_het=2 → H=3 (3 alleles, one missing cross)',
     ctx._hlabelImpliedHFromHetCount(2) === 3);
  ok('n_het=3 → H=3 (3 alleles, all crosses)',
     ctx._hlabelImpliedHFromHetCount(3) === 3);
  ok('n_het=4 → H=4 (suspect, beyond ceiling)',
     ctx._hlabelImpliedHFromHetCount(4) === 4);
  ok('n_het=6 → H=4 (still capped)',
     ctx._hlabelImpliedHFromHetCount(6) === 4);
  ok('n_het=-1 → H=1 (defensive)',
     ctx._hlabelImpliedHFromHetCount(-1) === 1);
  ok('n_het=null → H=1 (defensive)',
     ctx._hlabelImpliedHFromHetCount(null) === 1);
  ok('n_het=NaN → H=1 (defensive)',
     ctx._hlabelImpliedHFromHetCount(NaN) === 1);
}

// ============================================================================
// 5. _hlabelRegimeConsistency — (n_hom, implied_H) → outcome
// ============================================================================
console.log('\n=== 5. _hlabelRegimeConsistency ===');

{
  const ctx = vm.createContext({ Number });
  vm.runInContext(extractFnSrc('_hlabelRegimeConsistency'), ctx);
  ok('(2, 2) → clean',              ctx._hlabelRegimeConsistency(2, 2) === 'clean');
  ok('(3, 3) → clean',              ctx._hlabelRegimeConsistency(3, 3) === 'clean');
  ok('(2, 3) → partial_obs_minor',  ctx._hlabelRegimeConsistency(2, 3) === 'partial_obs_minor');
  ok('(1, 3) → partial_obs_major',  ctx._hlabelRegimeConsistency(1, 3) === 'partial_obs_major');
  ok('(0, 3) → suspect',            ctx._hlabelRegimeConsistency(0, 3) === 'suspect');
  ok('(4, 3) → suspect_extra_hom',  ctx._hlabelRegimeConsistency(4, 3) === 'suspect_extra_hom');
  ok('(5, 3) → suspect_extra_hom',  ctx._hlabelRegimeConsistency(5, 3) === 'suspect_extra_hom');
  ok('(null, 3) → unknown (defensive)',
     ctx._hlabelRegimeConsistency(null, 3) === 'unknown');
  ok('(2, null) → unknown (defensive)',
     ctx._hlabelRegimeConsistency(2, null) === 'unknown');
}

// ============================================================================
// 6. _classifyHLabelBands — full sandbox with mocked dosage rates
// ============================================================================
console.log('\n=== 6. _classifyHLabelBands sandbox ===');

// Build a sandbox that has the function and all its dependencies
function makeClassifierCtx(opts) {
  opts = opts || {};
  const stateStub = opts.state || {
    data: {
      chrom: 'LG28',
      n_samples: 0,
      samples: [],
      l2_envelopes: [],
    },
    candidate: null,
    candidateList: [],
    k: 3,
  };

  // Dosage rates supplied directly via mock _computeHetRateForL2.
  // Map: l2idx → Float32Array[n_samples] of rates.
  const rateMap = opts.rateMap || {};
  // PC1 values per L2 (sandbox), if needed.
  const pc1Map = opts.pc1Map || {};

  const ctx = {
    state: stateStub,
    window: undefined,
    Number, Array, Math, JSON, Float32Array, Int8Array, Int32Array, Set,
    isFinite, parseInt, parseFloat,
    console: { warn: () => {}, log: () => {} },
    setTimeout: () => {},
    // Mock the upstream het-rate function
    _computeHetRateForL2: function (l2idx) {
      if (rateMap[l2idx]) return rateMap[l2idx];
      // Default: all-NaN
      const arr = new Float32Array(stateStub.data.samples.length);
      for (let i = 0; i < arr.length; i++) arr[i] = NaN;
      return arr;
    },
    getL2Cluster: function (l2idx) {
      if (pc1Map[l2idx]) return { xs: pc1Map[l2idx] };
      return null;
    },
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  // Need all the constants and helpers
  vm.runInContext(`
    const _HLABEL_HET_FRACTION_HIGH = 0.70;
    const _HLABEL_HET_FRACTION_LOW = 0.30;
    const _HLABEL_MIN_DOSAGE_N = 5;
    const _HLABEL_SAMPLE_HET_RATE_LO = 0.30;
    const _HLABEL_SAMPLE_HET_RATE_HI = 0.70;
    ${extractFnSrc('_hlabelDosageStateFromRate')}
    ${extractFnSrc('_hlabelImpliedHFromHetCount')}
    ${extractFnSrc('_hlabelRegimeConsistency')}
    ${extractFnSrc('_classifyHLabelBands')}
  `, ctx);
  return { ctx, state: stateStub };
}

// 6a. Clean K=3, all 3 classes (HOM/HET/HOM), good agreement
{
  console.log('\n--- 6a. K=3 clean (HOM/HET/HOM pattern) ---');
  // 18 samples: 6 in band 0 (HOM_REF), 6 in band 1 (HET), 6 in band 2 (HOM_INV)
  // Band sizes >= _HLABEL_MIN_DOSAGE_N (=5) so classification proceeds.
  // Dosage rates: band 0 rates ~0.05 (HOM), band 1 rates ~0.50 (HET), band 2 rates ~0.05 (HOM)
  const rates = new Float32Array([
    0.04, 0.05, 0.06, 0.04, 0.05, 0.06,   // band 0 — 6 HOM
    0.50, 0.51, 0.49, 0.50, 0.51, 0.49,   // band 1 — 6 HET
    0.04, 0.05, 0.06, 0.04, 0.05, 0.06,   // band 2 — 6 HOM
  ]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 18,
        samples: Array.from({length: 18}, (_, i) => ({ cga: 'CGA' + i.toString().padStart(3, '0') })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 3,
    },
    rateMap: { 0: rates },
    pc1Map: { 0: new Float32Array([
      -2,-2,-2, -2,-2,-2,
       0, 0, 0,  0, 0, 0,
       2, 2, 2,  2, 2, 2,
    ])},
  });
  const cand = {
    id: 'cand_LG28_clean',
    K: 3,
    ref_l2: 0,
    locked_labels: new Int8Array([
      0,0,0, 0,0,0,
      1,1,1, 1,1,1,
      2,2,2, 2,2,2,
    ]),
  };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('returns object',                  cls && typeof cls === 'object');
  ok('K=3',                             cls.K === 3);
  ok('vocab_mode = observable_allele_v1', cls.vocab_mode === 'observable_allele_v1');
  ok('band_counts.n_hom = 2',           cls.band_counts.n_hom === 2);
  ok('band_counts.n_het = 1',           cls.band_counts.n_het === 1);
  ok('band_counts.n_ambiguous = 0',     cls.band_counts.n_ambiguous === 0);
  ok('band_counts.n_no_dosage = 0',     cls.band_counts.n_no_dosage === 0);
  ok('implied_regime.implied_H = 2',    cls.implied_regime.implied_H === 2);
  ok('implied_regime.n_inferred_systems = 1',
     cls.implied_regime.n_inferred_systems === 1);
  ok("implied_regime.consistency = 'clean'",
     cls.implied_regime.consistency === 'clean');
  ok('band 0 classification HOM',        cls.bands[0].classification === 'HOM');
  ok('band 1 classification HET',        cls.bands[1].classification === 'HET');
  ok('band 2 classification HOM',        cls.bands[2].classification === 'HOM');
  ok('band 0 het_fraction near 0',       cls.bands[0].het_fraction === 0);
  ok('band 1 het_fraction = 1.0',        cls.bands[1].het_fraction === 1);
  ok('band 0 median_pc1 = -2',           cls.bands[0].median_pc1 === -2);
  ok('band 1 median_pc1 = 0',            cls.bands[1].median_pc1 === 0);
  ok('band 2 median_pc1 = 2',            cls.bands[2].median_pc1 === 2);
  ok('summary.n_agree = 18 (all match)', cls.summary.n_agree === 18);
  ok('summary.n_mismatch = 0',           cls.summary.n_mismatch === 0);
  ok('summary.agree_fraction = 1.0',     cls.summary.agree_fraction === 1);
  ok('cached on candidate.h_classification',
     cand.h_classification === cls);
}

// 6b. K=4 partial obs: HOM/HET/HOM/HET pattern → 2 hom + 2 het → H=3 partial_obs_minor
{
  console.log('\n--- 6b. K=4 partial_obs_minor (HOM/HET/HOM/HET) ---');
  // 24 samples, 6 per band (above min_dosage_n=5)
  // band 0: HOM rates, band 1: HET rates, band 2: HOM rates, band 3: HET rates
  const rates = new Float32Array([
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 0 — HOM
    0.50, 0.50, 0.50, 0.50, 0.50, 0.50,   // band 1 — HET
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 2 — HOM
    0.50, 0.50, 0.50, 0.50, 0.50, 0.50,   // band 3 — HET
  ]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 24,
        samples: Array.from({length: 24}, (_, i) => ({ cga: 'CGA' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 4,
    },
    rateMap: { 0: rates },
  });
  const cand = {
    id: 'cand_partial', K: 4, ref_l2: 0,
    locked_labels: new Int8Array([
      0,0,0, 0,0,0,
      1,1,1, 1,1,1,
      2,2,2, 2,2,2,
      3,3,3, 3,3,3,
    ]),
  };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('n_hom = 2',                       cls.band_counts.n_hom === 2);
  ok('n_het = 2',                       cls.band_counts.n_het === 2);
  ok('implied_H = 3',                   cls.implied_regime.implied_H === 3);
  ok('n_inferred_systems = 2',          cls.implied_regime.n_inferred_systems === 2);
  ok("consistency = 'partial_obs_minor'",
     cls.implied_regime.consistency === 'partial_obs_minor');
}

// 6c. Suspect_extra_hom: 3 hom + 1 het at K=4
{
  console.log('\n--- 6c. K=4 suspect_extra_hom (3 hom + 1 het) ---');
  const rates = new Float32Array([
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 0 — HOM
    0.50, 0.50, 0.50, 0.50, 0.50, 0.50,   // band 1 — HET
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 2 — HOM
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 3 — HOM
  ]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 24,
        samples: Array.from({length: 24}, (_, i) => ({ cga: 'CGA' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 4,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_X', K: 4, ref_l2: 0,
    locked_labels: new Int8Array([
      0,0,0, 0,0,0,
      1,1,1, 1,1,1,
      2,2,2, 2,2,2,
      3,3,3, 3,3,3,
    ]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('n_hom = 3',                       cls.band_counts.n_hom === 3);
  ok('n_het = 1',                       cls.band_counts.n_het === 1);
  ok('implied_H from n_het=1 is 2',     cls.implied_regime.implied_H === 2);
  ok("consistency = 'suspect_extra_hom' (3 hom > 2 implied)",
     cls.implied_regime.consistency === 'suspect_extra_hom');
}

// 6d. Ambiguous band — fraction between thresholds
{
  console.log('\n--- 6d. Ambiguous band (het_fraction = 0.5) ---');
  // 6 samples in band 0: 3 HET (rate 0.5), 3 HOM (rate 0.05) → het_fraction = 0.5
  const rates = new Float32Array([0.50, 0.50, 0.50, 0.05, 0.05, 0.05]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 6,
        samples: Array.from({length: 6}, (_, i) => ({ cga: 'CGA' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 1,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_ambig', K: 1, ref_l2: 0,
    locked_labels: new Int8Array([0, 0, 0, 0, 0, 0]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('1 ambiguous band',                cls.band_counts.n_ambiguous === 1);
  ok('band 0 classification = AMBIGUOUS',
     cls.bands[0].classification === 'AMBIGUOUS');
  ok('het_fraction = 0.5',              cls.bands[0].het_fraction === 0.5);
}

// 6e. NO_DOSAGE — too few non-NA samples (below min_n)
{
  console.log('\n--- 6e. NO_DOSAGE (below min_n=5) ---');
  // 4 samples in band 0, all dosage HET → met threshold but n<5 → NO_DOSAGE
  const rates = new Float32Array([0.50, 0.50, 0.50, 0.50]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 4,
        samples: Array.from({length: 4}, (_, i) => ({ cga: 'CGA' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 1,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_smol', K: 1, ref_l2: 0,
    locked_labels: new Int8Array([0, 0, 0, 0]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('classification = NO_DOSAGE (n<min)',
     cls.bands[0].classification === 'NO_DOSAGE');
  ok('n_no_dosage count = 1',           cls.band_counts.n_no_dosage === 1);
  // het_fraction is still computed (4/4 = 1.0) but classification respects min_n
  ok('het_fraction still computed',     cls.bands[0].het_fraction === 1);
}

// 6f. NO_DOSAGE — all-NaN dosage rates
{
  console.log('\n--- 6f. NO_DOSAGE (all NaN rates) ---');
  const rates = new Float32Array([NaN, NaN, NaN, NaN, NaN, NaN]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 6,
        samples: Array.from({length: 6}, (_, i) => ({ cga: 'CGA' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 1,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_nan', K: 1, ref_l2: 0,
    locked_labels: new Int8Array([0, 0, 0, 0, 0, 0]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('classification = NO_DOSAGE (all NaN)',
     cls.bands[0].classification === 'NO_DOSAGE');
  ok('het_fraction = NaN',              isNaN(cls.bands[0].het_fraction));
  ok('n_dosage_non_na = 0',             cls.bands[0].n_dosage_non_na === 0);
}

// 6g. Per-sample verdicts — MISMATCH cases
{
  console.log('\n--- 6g. Per-sample verdicts ---');
  // 18 samples (6 per band, above min_n=5).
  // Band 0: 5 HOM + 1 HET (rate 0.50) → het_frac 0.167 → HOM
  //   sample 5 is the dosage-HET one → MISMATCH_FALSE_HOM
  // Band 1: 5 HET + 1 HOM (rate 0.05) → het_frac 0.833 → HET
  //   sample 11 is the dosage-HOM one → MISMATCH_FALSE_HET
  // Band 2: 6 HOM → clean HOM, no mismatches
  const rates = new Float32Array([
    0.05, 0.05, 0.05, 0.05, 0.05, 0.50,   // band 0 (HOM): sample 5 is HET-by-dosage
    0.50, 0.50, 0.50, 0.50, 0.50, 0.05,   // band 1 (HET): sample 11 is HOM-by-dosage
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   // band 2 (HOM): all HOM
  ]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 18,
        samples: Array.from({length: 18}, (_, i) => ({ cga: 'CGA' + i.toString().padStart(3, '0') })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 3,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_mix', K: 3, ref_l2: 0,
    locked_labels: new Int8Array([
      0,0,0, 0,0,0,
      1,1,1, 1,1,1,
      2,2,2, 2,2,2,
    ]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('band 0 classification HOM (1/6 = 0.167)',
     cls.bands[0].classification === 'HOM');
  ok('band 1 classification HET (5/6 = 0.833)',
     cls.bands[1].classification === 'HET');
  ok('band 2 classification HOM (0/6 = 0.000)',
     cls.bands[2].classification === 'HOM');
  ok('per_sample_verdicts has 18 entries',
     cls.per_sample_verdicts.length === 18);

  // Find sample 5's verdict (HOM band, dosage HET → MISMATCH_FALSE_HOM)
  const v5 = cls.per_sample_verdicts.find(v => v.sample_idx === 5);
  ok('sample 5 verdict = MISMATCH_FALSE_HOM',
     v5 && v5.verdict === 'MISMATCH_FALSE_HOM');
  ok('sample 5 dosage = HET',          v5 && v5.dosage === 'HET');

  // Find sample 11's verdict (HET band, dosage HOM → MISMATCH_FALSE_HET)
  const v11 = cls.per_sample_verdicts.find(v => v.sample_idx === 11);
  ok('sample 11 verdict = MISMATCH_FALSE_HET',
     v11 && v11.verdict === 'MISMATCH_FALSE_HET');

  ok('summary.n_agree = 16', cls.summary.n_agree === 16);
  ok('summary.n_mismatch = 2', cls.summary.n_mismatch === 2);
}

// 6h. Per-sample verdicts with larger band populations
{
  console.log('\n--- 6h. Per-sample verdicts (large bands) ---');
  // 24 samples in K=3 — 8 per band. Band 0 mostly HOM (1/8 het = 0.125),
  // band 1 mostly HET (7/8 = 0.875), band 2 mostly HOM (1/8 = 0.125).
  const rates = new Float32Array([
    0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.50,   // band 0 → HOM (sample 7 is mismatch)
    0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.05,   // band 1 → HET (sample 15 is mismatch)
    0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,   // band 2 → HOM (all match)
  ]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 24,
        samples: Array.from({length: 24}, (_, i) => ({ cga: 'CGA' + i.toString().padStart(3, '0') })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 3,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_big', K: 3, ref_l2: 0,
    locked_labels: new Int8Array([
      0,0,0,0,0,0,0,0,
      1,1,1,1,1,1,1,1,
      2,2,2,2,2,2,2,2,
    ])};
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('band 0 → HOM (het_frac=0.125)',  cls.bands[0].classification === 'HOM');
  ok('band 1 → HET (het_frac=0.875)',  cls.bands[1].classification === 'HET');
  ok('band 2 → HOM (het_frac=0.0)',    cls.bands[2].classification === 'HOM');
  ok('per_sample_verdicts has 24 entries',
     cls.per_sample_verdicts.length === 24);

  // Find sample 7's verdict (HOM band, dosage HET → MISMATCH_FALSE_HOM)
  const v7 = cls.per_sample_verdicts.find(v => v.sample_idx === 7);
  ok('sample 7 verdict = MISMATCH_FALSE_HOM',
     v7 && v7.verdict === 'MISMATCH_FALSE_HOM');
  ok('sample 7 dosage = HET',          v7 && v7.dosage === 'HET');
  ok('sample 7 sample_id populated from .cga',
     v7 && v7.sample_id === 'CGA007');

  // Find sample 15's verdict (HET band, dosage HOM → MISMATCH_FALSE_HET)
  const v15 = cls.per_sample_verdicts.find(v => v.sample_idx === 15);
  ok('sample 15 verdict = MISMATCH_FALSE_HET',
     v15 && v15.verdict === 'MISMATCH_FALSE_HET');

  ok('summary.n_agree = 22',           cls.summary.n_agree === 22);
  ok('summary.n_mismatch = 2',         cls.summary.n_mismatch === 2);
  ok('summary.agree_fraction ≈ 22/24', Math.abs(cls.summary.agree_fraction - 22/24) < 1e-6);
}

// 6i. Cache hit — second call returns cached
{
  console.log('\n--- 6i. Caching ---');
  const rates = new Float32Array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]);
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 6,
        samples: Array.from({length: 6}, (_, i) => ({ cga: 'C' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 1,
    },
    rateMap: { 0: rates },
  });
  const cand = { id: 'cand_cache', K: 1, ref_l2: 0,
    locked_labels: new Int8Array([0, 0, 0, 0, 0, 0]) };
  const cls1 = sb.ctx._classifyHLabelBands(cand);
  ok('first call returns object', cls1 != null);
  ok('cached on candidate', cand.h_classification === cls1);
  // Mutate the cache — second call should return the mutated cache
  cand.h_classification.K = 999;
  const cls2 = sb.ctx._classifyHLabelBands(cand);
  ok('second call returns same cached object',
     cls2.K === 999, 'cache wasn\'t hit');
  // Force=true bypasses cache
  const cls3 = sb.ctx._classifyHLabelBands(cand, { force: true });
  ok('force:true recomputes (K back to 1)',
     cls3.K === 1);
}

// 6j. Null candidate / missing locked_labels → null
{
  console.log('\n--- 6j. Null inputs ---');
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 0, samples: [], l2_envelopes: [],
      },
      candidate: null, candidateList: [], k: 3,
    },
  });
  ok('null candidate returns null',
     sb.ctx._classifyHLabelBands(null) === null);
  ok('missing locked_labels returns null',
     sb.ctx._classifyHLabelBands({ K: 3 }) === null);
  ok('empty samples returns null',
     sb.ctx._classifyHLabelBands({ K: 3, locked_labels: new Int8Array([]) }) === null);
}

// 6k. No dosage rates available (no _computeHetRateForL2 → all NA) — but candidate is valid
{
  console.log('\n--- 6k. No dosage layer (graceful fallback) ---');
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 6,
        samples: Array.from({length: 6}, (_, i) => ({ cga: 'C' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null, candidateList: [], k: 1,
    },
    // No rateMap — sandbox returns all-NaN
  });
  const cand = { id: 'cand_nodose', K: 1, ref_l2: 0,
    locked_labels: new Int8Array([0, 0, 0, 0, 0, 0]) };
  const cls = sb.ctx._classifyHLabelBands(cand);
  ok('returns object even without dosage data',
     cls && typeof cls === 'object');
  ok('all bands NO_DOSAGE',             cls.band_counts.n_no_dosage === 1);
  ok('no per-sample verdicts',          cls.per_sample_verdicts.length === 0);
  ok('summary.agree_fraction is NaN',   isNaN(cls.summary.agree_fraction));
}

// ============================================================================
// 7. _hlabelBuildRegimeHistogramTSV — TSV row construction
// ============================================================================
console.log('\n=== 7. _hlabelBuildRegimeHistogramTSV ===');

{
  // Build a sandbox with multiple candidates pre-classified
  const sb = makeClassifierCtx({
    state: {
      data: {
        chrom: 'LG28', n_samples: 9,
        samples: Array.from({length: 9}, (_, i) => ({ cga: 'C' + i })),
        l2_envelopes: [{ start_bp: 0, end_bp: 1e6 }],
      },
      candidate: null,
      candidateList: [
        // Candidate 1 — pre-classified
        {
          id: 'cand_a', K: 3, ref_l2: 0,
          locked_labels: new Int8Array([0,0,0, 1,1,1, 2,2,2]),
          h_classification: {
            K: 3,
            band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
            implied_regime: { implied_H: 2, n_inferred_systems: 1, consistency: 'clean' },
            summary: { n_agree: 9, n_mismatch: 0, agree_fraction: 1.0 },
          },
        },
        // Candidate 2 — pre-classified (different regime)
        {
          id: 'cand_b', K: 4, ref_l2: 0,
          locked_labels: new Int8Array([0,0,0, 1,1,1, 2,2,2, 3,3,3]),
          h_classification: {
            K: 4,
            band_counts: { n_hom: 2, n_het: 2, n_ambiguous: 0, n_no_dosage: 0 },
            implied_regime: { implied_H: 3, n_inferred_systems: 2, consistency: 'partial_obs_minor' },
            summary: { n_agree: 11, n_mismatch: 1, agree_fraction: 11/12 },
          },
        },
        // Candidate 3 — NOT pre-classified, has runMissing flag check
        {
          id: 'cand_c', K: 3, ref_l2: 0,
          locked_labels: new Int8Array([0,0,0, 1,1,1, 2,2,2]),
        },
      ],
      k: 3,
    },
    rateMap: {
      0: new Float32Array([0.05,0.05,0.05, 0.50,0.50,0.50, 0.05,0.05,0.05]),
    },
  });
  vm.runInContext(extractFnSrc('_hlabelBuildRegimeHistogramTSV'), sb.ctx);

  // Without runMissing: only 2 rows (cand_a + cand_b)
  const tsv1 = sb.ctx._hlabelBuildRegimeHistogramTSV();
  const lines1 = tsv1.trim().split('\n');
  ok('header line present',                lines1[0].startsWith('candidate_id\tK\tn_hom'));
  ok('without runMissing: 1 header + 2 rows = 3 lines',
     lines1.length === 3);
  ok('cand_a row contains expected fields',
     lines1[1].startsWith('cand_a\t3\t2\t1\t0\t0\t2\t1\tclean\t1.0000\t9'));
  ok('cand_b row contains expected fields',
     lines1[2].startsWith('cand_b\t4\t2\t2\t0\t0\t3\t2\tpartial_obs_minor\t'));

  // With runMissing: 3 rows (cand_c gets classified)
  const tsv2 = sb.ctx._hlabelBuildRegimeHistogramTSV({ runMissing: true });
  const lines2 = tsv2.trim().split('\n');
  ok('with runMissing: 1 header + 3 rows = 4 lines',
     lines2.length === 4);
  ok('cand_c classified on the fly', /^cand_c\t3\t/.test(lines2[3]));
}

// ============================================================================
// 8. _hlabelExportRegimeHistogramTSV — sandbox returns false, real Blob returns true
// ============================================================================
console.log('\n=== 8. _hlabelExportRegimeHistogramTSV ===');

// 8a. Sandbox (no Blob) → false
{
  const sb = makeClassifierCtx({
    state: {
      data: { chrom: 'LG28', samples: [{cga:'X'}], l2_envelopes: [], n_samples: 1 },
      candidate: null,
      candidateList: [{ id: 'c', K: 1, ref_l2: 0, locked_labels: new Int8Array([0]) }],
      k: 1,
    },
  });
  vm.runInContext(extractFnSrc('_hlabelBuildRegimeHistogramTSV'), sb.ctx);
  vm.runInContext(extractFnSrc('_hlabelExportRegimeHistogramTSV'), sb.ctx);
  // Override: Blob undefined in this sandbox
  ok('sandbox returns false (no Blob)',
     sb.ctx._hlabelExportRegimeHistogramTSV() === false);
}

// 8b. With Blob mock → true, correct filename pattern
{
  const downloads = [];
  const stateStub = {
    data: { chrom: 'LG28', samples: [{cga:'A'},{cga:'B'},{cga:'C'},{cga:'D'},{cga:'E'},{cga:'F'}], l2_envelopes: [], n_samples: 6 },
    candidate: null,
    candidateList: [{
      id: 'cand_X', K: 1, ref_l2: 0,
      locked_labels: new Int8Array([0,0,0,0,0,0]),
      h_classification: {
        K: 1,
        band_counts: { n_hom: 1, n_het: 0, n_ambiguous: 0, n_no_dosage: 0 },
        implied_regime: { implied_H: 1, n_inferred_systems: 0, consistency: 'clean' },
        summary: { n_agree: 6, n_mismatch: 0, agree_fraction: 1.0 },
      },
    }],
    k: 1,
  };
  const ctx = {
    state: stateStub,
    window: undefined,
    Number, Array, Math, Float32Array, Int8Array,
    isFinite,
    console: { warn: () => {} },
    setTimeout: () => {},
    Blob: function (parts, opts) { return { _parts: parts, _opts: opts }; },
    URL: { createObjectURL: () => 'blob://x', revokeObjectURL: () => {} },
    document: {
      body: { appendChild() {}, removeChild() {} },
      createElement(tag) {
        const el = { tagName: tag };
        el.click = function () {
          downloads.push({ download: el.download, href: el._href });
        };
        Object.defineProperty(el, 'href', {
          set(v) { el._href = v; }, get() { return el._href; },
        });
        return el;
      },
    },
    _classifyHLabelBands: () => null,   // stub; cand has cached
  };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(extractFnSrc('_hlabelBuildRegimeHistogramTSV'), ctx);
  vm.runInContext(extractFnSrc('_hlabelExportRegimeHistogramTSV'), ctx);
  const result = ctx._hlabelExportRegimeHistogramTSV();
  ok('with Blob mock: returns true', result === true);
  ok('one download triggered',       downloads.length === 1);
  ok('filename pattern h_label_regime_histogram.<chrom>.tsv',
     /^h_label_regime_histogram\.LG28\.tsv$/.test(downloads[0].download || ''));
}

// ============================================================================
// 9. _invalidateHLabelClassificationCache
// ============================================================================
console.log('\n=== 9. _invalidateHLabelClassificationCache ===');

{
  const stateStub = {
    candidate: { h_classification: { K: 3, _stale: true } },
    candidateList: [
      { h_classification: { K: 4, _stale: true } },
      { /* no h_classification */ },
      { h_classification: { K: 6, _stale: true } },
    ],
  };
  const ctx = { state: stateStub, window: undefined, console: { warn: () => {} } };
  ctx.window = ctx;
  vm.createContext(ctx);
  vm.runInContext(extractFnSrc('_invalidateHLabelClassificationCache'), ctx);
  ctx._invalidateHLabelClassificationCache();
  ok('focused candidate cache cleared',
     stateStub.candidate.h_classification === null);
  ok('list candidate 0 cache cleared',
     stateStub.candidateList[0].h_classification === null);
  ok('list candidate 1 untouched (no cache present)',
     stateStub.candidateList[1].h_classification === undefined);
  ok('list candidate 2 cache cleared',
     stateStub.candidateList[2].h_classification === null);
}

// ============================================================================
// 10. G-panel karyotype tab — Slice 1 additions present
// ============================================================================
console.log('\n=== 10. Karyotype tab additions (source-level) ===');

ok('renderer calls _classifyHLabelBands',
   /_classifyHLabelBands\(c\)/.test(html));
ok('renderer surfaces "Implied regime" header',
   /Implied regime · empirical-first classifier/.test(html));
ok('renderer shows "no dosage" hint when classifier returns null',
   /pending dosage_chunks layer/.test(html));
ok('per-band rows have classification chip column (6 columns now)',
   /grid-template-columns:\s*18px 80px 60px 90px 1fr auto/.test(html));
ok('chip styles for HET-class',
   /HET-class band: dosage-HET fraction/.test(html));
ok('chip styles for HOM-class',
   /HOM-class band: dosage-HET fraction/.test(html));
ok('chip styles for AMBIGUOUS',
   /AMBIGUOUS: dosage-HET fraction/.test(html));
ok('chip styles for NO_DOSAGE',
   /NO_DOSAGE.*non-NA samples|NO_DOSAGE.*not loaded/.test(html));
ok('calibration export button present',
   /id="gpHlabelCalibrationBtn"/.test(html));
ok('calibration button label includes "calibration TSV"',
   /calibration TSV \(regime histogram\)/.test(html));
ok('calibration button wired in renderer',
   /gpHlabelCalibrationBtn[\s\S]{0,500}_hlabelExportRegimeHistogramTSV/.test(html));

// ============================================================================
// 11. Cache invalidation hooked to dosage rotation
// ============================================================================
console.log('\n=== 11. Cache invalidation hook ===');

ok('_invalidateHetRateCache calls _invalidateHLabelClassificationCache',
   /_invalidateHetRateCache[\s\S]*?_invalidateHLabelClassificationCache/.test(html));

// ============================================================================
// Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
