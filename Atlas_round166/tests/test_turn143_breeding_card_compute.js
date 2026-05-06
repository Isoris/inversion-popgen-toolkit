// =============================================================================
// turn 143 — Breeding-readiness card, Turn B (computation layer)
// =============================================================================
// SPEC: specs_new_turn131/SPEC_per_candidate_breeding_readiness_card.md
//
// Builds on Turn A's cohort_diversity loader (turn 142). Adds five
// pure helpers + one builder:
//   - _breedingCardKaryotypePerSample(cand)  — REF/HET/INV from h_class + PC1
//   - _wilcoxonRankSumP(a, b)                 — two-sided Mann–Whitney U
//                                              (tie + continuity correction)
//   - _summarizeFROHGroup(values)             — n / mean / median / sd / IQR
//   - _perArrangementFROH(cand)               — F_ROH stats by karyotype
//   - _arrangementMAF(counts)                 — minor allele frequency
//   - _carrierByK8Table(cand)                 — K8 × karyotype contingency
//   - _generatePairingAdvice(card)            — SPEC §2 rule engine
//   - _buildBreedingCard(cand)                — pure data builder
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
function approx(a, b, tol) { return Math.abs(a - b) < (tol || 1e-6); }

// ============================================================================
// Sandbox builder
// ============================================================================
function buildSandbox() {
  const re = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?window\._buildBreedingCard = _buildBreedingCard;\s*\n}/;
  const m = html.match(re);
  if (!m) throw new Error('breeding-card region not found');
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  if (!normal) throw new Error('normalCDF region not found');
  const _ls = new Map();
  const ctx = vm.createContext({
    state: { data: null, cohortDiversity: null },
    window: {},
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
    console,
  });
  vm.runInContext(normal[0] + '\n' + m[0], ctx);
  ctx.__ls = _ls;
  return ctx;
}

// Standard K=3 candidate fixture used across many tests.
function k3Fixture(refF, hetF, invF) {
  refF = refF || [0.20, 0.21, 0.19];
  hetF = hetF || [0.24, 0.26, 0.25, 0.27, 0.23];
  invF = invF || [0.31, 0.33, 0.32];
  const samples = [];
  const labels = [];
  let id = 1;
  for (const f of refF) { samples.push({ sample_id: 'CGA' + String(id).padStart(3,'0'), k8: 'K' + ((id%3)+1), f_roh: f }); labels.push(0); id++; }
  for (const f of hetF) { samples.push({ sample_id: 'CGA' + String(id).padStart(3,'0'), k8: 'K' + ((id%3)+1), f_roh: f }); labels.push(1); id++; }
  for (const f of invF) { samples.push({ sample_id: 'CGA' + String(id).padStart(3,'0'), k8: 'K' + ((id%3)+1), f_roh: f }); labels.push(2); id++; }
  const chromSamples = samples.map(s => ({ cga: s.sample_id }));
  const cand = {
    id: 'LG28_test',
    chrom: 'LG28',
    start_bp: 14000000,
    end_bp:   16500000,
    K: 3,
    tier: 'Tier 1',
    confidence: 'high',
    locked_labels: new Int8Array(labels),
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2.0 },
        { band_idx: 1, classification: 'HET', median_pc1:  0.0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  2.0 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'clean' },
    },
  };
  return { samples, chromSamples, cand };
}

// ============================================================================
// 1. Source-level definitions
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_breedingCardKaryotypePerSample defined',  /function _breedingCardKaryotypePerSample\(/.test(html));
ok('_wilcoxonRankSumP defined',                 /function _wilcoxonRankSumP\(/.test(html));
ok('_summarizeFROHGroup defined',               /function _summarizeFROHGroup\(/.test(html));
ok('_perArrangementFROH defined',               /function _perArrangementFROH\(/.test(html));
ok('_arrangementMAF defined',                   /function _arrangementMAF\(/.test(html));
ok('_carrierByK8Table defined',                 /function _carrierByK8Table\(/.test(html));
ok('_generatePairingAdvice defined',            /function _generatePairingAdvice\(/.test(html));
ok('_buildBreedingCard defined',                /function _buildBreedingCard\(/.test(html));

ok('_BREEDING_FROH_P_THRESHOLD defined',         /const _BREEDING_FROH_P_THRESHOLD\s*=\s*0\.05/.test(html));
ok('_BREEDING_MAF_IMBALANCE_THRESH defined',     /const _BREEDING_MAF_IMBALANCE_THRESH\s*=\s*0\.10/.test(html));
ok('_BREEDING_RECOMBINANT_THRESH defined',       /const _BREEDING_RECOMBINANT_THRESH\s*=\s*0\.20/.test(html));

// ============================================================================
// 2. Window exports
// ============================================================================
console.log('\n=== 2. Window exports ===');

ok('window._breedingCardKaryotypePerSample',
   /window\._breedingCardKaryotypePerSample\s*=\s*_breedingCardKaryotypePerSample/.test(html));
ok('window._wilcoxonRankSumP',
   /window\._wilcoxonRankSumP\s*=\s*_wilcoxonRankSumP/.test(html));
ok('window._summarizeFROHGroup',
   /window\._summarizeFROHGroup\s*=\s*_summarizeFROHGroup/.test(html));
ok('window._perArrangementFROH',
   /window\._perArrangementFROH\s*=\s*_perArrangementFROH/.test(html));
ok('window._arrangementMAF',
   /window\._arrangementMAF\s*=\s*_arrangementMAF/.test(html));
ok('window._carrierByK8Table',
   /window\._carrierByK8Table\s*=\s*_carrierByK8Table/.test(html));
ok('window._generatePairingAdvice',
   /window\._generatePairingAdvice\s*=\s*_generatePairingAdvice/.test(html));
ok('window._buildBreedingCard',
   /window\._buildBreedingCard\s*=\s*_buildBreedingCard/.test(html));

// ============================================================================
// 3. Karyotype derivation behavioral
// ============================================================================
console.log('\n=== 3. Karyotype derivation ===');

const ctx = buildSandbox();
const W = ctx.window;

// 3.1 K=3 clean
{
  const { cand } = k3Fixture();
  const k = W._breedingCardKaryotypePerSample(cand);
  ok('K=3 returns object',                       k !== null);
  ok('K=3 karyotype length',                     k.karyotype.length === cand.locked_labels.length);
  ok('K=3 first 3 are HOM_REF',
     k.karyotype.slice(0, 3).every(v => v === 'HOM_REF'));
  ok('K=3 next 5 are HET',
     k.karyotype.slice(3, 8).every(v => v === 'HET'));
  ok('K=3 last 3 are HOM_INV',
     k.karyotype.slice(8).every(v => v === 'HOM_INV'));
  ok('K=3 counts.HOM_REF',                       k.counts.HOM_REF === 3);
  ok('K=3 counts.HET',                            k.counts.HET === 5);
  ok('K=3 counts.HOM_INV',                       k.counts.HOM_INV === 3);
  ok('K=3 counts.n_classified',                  k.counts.n_classified === 11);
  ok('K=3 hom_band_idx_by_role.REF',             k.hom_band_idx_by_role.REF === 0);
  ok('K=3 hom_band_idx_by_role.INV',             k.hom_band_idx_by_role.INV === 2);
  ok('K=3 not multi_haplotype',                   k.multi_haplotype === false);
  ok('K=3 not ambiguous_role',                    k.ambiguous_role === false);
  ok('K=3 classification_summary',                k.classification_summary === 'K=3 · 2H+1T');
  ok('K=3 consistency',                           k.consistency === 'clean');
}

// 3.2 K=4 multi-haplotype (3 HOM + 1 HET)
{
  const cand = {
    K: 4,
    locked_labels: new Int8Array([0,0, 1,1, 2,2, 3,3]),
    h_classification: {
      K: 4,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -3.0 },
        { band_idx: 1, classification: 'HOM', median_pc1:  0.0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  3.0 },
        { band_idx: 3, classification: 'HET', median_pc1:  1.5 },
      ],
      band_counts: { n_hom: 3, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const k = W._breedingCardKaryotypePerSample(cand);
  ok('K=4 multi: REF = leftmost-PC1 band 0',     k.hom_band_idx_by_role.REF === 0);
  ok('K=4 multi: INV = rightmost-PC1 band 2',    k.hom_band_idx_by_role.INV === 2);
  ok('K=4 multi: MID contains band 1',           k.hom_band_idx_by_role.MID.length === 1 &&
                                                  k.hom_band_idx_by_role.MID[0] === 1);
  ok('K=4 multi: counts.HOM_MID',                k.counts.HOM_MID === 2);
  ok('K=4 multi: multi_haplotype = true',         k.multi_haplotype === true);
}

// 3.3 Single-HOM ambiguous role
{
  const cand = {
    K: 2,
    locked_labels: new Int8Array([0, 0, 1]),
    h_classification: {
      K: 2,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -1.0 },
        { band_idx: 1, classification: 'HET', median_pc1:  1.0 },
      ],
      band_counts: { n_hom: 1, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const k = W._breedingCardKaryotypePerSample(cand);
  ok('single-HOM: ambiguous_role = true',         k.ambiguous_role === true);
  ok('single-HOM: REF/INV both null',             k.hom_band_idx_by_role.REF === null &&
                                                  k.hom_band_idx_by_role.INV === null);
  ok('single-HOM: counts.HOM_AMBIGUOUS_ROLE = 2', k.counts.HOM_AMBIGUOUS_ROLE === 2);
  ok('single-HOM: counts.HOM_REF = 0',            k.counts.HOM_REF === 0);
  ok('single-HOM: karyotype values are HOM',
     k.karyotype.slice(0, 2).every(v => v === 'HOM'));
}

// 3.4 AMBIGUOUS / NO_DOSAGE bands → null karyotype
{
  const cand = {
    K: 4,
    locked_labels: new Int8Array([0, 1, 2, 3]),
    h_classification: {
      K: 4,
      bands: [
        { band_idx: 0, classification: 'HOM',         median_pc1: -1.0 },
        { band_idx: 1, classification: 'HET',         median_pc1:  0.0 },
        { band_idx: 2, classification: 'AMBIGUOUS',   median_pc1:  0.5 },
        { band_idx: 3, classification: 'NO_DOSAGE',   median_pc1:  1.0 },
      ],
      band_counts: { n_hom: 1, n_het: 1, n_ambiguous: 1, n_no_dosage: 1 },
      implied_regime: { consistency: 'unknown' },
    },
  };
  const k = W._breedingCardKaryotypePerSample(cand);
  ok('mixed: si=2 AMBIGUOUS → null',              k.karyotype[2] === null);
  ok('mixed: si=3 NO_DOSAGE → null',              k.karyotype[3] === null);
  ok('mixed: counts.AMBIGUOUS = 1',               k.counts.AMBIGUOUS === 1);
  ok('mixed: counts.NO_DOSAGE = 1',               k.counts.NO_DOSAGE === 1);
  ok('mixed: counts.n_unclassified = 2',          k.counts.n_unclassified === 2);
}

// 3.5 Defensive returns
ok('null cand → null',                            W._breedingCardKaryotypePerSample(null) === null);
ok('cand without locked_labels → null',           W._breedingCardKaryotypePerSample({}) === null);
ok('cand without h_classification → null',
   W._breedingCardKaryotypePerSample({ locked_labels: new Int8Array([0,1]) }) === null);

// ============================================================================
// 4. Wilcoxon rank-sum behavioral
// ============================================================================
console.log('\n=== 4. Wilcoxon rank-sum ===');

// 4.1 Perfect separation 5v5 — hand-verified expected values
{
  const r = W._wilcoxonRankSumP([1,2,3,4,5], [6,7,8,9,10]);
  ok('5v5 perfect: n_a=5',                        r.n_a === 5);
  ok('5v5 perfect: n_b=5',                        r.n_b === 5);
  ok('5v5 perfect: R_a=15',                       r.R_a === 15);
  ok('5v5 perfect: U_a=0',                        r.U_a === 0);
  ok('5v5 perfect: mu=12.5',                      r.mu === 12.5);
  ok('5v5 perfect: sigma2 ≈ 22.917',              approx(r.sigma2, 22.9166666, 1e-5));
  ok('5v5 perfect: z ≈ 2.5067',                   approx(r.z, 2.5067, 1e-3));
  ok('5v5 perfect: p ≈ 0.0122',                   approx(r.p_two_sided, 0.0122, 1e-4));
  ok('5v5 perfect: direction=a_lower',            r.direction === 'a_lower');
  ok('5v5 perfect: n_tie_groups=0',               r.n_tie_groups === 0);
  ok('5v5 perfect: tie_correction_factor=0',      r.tie_correction_factor === 0);
}

// 4.2 Symmetry — swap groups → same p, opposite direction
{
  const r1 = W._wilcoxonRankSumP([1,2,3,4,5], [6,7,8,9,10]);
  const r2 = W._wilcoxonRankSumP([6,7,8,9,10], [1,2,3,4,5]);
  ok('symmetry: p(a,b) === p(b,a)',
     Math.abs(r1.p_two_sided - r2.p_two_sided) < 1e-12);
  ok('symmetry: direction flips',
     r1.direction === 'a_lower' && r2.direction === 'a_higher');
  ok('symmetry: U_a + U_b = n_a*n_b',
     r1.U_a + r2.U_a === r1.n_a * r1.n_b);
}

// 4.3 Full ties — degenerate case
{
  const r = W._wilcoxonRankSumP([5,5,5], [5,5,5]);
  ok('full ties: U_a == mu',                      r.U_a === r.mu);
  ok('full ties: sigma=0',                        r.sigma === 0);
  ok('full ties: z=NaN',                          Number.isNaN(r.z));
  ok('full ties: p=NaN',                          Number.isNaN(r.p_two_sided));
  ok('full ties: direction=equal',                r.direction === 'equal');
  ok('full ties: n_tie_groups=1',                 r.n_tie_groups === 1);
}

// 4.4 Mixed ties — tie correction applied
{
  const r = W._wilcoxonRankSumP([1,1,2,3], [1,2,3,4]);
  ok('mixed ties: n_tie_groups=3',                r.n_tie_groups === 3);
  ok('mixed ties: tie_correction > 0',            r.tie_correction_factor > 0);
  // Verify: with a values [1,1,2,3], R_a is sum of ranks of those values
  // in the pooled sorted [1,1,1,2,2,3,3,4]. The three 1's get rank 2
  // each, the two 2's get rank 4.5 each, the two 3's get rank 6.5 each,
  // 4 gets rank 8. Group a values (1,1,2,3) get ranks (2, 2, 4.5, 6.5)
  // = sum 15. So R_a should be 15.
  ok('mixed ties: R_a = 15 (hand-verified)',      r.R_a === 15);
  ok('mixed ties: U_a = 5',                       r.U_a === 5);
}

// 4.5 Realistic large-n with shift
{
  const a = [], b = [];
  for (let i = 0; i < 50; i++) {
    a.push(0.20 + Math.sin(i * 0.7) * 0.03);     // ~0.20
    b.push(0.30 + Math.sin(i * 0.3) * 0.03);     // ~0.30 (shift +0.10)
  }
  const r = W._wilcoxonRankSumP(a, b);
  ok('50v50 shift+0.10: significant',             r.p_two_sided < 1e-10);
  ok('50v50 shift+0.10: direction=a_lower',       r.direction === 'a_lower');
}

// 4.6 Edge cases
ok('null,null → null',                            W._wilcoxonRankSumP(null, null) === null);
ok('non-array → null',                            W._wilcoxonRankSumP(5, 6) === null);
ok('empty arrays → null',                         W._wilcoxonRankSumP([], []) === null);
ok('one empty → null',                            W._wilcoxonRankSumP([1,2], []) === null);
ok('all-NaN → null',                              W._wilcoxonRankSumP([NaN, NaN], [NaN, NaN]) === null);
ok('mix NaN+finite drops NaN',
   W._wilcoxonRankSumP([1, NaN, 2], [3, 4, NaN]).n_a === 2);

// 4.7 1v1 case — sigma should be 0 (only one ordering possible after avg)
{
  const r = W._wilcoxonRankSumP([1], [2]);
  // n=1+1, no ties: mu = 0.5, sigma2 = 1*1*(2+1)/12 = 0.25, sigma=0.5
  // U_a = ranks(1) - 1*2/2 = 1 - 1 = 0 → diff = -0.5 → |diff|=0.5 →
  // continuity-corrected returns p=1.
  ok('1v1: p = 1 (continuity correction)',        r.p_two_sided === 1);
}

// ============================================================================
// 5. _summarizeFROHGroup
// ============================================================================
console.log('\n=== 5. _summarizeFROHGroup ===');

{
  const s = W._summarizeFROHGroup([0.20, 0.21, 0.19]);
  ok('summarize n=3: n',                           s.n === 3);
  ok('summarize n=3: mean ≈ 0.2',                  approx(s.mean, 0.20, 1e-9));
  ok('summarize n=3: median = 0.20',               s.median === 0.20);
  ok('summarize n=3: sd > 0',                      s.sd > 0);
  ok('summarize n=3: min = 0.19',                  s.min === 0.19);
  ok('summarize n=3: max = 0.21',                  s.max === 0.21);
  ok('summarize n=3: q1 ≈ 0.195',                  approx(s.q1, 0.195, 1e-9));
  ok('summarize n=3: q3 ≈ 0.205',                  approx(s.q3, 0.205, 1e-9));
}
{
  const s = W._summarizeFROHGroup([]);
  ok('summarize empty: n=0',                       s.n === 0);
  ok('summarize empty: mean=NaN',                  Number.isNaN(s.mean));
}
{
  const s = W._summarizeFROHGroup([0.5]);
  ok('summarize n=1: n=1',                         s.n === 1);
  ok('summarize n=1: sd=0',                        s.sd === 0);
  ok('summarize n=1: mean=median=0.5',             s.mean === 0.5 && s.median === 0.5);
}
{
  const s = W._summarizeFROHGroup([0.1, NaN, 0.2, undefined, 0.3]);
  ok('summarize drops non-finite: n=3',            s.n === 3);
  ok('summarize after drop: mean ≈ 0.2',           approx(s.mean, 0.2, 1e-9));
}

// ============================================================================
// 6. _perArrangementFROH
// ============================================================================
console.log('\n=== 6. _perArrangementFROH ===');

const burdenCtx = buildSandbox();
{
  const fx = k3Fixture();
  burdenCtx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  burdenCtx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples,
  });
  const r = burdenCtx.window._perArrangementFROH(fx.cand);
  ok('available=true',                              r.available === true);
  ok('reason=null',                                  r.reason === null);
  ok('groups.HOM_REF.n=3',                           r.groups.HOM_REF.n === 3);
  ok('groups.HET.n=5',                                r.groups.HET.n === 5);
  ok('groups.HOM_INV.n=3',                           r.groups.HOM_INV.n === 3);
  ok('all_carriers.n=8 (HET+INV)',                   r.all_carriers.n === 8);
  ok('groups.HOM_REF.mean=0.20',                     approx(r.groups.HOM_REF.mean, 0.20, 1e-9));
  ok('groups.HOM_INV.mean=0.32',                     approx(r.groups.HOM_INV.mean, 0.32, 1e-9));
  ok('delta_mean ≈ 0.12',                            approx(r.delta_mean_inv_minus_ref, 0.12, 1e-9));
  ok('wilcoxon_ref_vs_inv exists',                   r.wilcoxon_ref_vs_inv !== null);
  ok('wilcoxon_ref_vs_inv direction=a_lower',        r.wilcoxon_ref_vs_inv.direction === 'a_lower');
  ok('wilcoxon_ref_vs_carriers exists',              r.wilcoxon_ref_vs_carriers !== null);
  ok('n_resolved=11',                                 r.n_resolved === 11);
  ok('n_unresolved=0',                                r.n_unresolved === 0);
  ok('values not included by default',
     !('values' in r.groups.HOM_REF));

  const r2 = burdenCtx.window._perArrangementFROH(fx.cand, { include_values: true });
  ok('include_values=true → values present',
     Array.isArray(r2.groups.HOM_REF.values) && r2.groups.HOM_REF.values.length === 3);
}

// 6.1 Stub when no cohort
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: [{ cga: 'CGA001' }] };
  const fx = k3Fixture();
  const r = c.window._perArrangementFROH(fx.cand);
  ok('no cohort: available=false',                   r.available === false);
  ok('no cohort: reason=no_cohort_diversity',        r.reason === 'no_cohort_diversity');
  // Karyotype summary still present for the renderer
  ok('no cohort: karyotype_summary set',             r.karyotype_summary !== null);
}

// 6.2 Stub when no h_classification
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: [{ cga: 'CGA001' }] };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [{ sample_id: 'CGA001', f_roh: 0.2 }],
  });
  const r = c.window._perArrangementFROH({ K: 3, locked_labels: new Int8Array([0]) });
  ok('no h_class: available=false',                  r.available === false);
  ok('no h_class: reason=no_h_classification',       r.reason === 'no_h_classification');
}

// 6.3 No candidate
{
  const c = buildSandbox();
  const r = c.window._perArrangementFROH(null);
  ok('null cand: available=false',                   r.available === false);
  ok('null cand: reason=no_candidate',               r.reason === 'no_candidate');
}

// 6.4 Samples without CGA → counted in n_unresolved
{
  const c = buildSandbox();
  c.state.data = {
    chrom: 'LG28',
    samples: [{ cga: 'CGA001' }, { cga: 'CGA002' }, { ind: 'Ind3' }],   // 3rd has no cga
  };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', f_roh: 0.20, k8: 'K1' },
      { sample_id: 'CGA002', f_roh: 0.25, k8: 'K1' },
    ],
  });
  const cand = {
    K: 3,
    locked_labels: new Int8Array([0, 1, 2]),
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2 },
        { band_idx: 1, classification: 'HET', median_pc1:  0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  2 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const r = c.window._perArrangementFROH(cand);
  ok('missing CGA on sample → n_unresolved bumps', r.n_unresolved === 1);
  ok('missing CGA: n_resolved counts the resolved 2', r.n_resolved === 2);
}

// ============================================================================
// 7. _arrangementMAF
// ============================================================================
console.log('\n=== 7. _arrangementMAF ===');

ok('balanced (10,10,10): maf=0.5',
   approx(W._arrangementMAF({HOM_REF:10, HET:10, HOM_INV:10}).maf, 0.5));
ok('balanced: minor_allele=tie',
   W._arrangementMAF({HOM_REF:10, HET:10, HOM_INV:10}).minor_allele === 'tie');
{
  const m = W._arrangementMAF({HOM_REF:50, HET:6, HOM_INV:0});
  ok('skewed (50,6,0): maf ≈ 0.054',                approx(m.maf, 6/112, 1e-9));
  ok('skewed: minor_allele=INV',                     m.minor_allele === 'INV');
  ok('skewed: p_ref ≈ 0.946',                        approx(m.p_ref, 106/112, 1e-9));
}
{
  const m = W._arrangementMAF({HOM_REF:0, HET:6, HOM_INV:50});
  ok('skewed reverse: minor_allele=REF',             m.minor_allele === 'REF');
}
ok('null counts → null',                              W._arrangementMAF(null) === null);
ok('all-zero counts → null',                         W._arrangementMAF({HOM_REF:0, HET:0, HOM_INV:0}) === null);

// ============================================================================
// 8. _carrierByK8Table
// ============================================================================
console.log('\n=== 8. _carrierByK8Table ===');

const k8Ctx = buildSandbox();
{
  const fx = k3Fixture();
  k8Ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  k8Ctx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples,
  });
  const t = k8Ctx.window._carrierByK8Table(fx.cand);
  ok('K8 table available',                          t.available === true);
  ok('K8 clusters K1,K2,K3 sorted',                 JSON.stringify(t.k8_clusters) === '["K1","K2","K3"]');
  ok('K8 columns: HOM_REF/HET/HOM_INV',
     JSON.stringify(t.karyotypes) === '["HOM_REF","HET","HOM_INV"]');
  ok('K8 totals match cohort sizes',
     t.totals.HOM_REF === 3 && t.totals.HET === 5 && t.totals.HOM_INV === 3);
  ok('K8 row count = 3',                             t.rows.length === 3);
  // For each row, frac sums to 1 (within float tolerance) when all classified
  const r0 = t.rows[0];
  const fracSum = r0.frac.HOM_REF + r0.frac.HET + r0.frac.HOM_INV;
  ok('K8 row frac sums to 1',                        approx(fracSum, 1, 1e-9));
}

// 8.1 K8 column adapts to multi_haplotype
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: Array.from({length: 8}, (_, i) => ({ cga: 'CGA' + String(i+1).padStart(3,'0') })) };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1, samples:
      Array.from({length: 8}, (_, i) => ({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K1', f_roh: 0.25 })),
  });
  const cand = {
    K: 4,
    locked_labels: new Int8Array([0,0, 1,1, 2,2, 3,3]),
    h_classification: {
      K: 4,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -3 },
        { band_idx: 1, classification: 'HOM', median_pc1:  0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  3 },
        { band_idx: 3, classification: 'HET', median_pc1:  1.5 },
      ],
      band_counts: { n_hom: 3, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const t = c.window._carrierByK8Table(cand);
  ok('multi-hap: HOM_MID column present',           t.karyotypes.includes('HOM_MID'));
}

// 8.2 K8 column adapts to ambiguous_role
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }] };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', k8: 'K1', f_roh: 0.2 },
      { sample_id: 'CGA002', k8: 'K1', f_roh: 0.2 },
      { sample_id: 'CGA003', k8: 'K2', f_roh: 0.3 },
    ],
  });
  const cand = {
    K: 2,
    locked_labels: new Int8Array([0, 0, 1]),
    h_classification: {
      K: 2,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -1 },
        { band_idx: 1, classification: 'HET', median_pc1:  1 },
      ],
      band_counts: { n_hom: 1, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const t = c.window._carrierByK8Table(cand);
  ok('ambig-role: HOM column present',              t.karyotypes.includes('HOM'));
}

// 8.3 Stubs
ok('K8 stub: no cohort',                             buildSandbox().window._carrierByK8Table(k3Fixture().cand).available === false);

// ============================================================================
// 9. _generatePairingAdvice — rule firing
// ============================================================================
console.log('\n=== 9. Pairing advice rule firing ===');

// 9.1 F_ROH asymmetric — synthetic large-n with clear separation
{
  const c = buildSandbox();
  const n = 90;
  c.state.data = { chrom: 'LG28', samples: Array.from({length: n}, (_, i) => ({ cga: 'CGA' + String(i+1).padStart(3,'0') })) };
  const samps = [];
  for (let i = 0; i < 30; i++) samps.push({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K1', f_roh: 0.20 + Math.sin(i*0.7)*0.02 });
  for (let i = 30; i < 60; i++) samps.push({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K2', f_roh: 0.25 + Math.sin(i*0.7)*0.02 });
  for (let i = 60; i < 90; i++) samps.push({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K3', f_roh: 0.32 + Math.sin(i*0.7)*0.02 });
  c.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: samps });
  const labels = new Int8Array(n);
  for (let i = 0; i < 30; i++) labels[i] = 0;
  for (let i = 30; i < 60; i++) labels[i] = 1;
  for (let i = 60; i < 90; i++) labels[i] = 2;
  const cand = {
    id: 'big', chrom: 'LG28', start_bp: 14e6, end_bp: 16.5e6, K: 3,
    locked_labels: labels,
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2 },
        { band_idx: 1, classification: 'HET', median_pc1:  0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  2 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'clean' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  const fr = card.advice.find(a => a.kind === 'froh_asymmetric');
  ok('FROH asym fires when INV>>REF',                fr !== undefined);
  ok('FROH asym severity=strong',                    fr && fr.severity === 'strong');
  ok('FROH asym evidence has p_value',               fr && fr.evidence && Number.isFinite(fr.evidence.p_value));
  ok('FROH asym evidence has delta',                 fr && Number.isFinite(fr.evidence.delta_mean_inv_minus_ref));
  ok('FROH asym evidence delta > 0',                 fr && fr.evidence.delta_mean_inv_minus_ref > 0);
  // No default advice when actionable rule fired
  const def = card.advice.find(a => a.kind === 'default');
  ok('no default when actionable rule fires',        def === undefined);
}

// 9.2 F_ROH inverted asymmetry — direction matters
{
  const c = buildSandbox();
  const n = 60;
  c.state.data = { chrom: 'LG28', samples: Array.from({length: n}, (_, i) => ({ cga: 'CGA' + String(i+1).padStart(3,'0') })) };
  const samps = [];
  // Reverse: REF has HIGH F_ROH, INV has LOW
  for (let i = 0; i < 30; i++) samps.push({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K1', f_roh: 0.32 + Math.sin(i*0.7)*0.02 });
  for (let i = 30; i < 60; i++) samps.push({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K2', f_roh: 0.20 + Math.sin(i*0.7)*0.02 });
  c.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: samps });
  const labels = new Int8Array(n);
  for (let i = 0; i < 30; i++) labels[i] = 0;
  for (let i = 30; i < 60; i++) labels[i] = 2;       // INV in band 2
  const cand = {
    id: 'inv', chrom: 'LG28', start_bp: 14e6, end_bp: 16.5e6, K: 3,
    locked_labels: labels,
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2 },
        { band_idx: 1, classification: 'HET', median_pc1:  0 },    // empty in this synthetic
        { band_idx: 2, classification: 'HOM', median_pc1:  2 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'clean' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  const inv = card.advice.find(a => a.kind === 'froh_inverted_asymmetry');
  ok('FROH inverted asym fires when REF>>INV',       inv !== undefined);
  ok('FROH inverted asym severity=info',             inv && inv.severity === 'info');
}

// 9.3 MAF imbalance fires
{
  const c = buildSandbox();
  const n = 56;
  c.state.data = { chrom: 'LG28', samples: Array.from({length: n}, (_, i) => ({ cga: 'CGA' + String(i+1).padStart(3,'0') })) };
  const samps = Array.from({length: n}, (_, i) => ({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K1', f_roh: 0.25 }));
  c.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: samps });
  const labels = new Int8Array(n);
  for (let i = 0; i < 50; i++) labels[i] = 0;
  for (let i = 50; i < 56; i++) labels[i] = 1;   // 6 HET
  // bands have 2 HOM positions; band 2 is empty
  const cand = {
    id: 'imb', chrom: 'LG28', start_bp: 14e6, end_bp: 16.5e6, K: 3,
    locked_labels: labels,
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2 },
        { band_idx: 1, classification: 'HET', median_pc1:  0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  2 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'clean' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  const imb = card.advice.find(a => a.kind === 'maf_imbalanced');
  ok('MAF imbalance fires at maf<0.10',              imb !== undefined);
  ok('MAF imbalance severity=warn',                  imb && imb.severity === 'warn');
  ok('MAF imbalance evidence has minor_allele',      imb && imb.evidence.minor_allele === 'INV');
}

// 9.4 Multi-haplotype caveat fires
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: Array.from({length: 8}, (_, i) => ({ cga: 'CGA' + String(i+1).padStart(3,'0') })) };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: Array.from({length: 8}, (_, i) => ({ sample_id: 'CGA' + String(i+1).padStart(3,'0'), k8: 'K1', f_roh: 0.25 })),
  });
  const cand = {
    K: 4,
    locked_labels: new Int8Array([0,0, 1,1, 2,2, 3,3]),
    h_classification: {
      K: 4,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -3 },
        { band_idx: 1, classification: 'HOM', median_pc1:  0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  3 },
        { band_idx: 3, classification: 'HET', median_pc1:  1.5 },
      ],
      band_counts: { n_hom: 3, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  const mh = card.advice.find(a => a.kind === 'multi_haplotype_caveat');
  ok('multi_haplotype caveat fires',                 mh !== undefined);
  ok('multi_haplotype severity=info',                mh && mh.severity === 'info');
}

// 9.5 Ambiguous-role caveat fires
{
  const c = buildSandbox();
  c.state.data = { chrom: 'LG28', samples: [{ cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' }] };
  c.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [{sample_id:'CGA001',k8:'K1',f_roh:0.2},{sample_id:'CGA002',k8:'K1',f_roh:0.2},{sample_id:'CGA003',k8:'K1',f_roh:0.3}],
  });
  const cand = {
    K: 2,
    locked_labels: new Int8Array([0, 0, 1]),
    h_classification: {
      K: 2,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -1 },
        { band_idx: 1, classification: 'HET', median_pc1:  1 },
      ],
      band_counts: { n_hom: 1, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'partial_obs_minor' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  const amb = card.advice.find(a => a.kind === 'ambiguous_role_caveat');
  ok('ambiguous_role caveat fires',                  amb !== undefined);
}

// 9.6 data_pending always emitted
{
  const c = buildSandbox();
  const fx = k3Fixture();
  c.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  c.window._storeCohortDiversity({tool:'cohort_diversity_v1',schema_version:1,samples:fx.samples});
  const card = c.window._buildBreedingCard(fx.cand);
  const pending = card.advice.filter(a => a.kind === 'data_pending');
  ok('two data_pending items emitted',               pending.length === 2);
  ok('damaging_load placeholder',
     pending.some(a => a.evidence && a.evidence.needs === 'damaging_load_v1'));
  ok('recombinant placeholder',
     pending.some(a => a.evidence && a.evidence.needs === 'recombinant_dosage_changepoint_detector'));
}

// 9.7 Default rule fires only when no actionable rules
{
  // Use the original k3Fixture which has insufficient n for significance
  const c = buildSandbox();
  const fx = k3Fixture();
  c.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  c.window._storeCohortDiversity({tool:'cohort_diversity_v1',schema_version:1,samples:fx.samples});
  const card = c.window._buildBreedingCard(fx.cand);
  // Wilcoxon p≈0.0518 > 0.05 → no FROH rule. MAF balanced → no MAF rule.
  // Default should fire.
  const def = card.advice.find(a => a.kind === 'default');
  ok('default rule fires when no actionable rules', def !== undefined);
}

// ============================================================================
// 10. _buildBreedingCard end-to-end
// ============================================================================
console.log('\n=== 10. _buildBreedingCard end-to-end ===');

const cardCtx = buildSandbox();
const fxCard = k3Fixture();
cardCtx.state.data = { chrom: 'LG28', samples: fxCard.chromSamples };
cardCtx.window._storeCohortDiversity({tool:'cohort_diversity_v1',schema_version:1,samples:fxCard.samples});
const card = cardCtx.window._buildBreedingCard(fxCard.cand);

ok('schema = breeding_readiness_card_v1',           card.schema === 'breeding_readiness_card_v1');
ok('generated_at is ISO string',                     /^\d{4}-\d{2}-\d{2}T/.test(card.generated_at));
ok('candidate.id preserved',                         card.candidate.id === 'LG28_test');
ok('candidate.span_bp computed',                    card.candidate.span_bp === 2500000);
ok('candidate.span_mb = 2.5',                       card.candidate.span_mb === 2.5);
ok('candidate.tier preserved',                       card.candidate.tier === 'Tier 1');
ok('karyotype.summary present',                      card.karyotype.summary.HOM_REF === 3);
ok('karyotype.classification_summary',               card.karyotype.classification_summary === 'K=3 · 2H+1T');
ok('burden.available=true',                          card.burden.available === true);
ok('ancestry.available=true',                        card.ancestry.available === true);
ok('maf computed',                                   card.maf !== null && card.maf.maf > 0);
ok('coverage.n_total=11',                            card.coverage.n_total === 11);
ok('coverage.n_resolved=11',                         card.coverage.n_resolved === 11);
ok('advice is non-empty array',                      Array.isArray(card.advice) && card.advice.length > 0);
ok('atlas_url null in vm sandbox',                   card.atlas_url === null);

// Card is JSON-serializable
let jsonable = false;
try { JSON.stringify(card); jsonable = true; } catch (_) {}
ok('card is JSON-serializable',                      jsonable);

// 10.1 No candidate
ok('null cand → null',                              cardCtx.window._buildBreedingCard(null) === null);

// 10.2 Card with no h_classification (no dosage data)
{
  const candNoHC = {
    id: 'no_hc', chrom: 'LG28', start_bp: 1e6, end_bp: 2e6, K: 3,
    locked_labels: new Int8Array([0, 1, 2]),
  };
  const c = cardCtx.window._buildBreedingCard(candNoHC);
  ok('no h_class: schema present',                   c.schema === 'breeding_readiness_card_v1');
  ok('no h_class: candidate fields populated',       c.candidate.id === 'no_hc');
  ok('no h_class: karyotype.available=false',        c.karyotype.available === false);
  ok('no h_class: burden.available=false',           c.burden.available === false);
  ok('no h_class: ancestry.available=false',         c.ancestry.available === false);
  ok('no h_class: still JSON-serializable',
     (() => { try { JSON.stringify(c); return true; } catch (_) { return false; } })());
}

// 10.3 No cohort_diversity
{
  const c2 = buildSandbox();
  c2.state.data = { chrom: 'LG28', samples: fxCard.chromSamples };
  // Don't load cohort_diversity
  const c = c2.window._buildBreedingCard(fxCard.cand);
  ok('no cohort: candidate populated',               c.candidate.id === 'LG28_test');
  ok('no cohort: karyotype populated',               c.karyotype.summary !== undefined &&
                                                      c.karyotype.summary.HOM_REF === 3);
  ok('no cohort: burden.available=false',            c.burden.available === false);
  ok('no cohort: burden.reason=no_cohort_diversity', c.burden.reason === 'no_cohort_diversity');
  ok('no cohort: ancestry.available=false',          c.ancestry.available === false);
  // Advice still emits caveats + data_pending + default
  ok('no cohort: advice still emitted',              c.advice.length >= 3);
}

// ============================================================================
// 11. Full pipeline against the real fixture (json/cohort_diversity_v1.json)
// ============================================================================
console.log('\n=== 11. Real fixture pipeline ===');

const FIXTURE_PATH = path.resolve(__dirname, '..', 'json', 'cohort_diversity_v1.json');
let fixture = null;
try { fixture = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf8')); } catch (_) {}

if (fixture) {
  const c = buildSandbox();
  // Synthetic 226-sample chrom matching the fixture's CGA ids.
  const chromSamples = fixture.samples.map(s => ({ cga: s.sample_id }));
  c.state.data = { chrom: 'LG28', samples: chromSamples };
  c.window._storeCohortDiversity(fixture);

  // Build a synthetic candidate that splits the 226 cohort 60/106/60 — the
  // actual ~2.89 Mb LG28 inversion karyotype reported in the user memories
  // (60 REF, 106 HET, 60 INV). We can't recover the true assignment here
  // since we don't have the chromosome's real K-means labels — just
  // verify the pipeline runs end-to-end on real-sized n.
  const labels = new Int8Array(226);
  for (let i = 0; i < 60; i++) labels[i] = 0;
  for (let i = 60; i < 166; i++) labels[i] = 1;
  for (let i = 166; i < 226; i++) labels[i] = 2;

  const cand = {
    id: 'LG28_real_n_test', chrom: 'LG28',
    start_bp: 13800000, end_bp: 16700000, K: 3, tier: 'Tier 1',
    locked_labels: labels,
    h_classification: {
      K: 3,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -2.0 },
        { band_idx: 1, classification: 'HET', median_pc1:  0.0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  2.0 },
      ],
      band_counts: { n_hom: 2, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'clean' },
    },
  };
  const card = c.window._buildBreedingCard(cand);
  ok('real-n: schema present',                       card.schema === 'breeding_readiness_card_v1');
  ok('real-n: 226 samples in coverage',              card.coverage.n_total === 226);
  ok('real-n: all 226 resolved',                     card.coverage.n_resolved === 226);
  ok('real-n: karyotype counts to 226',              card.karyotype.summary.HOM_REF +
                                                       card.karyotype.summary.HET +
                                                       card.karyotype.summary.HOM_INV === 226);
  ok('real-n: burden has wilcoxon ref vs inv',
     card.burden.wilcoxon_ref_vs_inv && Number.isFinite(card.burden.wilcoxon_ref_vs_inv.p_two_sided));
  ok('real-n: ancestry rows match # of K=8 clusters',
     card.ancestry.rows.length >= 1 && card.ancestry.rows.length <= 8);
  ok('real-n: card JSON-serializable',
     (() => { try { JSON.stringify(card); return true; } catch (_) { return false; } })());

  // Print the actual computed Wilcoxon for inspection (not assertion):
  console.log('  [info] real-n Wilcoxon REF vs INV p =',
              card.burden.wilcoxon_ref_vs_inv.p_two_sided);
  console.log('  [info] real-n delta_mean_inv_minus_ref =',
              card.burden.delta_mean_inv_minus_ref);
}

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
