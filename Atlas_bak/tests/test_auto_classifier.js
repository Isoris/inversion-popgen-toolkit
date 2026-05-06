// Tests for turn 2f — auto-classifier + vocabulary system.
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
  url: 'http://localhost/',
});
const w = dom.window;

function run() {
  const failures = [];
  let testNum = 0;
  function t(name, fn) {
    testNum++;
    try { fn(); console.log(`  PASS [${testNum}] ${name}`); }
    catch (e) { failures.push({ name, err: e.message }); console.log(`  FAIL [${testNum}] ${name}: ${e.message}`); }
  }
  function eq(a, b, m) { if (a !== b) throw new Error(`${m||''} expected ${JSON.stringify(b)}, got ${JSON.stringify(a)}`); }

  if (typeof w.autoClassifyCandidate !== 'function'
      || typeof w.applyAutoClassificationToCandidate !== 'function'
      || typeof w.getHaplotypeVocab !== 'function'
      || typeof w.setHaplotypeVocab !== 'function'
      || typeof w.autoSelectVocabForCandidate !== 'function') {
    console.log('  FAIL: turn 2f functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2f functions exposed');

  function clearLS() { try { w.localStorage.clear(); } catch (_) {} }

  // Build a synthetic state.data with controllable per-window PC1 distributions
  function makeState(K, perBandSpec, n_windows) {
    // perBandSpec[k] = { n_samples, centroid, drift }
    // - n_samples: how many fish in band k
    // - centroid: PC1 mean for fish in band k
    // - drift: per-fish per-window standard deviation (het proxy)
    n_windows = n_windows || 30;
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.candidates = {};
    _state.candidates_detailed = {};
    _state.activeMode = 'default';
    _state.haplotypeVocabs = {};

    // Build labels and synthetic PC1 data
    const labels = [];
    const sampleBand = [];
    for (let k = 0; k < K; k++) {
      for (let i = 0; i < perBandSpec[k].n_samples; i++) {
        labels.push(k);
        sampleBand.push(k);
      }
    }
    const n_samples = labels.length;
    // Build per-window PC1 with noise
    const windows = [];
    let rngState = 12345;
    function rand() {
      // Deterministic LCG for reproducibility
      rngState = (rngState * 1103515245 + 12345) & 0x7fffffff;
      return rngState / 0x7fffffff;
    }
    for (let w_idx = 0; w_idx < n_windows; w_idx++) {
      const pc1 = new Float32Array(n_samples);
      for (let s = 0; s < n_samples; s++) {
        const k = sampleBand[s];
        const c = perBandSpec[k].centroid;
        const drift = perBandSpec[k].drift;
        // Box-Muller for gaussian noise
        const u1 = rand(), u2 = rand();
        const noise = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        pc1[s] = c + drift * noise;
      }
      windows.push({ pca: { pc1 } });
    }
    _state.data = {
      chrom: 'TEST',
      n_windows,
      windows,
      samples: Array.from({ length: n_samples }, (_, i) => ({ id: `S${i}` })),
    };
    return { _state, labels, n_samples, n_windows };
  }

  // ============================================================
  // autoSelectVocabForCandidate
  // ============================================================
  t('vocab: K=3 -> standard', () => {
    eq(w.autoSelectVocabForCandidate({ K: 3 }), 'standard');
  });
  t('vocab: K=2 -> standard', () => {
    eq(w.autoSelectVocabForCandidate({ K: 2 }), 'standard');
  });
  t('vocab: K=4 -> multi2', () => {
    eq(w.autoSelectVocabForCandidate({ K: 4 }), 'multi2');
  });
  t('vocab: K=6 -> multi3', () => {
    eq(w.autoSelectVocabForCandidate({ K: 6 }), 'multi3');
  });
  t('vocab: K=8 -> free', () => {
    eq(w.autoSelectVocabForCandidate({ K: 8 }), 'free');
  });

  // ============================================================
  // get/set vocab persistence
  // ============================================================
  t('vocab: set persists in localStorage', () => {
    clearLS();
    const c = { id: 'cand_v', chrom: 'LG28', K: 3 };
    w.setHaplotypeVocab(c, 'binary');
    eq(w.getHaplotypeVocab(c), 'binary');
    // Reload state to simulate fresh page
    delete w.state.haplotypeVocabs;
    eq(w.getHaplotypeVocab(c), 'binary');   // restored from LS
  });

  t('vocab: getHaplotypeVocab falls back to auto-select if not set', () => {
    clearLS();
    delete w.state.haplotypeVocabs;
    const c6 = { id: 'cand_v6', chrom: 'LG28', K: 6 };
    eq(w.getHaplotypeVocab(c6), 'multi3');
  });

  // ============================================================
  // _candPerBandStats
  // ============================================================
  t('per-band stats: returns K entries with n, centroid, mean_sigma', () => {
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.1 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },   // het: high drift
      { n_samples: 30, centroid: 1,  drift: 0.1 },
    ]);
    const c = { id: 'cps', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['cps'] = c;
    const stats = w._candPerBandStats(c);
    eq(stats.length, 3);
    eq(stats[0].n, 30);
    eq(stats[1].n, 30);
    eq(stats[2].n, 30);
    // band 0 centroid should be near -1
    if (Math.abs(stats[0].centroid_pc1 - (-1)) > 0.2) throw new Error(`band 0 centroid wrong: ${stats[0].centroid_pc1}`);
    if (Math.abs(stats[2].centroid_pc1 - 1) > 0.2) throw new Error(`band 2 centroid wrong: ${stats[2].centroid_pc1}`);
    // Mean sigma should be highest in band 1 (high drift)
    if (stats[1].mean_sigma < stats[0].mean_sigma) throw new Error('band 1 should have highest sigma');
    if (stats[1].mean_sigma < stats[2].mean_sigma) throw new Error('band 1 should have highest sigma');
  });

  // ============================================================
  // autoClassifyCandidate — vocab=standard with clear het signal
  // ============================================================
  t('classify: K=3 standard vocab with clear het -> HOM_REF / HET / HOM_INV high confidence', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.05 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },
      { n_samples: 30, centroid: 1,  drift: 0.05 },
    ]);
    const c = { id: 'std', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['std'] = c;
    const r = w.autoClassifyCandidate(c);
    eq(r.vocab, 'standard');
    eq(r.per_band[0].label, 'HOM_REF');
    eq(r.per_band[1].label, 'HET');
    eq(r.per_band[2].label, 'HOM_INV');
    eq(r.per_band[0].confidence, 'high');
    eq(r.per_band[1].confidence, 'high');
    eq(r.per_band[2].confidence, 'high');
  });

  t('classify: K=3 standard vocab with NO het signal -> medium/low confidence', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.1 },
      { n_samples: 30, centroid: 0,  drift: 0.1 },
      { n_samples: 30, centroid: 1,  drift: 0.1 },
    ]);
    const c = { id: 'std_nohet', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['std_nohet'] = c;
    const r = w.autoClassifyCandidate(c);
    eq(r.vocab, 'standard');
    // Without clear het, confidence is medium for outer bands and low for the middle
    if (r.per_band[1].confidence === 'high') throw new Error('middle band should not be high without het signal');
  });

  t('classify: small band (n<5) -> RECOMBINANT', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.1 },
      { n_samples: 2,  centroid: 0,  drift: 0.1 },   // tiny band
      { n_samples: 30, centroid: 1,  drift: 0.1 },
    ]);
    const c = { id: 'rec', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['rec'] = c;
    const r = w.autoClassifyCandidate(c);
    eq(r.per_band[1].label, 'RECOMBINANT');
    if (r.per_band[1].confidence !== 'medium') throw new Error('expected medium confidence for recombinant');
  });

  // ============================================================
  // autoClassifyCandidate — vocab=binary
  // ============================================================
  t('classify: vocab binary -> AA / AB / BB', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.05 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },
      { n_samples: 30, centroid: 1,  drift: 0.05 },
    ]);
    const c = { id: 'bin', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['bin'] = c;
    w.setHaplotypeVocab(c, 'binary');
    const r = w.autoClassifyCandidate(c);
    eq(r.vocab, 'binary');
    eq(r.per_band[0].label, 'AA');
    eq(r.per_band[1].label, 'AB');
    eq(r.per_band[2].label, 'BB');
  });

  // ============================================================
  // autoClassifyCandidate — vocab=multi3 K=6
  // ============================================================
  t('classify: K=6 multi3 vocab returns 6 labels with low confidence', () => {
    clearLS();
    const spec = [];
    for (let k = 0; k < 6; k++) {
      spec.push({ n_samples: 20, centroid: k * 0.5, drift: 0.05 });
    }
    const { _state, labels } = makeState(6, spec);
    const c = { id: 'm3', chrom: 'TEST', K: 6, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['m3'] = c;
    w.setHaplotypeVocab(c, 'multi3');
    const r = w.autoClassifyCandidate(c);
    eq(r.vocab, 'multi3');
    eq(r.per_band.length, 6);
    // All labels should be from the multi3 vocab (or RECOMBINANT)
    const validLabels = ['H1/H1','H1/H2','H1/H3','H2/H2','H2/H3','H3/H3','RECOMBINANT'];
    for (const pb of r.per_band) {
      if (!validLabels.includes(pb.label)) throw new Error(`unexpected label: ${pb.label}`);
    }
    // K=6 multi-haplotype is genuinely ambiguous; most should be 'low' or 'medium'
    const highCount = r.per_band.filter(pb => pb.confidence === 'high').length;
    if (highCount > 1) throw new Error(`expected at most 1 high confidence in K=6, got ${highCount}`);
  });

  // ============================================================
  // applyAutoClassificationToCandidate
  // ============================================================
  t('apply: fills empty bands, leaves user input alone', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.05 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },
      { n_samples: 30, centroid: 1,  drift: 0.05 },
    ]);
    const c = { id: 'apply_test', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['apply_test'] = c;
    // User pre-set band 1
    w.setHaplotypeLabel(c, 1, 'CUSTOM_USER');
    const r = w.applyAutoClassificationToCandidate(c);
    // Bands 0 and 2 should be auto-filled
    eq(c.haplotype_labels[0], 'HOM_REF');
    eq(c.haplotype_labels[2], 'HOM_INV');
    // Band 1 stays untouched
    eq(c.haplotype_labels[1], 'CUSTOM_USER');
    if (r.applied !== 2) throw new Error(`expected 2 applied, got ${r.applied}`);
  });

  t('apply: overwriteUserInput=true overwrites existing labels', () => {
    clearLS();
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.05 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },
      { n_samples: 30, centroid: 1,  drift: 0.05 },
    ]);
    const c = { id: 'apply_ov', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['apply_ov'] = c;
    w.setHaplotypeLabel(c, 0, 'OLD');
    w.applyAutoClassificationToCandidate(c, { overwriteUserInput: true });
    eq(c.haplotype_labels[0], 'HOM_REF');
  });

  t('apply: confidenceFloor="high" only applies high-confidence labels', () => {
    clearLS();
    // K=3 with NO clear het — only outer bands are medium, middle is low
    const { _state, labels } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.1 },
      { n_samples: 30, centroid: 0,  drift: 0.1 },
      { n_samples: 30, centroid: 1,  drift: 0.1 },
    ]);
    const c = { id: 'apply_high', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['apply_high'] = c;
    const r = w.applyAutoClassificationToCandidate(c, { confidenceFloor: 'high' });
    // None should be applied since no band reached 'high' without het
    eq(r.applied, 0);
  });

  // ============================================================
  // applyAutoClassificationToAllCandidates
  // ============================================================
  t('apply-all: walks every candidate', () => {
    clearLS();
    const { _state } = makeState(3, [
      { n_samples: 30, centroid: -1, drift: 0.05 },
      { n_samples: 30, centroid: 0,  drift: 0.5 },
      { n_samples: 30, centroid: 1,  drift: 0.05 },
    ]);
    // Replace synthetic data's single candidate with two candidates sharing same labels
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    _state.candidates = {
      'A': { id: 'A', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 },
      'B': { id: 'B', chrom: 'TEST', K: 3, locked_labels: labels, start_w: 0, end_w: 29 },
    };
    const r = w.applyAutoClassificationToAllCandidates();
    eq(Object.keys(r.perCandidate).length, 2);
    eq(r.totalApplied, 6);   // 3 bands × 2 candidates
    eq(_state.candidates['A'].haplotype_labels[0], 'HOM_REF');
    eq(_state.candidates['B'].haplotype_labels[2], 'HOM_INV');
  });

  // ============================================================
  // Realistic LG28-style scenario
  // ============================================================
  t('end-to-end: K=6 LG28-style auto-fill produces sensible labels', () => {
    clearLS();
    // K=6 with the structure: 2 stacked systems, each with 3 (HOM/HET/HOM) bands
    // We won't get them perfect but the classifier should produce SOMETHING
    // and the outputs should be from the multi3 vocabulary.
    const spec = [];
    // First system (bands 0, 1, 2 by PC1)
    spec.push({ n_samples: 40, centroid: -1.5, drift: 0.05 });
    spec.push({ n_samples: 38, centroid: -0.5, drift: 0.4  });   // HET-like
    spec.push({ n_samples: 38, centroid:  0.5, drift: 0.05 });
    // Second system (bands 3, 4, 5 by PC1)
    spec.push({ n_samples: 38, centroid:  1.5, drift: 0.05 });
    spec.push({ n_samples: 36, centroid:  2.5, drift: 0.05 });
    spec.push({ n_samples: 36, centroid:  3.5, drift: 0.05 });
    const { _state, labels } = makeState(6, spec);
    const c = { id: 'lg28', chrom: 'C_gar_LG28', K: 6, locked_labels: labels, start_w: 0, end_w: 29 };
    _state.candidates['lg28'] = c;
    const r = w.autoClassifyCandidate(c);
    eq(r.vocab, 'multi3');
    eq(r.per_band.length, 6);
    // All labels should be valid multi3 vocab values
    const validLabels = ['H1/H1','H1/H2','H1/H3','H2/H2','H2/H3','H3/H3','RECOMBINANT'];
    for (const pb of r.per_band) {
      if (!validLabels.includes(pb.label)) throw new Error(`bad label ${pb.label}`);
    }
  });

  return failures;
}

setTimeout(() => {
  const failures = run();
  if (failures.length > 0) {
    console.log(`\n${failures.length} test(s) failed`);
    failures.forEach(f => console.log(`  - ${f.name}: ${f.err}`));
    process.exit(1);
  } else {
    console.log('\nAll tests passed');
    process.exit(0);
  }
}, 500);
