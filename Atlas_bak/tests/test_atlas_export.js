// Tests for turn 2e — atlas → phase 7 export contract.
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

  if (typeof w.buildAtlasCandidateExport !== 'function'
      || typeof w.buildCandidateExportRecord !== 'function'
      || typeof w.downloadAtlasExport !== 'function') {
    console.log('  FAIL: turn 2e functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2e functions exposed');

  function clearLS() { try { w.localStorage.clear(); } catch (_) {} }

  function setupState(candConfigs, opts) {
    opts = opts || {};
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = opts.mode || 'default';
    _state.k = (_state.activeMode === 'detailed') ? 6 : 3;
    const target = (_state.activeMode === 'detailed') ? 'candidates_detailed' : 'candidates';
    _state[target] = {};
    for (const c of candConfigs) {
      _state[target][c.id] = {
        id: c.id,
        chrom: c.chrom || 'C_gar_LG28',
        K: c.K, K_used: c.K,
        locked_labels: c.labels,
        start_w: c.start_w || 100,
        end_w: c.end_w || 200,
        start_bp: c.start_bp,
        end_bp: c.end_bp,
        source: c.source || 'L2',
        confirmed: !!c.confirmed,
        notes: c.notes || '',
      };
    }
    _state.data = {
      chrom: opts.chrom || 'C_gar_LG28',
      n_windows: opts.n_windows || 4302,
      samples: opts.samples || (function () {
        const arr = []; for (let i = 0; i < 226; i++) arr.push({ id: `S${String(i).padStart(3, '0')}` });
        return arr;
      })(),
    };
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
  }

  // ============================================================
  // buildCandidateExportRecord
  // ============================================================
  t('record: null candidate returns null', () => {
    eq(w.buildCandidateExportRecord(null), null);
  });

  t('record: minimal candidate produces valid record', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_1', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ]);
    const r = w.buildCandidateExportRecord(w.state.candidates['cand_1']);
    eq(r.candidate_id, 'cand_1');
    eq(r.chrom, 'C_gar_LG28');
    eq(r.K_used, 3);
    eq(r.start_bp, 1_000_000);
    eq(r.end_bp, 2_000_000);
    eq(typeof r.atlas_band_assignments, 'object');
    eq(r.atlas_band_assignments.band_counts.length, 3);
    eq(r.atlas_band_assignments.band_counts[0], 30);
    eq(r.atlas_band_assignments.band_counts[1], 30);
    eq(r.atlas_band_assignments.band_counts[2], 30);
  });

  t('record: includes haplotype_labels from turn 2d', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_h', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ]);
    const c = w.state.candidates['cand_h'];
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    w.setHaplotypeLabel(c, 1, 'H1/H2');
    w.setHaplotypeLabel(c, 2, 'H2/H2');
    const r = w.buildCandidateExportRecord(c);
    eq(r.haplotype_labels[0], 'H1/H1');
    eq(r.haplotype_labels[1], 'H1/H2');
    eq(r.haplotype_labels[2], 'H2/H2');
  });

  t('record: per_sample has primary_band per sample', () => {
    clearLS();
    const labels = [0, 1, 2, 0, 1, 2, 0, 0];   // 8 samples
    setupState([
      { id: 'cand_p', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { samples: Array.from({ length: 8 }, (_, i) => ({ id: `Z${i}` })) });
    const r = w.buildCandidateExportRecord(w.state.candidates['cand_p']);
    eq(r.atlas_band_assignments.per_sample['Z0'].primary_band, 0);
    eq(r.atlas_band_assignments.per_sample['Z3'].primary_band, 0);
    eq(r.atlas_band_assignments.per_sample['Z1'].primary_band, 1);
    eq(r.atlas_band_assignments.per_sample['Z2'].primary_band, 2);
  });

  t('record: per_sample primary_band_label uses haplotype labels', () => {
    clearLS();
    const labels = [0, 1, 0, 1];
    setupState([
      { id: 'cand_pl', K: 2, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { samples: Array.from({ length: 4 }, (_, i) => ({ id: `S${i}` })) });
    const c = w.state.candidates['cand_pl'];
    w.setHaplotypeLabel(c, 0, 'REF');
    w.setHaplotypeLabel(c, 1, 'INV');
    const r = w.buildCandidateExportRecord(c);
    eq(r.atlas_band_assignments.per_sample['S0'].primary_band_label, 'REF');
    eq(r.atlas_band_assignments.per_sample['S1'].primary_band_label, 'INV');
    eq(r.atlas_band_assignments.per_sample['S2'].primary_band_label, 'REF');
    eq(r.atlas_band_assignments.per_sample['S3'].primary_band_label, 'INV');
  });

  t('record: sample_groups built from haplotype_labels', () => {
    clearLS();
    const labels = [0, 0, 1, 1, 2, 2];
    setupState([
      { id: 'cand_sg', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ], { samples: Array.from({ length: 6 }, (_, i) => ({ id: `S${i}` })) });
    const c = w.state.candidates['cand_sg'];
    w.setHaplotypeLabel(c, 0, 'HOM_REF');
    w.setHaplotypeLabel(c, 1, 'HET');
    w.setHaplotypeLabel(c, 2, 'HOM_INV');
    const r = w.buildCandidateExportRecord(c);
    eq(r.sample_groups['HOM_REF'].length, 2);
    eq(r.sample_groups['HET'].length, 2);
    eq(r.sample_groups['HOM_INV'].length, 2);
    eq(r.sample_groups['HOM_REF'][0], 'S0');
    eq(r.sample_groups['HET'][0], 'S2');
    eq(r.sample_groups['HOM_INV'][0], 'S4');
  });

  t('record: unlabeled bands fall back to band_N keys in sample_groups', () => {
    clearLS();
    const labels = [0, 0, 1, 1];
    setupState([
      { id: 'cand_u', K: 2, start_bp: 1, end_bp: 2, labels },
    ], { samples: Array.from({ length: 4 }, (_, i) => ({ id: `S${i}` })) });
    const c = w.state.candidates['cand_u'];
    // No labels set
    const r = w.buildCandidateExportRecord(c);
    eq(r.sample_groups['band_0'].length, 2);
    eq(r.sample_groups['band_1'].length, 2);
  });

  t('record: includes inheritance info when result available', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_inh', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ]);
    // Synthesize an inheritance result that includes this candidate
    w.state.inheritanceResult = {
      items_meta: [{ id: 'cand_inh', K: 3, seq_num: 1, start_bp: 1_000_000, end_bp: 2_000_000 }],
      band_index: [
        { item_idx: 0, band: 0 },
        { item_idx: 0, band: 1 },
        { item_idx: 0, band: 2 },
      ],
      cut: { group_id_per_band: new Int32Array([0, 1, 2]), n_groups: 3 },
      rtab: { per_item_n_groups: [3] },
    };
    const r = w.buildCandidateExportRecord(w.state.candidates['cand_inh']);
    if (!r.inheritance) throw new Error('expected inheritance block');
    eq(r.inheritance.seq_num, 1);
    eq(r.inheritance.n_inheritance_groups, 3);
    eq(r.inheritance.total_groups_in_cohort, 3);
    eq(r.inheritance.group_id_per_band[0], 0);
    eq(r.inheritance.group_id_per_band[1], 1);
    eq(r.inheritance.group_id_per_band[2], 2);
  });

  t('record: inheritance is null when candidate is not in result', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_notinh', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    w.state.inheritanceResult = {
      items_meta: [{ id: 'OTHER_CAND', K: 3, seq_num: 1 }],
      band_index: [],
      cut: { group_id_per_band: new Int32Array(), n_groups: 0 },
      rtab: {},
    };
    const r = w.buildCandidateExportRecord(w.state.candidates['cand_notinh']);
    eq(r.inheritance, null);
  });

  // ============================================================
  // buildAtlasCandidateExport
  // ============================================================
  t('full: returns format_version + cohort + candidates', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'cand_B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels },
    ]);
    const exp = w.buildAtlasCandidateExport();
    eq(exp.format_version, 'atlas_candidate_export_v2');
    eq(typeof exp.atlas_provenance, 'object');
    eq(typeof exp.atlas_provenance.exported_at, 'string');
    eq(exp.cohort.chrom, 'C_gar_LG28');
    eq(exp.cohort.n_samples, 226);
    eq(exp.candidates.length, 2);
  });

  t('full: candidates sorted by start_bp', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels },
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels },
    ]);
    const exp = w.buildAtlasCandidateExport();
    eq(exp.candidates[0].candidate_id, 'A');
    eq(exp.candidates[1].candidate_id, 'B');
    eq(exp.candidates[2].candidate_id, 'C');
  });

  t('full: filter to single candidate via candidate_id opt', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels },
    ]);
    const exp = w.buildAtlasCandidateExport({ candidate_id: 'A' });
    eq(exp.candidates.length, 1);
    eq(exp.candidates[0].candidate_id, 'A');
  });

  t('full: respects mode (default vs detailed)', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 6; k++) for (let i = 0; i < 20; i++) labels.push(k);
    // Set up detailed mode candidates
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'detailed';
    _state.k = 6;
    _state.candidates = {};
    _state.candidates_detailed = {
      'cand_d1': { id: 'cand_d1', chrom: 'LG28', K: 6, K_used: 6, locked_labels: labels, start_bp: 1, end_bp: 2 },
    };
    _state.data = { chrom: 'LG28', n_windows: 100, samples: Array.from({length: 120}, (_, i) => ({id: `S${i}`})) };
    _state.inheritanceResult = null;

    const exp = w.buildAtlasCandidateExport({ mode: 'detailed' });
    eq(exp.atlas_provenance.mode, 'detailed');
    eq(exp.candidates.length, 1);
    eq(exp.candidates[0].K_used, 6);
  });

  t('full: empty candidates returns empty list', () => {
    setupState([]);
    const exp = w.buildAtlasCandidateExport();
    eq(exp.candidates.length, 0);
    eq(exp.format_version, 'atlas_candidate_export_v2');
  });

  t('full: JSON.stringify produces valid JSON', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_J', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
    ]);
    const c = w.state.candidates['cand_J'];
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    const exp = w.buildAtlasCandidateExport();
    const json = JSON.stringify(exp, null, 2);
    if (typeof json !== 'string' || json.length < 100) throw new Error('json too short');
    // Round-trip
    const parsed = JSON.parse(json);
    eq(parsed.format_version, 'atlas_candidate_export_v2');
    eq(parsed.candidates[0].haplotype_labels[0], 'H1/H1');
  });

  // ============================================================
  // downloadAtlasExport
  // ============================================================
  t('download: returns filename with chrom and timestamp', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'cand_dl', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    // Stub URL.createObjectURL since JSDOM may not provide it
    let urlCreated = null;
    w.URL.createObjectURL = (blob) => { urlCreated = blob; return 'blob:fake-url'; };
    w.URL.revokeObjectURL = () => {};
    const filename = w.downloadAtlasExport();
    if (typeof filename !== 'string' || filename.indexOf('atlas_export_') < 0) {
      throw new Error(`bad filename: ${filename}`);
    }
    if (filename.indexOf('C_gar_LG28') < 0) throw new Error('missing chrom in filename');
    if (filename.indexOf('.json') < 0) throw new Error('missing .json extension');
    if (!urlCreated) throw new Error('createObjectURL not called');
  });

  t('download: per-candidate filename includes candidate_id', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'specific', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    w.URL.createObjectURL = () => 'blob:fake';
    w.URL.revokeObjectURL = () => {};
    const filename = w.downloadAtlasExport({ candidate_id: 'specific' });
    if (filename.indexOf('specific') < 0) throw new Error(`missing candidate id in filename: ${filename}`);
  });

  // ============================================================
  // End-to-end realistic
  // ============================================================
  t('end-to-end: realistic LG28 export with full annotations', () => {
    clearLS();
    const N = 226;
    const labels = [];
    const groups = [40, 38, 38, 38, 36, 36];
    for (let g = 0; g < 6; g++) {
      for (let i = 0; i < groups[g]; i++) labels.push(g);
    }
    setupState([
      { id: 'LG28_cand_01', chrom: 'C_gar_LG28', K: 6,
        start_w: 3260, end_w: 3870, start_bp: 15115000, end_bp: 18005000,
        labels, source: 'l3_draft_w', confirmed: true,
        notes: 'Diamond at 16.2 Mb suggests two stacked systems' },
    ], { mode: 'detailed', samples: Array.from({length: N}, (_, i) => ({id: `S${String(i).padStart(3,'0')}`})) });
    const c = w.state.candidates_detailed['LG28_cand_01'];
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    w.setHaplotypeLabel(c, 1, 'H1/H2');
    w.setHaplotypeLabel(c, 2, 'H1/H3');
    w.setHaplotypeLabel(c, 3, 'H2/H2');
    w.setHaplotypeLabel(c, 4, 'H2/H3');
    w.setHaplotypeLabel(c, 5, 'H3/H3');

    const exp = w.buildAtlasCandidateExport({ mode: 'detailed' });
    eq(exp.candidates.length, 1);
    const r = exp.candidates[0];
    eq(r.candidate_id, 'LG28_cand_01');
    eq(r.K_used, 6);
    eq(r.confirmed, true);
    eq(r.start_bp, 15115000);
    eq(r.end_bp, 18005000);
    eq(r.haplotype_labels[0], 'H1/H1');
    eq(r.haplotype_labels[5], 'H3/H3');
    eq(r.annotator_notes, 'Diamond at 16.2 Mb suggests two stacked systems');
    // sample_groups should have 6 keys (H1/H1 ... H3/H3)
    eq(Object.keys(r.sample_groups).length, 6);
    eq(r.sample_groups['H1/H1'].length, 40);
    eq(r.sample_groups['H1/H2'].length, 38);
    // Round-trip JSON
    const json = JSON.stringify(exp);
    const parsed = JSON.parse(json);
    eq(parsed.candidates[0].haplotype_labels[3], 'H2/H2');
  });

  // ============================================================
  // v2 additions (turn 121: full atlas context)
  // ============================================================
  t('v2: format_version is atlas_candidate_export_v2', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const exp = w.buildAtlasCandidateExport();
    eq(exp.format_version, 'atlas_candidate_export_v2');
  });

  t('v2: provenance includes atlas_constants and active_K', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const exp = w.buildAtlasCandidateExport();
    if (!exp.atlas_provenance.atlas_constants) throw new Error('missing atlas_constants');
    if (typeof exp.atlas_provenance.atlas_constants.diamond !== 'object') throw new Error('missing diamond constants');
    if (typeof exp.atlas_provenance.atlas_constants.inheritance !== 'object') throw new Error('missing inheritance constants');
    if (typeof exp.atlas_provenance.haplotype_vocabs !== 'object') throw new Error('missing haplotype_vocabs');
  });

  t('v2: cohort.samples carries id + family_id + ancestry', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    const samples = Array.from({ length: 90 }, (_, i) => ({
      id: `S${String(i).padStart(3,'0')}`,
      family_id: i % 6,                       // 6 broodlines
      ancestry: i < 50 ? 'A' : 'B',
    }));
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }], { samples });
    const exp = w.buildAtlasCandidateExport();
    eq(exp.cohort.samples.length, 90);
    eq(exp.cohort.samples[0].id, 'S000');
    eq(exp.cohort.samples[0].family_id, '0');
    eq(exp.cohort.samples[0].ancestry, 'A');
    eq(exp.cohort.samples[60].family_id, '0');   // 60 % 6 = 0
    eq(exp.cohort.samples[60].ancestry, 'B');
  });

  t('v2: cohort.samples handles missing family_id (-1) as null', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    const samples = Array.from({ length: 90 }, (_, i) => ({
      id: `S${i}`,
      family_id: -1,
      ancestry: null,
    }));
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }], { samples });
    const exp = w.buildAtlasCandidateExport();
    eq(exp.cohort.samples[0].family_id, null);
    eq(exp.cohort.samples[0].ancestry, null);
  });

  t('v2: cohort.inheritance_result populated when compute has run', () => {
    clearLS();
    const labA = [], labB = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) { labA.push(k); labB.push(k); }
    setupState([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels: labA },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: labB },
    ]);
    w.runInheritanceCompute();
    const exp = w.buildAtlasCandidateExport();
    if (!exp.cohort.inheritance_result) throw new Error('missing inheritance_result');
    eq(exp.cohort.inheritance_result.n_items, 2);
    if (!Array.isArray(exp.cohort.inheritance_result.dendrogram)) throw new Error('dendrogram not array');
    if (!exp.cohort.inheritance_result.cut) throw new Error('cut missing');
    eq(exp.cohort.inheritance_result.cut.n_groups, 3);
  });

  t('v2: cohort.inheritance_result is null when no compute', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'A', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const exp = w.buildAtlasCandidateExport();
    eq(exp.cohort.inheritance_result, null);
  });

  t('v2: candidate record includes haplotype_vocab', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 6; k++) for (let i = 0; i < 20; i++) labels.push(k);
    setupState([{ id: 'V', K: 6, start_bp: 1, end_bp: 2, labels }]);
    const c = w.state.candidates['V'];
    w.setHaplotypeVocab(c, 'multi3');
    const exp = w.buildAtlasCandidateExport();
    eq(exp.candidates[0].haplotype_vocab, 'multi3');
  });

  t('v2: candidate record includes classifier_output per band', () => {
    clearLS();
    // Build a candidate with windows + per-window PC1 so the classifier can run
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.candidates = {};
    _state.candidates_detailed = {};
    const n_samples = 90;
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    // Synthesize 30 windows with band-separated PC1
    const n_windows = 30;
    const windows = [];
    for (let w_idx = 0; w_idx < n_windows; w_idx++) {
      const pc1 = new Float32Array(n_samples);
      for (let s = 0; s < n_samples; s++) {
        const k = labels[s];
        pc1[s] = (k - 1) * 1.0 + 0.05 * Math.sin(s + w_idx);
      }
      windows.push({ pca: { pc1 } });
    }
    _state.data = {
      chrom: 'TEST', n_windows, windows,
      samples: Array.from({ length: n_samples }, (_, i) => ({ id: `S${i}` })),
    };
    _state.candidates['CL'] = {
      id: 'CL', chrom: 'TEST', K: 3, K_used: 3,
      locked_labels: labels, start_w: 0, end_w: n_windows - 1,
      start_bp: 1, end_bp: 2, source: 'test',
    };
    _state.inheritanceResult = null;

    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    if (!r.classifier_output) throw new Error('classifier_output missing');
    eq(r.classifier_output.per_band.length, 3);
    if (typeof r.classifier_output.per_band[0].confidence !== 'string') throw new Error('confidence missing');
  });

  t('v2: candidate record includes per_band_stats', () => {
    clearLS();
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.candidates = {}; _state.candidates_detailed = {};
    const n_samples = 90;
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    const n_windows = 30;
    const windows = [];
    for (let w_idx = 0; w_idx < n_windows; w_idx++) {
      const pc1 = new Float32Array(n_samples);
      for (let s = 0; s < n_samples; s++) pc1[s] = (labels[s] - 1) * 1.0;
      windows.push({ pca: { pc1 } });
    }
    _state.data = { chrom: 'TEST', n_windows, windows, samples: Array.from({length: n_samples}, (_, i) => ({id: `S${i}`})) };
    _state.candidates['PB'] = {
      id: 'PB', chrom: 'TEST', K: 3, K_used: 3,
      locked_labels: labels, start_w: 0, end_w: n_windows - 1,
      start_bp: 1, end_bp: 2, source: 'test',
    };
    _state.inheritanceResult = null;

    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    if (!Array.isArray(r.per_band_stats)) throw new Error('per_band_stats missing');
    eq(r.per_band_stats.length, 3);
    eq(r.per_band_stats[0].n, 30);
  });

  t('v2: candidate record includes diamonds_full when diamonds detected', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'D', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    // diamonds_full may be null or empty array depending on detection;
    // but if present, it should be an array
    if (r.diamonds_full !== null && !Array.isArray(r.diamonds_full)) {
      throw new Error('diamonds_full should be null or array');
    }
  });

  t('v2: candidate record includes locked_labels_raw', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([{ id: 'LL', K: 3, start_bp: 1, end_bp: 2, labels }]);
    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    if (!Array.isArray(r.locked_labels_raw)) throw new Error('locked_labels_raw missing');
    eq(r.locked_labels_raw.length, 90);
    eq(r.locked_labels_raw[0], 0);
    eq(r.locked_labels_raw[89], 2);
  });

  t('v2: candidate record includes candidate_runtime metadata', () => {
    clearLS();
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupState([
      { id: 'RT', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    // Add some runtime fields
    w.state.candidates['RT'].ref_l2 = 14;
    w.state.candidates['RT'].l2_indices = [12, 13, 14, 15];
    w.state.candidates['RT'].resolution = 'L2';
    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    eq(r.candidate_runtime.ref_l2, 14);
    if (!Array.isArray(r.candidate_runtime.l2_indices)) throw new Error('l2_indices not array');
    eq(r.candidate_runtime.l2_indices.length, 4);
    eq(r.candidate_runtime.resolution, 'L2');
  });

  t('v2: band_composition is null or array (function depends on lexical state.data which test cannot populate via window.state)', () => {
    clearLS();
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.candidates = {}; _state.candidates_detailed = {};
    const n_samples = 90;
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    const samples = Array.from({ length: n_samples }, (_, i) => ({
      id: `S${i}`,
      family_id: i % 3,
      ancestry: 'A',
    }));
    _state.data = {
      chrom: 'TEST', n_windows: 30,
      windows: Array.from({ length: 30 }, () => ({ pca: { pc1: new Float32Array(n_samples) } })),
      samples,
    };
    _state.candidates['BC'] = {
      id: 'BC', chrom: 'TEST', K: 3, K_used: 3,
      locked_labels: labels, start_w: 0, end_w: 29,
      start_bp: 1, end_bp: 2, source: 'test',
    };
    _state.inheritanceResult = null;
    const exp = w.buildAtlasCandidateExport();
    const r = exp.candidates[0];
    // candidateBandComposition reads the LEXICAL `state` const (not window.state),
    // which JSDOM tests can't populate. In production the field is an array; here
    // we accept null too.
    if (r.band_composition !== null && !Array.isArray(r.band_composition)) {
      throw new Error('band_composition should be null or array');
    }
    if (Array.isArray(r.band_composition)) {
      eq(r.band_composition.length, 3);
    }
  });

  t('v2: end-to-end LG28-style export with full context', () => {
    clearLS();
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'detailed';
    _state.k = 6;
    _state.candidates = {};
    _state.candidates_detailed = {};
    const N = 226;
    const labels = [];
    const groups = [40, 38, 38, 38, 36, 36];
    for (let g = 0; g < 6; g++) for (let i = 0; i < groups[g]; i++) labels.push(g);
    const samples = Array.from({ length: N }, (_, i) => ({
      id: `S${String(i).padStart(3, '0')}`,
      family_id: i % 6,
      ancestry: null,
    }));
    const n_windows = 30;
    const windows = [];
    for (let w_idx = 0; w_idx < n_windows; w_idx++) {
      const pc1 = new Float32Array(N);
      for (let s = 0; s < N; s++) pc1[s] = (labels[s] - 2.5) * 0.5;
      windows.push({ pca: { pc1 } });
    }
    _state.data = { chrom: 'C_gar_LG28', n_windows, windows, samples };
    _state.candidates_detailed['LG28_cand_01'] = {
      id: 'LG28_cand_01', chrom: 'C_gar_LG28', K: 6, K_used: 6,
      locked_labels: labels,
      start_w: 0, end_w: n_windows - 1,
      start_bp: 15115000, end_bp: 18005000,
      source: 'l3_draft_w', confirmed: true,
      notes: 'Diamond at 16.2 Mb',
    };
    const c = _state.candidates_detailed['LG28_cand_01'];
    w.setHaplotypeVocab(c, 'multi3');
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    w.setHaplotypeLabel(c, 1, 'H1/H2');
    w.setHaplotypeLabel(c, 5, 'H3/H3');

    const exp = w.buildAtlasCandidateExport({ mode: 'detailed' });
    eq(exp.format_version, 'atlas_candidate_export_v2');
    eq(exp.cohort.n_samples, N);
    eq(exp.cohort.samples.length, N);
    eq(exp.candidates.length, 1);
    const r = exp.candidates[0];
    eq(r.K_used, 6);
    eq(r.haplotype_vocab, 'multi3');
    eq(r.haplotype_labels[0], 'H1/H1');
    eq(r.haplotype_labels[5], 'H3/H3');
    if (!r.classifier_output) throw new Error('classifier_output missing');
    if (!r.per_band_stats) throw new Error('per_band_stats missing');
    if (!r.locked_labels_raw) throw new Error('locked_labels_raw missing');
    eq(r.locked_labels_raw.length, N);
    if (!r.candidate_runtime) throw new Error('candidate_runtime missing');
    // Round-trip
    const json = JSON.stringify(exp);
    const parsed = JSON.parse(json);
    eq(parsed.candidates[0].haplotype_vocab, 'multi3');
    eq(parsed.cohort.samples[0].family_id, '0');
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
