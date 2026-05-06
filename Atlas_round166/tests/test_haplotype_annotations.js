// Tests for turn 2d — haplotype annotations panel.
const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.resolve(__dirname, 'atlas.html'), 'utf8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: new (require('jsdom').VirtualConsole)(),
  url: 'http://localhost/',   // needed for localStorage
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

  if (typeof w.candidateHaplotypeAnnotationsHtml !== 'function'
      || typeof w._wireCandidateHaplotypeAnnotations !== 'function'
      || typeof w.setHaplotypeLabel !== 'function'
      || typeof w.loadHaplotypeLabels !== 'function'
      || typeof w._inheritanceSuggestionsForCandidate !== 'function') {
    console.log('  FAIL: turn 2d functions not exposed');
    return [{ name: 'fns exposed', err: 'missing' }];
  }
  console.log('  All turn 2d functions exposed');

  // Clear localStorage between tests
  function clearLS() { try { w.localStorage.clear(); } catch (_) {} }

  // ============================================================
  // _hapLabelLSKey
  // ============================================================
  t('LS key: includes chrom and candidate id', () => {
    const k = w._hapLabelLSKey('LG28', 'cand_1');
    if (k.indexOf('LG28') < 0) throw new Error(`missing chrom: ${k}`);
    if (k.indexOf('cand_1') < 0) throw new Error(`missing id: ${k}`);
  });

  t('LS key: handles missing chrom/id', () => {
    const k = w._hapLabelLSKey(null, null);
    if (typeof k !== 'string' || k.length === 0) throw new Error(`bad key: ${k}`);
  });

  // ============================================================
  // loadHaplotypeLabels / setHaplotypeLabel
  // ============================================================
  t('load: returns empty {} for fresh candidate', () => {
    clearLS();
    const c = { id: 'fresh', chrom: 'LG28', K: 3 };
    const labels = w.loadHaplotypeLabels(c);
    eq(typeof labels, 'object');
    eq(Object.keys(labels).length, 0);
  });

  t('set: writes and persists', () => {
    clearLS();
    const c = { id: 'cand_X', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    w.setHaplotypeLabel(c, 1, 'H1/H2');
    w.setHaplotypeLabel(c, 2, 'H2/H2');
    eq(c.haplotype_labels[0], 'H1/H1');
    eq(c.haplotype_labels[1], 'H1/H2');
    eq(c.haplotype_labels[2], 'H2/H2');
    // Persisted in localStorage
    const key = w._hapLabelLSKey(c.chrom, c.id);
    const raw = w.localStorage.getItem(key);
    if (!raw) throw new Error('localStorage not persisted');
    const parsed = JSON.parse(raw);
    eq(parsed[0], 'H1/H1');
  });

  t('set: trim whitespace', () => {
    clearLS();
    const c = { id: 'cand_T', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(c, 0, '  H1/H1   ');
    eq(c.haplotype_labels[0], 'H1/H1');
  });

  t('set: empty value clears label', () => {
    clearLS();
    const c = { id: 'cand_E', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(c, 0, 'H1/H1');
    eq(c.haplotype_labels[0], 'H1/H1');
    w.setHaplotypeLabel(c, 0, '');
    if (c.haplotype_labels[0] !== undefined) throw new Error('expected deletion');
  });

  t('load: restores from localStorage on reload', () => {
    clearLS();
    const c1 = { id: 'cand_R', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(c1, 0, 'H1/H1');
    w.setHaplotypeLabel(c1, 1, 'H1/H2');
    // Simulate fresh candidate object (no haplotype_labels field)
    const c2 = { id: 'cand_R', chrom: 'LG28', K: 3 };
    const labels = w.loadHaplotypeLabels(c2);
    eq(labels[0], 'H1/H1');
    eq(labels[1], 'H1/H2');
  });

  t('load: separate candidates have separate stores', () => {
    clearLS();
    const cA = { id: 'cand_A', chrom: 'LG28', K: 3 };
    const cB = { id: 'cand_B', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(cA, 0, 'AA');
    w.setHaplotypeLabel(cB, 0, 'XX');
    const cA2 = { id: 'cand_A', chrom: 'LG28', K: 3 };
    const cB2 = { id: 'cand_B', chrom: 'LG28', K: 3 };
    eq(w.loadHaplotypeLabels(cA2)[0], 'AA');
    eq(w.loadHaplotypeLabels(cB2)[0], 'XX');
  });

  t('load: same candidate id on different chromosomes are separate', () => {
    clearLS();
    const cLG28 = { id: 'cand_1', chrom: 'LG28', K: 3 };
    const cLG12 = { id: 'cand_1', chrom: 'LG12', K: 3 };
    w.setHaplotypeLabel(cLG28, 0, 'LG28-band0');
    w.setHaplotypeLabel(cLG12, 0, 'LG12-band0');
    const cLG28_fresh = { id: 'cand_1', chrom: 'LG28', K: 3 };
    eq(w.loadHaplotypeLabels(cLG28_fresh)[0], 'LG28-band0');
  });

  // ============================================================
  // _inheritanceSuggestionsForCandidate
  // ============================================================
  t('suggestions: empty when no inheritance result', () => {
    if (!w.state) w.state = {};
    w.state.inheritanceResult = null;
    const c = { id: 'cand_1', chrom: 'LG28', K: 3 };
    const s = w._inheritanceSuggestionsForCandidate(c);
    eq(Object.keys(s).length, 0);
  });

  t('suggestions: returns group ids when candidate is in result', () => {
    if (!w.state) w.state = {};
    // Synthesize an inheritance result
    w.state.inheritanceResult = {
      items_meta: [
        { id: 'cand_A', K: 3, seq_num: 1 },
        { id: 'cand_B', K: 3, seq_num: 2 },
      ],
      band_index: [
        { item_idx: 0, band: 0 },
        { item_idx: 0, band: 1 },
        { item_idx: 0, band: 2 },
        { item_idx: 1, band: 0 },
        { item_idx: 1, band: 1 },
        { item_idx: 1, band: 2 },
      ],
      cut: {
        group_id_per_band: new Int32Array([0, 1, 2, 0, 1, 2]),
        n_groups: 3,
      },
      rtab: {},
    };
    const cA = { id: 'cand_A', chrom: 'LG28', K: 3 };
    const sA = w._inheritanceSuggestionsForCandidate(cA);
    eq(sA[0], 0); eq(sA[1], 1); eq(sA[2], 2);
    const cB = { id: 'cand_B', chrom: 'LG28', K: 3 };
    const sB = w._inheritanceSuggestionsForCandidate(cB);
    eq(sB[0], 0); eq(sB[1], 1); eq(sB[2], 2);
  });

  t('suggestions: candidate not in result returns empty', () => {
    if (!w.state) w.state = {};
    w.state.inheritanceResult = {
      items_meta: [{ id: 'cand_A', K: 3, seq_num: 1 }],
      band_index: [{ item_idx: 0, band: 0 }],
      cut: { group_id_per_band: new Int32Array([0]), n_groups: 1 },
      rtab: {},
    };
    const c = { id: 'cand_NOT_IN_RESULT', chrom: 'LG28', K: 3 };
    const s = w._inheritanceSuggestionsForCandidate(c);
    eq(Object.keys(s).length, 0);
  });

  // ============================================================
  // candidateHaplotypeAnnotationsHtml
  // ============================================================
  t('html: returns empty string for null candidate', () => {
    eq(w.candidateHaplotypeAnnotationsHtml(null), '');
  });

  t('html: returns empty string for K=0 candidate', () => {
    const c = { id: 'cand_0', chrom: 'LG28', K: 0 };
    eq(w.candidateHaplotypeAnnotationsHtml(c), '');
  });

  t('html: K=3 produces 3 input rows', () => {
    clearLS();
    const c = { id: 'cand_h', chrom: 'LG28', K: 3 };
    const html = w.candidateHaplotypeAnnotationsHtml(c);
    if (html.indexOf('hapLabelsSection') < 0) throw new Error('no section id');
    const inputs = (html.match(/data-band-idx=/g) || []).length;
    // Each band has BOTH a picker and a custom input → 2 elements per band
    // Plus the vocab picker doesn't carry data-band-idx
    if (inputs < 3) throw new Error(`expected at least 3, got ${inputs}`);
  });

  t('html: K=6 produces 6 input rows', () => {
    clearLS();
    const c = { id: 'cand_h6', chrom: 'LG28', K: 6 };
    const html = w.candidateHaplotypeAnnotationsHtml(c);
    const inputs = (html.match(/data-band-idx=/g) || []).length;
    if (inputs < 6) throw new Error(`expected at least 6, got ${inputs}`);
  });

  t('html: prefills user-typed values', () => {
    clearLS();
    const c = { id: 'cand_p', chrom: 'LG28', K: 3 };
    w.setHaplotypeLabel(c, 0, 'HOM_REF');
    w.setHaplotypeLabel(c, 1, 'HET');
    const html = w.candidateHaplotypeAnnotationsHtml(c);
    // Now rendered as <option ... selected>HOM_REF</option> in dropdowns
    if (html.indexOf('selected>HOM_REF') < 0) throw new Error('did not prefill band 0');
    if (html.indexOf('selected>HET') < 0) throw new Error('did not prefill band 1');
  });

  t('html: shows inheritance suggestions when result exists', () => {
    clearLS();
    if (!w.state) w.state = {};
    w.state.inheritanceResult = {
      items_meta: [{ id: 'cand_S', K: 3, seq_num: 1 }],
      band_index: [
        { item_idx: 0, band: 0 },
        { item_idx: 0, band: 1 },
        { item_idx: 0, band: 2 },
      ],
      cut: { group_id_per_band: new Int32Array([0, 1, 2]), n_groups: 3 },
      rtab: {},
    };
    const c = { id: 'cand_S', chrom: 'LG28', K: 3 };
    const html = w.candidateHaplotypeAnnotationsHtml(c);
    // Compact form: 'inh0' / 'inh1' / 'inh2'
    if (html.indexOf('inh0') < 0) throw new Error('missing suggestion 0');
    if (html.indexOf('inh1') < 0) throw new Error('missing suggestion 1');
    if (html.indexOf('inh2') < 0) throw new Error('missing suggestion 2');
  });

  t('html: escapes quotes in user input to prevent attribute breakage', () => {
    clearLS();
    const c = { id: 'cand_Q', chrom: 'LG28', K: 3 };
    w.setHaplotypeVocab(c, 'free');   // free vocab uses text input
    w.setHaplotypeLabel(c, 0, 'with "quotes"');
    const html = w.candidateHaplotypeAnnotationsHtml(c);
    if (html.indexOf('value="with "quotes""') >= 0) {
      throw new Error('raw quotes break attribute');
    }
    if (html.indexOf('with &quot;quotes&quot;') < 0) {
      throw new Error('quotes not escaped');
    }
  });

  // ============================================================
  // _wireCandidateHaplotypeAnnotations
  // ============================================================
  t('wire: picker change writes to candidate object and localStorage', () => {
    clearLS();
    const c = { id: 'cand_W', chrom: 'LG28', K: 3 };
    // Build the DOM for this candidate
    const container = w.document.createElement('div');
    container.innerHTML = w.candidateHaplotypeAnnotationsHtml(c);
    w.document.body.appendChild(container);
    w._wireCandidateHaplotypeAnnotations(c);

    // For 'standard' vocab (default for K=3), it's a select element
    const picker = w.document.querySelector('.hap-label-picker[data-band-idx="0"]');
    if (!picker) throw new Error('picker not found');
    picker.value = 'HOM_REF';
    picker.dispatchEvent(new w.Event('change', { bubbles: true }));

    eq(c.haplotype_labels[0], 'HOM_REF');
    const key = w._hapLabelLSKey(c.chrom, c.id);
    const raw = w.localStorage.getItem(key);
    if (!raw) throw new Error('not persisted');

    // Cleanup
    container.remove();
  });

  t('wire: free vocab uses text input', () => {
    clearLS();
    const c = { id: 'cand_FREE', chrom: 'LG28', K: 3 };
    w.setHaplotypeVocab(c, 'free');
    const container = w.document.createElement('div');
    container.innerHTML = w.candidateHaplotypeAnnotationsHtml(c);
    w.document.body.appendChild(container);
    w._wireCandidateHaplotypeAnnotations(c);

    const inp = w.document.querySelector('.hap-label-free[data-band-idx="0"]');
    if (!inp) throw new Error('free input not found');
    inp.value = 'custom_label';
    inp.dispatchEvent(new w.Event('change', { bubbles: true }));
    eq(c.haplotype_labels[0], 'custom_label');
    container.remove();
  });

  t('wire: handles null candidate without crashing', () => {
    w._wireCandidateHaplotypeAnnotations(null);
  });

  t('wire: handles missing DOM section without crashing', () => {
    const c = { id: 'cand_NoDOM', chrom: 'LG28', K: 3 };
    // Don't render the HTML — section won't exist
    w._wireCandidateHaplotypeAnnotations(c);
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
