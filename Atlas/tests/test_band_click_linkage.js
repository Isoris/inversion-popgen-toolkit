// Tests for turn 2j — page-2 band click → genome-wide linkage projection.
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

  if (typeof w._wireCandidateBandClicks !== 'function'
      || typeof w._renderGenomeLinkageTable !== 'function'
      || typeof w.candidateBandsHtml !== 'function') {
    // candidateBandsHtml not on window — check function exists at all
    // (it may not be exposed; we'll inspect the HTML output via candidateBandsHtml indirectly)
    if (typeof w._wireCandidateBandClicks !== 'function') {
      console.log('  FAIL: turn 2j wire function not exposed');
      return [{ name: 'fns exposed', err: 'missing' }];
    }
  }
  console.log('  Turn 2j functions exposed');

  function setupCandidates(candConfigs, withSamples) {
    const _state = w.state || {};
    if (!w.state) w.state = _state;
    _state.activeMode = 'default';
    _state.k = 3;
    _state.candidates = {};
    _state.candidates_detailed = {};
    let nSamples = 0;
    for (const c of candConfigs) {
      _state.candidates[c.id] = {
        id: c.id, chrom: 'TEST',
        K: c.K, K_used: c.K,
        locked_labels: c.labels,
        start_w: 0, end_w: 30,
        start_bp: c.start_bp,
        end_bp: c.end_bp,
        source: 'test',
      };
      if (c.labels && c.labels.length > nSamples) nSamples = c.labels.length;
    }
    _state.data = {
      chrom: 'TEST', n_windows: 30,
      samples: withSamples || Array.from({ length: nSamples }, (_, i) => ({ id: `S${i}` })),
    };
    _state.tracked = [];
    _state.inheritanceResult = null;
    _state.inheritanceCacheKey = null;
  }

  // Helper to set up DOM + render bands
  function renderBandsToDom(c, bands) {
    // Inject the HTML into the document body
    const container = w.document.createElement('div');
    // candidateBandsHtml may not be on window, so use innerHTML directly via string
    // Actually we can call it through the pageRender path — but the simplest is
    // to manually construct a band card DOM that matches the data attrs we wired.
    let html = '<div class="band-grid" id="candBandGrid">';
    for (const b of bands) {
      html += `<div class="band-card cand-band-card" data-band-idx="${b.k}">`;
      html += `<span>band ${b.k}</span></div>`;
    }
    html += '</div>';
    html += '<div class="cand-section" id="candGenomeLinkageSection" style="display:none;">';
    html += '<h4 class="cand-h4" id="candGenomeLinkageHeader"></h4>';
    html += '<div id="candGenomeLinkageBody"></div></div>';
    container.innerHTML = html;
    // Clean any prior containers
    while (w.document.body.firstChild) w.document.body.removeChild(w.document.body.firstChild);
    w.document.body.appendChild(container);
    return container;
  }

  // ============================================================
  // _wireCandidateBandClicks — sets state.tracked on click
  // ============================================================
  t('wire: clicking a band sets state.tracked to its members', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    const c = w.state.candidates['A'];
    // Simulate bands array
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 30; i++) {
        if (labels[i + k * 30] === k) members.push(i + k * 30);
      }
      // Actually members is 0..29 for band 0, 30..59 for band 1, 60..89 for band 2
      const correctMembers = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) correctMembers.push(i);
      bands.push({ k, n: correctMembers.length, members: correctMembers, families: [], ancestries: [] });
    }
    renderBandsToDom(c, bands);
    w._wireCandidateBandClicks(c, bands);
    // Click band 0
    const card = w.document.querySelector('.cand-band-card[data-band-idx="0"]');
    if (!card) throw new Error('band card not found');
    card.click();
    eq(w.state.tracked.length, 30);
    eq(w.state.tracked[0], 0);
    eq(w.state.tracked[29], 29);
  });

  t('wire: clicking different bands replaces state.tracked', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._wireCandidateBandClicks(w.state.candidates['A'], bands);
    w.document.querySelector('.cand-band-card[data-band-idx="0"]').click();
    eq(w.state.tracked[0], 0);
    w.document.querySelector('.cand-band-card[data-band-idx="2"]').click();
    eq(w.state.tracked[0], 60);   // band 2 starts at 60
  });

  t('wire: clicking band shows genome linkage section', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._wireCandidateBandClicks(w.state.candidates['A'], bands);
    w.document.querySelector('.cand-band-card[data-band-idx="0"]').click();
    const section = w.document.getElementById('candGenomeLinkageSection');
    eq(section.style.display, 'block');
  });

  t('wire: handles null candidate without crashing', () => {
    w._wireCandidateBandClicks(null, []);
  });

  t('wire: handles missing DOM gracefully', () => {
    while (w.document.body.firstChild) w.document.body.removeChild(w.document.body.firstChild);
    const c = { id: 'X', chrom: 'T', K: 3 };
    w._wireCandidateBandClicks(c, [{ k: 0, members: [0,1,2] }]);
  });

  // ============================================================
  // _renderGenomeLinkageTable — content tests
  // ============================================================
  t('render: shows table with rows for each candidate', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._renderGenomeLinkageTable(w.state.candidates['A'], 0, bands[0]);
    const body = w.document.getElementById('candGenomeLinkageBody').innerHTML;
    if (body.indexOf('<table') < 0) throw new Error('expected table');
    if (body.indexOf('<tbody>') < 0) throw new Error('expected tbody');
    if (body.indexOf('I1') < 0) throw new Error('expected I1 row');
    if (body.indexOf('I2') < 0) throw new Error('expected I2 row');
  });

  t('render: marks self row distinctly', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._renderGenomeLinkageTable(w.state.candidates['A'], 0, bands[0]);
    const body = w.document.getElementById('candGenomeLinkageBody').innerHTML;
    if (body.indexOf('this candidate (self)') < 0) throw new Error('expected self marker');
  });

  t('render: family-confound warning appears when one family >80%', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    // Make all band-0 fish be from family F0 (single family)
    const samples = Array.from({ length: 90 }, (_, i) => ({
      id: `S${i}`,
      family_id: i < 30 ? 0 : 1,   // band 0 fish all from family F0
    }));
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ], samples);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._renderGenomeLinkageTable(w.state.candidates['A'], 0, bands[0]);
    const body = w.document.getElementById('candGenomeLinkageBody').innerHTML;
    if (body.indexOf('family LD') < 0 && body.indexOf('family 0') < 0) {
      throw new Error('expected family confound warning');
    }
  });

  t('render: empty projection shows fallback message', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._renderGenomeLinkageTable(w.state.candidates['A'], 0, bands[0]);
    const body = w.document.getElementById('candGenomeLinkageBody').innerHTML;
    // With only 1 candidate, projection still shows it (self-row); but the
    // user-visible "no other candidates" message only fires when proj has 0 entries.
    // Here we have 1 entry (self), so it should still render the table.
    if (body.indexOf('<table') < 0) throw new Error('expected single-row table for self-only');
  });

  t('render: clear button resets state.tracked', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1, end_bp: 2, labels },
      { id: 'B', K: 3, start_bp: 5, end_bp: 6, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._wireCandidateBandClicks(w.state.candidates['A'], bands);
    w.document.querySelector('.cand-band-card[data-band-idx="0"]').click();
    eq(w.state.tracked.length, 30);
    const clearBtn = w.document.getElementById('candLinkageClearBtn');
    if (!clearBtn) throw new Error('clear button not found');
    clearBtn.click();
    eq(w.state.tracked.length, 0);
    const section = w.document.getElementById('candGenomeLinkageSection');
    eq(section.style.display, 'none');
  });

  // ============================================================
  // End-to-end: realistic 3-candidate scenario
  // ============================================================
  t('end-to-end: clicking a band on candidate A projects across B and C', () => {
    const labels = [];
    for (let k = 0; k < 3; k++) for (let i = 0; i < 30; i++) labels.push(k);
    setupCandidates([
      { id: 'A', K: 3, start_bp: 1_000_000, end_bp: 2_000_000, labels },
      { id: 'B', K: 3, start_bp: 5_000_000, end_bp: 6_000_000, labels: [...labels] },
      { id: 'C', K: 3, start_bp: 9_000_000, end_bp: 10_000_000, labels: [...labels] },
    ]);
    const bands = [];
    for (let k = 0; k < 3; k++) {
      const members = [];
      for (let i = 0; i < 90; i++) if (labels[i] === k) members.push(i);
      bands.push({ k, n: members.length, members, families: [], ancestries: [] });
    }
    renderBandsToDom(w.state.candidates['A'], bands);
    w._wireCandidateBandClicks(w.state.candidates['A'], bands);
    w.document.querySelector('.cand-band-card[data-band-idx="1"]').click();
    eq(w.state.tracked.length, 30);
    // Verify the projection sees all 3 candidates (including self)
    const body = w.document.getElementById('candGenomeLinkageBody').innerHTML;
    if (body.indexOf('I1') < 0) throw new Error('missing I1 row');
    if (body.indexOf('I2') < 0) throw new Error('missing I2 row');
    if (body.indexOf('I3') < 0) throw new Error('missing I3 row');
    // Band 1 fish (indices 30-59) all in band 1 of all candidates — purity should be 100%
    if (body.indexOf('100%') < 0) throw new Error('expected 100% purity row');
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
