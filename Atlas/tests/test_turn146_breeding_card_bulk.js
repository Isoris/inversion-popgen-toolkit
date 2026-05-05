// =============================================================================
// turn 146 — Breeding-readiness card, Turn D (bulk catalogue export)
// =============================================================================
// SPEC: specs_new_turn131/SPEC_per_candidate_breeding_readiness_card.md §4
// Slice 3 (bulk PDF export from the catalogue page).
//
// Builds on Turn C's renderer (turn 144). Adds:
//   - _filterCandsForBreedingExport(cands, tierMode)
//   - _buildBreedingCardsCombinedHTML(cards, opts)
//   - _buildBreedingCardsJSONBundle(cards, opts)
//   - _exportBreedingCardsHTML(opts)             — DOM dispatcher
//   - _exportBreedingCardsJSON(opts)             — DOM dispatcher
//   - _wireCatalogueBreedingExportBtns()         — toolbar wiring
//   - Catalogue-page UI: 🧬 breeding cards group with tier dropdown +
//     HTML / JSON buttons
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

// Build a sandbox containing all the helpers defined in turn 142 + 143 +
// 144 + 146. This regex span matches from the Turn-A loader banner through
// the end of Turn D's last window export.
function buildSandbox(extra) {
  const re = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 146 — Breeding-readiness card, Turn D \(bulk catalogue export\)[\s\S]*?window\._wireCatalogueBreedingExportBtns = _wireCatalogueBreedingExportBtns;\s*\n[\s\S]*?\n\}/;
  const m = html.match(re);
  if (!m) throw new Error('combined Turn-B/C/D region not found');
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  if (!normal) throw new Error('normalCDF region not found');
  const _ls = new Map();
  const sandbox = Object.assign({
    state: { data: null, cohortDiversity: null, candidateList: [] },
    window: {},
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
    console,
  }, extra || {});
  const ctx = vm.createContext(sandbox);
  vm.runInContext(normal[0] + '\n' + m[0], ctx);
  ctx.__ls = _ls;
  return ctx;
}

// Build a synthetic candidate that matches the shape Turn B's
// _breedingCardKaryotypePerSample expects.
function makeCand(opts) {
  const o = opts || {};
  const labels = [];
  const refF = o.refF || [0.20, 0.21, 0.19];
  const hetF = o.hetF || [0.24, 0.26, 0.25];
  const invF = o.invF || [0.31, 0.33, 0.32];
  for (const _ of refF) labels.push(0);
  for (const _ of hetF) labels.push(1);
  for (const _ of invF) labels.push(2);
  return {
    id:        o.id || 'TEST',
    chrom:     o.chrom || 'LG28',
    start_bp:  o.start_bp || 14_000_000,
    end_bp:    o.end_bp   || 16_500_000,
    K: 3,
    tier:        o.tier || 'Tier 1',
    confidence:  o.confidence || 'high',
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
}

// ============================================================================
// 1. Source-level definitions
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');
ok('_filterCandsForBreedingExport defined',
   /function _filterCandsForBreedingExport\(/.test(html));
ok('_buildBreedingCardsCombinedHTML defined',
   /function _buildBreedingCardsCombinedHTML\(/.test(html));
ok('_buildBreedingCardsJSONBundle defined',
   /function _buildBreedingCardsJSONBundle\(/.test(html));
ok('_exportBreedingCardsHTML defined',
   /function _exportBreedingCardsHTML\(/.test(html));
ok('_exportBreedingCardsJSON defined',
   /function _exportBreedingCardsJSON\(/.test(html));
ok('_wireCatalogueBreedingExportBtns defined',
   /function _wireCatalogueBreedingExportBtns\(/.test(html));
ok('_BREEDING_EXPORT_TIER_MODES const defined',
   /const _BREEDING_EXPORT_TIER_MODES\s*=\s*\{/.test(html));
ok('Turn D banner present',
   /turn 146 — Breeding-readiness card, Turn D \(bulk catalogue export\)/.test(html));

// ============================================================================
// 2. Window exports
// ============================================================================
console.log('\n=== 2. Window exports ===');
const wins = [
  '_filterCandsForBreedingExport',
  '_BREEDING_EXPORT_TIER_MODES',
  '_buildBreedingCardsCombinedHTML',
  '_buildBreedingCardsJSONBundle',
  '_exportBreedingCardsHTML',
  '_exportBreedingCardsJSON',
  '_breedingExportTrigger',
  '_gatherBreedingCardsFromState',
  '_wireCatalogueBreedingExportBtns',
];
for (const w of wins) {
  const re = new RegExp('window\\.' + w + '\\s*=\\s*' + w);
  ok('window.' + w, re.test(html));
}

// ============================================================================
// 3. Catalogue-page UI shell
// ============================================================================
console.log('\n=== 3. Catalogue-page UI shell ===');
ok('🧬 breeding cards label present',
   /<span class="export-group-label">🧬 breeding cards<\/span>/.test(html));
ok('tier select #catBreedingTierSel present',
   /id="catBreedingTierSel"/.test(html));
ok('tier select has 4 options',
   (() => {
     // Anchor to the specific dropdown — `value="all"` and others appear
     // elsewhere in the atlas.
     const m = html.match(/id="catBreedingTierSel"[\s\S]*?<\/select>/);
     if (!m) return false;
     return (m[0].match(/<option value="(tier_1_2|tier_1|tier_1_2_3|all)"/g) || []).length === 4;
   })());
ok('Tier 1+2 is the default option',
   /<option value="tier_1_2" selected>Tier 1\+2<\/option>/.test(html));
ok('HTML export button #catExportBreedingHTML present',
   /id="catExportBreedingHTML"/.test(html));
ok('JSON export button #catExportBreedingJSON present',
   /id="catExportBreedingJSON"/.test(html));
ok('export-group select CSS rule added',
   /\.export-group select\.export-fmt/.test(html));
ok('select sits inside an export-group span',
   /<span class="export-group"[\s\S]{0,600}?id="catBreedingTierSel"/.test(html));

// ============================================================================
// 4. _filterCandsForBreedingExport — pure tier filter
// ============================================================================
console.log('\n=== 4. _filterCandsForBreedingExport ===');
{
  const ctx = buildSandbox();
  const filt = ctx.window._filterCandsForBreedingExport;

  // 4.1 Default mode = tier_1_2
  {
    const cands = [
      { id: 'a', tier: 'Tier 1' },
      { id: 'b', tier: 'Tier 2' },
      { id: 'c', tier: 'Tier 3' },
      { id: 'd', tier: 'Tier 4' },
      { id: 'e' },                  // no tier
      null,
    ];
    const r = filt(cands);
    ok('default mode kept Tier 1+2 only',
       r.kept.length === 2 && r.kept[0].id === 'a' && r.kept[1].id === 'b');
    ok('default mode dropped Tier 3/4 (off-tier)',
       r.dropped_off_tier === 2);
    ok('default mode dropped untiered (no_tier)',
       r.dropped_no_tier === 1);
    ok('default mode label = tier_1_2',
       r.mode === 'tier_1_2');
    ok('default accepted_tiers = ["Tier 1","Tier 2"]',
       JSON.stringify(r.accepted_tiers) === '["Tier 1","Tier 2"]');
  }

  // 4.2 tier_1 mode
  {
    const cands = [
      { id: 'a', tier: 'Tier 1' },
      { id: 'b', tier: 'Tier 2' },
      { id: 'c', tier: 'Tier 1' },
    ];
    const r = filt(cands, 'tier_1');
    ok('tier_1 mode kept only Tier 1',
       r.kept.length === 2 && r.kept.every(c => c.tier === 'Tier 1'));
    ok('tier_1 mode dropped Tier 2 (off-tier)',
       r.dropped_off_tier === 1);
  }

  // 4.3 tier_1_2_3 mode
  {
    const cands = [
      { id: 'a', tier: 'Tier 1' },
      { id: 'b', tier: 'Tier 2' },
      { id: 'c', tier: 'Tier 3' },
      { id: 'd', tier: 'Tier 4' },
    ];
    const r = filt(cands, 'tier_1_2_3');
    ok('tier_1_2_3 mode kept Tier 1+2+3', r.kept.length === 3);
    ok('tier_1_2_3 mode dropped only Tier 4',
       r.dropped_off_tier === 1 && r.dropped_no_tier === 0);
  }

  // 4.4 all mode includes everything
  {
    const cands = [
      { id: 'a', tier: 'Tier 1' },
      { id: 'b', tier: 'Tier 4' },
      { id: 'c' },                  // no tier
    ];
    const r = filt(cands, 'all');
    ok('all mode kept everything', r.kept.length === 3);
    ok('all mode dropped 0',
       r.dropped_off_tier === 0 && r.dropped_no_tier === 0);
    ok('all mode accepted_tiers = null', r.accepted_tiers === null);
  }

  // 4.5 unknown tierMode falls back to default
  {
    const r = filt([{ id: 'a', tier: 'Tier 1' }, { id: 'b', tier: 'Tier 3' }],
                   'totally-bogus-mode');
    ok('unknown mode falls back to tier_1_2', r.mode === 'tier_1_2');
    ok('unknown mode kept Tier 1, dropped Tier 3',
       r.kept.length === 1 && r.dropped_off_tier === 1);
  }

  // 4.6 Custom tier array
  {
    const r = filt(
      [{ id: 'a', tier: 'Tier 1' },
       { id: 'b', tier: 'Tier 4' },
       { id: 'c', tier: 'Tier 2' }],
      ['Tier 1', 'Tier 4']);
    ok('custom array tier filter kept Tier 1+4',
       r.kept.length === 2 && r.kept.map(c => c.tier).sort().join(',') === 'Tier 1,Tier 4');
    ok('custom array filter mode = "custom"', r.mode === 'custom');
  }

  // 4.7 Defensive: non-array input
  {
    const r = filt(null);
    ok('null input returns empty kept', Array.isArray(r.kept) && r.kept.length === 0);
  }
}

// ============================================================================
// 5. _buildBreedingCardsJSONBundle — pure envelope builder
// ============================================================================
console.log('\n=== 5. _buildBreedingCardsJSONBundle ===');
{
  const ctx = buildSandbox();
  const fxA = makeCand({ id: 'A', tier: 'Tier 1' });
  const fxB = makeCand({ id: 'B', tier: 'Tier 2' });
  ctx.state.data = { chrom: 'LG28', samples: [
    { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
    { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
    { cga: 'CGA007' }, { cga: 'CGA008' }, { cga: 'CGA009' },
  ] };
  ctx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', k8: 'K1', f_roh: 0.20 },
      { sample_id: 'CGA002', k8: 'K1', f_roh: 0.21 },
      { sample_id: 'CGA003', k8: 'K1', f_roh: 0.19 },
      { sample_id: 'CGA004', k8: 'K1', f_roh: 0.24 },
      { sample_id: 'CGA005', k8: 'K1', f_roh: 0.26 },
      { sample_id: 'CGA006', k8: 'K1', f_roh: 0.25 },
      { sample_id: 'CGA007', k8: 'K1', f_roh: 0.31 },
      { sample_id: 'CGA008', k8: 'K1', f_roh: 0.33 },
      { sample_id: 'CGA009', k8: 'K1', f_roh: 0.32 },
    ],
  });
  const cardA = ctx.window._buildBreedingCard(fxA);
  const cardB = ctx.window._buildBreedingCard(fxB);
  const bundle = ctx.window._buildBreedingCardsJSONBundle([cardA, cardB], {
    chrom: 'LG28',
    filter_meta: { mode: 'tier_1_2', accepted_tiers: ['Tier 1', 'Tier 2'],
                   dropped_no_tier: 0, dropped_off_tier: 0, n_total: 2 },
  });
  ok('schema = breeding_readiness_card_bundle_v1',
     bundle.schema === 'breeding_readiness_card_bundle_v1');
  ok('generated_at is ISO',                      /^\d{4}-\d{2}-\d{2}T/.test(bundle.generated_at));
  ok('chrom = LG28',                             bundle.chrom === 'LG28');
  ok('n_cards = 2',                              bundle.n_cards === 2);
  ok('filter.mode = tier_1_2',                   bundle.filter.mode === 'tier_1_2');
  ok('cards array has 2 entries',                Array.isArray(bundle.cards) && bundle.cards.length === 2);
  ok('cards[0] is the EXACT _buildBreedingCard output (deep equal)',
     JSON.stringify(bundle.cards[0]) === JSON.stringify(cardA));
  ok('cards[1] preserved order',                 bundle.cards[1].candidate.id === 'B');
  // JSON-serializable as a whole
  ok('bundle JSON.stringify roundtrips',
     (() => { try { JSON.parse(JSON.stringify(bundle)); return true; } catch (_) { return false; } })());

  // Empty input
  const empty = ctx.window._buildBreedingCardsJSONBundle([], {});
  ok('empty bundle: n_cards=0',                  empty.n_cards === 0);
  ok('empty bundle: cards is []',                Array.isArray(empty.cards) && empty.cards.length === 0);
  ok('empty bundle: chrom is null',              empty.chrom === null);

  // Non-array input
  const nonArr = ctx.window._buildBreedingCardsJSONBundle(null, {});
  ok('null cards → n_cards = 0',                 nonArr.n_cards === 0);
}

// ============================================================================
// 6. _buildBreedingCardsCombinedHTML — pure combined-doc builder
// ============================================================================
console.log('\n=== 6. _buildBreedingCardsCombinedHTML ===');
{
  const ctx = buildSandbox();
  const fxA = makeCand({ id: 'CAND_A', tier: 'Tier 1' });
  const fxB = makeCand({ id: 'CAND_B', tier: 'Tier 2' });
  ctx.state.data = { chrom: 'LG28', samples: [
    { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
    { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
    { cga: 'CGA007' }, { cga: 'CGA008' }, { cga: 'CGA009' },
  ] };
  ctx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', k8: 'K1', f_roh: 0.20 },
      { sample_id: 'CGA002', k8: 'K1', f_roh: 0.21 },
      { sample_id: 'CGA003', k8: 'K1', f_roh: 0.19 },
      { sample_id: 'CGA004', k8: 'K1', f_roh: 0.24 },
      { sample_id: 'CGA005', k8: 'K1', f_roh: 0.26 },
      { sample_id: 'CGA006', k8: 'K1', f_roh: 0.25 },
      { sample_id: 'CGA007', k8: 'K1', f_roh: 0.31 },
      { sample_id: 'CGA008', k8: 'K1', f_roh: 0.33 },
      { sample_id: 'CGA009', k8: 'K1', f_roh: 0.32 },
    ],
  });
  const cardA = ctx.window._buildBreedingCard(fxA);
  const cardB = ctx.window._buildBreedingCard(fxB);
  const out = ctx.window._buildBreedingCardsCombinedHTML([cardA, cardB], {
    chrom: 'LG28',
    filter_meta: { mode: 'tier_1_2', accepted_tiers: ['Tier 1', 'Tier 2'],
                   dropped_no_tier: 0, dropped_off_tier: 0, n_total: 2 },
  });
  ok('returns string',                            typeof out === 'string' && out.length > 1000);
  ok('starts with !doctype',                      /^<!doctype html>/i.test(out));
  ok('ends with </html>',                         /<\/html>\s*$/.test(out.trim()));
  ok('has charset utf-8',                         /<meta charset="utf-8"/.test(out));
  ok('title contains LG28',                       /<title>[^<]*LG28/.test(out));
  ok('title shows n=2',                           /<title>[^<]*n=2/.test(out));
  ok('embeds @page A4',                           /@page\s*\{[^}]*size:\s*A4/.test(out));
  ok('embeds page-break-after CSS for cards',     /\.brc-print-page[\s\S]{0,200}?page-break-after:\s*always/.test(out));
  ok('TOC table present',                         /class="brc-bundle-toc"/.test(out));
  ok('TOC has both candidate IDs',
     out.indexOf('CAND_A') >= 0 && out.indexOf('CAND_B') >= 0);
  ok('TOC has filter mode shown',                 /<b>Filter:<\/b>\s*tier_1_2/.test(out));
  ok('TOC has tiers included list',               /<b>Tiers included:<\/b>\s*Tier 1,\s*Tier 2/.test(out));
  ok('TOC has total count',                       /<b>Total candidates:<\/b>\s*2/.test(out));
  ok('TOC has bundle count',                      /<b>In bundle:<\/b>\s*2/.test(out));
  ok('TOC has generated timestamp',               /<b>Generated:<\/b>/.test(out));
  ok('exactly two .brc-print-page sections',
     (out.match(/class="brc-print-page"/g) || []).length === 2);
  ok('first section has anchor brc-card-1',       /id="brc-card-1"/.test(out));
  ok('second section has anchor brc-card-2',      /id="brc-card-2"/.test(out));
  ok('TOC links use the anchor href',
     /href="#brc-card-1"/.test(out) && /href="#brc-card-2"/.test(out));
  ok('contains an actual rendered .brc-card from cardA',
     out.indexOf('CAND_A') >= 0);
  ok('hint banner present by default',            /class="brc-print-hint/.test(out));
  ok('auto-print OFF by default (no script tag)',
     out.indexOf('window.print()') === -1);

  // Empty bundle
  const emptyOut = ctx.window._buildBreedingCardsCombinedHTML([], {
    chrom: 'LG28',
    filter_meta: { mode: 'tier_1_2', accepted_tiers: ['Tier 1', 'Tier 2'],
                   dropped_no_tier: 0, dropped_off_tier: 5, n_total: 5 },
  });
  ok('empty bundle: still returns full HTML doc', /^<!doctype html>/i.test(emptyOut));
  ok('empty bundle: warns no candidates matched', /No candidates matched/.test(emptyOut));
  ok('empty bundle: 0 .brc-print-page sections',
     (emptyOut.match(/class="brc-print-page"/g) || []).length === 0);

  // auto_print=true → window.print() in script
  const autoOut = ctx.window._buildBreedingCardsCombinedHTML([cardA], {
    chrom: 'LG28', auto_print: true,
  });
  ok('auto_print=true emits window.print()',      autoOut.indexOf('window.print()') >= 0);
  ok('auto_print=true keeps hint by default',     /class="brc-print-hint/.test(autoOut));

  // show_hint=false drops the banner DIV (the CSS class name still
  // appears in the embedded stylesheet, which is fine — what matters
  // is that the user-visible <div> isn't rendered).
  const noHintOut = ctx.window._buildBreedingCardsCombinedHTML([cardA], {
    chrom: 'LG28', show_hint: false,
  });
  ok('show_hint=false drops the hint',            noHintOut.indexOf('<div class="brc-print-hint') === -1);
}

// ============================================================================
// 7. _gatherBreedingCardsFromState — DOM-light orchestrator
// ============================================================================
console.log('\n=== 7. _gatherBreedingCardsFromState ===');
{
  // Sandbox with state pre-populated. Mock alert so we can capture it.
  const _alerts = [];
  const ctx = buildSandbox({ alert: (msg) => _alerts.push(msg) });

  // 7.1 Empty candidateList
  ctx.state.data = { chrom: 'LG28', samples: [] };
  ctx.state.candidateList = [];
  let r = ctx.window._gatherBreedingCardsFromState();
  ok('empty candidateList → null + alert',         r === null && _alerts.length >= 1);
  ok('empty alert mentions promote',               /Promote at least one/.test(_alerts[0]));

  // 7.2 Off-tier filter result
  _alerts.length = 0;
  ctx.state.candidateList = [
    makeCand({ id: 'X', tier: 'Tier 4' }),
    makeCand({ id: 'Y', tier: 'Tier 3' }),
  ];
  r = ctx.window._gatherBreedingCardsFromState({ tierMode: 'tier_1_2' });
  ok('off-tier filter → null + alert',             r === null && _alerts.length >= 1);
  ok('off-tier alert mentions tier dropdown',
     /tier filter|tier dropdown/.test(_alerts[0]));

  // 7.3 Happy path — 2 in-tier candidates
  ctx.state.candidateList = [
    makeCand({ id: 'A', tier: 'Tier 1' }),
    makeCand({ id: 'B', tier: 'Tier 2' }),
    makeCand({ id: 'C', tier: 'Tier 3' }),
  ];
  ctx.state.data = { chrom: 'LG28', samples: [
    { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
    { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
    { cga: 'CGA007' }, { cga: 'CGA008' }, { cga: 'CGA009' },
  ] };
  ctx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', k8: 'K1', f_roh: 0.20 },
      { sample_id: 'CGA002', k8: 'K1', f_roh: 0.21 },
      { sample_id: 'CGA003', k8: 'K1', f_roh: 0.19 },
      { sample_id: 'CGA004', k8: 'K1', f_roh: 0.24 },
      { sample_id: 'CGA005', k8: 'K1', f_roh: 0.26 },
      { sample_id: 'CGA006', k8: 'K1', f_roh: 0.25 },
      { sample_id: 'CGA007', k8: 'K1', f_roh: 0.31 },
      { sample_id: 'CGA008', k8: 'K1', f_roh: 0.33 },
      { sample_id: 'CGA009', k8: 'K1', f_roh: 0.32 },
    ],
  });
  r = ctx.window._gatherBreedingCardsFromState({ tierMode: 'tier_1_2' });
  ok('happy path returns object',                 r !== null && typeof r === 'object');
  ok('cards.length = 2 (Tier 1+2 only)',          r.cards.length === 2);
  ok('chrom = LG28',                              r.chrom === 'LG28');
  ok('n_input = 3 (all candidates)',              r.n_input === 3);
  ok('filter_meta.dropped_off_tier = 1 (Tier 3)', r.filter_meta.dropped_off_tier === 1);
  ok('filter_meta.mode = tier_1_2',               r.filter_meta.mode === 'tier_1_2');

  // 7.4 'all' tier mode includes Tier 3 too
  r = ctx.window._gatherBreedingCardsFromState({ tierMode: 'all' });
  ok('all-tier mode: 3 cards',                    r.cards.length === 3);
  ok('all-tier mode: 0 dropped',
     r.filter_meta.dropped_off_tier === 0 && r.filter_meta.dropped_no_tier === 0);
}

// ============================================================================
// 8. End-to-end DOM dispatcher with mocked Blob + URL
// ============================================================================
console.log('\n=== 8. _exportBreedingCardsHTML / JSON dispatcher ===');
{
  const _downloads = [];
  const _alerts = [];
  // Mock the bits the dispatcher touches: Blob, URL, document.createElement,
  // appendChild, click. We also want to read what the Blob contained.
  const sandbox = {
    state: { data: null, cohortDiversity: null, candidateList: [] },
    window: {},
    document: {
      readyState: 'complete',
      addEventListener: () => {},
      body: {
        appendChild: () => {},
        removeChild: () => {},
      },
      createElement: (tag) => {
        const el = {
          tagName: tag.toUpperCase(),
          style: {}, dataset: {},
          set href(v) { this._href = v; },
          get href() { return this._href; },
          download: '',
          click: function () {
            _downloads.push({ href: this._href, download: this.download,
                              capturedBody: this.__capturedBody });
          },
          appendChild: () => {},
        };
        return el;
      },
    },
    Blob: class {
      constructor(parts, opts) {
        this.parts = parts;
        this.type  = (opts && opts.type) || '';
        this.size  = parts && parts[0] ? parts[0].length : 0;
        // capture body for the next createElement('a').click() call
        this.__body = parts && parts[0] || '';
      }
    },
    URL: {
      _last: null,
      createObjectURL: function (b) {
        // We carry the Blob body via the URL→createElement('a').__capturedBody hop.
        this._last = b;
        return 'blob:' + Math.random();
      },
      revokeObjectURL: () => {},
    },
    setTimeout: setTimeout,
    alert: (msg) => _alerts.push(msg),
    console,
    localStorage: { getItem: () => null, setItem: () => {}, removeItem: () => {} },
  };
  // Monkey-patch createElement to wire the captured body onto the next <a>
  const origCreate = sandbox.document.createElement;
  sandbox.document.createElement = (tag) => {
    const el = origCreate(tag);
    if (tag === 'a') {
      // Read the most recently created Blob body off URL._last
      Object.defineProperty(el, '__capturedBody', {
        get: () => (sandbox.URL._last && sandbox.URL._last.__body) || '',
      });
    }
    return el;
  };

  const ctx = vm.createContext(sandbox);
  const reBoth = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 146 — Breeding-readiness card, Turn D \(bulk catalogue export\)[\s\S]*?window\._wireCatalogueBreedingExportBtns = _wireCatalogueBreedingExportBtns;\s*\n[\s\S]*?\n\}/;
  const m = html.match(reBoth);
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  vm.runInContext(normal[0] + '\n' + m[0], ctx);

  // Set up state: 3 candidates, Tier 1, Tier 2, Tier 3
  ctx.state.data = { chrom: 'LG28', samples: [
    { cga: 'CGA001' }, { cga: 'CGA002' }, { cga: 'CGA003' },
    { cga: 'CGA004' }, { cga: 'CGA005' }, { cga: 'CGA006' },
    { cga: 'CGA007' }, { cga: 'CGA008' }, { cga: 'CGA009' },
  ] };
  ctx.state.candidateList = [
    makeCand({ id: 'A', tier: 'Tier 1' }),
    makeCand({ id: 'B', tier: 'Tier 2' }),
    makeCand({ id: 'C', tier: 'Tier 3' }),
  ];
  ctx.window._storeCohortDiversity({
    tool: 'cohort_diversity_v1', schema_version: 1,
    samples: [
      { sample_id: 'CGA001', k8: 'K1', f_roh: 0.20 },
      { sample_id: 'CGA002', k8: 'K1', f_roh: 0.21 },
      { sample_id: 'CGA003', k8: 'K1', f_roh: 0.19 },
      { sample_id: 'CGA004', k8: 'K1', f_roh: 0.24 },
      { sample_id: 'CGA005', k8: 'K1', f_roh: 0.26 },
      { sample_id: 'CGA006', k8: 'K1', f_roh: 0.25 },
      { sample_id: 'CGA007', k8: 'K1', f_roh: 0.31 },
      { sample_id: 'CGA008', k8: 'K1', f_roh: 0.33 },
      { sample_id: 'CGA009', k8: 'K1', f_roh: 0.32 },
    ],
  });

  // 8.1 HTML export, default tier_1_2 → 2 candidates
  _downloads.length = 0;
  const ok1 = ctx.window._exportBreedingCardsHTML({ tierMode: 'tier_1_2' });
  ok('HTML export returns true',                  ok1 === true);
  ok('HTML export triggered exactly 1 download',  _downloads.length === 1);
  ok('HTML download filename ends in .html',
     _downloads[0] && /_breeding_cards_tier_1_2\.html$/.test(_downloads[0].download));
  ok('HTML download filename starts with chrom',
     _downloads[0] && _downloads[0].download.indexOf('LG28_') === 0);
  ok('HTML body starts with !doctype',
     _downloads[0] && /^<!doctype html>/i.test(_downloads[0].capturedBody));
  ok('HTML body contains both candidate IDs',
     _downloads[0] &&
     _downloads[0].capturedBody.indexOf('A') >= 0 &&
     _downloads[0].capturedBody.indexOf('B') >= 0);
  ok('HTML body does NOT contain Tier 3 candidate C',
     _downloads[0] &&
     /id="brc-card-1"[\s\S]*?CAND.*A|>A</.test(_downloads[0].capturedBody) &&
     // C is the off-tier id; its presence elsewhere (e.g. in totals) is OK,
     // but the rendered card shouldn't be there.
     (_downloads[0].capturedBody.match(/id="brc-card-/g) || []).length === 2);

  // 8.2 JSON export, all tiers → 3 candidates
  _downloads.length = 0;
  const ok2 = ctx.window._exportBreedingCardsJSON({ tierMode: 'all' });
  ok('JSON export returns true',                  ok2 === true);
  ok('JSON export triggered 1 download',          _downloads.length === 1);
  ok('JSON filename ends in .json',
     _downloads[0] && /_breeding_cards_all\.json$/.test(_downloads[0].download));
  // Parse the captured JSON body
  let parsed = null;
  try { parsed = JSON.parse(_downloads[0].capturedBody); } catch (_) {}
  ok('JSON body parses',                          parsed !== null);
  ok('JSON has schema = breeding_readiness_card_bundle_v1',
     parsed && parsed.schema === 'breeding_readiness_card_bundle_v1');
  ok('JSON has n_cards = 3',                      parsed && parsed.n_cards === 3);
  ok('JSON has filter.mode = all',                parsed && parsed.filter && parsed.filter.mode === 'all');
  ok('JSON cards in original order',
     parsed && parsed.cards.map(c => c.candidate.id).join(',') === 'A,B,C');

  // 8.3 Empty path: tier_1 with no Tier 1 candidates → null + alert
  _alerts.length = 0;
  ctx.state.candidateList = [
    makeCand({ id: 'B', tier: 'Tier 2' }),
    makeCand({ id: 'C', tier: 'Tier 3' }),
  ];
  _downloads.length = 0;
  const ok3 = ctx.window._exportBreedingCardsHTML({ tierMode: 'tier_1' });
  ok('empty result: HTML export returns false',   ok3 === false);
  ok('empty result: no download triggered',       _downloads.length === 0);
  ok('empty result: alert was raised',            _alerts.length >= 1);

  // 8.4 Empty candidateList → null + alert
  _alerts.length = 0;
  _downloads.length = 0;
  ctx.state.candidateList = [];
  const ok4 = ctx.window._exportBreedingCardsJSON({ tierMode: 'tier_1_2' });
  ok('no candidates: JSON export returns false',  ok4 === false);
  ok('no candidates: no download triggered',      _downloads.length === 0);
  ok('no candidates: alert mentions promote',
     _alerts.length >= 1 && /promote/i.test(_alerts[0]));
}

// ============================================================================
// 9. _wireCatalogueBreedingExportBtns — DOM idempotence
// ============================================================================
console.log('\n=== 9. _wireCatalogueBreedingExportBtns ===');
{
  // Tiny DOM mock — only need getElementById, addEventListener, dataset.
  const elements = new Map();
  function makeEl(id) {
    const el = {
      id, dataset: {}, value: 'tier_1_2',
      _listeners: {},
      addEventListener(type, fn) {
        if (!this._listeners[type]) this._listeners[type] = [];
        this._listeners[type].push(fn);
      },
      _fire(type, evt) {
        (this._listeners[type] || []).forEach(fn => fn(evt || {}));
      },
    };
    elements.set(id, el);
    return el;
  }
  const sel  = makeEl('catBreedingTierSel');
  const btnH = makeEl('catExportBreedingHTML');
  const btnJ = makeEl('catExportBreedingJSON');

  const _ls = new Map();
  const sandbox = {
    state: { data: null, candidateList: [] },
    window: {},
    document: {
      readyState: 'complete',
      addEventListener: () => {},
      getElementById: (id) => elements.get(id) || null,
      body: { appendChild: () => {}, removeChild: () => {} },
      createElement: () => ({ style: {}, click: () => {}, dataset: {} }),
    },
    Blob: class { constructor() {} },
    URL: { createObjectURL: () => 'blob:x', revokeObjectURL: () => {} },
    setTimeout: setTimeout,
    alert: () => {},
    console,
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
  };
  const ctx = vm.createContext(sandbox);
  const reBoth = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 146 — Breeding-readiness card, Turn D \(bulk catalogue export\)[\s\S]*?window\._wireCatalogueBreedingExportBtns = _wireCatalogueBreedingExportBtns;\s*\n[\s\S]*?\n\}/;
  const m = html.match(reBoth);
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  vm.runInContext(normal[0] + '\n' + m[0], ctx);

  // Pre-seed localStorage with a saved tier mode
  _ls.set('pca_scrubber_v3.breeding_export_tier', 'tier_1');
  ctx.window._wireCatalogueBreedingExportBtns();
  ok('saved tier value restored to dropdown',     sel.value === 'tier_1');
  ok('select marked _wired=1',                    sel.dataset._wired === '1');
  ok('HTML btn marked _wired=1',                  btnH.dataset._wired === '1');
  ok('JSON btn marked _wired=1',                  btnJ.dataset._wired === '1');
  ok('select has at least 1 change listener',     (sel._listeners.change || []).length === 1);
  ok('HTML btn has 1 click listener',             (btnH._listeners.click || []).length === 1);
  ok('JSON btn has 1 click listener',             (btnJ._listeners.click || []).length === 1);

  // Idempotent — second call doesn't add duplicate listeners
  ctx.window._wireCatalogueBreedingExportBtns();
  ok('idempotent: HTML btn still 1 listener',     (btnH._listeners.click || []).length === 1);
  ok('idempotent: JSON btn still 1 listener',     (btnJ._listeners.click || []).length === 1);
  ok('idempotent: select still 1 listener',       (sel._listeners.change || []).length === 1);

  // Changing the dropdown persists the value
  sel.value = 'tier_1_2_3';
  sel._fire('change', { target: sel });
  ok('selecting persists to localStorage',
     _ls.get('pca_scrubber_v3.breeding_export_tier') === 'tier_1_2_3');

  // Bogus saved value falls back to default (no restore)
  _ls.set('pca_scrubber_v3.breeding_export_tier', 'totally-bogus-mode');
  const sel2 = makeEl('catBreedingTierSel'); sel2.value = 'tier_1_2';
  const btnH2 = makeEl('catExportBreedingHTML');
  const btnJ2 = makeEl('catExportBreedingJSON');
  // Re-call won't re-restore (it's idempotent on dataset._wired) — so we
  // instead re-run on fresh elements
  ctx.window._wireCatalogueBreedingExportBtns();
  ok('bogus saved value: defaults preserved',     sel2.value === 'tier_1_2');
}

// ============================================================================
// 10. Real-fixture render — 226-sample LG28
// ============================================================================
console.log('\n=== 10. Real-fixture render ===');
{
  const FIXTURE_PATH = path.resolve(__dirname, '..', 'json', 'cohort_diversity_v1.json');
  let fixture = null;
  try { fixture = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf8')); } catch (_) {}
  if (fixture) {
    const ctx = buildSandbox();
    const chromSamples = fixture.samples.map(s => ({ cga: s.sample_id }));
    ctx.state.data = { chrom: 'LG28', samples: chromSamples };
    ctx.window._storeCohortDiversity(fixture);
    const labels = new Int8Array(226);
    for (let i = 0; i < 60; i++)   labels[i] = 0;
    for (let i = 60; i < 166; i++) labels[i] = 1;
    for (let i = 166; i < 226; i++) labels[i] = 2;
    const cand1 = {
      id: 'LG28_inv1', chrom: 'LG28',
      start_bp: 13_800_000, end_bp: 16_700_000, K: 3,
      tier: 'Tier 1', confidence: 'high',
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
    const cand2 = Object.assign({}, cand1, {
      id: 'LG28_inv2', start_bp: 28_000_000, end_bp: 30_500_000,
      tier: 'Tier 2',
    });
    ctx.state.candidateList = [cand1, cand2];
    const r = ctx.window._gatherBreedingCardsFromState({ tierMode: 'tier_1_2' });
    ok('real-n: gathered 2 cards',                r && r.cards.length === 2);
    const out = ctx.window._buildBreedingCardsCombinedHTML(r.cards, {
      chrom: r.chrom, filter_meta: r.filter_meta,
    });
    ok('real-n: combined HTML > 20KB',            out.length > 20_000);
    ok('real-n: TOC has both inv1 + inv2',
       out.indexOf('LG28_inv1') >= 0 && out.indexOf('LG28_inv2') >= 0);
    ok('real-n: 2 .brc-print-page sections',
       (out.match(/class="brc-print-page"/g) || []).length === 2);
    ok('real-n: contains 60 (HOM_REF count)',     />60</.test(out));
    ok('real-n: contains 106 (HET count)',        />106</.test(out));
    const json = ctx.window._buildBreedingCardsJSONBundle(r.cards, {
      chrom: r.chrom, filter_meta: r.filter_meta,
    });
    ok('real-n: JSON bundle has 2 cards',         json.n_cards === 2);
    ok('real-n: JSON bundle is JSON.stringifiable',
       (() => { try { JSON.parse(JSON.stringify(json)); return true; } catch (_) { return false; } })());
  } else {
    ok('real-n fixture loaded', false, 'fixture not found');
  }
}

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
