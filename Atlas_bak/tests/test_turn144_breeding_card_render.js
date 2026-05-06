// =============================================================================
// turn 144 — Breeding-readiness card, Turn C (render layer)
// =============================================================================
// SPEC: specs_new_turn131/SPEC_per_candidate_breeding_readiness_card.md
//
// Builds on Turn B's _buildBreedingCard data. Adds five render helpers + a
// DOM dispatcher + a toolbar wirer:
//
//   - _renderBreedingCardFROHBarsSVG(burden, opts) — pure SVG bar chart
//   - _renderBreedingCardHTML(card, opts)          — pure body-HTML string
//   - _breedingCardPrintCSS()                       — pure print stylesheet
//   - _breedingCardPrintHTML(card, opts)           — full <!doctype> page
//   - openBreedingCardPrintWindow(cand, opts)      — DOM dispatcher
//   - _wireCandKaryoPrintCardBtn(cand)             — toolbar wiring
//
// Tests cover: source presence, window exports, structural HTML / SVG
// shapes, the print stylesheet, the auto-print wrapping, real-fixture
// rendering, and DOM dispatcher behavior under mocked window/document.
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
// Sandbox builder — same approach as test_turn143_breeding_card_compute.js,
// but extends the regex to swallow the Turn C block too.
// ============================================================================
function buildSandbox(extra) {
  const re = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 144 — Breeding-readiness card, Turn C \(render layer\)[\s\S]*?window\._wireCandKaryoPrintCardBtn = _wireCandKaryoPrintCardBtn;\s*\n}/;
  const m = html.match(re);
  if (!m) throw new Error('combined Turn-B/C region not found');
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  if (!normal) throw new Error('normalCDF region not found');
  const _ls = new Map();
  const sandbox = Object.assign({
    state: { data: null, cohortDiversity: null },
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

// Standard K=3 candidate fixture (same as Turn B)
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

// Build a card in the sandbox
function buildCard(ctx, cand) {
  return ctx.window._buildBreedingCard(cand);
}

// ============================================================================
// 1. Source-level definitions
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');
ok('_brHtmlEsc defined',                     /function _brHtmlEsc\(/.test(html));
ok('_renderBreedingCardFROHBarsSVG defined', /function _renderBreedingCardFROHBarsSVG\(/.test(html));
ok('_renderBreedingCardHTML defined',        /function _renderBreedingCardHTML\(/.test(html));
ok('_breedingCardPrintCSS defined',          /function _breedingCardPrintCSS\(/.test(html));
ok('_breedingCardPrintHTML defined',         /function _breedingCardPrintHTML\(/.test(html));
ok('openBreedingCardPrintWindow defined',    /function openBreedingCardPrintWindow\(/.test(html));
ok('_wireCandKaryoPrintCardBtn defined',     /function _wireCandKaryoPrintCardBtn\(/.test(html));
ok('_BR_KARYO_COLORS defined',               /const _BR_KARYO_COLORS\s*=/.test(html));
ok('_BR_SEVERITY_COLORS defined',            /const _BR_SEVERITY_COLORS\s*=/.test(html));
ok('Turn C banner present',
   /turn 144 — Breeding-readiness card, Turn C \(render layer\)/.test(html));

// ============================================================================
// 2. Window exports
// ============================================================================
console.log('\n=== 2. Window exports ===');
ok('window._renderBreedingCardFROHBarsSVG',
   /window\._renderBreedingCardFROHBarsSVG\s*=\s*_renderBreedingCardFROHBarsSVG/.test(html));
ok('window._renderBreedingCardHTML',
   /window\._renderBreedingCardHTML\s*=\s*_renderBreedingCardHTML/.test(html));
ok('window._breedingCardPrintCSS',
   /window\._breedingCardPrintCSS\s*=\s*_breedingCardPrintCSS/.test(html));
ok('window._breedingCardPrintHTML',
   /window\._breedingCardPrintHTML\s*=\s*_breedingCardPrintHTML/.test(html));
ok('window.openBreedingCardPrintWindow',
   /window\.openBreedingCardPrintWindow\s*=\s*openBreedingCardPrintWindow/.test(html));
ok('window._wireCandKaryoPrintCardBtn',
   /window\._wireCandKaryoPrintCardBtn\s*=\s*_wireCandKaryoPrintCardBtn/.test(html));

// ============================================================================
// 3. _brHtmlEsc behavior
// ============================================================================
console.log('\n=== 3. _brHtmlEsc behavior ===');
{
  const ctx = buildSandbox();
  const esc = ctx.window._renderBreedingCardHTML;
  // We don't expose _brHtmlEsc on window directly — test via output observation.
  // Instead, exercise via card.candidate.id containing HTML chars.
  const cand = k3Fixture().cand;
  cand.id = '<script>alert(1)</script>';
  ctx.state.data = { chrom: 'LG28', samples: k3Fixture().chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: k3Fixture().samples });
  const card = ctx.window._buildBreedingCard(cand);
  const html2 = ctx.window._renderBreedingCardHTML(card);
  ok('escapes < in id',   html2.indexOf('<script>alert') === -1);
  ok('encodes &lt; in id', html2.indexOf('&lt;script&gt;alert(1)&lt;/script&gt;') >= 0);
}

// ============================================================================
// 4. _renderBreedingCardHTML structural shape
// ============================================================================
console.log('\n=== 4. _renderBreedingCardHTML structural shape ===');
{
  const ctx = buildSandbox();
  const fx = k3Fixture();
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });
  const card = ctx.window._buildBreedingCard(fx.cand);
  const out = ctx.window._renderBreedingCardHTML(card);
  ok('returns string',                    typeof out === 'string' && out.length > 100);
  ok('has .brc-card root',                /class="brc-card"/.test(out));
  ok('has .brc-header',                   /class="brc-header"/.test(out));
  ok('has .brc-title',                    /class="brc-title"/.test(out));
  ok('title contains candidate id',       out.indexOf('LG28_test') >= 0);
  ok('chrom in meta',                     out.indexOf('LG28') >= 0);
  ok('span_mb shown',                     /span 2\.50 Mb/.test(out));
  ok('coords shown in Mb',                /14\.00–16\.50 Mb/.test(out));
  ok('K=3 shown',                         /K=3/.test(out));
  ok('Tier 1 shown',                      out.indexOf('Tier 1') >= 0);
  ok('confidence shown',                  /confidence: high/.test(out));
  ok('has .brc-section.brc-karyo',        /class="brc-section brc-karyo"/.test(out));
  ok('has .brc-chip-row',                 /class="brc-chip-row"/.test(out));
  ok('REF chip present',                  /HOM_REF/.test(out));
  ok('HET chip present',                  out.indexOf('HET') >= 0);
  ok('INV chip present',                  /HOM_INV/.test(out));
  ok('classification_summary present',    /K=3 · 2H\+1T/.test(out));
  ok('has .brc-section.brc-froh',         /class="brc-section brc-froh"/.test(out));
  ok('embeds an SVG inside froh section', /brc-svg-wrap/.test(out) && /<svg/.test(out));
  ok('has summary table',                 /class="brc-table"/.test(out));
  ok('table headers present',
     /<th>Arrangement<\/th>/.test(out) && /<th>n<\/th>/.test(out) && /<th>mean<\/th>/.test(out));
  ok('Wilcoxon line present',             /class="brc-wilcoxon"/.test(out));
  ok('Wilcoxon shows p =',                /<b>p = /.test(out));
  ok('has .brc-section.brc-ancestry',     /class="brc-section brc-ancestry"/.test(out));
  ok('ancestry K8 row present',           /<b>K1<\/b>|<b>K2<\/b>|<b>K3<\/b>/.test(out));
  ok('ancestry totals row present',       /class="brc-totals"/.test(out));
  ok('has .brc-section.brc-maf',          /class="brc-section brc-maf"/.test(out));
  ok('p(REF) shown',                      out.indexOf('p(REF)') >= 0);
  ok('p(INV) shown',                      out.indexOf('p(INV)') >= 0);
  ok('MAF label shown',                   /<b>MAF<\/b>/.test(out));
  ok('has .brc-section.brc-advice',       /class="brc-section brc-advice"/.test(out));
  ok('at least one advice item',          /class="brc-advice-item"/.test(out));
  ok('disclaimer present',                /class="brc-disclaimer"/.test(out));
  ok('disclaimer mentions human review',  /human review/.test(out));
  ok('footer present',                    /class="brc-footer"/.test(out));
  ok('footer shows schema',               /schema:/.test(out));
  ok('footer shows generated',            /generated:/.test(out));
}

// ============================================================================
// 5. _renderBreedingCardHTML — empty / fallback states
// ============================================================================
console.log('\n=== 5. _renderBreedingCardHTML fallback states ===');
{
  const ctx = buildSandbox();
  // null card
  const empty = ctx.window._renderBreedingCardHTML(null);
  ok('null card returns string',          typeof empty === 'string' && empty.length > 0);
  ok('null card uses brc-empty class',    /class="brc-card brc-empty"/.test(empty));
  ok('null card mentions No candidate',   /No candidate/.test(empty));
  // candidate without h_classification
  const fx = k3Fixture();
  const cand2 = {
    id: 'no_hc', chrom: 'LG28', start_bp: 1e6, end_bp: 2e6, K: 3,
    locked_labels: new Int8Array([0, 1, 2]),
  };
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });
  const card2 = ctx.window._buildBreedingCard(cand2);
  const html2 = ctx.window._renderBreedingCardHTML(card2);
  ok('no h_class: still renders',         typeof html2 === 'string' && html2.length > 100);
  ok('no h_class: header has id no_hc',   /no_hc/.test(html2));
  ok('no h_class: karyo empty msg',       /Karyotype not available/.test(html2));
  ok('no h_class: froh empty msg',
     /F<sub>ROH<\/sub> per-arrangement statistics/.test(html2) ||
     /F_ROH chart unavailable/.test(html2));
  ok('no h_class: advice still emitted',  /class="brc-advice-item"/.test(html2));
  // candidate with no cohort_diversity
  const ctx2 = buildSandbox();
  ctx2.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  // Don't load cohort_diversity
  const card3 = ctx2.window._buildBreedingCard(fx.cand);
  const html3 = ctx2.window._renderBreedingCardHTML(card3);
  ok('no cohort: still renders',           typeof html3 === 'string');
  ok('no cohort: shows no_cohort_diversity reason in body',
     /no_cohort_diversity/.test(html3));
  ok('no cohort: shows karyotype counts',  /class="brc-chip-row"/.test(html3));
  ok('no cohort: shows ancestry empty',    /Ancestry table not available/.test(html3));
}

// ============================================================================
// 6. _renderBreedingCardFROHBarsSVG structural shape
// ============================================================================
console.log('\n=== 6. _renderBreedingCardFROHBarsSVG structural shape ===');
{
  const ctx = buildSandbox();
  const fx = k3Fixture();
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });
  const card = ctx.window._buildBreedingCard(fx.cand);
  const svg = ctx.window._renderBreedingCardFROHBarsSVG(card.burden);
  ok('returns string',                     typeof svg === 'string' && svg.length > 50);
  ok('starts with <svg',                   svg.indexOf('<svg') === 0);
  ok('ends with </svg>',                   svg.trim().endsWith('</svg>'));
  ok('has xmlns',                          /xmlns="http:\/\/www\.w3\.org\/2000\/svg"/.test(svg));
  ok('has viewBox',                        /viewBox="0 0 \d+ \d+"/.test(svg));
  ok('has REF tick label',                 svg.indexOf('>REF<') >= 0);
  ok('has HET tick label',                 svg.indexOf('>HET<') >= 0);
  ok('has INV tick label',                 svg.indexOf('>INV<') >= 0);
  ok('contains REF n=3 label',             /n=3/.test(svg));
  ok('contains REF color (#3b6fb6)',       svg.indexOf('#3b6fb6') >= 0);
  ok('contains INV color (#d97a2c)',       svg.indexOf('#d97a2c') >= 0);
  ok('Wilcoxon p annotation present',      /REF vs INV: p =/.test(svg));
  ok('has at least 3 <rect>',              (svg.match(/<rect/g) || []).length >= 3);
  // Whisker lines
  ok('has whisker lines',                  (svg.match(/<line/g) || []).length >= 4);

  // Custom width/height
  const svg2 = ctx.window._renderBreedingCardFROHBarsSVG(card.burden, { width: 700, height: 250 });
  ok('honors custom width',                /viewBox="0 0 700 250"/.test(svg2) && /width="700"/.test(svg2));

  // Unavailable burden -> placeholder
  const svg3 = ctx.window._renderBreedingCardFROHBarsSVG({ available: false, reason: 'no_cohort_diversity' });
  ok('unavailable returns SVG',            svg3.indexOf('<svg') === 0);
  ok('unavailable shows reason',           /no_cohort_diversity/.test(svg3));

  // Null burden -> placeholder
  const svg4 = ctx.window._renderBreedingCardFROHBarsSVG(null);
  ok('null burden returns SVG placeholder', svg4.indexOf('<svg') === 0);
  ok('null burden shows unavailable',       /unavailable/.test(svg4));

  // Burden with HOM_MID -> 4 bars
  const candMid = {
    id: 'LG28_K4', chrom: 'LG28',
    start_bp: 1e6, end_bp: 2e6, K: 4,
    locked_labels: new Int8Array([0,0,1,1,1,2,2,3,3]),
    h_classification: {
      K: 4,
      bands: [
        { band_idx: 0, classification: 'HOM', median_pc1: -3.0 },
        { band_idx: 1, classification: 'HET', median_pc1:  0.0 },
        { band_idx: 2, classification: 'HOM', median_pc1:  1.0 },   // MID
        { band_idx: 3, classification: 'HOM', median_pc1:  3.0 },
      ],
      band_counts: { n_hom: 3, n_het: 1, n_ambiguous: 0, n_no_dosage: 0 },
      implied_regime: { consistency: 'multi_haplotype' },
    },
  };
  const ctx2 = buildSandbox();
  const samples4 = [];
  for (let i = 0; i < 9; i++) samples4.push({ sample_id: 'CGA' + (100+i), k8: 'K1', f_roh: 0.2 + i*0.01 });
  ctx2.state.data = { chrom: 'LG28', samples: samples4.map(s => ({ cga: s.sample_id })) };
  ctx2.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: samples4 });
  const cardMid = ctx2.window._buildBreedingCard(candMid);
  const svgMid = ctx2.window._renderBreedingCardFROHBarsSVG(cardMid.burden);
  ok('K=4 burden: MID label in SVG',       svgMid.indexOf('>MID<') >= 0);
  ok('K=4 burden: 4 bar slots',            (svgMid.match(/text-anchor="middle"[^>]*>(REF|HET|INV|MID)</g) || []).length >= 4);
}

// ============================================================================
// 7. _breedingCardPrintCSS contents
// ============================================================================
console.log('\n=== 7. _breedingCardPrintCSS contents ===');
{
  const ctx = buildSandbox();
  const css = ctx.window._breedingCardPrintCSS();
  ok('css is non-empty string',            typeof css === 'string' && css.length > 200);
  ok('has @page A4',                       /@page\s*\{[^}]*size:\s*A4/.test(css));
  ok('has 14mm margin',                    /margin:\s*14mm/.test(css));
  ok('targets .brc-card',                  /\.brc-card/.test(css));
  ok('targets .brc-section',               /\.brc-section/.test(css));
  ok('targets .brc-table',                 /\.brc-table/.test(css));
  ok('targets .brc-advice-item',           /\.brc-advice-item/.test(css));
  ok('font-size in pt',                    /font-size:\s*\d+(\.\d+)?pt/.test(css));
  ok('page-break-inside: avoid present',
     /page-break-inside:\s*avoid/.test(css) && /break-inside:\s*avoid/.test(css));
  ok('@media print rule present',          /@media print/.test(css));
  ok('hides .brc-no-print in @media print',
     /\.brc-no-print[^{]*\{[^}]*display:\s*none/.test(css));
}

// ============================================================================
// 8. _breedingCardPrintHTML wrapping
// ============================================================================
console.log('\n=== 8. _breedingCardPrintHTML wrapping ===');
{
  const ctx = buildSandbox();
  const fx = k3Fixture();
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });
  const card = ctx.window._buildBreedingCard(fx.cand);
  const out = ctx.window._breedingCardPrintHTML(card);
  ok('starts with !doctype',               /^<!doctype html>/i.test(out));
  ok('has <html>',                         /<html[^>]*>/.test(out));
  ok('has <meta charset',                  /<meta charset="utf-8"/.test(out));
  ok('title contains candidate id',        /<title>[^<]*LG28_test/.test(out));
  ok('title mentions breeding-readiness card', /breeding-readiness card<\/title>/.test(out));
  ok('inline <style> tag present',         /<style>[\s\S]*<\/style>/.test(out));
  ok('embeds @page A4',                    /@page\s*\{[^}]*size:\s*A4/.test(out));
  ok('embeds the body via brc-card',       /class="brc-card"/.test(out));
  ok('embeds Wilcoxon section',            /brc-wilcoxon/.test(out));
  ok('auto-print script included',
     /window\.addEventListener\("load"/.test(out) && /window\.print\(\)/.test(out));
  ok('print hint included by default',     /class="brc-print-hint/.test(out));

  // auto_print = false
  const out2 = ctx.window._breedingCardPrintHTML(card, { auto_print: false });
  ok('auto_print=false omits print()',     out2.indexOf('window.print()') === -1);
  ok('auto_print=false still has body',    /class="brc-card"/.test(out2));

  // show_hint = false
  const out3 = ctx.window._breedingCardPrintHTML(card, { show_hint: false });
  ok('show_hint=false omits hint',         /class="brc-print-hint/.test(out3) === false);

  // null card still wraps cleanly
  const out4 = ctx.window._breedingCardPrintHTML(null);
  ok('null card: returns valid HTML doc',  /^<!doctype html>/i.test(out4) && /<\/html>/.test(out4));
  ok('null card: title generic',           /breeding-readiness card<\/title>/.test(out4));
  ok('null card: body brc-empty',          /class="brc-card brc-empty"/.test(out4));
}

// ============================================================================
// 9. openBreedingCardPrintWindow with mocked DOM
// ============================================================================
console.log('\n=== 9. openBreedingCardPrintWindow with mocked DOM ===');
{
  // Build a sandbox where window.open + document.* are mocked to capture
  // what the dispatcher does.
  const mockWriteCalls = [];
  const mockDocument = {
    open:  () => { mockWriteCalls.push({ kind: 'open' }); },
    write: (s) => { mockWriteCalls.push({ kind: 'write', html: s }); },
    close: () => { mockWriteCalls.push({ kind: 'close' }); },
    location: { origin: 'http://localhost', pathname: '/atlas' },
    querySelector: () => null,
  };
  let openCalls = 0;
  const mockWindow = {};
  // Provide window.open that returns a fake child window
  function makeOpenedWindow() {
    return { document: { open: mockDocument.open, write: mockDocument.write, close: mockDocument.close } };
  }
  const _ls = new Map();
  const ctx = vm.createContext({
    state: { data: null, cohortDiversity: null },
    window: mockWindow,
    document: mockDocument,
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
    alert: (msg) => { ctx.__alerts && ctx.__alerts.push(msg); },
    console,
  });
  // Wire window.open AFTER context creation
  mockWindow.open = function () {
    openCalls++;
    return makeOpenedWindow();
  };
  ctx.__alerts = [];
  const reBoth = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 144 — Breeding-readiness card, Turn C \(render layer\)[\s\S]*?window\._wireCandKaryoPrintCardBtn = _wireCandKaryoPrintCardBtn;\s*\n}/;
  const m = html.match(reBoth);
  if (!m) throw new Error('combined region not found in mock-DOM test');
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  vm.runInContext(normal[0] + '\n' + m[0], ctx);

  // Set up data + candidate
  const fx = k3Fixture();
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });

  // Dispatch
  const result = ctx.window.openBreedingCardPrintWindow(fx.cand);
  ok('dispatch returns truthy on success',  !!result);
  ok('window.open called once',             openCalls === 1);
  ok('document.open called',                mockWriteCalls.some(c => c.kind === 'open'));
  ok('document.write called',               mockWriteCalls.some(c => c.kind === 'write'));
  ok('document.close called',               mockWriteCalls.some(c => c.kind === 'close'));
  const wrote = mockWriteCalls.find(c => c.kind === 'write');
  ok('written HTML starts with !doctype',   wrote && /^<!doctype html>/i.test(wrote.html));
  ok('written HTML contains brc-card',      wrote && /class="brc-card"/.test(wrote.html));
  ok('written HTML contains LG28_test',     wrote && /LG28_test/.test(wrote.html));

  // Null candidate -> alert + false
  ctx.__alerts = [];
  const r2 = ctx.window.openBreedingCardPrintWindow(null);
  ok('null cand returns false',             r2 === false);
  ok('null cand triggers alert',            ctx.__alerts.length === 1);

  // Popup blocked: window.open returns null
  mockWindow.open = function () { return null; };
  ctx.__alerts = [];
  const r3 = ctx.window.openBreedingCardPrintWindow(fx.cand);
  ok('popup blocked returns false',         r3 === false);
  ok('popup blocked triggers alert',        ctx.__alerts.length === 1 &&
                                             /popup/i.test(ctx.__alerts[0]));
}

// ============================================================================
// 10. _wireCandKaryoPrintCardBtn — DOM idempotence
// ============================================================================
console.log('\n=== 10. _wireCandKaryoPrintCardBtn DOM idempotence ===');
{
  // Mock just enough DOM to verify the wiring
  const elements = new Map();
  let idCounter = 1;
  function makeEl(tag) {
    const el = {
      tagName: tag.toUpperCase(),
      id: '', className: '', title: '', textContent: '',
      _children: [], _listeners: {},
      style: {}, dataset: {},
      appendChild: function (c) { this._children.push(c); c._parent = this; return c; },
      insertBefore: function (c, ref) {
        const idx = this._children.indexOf(ref);
        if (idx >= 0) this._children.splice(idx, 0, c);
        else this._children.push(c);
        c._parent = this;
        return c;
      },
      addEventListener: function (e, h) {
        if (!this._listeners[e]) this._listeners[e] = [];
        this._listeners[e].push(h);
      },
      querySelector: function (sel) {
        // Very narrow CSS support: handle '#xxx'
        if (sel.indexOf('#') === 0) {
          const wantId = sel.slice(1);
          const dfs = (node) => {
            if (node.id === wantId) return node;
            for (const ch of node._children) {
              const r = dfs(ch);
              if (r) return r;
            }
            return null;
          };
          return dfs(this);
        }
        return null;
      },
      getAttribute: function () { return null; },
    };
    return el;
  }

  const root  = makeEl('body');
  const pane  = makeEl('div');  pane.id = 'candKaryoPane';
  const tbar  = makeEl('div');  tbar.className = 'ck-toolbar';
  const expBtn = makeEl('button'); expBtn.id = 'ckExportTSV';
  const info   = makeEl('span');   info.id = 'ckInfo';
  pane.appendChild(tbar);
  tbar.appendChild(expBtn);
  tbar.appendChild(info);
  root.appendChild(pane);

  const docMock = {
    body: root,
    createElement: (tag) => makeEl(tag),
    querySelector: function (sel) {
      // Support '#candKaryoPane .ck-toolbar'
      if (sel === '#candKaryoPane .ck-toolbar') return tbar;
      return null;
    },
  };

  const _ls = new Map();
  const sandbox = {
    state: { data: null, cohortDiversity: null },
    window: {},
    document: docMock,
    localStorage: {
      getItem: (k) => _ls.has(k) ? _ls.get(k) : null,
      setItem: (k, v) => _ls.set(k, String(v)),
      removeItem: (k) => _ls.delete(k),
    },
    alert: () => {},
    console,
  };
  const ctx = vm.createContext(sandbox);
  const reBoth = /\/\/ turn 142 — cohort_diversity_v1 loader[\s\S]*?\/\/ turn 144 — Breeding-readiness card, Turn C \(render layer\)[\s\S]*?window\._wireCandKaryoPrintCardBtn = _wireCandKaryoPrintCardBtn;\s*\n}/;
  const m = html.match(reBoth);
  const normalRe = /\/\/ Abramowitz & Stegun 7\.1\.26 normal CDF approximation\nfunction normalCDF[\s\S]*?\n\}/;
  const normal = html.match(normalRe);
  vm.runInContext(normal[0] + '\n' + m[0], ctx);

  const cand = k3Fixture().cand;
  ctx.window._wireCandKaryoPrintCardBtn(cand);

  // Find the new button via querySelector on the toolbar
  const btn = tbar.querySelector('#ckPrintBreedingCard');
  ok('button created with id ckPrintBreedingCard', !!btn);
  ok('button class includes export',                btn && btn.className === 'export');
  ok('button title set',                            btn && btn.title.length > 0);
  ok('button text contains print breeding card',    btn && /print breeding card/i.test(btn.textContent));
  // It should be inserted before #ckInfo
  const idx = tbar._children.indexOf(btn);
  const idxInfo = tbar._children.indexOf(info);
  ok('button inserted before #ckInfo',              idx >= 0 && idx < idxInfo);

  // Idempotent — second call doesn't add a duplicate
  const before = tbar._children.length;
  ctx.window._wireCandKaryoPrintCardBtn(cand);
  ok('idempotent: child count unchanged',           tbar._children.length === before);

  // No-op when no toolbar present
  const sandbox2 = {
    state: { data: null, cohortDiversity: null },
    window: {},
    document: { querySelector: () => null, body: makeEl('body') },
    localStorage: { getItem: () => null, setItem: () => {}, removeItem: () => {} },
    console,
  };
  const ctx2 = vm.createContext(sandbox2);
  vm.runInContext(normal[0] + '\n' + m[0], ctx2);
  let threw = false;
  try { ctx2.window._wireCandKaryoPrintCardBtn(cand); } catch (_) { threw = true; }
  ok('no-toolbar: does not throw',                  !threw);
}

// ============================================================================
// 11. Real-fixture render
// ============================================================================
console.log('\n=== 11. Real-fixture render ===');
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
    const cand = {
      id: 'LG28_real_n', chrom: 'LG28',
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
    const card = ctx.window._buildBreedingCard(cand);
    const body = ctx.window._renderBreedingCardHTML(card);
    ok('real-n: body string length > 5KB',   body.length > 5000);
    ok('real-n: contains 60 (HOM_REF count)', />60</.test(body));
    ok('real-n: contains 106 (HET count)',    />106</.test(body));
    ok('real-n: at least 1 K8 row',           /<b>K1<\/b>|<b>K2<\/b>|<b>K3<\/b>|<b>K4<\/b>|<b>K5<\/b>|<b>K6<\/b>|<b>K7<\/b>|<b>K8<\/b>/.test(body));
    const full = ctx.window._breedingCardPrintHTML(card);
    ok('real-n: full doc starts !doctype',    /^<!doctype html>/i.test(full));
    ok('real-n: full doc has @page A4',       /@page\s*\{[^}]*size:\s*A4/.test(full));
    ok('real-n: full doc embeds 226-derived counts',
       />60</.test(full) && />106</.test(full));
  } else {
    ok('real-n fixture loaded',               false, 'fixture not found at ' + FIXTURE_PATH);
  }
}

// ============================================================================
// 12. Severity styling on advice items
// ============================================================================
console.log('\n=== 12. Severity styling on advice items ===');
{
  const ctx = buildSandbox();
  // Build a card whose burden makes Wilcoxon p < 0.05 with INV >> REF →
  // froh_asymmetric strong rule fires.
  const refF = [0.10, 0.11, 0.09, 0.12, 0.10];
  const hetF = [0.20, 0.21, 0.22];
  const invF = [0.40, 0.41, 0.42, 0.39, 0.40];
  const fx = k3Fixture(refF, hetF, invF);
  ctx.state.data = { chrom: 'LG28', samples: fx.chromSamples };
  ctx.window._storeCohortDiversity({ tool: 'cohort_diversity_v1', schema_version: 1, samples: fx.samples });
  const card = ctx.window._buildBreedingCard(fx.cand);
  const body = ctx.window._renderBreedingCardHTML(card);
  ok('strong severity in HTML body',        /STRONG/.test(body));
  ok('info severity in HTML body',          /INFO/.test(body));
  ok('advice item has border-left-color',   /border-left-color:\s*#d97a2c/.test(body));
  // strong palette uses #d97a2c (orange)
}

// ============================================================================
// 13. Atlas-side wiring: button hookup + tutorials placeholder
// ============================================================================
// These are source-level grep assertions on Inversion_atlas.html — no
// sandbox needed. They confirm the wiring + placeholder splices took.
console.log('\n=== 13. Atlas-side wiring + tutorials placeholder ===');
ok('renderCandidateKaryotype wires _wireCandKaryoPrintCardBtn',
   /try \{ _wireCandKaryoPrintCardBtn\(c\); \} catch/.test(html));
ok('wiring sits after expBtn handler',
   (() => {
     const expIdx = html.indexOf("expBtn.addEventListener('click', () => exportKaryotypeTSV");
     const wireIdx = html.indexOf('_wireCandKaryoPrintCardBtn(c)');
     return expIdx >= 0 && wireIdx > expIdx && (wireIdx - expIdx) < 800;
   })());
ok('tutorials section heading present',
   /<h4 id="tutorials-section">Tutorials<\/h4>/.test(html));
ok('tutorials grid container present',
   /id="tutorialsGrid"/.test(html));
ok('tutorials grid uses responsive auto-fill',
   /grid-template-columns:\s*repeat\(auto-fill,\s*minmax\(260px/.test(html));
ok('discover_first_inversion card present',
   /data-tutorial-id="discover_first_inversion"/.test(html));
ok('boundary_refinement card present',
   /data-tutorial-id="boundary_refinement"/.test(html));
ok('breeding_card tutorial card present',
   /data-tutorial-id="breeding_card"/.test(html));
ok('manuscript_bundle tutorial card present',
   /data-tutorial-id="manuscript_bundle"/.test(html));
ok('cli_inversion_toolkit card present',
   /data-tutorial-id="cli_inversion_toolkit"/.test(html));
ok('cli_te_density card present',
   /data-tutorial-id="cli_te_density"/.test(html));
ok('cli_step_d17_boundary card present',
   /data-tutorial-id="cli_step_d17_boundary"/.test(html));
ok('cli_module_2b_popstruct card present',
   /data-tutorial-id="cli_module_2b_popstruct"/.test(html));
ok('all 8 cards default to data-status="pending"',
   (html.match(/class="tutorial-card"[^>]*data-status="pending"/g) || []).length === 8);
ok('coming-soon pill rendered for each card',
   (html.match(/coming soon<\/span>/g) || []).length >= 8);
ok('authoring note present',
   /Authoring note:/.test(html) && /data-status="available"/.test(html));
ok('30-seconds reference present (UX pattern continuity)',
   /30 seconds to figure 3/.test(html));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
