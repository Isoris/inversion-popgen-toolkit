// Smoke test for the cross-species page (page16) — schema v2 with synteny.
// Loads the atlas in JSDOM, executes the inline script, and exercises both
// the breakpoints layer and the macro-synteny layer.

const { JSDOM } = require('jsdom');
const fs   = require('fs');

const ATLAS_PATH = '/home/claude/work/Inversion_atlas.html';
const html = fs.readFileSync(ATLAS_PATH, 'utf-8');

const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window } = dom;
const { document } = window;

function fail(msg) { console.error('[smoke] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[smoke] ok  :', msg); }

// ---- 1. Page16 DOM ----
const page16 = document.getElementById('page16');
if (!page16) fail('page16 div not found');
ok('page16 div present');

const tabBtn = document.querySelector('#tabBar button[data-page="page16"]');
if (!tabBtn) fail('page16 tab button missing');
ok('page16 tab button present: "' + (tabBtn.textContent || '').trim().replace(/\s+/g, ' ') + '"');

// All breakpoint-layer + synteny-layer helpers should be exposed
for (const name of [
  // breakpoint layer
  '_isCrossSpeciesJSON', '_storeCrossSpecies', '_persistCrossSpecies',
  '_restoreCrossSpecies', '_clearCrossSpecies', '_renderCrossSpeciesPage',
  '_renderCrossSpeciesToolbar', '_renderCrossSpeciesCatalogue', '_renderCrossSpeciesFocus',
  '_csFilteredSorted', '_csIsSpalaxFlagged', '_csGarAllTeFlankMean',
  '_crossSpeciesStep', '_crossSpeciesToggleEventFilter', '_crossSpeciesClearFilters',
  '_wireCrossSpeciesKeys',
  // synteny layer (v4 turn 75)
  '_csComputeSynteny', '_csSyntenyEdgesByChrom', '_csDistanceToNearestEdge',
  '_csInversionContexts', '_csPermutationTest', '_csBuildSankeySvg',
  '_csBuildSyntenySection', '_renderCrossSpeciesSynteny', '_csNaturalSortChroms',
  '_csMulberry32', '_csSeedFromString', '_csBinomialTailP',
]) {
  if (typeof window[name] !== 'function') fail(name + ' missing from window');
}
ok('all breakpoint + synteny helpers present');

// ---- 2. Synthesize a cs_breakpoints_v1 v2 JSON ----
// Layout designed to exercise:
//   - one_to_one: gar LG28 ↔ mac LG28 (mostly conserved, with the inversion)
//   - fission_in_target: gar LG07 → mac LG07 + LG10 (translocation)
//   - one_to_one + 2nd inversion: gar LG15 ↔ mac LG15 (inside a single block)
const goodJson = {
  tool: 'cross_species_breakpoints_v1',
  schema_version: 2,
  generated_at: '2026-04-30T00:00:00Z',
  species_query:  { name: 'Clarias gariepinus',  haplotype: 'fClaHyb_Gar_LG' },
  species_target: { name: 'Clarias macrocephalus', haplotype: 'fClaHyb_Mac_LG' },
  input_paf: { path: '/scratch/test.paf', sha256: 'deadbeef' },
  params: { min_mapq: 40, min_block_bp: 50000, cluster_radius_bp: 50000, flank_bp: 100000, classes_annotated: ['all_TE'] },
  n_breakpoints: 3,
  n_by_event_type: { inversion: 2, fission_or_fusion: 1 },
  breakpoints: [
    {
      id: 'cs_bp_0001', event_type: 'inversion',
      gar_chr: 'LG28', gar_pos_start: 3_000_000, gar_pos_end: 3_010_000, gar_pos_mb: 3.005,
      prev_block: { mac_chr: 'LG28', mac_start_bp: 0,         mac_end_bp: 3_000_000,  strand: '+', block_size_bp: 3_000_000, mapping_quality: 60 },
      next_block: { mac_chr: 'LG28', mac_start_bp: 3_000_000, mac_end_bp: 5_000_000,  strand: '-', block_size_bp: 2_000_000, mapping_quality: 60 },
      flanking_repeat_density_gar: { all_TE: { mean: 0.78, max: 0.91, n_windows: 4 } },
      flanking_repeat_density_mac: {
        prev: { mac_chr: 'LG28', anchor_bp: 3_000_000, by_class: { all_TE: { mean: 0.42, max: 0.50, n_windows: 4 } } },
        next: { mac_chr: 'LG28', anchor_bp: 3_000_000, by_class: { all_TE: { mean: 0.41, max: 0.49, n_windows: 4 } } },
      },
      manuscript_note: 'Spalax-style enrichment: gar_flank all_TE mean 0.78 vs chrom mean 0.45 (1.7x)',
      candidate_overlap: ['CAND001'],
    },
    {
      id: 'cs_bp_0002', event_type: 'inversion',
      gar_chr: 'LG28', gar_pos_start: 5_000_000, gar_pos_end: 5_010_000, gar_pos_mb: 5.005,
      prev_block: { mac_chr: 'LG28', mac_start_bp: 3_000_000, mac_end_bp: 5_000_000,  strand: '-', block_size_bp: 2_000_000, mapping_quality: 60 },
      next_block: { mac_chr: 'LG28', mac_start_bp: 5_000_000, mac_end_bp: 10_000_000, strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
      flanking_repeat_density_gar: { all_TE: { mean: 0.71, max: 0.88, n_windows: 4 } },
      flanking_repeat_density_mac: {
        prev: { mac_chr: 'LG28', anchor_bp: 5_000_000, by_class: { all_TE: { mean: 0.40, max: 0.48, n_windows: 4 } } },
        next: { mac_chr: 'LG28', anchor_bp: 5_000_000, by_class: { all_TE: { mean: 0.39, max: 0.46, n_windows: 4 } } },
      },
      candidate_overlap: ['CAND001'],
    },
    {
      id: 'cs_bp_0003', event_type: 'translocation_or_fission', event_type_refined: 'fission_or_fusion',
      gar_chr: 'LG07', gar_pos_start: 5_000_000, gar_pos_end: 5_010_000, gar_pos_mb: 5.005,
      prev_block: { mac_chr: 'LG07', mac_start_bp: 0,          mac_end_bp: 5_000_000,  strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
      next_block: { mac_chr: 'LG10', mac_start_bp: 20_000_000, mac_end_bp: 25_000_000, strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
      flanking_repeat_density_gar: { all_TE: { mean: 0.45, max: 0.51, n_windows: 4 } },
      flanking_repeat_density_mac: {
        prev: { mac_chr: 'LG07', anchor_bp: 5_000_000, by_class: { all_TE: { mean: 0.40, max: 0.45, n_windows: 4 } } },
        next: { mac_chr: 'LG10', anchor_bp: 20_000_000, by_class: { all_TE: { mean: 0.38, max: 0.43, n_windows: 4 } } },
      },
    },
  ],
  // v2 additions:
  n_synteny_blocks: 6,
  synteny_blocks: [
    // LG28: 3 blocks on mac LG28 (mostly 1:1, contains the inversion)
    { gar_chr: 'LG28', gar_start: 0,         gar_end: 3_000_000,  mac_chr: 'LG28', mac_start: 0,          mac_end: 3_000_000,  strand: '+', block_size_bp: 3_000_000, mapping_quality: 60 },
    { gar_chr: 'LG28', gar_start: 3_000_000, gar_end: 5_000_000,  mac_chr: 'LG28', mac_start: 3_000_000,  mac_end: 5_000_000,  strand: '-', block_size_bp: 2_000_000, mapping_quality: 60 },
    { gar_chr: 'LG28', gar_start: 5_000_000, gar_end: 10_000_000, mac_chr: 'LG28', mac_start: 5_000_000,  mac_end: 10_000_000, strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
    // LG07: 2 blocks on different mac chroms (fission in target)
    { gar_chr: 'LG07', gar_start: 0,         gar_end: 5_000_000,  mac_chr: 'LG07', mac_start: 0,          mac_end: 5_000_000,  strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
    { gar_chr: 'LG07', gar_start: 5_000_000, gar_end: 10_000_000, mac_chr: 'LG10', mac_start: 20_000_000, mac_end: 25_000_000, strand: '+', block_size_bp: 5_000_000, mapping_quality: 60 },
    // LG15: single mac LG15 block (1:1 conserved)
    { gar_chr: 'LG15', gar_start: 0,         gar_end: 12_000_000, mac_chr: 'LG15', mac_start: 0,          mac_end: 12_000_000, strand: '+', block_size_bp: 12_000_000, mapping_quality: 60 },
  ],
  chrom_lengths_query:  { LG07: 10_000_000, LG15: 12_000_000, LG28: 10_000_000 },
  chrom_lengths_target: { LG07: 5_000_000,  LG10: 25_000_000, LG15: 12_000_000, LG28: 10_000_000 },
};
if (!window._isCrossSpeciesJSON(goodJson)) fail('rejected a valid JSON');
if (!window._storeCrossSpecies(goodJson)) fail('_storeCrossSpecies returned false');
ok('stored v2 JSON with 6 synteny blocks');

// Activate page16 so the renderer's DOM lookups succeed
document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
page16.classList.add('active');
window._renderCrossSpeciesPage();
ok('_renderCrossSpeciesPage() did not throw');

// ---- 3. Synteny computation ----
const synteny = window._csComputeSynteny();
if (!synteny) fail('_csComputeSynteny returned null');
ok('synteny computed: ' + synteny.flows.length + ' flows');

// LG28 should be one_to_one (single mac partner LG28, that mac partner has only LG28 as gar)
if (synteny.gar_class.LG28 !== 'one_to_one')
  fail('expected LG28 → one_to_one, got: ' + synteny.gar_class.LG28);
ok('LG28 classified as one_to_one');

// LG07 should be fission_in_target (1 gar → 2 mac)
if (synteny.gar_class.LG07 !== 'fission_in_target')
  fail('expected LG07 → fission_in_target, got: ' + synteny.gar_class.LG07);
ok('LG07 classified as fission_in_target');

// LG15 should be one_to_one
if (synteny.gar_class.LG15 !== 'one_to_one')
  fail('expected LG15 → one_to_one, got: ' + synteny.gar_class.LG15);
ok('LG15 classified as one_to_one');

// ---- 4. Synteny edges ----
const edges = window._csSyntenyEdgesByChrom();
if (!edges || !edges.LG28) fail('synteny edges missing');
// LG28 has blocks at 0-3M, 3M-5M, 5M-10M → unique edges: 0, 3M, 5M, 10M
const lg28edges = edges.LG28;
if (lg28edges.length !== 4) fail('expected 4 unique LG28 edges, got ' + lg28edges.length + ': ' + JSON.stringify(lg28edges));
ok('LG28 edges: ' + lg28edges.map(e => e/1e6 + 'M').join(', '));

// Distance check
const d = window._csDistanceToNearestEdge('LG28', 4_500_000);   // halfway between 3M and 5M edges
if (d !== 500_000) fail('expected dist 500k from 4.5M to nearest edge, got ' + d);
ok('distance check: 4.5M → 500kb to nearest edge');

// ---- 5. Inversion contexts (synthesize a candidate list) ----
// We can't directly mutate state.candidateList from outside, but the helpers
// _csInversionContexts read from state.candidateList. JSDOM exposes functions
// but not the const state — so we test by directly invoking the underlying
// computation with a synthetic context list. The function operates on
// state.candidateList which is empty in this fresh JSDOM. Verify the
// "no candidates" path returns an empty array, then move on.
const ctxs = window._csInversionContexts();
if (!Array.isArray(ctxs)) fail('_csInversionContexts should return array (got ' + typeof ctxs + ')');
ok('_csInversionContexts returns array (n=' + ctxs.length + ', no candidate list set)');

// ---- 6. Sankey rendering ----
const sankey = window._csBuildSankeySvg();
if (!sankey || !sankey.includes('<svg')) fail('Sankey SVG missing');
if (!sankey.includes('LG28')) fail('Sankey missing LG28 label');
if (!sankey.includes('LG07')) fail('Sankey missing LG07 label');
// Class band labels should appear (one_to_one and fission_in_target both have flow)
if (!sankey.includes('1:1 conserved')) fail('Sankey missing one_to_one label');
if (!sankey.includes('fission in target')) fail('Sankey missing fission_in_target label');
ok('Sankey SVG renders with chrom + class labels');

// ---- 7. Synteny section in page DOM ----
const synSlot = document.getElementById('csSynteny');
if (!synSlot) fail('csSynteny div missing');
if (synSlot.style.display === 'none') fail('csSynteny hidden after render');
if (!synSlot.innerHTML.includes('Macro-synteny')) fail('synteny section title missing');
if (!synSlot.innerHTML.includes('1:1 conserved')) fail('synteny legend missing');
if (!synSlot.innerHTML.includes('Per-chromosome relationships')) fail('per-chrom table missing');
ok('synteny section rendered with Sankey + per-chrom table');

// ---- 8. Permutation test (no candidates → empty result, but should not throw) ----
const permResult = window._csPermutationTest();
if (!permResult) fail('permutation test returned null');
if (!permResult.global) fail('permutation test missing global');
if (permResult.global.n_candidates !== 0) fail('expected 0 candidates, got ' + permResult.global.n_candidates);
ok('permutation test runs cleanly with empty candidate list');

// ---- 8b. Permutation test with mock candidates ----
// Inject candidates directly into the module-scope `state` via the window
// object. In the atlas, `state` is a `const` declared at module top — JSDOM
// runs the script via runScripts:'dangerously', so the const is on the
// inline-script scope, not on window. But functions that close over `state`
// were also declared in that scope, so `_csComputeSynteny` etc. all read
// from the same const. We use window.eval() to mutate it from the inside.
window.eval(`
  state.candidateList = [
    // CAND_A: inside a single LG28 block (3M-5M region) — should be
    // "inside conserved 1:1 syntenic block" → recent within-lineage
    { id: 'CAND_A', chrom: 'LG28', start_bp: 3500000, end_bp: 4500000 },
    // CAND_B: spans the LG28 3M and 5M edges (one breakpoint at each) —
    // should be "spans 3 blocks" → ancient/fragile architecture reuse
    { id: 'CAND_B', chrom: 'LG28', start_bp: 2900000, end_bp: 5100000 },
    // CAND_C: on LG07 (fission chrom), straddles the 5M synteny edge
    // (where the gar→mac partner switches from LG07 to LG10)
    { id: 'CAND_C', chrom: 'LG07', start_bp: 4500000, end_bp: 5500000 },
    // CAND_D: deep inside LG15 (1:1 conserved) — should be inside-single-block
    { id: 'CAND_D', chrom: 'LG15', start_bp: 4000000, end_bp: 5000000 },
  ];
  // Invalidate caches so contexts get recomputed
  state._csInversionContextCache = null;
`);

const ctxsWithCands = window._csInversionContexts();
if (!Array.isArray(ctxsWithCands)) fail('_csInversionContexts returned non-array after candidate injection');
if (ctxsWithCands.length !== 4) fail('expected 4 contexts, got ' + ctxsWithCands.length);
ok('inversion contexts computed for 4 mock candidates');

// CAND_A: inside_single_block, chrom_class one_to_one
const cA = ctxsWithCands.find(c => c.candidate_id === 'CAND_A');
if (!cA.inside_single_block) fail('CAND_A should be inside single block, got spanning_n_blocks=' + cA.spanning_n_blocks);
if (cA.chrom_class !== 'one_to_one') fail('CAND_A chrom_class should be one_to_one');
if (!cA.interpretation.includes('within-lineage')) fail('CAND_A interp wrong: ' + cA.interpretation);
ok('CAND_A: inside conserved 1:1 block → ' + cA.interpretation);

// CAND_B: spans 3 blocks (overlaps 0-3M, 3M-5M, 5M-10M)
const cB = ctxsWithCands.find(c => c.candidate_id === 'CAND_B');
if (cB.spanning_n_blocks !== 3) fail('CAND_B should span 3 blocks, got ' + cB.spanning_n_blocks);
if (cB.inside_single_block) fail('CAND_B should not be inside single block');
if (!cB.interpretation.includes('coincides with a synteny edge'))
  fail('CAND_B should mention synteny edge coincidence: ' + cB.interpretation);
ok('CAND_B: spans 3 LG28 blocks → ' + cB.interpretation);

// CAND_C: on LG07 (fission), should mention fusion/fission in interp
const cC = ctxsWithCands.find(c => c.candidate_id === 'CAND_C');
if (cC.chrom_class !== 'fission_in_target') fail('CAND_C chrom_class wrong: ' + cC.chrom_class);
if (!cC.interpretation.includes('ancient rearrangement boundary'))
  fail('CAND_C should flag ancient rearrangement reuse: ' + cC.interpretation);
ok('CAND_C: on fission chrom → ' + cC.interpretation);

// CAND_D: inside LG15 single block
const cD = ctxsWithCands.find(c => c.candidate_id === 'CAND_D');
if (!cD.inside_single_block) fail('CAND_D should be inside single block');
if (cD.chrom_class !== 'one_to_one') fail('CAND_D chrom_class should be one_to_one');
ok('CAND_D: inside LG15 1:1 block → ' + cD.interpretation);

// Distances should be sensible
// CAND_A is at 3.5-4.5M with edges at {0, 3M, 5M, 10M}. Closest edge to start
// (3.5M) is 3M (dist=500kb); closest to end (4.5M) is 5M (dist=500kb).
if (cA.dist_min_edge_bp !== 500000) fail('CAND_A dist_min should be 500kb, got ' + cA.dist_min_edge_bp);
ok('CAND_A dist_min_edge = 500 kb (correct for breakpoints in middle of 3M-5M block)');

// Now run the permutation test with real candidates
const permWithCands = window._csPermutationTest({ n_perm: 1000 });
if (!permWithCands) fail('permutation test returned null with candidates');
if (permWithCands.global.n_valid !== 4) fail('expected n_valid=4, got ' + permWithCands.global.n_valid);
if (permWithCands.per_candidate.length !== 4) fail('expected 4 per-candidate results');
// Each candidate should have a p_value
for (const r of permWithCands.per_candidate) {
  if (r.p_value == null) fail('candidate ' + r.candidate_id + ' has null p_value');
  if (r.p_value < 0 || r.p_value > 1) fail('candidate ' + r.candidate_id + ' p_value out of [0,1]: ' + r.p_value);
}
ok('permutation test with 4 candidates: n_valid=' + permWithCands.global.n_valid +
   ', enrichment=' + (permWithCands.global.enrichment_ratio || 0).toFixed(2) +
   ', binomial_p=' + permWithCands.global.binomial_p.toFixed(3));

// Render the result HTML and check it doesn't throw
const permHtml = window._csBuildPermResultHtml(permWithCands);
if (!permHtml.includes('Global enrichment test')) fail('perm result HTML missing summary');
if (!permHtml.includes('CAND_')) fail('perm result HTML missing per-candidate rows');
ok('_csBuildPermResultHtml renders cleanly');

// Re-render synteny section with candidates loaded → should now show
// the inversion context table and Run button
window._renderCrossSpeciesSynteny();
const synSlotWithCands = document.getElementById('csSynteny');
if (!synSlotWithCands.innerHTML.includes('Inversion candidates in synteny context'))
  fail('inversion-context table missing after candidate injection');
if (!synSlotWithCands.innerHTML.includes('csPermRunBtn'))
  fail('run-permutation button missing after candidate injection');
ok('synteny section renders inversion table + perm button when candidates present');

// ---- 9. Permutation helpers ----
// Mulberry32 determinism
const r1 = window._csMulberry32(42);
const r2 = window._csMulberry32(42);
const a1 = []; const a2 = [];
for (let i = 0; i < 5; i++) { a1.push(r1()); a2.push(r2()); }
if (JSON.stringify(a1) !== JSON.stringify(a2)) fail('Mulberry32 not deterministic');
ok('Mulberry32 deterministic');
const r3 = window._csMulberry32(43);
const a3 = []; for (let i = 0; i < 5; i++) a3.push(r3());
if (JSON.stringify(a1) === JSON.stringify(a3)) fail('Mulberry32 should differ for different seeds');
ok('Mulberry32 differs across seeds');

// FNV hash determinism
if (window._csSeedFromString('LG28') !== window._csSeedFromString('LG28'))
  fail('seed hash not deterministic');
if (window._csSeedFromString('LG28') === window._csSeedFromString('LG07'))
  fail('seed hash should differ for different inputs');
ok('seed hash deterministic + differs across inputs');

// Binomial tail p-value sanity
// For n=10, p=0.5, k=10: prob = (1/2)^10 = 0.000977
const p10 = window._csBinomialTailP(10, 0.5, 10);
if (Math.abs(p10 - 0.0009765625) > 1e-6) fail('binomial p(10/10 successes, p=0.5) wrong: ' + p10);
ok('binomial tail p(10/10, p=0.5) = ' + p10.toFixed(6));
// For n=10, p=0.5, k=0: prob >= 0 = 1.0
const p0 = window._csBinomialTailP(10, 0.5, 0);
if (Math.abs(p0 - 1.0) > 1e-6) fail('binomial p(>=0) should be 1.0, got ' + p0);
ok('binomial tail p(>=0) = 1.0');

// ---- 10. Natural chrom sort ----
const sorted = window._csNaturalSortChroms(['LG10', 'LG02', 'LG28', 'LG01', 'LG15']);
if (sorted.join(',') !== 'LG01,LG02,LG10,LG15,LG28')
  fail('natural sort wrong: ' + sorted.join(','));
ok('natural chrom sort works');

// ---- 11. Filter / step nav (still working from breakpoint layer) ----
const cat = document.getElementById('csCatalogue');
const rows = cat.querySelectorAll('[data-cs-bp-id]');
if (rows.length !== 3) fail('expected 3 catalogue rows, got ' + rows.length);
ok('catalogue has 3 rows');
window._crossSpeciesToggleEventFilter('inversion');
if (cat.querySelectorAll('[data-cs-bp-id]').length !== 2) fail('inversion filter expected 2 rows');
ok('event filter still works');
window._crossSpeciesClearFilters();

// ---- 12. v1 backward compatibility ----
const v1Json = { ...goodJson, schema_version: 1 };
delete v1Json.synteny_blocks;
delete v1Json.n_synteny_blocks;
delete v1Json.chrom_lengths_query;
delete v1Json.chrom_lengths_target;
window._clearCrossSpecies();
window._storeCrossSpecies(v1Json);
window._renderCrossSpeciesPage();
const synSlotV1 = document.getElementById('csSynteny');
if (synSlotV1.style.display !== 'none')
  fail('synteny section should be hidden for v1 JSON');
ok('v1 JSON: synteny section correctly hidden');

// ---- 13. Restore v2 round-trip ----
window._clearCrossSpecies();
window._storeCrossSpecies(goodJson);
window._persistCrossSpecies();
window._clearCrossSpecies();   // wipes state and LS
// Re-set LS manually
window.localStorage.setItem('pca_scrubber_v3.crossSpecies.v1', JSON.stringify({
  ...goodJson, loaded_at: '2026-04-30T00:00:00Z',
}));
if (!window._restoreCrossSpecies()) fail('_restoreCrossSpecies failed for v2 LS');
window._renderCrossSpeciesPage();
const restoredSyn = window._csComputeSynteny();
if (!restoredSyn) fail('synteny missing after restore');
if (restoredSyn.flows.length !== 4) fail('expected 4 flows after restore, got ' + restoredSyn.flows.length);
ok('v2 restore round-trip preserves synteny_blocks');

// ---- 14. Help table ----
const helpPage = document.getElementById('page5');
if (!helpPage.innerHTML.includes('cross-species')) fail('help page missing cross-species row');
if (!helpPage.innerHTML.includes('15 marker panel')) fail('help table missing marker panel row (renumbered to 15)');
if (!helpPage.innerHTML.includes('16 help')) fail('help table not renumbered to "16 help"');
if (!helpPage.innerHTML.includes('14 stats profile')) fail('help table missing stats profile row (renumbered to 14)');
ok('help table updated (with stats profile + marker panel)');

console.log('\n[smoke] ALL CHECKS PASSED');
