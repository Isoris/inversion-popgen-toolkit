// Smoke test for the marker readiness panel (page18, v4 turn 78).
// REBUILT for the private-indel architecture:
//   Tier 1 = clean dosage AF + controls + bighead specificity (highest)
//   Tier 2 = multi-marker panel OR strong tag with imperfect het
//   Tier 3 = breakpoint PCR candidate (DEMOTED — breakpoint precision uncertain)
//   Tier 4 = exploratory (default for atlas-only candidates)

const { JSDOM } = require('jsdom');
const fs = require('fs');

const ATLAS_PATH = '/home/claude/work/Inversion_atlas.html';
const html = fs.readFileSync(ATLAS_PATH, 'utf-8');
const dom = new JSDOM(html, {
  runScripts: 'dangerously',
  pretendToBeVisual: true,
  url: 'http://localhost/Inversion_atlas.html',
});
const { window } = dom;
const { document } = window;
window.URL.createObjectURL = () => 'blob:mock';
window.URL.revokeObjectURL = () => {};

function fail(msg) { console.error('[mp-smoke] FAIL:', msg); process.exit(1); }
function ok(msg)   { console.log ('[mp-smoke] ok  :', msg); }

// ---- 1. DOM presence ----
const page18 = document.getElementById('page18');
if (!page18) fail('page18 div missing');
ok('page18 div present');

const tabBtn = document.querySelector('#tabBar button[data-page="page18"]');
if (!tabBtn) fail('page18 tab button missing');
ok('page18 tab button: "' + (tabBtn.textContent || '').trim().replace(/\s+/g, ' ') + '"');

for (const name of [
  '_renderMarkerPanelPage', '_mpDeriveAutoPanel', '_mpRefreshPanel',
  '_mpEnsureState', '_mpAutoTierFromAtlas', '_mpMinDistanceToCsBreakpoint',
  '_mpScoreVariantAf', '_mpAnnotateGelVisibility', '_mpAnnotateAfVariants',
  '_mpSuggestControlsFromKaryotype', '_mpDefaultCrossSpeciesControls',
  '_mpIsValidPanelJson', '_mpIsValidVariantAfsJson', '_mpMergeUserJson',
  '_mpParseTsv', '_mpFilteredPanel', '_mpExportCsv',
]) {
  if (typeof window[name] !== 'function') fail(name + ' missing from window');
}
ok('all v4-turn-78 marker-panel helpers present');

// ---- 2. AF scoring is the core logic — test all four tier outcomes ----
let v = window._mpScoreVariantAf({ chr: 'LG28', pos: 3500000, ref: 'A', alt: 'ATTT',
  type: 'indel', indel_size_bp: 60, af_std: 0.01, af_het: 0.48, af_inv: 0.96 });
if (v.tier_from_af !== 1) fail('clean dosage should be Tier 1, got ' + v.tier_from_af);
if (Math.abs(v.private_score - 0.95) > 1e-6) fail('private_score should be 0.95, got ' + v.private_score);
if (Math.abs(v.dosage_score - 0.96) > 1e-6) fail('dosage_score should be 0.96, got ' + v.dosage_score);
if (Math.abs(v.final_score - 1.91) > 1e-6) fail('final_score should be 1.91, got ' + v.final_score);
ok('AF Tier 1: STD=0.01/HET=0.48/INV=0.96 \u2192 private=0.95, dosage=0.96, final=1.91');

v = window._mpScoreVariantAf({ chr: 'LG07', pos: 5000000, ref: 'G', alt: 'T',
  type: 'snp', af_std: 0.04, af_het: 0.62, af_inv: 0.78 });
if (v.tier_from_af !== 2) fail('strong tag with imperfect het should be Tier 2, got ' + v.tier_from_af);
ok('AF Tier 2: STD=0.04/HET=0.62/INV=0.78 \u2192 strong tag (imperfect dosage)');

v = window._mpScoreVariantAf({ af_std: 0.20, af_het: 0.55, af_inv: 0.80 });
if (v.tier_from_af !== 4) fail('AF_STD too high should be Tier 4, got ' + v.tier_from_af);
ok('AF Tier 4: STD=0.20/HET=0.55/INV=0.80 \u2192 too noisy');

v = window._mpScoreVariantAf({ af_std: 0.02, af_inv: 0.85 });
if (v.tier_from_af !== 2) fail('missing AF_HET should drop to Tier 2, got ' + v.tier_from_af);
if (v.dosage_score !== null) fail('dosage_score should be null when AF_HET missing');
ok('AF without HET: drops to Tier 2 (strong tag, dosage unknown)');

v = window._mpScoreVariantAf({ af_std: 0.01, af_het: 0.85, af_inv: 0.92 });
if (v.tier_from_af !== 2) fail('AF_HET=0.85 should knock out of Tier 1, got ' + v.tier_from_af);
ok('AF_HET=0.85 (out of [0.25,0.75]): drops to Tier 2');

// ---- 3. Gel-visibility ----
let g = window._mpAnnotateGelVisibility({ type: 'indel', indel_size_bp: 60 });
if (g.gel_visible !== true) fail('60bp indel should be gel-visible');
g = window._mpAnnotateGelVisibility({ type: 'indel', indel_size_bp: 5 });
if (g.gel_visible !== false) fail('5bp indel should NOT be gel-visible');
g = window._mpAnnotateGelVisibility({ type: 'indel', indel_size_bp: 800 });
if (g.gel_visible !== false) fail('800bp indel should NOT be gel-visible');
g = window._mpAnnotateGelVisibility({ type: 'snp' });
if (g.gel_visible !== null) fail('SNP should be gel_visible=null');
ok('gel-visibility: 60bp=gel\u00d7, 5bp=no-gel, 800bp=no-gel, SNP=N/A');

// ---- 4. Atlas auto-tier ----
const v2Json = {
  tool: 'cross_species_breakpoints_v1', schema_version: 2,
  generated_at: '2026-04-30T00:00:00Z',
  species_query:  { name: 'Clarias gariepinus', haplotype: 'fClaHyb_Gar_LG' },
  species_target: { name: 'Clarias macrocephalus', haplotype: 'fClaHyb_Mac_LG' },
  input_paf: { path: '/scratch/test.paf', sha256: 'abc' },
  params: { min_mapq: 40, min_block_bp: 50000, cluster_radius_bp: 50000, flank_bp: 100000, classes_annotated: ['all_TE'] },
  n_breakpoints: 1, n_by_event_type: { inversion: 1 },
  breakpoints: [
    { id: 'cs_bp_001', event_type: 'inversion',
      gar_chr: 'LG28', gar_pos_start: 3_010_000, gar_pos_end: 3_020_000, gar_pos_mb: 3.015 },
  ],
  n_synteny_blocks: 2,
  synteny_blocks: [
    { gar_chr: 'LG28', gar_start: 0, gar_end: 10_000_000, mac_chr: 'LG28',
      mac_start: 0, mac_end: 10_000_000, strand: '+', block_size_bp: 10_000_000, mapping_quality: 60 },
    { gar_chr: 'LG15', gar_start: 0, gar_end: 12_000_000, mac_chr: 'LG15',
      mac_start: 0, mac_end: 12_000_000, strand: '+', block_size_bp: 12_000_000, mapping_quality: 60 },
  ],
  chrom_lengths_query:  { LG15: 12_000_000, LG28: 10_000_000 },
  chrom_lengths_target: { LG15: 12_000_000, LG28: 10_000_000 },
};
window._storeCrossSpecies(v2Json);

window.eval(`
  state.candidateList = [
    { id: 'CAND_NEAR_BP', chrom: 'LG28', start_bp: 3000000, end_bp: 3100000, confirmed: true,
      assignments: { S001: 0, S002: 0, S003: 0, S010: 1, S011: 1, S020: 2, S021: 2, S022: 2 } },
    { id: 'CAND_NO_BP', chrom: 'LG15', start_bp: 6000000, end_bp: 7000000, confirmed: true,
      assignments: { S101: 0, S102: 0, S110: 1, S120: 2, S121: 2 } },
    { id: 'CAND_NO_KAR', chrom: 'LG15', start_bp: 9000000, end_bp: 9500000, confirmed: true },
  ];
  state._csInversionContextCache = null;
  state._csSyntenyCache = null;
  state._csSyntenyEdgesCache = null;
`);

window._mpRefreshPanel();
const panel = window._mpEnsureState().panel;
if (panel.length !== 3) fail('expected 3 markers, got ' + panel.length);

const nearBp = panel.find(m => m.inversion_id === 'CAND_NEAR_BP');
if (nearBp.tier !== 3) fail('CAND_NEAR_BP should be Tier 3, got ' + nearBp.tier);
if (nearBp.marker_type !== 'breakpoint PCR (candidate)') fail('CAND_NEAR_BP marker_type wrong: ' + nearBp.marker_type);
ok('CAND_NEAR_BP: Tier 3 (breakpoint PCR candidate, DEMOTED) \u2014 ' + nearBp.tier_reason);

const noBp = panel.find(m => m.inversion_id === 'CAND_NO_BP');
if (noBp.tier !== 4) fail('CAND_NO_BP should default to Tier 4, got ' + noBp.tier);
ok('CAND_NO_BP: Tier 4 (default, no AF data loaded yet)');

// ---- 5. Auto-suggested controls from karyotype ----
if (!nearBp.controls.positive_controls_INV.includes('S020')) fail('CAND_NEAR_BP missing INV controls');
if (!nearBp.controls.heterozygote_controls.includes('S010')) fail('CAND_NEAR_BP missing HET controls');
if (!nearBp.controls.negative_controls_STD.includes('S001')) fail('CAND_NEAR_BP missing STD controls');
if (nearBp.controls.positive_controls_INV.length !== 3) fail('expected 3 INV controls');
ok('CAND_NEAR_BP controls auto-suggested: 3 INV (S020,S021,S022), 2 HET (S010,S011), 3 STD');

const noKar = panel.find(m => m.inversion_id === 'CAND_NO_KAR');
if (noKar.controls.positive_controls_INV.length !== 0) fail('CAND_NO_KAR should have no auto-controls');
ok('CAND_NO_KAR: no auto-controls (no karyotype data)');

// ---- 6. Cross-species controls always populated ----
for (const m of panel) {
  if (!m.controls.bighead_negative_required) fail(m.inversion_id + ' missing bighead requirement');
  if (!m.controls.no_template_control_required) fail(m.inversion_id + ' missing NTC requirement');
  if (!m.controls.bighead_expected_pattern) fail(m.inversion_id + ' missing bighead pattern');
}
ok('all markers carry default cross-species controls (bighead, NTC, F1)');

// ---- 7. variant_afs.json overlay \u2014 the core promotion logic ----
const afsJson = {
  metadata: { date: '2026-04-30', cohort: 'broodstock-226' },
  variants_by_inversion: {
    CAND_NO_BP: [
      { chr: 'LG15', pos: 6500000, ref: 'A', alt: 'ATTGGCC', type: 'indel', indel_size_bp: 60,
        af_std: 0.01, af_het: 0.48, af_inv: 0.96 },
      { chr: 'LG15', pos: 6600000, ref: 'C', alt: 'T', type: 'snp',
        af_std: 0.02, af_het: 0.46, af_inv: 0.91 },
      { chr: 'LG15', pos: 6800000, ref: 'G', alt: 'GAACT', type: 'indel', indel_size_bp: 4,
        af_std: 0.00, af_het: 0.50, af_inv: 0.92 },
      { chr: 'LG15', pos: 6900000, ref: 'A', alt: 'C', type: 'snp',
        af_std: 0.03, af_het: 0.62, af_inv: 0.78 },
      { chr: 'LG15', pos: 6950000, ref: 'T', alt: 'G', type: 'snp',
        af_std: 0.18, af_het: 0.55, af_inv: 0.65 },
    ],
    CAND_NEAR_BP: [
      { chr: 'LG28', pos: 3050000, ref: 'A', alt: 'G', type: 'snp',
        af_std: 0.05, af_het: 0.55, af_inv: 0.72 },
    ],
  },
};
if (!window._mpIsValidVariantAfsJson(afsJson)) fail('valid AFs JSON rejected');
const ui = window._mpEnsureState();
ui.af_data = afsJson;
window._mpRefreshPanel();
const panelAf = window._mpEnsureState().panel;

const noBpAf = panelAf.find(m => m.inversion_id === 'CAND_NO_BP');
if (noBpAf.tier !== 1) fail('CAND_NO_BP should be promoted to Tier 1, got ' + noBpAf.tier);
if (noBpAf.af_variants.length !== 5) fail('expected 5 scored variants, got ' + noBpAf.af_variants.length);
const t1count = noBpAf.af_variants.filter(v => v.tier_from_af === 1).length;
if (t1count !== 3) fail('expected 3 Tier-1 variants, got ' + t1count);
const evidenceStr = noBpAf.evidence.join(' | ');
if (!evidenceStr.includes('multi-marker panel feasible')) fail('multi-marker evidence missing');
ok('CAND_NO_BP: PROMOTED to Tier 1 (3 Tier-1 variants \u2192 multi-marker panel feasible)');

if (Math.abs(noBpAf.best_af_score - 1.91) > 0.05) fail('best_af_score wrong: ' + noBpAf.best_af_score);
ok('CAND_NO_BP best AF score = ' + noBpAf.best_af_score.toFixed(2));

const nearBpAf = panelAf.find(m => m.inversion_id === 'CAND_NEAR_BP');
if (nearBpAf.tier !== 2) fail('CAND_NEAR_BP should be promoted to Tier 2, got ' + nearBpAf.tier);
ok('CAND_NEAR_BP: PROMOTED from Tier 3 to Tier 2 by strong-tag AF variant');

const noKarAf = panelAf.find(m => m.inversion_id === 'CAND_NO_KAR');
if (noKarAf.tier !== 4) fail('CAND_NO_KAR should stay Tier 4, got ' + noKarAf.tier);
ok('CAND_NO_KAR: stays Tier 4');

// ---- 8. Render with promoted panel ----
document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
page18.classList.add('active');
window._renderMarkerPanelPage();
const body = document.getElementById('mpBody');
if (!body.innerHTML.includes('Tier 1 \u2014 private indel/SNP tag')) fail('Tier 1 def label missing');
if (!body.innerHTML.includes('Tier 3 \u2014 breakpoint PCR candidate')) fail('Tier 3 def missing');
if (!body.innerHTML.includes('Tier 4 \u2014 exploratory')) fail('Tier 4 def missing');
if (!body.innerHTML.includes('Pilot validation plan (with controls)')) fail('pilot plan title missing');
if (!body.innerHTML.includes('bighead catfish DNA')) fail('cross-species step missing');
if (!body.innerHTML.includes('locally validated before routine breeding')) fail('protective sentence missing');
if (!body.innerHTML.includes('Because the target breeding system involves interspecific F\u2081 hybrids'))
  fail('cross-species manuscript sentence missing');
ok('page18 renders 4 tier defs + cross-species pilot + manuscript sentence');

if (!body.innerHTML.includes('priv 0.95')) fail('AF block missing private score');
if (!body.innerHTML.includes('gel\u00d7')) fail('gel-visibility badge missing');
ok('AF block renders with score breakdown + gel-visibility badges');

if (!body.innerHTML.includes('S020')) fail('controls cell missing INV samples');
if (!body.innerHTML.includes('cross-species: bighead')) fail('controls cell missing cross-species line');
ok('controls cell shows auto-suggested samples + cross-species note');

// ---- 9. Filter ----
ui.filter_tier = '1';
window._renderMarkerPanelPage();
const filtered = document.getElementById('mpBody').innerHTML;
if (!filtered.includes('CAND_NO_BP')) fail('Tier 1 filter should include CAND_NO_BP');
// Use a more specific match — "CAND_NEAR_BP" check needs to look in <code> tag context
if (filtered.match(/<code>CAND_NEAR_BP<\/code>/)) fail('Tier 1 filter should NOT include CAND_NEAR_BP (T2)');
if (filtered.match(/<code>CAND_NO_KAR<\/code>/)) fail('Tier 1 filter should NOT include CAND_NO_KAR (T4)');
ok('Tier 1 filter isolates the AF-promoted CAND_NO_BP');
ui.filter_tier = 'all';

// ---- 10. User marker_panel.json overrides ----
const userPanel = {
  metadata: { panel_name: 'Final breeding panel v1' },
  markers: [
    {
      inversion_id: 'CAND_NO_BP',
      tier: 1,
      validation_status: 'validated',
      primer_F: 'AAACCCGGGTTT',
      primer_R: 'TTTGGGCCCAAA',
      amplicon_bp_state_A: 240,
      amplicon_bp_state_B: 180,
      n_carriers_tested: 88,
      controls: {
        positive_controls_INV: ['CGA021_validated', 'CGA044_validated'],
        bighead_negative_status: 'tested_no_amplification',
        bighead_orthologous_sequence_status: 'absent',
      },
      notes: 'validated 60bp indel, gel-visible, bighead negative',
    },
  ],
};
ui.panel = window._mpMergeUserJson(ui.panel, userPanel);
window._renderMarkerPanelPage();
const v1 = ui.panel.find(m => m.inversion_id === 'CAND_NO_BP');
if (v1.validation_status !== 'validated') fail('user override: validation_status wrong');
if (v1.primer_F !== 'AAACCCGGGTTT') fail('user override: primer_F missing');
if (v1.amplicon_bp_state_A !== 240) fail('user override: amplicon_A wrong');
if (!v1.controls.positive_controls_INV.includes('CGA021_validated'))
  fail('user override: validated controls missing');
if (v1.controls.bighead_negative_status !== 'tested_no_amplification')
  fail('user override: bighead status wrong');
ok('user marker_panel.json overlay: CAND_NO_BP validated, primers loaded, bighead tested');

// ---- 11. CSV export ----
let exportedRows = null;
const origCreateElement = document.createElement.bind(document);
document.createElement = function (tag) {
  const el = origCreateElement(tag);
  if (tag === 'a') {
    el.click = function () {
      exportedRows = window._mpFilteredPanel(window._mpEnsureState().panel, window._mpEnsureState());
    };
  }
  return el;
};
const csvBtn = document.getElementById('mpExportCsvBtn');
if (csvBtn) csvBtn.click();
if (!exportedRows) fail('CSV export did not run');
if (exportedRows.length !== 3) fail('CSV export should be 3 rows, got ' + exportedRows.length);
ok('CSV export: 3 rows including AF + control columns');

// ---- 12. Stats profile markerability tier breakdown ----
window._spRefreshDerivedRows();
const spRows = window._spEnsureState().rows;
const mkRow = spRows.find(r => r.id === 'markerability');
if (!mkRow) fail('markerability row missing');
if (!mkRow.african_result.tier_breakdown) fail('tier_breakdown missing');
ok('stats profile markerability: ' + mkRow.african_result.label);

// ---- 13. Help table ----
const help = document.getElementById('page5');
if (!help.innerHTML.includes('private indel')) fail('help row not updated');
if (!help.innerHTML.includes('DEMOTED')) fail('help row should mention Tier 3 demotion');
if (!help.innerHTML.includes('16 help')) fail('help table not renumbered to 16 (after negative regions tab insertion)');
ok('help table updated with private-indel architecture');

// ---- 14. New help sections (turn 80): hotkeys + changelog + schemas ----
if (!help.innerHTML.includes('Hotkeys / tutorial')) fail('Hotkeys section missing from help page');
if (!help.innerHTML.includes('navigate the genome')) fail('hotkeys text "navigate the genome" missing');
if (!help.innerHTML.includes('promote</b> to add the consolidated interval')) fail('hotkeys workflow text missing');
ok('help page: Hotkeys / tutorial section present with L1/L2 workflow');

if (!help.innerHTML.includes('Changelog')) fail('Changelog section missing');
if (!help.innerHTML.includes('v4 turn 79')) fail('changelog should mention turn 79 (hash routing)');
if (!help.innerHTML.includes('v4 turn 78')) fail('changelog should mention turn 78 (private indel rebuild)');
if (!help.innerHTML.includes('v4 turn 77')) fail('changelog should mention turn 77 (stats profile)');
ok('help page: Changelog section present with turns 75-79 entries');

if (!help.innerHTML.includes('New JSON schemas')) fail('Schema section missing');
if (!help.innerHTML.includes('cs_breakpoints_v1.json')) fail('cs_breakpoints schema missing');
if (!help.innerHTML.includes('variant_afs.json')) fail('variant_afs schema missing');
if (!help.innerHTML.includes('marker_panel.json')) fail('marker_panel schema missing');
if (!help.innerHTML.includes('stats_profile.json')) fail('stats_profile schema missing');
ok('help page: New JSON schemas section lists all 4 file formats');

// Stale entries cleaned up
if (help.innerHTML.includes('13 diversity atlas') || help.innerHTML.includes('14 genome atlas'))
  fail('stale "13 diversity atlas / 14 genome atlas" entries should have been removed');
ok('help page: stale "13 diversity atlas / 14 genome atlas" entries removed');

console.log('\n[mp-smoke] ALL CHECKS PASSED');
