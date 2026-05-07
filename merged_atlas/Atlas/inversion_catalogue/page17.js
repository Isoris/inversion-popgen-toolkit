// Atlas/inversion_catalogue/page17.js
// =============================================================================
// page17 — Statistical profile of inversion-associated genomic features
// (`<div id="page17">` contains `<div id="spBody">`)
//
// Source: legacy/Inversion_atlas.html lines 28420–29306 (JS body)
//                                    lines  8110–8127  (HTML shell)
// Public entry: renderStatsProfilePage()  (legacy: _renderStatsProfilePage line 29146)
//
// External dependencies — marked TODO_MISSING. These are referenced by the
// extracted body but not defined in this module. The merge chat resolves them
// (either promotes to shared/ or co-locates with the page that owns them):
//
//   TODO_MISSING(_csGetSyntenyBlocks) — likely legacy line ~? (cross-species
//                                       synteny lookup, owned by Batch 5
//                                       inversion_comparative)
//   TODO_MISSING(_csPermutationTest)  — likely legacy line ~? (cross-species
//                                       permutation test, owned by Batch 5
//                                       inversion_comparative)
//   TODO_MISSING(_esc)                — HTML-escape helper. High frequency
//                                       across many pages — strong candidate
//                                       for promotion to shared/dom_utils.js.
//   TODO_MISSING(_mpDeriveAutoPanel)  — exported by sibling page18.js in this
//                                       same batch; merge wires the import.
//
//   global `state`                    — atlas state singleton; reads:
//                                         state.crossSpecies
//                                         state.candidateList
//                                         state.data
//                                         state._statsProfile
//                                         state._markerPanel
//                                       The merge chat refactors to ctx-pattern
//                                       per shared/per_l2_cluster.js.
// =============================================================================

// Reuse global `state` until the merge chat refactors to contextFromState.
const state = (typeof window !== 'undefined' && window.state) ? window.state : {};

// ---------------------------------------------------------------------------
// VERBATIM extraction from legacy lines 28445–29306. Only modification: the
// outermost legacy-global function declarations remain `function foo()` so
// hoisting still works inside this module; we then export the public entry
// at the bottom.
// ---------------------------------------------------------------------------

// Effect-tag colors (subtle, professional). Tied to data-effect attrs in CSS.
const SP_EFFECT_COLORS = {
  'enriched':              '#3cc08a',
  'depleted':              '#f0a35e',
  'closer_than_expected':  '#3cc08a',
  'farther_than_expected': '#f0a35e',
  'higher_than_background':'#3cc08a',
  'lower_than_background': '#f0a35e',
  'common':                '#3cc08a',
  'rare':                  '#f0a35e',
  'not_significant':       '#9aa3ad',
  'not_tested':            '#6e7782',
  'african_only':          '#7ad3db',
  'cross_species':         '#b07cf7',
  // Orthogonal flag for "this row is supported by atlas data right now"
  'derived':               '#4fa3ff',
  'loaded':                '#3cc08a',
  'placeholder':           '#9aa3ad',
};

// Significance threshold for the "n significant" summary card. Same as the
// rest of the atlas's defaults.
const SP_ALPHA = 0.05;

// Default schema for the stats-profile table. Each row has a stable `id`
// (used to merge user-supplied overrides), a `category`, and a description
// of how to derive each species' result. The 13 rows mirror the screenshots
// from the user's "comparative inversion phenotype page" mockup.
const SP_DEFAULT_ROWS = [
  // ------ A. Breakpoint architecture ------------------------------------------
  {
    id: 'breakpoints_near_synteny_edges',
    category: 'Breakpoint architecture',
    statistic: 'Breakpoints near synteny edges',
    test: 'matched permutation (length-matched)',
    null_comparison: 'random intervals on same chrom, length matched',
    derive_from: 'cs_permutation_test',
    african_hint: 'derived from page 16 cross-species permutation test',
    bighead_hint: 'orthologous context — same synteny blocks',
    interpretation_default: 'Whether African inversion breakpoints preferentially occur near macro-syntenic transition zones.',
  },
  // ------ B. Fusion/fission context ------------------------------------------
  {
    id: 'fusion_fission_proximity',
    category: 'Breakpoint architecture',
    statistic: 'Proximity to fusion/fission boundaries',
    test: 'Fisher / matched permutation',
    null_comparison: 'random matched breakpoint positions',
    derive_from: 'cs_breakpoints_event_type',
    african_hint: 'count of inversions on chroms classified as fission_in_target or many_to_many (page 16)',
    bighead_hint: 'symmetric',
    interpretation_default: 'Subset of inversions may reuse deeper chromosome-evolution boundaries.',
  },
  // ------ C. Repeat context around breakpoints --------------------------------
  {
    id: 'repeat_density_flanks',
    category: 'Genomic composition',
    statistic: 'Repeat (all_TE) density flanking breakpoints',
    test: 'Wilcoxon rank-sum',
    null_comparison: 'matched random breakpoint windows',
    derive_from: 'cs_breakpoints_flank_density',
    african_hint: 'mean ±100kb all_TE on Gar haplotype (page 16 flanks)',
    bighead_hint: 'mean ±100kb all_TE on Mac haplotype (page 16 flanks)',
    interpretation_default: 'Repeats may contribute to local structural instability at breakpoints.',
  },
  // ------ D. Telomere / centromere / gap proximity ----------------------------
  {
    id: 'telomere_centromere_proximity',
    category: 'Genomic composition',
    statistic: 'Distance to telomere / centromere / scaffold gap',
    test: 'permutation / threshold enrichment',
    null_comparison: 'random matched intervals',
    derive_from: null,
    african_hint: 'requires telomere + centromere annotation track',
    bighead_hint: 'requires Mac telomere + centromere annotation track',
    interpretation_default: 'Breakpoint regions may occur near difficult or fragile chromosome architecture.',
  },
  // ------ E. Gene density inside inversions -----------------------------------
  {
    id: 'gene_density_inversions',
    category: 'Genomic composition',
    statistic: 'Gene density inside inversions',
    test: 'Wilcoxon / matched interval permutation',
    null_comparison: 'size-matched random intervals on same chromosomes',
    derive_from: null,
    african_hint: 'requires GFF/BED gene annotation',
    bighead_hint: 'requires Mac GFF/BED gene annotation',
    interpretation_default: 'Inversions capture gene-rich or gene-poor domains non-randomly.',
  },
  // ------ F. Functional cargo (GO/KEGG) ---------------------------------------
  {
    id: 'go_kegg_enrichment',
    category: 'Functional cargo',
    statistic: 'GO / KEGG / Pfam enrichment among inversion-internal genes',
    test: 'hypergeometric + BH correction',
    null_comparison: 'all annotated genes (background)',
    derive_from: null,
    african_hint: 'requires gene-set enrichment table (e.g. clusterProfiler output)',
    bighead_hint: 'requires Mac gene-set enrichment table',
    interpretation_default: 'Inversion cargo may be enriched for immune, stress, growth, reproduction, or metabolism functions.',
  },
  // ------ G. Immune/growth/stress/reproduction candidates ---------------------
  {
    id: 'breeding_candidate_genes',
    category: 'Functional cargo',
    statistic: 'Breeding-candidate genes inside inversions',
    test: 'Fisher / hypergeometric',
    null_comparison: 'matched random gene sets',
    derive_from: null,
    african_hint: 'requires curated gene list (immune/growth/stress/reproduction)',
    bighead_hint: 'requires curated gene list',
    interpretation_default: 'Possible breeding-relevant cargo. Treat as candidate, not causal.',
  },
  // ------ H. Heterozygosity inside inversions ---------------------------------
  {
    id: 'heterozygosity_inversions',
    category: 'Population variation',
    statistic: 'Heterozygosity / nucleotide diversity inside inversions',
    test: 'Wilcoxon / permutation',
    null_comparison: 'matched genomic intervals',
    derive_from: null,
    african_hint: 'requires per-window het from page 1 (state.data.diversity_tracks)',
    bighead_hint: 'African-only unless Mac population data available',
    african_only: true,
    interpretation_default: 'Inversions preserve divergent haplotype backgrounds or reduced diversity.',
  },
  // ------ I. ROH / inbreeding overlap -----------------------------------------
  {
    id: 'roh_overlap',
    category: 'Population variation',
    statistic: 'ROH / FROH overlap with inversions',
    test: 'permutation / Fisher / Wilcoxon',
    null_comparison: 'matched intervals',
    derive_from: null,
    african_hint: 'requires ROH BED from MODULE_3_ROH (ngsF-HMM)',
    bighead_hint: 'African-only',
    african_only: true,
    interpretation_default: 'Inversion haplotypes may overlap founder/inbreeding blocks.',
  },
  // ------ J. Deleterious burden -----------------------------------------------
  {
    id: 'deleterious_burden',
    category: 'Population variation',
    statistic: 'Deleterious-variant burden in inversion haplotypes',
    test: 'Wilcoxon / GLM (family-corrected)',
    null_comparison: 'carriers vs non-carriers, matched intervals',
    derive_from: null,
    african_hint: 'requires SIFT4G / VESM_650M output (MODULE_CONSERVATION)',
    bighead_hint: 'African-only',
    african_only: true,
    interpretation_default: 'Some structural haplotypes may carry different breeding risk.',
  },
  // ------ K. FST / dXY between inversion haplotypes ---------------------------
  {
    id: 'haplotype_differentiation',
    category: 'Population variation',
    statistic: 'FST / dXY between inversion structural states',
    test: 'permutation / bootstrap window comparison',
    null_comparison: 'matched intervals or randomized labels',
    derive_from: 'candidate_list_dosage_states',
    african_hint: 'derived from candidate list dosage groups (page 1) + region_popstats',
    bighead_hint: 'African-only',
    african_only: true,
    interpretation_default: 'Inversion states represent differentiated structural haplotypes.',
  },
  // ------ L. Markerability ----------------------------------------------------
  {
    id: 'markerability',
    category: 'Breeding utility',
    statistic: 'High-confidence inversions with breakpoint / tagging markers',
    test: 'descriptive (no p-value)',
    null_comparison: 'NA',
    derive_from: 'candidate_list_count',
    african_hint: 'count of confirmed candidates from page 3 catalogue',
    bighead_hint: 'NA',
    interpretation_default: 'Practical conversion of atlas candidates into marker-assisted breeding tools.',
  },
];

// =============================================================================
// Atlas-derivation helpers
// -----------------------------------------------------------------------------

// Run the cross-species permutation test (page 16) and convert the result
// into a `result` object usable by both the African and Bighead columns.
// Returns null if cs_breakpoints v2 + candidate list aren't both loaded.
function _spDeriveCsPermutation() {
  if (typeof _csPermutationTest !== 'function') return null;
  const blocks = (typeof _csGetSyntenyBlocks === 'function') ? _csGetSyntenyBlocks() : null;
  if (!blocks) return null;
  if (!Array.isArray(state.candidateList) || state.candidateList.length === 0) return null;
  let perm;
  try { perm = _csPermutationTest({ n_perm: 2000 }); }
  catch (e) { console.warn('[stats_profile] perm test failed:', e); return null; }
  if (!perm || !perm.global || perm.global.n_valid === 0) return null;
  const g = perm.global;
  const enr = g.enrichment_ratio;
  const p   = g.binomial_p;
  // Effect tag based on direction + significance
  let effect = 'not_significant';
  if (p != null && p < SP_ALPHA) {
    effect = (enr != null && enr > 1) ? 'closer_than_expected' : 'farther_than_expected';
  }
  const label = (g.near_edge_observed + ' / ' + g.n_valid +
    ' candidates near synteny edges (q=' + (g.quantile * 100).toFixed(0) + '%); ' +
    'enrichment ' + (enr != null ? enr.toFixed(2) : '\u2014') + ', ' +
    'binomial p=' + (p != null ? (p < 0.001 ? p.toExponential(2) : p.toFixed(3)) : '\u2014'));
  return {
    state: 'derived',
    effect,
    observed: g.near_edge_observed,
    expected: g.near_edge_expected,
    unit: 'candidates',
    fold_change: enr,
    p_value: p,
    q_value: null,
    n: g.n_valid,
    label,
  };
}

// Derive flanking-repeat-density Wilcoxon-ish proxy from the cs_breakpoints
// JSON. We don't actually run the Wilcoxon (no chrom-mean per-window
// distribution in browser) — instead we report the fraction of breakpoints
// flagged with the Spalax-style manuscript_note (flank all_TE >= 1.5x
// chrom mean, computed by STEP_CS01_extract_breakpoints.py).
function _spDeriveRepeatFlankSpalax(species) {
  if (!state.crossSpecies || !Array.isArray(state.crossSpecies.breakpoints)) return null;
  const bps = state.crossSpecies.breakpoints;
  if (bps.length === 0) return null;
  // The Spalax manuscript_note is computed on the GAR side — for the Mac
  // (target/bighead) we don't have a comparable chrom-mean baseline yet,
  // so report the simple mean flank density as a placeholder summary.
  if (species === 'african') {
    const flagged = bps.filter(b => b.manuscript_note && b.manuscript_note.indexOf('Spalax') >= 0).length;
    const frac = flagged / bps.length;
    let effect = 'not_significant';
    if (frac >= 0.20) effect = 'higher_than_background';
    else if (frac > 0.0 && frac < 0.05) effect = 'lower_than_background';
    return {
      state: 'derived',
      effect,
      observed: flagged,
      expected: null,
      unit: 'breakpoints',
      fold_change: null,
      p_value: null,
      q_value: null,
      n: bps.length,
      label: flagged + ' / ' + bps.length + ' breakpoints flagged Spalax-style ' +
             '(\u22651.5\u00D7 chrom-mean all_TE in \u00B1100kb flank)',
    };
  }
  // Mac/bighead side — report mean ±100kb flank density without a baseline.
  // The `flanking_repeat_density_mac` field has prev/next dicts with
  // by_class.all_TE.mean per side. We average the two and report.
  let nWith = 0, sum = 0;
  for (const b of bps) {
    const m = b.flanking_repeat_density_mac;
    if (!m || !m.prev || !m.next) continue;
    const pv = m.prev.by_class && m.prev.by_class.all_TE && m.prev.by_class.all_TE.mean;
    const nv = m.next.by_class && m.next.by_class.all_TE && m.next.by_class.all_TE.mean;
    if (pv == null && nv == null) continue;
    const vals = [pv, nv].filter(v => v != null);
    sum += vals.reduce((a, b) => a + b, 0) / vals.length;
    nWith += 1;
  }
  if (nWith === 0) return null;
  const mean = sum / nWith;
  return {
    state: 'derived',
    effect: 'not_tested',
    observed: mean,
    expected: null,
    unit: 'all_TE fraction',
    fold_change: null,
    p_value: null,
    q_value: null,
    n: nWith,
    label: 'mean Mac-side flank all_TE = ' + mean.toFixed(3) + ' (n=' + nWith +
           ' breakpoints with both flanks resolved); no chrom-mean baseline yet',
  };
}

// Derive fusion/fission proximity from cs_breakpoints event_type field.
function _spDeriveFusionFission() {
  if (!state.crossSpecies || !Array.isArray(state.crossSpecies.breakpoints)) return null;
  const bps = state.crossSpecies.breakpoints;
  if (bps.length === 0) return null;
  const ff = bps.filter(b => {
    const et = b.event_type_refined || b.event_type || '';
    return et === 'fission_or_fusion' || et === 'translocation_or_fission';
  }).length;
  const inv = bps.filter(b => (b.event_type_refined || b.event_type) === 'inversion').length;
  return {
    state: 'derived',
    effect: ff > 0 ? 'enriched' : 'not_significant',
    observed: ff,
    expected: null,
    unit: 'breakpoints',
    fold_change: null,
    p_value: null,
    q_value: null,
    n: bps.length,
    label: ff + ' / ' + bps.length + ' cross-species breakpoints classified ' +
           'as fission/fusion (' + inv + ' inversion breakpoints)',
  };
}

// Markerability — count promoted candidates from page 3.
// v4 turn 77: markerability derivation upgraded to a tier breakdown.
// Reads the marker panel (auto-tiered from candidate list + cs_breakpoints
// proximity on page 18) and reports the per-tier counts. The label is the
// manuscript-ready breakdown: "T1: N \u00b7 T2: N \u00b7 T3: N (M confirmed)".
function _spDeriveMarkerability() {
  const cands = state.candidateList || [];
  const confirmed = cands.filter(c => c && c.confirmed).length;
  // v4 turn 78: read the LIVE promoted panel (state._markerPanel.panel) when
  // available, so the tier breakdown reflects AF-derived promotions and
  // user-supplied overrides. Fall back to _mpDeriveAutoPanel() when no
  // marker-panel state exists yet, then to a flat count.
  let tierCounts = { 1: 0, 2: 0, 3: 0, 4: 0 };
  let panelAvailable = false;
  let usingLivePanel = false;
  let livePanel = null;
  if (state._markerPanel && Array.isArray(state._markerPanel.panel)) {
    livePanel = state._markerPanel.panel;
    usingLivePanel = true;
  } else if (typeof _mpDeriveAutoPanel === 'function') {
    try { livePanel = _mpDeriveAutoPanel(); }
    catch (e) { console.warn('[stats_profile] tier breakdown failed:', e); }
  }
  if (Array.isArray(livePanel)) {
    for (const m of livePanel) tierCounts[m.tier] = (tierCounts[m.tier] || 0) + 1;
    panelAvailable = livePanel.length > 0;
  }
  const totalAcrossTiers = tierCounts[1] + tierCounts[2] + tierCounts[3] + tierCounts[4];
  const label = panelAvailable
    ? ('Tier 1: ' + tierCounts[1] + ' \u00B7 Tier 2: ' + tierCounts[2] +
       ' \u00B7 Tier 3: ' + tierCounts[3] + ' \u00B7 Tier 4: ' + tierCounts[4] +
       ' (out of ' + totalAcrossTiers + ' markers' +
       (usingLivePanel ? ', live AF-aware promotion' : ', atlas-only auto-tier') +
       '; ' + confirmed + ' confirmed candidates) \u2014 see page 14 marker panel')
    : (confirmed + ' / ' + cands.length + ' confirmed candidates available ' +
       'as marker-assisted breeding targets');
  return {
    state: 'derived',
    effect: 'not_tested',
    observed: tierCounts[1],   // headline metric: count of Tier 1 markers
    expected: null,
    unit: 'Tier 1 markers',
    fold_change: null,
    p_value: null,
    q_value: null,
    n: cands.length,
    tier_breakdown: tierCounts,
    label,
  };
}

// Wire up the derivation pipeline. Returns a fresh array of row objects
// with `african_result` and `bighead_result` populated where possible.
function _spDeriveAllRows() {
  // Defensive deep clone of the schema so live edits don't mutate the const.
  const rows = SP_DEFAULT_ROWS.map(r => Object.assign({}, r));
  for (const r of rows) {
    let af = null, bh = null;
    if (r.derive_from === 'cs_permutation_test') {
      af = _spDeriveCsPermutation();
      // Bighead column is the orthologous context — same permutation test
      // applies symmetrically (the synteny edges are shared). For now we
      // surface "symmetric / orthologous" as a not_tested entry until we
      // have a Mac-cohort candidate list.
      if (af) {
        bh = {
          state: 'derived', effect: 'not_tested',
          observed: null, expected: null, unit: 'candidates',
          fold_change: null, p_value: null, q_value: null, n: null,
          label: 'orthologous context (same synteny edges); no Mac-cohort candidate list yet',
        };
      }
    } else if (r.derive_from === 'cs_breakpoints_event_type') {
      af = _spDeriveFusionFission();
      if (af) bh = { state: 'derived', effect: 'not_tested', label: 'symmetric \u2014 same cross-species breakpoint set', n: af.n };
    } else if (r.derive_from === 'cs_breakpoints_flank_density') {
      af = _spDeriveRepeatFlankSpalax('african');
      bh = _spDeriveRepeatFlankSpalax('bighead');
    } else if (r.derive_from === 'candidate_list_count') {
      af = _spDeriveMarkerability();
    }
    r.african_result = af || _spPlaceholder(r.african_hint);
    r.bighead_result = bh || _spPlaceholder(r.bighead_hint, r.african_only);
    r.interpretation = r.interpretation_default;
  }
  return rows;
}

function _spPlaceholder(hint, africanOnly) {
  return {
    state: 'placeholder',
    effect: africanOnly ? 'african_only' : 'not_tested',
    observed: null, expected: null, unit: null,
    fold_change: null, p_value: null, q_value: null, n: null,
    label: 'add this data \u2014 ' + (hint || 'no atlas-side derivation available'),
  };
}

// =============================================================================
// User-supplied JSON / TSV ingestion (overlays atlas-derived rows)
// -----------------------------------------------------------------------------

function _spIsValidProfileJson(parsed) {
  return parsed
      && typeof parsed === 'object'
      && parsed.metadata
      && Array.isArray(parsed.summary_rows);
}

// Merge user JSON over derived rows. User rows are matched by `id`; if no id
// match, the row is appended. Derived results survive when the user provides
// `null` or omits the result field, so loading a sparse JSON only overrides
// what's in it.
function _spMergeUserJson(derived, parsed) {
  const byId = new Map(derived.map(r => [r.id, r]));
  const meta = parsed.metadata || {};
  const speciesNames = {
    species_1: meta.species_1 || 'C. gariepinus',
    species_2: meta.species_2 || 'C. macrocephalus',
  };
  for (const ur of parsed.summary_rows) {
    const id = ur.id || _spSlug(ur.statistic || ur.category || 'row_' + Math.random());
    let target = byId.get(id);
    if (!target) {
      target = {
        id,
        category:    ur.category    || 'User-supplied',
        statistic:   ur.statistic   || '',
        test:        ur.test        || '',
        null_comparison: ur.null_comparison || '',
        african_result: null,
        bighead_result: null,
        interpretation: '',
      };
      derived.push(target);
      byId.set(id, target);
    }
    if (ur.african_result) {
      target.african_result = Object.assign({}, ur.african_result, { state: 'loaded' });
    }
    if (ur.bighead_result) {
      target.bighead_result = Object.assign({}, ur.bighead_result, { state: 'loaded' });
    }
    if (ur.interpretation) target.interpretation = ur.interpretation;
    if (ur.category)  target.category  = ur.category;
    if (ur.statistic) target.statistic = ur.statistic;
    if (ur.test)      target.test      = ur.test;
    if (ur.null_comparison) target.null_comparison = ur.null_comparison;
  }
  return { rows: derived, species: speciesNames, metadata: meta };
}

function _spSlug(s) {
  return String(s).toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_|_$/g, '');
}

// Parse a TSV with the documented column order. Any column may be empty.
function _spParseTsv(text) {
  const lines = text.split(/\r?\n/).filter(l => l.trim().length > 0);
  if (lines.length < 2) return null;
  const header = lines[0].split('\t').map(s => s.trim());
  const wanted = ['category','statistic','test','null_comparison',
                  'african_label','african_effect','african_observed','african_expected','african_p','african_q',
                  'bighead_label','bighead_effect','bighead_observed','bighead_expected','bighead_p','bighead_q',
                  'interpretation'];
  const idx = {};
  for (const w of wanted) idx[w] = header.indexOf(w);
  const summary_rows = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split('\t');
    const get = k => idx[k] >= 0 ? (parts[idx[k]] || '').trim() : '';
    const num = k => { const v = get(k); return (v === '' || v === 'NA') ? null : Number(v); };
    summary_rows.push({
      id: _spSlug(get('statistic') || get('category')),
      category:  get('category'),
      statistic: get('statistic'),
      test:      get('test'),
      null_comparison: get('null_comparison'),
      african_result: {
        label: get('african_label'),
        effect: get('african_effect') || 'not_tested',
        observed: num('african_observed'),
        expected: num('african_expected'),
        p_value: num('african_p'),
        q_value: num('african_q'),
      },
      bighead_result: {
        label: get('bighead_label'),
        effect: get('bighead_effect') || 'not_tested',
        observed: num('bighead_observed'),
        expected: num('bighead_expected'),
        p_value: num('bighead_p'),
        q_value: num('bighead_q'),
      },
      interpretation: get('interpretation'),
    });
  }
  return { metadata: {}, summary_rows };
}

// =============================================================================
// State container for the page (what's rendered + what filters are active)
// -----------------------------------------------------------------------------
function _spEnsureState() {
  if (!state._statsProfile) {
    state._statsProfile = {
      compact: false,         // detailed vs compact view
      filter_category: 'all',
      filter_significant_only: false,
      filter_scope: 'all',    // 'all' | 'african_only' | 'cross_species'
      species: { species_1: 'C. gariepinus', species_2: 'C. macrocephalus' },
      metadata: { analysis_name: 'Inversion statistical profile' },
      rows: null,
    };
  }
  return state._statsProfile;
}

function _spRefreshDerivedRows() {
  const ui = _spEnsureState();
  ui.rows = _spDeriveAllRows();
}

// =============================================================================
// Rendering
// -----------------------------------------------------------------------------
function _spRenderSummaryCards(rows) {
  const total = rows.length;
  let nSigAf = 0, nSigBh = 0, nBreedingRel = 0, nDerived = 0;
  for (const r of rows) {
    if (r.african_result && r.african_result.p_value != null && r.african_result.p_value < SP_ALPHA) nSigAf += 1;
    if (r.bighead_result && r.bighead_result.p_value != null && r.bighead_result.p_value < SP_ALPHA) nSigBh += 1;
    if (r.category === 'Breeding utility' || r.category === 'Functional cargo' || r.category === 'Population variation') {
      if (r.african_result && r.african_result.state === 'derived') nBreedingRel += 1;
      else if (r.african_result && r.african_result.state === 'loaded') nBreedingRel += 1;
    }
    if (r.african_result && (r.african_result.state === 'derived' || r.african_result.state === 'loaded')) nDerived += 1;
  }
  const card = (label, value, hint) =>
    '<div class="sp-card"><div class="sp-card-value">' + value + '</div>' +
    '<div class="sp-card-label">' + label + '</div>' +
    (hint ? '<div class="sp-card-hint">' + hint + '</div>' : '') + '</div>';
  return '<div class="sp-cards">' +
    card('Statistics tested',     total,        'rows in synthesis table') +
    // v4 turn 111: italicize species abbreviations per scientific convention.
    // Cgar = C. gariepinus, Cmac = C. macrocephalus. Both binomial-style.
    card('Significant in <i>Cgar</i>',   nSigAf,       'p &lt; ' + SP_ALPHA) +
    card('Significant in <i>Cmac</i>',   nSigBh,       'p &lt; ' + SP_ALPHA + ' (orthologous)') +
    card('Atlas-supported rows',  nDerived + ' / ' + total, 'derived live or loaded') +
    '</div>';
}

function _spFilteredRows(rows, ui) {
  return rows.filter(r => {
    if (ui.filter_category !== 'all' && r.category !== ui.filter_category) return false;
    if (ui.filter_significant_only) {
      const af = r.african_result, bh = r.bighead_result;
      const sig = (af && af.p_value != null && af.p_value < SP_ALPHA) ||
                  (bh && bh.p_value != null && bh.p_value < SP_ALPHA);
      if (!sig) return false;
    }
    if (ui.filter_scope === 'african_only' && !r.african_only) return false;
    if (ui.filter_scope === 'cross_species' && r.african_only) return false;
    return true;
  });
}

function _spRenderToolbar(rows, ui) {
  const cats = Array.from(new Set(rows.map(r => r.category)));
  const catOpts = ['<option value="all">all categories</option>']
    .concat(cats.map(c => '<option value="' + _esc(c) + '"' +
      (c === ui.filter_category ? ' selected' : '') + '>' + _esc(c) + '</option>')).join('');
  return '<div class="sp-toolbar">' +
    '<label class="sp-tb-label">View ' +
    '<select id="spViewMode">' +
      '<option value="detailed"' + (!ui.compact ? ' selected' : '') + '>detailed (6 cols)</option>' +
      '<option value="compact"'  + (ui.compact  ? ' selected' : '') + '>compact (4 cols)</option>' +
    '</select></label>' +
    '<label class="sp-tb-label">Category ' +
    '<select id="spFilterCat">' + catOpts + '</select></label>' +
    '<label class="sp-tb-label">Scope ' +
    '<select id="spFilterScope">' +
      '<option value="all"' + (ui.filter_scope === 'all' ? ' selected' : '') + '>all</option>' +
      '<option value="african_only"' + (ui.filter_scope === 'african_only' ? ' selected' : '') + '>African-only</option>' +
      '<option value="cross_species"' + (ui.filter_scope === 'cross_species' ? ' selected' : '') + '>cross-species</option>' +
    '</select></label>' +
    '<label class="sp-tb-label sp-tb-checkbox">' +
      '<input type="checkbox" id="spFilterSig"' + (ui.filter_significant_only ? ' checked' : '') + '> significant only' +
    '</label>' +
    '<div class="sp-tb-spacer"></div>' +
    '<button class="sp-tb-btn" id="spLoadJsonBtn" title="Load a stats_profile.json or .tsv to overlay rows the atlas can\'t derive itself.">load JSON / TSV\u2026</button>' +
    '<input type="file" id="spLoadInput" style="display:none" accept=".json,.tsv,.txt">' +
    '<button class="sp-tb-btn" id="spExportCsvBtn" title="Download the current view as a CSV.">export CSV</button>' +
    '<button class="sp-tb-btn" id="spExportSvgBtn" title="Download the table as a single-file SVG suitable for manuscript figures.">export SVG</button>' +
    '<button class="sp-tb-btn" id="spResetBtn" title="Drop user-supplied overlay; re-derive everything from atlas state.">reset</button>' +
    '</div>';
}

function _spEffectPill(effect) {
  if (!effect) effect = 'not_tested';
  const color = SP_EFFECT_COLORS[effect] || 'var(--ink-dim)';
  const label = String(effect).replace(/_/g, ' ');
  return '<span class="sp-effect-pill" style="color:' + color +
    ';border-color:' + color + ';" data-effect="' + _esc(effect) + '">' +
    _esc(label) + '</span>';
}

function _spStateBadge(s) {
  if (!s) s = 'placeholder';
  const color = SP_EFFECT_COLORS[s] || 'var(--ink-dim)';
  const label = (s === 'derived' ? 'atlas' : s === 'loaded' ? 'loaded' : 'add data');
  return '<span class="sp-state-badge" style="color:' + color +
    ';border-color:' + color + ';" data-state="' + _esc(s) + '">' +
    _esc(label) + '</span>';
}

function _spFmtP(p) {
  if (p == null) return '\u2014';
  if (p < 0.001) return p.toExponential(2);
  return p.toFixed(3);
}

function _spRenderResultCell(result, isAfrican) {
  if (!result) result = _spPlaceholder('no data');
  const eff = _spEffectPill(result.effect);
  const badge = _spStateBadge(result.state);
  const lbl = result.label ? '<div class="sp-cell-label">' + _esc(result.label) + '</div>' : '';
  const stats = (result.p_value != null || result.fold_change != null)
    ? '<div class="sp-cell-stats">' +
      (result.fold_change != null ? 'fold ' + result.fold_change.toFixed(2) + ' \u00B7 ' : '') +
      (result.p_value != null ? 'p=' + _spFmtP(result.p_value) : '') +
      (result.n != null ? ' \u00B7 n=' + result.n : '') +
      '</div>'
    : '';
  return '<div class="sp-cell-effect-row">' + eff + ' ' + badge + '</div>' + lbl + stats;
}

function _spRenderTable(rows, ui, species) {
  const filtered = _spFilteredRows(rows, ui);
  if (filtered.length === 0) {
    return '<div class="sp-empty">No rows match the current filters.</div>';
  }
  const isCompact = ui.compact;
  const sp1 = species.species_1, sp2 = species.species_2;
  // v4 turn 111: italicize species names per scientific convention.
  // species_1 / species_2 are typically "C. gariepinus" / "C. macrocephalus"
  // (full binomial or genus-abbrev + specific epithet) — both should be in
  // italics. We wrap _esc() output in <i>..</i> tags. _esc() guarantees no
  // HTML injection so wrapping is safe.
  const sp1Italic = '<i>' + _esc(sp1) + '</i>';
  const sp2Italic = '<i>' + _esc(sp2) + '</i>';
  const head = isCompact
    ? '<thead><tr>' +
        '<th>Statistic</th>' +
        '<th>' + sp1Italic + ' result</th>' +
        '<th>' + sp2Italic + ' result</th>' +
        '<th>Biological interpretation</th>' +
      '</tr></thead>'
    : '<thead><tr>' +
        '<th>Category</th><th>Statistic / test</th><th>Null comparison</th>' +
        '<th>' + sp1Italic + ' result</th>' +
        '<th>' + sp2Italic + ' result</th>' +
        '<th>Biological conclusion</th>' +
      '</tr></thead>';
  const body = filtered.map(r => {
    const af = _spRenderResultCell(r.african_result, true);
    const bh = _spRenderResultCell(r.bighead_result, false);
    if (isCompact) {
      return '<tr>' +
        '<td><div class="sp-stat-name">' + _esc(r.statistic) + '</div>' +
        '<div class="sp-stat-test">' + _esc(r.test || '') + '</div></td>' +
        '<td>' + af + '</td>' +
        '<td>' + bh + '</td>' +
        '<td class="sp-interp">' + _esc(r.interpretation || '') + '</td>' +
        '</tr>';
    }
    return '<tr>' +
      '<td><div class="sp-cat">' + _esc(r.category) + '</div></td>' +
      '<td><div class="sp-stat-name">' + _esc(r.statistic) + '</div>' +
        '<div class="sp-stat-test">' + _esc(r.test || '') + '</div></td>' +
      '<td class="sp-null">' + _esc(r.null_comparison || '') + '</td>' +
      '<td>' + af + '</td>' +
      '<td>' + bh + '</td>' +
      '<td class="sp-interp">' + _esc(r.interpretation || '') + '</td>' +
      '</tr>';
  }).join('');
  return '<table class="sp-table' + (isCompact ? ' sp-compact' : '') + '">' +
    head + '<tbody>' + body + '</tbody></table>';
}

// Render the whole page17 body. Idempotent — call any time atlas state
// changes (cs_breakpoints loaded, candidates added, etc).
function _renderStatsProfilePage() {
  const slot = document.getElementById('spBody');
  if (!slot) return;
  const ui = _spEnsureState();
  if (!ui.rows) _spRefreshDerivedRows();
  const cards = _spRenderSummaryCards(ui.rows);
  const toolbar = _spRenderToolbar(ui.rows, ui);
  const table = _spRenderTable(ui.rows, ui, ui.species);
  const legend = '<div class="sp-legend">' +
    '<div class="sp-legend-title">Effect tags:</div> ' +
    Object.keys(SP_EFFECT_COLORS).filter(k => k !== 'derived' && k !== 'loaded' && k !== 'placeholder')
      .map(_spEffectPill).join(' ') + '</div>';
  const methods = '<div class="sp-methods">' +
    '<b>Methods note.</b> Statistics are compared against chromosome-, size-, ' +
    'or feature-matched genomic backgrounds where possible. Population-based ' +
    'tests are performed in the African catfish broodstock panel; bighead ' +
    'catfish values refer to orthologous genomic context unless population ' +
    'resequencing data are available. Atlas-derived rows recompute live ' +
    'from <code>state.crossSpecies</code> and <code>state.candidateList</code>; ' +
    'rows requiring annotation tracks (gene density, GO/KEGG, ROH, deleterious ' +
    'burden, FST) accept a stats_profile.json or .tsv overlay.' +
    '</div>';
  slot.innerHTML = cards + toolbar + table + legend + methods;
  _spWireToolbar(slot);
}

function _spWireToolbar(slot) {
  const ui = _spEnsureState();
  const $vm = slot.querySelector('#spViewMode');
  if ($vm) $vm.addEventListener('change', () => {
    ui.compact = ($vm.value === 'compact');
    _renderStatsProfilePage();
  });
  const $fc = slot.querySelector('#spFilterCat');
  if ($fc) $fc.addEventListener('change', () => {
    ui.filter_category = $fc.value;
    _renderStatsProfilePage();
  });
  const $fs = slot.querySelector('#spFilterScope');
  if ($fs) $fs.addEventListener('change', () => {
    ui.filter_scope = $fs.value;
    _renderStatsProfilePage();
  });
  const $sig = slot.querySelector('#spFilterSig');
  if ($sig) $sig.addEventListener('change', () => {
    ui.filter_significant_only = $sig.checked;
    _renderStatsProfilePage();
  });
  const $load = slot.querySelector('#spLoadJsonBtn');
  const $input = slot.querySelector('#spLoadInput');
  if ($load && $input) {
    $load.addEventListener('click', () => $input.click());
    $input.addEventListener('change', e => {
      const f = e.target.files && e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = ev => _spIngestText(ev.target.result, f.name);
      r.readAsText(f);
      $input.value = '';
    });
  }
  const $csv = slot.querySelector('#spExportCsvBtn');
  if ($csv) $csv.addEventListener('click', () => _spExportCsv(ui));
  const $svg = slot.querySelector('#spExportSvgBtn');
  if ($svg) $svg.addEventListener('click', () => _spExportSvg(slot));
  const $reset = slot.querySelector('#spResetBtn');
  if ($reset) $reset.addEventListener('click', () => {
    _spRefreshDerivedRows();
    _renderStatsProfilePage();
  });
}

function _spIngestText(text, filename) {
  const ui = _spEnsureState();
  let parsed = null;
  if (filename && /\.tsv$|\.txt$/i.test(filename)) {
    parsed = _spParseTsv(text);
  } else {
    try { parsed = JSON.parse(text); } catch (e) {
      console.warn('[stats_profile] JSON parse failed:', e);
      try { parsed = _spParseTsv(text); } catch (_) { parsed = null; }
    }
  }
  if (!parsed) { console.warn('[stats_profile] could not parse', filename); return; }
  if (!_spIsValidProfileJson(parsed) && !Array.isArray(parsed.summary_rows)) {
    console.warn('[stats_profile] file does not look like a stats_profile schema');
    return;
  }
  _spRefreshDerivedRows();
  const merged = _spMergeUserJson(ui.rows, parsed);
  ui.rows = merged.rows;
  ui.species = merged.species;
  if (parsed.metadata) ui.metadata = Object.assign({}, ui.metadata, parsed.metadata);
  _renderStatsProfilePage();
}

// CSV export — current filtered view, full 17 columns regardless of compact mode
function _spExportCsv(ui) {
  const rows = _spFilteredRows(ui.rows, ui);
  const sp1 = ui.species.species_1, sp2 = ui.species.species_2;
  const header = ['category','statistic','test','null_comparison',
    sp1 + '_label', sp1 + '_effect', sp1 + '_observed', sp1 + '_expected', sp1 + '_p',
    sp2 + '_label', sp2 + '_effect', sp2 + '_observed', sp2 + '_expected', sp2 + '_p',
    'interpretation'];
  const csvCell = v => {
    if (v == null) return '';
    const s = String(v);
    if (/[",\n]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
    return s;
  };
  const lines = [header.map(csvCell).join(',')];
  for (const r of rows) {
    const af = r.african_result || {}, bh = r.bighead_result || {};
    lines.push([r.category, r.statistic, r.test, r.null_comparison,
      af.label, af.effect, af.observed, af.expected, af.p_value,
      bh.label, bh.effect, bh.observed, bh.expected, bh.p_value,
      r.interpretation].map(csvCell).join(','));
  }
  const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'stats_profile.csv';
  document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 250);
}

// SVG export — wraps the current rendered HTML table inside a foreignObject
// for a single-file figure-quality export.
function _spExportSvg(slot) {
  const tbl = slot.querySelector('table.sp-table');
  if (!tbl) return;
  // Get dimensions from the rendered table
  const w = Math.max(900, tbl.scrollWidth || tbl.offsetWidth || 1200);
  const h = Math.max(300, tbl.scrollHeight || tbl.offsetHeight || 600);
  // Minimal inline CSS — enough to keep the export self-contained
  const inlineCss =
    'body{margin:0;font-family:ui-monospace,monospace;color:#1c1c1f;background:#fff;}' +
    'table{border-collapse:collapse;width:100%;font-size:11px;}' +
    'th,td{padding:7px 9px;border-bottom:1px solid #e0e2e6;vertical-align:top;text-align:left;}' +
    'th{font-weight:600;background:#f4f5f7;}' +
    '.sp-effect-pill,.sp-state-badge{display:inline-block;padding:1px 7px;border:1px solid;border-radius:9px;font-size:10px;margin-right:4px;}' +
    '.sp-cell-label{font-size:10.5px;color:#454952;margin-top:3px;line-height:1.4;}' +
    '.sp-cell-stats{font-size:10px;color:#6e7782;margin-top:2px;}' +
    '.sp-cat{font-size:10px;text-transform:uppercase;letter-spacing:0.05em;color:#6e7782;}' +
    '.sp-stat-name{font-weight:500;}' +
    '.sp-stat-test{font-size:10px;color:#6e7782;font-style:italic;margin-top:2px;}' +
    '.sp-null{font-size:10px;color:#6e7782;}' +
    '.sp-interp{font-size:10.5px;line-height:1.4;}';
  const xml = '<svg xmlns="http://www.w3.org/2000/svg" width="' + w + '" height="' + h + '" viewBox="0 0 ' + w + ' ' + h + '">' +
    '<foreignObject width="100%" height="100%">' +
    '<div xmlns="http://www.w3.org/1999/xhtml">' +
    '<style>' + inlineCss + '</style>' +
    tbl.outerHTML +
    '</div></foreignObject></svg>';
  const blob = new Blob([xml], { type: 'image/svg+xml;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'stats_profile.svg';
  document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 250);
}

// ---------------------------------------------------------------------------
// End verbatim extraction.
// ---------------------------------------------------------------------------

// Public entry — alias the underscore-prefixed legacy name to a clean export.
export function renderStatsProfilePage() {
  return _renderStatsProfilePage();
}

// Constants other modules may want to read.
export { SP_EFFECT_COLORS, SP_ALPHA, SP_DEFAULT_ROWS };

// Internal helpers re-exported for the merge chat / tests. Treat as
// implementation detail — public API is renderStatsProfilePage().
export {
  _spEnsureState,
  _spRefreshDerivedRows,
  _spDeriveAllRows,
  _spDeriveCsPermutation,
  _spDeriveRepeatFlankSpalax,
  _spDeriveFusionFission,
  _spDeriveMarkerability,
  _spIsValidProfileJson,
  _spMergeUserJson,
  _spParseTsv,
  _spExportCsv,
  _spExportSvg,
  _spIngestText,
};

export const __MODULE_ID__ = 'inversion_catalogue/page17';
