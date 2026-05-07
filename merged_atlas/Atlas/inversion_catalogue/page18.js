// Atlas/inversion_catalogue/page18.js
// =============================================================================
// page18 — Marker readiness panel (private-indel architecture)
// (`<div id="page18">` contains `<div id="mpBody">`)
//
// Source: legacy/Inversion_atlas.html lines 29307–30160 (JS body)
//                                    lines  8134–8157  (HTML shell)
// Public entry: renderMarkerPanelPage()  (legacy: _renderMarkerPanelPage line 30045)
//
// Tier hierarchy (manuscript spec):
//   Tier 1 — private indel/SNP tag with clean dosage (HIGHEST)
//             AF_STD <= 0.02 AND AF_HET in [0.25, 0.75] AND AF_INV >= 0.80
//   Tier 2 — multi-marker haplotype panel OR strong-tag with imperfect het
//   Tier 3 — breakpoint PCR candidate (DEMOTED)
//   Tier 4 — exploratory (LOWEST)
//
// External dependencies:
//   TODO_MISSING(_esc)  — HTML-escape helper (high frequency, candidate for
//                         shared/dom_utils.js promotion)
//   global `state`      — reads:
//                           state.candidateList
//                           state.crossSpecies
//                           state._markerPanel
//                           state.markerThresholds
// =============================================================================

const state = (typeof window !== 'undefined' && window.state) ? window.state : {};

// ---------------------------------------------------------------------------
// VERBATIM extraction from legacy lines 29307–30160.
// ---------------------------------------------------------------------------


// =============================================================================
// v4 turn 78 — Marker readiness panel (page18) — REBUILT for private-indel
//              architecture. Breakpoint PCR demoted to Tier 3 because
//              breakpoint precision is uncertain even on the cs_breakpoint
//              catalogue.
// -----------------------------------------------------------------------------
// Tier hierarchy (manuscript spec):
//
//   Tier 1 — private indel/SNP tag with clean dosage (HIGHEST confidence)
//     AF_STD <= 0.02 AND AF_HET in [0.25, 0.75] AND AF_INV >= 0.80
//     Plus: positive controls in INV/INV and HET groups, negative controls
//           in STD/STD group, and bighead specificity confirmed.
//     Indel size 20-300 bp preferred for cheap gel-visible PCR.
//
//   Tier 2 — multi-marker haplotype panel OR strong-tag with imperfect het
//     >=3 linked private markers across the inversion (left/middle/right) OR
//     AF_STD <= 0.05 AND AF_INV >= 0.70 with controls but het not perfectly
//     dosage-shaped.
//
//   Tier 3 — breakpoint PCR candidate (DEMOTED from prior architecture)
//     Breakpoint within +/-50 kb of a cs_breakpoint AND wet-lab validation
//     pending. Lower confidence because the cs_breakpoint catalogue itself
//     can mis-localise breakpoints by 10s of kb.
//
//   Tier 4 — exploratory (LOWEST)
//     Association exists but controls incomplete or weak. For research only.
//
// AF scoring runs live in the browser:
//   private_score = AF_INV - AF_STD                 (0 to 1; higher = better)
//   dosage_score  = 1 - |AF_HET - 0.5| * 2          (0 to 1; 1 = perfect 50%)
//   final_score   = private_score + dosage_score    (0 to 2; >=1.7 = clean)
//
// Inputs:
//   variant_afs.json — per-inversion-per-variant AF table (atlas computes
//                      tier + score live)
//   marker_panel.json — user-curated overlay with primers, controls,
//                       cross-species notes (overrides tier)
//
// Auto-suggested controls come from state.candidateList karyotype data
// (each candidate carries assignments[sample] = 0|1|2 from K-means PCA on
// page 1). When that data is present, the panel proposes positive/negative
// control samples per group; user JSON can override or add cross-species
// controls (bighead negative, F1 hybrid expected pattern, no-template).
// =============================================================================

// Tier definitions. Tier 1 is now the private-indel tag, not breakpoint PCR.
const MP_TIER_DEFS = {
  1: {
    label: 'Tier 1 — private indel/SNP tag',
    short: 'Tier 1',
    color: '#3cc08a',
    description: 'Clean dosage signature: AF_STD\u00A0\u22640.02, AF_HET\u00A0\u2208\u00A0[0.25,\u00A00.75], AF_INV\u00A0\u22650.80. Positive controls (INV/INV and HET) and negative controls (STD/STD) identified, plus bighead specificity confirmed.',
    expected_use: 'direct genotyping (carrier vs non-carrier)',
  },
  2: {
    label: 'Tier 2 — multi-marker panel / strong tag',
    short: 'Tier 2',
    color: '#7ad3db',
    description: 'Either \u22653 linked private markers across the inversion (left/mid/right), OR a single strong tag (AF_STD\u00A0\u22640.05, AF_INV\u00A0\u22650.70) with controls but het not perfectly dosage-shaped.',
    expected_use: 'population screening (validate per batch)',
  },
  3: {
    label: 'Tier 3 — breakpoint PCR candidate',
    short: 'Tier 3',
    color: '#f0a35e',
    description: 'Breakpoint within \u00B150 kb of a cs_breakpoint, controls identified, wet-lab validation pending. Lower confidence than private indels because breakpoint precision itself is uncertain.',
    expected_use: 'wet-lab pilot validation',
  },
  4: {
    label: 'Tier 4 — exploratory',
    short: 'Tier 4',
    color: '#9aa3ad',
    description: 'Association exists but controls incomplete or weak. For research-only follow-up; do not use for breeding decisions.',
    expected_use: 'research follow-up only',
  },
};

const MP_VALIDATION_DEFS = {
  'validated':         { label: 'validated',         color: '#3cc08a' },
  'pilot_tested':      { label: 'pilot tested',      color: '#7ad3db' },
  'designed_only':     { label: 'designed only',     color: '#f0a35e' },
  'auto_assigned':     { label: 'auto-assigned',     color: '#9aa3ad' },
  'untested':          { label: 'untested',          color: '#6e7782' },
};

// Private-allele AF thresholds. Constants so the user can later override
// via state.markerThresholds if needed. Two thresholds: relaxed (Tier 2
// strong tag) and strict (Tier 1).
const MP_AF_THRESHOLDS = {
  tier1: { af_std_max: 0.02, af_het_min: 0.25, af_het_max: 0.75, af_inv_min: 0.80 },
  tier2: { af_std_max: 0.05, af_inv_min: 0.70 },
};

// Indel size band that's gel-visible on cheap agarose PCR
const MP_INDEL_GEL_VISIBLE_BP = { min: 20, max: 300 };

// Breakpoint PCR proximity radius (Tier 3 candidate)
const MP_TIER3_BREAKPOINT_RADIUS_BP = 50000;

// =============================================================================
// AF scoring — the core private-allele logic
// -----------------------------------------------------------------------------

// Score a single variant given its allele frequencies in three karyotype
// groups (STD/STD, STD/INV, INV/INV). Returns the manuscript-spec scores
// plus the assigned tier.
//
// Input: { af_std, af_het, af_inv, indel_size_bp (optional, +/-) }
// Output: { private_score, dosage_score, final_score, tier_from_af, ...inputs }
function _mpScoreVariantAf(v) {
  const af_std = (typeof v.af_std === 'number') ? v.af_std : null;
  const af_het = (typeof v.af_het === 'number') ? v.af_het : null;
  const af_inv = (typeof v.af_inv === 'number') ? v.af_inv : null;
  if (af_std == null || af_inv == null) {
    return Object.assign({}, v, {
      private_score: null, dosage_score: null, final_score: null,
      tier_from_af: 4, tier_reason: 'AF table incomplete (need af_std + af_inv)',
    });
  }
  const private_score = af_inv - af_std;
  const dosage_score = (af_het == null) ? null : (1 - Math.abs(af_het - 0.5) * 2);
  const final_score = (dosage_score != null) ? (private_score + dosage_score) : private_score;

  // Tier 1 requires clean dosage AND strict AF separation
  const t1 = MP_AF_THRESHOLDS.tier1;
  const meetsT1 = (af_std <= t1.af_std_max)
               && (af_inv >= t1.af_inv_min)
               && (af_het != null)
               && (af_het >= t1.af_het_min)
               && (af_het <= t1.af_het_max);
  const t2 = MP_AF_THRESHOLDS.tier2;
  const meetsT2Strong = (af_std <= t2.af_std_max) && (af_inv >= t2.af_inv_min);

  let tier_from_af, tier_reason;
  if (meetsT1) {
    tier_from_af = 1;
    tier_reason = 'AF clean: STD ' + af_std.toFixed(2) +
                  ', HET ' + af_het.toFixed(2) +
                  ', INV ' + af_inv.toFixed(2);
  } else if (meetsT2Strong) {
    tier_from_af = 2;
    tier_reason = 'strong tag: STD ' + af_std.toFixed(2) +
                  ', INV ' + af_inv.toFixed(2) +
                  (af_het != null ? ', HET ' + af_het.toFixed(2) + ' (imperfect dosage)' :
                                    ' (HET not measured)');
  } else {
    tier_from_af = 4;
    tier_reason = 'AF separation insufficient: STD ' + af_std.toFixed(2) +
                  ', INV ' + af_inv.toFixed(2);
  }
  return Object.assign({}, v, {
    private_score, dosage_score, final_score, tier_from_af, tier_reason,
  });
}

// Apply gel-visibility annotation to an indel marker. Returns the variant
// with `gel_visible: bool` added. SNPs return null (not gel-visible by PCR).
function _mpAnnotateGelVisibility(v) {
  const t = (v.type || '').toLowerCase();
  if (t !== 'indel' && t !== 'insertion' && t !== 'deletion') {
    return Object.assign({}, v, { gel_visible: null });
  }
  const size = Math.abs(v.indel_size_bp || 0);
  const gel_visible = (size >= MP_INDEL_GEL_VISIBLE_BP.min) &&
                      (size <= MP_INDEL_GEL_VISIBLE_BP.max);
  return Object.assign({}, v, { gel_visible });
}

// =============================================================================
// Auto-suggested positive/negative controls from candidate karyotype state
// -----------------------------------------------------------------------------

// For a candidate in state.candidateList that has karyotype assignments
// (cand.assignments = { sample_id: 0|1|2 }), return up to N sample IDs from
// each karyotype class. Used to pre-populate the controls column on every
// marker for that inversion.
//
// Output: { positive_controls_INV: [...], heterozygote_controls: [...],
//           negative_controls_STD: [...], n_per_class: {...} }
function _mpSuggestControlsFromKaryotype(cand, nPerClass) {
  if (!cand || typeof cand !== 'object') return null;
  nPerClass = nPerClass || 3;
  // Karyotype assignments may live under several names depending on which
  // page promoted the candidate. Check the common ones.
  const assn = cand.assignments || cand.karyotype || cand.assignment || null;
  if (!assn || typeof assn !== 'object') return null;
  const groups = { 0: [], 1: [], 2: [] };
  for (const sample in assn) {
    const k = assn[sample];
    if (k === 0 || k === 1 || k === 2) groups[k].push(sample);
  }
  // Sample N per group; if user has confidence scores attached, sort by them
  const sample = arr => arr.slice(0, nPerClass);
  return {
    negative_controls_STD: sample(groups[0]),
    heterozygote_controls: sample(groups[1]),
    positive_controls_INV: sample(groups[2]),
    n_per_class: { 0: groups[0].length, 1: groups[1].length, 2: groups[2].length },
  };
}

// Default cross-species control template — added to every marker's
// `controls` block so the manuscript-required cross-species checks are
// always visible. Users override via marker_panel.json.
function _mpDefaultCrossSpeciesControls() {
  return {
    bighead_negative_required: true,
    bighead_negative_status: 'not_yet_tested',
    bighead_orthologous_sequence_status: 'unknown',
    bighead_expected_pattern: 'no amplification, OR clearly different size band, OR non-INV allele',
    f1_hybrid_expected_pattern: 'African paternal INV allele + bighead maternal allele (codominant)',
    no_template_control_required: true,
  };
}

// =============================================================================
// Auto-tier from atlas state
// -----------------------------------------------------------------------------

// Find the closest cs_breakpoint to a candidate's endpoints. Returns
// Infinity if no cs_breakpoints are loaded.
function _mpMinDistanceToCsBreakpoint(chrom, startBp, endBp) {
  if (!state.crossSpecies || !Array.isArray(state.crossSpecies.breakpoints)) return Infinity;
  let best = Infinity;
  for (const bp of state.crossSpecies.breakpoints) {
    if (bp.gar_chr !== chrom) continue;
    const bpPos = (typeof bp.gar_pos_start === 'number')
      ? bp.gar_pos_start
      : (typeof bp.gar_pos_mb === 'number' ? bp.gar_pos_mb * 1e6 : null);
    if (bpPos == null) continue;
    const d = Math.min(Math.abs(bpPos - startBp), Math.abs(bpPos - endBp));
    if (d < best) best = d;
  }
  return best;
}

// Auto-tier a candidate based on what the atlas knows about it. This is the
// "starting tier" before AF data or user uploads override it.
//
// Tier 4 is the default: an inversion call with no quantitative evidence
// for marker design. Promote to Tier 3 only when a cs_breakpoint is nearby
// (the candidate could become a breakpoint PCR target with wet-lab work).
// Tier 1 and Tier 2 require AF data — they're never assigned by atlas
// state alone.
function _mpAutoTierFromAtlas(cand) {
  const start = cand.start_bp != null ? cand.start_bp : cand.startBp;
  const end   = cand.end_bp   != null ? cand.end_bp   : cand.endBp;
  const distBp = _mpMinDistanceToCsBreakpoint(cand.chrom, start, end);
  const evidence = [];
  if (Number.isFinite(distBp) && distBp <= MP_TIER3_BREAKPOINT_RADIUS_BP) {
    evidence.push('cs_breakpoint within ' + (distBp < 1000 ? distBp + ' bp' :
                                             (distBp / 1000).toFixed(1) + ' kb'));
    return {
      tier: 3, marker_type: 'breakpoint PCR (candidate)',
      evidence, dist_to_bp_bp: distBp,
      tier_reason: 'cs_breakpoint nearby; wet-lab validation pending',
    };
  }
  evidence.push('no cs_breakpoint within \u00B1' +
                (MP_TIER3_BREAKPOINT_RADIUS_BP / 1000) + ' kb; no AF data loaded');
  return {
    tier: 4, marker_type: 'exploratory',
    evidence, dist_to_bp_bp: Number.isFinite(distBp) ? distBp : null,
    tier_reason: 'no quantitative evidence; load variant_afs.json to score',
  };
}

// Build the full marker panel from state.candidateList. Each candidate
// becomes one row; AF-derived sub-rows are added when variant_afs.json
// has been loaded for that candidate.
function _mpDeriveAutoPanel() {
  const cands = (state.candidateList || []).filter(c =>
    c && c.chrom && (Number.isFinite(c.start_bp) || Number.isFinite(c.startBp)));
  return cands.map(c => {
    const start = c.start_bp != null ? c.start_bp : c.startBp;
    const end   = c.end_bp   != null ? c.end_bp   : c.endBp;
    const t = _mpAutoTierFromAtlas(c);
    const id = c.id || c.candidate_id || ('INV_' + Math.random().toString(36).slice(2, 7));
    const controls = _mpSuggestControlsFromKaryotype(c) || {
      positive_controls_INV: [], heterozygote_controls: [], negative_controls_STD: [],
      n_per_class: { 0: 0, 1: 0, 2: 0 },
    };
    return {
      inversion_id: id,
      chrom: c.chrom,
      start_bp: start,
      end_bp: end,
      len_bp: (Number.isFinite(start) && Number.isFinite(end)) ? (end - start) : null,
      confirmed: c.confirmed === true,
      tier: t.tier,
      marker_type: t.marker_type,
      expected_use: MP_TIER_DEFS[t.tier].expected_use,
      tier_reason: t.tier_reason,
      evidence: t.evidence,
      dist_to_bp_bp: t.dist_to_bp_bp,
      validation_status: 'auto_assigned',
      // AF-scored variants (filled by _mpAnnotateAfVariants when
      // variant_afs.json is loaded for this inversion):
      af_variants: [],
      best_af_score: null,
      // Controls — auto-suggested from karyotype, user-overridable
      controls: Object.assign(controls, _mpDefaultCrossSpeciesControls()),
      // User-supplied fields (populated by marker_panel.json):
      primer_F: null, primer_R: null,
      amplicon_bp_state_A: null, amplicon_bp_state_B: null,
      n_carriers_tested: null,
      notes: null,
      source: 'derived',
    };
  });
}

// Annotate the panel with AF-scored variants from a variant_afs.json upload.
// The JSON has shape:
//   { metadata: {...}, variants_by_inversion: { <inv_id>: [<variant>, ...] } }
// where each variant has chr, pos, ref, alt, type, indel_size_bp, af_std,
// af_het, af_inv. We score each, attach to the matching panel row, and
// re-evaluate the row's overall tier (the AF tier overrides the atlas tier
// when AF tier is better — i.e. lower number).
function _mpAnnotateAfVariants(panel, afJson) {
  if (!afJson || !afJson.variants_by_inversion) return panel;
  const byInv = afJson.variants_by_inversion;
  for (const row of panel) {
    const variants = byInv[row.inversion_id];
    if (!Array.isArray(variants)) continue;
    const scored = variants.map(v => _mpAnnotateGelVisibility(_mpScoreVariantAf(v)));
    // Sort by final_score desc (best markers first)
    scored.sort((a, b) => (b.final_score || -Infinity) - (a.final_score || -Infinity));
    row.af_variants = scored;
    row.best_af_score = scored.length > 0 ? scored[0].final_score : null;
    // Re-tier based on best AF score:
    const best = scored[0];
    if (best && best.tier_from_af < row.tier) {
      // Promote: AF evidence is stronger than the atlas-derived tier
      row.evidence = (row.evidence || []).concat([
        'AF-derived tier ' + best.tier_from_af + ' (' + best.tier_reason + ')',
      ]);
      row.tier = best.tier_from_af;
      row.marker_type = (best.tier_from_af === 1) ? 'private indel/SNP tag'
                      : (best.tier_from_af === 2) ? 'multi-marker / strong tag'
                      : row.marker_type;
      row.expected_use = MP_TIER_DEFS[row.tier].expected_use;
      row.tier_reason = best.tier_reason;
    }
    // If multiple Tier-1 markers exist across the inversion (left/mid/right),
    // explicitly note that this inversion has a multi-marker panel option.
    const t1Count = scored.filter(s => s.tier_from_af === 1).length;
    if (t1Count >= 3) {
      row.evidence = (row.evidence || []).concat([
        t1Count + ' Tier-1 private markers \u2014 multi-marker panel feasible',
      ]);
    }
  }
  return panel;
}

// =============================================================================
// User-supplied marker_panel.json overlay — same shape as before but with
// expanded controls block and tier 4
// -----------------------------------------------------------------------------

function _mpIsValidPanelJson(parsed) {
  return parsed && typeof parsed === 'object' && Array.isArray(parsed.markers);
}

function _mpIsValidVariantAfsJson(parsed) {
  return parsed && typeof parsed === 'object' &&
         parsed.variants_by_inversion &&
         typeof parsed.variants_by_inversion === 'object';
}

function _mpMergeUserJson(autoPanel, parsed) {
  const byId = new Map(autoPanel.map(r => [r.inversion_id, r]));
  for (const um of parsed.markers) {
    const id = um.inversion_id || um.id;
    if (!id) continue;
    let target = byId.get(id);
    if (!target) {
      target = {
        inversion_id: id,
        chrom: um.chrom || '?',
        start_bp: um.start_bp != null ? um.start_bp : null,
        end_bp:   um.end_bp   != null ? um.end_bp   : null,
        len_bp:   null,
        confirmed: false,
        tier: 4, marker_type: 'exploratory',
        expected_use: MP_TIER_DEFS[4].expected_use,
        tier_reason: 'user-supplied marker (no atlas candidate match)',
        evidence: ['user-supplied (not in current atlas candidate list)'],
        dist_to_bp_bp: null,
        validation_status: 'untested',
        af_variants: [],
        best_af_score: null,
        controls: Object.assign({
          positive_controls_INV: [], heterozygote_controls: [],
          negative_controls_STD: [],
          n_per_class: { 0: 0, 1: 0, 2: 0 },
        }, _mpDefaultCrossSpeciesControls()),
        primer_F: null, primer_R: null,
        amplicon_bp_state_A: null, amplicon_bp_state_B: null,
        n_carriers_tested: null,
        notes: null,
        source: 'user',
      };
      autoPanel.push(target); byId.set(id, target);
    }
    if (typeof um.tier === 'number' && um.tier >= 1 && um.tier <= 4) {
      if (target.tier !== um.tier) {
        target.evidence = (target.evidence || []).concat(
          ['auto-tier was ' + target.tier + ', overridden to ' + um.tier + ' by user panel']);
      }
      target.tier = um.tier;
      target.marker_type = um.marker_type || MP_TIER_DEFS[um.tier].label.split('—')[1].trim();
      target.expected_use = um.expected_use || MP_TIER_DEFS[um.tier].expected_use;
    }
    if (um.marker_type)  target.marker_type  = um.marker_type;
    if (um.expected_use) target.expected_use = um.expected_use;
    if (um.primer_F)     target.primer_F     = um.primer_F;
    if (um.primer_R)     target.primer_R     = um.primer_R;
    if (um.amplicon_bp_state_A != null) target.amplicon_bp_state_A = um.amplicon_bp_state_A;
    if (um.amplicon_bp_state_B != null) target.amplicon_bp_state_B = um.amplicon_bp_state_B;
    if (um.n_carriers_tested  != null) target.n_carriers_tested  = um.n_carriers_tested;
    if (um.validation_status) target.validation_status = um.validation_status;
    if (um.notes) target.notes = um.notes;
    // Controls: user overrides per-field, others inherit auto-suggested
    if (um.controls && typeof um.controls === 'object') {
      target.controls = Object.assign({}, target.controls, um.controls);
    }
    if (target.source !== 'user') target.source = 'overlaid';
  }
  return autoPanel;
}

function _mpParseTsv(text) {
  const lines = text.split(/\r?\n/).filter(l => l.trim().length > 0);
  if (lines.length < 2) return null;
  const header = lines[0].split('\t').map(s => s.trim());
  const idx = name => header.indexOf(name);
  const markers = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split('\t');
    const get = name => idx(name) >= 0 ? (parts[idx(name)] || '').trim() : '';
    const num = name => { const v = get(name); return (v === '' || v === 'NA') ? null : Number(v); };
    markers.push({
      inversion_id:        get('inversion_id'),
      chrom:               get('chrom'),
      start_bp:            num('start_bp'),
      end_bp:              num('end_bp'),
      tier:                num('tier'),
      marker_type:         get('marker_type'),
      expected_use:        get('expected_use'),
      validation_status:   get('validation_status'),
      primer_F:            get('primer_F'),
      primer_R:            get('primer_R'),
      amplicon_bp_state_A: num('amplicon_bp_state_A'),
      amplicon_bp_state_B: num('amplicon_bp_state_B'),
      n_carriers_tested:   num('n_carriers_tested'),
      notes:               get('notes'),
    });
  }
  return { metadata: {}, markers };
}

// =============================================================================
// Page state
// -----------------------------------------------------------------------------
function _mpEnsureState() {
  if (!state._markerPanel) {
    state._markerPanel = {
      filter_tier: 'all',
      filter_validation: 'all',
      filter_chrom: 'all',
      panel: null,
      af_data: null,         // last-loaded variant_afs.json
      metadata: { panel_name: 'Inversion marker readiness panel' },
    };
  }
  return state._markerPanel;
}

function _mpRefreshPanel() {
  const ui = _mpEnsureState();
  ui.panel = _mpDeriveAutoPanel();
  // Re-apply AF annotation if we previously loaded a variant_afs.json
  if (ui.af_data) ui.panel = _mpAnnotateAfVariants(ui.panel, ui.af_data);
}

// =============================================================================
// Rendering
// -----------------------------------------------------------------------------
function _mpFmtBp(bp) {
  if (bp == null) return '\u2014';
  if (bp >= 1e6)  return (bp / 1e6).toFixed(2) + ' Mb';
  if (bp >= 1000) return (bp / 1000).toFixed(1) + ' kb';
  return bp + ' bp';
}

function _mpRenderTierDefsCard() {
  const items = [1, 2, 3, 4].map(t => {
    const def = MP_TIER_DEFS[t];
    return '<div class="mp-tier-def" data-tier="' + t + '">' +
      '<div class="mp-tier-pill" style="color:' + def.color + ';border-color:' + def.color + ';">' +
      _esc(def.label) + '</div>' +
      '<div class="mp-tier-desc">' + _esc(def.description) + '</div>' +
      '<div class="mp-tier-use">expected use: <b>' + _esc(def.expected_use) + '</b></div>' +
      '</div>';
  }).join('');
  return '<div class="mp-tier-defs">' + items + '</div>';
}

function _mpRenderSummaryCards(panel) {
  const counts = { 1: 0, 2: 0, 3: 0, 4: 0 };
  let nValidated = 0, nUntested = 0, nWithControls = 0;
  for (const m of panel) {
    counts[m.tier] = (counts[m.tier] || 0) + 1;
    if (m.validation_status === 'validated' || m.validation_status === 'pilot_tested') nValidated += 1;
    if (m.validation_status === 'untested' || m.validation_status === 'auto_assigned') nUntested += 1;
    if (m.controls &&
        m.controls.positive_controls_INV && m.controls.positive_controls_INV.length > 0 &&
        m.controls.negative_controls_STD && m.controls.negative_controls_STD.length > 0) {
      nWithControls += 1;
    }
  }
  const card = (label, value, hint, color) =>
    '<div class="mp-card"' + (color ? ' style="border-left:3px solid ' + color + '"' : '') + '>' +
    '<div class="mp-card-value">' + value + '</div>' +
    '<div class="mp-card-label">' + label + '</div>' +
    (hint ? '<div class="mp-card-hint">' + hint + '</div>' : '') + '</div>';
  return '<div class="mp-cards">' +
    card('Tier 1 \u2014 private tag', counts[1], 'clean dosage \u00b7 controls confirmed', MP_TIER_DEFS[1].color) +
    card('Tier 2 \u2014 strong/multi',  counts[2], '\u22653 linked OR strong tag', MP_TIER_DEFS[2].color) +
    card('Tier 3 \u2014 breakpoint PCR', counts[3], 'wet-lab pending',                MP_TIER_DEFS[3].color) +
    card('Tier 4 \u2014 exploratory',   counts[4], 'controls weak \u00b7 research only', MP_TIER_DEFS[4].color) +
    card('Controls identified', nWithControls + ' / ' + panel.length,
         'pos + neg samples present',                                         '#7ad3db') +
    card('Validation', nValidated + ' / ' + panel.length,
         'validated or pilot-tested (' + nUntested + ' awaiting)',            '#7ad3db') +
    '</div>';
}

function _mpRenderToolbar(panel, ui) {
  const chroms = Array.from(new Set(panel.map(m => m.chrom).filter(c => c)));
  chroms.sort((a, b) => {
    const na = String(a).match(/(\d+)/), nb = String(b).match(/(\d+)/);
    if (na && nb) return parseInt(na[1], 10) - parseInt(nb[1], 10);
    return String(a).localeCompare(String(b));
  });
  const chromOpts = ['<option value="all">all chroms</option>']
    .concat(chroms.map(c => '<option value="' + _esc(c) + '"' +
      (c === ui.filter_chrom ? ' selected' : '') + '>' + _esc(c) + '</option>')).join('');
  return '<div class="mp-toolbar">' +
    '<label class="mp-tb-label">Tier ' +
    '<select id="mpFilterTier">' +
      '<option value="all"' + (ui.filter_tier === 'all' ? ' selected' : '') + '>all tiers</option>' +
      '<option value="1"' + (ui.filter_tier === '1' ? ' selected' : '') + '>Tier 1 only</option>' +
      '<option value="2"' + (ui.filter_tier === '2' ? ' selected' : '') + '>Tier 2 only</option>' +
      '<option value="3"' + (ui.filter_tier === '3' ? ' selected' : '') + '>Tier 3 only</option>' +
      '<option value="4"' + (ui.filter_tier === '4' ? ' selected' : '') + '>Tier 4 only</option>' +
    '</select></label>' +
    '<label class="mp-tb-label">Validation ' +
    '<select id="mpFilterValidation">' +
      '<option value="all"' + (ui.filter_validation === 'all' ? ' selected' : '') + '>any status</option>' +
      Object.keys(MP_VALIDATION_DEFS).map(k =>
        '<option value="' + k + '"' + (ui.filter_validation === k ? ' selected' : '') + '>' +
        MP_VALIDATION_DEFS[k].label + '</option>').join('') +
    '</select></label>' +
    '<label class="mp-tb-label">Chrom ' +
    '<select id="mpFilterChrom">' + chromOpts + '</select></label>' +
    '<div class="mp-tb-spacer"></div>' +
    '<button class="mp-tb-btn" id="mpLoadAfBtn" title="Upload variant_afs.json with per-variant AF in STD/HET/INV groups. Atlas computes private_score, dosage_score, gel-visibility, and assigns tiers.">load variant AFs\u2026</button>' +
    '<button class="mp-tb-btn" id="mpLoadJsonBtn" title="Upload marker_panel.json with primers, controls, validation status. Overrides auto-tiering.">load marker panel\u2026</button>' +
    '<input type="file" id="mpLoadInput" style="display:none" accept=".json,.tsv,.txt">' +
    '<button class="mp-tb-btn" id="mpExportCsvBtn" title="Download the marker panel as CSV.">export CSV</button>' +
    '<button class="mp-tb-btn" id="mpResetBtn" title="Drop all overlays; re-derive from atlas state only.">reset</button>' +
    '</div>';
}

function _mpFilteredPanel(panel, ui) {
  return panel.filter(m => {
    if (ui.filter_tier !== 'all' && String(m.tier) !== ui.filter_tier) return false;
    if (ui.filter_validation !== 'all' && m.validation_status !== ui.filter_validation) return false;
    if (ui.filter_chrom !== 'all' && m.chrom !== ui.filter_chrom) return false;
    return true;
  });
}

function _mpTierPill(tier) {
  const def = MP_TIER_DEFS[tier];
  if (!def) return '\u2014';
  return '<span class="mp-tier-badge" style="color:' + def.color + ';border-color:' + def.color + ';">' +
    _esc(def.short) + '</span>';
}

function _mpValidationPill(status) {
  const def = MP_VALIDATION_DEFS[status] || MP_VALIDATION_DEFS['untested'];
  return '<span class="mp-val-pill" style="color:' + def.color + ';border-color:' + def.color + ';">' +
    _esc(def.label) + '</span>';
}

function _mpRenderControlsCell(controls) {
  if (!controls) return '<span class="mp-no-ctrls">\u2014</span>';
  const fmt = arr => Array.isArray(arr) && arr.length
    ? arr.slice(0, 3).map(s => '<code>' + _esc(s) + '</code>').join(', ') +
      (arr.length > 3 ? ' \u00b7 +' + (arr.length - 3) + ' more' : '')
    : '<span class="mp-no-ctrls">none</span>';
  const bgStatus = controls.bighead_negative_status || 'unknown';
  return '<div class="mp-ctrl-block">' +
    '<div><span class="mp-ctrl-lbl">+ INV/INV:</span> ' + fmt(controls.positive_controls_INV) + '</div>' +
    '<div><span class="mp-ctrl-lbl">+ HET:</span> ' + fmt(controls.heterozygote_controls) + '</div>' +
    '<div><span class="mp-ctrl-lbl">- STD/STD:</span> ' + fmt(controls.negative_controls_STD) + '</div>' +
    '<div class="mp-ctrl-xs">cross-species: bighead ' + _esc(bgStatus) + ' \u00b7 NTC required</div>' +
    '</div>';
}

function _mpRenderAfBlock(m) {
  if (!Array.isArray(m.af_variants) || m.af_variants.length === 0) {
    return '<span class="mp-no-af">load variant_afs.json</span>';
  }
  // Show top 3 by final_score
  const top = m.af_variants.slice(0, 3);
  const rows = top.map(v => {
    const fs = (v.final_score != null) ? v.final_score.toFixed(2) : '\u2014';
    const ps = (v.private_score != null) ? v.private_score.toFixed(2) : '\u2014';
    const ds = (v.dosage_score != null) ? v.dosage_score.toFixed(2) : '\u2014';
    const t  = v.tier_from_af;
    const tcolor = MP_TIER_DEFS[t] ? MP_TIER_DEFS[t].color : '#9aa3ad';
    const gel = (v.gel_visible === true) ? ' \u00b7 <span class="mp-af-gel">gel\u00d7</span>' :
                (v.gel_visible === false) ? ' \u00b7 <span class="mp-af-nogel">no-gel</span>' : '';
    const indelInfo = (v.indel_size_bp != null) ? ' (' + (v.indel_size_bp > 0 ? '+' : '') + v.indel_size_bp + ' bp)' : '';
    return '<div class="mp-af-row">' +
      '<code>' + _esc(v.chr || m.chrom) + ':' + (v.pos != null ? v.pos.toLocaleString() : '?') + '</code> ' +
      _esc(v.type || 'snp') + indelInfo +
      ' \u00b7 STD ' + (v.af_std != null ? v.af_std.toFixed(2) : '\u2014') +
      ' / HET '       + (v.af_het != null ? v.af_het.toFixed(2) : '\u2014') +
      ' / INV '       + (v.af_inv != null ? v.af_inv.toFixed(2) : '\u2014') +
      ' \u00b7 score <b style="color:' + tcolor + ';">' + fs + '</b>' +
      ' (priv ' + ps + ', dosage ' + ds + ')' + gel +
      '</div>';
  }).join('');
  const more = m.af_variants.length > 3
    ? '<div class="mp-af-more">+' + (m.af_variants.length - 3) + ' more variants scored</div>'
    : '';
  return '<div class="mp-af-block">' + rows + more + '</div>';
}

function _mpRenderTable(panel, ui) {
  const filtered = _mpFilteredPanel(panel, ui);
  if (filtered.length === 0) {
    return '<div class="mp-empty">No markers match the current filters. ' +
           'Promote candidates on page 2 (focus) or page 3 (catalogue) to populate the panel, ' +
           'or upload a <code>marker_panel.json</code> / <code>variant_afs.json</code>.</div>';
  }
  const sorted = filtered.slice().sort((a, b) => {
    if (a.tier !== b.tier) return a.tier - b.tier;
    const na = String(a.chrom).match(/(\d+)/), nb = String(b.chrom).match(/(\d+)/);
    if (na && nb) {
      const d = parseInt(na[1], 10) - parseInt(nb[1], 10);
      if (d !== 0) return d;
    }
    return (a.start_bp || 0) - (b.start_bp || 0);
  });
  const head = '<thead><tr>' +
    '<th>inversion_id</th>' +
    '<th>chrom \u00b7 position</th>' +
    '<th>tier</th>' +
    '<th>marker type</th>' +
    '<th>top scored variants (AF)</th>' +
    '<th>controls</th>' +
    '<th>validation</th>' +
    '<th>tier reason</th>' +
    '</tr></thead>';
  const rows = sorted.map(m => {
    const pos = (m.start_bp != null && m.end_bp != null)
      ? (_mpFmtBp(m.start_bp) + '\u2013' + _mpFmtBp(m.end_bp) +
         '<br><span class="mp-len">len ' + _mpFmtBp(m.len_bp) + '</span>')
      : '\u2014';
    return '<tr data-tier="' + m.tier + '">' +
      '<td><code>' + _esc(m.inversion_id) + '</code></td>' +
      '<td><code>' + _esc(m.chrom) + '</code><br>' + pos + '</td>' +
      '<td>' + _mpTierPill(m.tier) + '</td>' +
      '<td>' + _esc(m.marker_type) + '<br><span class="mp-use">' + _esc(m.expected_use) + '</span></td>' +
      '<td>' + _mpRenderAfBlock(m) + '</td>' +
      '<td>' + _mpRenderControlsCell(m.controls) + '</td>' +
      '<td>' + _mpValidationPill(m.validation_status) + '</td>' +
      '<td class="mp-evidence">' + _esc(m.tier_reason || '') + '</td>' +
      '</tr>';
  }).join('');
  return '<table class="mp-table">' + head + '<tbody>' + rows + '</tbody></table>';
}

function _mpRenderPilotPlan() {
  return '<div class="mp-pilot">' +
    '<div class="mp-pilot-title">Pilot validation plan (with controls)</div>' +
    '<div class="mp-pilot-body">' +
      'Before scaling up to a full breeding panel, validate a small pilot subset of 10\u201320 ' +
      'low-cost markers across the highest-confidence inversions, with positive AND negative ' +
      'controls in every PCR run. The 5-step sequence:' +
    '</div>' +
    '<ol class="mp-pilot-steps">' +
      '<li>Run each marker on positive controls (INV/INV and HET samples from the resequencing dataset).</li>' +
      '<li>Run the same marker on negative controls (STD/STD samples) and a no-template control (water).</li>' +
      '<li>Run on bighead catfish DNA to confirm cross-species specificity \u2014 the marker should NOT amplify the bighead allele in the same way.</li>' +
      '<li>Test on 20\u201330 unrelated broodstock from a different batch; keep only markers with clean genotype separation across all controls.</li>' +
      '<li>Then expand to larger breeding populations.</li>' +
    '</ol>' +
    '<div class="mp-pilot-footer">' +
      'Tier 1 markers are the recommended starting set. A marker without positive AND negative controls ' +
      'cannot be Tier 1 \u2014 the controls column on every row reflects this requirement.' +
    '</div>' +
    '</div>';
}

function _mpRenderMethods() {
  return '<div class="mp-methods">' +
    '<b>Methods note.</b> Markers are scored from per-variant allele frequencies in three ' +
    'karyotype groups (STD/STD, STD/INV, INV/INV) inferred from K-means PCA dosage. ' +
    'A variant is assigned <b>Tier 1</b> when it shows clean dosage ' +
    '(AF_STD\u00A0\u22640.02, AF_HET\u00A0\u2208\u00A0[0.25,\u00A00.75], AF_INV\u00A0\u22650.80) ' +
    'AND positive/negative controls are identified in all three groups AND bighead specificity is ' +
    'confirmed; <b>Tier 2</b> requires either \u22653 linked private markers across the inversion or ' +
    'a strong tag (AF_STD\u00A0\u22640.05, AF_INV\u00A0\u22650.70) with controls but imperfect het; ' +
    '<b>Tier 3</b> covers breakpoint PCR candidates within \u00B150 kb of a cs_breakpoint, controls ' +
    'identified, wet-lab validation pending; <b>Tier 4</b> covers exploratory markers with ' +
    'incomplete or weak controls. ' +
    '<b>Because the target breeding system involves interspecific F\u2081 hybrids, candidate markers ' +
    'were additionally screened against the bighead catfish orthologous genome to distinguish African ' +
    'inversion tags from maternal-species sequence differences.</b> ' +
    'Candidate markers are provided as a starting resource and should be locally validated before ' +
    'routine breeding deployment. No marker on this panel constitutes a guarantee of genotype-' +
    'phenotype association without independent validation in the target population. ' +
    '<br><br>' +
    'See also: <a href="Population_atlas.html#page8" style="color:#7ad3db; text-decoration:none;">' +
    'Population Atlas \u2192 Breeding</a> for sample-level highlights of marker-validation control ' +
    'samples and other broodstock flagged for committee review.' +
    '</div>';
}

function _renderMarkerPanelPage() {
  const slot = document.getElementById('mpBody');
  if (!slot) return;
  const ui = _mpEnsureState();
  if (!ui.panel) _mpRefreshPanel();
  const tierDefs = _mpRenderTierDefsCard();
  const cards = _mpRenderSummaryCards(ui.panel);
  const toolbar = _mpRenderToolbar(ui.panel, ui);
  const table = _mpRenderTable(ui.panel, ui);
  const pilot = _mpRenderPilotPlan();
  const methods = _mpRenderMethods();
  slot.innerHTML = tierDefs + cards + toolbar + table + pilot + methods;
  _mpWireToolbar(slot);
}

function _mpWireToolbar(slot) {
  const ui = _mpEnsureState();
  const $ft = slot.querySelector('#mpFilterTier');
  if ($ft) $ft.addEventListener('change', () => { ui.filter_tier = $ft.value; _renderMarkerPanelPage(); });
  const $fv = slot.querySelector('#mpFilterValidation');
  if ($fv) $fv.addEventListener('change', () => { ui.filter_validation = $fv.value; _renderMarkerPanelPage(); });
  const $fc = slot.querySelector('#mpFilterChrom');
  if ($fc) $fc.addEventListener('change', () => { ui.filter_chrom = $fc.value; _renderMarkerPanelPage(); });
  const $loadAf = slot.querySelector('#mpLoadAfBtn');
  const $load   = slot.querySelector('#mpLoadJsonBtn');
  const $input  = slot.querySelector('#mpLoadInput');
  if ($loadAf && $input) {
    $loadAf.addEventListener('click', () => { ui._loadKind = 'af'; $input.click(); });
  }
  if ($load && $input) {
    $load.addEventListener('click', () => { ui._loadKind = 'panel'; $input.click(); });
  }
  if ($input) {
    $input.addEventListener('change', e => {
      const f = e.target.files && e.target.files[0];
      if (!f) return;
      const r = new FileReader();
      r.onload = ev => _mpIngestText(ev.target.result, f.name, ui._loadKind || 'auto');
      r.readAsText(f);
      $input.value = '';
    });
  }
  const $csv = slot.querySelector('#mpExportCsvBtn');
  if ($csv) $csv.addEventListener('click', () => _mpExportCsv(ui));
  const $reset = slot.querySelector('#mpResetBtn');
  if ($reset) $reset.addEventListener('click', () => {
    ui.af_data = null;
    _mpRefreshPanel();
    _renderMarkerPanelPage();
  });
}

function _mpIngestText(text, filename, kind) {
  const ui = _mpEnsureState();
  let parsed = null;
  if (filename && /\.tsv$|\.txt$/i.test(filename)) {
    parsed = _mpParseTsv(text);
  } else {
    try { parsed = JSON.parse(text); } catch (_) {
      try { parsed = _mpParseTsv(text); } catch (__) { parsed = null; }
    }
  }
  if (!parsed) { console.warn('[markerPanel] could not parse', filename); return; }
  // Auto-detect: variant_afs or marker_panel
  const isAfs = _mpIsValidVariantAfsJson(parsed);
  const isPanel = _mpIsValidPanelJson(parsed);
  if (kind === 'af' || (kind === 'auto' && isAfs)) {
    if (!isAfs) { console.warn('[markerPanel] expected variant_afs schema'); return; }
    ui.af_data = parsed;
    _mpRefreshPanel();
    if (parsed.metadata) ui.metadata = Object.assign({}, ui.metadata, parsed.metadata);
  } else {
    if (!isPanel) { console.warn('[markerPanel] expected marker_panel schema'); return; }
    _mpRefreshPanel();
    ui.panel = _mpMergeUserJson(ui.panel, parsed);
    if (parsed.metadata) ui.metadata = Object.assign({}, ui.metadata, parsed.metadata);
  }
  _renderMarkerPanelPage();
}

function _mpExportCsv(ui) {
  const filtered = _mpFilteredPanel(ui.panel, ui);
  const header = ['inversion_id','chrom','start_bp','end_bp','tier','marker_type',
                  'expected_use','validation_status','tier_reason',
                  'best_af_score','n_af_variants','n_tier1_variants',
                  'positive_controls_INV','heterozygote_controls','negative_controls_STD',
                  'bighead_negative_status','bighead_orthologous_sequence_status',
                  'primer_F','primer_R','amplicon_bp_state_A','amplicon_bp_state_B',
                  'n_carriers_tested','dist_to_bp_bp','notes','source'];
  const csvCell = v => {
    if (v == null) return '';
    if (Array.isArray(v)) v = v.join('; ');
    const s = String(v);
    if (/[",\n]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
    return s;
  };
  const lines = [header.map(csvCell).join(',')];
  for (const m of filtered) {
    const c = m.controls || {};
    const t1Count = (m.af_variants || []).filter(v => v.tier_from_af === 1).length;
    lines.push([m.inversion_id, m.chrom, m.start_bp, m.end_bp, m.tier, m.marker_type,
      m.expected_use, m.validation_status, m.tier_reason,
      m.best_af_score != null ? m.best_af_score.toFixed(3) : null,
      (m.af_variants || []).length, t1Count,
      c.positive_controls_INV, c.heterozygote_controls, c.negative_controls_STD,
      c.bighead_negative_status, c.bighead_orthologous_sequence_status,
      m.primer_F, m.primer_R, m.amplicon_bp_state_A, m.amplicon_bp_state_B,
      m.n_carriers_tested, m.dist_to_bp_bp, m.notes, m.source].map(csvCell).join(','));
  }
  const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'marker_panel.csv';
  document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 250);
}

// ---------------------------------------------------------------------------
// End verbatim extraction.
// ---------------------------------------------------------------------------

// Public entry — alias the underscore-prefixed legacy name to a clean export.
export function renderMarkerPanelPage() {
  return _renderMarkerPanelPage();
}

// Internal helpers re-exported. _mpDeriveAutoPanel is read by sibling
// page17 (stats profile) so it MUST be exported. The merge chat may
// promote it (and _mpScoreVariantAf) to a shared module if other batches
// also need it.
export {
  _mpScoreVariantAf,
  _mpAnnotateGelVisibility,
  _mpSuggestControlsFromKaryotype,
  _mpDefaultCrossSpeciesControls,
  _mpMinDistanceToCsBreakpoint,
  _mpAutoTierFromAtlas,
  _mpDeriveAutoPanel,
  _mpAnnotateAfVariants,
  _mpIsValidPanelJson,
  _mpIsValidVariantAfsJson,
  _mpMergeUserJson,
  _mpParseTsv,
  _mpEnsureState,
  _mpRefreshPanel,
  _mpExportCsv,
  _mpIngestText,
};

export const __MODULE_ID__ = 'inversion_catalogue/page18';
