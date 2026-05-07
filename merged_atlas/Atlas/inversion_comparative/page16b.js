// page16b.js — Multi-species classification cockpit (page16b)
// Extracted verbatim from legacy/Inversion_atlas.html for batch 5 of the parallel migration.
//
// Source line ranges in legacy:
//   26026-28366 — dotplot mashmap IO + multi-species cockpit (init/state, lineage,
//                 architecture/age auto-suggest, classification IO, tree/center/right
//                 column rendering, TE fragility layer, karyotype-lineage layer,
//                 polarize-karyotype-event helpers, refinement chips,
//                 _renderMultiSpeciesPage())
//
// Note: the dotplot_mashmap_v1 IO loader (26026-26124) is co-located here because
// the multi-species page is the primary consumer; page16 reads state.dotplotMashmap
// as an optional overlay only. The merge chat may decide to promote the IO loader
// to shared/ if both pages end up depending on it equally.

// =============================================================================
// TODO_MISSING markers (frequency × symbol — likely legacy line)
// =============================================================================
// All references below are unresolved within this module; the merge chat
// resolves them by either promoting to shared/ (multi-page consumers) or by
// copying locally (single-page consumers).
//
//  TODO_MISSING(_esc)                  ×47   legacy line 13925   HTML escape helper, used everywhere
//  TODO_PROMOTE_TO_SHARED(_esc)        — used by every page, belongs in shared/
//
// External library globals (vendor):
//  TODO_MISSING(window.popgenDotplot)        external — popgen dotplot library (referenced by ms dotplot panel)
//
// Top-level state slots this module reads/writes (handled by shared/state.js):
//   READS: state.crossSpecies (page16-owned), state._crossSpeciesUI.active_id,
//          state._csDotplotPanel*, state._focalVsBg, state.repeatDensity
//   OWNS: state.dotplotMashmap, state.syntenyMultispecies, state.phyloTree,
//         state.dxyPerInversion, state.compTEFragility, state.karyotypeLineage,
//         state._multiSpeciesUI, state._msClassifications, state.classifications
//
// All of these are atlas state slots managed by shared/state.js's SLOT_REGISTRY.
//
// =============================================================================
// EXTERNAL STATE READS (not function calls — flagged for the merge chat to
// confirm the global `state` object exposes these slots):
//   state.crossSpecies               — written by page16.js (_storeCrossSpecies)
//   state.dotplotMashmap             — written by this file (_storeDotplotMashmap)
//   state.syntenyMultispecies        — written by this file (_storeSyntenyMultispecies)
//   state.phyloTree                  — written by this file (_storePhyloTree)
//   state.dxyPerInversion            — written by this file (_storeDxyPerInversion)
//   state.compTEFragility            — written by this file (_storeCompTEFragility)
//   state.karyotypeLineage           — written by this file (_storeKaryotypeLineage)
//   state._multiSpeciesUI            — UI scratch, written by _msInitState
//   state._msClassifications         — written by _msInitClassifications
//   state._crossSpeciesUI.active_id  — read for active breakpoint context
// All of these are atlas state slots managed by shared/state.js's SLOT_REGISTRY.
//
// =============================================================================


// =============================================================================
// PRELUDE — dotplot_mashmap_v1 constants (legacy 26016-26024)
// Re-anchored here from immediately above the extraction window so that the
// IO loader functions (_isDotplotMashmapJSON, _storeDotplotMashmap, ...,
// _clearDotplotMashmap) keep working without TODO_MISSING markers.
// =============================================================================

// -----------------------------------------------------------------------------
// dotplot_mashmap_v1.json loader. Recognized by tool === 'dotplot_mashmap_v1'
// and a resolutions[] array. Stored on state.dotplotMashmap; persisted via
// localStorage like cross_species. Loader is registered in _classifyJSONKind
// + the loadMultipleJSONs dispatch chain.
// -----------------------------------------------------------------------------
const DOTPLOT_MASHMAP_TOOL = 'dotplot_mashmap_v1';
const DOTPLOT_MASHMAP_LS_KEY = 'inversion_atlas.dotplot_mashmap';


// =============================================================================
// BLOCK A — full multi-species cockpit (legacy 26026-28366)
// =============================================================================
function _isDotplotMashmapJSON(data) {
  return !!(
    data &&
    data.tool === DOTPLOT_MASHMAP_TOOL &&
    typeof data.schema_version === 'number' &&
    Array.isArray(data.resolutions)
  );
}

function _storeDotplotMashmap(parsed) {
  if (!_isDotplotMashmapJSON(parsed)) return false;
  // Defensive copy — JSON round-trip ensures no shared references.
  const cleaned = {
    schema_version:       parsed.schema_version,
    tool:                 parsed.tool,
    generated_at:         parsed.generated_at || null,
    species_query:        parsed.species_query  || null,
    species_target:       parsed.species_target || null,
    chrom_lengths_query:  parsed.chrom_lengths_query  || {},
    chrom_lengths_target: parsed.chrom_lengths_target || {},
    resolutions:          JSON.parse(JSON.stringify(parsed.resolutions)),
    loaded_at:            new Date().toISOString(),
  };
  state.dotplotMashmap = cleaned;
  return true;
}

function _persistDotplotMashmap() {
  try {
    if (state.dotplotMashmap) {
      localStorage.setItem(DOTPLOT_MASHMAP_LS_KEY,
                           JSON.stringify(state.dotplotMashmap));
    } else {
      localStorage.removeItem(DOTPLOT_MASHMAP_LS_KEY);
    }
  } catch (_) { /* fail-soft — quota or disabled localStorage */ }
}

function _restoreDotplotMashmap() {
  try {
    const raw = localStorage.getItem(DOTPLOT_MASHMAP_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeDotplotMashmap(parsed);
  } catch (_) { return false; }
}

function _clearDotplotMashmap() {
  state.dotplotMashmap = null;
  state._csDotplotPanel = null;
  state._csDotplotPanelFp = null;
  _persistDotplotMashmap();
}

// =============================================================================
// turn 121: multi-species classification cockpit (page16b)
// -----------------------------------------------------------------------------
// Reads:
//   state.crossSpecies                   active bp from page 16
//   state._crossSpeciesUI.active_id      which bp is active
//   state.syntenyMultispecies            NEW layer (synteny_multispecies_v1)
//   state.phyloTree                      NEW layer (phylo_tree_v1)
//   state.dotplotMashmap                 existing layer (optional)
//
// State on state._multiSpeciesUI = { active_species: 'Cfus' | null }.
//
// Empty-state hierarchy (in order of severity, from "all good" to "nothing"):
//   1. Have synteny_multispecies + cs_breakpoint selected → full cockpit
//   2. Have cs_breakpoint selected, NO multispecies layer    → fallback tree, jump-to-page-16 hint
//   3. No cs_breakpoint selected                              → top-level empty hint
// =============================================================================

const MULTI_SPECIES_UI_LS_KEY = 'inversion_atlas.multiSpeciesUI.v1';
const SYNTENY_MULTISPECIES_LS_KEY = 'inversion_atlas.syntenyMultispecies.v1';
const PHYLO_TREE_LS_KEY = 'inversion_atlas.phyloTree.v1';

// -----------------------------------------------------------------------------
// Default reference species manifest. Used as the fallback when no
// phylo_tree_v1.json is loaded so the page is functional at first open.
// Mirrors catfish-synteny-toolkit/config/species_manifest.tsv but extended
// to the 9-species recommendation that comes out of Quentin's pre-filtering
// work (see prior chat: chromosome breakpoint analysis, 2026-04-23 + the
// 17-species Siluriformes dataset). The leaf order is rough taxonomy:
// Trichomycteridae outgroup, then Bagridae/Ictaluridae/Siluridae deep
// outgroups, then close Clarias context, then the two focal species.
// -----------------------------------------------------------------------------
const _MS_DEFAULT_SPECIES = [
  { id: 'Tros',  label: 'Trichomycterus rosablanca', tier: 'deep_outgroup',     focal: false },
  { id: 'Smer',  label: 'Silurus meridionalis',      tier: 'far_outgroup',      focal: false },
  { id: 'Tfulv', label: 'Tachysurus fulvidraco',     tier: 'far_outgroup',      focal: false },
  { id: 'Ipun',  label: 'Ictalurus punctatus',       tier: 'far_outgroup',      focal: false },
  { id: 'Hwyc',  label: 'Hemibagrus wyckioides',     tier: 'far_outgroup',      focal: false },
  { id: 'Phyp',  label: 'Pangasianodon hypophthalmus',tier: 'mid_outgroup',     focal: false },
  { id: 'Capus', label: 'Channallabes apus',         tier: 'clarias_context',   focal: false },
  { id: 'Cfus',  label: 'Clarias fuscus',            tier: 'clarias_context',   focal: false },
  { id: 'Cmac',  label: 'Clarias macrocephalus',     tier: 'core',              focal: true  },
  { id: 'Cgar',  label: 'Clarias gariepinus',        tier: 'core',              focal: true  },
];

function _msInitState() {
  if (!state._multiSpeciesUI || typeof state._multiSpeciesUI !== 'object') {
    state._multiSpeciesUI = { active_species: null };
    try {
      const raw = localStorage.getItem(MULTI_SPECIES_UI_LS_KEY);
      if (raw) {
        const parsed = JSON.parse(raw);
        if (parsed && typeof parsed === 'object') {
          state._multiSpeciesUI = {
            active_species: typeof parsed.active_species === 'string' ? parsed.active_species : null,
          };
        }
      }
    } catch (_) { /* fail-soft */ }
  }
  if (!state.syntenyMultispecies) state.syntenyMultispecies = null;
  if (!state.phyloTree) state.phyloTree = null;
  return state._multiSpeciesUI;
}

function _msPersistUI() {
  try {
    if (state._multiSpeciesUI) {
      localStorage.setItem(MULTI_SPECIES_UI_LS_KEY,
                           JSON.stringify(state._multiSpeciesUI));
    } else {
      localStorage.removeItem(MULTI_SPECIES_UI_LS_KEY);
    }
  } catch (_) { /* quota or disabled */ }
}

// -----------------------------------------------------------------------------
// JSON layer detectors + loaders
// -----------------------------------------------------------------------------
function _isSyntenyMultispeciesJSON(data) {
  return !!(
    data &&
    (data.tool === 'catfish_synteny_toolkit_v1' ||
     data.tool === 'synteny_multispecies_v1') &&
    typeof data.schema_version === 'number' &&
    Array.isArray(data.species) &&
    Array.isArray(data.synteny_blocks)
  );
}

function _storeSyntenyMultispecies(parsed) {
  if (!_isSyntenyMultispeciesJSON(parsed)) return false;
  state.syntenyMultispecies = {
    schema_version: parsed.schema_version,
    tool: parsed.tool,
    generated_at: parsed.generated_at || null,
    species: JSON.parse(JSON.stringify(parsed.species)),
    synteny_blocks: JSON.parse(JSON.stringify(parsed.synteny_blocks)),
    breakpoints_multilineage: Array.isArray(parsed.breakpoints_multilineage)
      ? JSON.parse(JSON.stringify(parsed.breakpoints_multilineage)) : [],
    loaded_at: new Date().toISOString(),
  };
  return true;
}

function _persistSyntenyMultispecies() {
  try {
    if (state.syntenyMultispecies) {
      localStorage.setItem(SYNTENY_MULTISPECIES_LS_KEY,
                           JSON.stringify(state.syntenyMultispecies));
    } else {
      localStorage.removeItem(SYNTENY_MULTISPECIES_LS_KEY);
    }
  } catch (_) { /* fail-soft */ }
}

function _restoreSyntenyMultispecies() {
  try {
    const raw = localStorage.getItem(SYNTENY_MULTISPECIES_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeSyntenyMultispecies(parsed);
  } catch (_) { return false; }
}

function _clearSyntenyMultispecies() {
  state.syntenyMultispecies = null;
  _persistSyntenyMultispecies();
}

function _isPhyloTreeJSON(data) {
  return !!(
    data &&
    data.tool === 'phylo_tree_v1' &&
    typeof data.schema_version === 'number' &&
    typeof data.newick === 'string' &&
    Array.isArray(data.species_set)
  );
}

function _storePhyloTree(parsed) {
  if (!_isPhyloTreeJSON(parsed)) return false;
  state.phyloTree = {
    schema_version: parsed.schema_version,
    tool: parsed.tool,
    species_set: JSON.parse(JSON.stringify(parsed.species_set)),
    newick: parsed.newick,
    calibration: parsed.calibration || null,
    support_values: parsed.support_values || null,
    provenance: parsed.provenance || null,
    loaded_at: new Date().toISOString(),
  };
  return true;
}

function _persistPhyloTree() {
  try {
    if (state.phyloTree) {
      localStorage.setItem(PHYLO_TREE_LS_KEY, JSON.stringify(state.phyloTree));
    } else {
      localStorage.removeItem(PHYLO_TREE_LS_KEY);
    }
  } catch (_) { /* fail-soft */ }
}

function _restorePhyloTree() {
  try {
    const raw = localStorage.getItem(PHYLO_TREE_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storePhyloTree(parsed);
  } catch (_) { return false; }
}

function _clearPhyloTree() {
  state.phyloTree = null;
  _persistPhyloTree();
}

// -----------------------------------------------------------------------------
// Active breakpoint helper — pulls the page-16 selected breakpoint object.
// Returns null if no breakpoint is selected or no cross-species data loaded.
// -----------------------------------------------------------------------------
function _msGetActiveBreakpoint() {
  if (!state.crossSpecies || !state._crossSpeciesUI) return null;
  const id = state._crossSpeciesUI.active_id;
  if (!id) return null;
  const bps = state.crossSpecies.breakpoints || [];
  for (let i = 0; i < bps.length; i++) {
    if (bps[i] && bps[i].id === id) return bps[i];
  }
  return null;
}

// -----------------------------------------------------------------------------
// Get the currently-effective species list. Order:
//   1. synteny_multispecies.species (if loaded)
//   2. _MS_DEFAULT_SPECIES fallback
// Returned list is shaped uniformly: [{id, label, focal}].
// -----------------------------------------------------------------------------
function _msGetEffectiveSpeciesList() {
  const sm = state.syntenyMultispecies;
  if (sm && Array.isArray(sm.species) && sm.species.length > 0) {
    return sm.species.map(sp => ({
      id:    sp.id || sp.name || '',
      label: sp.label || sp.name || sp.id || '',
      tier:  sp.tier || sp.role || null,
      focal: sp.focal === true || sp.role === 'focal' ||
             sp.role === 'sister' || sp.id === 'Cgar' || sp.id === 'Cmac',
    })).filter(sp => sp.id);
  }
  return _MS_DEFAULT_SPECIES.map(sp => ({ ...sp }));
}

// -----------------------------------------------------------------------------
// Lineage distribution lookup for a breakpoint. Reads
// synteny_multispecies.breakpoints_multilineage when present; otherwise
// returns null.
// -----------------------------------------------------------------------------
function _msGetLineageDistribution(bp) {
  if (!bp) return null;
  const sm = state.syntenyMultispecies;
  if (!sm || !Array.isArray(sm.breakpoints_multilineage)) return null;
  // Match either by exact bp_id or by gar_chr + position window
  const bpId = bp.id;
  const garChr = bp.gar_chr;
  const garMid = bp.gar_pos_mb != null ? bp.gar_pos_mb * 1e6 :
                 (bp.gar_start && bp.gar_end ? (bp.gar_start + bp.gar_end) / 2 : null);
  for (const mlbp of sm.breakpoints_multilineage) {
    if (!mlbp) continue;
    if (mlbp.bp_id && mlbp.bp_id === bpId) return mlbp.lineage_distribution || null;
    if (garChr && garMid != null && mlbp.position) {
      const matchChr = mlbp.position.chrom === garChr ||
                       mlbp.position.chrom === ('C_gar_' + garChr) ||
                       ('C_gar_' + mlbp.position.chrom) === garChr;
      if (matchChr && mlbp.position.pos_bp != null &&
          Math.abs(mlbp.position.pos_bp - garMid) < 100000) {
        return mlbp.lineage_distribution || null;
      }
    }
  }
  return null;
}

// -----------------------------------------------------------------------------
// Architecture-class auto-suggest from lineage_distribution.
// Schema (Quentin's two-layer framework, recorded in atlas page 5):
//   A = simple inversion       (boundary present in Cgar+Cmac, absent elsewhere or unknown)
//   B = synteny-boundary inv.  (breakpoint sits at conserved synteny block edge)
//   C = fusion/fission-assoc.  (sister species shows boundary_present + chrom-context change)
//   D = terminal translocation (boundary near chrom end with reversed attachment)
//   E = recurrent rearr. hotspot (boundary present in >=3 lineages)
//   F = ambiguous / assembly-risk (mixed signal, low data)
// Returns { class, confidence, rationale, n_supporting }.
// -----------------------------------------------------------------------------
function _msAutoSuggestArchitecture(bp, lineageDist) {
  if (!bp) return { class: null, confidence: 'unknown', rationale: 'No active breakpoint.', n_supporting: 0 };
  if (!lineageDist) {
    // No multi-species layer loaded — best we can do is base it on the bp event_type.
    const eventType = bp.event_type_refined || bp.event_type || '';
    if (eventType.indexOf('translocation') >= 0) {
      return {
        class: 'D',
        confidence: 'low',
        rationale: 'Event type is <b>translocation</b>; without multi-species data we cannot ' +
                   'distinguish terminal vs internal. Load <code>synteny_multispecies_v1.json</code> ' +
                   'to refine.',
        n_supporting: 0,
      };
    }
    if (eventType === 'fission' || eventType === 'fusion' || eventType.indexOf('fis') >= 0) {
      return {
        class: 'C',
        confidence: 'low',
        rationale: 'Event type is <b>' + _esc(eventType) + '</b>; provisional Class C ' +
                   '(fusion/fission-associated). Multi-species data would confirm chromosome-context ' +
                   'change in sister species.',
        n_supporting: 0,
      };
    }
    if (eventType === 'inversion') {
      return {
        class: 'A',
        confidence: 'low',
        rationale: 'Event type is <b>inversion</b>; provisional Class A (simple inversion). ' +
                   'Multi-species evidence would re-classify as B (synteny-boundary) or ' +
                   'E (recurrent hotspot) if applicable.',
        n_supporting: 0,
      };
    }
    return {
      class: 'F',
      confidence: 'unknown',
      rationale: 'Event type unclear and no multi-species data loaded. Class assignment deferred.',
      n_supporting: 0,
    };
  }
  // We have lineage_distribution. Count "boundary_present" support across species
  // and check for chrom-context changes.
  const states = Object.values(lineageDist);
  const nPresent = states.filter(s => s === 'boundary_present').length;
  const nAbsent  = states.filter(s => s === 'boundary_absent_internal_to_block').length;
  const nFission = states.filter(s => typeof s === 'string' && s.indexOf('fission') >= 0).length;
  const nFusion  = states.filter(s => typeof s === 'string' && s.indexOf('fusion')  >= 0).length;
  const nUnknown = states.filter(s => !s || s === 'unknown').length;
  const nTotal = states.length;
  // Class E — recurrent rearrangement hotspot — strongest claim, requires 3+ lineages
  if (nPresent >= 3) {
    return {
      class: 'E',
      confidence: 'medium',
      rationale: '<b>Boundary present in ' + nPresent + '/' + nTotal + ' species.</b> ' +
                 'Pattern is consistent with a <b>recurrent rearrangement hotspot</b> — the ' +
                 'breakpoint substrate is reused across the catfish phylogeny. Strongest ' +
                 'mechanism class. Pair with TE density at the homologous region in each ' +
                 'species to support the fragility hypothesis.',
      n_supporting: nPresent,
    };
  }
  // Class C — fusion/fission with sister-species support
  if (nFission + nFusion >= 1 && nPresent >= 1) {
    return {
      class: 'C',
      confidence: 'medium',
      rationale: '<b>' + (nFission + nFusion) + '</b> species shows fusion/fission at this ' +
                 'position; another <b>' + nPresent + '</b> shows the boundary as present. ' +
                 'Class C — <b>fusion/fission-associated</b>. Polarize with the species tree ' +
                 'to determine on which lineage the karyotype event occurred.',
      n_supporting: nPresent + nFission + nFusion,
    };
  }
  // Class B — synteny-boundary present in 2 species (focal + 1)
  if (nPresent === 2) {
    return {
      class: 'B',
      confidence: 'low',
      rationale: '<b>Boundary present in ' + nPresent + '/' + nTotal + ' species.</b> ' +
                 'Class B — <b>synteny-boundary inversion</b>. The breakpoint sits at a ' +
                 'conserved synteny block edge. Distinct from Class A (simple inversion ' +
                 'inside a syntenic block).',
      n_supporting: nPresent,
    };
  }
  // Class A — only the focal species shows the boundary
  if (nPresent <= 1) {
    return {
      class: 'A',
      confidence: 'low',
      rationale: 'Boundary appears only in the focal species. Class A — <b>simple inversion</b> ' +
                 'within a syntenic block, not at a conserved boundary.',
      n_supporting: nPresent,
    };
  }
  // Fallback — ambiguous mixed signal
  return {
    class: 'F',
    confidence: 'unknown',
    rationale: 'Mixed lineage distribution (' + nPresent + ' present, ' + nAbsent + ' absent, ' +
               nFission + ' fission, ' + nFusion + ' fusion, ' + nUnknown + ' unknown). ' +
               'Class F — ambiguous; manual classification recommended.',
    n_supporting: nPresent,
  };
}

// -----------------------------------------------------------------------------
// Phylo tree mini renderer. Vertical cladogram, leaves on the right.
// Pure SVG, no external library. Click-handler on each leaf updates
// state._multiSpeciesUI.active_species and triggers a re-render of the
// detail column.
//
// Two modes:
//   - Newick from state.phyloTree.newick (preferred)
//   - Fallback flat tree from _MS_DEFAULT_SPECIES (each leaf at same depth,
//     just for navigation; not a real phylogeny)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// turn 122: dxy_per_inversion_v1 layer.
// Provides per-candidate dXY (between standard and inverted haplotypes inside
// the inversion) + fold_elevation vs flank, which is the quantitative signal
// for the age-model auto-suggest. Computed by R/Python on cohort BEAGLE
// dosage matrices; atlas just consumes the JSON.
// -----------------------------------------------------------------------------
const DXY_PER_INVERSION_LS_KEY = 'inversion_atlas.dxyPerInversion.v1';

function _isDxyPerInversionJSON(data) {
  return !!(
    data &&
    data.tool === 'dxy_per_inversion_v1' &&
    typeof data.schema_version === 'number' &&
    Array.isArray(data.per_inversion)
  );
}

function _storeDxyPerInversion(parsed) {
  if (!_isDxyPerInversionJSON(parsed)) return false;
  state.dxyPerInversion = {
    schema_version: parsed.schema_version,
    tool: parsed.tool,
    generated_at: parsed.generated_at || null,
    params: parsed.params || null,
    per_inversion: JSON.parse(JSON.stringify(parsed.per_inversion)),
    loaded_at: new Date().toISOString(),
  };
  return true;
}

function _persistDxyPerInversion() {
  try {
    if (state.dxyPerInversion) {
      localStorage.setItem(DXY_PER_INVERSION_LS_KEY,
                           JSON.stringify(state.dxyPerInversion));
    } else {
      localStorage.removeItem(DXY_PER_INVERSION_LS_KEY);
    }
  } catch (_) { /* fail-soft */ }
}

function _restoreDxyPerInversion() {
  try {
    const raw = localStorage.getItem(DXY_PER_INVERSION_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeDxyPerInversion(parsed);
  } catch (_) { return false; }
}

function _clearDxyPerInversion() {
  state.dxyPerInversion = null;
  _persistDxyPerInversion();
}

// Fetch dXY entry for a breakpoint/candidate. Match by candidate_id first,
// then by chrom + position window. Returns null if no match.
function _msGetDxyForBreakpoint(bp) {
  if (!bp) return null;
  const dxy = state.dxyPerInversion;
  if (!dxy || !Array.isArray(dxy.per_inversion)) return null;
  // Atlas cs_breakpoints are usually 1 boundary; the dxy layer is per-INVERSION
  // (covers the whole inverted span). Match: bp.candidate_overlap[0] if
  // populated; fall back to bp.id; fall back to chrom + position window.
  const candIds = Array.isArray(bp.candidate_overlap) ? bp.candidate_overlap : [];
  for (const e of dxy.per_inversion) {
    if (!e) continue;
    if (e.candidate_id && (candIds.indexOf(e.candidate_id) >= 0 || e.candidate_id === bp.id)) {
      return e;
    }
  }
  // Position-based fallback: bp gar_chr/pos within e's chrom+span.
  const garChr = bp.gar_chr;
  const garMid = bp.gar_pos_mb != null ? bp.gar_pos_mb * 1e6 :
                 (bp.gar_start && bp.gar_end ? (bp.gar_start + bp.gar_end) / 2 : null);
  if (!garChr || garMid == null) return null;
  for (const e of dxy.per_inversion) {
    if (!e || !e.chrom) continue;
    const matchChr = e.chrom === garChr ||
                     e.chrom === ('C_gar_' + garChr) ||
                     ('C_gar_' + e.chrom) === garChr;
    if (!matchChr) continue;
    if (e.start_bp != null && e.end_bp != null) {
      if (garMid >= e.start_bp && garMid <= e.end_bp) return e;
    }
  }
  return null;
}

// -----------------------------------------------------------------------------
// turn 122: Age-model auto-suggest.
// Five age tags, from Quentin's two-layer framework:
//   YOUNG-POP            recent within-Cgar polymorphism, low dXY
//   OLD-POLY             old polymorphism, elevated dXY between arrangements
//   OLD-BP-YOUNG-INV     old fragile breakpoint, current inversion is young
//   LINEAGE-KARYO        species-lineage karyotype event (fission/fusion)
//   MULTI-AGE-HOTSPOT    same region reused with multiple rearrangement types
//
// Inputs:
//   bp           page-16 cs_breakpoint object
//   lineageDist  per-species boundary status (or null)
//   dxyEntry     per-inversion dXY entry (or null)
// Returns: { age_model, confidence, rationale, signals }
// -----------------------------------------------------------------------------
function _msAutoSuggestAgeModel(bp, lineageDist, dxyEntry, karyoVerdict) {
  if (!bp) {
    return { age_model: null, confidence: 'unknown',
             rationale: 'No active breakpoint.', signals: {} };
  }
  // Compute the comparative-evidence signals
  const signals = { n_present: 0, n_fission: 0, n_fusion: 0, n_total: 0,
                    fold_elevation: null, dxy_within: null,
                    has_mixed_event_types: false, n_rearrangement_types: 0,
                    karyo_verdict: null };
  if (lineageDist) {
    const states = Object.values(lineageDist);
    signals.n_total   = states.length;
    signals.n_present = states.filter(s => s === 'boundary_present').length;
    signals.n_fission = states.filter(s => typeof s === 'string' && s.indexOf('fission') >= 0).length;
    signals.n_fusion  = states.filter(s => typeof s === 'string' && s.indexOf('fusion')  >= 0).length;
    // mixed event types signal — different species show different rearrangement types
    const types = new Set();
    if (signals.n_present > 0) types.add('boundary_present');
    if (signals.n_fission > 0) types.add('fission');
    if (signals.n_fusion  > 0) types.add('fusion');
    signals.n_rearrangement_types = types.size;
    signals.has_mixed_event_types = types.size >= 2;
  }
  if (dxyEntry) {
    if (typeof dxyEntry.fold_elevation_inside_vs_flank === 'number') {
      signals.fold_elevation = dxyEntry.fold_elevation_inside_vs_flank;
    }
    if (typeof dxyEntry.dxy_within_inversion_ref_vs_inv === 'number') {
      signals.dxy_within = dxyEntry.dxy_within_inversion_ref_vs_inv;
    }
  }
  if (karyoVerdict && karyoVerdict.verdict) {
    signals.karyo_verdict = karyoVerdict.verdict;
  }
  const eventType = bp.event_type_refined || bp.event_type || '';
  const isFissionFusion = eventType === 'fission' || eventType === 'fusion' ||
                          eventType.indexOf('fis') >= 0 || eventType.indexOf('fus') >= 0;

  // RULE 0 — Karyotype-driven LINEAGE-KARYO override (turn 124)
  // Species-agnostic verdict keys (focal_lineage_fission, sister_lineage_fission)
  // are primary; legacy (cgar_/cmac_) are accepted for backward compat (turn 125).
  const lineageVerdicts = ['focal_lineage_fission', 'sister_lineage_fission',
                           'cgar_lineage_fission', 'cmac_lineage_fission',
                           'ancestral_split_both_retained'];
  if (lineageVerdicts.indexOf(signals.karyo_verdict) >= 0) {
    return {
      age_model: 'LINEAGE-KARYO',
      confidence: 'high',
      rationale: 'Karyotype-scale polarization: <b>' +
                 _esc(signals.karyo_verdict.replace(/_/g, ' ')) + '</b>. ' +
                 'The chromosome-context evidence (mashmap 1\u20131 / 1\u20132 ' +
                 'classes across catfish lineages) indicates this is a ' +
                 'species-evolution breakpoint, not just a population polymorphism. ' +
                 'Class <b>LINEAGE-KARYO</b>.',
      signals,
    };
  }
  if (signals.karyo_verdict === 'recurrent_fission_hotspot') {
    return {
      age_model: 'MULTI-AGE-HOTSPOT',
      confidence: 'high',
      rationale: 'Karyotype-scale polarization: <b>recurrent fission hotspot</b>. ' +
                 'Multiple lineages show 1\u20132 mappings with <b>different ' +
                 'target chromosome combinations</b> at this focal-species chromosome \u2014 ' +
                 'consistent with reuse of the same fragile region by independent ' +
                 'rearrangement events. Strongest mechanism class.',
      signals,
    };
  }
  // (no_karyotype_change does NOT block downstream rules; absence of
  // chromosome-scale rearrangement is consistent with intra-chromosomal
  // inversions, which is what most YOUNG-POP / OLD-POLY breakpoints are)

  // RULE 1 — MULTI-AGE-HOTSPOT (highest specificity, strongest claim)
  // Same region reused with multiple rearrangement types
  if (signals.has_mixed_event_types && signals.n_rearrangement_types >= 2 &&
      (signals.n_present + signals.n_fission + signals.n_fusion) >= 3) {
    return {
      age_model: 'MULTI-AGE-HOTSPOT',
      confidence: 'medium',
      rationale: '<b>Mixed rearrangement types</b> across ' +
                 (signals.n_present + signals.n_fission + signals.n_fusion) +
                 ' species (' + signals.n_present + ' boundary present, ' +
                 signals.n_fission + ' fission, ' + signals.n_fusion +
                 ' fusion). The fragile region appears reused by rearrangements ' +
                 'of <b>different ages and types</b> — strongest mechanism class.',
      signals,
    };
  }
  // RULE 2 — LINEAGE-KARYO (fission/fusion event polarized to a lineage)
  if (isFissionFusion && (signals.n_fission + signals.n_fusion) >= 1) {
    return {
      age_model: 'LINEAGE-KARYO',
      confidence: 'medium',
      rationale: 'Event type is <b>' + _esc(eventType) + '</b> and ' +
                 (signals.n_fission + signals.n_fusion) + ' comparative species ' +
                 'show a chromosome-context change. Class <b>LINEAGE-KARYO</b> — ' +
                 'this is mostly a species-evolution breakpoint, not just a ' +
                 'population polymorphism. Tree polarization (STEP_11 in the ' +
                 'synteny toolkit) can refine which lineage carries the event.',
      signals,
    };
  }
  // RULE 3 — OLD-BP-YOUNG-INV (≥3 species share boundary, dXY moderate or low)
  if (signals.n_present >= 3 &&
      (signals.fold_elevation == null || signals.fold_elevation < 1.5)) {
    return {
      age_model: 'OLD-BP-YOUNG-INV',
      confidence: signals.fold_elevation != null ? 'medium' : 'low',
      rationale: 'Breakpoint shared in <b>' + signals.n_present + '/' +
                 signals.n_total + '</b> species' +
                 (signals.fold_elevation != null
                   ? ' but dXY inside the inversion is <b>not strongly elevated</b> (fold = ' +
                     signals.fold_elevation.toFixed(2) + '\u00D7 vs flank).'
                   : '.') +
                 ' The fragile region appears <b>old</b> (reused across catfish ' +
                 'lineages) but the present-day inversion is likely <b>young</b>. ' +
                 'Class <b>OLD-BP-YOUNG-INV</b>.',
      signals,
    };
  }
  // RULE 4 — OLD-POLY (high dXY, shared across deep splits)
  if (signals.fold_elevation != null && signals.fold_elevation >= 1.5 &&
      signals.n_present >= 2) {
    return {
      age_model: 'OLD-POLY',
      confidence: 'medium',
      rationale: '<b>dXY elevated</b> inside the inversion (' +
                 signals.fold_elevation.toFixed(2) + '\u00D7 vs flank) AND ' +
                 'boundary shared in ' + signals.n_present + ' species. The ' +
                 'inversion itself is likely <b>old polymorphism</b> — ' +
                 'standard and inverted haplotypes have been diverging for a ' +
                 'long time. Class <b>OLD-POLY</b>.',
      signals,
    };
  }
  // RULE 5 — YOUNG-POP (low dXY, no comparative evidence)
  if (signals.fold_elevation != null && signals.fold_elevation < 1.2 &&
      signals.n_present <= 1) {
    return {
      age_model: 'YOUNG-POP',
      confidence: 'medium',
      rationale: '<b>dXY low</b> inside the inversion (' +
                 signals.fold_elevation.toFixed(2) + '\u00D7 vs flank) AND ' +
                 'boundary present in only ' + signals.n_present + ' species. ' +
                 'The inversion is likely <b>recent</b>, segregating within ' +
                 'the C. gariepinus population. Class <b>YOUNG-POP</b>.',
      signals,
    };
  }
  // Fallback — uncertain (mixed or insufficient evidence)
  let why = 'Insufficient evidence to commit to an age model.';
  if (signals.fold_elevation == null && !lineageDist) {
    why = 'No dXY layer and no multi-species lineage data loaded. ' +
          'Drop <code>dxy_per_inversion_v1.json</code> + ' +
          '<code>synteny_multispecies_v1.json</code> to populate.';
  } else if (signals.fold_elevation == null) {
    why = 'No dXY layer loaded. Drop <code>dxy_per_inversion_v1.json</code> ' +
          'to add quantitative inversion-age signal.';
  } else if (!lineageDist) {
    why = 'No multi-species lineage data loaded. Drop ' +
          '<code>synteny_multispecies_v1.json</code> to add comparative evidence.';
  }
  return {
    age_model: null,
    confidence: 'unknown',
    rationale: why,
    signals,
  };
}

// -----------------------------------------------------------------------------
// turn 122: per-breakpoint manual classification override.
// Stored on state.classifications keyed by bp_id. Persisted across sessions.
// Override always trumps the auto-suggest when shown in chips/exports.
// -----------------------------------------------------------------------------
const CLASSIFICATIONS_LS_KEY = 'inversion_atlas.classifications.v1';

function _msInitClassifications() {
  if (!state.classifications || typeof state.classifications !== 'object') {
    state.classifications = {};
    try {
      const raw = localStorage.getItem(CLASSIFICATIONS_LS_KEY);
      if (raw) {
        const parsed = JSON.parse(raw);
        if (parsed && typeof parsed === 'object') {
          state.classifications = parsed;
        }
      }
    } catch (_) {}
  }
  return state.classifications;
}

function _msPersistClassifications() {
  try {
    localStorage.setItem(CLASSIFICATIONS_LS_KEY,
                         JSON.stringify(state.classifications || {}));
  } catch (_) {}
}

function _msSetClassification(bpId, fields) {
  if (!bpId) return false;
  _msInitClassifications();
  const prev = state.classifications[bpId] || {};
  state.classifications[bpId] = {
    ...prev,
    ...fields,
    updated_at: new Date().toISOString(),
  };
  _msPersistClassifications();
  return true;
}

function _msGetClassification(bpId) {
  if (!bpId) return null;
  _msInitClassifications();
  return state.classifications[bpId] || null;
}

function _msClearClassification(bpId) {
  if (!bpId) return false;
  _msInitClassifications();
  if (state.classifications[bpId]) {
    delete state.classifications[bpId];
    _msPersistClassifications();
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
// turn 122: Manuscript-export TSV.
// One row per cs_breakpoint with both classification labels (architecture +
// age model) — manual override wins, falls back to auto-suggest. Includes
// supporting evidence counts so reviewers can audit the call.
// -----------------------------------------------------------------------------
function _msBuildClassificationTSV() {
  const cs = state.crossSpecies;
  const bps = (cs && Array.isArray(cs.breakpoints)) ? cs.breakpoints : [];
  const speciesList = _msGetEffectiveSpeciesList();
  const speciesIds = speciesList.map(sp => sp.id);
  const headers = [
    'bp_id', 'gar_chr', 'gar_pos_mb', 'event_type',
    'architecture_class', 'architecture_source', 'architecture_confidence',
    'age_model', 'age_source', 'age_confidence',
    'n_species_total', 'n_species_present', 'n_species_fission_fusion',
    'fold_elevation_dxy', 'dxy_within',
    'classification_notes',
  ];
  // Append per-species lineage status columns (one column per species_id)
  for (const sid of speciesIds) headers.push('lineage_' + sid);
  const lines = [headers.join('\t')];
  for (const bp of bps) {
    if (!bp || !bp.id) continue;
    const lineageDist = _msGetLineageDistribution(bp);
    const dxyEntry = _msGetDxyForBreakpoint(bp);
    const karyoVerdict = _msPolarizeKaryotypeEvent(bp.gar_chr);
    const arch = _msAutoSuggestArchitecture(bp, lineageDist);
    const age  = _msAutoSuggestAgeModel(bp, lineageDist, dxyEntry, karyoVerdict);
    const override = _msGetClassification(bp.id);
    const archClass = (override && override.architecture_class) || arch.class || '';
    const archSrc   = (override && override.architecture_class) ? 'manual' :
                      (arch.class ? 'auto' : '');
    const archConf  = (override && override.architecture_class)
                      ? (override.confidence || 'manual')
                      : arch.confidence;
    const ageModel  = (override && override.age_model) || age.age_model || '';
    const ageSrc    = (override && override.age_model) ? 'manual' :
                      (age.age_model ? 'auto' : '');
    const ageConf   = (override && override.age_model)
                      ? (override.confidence || 'manual')
                      : age.confidence;
    const notes = (override && override.notes) || '';
    const garMb = bp.gar_pos_mb != null ? bp.gar_pos_mb.toFixed(3) :
                  (bp.gar_start ? (bp.gar_start / 1e6).toFixed(3) : '');
    const row = [
      bp.id,
      bp.gar_chr || '',
      garMb,
      bp.event_type_refined || bp.event_type || '',
      archClass, archSrc, archConf,
      ageModel,  ageSrc,  ageConf,
      age.signals.n_total       != null ? age.signals.n_total       : '',
      age.signals.n_present     != null ? age.signals.n_present     : '',
      (age.signals.n_fission || 0) + (age.signals.n_fusion || 0),
      age.signals.fold_elevation != null ? age.signals.fold_elevation.toFixed(3) : '',
      age.signals.dxy_within     != null ? age.signals.dxy_within.toExponential(3) : '',
      notes.replace(/\t/g, ' ').replace(/\n/g, ' '),
    ];
    for (const sid of speciesIds) {
      const st = lineageDist ? (lineageDist[sid] || lineageDist[sid + ''] || '') : '';
      row.push(st);
    }
    lines.push(row.join('\t'));
  }
  return lines.join('\n') + '\n';
}

function _msDownloadClassificationTSV() {
  const tsv = _msBuildClassificationTSV();
  const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  const stamp = new Date().toISOString().replace(/[:.]/g, '-').replace('T', '_').slice(0, 19);
  a.href = url;
  a.download = 'multispecies_classification_' + stamp + '.tsv';
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 200);
}

function _msParseNewickToLeaves(newick, allowedIds) {
  // Minimal Newick parser → flat list of leaves with depth approximation.
  // We don't need full topology yet (only leaf order + relative depth) to
  // render a left-rooted cladogram. Returns [{id, depth}].
  if (!newick || typeof newick !== 'string') return [];
  const allowed = allowedIds ? new Set(allowedIds) : null;
  const leaves = [];
  let depth = 0;
  let i = 0;
  let token = '';
  let curDepth = 0;
  function flushToken() {
    if (token) {
      // Strip optional :branch_length
      const colon = token.indexOf(':');
      const id = (colon >= 0 ? token.substring(0, colon) : token).trim();
      if (id && (!allowed || allowed.has(id))) {
        leaves.push({ id, depth: curDepth });
      } else if (id) {
        // Still record with the "unknown" flag so the fallback still works.
        leaves.push({ id, depth: curDepth, _unknown: true });
      }
      token = '';
    }
  }
  while (i < newick.length) {
    const ch = newick.charAt(i);
    if (ch === '(') { depth++; curDepth = depth; flushToken(); }
    else if (ch === ')') { flushToken(); depth--; curDepth = depth; }
    else if (ch === ',') { flushToken(); }
    else if (ch === ';') { flushToken(); break; }
    else { token += ch; }
    i++;
  }
  flushToken();
  return leaves;
}

function _msRenderTreeSvg() {
  const slot = document.getElementById('msTree');
  if (!slot) return;
  const ui = _msInitState();
  const speciesList = _msGetEffectiveSpeciesList();
  const speciesById = {};
  for (const sp of speciesList) speciesById[sp.id] = sp;
  let leaves;
  let provenance;
  if (state.phyloTree && state.phyloTree.newick) {
    leaves = _msParseNewickToLeaves(state.phyloTree.newick, speciesList.map(s => s.id));
    if (leaves.length === 0) leaves = speciesList.map(s => ({ id: s.id, depth: 0 }));
    provenance = 'phylo_tree_v1.json';
  } else {
    leaves = speciesList.map(s => ({ id: s.id, depth: 0 }));
    provenance = 'reference fallback';
  }
  const N = leaves.length;
  const lineH = 22;
  const padTop = 8, padBot = 8, padL = 12, padR = 14;
  const labelXOffset = 12;  // gap between branch tip and label start
  const treeWidth = 70;     // px reserved for branch lines
  const minLabelWidth = 130;
  const W = padL + treeWidth + labelXOffset + minLabelWidth + padR;
  const H = padTop + padBot + N * lineH;
  const tipX = padL + treeWidth;
  // Build SVG
  const parts = [];
  parts.push('<div class="ms-tree-title">Catfish phylogeny <span style="float:right;color:var(--ink-dimmer)">' +
             '(' + _esc(provenance) + ')</span></div>');
  parts.push('<svg class="ms-tree-svg" viewBox="0 0 ' + W + ' ' + H +
             '" preserveAspectRatio="xMinYMin meet">');
  // Backbone vertical line at the root
  parts.push('<line class="ms-tree-branch" x1="' + padL + '" y1="' + (padTop + lineH/2) +
             '" x2="' + padL + '" y2="' + (padTop + (N - 0.5) * lineH) + '"></line>');
  for (let li = 0; li < N; li++) {
    const leaf = leaves[li];
    const sp = speciesById[leaf.id];
    if (!sp) continue;
    const y = padTop + (li + 0.5) * lineH;
    // Branch line: horizontal from backbone to tip
    parts.push('<line class="ms-tree-branch" x1="' + padL + '" y1="' + y +
               '" x2="' + tipX + '" y2="' + y + '"></line>');
    // Leaf marker (small filled circle at tip)
    const isFocal  = !!sp.focal;
    const isActive = !isFocal && ui.active_species === sp.id;
    const markerCls = ['ms-tree-leaf-marker'];
    if (isFocal) markerCls.push('ms-tree-leaf-focal');
    else if (isActive) markerCls.push('ms-tree-leaf-active');
    parts.push('<circle class="' + markerCls.join(' ') + '" cx="' + tipX + '" cy="' + y +
               '" r="3.5"></circle>');
    // Leaf label
    const labelCls = ['ms-tree-leaf-label'];
    if (isFocal) labelCls.push('ms-tree-leaf-focal');
    else if (isActive) labelCls.push('ms-tree-leaf-active');
    const prefix = isFocal ? '\u2605 ' : (isActive ? '\u25CF ' : '');
    parts.push('<text class="' + labelCls.join(' ') +
               '" data-species-id="' + _esc(sp.id) + '" x="' + (tipX + labelXOffset) +
               '" y="' + (y + 4) + '">' + _esc(prefix + sp.label) + '</text>');
  }
  parts.push('</svg>');
  parts.push('<div class="ms-tree-legend">' +
             '<b>\u2605</b> focal species (Cgar, Cmac) &middot; ' +
             '<b>\u25CF</b> active comparison species &middot; ' +
             'click any non-focal species to compare it to the active breakpoint.' +
             '</div>');
  slot.innerHTML = parts.join('');
  // Wire click handlers (delegate)
  const labels = slot.querySelectorAll('.ms-tree-leaf-label');
  labels.forEach(label => {
    label.addEventListener('click', (ev) => {
      ev.preventDefault();
      const id = label.getAttribute('data-species-id');
      const sp = speciesById[id];
      if (!sp) return;
      if (sp.focal) return;  // focal species are always-active anchors; click is a no-op
      ui.active_species = id;
      _msPersistUI();
      _msRenderTreeSvg();           // re-render to update marker styles
      _msRenderDetailColumn();      // update detail column for new species
    });
  });
}

// -----------------------------------------------------------------------------
// Active-breakpoint header. Reads from state.crossSpecies + UI.
// -----------------------------------------------------------------------------
function _msRenderActiveHeader() {
  const slot = document.getElementById('msActiveHeader');
  if (!slot) return;
  const bp = _msGetActiveBreakpoint();
  if (!bp) {
    slot.innerHTML = '<span class="ms-active-header-empty">' +
                     'No breakpoint selected on page 13 cross-species. ' +
                     '<a href="#" class="ms-jump-link" data-jump="page16">switch to page 13 \u2192</a>' +
                     '</span>';
    const link = slot.querySelector('.ms-jump-link');
    if (link) {
      link.addEventListener('click', (ev) => {
        ev.preventDefault();
        const btn = document.querySelector('#tabBar button[data-page="page16"]');
        if (btn) btn.click();
      });
    }
    return;
  }
  const eventType = bp.event_type_refined || bp.event_type || 'unknown';
  const garChr = bp.gar_chr || '?';
  const garMb = bp.gar_pos_mb != null ? bp.gar_pos_mb.toFixed(2) :
                (bp.gar_start ? (bp.gar_start / 1e6).toFixed(2) : '?');
  slot.innerHTML =
    '<span class="ms-active-bp-id">' + _esc(bp.id) + '</span>' +
    '<span class="ms-active-bp-coord">' +
      _esc(garChr) + ' @ ' + _esc(garMb) + ' Mb &middot; ' + _esc(eventType) +
    '</span>' +
    '<span style="flex:1 1 auto"></span>' +
    '<a href="#" class="ms-jump-link" data-jump="page16">view on page 13 \u2192</a>';
  const link = slot.querySelector('.ms-jump-link');
  if (link) {
    link.addEventListener('click', (ev) => {
      ev.preventDefault();
      const btn = document.querySelector('#tabBar button[data-page="page16"]');
      if (btn) btn.click();
    });
  }
}

// -----------------------------------------------------------------------------
// Center column: classification chips + lineage distribution table.
// -----------------------------------------------------------------------------
function _msRenderCenterColumn() {
  const empty = document.getElementById('msEmpty');
  const cls   = document.getElementById('msClassification');
  const lin   = document.getElementById('msLineageTable');
  if (!empty || !cls || !lin) return;
  const bp = _msGetActiveBreakpoint();
  if (!bp) {
    empty.style.display = 'block';
    cls.style.display = 'none';
    lin.style.display = 'none';
    return;
  }
  empty.style.display = 'none';
  // --- Pull all evidence
  const lineageDist = _msGetLineageDistribution(bp);
  const dxyEntry = _msGetDxyForBreakpoint(bp);
  const karyoVerdict = _msPolarizeKaryotypeEvent(bp.gar_chr);
  const arch = _msAutoSuggestArchitecture(bp, lineageDist);
  const age  = _msAutoSuggestAgeModel(bp, lineageDist, dxyEntry, karyoVerdict);
  const override = _msGetClassification(bp.id) || {};
  const effArchClass = override.architecture_class || arch.class || null;
  const effArchSrc   = override.architecture_class ? 'manual' :
                       (arch.class ? 'auto' : null);
  const effArchConf  = override.architecture_class
                       ? (override.confidence || 'manual')
                       : arch.confidence;
  const effAgeModel  = override.age_model || age.age_model || null;
  const effAgeSrc    = override.age_model ? 'manual' :
                       (age.age_model ? 'auto' : null);
  const effAgeConf   = override.age_model
                       ? (override.confidence || 'manual')
                       : age.confidence;
  // --- Classification block (architecture + age + override editor + export)
  cls.style.display = 'block';
  const archChipHtml = effArchClass
    ? '<span class="ms-cls-chip ms-cls-' + effArchClass + '">Class ' + effArchClass +
        '<span class="ms-cls-confidence">(' + _esc(effArchConf) +
        (effArchSrc === 'manual' ? ', manual' : '') + ')</span>' +
      '</span>'
    : '<span class="ms-cls-empty">No architecture call yet \u2014 ' + _esc(arch.rationale) + '</span>';
  const ageChipHtml = effAgeModel
    ? '<span class="ms-cls-chip ms-age-' + effAgeModel.replace(/[^A-Za-z0-9]/g, '_') + '">' +
        _esc(effAgeModel) +
        '<span class="ms-cls-confidence">(' + _esc(effAgeConf) +
        (effAgeSrc === 'manual' ? ', manual' : '') + ')</span>' +
      '</span>'
    : '<span class="ms-cls-empty">No age model yet</span>';
  // Architecture options
  const archOpts = ['', 'A', 'B', 'C', 'D', 'E', 'F'].map(c =>
    '<option value="' + c + '"' +
    (override.architecture_class === c ? ' selected' : '') +
    '>' + (c ? 'Class ' + c : '\u2014 auto \u2014') + '</option>').join('');
  // Age options
  const ageOpts = ['', 'YOUNG-POP', 'OLD-POLY', 'OLD-BP-YOUNG-INV',
                   'LINEAGE-KARYO', 'MULTI-AGE-HOTSPOT'].map(a =>
    '<option value="' + a + '"' +
    (override.age_model === a ? ' selected' : '') +
    '>' + (a || '\u2014 auto \u2014') + '</option>').join('');
  const confOpts = ['low', 'medium', 'high', 'manual'].map(c =>
    '<option value="' + c + '"' +
    (override.confidence === c ? ' selected' : '') +
    '>' + c + '</option>').join('');
  cls.innerHTML =
    '<div class="ms-cls-block">' +
      '<div class="ms-cls-title">Classification \u2014 architecture + age model</div>' +
      '<div class="ms-cls-chips">' + archChipHtml + ' ' + ageChipHtml + '</div>' +
      '<div class="ms-cls-rationale">' +
        '<b>Architecture rationale.</b> ' + arch.rationale +
      '</div>' +
      '<div class="ms-cls-rationale">' +
        '<b>Age-model rationale.</b> ' + age.rationale +
      '</div>' +
      '<div class="ms-cls-override">' +
        '<div class="ms-cls-override-title">Manual override</div>' +
        '<div class="ms-cls-override-grid">' +
          '<label class="ms-kv-label">Architecture</label>' +
          '<select id="msOverrideArch" class="ms-cls-select">' + archOpts + '</select>' +
          '<label class="ms-kv-label">Age model</label>' +
          '<select id="msOverrideAge" class="ms-cls-select">' + ageOpts + '</select>' +
          '<label class="ms-kv-label">Confidence</label>' +
          '<select id="msOverrideConf" class="ms-cls-select"><option value="">\u2014 default \u2014</option>' +
            confOpts + '</select>' +
          '<label class="ms-kv-label">Notes</label>' +
          '<textarea id="msOverrideNotes" class="ms-cls-textarea" rows="2" ' +
            'placeholder="optional manuscript-style note">' +
            _esc(override.notes || '') + '</textarea>' +
        '</div>' +
        '<div class="ms-cls-override-buttons">' +
          '<button id="msOverrideSaveBtn" class="ms-btn ms-btn-primary">save override</button> ' +
          '<button id="msOverrideClearBtn" class="ms-btn">clear override</button> ' +
          '<span class="ms-cls-override-meta">' +
            (override && override.updated_at
              ? 'last saved: ' + _esc(override.updated_at.replace('T', ' ').slice(0, 19))
              : 'no override saved') +
          '</span>' +
        '</div>' +
      '</div>' +
      '<div class="ms-cls-export">' +
        '<button id="msExportTSVBtn" class="ms-btn">' +
          '\u2B07 export classification TSV \u2014 all breakpoints' +
        '</button>' +
        '<span class="ms-cls-export-hint">one row per cs_breakpoint with both labels' +
          ' &middot; manual override wins, auto-suggest fallback</span>' +
      '</div>' +
    '</div>';
  // Wire override buttons
  const saveBtn = document.getElementById('msOverrideSaveBtn');
  const clearBtn = document.getElementById('msOverrideClearBtn');
  const exportBtn = document.getElementById('msExportTSVBtn');
  if (saveBtn) {
    saveBtn.addEventListener('click', () => {
      const fields = {};
      const a = document.getElementById('msOverrideArch');
      const g = document.getElementById('msOverrideAge');
      const c = document.getElementById('msOverrideConf');
      const n = document.getElementById('msOverrideNotes');
      if (a && a.value) fields.architecture_class = a.value; else fields.architecture_class = '';
      if (g && g.value) fields.age_model = g.value; else fields.age_model = '';
      if (c && c.value) fields.confidence = c.value; else fields.confidence = '';
      if (n) fields.notes = n.value || '';
      // Strip empty fields so they don't shadow the auto-suggest
      const cleaned = {};
      for (const k of Object.keys(fields)) {
        if (fields[k] !== '' && fields[k] != null) cleaned[k] = fields[k];
      }
      _msSetClassification(bp.id, cleaned);
      _msRenderCenterColumn();
    });
  }
  if (clearBtn) {
    clearBtn.addEventListener('click', () => {
      _msClearClassification(bp.id);
      _msRenderCenterColumn();
    });
  }
  if (exportBtn) {
    exportBtn.addEventListener('click', () => {
      try { _msDownloadClassificationTSV(); }
      catch (e) { console.warn('[multiSpecies] TSV export failed:',
                              e && e.message ? e.message : e); }
    });
  }
  // --- Lineage distribution table
  const ui = _msInitState();
  const speciesList = _msGetEffectiveSpeciesList();
  if (lineageDist) {
    const rows = speciesList.map(sp => {
      const st = lineageDist[sp.id] || lineageDist[sp.label] || 'unknown';
      const stCls = (st === 'boundary_present' || st === 'present') ? 'present' :
                    (st === 'boundary_absent_internal_to_block' || st === 'absent') ? 'absent' :
                    (typeof st === 'string' && st.indexOf('fission') >= 0) ? 'fission' :
                    (typeof st === 'string' && st.indexOf('fusion')  >= 0) ? 'fusion' :
                    'unknown';
      const focalStar = sp.focal ? '\u2605 ' : '';
      const isActive = ui.active_species === sp.id;
      const rowCls = isActive ? ' class="ms-lineage-row-active"' : '';
      // turn 124: optional karyotype-class column merge when karyotype data is loaded
      let karyoCol = '';
      if (state.karyotypeLineage) {
        const kentry = _msGetKaryotypeEntryForCgarChr(bp.gar_chr);
        const kcls = kentry && kentry.classes_by_species &&
                     kentry.classes_by_species[sp.id]
                     ? kentry.classes_by_species[sp.id].class : null;
        const kclsCss = kcls === '1-1' ? 'ms-karyo-class-1to1' :
                        kcls === '1-2' ? 'ms-karyo-class-1to2' :
                        kcls === '1-3' ? 'ms-karyo-class-1to3' :
                        kcls === '1-4+' ? 'ms-karyo-class-1to4plus' :
                        'ms-karyo-class-other';
        karyoCol = '<td><span class="ms-karyo-class-chip ' + kclsCss + '">' +
                   _esc(kcls || '\u2014') + '</span></td>';
      }
      return '<tr' + rowCls + '><td>' + focalStar + _esc(sp.label) + '</td>' +
             '<td><span class="ms-lineage-status-' + stCls + '">' + _esc(st) + '</span></td>' +
             karyoCol + '</tr>';
    }).join('');
    const karyoColHeader = state.karyotypeLineage ? '<th>karyo class</th>' : '';
    lin.style.display = 'block';
    lin.innerHTML =
      // turn 124: karyotype context panel sits FIRST — chromosome-scale
      // polarization, separate from breakpoint-scale boundary status.
      _msBuildKaryotypeContextHtml(bp) +
      '<div class="ms-cls-block">' +
        '<div class="ms-cls-title">Lineage distribution (breakpoint scale)' +
          (dxyEntry ? ' &middot; <span class="ms-cls-dxy-note">dXY ' +
            (typeof dxyEntry.fold_elevation_inside_vs_flank === 'number'
              ? dxyEntry.fold_elevation_inside_vs_flank.toFixed(2) + '\u00D7 vs flank'
              : 'loaded') + '</span>'
            : '') +
        '</div>' +
        '<table class="ms-lineage-table">' +
          '<thead><tr><th>Species</th><th>Boundary status</th>' + karyoColHeader + '</tr></thead>' +
          '<tbody>' + rows + '</tbody>' +
        '</table>' +
      '</div>' +
      // turn 123: comparative TE fragility strip — same column, beneath
      // the lineage table. Renders its own empty-state when no
      // comparative_te_breakpoint_fragility_v1.json is loaded.
      _msBuildTEFragilityStripHtml(bp);
  } else {
    lin.style.display = 'block';
    lin.innerHTML =
      // turn 124: even without breakpoint-scale boundary data, the
      // chromosome-scale karyotype context can still be shown.
      _msBuildKaryotypeContextHtml(bp) +
      '<div class="ms-cls-block">' +
        '<div class="ms-cls-title">Lineage distribution (breakpoint scale)</div>' +
        '<div class="ms-cls-empty">' +
        'Per-species boundary status unavailable. Load <code>synteny_multispecies_v1.json</code> ' +
        'to populate this table.' +
        '</div>' +
      '</div>' +
      // turn 123: TE fragility strip is independently meaningful even
      // without lineage_distribution. Show it here too.
      _msBuildTEFragilityStripHtml(bp);
  }
}

// -----------------------------------------------------------------------------
// Detail column: per-species view of how the active species relates to the
// active breakpoint. Reuses popgenDotplot for a 3-row (Cgar / active / Cmac)
// dotplot when synteny_multispecies provides per-species blocks; otherwise
// shows a compact KV summary.
// -----------------------------------------------------------------------------
function _msRenderDetailColumn() {
  const empty   = document.getElementById('msDetailEmpty');
  const header  = document.getElementById('msDetailHeader');
  const dotplot = document.getElementById('msDetailDotplot');
  const bound   = document.getElementById('msDetailBoundary');
  if (!empty || !header || !dotplot || !bound) return;
  const ui = _msInitState();
  const bp = _msGetActiveBreakpoint();
  const activeId = ui.active_species;
  if (!bp || !activeId) {
    empty.style.display = 'block';
    header.style.display = 'none';
    dotplot.style.display = 'none';
    bound.style.display = 'none';
    return;
  }
  const speciesList = _msGetEffectiveSpeciesList();
  const sp = speciesList.find(s => s.id === activeId);
  if (!sp) {
    empty.style.display = 'block';
    header.style.display = 'none';
    dotplot.style.display = 'none';
    bound.style.display = 'none';
    return;
  }
  empty.style.display = 'none';
  header.style.display = 'block';
  header.innerHTML =
    '<div class="ms-detail-title">' + _esc(sp.label) + '</div>' +
    '<div class="ms-detail-subtitle">vs Cgar &middot; vs Cmac at ' +
    _esc(bp.id) + '</div>';
  // --- Dotplot block (only if we have synteny_multispecies blocks for this species)
  const blocks = _msGetSyntenyBlocksForSpecies(sp.id, bp);
  if (blocks && (blocks.cgar_blocks.length > 0 || blocks.cmac_blocks.length > 0) &&
      typeof window !== 'undefined' && window.popgenDotplot &&
      typeof window.popgenDotplot.makeDotplotPanel === 'function') {
    dotplot.style.display = 'block';
    dotplot.innerHTML = '';
    const headerEl = document.createElement('div');
    headerEl.className = 'ms-detail-section-title';
    headerEl.textContent = 'Synteny dotplot vs Cgar';
    dotplot.appendChild(headerEl);
    try {
      const panel = window.popgenDotplot.makeDotplotPanel({
        data: {
          species_query:  { name: 'C. gariepinus', haplotype: 'fClaHyb_Gar_LG' },
          species_target: { name: sp.label, haplotype: sp.id },
          chrom_lengths_query:  blocks.cgar_chrom_lengths || {},
          chrom_lengths_target: blocks.target_chrom_lengths || {},
          wfmash_blocks: blocks.cgar_blocks,
          mashmap_resolutions: null,
        },
        defaultResolution: 'wfmash',
        miniSize: 240,
        enlargedSize: 600,
      });
      dotplot.appendChild(panel);
    } catch (e) {
      console.warn('[multiSpecies] dotplot render failed:',
                   e && e.message ? e.message : e);
      dotplot.style.display = 'none';
    }
  } else {
    dotplot.style.display = 'block';
    dotplot.innerHTML =
      '<div class="ms-detail-section-title">Synteny dotplot vs Cgar</div>' +
      '<div class="ms-detail-empty-msg">' +
      'No synteny blocks loaded for ' + _esc(sp.label) + '. Drop a ' +
      '<code>synteny_multispecies_v1.json</code> with this species included.' +
      '</div>';
  }
  // --- Boundary status (KV summary)
  const lineageDist = _msGetLineageDistribution(bp);
  const status = lineageDist ? (lineageDist[sp.id] || lineageDist[sp.label] || 'unknown') : null;
  bound.style.display = 'block';
  const kvParts = [
    '<div class="ms-detail-section-title">Boundary status at this breakpoint</div>',
  ];
  if (status) {
    const stCls = (status === 'boundary_present' || status === 'present') ? 'present' :
                  (status === 'boundary_absent_internal_to_block' || status === 'absent') ? 'absent' :
                  (typeof status === 'string' && status.indexOf('fission') >= 0) ? 'fission' :
                  (typeof status === 'string' && status.indexOf('fusion')  >= 0) ? 'fusion' :
                  'unknown';
    kvParts.push('<div class="ms-detail-kv">');
    kvParts.push('<span class="ms-kv-label">status</span>');
    kvParts.push('<span class="ms-kv-val ms-lineage-status-' + stCls + '">' +
                 _esc(status) + '</span>');
    if (blocks && blocks.target_chrom) {
      kvParts.push('<span class="ms-kv-label">homologous chrom</span>');
      kvParts.push('<span class="ms-kv-val">' + _esc(blocks.target_chrom) + '</span>');
    }
    if (blocks && blocks.target_orientation) {
      kvParts.push('<span class="ms-kv-label">orientation vs Cgar</span>');
      kvParts.push('<span class="ms-kv-val">' + _esc(blocks.target_orientation) + '</span>');
    }
    kvParts.push('</div>');
  } else {
    kvParts.push('<div class="ms-detail-empty-msg">' +
                 'Boundary status not available. Load <code>synteny_multispecies_v1.json</code> ' +
                 'with a <code>breakpoints_multilineage</code> block covering this breakpoint.' +
                 '</div>');
  }
  bound.innerHTML = kvParts.join('');
  // --- TE-fragility focal-vs-bg (turn 122). Reuses popgenFocalVsBg widget,
  // anchored on the active breakpoint (in Cgar coordinates). The widget
  // surfaces |Z| / theta_pi / fst at the breakpoint focal vs chrom background;
  // when the all_TE class is loaded via state.repeatDensity, it adds an
  // all_TE metric to the widget's metric pool — the comparative fragility
  // signal. The widget renders as empty-state when nothing matches.
  const fvb = document.getElementById('msDetailFocalVsBg');
  if (fvb) {
    if (typeof window !== 'undefined' && window.popgenFocalVsBg &&
        typeof window.popgenFocalVsBg.makeFocalVsBgPanel === 'function' &&
        bp && bp.gar_chr) {
      fvb.style.display = 'block';
      fvb.innerHTML = '';
      const titleEl = document.createElement('div');
      titleEl.className = 'ms-detail-section-title';
      titleEl.textContent = 'TE fragility \u2014 focal vs background';
      fvb.appendChild(titleEl);
      // Initialize the focal-vs-bg radius slot (shared with page 16's widget,
      // so the user's preference travels between pages).
      if (!state._focalVsBg || typeof state._focalVsBg !== 'object') {
        state._focalVsBg = { csRadiusBp: 2000000 };
      }
      const radius = state._focalVsBg.csRadiusBp || 2000000;
      const anchorBp = bp.gar_pos_mb != null ? bp.gar_pos_mb * 1e6 :
                       (bp.gar_start && bp.gar_end ? (bp.gar_start + bp.gar_end) / 2 :
                        (bp.gar_start || bp.gar_end || 0));
      const focal_lo_bp = Math.max(0, anchorBp - radius);
      const focal_hi_bp = anchorBp + radius;
      const anchor_marks_mb = [];
      if (Number.isFinite(bp.gar_pos_start)) anchor_marks_mb.push(bp.gar_pos_start / 1e6);
      if (Number.isFinite(bp.gar_pos_end) && bp.gar_pos_end !== bp.gar_pos_start) {
        anchor_marks_mb.push(bp.gar_pos_end / 1e6);
      }
      try {
        const panel = window.popgenFocalVsBg.makeFocalVsBgPanel({
          page_id:           'page16b',
          anchor_id:         bp.id + ':' + (sp ? sp.id : 'active'),
          anchor_type:       'cs_breakpoint',
          chrom:             bp.gar_chr,
          focal_lo_bp:       focal_lo_bp,
          focal_hi_bp:       focal_hi_bp,
          anchor_marks_mb:   anchor_marks_mb,
          radius_source:     'self',
          radius_default_bp: radius,
          get_radius_bp:     () => state._focalVsBg.csRadiusBp,
          set_radius_bp:     (r) => { state._focalVsBg.csRadiusBp = r; },
          on_radius_change:  () => {
            // Re-render so focal_lo/hi update, but only this column.
            try { _msRenderDetailColumn(); } catch (_) {}
          },
          cohort_selector:   null,
        });
        fvb.appendChild(panel);
      } catch (e) {
        console.warn('[multiSpecies] focal-vs-bg render failed:',
                     e && e.message ? e.message : e);
        fvb.style.display = 'none';
        fvb.innerHTML = '';
      }
    } else {
      fvb.style.display = 'none';
      fvb.innerHTML = '';
    }
  }
}

// -----------------------------------------------------------------------------
// Synteny block lookup for an (active species, active breakpoint) pair.
// Returns { cgar_blocks, cmac_blocks, target_chrom, target_orientation,
// cgar_chrom_lengths, target_chrom_lengths } if synteny_multispecies has
// data; otherwise null.
// -----------------------------------------------------------------------------
function _msGetSyntenyBlocksForSpecies(speciesId, bp) {
  const sm = state.syntenyMultispecies;
  if (!sm || !Array.isArray(sm.synteny_blocks) || !bp) return null;
  const cgarChr = bp.gar_chr;
  const cgarBlocks = [];
  let targetChrom = null;
  let targetOrientation = null;
  const cgarChromLengths = {};
  const targetChromLengths = {};
  for (const block of sm.synteny_blocks) {
    if (!block || !Array.isArray(block.members)) continue;
    let cgarMember = null, targetMember = null;
    for (const mb of block.members) {
      if (!mb) continue;
      if (mb.species === 'Cgar' || mb.species === 'C_gariepinus' ||
          mb.species === 'Clarias gariepinus') cgarMember = mb;
      if (mb.species === speciesId) targetMember = mb;
    }
    if (!cgarMember || !targetMember) continue;
    if (cgarChr && cgarMember.chrom !== cgarChr &&
        cgarMember.chrom !== ('C_gar_' + cgarChr) &&
        ('C_gar_' + cgarMember.chrom) !== cgarChr) continue;
    // Build wfmash-shaped block for popgenDotplot
    cgarBlocks.push({
      gar_chr:    cgarMember.chrom,
      gar_start:  cgarMember.start_bp,
      gar_end:    cgarMember.end_bp,
      mac_chr:    targetMember.chrom,
      mac_start:  targetMember.start_bp,
      mac_end:    targetMember.end_bp,
      strand:     (cgarMember.strand || '+') === (targetMember.strand || '+') ? '+' : '-',
      block_size_bp: Math.max(0, (cgarMember.end_bp || 0) - (cgarMember.start_bp || 0)),
      mapping_quality: block.support && block.support.n_anchors ? block.support.n_anchors : 60,
    });
    if (!targetChrom) targetChrom = targetMember.chrom;
    if (!targetOrientation) {
      targetOrientation = (cgarMember.strand || '+') === (targetMember.strand || '+') ? 'same' : 'inverted';
    }
    if (cgarMember.chrom_length_bp) cgarChromLengths[cgarMember.chrom] = cgarMember.chrom_length_bp;
    if (targetMember.chrom_length_bp) targetChromLengths[targetMember.chrom] = targetMember.chrom_length_bp;
  }
  return {
    cgar_blocks: cgarBlocks,
    cmac_blocks: [],   // Phase 2 — for now we focus on Cgar↔active comparison only
    target_chrom: targetChrom,
    target_orientation: targetOrientation,
    cgar_chrom_lengths: cgarChromLengths,
    target_chrom_lengths: targetChromLengths,
  };
}

// -----------------------------------------------------------------------------
// Top-level page renderer. Called from the tab-dispatcher activation hook.
// -----------------------------------------------------------------------------
function _renderMultiSpeciesPage() {
  _msInitState();
  _msRenderActiveHeader();
  _msRenderTreeSvg();
  _msRenderCenterColumn();
  _msRenderDetailColumn();
  // Hook the help-page jump link (same pattern as page 16)
  const link = document.getElementById('msClassificationFrameworkLink');
  if (link && !link._msWired) {
    link._msWired = true;
    link.addEventListener('click', (ev) => {
      ev.preventDefault();
      const btn = document.querySelector('#tabBar button[data-page="page5"]');
      if (btn) {
        btn.click();
        setTimeout(() => {
          const target = document.getElementById('bp-classification-section');
          if (target && target.scrollIntoView) {
            target.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }
        }, 100);
      }
    });
  }
}

// Expose for tests
if (typeof window !== 'undefined') {
  window._isSyntenyMultispeciesJSON = _isSyntenyMultispeciesJSON;
  window._isPhyloTreeJSON           = _isPhyloTreeJSON;
  window._isDxyPerInversionJSON     = _isDxyPerInversionJSON;
  window._storeSyntenyMultispecies  = _storeSyntenyMultispecies;
  window._storePhyloTree            = _storePhyloTree;
  window._storeDxyPerInversion      = _storeDxyPerInversion;
  window._msAutoSuggestArchitecture = _msAutoSuggestArchitecture;
  window._msAutoSuggestAgeModel     = _msAutoSuggestAgeModel;
  window._msGetEffectiveSpeciesList = _msGetEffectiveSpeciesList;
  window._msParseNewickToLeaves     = _msParseNewickToLeaves;
  window._msSetClassification       = _msSetClassification;
  window._msGetClassification       = _msGetClassification;
  window._msClearClassification     = _msClearClassification;
  window._msBuildClassificationTSV  = _msBuildClassificationTSV;
  window._msDownloadClassificationTSV = _msDownloadClassificationTSV;
  window._renderMultiSpeciesPage    = _renderMultiSpeciesPage;
  window._MS_DEFAULT_SPECIES        = _MS_DEFAULT_SPECIES;
}

// =============================================================================
// turn 123: comparative TE breakpoint fragility (page16b right column + center)
// -----------------------------------------------------------------------------
// Adds a NEW JSON layer comparative_te_breakpoint_fragility_v1 that carries
// per-species TE density at each homologous breakpoint position. Implements
// hypotheses 4 + 5 from chat 1b4d8e12:
//   H4: some Cgar polymorphic inversions reuse ancient synteny boundaries
//   H5: some homologous old-species regions are predicted polymorphic hotspots
//
// The comparative species cannot be tested for polymorphism (no resequencing),
// but their assemblies + EDTA TE annotations are loaded. So we report the
// homologous-region TE density as a fragility proxy and let the user/manuscript
// distinguish "confirmed polymorphic in Gar" from "candidate fragile region in
// species X" using the vocabulary discipline already in the manuscript draft.
// =============================================================================

const TE_FRAGILITY_LS_KEY = 'inversion_atlas.teFragility.v1';

function _isCompTEFragilityJSON(data) {
  return !!(
    data &&
    (data.tool === 'comparative_te_breakpoint_fragility_v1' ||
     data.tool === 'te_fragility_v1') &&
    typeof data.schema_version === 'number' &&
    Array.isArray(data.per_breakpoint_per_species)
  );
}

function _storeCompTEFragility(parsed) {
  if (!_isCompTEFragilityJSON(parsed)) return false;
  state.teFragility = {
    schema_version: parsed.schema_version,
    tool: parsed.tool,
    generated_at: parsed.generated_at || null,
    params: parsed.params ? JSON.parse(JSON.stringify(parsed.params)) : null,
    per_breakpoint_per_species: JSON.parse(JSON.stringify(parsed.per_breakpoint_per_species)),
    loaded_at: new Date().toISOString(),
  };
  return true;
}

function _persistCompTEFragility() {
  try {
    if (state.teFragility) {
      localStorage.setItem(TE_FRAGILITY_LS_KEY, JSON.stringify(state.teFragility));
    } else {
      localStorage.removeItem(TE_FRAGILITY_LS_KEY);
    }
  } catch (_) { /* fail-soft */ }
}

function _restoreCompTEFragility() {
  try {
    const raw = localStorage.getItem(TE_FRAGILITY_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeCompTEFragility(parsed);
  } catch (_) { return false; }
}

function _clearCompTEFragility() {
  state.teFragility = null;
  _persistCompTEFragility();
}

// -----------------------------------------------------------------------------
// Lookup: get TE fragility entry for an (active breakpoint, species) pair.
// Returns { focal_te_density, bg_te_density_chrom, fold_enrichment, percentile,
//           homologous_chrom, focal_window_bp, focal_lo_bp, focal_hi_bp }
// or null when no entry matches.
// Match order: bp_id exact, then gar_chr + gar_pos within 100kb.
// -----------------------------------------------------------------------------
function _msGetTEFragilityForBreakpoint(bp, speciesId) {
  if (!bp || !speciesId) return null;
  const tf = state.teFragility;
  if (!tf || !Array.isArray(tf.per_breakpoint_per_species)) return null;
  const bpId = bp.id;
  const garChr = bp.gar_chr;
  const garMid = bp.gar_pos_mb != null ? bp.gar_pos_mb * 1e6 :
                 (bp.gar_start && bp.gar_end ? (bp.gar_start + bp.gar_end) / 2 : null);
  for (const entry of tf.per_breakpoint_per_species) {
    if (!entry || entry.species !== speciesId) continue;
    if (entry.bp_id && entry.bp_id === bpId) return entry;
    if (garChr && garMid != null && entry.gar_chr) {
      const matchChr = entry.gar_chr === garChr ||
                       entry.gar_chr === ('C_gar_' + garChr) ||
                       ('C_gar_' + entry.gar_chr) === garChr;
      if (matchChr && entry.gar_pos_bp != null &&
          Math.abs(entry.gar_pos_bp - garMid) < 100000) {
        return entry;
      }
    }
  }
  return null;
}

// -----------------------------------------------------------------------------
// Build the per-species comparative TE fragility strip for the active
// breakpoint. Returns an HTML string; renders as a horizontal-ish small
// table: one row per species in the effective species list, columns =
// fold-enrichment chip + percentile + raw density. Empty-state hint when
// no JSON loaded.
// -----------------------------------------------------------------------------
function _msBuildTEFragilityStripHtml(bp) {
  if (!bp) return '';
  const tf = state.teFragility;
  if (!tf) {
    return '<div class="ms-cls-block">' +
             '<div class="ms-cls-title">Comparative TE fragility \u2014 per species</div>' +
             '<div class="ms-cls-empty">' +
             'No <code>comparative_te_breakpoint_fragility_v1.json</code> loaded. ' +
             'Drop one to populate per-species TE density at the homologous region. ' +
             'Pipeline: EDTA annotations \u2192 lift-over to homologous region \u2192 ' +
             'density inside focal radius vs chromosome background.' +
             '</div>' +
           '</div>';
  }
  const speciesList = _msGetEffectiveSpeciesList();
  const rows = [];
  let nWithData = 0;
  for (const sp of speciesList) {
    const entry = _msGetTEFragilityForBreakpoint(bp, sp.id);
    if (!entry) {
      rows.push(
        '<tr>' +
          '<td>' + (sp.focal ? '\u2605 ' : '') + _esc(sp.label) + '</td>' +
          '<td colspan="3" class="ms-tef-empty">no data</td>' +
        '</tr>'
      );
      continue;
    }
    nWithData++;
    const fold = (typeof entry.fold_enrichment === 'number')
                 ? entry.fold_enrichment : null;
    const pct  = (typeof entry.percentile === 'number')
                 ? entry.percentile : null;
    const dens = (typeof entry.focal_te_density === 'number')
                 ? entry.focal_te_density : null;
    // fold-enrichment colour band:
    //   < 1.0 → grey  (depleted)
    //   1.0 - 1.5 → neutral
    //   1.5 - 2.5 → orange (elevated, fragility-suggestive)
    //   > 2.5 → red (strongly elevated, candidate hotspot)
    let foldCls = 'ms-tef-fold-neutral';
    if (fold != null && fold < 1.0) foldCls = 'ms-tef-fold-depleted';
    else if (fold != null && fold >= 2.5) foldCls = 'ms-tef-fold-strong';
    else if (fold != null && fold >= 1.5) foldCls = 'ms-tef-fold-elevated';
    rows.push(
      '<tr>' +
        '<td>' + (sp.focal ? '\u2605 ' : '') + _esc(sp.label) + '</td>' +
        '<td><span class="ms-tef-fold-chip ' + foldCls + '">' +
            (fold != null ? fold.toFixed(2) + '\u00D7' : '\u2014') +
        '</span></td>' +
        '<td>' + (pct != null ? pct.toFixed(0) + '%ile' : '\u2014') + '</td>' +
        '<td>' + (dens != null ? (dens * 100).toFixed(1) + '%' : '\u2014') + '</td>' +
      '</tr>'
    );
  }
  if (nWithData === 0) {
    return '<div class="ms-cls-block">' +
             '<div class="ms-cls-title">Comparative TE fragility \u2014 per species</div>' +
             '<div class="ms-cls-empty">' +
             'TE fragility data is loaded but contains no entry matching this ' +
             'breakpoint (' + _esc(bp.id) + '). Verify the layer covers ' +
             _esc(bp.gar_chr || '?') + '.' +
             '</div>' +
           '</div>';
  }
  return '<div class="ms-cls-block">' +
           '<div class="ms-cls-title">Comparative TE fragility \u2014 per species</div>' +
           '<table class="ms-tef-table">' +
             '<thead><tr>' +
               '<th>Species</th>' +
               '<th>fold vs chrom bg</th>' +
               '<th>%ile</th>' +
               '<th>density</th>' +
             '</tr></thead>' +
             '<tbody>' + rows.join('') + '</tbody>' +
           '</table>' +
           '<div class="ms-tef-legend">' +
           '<span class="ms-tef-fold-chip ms-tef-fold-strong">\u22652.5\u00D7</span> ' +
           'candidate hotspot \u00B7 ' +
           '<span class="ms-tef-fold-chip ms-tef-fold-elevated">1.5\u20132.5\u00D7</span> ' +
           'elevated \u00B7 ' +
           '<span class="ms-tef-fold-chip ms-tef-fold-neutral">~1\u00D7</span> ' +
           'baseline ' +
           '<br><i>fold = TE density in focal radius vs chromosome-wide median. ' +
           'Fragility proxy only \u2014 not evidence of polymorphism in non-focal species.</i>' +
           '</div>' +
         '</div>';
}

// Expose
if (typeof window !== 'undefined') {
  window._isCompTEFragilityJSON       = _isCompTEFragilityJSON;
  window._storeCompTEFragility        = _storeCompTEFragility;
  window._msGetTEFragilityForBreakpoint = _msGetTEFragilityForBreakpoint;
  window._msBuildTEFragilityStripHtml = _msBuildTEFragilityStripHtml;
}

// =============================================================================
// turn 124: karyotype-context lineage layer (page16b chromosome-scale evidence)
// -----------------------------------------------------------------------------
// Source pipeline: mashmap one-to-one (1 Mb / 85% identity) all-vs-all between
// catfish genomes, summarized into per-(query_chr) karyotype class:
//   1-1  → query chrom maps cleanly to a single ref chrom
//   1-2  → query chrom maps to 2 distinct ref chroms (fission OR fusion source,
//          depending on polarity)
//   1-3 / 1-4+ → highly fragmented mapping, often artifacts of deep divergence
//
// The atlas integrates this at chromosome scale (NOT breakpoint scale). The
// breakpoint-scale wfmash refinement remains a separate layer.
//
// The polarization verdict is computed by comparing focal species (Cgar) to
// outgroup(s) and to the rest of the comparative panel. Five outcomes:
//   - cgar_lineage_fission          (Cgar 1-2, others 1-1)
//   - cmac_lineage_fission          (Cmac 1-2, others 1-1)
//   - ancestral_split_both_retained (most species 1-2)
//   - recurrent_fission_hotspot     (multiple species 1-2 with different targets)
//   - no_karyotype_change           (everything 1-1 at this resolution)
//   - unresolved                    (sparse coverage)
//
// HONEST LIMITATION: 1 Mb / 85% identity resolves chromosome-painting only,
// not breakpoint coordinates. Verdicts are at chromosome scale.
// =============================================================================

const KARYO_LINEAGE_LS_KEY = 'inversion_atlas.karyotypeLineage.v1';

function _isKaryotypeLineageJSON(data) {
  return !!(
    data &&
    (data.tool === 'karyotype_lineage_v1' ||
     data.tool === 'mashmap_karyotype_lineage_v1') &&
    typeof data.schema_version === 'number' &&
    typeof data.params === 'object' &&
    // accept either species-agnostic per_focal_chr OR legacy per_cgar_chr
    (Array.isArray(data.per_focal_chr) || Array.isArray(data.per_cgar_chr))
  );
}

function _storeKaryotypeLineage(parsed) {
  if (!_isKaryotypeLineageJSON(parsed)) return false;
  // Normalize to species-agnostic shape on store. The atlas always reads
  // `per_focal_chr` and `focal_chr` after normalization, even if the source
  // JSON used the legacy `per_cgar_chr` / `cgar_chr` field names.
  const rawEntries = Array.isArray(parsed.per_focal_chr)
                      ? parsed.per_focal_chr
                      : (parsed.per_cgar_chr || []);
  const normalizedEntries = rawEntries.map(e => {
    if (!e) return null;
    const copy = JSON.parse(JSON.stringify(e));
    // Promote legacy cgar_chr → focal_chr
    if (copy.focal_chr == null && copy.cgar_chr != null) {
      copy.focal_chr = copy.cgar_chr;
    }
    return copy;
  }).filter(Boolean);
  // Determine focal_species: explicit > inferred from data (the species id
  // that maps 1-1 to its own chrom in every entry)
  let focalSpecies = parsed.focal_species || null;
  if (!focalSpecies && normalizedEntries.length > 0) {
    // Try to infer from first entry: the species whose only target equals focal_chr
    const first = normalizedEntries[0];
    const cls = first.classes_by_species || {};
    for (const sp of Object.keys(cls)) {
      const e = cls[sp];
      if (e && e.class === '1-1' && Array.isArray(e.targets) &&
          e.targets.length === 1 && e.targets[0] === first.focal_chr) {
        focalSpecies = sp;
        break;
      }
    }
  }
  // Sister species: explicit list, or empty (means polarization is binary
  // focal-vs-rest)
  const sisterSpecies = Array.isArray(parsed.sister_species)
                         ? parsed.sister_species.slice() : [];
  const outgroupSpecies = Array.isArray(parsed.outgroup_species)
                          ? parsed.outgroup_species.slice() : [];
  state.karyotypeLineage = {
    schema_version: parsed.schema_version,
    tool: parsed.tool,
    generated_at: parsed.generated_at || null,
    params: JSON.parse(JSON.stringify(parsed.params)),
    focal_species: focalSpecies,
    sister_species: sisterSpecies,
    outgroup_species: outgroupSpecies,
    per_focal_chr: normalizedEntries,
    loaded_at: new Date().toISOString(),
  };
  return true;
}

function _persistKaryotypeLineage() {
  try {
    if (state.karyotypeLineage) {
      localStorage.setItem(KARYO_LINEAGE_LS_KEY,
                           JSON.stringify(state.karyotypeLineage));
    } else {
      localStorage.removeItem(KARYO_LINEAGE_LS_KEY);
    }
  } catch (_) { /* fail-soft */ }
}

function _restoreKaryotypeLineage() {
  try {
    const raw = localStorage.getItem(KARYO_LINEAGE_LS_KEY);
    if (!raw) return false;
    const parsed = JSON.parse(raw);
    return _storeKaryotypeLineage(parsed);
  } catch (_) { return false; }
}

function _clearKaryotypeLineage() {
  state.karyotypeLineage = null;
  _persistKaryotypeLineage();
}

// -----------------------------------------------------------------------------
// Lookup: karyotype entry for a given focal-species chromosome.
// Match priority: exact focal_chr, then prefix-tolerant match using the
// focal_species to derive a 'C_<short>_' style prefix when needed
// (handles 'C_gar_LG28' vs 'LG28' for Cgar focal, etc).
// -----------------------------------------------------------------------------
function _msGetKaryotypeEntryForFocalChr(focalChr) {
  if (!focalChr) return null;
  const kl = state.karyotypeLineage;
  if (!kl || !Array.isArray(kl.per_focal_chr)) return null;
  for (const entry of kl.per_focal_chr) {
    if (!entry || !entry.focal_chr) continue;
    if (entry.focal_chr === focalChr) return entry;
    // Generic prefix tolerance: e.g. 'C_gar_LG28' vs 'LG28'
    // We try stripping/adding common reference prefixes.
    const stripped = entry.focal_chr.replace(/^[A-Z][a-z]?_[a-z]+_/i, '');
    if (stripped === focalChr) return entry;
    if (entry.focal_chr.endsWith('_' + focalChr) ||
        entry.focal_chr.endsWith('|' + focalChr)) return entry;
  }
  return null;
}
// Backward-compat wrapper for the legacy name (used by older callers + tests)
function _msGetKaryotypeEntryForCgarChr(focalChr) {
  return _msGetKaryotypeEntryForFocalChr(focalChr);
}

// -----------------------------------------------------------------------------
// Polarize karyotype evolution for a given focal-species chromosome.
// SPECIES-AGNOSTIC: reads focal_species + sister_species + outgroup_species
// from the loaded JSON, and applies the polarization rules generically.
// Returns { verdict, confidence, rationale, signals, classes_by_species }.
// signals: { n_one_to_one, n_one_to_two, n_other, n_total,
//            focal_class, sister_classes (dict),
//            outgroup_classes (dict), outgroup_one_to_one, outgroup_one_to_two }
// Possible verdicts:
//   - focal_lineage_fission   (focal 1-2, sister + outgroup 1-1)
//   - sister_lineage_fission  (sister 1-2, focal + outgroup 1-1)
//   - ancestral_split_both_retained (≥3 species 1-2 with concordant targets)
//   - recurrent_fission_hotspot (≥3 species 1-2 with discordant targets)
//   - no_karyotype_change (all 1-1)
//   - unresolved (mixed sparse)
// -----------------------------------------------------------------------------
function _msPolarizeKaryotypeEvent(focalChr) {
  const entry = _msGetKaryotypeEntryForFocalChr(focalChr);
  if (!entry) {
    return {
      verdict: null,
      confidence: 'unknown',
      rationale: 'No karyotype lineage data for ' + (focalChr || '?') + '. ' +
                 'Load <code>karyotype_lineage_v1.json</code> to populate.',
      signals: {},
      classes_by_species: {},
    };
  }
  const classes = entry.classes_by_species || {};
  const kl = state.karyotypeLineage;
  const focalSp = (kl && kl.focal_species) || null;
  const sisterSp = (kl && Array.isArray(kl.sister_species)) ? kl.sister_species : [];
  const outgroups = (kl && Array.isArray(kl.outgroup_species)) ? kl.outgroup_species : [];
  const speciesIds = Object.keys(classes);
  const focalLabel = focalSp || 'focal';
  const sisterLabel = sisterSp.length > 0 ? sisterSp.join('/') : 'sister';
  const signals = {
    n_one_to_one: 0,
    n_one_to_two: 0,
    n_other: 0,
    n_total: speciesIds.length,
    focal_species: focalSp,
    focal_class: focalSp && classes[focalSp] ? classes[focalSp].class : null,
    sister_species: sisterSp,
    sister_classes: {},
    outgroup_classes: {},
    outgroup_one_to_one: 0,
    outgroup_one_to_two: 0,
    sister_one_to_two: 0,
    sister_one_to_one: 0,
  };
  for (const sp of speciesIds) {
    const c = classes[sp] && classes[sp].class;
    if (c === '1-1') signals.n_one_to_one++;
    else if (c === '1-2') signals.n_one_to_two++;
    else if (c) signals.n_other++;
    if (outgroups.indexOf(sp) >= 0) {
      signals.outgroup_classes[sp] = c;
      if (c === '1-1') signals.outgroup_one_to_one++;
      else if (c === '1-2') signals.outgroup_one_to_two++;
    }
    if (sisterSp.indexOf(sp) >= 0) {
      signals.sister_classes[sp] = c;
      if (c === '1-1') signals.sister_one_to_one++;
      else if (c === '1-2') signals.sister_one_to_two++;
    }
  }
  // RULE 1 — Recurrent hotspot: ≥3 species 1-2 but with different target sets
  if (signals.n_one_to_two >= 3) {
    const targetSets = [];
    for (const sp of speciesIds) {
      const c = classes[sp];
      if (c && c.class === '1-2' && Array.isArray(c.targets) && c.targets.length === 2) {
        targetSets.push(c.targets.slice().sort().join('|'));
      }
    }
    const uniqueSets = new Set(targetSets);
    if (uniqueSets.size >= 2) {
      return {
        verdict: 'recurrent_fission_hotspot',
        confidence: 'medium',
        rationale: '<b>' + signals.n_one_to_two + '/' + signals.n_total + ' species ' +
                   'show 1\u20132 mapping</b>, but the target chromosome sets <b>differ ' +
                   'between species</b> (' + uniqueSets.size + ' distinct target combinations). ' +
                   'This is consistent with a <b>recurrent fission/fusion hotspot</b> \u2014 ' +
                   'the same region in <b>' + _esc(focalLabel) + '</b> has been involved in ' +
                   'independent rearrangement events across multiple lineages.',
        signals,
        classes_by_species: classes,
      };
    }
    return {
      verdict: 'ancestral_split_both_retained',
      confidence: 'medium',
      rationale: '<b>' + signals.n_one_to_two + '/' + signals.n_total + ' species ' +
                 'share the 1\u20132 mapping with concordant targets</b>. The split is ' +
                 'consistent with an <b>ancestral fission</b> retained across most ' +
                 'lineages; <b>' + _esc(focalLabel) + '</b> (and any other 1\u20131 species) ' +
                 'carries a derived fusion.',
      signals,
      classes_by_species: classes,
    };
  }
  // RULE 2 — Focal-lineage event (focal 1-2, outgroup 1-1, most others 1-1)
  if (signals.focal_class === '1-2' &&
      signals.outgroup_one_to_one >= 1 &&
      signals.outgroup_one_to_two === 0 &&
      signals.n_one_to_one >= 2) {
    return {
      verdict: 'focal_lineage_fission',
      confidence: 'medium',
      rationale: '<b>' + _esc(focalLabel) + ' shows 1\u20132 mapping</b> while ' +
                 signals.outgroup_one_to_one + ' outgroup species and ' +
                 signals.n_one_to_one + ' total species are 1\u20131. The ' +
                 'rearrangement appears to have occurred on the <b>' +
                 _esc(focalLabel) + ' lineage</b> after divergence from the common ancestor.',
      signals,
      classes_by_species: classes,
    };
  }
  // RULE 3 — Sister-lineage event (sister 1-2, focal 1-1, outgroup 1-1)
  if (signals.sister_one_to_two >= 1 &&
      signals.focal_class === '1-1' &&
      signals.outgroup_one_to_one >= 1 &&
      signals.outgroup_one_to_two === 0) {
    // Identify which sister(s) carry the event
    const sisterCarriers = [];
    for (const sp of sisterSp) {
      const c = classes[sp] && classes[sp].class;
      if (c === '1-2') sisterCarriers.push(sp);
    }
    return {
      verdict: 'sister_lineage_fission',
      confidence: 'medium',
      rationale: '<b>' + sisterCarriers.map(_esc).join(', ') + ' shows 1\u20132 mapping</b> ' +
                 'while ' + _esc(focalLabel) + ' and ' + signals.outgroup_one_to_one +
                 ' outgroup species are 1\u20131. The rearrangement appears to have ' +
                 'occurred on the <b>' + sisterCarriers.map(_esc).join('/') + ' lineage</b> ' +
                 'after divergence from the common ancestor.',
      signals,
      classes_by_species: classes,
    };
  }
  // RULE 4 — All 1-1: no detectable karyotype change at this resolution
  if (signals.n_one_to_one === signals.n_total && signals.n_total >= 3) {
    return {
      verdict: 'no_karyotype_change',
      confidence: 'medium',
      rationale: 'All ' + signals.n_total + ' species show <b>1\u20131 mapping</b> for ' +
                 'this <b>' + _esc(focalLabel) + '</b> chromosome. ' +
                 '<b>No fusion/fission detected</b> at the 1 Mb / 85% identity scale. ' +
                 'Any breakpoint here is likely an <b>intra-chromosomal inversion</b>, ' +
                 'not a karyotype-level event.',
      signals,
      classes_by_species: classes,
    };
  }
  // FALLBACK — sparse / mixed signal
  return {
    verdict: 'unresolved',
    confidence: 'low',
    rationale: 'Karyotype signal is mixed: ' + signals.n_one_to_one + ' species 1\u20131, ' +
               signals.n_one_to_two + ' species 1\u20132, ' + signals.n_other + ' other. ' +
               'Insufficient evidence for a confident polarization verdict.',
    signals,
    classes_by_species: classes,
  };
}

// -----------------------------------------------------------------------------
// turn 126: wfmash refinement helpers
// -----------------------------------------------------------------------------
// Each classes_by_species[species_id] entry can carry these optional fields
// from a higher-resolution wfmash pass:
//   refined_by_wfmash: 'confirmed' | 'refuted' | 'refined' | 'failed' | null
//   wfmash_class:      the class string after refinement (when present)
//   wfmash_targets:    array of refined targets with coordinates
//   wfmash_n_blocks:   number of wfmash blocks supporting the call
//   wfmash_pct_identity: median identity across blocks
//
// SEMANTICS
//   confirmed → wfmash agrees with mashmap at higher resolution
//   refuted   → wfmash collapses 1-2 to 1-1 (or otherwise reduces fragmentation)
//   refined   → wfmash gives a *different* class than mashmap (e.g. mashmap 1-2,
//               wfmash splits to 1-3, OR mashmap 1-1, wfmash splits to 1-2)
//   failed    → wfmash attempted but produced no confident alignment
//   null      → not yet attempted
// -----------------------------------------------------------------------------

function _msGetEffectiveClassForCell(cell) {
  // Return the effective class — wfmash overrides mashmap when refined or refuted.
  if (!cell) return null;
  if (cell.refined_by_wfmash === 'refuted' || cell.refined_by_wfmash === 'refined') {
    if (cell.wfmash_class) return cell.wfmash_class;
  }
  return cell.class || null;
}

function _msGetEffectiveTargetsForCell(cell) {
  // Return the effective targets — wfmash refined targets override mashmap when present.
  if (!cell) return [];
  if ((cell.refined_by_wfmash === 'refuted' || cell.refined_by_wfmash === 'refined') &&
      Array.isArray(cell.wfmash_targets) && cell.wfmash_targets.length > 0) {
    return cell.wfmash_targets.map(t =>
      typeof t === 'string' ? t :
      (t && t.chrom ? t.chrom : ''));
  }
  return Array.isArray(cell.targets) ? cell.targets.slice() : [];
}

function _msSummarizeRefinement(entry) {
  // Across all species cells in a per-focal-chr entry, count refinement states.
  // Returns { confirmed, refuted, refined, failed, not_attempted, total }.
  const out = { confirmed: 0, refuted: 0, refined: 0, failed: 0,
                not_attempted: 0, total: 0 };
  if (!entry || !entry.classes_by_species) return out;
  const cls = entry.classes_by_species;
  for (const sp of Object.keys(cls)) {
    out.total++;
    const r = cls[sp] && cls[sp].refined_by_wfmash;
    if (r === 'confirmed') out.confirmed++;
    else if (r === 'refuted')  out.refuted++;
    else if (r === 'refined')  out.refined++;
    else if (r === 'failed')   out.failed++;
    else                       out.not_attempted++;
  }
  return out;
}

function _msAdjustConfidenceForRefinement(verdict, refinementSummary) {
  // Pure function that returns a possibly-adjusted verdict based on refinement.
  // Rules:
  //   - If verdict is null/unresolved → no change (refinement doesn't synthesize new info).
  //   - If ≥80% of cells confirmed and 0 refuted → boost confidence one tier (medium → high).
  //   - If ≥30% of cells refuted → drop one tier and append refutation hint.
  //   - Otherwise unchanged.
  if (!verdict || !verdict.verdict) return verdict;
  if (verdict.verdict === 'unresolved') return verdict;
  const s = refinementSummary;
  if (!s || s.total === 0) return verdict;
  const fracConfirmed = s.confirmed / s.total;
  const fracRefuted   = s.refuted / s.total;
  const fracTouched   = (s.confirmed + s.refuted + s.refined + s.failed) / s.total;
  const order = ['low', 'medium', 'high'];
  function boost(c, n) {
    const i = order.indexOf(c);
    if (i < 0) return c;
    const j = Math.min(order.length - 1, Math.max(0, i + n));
    return order[j];
  }
  // Need at least 50% of cells to have been touched by wfmash before we
  // adjust confidence based on refinement (otherwise the sample is too small).
  if (fracTouched < 0.5) return verdict;
  if (fracConfirmed >= 0.8 && s.refuted === 0) {
    return Object.assign({}, verdict, {
      confidence: boost(verdict.confidence, +1),
      rationale: verdict.rationale +
                 ' <span class="ms-refine-note">Boosted by wfmash confirmation across ' +
                 s.confirmed + '/' + s.total + ' species cells.</span>',
      _refinement: s,
    });
  }
  if (fracRefuted >= 0.3) {
    return Object.assign({}, verdict, {
      confidence: boost(verdict.confidence, -1),
      rationale: verdict.rationale +
                 ' <span class="ms-refine-note ms-refine-note-warn">Caution: wfmash refuted ' +
                 s.refuted + '/' + s.total + ' mashmap cells at higher resolution.</span>',
      _refinement: s,
    });
  }
  return Object.assign({}, verdict, { _refinement: s });
}

// -----------------------------------------------------------------------------
// Build the small refinement chip HTML for one cell.
// -----------------------------------------------------------------------------
function _msBuildRefinementChipHtml(cell) {
  if (!cell) return '';
  const r = cell.refined_by_wfmash;
  if (!r) return '<span class="ms-refine-chip ms-refine-chip-pending" title="not yet attempted">\u2014</span>';
  if (r === 'confirmed') {
    return '<span class="ms-refine-chip ms-refine-chip-confirmed" title="wfmash confirmed">' +
           '\u2713 wfmash</span>';
  }
  if (r === 'refuted') {
    return '<span class="ms-refine-chip ms-refine-chip-refuted" ' +
           'title="wfmash refuted at higher resolution">' +
           '\u2717 refuted' +
           (cell.wfmash_class ? ' \u2192 ' + _esc(cell.wfmash_class) : '') +
           '</span>';
  }
  if (r === 'refined') {
    return '<span class="ms-refine-chip ms-refine-chip-refined" ' +
           'title="wfmash refined to a different class">' +
           '\u279C ' + _esc(cell.wfmash_class || 'refined') + '</span>';
  }
  if (r === 'failed') {
    return '<span class="ms-refine-chip ms-refine-chip-failed" ' +
           'title="wfmash ran but produced no confident alignment">' +
           '\u26A0 failed</span>';
  }
  return '<span class="ms-refine-chip ms-refine-chip-pending">' + _esc(r) + '</span>';
}

// Polarize-with-refinement convenience: applies refinement adjustment after polarization.
function _msPolarizeKaryotypeEventWithRefinement(focalChr) {
  const verdict = _msPolarizeKaryotypeEvent(focalChr);
  if (!verdict || !verdict.verdict) return verdict;
  const entry = _msGetKaryotypeEntryForFocalChr(focalChr);
  const summary = _msSummarizeRefinement(entry);
  return _msAdjustConfidenceForRefinement(verdict, summary);
}

// -----------------------------------------------------------------------------
// Render the karyotype context panel HTML for a given breakpoint.
// Returns an HTML string.
// -----------------------------------------------------------------------------
function _msBuildKaryotypeContextHtml(bp) {
  if (!bp) return '';
  // bp.gar_chr is the focal-species chromosome the breakpoint is anchored on.
  // Even when the focal species is renamed (e.g. Cmac for the Cmac paper),
  // the page-16 cs_breakpoints schema uses gar_chr as the convention key —
  // so we read from bp.gar_chr regardless of which species is "focal".
  const focalChr = bp.gar_chr;
  if (!focalChr) return '';
  const kl = state.karyotypeLineage;
  if (!kl) {
    return '<div class="ms-cls-block">' +
             '<div class="ms-cls-title">Karyotype context (chromosome scale) \u2014 ' +
             _esc(focalChr) + '</div>' +
             '<div class="ms-cls-empty">' +
             'No <code>karyotype_lineage_v1.json</code> loaded. This layer comes from a ' +
             'mashmap one-to-one (1 Mb / 85% identity) all-vs-all scan that classifies ' +
             'each focal-species chromosome\u2019s mapping to each comparative species ' +
             '(1\u20131 / 1\u20132 / 1\u20133 / 1\u20134+). It polarizes ' +
             'fusion/fission events to specific lineages using outgroup species.' +
             '</div>' +
           '</div>';
  }
  const verdict = _msPolarizeKaryotypeEventWithRefinement(focalChr);
  if (!verdict.verdict) {
    return '<div class="ms-cls-block">' +
             '<div class="ms-cls-title">Karyotype context (chromosome scale) \u2014 ' +
             _esc(focalChr) + '</div>' +
             '<div class="ms-cls-empty">' + verdict.rationale + '</div>' +
           '</div>';
  }
  // Verdict chip color. Species-agnostic verdict keys are primary; legacy
  // cgar_/cmac_ keys are backward-compat aliases.
  const verdictColors = {
    focal_lineage_fission:  'ms-karyo-verdict-focal',
    sister_lineage_fission: 'ms-karyo-verdict-sister',
    cgar_lineage_fission:   'ms-karyo-verdict-focal',   // legacy alias
    cmac_lineage_fission:   'ms-karyo-verdict-sister',  // legacy alias
    ancestral_split_both_retained: 'ms-karyo-verdict-ancestral',
    recurrent_fission_hotspot:     'ms-karyo-verdict-hotspot',
    no_karyotype_change:           'ms-karyo-verdict-stable',
    unresolved:                    'ms-karyo-verdict-unresolved',
  };
  const verdictColorCls = verdictColors[verdict.verdict] || 'ms-karyo-verdict-unresolved';
  // Build the verdict label: the species-agnostic verdict is generic, so
  // substitute the focal species name into the label so the user sees
  // "Cgar lineage fission" rather than "focal lineage fission".
  let verdictLabel = verdict.verdict.replace(/_/g, ' ');
  const focalSp = kl.focal_species;
  if (focalSp && verdict.verdict === 'focal_lineage_fission') {
    verdictLabel = _esc(focalSp) + ' lineage fission';
  }
  if (verdict.verdict === 'sister_lineage_fission' &&
      Array.isArray(verdict.signals && verdict.signals.sister_species) &&
      verdict.signals.sister_species.length > 0) {
    // Pick the sister(s) that actually carry the event
    const sisterCarriers = [];
    const cls = verdict.classes_by_species || {};
    for (const sp of verdict.signals.sister_species) {
      if (cls[sp] && cls[sp].class === '1-2') sisterCarriers.push(sp);
    }
    if (sisterCarriers.length > 0) {
      verdictLabel = sisterCarriers.map(_esc).join('/') + ' lineage fission';
    }
  }
  // Per-species table rows
  const speciesList = _msGetEffectiveSpeciesList();
  const classes = verdict.classes_by_species || {};
  const outgroups = (kl && Array.isArray(kl.outgroup_species)) ? kl.outgroup_species : [];
  // turn 126: detect if this entry has any wfmash refinement at all so we
  // can decide whether to render the refinement column header.
  const refSummary = verdict._refinement || _msSummarizeRefinement(
    _msGetKaryotypeEntryForFocalChr(focalChr));
  const hasAnyRefinement = refSummary &&
    (refSummary.confirmed + refSummary.refuted +
     refSummary.refined + refSummary.failed) > 0;
  const rows = speciesList.map(sp => {
    const c = classes[sp.id];
    const mashmapCls = c && c.class ? c.class : '\u2014';
    const targets = c && Array.isArray(c.targets) ? c.targets.join(', ') : '';
    const focalStar = sp.focal ? '\u2605 ' : '';
    const isOutgroup = outgroups.indexOf(sp.id) >= 0;
    const outgroupTag = isOutgroup ? ' <span class="ms-karyo-outgroup-tag">outgroup</span>' : '';
    let classCls = 'ms-karyo-class-other';
    if (mashmapCls === '1-1') classCls = 'ms-karyo-class-1to1';
    else if (mashmapCls === '1-2') classCls = 'ms-karyo-class-1to2';
    else if (mashmapCls === '1-3') classCls = 'ms-karyo-class-1to3';
    else if (mashmapCls === '1-4+') classCls = 'ms-karyo-class-1to4plus';
    // turn 126: when wfmash refutes the mashmap call, strike through it
    const isRefuted = c && c.refined_by_wfmash === 'refuted';
    if (isRefuted) classCls += ' ms-karyo-class-refuted';
    const refineCol = hasAnyRefinement
      ? '<td class="ms-karyo-refine-col">' + _msBuildRefinementChipHtml(c) + '</td>'
      : '';
    return '<tr>' +
             '<td>' + focalStar + _esc(sp.label) + outgroupTag + '</td>' +
             '<td><span class="ms-karyo-class-chip ' + classCls + '">' +
                 _esc(mashmapCls) + '</span></td>' +
             '<td class="ms-karyo-targets">' + _esc(targets || '\u2014') + '</td>' +
             refineCol +
           '</tr>';
  }).join('');
  const refineHeader = hasAnyRefinement ? '<th>wfmash refinement</th>' : '';
  // Refinement summary line above the table
  let refineSummaryLine = '';
  if (hasAnyRefinement) {
    const parts = [];
    if (refSummary.confirmed) parts.push('<b>' + refSummary.confirmed + '</b> confirmed');
    if (refSummary.refuted)   parts.push('<b>' + refSummary.refuted + '</b> refuted');
    if (refSummary.refined)   parts.push('<b>' + refSummary.refined + '</b> refined');
    if (refSummary.failed)    parts.push('<b>' + refSummary.failed + '</b> failed');
    if (refSummary.not_attempted)
      parts.push('<b>' + refSummary.not_attempted + '</b> not yet attempted');
    refineSummaryLine =
      '<div class="ms-karyo-refine-summary">' +
        '<span class="ms-karyo-refine-summary-label">wfmash refinement:</span> ' +
        parts.join(' &middot; ') +
      '</div>';
  }
  return '<div class="ms-cls-block">' +
           '<div class="ms-cls-title">Karyotype context (chromosome scale) \u2014 ' +
           _esc(focalChr) + '</div>' +
           '<div class="ms-karyo-verdict-row">' +
             '<span class="ms-karyo-verdict-chip ' + verdictColorCls + '">' +
                 _esc(verdictLabel) +
                 '<span class="ms-cls-confidence">(' + _esc(verdict.confidence) +
                 ' confidence)</span>' +
             '</span>' +
           '</div>' +
           '<div class="ms-cls-rationale">' + verdict.rationale + '</div>' +
           refineSummaryLine +
           '<table class="ms-karyo-table">' +
             '<thead><tr><th>Species</th><th>class</th><th>targets</th>' +
               refineHeader + '</tr></thead>' +
             '<tbody>' + rows + '</tbody>' +
           '</table>' +
           '<div class="ms-karyo-legend">' +
             '<i>Mashmap (1 Mb / 85% id) gives chromosome-scale class. ' +
             'wfmash (kbp / higher id) refines: \u2713 confirmed, \u2717 refuted, ' +
             '\u279C refined to different class. ' +
             'Pair with breakpoint-scale wfmash synteny layer for sub-Mb resolution.</i>' +
           '</div>' +
         '</div>';
}

// Expose
if (typeof window !== 'undefined') {
  window._isKaryotypeLineageJSON      = _isKaryotypeLineageJSON;
  window._storeKaryotypeLineage       = _storeKaryotypeLineage;
  window._msGetKaryotypeEntryForFocalChr = _msGetKaryotypeEntryForFocalChr;
  window._msGetKaryotypeEntryForCgarChr = _msGetKaryotypeEntryForCgarChr;
  window._msPolarizeKaryotypeEvent    = _msPolarizeKaryotypeEvent;
  window._msPolarizeKaryotypeEventWithRefinement = _msPolarizeKaryotypeEventWithRefinement;
  window._msSummarizeRefinement       = _msSummarizeRefinement;
  window._msAdjustConfidenceForRefinement = _msAdjustConfidenceForRefinement;
  window._msGetEffectiveClassForCell  = _msGetEffectiveClassForCell;
  window._msGetEffectiveTargetsForCell= _msGetEffectiveTargetsForCell;
  window._msBuildRefinementChipHtml   = _msBuildRefinementChipHtml;
  window._msBuildKaryotypeContextHtml = _msBuildKaryotypeContextHtml;
}

