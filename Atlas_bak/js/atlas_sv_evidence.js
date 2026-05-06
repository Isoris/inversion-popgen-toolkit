// =============================================================================
// atlas_sv_evidence.js — SV Evidence & Boundaries page
// =============================================================================
//
// Per-candidate page surfacing every SV call (DELLY2 + Manta) within ±500 kb
// of the candidate's left and right boundary, scored against karyotype-group
// genotype counts (H1/H1, H1/H2, H2/H2), classified into pattern labels.
//
// Spec: specs_todo/SPEC_sv_evidence_page.md
// Mockup: specs_todo/_mockups/sv_evidence_page_mockup_1_canonical.png
//
// Implementation order (from spec §7):
//   step 1   — skeleton + tab + empty-state mount         [DONE]
//   step 2   — loader + main SV table (sort/filter/etc)   [DONE]
//   step 3   — locus track strip                          [DONE]
//   step 3.5 — cursor + hover + zoom + E/F/Enter/Esc      [DONE]
//   step 4   — right-rail boundary summary + legend       [DONE]
//              + table row annotations (L / R / ★)
//              + ↑/↓ table-row navigation
//              + Shift+↑/↓ now cycles SV-type filter
//              + interactive legend (FDR + pattern → filter)
//              + boundary-summary row click → filter table + locus
//   step 4.5 — region-select mode (drag → filter region)  [DONE]
//              toolbar Mode toggle: zoom | select
//              drag in select mode → translucent region rect
//              region applied as a filter to table + locus
//              Esc clears it; "clear" pill in readout
//   step 5   — UpSet panel (the redirect)                 [DONE]
//              rows = per-SV evidence types (left_SA, right_SA, …)
//              top bars = fish counts per evidence combination
//              click bar → selectedSamples populated; readout banner
//              folder-walk drag-drop via webkitGetAsEntry for per-cand folders
//              new sibling layer: sv_evidence_combinations_v1.json
//              see specs_todo/SPEC_sv_evidence_page__upset_redirect.md
//   step 6   — sample × SV dosage heatmap                  [DONE — this commit]
//              new sibling layer: sv_support_by_sample_v1.json
//              right-rail compact + main-area large overlay
//              cell hover → tooltip; cell click → highlight SV
//              also unlocks: dim non-selected SV glyphs when an UpSet
//              bar is active and supportLayer is loaded
//   step 7   — STEP_SV_GT_AGG + STEP_SV_EVID_COMB + STEP_SV_SUPPORT
//              Python emitters + per-candidate folder writer       [DONE]
//              see producers/ directory; STUBS for VCF parsing —
//              JSON shapes locked, ready to wire to cyvcf2/pysam
//
// Read-only page — does not edit boundaries (page11), karyotype groups
// (scrubber), or SV calls (MODULE_5A2).
//
// API (window.AtlasSVEvidence):
//   .init(opts)          — called once at atlas startup; mounts empty shell
//   .loadCandidate(cid)  — fetches sv_genotype_counts/<cid>.json, renders
//   .refresh()           — re-renders active candidate using state.* updates
//   .setFilters(filters) — programmatic filter override (URL-state)
//   .exportFilteredTSV() — emits filtered+sorted rows as TSV
//   ._state              — for debugging
//
// Lifecycle: tab-bar click handler in Inversion_atlas.html dispatches
//   `if (target === 'page_sv_evidence')` → AtlasSVEvidence.loadCandidate(cid).
// Scrubber's `karyotype_locks_changed` event triggers .refresh().
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.AtlasSVEvidence = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const LS_PREFIX = 'inversion_atlas.sv_evidence.';
  const LS_FILTERS = LS_PREFIX + 'filters'; // candId → filter state
  const LS_ANNOTATIONS = LS_PREFIX + 'annotations'; // candId → { sv_id: 'L'|'R'|'star' }

  // Default filter state (spec §3.3)
  const DEFAULT_FILTERS = {
    sv_type: 'All',           // All | BND | INV | DEL | DUP | Other
    quality: 'PASS',          // PASS | All
    zone: 'boundary_500kb',   // all | boundary_500kb | inversion_body | flanks
    min_samples: 5,
    show_only_associated: true,
    fdr_threshold: 0.05,      // 1e-6 | 1e-4 | 0.01 | 0.05 | 0.10
  };

  // FDR colour palette (spec §5.2)
  const FDR_COLOURS = [
    { max: 1e-6, color: '#bf616a', label: '< 1e-6' },
    { max: 1e-4, color: '#d08770', label: '1e-6 – 1e-4' },
    { max: 0.01, color: '#ebcb8b', label: '1e-4 – 0.01' },
    { max: 0.05, color: '#a3be8c', label: '0.01 – 0.05' },
    { max: Infinity, color: '#888888', label: '> 0.05' },
  ];

  // Pattern-label vocabulary (spec §2.3) — atlas does NOT recompute, only colour-codes.
  const PATTERN_LABELS = {
    canonical_breakpoint_marker: { color: '#7d4cdb', display: 'Canonical breakpoint marker' },
    dominant_presence_marker:    { color: '#d08770', display: 'Dominant / presence marker' },
    het_specific_marker:         { color: '#5e81ac', display: 'Het-specific marker' },
    sub_haplotype_marker:        { color: '#88c0d0', display: 'Sub-haplotype marker' },
    internal_linked_marker:      { color: '#ebcb8b', display: 'Internal linked marker' },
    uninformative:               { color: '#888888', display: 'Uninformative' },
  };

  // SV-type glyph + colour (spec §3.4)
  const SV_TYPE_STYLE = {
    BND:   { color: '#bf616a', glyph: '▼', display: 'BND' },
    INV:   { color: '#7d4cdb', glyph: '▼', display: 'INV' },
    DEL:   { color: '#d08770', glyph: '▼', display: 'DEL' },
    DUP:   { color: '#5e81ac', glyph: '▼', display: 'DUP' },
    Other: { color: '#888888', glyph: '▲', display: 'Other' },
  };

  // Karyotype group → display label (atlas-side scrubber locks → spec H-system)
  // spec §5.1: locked_labels carry HOMO_1 / HET / HOMO_2; remap to H1/H1, H1/H2, H2/H2.
  const KARYOTYPE_REMAP = {
    HOMO_1: 'H1/H1',
    HET:    'H1/H2',
    HOMO_2: 'H2/H2',
  };

  // Per-group accent colours (mockup 1, karyotype chips)
  const GROUP_COLOURS = {
    'H1/H1': '#3cc08a', // teal
    'H1/H2': '#f5a524', // amber
    'H2/H2': '#7d4cdb', // violet
  };

  // Default rows-per-page for the main table (spec §3.5)
  const PAGE_SIZES = [10, 25, 50, 100];
  const DEFAULT_PAGE_SIZE = 25;

  // Column registry for the Summary tab. Each column has:
  //   id: stable string id (used for sort state + TSV header)
  //   label: short header shown in the table
  //   group: visual group ('id' | 'loc' | 'gc-h1h1' | 'gc-h1h2' | 'gc-h2h2'
  //          | 'fisher' | 'meta'). Used for the two-row header (top row spans
  //          a group; second row shows individual columns).
  //   get(call): pluck the cell value from a sv_call row
  //   fmt(value, call): format the cell to display string
  //   sortKey(call): comparable value for sorting (numbers preferred for
  //          stable ascending/descending; null falls last)
  //   tsv(call): TSV-export string (no HTML, no formatting tricks)
  //   align: 'left' | 'right' | 'center' for cell alignment
  //
  // The registry is the single source of truth for the table; the renderer
  // iterates it for headers, body cells, sort, and TSV export.
  const TABLE_COLUMNS = (function () {
    const num = v => (v == null || !isFinite(v)) ? null : v;
    // Genotype-count cell accessor. When _state.gtCountsView is present
    // (an UpSet bar is active AND supportLayer is loaded), counts come
    // from the selection-scoped recount instead of the layer's full-cohort
    // counts. This is the step-6 deferred-feature unlock — same column
    // identity, just different denominator.
    const gc = (call, group, allele) => {
      const view = _state && _state.gtCountsView;
      if (view && call && call.sv_id && view[call.sv_id]) {
        const g = view[call.sv_id][group];
        return (g && typeof g[allele] === 'number') ? g[allele] : null;
      }
      const g = call && call.genotype_counts && call.genotype_counts[group];
      return (g && typeof g[allele] === 'number') ? g[allele] : null;
    };
    return [
      // Identity
      { id: 'sv_id', label: 'sv_id', group: 'id', align: 'left',
        get:     c => c.sv_id || '',
        fmt:     v => v || '—',
        sortKey: c => (c.sv_id || ''),
        tsv:     c => c.sv_id || '' },
      { id: 'sv_type', label: 'SV type', group: 'id', align: 'left',
        get:     c => c.sv_type || 'Other',
        fmt:     v => v || '—',
        sortKey: c => (c.sv_type || ''),
        tsv:     c => c.sv_type || '' },
      { id: 'zone', label: 'Zone', group: 'id', align: 'left',
        get:     c => c.zone || '',
        fmt:     v => _zoneDisplay(v),
        sortKey: c => (c.zone || ''),
        tsv:     c => c.zone || '' },
      // Locus
      { id: 'position_bp', label: 'Position (bp)', group: 'loc', align: 'right',
        get:     c => num(c.position_bp),
        fmt:     v => v == null ? '—' : v.toLocaleString(),
        sortKey: c => num(c.position_bp),
        tsv:     c => num(c.position_bp) == null ? '' : String(c.position_bp) },
      { id: 'distance_to_edge_bp', label: 'Dist. to edge (bp)', group: 'loc', align: 'right',
        get:     c => num(c.distance_to_edge_bp),
        fmt:     v => v == null ? '—' : (v >= 0 ? '+' : '') + v.toLocaleString(),
        sortKey: c => num(c.distance_to_edge_bp),
        tsv:     c => num(c.distance_to_edge_bp) == null ? '' : String(c.distance_to_edge_bp) },
      { id: 'n_samples', label: 'Samples (n)', group: 'loc', align: 'right',
        get:     c => num(c.n_samples_with_call),
        fmt:     v => v == null ? '—' : String(v),
        sortKey: c => num(c.n_samples_with_call),
        tsv:     c => num(c.n_samples_with_call) == null ? '' : String(c.n_samples_with_call) },
      // Genotype counts H1/H1
      { id: 'h1h1_AA',   label: 'AA',   group: 'gc-h1h1', align: 'right',
        get: c => gc(c,'H1/H1','AA'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H1','AA'), tsv: c => String(gc(c,'H1/H1','AA') || 0) },
      { id: 'h1h1_AB',   label: 'AB',   group: 'gc-h1h1', align: 'right',
        get: c => gc(c,'H1/H1','AB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H1','AB'), tsv: c => String(gc(c,'H1/H1','AB') || 0) },
      { id: 'h1h1_BB',   label: 'BB',   group: 'gc-h1h1', align: 'right',
        get: c => gc(c,'H1/H1','BB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H1','BB'), tsv: c => String(gc(c,'H1/H1','BB') || 0) },
      { id: 'h1h1_miss', label: 'miss', group: 'gc-h1h1', align: 'right',
        get: c => gc(c,'H1/H1','miss'), fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H1','miss'), tsv: c => String(gc(c,'H1/H1','miss') || 0) },
      // Genotype counts H1/H2
      { id: 'h1h2_AA',   label: 'AA',   group: 'gc-h1h2', align: 'right',
        get: c => gc(c,'H1/H2','AA'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H2','AA'), tsv: c => String(gc(c,'H1/H2','AA') || 0) },
      { id: 'h1h2_AB',   label: 'AB',   group: 'gc-h1h2', align: 'right',
        get: c => gc(c,'H1/H2','AB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H2','AB'), tsv: c => String(gc(c,'H1/H2','AB') || 0) },
      { id: 'h1h2_BB',   label: 'BB',   group: 'gc-h1h2', align: 'right',
        get: c => gc(c,'H1/H2','BB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H2','BB'), tsv: c => String(gc(c,'H1/H2','BB') || 0) },
      { id: 'h1h2_miss', label: 'miss', group: 'gc-h1h2', align: 'right',
        get: c => gc(c,'H1/H2','miss'), fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H1/H2','miss'), tsv: c => String(gc(c,'H1/H2','miss') || 0) },
      // Genotype counts H2/H2
      { id: 'h2h2_AA',   label: 'AA',   group: 'gc-h2h2', align: 'right',
        get: c => gc(c,'H2/H2','AA'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H2/H2','AA'), tsv: c => String(gc(c,'H2/H2','AA') || 0) },
      { id: 'h2h2_AB',   label: 'AB',   group: 'gc-h2h2', align: 'right',
        get: c => gc(c,'H2/H2','AB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H2/H2','AB'), tsv: c => String(gc(c,'H2/H2','AB') || 0) },
      { id: 'h2h2_BB',   label: 'BB',   group: 'gc-h2h2', align: 'right',
        get: c => gc(c,'H2/H2','BB'),   fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H2/H2','BB'), tsv: c => String(gc(c,'H2/H2','BB') || 0) },
      { id: 'h2h2_miss', label: 'miss', group: 'gc-h2h2', align: 'right',
        get: c => gc(c,'H2/H2','miss'), fmt: v => v == null ? '—' : v,
        sortKey: c => gc(c,'H2/H2','miss'), tsv: c => String(gc(c,'H2/H2','miss') || 0) },
      // H1/H1 vs H2/H2: OR / p / FDR
      { id: 'odds_ratio', label: 'OR (Fisher)', group: 'fisher', align: 'right',
        get:     c => num(c.fisher && c.fisher.odds_ratio),
        fmt:     v => v == null ? '—' : (v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2)),
        sortKey: c => num(c.fisher && c.fisher.odds_ratio),
        tsv:     c => { const v = num(c.fisher && c.fisher.odds_ratio); return v == null ? '' : String(v); } },
      { id: 'p_value', label: 'p value', group: 'fisher', align: 'right',
        get:     c => num(c.fisher && c.fisher.p_value),
        fmt:     v => v == null ? '—' : _fmtSci(v),
        sortKey: c => num(c.fisher && c.fisher.p_value),
        tsv:     c => { const v = num(c.fisher && c.fisher.p_value); return v == null ? '' : String(v); } },
      { id: 'fdr', label: 'FDR', group: 'fisher', align: 'right',
        get:     c => num(c.fisher && c.fisher.fdr_bh),
        fmt:     v => v == null ? '—' : _fmtSci(v),
        sortKey: c => num(c.fisher && c.fisher.fdr_bh),
        tsv:     c => { const v = num(c.fisher && c.fisher.fdr_bh); return v == null ? '' : String(v); } },
      // Pattern + notes
      { id: 'pattern', label: 'Pattern (label)', group: 'meta', align: 'left',
        get:     c => c.pattern_label || '',
        fmt:     v => (PATTERN_LABELS[v] && PATTERN_LABELS[v].display) || v || '—',
        sortKey: c => (c.pattern_label || ''),
        tsv:     c => c.pattern_label || '' },
      { id: 'notes', label: 'Notes', group: 'meta', align: 'left',
        get:     c => c.notes || '',
        fmt:     v => v || '—',
        sortKey: c => (c.notes || ''),
        tsv:     c => c.notes || '' },
    ];
  }());

  // Group header spec for the table's two-row header (mockup 1).
  const TABLE_GROUP_HEADERS = [
    { group: 'id',      label: '',                       cols: 3 },
    { group: 'loc',     label: '',                       cols: 3 },
    { group: 'gc-h1h1', label: 'H1/H1 (n=?)',            cols: 4 },
    { group: 'gc-h1h2', label: 'H1/H2 (n=?)',            cols: 4 },
    { group: 'gc-h2h2', label: 'H2/H2 (n=?)',            cols: 4 },
    { group: 'fisher',  label: 'H1/H1 vs H2/H2',         cols: 3 },
    { group: 'meta',    label: '',                       cols: 2 },
  ];

  // Display name for zone identifiers (handles the 5 spec zones + unknowns)
  function _zoneDisplay(z) {
    if (!z) return '—';
    return ({
      left_flank:     'Left flank',
      left_boundary:  'Left boundary',
      inversion_body: 'Inversion body',
      right_boundary: 'Right boundary',
      right_flank:    'Right flank',
    })[z] || z;
  }

  // Compact scientific formatter: 2.1e-9 → "2.1e-9"; 0.0041 → "0.0041";
  // 0.43 → "0.43". Avoids "4.1000000000000004e-3" output.
  function _fmtSci(v) {
    if (v == null || !isFinite(v)) return '—';
    if (v === 0) return '0';
    const abs = Math.abs(v);
    if (abs >= 0.001 && abs < 1) return v.toFixed(4).replace(/0+$/, '').replace(/\.$/, '');
    if (abs >= 1) return v.toFixed(2);
    // scientific
    const exp = Math.floor(Math.log10(abs));
    const mant = v / Math.pow(10, exp);
    const mantStr = (Math.abs(mant) >= 10) ? mant.toFixed(0)
                                            : mant.toFixed(1).replace(/\.0$/, '');
    return mantStr + 'e' + exp;
  }


  // ===========================================================================
  // Internal state — module-scoped, read by .refresh() / debug
  // ===========================================================================

  const _state = {
    mounted: false,
    activeCandidateId: null,
    layer: null,           // sv_genotype_counts_v1 JSON for active candidate
    layerLoading: false,
    layerError: null,
    filters: { ...DEFAULT_FILTERS },
    rootEl: null,
    // Table state — sort + pagination. Persisted in memory only (filters
    // persist to localStorage; table-view state resets on candidate switch
    // because column visibility doesn't change between candidates).
    sortColumnId: 'fdr',   // default sort: FDR ascending (spec §3.5)
    sortDirection: 'asc',  // 'asc' | 'desc'
    pageSize: DEFAULT_PAGE_SIZE,
    pageIndex: 0,          // 0-based
    highlightedSvId: null, // for row ↔ locus glyph crosstalk (step 3)
    // Cursor + sticky markers (step 3.5). cursor_bp follows the mouse on
    // hover and snaps to track ticks on click/keyboard. markerLeft/Right
    // are sticky pins that the user drops with E/F; they're VISUAL ONLY
    // (this page is read-only per spec §1 — boundary editing belongs to
    // page 11). View preset = which window of bp to show: 'default' |
    // 'left_close' | 'right_close' (spec §3.2).
    cursorBp: null,        // number (bp) | null
    markerLeftBp: null,    // bp pin from "E" key
    markerRightBp: null,   // bp pin from "F" key
    viewPreset: 'default', // current preset id
    hotkeysAttached: false,
    // Step 4 additions:
    highlightPattern: null, // pattern_label key (clicked legend) | null
    rowAnnotations: {},     // { sv_id: 'L' | 'R' | 'star' | null } per active candidate
    // Step 4.5 — region-select mode:
    selectMode: 'zoom',          // 'zoom' (default — click=zoom-via-cursor) | 'select' (drag=region)
    selection: null,             // { startBp, endBp } | null — sticky region filter
    _selectDrag: null,           // { startBp, currentBp } during an active drag
    // Step 5 — UpSet panel (the redirect):
    //   combinationsLayer = parsed sv_evidence_combinations_v1.json (sibling
    //     of sv_genotype_counts in the per-candidate folder).
    //   selectedSamples = a Set<sample_id> populated when the user clicks a
    //     bar in the UpSet panel; locus dims non-selected SV glyphs and the
    //     main table re-counts genotype columns within selected_samples.
    //   activeCombinationIndex = the index of the clicked bar (for visual
    //     "active bar" styling). null = no active selection.
    combinationsLayer: null,
    selectedSamples: null,       // Set<sample_id> | null
    activeCombinationIndex: null,
    // Step 6 — sample × SV dosage heatmap:
    //   supportLayer = parsed sv_support_by_sample_v1.json (third sibling
    //     in the per-candidate folder). Compact dosage matrix:
    //     rows = samples, cols = SVs, '0'=AA, '1'=AB, '2'=BB, '.'=miss.
    //   heatmapView = 'compact' (right-rail thumbnail) | 'large' (main-area
    //     overlay panel). Toggle via the "expand" button next to the heatmap.
    supportLayer: null,
    heatmapView: 'compact',
    heatmapHoverCell: null,      // { rowIdx, colIdx } during cell hover
    // Selection-scoped recount (step 6 deferred-feature unlock):
    //   When `selectedSamples` is set AND `supportLayer` is loaded, the
    //   table's H1/H1 / H1/H2 / H2/H2 columns recount within the selected
    //   fish only. The recount is cached here so we don't traverse the
    //   support matrix on every cell access. Set by _recomputeGtCountsView,
    //   cleared whenever selection or support changes. Shape:
    //     { sv_id: { 'H1/H1': {AA, AB, BB, miss}, 'H1/H2': {...}, 'H2/H2': {...} },
    //       _groupNs: { 'H1/H1': nSel, 'H1/H2': nSel, 'H2/H2': nSel } }
    //   `_groupNs` is the per-group sample count *within the selection*,
    //   used for the table's group-header `(n=…)` annotation.
    gtCountsView: null,
  };

  // ===========================================================================
  // Filter / sort / paginate engine
  // ===========================================================================

  // Apply the user's filter selections (sv_type / quality / zone / min_samples /
  // show_only_associated / fdr_threshold) to the layer's sv_calls.
  // Returns a new array; does not mutate the layer.
  function _applyFilters(calls, filters) {
    if (!Array.isArray(calls)) return [];
    const f = filters || DEFAULT_FILTERS;
    // Step 4.5 — sticky region selection. Read from module _state, not
    // filters, because selection is a transient locus-level filter (set
    // by drag, cleared by Esc), not a persisted filter form value.
    const sel = _state.selection;
    return calls.filter(c => {
      if (f.sv_type && f.sv_type !== 'All' && c.sv_type !== f.sv_type) return false;
      if (f.quality === 'PASS' && c.quality && c.quality !== 'PASS') return false;
      // zone
      if (f.zone === 'boundary_500kb') {
        if (c.zone !== 'left_boundary' && c.zone !== 'right_boundary') return false;
      } else if (f.zone === 'inversion_body') {
        if (c.zone !== 'inversion_body') return false;
      } else if (f.zone === 'flanks') {
        if (c.zone !== 'left_flank' && c.zone !== 'right_flank') return false;
      }
      // min_samples
      if (f.min_samples != null && (c.n_samples_with_call || 0) < f.min_samples) return false;
      // show_only_associated → require fdr_bh < threshold
      if (f.show_only_associated) {
        const fdr = c.fisher && c.fisher.fdr_bh;
        if (fdr == null || !isFinite(fdr) || fdr >= (f.fdr_threshold || 0.05)) return false;
      }
      // Region selection (step 4.5): keep only SVs whose position falls
      // inside [startBp, endBp]. Both endpoints inclusive.
      if (sel && sel.startBp != null && sel.endBp != null) {
        const pos = c.position_bp;
        if (pos == null || pos < sel.startBp || pos > sel.endBp) return false;
      }
      return true;
    });
  }

  // Sort by the active sort column (TABLE_COLUMNS entry by id). Stable
  // tiebreak on sv_id to prevent flicker when multiple rows share a key.
  function _sortRows(rows, columnId, direction) {
    const col = TABLE_COLUMNS.find(c => c.id === columnId) ||
                TABLE_COLUMNS.find(c => c.id === 'fdr');
    const dir = direction === 'desc' ? -1 : 1;
    const out = rows.slice();
    out.sort((a, b) => {
      const ka = col.sortKey(a);
      const kb = col.sortKey(b);
      // null/undefined falls last regardless of direction
      const aNull = (ka == null || (typeof ka === 'number' && !isFinite(ka)));
      const bNull = (kb == null || (typeof kb === 'number' && !isFinite(kb)));
      if (aNull && bNull) return (a.sv_id || '').localeCompare(b.sv_id || '');
      if (aNull) return 1;
      if (bNull) return -1;
      if (ka < kb) return -1 * dir;
      if (ka > kb) return  1 * dir;
      return (a.sv_id || '').localeCompare(b.sv_id || '');
    });
    return out;
  }

  // Slice a sorted+filtered row set into a single page.
  function _paginateRows(rows, pageSize, pageIndex) {
    const total = rows.length;
    const lastPage = Math.max(0, Math.ceil(total / pageSize) - 1);
    const idx = Math.min(Math.max(0, pageIndex), lastPage);
    const start = idx * pageSize;
    return {
      rows: rows.slice(start, start + pageSize),
      total: total,
      pageIndex: idx,
      lastPage: lastPage,
      pageSize: pageSize,
      startRank: total === 0 ? 0 : start + 1,
      endRank:   Math.min(total, start + pageSize),
    };
  }

  // Compose: filter → sort → paginate. Pure helper, no DOM. Used by the
  // renderer + by exportFilteredTSV.
  function _computeVisibleRows(layer, filters, sortColumnId, sortDir, pageSize, pageIndex) {
    if (!layer || !Array.isArray(layer.sv_calls)) {
      return { rows: [], total: 0, pageIndex: 0, lastPage: 0, pageSize, startRank: 0, endRank: 0,
               filteredAll: [] };
    }
    const filtered = _applyFilters(layer.sv_calls, filters);
    const sorted   = _sortRows(filtered, sortColumnId, sortDir);
    const paged    = _paginateRows(sorted, pageSize, pageIndex);
    paged.filteredAll = sorted; // full filtered+sorted set (for TSV export)
    return paged;
  }

  // ===========================================================================
  // Utilities
  // ===========================================================================

  function _esc(s) {
    if (s == null) return '';
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
                    .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }

  function _fmtBp(bp) {
    if (bp == null || !isFinite(bp)) return '—';
    if (Math.abs(bp) >= 1e6) return (bp / 1e6).toFixed(2) + ' Mb';
    if (Math.abs(bp) >= 1e3) return (bp / 1e3).toFixed(1) + ' kb';
    return bp + ' bp';
  }

  function _fmtMb(bp) {
    if (bp == null || !isFinite(bp)) return '—';
    return (bp / 1e6).toFixed(2) + ' Mb';
  }

  function _fdrColour(fdr) {
    if (fdr == null || !isFinite(fdr)) return '#888888';
    for (const tier of FDR_COLOURS) {
      if (fdr < tier.max) return tier.color;
    }
    return '#888888';
  }

  // Resolve the active candidate id from the global atlas state. Multiple
  // call sites use different keys (state.candidate, state.activeCandidateId,
  // state.candidateState[cid] — depends on which page set it). Probe all.
  function _resolveActiveCandidateId() {
    const state = (typeof window !== 'undefined') ? window.state : null;
    if (!state) return null;
    if (state.activeCandidateId) return state.activeCandidateId;
    if (state.candidate && state.candidate.id) return state.candidate.id;
    if (state.candidate && state.candidate.candidate_id) return state.candidate.candidate_id;
    return null;
  }

  // Resolve karyotype-group membership for the active candidate.
  // Returns {groupOfSample: Map<sampleId, "H1/H1"|"H1/H2"|"H2/H2">,
  //          counts: {H1/H1: n, H1/H2: n, H2/H2: n}, locked: bool}.
  function _resolveKaryotypeGroups(cid) {
    const state = (typeof window !== 'undefined') ? window.state : null;
    if (!state || !state.candidateState || !cid) {
      return { groupOfSample: new Map(), counts: {}, locked: false };
    }
    const cs = state.candidateState[cid] || {};
    const labels = cs.locked_labels || [];
    if (!Array.isArray(labels) || labels.length === 0) {
      return { groupOfSample: new Map(), counts: {}, locked: false };
    }
    const groupOfSample = new Map();
    const counts = { 'H1/H1': 0, 'H1/H2': 0, 'H2/H2': 0 };
    for (const r of labels) {
      const sid = r.sample_id || r.sid;
      const grp = KARYOTYPE_REMAP[r.label] || null;
      if (sid && grp) {
        groupOfSample.set(sid, grp);
        counts[grp] = (counts[grp] || 0) + 1;
      }
    }
    return { groupOfSample, counts, locked: true };
  }

  // Persist filter state per-candidate (spec §3.3, mirrors page-11 pattern).
  function _readFilters(cid) {
    if (!cid) return { ...DEFAULT_FILTERS };
    try {
      const raw = window.localStorage.getItem(LS_FILTERS);
      if (!raw) return { ...DEFAULT_FILTERS };
      const all = JSON.parse(raw);
      return { ...DEFAULT_FILTERS, ...(all[cid] || {}) };
    } catch (_) {
      return { ...DEFAULT_FILTERS };
    }
  }

  function _writeFilters(cid, filters) {
    if (!cid) return;
    try {
      const raw = window.localStorage.getItem(LS_FILTERS);
      const all = raw ? JSON.parse(raw) : {};
      all[cid] = filters;
      window.localStorage.setItem(LS_FILTERS, JSON.stringify(all));
    } catch (_) { /* localStorage full / disabled — silently skip */ }
  }

  // Read/write row annotations (step 4). Annotations are { sv_id: 'L'|'R'|'star' }
  // per candidate. They're a researcher's note-to-self ("this SV looks like
  // it's the actual left breakpoint") — not a structural change to the
  // candidate's boundaries (that's page 11). Persisted per-candidate in
  // localStorage so they survive page refreshes.
  function _readAnnotations(cid) {
    if (!cid) return {};
    try {
      const raw = window.localStorage.getItem(LS_ANNOTATIONS);
      if (!raw) return {};
      const all = JSON.parse(raw);
      return all[cid] || {};
    } catch (_) { return {}; }
  }

  function _writeAnnotations(cid, annos) {
    if (!cid) return;
    try {
      const raw = window.localStorage.getItem(LS_ANNOTATIONS);
      const all = raw ? JSON.parse(raw) : {};
      all[cid] = annos;
      window.localStorage.setItem(LS_ANNOTATIONS, JSON.stringify(all));
    } catch (_) { /* localStorage full / disabled — silently skip */ }
  }

  // Cycle the annotation for one SV: null → L → R → star → null.
  // Mutates _state.rowAnnotations and persists.
  function _cycleAnnotation(svId) {
    const cur = _state.rowAnnotations[svId] || null;
    const next = ({ null: 'L', L: 'R', R: 'star', star: null })[String(cur)];
    if (next == null) delete _state.rowAnnotations[svId];
    else              _state.rowAnnotations[svId] = next;
    _writeAnnotations(_state.activeCandidateId, _state.rowAnnotations);
  }

  // Compact visual chip for an annotation. Used in the table sv_id cell.
  // Pre-defined glyphs: L (green) = putative left breakpoint, R (orange) =
  // putative right, ★ (amber) = starred / interesting. Empty = no anno.
  function _renderAnnoChip(anno) {
    if (anno === 'L')    return '<span class="sv-anno sv-anno-L">L</span>';
    if (anno === 'R')    return '<span class="sv-anno sv-anno-R">R</span>';
    if (anno === 'star') return '<span class="sv-anno sv-anno-star">★</span>';
    return '<span class="sv-anno sv-anno-empty">·</span>';
  }

  // ===========================================================================
  // Layer loader — fetch json/sv_genotype_counts/<cid>.json
  // ===========================================================================

  function _loadLayer(cid) {
    _state.layerLoading = true;
    _state.layerError = null;
    _state.layer = null;
    const url = 'json/sv_genotype_counts/' + encodeURIComponent(cid) + '.json';
    return fetch(url)
      .then(resp => {
        if (!resp.ok) throw new Error('HTTP ' + resp.status);
        return resp.json();
      })
      .then(json => {
        // Light validation — reject if format_version doesn't match.
        if (json && json.format_version && json.format_version !== 'sv_genotype_counts_v1') {
          throw new Error('unexpected format_version: ' + json.format_version);
        }
        _state.layer = json;
        _state.layerLoading = false;
        return json;
      })
      .catch(err => {
        _state.layerError = err && err.message ? err.message : String(err);
        _state.layerLoading = false;
        _state.layer = null;
        return null;
      });
  }

  // ===========================================================================
  // Rendering — empty shell + dispatcher
  // ===========================================================================

  // Render the page-level shell (left rail / main / right rail boxes) into
  // _state.rootEl. Idempotent: builds DOM once, then subsequent renders only
  // patch in-place. For step 1 we ship the shell + empty-state messages only.
  function _renderShell() {
    if (!_state.rootEl) return;
    if (_state.rootEl.dataset.svShellMounted === '1') return;
    _state.rootEl.dataset.svShellMounted = '1';
    _state.rootEl.innerHTML = `
      <div class="sv-evidence-page" data-sv-page>
        <div class="sv-toolbar" data-sv-toolbar>
          <div class="sv-toolbar-locus" data-sv-locus-label>Locus: —</div>
          <div class="sv-toolbar-spacer"></div>
          <div class="sv-toolbar-presets">
            <label>View presets:</label>
            <select data-sv-preset>
              <option value="default">Default</option>
              <option value="left_close">Left boundary close-up</option>
              <option value="right_close">Right boundary close-up</option>
            </select>
          </div>
          <div class="sv-toolbar-mode">
            <label>Mode:</label>
            <button data-sv-mode="zoom" class="sv-mode-btn"
                    title="Zoom mode — click locus to set cursor + nearest SV">🔍 zoom</button>
            <button data-sv-mode="select" class="sv-mode-btn"
                    title="Select mode — drag on locus to filter table + locus to a region">◧ select</button>
          </div>
          <div class="sv-toolbar-zoom">
            <button data-sv-zoom="in"  title="Zoom in (tighten ±10% around cursor)">＋</button>
            <button data-sv-zoom="out" title="Zoom out (widen ±10%)">－</button>
            <button data-sv-zoom="fit" title="Fit (back to Default preset)">⛶</button>
            <button data-sv-zoom="reset" title="Reset cursor + markers">↺</button>
          </div>
        </div>
        <div class="sv-grid">
          <aside class="sv-leftrail" data-sv-leftrail>
            <div class="sv-section">
              <div class="sv-section-title">Candidate</div>
              <div data-sv-candidate-context class="sv-candidate-context">
                <div class="sv-empty-hint">No candidate active.</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Karyotype</div>
              <div data-sv-karyotype class="sv-karyotype-block">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Boundaries</div>
              <div data-sv-boundaries class="sv-boundaries-block">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Zones (±500 kb)</div>
              <div data-sv-zones class="sv-zones-block">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Filter SVs</div>
              <div data-sv-filters class="sv-filters-block">
                <div class="sv-empty-hint">Load a candidate to enable filters.</div>
              </div>
            </div>
          </aside>
          <main class="sv-main" data-sv-main>
            <div class="sv-locus-strip" data-sv-locus-strip>
              <div class="sv-empty-state" data-sv-empty-main>
                <div class="sv-empty-title">SV Evidence</div>
                <div class="sv-empty-body">
                  Promote a candidate (page 3 catalogue or page 1 local PCA |z|),
                  lock karyotype groups on the scrubber, then return to this tab.
                </div>
              </div>
            </div>
            <div class="sv-table-wrap" data-sv-table-wrap></div>
          </main>
          <aside class="sv-rightrail" data-sv-rightrail>
            <div class="sv-section">
              <div class="sv-section-title">Boundary summary</div>
              <div data-sv-boundary-summary class="sv-boundary-summary">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Legend</div>
              <div data-sv-legend class="sv-legend">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">UpSet (top 12 associated SVs)</div>
              <div data-sv-upset class="sv-upset">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
            <div class="sv-section">
              <div class="sv-section-title">Sample × SV heatmap</div>
              <div data-sv-heatmap class="sv-heatmap">
                <div class="sv-empty-hint">—</div>
              </div>
            </div>
          </aside>
        </div>
      </div>
    `;
    _injectStylesOnce();
    // Wire drag-drop on the outer page wrapper. Drop a JSON anywhere on
    // the page → ingest as the active layer (matches the rest of the
    // atlas's "drag-drop a JSON" convention).
    const wrap = _state.rootEl.querySelector('.sv-evidence-page');
    if (wrap) {
      wrap.addEventListener('dragover', (ev) => {
        ev.preventDefault();
        ev.stopPropagation();
        if (ev.dataTransfer) ev.dataTransfer.dropEffect = 'copy';
        wrap.classList.add('sv-drop-active');
      });
      wrap.addEventListener('dragleave', (ev) => {
        // Only clear if leaving the wrap (not a child). Cheap check via
        // relatedTarget containment.
        if (!wrap.contains(ev.relatedTarget)) {
          wrap.classList.remove('sv-drop-active');
        }
      });
      wrap.addEventListener('drop', _handlePageDrop);
    }

    // Wire toolbar: view-preset dropdown + zoom buttons + mode toggle.
    const presetEl = _state.rootEl.querySelector('[data-sv-preset]');
    if (presetEl) {
      presetEl.value = _state.viewPreset;
      presetEl.addEventListener('change', () => _setViewPreset(presetEl.value));
    }
    _state.rootEl.querySelectorAll('[data-sv-zoom]').forEach(btn => {
      btn.addEventListener('click', () => _onZoomButton(btn.dataset.svZoom));
    });
    _state.rootEl.querySelectorAll('[data-sv-mode]').forEach(btn => {
      btn.addEventListener('click', () => _setSelectMode(btn.dataset.svMode));
    });
    _refreshModeButtonStyles();
  }

  // Reflect the active select-mode in the toolbar buttons (highlights the
  // active one). Called after _renderShell wiring and after _setSelectMode.
  function _refreshModeButtonStyles() {
    if (!_state.rootEl) return;
    _state.rootEl.querySelectorAll('[data-sv-mode]').forEach(btn => {
      if (btn.dataset.svMode === _state.selectMode) btn.classList.add('sv-mode-active');
      else                                          btn.classList.remove('sv-mode-active');
    });
  }

  // Switch select mode. Clears any in-progress drag but preserves the
  // current sticky selection (so the user can switch back to zoom mode
  // without losing their region — Esc clears the region explicitly).
  function _setSelectMode(mode) {
    if (mode !== 'zoom' && mode !== 'select') return;
    _state.selectMode = mode;
    _state._selectDrag = null;
    _refreshModeButtonStyles();
    // Update cursor on the hit-rect so the mode is visually obvious.
    if (_state.rootEl) {
      const hit = _state.rootEl.querySelector('[data-sv-locus-hit]');
      if (hit) hit.style.cursor = (mode === 'select') ? 'col-resize' : 'crosshair';
    }
  }

  // Toolbar zoom buttons (step 3.5). +/- tighten/widen the window by 10%
  // around the cursor (or window centre if no cursor); fit returns to the
  // Default preset; reset clears cursor + markers + filters.
  function _onZoomButton(action) {
    const layer = _state.layer;
    const w = _state._lastWindow;
    if (!layer || !w) return;
    if (action === 'fit') { _setViewPreset('default'); return; }
    if (action === 'reset') {
      _state.cursorBp = null;
      _state.markerLeftBp = null;
      _state.markerRightBp = null;
      _state.highlightedSvId = null;
      _state.selection = null;       // step 4.5: clear region selection too
      _state._selectDrag = null;
      _state.filters = { ...DEFAULT_FILTERS };
      _writeFilters(_state.activeCandidateId, _state.filters);
      _state.viewPreset = 'default';
      const presetEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-preset]');
      if (presetEl) presetEl.value = 'default';
      const filEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-filters]');
      if (filEl) _renderFiltersBlock(filEl);
      _renderLocusStrip();
      _renderTable();
      return;
    }
    // in / out — local tighten/widen around cursor
    const span = w.hi - w.lo;
    const factor = (action === 'in') ? 0.8 : 1.25;
    const newSpan = span * factor;
    const centre = (_state.cursorBp != null) ? _state.cursorBp : (w.lo + span / 2);
    const newLo = centre - newSpan / 2;
    const newHi = centre + newSpan / 2;
    // Stash a custom window override on _state so _renderLocusStrip uses it
    _state.viewPreset = '_custom';
    _state._customWindow = { lo: newLo, hi: newHi };
    _renderLocusStrip();
  }

  // Update left-rail candidate context + karyotype block + boundaries +
  // zones using whatever combination of (candidate, locked_labels, layer)
  // is currently available. Always renders something — never throws.
  function _renderLeftRail() {
    if (!_state.rootEl) return;
    const cid = _state.activeCandidateId;
    const layer = _state.layer;
    const groups = _resolveKaryotypeGroups(cid);

    // --- Candidate context -------------------------------------------------
    const ctxEl = _state.rootEl.querySelector('[data-sv-candidate-context]');
    if (ctxEl) {
      if (!cid) {
        ctxEl.innerHTML = '<div class="sv-empty-hint">No candidate active.</div>';
      } else {
        const chrom = layer ? layer.chrom : (_resolveCandidateChrom(cid) || '—');
        const left  = layer ? layer.boundary_left_bp  : null;
        const right = layer ? layer.boundary_right_bp : null;
        const span  = (left != null && right != null) ? (right - left) : null;
        ctxEl.innerHTML = `
          <div class="sv-row"><span class="sv-key">Candidate</span>
                              <span class="sv-val sv-strong">${_esc(cid)}</span></div>
          <div class="sv-row"><span class="sv-key">Chromosome</span>
                              <span class="sv-val">${_esc(chrom)}</span></div>
          <div class="sv-row"><span class="sv-key">Region</span>
                              <span class="sv-val">${left != null && right != null
                                  ? _fmtMb(left) + ' – ' + _fmtMb(right) : '—'}</span></div>
          <div class="sv-row"><span class="sv-key">Span</span>
                              <span class="sv-val">${span != null ? _fmtMb(span) : '—'}</span></div>
        `;
      }
    }

    // --- Karyotype block ---------------------------------------------------
    const karyEl = _state.rootEl.querySelector('[data-sv-karyotype]');
    if (karyEl) {
      if (!groups.locked) {
        karyEl.innerHTML = '<div class="sv-empty-hint">' +
          'Lock karyotype groups on the scrubber (page 1 / page 4) first.' +
          '</div>';
      } else {
        const chips = ['H1/H1', 'H1/H2', 'H2/H2'].map(g => {
          const n = groups.counts[g] || 0;
          const c = GROUP_COLOURS[g];
          return `
            <button class="sv-karyo-chip" data-sv-group="${g}"
                    style="border-color:${c};color:${c};">
              <span class="sv-karyo-name">${g}</span>
              <span class="sv-karyo-n">n=${n}</span>
            </button>`;
        }).join('');
        karyEl.innerHTML = '<div class="sv-karyo-row">' + chips + '</div>';
      }
    }

    // --- Boundaries --------------------------------------------------------
    const bndEl = _state.rootEl.querySelector('[data-sv-boundaries]');
    if (bndEl) {
      if (!layer) {
        bndEl.innerHTML = '<div class="sv-empty-hint">' +
          (cid ? 'sv_genotype_counts layer not loaded' : '—') +
          '</div>';
      } else {
        bndEl.innerHTML = `
          <div class="sv-row"><span class="sv-key">Left</span>
                              <span class="sv-val">${_esc(layer.boundary_left_bp.toLocaleString())} bp</span></div>
          <div class="sv-row"><span class="sv-key">Right</span>
                              <span class="sv-val">${_esc(layer.boundary_right_bp.toLocaleString())} bp</span></div>
        `;
      }
    }

    // --- Zones (±500 kb) ---------------------------------------------------
    const zoneEl = _state.rootEl.querySelector('[data-sv-zones]');
    if (zoneEl) {
      if (!layer || !layer.zone_definitions_bp) {
        zoneEl.innerHTML = '<div class="sv-empty-hint">—</div>';
      } else {
        const z = layer.zone_definitions_bp;
        const zoneRows = [
          ['Left boundary',  z.left_boundary,  '#4fa3ff'],
          ['Inversion body', z.inversion_body, '#f5a524'],
          ['Right boundary', z.right_boundary, '#bf616a'],
          ['Left flank',     z.left_flank,     '#5a6472'],
          ['Right flank',    z.right_flank,    '#5a6472'],
        ];
        zoneEl.innerHTML = zoneRows
          .filter(([, range]) => Array.isArray(range))
          .map(([label, range, c]) => `
            <div class="sv-zone-row">
              <span class="sv-zone-swatch" style="background:${c}"></span>
              <span class="sv-zone-label">${_esc(label)}</span>
              <span class="sv-zone-range">${_fmtMb(range[0])} – ${_fmtMb(range[1])}</span>
            </div>
          `).join('');
      }
    }

    // --- Filters block -----------------------------------------------------
    // The filter block is enabled regardless of whether a layer is loaded —
    // changing filters before the layer arrives is harmless (they apply
    // when the layer lands). Only the "Load layer" button is unconditional.
    const filEl = _state.rootEl.querySelector('[data-sv-filters]');
    if (filEl) {
      _renderFiltersBlock(filEl);
    }
  }

  // Render the left-rail filter block (spec §3.3). The form is single-pane
  // (no "Apply"/"Reset" delays — every change is live, mirroring how the
  // rest of the atlas behaves) but a Reset button is still useful for
  // returning to defaults in one click.
  function _renderFiltersBlock(filEl) {
    if (!filEl) return;
    const F = _state.filters;
    filEl.innerHTML = `
      <div class="sv-filter-row">
        <label>SV type</label>
        <select data-sv-filter="sv_type">
          ${['All','BND','INV','DEL','DUP','Other'].map(o =>
            `<option value="${o}"${F.sv_type === o ? ' selected' : ''}>${o}</option>`).join('')}
        </select>
      </div>
      <div class="sv-filter-row">
        <label>Quality</label>
        <select data-sv-filter="quality">
          <option value="PASS"${F.quality === 'PASS' ? ' selected' : ''}>PASS</option>
          <option value="All"${F.quality === 'All'  ? ' selected' : ''}>All</option>
        </select>
      </div>
      <div class="sv-filter-row">
        <label>Zone</label>
        <select data-sv-filter="zone">
          <option value="all"             ${F.zone === 'all'             ? 'selected' : ''}>All</option>
          <option value="boundary_500kb"  ${F.zone === 'boundary_500kb'  ? 'selected' : ''}>Boundary ±500 kb</option>
          <option value="inversion_body"  ${F.zone === 'inversion_body'  ? 'selected' : ''}>Inversion body</option>
          <option value="flanks"          ${F.zone === 'flanks'          ? 'selected' : ''}>Flanks only</option>
        </select>
      </div>
      <div class="sv-filter-row">
        <label>Min samples</label>
        <input type="range" min="1" max="226" step="1"
               value="${F.min_samples}" data-sv-filter="min_samples">
        <span class="sv-filter-num" data-sv-filter-display="min_samples">${F.min_samples}</span>
      </div>
      <div class="sv-filter-row sv-filter-checkrow">
        <input type="checkbox" id="sv-fil-only-assoc"
               data-sv-filter="show_only_associated"
               ${F.show_only_associated ? 'checked' : ''}>
        <label for="sv-fil-only-assoc">Show only associated</label>
      </div>
      <div class="sv-filter-row">
        <label>FDR threshold</label>
        <select data-sv-filter="fdr_threshold">
          ${[1e-6, 1e-4, 0.01, 0.05, 0.10].map(v => {
            const sel = (Math.abs(F.fdr_threshold - v) < 1e-15) ? ' selected' : '';
            return `<option value="${v}"${sel}>${_fmtSci(v)}</option>`;
          }).join('')}
        </select>
      </div>
      <div class="sv-filter-actions">
        <button data-sv-action="reset-filters" class="sv-btn-secondary">Reset</button>
        <button data-sv-action="load-layer" class="sv-btn-primary">Load layer…</button>
      </div>
      <div class="sv-filter-hint">
        Layers ship as <code>sv_genotype_counts/&lt;cid&gt;.json</code>.
        Drop a file on this page or click <em>Load layer…</em>.
      </div>
    `;
    _wireFilterInputs(filEl);
  }

  // Wire the filter form's change events. Live-updates _state.filters and
  // re-renders the table without a debounce — the filter pipeline is cheap
  // (single linear pass over a bounded set; each candidate caps at a few
  // hundred SV calls).
  function _wireFilterInputs(filEl) {
    if (!filEl) return;
    const onFilterChange = (key, raw) => {
      let v = raw;
      if (key === 'min_samples')           v = parseInt(raw, 10) || 1;
      else if (key === 'fdr_threshold')    v = parseFloat(raw);
      else if (key === 'show_only_associated') v = !!raw;
      _state.filters = { ..._state.filters, [key]: v };
      _writeFilters(_state.activeCandidateId, _state.filters);
      _state.pageIndex = 0;
      // Live display update for the slider numeric.
      const disp = filEl.querySelector('[data-sv-filter-display="' + key + '"]');
      if (disp) disp.textContent = String(v);
      _renderTable();
      _renderLocusStrip();
    };
    filEl.querySelectorAll('select[data-sv-filter], input[data-sv-filter]').forEach(el => {
      const key = el.dataset.svFilter;
      const evt = (el.tagName === 'INPUT' && el.type === 'range') ? 'input' :
                  (el.tagName === 'INPUT' && el.type === 'checkbox') ? 'change' : 'change';
      el.addEventListener(evt, () => {
        const val = (el.type === 'checkbox') ? el.checked : el.value;
        onFilterChange(key, val);
      });
    });
    const reset = filEl.querySelector('[data-sv-action="reset-filters"]');
    if (reset) reset.addEventListener('click', () => {
      _state.filters = { ...DEFAULT_FILTERS };
      _writeFilters(_state.activeCandidateId, _state.filters);
      _state.pageIndex = 0;
      _renderFiltersBlock(filEl);
      _renderTable();
      _renderLocusStrip();
    });
    const load = filEl.querySelector('[data-sv-action="load-layer"]');
    if (load) load.addEventListener('click', _triggerLoadLayerFile);
  }

  // Update the toolbar locus label and main empty-state message based on
  // current loading / candidate / layer state.
  function _renderMain() {
    if (!_state.rootEl) return;
    const cid = _state.activeCandidateId;
    const layer = _state.layer;
    const groups = _resolveKaryotypeGroups(cid);

    // Locus label
    const locusEl = _state.rootEl.querySelector('[data-sv-locus-label]');
    if (locusEl) {
      if (!layer) {
        locusEl.textContent = 'Locus: —';
      } else {
        const z = layer.zone_definitions_bp || {};
        const lo = (z.left_flank  || [])[0];
        const hi = (z.right_flank || [])[1];
        const span = (lo != null && hi != null) ? (hi - lo) : null;
        locusEl.textContent = 'Locus: ' + (layer.chrom || '—') +
          ((lo != null && hi != null)
            ? '  ' + _fmtMb(lo) + ' – ' + _fmtMb(hi) + '  (' +
              (span != null ? (span / 1e6).toFixed(2) + ' Mb window' : '—') + ')'
            : '');
      }
    }

    // Main empty / partial message
    const emptyEl = _state.rootEl.querySelector('[data-sv-empty-main]');
    if (!emptyEl) return;

    // Reset the locus-strip layout mode each render. _renderLocusStrip() will
    // re-add the .sv-locus-has-content class when actually rendering tracks.
    const stripEl = _state.rootEl.querySelector('[data-sv-locus-strip]');
    if (stripEl) stripEl.classList.remove('sv-locus-has-content');

    if (_state.layerLoading) {
      emptyEl.innerHTML = `
        <div class="sv-empty-title">Loading…</div>
        <div class="sv-empty-body">Fetching <code>sv_genotype_counts/${_esc(cid)}.json</code>.</div>
      `;
      emptyEl.style.display = '';
    } else if (!cid) {
      emptyEl.innerHTML = `
        <div class="sv-empty-title">SV Evidence</div>
        <div class="sv-empty-body">
          Promote a candidate (page 3 catalogue or page 1 local PCA |z|),
          lock karyotype groups on the scrubber, then return to this tab.
        </div>
      `;
      emptyEl.style.display = '';
    } else if (_state.layerError) {
      emptyEl.innerHTML = `
        <div class="sv-empty-title">No layer for <code>${_esc(cid)}</code></div>
        <div class="sv-empty-body">
          Load <code>json/sv_genotype_counts/${_esc(cid)}.json</code> to populate this page.<br>
          <span class="sv-empty-detail">(${_esc(_state.layerError)})</span>
        </div>
      `;
      emptyEl.style.display = '';
    } else if (!groups.locked) {
      emptyEl.innerHTML = `
        <div class="sv-empty-title">Karyotype groups not locked</div>
        <div class="sv-empty-body">
          Lock H1/H1, H1/H2, H2/H2 on the scrubber (page 1 or page 4) for this
          candidate, then refresh.
        </div>
      `;
      emptyEl.style.display = '';
    } else if (layer) {
      // Layer + groups present — locus strip renders into [data-sv-locus-strip],
      // table renders into [data-sv-table-wrap]. Hide the empty-state card.
      emptyEl.style.display = 'none';
      _renderLocusStrip();
      _renderTable();
    } else {
      emptyEl.style.display = 'none';
    }
  }

  function _renderRightRail() {
    if (!_state.rootEl) return;
    const layer = _state.layer;

    // --- Boundary summary (spec §3.6) -------------------------------------
    // Two stacked tables: left boundary, right boundary. Each shows
    // SV type | # SVs | # Associated (FDR < 0.05). Click any row to
    // filter the main table + locus to that SV type + that boundary side.
    const bsEl = _state.rootEl.querySelector('[data-sv-boundary-summary]');
    if (bsEl) {
      if (!layer) {
        bsEl.innerHTML = '<div class="sv-empty-hint">—</div>';
      } else {
        bsEl.innerHTML = _renderBoundarySummaryHTML(layer);
        _wireBoundarySummary(bsEl);
      }
    }

    // --- Legend (spec §3.6) -----------------------------------------------
    // Two side-by-side blocks: FDR colour tiers, Pattern-label colours.
    // Each entry is interactive — clicking sets the corresponding filter.
    const legEl = _state.rootEl.querySelector('[data-sv-legend]');
    if (legEl) {
      if (!layer) {
        legEl.innerHTML = '<div class="sv-empty-hint">—</div>';
      } else {
        legEl.innerHTML = _renderLegendHTML();
        _wireLegend(legEl);
      }
    }

    // UpSet (step 5) — the redirect: rows = per-SV evidence types, top
    // bars = fish counts per evidence combination. Click a bar to populate
    // selected_samples; locus + table re-render scoped to that selection.
    // Reads _state.combinationsLayer (sibling JSON loaded via folder-walk).
    const upEl = _state.rootEl.querySelector('[data-sv-upset]');
    if (upEl) {
      _renderUpSetPanel(upEl);
    }

    // Heatmap (step 6) — sample × SV dosage matrix. Reads
    // _state.supportLayer (sv_support_by_sample_v1.json — third sibling
    // in the per-candidate folder).
    const hmEl = _state.rootEl.querySelector('[data-sv-heatmap]');
    if (hmEl) {
      _renderHeatmapPanel(hmEl);
    }
  }

  // Produce the boundary-summary HTML — two stacked tables.
  function _renderBoundarySummaryHTML(layer) {
    const bs = layer.boundary_summary || {};
    const blocks = [];
    [['left', 'Left boundary'], ['right', 'Right boundary']].forEach(([side, label]) => {
      const block = bs[side];
      if (!block || !block.by_sv_type) {
        blocks.push(`
          <div class="sv-bsblock">
            <div class="sv-bsblock-title sv-bsblock-${side}">${label}</div>
            <div class="sv-empty-hint">no summary in layer</div>
          </div>`);
        return;
      }
      const interval = Array.isArray(block.interval_bp)
        ? `(${_fmtMb(block.interval_bp[0])} – ${_fmtMb(block.interval_bp[1])})`
        : '';
      const types = ['BND', 'INV', 'DEL', 'DUP', 'Other'];
      let totN = 0, totA = 0;
      const rows = types.map(t => {
        const r = block.by_sv_type[t] || { n_total: 0, n_associated_fdr_lt_0_05: 0 };
        totN += r.n_total || 0;
        totA += r.n_associated_fdr_lt_0_05 || 0;
        const c = (SV_TYPE_STYLE[t] || SV_TYPE_STYLE.Other).color;
        return `
          <tr class="sv-bs-row" data-sv-bs-side="${side}" data-sv-bs-type="${t}"
              title="Click to filter table + locus to ${t} SVs at the ${label.toLowerCase()}">
            <td class="sv-bs-type" style="color:${c}">${t}</td>
            <td class="sv-bs-num">${r.n_total || 0}</td>
            <td class="sv-bs-num">${r.n_associated_fdr_lt_0_05 || 0}</td>
          </tr>`;
      }).join('');
      blocks.push(`
        <div class="sv-bsblock">
          <div class="sv-bsblock-title sv-bsblock-${side}">
            ${label} <span class="sv-bsblock-interval">${interval}</span>
          </div>
          <table class="sv-bs-table">
            <thead>
              <tr>
                <th>SV type</th>
                <th title="Total SV calls in this boundary zone"># SVs</th>
                <th title="SV calls with Fisher FDR &lt; 0.05"># Assoc.</th>
              </tr>
            </thead>
            <tbody>
              ${rows}
              <tr class="sv-bs-total">
                <td>Total</td>
                <td class="sv-bs-num">${totN}</td>
                <td class="sv-bs-num">${totA}</td>
              </tr>
            </tbody>
          </table>
        </div>`);
    });
    return blocks.join('');
  }

  // Wire boundary-summary row clicks → filter to SV type + boundary side.
  function _wireBoundarySummary(bsEl) {
    bsEl.querySelectorAll('[data-sv-bs-type]').forEach(tr => {
      tr.addEventListener('click', () => {
        const side = tr.dataset.svBsSide;
        const type = tr.dataset.svBsType;
        // Compose: zone = boundary_500kb (the spec's "Boundary ±500 kb"),
        // sv_type = clicked type. The user can read both from the left rail
        // dropdowns immediately afterward — they're set live.
        _state.filters = {
          ..._state.filters,
          sv_type: type,
          zone: 'boundary_500kb',
          show_only_associated: false, // show Total, not just Assoc.
        };
        _writeFilters(_state.activeCandidateId, _state.filters);
        _state.pageIndex = 0;
        // Optional: snap the view preset to the clicked side close-up so
        // the user sees the boundary they clicked.
        if (side === 'left')  _state.viewPreset = 'left_close';
        if (side === 'right') _state.viewPreset = 'right_close';
        _state._customWindow = null;
        const presetEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-preset]');
        if (presetEl) presetEl.value = _state.viewPreset;
        // Re-render the filter block so the dropdowns reflect the new state
        const filEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-filters]');
        if (filEl) _renderFiltersBlock(filEl);
        _renderLocusStrip();
        _renderTable();
      });
    });
  }

  // Render the legend HTML — two columns: FDR colours, Pattern labels.
  function _renderLegendHTML() {
    const fdrItems = FDR_COLOURS.map(t => `
      <div class="sv-legend-item" data-sv-legend-fdr="${t.max}"
           title="Click to set FDR threshold to ${t.label}">
        <span class="sv-legend-swatch" style="background:${t.color}"></span>
        <span class="sv-legend-text">${_esc(t.label)}</span>
      </div>`).join('');
    const patternItems = Object.keys(PATTERN_LABELS).map(key => {
      const pl = PATTERN_LABELS[key];
      return `
        <div class="sv-legend-item" data-sv-legend-pattern="${key}"
             title="Click to highlight ${pl.display} rows">
          <span class="sv-legend-swatch" style="background:${pl.color}"></span>
          <span class="sv-legend-text" style="color:${pl.color}">${_esc(pl.display)}</span>
        </div>`;
    }).join('');
    return `
      <div class="sv-legend-cols">
        <div class="sv-legend-col">
          <div class="sv-legend-coltitle">Association (FDR)</div>
          ${fdrItems}
        </div>
        <div class="sv-legend-col">
          <div class="sv-legend-coltitle">Pattern labels</div>
          ${patternItems}
        </div>
      </div>`;
  }

  // Wire legend interactions:
  //   - FDR swatch click → set fdr_threshold + show_only_associated=true
  //   - Pattern swatch click → toggle a "highlight pattern" filter (sticky;
  //     adds a CSS class to matching rows). Click again to clear.
  function _wireLegend(legEl) {
    legEl.querySelectorAll('[data-sv-legend-fdr]').forEach(it => {
      it.addEventListener('click', () => {
        const max = parseFloat(it.dataset.svLegendFdr);
        // Map the clicked tier to one of the dropdown's discrete thresholds.
        // The dropdown has 1e-6 / 1e-4 / 0.01 / 0.05 / 0.10 — pick the
        // closest one ≤ the swatch's upper bound. The "> 0.05" tier
        // becomes 0.10 (most permissive).
        const choices = [1e-6, 1e-4, 0.01, 0.05, 0.10];
        let pick = 0.05;
        for (const c of choices) if (c <= max + 1e-15) pick = c;
        if (!isFinite(max)) pick = 0.10;
        _state.filters = {
          ..._state.filters,
          fdr_threshold: pick,
          show_only_associated: true,
        };
        _writeFilters(_state.activeCandidateId, _state.filters);
        _state.pageIndex = 0;
        const filEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-filters]');
        if (filEl) _renderFiltersBlock(filEl);
        _renderLocusStrip();
        _renderTable();
      });
    });
    legEl.querySelectorAll('[data-sv-legend-pattern]').forEach(it => {
      it.addEventListener('click', () => {
        const key = it.dataset.svLegendPattern;
        _state.highlightPattern = (_state.highlightPattern === key) ? null : key;
        _renderTable();
        // Refresh legend so the active swatch shows a ring
        _renderRightRail();
      });
    });
  }

  // ===========================================================================
  // Step 5 — UpSet panel (the redirect from Quentin's step-4 chat turn)
  // ===========================================================================
  // Reads _state.combinationsLayer (sv_evidence_combinations_v1.json,
  // sibling of sv_genotype_counts in the per-candidate folder). Renders:
  //   - Top bars: one bar per combination (height ∝ intersection_size)
  //   - Dot matrix below: rows = evidence types, cols = combinations,
  //     filled dot if evidence type is in combination's `members`
  //   - Right-side mini-bars: per-evidence-type total carrier counts
  //   - Empty state when the layer is missing
  //
  // Click a top bar → _state.selectedSamples = Set(combination.samples).
  // The locus + table pick this up via _renderLocusStrip / _renderTable.
  // Click again on the active bar → deselect.
  // ===========================================================================

  // UpSet layout constants (kept here, not in the global LOCUS_TRACKS list,
  // because UpSet has its own grid math).
  const UPSET = {
    barH:        50,   // top bar area height (px)
    rowH:        14,   // one evidence-type row height
    dotR:         3.5, // filled dot radius
    colW:        18,   // one combination column width
    labelW:     104,   // left column for evidence-type labels
    setBarW:     56,   // right column for per-evidence-type set-size mini-bars
    gapTopBar:    4,   // space between top bars and dot matrix
    padX:         8,
    padY:         6,
  };

  // Validate a candidate sv_evidence_combinations_v1 object. Returns
  // { ok: true } or { ok: false, reason }. Cheap structural check; the
  // renderer is defensive against missing inner fields too.
  function _validateCombinationsLayer(obj) {
    if (!obj || typeof obj !== 'object') return { ok: false, reason: 'NOT_OBJECT' };
    if (obj.format_version !== 'sv_evidence_combinations_v1')
      return { ok: false, reason: 'BAD_FORMAT_VERSION' };
    if (!Array.isArray(obj.evidence_types) || obj.evidence_types.length === 0)
      return { ok: false, reason: 'NO_EVIDENCE_TYPES' };
    if (!Array.isArray(obj.combinations))
      return { ok: false, reason: 'NO_COMBINATIONS' };
    return { ok: true };
  }

  function _renderUpSetPanel(host) {
    if (!host) return;
    const layer = _state.layer;
    const cl = _state.combinationsLayer;

    if (!layer && !cl) {
      host.innerHTML = '<div class="sv-empty-hint">—</div>';
      return;
    }
    if (layer && !cl) {
      host.innerHTML = '<div class="sv-empty-hint sv-upset-emptyhint">' +
        'No <code>sv_evidence_combinations_v1.json</code> loaded.<br>' +
        '<span class="sv-readout-dim">Drop the candidate folder ' +
        '(or this single file) onto the page to populate.</span></div>';
      return;
    }

    const v = _validateCombinationsLayer(cl);
    if (!v.ok) {
      host.innerHTML = '<div class="sv-empty-hint sv-upset-error">' +
        'Invalid combinations layer: ' + _esc(v.reason) + '</div>';
      return;
    }

    const types = cl.evidence_types;
    // Show top N combinations (cap at 20 — anything more is noise in a
    // small right-rail panel; producer should pre-trim too).
    const combos = (cl.combinations || []).slice(0, 20);
    const totals = cl.per_evidence_totals || {};

    // Geometry
    const nCols   = combos.length || 1;
    const nRows   = types.length;
    const matrixH = nRows * UPSET.rowH;
    const matrixW = nCols * UPSET.colW;
    const totalW  = UPSET.padX * 2 + UPSET.labelW + matrixW + UPSET.setBarW + 6;
    const totalH  = UPSET.padY * 2 + UPSET.barH + UPSET.gapTopBar + matrixH + 14;

    // Bar scaling — max intersection_size, clamped to ≥ 1
    let maxI = 1;
    for (const c of combos) maxI = Math.max(maxI, c.intersection_size || 0);
    const barScale = (n) => Math.max(1, ((n || 0) / maxI) * UPSET.barH);

    // Set-size scaling for right-side mini-bars
    let maxS = 1;
    for (const t of types) {
      const v = (totals[t.id] && totals[t.id].n_samples) || 0;
      maxS = Math.max(maxS, v);
    }
    const setScale = (n) => Math.max(1, ((n || 0) / maxS) * (UPSET.setBarW - 8));

    // Build top-bars (one rect per combination)
    const topBars = combos.map((c, i) => {
      const x = UPSET.padX + UPSET.labelW + i * UPSET.colW + UPSET.colW * 0.25;
      const w = UPSET.colW * 0.5;
      const h = barScale(c.intersection_size);
      const y = UPSET.padY + (UPSET.barH - h);
      const isActive = (_state.activeCombinationIndex === i);
      const fill = isActive ? 'var(--accent, #f5a524)' : 'var(--accent-2, #4fa3ff)';
      const stroke = isActive ? 'var(--accent, #f5a524)' : 'transparent';
      return `
        <g class="sv-upset-bar-g" data-sv-upset-bar="${i}">
          <rect class="sv-upset-bar-hit"
                x="${(x - UPSET.colW * 0.5).toFixed(1)}"
                y="${UPSET.padY.toFixed(1)}"
                width="${UPSET.colW.toFixed(1)}"
                height="${(UPSET.barH + UPSET.gapTopBar + matrixH).toFixed(1)}"
                fill="transparent"
                style="cursor:pointer;"/>
          <rect x="${x.toFixed(1)}" y="${y.toFixed(1)}"
                width="${w.toFixed(1)}" height="${h.toFixed(1)}"
                fill="${fill}" stroke="${stroke}" stroke-width="1"/>
          <text x="${(x + w / 2).toFixed(1)}" y="${(y - 2).toFixed(1)}"
                text-anchor="middle" font-size="9.5"
                fill="var(--ink, #e7edf3)">${c.intersection_size}</text>
        </g>`;
    }).join('');

    // Build evidence-type row labels (left column)
    const rowLabels = types.map((t, r) => {
      const y = UPSET.padY + UPSET.barH + UPSET.gapTopBar + r * UPSET.rowH + UPSET.rowH * 0.65;
      const sideColor = t.side === 'left'  ? '#4fa3ff'
                      : t.side === 'right' ? '#bf616a'
                      : 'var(--ink-dim, #8a94a3)';
      return `
        <text x="${(UPSET.padX + UPSET.labelW - 6).toFixed(1)}"
              y="${y.toFixed(1)}"
              text-anchor="end" font-size="10"
              fill="${sideColor}">${_esc(t.label || t.id)}</text>`;
    }).join('');

    // Build dot matrix (rows × cols)
    const matrixCells = [];
    types.forEach((t, r) => {
      // Per-row guide line behind the dots
      const rowY = UPSET.padY + UPSET.barH + UPSET.gapTopBar + r * UPSET.rowH + UPSET.rowH * 0.5;
      matrixCells.push(`
        <line x1="${(UPSET.padX + UPSET.labelW).toFixed(1)}"
              x2="${(UPSET.padX + UPSET.labelW + matrixW).toFixed(1)}"
              y1="${rowY.toFixed(1)}" y2="${rowY.toFixed(1)}"
              stroke="var(--rule, #2a3242)" stroke-width="0.5"/>`);
      combos.forEach((c, i) => {
        const x = UPSET.padX + UPSET.labelW + i * UPSET.colW + UPSET.colW * 0.5;
        const filled = (c.members || []).indexOf(t.id) !== -1;
        const fill = filled ? 'var(--ink, #e7edf3)' : 'var(--ink-dimmer, #5a6472)';
        matrixCells.push(`
          <circle cx="${x.toFixed(1)}" cy="${rowY.toFixed(1)}"
                  r="${filled ? UPSET.dotR : (UPSET.dotR * 0.45)}"
                  fill="${fill}" opacity="${filled ? '1' : '0.5'}"/>`);
      });
      // Connector line within filled dots of a column — drawn per-column
      // below, so we don't repeat per-row here.
    });

    // Vertical connector lines linking filled dots within each column
    combos.forEach((c, i) => {
      const x = UPSET.padX + UPSET.labelW + i * UPSET.colW + UPSET.colW * 0.5;
      // Find first and last filled-row index in this column
      let first = -1, last = -1;
      types.forEach((t, r) => {
        if ((c.members || []).indexOf(t.id) !== -1) {
          if (first < 0) first = r;
          last = r;
        }
      });
      if (first >= 0 && last > first) {
        const y0 = UPSET.padY + UPSET.barH + UPSET.gapTopBar + first * UPSET.rowH + UPSET.rowH * 0.5;
        const y1 = UPSET.padY + UPSET.barH + UPSET.gapTopBar + last  * UPSET.rowH + UPSET.rowH * 0.5;
        matrixCells.push(`
          <line x1="${x.toFixed(1)}" x2="${x.toFixed(1)}"
                y1="${y0.toFixed(1)}" y2="${y1.toFixed(1)}"
                stroke="var(--ink, #e7edf3)" stroke-width="1.2"/>`);
      }
    });

    // Right-side per-evidence-type mini-bars (set sizes)
    const setBars = types.map((t, r) => {
      const n = (totals[t.id] && totals[t.id].n_samples) || 0;
      const y = UPSET.padY + UPSET.barH + UPSET.gapTopBar + r * UPSET.rowH + UPSET.rowH * 0.25;
      const w = setScale(n);
      const x = UPSET.padX + UPSET.labelW + matrixW + 6;
      return `
        <rect x="${x.toFixed(1)}" y="${y.toFixed(1)}"
              width="${w.toFixed(1)}" height="${(UPSET.rowH * 0.5).toFixed(1)}"
              fill="var(--ink-dim, #8a94a3)" opacity="0.55"/>
        <text x="${(x + w + 2).toFixed(1)}"
              y="${(y + UPSET.rowH * 0.5).toFixed(1)}"
              font-size="9" fill="var(--ink-dim, #8a94a3)">${n}</text>`;
    }).join('');

    // Header row above the panel
    const candId = cl.candidate_id || _state.activeCandidateId || '';
    const nSel = _state.selectedSamples ? _state.selectedSamples.size : 0;
    const headerHTML = `
      <div class="sv-upset-header">
        <span class="sv-readout-dim">UpSet · evidence patterns</span>
        ${candId ? `<span class="sv-readout-dim">· ${_esc(candId)}</span>` : ''}
        ${nSel > 0
          ? `<span class="sv-upset-pill" data-sv-action="upset-clear">
               selected: ${nSel} fish · clear
             </span>`
          : ''}
      </div>`;

    host.innerHTML = headerHTML + `
      <svg class="sv-upset-svg" width="${totalW}" height="${totalH}"
           viewBox="0 0 ${totalW} ${totalH}" preserveAspectRatio="xMinYMin meet">
        ${rowLabels}
        ${matrixCells.join('')}
        ${setBars}
        ${topBars}
      </svg>
      <div class="sv-upset-foothint sv-readout-dim">
        click a bar → highlight those fish in locus + table
      </div>`;

    _wireUpSetInteractions(host);
  }

  function _wireUpSetInteractions(host) {
    if (!host) return;
    // Bar click → set/toggle selectedSamples
    host.querySelectorAll('[data-sv-upset-bar]').forEach(g => {
      g.addEventListener('click', () => {
        const i = parseInt(g.dataset.svUpsetBar, 10);
        _onUpSetBarClick(i);
      });
    });
    // Clear-selection pill
    const pill = host.querySelector('[data-sv-action="upset-clear"]');
    if (pill) {
      pill.addEventListener('click', (ev) => {
        ev.stopPropagation();
        _clearSampleSelection();
      });
    }
  }

  function _onUpSetBarClick(i) {
    const cl = _state.combinationsLayer;
    if (!cl || !Array.isArray(cl.combinations) || !cl.combinations[i]) return;
    if (_state.activeCombinationIndex === i) {
      // Toggle off — second click clears
      _clearSampleSelection();
      return;
    }
    const samples = cl.combinations[i].samples || [];
    _state.selectedSamples = new Set(samples);
    _state.activeCombinationIndex = i;
    // Step 6 unlock: recompute the per-SV genotype counts within the
    // selected fish so the table re-counts within selection. No-op when
    // supportLayer isn't loaded; the table falls back to full-cohort
    // counts in that case.
    _recomputeGtCountsView();
    // Re-render everything that visualises selected_samples
    _renderRightRail();   // re-renders UpSet with active bar styling + pill
    _renderLocusStrip();  // dim non-selected SV glyphs
    _renderTable();       // re-count columns within selection (step 6 unlock)
  }

  function _clearSampleSelection() {
    _state.selectedSamples = null;
    _state.activeCombinationIndex = null;
    _state.gtCountsView = null;     // step 6 unlock — reset the recount
    _renderRightRail();
    _renderLocusStrip();
    _renderTable();
  }

  // Step 6 deferred-feature unlock: when an UpSet bar is clicked AND the
  // support layer is loaded, recount each SV's AA/AB/BB/miss counts within
  // the selected fish only, broken down by karyotype group. Result lives
  // in _state.gtCountsView; the table's gc() accessor reads from it.
  //
  // Algorithm:
  //   For each sample row in support_layer.dosage_compact whose ID is in
  //   selectedSamples, look up which karyotype group the row belongs to
  //   (via support.row_groups), then for each sv_id read the dosage char
  //   and tally AA/AB/BB/miss into out[sv_id][group][allele].
  //
  // Edge cases:
  //   - No support layer or no selection -> set view to null (the table
  //     falls back to call.genotype_counts).
  //   - SV in layer but not in support_layer.sv_ids -> not added to view;
  //     gc() falls back transparently.
  //   - Sample row not in any group range -> skipped (defensive; shouldn't
  //     happen with a well-formed layer).
  //
  // Returns the view object (or null) for testing convenience; also
  // assigns to _state.gtCountsView.
  function _recomputeGtCountsView() {
    const supp = _state.supportLayer;
    const sel  = _state.selectedSamples;
    if (!supp || !sel || sel.size === 0) {
      _state.gtCountsView = null;
      return null;
    }
    // Map sample-row index -> karyotype group label
    const rowGroup = new Array(supp.samples.length).fill(null);
    if (supp.row_groups) {
      for (const g of ['H1/H1', 'H1/H2', 'H2/H2']) {
        const r = supp.row_groups[g];
        if (Array.isArray(r) && r.length === 2) {
          for (let i = r[0]; i <= r[1] && i < rowGroup.length; i++) rowGroup[i] = g;
        }
      }
    }
    // Char -> allele bucket
    const allele = { '0': 'AA', '1': 'AB', '2': 'BB', '.': 'miss' };
    const view = {};
    const groupNs = { 'H1/H1': 0, 'H1/H2': 0, 'H2/H2': 0 };
    // Initialise view shells per SV
    for (const svId of supp.sv_ids) {
      view[svId] = {
        'H1/H1': { AA: 0, AB: 0, BB: 0, miss: 0 },
        'H1/H2': { AA: 0, AB: 0, BB: 0, miss: 0 },
        'H2/H2': { AA: 0, AB: 0, BB: 0, miss: 0 },
      };
    }
    for (let r = 0; r < supp.samples.length; r++) {
      const sid = supp.samples[r];
      if (!sel.has(sid)) continue;
      const g = rowGroup[r];
      if (!g) continue;
      groupNs[g]++;
      const row = supp.dosage_compact[r];
      if (typeof row !== 'string') continue;
      for (let c = 0; c < supp.sv_ids.length; c++) {
        const ch = row.charAt(c);
        const a = allele[ch];
        if (a) view[supp.sv_ids[c]][g][a]++;
      }
    }
    view._groupNs = groupNs;
    _state.gtCountsView = view;
    return view;
  }

  // ===========================================================================
  // Step 6 — Sample × SV dosage heatmap
  // ===========================================================================
  // Reads _state.supportLayer (sv_support_by_sample_v1.json, third sibling
  // in the per-candidate folder). Two views:
  //   - 'compact' — right-rail thumbnail under the UpSet panel (cells small)
  //   - 'large'   — full-width overlay panel above the locus strip
  //
  // The dosage matrix is encoded as one string per sample:
  //   '0' = AA (REF), '1' = AB (HET), '2' = BB (ALT/INV), '.' = missing.
  // Char-based parsing is fast (charCodeAt) and the producer side is
  // trivial in Python (''.join(map(...))).
  //
  // Cell hover → tooltip with (sample_id, sv_id, dosage). Cell click →
  // highlight that SV in the table + locus.
  // ===========================================================================

  const HEATMAP = {
    compact: {
      cellW: 1.5,    // each sample column is this many px wide
      cellH: 8,      // each SV row is this tall
      maxRows: 12,   // truncate rows if more SVs than this fit
      labelW: 86,
      gapY: 1,
    },
    large: {
      cellW: 3,
      cellH: 14,
      maxRows: 100,
      labelW: 110,
      gapY: 0,
    },
    colours: {
      // Dosage → CSS variable. Stable across light + dark themes.
      '0': 'var(--panel-3, #2a3242)',     // AA / REF: muted
      '1': 'var(--accent-2, #4fa3ff)',    // AB / HET: blue
      '2': 'var(--bad,     #e0555c)',     // BB / ALT (inversion): red
      '.': 'transparent',                  // missing: blank
    },
    // Light-mode AA needs a lighter shade than panel-3 for contrast
    coloursLightAA: '#cdd5df',
  };

  // Validate sv_support_by_sample_v1. Cheap structural; defensive renderers
  // tolerate empty inner blocks.
  function _validateSupportLayer(obj) {
    if (!obj || typeof obj !== 'object')               return { ok: false, reason: 'NOT_OBJECT' };
    if (obj.format_version !== 'sv_support_by_sample_v1')
      return { ok: false, reason: 'BAD_FORMAT_VERSION' };
    if (!Array.isArray(obj.samples) || obj.samples.length === 0)
      return { ok: false, reason: 'NO_SAMPLES' };
    if (!Array.isArray(obj.sv_ids) || obj.sv_ids.length === 0)
      return { ok: false, reason: 'NO_SV_IDS' };
    if (!Array.isArray(obj.dosage_compact) ||
        obj.dosage_compact.length !== obj.samples.length)
      return { ok: false, reason: 'DOSAGE_ROW_COUNT_MISMATCH' };
    return { ok: true };
  }

  // Look up dosage for (rowIdx, colIdx). Returns char '0'/'1'/'2'/'.' or null.
  function _supportCell(supp, rowIdx, colIdx) {
    if (!supp || !supp.dosage_compact) return null;
    const row = supp.dosage_compact[rowIdx];
    if (typeof row !== 'string') return null;
    const ch = row.charAt(colIdx);
    return ch || null;
  }

  function _renderHeatmapPanel(host) {
    if (!host) return;
    const supp = _state.supportLayer;
    if (!supp) {
      host.innerHTML = '<div class="sv-empty-hint sv-heatmap-emptyhint">' +
        'No <code>sv_support_by_sample_v1.json</code> loaded.<br>' +
        '<span class="sv-readout-dim">Drop the candidate folder to populate the heatmap.</span></div>';
      return;
    }
    const v = _validateSupportLayer(supp);
    if (!v.ok) {
      host.innerHTML = '<div class="sv-empty-hint sv-heatmap-error">' +
        'Invalid support layer: ' + _esc(v.reason) + '</div>';
      return;
    }

    // Compact view (the default for the right rail). The "expand" button
    // toggles to the large overlay.
    if (_state.heatmapView === 'large') {
      _renderHeatmapLarge();
      // Compact host shows a dim placeholder while large is open
      host.innerHTML = '<div class="sv-empty-hint sv-readout-dim">' +
        'Heatmap open in main area — close it to return.</div>';
      return;
    }
    host.innerHTML = _renderHeatmapHTML(supp, 'compact');
    _wireHeatmapInteractions(host, supp, 'compact');
  }

  // Build the heatmap SVG + header. `mode` ∈ {'compact', 'large'}.
  // Filters SVs by the active filter set so the heatmap stays in sync
  // with the table; samples ordered by row_groups (H1/H1 → H1/H2 → H2/H2).
  function _renderHeatmapHTML(supp, mode) {
    const cfg = HEATMAP[mode];
    const layer = _state.layer;

    // Filter SVs: keep only ones that pass the active filter set AND are
    // in the support matrix. Cap at cfg.maxRows for performance / sanity.
    let svIds = supp.sv_ids.slice(0);
    if (layer && Array.isArray(layer.sv_calls)) {
      const visibleSet = new Set(_applyFilters(layer.sv_calls, _state.filters)
                                   .map(c => c.sv_id));
      svIds = svIds.filter(id => visibleSet.has(id));
    }
    let truncated = false;
    if (svIds.length > cfg.maxRows) {
      svIds = svIds.slice(0, cfg.maxRows);
      truncated = true;
    }
    const colIdxByRow = svIds.map(id => supp.sv_ids.indexOf(id));

    const samples = supp.samples;
    const nSamples = samples.length;
    const nRows = svIds.length;

    // Group separators (vertical lines between karyotype groups)
    const groupSeps = [];
    if (supp.row_groups) {
      const order = ['H1/H1', 'H1/H2', 'H2/H2'];
      let cumul = 0;
      for (const g of order) {
        const range = supp.row_groups[g];
        if (Array.isArray(range) && range.length === 2) {
          cumul = range[1] + 1;
          if (cumul < nSamples) groupSeps.push(cumul);
        }
      }
    }

    const matrixW = nSamples * cfg.cellW;
    const matrixH = nRows * (cfg.cellH + cfg.gapY);
    const totalW = cfg.labelW + matrixW + 12;
    const totalH = matrixH + 24;

    // Build cells. One <rect> per non-AA cell (saves DOM nodes for the
    // common-mode 0 background — we render a single AA-coloured backdrop
    // rect first, then overlay HET/HOM/miss rects only).
    const cells = [];
    const aaCol = HEATMAP.colours['0'];

    for (let r = 0; r < nRows; r++) {
      const colIdx = colIdxByRow[r];
      const yRow = r * (cfg.cellH + cfg.gapY);
      // AA backdrop for the whole row
      cells.push(`<rect x="${cfg.labelW}" y="${yRow}"
                       width="${matrixW}" height="${cfg.cellH}"
                       fill="${aaCol}" data-sv-row="${r}"/>`);
      // Overlay HET / HOM / miss
      for (let s = 0; s < nSamples; s++) {
        const ch = _supportCell(supp, s, colIdx);
        if (!ch || ch === '0') continue;
        const fill = HEATMAP.colours[ch] || 'transparent';
        if (fill === 'transparent') continue;
        cells.push(`<rect x="${(cfg.labelW + s * cfg.cellW).toFixed(2)}"
                         y="${yRow}" width="${cfg.cellW}" height="${cfg.cellH}"
                         fill="${fill}"/>`);
      }
      // Row label (SV id)
      const isHi = (_state.highlightedSvId === svIds[r]);
      const labelClass = isHi ? 'sv-hm-rowlabel sv-hm-rowlabel-hi' : 'sv-hm-rowlabel';
      cells.push(
        `<text x="${(cfg.labelW - 4)}" y="${(yRow + cfg.cellH * 0.75)}"
              text-anchor="end" font-size="9.5"
              class="${labelClass}"
              data-sv-hm-svid="${_esc(svIds[r])}">${_esc(svIds[r])}</text>`
      );
    }
    // Group separators
    const sepLines = groupSeps.map(s => {
      const x = cfg.labelW + s * cfg.cellW;
      return `<line x1="${x.toFixed(2)}" x2="${x.toFixed(2)}" y1="0"
                    y2="${matrixH}" stroke="var(--ink-dim, #8a94a3)"
                    stroke-width="0.7" opacity="0.7"/>`;
    }).join('');
    // Group labels along the bottom
    const groupLabels = [];
    if (supp.row_groups) {
      const order = ['H1/H1', 'H1/H2', 'H2/H2'];
      for (const g of order) {
        const range = supp.row_groups[g];
        if (!Array.isArray(range)) continue;
        const xMid = cfg.labelW + ((range[0] + range[1]) / 2) * cfg.cellW;
        groupLabels.push(
          `<text x="${xMid.toFixed(2)}" y="${(matrixH + 12).toFixed(2)}"
                font-size="9" fill="var(--ink-dim, #8a94a3)"
                text-anchor="middle">${g} (n=${range[1] - range[0] + 1})</text>`
        );
      }
    }
    // Hover-cell highlight rect (drawn last so it overlays everything)
    let hoverRect = '';
    if (_state.heatmapHoverCell) {
      const { rowIdx, colIdx } = _state.heatmapHoverCell;
      // colIdx in this view's display order
      const dispRowIdx = svIds.indexOf(supp.sv_ids[colIdx]);
      if (dispRowIdx >= 0 && rowIdx >= 0 && rowIdx < nSamples) {
        const x = cfg.labelW + rowIdx * cfg.cellW;
        const y = dispRowIdx * (cfg.cellH + cfg.gapY);
        hoverRect = `<rect x="${x.toFixed(2)}" y="${y}"
                          width="${cfg.cellW}" height="${cfg.cellH}"
                          fill="none" stroke="var(--accent, #f5a524)"
                          stroke-width="1" pointer-events="none"/>`;
      }
    }

    // Active-cell tooltip line (sample id, sv id, dosage)
    let tipHtml = '';
    if (_state.heatmapHoverCell) {
      const { rowIdx, colIdx } = _state.heatmapHoverCell;
      const sId = samples[rowIdx];
      const svId = supp.sv_ids[colIdx];
      const ch = _supportCell(supp, rowIdx, colIdx);
      const dosLabel = { '0': 'AA (REF)', '1': 'AB (HET)', '2': 'BB (HOM-ALT)', '.': 'missing' }[ch] || '?';
      tipHtml = `<div class="sv-hm-tip">
        <code>${_esc(sId)}</code> · <code>${_esc(svId)}</code>
        · <strong>${dosLabel}</strong>
      </div>`;
    }

    // Header with title + optional "expand" button
    const expandBtn = mode === 'compact'
      ? `<button class="sv-hm-expand" data-sv-hm-action="expand"
                 title="Open large heatmap in main area">↗ expand</button>`
      : `<button class="sv-hm-expand" data-sv-hm-action="collapse"
                 title="Close large heatmap">close</button>`;

    const truncatedNote = truncated
      ? `<span class="sv-readout-dim">· showing top ${cfg.maxRows} of ${supp.sv_ids.length}</span>`
      : '';

    const legendInline = `
      <span class="sv-hm-leg">
        <span class="sv-hm-sw" style="background:${HEATMAP.colours['0']}"></span> AA
        <span class="sv-hm-sw" style="background:${HEATMAP.colours['1']}"></span> AB
        <span class="sv-hm-sw" style="background:${HEATMAP.colours['2']}"></span> BB
      </span>`;

    return `
      <div class="sv-hm-header">
        <span class="sv-readout-dim">Dosage heatmap</span>
        ${legendInline}
        ${truncatedNote}
        ${expandBtn}
      </div>
      ${tipHtml}
      <div class="sv-hm-svgwrap">
        <svg class="sv-hm-svg" width="${totalW}" height="${totalH}"
             viewBox="0 0 ${totalW} ${totalH}" preserveAspectRatio="xMinYMin meet"
             data-sv-hm-svg="1">
          ${cells.join('')}
          ${sepLines}
          ${groupLabels.join('')}
          ${hoverRect}
        </svg>
      </div>`;
  }

  // Render the large overlay panel (full main area). Clicked via expand
  // button. Closes on collapse / Esc / clicking the close button.
  function _renderHeatmapLarge() {
    if (!_state.rootEl) return;
    let overlay = _state.rootEl.querySelector('[data-sv-hm-overlay]');
    if (!overlay) {
      overlay = document.createElement('div');
      overlay.setAttribute('data-sv-hm-overlay', '1');
      overlay.className = 'sv-hm-overlay';
      _state.rootEl.querySelector('.sv-evidence-page').appendChild(overlay);
    }
    const supp = _state.supportLayer;
    if (!supp) {
      overlay.innerHTML = '<div class="sv-empty-hint">no support layer</div>';
      return;
    }
    overlay.innerHTML = _renderHeatmapHTML(supp, 'large');
    _wireHeatmapInteractions(overlay, supp, 'large');
  }

  // Wire hover + click for one heatmap host (compact OR large overlay).
  function _wireHeatmapInteractions(host, supp, mode) {
    const cfg = HEATMAP[mode];
    const svg = host.querySelector('[data-sv-hm-svg]');

    if (svg) {
      svg.addEventListener('mousemove', (ev) => {
        const rect = svg.getBoundingClientRect();
        const px = ev.clientX - rect.left;
        const py = ev.clientY - rect.top;
        if (px < cfg.labelW) { _heatmapClearHover(host); return; }
        const colInDisp = Math.floor((px - cfg.labelW) / cfg.cellW);
        const rowInDisp = Math.floor(py / (cfg.cellH + cfg.gapY));
        if (colInDisp < 0 || colInDisp >= supp.samples.length) {
          _heatmapClearHover(host); return;
        }
        // Translate display row → matrix col index (sv_ids index)
        // We need the same display order used in _renderHeatmapHTML.
        // Cheapest: re-derive svIds list. Or stash the mapping.
        // Actually, the SVG cells carry data-sv-row but not the SV index;
        // re-compute from the layer + filters + caps.
        let svIds = supp.sv_ids.slice(0);
        if (_state.layer && Array.isArray(_state.layer.sv_calls)) {
          const visibleSet = new Set(_applyFilters(_state.layer.sv_calls, _state.filters)
                                       .map(c => c.sv_id));
          svIds = svIds.filter(id => visibleSet.has(id));
        }
        if (svIds.length > cfg.maxRows) svIds = svIds.slice(0, cfg.maxRows);
        if (rowInDisp < 0 || rowInDisp >= svIds.length) {
          _heatmapClearHover(host); return;
        }
        const colIdx = supp.sv_ids.indexOf(svIds[rowInDisp]);
        // colIdx here is the sv_ids matrix index; rowInDisp's sample index
        // is colInDisp (because px is along the sample axis).
        _state.heatmapHoverCell = { rowIdx: colInDisp, colIdx: colIdx };
        // Re-render only the host (not the whole right rail) — partial
        // update that doesn't churn the UpSet panel.
        host.innerHTML = _renderHeatmapHTML(supp, mode);
        _wireHeatmapInteractions(host, supp, mode);
      });
      svg.addEventListener('mouseleave', () => _heatmapClearHover(host));

      svg.addEventListener('click', (ev) => {
        // Click on a row label or anywhere in a row → highlight that SV in
        // table + locus.
        const tgt = ev.target;
        if (tgt && tgt.dataset && tgt.dataset.svHmSvid) {
          _state.highlightedSvId = tgt.dataset.svHmSvid;
          _renderTable();
          _renderLocusStrip();
          _scrollTableRowIntoView(tgt.dataset.svHmSvid);
          return;
        }
        // Click on a cell → infer the SV by row position
        if (_state.heatmapHoverCell) {
          const { colIdx } = _state.heatmapHoverCell;
          const svId = supp.sv_ids[colIdx];
          if (svId) {
            _state.highlightedSvId = svId;
            _renderTable();
            _renderLocusStrip();
            _scrollTableRowIntoView(svId);
          }
        }
      });
    }
    // Expand / collapse buttons
    host.querySelectorAll('[data-sv-hm-action]').forEach(btn => {
      btn.addEventListener('click', (ev) => {
        ev.stopPropagation();
        const action = btn.dataset.svHmAction;
        if (action === 'expand') {
          _state.heatmapView = 'large';
          _renderRightRail();
          _renderHeatmapLarge();
        } else if (action === 'collapse') {
          _state.heatmapView = 'compact';
          // Remove the overlay
          const ov = _state.rootEl && _state.rootEl.querySelector('[data-sv-hm-overlay]');
          if (ov && ov.parentNode) ov.parentNode.removeChild(ov);
          _renderRightRail();
        }
      });
    });
  }

  function _heatmapClearHover(host) {
    if (_state.heatmapHoverCell == null) return;
    _state.heatmapHoverCell = null;
    const supp = _state.supportLayer;
    if (!supp) return;
    // Detect compact vs large by host's class
    const mode = host.classList && host.classList.contains('sv-hm-overlay') ? 'large' : 'compact';
    host.innerHTML = _renderHeatmapHTML(supp, mode);
    _wireHeatmapInteractions(host, supp, mode);
  }

  // ===========================================================================
  // Locus track strip (spec §3.4) — SVG-based, self-contained
  // ===========================================================================
  //
  // The mockup shows seven vertically stacked tracks sharing one bp x-axis:
  //   1. Zone bar (5 coloured segments)
  //   2. Coordinate axis (bp ticks every 0.5 Mb, dashed verticals at boundaries)
  //   3. Dosage heatmap                — placeholder until atlas-state wiring
  //   4. PCA local (PC1/PC2/PC3)       — placeholder until atlas-state wiring
  //   5. SV calls (one row per type)   — RENDERED FROM LAYER (this is the new track)
  //   6. Genes (GFF strand arrows)     — placeholder until atlas-state wiring
  //   7. Repeat density (mini bars)    — placeholder until atlas-state wiring
  //
  // Tracks 1, 2, 5 ship in step 3 because they're fully self-contained from
  // the sv_genotype_counts_v1 layer. Tracks 3/4/6/7 require data that lives
  // elsewhere in atlas state (state.candidatePCA, state.geneTrack,
  // state.repeatDensity, etc.) — for the local-first workflow Quentin can
  // turn them on when he has them; for now they render as labelled
  // placeholders rather than throwing or hiding.

  // Track configuration. h = pixel height, key = state key (or null for
  // self-contained). Order is the visual stack top-to-bottom.
  const LOCUS_TRACKS = [
    { id: 'zones',        label: 'Zones',                h: 22, kind: 'self' },
    { id: 'axis',         label: '',                     h: 22, kind: 'self' },
    { id: 'dosage',       label: 'Dosage (cap 200)',     h: 64, kind: 'state', stateKey: 'dosageMatrix' },
    { id: 'pca',          label: 'PCA (local)',          h: 56, kind: 'state', stateKey: 'candidatePCA' },
    { id: 'sv_calls',     label: 'SV calls (all types)', h: 60, kind: 'self' },
    { id: 'genes',        label: 'Genes (GFF)',          h: 28, kind: 'state', stateKey: 'geneTrack' },
    { id: 'repeats',      label: 'Repeat density',       h: 28, kind: 'state', stateKey: 'repeatDensity' },
  ];

  // Layout constants
  const LOCUS_PADDING_X = 12;
  const LOCUS_LABEL_W   = 110;
  const LOCUS_TRACK_GAP = 4;

  // Compute the [lo, hi] bp window for a view preset (spec §3.2).
  // Defensive against missing zone definitions; always returns finite bounds.
  // The '_custom' preset reads _state._customWindow (set by zoom in/out).
  function _windowForPreset(layer, preset, z) {
    z = z || (layer && layer.zone_definitions_bp) || {};
    const fullLo = (z.left_flank  && z.left_flank[0])  != null ? z.left_flank[0]
                    : ((layer && layer.boundary_left_bp  != null) ? layer.boundary_left_bp  - 500000 : 0);
    const fullHi = (z.right_flank && z.right_flank[1]) != null ? z.right_flank[1]
                    : ((layer && layer.boundary_right_bp != null) ? layer.boundary_right_bp + 500000 : 0);
    if (preset === '_custom' && _state._customWindow) {
      return { lo: _state._customWindow.lo, hi: _state._customWindow.hi, preset: preset };
    }
    if (preset === 'left_close' && layer && layer.boundary_left_bp != null) {
      return { lo: layer.boundary_left_bp  - 250000, hi: layer.boundary_left_bp  + 250000, preset: preset };
    }
    if (preset === 'right_close' && layer && layer.boundary_right_bp != null) {
      return { lo: layer.boundary_right_bp - 250000, hi: layer.boundary_right_bp + 250000, preset: preset };
    }
    return { lo: fullLo, hi: fullHi, preset: 'default' };
  }

  // Render the locus strip into the [data-sv-locus-strip] container.
  // Idempotent — full innerHTML rebuild every call (cheap; SVG is bounded).
  function _renderLocusStrip() {
    if (!_state.rootEl) return;
    const host = _state.rootEl.querySelector('[data-sv-locus-strip]');
    if (!host) return;
    const layer = _state.layer;
    if (!layer || !layer.zone_definitions_bp) {
      // Don't blow away the empty-state card; just leave it.
      return;
    }

    // Determine the x-axis domain from the active view preset (spec §3.2):
    //   default      — full window: left_flank.start → right_flank.end
    //   left_close   — boundary_left_bp ± 250 kb
    //   right_close  — boundary_right_bp ± 250 kb
    // If preset zones are missing, fall back to ± 500 kb of boundaries.
    const z = layer.zone_definitions_bp;
    const win = _windowForPreset(layer, _state.viewPreset, z);
    const lo = win.lo, hi = win.hi;
    if (!isFinite(lo) || !isFinite(hi) || hi <= lo) {
      host.innerHTML = '<div class="sv-empty-hint" style="padding:24px">' +
                       'Cannot compute locus window — boundary or zone metadata invalid.</div>';
      return;
    }

    // Width responds to container size; height = sum(track h) + gaps.
    const containerW = Math.max(640, host.clientWidth || 0);
    const drawW = containerW - LOCUS_PADDING_X * 2 - LOCUS_LABEL_W - 12;
    const totalH = LOCUS_TRACKS.reduce((a, t) => a + t.h + LOCUS_TRACK_GAP, 0);

    // bp → x pixel mapper, used by every track.
    const x = bp => ((bp - lo) / (hi - lo)) * drawW;
    // Inverse: x pixel → bp (used by hover/click cursor logic)
    const xInv = px => lo + (px / drawW) * (hi - lo);
    // Stash on _state so the hotkey handler can read it without re-deriving
    _state._lastWindow = { lo, hi, drawW, totalH, x, xInv };

    // Build each track group; stack them vertically.
    let yCursor = 0;
    const trackSvgs = LOCUS_TRACKS.map(t => {
      const svg = _renderLocusTrack(t, layer, lo, hi, drawW, t.h, x);
      const block = `
        <div class="sv-locus-track" style="height:${t.h}px;">
          <span class="sv-locus-track-label">${_esc(t.label)}</span>
          <div class="sv-locus-track-svg">${svg}</div>
        </div>`;
      yCursor += t.h + LOCUS_TRACK_GAP;
      return block;
    }).join('');

    // Boundary verticals — dashed lines at boundary_left_bp (blue) and
    // boundary_right_bp (red), spanning all tracks. Plus the cursor (pale
    // accent vertical) and sticky E/F markers if set. All rendered as a
    // single overlay SVG positioned absolute over the stack.
    //
    // The cursor + markers are interactive: the overlay SVG receives mouse
    // events (pointer-events:auto on the inner <rect> hit-target only;
    // glyphs in the underlying tracks still receive their own clicks
    // because the overlay's other elements have pointer-events:none).
    const xL = x(layer.boundary_left_bp);
    const xR = x(layer.boundary_right_bp);
    const cBp = _state.cursorBp;
    const cInside = (cBp != null && cBp >= lo && cBp <= hi);
    const xCur = cInside ? x(cBp) : null;
    const mL = _state.markerLeftBp;
    const mR = _state.markerRightBp;
    const xML = (mL != null && mL >= lo && mL <= hi) ? x(mL) : null;
    const xMR = (mR != null && mR >= lo && mR <= hi) ? x(mR) : null;

    // Selection rect (step 4.5) — drawn first so the boundary verticals
    // and cursor render on top. Two layers: sticky selection (committed)
    // and in-progress drag preview.
    const sel = _state.selection;
    const drag = _state._selectDrag;
    let selRect = '';
    if (sel && sel.startBp != null && sel.endBp != null
        && sel.endBp >= lo && sel.startBp <= hi) {
      const sx0 = x(Math.max(sel.startBp, lo));
      const sx1 = x(Math.min(sel.endBp,   hi));
      selRect = `
        <rect x="${sx0.toFixed(1)}" y="0"
              width="${(sx1 - sx0).toFixed(1)}" height="${totalH}"
              fill="var(--accent, #f5a524)" opacity="0.10"
              style="pointer-events:none;"/>
        <line x1="${sx0.toFixed(1)}" x2="${sx0.toFixed(1)}" y1="0" y2="${totalH}"
              stroke="var(--accent, #f5a524)" stroke-width="1" opacity="0.7"
              style="pointer-events:none;"/>
        <line x1="${sx1.toFixed(1)}" x2="${sx1.toFixed(1)}" y1="0" y2="${totalH}"
              stroke="var(--accent, #f5a524)" stroke-width="1" opacity="0.7"
              style="pointer-events:none;"/>`;
    }
    let dragRect = '';
    if (drag) {
      const dlo = Math.min(drag.startBp, drag.currentBp);
      const dhi = Math.max(drag.startBp, drag.currentBp);
      const dx0 = x(Math.max(dlo, lo));
      const dx1 = x(Math.min(dhi, hi));
      dragRect = `
        <rect x="${dx0.toFixed(1)}" y="0"
              width="${(dx1 - dx0).toFixed(1)}" height="${totalH}"
              fill="var(--accent, #f5a524)" opacity="0.20"
              stroke="var(--accent, #f5a524)" stroke-width="1"
              stroke-dasharray="2 2"
              style="pointer-events:none;"/>`;
    }

    const overlay = `
      <svg class="sv-locus-overlay" width="${drawW}" height="${totalH}"
           viewBox="0 0 ${drawW} ${totalH}" preserveAspectRatio="none">
        <!-- transparent hit-rect captures hover/click for cursor logic -->
        <rect class="sv-locus-hit" x="0" y="0" width="${drawW}" height="${totalH}"
              fill="transparent" data-sv-locus-hit="1"
              style="pointer-events:auto; cursor: ${_state.selectMode === 'select' ? 'col-resize' : 'crosshair'};"/>
        ${selRect}
        ${dragRect}
        <!-- inferred boundary verticals (always shown) -->
        <line x1="${xL.toFixed(1)}" x2="${xL.toFixed(1)}" y1="0" y2="${totalH}"
              stroke="#4fa3ff" stroke-width="1" stroke-dasharray="3 3" opacity="0.85"
              style="pointer-events:none;"/>
        <line x1="${xR.toFixed(1)}" x2="${xR.toFixed(1)}" y1="0" y2="${totalH}"
              stroke="#bf616a" stroke-width="1" stroke-dasharray="3 3" opacity="0.85"
              style="pointer-events:none;"/>
        <text x="${xL.toFixed(1)}" y="10" fill="#4fa3ff" font-size="9.5"
              text-anchor="start" dx="3" style="pointer-events:none;">${_fmtMb(layer.boundary_left_bp)} (left)</text>
        <text x="${xR.toFixed(1)}" y="10" fill="#bf616a" font-size="9.5"
              text-anchor="end" dx="-3" style="pointer-events:none;">${_fmtMb(layer.boundary_right_bp)} (right)</text>
        ${xML != null ? `
          <!-- E marker — left visual pin -->
          <line x1="${xML.toFixed(1)}" x2="${xML.toFixed(1)}" y1="0" y2="${totalH}"
                stroke="#3cc08a" stroke-width="1.5" opacity="0.9"
                style="pointer-events:none;"/>
          <text x="${xML.toFixed(1)}" y="${(totalH - 2).toFixed(1)}" fill="#3cc08a"
                font-size="9.5" text-anchor="start" dx="3"
                style="pointer-events:none; font-weight:600;">E ${_fmtMb(mL)}</text>` : ''}
        ${xMR != null ? `
          <!-- F marker — right visual pin -->
          <line x1="${xMR.toFixed(1)}" x2="${xMR.toFixed(1)}" y1="0" y2="${totalH}"
                stroke="#d08770" stroke-width="1.5" opacity="0.9"
                style="pointer-events:none;"/>
          <text x="${xMR.toFixed(1)}" y="${(totalH - 2).toFixed(1)}" fill="#d08770"
                font-size="9.5" text-anchor="end" dx="-3"
                style="pointer-events:none; font-weight:600;">F ${_fmtMb(mR)}</text>` : ''}
        ${xCur != null ? `
          <!-- live cursor -->
          <line x1="${xCur.toFixed(1)}" x2="${xCur.toFixed(1)}" y1="0" y2="${totalH}"
                stroke="var(--accent-2, #4fa3ff)" stroke-width="1" opacity="0.55"
                style="pointer-events:none;"/>` : ''}
      </svg>`;

    // Legend strip (matches mockup): "SV types: ▼ BND  ▼ INV  ▼ DEL  ▼ DUP  ▲ Other"
    const svLegend = ['BND','INV','DEL','DUP','Other'].map(t => {
      const s = SV_TYPE_STYLE[t];
      return `<span class="sv-locus-legend-item">
        <span class="sv-locus-legend-glyph" style="color:${s.color}">${s.glyph}</span>
        ${t}
      </span>`;
    }).join('');

    // Cursor readout — small floating panel at top-right of the strip.
    // Shows: cursor bp (Mb), nearest SV id, distance to nearest SV, hint
    // about hotkeys. When cursor is null, shows hotkey hint only.
    const readout = _renderLocusReadout(layer, lo, hi);

    host.innerHTML = `
      <div class="sv-locus-wrap" style="padding:8px ${LOCUS_PADDING_X}px;">
        <div class="sv-locus-readout">${readout}</div>
        <div class="sv-locus-stack" style="position:relative;">
          ${trackSvgs}
          <div class="sv-locus-overlay-wrap"
               style="position:absolute; left:${LOCUS_LABEL_W + 4}px; top:0;
                      width:${drawW}px; height:${totalH}px;">
            ${overlay}
          </div>
        </div>
        <div class="sv-locus-legend">
          <span class="sv-locus-legend-label">SV types:</span> ${svLegend}
        </div>
      </div>
    `;
    // Switch the locus-strip container from centered-empty-state mode to
    // stretched mode when rendering tracks (CSS class toggle).
    host.classList.add('sv-locus-has-content');

    // Wire SV-icon click + cursor hover/click handshake.
    _wireLocusInteractions(host);
  }

  // Render the small readout above the locus stack. Three lines:
  //   cursor: <Mb> bp  ·  nearest: <SV id> (Δ <bp>)
  //   markers: E <Mb> · F <Mb>  (omitted if neither set)
  //   hotkeys: ←/→ cursor · Shift+←/→ jump SV · ↑/↓ rows · Shift+↑/↓ SV type · E/F pin · Enter highlight · Esc clear
  function _renderLocusReadout(layer, lo, hi) {
    const cBp = _state.cursorBp;
    const lines = [];
    if (cBp == null) {
      lines.push('<span class="sv-readout-dim">cursor: —</span>' +
                 '<span class="sv-readout-sep">·</span>' +
                 '<span class="sv-readout-dim">click or hover the strip</span>');
    } else {
      const near = _nearestSV(cBp);
      const nearStr = near
        ? `nearest: <strong style="color:${(SV_TYPE_STYLE[near.sv.sv_type] || SV_TYPE_STYLE.Other).color}">${_esc(near.sv.sv_id)}</strong> (${near.sv.sv_type})`
          + ` <span class="sv-readout-dim">· Δ ${_fmtBp(near.dist)}</span>`
        : '<span class="sv-readout-dim">no SVs in window</span>';
      lines.push(`<span class="sv-readout-strong">cursor: ${_fmtMb(cBp)}</span>` +
                 '<span class="sv-readout-sep">·</span>' + nearStr);
    }
    if (_state.markerLeftBp != null || _state.markerRightBp != null) {
      const e = _state.markerLeftBp  != null ? `<span style="color:#3cc08a">E ${_fmtMb(_state.markerLeftBp)}</span>`  : '';
      const f = _state.markerRightBp != null ? `<span style="color:#d08770">F ${_fmtMb(_state.markerRightBp)}</span>` : '';
      lines.push(['<span class="sv-readout-dim">markers:</span>', e, f].filter(Boolean).join(' '));
    }
    // Selection info (step 4.5) — when a sticky region is active
    const sel = _state.selection;
    if (sel && sel.startBp != null && sel.endBp != null) {
      const span = sel.endBp - sel.startBp;
      lines.push(
        '<span class="sv-readout-dim">region:</span> ' +
        `<span style="color:var(--accent, #f5a524); font-weight:600;">${_fmtMb(sel.startBp)} – ${_fmtMb(sel.endBp)}</span>` +
        ` <span class="sv-readout-dim">(${(span / 1e6).toFixed(2)} Mb)</span> ` +
        '<button class="sv-clear-selection" data-sv-action="clear-selection"' +
        ' title="Clear region selection (also clears with Esc)">clear</button>'
      );
    }
    // Sample selection info (step 5) — when an UpSet bar is active
    if (_state.selectedSamples && _state.selectedSamples.size > 0) {
      const cl = _state.combinationsLayer;
      let comboDesc = '';
      if (_state.activeCombinationIndex != null && cl &&
          cl.combinations && cl.combinations[_state.activeCombinationIndex]) {
        const cm = cl.combinations[_state.activeCombinationIndex];
        comboDesc = ` <span class="sv-readout-dim">[${(cm.members || []).join(' + ')}]</span>`;
      }
      lines.push(
        '<span class="sv-readout-dim">samples:</span> ' +
        `<span style="color:var(--accent, #f5a524); font-weight:600;">` +
          `${_state.selectedSamples.size} fish</span>${comboDesc} ` +
        '<button class="sv-clear-selection" data-sv-action="clear-samples"' +
        ' title="Clear UpSet sample selection (also via UpSet panel)">clear</button>'
      );
    }
    lines.push('<span class="sv-readout-dim sv-readout-hotkeys">' +
               '←/→ cursor · Shift+←/→ jump SV · ↑/↓ rows · Shift+↑/↓ SV type · E/F pin · Enter highlight · Esc clear' +
               '</span>');
    return lines.join('<br>');
  }

  // Find the SV (within the active filter set) whose position_bp is closest
  // to the cursor. Returns { sv, dist } or null if no candidates.
  // Helper: does this SV have at least one carrier (HET or HOM-ALT) in the
  // current sample selection (UpSet bar click)? Returns false if no support
  // layer or no selection (so we don't dim anything by default). Returns
  // true when there's an active selection AND we can confirm zero carriers
  // in it — caller uses that to dim the glyph.
  function _hasNoCarrierInSelection(svId) {
    const sel = _state.selectedSamples;
    const supp = _state.supportLayer;
    if (!sel || !supp) return false;        // can't determine → don't dim
    const colIdx = (supp.sv_ids || []).indexOf(svId);
    if (colIdx < 0) return false;            // SV not in support matrix
    const samples = supp.samples || [];
    for (let r = 0; r < samples.length; r++) {
      if (!sel.has(samples[r])) continue;
      const ch = _supportCell(supp, r, colIdx);
      if (ch === '1' || ch === '2') return false;   // at least one carrier
    }
    return true;                              // confirmed: 0 carriers in sel
  }

  function _nearestSV(bp) {
    const layer = _state.layer;
    if (!layer || !Array.isArray(layer.sv_calls) || bp == null) return null;
    const visible = _applyFilters(layer.sv_calls, _state.filters);
    if (visible.length === 0) return null;
    let best = null, bestD = Infinity;
    for (const c of visible) {
      const d = Math.abs((c.position_bp || 0) - bp);
      if (d < bestD) { bestD = d; best = c; }
    }
    return best ? { sv: best, dist: bestD } : null;
  }

  // Dispatch one track. Returns SVG markup (not the wrapper div).
  function _renderLocusTrack(t, layer, lo, hi, drawW, h, x) {
    if (t.id === 'zones')    return _drawTrackZones(layer, lo, hi, drawW, h, x);
    if (t.id === 'axis')     return _drawTrackAxis(layer, lo, hi, drawW, h, x);
    if (t.id === 'sv_calls') return _drawTrackSVCalls(layer, lo, hi, drawW, h, x);
    if (t.kind === 'state')  return _drawTrackPlaceholder(t, drawW, h);
    return _drawTrackPlaceholder(t, drawW, h);
  }

  // -- Track 1: Zones (5 coloured segments) ----------------------------------
  function _drawTrackZones(layer, lo, hi, drawW, h, x) {
    const z = layer.zone_definitions_bp || {};
    const segs = [
      { id: 'left_flank',     label: 'Left flank',     range: z.left_flank,     fill: '#5a6472' },
      { id: 'left_boundary',  label: 'Left boundary',  range: z.left_boundary,  fill: '#4fa3ff' },
      { id: 'inversion_body', label: 'Inversion body', range: z.inversion_body, fill: '#f5a524' },
      { id: 'right_boundary', label: 'Right boundary', range: z.right_boundary, fill: '#bf616a' },
      { id: 'right_flank',    label: 'Right flank',    range: z.right_flank,    fill: '#5a6472' },
    ].filter(s => Array.isArray(s.range));

    const rects = segs.map(s => {
      const x0 = x(s.range[0]);
      const x1 = x(s.range[1]);
      const w  = Math.max(0.5, x1 - x0);
      // Bar: top half is the colour band, label sits above on the boundary
      // segments (left_boundary / right_boundary / inversion_body) only.
      const labelY = 9;
      const showLabel = (s.id === 'left_boundary' || s.id === 'right_boundary' ||
                         s.id === 'inversion_body');
      const labelTxt = showLabel ? `<text x="${(x0 + w/2).toFixed(1)}" y="${labelY}"
                                          text-anchor="middle" font-size="10"
                                          fill="${s.fill}" font-weight="600">${_esc(s.label)}</text>` : '';
      // Bar at h/2 → h
      return `
        ${labelTxt}
        <rect x="${x0.toFixed(1)}" y="${(h - 8).toFixed(1)}"
              width="${w.toFixed(1)}" height="6"
              fill="${s.fill}" opacity="0.9" rx="1.5" ry="1.5"/>`;
    }).join('');

    return `<svg width="${drawW}" height="${h}" viewBox="0 0 ${drawW} ${h}"
                 preserveAspectRatio="none">${rects}</svg>`;
  }

  // -- Track 2: Coordinate axis ---------------------------------------------
  function _drawTrackAxis(layer, lo, hi, drawW, h, x) {
    // Tick every 500 kb, snap to nearest 500 kb above lo
    const STEP = 500000;
    const t0 = Math.ceil(lo / STEP) * STEP;
    const ticks = [];
    for (let bp = t0; bp <= hi; bp += STEP) ticks.push(bp);

    const baseY = h - 4;
    const tickH = 5;

    const axisLine = `<line x1="0" x2="${drawW}" y1="${baseY}" y2="${baseY}"
                            stroke="var(--ink-dim, #8a94a3)" stroke-width="0.6"/>`;
    const tickMarks = ticks.map(bp => {
      const xp = x(bp);
      const txt = (bp / 1e6).toFixed(2) + ' Mb';
      return `
        <line x1="${xp.toFixed(1)}" x2="${xp.toFixed(1)}"
              y1="${baseY}" y2="${baseY - tickH}"
              stroke="var(--ink-dim, #8a94a3)" stroke-width="0.6"/>
        <text x="${xp.toFixed(1)}" y="${baseY - tickH - 2}"
              text-anchor="middle" font-size="9.5"
              fill="var(--ink-dim, #8a94a3)">${txt}</text>`;
    }).join('');

    return `<svg width="${drawW}" height="${h}" viewBox="0 0 ${drawW} ${h}"
                 preserveAspectRatio="none">${axisLine}${tickMarks}</svg>`;
  }

  // -- Track 5: SV calls (one row per SV type, glyph per call) --------------
  // Per spec §3.4: ▼ BND red, ▼ INV violet, ▼ DEL orange, ▼ DUP teal, ▲ Other grey.
  // Click a glyph → highlights the table row (handshake via _state.highlightedSvId).
  function _drawTrackSVCalls(layer, lo, hi, drawW, h, x) {
    const calls = Array.isArray(layer.sv_calls) ? layer.sv_calls : [];
    if (calls.length === 0) {
      return `<svg width="${drawW}" height="${h}"><text x="${drawW/2}" y="${h/2}"
              text-anchor="middle" font-size="10"
              fill="var(--ink-dimmer, #5a6472)">no SV calls in layer</text></svg>`;
    }

    // Five lanes (BND, INV, DEL, DUP, Other), each ~h/5 tall
    const lanes = ['BND', 'INV', 'DEL', 'DUP', 'Other'];
    const laneH = h / lanes.length;

    // Apply current filters here so the locus view stays in sync with the
    // table. This is essential per acceptance §9: the locus reflects the
    // active filter set.
    const visible = _applyFilters(calls, _state.filters);
    const visibleSet = new Set(visible.map(c => c.sv_id));

    const elements = calls.map(c => {
      const laneIdx = lanes.indexOf(c.sv_type);
      const li = laneIdx >= 0 ? laneIdx : (lanes.length - 1); // unknown → Other lane
      const yMid = li * laneH + laneH / 2;
      const xp = x(c.position_bp);
      // Cull off-domain calls
      if (xp < -8 || xp > drawW + 8) return '';
      const style = SV_TYPE_STYLE[c.sv_type] || SV_TYPE_STYLE.Other;
      const filteredOut = !visibleSet.has(c.sv_id);
      // Step 5+6: when a sample selection is active AND we have a support
      // layer, SVs whose carriers don't overlap selected_samples are dimmed.
      // (Without supportLayer we can't compute this, so glyphs stay 1.0.)
      const noCarrierInSel = _hasNoCarrierInSelection(c.sv_id);
      let opacity = 1.0;
      if (filteredOut)            opacity = 0.18;
      else if (noCarrierInSel)    opacity = 0.35;
      const isHi = (_state.highlightedSvId === c.sv_id);
      const stroke = isHi ? 'var(--accent-2, #4fa3ff)' : 'none';
      const strokeW = isHi ? 1.5 : 0;
      // FDR halo: small ring underneath the glyph for FDR-coloured rim
      const fdr = (c.fisher && c.fisher.fdr_bh);
      const fdrCol = _fdrColour(fdr);

      // ▼ → polygon downward triangle; ▲ → upward
      const tri = (style.glyph === '▲')
        ? `${xp - 4},${yMid + 3.5} ${xp + 4},${yMid + 3.5} ${xp},${yMid - 3.5}`
        : `${xp - 4},${yMid - 3.5} ${xp + 4},${yMid - 3.5} ${xp},${yMid + 3.5}`;

      // Title gives a hover tooltip
      const tip = `${c.sv_id} · ${c.sv_type} · ${(c.position_bp / 1e6).toFixed(3)} Mb` +
                  ` · FDR ${fdr != null ? _fmtSci(fdr) : '—'}` +
                  (c.distance_to_edge_bp != null
                    ? ` · ${c.distance_to_edge_bp >= 0 ? '+' : ''}${c.distance_to_edge_bp.toLocaleString()} bp from edge`
                    : '');

      return `
        <g class="sv-locus-glyph" data-sv-glyph="${_esc(c.sv_id)}"
           style="cursor:pointer;" opacity="${opacity}">
          <title>${_esc(tip)}</title>
          <circle cx="${xp.toFixed(1)}" cy="${yMid.toFixed(1)}" r="6.5"
                  fill="${fdrCol}" opacity="0.18"/>
          <polygon points="${tri}" fill="${style.color}"
                   stroke="${stroke}" stroke-width="${strokeW}"/>
        </g>`;
    }).join('');

    // Lane separators (faint horizontal lines)
    const sep = lanes.slice(1).map((_, i) => {
      const y = (i + 1) * laneH;
      return `<line x1="0" x2="${drawW}" y1="${y}" y2="${y}"
                    stroke="var(--rule, #2a3242)" stroke-width="0.3"/>`;
    }).join('');

    // Lane labels at far left (very small)
    const laneLabels = lanes.map((t, i) => {
      const y = i * laneH + laneH / 2 + 3;
      return `<text x="2" y="${y}" font-size="8.5" fill="var(--ink-dimmer, #5a6472)">${t}</text>`;
    }).join('');

    return `<svg width="${drawW}" height="${h}" viewBox="0 0 ${drawW} ${h}"
                 preserveAspectRatio="none">
              ${sep}${laneLabels}${elements}
            </svg>`;
  }

  // -- Tracks 3/4/6/7: state-bound placeholders ------------------------------
  // Render an "empty / data not loaded" stripe with the track label and a
  // hint about what state key would populate it. Honest about what's
  // missing instead of pretending the track is empty.
  function _drawTrackPlaceholder(t, drawW, h) {
    const hint = ({
      dosage:  'Drop the per-window dosage matrix to render',
      pca:     'Drop the per-window local PCA to render (PC1/PC2/PC3 lines)',
      genes:   'Load a GFF gene track to render strand arrows',
      repeats: 'Load repeat-density per-window data to render',
    })[t.id] || 'Data not loaded';
    return `<svg width="${drawW}" height="${h}" viewBox="0 0 ${drawW} ${h}"
                 preserveAspectRatio="none">
              <rect x="0" y="0" width="${drawW}" height="${h}"
                    fill="var(--panel-2, #1c2330)" opacity="0.4"/>
              <text x="${drawW/2}" y="${h/2 + 3}" text-anchor="middle"
                    font-size="9.5" fill="var(--ink-dimmer, #5a6472)"
                    font-style="italic">${_esc(hint)}</text>
            </svg>`;
  }

  // Wire glyph click → highlight the corresponding table row.
  function _wireLocusInteractions(host) {
    if (!host) return;
    host.querySelectorAll('[data-sv-glyph]').forEach(g => {
      g.addEventListener('click', (ev) => {
        ev.stopPropagation();
        const id = g.dataset.svGlyph;
        _state.highlightedSvId = id;
        // Re-render locus (so the click gets the active rim) AND table
        // (so the corresponding row gets the inset accent border + scrolls
        // into view).
        _renderLocusStrip();
        _renderTable();
        _scrollTableRowIntoView(id);
      });
    });

    // Hit-rect: hover updates cursor; click sets cursor + highlights nearest.
    // The rect is the only overlay element with pointer-events:auto, so glyph
    // clicks (in tracks below the overlay) still receive their handlers via
    // SVG hit-testing (the rect is transparent and lives in the overlay
    // SVG, but pointer-events:none is set on every other overlay element so
    // they don't shadow glyph clicks).
    const hit = host.querySelector('[data-sv-locus-hit]');
    if (hit) {
      hit.addEventListener('mousemove', _onLocusMove);
      hit.addEventListener('mouseleave', _onLocusLeave);
      hit.addEventListener('mousedown', _onLocusMouseDown);
      hit.addEventListener('mouseup',   _onLocusMouseUp);
      hit.addEventListener('click',     _onLocusClick);
    }

    // Clear-selection pill in the readout (step 4.5)
    const clearBtn = host.querySelector('[data-sv-action="clear-selection"]');
    if (clearBtn) {
      clearBtn.addEventListener('click', (ev) => {
        ev.stopPropagation();
        _state.selection = null;
        _state.pageIndex = 0;
        _renderLocusStrip();
        _renderTable();
      });
    }
    // Clear UpSet sample selection pill (step 5)
    const clearSamplesBtn = host.querySelector('[data-sv-action="clear-samples"]');
    if (clearSamplesBtn) {
      clearSamplesBtn.addEventListener('click', (ev) => {
        ev.stopPropagation();
        _clearSampleSelection();
      });
    }
  }

  // Map a mouse event on the hit-rect to a bp value using the cached x-inverse.
  function _eventToBp(ev) {
    const w = _state._lastWindow;
    if (!w) return null;
    // Use bounding rect so we get the px relative to the SVG, not the page.
    const rect = ev.target.getBoundingClientRect();
    const px = ev.clientX - rect.left;
    // Clamp to window so a drag that exits the rect doesn't yield NaN bp
    return Math.max(w.lo, Math.min(w.hi, w.xInv(px)));
  }

  function _onLocusMove(ev) {
    const bp = _eventToBp(ev);
    if (bp == null) return;
    _state.cursorBp = bp;
    // If we're in the middle of a select drag, extend the rectangle.
    if (_state.selectMode === 'select' && _state._selectDrag) {
      _state._selectDrag.currentBp = bp;
    }
    // Cheap re-render of just the readout + overlay, not the whole strip.
    // For now, full re-render is fine — bounded SVG, no perf hit at 1000 SVs.
    _renderLocusStrip();
  }

  function _onLocusLeave() {
    // Hover leave — keep cursor sticky (matches page-11 behaviour where
    // tracks remember the last cursor on mouseleave). Esc clears it.
    // If a drag is in flight, treat mouseleave like mouseup so the user
    // doesn't end up with a stuck drag-state when they release outside.
    if (_state._selectDrag) {
      _commitSelectDrag();
    }
  }

  // Step 4.5 — mousedown starts a drag in select mode. In zoom mode it's
  // a no-op; the click handler does the cursor work.
  function _onLocusMouseDown(ev) {
    if (_state.selectMode !== 'select') return;
    const bp = _eventToBp(ev);
    if (bp == null) return;
    ev.preventDefault();
    _state._selectDrag = { startBp: bp, currentBp: bp };
    // Show the cursor at the drag-start position
    _state.cursorBp = bp;
    _renderLocusStrip();
  }

  // mouseup: in select mode, commit the drag as a sticky selection if the
  // drag spans more than a trivial number of bp; otherwise treat it as a
  // click (clear selection, place cursor). In zoom mode, mouseup is a
  // no-op — the click handler runs separately.
  function _onLocusMouseUp(ev) {
    if (_state.selectMode !== 'select') return;
    if (!_state._selectDrag) return;
    ev.preventDefault();
    _commitSelectDrag();
  }

  // Convert the in-progress drag into a sticky selection. If the drag was
  // shorter than ~1 kb (effectively a click in select mode), clear the
  // selection and treat it as a cursor placement.
  function _commitSelectDrag() {
    const d = _state._selectDrag;
    _state._selectDrag = null;
    if (!d) return;
    const lo = Math.min(d.startBp, d.currentBp);
    const hi = Math.max(d.startBp, d.currentBp);
    const span = hi - lo;
    if (span < 1000) {
      // Effectively a click — clear any existing selection, place cursor
      _state.selection = null;
      _state.cursorBp = (lo + hi) / 2;
    } else {
      _state.selection = { startBp: lo, endBp: hi };
      // Reset to first page so the user sees the top of the new filtered set
      _state.pageIndex = 0;
    }
    _renderLocusStrip();
    _renderTable();
  }

  function _onLocusClick(ev) {
    // In select mode the click event also fires after mouseup; suppress it
    // so we don't double-handle (zoom-mode cursor placement on top of the
    // drag commit).
    if (_state.selectMode === 'select') return;
    const bp = _eventToBp(ev);
    if (bp == null) return;
    _state.cursorBp = bp;
    // Promote to nearest SV: highlight + scroll table.
    const near = _nearestSV(bp);
    if (near) {
      _state.highlightedSvId = near.sv.sv_id;
      _renderLocusStrip();
      _renderTable();
      _scrollTableRowIntoView(near.sv.sv_id);
    } else {
      _renderLocusStrip();
    }
  }

  // After the table re-renders, scroll the highlighted row into view.
  function _scrollTableRowIntoView(svId) {
    if (!svId || !_state.rootEl) return;
    const tr = _state.rootEl.querySelector('[data-sv-row="' + svId + '"]');
    if (tr && typeof tr.scrollIntoView === 'function') {
      try {
        tr.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      } catch (_) { /* old browsers without smooth */ }
    } else if (svId) {
      // Row not in current page — flip to its page. Find rank in
      // filtered+sorted set; compute pageIndex; re-render.
      const layer = _state.layer;
      if (!layer) return;
      const filtered = _applyFilters(layer.sv_calls, _state.filters);
      const sorted   = _sortRows(filtered, _state.sortColumnId, _state.sortDirection);
      const idx = sorted.findIndex(c => c.sv_id === svId);
      if (idx >= 0) {
        _state.pageIndex = Math.floor(idx / _state.pageSize);
        _renderTable();
        // Re-scroll after the new render
        requestAnimationFrame(() => {
          const tr2 = _state.rootEl.querySelector('[data-sv-row="' + svId + '"]');
          if (tr2 && typeof tr2.scrollIntoView === 'function') {
            try { tr2.scrollIntoView({ behavior: 'smooth', block: 'nearest' }); }
            catch (_) {}
          }
        });
      }
    }
  }

  // ===========================================================================
  // Page-scoped hotkeys (step 3.5) — mirrors page 11's _bndAttachHotkeys.
  // Only fires when the SV evidence page is .active AND no input/select has
  // focus AND no Ctrl/Meta/Alt is held. Keys:
  //   ←/→         move cursor by 0.5% of window span
  //   Shift+←/→   jump to next/prev SV glyph (filtered set)
  //   ↑/↓         cycle SV-type filter (All → BND → INV → DEL → DUP → Other)
  //   E           drop a green left-pin marker at the cursor
  //   F           drop an orange right-pin marker at the cursor
  //   Enter       highlight the nearest SV (same as click) + scroll table
  //   Esc         clear cursor + both markers + highlighted SV
  //   1/2/3       switch view preset (Default / Left close-up / Right close-up)
  // ===========================================================================

  const SV_TYPE_CYCLE = ['All', 'BND', 'INV', 'DEL', 'DUP', 'Other'];

  function _onLocusHotkey(ev) {
    // Only fire when the SV evidence page is the active tab
    const root = _state.rootEl
      ? _state.rootEl.closest('.page')
      : document.getElementById('page_sv_evidence');
    if (!root || !root.classList.contains('active')) return;

    // Don't hijack typing in form controls
    const ae = (typeof document !== 'undefined') ? document.activeElement : null;
    if (ae && /^(INPUT|TEXTAREA|SELECT|BUTTON)$/.test(ae.tagName)) {
      // Allow Esc to still fire even when filter inputs have focus —
      // Esc is universally "clear" and shouldn't be eaten by a select.
      if (ev.key !== 'Escape') return;
    }

    // Skip browser-shortcut combos
    if (ev.ctrlKey || ev.metaKey || ev.altKey) return;

    // Need a layer + window for any of these to be meaningful
    const layer = _state.layer;
    const win = _state._lastWindow;
    if (!layer || !win) return;

    const key = ev.key;

    if (key === 'ArrowLeft' || key === 'ArrowRight') {
      ev.preventDefault();
      const dir = key === 'ArrowRight' ? 1 : -1;
      if (ev.shiftKey) {
        _jumpCursorToNextSV(dir);
      } else {
        _stepCursor(dir);
      }
      return;
    }

    if (key === 'ArrowUp' || key === 'ArrowDown') {
      ev.preventDefault();
      if (ev.shiftKey) {
        // Shift+↑/↓ — cycle the SV-type filter
        _cycleSvTypeFilter(key === 'ArrowDown' ? 1 : -1);
      } else {
        // Plain ↑/↓ — step through table rows (Quentin's step-4 request).
        // Updates highlightedSvId, scrolls table, lights up locus glyph.
        _stepHighlightedRow(key === 'ArrowDown' ? 1 : -1);
      }
      return;
    }

    if (key === 'Enter') {
      ev.preventDefault();
      const near = _nearestSV(_state.cursorBp);
      if (near) {
        _state.highlightedSvId = near.sv.sv_id;
        _renderLocusStrip();
        _renderTable();
        _scrollTableRowIntoView(near.sv.sv_id);
      }
      return;
    }

    if (key === 'Escape') {
      ev.preventDefault();
      _state.cursorBp = null;
      _state.markerLeftBp = null;
      _state.markerRightBp = null;
      _state.highlightedSvId = null;
      _state.selection = null;       // step 4.5: also clear region selection
      _state._selectDrag = null;
      _state.selectedSamples = null; // step 5: also clear UpSet sample selection
      _state.activeCombinationIndex = null;
      _state.gtCountsView = null;    // step 6: also clear the recount view
      _state.heatmapHoverCell = null;
      // Close large heatmap if open
      if (_state.heatmapView === 'large') {
        _state.heatmapView = 'compact';
        const ov = _state.rootEl && _state.rootEl.querySelector('[data-sv-hm-overlay]');
        if (ov && ov.parentNode) ov.parentNode.removeChild(ov);
      }
      _state.pageIndex = 0;
      _renderLocusStrip();
      _renderTable();
      _renderRightRail();
      return;
    }

    const k = key.toLowerCase();
    if (k === 'e') {
      ev.preventDefault();
      if (_state.cursorBp != null) {
        _state.markerLeftBp = _state.cursorBp;
        _renderLocusStrip();
      }
      return;
    }
    if (k === 'f') {
      ev.preventDefault();
      if (_state.cursorBp != null) {
        _state.markerRightBp = _state.cursorBp;
        _renderLocusStrip();
      }
      return;
    }

    // View-preset switch via 1/2/3 (mirroring the page-11 hotkey-on-toolbar
    // muscle memory; the toolbar dropdown still works too).
    if (key === '1') { ev.preventDefault(); _setViewPreset('default');     return; }
    if (key === '2') { ev.preventDefault(); _setViewPreset('left_close');  return; }
    if (key === '3') { ev.preventDefault(); _setViewPreset('right_close'); return; }
  }

  function _stepCursor(dir) {
    const w = _state._lastWindow;
    if (!w) return;
    const span = w.hi - w.lo;
    const step = span * 0.005;       // 0.5% of window per arrow press
    const cur = (_state.cursorBp != null) ? _state.cursorBp : (w.lo + span / 2);
    let next = cur + dir * step;
    if (next < w.lo) next = w.lo;
    if (next > w.hi) next = w.hi;
    _state.cursorBp = next;
    _renderLocusStrip();
  }

  function _jumpCursorToNextSV(dir) {
    const layer = _state.layer;
    const w = _state._lastWindow;
    if (!layer || !w) return;
    // Only consider SVs within the current visible window AND active filters
    const visible = _applyFilters(layer.sv_calls, _state.filters)
      .filter(c => c.position_bp >= w.lo && c.position_bp <= w.hi)
      .sort((a, b) => (a.position_bp || 0) - (b.position_bp || 0));
    if (visible.length === 0) return;
    const cur = (_state.cursorBp != null) ? _state.cursorBp : (w.lo + (w.hi - w.lo) / 2);
    let target = null;
    if (dir > 0) {
      for (const c of visible) if (c.position_bp > cur + 1) { target = c; break; }
      if (!target) target = visible[0];                   // wrap
    } else {
      for (let i = visible.length - 1; i >= 0; i--) {
        if (visible[i].position_bp < cur - 1) { target = visible[i]; break; }
      }
      if (!target) target = visible[visible.length - 1];  // wrap
    }
    if (target) {
      _state.cursorBp = target.position_bp;
      _state.highlightedSvId = target.sv_id;
      _renderLocusStrip();
      _renderTable();
      _scrollTableRowIntoView(target.sv_id);
    }
  }

  // Step the highlighted row through the visible table rows in current
  // sort order. Wraps. If no row is highlighted yet, picks the first
  // visible row. Re-renders the table + locus, and scrolls the new row
  // into view. Crosses page boundaries (paginates) when necessary.
  function _stepHighlightedRow(dir) {
    const layer = _state.layer;
    if (!layer || !Array.isArray(layer.sv_calls) || layer.sv_calls.length === 0) return;
    const filtered = _applyFilters(layer.sv_calls, _state.filters);
    if (filtered.length === 0) return;
    const sorted = _sortRows(filtered, _state.sortColumnId, _state.sortDirection);
    let idx = sorted.findIndex(c => c.sv_id === _state.highlightedSvId);
    if (idx < 0) {
      idx = (dir > 0) ? 0 : (sorted.length - 1);
    } else {
      idx = (idx + dir + sorted.length) % sorted.length;
    }
    const next = sorted[idx];
    _state.highlightedSvId = next.sv_id;
    // Move cursor to the new row's position so the locus tracks it
    if (next.position_bp != null) _state.cursorBp = next.position_bp;
    // Flip to the new row's page if it's outside the current page
    const onPage = Math.floor(idx / _state.pageSize);
    if (onPage !== _state.pageIndex) _state.pageIndex = onPage;
    _renderTable();
    _renderLocusStrip();
    _scrollTableRowIntoView(next.sv_id);
  }

  function _cycleSvTypeFilter(dir) {
    const cur = SV_TYPE_CYCLE.indexOf(_state.filters.sv_type);
    const idx = ((cur >= 0 ? cur : 0) + dir + SV_TYPE_CYCLE.length) % SV_TYPE_CYCLE.length;
    _state.filters = { ..._state.filters, sv_type: SV_TYPE_CYCLE[idx] };
    _writeFilters(_state.activeCandidateId, _state.filters);
    _state.pageIndex = 0;
    // Re-render filter form so the SELECT shows the new value
    const filEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-filters]');
    if (filEl) _renderFiltersBlock(filEl);
    _renderLocusStrip();
    _renderTable();
  }

  function _setViewPreset(preset) {
    if (preset !== 'default' && preset !== 'left_close' && preset !== 'right_close') return;
    _state.viewPreset = preset;
    _state._customWindow = null;
    // If the cursor is now outside the new window, clip it to the new bounds
    // (rather than dropping it — clipping preserves user intent).
    const layer = _state.layer;
    if (layer) {
      const w = _windowForPreset(layer, preset);
      if (_state.cursorBp != null) {
        if (_state.cursorBp < w.lo) _state.cursorBp = w.lo;
        if (_state.cursorBp > w.hi) _state.cursorBp = w.hi;
      }
    }
    // Reflect in the toolbar dropdown
    const presetEl = _state.rootEl && _state.rootEl.querySelector('[data-sv-preset]');
    if (presetEl) presetEl.value = preset;
    _renderLocusStrip();
  }

  function _attachLocusHotkeys() {
    if (_state.hotkeysAttached) return;
    if (typeof document === 'undefined') return;
    document.addEventListener('keydown', _onLocusHotkey);
    _state.hotkeysAttached = true;
  }

  function _detachLocusHotkeys() {
    if (!_state.hotkeysAttached) return;
    if (typeof document === 'undefined') return;
    document.removeEventListener('keydown', _onLocusHotkey);
    _state.hotkeysAttached = false;
  }


  // Open a hidden file input and read the chosen JSON. Mirrors the regime
  // import pattern at line ~12618 of Inversion_atlas.html.
  function _triggerLoadLayerFile() {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = 'application/json,.json';
    input.style.display = 'none';
    input.addEventListener('change', () => {
      const file = input.files && input.files[0];
      if (input.parentNode) input.parentNode.removeChild(input);
      if (!file) return;
      _readAndIngestFile(file);
    });
    document.body.appendChild(input);
    input.click();
  }

  // Drag-drop handler. Wired to the page's outer .sv-evidence-page element
  // in _renderShell so dropping anywhere on the page works.
  function _handlePageDrop(ev) {
    ev.preventDefault();
    ev.stopPropagation();
    if (_state.rootEl) {
      const wrap = _state.rootEl.querySelector('.sv-evidence-page');
      if (wrap) wrap.classList.remove('sv-drop-active');
    }
    const dt = ev.dataTransfer;
    if (!dt) return;
    // Try the items API (lets us detect folder drops via webkitGetAsEntry).
    // Fall back to dt.files for environments without it.
    if (dt.items && dt.items.length > 0 && typeof dt.items[0].webkitGetAsEntry === 'function') {
      const entries = [];
      for (let i = 0; i < dt.items.length; i++) {
        const e = dt.items[i].webkitGetAsEntry();
        if (e) entries.push(e);
      }
      _ingestEntries(entries);
      return;
    }
    if (dt.files && dt.files.length > 0) {
      // No folder support — single-file drop. If multiple files were dropped,
      // ingest them all.
      const promises = [];
      for (let i = 0; i < dt.files.length; i++) {
        promises.push(_readFileAsText(dt.files[i]).then(_ingestJsonText));
      }
      Promise.all(promises).then(() => refresh());
    }
  }

  // Walk a list of FileSystemEntry objects (single-files or directories).
  // Collect all *.json files and ingest them in dispatch order. Folder
  // drops are the canonical local-first path: drop a candidate folder,
  // populate every layer the atlas knows.
  function _ingestEntries(entries) {
    const filePromises = [];
    function visit(entry) {
      if (entry.isFile) {
        if (!/\.json$/i.test(entry.name)) return;
        filePromises.push(new Promise((resolve) => {
          entry.file(f => resolve(_readFileAsText(f).then(_ingestJsonText)));
        }));
      } else if (entry.isDirectory) {
        const reader = entry.createReader();
        // readEntries is paginated — drain it
        function drain() {
          reader.readEntries((ents) => {
            if (!ents.length) return;
            ents.forEach(visit);
            drain();
          });
        }
        drain();
      }
    }
    entries.forEach(visit);
    // Wait a tick for folder traversal to schedule, then await all reads.
    setTimeout(() => {
      Promise.all(filePromises).then(() => refresh());
    }, 50);
  }

  function _readFileAsText(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result);
      reader.onerror = () => reject(reader.error || new Error('read failed'));
      reader.readAsText(file);
    });
  }

  // Dispatch by format_version. Sets the matching _state slot. Idempotent;
  // a second drop overwrites the same slot.
  function _ingestJsonText(text) {
    let parsed;
    try {
      parsed = JSON.parse(text);
    } catch (e) {
      _state.layerError = 'JSON parse failed: ' + (e && e.message ? e.message : e);
      return;
    }
    if (!parsed || typeof parsed !== 'object') return;
    const fv = parsed.format_version;
    if (fv === 'sv_genotype_counts_v1') {
      if (parsed.candidate_id) {
        _state.activeCandidateId = parsed.candidate_id;
        _state.filters = _readFilters(_state.activeCandidateId);
        _state.rowAnnotations = _readAnnotations(_state.activeCandidateId);
      }
      _state.layer = parsed;
      _state.layerError = null;
      _state.layerLoading = false;
      _state.pageIndex = 0;
      return;
    }
    if (fv === 'sv_evidence_combinations_v1') {
      const v = _validateCombinationsLayer(parsed);
      if (!v.ok) {
        _state.layerError = 'invalid combinations layer: ' + v.reason;
        return;
      }
      _state.combinationsLayer = parsed;
      // Don't auto-clear an in-flight selection — the user might be
      // mid-investigation. They can press Esc to clear.
      return;
    }
    if (fv === 'sv_support_by_sample_v1') {
      const v = _validateSupportLayer(parsed);
      if (!v.ok) {
        _state.layerError = 'invalid support layer: ' + v.reason;
        return;
      }
      _state.supportLayer = parsed;
      return;
    }
    _state.layerError = 'unknown format_version: "' + fv + '"';
  }

  // Backwards-compatible single-file ingest (used by the old API and tests).
  function _readAndIngestFile(file) {
    if (!file) return;
    _readFileAsText(file).then(text => {
      _ingestJsonText(text);
      refresh();
    }).catch(e => {
      _state.layerError = 'file read failed: ' + (e && e.message ? e.message : e);
      _state.layer = null;
      _state.layerLoading = false;
      refresh();
    });
  }

  // ===========================================================================
  // Table renderer (spec §3.5)
  // ===========================================================================

  // Render the main SV table into the .sv-table-wrap container. Replaces
  // the step-1 "Ready (skeleton)" hint when a layer is present. Idempotent:
  // safe to call after every filter / sort / page-index change.
  function _renderTable() {
    if (!_state.rootEl) return;
    const wrap  = _state.rootEl.querySelector('[data-sv-table-wrap]');
    const empty = _state.rootEl.querySelector('[data-sv-empty-main]');
    if (!wrap) return;

    const layer = _state.layer;
    const groups = _resolveKaryotypeGroups(_state.activeCandidateId);

    // No layer → table area stays empty; the empty-state card in
    // .sv-locus-strip carries the user-facing message (see _renderMain).
    if (!layer) {
      wrap.innerHTML = '';
      return;
    }

    // Hide the main empty-state card when the table is going to render.
    if (empty) empty.style.display = 'none';

    const visible = _computeVisibleRows(
      layer, _state.filters,
      _state.sortColumnId, _state.sortDirection,
      _state.pageSize, _state.pageIndex
    );

    // Rebuild table HTML. Table is lightweight (<= page_size rows × ~22
    // columns) so a full innerHTML rebuild is fine — no virtual scroll
    // needed.
    wrap.innerHTML = `
      <div class="sv-tabletop">
        <div class="sv-tablebar">
          <span class="sv-tablebar-label">SVs near boundaries (±500 kb)</span>
          <div class="sv-tablebar-tabs">
            <button class="sv-tab active" data-sv-table-tab="summary">Summary table</button>
            <button class="sv-tab" data-sv-table-tab="genotype" disabled
                    title="Step 2: separate genotype-counts tab uses the same dataset; defer to step 3+">Genotype counts by group</button>
            <button class="sv-tab" data-sv-table-tab="upset" disabled
                    title="Step 5: UpSet panel">UpSet (combinations)</button>
            <button class="sv-tab" data-sv-table-tab="heatmap" disabled
                    title="Step 6: sample × SV heatmap">Sample × SV heatmap</button>
          </div>
          <div class="sv-tablebar-spacer"></div>
          <button class="sv-btn-secondary" data-sv-action="export-tsv"
                  title="Export the currently filtered + sorted rows as TSV">⤓ Export table (TSV)</button>
        </div>
      </div>
      <div class="sv-table-scroll">
        ${_renderTableHTML(visible, groups)}
      </div>
      <div class="sv-table-footer">
        <div class="sv-pager">
          <span>Rows per page</span>
          <select data-sv-action="page-size">
            ${PAGE_SIZES.map(s =>
              `<option value="${s}"${s === _state.pageSize ? ' selected' : ''}>${s}</option>`).join('')}
          </select>
          <span class="sv-pager-range">
            ${visible.startRank}–${visible.endRank} of ${visible.total}
          </span>
          ${_renderPagerButtons(visible)}
        </div>
        <div class="sv-table-legend-line">
          AA = 0/0 (ref/ref) &nbsp;&nbsp; AB = 0/1 (het) &nbsp;&nbsp;
          BB = 1/1 (alt/alt) &nbsp;&nbsp; OR = odds ratio
        </div>
      </div>
    `;

    _wireTableInteractions(wrap);
  }

  // Build the inner <table> HTML. Two-row header: top row spans groups,
  // second row carries individual columns. Body row per visible SV call.
  function _renderTableHTML(visible, groups) {
    if (!visible || visible.total === 0) {
      return `
        <div class="sv-table-empty">
          <strong>No SV calls match current filters.</strong><br>
          <span class="sv-empty-detail">
            ${_state.layer && _state.layer.sv_calls.length
              ? _state.layer.sv_calls.length + ' total in layer; relax filters to see more.'
              : 'No SV calls in layer.'}
          </span>
        </div>`;
    }

    // Top header — group spans. Patch in n=? from groups.counts.
    // When _state.gtCountsView is active (UpSet selection + supportLayer),
    // show the in-selection per-group sample counts instead so the user
    // sees the right denominator at a glance. A "(in sel)" tag flags the
    // mode. The full-cohort label is recovered the moment Esc is pressed.
    const view = _state.gtCountsView;
    const gnSel = view && view._groupNs;
    const _gn = (g, full) => {
      if (gnSel && typeof gnSel[g] === 'number') return `${gnSel[g]}/${full} sel`;
      return String(full);
    };
    const topHeaders = TABLE_GROUP_HEADERS.map(g => {
      let label = g.label;
      if (g.group === 'gc-h1h1') label = `H1/H1 (n=${_gn('H1/H1', groups.counts['H1/H1'] || 0)})`;
      if (g.group === 'gc-h1h2') label = `H1/H2 (n=${_gn('H1/H2', groups.counts['H1/H2'] || 0)})`;
      if (g.group === 'gc-h2h2') label = `H2/H2 (n=${_gn('H2/H2', groups.counts['H2/H2'] || 0)})`;
      // 3 group headers carry borders; rest are blank
      const isGroup = g.group.startsWith('gc-') || g.group === 'fisher';
      return `<th class="sv-th-group${isGroup ? ' sv-th-group-divider' : ''}"
                  colspan="${g.cols}">${_esc(label)}</th>`;
    }).join('');

    // Second row — individual column headers
    const colHeaders = TABLE_COLUMNS.map(col => {
      const isSorted = (_state.sortColumnId === col.id);
      const arrow = isSorted ? (_state.sortDirection === 'asc' ? ' ↑' : ' ↓') : '';
      return `<th class="sv-th sv-th-col sv-th-${col.align}"
                  data-sv-sort-col="${col.id}"
                  title="Click to sort by ${_esc(col.label)}">${_esc(col.label)}${arrow}</th>`;
    }).join('');

    // Body
    const bodyRows = visible.rows.map(call => {
      const cells = TABLE_COLUMNS.map(col => {
        const v = col.get(call);
        const txt = col.fmt(v, call);
        let cellStyle = '';
        let cellClass = 'sv-td sv-td-' + col.align;
        // FDR cell: background-coloured per spec §5.2
        if (col.id === 'fdr') {
          const c = _fdrColour(v);
          cellStyle = `background:${c};color:#0e1116;font-weight:600;`;
        }
        // Pattern cell: text-coloured per spec §2.3 / §5.3
        if (col.id === 'pattern') {
          const pl = PATTERN_LABELS[call.pattern_label];
          if (pl) cellStyle = `color:${pl.color};font-weight:600;`;
        }
        // SV type cell: dot + label
        if (col.id === 'sv_type') {
          const st = SV_TYPE_STYLE[call.sv_type] || SV_TYPE_STYLE.Other;
          return `<td class="${cellClass}"><span class="sv-typedot" style="color:${st.color}">${st.glyph}</span> ${_esc(txt)}</td>`;
        }
        // Zone cell: coloured per zone
        if (col.id === 'zone') {
          const z = call.zone;
          const zc = ({
            left_boundary: '#4fa3ff', right_boundary: '#bf616a',
            inversion_body: '#f5a524', left_flank: '#888', right_flank: '#888'
          })[z] || '#888';
          return `<td class="${cellClass}" style="color:${zc}">${_esc(txt)}</td>`;
        }
        // sv_id cell: annotation chip + sv_id button. Annotation cycles
        // through null → L (putative left breakpoint) → R (putative right)
        // → ★ (starred) → null on click. Stored per-candidate in localStorage.
        if (col.id === 'sv_id') {
          const anno = _state.rowAnnotations[call.sv_id] || null;
          const chip = _renderAnnoChip(anno);
          const bold = (_state.highlightedSvId === call.sv_id) ? ' sv-row-hi' : '';
          // Pattern-highlight class (legend click). Adds a soft tint.
          const phClass = (_state.highlightPattern && call.pattern_label === _state.highlightPattern)
            ? ' sv-row-pattern-hi' : '';
          return `<td class="${cellClass}${bold}${phClass}">` +
            `<button class="sv-anno-btn" data-sv-toggle-anno="${_esc(call.sv_id)}"
                     title="Click to cycle annotation (L → R → ★ → none)">${chip}</button>` +
            `<button class="sv-id-btn" data-sv-row-id="${_esc(call.sv_id)}"
                     title="Click to highlight this SV in the locus">${_esc(txt)}</button></td>`;
        }
        return `<td class="${cellClass}" style="${cellStyle}">${_esc(txt)}</td>`;
      }).join('');
      return `<tr class="sv-tr" data-sv-row="${_esc(call.sv_id)}">${cells}</tr>`;
    }).join('');

    return `
      <table class="sv-table">
        <thead>
          <tr class="sv-thead-top">${topHeaders}</tr>
          <tr class="sv-thead-cols">${colHeaders}</tr>
        </thead>
        <tbody>${bodyRows}</tbody>
      </table>`;
  }

  function _renderPagerButtons(visible) {
    const cur = visible.pageIndex;
    const last = visible.lastPage;
    // For very long lists, build a windowed page list with ellipses.
    const pages = [];
    if (last <= 6) {
      for (let i = 0; i <= last; i++) pages.push(i);
    } else {
      pages.push(0);
      if (cur > 2) pages.push('…');
      for (let i = Math.max(1, cur - 1); i <= Math.min(last - 1, cur + 1); i++) {
        pages.push(i);
      }
      if (cur < last - 2) pages.push('…');
      pages.push(last);
    }
    const prev = `<button class="sv-pager-btn" data-sv-page="prev" ${cur === 0 ? 'disabled' : ''}>‹</button>`;
    const next = `<button class="sv-pager-btn" data-sv-page="next" ${cur >= last ? 'disabled' : ''}>›</button>`;
    const items = pages.map(p => {
      if (p === '…') return '<span class="sv-pager-ellipsis">…</span>';
      const active = (p === cur) ? ' sv-pager-active' : '';
      return `<button class="sv-pager-btn${active}" data-sv-page="${p}">${p + 1}</button>`;
    }).join('');
    return prev + items + next;
  }

  function _wireTableInteractions(wrap) {
    if (!wrap) return;
    // Sort by column header click
    wrap.querySelectorAll('[data-sv-sort-col]').forEach(th => {
      th.addEventListener('click', () => {
        const colId = th.dataset.svSortCol;
        if (_state.sortColumnId === colId) {
          _state.sortDirection = (_state.sortDirection === 'asc') ? 'desc' : 'asc';
        } else {
          _state.sortColumnId = colId;
          _state.sortDirection = 'asc';
        }
        _state.pageIndex = 0;
        _renderTable();
      });
    });

    // Pagination
    wrap.querySelectorAll('[data-sv-page]').forEach(btn => {
      btn.addEventListener('click', () => {
        const v = btn.dataset.svPage;
        const last = Math.max(0, Math.ceil(_filteredCount() / _state.pageSize) - 1);
        if (v === 'prev') _state.pageIndex = Math.max(0, _state.pageIndex - 1);
        else if (v === 'next') _state.pageIndex = Math.min(last, _state.pageIndex + 1);
        else _state.pageIndex = parseInt(v, 10) || 0;
        _renderTable();
      });
    });
    const ps = wrap.querySelector('[data-sv-action="page-size"]');
    if (ps) ps.addEventListener('change', () => {
      _state.pageSize = parseInt(ps.value, 10) || DEFAULT_PAGE_SIZE;
      _state.pageIndex = 0;
      _renderTable();
    });

    // Annotation cycle (step 4) — replaces starring; click cycles
    // null → L → R → ★ → null, persisted per-candidate.
    wrap.querySelectorAll('[data-sv-toggle-anno]').forEach(btn => {
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        _cycleAnnotation(btn.dataset.svToggleAnno);
        _renderTable();
      });
    });

    // sv_id button click → highlight that SV in locus + table (without
    // changing annotation). Stops propagation so the row-click handler
    // below doesn't fire twice.
    wrap.querySelectorAll('[data-sv-row-id]').forEach(btn => {
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        _state.highlightedSvId = btn.dataset.svRowId;
        _renderTable();
        _renderLocusStrip();
      });
    });

    // Row click → highlight + sync locus glyph.
    wrap.querySelectorAll('[data-sv-row]').forEach(tr => {
      tr.addEventListener('click', () => {
        _state.highlightedSvId = tr.dataset.svRow;
        _renderTable();
        _renderLocusStrip();
      });
    });

    // TSV export
    const exp = wrap.querySelector('[data-sv-action="export-tsv"]');
    if (exp) exp.addEventListener('click', _downloadFilteredTSV);
  }

  function _filteredCount() {
    if (!_state.layer) return 0;
    return _applyFilters(_state.layer.sv_calls, _state.filters).length;
  }

  // ===========================================================================
  // TSV export
  // ===========================================================================

  // Build TSV text from the currently-filtered, currently-sorted full set
  // (NOT the current page). Header order matches TABLE_COLUMNS.
  function _buildFilteredTSV() {
    const layer = _state.layer;
    if (!layer) return 'sv_id\tSV type\tZone\tPosition (bp)\tFDR\n';
    const groups = _resolveKaryotypeGroups(_state.activeCandidateId);
    const visible = _computeVisibleRows(
      layer, _state.filters,
      _state.sortColumnId, _state.sortDirection,
      Number.MAX_SAFE_INTEGER, 0
    );
    const rows = visible.filteredAll || [];

    // Header — re-uses TABLE_COLUMNS labels but with group annotation in
    // genotype-count headers so the TSV is interpretable standalone.
    // Append the user's row-annotation column at the end so the export
    // captures the (L)/(R)/star marks researchers added during review.
    const header = TABLE_COLUMNS.map(col => {
      if (col.id.startsWith('h1h1_')) return 'H1/H1 (n=' + (groups.counts['H1/H1'] || 0) + ') ' + col.label;
      if (col.id.startsWith('h1h2_')) return 'H1/H2 (n=' + (groups.counts['H1/H2'] || 0) + ') ' + col.label;
      if (col.id.startsWith('h2h2_')) return 'H2/H2 (n=' + (groups.counts['H2/H2'] || 0) + ') ' + col.label;
      return col.label;
    });
    header.push('User annotation');
    const lines = [header.join('\t')];

    // Provenance comment lines (TSV readers like R's read.table skip lines
    // starting with #). Useful for the manuscript supplement.
    const annoCount = Object.keys(_state.rowAnnotations || {}).length;
    const meta = [
      '# atlas_sv_evidence TSV export',
      '# candidate_id: '   + (layer.candidate_id || ''),
      '# chrom: '          + (layer.chrom || ''),
      '# boundary_left: '  + (layer.boundary_left_bp != null ? layer.boundary_left_bp : ''),
      '# boundary_right: ' + (layer.boundary_right_bp != null ? layer.boundary_right_bp : ''),
      '# filters: '        + JSON.stringify(_state.filters),
      '# sort: '           + _state.sortColumnId + ' ' + _state.sortDirection,
      '# n_rows: '         + rows.length,
      '# n_annotations: '  + annoCount,
      '# exported: '       + new Date().toISOString(),
    ].join('\n');

    for (const call of rows) {
      const annoCols = TABLE_COLUMNS.map(col => col.tsv(call));
      annoCols.push(_state.rowAnnotations[call.sv_id] || '');
      lines.push(annoCols.join('\t'));
    }
    return meta + '\n' + lines.join('\n') + '\n';
  }

  function _downloadFilteredTSV() {
    const tsv = _buildFilteredTSV();
    const cid = _state.activeCandidateId || 'unknown';
    const ts  = new Date().toISOString().slice(0, 10);
    const filename = 'sv_evidence_' + cid + '_' + ts + '.tsv';
    try {
      const blob = new Blob([tsv], { type: 'text/tab-separated-values;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = filename;
      document.body.appendChild(a);
      a.click();
      setTimeout(() => {
        if (a.parentNode) a.parentNode.removeChild(a);
        URL.revokeObjectURL(url);
      }, 100);
    } catch (e) {
      console.warn('[svEvidence] TSV download failed:', e && e.message ? e.message : e);
    }
  }

  // Pull a chrom hint from the global candidate registry for the loading
  // empty-state, so the left rail isn't completely blank during the fetch.
  function _resolveCandidateChrom(cid) {
    const state = (typeof window !== 'undefined') ? window.state : null;
    if (!state) return null;
    const list = state.candidateList || (state.atlas_catalogue && state.atlas_catalogue.candidates) || [];
    for (const c of list) {
      if (!c) continue;
      if (c.id === cid || c.candidate_id === cid) return c.chrom || null;
    }
    return null;
  }

  // ===========================================================================
  // Public API
  // ===========================================================================

  function init(opts) {
    opts = opts || {};
    const root = opts.root || document.getElementById('sv_evidence_root') ||
                 document.getElementById('page_sv_evidence');
    if (!root) {
      // Page DOM not yet present — caller will retry on page activation.
      return false;
    }
    _state.rootEl = root;
    _renderShell();
    _renderLeftRail();
    _renderMain();
    _renderRightRail();
    _state.mounted = true;
    // Listen for scrubber karyotype lock changes.
    if (typeof window !== 'undefined' && window.addEventListener) {
      window.addEventListener('karyotype_locks_changed', () => {
        try { refresh(); } catch (_) {}
      });
    }
    return true;
  }

  function loadCandidate(cid) {
    if (!_state.rootEl) {
      // First-time mount — build the shell so empty-states render even
      // before any data arrives.
      const ok = init();
      if (!ok) return Promise.resolve(null);
    }
    if (!cid) {
      cid = _resolveActiveCandidateId();
    }
    _state.activeCandidateId = cid || null;
    _state.filters = _readFilters(_state.activeCandidateId);
    _state.rowAnnotations = _readAnnotations(_state.activeCandidateId);

    // No candidate → empty state.
    if (!_state.activeCandidateId) {
      _state.layer = null;
      _state.layerLoading = false;
      _state.layerError = null;
      _renderLeftRail();
      _renderMain();
      _renderRightRail();
      return Promise.resolve(null);
    }

    // Render loading state, then fetch.
    _state.layerLoading = true;
    _state.layerError = null;
    _state.layer = null;
    _renderLeftRail();
    _renderMain();
    _renderRightRail();

    return _loadLayer(_state.activeCandidateId).then(json => {
      _renderLeftRail();
      _renderMain();
      _renderRightRail();
      return json;
    });
  }

  function refresh() {
    if (!_state.rootEl) return;
    _renderLeftRail();
    _renderMain();
    _renderRightRail();
  }

  function setFilters(filters) {
    _state.filters = { ...DEFAULT_FILTERS, ...(filters || {}) };
    _writeFilters(_state.activeCandidateId, _state.filters);
    refresh();
  }

  function exportFilteredTSV() {
    return _buildFilteredTSV();
  }

  // ===========================================================================
  // CSS — injected once on first init()
  // ===========================================================================

  const STYLES = `
    /* Page container — fills the available .page area. */
    .sv-evidence-page {
      display: flex; flex-direction: column;
      height: 100%; min-height: 0;
      background: var(--bg, #0e1116); color: var(--ink, #e7edf3);
      font-family: var(--ui, system-ui, -apple-system, sans-serif);
      font-size: 12px;
    }

    /* Top toolbar */
    .sv-evidence-page .sv-toolbar {
      display: flex; align-items: center; gap: 12px;
      padding: 6px 12px;
      background: var(--panel, #151a22);
      border-bottom: 1px solid var(--rule, #2a3242);
      flex: 0 0 auto;
    }
    .sv-evidence-page .sv-toolbar-locus {
      font-weight: 600; color: var(--ink, #e7edf3);
    }
    .sv-evidence-page .sv-toolbar-spacer { flex: 1 1 auto; }
    .sv-evidence-page .sv-toolbar-presets,
    .sv-evidence-page .sv-toolbar-mode,
    .sv-evidence-page .sv-toolbar-zoom {
      display: inline-flex; align-items: center; gap: 4px;
      color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-toolbar-presets select,
    .sv-evidence-page .sv-toolbar-zoom button,
    .sv-evidence-page .sv-mode-btn {
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
      border: 1px solid var(--rule, #2a3242); border-radius: 3px;
      padding: 2px 8px; font-size: 11px; cursor: pointer;
    }
    .sv-evidence-page .sv-mode-btn.sv-mode-active {
      background: var(--accent-2, #4fa3ff); color: #0e1116;
      border-color: var(--accent-2, #4fa3ff); font-weight: 600;
    }
    .sv-evidence-page .sv-toolbar-zoom button[disabled],
    .sv-evidence-page .sv-toolbar-presets select[disabled] {
      opacity: 0.45; cursor: not-allowed;
    }

    /* Clear-selection pill in the locus readout (step 4.5) */
    .sv-evidence-page .sv-clear-selection {
      display: inline-block;
      padding: 1px 6px; margin-left: 4px;
      font-size: 10px; font-weight: 600;
      background: var(--panel-2, #1c2330); color: var(--accent, #f5a524);
      border: 1px solid var(--accent, #f5a524); border-radius: 2px;
      cursor: pointer;
    }
    .sv-evidence-page .sv-clear-selection:hover {
      background: var(--accent, #f5a524); color: #0e1116;
    }

    /* Main 3-column grid */
    .sv-evidence-page .sv-grid {
      display: grid;
      grid-template-columns: 220px 1fr 300px;
      gap: 0;
      flex: 1 1 auto; min-height: 0;
      overflow: hidden;
    }
    .sv-evidence-page .sv-leftrail,
    .sv-evidence-page .sv-rightrail {
      background: var(--panel, #151a22);
      border-right: 1px solid var(--rule, #2a3242);
      overflow-y: auto; overflow-x: hidden;
      padding: 8px 10px;
    }
    .sv-evidence-page .sv-rightrail {
      border-right: none; border-left: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-main {
      display: flex; flex-direction: column;
      min-width: 0; min-height: 0;
      background: var(--bg, #0e1116);
      overflow: auto;
    }

    /* Sections in the rails */
    .sv-evidence-page .sv-section { margin-bottom: 14px; }
    .sv-evidence-page .sv-section-title {
      font-size: 10.5px; font-weight: 700; letter-spacing: 0.06em;
      text-transform: uppercase; color: var(--ink-dimmer, #5a6472);
      margin-bottom: 6px; padding-bottom: 3px;
      border-bottom: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-empty-hint {
      font-size: 11px; color: var(--ink-dim, #8a94a3); font-style: italic;
      padding: 2px 0;
    }

    /* Key/value rows in the rails */
    .sv-evidence-page .sv-row {
      display: flex; justify-content: space-between; align-items: baseline;
      padding: 2px 0; font-size: 11.5px;
    }
    .sv-evidence-page .sv-key   { color: var(--ink-dim, #8a94a3); }
    .sv-evidence-page .sv-val   { color: var(--ink, #e7edf3); font-variant-numeric: tabular-nums; }
    .sv-evidence-page .sv-strong{ font-weight: 600; }

    /* Karyotype chips */
    .sv-evidence-page .sv-karyo-row {
      display: flex; gap: 4px; flex-wrap: wrap;
    }
    .sv-evidence-page .sv-karyo-chip {
      flex: 1 1 0; min-width: 50px;
      display: flex; flex-direction: column; align-items: center;
      gap: 1px; padding: 4px 6px;
      background: var(--panel-2, #1c2330);
      border: 1.5px solid var(--rule, #2a3242);
      border-radius: 4px;
      color: var(--ink, #e7edf3);
      font-size: 11px; font-weight: 600;
      cursor: pointer; transition: background 80ms;
    }
    .sv-evidence-page .sv-karyo-chip:hover { background: var(--panel-3, #2a3242); }
    .sv-evidence-page .sv-karyo-name { font-family: var(--mono, monospace); font-size: 11.5px; }
    .sv-evidence-page .sv-karyo-n    { font-size: 10px; color: var(--ink-dim, #8a94a3); font-weight: 400; }

    /* Zone rows */
    .sv-evidence-page .sv-zone-row {
      display: flex; align-items: center; gap: 6px;
      padding: 2px 0; font-size: 11px;
    }
    .sv-evidence-page .sv-zone-swatch {
      display: inline-block; width: 9px; height: 9px; border-radius: 2px;
      flex: 0 0 auto;
    }
    .sv-evidence-page .sv-zone-label { flex: 1 1 auto; color: var(--ink, #e7edf3); }
    .sv-evidence-page .sv-zone-range { color: var(--ink-dim, #8a94a3); font-variant-numeric: tabular-nums; }

    /* Main empty-state card */
    .sv-evidence-page .sv-locus-strip {
      flex: 0 0 auto;
      min-height: 200px;
      display: flex; align-items: center; justify-content: center;
      padding: 24px;
    }
    .sv-evidence-page .sv-empty-state {
      max-width: 540px; text-align: center;
      padding: 24px 28px;
      background: var(--panel, #151a22);
      border: 1px dashed var(--rule, #2a3242);
      border-radius: 6px;
    }
    .sv-evidence-page .sv-empty-title {
      font-size: 14px; font-weight: 600; color: var(--ink, #e7edf3);
      margin-bottom: 8px;
    }
    .sv-evidence-page .sv-empty-body {
      font-size: 12px; line-height: 1.6; color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-empty-body code {
      font-family: var(--mono, monospace); font-size: 11px;
      background: var(--panel-2, #1c2330); padding: 1px 5px; border-radius: 3px;
      color: var(--ink, #e7edf3);
    }
    .sv-evidence-page .sv-empty-detail {
      font-size: 11px; color: var(--ink-dimmer, #5a6472); font-style: italic;
    }
    .sv-evidence-page .sv-table-wrap {
      flex: 1 1 auto; min-height: 0;
      padding: 0 12px 12px;
    }

    /* Locus track strip (step 3) */
    .sv-evidence-page .sv-locus-strip.sv-locus-has-content {
      align-items: stretch; justify-content: flex-start;
      min-height: 0; padding: 0;
    }
    .sv-evidence-page .sv-locus-wrap { width: 100%; }
    .sv-evidence-page .sv-locus-stack {
      display: flex; flex-direction: column;
      gap: 4px;
    }
    .sv-evidence-page .sv-locus-track {
      display: flex; align-items: center; gap: 6px;
      width: 100%;
    }
    .sv-evidence-page .sv-locus-track-label {
      flex: 0 0 110px; padding-right: 6px;
      font-size: 10.5px; color: var(--ink-dim, #8a94a3);
      text-align: right;
    }
    .sv-evidence-page .sv-locus-track-svg {
      flex: 1 1 auto; min-width: 0;
      line-height: 0; /* SVG inline gap fix */
    }
    .sv-evidence-page .sv-locus-track-svg svg {
      display: block; width: 100%; height: 100%;
    }
    .sv-evidence-page .sv-locus-overlay {
      display: block; width: 100%; height: 100%;
    }
    .sv-evidence-page .sv-locus-glyph:hover { opacity: 1 !important; }
    .sv-evidence-page .sv-locus-legend {
      display: flex; align-items: center; gap: 12px;
      padding: 6px 0 0 116px;
      font-size: 11px; color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-locus-legend-label { color: var(--ink-dimmer, #5a6472); }
    .sv-evidence-page .sv-locus-legend-item {
      display: inline-flex; align-items: center; gap: 4px;
    }
    .sv-evidence-page .sv-locus-legend-glyph { font-size: 13px; }

    /* Cursor readout (step 3.5) */
    .sv-evidence-page .sv-locus-readout {
      padding: 4px 0 8px 116px;
      font-size: 11px; line-height: 1.45;
      color: var(--ink, #e7edf3);
      font-variant-numeric: tabular-nums;
    }
    .sv-evidence-page .sv-readout-strong { font-weight: 600; }
    .sv-evidence-page .sv-readout-dim { color: var(--ink-dim, #8a94a3); }
    .sv-evidence-page .sv-readout-sep {
      color: var(--ink-dimmer, #5a6472); margin: 0 6px;
    }
    .sv-evidence-page .sv-readout-hotkeys {
      font-size: 10px; color: var(--ink-dimmer, #5a6472);
    }
    /* Hit-rect cursor */
    .sv-evidence-page .sv-locus-hit { cursor: crosshair; }

    /* Filters block */
    .sv-evidence-page .sv-filter-row {
      display: flex; align-items: center; gap: 6px;
      padding: 3px 0; font-size: 11px;
    }
    .sv-evidence-page .sv-filter-row label {
      flex: 0 0 84px; color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-filter-row select,
    .sv-evidence-page .sv-filter-row input[type="range"] {
      flex: 1 1 auto;
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
      border: 1px solid var(--rule, #2a3242); border-radius: 3px;
      padding: 2px 6px; font-size: 11px;
    }
    .sv-evidence-page .sv-filter-row input[type="range"] { padding: 0; }
    .sv-evidence-page .sv-filter-row .sv-filter-num {
      flex: 0 0 28px; text-align: right;
      font-variant-numeric: tabular-nums; color: var(--ink, #e7edf3);
    }
    .sv-evidence-page .sv-filter-checkrow { gap: 8px; }
    .sv-evidence-page .sv-filter-checkrow label {
      flex: 1 1 auto; color: var(--ink, #e7edf3);
    }
    .sv-evidence-page .sv-filter-actions {
      display: flex; gap: 6px; margin-top: 8px;
    }
    .sv-evidence-page .sv-btn-primary,
    .sv-evidence-page .sv-btn-secondary {
      flex: 1 1 0;
      padding: 4px 8px; font-size: 11px;
      border-radius: 3px; cursor: pointer;
      border: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-btn-primary {
      background: var(--accent-2, #4fa3ff); color: #0e1116; border-color: var(--accent-2, #4fa3ff);
      font-weight: 600;
    }
    .sv-evidence-page .sv-btn-secondary {
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
    }
    .sv-evidence-page .sv-btn-secondary:hover { background: var(--panel-3, #2a3242); }
    .sv-evidence-page .sv-filter-hint {
      margin-top: 8px;
      font-size: 10.5px; line-height: 1.4;
      color: var(--ink-dimmer, #5a6472);
    }
    .sv-evidence-page .sv-filter-hint code {
      font-family: var(--mono, monospace); font-size: 10px;
      background: var(--panel-2, #1c2330); padding: 0 4px; border-radius: 2px;
      color: var(--ink-dim, #8a94a3);
    }

    /* Table */
    .sv-evidence-page .sv-tabletop { padding: 8px 0 4px; }
    .sv-evidence-page .sv-tablebar {
      display: flex; align-items: center; gap: 8px;
      padding: 6px 8px;
      background: var(--panel, #151a22);
      border: 1px solid var(--rule, #2a3242);
      border-radius: 4px 4px 0 0;
      border-bottom: none;
    }
    .sv-evidence-page .sv-tablebar-label {
      font-weight: 600; color: var(--ink, #e7edf3); font-size: 11.5px;
    }
    .sv-evidence-page .sv-tablebar-tabs { display: flex; gap: 0; }
    .sv-evidence-page .sv-tablebar-spacer { flex: 1 1 auto; }
    .sv-evidence-page .sv-tab {
      background: none; border: none; border-bottom: 2px solid transparent;
      color: var(--ink-dim, #8a94a3);
      padding: 4px 10px; font-size: 11px; cursor: pointer;
    }
    .sv-evidence-page .sv-tab.active {
      color: var(--ink, #e7edf3);
      border-bottom-color: var(--accent-2, #4fa3ff);
    }
    .sv-evidence-page .sv-tab[disabled] { opacity: 0.5; cursor: not-allowed; }
    .sv-evidence-page .sv-table-scroll {
      overflow: auto;
      border: 1px solid var(--rule, #2a3242);
      max-height: 56vh;
      background: var(--panel, #151a22);
    }
    .sv-evidence-page .sv-table {
      width: 100%; border-collapse: collapse;
      font-size: 11px; font-variant-numeric: tabular-nums;
    }
    .sv-evidence-page .sv-table thead {
      position: sticky; top: 0; z-index: 1;
      background: var(--panel, #151a22);
    }
    .sv-evidence-page .sv-th-group {
      padding: 4px 6px; text-align: center;
      font-size: 10px; font-weight: 600;
      color: var(--ink-dim, #8a94a3);
      border-bottom: 1px solid var(--rule, #2a3242);
      background: var(--panel-2, #1c2330);
    }
    .sv-evidence-page .sv-th-group-divider { border-left: 1px solid var(--rule, #2a3242); }
    .sv-evidence-page .sv-th-col {
      padding: 4px 6px;
      font-size: 10.5px; font-weight: 600;
      color: var(--ink, #e7edf3);
      border-bottom: 1px solid var(--rule, #2a3242);
      cursor: pointer; user-select: none;
      white-space: nowrap;
    }
    .sv-evidence-page .sv-th-col:hover { background: var(--panel-2, #1c2330); }
    .sv-evidence-page .sv-th-left   { text-align: left; }
    .sv-evidence-page .sv-th-right  { text-align: right; }
    .sv-evidence-page .sv-th-center { text-align: center; }
    .sv-evidence-page .sv-td {
      padding: 3px 6px; white-space: nowrap;
      border-bottom: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-td-left   { text-align: left; }
    .sv-evidence-page .sv-td-right  { text-align: right; }
    .sv-evidence-page .sv-td-center { text-align: center; }
    .sv-evidence-page .sv-tr:hover {
      background: var(--panel-2, #1c2330); cursor: pointer;
    }
    .sv-evidence-page .sv-row-hi {
      box-shadow: inset 3px 0 0 0 var(--accent-2, #4fa3ff);
    }
    .sv-evidence-page .sv-typedot {
      font-size: 10px; vertical-align: middle;
    }
    .sv-evidence-page .sv-id-btn {
      background: none; border: none; padding: 0; color: var(--accent-2, #4fa3ff);
      cursor: pointer; font: inherit;
    }
    .sv-evidence-page .sv-id-btn:hover { text-decoration: underline; }
    /* Row annotation chips (step 4) — replaces the old single-star icon.
       Cycle: ·  → L → R → ★ → ·  on click. */
    .sv-evidence-page .sv-anno-btn {
      background: none; border: none; padding: 0 4px 0 0;
      cursor: pointer; font: inherit;
      vertical-align: middle;
    }
    .sv-evidence-page .sv-anno {
      display: inline-block; min-width: 14px;
      padding: 0 3px;
      font-family: var(--mono, monospace); font-size: 10.5px; font-weight: 700;
      text-align: center; border-radius: 2px;
      vertical-align: middle;
    }
    .sv-evidence-page .sv-anno-empty {
      color: var(--ink-dimmer, #5a6472); font-weight: 400;
    }
    .sv-evidence-page .sv-anno-L {
      color: #fff; background: #3cc08a;
    }
    .sv-evidence-page .sv-anno-R {
      color: #fff; background: #d08770;
    }
    .sv-evidence-page .sv-anno-star {
      color: var(--accent, #f5a524); background: rgba(245,165,36,0.12);
    }
    .sv-evidence-page .sv-row-pattern-hi {
      box-shadow: inset 3px 0 0 0 var(--accent, #f5a524);
      background: rgba(245,165,36,0.04);
    }

    /* Right-rail boundary summary (step 4) */
    .sv-evidence-page .sv-bsblock { margin-bottom: 12px; }
    .sv-evidence-page .sv-bsblock-title {
      font-size: 11.5px; font-weight: 600;
      margin-bottom: 4px; padding-bottom: 2px;
      display: flex; align-items: baseline; gap: 6px;
      border-bottom: 1.5px solid currentColor;
    }
    .sv-evidence-page .sv-bsblock-left  { color: #4fa3ff; }
    .sv-evidence-page .sv-bsblock-right { color: #bf616a; }
    .sv-evidence-page .sv-bsblock-interval {
      font-size: 10px; font-weight: 400; color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-bs-table {
      width: 100%; border-collapse: collapse;
      font-size: 11px; font-variant-numeric: tabular-nums;
    }
    .sv-evidence-page .sv-bs-table th {
      text-align: left; font-weight: 600;
      color: var(--ink-dim, #8a94a3); font-size: 10px;
      padding: 2px 4px;
      border-bottom: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-bs-table th:nth-child(2),
    .sv-evidence-page .sv-bs-table th:nth-child(3),
    .sv-evidence-page .sv-bs-table .sv-bs-num { text-align: right; }
    .sv-evidence-page .sv-bs-row {
      cursor: pointer;
    }
    .sv-evidence-page .sv-bs-row:hover {
      background: var(--panel-2, #1c2330);
    }
    .sv-evidence-page .sv-bs-row td {
      padding: 2px 4px;
    }
    .sv-evidence-page .sv-bs-type { font-weight: 600; }
    .sv-evidence-page .sv-bs-total {
      font-weight: 600; color: var(--ink, #e7edf3);
      border-top: 1px solid var(--rule, #2a3242);
    }
    .sv-evidence-page .sv-bs-total td { padding: 3px 4px; }

    /* Right-rail legend (step 4) */
    .sv-evidence-page .sv-legend-cols {
      display: flex; gap: 12px;
    }
    .sv-evidence-page .sv-legend-col { flex: 1 1 0; min-width: 0; }
    .sv-evidence-page .sv-legend-coltitle {
      font-size: 10px; font-weight: 700; letter-spacing: 0.04em;
      text-transform: uppercase;
      color: var(--ink-dimmer, #5a6472);
      margin-bottom: 4px;
    }
    .sv-evidence-page .sv-legend-item {
      display: flex; align-items: center; gap: 5px;
      padding: 2px 0; font-size: 10.5px;
      cursor: pointer; border-radius: 2px;
    }
    .sv-evidence-page .sv-legend-item:hover { background: var(--panel-2, #1c2330); }
    .sv-evidence-page .sv-legend-swatch {
      display: inline-block; width: 9px; height: 9px;
      border-radius: 50%; flex: 0 0 auto;
    }
    .sv-evidence-page .sv-legend-text { white-space: nowrap; }

    /* Step 5 — UpSet panel */
    .sv-evidence-page .sv-upset { font-size: 11px; }
    .sv-evidence-page .sv-upset-header {
      display: flex; gap: 6px; align-items: center;
      margin-bottom: 4px; flex-wrap: wrap;
    }
    .sv-evidence-page .sv-upset-pill {
      margin-left: auto;
      padding: 1px 6px;
      font-size: 10px; font-weight: 600;
      background: var(--accent, #f5a524); color: #0e1116;
      border-radius: 2px; cursor: pointer;
    }
    .sv-evidence-page .sv-upset-pill:hover { opacity: 0.85; }
    .sv-evidence-page .sv-upset-svg {
      display: block;
      max-width: 100%;
    }
    .sv-evidence-page .sv-upset-bar-g .sv-upset-bar-hit:hover ~ rect {
      filter: brightness(1.15);
    }
    .sv-evidence-page .sv-upset-foothint {
      margin-top: 4px;
      font-size: 9.5px;
    }
    .sv-evidence-page .sv-upset-emptyhint code {
      background: var(--panel-2, #1c2330);
      padding: 1px 4px; border-radius: 2px;
      font-family: var(--mono, monospace);
      font-size: 10px;
    }
    .sv-evidence-page .sv-upset-error {
      color: var(--bad, #e0555c);
    }

    /* Step 6 — heatmap CSS */
    .sv-evidence-page .sv-heatmap { font-size: 11px; }
    .sv-evidence-page .sv-hm-header {
      display: flex; align-items: center; gap: 8px;
      margin-bottom: 4px; flex-wrap: wrap;
      font-size: 10.5px;
    }
    .sv-evidence-page .sv-hm-leg {
      display: inline-flex; align-items: center; gap: 4px;
      color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-hm-sw {
      display: inline-block; width: 9px; height: 9px;
      border-radius: 1px; vertical-align: middle;
    }
    .sv-evidence-page .sv-hm-expand {
      margin-left: auto;
      padding: 1px 6px;
      font-size: 10px; font-weight: 500;
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
      border: 1px solid var(--rule, #2a3242); border-radius: 2px;
      cursor: pointer;
    }
    .sv-evidence-page .sv-hm-expand:hover { background: var(--accent-2, #4fa3ff); color: #0e1116; }
    .sv-evidence-page .sv-hm-svgwrap { max-width: 100%; overflow-x: auto; }
    .sv-evidence-page .sv-hm-svg { display: block; }
    .sv-evidence-page .sv-hm-rowlabel {
      font-family: var(--mono, monospace);
      fill: var(--ink-dim, #8a94a3);
      cursor: pointer;
    }
    .sv-evidence-page .sv-hm-rowlabel-hi {
      fill: var(--accent-2, #4fa3ff); font-weight: 700;
    }
    .sv-evidence-page .sv-hm-rowlabel:hover { fill: var(--ink, #e7edf3); }
    .sv-evidence-page .sv-hm-tip {
      padding: 3px 6px;
      background: var(--panel-2, #1c2330);
      border-left: 2px solid var(--accent, #f5a524);
      font-family: var(--mono, monospace);
      font-size: 10.5px;
      margin-bottom: 4px;
    }
    .sv-evidence-page .sv-hm-tip code { color: var(--accent-2, #4fa3ff); }
    .sv-evidence-page .sv-hm-emptyhint code {
      background: var(--panel-2, #1c2330);
      padding: 1px 4px; border-radius: 2px;
      font-family: var(--mono, monospace);
      font-size: 10px;
    }
    .sv-evidence-page .sv-hm-error { color: var(--bad, #e0555c); }
    /* Large overlay panel — sits on top of the locus area */
    .sv-evidence-page .sv-hm-overlay {
      position: absolute; top: 56px; left: 232px; right: 312px;
      max-height: calc(100vh - 140px);
      overflow: auto;
      background: var(--panel, #151a22);
      border: 1px solid var(--rule, #2a3242);
      border-radius: 4px;
      padding: 12px 16px;
      box-shadow: 0 4px 24px rgba(0,0,0,0.5);
      z-index: 100;
    }
    .sv-evidence-page .sv-table-empty {
      padding: 32px 24px; text-align: center;
      color: var(--ink-dim, #8a94a3);
    }

    /* Footer + pager */
    .sv-evidence-page .sv-table-footer {
      display: flex; flex-direction: column; gap: 6px;
      padding: 8px 4px;
      font-size: 11px;
      color: var(--ink-dim, #8a94a3);
    }
    .sv-evidence-page .sv-pager {
      display: flex; align-items: center; gap: 8px; flex-wrap: wrap;
    }
    .sv-evidence-page .sv-pager select {
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
      border: 1px solid var(--rule, #2a3242); border-radius: 3px;
      padding: 1px 4px; font-size: 11px;
    }
    .sv-evidence-page .sv-pager-range { margin-right: 6px; }
    .sv-evidence-page .sv-pager-btn {
      background: var(--panel-2, #1c2330); color: var(--ink, #e7edf3);
      border: 1px solid var(--rule, #2a3242); border-radius: 3px;
      padding: 1px 7px; font-size: 11px; cursor: pointer;
      min-width: 22px;
    }
    .sv-evidence-page .sv-pager-btn:hover { background: var(--panel-3, #2a3242); }
    .sv-evidence-page .sv-pager-btn[disabled] { opacity: 0.4; cursor: not-allowed; }
    .sv-evidence-page .sv-pager-btn.sv-pager-active {
      background: var(--accent-2, #4fa3ff); color: #0e1116;
      border-color: var(--accent-2, #4fa3ff); font-weight: 600;
    }
    .sv-evidence-page .sv-pager-ellipsis { color: var(--ink-dimmer, #5a6472); padding: 0 2px; }
    .sv-evidence-page .sv-table-legend-line {
      font-size: 10.5px; color: var(--ink-dimmer, #5a6472);
    }

    /* Drag-drop active state — overlay across the whole page wrapper */
    .sv-evidence-page.sv-drop-active::after {
      content: 'Drop sv_genotype_counts JSON to load…';
      position: absolute; inset: 0;
      display: flex; align-items: center; justify-content: center;
      background: rgba(79,163,255,0.15);
      color: var(--ink, #e7edf3);
      border: 3px dashed var(--accent-2, #4fa3ff);
      font-size: 16px; font-weight: 600;
      pointer-events: none;
      z-index: 10;
    }
    .sv-evidence-page { position: relative; }
  `;

  function _injectStylesOnce() {
    if (document.getElementById('atlas-sv-evidence-styles')) return;
    const tag = document.createElement('style');
    tag.id = 'atlas-sv-evidence-styles';
    tag.textContent = STYLES;
    document.head.appendChild(tag);
  }

  // ===========================================================================
  // Exports
  // ===========================================================================

  return {
    init,
    loadCandidate,
    refresh,
    setFilters,
    exportFilteredTSV,
    // step 3.5: cursor + hotkey API
    attachHotkeys: _attachLocusHotkeys,
    detachHotkeys: _detachLocusHotkeys,
    setViewPreset: _setViewPreset,
    setCursorBp: function (bp) { _state.cursorBp = bp; _renderLocusStrip(); },
    // step 4.5: select-mode API
    setSelectMode: _setSelectMode,
    setSelection: function (startBp, endBp) {
      if (startBp == null || endBp == null) {
        _state.selection = null;
      } else {
        _state.selection = { startBp: Math.min(startBp, endBp),
                             endBp:   Math.max(startBp, endBp) };
      }
      _state.pageIndex = 0;
      _renderLocusStrip();
      _renderTable();
    },
    // step 5: UpSet panel + folder-walk + sample selection API
    onUpSetBarClick: _onUpSetBarClick,
    clearSampleSelection: _clearSampleSelection,
    setCombinationsLayer: function (obj) {
      const v = _validateCombinationsLayer(obj);
      if (!v.ok) { _state.layerError = 'invalid combinations layer: ' + v.reason; return false; }
      _state.combinationsLayer = obj;
      _renderRightRail();
      return true;
    },
    // step 6: support layer + heatmap API
    setSupportLayer: function (obj) {
      const v = _validateSupportLayer(obj);
      if (!v.ok) { _state.layerError = 'invalid support layer: ' + v.reason; return false; }
      _state.supportLayer = obj;
      // If a sample selection is already active, the recount becomes
      // available now that support is present. Recompute and re-render
      // the table so its counts switch to the in-selection view.
      if (_state.selectedSamples && _state.selectedSamples.size > 0) {
        _recomputeGtCountsView();
      }
      _renderRightRail();
      _renderLocusStrip();
      _renderTable();
      return true;
    },
    setHeatmapView: function (mode) {
      if (mode !== 'compact' && mode !== 'large') return;
      _state.heatmapView = mode;
      _renderRightRail();
      if (mode === 'large') _renderHeatmapLarge();
      else {
        const ov = _state.rootEl && _state.rootEl.querySelector('[data-sv-hm-overlay]');
        if (ov && ov.parentNode) ov.parentNode.removeChild(ov);
      }
    },
    _state,
    // Constants exposed for the upcoming step-2 table renderer + tests.
    _const: {
      DEFAULT_FILTERS,
      FDR_COLOURS,
      PATTERN_LABELS,
      SV_TYPE_STYLE,
      KARYOTYPE_REMAP,
      GROUP_COLOURS,
      TABLE_COLUMNS,
      TABLE_GROUP_HEADERS,
      PAGE_SIZES,
      DEFAULT_PAGE_SIZE,
      LOCUS_TRACKS,
      SV_TYPE_CYCLE,
      UPSET,
      HEATMAP,
    },
    _internals: {
      _resolveActiveCandidateId,
      _resolveKaryotypeGroups,
      _readFilters, _writeFilters,
      _readAnnotations, _writeAnnotations, _cycleAnnotation,
      _fdrColour, _fmtBp, _fmtMb, _fmtSci,
      _zoneDisplay,
      // step 2:
      _applyFilters, _sortRows, _paginateRows, _computeVisibleRows,
      _buildFilteredTSV,
      // step 3:
      _drawTrackZones, _drawTrackAxis, _drawTrackSVCalls, _drawTrackPlaceholder,
      // step 3.5:
      _windowForPreset, _nearestSV, _stepCursor, _jumpCursorToNextSV,
      _cycleSvTypeFilter, _setViewPreset, _onZoomButton,
      // step 4:
      _renderAnnoChip, _renderBoundarySummaryHTML, _renderLegendHTML,
      _stepHighlightedRow,
      // step 4.5:
      _setSelectMode, _commitSelectDrag, _refreshModeButtonStyles,
      // step 5:
      _validateCombinationsLayer, _renderUpSetPanel,
      _onUpSetBarClick, _clearSampleSelection,
      _ingestJsonText,
      // step 6:
      _validateSupportLayer, _supportCell, _renderHeatmapPanel,
      _renderHeatmapHTML, _renderHeatmapLarge, _hasNoCarrierInSelection,
      _recomputeGtCountsView,
    },
  };
}));
