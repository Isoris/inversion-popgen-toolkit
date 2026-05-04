# Atlas audit — full integration spec for the six-axis classification

**Source of truth read:** `Inversion_atlas.html` (~61 570 lines, atlas v4.t124 by inspection).

**Scope of audit:** every page that today carries classification semantics, every JSON layer the atlas consumes, every export schema the atlas emits, and every place where the **six-axis classification** (Evidence / Structure / Mechanism / History / Population / Function & Application) needs to slot in **without breaking the existing read-only contract** between cluster-side R pipelines and the browser-side scrubber.

The atlas explicitly distinguishes WORKFLOW (the Q1–Q7 question set, §13) from OUTPUT CATALOG (the 14-axis `final_classification`, §19). The audit respects that split: **the six-axis framework belongs at §19, NOT §13.** §13 is the implementation backend; §19 is what the user sees and what the manuscript reports.

---

## 0. What the atlas already does well — DO NOT TOUCH

These pieces of the atlas are already aligned with the six-axis framework even though they predate it. Leave them alone; they only need re-labelling at the most superficial level.

### 0.1 Page 4 (`page4`) Tier sub-view — THE MAIN INTEGRATION TARGET

- Lines 49060–49200 define the 14-axis schema as a JS constant `_TIER_AXES` and consume `state.data.final_classification[<cid>]` (an object keyed by candidate_id, values are per-axis category labels).
- The render is read-only by design (`scrubber never writes to final_classification`, line 8096).
- Empty-state UI is already wired (line 49200+: `layerStatus.textContent = 'final_classification: loaded · ...'`).

**This is the surface that needs to migrate from 14 axes to 6.** Migration is mechanical because the rendering is data-driven from `_TIER_AXES`. See §3 below.

### 0.2 Page 16b (`page16b`) Multi-species cockpit — THE AGE/HISTORY ENGINE

- Function `_msAutoSuggestArchitecture` (line 21626) classifies breakpoints into A–F (simple inversion / synteny-boundary / fusion-fission / terminal translocation / recurrent hotspot / ambiguous).
- Function `_msAutoSuggestAgeModel` (line 21850) classifies inversions into 5 age models: `YOUNG-POP`, `OLD-POLY`, `OLD-BP-YOUNG-INV`, `LINEAGE-KARYO`, `MULTI-AGE-HOTSPOT`.
- Already consumes `dxy_per_inversion_v1.json` (line ~21880 — `signals.fold_elevation`, `signals.dxy_within`) and `synteny_multispecies_v1.json` + `phylo_tree_v1.json`.
- Already exports per-bp TSV with columns `architecture_class, architecture_source, architecture_confidence, age_model, age_source, age_confidence` (line 22096).

**This is the most sophisticated classification surface in the atlas already.** It just needs to feed its outputs into the unified six-axis structure rather than living as an isolated cockpit. See §6 below.

### 0.3 Page 17 (`page17`) Stats profile — already category-organised

- Lines 23475–23625 define `SP_DEFAULT_ROWS` with 12 statistics under 5 categories: Breakpoint architecture / Genomic composition / Functional cargo / Population variation / Breeding utility.
- Each row already declares `test`, `null_comparison`, `derive_from`, and a `african_hint` / `bighead_hint` placeholder for cross-species rendering.

**This page is already 90% aligned with the six-axis framework — it just uses different category names.** Renaming + 2 added categories closes the gap. See §7 below.

### 0.4 Page 18 (`page18`) Marker readiness — already six-axis Function & Application

- Lines 7188–7194 + 4595–4604 define Tier 1–4 with explicit AF thresholds. This **is** Axis 6.5 (Application / marker readiness) of the six-axis scheme.
- Already consumes `variant_afs.json` + optional `marker_panel.json` overlay.

**No migration needed; just relabel the page header to "Axis 6 · Application — marker readiness" and document the alignment.**

### 0.5 Catalogue (`page3`, `CAT_COLS_BASE` line 58494)

- Already a flat per-candidate spreadsheet with sortable/filterable columns. New axis columns added here are read by the bundle export automatically (line 47100+).

**This is the right place to add `axis_*_summary` columns.** See §4 below.

---

## 1. The atlas's current classification surfaces, mapped to the six axes

| # | Atlas surface | Code location | What it carries | Six-axis mapping | Status |
|---|---|---|---|---|---|
| 1 | Existence layers A–D | `_TIER_AXES[0..3]` (line 49093) | pass/fail per layer (Local PCA, SV, GHSL, Fisher) | **Axis 1 — Evidence (1.1–1.4)** | already there, rename only |
| 2 | Boundary quality | `_TIER_AXES[4]` | sharp / fuzzy | **Axis 2 — Structure (2.5)** | rename |
| 3 | Group validation | `_TIER_AXES[5]` | NONE / UNCERTAIN / SUPPORTED / VALIDATED / SUSPECT | **Axis 1 — Evidence (1.2 cluster quality)** | split — currently treated as its own axis but biologically a sub-axis of Evidence |
| 4 | Internal structure | `_TIER_AXES[6]` | clean / gradient / composite_undecomposed | **Axis 5 — Population dynamics (5.5 within-arrangement substructure)** | moved |
| 5 | Recombinant class | `_TIER_AXES[7]` | none / gene_conversion / double_crossover / mixed | **Axis 3 — Mechanism (3.4 gene-flux signature)** | moved |
| 6 | Family linkage | `_TIER_AXES[8]` | multi_family / few_family / single_family / pca_family_confounded | **Axis 5 — Population dynamics (5.3 family confound)** | moved |
| 7 | Polymorphism class | `_TIER_AXES[9]` | cohort_wide / lineage_restricted / family_restricted | **Axis 4 — History (4.3 Type I/II)** + **Axis 5 (5.4 cline)** | split |
| 8 | Mechanism class | `_TIER_AXES[10]` | NAHR / NHEJ / MMBIR | **Axis 3 — Mechanism (3.1 breakpoint substrate)** | rename |
| 9 | Age class | `_TIER_AXES[11]` | young / intermediate / ancient | **Axis 4 — History (4.1 age)** | rename |
| 10 | Burden class | `_TIER_AXES[12]` | enriched / neutral / depleted | **Axis 6 — Function & Application (6.1 cargo + load)** | moved |
| 11 | Confidence tier | `_TIER_AXES[13]` | T1 / T2 / T3 / T4 | **summary across Axis 1** | keep as composite |
| 12 | Age models (page 16b) | `_msAutoSuggestAgeModel` (line 21850) | YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV / LINEAGE-KARYO / MULTI-AGE-HOTSPOT | **Axis 4 — History (4.1+4.4 recurrence)** | richer than `_TIER_AXES[11]` — promote |
| 13 | Architecture class A–F (page 16b) | `_msAutoSuggestArchitecture` (line 21626) | simple / synteny-boundary / fusion-fission / terminal / recurrent / ambiguous | **Axis 2 — Structure (2.4 topology)** | promote |
| 14 | Regime topology (regime registry) | `REGIME_AXIS_TOPOLOGY` (line 11532) | one_axis_3band / one_axis_nested / two_axes_independent / two_axes_correlated / multiband_stable / artifact_suspect | **Axis 1 — Evidence (1.2 cluster quality)** | keep — this is cluster-quality vocabulary |
| 15 | Stats-profile categories (page 17) | `SP_DEFAULT_ROWS` (line 23475) | Breakpoint architecture / Genomic composition / Functional cargo / Population variation / Breeding utility | **rename to Axes 2–6** | rename + extend |
| 16 | Marker tier 1–4 (page 18) | line 7188 | Private indel / multi-marker / breakpoint PCR demoted / exploratory | **Axis 6 — Function & Application (6.5)** | rename |
| 17 | Negative regions (page 6) | `negative_regions.json` | inversion_candidate / complex_candidate / no_detectable_inversion_high_confidence / no_detectable_inversion_low_power / not_testable | **Axis 1 — Evidence (1.1–1.4 inverted = "all layers fail")** | keep — orthogonal use of same axis |
| 18 | Cross-species breakpoints (page 13) | `cs_breakpoints_v1.json` | event_type per bp + TE flank density | feeds **Axis 2 (2.4 topology), Axis 3 (3.2 substrate), Axis 4 (4.4 recurrence)** | keep — multi-axis source |
| 19 | dxy_per_inversion_v1 | layer | `fold_elevation_inside_vs_flank`, `dxy_within_inversion_ref_vs_inv` | **Axis 4 — History (4.1)** | keep — feeds age engine |
| 20 | Variant AFs (page 18) | `variant_afs.json` | per-variant AF in STD/HET/INV groups | **Axis 6 — Function & Application (6.5)** | keep |

**Bottom line**: 18 of 20 surfaces already exist and just need relabelling or moving to the right axis. Two genuinely new things are needed: (a) the **six-card layout on page 4** that replaces the 14-axis grid, and (b) the **back-compat alias layer in `_TIER_AXES`** that makes existing R-pipeline outputs keep working.

---

## 2. Critical findings from the audit

### 2.1 The atlas's vocabulary is right, the categorisation is wrong

The 14-axis schema is biologically mostly correct, but it **mixes evidence-level claims with biology-level claims** in a way the six-axis framework explicitly forbids (Rule 2 from the post-classification spec: "stratify features by biological level").

Specifically:
- "Group validation" (axis 5) is an **Evidence-level** claim (does the K=3 PCA cluster hold up?), not a separate axis.
- "Internal structure" (axis 7) is a **Population-level** claim (is there sub-arrangement diversity within HOMO_1 or HOMO_2?), not a Group claim.
- "Recombinant class" (axis 8) is a **Mechanism-level** claim (gene flux mechanism), not a Group claim.
- "Polymorphism class" (axis 10) packs **two distinct things** — (i) is the inversion segregating cohort-wide vs lineage-restricted (a **History** claim, Type I vs Type II from Faria 2019), and (ii) is it cline-correlated (a **Population** claim from Wellenreuther 2018).

The migration to six axes resolves these confusions.

### 2.2 The atlas has TWO age signals that don't talk to each other

- `_TIER_AXES[11].age_class` = `young / intermediate / ancient / unknown` — coarse, R-pipeline-emitted, currently **not yet shipped** ("renders empty until final_classification ships").
- `_msAutoSuggestAgeModel` on page 16b = `YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV / LINEAGE-KARYO / MULTI-AGE-HOTSPOT` — fine-grained, browser-computed live, **already working** when `dxy_per_inversion_v1` + `synteny_multispecies_v1` are loaded.

These need to merge. The page-16b version is richer (5 categories vs 3) and grounded in the cross-species evidence that the manuscript anyway uses. Recommendation: **the page-16b age model becomes the canonical Axis 4 value**; the legacy `age_class` becomes a coarse alias derived from it.

### 2.3 The 14-axis schema is missing Axis 6 (Function & Application) almost entirely

The 14 axes have **one** function-related slot (`burden_class` = enriched/neutral/depleted) and zero application/marker slots. The marker tier on page 18 is a separate page with no link back to page 4's tier view.

Recommendation: **add Axis 6 sub-axes explicitly** (cargo class, GO/KEGG enrichment, candidate-gene presence, trait association, marker tier) to the canonical schema. The six-card layout makes this natural — Axis 6 becomes one of the six cards.

### 2.4 The atlas has NO Mechanism-level recombination/gene-flux readout

Axis 3 of the six-axis framework includes:
- 3.1 Breakpoint substrate (NAHR / NHEJ / MMBIR / unknown) — **the atlas has this** as `mechanism_class`.
- 3.2 Repeat-class enrichment at breakpoints — **the atlas has this** via cs_breakpoints flank density.
- 3.3 Recurrent vs single origin — **the atlas has this** via Twisst-style cross-species reuse on page 16b.
- 3.4 Recombination-suppression profile + gene-flux signature — **the atlas has `recombinant_class`** but it's currently classified as a "groups" axis, not a Mechanism axis.

Recommendation: move `recombinant_class` into Axis 3 as `axis3_recombination_class`; rename to be unambiguous. Also add `axis3_dxy_profile_shape` (flat/linear/U-shape from §1.6 of the previous review).

### 2.5 The atlas does NOT yet consume any of the TIER-1 missing tests from the previous review

Recall the six TIER-1 gaps:
1. Breakpoint-restricted dXY + U-shape profile classification
2. Twisst topology weighting per inversion
3. ABC coalescent age estimate (fastsimcoal2/msprime)
4. Matched-null GO + cargo enrichment
5. R_xy formal load ratio with CI
6. pyrho recombination-rate landscape

None of these have a JSON layer slot in the atlas. They need to be added. See §5 below for the per-test integration spec.

### 2.6 The bundle export (`atlas_candidate_export_v2`, line 47097) does NOT carry the axis classification

Currently the export carries: candidate coordinates + per-sample band assignments + haplotype labels + atlas analytics (diamond, sigma, inheritance). It does **not** carry any of the six-axis classification, because that is consumed read-only from `final_classification` on page 4 and the export was designed before the tier sub-view.

Recommendation: **bump to `atlas_candidate_export_v3`** that includes a top-level `classification` block with the six-axis values (auto-extracted from `state.data.final_classification[<cid>]` at export time). See §8 below for the schema.

### 2.7 The catalogue does NOT show any axis at a glance

Currently `CAT_COLS_BASE` (line 58494) shows: chrom, span, K, silhouette, coherence, family-purity, σ-verdict, diamond count, power, regimes. Only `regimes` is an explicit classification column, and it's user-curated, not auto-derived.

Recommendation: add **six axis-summary columns** (one per axis) that read from `final_classification` at row-render time. See §4 below.

### 2.8 The atlas is not species-portable yet

Hard-coded references throughout to *C. gariepinus* / *C. macrocephalus* (page 16, 16b, 18). The six-axis framework is portable; the atlas is not.

Recommendation: define a top-level `state.species` object with `focal_species_label`, `focal_species_short` (default `"Cgar"`), `outgroup_species_label`, `outgroup_species_short` (default `"Cmac"`). All page text reads from these. The bundle import auto-fills them from a `species_metadata.json` layer.

---

## 3. Page 4 (Karyotype/Tier) — the central refactor

### 3.1 Replace the 14-axis grid with a six-card layout

**Current** (line 49093, `_TIER_AXES`): a flat 14-row pill list grouped into 6 visual sections (existence / boundary / groups / population / biology / tier).

**Target**: six **Axis Cards** (one per axis) each showing a headline value + 2–4 sub-values + the underlying tests + their pass/fail + a "see source layer" link. Each card collapses by default to one-line, expands on click to show sub-values and supporting tests.

#### 3.1.1 The six cards

```
+- AXIS 1 - EVIDENCE -----------------------------------------+
| Layer: A.B.C.D OK,OK,OK,-  |  Group validation: VALIDATED   |
| Cluster quality: silhouette 0.62, K=3 stable                |
| Confidence tier: T1                                         |
+-------------------------------------------------------------+
+- AXIS 2 - STRUCTURE ----------------------------------------+
| Size: 4.2 Mb (medium), Position: paracentric, Sub-tel: yes  |
| Topology: simple (page-16b class A), Boundary: sharp        |
| coords: chr28:21,400,000-25,600,000                         |
+-------------------------------------------------------------+
+- AXIS 3 - MECHANISM ----------------------------------------+
| Substrate: NAHR (LTR retrotransposon flanks), TE fold 2.4x  |
| Recombination: suppressed in HET, rho trough confirmed pyrho|
| Gene flux: U-shape dXY profile (gene conversion ~mid)       |
+-------------------------------------------------------------+
+- AXIS 4 - HISTORY ------------------------------------------+
| Age: OLD-POLY (page-16b), dXY 1.8x flank, Tajima D HOMO_1+1.3|
| Polarisation: derived (Cmac=ancestral, Dollo), ABC: 0.8 Ma  |
| Recurrence: shared with C. mac at boundary L (Twisst 0.78)  |
+-------------------------------------------------------------+
+- AXIS 5 - POPULATION DYNAMICS ------------------------------+
| Frequency: 0.32 (intermediate), HWE: het excess (p<0.01)    |
| Family confound: multi_family (NAToRA-pruned re-test OK)    |
| Substructure: BB has 2 sub-clusters (silhouette + gap k=2)  |
+-------------------------------------------------------------+
+- AXIS 6 - FUNCTION & APPLICATION ---------------------------+
| Cargo: 18 genes, GO enriched: immune, oxidative-stress      |
| Candidate gene: tlr5, Trait assoc: pending phenotypes       |
| Marker tier: T1 (private indel rs..., controls OK)          |
+-------------------------------------------------------------+
```

#### 3.1.2 Concrete code change — `_TIER_AXES` becomes `_AXIS_CARDS`

Replace the flat list at line 49093 with a hierarchical structure:

```javascript
// six-axis classification cards (replaces _TIER_AXES, scrubber-side render
// of state.data.final_classification[<cid>]). Each card declares its
// sub-axes; renderer walks the cards and produces collapsible UI. The
// order of cards is fixed for cross-candidate visual consistency.
const _AXIS_CARDS = [
  { axis: 1, label: 'Evidence', short: 'EVID', color: '#4fa3ff',
    sub: [
      { id: 'axis1_layer_a', label: 'Layer A - Local PCA',
        cats: ['pass','fail','unknown'],
        // back-compat alias: writes from existing final_classification.existence_layer_a
        legacy_alias: 'existence_layer_a' },
      { id: 'axis1_layer_b', label: 'Layer B - SV callers',
        cats: ['pass','fail','unknown'], legacy_alias: 'existence_layer_b' },
      { id: 'axis1_layer_c', label: 'Layer C - GHSL haplotype',
        cats: ['pass','fail','unknown'], legacy_alias: 'existence_layer_c' },
      { id: 'axis1_layer_d', label: 'Layer D - Fisher genotype-bp',
        cats: ['pass','fail','unknown'], legacy_alias: 'existence_layer_d' },
      { id: 'axis1_cluster_quality', label: 'Cluster quality',
        cats: ['VALIDATED','SUPPORTED','UNCERTAIN','SUSPECT','NONE'],
        legacy_alias: 'group_validation' },
      { id: 'axis1_regime_topology', label: 'Regime topology',
        cats: ['one_axis_3band','one_axis_nested','two_axes_independent',
               'two_axes_correlated','multiband_stable','artifact_suspect','unspecified'],
        // populated from state.regimeRegistry, not from final_classification
        source: 'regime_registry' },
      { id: 'axis1_confidence_tier', label: 'Confidence tier (composite)',
        cats: ['T1','T2','T3','T4','SV_only','Layer_C_only','unknown'],
        legacy_alias: 'confidence_tier', composite: true },
    ]},
  { axis: 2, label: 'Structure', short: 'STR', color: '#8a94a3',
    sub: [
      { id: 'axis2_size_class', label: 'Size class',
        cats: ['small_lt1Mb','intermediate_1to5Mb','large_5to10Mb','supergene_gte10Mb'],
        derive_from: 'span_bp' },         // browser-computed
      { id: 'axis2_centromere_position', label: 'Centromere position',
        cats: ['paracentric','pericentric','unknown'] },
      { id: 'axis2_arm_position', label: 'Arm position',
        cats: ['sub_telomeric','interstitial','peri_centromeric','unknown'] },
      { id: 'axis2_topology', label: 'Topology',
        cats: ['simple_A','synteny_boundary_B','fusion_fission_C',
               'terminal_translocation_D','recurrent_hotspot_E','ambiguous_F'],
        // promote page-16b _msAutoSuggestArchitecture into the canonical schema
        derive_from: 'architecture_class_msAutoSuggest' },
      { id: 'axis2_boundary_quality', label: 'Boundary quality',
        cats: ['sharp','fuzzy','unknown'], legacy_alias: 'boundary_quality' },
    ]},
  { axis: 3, label: 'Mechanism', short: 'MECH', color: '#7ad3db',
    sub: [
      { id: 'axis3_breakpoint_substrate', label: 'Breakpoint substrate',
        cats: ['NAHR_TE','NAHR_segdup','NHEJ','FoSTeS_MMBIR','unknown'],
        legacy_alias: 'mechanism_class' },
      { id: 'axis3_te_enrichment', label: 'TE enrichment at breakpoints',
        cats: ['enriched_LTR','enriched_LINE','enriched_DNA_TE','enriched_tandem',
               'no_enrichment','unknown'],
        // populated from cs_breakpoints flank density - already in atlas
        source: 'cs_breakpoints_flank_density' },
      { id: 'axis3_recombination_class', label: 'Recombination/gene-flux',
        cats: ['fully_suppressed','partial_with_conversion','breakpoint_leaky',
               'gene_conversion','double_crossover','mixed','unknown'],
        legacy_alias: 'recombinant_class' },
      { id: 'axis3_dxy_profile_shape', label: 'dXY profile shape',
        cats: ['flat','linear_decay','U_shape','unknown'],
        // NEW - needs cluster-side computation (TIER-1 missing test #1)
        source: 'dxy_breakpoint_profile_v1' },
    ]},
  { axis: 4, label: 'History', short: 'HIST', color: '#3cc08a',
    sub: [
      { id: 'axis4_age_model', label: 'Age model',
        cats: ['YOUNG-POP','OLD-POLY','OLD-BP-YOUNG-INV',
               'LINEAGE-KARYO','MULTI-AGE-HOTSPOT','unknown'],
        // promote page-16b age engine to canonical
        derive_from: 'msAutoSuggestAgeModel' },
      { id: 'axis4_age_estimate_my', label: 'Age estimate (My, ABC posterior)',
        cats: null, type: 'numeric_with_ci',
        // NEW - TIER-1 missing test #3
        source: 'abc_age_v1' },
      { id: 'axis4_polarisation', label: 'Polarisation',
        cats: ['derived','ancestral','unresolved','unknown'],
        source: 'synteny_dollo' },        // already in registry
      { id: 'axis4_polymorphism_type', label: 'Polymorphism type',
        cats: ['type_I_transient','type_II_balanced','fixed_difference','unknown'],
        // Faria 2019 Type I/II - derived from age + frequency
        legacy_alias_partial: 'polymorphism_class' },
      { id: 'axis4_recurrence', label: 'Recurrence / cross-species sharing',
        cats: ['focal_only','sister_shared','clade_shared','recurrent_hotspot','unknown'],
        // populated from page-16b lineage_distribution + Twisst weight
        source: 'synteny_multispecies_v1' },
    ]},
  { axis: 5, label: 'Population dynamics', short: 'POP', color: '#f5a524',
    sub: [
      { id: 'axis5_frequency_class', label: 'Frequency class',
        cats: ['rare_lt5pct','intermediate_5to50pct','common_gt50pct','unknown'],
        derive_from: 'cohort_af' },        // browser-computed
      { id: 'axis5_hwe_state', label: 'HWE',
        cats: ['het_excess','het_deficit','equilibrium','unknown'] },
      { id: 'axis5_family_confound', label: 'Family / pedigree confound',
        cats: ['multi_family','few_family','single_family',
               'pca_family_confounded','cleared_via_ngsRelate','unknown'],
        legacy_alias: 'family_linkage' },
      { id: 'axis5_cline', label: 'Geographic / environmental cline',
        cats: ['present','absent','parallel','not_applicable_single_population','unknown'],
        // for current cohort = not_applicable_single_population
      },
      { id: 'axis5_substructure', label: 'Within-arrangement substructure',
        cats: ['monolithic','two_subclusters','three_or_more','unknown'],
        legacy_alias: 'internal_structure' },
    ]},
  { axis: 6, label: 'Function & Application', short: 'FUNC', color: '#e0555c',
    sub: [
      { id: 'axis6_cargo_class', label: 'Cargo class',
        cats: ['gene_rich','gene_poor','contains_immune_MHC',
               'contains_reproductive','contains_growth_metabolic','unknown'] },
      { id: 'axis6_go_enrichment', label: 'GO/KEGG enrichment (matched null)',
        cats: ['significant_immune','significant_metabolic','significant_other',
               'no_enrichment','unknown'],
        // NEW - TIER-1 missing test #4 (matched-null permutation)
        source: 'cargo_enrichment_v1' },
      { id: 'axis6_candidate_genes', label: 'Candidate breeding genes',
        cats: null, type: 'gene_list_short',
        source: 'candidate_gene_list' },
      { id: 'axis6_trait_association', label: 'Trait association (GLMM)',
        cats: ['significant','suggestive','none','no_phenotype_data','unknown'],
        source: 'trait_glmm_v1' },         // optional, only if phenotypes
      { id: 'axis6_load_class', label: 'Deleterious load (R_xy)',
        cats: ['enriched_HOMO_1','enriched_HOMO_2','symmetric','depleted','unknown'],
        legacy_alias: 'burden_class',
        // NEW - TIER-1 missing test #5 (R_xy with CI) refines this
        source: 'rxy_load_v1' },
      { id: 'axis6_marker_tier', label: 'Marker tier',
        cats: ['T1_private_indel','T2_multi_marker','T3_breakpoint_pcr','T4_exploratory','unknown'],
        // already in atlas on page 18
        source: 'marker_panel' },
    ]},
];
```

#### 3.1.3 Back-compat: do NOT delete `_TIER_AXES`

The R-pipeline writes `final_classification.<cid>.{existence_layer_a, ...}` per the **current** schema. Removing `_TIER_AXES` would break that contract.

**Solution**: keep `_TIER_AXES` in place as a **legacy alias map**. Add a thin translation function:

```javascript
// Translate a final_classification entry (legacy schema) -> axis-keyed view.
// Called once per render in renderCandidateTier(). The R-pipeline can
// continue to write the legacy schema unchanged; the scrubber renders the
// six-card layout from it.
function _translateLegacyToAxes(candEntry) {
  if (!candEntry) return null;
  const out = { axis1: {}, axis2: {}, axis3: {}, axis4: {}, axis5: {}, axis6: {} };
  for (const card of _AXIS_CARDS) {
    for (const sub of card.sub) {
      // Direct legacy alias
      if (sub.legacy_alias && candEntry[sub.legacy_alias] != null) {
        out[`axis${card.axis}`][sub.id] = candEntry[sub.legacy_alias];
        continue;
      }
      // New keys (axis-prefixed) - read directly if present
      if (candEntry[sub.id] != null) {
        out[`axis${card.axis}`][sub.id] = candEntry[sub.id];
        continue;
      }
      // Derived from atlas state (no R-pipeline contribution needed)
      if (sub.derive_from === 'span_bp')  out[`axis${card.axis}`][sub.id] = _classifySize(state.candidate);
      if (sub.derive_from === 'cohort_af') out[`axis${card.axis}`][sub.id] = _classifyAF(state.candidate);
      if (sub.derive_from === 'msAutoSuggestAgeModel') {
        const r = _msAutoSuggestAgeModel(/* current bp lookup */);
        out[`axis${card.axis}`][sub.id] = r ? r.age_model : 'unknown';
      }
      // No data - leave undefined; renderer shows "unknown" pill
    }
  }
  return out;
}
```

This means **the R-pipeline can keep emitting the existing schema** until it migrates voluntarily. The scrubber gracefully accepts both the legacy keys and the new axis-prefixed keys.

#### 3.1.4 New `_axisCardValueColor()` replaces `_tierAxisValueColor()`

The existing color logic (line 49150) is correct for the legacy axes. Extend with the new sub-axis vocabularies. No structural change - just more cases in the switch.

### 3.2 The empty-state UI on page 4

Currently page 4 says "Tier renders empty until the cluster-side `final_classification.json` layer is loaded".

After migration: **each card renders independently** based on what data is available. So even with zero `final_classification` loaded, the card UI still renders Axis 2 (size from `span_bp`, computed live), Axis 5 frequency (from cohort AF, computed live), and Axis 6 marker tier (from page 18 / `variant_afs.json` if loaded). The cards that depend on R-pipeline output show "data not yet available" placeholders.

This is a **UX improvement**, not a regression - the user can see partial classification immediately, and only the truly R-pipeline-dependent axes (most of Axis 1 and Axis 4) wait for the layer.

---

## 4. Catalogue (page 3) — add six axis-summary columns

### 4.1 New columns on `CAT_COLS_BASE` (line 58494)

Add six columns to the catalogue, each a **short categorical summary** read from the per-candidate `final_classification` entry:

```javascript
// Append to CAT_COLS_BASE, in this order, after 'verdict' and before 'diamond':
  { key: 'axis1_tier',     label: 'Evid (T)',  type: 'str',
    derive: c => _axisGet(c, 1, 'axis1_confidence_tier') || '-',
    tooltip: 'Axis 1 - Evidence tier composite. T1=Layer A+B+C+D pass + VALIDATED groups. T4=evidence weakest.' },
  { key: 'axis2_size',     label: 'Size',      type: 'str',
    derive: c => _classifySize(c),    // browser-computed from span_bp
    tooltip: 'Axis 2 - Structure size class. small <1 Mb / intermediate 1-5 Mb / large 5-10 Mb / supergene >=10 Mb.' },
  { key: 'axis3_mech',     label: 'Mech',      type: 'str',
    derive: c => _axisGet(c, 3, 'axis3_breakpoint_substrate') || '-',
    tooltip: 'Axis 3 - Mechanism: NAHR-TE / NAHR-segdup / NHEJ / FoSTeS-MMBIR / unknown.' },
  { key: 'axis4_age',      label: 'Age',       type: 'str',
    derive: c => _axisGet(c, 4, 'axis4_age_model') || '-',
    tooltip: 'Axis 4 - Age model: YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV / LINEAGE-KARYO / MULTI-AGE-HOTSPOT.' },
  { key: 'axis5_freq',     label: 'Freq',      type: 'str',
    derive: c => _classifyAF(c),
    tooltip: 'Axis 5 - Population frequency class: rare <5% / intermediate 5-50% / common >50%.' },
  { key: 'axis6_marker',   label: 'Marker',    type: 'str',
    derive: c => _axisGet(c, 6, 'axis6_marker_tier') || '-',
    tooltip: 'Axis 6 - Marker tier: T1 private indel / T2 multi-marker / T3 breakpoint PCR / T4 exploratory.' },
```

### 4.2 Catalogue export already supports new columns

Line 58576: `let cols = [id, ...CAT_COLS_BASE.slice(0, 8), conservationCol, ...CAT_COLS_BASE.slice(8)];`

The TSV/Markdown export iterates `getCatCols()`, so new columns are exported automatically. No additional code change needed for export.

### 4.3 Filter chips for the catalogue

Add 6 filter chips (one per axis) to the catalogue toolbar so users can filter on, e.g. "show all OLD-POLY · large · NAHR_TE" candidates. Implementation: extend the existing `catState.filterX` pattern (already used for `verdict`, `regimes`).

---

## 5. The six TIER-1 missing tests — atlas integration spec

Each of these needs (i) a JSON layer name with schema, (ii) a recogniser in `_identifyRegistryArtifact` (line 50157), (iii) a state slot, (iv) at least one render surface.

### 5.1 Breakpoint-restricted dXY + U-shape profile

- **Layer name:** `dxy_breakpoint_profile_v1.json` (extends existing `dxy_per_inversion_v1`).
- **Schema:**
  ```json
  {
    "format_version": "dxy_breakpoint_profile_v1",
    "candidate_id": "...",
    "windows": [{"window_bp_mid": 21500000, "dist_from_bp_kb": 100, "dxy": 0.0042, "in_breakpoint_band": true}],
    "summary": {
      "dxy_left_bp": 0.0048, "dxy_right_bp": 0.0046,
      "dxy_centre": 0.0021,
      "fold_centre_vs_bp": 0.45,
      "profile_shape": "U_shape",
      "profile_shape_confidence": "high"
    }
  }
  ```
- **State slot:** `state.data.dxy_breakpoint_profile[<cid>]`.
- **Render surfaces:**
  - Page 4 Axis 3 card -> `axis3_dxy_profile_shape` pill.
  - Page 6 popstats -> new track "dXY (per arrangement)" overlay.
  - Page 8 windows table -> new column `dxy_to_bp`.
  - Page 17 stats profile -> row in Axis 3 category (rename "Repeat context" to "Mechanism").
- **Loader hook:** add to `_identifyRegistryArtifact` filename regex `/dxy_breakpoint_profile_v\d+\.json$/`.

### 5.2 Twisst topology weighting per inversion

- **Layer name:** `twisst_per_inversion_v1.json`.
- **Schema:**
  ```json
  {
    "format_version": "twisst_per_inversion_v1",
    "candidate_id": "...",
    "taxa": ["HOMO_1", "HOMO_2", "Cmac"],
    "topologies": ["((HOMO_1,HOMO_2),Cmac)", "((HOMO_1,Cmac),HOMO_2)", "((HOMO_2,Cmac),HOMO_1)"],
    "windows": [{"start_bp": 21500000, "end_bp": 21600000, "weights": [0.78, 0.12, 0.10]}],
    "summary": {
      "monophyly_weight_mean": 0.78,
      "monophyly_weight_centre": 0.85,
      "monophyly_weight_breakpoints": 0.45,
      "monophyletic": true
    }
  }
  ```
- **Render surfaces:**
  - Page 4 Axis 4 card -> `axis4_recurrence` and a small inline ribbon showing topology weights along the inversion.
  - Page 16b -> augment the existing per-species ribbon view with the Twisst weight track.
  - Page 17 stats profile -> new row "Arrangement monophyly (Twisst)" in Axis 4 category.

### 5.3 ABC coalescent age estimate

- **Layer name:** `abc_age_v1.json` (per-candidate, flagship inversions only).
- **Schema:**
  ```json
  {
    "format_version": "abc_age_v1",
    "candidate_id": "...",
    "method": "fastsimcoal2",
    "model_specification": "...",
    "age_my_point": 0.82,
    "age_my_lo95": 0.42, "age_my_hi95": 1.51,
    "age_generations_point": 4100000,
    "additional_params": {"Ne_inverted": 12000, "Ne_standard": 28000, "migration_rate": null},
    "fit_diagnostic": {"posterior_log_prob": -1234.5, "n_simulations": 100000}
  }
  ```
- **Render surfaces:**
  - Page 4 Axis 4 card -> `axis4_age_estimate_my` numeric pill with CI bracket.
  - Page 17 stats profile -> new row "Inversion age (ABC)" in Axis 4 category.
  - Bundle export -> carries through.

### 5.4 Matched-null GO + cargo enrichment

- **Layer name:** `cargo_enrichment_v1.json`.
- **Schema:**
  ```json
  {
    "format_version": "cargo_enrichment_v1",
    "candidate_id": "...",
    "null_model": "size_and_recombination_class_matched",
    "n_null_replicates": 1000,
    "gene_density": {"observed": 5.2, "null_mean": 4.8, "null_p": 0.32},
    "go_enrichment": [
      {"go_id": "GO:0006955", "go_name": "immune response", "n_observed": 4, "fdr": 0.003, "passes_matched_null": true}
    ],
    "kegg_enrichment": [],
    "candidate_breeding_genes": ["tlr5", "il10", "cyp19a1"],
    "candidate_gene_categories": ["immunity", "reproduction"]
  }
  ```
- **Render surfaces:**
  - Page 4 Axis 6 card -> `axis6_go_enrichment` and `axis6_candidate_genes`.
  - Page 17 stats profile -> already has `go_kegg_enrichment` row but currently null-derive - populate from this layer.

### 5.5 R_xy formal load ratio with CI

- **Layer name:** `rxy_load_v1.json`.
- **Schema:**
  ```json
  {
    "format_version": "rxy_load_v1",
    "candidate_id": "...",
    "predictor": "VESM_650M",
    "groups": ["HOMO_1", "HOMO_2"],
    "rxy_point": 1.34,
    "rxy_lo95": 1.08, "rxy_hi95": 1.62,
    "n_jackknife": 50,
    "interpretation": "HOMO_1 carries higher deleterious load",
    "n_deleterious_HOMO_1": 142, "n_deleterious_HOMO_2": 98,
    "n_neutral_HOMO_1": 1240, "n_neutral_HOMO_2": 1180
  }
  ```
- **Render surfaces:**
  - Page 4 Axis 6 card -> `axis6_load_class` pill with CI tooltip.
  - Page 17 stats profile -> existing `deleterious_burden` row populates from this.

### 5.6 pyrho recombination-rate landscape

- **Layer name:** `pyrho_rho_v1.json` (chromosome-wide).
- **Schema:**
  ```json
  {
    "format_version": "pyrho_rho_v1",
    "chrom": "C_gar_LG28",
    "windows": [{"start_bp": 0, "end_bp": 50000, "rho_per_bp": 0.0023}],
    "summary_per_candidate": {"<cid>": {"rho_mean_inside": 0.00018, "rho_mean_flank": 0.00420, "fold_reduction": 23.3, "trough_present": true}}
  }
  ```
- **Render surfaces:**
  - Page 6 popstats -> new track "rho (population recombination rate)".
  - Page 4 Axis 3 card -> small visual chip "rho trough confirmed (23x)".
  - Page 17 stats profile -> new row "Recombination rate (pyrho)" in Axis 3 category.

---

## 6. Page 16b promotion — make it the canonical Axis 4 engine

### 6.1 Currently

Page 16b (multi-species cockpit) is a **local cockpit** that displays per-breakpoint architecture/age classification but does **NOT** write to `final_classification`. Its output is exported as a TSV (line 22090) - a one-off file users download manually.

### 6.2 Recommended change

When page 16b's `_msAutoSuggestArchitecture` and `_msAutoSuggestAgeModel` produce a result with confidence >= medium, **mirror that result into `state.data.final_classification[<cid>]`** under the new axis-keyed schema:

```javascript
function _msPushClassificationToAxes(bp, archResult, ageResult) {
  // bp.id maps to candidate_id when bp is a polymorphic-cohort breakpoint
  // (the star flag in catalogue means cs_breakpoint coincides with cohort cand).
  const cid = _msResolveCandidateId(bp);
  if (!cid) return;
  if (!state.data.final_classification) state.data.final_classification = {};
  if (!state.data.final_classification[cid]) state.data.final_classification[cid] = {};
  const fc = state.data.final_classification[cid];
  // Browser-derived axis values are written under axis-prefixed keys to
  // distinguish them from R-pipeline-emitted legacy keys.
  if (archResult.class && archResult.confidence !== 'unknown') {
    fc.axis2_topology = _archClassToCanonical(archResult.class);  // 'simple_A' etc
    fc.axis2_topology_source = 'browser_msAutoSuggest';
    fc.axis2_topology_confidence = archResult.confidence;
  }
  if (ageResult.age_model) {
    fc.axis4_age_model = ageResult.age_model;
    fc.axis4_age_model_source = 'browser_msAutoSuggest';
    fc.axis4_age_model_confidence = ageResult.confidence;
    fc.axis4_age_model_signals = ageResult.signals;
  }
  // Trigger page 4 re-render if active candidate matches
  if (state.candidate && state.candidate.id === cid && state.activePage === 'page4') {
    if (typeof renderCandidateTier === 'function') renderCandidateTier();
  }
}
```

This keeps page 16b's user-facing flow unchanged but makes its output visible on page 4 (the central tier view) and exports cleanly.

### 6.3 Manual override propagation

The `_msGetClassification(bp.id)` override (line 22115) is already user-curated. Mirror those overrides into `final_classification` too with `_source: 'manual'`.

---

## 7. Page 17 (stats profile) — rename categories to axes

### 7.1 Current categories (line 23475)

- Breakpoint architecture
- Genomic composition
- Functional cargo
- Population variation
- Breeding utility

### 7.2 Target six-axis categories

Rename + reassign rows. **Existing rows do NOT change semantically**, only category-tag is updated:

| Old category | New axis | Rows that move |
|---|---|---|
| Breakpoint architecture | **Axis 3 - Mechanism** | breakpoints_near_synteny_edges -> 3.3 recurrence; fusion_fission_proximity -> 3.3 |
| Genomic composition | **Axis 2 - Structure** + **Axis 3** | repeat_density_flanks -> Axis 3.2; telomere_centromere_proximity -> Axis 2.3; gene_density_inversions -> Axis 2 (descriptor) |
| Functional cargo | **Axis 6 - Function & Application** | go_kegg_enrichment, breeding_candidate_genes |
| Population variation | **Axis 5 + Axis 6** | heterozygosity_inversions -> Axis 5; roh_overlap -> Axis 5; deleterious_burden -> Axis 6.1; haplotype_differentiation -> Axis 5 (frequency-related) or Axis 4 (history-related, fold_elevation) |
| Breeding utility | **Axis 6** | markerability |

### 7.3 New rows to ADD per the 6 missing tests

```javascript
// Append to SP_DEFAULT_ROWS (line 23475):
  { id: 'dxy_profile_shape', category: 'Axis 3 - Mechanism',
    statistic: 'dXY profile shape (U-shape / flat)',
    test: 'breakpoint-distance regression + AIC',
    derive_from: 'dxy_breakpoint_profile_v1',
    interpretation_default: 'U-shape = old + gene flux; flat = young or no flux.' },
  { id: 'twisst_monophyly', category: 'Axis 4 - History',
    statistic: 'Arrangement monophyly weight (Twisst)',
    test: 'topology weighting',
    derive_from: 'twisst_per_inversion_v1',
    interpretation_default: 'Higher = arrangements more cleanly monophyletic = older or stronger recombination suppression.' },
  { id: 'abc_age', category: 'Axis 4 - History',
    statistic: 'Inversion age (ABC posterior, My)',
    test: 'fastsimcoal2 + ABC',
    derive_from: 'abc_age_v1',
    interpretation_default: 'Posterior on inversion age. Wide CI under model uncertainty.' },
  { id: 'pyrho_trough', category: 'Axis 3 - Mechanism',
    statistic: 'Recombination-rate trough (pyrho rho)',
    test: 'fold-reduction inside vs flank',
    derive_from: 'pyrho_rho_v1',
    interpretation_default: 'Strong trough confirms recombination suppression in heterokaryotypes.' },
  { id: 'rxy_load_ratio', category: 'Axis 6 - Function & Application',
    statistic: 'Deleterious load R_xy ratio (HOMO_1 vs HOMO_2)',
    test: 'Do 2015 R_xy with jackknife CI',
    derive_from: 'rxy_load_v1',
    interpretation_default: 'R_xy > 1 -> HOMO_1 carries higher deleterious load; CI excludes 1 = significant.' },
  { id: 'cargo_matched_null', category: 'Axis 6 - Function & Application',
    statistic: 'GO/KEGG enrichment vs matched null',
    test: 'permutation against size+recombination-matched intervals',
    derive_from: 'cargo_enrichment_v1',
    interpretation_default: 'Enrichment robust to size/recombination-class confound.' },
```

### 7.4 Filter pillchip toolbar update

Page 17 currently has `filter_category: 'all' | <category>`. Update the chip set to: `all / Axis 1 / Axis 2 / Axis 3 / Axis 4 / Axis 5 / Axis 6`.

---

## 8. Bundle export — `atlas_candidate_export_v3`

### 8.1 Current export (line 47097, `_EXPORT_FORMAT_VERSION = 'atlas_candidate_export_v2'`)

Carries: candidate coordinates, per-sample band assignments, haplotype labels, atlas analytics (diamond, sigma, inheritance), atlas constants. Does **not** carry classification.

### 8.2 v3 schema additions

```javascript
const _EXPORT_FORMAT_VERSION = 'atlas_candidate_export_v3';

// In buildCandidateExportRecord(c), append after sigmaSummary (line ~47290):

  // Six-axis classification block (NEW in v3). Auto-extracts from
  // state.data.final_classification[cid] using _translateLegacyToAxes and
  // adds browser-derived axis values (size, frequency, age model from
  // page-16b auto-suggest).
  let classification = null;
  if (state.data && state.data.final_classification) {
    const cid = c.id || (c.ref_l2 != null ? `ref_l2_${c.ref_l2}` : null);
    const candEntry = (cid != null) ?
      (state.data.final_classification[cid] || state.data.final_classification[c.id] || null) : null;
    if (candEntry) {
      classification = _translateLegacyToAxes(candEntry);
      classification._format_version = 'six_axis_v1';
      classification._extracted_at = new Date().toISOString();
    }
  }
  // Page-16b auto-suggest mirror (always available if cs_breakpoints loaded)
  let auto_suggest = null;
  if (typeof _msAutoSuggestArchitecture === 'function' && state.crossSpecies) {
    const bp = _msResolveBpForCandidate(c);
    if (bp) {
      const arch = _msAutoSuggestArchitecture(bp, _msGetLineageDistribution(bp));
      const age  = _msAutoSuggestAgeModel(bp,
                     _msGetLineageDistribution(bp),
                     _msGetDxyForBreakpoint(bp),
                     _msPolarizeKaryotypeEvent(bp.gar_chr));
      auto_suggest = { architecture_class: arch.class, architecture_confidence: arch.confidence,
                       age_model: age.age_model, age_confidence: age.confidence,
                       age_signals: age.signals };
    }
  }
  // Marker tier (page 18) - mirror if computed
  let marker_summary = null;
  if (typeof _mpResolveTierForCandidate === 'function') {
    marker_summary = _mpResolveTierForCandidate(c) || null;
  }

  return {
    // ... existing fields unchanged ...
    classification: classification,             // NEW v3
    auto_suggest:   auto_suggest,               // NEW v3
    marker_summary: marker_summary,             // NEW v3
  };
```

### 8.3 Top-level metadata

```javascript
// In buildAtlasCandidateExport(opts), add to top-level export object:
  classification_summary: {
    schema:        'six_axis_v1',
    n_classified:  /* count of candidates with non-null classification block */,
    axis_coverage: {
      axis1: /* fraction of candidates with axis1 sub-axes populated */,
      axis2: 0, axis3: 0, axis4: 0, axis5: 0, axis6: 0,
    },
  },
  species_metadata: {
    focal_species_label:  state.species && state.species.focal_species_label  || 'Clarias gariepinus',
    focal_species_short:  state.species && state.species.focal_species_short  || 'Cgar',
    outgroup_species_label: state.species && state.species.outgroup_species_label || 'Clarias macrocephalus',
    outgroup_species_short: state.species && state.species.outgroup_species_short || 'Cmac',
  },
```

### 8.4 The bundle is now species-portable

For *C. macrocephalus* future paper: change `state.species` -> re-export -> all per-candidate panels regenerate with correct species labels. No code change.

---

## 9. New optional layer: `species_metadata.json`

Tiny layer that lets the user swap species without touching code.

```json
{
  "format_version": "species_metadata_v1",
  "focal_species_label": "Clarias gariepinus",
  "focal_species_short": "Cgar",
  "focal_species_common_name": "African catfish",
  "outgroup_species_label": "Clarias macrocephalus",
  "outgroup_species_short": "Cmac",
  "outgroup_species_common_name": "bighead catfish",
  "additional_outgroups": [
    {"label": "...", "short": "..."}
  ],
  "phenotype_available": true,
  "phenotypes_list": ["body_weight_g", "growth_rate_g_per_day", "survival_days"],
  "cohort_design": "hatchery_pedigreed",
  "n_samples": 226,
  "coverage_x": 9,
  "reference_genome": "fClaHyb_Gar_LG.fa"
}
```

Loader recogniser pattern: `/species_metadata_v\d+\.json$/`. Stored at `state.species`. Read by all page-render functions that currently hardcode "C. gariepinus" / "Cgar".

Pages that need to switch from hardcoded to `state.species` lookups:
- Page 13 cross-species header (line 6982 onwards)
- Page 16b multi-species cockpit
- Page 17 stats profile (african_hint / bighead_hint columns currently hard-name "African" / "Bighead")
- Page 18 marker readiness (the protective methods sentence)

---

## 10. Help page (page 5) update

Lines 7204–7400 are the help reference. Current help references "14-axis classification" in 3 places. Update to:
- Replace "14-axis classification view" -> "six-axis classification cards (Evidence / Structure / Mechanism / History / Population / Function & Application)"
- Add a new help-table row explaining the six axes mapping to manuscript Result blocks.
- Add a glossary entry for `axis1_*` ... `axis6_*` keys.

Also: the Q1–Q7 documentation in §13 of the help stays as-is (it's the IMPLEMENTATION layer); just add a paragraph noting "the six axes are the OUTPUT view; Q1–Q7 are the question-driven workflow that produces the data".

---

## 11. Migration order — sequenced so nothing breaks mid-flight

Each step is independently shippable, leaves the atlas working, and progresses the migration.

### Step 1 — Mechanical (no behaviour change, ~1 day)
- Add `_AXIS_CARDS` constant alongside `_TIER_AXES` (do not remove `_TIER_AXES`).
- Add `_translateLegacyToAxes()` helper. Test on existing fixture data — output must round-trip.
- Add `_classifySize()` and `_classifyAF()` browser-side helpers.
- Add `state.species` default object with current Cgar/Cmac values. No reads from it yet.

### Step 2 — Page 4 dual-render (~1–2 days, opt-in)
- Add a new toggle on page 4 sub-view bar: `Karyotype | Tier (legacy) | Axis cards (new)`.
- Wire `Axis cards` mode to `_AXIS_CARDS` + `_translateLegacyToAxes`.
- Default = `Tier (legacy)` for back-compat; users opt in via the new button.
- Test that data flows: load existing `final_classification.json` -> both Tier and Axis cards show same content.

### Step 3 — Page 16b auto-mirror (~1 day)
- Implement `_msPushClassificationToAxes`.
- Hook into the page-16b auto-suggest result paths (line 22110+).
- Page 4 Axis 4 card now populates from page-16b auto-suggest even with no R-pipeline output.

### Step 4 — Catalogue axis columns (~1 day)
- Add the 6 axis-summary columns to `CAT_COLS_BASE`.
- Catalogue export TSV picks them up automatically.

### Step 5 — Stats profile rename + new rows (~1 day)
- Rename `category` field on existing `SP_DEFAULT_ROWS` to use Axis 1–6 vocabulary.
- Add the 6 new rows for the missing TIER-1 tests.
- Update filter chip labels.
- Old user-supplied `stats_profile.json` overlays continue to work (the `category` field is just a label).

### Step 6 — Loader recognisers for the 6 new layers (~1–2 days)
- Add filename regex matches in `_identifyRegistryArtifact`.
- Add light schema validators (required-field checks) per the schemas in §5.
- Add state slots: `state.data.dxy_breakpoint_profile`, `state.data.twisst_per_inversion`, `state.data.abc_age`, `state.data.cargo_enrichment`, `state.data.rxy_load`, `state.data.pyrho_rho`.
- Layers do nothing until the corresponding axis-card sub-axis reads them — graceful degradation.

### Step 7 — Bundle export v3 (~1 day)
- Bump `_EXPORT_FORMAT_VERSION`.
- Add `classification`, `auto_suggest`, `marker_summary` to per-candidate record.
- Add `classification_summary` and `species_metadata` to top-level.
- Old v2 readers continue to work because v3 is additive.

### Step 8 — Default to Axis cards (~1 day)
- After 2–3 weeks of dual-render testing, flip the default mode on page 4 from `Tier (legacy)` to `Axis cards`.
- Keep `Tier (legacy)` as opt-in for one further release cycle.

### Step 9 — Retire (~later, after manuscript submission)
- Once R-pipeline emits axis-keyed schema natively, remove `_TIER_AXES` and the legacy tier render.
- Keep `_translateLegacyToAxes` as a back-compat reader for any user with old `final_classification.json` files.

---

## 12. Per-tab change summary

| Tab | Change required | Effort |
|---|---|---|
| 1 local PCA \|z\| (page1) | None | — |
| 2 local PCA θπ (page12) | None | — |
| 2b local PCA GHSL (page15) | None | — |
| 3 candidate focus (page2) | Update the 18-block chip row tooltip to mention axis derivation | trivial |
| 4 boundaries (page11) | None | — |
| 5 catalogue (page3) | Add 6 axis-summary columns (§4) | small |
| 6 negative regions (page19) | None — already serves Axis 1 inverted | — |
| 7 karyotype/tier (page4) | **MAIN refactor** — six-card layout (§3) | medium |
| 8 popstats (page6) | Add 2 new tracks (dXY profile, ρ/pyrho) (§5.1, §5.6) | small |
| 9 ancestry (page7) | None | — |
| 10 windows (page8) | Add `dxy_to_bp` column (§5.1) | trivial |
| 11 confirmed (page9) | None | — |
| 11b annotation (page21) | None | — |
| 12 markers (page10) | None — legacy phase-13 layer, leave as-is | — |
| 13 cross-species (page16) | Add Axis 4 recurrence chip per breakpoint | small |
| 13b multi-species (page16b) | **Promote to canonical engine** — auto-mirror to final_classification (§6) | medium |
| 14 stats profile (page17) | Rename categories, add 6 new rows (§7) | small |
| 15 marker panel (page18) | Add header tag "Axis 6 · Application", no other change | trivial |
| 15b overview (page_overview) | Inherit catalogue axis columns automatically | — |
| 16 help (page5) | Update help text to reflect six-axis (§10) | small |

**Total estimated effort:** 7–10 days of focused work, sequenced as in §11.

---

## 13. The integration matrix — each test ↔ each tab

This is the single table to keep open while building. Each cell tells you whether tab consumes test output, derives data for it, or shows it.

| Tab \\ Layer | dxy_profile | twisst | abc_age | cargo_enrich | rxy_load | pyrho_rho | species_meta |
|---|---|---|---|---|---|---|---|
| 5 catalogue (axis cols) | derives Axis 3 col | — | derives Axis 4 col | derives Axis 6 col | derives Axis 6 col | — | reads species |
| 6 popstats | new track | — | — | — | — | new track | — |
| 7 karyotype/tier | Axis 3 card | Axis 4 card | Axis 4 card | Axis 6 card | Axis 6 card | Axis 3 card | reads species |
| 8 windows | new column | — | — | — | — | — | — |
| 13 cross-species | — | — | — | — | — | — | reads species |
| 13b multi-species | feeds Axis 4 | feeds Axis 4 | feeds Axis 4 | — | — | — | reads species |
| 14 stats profile | new row | new row | new row | new row | new row | new row | reads species |
| 15 marker panel | — | — | — | — | — | — | reads species |
| Bundle export v3 | classification.axis3 | classification.axis4 | classification.axis4 | classification.axis6 | classification.axis6 | classification.axis3 | top-level |

---

## 14. What the cluster-side R pipeline must emit (post-migration target)

This is the contract the cluster-side R pipeline (`characterize_candidate.R` + `classify_inversions.R`) commits to. **Today, it emits the legacy 14-axis schema.** **Tomorrow, it should emit the axis-keyed schema directly.** The translation step in §3.1.3 handles the transition.

```json
{
  "format_version": "final_classification_six_axis_v1",
  "<candidate_id>": {
    "axis1_layer_a": "pass",
    "axis1_layer_b": "pass",
    "axis1_layer_c": "pass",
    "axis1_layer_d": "fail",
    "axis1_cluster_quality": "VALIDATED",
    "axis1_regime_topology": "one_axis_3band",
    "axis1_confidence_tier": "T1",

    "axis2_size_class": "intermediate_1to5Mb",
    "axis2_centromere_position": "paracentric",
    "axis2_arm_position": "sub_telomeric",
    "axis2_topology": "simple_A",
    "axis2_boundary_quality": "sharp",

    "axis3_breakpoint_substrate": "NAHR_TE",
    "axis3_te_enrichment": "enriched_LTR",
    "axis3_recombination_class": "fully_suppressed",
    "axis3_dxy_profile_shape": "U_shape",

    "axis4_age_model": "OLD-POLY",
    "axis4_age_estimate_my": 0.82,
    "axis4_age_estimate_my_lo95": 0.42,
    "axis4_age_estimate_my_hi95": 1.51,
    "axis4_polarisation": "derived",
    "axis4_polymorphism_type": "type_II_balanced",
    "axis4_recurrence": "sister_shared",

    "axis5_frequency_class": "intermediate_5to50pct",
    "axis5_hwe_state": "het_excess",
    "axis5_family_confound": "multi_family",
    "axis5_cline": "not_applicable_single_population",
    "axis5_substructure": "two_subclusters",

    "axis6_cargo_class": "contains_immune_MHC",
    "axis6_go_enrichment": "significant_immune",
    "axis6_candidate_genes": ["tlr5", "il10"],
    "axis6_trait_association": "no_phenotype_data",
    "axis6_load_class": "enriched_HOMO_1",
    "axis6_marker_tier": "T1_private_indel",

    "_legacy_aliases_emitted": false,
    "_provenance": {
      "characterize_candidate_version": "vX.Y",
      "classify_inversions_version": "vX.Y",
      "emitted_at": "2026-05-03T00:00:00Z"
    }
  }
}
```

Compatibility note: if the R pipeline still emits legacy keys, the scrubber's `_translateLegacyToAxes` reads them. R-pipeline migration is decoupled from atlas migration — they can be done in either order.

---

## 15. Final recommendation summary

The atlas is **closer to the six-axis framework than it looks** — most of the heavy lifting is already done in code, just classified under different vocabulary. The migration is mostly **renaming and re-grouping**, with two genuinely new things:

1. **The six-card layout on page 4** (replacing the 14-axis grid).
2. **Six new optional JSON layers** for the TIER-1 missing tests (one each for breakpoint-restricted dXY, Twisst, ABC age, matched-null cargo, R_xy load, pyrho ρ).

Everything else is data-flow plumbing — back-compat aliases, browser-derived sub-axes that don't need any cluster-side change, and a bundle export bump.

**Lowest-risk first step:** Add `_AXIS_CARDS` and `_translateLegacyToAxes` (Step 1 above). The atlas keeps rendering exactly as today, and you have the migration scaffolding in place. Everything afterwards is incremental.

**Highest-value first user-visible step:** Page 4 Axis-cards mode behind a toggle (Step 2). Reviewers and collaborators can see the new layout immediately without any cluster-side change because page-16b auto-suggest + browser-derived size/frequency already populate Axes 2, 4, 5 in part.

**The species portability win is free**: once `state.species` is consulted instead of hardcoded "Cgar"/"Cmac" strings (Step 1, also), the same atlas serves the *C. macrocephalus* paper with a one-line config change.

End of audit.
