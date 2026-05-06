# SPEC — Hypothesis registry, multi-species support, classification UI

**Status**: forward-looking spec. None of this is implemented. Scoped end of
turn 113-tail (2026-05-02) to provide a working blueprint for whichever
fresh chat picks up the work.

**This spec is one of five**. See `SPEC_OVERVIEW_multispecies_architecture.md`
for how it composes with:
- `SPEC_comparative_te_breakpoint_fragility.md` — TE density across species
- `SPEC_phylogenetic_tree_integration.md` — branch placement of events
- `SPEC_busco_anchors.md` — deep-divergence protein anchors
This document focuses on the registry/UI side; layer-specific schemas live
in their own specs.

**Scope owners** (the people / sessions this spec is for):
- Future Claude sessions implementing the hypothesis registry
- Future Quentin reviewing scope before authorizing work
- Reviewers looking for what the atlas could become beyond this manuscript

**Reading order**:
1. This spec (forward-looking)
2. `HANDOFF_2026-05-01_session_summary.md` (what's already done)
3. `ATLAS_REFERENCE_for_phase4b_synthesis.md` (atlas internals reference)
4. Atlas page 5 — "Breakpoint classification framework" section (the schema
   this spec operationalizes)

---

## Why this spec exists

The atlas currently does **discovery** (find candidates, assign K-band labels)
and **display** (cross-species breakpoints, candidate focus, contingency).
What it does NOT yet do:

1. **Classify** each cs-breakpoint or candidate against the two-layer schema
   (architecture A–F, age YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV /
   LINEAGE-KARYO / MULTI-AGE-HOTSPOT)
2. **Score evidence** for each classification (which evidence layers
   support / refute the assigned class?)
3. **Track hypotheses** as first-class objects (Model A vs B vs C for the
   LG28 case, with predictions per evidence layer)
4. **Consume new layer types** that don't yet exist: 10-species synteny,
   miniprot anchors, dXY age signal

This spec describes a coherent design that addresses all four together,
because they're the same architectural problem: the atlas needs a
**hypothesis-and-evidence object model**.

---

## Sequencing — multi-chat plan

| Chat | Scope | Effort | Blocks |
|------|-------|--------|--------|
| 115a | Hypothesis registry data model + minimal page-14 UI | 2-3 hr | none |
| 115b | Per-breakpoint classification UI on page 16 | 2 hr | 115a |
| 115c | dXY age-signal evidence layer (`dxy_per_inversion_v1`) | 2 hr | 115a |
| 115d | 9-11 species synteny layer (`synteny_multispecies_v1`) | 3 hr | toolkit run on real data |
| 115e | Miniprot anchors layer (`miniprot_anchors_v1`) | 2 hr | miniprot pipeline output |
| 115f | Auto-suggest classification from evidence layers | 2-3 hr | 115b through 115e |
| 115g | Hypothesis-scoring engine + manuscript table export | 2 hr | all above |

Each chat ships independently. Earlier chats produce useful artifacts even if
later chats never happen. **Critical prerequisite for 115d**: the
catfish-synteny-toolkit (~3700 lines, shipped earlier as tarball) has never
been run on real LANTA data. Its output schema is theoretical until then.
Don't let chat 115d start until Quentin has produced one real example output
JSON. Otherwise the layer reader is built against a guess and rewritten later.

---

## Data model

### `Hypothesis` object

```jsonc
{
  "hypothesis_id": "lg28_model_b_breakpoint_reuse",
  "name": "LG28: recurrent breakpoint reuse",
  "applies_to": {
    "kind": "candidate",                  // "candidate" | "cs_breakpoint" | "region"
    "id": "cand_lg28_15_18",              // FK to candidate or cs_breakpoint
  },
  "architecture_class": "E",              // A..F per page-5 schema
  "age_model": "OLD-BP/YOUNG-INV",        // 5 tags per page-5 schema
  "confidence": "medium",                 // low / medium / high
  "predictions": [
    {
      "layer": "te_density",
      "expected": "high_all_te_at_breakpoint",
      "weight": 1.0,
      "evidence_status": "supports"       // supports | neutral | refutes | not_yet_tested
    },
    {
      "layer": "synteny_multispecies",
      "expected": "breakpoint_appears_in_multiple_lineages",
      "weight": 2.0,
      "evidence_status": "not_yet_tested"
    },
    {
      "layer": "dxy_per_inversion",
      "expected": "dxy_inside_inversion_NOT_extremely_old",
      "weight": 1.0,
      "evidence_status": "not_yet_tested"
    }
  ],
  "evidence_score": 1.0,                  // sum of weight * (supports=1, refutes=-1, neutral=0)
  "evidence_score_max": 4.0,              // sum of weights (max possible if all supports)
  "evidence_score_normalized": 0.25,      // score / score_max in [-1, 1]
  "manuscript_text": "...",               // free-text; populated from page-5 framing
  "annotator": "quentin",
  "created_at": "2026-05-02T12:00Z",
  "updated_at": "2026-05-02T12:00Z"
}
```

### `Classification` object — simpler, per-breakpoint

```jsonc
{
  "applies_to": {"kind": "cs_breakpoint", "id": "cs_bp_0001"},
  "architecture_class": "C",
  "age_model": "LINEAGE-KARYO",
  "confidence": "high",
  "annotator_notes": "Cgar LG28 fission/fusion onto Cmac LG01; clean break.",
  "annotator": "quentin",
  "updated_at": "2026-05-02T12:00Z"
}
```

A classification is a lightweight per-object record. A hypothesis is a richer
object that bundles **multiple** classifications + predictions + evidence
scoring. Most users will only need classifications; hypotheses are for the
manuscript / methods comparisons.

### Storage

- localStorage keys:
  - `inversion_atlas.classifications` — array of Classification objects
  - `inversion_atlas.hypotheses` — array of Hypothesis objects
- Export / import: same JSON shape, downloaded via "Export classifications"
  / "Export hypotheses" buttons on page 14
- IndexedDB: not needed at this scale (cohort has < 200 candidates × < 5
  hypotheses each = bounded)
- Per-cohort scoping: classifications are loaded/saved alongside the precomp
  JSON. Switching to a different cohort loads that cohort's classifications.

---

## Chat 115a — Hypothesis registry data model + page-14 UI

### Deliverables

1. State slots: `state.classifications[]`, `state.hypotheses[]`
2. Helpers: `_addClassification(c)`, `_updateClassification(id, patch)`,
   `_deleteClassification(id)`, `_addHypothesis(h)`, etc.
3. Persistence: `_persistClassifications()` + `_restoreClassifications()`
   (mirrors existing `_persistCrossSpecies` pattern)
4. Page 14 (currently "stats profile") gets a new collapsible section
   **"Hypothesis registry"** above the existing stats-profile rows.
   Shows: list of hypotheses, per-hypothesis evidence score bar, "Add new
   hypothesis" button.
5. **Pre-registered defaults for the LG28 case** — three hypotheses
   (Models A/B/C from the breakpoint-reuse framing), populated on first
   load:
   - `lg28_model_a_terminal_translocation` — Class D, MULTI-AGE-HOTSPOT
   - `lg28_model_b_breakpoint_reuse` — Class E, OLD-BP/YOUNG-INV
   - `lg28_model_c_ancestral_attachment_then_fission` — Class C, LINEAGE-KARYO

### Tests

- Add / update / delete / list hypotheses
- Persistence round-trip (save → reload → matches)
- Pre-registered defaults appear on fresh state
- LocalStorage quota guard (gracefully fail if quota exceeded)

### NOT in scope this chat

- Auto-suggest classification (chat 115f)
- Evidence scoring against real layer data (chat 115g)
- Manuscript table export (chat 115g)

---

## Chat 115b — Per-breakpoint classification UI on page 16

### Deliverables

1. In the page-16 focus panel (`_renderCrossSpeciesFocus`), add a
   classification editor below the breakpoint header:
   ```
   [Architecture] [A ▾] [Age model] [YOUNG-POP ▾] [Confidence] [low ▾]
   [Notes: free-text textarea]
   [Save]  [Reset]
   ```
2. Architecture and age dropdowns populated from page-5 schema (six A–F
   options + five age tags).
3. On Save: write a `Classification` object to `state.classifications[]`
   keyed by breakpoint ID.
4. Display the saved classification in the catalogue list (left column)
   as a small chip — `[D · OLD-BP]` next to the breakpoint name.

### Tests

- Classification persists across atlas reload
- Multiple breakpoints can carry independent classifications
- Catalogue chips render when classification exists
- Reset clears the editor without deleting saved classification

### Visual style

- Match existing `.cs-chip` style for the catalogue chip
- Architecture-class color coding: A=blue, B=teal, C=orange, D=purple,
  E=red, F=grey (matches page 5 emphasis)

---

## Chat 115c — dXY age-signal evidence layer

### Why this layer matters

Age model classification (Layer 2) needs a quantitative signal for
"is this inversion old or young?" The standard signal is dXY between
standard and inverted haplotypes inside the inversion (high = old, low =
young). Currently the atlas has no dXY layer. Adding one unlocks the
YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV distinction.

### Layer schema

```jsonc
// dxy_per_inversion_v1.json
{
  "tool": "dxy_per_inversion_v1",
  "schema_version": 1,
  "generated_at": "...",
  "params": { "min_carriers_per_class": 5, "windowed": false },
  "per_inversion": [
    {
      "candidate_id": "cand_lg28_15_18",
      "n_hom_ref": 60, "n_het": 106, "n_hom_inv": 60,
      "dxy_within_inversion_ref_vs_inv": 0.0042,
      "dxy_flank_left_ref_vs_inv": 0.0019,
      "dxy_flank_right_ref_vs_inv": 0.0021,
      "fold_elevation_inside_vs_flank": 2.10,
      "interpretation_default": "elevated_consistent_with_old_inversion"
    }
  ]
}
```

### Atlas-side

- Layer key: `dxy_per_inversion`
- Reader: parse and store under `state.layers.dxy_per_inversion`
- Display: per-candidate row on page 7 (karyotype/tier) gets a new column
- Hypothesis hook: hypotheses with `predictions[].layer === "dxy_per_inversion"`
  evaluate themselves against the loaded data; evidence_status changes
  from `not_yet_tested` to `supports` / `refutes` based on threshold
- Default thresholds: `fold_elevation > 1.5` ⇒ supports OLD-POLY;
  `fold_elevation < 1.2` ⇒ supports YOUNG-POP;
  `1.2 ≤ fold_elevation ≤ 1.5` ⇒ neutral

### NOT in scope this chat

- Computing dXY (that's R/Python pipeline work, not atlas)
- Per-window dXY (only per-inversion summary)

---

## Chat 115d — 9-11 species synteny layer

### Prerequisite

Quentin runs `catfish-synteny-toolkit` on real LANTA data and produces one
example `synteny_multispecies_v1.json`. Until then, this chat doesn't start.

### Provisional layer schema (subject to revision after first real run)

```jsonc
{
  "tool": "catfish_synteny_toolkit_v1",
  "schema_version": 1,
  "species": [
    { "name": "Clarias gariepinus", "haplotype": "fClaHyb_Gar_LG", "role": "focal" },
    { "name": "Clarias macrocephalus", "haplotype": "fClaHyb_Mac_LG", "role": "sister" },
    { "name": "Heterobranchus longifilis", "role": "near_outgroup" },
    { "name": "Trichomycterus rosablanca", "role": "deep_outgroup" }
    // ... 7-11 total
  ],
  "synteny_blocks": [
    {
      "block_id": "blk_0001",
      "members": [
        { "species": "Cgar", "chrom": "C_gar_LG28", "start_bp": 14800000, "end_bp": 18200000, "strand": "+" },
        { "species": "Cmac", "chrom": "C_mac_LG01", "start_bp": 8100000, "end_bp": 11300000, "strand": "-" }
      ],
      "support": { "wfmash": true, "miniprot": true, "n_anchors": 142 }
    }
  ],
  "breakpoints_multilineage": [
    {
      "bp_id": "mlbp_0001",
      "position": { "species": "Cgar", "chrom": "C_gar_LG28", "pos_bp": 15115000 },
      "lineage_distribution": {
        "Cgar": "boundary_present",
        "Cmac": "boundary_absent_internal_to_block",
        "Hetero": "boundary_present"
      },
      "interpretation_default": "boundary_reuse_across_lineages"
    }
  ]
}
```

### Atlas-side

- Layer key: `synteny_multispecies`
- New page-16 view mode: "multi-species" (toggle button) that switches the
  ribbon plot from 2-species to all loaded species
- Hypothesis hook: hypotheses with predictions on this layer evaluate
  against the lineage_distribution field
- Architecture class auto-suggest:
  - Class B (synteny-boundary inversion) ⇒ if breakpoint sits at edge of
    a multi-species block boundary
  - Class C (fusion/fission-associated) ⇒ if breakpoint colocates with
    a `boundary_present` in a sister species
  - Class E (recurrent hotspot) ⇒ if `lineage_distribution` shows the
    breakpoint reused across ≥3 lineages

---

## Chat 115e — Miniprot anchors layer

### Why this layer matters

Miniprot uses protein-vs-genome alignment (slow-evolving signal) and
survives across deeper divergences than wfmash (DNA-vs-DNA). For very
distant outgroups, miniprot anchors are the only signal that gives
synteny information. Required for proper Layer-2 age tagging
(distinguishing "young inversion within Cgar" from "old breakpoint
predating Cgar–Cmac split").

### Provisional layer schema

```jsonc
{
  "tool": "miniprot_anchors_v1",
  "schema_version": 1,
  "reference_species": "Cgar",
  "query_proteomes": ["Cmac", "Hetero", "Trichomycterus"],
  "anchors": [
    {
      "anchor_id": "ma_000001",
      "ref_chrom": "C_gar_LG28",
      "ref_start_bp": 14802340, "ref_end_bp": 14802980,
      "protein_id": "Cmac.gene12345",
      "query_species": "Cmac",
      "query_chrom": "C_mac_LG01",
      "query_start_bp": 8101200, "query_end_bp": 8101840,
      "identity": 0.87,
      "is_collinear_with_neighbors": true
    }
  ]
}
```

### Atlas-side

- Layer key: `miniprot_anchors`
- Render: thin tick marks on page-16 ribbons (one per anchor, color by
  query species)
- Hypothesis hook: hypotheses with `predictions[].layer === "miniprot_anchors"`
  evaluate against anchor density inside vs flanking the breakpoint

---

## Chat 115f — Auto-suggest classification from evidence layers

### Behavior

When a user opens the classification editor (chat 115b) for a breakpoint
and presses "Auto-suggest", the atlas:

1. Reads all evidence layers currently loaded
2. For each architecture class A–F, computes a fit score from layer
   evidence:
   - A (simple inversion): high if local PCA bands strong, low if synteny
     boundary present
   - B (synteny-boundary): high if multi-species synteny shows a block
     edge here
   - C (fusion/fission): high if cs_breakpoints carries a fission/fusion
     event at this position
   - D (terminal translocation): high if breakpoint is within 5% of a
     chromosome end AND homologous segment appears on a different
     chromosome in another species
   - E (recurrent hotspot): high if multiple rearrangement signals
     colocate (≥2 of: high TE flank, fission/fusion, multi-species
     boundary reuse)
   - F (ambiguous): high if assembly gap nearby OR low mappability OR
     only one breakpoint visible
3. For each age model, similar fit logic against dxy + multi-species
   layer
4. Returns top architecture + top age model + confidence score
5. Pre-fills the editor; user can accept or override

### Critical caveat

Auto-suggest is **assistive, not authoritative**. The annotator (Quentin)
always has final say. The score is shown but not the only thing surfaced —
the atlas should display *why* each suggestion was made (which evidence
layers contributed) so the annotator can audit.

### Visual style

- Auto-suggest button next to architecture / age dropdowns
- Result shown as: `[D · OLD-BP/YOUNG-INV · 0.72] [accept] [override]`
- Hover the score → tooltip showing layer-by-layer breakdown

---

## Chat 115g — Hypothesis-scoring engine + manuscript table export

### Behavior

For each hypothesis in `state.hypotheses[]`:

1. Iterate over `predictions[]`
2. For each prediction, check if its layer is loaded; if so, evaluate
   the prediction (`expected` value vs actual layer data)
3. Update `evidence_status` on each prediction
4. Compute `evidence_score` = Σ weight × (supports=+1, refutes=−1,
   neutral=0)
5. Compute `evidence_score_normalized` = score / max possible
6. Sort hypotheses by normalized score; flag top one per
   `applies_to.id` as the leading model

### Page-14 UI

- Table view: one row per hypothesis, columns:
  - hypothesis name | applies to | architecture | age | evidence score | rank
- Per-hypothesis expand: shows prediction-by-prediction evidence breakdown
- Filters: by candidate, by class, by score threshold
- Export: TSV (data) + Markdown (manuscript-ready paragraph format)

### Manuscript export format

For each candidate, generate a paragraph like:

> **Candidate cand_lg28_15_18** (LG28 15.115–18.005 Mb): three competing
> models were evaluated. Model B (recurrent breakpoint reuse, Class E,
> age OLD-BP/YOUNG-INV) achieved the highest evidence score
> (normalized = 0.72), supported by elevated TE flanking density and
> presence of the breakpoint at a synteny boundary shared with
> *C. macrocephalus*. Model A (terminal translocation with reversed
> attachment) and Model C (ancestral fission) scored lower
> (0.18 and −0.34 respectively).

---

## Open design questions (decisions not yet made)

1. **Confidence calibration**: should `confidence` field be derived
   automatically from `evidence_score_normalized`, or stay user-set?
   Recommendation: user-set, but show the normalized score as a hint.

2. **Multi-annotator support**: only Quentin annotates today. If
   collaborators come in, should classifications carry annotator name
   + version history? Probably yes, but that's chat 116+.

3. **Cross-cohort hypothesis sharing**: can a hypothesis registered for
   Cgar be evaluated against Cmac data when the Cmac cohort loads?
   Probably yes via `applies_to.kind = "region"` (genomic coordinates
   instead of candidate IDs), but defer this until the Cmac wild cohort
   work begins.

4. **Versioning the classification schema**: as the architecture classes
   or age tags evolve, classifications need a `schema_version` field.
   Initial: schema_version: 1 (the page-5 schema as of turn 113-tail).

5. **Pre-registered hypotheses**: should the LG28 default hypotheses
   ship with the atlas (immutable, atlas-author-curated) or be created
   by Quentin on first load (mutable, user-curated)? Recommendation:
   user-curated. The atlas can offer a "create LG28 default hypotheses"
   button on first load but doesn't auto-create them, so users on
   different species don't see catfish-specific defaults.

---

## Risks the implementing chat should watch for

1. **State pollution across page switches**. Current atlas has multiple
   `state.candidates` / `state.candidates_detailed` registries (turn 87
   work). Hypotheses must NOT silently overwrite each other when
   switching K-modes. Use `assertSameMode` patterns from existing code.

2. **Storage quota**. localStorage limits to 5–10 MB per origin. With
   200 candidates × 5 hypotheses each × moderately-sized predictions
   each, expect ~500 KB. Well under quota, but the user might also have
   large precomp JSONs cached. Test quota-exceeded path.

3. **Cross-page rendering performance**. Auto-suggest runs across all
   loaded layers per breakpoint. With 100+ breakpoints, this is O(N×L)
   where N=breakpoints, L=layers. Cache aggressively; recompute only
   on layer change.

4. **Breakpoint ID stability**. cs-breakpoints get IDs like
   `cs_bp_0001` from STEP_CS01. If the user re-runs STEP_CS01 with
   different parameters, IDs may shift. Classifications keyed by ID
   would silently misalign. Mitigation: when loading a new
   `cs_breakpoints_v1.json`, hash the breakpoints by chrom+pos and
   warn if any classifications no longer match.

---

## Connection to other documents

- Page 5 (help) — "Breakpoint classification framework" section: the
  authoritative schema this spec operationalizes
- `HANDOFF_2026-05-01_session_summary.md` — what was shipped tonight
  (turn 113-tail through 114c-partial), including the existing
  cs-breakpoint cross-page integration the hypothesis registry will
  build on top of
- `ATLAS_REFERENCE_for_phase4b_synthesis.md` — how the existing atlas
  state model works, K-mode handling, schema-v2 layer machinery
- `breakpoint_pipeline/docs/HANDOFF.md` (parallel breakpoint-pipeline
  chat) — the multi-haplotype methodology question that determines
  whether classifications are binary (REF/HET/INV) or multi-class
  (H1/H1...H3/H3). The hypothesis registry should accommodate both
  via a flexible JSON schema.

---

## End

This spec is concrete enough to implement, abstract enough to revise.
The implementing chat should:

1. Read this whole spec
2. Read the relevant section of `HANDOFF_2026-05-01_session_summary.md`
3. Read the page-5 "Breakpoint classification framework" in the atlas
   itself (open `Inversion_atlas.html`, click tab 16, look at the schema)
4. Pick one chat (115a is the obvious starting point)
5. Build, test, ship as one cleanly-scoped turn

Don't try to do all of 115a–115g in one session. The lesson from turn
114 was that 4 small turns ship cleanly while one big turn doesn't.
