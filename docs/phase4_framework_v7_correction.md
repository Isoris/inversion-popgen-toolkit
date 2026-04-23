# Phase 4 framework v7 — corrections after reading `breakpoint_pipeline/`

**This is a correction doc, not a new framework.** After reading the
`breakpoint_pipeline/` tar + `SESSION_AUDIT.md` + `METHODOLOGY.md` +
`WORKFLOW_DIAGRAM.md`, two things are clear:

1. My v6 interior-structure diagnostic was duplicating logic the user
   had already designed and partly built.
2. I was reaching for a DUP-TRP/INV-DUP template from human cancer
   genetics that doesn't fit catfish hatchery inversion data.

This doc fixes both.

---

## The DUP-TRP/INV-DUP framing is wrong for this data

### Why I proposed it
The LG28 dosage heatmap shows interior banding → I reached for a
known structural-variant architecture that produces banded dosage
patterns → DUP-TRP/INV-DUP (Carvalho et al.).

### Why it doesn't fit
- **Size range.** Published DUP-TRP/INV-DUP events are typically
  20–500 kb in human cancer genomics. The LG28 candidate is 2.89 Mb.
  Hatchery polymorphic inversions are predominantly simple balanced
  inversions at the 100-kb to low-Mb scale.
- **Biological context.** DUP-TRP/INV-DUP is a copy-number
  rearrangement involving triplications flanked by duplications,
  typically arising in tumour genomes via BIR-mediated mechanisms.
  It is not the expected architecture for segregating hatchery
  polymorphisms under standard broodstock management.
- **Data type.** DUP-TRP/INV-DUP is diagnosed by aCGH copy-number
  patterns across 3+ dosage levels. Our dosage is BEAGLE posterior
  expected count, scaled 0–2 for a biallelic polymorphism. The
  biology we measure is allele-frequency variation at an inversion
  polymorphism, not copy-number variation.

**Conclusion:** DUP-TRP/INV-DUP is out-of-scope for this project.
Any candidate whose interior structure genuinely requires that
vocabulary to describe should be **gated out of the main analysis**,
not reconstructed.

### The honest gating rule

Add a terminal label:

  `complex_rearrangement_out_of_scope`

Triggered when interior diagnostic hits either:
- `q2_interior_structure_class = complex_unresolved` AND
  `q2_interior_n_changepoints >= 3`
- OR carrier heatmap shows three or more dosage levels distinct
  from 0/1/2 on BEAGLE scale (this is a red flag for copy-number
  rearrangement rather than balanced inversion)

Candidates receiving this label are:
- Retained in the catalog for transparency
- Excluded from mechanism, age, and cross-species interpretation
- Flagged in manuscript Table S1 with the note "complex architecture
  not consistent with simple inversion model; excluded from further
  interpretation"
- Not forced into `supported_balanced_inversion_*` labels

**This is the honest handling.** Reconstruction of complex
rearrangements from ~9× short-read data is not realistic for this
manuscript. The user was right: this is a gating decision, not a
reconstruction problem.

### The manuscript framing for complex candidates

"Of X candidates in our catalog, Y were flagged as having
interior architecture inconsistent with a simple balanced-inversion
model. These candidates are retained as Table Sz and excluded from
mechanism, age, and conservation analyses. Characterization of
complex rearrangements would require long-read haplotype-resolved
assembly per carrier, which is beyond the scope of this study."

One sentence. Honest. No reviewer will argue with it.

---

## The existing breakpoint_pipeline supersedes my v6 interior work

### What the user already has (confirmed from file reading)

`breakpoint_pipeline/` (7 scripts, ~3800 lines total):

| Script | Purpose |
|---|---|
| `01_dosage_signal.R` | informative markers, core block, block boundary extension (adapted from C01i) |
| `02_ancestral_fragments.R` | per-carrier ancestral-fragment scan + KDE mode + bootstrap CI (the primary breakpoint signal, weight 3.0) |
| `03_consensus_merge.R` | weighted-median consensus across 7 method sources, SV at weight 0.5 |
| `04_diagnostic_figure.R` | multi-track per-candidate PDF (7 tracks including fragment boundary RUG + histogram inset) |
| `05_four_encoding_diagnostic.R` | parallel sample comparison under 4 encodings (minor / major / 012 / polarized) — clustering robustness test |
| `06_regime_stream_graph.R` | regime / ancestry stream-graph visualization |
| `07_breakpoint_evidence_pileup.py` | per-sample stacked read evidence (split reads, discordant pairs, coverage) — Panel D |

The design document `METHODOLOGY.md` is explicit:

> **C1.** The union of all per-sample ancestral fragment boundaries
> on one side of the inversion clusters around the true breakpoint on
> that side. Samples whose ancestors never recombined vote for the
> true breakpoint position. Samples whose ancestors recombined inside
> the inversion vote for a position interior to the true breakpoint.

This is the **primary signal**. Fragment-distribution **mode** is the
breakpoint. Fragment-distribution **tail** quantifies recombinant
history. Fragment-distribution **spread** gives the CI.

### What my v6 `interior_structure_diagnostic.py` was doing

I was re-deriving the same signal using changepoints on per-group mean
dosage + Hartigan dip on HOM1 intragroup correlations. This is a
**coarser, less principled** version of the same diagnostic.

The user's approach is better because:
- It operates at SNP resolution, not window resolution (~230 bp
  spacing on LG28, vs my arbitrary position bins).
- It gives per-sample fragments explicitly — each sample has a
  documented left/right boundary, not just a group summary.
- The recombinant-tail distribution IS the "interior complexity"
  signal. No separate Hartigan dip needed — if the left-side
  fragment distribution has a bimodal shape, there are two
  sub-populations of ancestral fragments.
- It has native CI via bootstrap, which my version doesn't.

### What v6 should actually do

**Nothing new on the interior.** Use the existing breakpoint_pipeline
outputs. The keys that come for free from the existing `02_ancestral_fragments.R`:

```
q2_fragment_left_mode_bp              from summarize_fragment_side()
q2_fragment_left_ci_low_bp            bootstrap 2.5%
q2_fragment_left_ci_high_bp           bootstrap 97.5%
q2_fragment_left_ci_width_kb          → precision class derived from this
q2_fragment_right_mode_bp             analogous
q2_fragment_right_ci_low_bp
q2_fragment_right_ci_high_bp
q2_fragment_right_ci_width_kb
q2_fragment_mad_kb_left               recombinant-tail breadth
q2_fragment_mad_kb_right
q2_n_fragment_carriers                carrier count used
q2_fragment_left_distribution_shape   derived: unimodal | bimodal | right-skewed
q2_fragment_right_distribution_shape
q2_fragment_interior_recomb_frequency fraction of carriers with fragment
                                       extending less than 80% of block width
```

Two derived summary keys:

```
q2_interior_class                     clean_simple |
                                       edge_recombinants_only |
                                       bimodal_boundary_signal |
                                       deep_interior_recombinants |
                                       complex_rearrangement_out_of_scope

q2_breakpoint_precision_class         snp_resolution (CI < 5 kb)
                                       sub_5kb (CI 5-20 kb)
                                       sub_30kb (CI 20-100 kb)
                                       coarse (CI > 100 kb)
```

All derived from the existing breakpoint_pipeline outputs. No new
script.

### What I should delete from v6

- `interior_structure_diagnostic.py` — SUPERSEDED by the user's
  `02_ancestral_fragments.R`. Delete.
- `interior_structure.schema.json` — keep the schema block-type name
  but point its `source_script` at `02_ancestral_fragments.R` and
  redefine the properties to match the fragment-distribution output.

### What I should add to v6 instead

A **thin bridge script** that reads the three existing breakpoint_pipeline
output TSVs and writes the derived v6 keys into the registry:

- `bp_pipeline_bridge.py` — reads `candidate_breakpoints_consensus.tsv`,
  `candidate_ancestral_fragments.tsv`, `candidate_breakpoints_per_method.tsv`
  and writes the 13 derived keys listed above.

~150 lines. No new biology, no new statistics.

---

## Updated structural-class assigner

`assign_structural_class_v6.py` needs two changes:

1. Add the gate label `complex_rearrangement_out_of_scope`. Fires
   when `q2_interior_class = complex_rearrangement_out_of_scope`
   is present. Terminal — no further labelling.

2. Rename `supported_balanced_inversion_with_interior_complexity`
   to `supported_balanced_inversion_with_edge_recombinants` —
   honest about what the interior signal actually means. Fragment
   distributions with wide spread on one side indicate edge
   recombinants, not necessarily nested inversions.

3. The interior class feeds into the suffix rather than being a
   terminal override:
   - `clean_simple` → no suffix modifier
   - `edge_recombinants_only` → no suffix modifier (edges are normal
     in segregating inversions — this is a feature not a bug)
   - `bimodal_boundary_signal` → `_with_partial_interior_recombination`
     (two sub-populations of fragments — one side shows historical
     recombinants pushed deep into the interior)
   - `deep_interior_recombinants` → `_with_deep_recombinants` (whole
     fragment distribution is wide; may indicate gene conversion
     tracts or historical admixture within the inversion)
   - `complex_rearrangement_out_of_scope` → terminal label, overrides
     everything else

Updated label list (16 labels, alphabetical):

```
candidate_weak_evidence
complex_rearrangement_out_of_scope                        ← NEW terminal
complex_unresolved_locus
contradicted_candidate
diffuse_region_of_interest
family_confounded_locus
single_family_polymorphism
supported_balanced_inversion_NAHR_like_supported_by_assembly
supported_balanced_inversion_NAHR_like_hypothesis
supported_balanced_inversion_NHEJ_like_supported_by_assembly
supported_balanced_inversion_NHEJ_like_hypothesis
supported_balanced_inversion_simple
supported_balanced_inversion_with_substrate_mechanism_unresolved
supported_balanced_inversion_with_deep_recombinants       ← renamed
supported_balanced_inversion_with_partial_interior_recombination  ← NEW
supported_fixed_inversion_SV_only
supported_locus_unreconstructed
supported_nested_inversion
supported_shared_between_species_inversion
```

The suffix ordering now matches biological reality: `simple <
with_edge_recombinants < with_partial_interior_recombination <
with_deep_recombinants < complex_out_of_scope`.

---

## What v6 Kuang-style annotations still need (unchanged from previous doc)

These still come from a different layer (boundary-region annotation
of reference sequence) and don't overlap with the breakpoint_pipeline:

- LOESS-smoothed repeat density ±100 kb
- Dyad symmetry at boundaries
- Centromeric / satellite hits at boundaries
- StainedGlass-style internal alignment
- Mappability class per boundary zone

~300 lines R. Still TODO, not affected by this correction.

---

## What v6 cross-species still needs (unchanged)

`cross_species_bridge_v6.py` is correct as written. Primary = overlap,
secondary = flank coherence, tertiary = tree polarization, Dollo as
cross-check. Unchanged by this correction.

---

## Session mandate — aligned with user's D19

From `SESSION_AUDIT.md` §D19:

> user's explicit direction for next chat: "wire everything to the
> registries and api so it collects all genome locations and help each
> script to work together so our results is all automatic bc otherwise
> its too tired."
>
> **The trap to avoid**: "let's improve X while wiring." No. Every
> script stays as-is. Wiring is mechanical and must finish in one
> session.

**I was walking into this trap.** The v6 interior-structure diagnostic
was a new module while the user was asking for wiring. Correcting now.

**What this session's v6 delivery should actually produce:**

1. A bridge script `bp_pipeline_bridge.py` that reads the 3 existing
   breakpoint_pipeline TSVs and writes 13 derived keys + a
   block `fragment_distribution.json`.
2. A small patch to `assign_structural_class_v6.py` adding the 3 new
   labels (`complex_rearrangement_out_of_scope`, the two renamed
   interior labels).
3. An updated schema `fragment_distribution.schema.json`.
4. No new interior-structure-diagnostic logic — the user's pipeline
   already does it correctly.

---

## Toolkit audit — one update

The previous toolkit audit (`docs/toolkit_audit.md` in v6) listed
items to merge/wire/simplify. Add one more item:

**Merge candidate 5: v6 `interior_structure_diagnostic.py`
(proposed) → DELETE in favour of `breakpoint_pipeline/02_ancestral_fragments.R`**

Don't ship v6's interior diagnostic. Use what already exists.

---

## Phase 7 cargo observation

The uploaded `phase_7_cargo` tar is a **separate module** (MODULE_6_Cargo
per the INTEGRATION_NOTES) for burden scoring, gene annotation, breeding
plots, and population-level variant tables. Not directly Phase 4 work.
I note it exists and that v6 structural-class assigner should eventually
feed into `phase_7_cargo` via `results_registry` — but that's the wiring
session user wants, not this design doc.

Breeding scripts (`BREEDING_A/C/D`) are downstream consumers of the
inversion catalog. They need the final `q_overall_structural_class`
labels to work. When wiring happens, feeding v6's `final_label.json`
into the `BREEDING_A_broodstock_compatibility.R` input path is one of
the concrete integrations to make.

---

## Summary table — what the user already built, what I was proposing,
## what should ship

| Capability | User's existing | v6 (proposed) | Ship in this session |
|---|---|---|---|
| Breakpoint refinement at SNP resolution | `breakpoint_pipeline/01+02+03` | N/A | user's — unchanged |
| Per-carrier fragment distribution | `02_ancestral_fragments.R` | N/A | user's — unchanged |
| Fragment-mode breakpoint + bootstrap CI | `02_ancestral_fragments.R` | N/A | user's — unchanged |
| Consensus across 7 methods | `03_consensus_merge.R` | N/A | user's — unchanged |
| Interior complexity signal | fragment distribution shape + recombinant tail | my proposed diagnostic | DELETE my version; derive from user's |
| DUP-TRP/INV-DUP reconstruction | N/A | my proposed | DELETE; gate as out-of-scope |
| 4-encoding robustness | `05_four_encoding_diagnostic.R` | N/A | user's — unchanged |
| Per-sample read pileup Panel D | `07_breakpoint_evidence_pileup.py` | N/A | user's — unchanged |
| Kuang-style boundary annotations | N/A | TODO | ship after wiring |
| Cross-species bridge (v6) | N/A | `cross_species_bridge_v6.py` | unchanged from previous |
| Assembled-junction parser (v5) | N/A | `cheat29b_assembled_junction.py` | unchanged from v5 |
| Single-sided BND (v5) | N/A | `bnd_sided_support.py` | unchanged from v5 |
| Final label synthesis | N/A | `assign_structural_class_v6.py` | update with new labels + gating |
| Bridge to breakpoint_pipeline | N/A | needs new script | **ship this session** |

Four things to ship:

1. `bp_pipeline_bridge.py` — reads user's 3 TSVs, writes 13 derived
   keys.
2. `fragment_distribution.schema.json` — new schema for the
   bridge output.
3. Updated `assign_structural_class_v6.py` — 3 new labels, gating.
4. Updated `phase4_framework_v6.md` — reflects everything above.

Nothing else new this session.
