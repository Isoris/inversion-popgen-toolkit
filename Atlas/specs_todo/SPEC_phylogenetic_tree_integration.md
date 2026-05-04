# SPEC — Phylogenetic tree integration

**Status**: forward-looking spec. Tree exists (per Quentin), atlas integration not yet built.

**Reading order**: this spec → `SPEC_OVERVIEW_multispecies_architecture.md` for how it fits the multi-species evidence stack.

---

## Purpose

A phylogenetic tree of catfish species exists. The atlas should use it for two things:

1. **Visualize species relationships** alongside cross-species evidence on page 16 / page 14
2. **Place rearrangement events on branches** — when a structural rearrangement is detected as shared between species, the tree tells you whether it's ancestral, lineage-specific, or convergent

The tree is the connective tissue that turns "species X has feature Y" into "feature Y emerged on branch Z." Without it, comparative evidence is just a list; with it, comparative evidence becomes evolutionary inference.

## What we have

Quentin has a phylogenetic tree (Newick, presumably from a multi-locus or whole-genome ML/Bayesian analysis). Exact provenance not documented in this session — needs to be captured when the tree is loaded.

## What we don't have

- The tree's input alignment / locus set (manuscript needs this)
- Branch support values explicitly recorded (may be in the Newick already)
- Calibrated divergence times (may exist; check)

These don't block atlas integration but they do affect what claims the tree can support.

---

## Atlas-side data model

```jsonc
// phylo_tree_v1.json
{
  "tool": "phylo_tree_v1",
  "schema_version": 1,
  "species_set": ["C_gariepinus", "C_macrocephalus", "Heterobranchus_longifilis", "...", "Trichomycterus_rosablanca"],
  "newick": "(((C_gariepinus:0.05,C_macrocephalus:0.06):0.10,Heterobranchus_longifilis:0.18):0.30,Trichomycterus_rosablanca:0.85);",
  "calibration": {
    "unit": "substitutions_per_site",  // or "Mya" if calibrated
    "calibrated": false,
    "node_ages_mya": null              // optional dict if calibrated
  },
  "support_values": {
    "type": "bootstrap_proportion",    // or "posterior_probability" / null
    "encoded_in_newick": true
  },
  "provenance": {
    "method": "iqtree2 / mrbayes / etc.",
    "alignment": "single_copy_busco_concatenated_or_other",
    "n_loci": 245,
    "model": "GTR+G+I",
    "reference_paper": "Quentin's analysis or external source"
  }
}
```

Atlas reads this as a layer (key: `phylo_tree`).

---

## Page integration

### Page 16 (cross-species breakpoints)

- Render a small tree mini-view in the focus panel header showing only the species relevant to the current breakpoint
- Highlight branches where the breakpoint has supporting evidence vs where it's absent
- Lineage-specific events show as a single highlighted branch
- Ancestral events show as the deepest branch where the feature is present
- Convergent events (rare but informative) show as multiple non-monophyletic branches

### Page 14 (stats profile / hypothesis registry)

- For each candidate, the tree appears as part of the evidence display
- Architecture class auto-suggest (see `SPEC_hypothesis_registry_and_multispecies.md` chat 115f) reads tree topology to distinguish:
  - LINEAGE_KARYO — fixed in one lineage
  - OLD_POLY — shared across deep splits
  - PREDICTED_POLYMORPHIC_HOTSPOT — homologous fragility on multiple branches

### New page (or part of page 14): Tree viewer

- Standalone tree viewer panel, expandable
- Click a species → highlight its evidence rows in cs_breakpoints + classifications
- Click a node → see what events are inferred on that branch

---

## Ancestral state reconstruction

When the tree + per-species presence/absence of a feature is loaded, the atlas can infer where on the tree the feature emerged. Two methods:

1. **Parsimony** — simple, fast, fine for binary features (boundary present/absent, TE-rich yes/no). Use for first-pass inference.
2. **ML / Bayesian** — better for continuous features (TE density values). Use only if needed for manuscript-grade claims.

For atlas v1, parsimony is enough. ML reconstruction can be left to external tools (ape, phytools) that the user runs offline and feeds back via a `phylo_ancestral_states_v1.json` layer.

---

## What the tree DOES enable

- Distinguishing **ancestral** rearrangements from **lineage-specific** rearrangements
- Placing **breakpoint reuse** events on branches (so the manuscript can claim "this breakpoint was reused on the C. gariepinus branch independently of the breakpoint reuse on the Pangasius branch" if topology supports it)
- Estimating **age** in relative terms (deeper branch = older event)
- Calibrating **age in absolute terms** (Mya) IF the tree is calibrated

## What the tree does NOT enable

- Polymorphism inference in non-resequenced species — tree is per-species, not per-individual
- Direct breakpoint coordinate inference for ancestral genomes — that needs ancestral-genome reconstruction, which is much harder
- Resolving rapid radiations or ILS-confounded branches without explicit gene-tree analysis

---

## Connection to BUSCO anchors spec

The single-copy BUSCO protein alignment is one of the standard inputs for building a catfish phylogenetic tree. So the BUSCO faa Quentin already has likely IS the alignment behind the tree (or could be).

If the existing tree was built from BUSCOs:
- The same BUSCOs serve as cross-species anchors (see `SPEC_busco_anchors.md`)
- One input, two evidence layers
- This is unusually clean and worth exploiting

If the existing tree was built from something else:
- BUSCOs can still be used as anchors
- The tree provenance just needs to be documented separately

---

## Atlas-side implementation

Suggested chat scope for the phylogenetic tree work:
- **Chat 115h**: load + render tree, basic viewer panel, lineage highlighting
- **Chat 115i**: ancestral state reconstruction (parsimony only) for binary features
- **Chat 115j**: tree + ML ancestral states for continuous features (deferred unless manuscript-needed)

These come AFTER the comparative TE layer + hypothesis registry, since the tree is most useful when there's evidence to display on it.

Effort: ~1-2 days per chat. Not a major architectural lift; mostly UI work + a small phylo library (use `phylocanvas` or similar in-browser, or pre-compute in Python and ship JSON).

---

## Manuscript implication

With the tree + comparative evidence layers:

> *"By placing inversion-associated breakpoints on the catfish phylogeny, we identified [N] breakpoints that appear ancestral (shared across multiple lineages) and [N] that appear specific to the C. gariepinus lineage. Lineage-specific breakpoints in repeat-rich regions are consistent with the recurrent-rearrangement-hotspot model: the same fragile substrate is reused independently in different lineages."*

Don't claim more than the topology supports. With ~10 species and limited deep branches, you can confidently distinguish "C. gariepinus + C. macrocephalus" from "all catfish" but not much finer. Calibrate language to data resolution.
