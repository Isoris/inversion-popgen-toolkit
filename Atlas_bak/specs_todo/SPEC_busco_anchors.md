# SPEC — BUSCO single-copy proteins as cross-species anchors

**Status**: forward-looking spec. BUSCO faa exists per Quentin, atlas integration not yet built.

**Reading order**: this spec → `SPEC_OVERVIEW_multispecies_architecture.md` for how it composes with miniprot, wfmash, comparative TE, and the phylogenetic tree.

---

## Purpose

BUSCO single-copy protein alignment is a complementary cross-species anchor to miniprot. Where miniprot uses a query species' full proteome to align against another species' genome, BUSCO uses a curated set of universally single-copy genes that are highly conserved across vertebrates / actinopterygians.

Why two protein-anchor signals?

- **Miniprot**: dense, comprehensive, but slow and fails on highly divergent species
- **BUSCO**: sparse (typically 3000-5000 genes for actinopterygians), survives at much deeper divergences, already standardized

For the catfish phylogeny work, BUSCO is likely already the input alignment. For breakpoint analysis, BUSCO anchors are sparser than miniprot but more consistent across distant species.

## What we have

Quentin has BUSCO faa (`.faa` = protein FASTA) files for the catfish species. These were likely generated as part of genome annotation QC (BUSCO completeness scoring) and as the basis for the phylogenetic tree (concatenated single-copy BUSCOs is the standard input for catfish phylogenomics).

## Why BUSCO is useful even though it's sparse at breakpoints

The honest concern: BUSCO genes are housekeeping / highly-conserved, so they tend to live in **conserved gene-rich euchromatin**, NOT in repeat-rich breakpoint regions. So BUSCO anchors directly AT a breakpoint will be rare.

But BUSCO is still useful because:

1. **It anchors the chromosome flanks**, which are usually gene-rich. If chromosome A in species 1 has BUSCO genes X, Y, Z and chromosome B in species 2 has BUSCO genes X, Y, Z — those chromosomes are homologous, even if the breakpoint between them is repeat-rich and has no BUSCO genes inside it.
2. **It survives deep divergence**. Where wfmash gives up on a Cgar-vs-Trichomycterus comparison (~150 Mya divergence), BUSCO still works.
3. **It gives a consistent yardstick**. Same gene set across all species means cross-species comparisons are apples-to-apples.

So BUSCO is not a breakpoint-localization tool. It's a **chromosome-homology backbone** tool. Used together with miniprot (denser, less deep) and wfmash (densest, shallowest), it covers the full divergence range.

---

## Layer schema

```jsonc
// busco_anchors_v1.json
{
  "tool": "busco_anchors_v1",
  "schema_version": 1,
  "lineage_dataset": "actinopterygii_odb10",  // or whichever BUSCO lineage was used
  "n_buscos_in_dataset": 3640,
  "species": [
    {
      "species": "C_gariepinus",
      "assembly": "fClaHyb_Gar_LG",
      "n_complete_single_copy": 3415,
      "completeness_pct": 93.8,
      "anchors": [
        {
          "busco_id": "10004at7898",
          "chrom": "C_gar_LG28",
          "start_bp": 14802340,
          "end_bp": 14805100,
          "strand": "+",
          "gene_id": "internal_gene_xyz",
          "protein_length_aa": 412
        },
        ...
      ]
    },
    {
      "species": "C_macrocephalus",
      "assembly": "fClaHyb_Mac_LG",
      "n_complete_single_copy": 3398,
      "completeness_pct": 93.4,
      "anchors": [...]
    },
    ...
  ],
  "homology_pairs": [
    // optional: precomputed cross-species pairings
    {
      "busco_id": "10004at7898",
      "species_locations": {
        "C_gariepinus": { "chrom": "C_gar_LG28", "pos_bp": 14803720 },
        "C_macrocephalus": { "chrom": "C_mac_LG01", "pos_bp": 8101520 },
        "Heterobranchus_longifilis": { "chrom": "Het_LG14", "pos_bp": 12018800 }
      }
    },
    ...
  ]
}
```

The `homology_pairs` block is optional. If absent, the atlas can compute cross-species BUSCO pairings on the fly (one BUSCO ID = one homology relationship). If present, it's pre-cached.

---

## How the atlas uses this

### Page 16 (cross-species breakpoints)

- BUSCO anchors render as small ticks on the ribbon plot (one per BUSCO, color by species)
- Sparse near breakpoints (because breakpoints are repeat-rich and depleted of BUSCOs) — this is informative: a chromosome region with NO BUSCO anchors AND high TE density is a strong fragility candidate
- Dense in chromosome flanks — these confirm chromosome homology even when miniprot/wfmash fail

### Page 14 (hypothesis registry)

- BUSCO-based homology backbone can confirm chromosome-context for each candidate
- Architecture class auto-suggest (Class C = fusion/fission-associated) uses BUSCO topology: if Cgar LG28 has BUSCOs A-B-C-D-E and Cmac has A-B-C on LG01 + D-E on LG06, that's a clear fission, anchored by the BUSCO backbone
- Class B (synteny-boundary inversion) similarly: BUSCOs on either side of the breakpoint are syntenic across species

### Page 13/16 cross-species view

- Toggle "show BUSCO backbone" — overlays BUSCO anchors as a thin track at the top of the ribbon plot
- Useful sanity check: if BUSCO anchors are properly collinear in the comparison, the assemblies are reliable in this region; if BUSCOs are scrambled, the assembly may have errors

---

## Pipeline-side: extracting BUSCO anchors per species

Quentin already has the faa files. Adding genomic coordinates needs:

1. Map each BUSCO faa to its source genome — usually the BUSCO output directory has this info (`full_table.tsv` from BUSCO has chromosome + start + end per BUSCO ID)
2. Extract the per-species table
3. Standardize chromosome names across species (often inconsistent: `LG28`, `chr28`, `C_gar_LG28`)
4. Pair BUSCOs across species by `busco_id` (same ID = same gene = homologous)
5. Output the JSON above

This is straightforward — BUSCO output is well-structured. ~1 day of script work.

---

## What BUSCO does NOT solve

- **Direct breakpoint localization** — BUSCOs are gene-rich, breakpoints are repeat-rich, they don't overlap
- **Polymorphism in non-resequenced species** — BUSCO is per-assembly, not per-individual
- **Resolving very recent rearrangements** — at < 10 Mya divergence, miniprot and wfmash work better; BUSCO's strength is at deeper divergences

## What BUSCO uniquely enables

- **Phylogenetic tree input** (probably already used)
- **Chromosome-homology backbone at deep divergences** where DNA-level alignment fails
- **Sanity-checking miniprot / wfmash results** — if BUSCOs disagree with miniprot on chromosome homology, something is wrong with one of them
- **Standardized comparison** across many species without re-running expensive aligners

---

## Atlas-side implementation

Suggested chat scope:
- **Chat 115k**: ingest BUSCO faa + position info, build the JSON layer
- **Chat 115l**: render BUSCO ticks on page 16 ribbon plot
- **Chat 115m**: BUSCO-based homology backbone for architecture class auto-suggest

These come AFTER comparative TE + miniprot anchors are wired, since BUSCO complements rather than replaces them.

Effort: ~1 day per chat for the pipeline + atlas side.

---

## Connection to phylogenetic tree

If the existing catfish tree was built from concatenated single-copy BUSCOs (likely), then:
- The BUSCO set used for the tree is the same as the BUSCO set for anchors
- Tree provenance documentation (`SPEC_phylogenetic_tree_integration.md`) and BUSCO anchor documentation share inputs
- Worth checking with Quentin: was the tree built from BUSCOs? If yes, document the locus list once and reference it from both specs.

---

## Manuscript implication

> *"To anchor cross-species comparisons across deep catfish divergences (up to ~150 Mya to the deepest outgroup), we used [N] single-copy BUSCO genes from the actinopterygii_odb10 dataset as homology anchors. BUSCO anchors confirm chromosome-scale homology even where DNA-level alignment loses sensitivity, and their depletion at inversion breakpoints (gene-poor, repeat-rich regions) provides additional evidence that breakpoints sit in evolutionarily distinct chromatin contexts from the surrounding gene-dense regions."*

This is a clean methods-section paragraph that doesn't overclaim.
